#include "search.h"

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/random/random.h"
#include "absl/status/statusor.h"
#include "core/kmer.h"
#include "core/kmer_counter.h"
#include "core/kmer_set.h"
#include "flags.h"
#include "io.h"
#include "log.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(int, k, 15, "the length of kmers");
ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(int, workers, 1, "number of workers");
ABSL_FLAG(int, n, 1, "number of iterations");
ABSL_FLAG(std::string, decompressor, "", "decompressor for input files");
ABSL_FLAG(int, cutoff, 1, "cutoff threshold");

template <int K, int N, typename KeyType>
void Main(const std::string& file1, const std::string& file2) {
  InitDefaultLogger();
  if (absl::GetFlag(FLAGS_debug)) EnableDebugLogs();

  const int n_workers = absl::GetFlag(FLAGS_workers);
  const std::string decompressor = absl::GetFlag(FLAGS_decompressor);

  std::vector<std::string> reads1;
  std::vector<std::string> reads2;

  {
    // Get reads from a FASTA file.
    const auto Read = [&](const std::string& file) {
      spdlog::info("reading file: {}", file);

      absl::StatusOr<std::vector<std::string>> statusor =
          ReadFASTAFile(file, decompressor, n_workers);

      if (!statusor.ok()) {
        spdlog::error("failed to read file: {}", statusor.status().ToString());
        std::exit(1);
      }

      return std::move(statusor).value();
    };

    reads1 = Read(file1);
    reads2 = Read(file2);
  }

  const std::int64_t n_reads = reads1.size();

  if (n_reads != static_cast<std::int64_t>(reads2.size())) {
    spdlog::error("invalid file");
    std::exit(1);
  }

  spdlog::info("n_reads = {}", n_reads);

  KmerSet<K, N, KeyType> kmer_set;

  {
    const int cutoff = absl::GetFlag(FLAGS_cutoff);

    // Constructs a KmerSet from reads.
    const auto GetKmerSet = [&](const std::vector<std::string>& reads) {
      KmerSet<K, N, KeyType> kmer_set;

      const KmerCounter<K, N, KeyType> kmer_counter =
          KmerCounter<K, N, KeyType>::FromReads(reads, false, n_workers);

      std::tie(kmer_set, std::ignore) =
          kmer_counter.ToKmerSet(cutoff, n_workers);

      return kmer_set;
    };

    kmer_set = GetKmerSet(reads1).Add(GetKmerSet(reads2), n_workers);
  }

  spdlog::info("kmer_set.Size() = {}", kmer_set.Size());

  spdlog::info("constructing kmer_graph");
  const KmerGraph<K> kmer_graph = ConstructKmerGraph(kmer_set, n_workers);
  spdlog::info("constructed kmer_graph");

  absl::InsecureBitGen bitgen;

  int n = absl::GetFlag(FLAGS_n);

  while (n--) {
    spdlog::info("n = {}", n);

    Kmer<K> start, goal;

    // Finds start and goal randomly.
    while (true) {
      std::int64_t j = absl::Uniform(bitgen, 0, n_reads);

      spdlog::info("finding goal and start: j = {}", j);

      const std::string& read1 = reads1[j];
      const std::string& read2 = reads2[j];

      // Finds a kmer in "read" that exists in "kmer_graph".
      const auto FindKmer = [&](const std::string& read) {
        for (std::size_t i = 0; i < read.length() - K + 1; i++) {
          const std::string s = read.substr(i, K);

          // "N" should not be contained.
          if (s.find('N') == std::string::npos) continue;

          const Kmer<K> kmer(s);

          // "kmer" should be in the graph.
          if (kmer_graph.ids.find(kmer) == kmer_graph.ids.end()) continue;

          return std::make_optional(kmer);
        }

        return std::optional<Kmer<K>>();
      };

      {
        std::optional<Kmer<K>> optional = FindKmer(read1);
        if (!optional.has_value()) continue;
        start = std::move(optional).value();
      }

      {
        std::optional<Kmer<K>> optional = FindKmer(read2);
        if (!optional.has_value()) continue;
        goal = std::move(optional).value();
      }

      break;
    }

    spdlog::info("start = {}, goal = {}", start.String(), goal.String());

    {
      spdlog::info("executing dijkstra");
      const SearchResult result = DijkstraSearch(kmer_graph, start, goal);
      spdlog::info("executed dijkstra");
      spdlog::info("found = {}, distance = {}, visited_nodes = {}",
                   result.found, result.distance, result.visited_nodes);
      std::cout << result.distance << ' ' << result.visited_nodes;
    }

    {
      spdlog::info("executing fast A*");
      const SearchResult result = AStarSearch<K>(kmer_graph, start, goal);
      spdlog::info("executed fast A*");
      spdlog::info("found = {}, distance = {}, visited_nodes = {}",
                   result.found, result.distance, result.visited_nodes);
      std::cout << ' ' << result.visited_nodes;
    }

    std::cout << std::endl;
  }
}

int main(int argc, char** argv) {
  const std::vector<std::string> files = ParseFlags(argc, argv);

  const int k = absl::GetFlag(FLAGS_k);

  switch (k) {
    case 15:
      Main<15, 14, std::uint16_t>(files[0], files[1]);
      break;
    case 19:
      Main<19, 10, std::uint32_t>(files[0], files[1]);
      break;
    case 23:
      Main<23, 14, std::uint32_t>(files[0], files[1]);
      break;
    default:
      spdlog::error("unsupported k value");
      std::exit(1);
  }
}

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
#include "search.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(int, k1, 23, "the length of kmers");
ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(int, workers, 1, "number of workers");
ABSL_FLAG(int, n, 1, "number of iterations");
ABSL_FLAG(std::string, decompressor, "", "decompressor for input files");
ABSL_FLAG(int, cutoff, 1, "cutoff threshold");
ABSL_FLAG(int, k2, 9, "the length of kmers used for A* search");

template <int K1, int N1, typename KeyType1, int K2, int N2, typename KeyType2>
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

  KmerSet<K1, N1, KeyType1> kmer_set;

  {
    const int cutoff = absl::GetFlag(FLAGS_cutoff);

    // Constructs a KmerSet from reads.
    const auto GetKmerSet = [&](const std::vector<std::string>& reads) {
      KmerSet<K1, N1, KeyType1> kmer_set;

      const KmerCounter<K1, N1, KeyType1> kmer_counter =
          KmerCounter<K1, N1, KeyType1>::FromReads(reads, false, n_workers);

      std::tie(kmer_set, std::ignore) =
          kmer_counter.ToKmerSet(cutoff, n_workers);

      return kmer_set;
    };

    kmer_set = GetKmerSet(reads1).Add(GetKmerSet(reads2), n_workers);
  }

  spdlog::info("kmer_set.Size() = {}", kmer_set.Size());

  absl::flat_hash_map<std::pair<Kmer<K2>, Kmer<K2>>, std::int64_t>
      all_pair_distances;
  {
    spdlog::info("constructing small_kmer_set");
    KmerSet<K2, N2, KeyType2> small_kmer_set =
        kmer_set.template Extract<K2, N2, KeyType2>(n_workers);
    spdlog::info("constructed small_kmer_set");

    spdlog::info("constructing all_pair_distances");
    all_pair_distances = GetAllPairDistances(small_kmer_set, n_workers);
    spdlog::info("constructed all_pair_distances");
  }

  spdlog::info("constructing kmer_graph");
  const KmerGraph<K1> kmer_graph = ConstructKmerGraph(kmer_set, n_workers);
  spdlog::info("constructed kmer_graph");

  absl::InsecureBitGen bitgen;

  int n = absl::GetFlag(FLAGS_n);

  while (n--) {
    spdlog::info("n = {}", n);

    Kmer<K1> start, goal;

    // Finds start and goal randomly.
    while (true) {
      std::int64_t j = absl::Uniform(bitgen, 0, n_reads);

      spdlog::info("finding goal and start: j = {}", j);

      const std::string& read1 = reads1[j];
      const std::string& read2 = reads2[j];

      // Finds a kmer in "read" that exists in "kmer_graph".
      const auto FindKmer = [&](const std::string& read) {
        for (std::size_t i = 0; i < read.length() - K1 + 1; i++) {
          const std::string s = read.substr(i, K1);

          // "N1" should not be contained.
          if (s.find('N') == std::string::npos) continue;

          const Kmer<K1> kmer(s);

          // "kmer" should be in the graph.
          if (kmer_graph.ids.find(kmer) == kmer_graph.ids.end()) continue;

          return std::make_optional(kmer);
        }

        return std::optional<Kmer<K1>>();
      };

      {
        std::optional<Kmer<K1>> optional = FindKmer(read1);
        if (!optional.has_value()) continue;
        start = std::move(optional).value();
      }

      {
        std::optional<Kmer<K1>> optional = FindKmer(read2);
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
      spdlog::info("executing A*");
      const SearchResult result =
          AStarSearch<K1, K2>(kmer_graph, start, goal, all_pair_distances);
      spdlog::info("executed A*");
      spdlog::info("found = {}, distance = {}, visited_nodes = {}",
                   result.found, result.distance, result.visited_nodes);
      std::cout << ' ' << result.visited_nodes;
    }

    std::cout << std::endl;
  }
}

int main(int argc, char** argv) {
  const std::vector<std::string> files = ParseFlags(argc, argv);

  if (files.size() != 2) {
    spdlog::error("invalid number of arguments");
    return 1;
  }

  const int k1 = absl::GetFlag(FLAGS_k1);
  const int k2 = absl::GetFlag(FLAGS_k2);

  if (k1 == 23) {
    const int K1 = 23;
    const int N1 = 14;
    using KeyType1 = std::uint32_t;

    if (k2 == 9) {
      Main<K1, N1, KeyType1, 9, 10, std::uint8_t>(files[0], files[1]);
      return 0;
    }

    if (k2 == 10) {
      Main<K1, N1, KeyType1, 10, 12, std::uint8_t>(files[0], files[1]);
      return 0;
    }

    if (k2 == 11) {
      Main<K1, N1, KeyType1, 11, 14, std::uint8_t>(files[0], files[1]);
      return 0;
    }

    if (k2 == 12) {
      Main<K1, N1, KeyType1, 12, 8, std::uint16_t>(files[0], files[1]);
      return 0;
    }
  }

  spdlog::error("unsupported k1 and k2 values");

  return 1;
}

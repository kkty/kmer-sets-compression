#include "search.h"

#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/random/random.h"
#include "absl/status/statusor.h"
#include "core/kmer.h"
#include "core/kmer_counter.h"
#include "core/kmer_set.h"
#include "core/range.h"
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

template <int K, typename KeyType>
void Main(const std::string& file1, const std::string& file2) {
  InitDefaultLogger();
  if (absl::GetFlag(FLAGS_debug)) EnableDebugLogs();

  const int n_workers = absl::GetFlag(FLAGS_workers);
  const std::string decompressor = absl::GetFlag(FLAGS_decompressor);

  std::vector<std::string> reads1;
  std::vector<std::string> reads2;

  {
    spdlog::info("reading file: {}", file1);

    absl::StatusOr<std::vector<std::string>> statusor =
        ReadFASTAFile(file1, decompressor, n_workers);

    if (!statusor.ok()) {
      spdlog::error("failed to read file: {}", statusor.status().ToString());
      std::exit(1);
    }

    reads1 = std::move(statusor).value();
  }

  {
    spdlog::info("reading file: {}", file2);

    absl::StatusOr<std::vector<std::string>> statusor =
        ReadFASTAFile(file2, decompressor, n_workers);

    if (!statusor.ok()) {
      spdlog::error("failed to read file: {}", statusor.status().ToString());
      std::exit(1);
    }

    reads2 = std::move(statusor).value();
  }

  const int64_t n_reads = reads1.size();

  if (n_reads != (int64_t)reads2.size()) {
    spdlog::error("invalid file");
    std::exit(1);
  }

  spdlog::info("n_reads = {}", n_reads);

  KmerSet<K, KeyType> kmer_set;

  {
    const int cutoff = absl::GetFlag(FLAGS_cutoff);

    KmerSet<K, KeyType> kmer_set1;
    KmerSet<K, KeyType> kmer_set2;

    const KmerCounter<K, KeyType> kmer_counter1 =
        KmerCounter<K, KeyType>::FromReads(reads1, false, n_workers);
    const KmerCounter<K, KeyType> kmer_counter2 =
        KmerCounter<K, KeyType>::FromReads(reads2, false, n_workers);

    std::tie(kmer_set1, std::ignore) =
        kmer_counter1.ToKmerSet(cutoff, n_workers);
    std::tie(kmer_set2, std::ignore) =
        kmer_counter2.ToKmerSet(cutoff, n_workers);

    kmer_set = std::move(kmer_set1);
    kmer_set.Add(kmer_set2, n_workers);
  }

  spdlog::info("kmer_set.Size() = {}", kmer_set.Size());

  spdlog::info("constructing kmer_graph");
  const KmerGraph<K> kmer_graph = ConstructKmerGraph(kmer_set, n_workers);
  spdlog::info("constructed kmer_graph");

  absl::InsecureBitGen bitgen;

  int n = absl::GetFlag(FLAGS_n);

  while (true) {
    Kmer<K> start, goal;

    // Finds start and goal.
    {
      int64_t j = absl::Uniform(bitgen, 0, n_reads);
      const std::string& read1 = reads1[j];
      const std::string& read2 = reads2[j];

      {
        bool found = false;

        for (size_t i = 0; i < read1.length() - K + 1; i++) {
          const std::string s = read1.substr(i, K);

          // "N" should not be contained.
          if (s.find("N") == std::string::npos) continue;

          const Kmer<K> kmer(s);

          // "kmer" should be in the graph.
          if (kmer_graph.ids.find(kmer) == kmer_graph.ids.end()) continue;

          found = true;
          start = kmer;
          break;
        }

        if (!found) {
          continue;
        }
      }

      {
        bool found = false;

        for (size_t i = 0; i < read2.length() - K + 1; i++) {
          const std::string s = read2.substr(i, K);

          // "N" should not be contained.
          if (s.find("N") == std::string::npos) continue;

          const Kmer<K> kmer(s);

          // "kmer" should be in the graph.
          if (kmer_graph.ids.find(kmer) == kmer_graph.ids.end()) continue;

          found = true;
          goal = kmer;
          break;
        }

        if (!found) {
          continue;
        }
      }
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

    for (int l = 2; l < 10; l++) {
      spdlog::info("l = {}", l);
      spdlog::info("executing A*");
      const SearchResult result = AStarSearch<K>(kmer_graph, start, goal, l);
      spdlog::info("executed A*");
      spdlog::info("found = {}, distance = {}, visited_nodes = {}",
                   result.found, result.distance, result.visited_nodes);
      std::cout << ' ' << result.visited_nodes;
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

    n--;
  }
}

int main(int argc, char** argv) {
  const std::vector<std::string> files = ParseFlags(argc, argv);

  const int k = absl::GetFlag(FLAGS_k);

  switch (k) {
    case 15:
      // 15 * 2 - 16 = 14 bits are used to select buckets.
      Main<15, uint16_t>(files[0], files[1]);
      break;
    case 19:
      // 19 * 2 - 16 = 22 bits are used to select buckets.
      Main<19, uint16_t>(files[0], files[1]);
      break;
    case 23:
      // 23 * 2 - 32 = 14 bits are used to select buckets.
      Main<23, uint32_t>(files[0], files[1]);
      break;
    case 27:
      // 27 * 2 - 32 = 22 bits are used to select buckets.
      Main<27, uint32_t>(files[0], files[1]);
      break;
    default:
      spdlog::error("unsupported k value");
      std::exit(1);
  }
}

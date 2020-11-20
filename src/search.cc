#include "search.h"

#include <iostream>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "core/kmer.h"
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

template <int K, typename KeyType>
void Main(const std::string& file_name) {
  InitDefaultLogger();
  if (absl::GetFlag(FLAGS_debug)) EnableDebugLogs();

  const int n_workers = absl::GetFlag(FLAGS_workers);

  const KmerSet<K, KeyType> kmer_set =
      GetKmerSetFromCompressedKmersFile<K, KeyType>(
          file_name, absl::GetFlag(FLAGS_decompressor), false, n_workers);

  spdlog::info("kmer_set.Size() = {}", kmer_set.Size());

  spdlog::info("constructing kmer_graph");
  const KmerGraph<K> kmer_graph = ConstructKmerGraph(kmer_set, n_workers);
  spdlog::info("constructed kmer_graph");

  for (int i = 0; i < absl::GetFlag(FLAGS_n); i++) {
    Kmer<K> start, goal;
    std::tie(start, goal) = SampleConnectedPair(kmer_graph, 100);

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

    std::cout << std::endl;
  }
}

int main(int argc, char** argv) {
  const std::string file_name = ParseFlags(argc, argv)[0];

  const int k = absl::GetFlag(FLAGS_k);

  switch (k) {
    case 15:
      // 15 * 2 - 16 = 14 bits are used to select buckets.
      Main<15, uint16_t>(file_name);
      break;
    case 19:
      // 19 * 2 - 16 = 22 bits are used to select buckets.
      Main<19, uint16_t>(file_name);
      break;
    case 23:
      // 23 * 2 - 32 = 14 bits are used to select buckets.
      Main<23, uint32_t>(file_name);
      break;
    case 27:
      // 27 * 2 - 32 = 22 bits are used to select buckets.
      Main<27, uint32_t>(file_name);
      break;
    default:
      spdlog::error("unsupported k value");
      std::exit(1);
  }
}

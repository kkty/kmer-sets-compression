#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "core/graph.h"
#include "core/kmer.h"
#include "core/kmer_set.h"
#include "core/kmer_set_set.h"
#include "flags.h"
#include "io.h"
#include "log.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(int, k, 15, "the length of kmers");
ABSL_FLAG(bool, similarity, false, "calculate Jaccard similarity");
ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(std::string, decompressor, "",
          "specify decompressor for input files");
ABSL_FLAG(bool, canonical, false, "count canonical k-mers");
ABSL_FLAG(int, workers, 1, "number of workers");
ABSL_FLAG(bool, parallel_input, false, "read files in parallel");
ABSL_FLAG(int, iteration, 1, "number of iterations for KmerSetSet");
ABSL_FLAG(int, recursion, 1, "recursion limit for KmerSetSetMM");
ABSL_FLAG(bool, check, false, "check if k-mer sets can be reconstructed");
ABSL_FLAG(bool, mm, false, "use maximum weighted matching algorithm");
ABSL_FLAG(bool, approximate_matching, false,
          "use approximate algorithm for matching");
ABSL_FLAG(bool, approximate_weights, false,
          "use approximate weights for the matching algorithm");
ABSL_FLAG(bool, approximate_graph, false,
          "use an approximate (sparse) graph for the matching algorithm");
ABSL_FLAG(bool, partial_matching, false,
          "make mates partially in the matching algorithm");
ABSL_FLAG(std::string, out, "", "path to save dumped file");
ABSL_FLAG(std::string, out_graph, "", "path to save dumped DOT file");
ABSL_FLAG(std::string, out_folder, "", "folder to save dumped files");
ABSL_FLAG(std::string, out_extension, "bin", "extension for output files");
ABSL_FLAG(std::string, compressor, "", "program to compress dumped file");

template <int K, typename KeyType>
void Main(const std::vector<std::string>& files) {
  InitDefaultLogger();

  if (absl::GetFlag(FLAGS_debug)) EnableDebugLogs();

  const int n_workers = absl::GetFlag(FLAGS_workers);
  const bool canonical = absl::GetFlag(FLAGS_canonical);

  const int n_datasets = files.size();

  std::vector<KmerSet<K, KeyType>> kmer_sets(n_datasets);

  // Reads the ith file and constructs kmer_sets[i].
  const auto ReadFile = [&](int i, int n_workers) {
    const std::string& file = files[i];
    const std::string decompressor = absl::GetFlag(FLAGS_decompressor);

    spdlog::info("reading {}", file);

    absl::StatusOr<KmerSet<K, KeyType>> statusor =
        GetKmerSetFromCompressedKmersFile<K, KeyType>(file, decompressor,
                                                      canonical, n_workers);

    if (!statusor.ok()) {
      spdlog::error("failed to read file: {}", statusor.status().ToString());
      std::exit(1);
    }

    kmer_sets[i] = std::move(statusor).value();

    spdlog::info("finished reading {}", file);
  };

  // Reading files in parallel may use a lot of memory.
  if (absl::GetFlag(FLAGS_parallel_input)) {
    boost::asio::thread_pool pool(n_workers);

    for (int i = 0; i < n_datasets; i++) {
      boost::asio::post(pool, [&, i] { ReadFile(i, 1); });
    }

    pool.join();
  } else {
    for (int i = 0; i < n_datasets; i++) {
      ReadFile(i, n_workers);
    }
  }

  {
    int64_t total_size = 0;

    for (int i = 0; i < n_datasets; i++) {
      int64_t kmer_set_size = kmer_sets[i].Size();
      spdlog::info("i = {}, kmer_set_size = {}", i, kmer_set_size);
      total_size += kmer_set_size;
    }

    spdlog::info("total_size = {}", total_size);
  }

  if (absl::GetFlag(FLAGS_similarity)) {
    for (int i = 0; i < n_datasets; i++) {
      for (int j = i + 1; j < n_datasets; j++) {
        double similarity = kmer_sets[i].Similarity(kmer_sets[j], n_workers);
        spdlog::info("i = {}, j = {}, similarity = {}", i, j, similarity);
      }
    }
  }

  if (absl::GetFlag(FLAGS_mm)) {
    spdlog::info("constructing kmer_set_set_mm");

    KmerSetSetMM<K, KeyType> kmer_set_set_mm = [&] {
      if (absl::GetFlag(FLAGS_check)) {
        return KmerSetSetMM<K, KeyType>(
            kmer_sets, absl::GetFlag(FLAGS_recursion),
            absl::GetFlag(FLAGS_approximate_matching),
            absl::GetFlag(FLAGS_approximate_weights),
            absl::GetFlag(FLAGS_approximate_graph),
            absl::GetFlag(FLAGS_partial_matching), n_workers);
      } else {
        // We can move kmer_sets.
        return KmerSetSetMM<K, KeyType>(
            std::move(kmer_sets), absl::GetFlag(FLAGS_recursion),
            absl::GetFlag(FLAGS_approximate_matching),
            absl::GetFlag(FLAGS_approximate_weights),
            absl::GetFlag(FLAGS_approximate_graph),
            absl::GetFlag(FLAGS_partial_matching), n_workers);
      }
    }();

    spdlog::info("constructed kmer_set_set_mm");

    spdlog::info("kmer_set_set_mm.Cost() = {}", kmer_set_set_mm.Cost());

    const std::string out_file = absl::GetFlag(FLAGS_out);

    if (out_file != "") {
      const absl::Status status = kmer_set_set_mm.Dump(
          out_file, absl::GetFlag(FLAGS_compressor), canonical, n_workers);

      if (!status.ok()) {
        spdlog::error("failed to dump kmer_set_mm: {}", status.ToString());
        std::exit(1);
      }
    }

    if (absl::GetFlag(FLAGS_check)) {
      spdlog::info("dumping");
      std::vector<std::string> dumped =
          kmer_set_set_mm.Dump(absl::GetFlag(FLAGS_canonical), n_workers);
      spdlog::info("dumped");

      spdlog::info("loading");
      std::unique_ptr<KmerSetSetMM<K, KeyType>> loaded =
          std::unique_ptr<KmerSetSetMM<K, KeyType>>(
              KmerSetSetMM<K, KeyType>::Load(std::move(dumped),
                                             absl::GetFlag(FLAGS_canonical),
                                             n_workers));
      spdlog::info("loaded");

      for (int i = 0; i < n_datasets; i++) {
        spdlog::info("kmer_sets[{}].Equals(loaded->Get({})) = {}", i, i,
                     kmer_sets[i].Equals(loaded->Get(i, n_workers), n_workers));
      }
    }
  } else {
    spdlog::info("constructing kmer_set_set");

    KmerSetSet<K, KeyType> kmer_set_set;
    const bool check = absl::GetFlag(FLAGS_check);

    // If --check is not specified, it is OK to invalidate kmer_sets.
    {
      const int n_iterations = absl::GetFlag(FLAGS_iteration);

      if (check) {
        kmer_set_set =
            KmerSetSet<K, KeyType>(kmer_sets, n_iterations, n_workers);
      } else {
        kmer_set_set = KmerSetSet<K, KeyType>(std::move(kmer_sets),
                                              n_iterations, n_workers);
      }
    }

    spdlog::info("constructed kmer_set_set");

    // If --out-graph is specified, dumps a DOT file.
    {
      const std::string out_graph = absl::GetFlag(FLAGS_out_graph);

      if (!out_graph.empty()) {
        spdlog::info("dumping graph");

        const absl::Status status = kmer_set_set.DumpGraph(out_graph);

        if (!status.ok()) {
          spdlog::error("failed to dump graph: {}", status.ToString());
        }

        spdlog::info("dumped graph");
      }
    }

    // If --out_file or --out_folder is specified, dumps the structure to a file
    // or files in a folder.
    {
      const std::string out_file = absl::GetFlag(FLAGS_out);
      const std::string out_folder = absl::GetFlag(FLAGS_out_folder);
      const std::string compressor = absl::GetFlag(FLAGS_compressor);

      // If --check is not specified, it is OK to invalidate kmer_set_set.
      const bool clear = !check;

      absl::Status status;

      if (!out_file.empty()) {
        status = kmer_set_set.Dump(out_file, compressor, canonical, clear,
                                   n_workers);

      } else if (!out_folder.empty()) {
        const std::string out_extension = absl::GetFlag(FLAGS_out_extension);

        status = kmer_set_set.DumpToFolder(
            out_folder, compressor, out_extension, canonical, clear, n_workers);
      }

      if (!status.ok()) {
        spdlog::error("failed to dump kmer_set_set: {}", status.ToString());
        std::exit(1);
      }
    }

    if (check) {
      spdlog::info("dumping kmer_set_set");
      std::vector<std::string> dumped =
          kmer_set_set.Dump(canonical, true, n_workers);
      spdlog::info("dumped kmer_set_set");

      spdlog::info("loading");
      KmerSetSet<K, KeyType> loaded =
          KmerSetSet<K, KeyType>::Load(std::move(dumped), canonical, n_workers);
      spdlog::info("loaded");

      for (int i = 0; i < n_datasets; i++) {
        spdlog::info("kmer_sets[{}].Equals(loaded.Get({})) = {}", i, i,
                     kmer_sets[i].Equals(loaded.Get(i, n_workers), n_workers));
      }
    }
  }
}

int main(int argc, char** argv) {
  const std::vector<std::string> files = ParseFlags(argc, argv);

  const int k = absl::GetFlag(FLAGS_k);

  switch (k) {
    case 15:
      // 15 * 2 - 16 = 14 bits are used to select buckets.
      Main<15, uint16_t>(files);
      break;
    case 19:
      // 19 * 2 - 16 = 22 bits are used to select buckets.
      Main<19, uint16_t>(files);
      break;
    case 23:
      // 23 * 2 - 32 = 14 bits are used to select buckets.
      Main<23, uint32_t>(files);
      break;
    case 27:
      // 27 * 2 - 32 = 22 bits are used to select buckets.
      Main<27, uint32_t>(files);
      break;
    default:
      spdlog::error("unsupported k value");
      std::exit(1);
  }
}

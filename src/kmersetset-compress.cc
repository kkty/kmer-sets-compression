#include <algorithm>
#include <cstdlib>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "core/kmer_set.h"
#include "core/kmer_set_immutable.h"
#include "core/kmer_set_set.h"
#include "flags.h"
#include "io.h"
#include "log.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(int, k, 15, "the length of kmers");
ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(std::string, decompressor, "",
          "specify decompressor for input files");
ABSL_FLAG(bool, canonical, false, "count canonical k-mers");
ABSL_FLAG(int, workers, 1, "number of workers");
ABSL_FLAG(bool, parallel_input, false, "read files in parallel");
ABSL_FLAG(int, iteration, 1, "number of iterations for KmerSetSet");
ABSL_FLAG(bool, check, false, "check if k-mer sets can be reconstructed");
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
  const bool check = absl::GetFlag(FLAGS_check);

  const int n_datasets = files.size();

  std::vector<KmerSetImmutable<K, KeyType>> kmer_sets_immutable(n_datasets);

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

    KmerSet<K, KeyType> kmer_set = std::move(statusor).value();

    kmer_sets_immutable[i] = KmerSetImmutable<K, KeyType>(kmer_set, n_workers);

    if (check) {
      if (kmer_set.Equals(kmer_sets_immutable[i].ToKmerSet(n_workers),
                          n_workers)) {
        spdlog::info("kmer_sets_immutable[{}] -> kmer_set: ok", i);
      } else {
        spdlog::error("kmer_sets_immutable[{}] -> kmer_set: failed", i);
        std::exit(1);
      }
    }

    spdlog::info("finished reading {}", file);
  };

  {
    // Reading files in parallel may use a lot of memory.
    const bool parallel_input = absl::GetFlag(FLAGS_parallel_input);

    if (parallel_input) {
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
  }

  {
    int64_t total_size = 0;

    for (int i = 0; i < n_datasets; i++) {
      int64_t size = kmer_sets_immutable[i].Size();
      spdlog::info("i = {}, size = {}", i, size);
      total_size += size;
    }

    spdlog::info("total_size = {}", total_size);
  }

  spdlog::info("constructing kmer_set_set");

  KmerSetSet<K, KeyType> kmer_set_set;

  {
    const int n_iterations = absl::GetFlag(FLAGS_iteration);

    // If --check is specified, kmer_sets_immutable should be copied.
    if (check) {
      kmer_set_set =
          KmerSetSet<K, KeyType>(kmer_sets_immutable, n_iterations, n_workers);
    } else {
      kmer_set_set = KmerSetSet<K, KeyType>(std::move(kmer_sets_immutable),
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
      status =
          kmer_set_set.Dump(out_file, compressor, canonical, clear, n_workers);

    } else if (!out_folder.empty()) {
      const std::string out_extension = absl::GetFlag(FLAGS_out_extension);

      status = kmer_set_set.DumpToFolder(out_folder, compressor, out_extension,
                                         canonical, clear, n_workers);
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
      // NOLINTNEXTLINE
      if (kmer_sets_immutable[i].ToKmerSet(n_workers).Equals(
              loaded.Get(i, n_workers), n_workers)) {
        spdlog::info("loaded.Get({}, n_workers): ok", i);
      } else {
        spdlog::error("loaded.Get({}, n_workers): failed", i);
        std::exit(1);
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

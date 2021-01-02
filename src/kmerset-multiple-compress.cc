#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/usage.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_format.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "core/kmer_set.h"
#include "core/kmer_set_compact.h"
#include "core/kmer_set_immutable.h"
#include "core/kmer_set_set.h"
#include "flags.h"
#include "log.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(int, k, 15, GetFlagMessage("k"));
ABSL_FLAG(bool, debug, false, GetFlagMessage("debug"));
ABSL_FLAG(std::string, decompressor, "", GetFlagMessage("decompressor"));
ABSL_FLAG(int, workers, 1, GetFlagMessage("workers"));
ABSL_FLAG(bool, canonical, true, GetFlagMessage("canonical"));
ABSL_FLAG(std::string, compressor, "", GetFlagMessage("compressor"));

ABSL_FLAG(bool, parallel_input, false, "read files in parallel");
ABSL_FLAG(std::string, out, "", "directory path to save dumped files");
ABSL_FLAG(std::string, extension, "txt", "extension for output files");
ABSL_FLAG(std::string, out_graph, "", "path to save dumped DOT file");

template <int K, int N, typename KeyType>
void Main(const std::vector<std::string>& files) {
  InitDefaultLogger();

  if (absl::GetFlag(FLAGS_debug)) EnableDebugLogs();

  const int n_workers = absl::GetFlag(FLAGS_workers);
  const bool canonical = absl::GetFlag(FLAGS_canonical);

  const int n_datasets = files.size();

  std::vector<KmerSetImmutable<K, N, KeyType>> kmer_sets_immutable(n_datasets);

  // Reads the ith file and constructs kmer_sets_immutable[i].
  const auto ReadFile = [&](int i, int n_workers) {
    const std::string& file = files[i];
    const std::string decompressor = absl::GetFlag(FLAGS_decompressor);

    spdlog::info("reading {}", file);

    spdlog::info("constructing kmer_set_compact");

    KmerSetCompact<K, N, KeyType> kmer_set_compact;

    {
      absl::StatusOr<KmerSetCompact<K, N, KeyType>> statusor =
          KmerSetCompact<K, N, KeyType>::Load(file, decompressor);

      if (!statusor.ok()) {
        spdlog::error("failed to read file: {}", statusor.status().ToString());
        std::exit(1);
      }

      kmer_set_compact = std::move(statusor).value();
    }

    spdlog::info("constructed kmer_set_compact");

    spdlog::info("constructing kmers");

    std::vector<Kmer<K>> kmers = kmer_set_compact.ToKmers(canonical, n_workers);

    spdlog::info("constructed kmers");

    spdlog::info("kmers.size() = {}", kmers.size());

    spdlog::info("constructing kmer_sets_immutable[{}]", i);

    kmer_sets_immutable[i] = KmerSetImmutable<K, N, KeyType>(kmers, n_workers);

    spdlog::info("constructed kmer_sets_immutable[{}]", i);

    spdlog::info("kmer_sets_immutable[{}].Bytes() = {}", i,
                 kmer_sets_immutable[i].Bytes());

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
    std::int64_t total_size = 0;

    for (int i = 0; i < n_datasets; i++) {
      std::int64_t size = kmer_sets_immutable[i].Size();
      spdlog::info("i = {}, size = {}", i, size);
      total_size += size;
    }

    spdlog::info("total_size = {}", total_size);
  }

  spdlog::info("constructing kmer_set_set");

  KmerSetSet<K, N, KeyType> kmer_set_set = KmerSetSet<K, N, KeyType>(
      std::move(kmer_sets_immutable), canonical, n_workers);

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

  // If --out is specified, dumps the structure to files in a directory.
  {
    const std::string out = absl::GetFlag(FLAGS_out);
    const std::string compressor = absl::GetFlag(FLAGS_compressor);
    const std::string extension = absl::GetFlag(FLAGS_extension);

    absl::Status status;

    if (!out.empty()) {
      status = kmer_set_set.Dump(out, compressor, extension, canonical, true,
                                 n_workers);
    }

    if (!status.ok()) {
      spdlog::error("failed to dump kmer_set_set: {}", status.ToString());
      std::exit(1);
    }
  }
}

int main(int argc, char** argv) {
  absl::SetProgramUsageMessage(
      absl::StrFormat("Compresses multiple k-mer sets. Usage: %s [options] "
                      "<paths to file> <path to file> ...",
                      argv[0]));

  const std::vector<std::string> files = ParseFlags(argc, argv);

  const int k = absl::GetFlag(FLAGS_k);

  switch (k) {
    case 15:
      Main<15, 14, std::uint16_t>(files);
      break;
    case 19:
      Main<19, 10, std::uint32_t>(files);
      break;
    case 23:
      Main<23, 14, std::uint32_t>(files);
      break;
    default:
      spdlog::error("unsupported k value");
      std::exit(1);
  }
}

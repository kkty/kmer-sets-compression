#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/status/statusor.h"
#include "core/graph.h"
#include "core/kmer.h"
#include "core/kmer_counter.h"
#include "core/kmer_set.h"
#include "core/kmer_set_compressed.h"
#include "core/kmer_set_set.h"
#include "log.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(int, k, 15, "the length of kmers");
ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(std::string, decompressor, "cat",
          "specify decompressor for FASTQ files");
ABSL_FLAG(std::string, compressor, "cat", "specify compressor for output");
ABSL_FLAG(bool, canonical, false, "count canonical k-mers");
ABSL_FLAG(int, cutoff, 1, "cut off threshold");
ABSL_FLAG(int, workers, 1, "number of workers");
ABSL_FLAG(std::string, in, "", "input file name");
ABSL_FLAG(std::string, out, "", "output file name");

template <int K, typename KeyType>
void Main() {
  InitDefaultLogger();
  if (absl::GetFlag(FLAGS_debug)) EnableDebugLogs();

  const std::string file_name = absl::GetFlag(FLAGS_in);

  const int n_workers = absl::GetFlag(FLAGS_workers);
  const bool canonical = absl::GetFlag(FLAGS_canonical);

  spdlog::info("constructing kmer_counter");

  const KmerCounter<K, KeyType> kmer_counter = [&] {
    const std::string decompressor = absl::GetFlag(FLAGS_decompressor);

    absl::StatusOr<KmerCounter<K, KeyType>> status_or =
        decompressor != "" ? KmerCounter<K, KeyType>::FromFASTQ(
                                 file_name, decompressor, canonical, n_workers)
                           : KmerCounter<K, KeyType>::FromFASTQ(
                                 file_name, canonical, n_workers);

    if (!status_or.ok()) {
      spdlog::error("failed to parse FASTQ file: {}",
                    status_or.status().ToString());
    }

    return *status_or;
  }();

  spdlog::info("constructed kmer_counter");
  spdlog::info("constructing kmer_set");

  const auto [kmer_set, cutoff_count] =
      kmer_counter.ToSet(absl::GetFlag(FLAGS_cutoff), n_workers);

  spdlog::info("constructed kmer_set");
  spdlog::info("kmer_set.Size() = {}", kmer_set.Size());
  spdlog::info("kmer_set.Hash() = {}", kmer_set.Hash(n_workers));

  spdlog::info("constructing kmer_set_compressed");
  const KmerSetCompressed<K, KeyType> kmer_set_compressed =
      KmerSetCompressed<K, KeyType>::FromKmerSet(kmer_set, canonical,
                                                 n_workers);
  spdlog::info("constructed kmer_set_compressed");

  spdlog::info("kmer_set_compressed.Size() = {}",
               kmer_set_compressed.Size(n_workers));

  std::string output_file_name;

  kmer_set_compressed.Dump(absl::GetFlag(FLAGS_out),
                           absl::GetFlag(FLAGS_compressor));
}

int main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);

  const int k = absl::GetFlag(FLAGS_k);

  switch (k) {
    case 15:
      // 15 * 2 - 16 = 14 bits are used to select buckets.
      Main<15, uint16_t>();
      break;
    case 19:
      // 19 * 2 - 16 = 22 bits are used to select buckets.
      Main<19, uint16_t>();
      break;
    case 23:
      // 23 * 2 - 32 = 14 bits are used to select buckets.
      Main<23, uint32_t>();
      break;
    case 27:
      // 27 * 2 - 32 = 22 bits are used to select buckets.
      Main<27, uint32_t>();
      break;
    default:
      spdlog::error("unsupported k value");
      std::exit(1);
  }
}

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
#include "graph.h"
#include "kmer.h"
#include "kmer_counter.h"
#include "kmer_set.h"
#include "kmer_set_compressed.h"
#include "kmer_set_set.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(std::string, decompressor, "cat",
          "specify decompressor for FASTQ files");
ABSL_FLAG(std::string, compressor, "cat", "specify compressor for output");
ABSL_FLAG(bool, canonical, false, "count canonical k-mers");
ABSL_FLAG(int, cutoff, 1, "cut off threshold");
ABSL_FLAG(int, workers, 1, "number of workers");
ABSL_FLAG(std::string, in, "", "input file name");
ABSL_FLAG(std::string, out, "", "output file name");

int main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);

  const std::string file_name = absl::GetFlag(FLAGS_in);

  if (absl::GetFlag(FLAGS_debug)) spdlog::set_level(spdlog::level::debug);

  const int K = 17;
  const int n_workers = absl::GetFlag(FLAGS_workers);
  const bool canonical = absl::GetFlag(FLAGS_canonical);
  using KeyType = uint16_t;

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

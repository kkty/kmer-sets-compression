#include <algorithm>
#include <cstdlib>
#include <string>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/status/statusor.h"
#include "core/kmer_counter.h"
#include "core/kmer_set.h"
#include "core/kmer_set_compact.h"
#include "flags.h"
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
ABSL_FLAG(bool, check, false,
          "does compression & decompression to see if it is working");
ABSL_FLAG(std::string, out, "", "output file name");

template <int K, typename KeyType>
void Main(const std::string& file_name) {
  InitDefaultLogger();

  const bool debug = absl::GetFlag(FLAGS_debug);
  const std::string decompressor = absl::GetFlag(FLAGS_decompressor);
  const std::string compressor = absl::GetFlag(FLAGS_compressor);
  const bool canonical = absl::GetFlag(FLAGS_canonical);
  const int cutoff = absl::GetFlag(FLAGS_cutoff);
  const int n_workers = absl::GetFlag(FLAGS_workers);
  const bool check = absl::GetFlag(FLAGS_check);
  const std::string out = absl::GetFlag(FLAGS_out);

  if (debug) EnableDebugLogs();

  KmerSet<K, KeyType> kmer_set;

  {
    KmerCounter<K, KeyType> kmer_counter;

    {
      spdlog::info("constructing kmer_counter");

      absl::StatusOr<KmerCounter<K, KeyType>> statusor =
          KmerCounter<K, KeyType>::FromFASTA(file_name, decompressor, canonical,
                                             n_workers);

      if (!statusor.ok()) {
        spdlog::error("failed to parse FASTA file: {}",
                      statusor.status().ToString());
        std::exit(1);
      }

      kmer_counter = std::move(statusor).value();

      spdlog::info("constructed kmer_counter");
    }

    spdlog::info("constructing kmer_set");

    int64_t cutoff_count;
    std::tie(kmer_set, cutoff_count) =
        kmer_counter.ToKmerSet(cutoff, n_workers);

    spdlog::info("constructed kmer_set");
    spdlog::info("cutoff_count = {}", cutoff_count);
    spdlog::info("kmer_set.Size() = {}", kmer_set.Size());
    spdlog::info("kmer_set.Hash() = {}", kmer_set.Hash(n_workers));
  }

  spdlog::info("constructing kmer_set_compact");
  const KmerSetCompact<K, KeyType> kmer_set_compact =
      KmerSetCompact<K, KeyType>::FromKmerSet(kmer_set, canonical, n_workers);
  spdlog::info("constructed kmer_set_compact");

  spdlog::info("kmer_set_compact.Size() = {}",
               kmer_set_compact.Size(n_workers));

  if (check) {
    const KmerSet<K, KeyType> decompressed =
        kmer_set_compact.ToKmerSet(canonical, n_workers);
    if (!kmer_set.Equals(decompressed, n_workers)) {
      spdlog::error("decompressed is not equal to kmer_set");
      std::exit(1);
    }
  }

  if (!out.empty()) {
    absl::Status status = kmer_set_compact.Dump(out, compressor);

    if (!status.ok()) {
      spdlog::error("failed to dump kmer_set_compact: {}", status.ToString());
      std::exit(1);
    }
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

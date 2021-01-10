#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <utility>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/flags/usage.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_format.h"
#include "core/kmer_counter.h"
#include "core/kmer_set.h"
#include "core/kmer_set_compact.h"
#include "flags.h"
#include "log.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(int, k, 15, GetFlagMessage("k"));
ABSL_FLAG(bool, debug, false, GetFlagMessage("debug"));
ABSL_FLAG(std::string, decompressor, "", GetFlagMessage("decompressor"));
ABSL_FLAG(std::string, compressor, "", GetFlagMessage("compressor"));
ABSL_FLAG(int, workers, 1, GetFlagMessage("workers"));
ABSL_FLAG(bool, canonical, true, GetFlagMessage("canonical"));

ABSL_FLAG(int, cutoff, 1,
          "ignore k-mers that appear less often than this value");
ABSL_FLAG(bool, check, false,
          "does compression & decompression to see if it is working correctly");
ABSL_FLAG(std::string, out, "", "output file name");

template <int K, int N, typename KeyType>
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

  KmerSet<K, N, KeyType> kmer_set;

  {
    KmerCounter<K, N, KeyType> kmer_counter;

    {
      spdlog::info("constructing kmer_counter");

      absl::StatusOr<KmerCounter<K, N, KeyType>> statusor =
          KmerCounter<K, N, KeyType>::FromFASTA(file_name, decompressor,
                                                canonical, n_workers);

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
  const KmerSetCompact<K, N, KeyType> kmer_set_compact =
      KmerSetCompact<K, N, KeyType>::FromKmerSet(kmer_set, canonical, true,
                                                 n_workers);
  spdlog::info("constructed kmer_set_compact");

  spdlog::info("kmer_set_compact.Size() = {}",
               kmer_set_compact.Size(n_workers));

  if (check) {
    const KmerSet<K, N, KeyType> decompressed =
        kmer_set_compact.ToKmerSet(canonical, n_workers);

    if (kmer_set.Equals(decompressed, n_workers)) {
      spdlog::info("kmer_set_compact -> KmerSet: ok");
    } else {
      spdlog::error("kmer_set_compact -> KmerSet: failed");
      std::exit(1);
    }
  }

  if (!out.empty()) {
    absl::Status status = kmer_set_compact.Dump(out, compressor, n_workers);

    if (!status.ok()) {
      spdlog::error("failed to dump kmer_set_compact: {}", status.ToString());
      std::exit(1);
    }
  }
}

int main(int argc, char** argv) {
  absl::SetProgramUsageMessage(
      absl::StrFormat("Reads a FASTA file and constructs a set of k-mers. "
                      "Usage: %s [options] <path to file>",
                      argv[0]));

  const std::vector<std::string> args = ParseFlags(argc, argv);

  if (args.size() != 1) {
    spdlog::error("file name should be provided");
    std::exit(1);
  }

  const std::string file_name = args.front();

  const int k = absl::GetFlag(FLAGS_k);

  switch (k) {
    case 15:
      Main<15, 14, std::uint16_t>(file_name);
      break;
    case 19:
      Main<19, 10, std::uint32_t>(file_name);
      break;
    case 23:
      Main<23, 14, std::uint32_t>(file_name);
      break;
    default:
      spdlog::error("unsupported k value");
      std::exit(1);
  }
}

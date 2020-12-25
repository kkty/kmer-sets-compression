#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/usage.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_format.h"
#include "core/kmer_set.h"
#include "core/kmer_set_compact.h"
#include "flags.h"
#include "log.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(int, k, 15, GetFlagMessage("k"));
ABSL_FLAG(bool, debug, false, GetFlagMessage("debug"));
ABSL_FLAG(std::string, decompressor, "", GetFlagMessage("decompressor"));
ABSL_FLAG(int, workers, 1, GetFlagMessage("workers"));
ABSL_FLAG(bool, canonical, true, GetFlagMessage("canonical"));

template <int K, int N, typename KeyType>
void Main(const std::vector<std::string>& files) {
  InitDefaultLogger();

  if (absl::GetFlag(FLAGS_debug)) EnableDebugLogs();

  const int n_workers = absl::GetFlag(FLAGS_workers);
  const bool canonical = absl::GetFlag(FLAGS_canonical);
  const std::string decompressor = absl::GetFlag(FLAGS_decompressor);

  for (size_t i = 0; i < files.size(); i++) {
    const std::string& file_name = files[i];
    spdlog::info("processing: i = {}, file_name = {}", i, file_name);

    KmerSetCompact<K, N, KeyType> kmer_set_compact;

    {
      spdlog::info("constructing kmer_set_compact");

      absl::StatusOr<KmerSetCompact<K, N, KeyType>> statusor =
          KmerSetCompact<K, N, KeyType>::Load(file_name, decompressor);

      if (!statusor.ok()) {
        spdlog::error("failed to load kmer_set_compact: {}",
                      statusor.status().ToString());
        std::exit(1);
      }

      kmer_set_compact = std::move(statusor).value();

      spdlog::info("constructed kmer_set_compact");
    }

    spdlog::info("constructing kmer_set");
    const KmerSet<K, N, KeyType> kmer_set =
        kmer_set_compact.ToKmerSet(canonical, n_workers);
    spdlog::info("constructed kmer_set");

    int64_t size = kmer_set.Size();
    size_t hash = kmer_set.Hash(n_workers);

    spdlog::info("size = {}", size);
    spdlog::info("hash = {}", hash);

    std::cout << i << '\t' << file_name << '\t' << size << '\t' << hash
              << std::endl;

    spdlog::info("processed: i = {}, file_name = {}", i, file_name);
  }
}

int main(int argc, char** argv) {
  absl::SetProgramUsageMessage(
      absl::StrFormat("Prints the metadata of a k-mer set. Usage: %s [options] "
                      "<path to file>",
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

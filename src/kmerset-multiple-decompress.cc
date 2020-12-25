#include <cstdint>
#include <cstdlib>
#include <utility>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/usage.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_format.h"
#include "core/kmer_set.h"
#include "core/kmer_set_set.h"
#include "flags.h"
#include "log.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(int, k, 15, GetFlagMessage("k"));
ABSL_FLAG(bool, debug, false, GetFlagMessage("debug"));
ABSL_FLAG(std::string, decompressor, "", GetFlagMessage("decompressor"));
ABSL_FLAG(int, workers, 1, GetFlagMessage("workers"));
ABSL_FLAG(bool, canonical, true, GetFlagMessage("canonical"));

ABSL_FLAG(std::string, extension, "txt", "extension of files in folder");

template <int K, int N, typename KeyType>
void Main(const std::string& directory_name) {
  InitDefaultLogger();

  if (absl::GetFlag(FLAGS_debug)) EnableDebugLogs();

  const int n_workers = absl::GetFlag(FLAGS_workers);
  const bool canonical = absl::GetFlag(FLAGS_canonical);
  const std::string decompressor = absl::GetFlag(FLAGS_decompressor);
  const std::string extension = absl::GetFlag(FLAGS_extension);

  KmerSetSetReader<K, N, KeyType> kmer_set_set_reader;

  spdlog::info("loading kmer_set_set_reader");

  {
    absl::StatusOr<KmerSetSetReader<K, N, KeyType>> statusor =
        KmerSetSetReader<K, N, KeyType>::FromDirectory(
            directory_name, extension, decompressor, canonical);

    if (!statusor.ok()) {
      spdlog::error("failed to load data: {}", statusor.status().ToString());
      std::exit(1);
    }

    kmer_set_set_reader = std::move(statusor).value();
  }

  spdlog::info("loaded kmer_set_set_reader");

  spdlog::info("kmer_set_set_reader.Size() = {}", kmer_set_set_reader.Size());

  for (int i = 0; i < kmer_set_set_reader.Size(); i++) {
    KmerSet<K, N, KeyType> kmer_set;

    spdlog::info("constructing kmer_set: i = {}", i);

    {
      absl::StatusOr<KmerSet<K, N, KeyType>> statusor =
          kmer_set_set_reader.Get(i, n_workers);

      if (!statusor.ok()) {
        spdlog::error("failed to construct kmer_set: {}",
                      statusor.status().ToString());
        std::exit(1);
      }

      kmer_set = std::move(statusor).value();
    }

    spdlog::info("constructed kmer_set: i = {}", i);

    spdlog::info("kmer_set.Hash() = {}", kmer_set.Hash(n_workers));
    spdlog::info("kmer_set.Size() = {}", kmer_set.Size());
  }
}

int main(int argc, char** argv) {
  absl::SetProgramUsageMessage(
      absl::StrFormat("Decompresses the output of \"kmerset-compress\". Usage: "
                      "%s [options] <path to directory>",
                      argv[0]));

  std::string directory_name;

  {
    const std::vector<std::string> args = ParseFlags(argc, argv);

    if (args.size() != 1) {
      spdlog::error("invalid argument");
      std::exit(1);
    }

    directory_name = args.front();
  }

  const int k = absl::GetFlag(FLAGS_k);

  switch (k) {
    case 15:
      Main<15, 14, std::uint16_t>(directory_name);
      break;
    case 19:
      Main<19, 10, std::uint32_t>(directory_name);
      break;
    case 23:
      Main<23, 14, std::uint32_t>(directory_name);
      break;
    default:
      spdlog::error("unsupported k value");
      std::exit(1);
  }
}

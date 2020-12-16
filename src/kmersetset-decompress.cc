#include <algorithm>
#include <cstdlib>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/status/statusor.h"
#include "core/kmer_set.h"
#include "core/kmer_set_set.h"
#include "flags.h"
#include "log.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(int, k, 15, "the length of kmers");
ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(std::string, decompressor, "",
          "specify decompressor for input files");
ABSL_FLAG(bool, canonical, false, "count canonical k-mers");
ABSL_FLAG(int, workers, 1, "number of workers");
ABSL_FLAG(bool, folder, false, "load from files in a folder");
ABSL_FLAG(std::string, extension, "txt", "extension of files in folder");

template <int K, int N, typename KeyType>
void Main(const std::string& path) {
  InitDefaultLogger();

  if (absl::GetFlag(FLAGS_debug)) EnableDebugLogs();

  const int n_workers = absl::GetFlag(FLAGS_workers);
  const bool canonical = absl::GetFlag(FLAGS_canonical);
  const std::string decompressor = absl::GetFlag(FLAGS_decompressor);
  const bool folder = absl::GetFlag(FLAGS_folder);

  if (folder) {
    const std::string extension = absl::GetFlag(FLAGS_extension);

    KmerSetSetReader<K, N, KeyType> kmer_set_set_reader;

    spdlog::info("loading kmer_set_set_reader");

    {
      absl::StatusOr<KmerSetSetReader<K, N, KeyType>> statusor =
          KmerSetSetReader<K, N, KeyType>::FromFolder(path, extension,
                                                      decompressor, canonical);

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
  } else {
    KmerSetSet<K, N, KeyType> kmer_set_set;

    {
      spdlog::info("loading kmer_set_set");

      absl::StatusOr<KmerSetSet<K, N, KeyType>> statusor =
          KmerSetSet<K, N, KeyType>::Load(path, decompressor, canonical,
                                          n_workers);

      if (!statusor.ok()) {
        spdlog::error("failed to load data: {}", statusor.status().ToString());
        std::exit(1);
      }

      kmer_set_set = std::move(statusor).value();

      spdlog::info("loaded kmer_set_set");
    }

    for (int i = 0; i < kmer_set_set.Size(); i++) {
      spdlog::info("i = {}", i);
      const KmerSet<K, N, KeyType> kmer_set = kmer_set_set.Get(i, n_workers);
      spdlog::info("kmer_set.Size() = {}", kmer_set.Size());
      spdlog::info("kmer_set.Hash() = {}", kmer_set.Hash(n_workers));
    }
  }
}

int main(int argc, char** argv) {
  std::string path;

  {
    const std::vector<std::string> args = ParseFlags(argc, argv);

    if (args.size() != 1) {
      spdlog::error("invalid argument");
      std::exit(1);
    }

    path = args.front();
  }

  const int k = absl::GetFlag(FLAGS_k);

  switch (k) {
    case 15:
      Main<15, 14, uint16_t>(path);
      break;
    case 19:
      Main<19, 10, uint32_t>(path);
      break;
    case 23:
      Main<23, 14, uint32_t>(path);
      break;
    default:
      spdlog::error("unsupported k value");
      std::exit(1);
  }
}

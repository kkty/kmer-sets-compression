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

template <int K, typename KeyType>
void Main(const std::string& file) {
  InitDefaultLogger();

  if (absl::GetFlag(FLAGS_debug)) EnableDebugLogs();

  const int n_workers = absl::GetFlag(FLAGS_workers);
  const bool canonical = absl::GetFlag(FLAGS_canonical);
  const std::string decompressor = absl::GetFlag(FLAGS_decompressor);

  KmerSetSet<K, KeyType> kmer_set_set;

  {
    absl::StatusOr<KmerSetSet<K, KeyType>> statusor =
        KmerSetSet<K, KeyType>::Load(file, decompressor, canonical, n_workers);

    if (!statusor.ok()) {
      spdlog::error("failed to load data: {}", statusor.status().ToString());
      std::exit(1);
    }

    kmer_set_set = std::move(statusor).value();
  }

  for (int i = 0; i < kmer_set_set.Size(); i++) {
    spdlog::info("i = {}", i);
    const KmerSet<K, KeyType> kmer_set = kmer_set_set.Get(i, n_workers);
    spdlog::info("kmer_set.Size() = {}", kmer_set.Size());
    spdlog::info("kmer_set.Hash() = {}", kmer_set.Hash(n_workers));
  }
}

int main(int argc, char** argv) {
  std::string file;

  {
    const std::vector<std::string> args = ParseFlags(argc, argv);

    if (args.size() != 1) {
      spdlog::error("invalid argument");
      std::exit(1);
    }

    file = args.front();
  }

  const int k = absl::GetFlag(FLAGS_k);

  switch (k) {
    case 15:
      // 15 * 2 - 16 = 14 bits are used to select buckets.
      Main<15, uint16_t>(file);
      break;
    case 19:
      // 19 * 2 - 16 = 22 bits are used to select buckets.
      Main<19, uint16_t>(file);
      break;
    case 23:
      // 23 * 2 - 32 = 14 bits are used to select buckets.
      Main<23, uint32_t>(file);
      break;
    case 27:
      // 27 * 2 - 32 = 22 bits are used to select buckets.
      Main<27, uint32_t>(file);
      break;
    default:
      spdlog::error("unsupported k value");
      std::exit(1);
  }
}

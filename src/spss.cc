#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "core/kmer.h"
#include "core/kmer_set.h"
#include "core/kmer_set_compact.h"
#include "core/unitigs.h"
#include "log.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(int, k, 15, "the length of kmers");
ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(std::string, decompressor, "cat", "specify decompressor");
ABSL_FLAG(int, workers, 1, "number of workers");
ABSL_FLAG(int, buckets, 1, "number of buckets for SPSSS calculation");

template <int K, typename KeyType>
void Main(const std::string& file_name) {
  InitDefaultLogger();

  if (absl::GetFlag(FLAGS_debug)) EnableDebugLogs();

  const int n_workers = absl::GetFlag(FLAGS_workers);
  const int n_buckets = absl::GetFlag(FLAGS_buckets);

  spdlog::info("constructing kmer_set_compact");
  const KmerSetCompact<K, KeyType> kmer_set_compact =
      KmerSetCompact<K, KeyType>::Load(file_name,
                                       absl::GetFlag(FLAGS_decompressor));
  spdlog::info("constructed kmer_set_compact");

  spdlog::info("constructing kmer_set");
  const KmerSet<K, KeyType> kmer_set =
      kmer_set_compact.ToKmerSet(true, n_workers);
  spdlog::info("constructed kmer_set");

  spdlog::info("constructing spss");
  std::vector<std::string> spss =
      GetSPSSCanonical(kmer_set, n_workers, n_buckets);
  spdlog::info("constructed spss");

  int64_t total_size = 0;
  for (const std::string& s : spss) total_size += s.length();
  spdlog::info("total_size = {}", total_size);
}

int main(int argc, char** argv) {
  const std::string file_name = absl::ParseCommandLine(argc, argv)[1];

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
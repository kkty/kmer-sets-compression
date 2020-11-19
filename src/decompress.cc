#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "core/kmer.h"
#include "core/kmer_set.h"
#include "core/kmer_set_compressed.h"
#include "log.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(int, k, 15, "the length of kmers");
ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(std::string, decompressor, "cat", "specify decompressor");
ABSL_FLAG(int, workers, 1, "number of workers");
ABSL_FLAG(bool, canonical, false, "use canonical k-mers");

template <int K, typename KeyType>
void Main(const std::string& file_name) {
  InitDefaultLogger();

  if (absl::GetFlag(FLAGS_debug)) EnableDebugLogs();

  const int n_workers = absl::GetFlag(FLAGS_workers);

  spdlog::info("constructing kmer_set_compressed");
  const KmerSetCompressed<K, KeyType> kmer_set_compressed =
      KmerSetCompressed<K, KeyType>::Load(file_name,
                                          absl::GetFlag(FLAGS_decompressor));
  spdlog::info("constructed kmer_set_compressed");

  spdlog::info("constructing kmer_set");
  const KmerSet<K, KeyType> kmer_set =
      kmer_set_compressed.ToKmerSet(absl::GetFlag(FLAGS_canonical), n_workers);
  spdlog::info("constructed kmer_set");

  spdlog::info("kmer_set.Size() = {}", kmer_set.Size());
  spdlog::info("kmer_set.Hash() = {}", kmer_set.Hash(n_workers));
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

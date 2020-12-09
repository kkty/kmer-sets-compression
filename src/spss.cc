#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "core/kmer_set.h"
#include "core/kmer_set_compact.h"
#include "core/unitigs.h"
#include "log.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(int, k, 15, "the length of kmers");
ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(std::string, decompressor, "", "specify decompressor");
ABSL_FLAG(int, workers, 1, "number of workers");
ABSL_FLAG(int, buckets, 1, "number of buckets for SPSS calculation");

template <int K, int N, typename KeyType>
void Main(const std::string& file_name) {
  InitDefaultLogger();

  if (absl::GetFlag(FLAGS_debug)) EnableDebugLogs();

  const int n_workers = absl::GetFlag(FLAGS_workers);
  const int n_buckets = absl::GetFlag(FLAGS_buckets);

  KmerSetCompact<K, N, KeyType> kmer_set_compact;

  {
    spdlog::info("constructing kmer_set_compact");

    const std::string decompressor = absl::GetFlag(FLAGS_decompressor);

    absl::StatusOr<KmerSetCompact<K, N, KeyType>> statusor =
        KmerSetCompact<K, N, KeyType>::Load(file_name, decompressor);

    if (!statusor.ok()) {
      spdlog::error("failed to load kmer_set_compact from file: {}",
                    statusor.status().ToString());
      std::exit(1);
    }

    kmer_set_compact = std::move(statusor).value();

    spdlog::info("constructed kmer_set_compact");
  }

  spdlog::info("constructing kmer_set");
  const KmerSet<K, N, KeyType> kmer_set =
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
      Main<15, 14, uint16_t>(file_name);
      break;
    case 19:
      Main<19, 10, uint32_t>(file_name);
      break;
    case 23:
      Main<23, 14, uint32_t>(file_name);
      break;
    default:
      spdlog::error("unsupported k value");
      std::exit(1);
  }
}

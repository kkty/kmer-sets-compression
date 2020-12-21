#include "core/spss.h"

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "core/kmer_set.h"
#include "io.h"
#include "log.h"
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

ABSL_FLAG(int, k, 15, "the length of kmers");
ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(std::string, decompressor, "", "specify decompressor");
ABSL_FLAG(int, workers, 1, "number of workers");
ABSL_FLAG(int, buckets, 1, "number of buckets for SPSS calculation");
ABSL_FLAG(int, repests, 1, "number of repeats");

template <int K, int N, typename KeyType>
void Main(const std::string& file_name) {
  InitDefaultLogger();

  if (absl::GetFlag(FLAGS_debug)) EnableDebugLogs();

  const int n_workers = absl::GetFlag(FLAGS_workers);
  const int n_buckets = absl::GetFlag(FLAGS_buckets);
  const std::string decompressor = absl::GetFlag(FLAGS_decompressor);
  const int n_repeats = absl::GetFlag(FLAGS_repests);

  KmerSet<K, N, KeyType> kmer_set;

  {
    absl::StatusOr<KmerSet<K, N, KeyType>> statusor =
        GetKmerSetFromFile<K, N, KeyType>(file_name, decompressor, true,
                                          n_workers);

    kmer_set = std::move(statusor).value();
  }

  spdlog::info("kmer_set.Size() = {}", kmer_set.Size());
  spdlog::info("kmer_set.Hash() = {}", kmer_set.Hash(n_workers));

  for (int i = 0; i < n_repeats; i++) {
    for (bool fast : std::vector<bool>{false, true}) {
      spdlog::info("fast = {}", fast);

      std::vector<std::string> spss;

      {
        spdlog::info("constructing spss");
        spdlog::stopwatch sw;

        spss = GetSPSSCanonical(kmer_set, fast, n_workers, n_buckets);

        spdlog::info("constructed spss: elapsed = {}", sw);

        std::cout << sw.elapsed().count() << ' ';
      }

      {
        std::int64_t total_size = 0;
        for (const std::string& s : spss) total_size += s.length();

        spdlog::info("total_size = {}", total_size);

        std::cout << total_size << ' ';
      }

      {
        spdlog::info("reconstructing");
        spdlog::stopwatch sw;

        const KmerSet<K, N, KeyType> reconstructed =
            GetKmerSetFromSPSS<K, N, KeyType>(spss, true, n_workers);

        spdlog::info("reconstructed: elapsed = {}", sw);

        std::cout << sw.elapsed().count() << ' ';

        const bool is_equal = kmer_set.Equals(reconstructed, n_workers);

        spdlog::info("is_equal = {}", is_equal);

        std::cout << is_equal << ' ';
      }
    }

    std::cout << std::endl;
  }
}

int main(int argc, char** argv) {
  const std::string file_name = absl::ParseCommandLine(argc, argv)[1];

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

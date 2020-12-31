#include "core/spss.h"

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/flags/usage.h"
#include "absl/strings/str_format.h"
#include "core/kmer.h"
#include "core/kmer_set.h"
#include "flags.h"
#include "io.h"
#include "log.h"
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

ABSL_FLAG(int, k, 15, GetFlagMessage("k"));
ABSL_FLAG(bool, debug, false, GetFlagMessage("debug"));
ABSL_FLAG(std::string, decompressor, "", GetFlagMessage("decompressor"));
ABSL_FLAG(int, workers, 1, GetFlagMessage("workers"));

ABSL_FLAG(int, buckets, 1, "number of buckets for SPSS calculation");
ABSL_FLAG(int, repeats, 1, "number of repeats");

template <int K, int N, typename KeyType>
void Main(const std::string& file_name) {
  InitDefaultLogger();

  if (absl::GetFlag(FLAGS_debug)) EnableDebugLogs();

  const int n_workers = absl::GetFlag(FLAGS_workers);
  const int n_buckets = absl::GetFlag(FLAGS_buckets);
  const std::string decompressor = absl::GetFlag(FLAGS_decompressor);
  const int n_repeats = absl::GetFlag(FLAGS_repeats);

  KmerSet<K, N, KeyType> kmer_set;

  {
    absl::StatusOr<KmerSet<K, N, KeyType>> statusor =
        GetKmerSetFromFile<K, N, KeyType>(file_name, decompressor, true,
                                          n_workers);

    kmer_set = std::move(statusor).value();
  }

  spdlog::info("kmer_set.Size() = {}", kmer_set.Size());
  spdlog::info("kmer_set.Hash() = {}", kmer_set.Hash(n_workers));

  std::int64_t weight = GetSPSSWeightCanonical(kmer_set, n_workers);

  spdlog::info("weight = {}", weight);

  spdlog::info("constructing unitigs");

  const std::vector<std::string> unitigs =
      GetUnitigsCanonical(kmer_set, n_workers);

  spdlog::info("constructed unitigs");

  spdlog::info("constructing prefixes and suffixes");

  absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>> prefixes =
      GetPrefixesFromUnitigs<K>(unitigs, n_workers);

  absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>> suffixes =
      GetSuffixesFromUnitigs<K>(unitigs, n_workers);

  spdlog::info("constructed prefixes and suffixes");

  for (int i = 0; i < n_repeats; i++) {
    for (bool fast : std::vector<bool>{false, true}) {
      spdlog::info("fast = {}", fast);

      std::vector<std::string> spss;

      {
        spdlog::info("constructing spss");
        spdlog::stopwatch sw;

        spss = GetSPSSCanonical<K, N, KeyType>(unitigs, prefixes, suffixes,
                                               fast, n_workers, n_buckets);

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
  absl::SetProgramUsageMessage(
      absl::StrFormat("Runs a benchmark for SPSS construction using a single "
                      "k-mer set. Usage: %s [options] <path to file>",
                      argv[0]));

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

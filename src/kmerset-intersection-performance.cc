#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/usage.h"
#include "absl/strings/str_format.h"
#include "core/kmer_set.h"
#include "core/kmer_set_immutable.h"
#include "flags.h"
#include "io.h"
#include "log.h"
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

ABSL_FLAG(int, k, 15, GetFlagMessage("k"));
ABSL_FLAG(bool, debug, false, GetFlagMessage("debug"));
ABSL_FLAG(int, workers, 1, GetFlagMessage("workers"));
ABSL_FLAG(std::string, decompressor, "", GetFlagMessage("decompressor"));
ABSL_FLAG(bool, canonical, true, GetFlagMessage("canonical"));

ABSL_FLAG(int, repeats, 100, "number of repeats");

template <int K, int N, typename KeyType>
void Main(const std::string& file1, const std::string& file2) {
  InitDefaultLogger();

  const bool debug = absl::GetFlag(FLAGS_debug);
  const int n_workers = absl::GetFlag(FLAGS_workers);
  const int n_repeats = absl::GetFlag(FLAGS_repeats);
  const bool canonical = absl::GetFlag(FLAGS_canonical);
  const std::string decompressor = absl::GetFlag(FLAGS_decompressor);

  if (debug) EnableDebugLogs();

  std::vector<KmerSetImmutable<K, N, KeyType>> kmer_sets_immutable;
  std::vector<KmerSetImmutableHashSet<K>> kmer_sets_immutable_hash_set;
  std::vector<KmerSetImmutableVector<K>> kmer_sets_immutable_vector;

  const auto ReadFile = [&](const std::string& file) {
    absl::StatusOr<KmerSet<K, N, KeyType>> statusor =
        GetKmerSetFromFile<K, N, KeyType>(file, decompressor, canonical,
                                          n_workers);

    if (!statusor.ok()) {
      spdlog::error("failed to read file: {}", statusor.status().ToString());
      std::exit(1);
    }

    KmerSet<K, N, KeyType> kmer_set = std::move(statusor).value();

    KmerSetImmutable<K, N, KeyType> kmer_set_immutable(kmer_set, n_workers);

    spdlog::info("kmer_set_immutable.Size() = {}", kmer_set_immutable.Size());
    spdlog::info("kmer_set_immutable.Bytes() = {}", kmer_set_immutable.Bytes());

    kmer_sets_immutable.push_back(std::move(kmer_set_immutable));

    kmer_sets_immutable_hash_set.push_back(
        KmerSetImmutableHashSet<K>::FromKmerSet(kmer_set, n_workers));

    kmer_sets_immutable_vector.push_back(
        KmerSetImmutableVector<K>::FromKmerSet(kmer_set, n_workers));
  };

  ReadFile(file1);
  ReadFile(file2);

  for (int i = 0; i < n_repeats; i++) {
    {
      spdlog::info("intersecting kmer_sets_immutable");

      spdlog::stopwatch sw;

      std::int64_t size = kmer_sets_immutable[0]
                              .Intersection(kmer_sets_immutable[1], n_workers)
                              .Size();

      spdlog::info("size = {}, elapsed = {}", size, sw);

      std::cout << sw.elapsed().count() << ' ';
    }

    {
      spdlog::info("intersecting kmer_sets_immutable_hash_set");

      spdlog::stopwatch sw;

      std::int64_t size = kmer_sets_immutable_hash_set[0]
                              .Intersection(kmer_sets_immutable_hash_set[1])
                              .Size();

      spdlog::info("size = {}, elapsed = {}", size, sw);

      std::cout << sw.elapsed().count() << ' ';
    }

    {
      spdlog::info("intersecting kmer_sets_immutable_vector");

      spdlog::stopwatch sw;

      std::int64_t size = kmer_sets_immutable_vector[0]
                              .Intersection(kmer_sets_immutable_vector[1])
                              .Size();

      spdlog::info("size = {}, elapsed = {}", size, sw);

      std::cout << sw.elapsed().count() << ' ';
    }

    {
      spdlog::info("adding kmer_sets_immutable");

      spdlog::stopwatch sw;

      std::int64_t size =
          kmer_sets_immutable[0].Add(kmer_sets_immutable[1], n_workers).Size();

      spdlog::info("size = {}, elapsed = {}", size, sw);

      std::cout << sw.elapsed().count() << ' ';
    }

    {
      spdlog::info("adding kmer_sets_immutable_hash_set");

      spdlog::stopwatch sw;

      std::int64_t size = kmer_sets_immutable_hash_set[0]
                              .Add(kmer_sets_immutable_hash_set[1])
                              .Size();

      spdlog::info("size = {}, elapsed = {}", size, sw);

      std::cout << sw.elapsed().count() << ' ';
    }

    {
      spdlog::info("adding kmer_sets_immutable_vector");

      spdlog::stopwatch sw;

      std::int64_t size = kmer_sets_immutable_vector[0]
                              .Add(kmer_sets_immutable_vector[1])
                              .Size();

      spdlog::info("size = {}, elapsed = {}", size, sw);

      std::cout << sw.elapsed().count() << ' ';
    }

    {
      spdlog::info("subtracting kmer_sets_immutable");

      spdlog::stopwatch sw;

      std::int64_t size =
          kmer_sets_immutable[0].Sub(kmer_sets_immutable[1], n_workers).Size();

      spdlog::info("size = {}, elapsed = {}", size, sw);

      std::cout << sw.elapsed().count() << ' ';
    }

    {
      spdlog::info("subtracting kmer_sets_immutable_hash_set");

      spdlog::stopwatch sw;

      std::int64_t size = kmer_sets_immutable_hash_set[0]
                              .Sub(kmer_sets_immutable_hash_set[1])
                              .Size();

      spdlog::info("size = {}, elapsed = {}", size, sw);

      std::cout << sw.elapsed().count() << ' ';
    }

    {
      spdlog::info("subtracting kmer_sets_immutable_vector");

      spdlog::stopwatch sw;

      std::int64_t size = kmer_sets_immutable_vector[0]
                              .Sub(kmer_sets_immutable_vector[1])
                              .Size();

      spdlog::info("size = {}, elapsed = {}", size, sw);

      std::cout << sw.elapsed().count() << ' ';
    }

    std::cout << std::endl;
  }
}

int main(int argc, char** argv) {
  absl::SetProgramUsageMessage(
      absl::StrFormat("Runs a benchmark for intersecting 2 k-mer sets. Usage: "
                      "%s [options] <path to file> <path to file>",
                      argv[0]));

  const std::vector<std::string> files = ParseFlags(argc, argv);

  const int k = absl::GetFlag(FLAGS_k);

  switch (k) {
    case 15:
      Main<15, 14, std::uint16_t>(files[0], files[1]);
      break;
    case 19:
      Main<19, 10, std::uint32_t>(files[0], files[1]);
      break;
    case 23:
      Main<23, 14, std::uint32_t>(files[0], files[1]);
      break;
    default:
      spdlog::error("unsupported k value");
      std::exit(1);
  }
}

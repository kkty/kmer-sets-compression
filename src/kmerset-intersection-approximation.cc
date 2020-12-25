#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
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

ABSL_FLAG(int, k, 15, GetFlagMessage("k"));
ABSL_FLAG(bool, debug, false, GetFlagMessage("debug"));
ABSL_FLAG(int, workers, 1, GetFlagMessage("workers"));
ABSL_FLAG(std::string, decompressor, "", GetFlagMessage("decompressor"));
ABSL_FLAG(bool, canonical, true, GetFlagMessage("canonical"));

ABSL_FLAG(double, factor, 0.1, "factor to decide the number of buckets to use");
ABSL_FLAG(int, repeats, 100, "number of repeats");

template <int K, int N, typename KeyType>
void Main(const std::string& file1, const std::string& file2) {
  InitDefaultLogger();

  const bool debug = absl::GetFlag(FLAGS_debug);
  const int n_workers = absl::GetFlag(FLAGS_workers);
  const int n_repeats = absl::GetFlag(FLAGS_repeats);
  const bool canonical = absl::GetFlag(FLAGS_canonical);
  const double factor = absl::GetFlag(FLAGS_factor);
  const std::string decompressor = absl::GetFlag(FLAGS_decompressor);

  if (debug) EnableDebugLogs();

  std::vector<KmerSetImmutable<K, N, KeyType>> kmer_sets_immutable;

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
  };

  ReadFile(file1);
  ReadFile(file2);

  const std::int64_t exact_size = kmer_sets_immutable[0].IntersectionSize(
      kmer_sets_immutable[1], n_workers);

  spdlog::info("exact_size = {}", exact_size);

  std::cout << std::setprecision(5);

  for (int i = 0; i < n_repeats; i++) {
    spdlog::info("intersecting kmer_sets_immutable");

    const std::int64_t size = kmer_sets_immutable[0].IntersectionSizeEstimate(
        kmer_sets_immutable[1], factor);

    const double error = static_cast<double>(size - exact_size) / exact_size;

    spdlog::info("size = {}, error = {}", size, error);

    std::cout << exact_size << ' ' << size << ' ' << error << std::endl;
  }
}

int main(int argc, char** argv) {
  absl::SetProgramUsageMessage(absl::StrFormat(
      "Evaluates the accuracy of intersection size estimation, "
      "using 2 k-mer sets. Usage: %s [options] <path to file> <path to file>",
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

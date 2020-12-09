#include <algorithm>
#include <cstdlib>
#include <string>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/status/statusor.h"
#include "core/kmer_set.h"
#include "core/kmer_set_immutable.h"
#include "flags.h"
#include "io.h"
#include "log.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(int, k, 15, "the length of kmers");
ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(std::string, decompressor, "",
          "specify decompressor for input files");
ABSL_FLAG(bool, canonical, false, "count canonical k-mers");
ABSL_FLAG(int, workers, 1, "number of workers");

template <int K, int N, typename KeyType>
void Main(const std::vector<std::string>& files) {
  InitDefaultLogger();

  const bool debug = absl::GetFlag(FLAGS_debug);
  const std::string decompressor = absl::GetFlag(FLAGS_decompressor);
  const bool canonical = absl::GetFlag(FLAGS_canonical);
  const int n_workers = absl::GetFlag(FLAGS_workers);

  if (debug) EnableDebugLogs();

  const int n_files = files.size();

  std::vector<KmerSetImmutable<K, N, KeyType>> kmer_sets_immutable(n_files);

  for (int i = 0; i < n_files; i++) {
    const std::string& file = files[i];

    spdlog::info("reading file: {}", file);

    absl::StatusOr<KmerSet<K, N, KeyType>> statusor =
        GetKmerSetFromCompressedKmersFile<K, N, KeyType>(file, decompressor,
                                                         canonical, n_workers);

    if (!statusor.ok()) {
      spdlog::error("failed to read file: {}", statusor.status().ToString());
      std::exit(1);
    }

    spdlog::info("finished reading file: {}", file);

    const KmerSet<K, N, KeyType> kmer_set = std::move(statusor).value();

    spdlog::info("kmer_set.Size() = {}", kmer_set.Size());
    spdlog::info("kmer_set.Hash() = {}", kmer_set.Hash(n_workers));

    KmerSetImmutable<K, N, KeyType> kmer_set_immutable(kmer_set, n_workers);

    spdlog::info("kmer_set_immutable.Size() = {}", kmer_set_immutable.Size());
    spdlog::info("kmer_set_immutable.Bytes() = {}", kmer_set_immutable.Bytes());

    kmer_sets_immutable[i] = std::move(kmer_set_immutable);
  }

  {
    int64_t total_size = 0;

    for (int i = 0; i < n_files; i++) {
      total_size += kmer_sets_immutable[i].Size();
    }

    spdlog::info("total_size = {}", total_size);
  }

  // Lists up all the k-mers that appear at least once.
  {
    KmerSetImmutable<K, N, KeyType> all;

    for (int i = 0; i < n_files; i++) {
      all = all.Add(kmer_sets_immutable[i], n_workers);
    }

    spdlog::info("all.Size() = {}", all.Size());
    spdlog::info("all.Bytes() = {}", all.Bytes());
  }

  for (int i = 0; i < n_files; i++) {
    for (int j = i + 1; j < n_files; j++) {
      spdlog::info("i = {}, j = {}", i, j);

      spdlog::info("constructing intersection");
      const KmerSetImmutable<K, N, KeyType> intersection =
          kmer_sets_immutable[i].Intersection(kmer_sets_immutable[j],
                                              n_workers);
      spdlog::info("constructed intersection");

      spdlog::info("constructing sub_i");
      const KmerSetImmutable<K, N, KeyType> sub_i =
          kmer_sets_immutable[i].Sub(intersection, n_workers);
      spdlog::info("constructed sub_i");

      spdlog::info("constructing sub_j");
      const KmerSetImmutable<K, N, KeyType> sub_j =
          kmer_sets_immutable[j].Sub(intersection, n_workers);
      spdlog::info("constructed sub_j");

      spdlog::info("intersection.Size() = {}", intersection.Size());
      spdlog::info("sub_i.Size() = {}", sub_i.Size());
      spdlog::info("sub_j.Size() = {}", sub_j.Size());

      const double similarity =
          static_cast<double>(intersection.Size()) /
          (sub_i.Size() + sub_j.Size() + intersection.Size());

      spdlog::info("similarity = {}", similarity);
    }
  }
}

int main(int argc, char** argv) {
  const std::vector<std::string> files = ParseFlags(argc, argv);

  const int k = absl::GetFlag(FLAGS_k);

  switch (k) {
    case 15:
      Main<15, 14, uint16_t>(files);
      break;
    case 19:
      Main<19, 10, uint32_t>(files);
      break;
    case 23:
      Main<23, 14, uint32_t>(files);
      break;
    default:
      spdlog::error("unsupported k value");
      std::exit(1);
  }
}

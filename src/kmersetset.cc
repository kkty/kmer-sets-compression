#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/status/statusor.h"
#include "flags.h"
#include "core/graph.h"
#include "core/kmer.h"
#include "core/kmer_counter.h"
#include "core/kmer_set.h"
#include "core/kmer_set_compressed.h"
#include "core/kmer_set_set.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(int, k, 15, "the length of kmers");
ABSL_FLAG(std::string, type, "fastq", "input file type (fastq or ckmers)");
ABSL_FLAG(bool, similarity, false, "calculate Jaccard similarity");
ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(std::string, decompressor, "",
          "specify decompressor for input files");
ABSL_FLAG(bool, canonical, false, "count canonical k-mers");
ABSL_FLAG(int, cutoff, 1, "cut off threshold");
ABSL_FLAG(int, workers, 1, "number of workers");
ABSL_FLAG(int, recursion, 1, "recursion limit for KmerSetSet");
ABSL_FLAG(bool, check, false, "check if k-mer sets can be reconstructed");
ABSL_FLAG(bool, nj, false, "use neighbor joining algorithm");
ABSL_FLAG(bool, mm, false, "use maximum weighted matching algorithm");
ABSL_FLAG(bool, approximate_matching, false,
          "use approximate algorithm for matching");
ABSL_FLAG(bool, approximate_weights, false,
          "use approximate weights for the matching algorithm");
ABSL_FLAG(bool, approximate_graph, false,
          "use an approximate (sparse) graph for the matching algorithm");

template <int K, typename KeyType>
KmerSet<K, KeyType> GetKmerSetFromCompressedKmersFile(
    const std::string& file_name, const std::string& decompressor,
    bool canonical, int n_workers) {
  spdlog::info("constructing kmer_set_compressed");

  const KmerSetCompressed<K, KeyType> kmer_set_compressed =
      KmerSetCompressed<K, KeyType>::Load(file_name, decompressor);

  spdlog::info("constructed kmer_set_compressed");

  spdlog::info("constructing kmer_set");

  const KmerSet<K, KeyType> kmer_set =
      kmer_set_compressed.ToKmerSet(canonical, n_workers);

  spdlog::info("constructed kmer_set");

  return kmer_set;
}

template <int K, typename KeyType>
KmerSet<K, KeyType> GetKmerSetFromFASTQFile(const std::string& file_name,
                                            const std::string& decompressor,
                                            bool canonical, int cutoff,
                                            int n_workers) {
  spdlog::info("constructing kmer_counter");

  const KmerCounter<K, KeyType> kmer_counter = [&] {
    absl::StatusOr<KmerCounter<K, KeyType>> statusor =
        decompressor != "" ? KmerCounter<K, KeyType>::FromFASTQ(
                                 file_name, decompressor, canonical, n_workers)
                           : KmerCounter<K, KeyType>::FromFASTQ(
                                 file_name, canonical, n_workers);

    if (!statusor.ok()) {
      spdlog::error("failed to parse FASTQ file: {}",
                    statusor.status().ToString());
      std::exit(1);
    }

    return *statusor;
  }();

  spdlog::info("constructed kmer_counter");

  spdlog::info("constructing kmer_set");

  KmerSet<K, KeyType> kmer_set;
  int64_t cutoff_count;
  std::tie(kmer_set, cutoff_count) = kmer_counter.ToSet(cutoff, n_workers);

  spdlog::info("constructed kmer_set");

  return kmer_set;
}

template <int K, typename KeyType>
void Main(const std::vector<std::string>& files) {
  spdlog::set_default_logger(spdlog::stderr_color_mt("default"));

  if (absl::GetFlag(FLAGS_debug)) spdlog::set_level(spdlog::level::debug);

  const int n_workers = absl::GetFlag(FLAGS_workers);
  const std::string decompressor = absl::GetFlag(FLAGS_decompressor);
  const bool canonical = absl::GetFlag(FLAGS_canonical);
  const int cutoff = absl::GetFlag(FLAGS_cutoff);

  const int n_datasets = files.size();

  std::vector<KmerSet<K, KeyType>> kmer_sets(n_datasets);

  for (int i = 0; i < n_datasets; i++) {
    const std::string& file = files[i];

    spdlog::info("file = {}", file);

    if (absl::GetFlag(FLAGS_type) == "fastq") {
      kmer_sets[i] = GetKmerSetFromFASTQFile<K, KeyType>(
          file, decompressor, canonical, cutoff, n_workers);
    } else {
      kmer_sets[i] = GetKmerSetFromCompressedKmersFile<K, KeyType>(
          file, decompressor, canonical, n_workers);
    }
  }

  for (int i = 0; i < n_datasets; i++) {
    int64_t kmer_set_size = kmer_sets[i].Size();
    spdlog::info("i = {}, kmer_set_size = {}", i, kmer_set_size);
  }

  if (absl::GetFlag(FLAGS_similarity)) {
    for (int i = 0; i < n_datasets; i++) {
      for (int j = i + 1; j < n_datasets; j++) {
        double similarity = kmer_sets[i].Similarity(kmer_sets[j], n_workers);
        spdlog::info("i = {}, j = {}, similarity = {}", i, j, similarity);
      }
    }
  }

  const auto cost_function = [](const KmerSet<K, KeyType>& lhs,
                                const KmerSet<K, KeyType>& rhs, int n_workers) {
    return lhs.Diff(rhs, n_workers);
  };

  {
    int64_t total_cost = 0;
    for (int i = 0; i < n_datasets; i++) {
      int64_t cost =
          cost_function(KmerSet<K, KeyType>(), kmer_sets[i], n_workers);
      spdlog::info("i = {}, cost = {}", i, cost);
      total_cost += cost;
    }
    spdlog::info("total_cost = {}", total_cost);
  }

  if (absl::GetFlag(FLAGS_nj)) {
    spdlog::info("constructing kmer_set_set_nj");

    KmerSetSetNJ<K, KeyType, decltype(cost_function)> kmer_set_set_nj(
        kmer_sets, absl::GetFlag(FLAGS_recursion), cost_function, n_workers);
    spdlog::info("constructed kmer_set_set_nj");

    spdlog::info("kmer_set_set_nj.Cost() = {}", kmer_set_set_nj.Cost());

    if (absl::GetFlag(FLAGS_check)) {
      for (int i = 0; i < n_datasets; i++) {
        spdlog::info(
            "kmer_sets[{}].Equals(kmer_set_set_nj.Get({})) = {}", i, i,
            kmer_sets[i].Equals(kmer_set_set_nj.Get(i, n_workers), n_workers));
      }
    }
  } else if (absl::GetFlag(FLAGS_mm)) {
    spdlog::info("constructing kmer_set_set_mm");

    KmerSetSetMM<K, KeyType, decltype(cost_function)> kmer_set_set_mm(
        kmer_sets, absl::GetFlag(FLAGS_recursion),
        absl::GetFlag(FLAGS_approximate_matching),
        absl::GetFlag(FLAGS_approximate_weights),
        absl::GetFlag(FLAGS_approximate_graph), cost_function, n_workers);

    spdlog::info("constructed kmer_set_set_mm");

    spdlog::info("kmer_set_set_mm.Size() = {}", kmer_set_set_mm.Size());
    spdlog::info("kmer_set_set_mm.Cost() = {}", kmer_set_set_mm.Cost());

    if (absl::GetFlag(FLAGS_check)) {
      spdlog::info("dumping");
      std::vector<std::string> dumped =
          kmer_set_set_mm.Dump(absl::GetFlag(FLAGS_canonical), n_workers);
      spdlog::info("dumped");

      spdlog::info("loading");
      std::unique_ptr<KmerSetSetMM<K, KeyType, decltype(cost_function)>>
          loaded = std::unique_ptr<
              KmerSetSetMM<K, KeyType, decltype(cost_function)>>(
              KmerSetSetMM<K, KeyType, decltype(cost_function)>::Load(
                  std::move(dumped), absl::GetFlag(FLAGS_canonical),
                  n_workers));
      spdlog::info("loaded");

      for (int i = 0; i < n_datasets; i++) {
        spdlog::info("kmer_sets[{}].Equals(loaded->Get({})) = {}", i, i,
                     kmer_sets[i].Equals(loaded->Get(i, n_workers), n_workers));
      }
    }
  } else {
    spdlog::info("constructing kmer_set_set");

    KmerSetSet<K, KeyType, decltype(cost_function)> kmer_set_set(
        kmer_sets, absl::GetFlag(FLAGS_recursion), cost_function, n_workers);
    spdlog::info("constructed kmer_set_set");

    spdlog::info("kmer_set_set.Size() = {}", kmer_set_set.Size());
    spdlog::info("kmer_set_set.Cost() = {}", kmer_set_set.Cost());

    if (absl::GetFlag(FLAGS_check)) {
      for (int i = 0; i < n_datasets; i++) {
        spdlog::info(
            "kmer_sets[{}].Equals(kmer_set_set.Get({})) = {}", i, i,
            kmer_sets[i].Equals(kmer_set_set.Get(i, n_workers), n_workers));
      }
    }
  }
}

int main(int argc, char** argv) {
  const std::vector<std::string> files = ParseFlags(argc, argv);

  const int k = absl::GetFlag(FLAGS_k);

  switch (k) {
    case 15:
      // 15 * 2 - 16 = 14 bits are used to select buckets.
      Main<15, uint16_t>(files);
      break;
    case 19:
      // 19 * 2 - 16 = 22 bits are used to select buckets.
      Main<19, uint16_t>(files);
      break;
    case 23:
      // 23 * 2 - 32 = 14 bits are used to select buckets.
      Main<23, uint32_t>(files);
      break;
    case 27:
      // 27 * 2 - 32 = 22 bits are used to select buckets.
      Main<27, uint32_t>(files);
      break;
    default:
      spdlog::error("unsupported k value");
      std::exit(1);
  }
}

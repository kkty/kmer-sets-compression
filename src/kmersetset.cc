#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/status/statusor.h"
#include "graph.h"
#include "kmer.h"
#include "kmer_counter.h"
#include "kmer_set.h"
#include "kmer_set_compact.h"
#include "kmer_set_set.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(bool, similarity, false, "calculate Jaccard similarity");
ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(std::string, decompressor, "",
          "specify decompressor for FASTQ files");
ABSL_FLAG(bool, canonical, false, "count canonical k-mers");
ABSL_FLAG(int, cutoff, 1, "cut off threshold");
ABSL_FLAG(int, workers, 1, "number of workers");
ABSL_FLAG(int, recursion, 1, "recursion limit for KmerSetSet");
ABSL_FLAG(bool, check, false, "check if k-mer sets can be reconstructed");
ABSL_FLAG(bool, nj, false, "use neighbor joining algorithm");
ABSL_FLAG(bool, mm, false, "use maximum weighted matching algorithm");
ABSL_FLAG(bool, amm, false,
          "use approximate maximum weighted matching algorithm");

int main(int argc, char** argv) {
  spdlog::set_default_logger(spdlog::stderr_color_mt("default"));

  // List of FASTQ file names.
  const std::vector<std::string> files = [&] {
    std::vector<std::string> v;

    for (const auto arg : absl::ParseCommandLine(argc, argv)) {
      v.push_back(arg);
    }

    // Program name.
    v.erase(v.begin());
    return v;
  }();

  if (absl::GetFlag(FLAGS_debug)) spdlog::set_level(spdlog::level::debug);

  const int n_workers = absl::GetFlag(FLAGS_workers);

  // 15 * 2 - 16 = 14 bits will be used to select a bucket.
  const int K = 15;
  using KeyType = uint16_t;

  const int n_datasets = files.size();

  std::vector<KmerSet<K, KeyType>> kmer_sets(n_datasets);

  for (int i = 0; i < n_datasets; i++) {
    const std::string& file = files[i];

    spdlog::info("constructing kmer_counter for {}", file);

    const KmerCounter<K, KeyType> kmer_counter = [&] {
      const std::string decompressor = absl::GetFlag(FLAGS_decompressor);
      const bool canonical = absl::GetFlag(FLAGS_canonical);

      absl::StatusOr<KmerCounter<K, KeyType>> status_or =
          decompressor != ""
              ? KmerCounter<K, KeyType>::FromFASTQ(file, decompressor,
                                                   canonical, n_workers)
              : KmerCounter<K, KeyType>::FromFASTQ(file, canonical, n_workers);

      if (!status_or.ok()) {
        spdlog::error("failed to parse FASTQ file: {}",
                      status_or.status().ToString());
      }

      return *status_or;
    }();

    spdlog::info("constructed kmer_counter for {}", file);
    spdlog::info("constructing kmer_set for {}", file);

    const auto [kmer_set, cutoff_count] =
        kmer_counter.ToSet(absl::GetFlag(FLAGS_cutoff), n_workers);

    spdlog::info("constructed kmer_set for {}", file);
    spdlog::info("cutoff_count = {}", cutoff_count);

    kmer_sets[i] = std::move(kmer_set);
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

  const auto cost_function =
      [canonical = absl::GetFlag(FLAGS_canonical)](
          const KmerSet<K, KeyType>& lhs, const KmerSet<K, KeyType>& rhs,
          int n_workers) { return lhs.Diff(rhs, n_workers); };

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
        kmer_sets, absl::GetFlag(FLAGS_recursion), cost_function, n_workers);

    spdlog::info("constructed kmer_set_set_mm");

    spdlog::info("kmer_set_set_mm.Size() = {}", kmer_set_set_mm.Size());
    spdlog::info("kmer_set_set_mm.Cost() = {}", kmer_set_set_mm.Cost());

    if (absl::GetFlag(FLAGS_check)) {
      for (int i = 0; i < n_datasets; i++) {
        spdlog::info(
            "kmer_sets[{}].Equals(kmer_set_set_mm.Get({})) = {}", i, i,
            kmer_sets[i].Equals(kmer_set_set_mm.Get(i, n_workers), n_workers));
      }
    }
  } else if (absl::GetFlag(FLAGS_amm)) {
    spdlog::info("constructing kmer_set_set_amm");

    KmerSetSetAMM<K, KeyType, decltype(cost_function)> kmer_set_set_amm(
        kmer_sets, absl::GetFlag(FLAGS_recursion), cost_function, n_workers);

    spdlog::info("constructed kmer_set_set_amm");

    spdlog::info("kmer_set_set_amm.Size() = {}", kmer_set_set_amm.Size());
    spdlog::info("kmer_set_set_amm.Cost() = {}", kmer_set_set_amm.Cost());

    if (absl::GetFlag(FLAGS_check)) {
      for (int i = 0; i < n_datasets; i++) {
        spdlog::info(
            "kmer_sets[{}].Equals(kmer_set_set_amm.Get({})) = {}", i, i,
            kmer_sets[i].Equals(kmer_set_set_amm.Get(i, n_workers), n_workers));
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

  return 0;
}

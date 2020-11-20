#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/status/statusor.h"
#include "core/graph.h"
#include "core/kmer.h"
#include "core/kmer_counter.h"
#include "core/kmer_set.h"
#include "core/kmer_set_compressed.h"
#include "core/kmer_set_set.h"
#include "flags.h"
#include "io.h"
#include "log.h"
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
ABSL_FLAG(bool, mm, false, "use maximum weighted matching algorithm");
ABSL_FLAG(bool, approximate_matching, false,
          "use approximate algorithm for matching");
ABSL_FLAG(bool, approximate_weights, false,
          "use approximate weights for the matching algorithm");
ABSL_FLAG(bool, approximate_graph, false,
          "use an approximate (sparse) graph for the matching algorithm");
ABSL_FLAG(std::string, out, "", "path to save dumped file");
ABSL_FLAG(std::string, compressor, "", "program to compress dumped file");

template <int K, typename KeyType>
void Main(const std::vector<std::string>& files) {
  InitDefaultLogger();

  if (absl::GetFlag(FLAGS_debug)) EnableDebugLogs();

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

  int64_t total_size = 0;
  for (int i = 0; i < n_datasets; i++) {
    int64_t kmer_set_size = kmer_sets[i].Size();
    spdlog::info("i = {}, kmer_set_size = {}", i, kmer_set_size);
    total_size += kmer_set_size;
  }

  spdlog::info("total_size = {}", total_size);

  if (absl::GetFlag(FLAGS_similarity)) {
    for (int i = 0; i < n_datasets; i++) {
      for (int j = i + 1; j < n_datasets; j++) {
        double similarity = kmer_sets[i].Similarity(kmer_sets[j], n_workers);
        spdlog::info("i = {}, j = {}, similarity = {}", i, j, similarity);
      }
    }
  }

  if (absl::GetFlag(FLAGS_mm)) {
    spdlog::info("constructing kmer_set_set_mm");

    KmerSetSetMM<K, KeyType> kmer_set_set_mm = [&] {
      if (absl::GetFlag(FLAGS_check)) {
        return KmerSetSetMM<K, KeyType>(
            kmer_sets, absl::GetFlag(FLAGS_recursion),
            absl::GetFlag(FLAGS_approximate_matching),
            absl::GetFlag(FLAGS_approximate_weights),
            absl::GetFlag(FLAGS_approximate_graph), n_workers);
      } else {
        // We can move kmer_sets.
        return KmerSetSetMM<K, KeyType>(
            std::move(kmer_sets), absl::GetFlag(FLAGS_recursion),
            absl::GetFlag(FLAGS_approximate_matching),
            absl::GetFlag(FLAGS_approximate_weights),
            absl::GetFlag(FLAGS_approximate_graph), n_workers);
      }
    }();

    spdlog::info("constructed kmer_set_set_mm");

    spdlog::info("kmer_set_set_mm.Size() = {}", kmer_set_set_mm.Size());
    spdlog::info("kmer_set_set_mm.Cost() = {}", kmer_set_set_mm.Cost());

    const std::string out_file = absl::GetFlag(FLAGS_out);

    if (out_file != "") {
      kmer_set_set_mm.Dump(out_file, absl::GetFlag(FLAGS_compressor), canonical,
                           n_workers);
    }

    if (absl::GetFlag(FLAGS_check)) {
      spdlog::info("dumping");
      std::vector<std::string> dumped =
          kmer_set_set_mm.Dump(absl::GetFlag(FLAGS_canonical), n_workers);
      spdlog::info("dumped");

      spdlog::info("loading");
      std::unique_ptr<KmerSetSetMM<K, KeyType>> loaded =
          std::unique_ptr<KmerSetSetMM<K, KeyType>>(
              KmerSetSetMM<K, KeyType>::Load(std::move(dumped),
                                             absl::GetFlag(FLAGS_canonical),
                                             n_workers));
      spdlog::info("loaded");

      for (int i = 0; i < n_datasets; i++) {
        spdlog::info("kmer_sets[{}].Equals(loaded->Get({})) = {}", i, i,
                     kmer_sets[i].Equals(loaded->Get(i, n_workers), n_workers));
      }
    }
  } else {
    spdlog::info("constructing kmer_set_set");

    KmerSetSet<K, KeyType> kmer_set_set(
        kmer_sets, absl::GetFlag(FLAGS_recursion), n_workers);
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

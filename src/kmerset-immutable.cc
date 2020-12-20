#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "absl/flags/flag.h"
#include "core/kmer_set.h"
#include "core/kmer_set_immutable.h"
#include "flags.h"
#include "log.h"
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

ABSL_FLAG(int, k, 15, "the length of kmers");
ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(int, workers, 1, "number of workers");
ABSL_FLAG(int, kmers, 10000000, "number of kmers");
ABSL_FLAG(int, repeats, 100, "number of repeats");

template <int K, int N, typename KeyType>
void Main() {
  InitDefaultLogger();

  const bool debug = absl::GetFlag(FLAGS_debug);
  const int n_workers = absl::GetFlag(FLAGS_workers);
  const int n_kmers = absl::GetFlag(FLAGS_kmers);
  const int n_repeats = absl::GetFlag(FLAGS_repeats);

  if (debug) EnableDebugLogs();

  std::vector<KmerSetImmutable<K, N, KeyType>> kmer_sets_immutable;
  std::vector<KmerSetImmutableHashSet<K>> kmer_sets_immutable_hash_set;
  std::vector<KmerSetImmutableVector<K>> kmer_sets_immutable_vector;

  for (int i = 0; i < 2; i++) {
    std::vector<Kmer<K>> v;
    for (int j = 0; j < n_kmers; j++) {
      v.push_back(GetRandomKmer<K>());
    }

    KmerSet<K, N, KeyType> kmer_set(v, n_workers);

    kmer_sets_immutable.push_back(
        KmerSetImmutable<K, N, KeyType>(kmer_set, n_workers));

    kmer_sets_immutable_hash_set.push_back(
        KmerSetImmutableHashSet<K>::FromKmerSet(kmer_set, n_workers));

    kmer_sets_immutable_vector.push_back(
        KmerSetImmutableVector<K>::FromKmerSet(kmer_set, n_workers));

    spdlog::info("kmer_sets_immutable[{}].Size() = {}", i,
                 kmer_sets_immutable[i].Size());
    spdlog::info("kmer_sets_immutable[{}].Bytes() = {}", i,
                 kmer_sets_immutable[i].Bytes());
  }

  {
    spdlog::info("intersecting kmer_sets_immutable");
    spdlog::stopwatch sw;
    std::int64_t total_size = 0;
    for (int i = 0; i < n_repeats; i++) {
      total_size += kmer_sets_immutable[0]
                        .Intersection(kmer_sets_immutable[1], n_workers)
                        .Size();
    }
    spdlog::info("total_size = {}, elapsed = {}", total_size, sw);
    std::cout << sw.elapsed().count() << ' ';
  }

  {
    spdlog::info("intersecting kmer_sets_immutable_hash_set");
    spdlog::stopwatch sw;
    std::int64_t total_size = 0;
    for (int i = 0; i < n_repeats; i++) {
      total_size += kmer_sets_immutable_hash_set[0]
                        .Intersection(kmer_sets_immutable_hash_set[1])
                        .Size();
    }
    spdlog::info("total_size = {}, elapsed = {}", total_size, sw);
    std::cout << sw.elapsed().count() << ' ';
  }

  {
    spdlog::info("intersecting kmer_sets_immutable_vector");
    spdlog::stopwatch sw;
    std::int64_t total_size = 0;
    for (int i = 0; i < n_repeats; i++) {
      total_size += kmer_sets_immutable_vector[0]
                        .Intersection(kmer_sets_immutable_vector[1])
                        .Size();
    }
    spdlog::info("total_size = {}, elapsed = {}", total_size, sw);
    std::cout << sw.elapsed().count() << '\n';
  }

  {
    spdlog::info("adding kmer_sets_immutable");
    spdlog::stopwatch sw;
    std::int64_t total_size = 0;
    for (int i = 0; i < n_repeats; i++) {
      total_size +=
          kmer_sets_immutable[0].Add(kmer_sets_immutable[1], n_workers).Size();
    }
    spdlog::info("total_size = {}, elapsed = {}", total_size, sw);
    std::cout << sw.elapsed().count() << ' ';
  }

  {
    spdlog::info("adding kmer_sets_immutable_hash_set");
    spdlog::stopwatch sw;
    std::int64_t total_size = 0;
    for (int i = 0; i < n_repeats; i++) {
      total_size += kmer_sets_immutable_hash_set[0]
                        .Add(kmer_sets_immutable_hash_set[1])
                        .Size();
    }
    spdlog::info("total_size = {}, elapsed = {}", total_size, sw);
    std::cout << sw.elapsed().count() << ' ';
  }

  {
    spdlog::info("adding kmer_sets_immutable_vector");
    spdlog::stopwatch sw;
    std::int64_t total_size = 0;
    for (int i = 0; i < n_repeats; i++) {
      total_size += kmer_sets_immutable_vector[0]
                        .Add(kmer_sets_immutable_vector[1])
                        .Size();
    }
    spdlog::info("total_size = {}, elapsed = {}", total_size, sw);
    std::cout << sw.elapsed().count() << '\n';
  }

  {
    spdlog::info("subtracting kmer_sets_immutable");
    spdlog::stopwatch sw;
    std::int64_t total_size = 0;
    for (int i = 0; i < n_repeats; i++) {
      total_size +=
          kmer_sets_immutable[0].Sub(kmer_sets_immutable[1], n_workers).Size();
    }
    spdlog::info("total_size = {}, elapsed = {}", total_size, sw);
    std::cout << sw.elapsed().count() << ' ';
  }

  {
    spdlog::info("subtracting kmer_sets_immutable_hash_set");
    spdlog::stopwatch sw;
    std::int64_t total_size = 0;
    for (int i = 0; i < n_repeats; i++) {
      total_size += kmer_sets_immutable_hash_set[0]
                        .Sub(kmer_sets_immutable_hash_set[1])
                        .Size();
    }
    spdlog::info("total_size = {}, elapsed = {}", total_size, sw);
    std::cout << sw.elapsed().count() << ' ';
  }

  {
    spdlog::info("subtracting kmer_sets_immutable_vector");
    spdlog::stopwatch sw;
    std::int64_t total_size = 0;
    for (int i = 0; i < n_repeats; i++) {
      total_size += kmer_sets_immutable_vector[0]
                        .Sub(kmer_sets_immutable_vector[1])
                        .Size();
    }
    spdlog::info("total_size = {}, elapsed = {}", total_size, sw);
    std::cout << sw.elapsed().count() << '\n';
  }
}

int main(int argc, char** argv) {
  ParseFlags(argc, argv);

  const int k = absl::GetFlag(FLAGS_k);

  switch (k) {
    case 15:
      Main<15, 14, uint16_t>();
      break;
    case 19:
      Main<19, 10, uint32_t>();
      break;
    case 23:
      Main<23, 14, uint32_t>();
      break;
    default:
      spdlog::error("unsupported k value");
      std::exit(1);
  }
}

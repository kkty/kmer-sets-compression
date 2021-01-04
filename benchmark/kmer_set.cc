#include <cstdint>

#include "benchmark/benchmark.h"
#include "core/kmer_set_immutable.h"
#include "random.h"

void Benchmark_KmerSet_Find(benchmark::State& state) {
  const int K = 11;
  const int N = 14;
  using KeyType = std::uint8_t;

  KmerSet<K, N, KeyType> kmer_set =
      GetRandomKmerSet<K, N, KeyType>(1'000'000, true);

  for (auto _ : state) {
    benchmark::DoNotOptimize(kmer_set.Find(state.range(0)));
  }
}

BENCHMARK(Benchmark_KmerSet_Find)->RangeMultiplier(2)->Range(1, 8);

BENCHMARK_MAIN();

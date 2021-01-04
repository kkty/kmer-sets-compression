#include "core/kmer_counter.h"

#include <cstdint>

#include "benchmark/benchmark.h"
#include "random.h"

void Benchmark_KmerCounter_ToKmerSet(benchmark::State& state) {
  const int K = 11;
  const int N = 14;
  using KeyType = std::uint8_t;

  KmerCounter<K, N, KeyType> kmer_counter =
      GetRandomKmerCounter<K, N, KeyType>(1'000'000, true);

  for (auto _ : state) {
    benchmark::DoNotOptimize(kmer_counter.ToKmerSet(2, state.range(0)));
  }
}

BENCHMARK(Benchmark_KmerCounter_ToKmerSet)->RangeMultiplier(2)->Range(1, 8);

BENCHMARK_MAIN();

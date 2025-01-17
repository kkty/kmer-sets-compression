#include "core/spss.h"

#include <cstdint>
#include <string>
#include <thread>
#include <vector>

#include "benchmark/benchmark.h"
#include "core/kmer_set.h"
#include "random.h"

void Benchmark_GetUnitigsCanonical(benchmark::State& state) {
  const int K = 11;
  const int N = 14;
  using KeyType = std::uint8_t;

  KmerSet<K, N, KeyType> kmer_set =
      GetRandomKmerSet<K, N, KeyType>(1'000'000, true);

  for (auto _ : state) {
    benchmark::DoNotOptimize(GetUnitigsCanonical(kmer_set, state.range(0)));
  }
}

BENCHMARK(Benchmark_GetUnitigsCanonical)->RangeMultiplier(2)->Range(1, 8);

void Benchmark_GetSPSSCanonical(benchmark::State& state) {
  const int K = 11;
  const int N = 14;
  using KeyType = std::uint8_t;

  KmerSet<K, N, KeyType> kmer_set =
      GetRandomKmerSet<K, N, KeyType>(1'000'000, true);

  for (auto _ : state) {
    benchmark::DoNotOptimize(GetSPSSCanonical(kmer_set, true, state.range(0)));
  }
}

BENCHMARK(Benchmark_GetSPSSCanonical)->RangeMultiplier(2)->Range(1, 8);

void Benchmark_GetKmerSetFromSPSS(benchmark::State& state) {
  const int K = 11;
  const int N = 14;
  using KeyType = std::uint8_t;

  KmerSet<K, N, KeyType> kmer_set =
      GetRandomKmerSet<K, N, KeyType>(1'000'000, true);

  std::vector<std::string> spss =
      GetSPSSCanonical(kmer_set, true, std::thread::hardware_concurrency());

  for (auto _ : state) {
    benchmark::DoNotOptimize(
        GetKmerSetFromSPSS<K, N, KeyType>(spss, true, state.range(0)));
  }
}

BENCHMARK(Benchmark_GetKmerSetFromSPSS)->RangeMultiplier(2)->Range(1, 8);

BENCHMARK_MAIN();

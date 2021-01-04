#include "core/kmer.h"

#include <string>

#include "benchmark/benchmark.h"
#include "random.h"

void Benchmark_Kmer_FromString(benchmark::State& state) {
  const int K = 15;

  std::string s = GetRandomKmer<K>().String();

  for (auto _ : state) {
    benchmark::DoNotOptimize(Kmer<K>(s));
  }
}

BENCHMARK(Benchmark_Kmer_FromString);

void Benchmark_Kmer_Complement(benchmark::State& state) {
  const int K = 15;

  Kmer<K> kmer = GetRandomKmer<K>();

  for (auto _ : state) {
    benchmark::DoNotOptimize(kmer.Complement());
  }
}

BENCHMARK(Benchmark_Kmer_Complement);

BENCHMARK_MAIN();

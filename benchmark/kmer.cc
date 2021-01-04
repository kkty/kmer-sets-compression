#include "core/kmer.h"

#include <string>

#include "benchmark/benchmark.h"

void Benchmark_Kmer_FromString(benchmark::State& state) {
  const int K = 15;

  std::string s("ACGTACGTACGTACG");

  for (auto _ : state) {
    benchmark::DoNotOptimize(Kmer<K>(s));
  }
}

BENCHMARK(Benchmark_Kmer_FromString);

void Benchmark_Kmer_Complement(benchmark::State& state) {
  const int K = 15;

  Kmer<K> kmer("ACGTACGTACGTACG");

  for (auto _ : state) {
    benchmark::DoNotOptimize(kmer.Complement());
  }
}

BENCHMARK(Benchmark_Kmer_Complement);

void Benchmark_Kmer_Next(benchmark::State& state) {
  const int K = 15;

  Kmer<K> kmer("ACGTACGTACGTACG");

  for (auto _ : state) {
    benchmark::DoNotOptimize(kmer.Next('T'));
  }
}

BENCHMARK(Benchmark_Kmer_Next);

BENCHMARK_MAIN();

#include "core/int_set.h"

#include <vector>

#include "benchmark/benchmark.h"
#include "random.h"

void Benchmark_IntSet(benchmark::State& state) {
  std::vector<int> v = GetRandomInts(1'000'000, true, true);

  for (auto _ : state) {
    IntSet<int> int_set(v);
    benchmark::DoNotOptimize(int_set.Decode());
  }
}

BENCHMARK(Benchmark_IntSet);

void Benchmark_IntSet_Intersection(benchmark::State& state) {
  std::vector<int> v1 = GetRandomInts(1'000'000, true, true);
  std::vector<int> v2 = GetRandomInts(1'000'000, true, true);

  IntSet<int> int_set1(v1);
  IntSet<int> int_set2(v2);

  for (auto _ : state) {
    benchmark::DoNotOptimize(int_set1.Intersection(int_set2));
  }
}

BENCHMARK(Benchmark_IntSet_Intersection);

BENCHMARK_MAIN();

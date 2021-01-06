#include "core/int_set.h"

#include <cstdint>
#include <vector>

#include "benchmark/benchmark.h"
#include "random.h"

void Benchmark_IntSet_Int(benchmark::State& state) {
  std::vector<int> v = GetRandomInts(1'000'000, true, true);

  for (auto _ : state) {
    IntSet<int> int_set(v);
    benchmark::DoNotOptimize(int_set.Decode());
  }
}

BENCHMARK(Benchmark_IntSet_Int);

void Benchmark_IntSet_Uint32(benchmark::State& state) {
  std::vector<std::uint32_t> v =
      GetRandomInts<std::uint32_t>(1'000'000, true, true);

  for (auto _ : state) {
    IntSet<std::uint32_t> int_set(v);
    benchmark::DoNotOptimize(int_set.Decode());
  }
}

BENCHMARK(Benchmark_IntSet_Uint32);

void Benchmark_IntSet_Intersection_Int(benchmark::State& state) {
  std::vector<int> v1 = GetRandomInts(1'000'000, true, true);
  std::vector<int> v2 = GetRandomInts(1'000'000, true, true);

  IntSet<int> int_set1(v1);
  IntSet<int> int_set2(v2);

  for (auto _ : state) {
    benchmark::DoNotOptimize(int_set1.Intersection(int_set2));
  }
}

BENCHMARK(Benchmark_IntSet_Intersection_Int);

void Benchmark_IntSet_Intersection_Uint32(benchmark::State& state) {
  std::vector<std::uint32_t> v1 =
      GetRandomInts<std::uint32_t>(1'000'000, true, true);
  std::vector<std::uint32_t> v2 =
      GetRandomInts<std::uint32_t>(1'000'000, true, true);

  IntSet<std::uint32_t> int_set1(v1);
  IntSet<std::uint32_t> int_set2(v2);

  for (auto _ : state) {
    benchmark::DoNotOptimize(int_set1.Intersection(int_set2));
  }
}

BENCHMARK(Benchmark_IntSet_Intersection_Uint32);

BENCHMARK_MAIN();

#include "core/parallel_disjoint_set.h"

#include <atomic>
#include <thread>
#include <tuple>
#include <vector>

#include "absl/random/random.h"
#include "core/range.h"
#include "gtest/gtest.h"

class DisjointSet {
 public:
  DisjointSet(int size) : parents_(size) {}

  int Find(int i) {
    if (parents_[i] == i) return i;
    return Find(parents_[i]);
  }

  void Unite(int x, int y) {
    x = Find(x);
    y = Find(y);
    parents_[x] = y;
  }

  bool IsSame(int x, int y) { return Find(x) == Find(y); }

 private:
  std::vector<int> parents_;
};

TEST(ParallelDisjointSet, Random) {
  const int n_workers = 8;
  const int n = 10000;
  const int m = 2500;

  std::vector<std::pair<int, int>> pairs;

  {
    absl::InsecureBitGen bitgen;
    for (int i = 0; i < m; i++) {
      int x = absl::Uniform(bitgen, 0, n);
      int y = absl::Uniform(bitgen, 0, n);
      pairs.emplace_back(x, y);
    }
  }

  DisjointSet disjoint_set(n);

  for (const std::pair<int, int>& p : pairs) {
    disjoint_set.Unite(p.first, p.second);
  }

  ParallelDisjointSet parallel_disjoint_set(n);

  {
    std::vector<std::thread> threads;

    for (const Range& range : Range(0, pairs.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        for (int i : range) {
          parallel_disjoint_set.Unite(pairs[i].first, pairs[i].second);
        }
      });
    }

    for (std::thread& t : threads) t.join();
  }

  {
    std::vector<std::thread> threads;
    std::atomic_int count = 0;

    for (const Range& range : Range(0, n).Split(n_workers)) {
      threads.emplace_back([&, range] {
        for (int i : range) {
          for (int j = 0; j < n; j++) {
            if (disjoint_set.IsSame(i, j) ==
                parallel_disjoint_set.IsSame(i, j)) {
              count += 1;
            }
          }
        }
      });
    }

    for (std::thread& t : threads) t.join();

    ASSERT_EQ(count, n * n);
  }
}

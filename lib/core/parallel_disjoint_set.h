#ifndef PARALLEL_DISJOINT_SET_H_
#define PARALLEL_DISJOINT_SET_H_

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <vector>

// ParallelDisjointSet is an implementation of the disjoint set data structure
// that can be used by multiple threads.
// It is based on "Anderson, R. J., & Woll, H. (1991, January). Wait-free
// parallel algorithms for the union-find problem. In Proceedings of the
// twenty-third annual ACM symposium on Theory of computing (pp. 370-380)".
class ParallelDisjointSet {
 public:
  ParallelDisjointSet(int size) : a_(size) {
    for (int i = 0; i < size; i++) {
      a_[i] = i;
    }
  }

  int Find(int x) {
    int y = x;

    while (x != GetNext(x)) {
      x = GetNext(x);
    }

    while (LessThan(y, x)) {
      std::uint64_t expected =
          (static_cast<std::uint64_t>(GetRank(y)) << 32) + GetNext(y);
      std::uint64_t desired = ((expected >> 32) << 32) + x;
      a_[y].compare_exchange_weak(expected, desired);
      y = GetNext(y);
    }

    return x;
  }

  bool IsSame(int x, int y) {
    while (true) {
      x = Find(x);
      y = Find(y);
      if (x == y) return true;
      if (GetNext(x) == x) return false;
    }
  }

  void Unite(int x, int y) {
    while (true) {
      x = Find(x);
      y = Find(y);

      if (x == y) return;

      int rank_x = GetRank(x);
      int rank_y = GetRank(y);

      if (rank_x > rank_y || (rank_x == rank_y && x > y)) {
        std::swap(x, y);
        std::swap(rank_x, rank_y);
      }

      if (!UpdateRoot(x, rank_x, y, rank_x)) {
        continue;
      }

      if (rank_x == rank_y) {
        UpdateRoot(y, rank_y, y, rank_y + 1);
      }

      break;
    }
  }

 private:
  int GetRank(int i) const { return a_[i] >> 32; }

  int GetNext(int i) const { return (a_[i] << 32) >> 32; }

  bool UpdateRoot(int x, int old_rank, int y, int new_rank) {
    std::uint64_t old = a_[x];
    if ((old << 32) >> 32 != x || old >> 32 != old_rank) return false;
    std::uint64_t updated = (static_cast<std::uint64_t>(new_rank) << 32) + y;
    return a_[x].compare_exchange_weak(old, updated);
  }

  bool LessThan(int x, int y) const {
    int rank_x = GetRank(x);
    int rank_y = GetRank(y);

    if (rank_x < rank_y) return true;
    if (rank_x > rank_y) return false;

    return x < y;
  }

  std::vector<std::atomic_uint64_t> a_;
};

#endif

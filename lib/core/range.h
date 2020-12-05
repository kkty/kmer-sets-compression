#ifndef CORE_RANGE_H_
#define CORE_RANGE_H_

#include <cstdlib>
#include <vector>

struct Range {
  Range(int64_t begin, int64_t end) : begin(begin), end(end) {}

  int64_t Size() const { return end - begin; }

  // Splits the range to n distinct ranges. If n is smaller than Size(), some
  // resulting ranges have a size of zero.
  std::vector<Range> Split(int64_t n) const {
    const int64_t small_chunk_size = Size() / n;
    const int64_t large_chunk_size = small_chunk_size + 1;

    const int64_t large_chunk_n = Size() - small_chunk_size * n;
    const int64_t small_chunk_n = n - large_chunk_n;

    std::vector<Range> ranges;
    ranges.reserve(n);

    for (int64_t i = 0; i < small_chunk_n; i++) {
      if (i == 0) {
        ranges.emplace_back(begin, begin + small_chunk_size);
      } else {
        ranges.emplace_back(ranges.back().end,
                            ranges.back().end + small_chunk_size);
      }
    }

    for (int64_t i = 0; i < large_chunk_n; i++) {
      ranges.emplace_back(ranges.back().end,
                          ranges.back().end + large_chunk_size);
    }

    return ranges;
  }

  // Executes a function for each integer in the range.
  // Range(...).ForEach([&](int64_t i) { ... });
  template <typename F>
  void ForEach(F f) const {
    for (int64_t i = begin; i < end; i++) {
      f(i);
    }
  }

  const int64_t begin;
  const int64_t end;
};

#endif
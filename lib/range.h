#ifndef RANGE_H_
#define RANGE_H_

#include <cstdlib>
#include <vector>

struct Range {
  Range(int64_t begin, int64_t end) : begin(begin), end(end) {}

  int64_t Size() const { return end - begin; }

  std::vector<Range> Split(int64_t n) const {
    const int small_chunk_size = Size() / n;
    const int large_chunk_size = small_chunk_size + 1;

    const int large_chunk_n = Size() - small_chunk_size * n;
    const int small_chunk_n = n - large_chunk_n;

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

  const int64_t begin;
  const int64_t end;
};

#endif
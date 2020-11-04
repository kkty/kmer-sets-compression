#ifndef RANGE_H_
#define RANGE_H_

#include <cstdlib>
#include <vector>

std::vector<std::pair<int64_t, int64_t>> SplitRange(int64_t begin, int64_t end,
                                                    int64_t n) {
  const int small_chunk_size = (end - begin) / n;
  const int large_chunk_size = small_chunk_size + 1;

  const int large_chunk_n = (end - begin) - small_chunk_size * n;
  const int small_chunk_n = n - large_chunk_n;

  std::vector<std::pair<int64_t, int64_t>> ranges;
  ranges.reserve(n);

  for (int64_t i = 0; i < small_chunk_n; i++) {
    if (i == 0) {
      ranges.emplace_back(begin, begin + small_chunk_size);
    } else {
      ranges.emplace_back(ranges.back().second,
                          ranges.back().second + small_chunk_size);
    }
  }

  for (int64_t i = 0; i < large_chunk_n; i++) {
    ranges.emplace_back(ranges.back().second,
                        ranges.back().second + large_chunk_size);
  }

  return ranges;
}

#endif
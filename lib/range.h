#ifndef RANGE_H_
#define RANGE_H_

#include <cstdlib>
#include <vector>

std::vector<std::pair<int64_t, int64_t>> SplitRange(int64_t begin, int64_t end,
                                                    int64_t n) {
  int chunk_size = (end - begin + n - 1) / n;

  std::vector<std::pair<int64_t, int64_t>> ranges;
  for (int64_t i = 0; i < n; i++) {
    ranges.emplace_back(i * chunk_size, std::min((i + 1) * chunk_size, end));
  }

  return ranges;
}

#endif
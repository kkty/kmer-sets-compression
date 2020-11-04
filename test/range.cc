#include "range.h"

#include "gtest/gtest.h"

TEST(SplitRange, Matrix) {
  for (int begin = 0; begin < 100; begin++) {
    for (int end = begin; end < 100; end++) {
      for (int n = 1; n < 100; n++) {
        const auto ranges = SplitRange(begin, end, n);

        ASSERT_EQ(ranges.size(), n);
        ASSERT_EQ(ranges.front().first, begin);
        for (int i = 0; i < n - 1; i++)
          ASSERT_EQ(ranges[i].second, ranges[i + 1].first);
        ASSERT_EQ(ranges.back().second, end);
      }
    }
  }
}

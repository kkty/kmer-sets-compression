#include "range.h"

#include "gtest/gtest.h"

TEST(Range, Split) {
  for (int begin = 0; begin < 100; begin++) {
    for (int end = begin; end < 100; end++) {
      for (int n = 1; n < 100; n++) {
        const auto ranges = Range(begin, end).Split(n);

        ASSERT_EQ(ranges.size(), n);
        ASSERT_EQ(ranges.front().begin, begin);
        for (int i = 0; i < n - 1; i++)
          ASSERT_EQ(ranges[i].begin, ranges[i + 1].end);
        ASSERT_EQ(ranges.back().end, end);
      }
    }
  }
}

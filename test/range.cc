#include "core/range.h"

#include "gtest/gtest.h"

TEST(Range, Split) {
  for (int begin = 0; begin < 100; begin++) {
    for (int end = begin; end < 100; end++) {
      for (int n = 1; n < 100; n++) {
        const auto ranges = Range(begin, end).Split(n);

        ASSERT_EQ(ranges.size(), n);
        ASSERT_EQ(ranges.front().Begin(), begin);
        for (int i = 0; i < n - 1; i++)
          ASSERT_EQ(ranges[i].End(), ranges[i + 1].Begin());
        ASSERT_EQ(ranges.back().End(), end);
      }
    }
  }
}

TEST(Range, ForLoop) {
  int sum = 0;
  for (int i : Range(0, 10)) {
    sum += i;
  }
  ASSERT_EQ(sum, 45);
}

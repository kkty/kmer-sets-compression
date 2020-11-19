#include "core/suffix_array.h"

#include <vector>

#include "gtest/gtest.h"

TEST(SuffixArray, Find) {
  SuffixArray sa("mississippi", 1);

  {
    std::vector<int64_t> v = sa.Find("m");
    ASSERT_EQ(v.size(), 1);
    ASSERT_EQ(v[0], 0);
  }

  {
    std::vector<int64_t> v = sa.Find("ss");
    ASSERT_EQ(v.size(), 2);
    ASSERT_EQ(v[0], 2);
    ASSERT_EQ(v[1], 5);
  }

  {
    std::vector<int64_t> v = sa.Find("sss");
    ASSERT_EQ(v.size(), 0);
  }
}

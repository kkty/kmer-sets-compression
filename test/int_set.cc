#include "core/int_set.h"

#include <algorithm>
#include <limits>
#include <vector>

#include "absl/random/random.h"
#include "gtest/gtest.h"

TEST(IntSet, MoveConstructor) {
  std::vector<int> v{1, 3, 5, 7, 9};
  IntSet<int> int_set(v);
  IntSet<int> moved = std::move(int_set);
  ASSERT_EQ(moved.Decode(), v);
  ASSERT_EQ(int_set.Size(), 0);  // NOLINT
}

TEST(IntSet, MoveAssignment) {
  std::vector<int> v{1, 3, 5, 7, 9};
  IntSet<int> int_set(v);
  IntSet<int> moved;
  moved = std::move(int_set);
  ASSERT_EQ(moved.Decode(), v);
  ASSERT_EQ(int_set.Size(), 0);  // NOLINT
}

TEST(IntSet, CopyConstructor) {
  std::vector<int> v{1, 3, 5, 7, 9};
  IntSet<int> int_set(v);
  IntSet<int> copied(int_set);  // NOLINT
  ASSERT_EQ(copied.Decode(), v);
}

TEST(IntSet, CopyAssignment) {
  std::vector<int> v{1, 3, 5, 7, 9};
  IntSet<int> int_set(v);
  IntSet<int> copied;
  copied = int_set;
  ASSERT_EQ(copied.Decode(), v);
}

TEST(IntSet, EncodeAndDecodeRandom) {
  absl::InsecureBitGen bitgen;

  const int n = 1000000;
  std::vector<int> v;
  v.reserve(n);

  for (int i = 0; i < n; i++) {
    v.push_back(absl::Uniform<uint32_t>(absl::IntervalClosedClosed, bitgen,
                                        std::numeric_limits<int>::min(),
                                        std::numeric_limits<int>::max()));
  }

  std::sort(v.begin(), v.end());

  {
    auto it = std::unique(v.begin(), v.end());
    v.resize(std::distance(v.begin(), it));
  }

  IntSet int_set(v);

  ASSERT_EQ(int_set.Size(), v.size());

  std::vector<int> decoded = int_set.Decode();

  ASSERT_EQ(v, decoded);
}

TEST(IntSet, IntersectionAddSub) {
  std::vector<int> v1{1, 3, 5, 7, 9};
  std::vector<int> v2{0, 1, 5, 6, 9, 10};

  {
    std::vector<int> expected{1, 5, 9};
    ASSERT_EQ(IntSet(v1).Intersection(IntSet(v2)).Decode(), expected);
    ASSERT_EQ(IntSet(v1).Intersection(IntSet(v2)).Decode(), expected);
  }

  {
    std::vector<int> expected{0, 1, 3, 5, 6, 7, 9, 10};
    ASSERT_EQ(IntSet(v1).Add(IntSet(v2)).Decode(), expected);
  }

  {
    std::vector<int> expected{3, 7};
    ASSERT_EQ(IntSet(v1).Sub(IntSet(v2)).Decode(), expected);
  }
}

TEST(IntSet, IntersectionAddSubRandom) {
  absl::InsecureBitGen bitgen;

  const auto GetData = [&] {
    std::vector<int> v;

    const int n = absl::Uniform(bitgen, 0, 10000);
    const int m = absl::Uniform(bitgen, 0, 10000);

    v.reserve(n);

    for (int i = 0; i < n; i++) {
      v.push_back(absl::Uniform(absl::IntervalClosedClosed, bitgen, 0, m));
    }

    std::sort(v.begin(), v.end());

    // Removes duplicates.
    auto it = std::unique(v.begin(), v.end());
    v.resize(std::distance(v.begin(), it));

    return v;
  };

  std::vector<int> v1 = GetData();
  std::vector<int> v2 = GetData();

  {
    std::vector<int> expected;
    expected.insert(expected.end(), v1.begin(), v1.end());
    expected.insert(expected.end(), v2.begin(), v2.end());
    std::sort(expected.begin(), expected.end());
    auto it = std::unique(expected.begin(), expected.end());
    expected.resize(std::distance(expected.begin(), it));

    ASSERT_EQ(IntSet(v1).Add(IntSet(v2)).Decode(), expected);
  }

  {
    std::vector<int> expected;

    for (int i : v1) {
      if (std::find(v2.begin(), v2.end(), i) != v2.end()) {
        expected.push_back(i);
      }
    }

    ASSERT_EQ(IntSet(v1).Intersection(IntSet(v2)).Decode(), expected);
  }

  {
    std::vector<int> expected;

    for (int i : v1) {
      if (std::find(v2.begin(), v2.end(), i) == v2.end()) {
        expected.push_back(i);
      }
    }

    ASSERT_EQ(IntSet(v1).Sub(IntSet(v2)).Decode(), expected);
  }
}

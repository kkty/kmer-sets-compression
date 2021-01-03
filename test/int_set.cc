#include "core/int_set.h"

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/random/random.h"
#include "gtest/gtest.h"

std::vector<int> GetRandomInts(int n, bool is_unique = false,
                               bool is_sorted = false,
                               int min = std::numeric_limits<int>::min(),
                               int max = std::numeric_limits<int>::max()) {
  absl::InsecureBitGen bitgen;

  std::vector<int> v;

  if (is_unique) {
    absl::flat_hash_set<int> s;

    while (static_cast<int>(s.size()) < n) {
      s.insert(absl::Uniform(absl::IntervalClosed, bitgen, min, max));
    }

    v.insert(v.end(), s.begin(), s.end());
  } else {
    for (int i = 0; i < n; i++) {
      v.push_back(absl::Uniform(absl::IntervalClosed, bitgen, min, max));
    }
  }

  if (is_sorted) {
    std::sort(v.begin(), v.end());
  }

  return v;
}

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

  const int n = 100000;
  std::vector<int> v = GetRandomInts(n, true, true);

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
    const int n = absl::Uniform(absl::IntervalClosed, bitgen, 0, 10000);
    const int m = absl::Uniform(absl::IntervalClosed, bitgen, 2 * n, 8 * n);

    return GetRandomInts(n, true, true, 0, m);
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

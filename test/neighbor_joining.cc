#include "core/neighbor_joining.h"

#include <string>
#include <unordered_map>

#include "gtest/gtest.h"

using namespace neighbor_joining;

TEST(NeighborJoining, Execute) {
  // Using an example in https://en.wikipedia.org/wiki/Neighbor_joining.

  using DistanceType = double;

  SymmetricMatrix<DistanceType> distances;

  distances.Set(0, 0, 0);
  distances.Set(0, 1, 5);
  distances.Set(0, 2, 9);
  distances.Set(0, 3, 9);
  distances.Set(0, 4, 8);
  distances.Set(1, 1, 0);
  distances.Set(1, 2, 10);
  distances.Set(1, 3, 10);
  distances.Set(1, 4, 9);
  distances.Set(2, 2, 0);
  distances.Set(2, 3, 8);
  distances.Set(2, 4, 7);
  distances.Set(3, 3, 0);
  distances.Set(3, 4, 3);
  distances.Set(4, 4, 0);

  const Result<DistanceType> result = Execute<DistanceType>(distances, 5);

  ASSERT_EQ(result.Neighbors(0), std::vector<int>{5});
  ASSERT_EQ(result.Neighbors(1), std::vector<int>{5});
  ASSERT_EQ(result.Neighbors(2), std::vector<int>{6});
  ASSERT_EQ(result.Neighbors(3), std::vector<int>{7});
  ASSERT_EQ(result.Neighbors(4), std::vector<int>{7});
  ASSERT_EQ(result.Neighbors(5), (std::vector<int>{0, 1, 6}));
  ASSERT_EQ(result.Neighbors(6), (std::vector<int>{2, 5, 7}));
  ASSERT_EQ(result.Neighbors(7), (std::vector<int>{3, 4, 6}));

  // Distances between adjacent nodes.
  ASSERT_DOUBLE_EQ(result.Distance(0, 5), 2);
  ASSERT_DOUBLE_EQ(result.Distance(1, 5), 3);
  ASSERT_DOUBLE_EQ(result.Distance(2, 6), 4);
  ASSERT_DOUBLE_EQ(result.Distance(3, 7), 2);
  ASSERT_DOUBLE_EQ(result.Distance(4, 7), 1);
  ASSERT_DOUBLE_EQ(result.Distance(5, 6), 3);
  ASSERT_DOUBLE_EQ(result.Distance(6, 7), 2);

  // Distances between non-adjacent nodes.
  ASSERT_DOUBLE_EQ(result.Distance(1, 3), 10);
}

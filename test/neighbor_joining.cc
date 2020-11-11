#include "neighbor_joining.h"

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

  ASSERT_EQ(result.parent.find(0)->second, 5);
  ASSERT_EQ(result.parent.find(1)->second, 5);
  ASSERT_EQ(result.parent.find(2)->second, 6);
  ASSERT_EQ(result.parent.find(3)->second, 7);
  ASSERT_EQ(result.parent.find(4)->second, 7);
  ASSERT_EQ(result.parent.find(5)->second, 6);
  ASSERT_EQ(result.parent.find(6)->second, 7);

  ASSERT_DOUBLE_EQ(result.distances.Get(0, 5), 2);
  ASSERT_DOUBLE_EQ(result.distances.Get(1, 5), 3);
  ASSERT_DOUBLE_EQ(result.distances.Get(2, 6), 4);
  ASSERT_DOUBLE_EQ(result.distances.Get(3, 7), 2);
  ASSERT_DOUBLE_EQ(result.distances.Get(4, 7), 1);
  ASSERT_DOUBLE_EQ(result.distances.Get(5, 6), 3);
  ASSERT_DOUBLE_EQ(result.distances.Get(6, 7), 2);
}

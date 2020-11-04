#include "graph.h"

#include "gtest/gtest.h"

TEST(UnionFind, UniteAndFind) {
  UnionFind uf;
  uf.Unite(1, 2);
  uf.Unite(2, 3);
  uf.Unite(4, 5);

  ASSERT_TRUE(uf.Same(1, 2));
  ASSERT_TRUE(uf.Same(1, 3));
  ASSERT_TRUE(uf.Same(2, 3));
  ASSERT_TRUE(uf.Same(4, 5));

  ASSERT_FALSE(uf.Same(1, 4));
  ASSERT_FALSE(uf.Same(2, 4));
  ASSERT_FALSE(uf.Same(3, 4));
  ASSERT_FALSE(uf.Same(1, 5));
  ASSERT_FALSE(uf.Same(2, 5));
  ASSERT_FALSE(uf.Same(3, 5));
}

TEST(BidirectionalGraph, MST) {
  BidirectionalGraph g;

  for (int i = 1; i <= 4; i++) {
    for (int j = i + 1; j <= 5; j++) {
      if ((i == 1 && j == 2) || (i == 1 && j == 3) || (i == 2 && j == 4) ||
          (i == 2 && j == 5))
        g.AddEdge(i, j, 100);
      else
        g.AddEdge(i, j, 1000);
    }
  }

  const auto [cost, tree] = g.MST(1);

  ASSERT_EQ(cost, 400);
  ASSERT_EQ(tree.Root(), 1);
  ASSERT_EQ(tree.Size(), 5);
  ASSERT_EQ(tree.Parent(2), 1);
  ASSERT_EQ(tree.Parent(3), 1);
  ASSERT_EQ(tree.Parent(4), 2);
  ASSERT_EQ(tree.Parent(5), 2);
}
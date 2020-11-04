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
  g.AddEdge(1, 2, 100);
  g.AddEdge(2, 3, 100);
  g.AddEdge(3, 4, 100);
  g.AddEdge(1, 4, 1);
  const auto [cost, tree] = g.MST(1);
  ASSERT_EQ(cost, 201);
  ASSERT_EQ(tree.Root(), 1);
  ASSERT_EQ(tree.Parent(2), 1);
  ASSERT_EQ(tree.Parent(3), 2);
  ASSERT_EQ(tree.Parent(4), 1);
}
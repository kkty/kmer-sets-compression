#ifndef GRAPH_H_
#define GRAPH_H_

#include <queue>
#include <tuple>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"

class UnionFind {
 public:
  int64_t Find(int64_t x) {
    if (parent_.find(x) == parent_.end()) return x;

    return parent_[x] = Find(parent_[x]);
  }

  void Unite(int64_t x, int64_t y) {
    x = Find(x);
    y = Find(y);

    if (x == y) return;

    parent_[x] = y;
  }

  bool Same(int64_t x, int64_t y) { return Find(x) == Find(y); }

 private:
  absl::flat_hash_map<int64_t, int64_t> parent_;
};

class Tree {
 public:
  Tree(int64_t root, absl::flat_hash_map<int64_t, int64_t> parent,
       absl::flat_hash_map<int64_t, std::vector<int64_t>> children)
      : root_(root),
        parent_(std::move(parent)),
        children_(std::move(children)) {}

  int64_t Root() const { return root_; }

  int64_t Parent(int64_t n) const { return parent_.find(n)->second; }

  std::vector<int64_t> Children(int64_t n) const {
    return children_.find(n)->second;
  }

 private:
  int64_t root_;
  absl::flat_hash_map<int64_t, int64_t> parent_;
  absl::flat_hash_map<int64_t, std::vector<int64_t>> children_;
};

class BidirectionalGraph {
 public:
  void AddEdge(int64_t x, int64_t y, int64_t cost) {
    edges_[x][y] = cost;
    edges_[y][x] = cost;
  }

  std::pair<int64_t, Tree> MST(int64_t root) const {
    std::vector<std::pair<int64_t, int64_t>> edges;
    for (const auto& [from, m] : edges_) {
      for (const auto& [to, _] : m) {
        edges.emplace_back(from, to);
      }
    }

    // get_cost(from, to) is similar to edges_[from][to] but does not break
    // const-ness.
    const auto get_cost = [&](const int64_t from, const int64_t to) {
      const auto& m = edges_.find(from)->second;
      return m.find(to)->second;
    };

    std::sort(edges.begin(), edges.end(),
              [&](const std::pair<int64_t, int64_t>& x,
                  const std::pair<int64_t, int64_t>& y) {
                return get_cost(x.first, x.second) <
                       get_cost(y.first, y.second);
              });

    UnionFind uf;

    int64_t cost = 0;
    absl::flat_hash_map<int64_t, std::vector<int64_t>> selected_edges;

    for (const std::pair<int64_t, int64_t>& edge : edges) {
      if (uf.Same(edge.first, edge.second)) continue;

      cost += get_cost(edge.first, edge.second);

      selected_edges[edge.first].push_back(edge.second);
      selected_edges[edge.second].push_back(edge.first);

      uf.Unite(edge.first, edge.second);
    }

    absl::flat_hash_map<int64_t, int64_t> parent;
    absl::flat_hash_map<int64_t, std::vector<int64_t>> children;

    // BFS to construct Tree.
    {
      std::queue<int64_t> queue;
      absl::flat_hash_set<int64_t> visited;
      queue.push(root);
      visited.insert(root);

      while (!queue.empty()) {
        int64_t n = queue.front();
        queue.pop();

        for (const int64_t to : selected_edges[n]) {
          if (visited.find(to) == visited.end()) {
            visited.insert(to);
            queue.push(to);
            parent[to] = n;
            children[n].push_back(to);
          }
        }
      }
    }

    Tree tree(root, std::move(parent), std::move(children));

    return std::make_pair(cost, tree);
  }

 private:
  absl::flat_hash_map<int64_t, absl::flat_hash_map<int64_t, int64_t>> edges_;
};

#endif
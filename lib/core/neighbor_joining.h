#ifndef NEIGHBOR_JOINING_H_
#define NEIGHBOR_JOINING_H_

#include <algorithm>
#include <limits>
#include <queue>
#include <tuple>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "core/graph.h"

namespace neighbor_joining {

// Symmetric matrix whose items can be empty.
template <typename ValueType>
class SymmetricMatrix {
 public:
  void Set(int i, int j, ValueType distance) {
    if (i > j) std::swap(i, j);

    values_[std::make_pair(i, j)] = distance;
  }

  // Before calling Get(i, j), Set(i, j) or Set(j, i) should be called.
  ValueType Get(int i, int j) const {
    if (i > j) std::swap(i, j);

    return values_.find(std::make_pair(i, j))->second;
  }

  // Finds a minimum element from the ones whose value is set.
  // There should be at least one element in the matrix.
  std::tuple<ValueType, int, int> Min() const {
    std::vector<std::tuple<ValueType, int, int>> candidates;

    for (const auto& [key, value] : values_) {
      candidates.emplace_back(value, key.first, key.second);
    }

    return *std::min_element(candidates.begin(), candidates.end());
  }

 private:
  absl::flat_hash_map<std::pair<int, int>, ValueType> values_;
};

// For each i, j in active, the distance between i and j should be
// defined.
template <typename DistanceType>
SymmetricMatrix<DistanceType> GetQ(
    const SymmetricMatrix<DistanceType>& distances, std::vector<int>& active) {
  SymmetricMatrix<DistanceType> q;

  for (size_t i = 0; i < active.size(); i++) {
    for (size_t j = i + 1; j < active.size(); j++) {
      DistanceType sum1 = [&] {
        DistanceType sum = 0;
        for (size_t k = 0; k < active.size(); k++)
          sum += distances.Get(active[i], active[k]);
        return sum;
      }();

      DistanceType sum2 = [&] {
        DistanceType sum = 0;
        for (size_t k = 0; k < active.size(); k++)
          sum += distances.Get(active[j], active[k]);
        return sum;
      }();

      q.Set(active[i], active[j],
            (active.size() - 2) * distances.Get(active[i], active[j]) - sum1 -
                sum2);
    }
  }

  return q;
}

// Tree structure obtained by neighbor joining.
template <typename DistanceType>
class Result {
 public:
  Result(absl::flat_hash_map<int, std::vector<int>> adjacency_list,
         SymmetricMatrix<DistanceType> distances)
      : adjacency_list_(std::move(adjacency_list)),
        distances_(std::move(distances)) {}

  // Returns a list of nodes.
  std::vector<int> Nodes() const {
    std::vector<int> nodes;
    for (auto it = adjacency_list_.cbegin(); it != adjacency_list_.cend();
         it++) {
      nodes.push_back(it->first);
    }
    return nodes;
  }

  // Returns the neighbors of "node".
  std::vector<int> Neighbors(int node) const {
    std::vector<int> neighbors = adjacency_list_.find(node)->second;
    std::sort(neighbors.begin(), neighbors.end());
    return neighbors;
  }

  // Returns the distances to each node from "root".
  absl::flat_hash_map<int, DistanceType> Distances(int root) const {
    absl::flat_hash_map<int, DistanceType> distances;

    std::queue<int> queue;
    queue.push(root);
    distances[root] = 0;

    while (!queue.empty()) {
      int node = queue.front();
      queue.pop();

      for (int to : adjacency_list_.find(node)->second) {
        if (distances.find(to) == distances.end()) {
          distances[to] = distances[node] + distances_.Get(node, to);
          queue.push(to);
        }
      }
    }

    return distances;
  }

  // Returns the distance from "from" to "to".
  DistanceType Distance(int from, int to) const { return Distances(from)[to]; }

  // Returns the all-pair distances.
  SymmetricMatrix<DistanceType> AllDistances() const {
    SymmetricMatrix<DistanceType> all_distances;

    for (int root : Nodes()) {
      absl::flat_hash_map<int, DistanceType> distances = Distances(root);
      for (auto it = distances.cbegin(); it != distances.cend(); it++) {
        all_distances.Set(root, it->first, it->second);
      }
    }

    return all_distances;
  }

 private:
  absl::flat_hash_map<int, std::vector<int>> adjacency_list_;
  SymmetricMatrix<DistanceType> distances_;
};

// Executes the neighbor joining algorithm.
// For details, please refer to https://en.wikipedia.org/wiki/Neighbor_joining.
//
// n is the number of initial nodes. For 0 <= i, j < n, the distances between
// i and j should be defined.
template <typename DistanceType>
Result<DistanceType> Execute(SymmetricMatrix<DistanceType> distances, int n) {
  absl::flat_hash_map<int, std::vector<int>> adjacency_list;

  if (n == 1) {
    return {adjacency_list, distances};
  }

  std::vector<int> active;
  for (int i = 0; i < n; i++) active.push_back(i);

  while (active.size() > 2) {
    const SymmetricMatrix<DistanceType> q = GetQ(distances, active);

    int f, g;
    std::tie(std::ignore, f, g) = q.Min();
    int u = n++;
    adjacency_list[f].push_back(u);
    adjacency_list[g].push_back(u);
    adjacency_list[u].push_back(f);
    adjacency_list[u].push_back(g);

    // Joins "f" and "g" and creates "u".

    distances.Set(u, u, 0);

    // Calculating the distances of (f, u) and (g, u).
    {
      DistanceType sum1 = 0;
      for (int i : active) sum1 += distances.Get(f, i);

      DistanceType sum2 = 0;
      for (int i : active) sum2 += distances.Get(g, i);

      distances.Set(
          f, u,
          distances.Get(f, g) / 2 + (sum1 - sum2) / (2 * (active.size() - 2)));

      distances.Set(g, u, distances.Get(f, g) - distances.Get(f, u));
    }

    // Calculating the distances of (i, u) for each i.
    for (int i : active) {
      if (i == f || i == g) continue;

      distances.Set(
          u, i,
          (distances.Get(f, i) + distances.Get(g, i) - distances.Get(f, g)) /
              2);
    }

    // Removes "f" and "g" from "active" and adds "u" to it.
    active.erase(std::find(active.begin(), active.end(), f));
    active.erase(std::find(active.begin(), active.end(), g));
    active.push_back(u);
  }

  assert(active.size() == 2);

  adjacency_list[active[0]].push_back(active[1]);
  adjacency_list[active[1]].push_back(active[0]);

  return {adjacency_list, distances};
}

}  // namespace neighbor_joining

#endif
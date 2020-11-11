#ifndef NEIGHBOR_JOINING_H_
#define NEIGHBOR_JOINING_H_

#include <algorithm>
#include <limits>
#include <tuple>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "graph.h"

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

template <typename DistanceType>
struct Result {
  int root;
  absl::flat_hash_map<int, int> parent;
  absl::flat_hash_map<int, std::vector<int>> children;
  SymmetricMatrix<DistanceType> distances;
};

// Executes the neighbor joining algorithm.
// For details, please refer to https://en.wikipedia.org/wiki/Neighbor_joining.
//
// n is the number of initial nodes. For 0 <= i, j < n, the distances between
// i and j should be set in "distances".
template <typename DistanceType>
Result<DistanceType> Execute(SymmetricMatrix<DistanceType> distances, int n) {
  absl::flat_hash_map<int, int> parent;
  absl::flat_hash_map<int, std::vector<int>> children;

  if (n == 1) {
    return {0, parent, children, distances};
  }

  std::vector<int> active;
  for (int i = 0; i < n; i++) active.push_back(i);

  while (active.size() > 2) {
    const SymmetricMatrix<DistanceType> q = GetQ(distances, active);

    int f, g;
    std::tie(std::ignore, f, g) = q.Min();
    int u = n++;
    parent[f] = parent[g] = u;
    children[u] = std::vector<int>{f, g};

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

  // Whichever can be the root. Here, the one with larger index will be the
  // root.
  parent[active[0]] = active[1];
  children[active[1]] = std::vector<int>{active[0]};

  return {active[1], parent, children, distances};
}

}  // namespace neighbor_joining

#endif
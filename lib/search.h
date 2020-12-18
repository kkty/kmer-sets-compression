#ifndef SEARCH_H_
#define SEARCH_H_

#include <algorithm>
#include <cstdint>
#include <mutex>
#include <optional>
#include <queue>
#include <stack>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/random/random.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "core/kmer.h"
#include "core/kmer_set.h"
#include "core/range.h"

// A graph where each node represents a kmer and each edge has its distance.
template <int K>
struct KmerGraph {
  // edges[i] is a list of edges that goes out from node i.
  // The first element of the pair is for node ids.
  // The second element of the pair is for distances.
  std::vector<std::vector<std::pair<std::int64_t, std::int64_t>>> edges;

  std::vector<Kmer<K>> kmers;                      // id to kmer
  absl::flat_hash_map<Kmer<K>, std::int64_t> ids;  // kmer to id
};

// Finds branching kmers (kmers with 2 or more incoming edges or 2 or more
// outgoing edges) in a KmerSet and constructs a graph where each node
// represents a branching kmer.
template <int K, int N, typename KeyType>
KmerGraph<K> ConstructKmerGraph(const KmerSet<K, N, KeyType>& kmer_set,
                                int n_workers) {
  const auto GetNexts = [&](const Kmer<K>& kmer) {
    std::vector<Kmer<K>> v;

    for (const Kmer<K>& next : kmer.Nexts()) {
      if (next != kmer && kmer_set.Contains(next)) v.push_back(next);
    }

    return v;
  };

  const auto GetPrevs = [&](const Kmer<K>& kmer) {
    std::vector<Kmer<K>> v;

    for (const Kmer<K>& prev : kmer.Prevs()) {
      if (prev != kmer && kmer_set.Contains(prev)) v.push_back(prev);
    }

    return v;
  };

  // Finds branching kmers.
  std::vector<Kmer<K>> kmers = kmer_set.Find(
      [&](const Kmer<K>& kmer) {
        return GetNexts(kmer).size() >= 2 || GetPrevs(kmer).size() >= 2;
      },
      n_workers);

  const std::int64_t n_kmers = kmers.size();

  absl::flat_hash_map<Kmer<K>, std::int64_t> ids;  // kmer to id

  for (std::int64_t i = 0; i < n_kmers; i++) ids[kmers[i]] = i;

  std::vector<std::vector<std::pair<std::int64_t, std::int64_t>>> edges(
      n_kmers);
  std::mutex mu;

  {
    boost::asio::thread_pool pool(n_workers);

    for (std::int64_t i = 0; i < n_kmers; i++) {
      boost::asio::post(pool, [&, i] {
        const Kmer<K>& start = kmers[i];
        absl::flat_hash_map<Kmer<K>, std::int64_t> distances;
        distances[start] = 0;

        std::queue<Kmer<K>> queue;
        queue.push(start);

        while (!queue.empty()) {
          const Kmer<K> kmer = queue.front();
          queue.pop();

          for (const Kmer<K>& next : GetNexts(kmer)) {
            if (distances.find(next) == distances.end()) {
              distances[next] = distances[kmer] + 1;

              if (ids.find(next) != ids.end()) {
                std::lock_guard lck(mu);
                edges[ids[start]].emplace_back(ids[next], distances[next]);
              } else {
                queue.push(next);
              }
            }
          }
        }
      });
    }

    pool.join();
  }

  return {edges, kmers, ids};
}

struct SearchResult {
  bool found;
  std::int64_t distance;
  std::int64_t visited_nodes;
};

template <int K>
SearchResult DijkstraSearch(const KmerGraph<K>& g, const Kmer<K>& start,
                            const Kmer<K>& goal) {
  const std::int64_t start_id = g.ids.find(start)->second;
  const std::int64_t goal_id = g.ids.find(goal)->second;

  absl::flat_hash_map<std::int64_t, std::int64_t> distances;
  distances[start_id] = 0;
  std::priority_queue<std::pair<std::int64_t, std::int64_t>,
                      std::vector<std::pair<std::int64_t, std::int64_t>>,
                      std::greater<>>
      queue;
  queue.push({distances[start_id], start_id});

  while (!queue.empty()) {
    int from = queue.top().second;
    queue.pop();

    if (from == goal_id) {
      return {true, distances[from],
              static_cast<std::int64_t>(distances.size())};
    }

    for (const std::pair<std::int64_t, std::int64_t>& edge : g.edges[from]) {
      std::int64_t to, cost;
      std::tie(to, cost) = edge;

      if (distances.find(to) == distances.end() ||
          distances[to] > distances[from] + cost) {
        distances[to] = distances[from] + cost;
        queue.push({distances[to], to});
      }
    }
  }

  return {false, 0, 0};
}

// Calculates the all-pair distances among the branching nodes.
template <int K, int N, typename KeyType>
absl::flat_hash_map<std::pair<Kmer<K>, Kmer<K>>, std::int64_t>
GetDistancesAmongBranchingKmers(const KmerSet<K, N, KeyType>& kmer_set,
                                int n_workers) {
  const auto GetNexts = [&](const Kmer<K>& kmer) {
    std::vector<Kmer<K>> v;
    for (const Kmer<K>& next : kmer.Nexts()) {
      if (kmer_set.Contains(next)) {
        v.push_back(kmer);
      }
    }
    return v;
  };

  const auto GetPrevs = [&](const Kmer<K>& kmer) {
    std::vector<Kmer<K>> v;
    for (const Kmer<K>& prev : kmer.Prevs()) {
      if (kmer_set.Contains(prev)) {
        v.push_back(kmer);
      }
    }
    return v;
  };

  absl::flat_hash_set<Kmer<K>> branching_kmers;

  {
    const std::vector<Kmer<K>> v = kmer_set.Find(
        [&](const Kmer<K>& kmer) {
          return GetNexts(kmer).size() >= 2 || GetPrevs(kmer).size() >= 2;
        },
        n_workers);

    branching_kmers.reserve(v.size());
    branching_kmers.insert(v.begin(), v.end());
  }

  absl::flat_hash_map<Kmer<K>, absl::flat_hash_map<Kmer<K>, std::int64_t>>
      distances;

  for (const Kmer<K>& kmer : branching_kmers) {
    for (const Kmer<K>& next : kmer.Nexts()) {
      if (!kmer_set.Contains(next)) continue;

      Kmer<K> current = next;
      int distance = 1;
      while (true) {
        const std::vector<Kmer<K>> nexts = GetNexts(current);
        const std::vector<Kmer<K>> prevs = GetPrevs(current);

        if (nexts.size() >= 2 || prevs.size() >= 2) {
          distances[kmer][current] = distance;
          break;
        }

        if (nexts.empty()) break;

        current = nexts.front();
        distance += 1;
      }
    }
  }

  for (auto it = distances.begin(); it != distances.end(); ++it) {
    const Kmer<K> start = it->first;

    // Distances from "start".
    absl::flat_hash_map<Kmer<K>, std::int64_t>& d = distances[start];

    using Item = std::pair<Kmer<K>, std::int64_t>;
    const auto cmp = [](const Item& lhs, const Item& rhs) {
      return lhs.second > rhs.second;
    };
    std::priority_queue<Item, std::vector<Item>, decltype(cmp)> pq(cmp);

    for (std::pair<const Kmer<K>, std::int64_t>& p : d) {
      pq.emplace(p.first, p.second);
    }

    while (!pq.empty()) {
      const Kmer<K> from = pq.top().first;
      pq.pop();

      for (std::pair<const Kmer<K>, std::int64_t>& p : distances[from]) {
        const Kmer<K> to = p.first;
        const std::int64_t distance = p.second;

        if (d.find(to) == d.end() || d[to] > d[from] + distance) {
          d[to] = d[from] + distance;
          pq.emplace(to, d[to]);
        }
      }
    }
  }

  absl::flat_hash_map<std::pair<Kmer<K>, Kmer<K>>, std::int64_t> flat;
  for (std::pair<const Kmer<K>, absl::flat_hash_map<Kmer<K>, std::int64_t>>&
           p1 : distances) {
    for (std::pair<const Kmer<K>, std::int64_t>& p2 : p1.second) {
      flat[std::make_pair(p1.first, p2.first)] = p2.second;
    }
  }

  return flat;
}

// Executes A* search using the all-pair shortest distances among branching
// K2-kmers.
template <int K1, int K2, int N2, typename KeyType2>
SearchResult AStarSearch(
    const KmerGraph<K1>& g, const Kmer<K1>& start, const Kmer<K1>& goal,
    const KmerSet<K2, N2, KeyType2>& kmer_set,
    const absl::flat_hash_map<std::pair<Kmer<K2>, Kmer<K2>>, std::int64_t>&
        distances) {
  const auto GetNexts = [&](const Kmer<K2>& kmer) {
    std::vector<Kmer<K2>> v;
    for (const Kmer<K2>& next : kmer.Nexts()) {
      if (kmer_set.Contains(next)) {
        v.push_back(kmer);
      }
    }
    return v;
  };

  const auto GetPrevs = [&](const Kmer<K2>& kmer) {
    std::vector<Kmer<K2>> v;
    for (const Kmer<K2>& prev : kmer.Prevs()) {
      if (kmer_set.Contains(prev)) {
        v.push_back(kmer);
      }
    }
    return v;
  };

  const auto h = [&](int i) {
    std::string s = g.kmers[i].String();
    std::string s_goal = goal.String();

    std::int64_t estimate = K1;

    // Considers the overlap of the suffix of s and the prefix of s_goal.
    for (int j = 0; j < K1; j++) {
      // If the (K1-j)-suffix of s is equal to the (K1-j)-prefix of s_goal,
      // there may be a path whose distance is j that goes from s to s_goal.
      if (s.substr(j, K1 - j) == s_goal.substr(0, K1 - j)) {
        estimate = j;
        break;
      }
    }

    for (int j = 0; j < K1 - K2 + 1; j++) {
      Kmer<K2> from = Kmer<K2>(s.substr(j, K2));
      Kmer<K2> to = Kmer<K2>(s_goal.substr(j, K2));

      // Calculates the distance from "from" to "to" if possible.

      int distance = 0;

      // Moves "from" so that it becomes a branching node.
      {
        bool is_valid = true;

        while (true) {
          if (GetNexts(from).size() >= 2 || GetPrevs(from).size() >= 2) {
            break;
          }

          const std::vector<Kmer<K2>> nexts = GetNexts(from);

          if (nexts.empty()) {
            is_valid = false;
            break;
          }

          from = nexts.front();
          distance += 1;
        }

        if (!is_valid) continue;
      }

      // Moves "to" so that it becomes a branching node.
      {
        bool is_valid = true;

        while (true) {
          if (GetNexts(to).size() >= 2 || GetPrevs(to).size() >= 2) {
            break;
          }

          const std::vector<Kmer<K2>> prevs = GetPrevs(to);

          if (prevs.empty()) {
            is_valid = false;
            break;
          }

          to = prevs.front();
          distance += 1;
        }

        if (!is_valid) continue;
      }

      const auto it = distances.find(std::make_pair(from, to));

      if (it != distances.end()) {
        estimate = std::max(estimate, distance + it->second);
      }
    }

    return estimate;
  };

  const std::int64_t start_id = g.ids.find(start)->second;
  const std::int64_t goal_id = g.ids.find(goal)->second;

  // f[i] will be distance(start_id, i) + h(i)
  absl::flat_hash_map<std::int64_t, std::int64_t> f;
  f[start_id] = h(start_id);
  std::priority_queue<std::pair<std::int64_t, std::int64_t>,
                      std::vector<std::pair<std::int64_t, std::int64_t>>,
                      std::greater<>>
      queue;
  queue.emplace(f[start_id], start_id);

  while (!queue.empty()) {
    int from = queue.top().second;
    queue.pop();

    if (from == goal_id) {
      return {true, f[from], static_cast<std::int64_t>(f.size())};
    }

    for (const std::pair<std::int64_t, std::int64_t>& edge : g.edges[from]) {
      std::int64_t to, cost;
      std::tie(to, cost) = edge;

      if (f.find(to) == f.end() || f[to] > f[from] - h(from) + h(to) + cost) {
        f[to] = f[from] - h(from) + h(to) + cost;
        queue.push({f[to], to});
      }
    }
  }

  return {false, 0, 0};
}

template <int K>
SearchResult AStarSearch(const KmerGraph<K>& g, const Kmer<K>& start,
                         const Kmer<K>& goal) {
  const auto h = [&](int i) {
    const std::string s = g.kmers[i].String();
    const std::string s_goal = goal.String();

    for (int j = 0; j < K; j++) {
      // If the (K-j)-suffix of s is equal to the (K-j)-prefix of s_goal, there
      // may be a path whose distance is j that goes from s to s_goal.
      if (s.substr(j, K - j) == s_goal.substr(0, K - j)) return j;
    }

    return K;
  };

  const std::int64_t start_id = g.ids.find(start)->second;
  const std::int64_t goal_id = g.ids.find(goal)->second;

  // f[i] will be distance(start_id, i) + h(i)
  absl::flat_hash_map<std::int64_t, std::int64_t> f;
  f[start_id] = h(start_id);
  std::priority_queue<std::pair<std::int64_t, std::int64_t>,
                      std::vector<std::pair<std::int64_t, std::int64_t>>,
                      std::greater<>>
      queue;
  queue.push({f[start_id], start_id});

  while (!queue.empty()) {
    int from = queue.top().second;
    queue.pop();

    if (from == goal_id) {
      return {true, f[from], static_cast<int64_t>(f.size())};
    }

    for (const std::pair<std::int64_t, std::int64_t>& edge : g.edges[from]) {
      std::int64_t to, cost;
      std::tie(to, cost) = edge;

      if (f.find(to) == f.end() || f[to] > f[from] - h(from) + h(to) + cost) {
        f[to] = f[from] - h(from) + h(to) + cost;
        queue.push({f[to], to});
      }
    }
  }

  return {false, 0, 0};
}

#endif
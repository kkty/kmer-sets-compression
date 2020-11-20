#ifndef SEARCH_H_
#define SEARCH_H_

#include <mutex>
#include <queue>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/random/random.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "core/graph.h"
#include "core/kmer.h"
#include "core/kmer_set.h"
#include "core/range.h"

// A graph where each node represents a kmer and each edge has its distance.
template <int K>
struct KmerGraph {
  // edges[i] is a list of edges that goes out from node i.
  // The first element of the pair is for node ids.
  // The second element of the pair is for distances.
  std::vector<std::vector<std::pair<int64_t, int64_t>>> edges;

  std::vector<Kmer<K>> kmers;                 // id to kmer
  absl::flat_hash_map<Kmer<K>, int64_t> ids;  // kmer to id
};

// Randomly finds p such that there is a path from p.first to p.second in
// g and the path consists of n or more edges.
template <int K>
std::pair<Kmer<K>, Kmer<K>> SampleConnectedPair(const KmerGraph<K>& g, int n) {
  absl::InsecureBitGen bitgen;

  while (true) {
    int64_t start = absl::Uniform(bitgen, (int64_t)0, (int64_t)g.edges.size());

    // DFS from start.

    std::stack<int64_t> stack;
    absl::flat_hash_map<int64_t, int64_t> distances;

    stack.push(start);
    distances[start] = 0;

    while (!stack.empty()) {
      int64_t current = stack.top();
      stack.pop();

      for (const std::pair<int64_t, int64_t>& edge : g.edges[current]) {
        int64_t to;
        std::tie(to, std::ignore) = edge;

        if (distances.find(to) == distances.end()) {
          distances[to] = distances[current] + 1;
          if (distances[to] >= n) return {g.kmers[start], g.kmers[to]};
          stack.push(to);
        }
      }
    }
  }
}

template <int K, typename KeyType>
KmerGraph<K> ConstructKmerGraph(const KmerSet<K, KeyType>& kmer_set,
                                int n_workers) {
  const auto get_nexts = [&](const Kmer<K>& kmer) {
    std::vector<Kmer<K>> v;

    for (const Kmer<K>& next : kmer.Nexts()) {
      if (next != kmer && kmer_set.Contains(next)) v.push_back(next);
    }

    return v;
  };

  const auto get_prevs = [&](const Kmer<K>& kmer) {
    std::vector<Kmer<K>> v;

    for (const Kmer<K>& prev : kmer.Prevs()) {
      if (prev != kmer && kmer_set.Contains(prev)) v.push_back(prev);
    }

    return v;
  };

  // Finds branching kmers.
  std::vector<Kmer<K>> kmers = kmer_set.Find(
      [&](const Kmer<K>& kmer) {
        return get_nexts(kmer).size() >= 2 || get_prevs(kmer).size() >= 2;
      },
      n_workers);

  const int64_t n_kmers = kmers.size();

  absl::flat_hash_map<Kmer<K>, int64_t> ids;  // kmer to id

  for (int64_t i = 0; i < n_kmers; i++) ids[kmers[i]] = i;

  std::vector<std::vector<std::pair<int64_t, int64_t>>> edges(n_kmers);
  std::mutex mu;

  {
    boost::asio::thread_pool pool(n_workers);

    for (int64_t i = 0; i < n_kmers; i++) {
      boost::asio::post(pool, [&, i] {
        const Kmer<K>& start = kmers[i];
        absl::flat_hash_map<Kmer<K>, int64_t> distances;
        distances[start] = 0;

        std::queue<Kmer<K>> queue;
        queue.push(start);

        while (!queue.empty()) {
          const Kmer<K> kmer = queue.front();
          queue.pop();

          for (const Kmer<K>& next : get_nexts(kmer)) {
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
  int64_t distance;
  int64_t visited_nodes;
};

template <int K>
SearchResult DijkstraSearch(const KmerGraph<K>& g, const Kmer<K>& start,
                            const Kmer<K>& goal) {
  const int64_t start_id = g.ids.find(start)->second;
  const int64_t goal_id = g.ids.find(goal)->second;

  absl::flat_hash_map<int64_t, int64_t> distances;
  distances[start_id] = 0;
  std::priority_queue<std::pair<int64_t, int64_t>,
                      std::vector<std::pair<int64_t, int64_t>>,
                      std::greater<std::pair<int64_t, int64_t>>>
      queue;
  queue.push({distances[start_id], start_id});

  while (!queue.empty()) {
    int from = queue.top().second;
    queue.pop();

    if (from == goal_id) {
      return {true, distances[from], (int64_t)distances.size()};
    }

    for (const std::pair<int64_t, int64_t>& edge : g.edges[from]) {
      int64_t to, cost;
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

template <int K, int L>
SearchResult AStarSearch(const KmerGraph<K>& g, const Kmer<K>& start,
                         const Kmer<K>& goal) {
  const std::vector<Kmer<K>>& kmers = g.kmers;
  const absl::flat_hash_map<Kmer<K>, int64_t>& ids = g.ids;

  // Returns the set of l-mers present in "kmer".
  const auto get_lmers = [](const Kmer<K>& kmer) {
    absl::flat_hash_set<Kmer<L>> lmers;

    for (int i = 0; i < K - L + 1; i++) {
      Kmer<L> lmer;

      for (int j = 0; j < L; j++) {
        lmer.Set(j, kmer.Get(i + j));
      }

      lmers.insert(lmer);
    }

    return lmers;
  };

  // The set of l-mers present in "goal".
  // It is used to calculate h(n).
  const auto lmers_goal = get_lmers(goal);

  // h(i) gives an approximate of distance(i, goal_id).
  const auto h = [&](int i) {
    const auto lmers = get_lmers(kmers[i]);

    int distance = 0;

    for (const auto& lmer_goal : lmers_goal) {
      if (lmers.find(lmer_goal) == lmers.end()) distance += 1;
    }

    for (const auto& lmer : lmers) {
      if (lmers_goal.find(lmer) == lmers_goal.end()) distance += 1;
    }

    return distance / 2;
  };

  const int64_t start_id = ids.find(start)->second;
  const int64_t goal_id = ids.find(goal)->second;

  // f[i] will be distance(start_id, i) + h(i)
  absl::flat_hash_map<int64_t, int64_t> f;
  f[start_id] = h(start_id);
  std::priority_queue<std::pair<int64_t, int64_t>,
                      std::vector<std::pair<int64_t, int64_t>>,
                      std::greater<std::pair<int64_t, int64_t>>>
      queue;
  queue.push({f[start_id], start_id});

  while (!queue.empty()) {
    int from = queue.top().second;
    queue.pop();

    if (from == goal_id) {
      return {true, f[from], (int64_t)f.size()};
    }

    for (const std::pair<int64_t, int64_t>& edge : g.edges[from]) {
      int64_t to, cost;
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
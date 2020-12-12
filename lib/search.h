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
      return {true, distances[from], (std::int64_t)distances.size()};
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

// Calculates all-pair shortest distances in dBG.
// If GetAllPairDistances(...)[{kmer1, kmer2}] = i, the distance from kmer1 to
// kmer2 is i.
template <int K, int N, typename KeyType>
absl::flat_hash_map<std::pair<Kmer<K>, Kmer<K>>, std::int64_t>
GetAllPairDistances(const KmerSet<K, N, KeyType>& kmer_set,
                    std::optional<int> max_distance, int n_workers) {
  const std::vector<Kmer<K>> kmers = kmer_set.Find(n_workers);

  absl::flat_hash_map<std::pair<Kmer<K>, Kmer<K>>, std::int64_t> distances;

  std::vector<std::thread> threads;
  std::mutex mu;

  for (const Range& range : Range(0, kmers.size()).Split(n_workers)) {
    threads.emplace_back([&, range] {
      absl::flat_hash_map<std::pair<Kmer<K>, Kmer<K>>, std::int64_t>
          buf_distances;

      for (std::int64_t i : range) {
        const Kmer<K>& start = kmers[i];

        absl::flat_hash_map<Kmer<K>, std::int64_t> d;
        d[start] = 0;

        std::queue<Kmer<K>> queue;
        queue.push(start);

        while (!queue.empty()) {
          const Kmer<K> current = queue.front();
          queue.pop();

          for (const Kmer<K>& next : current.Nexts()) {
            // If "next" is not in kmer_set, skips it.
            if (!kmer_set.Contains(next)) continue;

            // If "next" has been visited already, skips it.
            if (d.find(next) != d.end()) continue;

            // If "next" is too far away from "start", skips it.
            if (max_distance.has_value() &&
                d[current] + 1 > max_distance.value())
              continue;

            d[next] = d[current] + 1;
            queue.push(next);
          }
        }

        for (const std::pair<const Kmer<K>, std::int64_t>& p : d) {
          buf_distances[{start, p.first}] = p.second;
        }
      }

      std::lock_guard lck(mu);
      distances.insert(buf_distances.begin(), buf_distances.end());
    });
  }

  for (std::thread& t : threads) t.join();

  return distances;
}

// Executes A* search using the all-pair shortest distances for the dBG
// represented by L-mers.
template <int K, int L>
SearchResult AStarSearch(const KmerGraph<K>& g, const Kmer<K>& start,
                         const Kmer<K>& goal,
                         const absl::flat_hash_map<std::pair<Kmer<L>, Kmer<L>>,
                                                   std::int64_t>& distances) {
  const auto h = [&](int i) {
    std::int64_t max = 0;

    std::string s = g.kmers[i].String();
    std::string s_goal = goal.String();

    for (int j = 0; j < K - L + 1; j++) {
      max =
          std::max(max, distances
                            .find(std::make_pair(Kmer<L>(s.substr(j, L)),
                                                 Kmer<L>(s_goal.substr(j, L))))
                            ->second

          );
    }

    return max;
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
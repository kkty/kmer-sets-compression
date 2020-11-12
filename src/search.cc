#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <optional>
#include <queue>
#include <string>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/status/status.h"
#include "kmer.h"
#include "kmer_counter.h"
#include "kmer_set.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(bool, canonical, false, "count canonical k-mers");
ABSL_FLAG(int, workers, 1, "number of workers");
ABSL_FLAG(int, seed, 0, "seed for rand()");
ABSL_FLAG(int, cutoff, 1, "cut off threshold");
ABSL_FLAG(int, n, 1, "number of iterations");

template <int K>
struct BFSResult {
  absl::flat_hash_map<Kmer<K>, int64_t> distances;
  int64_t reached_nodes;
  int64_t max_distance;
  double avg_distance;
};

// BFS from "start". "limit" is used to specify the search space.
template <int K, typename KeyType>
BFSResult<K> BFS(const KmerSet<K, KeyType>& kmer_set, const Kmer<K>& start,
                 const int limit) {
  if (!kmer_set.Contains(start)) {
    spdlog::error("node not found: {}", start.String());
    std::exit(1);
  }

  absl::flat_hash_map<Kmer<K>, int64_t> distances;
  distances[start] = 0;

  const auto is_visited = [&](const Kmer<K>& node) {
    return distances.find(node) != distances.end();
  };

  std::queue<Kmer<K>> queue;
  queue.push(start);

  int64_t max_distance = 0;
  int64_t sum_distance = 0;

  while (!queue.empty()) {
    const Kmer<K> from = queue.front();
    queue.pop();

    for (const auto& to : from.Nexts()) {
      if (kmer_set.Contains(to) && !is_visited(to) && distances[from] < limit) {
        distances[to] = distances[from] + 1;
        queue.push(to);
        max_distance = std::max(max_distance, distances[to]);
        sum_distance += distances[to];
      }
    }
  }

  int64_t visited_nodes = distances.size();

  return {
      distances,
      visited_nodes,
      max_distance,
      (double)sum_distance / distances.size(),
  };

  spdlog::info("reached {} nodes", distances.size());
  spdlog::info("max_distance = {}", max_distance);

  double avg_distance = (double)sum_distance / distances.size();
  spdlog::info("avg_distance = {:.1f}", avg_distance);
}

template <int K>
struct AStarSearchResult {
  std::vector<Kmer<K>> path;
  int visited_nodes;
};

// A* search from "start" to "goal".
// Variable names, etc. follow https://ja.wikipedia.org/wiki/A*.
template <int K, typename KeyType, int L>
std::optional<AStarSearchResult<K>> AStarSearch(
    const KmerSet<K, KeyType>& kmer_set, const Kmer<K>& start,
    const Kmer<K>& goal) {
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

  // h(n) gives an approximate of distance(n, goal).
  const auto h = [&](const Kmer<K>& n) {
    const auto lmers_n = get_lmers(n);

    int distance = 0;

    for (const auto& lmer_goal : lmers_goal) {
      if (lmers_n.find(lmer_goal) == lmers_n.end()) distance += 1;
    }

    for (const auto& lmer_n : lmers_n) {
      if (lmers_goal.find(lmer_n) == lmers_goal.end()) distance += 1;
    }

    return distance;
  };

  // f[n] = h(n) + g(n).
  absl::flat_hash_map<Kmer<K>, int> f;

  // g(n) gives an approximate of distance(start, n).
  const auto g = [&](const Kmer<K>& n) { return f[n] - h(n); };

  const auto cmp = [&](const Kmer<K>& lhs, const Kmer<K>& rhs) {
    return f[lhs] > f[rhs];
  };

  std::priority_queue<Kmer<K>, std::vector<Kmer<K>>, decltype(cmp)> open(cmp);

  absl::flat_hash_set<Kmer<K>> close;

  // This will be used to re-construct the path from "start" to "goal".
  absl::flat_hash_map<Kmer<K>, Kmer<K>> parents;

  f[start] = h(start);
  open.emplace(start);

  while (!open.empty()) {
    spdlog::debug("open.size() = {}", open.size());
    spdlog::debug("close.size() = {}", close.size());

    const Kmer<K> n = open.top();
    open.pop();

    spdlog::debug("n = {}", n.String());

    if (goal == n) {
      // Re-constrcut the path.
      std::vector<Kmer<K>> path;
      path.push_back(n);
      while (parents.find(path.back()) != parents.end()) {
        path.push_back(parents[path.back()]);
      }
      std::reverse(path.begin(), path.end());
      int visited_nodes = f.size();
      spdlog::info("path.size() = {}", path.size());
      return {{path, visited_nodes}};
    }

    for (const auto& to : n.Nexts()) {
      if (kmer_set.Contains(to)) {
        const int new_f = g(n) + 1 + h(to);

        spdlog::debug("to = {}", to.String());
        spdlog::debug("new_f = {}", new_f);

        const auto it = f.find(to);
        if (it == f.end() || new_f < it->second) {
          f[to] = new_f;
          open.push(to);
          parents[to] = n;
        }
      }
    }
  }

  // Paths were not found.
  return {};
}

int main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);

  if (absl::GetFlag(FLAGS_debug)) {
    spdlog::set_level(spdlog::level::debug);
  }

  std::srand(absl::GetFlag(FLAGS_seed));
  std::ios_base::sync_with_stdio(false);

  const int K = 15;
  const int n_workers = absl::GetFlag(FLAGS_workers);
  using KeyType = uint16_t;

  spdlog::info("constructing kmer_counter");

  KmerCounter<K, KeyType> kmer_counter;

  {
    const absl::Status status = kmer_counter.FromFASTQ(
        std::cin, absl::GetFlag(FLAGS_canonical), n_workers);

    if (!status.ok()) {
      spdlog::error("failed to parse FASTQ file");
      std::exit(1);
    }
  }

  spdlog::info("constructed kmer_counter");
  spdlog::info("kmer_counter.Size() = {}", kmer_counter.Size());

  spdlog::info("constructing kmer_set");

  KmerSet<K, KeyType> kmer_set;
  int64_t cutoff_count;
  std::tie(kmer_set, cutoff_count) =
      kmer_counter.ToSet(absl::GetFlag(FLAGS_cutoff), n_workers);

  spdlog::info("constructed kmer_set");
  spdlog::info("cutoff_count = {}", cutoff_count);

  for (int i = 0; i < absl::GetFlag(FLAGS_n); i++) {
    const auto start = kmer_set.Sample();

    spdlog::info("searching: start = {}", start.String());

    const BFSResult bfs_result = BFS(kmer_set, start, 3000);

    spdlog::info("searched: start = {}", start.String());
    spdlog::info("reached_nodes = {}", bfs_result.reached_nodes);
    spdlog::info("max_distance = {}", bfs_result.max_distance);
    spdlog::info("avg_distance = {:.1f}", bfs_result.avg_distance);

    std::map<int, int> distance_count;

    for (const auto& p : bfs_result.distances) {
      distance_count[p.second] += 1;
      spdlog::debug("distances[{}] = {}", p.first.String(), p.second);
    }

    // distance_count_sum[i] = distance_count[0] + ... + distance_count[i];
    std::map<int, int> distance_count_sum;

    {
      int sum = 0;

      for (const auto& p : distance_count) {
        sum += p.second;
        distance_count_sum[p.first] = sum;
      }
    }

    // Randomly selects the goal which was found in the BFS.
    const Kmer<K> goal = [&] {
      std::vector<Kmer<K>> nodes;

      for (const auto& p : bfs_result.distances) {
        nodes.push_back(p.first);
      }

      return nodes[rand() % nodes.size()];
    }();

    spdlog::info("goal = {}", goal.String());
    spdlog::info("distances[{}] = {}", goal.String(),
                 bfs_result.distances.find(goal)->second);
    spdlog::info(
        "distance_count_sum[{}] = {}",
        bfs_result.distances.find(goal)->second - 1,
        distance_count_sum[bfs_result.distances.find(goal)->second - 1]);
    spdlog::info("distance_count_sum[{}] = {}",
                 bfs_result.distances.find(goal)->second,
                 distance_count_sum[bfs_result.distances.find(goal)->second]);

    for (int L = 1; L <= 9; L++) {
      spdlog::info("starting A* search from {} to {}: L = {}", start.String(),
                   goal.String(), L);

      const auto a_star_search_result = [&] {
        if (L == 1) return AStarSearch<K, KeyType, 1>(kmer_set, start, goal);
        if (L == 2) return AStarSearch<K, KeyType, 2>(kmer_set, start, goal);
        if (L == 3) return AStarSearch<K, KeyType, 3>(kmer_set, start, goal);
        if (L == 4) return AStarSearch<K, KeyType, 4>(kmer_set, start, goal);
        if (L == 5) return AStarSearch<K, KeyType, 5>(kmer_set, start, goal);
        if (L == 6) return AStarSearch<K, KeyType, 6>(kmer_set, start, goal);
        if (L == 7) return AStarSearch<K, KeyType, 7>(kmer_set, start, goal);
        if (L == 8) return AStarSearch<K, KeyType, 8>(kmer_set, start, goal);
        return AStarSearch<K, KeyType, 9>(kmer_set, start, goal);
      }();

      if (!a_star_search_result) {
        spdlog::error("A* search failed");
      } else {
        spdlog::info("finished A* search from {} to {}", start.String(),
                     goal.String());
        spdlog::info("visited_nodes = {}",
                     (*a_star_search_result).visited_nodes);
        spdlog::info("distance = {}", (*a_star_search_result).path.size() - 1);
      }
    }
  }
}

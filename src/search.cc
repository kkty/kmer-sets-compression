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
#include "kmer.h"
#include "kmer_counter.h"
#include "kmer_set.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(bool, dump, false, "enable dumping kmers");

template <int K>
struct BFSResult {
  absl::flat_hash_map<Kmer<K>, int64_t> distances;
  int64_t reached_nodes;
  int64_t max_distance;
  double avg_distance;
};

// BFS from "start". "limit" is used to specify the search space.
template <int K, int B>
BFSResult<K> BFS(const KmerSet<K, B>& kmer_set, const Kmer<K>& start,
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
template <int K, int B, int L>
std::optional<AStarSearchResult<K>> AStarSearch(const KmerSet<K, B>& kmer_set,
                                                const Kmer<K>& start,
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

  std::srand(std::time(nullptr));
  std::ios_base::sync_with_stdio(false);

  const int K = 31;
  const int B = 6;

  spdlog::info("constructing kmer_counter");

  KmerCounter<K, B> kmer_counter;
  kmer_counter.FromFASTQ(std::cin);

  spdlog::info("constructed kmer_counter");
  spdlog::info("kmer_counter.Size() = {}", kmer_counter.Size());

  spdlog::info("constructing kmer_set");

  const auto [kmer_set, cut_off_count] = kmer_counter.Set(2);

  spdlog::info("constructed kmer_set");
  spdlog::info("cut_off_count = {}", cut_off_count);

  const auto start = kmer_set.Sample();

  spdlog::info("searching: start = {}", start.String());

  const BFSResult bfs_result = BFS(kmer_set, start, 1000);

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
  spdlog::info("distance_count_sum[{}] = {}",
               bfs_result.distances.find(goal)->second,
               distance_count_sum[bfs_result.distances.find(goal)->second]);
  spdlog::info("starting A* search from {} to {}", start.String(),
               goal.String());

  const auto a_star_search_result = AStarSearch<K, B, 5>(kmer_set, start, goal);

  if (!a_star_search_result) {
    spdlog::error("A* search failed");
  } else {
    spdlog::info("finished A* search from {} to {}", start.String(),
                 goal.String());
    spdlog::info("visited_nodes = {}", (*a_star_search_result).visited_nodes);
    spdlog::info("distance = {}", (*a_star_search_result).path.size() - 1);
  }

  if (absl::GetFlag(FLAGS_dump)) {
    spdlog::info("dumping kmer_set");
    std::ofstream kmer_dump_file;
    kmer_dump_file.open("kmer_dump.txt");
    kmer_set.Dump(kmer_dump_file);
    kmer_dump_file.close();
    spdlog::info("dumped kmer_set");
  }
}
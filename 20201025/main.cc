#include <algorithm>
#include <array>
#include <atomic>
#include <bitset>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/str_split.h"
#include "absl/time/time.h"
#include "omp.h"
#include "spdlog/spdlog.h"

template <int K>
std::bitset<K * 2> GetCompactKmer(const std::string& kmer_s) {
  std::bitset<K * 2> kmer_b;

  // Sets the ith bit from the left.
  const auto set = [&kmer_b](int i) { kmer_b.set(K * 2 - 1 - i); };

  for (int j = 0; j < K; j++) {
    switch (kmer_s[j]) {
      case 'A':
        break;
      case 'C':
        set(j * 2);
        break;
      case 'G':
        set(j * 2 + 1);
        break;
      case 'T':
        set(j * 2);
        set(j * 2 + 1);
        break;
      default:
        spdlog::error("unknown character: {}", kmer_s[j]);
        std::exit(1);
    }
  }

  return kmer_b;
}

// Use (1 << B) buckets internally.
template <int K, int B>
class KmerSet {
 public:
  KmerSet(std::istream& is, const int64_t lines_approximate,
          const int64_t kmers_approximate) {
    std::atomic_bool done = false;
    std::vector<std::thread> threads;
    std::queue<std::string> queue;
    std::mutex queue_mutex;
    std::array<std::mutex, 1 << B> buckets_mutexes;

    for (auto& bucket : buckets_) {
      bucket.reserve(kmers_approximate >> B);
    }

    std::vector<std::string> lines;
    lines.reserve(lines_approximate);

    {
      std::array<std::string, 4> buf;

      // Reads 4 lines from "is" and pushes to "buf".
      auto consume = [&is, &buf]() {
        for (int i = 0; i < 4; i++) {
          if (!std::getline(is, buf[i])) return false;
        }

        return true;
      };

      while (consume() && buf[0].length() && buf[0][0] == '@') {
        lines.push_back(std::move(buf[1]));
      }
    }

    const int64_t n_lines = lines.size();

#pragma omp parallel for
    for (int i = 0; i < n_lines; i++) {
      const auto& line = lines[i];
      std::vector<std::string> fragments = absl::StrSplit(line, "N");
      for (const auto& fragment : fragments) {
        for (int i = 0; i + K <= fragment.length(); i++) {
          const std::string kmer_s = fragment.substr(i, K);
          std::bitset<K* 2> kmer_b = GetCompactKmer<K>(kmer_s);

          int bucket;
          std::bitset<K * 2 - B> key;
          std::tie(bucket, key) = GetBucketAndKey(kmer_b);

          std::lock_guard lck{buckets_mutexes[bucket]};
          buckets_[bucket].insert(key);
        }
      }
    }
  }

  // Loads from a dumped data.
  KmerSet(std::istream& is) {
    int k, b;
    is >> k >> b;
    if (k != K || b != B) {
      spdlog::error("incompatible data");
      return;
    }

    for (int i = 0; i < (1 << B); i++) {
      int64_t size;
      is >> size;
      buckets_[i].reserve(size);

      for (int64_t j = 0; j < size; j++) {
        std::string kmer_s;
        is >> kmer_s;
        std::bitset<K * 2 - B> kmer_b(kmer_s);
        buckets_[i].insert(std::move(kmer_b));
      }
    }
  }

  void Dump(std::ostream& os) const {
    os << K << ' ' << B << '\n';

    for (int i = 0; i < (1 << B); i++) {
      os << buckets_[i].size() << '\n';
      for (const auto& key : buckets_[i]) {
        os << key.to_string() << '\n';
      }
    }
  }

  int Size() const {
    int sum = 0;
    for (const auto& bucket : buckets_) {
      sum += bucket.size();
    }
    return sum;
  }

  bool Contains(const std::bitset<K * 2>& kmer) const {
    int bucket;
    std::bitset<K * 2 - B> key;
    std::tie(bucket, key) = GetBucketAndKey(kmer);
    return buckets_[bucket].find(key) != buckets_[bucket].end();
  }

  std::bitset<K * 2> Sample() const {
    const int bucket = rand() % (1 << B);
    const int64_t n = rand() % buckets_[bucket].size();

    auto it = buckets_[bucket].begin();
    for (int64_t i = 0; i < n; i++) it++;

    const auto key = *it;

    std::bitset<K * 2> kmer;

    for (int i = 0; i < K * 2 - B; i++) kmer.set(i, key[i]);
    {
      std::bitset<B> bucket_b(bucket);
      for (int i = 0; i < B; i++) kmer.set(i + K * 2 - B, bucket_b[i]);
    }

    return kmer;
  }

 private:
  std::array<absl::flat_hash_set<std::bitset<K * 2 - B>>, 1 << B> buckets_;

  static std::pair<int, std::bitset<K * 2 - B>> GetBucketAndKey(
      const std::bitset<K * 2>& kmer) {
    const int key_bits = K * 2 - B;
    const int bucket = kmer.to_ullong() >> key_bits;
    const std::bitset<K * 2 - B> key{kmer.to_ullong() % (1ull << key_bits)};
    return {bucket, key};
  }
};

template <int K>
struct BFSResult {
  absl::flat_hash_map<std::bitset<K * 2>, int64_t> distances;
  int64_t reached_nodes;
  int64_t max_distance;
  double avg_distance;
};

// BFS from "start". "limit" is used to specify the search space.
template <int K, int B>
BFSResult<K> BFS(const KmerSet<K, B>& kmer_set, const std::bitset<K * 2>& start,
                 const int limit) {
  if (!kmer_set.Contains(start)) {
    spdlog::error("node not found: {}", start.to_string());
    std::exit(1);
  }

  absl::flat_hash_map<std::bitset<K * 2>, int64_t> distances;
  distances[start] = 0;

  const auto is_visited = [&](const std::bitset<K * 2>& node) {
    return distances.find(node) != distances.end();
  };

  std::queue<std::bitset<K * 2>> queue;
  queue.push(start);

  int64_t max_distance = 0;
  int64_t sum_distance = 0;

  while (!queue.empty()) {
    const std::bitset<K* 2> from = queue.front();
    queue.pop();

    for (const char c : std::vector<char>{'A', 'C', 'G', 'T'}) {
      std::bitset<K* 2> to = from << 2;

      switch (c) {
        case 'A':
          break;
        case 'C':
          to.set(1);
          break;
        case 'G':
          to.set(0);
          break;
        case 'T':
          to.set(0);
          to.set(1);
      }

      if (kmer_set.Contains(to) && !is_visited(to) && distances[from] < limit) {
        distances[to] = distances[from] + 1;
        queue.push(to);
        max_distance = std::max(max_distance, distances[to]);
        sum_distance += distances[to];
      }
    }
  }

  return {
      distances,
      distances.size(),
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
  std::vector<std::bitset<K * 2>> path;
  int visited_nodes;
};

// A* search from "start" to "goal".
// Variable names, etc. follow https://ja.wikipedia.org/wiki/A*.
template <int K, int B, int L>
std::optional<AStarSearchResult<K>> AStarSearch(
    const KmerSet<K, B>& kmer_set, const std::bitset<K * 2>& start,
    const std::bitset<K * 2>& goal) {
  // Returns the set of l-mers present in "kmer".
  const auto get_lmers = [](const std::bitset<K * 2>& kmer) {
    // Returns the ith bit of "kmer" from left.
    const auto get = [&](int i) { return kmer[K * 2 - 1 - i]; };

    absl::flat_hash_set<std::bitset<L * 2>> lmers;

    for (int i = 0; i < K - L + 1; i++) {
      std::bitset<L * 2> lmer;

      // Sets the ith bit of "lmer" from left.
      const auto set = [&](int i, bool v) { lmer[L * 2 - 1 - i] = v; };

      for (int j = 0; j < L * 2; j++) {
        set(j, get(i * 2 + j));
      }

      lmers.insert(lmer);
    }

    return lmers;
  };

  // The set of l-mers present in "goal".
  // It is used to calculate h(n).
  const auto lmers_goal = get_lmers(goal);

  // h(n) gives an approximate of distance(n, goal).
  const auto h = [&](const std::bitset<K * 2>& n) {
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
  absl::flat_hash_map<std::bitset<K * 2>, int> f;

  // g(n) gives an approximate of distance(start, n).
  const auto g = [&](const std::bitset<K * 2>& n) { return f[n] - h(n); };

  const auto cmp = [&](const std::bitset<K * 2>& lhs,
                       const std::bitset<K * 2>& rhs) {
    return f[lhs] > f[rhs];
  };

  std::priority_queue<std::bitset<K * 2>, std::vector<std::bitset<K * 2>>,
                      decltype(cmp)>
      open(cmp);

  absl::flat_hash_set<std::bitset<K * 2>> close;

  // This will be used to re-construct the path from "start" to "goal".
  absl::flat_hash_map<std::bitset<K * 2>, std::bitset<K * 2>> parents;

  f[start] = h(start);
  open.emplace(start);

  while (!open.empty()) {
    spdlog::debug("open.size() = {}", open.size());
    spdlog::debug("close.size() = {}", close.size());

    const std::bitset<K* 2> n = open.top();
    open.pop();

    spdlog::debug("n = {}", n.to_string());

    if (goal == n) {
      // Re-constrcut the path.
      std::vector<std::bitset<K * 2>> path;
      path.push_back(n);
      while (parents.find(path.back()) != parents.end()) {
        path.push_back(parents[path.back()]);
      }
      std::reverse(path.begin(), path.end());
      int visited_nodes = f.size();
      spdlog::info("path.size() = {}", path.size());
      return {{path, visited_nodes}};
    }

    for (const char c : std::vector<char>{'A', 'C', 'G', 'T'}) {
      std::bitset<K* 2> to = n << 2;

      switch (c) {
        case 'A':
          break;
        case 'C':
          to.set(1);
          break;
        case 'G':
          to.set(0);
          break;
        case 'T':
          to.set(0);
          to.set(1);
      }

      if (kmer_set.Contains(to)) {
        const int new_f = g(n) + 1 + h(to);

        spdlog::debug("to = {}", to.to_string());
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

int main() {
  // spdlog::set_level(spdlog::level::debug);
  std::srand(std::time(nullptr));
  std::ios_base::sync_with_stdio(false);

  const int K = 31;
  const int B = 6;

  spdlog::info("constructing kmer_set");
  KmerSet<K, B> kmer_set(std::cin, 100'000'000, 500'000'000);
  spdlog::info("constructed kmer_set");
  spdlog::info("kmer_set.Size() = {}", kmer_set.Size());

  const auto start = kmer_set.Sample();

  spdlog::info("searching: start = {}", start.to_string());

  const BFSResult bfs_result = BFS(kmer_set, start, 1000);

  spdlog::info("searched: start = {}", start.to_string());
  spdlog::info("reached_nodes = {}", bfs_result.reached_nodes);
  spdlog::info("max_distance = {}", bfs_result.max_distance);
  spdlog::info("avg_distance = {:.1f}", bfs_result.avg_distance);

  std::map<int, int> distance_count;

  for (const auto& p : bfs_result.distances) {
    distance_count[p.second] += 1;
    spdlog::debug("distances[{}] = {}", p.first.to_string(), p.second);
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

  const std::bitset<K* 2> goal = [&] {
    std::vector<std::bitset<K * 2>> nodes;

    for (const auto& p : bfs_result.distances) {
      nodes.push_back(p.first);
    }

    return nodes[rand() % nodes.size()];
  }();

  spdlog::info("goal = {}", goal.to_string());
  spdlog::info("distances[{}] = {}", goal.to_string(),
               bfs_result.distances.find(goal)->second);
  spdlog::info("distance_count_sum[{}] = {}",
               bfs_result.distances.find(goal)->second,
               distance_count_sum[bfs_result.distances.find(goal)->second]);
  spdlog::info("starting A* search from {} to {}", start.to_string(),
               goal.to_string());

  const auto a_star_search_result = AStarSearch<K, B, 5>(kmer_set, start, goal);

  if (!a_star_search_result) {
    spdlog::error("A* search failed");
  } else {
    spdlog::info("finished A* search from {} to {}", start.to_string(),
                 goal.to_string());
    spdlog::info("visited_nodes = {}", (*a_star_search_result).visited_nodes);
    spdlog::info("distance = {}", (*a_star_search_result).path.size() - 1);
  }

  {
    spdlog::info("dumping kmer_set");
    std::ofstream kmer_dump_file;
    kmer_dump_file.open("kmer_dump.txt");
    kmer_set.Dump(kmer_dump_file);
    kmer_dump_file.close();
    spdlog::info("dumped kmer_set");
  }

  {
    spdlog::info("loading kmer_set_copy");
    std::ifstream kmer_dump_file{"kmer_dump.txt"};
    KmerSet<K, B> kmer_set_copy(kmer_dump_file);
    spdlog::info("loaded kmer_set_copy");
    spdlog::info("kmer_set_copy.Size() = {}", kmer_set_copy.Size());
  }
}

#include <algorithm>
#include <array>
#include <atomic>
#include <bitset>
#include <cstdlib>
#include <fstream>
#include <iostream>
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

  std::bitset<K * 2> Sample() {
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

struct SearchResult {
  int64_t reached_nodes;
  int64_t max_distance;
  double avg_distance;
};

// BFS from "start". "limit" is used to specify the search space.
template <int K, int B>
SearchResult Search(const KmerSet<K, B>& kmer_set,
                    const std::bitset<K * 2>& start, const int limit) {
  if (!kmer_set.Contains(start)) {
    spdlog::error("node not found: {}", start.to_string());
    std::exit(1);
  }

  absl::flat_hash_map<std::bitset<K * 2>, int64_t> distances;
  distances[start] = 0;

  const auto is_visited = [&](const std::bitset<K * 2>& node) {
    auto it = distances.find(node);
    return it != distances.end();
  };

  std::queue<std::bitset<K * 2>> queue;
  queue.push(start);

  int64_t max_distance = 0;
  int64_t sum_distance = 0;

  while (!queue.empty()) {
    const auto& from = queue.front();
    queue.pop();

    for (const char c : "ACGT") {
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
      distances.size(),
      max_distance,
      (double)sum_distance / distances.size(),
  };

  spdlog::info("reached {} nodes", distances.size());
  spdlog::info("max_distance = {}", max_distance);

  double avg_distance = (double)sum_distance / distances.size();
  spdlog::info("avg_distance = {:.1f}", avg_distance);
}

int main() {
  std::ios_base::sync_with_stdio(false);

  spdlog::info("constructing kmer_set");
  KmerSet<21, 6> kmer_set(std::cin, 100'000'000, 500'000'000);
  spdlog::info("constructed kmer_set");
  spdlog::info("kmer_set.Size() = {}", kmer_set.Size());

  std::mutex mu;

  std::vector<std::thread> threads;

  for (int i = 0; i < std::thread::hardware_concurrency(); i++) {
    threads.emplace_back([&] {
      const auto& start = kmer_set.Sample();

      {
        std::lock_guard l{mu};
        spdlog::info("searching: start = {}", start.to_string());
      }

      auto result = Search(kmer_set, start, 1000);

      {
        std::lock_guard l{mu};
        spdlog::info("searched: start = {}", start.to_string());
        spdlog::info("reached_nodes = {}", result.reached_nodes);
        spdlog::info("max_distance = {}", result.max_distance);
        spdlog::info("avg_distance = {:.1f}", result.avg_distance);
      }
    });
  }

  for (auto& thread : threads) {
    thread.join();
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
    KmerSet<21, 6> kmer_set_copy(kmer_dump_file);
    spdlog::info("loaded kmer_set_copy");
    spdlog::info("kmer_set_copy.Size() = {}", kmer_set_copy.Size());
  }
}

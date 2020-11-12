#ifndef KMER_SET_COMPACT_H_
#define KMER_SET_COMPACT_H_

#include <algorithm>
#include <atomic>
#include <cassert>
#include <mutex>
#include <thread>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "kmer_set.h"
#include "range.h"
#include "suffix_array.h"

// "ACCG", "CCGT", "CGTT" -> "ACCGTT"
template <int K>
std::string GetUnitigFromKmers(const std::vector<Kmer<K>>& kmers) {
  std::string s;
  s.reserve(K + kmers.size() - 1);
  s += kmers[0].String();
  for (int64_t i = 1; i < (int64_t)kmers.size(); i++) {
    s += kmers[i].String()[K - 1];
  }
  return s;
}

template <int K, typename KeyType>
std::vector<std::string> GetUnitigsCanonical(
    const KmerSet<K, KeyType>& kmer_set, int n_workers) {
  // Consider a bi-directional graph.

  // Returns neighboring k-mers. Complements are considered.
  const auto get_neighbors = [&](const Kmer<K>& kmer) {
    absl::flat_hash_set<Kmer<K>> neighbors;

    for (const Kmer<K>& next : kmer.Nexts()) {
      const Kmer<K> next_canonical = next.Canonical();
      if (kmer != next_canonical && kmer_set.Contains(next_canonical))
        neighbors.insert(next_canonical);
    }

    for (const Kmer<K>& prev : kmer.Prevs()) {
      const Kmer<K> prev_canonical = prev.Canonical();
      if (kmer != prev_canonical && kmer_set.Contains(prev_canonical))
        neighbors.insert(prev_canonical);
    }

    return neighbors;
  };

  // Terminal nodes of simple paths.
  // Every node in a simple path has degree of 0, 1, or 2.
  const std::vector<Kmer<K>> terminal_kmers = kmer_set.Find(
      [&](const Kmer<K>& kmer) {
        const auto neighbors = get_neighbors(kmer);

        if (neighbors.size() >= 3) return false;
        if (neighbors.size() <= 1) return true;

        for (const Kmer<K>& neighbor : neighbors)
          if (get_neighbors(neighbor).size() >= 3) return true;

        return false;
      },
      n_workers);

  // If "circular" is false, finds a simple path from "start".
  // "start" should be one of "terminal_kmers" in this case.
  // If "circular" is true, finds a path in a non-branching loop.
  // "start" should belong to a non-branching loop in this case.
  const auto find_path = [&](const Kmer<K>& start, bool circular) {
    std::vector<Kmer<K>> path;

    Kmer<K> current = start;

    while (true) {
      path.push_back(current);

      std::vector<Kmer<K>> unvisited_neighbors;

      for (const Kmer<K>& neighbor : get_neighbors(current)) {
        if (get_neighbors(neighbor).size() <= 2 &&
            std::find(path.begin(), path.end(), neighbor) == path.end())
          unvisited_neighbors.push_back(neighbor);
      }

      if (!circular) {
        assert(unvisited_neighbors.size() <= 1);
      } else {
        circular = false;
      }

      if (unvisited_neighbors.size() == 0) break;

      const Kmer<K> unvisited_neighbor = unvisited_neighbors[0];

      current = unvisited_neighbor;
    }

    return path;
  };

  // Finds unitigs in a path.
  // As k-mers in a path can either be used as it is or as its complement,
  // a path can be split into multiple segments, leading to multiple unitigs.
  const auto get_unitigs = [](const std::vector<Kmer<K>>& path) {
    // Returns true if the (K-1)-suffix of "lhs" is equal to the (K-1)-prefix of
    // "rhs".
    const auto is_joinable = [](const Kmer<K>& lhs, const Kmer<K>& rhs) {
      for (int i = 0; i < K - 1; i++) {
        if (lhs.Get(i + 1) != rhs.Get(i)) return false;
      }

      return true;
    };

    std::vector<std::string> unitigs;
    std::vector<std::vector<Kmer<K>>> joinable_paths;

    for (const Kmer<K>& kmer : path) {
      if (joinable_paths.size() &&
          is_joinable(joinable_paths.back().back(), kmer)) {
        joinable_paths.back().push_back(kmer);
      } else if (joinable_paths.size() &&
                 is_joinable(joinable_paths.back().back(), kmer.Complement())) {
        joinable_paths.back().push_back(kmer.Complement());
      } else {
        joinable_paths.push_back(std::vector<Kmer<K>>{kmer});
      }
    }

    for (const std::vector<Kmer<K>>& path : joinable_paths) {
      unitigs.push_back(GetUnitigFromKmers(path));
    }

    return unitigs;
  };

  std::vector<std::string> unitigs;
  absl::flat_hash_set<Kmer<K>> visited;

  std::vector<std::thread> threads;
  std::mutex mu_unitigs;
  std::mutex mu_visited;

  for (const Range& range : Range(0, terminal_kmers.size()).Split(n_workers)) {
    threads.emplace_back([&, range] {
      std::vector<std::string> buf_unitigs;
      absl::flat_hash_set<Kmer<K>> buf_visited;

      range.ForEach([&](int i) {
        const Kmer<K> terminal_kmer = terminal_kmers[i];

        std::vector<Kmer<K>> path = find_path(terminal_kmer, false);

        // Not to process the same path for multiple times.
        if (path[0] < path[path.size() - 1]) return;

        for (const Kmer<K>& kmer : path) buf_visited.insert(kmer);

        for (const std::string& unitig : get_unitigs(path)) {
          buf_unitigs.push_back(unitig);
        }
      });

      std::vector<std::thread> threads;

      threads.emplace_back([&] {
        std::lock_guard _{mu_unitigs};
        unitigs.reserve(unitigs.size() + buf_unitigs.size());
        for (const std::string& unitig : buf_unitigs) unitigs.push_back(unitig);
      });

      threads.emplace_back([&] {
        std::lock_guard _{mu_visited};
        visited.reserve(visited.size() + buf_visited.size());
        for (const Kmer<K>& kmer : buf_visited) visited.insert(kmer);
      });

      for (std::thread& thread : threads) thread.join();
    });
  }

  for (std::thread& thread : threads) thread.join();

  // Consider non-branching loops.
  // This is hard to parallelize.

  const std::vector<Kmer<K>> nodes_in_non_branching_loop = kmer_set.Find(
      [&](const Kmer<K>& kmer) {
        return get_neighbors(kmer).size() == 2 &&
               visited.find(kmer) == visited.end();
      },
      n_workers);

  for (const Kmer<K>& kmer : nodes_in_non_branching_loop) {
    if (visited.find(kmer) != visited.end()) continue;

    std::vector<Kmer<K>> path = find_path(kmer, true);

    for (const Kmer<K>& kmer : path) {
      visited.insert(kmer);
    }

    for (const std::string& unitig : get_unitigs(path)) {
      unitigs.push_back(unitig);
    }
  }

  const std::vector<Kmer<K>> branching_kmers = kmer_set.Find(
      [&](const Kmer<K>& kmer) { return get_neighbors(kmer).size() >= 3; },
      n_workers);

  for (const Kmer<K>& kmer : branching_kmers) {
    unitigs.push_back(kmer.String());
  }

  return unitigs;
}

template <int K, typename KeyType>
std::vector<std::string> GetUnitigs(const KmerSet<K, KeyType>& kmer_set,
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

  const auto start_kmers = kmer_set.Find(
      [&](const Kmer<K>& kmer) {
        const auto prevs = get_prevs(kmer);

        // If the k-mer has no incoming edges.
        if (prevs.size() == 0) return true;

        // If the k-mer has multiple incoming edges.
        if (prevs.size() >= 2) return true;

        // There is an edge from "prev" to "kmer".
        const Kmer<K> prev = prevs[0];

        const auto prev_nexts = get_nexts(prev);

        // If "prev" has multiple outgoing edges.
        if (prev_nexts.size() >= 2) return true;

        return false;
      },
      n_workers);

  const auto end_kmers = [&] {
    const auto v = kmer_set.Find(
        [&](const Kmer<K>& kmer) {
          const auto nexts = get_nexts(kmer);

          // If the k-mer has no outgoing edges.
          if (nexts.size() == 0) return true;

          // If the k-mer has multiple outgoing edges.
          if (nexts.size() >= 2) return true;

          // There is an edge from "kmer" to "next".
          const Kmer<K> next = nexts[0];

          const auto next_prevs = get_prevs(next);

          // If "next" has multiple incoming edges.
          if (next_prevs.size() >= 2) return true;

          return false;
        },
        n_workers);

    const absl::flat_hash_set<Kmer<K>> s(v.begin(), v.end());

    return s;
  }();

  std::vector<std::string> unitigs;
  KmerSet<K, KeyType> visited;

  std::mutex mu_unitigs;
  std::mutex mu_visited;

  std::vector<std::thread> threads;

  for (const Range& range : Range(0, start_kmers.size()).Split(n_workers)) {
    threads.emplace_back([&, range] {
      std::vector<std::string> buf_unitigs;
      KmerSet<K, KeyType> buf_visited;

      range.ForEach([&](int64_t i) {
        const Kmer<K>& start_kmer = start_kmers[i];
        std::vector<Kmer<K>> path;

        Kmer<K> current = start_kmer;
        while (true) {
          buf_visited.Add(current);
          path.push_back(current);
          if (end_kmers.find(current) != end_kmers.end()) break;
          current = get_nexts(current)[0];
        }

        buf_unitigs.push_back(GetUnitigFromKmers(path));
      });

      std::vector<std::thread> threads;

      threads.emplace_back([&] {
        std::lock_guard _{mu_visited};
        visited.Add(buf_visited, n_workers);
      });

      threads.emplace_back([&] {
        std::lock_guard _{mu_unitigs};
        unitigs.reserve(unitigs.size() + buf_unitigs.size());
        for (const auto& unitig : buf_unitigs) unitigs.push_back(unitig);
      });

      for (std::thread& thread : threads) thread.join();
    });
  }

  for (std::thread& thread : threads) thread.join();

  // Consider loops where every no has 1 incoming edge and 1 outgoing edge.
  {
    const std::vector<Kmer<K>> not_visited = kmer_set.Find(
        [&](const Kmer<K>& kmer) { return !visited.Contains(kmer); },
        n_workers);

    for (const Kmer<K>& kmer : not_visited) {
      if (visited.Contains(kmer)) continue;

      Kmer<K> current = kmer;

      std::vector<Kmer<K>> path;

      while (!visited.Contains(current)) {
        path.push_back(current);
        visited.Add(current);
        current = get_nexts(current)[0];
      }

      unitigs.push_back(GetUnitigFromKmers(path));
    }
  }

  return unitigs;
}

// KmerSetSetCompact holds a set of k-mers in less space.
template <int K>
class KmerSetCompact {
 public:
  template <typename KeyType>
  KmerSetCompact(const KmerSet<K, KeyType>& kmer_set, bool canonical,
                 int n_workers)
      : canonical_(canonical) {
    const std::vector<std::string> unitigs =
        canonical ? GetUnitigsCanonical(kmer_set, n_workers)
                  : GetUnitigs(kmer_set, n_workers);

    std::string concatenated;

    // Pre-calculate the required length.
    {
      int64_t sum = 0;
      for (const std::string& unitig : unitigs) sum += unitig.length();
      concatenated.reserve(sum);
    }

    for (const std::string& unitig : unitigs) {
      boundary_.push_back(concatenated.length());
      concatenated += unitig;
    }

    sa_ = SuffixArray(concatenated, n_workers);
  }

  bool Contains(const Kmer<K>& kmer) const {
    const auto contains = [&](const Kmer<K>& kmer) {
      std::vector<int64_t> v = sa_.Find(kmer.String());

      // Returns the largest j such that boundary_[j] <= i;
      const auto f = [&](int64_t i) {
        // begin == -1 or boundary_[begin] <= i
        int64_t begin = -1;

        // end == boundary_.size() or boundary_[end] > i
        int64_t end = boundary_.size();

        while (end - begin > 1) {
          int64_t mid = (begin + end) / 2;

          if (boundary_[mid] <= i)
            begin = mid;
          else
            end = mid;
        }

        return begin;
      };

      return std::any_of(v.begin(), v.end(),
                         [&](int64_t i) { return f(i) == f(i + K - 1); });
    };

    if (canonical_) {
      return contains(kmer) && contains(kmer.Complement());
    } else {
      return contains(kmer);
    }
  }

  // Returns the k-mers that match the condition.
  template <typename Pred>
  std::vector<Kmer<K>> Find(Pred&& pred, int n_workers) const {
    std::vector<Kmer<K>> kmers;
    std::mutex mu;

    std::vector<std::thread> threads;

    for (const Range& range : Range(0, boundary_.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<Kmer<K>> buf;

        range.ForEach([&](int64_t i) {
          const int64_t begin = boundary_[i];
          const int64_t end = (i == (int64_t)boundary_.size() - 1)
                                  ? sa_.String().length()
                                  : boundary_[i + 1];
          for (int64_t j = 0; begin + j + K - 1 < end; j++) {
            const Kmer<K> kmer{sa_.String().substr(begin + j, K)};
            if (pred(kmer)) buf.push_back(canonical_ ? kmer.Canonical() : kmer);
          }
        });

        std::lock_guard _{mu};
        kmers.reserve(kmers.size() + buf.size());
        for (const Kmer<K>& kmer : buf) kmers.push_back(kmer);
      });
    }

    for (std::thread& thread : threads) thread.join();

    return kmers;
  }

  // Return all the k-mers.
  std::vector<Kmer<K>> Find(int n_workers) const {
    return Find([&](const Kmer<K>&) { return true; }, n_workers);
  }

  int64_t Size() const { return sa_.String().length(); }

 private:
  bool canonical_;

  SuffixArray sa_;

  std::vector<int64_t> boundary_;
};

#endif
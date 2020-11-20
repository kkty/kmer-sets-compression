#ifndef CORE_UNITIGS_H_
#define CORE_UNITIGS_H_

#include <cassert>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "core/kmer_set.h"
#include "core/range.h"

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

// Constructs unitigs from a kmer set where canonical kmers are stored.
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

        // Here, neighbors.size() is 2.

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

    // Splits the path into multiple paths each of which can be made into a
    // unitig by using "GetUnitigFromKmers()".
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

    std::vector<std::string> unitigs;

    for (const std::vector<Kmer<K>>& path : joinable_paths) {
      unitigs.push_back(GetUnitigFromKmers(path));
    }

    return unitigs;
  };

  std::vector<std::string> unitigs;
  absl::flat_hash_set<Kmer<K>> visited;

  // For each kmer in terminal_kmers, finds the unitig from it.
  {
    std::vector<std::thread> threads;
    std::mutex mu_unitigs;
    std::mutex mu_visited;

    for (const Range& range :
         Range(0, terminal_kmers.size()).Split(n_workers)) {
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

        // 3 * n_workers threads are created in total, but no more than
        // n_workers threads are active at any point in time.

        std::vector<std::thread> threads;

        threads.emplace_back([&] {
          std::lock_guard lck(mu_unitigs);
          for (std::string& unitig : buf_unitigs)
            unitigs.push_back(std::move(unitig));
        });

        threads.emplace_back([&] {
          std::lock_guard lck(mu_visited);
          for (const Kmer<K>& kmer : buf_visited) visited.insert(kmer);
        });

        for (std::thread& t : threads) t.join();
      });
    }

    for (std::thread& t : threads) t.join();
  }

  // Considers non-branching loops.
  // This is hard to parallelize.
  {
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
  }

  // Kmers with branches are themselves treated as unitigs.
  {
    const std::vector<Kmer<K>> branching_kmers = kmer_set.Find(
        [&](const Kmer<K>& kmer) { return get_neighbors(kmer).size() >= 3; },
        n_workers);

    for (const Kmer<K>& kmer : branching_kmers) {
      unitigs.push_back(kmer.String());
    }
  }

  return unitigs;
}

// Constructs unitigs from a kmer set.
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

  // Kmers where a unitig starts.
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

  // Kmers where a unitig ends.
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

  // For each kmer in start_kmers, finds the unitig from it.
  {
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

        // 3 * n_workers threads are created in total, but no more than
        // n_workers threads are active at any point in time.

        std::vector<std::thread> threads;

        threads.emplace_back([&] {
          std::lock_guard lck(mu_visited);
          visited.Add(buf_visited, 1);
        });

        threads.emplace_back([&] {
          std::lock_guard lck(mu_unitigs);

          for (std::string& unitig : buf_unitigs)
            unitigs.push_back(std::move(unitig));
        });

        for (std::thread& t : threads) t.join();
      });
    }

    for (std::thread& t : threads) t.join();
  }

  // Considers loops where every node has 1 incoming edge and 1 outgoing edge.
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

#endif
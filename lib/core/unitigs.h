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
    assert(kmers[i - 1].String().substr(1, K - 1) ==
           kmers[i].String().substr(0, K - 1));

    s += kmers[i].String()[K - 1];
  }

  return s;
}

// Constructs unitigs from a set of canonical kmers.
template <int K, typename KeyType>
std::vector<std::string> GetUnitigsCanonical(
    const KmerSet<K, KeyType>& kmer_set, int n_workers) {
  // If the second element is true, the edge connects the same side of two kmers.
  using Neighbor = std::pair<Kmer<K>, bool>;

  // Lists neighbors of the right side of "kmer".
  const auto GetNeighborsRight = [&](const Kmer<K>& kmer) {
    std::vector<Neighbor> neighbors;

    for (const Kmer<K>& next : kmer.Nexts()) {
      if (kmer != next && kmer_set.Contains(next)) {
        neighbors.emplace_back(next, false);
      }

      const Kmer<K> next_complement = next.Complement();

      if (kmer != next_complement && kmer_set.Contains(next_complement)) {
        neighbors.emplace_back(next_complement, true);
      }
    }

    return neighbors;
  };

  // Lists neighbors of the left side of "kmer".
  const auto GetNeighborsLeft = [&](const Kmer<K>& kmer) {
    std::vector<Neighbor> neighbors;

    for (const Kmer<K>& prev : kmer.Prevs()) {
      if (kmer != prev && kmer_set.Contains(prev)) {
        neighbors.emplace_back(prev, false);
      }

      const Kmer<K> prev_complement = prev.Complement();

      if (kmer != prev_complement && kmer_set.Contains(prev_complement)) {
        neighbors.emplace_back(prev_complement, true);
      }
    }

    return neighbors;
  };

  // Returns true if the kmer has no mates on the left side.
  const auto IsTerminalLeft = [&](const Kmer<K>& kmer) {
    const std::vector<Neighbor> neighbors = GetNeighborsLeft(kmer);

    if (neighbors.size() != 1) return true;

    Kmer<K> neighbor;
    bool is_same_side;

    std::tie(neighbor, is_same_side) = neighbors.front();

    if (is_same_side) {
      if (GetNeighborsLeft(neighbor).size() != 1) return true;
    } else {
      if (GetNeighborsRight(neighbor).size() != 1) return true;
    }

    return false;
  };

  // Returns true if the kmer has no mates on the right side.
  const auto IsTerminalRight = [&](const Kmer<K>& kmer) {
    const std::vector<Neighbor> neighbors = GetNeighborsRight(kmer);

    if (neighbors.size() != 1) return true;

    Kmer<K> neighbor;
    bool is_same_side;

    std::tie(neighbor, is_same_side) = neighbors.front();

    if (is_same_side) {
      if (GetNeighborsRight(neighbor).size() != 1) return true;
    } else {
      if (GetNeighborsLeft(neighbor).size() != 1) return true;
    }

    return false;
  };

  // kmers that has no mates on the left side, but has a mate on the right side.
  const std::vector<Kmer<K>> terminals_left = kmer_set.Find(
      [&](const Kmer<K>& kmer) {
        return IsTerminalLeft(kmer) && !IsTerminalRight(kmer);
      },
      n_workers);

  // kmers that has no mates on the right side, but has a mate on the left side.
  const std::vector<Kmer<K>> terminals_right = kmer_set.Find(
      [&](const Kmer<K>& kmer) {
        return !IsTerminalLeft(kmer) && IsTerminalRight(kmer);
      },
      n_workers);

  // kmers that has no mates on the both sides.
  const std::vector<Kmer<K>> terminals_both = kmer_set.Find(
      [&](const Kmer<K>& kmer) {
        return IsTerminalLeft(kmer) && IsTerminalRight(kmer);
      },
      n_workers);

  // If "is_right_side" is true, finds a path from the right side of "start".
  // If "is_right_side" is false, finds a path from the left side of "start".
  const auto FindPath = [&](Kmer<K> start, bool is_right_side) {
    Kmer<K> current = start;
    std::vector<Kmer<K>> path;

    while (true) {
      path.push_back(is_right_side ? current : current.Complement());

      if (is_right_side) {
        if (IsTerminalRight(current)) break;
      } else {
        if (IsTerminalLeft(current)) break;
      }

      std::vector<Neighbor> neighbors = is_right_side
                                            ? GetNeighborsRight(current)
                                            : GetNeighborsLeft(current);

      assert(neighbors.size() == 1);

      bool is_same_side;

      std::tie(current, is_same_side) = neighbors.front();

      if (is_same_side) is_right_side = !is_right_side;
    }

    return path;
  };

  std::vector<std::string> unitigs;
  absl::flat_hash_set<Kmer<K>> visited;

  const auto MoveFromBuffer = [&](std::mutex& mu_unitigs,
                                  std::mutex& mu_visited,
                                  std::vector<std::string>& buf_unitigs,
                                  std::vector<Kmer<K>>& buf_visited) {
    bool done_unitigs = false;
    bool done_visited = false;

    while (!done_unitigs || !done_visited) {
      if (!done_unitigs && mu_unitigs.try_lock()) {
        for (std::string& unitig : buf_unitigs) {
          unitigs.push_back(std::move(unitig));
        }

        mu_unitigs.unlock();
        done_unitigs = true;
      }

      if (!done_visited && mu_visited.try_lock()) {
        for (const Kmer<K>& kmer : buf_visited) {
          visited.insert(kmer);
        }

        mu_visited.unlock();
        done_visited = true;
      }
    }
  };

  // Processes kmers in terminals_both.
  {
    std::vector<std::thread> threads;
    std::mutex mu_unitigs;
    std::mutex mu_visited;

    for (const Range& range :
         Range(0, terminals_both.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<std::string> buf_unitigs;
        std::vector<Kmer<K>> buf_visited;

        range.ForEach([&](int64_t i) {
          const Kmer<K>& kmer = terminals_both[i];
          buf_unitigs.push_back(kmer.String());
          buf_visited.push_back(kmer);
        });

        MoveFromBuffer(mu_unitigs, mu_visited, buf_unitigs, buf_visited);
      });
    }

    for (std::thread& t : threads) t.join();
  }

  // Processes paths from terminals_left.
  {
    std::vector<std::thread> threads;
    std::mutex mu_unitigs;
    std::mutex mu_visited;

    for (const Range& range :
         Range(0, terminals_left.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<std::string> buf_unitigs;
        std::vector<Kmer<K>> buf_visited;

        range.ForEach([&](int64_t i) {
          std::vector<Kmer<K>> path = FindPath(terminals_left[i], true);

          if (path.front().Canonical() < path.back().Canonical()) return;

          for (const Kmer<K>& kmer : path)
            buf_visited.push_back(kmer.Canonical());

          buf_unitigs.push_back(GetUnitigFromKmers(path));
        });

        MoveFromBuffer(mu_unitigs, mu_visited, buf_unitigs, buf_visited);
      });
    }

    for (std::thread& t : threads) t.join();
  }

  // Processes paths from terminals_right.
  {
    std::vector<std::thread> threads;
    std::mutex mu_unitigs;
    std::mutex mu_visited;

    for (const Range& range :
         Range(0, terminals_right.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<std::string> buf_unitigs;
        std::vector<Kmer<K>> buf_visited;

        range.ForEach([&](int64_t i) {
          std::vector<Kmer<K>> path = FindPath(terminals_right[i], false);

          if (path.front().Canonical() < path.back().Canonical()) return;

          for (const Kmer<K>& kmer : path)
            buf_visited.push_back(kmer.Canonical());

          buf_unitigs.push_back(GetUnitigFromKmers(path));
        });

        MoveFromBuffer(mu_unitigs, mu_visited, buf_unitigs, buf_visited);
      });
    }

    for (std::thread& t : threads) t.join();
  }

  // Considers non-branching loops.

  std::vector<Kmer<K>> not_visited = kmer_set.Find(
      [&](const Kmer<K>& kmer) { return visited.find(kmer) == visited.end(); },
      n_workers);

  for (const Kmer<K>& start : not_visited) {
    if (visited.find(start) != visited.end()) continue;

    bool is_right_side = true;
    Kmer<K> current = start;
    std::vector<Kmer<K>> path;

    while (visited.find(current) == visited.end()) {
      visited.insert(current);
      path.push_back(is_right_side ? current : current.Complement());

      std::vector<Neighbor> neighbors = is_right_side
                                            ? GetNeighborsRight(current)
                                            : GetNeighborsLeft(current);
      assert(neighbors.size() == 1);
      bool is_same_side;
      std::tie(current, is_same_side) = neighbors.front();
      if (is_same_side) is_right_side = !is_right_side;
    }

    unitigs.push_back(GetUnitigFromKmers(path));
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
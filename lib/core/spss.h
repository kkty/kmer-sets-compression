#ifndef CORE_SPSS_H_
#define CORE_SPSS_H_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "boost/sort/sort.hpp"
#include "core/kmer_set.h"
#include "core/range.h"
#include "parallel_disjoint_set.h"
#include "spdlog/spdlog.h"

namespace internal {

// "ACCG", "CCGT", "CGTT" -> "ACCGTT"
template <int K>
std::string ConcatenateKmers(const std::vector<Kmer<K>>& kmers) {
  std::string s;
  s.reserve(K + kmers.size() - 1);

  s += kmers[0].String();

  for (std::int64_t i = 1; i < static_cast<std::int64_t>(kmers.size()); i++) {
    assert(kmers[i - 1].String().substr(1, K - 1) ==
           kmers[i].String().substr(0, K - 1));

    s += kmers[i].String()[K - 1];
  }

  return s;
}

// Returns the complement of s. The complement of s can be constructed by
// reversing s and replacing A, C, G and T with T, G, C and A respectively.
std::string Complement(std::string s) {
  std::reverse(s.begin(), s.end());

  for (char& c : s) {
    switch (c) {
      case 'A':
        c = 'T';
        break;
      case 'C':
        c = 'G';
        break;
      case 'G':
        c = 'C';
        break;
      case 'T':
        c = 'A';
        break;
      default:
        assert(false);
    }
  }

  return s;
}

}  // namespace internal

// Constructs unitigs from a kmer set.
template <int K, int N, typename KeyType>
std::vector<std::string> GetUnitigs(const KmerSet<K, N, KeyType>& kmer_set,
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

  // Kmers where a unitig starts.
  const auto start_kmers = kmer_set.Find(
      [&](const Kmer<K>& kmer) {
        const auto prevs = GetPrevs(kmer);

        // If the k-mer has no incoming edges.
        if (prevs.size() == 0) return true;

        // If the k-mer has multiple incoming edges.
        if (prevs.size() >= 2) return true;

        // There is an edge from "prev" to "kmer".
        const Kmer<K> prev = prevs[0];

        const auto prev_nexts = GetNexts(prev);

        // If "prev" has multiple outgoing edges.
        if (prev_nexts.size() >= 2) return true;

        return false;
      },
      n_workers);

  // Kmers where a unitig ends.
  const auto end_kmers = [&] {
    const auto v = kmer_set.Find(
        [&](const Kmer<K>& kmer) {
          const auto nexts = GetNexts(kmer);

          // If the k-mer has no outgoing edges.
          if (nexts.size() == 0) return true;

          // If the k-mer has multiple outgoing edges.
          if (nexts.size() >= 2) return true;

          // There is an edge from "kmer" to "next".
          const Kmer<K> next = nexts[0];

          const auto next_prevs = GetPrevs(next);

          // If "next" has multiple incoming edges.
          if (next_prevs.size() >= 2) return true;

          return false;
        },
        n_workers);

    const absl::flat_hash_set<Kmer<K>> s(v.begin(), v.end());

    return s;
  }();

  std::vector<std::string> unitigs;
  KmerSet<K, N, KeyType> visited;

  // For each kmer in start_kmers, finds the unitig from it.
  {
    std::mutex mu_unitigs;
    std::mutex mu_visited;

    std::vector<std::thread> threads;

    for (const Range& range : Range(0, start_kmers.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<std::string> buf_unitigs;
        KmerSet<K, N, KeyType> buf_visited;

        for (std::int64_t i : range) {
          const Kmer<K>& start_kmer = start_kmers[i];
          std::vector<Kmer<K>> path;

          Kmer<K> current = start_kmer;
          while (true) {
            buf_visited.Add(current);
            path.push_back(current);
            if (end_kmers.find(current) != end_kmers.end()) break;
            current = GetNexts(current)[0];
          }

          buf_unitigs.push_back(internal::ConcatenateKmers(path));
        }

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
        current = GetNexts(current)[0];
      }

      unitigs.push_back(internal::ConcatenateKmers(path));
    }
  }

  return unitigs;
}

// Constructs unitigs from a set of canonical kmers.
template <int K, int N, typename KeyType>
std::vector<std::string> GetUnitigsCanonical(
    const KmerSet<K, N, KeyType>& kmer_set, int n_workers) {
  // If the second element is true, the edge connects the same side of two
  // kmers.
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

  spdlog::debug("constructing terminals_left");

  // kmers that has no mates on the left side, but has a mate on the right side.
  std::vector<Kmer<K>> terminals_left = kmer_set.Find(
      [&](const Kmer<K>& kmer) { return IsTerminalLeft(kmer); }, n_workers);

  spdlog::debug("constructed terminals_left");

  spdlog::debug("sorting terminals_left");

  if (n_workers == 1) {
    std::sort(terminals_left.begin(), terminals_left.end());
  } else {
    boost::sort::block_indirect_sort(terminals_left.begin(),
                                     terminals_left.end(), n_workers);
  }

  spdlog::debug("sorted terminals_left");

  spdlog::debug("constructing terminals_right");

  // kmers that has no mates on the right side, but has a mate on the left side.
  std::vector<Kmer<K>> terminals_right = kmer_set.Find(
      [&](const Kmer<K>& kmer) { return IsTerminalRight(kmer); }, n_workers);

  spdlog::debug("constructed terminals_right");

  spdlog::debug("sorting terminals_right");

  if (n_workers == 1) {
    std::sort(terminals_right.begin(), terminals_right.end());
  } else {
    boost::sort::block_indirect_sort(terminals_right.begin(),
                                     terminals_right.end(), n_workers);
  }

  spdlog::debug("sorted terminals_right");

  spdlog::debug("constructing terminals_both");

  // kmers that has no mates on the both sides.
  std::vector<Kmer<K>> terminals_both;

  std::set_intersection(terminals_left.begin(), terminals_left.end(),
                        terminals_right.begin(), terminals_right.end(),
                        std::back_inserter(terminals_both));

  spdlog::debug("constructed terminals_both");

  {
    spdlog::debug("updating terminals_left");

    std::vector<Kmer<K>> buf;
    buf.reserve(terminals_left.size() - terminals_both.size());

    std::set_difference(terminals_left.begin(), terminals_left.end(),
                        terminals_both.begin(), terminals_both.end(),
                        std::back_inserter(buf));

    buf.swap(terminals_left);

    spdlog::debug("updated terminals_left");
  }

  {
    spdlog::debug("updating terminals_right");

    std::vector<Kmer<K>> buf;
    buf.reserve(terminals_right.size() - terminals_both.size());

    std::set_difference(terminals_right.begin(), terminals_right.end(),
                        terminals_both.begin(), terminals_both.end(),
                        std::back_inserter(buf));

    buf.swap(terminals_right);

    spdlog::debug("updated terminals_right");
  }

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

  unitigs.reserve(terminals_both.size() +
                  (terminals_left.size() + terminals_right.size()) / 2);

  visited.reserve(kmer_set.Size());

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
        visited.insert(buf_visited.begin(), buf_visited.end());

        mu_visited.unlock();
        done_visited = true;
      }
    }
  };

  spdlog::debug("processing kmers in terminals_both");

  {
    std::vector<std::thread> threads;
    std::mutex mu_unitigs;
    std::mutex mu_visited;

    for (const Range& range :
         Range(0, terminals_both.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<std::string> buf_unitigs;
        std::vector<Kmer<K>> buf_visited;

        for (std::int64_t i : range) {
          const Kmer<K>& kmer = terminals_both[i];

          if (n_workers == 1) {
            unitigs.push_back(kmer.String());
            visited.insert(kmer);
          } else {
            buf_unitigs.push_back(kmer.String());
            buf_visited.push_back(kmer);
          }
        }

        if (n_workers != 1) {
          MoveFromBuffer(mu_unitigs, mu_visited, buf_unitigs, buf_visited);
        }
      });
    }

    for (std::thread& t : threads) t.join();
  }

  spdlog::debug("processed kmers in terminals_both");

  spdlog::debug("processing paths from terminals_left");

  {
    std::vector<std::thread> threads;
    std::mutex mu_unitigs;
    std::mutex mu_visited;

    for (const Range& range :
         Range(0, terminals_left.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<std::string> buf_unitigs;
        std::vector<Kmer<K>> buf_visited;

        for (std::int64_t i : range) {
          std::vector<Kmer<K>> path = FindPath(terminals_left[i], true);

          if (path.front().Canonical() < path.back().Canonical()) continue;

          for (const Kmer<K>& kmer : path) {
            if (n_workers == 1) {
              visited.insert(kmer.Canonical());
            } else {
              buf_visited.push_back(kmer.Canonical());
            }
          }

          if (n_workers == 1) {
            unitigs.push_back(internal::ConcatenateKmers(path));
          } else {
            buf_unitigs.push_back(internal::ConcatenateKmers(path));
          }
        }

        if (n_workers != 1) {
          MoveFromBuffer(mu_unitigs, mu_visited, buf_unitigs, buf_visited);
        }
      });
    }

    for (std::thread& t : threads) t.join();
  }

  spdlog::debug("processed paths from terminals_left");

  spdlog::debug("processing paths from terminals_right");

  {
    std::vector<std::thread> threads;
    std::mutex mu_unitigs;
    std::mutex mu_visited;

    for (const Range& range :
         Range(0, terminals_right.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<std::string> buf_unitigs;
        std::vector<Kmer<K>> buf_visited;

        for (std::int64_t i : range) {
          std::vector<Kmer<K>> path = FindPath(terminals_right[i], false);

          if (path.front().Canonical() < path.back().Canonical()) continue;

          for (const Kmer<K>& kmer : path) {
            if (n_workers == 1) {
              visited.insert(kmer.Canonical());
            } else {
              buf_visited.push_back(kmer.Canonical());
            }
          }

          if (n_workers == 1) {
            unitigs.push_back(internal::ConcatenateKmers(path));
          } else {
            buf_unitigs.push_back(internal::ConcatenateKmers(path));
          }
        }

        if (n_workers != 1) {
          MoveFromBuffer(mu_unitigs, mu_visited, buf_unitigs, buf_visited);
        }
      });
    }

    for (std::thread& t : threads) t.join();
  }

  spdlog::debug("processed paths from terminals_right");

  spdlog::debug("processing non-branching loops");

  std::vector<Kmer<K>> not_visited = kmer_set.Find(
      [&](const Kmer<K>& kmer) { return visited.find(kmer) == visited.end(); },
      n_workers, kmer_set.Size() - visited.size());

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

    unitigs.push_back(internal::ConcatenateKmers(path));
  }

  spdlog::debug("processed non-branching loops");

  return unitigs;
}

// For each kmer s, finds a list of integers v such that s is the prefix of
// unitigs[v[i]] for each i.
template <int K>
absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>> GetPrefixesFromUnitigs(
    const std::vector<std::string>& unitigs, int n_workers) {
  const std::int64_t n = unitigs.size();

  absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>> prefixes;

  {
    std::vector<std::thread> threads;
    std::mutex mu;

    for (const Range& range : Range(0, n).Split(n_workers)) {
      threads.emplace_back([&, range] {
        absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>> buf;

        for (std::int64_t i : range) {
          const std::string& unitig = unitigs[i];
          const Kmer<K> prefix(unitig.substr(0, K));
          buf[prefix].push_back(i);
        }

        // Moves from buffers.

        if (n_workers == 1) {
          prefixes = std::move(buf);
        } else {
          std::lock_guard lck(mu);
          prefixes.insert(buf.begin(), buf.end());
        }
      });
    }

    for (std::thread& t : threads) t.join();
  }

  return prefixes;
}

// For each kmer s, finds a list of integers v such that s is the suffix of
// unitigs[v[i]] for each i.
template <int K>
absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>> GetSuffixesFromUnitigs(
    const std::vector<std::string>& unitigs, int n_workers) {
  const std::int64_t n = unitigs.size();

  absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>> suffixes;

  {
    std::vector<std::thread> threads;
    std::mutex mu;

    for (const Range& range : Range(0, n).Split(n_workers)) {
      threads.emplace_back([&, range] {
        absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>> buf;

        for (std::int64_t i : range) {
          const std::string& unitig = unitigs[i];
          const Kmer<K> suffix(unitig.substr(unitig.length() - K, K));
          buf[suffix].push_back(i);
        }

        // Moves from buffers.

        if (n_workers == 1) {
          suffixes = std::move(buf);
        } else {
          std::lock_guard lck(mu);
          suffixes.insert(buf.begin(), buf.end());
        }
      });
    }

    for (std::thread& t : threads) t.join();
  }

  return suffixes;
}

template <int K, int N, typename KeyType>
std::vector<std::string> GetSPSS(
    const std::vector<std::string>& unitigs,
    const absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>>& prefixes,
    int n_workers, int n_buckets = 64) {
  const std::int64_t n = unitigs.size();

  // Considers a directed graph where each node represents a unitig.

  // Returns a list of outgoing edges from i.
  const auto GetEdgesOut = [&](std::int64_t i) {
    std::vector<std::int64_t> edges;

    const std::string& unitig = unitigs[i];

    const Kmer<K> suffix(unitig.substr(unitig.length() - K, K));

    for (const Kmer<K>& suffix_next : suffix.Nexts()) {
      auto it = prefixes.find(suffix_next);

      if (it == prefixes.end()) continue;

      for (std::int64_t j : it->second) {
        if (i == j) continue;

        edges.push_back(j);
      }
    }

    return edges;
  };

  // edge_in[i] is the selected incoming edge to i.
  absl::flat_hash_map<std::int64_t, std::int64_t> edge_in;

  // edge_out[i] is the selected outgoing edge from i.
  absl::flat_hash_map<std::int64_t, std::int64_t> edge_out;

  {
    spdlog::debug("constructing edge_in and edge_out");

    std::vector<std::thread> threads;
    std::vector<std::mutex> mus(n_buckets);
    std::vector<absl::flat_hash_map<std::int64_t, std::int64_t>> buf_edge_in(
        n_buckets);
    std::vector<absl::flat_hash_map<std::int64_t, std::int64_t>> buf_edge_out(
        n_buckets);

    // Acquires locks for node i and j.
    const auto AcquireLock = [&](std::int64_t i, std::int64_t j) {
      int bucket_i = i % n_buckets;
      int bucket_j = j % n_buckets;

      if (bucket_i == bucket_j) {
        mus[bucket_i].lock();
        return;
      }

      mus[std::min(bucket_i, bucket_j)].lock();
      mus[std::max(bucket_i, bucket_j)].lock();
    };

    // Releases locks for node i and j.
    const auto ReleaseLock = [&](std::int64_t i, std::int64_t j) {
      int bucket_i = i % n_buckets;
      int bucket_j = j % n_buckets;

      if (bucket_i == bucket_j) {
        mus[bucket_i].unlock();
        return;
      }

      mus[std::max(bucket_i, bucket_j)].unlock();
      mus[std::min(bucket_i, bucket_j)].unlock();
    };

    // Returns true if there is an incoming edge to i.
    const auto HasEdgeIn = [&](std::int64_t i) {
      int bucket = i % n_buckets;
      return buf_edge_in[bucket].find(i) != buf_edge_in[bucket].end();
    };

    // Returns true if there is an outgoing edge from i.
    const auto HasEdgeOut = [&](std::int64_t i) {
      int bucket = i % n_buckets;
      return buf_edge_out[bucket].find(i) != buf_edge_out[bucket].end();
    };

    // Adds an incoming edge to i.
    const auto AddEdgeIn = [&](std::int64_t i, std::int64_t j) {
      int bucket_i = i % n_buckets;
      buf_edge_in[bucket_i][i] = j;
    };

    // Add an outgoing edge from i.
    const auto AddEdgeOut = [&](std::int64_t i, std::int64_t j) {
      int bucket_i = i % n_buckets;
      buf_edge_out[bucket_i][i] = j;
    };

    for (const Range& range : Range(0, n).Split(n_workers)) {
      threads.emplace_back([&, range] {
        for (std::int64_t i : range) {
          for (std::int64_t j : GetEdgesOut(i)) {
            AcquireLock(i, j);

            if (!HasEdgeOut(i) && !HasEdgeIn(j)) {
              AddEdgeOut(i, j);
              AddEdgeIn(j, i);
            }

            ReleaseLock(i, j);
          }
        }
      });
    }

    for (std::thread& t : threads) t.join();

    spdlog::debug("moving from buffers");

    {
      boost::asio::thread_pool pool(n_workers);

      boost::asio::post(pool, [&] {
        std::int64_t size = 0;
        for (int i = 0; i < n_buckets; i++) {
          size += buf_edge_in[i].size();
        }
        edge_in.reserve(size);
        for (int i = 0; i < n_buckets; i++) {
          edge_in.insert(buf_edge_in[i].begin(), buf_edge_in[i].end());
        }
      });

      boost::asio::post(pool, [&] {
        std::int64_t size = 0;
        for (int i = 0; i < n_buckets; i++) {
          size += buf_edge_out[i].size();
        }
        edge_out.reserve(size);
        for (int i = 0; i < n_buckets; i++) {
          edge_out.insert(buf_edge_out[i].begin(), buf_edge_out[i].end());
        }
      });

      pool.join();
    }

    spdlog::debug("constructed edge_in and edge_out");
  }

  {
    ParallelDisjointSet disjoint_set(n);

    spdlog::debug("constructing disjoint_set");

    {
      std::vector<std::thread> threads;

      for (const Range& range : Range(0, n).Split(n_workers)) {
        threads.emplace_back([&, range] {
          for (std::int64_t i : range) {
            auto it = edge_out.find(i);
            if (it == edge_out.end()) continue;
            disjoint_set.Unite(i, it->second);
          }
        });
      }

      for (std::thread& t : threads) t.join();
    }

    spdlog::debug("constructed disjoint_set");

    spdlog::debug("removing loops");

    {
      absl::flat_hash_set<std::int64_t> groups;
      absl::flat_hash_set<std::int64_t> groups_with_terminals;

      spdlog::debug("constructing groups and groups_with_terminals");

      std::vector<std::thread> threads;
      std::mutex mu_groups;
      std::mutex mu_groups_with_terminals;

      for (const Range& range : Range(0, n).Split(n_workers)) {
        threads.emplace_back([&, range] {
          absl::flat_hash_set<std::int64_t> buf_groups;
          absl::flat_hash_set<std::int64_t> buf_groups_with_terminals;

          for (std::int64_t i : range) {
            int group = disjoint_set.Find(i);

            buf_groups.insert(group);

            if (edge_out.find(i) == edge_out.end()) {
              buf_groups_with_terminals.insert(group);
            }
          }

          if (n_workers == 1) {
            groups = std::move(buf_groups);
            groups_with_terminals = std::move(buf_groups_with_terminals);
          } else {
            {
              std::lock_guard lck(mu_groups);
              groups.insert(buf_groups.begin(), buf_groups.end());
            }

            {
              std::lock_guard lck(mu_groups_with_terminals);
              groups_with_terminals.insert(buf_groups_with_terminals.begin(),
                                           buf_groups_with_terminals.end());
            }
          }
        });
      }

      for (std::thread& t : threads) t.join();

      spdlog::debug("constructed groups and groups_with_terminals");

      for (std::int64_t i : groups) {
        if (groups_with_terminals.find(i) == groups_with_terminals.end()) {
          std::int64_t j = edge_out[i];
          edge_out.erase(i);
          edge_in.erase(j);
        }
      }
    }

    spdlog::debug("removed loops");
  }

  std::vector<std::int64_t> starts;

  spdlog::debug("constructing starts");

  {
    std::vector<std::thread> threads;
    std::mutex mu;

    for (const Range& range : Range(0, n).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<std::int64_t> buf;

        for (std::int64_t i : range) {
          if (edge_in.find(i) == edge_in.end()) {
            buf.push_back(i);
          }
        }

        {
          std::lock_guard lck(mu);
          starts.insert(starts.end(), buf.begin(), buf.end());
        }
      });
    }

    for (std::thread& t : threads) t.join();
  }

  spdlog::debug("constructed starts");

  std::vector<std::string> spss;

  spdlog::debug("constructing spss");

  {
    std::vector<std::thread> threads;
    std::mutex mu;

    for (const Range& range : Range(0, starts.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<std::string> buf;

        for (std::int64_t i : range) {
          std::int64_t start = starts[i];

          std::vector<std::int64_t> path;

          std::int64_t current = start;
          while (true) {
            path.push_back(current);

            auto it = edge_out.find(current);
            if (it == edge_out.end()) break;

            current = it->second;
          }

          std::string s = unitigs[path[0]];
          for (std::size_t j = 1; j < path.size(); j++) {
            s += unitigs[path[j]].substr(K - 1,
                                         unitigs[path[j]].length() - (K - 1));
          }

          buf.push_back(std::move(s));
        }

        {
          std::lock_guard lck(mu);
          spss.insert(spss.end(), buf.begin(), buf.end());
        }
      });
    }

    for (std::thread& t : threads) t.join();
  }

  spdlog::debug("constructed spss");

  return spss;
}

// Constructs a small-weight SPSS from a set of kmers.
template <int K, int N, typename KeyType>
std::vector<std::string> GetSPSS(const KmerSet<K, N, KeyType>& kmer_set,
                                 int n_workers, int n_buckets = 64) {
  spdlog::debug("constructing unitigs");

  const std::vector<std::string> unitigs = GetUnitigs(kmer_set, n_workers);

  spdlog::debug("constructed unitigs");

  spdlog::debug("constructing prefixes");

  absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>> prefixes =
      GetPrefixesFromUnitigs<K>(unitigs, n_workers);

  spdlog::debug("constructed prefixes");

  return GetSPSS<K, N, KeyType>(unitigs, prefixes, n_workers, n_buckets);
}

// Constructs a small-weight SPSS from unitigs.
template <int K, int N, typename KeyType>
std::vector<std::string> GetSPSSCanonical(
    const std::vector<std::string>& unitigs,
    const absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>>& prefixes,
    const absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>>& suffixes,
    bool fast, int n_workers, int n_buckets = 64) {
  const std::int64_t n = unitigs.size();

  // Considers a graph where node i represents unitigs[i].
  // Each node has two sides just as the one considered in
  // GetUnitigsCanonical().

  // If the second element is true, the edge connects the same side of two
  // nodes.
  using Edge = std::pair<std::int64_t, bool>;

  // Returns all the edges incident to the right side of ith node.
  // Self loops are ignored.
  const auto GetEdgesRight = [&](std::int64_t i) {
    std::vector<Edge> edges;

    const std::string& unitig = unitigs[i];

    const Kmer<K> suffix(unitig.substr(unitig.length() - K, K));

    for (const Kmer<K>& suffix_next : suffix.Nexts()) {
      {
        auto it = prefixes.find(suffix_next);

        if (it != prefixes.end()) {
          for (std::int64_t j : it->second) {
            if (i == j) continue;

            // There is an edge from the right side from i to the left
            // side of j.

            edges.emplace_back(j, false);
          }
        }
      }

      const Kmer<K> suffix_next_complement = suffix_next.Complement();

      {
        auto it = suffixes.find(suffix_next_complement);

        if (it != suffixes.end()) {
          for (std::int64_t j : it->second) {
            if (i == j) continue;

            // There is an edge from the right side from i to the right side
            // of j.

            edges.emplace_back(j, true);
          }
        }
      }
    }

    return edges;
  };

  // Returns all the edges incident to the left side of ith node.
  // Self loops are ignored.
  const auto GetEdgesLeft = [&](std::int64_t i) {
    std::vector<Edge> edges;

    const std::string& unitig = unitigs[i];

    const Kmer<K> prefix(unitig.substr(0, K));

    for (const Kmer<K>& prefix_prev : prefix.Prevs()) {
      {
        auto it = suffixes.find(prefix_prev);

        if (it != suffixes.end()) {
          for (std::int64_t j : it->second) {
            if (i == j) continue;

            // There is an edge from the left side of i to the right side of
            // j.

            edges.emplace_back(j, false);
          }
        }
      }

      const Kmer<K> prefix_prev_complement = prefix_prev.Complement();

      {
        auto it = prefixes.find(prefix_prev_complement);

        if (it != prefixes.end()) {
          for (std::int64_t j : it->second) {
            if (i == j) continue;

            // There is an edge from the left side of i to the left side of
            // j.

            edges.emplace_back(j, true);
          }
        }
      }
    }

    return edges;
  };

  // edge_left[i] is the selected edge incident to the left side of i.
  absl::flat_hash_map<std::int64_t, Edge> edge_left;

  // edge_right[i] is the selected edge incident to the right side of i.
  absl::flat_hash_map<std::int64_t, Edge> edge_right;

  // If the second element of a pair is true, the compliment of the unitig
  // should be considered in concatenation.
  using Path = std::vector<std::pair<std::int64_t, bool>>;

  // If is_right_side is true, finds a path from the right side of start.
  // If is_right_side is false, finds a path from the left side of start.
  const auto FindPath = [&](std::int64_t start, bool is_right_side) {
    Path path;

    std::int64_t current = start;

    while (true) {
      bool is_same_side;

      if (is_right_side) {
        path.emplace_back(current, false);
        if (edge_right.find(current) == edge_right.end()) break;
        std::tie(current, is_same_side) = edge_right[current];
      } else {
        path.emplace_back(current, true);
        if (edge_left.find(current) == edge_left.end()) break;
        std::tie(current, is_same_side) = edge_left[current];
      }

      if (is_same_side) {
        is_right_side = !is_right_side;
      }
    }

    return path;
  };

  // Constructs a string by concatenating the unitigs in the path.
  const auto GetStringFromPath = [&](const Path& path) {
    assert(!path.empty());

    std::string s;
    bool is_first = true;

    for (const std::pair<std::int64_t, bool>& p : path) {
      if (is_first) {
        s += p.second ? internal::Complement(unitigs[p.first])
                      : unitigs[p.first];
        is_first = false;
      } else {
        s += p.second ? internal::Complement(unitigs[p.first])
                            .substr(K - 1, unitigs[p.first].length() - (K - 1))
                      : unitigs[p.first].substr(
                            K - 1, unitigs[p.first].length() - (K - 1));
      }
    }

    return s;
  };

  if (!fast) {
    // Constructs SPSS with one thread.

    spdlog::debug("constructing start, edge_left, and edge_right");

    for (std::int64_t i = 0; i < n; i++) {
      // Skips if it has already been visited.
      if (edge_left.find(i) != edge_left.end() ||
          edge_right.find(i) != edge_right.end())
        continue;

      std::int64_t current = i;

      bool is_right_side;

      {
        std::vector<Edge> edges_right = GetEdgesRight(current);
        std::vector<Edge> edges_left = GetEdgesLeft(current);

        // Skips if the node is isolated.
        if (edges_right.empty() && edges_left.empty()) continue;

        is_right_side = !edges_right.empty();
      }

      while (true) {
        if (is_right_side) {
          // Stops if the right side of the node was already used.
          if (edge_right.find(current) != edge_right.end()) break;

          std::vector<Edge> edges_right = GetEdgesRight(current);

          // Stops if there is no edges on the right side.
          if (edges_right.empty()) break;

          std::int64_t next;
          bool is_same_side;
          bool found = false;

          for (auto it = edges_right.begin(); it != edges_right.end(); ++it) {
            std::tie(next, is_same_side) = *it;

            // Skips if it creates a loop.
            if (next == i) continue;

            // Skips if the side of the next node was already used.
            if (is_same_side) {
              if (edge_right.find(next) != edge_right.end()) continue;
            } else {
              if (edge_left.find(next) != edge_left.end()) continue;
            }

            found = true;
            break;
          }

          if (!found) break;

          // Extends the path to next.

          edge_right[current] = std::make_pair(next, is_same_side);

          if (is_same_side) {
            edge_right[next] = std::make_pair(current, is_same_side);

            is_right_side = !is_right_side;
          } else {
            edge_left[next] = std::make_pair(current, is_same_side);
          }

          current = next;
        } else {
          // Stops if the left side of the node was already used.
          if (edge_left.find(current) != edge_left.end()) break;

          std::vector<Edge> edges_left = GetEdgesLeft(current);

          // Stops if there is no edges on the left side.
          if (edges_left.empty()) break;

          std::int64_t next;
          bool is_same_side;
          bool found = false;

          for (auto it = edges_left.begin(); it != edges_left.end(); ++it) {
            std::tie(next, is_same_side) = *it;

            // Skips if it creates a loop.
            if (next == i) continue;

            // Skips if the side of the next node was already used.
            if (is_same_side) {
              if (edge_left.find(next) != edge_left.end()) continue;
            } else {
              if (edge_right.find(next) != edge_right.end()) continue;
            }

            found = true;
            break;
          }

          if (!found) break;

          // Extends the path to next.

          edge_left[current] = std::make_pair(next, is_same_side);

          if (is_same_side) {
            edge_left[next] = std::make_pair(current, is_same_side);

            is_right_side = !is_right_side;
          } else {
            edge_right[next] = std::make_pair(current, is_same_side);
          }

          current = next;
        }
      }
    }

    spdlog::debug("constructed start, edge_left, and edge_right");

    std::vector<std::string> spss;

    spdlog::debug("constructing spss");

    for (std::int64_t i = 0; i < n; i++) {
      const bool has_left = edge_left.find(i) != edge_left.end();
      const bool has_right = edge_right.find(i) != edge_right.end();

      if (has_left && has_right) continue;

      Path path;

      if (has_right) {
        path = FindPath(i, true);
      } else {
        path = FindPath(i, false);
      }

      if (path.front().first <= path.back().first) {
        spss.push_back(GetStringFromPath(path));
      }
    }

    spdlog::debug("constructed spss");

    return spss;
  }

  {
    spdlog::debug("constructing edge_left and edge_right");

    // Nodes are divided into n_buckets buckets to allow parallel processing.

    std::vector<std::thread> threads;
    std::vector<std::mutex> mus(n_buckets);
    std::vector<absl::flat_hash_map<std::int64_t, Edge>> buf_edge_left(
        n_buckets);
    std::vector<absl::flat_hash_map<std::int64_t, Edge>> buf_edge_right(
        n_buckets);

    // Acquires locks for node i and j.
    const auto AcquireLock = [&](std::int64_t i, std::int64_t j) {
      int bucket_i = i % n_buckets;
      int bucket_j = j % n_buckets;

      if (bucket_i == bucket_j) {
        mus[bucket_i].lock();
        return;
      }

      mus[std::min(bucket_i, bucket_j)].lock();
      mus[std::max(bucket_i, bucket_j)].lock();
    };

    // Releases locks for node i and j.
    const auto ReleaseLock = [&](std::int64_t i, std::int64_t j) {
      int bucket_i = i % n_buckets;
      int bucket_j = j % n_buckets;

      if (bucket_i == bucket_j) {
        mus[bucket_i].unlock();
        return;
      }

      mus[std::max(bucket_i, bucket_j)].unlock();
      mus[std::min(bucket_i, bucket_j)].unlock();
    };

    // Returns true if an edge is incident to the left side of i.
    const auto HasEdgeLeft = [&](std::int64_t i) {
      int bucket = i % n_buckets;
      return buf_edge_left[bucket].find(i) != buf_edge_left[bucket].end();
    };

    // Returns true if an edge is incident to the right side of i.
    const auto HasEdgeRight = [&](std::int64_t i) {
      int bucket = i % n_buckets;
      return buf_edge_right[bucket].find(i) != buf_edge_right[bucket].end();
    };

    // Adds an edge incident to the left side of i.
    const auto AddEdgeLeft = [&](std::int64_t i, std::int64_t j,
                                 bool is_same_side) {
      int bucket_i = i % n_buckets;
      buf_edge_left[bucket_i][i] = {j, is_same_side};
    };

    // Adds an edge incident to the right side of i.
    const auto AddEdgeRight = [&](std::int64_t i, std::int64_t j,
                                  bool is_same_side) {
      int bucket_i = i % n_buckets;
      buf_edge_right[bucket_i][i] = {j, is_same_side};
    };

    for (const Range& range : Range(0, n).Split(n_workers)) {
      threads.emplace_back([&, range] {
        for (std::int64_t i : range) {
          for (const Edge& edge : GetEdgesRight(i)) {
            std::int64_t j;
            bool is_same_side;
            std::tie(j, is_same_side) = edge;

            AcquireLock(i, j);

            if (is_same_side) {
              if (!HasEdgeRight(i) && !HasEdgeRight(j)) {
                AddEdgeRight(i, j, true);
                AddEdgeRight(j, i, true);
              }
            } else {
              if (!HasEdgeRight(i) && !HasEdgeLeft(j)) {
                AddEdgeRight(i, j, false);
                AddEdgeLeft(j, i, false);
              }
            }

            ReleaseLock(i, j);
          }

          for (const Edge& edge : GetEdgesLeft(i)) {
            std::int64_t j;
            bool is_same_side;
            std::tie(j, is_same_side) = edge;

            AcquireLock(i, j);

            if (is_same_side) {
              if (!HasEdgeLeft(i) && !HasEdgeLeft(j)) {
                AddEdgeLeft(i, j, true);
                AddEdgeLeft(j, i, true);
              }
            } else {
              if (!HasEdgeLeft(i) && !HasEdgeRight(j)) {
                AddEdgeLeft(i, j, false);
                AddEdgeRight(j, i, false);
              }
            }

            ReleaseLock(i, j);
          }
        }
      });
    }

    for (std::thread& t : threads) t.join();

    spdlog::debug("moving from buffers");

    // Moves from buffers.
    {
      boost::asio::thread_pool pool(n_workers);

      boost::asio::post(pool, [&] {
        {
          std::int64_t size = 0;
          for (int i = 0; i < n_buckets; i++) {
            size += buf_edge_left[i].size();
          }
          edge_left.reserve(size);
        }

        for (int i = 0; i < n_buckets; i++) {
          edge_left.insert(buf_edge_left[i].begin(), buf_edge_left[i].end());
        }
      });

      boost::asio::post(pool, [&] {
        {
          std::int64_t size = 0;
          for (int i = 0; i < n_buckets; i++) {
            size += buf_edge_right[i].size();
          }
          edge_right.reserve(size);
        }

        for (int i = 0; i < n_buckets; i++) {
          edge_right.insert(buf_edge_right[i].begin(), buf_edge_right[i].end());
        }
      });

      pool.join();
    }

    spdlog::debug("constructed edge_left and edge_right");
  }

  {
    ParallelDisjointSet disjoint_set(n);

    spdlog::debug("constructing disjoint_set");

    {
      std::vector<std::thread> threads;

      for (const Range& range : Range(0, n).Split(n_workers)) {
        threads.emplace_back([&, range] {
          for (std::int64_t i : range) {
            {
              auto it = edge_left.find(i);

              if (it != edge_left.end()) {
                disjoint_set.Unite(i, it->second.first);
              }
            }

            {
              auto it = edge_right.find(i);

              if (it != edge_right.end()) {
                disjoint_set.Unite(i, it->second.first);
              }
            }
          }
        });
      }

      for (std::thread& t : threads) t.join();
    }

    spdlog::debug("constructed disjoint_set");

    spdlog::debug("removing loops");

    {
      absl::flat_hash_set<int> groups;
      absl::flat_hash_set<int> groups_with_terminals;

      spdlog::debug("constructing groups and groups_with_terminals");

      std::vector<std::thread> threads;
      std::mutex mu_groups;
      std::mutex mu_groups_with_terminals;

      for (const Range& range : Range(0, n).Split(n_workers)) {
        threads.emplace_back([&, range] {
          absl::flat_hash_set<int> buf_groups;
          absl::flat_hash_set<int> buf_groups_with_terminals;

          for (std::int64_t i : range) {
            int group = disjoint_set.Find(i);

            buf_groups.insert(group);

            if (edge_left.find(i) == edge_left.end() ||
                edge_right.find(i) == edge_right.end()) {
              buf_groups_with_terminals.insert(group);
            }
          }

          if (n_workers == 1) {
            groups = std::move(buf_groups);
            groups_with_terminals = std::move(buf_groups_with_terminals);
          } else {
            {
              std::lock_guard lck(mu_groups);
              groups.insert(buf_groups.begin(), buf_groups.end());
            }

            {
              std::lock_guard lck(mu_groups_with_terminals);
              groups_with_terminals.insert(buf_groups_with_terminals.begin(),
                                           buf_groups_with_terminals.end());
            }
          }
        });
      }

      for (std::thread& t : threads) t.join();

      spdlog::debug("constructed groups and groups_with_terminals");

      for (int i : groups) {
        if (groups_with_terminals.find(i) == groups_with_terminals.end()) {
          auto it = edge_left.find(i);
          assert(it != edge_left.end());

          int j;
          bool is_same_side;
          std::tie(j, is_same_side) = it->second;

          edge_left.erase(i);

          if (is_same_side) {
            edge_left.erase(j);
          } else {
            edge_right.erase(j);
          }
        }
      }
    }

    spdlog::debug("removed loops");
  }

  // Nodes without edges on the left side.
  std::vector<std::int64_t> terminals_left;
  // Nodes without edges on the right side.
  std::vector<std::int64_t> terminals_right;
  // Nodes without edges on both sides.
  std::vector<std::int64_t> terminals_both;

  {
    spdlog::debug(
        "constructing terminals_left, terminals_right, and terminals_both");

    std::vector<std::thread> threads;
    std::mutex mu_terminals_left;
    std::mutex mu_terminals_right;
    std::mutex mu_terminals_both;

    for (const Range& range : Range(0, n).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<std::int64_t> buf_terminals_left;
        std::vector<std::int64_t> buf_terminals_right;
        std::vector<std::int64_t> buf_terminals_both;

        for (std::int64_t i : range) {
          const bool has_left = edge_left.find(i) != edge_left.end();
          const bool has_right = edge_right.find(i) != edge_right.end();

          if (!has_left && !has_right) {
            buf_terminals_both.push_back(i);
          } else if (!has_left) {
            buf_terminals_left.push_back(i);
          } else if (!has_right) {
            buf_terminals_right.push_back(i);
          }
        }

        // Moves from buffers.

        if (n_workers == 1) {
          terminals_left = std::move(buf_terminals_left);
          terminals_right = std::move(buf_terminals_right);
          terminals_both = std::move(buf_terminals_both);
        } else {
          bool done_terminals_left = false;
          bool done_terminals_right = false;
          bool done_terminals_both = false;

          while (!done_terminals_left || !done_terminals_right ||
                 !done_terminals_both) {
            if (!done_terminals_left && mu_terminals_left.try_lock()) {
              terminals_left.insert(terminals_left.end(),
                                    buf_terminals_left.begin(),
                                    buf_terminals_left.end());
              mu_terminals_left.unlock();
              done_terminals_left = true;
            }

            if (!done_terminals_right && mu_terminals_right.try_lock()) {
              terminals_right.insert(terminals_right.end(),
                                     buf_terminals_right.begin(),
                                     buf_terminals_right.end());
              mu_terminals_right.unlock();
              done_terminals_right = true;
            }

            if (!done_terminals_both && mu_terminals_both.try_lock()) {
              terminals_both.insert(terminals_both.end(),
                                    buf_terminals_both.begin(),
                                    buf_terminals_both.end());
              mu_terminals_both.unlock();
              done_terminals_both = true;
            }
          }
        }
      });
    }

    for (std::thread& t : threads) t.join();

    spdlog::debug(
        "constructed terminals_left, terminals_right, and terminals_both");
  }

  std::vector<std::string> spss;

  {
    spdlog::debug("constructing spss");

    // Considers paths from terminals_left.
    {
      std::vector<std::thread> threads;
      std::mutex mu_spss;

      for (const Range& range :
           Range(0, terminals_left.size()).Split(n_workers)) {
        threads.emplace_back([&, range] {
          std::vector<std::string> buf_spss;

          for (std::int64_t i : range) {
            Path path = FindPath(terminals_left[i], true);

            if (path.front().first > path.back().first) continue;

            if (n_workers == 1) {
              spss.push_back(GetStringFromPath(path));
            } else {
              buf_spss.push_back(GetStringFromPath(path));
            }
          }

          if (n_workers != 1) {
            std::lock_guard lck(mu_spss);
            spss.insert(spss.end(), buf_spss.begin(), buf_spss.end());
          }
        });
      }

      for (std::thread& t : threads) t.join();
    }

    // Considers paths from terminals_right.
    {
      std::vector<std::thread> threads;
      std::mutex mu_spss;

      for (const Range& range :
           Range(0, terminals_right.size()).Split(n_workers)) {
        threads.emplace_back([&, range] {
          std::vector<std::string> buf_spss;

          for (std::int64_t i : range) {
            Path path = FindPath(terminals_right[i], false);

            if (path.front().first > path.back().first) continue;

            if (n_workers == 1) {
              spss.push_back(GetStringFromPath(path));
            } else {
              buf_spss.push_back(GetStringFromPath(path));
            }
          }

          if (n_workers != 1) {
            std::lock_guard lck(mu_spss);
            spss.insert(spss.end(), buf_spss.begin(), buf_spss.end());
          }
        });
      }

      for (std::thread& t : threads) t.join();
    }

    // Considers nodes in terminal_both.
    {
      std::vector<std::thread> threads;
      std::mutex mu_spss;

      for (const Range& range :
           Range(0, terminals_both.size()).Split(n_workers)) {
        threads.emplace_back([&, range] {
          std::vector<std::string> buf_spss;

          for (std::int64_t i : range) {
            if (n_workers == 1) {
              spss.push_back(unitigs[terminals_both[i]]);
            } else {
              buf_spss.push_back(unitigs[terminals_both[i]]);
            }
          }

          if (n_workers != 1) {
            std::lock_guard lck(mu_spss);
            spss.insert(spss.end(), buf_spss.begin(), buf_spss.end());
          }
        });
      }

      for (std::thread& t : threads) t.join();
    }

    spdlog::debug("constructed spss");
  }

  return spss;
}

// Constructs a small-weight SPSS from a set of canonical kmers.
template <int K, int N, typename KeyType>
std::vector<std::string> GetSPSSCanonical(
    const KmerSet<K, N, KeyType>& kmer_set, bool fast, int n_workers,
    int n_buckets = 64) {
  spdlog::debug("constructing unitigs");

  const std::vector<std::string> unitigs =
      GetUnitigsCanonical(kmer_set, n_workers);

  spdlog::debug("constructed unitigs");

  spdlog::debug("constructing prefixes and suffixes");

  absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>> prefixes =
      GetPrefixesFromUnitigs<K>(unitigs, n_workers);

  absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>> suffixes =
      GetSuffixesFromUnitigs<K>(unitigs, n_workers);

  spdlog::debug("constructed prefixes and suffixes");

  return GetSPSSCanonical<K, N, KeyType>(unitigs, prefixes, suffixes, fast,
                                         n_workers, n_buckets);
}

// Given a set of kmers, returns the lower bound of SPSS weight.
// If factor is 0.8, 80% of kmers are used to estimate the correct value.
template <int K, int N, typename KeyType>
std::int64_t GetSPSSWeight(const KmerSet<K, N, KeyType>& kmer_set,
                           int n_workers, float factor = 1.0) {
  std::vector<Kmer<K>> kmers;

  {
    absl::InsecureBitGen bitgen;

    if (factor == 1.0) {
      kmers = kmer_set.Find(n_workers);
    } else {
      kmers = kmer_set.Find(
          [&](const Kmer<K>&) {
            return absl::Uniform(bitgen, 0.0, 1.0) < factor;
          },
          n_workers,
          static_cast<std::int64_t>(static_cast<float>(kmer_set.Size()) *
                                    factor));
    }
  }

  std::atomic_int64_t n_iso = 0;
  std::atomic_int64_t n_dead = 0;
  std::atomic_int64_t n_sp = 0;

  const auto GetNeighborsOut = [&](const Kmer<K>& kmer) {
    std::vector<Kmer<K>> neighbors;

    for (const Kmer<K>& next : kmer.Nexts()) {
      if (kmer_set.Contains(next)) {
        neighbors.push_back(next);
      }
    }

    return neighbors;
  };

  const auto GetNeighborsIn = [&](const Kmer<K>& kmer) {
    std::vector<Kmer<K>> neighbors;

    for (const Kmer<K>& prev : kmer.Prevs()) {
      if (kmer_set.Contains(prev)) {
        neighbors.push_back(prev);
      }
    }

    return neighbors;
  };

  std::vector<std::thread> threads;

  for (const Range& range : Range(0, kmers.size()).Split(n_workers)) {
    for (std::int64_t i : range) {
      const Kmer<K>& kmer = kmers[i];

      std::vector<Kmer<K>> neighbors_in = GetNeighborsIn(kmer);
      std::vector<Kmer<K>> neighbors_out = GetNeighborsOut(kmer);

      if (neighbors_in.size() == 0 && neighbors_out.size() == 0) {
        n_iso += 1;
        continue;
      }

      if (neighbors_in.size() == 0 || neighbors_out.size() == 0) {
        n_dead += 1;
      }

      int n_special_neighbors_in = 0;
      int n_special_neighbors_out = 0;

      for (const Kmer<K>& neighbor : neighbors_in) {
        if (GetNeighborsOut(neighbor).size() == 1) {
          n_special_neighbors_in += 1;
        }
      }

      for (const Kmer<K>& neighbor : neighbors_out) {
        if (GetNeighborsIn(neighbor).size() == 1) {
          n_special_neighbors_out += 1;
        }
      }

      n_sp += std::max(0, n_special_neighbors_in - 1);
      n_sp += std::max(0, n_special_neighbors_out - 1);
    }
  }

  for (std::thread& t : threads) t.join();

  std::int64_t weight =
      ((n_dead + n_sp + 1) / 2 + n_iso) * (K - 1) + kmers.size();

  if (factor != 1.0) {
    weight = static_cast<std::int64_t>(static_cast<float>(weight) / factor);
  }

  return weight;
}

// Given a set of canonical kmers, returns the lower bound of the SPSS weight.
// If factor is 0.8, 80% of kmers are used to estimate the correct value.
template <int K, int N, typename KeyType>
std::int64_t GetSPSSWeightCanonical(const KmerSet<K, N, KeyType>& kmer_set,
                                    int n_workers, float factor = 1.0) {
  std::vector<Kmer<K>> kmers;
  absl::InsecureBitGen bitgen;

  if (factor == 1.0) {
    kmers = kmer_set.Find(n_workers);
  } else {
    kmers = kmer_set.Find(
        [&](const Kmer<K>&) {
          return absl::Uniform(bitgen, 0.0, 1.0) < factor;
        },
        n_workers,
        static_cast<std::int64_t>(static_cast<float>(kmer_set.Size()) *
                                  factor));
  }

  std::atomic_int64_t n_iso = 0;
  std::atomic_int64_t n_dead = 0;
  std::atomic_int64_t n_sp = 0;

  using Neighbor = std::pair<Kmer<K>, bool>;

  const auto GetNeighborsLeft = [&](const Kmer<K>& kmer) {
    std::vector<Neighbor> neighbors;

    for (const Kmer<K>& prev : kmer.Prevs()) {
      if (kmer_set.Contains(prev)) {
        neighbors.emplace_back(prev, false);
      }

      const Kmer<K> prev_complement = prev.Complement();

      if (kmer_set.Contains(prev_complement)) {
        neighbors.emplace_back(prev_complement, true);
      }
    }

    return neighbors;
  };

  const auto GetNeighborsRight = [&](const Kmer<K>& kmer) {
    std::vector<Neighbor> neighbors;

    for (const Kmer<K>& next : kmer.Nexts()) {
      if (kmer_set.Contains(next)) {
        neighbors.emplace_back(next, false);
      }

      const Kmer<K> next_complement = next.Complement();

      if (kmer_set.Contains(next_complement)) {
        neighbors.emplace_back(next_complement, true);
      }
    }

    return neighbors;
  };

  std::vector<std::thread> threads;

  for (const Range& range : Range(0, kmers.size()).Split(n_workers)) {
    threads.emplace_back([&, range] {
      for (std::int64_t i : range) {
        const Kmer<K>& kmer = kmers[i];

        std::vector<Neighbor> neighbors_left = GetNeighborsLeft(kmer);
        std::vector<Neighbor> neighbors_right = GetNeighborsRight(kmer);

        if (neighbors_left.size() + neighbors_right.size() == 0) {
          n_iso += 1;
          continue;
        }

        if (neighbors_left.size() == 0 || neighbors_right.size() == 0) {
          n_dead += 1;
        }

        int n_special_neighbors_left = 0;
        int n_special_neighbors_right = 0;

        for (const Neighbor& neighbor : neighbors_left) {
          if (neighbor.second) {
            if (GetNeighborsLeft(neighbor.first).size() == 1) {
              n_special_neighbors_left += 1;
            }
          } else {
            if (GetNeighborsRight(neighbor.first).size() == 1) {
              n_special_neighbors_left += 1;
            }
          }
        }

        n_sp += std::max(0, n_special_neighbors_left - 1);

        for (const Neighbor& neighbor : neighbors_right) {
          if (neighbor.second) {
            if (GetNeighborsRight(neighbor.first).size() == 1) {
              n_special_neighbors_right += 1;
            }
          } else {
            if (GetNeighborsLeft(neighbor.first).size() == 1) {
              n_special_neighbors_right += 1;
            }
          }
        }

        n_sp += std::max(0, n_special_neighbors_right - 1);
      }
    });
  }

  for (std::thread& t : threads) t.join();

  std::int64_t weight =
      ((n_dead + n_sp + 1) / 2 + n_iso) * (K - 1) + kmers.size();

  if (factor != 1.0) {
    weight = static_cast<std::int64_t>(static_cast<float>(weight) / factor);
  }

  return weight;
}

// Reads SPSS and returns a list of kmers.
template <int K>
std::vector<Kmer<K>> GetKmersFromSPSS(const std::vector<std::string>& spss,
                                      bool canonical, int n_workers) {
  std::vector<Kmer<K>> kmers;
  std::vector<std::thread> threads;
  std::mutex mu;

  if (n_workers != 1) {
    std::int64_t size = 0;

    for (const std::string& s : spss) {
      size += s.length() - K + 1;
    }

    kmers.reserve(size);
  }

  for (const Range& range : Range(0, spss.size()).Split(n_workers)) {
    threads.emplace_back([&, range] {
      std::vector<Kmer<K>> buf;

      {
        std::int64_t size = 0;

        for (std::int64_t i : range) {
          const std::string& s = spss[i];
          size += s.length() - K + 1;
        }

        buf.reserve(size);
      }

      for (std::int64_t i : range) {
        const std::string& s = spss[i];

        for (int j = 0; j < static_cast<int>(s.length()) - K + 1; j++) {
          Kmer<K> kmer(s.substr(j, K));
          if (canonical) kmer = kmer.Canonical();
          buf.push_back(kmer);
        }
      }

      if (n_workers == 1) {
        kmers = std::move(buf);
      } else {
        std::lock_guard lck(mu);
        kmers.insert(kmers.end(), buf.begin(), buf.end());
      }
    });
  }

  for (std::thread& t : threads) t.join();

  return kmers;
}

// Reads SPSS and returns the corresponding kmer set.
template <int K, int N, typename KeyType>
KmerSet<K, N, KeyType> GetKmerSetFromSPSS(const std::vector<std::string>& spss,
                                          bool canonical, int n_workers) {
  std::vector<std::thread> threads;

  KmerSet<K, N, KeyType> kmer_set;
  std::mutex mu;
  int done_count = 0;

  for (const Range& range : Range(0, spss.size()).Split(n_workers)) {
    threads.emplace_back([&, range] {
      KmerSet<K, N, KeyType> buf;

      for (std::int64_t i : range) {
        const std::string& s = spss[i];
        for (int j = 0; j < static_cast<int>(s.length()) - K + 1; j++) {
          Kmer<K> kmer(s.substr(j, K));
          if (canonical) kmer = kmer.Canonical();
          buf.Add(kmer);
        }
      }

      std::lock_guard lck(mu);
      kmer_set.Add(buf, 1 + done_count);
      done_count += 1;
    });
  }

  for (std::thread& t : threads) t.join();

  return kmer_set;
}

#endif
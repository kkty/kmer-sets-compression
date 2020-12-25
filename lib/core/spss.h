#ifndef CORE_SPSS_H_
#define CORE_SPSS_H_

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "core/kmer_set.h"
#include "core/range.h"
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
  const std::vector<Kmer<K>> terminals_left = kmer_set.Find(
      [&](const Kmer<K>& kmer) {
        return IsTerminalLeft(kmer) && !IsTerminalRight(kmer);
      },
      n_workers);

  spdlog::debug("constructed terminals_left");

  spdlog::debug("constructing terminals_right");

  // kmers that has no mates on the right side, but has a mate on the left side.
  const std::vector<Kmer<K>> terminals_right = kmer_set.Find(
      [&](const Kmer<K>& kmer) {
        return !IsTerminalLeft(kmer) && IsTerminalRight(kmer);
      },
      n_workers);

  spdlog::debug("constructed terminals_right");

  spdlog::debug("constructing terminals_both");

  // kmers that has no mates on the both sides.
  const std::vector<Kmer<K>> terminals_both = kmer_set.Find(
      [&](const Kmer<K>& kmer) {
        return IsTerminalLeft(kmer) && IsTerminalRight(kmer);
      },
      n_workers);

  spdlog::debug("constructed terminals_both");

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

    unitigs.push_back(internal::ConcatenateKmers(path));
  }

  spdlog::debug("processed non-branching loops");

  return unitigs;
}

// Constructs a small-weight SPSS from unitigs.
template <int K, int N, typename KeyType>
std::vector<std::string> GetSPSSCanonical(
    const std::vector<std::string>& unitigs, bool fast, int n_workers,
    int n_buckets = 64) {
  const std::int64_t n = unitigs.size();

  // Considers a graph where node i represents unitigs[i].
  // Each node has two sides just as the one considered in
  // GetUnitigsCanonical().

  // If i is in prefixes[kmer], the unitigs.substr(0, K) == kmer.String().
  absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>> prefixes;

  // Similar to prefixes, but suffixes are considered.
  absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>> suffixes;

  {
    spdlog::debug("constructing prefixes and suffixes");

    std::vector<std::thread> threads;
    std::mutex mu_prefixes;
    std::mutex mu_suffixes;

    for (const Range& range : Range(0, n).Split(n_workers)) {
      threads.emplace_back([&, range] {
        absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>> buf_prefixes;
        absl::flat_hash_map<Kmer<K>, std::vector<std::int64_t>> buf_suffixes;

        for (std::int64_t i : range) {
          const std::string& unitig = unitigs[i];
          const Kmer<K> prefix(unitig.substr(0, K));
          const Kmer<K> suffix(unitig.substr(unitig.length() - K, K));
          buf_prefixes[prefix].push_back(i);
          buf_suffixes[suffix].push_back(i);
        }

        // Moves from buffers.

        if (n_workers == 1) {
          prefixes = std::move(buf_prefixes);
          suffixes = std::move(buf_suffixes);
        } else {
          bool done_prefixes = false;
          bool done_suffixes = false;

          while (!done_prefixes || !done_suffixes) {
            if (!done_prefixes && mu_prefixes.try_lock()) {
              prefixes.insert(buf_prefixes.begin(), buf_prefixes.end());
              mu_prefixes.unlock();
              done_prefixes = true;
            }

            if (!done_suffixes && mu_suffixes.try_lock()) {
              suffixes.insert(buf_suffixes.begin(), buf_suffixes.end());
              mu_suffixes.unlock();
              done_suffixes = true;
            }
          }
        }
      });
    }

    for (std::thread& t : threads) t.join();

    spdlog::debug("constructed prefixes and suffixes");
  }

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
      if (prefixes.find(suffix_next) != prefixes.end()) {
        for (std::int64_t j : prefixes[suffix_next]) {
          if (i == j) continue;

          // There is an edge from the right side from i to the left
          // side of j.

          edges.emplace_back(j, false);
        }
      }

      const Kmer<K> suffix_next_complement = suffix_next.Complement();

      if (suffixes.find(suffix_next_complement) != suffixes.end()) {
        for (std::int64_t j : suffixes[suffix_next_complement]) {
          if (i == j) continue;

          // There is an edge from the right side from i to the right side
          // of j.

          edges.emplace_back(j, true);
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
      if (suffixes.find(prefix_prev) != suffixes.end()) {
        for (std::int64_t j : suffixes[prefix_prev]) {
          if (i == j) continue;

          // There is an edge from the left side of i to the right side of
          // j.

          edges.emplace_back(j, false);
        }
      }

      const Kmer<K> prefix_prev_complement = prefix_prev.Complement();

      if (prefixes.find(prefix_prev_complement) != prefixes.end()) {
        for (std::int64_t j : prefixes[prefix_prev_complement]) {
          if (i == j) continue;

          // There is an edge from the left side of i to the left side of
          // j.

          edges.emplace_back(j, true);
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
    // NOLINTNEXTLINE
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
  absl::flat_hash_set<std::int64_t> visited;

  {
    spdlog::debug("constructing spss and visited");

    const auto MoveFromBuffer = [&](std::mutex& mu_spss, std::mutex& mu_visited,
                                    std::vector<std::string>& buf_spss,
                                    std::vector<std::int64_t>& buf_visited) {
      bool done_spss = false;
      bool done_visited = false;

      while (!done_spss || !done_visited) {
        if (!done_spss && mu_spss.try_lock()) {
          for (std::string& s : buf_spss) spss.push_back(std::move(s));
          mu_spss.unlock();
          done_spss = true;
        }

        if (!done_visited && mu_visited.try_lock()) {
          visited.insert(buf_visited.begin(), buf_visited.end());
          mu_visited.unlock();
          done_visited = true;
        }
      }
    };

    // Considers paths from terminals_left.
    {
      std::vector<std::thread> threads;
      std::mutex mu_spss;
      std::mutex mu_visited;

      for (const Range& range :
           Range(0, terminals_left.size()).Split(n_workers)) {
        threads.emplace_back([&, range] {
          std::vector<std::string> buf_spss;
          std::vector<std::int64_t> buf_visited;

          for (std::int64_t i : range) {
            Path path = FindPath(terminals_left[i], true);

            if (path.front().first > path.back().first) continue;

            for (std::size_t j = 0; j < path.size(); j++) {
              if (n_workers == 1) {
                visited.insert(path[j].first);
              } else {
                buf_visited.push_back(path[j].first);
              }
            }

            if (n_workers == 1) {
              spss.push_back(GetStringFromPath(path));
            } else {
              buf_spss.push_back(GetStringFromPath(path));
            }
          }

          if (n_workers != 1) {
            MoveFromBuffer(mu_spss, mu_visited, buf_spss, buf_visited);
          }
        });
      }

      for (std::thread& t : threads) t.join();
    }

    // Considers paths from terminals_right.
    {
      std::vector<std::thread> threads;
      std::mutex mu_spss;
      std::mutex mu_visited;

      for (const Range& range :
           Range(0, terminals_right.size()).Split(n_workers)) {
        threads.emplace_back([&, range] {
          std::vector<std::string> buf_spss;
          std::vector<std::int64_t> buf_visited;

          for (std::int64_t i : range) {
            Path path = FindPath(terminals_right[i], false);

            if (path.front().first > path.back().first) continue;

            for (std::size_t j = 0; j < path.size(); j++) {
              if (n_workers == 1) {
                visited.insert(path[j].first);
              } else {
                buf_visited.push_back(path[j].first);
              }
            }

            if (n_workers == 1) {
              spss.push_back(GetStringFromPath(path));
            } else {
              buf_spss.push_back(GetStringFromPath(path));
            }
          }

          if (n_workers != 1) {
            MoveFromBuffer(mu_spss, mu_visited, buf_spss, buf_visited);
          }
        });
      }

      for (std::thread& t : threads) t.join();
    }

    // Considers nodes in terminal_both.
    {
      std::vector<std::thread> threads;
      std::mutex mu_spss;
      std::mutex mu_visited;

      for (const Range& range :
           Range(0, terminals_both.size()).Split(n_workers)) {
        threads.emplace_back([&, range] {
          std::vector<std::string> buf_spss;
          std::vector<std::int64_t> buf_visited;

          for (std::int64_t i : range) {
            if (n_workers == 1) {
              visited.insert(terminals_both[i]);
              spss.push_back(unitigs[terminals_both[i]]);
            } else {
              buf_visited.push_back(terminals_both[i]);
              buf_spss.push_back(unitigs[terminals_both[i]]);
            }
          }

          if (n_workers != 1) {
            MoveFromBuffer(mu_spss, mu_visited, buf_spss, buf_visited);
          }
        });
      }

      for (std::thread& t : threads) t.join();
    }

    spdlog::debug("constructed spss and visited");
  }

  {
    spdlog::debug("processing non-branching loops");

    std::vector<std::int64_t> not_visited;

    // Constructs "not_visited".
    {
      std::vector<std::thread> threads;
      std::mutex mu;

      for (const Range& range : Range(0, n).Split(n_workers)) {
        threads.emplace_back([&, range] {
          std::vector<std::int64_t> buf;

          for (std::int64_t i : range) {
            if (visited.find(i) == visited.end()) {
              buf.push_back(i);
            }
          }

          std::lock_guard lck(mu);
          not_visited.insert(not_visited.end(), buf.begin(), buf.end());
        });
      }

      for (std::thread& t : threads) t.join();
    }

    for (std::int64_t start : not_visited) {
      if (visited.find(start) != visited.end()) continue;

      Path path;
      std::int64_t current = start;
      bool is_right_side = true;

      while (visited.find(current) == visited.end()) {
        visited.insert(current);
        path.emplace_back(current, !is_right_side);

        bool is_same_side;
        std::tie(current, is_same_side) =
            is_right_side ? edge_right[current] : edge_left[current];

        if (is_same_side) is_right_side = !is_right_side;
      }

      spss.push_back(GetStringFromPath(path));
    }

    spdlog::debug("processed non-branching loops");
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

  return GetSPSSCanonical<K, N, KeyType>(unitigs, fast, n_workers, n_buckets);
}

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
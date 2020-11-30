#ifndef CORE_UNITIGS_H_
#define CORE_UNITIGS_H_

#include <cassert>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "core/kmer_set.h"
#include "core/range.h"
#include "spdlog/spdlog.h"

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

// Returns the complement of s. The complement of s can be constructed by
// reversing s and replacing A, C, G and T with T, G, C and A respectively.
std::string Complement(std::string s) {
  std::reverse(s.begin(), s.end());

  for (size_t i = 0; i < s.length(); i++) {
    switch (s[i]) {
      case 'A':
        s[i] = 'T';
        break;
      case 'C':
        s[i] = 'G';
        break;
      case 'G':
        s[i] = 'C';
        break;
      case 'T':
        s[i] = 'A';
        break;
    }
  }

  return s;
}

// Constructs unitigs from a set of canonical kmers.
template <int K, typename KeyType>
std::vector<std::string> GetUnitigsCanonical(
    const KmerSet<K, KeyType>& kmer_set, int n_workers) {
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

// Constructs a small-weight SPSS from a set of canonical kmers.
template <int K, typename KeyType>
std::vector<std::string> GetSPSSCanonical(const KmerSet<K, KeyType>& kmer_set,
                                          int n_workers, int n_buckets = 64) {
  spdlog::debug("constructing unitigs");
  const std::vector<std::string> unitigs = GetUnitigs(kmer_set, n_workers);
  spdlog::debug("constructed unitigs");

  const int64_t n = unitigs.size();

  // Considers a graph where node i represents unitigs[i].
  // Each node has two sides just as the one considered in
  // GetUnitigsCanonical().

  // If i is in prefixes[kmer], the unitigs.substr(0, K) == kmer.String().
  absl::flat_hash_map<Kmer<K>, std::vector<int64_t>> prefixes;

  // Similar to prefixes, but suffixes are considered.
  absl::flat_hash_map<Kmer<K>, std::vector<int64_t>> suffixes;

  {
    spdlog::debug("constructing prefixes and suffixes");

    std::vector<std::thread> threads;
    std::mutex mu_prefixes;
    std::mutex mu_suffixes;

    for (const Range& range : Range(0, n).Split(n_workers)) {
      threads.emplace_back([&, range] {
        absl::flat_hash_map<Kmer<K>, std::vector<int64_t>> buf_prefixes;
        absl::flat_hash_map<Kmer<K>, std::vector<int64_t>> buf_suffixes;

        range.ForEach([&](int64_t i) {
          const std::string& unitig = unitigs[i];
          const Kmer<K> prefix(unitig.substr(0, K));
          const Kmer<K> suffix(unitig.substr(unitig.length() - K, K));
          buf_prefixes[prefix].push_back(i);
          buf_suffixes[suffix].push_back(i);
        });

        // Moves from buffers.

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
      });
    }

    for (std::thread& t : threads) t.join();

    spdlog::debug("constructed prefixes and suffixes");
  }

  // If the second element is true, the edge connects the same side of two
  // nodes.
  using Edge = std::pair<int64_t, bool>;

  // edge_left[i] is the edge incident to the left side of i.
  absl::flat_hash_map<int64_t, Edge> edges_left;

  // edge_right[i] is the edge incident to the right side of i.
  absl::flat_hash_map<int64_t, Edge> edges_right;

  {
    spdlog::debug("constructing edges_left and edges_right");

    // Nodes are divided into n_buckets buckets to allow parallel processing.

    std::vector<std::thread> threads;
    std::vector<std::mutex> mus(n_buckets);
    std::vector<absl::flat_hash_map<int64_t, Edge>> buf_edges_left(n_buckets);
    std::vector<absl::flat_hash_map<int64_t, Edge>> buf_edges_right(n_buckets);

    // Acquires locks for node i and j.
    const auto AcquireLock = [&](int64_t i, int64_t j) {
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
    const auto ReleaseLock = [&](int64_t i, int64_t j) {
      int bucket_i = i % n_buckets;
      int bucket_j = j % n_buckets;

      if (bucket_i == bucket_j) {
        mus[bucket_i].unlock();
        return;
      }

      mus[std::max(bucket_i, bucket_j)].unlock();
      mus[std::min(bucket_i, bucket_j)].unlock();
    };

    // Returns true if an edge is indicent to the left side of i.
    const auto HasEdgeLeft = [&](int64_t i) {
      int bucket = i % n_buckets;
      return buf_edges_left[bucket].find(i) != buf_edges_left[bucket].end();
    };

    // Returns true if an edge is indicent to the right side of i.
    const auto HasEdgeRight = [&](int64_t i) {
      int bucket = i % n_buckets;
      return buf_edges_right[bucket].find(i) != buf_edges_right[bucket].end();
    };

    // Adds an edge incident to the left side of i.
    const auto AddEdgeLeft = [&](int64_t i, int64_t j, bool is_same_side) {
      int bucket_i = i % n_buckets;
      buf_edges_left[bucket_i][i] = {j, is_same_side};
    };

    // Adds an edge incident to the right side of i.
    const auto AddEdgeRight = [&](int64_t i, int64_t j, bool is_same_side) {
      int bucket_i = i % n_buckets;
      buf_edges_right[bucket_i][i] = {j, is_same_side};
    };

    for (const Range& range : Range(0, n).Split(n_workers)) {
      threads.emplace_back([&, range] {
        range.ForEach([&](int64_t i) {
          const std::string& unitig = unitigs[i];

          const Kmer<K> prefix(unitig.substr(0, K));
          const Kmer<K> suffix(unitig.substr(unitig.length() - K, K));

          for (const Kmer<K>& suffix_next : suffix.Nexts()) {
            if (prefixes.find(suffix_next) != prefixes.end()) {
              for (int64_t j : prefixes[suffix_next]) {
                if (i == j) continue;

                // There is an edge from the right side from i to the left
                // side of j.

                AcquireLock(i, j);

                if (!HasEdgeRight(i) && !HasEdgeLeft(j)) {
                  AddEdgeRight(i, j, false);
                  AddEdgeLeft(j, i, false);
                }

                ReleaseLock(i, j);
              }
            }

            const Kmer<K> suffix_next_complement = suffix_next.Complement();

            if (suffixes.find(suffix_next_complement) != suffixes.end()) {
              for (int64_t j : suffixes[suffix_next_complement]) {
                if (i == j) continue;

                // There is an edge from the right side from i to the right side
                // of j.

                AcquireLock(i, j);

                if (!HasEdgeRight(i) && !HasEdgeRight(j)) {
                  AddEdgeRight(i, j, true);
                  AddEdgeRight(j, i, true);
                }

                ReleaseLock(i, j);
              }
            }
          }

          for (const Kmer<K>& prefix_prev : prefix.Prevs()) {
            if (suffixes.find(prefix_prev) != suffixes.end()) {
              for (int64_t j : suffixes[prefix_prev]) {
                if (i == j) continue;

                // There is an edge from the left side of i to the right side of
                // j.

                AcquireLock(i, j);

                if (!HasEdgeLeft(i) && !HasEdgeRight(j)) {
                  AddEdgeLeft(i, j, false);
                  AddEdgeRight(j, i, false);
                }

                ReleaseLock(i, j);
              }
            }

            const Kmer<K> prefix_prev_complement = prefix_prev.Complement();

            if (prefixes.find(prefix_prev_complement) != prefixes.end()) {
              for (int64_t j : prefixes[prefix_prev_complement]) {
                if (i == j) continue;

                // There is an edge from the left side of i to the left side of
                // j.

                AcquireLock(i, j);

                if (!HasEdgeLeft(i) && !HasEdgeLeft(j)) {
                  AddEdgeLeft(i, j, true);
                  AddEdgeLeft(j, i, true);
                }

                ReleaseLock(i, j);
              }
            }
          }
        });
      });
    }

    for (std::thread& t : threads) t.join();

    // Moves from buffers.
    {
      boost::asio::thread_pool pool(n_workers);

      boost::asio::post(pool, [&] {
        for (int i = 0; i < n_buckets; i++) {
          edges_left.insert(buf_edges_left[i].begin(), buf_edges_left[i].end());
        }
      });

      boost::asio::post(pool, [&] {
        for (int i = 0; i < n_buckets; i++) {
          edges_right.insert(buf_edges_right[i].begin(),
                             buf_edges_right[i].end());
        }
      });

      pool.join();
    }

    spdlog::debug("constructed edges_left and edges_right");
  }

  // Nodes without edges on the left side.
  std::vector<int64_t> terminals_left;
  // Nodes without edges on the right side.
  std::vector<int64_t> terminals_right;
  // Nodes without edges on both sides.
  std::vector<int64_t> terminals_both;

  {
    spdlog::debug(
        "constructing terminals_left, terminals_right, and terminals_both");

    std::vector<std::thread> threads;
    std::mutex mu_terminals_left;
    std::mutex mu_terminals_right;
    std::mutex mu_terminals_both;

    for (const Range& range : Range(0, n).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<int64_t> buf_terminals_left;
        std::vector<int64_t> buf_terminals_right;
        std::vector<int64_t> buf_terminals_both;

        range.ForEach([&](int64_t i) {
          const bool has_left = edges_left.find(i) != edges_left.end();
          const bool has_right = edges_right.find(i) != edges_right.end();

          if (!has_left && !has_right) {
            buf_terminals_both.push_back(i);
          } else if (!has_left) {
            buf_terminals_left.push_back(i);
          } else if (!has_right) {
            buf_terminals_right.push_back(i);
          }
        });

        // Moves from buffers.

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
      });
    }

    for (std::thread& t : threads) t.join();

    spdlog::debug(
        "constructed terminals_left, terminals_right, and terminals_both");
  }

  // If the second pair of a pair is true, the compliment of the unitig should
  // be considered when concatenating.
  using Path = std::vector<std::pair<int64_t, bool>>;

  // If is_right_side is true, finds a path from the right side of start.
  // If is_right_side is false, finds a path from the left side of start.
  const auto FindPath = [&](int64_t start, bool is_right_side) {
    Path path;

    int64_t current = start;

    while (true) {
      bool is_same_side;

      if (is_right_side) {
        path.emplace_back(current, false);
        if (edges_right.find(current) == edges_right.end()) break;
        std::tie(current, is_same_side) = edges_right[current];
      } else {
        path.emplace_back(current, true);
        if (edges_left.find(current) == edges_left.end()) break;
        std::tie(current, is_same_side) = edges_left[current];
      }

      if (is_same_side) is_right_side = !is_right_side;
    }

    return path;
  };

  const auto GetStringFromPath = [&](const Path& path) {
    assert(path.size() > 0);

    std::string s;
    bool is_first = true;

    for (const std::pair<int64_t, bool>& p : path) {
      if (is_first) {
        s += p.second ? Complement(unitigs[p.first]) : unitigs[p.first];
        is_first = false;
      } else {
        s += p.second ? Complement(unitigs[p.first])
                            .substr(K - 1, unitigs[p.first].length() - (K - 1))
                      : unitigs[p.first].substr(
                            K - 1, unitigs[p.first].length() - (K - 1));
      }
    }

    return s;
  };

  std::vector<std::string> spss;
  absl::flat_hash_set<int64_t> visited;

  {
    spdlog::debug("constructing spss and visited");

    const auto MoveFromBuffer = [&](std::mutex& mu_spss, std::mutex& mu_visited,
                                    std::vector<std::string>& buf_spss,
                                    std::vector<int64_t>& buf_visited) {
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
          std::vector<int64_t> buf_visited;

          range.ForEach([&](int64_t i) {
            Path path = FindPath(terminals_left[i], true);

            if (path.front().first > path.back().first) return;

            for (size_t j = 0; j < path.size(); j++) {
              buf_visited.push_back(path[j].first);
            }

            buf_spss.push_back(GetStringFromPath(path));
          });

          MoveFromBuffer(mu_spss, mu_visited, buf_spss, buf_visited);
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
          std::vector<int64_t> buf_visited;

          range.ForEach([&](int64_t i) {
            Path path = FindPath(terminals_right[i], false);

            if (path.front().first > path.back().first) return;

            for (size_t j = 0; j < path.size(); j++) {
              buf_visited.push_back(path[j].first);
            }

            buf_spss.push_back(GetStringFromPath(path));
          });

          MoveFromBuffer(mu_spss, mu_visited, buf_spss, buf_visited);
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
          std::vector<int64_t> buf_visited;

          range.ForEach([&](int64_t i) {
            buf_visited.push_back(terminals_both[i]);
            buf_spss.push_back(unitigs[terminals_both[i]]);
          });

          MoveFromBuffer(mu_spss, mu_visited, buf_spss, buf_visited);
        });
      }

      for (std::thread& t : threads) t.join();
    }

    spdlog::debug("constructed spss and visited");
  }

  {
    spdlog::debug("processing non-branching loops");

    std::vector<int64_t> not_visited;

    // Constructs "not_visited".
    {
      std::vector<std::thread> threads;
      std::mutex mu;

      for (const Range& range : Range(0, n).Split(n_workers)) {
        threads.emplace_back([&, range] {
          std::vector<int64_t> buf;

          range.ForEach([&](int64_t i) {
            if (visited.find(i) == visited.end()) {
              buf.push_back(i);
            }
          });

          std::lock_guard lck(mu);
          not_visited.insert(not_visited.end(), buf.begin(), buf.end());
        });
      }

      for (std::thread& t : threads) t.join();
    }

    for (int64_t start : not_visited) {
      if (visited.find(start) != visited.end()) continue;

      Path path;
      int64_t current = start;
      bool is_right_side = true;

      while (visited.find(current) == visited.end()) {
        visited.insert(current);
        path.emplace_back(current, !is_right_side);

        bool is_same_side;
        std::tie(current, is_same_side) =
            is_right_side ? edges_right[current] : edges_left[current];

        if (is_same_side) is_right_side = !is_right_side;
      }

      spss.push_back(GetStringFromPath(path));
    }

    spdlog::debug("processed non-branching loops");
  }

  return spss;
}

// Constructs unitigs from a kmer set.
template <int K, typename KeyType>
std::vector<std::string> GetUnitigs(const KmerSet<K, KeyType>& kmer_set,
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
            current = GetNexts(current)[0];
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
        current = GetNexts(current)[0];
      }

      unitigs.push_back(GetUnitigFromKmers(path));
    }
  }

  return unitigs;
}

#endif
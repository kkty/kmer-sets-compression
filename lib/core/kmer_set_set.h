#ifndef CORE_KMER_SET_SET_H_
#define CORE_KMER_SET_SET_H_

#include <algorithm>
#include <cmath>
#include <mutex>
#include <optional>
#include <queue>
#include <sstream>
#include <thread>
#include <tuple>
#include <variant>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/random/random.h"
#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "core/graph.h"
#include "core/io.h"
#include "core/kmer.h"
#include "core/kmer_set.h"
#include "core/kmer_set_compact.h"
#include "core/kmer_set_compressed.h"
#include "core/neighbor_joining.h"
#include "core/range.h"
#include "lemon/list_graph.h"
#include "lemon/matching.h"
#include "spdlog/spdlog.h"

// KmerSetSet can be used to represent multiple k-mer sets in less space.
// It is a recursive structure. That is, a KmerSetSet can contain another
// KmerSetSet internally.
template <int K, typename KeyType, typename CostFunctionType>
class KmerSetSet {
 public:
  // Not to have a recursive structure, set "recursion_limit" to 0.
  KmerSetSet(std::vector<KmerSet<K, KeyType>> kmer_sets, int recursion_limit,
             CostFunctionType cost_function, int n_workers)
      : n_(kmer_sets.size()) {
    // kmer_sets[n_] is an empty set.
    kmer_sets.push_back(KmerSet<K, KeyType>());

    std::tie(cost_, tree_) =
        ConstructMST(kmer_sets, n_, cost_function, n_workers);

    spdlog::debug("constructed MST: cost_ = {}", cost_);

    std::vector<KmerSet<K, KeyType>> diffs;

    for (int i = 0; i < n_; i++) {
      const int parent = tree_->Parent(i);

      diff_table_[i][parent] = diffs.size();
      diffs.push_back(Sub(kmer_sets[parent], kmer_sets[i], n_workers));

      diff_table_[parent][i] = diffs.size();
      diffs.push_back(Sub(kmer_sets[i], kmer_sets[parent], n_workers));
    }

    if (recursion_limit == 0) {
      diffs_ = std::move(diffs);
    } else {
      diffs_ = new KmerSetSet(std::move(diffs), recursion_limit - 1,
                              cost_function, n_workers);
    }
  }

  ~KmerSetSet() {
    if (!IsTerminal()) delete std::get<KmerSetSet*>(diffs_);
  }

  // Returns the ith k-mer set.
  KmerSet<K, KeyType> Get(int i, int n_workers) const {
    std::vector<int> path;

    while (true) {
      path.push_back(i);
      if (i == n_) break;
      i = tree_->Parent(i);
    }

    std::reverse(path.begin(), path.end());

    KmerSet<K, KeyType> kmer_set;

    // Traverses the path and reconstructs the set.
    for (size_t i = 0; i < path.size() - 1; i++) {
      if (IsTerminal()) {
        const std::vector<KmerSet<K, KeyType>>& diffs =
            std::get<std::vector<KmerSet<K, KeyType>>>(diffs_);

        kmer_set
            .Add(diffs[diff_table_.find(path[i])
                           ->second.find(path[i + 1])
                           ->second],
                 n_workers)
            .Sub(diffs[diff_table_.find(path[i + 1])
                           ->second.find(path[i])
                           ->second],
                 n_workers);
      } else {
        KmerSetSet* diffs = std::get<KmerSetSet*>(diffs_);

        KmerSet<K, KeyType> add;
        KmerSet<K, KeyType> sub;

        const auto calculate_add = [&](int n_workers) {
          add = diffs->Get(
              diff_table_.find(path[i])->second.find(path[i + 1])->second,
              n_workers);
        };

        const auto calculate_sub = [&](int n_workers) {
          sub = diffs->Get(
              diff_table_.find(path[i + 1])->second.find(path[i])->second,
              n_workers);
        };

        if (n_workers == 1) {
          calculate_add(1);
          calculate_sub(1);
        } else {
          std::vector<std::thread> threads;

          threads.emplace_back([&] { calculate_add(n_workers / 2); });

          threads.emplace_back(
              [&] { calculate_sub(n_workers - n_workers / 2); });

          for (std::thread& t : threads) t.join();
        }

        kmer_set.Add(add, n_workers).Sub(sub, n_workers);
      }
    }

    return kmer_set;
  }

  // Returns the number of k-mers stored internally.
  int64_t Size() const {
    if (!IsTerminal()) {
      return std::get<KmerSetSet*>(diffs_)->Size();
    }

    int64_t size = 0;

    for (const KmerSet<K, KeyType>& diff :
         std::get<std::vector<KmerSet<K, KeyType>>>(diffs_)) {
      size += diff.Size();
    }
    return size;
  }

  int64_t Cost() const {
    if (IsTerminal()) return cost_;

    return std::get<KmerSetSet*>(diffs_)->Cost();
  }

  // Uses another cost function to get the cost of the entire structure.
  // For each internal kmer set, cost_function(empty_kmer_set, kmer_set,
  // n_workers) is calculated and their sum is returned.
  template <typename F>
  int64_t Cost(F cost_function, int n_workers) const {
    if (!IsTerminal()) {
      return std::get<KmerSetSet*>(diffs_)->Cost(cost_function, n_workers);
    }

    int64_t cost = 0;

    const std::vector<KmerSet<K, KeyType>>& diffs =
        std::get<std::vector<KmerSet<K, KeyType>>>(diffs_);

    for (const KmerSet<K, KeyType>& diff : diffs) {
      cost += cost_function(KmerSet<K, KeyType>(), diff, n_workers);
    }

    return cost;
  }

 private:
  // Returns true if it is the end of recursion.
  bool IsTerminal() const {
    return std::holds_alternative<std::vector<KmerSet<K, KeyType>>>(diffs_);
  }

  static std::pair<int64_t, Tree> ConstructMST(
      const std::vector<KmerSet<K, KeyType>>& kmer_sets, int root,
      CostFunctionType cost_function, int n_workers) {
    BidirectionalGraph g;

    // Adds edges.
    {
      std::vector<std::pair<int, int>> pairs;

      for (size_t i = 0; i < kmer_sets.size(); i++) {
        for (size_t j = i + 1; j < kmer_sets.size(); j++) {
          pairs.emplace_back(i, j);
        }
      }

      std::vector<std::thread> threads;
      std::mutex mu;

      for (const Range& range : Range(0, pairs.size()).Split(n_workers)) {
        threads.emplace_back([&, range] {
          range.ForEach([&](int i) {
            std::pair<int, int> pair = pairs[i];

            int64_t cost =
                cost_function(kmer_sets[pair.first], kmer_sets[pair.second], 1);

            std::lock_guard _{mu};

            g.AddEdge(pair.first, pair.second, cost);
          });
        });
      }

      for (std::thread& t : threads) t.join();
    }

    return g.MST(root);
  }

  int n_;

  int64_t cost_;

  // Tree does not have a default constructor.
  std::optional<Tree> tree_;

  // If diffs_ is of type std::vector<KmerSet<K, KeyType>>,
  // diffs_[diff_table_[i][j]] represents the diff from ith data to jth data.
  // If diffs_ is of type KmerSetSet*, diffs_->Get(diff_table_[i][j]) represents
  // the diff from ith data to jth data.
  std::variant<std::vector<KmerSet<K, KeyType>>, KmerSetSet*> diffs_;
  absl::flat_hash_map<int, absl::flat_hash_map<int, int>> diff_table_;
};

// KmerSetSetNJ is similar to KmerSetSet, but it uses the neighbor joining
// algorithm internally.
template <int K, typename KeyType, typename CostFunctionType>
class KmerSetSetNJ {
 public:
  KmerSetSetNJ(std::vector<KmerSet<K, KeyType>> kmer_sets, int recursion_limit,
               CostFunctionType cost_function, int n_workers)
      : n_(kmer_sets.size()) {
    for (int i = 0; i < n_; i++) {
      spdlog::debug("kmer_sets[{}].Size() = {}", i, kmer_sets[i].Size());
    }

    // kmer_sets[n_] is an empty set.
    kmer_sets.push_back(KmerSet<K, KeyType>());

    neighbor_joining::SymmetricMatrix<float> all_distances;

    {
      std::vector<std::pair<int, int>> pairs;

      for (int i = 0; i < n_ + 1; i++) {
        for (int j = i; j < n_ + 1; j++) {
          pairs.emplace_back(i, j);
        }
      }

      std::vector<std::thread> threads;
      std::mutex mu;

      for (const Range& range : Range(0, pairs.size()).Split(n_workers)) {
        threads.emplace_back([&, range] {
          range.ForEach([&](int i) {
            const std::pair<int, int>& pair = pairs[i];

            int64_t cost = pair.first == pair.second
                               ? 0
                               : cost_function(kmer_sets[pair.first],
                                               kmer_sets[pair.second], 1);

            std::lock_guard _{mu};
            all_distances.Set(pair.first, pair.second, cost);
          });
        });
      }

      for (std::thread& thread : threads) thread.join();
    }

    const neighbor_joining::Result<float> result =
        neighbor_joining::Execute(all_distances, n_ + 1);

    absl::flat_hash_map<int, KmerSet<K, KeyType>> new_kmer_sets;

    for (int node : result.Nodes()) {
      // Finds newly added nodes.
      if (node > n_) {
        // Composes the kmer set corresponding to the node.
        KmerCounter<K, KeyType, float> kmer_counter;

        float denominator = 0;

        float distance_avg = 0;
        for (const auto& p : result.Distances(node)) {
          distance_avg += p.second;
        }
        distance_avg /= result.Nodes().size();

        for (const auto& [another_node, distance] : result.Distances(node)) {
          if (another_node <= n_) {
            KmerCounter<K, KeyType, float> another_kmer_counter =
                KmerCounter<K, KeyType, float>::FromSet(kmer_sets[another_node],
                                                        n_workers);
            float factor = std::pow(distance_avg / distance, 2);
            another_kmer_counter.Multiply(factor, n_workers);
            denominator += factor;
            kmer_counter.Add(another_kmer_counter, n_workers);
          }
        }

        KmerSet<K, KeyType> kmer_set;
        std::tie(kmer_set, std::ignore) =
            kmer_counter.ToSet(denominator / 2, n_workers);

        spdlog::debug("kmer_set.Size() = {}", kmer_set.Size());

        new_kmer_sets[node] = std::move(kmer_set);
      }
    }

    std::vector<KmerSet<K, KeyType>> diffs;

    // This will be used to determine which node is the parent of two adjacent
    // nodes.
    absl::flat_hash_map<int, float> distances_from_root = result.Distances(n_);

    for (int node : result.Nodes()) {
      if (node == n_) continue;
      for (int neighbor : result.Neighbors(node)) {
        if (distances_from_root[neighbor] < distances_from_root[node]) {
          const int child = node;
          const int parent = neighbor;

          parent_[child] = parent;

          const KmerSet<K, KeyType>& parent_kmer_set =
              parent <= n_ ? kmer_sets[parent] : new_kmer_sets[parent];
          const KmerSet<K, KeyType>& child_kmer_set =
              child <= n_ ? kmer_sets[child] : new_kmer_sets[child];

          cost_ += cost_function(parent_kmer_set, child_kmer_set, n_workers);

          {
            KmerSet<K, KeyType> sub =
                Sub(parent_kmer_set, child_kmer_set, n_workers);
            if (sub.Size() == 0) {
              diff_table_[child][parent] = -1;
            } else {
              diff_table_[child][parent] = diffs.size();
              diffs.push_back(std::move(sub));
            }
          }

          {
            KmerSet<K, KeyType> sub =
                Sub(child_kmer_set, parent_kmer_set, n_workers);
            if (sub.Size() == 0) {
              diff_table_[parent][child] = -1;
            } else {
              diff_table_[parent][child] = diffs.size();
              diffs.push_back(std::move(sub));
            }
          }
        }
      }
    }

    spdlog::debug("cost_ = {}", cost_);

    if (recursion_limit == 0) {
      diffs_ = std::move(diffs);
    } else {
      diffs_ = new KmerSetSetNJ(std::move(diffs), recursion_limit - 1,
                                cost_function, n_workers);
    }
  }

  ~KmerSetSetNJ() {
    if (!IsTerminal()) delete std::get<KmerSetSetNJ*>(diffs_);
  }

  int64_t Cost() const {
    if (IsTerminal()) return cost_;

    return std::get<KmerSetSetNJ*>(diffs_)->Cost();
  }

  // Reconstructs the ith kmer set.
  KmerSet<K, KeyType> Get(int i, int n_workers) const {
    std::vector<int> path;

    while (true) {
      path.push_back(i);
      if (i == n_) break;
      i = parent_.find(i)->second;
    }

    std::reverse(path.begin(), path.end());

    KmerSet<K, KeyType> kmer_set;

    for (size_t i = 0; i < path.size() - 1; i++) {
      if (IsTerminal()) {
        const std::vector<KmerSet<K, KeyType>>& diffs =
            std::get<std::vector<KmerSet<K, KeyType>>>(diffs_);

        {
          int pos = diff_table_.find(path[i])->second.find(path[i + 1])->second;
          if (pos != -1) {
            kmer_set.Add(diffs[pos], n_workers);
          }
        }

        {
          int pos = diff_table_.find(path[i + 1])->second.find(path[i])->second;
          if (pos != -1) {
            kmer_set.Sub(diffs[pos], n_workers);
          }
        }

      } else {
        KmerSetSetNJ* diffs = std::get<KmerSetSetNJ*>(diffs_);

        KmerSet<K, KeyType> add;
        KmerSet<K, KeyType> sub;

        const auto calculate_add = [&](int n_workers) {
          int pos = diff_table_.find(path[i])->second.find(path[i + 1])->second;
          if (pos != -1) add = diffs->Get(pos, n_workers);
        };

        const auto calculate_sub = [&](int n_workers) {
          int pos = diff_table_.find(path[i + 1])->second.find(path[i])->second;
          if (pos != -1) sub = diffs->Get(pos, n_workers);
        };

        if (n_workers == 1) {
          calculate_add(1);
          calculate_sub(1);
        } else {
          std::vector<std::thread> threads;

          threads.emplace_back([&] { calculate_add(n_workers / 2); });

          threads.emplace_back(
              [&] { calculate_sub(n_workers - n_workers / 2); });

          for (std::thread& thread : threads) thread.join();
        }

        kmer_set.Add(add, n_workers).Sub(sub, n_workers);
      }
    }

    return kmer_set;
  }

 private:
  // Returns true if it is the end of recursion.
  bool IsTerminal() const {
    return std::holds_alternative<std::vector<KmerSet<K, KeyType>>>(diffs_);
  }

  int n_;

  int64_t cost_ = 0;

  // parent_[i] is the parent of i.
  absl::flat_hash_map<int, int> parent_;

  std::variant<std::vector<KmerSet<K, KeyType>>, KmerSetSetNJ*> diffs_;

  // If diff_table_[i][j] is -1, it means that there is no addition from i to j.
  absl::flat_hash_map<int, absl::flat_hash_map<int, int>> diff_table_;
};

// KmerSetMM is similar to KmerSetSet, but it uses maximum weighted matching.
// If "approximate_matching" is true, it uses an approximate matching algorithm.
// If "approximate_weights" is true, it uses approximate weights for the
// matching algorithm. If "approximate_graph" is true, it uses a aparse graph
// for the matching algorithm.
template <int K, typename KeyType, typename CostFunctionType>
class KmerSetSetMM {
 public:
  KmerSetSetMM(std::vector<KmerSet<K, KeyType>> kmer_sets, int recursion_limit,
               bool approximate_matching, bool approximate_weights,
               bool approximate_graph, CostFunctionType cost_function,
               int n_workers) {
    const int n = kmer_sets.size();

    // If mates[i] is j (and mates[j] is i), we have a match between
    // kmer_sets[i] and kmer_sets[j].
    absl::flat_hash_map<int, int> mates;

    spdlog::debug("calculating mates");

    // The first two are for nodes, and the last one is for weights.
    std::vector<std::tuple<int, int, int64_t>> edges;

    // Calculates edges. The topology is affected by "approximate_graph" and the
    // weights are affected by "approximate_weights".
    {
      absl::InsecureBitGen bitgen;
      boost::asio::thread_pool pool(n_workers);
      std::mutex mu;

      for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
          if (approximate_graph && absl::Bernoulli(bitgen, 0.9)) continue;

          boost::asio::post(pool, [&, i, j] {
            int64_t weight =
                approximate_weights
                    ? kmer_sets[i].CommonEstimate(kmer_sets[j], 0.1)
                    : kmer_sets[i].Common(kmer_sets[j], 1);

            std::lock_guard lck(mu);
            edges.emplace_back(i, j, weight);
          });
        }
      }

      pool.join();
    }

    if (approximate_matching) {
      // Sorts "edges" so that the heaviest edge comes first.
      std::sort(edges.begin(), edges.end(),
                [](const std::tuple<int, int, int64_t>& lhs,
                   const std::tuple<int, int, int64_t>& rhs) {
                  return std::get<2>(lhs) > std::get<2>(rhs);
                });

      for (const std::tuple<int, int, int64_t>& edge : edges) {
        int i, j;
        std::tie(i, j, std::ignore) = edge;

        if (mates.find(i) == mates.end() && mates.find(j) == mates.end()) {
          mates[i] = j;
          mates[j] = i;
        }
      }
    } else {
      // Boost Graph Library was used at first here, but there seemed to be a
      // bug. Therefore, lemon library is used.

      lemon::ListGraph g;

      std::vector<lemon::ListGraph::Node> nodes(n);  // id to node.
      lemon::ListGraph::NodeMap<int> ids(g);         // node to id.

      // Adds nodes.
      for (int i = 0; i < n; i++) {
        lemon::ListGraph::Node node = g.addNode();
        ids[node] = i;
        nodes[i] = node;
      }

      lemon::ListGraph::EdgeMap<int64_t> weights(g);

      // Adds edges and sets their weights.
      for (const std::tuple<int, int, int64_t>& edge : edges) {
        int i, j;
        int64_t weight;
        std::tie(i, j, weight) = edge;

        lemon::ListGraph::Edge e = g.addEdge(nodes[i], nodes[j]);
        weights[e] = weight;
      }

      lemon::MaxWeightedMatching<lemon::ListGraph,
                                 lemon::ListGraph::EdgeMap<int64_t>>
          matching(g, weights);

      matching.run();

      for (int i = 0; i < n; i++) {
        lemon::ListGraph::Node mate = matching.mate(nodes[i]);

        if (mate != lemon::INVALID) {
          mates[i] = ids[mate];
        }
      }
    }

    spdlog::debug("calculated mates");

    std::vector<KmerSet<K, KeyType>> diffs;

    // Calculates parents, diffs, and diff_table_.
    // parents are calculated in the main thread, and the others are calculated
    // in sub threads.
    {
      int next_parent_id = n;

      std::mutex mu;
      boost::asio::thread_pool pool(n_workers);

      for (int i = 0; i < n; i++) {
        if (mates.find(i) == mates.end()) {
          parents_[i] = -1;

          boost::asio::post(pool, [&, i] {
            KmerSet<K, KeyType> kmer_set = kmer_sets[i];

            std::lock_guard lck(mu);
            diff_table_[-1][i] = diffs.size();
            diffs.push_back(std::move(kmer_set));
          });
        } else {
          int mate = mates[i];

          if (i < mate) {
            int parent = next_parent_id++;

            parents_[i] = parent;
            parents_[mate] = parent;

            // -1 is an empty set.
            parents_[parent] = -1;

            boost::asio::post(pool, [&, i, mate, parent] {
              KmerSet<K, KeyType> sub = Sub(kmer_sets[i], kmer_sets[mate], 1);

              std::lock_guard lck(mu);
              diff_table_[parent][i] = diffs.size();
              diffs.push_back(std::move(sub));
            });

            boost::asio::post(pool, [&, i, mate, parent] {
              KmerSet<K, KeyType> sub = Sub(kmer_sets[mate], kmer_sets[i], 1);

              std::lock_guard lck(mu);
              diff_table_[parent][mate] = diffs.size();
              diffs.push_back(std::move(sub));
            });

            boost::asio::post(pool, [&, i, mate, parent] {
              KmerSet<K, KeyType> intersection =
                  Intersection(kmer_sets[i], kmer_sets[mate], 1);

              std::lock_guard lck(mu);
              diff_table_[-1][parent] = diffs.size();
              diffs.push_back(std::move(intersection));
            });
          }
        }
      }

      pool.join();
    }

    // kmer_sets is no longer needed.
    std::vector<KmerSet<K, KeyType>>().swap(kmer_sets);

    // Calculates cost.
    {
      std::atomic_int64_t cost = 0;

      boost::asio::thread_pool pool(n_workers);

      for (size_t i = 0; i < diffs.size(); i++) {
        boost::asio::post(pool, [&, i] {
          cost += cost_function(KmerSet<K, KeyType>(), diffs[i], 1);
        });
      }

      pool.join();

      cost_ = cost;
    }

    spdlog::debug("cost_ = {}", cost_);

    if (recursion_limit == 0) {
      diffs_ = std::move(diffs);
    } else {
      diffs_ = new KmerSetSetMM(std::move(diffs), recursion_limit - 1,
                                approximate_matching, approximate_weights,
                                approximate_graph, std::move(cost_function),
                                n_workers);
    }
  }

  ~KmerSetSetMM() {
    if (!IsTerminal()) delete std::get<KmerSetSetMM*>(diffs_);
  }

  // Returns the number of k-mers stored internally.
  int64_t Size() const {
    if (IsTerminal()) {
      int64_t sum = 0;
      for (const KmerSet<K, KeyType>& kmer_set :
           std::get<std::vector<KmerSet<K, KeyType>>>(diffs_)) {
        sum += kmer_set.Size();
      }
      return sum;
    } else {
      return std::get<KmerSetSetMM*>(diffs_)->Size();
    }
  }

  int64_t Cost() const {
    if (IsTerminal()) return cost_;
    return std::get<KmerSetSetMM*>(diffs_)->Cost();
  }

  // Reconstructs the ith kmer set.
  KmerSet<K, KeyType> Get(int i, int n_workers) const {
    // Returns the diff from "from" to "to".
    const auto get_diff = [&](int from, int to) {
      int pos = diff_table_.find(from)->second.find(to)->second;

      if (IsTerminal()) {
        return std::get<std::vector<KmerSet<K, KeyType>>>(diffs_)[pos];
      } else {
        return std::get<KmerSetSetMM*>(diffs_)->Get(pos, n_workers);
      }
    };

    int parent = parents_.find(i)->second;

    if (parent == -1) {
      return get_diff(-1, i);
    }

    return Add(get_diff(-1, parent), get_diff(parent, i), n_workers);
  }

  // Dumps data to a vector of strings.
  // "canonical" will be used to compress kmer sets with KmerSetCompressed.
  std::vector<std::string> Dump(
      bool canonical, int n_workers,
      std::vector<std::string> v = std::vector<std::string>()) const {
    // Dumps cost_ to a line.
    {
      std::stringstream ss;
      ss << cost_;
      v.push_back(ss.str());
    }

    // Dumps parents_ to a line.
    // Format: "size key value key value ..."
    {
      std::stringstream ss;
      ss << parents_.size();

      for (const auto& p : parents_) {
        ss << ' ' << p.first << ' ' << p.second;
      }

      v.push_back(ss.str());
    }

    // Dumps diff_table_ to a line.
    // Format: "size key1 key2 value key1 key2 value ..."
    {
      std::vector<std::tuple<int, int, int>> flatten;
      for (const auto& p1 : diff_table_) {
        for (const auto& p2 : p1.second) {
          flatten.emplace_back(p1.first, p2.first, p2.second);
        }
      }

      std::stringstream ss;
      ss << flatten.size();

      for (const auto& t : flatten) {
        ss << ' ' << std::get<0>(t) << ' ' << std::get<1>(t) << ' '
           << std::get<2>(t);
      }

      v.push_back(ss.str());
    }

    if (IsTerminal()) {
      // Section delimiter.
      v.push_back("");

      const std::vector<KmerSet<K, KeyType>>& diffs =
          std::get<std::vector<KmerSet<K, KeyType>>>(diffs_);

      // Dumps diffs to 1 + diffs_.size() lines.
      // Each KmerSet is converted to KmerSetCompressed and dumped to one line.

      {
        std::stringstream ss;
        ss << diffs.size();
        v.push_back(ss.str());
      }

      // This will later be appended to "v".
      std::vector<std::string> buf(diffs.size());

      boost::asio::thread_pool pool(n_workers);

      for (size_t i = 0; i < diffs.size(); i++) {
        boost::asio::post(pool, [&, i] {
          const KmerSetCompressed<K, KeyType> kmer_set_compressed =
              KmerSetCompressed<K, KeyType>::FromKmerSet(diffs[i], canonical,
                                                         1);

          buf[i] = absl::StrJoin(kmer_set_compressed.Dump(), " ");
        });
      }

      pool.join();

      v.reserve(v.size() + buf.size());
      for (std::string& s : buf) v.push_back(std::move(s));

      return v;
    } else {
      // Recursive
      return std::get<KmerSetSetMM*>(diffs_)->Dump(canonical, n_workers,
                                                   std::move(v));
    }
  }

  // Dumps data to a file.
  void Dump(const std::string& file_name, const std::string& compressor,
            bool canonical, int n_workers) const {
    WriteLines(file_name, compressor, Dump(canonical, n_workers));
  }

  // Loads data from a vector of strings.
  static KmerSetSetMM* Load(const std::vector<std::string>& lines,
                            bool canonical, int n_workers, int64_t i = 0) {
    int64_t cost;
    absl::flat_hash_map<int, int> parents;
    absl::flat_hash_map<int, absl::flat_hash_map<int, int>> diff_table;

    {
      std::stringstream ss(lines[i]);
      ss >> cost;
    }

    {
      std::stringstream ss(lines[i + 1]);
      int64_t size;
      ss >> size;

      for (int64_t i = 0; i < size; i++) {
        int key, value;
        ss >> key >> value;
        parents[key] = value;
      }
    }

    {
      std::stringstream ss(lines[i + 2]);
      int64_t size;
      ss >> size;

      for (int64_t i = 0; i < size; i++) {
        int key1, key2, value;
        ss >> key1 >> key2 >> value;
        diff_table[key1][key2] = value;
      }
    }

    if (lines[i + 3] == "") {
      // The end of recursion.

      int64_t size;
      {
        std::stringstream ss(lines[i + 4]);
        ss >> size;
      }

      std::vector<KmerSet<K, KeyType>> diffs(size);

      boost::asio::thread_pool pool(n_workers);

      for (int64_t j = 0; j < size; j++) {
        boost::asio::post(pool, [&, j] {
          const std::string& line = lines[i + 5 + j];

          if (line.empty()) return;

          const KmerSetCompressed<K, KeyType> kmer_set_compressed =
              KmerSetCompressed<K, KeyType>::Load(absl::StrSplit(line, ' '));

          diffs[j] = kmer_set_compressed.ToKmerSet(canonical, 1);
        });
      }

      pool.join();

      return new KmerSetSetMM(cost, parents, diff_table, std::move(diffs));
    } else {
      return new KmerSetSetMM(cost, parents, diff_table,
                              Load(lines, canonical, n_workers, i + 3));
    }
  }

  static KmerSetSetMM* Load(const std::string& file_name,
                            const std::string& decompressor, bool canonical,
                            int n_workers) {
    return Load(ReadLines(file_name, decompressor), canonical, n_workers);
  }

 private:
  KmerSetSetMM(
      int64_t cost, absl::flat_hash_map<int, int> parents,
      absl::flat_hash_map<int, absl::flat_hash_map<int, int>> diff_table,
      std::vector<KmerSet<K, KeyType>> diffs)
      : cost_(cost),
        parents_(std::move(parents)),
        diff_table_(std::move(diff_table)),
        diffs_(std::move(diffs)) {}

  KmerSetSetMM(
      int64_t cost, absl::flat_hash_map<int, int> parents,
      absl::flat_hash_map<int, absl::flat_hash_map<int, int>> diff_table,
      KmerSetSetMM* diffs)
      : cost_(cost),
        parents_(std::move(parents)),
        diff_table_(std::move(diff_table)),
        diffs_(diffs) {}

  // Returns true if it is the end of recursion.
  bool IsTerminal() const {
    return std::holds_alternative<std::vector<KmerSet<K, KeyType>>>(diffs_);
  }

  int64_t cost_ = 0;

  // If parents_[i] == -1, the parent of i is an empty set.
  absl::flat_hash_map<int, int> parents_;

  // If diffs_ is a vector, diffs_[diff_table[i][j]] is a diff from i to j.
  // If diffs_ is a KmerSetSetMM*, diffs_.Get(diff_table[i][j]) is a diff from
  // i to j.
  absl::flat_hash_map<int, absl::flat_hash_map<int, int>> diff_table_;
  std::variant<std::vector<KmerSet<K, KeyType>>, KmerSetSetMM*> diffs_;
};

#endif

#ifndef KMER_SET_SET_H_
#define KMER_SET_SET_H_

#include <algorithm>
#include <cmath>
#include <mutex>
#include <optional>
#include <thread>
#include <tuple>
#include <variant>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "graph.h"
#include "kmer.h"
#include "kmer_set.h"
#include "kmer_set_compact.h"
#include "lemon/list_graph.h"
#include "lemon/matching.h"
#include "neighbor_joining.h"
#include "range.h"
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
        std::vector<KmerSet<K, KeyType>> diffs =
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

    std::vector<KmerSet<K, KeyType>> diffs =
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
        std::vector<KmerSet<K, KeyType>> diffs =
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
template <int K, typename KeyType, typename CostFunctionType>
class KmerSetSetMM {
 public:
  KmerSetSetMM(std::vector<KmerSet<K, KeyType>> kmer_sets, int recursion_limit,
               CostFunctionType cost_function, int n_workers)
      : n_(kmer_sets.size()) {
    using Graph = lemon::ListGraph;

    Graph g;

    std::vector<Graph::Node> nodes(n_);
    Graph::NodeMap<int> ids(g);
    for (int i = 0; i < n_; i++) {
      Graph::Node node = g.addNode();
      ids[node] = i;
      nodes[i] = node;
    }

    Graph::EdgeMap<int64_t> weights(g);

    // Adds edges and sets their weights.
    {
      boost::asio::thread_pool pool(n_workers);
      std::mutex mu;

      for (int i = 0; i < n_; i++) {
        for (int j = i + 1; j < n_; j++) {
          boost::asio::post(pool, [&, i, j] {
            int64_t weight = kmer_sets[i].Common(kmer_sets[j], 1);

            std::lock_guard lck(mu);
            Graph::Edge edge = g.addEdge(nodes[i], nodes[j]);
            weights[edge] = weight;
          });
        }
      }

      pool.join();
    }

    lemon::MaxWeightedMatching<Graph, Graph::EdgeMap<int64_t>> matching(
        g, weights);

    matching.run();

    std::vector<KmerSet<K, KeyType>> diffs;

    // Calculates parents, diffs, and diff_table_.
    // parents are calculated in the main thread, and the others are calculated
    // in sub threads.
    {
      int next_parent_id = n_;

      std::mutex mu;
      boost::asio::thread_pool pool(n_workers);

      for (int i = 0; i < n_; i++) {
        lemon::ListGraph::Node mate = matching.mate(nodes[i]);

        if (mate == lemon::INVALID) {
          parents_[i] = -1;

          boost::asio::post(pool, [&, i] {
            KmerSet<K, KeyType> kmer_set = kmer_sets[i];

            std::lock_guard lck(mu);
            diff_table_[-1][i] = diffs.size();
            diffs.push_back(std::move(kmer_set));
          });
        } else {
          int mate_i = ids[mate];

          if (i < mate_i) {
            int parent = next_parent_id++;

            parents_[i] = parent;
            parents_[mate_i] = parent;

            // -1 is an empty set.
            parents_[parent] = -1;

            boost::asio::post(pool, [&, i, mate_i, parent] {
              KmerSet<K, KeyType> sub =
                  Sub(kmer_sets[i], kmer_sets[mate_i], 1);

              std::lock_guard lck(mu);
              diff_table_[parent][i] = diffs.size();
              diffs.push_back(std::move(sub));
            });

            boost::asio::post(pool, [&, i, mate_i, parent] {
              KmerSet<K, KeyType> sub =
                  Sub(kmer_sets[mate_i], kmer_sets[i], 1);

              std::lock_guard lck(mu);
              diff_table_[parent][mate_i] = diffs.size();
              diffs.push_back(std::move(sub));
            });

            boost::asio::post(pool, [&, i, mate_i, parent] {
              KmerSet<K, KeyType> intersection =
                  Intersection(kmer_sets[i], kmer_sets[mate_i], 1);

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
                                cost_function, n_workers);
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

 private:
  // Returns true if it is the end of recursion.
  bool IsTerminal() const {
    return std::holds_alternative<std::vector<KmerSet<K, KeyType>>>(diffs_);
  }

  int n_;

  int64_t cost_ = 0;

  // If parents_[i] == -1, the parent of i is an empty set.
  absl::flat_hash_map<int, int> parents_;

  // If diffs_ is a vector, diffs_[diff_table[i][j]] is a diff set from i to j.
  absl::flat_hash_map<int, absl::flat_hash_map<int, int>> diff_table_;

  std::variant<std::vector<KmerSet<K, KeyType>>, KmerSetSetMM*> diffs_;
};

#endif
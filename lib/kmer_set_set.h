#ifndef KMER_SET_SET_H_
#define KMER_SET_SET_H_

#include <algorithm>
#include <optional>
#include <thread>
#include <variant>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "graph.h"
#include "kmer.h"
#include "kmer_set.h"
#include "kmer_set_compact.h"
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

    bool use_recursion;

    if (recursion_limit == 0) {
      use_recursion = false;
    } else {
      bool improve_by_recursion;

      {
        diffs.push_back(KmerSet<K, KeyType>());

        // The same MST can be constructed twice.
        // This is for the sake of simplicity of implementation.
        improve_by_recursion =
            ConstructMST(diffs, diffs.size() - 1, cost_function, n_workers)
                .first < cost_;

        diffs.pop_back();
      }

      spdlog::debug("improve_by_recursion = {}", improve_by_recursion);

      use_recursion = improve_by_recursion;
    }

    if (use_recursion) {
      diffs_ = new KmerSetSet(std::move(diffs), recursion_limit - 1,
                              cost_function, n_workers);
    } else {
      diffs_ = diffs;
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

          for (std::thread& thread : threads) thread.join();
        }

        kmer_set.Add(add, n_workers).Sub(sub, n_workers);
      }
    }

    return kmer_set;
  }

  // Returns the number of k-mers stored internally.
  int64_t Size() const {
    int64_t size = 0;

    if (IsTerminal()) {
      std::vector<KmerSet<K, KeyType>> diffs =
          std::get<std::vector<KmerSet<K, KeyType>>>(diffs_);

      for (const auto& diff : diffs) {
        size += diff.Size();
      }
    } else {
      KmerSetSet* diffs = std::get<KmerSetSet*>(diffs_);

      return diffs->Size();
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

  // Returns the level of recursion.
  // If it is not a recursive structure, 0 is returned.
  int Depth() const {
    int depth = 0;

    const KmerSetSet* current = this;

    while (!current->IsTerminal()) {
      current = std::get<KmerSetSet*>(current->diffs_);
      depth += 1;
    }

    return depth;
  }

 private:
  bool IsTerminal() const {
    return std::holds_alternative<std::vector<KmerSet<K, KeyType>>>(diffs_);
  }

  static std::pair<int64_t, Tree> ConstructMST(
      const std::vector<KmerSet<K, KeyType>>& kmer_sets, int root,
      CostFunctionType cost_function, int n_workers) {
    std::vector<std::pair<int, int>> pairs;

    for (size_t i = 0; i < kmer_sets.size(); i++) {
      for (size_t j = i + 1; j < kmer_sets.size(); j++) {
        pairs.emplace_back(i, j);
      }
    }

    BidirectionalGraph g;
    std::vector<std::thread> threads;
    std::mutex mu;

    for (const Range& range : Range(0, pairs.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        for (int i = range.begin; i < range.end; i++) {
          std::pair<int, int> pair = pairs[i];

          int64_t cost =
              cost_function(kmer_sets[pair.first], kmer_sets[pair.second], 1);

          std::lock_guard lck(mu);

          g.AddEdge(pair.first, pair.second, cost);
        }
      });
    }

    for (std::thread& thread : threads) thread.join();

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

#endif
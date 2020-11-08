#ifndef KMER_SET_SET_H_
#define KMER_SET_SET_H_

#include <algorithm>
#include <optional>
#include <variant>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "graph.h"
#include "kmer.h"
#include "kmer_set.h"
#include "kmer_set_compact.h"
#include "spdlog/spdlog.h"

template <int K, int B>
std::pair<int64_t, Tree> ConstructMST(
    const std::vector<KmerSet<K, B>>& kmer_sets, int root) {
  BidirectionalGraph g;

  for (size_t i = 0; i < kmer_sets.size(); i++) {
    for (size_t j = i + 1; j < kmer_sets.size(); j++) {
      g.AddEdge(i, j,
                (kmer_sets[i] - kmer_sets[j]).Size() +
                    (kmer_sets[j] - kmer_sets[i]).Size());
    }
  }

  return g.MST(root);
}

// KmerSetSet can be used to represent multiple k-mer sets in less space.
// It is a recursive structure. That is, a KmerSetSet can contain another
// KmerSetSet internally.
template <int K, int B>
class KmerSetSet {
 public:
  KmerSetSet(std::vector<KmerSet<K, B>> kmer_sets) : n_(kmer_sets.size()) {
    // kmer_sets[n_] is an empty set.
    kmer_sets.push_back(KmerSet<K, B>());

    int64_t cost;
    std::tie(cost, tree_) = ConstructMST(kmer_sets, n_);

    spdlog::debug("constructed MST: cost = {}", cost);

    std::vector<KmerSet<K, B>> diffs;

    for (int i = 0; i < n_; i++) {
      const int parent = tree_->Parent(i);

      diff_table_[i][parent] = diffs.size();
      diffs.push_back(kmer_sets[parent] - kmer_sets[i]);

      diff_table_[parent][i] = diffs.size();
      diffs.push_back(kmer_sets[i] - kmer_sets[parent]);
    }

    bool improve_by_recursion;

    {
      diffs.push_back(KmerSet<K, B>());

      // The same MST can be constructed twice.
      // This is for the sake of simplicity of implementation.
      improve_by_recursion = ConstructMST(diffs, diffs.size() - 1).first < cost;

      diffs.pop_back();
    }

    spdlog::debug("improve_by_recursion = {}", improve_by_recursion);

    if (improve_by_recursion) {
      diffs_ = new KmerSetSet(diffs);
    } else {
      diffs_ = diffs;
    }
  }

  ~KmerSetSet() {
    KmerSetSet* current = this;

    std::vector<KmerSetSet*> ptrs;

    while (!current->IsTerminal()) {
      KmerSetSet* next = std::get<KmerSetSet*>(current->diffs_);
      ptrs.push_back(next);
      current = next;
    }

    for (KmerSetSet* ptr : ptrs) delete ptr;
  }

  // Returns the ith k-mer set.
  KmerSet<K, B> Get(int i) const {
    std::vector<int> path;

    while (true) {
      path.push_back(i);
      if (i == n_) break;
      i = tree_->Parent(i);
    }

    std::reverse(path.begin(), path.end());

    KmerSet<K, B> kmer_set;

    for (size_t i = 0; i < path.size() - 1; i++) {
      if (IsTerminal()) {
        std::vector<KmerSet<K, B>> diffs =
            std::get<std::vector<KmerSet<K, B>>>(diffs_);

        kmer_set =
            kmer_set +
            diffs[diff_table_.find(path[i])->second.find(path[i + 1])->second] -
            diffs[diff_table_.find(path[i + 1])->second.find(path[i])->second];
      } else {
        KmerSetSet<K, B>* diffs = std::get<KmerSetSet<K, B>*>(diffs_);

        kmer_set =
            kmer_set +
            diffs->Get(
                diff_table_.find(path[i])->second.find(path[i + 1])->second) -
            diffs->Get(
                diff_table_.find(path[i + 1])->second.find(path[i])->second);
      }
    }

    return kmer_set;
  }

  int64_t Size() const {
    int64_t size = 0;

    if (IsTerminal()) {
      std::vector<KmerSet<K, B>> diffs =
          std::get<std::vector<KmerSet<K, B>>>(diffs_);

      for (const auto& diff : diffs) {
        size += diff.Size();
      }
    } else {
      KmerSetSet<K, B>* diffs = std::get<KmerSetSet<K, B>*>(diffs_);

      return diffs->Size();
    }

    return size;
  }

  int Depth() const {
    int depth = 0;

    const KmerSetSet* current = this;

    while (!current->IsTerminal()) {
      current = std::get<KmerSetSet<K, B>*>(current->diffs_);
      depth += 1;
    }

    return depth;
  }

 private:
  bool IsTerminal() const {
    return std::holds_alternative<std::vector<KmerSet<K, B>>>(diffs_);
  }

  int n_;

  // Tree does not have a default constructor.
  std::optional<Tree> tree_;

  // If diffs_ is of type std::vector<KmerSet<K, B>>,
  // diffs_[diff_table_[i][j]] represents the diff from ith data to jth data.
  // If diffs_ is of type KmerSetSet*, diffs_->Get(diff_table_[i][j]) represents
  // the diff from ith data to jth data.
  std::variant<std::vector<KmerSet<K, B>>, KmerSetSet*> diffs_;
  absl::flat_hash_map<int, absl::flat_hash_map<int, int>> diff_table_;
};

#endif
#ifndef KMER_SET_SET_H_
#define KMER_SET_SET_H_

#include <algorithm>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "graph.h"
#include "kmer.h"
#include "kmer_set.h"
#include "kmer_set_compact.h"

// KmerSetSet can be used to represent multiple k-mer sets in less space.
template <int K, int B>
class KmerSetSet {
 public:
  KmerSetSet(std::vector<KmerSet<K, B>> kmer_sets) : n_(kmer_sets.size()) {
    // kmer_sets[n_] is an empty set.
    kmer_sets.push_back(KmerSet<K, B>());

    BidirectionalGraph g;

    for (size_t i = 0; i < kmer_sets.size(); i++) {
      for (size_t j = i + 1; j < kmer_sets.size(); j++) {
        g.AddEdge(i, j,
                  (kmer_sets[i] - kmer_sets[j]).Size() +
                      (kmer_sets[j] - kmer_sets[i]).Size());
      }
    }

    tree_ = g.MST(kmer_sets.size() - 1).second;

    for (int i = 0; i < n_; i++) {
      const int parent = tree_->Parent(i);
      diff_[i][parent] = kmer_sets[parent] - kmer_sets[i];
      diff_[parent][i] = kmer_sets[i] - kmer_sets[parent];
    }
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
      kmer_set = kmer_set +
                 diff_.find(path[i])->second.find(path[i + 1])->second -
                 diff_.find(path[i + 1])->second.find(path[i])->second;
    }

    return kmer_set;
  }

  int64_t Size() const {
    int64_t size = 0;

    for (const auto& p1 : diff_) {
      for (const auto& p2 : p1.second) {
        size += p2.second.Size();
      }
    }

    return size;
  }

  std::vector<KmerSet<K, B>> Diffs() const {
    std::vector<KmerSet<K, B>> diffs;

    for (const auto& p1 : diff_) {
      for (const auto& p2 : p1.second) {
        diffs.push_back(p2.second);
      }
    }

    return diffs;
  }

 private:
  int n_;

  // Tree does not have a default constructor.
  std::optional<Tree> tree_;

  absl::flat_hash_map<int, absl::flat_hash_map<int, KmerSet<K, B>>> diff_;
};

#endif
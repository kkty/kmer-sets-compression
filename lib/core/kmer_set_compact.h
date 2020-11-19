#ifndef CORE_KMER_SET_COMPACT_H_
#define CORE_KMER_SET_COMPACT_H_

#include <algorithm>
#include <atomic>
#include <cassert>
#include <mutex>
#include <thread>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "core/kmer_set.h"
#include "core/range.h"
#include "core/suffix_array.h"
#include "core/unitigs.h"

// KmerSetCompact holds a set of k-mers in less space.
template <int K>
class KmerSetCompact {
 public:
  template <typename KeyType>
  KmerSetCompact(const KmerSet<K, KeyType>& kmer_set, bool canonical,
                 int n_workers)
      : canonical_(canonical) {
    const std::vector<std::string> unitigs =
        canonical ? GetUnitigsCanonical(kmer_set, n_workers)
                  : GetUnitigs(kmer_set, n_workers);

    std::string concatenated;

    // Pre-calculate the required length.
    {
      int64_t sum = 0;
      for (const std::string& unitig : unitigs) sum += unitig.length();
      concatenated.reserve(sum);
    }

    for (const std::string& unitig : unitigs) {
      boundary_.push_back(concatenated.length());
      concatenated += unitig;
    }

    sa_ = SuffixArray(concatenated, n_workers);
  }

  bool Contains(const Kmer<K>& kmer) const {
    const auto contains = [&](const Kmer<K>& kmer) {
      std::vector<int64_t> v = sa_.Find(kmer.String());

      // Returns the largest j such that boundary_[j] <= i;
      const auto f = [&](int64_t i) {
        // begin == -1 or boundary_[begin] <= i
        int64_t begin = -1;

        // end == boundary_.size() or boundary_[end] > i
        int64_t end = boundary_.size();

        while (end - begin > 1) {
          int64_t mid = (begin + end) / 2;

          if (boundary_[mid] <= i)
            begin = mid;
          else
            end = mid;
        }

        return begin;
      };

      return std::any_of(v.begin(), v.end(),
                         [&](int64_t i) { return f(i) == f(i + K - 1); });
    };

    if (canonical_) {
      return contains(kmer) && contains(kmer.Complement());
    } else {
      return contains(kmer);
    }
  }

  // Returns the k-mers that match the condition.
  template <typename Pred>
  std::vector<Kmer<K>> Find(Pred&& pred, int n_workers) const {
    std::vector<Kmer<K>> kmers;
    std::mutex mu;

    std::vector<std::thread> threads;

    for (const Range& range : Range(0, boundary_.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<Kmer<K>> buf;

        range.ForEach([&](int64_t i) {
          const int64_t begin = boundary_[i];
          const int64_t end = (i == (int64_t)boundary_.size() - 1)
                                  ? sa_.String().length()
                                  : boundary_[i + 1];
          for (int64_t j = 0; begin + j + K - 1 < end; j++) {
            const Kmer<K> kmer{sa_.String().substr(begin + j, K)};
            if (pred(kmer)) buf.push_back(canonical_ ? kmer.Canonical() : kmer);
          }
        });

        std::lock_guard lck(mu);
        for (const Kmer<K>& kmer : buf) kmers.push_back(kmer);
      });
    }

    for (std::thread& thread : threads) thread.join();

    return kmers;
  }

  // Return all the k-mers.
  std::vector<Kmer<K>> Find(int n_workers) const {
    return Find([&](const Kmer<K>&) { return true; }, n_workers);
  }

  int64_t Size() const { return sa_.String().length(); }

 private:
  bool canonical_;

  SuffixArray sa_;

  std::vector<int64_t> boundary_;
};

#endif
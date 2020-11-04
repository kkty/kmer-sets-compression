#ifndef KMER_SET_COMPACT_H_
#define KMER_SET_COMPACT_H_

#include <algorithm>
#include <atomic>
#include <mutex>
#include <thread>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "kmer_set.h"
#include "range.h"
#include "suffix_array.h"

template <int K>
std::string GetUnitigFromKmers(const std::vector<Kmer<K>>& kmers) {
  std::string s;
  s.reserve(K + kmers.size() - 1);
  s += kmers[0].String();
  for (int64_t i = 1; i < (int64_t)kmers.size(); i++) {
    s += kmers[i].String()[K - 1];
  }
  return s;
}

template <int K>
class KmerSetCompact {
 public:
  template <int B>
  KmerSetCompact(const KmerSet<K, B>& kmer_set) {
    const auto start_kmers = kmer_set.Find([&](const Kmer<K>& kmer) {
      const auto prevs = kmer.Prevs();

      // If the k-mer has no incoming edges.
      if (std::none_of(prevs.begin(), prevs.end(), [&](const Kmer<K>& prev) {
            return kmer_set.Contains(prev) && prev != kmer;
          }))
        return true;

      // If the k-mer has multiple incoming edges.
      if (std::count_if(prevs.begin(), prevs.end(), [&](const Kmer<K>& prev) {
            return kmer_set.Contains(prev) && prev != kmer;
          }) >= 2)
        return true;

      // There is an edge from "prev" to "kmer".
      const auto prev =
          *std::find_if(prevs.begin(), prevs.end(), [&](const Kmer<K>& prev) {
            return kmer_set.Contains(prev) && prev != kmer;
          });

      const auto prev_nexts = prev.Nexts();

      // If "prev" has multiple outgoing edges.
      if (std::count_if(prev_nexts.begin(), prev_nexts.end(),
                        [&](const Kmer<K>& prev_next) {
                          return kmer_set.Contains(prev_next) &&
                                 prev_next != prev;
                        }) >= 2)
        return true;

      return false;
    });

    const auto end_kmers = [&] {
      const auto v = kmer_set.Find([&](const Kmer<K>& kmer) {
        const auto nexts = kmer.Nexts();

        // If the k-mer has no outgoing edges.
        if (std::none_of(nexts.begin(), nexts.end(), [&](const Kmer<K>& next) {
              return kmer_set.Contains(next) && next != kmer;
            }))
          return true;

        // If the k-mer has multiple outgoing edges.
        if (std::count_if(nexts.begin(), nexts.end(), [&](const Kmer<K>& next) {
              return kmer_set.Contains(next) && next != kmer;
            }) >= 2)
          return true;

        // There is an edge from "kmer" to "next".
        const auto next =
            *std::find_if(nexts.begin(), nexts.end(), [&](const Kmer<K>& next) {
              return kmer_set.Contains(next) && next != kmer;
            });

        const auto next_prevs = next.Prevs();

        // If "next" has multiple incoming edges.
        if (std::count_if(next_prevs.begin(), next_prevs.end(),
                          [&](const Kmer<K>& next_prev) {
                            return kmer_set.Contains(next_prev) &&
                                   next_prev != next;
                          }) >= 2)
          return true;

        return false;
      });

      const absl::flat_hash_set<Kmer<K>> s(v.begin(), v.end());

      return s;
    }();

    std::vector<std::string> unitigs;
    KmerSet<K, B> visited;

    std::mutex mu_unitigs;
    std::mutex mu_visited;

    std::vector<std::thread> threads;

    for (const auto& range : SplitRange(0, start_kmers.size(),
                                        std::thread::hardware_concurrency())) {
      threads.emplace_back([&, range] {
        std::vector<std::string> buf_unitigs;
        KmerSet<K, B> buf_visited;

        for (int64_t i = range.first; i < range.second; i++) {
          const Kmer<K>& start_kmer = start_kmers[i];
          std::vector<Kmer<K>> path;

          Kmer<K> current = start_kmer;
          while (true) {
            buf_visited.Add(current);
            path.push_back(current);
            if (end_kmers.find(current) != end_kmers.end()) break;
            const auto nexts = current.Nexts();
            current = *std::find_if(
                nexts.begin(), nexts.end(), [&](const Kmer<K>& next) {
                  return kmer_set.Contains(next) && next != current;
                });
          }

          buf_unitigs.push_back(GetUnitigFromKmers(path));
        }

        {
          std::lock_guard _{mu_visited};
          visited += buf_visited;
        }

        {
          std::lock_guard _{mu_unitigs};
          unitigs.reserve(unitigs.size() + buf_unitigs.size());
          for (const auto& unitig : buf_unitigs) unitigs.push_back(unitig);
        }
      });
    }

    for (std::thread& thread : threads) thread.join();

    // Consider loops where every no has 1 incoming edge and 1 outgoing edge.
    {
      const std::vector<Kmer<K>> not_visited = kmer_set.Find(
          [&](const Kmer<K>& kmer) { return !visited.Contains(kmer); });

      for (const Kmer<K>& kmer : not_visited) {
        if (visited.Contains(kmer)) continue;

        Kmer<K> current = kmer;

        std::vector<Kmer<K>> path;

        while (!visited.Contains(current)) {
          path.push_back(current);
          visited.Add(current);
          const auto nexts = current.Nexts();
          current = *std::find_if(
              nexts.begin(), nexts.end(),
              [&](const Kmer<K>& next) { return kmer_set.Contains(next); });
        }

        unitigs.push_back(GetUnitigFromKmers(path));
      }
    }

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

    sa_ = SuffixArray(concatenated);
  }

  bool Contains(const Kmer<K>& kmer) {
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
  }

  // Returns the k-mers that match the condition.
  template <typename Pred>
  std::vector<Kmer<K>> Find(Pred&& pred) const {
    std::vector<Kmer<K>> kmers;
    std::mutex mu;

    std::vector<std::thread> threads;

    for (const auto& range :
         SplitRange(0, boundary_.size(), std::thread::hardware_concurrency())) {
      threads.emplace_back([&] {
        std::vector<Kmer<K>> buf;

        for (int64_t i = range.first; i < range.second; i++) {
          const int64_t begin = boundary_[i];
          const int64_t end = (i == (int64_t)boundary_.size() - 1)
                                  ? sa_.String().length()
                                  : boundary_[i + 1];
          for (int64_t j = 0; begin + j + K - 1 < end; j++) {
            const Kmer<K> kmer{sa_.String().substr(begin + j, K)};
            if (pred(kmer)) buf.push_back(kmer);
          }
        }

        std::lock_guard _{mu};
        kmers.reserve(kmers.size() + buf.size());
        for (const Kmer<K>& kmer : buf) kmers.push_back(kmer);
      });
    }

    for (std::thread& thread : threads) thread.join();

    return kmers;
  }

  // Return all the k-mers.
  std::vector<Kmer<K>> Find() const {
    return Find([&](const Kmer<K>&) { return true; });
  }

  int64_t Size() const { return sa_.String().length(); }

 private:
  SuffixArray sa_;
  std::vector<int64_t> boundary_;
};

#endif
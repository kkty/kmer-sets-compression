#ifndef CORE_KMER_SET_COMPACT_H_
#define CORE_KMER_SET_COMPACT_H_

#include <algorithm>
#include <atomic>
#include <cassert>
#include <mutex>
#include <thread>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "core/kmer.h"
#include "core/kmer_set.h"
#include "core/range.h"
#include "core/unitigs.h"

// Represents a set of kmers with SPSS, thereby reducing space.
template <int K, typename KeyType>
class KmerSetCompact {
 public:
  KmerSetCompact() = default;

  // Constructs KmerSetCompact from a kmer set. If the kmer set is for
  // storing canonical kmers, "canonical" should be true.
  static KmerSetCompact FromKmerSet(const KmerSet<K, KeyType>& kmer_set,
                                    bool canonical, int n_workers) {
    std::vector<std::string> spss;

    if (canonical) {
      spss = GetSPSSCanonical<K, KeyType>(kmer_set, n_workers);
    } else {
      spss = GetUnitigs<K, KeyType>(kmer_set, n_workers);
    }

    return KmerSetCompact(std::move(spss));
  }

  // Constructs a KmerSet.
  // If we are considering canonical kmers, "canonical" should be true.
  KmerSet<K, KeyType> ToKmerSet(bool canonical, int n_workers) const {
    std::vector<Kmer<K>> kmers;
    std::vector<std::thread> threads;
    std::mutex mu;

    for (const Range& range : Range(0, spss_.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<Kmer<K>> buf;

        range.ForEach([&](int i) {
          const std::string s = spss_[i];
          for (size_t i = 0; i < s.size() + 1 - K; i++) {
            const Kmer<K> kmer(s.substr(i, K));
            buf.push_back(canonical ? kmer.Canonical() : kmer);
          }
        });

        std::lock_guard lck(mu);
        for (Kmer<K>& kmer : buf) kmers.push_back(std::move(kmer));
      });
    }

    for (std::thread& t : threads) t.join();

    return KmerSet<K, KeyType>(kmers, n_workers);
  }

  // Dumps data to a file.
  // Dumped data can be read by Load().
  absl::Status Dump(const std::string& file_name,
                    const std::string& compressor) const {
    return WriteLines(file_name, compressor, spss_);
  }

  // Dumps data to a vector of strings.
  // Dumped data can be read by Load().
  // Each string does not contain whitespaces or line breaks.
  std::vector<std::string> Dump() const { return spss_; }

  // Loads data from a file.
  static absl::StatusOr<KmerSetCompact> Load(const std::string& file_name,
                                             const std::string& decompressor) {
    std::vector<std::string> lines;

    {
      absl::StatusOr<std::vector<std::string>> statusor = ReadLines(file_name, decompressor);

      if (!statusor.ok()) {
        return statusor.status();
      }

      lines = std::move(statusor).value();
    }

    return KmerSetCompact(std::move(lines));
  }

  // Loads data from a vector of strings.
  static KmerSetCompact Load(std::vector<std::string> v) {
    return KmerSetCompact(std::move(v));
  }

  // Returns the total length of stored strings.
  int64_t Size(int n_workers) const {
    std::vector<std::thread> threads;
    std::atomic_int64_t size = 0;

    for (const Range& range : Range(0, spss_.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        range.ForEach([&](int i) { size += spss_[i].length(); });
      });
    }

    for (std::thread& t : threads) t.join();

    return size;
  }

 private:
  KmerSetCompact(std::vector<std::string> spss) : spss_(std::move(spss)){};

  std::vector<std::string> spss_;
};

#endif
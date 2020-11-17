#ifndef KMER_SET_COMPRESSED_H_
#define KMER_SET_COMPRESSED_H_

#include <algorithm>
#include <atomic>
#include <cassert>
#include <mutex>
#include <thread>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "kmer_set.h"
#include "range.h"
#include "suffix_array.h"
#include "unitigs.h"

// Represents a set of kmers with unitigs, thereby reducing space.
template <int K, typename KeyType>
class KmerSetCompressed {
 public:
  static KmerSetCompressed FromKmerSet(const KmerSet<K, KeyType>& kmer_set,
                                       bool canonical, int n_workers) {
    std::vector<std::string> unitigs;

    if (canonical) {
      unitigs = GetUnitigsCanonical<K, KeyType>(kmer_set, n_workers);
    } else {
      unitigs = GetUnitigs<K, KeyType>(kmer_set, n_workers);
    }

    return KmerSetCompressed(std::move(unitigs), n_workers);
  }

  KmerSet<K, KeyType> ToKmerSet(int n_workers) const {
    std::vector<std::string> kmers;
    std::vector<std::thread> threads;
    std::mutex mu;

    for (const Range& range : Range(0, unitigs_.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<std::string> buf;

        range.ForEach([&](int i) {
          const std::string s = unitigs_[i];
          for (size_t i = 0; i < s.size() + 1 - K; i++) {
            buf.push_back(s.substr(i, K));
          }
        });

        std::lock_guard lck(mu);
        for (std::string& s : buf) kmers.push_back(std::move(s));
      });
    }

    for (std::thread& t : threads) t.join();

    return KmerSet<K, KeyType>(kmers, n_workers);
  }

  // Dumps unitigs to a file.
  void Dump(const std::string& file_name, const std::string& compressor) const {
    WriteLines(file_name, compressor, unitigs_);
  }

  // Loads unitigs from a file.
  static KmerSetCompressed Load(const std::string& file_name,
                                const std::string& decompressor,
                                int n_workers) {
    return KmerSetCompressed(ReadLines(file_name, decompressor), n_workers);
  }

  int64_t Size() const { return size_; }

 private:
  KmerSetCompressed(std::vector<std::string> unitigs, int n_workers)
      : unitigs_(std::move(unitigs)) {
    std::vector<std::thread> threads;
    std::atomic_int64_t size = 0;

    for (const Range& range : Range(0, unitigs_.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        range.ForEach([&](int i) { size += unitigs_[i].length(); });
      });
    }

    for (std::thread& t : threads) t.join();

    size_ = size;
  };

  int64_t size_ = 0;

  std::vector<std::string> unitigs_;
};

#endif
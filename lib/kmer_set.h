#ifndef KMER_SET_H_
#define KMER_SET_H_

#include <array>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/strings/str_split.h"
#include "kmer.h"
#include "range.h"

// The first N bits are used to get the bucket ID.
template <int K, int N>
std::pair<int, std::bitset<K * 2 - N>> GetBucketAndKeyFromKmer(
    const Kmer<K>& kmer) {
  const auto& bits = kmer.Bits();
  const int key_bits = K * 2 - N;
  const int bucketID = bits.to_ullong() >> key_bits;
  const std::bitset<K * 2 - N> key{bits.to_ullong() % (1ull << key_bits)};
  return {bucketID, key};
}

template <int K, int N>
Kmer<K> GetKmerFromBucketAndKey(int bucketID,
                                const std::bitset<K * 2 - N>& key) {
  const int key_bits = K * 2 - N;
  std::bitset<K * 2> bits{((uint64_t)bucketID << key_bits) + key.to_ullong()};
  return Kmer<K>{bits};
}

// Use (1 << B) buckets internally.
template <int K, int B>
class KmerSet {
 public:
  // Loads from dumped data.
  // Returns "true" if successful.
  bool Load(std::istream& is) {
    int k, b;
    is >> k >> b;

    if (k != K || b != B) {
      return false;
    }

    for (int i = 0; i < (1 << B); i++) {
      int64_t size;
      is >> size;
      buckets_[i].reserve(size);

      for (int64_t j = 0; j < size; j++) {
        std::string kmer_s;
        is >> kmer_s;
        std::bitset<K * 2 - B> kmer_b(kmer_s);
        buckets_[i].insert(std::move(kmer_b));
      }
    }

    return true;
  }

  void Dump(std::ostream& os) const {
    os << K << ' ' << B << '\n';

    for (int i = 0; i < (1 << B); i++) {
      os << buckets_[i].size() << '\n';
      for (const auto& key : buckets_[i]) {
        os << key.to_string() << '\n';
      }
    }
  }

  int Size() const {
    int sum = 0;
    for (const auto& bucket : buckets_) {
      sum += bucket.size();
    }
    return sum;
  }

  void Add(const Kmer<K>& kmer) {
    const auto [bucket, key] = GetBucketAndKeyFromKmer<K, B>(kmer);
    buckets_[bucket].insert(key);
  }

  void Add(int bucket, const std::bitset<K * 2 - B>& key) {
    buckets_[bucket].insert(key);
  }

  bool Contains(const Kmer<K>& kmer) const {
    const auto [bucket, key] = GetBucketAndKeyFromKmer<K, B>(kmer);
    return buckets_[bucket].find(key) != buckets_[bucket].end();
  }

  // Fetches one k-mer randomly.
  Kmer<K> Sample() const {
    const int bucket = rand() % (1 << B);
    const int64_t n = rand() % buckets_[bucket].size();

    auto it = buckets_[bucket].begin();
    for (int64_t i = 0; i < n; i++) it++;

    const auto key = *it;

    return GetKmerFromBucketAndKey<K, B>(bucket, key);
  }

  template <typename Pred>
  std::vector<Kmer<K>> Find(Pred&& pred) {
    std::vector<Kmer<K>> kmers;
    std::mutex mu;

    ForEachBucket([&](const Bucket& bucket, int bucketID) {
      std::vector<Kmer<K>> buf;

      for (const auto& key : bucket) {
        const Kmer<K> kmer = GetKmerFromBucketAndKey<K, B>(bucketID, key);
        if (pred(kmer)) buf.push_back(kmer);
      }

      std::lock_guard _{mu};
      kmers.reserve(kmers.size() + buf.size());
      for (const auto& kmer : buf) {
        kmers.push_back(kmer);
      }
    });

    return kmers;
  }

 private:
  using Bucket = absl::flat_hash_set<std::bitset<K * 2 - B>>;

  std::array<Bucket, 1 << B> buckets_;

  template <typename F>
  void ForEachBucket(F&& f) {
    std::vector<std::thread> threads;

    for (const auto& range :
         SplitRange(0, 1 << B, std::thread::hardware_concurrency())) {
      threads.emplace_back([&] {
        for (int i = range.first; i < range.second; i++) {
          const Bucket& bucket = buckets_[i];
          f(bucket, i);
        }
      });
    }

    for (auto& thread : threads) thread.join();
  }
};

#endif

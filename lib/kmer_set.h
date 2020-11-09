#ifndef KMER_SET_H_
#define KMER_SET_H_

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
template <int K, typename KeyType>
std::pair<int, KeyType> GetBucketAndKeyFromKmer(const Kmer<K>& kmer) {
  const auto& bits = kmer.Bits();
  const int key_bits = sizeof(KeyType) * 8;
  const int bucket_id = bits.to_ullong() >> key_bits;
  const KeyType key = bits.to_ullong() % (1ull << key_bits);
  return {bucket_id, key};
}

template <int K, typename KeyType>
Kmer<K> GetKmerFromBucketAndKey(int bucketID, KeyType key) {
  const int key_bits = sizeof(KeyType) * 8;
  std::bitset<K * 2> bits{((uint64_t)bucketID << key_bits) + (uint64_t)key};
  return Kmer<K>{bits};
}

template <int K, typename KeyType>
class KmerSet {
 public:
  KmerSet() : buckets_(kBucketsNum) {}

  int64_t Size() const {
    int64_t sum = 0;
    for (const auto& bucket : buckets_) {
      sum += bucket.size();
    }
    return sum;
  }

  void Add(const Kmer<K>& kmer) {
    const auto [bucket, key] = GetBucketAndKeyFromKmer<K, KeyType>(kmer);
    buckets_[bucket].insert(key);
  }

  void Remove(const Kmer<K>& kmer) {
    const auto [bucket, key] = GetBucketAndKeyFromKmer<K, KeyType>(kmer);
    buckets_[bucket].erase(key);
  }

  bool Contains(const Kmer<K>& kmer) const {
    const auto [bucket, key] = GetBucketAndKeyFromKmer<K, KeyType>(kmer);
    return buckets_[bucket].find(key) != buckets_[bucket].end();
  }

  // Fetches one k-mer randomly.
  Kmer<K> Sample() const {
    const int bucket = rand() % kBucketsNum;
    const int64_t n = rand() % buckets_[bucket].size();

    auto it = buckets_[bucket].begin();
    for (int64_t i = 0; i < n; i++) it++;

    const auto key = *it;

    return GetKmerFromBucketAndKey<K, KeyType>(bucket, key);
  }

  // Finds all the k-mers that match the condition.
  template <typename Pred>
  std::vector<Kmer<K>> Find(Pred&& pred, int n_workers) const {
    std::vector<Kmer<K>> kmers;
    std::mutex mu;

    ForEachBucket(
        [&](const Bucket& bucket, int bucket_id) {
          std::vector<Kmer<K>> buf;

          for (const auto& key : bucket) {
            const Kmer<K> kmer =
                GetKmerFromBucketAndKey<K, KeyType>(bucket_id, key);
            if (pred(kmer)) buf.push_back(kmer);
          }

          std::lock_guard _{mu};
          kmers.reserve(kmers.size() + buf.size());
          for (const auto& kmer : buf) {
            kmers.push_back(kmer);
          }
        },
        n_workers);

    return kmers;
  }

  // Returns all the k-mers.
  std::vector<Kmer<K>> Find(int n_workers) const {
    return Find([&](const Kmer<K>&) { return true; }, n_workers);
  }

  KmerSet& Add(const KmerSet& other, int n_workers) {
    other.ForEachBucket(
        [&](const Bucket& other_bucket, int bucket_id) {
          for (const auto& key : other_bucket) {
            buckets_[bucket_id].insert(key);
          }
        },
        n_workers);

    return *this;
  }

  KmerSet& Sub(const KmerSet& other, int n_workers) {
    other.ForEachBucket(
        [&](const Bucket& other_bucket, int bucket_id) {
          for (const auto& key : other_bucket) {
            buckets_[bucket_id].erase(key);
          }
        },
        n_workers);

    return *this;
  }

  static constexpr int kBucketsNumFactor = 2 * K - sizeof(KeyType) * 8;
  static constexpr int kBucketsNum = 1 << (2 * K - sizeof(KeyType) * 8);

 private:
  using Bucket = absl::flat_hash_set<KeyType>;

  std::vector<Bucket> buckets_;

  void Add(int bucket, KeyType key) { buckets_[bucket].insert(key); }

  template <typename F>
  void ForEachBucket(F&& f, int n_workers) const {
    std::vector<std::thread> threads;

    for (const Range& range : Range(0, kBucketsNum).Split(n_workers)) {
      threads.emplace_back([&, range] {
        for (int i = range.begin; i < range.end; i++) {
          const Bucket& bucket = buckets_[i];
          f(bucket, i);
        }
      });
    }

    for (auto& thread : threads) thread.join();
  }

  template <int, typename>
  friend class KmerCounter;
};

template <int K, typename KeyType>
KmerSet<K, KeyType> Add(KmerSet<K, KeyType> lhs, const KmerSet<K, KeyType>& rhs,
                        int n_workers) {
  return lhs.Add(rhs, n_workers);
}

template <int K, typename KeyType>
KmerSet<K, KeyType> Sub(KmerSet<K, KeyType> lhs, const KmerSet<K, KeyType>& rhs,
                        int n_workers) {
  return lhs.Sub(rhs, n_workers);
}

template <int K, typename KeyType>
KmerSet<K, KeyType> Intersection(KmerSet<K, KeyType> lhs,
                                 const KmerSet<K, KeyType>& rhs,
                                 int n_workers) {
  return lhs.Sub(Sub(lhs, rhs, n_workers), n_workers);
}

template <int K, typename KeyType>
bool Equals(const KmerSet<K, KeyType>& lhs, const KmerSet<K, KeyType>& rhs,
            int n_workers) {
  return Sub(lhs, rhs, n_workers).Size() + Sub(rhs, lhs, n_workers).Size() == 0;
}

#endif

#ifndef CORE_KMER_SET_H_
#define CORE_KMER_SET_H_

#include <cstddef>
#include <cstdint>
#include <mutex>
#include <tuple>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/random/random.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "core/kmer.h"
#include "core/range.h"

// The first N bits are used to select buckets.
// The last K * 2 - N bits are used as keys.
template <int K, int N, typename KeyType>
std::pair<int, KeyType> GetBucketAndKeyFromKmer(const Kmer<K>& kmer) {
  const int n_key_bits = K * 2 - N;
  static_assert(n_key_bits <= sizeof(KeyType) * 8);

  const std::uint64_t bits = kmer.Bits();
  const int bucket_id = bits >> n_key_bits;
  const KeyType key = bits % (static_cast<std::uint64_t>(1) << n_key_bits);
  return std::make_pair(bucket_id, key);
}

// Inverse of GetBucketAndKeyFromKmer().
template <int K, int N, typename KeyType>
Kmer<K> GetKmerFromBucketAndKey(int bucket_id, KeyType key) {
  const int n_key_bits = K * 2 - N;
  static_assert(n_key_bits <= sizeof(KeyType) * 8);

  std::uint64_t data = (static_cast<std::uint64_t>(bucket_id) << n_key_bits) +
                       static_cast<std::uint64_t>(key);

  return Kmer<K>(data);
}

// KmerSet can be used to hold a unique set of kmers.
//
// The data is divided to 1 << N buckets and keys of 2 * K - N bits are saved to
// each bucket. The bucket id and the key for a kmer can be calculated with
// GetBucketAndKeyFromKmer().
//
// All the const methods are thread-safe, and it is also safe to add / remove
// kmers from multiple threads if they belong to different buckets.
//
// As data is divided into multiple buckets, some operations can be done
// efficiently in parallel. E.g., it supports multi-threaded operations for
// merging 2 KmerSets and finding all the kmers that match a condition.
template <int K, int N, typename KeyType>
class KmerSet {
  static_assert(2 * K - N <= sizeof(KeyType) * 8);

 public:
  KmerSet() : buckets_(kBucketsNum) {}

  // Returns the total number of stored kmers.
  std::int64_t Size() const {
    std::int64_t sum = 0;
    for (const Bucket& bucket : buckets_) {
      sum += bucket.size();
    }
    return sum;
  }

  // Clears up all the elements in the set.
  void Clear() {
    for (Bucket& bucket : buckets_) {
      Bucket().swap(bucket);
    }
  }

  // Adds a single kmer.
  void Add(const Kmer<K>& kmer) {
    int bucket;
    KeyType key;
    std::tie(bucket, key) = GetBucketAndKeyFromKmer<K, N, KeyType>(kmer);

    buckets_[bucket].insert(key);
  }

  // Removes a single kmer.
  void Remove(const Kmer<K>& kmer) {
    int bucket;
    KeyType key;
    std::tie(bucket, key) = GetBucketAndKeyFromKmer<K, N, KeyType>(kmer);

    buckets_[bucket].erase(key);
  }

  // Returns true if the kmer is contained in the set.
  bool Contains(const Kmer<K>& kmer) const {
    int bucket;
    KeyType key;
    std::tie(bucket, key) = GetBucketAndKeyFromKmer<K, N, KeyType>(kmer);

    return buckets_[bucket].find(key) != buckets_[bucket].end();
  }

  // Allocates memory to accommodate n kmers.
  void Reserve(std::int64_t n) {
    for (int i = 0; i < kBucketsNum; i++) {
      buckets_[i].reserve(n / kBucketsNum);
    }
  }

  // Finds all the k-mers that match the condition.
  // If an estimate of the resulting size is given, it is used for speed-up.
  template <typename PredType>
  std::vector<Kmer<K>> Find(PredType pred, int n_workers,
                            std::int64_t estimated_size = 0) const {
    std::vector<Kmer<K>> kmers;

    if (estimated_size > 0) {
      kmers.reserve(estimated_size);
    }

    std::mutex mu;

    ForEachBucket(
        [&](const Bucket& bucket, int bucket_id) {
          std::vector<Kmer<K>> buf;

          if (n_workers != 1) {
            buf.reserve(bucket.size());
          }

          for (const KeyType& key : bucket) {
            const Kmer<K> kmer =
                GetKmerFromBucketAndKey<K, N, KeyType>(bucket_id, key);

            if (pred(kmer)) {
              if (n_workers == 1) {
                kmers.push_back(kmer);
              } else {
                buf.push_back(kmer);
              }
            }
          }

          if (n_workers != 1) {
            std::lock_guard lck(mu);
            kmers.insert(kmers.end(), buf.begin(), buf.end());
          }
        },
        n_workers);

    return kmers;
  }

  // Returns all the k-mers.
  std::vector<Kmer<K>> Find(int n_workers) const {
    return Find([](const Kmer<K>&) { return true; }, n_workers, Size());
  }

  // Adds kmers that exist in "other".
  KmerSet& Add(const KmerSet& other, int n_workers) {
    other.ForEachBucket(
        [&](const Bucket& other_bucket, int bucket_id) {
          for (const KeyType& key : other_bucket) {
            buckets_[bucket_id].insert(key);
          }
        },
        n_workers);

    return *this;
  }

  // Removes kmers that exist in "other".
  KmerSet& Sub(const KmerSet& other, int n_workers) {
    other.ForEachBucket(
        [&](const Bucket& other_bucket, int bucket_id) {
          for (const KeyType& key : other_bucket) {
            buckets_[bucket_id].erase(key);
          }
        },
        n_workers);

    return *this;
  }

  // Returns the number of kmers which is in only one of the
  // two sets.
  std::int64_t Diff(const KmerSet& other, int n_workers) const {
    std::atomic_int64_t count = 0;

    other.ForEachBucket(
        [&](const Bucket& other_bucket, int bucket_id) {
          for (const KeyType& key : other_bucket) {
            if (buckets_[bucket_id].find(key) == buckets_[bucket_id].end())
              count += 1;
          }
        },
        n_workers);

    ForEachBucket(
        [&](const Bucket& bucket, int bucket_id) {
          for (const KeyType& key : bucket) {
            if (other.buckets_[bucket_id].find(key) ==
                other.buckets_[bucket_id].end())
              count += 1;
          }
        },
        n_workers);

    return count;
  }

  // Returns true if two sets are the same.
  bool Equals(const KmerSet& other, int n_workers) const {
    return Diff(other, n_workers) == 0;
  }

  // Returns the hash value.
  // If two KmerSets contain the same set of kmers, they will produce the same
  // hash value.
  std::size_t Hash(int n_workers) const {
    std::size_t hash = 0;
    std::mutex mu;

    ForEachBucket(
        [&](const Bucket& bucket, int bucket_id) {
          std::size_t buf = 0;

          for (KeyType key : bucket) {
            const Kmer<K> kmer =
                GetKmerFromBucketAndKey<K, N, KeyType>(bucket_id, key);
            buf ^= kmer.Hash();
          }

          std::lock_guard lck(mu);
          hash ^= buf;
        },
        n_workers);

    return hash;
  }

 private:
  using Bucket = absl::flat_hash_set<KeyType>;

  static constexpr int kBucketsNum = 1u << N;

  std::vector<Bucket> buckets_;

  void Add(int bucket, KeyType key) { buckets_[bucket].insert(key); }

  // Executes a function for each bucket (in parallel).
  //
  // Example:
  //   ForEachBucket([&](const Bucket& bucket, int bucket_id) { ... },
  //   n_workers);
  template <typename F>
  void ForEachBucket(F f, int n_workers) const {
    if (n_workers == 1) {
      for (int i = 0; i < kBucketsNum; i++) {
        f(buckets_[i], i);
      }

      return;
    }

    boost::asio::thread_pool pool(n_workers);

    for (const Range& range :
         Range(0, kBucketsNum).Split(n_workers * n_workers)) {
      boost::asio::post(pool, [&, range] {
        for (int i : range) {
          f(buckets_[i], i);
        }
      });
    }

    pool.join();
  }

  template <int, int, typename, typename>
  friend class KmerCounter;
};

// Calculates the union of two KmerSets.
template <int K, int N, typename KeyType>
KmerSet<K, N, KeyType> Add(KmerSet<K, N, KeyType> lhs,
                           const KmerSet<K, N, KeyType>& rhs, int n_workers) {
  return lhs.Add(rhs, n_workers);
}

// Calculates the difference between two KmerSets.
template <int K, int N, typename KeyType>
KmerSet<K, N, KeyType> Sub(KmerSet<K, N, KeyType> lhs,
                           const KmerSet<K, N, KeyType>& rhs, int n_workers) {
  return lhs.Sub(rhs, n_workers);
}

// Calculates the intersection between two KmerSets.
template <int K, int N, typename KeyType>
KmerSet<K, N, KeyType> Intersection(KmerSet<K, N, KeyType> lhs,
                                    const KmerSet<K, N, KeyType>& rhs,
                                    int n_workers) {
  return lhs.Sub(Sub(lhs, rhs, n_workers), n_workers);
}

#endif

#ifndef CORE_KMER_SET_H_
#define CORE_KMER_SET_H_

#include <cstddef>
#include <cstdint>
#include <mutex>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/random/random.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "core/io.h"
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

template <int K, int N, typename KeyType>
Kmer<K> GetKmerFromBucketAndKey(int bucket_id, KeyType key) {
  const int n_key_bits = K * 2 - N;
  static_assert(n_key_bits <= sizeof(KeyType) * 8);

  std::uint64_t data = (static_cast<std::uint64_t>(bucket_id) << n_key_bits) +
                       static_cast<std::uint64_t>(key);

  return Kmer<K>(data);
}

// KmerSet can be used to hold a unique set of kmers.
// As data are divided into multiple buckets, some operations can be done
// efficiently in parallel. E.g., it supports multi-threaded operations for
// merging 2 KmerSets and finding all the kmers that match a condition.
template <int K, int N, typename KeyType>
class KmerSet {
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

  // Finds all the k-mers that match the condition.
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

  // Returns the number of kmers that is in the both of two sets.
  std::int64_t Common(const KmerSet& other, int n_workers) const {
    std::atomic_int64_t count = 0;

    other.ForEachBucket(
        [&](const Bucket& other_bucket, int bucket_id) {
          for (const KeyType& key : other_bucket) {
            if (buckets_[bucket_id].find(key) != buckets_[bucket_id].end())
              count += 1;
          }
        },
        n_workers);

    return count;
  }

  // Returns the Jaccard similarity of two sets.
  double Similarity(const KmerSet& other, int n_workers) const {
    std::int64_t diff = Diff(other, n_workers);
    std::int64_t common = Common(other, n_workers);

    return static_cast<double>(common) / static_cast<double>(diff + common);
  }

  // Returns true if two sets are the same.
  bool Equals(const KmerSet& other, int n_workers) const {
    return Diff(other, n_workers) == 0;
  }

  // Extracts K2-mers from the K-mer set.
  // K2 should be equal to or less than K.
  template <int K2, int N2, typename KeyType2>
  KmerSet<K2, N2, KeyType2> Extract(int n_workers) const {
    std::vector<Kmer<K>> kmers = Find(n_workers);
    KmerSet<K2, N2, KeyType2> extracted;

    std::vector<std::thread> threads;
    int done = 0;
    std::mutex mu;

    for (const Range& range : Range(0, kmers.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        KmerSet<K2, N2, KeyType2> buf;

        for (std::int64_t i : range) {
          const std::string s = kmers[i].String();
          for (int j = 0; j < (int)s.size() - K2 + 1; j++) {
            buf.Add(Kmer<K2>(s.substr(j, K2)));
          }
        }

        std::lock_guard lck(mu);
        extracted.Add(buf, 1 + done);
        done += 1;
      });
    }

    for (std::thread& t : threads) t.join();

    return extracted;
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

  // Executes a function for each bucket (in parallel.)
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

template <int K, int N, typename KeyType>
KmerSet<K, N, KeyType> Add(KmerSet<K, N, KeyType> lhs,
                           const KmerSet<K, N, KeyType>& rhs, int n_workers) {
  return lhs.Add(rhs, n_workers);
}

template <int K, int N, typename KeyType>
KmerSet<K, N, KeyType> Sub(KmerSet<K, N, KeyType> lhs,
                           const KmerSet<K, N, KeyType>& rhs, int n_workers) {
  return lhs.Sub(rhs, n_workers);
}

template <int K, int N, typename KeyType>
KmerSet<K, N, KeyType> Intersection(KmerSet<K, N, KeyType> lhs,
                                    const KmerSet<K, N, KeyType>& rhs,
                                    int n_workers) {
  return lhs.Sub(Sub(lhs, rhs, n_workers), n_workers);
}

#endif

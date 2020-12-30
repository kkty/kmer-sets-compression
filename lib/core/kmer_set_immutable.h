#ifndef CORE_KMER_SET_IMMUTABLE_H_
#define CORE_KMER_SET_IMMUTABLE_H_

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/random/random.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "core/int_set.h"
#include "core/kmer.h"
#include "core/kmer_set.h"

// KmerSetImmutable can be used to hold an immutable set of kmers.
// It is not possible to add or remove kmers, but fast operations are supported
// for adding two KmerSetImmutables, intersecting two KmerSetImmutables, and
// subtracting one KmerSetImmutable from another.
template <int K, int N, typename KeyType>
class KmerSetImmutable {
 public:
  KmerSetImmutable() : buckets_(kBucketsNum) {}

  // Constructs KmerSetImmutable from a list of kmers.
  KmerSetImmutable(const std::vector<Kmer<K>>& kmers, int n_workers)
      : buckets_(kBucketsNum) {
    std::vector<std::thread> threads;
    std::vector<std::mutex> mus(kBucketsNum);

    for (const Range& range : Range(0, kmers.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<std::vector<KeyType>> buf(kBucketsNum);

        for (std::int64_t i : range) {
          int bucket;
          KeyType key;
          std::tie(bucket, key) =
              GetBucketAndKeyFromKmer<K, N, KeyType>(kmers[i]);

          buf[bucket].push_back(key);
        }

        std::vector<IntSet<KeyType>> int_sets(kBucketsNum);

        // Moves data from buf to int_sets.

        for (int i = 0; i < kBucketsNum; i++) {
          std::sort(buf[i].begin(), buf[i].end());
          int_sets[i] = IntSet<KeyType>(buf[i]);
          std::vector<KeyType>().swap(buf[i]);
        }

        std::vector<std::vector<KeyType>>().swap(buf);

        // Moves from int_sets.

        if (n_workers == 1) {
          buckets_ = std::move(int_sets);
        } else {
          std::atomic_int done_count = 0;
          while (done_count < kBucketsNum) {
            for (int i = 0; i < kBucketsNum; i++) {
              if (mus[i].try_lock()) {
                buckets_[i] = buckets_[i].Add(int_sets[i]);
                mus[i].unlock();
                done_count += 1;
              }
            }
          }
        }
      });
    }

    for (std::thread& t : threads) t.join();
  }

  // Constructs KmerSetImmutable from a KmerSet.
  KmerSetImmutable(const KmerSet<K, N, KeyType>& kmer_set, int n_workers)
      : buckets_(kBucketsNum) {
    boost::asio::thread_pool pool(n_workers);

    for (const Range& range :
         Range(0, kBucketsNum).Split(n_workers * n_workers)) {
      boost::asio::post(pool, [&, range] {
        for (int i : range) {
          std::vector<KeyType> v(kmer_set.buckets_[i].begin(),
                                 kmer_set.buckets_[i].end());
          std::sort(v.begin(), v.end());
          buckets_[i] = IntSet<KeyType>(v);
        }
      });
    }

    pool.join();
  }

  // Returns the number of bytes used by buckets.
  std::int64_t Bytes() const {
    std::int64_t bytes = 0;
    for (const Bucket& bucket : buckets_) {
      bytes += bucket.Bytes();
    }
    return bytes;
  }

  void Clear() {
    for (int i = 0; i < kBucketsNum; i++) {
      buckets_[i] = IntSet<KeyType>();
    }
  }

  // Returns the total number of stored kmers.
  std::int64_t Size() const {
    std::int64_t sum = 0;
    for (const Bucket& bucket : buckets_) {
      sum += bucket.Size();
    }
    return sum;
  }

  // Reconstructs a KmerSet.
  KmerSet<K, N, KeyType> ToKmerSet(int n_workers) const {
    KmerSet<K, N, KeyType> kmer_set;

    boost::asio::thread_pool pool(n_workers);

    for (const Range& range :
         Range(0, kBucketsNum).Split(n_workers * n_workers)) {
      boost::asio::post(pool, [&, range] {
        for (int i : range) {
          std::vector<KeyType> v = buckets_[i].Decode();
          kmer_set.buckets_[i] =
              absl::flat_hash_set<KeyType>(v.begin(), v.end());
        }
      });
    }

    pool.join();

    return kmer_set;
  }

  // Reconstructs a list of kmers.
  std::vector<Kmer<K>> ToKmers(int n_workers) const {
    std::vector<Kmer<K>> kmers;

    boost::asio::thread_pool pool(n_workers);
    std::mutex mu;

    for (const Range& range :
         Range(0, kBucketsNum).Split(n_workers * n_workers)) {
      boost::asio::post(pool, [&, range] {
        std::vector<Kmer<K>> buf;

        for (int i : range) {
          std::vector<KeyType> v = buckets_[i].Decode();

          for (KeyType key : v) {
            buf.push_back(GetKmerFromBucketAndKey<K, N, KeyType>(i, key));
          }
        }

        std::lock_guard lck(mu);
        kmers.insert(kmers.end(), buf.begin(), buf.end());
      });
    }

    pool.join();

    return kmers;
  }

  std::int64_t IntersectionSize(const KmerSetImmutable& other,
                                int n_workers) const {
    std::atomic_int64_t count = 0;

    boost::asio::thread_pool pool(n_workers);

    for (const Range& range :
         Range(0, kBucketsNum).Split(n_workers * n_workers)) {
      boost::asio::post(pool, [&, range] {
        for (int i : range) {
          count += buckets_[i].IntersectionSize(other.buckets_[i]);
        }
      });
    }

    pool.join();

    return count;
  }

  // Returns an estimate of the intersection size of two KmerSetImmutables.
  // "n_buckets" buckets are used for estimation.
  std::int64_t IntersectionSizeEstimate(const KmerSetImmutable& other,
                                        int n_buckets) const {
    std::int64_t count = 0;

    absl::InsecureBitGen bitgen;

    for (int i = 0; i < n_buckets; i++) {
      int bucket_id = absl::Uniform(bitgen, 0, kBucketsNum);

      count += buckets_[bucket_id].IntersectionSize(other.buckets_[bucket_id]);
    }

    return count * kBucketsNum / n_buckets;
  }

  // Returns an estimate of the intersection size of two KmerSetImmutables.
  // "(n_buckets_factor * 100)" % of buckets are used for estimation.
  std::int64_t IntersectionSizeEstimate(const KmerSetImmutable& other,
                                        double n_buckets_factor) const {
    return IntersectionSizeEstimate(
        other, static_cast<int>(kBucketsNum * n_buckets_factor));
  }

  KmerSetImmutable Intersection(const KmerSetImmutable& other,
                                int n_workers) const {
    KmerSetImmutable kmer_set_immutable;

    boost::asio::thread_pool pool(n_workers);

    for (const Range& range :
         Range(0, kBucketsNum).Split(n_workers * n_workers)) {
      boost::asio::post(pool, [&, range] {
        for (int i : range) {
          kmer_set_immutable.buckets_[i] =
              buckets_[i].Intersection(other.buckets_[i]);
        }
      });
    }

    pool.join();

    return kmer_set_immutable;
  }

  // Adds two KmerSetImmutables.
  KmerSetImmutable Add(const KmerSetImmutable& other, int n_workers) const {
    KmerSetImmutable kmer_set_immutable;

    boost::asio::thread_pool pool(n_workers);

    for (const Range& range :
         Range(0, kBucketsNum).Split(n_workers * n_workers)) {
      boost::asio::post(pool, [&, range] {
        for (int i : range) {
          kmer_set_immutable.buckets_[i] = buckets_[i].Add(other.buckets_[i]);
        }
      });
    }

    pool.join();

    return kmer_set_immutable;
  }

  // Subtracts two KmerSetImmutables.
  KmerSetImmutable Sub(const KmerSetImmutable& other, int n_workers) const {
    KmerSetImmutable kmer_set_immutable;

    boost::asio::thread_pool pool(n_workers);

    for (const Range& range :
         Range(0, kBucketsNum).Split(n_workers * n_workers)) {
      boost::asio::post(pool, [&, range] {
        for (int i : range) {
          kmer_set_immutable.buckets_[i] = buckets_[i].Sub(other.buckets_[i]);
        }
      });
    }

    pool.join();

    return kmer_set_immutable;
  }

 private:
  using Bucket = IntSet<KeyType>;

  static constexpr int kBucketsNum = 1u << N;

  std::vector<Bucket> buckets_;
};

// KmerSetImmutableHashSet is similar to KmerSetImmutable, but it is backed by a
// simple hash set and does not support parallelism.
template <int K>
class KmerSetImmutableHashSet {
 public:
  template <int N, typename KeyType>
  static KmerSetImmutableHashSet FromKmerSet(
      const KmerSet<K, N, KeyType>& kmer_set, int n_workers) {
    std::vector<Kmer<K>> v = kmer_set.Find(n_workers);
    return KmerSetImmutableHashSet(v);
  }

  KmerSetImmutableHashSet Add(const KmerSetImmutableHashSet& other) const {
    std::vector<Kmer<K>> v(data_.begin(), data_.end());

    for (const Kmer<K>& kmer : other.data_) {
      v.push_back(kmer);
    }

    return KmerSetImmutableHashSet(v);
  }

  KmerSetImmutableHashSet Sub(const KmerSetImmutableHashSet& other) const {
    std::vector<Kmer<K>> v;

    for (const Kmer<K>& kmer : data_) {
      if (other.data_.find(kmer) == other.data_.end()) {
        v.push_back(kmer);
      }
    }

    return KmerSetImmutableHashSet(v);
  }

  KmerSetImmutableHashSet Intersection(
      const KmerSetImmutableHashSet& other) const {
    std::vector<Kmer<K>> v;

    for (const Kmer<K>& kmer : data_) {
      if (other.data_.find(kmer) != other.data_.end()) {
        v.push_back(kmer);
      }
    }

    return KmerSetImmutableHashSet(v);
  }

  std::int64_t Size() const { return data_.size(); }

 private:
  explicit KmerSetImmutableHashSet(const std::vector<Kmer<K>>& v) {
    data_.reserve(v.size());
    data_.insert(v.begin(), v.end());
  }

  absl::flat_hash_set<Kmer<K>> data_;
};

// KmerSetImmutableVector is similar to KmerSetImmutable, but it is backed by a
// simple vector and does not support parallelism.
template <int K>
class KmerSetImmutableVector {
 public:
  template <int N, typename KeyType>
  static KmerSetImmutableVector FromKmerSet(
      const KmerSet<K, N, KeyType>& kmer_set, int n_workers) {
    std::vector<Kmer<K>> data = kmer_set.Find(n_workers);

    std::sort(data.begin(), data.end());

    return KmerSetImmutableVector(std::move(data));
  }

  KmerSetImmutableVector Add(const KmerSetImmutableVector& other) const {
    std::vector<Kmer<K>> data;

    data.reserve(std::max(data_.size(), other.data_.size()));

    std::set_union(data_.begin(), data_.end(), other.data_.begin(),
                   other.data_.end(), std::back_inserter(data));

    return KmerSetImmutableVector(std::move(data));
  }

  KmerSetImmutableVector Sub(const KmerSetImmutableVector& other) const {
    std::vector<Kmer<K>> data;

    std::set_difference(data_.begin(), data_.end(), other.data_.begin(),
                        other.data_.end(), std::back_inserter(data));

    return KmerSetImmutableVector(std::move(data));
  }

  KmerSetImmutableVector Intersection(
      const KmerSetImmutableVector& other) const {
    std::vector<Kmer<K>> data;

    std::set_intersection(data_.begin(), data_.end(), other.data_.begin(),
                          other.data_.end(), std::back_inserter(data));

    return KmerSetImmutableVector(std::move(data));
  }

  std::int64_t Size() const { return data_.size(); }

 private:
  explicit KmerSetImmutableVector(std::vector<Kmer<K>> data)
      : data_(std::move(data)) {}

  std::vector<Kmer<K>> data_;
};

#endif

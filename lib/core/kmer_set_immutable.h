#ifndef CORE_KMER_SET_IMMUTABLE_H_
#define CORE_KMER_SET_IMMUTABLE_H_

#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/random/random.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "core/int_set.h"

// KmerSetImmutable can be used to hold an immutable set of kmers.
// It is not possible to add or remove kmers, but fast operations are supported
// for adding two KmerSetImmutables, intersecting two KmerSetImmutables, and
// subtracting one KmerSetImmutable from another.
template <int K, int N, typename KeyType>
class KmerSetImmutable {
 public:
  KmerSetImmutable() : buckets_(kBucketsNum) {}

  KmerSetImmutable(const KmerSet<K, N, KeyType>& kmer_set, int n_workers)
      : buckets_(kBucketsNum) {
    boost::asio::thread_pool pool(n_workers);

    for (int i = 0; i < kBucketsNum; i++) {
      boost::asio::post(pool, [&, i] {
        std::vector<KeyType> v(kmer_set.buckets_[i].begin(),
                               kmer_set.buckets_[i].end());
        std::sort(v.begin(), v.end());
        buckets_[i] = IntSet<KeyType>(v);
      });
    }

    pool.join();
  }

  void Clear() {
    for (int i = 0; i < kBucketsNum; i++) {
      buckets_[i] = IntSet<KeyType>();
    }
  }

  // Returns the total number of stored kmers.
  int64_t Size() const {
    int64_t sum = 0;
    for (const Bucket& bucket : buckets_) {
      sum += bucket.Size();
    }
    return sum;
  }

  // Reconstructs a KmerSet.
  KmerSet<K, N, KeyType> ToKmerSet(int n_workers) const {
    KmerSet<K, N, KeyType> kmer_set;

    boost::asio::thread_pool pool(n_workers);

    for (int i = 0; i < kBucketsNum; i++) {
      boost::asio::post(pool, [&, i] {
        std::vector<KeyType> v = buckets_[i].Decode();
        kmer_set.buckets_[i] = absl::flat_hash_set<KeyType>(v.begin(), v.end());
      });
    }

    pool.join();

    return kmer_set;
  }

  // Returns an estimate of the intersection size of two KmerSetImmutables.
  // "n_buckets" buckets are used for estimation.
  int64_t IntersectionSizeEstimate(const KmerSetImmutable& other,
                                   int n_buckets) const {
    int64_t count = 0;

    absl::InsecureBitGen bitgen;

    for (int i = 0; i < n_buckets; i++) {
      int bucket_id = absl::Uniform(bitgen, 0, kBucketsNum);

      count += buckets_[bucket_id].IntersectionSize(other.buckets_[bucket_id]);
    }

    return count * kBucketsNum / n_buckets;
  }

  // Returns an estimate of the intersection size of two KmerSetImmutables.
  // "(n_buckets_factor * 100)" % of buckets are used for estimation.
  int64_t IntersectionSizeEstimate(const KmerSetImmutable& other,
                                   double n_buckets_factor) const {
    return IntersectionSizeEstimate(
        other, static_cast<int>(kBucketsNum * n_buckets_factor));
  }

  KmerSetImmutable Intersection(const KmerSetImmutable& other,
                                int n_workers) const {
    KmerSetImmutable kmer_set_immutable;

    boost::asio::thread_pool pool(n_workers);

    for (int i = 0; i < kBucketsNum; i++) {
      boost::asio::post(pool, [&, i] {
        kmer_set_immutable.buckets_[i] =
            buckets_[i].Intersection(other.buckets_[i]);
      });
    }

    pool.join();

    return kmer_set_immutable;
  }

  // Adds two KmerSetImmutables.
  KmerSetImmutable Add(const KmerSetImmutable& other, int n_workers) const {
    KmerSetImmutable kmer_set_immutable;

    boost::asio::thread_pool pool(n_workers);

    for (int i = 0; i < kBucketsNum; i++) {
      boost::asio::post(pool, [&, i] {
        kmer_set_immutable.buckets_[i] = buckets_[i].Add(other.buckets_[i]);
      });
    }

    pool.join();

    return kmer_set_immutable;
  }

  // Subtracts two KmerSetImmutables.
  KmerSetImmutable Sub(const KmerSetImmutable& other, int n_workers) const {
    KmerSetImmutable kmer_set_immutable;

    boost::asio::thread_pool pool(n_workers);

    for (int i = 0; i < kBucketsNum; i++) {
      boost::asio::post(pool, [&, i] {
        kmer_set_immutable.buckets_[i] = buckets_[i].Sub(other.buckets_[i]);
      });
    }

    pool.join();

    return kmer_set_immutable;
  }

 private:
  using Bucket = IntSet<KeyType>;

  static constexpr int kBucketsNum = 1ull << (2 * K - sizeof(KeyType) * 8);

  std::vector<Bucket> buckets_;
};

#endif

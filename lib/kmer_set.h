#ifndef KMER_SET_H_
#define KMER_SET_H_

#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "io.h"
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

  KmerSet(const std::vector<Kmer<K>>& kmers, int n_workers)
      : buckets_(kBucketsNum) {
    std::vector<std::thread> threads;
    std::vector<std::mutex> mus(kBucketsNum);

    for (const Range& range : Range(0, kmers.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<Bucket> buf(kBucketsNum);

        range.ForEach([&](int64_t i) {
          const Kmer<K>& kmer = kmers[i];
          const auto [bucket, key] = GetBucketAndKeyFromKmer<K, KeyType>(kmer);

          buf[bucket].insert(key);
        });

        std::vector<bool> done(kBucketsNum);
        int done_count = 0;

        while (done_count < kBucketsNum) {
          for (int i = 0; i < kBucketsNum; i++) {
            if (done[i]) continue;

            if (mus[i].try_lock()) {
              for (const auto& key : buf[i]) buckets_[i].insert(key);

              mus[i].unlock();

              done[i] = true;
              done_count += 1;
            }
          }
        }
      });
    }

    for (auto& t : threads) t.join();
  }

  int64_t Size() const {
    int64_t sum = 0;
    for (const auto& bucket : buckets_) {
      sum += bucket.size();
    }
    return sum;
  }

  // Adds a single kmer.
  void Add(const Kmer<K>& kmer) {
    const auto [bucket, key] = GetBucketAndKeyFromKmer<K, KeyType>(kmer);
    buckets_[bucket].insert(key);
  }

  // Removes a single kmer.
  void Remove(const Kmer<K>& kmer) {
    const auto [bucket, key] = GetBucketAndKeyFromKmer<K, KeyType>(kmer);
    buckets_[bucket].erase(key);
  }

  // Returns true if the kmer is contained in the set.
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
  std::vector<Kmer<K>> Find(Pred pred, int n_workers) const {
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

          std::lock_guard lck(mu);
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

  // Adds kmers that exist in "other".
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

  // Removes kmers that exist in "other".
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

  // Returns the number of kmers which is in only one of the
  // two sets.
  int64_t Diff(const KmerSet& other, int n_workers) const {
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
  int64_t Common(const KmerSet& other, int n_workers) const {
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

  // Approximates the number of common kmers by using some buckets.
  int64_t CommonEstimate(const KmerSet& other, int n_buckets) const {
    int64_t count = 0;

    for (int i = 0; i < n_buckets; i++) {
      for (const KeyType& key : other.buckets_[i]) {
        if (buckets_[i].find(key) != buckets_[i].end()) count += 1;
      }
    }

    return count * kBucketsNum / n_buckets;
  }

  // Approximates the number of common kmers by using some buckets.
  // If "n_buckets_factor" is 0.3, 30% of buckets are used.
  int64_t CommonEstimate(const KmerSet& other, double n_buckets_factor) const {
    return CommonEstimate(other,
                          static_cast<int>(kBucketsNum * n_buckets_factor));
  }

  // Returns the Jaccard similarity of two sets.
  double Similarity(const KmerSet& other, int n_workers) const {
    int64_t diff = Diff(other, n_workers);
    int64_t common = Common(other, n_workers);

    return (double)common / (diff + common);
  }

  // Estimates the difference between two sets using some buckets.
  // n_buckets should be no more than kBucketsNum.
  int64_t DiffEstimate(const KmerSet& other, int n_buckets) const {
    int64_t count = 0;

    for (int i = 0; i < n_buckets; i++) {
      const Bucket& bucket = buckets_[i];
      const Bucket& other_bucket = other.buckets_[i];

      for (const KeyType& key : bucket) {
        if (other_bucket.find(key) == other_bucket.end()) count += 1;
      }

      for (const KeyType& key : other_bucket) {
        if (bucket.find(key) == bucket.end()) count += 1;
      }
    }

    return count * kBucketsNum / n_buckets;
  }

  // Returns true if two sets are the same.
  bool Equals(const KmerSet& other, int n_workers) const {
    return Diff(other, n_workers) == 0;
  }

  // Dumps kmers to a file.
  void Dump(const std::string& file_name, const std::string& compressor,
            int n_workers) const {
    std::vector<Kmer<K>> kmers = Find(n_workers);

    std::vector<std::string> lines(kmers.size());

    std::vector<std::thread> threads;

    for (const Range& range : Range(0, kmers.size()).Split(1)) {
      threads.emplace_back([&, range] {
        range.ForEach([&](int64_t i) { lines[i] = kmers[i].String(); });
      });
    }

    for (auto& t : threads) t.join();

    WriteLines(file_name, compressor, lines);
  }

  // Loads kmers from a file.
  static KmerSet Load(const std::string& file_name,
                      const std::string& decompressor, int n_workers) {
    std::vector<std::string> lines = ReadLines(file_name, decompressor);
    std::vector<Kmer<K>> kmers(lines.size());

    std::vector<std::thread> threads;

    for (const Range& range : Range(0, lines.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        range.ForEach([&](int i) {
          kmers[i] = Kmer<K>(lines[i]);

          // lines[i] is not needed anymore.
          std::string().swap(lines[i]);
        });
      });
    }

    for (std::thread& t : threads) t.join();

    return KmerSet(kmers, n_workers);
  }

  static constexpr int kBucketsNumFactor = 2 * K - sizeof(KeyType) * 8;
  static constexpr int kBucketsNum = 1 << (2 * K - sizeof(KeyType) * 8);

 private:
  using Bucket = absl::flat_hash_set<KeyType>;

  std::vector<Bucket> buckets_;

  void Add(int bucket, KeyType key) { buckets_[bucket].insert(key); }

  // Execute a function for each bucket.
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

    for (int i = 0; i < kBucketsNum; i++) {
      boost::asio::post(pool, [&, i] { f(buckets_[i], i); });
    }

    pool.join();
  }

  template <int, typename, typename>
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

#endif

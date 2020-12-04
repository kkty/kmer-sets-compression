#ifndef CORE_KMER_SET_H_
#define CORE_KMER_SET_H_

#include <bitset>
#include <mutex>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/random/random.h"
#include "absl/status/statusor.h"
#include "absl/status/status.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "core/io.h"
#include "core/kmer.h"
#include "core/range.h"

// The last sizeof(KeyType) * 8 bits are used as keys.
// The other bits are used to select buckets.
template <int K, typename KeyType>
std::pair<int, KeyType> GetBucketAndKeyFromKmer(const Kmer<K>& kmer) {
  const std::bitset<2 * K> bits = kmer.Bits();
  const int key_bits = sizeof(KeyType) * 8;
  const int bucket_id = bits.to_ullong() >> key_bits;
  const KeyType key = bits.to_ullong() % (1ull << key_bits);
  return {bucket_id, key};
}

template <int K, typename KeyType>
Kmer<K> GetKmerFromBucketAndKey(int bucket_id, KeyType key) {
  const int key_bits = sizeof(KeyType) * 8;
  std::bitset<K * 2> bits{((uint64_t)bucket_id << key_bits) + (uint64_t)key};
  return Kmer<K>{bits};
}

// KmerSet can be used to hold a unique set of kmers.
// As data are divided into multiple buckets, some operations can be done
// efficiently in parallel. E.g., it supports multi-threaded operations for
// merging 2 KmerSets and finding all the kmers that match a condition.
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
          int bucket;
          KeyType key;
          std::tie(bucket, key) = GetBucketAndKeyFromKmer<K, KeyType>(kmers[i]);
          buf[bucket].insert(key);
        });

        // Moves data from buffer.

        std::vector<bool> done(kBucketsNum);
        int done_count = 0;

        while (done_count < kBucketsNum) {
          for (int i = 0; i < kBucketsNum; i++) {
            if (done[i]) continue;

            if (mus[i].try_lock()) {
              for (const KeyType& key : buf[i]) buckets_[i].insert(key);

              mus[i].unlock();

              done[i] = true;
              done_count += 1;
            }
          }
        }
      });
    }

    for (std::thread& t : threads) t.join();
  }

  // Returns the total number of stored kmers.
  int64_t Size() const {
    int64_t sum = 0;
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
    std::tie(bucket, key) = GetBucketAndKeyFromKmer<K, KeyType>(kmer);

    buckets_[bucket].insert(key);
  }

  // Removes a single kmer.
  void Remove(const Kmer<K>& kmer) {
    int bucket;
    KeyType key;
    std::tie(bucket, key) = GetBucketAndKeyFromKmer<K, KeyType>(kmer);

    buckets_[bucket].erase(key);
  }

  // Returns true if the kmer is contained in the set.
  bool Contains(const Kmer<K>& kmer) const {
    int bucket;
    KeyType key;
    std::tie(bucket, key) = GetBucketAndKeyFromKmer<K, KeyType>(kmer);

    return buckets_[bucket].find(key) != buckets_[bucket].end();
  }

  // Finds all the k-mers that match the condition.
  template <typename PredType>
  std::vector<Kmer<K>> Find(PredType pred, int n_workers) const {
    std::vector<Kmer<K>> kmers;
    std::mutex mu;

    ForEachBucket(
        [&](const Bucket& bucket, int bucket_id) {
          std::vector<Kmer<K>> buf;

          for (const KeyType& key : bucket) {
            const Kmer<K> kmer =
                GetKmerFromBucketAndKey<K, KeyType>(bucket_id, key);
            if (pred(kmer)) buf.push_back(kmer);
          }

          std::lock_guard lck(mu);
          kmers.insert(kmers.end(), buf.begin(), buf.end());
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

    absl::InsecureBitGen bitgen;

    for (int i = 0; i < n_buckets; i++) {
      int bucket_id = absl::Uniform(bitgen, 0, kBucketsNum);

      for (const KeyType& key : other.buckets_[bucket_id]) {
        if (buckets_[bucket_id].find(key) != buckets_[bucket_id].end())
          count += 1;
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

    absl::InsecureBitGen bitgen;

    for (int i = 0; i < n_buckets; i++) {
      int bucket_id = absl::Uniform(bitgen, 0, kBucketsNum);

      const Bucket& bucket = buckets_[bucket_id];
      const Bucket& other_bucket = other.buckets_[bucket_id];

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

  // Extracts K2-mers from the K-mer set.
  // K2 should be equal to or less than K.
  template <int K2, typename KeyType2>
  KmerSet<K2, KeyType2> Extract(int n_workers) const {
    std::vector<Kmer<K>> kmers = Find(n_workers);
    KmerSet<K2, KeyType2> extracted;

    std::vector<std::thread> threads;
    int done = 0;
    std::mutex mu;

    for (const Range& range : Range(0, kmers.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        KmerSet<K2, KeyType2> buf;

        range.ForEach([&](int64_t i) {
          const std::string s = kmers[i].String();
          for (int j = 0; j < (int)s.size() - K2 + 1; j++) {
            buf.Add(Kmer<K2>(s.substr(j, K2)));
          }
        });

        std::lock_guard lck(mu);
        extracted.Add(buf, 1 + done);
        done += 1;
      });
    }

    for (std::thread& t : threads) t.join();

    return extracted;
  }

  // Dumps kmers to a file.
  absl::Status Dump(const std::string& file_name, const std::string& compressor,
            int n_workers) const {
    std::vector<Kmer<K>> kmers = Find(n_workers);

    std::vector<std::string> lines(kmers.size());

    std::vector<std::thread> threads;

    for (const Range& range : Range(0, kmers.size()).Split(1)) {
      threads.emplace_back([&, range] {
        range.ForEach([&](int64_t i) { lines[i] = kmers[i].String(); });
      });
    }

    for (std::thread& t : threads) t.join();

    return WriteLines(file_name, compressor, lines);
  }

  // Returns the hash value.
  // If two KmerSets contain the same set of kmers, they will produce the same
  // hash value.
  size_t Hash(int n_workers) const {
    size_t hash = 0;
    std::mutex mu;

    ForEachBucket(
        [&](const Bucket& bucket, int bucket_id) {
          size_t buf = 0;

          for (KeyType key : bucket) {
            const Kmer<K> kmer =
                GetKmerFromBucketAndKey<K, KeyType>(bucket_id, key);
            buf ^= kmer.Hash();
          }

          std::lock_guard lck(mu);
          hash ^= buf;
        },
        n_workers);

    return hash;
  }

  // Loads kmers from a file.
  static absl::StatusOr<KmerSet> Load(const std::string& file_name,
                                      const std::string& decompressor,
                                      int n_workers) {
    std::vector<std::string> lines;

    {
      absl::StatusOr<std::vector<std::string>> statusor =
          ReadLines(file_name, decompressor);

      if (!statusor.ok()) {
        return statusor.status();
      }

      lines = std::move(statusor).value();
    }

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

 private:
  using Bucket = absl::flat_hash_set<KeyType>;

  static constexpr int kBucketsNum = 1 << (2 * K - sizeof(KeyType) * 8);

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

    for (int i = 0; i < kBucketsNum; i++) {
      boost::asio::post(pool, [&, i] { f(buckets_[i], i); });
    }

    pool.join();
  }

  template <int, typename, typename>
  friend class KmerCounter;
};

template <int K, typename KeyType>
KmerSet<K, KeyType> Add(const std::vector<KmerSet<K, KeyType>>& kmer_sets,
                        int begin, int end, int n_workers) {
  if (end - begin == 1) return kmer_sets[begin];

  const int mid = (begin + end) / 2;

  if (n_workers == 1) {
    return Add(Add(kmer_sets, begin, mid, 1), Add(kmer_sets, mid, end, 1), 1);
  }

  KmerSet<K, KeyType> lhs;
  KmerSet<K, KeyType> rhs;

  std::vector<std::thread> threads;

  threads.emplace_back(
      [&] { lhs = Add(kmer_sets, begin, mid, n_workers / 2); });

  threads.emplace_back(
      [&] { rhs = Add(kmer_sets, mid, end, n_workers - n_workers / 2); });

  for (std::thread& t : threads) t.join();

  return Add(std::move(lhs), rhs, n_workers);
}

template <int K, typename KeyType>
KmerSet<K, KeyType> Add(const std::vector<KmerSet<K, KeyType>>& kmer_sets,
                        int n_workers) {
  return Add(kmer_sets, n_workers, 0, kmer_sets.size(), n_workers);
}

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

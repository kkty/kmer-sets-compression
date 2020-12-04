#ifndef CORE_KMER_COUNTER_H_
#define CORE_KMER_COUNTER_H_

#include <atomic>
#include <limits>
#include <mutex>
#include <random>
#include <string>
#include <thread>
#include <tuple>
#include <type_traits>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_split.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "core/io.h"
#include "core/kmer.h"
#include "core/kmer_set.h"
#include "core/range.h"

// Adds two numerical values. If the result exceeds the max value for the type,
// the latter is returned.
template <typename T>
T AddWithMax(T x, T y) {
  if (std::is_integral<T>::value) {
    T max = std::numeric_limits<T>::max();
    return std::min((int64_t)max, (int64_t)x + (int64_t)y);
  }

  return x + y;
}

// Multiplies two numerical values. If the result exceeds the max value for the
// type, the latter is returned.
template <typename T>
T MultiplyWithMax(T x, T y) {
  if (std::is_integral<T>::value) {
    T max = std::numeric_limits<T>::max();
    return std::min((int64_t)max, (int64_t)x * (int64_t)y);
  }

  return x * y;
}

// KmerCounter can be used to count kmers.
template <int K, typename KeyType, typename ValueType = uint8_t>
class KmerCounter {
 public:
  KmerCounter() : buckets_(kBucketsNum) {}

  // Returns the number of distinct k-mers.
  int64_t Size() const {
    int64_t sum = 0;
    for (const Bucket& bucket : buckets_) {
      sum += bucket.size();
    }
    return sum;
  }

  static KmerCounter FromReads(std::vector<std::string> reads, bool canonical,
                               int n_workers) {
    KmerCounter kmer_counter;

    std::vector<std::mutex> buckets_mutexes(kBucketsNum);
    std::vector<std::thread> threads;

    for (const Range& range : Range(0, reads.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        std::vector<Bucket> buf(kBucketsNum);

        range.ForEach([&](int64_t i) {
          std::string& read = reads[i];
          std::vector<std::string> fragments = absl::StrSplit(read, "N");

          for (const std::string& fragment : fragments) {
            for (size_t i = 0; i + K <= fragment.length(); i++) {
              const Kmer<K> kmer(fragment.substr(i, K));

              int bucket;
              KeyType key;

              std::tie(bucket, key) = GetBucketAndKeyFromKmer<K, KeyType>(
                  canonical ? kmer.Canonical() : kmer);

              buf[bucket][key] = AddWithMax<ValueType>(buf[bucket][key], 1);
            }
          }

          // "read" is not needed anymore.
          std::string().swap(read);
        });

        std::vector<bool> done(kBucketsNum);
        int done_count = 0;

        while (done_count < kBucketsNum) {
          for (int i = 0; i < kBucketsNum; i++) {
            if (done[i]) continue;

            std::mutex& mu = buckets_mutexes[i];

            if (mu.try_lock()) {
              Bucket& bucket = kmer_counter.buckets_[i];

              for (const std::pair<const KeyType, ValueType>& p : buf[i]) {
                bucket[p.first] = AddWithMax(bucket[p.first], p.second);
              }

              mu.unlock();
              done[i] = true;
              done_count += 1;
            }
          }
        }
      });
    }

    for (std::thread& thread : threads) thread.join();

    return kmer_counter;
  }

  // Counts kmers in a FASTA file.
  static absl::StatusOr<KmerCounter> FromFASTA(const std::string& file_name,
                                               bool canonical, int n_workers) {
    std::vector<std::string> lines;

    {
      absl::StatusOr<std::vector<std::string>> statusor = ReadLines(file_name);

      if (!statusor.ok()) {
        return statusor.status();
      }

      lines = std::move(statusor).value();
    }

    return FromFASTA(std::move(lines), canonical, n_workers);
  }

  // Counts kmers in a FASTA file. "decompressor" is used to decompress the
  // file.
  static absl::StatusOr<KmerCounter> FromFASTA(const std::string& file_name,
                                               const std::string& decompressor,
                                               bool canonical, int n_workers) {
    std::vector<std::string> lines;

    {
      absl::StatusOr<std::vector<std::string>> statusor =
          ReadLines(file_name, decompressor);

      if (!statusor.ok()) {
        return statusor.status();
      }

      lines = std::move(statusor).value();
    }

    return FromFASTA(std::move(lines), canonical, n_workers);
  }

  // Counts kmers in a FASTA data.
  static absl::StatusOr<KmerCounter> FromFASTA(std::vector<std::string> lines,
                                               bool canonical, int n_workers) {
    if (lines.size() % 2 != 0) {
      return absl::FailedPreconditionError(
          "FASTA files should have an even number of lines");
    }

    std::vector<std::string> reads;
    reads.reserve(lines.size() / 2);

    for (size_t i = 0; i < lines.size(); i++) {
      if (i % 2 == 0 && lines[i][0] != '>') {
        return absl::FailedPreconditionError("invalid FASTA file");
      }

      if (i % 2 == 1) {
        reads.push_back(std::move(lines[i]));
      } else {
        std::string().swap(lines[i]);
      }
    }

    // "lines" is not needed anymore.
    std::vector<std::string>().swap(lines);

    return FromReads(std::move(reads), canonical, n_workers);
  }

  // Returns a KmerSet, ignoring ones that appear less often.
  // The number of ignored k-mers is also returned.
  std::pair<KmerSet<K, KeyType>, int64_t> ToKmerSet(ValueType cutoff,
                                                    int n_workers) const {
    KmerSet<K, KeyType> set;
    std::mutex mu;

    std::vector<std::thread> threads;
    std::atomic_int64_t cutoff_count = 0;

    ForEachBucket(
        [&](const Bucket& bucket, int bucket_id) {
          std::vector<KeyType> buf;

          for (const std::pair<const KeyType, ValueType>& p : bucket) {
            KeyType key;
            ValueType count;
            std::tie(key, count) = p;

            if (count < cutoff) {
              cutoff_count += 1;
              continue;
            }

            buf.push_back(key);
          }

          {
            std::lock_guard lck(mu);
            for (KeyType key : buf) {
              set.buckets_[bucket_id].insert(key);
            }
          }
        },
        n_workers);

    return {set, (int64_t)cutoff_count};
  }

  // Returns the count for a kmer.
  ValueType Get(const Kmer<K>& kmer) const {
    int bucket;
    KeyType key;
    std::tie(bucket, key) = GetBucketAndKeyFromKmer<K, KeyType>(kmer);

    auto it = buckets_[bucket].find(key);
    if (it == buckets_[bucket].end()) return 0;
    return it->second;
  }

  // Increments the count for a kmer.
  KmerCounter& Add(const Kmer<K>& kmer, ValueType v) {
    int bucket;
    KeyType key;
    std::tie(bucket, key) = GetBucketAndKeyFromKmer<K, KeyType>(kmer);

    buckets_[bucket][key] = AddWithMax(buckets_[bucket][key], v);
    return *this;
  }

  KmerCounter& Add(const KmerCounter& other, int n_workers) {
    other.ForEachBucket(
        [&](const Bucket& other_bucket, int bucket_id) {
          for (const std::pair<const KeyType, ValueType>& p : other_bucket) {
            KeyType key;
            ValueType value;
            std::tie(key, value) = p;

            buckets_[bucket_id][key] =
                AddWithMax(buckets_[bucket_id][key], value);
          }
        },
        n_workers);

    return *this;
  }

  KmerCounter& Multiply(ValueType v, int n_workers) {
    ForEachBucket(
        [&](const Bucket& bucket, int bucket_id) {
          for (const std::pair<const KeyType, ValueType>& p : bucket) {
            KeyType key;
            ValueType value;
            std::tie(key, value) = p;

            buckets_[bucket_id][key] = MultiplyWithMax(value, v);
          }
        },
        n_workers);

    return *this;
  }

  // Constructs a KmerCounter from a KmerSet. The counts are set to 1.
  static KmerCounter FromKmerSet(const KmerSet<K, KeyType>& kmer_set,
                                 int n_workers) {
    KmerCounter kmer_counter;

    kmer_set.ForEachBucket(
        [&](const typename KmerSet<K, KeyType>::Bucket& bucket, int bucket_id) {
          kmer_counter.buckets_[bucket_id].reserve(bucket.size());

          for (const KeyType& key : bucket) {
            kmer_counter.buckets_[bucket_id][key] = 1;
          }
        },
        n_workers);

    return kmer_counter;
  }

 private:
  static constexpr int kBucketsNum = 1 << (2 * K - sizeof(KeyType) * 8);

  using Bucket = absl::flat_hash_map<KeyType, ValueType>;

  std::vector<Bucket> buckets_;

  // Dispatches processing for each bucket.
  // Usage:
  //   ForEachBucket([&](const Bucket& bucket, int bucket_id){ ... }, 4);
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
};

#endif

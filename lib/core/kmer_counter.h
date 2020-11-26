#ifndef CORE_KMER_COUNTER_H_
#define CORE_KMER_COUNTER_H_

#include <atomic>
#include <limits>
#include <mutex>
#include <random>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

#include "absl/container/flat_hash_map.h"
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
  T max = std::numeric_limits<T>::max();

  if (std::is_floating_point<T>::value) {
    return std::min((double)max, (double)x + (double)y);
  } else {
    return std::min((int64_t)max, (int64_t)x + (int64_t)y);
  }
}

// Multiplies two numerical values. If the result exceeds the max value for the
// type, the latter is returned.
template <typename T>
T MultiplyWithMax(T x, T y) {
  T max = std::numeric_limits<T>::max();

  if (std::is_floating_point<T>::value) {
    return std::min((double)max, (double)x * (double)y);
  } else {
    return std::min((int64_t)max, (int64_t)x * (int64_t)y);
  }
}

// KmerCounter can be used to count kmers.
template <int K, typename KeyType, typename ValueType = uint8_t>
class KmerCounter {
 public:
  KmerCounter() : buckets_(kBucketsNum) {}

  // Returns the number of distinct k-mers.
  int64_t Size() const {
    int64_t sum = 0;
    for (const auto& bucket : buckets_) {
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

          for (const auto& fragment : fragments) {
            for (size_t i = 0; i + K <= fragment.length(); i++) {
              const Kmer<K> kmer(fragment.substr(i, K));

              const auto [bucket, key] = GetBucketAndKeyFromKmer<K, KeyType>(
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

            auto& mu = buckets_mutexes[i];

            if (mu.try_lock()) {
              auto& bucket = kmer_counter.buckets_[i];

              for (const auto& [key, count] : buf[i]) {
                bucket[key] = AddWithMax(bucket[key], count);
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

  static absl::StatusOr<KmerCounter> FromFASTQ(const std::string& file_name,
                                               bool canonical, int n_workers) {
    return FromFASTQ(ReadLines(file_name), canonical, n_workers);
  }

  static absl::StatusOr<KmerCounter> FromFASTQ(const std::string& file_name,
                                               const std::string& decompressor,
                                               bool canonical, int n_workers) {
    return FromFASTQ(ReadLines(file_name, decompressor), canonical, n_workers);
  }

  // Counts kmers in a FASTQ data.
  static absl::StatusOr<KmerCounter> FromFASTQ(std::vector<std::string> lines,
                                               bool canonical, int n_workers) {
    if (lines.size() % 4 != 0)
      return absl::UnknownError("there should be 4 * N lines");

    std::vector<std::string> reads;
    reads.reserve(lines.size() / 4);

    for (size_t i = 0; i < lines.size(); i++) {
      if (i % 4 == 0 && lines[i][0] != '@')
        return absl::UnknownError("the line should start with '@'");

      if (i % 4 == 1)
        reads.push_back(std::move(lines[i]));
      else
        std::string().swap(lines[i]);
    }

    // "lines" is not needed anymore.
    std::vector<std::string>().swap(lines);

    return FromReads(std::move(reads), canonical, n_workers);
  }

  // Counts kmers in a FASTA file.
  static absl::StatusOr<KmerCounter> FromFASTA(const std::string& file_name,
                                               bool canonical, int n_workers) {
    return FromFASTA(ReadLines(file_name), canonical, n_workers);
  }

  // Counts kmers in a FASTA file. "decompressor" is used to decompress the
  // file.
  static absl::StatusOr<KmerCounter> FromFASTA(const std::string& file_name,
                                               const std::string& decompressor,
                                               bool canonical, int n_workers) {
    return FromFASTA(ReadLines(file_name, decompressor), canonical, n_workers);
  }

  // Counts kmers in a FASTA data.
  static absl::StatusOr<KmerCounter> FromFASTA(std::vector<std::string> lines,
                                               bool canonical, int n_workers) {
    if (lines.size() % 2 != 0)
      return absl::UnknownError("there should be 2 * N lines");

    std::vector<std::string> reads;
    reads.reserve(lines.size() / 2);

    for (size_t i = 0; i < lines.size(); i++) {
      if (i % 2 == 0 && lines[i][0] != '>')
        return absl::UnknownError("the line should start with '>'");

      if (i % 2 == 1)
        reads.push_back(std::move(lines[i]));
      else
        std::string().swap(lines[i]);
    }

    // "lines" is not needed anymore.
    std::vector<std::string>().swap(lines);

    return FromReads(std::move(reads), canonical, n_workers);
  }

  // Returns a KmerSet, ignoring ones that appear less often.
  // The number of ignored k-mers is also returned.
  std::pair<KmerSet<K, KeyType>, int64_t> ToSet(ValueType cutoff,
                                                int n_workers) const {
    KmerSet<K, KeyType> set;
    std::mutex mu;

    std::vector<std::thread> threads;
    std::atomic_int64_t cutoff_count = 0;

    ForEachBucket(
        [&](const Bucket& bucket, int bucket_id) {
          std::vector<KeyType> buf;

          for (const auto& [key, count] : bucket) {
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

  ValueType Get(const Kmer<K>& kmer) const {
    const auto [bucket, key] = GetBucketAndKeyFromKmer<K, KeyType>(kmer);
    auto it = buckets_[bucket].find(key);
    if (it == buckets_[bucket].end()) return 0;
    return it->second;
  }

  KmerCounter& Add(const Kmer<K>& kmer, ValueType count) {
    const auto [bucket, key] = GetBucketAndKeyFromKmer<K, KeyType>(kmer);
    buckets_[bucket][key] = AddWithMax(buckets_[bucket][key], count);
    return *this;
  }

  KmerCounter& Add(const KmerCounter& other, int n_workers) {
    other.ForEachBucket(
        [&](const Bucket& other_bucket, int bucket_id) {
          for (const auto& [key, value] : other_bucket) {
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
          for (const auto& [key, value] : bucket) {
            buckets_[bucket_id][key] = MultiplyWithMax(value, v);
          }
        },
        n_workers);

    return *this;
  }

  static KmerCounter FromSet(const KmerSet<K, KeyType>& kmer_set,
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

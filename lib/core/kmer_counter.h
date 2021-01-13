#ifndef CORE_KMER_COUNTER_H_
#define CORE_KMER_COUNTER_H_

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <mutex>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
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
    return std::min(
        static_cast<std::int64_t>(max),
        static_cast<std::int64_t>(x) + static_cast<std::int64_t>(y));
  }

  return x + y;
}

// KmerCounter can be used to count kmers.
//
// The data is divided to 1 << N buckets and keys of 2 * K - N bits (in the form
// of KeyType) and values of ValueType are saved to each bucket.
//
// As data is divided into multiple buckets, some operations are performed
// efficiently in parallel, including the conversion to KmerSet and the counting
// of kmers in FASTA files.
template <int K, int N, typename KeyType, typename ValueType = std::uint8_t>
class KmerCounter {
 public:
  KmerCounter() : buckets_(kBucketsNum) {}

  // Returns the number of distinct k-mers.
  std::int64_t Size() const {
    std::int64_t sum = 0;
    for (const Bucket& bucket : buckets_) {
      sum += bucket.size();
    }
    return sum;
  }

  // Constructs KmerCounter from a list of reads.
  // Each read should only contain 'A', 'C', 'G', 'T', or 'N'.
  static KmerCounter FromReads(std::vector<std::string> reads, bool canonical,
                               int n_workers) {
    KmerCounter kmer_counter;

    std::vector<std::mutex> buckets_mutexes(kBucketsNum);
    boost::asio::thread_pool pool(n_workers);

    for (const Range& range :
         Range(0, reads.size()).Split(n_workers * n_workers)) {
      boost::asio::post(pool, [&, range] {
        std::vector<Bucket> buf(kBucketsNum);

        for (std::int64_t i : range) {
          std::string& read = reads[i];
          std::vector<std::string> fragments = absl::StrSplit(read, 'N');

          for (const std::string& fragment : fragments) {
            for (std::size_t j = 0; j + K <= fragment.length(); j++) {
              const Kmer<K> kmer(fragment.substr(j, K));

              int bucket;
              KeyType key;

              std::tie(bucket, key) = GetBucketAndKeyFromKmer<K, N, KeyType>(
                  canonical ? kmer.Canonical() : kmer);

              buf[bucket][key] = AddWithMax<ValueType>(buf[bucket][key], 1);
            }
          }

          // "read" is not needed anymore.
          std::string().swap(read);
        }

        if (n_workers == 1) {
          kmer_counter.buckets_ = std::move(buf);
          return;
        }

        // Moves from buffer.

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

    pool.join();

    return kmer_counter;
  }

  // Counts kmers in a FASTA file. "decompressor" is used to decompress the
  // file.
  // Example:
  //   ... = KmerCounter::FromFASTA("foo.fasta", "", canonical, n_workers)
  //   ... = KmerCounter::FromFASTA("foo.fasta.bz2", "bzip2 -d", canonical,
  //                                n_workers)
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

    std::vector<std::string> reads(lines.size() / 2);
    boost::asio::thread_pool pool(n_workers);
    std::atomic_bool is_valid = true;

    for (const Range& range :
         Range(0, lines.size()).Split(n_workers * n_workers)) {
      boost::asio::post(pool, [&, range] {
        for (int i : range) {
          std::string& line = lines[i];

          if (i % 2 == 0) {
            if (line.empty() || line[0] != '>') {
              is_valid = false;
              return;
            }

            std::string().swap(line);
          } else {
            for (char c : line) {
              if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') {
                is_valid = false;
                return;
              }
            }

            reads[i / 2] = std::move(line);
          }
        }
      });
    }

    pool.join();

    if (!is_valid) {
      return absl::FailedPreconditionError("invalid FASTA file");
    }

    // "lines" is not needed anymore.
    std::vector<std::string>().swap(lines);

    return FromReads(std::move(reads), canonical, n_workers);
  }

  // Returns a KmerSet, ignoring ones that appear not often.
  // The number of ignored k-mers is also returned.
  std::pair<KmerSet<K, N, KeyType>, std::int64_t> ToKmerSet(
      ValueType cutoff, int n_workers) const {
    KmerSet<K, N, KeyType> set;

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

          for (KeyType key : buf) {
            set.Add(GetKmerFromBucketAndKey<K, N, KeyType>(bucket_id, key));
          }
        },
        n_workers);

    return std::make_pair(set, static_cast<std::int64_t>(cutoff_count));
  }

  // Returns the count for a kmer.
  ValueType Get(const Kmer<K>& kmer) const {
    int bucket;
    KeyType key;
    std::tie(bucket, key) = GetBucketAndKeyFromKmer<K, N, KeyType>(kmer);

    auto it = buckets_[bucket].find(key);
    if (it == buckets_[bucket].end()) return 0;
    return it->second;
  }

  // Increments the count for a kmer by "v".
  KmerCounter& Add(const Kmer<K>& kmer, ValueType v) {
    int bucket;
    KeyType key;
    std::tie(bucket, key) = GetBucketAndKeyFromKmer<K, N, KeyType>(kmer);

    buckets_[bucket][key] = AddWithMax(buckets_[bucket][key], v);
    return *this;
  }

 private:
  static constexpr int kBucketsNum = 1 << N;

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
};

#endif

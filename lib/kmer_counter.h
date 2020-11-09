#ifndef KMER_COUNTER_H_
#define KMER_COUNTER_H_

#include <array>
#include <atomic>
#include <mutex>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/status/status.h"
#include "absl/strings/str_split.h"
#include "kmer.h"
#include "kmer_set.h"
#include "range.h"

template <int K, typename KeyType>
class KmerCounter {
 public:
  // Returns the number of distinct k-mers.
  int64_t Size() const {
    int64_t sum = 0;
    for (const auto& bucket : buckets_) {
      sum += bucket.size();
    }
    return sum;
  }

  absl::Status FromFASTQ(std::istream& is, bool canonical) {
    std::vector<std::string> lines;

    {
      const absl::Status error =
          absl::UnknownError("failed to parse FASTQ file");

      std::string s;

      while (true) {
        // EOF
        if (!std::getline(is, s)) break;

        if (s.length() == 0 || s[0] != '@') return error;

        if (!std::getline(is, s)) return error;

        lines.push_back(s);

        if (!std::getline(is, s)) return error;
        if (!std::getline(is, s)) return error;
      }
    }

    std::array<std::mutex, kBucketsNum> buckets_mutexes;

    std::vector<std::thread> threads;

    for (const Range& range :
         Range(0, lines.size()).Split(std::thread::hardware_concurrency())) {
      threads.emplace_back([&, range] {
        std::array<Bucket, kBucketsNum> buf;

        for (int64_t i = range.begin; i < range.end; i++) {
          const auto& line = lines[i];
          std::vector<std::string> fragments = absl::StrSplit(line, "N");
          for (const auto& fragment : fragments) {
            for (int i = 0; i + K <= (int)fragment.length(); i++) {
              const Kmer<K> kmer(fragment.substr(i, K));
              const auto [bucket, key] = GetBucketAndKeyFromKmer<K, KeyType>(
                  canonical ? kmer.Canonical() : kmer);
              buf[bucket][key] += 1;
            }
          }
        }

        std::vector<bool> done(kBucketsNum);
        int done_count = 0;

        while (done_count < kBucketsNum) {
          for (int i = 0; i < kBucketsNum; i++) {
            if (done[i]) continue;

            auto& mu = buckets_mutexes[i];

            if (mu.try_lock()) {
              auto& bucket = buckets_[i];

              bucket.reserve(bucket.size() + buf[i].size());
              for (const auto& [key, count] : buf[i]) {
                bucket[key] += count;
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

    return absl::OkStatus();
  }

  // Returns a KmerSet, ignoring ones that appear less often.
  // The number of ignored k-mers is also returned.
  std::pair<KmerSet<K, KeyType>, int64_t> Set(int cutoff) const {
    KmerSet<K, KeyType> set;
    std::mutex mu;

    std::vector<std::thread> threads;
    std::atomic_int64_t cutoff_count = 0;

    ForEachBucket([&](const Bucket& bucket, int bucketID) {
      std::vector<KeyType> buf;

      for (const auto& [key, count] : bucket) {
        if (count < cutoff) {
          cutoff_count += 1;
          continue;
        }

        buf.push_back(key);
      }

      {
        std::lock_guard _{mu};
        set.buckets_[bucketID].reserve(set.buckets_[bucketID].size() +
                                       buf.size());
        for (KeyType key : buf) {
          set.buckets_[bucketID].insert(key);
        }
      }
    });

    return {set, (int64_t)cutoff_count};
  }

 private:
  static constexpr int kBucketsNumFactor = 2 * K - sizeof(KeyType) * 8;
  static constexpr int kBucketsNum = 1 << (2 * K - sizeof(KeyType) * 8);

  using Bucket = absl::flat_hash_map<KeyType, int>;

  std::array<Bucket, kBucketsNum> buckets_;

  // Dispatches processing for each bucket.
  // Usage:
  //   ForEachBucket([&](const Bucket& bucket, int bucketID){ ... });
  template <typename F>
  void ForEachBucket(F&& f) const {
    std::vector<std::thread> threads;

    for (const Range& range :
         Range(0, kBucketsNum).Split(std::thread::hardware_concurrency())) {
      threads.emplace_back([&, range] {
        for (int i = range.begin; i < range.end; i++) {
          const Bucket& bucket = buckets_[i];
          f(bucket, i);
        }
      });
    }

    for (auto& thread : threads) thread.join();
  }
};

#endif

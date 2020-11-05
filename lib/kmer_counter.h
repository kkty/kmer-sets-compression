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
#include "absl/strings/str_split.h"
#include "kmer.h"
#include "kmer_set.h"
#include "range.h"

// Use (1 << B) buckets internally.
template <int K, int B>
class KmerCounter {
 public:
  // Find k-mers in a FASTQ file.
  void FromFASTQ(std::istream& is, bool canonical = false) {
    std::vector<std::string> lines;

    {
      std::array<std::string, 4> buf;

      // Reads 4 lines from "is" and pushes to "buf".
      auto consume = [&is, &buf]() {
        for (int i = 0; i < 4; i++) {
          if (!std::getline(is, buf[i])) return false;
        }

        return true;
      };

      while (consume() && buf[0].length() && buf[0][0] == '@') {
        lines.push_back(std::move(buf[1]));
      }
    }

    std::array<std::mutex, 1 << B> buckets_mutexes;

    std::vector<std::thread> threads;

    const auto ranges =
        SplitRange(0, lines.size(), std::thread::hardware_concurrency());

    for (const auto& range : ranges) {
      threads.emplace_back([&, range] {
        std::array<Bucket, 1 << B> buf;

        for (int64_t i = range.first; i < range.second; i++) {
          const auto& line = lines[i];
          std::vector<std::string> fragments = absl::StrSplit(line, "N");
          for (const auto& fragment : fragments) {
            for (int i = 0; i + K <= (int)fragment.length(); i++) {
              const Kmer<K> kmer(fragment.substr(i, K));
              const auto [bucket, key] = GetBucketAndKeyFromKmer<K, B>(
                  canonical ? kmer.Canonical() : kmer);
              buf[bucket][key] += 1;
            }
          }
        }

        std::vector<int> v;
        v.reserve(1 << B);
        for (int i = 0; i < (1 << B); i++) v.push_back(i);

        {
          std::default_random_engine e(range.first);
          std::shuffle(v.begin(), v.end(), e);
        }

        for (int i : v) {
          std::lock_guard _{buckets_mutexes[i]};
          auto& bucket = buckets_[i];
          bucket.reserve(bucket.size() + buf[i].size());
          for (const auto& [key, count] : buf[i]) {
            bucket[key] += count;
          }
        }
      });
    }

    for (std::thread& thread : threads) thread.join();
  }

  int Size() const {
    int sum = 0;
    for (const auto& bucket : buckets_) {
      sum += bucket.size();
    }
    return sum;
  }

  std::pair<KmerSet<K, B>, int64_t> Set(int cutoff) const {
    KmerSet<K, B> set;
    std::vector<std::thread> threads;
    std::atomic_int64_t cutoff_count = 0;

    for (const auto& range :
         SplitRange(0, 1 << B, std::thread::hardware_concurrency())) {
      threads.emplace_back([&, range] {
        for (int i = range.first; i < range.second; i++) {
          const Bucket& bucket = buckets_[i];
          for (const auto& [key, count] : bucket) {
            if (count < cutoff) {
              cutoff_count += 1;
              continue;
            }

            set.Add(i, key);
          }
        }
      });
    }

    for (std::thread& thread : threads) thread.join();

    return {set, (int64_t)cutoff_count};
  }

 private:
  using Bucket = absl::flat_hash_map<std::bitset<K * 2 - B>, int>;

  std::array<Bucket, 1 << B> buckets_;
};

#endif

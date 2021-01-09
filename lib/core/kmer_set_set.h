#ifndef CORE_KMER_SET_SET_H_
#define CORE_KMER_SET_SET_H_

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <mutex>
#include <optional>
#include <queue>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/random/random.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "core/io.h"
#include "core/kmer.h"
#include "core/kmer_set.h"
#include "core/kmer_set_compact.h"
#include "core/range.h"
#include "core/spss.h"
#include "spdlog/spdlog.h"

namespace internal {

using AdjacencyList = absl::flat_hash_map<int, std::vector<int>>;

// Serializes an AdjacencyList to a string. The string only contains digits and
// whitespaces.
std::string SerializeAdjacencyList(const AdjacencyList& adjacency_list) {
  std::stringstream ss;
  ss << adjacency_list.size();

  for (const std::pair<const int, std::vector<int>>& p : adjacency_list) {
    ss << ' ' << p.first;
    ss << ' ' << p.second.size();
    for (int i : p.second) ss << ' ' << i;
  }

  return ss.str();
}

// Deserializes an AdjacencyList from a string which was produced by
// SerializeAdjacencyList().
AdjacencyList DeserializeAdjacencyList(const std::string& s) {
  std::stringstream ss(s);

  AdjacencyList adjacency_list;

  std::size_t size;
  ss >> size;
  adjacency_list.reserve(size);

  for (std::size_t i = 0; i < size; i++) {
    int key;
    ss >> key;

    std::size_t value_size;
    ss >> value_size;

    std::vector<int> value(value_size);
    for (std::size_t j = 0; j < value_size; j++) {
      ss >> value[j];
    }

    adjacency_list[key] = std::move(value);
  }

  return adjacency_list;
}

}  // namespace internal

// KmerSetSet can be used to represent multiple kmer sets efficiently.
template <int K, int N, typename KeyType>
class KmerSetSet {
 public:
  KmerSetSet() = default;

  KmerSetSet(std::vector<KmerSetCompact<K, N, KeyType>> kmer_sets_compact,
             bool canonical, int n_workers)
      : kmer_sets_compact_(std::move(kmer_sets_compact)) {
    // Considers a non-directed complete graph where ith node represents
    // kmer_sets_compact_[i].
    // If the edge between ith node and jth node has large weight,
    // kmer_sets_compact_[i] and kmer_sets_compact_[j] have many common kmers.

    // To calculate weights among nodes efficiently, we use "SampledKmerSet"s.
    using SampledKmerSet = std::vector<std::vector<KeyType>>;

    // "bucket_ids" contains random numbers from 0 to (1 << N) - 1.
    // Elements in the ith element of a SampledKmerSet corresponds to
    // bucket_ids[i].
    std::vector<int> bucket_ids;

    {
      absl::InsecureBitGen bitgen;

      for (int i = 0; i < (1 << N); i++) {
        // Adds i to bucket_ids with the probability of 5%.
        if (absl::Uniform(bitgen, 0, 20) == 0) {
          bucket_ids.push_back(i);
        }
      }
    }

    int n_buckets = bucket_ids.size();

    // sampled_kmer_sets[i] holds a SampledKmerSet made from
    // kmer_sets_compact_[i].
    std::vector<SampledKmerSet> sampled_kmer_sets;

    spdlog::debug("constructing initial sampled_kmer_sets");

    {
      const int n = kmer_sets_compact_.size();

      sampled_kmer_sets.resize(n);

      boost::asio::thread_pool pool(n_workers);

      for (int i = 0; i < n; i++) {
        boost::asio::post(pool, [&, i] {
          sampled_kmer_sets[i] =
              kmer_sets_compact_[i].GetSampledKmerSet(bucket_ids, canonical, 1);
        });
      }

      pool.join();
    }

    spdlog::debug("constructed initial sampled_kmer_sets");

    // Calculates the weight of the edge between ith node and jth node.
    const auto GetEdgeWeight = [&](int i, int j) {
      std::int64_t count = 0;

      for (int k = 0; k < n_buckets; k++) {
        const std::vector<KeyType>& bucket_i = sampled_kmer_sets[i][k];
        const std::vector<KeyType>& bucket_j = sampled_kmer_sets[j][k];

        // Calculates the intersection size of bucket_i and bucket_j.

        auto it_i = bucket_i.begin();
        auto it_j = bucket_j.begin();

        while (it_i != bucket_i.end() && it_j != bucket_j.end()) {
          if (*it_i < *it_j) {
            it_i++;
          } else if (*it_i > *it_j) {
            it_j++;
          } else {
            count += 1;
            it_i += 1;
            it_j += 1;
          }
        }
      }

      return count;
    };

    // weights[{i, j}] (i < j) is the weight between the ith node and jth node.
    absl::flat_hash_map<std::pair<int, int>, std::int64_t> weights;

    spdlog::debug("calculating initial weights");

    {
      const int n = kmer_sets_compact_.size();

      boost::asio::thread_pool pool(n_workers);
      std::mutex mu;

      std::vector<std::pair<int, int>> pairs;

      for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
          pairs.emplace_back(i, j);
        }
      }

      for (const Range& range :
           Range(0, pairs.size()).Split(n_workers * n_workers)) {
        boost::asio::post(pool, [&, range] {
          for (int i : range) {
            int j, k;
            std::tie(j, k) = pairs[i];
            std::int64_t weight = GetEdgeWeight(j, k);
            std::lock_guard lck(mu);
            weights[std::make_pair(j, k)] = weight;
          }
        });
      }

      pool.join();
    }

    spdlog::debug("calculated initial weights");

    // For each i, calculates the size of kmer_sets_compact_[i] and sums them
    // up. This value will be updated in the main loop accordingly.
    std::atomic_int64_t total_size = 0;

    {
      boost::asio::thread_pool pool(n_workers);
      for (std::size_t i = 0; i < kmer_sets_compact_.size(); i++) {
        boost::asio::post(
            pool, [&, i] { total_size += kmer_sets_compact_[i].Size(); });
      }
      pool.join();
    }

    spdlog::debug("total_size = {}", total_size);

    // Returns an approximate of the final file size.
    const auto GetExpectedFinalFileSize = [&] {
      std::atomic_int64_t size = 0;

      boost::asio::thread_pool pool(n_workers);

      for (int i = 0; i < static_cast<int>(kmer_sets_compact_.size()); i++) {
        boost::asio::post(pool, [&, i] {
          if (canonical) {
            size += GetSPSSWeightCanonical(
                kmer_sets_compact_[i].ToKmerSet(canonical, 1), 1, 0.1);
          } else {
            size += GetSPSSWeight(kmer_sets_compact_[i].ToKmerSet(canonical, 1),
                                  1, 0.1);
          }
        });
      }

      pool.join();

      return static_cast<std::int64_t>(size);
    };

    spdlog::debug("calculating expected_final_file_size");

    std::int64_t expected_final_file_size = GetExpectedFinalFileSize();

    spdlog::debug("calculated expected_final_file_size");
    spdlog::debug("expected_final_file_size = {}", expected_final_file_size);

    // The interval with which expected_final_file_size is updated.
    const int interval = kmer_sets_compact_.size() / 4 + 1;

    // The threshold to decide when to stop the iterations.
    // If the (ratio of) decrease of expected_final_file_size during an interval
    // is less than this number, we break the loop.
    const float improvement_threshold =
        0.1 * interval / kmer_sets_compact_.size();

    spdlog::debug("interval = {}, improvement_threshold = {}", interval,
                  improvement_threshold);

    for (int i = 0;; i++) {
      spdlog::debug("i = {}", i);

      if (i > 0 && i % interval == 0) {
        // Updates expected_final_file_size.
        // If the improvement is less than the threshold, breaks the loop.

        spdlog::debug("updating expected_final_file_size");

        const std::int64_t updated = GetExpectedFinalFileSize();

        const float improvement =
            static_cast<float>(expected_final_file_size - updated) /
            expected_final_file_size;

        spdlog::debug(
            "expected_final_file_size = {}, updated = {}, improvement = {}",
            expected_final_file_size, updated, improvement);

        if (improvement <= improvement_threshold) {
          break;
        }

        expected_final_file_size = updated;

        spdlog::debug("updated expected_final_file_size");
      }

      const int n = kmer_sets_compact_.size();

      // Finds i, j such that the weight between ith node and jth node is
      // maximized.
      std::int64_t weight = 0;
      int j, k;
      for (const std::pair<const std::pair<int, int>, std::int64_t>& p :
           weights) {
        if (p.second > weight) {
          std::tie(j, k) = p.first;
          weight = p.second;
        }
      }

      // If there are no edges with positive weights, stops the iteration.
      if (weight == 0) {
        spdlog::debug("found no positive weights");
        break;
      }

      spdlog::debug("j = {}, k = {}, weight = {}", j, k, weight);

      const std::int64_t original_size =
          kmer_sets_compact_[j].Size() + kmer_sets_compact_[k].Size();

      spdlog::debug(
          "updating kmer_sets_compact_, sampled_kmer_sets_, and children_");

      {
        KmerSet<K, N, KeyType> kmer_set_j =
            kmer_sets_compact_[j].ToKmerSet(canonical, n_workers);

        KmerSet<K, N, KeyType> kmer_set_k =
            kmer_sets_compact_[k].ToKmerSet(canonical, n_workers);

        {
          const KmerSet<K, N, KeyType> kmer_set_n =
              Intersection(kmer_set_j, kmer_set_k, n_workers);

          kmer_set_j.Sub(kmer_set_n, n_workers);
          kmer_set_k.Sub(kmer_set_n, n_workers);

          kmer_sets_compact_.push_back(
              KmerSetCompact<K, N, KeyType>::FromKmerSet(kmer_set_n, canonical,
                                                         true, n_workers));

          sampled_kmer_sets.push_back(kmer_sets_compact_[n].GetSampledKmerSet(
              bucket_ids, canonical, n_workers));
        }

        kmer_sets_compact_[j] = KmerSetCompact<K, N, KeyType>::FromKmerSet(
            kmer_set_j, canonical, true, n_workers);

        sampled_kmer_sets[j] = kmer_sets_compact_[j].GetSampledKmerSet(
            bucket_ids, canonical, n_workers);

        kmer_sets_compact_[k] = KmerSetCompact<K, N, KeyType>::FromKmerSet(
            kmer_set_k, canonical, true, n_workers);

        sampled_kmer_sets[k] = kmer_sets_compact_[k].GetSampledKmerSet(
            bucket_ids, canonical, n_workers);

        children_[j].push_back(n);
        children_[k].push_back(n);
      }

      spdlog::debug(
          "updated kmer_sets_compact_, sampled_kmer_sets, and children_");

      // The change in the number of kmers that are stored in
      // kmer_sets_compact_. It should be negative.
      const std::int64_t size_diff =
          kmer_sets_compact_[n].Size() + kmer_sets_compact_[j].Size() +
          kmer_sets_compact_[k].Size() - original_size;

      total_size += size_diff;
      spdlog::debug("size_diff = {}, total_size = {}", size_diff, total_size);

      // Updates the graph.
      // Weights of edges that are incident to jth node, kth node, or nth node
      // should be modified / added.
      {
        spdlog::debug("updating weights");

        std::vector<std::pair<int, int>> pairs;

        for (int l = 0; l < n; l++) {
          if (j == l) continue;

          pairs.emplace_back(std::min(j, l), std::max(j, l));
        }

        for (int l = 0; l < n; l++) {
          if (k == l) continue;

          pairs.emplace_back(std::min(k, l), std::max(k, l));
        }

        for (int l = 0; l < n; l++) {
          pairs.emplace_back(l, n);
        }

        boost::asio::thread_pool pool(n_workers);
        std::mutex mu;

        for (const Range& range :
             Range(0, pairs.size()).Split(n_workers * n_workers)) {
          boost::asio::post(pool, [&, range] {
            for (int j : range) {
              int k, l;
              std::tie(k, l) = pairs[j];
              std::int64_t weight = GetEdgeWeight(k, l);
              std::lock_guard lck(mu);
              weights[std::make_pair(k, l)] = weight;
            }
          });
        }

        pool.join();

        spdlog::debug("updated weights");
      }
    }
  }

  // Returns the number of kmer sets in the structure.
  int Size() const { return kmer_sets_compact_.size(); }

  // Reconstructs the ith kmer set.
  KmerSet<K, N, KeyType> Get(int i, bool canonical, int n_workers) const {
    KmerSet<K, N, KeyType> kmer_set;

    std::queue<int> queue;
    queue.push(i);

    while (!queue.empty()) {
      int current = queue.front();
      queue.pop();

      kmer_set.Add(kmer_sets_compact_[current].ToKmerSet(canonical, n_workers),
                   n_workers);
      auto it = children_.find(current);
      if (it != children_.end()) {
        for (int child : it->second) {
          queue.push(child);
        }
      }
    }

    return kmer_set;
  }

  // Dumps data to files in a directory.
  // Example:
  //   ... = Dump("./dir", "bzip2", "bz2", n_workers);
  absl::Status Dump(const std::string& directory_name,
                    const std::string& compressor, const std::string& extension,
                    int n_workers) {
    try {
      std::filesystem::create_directories(directory_name);
    } catch (...) {
      return absl::InternalError("failed to create a directory");
    }

    // Dumps children_ and kmer_sets_compact_.size().
    {
      std::vector<std::string> v;

      v.push_back(internal::SerializeAdjacencyList(children_));

      {
        std::stringstream ss;
        ss << kmer_sets_compact_.size();
        v.push_back(ss.str());
      }

      const std::string file_name =
          (std::filesystem::path(directory_name) /
           std::filesystem::path(absl::StrFormat("meta.%s", extension)))
              .string();

      absl::Status status = WriteLines(file_name, compressor, v);

      if (!status.ok()) {
        return status;
      }
    }

    // Dumps kmer_sets_compact_.

    boost::asio::thread_pool pool(n_workers);
    std::atomic_int fail_count = 0;

    for (std::size_t i = 0; i < kmer_sets_compact_.size(); i++) {
      boost::asio::post(pool, [&, i] {
        spdlog::debug("dumping kmer set: i = {}", i);

        const KmerSetCompact<K, N, KeyType>& kmer_set_compact =
            kmer_sets_compact_[i];

        const std::string file_name = absl::StrFormat("%d.%s", i, extension);

        const absl::Status status =
            kmer_set_compact.Dump((std::filesystem::path(directory_name) /
                                   std::filesystem::path(file_name))
                                      .string(),
                                  compressor, 1);

        if (!status.ok()) {
          spdlog::error("failed to write data to a file: {}",
                        status.ToString());
          fail_count += 1;
        }

        spdlog::debug("dumped kmer set: i = {}", i);
      });
    }

    pool.join();

    if (fail_count > 0) {
      return absl::InternalError(
          absl::StrFormat("failed to write %d files", fail_count));
    }

    return absl::OkStatus();
  }

  // Dumps the graph structure in DOT format.
  absl::Status DumpGraph(const std::string& file_name) const {
    std::vector<std::string> lines;

    lines.emplace_back("digraph G {");

    for (const std::pair<const int, std::vector<int>>& p : children_) {
      for (int i : p.second) {
        lines.push_back(absl::StrFormat("v%d -> v%d", p.first, i));
      }
    }

    lines.emplace_back("}");

    return WriteLines(file_name, "", lines);
  }

  // Loads data from files in a directory.
  static absl::StatusOr<KmerSetSet> Load(const std::string& directory_name,
                                         const std::string& decompressor,
                                         const std::string& extension,
                                         int n_workers) {
    absl::flat_hash_map<int, std::vector<int>> children;
    int n;  // Number of kmer sets.

    // Loads meta data.
    {
      const std::string file_name =
          (std::filesystem::path(directory_name) /
           std::filesystem::path(absl::StrFormat("meta.%s", extension)))
              .string();

      absl::StatusOr<std::vector<std::string>> statusor =
          ReadLines(file_name, decompressor);

      if (!statusor.ok()) {
        return statusor.status();
      }

      std::vector<std::string> lines = std::move(statusor).value();

      children = internal::DeserializeAdjacencyList(lines[0]);

      {
        std::stringstream ss(lines[1]);
        ss >> n;
      }
    }

    std::vector<KmerSetCompact<K, N, KeyType>> kmer_sets_compact(n);

    boost::asio::thread_pool pool(n_workers);
    std::atomic_int n_fail = 0;

    for (int i = 0; i < n; i++) {
      boost::asio::post(pool, [&, i] {
        const std::string file_name =
            (std::filesystem::path(directory_name) /
             std::filesystem::path(absl::StrFormat("%d.%s", i, extension)))
                .string();

        absl::StatusOr<KmerSetCompact<K, N, KeyType>> statusor =
            KmerSetCompact<K, N, KeyType>::Load(file_name, decompressor,
                                                n_workers);

        if (!statusor.ok()) {
          spdlog::error("failed to load {}: {}", file_name,
                        statusor.status().ToString());
          n_fail += 1;
          return;
        }

        kmer_sets_compact[i] = std::move(statusor).value();
      });
    }

    pool.join();

    if (n_fail > 0) {
      return absl::InternalError(
          absl::StrFormat("failed to dump %d files", n_fail));
    }

    return KmerSetSet(std::move(children), std::move(kmer_sets_compact));
  }

 private:
  KmerSetSet(absl::flat_hash_map<int, std::vector<int>> children,
             std::vector<KmerSetCompact<K, N, KeyType>> kmer_sets_compact)
      : children_(std::move(children)),
        kmer_sets_compact_(std::move(kmer_sets_compact)) {}

  absl::flat_hash_map<int, std::vector<int>> children_;
  std::vector<KmerSetCompact<K, N, KeyType>> kmer_sets_compact_;
};

// KmerSetSetReader can be used to reconstruct kmer sets from files in a
// directory.
template <int K, int N, typename KeyType>
class KmerSetSetReader {
 public:
  KmerSetSetReader() = default;

  static absl::StatusOr<KmerSetSetReader> FromDirectory(
      std::string directory_name, std::string extension,
      std::string decompressor, bool canonical) {
    std::vector<std::string> v;

    {
      const std::string file_name =
          (std::filesystem::path(directory_name) /
           std::filesystem::path(absl::StrFormat("meta.%s", extension)))
              .string();

      absl::StatusOr<std::vector<std::string>> statusor =
          ReadLines(file_name, decompressor);

      if (!statusor.ok()) {
        return statusor.status();
      }

      v = std::move(statusor).value();
    }

    absl::flat_hash_map<int, std::vector<int>> children =
        internal::DeserializeAdjacencyList(v[0]);

    int size;
    {
      std::stringstream ss(v[1]);
      ss >> size;
    }

    return KmerSetSetReader(directory_name, extension, decompressor, canonical,
                            children, size);
  }

  // Returns the number of managed kmer sets.
  int Size() const { return size_; }

  // Loads data from files and reconstructs the ith kmer set.
  absl::StatusOr<KmerSet<K, N, KeyType>> Get(int i, int n_workers) const {
    std::vector<int> ids;

    // BFS from i.
    {
      std::queue<int> queue;
      queue.push(i);

      while (!queue.empty()) {
        int current = queue.front();
        queue.pop();

        ids.push_back(current);

        auto it = children_.find(current);

        if (it == children_.end()) continue;

        for (int child : it->second) {
          queue.push(child);
        }
      }
    }

    KmerSet<K, N, KeyType> to_return;

    {
      boost::asio::thread_pool pool(n_workers);
      int n_done = 0;
      int n_fail = 0;
      std::mutex mu;

      for (int id : ids) {
        boost::asio::post(pool, [&, id] {
          spdlog::debug("loading kmer_set: id = {}", id);

          const std::string file_name =
              (std::filesystem::path(directory_name_) /
               std::filesystem::path(absl::StrFormat("%d.%s", id, extension_)))
                  .string();

          spdlog::debug("file_name = {}", file_name);

          absl::StatusOr<KmerSetCompact<K, N, KeyType>> statusor =
              KmerSetCompact<K, N, KeyType>::Load(file_name, decompressor_, 1);

          if (!statusor.ok()) {
            spdlog::debug("failed to load kmer_set: {}",
                          statusor.status().ToString());

            std::lock_guard lck(mu);

            n_fail += 1;
            return;
          }

          KmerSet<K, N, KeyType> kmer_set =
              std::move(statusor).value().ToKmerSet(canonical_, 1);

          spdlog::debug("loaded kmer_set: id = {}", id);

          std::lock_guard lck(mu);

          spdlog::debug("adding kmer_set: id = {}", id);

          to_return.Add(kmer_set, 1 + n_done);

          spdlog::debug("added kmer_set: id = {}", id);

          n_done += 1;
        });
      }

      pool.join();

      if (n_fail > 0) {
        const std::string message =
            absl::StrFormat("failed to load data from %d files", n_fail);
        return absl::InternalError(message);
      }
    }

    return to_return;
  }

 private:
  KmerSetSetReader(std::string directory_name, std::string extension,
                   std::string decompressor, bool canonical,
                   absl::flat_hash_map<int, std::vector<int>> children,
                   int size)
      : directory_name_(std::move(directory_name)),
        extension_(std::move(extension)),
        decompressor_(std::move(decompressor)),
        canonical_(canonical),
        children_(std::move(children)),
        size_(size) {}

  std::string directory_name_;
  std::string extension_;
  std::string decompressor_;
  bool canonical_;
  absl::flat_hash_map<int, std::vector<int>> children_;
  int size_ = 0;
};

#endif

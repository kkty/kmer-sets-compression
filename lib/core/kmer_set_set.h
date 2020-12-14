#ifndef CORE_KMER_SET_SET_H_
#define CORE_KMER_SET_SET_H_

#include <algorithm>
#include <atomic>
#include <cmath>
#include <mutex>
#include <optional>
#include <queue>
#include <sstream>
#include <thread>
#include <tuple>
#include <variant>
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
#include "boost/filesystem.hpp"
#include "core/io.h"
#include "core/kmer.h"
#include "core/kmer_set.h"
#include "core/kmer_set_compact.h"
#include "core/kmer_set_immutable.h"
#include "core/range.h"
#include "core/unitigs.h"
#include "spdlog/spdlog.h"

// KmerSetSet can be used to represent multiple kmer sets efficiently.
template <int K, int N, typename KeyType>
class KmerSetSet {
 public:
  KmerSetSet() = default;

  KmerSetSet(std::vector<KmerSetImmutable<K, N, KeyType>> kmer_sets_immutable,
             int n_iterations, int n_workers)
      : kmer_sets_immutable_(std::move(kmer_sets_immutable)) {
    // Considers a non-directed complete graph where ith node represents
    // kmer_sets_immutable_[i].

    // Calculates the weight of the edge between ith node and jth node.
    const auto GetEdgeWeight = [&](int i, int j) {
      return kmer_sets_immutable_[i].IntersectionSizeEstimate(
          kmer_sets_immutable_[j], 0.1);
    };

    // weights[{i, j}] (i < j) is the weight between the ith node and jth node.
    absl::flat_hash_map<std::pair<int, int>, int64_t> weights;

    spdlog::debug("calculating initial weights");
    {
      const int n = kmer_sets_immutable_.size();
      boost::asio::thread_pool pool(n_workers);
      std::mutex mu;

      for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
          boost::asio::post(pool, [&, i, j] {
            spdlog::debug("calculating weight: i = {}, j = {}", i, j);
            int64_t weight = GetEdgeWeight(i, j);
            spdlog::debug("calculated weight: i = {}, j = {}", i, j);

            std::lock_guard lck(mu);
            weights[{i, j}] = weight;
          });
        }
      }

      pool.join();
    }
    spdlog::debug("calculated initial weights");

    // For each i, calculates the size of kmer_sets_immutable_[i] and sums them
    // up. This value will be updated in the main loop accordingly.
    std::atomic_int64_t total_size = 0;
    {
      boost::asio::thread_pool pool(n_workers);
      for (size_t i = 0; i < kmer_sets_immutable_.size(); i++) {
        boost::asio::post(
            pool, [&, i] { total_size += kmer_sets_immutable_[i].Size(); });
      }
      pool.join();
    }
    spdlog::debug("total_size = {}", total_size);

    while (n_iterations--) {
      spdlog::debug("n_iterations = {}", n_iterations);

      const int n = kmer_sets_immutable_.size();

      // Finds i, j such that the weight between ith node and jth node is
      // maximized.
      int64_t weight = 0;
      int i, j;
      for (const std::pair<const std::pair<int, int>, int64_t>& p : weights) {
        if (p.second > weight) {
          std::tie(i, j) = p.first;
          weight = p.second;
        }
      }

      // If there are no edges with positive weights, stops the iteration.
      if (weight == 0) break;

      spdlog::debug("i = {}, j = {}, weight = {}", i, j, weight);

      // Constructs kmer_set by intersecting kmer_sets_immutable_[i] and
      // kmer_sets_immutable_[j].
      KmerSetImmutable<K, N, KeyType> kmer_set_immutable =
          kmer_sets_immutable_[i].Intersection(kmer_sets_immutable_[j],
                                               n_workers);

      const int64_t original_size =
          kmer_sets_immutable_[i].Size() + kmer_sets_immutable_[j].Size();

      // Moves kmer_set to kmer_sets_immutable_[n].
      kmer_sets_immutable_.push_back(std::move(kmer_set_immutable));

      // Updates kmer_sets_immutable_[i] and kmer_sets_immutable_[j] and adds
      // metadata so that the original kmer_sets_immutable_[i] and
      // kmer_sets_immutable_[j] can be reconstructed.
      kmer_sets_immutable_[i] =
          kmer_sets_immutable_[i].Sub(kmer_sets_immutable_[n], n_workers);
      kmer_sets_immutable_[j] =
          kmer_sets_immutable_[j].Sub(kmer_sets_immutable_[n], n_workers);

      children_[i].push_back(n);
      children_[j].push_back(n);

      // The change in the number of kmers that are stored in
      // kmer_sets_immutable_. It should be negative.
      const int64_t size_diff = kmer_sets_immutable_[n].Size() +
                                kmer_sets_immutable_[i].Size() +
                                kmer_sets_immutable_[j].Size() - original_size;

      total_size += size_diff;
      spdlog::debug("size_diff = {}, total_size = {}", size_diff, total_size);

      // Updates the graph.
      // Weights of edges that are incident to ith node, jth node, or nth node
      // should be modified / added.
      {
        spdlog::debug("updating weights");

        boost::asio::thread_pool pool(n_workers);
        std::mutex mu;

        for (int k = 0; k < n; k++) {
          if (i == k) continue;

          boost::asio::post(pool, [&, k] {
            int64_t weight = GetEdgeWeight(i, k);
            std::lock_guard lck(mu);
            weights[{std::min(i, k), std::max(i, k)}] = weight;
          });
        }

        for (int k = 0; k < n; k++) {
          if (j == k) continue;

          boost::asio::post(pool, [&, k] {
            int64_t weight = GetEdgeWeight(j, k);
            std::lock_guard lck(mu);
            weights[{std::min(j, k), std::max(j, k)}] = weight;
          });
        }

        for (int k = 0; k < n; k++) {
          boost::asio::post(pool, [&, k] {
            int64_t weight = GetEdgeWeight(k, n);
            std::lock_guard lck(mu);
            weights[{k, n}] = weight;
          });
        }

        pool.join();

        spdlog::debug("updated weights");
      }
    }
  }

  // Returns the number of kmer sets in the structure.
  int Size() const { return kmer_sets_immutable_.size(); }

  // Reconstructs the ith kmer set.
  KmerSet<K, N, KeyType> Get(int i, int n_workers) const {
    KmerSet<K, N, KeyType> kmer_set;

    std::queue<int> queue;
    queue.push(i);

    while (!queue.empty()) {
      int current = queue.front();
      queue.pop();

      kmer_set.Add(kmer_sets_immutable_[current].ToKmerSet(n_workers),
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

  // Dumps the structure to a vector of strings.
  // If "clear" is true, the contents of the structure will be invalidated.
  std::vector<std::string> Dump(bool canonical, bool clear, int n_workers) {
    std::vector<std::string> v;

    // Dumps "children_".

    {
      std::stringstream ss;
      ss << children_.size();
      v.push_back(ss.str());
    }

    for (const std::pair<const int, std::vector<int>>& p : children_) {
      std::stringstream ss;
      ss << p.first << ' ' << p.second.size();
      for (int child : p.second) {
        ss << ' ' << child;
      }
      v.push_back(ss.str());
    }

    if (clear) {
      absl::flat_hash_map<int, std::vector<int>>().swap(children_);
    }

    // Dumps "kmer_sets_immutable_".

    {
      std::stringstream ss;
      ss << kmer_sets_immutable_.size();
      v.push_back(ss.str());
    }

    const int n = v.size();

    v.resize(n + kmer_sets_immutable_.size());

    {
      boost::asio::thread_pool pool(n_workers);

      for (size_t i = 0; i < kmer_sets_immutable_.size(); i++) {
        boost::asio::post(pool, [&, i] {
          spdlog::debug("dumping kmer set: i = {}", i);

          const KmerSetCompact<K, N, KeyType> kmer_set_compact =
              KmerSetCompact<K, N, KeyType>::FromKmerSet(
                  kmer_sets_immutable_[i].ToKmerSet(1), canonical, 1);

          v[n + i] = absl::StrJoin(kmer_set_compact.Dump(), " ");

          if (clear) {
            kmer_sets_immutable_[i].Clear();
          }

          spdlog::debug("dumped kmer set: i = {}", i);
        });
      }

      pool.join();
    }

    return v;
  }

  // Dumps data to files in a folder.
  absl::Status DumpToFolder(const std::string& folder_name,
                            const std::string& compressor,
                            const std::string& extension, bool canonical,
                            bool clear, int n_workers) {
    if (!boost::filesystem::exists(folder_name)) {
      boost::filesystem::create_directories(folder_name);
    }

    {
      std::vector<std::string> v;

      {
        std::stringstream ss;
        ss << children_.size();
        v.push_back(ss.str());
      }

      for (const std::pair<const int, std::vector<int>>& p : children_) {
        std::stringstream ss;
        ss << p.first << ' ' << p.second.size();
        for (int child : p.second) {
          ss << ' ' << child;
        }
        v.push_back(ss.str());
      }

      {
        std::stringstream ss;
        ss << kmer_sets_immutable_.size();
        v.push_back(ss.str());
      }

      const std::string file_name =
          (boost::filesystem::path(folder_name) /
           boost::filesystem::path(absl::StrFormat("meta.%s", extension)))
              .string();

      absl::Status status = WriteLines(file_name, compressor, v);

      if (!status.ok()) {
        return status;
      }
    }

    if (clear) {
      absl::flat_hash_map<int, std::vector<int>>().swap(children_);
    }

    boost::asio::thread_pool pool(n_workers);
    std::atomic_int fail_count = 0;

    for (size_t i = 0; i < kmer_sets_immutable_.size(); i++) {
      boost::asio::post(pool, [&, i] {
        spdlog::debug("dumping kmer set: i = {}", i);

        const KmerSetCompact<K, N, KeyType> kmer_set_compact =
            KmerSetCompact<K, N, KeyType>::FromKmerSet(
                kmer_sets_immutable_[i].ToKmerSet(1), canonical, 1);
        if (clear) {
          kmer_sets_immutable_[i].Clear();
        }

        const std::string file_name = absl::StrFormat("%d.%s", i, extension);

        const absl::Status status =
            WriteLines((boost::filesystem::path(folder_name) /
                        boost::filesystem::path(file_name))
                           .string(),
                       compressor, kmer_set_compact.Dump());

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
          absl::StrFormat("failed to write %d data", fail_count));
    }

    return absl::OkStatus();
  }

  // Dumps data to a file.
  absl::Status Dump(const std::string& file_name, const std::string& compressor,
                    bool canonical, bool clear, int n_workers) {
    return WriteLines(file_name, compressor, Dump(canonical, clear, n_workers));
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

  // Loads data from a vector of strings.
  static KmerSetSet Load(std::vector<std::string> v, bool canonical,
                         int n_workers) {
    absl::flat_hash_map<int, std::vector<int>> children;
    std::vector<KmerSetImmutable<K, N, KeyType>> kmer_sets_immutable;

    // Loads "children".

    int children_size;
    {
      std::stringstream ss(v[0]);
      ss >> children_size;
    }
    children.reserve(children_size);

    for (int i = 0; i < children_size; i++) {
      std::stringstream ss(v[i + 1]);
      int key;
      int size;
      ss >> key >> size;
      for (int j = 0; j < size; j++) {
        int value;
        ss >> value;
        children[key].push_back(value);
      }
    }

    // Loads "kmer_sets".

    // We are concerned with {v[offset] ... v[v.size() - 1]}.
    const int offset = children_size + 1;

    int kmer_sets_immutable_size;
    {
      std::stringstream ss(v[offset]);
      ss >> kmer_sets_immutable_size;
    }

    kmer_sets_immutable.resize(kmer_sets_immutable_size);

    {
      boost::asio::thread_pool pool(n_workers);

      for (int i = 0; i < kmer_sets_immutable_size; i++) {
        boost::asio::post(pool, [&, i] {
          std::string& s = v[offset + 1 + i];
          KmerSetCompact<K, N, KeyType> kmer_set_compact =
              KmerSetCompact<K, N, KeyType>::Load(absl::StrSplit(s, ' '));
          kmer_sets_immutable[i] = KmerSetImmutable<K, N, KeyType>(
              kmer_set_compact.ToKmerSet(canonical, n_workers), 1);
          std::string().swap(s);
        });
      }

      pool.join();
    }

    return KmerSetSet(children, kmer_sets_immutable);
  }

  static absl::StatusOr<KmerSetSet> Load(const std::string& file_name,
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

    return Load(std::move(lines), canonical, n_workers);
  }

 private:
  KmerSetSet(absl::flat_hash_map<int, std::vector<int>> children,
             std::vector<KmerSetImmutable<K, N, KeyType>> kmer_sets_immutable)
      : children_(std::move(children)),
        kmer_sets_immutable_(std::move(kmer_sets_immutable)) {}

  absl::flat_hash_map<int, std::vector<int>> children_;
  std::vector<KmerSetImmutable<K, N, KeyType>> kmer_sets_immutable_;
};

#endif

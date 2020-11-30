#ifndef CORE_KMER_SET_SET_H_
#define CORE_KMER_SET_SET_H_

#include <algorithm>
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
#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "core/graph.h"
#include "core/io.h"
#include "core/kmer.h"
#include "core/kmer_set.h"
#include "core/kmer_set_compressed.h"
#include "core/range.h"
#include "core/unitigs.h"
#include "lemon/list_graph.h"
#include "lemon/matching.h"
#include "spdlog/spdlog.h"

// KmerSetSet can be used to represent multiple kmer sets efficiently.
template <int K, typename KeyType>
class KmerSetSet {
 public:
  KmerSetSet() = default;

  KmerSetSet(std::vector<KmerSet<K, KeyType>> kmer_sets, int n_iterations,
             int n_workers)
      : kmer_sets_(std::move(kmer_sets)) {
    // Considers a non-directed complete graph where ith node represents
    // kmer_sets_[i].

    // Calculates the weight of the edge between ith node and jth node.
    const auto GetEdgeWeight = [&](int i, int j) {
      return kmer_sets_[i].CommonEstimate(kmer_sets_[j], 0.1);
    };

    // weights[{i, j}] (i < j) is the weight between the ith node and jth node.
    absl::flat_hash_map<std::pair<int, int>, int64_t> weights;

    spdlog::debug("calculating initial weights");
    {
      const int n = kmer_sets_.size();
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

    // For each i, calculates the size of kmer_sets_[i] and sums them up.
    // This value will be updated in the main loop accordingly.
    std::atomic_int64_t total_size = 0;
    {
      boost::asio::thread_pool pool(n_workers);
      for (size_t i = 0; i < kmer_sets_.size(); i++) {
        boost::asio::post(pool, [&, i] { total_size += kmer_sets_[i].Size(); });
      }
      pool.join();
    }
    spdlog::debug("total_size = {}", total_size);

    while (n_iterations--) {
      spdlog::debug("n_iterations = {}", n_iterations);

      const int n = kmer_sets_.size();

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

      // Constructs kmer_set by intersecting kmer_sets_[i] and kmer_sets_[j].
      KmerSet<K, KeyType> kmer_set =
          Intersection(kmer_sets_[i], kmer_sets_[j], n_workers);

      const int64_t original_size = kmer_sets_[i].Size() + kmer_sets_[j].Size();

      // Moves kmer_set to kmer_sets_[n].
      kmer_sets_.push_back(std::move(kmer_set));

      // Updates kmer_sets_[i] and kmer_sets_[j] and adds metadata so that the
      // original kmer_sets_[i] and kmer_sets_[j] can be reconstructed.
      kmer_sets_[i].Sub(kmer_sets_[n], n_workers);
      kmer_sets_[j].Sub(kmer_sets_[n], n_workers);
      children_[i].push_back(n);
      children_[j].push_back(n);

      // The change in the number of kmers that are stored in kmer_sets_.
      // It should be negative.
      const int64_t size_diff = kmer_sets_[n].Size() + kmer_sets_[i].Size() +
                                kmer_sets_[j].Size() - original_size;

      total_size += size_diff;
      spdlog::debug("size_diff = {}, total_size = {}", size_diff, total_size);

      // Updates the graph.
      // Weights of edges that are incident to ith node, jth node, or nth node
      // should be modified / added.
      {
        spdlog::debug("updating weights");

        boost::asio::thread_pool pool(n_workers);
        std::mutex mu;

        for (int k = i + 1; k < n; k++) {
          boost::asio::post(pool, [&, k] {
            int64_t weight = GetEdgeWeight(i, k);
            std::lock_guard lck(mu);
            weights[{i, k}] = weight;
          });
        }

        for (int k = j + 1; k < n; k++) {
          boost::asio::post(pool, [&, k] {
            int64_t weight = GetEdgeWeight(j, k);
            std::lock_guard lck(mu);
            weights[{j, k}] = weight;
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

  // Reconstructs the ith kmer set.
  KmerSet<K, KeyType> Get(int i, int n_workers) const {
    KmerSet<K, KeyType> kmer_set;

    std::queue<int> queue;
    queue.push(i);

    while (!queue.empty()) {
      int i = queue.front();
      queue.pop();

      kmer_set.Add(kmer_sets_[i], n_workers);
      auto it = children_.find(i);
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

    // Dumps "kmer_sets_".

    {
      std::stringstream ss;
      ss << kmer_sets_.size();
      v.push_back(ss.str());
    }

    const int n = v.size();

    v.resize(n + kmer_sets_.size());

    {
      boost::asio::thread_pool pool(n_workers);

      for (size_t i = 0; i < kmer_sets_.size(); i++) {
        boost::asio::post(pool, [&, i] {
          spdlog::debug("dumping kmer set: i = {}", i);

          const KmerSetCompressed<K, KeyType> kmer_set_compressed =
              KmerSetCompressed<K, KeyType>::FromKmerSet(kmer_sets_[i],
                                                         canonical, 1);

          v[n + i] = absl::StrJoin(kmer_set_compressed.Dump(), " ");

          if (clear) {
            kmer_sets_[i].Clear();
          }

          spdlog::debug("dumped kmer set: i = {}", i);
        });
      }

      pool.join();
    }

    return v;
  }

  // Dumps data to a file.
  void Dump(const std::string& file_name, const std::string& compressor,
            bool canonical, bool clear, int n_workers) {
    WriteLines(file_name, compressor, Dump(canonical, clear, n_workers));
  }

  // Loads data from a vector of strings.
  static KmerSetSet Load(std::vector<std::string> v, bool canonical,
                         int n_workers) {
    absl::flat_hash_map<int, std::vector<int>> children;
    std::vector<KmerSet<K, KeyType>> kmer_sets;

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

    int kmer_sets_size;
    {
      std::stringstream ss(v[offset]);
      ss >> kmer_sets_size;
    }

    kmer_sets.resize(kmer_sets_size);

    {
      boost::asio::thread_pool pool(n_workers);

      for (int i = 0; i < kmer_sets_size; i++) {
        boost::asio::post(pool, [&, i] {
          std::string& s = v[offset + 1 + i];
          KmerSetCompressed<K, KeyType> kmer_set_compressed =
              KmerSetCompressed<K, KeyType>::Load(absl::StrSplit(s, ' '));
          kmer_sets[i] = kmer_set_compressed.ToKmerSet(canonical, n_workers);
          std::string().swap(s);
        });
      }

      pool.join();
    }

    return KmerSetSet(children, kmer_sets);
  }

 private:
  KmerSetSet(absl::flat_hash_map<int, std::vector<int>> children,
             std::vector<KmerSet<K, KeyType>> kmer_sets)
      : children_(children), kmer_sets_(kmer_sets) {}

  absl::flat_hash_map<int, std::vector<int>> children_;
  std::vector<KmerSet<K, KeyType>> kmer_sets_;
};

// KmerSetMM is similar to KmerSetSet, but it uses maximum weighted matching.
// If "approximate_matching" is true, it uses an approximate matching algorithm.
// If "approximate_weights" is true, it uses approximate weights for the
// matching algorithm.
// If "approximate_graph" is true, it uses a aparse graph for the matching
// algorithm, cutting half of edges.
// If "partial_matching" is true, only half of the nodes are matched. The larger
// ones are selected.
// The same operation can be repeated multiple times,
// resulting in a recursive structure. "recursion_limit" can be used to
// configure how deep it goes.
template <int K, typename KeyType>
class KmerSetSetMM {
 public:
  KmerSetSetMM(const KmerSetSetMM&) = delete;
  KmerSetSetMM& operator=(const KmerSetSetMM&) = delete;
  KmerSetSetMM(KmerSetSetMM&&) = default;
  KmerSetSetMM& operator=(KmerSetSetMM&&) = default;

  KmerSetSetMM(std::vector<KmerSet<K, KeyType>> kmer_sets, int recursion_limit,
               bool approximate_matching, bool approximate_weights,
               bool approximate_graph, bool partial_matching, int n_workers) {
    const int n = kmer_sets.size();

    // If mates[i] is j (and mates[j] is i), we have a match between
    // kmer_sets[i] and kmer_sets[j].
    absl::flat_hash_map<int, int> mates;

    spdlog::debug("calculating mates");

    {
      // "partial_matching" affects what kmer sets (nodes) to consider.
      absl::flat_hash_set<int> nodes;

      if (partial_matching) {
        // Considers the larger kmer sets only.

        std::vector<int> v;
        for (int i = 0; i < n; i++) v.push_back(i);

        std::sort(v.begin(), v.end(), [&](int lhs, int rhs) {
          return kmer_sets[lhs].Size() > kmer_sets[rhs].Size();
        });

        for (int i = 0; i < n / 2; i++) nodes.insert(v[i]);
      } else {
        // Considers all the kmer sets.

        for (int i = 0; i < n; i++) nodes.insert(i);
      }

      // The first two are for nodes, and the last one is for weights.
      std::vector<std::tuple<int, int, int64_t>> edges;

      // Calculates edges. The topology is affected by "approximate_graph" and
      // the weights are affected by "approximate_weights".
      {
        absl::InsecureBitGen bitgen;
        boost::asio::thread_pool pool(n_workers);
        std::mutex mu;

        for (int i : nodes) {
          for (int j : nodes) {
            // Mekes sure i < j.
            if (i >= j) continue;

            if (approximate_graph && absl::Bernoulli(bitgen, 0.5)) continue;

            boost::asio::post(pool, [&, i, j] {
              int64_t weight =
                  approximate_weights
                      ? kmer_sets[i].CommonEstimate(kmer_sets[j], 0.1)
                      : kmer_sets[i].Common(kmer_sets[j], 1);

              std::lock_guard lck(mu);
              edges.emplace_back(i, j, weight);
            });
          }
        }

        pool.join();
      }

      if (approximate_matching) {
        // Executes approximate max-weighted matching.

        // Sorts "edges" so that the heaviest edge comes first.
        std::sort(edges.begin(), edges.end(),
                  [](const std::tuple<int, int, int64_t>& lhs,
                     const std::tuple<int, int, int64_t>& rhs) {
                    return std::get<2>(lhs) > std::get<2>(rhs);
                  });

        for (const std::tuple<int, int, int64_t>& edge : edges) {
          int i, j;
          std::tie(i, j, std::ignore) = edge;

          // If i and j do not have their mates, make a match.
          if (mates.find(i) == mates.end() && mates.find(j) == mates.end()) {
            mates[i] = j;
            mates[j] = i;
          }
        }
      } else {
        // Executes exact max-weighted matching.

        // Boost Graph Library was used at first here, but there seemed to be a
        // bug. Therefore, lemon library is used.

        lemon::ListGraph g;

        std::vector<lemon::ListGraph::Node> nodes(n);  // id to node.
        lemon::ListGraph::NodeMap<int> ids(g);         // node to id.

        // Adds nodes.
        for (int i = 0; i < n; i++) {
          lemon::ListGraph::Node node = g.addNode();
          ids[node] = i;
          nodes[i] = node;
        }

        lemon::ListGraph::EdgeMap<int64_t> weights(g);

        // Adds edges and sets their weights.
        for (const std::tuple<int, int, int64_t>& edge : edges) {
          int i, j;
          int64_t weight;
          std::tie(i, j, weight) = edge;

          lemon::ListGraph::Edge e = g.addEdge(nodes[i], nodes[j]);
          weights[e] = weight;
        }

        lemon::MaxWeightedMatching<lemon::ListGraph,
                                   lemon::ListGraph::EdgeMap<int64_t>>
            matching(g, weights);

        matching.run();

        for (int i = 0; i < n; i++) {
          lemon::ListGraph::Node mate = matching.mate(nodes[i]);

          if (mate != lemon::INVALID) {
            mates[i] = ids[mate];
          }
        }
      }
    }

    spdlog::debug("calculated mates");

    std::vector<KmerSet<K, KeyType>> diffs;

    // Calculates parents, diffs, and diff_table_.
    // parents are calculated in the main thread, and the others are
    // calculated in sub threads.
    {
      int next_parent_id = n;

      std::mutex mu;
      boost::asio::thread_pool pool(n_workers);

      for (int i = 0; i < n; i++) {
        if (mates.find(i) == mates.end()) {
          // If kmer_sets[i] has no mates, sets its parent to an empty set.
          parents_[i] = -1;

          boost::asio::post(pool, [&, i] {
            KmerSet<K, KeyType> kmer_set = kmer_sets[i];

            std::lock_guard lck(mu);
            diff_table_[-1][i] = diffs.size();
            diffs.push_back(std::move(kmer_set));
          });
        } else {
          // If kmer_sets[i] has a mate, creates a new node by intersecting
          // them and sets their parents to that node.
          int mate = mates[i];

          if (i < mate) {
            int parent = next_parent_id++;

            parents_[i] = parent;
            parents_[mate] = parent;

            // -1 is an empty set.
            parents_[parent] = -1;

            boost::asio::post(pool, [&, i, mate, parent] {
              KmerSet<K, KeyType> sub = Sub(kmer_sets[i], kmer_sets[mate], 1);

              std::lock_guard lck(mu);
              diff_table_[parent][i] = diffs.size();
              diffs.push_back(std::move(sub));
            });

            boost::asio::post(pool, [&, i, mate, parent] {
              KmerSet<K, KeyType> sub = Sub(kmer_sets[mate], kmer_sets[i], 1);

              std::lock_guard lck(mu);
              diff_table_[parent][mate] = diffs.size();
              diffs.push_back(std::move(sub));
            });

            boost::asio::post(pool, [&, i, mate, parent] {
              KmerSet<K, KeyType> intersection =
                  Intersection(kmer_sets[i], kmer_sets[mate], 1);

              std::lock_guard lck(mu);
              diff_table_[-1][parent] = diffs.size();
              diffs.push_back(std::move(intersection));
            });
          }
        }
      }

      pool.join();
    }

    // kmer_sets is no longer needed.
    std::vector<KmerSet<K, KeyType>>().swap(kmer_sets);

    // Calculates cost.
    {
      std::atomic_int64_t cost = 0;

      boost::asio::thread_pool pool(n_workers);

      for (size_t i = 0; i < diffs.size(); i++) {
        boost::asio::post(pool, [&, i] { cost += diffs[i].Size(); });
      }

      pool.join();

      cost_ = cost;
    }

    spdlog::debug("cost_ = {}", cost_);

    if (recursion_limit == 0) {
      diffs_ = std::move(diffs);
    } else {
      diffs_ = new KmerSetSetMM(std::move(diffs), recursion_limit - 1,
                                approximate_matching, approximate_weights,
                                approximate_graph, partial_matching, n_workers);
    }
  }

  ~KmerSetSetMM() {
    if (!IsTerminal()) delete std::get<KmerSetSetMM*>(diffs_);
  }

  int64_t Cost() const {
    if (IsTerminal()) return cost_;
    return std::get<KmerSetSetMM*>(diffs_)->Cost();
  }

  // Reconstructs the ith kmer set.
  KmerSet<K, KeyType> Get(int i, int n_workers) const {
    // Returns the diff from "from" to "to".
    const auto GetDiff = [&](int from, int to) {
      int pos = diff_table_.find(from)->second.find(to)->second;

      if (IsTerminal()) {
        return std::get<std::vector<KmerSet<K, KeyType>>>(diffs_)[pos];
      } else {
        return std::get<KmerSetSetMM*>(diffs_)->Get(pos, n_workers);
      }
    };

    int parent = parents_.find(i)->second;

    if (parent == -1) {
      return GetDiff(-1, i);
    }

    return Add(GetDiff(-1, parent), GetDiff(parent, i), n_workers);
  }

  // Dumps data to a vector of strings.
  // "canonical" will be used to compress kmer sets with KmerSetCompressed.
  // "v" is used internally to process recursions.
  std::vector<std::string> Dump(
      bool canonical, int n_workers,
      std::vector<std::string> v = std::vector<std::string>()) const {
    // Dumps cost_ to a string.
    {
      std::stringstream ss;
      ss << cost_;
      v.push_back(ss.str());
    }

    // Dumps parents_ to a string.
    // Format: "size key value key value ..."
    {
      std::stringstream ss;
      ss << parents_.size();

      for (const auto& p : parents_) {
        ss << ' ' << p.first << ' ' << p.second;
      }

      v.push_back(ss.str());
    }

    // Dumps diff_table_ to a string.
    // Format: "size key1 key2 value key1 key2 value ..."
    {
      std::vector<std::tuple<int, int, int>> flatten;
      for (const auto& p1 : diff_table_) {
        for (const auto& p2 : p1.second) {
          flatten.emplace_back(p1.first, p2.first, p2.second);
        }
      }

      std::stringstream ss;
      ss << flatten.size();

      for (const auto& t : flatten) {
        ss << ' ' << std::get<0>(t) << ' ' << std::get<1>(t) << ' '
           << std::get<2>(t);
      }

      v.push_back(ss.str());
    }

    if (IsTerminal()) {
      // Section delimiter.
      v.push_back("");

      const std::vector<KmerSet<K, KeyType>>& diffs =
          std::get<std::vector<KmerSet<K, KeyType>>>(diffs_);

      // Dumps diffs to 1 + diffs_.size() strings.
      // Each KmerSet is converted to KmerSetCompressed and dumped to a
      // string.

      {
        std::stringstream ss;
        ss << diffs.size();
        v.push_back(ss.str());
      }

      {
        std::vector<std::string> buf(diffs.size());

        boost::asio::thread_pool pool(n_workers);

        for (size_t i = 0; i < diffs.size(); i++) {
          boost::asio::post(pool, [&, i] {
            const KmerSetCompressed<K, KeyType> kmer_set_compressed =
                KmerSetCompressed<K, KeyType>::FromKmerSet(diffs[i], canonical,
                                                           1);

            buf[i] = absl::StrJoin(kmer_set_compressed.Dump(), " ");
          });
        }

        pool.join();

        v.reserve(v.size() + buf.size());
        for (std::string& s : buf) v.push_back(std::move(s));
      }

      return v;
    } else {
      // Recursive
      return std::get<KmerSetSetMM*>(diffs_)->Dump(canonical, n_workers,
                                                   std::move(v));
    }
  }

  // Dumps data to a file.
  void Dump(const std::string& file_name, const std::string& compressor,
            bool canonical, int n_workers) const {
    WriteLines(file_name, compressor, Dump(canonical, n_workers));
  }

  // Loads data from a vector of strings.
  // "i" represents the start index and is used internally to process
  // recursions.
  static KmerSetSetMM* Load(const std::vector<std::string>& v, bool canonical,
                            int n_workers, int64_t i = 0) {
    int64_t cost;
    absl::flat_hash_map<int, int> parents;
    absl::flat_hash_map<int, absl::flat_hash_map<int, int>> diff_table;

    // Loads cost.
    {
      std::stringstream ss(v[i]);
      ss >> cost;
    }

    // Loads parents.
    {
      std::stringstream ss(v[i + 1]);
      int64_t size;
      ss >> size;

      for (int64_t i = 0; i < size; i++) {
        int key, value;
        ss >> key >> value;
        parents[key] = value;
      }
    }

    // Loads diff_table.
    {
      std::stringstream ss(v[i + 2]);
      int64_t size;
      ss >> size;

      for (int64_t i = 0; i < size; i++) {
        int key1, key2, value;
        ss >> key1 >> key2 >> value;
        diff_table[key1][key2] = value;
      }
    }

    if (v[i + 3] == "") {
      // The end of recursion.

      // Loads diffs.

      int64_t size;
      {
        std::stringstream ss(v[i + 4]);
        ss >> size;
      }

      std::vector<KmerSet<K, KeyType>> diffs(size);

      boost::asio::thread_pool pool(n_workers);

      for (int64_t j = 0; j < size; j++) {
        boost::asio::post(pool, [&, j] {
          const std::string& line = v[i + 5 + j];

          if (line.empty()) {
            // diffs[j] should be an empty set.
            return;
          }

          const KmerSetCompressed<K, KeyType> kmer_set_compressed =
              KmerSetCompressed<K, KeyType>::Load(absl::StrSplit(line, ' '));

          diffs[j] = kmer_set_compressed.ToKmerSet(canonical, 1);
        });
      }

      pool.join();

      return new KmerSetSetMM(cost, parents, diff_table, std::move(diffs));
    } else {
      return new KmerSetSetMM(cost, parents, diff_table,
                              Load(v, canonical, n_workers, i + 3));
    }
  }

  static KmerSetSetMM* Load(const std::string& file_name,
                            const std::string& decompressor, bool canonical,
                            int n_workers) {
    return Load(ReadLines(file_name, decompressor), canonical, n_workers);
  }

 private:
  KmerSetSetMM(
      int64_t cost, absl::flat_hash_map<int, int> parents,
      absl::flat_hash_map<int, absl::flat_hash_map<int, int>> diff_table,
      std::vector<KmerSet<K, KeyType>> diffs)
      : cost_(cost),
        parents_(std::move(parents)),
        diff_table_(std::move(diff_table)),
        diffs_(std::move(diffs)) {}

  KmerSetSetMM(
      int64_t cost, absl::flat_hash_map<int, int> parents,
      absl::flat_hash_map<int, absl::flat_hash_map<int, int>> diff_table,
      KmerSetSetMM* diffs)
      : cost_(cost),
        parents_(std::move(parents)),
        diff_table_(std::move(diff_table)),
        diffs_(diffs) {}

  // Returns true if it is the end of recursion.
  bool IsTerminal() const {
    return std::holds_alternative<std::vector<KmerSet<K, KeyType>>>(diffs_);
  }

  int64_t cost_ = 0;

  // If parents_[i] == -1, the parent of i is an empty set.
  absl::flat_hash_map<int, int> parents_;

  // If diffs_ is a vector, diffs_[diff_table[i][j]] is a diff from i to j.
  // If diffs_ is a KmerSetSetMM*, diffs_.Get(diff_table[i][j]) is a diff from
  // i to j.
  absl::flat_hash_map<int, absl::flat_hash_map<int, int>> diff_table_;
  std::variant<std::vector<KmerSet<K, KeyType>>, KmerSetSetMM*> diffs_;
};

#endif

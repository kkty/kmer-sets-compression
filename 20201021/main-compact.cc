#include <glog/logging.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <bitset>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>
#include <mutex>
#include <queue>
#include <sparsehash/dense_hash_map>
#include <sparsehash/sparse_hash_map>
#include <sparsehash/sparse_hash_set>
#include <sstream>
#include <stack>
#include <string>
#include <vector>

std::vector<std::vector<bool>> ParseFASTQ(std::istream& is) {
  std::array<std::string, 4> buf;

  // Reads 4 lines from "is" and pushes to "buf".
  auto consume = [&is, &buf]() {
    for (int i = 0; i < 4; i++) {
      if (!std::getline(is, buf[i])) return false;
    }

    return true;
  };

  std::vector<std::string> reads_s;
  while (consume() && buf[0].length() && buf[0][0] == '@') {
    reads_s.push_back(buf[1]);
  }

  const int size = reads_s.size();

  LOG(INFO) << "size = " << size;

  std::vector<std::vector<bool>> reads_b(size);

  std::atomic_int count_a = 0;
  std::atomic_int count_c = 0;
  std::atomic_int count_g = 0;
  std::atomic_int count_t = 0;
  std::atomic_int count_unknown = 0;

#pragma omp parallel for
  for (int i = 0; i < size; i++) {
    const std::string& read_s = reads_s[i];
    std::vector<bool>& read_b = reads_b[i];

    read_b.reserve(read_s.length() * 3);

    bool skip = false;

    for (const char c : read_s) {
      switch (c) {
        case 'A':
          read_b.push_back(false);
          read_b.push_back(false);
          read_b.push_back(false);
          count_a++;
          break;
        case 'C':
          read_b.push_back(false);
          read_b.push_back(false);
          read_b.push_back(true);
          count_c++;
          break;
        case 'G':
          read_b.push_back(false);
          read_b.push_back(true);
          read_b.push_back(false);
          count_g++;
          break;
        case 'T':
          read_b.push_back(false);
          read_b.push_back(true);
          read_b.push_back(true);
          count_t++;
          break;
        default:
          // Unknown character.
          read_b.push_back(true);
          read_b.push_back(false);
          read_b.push_back(false);
          count_unknown++;
      }
    }
  }

  LOG(INFO) << "count_a = " << count_a;
  LOG(INFO) << "count_c = " << count_c;
  LOG(INFO) << "count_g = " << count_g;
  LOG(INFO) << "count_t = " << count_t;
  LOG(INFO) << "count_unknown = " << count_unknown;

  return reads_b;
}

std::vector<std::vector<bool>> ParseFASTQ(const std::string& path) {
  std::ifstream is{path};
  return ParseFASTQ(is);
}

template <int K, int8_t CUTOFF = 2>
std::vector<std::bitset<K * 2>> FindKmers(
    const std::vector<std::vector<bool>>& reads) {
  google::sparse_hash_map<std::bitset<K * 2>, int> kmers;

  {
    // Number of kmers, including duplicates;
    std::atomic_int64_t count = 0;

#pragma omp parallel for
    for (const auto& read : reads) {
      count += read.size() / 3 - K + 1;
    }

    // Rough estimation.
    kmers.resize(count / 10);
  }

  std::mutex mu;

  std::atomic_uint skipped = 0;

#pragma omp parallel for
  for (const auto& read : reads) {
    std::vector<std::bitset<K * 2>> kmers_buf;

    for (int i = 0; i + K <= read.size() / 3; i++) {
      std::bitset<K * 2> kmer;

      bool has_unkown = false;
      for (int j = 0; j < K; j++) {
        // Looking at the following.
        // - kmer[j * 2], kmer[j * 2 + 1]
        // - read[(i + j) * 3], read[(i + j) * 3 + 1], read[(i + j) * 3 + 2]

        if (read[(i + j) * 3]) {
          has_unkown = true;
          break;
        }

        kmer.set(j * 2, read[(i + j) * 3 + 1]);
        kmer.set(j * 2 + 1, read[(i + j) * 3 + 2]);
      }

      if (has_unkown) {
        skipped += 1;
        continue;
      }

      kmers_buf.push_back(kmer);
    }

    {
      std::unique_lock lck{mu};

      for (const auto& kmer : kmers_buf) kmers[kmer] += 1;
    }
  }

  LOG(INFO) << "skipped = " << skipped;

  std::vector<std::bitset<K * 2>> v;
  v.reserve(kmers.size());
  int cut = 0;

  for (const auto& p : kmers) {
    const auto& [kmer, count] = p;
    if (count < CUTOFF) {
      cut++;
      continue;
    }
    v.push_back(kmer);
  }

  LOG(INFO) << "cut = " << cut;

  v.shrink_to_fit();

  return v;
}

struct Node {
  Node* next_a = nullptr;
  Node* next_g = nullptr;
  Node* next_c = nullptr;
  Node* next_t = nullptr;
  int64_t index = -1;
};

template <int K>
struct Graph {
  google::sparse_hash_map<std::bitset<K * 2>, Node*> nodes;
  google::dense_hash_map<int64_t, int64_t> distances;
  int64_t edges;
};

template <int K>
Graph<K> ConstructGraph(const std::vector<std::bitset<K * 2>>& kmers) {
  google::sparse_hash_map<std::bitset<K * 2>, Node*> nodes;
  nodes.resize(kmers.size());

  for (const auto& kmer : kmers) {
    nodes[kmer] = nullptr;
  }

#pragma omp parallel for
  for (const auto& kmer : kmers) {
    nodes[kmer] = new Node;
  }

  LOG(INFO) << "adding edges";

  std::atomic_int64_t edges = 0;

#pragma omp parallel for
  for (const auto& from_s : kmers) {
    Node* const from_node = nodes[from_s];

    for (const char c : std::vector<char>{'A', 'C', 'G', 'T'}) {
      // Careful!
      std::bitset<K* 2> to_s = from_s >> 2;

      switch (c) {
        case 'A':
          break;
        case 'C':
          to_s.set(K * 2 - 1);
          break;
        case 'G':
          to_s.set(K * 2 - 2);
          break;
        case 'T':
          to_s.set(K * 2 - 1);
          to_s.set(K * 2 - 2);
          break;
      }

      auto it = nodes.find(to_s);
      if (it == nodes.end()) continue;

      edges += 1;

      Node* const to_node = it->second;

      switch (c) {
        case 'A':
          from_node->next_a = to_node;
          break;
        case 'C':
          from_node->next_c = to_node;
          break;
        case 'G':
          from_node->next_g = to_node;
          break;
        case 'T':
          from_node->next_t = to_node;
          break;
      }
    }
  }

  LOG(INFO) << "added edges";

  LOG(INFO) << "edges = " << edges;

  LOG(INFO) << "calculating index";

  int64_t next_index = 0;

  for (const auto& kmer : kmers) {
    if (nodes[kmer]->index != -1) continue;

    std::stack<Node*> q;
    q.push(nodes[kmer]);

    while (q.size()) {
      Node* const from = q.top();
      q.pop();

      from->index = next_index++;

      for (Node* to : std::array<Node*, 4>{from->next_a, from->next_c,
                                           from->next_g, from->next_t}) {
        if (!to || to->index != -1) continue;
        q.push(to);
      }
    }
  }

  LOG(INFO) << "calculated index";

  LOG(INFO) << "calculating distances";

  google::dense_hash_map<int64_t, int64_t> distances;
  distances.set_empty_key(-1);

  for (auto it = kmers.begin(); it != kmers.end(); ++it) {
    const auto& kmer = *it;

    Node* const from = nodes[kmer];

    for (Node* const to : std::vector<Node*>{from->next_a, from->next_c,
                                             from->next_g, from->next_t}) {
      if (!to) continue;
      distances[std::abs(to->index - from->index)] += 1;
    }
  }

  LOG(INFO) << "calculated distances";

  return Graph<K>{std::move(nodes), std::move(distances), edges};
}

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = 1;

  const int K = 31;

  std::vector<std::vector<bool>> reads;

  if (argc == 1) {
    reads = ParseFASTQ(std::cin);
  } else if (argc == 2) {
    std::string path{argv[1]};
    reads = ParseFASTQ(path);
  } else {
    LOG(INFO) << "invalid arguments";
    exit(1);
  }

  LOG(INFO) << "reads.size() = " << reads.size();

  auto kmers = FindKmers<K>(reads);
  LOG(INFO) << "kmers.size() = " << kmers.size();

  LOG(INFO) << "constructing the graph";
  auto graph = ConstructGraph<K>(kmers);
  LOG(INFO) << "constructed the graph";

  {
    auto& distances = graph.distances;
    int64_t edges = graph.edges;

    std::vector<int64_t> keys;
    for (auto const& p : distances) keys.push_back(p.first);
    std::sort(keys.begin(), keys.end());

    double sum = 0;
    std::cout << "i" << '\t' << "distances[i]" << '\t' << "distances[i] / edges"
              << '\t' << "sum" << std::endl;
    for (int64_t i : keys) {
      std::cout << i << '\t' << distances[i] << '\t'
                << (double)distances[i] / edges << '\t';
      sum += (double)distances[i] / edges;
      std::cout << sum << std::endl;
    }
  }
}
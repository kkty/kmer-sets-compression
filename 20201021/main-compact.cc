#include <omp.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <bitset>
#include <filesystem>
#include <fstream>
#include <google/dense_hash_map>
#include <google/sparse_hash_map>
#include <google/sparse_hash_set>
#include <iostream>
#include <mutex>
#include <queue>
#include <stack>
#include <string>
#include <vector>

std::vector<std::vector<bool>> ParseFASTQ(std::istream& is) {
  std::array<std::string, 4> lines;

  auto read_lines = [&]() {
    for (int i = 0; i < 4; i++) {
      if (!std::getline(is, lines[i])) return false;
    }
    return true;
  };

  std::vector<std::vector<bool>> reads;

  int skipped = 0;

  while (read_lines() && lines[0].length() && lines[0][0] == '@') {
    std::string read_s = lines[1];
    std::vector<bool> read_b;
    read_b.reserve(read_s.length() * 2);

    bool skip = false;

    for (const char c : read_s) {
      switch (c) {
        case 'A':
          read_b.push_back(false);
          read_b.push_back(false);
          break;
        case 'C':
          read_b.push_back(false);
          read_b.push_back(true);
          break;
        case 'G':
          read_b.push_back(true);
          read_b.push_back(false);
          break;
        case 'T':
          read_b.push_back(true);
          read_b.push_back(true);
          break;
        default:
          skip = true;
      }
    }

    if (skip) {
      skipped += 1;
      continue;
    }

    reads.push_back(read_b);
  }

  std::cerr << "skipped = " << skipped << std::endl;

  return reads;
}

std::vector<std::vector<bool>> ParseFASTQ(const std::filesystem::path& path) {
  std::ifstream is{path};
  return ParseFASTQ(is);
}

template <int K, int8_t CUTOFF = 2>
std::vector<std::bitset<K * 2>> FindKmers(
    const std::vector<std::vector<bool>>& reads) {
  google::sparse_hash_map<std::bitset<K * 2>, int> kmers;

  {
    std::atomic_int64_t count = 0;
#pragma omp parallel for
    for (const auto& read : reads) {
      count += read.size() / 2 - K + 1;
    }
    kmers.resize(count / 10);
  }

  std::mutex mu;

#pragma omp parallel for
  for (const auto& read : reads) {
    std::vector<std::bitset<K * 2>> kmers_buf;

    for (int i = 0; i + K <= read.size() / 2; i++) {
      std::bitset<K * 2> kmer;

      for (int j = 0; j < K * 2; j++) {
        kmer.set(j, read[i * 2 + j]);
      }

      kmers_buf.push_back(kmer);
    }

    {
      std::unique_lock lck{mu};

      for (const auto& kmer : kmers_buf) kmers[kmer] += 1;
    }
  }

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

  std::cerr << "cut = " << cut << std::endl;

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

  std::cerr << "adding edges" << std::endl;

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

  std::cerr << "added edges" << std::endl;

  std::cerr << "edges = " << edges << std::endl;

  std::cerr << "calculating index" << std::endl;

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

  std::cerr << "calculated index" << std::endl;

  std::cerr << "calculating distances" << std::endl;

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

  std::cerr << "calculated distances" << std::endl;

  return Graph<K>{std::move(nodes), std::move(distances), edges};
}

int main(int argc, char* argv[]) {
  const int K = 31;

  std::vector<std::vector<bool>> reads;

  if (argc == 1) {
    reads = ParseFASTQ(std::cin);
  } else if (argc == 2) {
    std::filesystem::path path{argv[1]};
    reads = ParseFASTQ(path);
  } else {
    std::cerr << "invalid arguments" << std::endl;
    exit(1);
  }

  std::cerr << "reads.size() = " << reads.size() << std::endl;

  auto kmers = FindKmers<K>(reads);
  std::cerr << "kmers.size() = " << kmers.size() << std::endl;

  std::cerr << "constructing the graph" << std ::endl;
  auto graph = ConstructGraph<K>(kmers);
  std::cerr << "constructed the graph" << std ::endl;

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
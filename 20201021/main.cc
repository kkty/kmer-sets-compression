#include <algorithm>
#include <array>
#include <atomic>
#include <execution>
#include <filesystem>
#include <fstream>
#include <google/dense_hash_map>
#include <google/dense_hash_set>
#include <iostream>
#include <mutex>
#include <optional>
#include <queue>
#include <stack>
#include <string>
#include <vector>

std::vector<std::string> parse_fastq(std::istream& is) {
  std::array<std::string, 4> lines;

  auto read_lines = [&]() {
    for (int i = 0; i < 4; i++) {
      if (!std::getline(is, lines[i])) return false;
    }
    return true;
  };

  std::vector<std::string> reads;

  while (read_lines() && lines[0].length() && lines[0][0] == '@') {
    reads.push_back(lines[1]);
  }

  return reads;
}

std::vector<std::string> parse_fastq(std::filesystem::path path) {
  std::ifstream is{path};
  return parse_fastq(is);
}

template <int K>
std::vector<std::string> find_kmers(std::vector<std::string> ss) {
  google::dense_hash_set<std::string> kmers;
  kmers.set_empty_key("X");

  std::mutex mu;

  std::for_each(std::execution::par_unseq, ss.begin(), ss.end(),
                [&](const std::string& s) {
                  std::vector<std::string> kmers_buf;
                  kmers_buf.reserve(s.length() - K + 1);

                  for (int i = 0; i + K <= s.length(); i++) {
                    const std::string kmer = s.substr(i, K);

                    for (const char c : kmer) {
                      if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
                        if (c == 'N') {
                          return;
                        } else {
                          std::cerr << "unknown character: " << c << std::endl;
                          exit(1);
                        }
                      }
                    }

                    kmers_buf.push_back(kmer);
                  }

                  {
                    std::unique_lock lck{mu};
                    for (const std::string& kmer : kmers_buf)
                      kmers.insert(kmer);
                  }
                });

  std::vector<std::string> v;
  v.reserve(kmers.size());
  for (const std::string& s : kmers) {
    v.push_back(s);
  }

  return v;
}

struct Node {
  Node* next_a = nullptr;
  Node* next_g = nullptr;
  Node* next_c = nullptr;
  Node* next_t = nullptr;
  int index = -1;
};

struct Graph {
  Node* root;
  int size;
};

template <int K>
Graph construct_graph(const std::vector<std::string>& kmers) {
  google::dense_hash_map<std::string, Node*> nodes;
  nodes.set_empty_key("X");

  std::for_each(std::execution::seq, kmers.begin(), kmers.end(),
                [&](const std::string& s) { nodes[s] = nullptr; });

  std::for_each(std::execution::par_unseq, kmers.begin(), kmers.end(),
                [&](const std::string& s) { nodes[s] = new Node; });

  const int size = nodes.size();

  std::atomic_int edges = 0;

  std::for_each(std::execution::par_unseq, kmers.begin(), kmers.end(),
                [&](const std::string& from_s) {
                  Node* const from_node = nodes[from_s];
                  const std::string from_s_suffix = from_s.substr(1);

                  for (const char c : std::vector<char>{'A', 'C', 'G', 'T'}) {
                    const std::string to_s = from_s_suffix + c;
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
                });

  std::cerr << "edges = " << edges << std::endl;

  int next_index = 0;

  for (const std::string& kmer : kmers) {
    if (nodes[kmer]->index != -1) continue;

    std::stack<Node*> q;
    q.push(nodes[kmer]);

    while (q.size()) {
      Node* const from = q.top();
      q.pop();

      from->index = next_index++;

      for (Node* to : std::vector<Node*>{from->next_a, from->next_c,
                                         from->next_g, from->next_t}) {
        if (!to || to->index != -1) continue;
        q.push(to);
      }
    }
  }

  google::dense_hash_map<int, int> distances;
  distances.set_empty_key(-1);

  std::for_each(
      std::execution::seq, kmers.begin(), kmers.end(),
      [&](const std::string& kmer) {
        Node* const from = nodes[kmer];

        for (Node* const to : std::vector<Node*>{from->next_a, from->next_c,
                                                 from->next_g, from->next_t}) {
          if (!to) continue;
          distances[std::abs(to->index - from->index)] += 1;
        }
      });

  double sum = 0;
  for (int i = 0; i < 100; i++) {
    std::cerr << "distances[" << i << "] = " << distances[i] << '\t';
    std::cerr << "distances[" << i
              << "] / edges = " << (double)distances[i] / edges << '\t';
    sum += (double)distances[i] / edges;
    std::cerr << "sum = " << sum << std::endl;
  }

  return Graph{};
}

int main(int argc, char* argv[]) {
  const int K = 21;

  auto reads = parse_fastq(std::cin);
  std::cout << "reads.size() = " << reads.size() << std::endl;

  auto kmers = find_kmers<K>(reads);
  std::cout << "kmers.size() = " << kmers.size() << std::endl;

  std::cout << "constructing the graph" << std ::endl;
  auto graph = construct_graph<K>(kmers);
  std::cout << "constructed the graph" << std ::endl;
  std::cout << "graph.size = " << graph.size << std::endl;
}
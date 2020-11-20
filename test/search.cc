#include "search.h"

#include <algorithm>
#include <vector>

#include "core/kmer.h"
#include "core/kmer_set.h"
#include "gtest/gtest.h"

TEST(search, ConstructGraphCanonical) {
  const int K = 7;
  using KeyType = uint8_t;
  KmerSet<K, KeyType> kmer_set;

  kmer_set.Add(Kmer<K>("AAAAAAA"));
  kmer_set.Add(Kmer<K>("AAAAAAT"));  // branching

  kmer_set.Add(Kmer<K>("AAAAATC"));
  kmer_set.Add(Kmer<K>("AAAATCC"));
  kmer_set.Add(Kmer<K>("AAATCCC"));  // branching
  kmer_set.Add(Kmer<K>("AATCCCC"));
  kmer_set.Add(Kmer<K>("AATCCCG"));

  kmer_set.Add(Kmer<K>("AAAAATG"));
  kmer_set.Add(Kmer<K>("AAAATGG"));
  kmer_set.Add(Kmer<K>("AAATGGG"));
  kmer_set.Add(Kmer<K>("AATGGGG"));  // branching
  kmer_set.Add(Kmer<K>("ATGGGGG"));
  kmer_set.Add(Kmer<K>("ATGGGGC"));

  KmerGraph<K> g = ConstructKmerGraph(kmer_set, 1);

  {
    std::vector<std::pair<int64_t, int64_t>> edges =
        g.edges[g.ids[Kmer<K>("AAAAAAT")]];

    ASSERT_EQ(std::find_if(edges.begin(), edges.end(),
                           [&](const std::pair<int64_t, int64_t>& p) {
                             return p.first == g.ids[Kmer<K>("AAATCCC")];
                           })
                  ->second,
              3

    );
  }

  {
    std::vector<std::pair<int64_t, int64_t>> edges =
        g.edges[g.ids[Kmer<K>("AAAAAAT")]];

    ASSERT_EQ(std::find_if(edges.begin(), edges.end(),
                           [&](const std::pair<int64_t, int64_t>& p) {
                             return p.first == g.ids[Kmer<K>("AATGGGG")];
                           })
                  ->second,
              4

    );
  }

  {
    SearchResult result =
        DijkstraSearch(g, Kmer<K>("AAAAAAT"), Kmer<K>("AAATCCC"));

    ASSERT_TRUE(result.found);
    ASSERT_EQ(result.distance, 3);
  }

  {
    SearchResult result =
        AStarSearch<K>(g, Kmer<K>("AAAAAAT"), Kmer<K>("AAATCCC"), 3);

    ASSERT_TRUE(result.found);
    ASSERT_EQ(result.distance, 3);
  }
}
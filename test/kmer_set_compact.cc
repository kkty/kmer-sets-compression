#include "kmer_set_compact.h"

#include <vector>

#include "gtest/gtest.h"
#include "kmer.h"

TEST(KmerSetCompact, Contains) {
  const int K = 5;
  const int B = 3;

  std::vector<Kmer<K>> kmers;
  kmers.emplace_back("AAAAA");
  kmers.emplace_back("AAAAC");
  kmers.emplace_back("AAACC");
  kmers.emplace_back("AAACA");

  KmerSet<K, B> kmer_set;
  for (const Kmer<K>& kmer : kmers) {
    kmer_set.Add(kmer);
  }

  KmerSetCompact<K> kmer_set_compact{kmer_set};

  for (const Kmer<K>& kmer : kmers) {
    ASSERT_TRUE(kmer_set_compact.Contains(kmer));
  }
}
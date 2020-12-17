#include "core/kmer_set_compact.h"

#include <string>

#include "core/kmer.h"
#include "core/kmer_set.h"
#include "gtest/gtest.h"

TEST(KmerSetCompact, FromAndToKmerSet) {
  const int K = 9;
  const int N = 10;
  using KeyType = uint8_t;

  KmerSet<K, N, KeyType> kmer_set;

  // Adding random kmers.
  for (int i = 0; i < 10000; i++) {
    kmer_set.Add(GetRandomKmer<K>().Canonical());
  }

  // Adding loops.
  for (int i = 0; i < 1000; i++) {
    const Kmer<K> kmer = GetRandomKmer<K>();
    const std::string s = kmer.String() + kmer.String();
    for (int j = 0; j < K; j++) {
      kmer_set.Add(Kmer<K>(s.substr(j, K)).Canonical());
    }
  }

  KmerSetCompact<K, N, KeyType> compressed =
      KmerSetCompact<K, N, KeyType>::FromKmerSet(kmer_set, true, false, 1);

  KmerSet<K, N, KeyType> decompressed = compressed.ToKmerSet(true, 1);

  ASSERT_TRUE(kmer_set.Equals(decompressed, 1));
}
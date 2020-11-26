#include "core/kmer_set_compressed.h"

#include <string>

#include "core/kmer.h"
#include "core/kmer_set.h"
#include "gtest/gtest.h"

TEST(KmerSetCompressed, FromAndToKmerSet) {
  const int K = 9;
  using KeyType = uint8_t;

  KmerSet<K, KeyType> kmer_set;

  // Adding random kmers.
  for (int i = 0; i < 10000; i++) {
    kmer_set.Add(GetRandomKmer<K>().Canonical());
  }

  // Adding loops.
  for (int i = 0; i < 1000; i++) {
    const Kmer<K> kmer = GetRandomKmer<K>();
    const std::string s = kmer.String() + kmer.String();
    for (int i = 0; i < K; i++) {
      kmer_set.Add(Kmer<K>(s.substr(i, K)).Canonical());
    }
  }

  KmerSetCompressed<K, KeyType> compressed =
      KmerSetCompressed<K, KeyType>::FromKmerSet(kmer_set, true, 1);

  KmerSet<K, KeyType> decompressed = compressed.ToKmerSet(true, 1);

  ASSERT_TRUE(kmer_set.Equals(decompressed, 1));
}
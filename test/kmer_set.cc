#include "kmer_set.h"

#include <bitset>
#include <tuple>

#include "gtest/gtest.h"
#include "kmer.h"
#include "kmer_counter.h"

TEST(kmer_set, BucketAndKey) {
  const int K = 5;
  const int N = 2;
  Kmer<K> kmer("AGCTG");
  const auto [bucket, key] = GetBucketAndKeyFromKmer<K, N>(kmer);
  ASSERT_EQ(kmer, (GetKmerFromBucketAndKey<K, N>(bucket, key)));
}

TEST(kmer_set, Find) {
  const int K = 5;
  const int B = 3;
  KmerSet<K, B> kmer_set;
  kmer_set.Add(Kmer<K>("AAAAA"));
  kmer_set.Add(Kmer<K>("CCCCC"));

  {
    const auto kmers = kmer_set.Find(
        [](const Kmer<K>& kmer) { return kmer.String()[0] == 'A'; });
    ASSERT_EQ(kmers.size(), 1);
    ASSERT_EQ(kmers[0].String(), "AAAAA");
  }

  {
    const auto kmers = kmer_set.Find(
        [](const Kmer<K>& kmer) { return kmer.String()[1] == 'C'; });
    ASSERT_EQ(kmers.size(), 1);
    ASSERT_EQ(kmers[0].String(), "CCCCC");
  }
}
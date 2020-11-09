#include "kmer_set.h"

#include <bitset>
#include <tuple>

#include "gtest/gtest.h"
#include "kmer.h"
#include "kmer_counter.h"

TEST(KmerSet, BucketAndKey) {
  const int K = 5;
  using KeyType = uint8_t;
  Kmer<K> kmer("AGCTG");
  const auto [bucket, key] = GetBucketAndKeyFromKmer<K, KeyType>(kmer);
  ASSERT_EQ(kmer.String(),
            (GetKmerFromBucketAndKey<K, KeyType>(bucket, key)).String());
}

TEST(KmerSet, AddRemove) {
  const int K = 5;
  using KeyType = uint8_t;

  Kmer<K> kmer("AAAAA");

  KmerSet<K, KeyType> kmer_set;

  ASSERT_EQ(kmer_set.Size(), 0);
  ASSERT_FALSE(kmer_set.Contains(kmer));

  kmer_set.Add(kmer);

  ASSERT_EQ(kmer_set.Size(), 1);
  ASSERT_TRUE(kmer_set.Contains(kmer));

  kmer_set.Remove(kmer);

  ASSERT_EQ(kmer_set.Size(), 0);
  ASSERT_FALSE(kmer_set.Contains(kmer));
}

TEST(KmerSet, Find) {
  const int K = 5;
  using KeyType = uint8_t;

  KmerSet<K, KeyType> kmer_set;
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

TEST(KmetSet, Operators) {
  const int K = 5;
  using KeyType = uint8_t;

  KmerSet<K, KeyType> kmer_set1;
  KmerSet<K, KeyType> kmer_set2;

  kmer_set1.Add(Kmer<K>("AAAAA"));
  kmer_set1.Add(Kmer<K>("TTTTT"));
  kmer_set1.Add(Kmer<K>("CCCCC"));

  kmer_set2.Add(Kmer<K>("AAAAA"));
  kmer_set2.Add(Kmer<K>("TTTTT"));
  kmer_set2.Add(Kmer<K>("GGGGG"));

  ASSERT_EQ((kmer_set1 + kmer_set2).Size(), 4);
  ASSERT_EQ((kmer_set1 - kmer_set2).Size(), 1);
  ASSERT_EQ((kmer_set2 - kmer_set1).Size(), 1);
  ASSERT_EQ((kmer_set2 - kmer_set1).Size(), 1);
  ASSERT_EQ((kmer_set2 & kmer_set1).Size(), 2);
}

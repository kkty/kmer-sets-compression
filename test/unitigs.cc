#include "core/unitigs.h"

#include <string>
#include <vector>

#include "core/kmer.h"
#include "core/kmer_set.h"
#include "gtest/gtest.h"

TEST(Unitigs, GetUnitigs) {
  const int K = 9;
  const int N = 10;
  using KeyType = uint8_t;

  KmerSet<K, N, KeyType> kmer_set;

  // Adding random kmers.
  for (int i = 0; i < 10000; i++) {
    kmer_set.Add(GetRandomKmer<K>());
  }

  // Adding loops.
  for (int i = 0; i < 1000; i++) {
    const Kmer<K> kmer = GetRandomKmer<K>();
    const std::string s = kmer.String() + kmer.String();
    for (int j = 0; j < K; j++) {
      kmer_set.Add(Kmer<K>(s.substr(j, K)));
    }
  }

  std::vector<std::string> unitigs = GetUnitigs(kmer_set, 1);

  KmerSet<K, N, KeyType> reconstructed;

  for (const std::string& s : unitigs) {
    const int n = s.length();
    ASSERT_TRUE(n >= K);

    for (int i = 0; i < n - K + 1; i++) {
      const Kmer<K> kmer = Kmer<K>(s.substr(i, K));
      ASSERT_FALSE(reconstructed.Contains(kmer));
      reconstructed.Add(kmer);
    }
  }

  ASSERT_TRUE(kmer_set.Equals(reconstructed, 1));
}

TEST(Unitigs, GetUnitigsCanonical) {
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

  std::vector<std::string> unitigs = GetUnitigsCanonical(kmer_set, 1);

  KmerSet<K, N, KeyType> reconstructed;

  for (const std::string& s : unitigs) {
    const int n = s.length();
    ASSERT_TRUE(n >= K);

    for (int i = 0; i < n - K + 1; i++) {
      const Kmer<K> kmer = Kmer<K>(s.substr(i, K)).Canonical();
      ASSERT_FALSE(reconstructed.Contains(kmer));
      reconstructed.Add(kmer);
    }
  }

  ASSERT_TRUE(kmer_set.Equals(reconstructed, 1));
}

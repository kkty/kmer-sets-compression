#include "core/unitigs.h"

#include <string>
#include <vector>

#include "absl/random/random.h"
#include "core/kmer.h"
#include "core/kmer_set.h"
#include "gtest/gtest.h"

template <int K, int N, typename KeyType>
KmerSet<K, N, KeyType> GetTestData() {
  absl::InsecureBitGen bitgen;

  KmerSet<K, N, KeyType> kmer_set;

  for (int i = 0; i < absl::Uniform(bitgen, 0, 10000); i++) {
    std::string s;

    for (int j = 0; j < absl::Uniform(bitgen, 0, 100); j++) {
      s += GetRandomKmer<K>().String();
    }

    // Creates a loop.
    if (absl::Uniform(bitgen, 0, 2) == 0) {
      s += s;
    }

    for (int j = 0; j < static_cast<int>(s.size()) - K + 1; j++) {
      kmer_set.Add(Kmer<K>(s.substr(j, K)).Canonical());
    }
  }

  return kmer_set;
}

TEST(Unitigs, Complement) { ASSERT_EQ(Complement("ACGTT"), "AACGT"); }

TEST(Unitigs, GetUnitigsRandom) {
  const int K = 9;
  const int N = 10;
  using KeyType = uint8_t;

  KmerSet<K, N, KeyType> kmer_set = GetTestData<K, N, KeyType>();

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

TEST(Unitigs, GetUnitigsCanonicalRandom) {
  const int K = 9;
  const int N = 10;
  using KeyType = uint8_t;

  KmerSet<K, N, KeyType> kmer_set = GetTestData<K, N, KeyType>();

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

TEST(Unitigs, GetSPSSCanonicalRandom) {
  const int K = 9;
  const int N = 10;
  using KeyType = uint8_t;

  KmerSet<K, N, KeyType> kmer_set = GetTestData<K, N, KeyType>();

  std::vector<std::string> spss = GetSPSSCanonical(kmer_set, false, 1);

  KmerSet<K, N, KeyType> reconstructed;

  for (const std::string& s : spss) {
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

TEST(Unitigs, GetSPSSCanonicalFastRandom) {
  const int K = 9;
  const int N = 10;
  using KeyType = uint8_t;

  KmerSet<K, N, KeyType> kmer_set = GetTestData<K, N, KeyType>();

  std::vector<std::string> spss = GetSPSSCanonical(kmer_set, true, 1);

  KmerSet<K, N, KeyType> reconstructed;

  for (const std::string& s : spss) {
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

#include "core/spss.h"

#include <cstdint>
#include <string>
#include <vector>

#include "absl/random/random.h"
#include "core/kmer.h"
#include "core/kmer_set.h"
#include "gtest/gtest.h"
#include "random.h"

TEST(SPSS, Complement) { ASSERT_EQ(internal::Complement("ACGTT"), "AACGT"); }

TEST(SPSS, GetUnitigsRandom) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  KmerSet<K, N, KeyType> kmer_set = GetRandomKmerSet<K, N, KeyType>(
      absl::Uniform(absl::IntervalClosed, absl::InsecureBitGen(), 1, 1 << 16),
      false);

  std::vector<std::string> unitigs = GetUnitigs(kmer_set, n_workers);

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

TEST(SPSS, GetUnitigsCanonicalRandom) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  KmerSet<K, N, KeyType> kmer_set = GetRandomKmerSet<K, N, KeyType>(
      absl::Uniform(absl::IntervalClosed, absl::InsecureBitGen(), 1, 1 << 16),
      true);

  std::vector<std::string> unitigs = GetUnitigsCanonical(kmer_set, n_workers);

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

  ASSERT_TRUE(kmer_set.Equals(reconstructed, n_workers));
}

TEST(SPSS, GetSPSSRandom) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  KmerSet<K, N, KeyType> kmer_set = GetRandomKmerSet<K, N, KeyType>(
      absl::Uniform(absl::IntervalClosed, absl::InsecureBitGen(), 1, 1 << 16),
      false);

  std::vector<std::string> spss = GetSPSS(kmer_set, n_workers);

  KmerSet<K, N, KeyType> reconstructed;

  for (const std::string& s : spss) {
    const int n = s.length();
    ASSERT_TRUE(n >= K);

    for (int i = 0; i < n - K + 1; i++) {
      const Kmer<K> kmer = Kmer<K>(s.substr(i, K));
      ASSERT_FALSE(reconstructed.Contains(kmer));
      reconstructed.Add(kmer);
    }
  }

  ASSERT_TRUE(kmer_set.Equals(reconstructed, n_workers));
}

TEST(SPSS, GetSPSSCanonicalRandom) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  KmerSet<K, N, KeyType> kmer_set = GetRandomKmerSet<K, N, KeyType>(
      absl::Uniform(absl::IntervalClosed, absl::InsecureBitGen(), 1, 1 << 16),
      true);

  std::vector<std::string> spss = GetSPSSCanonical(kmer_set, false, n_workers);

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

  ASSERT_TRUE(kmer_set.Equals(reconstructed, n_workers));
}

TEST(SPSS, GetSPSSCanonicalFastRandom) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  KmerSet<K, N, KeyType> kmer_set = GetRandomKmerSet<K, N, KeyType>(
      absl::Uniform(absl::IntervalClosed, absl::InsecureBitGen(), 1, 1 << 16),
      true);

  std::vector<std::string> spss = GetSPSSCanonical(kmer_set, true, n_workers);

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

  ASSERT_TRUE(kmer_set.Equals(reconstructed, n_workers));
}

TEST(SPSS, GetKmerSetFromSPSSRandom) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  KmerSet<K, N, KeyType> kmer_set = GetRandomKmerSet<K, N, KeyType>(
      absl::Uniform(absl::IntervalClosed, absl::InsecureBitGen(), 1, 1 << 16),
      false);

  std::vector<std::string> spss = GetSPSS(kmer_set, n_workers);

  KmerSet<K, N, KeyType> reconstructed =
      GetKmerSetFromSPSS<K, N, KeyType>(spss, false, n_workers);

  ASSERT_TRUE(kmer_set.Equals(reconstructed, n_workers));
}

TEST(SPSS, GetKmerSetFromSPSSCanonicalRandom) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  KmerSet<K, N, KeyType> kmer_set = GetRandomKmerSet<K, N, KeyType>(
      absl::Uniform(absl::IntervalClosed, absl::InsecureBitGen(), 1, 1 << 16),
      true);

  std::vector<std::string> spss = GetSPSSCanonical(kmer_set, true, n_workers);

  KmerSet<K, N, KeyType> reconstructed =
      GetKmerSetFromSPSS<K, N, KeyType>(spss, true, n_workers);

  ASSERT_TRUE(kmer_set.Equals(reconstructed, n_workers));
}

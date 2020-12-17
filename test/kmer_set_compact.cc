#include "core/kmer_set_compact.h"

#include <string>

#include "absl/random/random.h"
#include "core/kmer.h"
#include "core/kmer_set.h"
#include "gtest/gtest.h"

template <int K, int N, typename KeyType>
KmerSet<K, N, KeyType> GetTestData(bool canonical) {
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
      if (canonical) {
        kmer_set.Add(Kmer<K>(s.substr(j, K)).Canonical());
      } else {
        kmer_set.Add(Kmer<K>(s.substr(j, K)));
      }
    }
  }

  return kmer_set;
}

TEST(KmerSetCompact, FromAndToKmerSetRandom) {
  const int K = 9;
  const int N = 10;
  using KeyType = uint8_t;

  // Canonical
  {
    KmerSet<K, N, KeyType> kmer_set = GetTestData<K, N, KeyType>(true);

    KmerSetCompact<K, N, KeyType> compressed =
        KmerSetCompact<K, N, KeyType>::FromKmerSet(kmer_set, true, false, 1);

    KmerSet<K, N, KeyType> decompressed = compressed.ToKmerSet(true, 1);

    ASSERT_TRUE(kmer_set.Equals(decompressed, 1));
  }

  // Non-canonical
  {
    {
      KmerSet<K, N, KeyType> kmer_set = GetTestData<K, N, KeyType>(false);

      KmerSetCompact<K, N, KeyType> compressed =
          KmerSetCompact<K, N, KeyType>::FromKmerSet(kmer_set, false, false, 1);

      KmerSet<K, N, KeyType> decompressed = compressed.ToKmerSet(false, 1);

      ASSERT_TRUE(kmer_set.Equals(decompressed, 1));
    }
  }
}

TEST(KmerSetCompact, FromAndToKmerSetFastRandom) {
  const int K = 9;
  const int N = 10;
  using KeyType = uint8_t;

  // Canonical
  {
    KmerSet<K, N, KeyType> kmer_set = GetTestData<K, N, KeyType>(true);

    KmerSetCompact<K, N, KeyType> compressed =
        KmerSetCompact<K, N, KeyType>::FromKmerSet(kmer_set, true, true, 1);

    KmerSet<K, N, KeyType> decompressed = compressed.ToKmerSet(true, 1);

    ASSERT_TRUE(kmer_set.Equals(decompressed, 1));
  }

  // Non-canonical
  {
    KmerSet<K, N, KeyType> kmer_set = GetTestData<K, N, KeyType>(false);

    KmerSetCompact<K, N, KeyType> compressed =
        KmerSetCompact<K, N, KeyType>::FromKmerSet(kmer_set, false, true, 1);

    KmerSet<K, N, KeyType> decompressed = compressed.ToKmerSet(false, 1);

    ASSERT_TRUE(kmer_set.Equals(decompressed, 1));
  }
}

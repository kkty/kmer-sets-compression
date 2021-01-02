#ifndef RANDOM_H_
#define RANDOM_H_

#include <string>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/random/random.h"
#include "core/kmer.h"
#include "core/kmer_set.h"
#include "core/kmer_set_immutable.h"
#include "core/kmer_set_set.h"

// Returns a randomly-generated kmer.
template <int K>
Kmer<K> GetRandomKmer() {
  absl::InsecureBitGen bitgen;
  Kmer<K> kmer;

  for (int i = 0; i < K; i++) {
    switch (absl::Uniform(bitgen, 0, 4)) {
      case 0:
        kmer.Set(i, 'A');
        break;
      case 1:
        kmer.Set(i, 'G');
        break;
      case 2:
        kmer.Set(i, 'C');
        break;
      case 3:
        kmer.Set(i, 'T');
        break;
    }
  }

  return kmer;
}

// Returns n randomly-generated kmers.
template <int K>
std::vector<Kmer<K>> GetRandomKmers(int n) {
  absl::flat_hash_set<Kmer<K>> kmers;

  while (static_cast<int>(kmers.size()) < n) {
    kmers.insert(GetRandomKmer<K>());
  }

  return std::vector<Kmer<K>>(kmers.begin(), kmers.end());
}

// Randomly constructs a KmerSet with n kmers.
template <int K, int N, typename KeyType>
KmerSet<K, N, KeyType> GetRandomKmerSet(int n, bool canonical) {
  absl::InsecureBitGen bitgen;

  absl::flat_hash_set<Kmer<K>> kmers;

  while (true) {
    std::string s;

    for (int j = 0; j < absl::Uniform(absl::IntervalClosed, bitgen, 1, 100);
         j++) {
      s += GetRandomKmer<K>().String();
    }

    // Creates a loop.
    if (absl::Uniform(bitgen, 0, 2) == 0) {
      s += s;
    }

    for (int j = 0; j < static_cast<int>(s.size()) - K + 1; j++) {
      Kmer<K> kmer(s.substr(j, K));
      if (canonical) kmer = kmer.Canonical();
      kmers.insert(kmer);
      if (static_cast<int>(kmers.size()) == n) break;
    }

    if (static_cast<int>(kmers.size()) == n) break;
  }

  KmerSet<K, N, KeyType> kmer_set;
  for (const Kmer<K>& kmer : kmers) kmer_set.Add(kmer);

  return kmer_set;
}

// Constructs n "KmerSetImmutable"s each with m kmers.
template <int K, int N, typename KeyType>
std::vector<KmerSetImmutable<K, N, KeyType>> GetRandomKmerSetsImmutable(
    int n, int m, int n_workers) {
  std::vector<KmerSetImmutable<K, N, KeyType>> kmer_sets_immutable(n);

  for (int i = 0; i < n; i++) {
    kmer_sets_immutable[i] = KmerSetImmutable<K, N, KeyType>(
        GetRandomKmerSet<K, N, KeyType>(m, true), n_workers);
  }

  return kmer_sets_immutable;
}

// Constructs a KmerSetSet with n kmer sets, each with m kmers.
template <int K, int N, typename KeyType>
KmerSetSet<K, N, KeyType> GetRandomKmerSetSet(int n, int m, int n_workers) {
  return KmerSetSet<K, N, KeyType>(
      GetRandomKmerSetsImmutable<K, N, KeyType>(n, m, n_workers), true,
      n_workers);
}

#endif
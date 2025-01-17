#ifndef RANDOM_H_
#define RANDOM_H_

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/random/random.h"
#include "core/kmer.h"
#include "core/kmer_counter.h"
#include "core/kmer_set.h"
#include "core/kmer_set_set.h"

// Returns a randomly-generated kmer.
template <int K>
Kmer<K> GetRandomKmer() {
  return Kmer<K>(absl::Uniform<std::uint64_t>(
      absl::IntervalClosed, absl::InsecureBitGen(), 0,
      std::numeric_limits<std::uint64_t>::max() >> (64 - K * 2)));
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

// Returns a random read.
template <int K>
std::string GetRandomRead() {
  std::string s;
  absl::InsecureBitGen bitgen;

  for (int j = 0; j < absl::Uniform(absl::IntervalClosed, bitgen, 1, 100);
       j++) {
    s += GetRandomKmer<K>().String();
  }

  // Creates a loop.
  if (absl::Uniform(bitgen, 0, 2) == 0) {
    s += s;
  }

  return s;
}

template <int K, int N, typename KeyType>
KmerCounter<K, N, KeyType> GetRandomKmerCounter(int n, bool canonical) {
  int counter = 0;

  KmerCounter<K, N, KeyType> kmer_counter;

  while (true) {
    std::string s = GetRandomRead<K>();

    for (int j = 0; j < static_cast<int>(s.size()) - K + 1; j++) {
      Kmer<K> kmer(s.substr(j, K));
      if (canonical) kmer = kmer.Canonical();
      kmer_counter.Add(kmer, 1);
      counter += 1;

      if (counter == n) break;
    }

    if (counter == n) break;
  }

  return kmer_counter;
}

// Randomly constructs a KmerSet with n kmers.
template <int K, int N, typename KeyType>
KmerSet<K, N, KeyType> GetRandomKmerSet(int n, bool canonical) {
  absl::InsecureBitGen bitgen;

  absl::flat_hash_set<Kmer<K>> kmers;

  while (true) {
    std::string s = GetRandomRead<K>();

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

template <int K, int N, typename KeyType>
KmerSetCompact<K, N, KeyType> GetRandomKmerSetCompact(int n, bool canonical,
                                                      int n_workers) {
  return KmerSetCompact<K, N, KeyType>::FromKmerSet(
      GetRandomKmerSet<K, N, KeyType>(n, canonical), canonical, true,
      n_workers);
}

// Constructs n KmerSetCompact each with m kmers.
template <int K, int N, typename KeyType>
std::vector<KmerSetCompact<K, N, KeyType>> GetRandomKmerSetsCompact(
    int n, int m, bool canonical, int n_workers) {
  std::vector<KmerSetCompact<K, N, KeyType>> kmer_sets_compact(n);

  for (int i = 0; i < n; i++) {
    kmer_sets_compact[i] =
        GetRandomKmerSetCompact<K, N, KeyType>(m, canonical, n_workers);
  }

  return kmer_sets_compact;
}

// Constructs a KmerSetSet with n kmer sets, each with m kmers.
template <int K, int N, typename KeyType>
KmerSetSet<K, N, KeyType> GetRandomKmerSetSet(int n, int m, bool canonical,
                                              int n_workers) {
  return KmerSetSet<K, N, KeyType>(
      GetRandomKmerSetsCompact<K, N, KeyType>(n, m, canonical, n_workers),
      canonical, n_workers);
}

#endif
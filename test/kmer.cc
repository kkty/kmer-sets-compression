#include "core/kmer.h"

#include <string>
#include <unordered_map>

#include "absl/container/flat_hash_map.h"
#include "gtest/gtest.h"

TEST(kmer, String) {
  const std::string s = "AGCTG";
  Kmer<5> kmer(s);
  ASSERT_EQ(kmer.String(), s);
}

TEST(kmer, Canonical) {
  ASSERT_EQ(Kmer<5>("AAAAT").Canonical().String(), "AAAAT");
  ASSERT_EQ(Kmer<5>("TTTTA").Canonical().String(), "TAAAA");
  ASSERT_EQ(Kmer<5>("CCCCG").Canonical().String(), "CCCCG");
  ASSERT_EQ(Kmer<5>("GGGGC").Canonical().String(), "GCCCC");
}

TEST(kmer, Complement) {
  Kmer<5> kmer("AGCTA");
  ASSERT_EQ(kmer.Complement().String(), "TAGCT");
}

TEST(kmer, Next) {
  Kmer<5> kmer("AGCTG");
  ASSERT_EQ(kmer.Next('C').String(), "GCTGC");
}

TEST(kmer, Prev) {
  Kmer<5> kmer("AGCTG");
  ASSERT_EQ(kmer.Prev('C').String(), "CAGCT");
}

TEST(kmer, UnorderedMapKey) {
  Kmer<5> kmer("AGCTG");
  std::unordered_map<Kmer<5>, int> m;
  m[kmer] = 1;
  ASSERT_EQ(m[kmer], 1);
}

TEST(kmer, FlatHashMapKey) {
  Kmer<5> kmer("AGCTG");
  absl::flat_hash_map<Kmer<5>, int> m;
  m[kmer] = 1;
  ASSERT_EQ(m[kmer], 1);
}

#include "core/kmer.h"

#include <string>
#include <unordered_map>

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

TEST(kmer, Hash) {
  std::hash<Kmer<5>> h;
  ASSERT_EQ(h(Kmer<5>("ACGTA")), h(Kmer<5>("ACGTA")));
}

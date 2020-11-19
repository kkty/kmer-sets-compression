#include "core/io.h"

#include <string>
#include <vector>

#include "gtest/gtest.h"

std::vector<std::string> GetTestData() {
  std::vector<std::string> lines;

  for (int i = 0; i < 1000000; i++) {
    if (i % 2)
      lines.push_back("foo");
    else
      lines.push_back("bar");
  }

  return lines;
}

TEST(io, WriteAndRead) {
  const std::string file_name = "/tmp/data";
  const std::vector<std::string> lines = GetTestData();

  WriteLines(file_name, lines);

  ASSERT_EQ(ReadLines(file_name), lines);
}

TEST(io, WriteAndReadCompressed) {
  const std::string file_name = "/tmp/data";
  const std::vector<std::string> lines = GetTestData();

  WriteLines(file_name, "gzip -c", lines);

  ASSERT_EQ(ReadLines(file_name, "gzip -d -c"), lines);
}

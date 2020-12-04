#include "core/io.h"

#include <string>
#include <vector>

#include "absl/status/statusor.h"
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

  {
    absl::Status status = WriteLines(file_name, lines);
    ASSERT_TRUE(status.ok());
  }

  {
    absl::StatusOr<std::vector<std::string>> statusor = ReadLines(file_name);

    ASSERT_TRUE(statusor.ok());
    ASSERT_EQ(statusor.value(), lines);
  }
}

TEST(io, WriteAndReadCompressed) {
  const std::string file_name = "/tmp/data.gz";
  const std::vector<std::string> lines = GetTestData();

  {
    absl::Status status = WriteLines(file_name, "gzip -c", lines);
    ASSERT_TRUE(status.ok());
  }

  {
    absl::StatusOr<std::vector<std::string>> statusor =
        ReadLines(file_name, "gzip -d -c");

    ASSERT_TRUE(statusor.ok());
    ASSERT_EQ(statusor.value(), lines);
  }
}

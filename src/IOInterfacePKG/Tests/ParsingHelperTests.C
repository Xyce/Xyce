#include <string>
#include <iostream>

#include <gtest/gtest.h>

#include <Xyce_config.h>
#include <N_IO_ParsingHelpers.h>

namespace {
TEST(IO_Parsing_Helpers, isAbsolutePath)
{
  ASSERT_EQ(true, Xyce::IO::isAbsolutePath("/path"));
  ASSERT_EQ(false, Xyce::IO::isAbsolutePath("path"));
  ASSERT_EQ(false, Xyce::IO::isAbsolutePath("./path"));
  ASSERT_EQ(false, Xyce::IO::isAbsolutePath("../path"));

  #ifdef HAVE_WINDOWS_H
    ASSERT_EQ(true, Xyce::IO::isAbsolutePath("\\path"));
    ASSERT_EQ(true, Xyce::IO::isAbsolutePath("\\\\path"));
    ASSERT_EQ(true, Xyce::IO::isAbsolutePath("C:\\path"));
    ASSERT_EQ(true, Xyce::IO::isAbsolutePath("c:\\path"));
  #else
    ASSERT_EQ(false, Xyce::IO::isAbsolutePath("\\path"));
    ASSERT_EQ(false, Xyce::IO::isAbsolutePath("\\\\path"));
    ASSERT_EQ(false, Xyce::IO::isAbsolutePath("C:\\path"));
    ASSERT_EQ(false, Xyce::IO::isAbsolutePath("c:\\path"));
  #endif

  // These are relative paths even on Windows
  ASSERT_EQ(false, Xyce::IO::isAbsolutePath("C:path"));
  ASSERT_EQ(false, Xyce::IO::isAbsolutePath("c:path"));
}

TEST(IO_Parsing_Helpers, hasWinDriveLetter)
{
  #ifdef HAVE_WINDOWS_H
  ASSERT_EQ(true,Xyce::IO::hasWinDriveLetter("c:"));
  ASSERT_EQ(true,Xyce::IO::hasWinDriveLetter("C:"));
  #else
  ASSERT_EQ(false,Xyce::IO::hasWinDriveLetter("c:"));
  ASSERT_EQ(false,Xyce::IO::hasWinDriveLetter("C:"));
  #endif
}

} // namespace

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

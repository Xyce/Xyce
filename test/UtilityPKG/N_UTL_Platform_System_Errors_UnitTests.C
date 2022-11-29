#include <gtest/gtest.h>
#include <N_UTL_Platform.h>
#include <string>

namespace
{
    TEST(UTL_Platform_System_Errors, hostname_Fail_test)
    {
      std::string error;
      std::string expectedHostName(EXPECTED_HOSTNAME);
      std::string testHostName = Xyce::hostname(error);

      #if defined(HAVE_WINDOWS_H)
      // This test is expected to fail on Windows right now; the mock function for getComputerNameEx is not implemented
      // EXPECT_TRUE(testHostName.empty());
      // EXPECT_EQ(error, "hostname: There was an error retriving the local DNS hostname, the function GetComputerNameEx set Last Error to 234");

      #elif defined(HAVE_GETHOSTNAME)
      EXPECT_TRUE(testHostName.empty());
      EXPECT_EQ(error, "hostname: The function gethostname() failed with errno 14");
      #endif
    }
}

int main(int argc, char *argv[]) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
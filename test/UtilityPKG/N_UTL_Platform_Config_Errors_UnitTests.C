#include <gtest/gtest.h>
#include <N_UTL_Platform.h>
#include <string>

namespace
{
    TEST(UTL_Platform_Config_Errors, hostname_function_not_found_test)
    {
      std::string error;
      std::string expectedHostName(EXPECTED_HOSTNAME);
      std::string testHostName = Xyce::hostname(error);
      EXPECT_TRUE(testHostName.empty());
      EXPECT_EQ(error, "hostname: No platform specification defined to retrieve the local hostname");
    }
}

int main(int argc, char *argv[]) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
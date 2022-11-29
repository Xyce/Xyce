#include <gtest/gtest.h>
#include <N_UTL_Platform.h>

#include <string>

namespace
{
    TEST(UTL_Platform, hostname_test)
    {
        std::string expectedHostName(EXPECTED_HOSTNAME);
        std::string testHostName = Xyce::hostname();
        ASSERT_EQ(testHostName, expectedHostName);
        ASSERT_NE(testHostName, "");
    }

    TEST(UTL_Platform, hostname_WithOutputArg_test)
    {
        std::string error;
        std::string expectedHostName(EXPECTED_HOSTNAME);
        std::string testHostName = Xyce::hostname(error);
        EXPECT_EQ(testHostName, expectedHostName);
        EXPECT_TRUE(error.empty());
    }
}

int main(int argc, char *argv[]) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
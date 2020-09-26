
// this is from  https://youtu.be/nbFXI9SDfbk
// 
// can build with: 
//
//    clang++ Test.C -I/opt/local/include -L/opt/local/lib -lgtest -lgtest_main -pthread
//

#include <iostream>
#include <gtest/gtest.h>


TEST ( TestName, Subtest_1) 
{
  EXPECT_EQ(1,2);
  //ASSERT_EQ(1,2);

  //ASSERT_FALSE(1 == 1);

  std::cout << "After assertion" <<std::endl;
}

TEST ( TestName2, Subtest_1)
{
  ASSERT_TRUE(1 == 2);
}



// EXPECT_EQ, ASSERT_EQ - non-fatal vs fatal failures

int main (int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}


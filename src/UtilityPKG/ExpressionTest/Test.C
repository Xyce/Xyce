//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
//   Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
//   NTESS, the U.S. Government retains certain rights in this software.
//
//   This file is part of the Xyce(TM) Parallel Electrical Simulator.
//
//   Xyce(TM) is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   Xyce(TM) is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with Xyce(TM).
//   If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------

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


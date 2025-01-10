//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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

#include <complex>
#include <algorithm>
#include <iterator>

#include "ast.h"
#include <newExpression.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_ExtendedString.h>

class testExpressionGroup : public Xyce::Util::baseExpressionGroup 
{ 
  public: 
    testExpressionGroup () : Xyce::Util::baseExpressionGroup()  {}; 
    ~testExpressionGroup () {}; 
}; 

//-------------------------------------------------------------------------------
// test values of binary operators 
//
#define PARSER_SIMPLE_TEST_MACRO(NAME,SUBNAME,STREXP, CPPEXP) \
TEST ( NAME, SUBNAME ) \
{ \
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() ); \
  Xyce::Util::newExpression testExpression(std::string(STREXP), testGroup); \
  testExpression.lexAndParseExpression(); \
  Xyce::Util::newExpression copyExpression(testExpression); \
  double result(0.0); \
  copyExpression.evaluateFunction(result); \
  EXPECT_EQ( (result-(CPPEXP)), 0.0); \
} 

// number by itself
TEST ( Double_Parser_Test, numval) 
{ 
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("1.0"), testGroup); 
  testExpression.lexAndParseExpression(); 
  Xyce::Util::newExpression copyExpression(testExpression); 
  double result(1.0); 
  copyExpression.evaluateFunction(result); 
  EXPECT_EQ( (result-(1.0)), 0.0); 
}

// binary operators
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, binaryAdd, "1.0+2.0", (1.0+2.0) )
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, binaryMinus, "1.0-2.0", (1.0-2.0) )
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, binaryMul, "4.0*3.0", (4.0*3.0) )
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, binaryDiv, "1.0/4.0", (1.0/4.0) )

int main (int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}



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
#include "newExpression.h"
#include "N_UTL_BreakPoint.h"

#include "value.h"

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


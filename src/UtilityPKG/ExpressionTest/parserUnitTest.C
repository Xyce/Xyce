
// this is from  https://youtu.be/nbFXI9SDfbk
//
// can build with:
//
//    clang++ Test.C -I/opt/local/include -L/opt/local/lib -lgtest -lgtest_main -pthread
//

#include <fstream>
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

#define OUTPUT_MACRO(NAME,SUBNAME) \
{ \
  char filename[ ] = "parserUnitTest.out"; \
  std::fstream outputFile; \
  outputFile.open(filename,  std::fstream::in | std::fstream::out | std::fstream::app ); \
  EXPECT_TRUE(outputFile.is_open()); \
  outputFile << std::endl << "TEST (" #NAME"," #SUBNAME ")" <<std::endl; \
  testExpression.dumpParseTree(outputFile); \
  outputFile.close();  } { \
  char filename[ ] = "parserUnitTest_codeGen.C"; \
  std::fstream outputFile; \
  outputFile.open(filename,  std::fstream::in | std::fstream::out | std::fstream::app ); \
  EXPECT_TRUE(outputFile.is_open()); \
  outputFile << std::endl << "// TEST (" #NAME"," #SUBNAME ")" <<std::endl; \
  outputFile << "{" <<std::endl; \
  testExpression.codeGen(outputFile); \
  outputFile << "}" <<std::endl; \
  outputFile.close(); \
}

#define OUTPUT_MACRO2(NAME,SUBNAME,EXPRNAME) \
{ \
  char filename[ ] = "parserUnitTest.out"; \
  std::fstream outputFile; \
  outputFile.open(filename,  std::fstream::in | std::fstream::out | std::fstream::app ); \
  EXPECT_TRUE(outputFile.is_open()); \
  outputFile << std::endl << "TEST (" #NAME"," #SUBNAME ")" <<std::endl; \
  EXPRNAME.dumpParseTree(outputFile); \
  outputFile.close();  } { \
  char filename[ ] = "parserUnitTest_codeGen.C"; \
  std::fstream outputFile; \
  outputFile.open(filename,  std::fstream::in | std::fstream::out | std::fstream::app ); \
  EXPECT_TRUE(outputFile.is_open()); \
  outputFile << std::endl << "// TEST (" #NAME"," #SUBNAME ")" <<std::endl; \
  outputFile << "{" <<std::endl; \
  EXPRNAME.codeGen(outputFile); \
  outputFile << "}" <<std::endl; \
  outputFile.close(); \
}

#define OUTPUT_MACRO3(NAME,SUBNAME) \
{ \
  char filename[ ] = "parserUnitTest.out"; \
  std::fstream outputFile; \
  outputFile.open(filename,  std::fstream::in | std::fstream::out | std::fstream::app ); \
  EXPECT_TRUE(outputFile.is_open()); \
  outputFile << std::endl << "TEST (" #NAME"," #SUBNAME ")" <<std::endl; \
  testExpression->dumpParseTree(outputFile); \
  outputFile.close();  } { \
  char filename[ ] = "parserUnitTest_codeGen.C"; \
  std::fstream outputFile; \
  outputFile.open(filename,  std::fstream::in | std::fstream::out | std::fstream::app ); \
  EXPECT_TRUE(outputFile.is_open()); \
  outputFile << std::endl << "// TEST (" #NAME"," #SUBNAME ")" <<std::endl; \
  outputFile << "{" <<std::endl; \
  testExpression->codeGen(outputFile); \
  outputFile << "}" <<std::endl; \
  outputFile.close(); \
}

//-------------------------------------------------------------------------------
// test values of binary operators
//
#define PARSER_SIMPLE_TEST_MACRO(NAME,SUBNAME,STREXP, CPPEXP) \
TEST ( NAME, SUBNAME ) \
{ \
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() ); \
  Xyce::Util::newExpression testExpression(std::string(STREXP), testGroup); \
  testExpression.lexAndParseExpression(); \
  double result(0.0); \
  testExpression.evaluateFunction(result); \
  EXPECT_NEAR( (result-(CPPEXP)), 0.0,1.0e-15); \
  Xyce::Util::newExpression copyExpression(testExpression); \
  copyExpression.evaluateFunction(result); \
  EXPECT_NEAR( (result-(CPPEXP)), 0.0,1.0e-15); \
  Xyce::Util::newExpression assignExpression; \
  assignExpression = testExpression; \
  assignExpression.evaluateFunction(result); \
  EXPECT_NEAR( (result-(CPPEXP)), 0.0,1.0e-15); \
  { char filename[ ] = "parserUnitTest.out"; \
  std::fstream outputFile; \
  outputFile.open(filename,  std::fstream::in | std::fstream::out | std::fstream::app ); \
  EXPECT_TRUE(outputFile.is_open()); \
  outputFile << std::endl << "TEST (" #NAME"," #SUBNAME ")" <<std::endl; \
  testExpression.dumpParseTree(outputFile); \
  outputFile.close(); } { \
  char filename[ ] = "parserUnitTest_codeGen.C"; \
  std::fstream outputFile; \
  outputFile.open(filename,  std::fstream::in | std::fstream::out | std::fstream::app ); \
  EXPECT_TRUE(outputFile.is_open()); \
  outputFile << std::endl << "// TEST (" #NAME"," #SUBNAME ")" <<std::endl; \
  outputFile << "{" <<std::endl; \
  testExpression.codeGen(outputFile); \
  outputFile << "}" <<std::endl; \
  outputFile.close();} \
}

// number by itself
TEST ( Double_Parser_Test, numval)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("1.0"), testGroup);
  testExpression.lexAndParseExpression();
  double result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_DOUBLE_EQ( (result-(1.0)), 0.0);

  Xyce::Util::newExpression copyExpression(testExpression); 
  copyExpression.evaluateFunction(result); 
  EXPECT_DOUBLE_EQ( (result-(1.0)), 0.0);

  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 
  assignExpression.evaluateFunction(result); 
  EXPECT_DOUBLE_EQ( (result-(1.0)), 0.0);
}

TEST ( Double_Parser_Test, numval2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("10.61E6"), testGroup);
  testExpression.lexAndParseExpression();
  double result(10.61E6);
  testExpression.evaluateFunction(result);
  EXPECT_DOUBLE_EQ( (result-(10.61E6)), 0.0);

  Xyce::Util::newExpression copyExpression(testExpression); 
  copyExpression.evaluateFunction(result); 
  EXPECT_DOUBLE_EQ( (result-(10.61E6)), 0.0);

  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 
  assignExpression.evaluateFunction(result); 
  EXPECT_DOUBLE_EQ( (result-(10.61E6)), 0.0);
}

// these next 3 tests are for parameters that happen to have the same name as single-character operators
TEST ( Double_Parser_Test, singleParam_R)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("R"), testGroup);
  testExpression.lexAndParseExpression();
  double result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_DOUBLE_EQ( result, 0.0);
}

TEST ( Double_Parser_Test, singleParam_M)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("M"), testGroup);
  testExpression.lexAndParseExpression();
  double result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_DOUBLE_EQ( result, 0.0);
}

TEST ( Double_Parser_Test, singleParam_P)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("P"), testGroup);
  testExpression.lexAndParseExpression();
  double result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_DOUBLE_EQ( result, 0.0);
}

// binary operators
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, binaryAdd, "1.0+2.0", (1.0+2.0) )
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, binaryMinus, "1.0-2.0", (1.0-2.0) )
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, binaryMul, "4.0*3.0", (4.0*3.0) )
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, binaryDiv, "1.0/4.0", (1.0/4.0) )

// simple precendence testing (via binary operators):
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, precedence1, "3.0*2.0+4.0", (3.0*2.0+4.0) )
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, precedence2, "5.0+4.0/2.0", (5.0+4.0/2.0) )
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, precedence3, "4.0*6.0/2.0", (4.0*6.0/2.0) )
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, precedence4, "4.0*(6.0/2.0)", (4.0*(6.0/2.0)) )
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, precedence5, "1.0/4.0*10.0", (1.0/4.0*10.0) )
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, precedence6, "1.0/(4.0*10.0)", (1.0/(4.0*10.0)) )

PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, unaryPlus, "+2.0", (2.0) )
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, unaryMinus, "-2.0", (-2.0) )

PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, phase, "P(1.0)", std::arg(1.0) )
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, real1, "Re(1.0)", 1.0 )
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, real2, "R(1.0)", 1.0 )
// Im cannot work along with IM for imaginary current.  disabling
//PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, imag1, "Im(1.0)", 0.0 )
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, imag2, "Img(1.0)", 0.0 )

PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, int1, "int(11.2423)", 11)
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, int2, "int(-11.2423)", -11)

PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, uramp1, "uramp(11.2423)", 11.2423)
PARSER_SIMPLE_TEST_MACRO ( Double_Parser_Test, uramp2, "uramp(-11.2423)", 0)

// std library functions
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, sqrt,  "sqrt(4.0)",  std::sqrt(4.0))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, exp,   "exp(0.5)", std::exp(0.5))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, abs,   "abs(-0.5)", std::abs(-0.5))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, sin,   "sin(0.5)", std::sin(0.5))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, cos,   "cos(0.5)", std::cos(0.5))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, acos,  "acos(0.5)", std::acos(0.5))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, acosh, "acosh(1.5)", std::acosh(1.5))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, asin,  "asin(0.5)", std::asin(0.5))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, asinh, "asinh(0.5)", std::asinh(0.5))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, atan,  "atan(0.5)", std::atan(0.5))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, atanh, "atanh(0.5)", std::atanh(0.5))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, cosh,  "cosh(0.5)", std::cosh(0.5))
//PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, log,   "log(0.5)", std::log(0.5))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, ln,   "ln(0.5)", std::log(0.5))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, log,   "log(0.5)", std::log10(0.5))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, log10, "log10(0.5)", std::log10(0.5))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, sinh,  "sinh(0.5)", std::sinh(0.5))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, tan,   "tan(0.5)", std::tan(0.5))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, tanh,  "tanh(0.5)", std::tanh(0.5))

PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, pow1,  "pow(2.0,3.0)", std::pow(2.0,3.0))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, pow2,  "2.0**3.0", std::pow(2.0,3.0))

PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, pwrs1,  "pwrs(2.0,3.0)", std::pow(2.0,3.0))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, pwrs2,  "pwrs(0.0,3.0)", 0.0);
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, pwrs3,  "pwrs(-2.0,3.0)", -std::pow(2.0,3.0))

PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, sign1,  "sign(-25,10.25)", 25)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, sign2,  "sign(15,-10.25)", -15)

// Hspice only:
//PARSER_SIMPLE_TEST_MACRO(Double_Parser_UnaryFunc_Test, pow3,  "2.0^3.0", std::pow(2.0,3.0))

// lower case metrix prefix/suffix tests
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, tera,  "3.0t", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, giga,  "5.0g", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, kilo,  "7.0k", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, mega,  "2.0meg", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, mega2,  "4.0x", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, micro,  "2.0u", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, nano,  "9.0n", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, pico,  "6.0p", 6.0e-12)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, femto,  "6.0f", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, mil,  "2.0mil", 2.0*(25.4e-6) )

PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, teraSec,  "3.0ts", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, gigaSec,  "5.0gs", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, kiloSec,  "7.0ks", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, megaSec,  "2.0megs", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, mega2Sec,  "4.0xs", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, microSec,  "2.0us", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, nanoSec,  "9.0ns", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, picoSec,  "6.0ps", 6.0e-12)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, femtoSec,  "6.0fs", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, milSec,  "2.0mils", 2.0*(25.4e-6) )

PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, teraVolt,  "3.0tv", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, gigaVolt,  "5.0gv", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, kiloVolt,  "7.0kv", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, megaVolt,  "2.0megv", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, mega2Volt,  "4.0xv", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, microVolt,  "2.0uv", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, nanoVolt,  "9.0nv", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, picoVolt,  "6.0pv", 6.0e-12)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, femtoVolt,  "6.0fv", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, milVolt,  "2.0milv", 2.0*(25.4e-6) )

PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, teraAmp,  "3.0ta", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, gigaAmp,  "5.0ga", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, kiloAmp,  "7.0ka", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, megaAmp,  "2.0mega", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, mega2Amp,  "4.0xa", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, microAmp,  "2.0ua", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, nanoAmp,  "9.0na", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, picoAmp,  "6.0pa", 6.0e-12)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, femtoAmp,  "6.0fa", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, milAmp,  "2.0mila", 2.0*(25.4e-6) )

// unit suffixes, which should be ignored, lower case
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, volt,  "7.0v", 7.0 )
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, amp,  "6.0a", 6.0 )
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, sec,  "5.0s", 5.0 )

// upper case metrix prefix/suffix tests
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, TERA,  "3.0T", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, GIGA,  "5.0G", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, KILO,  "7.0K", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, MEGA,  "2.0MEG", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, MEGA2,  "4.0X", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, MICRO,  "2.0U", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, NANO,  "9.0N", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, PICO,  "6.0P", 6.0e-12)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, FEMTO,  "6.0F", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, MIL,  "2.0MIL", 2.0*(25.4e-6) )

PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, TERASEC,  "3.0TS", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, GIGASEC,  "5.0GS", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, KILOSEC,  "7.0KS", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, MEGASEC,  "2.0MEGS", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, MEGA2SEC,  "4.0XS", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, MICROSEC,  "2.0US", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, NANOSEC,  "9.0NS", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, PICOSEC,  "6.0PS", 6.0E-12)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, FEMTOSEC,  "6.0FS", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, MILSEC,  "2.0MILS", 2.0*(25.4e-6) )

PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, TERAVOLT,  "3.0TV", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, GIGAVOLT,  "5.0GV", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, KILOVOLT,  "7.0KV", 7.0E+3)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, MEGAVOLT,  "2.0MEGV", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, MEGA2VOLT,  "4.0XV", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, MICROVOLT,  "2.0UV", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, NANOVOLT,  "9.0NV", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, PICOVOLT,  "6.0PV", 6.0e-12)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, FEMTOVOLT,  "6.0FV", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, MILVOLT,  "2.0MILV", 2.0*(25.4e-6) )

PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, TERAAMP,  "3.0TA", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, GIGAAMP,  "5.0GA", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, KILOAMP,  "7.0KA", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, MEGAAMP,  "2.0MEGA", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, MEGA2AMP,  "4.0XA", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, MICROAMP,  "2.0UA", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, NANOAMP,  "9.0NA", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, PICOAMP,  "6.0PA", 6.0e-12)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, FEMTOAMP,  "6.0FA", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, MILAMP,  "2.0MILA", 2.0*(25.4e-6) )

// unit suffixes, which should be ignored, upper case
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, VOLT,  "4.0V", 4.0 )
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, AMP,  "3.0A", 3.0 )
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, SEC,  "2.0S", 2.0 )

// lower case metrix prefix/suffix tests
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, sin_tera,  "sin(3.0t)", std::sin(3.0e+12))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, sin_giga,  "sin(5.0g)", std::sin(5.0e+9))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, sin_kilo,  "sin(7.0k)", std::sin(7.0e+3))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, sin_mega,  "sin(2.0meg)", std::sin(2.0e+6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, sin_mega2,  "sin(4.0x)", std::sin(4.0e+6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, sin_micro,  "sin(2.0u)", std::sin(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, sin_nano,  "sin(9.0n)", std::sin(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, sin_pico,  "sin(6.0p)", std::sin(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, sin_femto,  "sin(6.0f)", std::sin(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, sin_mil,  "sin(2.0mil)", std::sin(2.0*(25.4e-6)))

PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, exp_teraSec,  "exp(3.0e-12ts)", std::exp(3.0))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, exp_gigaSec,  "exp(5.0e-9gs)", std::exp(5.0))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, exp_kiloSec,  "exp(7.0e-3ks)", std::exp(7.0))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, exp_megaSec,  "exp(2.0e-6megs)", std::exp(2.0))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, exp_mega2Sec,  "exp(4.0e-6xs)", std::exp(4.0))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, exp_microSec,  "exp(2.0us)", std::exp(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, exp_nanoSec,  "exp(9.0ns)", std::exp(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, exp_picoSec,  "exp(6.0ps)", std::exp(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, exp_femtoSec,  "exp(6.0fs)", std::exp(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, exp_milSec,  "exp(2.0mils)", std::exp(2.0*(25.4e-6)))

PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, cos_teraVolt,  "cos(3.0tv)", std::cos(3.0e+12))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, cos_gigaVolt,  "cos(5.0gv)", std::cos(5.0e+9))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, cos_kiloVolt,  "cos(7.0kv)", std::cos(7.0e+3))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, cos_megaVolt,  "cos(2.0megv)", std::cos(2.0e+6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, cos_mega2Volt,  "cos(4.0xv)", std::cos(4.0e+6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, cos_microVolt,  "cos(2.0uv)", std::cos(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, cos_nanoVolt,  "cos(9.0nv)", std::cos(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, cos_picoVolt,  "cos(6.0pv)", std::cos(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, cos_femtoVolt,  "cos(6.0fv)", std::cos(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, cos_milVolt,  "cos(2.0milv)", std::cos(2.0*(25.4e-6)))

PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, tan_teraAmp,  "tan(3.0ta)", std::tan(3.0e+12))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, tan_gigaAmp,  "tan(5.0ga)", std::tan(5.0e+9))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, tan_kiloAmp,  "tan(7.0ka)", std::tan(7.0e+3))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, tan_megaAmp,  "tan(2.0mega)", std::tan(2.0e+6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, tan_mega2Amp,  "tan(4.0xa)", std::tan(4.0e+6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, tan_microAmp,  "tan(2.0ua)", std::tan(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, tan_nanoAmp,  "tan(9.0na)", std::tan(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, tan_picoAmp,  "tan(6.0pa)", std::tan(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, tan_femtoAmp,  "tan(6.0fa)", std::tan(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, tan_milAmp,  "tan(2.0mila)", std::tan(2.0*(25.4e-6)))

// unit suffixes, which should be ignored, lower case
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, sin_volt,  "sin(2.0v)", std::sin(2.0))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, sin_amp,  "sin(3.0a)", std::sin(3.0))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, sin_sec,  "sin(4.0s)", std::sin(4.0))

// upper case
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, SIN_TERA,  "SIN(3.0T)", std::sin(3.0e+12))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, SIN_GIGA,  "SIN(5.0G)", std::sin(5.0e+9))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, SIN_KILO,  "SIN(7.0K)", std::sin(7.0e+3))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, SIN_MEGA,  "SIN(2.0MEG)", std::sin(2.0e+6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, SIN_MEGA2,  "SIN(4.0X)", std::sin(4.0e+6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, SIN_MICRO,  "SIN(2.0U)", std::sin(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, SIN_NANO,  "SIN(9.0N)", std::sin(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, SIN_PICO,  "SIN(6.0P)", std::sin(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, SIN_FEMTO,  "SIN(6.0F)", std::sin(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, SIN_MIL,  "SIN(2.0MIL)", std::sin(2.0*(25.4e-6)))

PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, EXP_TERASEC,  "EXP(3.0E-12TS)", std::exp(3.0))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, EXP_GIGASEC,  "EXP(5.0E-9GS)", std::exp(5.0))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, EXP_KILOSEC,  "EXP(7.0E-3KS)", std::exp(7.0))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, EXP_MEGASEC,  "EXP(2.0E-6MEGS)", std::exp(2.0))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, EXP_MEGA2SEC,  "EXP(4.0E-6XS)", std::exp(4.0))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, EXP_MICROSEC,  "EXP(2.0US)", std::exp(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, EXP_NANOSEC,  "EXP(9.0NS)", std::exp(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, EXP_PICOSEC,  "EXP(6.0PS)", std::exp(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, EXP_FEMTOSEC,  "EXP(6.0FS)", std::exp(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, EXP_MILSEC,  "EXP(2.0MILS)", std::exp(2.0*(25.4e-6)))

PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, COS_TERAVOLT,  "COS(3.0TV)", std::cos(3.0e+12))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, COS_GIGAVOLT,  "COS(5.0GV)", std::cos(5.0e+9))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, COS_KILOVOLT,  "COS(7.0KV)", std::cos(7.0e+3))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, COS_MEGAVOLT,  "COS(2.0MEGV)", std::cos(2.0e+6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, COS_MEGA2VOLT,  "COS(4.0XV)", std::cos(4.0e+6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, COS_MICROVOLT,  "COS(2.0UV)", std::cos(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, COS_NANOVOLT,  "COS(9.0NV)", std::cos(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, COS_PICOVOLT,  "COS(6.0PV)", std::cos(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, COS_FEMTOVOLT,  "COS(6.0FV)", std::cos(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, COS_MILVOLT,  "COS(2.0MILV)", std::cos(2.0*(25.4e-6)))

PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, TAN_TERAAMP,  "TAN(3.0TA)", std::tan(3.0e+12))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, TAN_GIGAAMP,  "TAN(5.0GA)", std::tan(5.0e+9))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, TAN_KILOAMP,  "TAN(7.0KA)", std::tan(7.0e+3))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, TAN_MEGAAMP,  "TAN(2.0MEGA)", std::tan(2.0e+6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, TAN_MEGA2AMP,  "TAN(4.0XA)", std::tan(4.0e+6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, TAN_MICROAMP,  "TAN(2.0UA)", std::tan(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, TAN_NANOAMP,  "TAN(9.0NA)", std::tan(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, TAN_PICOAMP,  "TAN(6.0PA)", std::tan(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, TAN_FEMTOAMP,  "TAN(6.0FA)", std::tan(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, TAN_MILAMP,  "TAN(2.0MILA)", std::tan(2.0*(25.4e-6)))

// unit suffixes, which should be ignored, upper case
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, SIN_VOLT,  "SIN(5.0V)", std::sin(5.0))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, SIN_AMP,  "SIN(6.0A)", std::sin(6.0))
PARSER_SIMPLE_TEST_MACRO(Double_Parser_Suffix_Test, SIN_SEC,  "SIN(7.0S)", std::sin(7.0))


// source functions:
class timeDepExpressionGroup : public Xyce::Util::baseExpressionGroup
{
  public:
    timeDepExpressionGroup () : Xyce::Util::baseExpressionGroup(), time(0.0), freq(0.0)  {};
    ~timeDepExpressionGroup () {};
    virtual double getTime() { return time; };
    virtual double getFreq() { return freq; };
    void setTime(double t) { time = t; };
    void setFreq(double f) { freq = f; };
    double time;
    double freq;
};

TEST ( Double_Parser_SourceFunc_Test, pulse)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;
  Xyce::Util::newExpression testExpression(std::string("spice_pulse(0.0,1.0,0.0,10e-6,10e-6,0.1e-6,20.1e-6)"), testGroup);
  testExpression.lexAndParseExpression();
  double result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( (result-(0.0)), 0.0);

  double tr=10.0e-6, pw=0.1e-6;
  timeDepGroup->setTime(tr+0.5*pw);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( (result-(1.0)), 0.0);

  Xyce::Util::newExpression copyExpression(testExpression); 
  copyExpression.evaluateFunction(result); 
  EXPECT_EQ( (result-(1.0)), 0.0);

  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 
  assignExpression.evaluateFunction(result); 
  EXPECT_EQ( (result-(1.0)), 0.0);

  OUTPUT_MACRO(Double_Parser_SourceFunc_Test,pulse)
}

TEST ( Double_Parser_SourceFunc_Test, pulse_breakpoints)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;
  Xyce::Util::newExpression testExpression(std::string("spice_pulse(0.0,1.0,0.0,10e-6,10e-6,0.1e-6,20.1e-6)"), testGroup);
  testExpression.lexAndParseExpression();
  double result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( (result-(0.0)), 0.0);

  double tr=10.0e-6, pw=0.1e-6;
  timeDepGroup->setTime(tr+0.5*pw);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( (result-(1.0)), 0.0);

  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  testExpression.getBreakPoints(breakPointTimes);

  OUTPUT_MACRO(Double_Parser_SourceFunc_Test,pulse_breakpoints)
}

TEST ( Double_Parser_SourceFunc_Test, sin)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;
  Xyce::Util::newExpression testExpression(std::string("spice_sin(1.65,1.65,10000,0,0,-90)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  int numpoints=100;
  double v0(1.65), va(1.65), freq(10000), td(0.0), theta(0.0), phase(-90),time(0.0);
  double dt=(1.0/freq)*(1.0/static_cast<double>(numpoints));
  std::vector<double> refRes(numpoints), result(numpoints);
  std::vector<double> copyResult(numpoints), assignResult(numpoints);
  for (int ii=0;ii<numpoints;ii++,time+=dt)
  {
    timeDepGroup->setTime(time); 
    testExpression.evaluateFunction(result[ii]);
    copyExpression.evaluateFunction(copyResult[ii]);
    assignExpression.evaluateFunction(assignResult[ii]);
    refRes[ii] = v0 + va * std::sin(2.0*M_PI*((freq)*time + (phase)/360)) * std::exp( -(time*(theta)));
  }
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( copyResult, refRes);
  EXPECT_EQ( assignResult, refRes);

  OUTPUT_MACRO(Double_Parser_SourceFunc_Test,sin)
}

// This test is taken from the sources.cir Xyce regression test.  This is the "B2" source in that test.
// The test runs with Xyce if I comment out the B2 source.  So I created this unit test to track down the problem.
//
// The problem is that (as of this writing 4/23/2020) is that I wrote the spice 
// sources in the new expression library to have a hardwired number of comma-separated arguments.
// For the spice_sin source, the hardwired number is 6.   In the regression test, only 5 arguments
// are given.  if I add the 6th (for phase) it runs fine and appears to give the right answer.
//
// I know how to modify the spice sources to have a variable number of args.  That is now on my to-do list.
//
// This means that I need to modify the spice sources so they can have a variable number of arguments.
// I know how to do this;  I did it with POLY for example.
TEST ( Double_Parser_SourceFunc_Test, sin2)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;

  // Note that the original test used a parameter named "v2" for va. So, doing that here.
  //Xyce::Util::newExpression testExpression(std::string("spice_sin(v0,v2,frequency,td,theta)-2mv"), testGroup);  // original, without phase
  Xyce::Util::newExpression testExpression(std::string("spice_sin(v0,v2,frequency,td,theta,phase)-2mv"), testGroup); // works, but must have all 6 args
  testExpression.lexAndParseExpression();

  // setup the parameters
  Teuchos::RCP<Xyce::Util::newExpression> v0Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("-0.5"), testGroup));
  v0Expression->lexAndParseExpression();
  std::string v0Name = "V0";
  testExpression.attachParameterNode(v0Name,v0Expression);

  Teuchos::RCP<Xyce::Util::newExpression> v2Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("2"), testGroup));
  v2Expression->lexAndParseExpression();
  std::string v2Name = "V2";
  testExpression.attachParameterNode(v2Name,v2Expression);

  Teuchos::RCP<Xyce::Util::newExpression> freqExpression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("3.4e+7"), testGroup));
  freqExpression->lexAndParseExpression();
  std::string freqName = "FREQUENCY";
  testExpression.attachParameterNode(freqName,freqExpression);

  Teuchos::RCP<Xyce::Util::newExpression> tdExpression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.5ns"), testGroup));
  tdExpression->lexAndParseExpression();
  std::string tdName = "TD";
  testExpression.attachParameterNode(tdName,tdExpression);

  Teuchos::RCP<Xyce::Util::newExpression> thetaExpression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.1"), testGroup));
  thetaExpression->lexAndParseExpression();
  std::string thetaName = "THETA";
  testExpression.attachParameterNode(thetaName,thetaExpression);

  // not part of the original:
  Teuchos::RCP<Xyce::Util::newExpression> phaseExpression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.0"), testGroup));
  phaseExpression->lexAndParseExpression();
  std::string phaseName = "PHASE";
  testExpression.attachParameterNode(phaseName,phaseExpression);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  int numpoints=100;
  double v0(-0.5), va(2.0), freq(3.4e+7), td(0.5e-9), theta(0.1), phase(0.0),time(0.0);

  double maxtime = 100e-9;
  double dt=maxtime * (1.0/static_cast<double>(numpoints));

  std::vector<double> refRes(numpoints), result(numpoints);
  std::vector<double> copyResult(numpoints), assignResult(numpoints);
  for (int ii=0;ii<numpoints;ii++,time+=dt)
  {
    timeDepGroup->setTime(time); 
    testExpression.evaluateFunction(result[ii]);
    copyExpression.evaluateFunction(copyResult[ii]);
    assignExpression.evaluateFunction(assignResult[ii]);
    double time2=time-td;
    if (time2<=0.0)
    {
      refRes[ii] = v0 + va * std::sin (2.0*M_PI*((std::real(phase))/360)) ;
    }
    else
    {
      refRes[ii] = v0 + va * std::sin(2.0*M_PI*((freq)*(time2) + (phase)/360)) * std::exp( -((time2)*(theta)));
    }
    refRes[ii] -= 2.0e-3;
  }

  EXPECT_EQ( result, refRes);
  EXPECT_EQ( copyResult, refRes);
  EXPECT_EQ( assignResult, refRes);

  OUTPUT_MACRO(Double_Parser_SourceFunc_Test,sin)
}

TEST ( Double_Parser_SourceFunc_Test, exp)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = timeDepGroup;
  Xyce::Util::newExpression
    testExpression(std::string("spice_exp(1.1,2.0,2e-9,15e-9,5e-9,30e-9)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  int numpoints=100;
  double v1(1.1), v2(2.0), td1(2e-9), tau1(15e-9), td2(5e-9), tau2(30e-9), time(0.0);
  double dt=2*td2/static_cast<double>(numpoints);
  std::vector<double> refRes(numpoints), result(numpoints);
  std::vector<double> copyResult(numpoints), assignResult(numpoints);
  for (int ii=0;ii<numpoints;ii++,time+=dt)
  {
    timeDepGroup->setTime(time); 
    testExpression.evaluateFunction(result[ii]);
    copyExpression.evaluateFunction(copyResult[ii]);
    assignExpression.evaluateFunction(assignResult[ii]);
    if (time <= td1)  refRes[ii] = v1;
    else if (time <= td2 && time > td1) refRes[ii] = v1 + (v2-v1)*(1.0-std::exp(-(time-td1)/tau1));
    else refRes[ii] = v1 + (v2-v1)*(1.0-std::exp(-(time-td1)/tau1)) + (v1-v2)*(1.0-std::exp(-(time-td2)/tau2)) ;
  }
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( copyResult, refRes);
  EXPECT_EQ( assignResult, refRes);

  OUTPUT_MACRO(Double_Parser_SourceFunc_Test,exp)
}

TEST ( Double_Parser_SourceFunc_Test, sffm)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = timeDepGroup;
  Xyce::Util::newExpression
    testExpression(std::string("spice_sffm(-0.5,2.0,100e6,0.3,2.1e6)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  int numpoints=100;
  double v0(-0.5), va(2.0), fc(100e6), mdi(0.3), fs(2.1e6), time(0.0);
  double dt=(1.0/2.1e6) /  static_cast<double>(numpoints);
  std::vector<double> refRes(numpoints), result(numpoints);
  std::vector<double> copyResult(numpoints), assignResult(numpoints);
  for (int ii=0;ii<numpoints;ii++,time+=dt)
  {
    timeDepGroup->setTime(time); 
    testExpression.evaluateFunction(result[ii]);
    copyExpression.evaluateFunction(copyResult[ii]);
    assignExpression.evaluateFunction(assignResult[ii]);
    refRes[ii] = v0 + va * sin((2 * M_PI * fc * time) + mdi * sin (2 * M_PI * fs * time));
  }
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( copyResult, refRes);
  EXPECT_EQ( assignResult, refRes);
  OUTPUT_MACRO(Double_Parser_SourceFunc_Test,sffm)
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class solnExpressionGroup : public Xyce::Util::baseExpressionGroup
{
  public:
    solnExpressionGroup () :
      Xyce::Util::baseExpressionGroup(), Aval(0.0), Bval(0.0), Cval(0.0), R1val(0.0), 
         VBval(0.0), VCval(0.0), VEval(0.0), VLPval(0.0), VLNval(0.0)
  {};
    ~solnExpressionGroup () {};

  virtual bool getSolutionVal(const std::string & nodeName, double & retval )
  {
    std::string tmp = nodeName; Xyce::Util::toLower(tmp);
    if (tmp==std::string("a")) { retval = Aval; return true; }
    else if (tmp==std::string("b")) { retval = Bval; return true; }
    else if (tmp==std::string("c")) { retval = Cval; return true; }
    else if (tmp==std::string("r1")) { retval = R1val; return true; }

    else if (tmp==std::string("vb")) { retval = VBval; return true; }
    else if (tmp==std::string("vc")) { retval = VCval; return true; }
    else if (tmp==std::string("ve")) { retval = VEval; return true; }
    else if (tmp==std::string("vlp")) { retval = VLPval; return true; }
    else if (tmp==std::string("vln")) { retval = VLNval; return true; }
    else { return 0.0; return false; }
  }

  void setSoln(const std::string & nodeName, double val)
  {
    std::string tmp = nodeName; Xyce::Util::toLower(tmp);
    if (tmp==std::string("a")) { Aval = val; }
    else if (tmp==std::string("b")) { Bval = val; }
    else if (tmp==std::string("c")) { Cval = val; }
    else if (tmp==std::string("r1")) { R1val = val; }

    else if (tmp==std::string("vb")) { VBval = val; }
    else if (tmp==std::string("vc")) { VCval = val; }
    else if (tmp==std::string("ve")) { VEval = val; }
    else if (tmp==std::string("vlp")) { VLPval = val; }
    else if (tmp==std::string("vln")) { VLNval = val; }
  }

#if 0
  void addParam (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
  {
    std::string lowerName = name;
    Xyce::Util::toLower(lowerName);

    parameters_[lowerName] = exp;
  };

  bool getParam       (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
  {
    bool retval=true;
    std::string lowerName = name;
    Xyce::Util::toLower(lowerName);
    if (parameters_.find(lowerName) != parameters_.end()) { exp = parameters_[lowerName]; }
    else { retval = false; }
    return retval;
  }
#endif

  private:
    //std::unordered_map <std::string, Teuchos::RCP<Xyce::Util::newExpression> >  parameters_;

  double Aval, Bval, Cval, R1val;
  double VBval, VCval, VEval, VLPval, VLNval;
};

TEST ( Double_Parser_VoltSoln_Test, test0)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("V(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=3.0;
  double refRes = Aval;
  solnGroup->setSoln(std::string("A"),Aval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(Double_Parser_VoltSoln_Test, test0)
}

TEST ( Double_Parser_VoltSoln_Test, test1)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("10.0*V(A)+1.0"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=3.0;
  double refRes = 10*Aval+1.0;
  solnGroup->setSoln(std::string("A"),Aval);
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(Double_Parser_VoltSoln_Test, test1)
}

TEST ( Double_Parser_VoltSoln_Test, test2)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("12.0*V(A,B)+7.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=6.3, Bval=2.1;
  double refRes = 12.0*(Aval-Bval)+7.5;
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("B"),Bval);
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(Double_Parser_VoltSoln_Test, test2)
}

// testing out the complex operators.  As this file tests the use case where the 
// default type is "double" rather than "std::complex<double>" these operators don't 
// do much.  But they do something, and need to be included.
TEST ( Double_Parser_VoltSoln_Test, vr_test0)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("Vr (A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=3.0;
  double refRes = std::real(Aval);
  solnGroup->setSoln(std::string("A"),Aval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(Double_Parser_VoltSoln_Test, test0)
}

TEST ( Double_Parser_VoltSoln_Test, vi_test0)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("vI (A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=3.0;
  double refRes = std::imag(3.0);
  solnGroup->setSoln(std::string("A"),Aval);

#if 0
  testExpression.dumpParseTree(std::cout);
#endif

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(Double_Parser_VoltSoln_Test, test0)
}

TEST ( Double_Parser_VoltDeriv_Test, test1)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("12.3*V(A)*V(B)+7.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=6.3, Bval=2.1;
  double refRes = 12.3*Aval*Bval+7.5;
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("B"),Bval);
  std::vector<double> refDer;
  refDer.push_back(12.3*Bval);
  refDer.push_back(12.3*Aval);
  std::vector<double> derivs;

  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
  OUTPUT_MACRO(Double_Parser_VoltDeriv_Test, test1)
}

TEST ( Double_Parser_VoltDeriv_Test, test2)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*V(A)**2.0+7.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=6.3;
  double refRes = 20.0*Aval*Aval+7.5;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<double> refDer;
  refDer.push_back(20.0* ( (2.0/Aval)*std::pow(Aval,2.0)) );

  std::vector<double> derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0]-refDer[0], 0.0);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0]-refDer[0], 0.0);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0]-refDer[0], 0.0);
  OUTPUT_MACRO(Double_Parser_VoltDeriv_Test, test2)
}

TEST ( Double_Parser_VoltDeriv_Test, test3)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*V(A)*V(A)+7.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=6.3;
  double refRes = 20.0*Aval*Aval+7.5;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<double> refDer;
  refDer.push_back(20.0*Aval + 20.0*Aval);
  std::vector<double> derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
  OUTPUT_MACRO(Double_Parser_VoltDeriv_Test, test3)
}

TEST ( Double_Parser_VoltDeriv_Test, test4)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*V(A)*V(A)+7.5*V(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=6.3;
  double refRes = 20.0*Aval*Aval+7.5*Aval;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<double> refDer;
  refDer.push_back( 20.0*Aval + 20.0*Aval + 7.5 );
  std::vector<double> derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
  OUTPUT_MACRO(Double_Parser_VoltDeriv_Test, test4)
}

TEST ( Double_Parser_VoltDeriv_Test, test5)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*(V(A)**3.0)+7.5*V(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=6.3;
  double refRes = 20.0*std::pow(Aval,3.0)+7.5*Aval;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<double> refDer;
  refDer.push_back( 20.0*(3.0/Aval)*std::pow(Aval,3.0)+7.5 );
  std::vector<double> derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
  OUTPUT_MACRO(Double_Parser_VoltDeriv_Test, test5)
}

TEST ( Double_Parser_VoltDeriv_Test, test6)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*V(A)*V(A)*V(A)+7.5*V(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=6.3;
  double refRes = 20.0*std::pow(Aval,3.0)+7.5*Aval;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<double> refDer;
  refDer.push_back( 20.0*(3.0/Aval)*std::pow(Aval,3.0)+7.5 );
  std::vector<double> derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
  OUTPUT_MACRO(Double_Parser_VoltDeriv_Test, test6)
}


TEST ( Double_Parser_CurrSoln_Test, test1)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("17.2*I(R1)+8.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, R1val=3.0;
  double refRes = 17.2*R1val+8.5;
  solnGroup->setSoln(std::string("R1"),R1val);
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(Double_Parser_CurrSoln_Test, test1)
}

TEST ( Double_Parser_CurrDeriv_Test, test1)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("17.2*I(R1)+8.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, R1val=3.0;
  double refRes = 17.2*R1val+8.5;
  solnGroup->setSoln(std::string("R1"),R1val);

  std::vector<double> refDer;
  refDer.push_back( 17.2 );
  std::vector<double> derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
  OUTPUT_MACRO(Double_Parser_CurrDeriv_Test, test1)
}

TEST ( Double_Parser_CurrDeriv_Test, test2)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*I(R1)*I(R1)*I(R1)+7.5*I(R1)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, R1val=6.3;
  double refRes = 20.0*std::pow(R1val,3.0)+7.5*R1val;
  solnGroup->setSoln(std::string("R1"),R1val);

  std::vector<double> refDer;
  refDer.push_back( 20.0*(3.0/R1val)*std::pow(R1val,3.0)+7.5 );
  std::vector<double> derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
  OUTPUT_MACRO(Double_Parser_CurrDeriv_Test, test2)
}


//-------------------------------------------------------------------------------
class internalDevExpressionGroup : public Xyce::Util::baseExpressionGroup
{
  public:
    internalDevExpressionGroup () : Xyce::Util::baseExpressionGroup() {};
    ~internalDevExpressionGroup () {};

  void setInternalDeviceVar (const std::string & name, double val)
  {
    std::string lowerName = name;
    Xyce::Util::toLower(lowerName);
    invernalVars_[lowerName] = val;
  };

  bool getInternalDeviceVar       (const std::string & name, double & val)
  {
    bool retval=true;
    std::string lowerName = name;
    Xyce::Util::toLower(lowerName);
    if (invernalVars_.find(lowerName) != invernalVars_.end()) { val = invernalVars_[lowerName]; }
    else { retval = false; }
    return retval;
  }

  private:
    std::unordered_map <std::string, double> invernalVars_;
};

TEST ( Double_Parser_InternalDeviceVariable_Test, test1)
{
  Teuchos::RCP<internalDevExpressionGroup> intVarGroup = Teuchos::rcp(new internalDevExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = intVarGroup;
  Xyce::Util::newExpression testExpression(std::string("17.2*N(M3:GM)+8.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, M3GMval=3.0;
  double refRes = 17.2*M3GMval+8.5;
  intVarGroup->setInternalDeviceVar(std::string("M3:GM"),M3GMval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);

  //copyExpression.evaluateFunction(result); 
  //assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
 
  OUTPUT_MACRO(Double_Parser_CurrSoln_Test, test1)
}


//-------------------------------------------------------------------------------
// .func tests
class testExpressionGroupWithFuncSupport : public Xyce::Util::baseExpressionGroup
{
  public:
    testExpressionGroupWithFuncSupport () : Xyce::Util::baseExpressionGroup()  {};
    ~testExpressionGroupWithFuncSupport () {};

#if 0
    void addFunction (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
    {
      std::string lowerName = name;
      Xyce::Util::toLower(lowerName);

      functions_[lowerName] = exp;
    };

    bool getFunction (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
    {
      bool retval=true;

      std::string lowerName = name;
      Xyce::Util::toLower(lowerName);

      if (functions_.find(lowerName) != functions_.end()) { exp = functions_[lowerName]; }
      else { retval = false; }

      return retval;
    }
#endif

  private:
    //std::unordered_map <std::string, Teuchos::RCP<Xyce::Util::newExpression>  >  functions_;
};

//-------------------------------------------------------------------------------
TEST ( Double_Parser_Func_Test, test1)
{
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression will use the .func f1.
  Xyce::Util::newExpression testExpression(std::string("F1(2,3)"), testGroup);
  testExpression.lexAndParseExpression();

  // this expression is the RHS of a .func statement:  .func F1(A,B) {A+B}
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("A+B"), testGroup) );

  // I originally had this set up so that the calling code would manually set the 
  // vector of prototype function arguments, as well as the name of the function 
  // itself.  But in a code like Xyce, that isn't how it is likely to work. The 
  // function prototype F1(A,B) has to be parsed, and the appropriate information 
  // pulled out of it.  In Xyce, the old expression library is used to parse the 
  // prototype(LHS), so attempting same here.
  Xyce::Util::newExpression f1_LHS (std::string("F1(A,B)"), testGroup);
  f1_LHS.lexAndParseExpression();

  std::vector<std::string> f1ArgStrings ;
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);
  f1Expression->setFunctionArgStringVec (f1ArgStrings);
  // during lex/parse, this vector of arg strings will be compared to any
  // param classes.  If it finds them, then they will be placed in the
  // functionArgOpVec object, which is used below, in the call to "setFuncArgs".
  f1Expression->lexAndParseExpression();

  // now parse the function name from the prototype
  //std::string f1Name = "F1";
  std::string f1Name; 
  f1_LHS.getFuncPrototypeName(f1Name);

  testExpression.attachFunctionNode(f1Name, f1Expression);

#if 0
  testExpression.dumpParseTree(std::cout);
#endif

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result;
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, 5.0 );
  OUTPUT_MACRO(Double_Parser_Func_Test, test1)
}

//-------------------------------------------------------------------------------
TEST ( Double_Parser_Func_Test, test1_multipleLexParse)
{
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression will use the .func f1.
  Xyce::Util::newExpression testExpression(std::string("F1(2,3)"), testGroup);
  testExpression.lexAndParseExpression();

  // this expression is the RHS of a .func statement:  .func F1(A,B) {A+B}
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("A+B"), testGroup) );

  // I originally had this set up so that the calling code would manually set the 
  // vector of prototype function arguments, as well as the name of the function 
  // itself.  But in a code like Xyce, that isn't how it is likely to work. The 
  // function prototype F1(A,B) has to be parsed, and the appropriate information 
  // pulled out of it.  In Xyce, the old expression library is used to parse the 
  // prototype(LHS), so attempting same here.
  Xyce::Util::newExpression f1_LHS (std::string("F1(A,B)"), testGroup);
  f1_LHS.lexAndParseExpression();

  std::vector<std::string> f1ArgStrings ;
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);

// testing to see if I can do this twice, once b4 and once after the setting of function args.  (maybe not, maybe needs a clear)
  f1Expression->lexAndParseExpression();

  f1Expression->setFunctionArgStringVec (f1ArgStrings);
  // during lex/parse, this vector of arg strings will be compared to any
  // param classes.  If it finds them, then they will be placed in the
  // functionArgOpVec object, which is used below, in the call to "setFuncArgs".
  f1Expression->lexAndParseExpression();

  // now parse the function name from the prototype
  //std::string f1Name = "F1";
  std::string f1Name; 
  f1_LHS.getFuncPrototypeName(f1Name);

  testExpression.attachFunctionNode(f1Name, f1Expression);

#if 0
  testExpression.dumpParseTree(std::cout);
#endif

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result;
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, 5.0 );
  OUTPUT_MACRO(Double_Parser_Func_Test, test1)
}



// tests are taken from the "ternary_precedence.cir" Xyce regression test
TEST ( Double_Parser_ternary_precedence, simple)
{
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression is the RHS of a .func statement:  .func simple(X) {x>0?2*x :0}
  Teuchos::RCP<Xyce::Util::newExpression> simpleExpression
    = Teuchos::rcp(new Xyce::Util::newExpression(std::string("x>0?2*x :0"), testGroup) );

  std::vector<std::string> simpleArgStrings;
  Xyce::Util::newExpression simple_LHS (std::string("simple(X)"), testGroup);
  simple_LHS.lexAndParseExpression();
  simple_LHS.getFuncPrototypeArgStrings(simpleArgStrings);
  simpleExpression->setFunctionArgStringVec (simpleArgStrings);
  simpleExpression->lexAndParseExpression();

  // these expressions uses the .func simple.
  {
    Xyce::Util::newExpression simpleTrue(std::string("simple(4)"), testGroup);
    simpleTrue.lexAndParseExpression();

    std::string simpleName;
    simple_LHS.getFuncPrototypeName(simpleName);

    simpleTrue.attachFunctionNode(simpleName ,  simpleExpression);

    Xyce::Util::newExpression copySimpleTrue(simpleTrue); 
    Xyce::Util::newExpression assignSimpleTrue; 
    assignSimpleTrue = simpleTrue; 

    double result;
    simpleTrue.evaluateFunction(result);       EXPECT_EQ( result, 8.0 );
    copySimpleTrue.evaluateFunction(result);   EXPECT_EQ( result, 8.0 );
    assignSimpleTrue.evaluateFunction(result); EXPECT_EQ( result, 8.0 );
    OUTPUT_MACRO2(Double_Parser_ternary_precedence, simple, simpleTrue) 
  }
  {
    Xyce::Util::newExpression simpleFalse(std::string("simple(0)"), testGroup);
    simpleFalse.lexAndParseExpression();

    std::string simpleName;
    simple_LHS.getFuncPrototypeName(simpleName);
    simpleFalse.attachFunctionNode(simpleName ,  simpleExpression);

    Xyce::Util::newExpression copySimpleFalse(simpleFalse); 
    Xyce::Util::newExpression assignSimpleFalse; 
    assignSimpleFalse = simpleFalse; 

    double result;
    simpleFalse.evaluateFunction(result);       EXPECT_EQ( result, 0.0 );
    copySimpleFalse.evaluateFunction(result);   EXPECT_EQ( result, 0.0 );
    assignSimpleFalse.evaluateFunction(result); EXPECT_EQ( result, 0.0 );
    OUTPUT_MACRO2(Double_Parser_ternary_precedence, simple, simpleFalse) 
  }
}

TEST ( Double_Parser_ternary_precedence, precplus)
{
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression is the RHS of a .func statement:  .func precplus(X) {1+x>0?2*x :0+2}
  Teuchos::RCP<Xyce::Util::newExpression> precplusExpression
    = Teuchos::rcp(new Xyce::Util::newExpression(std::string("1+x>0?2*x :0+2"), testGroup) );
  std::vector<std::string> precplusArgStrings; 

  Xyce::Util::newExpression precplus_LHS (std::string("precplus(X)"), testGroup);
  precplus_LHS.lexAndParseExpression();
  precplus_LHS.getFuncPrototypeArgStrings(precplusArgStrings);
  precplusExpression->setFunctionArgStringVec (precplusArgStrings);
  precplusExpression->lexAndParseExpression();

  {
    Xyce::Util::newExpression precplusTrue(std::string("precplus(4)"), testGroup);
    precplusTrue.lexAndParseExpression();
    std::string precplusName;
    precplus_LHS.getFuncPrototypeName(precplusName);
    precplusTrue.attachFunctionNode(precplusName ,  precplusExpression);

    Xyce::Util::newExpression copyPrecplusTrue(precplusTrue); 
    Xyce::Util::newExpression assignPrecplusTrue; 
    assignPrecplusTrue = precplusTrue; 

    double result;
    precplusTrue.evaluateFunction(result);       EXPECT_EQ( result, 8.0 );
    copyPrecplusTrue.evaluateFunction(result);   EXPECT_EQ( result, 8.0 );
    assignPrecplusTrue.evaluateFunction(result); EXPECT_EQ( result, 8.0 );
    OUTPUT_MACRO2(Double_Parser_ternary_precedence, precplus, precplusTrue) 
  }
  {
    Xyce::Util::newExpression precplusFalse(std::string("precplus(-4)"), testGroup);
    precplusFalse.lexAndParseExpression();
    std::string precplusName;
    precplus_LHS.getFuncPrototypeName(precplusName);
    precplusFalse.attachFunctionNode(precplusName ,  precplusExpression);

    Xyce::Util::newExpression copyPrecplusFalse(precplusFalse); 
    Xyce::Util::newExpression assignPrecplusFalse; 
    assignPrecplusFalse = precplusFalse; 

    double result;
    precplusFalse.evaluateFunction(result);       EXPECT_EQ( result, 2.0 );
    copyPrecplusFalse.evaluateFunction(result);   EXPECT_EQ( result, 2.0 );
    assignPrecplusFalse.evaluateFunction(result); EXPECT_EQ( result, 2.0 );
    OUTPUT_MACRO2(Double_Parser_ternary_precedence, precplus, precplusFalse) 
  }
}

TEST ( Double_Parser_ternary_precedence, precplusparen)
{
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression is the RHS of a .func statement:  .func precplusparen(X) {(1+x)>0?(2*x) :(0+2)}
  Teuchos::RCP<Xyce::Util::newExpression> precplusparenExpression
    = Teuchos::rcp(new Xyce::Util::newExpression(std::string("(1+x)>0?(2*x) :(0+2)"), testGroup));

  std::vector<std::string> precplusparenArgStrings;
  Xyce::Util::newExpression precplusparen_LHS (std::string("precplusparen(X)"), testGroup);
  precplusparen_LHS.lexAndParseExpression();
  precplusparen_LHS.getFuncPrototypeArgStrings(precplusparenArgStrings);
  precplusparenExpression->setFunctionArgStringVec ( precplusparenArgStrings );
  precplusparenExpression->lexAndParseExpression();

  {
    Xyce::Util::newExpression precplusparenTrue(std::string("precplusparen(4)"), testGroup);
    precplusparenTrue.lexAndParseExpression();
    std::string precplusparenName;
    precplusparen_LHS.getFuncPrototypeName(precplusparenName);
    precplusparenTrue.attachFunctionNode(precplusparenName ,  precplusparenExpression);

    Xyce::Util::newExpression copyPrecplusparenTrue(precplusparenTrue); 
    Xyce::Util::newExpression assignPrecplusparenTrue; 
    assignPrecplusparenTrue = precplusparenTrue; 

    double result;
    precplusparenTrue.evaluateFunction(result);       EXPECT_EQ( result, 8.0 );
    copyPrecplusparenTrue.evaluateFunction(result);   EXPECT_EQ( result, 8.0 );
    assignPrecplusparenTrue.evaluateFunction(result); EXPECT_EQ( result, 8.0 );
    OUTPUT_MACRO2(Double_Parser_ternary_precedence, precplusparen, precplusparenTrue) 
  }
  {
    Xyce::Util::newExpression precplusparenFalse(std::string("precplusparen(-4)"), testGroup);
    precplusparenFalse.lexAndParseExpression();
    std::string precplusparenName;
    precplusparen_LHS.getFuncPrototypeName(precplusparenName);
    precplusparenFalse.attachFunctionNode(precplusparenName ,  precplusparenExpression);

    Xyce::Util::newExpression copyPrecplusparenFalse(precplusparenFalse); 
    Xyce::Util::newExpression assignPrecplusparenFalse; 
    assignPrecplusparenFalse = precplusparenFalse; 

    double result;
    precplusparenFalse.evaluateFunction(result);       EXPECT_EQ( result, 2.0 );
    copyPrecplusparenFalse.evaluateFunction(result);   EXPECT_EQ( result, 2.0 );
    assignPrecplusparenFalse.evaluateFunction(result); EXPECT_EQ( result, 2.0 );
    OUTPUT_MACRO2(Double_Parser_ternary_precedence, precplusparen, precplusparenFalse) 
  }
}

TEST ( Double_Parser_ternary_precedence, simpleif)
{
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression is the RHS of a .func statement:  .func simpleif(X) {if(x>0,2*x,0)}
  Teuchos::RCP<Xyce::Util::newExpression> simpleifExpression 
     = Teuchos::rcp(new Xyce::Util::newExpression(std::string("if(x>0,2*x,0)"), testGroup) );

  std::vector<std::string> simpleifArgStrings;

  Xyce::Util::newExpression simpleif_LHS (std::string("simpleif(X)"), testGroup);
  simpleif_LHS.lexAndParseExpression();
  simpleif_LHS.getFuncPrototypeArgStrings(simpleifArgStrings);
  simpleifExpression->setFunctionArgStringVec ( simpleifArgStrings );
  simpleifExpression->lexAndParseExpression();

  // these expressions uses the .func simpleif.
  {
    Xyce::Util::newExpression simpleifTrue(std::string("simpleif(4)"), testGroup);
    simpleifTrue.lexAndParseExpression();
    std::string simpleifName;
    simpleif_LHS.getFuncPrototypeName(simpleifName);
    simpleifTrue.attachFunctionNode(simpleifName ,  simpleifExpression);

    Xyce::Util::newExpression copySimpleifTrue(simpleifTrue); 
    Xyce::Util::newExpression assignSimpleifTrue; 
    assignSimpleifTrue = simpleifTrue; 

    double result;
    simpleifTrue.evaluateFunction(result);       EXPECT_EQ( result, 8.0 );
    copySimpleifTrue.evaluateFunction(result);   EXPECT_EQ( result, 8.0 );
    assignSimpleifTrue.evaluateFunction(result); EXPECT_EQ( result, 8.0 );
    OUTPUT_MACRO2(Double_Parser_ternary_precedence, simpleif, simpleifTrue) 
  }
  {
    Xyce::Util::newExpression simpleifFalse(std::string("simpleif(0)"), testGroup);
    simpleifFalse.lexAndParseExpression();
    std::string simpleifName;
    simpleif_LHS.getFuncPrototypeName(simpleifName);
    simpleifFalse.attachFunctionNode(simpleifName ,  simpleifExpression);

    Xyce::Util::newExpression copySimpleifFalse(simpleifFalse); 
    Xyce::Util::newExpression assignSimpleifFalse; 
    assignSimpleifFalse = simpleifFalse; 

    double result;
    simpleifFalse.evaluateFunction(result);       EXPECT_EQ( result, 0.0 );
    copySimpleifFalse.evaluateFunction(result);   EXPECT_EQ( result, 0.0 );
    assignSimpleifFalse.evaluateFunction(result); EXPECT_EQ( result, 0.0 );
    OUTPUT_MACRO2(Double_Parser_ternary_precedence, simpleif, simpleifFalse) 
  }
}

TEST ( Double_Parser_ternary_precedence, precplusif)
{
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression is the RHS of a .func statement:  .func precplusif(X) {if((1+x>0),2*x,0+2)}
  Teuchos::RCP<Xyce::Util::newExpression> precplusifExpression
     = Teuchos::rcp(new Xyce::Util::newExpression(std::string("if((1+x>0),2*x,0+2)"), testGroup));

  std::vector<std::string> precplusifArgStrings;

  Xyce::Util::newExpression precplusif_LHS (std::string("precplusif(X)"), testGroup);
  precplusif_LHS.lexAndParseExpression();
  precplusif_LHS.getFuncPrototypeArgStrings(precplusifArgStrings);
  precplusifExpression->setFunctionArgStringVec ( precplusifArgStrings );
  precplusifExpression->lexAndParseExpression();
  {
    Xyce::Util::newExpression precplusifTrue(std::string("precplusif(4)"), testGroup);
    precplusifTrue.lexAndParseExpression();
    std::string precplusifName;
    precplusif_LHS.getFuncPrototypeName(precplusifName);
    precplusifTrue.attachFunctionNode(precplusifName ,  precplusifExpression);

    Xyce::Util::newExpression copyPrecplusifTrue(precplusifTrue); 
    Xyce::Util::newExpression assignPrecplusifTrue; 
    assignPrecplusifTrue = precplusifTrue; 

    double result;
    precplusifTrue.evaluateFunction(result);       EXPECT_EQ( result, 8.0 );
    copyPrecplusifTrue.evaluateFunction(result);   EXPECT_EQ( result, 8.0 );
    assignPrecplusifTrue.evaluateFunction(result); EXPECT_EQ( result, 8.0 );
    OUTPUT_MACRO2(Double_Parser_ternary_precedence, precplusif, precplusifTrue) 
  }
  {
    Xyce::Util::newExpression precplusifFalse(std::string("precplusif(-4)"), testGroup);
    precplusifFalse.lexAndParseExpression();
    std::string precplusifName;
    precplusif_LHS.getFuncPrototypeName(precplusifName);
    precplusifFalse.attachFunctionNode(precplusifName ,  precplusifExpression);

    Xyce::Util::newExpression copyPrecplusifFalse(precplusifFalse); 
    Xyce::Util::newExpression assignPrecplusifFalse; 
    assignPrecplusifFalse = precplusifFalse; 

    double result;
    precplusifFalse.evaluateFunction(result);       EXPECT_EQ( result, 2.0 );
    copyPrecplusifFalse.evaluateFunction(result);   EXPECT_EQ( result, 2.0 );
    assignPrecplusifFalse.evaluateFunction(result); EXPECT_EQ( result, 2.0 );
    OUTPUT_MACRO2(Double_Parser_ternary_precedence, precplusif, precplusifFalse) 
  }
}

TEST ( Double_Parser_ternary_precedence, precplusparenif)
{
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression is the RHS of a .func statement:  .func precplusparenif(X) {if(((1+x)>0),(2*x),(0+2))}
  //Xyce::Util::newExpression precplusparenifExpression(std::string("if(((1+x)>0),(2*x),(0+2))"), testGroup);
  Teuchos::RCP<Xyce::Util::newExpression> precplusparenifExpression
     = Teuchos::rcp(new Xyce::Util::newExpression(std::string("if(((1+x)>0),(2*x),(0+2))"), testGroup));

  std::vector<std::string> precplusparenifArgStrings; 

  Xyce::Util::newExpression precplusparenif_LHS (std::string("precplusparenif(X)"), testGroup);
  precplusparenif_LHS.lexAndParseExpression();
  precplusparenif_LHS.getFuncPrototypeArgStrings(precplusparenifArgStrings);
  precplusparenifExpression->setFunctionArgStringVec ( precplusparenifArgStrings );
  precplusparenifExpression->lexAndParseExpression();
  {
    Xyce::Util::newExpression precplusparenifTrue(std::string("precplusparenif(4)"), testGroup);
    precplusparenifTrue.lexAndParseExpression();
    std::string precplusparenifName;
    precplusparenif_LHS.getFuncPrototypeName (precplusparenifName);
    precplusparenifTrue.attachFunctionNode(precplusparenifName ,  precplusparenifExpression);

    Xyce::Util::newExpression copyPrecplusparenifTrue(precplusparenifTrue); 
    Xyce::Util::newExpression assignPrecplusparenifTrue; 
    assignPrecplusparenifTrue = precplusparenifTrue; 

    double result;
    precplusparenifTrue.evaluateFunction(result);       EXPECT_EQ( result, 8.0 );
    copyPrecplusparenifTrue.evaluateFunction(result);   EXPECT_EQ( result, 8.0 );
    assignPrecplusparenifTrue.evaluateFunction(result); EXPECT_EQ( result, 8.0 );
    OUTPUT_MACRO2(Double_Parser_ternary_precedence, precplusparenif, precplusparenifTrue) 
  }
  {
    Xyce::Util::newExpression precplusparenifFalse(std::string("precplusparenif(-4)"), testGroup);
    precplusparenifFalse.lexAndParseExpression();
    std::string precplusparenifName;
    precplusparenif_LHS.getFuncPrototypeName (precplusparenifName);
    precplusparenifFalse.attachFunctionNode(precplusparenifName ,  precplusparenifExpression);

    Xyce::Util::newExpression copyPrecplusparenifFalse(precplusparenifFalse); 
    Xyce::Util::newExpression assignPrecplusparenifFalse; 
    assignPrecplusparenifFalse = precplusparenifFalse; 

    double result;
    precplusparenifFalse.evaluateFunction(result);       EXPECT_EQ( result, 2.0 );
    copyPrecplusparenifFalse.evaluateFunction(result);   EXPECT_EQ( result, 2.0 );
    assignPrecplusparenifFalse.evaluateFunction(result); EXPECT_EQ( result, 2.0 );
    OUTPUT_MACRO2(Double_Parser_ternary_precedence, precplusparenif, precplusparenifFalse) 
  }
}

TEST ( Double_Parser_Func_Test, longArgList)
{
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression is the RHS of a .func statement:  .func crazy(X,A,B,Y,E,C,D,F) {E}
  Teuchos::RCP<Xyce::Util::newExpression> crazyExpression
     = Teuchos::rcp(new Xyce::Util::newExpression(std::string("E"), testGroup));

  std::vector<std::string> crazyArgStrings;
  Xyce::Util::newExpression crazy_LHS (std::string("crazy(X,A,B,Y,E,C,D,F)"), testGroup);
  crazy_LHS.lexAndParseExpression();
  crazy_LHS.getFuncPrototypeArgStrings(crazyArgStrings); // should be std::vector<std::string> = { "X","A","B","Y","E","C","D","F" }
  crazyExpression->setFunctionArgStringVec ( crazyArgStrings );
  crazyExpression->lexAndParseExpression();
  {
    Xyce::Util::newExpression crazyTrue(std::string("crazy(0,0,0,0,4,0,0,0)"), testGroup);
    crazyTrue.lexAndParseExpression();
    std::string crazyName;
    crazy_LHS.getFuncPrototypeName (crazyName); // should be "crazy"
    crazyTrue.attachFunctionNode(crazyName ,  crazyExpression);

    Xyce::Util::newExpression copyCrazyTrue(crazyTrue); 
    Xyce::Util::newExpression assignCrazyTrue; 
    assignCrazyTrue = crazyTrue; 

    double result;
    crazyTrue.evaluateFunction(result);       EXPECT_EQ( result, 4.0 );
    copyCrazyTrue.evaluateFunction(result);   EXPECT_EQ( result, 4.0 );
    assignCrazyTrue.evaluateFunction(result); EXPECT_EQ( result, 4.0 );
    OUTPUT_MACRO2(Double_Parser_Func_Test, longArgList, crazyTrue) 
  }
  {
    Xyce::Util::newExpression crazyFalse(std::string("crazy(0,0,0,0,-4,0,0,0)"), testGroup);
    crazyFalse.lexAndParseExpression();
    std::string crazyName;
    crazy_LHS.getFuncPrototypeName (crazyName); // should be "crazy"
    crazyFalse.attachFunctionNode(crazyName ,  crazyExpression);

    Xyce::Util::newExpression copyCrazyFalse(crazyFalse); 
    Xyce::Util::newExpression assignCrazyFalse; 
    assignCrazyFalse = crazyFalse; 

    double result;
    crazyFalse.evaluateFunction(result);       EXPECT_EQ( result, -4.0 );
    copyCrazyFalse.evaluateFunction(result);   EXPECT_EQ( result, -4.0 );
    assignCrazyFalse.evaluateFunction(result); EXPECT_EQ( result, -4.0 );
    OUTPUT_MACRO2(Double_Parser_Func_Test, longArgList, crazyFalse) 
  }
}

//-------------------------------------------------------------------------------
// tests are taken from the "ifstatement.cir" Xyce regression test
class ifStatementExpressionGroup : public Xyce::Util::baseExpressionGroup
{
  public:
    ifStatementExpressionGroup () : Xyce::Util::baseExpressionGroup(), time(0.0)  {};
    ~ifStatementExpressionGroup () {};

    virtual double getTime() { return time; };
    void setTime(double t) { time = t; };

    virtual bool getSolutionVal(const std::string & nodeName, double & retval )
    {
      std::string tmp = nodeName; Xyce::Util::toLower(tmp);
      if (tmp==std::string("b2")) { retval = B2; return true; }
      else if (tmp==std::string("6")) { retval = v6; return true; }
      else if (tmp==std::string("7")) { retval = v7; return true; }
      else { return 0.0; return false; }
    }

    void setSoln(const std::string & nodeName, double val)
    {
      std::string tmp = nodeName; Xyce::Util::toLower(tmp);
      if (tmp==std::string("b2")) { B2 = val; }
      else if (tmp==std::string("6")) { v6 = val; }
      else if (tmp==std::string("7")) { v7 = val; }
    }

#if 0
    void addFunction (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
    {
      std::string lowerName = name;
      Xyce::Util::toLower(lowerName);

      functions_[lowerName] = exp;
    };

    bool getFunction (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
    {
      bool retval=true;

      std::string lowerName = name;
      Xyce::Util::toLower(lowerName);

      if (functions_.find(lowerName) != functions_.end()) { exp = functions_[lowerName]; }
      else { retval = false; }

      return retval;
    }
#endif

  private:
    //std::unordered_map <std::string, Teuchos::RCP<Xyce::Util::newExpression> >  functions_;

    double time;
    double B2;
    double v6;
    double v7;
};

// currently doesn't work
TEST ( Double_Parser_ifstatement, ifmin_ifmax_func)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  // The ifmin expression is the RHS of a .func statement:  .func ifmin (a,b) {if(a<b, a, b)}
  Teuchos::RCP<Xyce::Util::newExpression> ifmin
     = Teuchos::rcp(new Xyce::Util::newExpression(std::string("if(a<b, a, b)"), baseGroup));
  std::vector<std::string> ifminArgStrings; // = { std::string("a"), std::string("b") };
  Xyce::Util::newExpression ifmin_LHS (std::string("ifmin (a,b)"), baseGroup);
  ifmin_LHS.lexAndParseExpression();
  ifmin_LHS.getFuncPrototypeArgStrings(ifminArgStrings);
  ifmin->setFunctionArgStringVec ( ifminArgStrings );
  ifmin->lexAndParseExpression();

  // The ifmax expression is the RHS of a .func statement:  .func ifmax (a,b) {if(a>b, a, b)}
  Teuchos::RCP<Xyce::Util::newExpression> ifmax
     = Teuchos::rcp(new Xyce::Util::newExpression(std::string("if(a>b, a, b)"), baseGroup));
  std::vector<std::string> ifmaxArgStrings; // = { std::string("a"), std::string("b") };
  Xyce::Util::newExpression ifmax_LHS (std::string("ifmax (a,b)"), baseGroup);
  ifmax_LHS.lexAndParseExpression();
  ifmax_LHS.getFuncPrototypeArgStrings(ifmaxArgStrings);
  ifmax->setFunctionArgStringVec ( ifmaxArgStrings );
  ifmax->lexAndParseExpression();

  std::string ifmaxName;//="ifmax";
  std::string ifminName;//="ifmin";

  ifmax_LHS.getFuncPrototypeName (ifmaxName);
  ifmin_LHS.getFuncPrototypeName (ifminName);

  //ifGroup->addFunction( ifmaxName ,  ifmax);
  //ifGroup->addFunction( ifminName ,  ifmin);

  // these expressions uses the .func ifmin and ifmax
  {
    Xyce::Util::newExpression e3(std::string("ifmax(ifmin(-I(B2), 2.5), 1.5)"), baseGroup);
    e3.lexAndParseExpression();
    e3.attachFunctionNode(ifmaxName, ifmax);
    e3.attachFunctionNode(ifminName, ifmin);

    Xyce::Util::newExpression copy_e3(e3); 
    Xyce::Util::newExpression assign_e3; 
    assign_e3 = e3; 

    int numpoints=100;
    double tmax=1.0;
    double time=0.0;
    double dt=tmax/(numpoints-1);
    std::vector<double> refRes(numpoints), result(numpoints);
    std::vector<double> copyResult(numpoints), assignResult(numpoints);
    for (int ii=0;ii<numpoints;ii++,time+=dt)
    {
      double V0=2, VA=1, FREQ=1, mpi = M_PI;
      double v1 = 0.0;
      if (time <= 0) { v1 = (V0) + (VA) * std::sin (0.0); }
      else { v1 = (V0) + (VA) * std::sin (2.0*mpi*((FREQ)*time + (0.0)/360)) * std::exp( -(time*(0.0))); }
      double v2 = 2*v1;
      double b2 = -0.5*v2;

      ifGroup->setTime(time);
      ifGroup->setSoln(std::string("b2"),b2);
      e3.evaluateFunction(result[ii]);
      copy_e3.evaluateFunction(copyResult[ii]);
      assign_e3.evaluateFunction(assignResult[ii]);
      refRes[ii] = std::max(std::min(-b2,2.5),1.5);
    }
    EXPECT_EQ( result, refRes);
    EXPECT_EQ( copyResult, refRes);
    EXPECT_EQ( assignResult, refRes);
    OUTPUT_MACRO2(Double_Parser_ifstatement, ifmin_ifmax_func, e3) 
  }
}

TEST ( Double_Parser_ifstatement, simple_nested_func)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  // The doubleIt expression is the RHS of a .func statement:  .func doubleIt(a) {2*a)}
  Teuchos::RCP<Xyce::Util::newExpression> doubleIt
     = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2*a"), baseGroup));
  std::vector<std::string> doubleItArgStrings;
  Xyce::Util::newExpression doubleIt_LHS (std::string("doubleIt(a)"), baseGroup);
  doubleIt_LHS.lexAndParseExpression();
  doubleIt_LHS.getFuncPrototypeArgStrings(doubleItArgStrings);
  doubleIt->setFunctionArgStringVec ( doubleItArgStrings );
  doubleIt->lexAndParseExpression();

  // The tripleIt expression is the RHS of a .func statement:  .func tripleIt (a) {3*a)}
  Teuchos::RCP<Xyce::Util::newExpression> tripleIt
     = Teuchos::rcp(new Xyce::Util::newExpression(std::string("3*a"), baseGroup));
  std::vector<std::string> tripleItArgStrings;
  Xyce::Util::newExpression tripleIt_LHS (std::string("tripleIt(a)"), baseGroup);
  tripleIt_LHS.lexAndParseExpression();
  tripleIt_LHS.getFuncPrototypeArgStrings(tripleItArgStrings);
  tripleIt->setFunctionArgStringVec ( tripleItArgStrings );
  tripleIt->lexAndParseExpression();

  std::string doubleItName;//="doubleIt";
  std::string tripleItName;//="tripleIt";
  doubleIt_LHS.getFuncPrototypeName(doubleItName);
  tripleIt_LHS.getFuncPrototypeName(tripleItName);

  // these expressions uses the .func doubleIt and tripleIt
  {
    Xyce::Util::newExpression e3(std::string("tripleIt(doubleIt(-I(B2)))"), baseGroup);
    e3.lexAndParseExpression();
    e3.attachFunctionNode(tripleItName, tripleIt);
    e3.attachFunctionNode(doubleItName, doubleIt);

    Xyce::Util::newExpression copy_e3(e3); 
    Xyce::Util::newExpression assign_e3; 
    assign_e3 = e3; 

    double result, refRes, copyResult, assignResult;
    double b2 = -0.5;
    ifGroup->setSoln(std::string("b2"),b2);
    e3.evaluateFunction(result);
    copy_e3.evaluateFunction(copyResult);
    assign_e3.evaluateFunction(assignResult);
    refRes = -b2*2*3;

    EXPECT_EQ( result, refRes);
    EXPECT_EQ( copyResult, refRes);
    EXPECT_EQ( assignResult, refRes);
    OUTPUT_MACRO2(Double_Parser_ifstatement, simple_nested_func, e3) 
  }
}

TEST ( Double_Parser_ifstatement, min_max)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  // these expressions uses min and max AST nodes
  Xyce::Util::newExpression e4(std::string("max(min(-I(B2), 2.5), 1.5)"), baseGroup);
  e4.lexAndParseExpression();

  Xyce::Util::newExpression copy_e4(e4); 
  Xyce::Util::newExpression assign_e4; 
  assign_e4 = e4; 

  int numpoints=100;
  double tmax=1.0;
  double time=0.0;
  double dt=tmax/(numpoints-1);
  std::vector<double> refRes(numpoints), result(numpoints);
  std::vector<double> copyResult(numpoints), assignResult(numpoints);
  for (int ii=0;ii<numpoints;ii++,time+=dt)
  {
    double V0=2, VA=1, FREQ=1, mpi = M_PI;
    double v1 = 0.0;
    if (time <= 0) { v1 = (V0) + (VA) * std::sin (0.0); }
    else { v1 = (V0) + (VA) * std::sin (2.0*mpi*((FREQ)*time + (0.0)/360)) * std::exp( -(time*(0.0))); }
    double v2 = 2*v1;
    double b2 = -0.5*v2;

    ifGroup->setTime(time);
    ifGroup->setSoln(std::string("b2"),b2);
    e4.evaluateFunction(result[ii]);
    copy_e4.evaluateFunction(copyResult[ii]);
    assign_e4.evaluateFunction(assignResult[ii]);
    refRes[ii] = std::max(std::min(-b2,2.5),1.5);
  }
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( copyResult, refRes);
  EXPECT_EQ( assignResult, refRes);
  OUTPUT_MACRO2(Double_Parser_ifstatement, min_max,e4) 
}

TEST ( Double_Parser_ifstatement, limit)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  // these expressions uses limit AST node
  Xyce::Util::newExpression e5(std::string("limit(-I(B2),1.5,2.5)"), baseGroup);
  e5.lexAndParseExpression();

  Xyce::Util::newExpression copy_e5(e5); 
  Xyce::Util::newExpression assign_e5; 
  assign_e5 = e5; 

  int numpoints=100;
  double tmax=1.0;
  double time=0.0;
  double dt=tmax/(numpoints-1);
  std::vector<double> refRes(numpoints), result(numpoints);
  std::vector<double> copyResult(numpoints), assignResult(numpoints);
  for (int ii=0;ii<numpoints;ii++,time+=dt)
  {
    double V0=2, VA=1, FREQ=1, mpi = M_PI;
    double v1 = 0.0;
    if (time <= 0) { v1 = (V0) + (VA) * std::sin (0.0); }
    else { v1 = (V0) + (VA) * std::sin (2.0*mpi*((FREQ)*time + (0.0)/360)) * std::exp( -(time*(0.0))); }
    double v2 = 2*v1;
    double b2 = -0.5*v2;

    ifGroup->setTime(time);
    ifGroup->setSoln(std::string("b2"),b2);
    e5.evaluateFunction(result[ii]);
    copy_e5.evaluateFunction(copyResult[ii]);
    assign_e5.evaluateFunction(assignResult[ii]);
    refRes[ii] = std::max(std::min(-b2,2.5),1.5);
  }
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( copyResult, refRes);
  EXPECT_EQ( assignResult, refRes);
  OUTPUT_MACRO2(Double_Parser_ifstatement,limit, e5) 
}

TEST ( Double_Parser_ifstatement, or_true)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e8(std::string("IF(((V(6) > 1.5) | (V(7) < 1.5)), 3, 1)"), baseGroup);
  e8.lexAndParseExpression();

  Xyce::Util::newExpression copy_e8(e8); 
  Xyce::Util::newExpression assign_e8; 
  assign_e8 = e8; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),1.0);

  double result=0.0;
  e8.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e8.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e8.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, or_true,e8) 
}

TEST ( Double_Parser_ifstatement, or_false)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e9(std::string("IF(((V(6) > 1.5) | (V(7) > 1.5)), 3, 1)"), baseGroup);
  e9.lexAndParseExpression();

  Xyce::Util::newExpression copy_e9(e9); 
  Xyce::Util::newExpression assign_e9; 
  assign_e9 = e9; 

  ifGroup->setSoln(std::string("6"),1.0);
  ifGroup->setSoln(std::string("7"),1.0);

  double result=0.0;
  e9.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e9.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e9.evaluateFunction(result); EXPECT_EQ( result, 1.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, or_false,e9) 
}

TEST ( Double_Parser_ifstatement, hspice_or_true)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e8(std::string("IF(((V(6) > 1.5) || (V(7) < 1.5)), 3, 1)"), baseGroup);
  e8.lexAndParseExpression();

  Xyce::Util::newExpression copy_e8(e8); 
  Xyce::Util::newExpression assign_e8; 
  assign_e8 = e8; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),1.0);

  double result=0.0;
  e8.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e8.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e8.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, hspice_or_true,e8) 
}

TEST ( Double_Parser_ifstatement, hspice_or_false)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e9(std::string("IF(((V(6) > 1.5) || (V(7) > 1.5)), 3, 1)"), baseGroup);
  e9.lexAndParseExpression();

  Xyce::Util::newExpression copy_e9(e9); 
  Xyce::Util::newExpression assign_e9; 
  assign_e9 = e9; 

  ifGroup->setSoln(std::string("6"),1.0);
  ifGroup->setSoln(std::string("7"),1.0);

  double result=0.0;
  e9.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e9.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e9.evaluateFunction(result); EXPECT_EQ( result, 1.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, hspice_or_false,e9) 
}

TEST ( Double_Parser_ifstatement, and_true)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e8(std::string("IF(((V(6) > 1.5) & (V(7) < 1.5)), 3, 1)"), baseGroup);
  e8.lexAndParseExpression();

  Xyce::Util::newExpression copy_e8(e8); 
  Xyce::Util::newExpression assign_e8; 
  assign_e8 = e8; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),1.0);

  double result=0.0;
  e8.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e8.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e8.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, and_true,e8) 
}

TEST ( Double_Parser_ifstatement, and_false)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e9(std::string("IF(((V(6) > 1.5) & (V(7) > 1.5)), 3, 1)"), baseGroup);
  e9.lexAndParseExpression();

  Xyce::Util::newExpression copy_e9(e9); 
  Xyce::Util::newExpression assign_e9; 
  assign_e9 = e9; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),1.0);

  double result=0.0;
  e9.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e9.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e9.evaluateFunction(result); EXPECT_EQ( result, 1.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, and_false,e9) 
}

TEST ( Double_Parser_ifstatement, hspice_and_true)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e8(std::string("IF(((V(6) > 1.5) && (V(7) < 1.5)), 3, 1)"), baseGroup);
  e8.lexAndParseExpression();

  Xyce::Util::newExpression copy_e8(e8); 
  Xyce::Util::newExpression assign_e8; 
  assign_e8 = e8; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),1.0);

  double result=0.0;
  e8.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e8.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e8.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, hspice_and_true,e8) 
}

TEST ( Double_Parser_ifstatement, hspice_and_false)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e9(std::string("IF(((V(6) > 1.5) && (V(7) > 1.5)), 3, 1)"), baseGroup);
  e9.lexAndParseExpression();

  Xyce::Util::newExpression copy_e9(e9); 
  Xyce::Util::newExpression assign_e9; 
  assign_e9 = e9; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),1.0);

  double result=0.0;
  e9.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e9.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e9.evaluateFunction(result); EXPECT_EQ( result, 1.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, hspice_and_false,e9) 
}

TEST ( Double_Parser_ifstatement, xor_true)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e8(std::string("IF(((V(6) > 1.5) ^ (V(7) < 1.5)), 3, 1)"), baseGroup);
  e8.lexAndParseExpression();

  Xyce::Util::newExpression copy_e8(e8); 
  Xyce::Util::newExpression assign_e8; 
  assign_e8 = e8; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),1.0);

  double result=0.0;
  e8.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e8.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e8.evaluateFunction(result); EXPECT_EQ( result, 1.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, xor_true,e8) 
}

TEST ( Double_Parser_ifstatement, xor_false)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e9(std::string("IF(((V(6) > 1.5) ^ (V(7) > 1.5)), 3, 1)"), baseGroup);
  e9.lexAndParseExpression();

  Xyce::Util::newExpression copy_e9(e9); 
  Xyce::Util::newExpression assign_e9; 
  assign_e9 = e9; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),1.0);

  double result=0.0;
  e9.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e9.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e9.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, xor_false,e9) 
}

TEST ( Double_Parser_ifstatement, neq)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e10(std::string("IF((V(6) != V(7)), 3, 1)"), baseGroup);
  e10.lexAndParseExpression();

  Xyce::Util::newExpression copy_e10(e10); 
  Xyce::Util::newExpression assign_e10; 
  assign_e10 = e10; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),1.0);

  double result=0.0;
  e10.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e10.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e10.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, neq, e10) 
}

TEST ( Double_Parser_ifstatement, not)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e11(std::string("IF( ~(V(6) > V(7)), 3, 1)"), baseGroup);
  e11.lexAndParseExpression();

  Xyce::Util::newExpression copy_e11(e11); 
  Xyce::Util::newExpression assign_e11; 
  assign_e11 = e11; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),1.0);

  double result=0.0;
  e11.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e11.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e11.evaluateFunction(result); EXPECT_EQ( result, 1.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, not, e11) 
}

TEST ( Double_Parser_ifstatement, equiv)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e11(std::string("IF(V(6) == V(7), 3, 1)"), baseGroup);
  e11.lexAndParseExpression();

  Xyce::Util::newExpression copy_e11(e11); 
  Xyce::Util::newExpression assign_e11; 
  assign_e11 = e11; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),1.0);

  double result=0.0;
  e11.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e11.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e11.evaluateFunction(result); EXPECT_EQ( result, 1.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, equiv, e11) 
}

TEST ( Double_Parser_ifstatement, ge1)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e11(std::string("IF(V(6) >= V(7), 3, 1)"), baseGroup);
  e11.lexAndParseExpression();

  Xyce::Util::newExpression copy_e11(e11); 
  Xyce::Util::newExpression assign_e11; 
  assign_e11 = e11; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),1.0);

  double result=0.0;
  e11.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e11.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e11.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, ge1, e11) 
}

TEST ( Double_Parser_ifstatement, ge2)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e11(std::string("IF(V(6) >= V(7), 3, 1)"), baseGroup);
  e11.lexAndParseExpression();

  Xyce::Util::newExpression copy_e11(e11); 
  Xyce::Util::newExpression assign_e11; 
  assign_e11 = e11; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),2.0);

  double result=0.0;
  e11.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e11.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e11.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, ge2, e11) 
}

TEST ( Double_Parser_ifstatement, ge3)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e11(std::string("IF(V(6) >= V(7), 3, 1)"), baseGroup);
  e11.lexAndParseExpression();

  Xyce::Util::newExpression copy_e11(e11); 
  Xyce::Util::newExpression assign_e11; 
  assign_e11 = e11; 

  ifGroup->setSoln(std::string("6"),1.0);
  ifGroup->setSoln(std::string("7"),2.0);

  double result=0.0;
  e11.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e11.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e11.evaluateFunction(result); EXPECT_EQ( result, 1.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, ge3, e11) 
}


TEST ( Double_Parser_ifstatement, le1)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e11(std::string("IF(V(6) <= V(7), 3, 1)"), baseGroup);
  e11.lexAndParseExpression();

  Xyce::Util::newExpression copy_e11(e11); 
  Xyce::Util::newExpression assign_e11; 
  assign_e11 = e11; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),1.0);

  double result=0.0;
  e11.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e11.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e11.evaluateFunction(result); EXPECT_EQ( result, 1.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, le1, e11) 
}

TEST ( Double_Parser_ifstatement, le2)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e11(std::string("IF(V(6) <= V(7), 3, 1)"), baseGroup);
  e11.lexAndParseExpression();

  Xyce::Util::newExpression copy_e11(e11); 
  Xyce::Util::newExpression assign_e11; 
  assign_e11 = e11; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),2.0);

  double result=0.0;
  e11.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e11.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e11.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, le2, e11) 
}

TEST ( Double_Parser_ifstatement, le3)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Xyce::Util::newExpression e11(std::string("IF(V(6) <= V(7), 3, 1)"), baseGroup);
  e11.lexAndParseExpression();

  Xyce::Util::newExpression copy_e11(e11); 
  Xyce::Util::newExpression assign_e11; 
  assign_e11 = e11; 

  ifGroup->setSoln(std::string("6"),1.0);
  ifGroup->setSoln(std::string("7"),2.0);

  double result=0.0;
  e11.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e11.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e11.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(Double_Parser_ifstatement, le3, e11) 
}

// from "ifstatement.cir":
// Also test modulus, along with its precedence vs. + - * and /.
// Because of precedence, this should evaluate to 4.
TEST ( Double_Parser_modulus, test1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new testExpressionGroup() );

  // 2 + 6*5/2%4 - 1  = 2 + 30/2%4 -1 = 2 + 15%4 - 1 = 2 + 3 - 1 = 4
  Xyce::Util::newExpression p1(std::string("2 + 6*5/2%4 - 1"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 

  double result=0.0;
  p1.evaluateFunction(result);        EXPECT_EQ( result, 4.0);
  copy_p1.evaluateFunction(result);   EXPECT_EQ( result, 4.0);
  assign_p1.evaluateFunction(result); EXPECT_EQ( result, 4.0);
  OUTPUT_MACRO2(Double_Parser_modulus, test1, p1) 
}

TEST ( Double_Parser_modulus, test2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new testExpressionGroup() );

  Xyce::Util::newExpression p1(std::string("15%4"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 

  double result=0.0;
  p1.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_p1.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_p1.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(Double_Parser_modulus, test2, p1) 
}

//-------------------------------------------------------------------------------
// table tests
//
// adapted from break.cir
TEST ( Double_Parser_table_Test, break1)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = timeDepGroup;
  Xyce::Util::newExpression tableExpression(std::string("Table(time, 0, 0, 0.3, 0, 0.301, 2, 0.302, 2, 0.6, 1, 1, 1)"), grp);
  tableExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_tableExpression(tableExpression); 
  Xyce::Util::newExpression assign_tableExpression; 
  assign_tableExpression = tableExpression; 

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<double> refRes = { 0, 0, 2, 2, 1, 1 };
  std::vector<double> result(times.size(),0.0);
  std::vector<double> copyResult(times.size(),0.0);
  std::vector<double> assignResult(times.size(),0.0);

  for (int ii=0;ii<times.size();ii++) 
  { 
    timeDepGroup->setTime(times[ii]); 
    tableExpression.evaluateFunction(result[ii]); 
    copy_tableExpression.evaluateFunction(copyResult[ii]); 
    assign_tableExpression.evaluateFunction(assignResult[ii]); 
  }
  EXPECT_EQ(refRes,result);
  EXPECT_EQ(refRes,copyResult);
  EXPECT_EQ(refRes,assignResult);
  OUTPUT_MACRO2(Double_Parser_table_Test, break1, tableExpression) 
}

TEST ( Double_Parser_table_Test, break2)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = timeDepGroup;
  Xyce::Util::newExpression tableExpression(std::string("Table({time} 0, 0, 0.3, 0, 0.301, 2, 0.302, 2, 0.6, 1, 1, 1)"), grp);
  tableExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_tableExpression(tableExpression); 
  Xyce::Util::newExpression assign_tableExpression; 
  assign_tableExpression = tableExpression; 

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<double> refRes = { 0, 0, 2, 2, 1, 1 };
  std::vector<double> result(times.size(),0.0);
  std::vector<double> copyResult(times.size(),0.0);
  std::vector<double> assignResult(times.size(),0.0);

  for (int ii=0;ii<times.size();ii++) 
  { 
    timeDepGroup->setTime(times[ii]); 
    tableExpression.evaluateFunction(result[ii]); 
    copy_tableExpression.evaluateFunction(copyResult[ii]); 
    assign_tableExpression.evaluateFunction(assignResult[ii]); 
  }
  EXPECT_EQ(refRes,result);
  EXPECT_EQ(refRes,copyResult);
  EXPECT_EQ(refRes,assignResult);
  OUTPUT_MACRO2(Double_Parser_table_Test, break2, tableExpression) 
}

class tempDepExpressionGroup : public Xyce::Util::baseExpressionGroup
{
  public:
    tempDepExpressionGroup () : Xyce::Util::baseExpressionGroup(), temp(0.0)  {};
    ~tempDepExpressionGroup () {};
    virtual double getTemp() { return temp; };
    void setTemp(double t) { temp = t; };
    double temp;
};

// adapted from power_thermalres_gear.cir
TEST ( Double_Parser_table_Test, power_thermalres)
{
  Teuchos::RCP<tempDepExpressionGroup> tempDepGroup = Teuchos::rcp(new tempDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  grp = tempDepGroup;
  {
    Xyce::Util::newExpression resistivity(std::string("table(temp+273.15, 0, 0.5e-9, 100, 3e-9, 1000, 6.6e-8)"), grp);
    resistivity.lexAndParseExpression();

    Xyce::Util::newExpression copy_resistivity(resistivity); 
    Xyce::Util::newExpression assign_resistivity; 
    assign_resistivity = resistivity; 

    std::vector<double> temps = { 0-273.15, 100-273.15, 1000-273.15 };
    std::vector<double> refRes = { 0.5e-9, 3e-9, 6.6e-8 };
    std::vector<double> result(temps.size(),0.0);
    std::vector<double> copyResult(temps.size(),0.0);
    std::vector<double> assignResult(temps.size(),0.0);

    for (int ii=0;ii<temps.size();ii++) 
    { 
      tempDepGroup->setTemp(temps[ii]); 
      resistivity.evaluateFunction(result[ii]); 
      copy_resistivity.evaluateFunction(copyResult[ii]); 
      assign_resistivity.evaluateFunction(assignResult[ii]); 
    }
    EXPECT_EQ(refRes,result);
    EXPECT_EQ(refRes,copyResult);
    EXPECT_EQ(refRes,assignResult);
    OUTPUT_MACRO2(Double_Parser_table_Test, power_thermalres, resistivity) 
  }
  {
    Xyce::Util::newExpression heatcapacity(std::string("8.92e+3*table(temp+273.15, 0, 1, 1000, 1500)"), grp);
    heatcapacity.lexAndParseExpression();
    
    Xyce::Util::newExpression copy_heatcapacity(heatcapacity); 
    Xyce::Util::newExpression assign_heatcapacity; 
    assign_heatcapacity = heatcapacity; 

    std::vector<double> temps = { 0-273.15, 1000-273.15 };
    std::vector<double> refRes = { 8.92e+3, 8.92e+3*1500 };
    std::vector<double> result(temps.size(),0.0);
    std::vector<double> copyResult(temps.size(),0.0);
    std::vector<double> assignResult(temps.size(),0.0);

    for (int ii=0;ii<temps.size();ii++) 
    { 
      tempDepGroup->setTemp(temps[ii]); 
      heatcapacity.evaluateFunction(result[ii]); 
      copy_heatcapacity.evaluateFunction(copyResult[ii]); 
      assign_heatcapacity.evaluateFunction(assignResult[ii]); 
    }
    EXPECT_EQ(refRes,result);
    EXPECT_EQ(refRes,copyResult);
    EXPECT_EQ(refRes,assignResult);
    OUTPUT_MACRO2(Double_Parser_table_Test, power_thermalres, heatcapacity) 
  }
}

class Bsrc_C1_ExpressionGroup : public Xyce::Util::baseExpressionGroup
{
  public:
    Bsrc_C1_ExpressionGroup () :
      Xyce::Util::baseExpressionGroup(), time(0.0), ONEval(0.0), TWOval(0.0) {};
    ~Bsrc_C1_ExpressionGroup () {};

  virtual bool getSolutionVal(const std::string & nodeName, double & retval )
  {
    std::string tmp = nodeName; Xyce::Util::toLower(tmp);
    if (tmp==std::string("1")) { retval = ONEval; return true; }
    else if (tmp==std::string("2")) { retval = TWOval; return true; }
    else { return 0.0; return false; }
  }

  void setSoln(const std::string & nodeName, double val)
  {
    std::string tmp = nodeName; Xyce::Util::toLower(tmp);
    if (tmp==std::string("1")) { ONEval = val; }
    else if (tmp==std::string("2")) { TWOval = val; }
  }

  virtual double getTime() { return time; };
  void setTime(double t) { time = t; };

#if 0
  void addParam (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
  {
    std::string lowerName = name;
    Xyce::Util::toLower(lowerName);

    parameters_[lowerName] = exp;
  };

  bool getParam       (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
  {
    bool retval=true;
    std::string lowerName = name;
    Xyce::Util::toLower(lowerName);
    if (parameters_.find(lowerName) != parameters_.end()) { exp = parameters_[lowerName]; }
    else { retval = false; }
    return retval;
  }
#endif

  private:
    //std::unordered_map <std::string, Teuchos::RCP<Xyce::Util::newExpression> >  parameters_;
    double time, ONEval, TWOval;
};



class Bsrc_C1_ExpressionGroup_noparam : public Xyce::Util::baseExpressionGroup
{
  public:
    Bsrc_C1_ExpressionGroup_noparam () :
      Xyce::Util::baseExpressionGroup(), time(0.0), ONEval(0.0), TWOval(0.0) {};
    ~Bsrc_C1_ExpressionGroup_noparam () {};

  virtual bool getSolutionVal(const std::string & nodeName, double & retval )
  {
    std::string tmp = nodeName; Xyce::Util::toLower(tmp);
    if (tmp==std::string("1")) { retval = ONEval; return true; }
    else if (tmp==std::string("2")) { retval = TWOval; return true; }
    else { return 0.0; return false; }
  }

  void setSoln(const std::string & nodeName, double val)
  {
    std::string tmp = nodeName; Xyce::Util::toLower(tmp);
    if (tmp==std::string("1")) { ONEval = val; }
    else if (tmp==std::string("2")) { TWOval = val; }
  }

  virtual double getTime() { return time; };
  void setTime(double t) { time = t; };

  private:
    double time, ONEval, TWOval;
};

#if 0
// adapted from Bsrc_C1.cir.
// See the Double_Parser_Param_Test.test2  test below,
// which also tests the first-argument expression and uses
// some of the same machinery.
//
// the mechanics of the table seem to work here, but I haven't generated a good gold standard yet
// so, disabling for now
TEST ( Double_Parser_table_Test, Bsrc_C1_withoutParens)
{
  Bsrc_C1_ExpressionGroup grp;
  {
    // this is a nice test b/c it has curly braces around the first expression, which is one of the supported formats
    Xyce::Util::newExpression BE_Dig(std::string("TABLE({ V(2) * (V(1) + 30) / 60 } 0.0000000, 0,0.0312500, 0,0.0312813, 1,0.0625000, 1,0.0625313, 2,0.0937500, 2,0.0937813, 3,0.1250000, 3,0.1250313, 4,0.1562500, 4,0.1562813, 5,0.1875000, 5,0.1875313, 6,0.2187500, 6,0.2187813, 7,0.2500000, 7,0.2500313, 8,0.2812500, 8,0.2812813, 9,0.3125000, 9,0.3125313, 10,0.3437500, 10,0.3437813, 11,0.3750000, 11,0.3750313, 12,0.4062500, 12,0.4062813, 13,0.4375000, 13,0.4375313, 14,0.4687500, 14,0.4687813, 15,0.5000000, 15,0.5000313, 16,0.5312500, 16,0.5312813, 17,0.5625000, 17,0.5625313, 18,0.5937500, 18,0.5937813, 19,0.6250000, 19,0.6250313, 20,0.6562500, 20,0.6562813, 21,0.6875000, 21,0.6875313, 22,0.7187500, 22,0.7187813, 23,0.7500000, 23,0.7500313, 24,0.7812500, 24,0.7812813, 25,0.8125000, 25,0.8125313, 26,0.8437500, 26,0.8437813, 27,0.8750000, 27,0.8750313, 28,0.9062500, 28,0.9062813, 29,0.9375000, 29,0.9375313, 30,0.9687500, 30,0.9687813, 31,1.0000000, 31)"), grp);
    BE_Dig.lexAndParseExpression();

    Xyce::Util::newExpression v1exp(std::string("spice_sin(0, 20, 1k, -.25e-3, 0, 0)" ), grp); v1exp.lexAndParseExpression();
    Xyce::Util::newExpression v2exp(std::string("spice_pulse(0, 1, 0, 0.5us, 0.5us, 2us, 20us) " ), grp); v2exp.lexAndParseExpression();

    // in the original test,
    // V(1) is a sinewave that goes between +20 and -20
    // V(2) is a pulsed source that goes between 0 and 1.  PW is short.
    // The expression, V(2) * (V(1) + 30) / 60  has roughly the scaled shape of V(1) but is spiked/digitized.
    // When V(2) is zero, so is the expression.  When V(2) is 1, then expression is (V(1)+30)/60  max = 50/60, min=10/60
    //
    // The table itself is a digitized signal; like stairsteps, with 64 points.  It goes from 0 to 31 in increments of 1
    // The result is that the interpolated final output is really similar to the expression.
    int numpoints=2000;
    double tfinal = 0.0005;
    double dt = tfinal/(numpoints-1), time=0.0;
    std::vector<double> refRes(numpoints), result(numpoints);
    for (int ii=0;ii<numpoints;ii++,time+=dt)
    {
      grp.setTime(time);
      double v1Value(0.0),v2Value(0.0);
      v1exp.evaluateFunction(v1Value);
      v2exp.evaluateFunction(v2Value);
      grp.setSoln(std::string("1"),v1Value);
      grp.setSoln(std::string("2"),v2Value);
      BE_Dig.evaluateFunction(result[ii]);
      //std::cout.setf(std::ios::scientific);
      double firstExpVal = (v2Value * (v1Value + 30.0) / 60.0);
      //std::cout << time << "\t" << v1Value << "\t" << v2Value << "\t" << firstExpVal << "\t" << result[ii] <<std::endl;
    }

    EXPECT_EQ(refRes,result);
  }
}
#endif

TEST ( Double_Parser_table_Test, Bsrc_C1_pureArray)
{
  Teuchos::RCP<Bsrc_C1_ExpressionGroup> bsrc_C1_grp = Teuchos::rcp(new Bsrc_C1_ExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  grp = bsrc_C1_grp;
  {
    // this is a nice test b/c it has curly braces around the first expression, which is one of the supported formats
    Xyce::Util::newExpression BE_Dig(std::string("TABLE({ V(2) * (V(1) + 30) / 60 } 0.0000000, 0,0.0312500, 0,0.0312813, 1,0.0625000, 1,0.0625313, 2,0.0937500, 2,0.0937813, 3,0.1250000, 3,0.1250313, 4,0.1562500, 4,0.1562813, 5,0.1875000, 5,0.1875313, 6,0.2187500, 6,0.2187813, 7,0.2500000, 7,0.2500313, 8,0.2812500, 8,0.2812813, 9,0.3125000, 9,0.3125313, 10,0.3437500, 10,0.3437813, 11,0.3750000, 11,0.3750313, 12,0.4062500, 12,0.4062813, 13,0.4375000, 13,0.4375313, 14,0.4687500, 14,0.4687813, 15,0.5000000, 15,0.5000313, 16,0.5312500, 16,0.5312813, 17,0.5625000, 17,0.5625313, 18,0.5937500, 18,0.5937813, 19,0.6250000, 19,0.6250313, 20,0.6562500, 20,0.6562813, 21,0.6875000, 21,0.6875313, 22,0.7187500, 22,0.7187813, 23,0.7500000, 23,0.7500313, 24,0.7812500, 24,0.7812813, 25,0.8125000, 25,0.8125313, 26,0.8437500, 26,0.8437813, 27,0.8750000, 27,0.8750313, 28,0.9062500, 28,0.9062813, 29,0.9375000, 29,0.9375313, 30,0.9687500, 30,0.9687813, 31,1.0000000, 31)"), grp);
    BE_Dig.lexAndParseExpression();

    Xyce::Util::newExpression copy_BE_Dig(BE_Dig); 
    Xyce::Util::newExpression assign_BE_Dig; 
    assign_BE_Dig = BE_Dig; 

    Xyce::Util::newExpression BE_Dig_leftArg(std::string("V(2) * (V(1) + 30) / 60"),grp);
    BE_Dig_leftArg.lexAndParseExpression();

    std::vector<double> xa = { 0.0000000, 0.0312500, 0.0312813, 0.0625000, 0.0625313, 0.0937500, 0.0937813, 0.1250000, 0.1250313, 0.1562500, 0.1562813, 0.1875000, 0.1875313, 0.2187500, 0.2187813, 0.2500000, 0.2500313, 0.2812500, 0.2812813, 0.3125000, 0.3125313, 0.3437500, 0.3437813, 0.3750000, 0.3750313, 0.4062500, 0.4062813, 0.4375000, 0.4375313, 0.4687500, 0.4687813, 0.5000000, 0.5000313, 0.5312500, 0.5312813, 0.5625000, 0.5625313, 0.5937500, 0.5937813, 0.6250000, 0.6250313, 0.6562500, 0.6562813, 0.6875000, 0.6875313, 0.7187500, 0.7187813, 0.7500000, 0.7500313, 0.7812500, 0.7812813, 0.8125000, 0.8125313, 0.8437500, 0.8437813, 0.8750000, 0.8750313, 0.9062500, 0.9062813, 0.9375000, 0.9375313, 0.9687500, 0.9687813, 1.0000000 };

    std::vector<double> ya = { 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31 };

    Xyce::Util::newExpression BE_Dig_pureArray(BE_Dig_leftArg.getAst(),xa,ya,grp);
    BE_Dig_pureArray.lexAndParseExpression();

    Xyce::Util::newExpression v1exp(std::string("spice_sin(0, 20, 1k, -.25e-3, 0, 0)" ), grp); v1exp.lexAndParseExpression();
    Xyce::Util::newExpression v2exp(std::string("spice_pulse(0, 1, 0, 0.5us, 0.5us, 2us, 20us) " ), grp); v2exp.lexAndParseExpression();

    // in the original test,
    // V(1) is a sinewave that goes between +20 and -20
    // V(2) is a pulsed source that goes between 0 and 1.  PW is short.
    // The expression, V(2) * (V(1) + 30) / 60  has roughly the scaled shape of V(1) but is spiked/digitized.
    // When V(2) is zero, so is the expression.  When V(2) is 1, then expression is (V(1)+30)/60  max = 50/60, min=10/60
    //
    // The table itself is a digitized signal; like stairsteps, with 64 points.  It goes from 0 to 31 in increments of 1
    // The result is that the interpolated final output is really similar to the expression.
    int numpoints=10;
    double tfinal = 0.0005;
    double dt = tfinal/(numpoints-1), time=0.0;
    std::vector<double> refRes(numpoints), result(numpoints);
    std::vector<double> copyResult(numpoints), assignResult(numpoints);
    for (int ii=0;ii<numpoints;ii++,time+=dt)
    {
      bsrc_C1_grp->setTime(time);
      double v1Value(0.0),v2Value(0.0);
      v1exp.evaluateFunction(v1Value);
      v2exp.evaluateFunction(v2Value);
      bsrc_C1_grp->setSoln(std::string("1"),v1Value);
      bsrc_C1_grp->setSoln(std::string("2"),v2Value);
      BE_Dig.evaluateFunction(result[ii]);
      copy_BE_Dig.evaluateFunction(copyResult[ii]);
      assign_BE_Dig.evaluateFunction(assignResult[ii]);
      BE_Dig_pureArray.evaluateFunction(refRes[ii]);
    }

    EXPECT_EQ(refRes,result);
    EXPECT_EQ(refRes,copyResult);
    EXPECT_EQ(refRes,assignResult);
  }
}

#if 0
// commenting this out for now, as it doesn't work yet.  It needs to, as the 
// original Bsrc_C1.cir regression test specifies the table with parentheses 
// around each time,value pair.
TEST ( Double_Parser_table_Test, Bsrc_C1_pairsWithParens)
{
  Teuchos::RCP<Bsrc_C1_ExpressionGroup> bsrc_C1_grp = Teuchos::rcp(new Bsrc_C1_ExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  grp = bsrc_C1_grp;
  {
    // this is a nice test b/c it has curly braces around the first expression, which is one of the supported formats
    Xyce::Util::newExpression BE_Dig(std::string("TABLE { V(2) * (V(1) + 30) / 60 }  (0.0000000, 0)  (0.0312500, 0)  (0.0312813, 1)  (0.0625000, 1)  (0.0625313, 2)  (0.0937500, 2)  (0.0937813, 3)  (0.1250000, 3)  (0.1250313, 4)  (0.1562500, 4)  (0.1562813, 5)  (0.1875000, 5)  (0.1875313, 6)  (0.2187500, 6)  (0.2187813, 7)  (0.2500000, 7)  (0.2500313, 8)  (0.2812500, 8)  (0.2812813, 9)  (0.3125000, 9)  (0.3125313, 10)  (0.3437500, 10)  (0.3437813, 11)  (0.3750000, 11)  (0.3750313, 12)  (0.4062500, 12)  (0.4062813, 13)  (0.4375000, 13)  (0.4375313, 14)  (0.4687500, 14)  (0.4687813, 15)  (0.5000000, 15)  (0.5000313, 16)  (0.5312500, 16)  (0.5312813, 17)  (0.5625000, 17)  (0.5625313, 18)  (0.5937500, 18)  (0.5937813, 19)  (0.6250000, 19)  (0.6250313, 20)  (0.6562500, 20)  (0.6562813, 21)  (0.6875000, 21)  (0.6875313, 22)  (0.7187500, 22)  (0.7187813, 23)  (0.7500000, 23)  (0.7500313, 24)  (0.7812500, 24)  (0.7812813, 25)  (0.8125000, 25)  (0.8125313, 26)  (0.8437500, 26)  (0.8437813, 27)  (0.8750000, 27)  (0.8750313, 28)  (0.9062500, 28)  (0.9062813, 29)  (0.9375000, 29)  (0.9375313, 30)  (0.9687500, 30)  (0.9687813, 31)  (1.0000000, 31)"), grp);
    BE_Dig.lexAndParseExpression();

    Xyce::Util::newExpression copy_BE_Dig(BE_Dig); 
    Xyce::Util::newExpression assign_BE_Dig; 
    assign_BE_Dig = BE_Dig; 

    Xyce::Util::newExpression BE_Dig_leftArg(std::string("V(2) * (V(1) + 30) / 60"),grp);
    BE_Dig_leftArg.lexAndParseExpression();

    std::vector<double> xa = { 0.0000000, 0.0312500, 0.0312813, 0.0625000, 0.0625313, 0.0937500, 0.0937813, 0.1250000, 0.1250313, 0.1562500, 0.1562813, 0.1875000, 0.1875313, 0.2187500, 0.2187813, 0.2500000, 0.2500313, 0.2812500, 0.2812813, 0.3125000, 0.3125313, 0.3437500, 0.3437813, 0.3750000, 0.3750313, 0.4062500, 0.4062813, 0.4375000, 0.4375313, 0.4687500, 0.4687813, 0.5000000, 0.5000313, 0.5312500, 0.5312813, 0.5625000, 0.5625313, 0.5937500, 0.5937813, 0.6250000, 0.6250313, 0.6562500, 0.6562813, 0.6875000, 0.6875313, 0.7187500, 0.7187813, 0.7500000, 0.7500313, 0.7812500, 0.7812813, 0.8125000, 0.8125313, 0.8437500, 0.8437813, 0.8750000, 0.8750313, 0.9062500, 0.9062813, 0.9375000, 0.9375313, 0.9687500, 0.9687813, 1.0000000 };

    std::vector<double> ya = { 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31 };

    Xyce::Util::newExpression BE_Dig_pureArray(BE_Dig_leftArg.getAst(),xa,ya,grp);
    BE_Dig_pureArray.lexAndParseExpression();

    Xyce::Util::newExpression v1exp(std::string("spice_sin(0, 20, 1k, -.25e-3, 0, 0)" ), grp); v1exp.lexAndParseExpression();
    Xyce::Util::newExpression v2exp(std::string("spice_pulse(0, 1, 0, 0.5us, 0.5us, 2us, 20us) " ), grp); v2exp.lexAndParseExpression();

    // in the original test,
    // V(1) is a sinewave that goes between +20 and -20
    // V(2) is a pulsed source that goes between 0 and 1.  PW is short.
    // The expression, V(2) * (V(1) + 30) / 60  has roughly the scaled shape of V(1) but is spiked/digitized.
    // When V(2) is zero, so is the expression.  When V(2) is 1, then expression is (V(1)+30)/60  max = 50/60, min=10/60
    //
    // The table itself is a digitized signal; like stairsteps, with 64 points.  It goes from 0 to 31 in increments of 1
    // The result is that the interpolated final output is really similar to the expression.
    int numpoints=10;
    double tfinal = 0.0005;
    double dt = tfinal/(numpoints-1), time=0.0;
    std::vector<double> refRes(numpoints), result(numpoints);
    std::vector<double> copyResult(numpoints), assignResult(numpoints);
    for (int ii=0;ii<numpoints;ii++,time+=dt)
    {
      bsrc_C1_grp->setTime(time);
      double v1Value(0.0),v2Value(0.0);
      v1exp.evaluateFunction(v1Value);
      v2exp.evaluateFunction(v2Value);
      bsrc_C1_grp->setSoln(std::string("1"),v1Value);
      bsrc_C1_grp->setSoln(std::string("2"),v2Value);
      BE_Dig.evaluateFunction(result[ii]);
      copy_BE_Dig.evaluateFunction(copyResult[ii]);
      assign_BE_Dig.evaluateFunction(assignResult[ii]);
      BE_Dig_pureArray.evaluateFunction(refRes[ii]);
    }

    EXPECT_EQ(refRes,result);
    EXPECT_EQ(refRes,copyResult);
    EXPECT_EQ(refRes,assignResult);
  }
}
#endif

// this form of test1 doesn't rely on the group to resolve the parameter.
// Instead, it allows the user to attach it.
TEST ( Double_Parser_Param_Test, test1)
{
  Teuchos::RCP<testExpressionGroup> noparamGroup = Teuchos::rcp(new testExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = noparamGroup;

  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("2+3"), testGroup));
  p1Expression->lexAndParseExpression();
  std::string p1Name = "p1";

  Xyce::Util::newExpression testExpression(std::string("p1"), testGroup);
  testExpression.lexAndParseExpression();
  testExpression.attachParameterNode(p1Name,p1Expression);

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  double result;
  testExpression.evaluateFunction(result);        EXPECT_EQ( result, 5.0 );
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( result, 5.0 );
  OUTPUT_MACRO(Double_Parser_Param_Test, test1)
}

//-----------------------------------------------------------------------------
// This tests the use of solution variables inside a parameter.
// It is also derived from the Bsrc_C1 table test, only without the table.
// I've added the twist that p1 is multiplied by 2 in the final expression.
// This one works.

#if 0
TEST ( Double_Parser_Param_Test, test2)
{
  Teuchos::RCP<Bsrc_C1_ExpressionGroup> paramGroup = Teuchos::rcp(new Bsrc_C1_ExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = paramGroup;

  Xyce::Util::newExpression v1exp(std::string("spice_sin(0, 20, 1k, -.25e-3, 0, 0)" ), grp);            v1exp.lexAndParseExpression();
  Xyce::Util::newExpression v2exp(std::string("spice_pulse(0, 1, 0, 0.5us, 0.5us, 2us, 20us) " ), grp); v2exp.lexAndParseExpression();
  Xyce::Util::newExpression testExpression(std::string("2*p1"), grp);                                   testExpression.lexAndParseExpression();
  Teuchos::RCP<Xyce::Util::newExpression> p1exp
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("V(2) * (V(1) + 30) / 60" ), grp));
  p1exp->lexAndParseExpression();
  std::string p1Name="p1";
  paramGroup->addParam(p1Name, p1exp);

  testExpression.resolveExpression();

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  int numpoints=100;
  double tfinal = 0.0005;
  double dt = tfinal/(numpoints-1), time=0.0;

  std::vector<double> refRes(numpoints), result(numpoints);
  std::vector<double> copyResult(numpoints), assignResult(numpoints);
  for (int ii=0;ii<numpoints;ii++,time+=dt)
  {
    paramGroup->setTime(time);
    double v1Value(0.0),v2Value(0.0);
    v1exp.evaluateFunction(v1Value);
    v2exp.evaluateFunction(v2Value);
    paramGroup->setSoln(std::string("1"),v1Value);
    paramGroup->setSoln(std::string("2"),v2Value);
    testExpression.evaluateFunction(result[ii]);
    copy_testExpression.evaluateFunction(copyResult[ii]);
    assign_testExpression.evaluateFunction(assignResult[ii]);
    refRes[ii] = 2 * v2Value * (v1Value + 30) / 60;
  }
  EXPECT_EQ(refRes,result);
  EXPECT_EQ(refRes,copyResult);
  EXPECT_EQ(refRes,assignResult);
  OUTPUT_MACRO(Double_Parser_Param_Test, test2)
}
#endif


TEST ( Double_Parser_Param_Test, test2)
{
  Teuchos::RCP<Bsrc_C1_ExpressionGroup_noparam> paramGroup = Teuchos::rcp(new Bsrc_C1_ExpressionGroup_noparam() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = paramGroup;

  Xyce::Util::newExpression v1exp(std::string("spice_sin(0, 20, 1k, -.25e-3, 0, 0)" ), grp);            v1exp.lexAndParseExpression();
  Xyce::Util::newExpression v2exp(std::string("spice_pulse(0, 1, 0, 0.5us, 0.5us, 2us, 20us) " ), grp); v2exp.lexAndParseExpression();
  Xyce::Util::newExpression testExpression(std::string("2*p1"), grp);                                   testExpression.lexAndParseExpression();
  Teuchos::RCP<Xyce::Util::newExpression> p1exp
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("V(2) * (V(1) + 30) / 60" ), grp));
  p1exp->lexAndParseExpression();
  std::string p1Name="p1";
  testExpression.attachParameterNode(p1Name,p1exp);

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  int numpoints=100;
  double tfinal = 0.0005;
  double dt = tfinal/(numpoints-1), time=0.0;

  std::vector<double> refRes(numpoints), result(numpoints);
  std::vector<double> copyResult(numpoints), assignResult(numpoints);
  for (int ii=0;ii<numpoints;ii++,time+=dt)
  {
    paramGroup->setTime(time);
    double v1Value(0.0),v2Value(0.0);
    v1exp.evaluateFunction(v1Value);
    v2exp.evaluateFunction(v2Value);
    paramGroup->setSoln(std::string("1"),v1Value);
    paramGroup->setSoln(std::string("2"),v2Value);
    testExpression.evaluateFunction(result[ii]);
    copy_testExpression.evaluateFunction(copyResult[ii]);
    assign_testExpression.evaluateFunction(assignResult[ii]);
    refRes[ii] = 2 * v2Value * (v1Value + 30) / 60;
  }
  EXPECT_EQ(refRes,result);
  EXPECT_EQ(refRes,copyResult);
  EXPECT_EQ(refRes,assignResult);
  OUTPUT_MACRO(Double_Parser_Param_Test, test2)
}


TEST ( Double_Parser_calculus, ddx1)
{
  //Teuchos::RCP<testExpressionGroupWithParamSupport> paramGroup = Teuchos::rcp(new testExpressionGroupWithParamSupport() );
  Teuchos::RCP<testExpressionGroup> paramGroup = Teuchos::rcp(new testExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;
  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2+3"), testGroup));
  p1Expression->lexAndParseExpression();
  std::string p1Name="p1";

  Xyce::Util::newExpression ddxTest(std::string("ddx(2*p1,p1)"), testGroup); ddxTest.lexAndParseExpression();
  ddxTest.attachParameterNode(p1Name,p1Expression);

  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  double result;
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result, 2.0 );
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result, 2.0 );
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result, 2.0 );
}

TEST ( Double_Parser_calculus, ddx2)
{
  //Teuchos::RCP<testExpressionGroupWithParamSupport> paramGroup = Teuchos::rcp(new testExpressionGroupWithParamSupport() );
  Teuchos::RCP<testExpressionGroup> paramGroup = Teuchos::rcp(new testExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;
  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("2+3"), testGroup));
  p1Expression->lexAndParseExpression();
  std::string p1Name="p1";

  Xyce::Util::newExpression ddxTest(std::string("ddx(p1*p1,p1)"), testGroup); ddxTest.lexAndParseExpression();
  ddxTest.attachParameterNode(p1Name,p1Expression);

  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  double result;
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result, 2.0*(2+3) );
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result, 2.0*(2+3) );
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result, 2.0*(2+3) );
}

TEST ( Double_Parser_calculus, ddx3)
{
  //Teuchos::RCP<testExpressionGroupWithParamSupport> paramGroup = Teuchos::rcp(new testExpressionGroupWithParamSupport() );
  Teuchos::RCP<testExpressionGroup> paramGroup = Teuchos::rcp(new testExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;
  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2+3"), testGroup));
  p1Expression->lexAndParseExpression();
  std::string p1Name="p1";

  Xyce::Util::newExpression ddxTest(std::string("ddx(sin(p1),p1)"), testGroup); ddxTest.lexAndParseExpression();
  ddxTest.attachParameterNode(p1Name,p1Expression);

  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  double result;
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result, std::cos(2+3) );
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result, std::cos(2+3) );
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result, std::cos(2+3) );
}

TEST ( Double_Parser_calculus, ddx4)
{
  Teuchos::RCP<testExpressionGroup> paramGroup = Teuchos::rcp(new testExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;
  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2+3"), testGroup)); p1Expression->lexAndParseExpression();
  std::string p1Name="p1";

  Xyce::Util::newExpression ddxTest(std::string("ddx(sin(p1*p1),p1)"), testGroup); ddxTest.lexAndParseExpression();
  ddxTest.attachParameterNode(p1Name,p1Expression);

  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  double result;
  double p1 = 2+3;
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result, 2.0*p1*std::cos(p1*p1) );
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result, 2.0*p1*std::cos(p1*p1) );
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result, 2.0*p1*std::cos(p1*p1) );
}

TEST ( Double_Parser_calculus, ddx5)
{
  //Teuchos::RCP<testExpressionGroupWithParamSupport> paramGroup = Teuchos::rcp(new testExpressionGroupWithParamSupport() );
  Teuchos::RCP<testExpressionGroup> paramGroup = Teuchos::rcp(new testExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;
  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2+3"), testGroup)); 
  p1Expression->lexAndParseExpression();
  std::string p1Name="p1";

  Xyce::Util::newExpression ddxTest(std::string("ddx( pow(sin(p1*p1),3.0),p1)"), testGroup); ddxTest.lexAndParseExpression();
  ddxTest.attachParameterNode(p1Name,p1Expression);

  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  double result;
  double p1 = 2+3;
  double p1Sq = p1*p1;
  double refRes = 3.0*(2.0*p1*std::cos(p1Sq))/(std::sin(p1Sq))*std::pow(std::sin(p1Sq),3.0);
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result-refRes, 0.0 );
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result-refRes, 0.0 );
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result-refRes, 0.0 );
}

TEST ( Double_Parser_calculus, ddx5b)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression ddxTest(std::string("pow(sin(V(A)*V(A)),3.0)"), testGroup); ddxTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  double Aval=5.0;
  solnGroup->setSoln(std::string("A"),Aval);
  double result;
  std::vector<double> derivs;
  double p1 = 5;
  double p1Sq = p1*p1;
  double refRes = 3.0*(2.0*p1*std::cos(p1Sq))/(std::sin(p1Sq))*std::pow(std::sin(p1Sq),3.0);
  std::vector<double> refderivs = { refRes };

  ddxTest.evaluate(result,derivs);        EXPECT_EQ( derivs, refderivs );
  copy_ddxTest.evaluate(result,derivs);   EXPECT_EQ( derivs, refderivs );
  assign_ddxTest.evaluate(result,derivs); EXPECT_EQ( derivs, refderivs );
} 

TEST ( Double_Parser_calculus, ddx6)
{
  //Teuchos::RCP<testExpressionGroupWithParamSupport> paramGroup = Teuchos::rcp(new testExpressionGroupWithParamSupport() );
  Teuchos::RCP<testExpressionGroup> paramGroup = Teuchos::rcp(new testExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;
  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2+3"), testGroup)); 
  p1Expression->lexAndParseExpression();
  std::string p1Name="p1";

  Xyce::Util::newExpression ddxTest(std::string("ddx( pow(p1,3.0),p1)"), testGroup); ddxTest.lexAndParseExpression();
  ddxTest.attachParameterNode(p1Name,p1Expression);

  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  double result;
  double p1 = 2+3;
  double refRes = 3.0*p1*p1;
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result, refRes );
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result, refRes );
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result, refRes );
} 

TEST ( Double_Parser_calculus, ddx7)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression ddxTest(std::string("ddx(sin(v(a)),v(a))"), testGroup); ddxTest.lexAndParseExpression();
  
  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  solnGroup->setSoln(std::string("A"),5.0);
  double result;
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result, std::cos(5.0) );
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result, std::cos(5.0) );
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result, std::cos(5.0) );
}

#if 0
TEST ( Double_Parser_calculus, ddx8)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression ddxTest(std::string("ddx(sin(v(a,b)),v(a,b))"), testGroup); ddxTest.lexAndParseExpression();
  
  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  solnGroup->setSoln(std::string("A"),5.0);
  solnGroup->setSoln(std::string("B"),3.0);
  double result;
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result, std::cos(2.0) );
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result, std::cos(2.0) );
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result, std::cos(2.0) );
}
#endif

TEST ( Double_Parser_calculus, ddx9)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression ddxTest(std::string("ddx(sin(i(a)),i(a))"), testGroup); ddxTest.lexAndParseExpression();
  
  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  solnGroup->setSoln(std::string("A"),5.0);
  double result;
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result, std::cos(5.0) );
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result, std::cos(5.0) );
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result, std::cos(5.0) );
}

TEST ( Double_Parser_calculus, ddx10)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression ddxTest(std::string("ddx(pow(5.0,v(a)),v(a))"), testGroup); ddxTest.lexAndParseExpression();
  
  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  double Aval=2.0;
  solnGroup->setSoln(std::string("A"),Aval);
  double result;
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result, std::log(5.0)*std::pow(5.0,2.0) );
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result, std::log(5.0)*std::pow(5.0,2.0) );
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result, std::log(5.0)*std::pow(5.0,2.0) );
}

TEST ( Double_Parser_calculus, ddx11)
{
  //Teuchos::RCP<testExpressionGroupWithParamSupport> paramGroup = Teuchos::rcp(new testExpressionGroupWithParamSupport() );
  Teuchos::RCP<testExpressionGroup> paramGroup = Teuchos::rcp(new testExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;

  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2"), testGroup)); 
  p1Expression->lexAndParseExpression();
  std::string p1Name="p1";

  Xyce::Util::newExpression ddxTest(std::string("ddx( pow(sin(p1),p1),p1)"), testGroup); ddxTest.lexAndParseExpression();
  ddxTest.attachParameterNode(p1Name,p1Expression);
  
  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  double Aval=2.0;
  double result;
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result, std::pow(std::sin(Aval),Aval)*(Aval/std::tan(Aval) + std::log(sin(Aval))) );
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result, std::pow(std::sin(Aval),Aval)*(Aval/std::tan(Aval) + std::log(sin(Aval))) );
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result, std::pow(std::sin(Aval),Aval)*(Aval/std::tan(Aval) + std::log(sin(Aval))) );
}

TEST ( Double_Parser_calculus, simpleDerivs1 )
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression ddxTest(std::string("0.5*(V(B)-3.0)**2.0"), testGroup); 
  ddxTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  double Aval=2.5;
  solnGroup->setSoln(std::string("B"),Aval);
  double result;
  std::vector<double> derivs;
  double refRes = 1.25e-01; 
  std::vector<double> refderivs = { -0.5 };

  ddxTest.evaluate(result,derivs);        EXPECT_EQ( derivs, refderivs );
  copy_ddxTest.evaluate(result,derivs);   EXPECT_EQ( derivs, refderivs );
  assign_ddxTest.evaluate(result,derivs); EXPECT_EQ( derivs, refderivs );
}

class solnAndFuncExpressionGroup : public Xyce::Util::baseExpressionGroup
{
  public:
    solnAndFuncExpressionGroup () :
      Xyce::Util::baseExpressionGroup(), Aval_(0.0), Bval_(0.0), Cval_(0.0), R1val_(0.0)  {};
    ~solnAndFuncExpressionGroup () {};

  virtual bool getSolutionVal(const std::string & nodeName, double & retval )
  {
    std::string tmp = nodeName; Xyce::Util::toLower(tmp);
    if (tmp==std::string("a")) { retval = Aval_; return true; }
    else if (tmp==std::string("b")) { retval = Bval_; return true; }
    else if (tmp==std::string("c")) { retval = Cval_; return true; }
    else if (tmp==std::string("r1")) { retval = R1val_; return true; }
    else { return 0.0; return false; }
  }

  void setSoln(const std::string & nodeName, double val)
  {
    std::string tmp = nodeName; Xyce::Util::toLower(tmp);
    if (tmp==std::string("a")) { Aval_ = val; }
    else if (tmp==std::string("b")) { Bval_ = val; }
    else if (tmp==std::string("c")) { Cval_ = val; }
    else if (tmp==std::string("r1")) { R1val_ = val; }
  }

  private:
    double Aval_, Bval_, Cval_, R1val_;
};

// These tests (derivsThruFuncs?) tests if derivatives work thru expression arguments.
// At the time of test creation (2/21/2020), the answer was NO.
//
TEST ( Double_Parser_calculus, derivsThruFuncs1 )
{
  Teuchos::RCP<solnAndFuncExpressionGroup> solnFuncGroup = Teuchos::rcp(new solnAndFuncExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnFuncGroup;

  // this expression is the RHS of a .func statement:  .func F1(A,B) {A-B}
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression (std::string("A-B"), testGroup));
  std::vector<std::string> f1ArgStrings;

  Xyce::Util::newExpression f1_LHS (std::string("F1(A,B)"), testGroup);
  f1_LHS.lexAndParseExpression();
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);
  f1Expression->setFunctionArgStringVec ( f1ArgStrings );
  // during lex/parse, this vector of arg strings will be compared to any
  // param classes.  If it finds them, then they will be placed in the
  // functionArgOpVec object, which is used below, in the call to "setFuncArgs".
  f1Expression->lexAndParseExpression();

  std::string f1Name;// = "F1";
  f1_LHS.getFuncPrototypeName(f1Name);

  Xyce::Util::newExpression derivFuncTestExpr(std::string("0.5*(F1(V(B),3.0))**2.0"), testGroup); 
  derivFuncTestExpr.lexAndParseExpression();
  derivFuncTestExpr.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copy_derivFuncTestExpr(derivFuncTestExpr); 
  Xyce::Util::newExpression assign_derivFuncTestExpr; 
  assign_derivFuncTestExpr = derivFuncTestExpr; 

  double Bval=2.5;
  solnFuncGroup->setSoln(std::string("B"),Bval);
  double result;
  std::vector<double> derivs;
  double refRes = 1.25e-01; 
  std::vector<double> refderivs = { -0.5 };

  derivFuncTestExpr.evaluate(result,derivs);        EXPECT_EQ( derivs, refderivs );
  copy_derivFuncTestExpr.evaluate(result,derivs);   EXPECT_EQ( derivs, refderivs );
  assign_derivFuncTestExpr.evaluate(result,derivs); EXPECT_EQ( derivs, refderivs );
}

// there are a bunch of tests in this one, all testing if derivatives propagate thru
// a .FUNC in the right way.
TEST ( Double_Parser_calculus, derivsThruFuncs2 )
{
  Teuchos::RCP<solnAndFuncExpressionGroup> solnFuncGroup = Teuchos::rcp(new solnAndFuncExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnFuncGroup;

  // this expression is the LHS of a .func statement:  .func F1(A,B) {sin(A)*cos(B)}
  std::vector<std::string> f1ArgStrings;// = { std::string("A"), std::string("B") };
  Xyce::Util::newExpression f1_LHS (std::string("F1(A,B)"), testGroup);
  f1_LHS.lexAndParseExpression();
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);

  // this expression is the RHS of a .func statement:  .func F1(A,B) {sin(A)*cos(B)}
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("sin(A)*cos(B)"), testGroup));
  f1Expression->setFunctionArgStringVec ( f1ArgStrings );
  f1Expression->lexAndParseExpression();

  // this expression is the LHS of a .func statement:  .func F2(A,B) {A*B}
  std::vector<std::string> f2ArgStrings;// = { std::string("A"), std::string("B") };
  Xyce::Util::newExpression f2_LHS (std::string("F2(A,B)"), testGroup);
  f2_LHS.lexAndParseExpression();
  f2_LHS.getFuncPrototypeArgStrings(f2ArgStrings);

  // this expression is the RHS of a .func statement:  .func F2(A,B) {A*B}
  Teuchos::RCP<Xyce::Util::newExpression> f2Expression  = Teuchos::rcp(new Xyce::Util::newExpression (std::string("A*B"), testGroup));
  f2Expression->setFunctionArgStringVec ( f2ArgStrings );
  f2Expression->lexAndParseExpression();

  Xyce::Util::newExpression derivFuncTestExpr1(std::string("0.5*(F1(V(A),V(B)))**2.0"), testGroup); 
  Xyce::Util::newExpression derivFuncTestExpr2(std::string("0.5*(F2(sin(V(A)),cos(V(B))))**2.0"), testGroup); 

  std::string f1Name, f2Name;
  f1_LHS.getFuncPrototypeName(f1Name);
  f2_LHS.getFuncPrototypeName(f2Name);

  derivFuncTestExpr1.lexAndParseExpression();
  derivFuncTestExpr1.attachFunctionNode(f1Name, f1Expression);
  derivFuncTestExpr1.attachFunctionNode(f2Name, f2Expression);

  derivFuncTestExpr2.lexAndParseExpression();
  derivFuncTestExpr2.attachFunctionNode(f1Name, f1Expression);
  derivFuncTestExpr2.attachFunctionNode(f2Name, f2Expression);

  Xyce::Util::newExpression copy_derivFuncTestExpr1(derivFuncTestExpr1); 
  Xyce::Util::newExpression assign_derivFuncTestExpr1; 
  assign_derivFuncTestExpr1 = derivFuncTestExpr1; 

  Xyce::Util::newExpression copy_derivFuncTestExpr2(derivFuncTestExpr2); 
  Xyce::Util::newExpression assign_derivFuncTestExpr2; 
  assign_derivFuncTestExpr2 = derivFuncTestExpr2; 

  double Aval=0.45;
  double Bval=0.6;
  solnFuncGroup->setSoln(std::string("A"),Aval);
  solnFuncGroup->setSoln(std::string("B"),Bval);
  double result;
  double result2;
  std::vector<double> derivs;
  std::vector<double> derivs2;

  // analytic answer: seems to have a roundoff problem compared to computed result.
  // I can make the expression library match analytic result, but only for 
  // certain Aval, Bval values.
  double refRes;
  std::vector<double> refderivs;
  {
    double f1val = std::sin(Aval)*std::cos(Bval);
    refRes =  0.5*f1val*f1val;
    double df1_dA = +std::cos(Aval)*std::cos(Bval);
    double df1_dB = -std::sin(Aval)*std::sin(Bval);
    double dExp_df1 = f1val;
    double dExp_dA = df1_dA*f1val;
    double dExp_dB = df1_dB*f1val;
    refderivs = { dExp_dA, dExp_dB };
  }

  derivFuncTestExpr1.evaluate(result,derivs);        
  derivFuncTestExpr2.evaluate(result2,derivs2);        
  EXPECT_EQ( result, result2 );
  EXPECT_EQ( derivs.size(), derivs2.size() );
  if ( derivs.size() == derivs2.size() ) { EXPECT_EQ( derivs, derivs2 ); }

  EXPECT_EQ( result-refRes, 0.0 );
  EXPECT_EQ(derivs.size(),refderivs.size());
  if ( derivs.size() == refderivs.size() )
  {
    std::vector<double> derivDiffs = { (derivs[0]-refderivs[0]),  (derivs[1]-refderivs[1]) };
    EXPECT_EQ( derivDiffs[0], 0.0 );
    EXPECT_EQ( derivDiffs[1], 0.0 );
  }

  copy_derivFuncTestExpr1.evaluate(result,derivs);   
  EXPECT_EQ( result, result2 );
  EXPECT_EQ( derivs.size(), derivs2.size() );
  if ( derivs.size() == derivs2.size() ) { EXPECT_EQ( derivs, derivs2 ); }

  assign_derivFuncTestExpr1.evaluate(result,derivs); 
  EXPECT_EQ( result, result2 );
  EXPECT_EQ( derivs.size(), derivs2.size() );
  if ( derivs.size() == derivs2.size() ) { EXPECT_EQ( derivs, derivs2 ); }
}

TEST ( Double_Parser_calculus, derivsThruFuncs3 )
{
  Teuchos::RCP<solnAndFuncExpressionGroup> solnFuncGroup = Teuchos::rcp(new solnAndFuncExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnFuncGroup;

  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression (std::string("A*B"), testGroup));
  std::vector<std::string> f1ArgStrings = { std::string("A"), std::string("B") };
  f1Expression->setFunctionArgStringVec ( f1ArgStrings );
  f1Expression->lexAndParseExpression();

  Xyce::Util::newExpression derivFuncTestExpr1(std::string("V(A)*F1(V(A),V(B)*V(B))+3.0*V(B)"), testGroup); 

  derivFuncTestExpr1.lexAndParseExpression();
  derivFuncTestExpr1.attachFunctionNode(std::string("F1"), f1Expression);

  Xyce::Util::newExpression copy_derivFuncTestExpr1(derivFuncTestExpr1); 
  Xyce::Util::newExpression assign_derivFuncTestExpr1; 
  assign_derivFuncTestExpr1 = derivFuncTestExpr1; 

  double Aval=1.0;
  double Bval=2.0;
  solnFuncGroup->setSoln(std::string("A"),Aval);
  solnFuncGroup->setSoln(std::string("B"),Bval);
  double result;
  std::vector<double> derivs;

  derivFuncTestExpr1.evaluate(result,derivs);   

  double resRef = Aval*Aval*Bval*Bval+3.0*Bval;
  double dfdA = (2.0*Aval*Bval*Bval);
  double dfdB = (2.0*Aval*Aval*Bval + 3.0);

  std::vector<double> derivsRef = { dfdA, dfdB };

  EXPECT_EQ( result-resRef, 0.0 );

  std::vector<double> derivDiffs = { (derivs[0]-derivsRef[0]),  (derivs[1]-derivsRef[1]) };
  EXPECT_EQ( derivDiffs[0], 0.0 );
  EXPECT_EQ( derivDiffs[1], 0.0 );

}

TEST ( Double_Parser_calculus, derivsThruFuncs4 )
{
  Teuchos::RCP<solnAndFuncExpressionGroup> solnFuncGroup = Teuchos::rcp(new solnAndFuncExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnFuncGroup;

  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("A*B+100*V(C)"), testGroup));
  std::vector<std::string> f1ArgStrings = { std::string("A"), std::string("B") };
  f1Expression->setFunctionArgStringVec ( f1ArgStrings );
  f1Expression->lexAndParseExpression();

  Xyce::Util::newExpression derivFuncTestExpr1(std::string("F1(V(A),V(B))"), testGroup); 

  derivFuncTestExpr1.lexAndParseExpression();
  derivFuncTestExpr1.attachFunctionNode(std::string("F1"), f1Expression);

  Xyce::Util::newExpression copy_derivFuncTestExpr1(derivFuncTestExpr1); 
  Xyce::Util::newExpression assign_derivFuncTestExpr1; 
  assign_derivFuncTestExpr1 = derivFuncTestExpr1; 

  double Aval=17.0;
  double Bval=2.0;
  double Cval=3.0;
  solnFuncGroup->setSoln(std::string("A"),Aval);
  solnFuncGroup->setSoln(std::string("B"),Bval);
  solnFuncGroup->setSoln(std::string("C"),Cval);
  double result;
  std::vector<double> derivs;

  derivFuncTestExpr1.evaluate(result,derivs);   

  double resRef = (Aval*Bval+100*Cval);
  std::vector<double> derivsRef = { Bval, Aval, 100.0 };

  EXPECT_EQ( result,resRef);
  EXPECT_EQ( derivs,derivsRef);
}

TEST ( Double_Parser_calculus, derivsThruFuncs5 )
{
  Teuchos::RCP<solnAndFuncExpressionGroup> solnFuncGroup = Teuchos::rcp(new solnAndFuncExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnFuncGroup;

  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("A*B+100*V(A)"), testGroup));
  std::vector<std::string> f1ArgStrings = { std::string("A"), std::string("B") };
  f1Expression->setFunctionArgStringVec ( f1ArgStrings );
  f1Expression->lexAndParseExpression();

  //solnFuncGroup->addFunction(std::string("F1"), f1Expression);

  Xyce::Util::newExpression derivFuncTestExpr1(std::string("F1(V(A),V(B))"), testGroup); 

  derivFuncTestExpr1.lexAndParseExpression();
  derivFuncTestExpr1.attachFunctionNode(std::string("F1"), f1Expression);

  //derivFuncTestExpr1.dumpParseTree(std::cout);

  Xyce::Util::newExpression copy_derivFuncTestExpr1(derivFuncTestExpr1); 
  Xyce::Util::newExpression assign_derivFuncTestExpr1; 
  assign_derivFuncTestExpr1 = derivFuncTestExpr1; 

  double Aval=17.0;
  double Bval=2.0;
  solnFuncGroup->setSoln(std::string("A"),Aval);
  solnFuncGroup->setSoln(std::string("B"),Bval);
  double result;
  std::vector<double> derivs;

  derivFuncTestExpr1.evaluate(result,derivs);   

  double resRef = (Aval*Bval+100.0*Aval);
  std::vector<double> derivsRef = { (Bval+100.0), Aval };

  EXPECT_EQ( result,resRef);
  EXPECT_EQ( derivs,derivsRef);
}

//class solnAndParamExpressionGroup : public Xyce::Util::baseExpressionGroup
class solnExpressionGroup2 : public Xyce::Util::baseExpressionGroup
{
  public:
    solnExpressionGroup2 () :
      Xyce::Util::baseExpressionGroup(), Aval_(0.0), Bval_(0.0), Cval_(0.0), R1val_(0.0)  {};
    ~solnExpressionGroup2 () {};

  virtual bool getSolutionVal(const std::string & nodeName, double & retval )
  {
    std::string tmp = nodeName; Xyce::Util::toLower(tmp);
    if (tmp==std::string("a")) { retval = Aval_; return true; }
    else if (tmp==std::string("b")) { retval = Bval_; return true; }
    else if (tmp==std::string("c")) { retval = Cval_; return true; }
    else if (tmp==std::string("r1")) { retval = R1val_; return true; }
    else { return 0.0; return false; }
  }

  void setSoln(const std::string & nodeName, double val)
  {
    std::string tmp = nodeName; Xyce::Util::toLower(tmp);
    if (tmp==std::string("a")) { Aval_ = val; }
    else if (tmp==std::string("b")) { Bval_ = val; }
    else if (tmp==std::string("c")) { Cval_ = val; }
    else if (tmp==std::string("r1")) { R1val_ = val; }
  }

#if 0
  void addParam (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
  {
    std::string lowerName = name;
    Xyce::Util::toLower(lowerName);

    params_[lowerName] = exp;
  };

  bool getParam (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
  {
    bool retval=true;

    std::string lowerName = name;
    Xyce::Util::toLower(lowerName);

    if (params_.find(lowerName) != params_.end()) { exp = params_[lowerName]; }
    else { retval = false; }

    return retval;
  }
#endif

  private:
    //std::unordered_map <std::string, Teuchos::RCP<Xyce::Util::newExpression> >  params_;
    double Aval_, Bval_, Cval_, R1val_;
};

TEST ( Double_Parser_calculus, derivsThruParams1 )
{
  //Teuchos::RCP<solnAndParamExpressionGroup> solnParamGroup = Teuchos::rcp(new solnAndParamExpressionGroup() );
  Teuchos::RCP<solnExpressionGroup2> solnParamGroup = Teuchos::rcp(new solnExpressionGroup2() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnParamGroup;

  // this expression is the RHS of a .global_param statement:  .global_param P1 {V(A)-V(B)}
  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("V(A)-V(B)"), testGroup));
  p1Expression->lexAndParseExpression();

  std::string p1Name = "P1";
  //solnParamGroup->addParam( p1Name ,  p1Expression);

  Xyce::Util::newExpression derivParamTestExpr(std::string("0.5*P1"), testGroup); 
  derivParamTestExpr.lexAndParseExpression();
  derivParamTestExpr.attachParameterNode(p1Name , p1Expression);

  Xyce::Util::newExpression copy_derivParamTestExpr(derivParamTestExpr); 
  Xyce::Util::newExpression assign_derivParamTestExpr; 
  assign_derivParamTestExpr = derivParamTestExpr; 

  double Aval=7.5;
  double Bval=2.5;
  solnParamGroup->setSoln(std::string("A"),Aval);
  solnParamGroup->setSoln(std::string("B"),Bval);
  double result;
  std::vector<double> derivs;
  double refRes = 2.5;
  std::vector<double> refderivs = { 0.5, -0.5 };

  derivParamTestExpr.evaluate(result,derivs);        EXPECT_EQ( derivs, refderivs );
  copy_derivParamTestExpr.evaluate(result,derivs);   EXPECT_EQ( derivs, refderivs );
  assign_derivParamTestExpr.evaluate(result,derivs); EXPECT_EQ( derivs, refderivs );
}

TEST ( Double_Parser_floor, test1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression floorTest(std::string("floor(10.25)"), testGroup);
  floorTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_floorTest(floorTest); 
  Xyce::Util::newExpression assign_floorTest; 
  assign_floorTest = floorTest; 

  double result;
  floorTest.evaluateFunction(result);        EXPECT_EQ( result, 10);
  copy_floorTest.evaluateFunction(result);   EXPECT_EQ( result, 10);
  assign_floorTest.evaluateFunction(result); EXPECT_EQ( result, 10);
  OUTPUT_MACRO2(Double_Parser_floor, test1, floorTest) 
}

TEST ( Double_Parser_floor, test2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression floorTest(std::string("floor(-34.251)"), testGroup);
  floorTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_floorTest(floorTest); 
  Xyce::Util::newExpression assign_floorTest; 
  assign_floorTest = floorTest; 

  double result;
  floorTest.evaluateFunction(result);        EXPECT_EQ( result, -35);
  copy_floorTest.evaluateFunction(result);   EXPECT_EQ( result, -35);
  assign_floorTest.evaluateFunction(result); EXPECT_EQ( result, -35);
  OUTPUT_MACRO2(Double_Parser_floor, test2, floorTest) 
}

TEST ( Double_Parser_floor, test3)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression floorTest(std::string("floor(0.71)"), testGroup);
  floorTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_floorTest(floorTest); 
  Xyce::Util::newExpression assign_floorTest; 
  assign_floorTest = floorTest; 

  double result;
  floorTest.evaluateFunction(result);        EXPECT_EQ( result, 0);
  copy_floorTest.evaluateFunction(result);   EXPECT_EQ( result, 0);
  assign_floorTest.evaluateFunction(result); EXPECT_EQ( result, 0);
  OUTPUT_MACRO2(Double_Parser_floor, test3, floorTest) 
}

TEST ( Double_Parser_ceil, test1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression ceilTest(std::string("ceil(10.25)"), testGroup);
  ceilTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_ceilTest(ceilTest); 
  Xyce::Util::newExpression assign_ceilTest; 
  assign_ceilTest = ceilTest; 

  double result;
  ceilTest.evaluateFunction(result);        EXPECT_EQ( result, 11);
  copy_ceilTest.evaluateFunction(result);   EXPECT_EQ( result, 11);
  assign_ceilTest.evaluateFunction(result); EXPECT_EQ( result, 11);
  OUTPUT_MACRO2(Double_Parser_ceil, test1, ceilTest) 
}

TEST ( Double_Parser_ceil, test2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression ceilTest(std::string("ceil(-34.251)"), testGroup);
  ceilTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_ceilTest(ceilTest); 
  Xyce::Util::newExpression assign_ceilTest; 
  assign_ceilTest = ceilTest; 

  double result;
  ceilTest.evaluateFunction(result);        EXPECT_EQ( result, -34);
  copy_ceilTest.evaluateFunction(result);   EXPECT_EQ( result, -34);
  assign_ceilTest.evaluateFunction(result); EXPECT_EQ( result, -34);
  OUTPUT_MACRO2(Double_Parser_ceil, test2, ceilTest) 
}

TEST ( Double_Parser_ceil, test3)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression ceilTest(std::string("ceil(0.71)"), testGroup);
  ceilTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_ceilTest(ceilTest); 
  Xyce::Util::newExpression assign_ceilTest; 
  assign_ceilTest = ceilTest; 

  double result;
  ceilTest.evaluateFunction(result);        EXPECT_EQ( result, 1);
  copy_ceilTest.evaluateFunction(result);   EXPECT_EQ( result, 1);
  assign_ceilTest.evaluateFunction(result); EXPECT_EQ( result, 1);
  OUTPUT_MACRO2(Double_Parser_ceil, test3, ceilTest) 
}

TEST ( Double_Parser_specials, pi1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression piTest(std::string("pi"), testGroup);
  piTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_piTest(piTest); 
  Xyce::Util::newExpression assign_piTest; 
  assign_piTest = piTest; 

  double result;
  piTest.evaluateFunction(result);        EXPECT_EQ( result, M_PI);
  copy_piTest.evaluateFunction(result);   EXPECT_EQ( result, M_PI);
  assign_piTest.evaluateFunction(result); EXPECT_EQ( result, M_PI);
  OUTPUT_MACRO2(Double_Parser_specials, pi1, piTest) 
}

TEST ( Double_Parser_specials, pi2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression piTest(std::string("sin(pi)"), testGroup);
  piTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_piTest(piTest); 
  Xyce::Util::newExpression assign_piTest; 
  assign_piTest = piTest; 

  double result;
  piTest.evaluateFunction(result);        EXPECT_EQ( result, std::sin(M_PI));
  copy_piTest.evaluateFunction(result);   EXPECT_EQ( result, std::sin(M_PI));
  assign_piTest.evaluateFunction(result); EXPECT_EQ( result, std::sin(M_PI));
  OUTPUT_MACRO2(Double_Parser_specials, pi2, piTest) 
}

TEST ( Double_Parser_specials, time)
{
  Teuchos::RCP<timeDepExpressionGroup> timeGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = timeGroup;
  Xyce::Util::newExpression testExpression(std::string("time"),testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  timeGroup->setTime(1.0);
  double result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( (result-(1.0)), 0.0);
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( (result-(1.0)), 0.0);
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( (result-(1.0)), 0.0);
  OUTPUT_MACRO(Double_Parser_specials, time)
}

TEST ( Double_Parser_specials, freq)
{
  Teuchos::RCP<timeDepExpressionGroup> timeGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = timeGroup;
  Xyce::Util::newExpression testExpression(std::string("freq"),testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  timeGroup->setFreq(1.0);
  double result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( (result-(1.0)), 0.0);
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( (result-(1.0)), 0.0);
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( (result-(1.0)), 0.0);
  OUTPUT_MACRO(Double_Parser_specials, freq)
}

TEST ( Double_Parser_specials, temp)
{
  Teuchos::RCP<tempDepExpressionGroup> tempGroup = Teuchos::rcp(new tempDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = tempGroup;
  Xyce::Util::newExpression testExpression(std::string("temp"),testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  tempGroup->setTemp(1.0);
  double result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( (result-(1.0)), 0.0);
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( (result-(1.0)), 0.0);
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( (result-(1.0)), 0.0);
  OUTPUT_MACRO(Double_Parser_specials, temp)
}

// these next two tests are for the use case of a parameter that is named either "I" or "V".
// In an earlier implementation, the parser would get confused by this.  The string res*I*I, 
// which is used by the regression test BUG_645_SON/bug645son_Func.cir, would simply fail to 
// parse and wouldn't even issue an error.  It would proceed but then fail to compute anything.
//

TEST ( Double_Parser_Param_Test, I )
{
  Teuchos::RCP<testExpressionGroup> paramGroup = Teuchos::rcp(new testExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;

  Teuchos::RCP<Xyce::Util::newExpression> iExpression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("2+3"), testGroup));
  iExpression->lexAndParseExpression();
  std::string iName = "I";

  Teuchos::RCP<Xyce::Util::newExpression> resExpression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("4"), testGroup));
  resExpression->lexAndParseExpression();
  std::string resName = "RES";

  Xyce::Util::newExpression testExpression(std::string("RES*I*I"), testGroup);
  testExpression.lexAndParseExpression();
  testExpression.attachParameterNode(iName,iExpression);
  testExpression.attachParameterNode(resName,resExpression);

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  double result;
  testExpression.evaluateFunction(result);        EXPECT_EQ( result, (2+3)*(2+3)*4 );
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( result, (2+3)*(2+3)*4 );
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( result, (2+3)*(2+3)*4 );
  OUTPUT_MACRO(Double_Parser_Param_Test, I)
}

TEST ( Double_Parser_Param_Test, V )
{
  Teuchos::RCP<testExpressionGroup> paramGroup = Teuchos::rcp(new testExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;

  Teuchos::RCP<Xyce::Util::newExpression> vExpression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("2+3"), testGroup));
  vExpression->lexAndParseExpression();
  std::string vName = "V";

  Teuchos::RCP<Xyce::Util::newExpression> resExpression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("4"), testGroup));
  resExpression->lexAndParseExpression();
  std::string resName = "RES";

  Xyce::Util::newExpression testExpression(std::string("RES*V*V"), testGroup);
  testExpression.lexAndParseExpression();
  testExpression.attachParameterNode(vName,vExpression);
  testExpression.attachParameterNode(resName,resExpression);

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  double result;
  testExpression.evaluateFunction(result);        EXPECT_EQ( result, (2+3)*(2+3)*4 );
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( result, (2+3)*(2+3)*4 );
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( result, (2+3)*(2+3)*4 );
  OUTPUT_MACRO(Double_Parser_Param_Test, V)
}


TEST ( Double_Parser_ASCTH_Test, test0)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("cosh(V(A))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=-10.0;
  double refRes = std::cosh(-10.0);
  solnGroup->setSoln(std::string("A"),Aval);

  std::vector<double> derivs;
  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(Double_Parser_ASCTH_Test, test0)
}

TEST ( Double_Parser_ASCTH_Test, test1)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("acosh(cosh(V(A)))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=-10.0;
  double refRes = std::acosh(std::cosh(-10.0));
  solnGroup->setSoln(std::string("A"),Aval);

  std::vector<double> derivs;
  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  copyExpression.evaluate(result,derivs);   
  EXPECT_EQ( result, refRes);
  assignExpression.evaluate(result,derivs); 
  EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(Double_Parser_ASCTH_Test, test1)
}

TEST ( Double_Parser_ASCTH_Test, test2)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("acosh(cosh(V(A)))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=0.0;
  double refRes = std::acosh(std::cosh(0.0));
  solnGroup->setSoln(std::string("A"),Aval);

  // this double checks if the derivatives are NOT Nan.
  std::vector<double> derivs;
  std::vector<double> refderivs = {0.0};
  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(Double_Parser_ASCTH_Test, test2)
}

TEST ( Double_Parser_STP_Test, test1)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("stp(V(A))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=2.0;
  double refRes = 1.0;
  solnGroup->setSoln(std::string("A"),Aval);

  std::vector<double> derivs;
  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(Double_Parser_STP_Test, test1)
}

TEST ( Double_Parser_STP_Test, test2)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("stp(V(A))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=-2.0;
  double refRes = 0.0;
  solnGroup->setSoln(std::string("A"),Aval);

  std::vector<double> derivs;
  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(Double_Parser_STP_Test, test2)
}


TEST ( Double_Parser_atan2_Test, test1)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("atan2(V(A),V(B))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=0.5, Bval=0.25;
  double refRes = std::atan2(Aval,Bval);
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("B"),Bval);

  std::vector<double> derivs;

  double denom = Aval*Aval+Bval*Bval;
  std::vector<double> refderivs = { (Bval/denom), (-Aval/denom) };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(Double_Parser_atan2_Test, test1)
}


TEST ( Double_Parser_atan2_Test, test2)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("atan2(V(A),V(B))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=0.5, Bval=-0.25;
  double refRes = std::atan2(Aval,Bval);
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("B"),Bval);

  std::vector<double> derivs;

  double denom = Aval*Aval+Bval*Bval;
  std::vector<double> refderivs = { (Bval/denom), (-Aval/denom) };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(Double_Parser_atan2_Test, test2)
}

TEST ( Double_Parser_atan2_Test, test3)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("atan2(V(A),V(B))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=-0.5, Bval=0.25;
  double refRes = std::atan2(Aval,Bval);
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("B"),Bval);

  std::vector<double> derivs;

  double denom = Aval*Aval+Bval*Bval;
  std::vector<double> refderivs = { (Bval/denom), (-Aval/denom) };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(Double_Parser_atan2_Test, test3)
}

TEST ( Double_Parser_atan2_Test, test4)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("atan2(V(A),V(B))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=-0.5, Bval=-0.25;
  double refRes = std::atan2(Aval,Bval);
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("B"),Bval);

  std::vector<double> derivs;

  double denom = Aval*Aval+Bval*Bval;
  std::vector<double> refderivs = { (Bval/denom), (-Aval/denom) };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(Double_Parser_atan2_Test, test4)
}

TEST ( Double_Parser_TwoNodeDeriv_Test, test1)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("V(A,B)"), testGroup);
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, vAval=5.0, vBval=1.0;

  solnGroup->setSoln(std::string("A"),vAval);
  solnGroup->setSoln(std::string("B"),vBval);

  double refRes = (vAval-vBval);

  std::vector<double> derivs;
  std::vector<double> refderivs = { 1, -1 };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(Double_Parser_TwoNodeDeriv_Test, test1)
}

//
TEST ( Double_Parser_poly_Test, test1)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("POLY(1) V(A) 1.0 2.0"), testGroup);
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=5.0;
  solnGroup->setSoln(std::string("A"),Aval);
  double refRes = 1.0 + 2.0*Aval;

  std::vector<double> derivs;
  std::vector<double> refderivs = { 2.0 };

  testExpression.evaluate(result, derivs);
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(Double_Parser_poly_Test, test1)
}


// poly test taken from amp_mod.cir
//
// BAM = 0 + 0*V(10) + 0*V(20) + 0*V(10)^2 + 1*V(10)*V(20) 
//
TEST ( Double_Parser_poly_Test, test2)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("POLY(2) V(A) V(B) 0 0 0 0 1 "), testGroup);
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, vAval=5.0, vBval=1.0;

  solnGroup->setSoln(std::string("A"),vAval);
  solnGroup->setSoln(std::string("B"),vBval);

  double refRes = vAval*vBval;

  std::vector<double> derivs;
  std::vector<double> refderivs = {  vBval, vAval };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(Double_Parser_poly_Test, test2)
}


// poly test taken from polyg.cir
//
// G1 8 9 POLY(2) 8 9 1 0  0 0 0 0 1
//
// this translates to:  POLY(2) V(8,9) V(1) 0 0 0 0 1
//
// And this means:
// G1 = V(8,9)*V(1)
//
TEST ( Double_Parser_poly_Test, test3)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("POLY(2) V(A,B) V(C) 0 0 0 0 1 "), testGroup);
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, vAval=5.0, vBval=1.0, vCval=2.0;

  solnGroup->setSoln(std::string("A"),vAval);
  solnGroup->setSoln(std::string("B"),vBval);
  solnGroup->setSoln(std::string("C"),vCval);

  double refRes = 
    (vAval-vBval)*vCval;

  std::vector<double> derivs;
  std::vector<double> refderivs = { vCval, -vCval, (vAval-vBval) };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(Double_Parser_poly_Test, test3)
}

TEST ( Double_Parser_poly_Test, test4)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("POLY(1) V(A,B) 0 1  "), testGroup);
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, vAval=5.0, vBval=1.0;

  solnGroup->setSoln(std::string("A"),vAval);
  solnGroup->setSoln(std::string("B"),vBval);

  double refRes = (vAval-vBval);

  std::vector<double> derivs;
  std::vector<double> refderivs = { 1, -1 };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(Double_Parser_poly_Test, test4)
}

// POLY(5) I(VB) I(VC)  I(VE)  I(VLP) I(VLN)  0 10.61E6 -10E6 10E6 10E6 -10E6
TEST ( Double_Parser_poly_Test, test5)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;

  // with parens
  Xyce::Util::newExpression testExpression(std::string("POLY(5) I(VB) I(VC)  I(VE)  I(VLP) I(VLN)  0 (10.61E6) (-10E6) (10E6) (10E6) (-10E6)"), testGroup);
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0; 
  double VB=-1.0;
  double VC=-2.0;
  double VE=-3.0;
  double VLP=-4.0;
  double VLN=-5.0;

  solnGroup->setSoln(std::string("VB"),VB);
  solnGroup->setSoln(std::string("VC"),VC);
  solnGroup->setSoln(std::string("VE"),VE);
  solnGroup->setSoln(std::string("VLP"),VLP);
  solnGroup->setSoln(std::string("VLN"),VLN);

  double refRes = 0 +  (10.61E6)*(-1) +  (-10E6)*(-2) +  (10E6)*(-3) +  (10E6)*(-4) +  (-10E6)*(-5);

  std::vector<double> derivs;
  //std::vector<double> refderivs = { 1, -1 };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(Double_Parser_poly_Test, test5)
}

// POLY(5) I(VB) I(VC)  I(VE)  I(VLP) I(VLN)  0 10.61E6 -10E6 10E6 10E6 -10E6
TEST ( Double_Parser_poly_Test, test6)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;

  // without parens
  Xyce::Util::newExpression testExpression(std::string("POLY(5) I(VB) I(VC)  I(VE)  I(VLP) I(VLN)  0 10.61E6 -10E6 10E6 10E6 -10E6"), testGroup);
  //
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0; 
  double VB=-1.0;
  double VC=-2.0;
  double VE=-3.0;
  double VLP=-4.0;
  double VLN=-5.0;

  solnGroup->setSoln(std::string("VB"),VB);
  solnGroup->setSoln(std::string("VC"),VC);
  solnGroup->setSoln(std::string("VE"),VE);
  solnGroup->setSoln(std::string("VLP"),VLP);
  solnGroup->setSoln(std::string("VLN"),VLN);

  double refRes = 0 +  (10.61E6)*(-1) +  (-10E6)*(-2) +  (10E6)*(-3) +  (10E6)*(-4) +  (-10E6)*(-5);

  std::vector<double> derivs;
  //std::vector<double> refderivs = { 1, -1 };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(Double_Parser_poly_Test, test6)
}

// POLY(5) I(VB) I(VC)  I(VE)  I(VLP) I(VLN)  0 10.61E6 -10E6 10E6 10E6 -10E6
TEST ( Double_Parser_poly_Test, test7)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;

  // without parens
  Xyce::Util::newExpression testExpression(std::string("POLY(5) I(VB) I(VC)  I(VE)  I(VLP) I(VLN)  +0 +10.61E6 -10E6 +10E6 +10E6 -10E6"), testGroup);
  //
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0; 
  double VB=-1.0;
  double VC=-2.0;
  double VE=-3.0;
  double VLP=-4.0;
  double VLN=-5.0;

  solnGroup->setSoln(std::string("VB"),VB);
  solnGroup->setSoln(std::string("VC"),VC);
  solnGroup->setSoln(std::string("VE"),VE);
  solnGroup->setSoln(std::string("VLP"),VLP);
  solnGroup->setSoln(std::string("VLN"),VLN);

  double refRes = 0 +  (10.61E6)*(-1) +  (-10E6)*(-2) +  (10E6)*(-3) +  (10E6)*(-4) +  (-10E6)*(-5);

  std::vector<double> derivs;
  //std::vector<double> refderivs = { 1, -1 };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);

  OUTPUT_MACRO(Double_Parser_poly_Test, test7)
}

TEST ( Double_Parser_poly_Test, test8)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("POLY(1) V(A) {1.0} {2.0}"), testGroup);
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=5.0;
  solnGroup->setSoln(std::string("A"),Aval);
  double refRes = 1.0 + 2.0*Aval;

  std::vector<double> derivs;
  std::vector<double> refderivs = { 2.0 };

  testExpression.evaluate(result, derivs);
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);

  OUTPUT_MACRO(Double_Parser_poly_Test, test8)
}

TEST ( Double_Parser_poly_Test, test9)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("POLY(1) V(A) ((1.0)) ((2.0))"), testGroup);
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=5.0;
  solnGroup->setSoln(std::string("A"),Aval);
  double refRes = 1.0 + 2.0*Aval;

  std::vector<double> derivs;
  std::vector<double> refderivs = { 2.0 };

  testExpression.evaluate(result, derivs);
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);

  OUTPUT_MACRO(Double_Parser_poly_Test, test9)
}

#if 0
TEST ( Double_Parser_poly_Test, test10)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  //Xyce::Util::newExpression testExpression(std::string("POLY(1) V(A) 1.0 2.0"), testGroup);
  Xyce::Util::newExpression testExpression(std::string("POLY(1) V(A) p1 p2"), testGroup);
  testExpression.lexAndParseExpression();


  Teuchos::RCP<Xyce::Util::newExpression> p1Expression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("1.0"), testGroup));
  Teuchos::RCP<Xyce::Util::newExpression> p2Expression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("2.0"), testGroup));
  p1Expression->lexAndParseExpression();
  p2Expression->lexAndParseExpression();
  std::string p1Name = "p1";
  std::string p2Name = "p2";
  solnGroup->addParam(p1Name,p1Expression);
  solnGroup->addParam(p2Name,p2Expression);

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result=0.0, Aval=5.0;
  solnGroup->setSoln(std::string("A"),Aval);
  double refRes = 1.0 + 2.0*Aval;

  std::vector<double> derivs;
  std::vector<double> refderivs = { 2.0 };

  testExpression.evaluate(result, derivs);
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
}
#endif


#if 0
TEST ( Double_Parser_NestedFunc_Test, func_cir)
{
  Teuchos::RCP<solnAndFuncExpressionGroup> funcGroup = Teuchos::rcp(new solnAndFuncExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression will use the .func f.
  //Xyce::Util::newExpression testExpression(std::string("f(v(x), v(y), v(z))"), testGroup);
  Xyce::Util::newExpression testExpression(std::string("f(4.0, 2.0, 3.0)"), testGroup);
  testExpression.lexAndParseExpression();

  // need to set up 3 functions:  fy, diff and f.
  //
  // this expression is the RHS of .func fy(x,y) {y}
  Teuchos::RCP<Xyce::Util::newExpression> fyExpression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("y"), testGroup) );
  Xyce::Util::newExpression fy_LHS (std::string("fy(x,y)"), testGroup);
  fy_LHS.lexAndParseExpression();
  std::vector<std::string> fyArgStrings ;
  fy_LHS.getFuncPrototypeArgStrings(fyArgStrings);
  fyExpression->setFunctionArgStringVec (fyArgStrings);
  fyExpression->lexAndParseExpression();
  std::string fyName; 
  fy_LHS.getFuncPrototypeName(fyName);

  // this expression is the RHS of .func diff(x,y) {x-y}
  Teuchos::RCP<Xyce::Util::newExpression> diffExpression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("x-y"), testGroup) );
  Xyce::Util::newExpression diff_LHS (std::string("diff(x,y)"), testGroup);
  diff_LHS.lexAndParseExpression();
  std::vector<std::string> diffArgStrings ;
  diff_LHS.getFuncPrototypeArgStrings(diffArgStrings);
  diffExpression->setFunctionArgStringVec (diffArgStrings);
  diffExpression->lexAndParseExpression();
  std::string diffName; 
  diff_LHS.getFuncPrototypeName(diffName);

  // this expression is the RHS of .func f(x,y,z) {diff(y,z)-4+fy(z,x)**2}
  Teuchos::RCP<Xyce::Util::newExpression> fExpression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("diff(y,z)-4+fy(z,x)**2"), testGroup) );
  Xyce::Util::newExpression f_LHS (std::string("f(x,y,z)"), testGroup);
  f_LHS.lexAndParseExpression();
  std::vector<std::string> fArgStrings ;
  f_LHS.getFuncPrototypeArgStrings(fArgStrings);
  fExpression->setFunctionArgStringVec (fArgStrings);
  fExpression->lexAndParseExpression();
  std::string fName; 
  f_LHS.getFuncPrototypeName(fName);

  funcGroup->addFunction( fyName ,  fyExpression);
  funcGroup->addFunction( diffName ,  diffExpression);
  funcGroup->addFunction( fName ,  fExpression);

  // what happens if not resolved?  IT WILL PRODUCE THE WRONG ANSWER
  // Putting these resolutions in "random" order.  If the order matters, then all this breaks.
  testExpression.resolveExpression();
  diffExpression->resolveExpression();
  fExpression->resolveExpression();
  fyExpression->resolveExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  double result;
  // x,y,z = 4,2,3
  // refresult = diff(y,z)-4+fy(z,x)**2;
  // refresult = diff(2,3)-4+fy(3,4)**2;
  // refresult = (2-3)-4+4*4;
  // refresult = -5+16;
  // refresult = 11;
  double refresult = 11;
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refresult );
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refresult );
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refresult );

  OUTPUT_MACRO(Double_Parser_NestedFunc_Test, func_cir)
}
#endif

//-------------------------------------------------------------------------------
// below is another version of the nested function test.
// In this version, the group is not relied upon for resolving external functions.
//
// Rather, the calling code takes care of it.
//
// I have been conflicted about whether or not the group should handle finding external dependencies such as functions. 
//
// this helper function is part of the test.
//-------------------------------------------------------------------------------
void resolveFunctions (
  std::unordered_map <std::string, Teuchos::RCP<Xyce::Util::newExpression> >  & functions_,
    Teuchos::RCP<Xyce::Util::newExpression> & expToResolve
    )
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & funcOpVec = expToResolve->getFuncOpVec ();
  for (int ii=0;ii<funcOpVec.size();++ii)
  {
    Teuchos::RCP<Xyce::Util::newExpression> exp;

    std::string lowerName = funcOpVec[ii]->getName();
    Xyce::Util::toLower(lowerName);

    if (functions_.find(lowerName) != functions_.end()) 
    {
      exp = functions_[lowerName];

      funcOpVec[ii]->setNode(exp->getAst());

      Teuchos::RCP<funcOp<usedType> > tmpPtr = Teuchos::rcp_dynamic_cast<funcOp<usedType> > (funcOpVec[ii]);
      if ( !(Teuchos::is_null(tmpPtr)) )
      {
        tmpPtr->setFuncArgs(  exp->getFunctionArgOpVec() );
      }
    }
  }
}

//-------------------------------------------------------------------------------
TEST ( Double_Parser_NestedFunc_Test, func_cir_newResolution)
{
  Teuchos::RCP<solnAndFuncExpressionGroup> funcGroup = Teuchos::rcp(new solnAndFuncExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression will use the .func f.
  //Xyce::Util::newExpression testExpression(std::string("f(v(x), v(y), v(z))"), testGroup);
  Teuchos::RCP<Xyce::Util::newExpression> testExpression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("f(4.0, 2.0, 3.0)"), testGroup));
  testExpression->lexAndParseExpression();

  std::unordered_map <std::string, Teuchos::RCP<Xyce::Util::newExpression> >  functions_;
  // need to set up 3 functions:  fy, diff and f, and put them in the map
  //
  // this expression is the RHS of .func fy(x,y) {y}
  Teuchos::RCP<Xyce::Util::newExpression> fyExpression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("y"), testGroup) );
  Xyce::Util::newExpression fy_LHS (std::string("fy(x,y)"), testGroup);
  fy_LHS.lexAndParseExpression();
  std::vector<std::string> fyArgStrings ;
  fy_LHS.getFuncPrototypeArgStrings(fyArgStrings);
  fyExpression->setFunctionArgStringVec (fyArgStrings);
  fyExpression->lexAndParseExpression();
  std::string fyName; 
  fy_LHS.getFuncPrototypeName(fyName);
  std::string lowerName = fyName;
  Xyce::Util::toLower(lowerName);
  functions_[lowerName] = fyExpression;

  // this expression is the RHS of .func diff(x,y) {x-y}
  Teuchos::RCP<Xyce::Util::newExpression> diffExpression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("x-y"), testGroup) );
  Xyce::Util::newExpression diff_LHS (std::string("diff(x,y)"), testGroup);
  diff_LHS.lexAndParseExpression();
  std::vector<std::string> diffArgStrings ;
  diff_LHS.getFuncPrototypeArgStrings(diffArgStrings);
  diffExpression->setFunctionArgStringVec (diffArgStrings);
  diffExpression->lexAndParseExpression();
  std::string diffName; 
  diff_LHS.getFuncPrototypeName(diffName);
  lowerName = diffName;
  Xyce::Util::toLower(lowerName);
  functions_[lowerName] = diffExpression;

  // this expression is the RHS of .func f(x,y,z) {diff(y,z)-4+fy(z,x)**2}
  Teuchos::RCP<Xyce::Util::newExpression> fExpression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("diff(y,z)-4+fy(z,x)**2"), testGroup) );
  Xyce::Util::newExpression f_LHS (std::string("f(x,y,z)"), testGroup);
  f_LHS.lexAndParseExpression();
  std::vector<std::string> fArgStrings ;
  f_LHS.getFuncPrototypeArgStrings(fArgStrings);
  fExpression->setFunctionArgStringVec (fArgStrings);
  fExpression->lexAndParseExpression();
  std::string fName; 
  f_LHS.getFuncPrototypeName(fName);
  lowerName = fName;
  Xyce::Util::toLower(lowerName);
  functions_[lowerName] = fExpression;


// a resolution must happen here, but I don't want to use the "resolveExpression" function.
// I want to use a stand-alone function instead, where the (a) the "group" is not involved 
// and (b) the function lookup is not the group, or the new-expression classes responsibility. 
#if 0 
  funcGroup->addFunction( fyName ,  fyExpression);
  funcGroup->addFunction( diffName ,  diffExpression);
  funcGroup->addFunction( fName ,  fExpression);
  testExpression.resolveExpression();
  diffExpression->resolveExpression();
  fExpression->resolveExpression();
  fyExpression->resolveExpression();
#else
  resolveFunctions(functions_, diffExpression);
  resolveFunctions(functions_, fExpression);
  resolveFunctions(functions_, fyExpression);
  resolveFunctions(functions_, testExpression);
#endif

  double result;
  // x,y,z = 4,2,3
  // refresult = diff(y,z)-4+fy(z,x)**2;
  // refresult = diff(2,3)-4+fy(3,4)**2;
  // refresult = (2-3)-4+4*4;
  // refresult = -5+16;
  // refresult = 11;
  double refresult = 11;
  testExpression->evaluateFunction(result);   EXPECT_EQ( result, refresult );

  OUTPUT_MACRO3(Double_Parser_NestedFunc_Test, func_cir)
}


//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
void newResolveFunctions (
  std::unordered_map <std::string, Teuchos::RCP<Xyce::Util::newExpression> >  & functions_,
    Teuchos::RCP<Xyce::Util::newExpression> & expToResolve
    )
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & funcOpVec = expToResolve->getFuncOpVec ();
  for (int ii=0;ii<funcOpVec.size();++ii)
  {
    Teuchos::RCP<Xyce::Util::newExpression> exp;

    std::string upperName = funcOpVec[ii]->getName();
    Xyce::Util::toUpper(upperName);

    if (functions_.find(upperName) != functions_.end()) 
    {
      exp = functions_[upperName];
      expToResolve->attachFunctionNode(upperName, exp);
    }
  }
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
TEST ( Double_Parser_NestedFunc_Test, func_cir_newResolution2)
{
  Teuchos::RCP<solnAndFuncExpressionGroup> funcGroup = Teuchos::rcp(new solnAndFuncExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression will use the .func f.
  //Xyce::Util::newExpression testExpression(std::string("f(v(x), v(y), v(z))"), testGroup);
  Teuchos::RCP<Xyce::Util::newExpression> testExpression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("f(4.0, 2.0, 3.0)"), testGroup));
  testExpression->lexAndParseExpression();

  std::unordered_map <std::string, Teuchos::RCP<Xyce::Util::newExpression> >  functions_;
  // need to set up 3 functions:  fy, diff and f, and put them in the map
  //
  // this expression is the RHS of .func fy(x,y) {y}
  Teuchos::RCP<Xyce::Util::newExpression> fyExpression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("y"), testGroup) );
  Xyce::Util::newExpression fy_LHS (std::string("fy(x,y)"), testGroup);
  fy_LHS.lexAndParseExpression();
  std::vector<std::string> fyArgStrings ;
  fy_LHS.getFuncPrototypeArgStrings(fyArgStrings);
  fyExpression->setFunctionArgStringVec (fyArgStrings);
  fyExpression->lexAndParseExpression();
  std::string fyName; 
  fy_LHS.getFuncPrototypeName(fyName);
  std::string upperName = fyName;
  Xyce::Util::toUpper(upperName);
  functions_[upperName] = fyExpression;

  // this expression is the RHS of .func diff(x,y) {x-y}
  Teuchos::RCP<Xyce::Util::newExpression> diffExpression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("x-y"), testGroup) );
  Xyce::Util::newExpression diff_LHS (std::string("diff(x,y)"), testGroup);
  diff_LHS.lexAndParseExpression();
  std::vector<std::string> diffArgStrings ;
  diff_LHS.getFuncPrototypeArgStrings(diffArgStrings);
  diffExpression->setFunctionArgStringVec (diffArgStrings);
  diffExpression->lexAndParseExpression();
  std::string diffName; 
  diff_LHS.getFuncPrototypeName(diffName);
  upperName = diffName;
  Xyce::Util::toUpper(upperName);
  functions_[upperName] = diffExpression;

  // this expression is the RHS of .func f(x,y,z) {diff(y,z)-4+fy(z,x)**2}
  Teuchos::RCP<Xyce::Util::newExpression> fExpression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("diff(y,z)-4+fy(z,x)**2"), testGroup) );
  Xyce::Util::newExpression f_LHS (std::string("f(x,y,z)"), testGroup);
  f_LHS.lexAndParseExpression();
  std::vector<std::string> fArgStrings ;
  f_LHS.getFuncPrototypeArgStrings(fArgStrings);
  fExpression->setFunctionArgStringVec (fArgStrings);
  fExpression->lexAndParseExpression();
  std::string fName; 
  f_LHS.getFuncPrototypeName(fName);
  upperName = fName;
  Xyce::Util::toUpper(upperName);
  functions_[upperName] = fExpression;


// a resolution must happen here, but I don't want to use the "resolveExpression" function.
// I want to use a stand-alone function instead, where the (a) the "group" is not involved 
// and (b) the function lookup is not the group, or the new-expression classes responsibility. 
  newResolveFunctions(functions_, diffExpression);
  newResolveFunctions(functions_, fExpression);
  newResolveFunctions(functions_, fyExpression);
  newResolveFunctions(functions_, testExpression);

  double result;
  // x,y,z = 4,2,3
  // refresult = diff(y,z)-4+fy(z,x)**2;
  // refresult = diff(2,3)-4+fy(3,4)**2;
  // refresult = (2-3)-4+4*4;
  // refresult = -5+16;
  // refresult = 11;
  double refresult = 11;
  testExpression->evaluateFunction(result);   EXPECT_EQ( result, refresult );

  OUTPUT_MACRO3(Double_Parser_NestedFunc_Test, func_cir)
}

//-------------------------------------------------------------------------------
// The point of these next tests is to test out large numbers of 
// nested function calls, or global_params
//
// This is the simplest list of function calls you could imagine, as
// all but one are pass-thrus.
//
// As the number of funcs increases, the size of the AST traversal increases.
//
// evaluateFunction which excludes derivatives scales much better than 
// evaluate, which includes them.  The funcOp::dx function needs optimization.
//
// Also, it would probably be good if we could automatically prune or 
// shrink the tree.
//-------------------------------------------------------------------------------
TEST ( Double_Parser_NestedFunc_Test, 1000nest_no_deriv)
{
  Teuchos::RCP<solnAndFuncExpressionGroup> funcGroup = Teuchos::rcp(new solnAndFuncExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  int numFuncs=1000;

  // this expression will use the .func f.
  std::string testExprFunctionString = std::string("f") + std::to_string(numFuncs-1) + std::string("(V(A))");
  Teuchos::RCP<Xyce::Util::newExpression> testExpression  = Teuchos::rcp(new Xyce::Util::newExpression( testExprFunctionString , testGroup));
  testExpression->lexAndParseExpression();


  // need to set up vector of functions  fn, all of the form:
  // this expression is the RHS of .func fn(x) {x}
  std::vector< std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > >  dotFuncVector;
  // set up the first func in the list (which will be the final one evaluated)  f0(x) {2*x}
  
  std::string fnName = "F" + std::to_string(0);
  std::vector<std::string> fnArgStrings = {"X"} ;

  Teuchos::RCP<Xyce::Util::newExpression> fnExpression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2*x"), testGroup) );
  fnExpression->setFunctionArgStringVec (fnArgStrings);
  fnExpression->lexAndParseExpression();
  std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > fnPair(fnName,fnExpression);
  dotFuncVector.push_back(fnPair);

  // need to set up vector of functions  fn, all of the form:
  // this expression is the RHS of .func fn(x) {x}
  for (int ii=1;ii<numFuncs;ii++)
  {
    // set up the func
    std::string fnName = "F" + std::to_string(ii);
    std::vector<std::string> fnArgStrings = {"X"} ;

    std::string rhsFname = std::string("f") + std::to_string(ii-1) + std::string("(x)");
    Teuchos::RCP<Xyce::Util::newExpression> fnExpression  = Teuchos::rcp(new Xyce::Util::newExpression(rhsFname, testGroup));
    fnExpression->setFunctionArgStringVec (fnArgStrings);
    fnExpression->lexAndParseExpression();

    std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > fnPair(fnName,fnExpression);
    dotFuncVector.push_back(fnPair);
  }

  // do all the attachments
  for (int ii=1;ii<numFuncs;ii++)
  {
    std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > fnPairM1 = dotFuncVector[ii-1];
    std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > fnPair   = dotFuncVector[ii];

    fnPair.second->attachFunctionNode(fnPairM1.first, fnPairM1.second);
  }

  testExpression->attachFunctionNode(dotFuncVector[numFuncs-1].first, dotFuncVector[numFuncs-1].second);

  double result=0.0, Aval=5.0;
  funcGroup->setSoln(std::string("A"),Aval);

  std::vector<double> derivs;
  std::vector<double> refderivs = {2.0};
  double refresult = 10.0;
  testExpression->evaluateFunction(result);   EXPECT_EQ( result, refresult );

  OUTPUT_MACRO3 ( Double_Parser_NestedFunc_Test, 1000nest_no_deriv)
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
TEST ( Double_Parser_NestedFunc_Test, 200nest_with_deriv)
{
  Teuchos::RCP<solnAndFuncExpressionGroup> funcGroup = Teuchos::rcp(new solnAndFuncExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  int numFuncs=200;

  // this expression will use the .func f.
  std::string testExprFunctionString = std::string("f") + std::to_string(numFuncs-1) + std::string("(V(A))");
  Teuchos::RCP<Xyce::Util::newExpression> testExpression  = Teuchos::rcp(new Xyce::Util::newExpression( testExprFunctionString, testGroup));
  testExpression->lexAndParseExpression();


  // need to set up vector of functions  fn, all of the form:
  // this expression is the RHS of .func fn(x) {x}
  std::vector< std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > >  dotFuncVector;
  // set up the first func in the list (which will be the final one evaluated)  f0(x) {2*x}
  
  std::string fnName = "F" + std::to_string(0);
  std::vector<std::string> fnArgStrings = {"X"} ;

  Teuchos::RCP<Xyce::Util::newExpression> fnExpression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2*x"), testGroup) );
  fnExpression->setFunctionArgStringVec (fnArgStrings);
  fnExpression->lexAndParseExpression();
  std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > fnPair(fnName,fnExpression);
  dotFuncVector.push_back(fnPair);

  // need to set up vector of functions  fn, all of the form:
  // this expression is the RHS of .func fn(x) {x}
  for (int ii=1;ii<numFuncs;ii++)
  {
    // set up the func
    std::string fnName = "F" + std::to_string(ii);
    std::vector<std::string> fnArgStrings = {"X"} ;

    std::string rhsFname = std::string("f") + std::to_string(ii-1) + std::string("(x)");
    Teuchos::RCP<Xyce::Util::newExpression> fnExpression  = Teuchos::rcp(new Xyce::Util::newExpression(rhsFname, testGroup));
    fnExpression->setFunctionArgStringVec (fnArgStrings);
    fnExpression->lexAndParseExpression();

    std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > fnPair(fnName,fnExpression);
    dotFuncVector.push_back(fnPair);
  }

  // do all the attachments
  for (int ii=1;ii<numFuncs;ii++)
  {
    std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > fnPairM1 = dotFuncVector[ii-1];
    std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > fnPair   = dotFuncVector[ii];

    fnPair.second->attachFunctionNode(fnPairM1.first, fnPairM1.second);
  }

  testExpression->attachFunctionNode(dotFuncVector[numFuncs-1].first, dotFuncVector[numFuncs-1].second);

  double result=0.0, Aval=5.0;
  funcGroup->setSoln(std::string("A"),Aval);

  std::vector<double> derivs;
  std::vector<double> refderivs = {2.0};
  double refresult = 10.0;

  // the dx function (called under evaluate) has some bottlenecks
  testExpression->evaluate(result,derivs);   
  EXPECT_EQ( result, refresult );
  EXPECT_EQ( derivs, refderivs);

  OUTPUT_MACRO3(Double_Parser_NestedFunc_Test, 200nest_with_deriv)
}

//-------------------------------------------------------------------------------
TEST ( Double_Parser_NestedGlobalParam_Test, 1000nest_no_deriv)
{
  Teuchos::RCP<solnAndFuncExpressionGroup> funcGroup = Teuchos::rcp(new solnAndFuncExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  int numParams=1000;

  // this expression will use the parameter p1.
  std::string testExprFunctionString = std::string("p") + std::to_string(numParams-1) + std::string("*V(A)");
  Teuchos::RCP<Xyce::Util::newExpression> testExpression  = Teuchos::rcp(new Xyce::Util::newExpression( testExprFunctionString , testGroup));
  testExpression->lexAndParseExpression();

  // need to set up vector of global parameter
  // this expression is the RHS of .global_param p0 {2.0}
  std::vector< std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > >  dotGlobalParamVector;
  // set up the first param in the list (which will be the final one evaluated)  
  
  std::string parName = "P" + std::to_string(0);

  Teuchos::RCP<Xyce::Util::newExpression> parExpression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2.0"), testGroup) );
  parExpression->lexAndParseExpression();
  std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > parPair(parName,parExpression);
  dotGlobalParamVector.push_back(parPair);

  // need to set up vector of parameters:
  for (int ii=1;ii<numParams;ii++)
  {
    // set up the func
    std::string parName = "P" + std::to_string(ii);

    std::string rhsFname = std::string("p") + std::to_string(ii-1);
    Teuchos::RCP<Xyce::Util::newExpression> parExpression  = Teuchos::rcp(new Xyce::Util::newExpression(rhsFname, testGroup));
    parExpression->lexAndParseExpression();

    std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > parPair(parName,parExpression);
    dotGlobalParamVector.push_back(parPair);
  }

  // do all the attachments
  for (int ii=1;ii<numParams;ii++)
  {
    std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > parPairM1 = dotGlobalParamVector[ii-1];
    std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > parPair   = dotGlobalParamVector[ii];

    parPair.second->attachParameterNode(parPairM1.first, parPairM1.second);
  }

  testExpression->attachParameterNode(dotGlobalParamVector[numParams-1].first, dotGlobalParamVector[numParams-1].second);

  double result=0.0, Aval=5.0;
  funcGroup->setSoln(std::string("A"),Aval);

  std::vector<double> derivs;
  std::vector<double> refderivs = {2.0};
  double refresult = 10.0;
  testExpression->evaluateFunction(result);   EXPECT_EQ( result, refresult );

  OUTPUT_MACRO3(Double_Parser_NestedGlobalParam_Test, 1000nest_no_deriv)
}

int main (int argc, char **argv)
{
  {
  char filename[ ] = "parserUnitTest.out";
  std::fstream outputFile;
  outputFile.open(filename,  std::fstream::in | std::fstream::out | std::fstream::trunc);
  outputFile << "parserUnitTest output file.  This is to test the output functions" <<std::endl;
  outputFile.close();
  }

  {
  char filename[ ] = "parserUnitTest_codeGen.C";
  std::fstream outputFile;
  outputFile.open(filename,  std::fstream::in | std::fstream::out | std::fstream::trunc);
  outputFile << "// parserUnitTest code gen file.  This is to test the codegen functions" <<std::endl;
  outputFile.close();
  }

  testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}


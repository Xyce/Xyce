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

//-------------------------------------------------------------------------
//
// Purpose        : Parser level unit tests for the expression library,
//                  std::complex<double> data type version.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/xx/2019
//
//
//
//
//-------------------------------------------------------------------------
//
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
#include <random>

#include "ast.h"
#include <newExpression.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_ExtendedString.h>
#include <N_IO_fwd.h>

//-------------------------------------------------------------------------------
// group classes.  To use the new expression library, it must be passed a group
// object, usually at construction of the expression object.  This group object
// must be derived from the base expression group class.  To be useful, they
// must implement some of the virtual functions of that class.
//
// The groups in this file were originally written "as needed" and thus are
// partially (or entirely) redundant with each other.  They were originally
// scattered throughout this file, near the various tests that used them, 
// but the file has been rearranged to have them at the top.  These groups
// could probably be consolidated to have a smaller list.
//-------------------------------------------------------------------------------
class testExpressionGroup : public Xyce::Util::baseExpressionGroup
{
  public:
    testExpressionGroup () : Xyce::Util::baseExpressionGroup()  {};
    ~testExpressionGroup () {};
};

//-------------------------------------------------------------------------------
class solutionGroup : public Xyce::Util::baseExpressionGroup
{
  public:
    solutionGroup () : Xyce::Util::baseExpressionGroup(), time(0.0), freq(0.0), gmin(0.0) , temp(0.0), VT(0.0), timeStep(0.0), stepNumber(0)
  {};
    ~solutionGroup () {};

  virtual double getTime() { return time; };
  void setTime(double t) { time = t; };

  virtual double getFreq() { return freq; };
  void setFreq(double f) { freq = f; };

  virtual double getGmin() { return gmin; };
  void setGmin(double g) { gmin = g; };

  virtual double getTemp() { return temp; };
  void setTemp(double t) { temp = t; };

  virtual double getVT() { return VT; };
  void setVT(double t) { VT = t; };

  virtual double getTimeStep() { return timeStep; };
  void setTimeStep(double dt) { timeStep = dt; };

  void setStepNumber (unsigned int number) { stepNumber = number; }
  unsigned int getStepNumber () { return stepNumber; }

  void setSoln(const std::string & name, std::complex<double> val)
  {
    std::string lowerName = name;
    Xyce::Util::toLower(lowerName);
    internalVars_[lowerName] = val;
  };

  bool getCurrentVal(const std::string & deviceName, const std::string & designator, std::complex<double> & retval )
  {
    return getSolutionVal(deviceName, retval);
  }

  virtual bool getSolutionVal(const std::string & name, std::complex<double> & val )
  {
    bool retval=true;
    std::string lowerName = name;
    Xyce::Util::toLower(lowerName);
    if (internalVars_.find(lowerName) != internalVars_.end()) { val = internalVars_[lowerName]; }
    else { retval = false; val=-2.0e9; }
    return retval;
  }

  void setPower(const std::string & deviceName, std::complex<double> val) { return setSoln(deviceName, val); }
  bool getPower(const std::string & tag, const std::string & deviceName, std::complex<double> & retval) { return getSolutionVal(deviceName, retval); }

  private:
    double time;
    double freq;
    double gmin;
    double temp;
    double VT;
    double timeStep;
    unsigned int stepNumber;
    std::unordered_map <std::string, std::complex<double> > internalVars_;
};

// these are here to support some deprecated groups.  They were redundant 
// with solutionGroup, so I got rid of them, but didn't want to edit every single 
// test to update.  Most of these groups were written on the fly to support a 
// small handful of tests, but ultimately they all basically did the same 
// sorts of things.  So they are gone, replaced by typedefs.
typedef  solutionGroup solnExpressionGroup;
typedef  solutionGroup timeDepExpressionGroup;
typedef  solutionGroup currSolnExpressionGroup;
typedef  solutionGroup testExpressionGroupWithFuncSupport;
typedef  solutionGroup ifStatementExpressionGroup;
typedef  solutionGroup tempDepExpressionGroup;
typedef  solutionGroup Bsrc_C1_ExpressionGroup;
typedef  solutionGroup solnAndFuncExpressionGroup;
typedef  solutionGroup solnExpressionGroup2;
typedef  solutionGroup testExpressionGroupWithParamSupport;
typedef  solutionGroup sdtExpressionGroup;

//-------------------------------------------------------------------------------
class leadCurrentExpressionGroup : public Xyce::Util::baseExpressionGroup
{
  public:
    leadCurrentExpressionGroup () :
      Xyce::Util::baseExpressionGroup() {};
    ~leadCurrentExpressionGroup () {};

  bool getCurrentVal  ( const std::string & deviceName, const std::string & designator, std::complex<double> & retval )
  {
    bool retbool = true;
    retval=0.0;
    std::string des=designator;
    Xyce::Util::toUpper(des);

    if (des==std::string("IG"))
    {
      if (IGvalues.find(deviceName) != IGvalues.end()) { retval = IGvalues[deviceName]; }
      else { retbool = false; }
    }
    else if (des==std::string("ID"))
    {
      if (IDvalues.find(deviceName) != IDvalues.end()) { retval = IDvalues[deviceName]; }
      else { retbool = false; }
    }
    else if (des==std::string("IS"))
    {
      if (ISvalues.find(deviceName) != ISvalues.end()) { retval = ISvalues[deviceName]; }
      else { retbool = false; }
    }

    return retbool;
  }

  bool setCurrentVal  ( const std::string & deviceName, const std::string & designator, std::complex<double> setval )
  {
    std::string des=designator;
    Xyce::Util::toUpper(des);

    if (des == std::string("IG")) IGvalues[deviceName] = setval;
    else if(des == std::string("ID")) IDvalues[deviceName] = setval;
    else if(des == std::string("IS")) ISvalues[deviceName] = setval;

    return true;
  }

  private:
    std::unordered_map <std::string, std::complex<double> > IGvalues;
    std::unordered_map <std::string, std::complex<double> > IDvalues;
    std::unordered_map <std::string, std::complex<double> > ISvalues;
};

//-------------------------------------------------------------------------------
class internalDevExpressionGroup : public Xyce::Util::baseExpressionGroup
{
  public:
    internalDevExpressionGroup () : Xyce::Util::baseExpressionGroup() {};
    ~internalDevExpressionGroup () {};

  void setInternalDeviceVar (const std::string & name, std::complex<double> val)
  {
    std::string lowerName = name;
    Xyce::Util::toLower(lowerName);
    internalVars_[lowerName] = val;
  };

  bool getInternalDeviceVar       (const std::string & name, std::complex<double> & val)
  {
    bool retval=true;
    std::string lowerName = name;
    Xyce::Util::toLower(lowerName);
    if (internalVars_.find(lowerName) != internalVars_.end()) { val = internalVars_[lowerName]; }
    else { retval = false; }
    return retval;
  }

  private:
    std::unordered_map <std::string, std::complex<double> > internalVars_;
};

//-------------------------------------------------------------------------------
class noiseExpressionGroup : public Xyce::Util::baseExpressionGroup
{
  public:
    noiseExpressionGroup () : Xyce::Util::baseExpressionGroup(), inoise_(std::complex<double>(0.0,0.0)), onoise_(std::complex<double>(0.0,0.0)) {};
    ~noiseExpressionGroup () {};

  void setDnoNoiseDeviceVar (const std::vector<std::string> & names, std::complex<double> val)
  {
    std::vector<std::string> lowerNames;
    if ( !(names.empty()) )
    {
      lowerNames = names;
      for (int ii=0;ii<lowerNames.size();ii++) { Xyce::Util::toLower(lowerNames[ii]); }
      dnoDeviceVars_[lowerNames[0]] = val;
    }
  };

  void setDniNoiseDeviceVar (const std::vector<std::string> & names, std::complex<double> val)
  {
    std::vector<std::string> lowerNames;
    if ( !(names.empty()) )
    {
      lowerNames = names;
      for (int ii=0;ii<lowerNames.size();ii++) { Xyce::Util::toLower(lowerNames[ii]); }
      dniDeviceVars_[lowerNames[0]] = val;
    }
  };

  void setONoise (std::complex<double> val) { onoise_ = val; };
  void setINoise (std::complex<double> val) { inoise_ = val; };

  virtual bool getDnoNoiseDeviceVar(const std::vector<std::string> & deviceNames, std::complex<double> & val) 
  { 
    bool retval=true;
    std::vector<std::string> lowerNames;
    if ( !(deviceNames.empty()) )
    {
      lowerNames = deviceNames;
      for (int ii=0;ii<lowerNames.size();ii++) { Xyce::Util::toLower(lowerNames[ii]); }
    }
    if (dnoDeviceVars_.find(lowerNames[0]) != dnoDeviceVars_.end()) { val = dnoDeviceVars_[lowerNames[0]]; }
    else { retval = false; }
    return retval;
  }

  virtual bool getDniNoiseDeviceVar(const std::vector<std::string> & deviceNames, std::complex<double> & val) 
  { 
    bool retval=true;
    std::vector<std::string> lowerNames;
    if ( !(deviceNames.empty()) )
    {
      lowerNames = deviceNames;
      for (int ii=0;ii<lowerNames.size();ii++) { Xyce::Util::toLower(lowerNames[ii]); }
    }
    if (dniDeviceVars_.find(lowerNames[0]) != dniDeviceVars_.end()) { val = dniDeviceVars_[lowerNames[0]]; }
    else { retval = false; }
    return retval;
  }

  virtual bool getONoise(std::complex<double> & retval) { retval=onoise_; return true; }
  virtual bool getINoise(std::complex<double> & retval) { retval=inoise_; return true; }

  private:
    std::unordered_map <std::string, std::complex<double> > dnoDeviceVars_;
    std::unordered_map <std::string, std::complex<double> > dniDeviceVars_;
    std::complex<double> inoise_, onoise_;
};

//-------------------------------------------------------------------------------
class rfParamGroup : public Xyce::Util::baseExpressionGroup
{
  public:
    rfParamGroup () : Xyce::Util::baseExpressionGroup() 
    { };

    ~rfParamGroup () {};

  public:
    void setSparam( const std::vector< std::vector< std::complex<double> > > & sp ) { sparams_ = sp; }
    void setYparam( const std::vector< std::vector< std::complex<double> > > & yp ) { yparams_ = yp; }
    void setZparam( const std::vector< std::vector< std::complex<double> > > & zp ) { zparams_ = zp; }

    virtual bool getSparam (const std::vector<int> & args, std::complex<double> & retval ) 
    { 
      int i1=args[0]-1;
      int i2=args[1]-1;
      retval=0.0;
      if (i1<sparams_.size()) { retval=sparams_[i1][i2]; }
      return true; 
    }

    virtual bool getYparam (const std::vector<int> & args, std::complex<double> & retval ) 
    { 
      int i1=args[0]-1;
      int i2=args[1]-1;
      retval=0.0;
      if (i1<yparams_.size()) { retval=yparams_[i1][i2]; }
      return true; 
    }

    virtual bool getZparam (const std::vector<int> & args, std::complex<double> & retval ) 
    { 
      int i1=args[0]-1;
      int i2=args[1]-1;
      retval=0.0;
      if (i1<zparams_.size()) { retval=zparams_[i1][i2]; }
      return true; 
    }

  std::vector< std::vector< std::complex<double> > > sparams_;
  std::vector< std::vector< std::complex<double> > > yparams_;
  std::vector< std::vector< std::complex<double> > > zparams_;
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
  std::complex<double> result(0.0,0.0); \
  std::complex<double> cppExp = (CPPEXP); \
  testExpression.evaluateFunction(result); \
  std::complex<double> testVal = result-cppExp; \
  EXPECT_NEAR( (std::real(testVal)), 0.0, 2.0e-14); \
  EXPECT_NEAR( (std::imag(testVal)), 0.0, 2.0e-14); \
  Xyce::Util::newExpression copyExpression(testExpression); \
  copyExpression.evaluateFunction(result); \
  testVal = result-cppExp; \
  EXPECT_NEAR( (std::real(testVal)), 0.0, 2.0e-14); \
  EXPECT_NEAR( (std::imag(testVal)), 0.0, 2.0e-14); \
  Xyce::Util::newExpression assignExpression; \
  assignExpression = testExpression; \
  assignExpression.evaluateFunction(result); \
  testVal = result-cppExp; \
  EXPECT_NEAR( (std::real(testVal)), 0.0, 2.0e-14); \
  EXPECT_NEAR( (std::imag(testVal)), 0.0, 2.0e-14); \
}

#define PARSER_SIMPLE_TEST_MACRO2(NAME,SUBNAME,STREXP, CPPEXP) \
TEST ( NAME, SUBNAME ) \
{ \
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() ); \
  Xyce::Util::newExpression testExpression(std::string(STREXP), testGroup); \
  testExpression.lexAndParseExpression(); \
  std::complex<double> result(0.0,0.0); \
  testExpression.evaluateFunction(result); \
  EXPECT_DOUBLE_EQ( std::real(result), std::real(CPPEXP)); \
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(CPPEXP)); \
  Xyce::Util::newExpression copyExpression(testExpression); \
  copyExpression.evaluateFunction(result); \
  EXPECT_DOUBLE_EQ( std::real(result), std::real(CPPEXP)); \
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(CPPEXP)); \
  Xyce::Util::newExpression assignExpression; \
  assignExpression = testExpression; \
  assignExpression.evaluateFunction(result); \
  EXPECT_DOUBLE_EQ( std::real(result), std::real(CPPEXP)); \
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(CPPEXP)); \
}

//-------------------------------------------------------------------------------
// Small function to make the creation of .func related objects simpler.
//
// This function mimics how .func creation and resolution happen inside of Xyce.
//
// The lhs and rhs strings are from the original .func declaration.  
// They are inputs.  For example:  .func f1(a) {sqrt(a)}
//     lhs = "f1(a)"
//     rhs = "sqrt(a)"
//
// The lhs string is used to create a temporary local expression 
// object (new_LHS).  This object isn't needed outside this function.  It is
// used only to pull out the name "f1" and the arguments "a".  The name "f1" 
// is one of the outputs of this function, and the arguments "a" are needed 
// to fully resolve the rhs expression.  The rhs expression is the other
// output of this function.
//
// The "newName" and "newExpression" are the outputs, and they contain (in 
// the above example) "f1" and "sqrt(a)" .  These are the only two 
// objects needed to attach this .func to another expression.
//-------------------------------------------------------------------------------
void createFunc( std::string & lhs, std::string & rhs,
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> & testGroup,
  std::string & newName,
  Teuchos::RCP<Xyce::Util::newExpression> & newExpression
    )
{
  newExpression = Teuchos::rcp(new Xyce::Util::newExpression( rhs , testGroup));
  Xyce::Util::newExpression new_LHS (lhs, testGroup);
  new_LHS.lexAndParseExpression();

  std::vector<std::string> newArgStrings ;
  new_LHS.getFuncPrototypeArgStrings(newArgStrings);
  newExpression->setFunctionArgStringVec (newArgStrings);

  newExpression->lexAndParseExpression();

  new_LHS.getFuncPrototypeName(newName);
}

//-------------------------------------------------------------------------------
// number by itself
TEST ( ComplexParserTest, numval)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("1.0+2.0J"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0,0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result,std::complex<double>(1.0,2.0) );

  Xyce::Util::newExpression copyExpression(testExpression); 
  copyExpression.evaluateFunction(result); 
  EXPECT_EQ( result,std::complex<double>(1.0,2.0) );

  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 
  assignExpression.evaluateFunction(result); 
  EXPECT_EQ( result,std::complex<double>(1.0,2.0) );
}


TEST ( ComplexParserTest, numval2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("10.61E6"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result (10.61E6);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( (result-(10.61E6)), 0.0);

  Xyce::Util::newExpression copyExpression(testExpression); 
  copyExpression.evaluateFunction(result); 
  EXPECT_EQ( (result-(10.61E6)), 0.0);

  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 
  assignExpression.evaluateFunction(result); 
  EXPECT_EQ( (result-(10.61E6)), 0.0);
}

// these next 3 tests are for parameters that happen to have the same name as single-character operators
TEST ( ComplexParserTest, singleParam_R)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("R"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, 0.0);
}

TEST ( ComplexParserTest, singleParam_M)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("M"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, 0.0);
}

TEST ( ComplexParserTest, singleParam_P)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("P"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, 0.0);
}

TEST ( ComplexParserTest, singleParam_E)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("E"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, 0.0);
}

TEST ( ComplexParserTest, singleParam_J)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("J"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, 0.0);
}

// this test is to ensure that a parameter name can include a period "."
//INVALIDLINES.S2P
TEST ( ComplexParserTest, param_INVALIDLINES_dot_S2PJ)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("INVALIDLINES.S2P"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, 0.0);
}

// this test is to ensure that a parameter name can include a period "!"
// See issue 191 on gitlab-ex.
TEST ( ComplexParserTest, param_exclamationPoint)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("YCAP!CAP2F:R"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, 0.0);
}

// these next 3 tests are for the single-character operators, which are mostly relevant to complex numbers
// In some codes, R(number) means "real part" of number.
TEST ( ComplexParserTest, singleCharacter_Rop)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("R(1.0+2.0J)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, 1.0);
}

// make sure it also works with a space, as the lexing of R operator requires that the 
// left paren be part of the token, so whitespace is optionally handled via regex
TEST ( ComplexParserTest, singleCharacter_Rop2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("R (1.0+2.0J)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, 1.0);
}


// M(number) = abs of complex number.
TEST ( ComplexParserTest, singleCharacter_Mop)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("M(1.0+2.0J)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> input(1.0,2.0);
  std::complex<double> refresult(std::abs(input));
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, refresult);
}



// M(parameter) = abs of complex valued parameter
TEST ( ComplexParserTest, singleCharacter_Mop_Param)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("M(par1)"), testGroup);
  testExpression.lexAndParseExpression();

  // setup the parameter
  Teuchos::RCP<Xyce::Util::newExpression> par1Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("1.0+2.0J"), testGroup));
  par1Expression->lexAndParseExpression();
  std::string par1Name = "par1";
  testExpression.attachParameterNode(par1Name,par1Expression);

  std::complex<double> input(1.0,2.0);
  std::complex<double> refresult(std::abs(input));
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, refresult);
}

// M(parameter) = abs of complex valued expression
TEST ( ComplexParserTest, singleCharacter_Mop_Expr)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("M(par1 + par2*4.0)"), testGroup);
  testExpression.lexAndParseExpression();

  // setup parameter par1
  Teuchos::RCP<Xyce::Util::newExpression> par1Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("1.0+2.0J"), testGroup));
  par1Expression->lexAndParseExpression();
  std::string par1Name = "par1";
  testExpression.attachParameterNode(par1Name,par1Expression);

  // setup parameter par2
  Teuchos::RCP<Xyce::Util::newExpression> par2Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("1.0+2.0J"), testGroup));
  par2Expression->lexAndParseExpression();
  std::string par2Name = "par2";
  testExpression.attachParameterNode(par2Name,par2Expression);

  std::complex<double> input(1.0,2.0);
  std::complex<double> refresult(std::abs(input)*5.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, refresult);
}

// M(parameter) = abs of complex valued expression with a function
TEST ( ComplexParserTest, singleCharacter_Mop_Func)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("M(par1 + F1(4.0))"), testGroup);
  testExpression.lexAndParseExpression();

  // setup parameter par1
  Teuchos::RCP<Xyce::Util::newExpression> par1Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("1.0+2.0J"), testGroup));
  par1Expression->lexAndParseExpression();
  std::string par1Name = "par1";
  testExpression.attachParameterNode(par1Name,par1Expression);

  // setup function F1
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("A*(1.0+2.0J)"), testGroup));

  Xyce::Util::newExpression f1_LHS (std::string("F1(A)"), testGroup);
  f1_LHS.lexAndParseExpression();

  std::vector<std::string> f1ArgStrings ;
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);
  f1Expression->setFunctionArgStringVec (f1ArgStrings);
  f1Expression->lexAndParseExpression();
  std::string f1Name;
  f1_LHS.getFuncPrototypeName(f1Name);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  std::complex<double> input(1.0,2.0);
  std::complex<double> refresult(std::abs(input)*5.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, refresult);
}

// make sure it also works with a space, as the lexing of M operator requires that the 
// left paren be part of the token, so whitespace is optionally handled via regex
TEST ( ComplexParserTest, singleCharacter_Mop2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("M (1.0+2.0J)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> input(1.0,2.0);
  std::complex<double> refresult(std::abs(input));
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, refresult);
}

// P(number) = phase of number
TEST ( ComplexParserTest, singleCharacter_Pop)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("Ph(1.0+2.0J)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> input(1.0,2.0);
  std::complex<double> refresult(std::arg(input));
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, refresult);
}

// make sure it also works with a space, as the lexing of P operator requires that the 
// left paren be part of the token, so whitespace is optionally handled via regex
TEST ( ComplexParserTest, singleCharacter_Pop2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("Ph (1.0+2.0J)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> input(1.0,2.0);
  std::complex<double> refresult(std::arg(input));
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, refresult);
}


// other complex number oriented operators:
TEST ( ComplexParserTest, complexOps_REop)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("RE(2.0+3.0J)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, 2.0);
}

#if 0
// had to disable this, as IM conflicts with current magnitude operator, IM
TEST ( ComplexParserTest, complexOps_IMop)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("IM(2.0+3.0J)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, 3.0);
}
#endif

TEST ( ComplexParserTest, complexOps_IMGop)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("IMG(2.0+3.0J)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, 3.0);
}

TEST ( ComplexParserTest, complexOps_ABSop)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("ABS(2.0+3.0J)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, std::abs(std::complex<double>(2.0,3.0)));
}


// ABS(parameter) = abs of complex valued parameter
TEST ( ComplexParserTest, singleCharacter_ABSop_Param)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("ABS(par1)"), testGroup);
  testExpression.lexAndParseExpression();

  // setup the parameter
  Teuchos::RCP<Xyce::Util::newExpression> par1Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("1.0+2.0J"), testGroup));
  par1Expression->lexAndParseExpression();
  std::string par1Name = "par1";
  testExpression.attachParameterNode(par1Name,par1Expression);

  std::complex<double> input(1.0,2.0);
  std::complex<double> refresult(std::abs(input));
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, refresult);
}

// ABS(parameter) = abs of complex valued expression
TEST ( ComplexParserTest, singleCharacter_ABSop_Expr)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("ABS(par1 + par2*4.0)"), testGroup);
  testExpression.lexAndParseExpression();

  // setup parameter par1
  Teuchos::RCP<Xyce::Util::newExpression> par1Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("1.0+2.0J"), testGroup));
  par1Expression->lexAndParseExpression();
  std::string par1Name = "par1";
  testExpression.attachParameterNode(par1Name,par1Expression);

  // setup parameter par2
  Teuchos::RCP<Xyce::Util::newExpression> par2Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("1.0+2.0J"), testGroup));
  par2Expression->lexAndParseExpression();
  std::string par2Name = "par2";
  testExpression.attachParameterNode(par2Name,par2Expression);

  std::complex<double> input(1.0,2.0);
  std::complex<double> refresult(std::abs(input)*5.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, refresult);
}

// ABS(parameter) = abs of complex valued expression with a function
TEST ( ComplexParserTest, singleCharacter_ABSop_Func)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("ABS(par1 + F1(4.0))"), testGroup);
  testExpression.lexAndParseExpression();

  // setup parameter par1
  Teuchos::RCP<Xyce::Util::newExpression> par1Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("1.0+2.0J"), testGroup));
  par1Expression->lexAndParseExpression();
  std::string par1Name = "par1";
  testExpression.attachParameterNode(par1Name,par1Expression);

  // setup function F1
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("A*(1.0+2.0J)"), testGroup));

  Xyce::Util::newExpression f1_LHS (std::string("F1(A)"), testGroup);
  f1_LHS.lexAndParseExpression();

  std::vector<std::string> f1ArgStrings ;
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);
  f1Expression->setFunctionArgStringVec (f1ArgStrings);
  f1Expression->lexAndParseExpression();
  std::string f1Name;
  f1_LHS.getFuncPrototypeName(f1Name);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  std::complex<double> input(1.0,2.0);
  std::complex<double> refresult(std::abs(input)*5.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, refresult);
}



TEST ( ComplexParserTest, complexOps_DBop)
{
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("db(i(vb))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, VBval=std::complex<double>(3.0,2.0);
  double refRes = 20.0*std::log10(std::abs(VBval)); 
  solnGroup->setSoln(std::string("vb"),VBval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO( ComplexParserTest, complexOps_DBop)
}

TEST ( ComplexParserTest, complexOps_DBop2)
{
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("db((2.0+3.0J))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0;
  double refRes = 20.0*std::log10(std::abs(std::complex<double>(2.0,3.0))); 

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO( ComplexParserTest, complexOps_DBop)
}
// binary operators
PARSER_SIMPLE_TEST_MACRO(ComplexParserTest, binaryAdd, "(1.0+2.0J)+(3.0+4.0J)", (std::complex<double>(1.0,2.0)+std::complex<double>(3.0,4.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserTest, binaryMinus, "(1.0+2.0J)-(3.0+4.0J)", (std::complex<double>(1.0,2.0)-std::complex<double>(3.0,4.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserTest, binaryMul, "(8.0+2.0J)*(3.0+4.0J)", (std::complex<double>(8.0,2.0)*std::complex<double>(3.0,4.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserTest, binaryDiv, "(8.0+2.0J)/(3.0+4.0J)", (std::complex<double>(8.0,2.0)/std::complex<double>(3.0,4.0)))

// these "extra" tests all came from my needing to add another bunch of lexer rules for processing numbers.
// Mainly, numbers specified as 1.  or .2, where one side of the decimal point doesn't have a number at all.
// I couldn't think of a way to handle them with a single regex (or lexer rule).  
// The regression test that triggered this was Certification_Tests/BUG_80_SON/Bad_PWL_Source.   The test was
// not about this issue, but it contains a couple of expressions like this.
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, binaryAddExtra1, "{4.*PI*1.E-7}", (4*M_PI *1.0e-7) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, binaryAddExtra2, "{.4*PI*.1E-7}", (0.4*M_PI *0.1e-7) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, binaryAddExtra3, "{4.u*PI*1.E-7}", (4e-6*M_PI *1.0e-7) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, binaryAddExtra4, "{.4u*PI*.1E-7}", (0.4e-6*M_PI *0.1e-7) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, binaryAddExtra5, "{4.*PI*1.E-7u}", (4*M_PI *1.0e-7*1.0e-6) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, binaryAddExtra6, "{.4*PI*.1E-7u}", (0.4*M_PI *0.1e-7*1.0e-6) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, binaryAddExtra7, "{4.*PI*1.E-7meg}", (4*M_PI *1.0e-7*1.0e+6) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, binaryAddExtra8, "{.4*PI*.1E-7meg}", (0.4*M_PI *0.1e-7*1.0e+6) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, binaryAddExtra9, "{4.meg*PI*1.E-7}", (4e+6*M_PI *1.0e-7) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, binaryAddExtra10, "{.4meg*PI*.1E-7}", (0.4e+6*M_PI *0.1e-7) )

PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, binaryAddExtra11, "{4.*PI*1.E-7mil}", (4*M_PI *1.0e-7*(25.4e-6)) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, binaryAddExtra12, "{.4*PI*.1E-7mil}", (0.4*M_PI *0.1e-7*(25.4e-6)) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, binaryAddExtra13, "{4.mil*PI*1.E-7}", (4*(25.4e-6)*M_PI *1.0e-7) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, binaryAddExtra14, "{.4mil*PI*.1E-7}", (0.4*(25.4e-6)*M_PI *0.1e-7) )

// simple precendence testing (via binary operators):
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, precedence1, "3.0*2.0+4.0", (3.0*2.0+4.0) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, precedence2, "5.0+4.0/2.0", (5.0+4.0/2.0) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, precedence3, "4.0*6.0/2.0", (4.0*6.0/2.0) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, precedence4, "4.0*(6.0/2.0)", (4.0*(6.0/2.0)) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, precedence5, "1.0/4.0*10.0", (1.0/4.0*10.0) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, precedence6, "1.0/(4.0*10.0)", (1.0/(4.0*10.0)) )

PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, unaryPlus, "+2.0", (2.0) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, unaryMinus, "-2.0", (-2.0) )

PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, phase, "Ph(1.0)", std::arg(1.0) )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, real1, "Re(1.0)", 1.0 )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, real2, "R(1.0)", 1.0 )
// Im cannot work along with IM for imaginary current.  disabling
//PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, imag1, "Im(1.0)", 0.0 )
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, imag2, "Img(1.0)", 0.0 )

PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, int1, "int(11.2423)", 11)
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, int2, "int(-11.2423)", -11)
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, int3, "int  (11.2423)", 11)

PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, uramp1, "uramp(11.2423)", 11.2423)
PARSER_SIMPLE_TEST_MACRO ( ComplexParserTest, uramp2, "uramp(-11.2423)", 0)

// std library functions
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest,sqrt,"sqrt(4.0+3.0J)",  std::sqrt(std::complex<double>(4.0,3.0)))

// necessary to use std::abs here b/c this test often produces opposite sign (even though both call std::sqrt)
PARSER_SIMPLE_TEST_MACRO2(ComplexParserUnaryFuncTest,sqrtNeg, "abs(sqrt(-4.0))", std::abs(std::sqrt(std::complex<double>(-4.0,0.0))))

PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, exp,   "exp(0.5)", std::exp(std::complex<double>(0.5,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, abs,   "abs(-0.5)", std::abs(std::complex<double>(-0.5,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, sin,   "sin(0.5)", std::sin(std::complex<double>(0.5,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, cos,   "cos(0.5)", std::cos(std::complex<double>(0.5,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, acos,  "acos(0.5)", std::acos(std::complex<double>(0.5,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, acosh, "acosh(1.5)", std::acosh(std::complex<double>(1.5,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, asin,  "asin(0.5)", std::asin(std::complex<double>(0.5,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, asinh, "asinh(0.5)", std::asinh(std::complex<double>(0.5,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, atan,  "atan(0.5)", std::atan(std::complex<double>(0.5,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, atanh, "atanh(0.5)", std::atanh(std::complex<double>(0.5,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, cosh,  "cosh(0.5)", std::cosh(std::complex<double>(0.5,0.0)))
//PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, log,   "log(0.5)", std::log(std::complex<double>(0.5,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, log,   "log(0.5)", std::log10(std::complex<double>(0.5,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, log10, "log10(0.5)", std::log10(std::complex<double>(0.5,0.0)))

// use abs b/c this often returns opposite sign
PARSER_SIMPLE_TEST_MACRO2(ComplexParserUnaryFuncTest, log10neg, "abs(log10(-0.5))", std::abs(std::log10(std::complex<double>(-0.5,0.0))))

PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, sinh,  "sinh(0.5)", std::sinh(std::complex<double>(0.5,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, tan,   "tan(0.5)", std::tan(std::complex<double>(0.5,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, tanh,  "tanh(0.5)", std::tanh(std::complex<double>(0.5,0.0)))

PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, pow1,  "pow(2.0,3.0)", std::pow(std::complex<double>(2.0,0.0),std::complex<double>(3.0,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, pow2,  "2.0**3.0", std::pow(std::complex<double>(2.0,0.0),std::complex<double>(3.0,0.0)))

PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, abmpow1,  "pow(2.5,2.0)", std::pow(std::complex<double>(2.5,0.0),std::complex<double>(2.0,0.0)))

PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, abmpow2a,  "pow(2.5,-2.0)", std::pow(std::complex<double>(2.5,0.0),std::complex<double>(-2.0,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, abmpow2b,  "pow(2.5,-3.0)", std::pow(std::complex<double>(2.5,0.0),std::complex<double>(-3.0,0.0)))

PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, abmpow3a,  "pow(2.5,2.1)", std::pow(std::complex<double>(2.5,0.0),std::complex<double>(2.1,0.0)))

PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, abmpow3b,  "pow(-2.5,3.1)", std::pow(std::complex<double>(-2.5,-0.0),std::complex<double>(3.1,0.0)))



PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, pwrs1,  "pwrs(2.0,3.0)" , std::pow(std::complex<double>(2.0,0.0),std::complex<double>(3.0,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, pwrs2,  "pwrs(0.0,3.0)", 0.0)
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, pwrs3,  "pwrs(-2.0,3.0)", -std::pow(std::complex<double>(2.0,0.0),std::complex<double>(3.0,0.0)))

PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, sign1,  "sign(-25,10.25)", 25)
PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, sign2,  "sign(15,-10.25)", -15)

// Hspice only:
//PARSER_SIMPLE_TEST_MACRO(ComplexParserUnaryFuncTest, pow3,  "2.0^3.0", std::pow(2.0,3.0))

// lower case metrix prefix/suffix tests
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, tera,  "3.0t", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, giga,  "5.0g", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, kilo,  "7.0k", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, mili,  "2.0m", 2.0e-3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, mega,  "2.0meg", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, mega2,  "4.0x", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, micro,  "2.0u", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, nano,  "9.0n", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, pico,  "6.0p", 6.0e-12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, femto,  "6.0f", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, mil,  "2.0mil", 2.0*(25.4e-6) )

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, Sec,      "3.0s", 3.0)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, teraSec,  "3.0ts", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, gigaSec,  "5.0gs", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, kiloSec,  "7.0ks", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, miliSec,  "2.0ms", 2.0e-3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, megaSec,  "2.0megs", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, mega2Sec,  "4.0xs", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, microSec,  "2.0us", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, nanoSec,  "9.0ns", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, picoSec,  "6.0ps", 6.0e-12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, femtoSec,  "6.0fs", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, milSec,  "2.0mils", 2.0*(25.4e-6) )

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, teraMeter,  "3.0tm", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, gigaMeter,  "5.0gm", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, kiloMeter,  "7.0km", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, miliMeter,  "2.0mm", 2.0e-3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, megaMeter,  "2.0megm", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, mega2Meter,  "4.0xm", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, microMeter,  "2.0um", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, nanoMeter,  "9.0nm", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, picoMeter,  "6.0pm", 6.0e-12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, femtoMeter,  "6.0fm", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, milMeter,  "2.0milm", 2.0*(25.4e-6) )

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, Henry,      "3.0h", 3.0)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, teraHenry,  "3.0th", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, gigaHenry,  "5.0gh", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, kiloHenry,  "7.0kh", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, miliHenry,  "2.0mh", 2.0e-3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, megaHenry,  "2.0megh", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, mega2Henry,  "4.0xh", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, microHenry,  "2.0uh", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, nanoHenry,  "9.0nh", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, picoHenry,  "6.0ph", 6.0e-12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, femtoHenry,  "6.0fh", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, milHenry,  "2.0milh", 2.0*(25.4e-6) )

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, Hz,      "3.0hz", 3.0)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, teraHz,  "3.0thz", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, gigaHz,  "5.0ghz", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, kiloHz,  "7.0khz", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, miliHz,  "2.0mhz", 2.0e-3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, megaHz,  "2.0meghz", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, mega2Hz,  "4.0xhz", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, microHz,  "2.0uhz", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, nanoHz,  "9.0nhz", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, picoHz,  "6.0phz", 6.0e-12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, femtoHz,  "6.0fhz", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, milHz,  "2.0milhz", 2.0*(25.4e-6) )

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, Volt,      "3.0v", 3.0)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, teraVolt,  "3.0tv", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, gigaVolt,  "5.0gv", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, kiloVolt,  "7.0kv", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, miliVolt,  "2.0mv", 2.0e-3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, megaVolt,  "2.0megv", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, mega2Volt,  "4.0xv", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, microVolt,  "2.0uv", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, nanoVolt,  "9.0nv", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, picoVolt,  "6.0pv", 6.0e-12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, femtoVolt,  "6.0fv", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, milVolt,  "2.0milv", 2.0*(25.4e-6) )

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, Amp,     "3.0a", 3.0)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, teraAmp,  "3.0ta", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, gigaAmp,  "5.0ga", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, kiloAmp,  "7.0ka", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, miliAmp,  "2.0ma", 2.0e-3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, megaAmp,  "2.0mega", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, mega2Amp,  "4.0xa", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, microAmp,  "2.0ua", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, nanoAmp,  "9.0na", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, picoAmp,  "6.0pa", 6.0e-12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, femtoAmp,  "6.0fa", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, milAmp,  "2.0mila", 2.0*(25.4e-6) )

// unit suffixes, which should be ignored, lower case
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, volt,  "7.0v", 7.0 )
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, amp,  "6.0a", 6.0 )
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, sec,  "5.0s", 5.0 )
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, hz,  "5.0hz", 5.0 )
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, henry,  "5.0h", 5.0 )

// upper case metrix prefix/suffix tests
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, TERA,  "3.0T", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, GIGA,  "5.0G", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, KILO,  "7.0K", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MILI,  "2.0M", 2.0e-3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MEGA,  "2.0MEG", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MEGA2,  "4.0X", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MICRO,  "2.0U", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, NANO,  "9.0N", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, PICO,  "6.0P", 6.0e-12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, FEMTO,  "6.0F", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MIL,  "2.0MIL", 2.0*(25.4e-6) )

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, SEC,  "3.0S", 3.0)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, TERASEC,  "3.0TS", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, GIGASEC,  "5.0GS", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, KILOSEC,  "7.0KS", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MILISEC,  "2.0MS", 2.0e-3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MEGASEC,  "2.0MEGS", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MEGA2SEC,  "4.0XS", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MICROSEC,  "2.0US", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, NANOSEC,  "9.0NS", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, PICOSEC,  "6.0PS", 6.0E-12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, FEMTOSEC,  "6.0FS", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MILSEC,  "2.0MILS", 2.0*(25.4e-6) )

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, TERAMETER,  "3.0TM", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, GIGAMETER,  "5.0GM", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, KILOMETER,  "7.0KM", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MILIMETER,  "2.0MM", 2.0e-3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MEGAMETER,  "2.0MEGM", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MEGA2METER,  "4.0XM", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MICROMETER,  "2.0UM", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, NANOMETER,  "9.0NM", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, PICOMETER,  "6.0PM", 6.0E-12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, FEMTOMETER,  "6.0FM", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MILMETER,  "2.0MILM", 2.0*(25.4e-6) )

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, HENRY,  "3.0H", 3.0)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, TERAHENRY,  "3.0TH", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, GIGAHENRY,  "5.0GH", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, KILOHENRY,  "7.0KH", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MILIHENRY,  "2.0MH", 2.0e-3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MEGAHENRY,  "2.0MEGH", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MEGA2HENRY,  "4.0XH", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MICROHENRY,  "2.0UH", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, NANOHENRY,  "9.0NH", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, PICOHENRY,  "6.0PH", 6.0E-12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, FEMTOHENRY,  "6.0FH", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MILHENRY,  "2.0MILH", 2.0*(25.4e-6) )

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, VOLT,  "3.0V", 3.0)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, TERAVOLT,  "3.0TV", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, GIGAVOLT,  "5.0GV", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, KILOVOLT,  "7.0KV", 7.0E+3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MILIVOLT,  "2.0MV", 2.0e-3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MEGAVOLT,  "2.0MEGV", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MEGA2VOLT,  "4.0XV", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MICROVOLT,  "2.0UV", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, NANOVOLT,  "9.0NV", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, PICOVOLT,  "6.0PV", 6.0e-12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, FEMTOVOLT,  "6.0FV", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MILVOLT,  "2.0MILV", 2.0*(25.4e-6) )

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, AMP,  "3.0A", 3.0)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, TERAAMP,  "3.0TA", 3.0e+12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, GIGAAMP,  "5.0GA", 5.0e+9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, KILOAMP,  "7.0KA", 7.0e+3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MILIAMP,  "2.0MA", 2.0e-3)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MEGAAMP,  "2.0MEGA", 2.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MEGA2AMP,  "4.0XA", 4.0e+6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MICROAMP,  "2.0UA", 2.0e-6)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, NANOAMP,  "9.0NA", 9.0*1.0e-9)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, PICOAMP,  "6.0PA", 6.0e-12)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, FEMTOAMP,  "6.0FA", 6.0*1.0e-15)
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, MILAMP,  "2.0MILA", 2.0*(25.4e-6) )

// unit suffixes, which should be ignored, upper case
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, VOLT2,  "4.0V", 4.0 )
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, AMP2,  "3.0A", 3.0 )
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, SEC2,  "2.0S", 2.0 )
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, HZ2,  "2.0HZ", 2.0 )
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, HENRY2,  "2.0HZ", 2.0 )

// lower case metrix prefix/suffix tests
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, sin_tera,  "sin(3.0t)", std::sin(3.0e+12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, sin_giga,  "sin(5.0g)", std::sin(5.0e+9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, sin_kilo,  "sin(7.0k)", std::sin(7.0e+3))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, sin_mili,  "sin(2.0m)", std::sin(2.0e-3))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, sin_mega,  "sin(2.0meg)", std::sin(2.0e+6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, sin_mega2,  "sin(4.0x)", std::sin(4.0e+6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, sin_micro,  "sin(2.0u)", std::sin(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, sin_nano,  "sin(9.0n)", std::sin(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, sin_pico,  "sin(6.0p)", std::sin(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, sin_femto,  "sin(6.0f)", std::sin(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, sin_mil,  "sin(2.0mil)", std::sin(2.0*(25.4e-6)))

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_teraSec,  "exp(3.0e-12ts)", std::exp(3.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_gigaSec,  "exp(5.0e-9gs)", std::exp(5.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_kiloSec,  "exp(7.0e-3ks)", std::exp(7.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_miliSec,  "exp(2.0e+3ms)", std::exp(2.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_megaSec,  "exp(2.0e-6megs)", std::exp(2.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_mega2Sec,  "exp(4.0e-6xs)", std::exp(4.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_microSec,  "exp(2.0us)", std::exp(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_nanoSec,  "exp(9.0ns)", std::exp(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_picoSec,  "exp(6.0ps)", std::exp(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_femtoSec,  "exp(6.0fs)", std::exp(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_milSec,  "exp(2.0mils)", std::exp(2.0*(25.4e-6)))

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_teraMeter,  "exp(3.0e-12tm)", std::exp(3.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_gigaMeter,  "exp(5.0e-9gm)", std::exp(5.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_kiloMeter,  "exp(7.0e-3km)", std::exp(7.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_miliMeter,  "exp(2.0e+3mm)", std::exp(2.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_megaMeter,  "exp(2.0e-6megm)", std::exp(2.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_mega2Meter,  "exp(4.0e-6xm)", std::exp(4.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_microMeter,  "exp(2.0um)", std::exp(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_nanoMeter,  "exp(9.0nm)", std::exp(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_picoMeter,  "exp(6.0pm)", std::exp(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_femtoMeter,  "exp(6.0fm)", std::exp(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_milMeter,  "exp(2.0milm)", std::exp(2.0*(25.4e-6)))

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_teraHz,  "exp(3.0e-12thz)", std::exp(3.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_gigaHz,  "exp(5.0e-9ghz)", std::exp(5.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_kiloHz,  "exp(7.0e-3khz)", std::exp(7.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_miliHz,  "exp(2.0e+3mhz)", std::exp(2.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_megaHz,  "exp(2.0e-6meghz)", std::exp(2.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_mega2Hz,  "exp(4.0e-6xhz)", std::exp(4.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_microHz,  "exp(2.0uhz)", std::exp(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_nanoHz,  "exp(9.0nhz)", std::exp(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_picoHz,  "exp(6.0phz)", std::exp(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_femtoHz,  "exp(6.0fhz)", std::exp(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, exp_milHz,  "exp(2.0milhz)", std::exp(2.0*(25.4e-6)))

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, cos_teraVolt,  "cos(3.0tv)", std::cos(3.0e+12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, cos_gigaVolt,  "cos(5.0gv)", std::cos(5.0e+9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, cos_kiloVolt,  "cos(7.0kv)", std::cos(7.0e+3))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, cos_miliVolt,  "cos(2.0mv)", std::cos(2.0e-3))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, cos_megaVolt,  "cos(2.0megv)", std::cos(2.0e+6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, cos_mega2Volt,  "cos(4.0xv)", std::cos(4.0e+6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, cos_microVolt,  "cos(2.0uv)", std::cos(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, cos_nanoVolt,  "cos(9.0nv)", std::cos(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, cos_picoVolt,  "cos(6.0pv)", std::cos(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, cos_femtoVolt,  "cos(6.0fv)", std::cos(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, cos_milVolt,  "cos(2.0milv)", std::cos(2.0*(25.4e-6)))

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, tan_teraAmp,  "tan(3.0ta)", std::tan(3.0e+12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, tan_gigaAmp,  "tan(5.0ga)", std::tan(5.0e+9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, tan_kiloAmp,  "tan(7.0ka)", std::tan(7.0e+3))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, tan_miliAmp,  "tan(2.0m)", std::tan(2.0e-3))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, tan_megaAmp,  "tan(2.0mega)", std::tan(2.0e+6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, tan_mega2Amp,  "tan(4.0xa)", std::tan(std::complex<double>(4.0e+6,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, tan_microAmp,  "tan(2.0ua)", std::tan(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, tan_nanoAmp,  "tan(9.0na)", std::tan(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, tan_picoAmp,  "tan(6.0pa)", std::tan(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, tan_femtoAmp,  "tan(6.0fa)", std::tan(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, tan_milAmp,  "tan(2.0mila)", std::tan(2.0*(25.4e-6)))

// unit suffixes, which should be ignored, lower case
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, sin_volt,  "sin(2.0v)", std::sin(2.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, sin_amp,  "sin(3.0a)", std::sin(3.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, sin_sec,  "sin(4.0s)", std::sin(4.0))

// upper case
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, SIN_TERA,  "SIN(3.0T)", std::sin(3.0e+12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, SIN_GIGA,  "SIN(5.0G)", std::sin(5.0e+9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, SIN_KILO,  "SIN(7.0K)", std::sin(7.0e+3))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, SIN_MILI,  "SIN(2.0M)", std::sin(2.0e-3))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, SIN_MEGA,  "SIN(2.0MEG)", std::sin(2.0e+6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, SIN_MEGA2,  "SIN(4.0X)", std::sin(4.0e+6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, SIN_MICRO,  "SIN(2.0U)", std::sin(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, SIN_NANO,  "SIN(9.0N)", std::sin(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, SIN_PICO,  "SIN(6.0P)", std::sin(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, SIN_FEMTO,  "SIN(6.0F)", std::sin(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, SIN_MIL,  "SIN(2.0MIL)", std::sin(2.0*(25.4e-6)))

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_TERASEC,  "EXP(3.0E-12TS)", std::exp(3.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_GIGASEC,  "EXP(5.0E-9GS)", std::exp(5.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_KILOSEC,  "EXP(7.0E-3KS)", std::exp(7.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_MILISEC,  "EXP(2.0E+3MS)", std::exp(2.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_MEGASEC,  "EXP(2.0E-6MEGS)", std::exp(2.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_MEGA2SEC,  "EXP(4.0E-6XS)", std::exp(4.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_MICROSEC,  "EXP(2.0US)", std::exp(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_NANOSEC,  "EXP(9.0NS)", std::exp(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_PICOSEC,  "EXP(6.0PS)", std::exp(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_FEMTOSEC,  "EXP(6.0FS)", std::exp(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_MILSEC,  "EXP(2.0MILS)", std::exp(2.0*(25.4e-6)))

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_TERAMETER,  "EXP(3.0E-12TM)", std::exp(3.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_GIGAMETER,  "EXP(5.0E-9GM)", std::exp(5.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_KILOMETER,  "EXP(7.0E-3KM)", std::exp(7.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_MILIMETER,  "EXP(2.0E+3MM)", std::exp(2.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_MEGAMETER,  "EXP(2.0E-6MEGM)", std::exp(2.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_MEGA2METER,  "EXP(4.0E-6XM)", std::exp(4.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_MICROMETER,  "EXP(2.0UM)", std::exp(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_NANOMETER,  "EXP(9.0NM)", std::exp(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_PICOMETER,  "EXP(6.0PM)", std::exp(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_FEMTOMETER,  "EXP(6.0FM)", std::exp(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_MILMETER,  "EXP(2.0MILM)", std::exp(2.0*(25.4e-6)))

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_TERAHZ,  "EXP(3.0E-12THZ)", std::exp(3.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_GIGAHZ,  "EXP(5.0E-9GHZ)", std::exp(5.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_KILOHZ,  "EXP(7.0E-3KHZ)", std::exp(7.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_MILIHZ,  "EXP(2.0E+3MHZ)", std::exp(2.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_MEGAHZ,  "EXP(2.0E-6MEGHZ)", std::exp(2.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_MEGA2HZ,  "EXP(4.0E-6XHZ)", std::exp(4.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_MICROHZ,  "EXP(2.0UHZ)", std::exp(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_NANOHZ,  "EXP(9.0NHZ)", std::exp(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_PICOHZ,  "EXP(6.0PHZ)", std::exp(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_FEMTOHZ,  "EXP(6.0FHZ)", std::exp(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, EXP_MILHZ,  "EXP(2.0MILHZ)", std::exp(2.0*(25.4e-6)))

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, COS_TERAVOLT,  "COS(3.0TV)", std::cos(3.0e+12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, COS_GIGAVOLT,  "COS(5.0GV)", std::cos(5.0e+9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, COS_KILOVOLT,  "COS(7.0KV)", std::cos(7.0e+3))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, COS_MILIVOLT,  "COS(2.0MV)", std::cos(2.0e-3))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, COS_MEGAVOLT,  "COS(2.0MEGV)", std::cos(2.0e+6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, COS_MEGA2VOLT,  "COS(4.0XV)", std::cos(4.0e+6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, COS_MICROVOLT,  "COS(2.0UV)", std::cos(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, COS_NANOVOLT,  "COS(9.0NV)", std::cos(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, COS_PICOVOLT,  "COS(6.0PV)", std::cos(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, COS_FEMTOVOLT,  "COS(6.0FV)", std::cos(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, COS_MILVOLT,  "COS(2.0MILV)", std::cos(2.0*(25.4e-6)))

PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, TAN_TERAAMP,  "TAN(3.0TA)", std::tan(3.0e+12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, TAN_GIGAAMP,  "TAN(5.0GA)", std::tan(5.0e+9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, TAN_KILOAMP,  "TAN(7.0KA)", std::tan(7.0e+3))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, TAN_MILIAMP,  "TAN(2.0MA)", std::tan(2.0e-3))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, TAN_MEGAAMP,  "TAN(2.0MEGA)", std::tan(2.0e+6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, TAN_MEGA2AMP,  "TAN(4.0XA)", std::tan(std::complex<double>(4.0e+6,0.0)))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, TAN_MICROAMP,  "TAN(2.0UA)", std::tan(2.0e-6))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, TAN_NANOAMP,  "TAN(9.0NA)", std::tan(9.0*1.0e-9))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, TAN_PICOAMP,  "TAN(6.0PA)", std::tan(6.0e-12))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, TAN_FEMTOAMP,  "TAN(6.0FA)", std::tan(6.0*1.0e-15))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, TAN_MILAMP,  "TAN(2.0MILA)", std::tan(2.0*(25.4e-6)))

// unit suffixes, which should be ignored, upper case
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, SIN_VOLT,  "SIN(5.0V)", std::sin(5.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, SIN_AMP,  "SIN(6.0A)", std::sin(6.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, SIN_SEC,  "SIN(7.0S)", std::sin(7.0))
PARSER_SIMPLE_TEST_MACRO(ComplexParserSuffixTest, SIN_HZ,  "SIN(7.0HZ)", std::sin(7.0))

// the next 6 tests were all added when I was trying to narrow down the problem
// with the lmod_indmod regression test case, which was failing.
// The issue turned out to be that I didn't include Henrys in the units that 
// could be handled by the lexer.  So, it (sort of) lexed "10mH", but ignored the "H", and
// then threw away everything to the right of it rather than emitting an obvious error.
TEST ( ComplexParserTest, simpleExpression_lmod_indmod1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("10mH*2*(1+0.010*(90-27)+0.926e-4*(90-27)**2)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result(0.0);
  std::complex<double> refres(0.039950588);
  testExpression.evaluateFunction(result);   
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));
  copyExpression.evaluateFunction(result);   
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));
  assignExpression.evaluateFunction(result); 
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));

  OUTPUT_MACRO(ComplexParserTest, simpleExpression_lmod_indmod1)
}

TEST ( ComplexParserTest, simpleExpression_lmod_indmod2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("(90-27)**2)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result(0.0);
  std::complex<double> refres = 63*63;
  testExpression.evaluateFunction(result);

  testExpression.evaluateFunction(result);   
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));
  copyExpression.evaluateFunction(result);   
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));
  assignExpression.evaluateFunction(result); 
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));

  OUTPUT_MACRO(ComplexParserTest, simpleExpression_lmod_indmod2)
}

TEST ( ComplexParserTest, simpleExpression_lmod_indmod3)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("10mH"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result(0.0);
  std::complex<double> refres = 10.0*(1.0e-3);
  testExpression.evaluateFunction(result);

  testExpression.evaluateFunction(result);   
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));
  copyExpression.evaluateFunction(result);   
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));
  assignExpression.evaluateFunction(result); 
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));

  OUTPUT_MACRO(ComplexParserTest, simpleExpression_lmod_indmod3)
}

TEST ( ComplexParserTest, simpleExpression_lmod_indmod4)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("1+0.010*(90-27)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result(0.0);
  std::complex<double> refres = 1.63;
  testExpression.evaluateFunction(result);

  testExpression.evaluateFunction(result);   
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));
  copyExpression.evaluateFunction(result);   
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));
  assignExpression.evaluateFunction(result); 
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));

  OUTPUT_MACRO(ComplexParserTest, simpleExpression_lmod_indmod4)
}

TEST ( ComplexParserTest, simpleExpression_lmod_indmod5)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("0.926e-4*(90-27)**2)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result(0.0);
  std::complex<double> refres = (0.926e-4)*63*63;
  testExpression.evaluateFunction(result);

  testExpression.evaluateFunction(result);   
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));
  copyExpression.evaluateFunction(result);   
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));
  assignExpression.evaluateFunction(result); 
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));

  OUTPUT_MACRO(ComplexParserTest, simpleExpression_lmod_indmod5)
}

TEST ( ComplexParserTest, simpleExpression_lmod_indmod6)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("10mH*2"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result(0.0);
  std::complex<double> refres = 10.0*(1.0e-3)*2.0;
  testExpression.evaluateFunction(result);

  testExpression.evaluateFunction(result);   
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));
  copyExpression.evaluateFunction(result);   
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));
  assignExpression.evaluateFunction(result); 
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refres));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refres));

  OUTPUT_MACRO(ComplexParserTest, simpleExpression_lmod_indmod6)
}

//-------------------------------------------------------------------------------
// source functions:
//-------------------------------------------------------------------------------
// pulse
//-------------------------------------------------------------------------------
TEST ( ComplexParserSourceFuncTest, pulse)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;
  Xyce::Util::newExpression testExpression(std::string("spice_pulse(0.0,1.0,0.0,10e-6,10e-6,0.1e-6,20.1e-6)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
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
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserSourceFuncTest, pulse_breakpoints)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;
  Xyce::Util::newExpression testExpression(std::string("spice_pulse(0.0,1.0,0.0,10e-6,10e-6,0.1e-6,20.1e-6)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( (result-(0.0)), 0.0);

  double tr=10.0e-6, pw=0.1e-6;
  timeDepGroup->setTime(tr+0.5*pw);
  testExpression.evaluateFunction(result);
  EXPECT_EQ( (result-(1.0)), 0.0);

  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  testExpression.getBreakPoints(breakPointTimes);

  bool timeDependent = testExpression.getTimeDependent();
  EXPECT_EQ(timeDependent, true);

  OUTPUT_MACRO(ComplexParserSourceFuncTest,pulse_breakpoints)
}

//-------------------------------------------------------------------------------
// identical to the first pulse test, except that it goes thru a .func
//-------------------------------------------------------------------------------
TEST ( ComplexParserSourceFuncTest, pulse_func)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;
  Xyce::Util::newExpression testExpression(std::string("f1(1.0)"), testGroup);
  testExpression.lexAndParseExpression();

  // this expression is the RHS of a .func statement:  .func F1(A) {A*spice_pulse(...)}
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(
    std::string("A*spice_pulse(0.0,1.0,0.0,10e-6,10e-6,0.1e-6,20.1e-6)"), testGroup));

  Xyce::Util::newExpression f1_LHS (std::string("F1(A)"), testGroup);
  f1_LHS.lexAndParseExpression();

  std::vector<std::string> f1ArgStrings ;
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);
  f1Expression->setFunctionArgStringVec (f1ArgStrings);
  f1Expression->lexAndParseExpression();
  std::string f1Name;
  f1_LHS.getFuncPrototypeName(f1Name);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  std::complex<double> result(0.0);
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

  bool timeDependent = testExpression.getTimeDependent();
  bool copyTimeDependent = copyExpression.getTimeDependent();
  bool assignTimeDependent = assignExpression.getTimeDependent();

  EXPECT_EQ(timeDependent, true);
  EXPECT_EQ(copyTimeDependent, true);
  EXPECT_EQ(assignTimeDependent, true);

  OUTPUT_MACRO(ComplexParserSourceFuncTest,pulse_func)
}

//-------------------------------------------------------------------------------
// sin
//-------------------------------------------------------------------------------
TEST ( ComplexParserSourceFuncTest, sin)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;
  Xyce::Util::newExpression testExpression(std::string("spice_sin(1.65,1.65,10000,0,0,-90)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  int numpoints=100;
  std::complex<double> v0(1.65), va(1.65); 
  double freq(10000), td(0.0), theta(0.0), phase(-90),time(0.0);
  double dt=(1.0/freq)*(1.0/static_cast<double>(numpoints));
  std::vector<std::complex<double>> refRes(numpoints), result(numpoints);
  std::vector<std::complex<double>> copyResult(numpoints), assignResult(numpoints);

  for (int ii=0;ii<numpoints;ii++,time+=dt)
  {
    timeDepGroup->setTime(time); 
    testExpression.evaluateFunction(result[ii]);
    copyExpression.evaluateFunction(copyResult[ii]);
    assignExpression.evaluateFunction(assignResult[ii]);
    refRes[ii] = v0 + va * std::sin(2.0*M_PI*((freq)*time + (phase)/360)) * std::exp( -(time*(theta)));
  }
  for( auto i=0; i<result.size(); i++) 
  {
    EXPECT_NEAR( std::real(result[i]), std::real(refRes[i]), 1e-14);
    EXPECT_NEAR( std::imag(result[i]), std::imag(refRes[i]), 1e-14);
    EXPECT_NEAR( std::real(copyResult[i]), std::real(refRes[i]), 1e-14);
    EXPECT_NEAR( std::imag(copyResult[i]), std::imag(refRes[i]), 1e-14);
    EXPECT_NEAR( std::real(assignResult[i]), std::real(refRes[i]), 1e-14);
    EXPECT_NEAR( std::imag(assignResult[i]), std::imag(refRes[i]), 1e-14);
  }
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserSourceFuncTest, sin_3arg)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;
  Xyce::Util::newExpression testExpression(std::string("spice_sin(1.65,1.65,10000)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  int numpoints=100;
  std::complex<double> v0(1.65), va(1.65); 
  double freq(10000), td(0.0), theta(0.0), phase(0.0),time(0.0);
  double dt=(1.0/freq)*(1.0/static_cast<double>(numpoints));
  std::vector<std::complex<double>> refRes(numpoints), result(numpoints);
  std::vector<std::complex<double>> copyResult(numpoints), assignResult(numpoints);

  for (int ii=0;ii<numpoints;ii++,time+=dt)
  {
    timeDepGroup->setTime(time); 
    testExpression.evaluateFunction(result[ii]);
    copyExpression.evaluateFunction(copyResult[ii]);
    assignExpression.evaluateFunction(assignResult[ii]);
    refRes[ii] = v0 + va * std::sin(2.0*M_PI*((freq)*time + (phase)/360)) * std::exp( -(time*(theta)));
  }
  for( auto i=0; i<result.size(); i++) 
  {
    EXPECT_NEAR( std::real(result[i]), std::real(refRes[i]), 1e-14);
    EXPECT_NEAR( std::imag(result[i]), std::imag(refRes[i]), 1e-14);
    EXPECT_NEAR( std::real(copyResult[i]), std::real(refRes[i]), 1e-14);
    EXPECT_NEAR( std::imag(copyResult[i]), std::imag(refRes[i]), 1e-14);
    EXPECT_NEAR( std::real(assignResult[i]), std::real(refRes[i]), 1e-14);
    EXPECT_NEAR( std::imag(assignResult[i]), std::imag(refRes[i]), 1e-14);
  }
}

//-------------------------------------------------------------------------------
// same as sin test, but thru a .func
//-------------------------------------------------------------------------------
TEST ( ComplexParserSourceFuncTest, sin_func)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;

  Xyce::Util::newExpression testExpression(std::string("f1(1.0)"), testGroup);
  testExpression.lexAndParseExpression();

  // this expression is the RHS of a .func statement:  .func F1(A) {A*spice_sin(...)}
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(
    std::string("A*spice_sin(1.65,1.65,10000,0,0,-90)"), testGroup));

  Xyce::Util::newExpression f1_LHS (std::string("F1(A)"), testGroup);
  f1_LHS.lexAndParseExpression();

  std::vector<std::string> f1ArgStrings ;
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);
  f1Expression->setFunctionArgStringVec (f1ArgStrings);
  f1Expression->lexAndParseExpression();
  std::string f1Name;
  f1_LHS.getFuncPrototypeName(f1Name);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  //testExpression.dumpParseTree(std::cout);

  int numpoints=100;
  double v0(1.65), va(1.65), freq(10000), td(0.0), theta(0.0), phase(-90),time(0.0);
  double dt=(1.0/freq)*(1.0/static_cast<double>(numpoints));

  std::vector<std::complex<double> > refRes(numpoints), result(numpoints);
  std::vector<std::complex<double> > copyResult(numpoints), assignResult(numpoints);

  for (int ii=0;ii<numpoints;ii++,time+=dt)
  {
    timeDepGroup->setTime(time); 
    testExpression.evaluateFunction(result[ii]);
    copyExpression.evaluateFunction(copyResult[ii]);
    assignExpression.evaluateFunction(assignResult[ii]);
    refRes[ii] = v0 + va * std::sin(2.0*M_PI*((freq)*time + (phase)/360)) * std::exp( -(time*(theta)));

    EXPECT_NEAR( std::real(result[ii]), std::real(refRes[ii]), 1e-14);
    EXPECT_NEAR( std::imag(result[ii]), std::imag(refRes[ii]), 1e-14);
    EXPECT_NEAR( std::real(copyResult[ii]),std::real( refRes[ii]), 1e-14);
    EXPECT_NEAR( std::imag(copyResult[ii]), std::imag(refRes[ii]), 1e-14);
    EXPECT_NEAR( std::real(assignResult[ii]), std::real(refRes[ii]), 1e-14);
    EXPECT_NEAR( std::imag(assignResult[ii]), std::imag(refRes[ii]), 1e-14);
  }

  bool timeDependent = testExpression.getTimeDependent();
  bool copyTimeDependent = copyExpression.getTimeDependent();
  bool assignTimeDependent = assignExpression.getTimeDependent();

  EXPECT_EQ(timeDependent, true);
  EXPECT_EQ(copyTimeDependent, true);
  EXPECT_EQ(assignTimeDependent, true);

  OUTPUT_MACRO(ComplexParserSourceFuncTest,sin_func)
}

//-------------------------------------------------------------------------------
// This test is taken from the sources.cir Xyce regression test.  This is the "B2" source in that test.
// The test runs with Xyce if I comment out the B2 source.  So I created this unit test to track down the problem.
// The test ultimately did find the problem and it is now fixed.
TEST ( ComplexParserSourceFuncTest, sin2)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;

  // Note that the original test used a parameter named "v2" for va. So, doing that here.
  Xyce::Util::newExpression testExpression(std::string("spice_sin(v0,v2,frequency,td,theta)-2mv"), testGroup);  // original, without phase
  //Xyce::Util::newExpression testExpression(std::string("spice_sin(v0,v2,frequency,td,theta,phase)-2mv"), testGroup); // works, but must have all 6 args
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

  std::vector<std::complex<double> > refRes(numpoints), result(numpoints);
  std::vector<std::complex<double> > copyResult(numpoints), assignResult(numpoints);

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

    EXPECT_NEAR( std::real(result[ii]), std::real(refRes[ii]), 5e-13);
    EXPECT_NEAR( std::imag(result[ii]), std::imag(refRes[ii]), 5e-13);
    EXPECT_NEAR( std::real(copyResult[ii]), std::real(refRes[ii]), 5e-13);
    EXPECT_NEAR( std::imag(copyResult[ii]), std::imag(refRes[ii]), 5e-13);
    EXPECT_NEAR( std::real(assignResult[ii]), std::real(refRes[ii]), 5e-13);
    EXPECT_NEAR( std::imag(assignResult[ii]), std::imag(refRes[ii]), 5e-13);
  }

  bool timeDependent = testExpression.getTimeDependent();
  bool copyTimeDependent = copyExpression.getTimeDependent();
  bool assignTimeDependent = assignExpression.getTimeDependent();

  EXPECT_EQ(timeDependent, true);
  EXPECT_EQ(copyTimeDependent, true);
  EXPECT_EQ(assignTimeDependent, true);

  OUTPUT_MACRO(ComplexParserSourceFuncTest,sin2)
}

//-------------------------------------------------------------------------------
// exp
//-------------------------------------------------------------------------------
TEST ( ComplexParserSourceFuncTest, exp)
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
  std::complex<double> v1(1.1), v2(2.0); 
  double td1(2e-9), tau1(15e-9), td2(5e-9), tau2(30e-9), time(0.0);
  double dt=2*td2/static_cast<double>(numpoints);
  std::vector<std::complex<double> > refRes(numpoints), result(numpoints);
  std::vector<std::complex<double> > copyResult(numpoints), assignResult(numpoints);
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
  for( auto i=0; i<result.size(); i++) 
  {
    EXPECT_NEAR( std::real(result[i]), std::real(refRes[i]), 1e-14);
    EXPECT_NEAR( std::imag(result[i]), std::imag(refRes[i]), 1e-14);
    EXPECT_NEAR( std::real(copyResult[i]), std::real(refRes[i]), 1e-14);
    EXPECT_NEAR( std::imag(copyResult[i]), std::imag(refRes[i]), 1e-14);
    EXPECT_NEAR( std::real(assignResult[i]), std::real(refRes[i]), 1e-14);
    EXPECT_NEAR( std::imag(assignResult[i]), std::imag(refRes[i]), 1e-14);
  }
}

//-------------------------------------------------------------------------------
// same as exp test, but thru a .func
//-------------------------------------------------------------------------------
TEST ( ComplexParserSourceFuncTest, exp_func)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = timeDepGroup;
  Xyce::Util::newExpression testExpression(std::string("f1(1.0)"), testGroup);
  testExpression.lexAndParseExpression();

  // this expression is the RHS of a .func statement:  .func F1(A) {A*spice_exp(...)}
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("A*spice_exp(1.1,2.0,2e-9,15e-9,5e-9,30e-9)"), testGroup));

  Xyce::Util::newExpression f1_LHS (std::string("F1(A)"), testGroup);
  f1_LHS.lexAndParseExpression();

  std::vector<std::string> f1ArgStrings ;
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);
  f1Expression->setFunctionArgStringVec (f1ArgStrings);
  f1Expression->lexAndParseExpression();
  std::string f1Name;
  f1_LHS.getFuncPrototypeName(f1Name);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  int numpoints=100;
  double v1(1.1), v2(2.0), td1(2e-9), tau1(15e-9), td2(5e-9), tau2(30e-9), time(0.0);
  double dt=2*td2/static_cast<double>(numpoints);

  std::vector<std::complex<double> > refRes(numpoints), result(numpoints);
  std::vector<std::complex<double> > copyResult(numpoints), assignResult(numpoints);

  for (int ii=0;ii<numpoints;ii++,time+=dt)
  {
    timeDepGroup->setTime(time); 
    testExpression.evaluateFunction(result[ii]);
    copyExpression.evaluateFunction(copyResult[ii]);
    assignExpression.evaluateFunction(assignResult[ii]);
    if (time <= td1)  refRes[ii] = v1;
    else if (time <= td2 && time > td1) refRes[ii] = v1 + (v2-v1)*(1.0-std::exp(-(time-td1)/tau1));
    else refRes[ii] = v1 + (v2-v1)*(1.0-std::exp(-(time-td1)/tau1)) + (v1-v2)*(1.0-std::exp(-(time-td2)/tau2)) ;

    EXPECT_NEAR( std::real(result[ii]), std::real(refRes[ii]), 5.0e-13);
    EXPECT_NEAR( std::imag(result[ii]), std::imag(refRes[ii]), 5.0e-13);
    EXPECT_NEAR( std::real(copyResult[ii]), std::real(refRes[ii]), 5.0e-13);
    EXPECT_NEAR( std::imag(copyResult[ii]), std::imag(refRes[ii]), 5.0e-13);
    EXPECT_NEAR( std::real(assignResult[ii]), std::real(refRes[ii]), 5.0e-13);
    EXPECT_NEAR( std::imag(assignResult[ii]), std::imag(refRes[ii]), 5.0e-13);
  }

  OUTPUT_MACRO(ComplexParserSourceFuncTest,exp_func)
}

//-------------------------------------------------------------------------------
// sffm
//-------------------------------------------------------------------------------
TEST ( ComplexParserSourceFuncTest, sffm)
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
  std::complex<double> v0(-0.5), va(2.0); 
  double fc(100e6), mdi(0.3), fs(2.1e6), time(0.0);
  double dt=(1.0/2.1e6) /  static_cast<double>(numpoints);

  std::vector<std::complex<double> > refRes(numpoints), result(numpoints);
  std::vector<std::complex<double> > copyResult(numpoints), assignResult(numpoints);

  for (int ii=0;ii<numpoints;ii++,time+=dt)
  {
    timeDepGroup->setTime(time); 
    testExpression.evaluateFunction(result[ii]);
    copyExpression.evaluateFunction(copyResult[ii]);
    assignExpression.evaluateFunction(assignResult[ii]);
    refRes[ii] = v0 + va * sin((2 * M_PI * fc * time) + mdi * sin (2 * M_PI * fs * time));
    
    EXPECT_NEAR( std::real(result[ii]), std::real(refRes[ii]), 5.0e-13);
    EXPECT_NEAR( std::imag(result[ii]), std::imag(refRes[ii]), 5.0e-13);
    EXPECT_NEAR( std::real(copyResult[ii]), std::real(refRes[ii]), 5.0e-13);
    EXPECT_NEAR( std::imag(copyResult[ii]), std::imag(refRes[ii]), 5.0e-13);
    EXPECT_NEAR( std::real(assignResult[ii]), std::real(refRes[ii]), 5.0e-13);
    EXPECT_NEAR( std::imag(assignResult[ii]), std::imag(refRes[ii]), 5.0e-13);
  }
}

//-------------------------------------------------------------------------------
// same as sffm, but thru a .func
//-------------------------------------------------------------------------------
TEST ( ComplexParserSourceFuncTest, sffm_func)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = timeDepGroup;
  Xyce::Util::newExpression testExpression(std::string("f1(1.0)"), testGroup);
  testExpression.lexAndParseExpression();

  // this expression is the RHS of a .func statement:  .func F1(A) {A*spice_sffm(...)}
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(
    std::string("A*spice_sffm(-0.5,2.0,100e6,0.3,2.1e6)"),testGroup));

  Xyce::Util::newExpression f1_LHS (std::string("F1(A)"), testGroup);
  f1_LHS.lexAndParseExpression();

  std::vector<std::string> f1ArgStrings ;
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);
  f1Expression->setFunctionArgStringVec (f1ArgStrings);
  f1Expression->lexAndParseExpression();
  std::string f1Name;
  f1_LHS.getFuncPrototypeName(f1Name);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  int numpoints=100;
  double v0(-0.5), va(2.0), fc(100e6), mdi(0.3), fs(2.1e6), time(0.0);
  double dt=(1.0/2.1e6) /  static_cast<double>(numpoints);

  std::vector<std::complex<double> > refRes(numpoints), result(numpoints);
  std::vector<std::complex<double> > copyResult(numpoints), assignResult(numpoints);

  for (int ii=0;ii<numpoints;ii++,time+=dt)
  {
    timeDepGroup->setTime(time); 
    testExpression.evaluateFunction(result[ii]);
    copyExpression.evaluateFunction(copyResult[ii]);
    assignExpression.evaluateFunction(assignResult[ii]);
    refRes[ii] = v0 + va * sin((2 * M_PI * fc * time) + mdi * sin (2 * M_PI * fs * time));
    
    EXPECT_NEAR( std::real(result[ii]), std::real(refRes[ii]), 5.0e-13);
    EXPECT_NEAR( std::imag(result[ii]), std::imag(refRes[ii]), 5.0e-13);
    EXPECT_NEAR( std::real(copyResult[ii]), std::real(refRes[ii]), 5.0e-13);
    EXPECT_NEAR( std::imag(copyResult[ii]), std::imag(refRes[ii]), 5.0e-13);
    EXPECT_NEAR( std::real(assignResult[ii]), std::real(refRes[ii]), 5.0e-13);
    EXPECT_NEAR( std::imag(assignResult[ii]), std::imag(refRes[ii]), 5.0e-13);
  }

  OUTPUT_MACRO(ComplexParserSourceFuncTest,sffm_func)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserVoltSolnTest, test0)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("V(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0; 
  std::complex<double> Aval=std::complex<double> (3.0,2.0);
  std::complex<double> refRes = Aval;
  solnGroup->setSoln(std::string("A"),Aval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
}

TEST ( ComplexParserVoltSolnTest, test1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp_dynamic_cast<solnExpressionGroup>(testGroup);
  Xyce::Util::newExpression testExpression(std::string("10.0*V(A)+1.0"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result(0.0,0.0);
  std::complex<double> Aval(3.0,5.0);
  std::complex<double> refRes = 10.0*Aval+1.0;
  solnGroup->setSoln(std::string("A"),Aval);
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(ComplexParserVoltSolnTest, test1)
}

TEST ( ComplexParserVoltSolnTest, test2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp_dynamic_cast<solnExpressionGroup>(testGroup);
  Xyce::Util::newExpression testExpression(std::string("12.0*V(A,B)+7.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result(0.0,0.0);
  std::complex<double> Aval(6.3,2.7);
  std::complex<double> Bval(2.1,8.1);
  std::complex<double> refRes = 12.0*(Aval-Bval)+7.5;
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("B"),Bval);
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(ComplexParserVoltSolnTest, test2)
}

//
TEST ( ComplexParserVoltSolnTest, test3)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("12.0*V(0)+7.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result(0.0,0.0);
  solnGroup->setSoln(std::string("0"),3.0); // this should do nothing
  std::complex<double> refRes = std::complex<double>(7.5,0.0);
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(ComplexParserVoltSolnTest, test3)
}

TEST ( ComplexParserVoltSolnTest, test4)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("12.0*V(A,0)+7.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result(0.0,0.0);
  std::complex<double> Aval(6.3,2.7);
  std::complex<double> refRes = 12.0*(Aval)+7.5;
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("0"),std::complex<double>(7.0,1.0)); // this should do nothing
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(ComplexParserVoltSolnTest, test4)
}

TEST ( ComplexParserVoltSolnTest, test5)
{
  Xyce::Util::preprocessFilter.resize(Xyce::IO::PreprocessType::NUM_PREPROCESS, false);
  Xyce::Util::preprocessFilter[Xyce::IO::PreprocessType::REPLACE_GROUND] = true;

  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("12.0*V(gnd,A)+7.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result(0.0,0.0);
  std::complex<double> Aval(6.3,2.7);
  std::complex<double> refRes = -12.0*(Aval)+7.5;
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("0"),std::complex<double>(7.0,1.0)); // this should do nothing
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(ComplexParserVoltSolnTest, test5)

  Xyce::Util::preprocessFilter.clear(); // reset this for next test
}
//
TEST ( ComplexParserVoltSolnTest, test5b)
{
  Xyce::Util::preprocessFilter.resize(Xyce::IO::PreprocessType::NUM_PREPROCESS, false);
  Xyce::Util::preprocessFilter[Xyce::IO::PreprocessType::REPLACE_GROUND] = true;

  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("12.0*V(gnd!,A)+7.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result(0.0,0.0);
  std::complex<double> Aval(6.3,2.7);
  std::complex<double> refRes = -12.0*(Aval)+7.5;
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("0"),std::complex<double>(7.0,1.0)); // this should do nothing
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(ComplexParserVoltSolnTest, test5b)

  Xyce::Util::preprocessFilter.clear(); // reset this for next test
}

TEST ( ComplexParserVoltSolnTest, test5c)
{
  Xyce::Util::preprocessFilter.resize(Xyce::IO::PreprocessType::NUM_PREPROCESS, false);
  Xyce::Util::preprocessFilter[Xyce::IO::PreprocessType::REPLACE_GROUND] = true;

  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("12.0*V(ground,A)+7.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result(0.0,0.0);
  std::complex<double> Aval(6.3,2.7);
  std::complex<double> refRes = -12.0*(Aval)+7.5;
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("0"),std::complex<double>(7.0,1.0)); // this should do nothing
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(ComplexParserVoltSolnTest, test5c)

  Xyce::Util::preprocessFilter.clear(); // reset this for next test
}

// Test complex .PRINT operators for voltage
TEST ( ComplexParserVoltSolnTest, vr_test0)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("Vr (A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=std::complex<double>(3.0,2.0);
  std::complex<double>  refRes = std::real(Aval);
  solnGroup->setSoln(std::string("A"),Aval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO ( ComplexParserVoltSolnTest, vr_test0)
}

TEST ( ComplexParserVoltSolnTest, vi_test0)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("vI (A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=std::complex<double>(3.0,2.0);
  double refRes = std::imag(Aval);
  solnGroup->setSoln(std::string("A"),Aval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO ( ComplexParserVoltSolnTest, vi_test0)
}

TEST ( ComplexParserVoltSolnTest, vm_test0)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("VM(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=std::complex<double>(3.0,2.0);
  double refRes = std::abs(Aval);
  solnGroup->setSoln(std::string("A"),Aval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO ( ComplexParserVoltSolnTest, vm_test0)
}

TEST ( ComplexParserVoltSolnTest, vp_test0)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("vp(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=std::complex<double>(3.0,2.0);
  double refRes = std::arg(Aval);
  solnGroup->setSoln(std::string("A"),Aval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO ( ComplexParserVoltSolnTest, vp_test0)
}

TEST ( ComplexParserVoltSolnTest, vdb_test0)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("vdb(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=std::complex<double>(3.0,2.0);
  double refRes = 20.0*std::log10(std::abs(Aval)); 
  solnGroup->setSoln(std::string("A"),Aval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO ( ComplexParserVoltSolnTest, vdb_test0)
}

//-------------------------------------------------------------------------------
// weird character voltage node tests.
// to support bug 1034, the NODE needs to recognize a lot of weird characters: 
// from bug 1034: 
// ` ~ ! @ # $ % ^ & - _ + [ ] | \ < > . ? / 
//
// From the INVALID_CHARS/valid_chars_expressions.cir test:
//
//.print DC WIDTH=6 PRECISION=1 v(3) {v(1`)} {v(1~)} {v(1!)} {v(1@)} 
// + {v(1$)} {v(1%)} {v(1^)} {v(1&)} 
// + {v(1*)} {v(1-)} {v(1_)} {v(1+)} {v(1[)} 
// + {v(1])} {v(1|)} {v(1\)} {v(1<)} {v(1>)} {v(1.)} 
// + {v(1?)} {v(1/)} 
// + {v(1#)}
// *+ {v(#)}
// + {v(`)} {v(~)} {v(!)} {v(@)} 
// + {v($)} {v(%)} {v(^)} {v(&)} 
// + {v(*)} {v(-)} {v(_)} {v(+)} {v([)} 
// + {v(])} {v(|)} {v(\)} {v(<)} {v(>)} {v(.)} 
// + {v(?)} {v(/)} 
// + v(X1:1`) v(X1:1~) v(X1:1!) v(X1:1@) 
// + v(X1:1#) v(X1:1$) v(X1:1%) v(X1:1^) v(X1:1&) v(X1:1-)
// + v(X1:1_) v(X1:1+) v(X1:1[) v(X1:1]) v(X1:1|) v(X1:1\) v(X1:1<) 
// + v(X1:1>) v(X1:1.) v(X1:1?) v(X1:1/) 
// + v(X1:`) v(X1:~) v(X1:!) v(X1:@) v(X1:#) v(X1:$) v(X1:%) v(X1:^)
// + v(X1:&) v(X1:-) v(X1:_) v(X1:+) v(X1:[) v(X1:]) v(X1:|) v(X1:\) 
// + v(X1:<) v(X1:>) v(X1:.) v(X1:?) v(X1:/) 
// + {I(R3+)} {P(R3+)} {W(R3+)} {I(R-4)}
// + {IC(Q1+)} {IB(Q1+)} {IE(Q1+)}
//
//-------------------------------------------------------------------------------
TEST ( ComplexParserVoltSolnTest, weirdChar)
{
  Teuchos::RCP<solutionGroup> solnGroup = Teuchos::rcp(new solutionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;


  std::vector<std::string> validNodes = 
  {
    std::string("1`"),
    std::string("1~"),
    std::string("1!"),
    std::string("1@"),
    std::string("1$"),
    std::string("1%"),
    std::string("1^"),
    std::string("1&"),
    std::string("1*"),
    std::string("1-"),
    std::string("1_"),
    std::string("1+"),
    std::string("1["),
    std::string("1]"),
    std::string("1|"),
    std::string("1\\"),
    std::string("1<"),
    std::string("1>"),
    std::string("1."),
    std::string("1?"),
    std::string("1/"),
    std::string("1#"),
    std::string("`"),
    std::string("~"),
    std::string("!"),
    std::string("@"),
    std::string("$"),
    std::string("%"),
    std::string("^"),
    std::string("&"),
    std::string("*"),
    std::string("-"),
    std::string("_"),
    std::string("+"),
    std::string("["),
    std::string("]"),
    std::string("|"),
    std::string("\\"),
    std::string("<"),
    std::string(">"),
    std::string("."),
    std::string("?"),
    std::string("/")
  };

  for (int ii=0;ii< validNodes.size();ii++)
  {
    std::string nodeName = validNodes[ii];

    std::string expressionString = "v(" + nodeName + ")";

    Xyce::Util::newExpression testExpression(expressionString, testGroup);
    testExpression.lexAndParseExpression();

    Xyce::Util::newExpression copyExpression(testExpression); 
    Xyce::Util::newExpression assignExpression; 
    assignExpression = testExpression; 

    std::complex<double> result=std::complex<double>(0.0,0.0); 
    std::complex<double> Aval=std::complex<double>(3.0,4.0);
    std::complex<double> refRes = Aval;
    solnGroup->setSoln(nodeName,Aval);

    testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
    copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
    assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
    OUTPUT_MACRO(ComplexParserVoltSolnTest, weirdChar)
  }
}

//-------------------------------------------------------------------------------
// Test complex .PRINT operators for current
//-------------------------------------------------------------------------------
TEST ( ComplexParserCurrentSolnTest, ir_test0)
{
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("Ir (vb)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, VBval=std::complex<double>(3.0,2.0);
  std::complex<double>  refRes = std::real(VBval);
  solnGroup->setSoln(std::string("vb"),VBval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(ComplexParserCurrentSolnTest, ir_test0)
}

TEST ( ComplexParserCurrentSolnTest, ii_test0)
{
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("iI (vb)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, VBval=std::complex<double>(3.0,2.0);
  double refRes = std::imag(VBval);
  solnGroup->setSoln(std::string("VB"),VBval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(ComplexParserCurrentSolnTest, ii_test0)
}

TEST ( ComplexParserCurrentSolnTest, im_test0)
{
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("IM(vb)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, VBval=std::complex<double>(3.0,2.0);
  double refRes = std::abs(VBval);
  solnGroup->setSoln(std::string("vb"),VBval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(ComplexParserCurrentSolnTest, im_test0)
}

TEST ( ComplexParserVoltSolnTest, ip_test0)
{
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("ip(vb)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, VBval=std::complex<double>(3.0,2.0);
  double refRes = std::arg(VBval);
  solnGroup->setSoln(std::string("vb"),VBval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(ComplexParserVoltSolnTest, ip_test0)
}

TEST ( ComplexParserVoltSolnTest, idb_test0)
{
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("idb(vb)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, VBval=std::complex<double>(3.0,2.0);
  double refRes = 20.0*std::log10(std::abs(VBval)); 
  solnGroup->setSoln(std::string("vb"),VBval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(ComplexParserVoltSolnTest, idb_test0)
}

TEST ( ComplexParserVoltDerivTest, test1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp_dynamic_cast<solnExpressionGroup>(testGroup);
  Xyce::Util::newExpression testExpression(std::string("12.3*V(A)*V(B)+7.5"), testGroup);
  testExpression.lexAndParseExpression();
  //testExpression.resolveExpression(); // this won't do anything, but is the proper order

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result(0.0,0.0), Aval(6.3,2.7), Bval(2.1,8.1);
  std::complex<double> refRes = 12.3*Aval*Bval+7.5;
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("B"),Bval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back(12.3*Bval);
  refDer.push_back(12.3*Aval);
  std::vector<std::complex<double> > derivs;

  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
}

TEST ( ComplexParserVoltDerivTest, test2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp_dynamic_cast<solnExpressionGroup>(testGroup);
  Xyce::Util::newExpression testExpression(std::string("20.0*V(A)**2.0+7.5"), testGroup);
  testExpression.lexAndParseExpression();
  //testExpression.resolveExpression(); // this won't do anything, but is the proper order

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result(0.0,0.0), Aval(6.3,2.7), Bval(2.1,8.1);
  std::complex<double> refRes = 20.0*std::pow(Aval,2.0)+7.5;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back(20.0* ( (2.0/Aval)*std::pow(Aval,2.0)) );

  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result,derivs);   
  EXPECT_NEAR( std::real(result), std::real(refRes), 1e-13); EXPECT_NEAR( std::real(derivs[0]-refDer[0]), std::real(0.0), 1e-13);
  EXPECT_NEAR( std::imag(result), std::imag(refRes), 1e-13); EXPECT_NEAR( std::imag(derivs[0]-refDer[0]), std::imag(0.0), 1e-13);
  copyExpression.evaluate(result,derivs);
  EXPECT_NEAR( std::real(result), std::real(refRes), 1e-13); EXPECT_NEAR( std::real(derivs[0]-refDer[0]), std::real(0.0), 1e-13);
  EXPECT_NEAR( std::imag(result), std::imag(refRes), 1e-13); EXPECT_NEAR( std::imag(derivs[0]-refDer[0]), std::imag(0.0), 1e-13);
  assignExpression.evaluate(result,derivs);
  EXPECT_NEAR( std::real(result), std::real(refRes), 1e-13); EXPECT_NEAR( std::real(derivs[0]-refDer[0]), std::real(0.0), 1e-13);
  EXPECT_NEAR( std::imag(result), std::imag(refRes), 1e-13); EXPECT_NEAR( std::imag(derivs[0]-refDer[0]), std::imag(0.0), 1e-13);
}

TEST ( ComplexParserVoltDerivTest, test3)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp_dynamic_cast<solnExpressionGroup>(testGroup);
  Xyce::Util::newExpression testExpression(std::string("20.0*V(A)*V(A)+7.5"), testGroup);
  testExpression.lexAndParseExpression();
  //testExpression.resolveExpression(); // this won't do anything, but is the proper order

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result(0.0,0.0), Aval(6.3,2.7);
  std::complex<double> refRes = 20.0*Aval*Aval+7.5;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back(20.0*Aval + 20.0*Aval);
  std::vector<std::complex<double> > derivs;

  testExpression.evaluate(result,derivs);
  EXPECT_DOUBLE_EQ(std::real(result), std::real(refRes));
  EXPECT_DOUBLE_EQ(std::imag(result), std::imag(refRes));
  EXPECT_EQ(derivs, refDer);

  copyExpression.evaluate(result,derivs);
  EXPECT_EQ(result, refRes);
  EXPECT_EQ(derivs, refDer);

  assignExpression.evaluate(result,derivs);
  EXPECT_EQ(result, refRes);
  EXPECT_EQ(derivs, refDer);
}

TEST ( ComplexParserVoltDerivTest, test4)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp_dynamic_cast<solnExpressionGroup>(testGroup);
  Xyce::Util::newExpression testExpression(std::string("20.0*V(A)*V(A)+7.5*V(A)"), testGroup);
  testExpression.lexAndParseExpression();
  //testExpression.resolveExpression(); // this won't do anything, but is the proper order

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result(0.0,0.0), Aval(6.3,2.7);
  std::complex<double> refRes = 20.0*Aval*Aval+7.5*Aval;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back( 20.0*Aval + 20.0*Aval + 7.5 );
  std::vector<std::complex<double> > derivs;

  testExpression.evaluate(result,derivs);
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refRes));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refRes));
  EXPECT_EQ( derivs,refDer);

  copyExpression.evaluate(result,derivs);
  EXPECT_EQ( result, refRes );
  EXPECT_EQ( derivs, refDer );

  assignExpression.evaluate(result,derivs);
  EXPECT_EQ( result, refRes );
  EXPECT_EQ( derivs, refDer );
}

TEST ( ComplexParserVoltDerivTest, test5)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp_dynamic_cast<solnExpressionGroup>(testGroup);
  Xyce::Util::newExpression testExpression(std::string("20.0*(V(A)**3.0)+7.5*V(A)"), testGroup);
  testExpression.lexAndParseExpression();
  //testExpression.resolveExpression(); // this won't do anything, but is the proper order

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result(0.0,0.0), Aval(6.3,0.0);
  std::complex<double> refRes = 20.0*std::pow(Aval,3.0)+7.5*Aval;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back( 20.0*(3.0/Aval)*std::pow(Aval,3.0)+7.5 );
  std::vector<std::complex<double> > derivs;

  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
}

TEST ( ComplexParserVoltDerivTest, test6)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp_dynamic_cast<solnExpressionGroup>(testGroup);
  Xyce::Util::newExpression testExpression(std::string("20.0*V(A)*V(A)*V(A)+7.5*V(A)"), testGroup);
  testExpression.lexAndParseExpression();
  //testExpression.resolveExpression(); // this won't do anything, but is the proper order

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result(0.0,0.0), Aval(6.3,0.0);
  std::complex<double> refRes = 20.0*std::pow(Aval,3.0)+7.5*Aval;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back( 20.0*(3.0/Aval)*std::pow(Aval,3.0)+7.5 );
  std::vector<std::complex<double> > derivs;

  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs,refDer);
}

TEST ( ComplexParserVoltDerivTest, test7)
{
  Xyce::Util::preprocessFilter.resize(Xyce::IO::PreprocessType::NUM_PREPROCESS, false);
  Xyce::Util::preprocessFilter[Xyce::IO::PreprocessType::REPLACE_GROUND] = true;

  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*V(A)*V(A,gnd)*(-V(gnd,A))+7.5*V(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=6.3;
  std::complex<double> refRes = 20.0*std::pow(Aval,3.0)+7.5*Aval;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back( 20.0*(3.0/Aval)*std::pow(Aval,3.0)+7.5 );
  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  OUTPUT_MACRO(ComplexParserVoltDerivTest, test7)

  Xyce::Util::preprocessFilter.clear(); // reset this for next test
}

TEST ( ComplexParserVoltDerivTest, test7b)
{
  Xyce::Util::preprocessFilter.resize(Xyce::IO::PreprocessType::NUM_PREPROCESS, false);
  Xyce::Util::preprocessFilter[Xyce::IO::PreprocessType::REPLACE_GROUND] = true;

  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*V(A)*V(A,gnd!)*(-V(gnd!,A))+7.5*V(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=6.3;
  std::complex<double> refRes = 20.0*std::pow(Aval,3.0)+7.5*Aval;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back( 20.0*(3.0/Aval)*std::pow(Aval,3.0)+7.5 );
  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  OUTPUT_MACRO(ComplexParserVoltDerivTest, test7b)

  Xyce::Util::preprocessFilter.clear(); // reset this for next test
}

TEST ( ComplexParserVoltDerivTest, test7c)
{
  Xyce::Util::preprocessFilter.resize(Xyce::IO::PreprocessType::NUM_PREPROCESS, false);
  Xyce::Util::preprocessFilter[Xyce::IO::PreprocessType::REPLACE_GROUND] = true;

  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*V(A)*V(A,ground)*(-V(ground,A))+7.5*V(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=6.3;
  std::complex<double> refRes = 20.0*std::pow(Aval,3.0)+7.5*Aval;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back( 20.0*(3.0/Aval)*std::pow(Aval,3.0)+7.5 );
  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  OUTPUT_MACRO(ComplexParserVoltDerivTest, test7c)

  Xyce::Util::preprocessFilter.clear(); // reset this for next test
}

TEST ( ComplexParserVoltDerivTest, test8)
{
  Xyce::Util::preprocessFilter.resize(Xyce::IO::PreprocessType::NUM_PREPROCESS, false);
  Xyce::Util::preprocessFilter[Xyce::IO::PreprocessType::REPLACE_GROUND] = true;

  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*(V(A,gnd)**3.0)+7.5*V(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=6.3;
  std::complex<double> refRes = 20.0*std::pow(Aval,3.0)+7.5*Aval;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back( 20.0*(3.0/Aval)*std::pow(Aval,3.0)+7.5 );
  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  OUTPUT_MACRO(ComplexParserVoltDerivTest, test8)

  Xyce::Util::preprocessFilter.clear(); // reset this for next test
}

TEST ( ComplexParserVoltDerivTest, test8b)
{
  Xyce::Util::preprocessFilter.resize(Xyce::IO::PreprocessType::NUM_PREPROCESS, false);
  Xyce::Util::preprocessFilter[Xyce::IO::PreprocessType::REPLACE_GROUND] = true;

  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*(V(A,gnd!)**3.0)+7.5*V(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=6.3;
  std::complex<double> refRes = 20.0*std::pow(Aval,3.0)+7.5*Aval;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back( 20.0*(3.0/Aval)*std::pow(Aval,3.0)+7.5 );
  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  OUTPUT_MACRO(ComplexParserVoltDerivTest, test8b)

  Xyce::Util::preprocessFilter.clear(); // reset this for next test
}

TEST ( ComplexParserVoltDerivTest, test8c)
{
  Xyce::Util::preprocessFilter.resize(Xyce::IO::PreprocessType::NUM_PREPROCESS, false);
  Xyce::Util::preprocessFilter[Xyce::IO::PreprocessType::REPLACE_GROUND] = true;

  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*(V(A,ground)**3.0)+7.5*V(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=6.3;
  std::complex<double> refRes = 20.0*std::pow(Aval,3.0)+7.5*Aval;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back( 20.0*(3.0/Aval)*std::pow(Aval,3.0)+7.5 );
  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  OUTPUT_MACRO(ComplexParserVoltDerivTest, test8c)

  Xyce::Util::preprocessFilter.clear(); // reset this for next test
}

TEST ( ComplexParserVoltDerivTest, test9)
{
  Xyce::Util::preprocessFilter.resize(Xyce::IO::PreprocessType::NUM_PREPROCESS, false);
  Xyce::Util::preprocessFilter[Xyce::IO::PreprocessType::REPLACE_GROUND] = true;

  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*((-V(gnd,A))**3.0)+7.5*V(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=6.3;
  std::complex<double> refRes = 20.0*std::pow(Aval,3.0)+7.5*Aval;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back( 20.0*(3.0/Aval)*std::pow(Aval,3.0)+7.5 );
  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  OUTPUT_MACRO(ComplexParserVoltDerivTest, test9)

  Xyce::Util::preprocessFilter.clear(); // reset this for next test
}

TEST ( ComplexParserVoltDerivTest, test9b)
{
  Xyce::Util::preprocessFilter.resize(Xyce::IO::PreprocessType::NUM_PREPROCESS, false);
  Xyce::Util::preprocessFilter[Xyce::IO::PreprocessType::REPLACE_GROUND] = true;

  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*((-V(gnd!,A))**3.0)+7.5*V(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=6.3;
  std::complex<double> refRes = 20.0*std::pow(Aval,3.0)+7.5*Aval;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back( 20.0*(3.0/Aval)*std::pow(Aval,3.0)+7.5 );
  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  OUTPUT_MACRO(ComplexParserVoltDerivTest, test9b)

  Xyce::Util::preprocessFilter.clear(); // reset this for next test
}

TEST ( ComplexParserVoltDerivTest, test9c)
{
  Xyce::Util::preprocessFilter.resize(Xyce::IO::PreprocessType::NUM_PREPROCESS, false);
  Xyce::Util::preprocessFilter[Xyce::IO::PreprocessType::REPLACE_GROUND] = true;

  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*((-V(ground,A))**3.0)+7.5*V(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=6.3;
  std::complex<double> refRes = 20.0*std::pow(Aval,3.0)+7.5*Aval;
  solnGroup->setSoln(std::string("A"),Aval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back( 20.0*(3.0/Aval)*std::pow(Aval,3.0)+7.5 );
  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  OUTPUT_MACRO(ComplexParserVoltDerivTest, test9c)

  Xyce::Util::preprocessFilter.clear(); // reset this for next test
}

TEST ( ComplexParserVoltDerivTest, test10)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("12.3*V(A)/V(B)-7.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=6.3, Bval=2.1;
  std::complex<double> refRes = 12.3*Aval/Bval-7.5;
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("B"),Bval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back(12.3/Bval);
  refDer.push_back(-12.3*Aval/(Bval*Bval));
  std::vector<std::complex<double> > derivs;

  testExpression.evaluate(result,derivs);
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs[0],refDer[0]);
  EXPECT_EQ( derivs[1],refDer[1]);

  copyExpression.evaluate(result,derivs);
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs[0],refDer[0]);
  EXPECT_EQ( derivs[1],refDer[1]);

  assignExpression.evaluate(result,derivs);
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs[0],refDer[0]);
  EXPECT_EQ( derivs[1],refDer[1]);

  OUTPUT_MACRO(ComplexParserVoltDerivTest, test1)
}

TEST ( ComplexParserVoltDerivTest, test11)
{
  Xyce::Util::preprocessFilter.resize(Xyce::IO::PreprocessType::NUM_PREPROCESS, false);
  Xyce::Util::preprocessFilter[Xyce::IO::PreprocessType::REPLACE_GROUND] = true;

  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*((-V(sub1:gnd,sub1:A))**3.0)+7.5*V(sub1:A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=6.3;
  std::complex<double> refRes = 20.0*std::pow(Aval,3.0)+7.5*Aval;
  solnGroup->setSoln(std::string("sub1:A"),Aval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back( 20.0*(3.0/Aval)*std::pow(Aval,3.0)+7.5 );
  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  OUTPUT_MACRO(ComplexParserVoltDerivTest, test11)

  Xyce::Util::preprocessFilter.clear(); // reset this for next test
}

TEST ( ComplexParserVoltDerivTest, test11b)
{
  Xyce::Util::preprocessFilter.resize(Xyce::IO::PreprocessType::NUM_PREPROCESS, false);
  Xyce::Util::preprocessFilter[Xyce::IO::PreprocessType::REPLACE_GROUND] = true;

  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*((-V(sub0:sub1:ground,sub1:A))**3.0)+7.5*V(sub1:A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=6.3;
  std::complex<double> refRes = 20.0*std::pow(Aval,3.0)+7.5*Aval;
  solnGroup->setSoln(std::string("sub1:A"),Aval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back( 20.0*(3.0/Aval)*std::pow(Aval,3.0)+7.5 );
  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  OUTPUT_MACRO(ComplexParserVoltDerivTest, test11b)

  Xyce::Util::preprocessFilter.clear(); // reset this for next test
}

TEST ( ComplexParserVoltDerivTest, test11c)
{
  Xyce::Util::preprocessFilter.resize(Xyce::IO::PreprocessType::NUM_PREPROCESS, false);
  Xyce::Util::preprocessFilter[Xyce::IO::PreprocessType::REPLACE_GROUND] = true;

  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("20.0*((-V(sub0:sub1:sub2:gnd!,sub1:A))**3.0)+7.5*V(sub1:A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=6.3;
  std::complex<double> refRes = 20.0*std::pow(Aval,3.0)+7.5*Aval;
  solnGroup->setSoln(std::string("sub1:A"),Aval);
  std::vector<std::complex<double> > refDer;
  refDer.push_back( 20.0*(3.0/Aval)*std::pow(Aval,3.0)+7.5 );
  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs[0],refDer[0]);
  OUTPUT_MACRO(ComplexParserVoltDerivTest, test11c)

  Xyce::Util::preprocessFilter.clear(); // reset this for next test
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserCurrSolnTest, test1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp_dynamic_cast<currSolnExpressionGroup>(testGroup);
  Xyce::Util::newExpression testExpression(std::string("17.2*I(vb)+8.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result(0.0,0.0), VBval(3.0,1.0);
  std::complex<double> refRes = 17.2*VBval+8.5;
  solnGroup->setSoln(std::string("VB"),VBval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
}

TEST ( ComplexParserPowerTest, test1)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("17.2*P(R1)+8.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result=0.0, R1val=3.0;
  std::complex<double> refRes = 17.2*R1val+8.5;
  solnGroup->setPower(std::string("R1"),R1val);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(ComplexParserPowerTest, test1)
}

TEST ( ComplexParserPowerTest, test2)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("17.2*W(R1)+8.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result=0.0, R1val=3.0;
  std::complex<double> refRes = 17.2*R1val+8.5;
  solnGroup->setPower(std::string("R1"),R1val);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(ComplexParserPowerTest, test2)
}

TEST ( ComplexParserPowerTest, test3)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("1-W(YACC!ACC1)"), testGroup); 
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result=0.0, ACC1val=3.0;
  std::complex<double> refRes = 1.0-ACC1val;
  solnGroup->setPower(std::string("YACC!ACC1"),ACC1val);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(ComplexParserPowerTest, test2)
}

TEST ( ComplexParserCurrDerivTest, test1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp_dynamic_cast<currSolnExpressionGroup>(testGroup);
  Xyce::Util::newExpression testExpression(std::string("17.2*I(vb)+8.5"), testGroup);
  testExpression.lexAndParseExpression();
  //testExpression.resolveExpression(); // this won't do anything, but is the proper order

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result(0.0,0.0), VBval(3.0,1.0);
  std::complex<double> refRes = 17.2*VBval+8.5;
  solnGroup->setSoln(std::string("VB"),VBval);

  std::vector<std::complex<double> > refDer;
  refDer.push_back( 17.2 );
  std::vector<std::complex<double> > derivs;
 
  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
}

TEST ( ComplexParserCurrDerivTest, test2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp_dynamic_cast<currSolnExpressionGroup>(testGroup);
  Xyce::Util::newExpression testExpression(std::string("20.0*I(VB)*I(VB)*I(VB)+7.5*I(VB)"), testGroup);
  testExpression.lexAndParseExpression();
  //testExpression.resolveExpression(); // this won't do anything, but is the proper order

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result(0.0,0.0), VBval(6.3,0.0);
  std::complex<double> refRes = 20.0*std::pow(VBval,3.0)+7.5*VBval;
  solnGroup->setSoln(std::string("VB"),VBval);

  std::vector<std::complex<double> > refDer;
  refDer.push_back( 20.0*(3.0/VBval)*std::pow(VBval,3.0)+7.5 );
  std::vector<std::complex<double> > derivs;

  testExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
  copyExpression.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
  assignExpression.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserLeadCurrTest, test1)
{
  Teuchos::RCP<leadCurrentExpressionGroup> leadCurrentGroup = Teuchos::rcp(new leadCurrentExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = leadCurrentGroup;
  Xyce::Util::newExpression testExpression(std::string("17.2*IG(M1)+8.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, M1val=3.0;
  std::complex<double> refRes = 17.2*M1val+8.5;
  leadCurrentGroup->setCurrentVal(std::string("M1"),std::string("IG"),M1val);

  //testExpression.dumpParseTree(std::cout);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);

  OUTPUT_MACRO(ComplexParserCurrSolnTest, test1)
}

TEST ( ComplexParserLeadCurrTest, test2)
{
  Teuchos::RCP<leadCurrentExpressionGroup> leadCurrentGroup = Teuchos::rcp(new leadCurrentExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = leadCurrentGroup;
  Xyce::Util::newExpression testExpression(std::string("17.2*ID(M1)+8.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, M1val=3.0;
  std::complex<double> refRes = 17.2*M1val+8.5;
  leadCurrentGroup->setCurrentVal(std::string("M1"),std::string("ID"),M1val);

  //testExpression.dumpParseTree(std::cout);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);

  OUTPUT_MACRO(ComplexParserCurrSolnTest, test1)
}

TEST ( ComplexParserLeadCurrTest, test3)
{
  Teuchos::RCP<leadCurrentExpressionGroup> leadCurrentGroup = Teuchos::rcp(new leadCurrentExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = leadCurrentGroup;
  Xyce::Util::newExpression testExpression(std::string("17.2*IS(M1)+8.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, M1val=3.0;
  std::complex<double> refRes = 17.2*M1val+8.5;
  leadCurrentGroup->setCurrentVal(std::string("M1"),std::string("IS"),M1val);

  //testExpression.dumpParseTree(std::cout);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);

  OUTPUT_MACRO(ComplexParserCurrSolnTest, test1)
}

TEST ( ComplexParserInternalDeviceVariableTest, test1)
{
  Teuchos::RCP<internalDevExpressionGroup> intVarGroup = Teuchos::rcp(new internalDevExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = intVarGroup;
  Xyce::Util::newExpression testExpression(std::string("17.2*N(M3:GM)+8.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0; 
  std::complex<double> M3GMval=std::complex<double>(3.0,0.0);
  std::complex<double> refRes = 17.2*M3GMval+8.5;
  intVarGroup->setInternalDeviceVar(std::string("M3:GM"),M3GMval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
 
  OUTPUT_MACRO(ComplexParserCurrSolnTest, test1)
}


// Test complex .PRINT operators for N()
TEST ( ComplexParserInternalDeviceVariableTest, nr_test0)
{
  Teuchos::RCP<internalDevExpressionGroup> intVarGroup = Teuchos::rcp(new internalDevExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = intVarGroup;

  Xyce::Util::newExpression testExpression(std::string("Nr (A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=std::complex<double>(3.0,2.0);
  std::complex<double>  refRes = std::real(Aval);
  intVarGroup->setInternalDeviceVar(std::string("A"),Aval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO ( ComplexParserInternalDeviceVariableTest, nr_test0)
}

TEST ( ComplexParserInternalDeviceVariableTest, ni_test0)
{
  Teuchos::RCP<internalDevExpressionGroup> intVarGroup = Teuchos::rcp(new internalDevExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = intVarGroup;

  Xyce::Util::newExpression testExpression(std::string("nI (A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=std::complex<double>(3.0,2.0);
  double refRes = std::imag(Aval);
  intVarGroup->setInternalDeviceVar(std::string("A"),Aval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO ( ComplexParserInternalDeviceVariableTest, ni_test0)
}

TEST ( ComplexParserInternalDeviceVariableTest, nm_test0)
{
  Teuchos::RCP<internalDevExpressionGroup> intVarGroup = Teuchos::rcp(new internalDevExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = intVarGroup;

  Xyce::Util::newExpression testExpression(std::string("nM(A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=std::complex<double>(3.0,2.0);
  double refRes = std::abs(Aval);
  intVarGroup->setInternalDeviceVar(std::string("A"),Aval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO ( ComplexParserInternalDeviceVariableTest, nm_test0)
}

TEST ( ComplexParserInternalDeviceVariableTest, np_test0)
{
  Teuchos::RCP<internalDevExpressionGroup> intVarGroup = Teuchos::rcp(new internalDevExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = intVarGroup;

  Xyce::Util::newExpression testExpression(std::string("Np  (A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=std::complex<double>(3.0,2.0);
  double refRes = std::arg(Aval);
  intVarGroup->setInternalDeviceVar(std::string("A"),Aval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO ( ComplexParserInternalDeviceVariableTest, np_test0)
}

TEST ( ComplexParserInternalDeviceVariableTest, ndb_test0)
{
  Teuchos::RCP<internalDevExpressionGroup> intVarGroup = Teuchos::rcp(new internalDevExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = intVarGroup;

  Xyce::Util::newExpression testExpression(std::string("nDb  (A)"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=std::complex<double>(3.0,2.0);
  double refRes = 20.0*std::log10(std::abs(Aval)); 
  intVarGroup->setInternalDeviceVar(std::string("A"),Aval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO ( ComplexParserInternalDeviceVariableTest, ndb_test0)
}

TEST ( ComplexParserInternalDeviceVariableTest, testConflict)
{
  Teuchos::RCP<internalDevExpressionGroup> intVarGroup = Teuchos::rcp(new internalDevExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = intVarGroup;
  Xyce::Util::newExpression testExpression(std::string("N+17.2*N(M3:GM)+8.5"), testGroup);
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Teuchos::RCP<Xyce::Util::newExpression> nExpression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("-0.5"), testGroup));
  nExpression->lexAndParseExpression();
  std::string nName = "N";
  testExpression.attachParameterNode(nName,nExpression);

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, M3GMval=3.0;
  std::complex<double> refRes = 17.2*M3GMval+8.5 + (-0.5);
  intVarGroup->setInternalDeviceVar(std::string("M3:GM"),M3GMval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result); 
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
 
  OUTPUT_MACRO(ComplexParserInternalDeviceVariableTest, test1)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserNoiseTest, dno_test)
{
  Teuchos::RCP<noiseExpressionGroup> noiseVarGroup = Teuchos::rcp(new noiseExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = noiseVarGroup;
  Xyce::Util::newExpression testExpression(std::string("17.2*DNO(RES1)+8.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, RES1val=std::complex<double>(3.0,0.0);
  std::complex<double>  refRes = 17.2*RES1val+8.5;
  std::vector<std::string> nameVec;
  nameVec.push_back(std::string("RES1"));
  noiseVarGroup->setDnoNoiseDeviceVar(nameVec,RES1val);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result); 
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
 
  OUTPUT_MACRO(ComplexParserNoiseTest, dno_test)
}

TEST ( ComplexParserNoiseTest, dni_test)
{
  Teuchos::RCP<noiseExpressionGroup> noiseVarGroup = Teuchos::rcp(new noiseExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = noiseVarGroup;
  Xyce::Util::newExpression testExpression(std::string("17.2*DNI(RES1)+8.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, RES1val=std::complex<double>(3.0,0.0);
  std::complex<double>  refRes = 17.2*RES1val+8.5;
  std::vector<std::string> nameVec;
  nameVec.push_back(std::string("RES1"));
  noiseVarGroup->setDniNoiseDeviceVar(nameVec,RES1val);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result); 
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
 
  OUTPUT_MACRO(ComplexParserNoiseTest, dni_test)
}

TEST ( ComplexParserNoiseTest, onoise_test)
{
  Teuchos::RCP<noiseExpressionGroup> noiseVarGroup = Teuchos::rcp(new noiseExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = noiseVarGroup;
  Xyce::Util::newExpression testExpression(std::string("17.2*ONOISE+8.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, ONOISEval=std::complex<double>(3.0,0.0);
  std::complex<double>  refRes = 17.2*ONOISEval+8.5;
  noiseVarGroup->setONoise(ONOISEval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result); 
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
 
  OUTPUT_MACRO(ComplexParserNoiseTest, onoise_test)
}

TEST ( ComplexParserNoiseTest, inoise_test)
{
  Teuchos::RCP<noiseExpressionGroup> noiseVarGroup = Teuchos::rcp(new noiseExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = noiseVarGroup;
  Xyce::Util::newExpression testExpression(std::string("17.2*INOISE+8.5"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, INOISEval=std::complex<double>(3.0,0.0);
  std::complex<double>  refRes = 17.2*INOISEval+8.5;
  noiseVarGroup->setINoise(INOISEval);

  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result); 
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refRes);
 
  OUTPUT_MACRO(ComplexParserNoiseTest, inoise_test)
}

//-------------------------------------------------------------------------------
// .func tests
//-------------------------------------------------------------------------------
TEST ( ComplexParserFuncTest, test1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp_dynamic_cast<testExpressionGroupWithFuncSupport>(testGroup);

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

  //std::vector<std::string> f1ArgStrings = { std::string("A"), std::string("B") };
  f1Expression->setFunctionArgStringVec ( f1ArgStrings );
  f1Expression->lexAndParseExpression();

  //std::string f1Name = "F1";
  std::string f1Name; 
  f1_LHS.getFuncPrototypeName(f1Name);
  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result;
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, std::complex<double>(5.0,0.0) );
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, std::complex<double>(5.0,0.0) );
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, std::complex<double>(5.0,0.0) );
}

//-------------------------------------------------------------------------------
// this test is inspired by the BUG_547_SON/mb_orig.cir test case, which seemed
// to have trouble with the function name DC_AC.
TEST ( ComplexParserFuncTest, test_mb_orig)
{
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression will use the .func DM_AC.
  Xyce::Util::newExpression testExpression(std::string("DM_AC(2,3)"), testGroup);
  testExpression.lexAndParseExpression();

  // this expression is the RHS of a .func statement:  .func DM_AC(A,B) {A+B}
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("A+B"), testGroup) );

  // I originally had this set up so that the calling code would manually set the 
  // vector of prototype function arguments, as well as the name of the function 
  // itself.  But in a code like Xyce, that isn't how it is likely to work. The 
  // function prototype DM_AC(A,B) has to be parsed, and the appropriate information 
  // pulled out of it.  In Xyce, the old expression library is used to parse the 
  // prototype(LHS), so attempting same here.
  Xyce::Util::newExpression f1_LHS (std::string("DM_AC(A,B)"), testGroup);
  f1_LHS.lexAndParseExpression();

  std::vector<std::string> f1ArgStrings ;
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);
  f1Expression->setFunctionArgStringVec (f1ArgStrings);
  // during lex/parse, this vector of arg strings will be compared to any
  // param classes.  If it finds them, then they will be placed in the
  // functionArgOpVec object, which is used below, in the call to "setFuncArgs".
  f1Expression->lexAndParseExpression();

  // now parse the function name from the prototype
  std::string f1Name;
  f1_LHS.getFuncPrototypeName(f1Name);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result;
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, 5.0 );
  OUTPUT_MACRO(ComplexParserFuncTest, test_mb_orig)
}

//-------------------------------------------------------------------------------
// this test is inspired by the BUG_547_SON/mb_orig.cir test case, which seemed
// to have trouble with the function name DC_AC, among other things
TEST ( ComplexParserFuncTest, test_mb_orig2)
{
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression will use the .func DM_AC.
  Xyce::Util::newExpression testExpression(std::string("(0 + I(B2)*I(V2)*DM_AC(v(z)))"), testGroup);
  testExpression.lexAndParseExpression();

  // this expression is the RHS of a .func statement:  .func DM_AC(x) {TABLE(x,0,0)}
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("TABLE(x,0,0)"), testGroup) );

  // I originally had this set up so that the calling code would manually set the 
  // vector of prototype function arguments, as well as the name of the function 
  // itself.  But in a code like Xyce, that isn't how it is likely to work. The 
  // function prototype DM_AC(A,B) has to be parsed, and the appropriate information 
  // pulled out of it.  In Xyce, the old expression library is used to parse the 
  // prototype(LHS), so attempting same here.
  Xyce::Util::newExpression f1_LHS (std::string("DM_AC(x)"), testGroup);
  f1_LHS.lexAndParseExpression();

  std::vector<std::string> f1ArgStrings ;
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);
  f1Expression->setFunctionArgStringVec (f1ArgStrings);
  // during lex/parse, this vector of arg strings will be compared to any
  // param classes.  If it finds them, then they will be placed in the
  // functionArgOpVec object, which is used below, in the call to "setFuncArgs".
  f1Expression->lexAndParseExpression();

  // now parse the function name from the prototype
  std::string f1Name;
  f1_LHS.getFuncPrototypeName(f1Name);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  double b2 = 7.0;
  double v2 = 3.0;
  double vz = 1.0;
  funcGroup->setSoln(std::string("v2"),v2);
  funcGroup->setSoln(std::string("z"),vz);
  funcGroup->setSoln(std::string("b2"),b2);

  std::complex<double> result;
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, 0.0 );
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, 0.0 );
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, 0.0 );
  OUTPUT_MACRO(ComplexParserFuncTest, test_mb_orig2)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserFuncTest, test_underscoreName)
{
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression will use the .func f1.
  Xyce::Util::newExpression testExpression(std::string("_F1(2,3)"), testGroup);
  testExpression.lexAndParseExpression();

  // this expression is the RHS of a .func statement:  .func _F1(A,B) {A+B}
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("A+B"), testGroup) );

  // I originally had this set up so that the calling code would manually set the 
  // vector of prototype function arguments, as well as the name of the function 
  // itself.  But in a code like Xyce, that isn't how it is likely to work. The 
  // function prototype F1(A,B) has to be parsed, and the appropriate information 
  // pulled out of it.  In Xyce, the old expression library is used to parse the 
  // prototype(LHS), so attempting same here.
  Xyce::Util::newExpression f1_LHS (std::string("_F1(A,B)"), testGroup);
  f1_LHS.lexAndParseExpression();

  std::vector<std::string> f1ArgStrings ;
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);
  f1Expression->setFunctionArgStringVec (f1ArgStrings);
  // during lex/parse, this vector of arg strings will be compared to any
  // param classes.  If it finds them, then they will be placed in the
  // functionArgOpVec object, which is used below, in the call to "setFuncArgs".
  f1Expression->lexAndParseExpression();

  // now parse the function name from the prototype
  std::string f1Name;
  f1_LHS.getFuncPrototypeName(f1Name);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result;
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, 5.0 );
  OUTPUT_MACRO(ComplexParserFuncTest, test_underscoreName)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserFuncTest, test_poundSymbolName)
{
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression will use the .func f1.
  Xyce::Util::newExpression testExpression(std::string("#F1(2,3)"), testGroup);
  testExpression.lexAndParseExpression();

  // this expression is the RHS of a .func statement:  .func #F1(A,B) {A+B}
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("A+B"), testGroup) );

  // I originally had this set up so that the calling code would manually set the 
  // vector of prototype function arguments, as well as the name of the function 
  // itself.  But in a code like Xyce, that isn't how it is likely to work. The 
  // function prototype F1(A,B) has to be parsed, and the appropriate information 
  // pulled out of it.  In Xyce, the old expression library is used to parse the 
  // prototype(LHS), so attempting same here.
  Xyce::Util::newExpression f1_LHS (std::string("#F1(A,B)"), testGroup);
  f1_LHS.lexAndParseExpression();

  std::vector<std::string> f1ArgStrings ;
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);
  f1Expression->setFunctionArgStringVec (f1ArgStrings);
  // during lex/parse, this vector of arg strings will be compared to any
  // param classes.  If it finds them, then they will be placed in the
  // functionArgOpVec object, which is used below, in the call to "setFuncArgs".
  f1Expression->lexAndParseExpression();

  // now parse the function name from the prototype
  std::string f1Name;
  f1_LHS.getFuncPrototypeName(f1Name);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result;
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, 5.0 );
  OUTPUT_MACRO(ComplexParserFuncTest, test_poundSymbolName)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserFuncTest, test_atSymbolName)
{
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression will use the .func f1.
  Xyce::Util::newExpression testExpression(std::string("@F1(2,3)"), testGroup);
  testExpression.lexAndParseExpression();

  // this expression is the RHS of a .func statement:  .func @F1(A,B) {A+B}
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("A+B"), testGroup) );

  // I originally had this set up so that the calling code would manually set the 
  // vector of prototype function arguments, as well as the name of the function 
  // itself.  But in a code like Xyce, that isn't how it is likely to work. The 
  // function prototype F1(A,B) has to be parsed, and the appropriate information 
  // pulled out of it.  In Xyce, the old expression library is used to parse the 
  // prototype(LHS), so attempting same here.
  Xyce::Util::newExpression f1_LHS (std::string("@F1(A,B)"), testGroup);
  f1_LHS.lexAndParseExpression();

  std::vector<std::string> f1ArgStrings ;
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);
  f1Expression->setFunctionArgStringVec (f1ArgStrings);
  // during lex/parse, this vector of arg strings will be compared to any
  // param classes.  If it finds them, then they will be placed in the
  // functionArgOpVec object, which is used below, in the call to "setFuncArgs".
  f1Expression->lexAndParseExpression();

  // now parse the function name from the prototype
  std::string f1Name;
  f1_LHS.getFuncPrototypeName(f1Name);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result;
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, 5.0 );
  OUTPUT_MACRO(ComplexParserFuncTest, test_poundSymbolName)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserFuncTest, test_backtickName)
{
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression will use the .func f1.
  Xyce::Util::newExpression testExpression(std::string("`F1(2,3)"), testGroup);
  testExpression.lexAndParseExpression();

  // this expression is the RHS of a .func statement:  .func `F1(A,B) {A+B}
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("A+B"), testGroup) );

  // I originally had this set up so that the calling code would manually set the 
  // vector of prototype function arguments, as well as the name of the function 
  // itself.  But in a code like Xyce, that isn't how it is likely to work. The 
  // function prototype F1(A,B) has to be parsed, and the appropriate information 
  // pulled out of it.  In Xyce, the old expression library is used to parse the 
  // prototype(LHS), so attempting same here.
  Xyce::Util::newExpression f1_LHS (std::string("`F1(A,B)"), testGroup);
  f1_LHS.lexAndParseExpression();

  std::vector<std::string> f1ArgStrings ;
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);
  f1Expression->setFunctionArgStringVec (f1ArgStrings);
  // during lex/parse, this vector of arg strings will be compared to any
  // param classes.  If it finds them, then they will be placed in the
  // functionArgOpVec object, which is used below, in the call to "setFuncArgs".
  f1Expression->lexAndParseExpression();

  // now parse the function name from the prototype
  std::string f1Name;
  f1_LHS.getFuncPrototypeName(f1Name);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result;
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, 5.0 );
  OUTPUT_MACRO(ComplexParserFuncTest, test_backtickName)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserFuncTest, test1_multipleLexParse)
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

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result;
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, 5.0 );
  OUTPUT_MACRO(ComplexParserFuncTest, test1)
}

// tests are taken from the "ternary_precedence.cir" Xyce regression test
TEST ( ComplexParserTernary_precedence, simple)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp_dynamic_cast<testExpressionGroupWithFuncSupport>(testGroup);

  // this expression is the RHS of a .func statement:  .func simple(X) {x>0?2*x :0}
  Teuchos::RCP<Xyce::Util::newExpression> simpleExpression
    = Teuchos::rcp(new Xyce::Util::newExpression(std::string("x>0?2*x :0"), testGroup) );

  //std::vector<std::string> simpleArgStrings = { std::string("x") }; // CASE MATTERS!!!  EEK
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

    std::complex<double> result;
    simpleTrue.evaluateFunction(result);       EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
    copySimpleTrue.evaluateFunction(result);   EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
    assignSimpleTrue.evaluateFunction(result); EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
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

    std::complex<double> result;
    simpleFalse.evaluateFunction(result);       EXPECT_EQ( result, std::complex<double>(0.0,0.0) );
    copySimpleFalse.evaluateFunction(result);   EXPECT_EQ( result, std::complex<double>(0.0,0.0) );
    assignSimpleFalse.evaluateFunction(result); EXPECT_EQ( result, std::complex<double>(0.0,0.0) );
  }
}

TEST ( ComplexParserTernary_precedence, precplus)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp_dynamic_cast<testExpressionGroupWithFuncSupport>(testGroup);

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

    std::complex<double> result;
    precplusTrue.evaluateFunction(result);       EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
    copyPrecplusTrue.evaluateFunction(result);   EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
    assignPrecplusTrue.evaluateFunction(result); EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
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

    std::complex<double> result;
    precplusFalse.evaluateFunction(result);       EXPECT_EQ( result, std::complex<double>(2.0,0.0) );
    copyPrecplusFalse.evaluateFunction(result);   EXPECT_EQ( result, std::complex<double>(2.0,0.0) );
    assignPrecplusFalse.evaluateFunction(result); EXPECT_EQ( result, std::complex<double>(2.0,0.0) );
  }
}

TEST ( ComplexParserTernary_precedence, precplusparen)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp_dynamic_cast<testExpressionGroupWithFuncSupport>(testGroup);

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

    std::complex<double> result;
    precplusparenTrue.evaluateFunction(result);       EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
    copyPrecplusparenTrue.evaluateFunction(result);   EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
    assignPrecplusparenTrue.evaluateFunction(result); EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
    OUTPUT_MACRO2(ComplexParserTernary_precedence, precplusparen, precplusparenTrue) 
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

    std::complex<double> result;
    precplusparenFalse.evaluateFunction(result);       EXPECT_EQ( result, std::complex<double>(2.0,0.0) );
    copyPrecplusparenFalse.evaluateFunction(result);   EXPECT_EQ( result, std::complex<double>(2.0,0.0) );
    assignPrecplusparenFalse.evaluateFunction(result); EXPECT_EQ( result, std::complex<double>(2.0,0.0) );
    OUTPUT_MACRO2(ComplexParserTernary_precedence, precplusparen, precplusparenFalse) 
  }
}

TEST ( ComplexParserTernary_precedence, simpleif)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp_dynamic_cast<testExpressionGroupWithFuncSupport>(testGroup);

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

    std::complex<double> result;
    simpleifTrue.evaluateFunction(result);       EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
    copySimpleifTrue.evaluateFunction(result);   EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
    assignSimpleifTrue.evaluateFunction(result); EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
    OUTPUT_MACRO2(ComplexParserTernary_precedence, simpleif, simpleifTrue) 
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

    std::complex<double> result;
    simpleifFalse.evaluateFunction(result);       EXPECT_EQ( result, std::complex<double>(0.0,0.0) );
    copySimpleifFalse.evaluateFunction(result);   EXPECT_EQ( result, std::complex<double>(0.0,0.0) );
    assignSimpleifFalse.evaluateFunction(result); EXPECT_EQ( result, std::complex<double>(0.0,0.0) );
    OUTPUT_MACRO2(ComplexParserTernary_precedence, simpleif, simpleifFalse) 
  }
}

TEST ( ComplexParserTernary_precedence, precplusif)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp_dynamic_cast<testExpressionGroupWithFuncSupport>(testGroup);

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

    std::complex<double> result;
    precplusifTrue.evaluateFunction(result);       EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
    copyPrecplusifTrue.evaluateFunction(result);   EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
    assignPrecplusifTrue.evaluateFunction(result); EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
    OUTPUT_MACRO2(ComplexParserTernary_precedence, precplusif, precplusifTrue) 
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

    std::complex<double> result;
    precplusifFalse.evaluateFunction(result);       EXPECT_EQ( result, std::complex<double>(2.0,0.0) );
    copyPrecplusifFalse.evaluateFunction(result);   EXPECT_EQ( result, std::complex<double>(2.0,0.0) );
    assignPrecplusifFalse.evaluateFunction(result); EXPECT_EQ( result, std::complex<double>(2.0,0.0) );
    OUTPUT_MACRO2(ComplexParserTernary_precedence, precplusif, precplusifFalse) 
  }
}

TEST ( ComplexParserTernary_precedence, precplusparenif)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp_dynamic_cast<testExpressionGroupWithFuncSupport>(testGroup);

  // this expression is the RHS of a .func statement:  .func precplusparenif(X) {if(((1+x)>0),(2*x),(0+2))}
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

    std::complex<double> result;
    precplusparenifTrue.evaluateFunction(result);       EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
    copyPrecplusparenifTrue.evaluateFunction(result);   EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
    assignPrecplusparenifTrue.evaluateFunction(result); EXPECT_EQ( result, std::complex<double>(8.0,0.0) );
    OUTPUT_MACRO2(ComplexParserTernary_precedence, precplusparenif, precplusparenifTrue) 
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

    std::complex<double> result;
    precplusparenifFalse.evaluateFunction(result);       EXPECT_EQ( result, std::complex<double>(2.0,0.0) );
    copyPrecplusparenifFalse.evaluateFunction(result);   EXPECT_EQ( result, std::complex<double>(2.0,0.0) );
    assignPrecplusparenifFalse.evaluateFunction(result); EXPECT_EQ( result, std::complex<double>(2.0,0.0) );
    OUTPUT_MACRO2(ComplexParserTernary_precedence, precplusparenif, precplusparenifFalse) 
  }
}

TEST ( ComplexParserFuncTest, longArgList)
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

    std::complex<double> result;
    crazyTrue.evaluateFunction(result);       EXPECT_EQ( result, 4.0 );
    copyCrazyTrue.evaluateFunction(result);   EXPECT_EQ( result, 4.0 );
    assignCrazyTrue.evaluateFunction(result); EXPECT_EQ( result, 4.0 );
    OUTPUT_MACRO2(ComplexParserFuncTest, longArgList, crazyTrue) 
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

    std::complex<double> result;
    crazyFalse.evaluateFunction(result);       EXPECT_EQ( result, std::complex<double>(-4.0,0.0) );
    copyCrazyFalse.evaluateFunction(result);   EXPECT_EQ( result, std::complex<double>(-4.0,0.0) );
    assignCrazyFalse.evaluateFunction(result); EXPECT_EQ( result, std::complex<double>(-4.0,0.0) );
    OUTPUT_MACRO2(ComplexParserFuncTest, longArgList, crazyFalse) 
  }
}

//-------------------------------------------------------------------------------
// tests are taken from the "ifstatement.cir" Xyce regression test
//-------------------------------------------------------------------------------
TEST ( ComplexParserIfstatement, ifmin_ifmax_func)
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

  std::string ifmaxName;
  std::string ifminName;

  ifmax_LHS.getFuncPrototypeName (ifmaxName);
  ifmin_LHS.getFuncPrototypeName (ifminName);

  // these expressions uses the .func ifmin and ifmax
  {
    Xyce::Util::newExpression e3(std::string("ifmax(ifmin(-I(V2), 2.5), 1.5)"), baseGroup);
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
    std::vector<std::complex<double> > refRes(numpoints), result(numpoints);
    std::vector<std::complex<double> > copyResult(numpoints), assignResult(numpoints);
    for (int ii=0;ii<numpoints;ii++,time+=dt)
    {
      std::complex<double> V0=std::complex<double>(2.0,0.0); 
      std::complex<double> VA=std::complex<double>(1.0,0.0);
      double FREQ=1, mpi = M_PI;
      std::complex<double> v1 = std::complex<double>(0.0,0.0);
      if (time <= 0) { v1 = (V0) + (VA) * std::sin (0.0); }
      else { v1 = (V0) + (VA) * std::sin (2.0*mpi*((FREQ)*time + (0.0)/360)) * std::exp( -(time*(0.0))); }
      std::complex<double> v11 = 2.0*v1;
      std::complex<double> v2 = -0.5*v11;

      ifGroup->setTime(time);
      ifGroup->setSoln(std::string("v2"),v2);
      e3.evaluateFunction(result[ii]);
      copy_e3.evaluateFunction(copyResult[ii]);
      assign_e3.evaluateFunction(assignResult[ii]);
      refRes[ii] = std::max(std::min(-v2.real(),2.5),1.5);
    }
    EXPECT_EQ( result, refRes);
    EXPECT_EQ( copyResult, refRes);
    EXPECT_EQ( assignResult, refRes);
    OUTPUT_MACRO2(ComplexParserIfstatement, ifmin_ifmax_func, e3) 
  }
}

TEST ( ComplexParserIfstatement, simple_nested_func)
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

  std::string tripleItName;//="tripleIt";
  std::string doubleItName;//="doubleIt";
  doubleIt_LHS.getFuncPrototypeName(doubleItName);
  tripleIt_LHS.getFuncPrototypeName(tripleItName);

  // these expressions uses the .func doubleIt and tripleIt
  {
    Xyce::Util::newExpression e3(std::string("tripleIt(doubleIt(-I(V2)))"), baseGroup);
    e3.lexAndParseExpression();
    e3.attachFunctionNode(tripleItName, tripleIt);
    e3.attachFunctionNode(doubleItName, doubleIt);

    Xyce::Util::newExpression copy_e3(e3); 
    Xyce::Util::newExpression assign_e3; 
    assign_e3 = e3; 

    std::complex<double> result, refRes, copyResult, assignResult;
    std::complex<double> v2 = std::complex<double>(-0.5,1.5);
    ifGroup->setSoln(std::string("v2"),v2);
    e3.evaluateFunction(result);
    copy_e3.evaluateFunction(copyResult);
    assign_e3.evaluateFunction(assignResult);
    refRes = -v2*2.0*3.0;

    EXPECT_EQ( result, refRes);
    EXPECT_EQ( copyResult, refRes);
    EXPECT_EQ( assignResult, refRes);
  }
}

TEST ( ComplexParserIfstatement, min_max)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp_dynamic_cast<ifStatementExpressionGroup>(testGroup);

  // these expressions uses min and max AST nodes
  Xyce::Util::newExpression e4(std::string("max(min(-I(V2), 2.5), 1.5)"), testGroup);
  e4.lexAndParseExpression();

  Xyce::Util::newExpression copy_e4(e4); 
  Xyce::Util::newExpression assign_e4; 
  assign_e4 = e4; 

  int numpoints=100;
  double tmax=1.0;
  double time=0.0;
  double dt=tmax/(numpoints-1);
  std::vector<std::complex<double> > refRes(numpoints), result(numpoints);
  std::vector<std::complex<double> > copyResult(numpoints), assignResult(numpoints);
  for (int ii=0;ii<numpoints;ii++,time+=dt)
  {
    std::complex<double> V0=2, VA=1;
    double FREQ=1, mpi = M_PI;
    std::complex<double> v1 = 0.0;
    if (time <= 0) { v1 = (V0) + (VA) * std::sin (0.0); }
    else { v1 = (V0) + (VA) * std::sin (2.0*mpi*((FREQ)*time + (0.0)/360)) * std::exp( -(time*(0.0))); }
    std::complex<double> v11 = 2.0*v1;
    std::complex<double> v2 = -0.5*v11;

    ifGroup->setTime(time);
    ifGroup->setSoln(std::string("v2"),v2);
    e4.evaluateFunction(result[ii]);
    copy_e4.evaluateFunction(copyResult[ii]);
    assign_e4.evaluateFunction(assignResult[ii]);
    refRes[ii] = std::max(std::min(-v2.real(),2.5),1.5);
  }
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( copyResult, refRes);
  EXPECT_EQ( assignResult, refRes);
}

TEST ( ComplexParserIfstatement, limit)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp_dynamic_cast<ifStatementExpressionGroup>(testGroup);

  // these expressions uses limit AST node
  Xyce::Util::newExpression e5(std::string("limit(-I(V2),1.5,2.5)"), testGroup);
  e5.lexAndParseExpression();

  Xyce::Util::newExpression copy_e5(e5); 
  Xyce::Util::newExpression assign_e5; 
  assign_e5 = e5; 

  int numpoints=100;
  double tmax=1.0;
  double time=0.0;
  double dt=tmax/(numpoints-1);
  std::vector<std::complex<double> > refRes(numpoints), result(numpoints);
  std::vector<std::complex<double> > copyResult(numpoints), assignResult(numpoints);
  for (int ii=0;ii<numpoints;ii++,time+=dt)
  {
    std::complex<double> V0=2, VA=1;
    double FREQ=1, mpi = M_PI;
    std::complex<double> v1 = 0.0;
    if (time <= 0) { v1 = (V0) + (VA) * std::sin (0.0); }
    else { v1 = (V0) + (VA) * std::sin (2.0*mpi*((FREQ)*time + (0.0)/360)) * std::exp( -(time*(0.0))); }
    std::complex<double> v11 = 2.0*v1;
    std::complex<double> v2 = -0.5*v11;

    ifGroup->setTime(time);
    ifGroup->setSoln(std::string("v2"),v2);
    e5.evaluateFunction(result[ii]);
    copy_e5.evaluateFunction(copyResult[ii]);
    assign_e5.evaluateFunction(assignResult[ii]);
    refRes[ii] = std::max(std::min(-v2.real(),2.5),1.5);
  }
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( copyResult, refRes);
  EXPECT_EQ( assignResult, refRes);
}

TEST ( ComplexParserIfstatement, limit2)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  // these expressions uses limit AST node
  Xyce::Util::newExpression e5(std::string("limit(-2.0*I(V2),3.0,5.0)"), baseGroup);
  e5.lexAndParseExpression();

  Xyce::Util::newExpression copy_e5(e5); 
  Xyce::Util::newExpression assign_e5; 
  assign_e5 = e5; 

  std::vector<std::complex<double> > v2Vals(5,0.0);
  v2Vals[0] =  0.0;
  v2Vals[1] = -1.0;
  v2Vals[2] = -2.0;
  v2Vals[3] = -3.0;
  v2Vals[4] = -4.0;

  for(int ii=0;ii<v2Vals.size();++ii)
  {
    std::complex<double> v2 = v2Vals[ii];
    ifGroup->setSoln(std::string("v2"),v2);

    std::complex<double> result;
    std::complex<double> refRes = std::max(std::min((-2.0*std::real(v2)),5.0),3.0);
    std::vector<std::complex<double> > derivs;
    std::vector<std::complex<double> > refDerivs(1,0.0);

    // derivative is the derivative of the first limit argument, 
    // but only if it is between the bounds set by the other 2 args.
    if(ii==2) refDerivs[0] = -2.0;
    else      refDerivs[0] = 0.0;

    e5.evaluate(result,derivs);

    //std::cout << "v2 = " << v2 
      //<< "  refRes = " << refRes << "  refDerivs[0] = " << refDerivs[0] 
      //<< "  result = " << result << "     derivs[0] = " << derivs[0] 
      //<< std::endl;

    ASSERT_EQ( result, refRes);
    ASSERT_EQ( derivs, refDerivs);

    copy_e5.evaluate(result,derivs);
    ASSERT_EQ( result, refRes);
    ASSERT_EQ( derivs, refDerivs);

    assign_e5.evaluate(result,derivs);
    ASSERT_EQ( result, refRes);
    ASSERT_EQ( derivs, refDerivs);
  }

  OUTPUT_MACRO2(ComplexParserIfstatement,limit2, e5) 
}

TEST ( ComplexParserIfstatement, or_true)
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

  std::complex<double> result=0.0;
  e8.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e8.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e8.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(ComplexParserIfstatement, or_true,e8) 
}

TEST ( ComplexParserIfstatement, or_false)
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

  std::complex<double> result=0.0;
  e9.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e9.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e9.evaluateFunction(result); EXPECT_EQ( result, 1.0);
  OUTPUT_MACRO2(ComplexParserIfstatement, or_false,e9) 
}

TEST ( ComplexParserIfstatement, hspice_or_true)
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

  std::complex<double> result=0.0;
  e8.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e8.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e8.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(ComplexParserIfstatement, hspice_or_true,e8) 
}

TEST ( ComplexParserIfstatement, hspice_or_false)
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

  std::complex<double> result=0.0;
  e9.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e9.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e9.evaluateFunction(result); EXPECT_EQ( result, 1.0);
  OUTPUT_MACRO2(ComplexParserIfstatement, hspice_or_false,e9) 
}

TEST ( ComplexParserIfstatement, and_true)
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

  std::complex<double> result=0.0;
  e8.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e8.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e8.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(ComplexParserIfstatement, and_true,e8) 
}

TEST ( ComplexParserIfstatement, and_false)
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

  std::complex<double> result=0.0;
  e9.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e9.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e9.evaluateFunction(result); EXPECT_EQ( result, 1.0);
  OUTPUT_MACRO2(ComplexParserIfstatement, and_false,e9) 
}

TEST ( ComplexParserIfstatement, hspice_and_true)
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

  std::complex<double> result=0.0;
  e8.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e8.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e8.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(ComplexParserIfstatement, hspice_and_true,e8) 
}

TEST ( ComplexParserIfstatement, hspice_and_false)
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

  std::complex<double> result=0.0;
  e9.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e9.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e9.evaluateFunction(result); EXPECT_EQ( result, 1.0);
  OUTPUT_MACRO2(ComplexParserIfstatement, hspice_and_false,e9) 
}

TEST ( ComplexParserIfstatement, xor_true)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp_dynamic_cast<ifStatementExpressionGroup>(testGroup);

  Xyce::Util::newExpression e8(std::string("IF(((V(6) > 1.5) ^ (V(7) < 1.5)), 3, 1)"), testGroup);
  e8.lexAndParseExpression();

  Xyce::Util::newExpression copy_e8(e8); 
  Xyce::Util::newExpression assign_e8; 
  assign_e8 = e8; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),1.0);

  std::complex<double> result=0.0;
  e8.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e8.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e8.evaluateFunction(result); EXPECT_EQ( result, 1.0);
}

TEST ( ComplexParserIfstatement, xor_false)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp_dynamic_cast<ifStatementExpressionGroup>(testGroup);

  Xyce::Util::newExpression e9(std::string("IF(((V(6) > 1.5) ^ (V(7) > 1.5)), 3, 1)"), testGroup);
  e9.lexAndParseExpression();

  Xyce::Util::newExpression copy_e9(e9); 
  Xyce::Util::newExpression assign_e9; 
  assign_e9 = e9; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),1.0);

  std::complex<double> result=0.0;
  e9.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e9.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e9.evaluateFunction(result); EXPECT_EQ( result, 3.0);
}

TEST ( ComplexParserIfstatement, neq)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp_dynamic_cast<ifStatementExpressionGroup>(testGroup);

  Xyce::Util::newExpression e10(std::string("IF((V(6) != V(7)), 3, 1)"), testGroup);
  e10.lexAndParseExpression();

  Xyce::Util::newExpression copy_e10(e10); 
  Xyce::Util::newExpression assign_e10; 
  assign_e10 = e10; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),1.0);

  std::complex<double> result=0.0;
  e10.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e10.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e10.evaluateFunction(result); EXPECT_EQ( result, 3.0);
}

TEST ( ComplexParserIfstatement, not)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp_dynamic_cast<ifStatementExpressionGroup>(testGroup);

  Xyce::Util::newExpression e11(std::string("IF( ~(V(6) > V(7)), 3, 1)"), testGroup);
  e11.lexAndParseExpression();

  Xyce::Util::newExpression copy_e11(e11); 
  Xyce::Util::newExpression assign_e11; 
  assign_e11 = e11; 

  ifGroup->setSoln(std::string("6"),2.0);
  ifGroup->setSoln(std::string("7"),1.0);

  std::complex<double> result=0.0;
  e11.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e11.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e11.evaluateFunction(result); EXPECT_EQ( result, 1.0);
}

TEST ( ComplexParserIfstatement, equiv)
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

  std::complex<double> result=0.0;
  e11.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e11.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e11.evaluateFunction(result); EXPECT_EQ( result, 1.0);
  OUTPUT_MACRO2(ComplexParserIfstatement, equiv, e11) 
}

TEST ( ComplexParserIfstatement, ge1)
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

  std::complex<double> result=0.0;
  e11.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e11.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e11.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(ComplexParserIfstatement, ge1, e11) 
}

TEST ( ComplexParserIfstatement, ge2)
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

  std::complex<double> result=0.0;
  e11.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e11.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e11.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(ComplexParserIfstatement, ge2, e11) 
}

TEST ( ComplexParserIfstatement, ge3)
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

  std::complex<double> result=0.0;
  e11.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e11.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e11.evaluateFunction(result); EXPECT_EQ( result, 1.0);
  OUTPUT_MACRO2(ComplexParserIfstatement, ge3, e11) 
}


TEST ( ComplexParserIfstatement, le1)
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

  std::complex<double> result=0.0;
  e11.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_e11.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_e11.evaluateFunction(result); EXPECT_EQ( result, 1.0);
  OUTPUT_MACRO2(ComplexParserIfstatement, le1, e11) 
}

TEST ( ComplexParserIfstatement, le2)
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

  std::complex<double> result=0.0;
  e11.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e11.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e11.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(ComplexParserIfstatement, le2, e11) 
}

TEST ( ComplexParserIfstatement, le3)
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

  std::complex<double> result=0.0;
  e11.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_e11.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_e11.evaluateFunction(result); EXPECT_EQ( result, 3.0);
  OUTPUT_MACRO2(ComplexParserIfstatement, le3, e11) 
}

// from "ifstatement.cir":
// Also test modulus, along with its precedence vs. + - * and /.
// Because of precedence, this should evaluate to 4.
TEST ( ComplexParserModulus, test1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new testExpressionGroup() );

  // 2 + 6*5/2%4 - 1  = 2 + 30/2%4 -1 = 2 + 15%4 - 1 = 2 + 3 - 1 = 4
  Xyce::Util::newExpression p1(std::string("2 + 6*5/2%4 - 1"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 

  std::complex<double> result=0.0;
  p1.evaluateFunction(result);        EXPECT_EQ( result, 4.0);
  copy_p1.evaluateFunction(result);   EXPECT_EQ( result, 4.0);
  assign_p1.evaluateFunction(result); EXPECT_EQ( result, 4.0);
}

TEST ( ComplexParserModulus, test2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new testExpressionGroup() );

  Xyce::Util::newExpression p1(std::string("15%4"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 

  std::complex<double> result=0.0;
  p1.evaluateFunction(result);        EXPECT_EQ( result, 3.0);
  copy_p1.evaluateFunction(result);   EXPECT_EQ( result, 3.0);
  assign_p1.evaluateFunction(result); EXPECT_EQ( result, 3.0);
}

#if 1
// test for modulus w/real numbers
TEST ( ComplexParserModulus, test3)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new testExpressionGroup() );

  Xyce::Util::newExpression p1(std::string("15.2%4.1"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 

  std::complex<double> result=0.0;
  std::complex<double> refRes=15.2-4.1*std::floor(15.2/4.1);
  p1.evaluateFunction(result);        EXPECT_EQ( result, refRes);
  copy_p1.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assign_p1.evaluateFunction(result); EXPECT_EQ( result, refRes);
  OUTPUT_MACRO2(ComplexParserModulus, test3, p1) 
}

// test for modulus w/real numbers, with problem: second arg is zero.
TEST ( ComplexParserModulus, test4)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new testExpressionGroup() );

  Xyce::Util::newExpression p1(std::string("15.2%0.0"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 

  std::complex<double> result=0.0;
  std::complex<double> refRes=1e50;
  p1.evaluateFunction(result);        EXPECT_EQ( std::fabs(std::real(result)), refRes);
  copy_p1.evaluateFunction(result);   EXPECT_EQ( std::fabs(std::real(result)), refRes);
  assign_p1.evaluateFunction(result); EXPECT_EQ( std::fabs(std::real(result)), refRes);
  OUTPUT_MACRO2(ComplexParserModulus, test4, p1) 
}

// test for modulus w/real numbers, and derivatives
TEST ( ComplexParserModulus, test5)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = solnGroup;

  Xyce::Util::newExpression p1(std::string("V(A)%4.1"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 

  std::complex<double> result=0.0;
  std::complex<double> Aval=std::complex<double> (15.2,0.0);
  std::complex<double> refRes=15.2-4.1*std::floor(15.2/4.1);
  solnGroup->setSoln(std::string("A"),Aval);

  std::vector<std::complex<double> > refDer;
  refDer.push_back(1.0);
  std::vector<std::complex<double> > derivs;

  p1.evaluate(result,derivs);        EXPECT_EQ( result, refRes); EXPECT_EQ(derivs, refDer);
  copy_p1.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ(derivs, refDer);
  assign_p1.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ(derivs, refDer);
  OUTPUT_MACRO2(ComplexParserModulus, test5, p1) 
}

// test for modulus w/real numbers, and derivatives
TEST ( ComplexParserModulus, test6)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = solnGroup;

  Xyce::Util::newExpression p1(std::string("15.2%V(A)"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 


  std::complex<double> result=0.0;
  std::complex<double> Aval=std::complex<double> (4.1,0.0);
  std::complex<double> refRes=15.2-4.1*std::floor(15.2/4.1);
  solnGroup->setSoln(std::string("A"),Aval);

  std::vector<std::complex<double> > refDer;
  refDer.push_back(-std::floor(15.2/4.1));
  std::vector<std::complex<double> > derivs;

  p1.evaluate(result,derivs);        EXPECT_EQ( result, refRes); EXPECT_EQ(derivs, refDer);
  copy_p1.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ(derivs, refDer);
  assign_p1.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ(derivs, refDer);
  OUTPUT_MACRO2(ComplexParserModulus, test6, p1) 
}

// test for modulus w/real numbers, and derivatives
TEST ( ComplexParserModulus, test7)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = solnGroup;

  Xyce::Util::newExpression p1(std::string("V(A)%4.1"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 

  std::complex<double> result=0.0;
  std::complex<double> Aval=std::complex<double> (-15.2,0.0);
  std::complex<double> refRes=-15.2+4.1*std::floor(15.2/4.1);
  solnGroup->setSoln(std::string("A"),Aval);

  std::vector<std::complex<double> > refDer;
  refDer.push_back(1.0);
  std::vector<std::complex<double> > derivs;

  p1.evaluate(result,derivs);        EXPECT_EQ( result, refRes); EXPECT_EQ(derivs, refDer);
  copy_p1.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_EQ(derivs, refDer);
  assign_p1.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_EQ(derivs, refDer);
  OUTPUT_MACRO2(ComplexParserModulus, test7, p1) 
}

// test for modulus w/real numbers, and derivatives
TEST ( ComplexParserModulus, test8)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = solnGroup;

  Xyce::Util::newExpression p1(std::string("15.2%V(A)"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 

  std::complex<double> result=0.0;
  std::complex<double> Aval=std::complex<double> (-4.1,0.0);
  std::complex<double> refRes=  std::complex<double>(std::fmod(15.2,std::real(Aval)),0.0);
  solnGroup->setSoln(std::string("A"),Aval);

  double Aval2=std::real(Aval)*1.1;
  std::complex<double> refResDiff=std::complex<double>(std::fmod(15.2,std::real(Aval2)),0.0);
  std::complex<double> dr= (refRes-refResDiff)/(0.1*Aval); // doing a numerical deriv b/c I keep confusing myself

  std::vector<std::complex<double> > refDer;
  refDer.push_back(dr);
  std::vector<std::complex<double> > derivs;

  p1.evaluate(result,derivs);        EXPECT_EQ( result, refRes); EXPECT_DOUBLE_EQ(std::real(derivs[0]), std::real(refDer[0]));
  copy_p1.evaluate(result,derivs);   EXPECT_EQ( result, refRes); EXPECT_DOUBLE_EQ(std::real(derivs[0]), std::real(refDer[0]));
  assign_p1.evaluate(result,derivs); EXPECT_EQ( result, refRes); EXPECT_DOUBLE_EQ(std::real(derivs[0]), std::real(refDer[0]));
  OUTPUT_MACRO2(ComplexParserModulus, test8, p1) 
}
#endif

TEST ( ComplexParserFmod, fmod1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new testExpressionGroup() );

  Xyce::Util::newExpression p1(std::string("fmod(1.2,0.31)"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 

  std::complex<double> result=0.0;
  std::complex<double> refres=std::fmod(1.2,0.31); // should be 0.27
  p1.evaluateFunction(result);        EXPECT_EQ( result, refres);
  copy_p1.evaluateFunction(result);   EXPECT_EQ( result, refres);
  assign_p1.evaluateFunction(result); EXPECT_EQ( result, refres);
  OUTPUT_MACRO2(ComplexParserFmod, fmod1, p1) 
}

TEST ( ComplexParserFmod, fmod2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new testExpressionGroup() );

  Xyce::Util::newExpression p1(std::string("fmod(5.3,2)"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 

  std::complex<double> result=0.0;
  std::complex<double> refres=std::fmod(5.3,2); // should be 1.3
  p1.evaluateFunction(result);        EXPECT_EQ( result, refres);
  copy_p1.evaluateFunction(result);   EXPECT_EQ( result, refres);
  assign_p1.evaluateFunction(result); EXPECT_EQ( result, refres);
  OUTPUT_MACRO2(ComplexParserFmod, fmod2, p1) 
}

TEST ( ComplexParserFmod, fmod3)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new testExpressionGroup() );

  Xyce::Util::newExpression p1(std::string("fmod(18.5,4.2)"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1);
  Xyce::Util::newExpression assign_p1;
  assign_p1 = p1;

  std::complex<double> result=0.0;
  std::complex<double> refres=std::fmod(18.5,4.2);  // should be 1.7
  p1.evaluateFunction(result);        EXPECT_EQ( result, refres);
  copy_p1.evaluateFunction(result);   EXPECT_EQ( result, refres);
  assign_p1.evaluateFunction(result); EXPECT_EQ( result, refres);
  OUTPUT_MACRO2(ComplexParserFmod, fmod3, p1) 
}

TEST ( ComplexParserNint, nint1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new testExpressionGroup() );

  Xyce::Util::newExpression p1(std::string("nint(2.2)"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 

  std::complex<double> result=0.0;
  std::complex<double> refres=2.0;
  p1.evaluateFunction(result);        EXPECT_EQ( result, refres);
  copy_p1.evaluateFunction(result);   EXPECT_EQ( result, refres);
  assign_p1.evaluateFunction(result); EXPECT_EQ( result, refres);
  OUTPUT_MACRO2(ComplexParserNint, nint1, p1) 
}

TEST ( ComplexParserNint, nint2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new testExpressionGroup() );

  Xyce::Util::newExpression p1(std::string("nint(2.5)"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 

  std::complex<double> result=0.0;
  std::complex<double> refres=3.0;
  p1.evaluateFunction(result);        EXPECT_EQ( result, refres);
  copy_p1.evaluateFunction(result);   EXPECT_EQ( result, refres);
  assign_p1.evaluateFunction(result); EXPECT_EQ( result, refres);
  OUTPUT_MACRO2(ComplexParserNint, nint2, p1) 
}

TEST ( ComplexParserNint, nint3)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new testExpressionGroup() );

  Xyce::Util::newExpression p1(std::string("nint(2.7)"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 

  std::complex<double> result=0.0;
  std::complex<double> refres=3.0;
  p1.evaluateFunction(result);        EXPECT_EQ( result, refres);
  copy_p1.evaluateFunction(result);   EXPECT_EQ( result, refres);
  assign_p1.evaluateFunction(result); EXPECT_EQ( result, refres);
  OUTPUT_MACRO2(ComplexParserNint, nint3, p1) 
}

TEST ( ComplexParserNint, nint4)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new testExpressionGroup() );

  Xyce::Util::newExpression p1(std::string("nint(-2.2)"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 

  std::complex<double> result=0.0;
  std::complex<double> refres=-2.0;
  p1.evaluateFunction(result);        EXPECT_EQ( result, refres);
  copy_p1.evaluateFunction(result);   EXPECT_EQ( result, refres);
  assign_p1.evaluateFunction(result); EXPECT_EQ( result, refres);
  OUTPUT_MACRO2(ComplexParserNint, nint4, p1) 
}

TEST ( ComplexParserNint, nint5)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new testExpressionGroup() );

  Xyce::Util::newExpression p1(std::string("nint(-2.5)"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 

  std::complex<double> result=0.0;
  std::complex<double> refres=-3.0;
  p1.evaluateFunction(result);        EXPECT_EQ( result, refres);
  copy_p1.evaluateFunction(result);   EXPECT_EQ( result, refres);
  assign_p1.evaluateFunction(result); EXPECT_EQ( result, refres);
  OUTPUT_MACRO2(ComplexParserNint, nint5, p1) 
}

TEST ( ComplexParserNint, nint6)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new testExpressionGroup() );

  Xyce::Util::newExpression p1(std::string("nint(-2.7)"), grp);
  p1.lexAndParseExpression();

  Xyce::Util::newExpression copy_p1(p1); 
  Xyce::Util::newExpression assign_p1; 
  assign_p1 = p1; 

  std::complex<double> result=0.0;
  std::complex<double> refres=-3.0;
  p1.evaluateFunction(result);        EXPECT_EQ( result, refres);
  copy_p1.evaluateFunction(result);   EXPECT_EQ( result, refres);
  assign_p1.evaluateFunction(result); EXPECT_EQ( result, refres);
  OUTPUT_MACRO2(ComplexParserNint, nint6, p1) 
}

//-------------------------------------------------------------------------------
// table tests
//
// adapted from break.cir
TEST ( ComplexParserTableTest, break1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<timeDepExpressionGroup> timeGroup = Teuchos::rcp_dynamic_cast<timeDepExpressionGroup>(grp);
  Xyce::Util::newExpression tableExpression(std::string("Table(time, 0, 0, 0.3, 0, 0.301, 2, 0.302, 2, 0.6, 1, 1, 1)"), grp);
  tableExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_tableExpression(tableExpression); 
  Xyce::Util::newExpression assign_tableExpression; 
  assign_tableExpression = tableExpression; 

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<std::complex<double> > refRes = { 0, 0, 2, 2, 1, 1 };
  std::vector<std::complex<double> > result(times.size(),0.0);
  std::vector<std::complex<double> > copyResult(times.size(),0.0);
  std::vector<std::complex<double> > assignResult(times.size(),0.0);

  for (int ii=0;ii<times.size();ii++)  
  { 
    timeGroup->setTime(times[ii]); 
    tableExpression.evaluateFunction(result[ii]); 
    copy_tableExpression.evaluateFunction(copyResult[ii]); 
    assign_tableExpression.evaluateFunction(assignResult[ii]); 
  }
  EXPECT_EQ(refRes,result);
  EXPECT_EQ(refRes,copyResult);
  EXPECT_EQ(refRes,assignResult);
}

TEST ( ComplexParserTableTest, tablefile_break1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<timeDepExpressionGroup> timeGroup = Teuchos::rcp_dynamic_cast<timeDepExpressionGroup>(grp);
  Xyce::Util::newExpression tableExpression(std::string("tablefile(\"test1.dat\")"), grp);
  tableExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_tableExpression(tableExpression); 
  Xyce::Util::newExpression assign_tableExpression; 
  assign_tableExpression = tableExpression; 

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<std::complex<double> > refRes = { 0, 0, 2, 2, 1, 1 };
  std::vector<std::complex<double> > result(times.size(),0.0);
  std::vector<std::complex<double> > copyResult(times.size(),0.0);
  std::vector<std::complex<double> > assignResult(times.size(),0.0);

  for (int ii=0;ii<times.size();ii++)  
  { 
    timeGroup->setTime(times[ii]); 
    tableExpression.evaluateFunction(result[ii]); 
    copy_tableExpression.evaluateFunction(copyResult[ii]); 
    assign_tableExpression.evaluateFunction(assignResult[ii]); 
  }
  EXPECT_EQ(refRes,result);
  EXPECT_EQ(refRes,copyResult);
  EXPECT_EQ(refRes,assignResult);
}

// This test is for :
// (1) making sure filenames can include dashes 
// (2) making sure table files can include comments.  
// The file, test1-1.dat includes comments in it.
TEST ( ComplexParserTableTest, tablefile_dash_comment)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<timeDepExpressionGroup> timeGroup = Teuchos::rcp_dynamic_cast<timeDepExpressionGroup>(grp);
  Xyce::Util::newExpression tableExpression(std::string("tablefile(\"test1-1.dat\")"), grp);
  tableExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_tableExpression(tableExpression); 
  Xyce::Util::newExpression assign_tableExpression; 
  assign_tableExpression = tableExpression; 

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<std::complex<double> > refRes = { 0, 0, 2, 2, 1, 1 };
  std::vector<std::complex<double> > result(times.size(),0.0);
  std::vector<std::complex<double> > copyResult(times.size(),0.0);
  std::vector<std::complex<double> > assignResult(times.size(),0.0);

  for (int ii=0;ii<times.size();ii++)  
  { 
    timeGroup->setTime(times[ii]); 
    tableExpression.evaluateFunction(result[ii]); 
    copy_tableExpression.evaluateFunction(copyResult[ii]); 
    assign_tableExpression.evaluateFunction(assignResult[ii]); 
  }
  EXPECT_EQ(refRes,result);
  EXPECT_EQ(refRes,copyResult);
  EXPECT_EQ(refRes,assignResult);
}

TEST ( ComplexParserTableTest, tablefile_break1b)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = timeDepGroup;
  Xyce::Util::newExpression tableExpression(std::string("tablefile(\"./SubDir1/test1.dat\")"), grp);
  tableExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_tableExpression(tableExpression); 
  Xyce::Util::newExpression assign_tableExpression; 
  assign_tableExpression = tableExpression; 

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<std::complex<double> > refRes = { 0, 0, 2, 2, 1, 1 };
  std::vector<std::complex<double> > result(times.size(),0.0);
  std::vector<std::complex<double> > copyResult(times.size(),0.0);
  std::vector<std::complex<double> > assignResult(times.size(),0.0);

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
  OUTPUT_MACRO2(ComplexParserTableTest, tablefile_break1b, tableExpression) 
}

TEST ( ComplexParserTableTest, tablefile_break1c)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = timeDepGroup;
  Xyce::Util::newExpression tableExpression(std::string("tablefile(\"./Sub_Dir/1test_5.dat\")"), grp);
  tableExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_tableExpression(tableExpression); 
  Xyce::Util::newExpression assign_tableExpression; 
  assign_tableExpression = tableExpression; 

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<std::complex<double> > refRes = { 0, 0, 2, 2, 1, 1 };
  std::vector<std::complex<double> > result(times.size(),0.0);
  std::vector<std::complex<double> > copyResult(times.size(),0.0);
  std::vector<std::complex<double> > assignResult(times.size(),0.0);

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
  OUTPUT_MACRO2(ComplexParserTableTest, tablefile_break1c, tableExpression) 
}

TEST ( ComplexParserTableTest, break2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<timeDepExpressionGroup> timeGroup = Teuchos::rcp_dynamic_cast<timeDepExpressionGroup>(grp);
  Xyce::Util::newExpression tableExpression(std::string("Table({time} 0, 0, 0.3, 0, 0.301, 2, 0.302, 2, 0.6, 1, 1, 1)"), grp);
  tableExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_tableExpression(tableExpression); 
  Xyce::Util::newExpression assign_tableExpression; 
  assign_tableExpression = tableExpression; 

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<std::complex<double> > refRes = { 0, 0, 2, 2, 1, 1 };
  std::vector<std::complex<double> > result(times.size(),0.0);
  std::vector<std::complex<double> > copyResult(times.size(),0.0);
  std::vector<std::complex<double> > assignResult(times.size(),0.0);

  for (int ii=0;ii<times.size();ii++) 
  { 
    timeGroup->setTime(times[ii]); 
    tableExpression.evaluateFunction(result[ii]); 
    copy_tableExpression.evaluateFunction(copyResult[ii]); 
    assign_tableExpression.evaluateFunction(assignResult[ii]); 
  }
  EXPECT_EQ(refRes,result);
  EXPECT_EQ(refRes,copyResult);
  EXPECT_EQ(refRes,assignResult);
}

//-------------------------------------------------------------------------------
// fasttable tests
//
// adapted from break.cir
TEST ( ComplexParserfasttableTest, break1)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = timeDepGroup;
  Xyce::Util::newExpression tableExpression(std::string("fastTable(time, 0, 0, 0.3, 0, 0.301, 2, 0.302, 2, 0.6, 1, 1, 1)"), grp);
  tableExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_tableExpression(tableExpression); 
  Xyce::Util::newExpression assign_tableExpression; 
  assign_tableExpression = tableExpression; 

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<std::complex<double> > refRes = { 0, 0, 2, 2, 1, 1 };
  std::vector<std::complex<double> > result(times.size(),0.0);
  std::vector<std::complex<double> > copyResult(times.size(),0.0);
  std::vector<std::complex<double> > assignResult(times.size(),0.0);

  // don't compare everything for now, just check parsing.
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

  OUTPUT_MACRO2(ComplexParserTableTest, break1, tableExpression) 
}

// adapted from power_thermalres_gear.cir
TEST ( ComplexParserTableTest, power_thermalres)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = Teuchos::rcp(new tempDepExpressionGroup() );
  Teuchos::RCP<tempDepExpressionGroup> tempGroup = Teuchos::rcp_dynamic_cast<tempDepExpressionGroup>(grp);
  {
    Xyce::Util::newExpression resistivity(std::string("table(temp+273.15, 0, 0.5e-9, 100, 3e-9, 1000, 6.6e-8)"), grp);
    resistivity.lexAndParseExpression();

    Xyce::Util::newExpression copy_resistivity(resistivity); 
    Xyce::Util::newExpression assign_resistivity; 
    assign_resistivity = resistivity; 

    std::vector<double> temps = { 0-273.15, 100-273.15, 1000-273.15 };
    std::vector<std::complex<double> > refRes = { 0.5e-9, 3e-9, 6.6e-8 };
    std::vector<std::complex<double> > result(temps.size(),0.0);
    std::vector<std::complex<double> > copyResult(temps.size(),0.0);
    std::vector<std::complex<double> > assignResult(temps.size(),0.0);

    for (int ii=0;ii<temps.size();ii++) 
    { 
      tempGroup->setTemp(temps[ii]); 
      resistivity.evaluateFunction(result[ii]); 
      copy_resistivity.evaluateFunction(copyResult[ii]); 
      assign_resistivity.evaluateFunction(assignResult[ii]); 
    }
    EXPECT_EQ(refRes,result);
    EXPECT_EQ(refRes,copyResult);
    EXPECT_EQ(refRes,assignResult);
  }
  {
    Xyce::Util::newExpression heatcapacity(std::string("8.92e+3*table(temp+273.15, 0, 1, 1000, 1500)"), grp);
    heatcapacity.lexAndParseExpression();
    
    Xyce::Util::newExpression copy_heatcapacity(heatcapacity); 
    Xyce::Util::newExpression assign_heatcapacity; 
    assign_heatcapacity = heatcapacity; 

    std::vector<double> temps = { 0-273.15, 1000-273.15 };
    std::vector<std::complex<double> > refRes = { 8.92e+3, 8.92e+3*1500 };
    std::vector<std::complex<double> > result(temps.size(),0.0);
    std::vector<std::complex<double> > copyResult(temps.size(),0.0);
    std::vector<std::complex<double> > assignResult(temps.size(),0.0);

    for (int ii=0;ii<temps.size();ii++) 
    { 
      tempGroup->setTemp(temps[ii]); 
      heatcapacity.evaluateFunction(result[ii]); 
      copy_heatcapacity.evaluateFunction(copyResult[ii]); 
      assign_heatcapacity.evaluateFunction(assignResult[ii]); 
    }
    EXPECT_EQ(refRes,result);
    EXPECT_EQ(refRes,copyResult);
    EXPECT_EQ(refRes,assignResult);
  }
}

#if 0
// adapted from Bsrc_C1.cir.
// See the ComplexParserParam_Test.test2  test below,
// which also tests the first-argument expression and uses
// some of the same machinery.
//
// the mechanics of the table seem to work here, but I haven't generated a good gold standard yet
// so, disabling for now
TEST ( ComplexParserTableTest, Bsrc_C1_withoutParens)
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

TEST ( ComplexParserTableTest, Bsrc_C1_pureArray)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  grp = Teuchos::rcp(new Bsrc_C1_ExpressionGroup() );
  Teuchos::RCP<Bsrc_C1_ExpressionGroup> bsrc_C1_grp = Teuchos::rcp_dynamic_cast<Bsrc_C1_ExpressionGroup>(grp);
  {
    // this is a nice test b/c it has curly braces around the first expression, which is one of the supported formats
    Xyce::Util::newExpression BE_Dig(std::string("TABLE({ V(2) * (V(1) + 30) / 60 } 0.0000000, 0,0.0312500, 0,0.0312813, 1,0.0625000, 1,0.0625313, 2,0.0937500, 2,0.0937813, 3,0.1250000, 3,0.1250313, 4,0.1562500, 4,0.1562813, 5,0.1875000, 5,0.1875313, 6,0.2187500, 6,0.2187813, 7,0.2500000, 7,0.2500313, 8,0.2812500, 8,0.2812813, 9,0.3125000, 9,0.3125313, 10,0.3437500, 10,0.3437813, 11,0.3750000, 11,0.3750313, 12,0.4062500, 12,0.4062813, 13,0.4375000, 13,0.4375313, 14,0.4687500, 14,0.4687813, 15,0.5000000, 15,0.5000313, 16,0.5312500, 16,0.5312813, 17,0.5625000, 17,0.5625313, 18,0.5937500, 18,0.5937813, 19,0.6250000, 19,0.6250313, 20,0.6562500, 20,0.6562813, 21,0.6875000, 21,0.6875313, 22,0.7187500, 22,0.7187813, 23,0.7500000, 23,0.7500313, 24,0.7812500, 24,0.7812813, 25,0.8125000, 25,0.8125313, 26,0.8437500, 26,0.8437813, 27,0.8750000, 27,0.8750313, 28,0.9062500, 28,0.9062813, 29,0.9375000, 29,0.9375313, 30,0.9687500, 30,0.9687813, 31,1.0000000, 31)"), grp);
    BE_Dig.lexAndParseExpression();

    Xyce::Util::newExpression copy_BE_Dig(BE_Dig); 
    Xyce::Util::newExpression assign_BE_Dig; 
    assign_BE_Dig = BE_Dig; 

    Xyce::Util::newExpression BE_Dig_leftArg(std::string("V(2) * (V(1) + 30) / 60"),grp);
    BE_Dig_leftArg.lexAndParseExpression();

    std::vector<std::complex<double> > xa = { 0.0000000, 0.0312500, 0.0312813, 0.0625000, 0.0625313, 0.0937500, 0.0937813, 0.1250000, 0.1250313, 0.1562500, 0.1562813, 0.1875000, 0.1875313, 0.2187500, 0.2187813, 0.2500000, 0.2500313, 0.2812500, 0.2812813, 0.3125000, 0.3125313, 0.3437500, 0.3437813, 0.3750000, 0.3750313, 0.4062500, 0.4062813, 0.4375000, 0.4375313, 0.4687500, 0.4687813, 0.5000000, 0.5000313, 0.5312500, 0.5312813, 0.5625000, 0.5625313, 0.5937500, 0.5937813, 0.6250000, 0.6250313, 0.6562500, 0.6562813, 0.6875000, 0.6875313, 0.7187500, 0.7187813, 0.7500000, 0.7500313, 0.7812500, 0.7812813, 0.8125000, 0.8125313, 0.8437500, 0.8437813, 0.8750000, 0.8750313, 0.9062500, 0.9062813, 0.9375000, 0.9375313, 0.9687500, 0.9687813, 1.0000000 };

    std::vector<std::complex<double> > ya = { 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31 };

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
    std::vector<std::complex<double> > refRes(numpoints), result(numpoints);
    std::vector<std::complex<double> > copyResult(numpoints), assignResult(numpoints);
    for (int ii=0;ii<numpoints;ii++,time+=dt)
    {
      bsrc_C1_grp->setTime(time);
      std::complex<double> v1Value(0.0),v2Value(0.0);
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

TEST ( ComplexParserTableTest, Bsrc_C1_pairsWithParens)
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

    std::vector<std::complex<double> > xa = { 0.0000000, 0.0312500, 0.0312813, 0.0625000, 0.0625313, 0.0937500, 0.0937813, 0.1250000, 0.1250313, 0.1562500, 0.1562813, 0.1875000, 0.1875313, 0.2187500, 0.2187813, 0.2500000, 0.2500313, 0.2812500, 0.2812813, 0.3125000, 0.3125313, 0.3437500, 0.3437813, 0.3750000, 0.3750313, 0.4062500, 0.4062813, 0.4375000, 0.4375313, 0.4687500, 0.4687813, 0.5000000, 0.5000313, 0.5312500, 0.5312813, 0.5625000, 0.5625313, 0.5937500, 0.5937813, 0.6250000, 0.6250313, 0.6562500, 0.6562813, 0.6875000, 0.6875313, 0.7187500, 0.7187813, 0.7500000, 0.7500313, 0.7812500, 0.7812813, 0.8125000, 0.8125313, 0.8437500, 0.8437813, 0.8750000, 0.8750313, 0.9062500, 0.9062813, 0.9375000, 0.9375313, 0.9687500, 0.9687813, 1.0000000 };

    std::vector<std::complex<double> > ya = { 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31 };

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
    std::vector<std::complex<double> > refRes(numpoints), result(numpoints);
    std::vector<std::complex<double> > copyResult(numpoints), assignResult(numpoints);
    for (int ii=0;ii<numpoints;ii++,time+=dt)
    {
      bsrc_C1_grp->setTime(time);
      std::complex<double> v1Value(0.0),v2Value(0.0);
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


TEST ( ComplexParserTableTest, power_e_gear)
{
  Teuchos::RCP<Bsrc_C1_ExpressionGroup> bsrc_C1_grp = Teuchos::rcp(new Bsrc_C1_ExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  grp = bsrc_C1_grp;
  {
    Xyce::Util::newExpression eTable(std::string("TABLE {V(1,0)} = ( 0 , 1 ) ( 1 , 2 )"), grp);
    eTable.lexAndParseExpression();

    Xyce::Util::newExpression eTable_leftArg(std::string("V(1,0)"),grp);
    eTable_leftArg.lexAndParseExpression();

    // v1:
    std::vector<std::complex<double> > v1 = { -5.00000000e-01, -4.00000000e-01, -3.00000000e-01, -2.00000000e-01, -1.00000000e-01, 0.00000000e+00, 1.00000000e-01, 2.00000000e-01, 3.00000000e-01, 4.00000000e-01, 5.00000000e-01, 6.00000000e-01, 7.00000000e-01, 8.00000000e-01, 9.00000000e-01, 1.00000000e+00, 1.10000000e+00, 1.20000000e+00, 1.30000000e+00, 1.40000000e+00, 1.50000000e+00};

    // table output
    std::vector<std::complex<double> > refArray = {1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.10000000e+00, 1.20000000e+00, 1.30000000e+00, 1.40000000e+00, 1.50000000e+00, 1.60000000e+00, 1.70000000e+00, 1.80000000e+00, 1.90000000e+00, 2.00000000e+00, 2.00000000e+00, 2.00000000e+00, 2.00000000e+00, 2.00000000e+00, 2.00000000e+00 };

    int size = v1.size();

    std::vector<std::complex<double> > result;
    result.resize(size,0.0);

    for (int ii=0;ii<size;ii++)
    {
      bsrc_C1_grp->setSoln(std::string("1"),v1[ii]);
      eTable.evaluateFunction(result[ii]);
    }

    EXPECT_EQ(refArray,result);
  }
}

TEST ( ComplexParserTableTest, power_endcomma)
{
  Teuchos::RCP<Bsrc_C1_ExpressionGroup> bsrc_C1_grp = Teuchos::rcp(new Bsrc_C1_ExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  grp = bsrc_C1_grp;
  {
    Xyce::Util::newExpression eTable(std::string("TABLE(V(1,0),0,1,1,2,)"), grp);
    eTable.lexAndParseExpression();

    Xyce::Util::newExpression eTable_leftArg(std::string("V(1,0)"),grp);
    eTable_leftArg.lexAndParseExpression();

    // v1:
    std::vector<std::complex<double> > v1 = { -5.00000000e-01, -4.00000000e-01, -3.00000000e-01, -2.00000000e-01, -1.00000000e-01, 0.00000000e+00, 1.00000000e-01, 2.00000000e-01, 3.00000000e-01, 4.00000000e-01, 5.00000000e-01, 6.00000000e-01, 7.00000000e-01, 8.00000000e-01, 9.00000000e-01, 1.00000000e+00, 1.10000000e+00, 1.20000000e+00, 1.30000000e+00, 1.40000000e+00, 1.50000000e+00};

    // table output
    std::vector<std::complex<double> > refArray = {1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.10000000e+00, 1.20000000e+00, 1.30000000e+00, 1.40000000e+00, 1.50000000e+00, 1.60000000e+00, 1.70000000e+00, 1.80000000e+00, 1.90000000e+00, 2.00000000e+00, 2.00000000e+00, 2.00000000e+00, 2.00000000e+00, 2.00000000e+00, 2.00000000e+00 };

    int size = v1.size();

    std::vector<std::complex<double> > result;
    result.resize(size,0.0);

    for (int ii=0;ii<size;ii++)
    {
      bsrc_C1_grp->setSoln(std::string("1"),v1[ii]);
      eTable.evaluateFunction(result[ii]);
    }

    EXPECT_EQ(refArray,result);
  }
}

//-------------------------------------------------------------------------------
// spline tests
//
// adapted from break.cir
TEST ( ComplexParsersplineTest, break1)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = timeDepGroup;
  Xyce::Util::newExpression splineExpression(std::string("Spline(time, 0, 0, 0.3, 0, 0.301, 2, 0.302, 2, 0.6, 1, 1, 1)"), grp);
  splineExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_splineExpression(splineExpression); 
  Xyce::Util::newExpression assign_splineExpression; 
  assign_splineExpression = splineExpression; 

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<std::complex<double> > refRes = { 0, 0, 2, 2, 1, 1 };
  std::vector<std::complex<double> > result(times.size(),0.0);
  std::vector<std::complex<double> > copyResult(times.size(),0.0);
  std::vector<std::complex<double> > assignResult(times.size(),0.0);

  for (int ii=0;ii<times.size();ii++) 
  { 
    timeDepGroup->setTime(times[ii]); 
    splineExpression.evaluateFunction(result[ii]); 
    copy_splineExpression.evaluateFunction(copyResult[ii]); 
    assign_splineExpression.evaluateFunction(assignResult[ii]); 

    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(result[ii]));
    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(copyResult[ii]));
    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(assignResult[ii]));

  }
  OUTPUT_MACRO2(ComplexParsersplineTest, break1, splineExpression) 
}

TEST ( ComplexParsersplineTest, break2)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = timeDepGroup;
  Xyce::Util::newExpression splineExpression(std::string("Spline({time} 0, 0, 0.3, 0, 0.301, 2, 0.302, 2, 0.6, 1, 1, 1)"), grp);
  splineExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_splineExpression(splineExpression); 
  Xyce::Util::newExpression assign_splineExpression; 
  assign_splineExpression = splineExpression; 

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<std::complex<double> > refRes = { 0, 0, 2, 2, 1, 1 };
  std::vector<std::complex<double> > result(times.size(),0.0);
  std::vector<std::complex<double> > copyResult(times.size(),0.0);
  std::vector<std::complex<double> > assignResult(times.size(),0.0);

  for (int ii=0;ii<times.size();ii++) 
  { 
    timeDepGroup->setTime(times[ii]); 
    splineExpression.evaluateFunction(result[ii]); 
    copy_splineExpression.evaluateFunction(copyResult[ii]); 
    assign_splineExpression.evaluateFunction(assignResult[ii]); 

    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(result[ii]));
    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(copyResult[ii]));
    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(assignResult[ii]));
  }
  OUTPUT_MACRO2(ComplexParsersplineTest, break2, splineExpression) 
}

// adapted from power_thermalres_gear.cir
TEST ( ComplexParsersplineTest, power_thermalres)
{
  Teuchos::RCP<tempDepExpressionGroup> tempDepGroup = Teuchos::rcp(new tempDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  grp = tempDepGroup;
  {
    Xyce::Util::newExpression resistivity(std::string("spline(temp+273.15, 0, 0.5e-9, 100, 3e-9, 1000, 6.6e-8)"), grp);
    resistivity.lexAndParseExpression();

    Xyce::Util::newExpression copy_resistivity(resistivity); 
    Xyce::Util::newExpression assign_resistivity; 
    assign_resistivity = resistivity; 

    std::vector<double> temps = { 0-273.15, 100-273.15, 1000-273.15 };
    std::vector<std::complex<double> > refRes = { 0.5e-9, 3e-9, 6.6e-8 };
    std::vector<std::complex<double> > result(temps.size(),0.0);
    std::vector<std::complex<double> > copyResult(temps.size(),0.0);
    std::vector<std::complex<double> > assignResult(temps.size(),0.0);

    for (int ii=0;ii<temps.size();ii++) 
    { 
      tempDepGroup->setTemp(temps[ii]); 
      resistivity.evaluateFunction(result[ii]); 
      copy_resistivity.evaluateFunction(copyResult[ii]); 
      assign_resistivity.evaluateFunction(assignResult[ii]); 

      EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(result[ii]));
      EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(copyResult[ii]));
      EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(assignResult[ii]));
    }
    OUTPUT_MACRO2(ComplexParsersplineTest, power_thermalres, resistivity) 
  }
  {
    Xyce::Util::newExpression heatcapacity(std::string("8.92e+3*spline(temp+273.15, 0, 1, 1000, 1500)"), grp);
    heatcapacity.lexAndParseExpression();
    
    Xyce::Util::newExpression copy_heatcapacity(heatcapacity); 
    Xyce::Util::newExpression assign_heatcapacity; 
    assign_heatcapacity = heatcapacity; 

    std::vector<double> temps = { 0-273.15, 1000-273.15 };
    std::vector<std::complex<double> > refRes = { 8.92e+3, 8.92e+3*1500 };
    std::vector<std::complex<double> > result(temps.size(),0.0);
    std::vector<std::complex<double> > copyResult(temps.size(),0.0);
    std::vector<std::complex<double> > assignResult(temps.size(),0.0);

    for (int ii=0;ii<temps.size();ii++) 
    { 
      tempDepGroup->setTemp(temps[ii]); 
      heatcapacity.evaluateFunction(result[ii]); 
      copy_heatcapacity.evaluateFunction(copyResult[ii]); 
      assign_heatcapacity.evaluateFunction(assignResult[ii]); 

      EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(result[ii]));
      EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(copyResult[ii]));
      EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(assignResult[ii]));
    }
    OUTPUT_MACRO2(ComplexParsersplineTest, power_thermalres, heatcapacity) 
  }
}

#if 1
TEST ( ComplexParsersplineTest, akima1)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = timeDepGroup;
  Xyce::Util::newExpression splineExpression(std::string("akima (time, 0, 0, 0.3, 0, 0.301, 2, 0.302, 2, 0.6, 1, 1, 1)"), grp);
  splineExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_splineExpression(splineExpression); 
  Xyce::Util::newExpression assign_splineExpression; 
  assign_splineExpression = splineExpression; 

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<std::complex<double> > refRes = { 0, 0, 2, 2, 1, 1 };
  std::vector<std::complex<double> > result(times.size(),0.0);
  std::vector<std::complex<double> > copyResult(times.size(),0.0);
  std::vector<std::complex<double> > assignResult(times.size(),0.0);

  for (int ii=0;ii<times.size();ii++) 
  { 
    timeDepGroup->setTime(times[ii]); 
    splineExpression.evaluateFunction(result[ii]); 
    copy_splineExpression.evaluateFunction(copyResult[ii]); 
    assign_splineExpression.evaluateFunction(assignResult[ii]); 

    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(result[ii]));
    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(copyResult[ii]));
    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(assignResult[ii]));
  }
  OUTPUT_MACRO2(ComplexParsersplineTest, akima1, splineExpression) 
}

TEST ( ComplexParsersplineTest, cubic1)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = timeDepGroup;
  Xyce::Util::newExpression splineExpression(std::string("cubic (time, 0, 0, 0.3, 0, 0.301, 2, 0.302, 2, 0.6, 1, 1, 1)"), grp);
  splineExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_splineExpression(splineExpression); 
  Xyce::Util::newExpression assign_splineExpression; 
  assign_splineExpression = splineExpression; 

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<std::complex<double> > refRes = { 0, 0, 2, 2, 1, 1 };
  std::vector<std::complex<double> > result(times.size(),0.0);
  std::vector<std::complex<double> > copyResult(times.size(),0.0);
  std::vector<std::complex<double> > assignResult(times.size(),0.0);

  for (int ii=0;ii<times.size();ii++) 
  { 
    timeDepGroup->setTime(times[ii]); 
    splineExpression.evaluateFunction(result[ii]); 
    copy_splineExpression.evaluateFunction(copyResult[ii]); 
    assign_splineExpression.evaluateFunction(assignResult[ii]); 

    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(result[ii]));
    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(copyResult[ii]));
    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(assignResult[ii]));
  }
  OUTPUT_MACRO2(ComplexParsersplineTest, cubic1, splineExpression) 
}

TEST ( ComplexParsersplineTest, barycentric1)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = timeDepGroup;
  Xyce::Util::newExpression splineExpression(std::string("bli (time, 0, 0, 0.3, 0, 0.301, 2, 0.302, 2, 0.6, 1, 1, 1)"), grp);
  splineExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_splineExpression(splineExpression); 
  Xyce::Util::newExpression assign_splineExpression; 
  assign_splineExpression = splineExpression; 

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<std::complex<double> > refRes = { 0, 0, 2, 2, 1, 1 };
  std::vector<std::complex<double> > result(times.size(),0.0);
  std::vector<std::complex<double> > copyResult(times.size(),0.0);
  std::vector<std::complex<double> > assignResult(times.size(),0.0);

  for (int ii=0;ii<times.size();ii++) 
  { 
    timeDepGroup->setTime(times[ii]); 
    splineExpression.evaluateFunction(result[ii]); 
    copy_splineExpression.evaluateFunction(copyResult[ii]); 
    assign_splineExpression.evaluateFunction(assignResult[ii]); 

    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(result[ii]));
    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(copyResult[ii]));
    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(assignResult[ii]));
  }
  OUTPUT_MACRO2(ComplexParsersplineTest, barycentric1, splineExpression) 
}

TEST ( ComplexParsersplineTest, wodicka1)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = timeDepGroup;
  Xyce::Util::newExpression splineExpression(std::string("wodicka (time, 0, 0, 0.3, 0, 0.301, 2, 0.302, 2, 0.6, 1, 1, 1)"), grp);
  splineExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_splineExpression(splineExpression); 
  Xyce::Util::newExpression assign_splineExpression; 
  assign_splineExpression = splineExpression; 

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<std::complex<double> > refRes = { 0, 0, 2, 2, 1, 1 };
  std::vector<std::complex<double> > result(times.size(),0.0);
  std::vector<std::complex<double> > copyResult(times.size(),0.0);
  std::vector<std::complex<double> > assignResult(times.size(),0.0);

  for (int ii=0;ii<times.size();ii++) 
  { 
    timeDepGroup->setTime(times[ii]); 
    splineExpression.evaluateFunction(result[ii]); 
    copy_splineExpression.evaluateFunction(copyResult[ii]); 
    assign_splineExpression.evaluateFunction(assignResult[ii]); 

    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(result[ii]));
    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(copyResult[ii]));
    EXPECT_FLOAT_EQ(std::real(refRes[ii]),std::real(assignResult[ii]));
  }
  OUTPUT_MACRO2(ComplexParsersplineTest, wodicka1, splineExpression) 
}
#endif

// Schedule is like a table, but with no interpolation.
// If the schedule is (t0, dt0, t1, dt1, t2, dt2) 
// then the value is: 
// if time < t0      value = 0 
// if t0 < time < t1 value = dt0 
// if t1 < time < t2 value = dt1 
// if t2 < time      value = dt2
TEST ( ComplexParserScheduleTest, test1)
{
  Teuchos::RCP<timeDepExpressionGroup> timeGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  grp = timeGroup;
  {
    Xyce::Util::newExpression eSchedule(std::string("schedule( 0.5e-3, 0, 1.0e-3, 1.0e-6, 1.5e-3, 1.0e-4, 2.0e-3, 0 )"), grp);
    eSchedule.lexAndParseExpression();

    std::complex<double> result = 0.0;
    std::complex<double> reference = 1.0e-6;
    timeGroup->setTime(1.2e-3);
    eSchedule.evaluateFunction(result);

    EXPECT_EQ(reference,result);
  }
}



//-------------------------------------------------------------------------------
// .param tests
//-------------------------------------------------------------------------------
// this form of test1 doesn't rely on the group to resolve the parameter.
// Instead, it allows the user to attach it.
//-------------------------------------------------------------------------------
TEST ( ComplexParserParamTest, test1)
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

  std::complex<double> result;
  testExpression.evaluateFunction(result);        EXPECT_EQ( result, 5.0 );
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( result, 5.0 );
  OUTPUT_MACRO(ComplexParserParamTest, test1)
}

// this form of test1 doesn't rely on the group to resolve the parameter.
// Instead, it allows the user to attach it.
TEST ( ComplexParserParamTest, testE1)
{
  Teuchos::RCP<testExpressionGroup> noparamGroup = Teuchos::rcp(new testExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = noparamGroup;

  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("2+3"), testGroup));
  p1Expression->lexAndParseExpression();
  std::string p1Name = "E1";

  Xyce::Util::newExpression testExpression(std::string("E1"), testGroup);
  testExpression.lexAndParseExpression();
  testExpression.attachParameterNode(p1Name,p1Expression);

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  //testExpression.dumpParseTree(std::cout);

  std::complex<double>  result;
  testExpression.evaluateFunction(result);        EXPECT_EQ( result, 5.0 );
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( result, 5.0 );
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( result, 5.0 );
  OUTPUT_MACRO(ComplexParserParamTest, testE1)
}

TEST ( ComplexParserParamTest, test2)
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
  testExpression.attachParameterNode(p1Name,p1exp);

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  int numpoints=100;
  double tfinal = 0.0005;
  double dt = tfinal/(numpoints-1), time=0.0;

  std::vector<std::complex<double> > refRes(numpoints), result(numpoints);
  std::vector<std::complex<double> > copyResult(numpoints), assignResult(numpoints);
  for (int ii=0;ii<numpoints;ii++,time+=dt)
  {
    paramGroup->setTime(time);
    std::complex<double>  v1Value = std::complex<double>(0.0,0.0);
    std::complex<double>  v2Value = std::complex<double>(0.0,0.0);
    v1exp.evaluateFunction(v1Value);
    v2exp.evaluateFunction(v2Value);
    paramGroup->setSoln(std::string("1"),v1Value);
    paramGroup->setSoln(std::string("2"),v2Value);
    testExpression.evaluateFunction(result[ii]);
    copy_testExpression.evaluateFunction(copyResult[ii]);
    assign_testExpression.evaluateFunction(assignResult[ii]);
    refRes[ii] = 2.0 * v2Value * (v1Value + 30.0) / 60.0;

    ASSERT_NEAR(std::real(refRes[ii]),std::real(result[ii]),1.0e-15);
    ASSERT_NEAR(std::real(refRes[ii]),std::real(copyResult[ii]),1.0e-15);
    ASSERT_NEAR(std::real(refRes[ii]),std::real(assignResult[ii]),1.0e-15);
  }
  OUTPUT_MACRO(ComplexParserParamTest, test2)
}

TEST ( ComplexParserCalculus, ddx1)
{
  Teuchos::RCP<testExpressionGroupWithParamSupport> paramGroup = Teuchos::rcp(new testExpressionGroupWithParamSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;
  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2+3"), testGroup));
  p1Expression->lexAndParseExpression();
  std::string p1Name="p1";

  Xyce::Util::newExpression ddxTest(std::string("ddx(2*p1,p1)"), testGroup); ddxTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  std::complex<double> result;
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result, 2.0 );
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result, 2.0 );
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result, 2.0 );
}

TEST ( ComplexParserCalculus, ddx2)
{
  Teuchos::RCP<testExpressionGroupWithParamSupport> paramGroup = Teuchos::rcp(new testExpressionGroupWithParamSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;
  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("2+3"), testGroup));
  p1Expression->lexAndParseExpression();
  std::string p1Name="p1";

  Xyce::Util::newExpression ddxTest(std::string("ddx(p1*p1,p1)"), testGroup); ddxTest.lexAndParseExpression();
  ddxTest.attachParameterNode(p1Name,p1Expression);

  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  std::complex<double> result;
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result, 2.0*(2+3) );
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result, 2.0*(2+3) );
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result, 2.0*(2+3) );
}

TEST ( ComplexParserCalculus, ddx3)
{
  Teuchos::RCP<testExpressionGroupWithParamSupport> paramGroup = Teuchos::rcp(new testExpressionGroupWithParamSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;
  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2+3"), testGroup));
  p1Expression->lexAndParseExpression();
  std::string p1Name="p1";

  Xyce::Util::newExpression ddxTest(std::string("ddx(sin(p1),p1)"), testGroup); ddxTest.lexAndParseExpression();
  ddxTest.attachParameterNode(p1Name,p1Expression);

  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  std::complex<double> result;
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result, std::cos(2+3) );
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result, std::cos(2+3) );
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result, std::cos(2+3) );
}

TEST ( ComplexParserCalculus, ddx4)
{
  Teuchos::RCP<testExpressionGroupWithParamSupport> paramGroup = Teuchos::rcp(new testExpressionGroupWithParamSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;
  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2+3"), testGroup)); 
  p1Expression->lexAndParseExpression();
  std::string p1Name="p1";

  Xyce::Util::newExpression ddxTest(std::string("ddx(sin(p1*p1),p1)"), testGroup); ddxTest.lexAndParseExpression();
  ddxTest.attachParameterNode(p1Name,p1Expression);

  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  std::complex<double> result;
  std::complex<double> p1 = 2+3;
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result, 2.0*p1*std::cos(p1*p1) );
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result, 2.0*p1*std::cos(p1*p1) );
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result, 2.0*p1*std::cos(p1*p1) );
}

TEST ( ComplexParserCalculus, ddx5)
{
  Teuchos::RCP<testExpressionGroupWithParamSupport> paramGroup = Teuchos::rcp(new testExpressionGroupWithParamSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;
  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2+3"), testGroup)); 
  p1Expression->lexAndParseExpression();
  std::string p1Name="p1";

  Xyce::Util::newExpression ddxTest(std::string("ddx( pow(sin(p1*p1),3.0),p1)"), testGroup); ddxTest.lexAndParseExpression();

  ddxTest.attachParameterNode(p1Name,p1Expression);

  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  std::complex<double> result;
  std::complex<double> p1 = 2+3;
  std::complex<double> p1Sq = p1*p1;
  std::complex<double> refRes = 3.0*(2.0*p1*std::cos(p1Sq))/(std::sin(p1Sq))*std::pow(std::sin(p1Sq),3.0);

  ddxTest.evaluateFunction(result);
  EXPECT_NEAR( std::real(result-refRes), 0.0, 1.0e-15 );
  EXPECT_NEAR( std::imag(result-refRes), 0.0, 1.0e-15 );

  copy_ddxTest.evaluateFunction(result);
  EXPECT_EQ((result-refRes), 0.0);

  assign_ddxTest.evaluateFunction(result);
  EXPECT_EQ((result-refRes), 0.0);
}

TEST ( ComplexParserCalculus, ddx5b)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression ddxTest(std::string("pow(sin(V(A)*V(A)),3.0)"), testGroup); ddxTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  std::complex<double> Aval=5.0;
  solnGroup->setSoln(std::string("A"),Aval);
  std::complex<double> result;
  std::vector<std::complex<double> > derivs;
  std::complex<double> p1 = 5;
  std::complex<double> p1Sq = p1*p1;
  std::complex<double> refRes = 3.0*(2.0*p1*std::cos(p1Sq))/(std::sin(p1Sq))*std::pow(std::sin(p1Sq),3.0);
  std::vector<std::complex<double> > refderivs = { refRes };

  ddxTest.evaluate(result,derivs);
  EXPECT_EQ( derivs.size(), refderivs.size() );
  for(int i=0; i < derivs.size(); ++i) {
    EXPECT_NEAR(std::real(derivs[i]), std::real(refderivs[i]), 1.0e-15);
    EXPECT_NEAR(std::imag(derivs[i]), std::imag(refderivs[i]), 1.0e-15);
  }

  copy_ddxTest.evaluate(result,derivs);
  EXPECT_EQ( derivs, refderivs );

  assign_ddxTest.evaluate(result,derivs);
  EXPECT_EQ( derivs, refderivs );
}

TEST ( ComplexParserCalculus, ddx6)
{
  Teuchos::RCP<testExpressionGroup> paramGroup = Teuchos::rcp(new testExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;
  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2.0+3.0"), testGroup)); 
  p1Expression->lexAndParseExpression();
  std::string p1Name="p1";

  Xyce::Util::newExpression ddxTest(std::string("ddx( pow(p1,3.0),p1)"), testGroup); ddxTest.lexAndParseExpression();
  ddxTest.attachParameterNode(p1Name,p1Expression);

  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  std::complex<double>  result;
  std::complex<double>  p1 = 2.0+3.0;
  std::complex<double>  exponent = 3.0;
  std::complex<double>  refRes = exponent*p1*p1; // doesn't work very well.
  std::complex<double>  refRes2 = (exponent/p1)*std::pow(p1,exponent); // exactly what the expression library does; works much better

  ddxTest.evaluateFunction(result);        ASSERT_NEAR( std::real(result), std::real(refRes2), 1.0e-13);
  copy_ddxTest.evaluateFunction(result);   ASSERT_NEAR( std::real(result), std::real(refRes2), 1.0e-13);
  assign_ddxTest.evaluateFunction(result); ASSERT_NEAR( std::real(result), std::real(refRes2), 1.0e-13);
} 

TEST ( ComplexParserCalculus, ddx7)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression ddxTest(std::string("ddx(sin(v(a)),v(a))"), testGroup); ddxTest.lexAndParseExpression();
  
  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  solnGroup->setSoln(std::string("A"),5.0);
  std::complex<double> result;
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result, std::cos(5.0) );
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result, std::cos(5.0) );
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result, std::cos(5.0) );
}

// this test will produce an error message.  That is (as of this writing) the correct behavior.
// I am leaving it in the unit test suite to prove that it doesn't produce a core dump (it used to do this)
TEST ( ComplexParserCalculus, ddx8)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression ddxTest(std::string("ddx(sin(v(a,b)),v(a,b))"), testGroup); ddxTest.lexAndParseExpression();
  
  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  solnGroup->setSoln(std::string("A"),5.0);
  solnGroup->setSoln(std::string("B"),3.0);
  std::complex<double> result;
  std::complex<double> refRes=0.0;
  //std::complex<double> refRes=std::cos(2.0); // this would be correct if the type was supported
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result, refRes);
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result, refRes);
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result, refRes);
}

TEST ( ComplexParserCalculus, ddx9)
{
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression ddxTest(std::string("ddx(sin(i(vb)),i(vb))"), testGroup); ddxTest.lexAndParseExpression();
  
  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  solnGroup->setSoln(std::string("VB"),5.0);
  std::complex<double> result;
  ddxTest.evaluateFunction(result);        EXPECT_EQ( result, std::cos(5.0) );
  copy_ddxTest.evaluateFunction(result);   EXPECT_EQ( result, std::cos(5.0) );
  assign_ddxTest.evaluateFunction(result); EXPECT_EQ( result, std::cos(5.0) );
}

TEST ( ComplexParserCalculus, ddx10)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression ddxTest(std::string("ddx(pow(5.0,v(a)),v(a))"), testGroup); ddxTest.lexAndParseExpression();
  
  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  std::complex<double> Aval=2.0;
  solnGroup->setSoln(std::string("A"),Aval);
  std::complex<double> result;

  //std::complex<double>  refRes = std::complex<double>(std::log(5.0)*std::pow(5.0,2.0),0.0); // doesn't quite work
  std::complex<double>  refRes = std::log(5.0)*std::pow(5.0,Aval); // exactly what the expression library does; works much better

  ddxTest.evaluateFunction(result);
  ASSERT_NEAR(std::real(result - refRes), 0.0, 1.0e-14);
  ASSERT_NEAR(std::imag(result - refRes), 0.0, 1.0e-14);

  copy_ddxTest.evaluateFunction(result);
  ASSERT_EQ(result, refRes);

  assign_ddxTest.evaluateFunction(result);
  ASSERT_EQ(result, refRes);
}

TEST ( ComplexParserCalculus, ddx11)
{
  Teuchos::RCP<testExpressionGroupWithParamSupport> paramGroup = Teuchos::rcp(new testExpressionGroupWithParamSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;

  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2"), testGroup)); 
  p1Expression->lexAndParseExpression();
  std::string p1Name="p1";
  //paramGroup->addParam(p1Name,p1Expression);

  Xyce::Util::newExpression ddxTest(std::string("ddx( pow(sin(p1),p1),p1)"), testGroup); ddxTest.lexAndParseExpression();
  //ddxTest.resolveExpression();
  ddxTest.attachParameterNode(p1Name,p1Expression);
  
  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  std::complex<double> Aval=2.0;
  std::complex<double> result;
  ddxTest.evaluateFunction(result);
  std::complex<double> acceptedVal(std::pow(std::sin(Aval),Aval)*(Aval/std::tan(Aval) + std::log(sin(Aval))));   
  EXPECT_DOUBLE_EQ( result.real(), acceptedVal.real() );
  EXPECT_DOUBLE_EQ( result.imag(), acceptedVal.imag() );
  copy_ddxTest.evaluateFunction(result);
  EXPECT_DOUBLE_EQ( result.real(), acceptedVal.real() );
  EXPECT_DOUBLE_EQ( result.imag(), acceptedVal.imag() );
  //EXPECT_DOUBLE_EQ( result, std::pow(std::sin(Aval),Aval)*(Aval/std::tan(Aval) + std::log(sin(Aval))) );
  assign_ddxTest.evaluateFunction(result);
  EXPECT_DOUBLE_EQ( result.real(), acceptedVal.real() );
  EXPECT_DOUBLE_EQ( result.imag(), acceptedVal.imag() );
  //EXPECT_DOUBLE_EQ( result, std::pow(std::sin(Aval),Aval)*(Aval/std::tan(Aval) + std::log(sin(Aval))) );
}

//-------------------------------------------------------------------------------
// using ddx with a table thru a func.  
//
// When this test was created (8/5/2020) the expression that uses ddx just 
// returned a zero. ie, it didn't work.  Now it is fixed.  The problem was 
// with the dx function in the tableOp class.
//
// ddx seems to work fine thru .funcs, even though I don't seem to have any 
// unit tests for it other than this one.
//-------------------------------------------------------------------------------
TEST ( ComplexParserCalculus, ddx12)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;

  Teuchos::RCP<Xyce::Util::newExpression> poffExpression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("0.005"), testGroup)); 
  poffExpression->lexAndParseExpression();
  std::string poffName="Poff";

  Teuchos::RCP<Xyce::Util::newExpression> zspaceExpression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("0.0707"), testGroup)); 
  zspaceExpression->lexAndParseExpression();
  std::string zspaceName="Zspace";

  //.FUNC LMz1_3(x) {TABLE(x, ... )}
  Teuchos::RCP<Xyce::Util::newExpression> lmz1_3_Expression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("TABLE(x, -1.562200e-001,8.314111e-007, -1.504200e-001,8.855486e-007, -1.446200e-001,9.441836e-007, -1.388200e-001,1.007741e-006, -1.330200e-001,1.076705e-006, -1.272200e-001,1.151629e-006, -1.214200e-001,1.233096e-006, -1.156200e-001,1.321782e-006, -1.098200e-001,1.418414e-006, -1.040200e-001,1.523815e-006, -9.822000e-002,1.638897e-006, -9.242000e-002,1.764608e-006, -8.662000e-002,1.902105e-006, -8.082000e-002,2.052570e-006, -7.502000e-002,2.217342e-006, -6.922000e-002,2.397851e-006, -6.342000e-002,2.595707e-006, -5.762000e-002,2.812610e-006, -5.182000e-002,3.050351e-006, -4.602000e-002,3.310836e-006, -4.022000e-002,3.596021e-006, -3.442000e-002,3.907863e-006, -2.862000e-002,4.248147e-006, -2.282000e-002,4.618467e-006, -1.702000e-002,5.019898e-006, -1.122000e-002,5.452745e-006, -5.420000e-003,5.916007e-006, 3.800000e-004,6.407217e-006, 6.180000e-003,6.920922e-006, 1.198000e-002,7.449045e-006, 1.778000e-002,7.979794e-006, 2.358000e-002,8.497570e-006, 2.938000e-002,8.984073e-006, 3.518000e-002,9.419870e-006, 4.098000e-002,9.788100e-006, 4.678000e-002,1.007599e-005, 5.258000e-002,1.027653e-005, 5.838000e-002,1.038781e-005, 6.418000e-002,1.041027e-005, 6.998000e-002,1.034424e-005, 7.578000e-002,1.018858e-005, 8.158000e-002,9.944564e-006, 8.738000e-002,9.616062e-006, 9.318000e-002,9.212598e-006, 9.898000e-002,8.749304e-006, 1.047800e-001,8.244964e-006, 1.105800e-001,7.718794e-006, 1.163800e-001,7.187442e-006, 1.221800e-001,6.665222e-006, 1.279800e-001,6.161786e-006, 1.337800e-001,5.683901e-006, 1.395800e-001,5.235392e-006, 1.453800e-001,4.818026e-006, 1.511800e-001,4.431996e-006)"), testGroup)); 

  Xyce::Util::newExpression lmz1_3_LHS (std::string("LMz1_3(x)"), testGroup);
  lmz1_3_LHS.lexAndParseExpression();
  std::vector<std::string> lmz1_3_ArgStrings ;
  lmz1_3_LHS.getFuncPrototypeArgStrings(lmz1_3_ArgStrings);

  lmz1_3_Expression->setFunctionArgStringVec ( lmz1_3_ArgStrings );
  lmz1_3_Expression->lexAndParseExpression();

  std::string lmz1_3_Name;
  lmz1_3_LHS.getFuncPrototypeName(lmz1_3_Name);

   //.FUNC dLMzdz1_3(x) TABLE(x,
  Teuchos::RCP<Xyce::Util::newExpression> dlmzdz1_3_Expression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("TABLE(x, -0.15622, 0, -0.15332, 9.33405e-06, -0.14752, 1.01095e-05, -0.14172, 1.09582e-05, -0.13592, 1.18903e-05, -0.13012, 1.29179e-05, -0.12432, 1.4046e-05, -0.11852, 1.52907e-05, -0.11272, 1.66607e-05, -0.10692, 1.81726e-05, -0.10112, 1.98417e-05, -0.09532, 2.16743e-05, -0.08952, 2.37064e-05, -0.08372, 2.59422e-05, -0.07792, 2.8409e-05, -0.07212, 3.11222e-05, -0.06632, 3.41131e-05, -0.06052, 3.73971e-05, -0.05472, 4.09898e-05, -0.04892, 4.49112e-05, -0.04312, 4.91698e-05, -0.03732, 5.37659e-05, -0.03152, 5.86697e-05, -0.02572, 6.38483e-05, -0.01992, 6.92122e-05, -0.01412, 7.46288e-05, -0.00832, 7.98728e-05, -0.00252, 8.46914e-05, 0.00328, 8.85698e-05, 0.00908, 9.10557e-05, 0.01488, 9.15084e-05, 0.02068, 8.92717e-05, 0.02648, 8.38798e-05, 0.03228, 7.51374e-05, 0.03808, 6.34879e-05, 0.04388, 4.96362e-05, 0.04968, 3.45759e-05, 0.05548, 1.91862e-05, 0.06128, 3.87241e-06, 0.06708, -1.13845e-05, 0.07288, -2.68379e-05, 0.07868, -4.20717e-05, 0.08448, -5.66383e-05, 0.09028, -6.95628e-05, 0.09608, -7.98783e-05, 0.10188, -8.69552e-05, 0.10768, -9.0719e-05, 0.11348, -9.16124e-05, 0.11928, -9.00379e-05, 0.12508, -8.67993e-05, 0.13088, -8.2394e-05, 0.13668, -7.73291e-05, 0.14248, -7.19597e-05, 0.14828, -6.65569e-05, 0.15118, 0)"), testGroup)); 

  Xyce::Util::newExpression dlmzdz1_3_LHS (std::string("dLMzdz1_3(x)"), testGroup);
  dlmzdz1_3_LHS.lexAndParseExpression();
  std::vector<std::string> dlmzdz1_3_ArgStrings ;
  dlmzdz1_3_LHS.getFuncPrototypeArgStrings(dlmzdz1_3_ArgStrings);

  dlmzdz1_3_Expression->setFunctionArgStringVec ( dlmzdz1_3_ArgStrings );
  dlmzdz1_3_Expression->lexAndParseExpression();

  std::string dlmzdz1_3_Name;
  dlmzdz1_3_LHS.getFuncPrototypeName(dlmzdz1_3_Name);

  Xyce::Util::newExpression ddxTest(std::string("ddx(LMz1_3(v(b)+Poff-Zspace),v(b))"), testGroup); 
  ddxTest.lexAndParseExpression();
  ddxTest.attachFunctionNode(lmz1_3_Name, lmz1_3_Expression);
  ddxTest.attachParameterNode(poffName,poffExpression);
  ddxTest.attachParameterNode(zspaceName,zspaceExpression);
  
  Xyce::Util::newExpression baseline(std::string("dLMzdz1_3(v(b)+Poff-Zspace)"), testGroup); 
  baseline.lexAndParseExpression();
  baseline.attachFunctionNode(dlmzdz1_3_Name, dlmzdz1_3_Expression);
  baseline.attachParameterNode(poffName,poffExpression);
  baseline.attachParameterNode(zspaceName,zspaceExpression);

std::vector<double> xvals = 
  {-0.15622, -0.15332, -0.14752, -0.14172, -0.13592, -0.13012, -0.12432, -0.11852,
   -0.11272, -0.10692, -0.10112, -0.09532, -0.08952, -0.08372, -0.07792, -0.07212,
   -0.06632, -0.06052, -0.05472, -0.04892, -0.04312, -0.03732, -0.03152, -0.02572,
   -0.01992, -0.01412, -0.00832, -0.00252, 0.00328, 0.00908, 0.01488, 0.02068,
    0.02648, 0.03228, 0.03808, 0.04388, 0.04968, 0.05548, 0.06128, 0.06708,
    0.07288, 0.07868, 0.08448, 0.09028, 0.09608, 0.10188, 0.10768, 0.11348,
    0.11928, 0.12508, 0.13088, 0.13668, 0.14248, 0.14828, 0.15118 };

  int size = xvals.size();
  for(int ii=0;ii<size;ii++)
  {
    solnGroup->setSoln(std::string("b"), xvals[ii]);
    std::complex<double> result1;
    std::complex<double> result2;
    ddxTest.evaluateFunction(result1);        
    baseline.evaluateFunction(result2);        
    ASSERT_NEAR(std::real(result1), std::real(result2), std::abs(1.0e-4*std::real(result1)));
    ASSERT_NEAR(std::imag(result1), std::imag(result2), std::abs(1.0e-4*std::imag(result1)));
  }
}

TEST ( ComplexParserCalculus, simpleDerivs1 )
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression ddxTest(std::string("0.5*(V(B)-2.0-1.2J)**2.0"), testGroup); 
  ddxTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  std::complex<double> Aval=std::complex<double>(2.5,1.0);
  solnGroup->setSoln(std::string("B"),Aval);
  std::complex<double> result;
  std::vector<std::complex<double> > derivs;

  std::complex<double> diff = Aval-std::complex<double>(2.0,1.2);
  std::complex<double> refRes = 0.5*std::pow(diff,2.0);
  std::vector<std::complex<double> > refderivs = { 0.5*2.0/diff*std::pow(diff,2.0) };

  ddxTest.evaluate(result,derivs);        
  EXPECT_DOUBLE_EQ(std::real(result), std::real(refRes));
  EXPECT_DOUBLE_EQ(std::imag(result), std::imag(refRes));

  EXPECT_EQ(derivs.size(), refderivs.size());
  for(int i=0; i < derivs.size(); ++i) {
    EXPECT_DOUBLE_EQ(std::real(derivs[i]), std::real(refderivs[i]));
    EXPECT_DOUBLE_EQ(std::imag(derivs[i]), std::imag(refderivs[i]));
  }

  copy_ddxTest.evaluate(refRes,refderivs);   
  EXPECT_EQ(result,refRes);
  EXPECT_EQ(derivs, refderivs);

  assign_ddxTest.evaluate(refRes,refderivs); 
  EXPECT_EQ(result,refRes);
  EXPECT_EQ(derivs, refderivs);
}

#if 1
//-------------------------------------------------------------------------------
TEST ( ComplexParserCalculus, simpleDerivs2 )
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression ddxTest(std::string("1.5/(1.0+2.0*V(A)*V(B))"), testGroup); 
  ddxTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_ddxTest(ddxTest); 
  Xyce::Util::newExpression assign_ddxTest; 
  assign_ddxTest = ddxTest; 

  std::complex<double> Aval=2.0;
  std::complex<double> Bval=0.5;
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("B"),Bval);
  std::complex<double> result;
  std::vector< std::complex<double> > derivs;
  std::complex<double> refRes = 0.5;
  std::vector<std::complex<double> > refderivs = { -(1.5/9.0), -(6.0/9.0) };

  ddxTest.evaluate(result,derivs);        EXPECT_EQ( derivs, refderivs );
  copy_ddxTest.evaluate(result,derivs);   EXPECT_EQ( derivs, refderivs );
  assign_ddxTest.evaluate(result,derivs); EXPECT_EQ( derivs, refderivs );
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserCalculus, ddx13)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression ddxTest1(std::string("ddx(1.5/(1.0+2.0*V(A)*V(B)),V(A))"), testGroup); 
  Xyce::Util::newExpression ddxTest2(std::string("ddx(1.5/(1.0+2.0*V(A)*V(B)),V(B))"), testGroup); 
  ddxTest1.lexAndParseExpression();
  ddxTest2.lexAndParseExpression();

  std::complex<double> Aval=2.0;
  std::complex<double> Bval=0.5;
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("B"),Bval);
  std::complex<double> result;
  std::complex<double> refRes = -1.5/9.0;
  ddxTest1.evaluateFunction(result);    EXPECT_EQ( result, refRes );

  refRes = -2.0/3.0;
  ddxTest2.evaluateFunction(result);    EXPECT_EQ( result, refRes );
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserCalculus, ddx14)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression ddxTest1(std::string("ddx(1.5/(1.0+2.0*I(VA)*I(VB)),I(VA))"), testGroup); 
  Xyce::Util::newExpression ddxTest2(std::string("ddx(1.5/(1.0+2.0*I(VA)*I(VB)),I(VB))"), testGroup); 
  ddxTest1.lexAndParseExpression();
  ddxTest2.lexAndParseExpression();

  std::complex<double> Aval=2.0;
  std::complex<double> Bval=0.5;
  solnGroup->setSoln(std::string("VA"),Aval);
  solnGroup->setSoln(std::string("VB"),Bval);
  std::complex<double> result;
  std::complex<double> refRes = -1.5/9.0;
  ddxTest1.evaluateFunction(result);    EXPECT_EQ( result, refRes );

  refRes = -2.0/3.0;
  ddxTest2.evaluateFunction(result);    EXPECT_EQ( result, refRes );
}
#endif

//-------------------------------------------------------------------------------
// These tests (derivsThruFuncs?) tests if derivatives work thru expression arguments.
// At the time of test creation (2/21/2020), the answer was NO.
//
//-------------------------------------------------------------------------------
TEST ( ComplexParserCalculus, derivsThruFuncs1 )
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

  std::string f1Name;
  f1_LHS.getFuncPrototypeName(f1Name);

  Xyce::Util::newExpression derivFuncTestExpr(std::string("0.5*(F1(V(B),(2.0+0.8J)))**2.0"), testGroup); 
  derivFuncTestExpr.lexAndParseExpression();
  derivFuncTestExpr.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copy_derivFuncTestExpr(derivFuncTestExpr); 
  Xyce::Util::newExpression assign_derivFuncTestExpr; 
  assign_derivFuncTestExpr = derivFuncTestExpr; 

  std::complex<double> Bval=std::complex<double>(2.5,1.2);
  solnFuncGroup->setSoln(std::string("B"),Bval);

  std::complex<double> result;
  std::vector<std::complex<double> > derivs;

  std::complex<double> diff = Bval - std::complex<double>(2.0,0.8);

  std::complex<double> refRes = 0.5*std::pow(diff,2.0);
  std::vector<std::complex<double> > refderivs = { 0.5*2.0/diff*std::pow(diff,2.0) };

  derivFuncTestExpr.evaluate(result,derivs);        
  EXPECT_EQ( derivs, refderivs );

  copy_derivFuncTestExpr.evaluate(result,derivs);   
  EXPECT_EQ( derivs, refderivs );

  assign_derivFuncTestExpr.evaluate(result,derivs); 
  EXPECT_EQ( derivs, refderivs );
}

// there are a bunch of tests in this one, all testing if derivatives propagate thru
// a .FUNC in the right way.
TEST ( ComplexParserCalculus, derivsThruFuncs2 )
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
  Teuchos::RCP<Xyce::Util::newExpression> f2Expression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("A*B"), testGroup));
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

  std::complex<double>  Aval=std::complex<double>(0.45,0.0);
  std::complex<double>  Bval=std::complex<double>(0.6,0.0);
  solnFuncGroup->setSoln(std::string("A"),Aval);
  solnFuncGroup->setSoln(std::string("B"),Bval);
  std::complex<double>  result;
  std::complex<double>  result2;
  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > derivs2;

#if 0
  // analytic answer: seems to have a roundoff problem compared to computed result.
  // I can make the expression library match analytic result, but only for 
  // certain Aval, Bval values.
  std::complex<double>  refRes;
  std::vector<std::complex<double> > refderivs;
  {
    std::complex<double>  f1val = std::sin(Aval)*std::cos(Bval);
    refRes =  0.5*f1val*f1val;
    std::complex<double>  df1_dA = +std::cos(Aval)*std::cos(Bval);
    std::complex<double>  df1_dB = -std::sin(Aval)*std::sin(Bval);
    std::complex<double>  dExp_df1 = f1val;
    std::complex<double>  dExp_dA = df1_dA*f1val;
    std::complex<double>  dExp_dB = df1_dB*f1val;
    refderivs = { dExp_dA, dExp_dB };
  }
#endif

  derivFuncTestExpr1.evaluate(result,derivs);
  derivFuncTestExpr2.evaluate(result2,derivs2);
  EXPECT_EQ( result, result2 );
  EXPECT_EQ( derivs, derivs2 );

#if 0
  EXPECT_EQ( result-refRes, 0.0 );

  std::vector<std::complex<double> > derivDiffs = { (derivs[0]-refderivs[0]),  (derivs[1]-refderivs[1]) };
  EXPECT_EQ( derivDiffs[0], 0.0 );
  EXPECT_EQ( derivDiffs[1], 0.0 );
#endif

  copy_derivFuncTestExpr1.evaluate(result,derivs);   
  EXPECT_EQ( result, result2 );
  EXPECT_EQ( derivs, derivs2 );

  assign_derivFuncTestExpr1.evaluate(result,derivs); 
  EXPECT_EQ( result, result2 );
  EXPECT_EQ( derivs, derivs2 );
}

TEST ( ComplexParserCalculus, derivsThruFuncs3 )
{
  Teuchos::RCP<solnAndFuncExpressionGroup> solnFuncGroup = Teuchos::rcp(new solnAndFuncExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnFuncGroup;

  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression (std::string("A*B"), testGroup));
  std::vector<std::string> f1ArgStrings = { std::string("A"), std::string("B") };
  f1Expression->setFunctionArgStringVec ( f1ArgStrings );
  f1Expression->lexAndParseExpression();

  //solnFuncGroup->addFunction( std::string("F1"),  f1Expression);

  Xyce::Util::newExpression derivFuncTestExpr1(std::string("V(A)*F1(V(A),V(B)*V(B))+3.0*V(B)"), testGroup); 

  derivFuncTestExpr1.lexAndParseExpression();
  //derivFuncTestExpr1.resolveExpression(); 
  derivFuncTestExpr1.attachFunctionNode(std::string("F1"), f1Expression);

  Xyce::Util::newExpression copy_derivFuncTestExpr1(derivFuncTestExpr1); 
  Xyce::Util::newExpression assign_derivFuncTestExpr1; 
  assign_derivFuncTestExpr1 = derivFuncTestExpr1; 

  std::complex<double>  Aval=std::complex<double>(1.0,3.0);
  std::complex<double>  Bval=std::complex<double>(2.0,4.0);
  solnFuncGroup->setSoln(std::string("A"),Aval);
  solnFuncGroup->setSoln(std::string("B"),Bval);
  std::complex<double>  result;
  std::vector<std::complex<double> > derivs;

  derivFuncTestExpr1.evaluate(result,derivs);   

  std::complex<double>  resRef = Aval*Aval*Bval*Bval+3.0*Bval;
  std::complex<double>  dfdA = (2.0*Aval*Bval*Bval);
  std::complex<double>  dfdB = (2.0*Aval*Aval*Bval + 3.0);

  std::vector<std::complex<double> > derivsRef = { dfdA, dfdB };

  EXPECT_EQ( result-resRef, 0.0 );

  std::vector<std::complex<double> > derivDiffs = { (derivs[0]-derivsRef[0]),  (derivs[1]-derivsRef[1]) };
  EXPECT_EQ( derivDiffs[0], 0.0 );
  EXPECT_EQ( derivDiffs[1], 0.0 );

}

TEST ( ComplexParserCalculus, derivsThruFuncs4 )
{
  Teuchos::RCP<solnAndFuncExpressionGroup> solnFuncGroup = Teuchos::rcp(new solnAndFuncExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnFuncGroup;

  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("A*B+100*V(C)"), testGroup));
  std::vector<std::string> f1ArgStrings = { std::string("A"), std::string("B") };
  f1Expression->setFunctionArgStringVec ( f1ArgStrings );
  f1Expression->lexAndParseExpression();

  //solnFuncGroup->addFunction( std::string("F1"),  f1Expression);

  Xyce::Util::newExpression derivFuncTestExpr1(std::string("F1(V(A),V(B))"), testGroup); 

  derivFuncTestExpr1.lexAndParseExpression();
  //derivFuncTestExpr1.resolveExpression(); 
  derivFuncTestExpr1.attachFunctionNode(std::string("F1"), f1Expression);

  Xyce::Util::newExpression copy_derivFuncTestExpr1(derivFuncTestExpr1); 
  Xyce::Util::newExpression assign_derivFuncTestExpr1; 
  assign_derivFuncTestExpr1 = derivFuncTestExpr1; 

  std::complex<double> Aval=std::complex<double>(17.0,4.0);
  std::complex<double> Bval=std::complex<double>(2.0,1.0);
  std::complex<double> Cval=std::complex<double>(3.0,5.0);
  solnFuncGroup->setSoln(std::string("A"),Aval);
  solnFuncGroup->setSoln(std::string("B"),Bval);
  solnFuncGroup->setSoln(std::string("C"),Cval);
  std::complex<double> result;
  std::vector<std::complex<double> > derivs;

  derivFuncTestExpr1.evaluate(result,derivs);   

  std::complex<double> resRef = (Aval*Bval+100.0*Cval);
  std::vector<std::complex<double> > derivsRef = { Bval, Aval, 100.0 };

  EXPECT_EQ( result,resRef);
  EXPECT_EQ( derivs,derivsRef);
}

TEST ( ComplexParserCalculus, derivsThruFuncs5 )
{
  Teuchos::RCP<solnAndFuncExpressionGroup> solnFuncGroup = Teuchos::rcp(new solnAndFuncExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnFuncGroup;

  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression(std::string("A*B+100*V(A)"), testGroup));
  std::vector<std::string> f1ArgStrings = { std::string("A"), std::string("B") };
  f1Expression->setFunctionArgStringVec ( f1ArgStrings );
  f1Expression->lexAndParseExpression();

  //solnFuncGroup->addFunction( std::string("F1"),  f1Expression);

  Xyce::Util::newExpression derivFuncTestExpr1(std::string("F1(V(A),V(B))"), testGroup); 

  derivFuncTestExpr1.lexAndParseExpression();
  //derivFuncTestExpr1.resolveExpression(); 
  derivFuncTestExpr1.attachFunctionNode(std::string("F1"), f1Expression);

  //derivFuncTestExpr1.dumpParseTree(std::cout);

  Xyce::Util::newExpression copy_derivFuncTestExpr1(derivFuncTestExpr1); 
  Xyce::Util::newExpression assign_derivFuncTestExpr1; 
  assign_derivFuncTestExpr1 = derivFuncTestExpr1; 

  std::complex<double>  Aval=std::complex<double>(17.0,4.0);
  std::complex<double>  Bval=std::complex<double>(2.0,1.0);
  solnFuncGroup->setSoln(std::string("A"),Aval);
  solnFuncGroup->setSoln(std::string("B"),Bval);
  std::complex<double>  result;
  std::vector<std::complex<double> > derivs;

  derivFuncTestExpr1.evaluate(result,derivs);   

  std::complex<double>  resRef = (Aval*Bval+100.0*Aval);
  std::vector<std::complex<double> > derivsRef = { (Bval+100.0), Aval };

  EXPECT_EQ( result,resRef);
  EXPECT_EQ( derivs,derivsRef);
}

//
TEST ( ComplexParserCalculus, derivsThruFuncs6 )
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
  f1Expression->lexAndParseExpression();

  std::string f1Name;
  f1_LHS.getFuncPrototypeName(f1Name);

  Xyce::Util::newExpression derivFuncTestExpr(std::string("DDX(F1(V(B),3.0),V(B))"), testGroup); 
  derivFuncTestExpr.lexAndParseExpression();
  derivFuncTestExpr.attachFunctionNode(f1Name, f1Expression);

  //Xyce::Util::newExpression copy_derivFuncTestExpr(derivFuncTestExpr); 
  //Xyce::Util::newExpression assign_derivFuncTestExpr; 
  //assign_derivFuncTestExpr = derivFuncTestExpr; 

  double Bval=2.5;
  solnFuncGroup->setSoln(std::string("B"),Bval);
  std::complex<double> result;
  std::complex<double> refRes = 1.0;

  //derivFuncTestExpr.dumpParseTree(std::cout);

  derivFuncTestExpr.evaluateFunction(result);       EXPECT_EQ(result, refRes);

  //copy_derivFuncTestExpr.evaluate(result,derivs);   EXPECT_EQ( derivs, refderivs );
  //assign_derivFuncTestExpr.evaluate(result,derivs); EXPECT_EQ( derivs, refderivs );
}

//
TEST ( ComplexParserCalculus, derivsThruFuncs7 )
{
  Teuchos::RCP<solnAndFuncExpressionGroup> solnFuncGroup = Teuchos::rcp(new solnAndFuncExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnFuncGroup;

  // this expression is the RHS of a .global_param statement:  .global_param P1 {V(B)}
  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("V(B)"), testGroup));
  p1Expression->lexAndParseExpression();
  std::string p1Name = "P1";

  // this expression is the RHS of a .func statement:  .func F1(A,B) {A-B}
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  = Teuchos::rcp(new Xyce::Util::newExpression (std::string("A-B"), testGroup));
  std::vector<std::string> f1ArgStrings;

  Xyce::Util::newExpression f1_LHS (std::string("F1(A,B)"), testGroup);
  f1_LHS.lexAndParseExpression();
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);
  f1Expression->setFunctionArgStringVec ( f1ArgStrings );
  f1Expression->lexAndParseExpression();

  std::string f1Name;
  f1_LHS.getFuncPrototypeName(f1Name);

  Xyce::Util::newExpression derivFuncTestExpr(std::string("DDX(F1(P1,3.0),P1)"), testGroup); 
  derivFuncTestExpr.lexAndParseExpression();
  derivFuncTestExpr.attachFunctionNode(f1Name, f1Expression);
  derivFuncTestExpr.attachParameterNode(p1Name , p1Expression);

  Xyce::Util::newExpression copy_derivFuncTestExpr(derivFuncTestExpr); 
  Xyce::Util::newExpression assign_derivFuncTestExpr; 
  assign_derivFuncTestExpr = derivFuncTestExpr; 

  double Bval=2.5;
  solnFuncGroup->setSoln(std::string("B"),Bval);
  std::complex<double> result;
  std::complex<double> refRes = 1.0;

  //derivFuncTestExpr.dumpParseTree(std::cout);

  derivFuncTestExpr.evaluateFunction(result);        EXPECT_EQ(result, refRes);
  copy_derivFuncTestExpr.evaluateFunction(result);   EXPECT_EQ(result, refRes);
  assign_derivFuncTestExpr.evaluateFunction(result); EXPECT_EQ(result, refRes);
}

TEST ( ComplexParserCalculus, derivsThruParams1 )
{
  Teuchos::RCP<solnExpressionGroup2> solnParamGroup = Teuchos::rcp(new solnExpressionGroup2() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnParamGroup;

  // this expression is the RHS of a .global_param statement:  .global_param P1 {V(A)-V(B)}
  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("V(A)-V(B)"), testGroup));
  p1Expression->lexAndParseExpression();

  std::string p1Name = "P1";

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
  std::complex<double> result;
  std::vector<std::complex<double> > derivs;
  std::complex<double> refRes = 2.5;
  std::vector<std::complex<double> > refderivs = { 0.5, -0.5 };

  derivParamTestExpr.evaluate(result,derivs);        EXPECT_EQ( derivs, refderivs );
  copy_derivParamTestExpr.evaluate(result,derivs);   EXPECT_EQ( derivs, refderivs );
  assign_derivParamTestExpr.evaluate(result,derivs); EXPECT_EQ( derivs, refderivs );
}

TEST ( ComplexParserFloor, test1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression floorTest(std::string("floor(10.25)"), testGroup);
  floorTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_floorTest(floorTest); 
  Xyce::Util::newExpression assign_floorTest; 
  assign_floorTest = floorTest; 

  std::complex<double> result;
  floorTest.evaluateFunction(result);        EXPECT_EQ( result, 10.0);
  copy_floorTest.evaluateFunction(result);   EXPECT_EQ( result, 10.0);
  assign_floorTest.evaluateFunction(result); EXPECT_EQ( result, 10.0);
}

TEST ( ComplexParserFloor, test2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression floorTest(std::string("floor(-34.251)"), testGroup);
  floorTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_floorTest(floorTest); 
  Xyce::Util::newExpression assign_floorTest; 
  assign_floorTest = floorTest; 

  std::complex<double> result;
  floorTest.evaluateFunction(result);        EXPECT_EQ( result, -35.0);
  copy_floorTest.evaluateFunction(result);   EXPECT_EQ( result, -35.0);
  assign_floorTest.evaluateFunction(result); EXPECT_EQ( result, -35.0);
}

TEST ( ComplexParserFloor, test3)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression floorTest(std::string("floor(0.71)"), testGroup);
  floorTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_floorTest(floorTest); 
  Xyce::Util::newExpression assign_floorTest; 
  assign_floorTest = floorTest; 

  std::complex<double> result;
  floorTest.evaluateFunction(result);        EXPECT_EQ( result, 0.0);
  copy_floorTest.evaluateFunction(result);   EXPECT_EQ( result, 0.0);
  assign_floorTest.evaluateFunction(result); EXPECT_EQ( result, 0.0);
}

TEST ( ComplexParserCeil, test1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression ceilTest(std::string("ceil(10.25)"), testGroup);
  ceilTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_ceilTest(ceilTest); 
  Xyce::Util::newExpression assign_ceilTest; 
  assign_ceilTest = ceilTest; 

  std::complex<double> result;
  ceilTest.evaluateFunction(result);        EXPECT_EQ( result, 11.0);
  copy_ceilTest.evaluateFunction(result);   EXPECT_EQ( result, 11.0);
  assign_ceilTest.evaluateFunction(result); EXPECT_EQ( result, 11.0);
}

TEST ( ComplexParserCeil, test2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression ceilTest(std::string("ceil(-34.251)"), testGroup);
  ceilTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_ceilTest(ceilTest); 
  Xyce::Util::newExpression assign_ceilTest; 
  assign_ceilTest = ceilTest; 

  std::complex<double> result;
  ceilTest.evaluateFunction(result);        EXPECT_EQ( result, -34.0);
  copy_ceilTest.evaluateFunction(result);   EXPECT_EQ( result, -34.0);
  assign_ceilTest.evaluateFunction(result); EXPECT_EQ( result, -34.0);
}

TEST ( ComplexParserCeil, test3)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression ceilTest(std::string("ceil(0.71)"), testGroup);
  ceilTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_ceilTest(ceilTest); 
  Xyce::Util::newExpression assign_ceilTest; 
  assign_ceilTest = ceilTest; 

  std::complex<double> result;
  ceilTest.evaluateFunction(result);        EXPECT_EQ( result, 1.0);
  copy_ceilTest.evaluateFunction(result);   EXPECT_EQ( result, 1.0);
  assign_ceilTest.evaluateFunction(result); EXPECT_EQ( result, 1.0);
}

TEST ( ComplexParserSpecials, pi1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression piTest(std::string("pi"), testGroup);
  piTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_piTest(piTest); 
  Xyce::Util::newExpression assign_piTest; 
  assign_piTest = piTest; 

  std::complex<double> result;
  piTest.evaluateFunction(result);        EXPECT_EQ( result, M_PI);
  copy_piTest.evaluateFunction(result);   EXPECT_EQ( result, M_PI);
  assign_piTest.evaluateFunction(result); EXPECT_EQ( result, M_PI);
}

TEST ( ComplexParserSpecials, pi2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression piTest(std::string("sin(pi)"), testGroup);
  piTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_piTest(piTest); 
  Xyce::Util::newExpression assign_piTest; 
  assign_piTest = piTest; 

  std::complex<double> result;
  piTest.evaluateFunction(result);        EXPECT_EQ( result, std::sin(M_PI));
  copy_piTest.evaluateFunction(result);   EXPECT_EQ( result, std::sin(M_PI));
  assign_piTest.evaluateFunction(result); EXPECT_EQ( result, std::sin(M_PI));
}

TEST ( ComplexParserSpecials, time)
{
  Teuchos::RCP<timeDepExpressionGroup> timeGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = timeGroup;
  Xyce::Util::newExpression testExpression(std::string("time"),testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  timeGroup->setTime(1.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( (result-(1.0)), 0.0);
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( (result-(1.0)), 0.0);
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( (result-(1.0)), 0.0);

  bool timeDependent = testExpression.getTimeDependent();
  bool copyTimeDependent = copy_testExpression.getTimeDependent();
  bool assignTimeDependent = assign_testExpression.getTimeDependent();

  EXPECT_EQ(timeDependent, true);
  EXPECT_EQ(copyTimeDependent, true);
  EXPECT_EQ(assignTimeDependent, true);

  OUTPUT_MACRO(ComplexParserSpecials, time)
}

TEST ( ComplexParserSpecials, time2)
{
  Teuchos::RCP<timeDepExpressionGroup> timeGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = timeGroup;

  Teuchos::RCP<Xyce::Util::newExpression> tExpression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("time"), testGroup));
  tExpression->lexAndParseExpression();
  std::string tName = "T";

  Xyce::Util::newExpression testExpression(std::string("t"),testGroup);
  testExpression.lexAndParseExpression();
  testExpression.attachParameterNode(tName,tExpression);

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  timeGroup->setTime(1.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( (result-(1.0)), 0.0);
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( (result-(1.0)), 0.0);
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( (result-(1.0)), 0.0);

  bool timeDependent = testExpression.getTimeDependent();
  bool copyTimeDependent = copy_testExpression.getTimeDependent();
  bool assignTimeDependent = assign_testExpression.getTimeDependent();

  EXPECT_EQ(timeDependent, true);
  EXPECT_EQ(copyTimeDependent, true);
  EXPECT_EQ(assignTimeDependent, true);

  OUTPUT_MACRO(ComplexParserSpecials, time2)
}

TEST ( ComplexParserSpecials, freq)
{
  Teuchos::RCP<timeDepExpressionGroup> timeGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = timeGroup;
  Xyce::Util::newExpression testExpression(std::string("freq"),testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  timeGroup->setFreq(1.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( (result-(1.0)), 0.0);
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( (result-(1.0)), 0.0);
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( (result-(1.0)), 0.0);

  bool freqDependent = testExpression.getFreqDependent();
  bool copyFreqDependent = copy_testExpression.getFreqDependent();
  bool assignFreqDependent = assign_testExpression.getFreqDependent();

  EXPECT_EQ(freqDependent, true);
  EXPECT_EQ(copyFreqDependent, true);
  EXPECT_EQ(assignFreqDependent, true);

  OUTPUT_MACRO(ComplexParserSpecials, freq)
}


TEST (ComplexParserSpecials, freq2)
{
  Teuchos::RCP<timeDepExpressionGroup> timeGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = timeGroup;

  Teuchos::RCP<Xyce::Util::newExpression> fExpression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("freq"), testGroup));
  fExpression->lexAndParseExpression();
  std::string fName = "F";

  Xyce::Util::newExpression testExpression(std::string("f"),testGroup);
  testExpression.lexAndParseExpression();
  testExpression.attachParameterNode(fName,fExpression);

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  timeGroup->setFreq(1.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( (result-(1.0)), 0.0);
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( (result-(1.0)), 0.0);
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( (result-(1.0)), 0.0);

  bool freqDependent = testExpression.getFreqDependent();
  bool copyFreqDependent = copy_testExpression.getFreqDependent();
  bool assignFreqDependent = assign_testExpression.getFreqDependent();

  EXPECT_EQ(freqDependent, true);
  EXPECT_EQ(copyFreqDependent, true);
  EXPECT_EQ(assignFreqDependent, true);

  OUTPUT_MACRO(ComplexParserSpecials, freq2)
}

TEST ( ComplexParserSpecials, hertz)
{
  Teuchos::RCP<timeDepExpressionGroup> timeGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = timeGroup;
  Xyce::Util::newExpression testExpression(std::string("hertz"),testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  timeGroup->setFreq(1.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( (result-(1.0)), 0.0);
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( (result-(1.0)), 0.0);
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( (result-(1.0)), 0.0);

  bool freqDependent = testExpression.getFreqDependent();
  bool copyFreqDependent = copy_testExpression.getFreqDependent();
  bool assignFreqDependent = assign_testExpression.getFreqDependent();

  EXPECT_EQ(freqDependent, true);
  EXPECT_EQ(copyFreqDependent, true);
  EXPECT_EQ(assignFreqDependent, true);

  OUTPUT_MACRO(ComplexParserSpecials, hertz)
}

TEST ( ComplexParserSpecials, hertz2)
{
  Teuchos::RCP<timeDepExpressionGroup> timeGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = timeGroup;

  Teuchos::RCP<Xyce::Util::newExpression> fExpression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("hertz"), testGroup));
  fExpression->lexAndParseExpression();
  std::string fName = "F";

  Xyce::Util::newExpression testExpression(std::string("f"),testGroup);
  testExpression.lexAndParseExpression();
  testExpression.attachParameterNode(fName,fExpression);

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  timeGroup->setFreq(1.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( (result-(1.0)), 0.0);
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( (result-(1.0)), 0.0);
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( (result-(1.0)), 0.0);

  bool freqDependent = testExpression.getFreqDependent();  
  bool copyFreqDependent = copy_testExpression.getFreqDependent();
  bool assignFreqDependent = assign_testExpression.getFreqDependent();

  EXPECT_EQ(freqDependent, true);
  EXPECT_EQ(copyFreqDependent, true);
  EXPECT_EQ(assignFreqDependent, true);

  OUTPUT_MACRO(ComplexParserSpecials, hertz2)
}

TEST ( ComplexParserSpecials, gmin)
{
  Teuchos::RCP<timeDepExpressionGroup> timeGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = timeGroup;
  Xyce::Util::newExpression testExpression(std::string("gmin"),testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  timeGroup->setGmin(1.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( (result-(1.0)), 0.0);
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( (result-(1.0)), 0.0);
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( (result-(1.0)), 0.0);

  bool gminDependent = testExpression.getGminDependent();
  bool copyGminDependent = copy_testExpression.getGminDependent();
  bool assignGminDependent = assign_testExpression.getGminDependent();

  EXPECT_EQ(gminDependent, true);
  EXPECT_EQ(copyGminDependent, true);
  EXPECT_EQ(assignGminDependent, true);

  OUTPUT_MACRO(ComplexParserSpecials, gmin)
}

TEST ( ComplexParserSpecials, gmin2)
{
  Teuchos::RCP<timeDepExpressionGroup> timeGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = timeGroup;

  Teuchos::RCP<Xyce::Util::newExpression> gExpression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("gmin"), testGroup));
  gExpression->lexAndParseExpression();
  std::string gName = "G";

  Xyce::Util::newExpression testExpression(std::string("g"),testGroup);
  testExpression.lexAndParseExpression();
  testExpression.attachParameterNode(gName,gExpression);

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  timeGroup->setGmin(1.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( (result-(1.0)), 0.0);
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( (result-(1.0)), 0.0);
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( (result-(1.0)), 0.0);

  bool gminDependent = testExpression.getGminDependent();  
  bool copyGminDependent = copy_testExpression.getGminDependent();
  bool assignGminDependent = assign_testExpression.getGminDependent();

  EXPECT_EQ(gminDependent, true);
  EXPECT_EQ(copyGminDependent, true);
  EXPECT_EQ(assignGminDependent, true);

  OUTPUT_MACRO(ComplexParserSpecials, gmin2)
}

TEST (ComplexParserSpecials, temp)
{
  Teuchos::RCP<tempDepExpressionGroup> tempGroup = Teuchos::rcp(new tempDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = tempGroup;
  Xyce::Util::newExpression testExpression(std::string("temp"),testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  tempGroup->setTemp(1.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( (result-(1.0)), 0.0);
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( (result-(1.0)), 0.0);
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( (result-(1.0)), 0.0);

  bool tempDependent = testExpression.getTempDependent();  
  bool copyTempDependent = copy_testExpression.getTempDependent();
  bool assignTempDependent = assign_testExpression.getTempDependent();

  EXPECT_EQ(tempDependent, true);
  EXPECT_EQ(copyTempDependent, true);
  EXPECT_EQ(assignTempDependent, true);

  OUTPUT_MACRO(ComplexParserSpecials, temp)
}

TEST (ComplexParserSpecials, temp2)
{
  Teuchos::RCP<tempDepExpressionGroup> tempGroup = Teuchos::rcp(new tempDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = tempGroup;

  Teuchos::RCP<Xyce::Util::newExpression> tExpression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("temp"), testGroup));
  tExpression->lexAndParseExpression();
  std::string tName = "T";

  Xyce::Util::newExpression testExpression(std::string("t"),testGroup);
  testExpression.lexAndParseExpression();
  testExpression.attachParameterNode(tName,tExpression);

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  tempGroup->setTemp(1.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( (result-(1.0)), 0.0);
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( (result-(1.0)), 0.0);
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( (result-(1.0)), 0.0);

  bool tempDependent = testExpression.getTempDependent();  
  bool copyTempDependent = copy_testExpression.getTempDependent();
  bool assignTempDependent = assign_testExpression.getTempDependent();

  EXPECT_EQ(tempDependent, true);
  EXPECT_EQ(copyTempDependent, true);
  EXPECT_EQ(assignTempDependent, true);

  OUTPUT_MACRO(ComplexParserSpecials, temp2)
}


TEST ( ComplexParserSpecials, temper)
{
  Teuchos::RCP<tempDepExpressionGroup> tempGroup = Teuchos::rcp(new tempDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = tempGroup;
  Xyce::Util::newExpression testExpression(std::string("temper"),testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  tempGroup->setTemp(1.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( (result-(1.0)), 0.0);
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( (result-(1.0)), 0.0);
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( (result-(1.0)), 0.0);

  bool tempDependent = testExpression.getTempDependent();  
  bool copyTempDependent = copy_testExpression.getTempDependent();
  bool assignTempDependent = assign_testExpression.getTempDependent();

  EXPECT_EQ(tempDependent, true);
  EXPECT_EQ(copyTempDependent, true);
  EXPECT_EQ(assignTempDependent, true);

  OUTPUT_MACRO(ComplexParserSpecials, temper)
}

TEST ( ComplexParserSpecials, temper2)
{
  Teuchos::RCP<tempDepExpressionGroup> tempGroup = Teuchos::rcp(new tempDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = tempGroup;

  Teuchos::RCP<Xyce::Util::newExpression> tExpression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("temper"), testGroup));
  tExpression->lexAndParseExpression();
  std::string tName = "T";

  Xyce::Util::newExpression testExpression(std::string("t"),testGroup);
  testExpression.lexAndParseExpression();
  testExpression.attachParameterNode(tName,tExpression);

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  tempGroup->setTemp(1.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( (result-(1.0)), 0.0);
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( (result-(1.0)), 0.0);
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( (result-(1.0)), 0.0);

  bool tempDependent = testExpression.getTempDependent();  
  bool copyTempDependent = copy_testExpression.getTempDependent();
  bool assignTempDependent = assign_testExpression.getTempDependent();

  EXPECT_EQ(tempDependent, true);
  EXPECT_EQ(copyTempDependent, true);
  EXPECT_EQ(assignTempDependent, true);

  OUTPUT_MACRO(ComplexParserSpecials, temper2)
}


//
TEST (ComplexParserSpecials, vt)
{
  Teuchos::RCP<tempDepExpressionGroup> vtGroup = Teuchos::rcp(new tempDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = vtGroup;
  Xyce::Util::newExpression testExpression(std::string("vt"),testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  vtGroup->setVT(1.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( (result-(1.0)), 0.0);
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( (result-(1.0)), 0.0);
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( (result-(1.0)), 0.0);

  bool vtDependent = testExpression.getVTDependent();  
  bool copyVTDependent = copy_testExpression.getVTDependent();
  bool assignVTDependent = assign_testExpression.getVTDependent();

  EXPECT_EQ(vtDependent, true);
  EXPECT_EQ(copyVTDependent, true);
  EXPECT_EQ(assignVTDependent, true);

  OUTPUT_MACRO(ComplexParserSpecials, vt)
}

TEST (ComplexParserSpecials, vt2)
{
  Teuchos::RCP<tempDepExpressionGroup> vtGroup = Teuchos::rcp(new tempDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = vtGroup;

  Teuchos::RCP<Xyce::Util::newExpression> tExpression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("vt"), testGroup));
  tExpression->lexAndParseExpression();
  std::string tName = "T";

  Xyce::Util::newExpression testExpression(std::string("t"),testGroup);
  testExpression.lexAndParseExpression();
  testExpression.attachParameterNode(tName,tExpression);

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  vtGroup->setVT(1.0);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( (result-(1.0)), 0.0);
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( (result-(1.0)), 0.0);
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( (result-(1.0)), 0.0);

  bool vtDependent = testExpression.getVTDependent();  
  bool copyVTDependent = copy_testExpression.getVTDependent();
  bool assignVTDependent = assign_testExpression.getVTDependent();

  EXPECT_EQ(vtDependent, true);
  EXPECT_EQ(copyVTDependent, true);
  EXPECT_EQ(assignVTDependent, true);

  OUTPUT_MACRO(ComplexParserSpecials, vt2)
}

TEST (  ComplexParserSpecials, CtoK1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("CtoK"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  std::complex<double> refRes(CONSTCtoK);
  testExpression.evaluateFunction(result);
  EXPECT_EQ(result, refRes);

  Xyce::Util::newExpression copyExpression(testExpression); 
  copyExpression.evaluateFunction(result); 
  EXPECT_EQ(result, refRes);

  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 
  assignExpression.evaluateFunction(result); 
  EXPECT_EQ(result, refRes);
}

TEST ( ComplexParserSpecials, exp1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression expTest(std::string("exp"), testGroup);
  expTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_expTest(expTest); 
  Xyce::Util::newExpression assign_expTest; 
  assign_expTest = expTest; 

  std::complex<double> result;
  expTest.evaluateFunction(result);        EXPECT_EQ( result, std::exp(1.0));
  copy_expTest.evaluateFunction(result);   EXPECT_EQ( result, std::exp(1.0));
  assign_expTest.evaluateFunction(result); EXPECT_EQ( result, std::exp(1.0));
  OUTPUT_MACRO2(ComplexParserSpecials, exp1, expTest) 
}

TEST ( ComplexParserSpecials, exp2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression expTest(std::string("5.0*exp"), testGroup);
  expTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_expTest(expTest); 
  Xyce::Util::newExpression assign_expTest; 
  assign_expTest = expTest; 

  std::complex<double> result;
  expTest.evaluateFunction(result);        EXPECT_EQ( result, (5.0*std::exp(1.0)));
  copy_expTest.evaluateFunction(result);   EXPECT_EQ( result, (5.0*std::exp(1.0)));
  assign_expTest.evaluateFunction(result); EXPECT_EQ( result, (5.0*std::exp(1.0)));
  OUTPUT_MACRO2(ComplexParserSpecials, exp2, expTest) 
}

TEST ( ComplexParserSpecials, exp3)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression expTest(std::string("exp*exp(2.0)"), testGroup);
  expTest.lexAndParseExpression();

  Xyce::Util::newExpression copy_expTest(expTest); 
  Xyce::Util::newExpression assign_expTest; 
  assign_expTest = expTest; 

  std::complex<double> result;
  expTest.evaluateFunction(result);        EXPECT_EQ( result, (std::exp(1.0)*std::exp(2.0)));
  copy_expTest.evaluateFunction(result);   EXPECT_EQ( result, (std::exp(1.0)*std::exp(2.0)));
  assign_expTest.evaluateFunction(result); EXPECT_EQ( result, (std::exp(1.0)*std::exp(2.0)));
  OUTPUT_MACRO2(ComplexParserSpecials, exp3, expTest) 
}




//


// these next two tests are for the use case of a parameter that is named either "I" or "V".
// In an earlier implementation, the parser would get confused by this.  The string res*I*I, 
// which is used by the regression test BUG_645_SON/bug645son_Func.cir, would simply fail to 
// parse and wouldn't even issue an error.  It would proceed but then fail to compute anything.
//

TEST ( ComplexParserParamTest, I )
{
  Teuchos::RCP<testExpressionGroupWithParamSupport> paramGroup = Teuchos::rcp(new testExpressionGroupWithParamSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;

  Teuchos::RCP<Xyce::Util::newExpression> iExpression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("2+3J"), testGroup));
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

  std::complex<double> iVal(2.0,3.0);
  std::complex<double> resVal(4.0,0.0);

  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( result, iVal*iVal*resVal );
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( result, iVal*iVal*resVal );
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( result, iVal*iVal*resVal );
  OUTPUT_MACRO(ComplexParserParamTest, I)
}

TEST ( ComplexParserParamTest, V )
{
  Teuchos::RCP<testExpressionGroupWithParamSupport> paramGroup = Teuchos::rcp(new testExpressionGroupWithParamSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = paramGroup;

  Teuchos::RCP<Xyce::Util::newExpression> vExpression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("2+3J"), testGroup));
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

  std::complex<double> vVal(2.0,3.0);
  std::complex<double> resVal(4.0,0.0);

  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);        EXPECT_EQ( result, vVal*vVal*resVal );
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( result, vVal*vVal*resVal );
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( result, vVal*vVal*resVal );
  OUTPUT_MACRO(ComplexParserParamTest, V)
}

// parameter node replacement tests
TEST ( ComplexParserParamTest, replace1)
{
  Teuchos::RCP<testExpressionGroup> noparamGroup = Teuchos::rcp(new testExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = noparamGroup;

  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("2+3"), testGroup));
  p1Expression->lexAndParseExpression();
  std::string p1Name = "p1";

  Xyce::Util::newExpression testExpression(std::string("p1"), testGroup);
  testExpression.lexAndParseExpression();
  testExpression.replaceParameterNode(p1Name,p1Expression);

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  std::complex<double> result;
  testExpression.evaluateFunction(result);        EXPECT_EQ( std::real(result), 5.0 );
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( std::real(result), 5.0 );
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( std::real(result), 5.0 );
  OUTPUT_MACRO(ComplexParserParamTest, replace1)
}

TEST ( ComplexParserParamTest, replace2)
{
  Teuchos::RCP<testExpressionGroup> noparamGroup = Teuchos::rcp(new testExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = noparamGroup;

  Teuchos::RCP<Xyce::Util::newExpression> p1Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("2+3"), testGroup));
  p1Expression->lexAndParseExpression();
  std::string p1Name = "p1";

  Teuchos::RCP<Xyce::Util::newExpression> p0Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("p1"), testGroup));
  p0Expression->lexAndParseExpression();
  p0Expression->attachParameterNode(p1Name,p1Expression);
  std::string p0Name = "p0";

  Xyce::Util::newExpression testExpression(std::string("p0"), testGroup);
  testExpression.lexAndParseExpression();
  testExpression.replaceParameterNode(p0Name,p0Expression);

  Xyce::Util::newExpression copy_testExpression(testExpression); 
  Xyce::Util::newExpression assign_testExpression; 
  assign_testExpression = testExpression; 

  std::complex<double> result;
  testExpression.evaluateFunction(result);        EXPECT_EQ( std::real(result), 5.0 );
  copy_testExpression.evaluateFunction(result);   EXPECT_EQ( std::real(result), 5.0 );
  assign_testExpression.evaluateFunction(result); EXPECT_EQ( std::real(result), 5.0 );
  OUTPUT_MACRO(ComplexParserParamTest, replace2)
}

TEST ( ComplexParserASCTHTest, test0)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("cosh(V(A))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=-10.0;
  std::complex<double> refRes = std::cosh(Aval);
  solnGroup->setSoln(std::string("A"),Aval);

  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result, derivs);   
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refRes));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refRes));

  copyExpression.evaluate(result, derivs);   
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refRes));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refRes));

  assignExpression.evaluate(result, derivs); 
  EXPECT_DOUBLE_EQ( std::real(result), std::real(refRes));
  EXPECT_DOUBLE_EQ( std::imag(result), std::imag(refRes));
}

TEST ( ComplexParserASCTHTest, test1)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("acosh(cosh(V(A)))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=-10.0;
  std::complex<double> refRes = std::acosh(std::cosh(Aval));
  solnGroup->setSoln(std::string("A"),Aval);

  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  copyExpression.evaluate(result,derivs);   
  EXPECT_EQ( result, refRes);
  assignExpression.evaluate(result,derivs); 
  EXPECT_EQ( result, refRes);
}

TEST ( ComplexParserASCTHTest, test2)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("acosh(cosh(V(A)))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=0.0;
  std::complex<double> refRes = std::acosh(std::cosh(0.0));
  solnGroup->setSoln(std::string("A"),Aval);

  // this double checks if the derivatives are NOT Nan.
  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  for( auto i=0; i<derivs.size(); i++ )
  {
    EXPECT_FALSE( std::isnan( derivs[i].real() ));
    EXPECT_FALSE( std::isnan( derivs[i].imag() ));
  }
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  for( auto i=0; i<derivs.size(); i++ )
  {
    EXPECT_FALSE( std::isnan( derivs[i].real() ));
    EXPECT_FALSE( std::isnan( derivs[i].imag() ));
  }
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  for( auto i=0; i<derivs.size(); i++ )
  {
    EXPECT_FALSE( std::isnan( derivs[i].real() ));
    EXPECT_FALSE( std::isnan( derivs[i].imag() ));
  }
}


TEST ( ComplexParserSTPTest, test1)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("stp(V(A))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=2.0;
  std::complex<double> refRes = 1.0;
  solnGroup->setSoln(std::string("A"),Aval);

  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(ComplexParserSTPTest, test1)
}

TEST ( ComplexParserSTPTest, test2)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("stp(V(A))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=-2.0;
  std::complex<double>  refRes = 0.0;
  solnGroup->setSoln(std::string("A"),Aval);

  std::vector<std::complex<double> > derivs;
  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  OUTPUT_MACRO(ComplexParserSTPTest, test2)
}


TEST ( ComplexParseratan2Test, test1)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("atan2(V(A),V(B))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=0.5, Bval=0.25;
  std::complex<double>  refRes = std::atan2(std::real(Aval),std::real(Bval));
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("B"),Bval);

  std::vector<std::complex<double> > derivs;

  std::complex<double>  denom = Aval*Aval+Bval*Bval;
  std::vector<std::complex<double> > refderivs = { (Bval/denom), (-Aval/denom) };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(ComplexParseratan2Test, test1)
}


TEST ( ComplexParseratan2Test, test2)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("atan2(V(A),V(B))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=0.5, Bval=-0.25;
  std::complex<double>  refRes = std::atan2(std::real(Aval),std::real(Bval));
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("B"),Bval);

  std::vector<std::complex<double> > derivs;

  std::complex<double>  denom = Aval*Aval+Bval*Bval;
  std::vector<std::complex<double> > refderivs = { (Bval/denom), (-Aval/denom) };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(ComplexParseratan2Test, test2)
}

TEST ( ComplexParseratan2Test, test3)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("atan2(V(A),V(B))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=-0.5, Bval=0.25;
  std::complex<double>  refRes = std::atan2(std::real(Aval),std::real(Bval));
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("B"),Bval);

  std::vector<std::complex<double> > derivs;

  std::complex<double>  denom = Aval*Aval+Bval*Bval;
  std::vector<std::complex<double> > refderivs = { (Bval/denom), (-Aval/denom) };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(ComplexParseratan2Test, test3)
}

TEST ( ComplexParseratan2Test, test4)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("atan2(V(A),V(B))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=-0.5, Bval=-0.25;
  std::complex<double>  refRes = std::atan2(std::real(Aval),std::real(Bval));
  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("B"),Bval);

  std::vector<std::complex<double> > derivs;

  std::complex<double>  denom = Aval*Aval+Bval*Bval;
  std::vector<std::complex<double> > refderivs = { (Bval/denom), (-Aval/denom) };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(ComplexParseratan2Test, test4)
}




TEST ( ComplexParserTwoNodeDerivTest, test1)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("V(A,B)"), testGroup);
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, vAval=5.0, vBval=1.0;

  solnGroup->setSoln(std::string("A"),vAval);
  solnGroup->setSoln(std::string("B"),vBval);

  std::complex<double> refRes = (vAval-vBval);

  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refderivs = { 1, -1 };

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

//
TEST ( ComplexParserpolyTest, test1)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("POLY(1) V(A) 1.0 2.0"), testGroup);
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=5.0;
  solnGroup->setSoln(std::string("A"),Aval);
  std::complex<double>  refRes = 1.0 + 2.0*Aval;

  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refderivs = { 2.0 };

  testExpression.evaluate(result, derivs);
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(ComplexParserpolyTest, test1)
}


// poly test taken from amp_mod.cir
//
// BAM = 0 + 0*V(10) + 0*V(20) + 0*V(10)^2 + 1*V(10)*V(20) 
//
TEST ( ComplexParserpolyTest, test2)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("POLY(2) V(A) V(B) 0 0 0 0 1 "), testGroup);
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, vAval=5.0, vBval=1.0;

  solnGroup->setSoln(std::string("A"),vAval);
  solnGroup->setSoln(std::string("B"),vBval);

  std::complex<double>  refRes = vAval*vBval;

  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refderivs = {  vBval, vAval };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(ComplexParserpolyTest, test2)
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
TEST ( ComplexParserpolyTest, test3)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("POLY(2) V(A,B) V(C) 0 0 0 0 1 "), testGroup);
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, vAval=5.0, vBval=1.0, vCval=2.0;

  solnGroup->setSoln(std::string("A"),vAval);
  solnGroup->setSoln(std::string("B"),vBval);
  solnGroup->setSoln(std::string("C"),vCval);

  std::complex<double>  refRes = 
    (vAval-vBval)*vCval;

  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refderivs = { vCval, -vCval, (vAval-vBval) };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(ComplexParserpolyTest, test3)
}

TEST ( ComplexParserpolyTest, test4)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("POLY(1) V(A,B) 0 1  "), testGroup);
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, vAval=5.0, vBval=1.0;

  solnGroup->setSoln(std::string("A"),vAval);
  solnGroup->setSoln(std::string("B"),vBval);

  std::complex<double>  refRes = (vAval-vBval);

  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refderivs = { 1, -1 };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(ComplexParserpolyTest, test4)
}

// POLY(5) I(VB) I(VC)  I(VE)  I(VLP) I(VLN)  0 10.61E6 -10E6 10E6 10E6 -10E6
TEST ( ComplexParserpolyTest, test5)
{
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;

  // with parens
  Xyce::Util::newExpression testExpression(std::string("POLY(5) I(VB) I(VC)  I(VE)  I(VLP) I(VLN)  0 (10.61E6) (-10E6) (10E6) (10E6) (-10E6)"), testGroup);
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=std::complex<double>(0.0,0.0); 
  std::complex<double>  VB=std::complex<double>(-1.0,0.0);
  std::complex<double>  VC=std::complex<double>(-2.0,0.0);
  std::complex<double>  VE=std::complex<double>(-3.0,0.0);
  std::complex<double>  VLP=std::complex<double>(-4.0,0.0);
  std::complex<double>  VLN=std::complex<double>(-5.0,0.0);

  solnGroup->setSoln(std::string("VB"),VB);
  solnGroup->setSoln(std::string("VC"),VC);
  solnGroup->setSoln(std::string("VE"),VE);
  solnGroup->setSoln(std::string("VLP"),VLP);
  solnGroup->setSoln(std::string("VLN"),VLN);

  std::complex<double>  refRes = 0 +  (10.61E6)*(-1) +  (-10E6)*(-2) +  (10E6)*(-3) +  (10E6)*(-4) +  (-10E6)*(-5);

  std::vector<std::complex<double> > derivs;
  //std::vector<std::complex<double> > refderivs = { 1, -1 };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(ComplexParserpolyTest, test5)
}

// POLY(5) I(VB) I(VC)  I(VE)  I(VLP) I(VLN)  0 10.61E6 -10E6 10E6 10E6 -10E6
TEST ( ComplexParserpolyTest, test6)
{
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;

  // without parens
  Xyce::Util::newExpression testExpression(std::string("POLY(5) I(VB) I(VC)  I(VE)  I(VLP) I(VLN)  0 10.61E6 -10E6 10E6 10E6 -10E6"), testGroup);
  //
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0; 
  std::complex<double>  VB=-1.0;
  std::complex<double>  VC=-2.0;
  std::complex<double>  VE=-3.0;
  std::complex<double>  VLP=-4.0;
  std::complex<double>  VLN=-5.0;

  solnGroup->setSoln(std::string("VB"),VB);
  solnGroup->setSoln(std::string("VC"),VC);
  solnGroup->setSoln(std::string("VE"),VE);
  solnGroup->setSoln(std::string("VLP"),VLP);
  solnGroup->setSoln(std::string("VLN"),VLN);

  std::complex<double>  refRes = 0 +  (10.61E6)*(-1) +  (-10E6)*(-2) +  (10E6)*(-3) +  (10E6)*(-4) +  (-10E6)*(-5);

  std::vector<std::complex<double> > derivs;
  //std::vector<std::complex<double> > refderivs = { 1, -1 };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);
  OUTPUT_MACRO(ComplexParserpolyTest, test6)
}

// POLY(5) I(VB) I(VC)  I(VE)  I(VLP) I(VLN)  0 10.61E6 -10E6 10E6 10E6 -10E6
TEST ( ComplexParserpolyTest, test7)
{
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;

  // without parens
  Xyce::Util::newExpression testExpression(std::string("POLY(5) I(VB) I(VC)  I(VE)  I(VLP) I(VLN)  +0 +10.61E6 -10E6 +10E6 +10E6 -10E6"), testGroup);
  //
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0; 
  std::complex<double>  VB=-1.0;
  std::complex<double>  VC=-2.0;
  std::complex<double>  VE=-3.0;
  std::complex<double>  VLP=-4.0;
  std::complex<double>  VLN=-5.0;

  solnGroup->setSoln(std::string("VB"),VB);
  solnGroup->setSoln(std::string("VC"),VC);
  solnGroup->setSoln(std::string("VE"),VE);
  solnGroup->setSoln(std::string("VLP"),VLP);
  solnGroup->setSoln(std::string("VLN"),VLN);

  std::complex<double>  refRes = 0 +  (10.61E6)*(-1) +  (-10E6)*(-2) +  (10E6)*(-3) +  (10E6)*(-4) +  (-10E6)*(-5);

  std::vector<std::complex<double> > derivs;
  //std::vector<std::complex<double> > refderivs = { 1, -1 };

  testExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  //EXPECT_EQ( derivs, refderivs);

  OUTPUT_MACRO(ComplexParserpolyTest, test7)
}

TEST ( ComplexParserpolyTest, test8)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("POLY(1) V(A) {1.0} {2.0}"), testGroup);
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=5.0;
  solnGroup->setSoln(std::string("A"),Aval);
  std::complex<double>  refRes = 1.0 + 2.0*Aval;

  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refderivs = { 2.0 };

  testExpression.evaluate(result, derivs);
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);

  OUTPUT_MACRO(ComplexParserpolyTest, test8)
}

TEST ( ComplexParserpolyTest, test9)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("POLY(1) V(A) ((1.0)) ((2.0))"), testGroup);
  testExpression.lexAndParseExpression();

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=5.0;
  solnGroup->setSoln(std::string("A"),Aval);
  std::complex<double>  refRes = 1.0 + 2.0*Aval;

  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refderivs = { 2.0 };

  testExpression.evaluate(result, derivs);
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  copyExpression.evaluate(result, derivs);   
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);
  assignExpression.evaluate(result, derivs); 
  EXPECT_EQ( result, refRes);
  EXPECT_EQ( derivs, refderivs);

  OUTPUT_MACRO(ComplexParserpolyTest, test9)
}

TEST ( ComplexParserpolyTest, test10)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
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

  testExpression.attachParameterNode(p1Name,p1Expression);
  testExpression.attachParameterNode(p2Name,p2Expression);

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=5.0;
  solnGroup->setSoln(std::string("A"),Aval);
  std::complex<double>  refRes = 1.0 + 2.0*Aval;

  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refderivs = { 2.0 };

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

// poly test with parameters and a function
// This test is named "test10a" b/c it is nearly identical to test10, other than use of .func.
TEST ( ComplexParserpolyTest, test10a)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  //Xyce::Util::newExpression testExpression(std::string("POLY(1) V(A) p1 p2"), testGroup);
  Xyce::Util::newExpression testExpression(std::string("f1(V(A))"), testGroup);
  testExpression.lexAndParseExpression();

  // this expression is the RHS of a .func statement:  .func F1(A) {POLY ...}
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression  
    = Teuchos::rcp(new Xyce::Util::newExpression(std::string("POLY(1) X p1 p2"), testGroup));

  Xyce::Util::newExpression f1_LHS (std::string("F1(X)"), testGroup);
  f1_LHS.lexAndParseExpression();

  std::vector<std::string> f1ArgStrings ;
  f1_LHS.getFuncPrototypeArgStrings(f1ArgStrings);
  f1Expression->setFunctionArgStringVec (f1ArgStrings);
  f1Expression->lexAndParseExpression();
  std::string f1Name;
  f1_LHS.getFuncPrototypeName(f1Name);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  Teuchos::RCP<Xyce::Util::newExpression> p1Expression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("1.0"), testGroup));
  Teuchos::RCP<Xyce::Util::newExpression> p2Expression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("2.0"), testGroup));
  p1Expression->lexAndParseExpression();
  p2Expression->lexAndParseExpression();

  std::string p1Name = "p1";
  std::string p2Name = "p2";

  f1Expression->attachParameterNode(p1Name,p1Expression);
  f1Expression->attachParameterNode(p2Name,p2Expression);

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, Aval=5.0;
  solnGroup->setSoln(std::string("A"),Aval);
  std::complex<double> refRes = 1.0 + 2.0*Aval;

  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refderivs = { 2.0 };

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

// single input, 4th order POLY test, with parameters
TEST ( ComplexParserpolyTest, test11)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("3.0*(POLY(1) V(A) p1 p2 p3 p4 p5)"), testGroup);
  testExpression.lexAndParseExpression();

  Teuchos::RCP<Xyce::Util::newExpression> p1Expression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("3.0"), testGroup));

  Teuchos::RCP<Xyce::Util::newExpression> p2Expression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("2.0"), testGroup));

  Teuchos::RCP<Xyce::Util::newExpression> p3Expression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("1.0"), testGroup));

  Teuchos::RCP<Xyce::Util::newExpression> p4Expression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("3.0"), testGroup));

  Teuchos::RCP<Xyce::Util::newExpression> p5Expression
    = Teuchos::rcp(new Xyce::Util::newExpression (std::string("4.0"), testGroup));

  p1Expression->lexAndParseExpression();
  p2Expression->lexAndParseExpression();
  p3Expression->lexAndParseExpression();
  p4Expression->lexAndParseExpression();
  p5Expression->lexAndParseExpression();

  std::string p1Name = "p1";
  std::string p2Name = "p2";
  std::string p3Name = "p3";
  std::string p4Name = "p4";
  std::string p5Name = "p5";

  testExpression.attachParameterNode(p1Name,p1Expression);
  testExpression.attachParameterNode(p2Name,p2Expression);
  testExpression.attachParameterNode(p3Name,p3Expression);
  testExpression.attachParameterNode(p4Name,p4Expression);
  testExpression.attachParameterNode(p5Name,p5Expression);

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double>  result=0.0, Aval=5.0;
  solnGroup->setSoln(std::string("A"),Aval);

  std::complex<double>  refRes = 3.0*((3.0+2.0*Aval+std::pow(Aval,2.0)+3.0*std::pow(Aval,3.0)+4.0*std::pow(Aval,4.0)));

  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refderivs = { 3.0*(( 2.0+ 2.0*Aval+ 9.0*std::pow(Aval,2.0)+ 4.0*4.0*std::pow(Aval,3.0))) };

  testExpression.evaluate(result, derivs);
  EXPECT_FLOAT_EQ( std::real(result), std::real(refRes));
  EXPECT_FLOAT_EQ( std::real(derivs[0]), std::real(refderivs[0]));

  copyExpression.evaluate(result, derivs);   
  EXPECT_FLOAT_EQ( std::real(result), std::real(refRes));
  EXPECT_FLOAT_EQ( std::real(derivs[0]), std::real(refderivs[0]));

  assignExpression.evaluate(result, derivs); 
  EXPECT_FLOAT_EQ( std::real(result), std::real(refRes));
  EXPECT_FLOAT_EQ( std::real(derivs[0]), std::real(refderivs[0]));
}

// two input, 3rd order POLY test
TEST ( ComplexParserpolyTest, test12)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;
  Xyce::Util::newExpression testExpression(std::string("POLY(2) V(A) V(B) C0 C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C12 C13 C14 C15 "), testGroup);
  testExpression.lexAndParseExpression();

  // originals:
   std::vector<double> C = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15};

  std::vector < Teuchos::RCP<Xyce::Util::newExpression> > C_param_expressions;

  C_param_expressions.push_back( Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.1"), testGroup)) );
  C_param_expressions.push_back( Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.2"), testGroup)) );
  C_param_expressions.push_back( Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.3"), testGroup)) );
  C_param_expressions.push_back( Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.4"), testGroup)) );
  C_param_expressions.push_back( Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.5"), testGroup)) );
  C_param_expressions.push_back( Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.6"), testGroup)) );
  C_param_expressions.push_back( Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.7"), testGroup)) );
  C_param_expressions.push_back( Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.8"), testGroup)) );
  C_param_expressions.push_back( Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.9"), testGroup)) );
  C_param_expressions.push_back( Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.10"), testGroup)) );
  C_param_expressions.push_back( Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.11"), testGroup)) );
  C_param_expressions.push_back( Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.12"), testGroup)) );
  C_param_expressions.push_back( Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.13"), testGroup)) );
  C_param_expressions.push_back( Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.14"), testGroup)) );
  C_param_expressions.push_back( Teuchos::rcp(new Xyce::Util::newExpression (std::string("0.15"), testGroup)) );

  for (int ii=0;ii<C_param_expressions.size();ii++) { C_param_expressions[ii]->lexAndParseExpression(); }

  for (int ii=0;ii<C_param_expressions.size();ii++) 
  { 
    std::string number = std::to_string(ii);
    std::string name = "C" + number;
    testExpression.attachParameterNode(name, C_param_expressions[ii]);
  }

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0; 
  std::vector<std::complex<double> > V = {5.0, 2.0};
  solnGroup->setSoln(std::string("A"),V[0]);
  solnGroup->setSoln(std::string("B"),V[1]);

  std::complex<double> refRes= C[0] + C[1]*V[0] + C[2]*V[1] + C[3]*V[0]*V[0] + C[4]*V[0]*V[1] + C[5]*V[0]*V[1] + C[6]*V[1]*V[1] + ( C[7] * V[0]*V[0]*V[0] + C[8] * V[0]*V[0]*V[1] + C[9] * V[0]*V[1]*V[0] + C[10] * V[0]*V[1]*V[1] + C[11] * V[1]*V[0]*V[0] + C[12] * V[1]*V[0]*V[1] + C[13] * V[1]*V[1]*V[0] + C[14] * V[1]*V[1]*V[1] );



  // deriv w.r.t. V[0]
  std::complex<double> deriv0 =  C[1] + 2.0*C[3]*V[0] + C[4]*V[1] + C[5]*V[1] + ( 3.0* C[7] * V[0]*V[0] + 2.0*C[8] * V[0]*V[1] + 2.0*C[9] * V[1]*V[0] + C[10] * V[1]*V[1] + 2.0*C[11] * V[1]*V[0] + C[12] * V[1]*V[1] + C[13] * V[1]*V[1] );

  // deriv w.r.t. V[1]
  std::complex<double> deriv1 = + C[2] + C[4]*V[0] + C[5]*V[0] + 2.0*C[6]*V[1] + ( C[8] * V[0]*V[0] + C[9] * V[0]*V[0] + 2.0*C[10] * V[0]*V[1] + C[11] * V[0]*V[0] + 2.0*C[12] * V[1]*V[0] + 2.0*C[13] * V[1]*V[0] + 3.0*C[14] * V[1]*V[1] );


  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refderivs = {deriv0, deriv1};

  testExpression.evaluate(result, derivs);
  EXPECT_FLOAT_EQ(std::real(result),std::real(refRes));
  EXPECT_FLOAT_EQ(std::real(derivs[0]),std::real(refderivs[0]));
  EXPECT_FLOAT_EQ(std::real(derivs[1]),std::real(refderivs[1]));

  copyExpression.evaluate(result, derivs);   
  EXPECT_FLOAT_EQ(std::real(result),std::real(refRes));
  EXPECT_FLOAT_EQ(std::real(derivs[0]),std::real(refderivs[0]));
  EXPECT_FLOAT_EQ(std::real(derivs[1]),std::real(refderivs[1]));

  assignExpression.evaluate(result, derivs); 
  EXPECT_FLOAT_EQ(std::real(result),std::real(refRes));
  EXPECT_FLOAT_EQ(std::real(derivs[0]),std::real(refderivs[0]));
  EXPECT_FLOAT_EQ(std::real(derivs[1]),std::real(refderivs[1]));
}

void setupCoordinates ( int order, int localIndex, int size, std::vector<int> & coords)
{
  coords.resize(order,0);
  int p=size;
  int d=order;
  int indx = localIndex;
  for (int k = d-1; k >= 0; k--) 
  {
    coords[k] = indx % (p);
    indx /= p;
  }
}

// three input, 3rd order POLY test
TEST ( ComplexParserpolyTest, test13)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;

  std::string exprString = std::string("POLY(3) V(A) V(B) V(C)");
  int numCoefs = 40;
  std::vector<std::complex<double> > C(numCoefs);
  std::vector < Teuchos::RCP<Xyce::Util::newExpression> > C_param_expressions;
  for (int ii=0;ii<numCoefs;ii++)
  {
    std::string number = std::to_string(ii);
    std::string name = " C" + number;
    exprString += name;
    C[ii] =  (static_cast<std::complex<double> >(ii+1))*0.1;

    std::string paramNumber = std::to_string(std::real(C[ii]));
    C_param_expressions.push_back( Teuchos::rcp(new Xyce::Util::newExpression (paramNumber, testGroup)) );
  }
  Xyce::Util::newExpression testExpression(exprString, testGroup);
  testExpression.lexAndParseExpression();

  for (int ii=0;ii<C_param_expressions.size();ii++) { C_param_expressions[ii]->lexAndParseExpression(); }

  for (int ii=0;ii<C_param_expressions.size();ii++) 
  { 
    std::string number = std::to_string(ii);
    std::string name = "C" + number;
    testExpression.attachParameterNode(name, C_param_expressions[ii]);
  }

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0; 
  std::vector<std::complex<double> > V = {5.0, 2.0, 3.0};
  solnGroup->setSoln(std::string("A"),V[0]);
  solnGroup->setSoln(std::string("B"),V[1]);
  solnGroup->setSoln(std::string("C"),V[2]);

  std::complex<double> refRes=0.0;
  {
    int currentOrderSize=1;
    int currentMaxIndex=currentOrderSize;
    std::vector<int> orderMinIndex(1,0);
    std::vector<int> orderMaxIndex(1,currentMaxIndex);
    while(currentMaxIndex < C.size())
    {
      currentOrderSize *= V.size();
      orderMinIndex.push_back(currentMaxIndex);
      currentMaxIndex += currentOrderSize;
      orderMaxIndex.push_back(currentMaxIndex);
    }

    refRes=C[0];
    int maxOrderPlus1 = orderMaxIndex.size();
    for (int iorder=1;iorder<maxOrderPlus1;iorder++)
    {
      int localIndex=0;
      for(int ii=orderMinIndex[iorder];ii<orderMaxIndex[iorder];++ii,++localIndex)
      {
        if(ii>=C.size()) break;
        std::complex<double> newval = C[ii];
        std::vector<int> coords;
        setupCoordinates(iorder, localIndex, V.size(), coords);
        for (int ie=0;ie<coords.size();ie++) { newval *= V[coords[ie]]; }
        refRes += newval;
      }
    }
  }

  testExpression.evaluateFunction(result);
  EXPECT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);
  EXPECT_EQ( result, refRes);
  assignExpression.evaluateFunction(result);
  EXPECT_EQ( result, refRes);
}


// four input, fourth order POLY test.  
//
// This function builds the POLY expression string, assigns values to parameters, 
// and evaluates the reference polynomial via arrays, which are sized based on 
// the number of variables and the order.   So this test can be modified to have 
// different numbers for maxOrder and numVars, and it will still be a valid test.
//
// This will have 341 total coefs, based on the math below
//
//ii=0 numCoefs=1   4^0 = +1
//ii=1 numCoefs=5  +4^1 = +4
//ii=2 numCoefs=21 +4^2 = +16 
//ii=3 numCoefs=85 +4^3 = +64
//ii=4 numCoefs=341 +4^4 = +256
TEST ( ComplexParserpolyTest, test14)
{
  Teuchos::RCP<solnExpressionGroup> solnGroup = Teuchos::rcp(new solnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;

  int maxOrder=4, numVars=4;
  //int maxOrder=1, numVars=4; // easy to debug, but not the full test
  std::string exprString = std::string("POLY(");
  std::string numberVars = std::to_string(numVars);
  exprString += numberVars + ")";
  for (int ii=0;ii<numVars;ii++)
  {
    char letter = 'A' + ii;
    exprString += " V(" + std::string(1,letter) + ")";
  }

  int numCoefs = 0;
  for (int ii=0;ii<=maxOrder;ii++) { numCoefs += std::pow(numVars,ii); }

  std::vector<std::complex<double> > C(numCoefs);
  std::vector < Teuchos::RCP<Xyce::Util::newExpression> > C_param_expressions;
  for (int ii=0;ii<numCoefs;ii++)
  {
    std::string name = " C" + std::to_string(ii);
    exprString += name;
    double realCvalue = 10.0-((static_cast<double>(ii+1))*0.05);
    C[ii] =  std::complex<double> (realCvalue , 0.0);
    std::string paramNumber = std::to_string(realCvalue);
    C_param_expressions.push_back( Teuchos::rcp(new Xyce::Util::newExpression (paramNumber, testGroup)) );
  }
  Xyce::Util::newExpression testExpression(exprString, testGroup);
  testExpression.lexAndParseExpression();

  for (int ii=0;ii<C_param_expressions.size();ii++) { C_param_expressions[ii]->lexAndParseExpression(); }

  for (int ii=0;ii<C_param_expressions.size();ii++) 
  { 
    std::string name = "C" + std::to_string(ii);
    testExpression.attachParameterNode(name, C_param_expressions[ii]);
  }

  //testExpression.dumpParseTree(std::cout);
  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0; 
  std::vector<std::complex<double> > V(numVars);
  for(int ii=0;ii<numVars;ii++) 
  { 
    V[ii] = 1.0+static_cast<double>(ii);
    char letter = 'A' + ii;
    solnGroup->setSoln(std::string(1,letter),V[ii]);
  }

  std::complex<double> refRes=0.0;
  {
    int currentOrderSize=1;
    int currentMaxIndex=currentOrderSize;
    std::vector<int> orderMinIndex(1,0);
    std::vector<int> orderMaxIndex(1,currentMaxIndex);
    while(currentMaxIndex < C.size())
    {
      currentOrderSize *= V.size();
      orderMinIndex.push_back(currentMaxIndex);
      currentMaxIndex += currentOrderSize;
      orderMaxIndex.push_back(currentMaxIndex);
    }

    refRes=C[0];
    int maxOrderPlus1 = orderMaxIndex.size();
    for (int iorder=1;iorder<maxOrderPlus1;iorder++)
    {
      int localIndex=0;
      for(int ii=orderMinIndex[iorder];ii<orderMaxIndex[iorder];++ii,++localIndex)
      {
        if(ii>=C.size()) break;
        std::complex<double> newval = C[ii];
        std::vector<int> coords;
        setupCoordinates(iorder, localIndex, V.size(), coords);
        for (int ie=0;ie<coords.size();ie++) { newval *= V[coords[ie]]; }
        refRes += newval;
      }
    }
  }


  testExpression.evaluateFunction(result);
  ASSERT_EQ( result, refRes);
  copyExpression.evaluateFunction(result);
  ASSERT_EQ( result, refRes);
  assignExpression.evaluateFunction(result);
  ASSERT_EQ( result, refRes);
}

#if 0
TEST ( ComplexParserNestedFuncTest, func_cir)
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

  std::complex<double>  result;
  // x,y,z = 4,2,3
  // refresult = diff(y,z)-4+fy(z,x)**2;
  // refresult = diff(2,3)-4+fy(3,4)**2;
  // refresult = (2-3)-4+4*4;
  // refresult = -5+16;
  // refresult = 11;
  std::complex<double>  refresult = 11;
  testExpression.evaluateFunction(result);   EXPECT_EQ( result, refresult );
  copyExpression.evaluateFunction(result);   EXPECT_EQ( result, refresult );
  assignExpression.evaluateFunction(result); EXPECT_EQ( result, refresult );

  OUTPUT_MACRO(ComplexParserNestedFuncTest, func_cir)
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
TEST ( ComplexParserNestedFuncTest, func_cir_newResolution)
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

  std::complex<double>  result;
  // x,y,z = 4,2,3
  // refresult = diff(y,z)-4+fy(z,x)**2;
  // refresult = diff(2,3)-4+fy(3,4)**2;
  // refresult = (2-3)-4+4*4;
  // refresult = -5+16;
  // refresult = 11;
  std::complex<double>  refresult = 11;
  testExpression->evaluateFunction(result);   EXPECT_FLOAT_EQ( std::real(result), std::real(refresult) );

  OUTPUT_MACRO3(ComplexParserNestedFuncTest, func_cir)
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
TEST ( ComplexParserNestedFuncTest, func_cir_newResolution2)
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

  std::complex<double>  result;
  // x,y,z = 4,2,3
  // refresult = diff(y,z)-4+fy(z,x)**2;
  // refresult = diff(2,3)-4+fy(3,4)**2;
  // refresult = (2-3)-4+4*4;
  // refresult = -5+16;
  // refresult = 11;
  std::complex<double>  refresult = 11;
  testExpression->evaluateFunction(result);   EXPECT_FLOAT_EQ( std::real(result), std::real(refresult) );

  OUTPUT_MACRO3(ComplexParserNestedFuncTest, func_cir)
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
TEST ( ComplexParserNestedFuncTest, 200nest_no_deriv)
{
  Teuchos::RCP<solnAndFuncExpressionGroup> funcGroup = Teuchos::rcp(new solnAndFuncExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  int numFuncs=200;

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

  std::complex<double>  result=0.0, Aval=5.0;
  funcGroup->setSoln(std::string("A"),Aval);

  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refderivs = {2.0};
  std::complex<double>  refresult = 10.0;
  testExpression->evaluateFunction(result);   EXPECT_EQ( result, refresult );

  OUTPUT_MACRO3 ( ComplexParserNestedFuncTest, 200nest_no_deriv)
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
TEST ( ComplexParserNestedFuncTest, 200nest_with_deriv)
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

  std::complex<double>  result=0.0, Aval=5.0;
  funcGroup->setSoln(std::string("A"),Aval);

  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refderivs = {2.0};
  std::complex<double>  refresult = 10.0;

  // the dx function (called under evaluate) has some bottlenecks
  testExpression->evaluate(result,derivs);   
  EXPECT_EQ( result, refresult );
  EXPECT_EQ( derivs, refderivs);

  OUTPUT_MACRO3(ComplexParserNestedFuncTest, 200nest_with_deriv)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserNestedGlobalParamTest, 1000nest_no_deriv)
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

  std::complex<double>  result=0.0, Aval=5.0;
  funcGroup->setSoln(std::string("A"),Aval);

  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refderivs = {2.0};
  std::complex<double>  refresult = 10.0;
  testExpression->evaluateFunction(result);   EXPECT_EQ( result, refresult );

  OUTPUT_MACRO3(ComplexParserNestedGlobalParamTest, 1000nest_no_deriv)
}

//-------------------------------------------------------------------------------
// SDT tests
//-------------------------------------------------------------------------------
#if 1
#define NUM_SDT_STEPS 101
#define NUM_SDT_STEPS2 1001
#else
#define NUM_SDT_STEPS 11
#define NUM_SDT_STEPS2 11
#endif
TEST ( ComplexParserIntegralTest, sdt1)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  Xyce::Util::newExpression testExpression(std::string("SDT(V(A))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result = 0.0, refRes = 0.0;
  double time=0.0, finalTime=1.0;
  int numSteps = NUM_SDT_STEPS;
  double dt = finalTime/(numSteps-1);
  for (int ii=0;ii<numSteps;ii++)
  {
    std::complex<double> Aval=time;
    sdtGroup->setSoln(std::string("A"),Aval);
    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt);
    refRes = time*time*0.5;

    testExpression.evaluateFunction(result);   ASSERT_NEAR( std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    copyExpression.evaluateFunction(result);   ASSERT_NEAR( std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    assignExpression.evaluateFunction(result); ASSERT_NEAR( std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));

    time += dt;

    // only one of these "clear" calls is necessary. but test calling all to make sure nothing broken  
    // Also, the best way is using the "static" function call:
    //     Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    // but that is sufficiently tested elsewhere
    //
    testExpression.clearProcessSuccessfulTimeStepMap(); 
    copyExpression.clearProcessSuccessfulTimeStepMap();
    assignExpression.clearProcessSuccessfulTimeStepMap();

    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserIntegralTest, sdt1)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserIntegralTest, sdt2)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  Xyce::Util::newExpression testExpression(std::string("SDT (3.0*cos(time))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::vector<std::complex<double> > derivs;
  std::complex<double> result = 0.0, refRes = 0.0;
  double time=0.0, finalTime=1.0;
  int numSteps = NUM_SDT_STEPS2;
  double dt = finalTime/(numSteps-1);
  for (int ii=0;ii<numSteps;ii++)
  {
    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt);
    refRes = 3.0*std::sin(time);

    testExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_EQ(derivs.size(),0);

    copyExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_EQ(derivs.size(),0);

    assignExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result),std::real(refRes),std::abs(1.0e-4*std::real(result)));
    ASSERT_EQ(derivs.size(),0);

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserIntegralTest, sdt2)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserIntegralTest, sdt3)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  Xyce::Util::newExpression testExpression(std::string("SDT(v(a))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refDerivs;
  std::complex<double> result = 0.0, refRes = 0.0;
  double time=0.0, finalTime=1.0;
  int numSteps = NUM_SDT_STEPS;
  double dt = finalTime/(numSteps-1);
  refDerivs.resize(1,0.5*dt);
  for (int ii=0;ii<numSteps;ii++)
  {
    std::complex<double> Aval=time;
    sdtGroup->setSoln(std::string("A"),Aval);
    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt); 

    if (ii!=0) { refDerivs[0] = 0.5*dt; }
    else       { refDerivs[0] = 0.0; }

    refRes = time*time*0.5;

    testExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));
    copyExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));
    assignExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserIntegralTest, sdt3)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserIntegralTest, sdt4)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  Xyce::Util::newExpression testExpression(std::string("SDT(2.0*v(a))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refDerivs;
  std::complex<double> result = 0.0, refRes = 0.0; 
  double time=0.0, finalTime=1.0;
  int numSteps = NUM_SDT_STEPS;
  double dt = finalTime/(numSteps-1);
  refDerivs.resize(1,0.5*dt);

  std::complex<double> Aval=1.0;
  sdtGroup->setSoln(std::string("A"),Aval);

  for (int ii=0;ii<numSteps;ii++)
  {
    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt); 

    if (ii!=0) { refDerivs[0] = dt;  }
    else       { refDerivs[0] = 0.0; }

    refRes = 2.0*time;

    testExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));

    copyExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));

    assignExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserIntegralTest, sdt4)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserIntegralTest, sdt5)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  Xyce::Util::newExpression testExpression(std::string("F1(V(A))"), testGroup);
  testExpression.lexAndParseExpression();

  // .func F1(A) {sdt(A)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(A)");
  std::string rhs=std::string("sdt(A)");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;
  //testExpression.dumpParseTree(std::cout);

  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refDerivs;
  std::complex<double> result = 0.0, refRes = 0.0;
  double time=0.0, finalTime=1.0;
  int numSteps = NUM_SDT_STEPS;
  double dt = finalTime/(numSteps-1);

  refDerivs.push_back(0.5*dt);

  for (int ii=0;ii<numSteps;ii++)
  {
    std::complex<double> Aval=time;
    sdtGroup->setSoln(std::string("A"),Aval);
    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt);
    refRes = time*time*0.5;

    if (ii!=0) { refDerivs[0] = 0.5*dt;  }
    else       { refDerivs[0] = 0.0; }

    testExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));

    copyExpression.evaluate(result,derivs);
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));

    assignExpression.evaluate(result,derivs);
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserIntegralTest, sdt5)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserIntegralTest, sdt6)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  Xyce::Util::newExpression testExpression(std::string("F1(V(A))+F1(V(A))"), testGroup);
  testExpression.lexAndParseExpression();

  // .func F1(A) {sdt(A)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(A)");
  std::string rhs=std::string("sdt(A)");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refDerivs;
  std::complex<double> result = 0.0, refRes = 0.0;
  double time=0.0, finalTime=1.0;
  int numSteps = NUM_SDT_STEPS;
  double dt = finalTime/(numSteps-1);
  refDerivs.push_back(dt);

  for (int ii=0;ii<numSteps;ii++)
  {
    double Aval=time;
    sdtGroup->setSoln(std::string("A"),Aval);
    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt);

    if (ii!=0) { refDerivs[0] = dt;  }
    else       { refDerivs[0] = 0.0; }

    refRes = time*time;

    testExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));

    copyExpression.evaluate(result,derivs);
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));

    assignExpression.evaluate(result,derivs);
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserIntegralTest, sdt6)
}

//-------------------------------------------------------------------------------
// this point of this test is to make sure that when SDT is inside of a .FUNC, 
// and that .FUNC is called more than once, that the integral informmation doesn't 
// get mangled.  Each call to the .FUNC should carry its own unique SDT state.
TEST ( ComplexParserIntegralTest, sdt7)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  Xyce::Util::newExpression testExpression(std::string("F1(V(A))+F1(V(B))"), testGroup);
  testExpression.lexAndParseExpression();

  // .func F1(A) {sdt(A)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(A)");
  std::string rhs=std::string("sdt(A)");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  //testExpression.dumpParseTree(std::cout);

  std::vector<std::complex<double> > derivs(2);
  std::vector<std::complex<double> > refDerivs(2);
  std::complex<double> result = 0.0, refRes = 0.0;
  double time=0.0, finalTime=1.0;
  int numSteps = NUM_SDT_STEPS2;
  double dt = finalTime/(numSteps-1);

  for (int ii=0;ii<numSteps;ii++)
  {
    double Aval=time;
    double Bval=3.0*std::cos(time);
    sdtGroup->setSoln(std::string("A"),Aval);
    sdtGroup->setSoln(std::string("B"),Bval);
    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt);

    double Aint = time*time*0.5;
    double Bint = 3.0*std::sin(time);
    refRes = Aint + Bint;

    if (ii!=0)
    {
      refDerivs[0] = 0.5*dt;
      refDerivs[1] = 0.5*dt;
    }
    else
    {
      refDerivs[0] = 0.0;
      refDerivs[1] = 0.0;
    }

    testExpression.evaluate(result,derivs);
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));
    ASSERT_NEAR(std::real(derivs[1]), std::real(refDerivs[1]), std::abs(1.0e-4*std::real(derivs[1])));

    copyExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));
    ASSERT_NEAR(std::real(derivs[1]), std::real(refDerivs[1]), std::abs(1.0e-4*std::real(derivs[1])));

    assignExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));
    ASSERT_NEAR(std::real(derivs[1]), std::real(refDerivs[1]), std::abs(1.0e-4*std::real(derivs[1])));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserIntegralTest, sdt7)
}

//-------------------------------------------------------------------------------
// this test is similar to sdt7, except that the SDT operators 
// are behind 2 layers of funcs instead of 1.
TEST ( ComplexParserIntegralTest, sdt8)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  Xyce::Util::newExpression testExpression(std::string("F1(V(A))+F1(V(B))"), testGroup);
  testExpression.lexAndParseExpression();

  // .func F1(A) {f2(a)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(A)");
  std::string rhs=std::string("f2(a)");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  // .func F2(A) {sdt(a)}
  std::string f2Name;
  Teuchos::RCP<Xyce::Util::newExpression> f2Expression;
  lhs=std::string("F2(A)");
  rhs=std::string("sdt(a)");
  createFunc(lhs,rhs,testGroup, f2Name,f2Expression);

  f1Expression->attachFunctionNode(f2Name, f2Expression);
  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  //testExpression.dumpParseTree(std::cout);

  std::vector<std::complex<double> > derivs(2);
  std::vector<std::complex<double> > refDerivs(2);
  std::complex<double> result = 0.0, refRes = 0.0;
  double time=0.0, finalTime=1.0;
  int numSteps = NUM_SDT_STEPS2;
  double dt = finalTime/(numSteps-1);
  for (int ii=0;ii<numSteps;ii++)
  {
    double Aval=time;
    double Bval=3.0*std::cos(time);
    sdtGroup->setSoln(std::string("A"),Aval);
    sdtGroup->setSoln(std::string("B"),Bval);
    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt);

    double Aint = time*time*0.5;
    double Bint = 3.0*std::sin(time);
    refRes = Aint + Bint;

    if (ii!=0)
    {
      refDerivs[0] = 0.5*dt;
      refDerivs[1] = 0.5*dt;
    }
    else
    {
      refDerivs[0] = 0.0;
      refDerivs[1] = 0.0;
    }

    testExpression.evaluate(result,derivs);
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));
    ASSERT_NEAR(std::real(derivs[1]), std::real(refDerivs[1]), std::abs(1.0e-4*std::real(derivs[1])));

    copyExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));
    ASSERT_NEAR(std::real(derivs[1]), std::real(refDerivs[1]), std::abs(1.0e-4*std::real(derivs[1])));

    assignExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));
    ASSERT_NEAR(std::real(derivs[1]), std::real(refDerivs[1]), std::abs(1.0e-4*std::real(derivs[1])));


    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserIntegralTest, sdt8)
}

//-------------------------------------------------------------------------------
// this test is similar to sdt7, except that the SDT operators 
// are behind 3 layers of funcs instead of 1.
TEST ( ComplexParserIntegralTest, sdt9)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  Xyce::Util::newExpression testExpression(std::string("F1(V(A))+F1(V(B))"), testGroup);
  testExpression.lexAndParseExpression();

  // .func F1(A) {f2(a)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(A)");
  std::string rhs=std::string("f2(a)");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  // .func F2(A) {f3(a)}
  std::string f2Name;
  Teuchos::RCP<Xyce::Util::newExpression> f2Expression;
  lhs=std::string("F2(A)");
  rhs=std::string("f3(a)");
  createFunc(lhs,rhs,testGroup, f2Name,f2Expression);

  // .func F3(A) {sdt(A)}
  std::string f3Name;
  Teuchos::RCP<Xyce::Util::newExpression> f3Expression;
  lhs=std::string("F3(A)");
  rhs=std::string("sdt(A)");
  createFunc(lhs,rhs,testGroup, f3Name,f3Expression);

  f2Expression->attachFunctionNode(f3Name, f3Expression);
  f1Expression->attachFunctionNode(f2Name, f2Expression);
  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  //testExpression.dumpParseTree(std::cout);

  std::vector<std::complex<double> > derivs(2);
  std::vector<std::complex<double> > refDerivs(2);
  std::complex<double> result = 0.0, refRes = 0.0;
  double time=0.0, finalTime=1.0;
  int numSteps = NUM_SDT_STEPS2;
  double dt = finalTime/(numSteps-1);
  for (int ii=0;ii<numSteps;ii++)
  {
    double Aval=time;
    double Bval=3.0*std::cos(time);
    sdtGroup->setSoln(std::string("A"),Aval);
    sdtGroup->setSoln(std::string("B"),Bval);
    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt);

    double Aint = time*time*0.5;
    double Bint = 3.0*std::sin(time);
    refRes = Aint + Bint;

    if (ii!=0)
    {
      refDerivs[0] = 0.5*dt;
      refDerivs[1] = 0.5*dt;
    }
    else
    {
      refDerivs[0] = 0.0;
      refDerivs[1] = 0.0;
    }

    testExpression.evaluate(result,derivs);
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));
    ASSERT_NEAR(std::real(derivs[1]), std::real(refDerivs[1]), std::abs(1.0e-4*std::real(derivs[1])));

    copyExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));
    ASSERT_NEAR(std::real(derivs[1]), std::real(refDerivs[1]), std::abs(1.0e-4*std::real(derivs[1])));

    assignExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));
    ASSERT_NEAR(std::real(derivs[1]), std::real(refDerivs[1]), std::abs(1.0e-4*std::real(derivs[1])));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserIntegralTest, sdt9)
}

//-------------------------------------------------------------------------------
// this test is similar to sdt7, except that the SDT operators 
// are behind 3 layers of funcs instead of 1.
// Also, the f1 layer calls sdt directly, which doubles the size of the answer
TEST ( ComplexParserIntegralTest, sdt10)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  Xyce::Util::newExpression testExpression(std::string("F1(V(A))+F1(V(B))"), testGroup);
  testExpression.lexAndParseExpression();

  // .func F1(A) {f2(a)+sdt(a)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(A)");
  std::string rhs=std::string("f2(a)+sdt(a)");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  // .func F2(A) {f3(a)}
  std::string f2Name;
  Teuchos::RCP<Xyce::Util::newExpression> f2Expression;
  lhs=std::string("F2(A)");
  rhs=std::string("f3(a)");
  createFunc(lhs,rhs,testGroup, f2Name,f2Expression);

  // .func F3(A) {sdt(A)}
  std::string f3Name;
  Teuchos::RCP<Xyce::Util::newExpression> f3Expression;
  lhs=std::string("F3(A)");
  rhs=std::string("sdt(A)");
  createFunc(lhs,rhs,testGroup, f3Name,f3Expression);

  f2Expression->attachFunctionNode(f3Name, f3Expression);
  f1Expression->attachFunctionNode(f2Name, f2Expression);
  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  //testExpression.dumpParseTree(std::cout);

  std::vector<std::complex<double> > derivs(2);
  std::vector<std::complex<double> > refDerivs(2);
  std::complex<double> result = 0.0, refRes = 0.0;
  double time=0.0, finalTime=1.0;
  int numSteps = NUM_SDT_STEPS2;
  double dt = finalTime/(numSteps-1);
  for (int ii=0;ii<numSteps;ii++)
  {
    double Aval=time;
    double Bval=3.0*std::cos(time);
    sdtGroup->setSoln(std::string("A"),Aval);
    sdtGroup->setSoln(std::string("B"),Bval);
    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt);

    double Aint = time*time*0.5;
    double Bint = 3.0*std::sin(time);
    refRes = Aint + Bint;
    refRes *= 2.0;

    if (ii!=0)
    {
      refDerivs[0] = dt;
      refDerivs[1] = dt;
    }
    else
    {
      refDerivs[0] = 0.0;
      refDerivs[1] = 0.0;
    }

    testExpression.evaluate(result,derivs);
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));
    ASSERT_NEAR(std::real(derivs[1]), std::real(refDerivs[1]), std::abs(1.0e-4*std::real(derivs[1])));

    copyExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));
    ASSERT_NEAR(std::real(derivs[1]), std::real(refDerivs[1]), std::abs(1.0e-4*std::real(derivs[1])));

    assignExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));
    ASSERT_NEAR(std::real(derivs[1]), std::real(refDerivs[1]), std::abs(1.0e-4*std::real(derivs[1])));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserIntegralTest, sdt10)
}

//-------------------------------------------------------------------------------
// this test is similar to sdt7, except that the SDT operators 
// are behind 3 layers of funcs instead of 1.
// Also, the f1 layer calls sdt directly, which doubles the size of the answer
TEST ( ComplexParserIntegralTest, sdt11)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  Xyce::Util::newExpression testExpression(std::string("F1(V(A))+F1(V(B))"), testGroup);
  testExpression.lexAndParseExpression();

  // .func F1(A) {f2(a)+sdt(a)+sdt(v(b))}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(A)");
  std::string rhs=std::string("f2(a)+sdt(a)+sdt(v(b))");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  // .func F2(A) {f3(a)*2.0+sdt(a)}
  std::string f2Name;
  Teuchos::RCP<Xyce::Util::newExpression> f2Expression;
  lhs=std::string("F2(A)");
  rhs=std::string("f3(a)*2.0+sdt(a)");
  createFunc(lhs,rhs,testGroup, f2Name,f2Expression);

  // .func F3(A) {sdt(A)}
  std::string f3Name;
  Teuchos::RCP<Xyce::Util::newExpression> f3Expression;
  lhs=std::string("F3(A)");
  rhs=std::string("sdt(a)");
  createFunc(lhs,rhs,testGroup, f3Name,f3Expression);

  f2Expression->attachFunctionNode(f3Name, f3Expression);
  f1Expression->attachFunctionNode(f2Name, f2Expression);
  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  //testExpression.dumpParseTree(std::cout);

  std::vector<std::complex<double> > derivs(2);
  std::vector<std::complex<double> > refDerivs(2);
  std::complex<double> result = 0.0, refRes = 0.0;
  double time=0.0, finalTime=1.0;
  int numSteps = NUM_SDT_STEPS2;
  double dt = finalTime/(numSteps-1);
  for (int ii=0;ii<numSteps;ii++)
  {
    double Aval=time;
    double Bval=3.0*std::cos(time);
    sdtGroup->setSoln(std::string("A"),Aval);
    sdtGroup->setSoln(std::string("B"),Bval);
    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt);

    double Aint = time*time*0.5;
    double Bint = 3.0*std::sin(time);
    refRes = 4.0*Aint + 6.0*Bint;

    if (ii!=0)
    {
      refDerivs[0] = 2.0*dt;
      refDerivs[1] = 3.0*dt;
    }
    else
    {
      refDerivs[0] = 0.0;
      refDerivs[1] = 0.0;
    }

    testExpression.evaluate(result,derivs);
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));
    ASSERT_NEAR(std::real(derivs[1]), std::real(refDerivs[1]), std::abs(1.0e-4*std::real(derivs[1])));

    copyExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));
    ASSERT_NEAR(std::real(derivs[1]), std::real(refDerivs[1]), std::abs(1.0e-4*std::real(derivs[1])));

    assignExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));
    ASSERT_NEAR(std::real(derivs[1]), std::real(refDerivs[1]), std::abs(1.0e-4*std::real(derivs[1])));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserIntegralTest, sdt11)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserIntegralTest, sdt12)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  Xyce::Util::newExpression testExpression(std::string("SDT (F1(V(A)))"), testGroup);
  testExpression.lexAndParseExpression();

  // .func F1(B) {(3.0*cos(B))}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(B)"), rhs=std::string("(3.0*cos(B))");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::vector<std::complex<double> > derivs(1);
  std::vector<std::complex<double> > refDerivs(1);
  std::complex<double> result = 0.0, refRes = 0.0;
  double time=0.0, finalTime=1.0;
  int numSteps = NUM_SDT_STEPS2;
  double dt = finalTime/(numSteps-1);
  for (int ii=0;ii<numSteps;ii++)
  {
    double Aval=time;
    sdtGroup->setSoln(std::string("A"),Aval);
    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt);

    refRes = 3.0*std::sin(time);

    if (ii!=0)
    {
      refDerivs[0] = -3.0*0.5*dt*sin(Aval);
    }
    else
    {
      refDerivs[0] = 0.0;
    }

    testExpression.evaluate(result,derivs);
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));

    copyExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));

    assignExpression.evaluate(result,derivs);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    ASSERT_NEAR(std::real(derivs[0]), std::real(refDerivs[0]), std::abs(1.0e-4*std::real(derivs[0])));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserIntegralTest, sdt12)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserIntegralTest, sdt13)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  Xyce::Util::newExpression testExpression(std::string("SDT (F1(V(A))*F1(V(A)))"), testGroup);
  testExpression.lexAndParseExpression();

  // .func F1(B) {(3.0*cos(B))}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(B)"), rhs=std::string("(3.0*cos(B))");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result = 0.0, refRes = 0.0; 
  double time=0.0, finalTime=1.0;
  int numSteps = NUM_SDT_STEPS2;
  double dt = finalTime/(numSteps-1);
  for (int ii=0;ii<numSteps;ii++)
  {
    std::complex<double> Aval=time;
    sdtGroup->setSoln(std::string("A"),Aval);
    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt);
    refRes = 9.0*(0.5*time  + 0.25 * std::sin(2.0*time));
 
    testExpression.evaluateFunction(result);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));

    copyExpression.evaluateFunction(result);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));

    assignExpression.evaluateFunction(result);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserIntegralTest, sdt13)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserIntegralTest, sdt14)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  Xyce::Util::newExpression testExpression(std::string("SDT (F1(V(A))*F1(V(A)))"), testGroup);
  testExpression.lexAndParseExpression();

  // .func F1(A) {f2(A)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(A)"), rhs=std::string("f2(a)");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  // .func F2(B) {(3.0*cos(B))}
  std::string f2Name;
  Teuchos::RCP<Xyce::Util::newExpression> f2Expression;
  lhs=std::string("F2(B)"); rhs=std::string("(3.0*cos(B))");
  createFunc(lhs,rhs,testGroup, f2Name,f2Expression);

  f1Expression->attachFunctionNode(f2Name, f2Expression);
  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result = 0.0, refRes = 0.0; 
  double time=0.0, finalTime=1.0;
  int numSteps = NUM_SDT_STEPS2;
  double dt = finalTime/(numSteps-1);

  for (int ii=0;ii<numSteps;ii++)
  {
    double Aval=time;
    sdtGroup->setSoln(std::string("A"),Aval);
    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt);
    refRes = 9.0*(0.5*time  + 0.25 * std::sin(2.0*time));
 
    testExpression.evaluateFunction(result);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));

    copyExpression.evaluateFunction(result);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));

    assignExpression.evaluateFunction(result);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserIntegralTest, sdt14)
}


//-------------------------------------------------------------------------------
// The point of this test is to make sure that one function call doesn't 
// corrupt another one.
//
// testExpression and testExpression3 are functionally identical.
//
// But, testExpression2 calls the same function as testExpression 
// but using different inputs.  If state is managed carefully, then this should
// not be a problem.
//
// In this case, the argument to sdt is a function argument parameter.  So, a 
// relatively simple scenario
//-------------------------------------------------------------------------------
TEST ( ComplexParserIntegralTest, sdt15)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  Xyce::Util::newExpression testExpression(std::string("f1(V(A))"), testGroup);
  testExpression.lexAndParseExpression();
  Xyce::Util::newExpression testExpression2(std::string("f1(V(B))"), testGroup);
  testExpression2.lexAndParseExpression();
  Xyce::Util::newExpression testExpression3(std::string("f3(V(A))"), testGroup);
  testExpression3.lexAndParseExpression();

  // .func F1(A) {f2(A)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(A)"), rhs=std::string("f2(a)");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  // .func F2(B) {sdt(B)}
  std::string f2Name;
  Teuchos::RCP<Xyce::Util::newExpression> f2Expression;
  lhs=std::string("F2(B)"); rhs=std::string("sdt(B)");
  createFunc(lhs,rhs,testGroup, f2Name,f2Expression);

  // .func F3(A) {f4(A)}
  std::string f3Name;
  Teuchos::RCP<Xyce::Util::newExpression> f3Expression;
  lhs=std::string("F3(A)"); rhs=std::string("f4(a)");
  createFunc(lhs,rhs,testGroup, f3Name,f3Expression);

  // .func F4(B) {sdt(B)}
  std::string f4Name;
  Teuchos::RCP<Xyce::Util::newExpression> f4Expression;
  lhs=std::string("F4(B)"); rhs=std::string("sdt(B)");
  createFunc(lhs,rhs,testGroup, f4Name,f4Expression);


  f1Expression->attachFunctionNode(f2Name, f2Expression);
  testExpression.attachFunctionNode(f1Name, f1Expression);
  testExpression2.attachFunctionNode(f1Name, f1Expression);


  f3Expression->attachFunctionNode(f4Name, f4Expression);
  testExpression3.attachFunctionNode(f3Name, f3Expression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result = 0.0, refRes = 0.0, refRes2 = 0.0; 
  double time=0.0;
  double finalTime=1.25000000e-07;
  int numSteps = NUM_SDT_STEPS2;
  double dt = finalTime/(numSteps-1);
  for (int ii=0;ii<numSteps;ii++)
  {
    double mpi = M_PI;
    double freq = 4000e3;
    double Aval= 1.0e3 * std::sin (2.0*mpi*((std::real(freq))*std::real(time) )) ;
    sdtGroup->setSoln(std::string("A"),Aval);

    double Bval = time;
    sdtGroup->setSoln(std::string("B"),Bval);

    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt);
 
    testExpression.evaluateFunction(result);   
    testExpression2.evaluateFunction(refRes);   
    testExpression3.evaluateFunction(refRes2);
    ASSERT_EQ( result, refRes2);

    copyExpression.evaluateFunction(result);   
    ASSERT_EQ( result, refRes2);

    assignExpression.evaluateFunction(result);   
    ASSERT_EQ( result, refRes2);

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    testExpression2.processSuccessfulTimeStep();
    testExpression3.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserIntegralTest, sdt15)
}

//-------------------------------------------------------------------------------
// The point of this test is to make sure that one function call doesn't 
// corrupt another one.
//
// testExpression and testExpression3 are functionally identical.
//
// But, testExpression2 calls the same function as testExpression 
// but using different inputs.  If state is managed carefully, then this should
// not be a problem.
//
// This test is similar, but different from sdt15, in that the argument to sdt
// is an expression. (B*B) instead of (B).  This was a little harder to manage.
//-------------------------------------------------------------------------------
TEST ( ComplexParserIntegralTest, sdt16)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  Xyce::Util::newExpression testExpression(std::string("f1(V(A))"), testGroup);
  testExpression.lexAndParseExpression();
  Xyce::Util::newExpression testExpression2(std::string("f1(V(B))"), testGroup);
  testExpression2.lexAndParseExpression();
  Xyce::Util::newExpression testExpression3(std::string("f3(V(A))"), testGroup);
  testExpression3.lexAndParseExpression();

  // .func F1(A) {f2(A)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(A)"), rhs=std::string("f2(a)");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  // .func F2(B) {sdt(B*B)}
  std::string f2Name;
  Teuchos::RCP<Xyce::Util::newExpression> f2Expression;
  lhs=std::string("F2(B)"); rhs=std::string("sdt(B*B)");
  createFunc(lhs,rhs,testGroup, f2Name,f2Expression);

  // .func F3(A) {f4(A)}
  std::string f3Name;
  Teuchos::RCP<Xyce::Util::newExpression> f3Expression;
  lhs=std::string("F3(A)"); rhs=std::string("f4(a)");
  createFunc(lhs,rhs,testGroup, f3Name,f3Expression);

  // .func F4(B) {sdt(B*B)}
  std::string f4Name;
  Teuchos::RCP<Xyce::Util::newExpression> f4Expression;
  lhs=std::string("F4(B)"); rhs=std::string("sdt(B*B)");
  createFunc(lhs,rhs,testGroup, f4Name,f4Expression);

  f1Expression->attachFunctionNode(f2Name, f2Expression);
  testExpression.attachFunctionNode(f1Name, f1Expression);
  testExpression2.attachFunctionNode(f1Name, f1Expression);

  f3Expression->attachFunctionNode(f4Name, f4Expression);
  testExpression3.attachFunctionNode(f3Name, f3Expression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result = 0.0, refRes = 0.0, refRes2 = 0.0; 
  double time = 0.0;
  double finalTime=1.25000000e-07;
  int numSteps = NUM_SDT_STEPS2;
  double dt = finalTime/(numSteps-1);
  for (int ii=0;ii<numSteps;ii++)
  {
    double mpi = M_PI;
    double freq = 4000e3;
    double Aval= 1.0e3 * std::sin (2.0*mpi*((std::real(freq))*std::real(time) )) ;
    sdtGroup->setSoln(std::string("A"),Aval);

    double Bval = time;
    sdtGroup->setSoln(std::string("B"),Bval);

    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt);
 
    testExpression.evaluateFunction(result);   
    testExpression2.evaluateFunction(refRes);   
    testExpression3.evaluateFunction(refRes2);
    
    ASSERT_EQ( result, refRes2);

    copyExpression.evaluateFunction(result);   
    ASSERT_EQ( result, refRes2);

    assignExpression.evaluateFunction(result);   
    ASSERT_EQ( result, refRes2);

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    testExpression2.processSuccessfulTimeStep();
    testExpression3.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserIntegralTest, sdt16)
}

//-------------------------------------------------------------------------------
// The point of this test is to make sure that one function call doesn't 
// corrupt another one.
//
// testExpression and testExpression3 are functionally identical.
//
// But, testExpression2 calls the same function as testExpression 
// but using different inputs.  If state is managed carefully, then this should
// not be a problem.
//
// This is similar, but slightly more complicated than sdt16, in that the 
// functions have 2 arguments rather than 1.  Among other issues, having 
// two arguments introduces the possible bug of having (A,B) be treated 
// the same as (B,A), in terms of  the state key. (the key that determines which
// state to use).
//-------------------------------------------------------------------------------
TEST ( ComplexParserIntegralTest, sdt17)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  Xyce::Util::newExpression testExpression(std::string("f1(V(B),V(A))"), testGroup);
  testExpression.lexAndParseExpression();
  Xyce::Util::newExpression testExpression2(std::string("f1(-V(A),V(C))"), testGroup);
  testExpression2.lexAndParseExpression();
  Xyce::Util::newExpression testExpression3(std::string("f3(V(B),V(A))"), testGroup);
  testExpression3.lexAndParseExpression();

  // .func F1(A) {f2(A)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(A,B)"), rhs=std::string("f2(a,b)");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  // .func F2(B) {sdt(B*B)}
  std::string f2Name;
  Teuchos::RCP<Xyce::Util::newExpression> f2Expression;
  lhs=std::string("F2(A,B)"); rhs=std::string("sdt(A+B*B)");
  createFunc(lhs,rhs,testGroup, f2Name,f2Expression);

  // .func F3(A) {f4(A)}
  std::string f3Name;
  Teuchos::RCP<Xyce::Util::newExpression> f3Expression;
  lhs=std::string("F3(A,B)"); rhs=std::string("f4(a,b)");
  createFunc(lhs,rhs,testGroup, f3Name,f3Expression);

  // .func F4(B) {sdt(B*B)}
  std::string f4Name;
  Teuchos::RCP<Xyce::Util::newExpression> f4Expression;
  lhs=std::string("F4(A,B)"); rhs=std::string("sdt(A+B*B)");
  createFunc(lhs,rhs,testGroup, f4Name,f4Expression);

  f1Expression->attachFunctionNode(f2Name, f2Expression);
  testExpression.attachFunctionNode(f1Name, f1Expression);
  testExpression2.attachFunctionNode(f1Name, f1Expression);

  f3Expression->attachFunctionNode(f4Name, f4Expression);
  testExpression3.attachFunctionNode(f3Name, f3Expression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result = 0.0, refRes = 0.0, refRes2=0.0; 
  double time=0.0;
  double finalTime=1.25000000e-07;
  int numSteps = NUM_SDT_STEPS2;
  double dt = finalTime/(numSteps-1);
  for (int ii=0;ii<numSteps;ii++)
  {
    double mpi = M_PI;
    double freq = 4000e3;
    std::complex<double> Aval= 1.0e3 * std::sin (2.0*mpi*((std::real(freq))*std::real(time) )) ;
    std::complex<double> Bval = time;
    std::complex<double> Cval= 1.0e3 * std::cos (2.0*mpi*((std::real(freq))*std::real(time) )) ;
    sdtGroup->setSoln(std::string("A"),Aval);
    sdtGroup->setSoln(std::string("B"),Bval);
    sdtGroup->setSoln(std::string("C"),Cval);

    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt);
 
    testExpression.evaluateFunction(result);   
    testExpression2.evaluateFunction(refRes);   
    testExpression3.evaluateFunction(refRes2);
    
    ASSERT_EQ( result, refRes2);

    copyExpression.evaluateFunction(result);   
    ASSERT_EQ( result, refRes2);

    assignExpression.evaluateFunction(result);   
    ASSERT_EQ( result, refRes2);

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    testExpression2.processSuccessfulTimeStep();
    testExpression3.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserIntegralTest, sdt17)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserIntegralTest, sdt_100nest_no_deriv)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;

  int numFuncs=100;

  // this expression will use the function f1.
  Xyce::Util::newExpression testExpression(std::string("F1(V(A))+F1(V(B))"), testGroup);
  testExpression.lexAndParseExpression();

  // need to set up vector of funcs
  std::vector< std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > >  dotFuncVector;
  // set up the first func in the list (which will be the final one evaluated)  
 
  // need to set up vector of funcs
  for (int ii=1;ii<numFuncs+1;ii++)
  { 
    // .func Fii(A) {fii+1(A)}
    std::string fii_name_args    = std::string("f") + std::to_string(ii) + std::string("(a)");   // lhs
    std::string fiiP1_name_args  = std::string("f") + std::to_string(ii+1) + std::string("(a)"); // rhs
    std::string fiiName;
    Teuchos::RCP<Xyce::Util::newExpression> fiiExpression;
    createFunc(fii_name_args,fiiP1_name_args,testGroup, fiiName,fiiExpression);

    std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > funcPair(fiiName,fiiExpression);
    dotFuncVector.push_back(funcPair);
  }
  
  { 
    // .func Fii(A) {sdt(a)}
    int ii=numFuncs;
    std::string fiiName;
    Teuchos::RCP<Xyce::Util::newExpression> fiiExpression;
    std::string lhs = std::string("f") + std::to_string(ii+1) + std::string("(a)");
    std::string rhs = std::string("sdt(a)");
    createFunc(lhs,rhs,testGroup, fiiName,fiiExpression);

    std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > funcPair(fiiName,fiiExpression);
    dotFuncVector.push_back(funcPair);
  }

  // do all the attachments
  for (int ii=0;ii<dotFuncVector.size()-1;ii++)
  {
    std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > funcPairP1 = dotFuncVector[ii+1];
    std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > funcPair   = dotFuncVector[ii];

    funcPair.second->attachFunctionNode(funcPairP1.first, funcPairP1.second);
  }

  testExpression.attachFunctionNode(dotFuncVector[0].first, dotFuncVector[0].second);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  //testExpression.dumpParseTree(std::cout);

  std::complex<double> result = 0.0, refRes = 0.0; 
  double time=0.0, finalTime=1.0;
  int numSteps = NUM_SDT_STEPS2;
  double dt = finalTime/(numSteps-1);

  for (int ii=0;ii<numSteps;ii++)
  {
    std::complex<double> Aval=time;
    std::complex<double> Bval=3.0*std::cos(time);
    sdtGroup->setSoln(std::string("A"),Aval);
    sdtGroup->setSoln(std::string("B"),Bval);
    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt);
    refRes = time*time*0.5 + 3.0*std::sin(time);

    testExpression.evaluateFunction(result);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    
    copyExpression.evaluateFunction(result);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    
    assignExpression.evaluateFunction(result);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  //OUTPUT_MACRO3(ComplexParserIntegralTest, sdt_1000nest_no_deriv)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserIntegralTest, sdt_100nest_no_deriv2)
{
  Teuchos::RCP<sdtExpressionGroup> sdtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = sdtGroup;


  // this expression will use the function f1.
  Xyce::Util::newExpression testExpression(std::string("SDT(F1(V(A))+F1(V(B)))"), testGroup);
  testExpression.lexAndParseExpression();

  // need to set up vector of funcs
  std::vector< std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > >  dotFuncVector;
 
  // need to set up vector of funcs
  int numFuncs=100;
  for (int ii=1;ii<numFuncs+1;ii++)
  { 
    // .func Fii(A) {fii+1(A)}
    std::string fii_name_args    = std::string("f") + std::to_string(ii) + std::string("(a)");   // lhs
    std::string fiiP1_name_args  = std::string("f") + std::to_string(ii+1) + std::string("(a)"); // rhs
    std::string fiiName;
    Teuchos::RCP<Xyce::Util::newExpression> fiiExpression;
    createFunc(fii_name_args,fiiP1_name_args,testGroup, fiiName,fiiExpression);

    std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > funcPair(fiiName,fiiExpression);
    dotFuncVector.push_back(funcPair);
  }
  
  { 
    // .func Fii(A) {a*2.0}
    int ii=numFuncs;
    std::string fiiName;
    Teuchos::RCP<Xyce::Util::newExpression> fiiExpression;
    std::string lhs = std::string("f") + std::to_string(ii+1) + std::string("(a)");
    std::string rhs = std::string("a*2.0");
    createFunc(lhs,rhs,testGroup, fiiName,fiiExpression);

    std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > funcPair(fiiName,fiiExpression);
    dotFuncVector.push_back(funcPair);
  }

  // do all the attachments
  for (int ii=0;ii<dotFuncVector.size()-1;ii++)
  {
    std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > funcPairP1 = dotFuncVector[ii+1];
    std::pair<std::string, Teuchos::RCP<Xyce::Util::newExpression> > funcPair   = dotFuncVector[ii];

    funcPair.second->attachFunctionNode(funcPairP1.first, funcPairP1.second);
  }

  testExpression.attachFunctionNode(dotFuncVector[0].first, dotFuncVector[0].second);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  //testExpression.dumpParseTree(std::cout);

  std::complex<double> result = 0.0, refRes = 0.0; 
  double time=0.0, finalTime=1.0;
  int numSteps = NUM_SDT_STEPS2;
  double dt = finalTime/(numSteps-1);

  for (int ii=0;ii<numSteps;ii++)
  {
    std::complex<double> Aval=time;
    std::complex<double> Bval=3.0*std::cos(time);
    sdtGroup->setSoln(std::string("A"),Aval);
    sdtGroup->setSoln(std::string("B"),Bval);
    sdtGroup->setTime(time);
    sdtGroup->setStepNumber(ii);
    sdtGroup->setTimeStep(dt);

    refRes = 2.0*(time*time*0.5 + 3.0*std::sin(time));

    testExpression.evaluateFunction(result);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    
    copyExpression.evaluateFunction(result);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));
    
    assignExpression.evaluateFunction(result);   
    ASSERT_NEAR(std::real(result), std::real(refRes), std::abs(1.0e-4*std::real(result)));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  //OUTPUT_MACRO3(ComplexParserIntegralTest, sdt_1000nest_no_deriv)
}


//-------------------------------------------------------------------------------
// DDT tests
//-------------------------------------------------------------------------------
TEST ( ComplexParserDerivativeTest, ddt1)
{
  Teuchos::RCP<sdtExpressionGroup> ddtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = ddtGroup;

  Xyce::Util::newExpression testExpression(std::string("DDT(V(A))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result = 0.0;
  std::complex<double> refRes = 0.0;

  double time=0.0;
  double finalTime=1.0;

  int numSteps = 101;
  double dt = finalTime/(numSteps-1);

  for (int ii=0;ii<numSteps;ii++)
  {
    std::complex<double> Aval=time;
    ddtGroup->setSoln(std::string("A"),Aval);
    ddtGroup->setTime(time);
    ddtGroup->setStepNumber(ii);
    ddtGroup->setTimeStep(dt);

    if (ii>0) {refRes = 1.0;}
    else  {refRes = 0.0;}

    testExpression.evaluateFunction(result);   
    ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));
    
    copyExpression.evaluateFunction(result);   
    ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));
    
    assignExpression.evaluateFunction(result);   
    ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserDerivativeTest, ddt1)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserDerivativeTest, ddt2)
{
  Teuchos::RCP<sdtExpressionGroup> ddtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = ddtGroup;

  Xyce::Util::newExpression testExpression(std::string("DDT (3.0*sin(time))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result = 0.0;
  std::complex<double> refRes = 0.0;
  double time=0.0, finalTime=0.01;
  int numSteps = 1001;
  double dt = finalTime/(numSteps-1);

  for (int ii=0;ii<numSteps;ii++)
  {
    ddtGroup->setTime(time);
    ddtGroup->setStepNumber(ii);
    ddtGroup->setTimeStep(dt);

    if (ii > 0) { refRes = 3.0*std::cos(time); }
    else { refRes = 0.0; }

    testExpression.evaluateFunction(result);   
    ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));
    
    copyExpression.evaluateFunction(result);   
    ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));
    
    assignExpression.evaluateFunction(result);   
    ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserDerivativeTest, ddt2)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserDerivativeTest, ddt3)
{
  Teuchos::RCP<sdtExpressionGroup> ddtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = ddtGroup;

  Xyce::Util::newExpression testExpression(std::string("DDT(v(a))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result = 0.0;
  std::complex<double> refRes = 0.0;
  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refDerivs;

  double time=0.0;
  double finalTime=1.0;

  int numSteps = 1001;
  double dt = finalTime/(numSteps-1);
  refDerivs.resize(1,0.5*dt);

  for (int ii=0;ii<numSteps;ii++)
  {
    std::complex<double> Aval=time;
    ddtGroup->setSoln(std::string("A"),Aval);
    ddtGroup->setTime(time);
    ddtGroup->setStepNumber(ii);
    ddtGroup->setTimeStep(dt); 

    testExpression.evaluate(result,derivs);   

    if (ii>0)
    {
      refDerivs[0] = 1.0/dt;
      refRes = 1.0;

      ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));
      ASSERT_FLOAT_EQ( std::real(derivs[0]), std::real(refDerivs[0]) );
    }

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserDerivativeTest, ddt3)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserDerivativeTest, ddt4)
{
  Teuchos::RCP<sdtExpressionGroup> ddtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = ddtGroup;

  Xyce::Util::newExpression testExpression(std::string("DDT(2.0*v(a))"), testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result = 0.0; 
  std::complex<double> refRes = 0.0;
  std::vector<std::complex<double> > derivs;
  std::vector<std::complex<double> > refDerivs;

  double time=0.0;
  double finalTime=1.0;

  int numSteps = 101;
  double dt = finalTime/(numSteps-1);
  refDerivs.resize(1,0.5*dt);


  for (int ii=0;ii<numSteps;ii++)
  {
    std::complex<double> Aval=2*time;
    ddtGroup->setSoln(std::string("A"),Aval);
    ddtGroup->setTime(time);
    ddtGroup->setStepNumber(ii);
    ddtGroup->setTimeStep(dt); 
    testExpression.evaluate(result,derivs);   

    if (ii>0)
    {
      refDerivs[0] = 2.0/dt;
      refRes = 4.0;
      ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));
      ASSERT_FLOAT_EQ( std::real(derivs[0]), std::real(refDerivs[0]) );
    }

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserDerivativeTest, ddt4)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserDerivativeTest, ddt5)
{
  Teuchos::RCP<sdtExpressionGroup> ddtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = ddtGroup;

  Xyce::Util::newExpression testExpression(std::string("F1(V(A))"), testGroup);
  testExpression.lexAndParseExpression();

  // .func F1(A) {ddt(A)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(A)"), rhs=std::string("ddt(a)");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result = 0.0; 
  std::complex<double> refRes = 0.0;
  double time=0.0, finalTime=1.0;
  int numSteps = 101;
  double dt = finalTime/(numSteps-1);

  for (int ii=0;ii<numSteps;ii++)
  {
    std::complex<double> Aval=time;
    ddtGroup->setSoln(std::string("A"),Aval);
    ddtGroup->setTime(time);
    ddtGroup->setStepNumber(ii);
    ddtGroup->setTimeStep(dt);

    if (ii>0) {refRes = 1.0;}
    else { refRes = 0.0; }

    testExpression.evaluateFunction(result);     ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));
    copyExpression.evaluateFunction(result);     ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));
    assignExpression.evaluateFunction(result);   ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserDerivativeTest, ddt5)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserDerivativeTest, ddt6)
{
  Teuchos::RCP<sdtExpressionGroup> ddtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = ddtGroup;

  Xyce::Util::newExpression testExpression(std::string("F1(V(A))+F1(V(A))"), testGroup);
  testExpression.lexAndParseExpression();

  // .func F1(A) {ddt(A)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(A)");
  std::string rhs=std::string("ddt(a)");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result = 0.0; 
  std::complex<double> refRes = 0.0;
  double time=0.0, finalTime=1.0;
  int numSteps = 101;
  double dt = finalTime/(numSteps-1);

  for (int ii=0;ii<numSteps;ii++)
  {
    std::complex<double> Aval=time;
    ddtGroup->setSoln(std::string("A"),Aval);
    ddtGroup->setTime(time);
    ddtGroup->setStepNumber(ii);
    ddtGroup->setTimeStep(dt);

    if (ii>0) {refRes = 2.0;}
    else { refRes = 0.0;}

    testExpression.evaluateFunction(result);     ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));
    copyExpression.evaluateFunction(result);     ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));
    assignExpression.evaluateFunction(result);   ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserDerivativeTest, ddt6)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserDerivativeTest, ddt7)
{
  Teuchos::RCP<sdtExpressionGroup> ddtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = ddtGroup;

  Xyce::Util::newExpression testExpression(std::string("F1(V(A))+F1(V(B))"), testGroup);
  testExpression.lexAndParseExpression();

  // .func F1(A) {ddt(A)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(A)");
  std::string rhs=std::string("ddt(a)");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  //testExpression.dumpParseTree(std::cout);

  std::complex<double> result = 0.0; 
  std::complex<double> refRes = 0.0;
  double time=0.0, finalTime=1.0;
  int numSteps = 1001;
  double dt = finalTime/(numSteps-1);

  for (int ii=0;ii<numSteps;ii++)
  {
    std::complex<double> Aval=time;
    std::complex<double> Bval=3.0*time;
    ddtGroup->setSoln(std::string("A"),Aval);
    ddtGroup->setSoln(std::string("B"),Bval);
    ddtGroup->setTime(time);
    ddtGroup->setStepNumber(ii);
    ddtGroup->setTimeStep(dt);
    if (ii>0) { refRes = 1.0 + 3.0; }
    else { refRes = 0.0; }

    testExpression.evaluateFunction(result);     ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));
    copyExpression.evaluateFunction(result);     ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));
    assignExpression.evaluateFunction(result);   ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));

    time += dt;
    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserDerivativeTest, ddt7)
}

//-------------------------------------------------------------------------------
// this test is similar to ddt7, except that the DDT operators 
// are behind 2 layers of funcs instead of 1.
TEST ( ComplexParserDerivativeTest, ddt8)
{
  Teuchos::RCP<sdtExpressionGroup> ddtGroup = Teuchos::rcp(new sdtExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = ddtGroup;

  Xyce::Util::newExpression testExpression(std::string("F1(V(A))+F1(V(B))"), testGroup);
  testExpression.lexAndParseExpression();

  // .func F1(A) {f2(A)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(A)");
  std::string rhs=std::string("f2(a)");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  // .func F2(A) {ddt(A)}
  std::string f2Name;
  Teuchos::RCP<Xyce::Util::newExpression> f2Expression;
  lhs=std::string("F2(A)");
  rhs=std::string("ddt(a)");
  createFunc(lhs,rhs,testGroup, f2Name,f2Expression);

  f1Expression->attachFunctionNode(f2Name, f2Expression);
  testExpression.attachFunctionNode(f1Name, f1Expression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  //testExpression.dumpParseTree(std::cout);

  std::complex<double> result = 0.0; 
  std::complex<double> refRes = 0.0;
  double time=0.0, finalTime=1.0;
  int numSteps = 1001;
  double dt = finalTime/(numSteps-1);
  for (int ii=0;ii<numSteps;ii++)
  {
    std::complex<double> Aval=time;
    std::complex<double> Bval=3.0*time;
    ddtGroup->setSoln(std::string("A"),Aval);
    ddtGroup->setSoln(std::string("B"),Bval);
    ddtGroup->setTime(time);
    ddtGroup->setStepNumber(ii);
    ddtGroup->setTimeStep(dt);

    if (ii>0) { refRes = 1.0 + 3.0; }
    else { refRes = 0.0; }

    testExpression.evaluateFunction(result);     ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));
    copyExpression.evaluateFunction(result);     ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));
    assignExpression.evaluateFunction(result);   ASSERT_FLOAT_EQ( std::real(result), std::real(refRes));

    time += dt;

    Xyce::Util::newExpression::clearProcessSuccessfulTimeStepMap();
    testExpression.processSuccessfulTimeStep();
    copyExpression.processSuccessfulTimeStep();
    assignExpression.processSuccessfulTimeStep();
  }

  OUTPUT_MACRO(ComplexParserDerivativeTest, ddt8)
}

//-------------------------------------------------------------------------------
// breakpoint function testing
//
// The two stp tests, below will cause the function computeBreakPoint to be 
// called in ast.h.  That function evaluates a Newton loop to obtain the next
// breakpoint.  However, for the STP operator, in these tests, this is a linear
// problem.  So, it won't use the bottom part of the function, where it continues
// for subsequent iterations.
//-------------------------------------------------------------------------------
TEST ( ComplexParserBreakpointTest, stp1)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;

  // this expression will use the .func stpTest.
  Xyce::Util::newExpression testExpression(std::string("stpTest(time)*0.5"), testGroup);
  testExpression.lexAndParseExpression();

  // .func stpTest(t) {t-0.5}
  std::string stpTestName;
  Teuchos::RCP<Xyce::Util::newExpression> stpTestExpression;
  std::string lhs=std::string("stpTest(t)");
  std::string rhs=std::string("stp(T-0.5)");
  createFunc(lhs,rhs,testGroup,stpTestName,stpTestExpression);

  testExpression.attachFunctionNode(stpTestName, stpTestExpression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  //testExpression.dumpParseTree(std::cout);

  timeDepGroup->setTime(0.4);
  std::complex<double> result = 0.0;

  {
  testExpression.evaluateFunction(result);
  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  testExpression.getBreakPoints(breakPointTimes);
  int numBp = breakPointTimes.size();
  std::vector<double> bpVals(1,0.0);
  if (numBp > 0) { bpVals[0] = breakPointTimes[0].value(); }
  EXPECT_EQ( numBp, 1 ); EXPECT_EQ( bpVals[0], 0.5 );
  }

  {
  copyExpression.evaluateFunction(result);
  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  copyExpression.getBreakPoints(breakPointTimes);
  int numBp = breakPointTimes.size();
  std::vector<double> bpVals(1,0.0);
  if (numBp > 0) { bpVals[0] = breakPointTimes[0].value(); }
  EXPECT_EQ( numBp, 1 ); EXPECT_EQ( bpVals[0], 0.5 );
  }

  {
  assignExpression.evaluateFunction(result);
  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  assignExpression.getBreakPoints(breakPointTimes);
  int numBp = breakPointTimes.size();
  std::vector<double> bpVals(1,0.0);
  if (numBp > 0) { bpVals[0] = breakPointTimes[0].value(); }
  EXPECT_EQ( numBp, 1 ); EXPECT_EQ( bpVals[0], 0.5 );
  }

  OUTPUT_MACRO(ComplexParserBreakpointTest, stp1)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserBreakpointTest, stp2)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;

  // this expression will use the .func f1.
  Xyce::Util::newExpression testExpression(std::string("f1(1.0)*stp(time-0.5)"), testGroup);
  testExpression.lexAndParseExpression();

  // .func f1(x) {5.0*x}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("f1(x)");
  std::string rhs=std::string("5.0*x");
  createFunc(lhs,rhs,testGroup,f1Name,f1Expression);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  timeDepGroup->setTime(0.4);
  std::complex<double> result;
  testExpression.evaluateFunction(result);

  //testExpression.dumpParseTree(std::cout);

  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  testExpression.getBreakPoints(breakPointTimes);
  int numBp = breakPointTimes.size();
  std::vector<double> bpVals(1,0.0);
  if (numBp > 0) { bpVals[0] = breakPointTimes[0].value(); }

  EXPECT_EQ( numBp, 1 );
  EXPECT_EQ( bpVals[0], 0.5 );

  OUTPUT_MACRO(ComplexParserBreakpointTest, stp2)
}

//-------------------------------------------------------------------------------
// limitOp breakpoint tests.  These are very similar to the STP breakpoint tests
// as they rely on the same computeBreakpoints function.
//-------------------------------------------------------------------------------
TEST ( ComplexParserBreakpointTest, limit1)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;

  // this expression will use the .func limitTest.
  Xyce::Util::newExpression testExpression(std::string("limitTest(time)*0.5"), testGroup);
  testExpression.lexAndParseExpression();

  // .func limitTest(t) {t-0.5}
  std::string limitTestName;
  Teuchos::RCP<Xyce::Util::newExpression> limitTestExpression;
  std::string lhs=std::string("limitTest(t)");
  std::string rhs=std::string("limit(T-0.5,0,1)");
  createFunc(lhs,rhs,testGroup,limitTestName,limitTestExpression);

  testExpression.attachFunctionNode(limitTestName, limitTestExpression);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  //testExpression.dumpParseTree(std::cout);

  timeDepGroup->setTime(0.4);
  std::complex<double> result = 0.0;

  {
  testExpression.evaluateFunction(result);
  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  testExpression.getBreakPoints(breakPointTimes);
  int numBp = breakPointTimes.size();
  std::vector<double> bpVals(numBp,0.0);
  for (int ii=0;ii<numBp;ii++) { bpVals[ii] = breakPointTimes[ii].value(); }
  EXPECT_EQ( numBp, 2 ); EXPECT_EQ( bpVals[0], 0.5 ); EXPECT_EQ( bpVals[1], 1.5 );
  }

  {
  copyExpression.evaluateFunction(result);
  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  copyExpression.getBreakPoints(breakPointTimes);
  int numBp = breakPointTimes.size();
  std::vector<double> bpVals(numBp,0.0);
  for (int ii=0;ii<numBp;ii++) { bpVals[ii] = breakPointTimes[ii].value(); }
  EXPECT_EQ( numBp, 2 ); EXPECT_EQ( bpVals[0], 0.5 ); EXPECT_EQ( bpVals[1], 1.5 );
  }

  {
  assignExpression.evaluateFunction(result);
  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  assignExpression.getBreakPoints(breakPointTimes);
  int numBp = breakPointTimes.size();
  std::vector<double> bpVals(numBp,0.0);
  for (int ii=0;ii<numBp;ii++) { bpVals[ii] = breakPointTimes[ii].value(); }
  EXPECT_EQ( numBp, 2 ); EXPECT_EQ( bpVals[0], 0.5 ); EXPECT_EQ( bpVals[1], 1.5 );
  }

  OUTPUT_MACRO(ComplexParserBreakpointTest, limit1)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserBreakpointTest, limit2)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;

  // this expression will use the .func f1.
  Xyce::Util::newExpression testExpression(std::string("f1(1.0)*limit(time-0.5,0,1)"), testGroup);
  testExpression.lexAndParseExpression();

  // .func f1(x) {5.0*x}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("f1(x)");
  std::string rhs=std::string("5.0*x");
  createFunc(lhs,rhs,testGroup,f1Name,f1Expression);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  timeDepGroup->setTime(0.4);
  std::complex<double> result = 0.0;
  testExpression.evaluateFunction(result);

  //testExpression.dumpParseTree(std::cout);

  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  testExpression.getBreakPoints(breakPointTimes);
  int numBp = breakPointTimes.size();
  std::vector<double> bpVals(numBp,0.0);
  for (int ii=0;ii<numBp;ii++) { bpVals[ii] = breakPointTimes[ii].value(); }
  EXPECT_EQ( numBp, 2);
  EXPECT_EQ( bpVals[0], 0.5 );
  EXPECT_EQ( bpVals[1], 1.5 );

  OUTPUT_MACRO(ComplexParserBreakpointTest, limit2)
}

//-------------------------------------------------------------------------------
// Breakpoint test inspired by ABM_BREAK/breaks.cir:
//
// V1   1  0  PULSE(0 3 0.01 1ms 1ms 0.05 0.15)
// B2   2  0  V = {Table(time, 0, 0, 0.3, 0, 0.301, 2, 0.302, 2, 0.6, 1, 1, 1)}
// B3   3  0  V = {v(2) + v(1) * if(abs(sin(5*PI*time)) > 0.9, (abs(sin(5*PI*time))-0.9)*10, 0)}
//
//
TEST ( ComplexParserBreakpointTest, abm_breaks1)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;

  // this expression is inspired by B3, above.
  // It forces the breakpointing function to call computeBreakPoint, 
  // and to do a multi-iteration Newton solve to obtain the breakpoint.  This
  // is necessary b/c "sin" is a nonlinear function.
  Xyce::Util::newExpression testExpression(std::string("{2.0 + 3.0 * if(abs(sin(5*PI*time)) > 0.9, (abs(sin(5*PI*time))-0.9)*10, 0)}"), testGroup);
  testExpression.lexAndParseExpression();

  // do some math.   sin(5*pi*time) = 0.9 ->  sin^(-1)(0.9) = 5*pi*time   ->     time = sin^(-1)(0.9)/(5*pi) = firstBP
  double firstBP = std::asin(0.9)/(5.0*M_PI);
  double seconBP = (M_PI - std::asin(0.9))/(5.0*M_PI);
  //std::cout << "firstBP = " << firstBP << std::endl;
  //std::cout << "seconBP = " << seconBP << std::endl;

  //testExpression.dumpParseTree(std::cout);

  std::complex<double> result;
  timeDepGroup->setTime(0.01);
  testExpression.evaluateFunction(result); 
  // This "evaluateFunction" call is to force the full AST setup, 
  // so that the arrays related to breakpoints are correct.  
  // This test does not check "result"

  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  testExpression.getBreakPoints(breakPointTimes);

  int numBp = breakPointTimes.size();
  if (numBp > 0) 
  { 
    double val = breakPointTimes[0].value(); 
    EXPECT_DOUBLE_EQ( val, firstBP );
  }

  timeDepGroup->setTime(seconBP-0.01);
  testExpression.evaluateFunction(result);
  breakPointTimes.clear();
  testExpression.getBreakPoints(breakPointTimes);
  if (numBp > 0)  
  { 
    double val = breakPointTimes[0].value(); 
    EXPECT_DOUBLE_EQ( val, seconBP );
  }

  OUTPUT_MACRO(ComplexParserBreakpointTest, abm_breaks1)
}

//-------------------------------------------------------------------------------
// made up test, for the equiv operator
TEST ( ComplexParserBreakpointTest, timeSquared1)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;

  // It forces the breakpointing function to call computeBreakPoint, 
  // and to do a multi-iteration Newton solve to obtain the breakpoint.  This
  // is necessary b/c "time*time" is a nonlinear function.
  Xyce::Util::newExpression testExpression(std::string("{2.0 + 3.0 * if((time*time == 4.0 ), 2.0, 1.0)}"), testGroup);
  testExpression.lexAndParseExpression();

  double firstBP = 2.0;

  //testExpression.dumpParseTree(std::cout);

  std::complex<double> result;
  timeDepGroup->setTime(0.01);
  testExpression.evaluateFunction(result); 
  // This "evaluateFunction" call is to force the full AST setup, 
  // so that the arrays related to breakpoints are correct.  
  // This test does not check "result"

  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  testExpression.getBreakPoints(breakPointTimes);

  int numBp = breakPointTimes.size();
  if (numBp == 1) 
  { 
    double val = breakPointTimes[0].value(); 
    EXPECT_DOUBLE_EQ( val, firstBP );
  }
  else
  {
    EXPECT_EQ( numBp, 1);
  }

  OUTPUT_MACRO(ComplexParserBreakpointTest, timeSquared1)
}

//-------------------------------------------------------------------------------
// made up test, for the equiv operator, which uses it in a .func
TEST ( ComplexParserBreakpointTest, timeSquared2)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = timeDepGroup;

  // this expression will use the .func f1.
  // It forces the breakpointing function to call computeBreakPoint, 
  // and to do a multi-iteration Newton solve to obtain the breakpoint.  This
  // is necessary b/c "time*time" is a nonlinear function.
  Xyce::Util::newExpression testExpression(std::string("{2.0 + 3.0 * f1(time)}"), testGroup);
  testExpression.lexAndParseExpression();

  // .func F1(x) {(if((x*x == 4.0 ), 2.0, 1.0)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(x)");
  std::string rhs=std::string("if((x*x == 4.0 ), 2.0, 1.0)");
  createFunc(lhs,rhs,testGroup,f1Name,f1Expression);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  double firstBP = 2.0;

  //testExpression.dumpParseTree(std::cout);

  std::complex<double> result;
  timeDepGroup->setTime(0.01);
  testExpression.evaluateFunction(result); 
  // This "evaluateFunction" call is to force the full AST setup, 
  // so that the arrays related to breakpoints are correct.  
  // This test does not check "result"

  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  testExpression.getBreakPoints(breakPointTimes);

  int numBp = breakPointTimes.size();
  if (numBp == 1) 
  { 
    double val = breakPointTimes[0].value(); 
    EXPECT_DOUBLE_EQ( val, firstBP );
  }
  else
  {
    EXPECT_EQ( numBp, 1);
  }

  OUTPUT_MACRO(ComplexParserBreakpointTest, timeSquared2)
}

//-------------------------------------------------------------------------------
// breakpoint test, for table source
TEST ( ComplexParserBreakpointTest, tableBreakPoint)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = timeDepGroup;
  Xyce::Util::newExpression tableExpression(std::string("Table({time} 0, 0, 0.3, 0, 0.301, 2, 0.302, 2, 0.6, 1, 1, 1)"), grp);
  tableExpression.lexAndParseExpression();

  std::complex<double> result=0.0;
  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  timeDepGroup->setTime(0.0);
  tableExpression.evaluateFunction(result); 
  tableExpression.getBreakPoints(breakPointTimes);

  // there are 6 entries in the table, but the tableOp breakpoint function only returns the first 5
  int size = breakPointTimes.size();
  EXPECT_EQ(size,5);

  if (size==5)
  {
    std::vector<double> refTimes = {0, 0.3, 0.301, 0.302, 0.6};
    for(int ii=0;ii<size;ii++)
    {
      EXPECT_DOUBLE_EQ( refTimes[ii], breakPointTimes[ii].value() );
    }
  }

  OUTPUT_MACRO2(ComplexParserBreakpointTest, tableBreakPoint, tableExpression) 
}

//-------------------------------------------------------------------------------
// breakpoint test, for table source
TEST ( ComplexParserBreakpointTest, tableBreakPoint2)
{
  Teuchos::RCP<timeDepExpressionGroup> timeDepGroup = Teuchos::rcp(new timeDepExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp = timeDepGroup;
  Xyce::Util::newExpression testExpression(std::string("f1(2.0)"), grp);
  testExpression.lexAndParseExpression();

  // .func F1(x) {(if((x*x == 4.0 ), 2.0, 1.0)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(x)");
  std::string rhs=std::string("X*Table({time} 0, 0, 0.3, 0, 0.301, 2, 0.302, 2, 0.6, 1, 1, 1)");
  createFunc(lhs,rhs,grp,f1Name,f1Expression);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  std::complex<double> result=0.0;
  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  timeDepGroup->setTime(0.0);
  testExpression.evaluateFunction(result); 
  testExpression.getBreakPoints(breakPointTimes);

  // there are 6 entries in the table, but the tableOp breakpoint function only returns the first 5
  int size = breakPointTimes.size();
  EXPECT_EQ(size,5);

  if (size==5)
  {
    std::vector<double> refTimes = {0, 0.3, 0.301, 0.302, 0.6};
    for(int ii=0;ii<size;ii++)
    {
      EXPECT_DOUBLE_EQ( refTimes[ii], breakPointTimes[ii].value() );
    }
  }

  OUTPUT_MACRO(ComplexParserBreakpointTest, tableBreakPoint2)
}


//-------------------------------------------------------------------------------
// Testing out error trapping.  This needs some work.

//-------------------------------------------------------------------------------
// For this test, see invalid math error message (See Message/Function/invalid_math_operator.cir, which has 'b =# n' )
//
//  As of this writing (6/27/2020) the new expression library doesn't believe this is an error.  It incorrectly passes.
TEST ( ComplexParserErrorTest, invalidMath)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  //Xyce::Util::newExpression testExpression(std::string("b # n"), testGroup); // 
  Xyce::Util::newExpression testExpression(std::string("b =# n"), testGroup); // this is what is in the regression test.  it now exits with error, but it was necessary to create a special rule in the parser to trap for it.
  //Xyce::Util::newExpression testExpression(std::string("b * "), testGroup); // this exits with an error
  //Xyce::Util::newExpression testExpression(std::string("sin()"), testGroup); // this exits with an error
  testExpression.lexAndParseExpression();
  //testExpression.dumpParseTree(std::cout);
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
//  ASSERT_EQ( (result-(1.0)), 0.0);
}

//-------------------------------------------------------------------------------
// For this test see bug 850 SON, which has the following bad_user_defined_func.cir:
// 
//  * User defined functions (udf) from the bug report for SON Bug 850.
//  * These functions were called UGAUSS in the bug report.
//  .FUNC udfA(mu,sigma) {exp(mu/sigma)}       ; User Defined Func. A
//  .func udfB(mu,sigma,t) {exp((t-mu)/sigma)} ; User Defined Func. B
//
//  VTEST    1    0     10
//  R1       1    2     50
//
//  * All of these B-source statements should produce a netlist parsing error.
//  * The inline comments give their behavior with Xyce 6.7.
//  
//  * udfA should have two arguments, but is given three arguments
//  B1 2 0     I = {(V(2)+udfA(10m,1u,23))/50} ; segfault with Func. A
//  B2 2 0     I = {udfA(10m,1u,23)/50}   ; runs fine with Func. A, but shouldn't
//  B3 2 0     I = {udfA(10m,1u,23)+V(2)} ; segfaults with Func. A.
//  B4 2 0     I = {udfA(10m,1u,23)}      ; correctly aborts with Func. A 
//  
//  * udfB should have three arguments, and is only given two arguments
//  B5 2 0     I = {(V(2)+udfB(10m,1u))/50}  ; segfault with Func. B
//  B6 2 0     I = {(udfB(10m,1u)+V(2))/50}  ; segfault with Func. B  
//  B7 2 0     I = {(udfB(20m,1u))/50}    ; runs fine with Func. B, but shouldn't
//
//  As of this writing (6/27/2020) the new expression library DOES believe these 
//  are errors.  However, as I haven't yet figured out how to properly do unit 
//  tests on error messages, I've hacked the code below to (rightfully) pass.
//-------------------------------------------------------------------------------
TEST ( ComplexParserErrorTest,  bad_user_defined_func )
{
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression will use the .func udfA.
  //Xyce::Util::newExpression testExpression(std::string("udfA(10m,1u,23)/50"), testGroup); // this, correctly fails
  Xyce::Util::newExpression testExpression(std::string("udfA(10p,1u)/50"), testGroup); // this passes, correctly.
  testExpression.lexAndParseExpression();

  // .func udfA(mu,sigma) {exp(mu/sigma)}
  std::string udfAName;
  Teuchos::RCP<Xyce::Util::newExpression> udfAExpression;
  std::string lhs=std::string("udfA(mu,sigma)");
  std::string rhs=std::string("exp(mu/sigma)");
  createFunc(lhs,rhs,testGroup, udfAName,udfAExpression);

  testExpression.attachFunctionNode(udfAName, udfAExpression);

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result;
  testExpression.evaluateFunction(result);   
  EXPECT_NEAR( std::real(result), std::real(std::exp( 10.0e-12/1.0e-6 )/50.0), 1e-12 );
  EXPECT_NEAR( std::imag(result), std::imag(std::exp( 10.0e-12/1.0e-6 )/50.0), 1e-12 );
  copyExpression.evaluateFunction(result);   
  EXPECT_NEAR( std::real(result), std::real(std::exp( 10.0e-12/1.0e-6 )/50.0), 1e-12 );
  EXPECT_NEAR( std::imag(result), std::imag(std::exp( 10.0e-12/1.0e-6 )/50.0), 1e-12 );
  assignExpression.evaluateFunction(result); 
  EXPECT_NEAR( std::real(result), std::real(std::exp( 10.0e-12/1.0e-6 )/50.0), 1e-12 );
  EXPECT_NEAR( std::imag(result), std::imag(std::exp( 10.0e-12/1.0e-6 )/50.0), 1e-12 );
  OUTPUT_MACRO(ComplexParserErrorTest,  bad_user_defined_func)
}

//-------------------------------------------------------------------------------
// This test is to make sure new expression library doesn't crash when you
// attach an invalid function.  The invalid function is the udfA function, which
// has a syntax error.
//
// This test does NOT compare against a test value.
//
// The test succeeds if the testExpression.evaluateFunction function can be
// called without causing a seg fault.  It should also test to ensure that the
// rhs expression triggers a syntax error, but I haven't figured out how to test
// that yet.
//
// Good news: no segfault.
// Bad news:  (as of 9/1/2020) no syntax error for extra paren at the end.
//-------------------------------------------------------------------------------
TEST ( ComplexParserErrorTest,  bad_user_defined_func2 )
{
  Teuchos::RCP<testExpressionGroupWithFuncSupport> funcGroup = Teuchos::rcp(new testExpressionGroupWithFuncSupport() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = funcGroup;

  // this expression will use the .func udfA.
  Xyce::Util::newExpression testExpression(std::string("udfA(10p,1u)/50"), testGroup);
  testExpression.lexAndParseExpression();

  // The correct function declaration would be: .func udfA(mu,sigma) {exp(mu/sigma)}
  // But I am using an incorrect one:           .func udfA(mu,sigma) {exp(mu/sigma))} (extra paren)
  // But I am using an incorrect one:           .func udfA(mu,sigma) {exp(mu/sigma} (missing paren)
  std::string udfAName;
  Teuchos::RCP<Xyce::Util::newExpression> udfAExpression;
  std::string lhs=std::string("udfA(mu,sigma)");
  //std::string rhs=std::string("exp(mu/sigma))");
  std::string rhs=std::string("exp(mu/sigma");
  createFunc(lhs,rhs,testGroup, udfAName,udfAExpression);

  testExpression.attachFunctionNode(udfAName, udfAExpression);

  //testExpression.dumpParseTree(std::cout);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  std::complex<double> result;
  testExpression.evaluateFunction(result);
  copyExpression.evaluateFunction(result);
  assignExpression.evaluateFunction(result);
  OUTPUT_MACRO(ComplexParserErrorTest,  bad_user_defined_func2)
}

//-------------------------------------------------------------------------------
// From: SENS/improperObjFormat.cir
//
//  objfunc={{I(VM),V(3)*V(3)}
//
//-------------------------------------------------------------------------------
TEST ( ComplexParserErrorTest,  bad_obj_func)
{
  Teuchos::RCP<solutionGroup> solGroup = Teuchos::rcp(new solutionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solGroup;

  //Xyce::Util::newExpression testExpression(std::string("POW(I(VM),V(3),V(3))"), testGroup);  // this should fail, and does
  //Xyce::Util::newExpression testExpression(std::string("(I(VM),V(3)*V(3))"), testGroup);  // this should fail in parsing, and does
  //Xyce::Util::newExpression testExpression(std::string("{I(VM),V(3)*V(3)}"), testGroup);  // this should fail in parsing, and does
  //Xyce::Util::newExpression testExpression(std::string("I(VM),V(3)*V(3)"), testGroup);  // this should fail in parsing, but does not
  Xyce::Util::newExpression testExpression(std::string("I(VM)*V(3)*V(3)"), testGroup);
  testExpression.lexAndParseExpression();

  std::complex<double> vmValue=2.0;
  std::complex<double> threeValue=3.0;

  solGroup->setSoln(std::string("VM"), vmValue);
  solGroup->setSoln(std::string("3"), threeValue);

  //Xyce::Util::newExpression copyExpression(testExpression);
  //Xyce::Util::newExpression assignExpression;
  //assignExpression = testExpression;

  //testExpression.dumpParseTree(std::cout);

  std::complex<double> result;
  testExpression.evaluateFunction(result);   ASSERT_EQ( result, 18.0 );
  //copyExpression.evaluateFunction(result);   ASSERT_EQ( result, 1.0 );
  //assignExpression.evaluateFunction(result); ASSERT_EQ( result, 1.0 );
  OUTPUT_MACRO(ComplexParserErrorTest,  bad_obj_func)
}

//-------------------------------------------------------------------------------
// ERK. 6/30/2020.
//
// This failed to parse in a recent circuit:  '1e-9*(w==2.0u)+1e9*(w!=2.0u)'
//
// The following tests are about this.  It appears that the "==" operator works
// but the "!=" operator does not, as of this writing.  The reason, as it turns
// out is because I added "!" to the allowed characters for parameter strings
// So, w!=2.0 can be interpretted as a parameter named "w!" equals 2.0.  This doesn't
// happen if there is a space between the w and the "!".  The fix appears to be
// to just exclude "!" from the TOK_ID in the lexer file.
//-------------------------------------------------------------------------------
TEST ( ComplexParserinlineComp, equiv)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Teuchos::RCP<Xyce::Util::newExpression> wExpression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2u"), baseGroup));
  wExpression->lexAndParseExpression();
  std::string wName="w";

  Xyce::Util::newExpression e11(std::string("1e-9*(w==2.0u)"), baseGroup);
  e11.lexAndParseExpression();
  e11.attachParameterNode(wName,wExpression);

  Xyce::Util::newExpression copy_e11(e11);
  Xyce::Util::newExpression assign_e11;
  assign_e11 = e11;

  std::complex<double> result=0.0;
  e11.evaluateFunction(result);        ASSERT_EQ( result, 1.0e-9);
  copy_e11.evaluateFunction(result);   ASSERT_EQ( result, 1.0e-9);
  assign_e11.evaluateFunction(result); ASSERT_EQ( result, 1.0e-9);

  OUTPUT_MACRO2(ComplexParserTernary_precedence, equiv, e11)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserinlineComp, notEquiv)
{
  Teuchos::RCP<ifStatementExpressionGroup> ifGroup = Teuchos::rcp(new ifStatementExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroup = ifGroup;

  Teuchos::RCP<Xyce::Util::newExpression> wExpression = Teuchos::rcp(new Xyce::Util::newExpression(std::string("2"), baseGroup));
  wExpression->lexAndParseExpression();
  std::string wName="w";

  Xyce::Util::newExpression e11(std::string("1e9*(w!=2.0u)"), baseGroup);
  e11.lexAndParseExpression();
  e11.attachParameterNode(wName,wExpression);

  Xyce::Util::newExpression copy_e11(e11);
  Xyce::Util::newExpression assign_e11;
  assign_e11 = e11;

  std::complex<double> result=0.0;
  e11.evaluateFunction(result);        ASSERT_EQ( result, 1.0e9);
  copy_e11.evaluateFunction(result);   ASSERT_EQ( result, 1.0e9);
  assign_e11.evaluateFunction(result); ASSERT_EQ( result, 1.0e9);

  OUTPUT_MACRO2(ComplexParserTernary_precedence, notEquiv, e11)
}


// This test emits the following error messages:
//
// Netlist error: Function or variable W(YACC!ACC1) is not defined
// Netlist error: Function or variable P(K1) is not defined
// Netlist error: Function or variable W(K2) is not defined
// Netlist error: Function or variable P(OLINE1) is not defined
// Netlist error: Function or variable P(UAND1) is not defined
// Netlist error: Function or variable P(YOR!OR1) is not defined
//
TEST ( ComplexParserErrorTest,  power_unsupported_devices)
{
  Teuchos::RCP<solutionGroup> solGroup = Teuchos::rcp(new solutionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solGroup;

  Xyce::Util::newExpression testExpression(std::string("{W(YACC!ACC1)}"), testGroup);
  //Xyce::Util::newExpression testExpression(std::string("{P(K1)}"), testGroup);
  testExpression.lexAndParseExpression();

  //Xyce::Util::newExpression copyExpression(testExpression);
  //Xyce::Util::newExpression assignExpression;
  //assignExpression = testExpression;

  //testExpression.dumpParseTree(std::cout);

  std::complex<double> result (0,0); // this number is the "can't find it" number
  std::complex<double> reference(-2e+09,0); // this number is the "can't find it" number
  testExpression.evaluateFunction(result);   ASSERT_EQ( result, reference );
  //copyExpression.evaluateFunction(result);   ASSERT_EQ( result, 1.0 );
  //assignExpression.evaluateFunction(result); ASSERT_EQ( result, 1.0 );
  OUTPUT_MACRO(ComplexParserErrorTest,  bad_obj_func)
}


//-------------------------------------------------------------------------------
//
// The point of this function is just to get thru the lex/parse phase w/o a
// syntax error.  So, I don't care about the result.
//
// The issue here (prompted by bug 28) is that the name "int" is used in that
// test as a node name.  And, that node name is initially processed as an
// expression.  The parser was interpretting "int" as the operator, when it
// needed to be perceiving it as an unresolved string.
//
// Note to self: I changed the lexer and parser so that the "int" operator
// token (TOK_INT) now includes the left paren.  That way, if the left paren
// is not present, it it tokenized as the generic TOK_WORD.  Given that the
// Xyce parser is used in this way, to process node names that ultimately are
// not really expressions, I should probably update a lot of the
// operators in this manner to force them to include the left paren in their
// token.
//
//-------------------------------------------------------------------------------
TEST ( Complex_Parsing_Syntax, bug28_1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("Int"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result;
  testExpression.evaluateFunction(result);
}

#if 0
// Similar to the above with "int", theoretically a user could attempt to
// use operator names as parameter names.  For example (below) sdt.
//
// Q:  Should I fix this for all operators that take an argument inside of parens?
TEST ( Complex_Parsing_Syntax, bug28_2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("SDT"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result;
  testExpression.evaluateFunction(result);
}
#endif

//-------------------------------------------------------------------------------
// tests for random operators
//-------------------------------------------------------------------------------
TEST ( ComplexParserRandom, agauss0)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("agauss(1.0,0.1,1.0)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  ASSERT_EQ( result, 1.0);

  OUTPUT_MACRO(ComplexParserRandom, agauss0)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserRandom, agauss1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("agauss(1.0,0.1,1.0)"), testGroup);
  testExpression.lexAndParseExpression();

  std::complex<double> result1(0.0);
  std::complex<double> result2(0.0);
  testExpression.evaluateFunction(result1);
  testExpression.evaluateFunction(result2);

  ASSERT_EQ( result1, result2); 

  OUTPUT_MACRO(ComplexParserRandom, agauss1)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserRandom, agauss2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("agauss(1.0,0.1)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  ASSERT_EQ( result, 1.0);

  OUTPUT_MACRO(ComplexParserRandom, agauss2)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserRandom, agauss3)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("agauss(1.0,0.1)"), testGroup);
  testExpression.lexAndParseExpression();

  std::complex<double> result1(0.0);
  std::complex<double> result2(0.0);
  testExpression.evaluateFunction(result1);
  testExpression.evaluateFunction(result2);

  ASSERT_EQ( result1, result2); 

  OUTPUT_MACRO(ComplexParserRandom, agauss3)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserRandom, agauss1_func)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("f1(1.0)"), testGroup);
  testExpression.lexAndParseExpression();

  // .func F1(A) {A*agauss(1.0,0.1,1.0)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(A)");
  std::string rhs=std::string("A*agauss(1.0,0.1,1.0)");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  std::complex<double> result1(0.0);
  std::complex<double> result2(0.0);
  testExpression.evaluateFunction(result1);
  testExpression.evaluateFunction(result2);

  ASSERT_EQ( result1, result2); 

  OUTPUT_MACRO(ComplexParserRandom, agauss1_func)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserRandom, agauss2_func)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("f1(1.0)"), testGroup);
  testExpression.lexAndParseExpression();

  // .func F1(A) {A*agauss(1.0,0.1,1.0)}
  std::string f1Name;
  Teuchos::RCP<Xyce::Util::newExpression> f1Expression;
  std::string lhs=std::string("F1(A)");
  std::string rhs=std::string("A*agauss(1.0,0.1)");
  createFunc(lhs,rhs,testGroup, f1Name,f1Expression);

  testExpression.attachFunctionNode(f1Name, f1Expression);

  std::complex<double> result1(0.0);
  std::complex<double> result2(0.0);
  testExpression.evaluateFunction(result1);
  testExpression.evaluateFunction(result2);

  ASSERT_EQ( result1, result2); 

  OUTPUT_MACRO(ComplexParserRandom, agauss2_func)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserRandom, gauss0)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("gauss(1.0,0.1,1.0)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  ASSERT_EQ( result, 1.0);

  OUTPUT_MACRO(ComplexParserRandom, gauss0)
}

TEST ( ComplexParserRandom, gauss1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("gauss(1.0,0.1)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  ASSERT_EQ( result, 1.0);

  OUTPUT_MACRO(ComplexParserRandom, gauss1)
}

//-------------------------------------------------------------------------------
TEST ( ComplexParserRandom, aunif0)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("aunif(1.0,0.1)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  ASSERT_EQ( result, 1.0);

  OUTPUT_MACRO(ComplexParserRandom, aunif0)
}

TEST ( ComplexParserRandom, unif0)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("unif(1.0,0.1)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  ASSERT_EQ( result, 1.0);

  OUTPUT_MACRO(ComplexParserRandom, unif0)
}

TEST ( Double_Parser_Random, rand0)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("rand()"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  ASSERT_EQ( result, 0.5);

  OUTPUT_MACRO(Double_Parser_Random, rand0)
}

TEST ( Double_Parser_Random, limit0)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("limit(1.0,0.5)"), testGroup);
  testExpression.lexAndParseExpression();
  std::complex<double> result(0.0);
  testExpression.evaluateFunction(result);
  ASSERT_EQ( result, 1.0);

  OUTPUT_MACRO(Double_Parser_Random, limit0)
}

// .print operators for S-params, Y-params, and Z-params
// from test Xyce_Regression/Netlists/Output/SPARAMS/sparams-ts2-dataFormat.cir
//.PRINT AC FILE=sparams-ts2-dataFormat.cir.FD.exp.prn
//+ {S(1,1)} {SR(1,1)} {SI(1,1)} {SM(1,1)} {SP(1,1)} {SDB(1,1)}
//+ {S(1,2)} {SR(1,2)} {SI(1,2)} {SM(1,2)} {SP(1,2)} {SDB(1,2)}
//+ {S(2,1)} {SR(2,1)} {SI(2,1)} {SM(2,1)} {SP(2,1)} {SDB(2,1)}
//+ {S(2,2)} {SR(2,2)} {SI(2,2)} {SM(2,2)} {SP(2,2)} {SDB(2,2)}
//+ {Y(1,1)} {YR(1,1)} {YI(1,1)} {YM(1,1)} {YP(1,1)} {YDB(1,1)}
//+ {Y(1,2)} {YR(1,2)} {YI(1,2)} {YM(1,2)} {YP(1,2)} {YDB(1,2)}
//+ {Y(2,1)} {YR(2,1)} {YI(2,1)} {YM(2,1)} {YP(2,1)} {YDB(2,1)}
//+ {Y(2,2)} {YR(2,2)} {YI(2,2)} {YM(2,2)} {YP(2,2)} {YDB(2,2)}
//+ {Z(1,1)} {ZR(1,1)} {ZI(1,1)} {ZM(1,1)} {ZP(1,1)} {ZDB(1,1)}
//+ {Z(1,2)} {ZR(1,2)} {ZI(1,2)} {ZM(1,2)} {ZP(1,2)} {ZDB(1,2)}
//+ {Z(2,1)} {ZR(2,1)} {ZI(2,1)} {ZM(2,1)} {ZP(2,1)} {ZDB(2,1)}
//+ {Z(2,2)} {ZR(2,2)} {ZI(2,2)} {ZM(2,2)} {ZP(2,2)} {ZDB(2,2)}


//SR(1,1)	-9.99200320e-01
//SI(1,1)	-6.36110782e-06
//SR(1,2)	-3.69390430e-10
//SI(1,2)	-1.27273037e-06
//SR(2,1)	-3.69390430e-10
//SI(2,1)	-1.27273037e-06
//SR(2,2)	5.99999973e-01
//SI(2,2)	-1.59093567e-04

//YR(1,1)	4.99968336e+01
//YI(1,1)	3.97862158e-01
//YR(1,2)	-3.08890787e-07
//YI(1,2)	3.97862762e-05
//YR(2,1)	-3.08890787e-07
//YI(2,1)	3.97862762e-05
//YR(2,2)	5.00000014e-03
//YI(2,2)	2.48583681e-06

//ZR(1,1)	2.00000000e-02
//ZI(1,1)	-1.59154943e-04
//ZR(1,2)	-1.09999962e-07
//ZI(1,2)	-1.59154878e-04
//ZR(2,1)	-1.09999962e-07
//ZI(2,1)	-1.59154878e-04
//ZR(2,2)	1.99999944e+02
//ZI(2,2)	-9.94334503e-02

TEST ( ComplexParserSparamTest, sparam1)
{
  Teuchos::RCP<rfParamGroup> rfGroup = Teuchos::rcp(new rfParamGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = rfGroup;

  Xyce::Util::newExpression testExpression_s(std::string("s(1,1)"), testGroup);
  Xyce::Util::newExpression testExpression_sr(std::string("sr(1,1)"), testGroup);
  Xyce::Util::newExpression testExpression_si(std::string("si(1,1)"), testGroup);
  Xyce::Util::newExpression testExpression_sm(std::string("sm(1,1)"), testGroup);
  Xyce::Util::newExpression testExpression_sdb(std::string("sdb(1,1)"), testGroup);

  testExpression_s.lexAndParseExpression();
  testExpression_sr.lexAndParseExpression();
  testExpression_si.lexAndParseExpression();
  testExpression_sm.lexAndParseExpression();
  testExpression_sdb.lexAndParseExpression();


  Xyce::Util::newExpression copyExpr_s(testExpression_s);Xyce::Util::newExpression assignExpr_s;assignExpr_s=testExpression_s; 
  Xyce::Util::newExpression copyExpr_sr(testExpression_sr);Xyce::Util::newExpression assignExpr_sr;assignExpr_sr=testExpression_sr; 
  Xyce::Util::newExpression copyExpr_si(testExpression_si);Xyce::Util::newExpression assignExpr_si;assignExpr_si=testExpression_si; 
  Xyce::Util::newExpression copyExpr_sm(testExpression_sm);Xyce::Util::newExpression assignExpr_sm;assignExpr_sm=testExpression_sm; 
  Xyce::Util::newExpression copyExpr_sdb(testExpression_sdb);Xyce::Util::newExpression assignExpr_sdb;assignExpr_sdb=testExpression_sdb; 
  
  std::vector< std::vector< std::complex<double> > > spars(2);
  spars[0] = 
  {
    std::complex<double>( -9.99200320e-01, -6.36110782e-06), // SR(1,1)	SI(1,1)	
    std::complex<double>( -3.69390430e-10, -1.27273037e-06) // SR(1,2)	SI(1,2)	
  };
  spars[1] = {
    std::complex<double>( -3.69390430e-10, -1.27273037e-06), // SR(2,1)	SI(2,1)	
    std::complex<double>( 5.99999973e-01, -1.59093567e-04) // SR(2,2)	SI(2,2)	
  };

  rfGroup->setSparam(spars);

  std::complex<double>  result=0.0;

  testExpression_s.evaluateFunction(result);   ASSERT_EQ( result, spars[0][0]);
  copyExpr_s.evaluateFunction(result);   ASSERT_EQ( result, spars[0][0]);
  assignExpr_s.evaluateFunction(result); ASSERT_EQ( result, spars[0][0]);


  testExpression_sr.evaluateFunction(result);   ASSERT_EQ( result, std::real(spars[0][0]));
  copyExpr_sr.evaluateFunction(result);   ASSERT_EQ( result, std::real(spars[0][0]));
  assignExpr_sr.evaluateFunction(result); ASSERT_EQ( result, std::real(spars[0][0]));

  testExpression_si.evaluateFunction(result);   ASSERT_EQ( result, std::imag(spars[0][0]));
  copyExpr_si.evaluateFunction(result);   ASSERT_EQ( result, std::imag(spars[0][0]));
  assignExpr_si.evaluateFunction(result); ASSERT_EQ( result, std::imag(spars[0][0]));

  testExpression_sm.evaluateFunction(result);   ASSERT_EQ( result, std::abs(spars[0][0]));
  copyExpr_sm.evaluateFunction(result);   ASSERT_EQ( result, std::abs(spars[0][0]));
  assignExpr_sm.evaluateFunction(result); ASSERT_EQ( result, std::abs(spars[0][0]));


  testExpression_sdb.evaluateFunction(result);   ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(spars[0][0])) );
  copyExpr_sdb.evaluateFunction(result);   ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(spars[0][0])) );
  assignExpr_sdb.evaluateFunction(result); ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(spars[0][0])));
  
  //OUTPUT_MACRO ( ComplexParserSparamTest, sparam1)
}


TEST ( ComplexParserSparamTest, sparam2)
{
  Teuchos::RCP<rfParamGroup> rfGroup = Teuchos::rcp(new rfParamGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = rfGroup;

  Xyce::Util::newExpression testExpression_s(std::string("S (2,1)"), testGroup);
  Xyce::Util::newExpression testExpression_sr(std::string("sr (2,1)"), testGroup);
  Xyce::Util::newExpression testExpression_si(std::string("sI (2,1)"), testGroup);
  Xyce::Util::newExpression testExpression_sm(std::string("SM (2,1)"), testGroup);
  Xyce::Util::newExpression testExpression_sdb(std::string("sdb (2,1)"), testGroup);

  testExpression_s.lexAndParseExpression();
  testExpression_sr.lexAndParseExpression();
  testExpression_si.lexAndParseExpression();
  testExpression_sm.lexAndParseExpression();
  testExpression_sdb.lexAndParseExpression();


  Xyce::Util::newExpression copyExpr_s(testExpression_s);Xyce::Util::newExpression assignExpr_s;assignExpr_s=testExpression_s; 
  Xyce::Util::newExpression copyExpr_sr(testExpression_sr);Xyce::Util::newExpression assignExpr_sr;assignExpr_sr=testExpression_sr; 
  Xyce::Util::newExpression copyExpr_si(testExpression_si);Xyce::Util::newExpression assignExpr_si;assignExpr_si=testExpression_si; 
  Xyce::Util::newExpression copyExpr_sm(testExpression_sm);Xyce::Util::newExpression assignExpr_sm;assignExpr_sm=testExpression_sm; 
  Xyce::Util::newExpression copyExpr_sdb(testExpression_sdb);Xyce::Util::newExpression assignExpr_sdb;assignExpr_sdb=testExpression_sdb; 
  
  std::vector< std::vector< std::complex<double> > > spars(2);
  spars[0] = 
  {
    std::complex<double>( -9.99200320e-01, -6.36110782e-06), // SR(1,1)	SI(1,1)	
    std::complex<double>( -3.69390430e-10, -1.27273037e-06) // SR(1,2)	SI(1,2)	
  };
  spars[1] = {
    std::complex<double>( -3.69390430e-10, -1.27273037e-06), // SR(2,1)	SI(2,1)	
    std::complex<double>( 5.99999973e-01, -1.59093567e-04) // SR(2,2)	SI(2,2)	
  };

  rfGroup->setSparam(spars);

  std::complex<double>  result=0.0;

  testExpression_s.evaluateFunction(result);   ASSERT_EQ( result, spars[1][0]);
  copyExpr_s.evaluateFunction(result);   ASSERT_EQ( result, spars[1][0]);
  assignExpr_s.evaluateFunction(result); ASSERT_EQ( result, spars[1][0]);


  testExpression_sr.evaluateFunction(result);   ASSERT_EQ( result, std::real(spars[1][0]));
  copyExpr_sr.evaluateFunction(result);   ASSERT_EQ( result, std::real(spars[1][0]));
  assignExpr_sr.evaluateFunction(result); ASSERT_EQ( result, std::real(spars[1][0]));

  testExpression_si.evaluateFunction(result);   ASSERT_EQ( result, std::imag(spars[1][0]));
  copyExpr_si.evaluateFunction(result);   ASSERT_EQ( result, std::imag(spars[1][0]));
  assignExpr_si.evaluateFunction(result); ASSERT_EQ( result, std::imag(spars[1][0]));

  testExpression_sm.evaluateFunction(result);   ASSERT_EQ( result, std::abs(spars[1][0]));
  copyExpr_sm.evaluateFunction(result);   ASSERT_EQ( result, std::abs(spars[1][0]));
  assignExpr_sm.evaluateFunction(result); ASSERT_EQ( result, std::abs(spars[1][0]));


  testExpression_sdb.evaluateFunction(result);   ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(spars[1][0])) );
  copyExpr_sdb.evaluateFunction(result);   ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(spars[1][0])) );
  assignExpr_sdb.evaluateFunction(result); ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(spars[1][0])));
  
  //OUTPUT_MACRO ( ComplexParserSparamTest, sparam2)
}

TEST ( ComplexParserYparamTest, yparam1)
{
  Teuchos::RCP<rfParamGroup> rfGroup = Teuchos::rcp(new rfParamGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = rfGroup;

  Xyce::Util::newExpression testExpression_s(std::string("y(1,1)"), testGroup);
  Xyce::Util::newExpression testExpression_sr(std::string("yr(1,1)"), testGroup);
  Xyce::Util::newExpression testExpression_si(std::string("yi(1,1)"), testGroup);
  Xyce::Util::newExpression testExpression_sm(std::string("ym(1,1)"), testGroup);
  Xyce::Util::newExpression testExpression_sdb(std::string("ydb(1,1)"), testGroup);

  testExpression_s.lexAndParseExpression();
  testExpression_sr.lexAndParseExpression();
  testExpression_si.lexAndParseExpression();
  testExpression_sm.lexAndParseExpression();
  testExpression_sdb.lexAndParseExpression();


  Xyce::Util::newExpression copyExpr_s(testExpression_s);Xyce::Util::newExpression assignExpr_s;assignExpr_s=testExpression_s; 
  Xyce::Util::newExpression copyExpr_sr(testExpression_sr);Xyce::Util::newExpression assignExpr_sr;assignExpr_sr=testExpression_sr; 
  Xyce::Util::newExpression copyExpr_si(testExpression_si);Xyce::Util::newExpression assignExpr_si;assignExpr_si=testExpression_si; 
  Xyce::Util::newExpression copyExpr_sm(testExpression_sm);Xyce::Util::newExpression assignExpr_sm;assignExpr_sm=testExpression_sm; 
  Xyce::Util::newExpression copyExpr_sdb(testExpression_sdb);Xyce::Util::newExpression assignExpr_sdb;assignExpr_sdb=testExpression_sdb; 
  
  std::vector< std::vector< std::complex<double> > > ypars(2);
  ypars[0] = 
  {
    std::complex<double>(4.99968336e+01, 3.97862158e-01), //YR(1,1)	YI(1,1)	
    std::complex<double>(-3.08890787e-07, 3.97862762e-05) // YR(1,2)	YI(1,2)	
  };
  ypars[1] =   
  {
    std::complex<double>(-3.08890787e-07, 3.97862762e-05), //YR(2,1)	YI(2,1)	
    std::complex<double>(5.00000014e-03, 2.48583681e-06) //YR(2,2)	YI(2,2)	
  };  
  rfGroup->setYparam(ypars);

  std::complex<double>  result=0.0;

  testExpression_s.evaluateFunction(result);   ASSERT_EQ( result, ypars[0][0]);
  copyExpr_s.evaluateFunction(result);   ASSERT_EQ( result, ypars[0][0]);
  assignExpr_s.evaluateFunction(result); ASSERT_EQ( result, ypars[0][0]);


  testExpression_sr.evaluateFunction(result);   ASSERT_EQ( result, std::real(ypars[0][0]));
  copyExpr_sr.evaluateFunction(result);   ASSERT_EQ( result, std::real(ypars[0][0]));
  assignExpr_sr.evaluateFunction(result); ASSERT_EQ( result, std::real(ypars[0][0]));

  testExpression_si.evaluateFunction(result);   ASSERT_EQ( result, std::imag(ypars[0][0]));
  copyExpr_si.evaluateFunction(result);   ASSERT_EQ( result, std::imag(ypars[0][0]));
  assignExpr_si.evaluateFunction(result); ASSERT_EQ( result, std::imag(ypars[0][0]));

  testExpression_sm.evaluateFunction(result);   ASSERT_EQ( result, std::abs(ypars[0][0]));
  copyExpr_sm.evaluateFunction(result);   ASSERT_EQ( result, std::abs(ypars[0][0]));
  assignExpr_sm.evaluateFunction(result); ASSERT_EQ( result, std::abs(ypars[0][0]));


  testExpression_sdb.evaluateFunction(result);   ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(ypars[0][0])) );
  copyExpr_sdb.evaluateFunction(result);   ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(ypars[0][0])) );
  assignExpr_sdb.evaluateFunction(result); ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(ypars[0][0])));
  
  //OUTPUT_MACRO ( ComplexParserYparamTest, yparam1)
}

TEST ( ComplexParserYparamTest, yparam2)
{
  Teuchos::RCP<rfParamGroup> rfGroup = Teuchos::rcp(new rfParamGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = rfGroup;

  Xyce::Util::newExpression testExpression_s(std::string("y (1,2)"), testGroup);
  Xyce::Util::newExpression testExpression_sr(std::string("yr (1,2)"), testGroup);
  Xyce::Util::newExpression testExpression_si(std::string("Yi (1,2)"), testGroup);
  Xyce::Util::newExpression testExpression_sm(std::string("ym (1,2)"), testGroup);
  Xyce::Util::newExpression testExpression_sdb(std::string("YDB (1,2)"), testGroup);

  testExpression_s.lexAndParseExpression();
  testExpression_sr.lexAndParseExpression();
  testExpression_si.lexAndParseExpression();
  testExpression_sm.lexAndParseExpression();
  testExpression_sdb.lexAndParseExpression();


  Xyce::Util::newExpression copyExpr_s(testExpression_s);Xyce::Util::newExpression assignExpr_s;assignExpr_s=testExpression_s; 
  Xyce::Util::newExpression copyExpr_sr(testExpression_sr);Xyce::Util::newExpression assignExpr_sr;assignExpr_sr=testExpression_sr; 
  Xyce::Util::newExpression copyExpr_si(testExpression_si);Xyce::Util::newExpression assignExpr_si;assignExpr_si=testExpression_si; 
  Xyce::Util::newExpression copyExpr_sm(testExpression_sm);Xyce::Util::newExpression assignExpr_sm;assignExpr_sm=testExpression_sm; 
  Xyce::Util::newExpression copyExpr_sdb(testExpression_sdb);Xyce::Util::newExpression assignExpr_sdb;assignExpr_sdb=testExpression_sdb; 
  
  std::vector< std::vector< std::complex<double> > > ypars(2);
  ypars[0] = 
  {
    std::complex<double>(4.99968336e+01, 3.97862158e-01), //YR(1,1)	YI(1,1)	
    std::complex<double>(-3.08890787e-07, 3.97862762e-05) // YR(1,2)	YI(1,2)	
  };
  ypars[1] =   
  {
    std::complex<double>(-3.08890787e-07, 3.97862762e-05), //YR(2,1)	YI(2,1)	
    std::complex<double>(5.00000014e-03, 2.48583681e-06) //YR(2,2)	YI(2,2)	
  };  
  rfGroup->setYparam(ypars);

  std::complex<double>  result=0.0;

  testExpression_s.evaluateFunction(result);   ASSERT_EQ( result, ypars[0][1]);
  copyExpr_s.evaluateFunction(result);   ASSERT_EQ( result, ypars[0][1]);
  assignExpr_s.evaluateFunction(result); ASSERT_EQ( result, ypars[0][1]);


  testExpression_sr.evaluateFunction(result);   ASSERT_EQ( result, std::real(ypars[0][1]));
  copyExpr_sr.evaluateFunction(result);   ASSERT_EQ( result, std::real(ypars[0][1]));
  assignExpr_sr.evaluateFunction(result); ASSERT_EQ( result, std::real(ypars[0][1]));

  testExpression_si.evaluateFunction(result);   ASSERT_EQ( result, std::imag(ypars[0][1]));
  copyExpr_si.evaluateFunction(result);   ASSERT_EQ( result, std::imag(ypars[0][1]));
  assignExpr_si.evaluateFunction(result); ASSERT_EQ( result, std::imag(ypars[0][1]));

  testExpression_sm.evaluateFunction(result);   ASSERT_EQ( result, std::abs(ypars[0][1]));
  copyExpr_sm.evaluateFunction(result);   ASSERT_EQ( result, std::abs(ypars[0][1]));
  assignExpr_sm.evaluateFunction(result); ASSERT_EQ( result, std::abs(ypars[0][1]));


  testExpression_sdb.evaluateFunction(result);   ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(ypars[0][1])) );
  copyExpr_sdb.evaluateFunction(result);   ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(ypars[0][1])) );
  assignExpr_sdb.evaluateFunction(result); ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(ypars[0][1])));
  
  //OUTPUT_MACRO ( ComplexParserYparamTest, yparam2)
}

TEST ( ComplexParserZparamTest, zparam1)
{
  Teuchos::RCP<rfParamGroup> rfGroup = Teuchos::rcp(new rfParamGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = rfGroup;

  Xyce::Util::newExpression testExpression_s(std::string("z(1,1)"), testGroup);
  Xyce::Util::newExpression testExpression_sr(std::string("zr(1,1)"), testGroup);
  Xyce::Util::newExpression testExpression_si(std::string("zi(1,1)"), testGroup);
  Xyce::Util::newExpression testExpression_sm(std::string("zm(1,1)"), testGroup);
  Xyce::Util::newExpression testExpression_sdb(std::string("zdb(1,1)"), testGroup);

  testExpression_s.lexAndParseExpression();
  testExpression_sr.lexAndParseExpression();
  testExpression_si.lexAndParseExpression();
  testExpression_sm.lexAndParseExpression();
  testExpression_sdb.lexAndParseExpression();


  Xyce::Util::newExpression copyExpr_s(testExpression_s);Xyce::Util::newExpression assignExpr_s;assignExpr_s=testExpression_s; 
  Xyce::Util::newExpression copyExpr_sr(testExpression_sr);Xyce::Util::newExpression assignExpr_sr;assignExpr_sr=testExpression_sr; 
  Xyce::Util::newExpression copyExpr_si(testExpression_si);Xyce::Util::newExpression assignExpr_si;assignExpr_si=testExpression_si; 
  Xyce::Util::newExpression copyExpr_sm(testExpression_sm);Xyce::Util::newExpression assignExpr_sm;assignExpr_sm=testExpression_sm; 
  Xyce::Util::newExpression copyExpr_sdb(testExpression_sdb);Xyce::Util::newExpression assignExpr_sdb;assignExpr_sdb=testExpression_sdb; 
  
  std::vector< std::vector< std::complex<double> > > zpars(2);
  zpars[0] = 
  {
    std::complex<double>(4.99968336e+01, 3.97862158e-01), //ZR(1,1)	ZI(1,1)	
    std::complex<double>(-3.08890787e-07, 3.97862762e-05) // ZR(1,2)	ZI(1,2)	
  };
  zpars[1] =   
  {
    std::complex<double>(-3.08890787e-07, 3.97862762e-05), //ZR(2,1)	ZI(2,1)	
    std::complex<double>(5.00000014e-03, 2.48583681e-06) //ZR(2,2)	ZI(2,2)	
  }; 
  rfGroup->setZparam(zpars);

  std::complex<double>  result=0.0;

  testExpression_s.evaluateFunction(result);   ASSERT_EQ( result, zpars[0][0]);
  copyExpr_s.evaluateFunction(result);   ASSERT_EQ( result, zpars[0][0]);
  assignExpr_s.evaluateFunction(result); ASSERT_EQ( result, zpars[0][0]);


  testExpression_sr.evaluateFunction(result);   ASSERT_EQ( result, std::real(zpars[0][0]));
  copyExpr_sr.evaluateFunction(result);   ASSERT_EQ( result, std::real(zpars[0][0]));
  assignExpr_sr.evaluateFunction(result); ASSERT_EQ( result, std::real(zpars[0][0]));

  testExpression_si.evaluateFunction(result);   ASSERT_EQ( result, std::imag(zpars[0][0]));
  copyExpr_si.evaluateFunction(result);   ASSERT_EQ( result, std::imag(zpars[0][0]));
  assignExpr_si.evaluateFunction(result); ASSERT_EQ( result, std::imag(zpars[0][0]));

  testExpression_sm.evaluateFunction(result);   ASSERT_EQ( result, std::abs(zpars[0][0]));
  copyExpr_sm.evaluateFunction(result);   ASSERT_EQ( result, std::abs(zpars[0][0]));
  assignExpr_sm.evaluateFunction(result); ASSERT_EQ( result, std::abs(zpars[0][0]));


  testExpression_sdb.evaluateFunction(result);   ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(zpars[0][0])) );
  copyExpr_sdb.evaluateFunction(result);   ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(zpars[0][0])) );
  assignExpr_sdb.evaluateFunction(result); ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(zpars[0][0])));
  
  //OUTPUT_MACRO ( ComplexParserZparamTest, zparam1)
}

TEST ( ComplexParserZparamTest, zparam2)
{
  Teuchos::RCP<rfParamGroup> rfGroup = Teuchos::rcp(new rfParamGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = rfGroup;

  Xyce::Util::newExpression testExpression_s(std::string("z (2,2)"), testGroup);
  Xyce::Util::newExpression testExpression_sr(std::string("ZR (2,2)"), testGroup);
  Xyce::Util::newExpression testExpression_si(std::string("zi (2,2)"), testGroup);
  Xyce::Util::newExpression testExpression_sm(std::string("ZM (2,2)"), testGroup);
  Xyce::Util::newExpression testExpression_sdb(std::string("zdb (2,2)"), testGroup);

  testExpression_s.lexAndParseExpression();
  testExpression_sr.lexAndParseExpression();
  testExpression_si.lexAndParseExpression();
  testExpression_sm.lexAndParseExpression();
  testExpression_sdb.lexAndParseExpression();


  Xyce::Util::newExpression copyExpr_s(testExpression_s);Xyce::Util::newExpression assignExpr_s;assignExpr_s=testExpression_s; 
  Xyce::Util::newExpression copyExpr_sr(testExpression_sr);Xyce::Util::newExpression assignExpr_sr;assignExpr_sr=testExpression_sr; 
  Xyce::Util::newExpression copyExpr_si(testExpression_si);Xyce::Util::newExpression assignExpr_si;assignExpr_si=testExpression_si; 
  Xyce::Util::newExpression copyExpr_sm(testExpression_sm);Xyce::Util::newExpression assignExpr_sm;assignExpr_sm=testExpression_sm; 
  Xyce::Util::newExpression copyExpr_sdb(testExpression_sdb);Xyce::Util::newExpression assignExpr_sdb;assignExpr_sdb=testExpression_sdb; 
  
  std::vector< std::vector< std::complex<double> > > zpars(2);
  zpars[0] = 
  {
    std::complex<double>(4.99968336e+01, 3.97862158e-01), //ZR(1,1)	ZI(1,1)	
    std::complex<double>(-3.08890787e-07, 3.97862762e-05) // ZR(1,2)	ZI(1,2)	
  };
  zpars[1] =   
  {
    std::complex<double>(-3.08890787e-07, 3.97862762e-05), //ZR(2,1)	ZI(2,1)	
    std::complex<double>(5.00000014e-03, 2.48583681e-06) //ZR(2,2)	ZI(2,2)	
  }; 
  rfGroup->setZparam(zpars);

  std::complex<double>  result=0.0;

  testExpression_s.evaluateFunction(result);   ASSERT_EQ( result, zpars[1][1]);
  copyExpr_s.evaluateFunction(result);   ASSERT_EQ( result, zpars[1][1]);
  assignExpr_s.evaluateFunction(result); ASSERT_EQ( result, zpars[1][1]);


  testExpression_sr.evaluateFunction(result);   ASSERT_EQ( result, std::real(zpars[1][1]));
  copyExpr_sr.evaluateFunction(result);   ASSERT_EQ( result, std::real(zpars[1][1]));
  assignExpr_sr.evaluateFunction(result); ASSERT_EQ( result, std::real(zpars[1][1]));

  testExpression_si.evaluateFunction(result);   ASSERT_EQ( result, std::imag(zpars[1][1]));
  copyExpr_si.evaluateFunction(result);   ASSERT_EQ( result, std::imag(zpars[1][1]));
  assignExpr_si.evaluateFunction(result); ASSERT_EQ( result, std::imag(zpars[1][1]));

  testExpression_sm.evaluateFunction(result);   ASSERT_EQ( result, std::abs(zpars[1][1]));
  copyExpr_sm.evaluateFunction(result);   ASSERT_EQ( result, std::abs(zpars[1][1]));
  assignExpr_sm.evaluateFunction(result); ASSERT_EQ( result, std::abs(zpars[1][1]));


  testExpression_sdb.evaluateFunction(result);   ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(zpars[1][1])) );
  copyExpr_sdb.evaluateFunction(result);   ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(zpars[1][1])) );
  assignExpr_sdb.evaluateFunction(result); ASSERT_FLOAT_EQ( std::real(result), 20.0*std::log10( std::abs(zpars[1][1])));
  
  //OUTPUT_MACRO ( ComplexParserZparamTest, zparam2)
}


// testing the "getIsComplex" function
TEST ( ComplexParserTestCmplxBoolean, isComplex1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("1.0+2.0J"), testGroup);
  testExpression.lexAndParseExpression();

  bool isComplex = testExpression.getIsComplex ();
  ASSERT_TRUE (isComplex);

  Xyce::Util::newExpression copyExpression(testExpression);
  isComplex = copyExpression.getIsComplex ();
  ASSERT_TRUE (isComplex);

  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;
  isComplex = assignExpression.getIsComplex ();
  ASSERT_TRUE (isComplex);
}

// testing the "getIsComplex" function
TEST ( ComplexParserTestCmplxBoolean, isComplex2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("Re(1.0+2.0J)"), testGroup);
  testExpression.lexAndParseExpression();

  bool isComplex = testExpression.getIsComplex ();
  ASSERT_FALSE (isComplex);

  Xyce::Util::newExpression copyExpression(testExpression);
  isComplex = copyExpression.getIsComplex ();
  ASSERT_FALSE (isComplex);

  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;
  isComplex = assignExpression.getIsComplex ();
  ASSERT_FALSE (isComplex);
}

// testing the "getIsComplex" function
TEST ( ComplexParserTestCmplxBoolean, isComplex3)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("Img(1.0+2.0J)"), testGroup);
  testExpression.lexAndParseExpression();

  bool isComplex = testExpression.getIsComplex ();
  ASSERT_FALSE (isComplex);

  Xyce::Util::newExpression copyExpression(testExpression);
  isComplex = copyExpression.getIsComplex ();
  ASSERT_FALSE (isComplex);

  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;
  isComplex = assignExpression.getIsComplex ();
  ASSERT_FALSE (isComplex);
}

// testing the "getIsComplex" function
TEST ( ComplexParserTestCmplxBoolean, isComplex4)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("Ph(1.0+2.0J)"), testGroup);
  testExpression.lexAndParseExpression();

  bool isComplex = testExpression.getIsComplex ();
  ASSERT_FALSE (isComplex);

  Xyce::Util::newExpression copyExpression(testExpression);
  isComplex = copyExpression.getIsComplex ();
  ASSERT_FALSE (isComplex);

  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;
  isComplex = assignExpression.getIsComplex ();
  ASSERT_FALSE (isComplex);
}

// testing the "getIsComplex" function
TEST ( ComplexParserTestCmplxBoolean, isComplex5)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("2.0"), testGroup);
  testExpression.lexAndParseExpression();

  bool isComplex = testExpression.getIsComplex ();
  ASSERT_FALSE (isComplex);

  Xyce::Util::newExpression copyExpression(testExpression);
  isComplex = copyExpression.getIsComplex ();
  ASSERT_FALSE (isComplex);

  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;
  isComplex = assignExpression.getIsComplex ();
  ASSERT_FALSE (isComplex);
}


// testing the "getIsComplex" function
// this expression is complex because while Ph(1.0+2.0J) is real, par1 is complex.
TEST ( ComplexParserTestCmplxBoolean, isComplex6)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("Ph(1.0+2.0J)*par1"), testGroup);
  testExpression.lexAndParseExpression();

  // setup the parameter
  Teuchos::RCP<Xyce::Util::newExpression> par1Expression = Teuchos::rcp(new Xyce::Util::newExpression (std::string("2.0-3.0J"), testGroup));
  par1Expression->lexAndParseExpression();
  std::string par1Name = "par1";
  testExpression.attachParameterNode(par1Name,par1Expression);

  bool isComplex = testExpression.getIsComplex ();
  ASSERT_TRUE (isComplex);

  Xyce::Util::newExpression copyExpression(testExpression);
  isComplex = copyExpression.getIsComplex ();
  ASSERT_TRUE (isComplex);

  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;
  isComplex = assignExpression.getIsComplex ();
  ASSERT_TRUE (isComplex);
}



// testing the "getVoltNameVec" function
TEST ( ComplexParserTest_getVoltNames, voltNames1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("V(A) + V(B) + V(A)"), testGroup);
  testExpression.lexAndParseExpression();

  int correctNumberNodes=2;

  {
  const std::vector<std::string> & voltNames = testExpression.getVoltNameVec();
  int numVoltageNodes= voltNames.size();
  ASSERT_EQ (numVoltageNodes, correctNumberNodes);
  }

  {
  Xyce::Util::newExpression copyExpression(testExpression);
  const std::vector<std::string> & voltNames = copyExpression.getVoltNameVec();
  int numVoltageNodes= voltNames.size();
  ASSERT_EQ (numVoltageNodes, correctNumberNodes);
  }

  {
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;
  const std::vector<std::string> & voltNames = assignExpression.getVoltNameVec();
  int numVoltageNodes= voltNames.size();
  ASSERT_EQ (numVoltageNodes, correctNumberNodes);
  }
}

// testing the "getVoltNameVec" function
TEST ( ComplexParserTest_getVoltNames, voltNames2)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("V(A) + V(0) + V(A)"), testGroup);
  testExpression.lexAndParseExpression();

  int correctNumberNodes=1;

  {
  const std::vector<std::string> & voltNames = testExpression.getVoltNameVec();
  int numVoltageNodes= voltNames.size();
  ASSERT_EQ (numVoltageNodes, correctNumberNodes);
  }

  {
  Xyce::Util::newExpression copyExpression(testExpression);
  const std::vector<std::string> & voltNames = copyExpression.getVoltNameVec();
  int numVoltageNodes= voltNames.size();
  ASSERT_EQ (numVoltageNodes, correctNumberNodes);
  }

  {
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;
  const std::vector<std::string> & voltNames = assignExpression.getVoltNameVec();
  int numVoltageNodes= voltNames.size();
  ASSERT_EQ (numVoltageNodes, correctNumberNodes);
  }
}

// disabling this test until it works
TEST ( ComplexParserTest_getVoltNames, voltNames3)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );
  Xyce::Util::newExpression testExpression(std::string("V(A) + V(Y) + V(Z)"), testGroup);
  testExpression.lexAndParseExpression();

  // the replaceName function is called when replacing a voltage node name with an alias 
  // (due to a subcircuit call)
  // The point of this test is to make sure that when nodes are replaced with redundant names, 
  // that the internal data structures of expression will properly update.
  std::string old_name1 = "Y";
  std::string old_name2 = "Z";
  std::string new_name1 = "A";
  std::string new_name2 = "A";
  testExpression.replaceName(old_name1, new_name1);
  testExpression.replaceName(old_name2, new_name2);

  int correctNumberNodes=1;
  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  {
  testExpression.setupVariousAstArrays();
  const std::vector<std::string> & voltNames = testExpression.getVoltNameVec();
  int numVoltageNodes= voltNames.size();
  ASSERT_EQ (numVoltageNodes, correctNumberNodes);
  //for (int ii=0;ii<numVoltageNodes; ii++) { std::cout << "voltNames3.  voltNames["<<ii<<"] = " << voltNames[ii] <<std::endl; }
  }

  {
  copyExpression.setupVariousAstArrays();
  const std::vector<std::string> & voltNames = copyExpression.getVoltNameVec();
  int numVoltageNodes= voltNames.size();
  ASSERT_EQ (numVoltageNodes, correctNumberNodes);
  }

  {
  assignExpression.setupVariousAstArrays();
  const std::vector<std::string> & voltNames = assignExpression.getVoltNameVec();
  int numVoltageNodes= voltNames.size();
  ASSERT_EQ (numVoltageNodes, correctNumberNodes);
  }
}

TEST ( ComplexParserTest_getVoltNames, voltNames4)
{
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;

  Xyce::Util::newExpression testExpression(std::string("V(A) + V(Y) + V(Z)"), testGroup);
  testExpression.lexAndParseExpression();

  // the replaceName function is called when replacing a voltage node name with an alias 
  // (due to a subcircuit call)
  // The point of this test is to make sure that when nodes are replaced with redundant names, 
  // that the internal data structures of expression will properly update.
  std::string old_name1 = "Y";
  std::string old_name2 = "Z";
  std::string new_name1 = "A";
  std::string new_name2 = "A";
  testExpression.replaceName(old_name1, new_name1);
  testExpression.replaceName(old_name2, new_name2);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  // After replacement, there should be a single derivative, value = 3;
  // This is b/c after replacement, the expression is equivalent to {3*V(A)}
  int correctNumberDerivs=1;

  std::complex<double> result(0.0,0.0), Aval(1.0,0.0);
  std::complex<double> refRes = 3.0*Aval;

  solnGroup->setSoln(std::string("A"),Aval);

  std::vector<std::complex<double> > refDer;
  refDer.push_back(3.0);

  std::vector<std::complex<double> > derivs;

  testExpression.evaluate(result,derivs);   
  ASSERT_EQ (correctNumberDerivs,derivs.size());
  EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);

  copyExpression.evaluate(result,derivs);   
  ASSERT_EQ (correctNumberDerivs,derivs.size());
  EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);

  assignExpression.evaluate(result,derivs); 
  ASSERT_EQ (correctNumberDerivs,derivs.size());
  EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
}

TEST ( ComplexParserTest_getVoltNames, voltNames5)
{
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;

  Xyce::Util::newExpression testExpression(std::string("8.0*I(vb) + V(A) + V(Y) + V(Z)"), testGroup);
  testExpression.lexAndParseExpression();

  // the replaceName function is called when replacing a voltage node name with an alias 
  // (due to a subcircuit call)
  // The point of this test is to make sure that when nodes are replaced with redundant names, 
  // that the internal data structures of expression will properly update.
  std::string old_name1 = "Y";
  std::string old_name2 = "Z";
  std::string new_name1 = "A";
  std::string new_name2 = "A";
  testExpression.replaceName(old_name1, new_name1);
  testExpression.replaceName(old_name2, new_name2);

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  // After replacement, there should be two derivatives, value = 3 and 8
  // This is b/c after replacement, the expression is equivalent to {8*I(vb) + 3*V(A)}
  //
  // Derivative array from expression library are in order of voltages first, currents second.
  int correctNumberDerivs=2;

  std::complex<double> result(0.0,0.0), Aval(1.0,0.0);
  std::complex<double> VBval=std::complex<double>(8.0,0.0);
  solnGroup->setSoln(std::string("vb"),VBval);
  std::complex<double> refRes = 3.0*Aval + 8.0*VBval;

  solnGroup->setSoln(std::string("A"),Aval);

  std::vector<std::complex<double> > refDer;
  refDer.push_back(3.0);
  refDer.push_back(8.0);

  std::vector<std::complex<double> > derivs;

  testExpression.evaluate(result,derivs);   
  ASSERT_EQ (correctNumberDerivs,derivs.size());
  EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);

  copyExpression.evaluate(result,derivs);   
  ASSERT_EQ (correctNumberDerivs,derivs.size());
  EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);

  assignExpression.evaluate(result,derivs); 
  ASSERT_EQ (correctNumberDerivs,derivs.size());
  EXPECT_EQ( result, refRes); EXPECT_EQ( derivs, refDer);
}

TEST ( ComplexParserTest_getVoltNames, voltNames6)
{
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;

  Xyce::Util::newExpression testExpression(std::string("V(A) + V(Y) + V(Z)"), testGroup);
  testExpression.lexAndParseExpression();

  // the replaceName function is called when replacing a voltage node name with an alias 
  // (due to a subcircuit call)
  // The point of this test is to make sure that when nodes are replaced with redundant names, 
  // that the internal data structures of expression will properly update.
  std::string old_name1 = "Y";
  std::string old_name2 = "Z";
  std::string new_name1 = "A";
  std::string new_name2 = "A";
  testExpression.replaceName(old_name1, new_name1);
  testExpression.replaceName(old_name2, new_name2);

  std::string gnd_name = "0";
  testExpression.replaceName(new_name1, gnd_name); // now everything is ground

  Xyce::Util::newExpression copyExpression(testExpression);
  Xyce::Util::newExpression assignExpression;
  assignExpression = testExpression;

  // After replacement, there should be a zero derivatives.
  // This is b/c after replacement, the expression is equivalent to {3*V(0)}
  int correctNumberDerivs=0;

  std::complex<double> result(0.0,0.0); 
  std::complex<double> Aval(1.0,0.0);
  std::complex<double> GNDval(0.0,0.0);
  std::complex<double> refRes = 0.0;

  solnGroup->setSoln(std::string("A"),Aval);
  solnGroup->setSoln(std::string("0"),GNDval);

  std::vector<std::complex<double> > derivs;

  testExpression.evaluate(result,derivs);   
  ASSERT_EQ (correctNumberDerivs,derivs.size());
  EXPECT_EQ( result, refRes); 

  copyExpression.evaluate(result,derivs);   
  ASSERT_EQ (correctNumberDerivs,derivs.size());
  EXPECT_EQ( result, refRes); 

  assignExpression.evaluate(result,derivs); 
  ASSERT_EQ (correctNumberDerivs,derivs.size());
  EXPECT_EQ( result, refRes); 
}


// testing the bracket operator in device names.  Test1 passes if it dosn't crash during parsing.
TEST ( ComplexParserDevNameWeirdChars, test1)
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  testGroup = Teuchos::rcp(new testExpressionGroup() );

  // original:
  std::string testExpr = "I(V_VDD)*V(_VDD)+I(V_VSS)*V(_VSS)+I(V_Q[3])*V(_Q[3])+I(V_Q[2])*V(_Q[2])+I(V_Q[1])*V(_Q[1])+I(V_Q[0])*V(_Q[0])+I(V_A[6])*V(_A[6])+I(V_A[5])*V(_A[5])+I(V_A[4])*V(_A[4])+I(V_A[3])*V(_A[3])+I(V_A[2])*V(_A[2])+I(V_A[1])*V(_A[1])+I(V_A[0])*V(_A[0])+I(V_E)*V(_E)+I(V_D[3])*V(_D[3])+I(V_D[2])*V(_D[2])+I(V_D[1])*V(_D[1])+I(V_D[0])*V(_D[0])+I(V_W[3])*V(_W[3])+I(V_W[2])*V(_W[2])+I(V_W[1])*V(_W[1])+I(V_W[0])*V(_W[0])";

  Xyce::Util::newExpression testExpression(testExpr, testGroup);

  testExpression.lexAndParseExpression();

}

// another test of bracket operator. This one looks at the answer.
TEST ( ComplexParserDevNameWeirdChars, test2)
{
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;

  // original:
  std::string testExpr = "17.2*I(V_Q[3])+8.5";
  Xyce::Util::newExpression testExpression(testExpr, testGroup);
  testExpression.lexAndParseExpression();

  Xyce::Util::newExpression copyExpression(testExpression); 
  Xyce::Util::newExpression assignExpression; 
  assignExpression = testExpression; 

  std::complex<double> result=0.0, V1val=3.0;
  std::complex<double> refRes = 17.2*V1val+8.5;
  solnGroup->setSoln(std::string("V_Q[3]"),V1val);
  testExpression.evaluateFunction(result);   EXPECT_EQ( std::real(result), std::real(refRes));
  copyExpression.evaluateFunction(result);   EXPECT_EQ( std::real(result), std::real(refRes));
  assignExpression.evaluateFunction(result); EXPECT_EQ( std::real(result), std::real(refRes));
}

TEST ( ComplexParserDevNameWeirdChars, test3)
{
  Teuchos::RCP<currSolnExpressionGroup> solnGroup = Teuchos::rcp(new currSolnExpressionGroup() );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> testGroup = solnGroup;


  std::vector<std::string> validDevs= 
  {
    std::string("1`"),
    std::string("1~"),
    std::string("1!"),
    std::string("1@"),
    std::string("1$"),
    std::string("1%"),
    std::string("1^"),
    std::string("1&"),
    std::string("1*"),
    std::string("1-"),
    std::string("1_"),
    std::string("1+"),
    std::string("1["),
    std::string("1]"),
    std::string("1|"),
    std::string("1\\"),
    std::string("1<"),
    std::string("1>"),
    std::string("1."),
    std::string("1?"),
    std::string("1/"),
    std::string("1#"),
    std::string("`"),
    std::string("~"),
    std::string("!"),
    std::string("@"),
    std::string("$"),
    std::string("%"),
    std::string("^"),
    std::string("&"),
    std::string("*"),
    std::string("-"),
    std::string("_"),
    std::string("+"),
    std::string("["),
    std::string("]"),
    std::string("|"),
    std::string("\\"),
    std::string("<"),
    std::string(">"),
    std::string("."),
    std::string("?"),
    std::string("/")
  };

  for (int ii=0;ii< validDevs.size();ii++)
  {
    std::string devName = "V" + validDevs[ii];

    std::string expressionString = "I(" + devName + ")";

    Xyce::Util::newExpression testExpression(expressionString, testGroup);
    testExpression.lexAndParseExpression();

    Xyce::Util::newExpression copyExpression(testExpression); 
    Xyce::Util::newExpression assignExpression; 
    assignExpression = testExpression; 

    std::complex<double> result=0.0, Aval=3.0;
    std::complex<double> refRes = Aval;
    solnGroup->setSoln(devName,Aval);

    testExpression.evaluateFunction(result);   EXPECT_EQ( std::real(result), std::real(refRes));
    copyExpression.evaluateFunction(result);   EXPECT_EQ( std::real(result), std::real(refRes));
    assignExpression.evaluateFunction(result); EXPECT_EQ( std::real(result), std::real(refRes));
    OUTPUT_MACRO(Double_Parser_VoltSolnTest, weirdChar)
  }
}

int main (int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}


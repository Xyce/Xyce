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

#include <Teuchos_RCP.hpp>

#include "ast.h"
#include <N_UTL_BreakPoint.h>
#include <N_UTL_ExtendedString.h>

using namespace Teuchos;

typedef std::complex<double> cmplx;

//ASSERT_TRUE(1 == 1);

//-------------------------------------------------------------------------------
// test values of binary operators 
//
#define AST_BINARY_OP_TEST_MACRO(NAME,OP,CPPOP, VAL1, VAL2) \
TEST ( NAME, OP ) \
{ \
  RCP<astNode<double> > arg1 = rcp(new numval<double> (VAL1)); \
  RCP<astNode<double> > arg2 = rcp(new numval<double> (VAL2)); \
  RCP<OP<double> > OP_1 = rcp(new OP<double>(arg1,arg2)); \
  EXPECT_DOUBLE_EQ(OP_1->val(),(VAL1 CPPOP VAL2)); \
}

#define AST_BINARY_OP_TEST_MACRO_CMPLX(NAME,OP,CPPOP, VAL1, VAL2) \
TEST ( NAME, OP ) \
{ \
  RCP<astNode<cmplx> > arg1 = rcp(new numval<cmplx> (VAL1)); \
  RCP<astNode<cmplx> > arg2 = rcp(new numval<cmplx> (VAL2)); \
  RCP<OP<cmplx> > OP_1 = rcp(new OP<cmplx>(arg1,arg2)); \
  EXPECT_DOUBLE_EQ(std::real(OP_1->val()),std::real(VAL1 CPPOP VAL2)); \
  EXPECT_DOUBLE_EQ(std::imag(OP_1->val()),std::imag(VAL1 CPPOP VAL2)); \
}

AST_BINARY_OP_TEST_MACRO(Double_Binary_Ast_Op,binaryAddOp,+,1.0,2.0) 
AST_BINARY_OP_TEST_MACRO(Double_Binary_Ast_Op,binaryMinusOp,-,1.0,2.0) 
AST_BINARY_OP_TEST_MACRO(Double_Binary_Ast_Op,binaryMulOp,*,1.0,2.0) 
AST_BINARY_OP_TEST_MACRO(Double_Binary_Ast_Op,binaryDivOp,/,1.0,2.0) 

AST_BINARY_OP_TEST_MACRO_CMPLX(Complex_Binary_Ast_Op,binaryAddOp,+,cmplx(1.0,0.5),cmplx(2.0,1.0)) 
AST_BINARY_OP_TEST_MACRO_CMPLX(Complex_Binary_Ast_Op,binaryMinusOp,-,cmplx(1.0,0.5),cmplx(2.0,1.0)) 
AST_BINARY_OP_TEST_MACRO_CMPLX(Complex_Binary_Ast_Op,binaryMulOp,*,cmplx(1.0,0.5),cmplx(2.0,1.0)) 
AST_BINARY_OP_TEST_MACRO_CMPLX(Complex_Binary_Ast_Op,binaryDivOp,/,cmplx(1.0,0.5),cmplx(2.0,1.0)) 

//-------------------------------------------------------------------------------
// test values of unary std functions 
//
#define AST_OP_TEST_MACRO(TYPE,NAME,OP,CPPFUNC, VAL) \
TEST ( NAME, OP )  \
{  \
  RCP<astNode<TYPE> > arg1 = rcp(new numval<TYPE> (VAL)); \
  RCP<OP<TYPE> > OP_1 = rcp(new OP<TYPE>(arg1)); \
  EXPECT_DOUBLE_EQ(OP_1->val(), CPPFUNC(VAL));  \
}

#define AST_OP_TEST_CMPLX_MACRO(TYPE,NAME,OP,CPPFUNC, VAL) \
TEST ( NAME, OP )  \
{  \
  RCP<astNode<TYPE> > arg1 = rcp(new numval<TYPE> (VAL)); \
  RCP<OP<TYPE> > OP_1 = rcp(new OP<TYPE>(arg1)); \
  EXPECT_DOUBLE_EQ(std::real(OP_1->val()), std::real(CPPFUNC(VAL)));  \
  EXPECT_DOUBLE_EQ(std::imag(OP_1->val()), std::imag(CPPFUNC(VAL)));  \
}

AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_Ops, sqrtOp,  std::sqrt,  4.0)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_Ops, expOp, std::exp, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_Ops, absOp, std::abs, -0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_Ops, sinOp, std::sin, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_Ops, cosOp, std::cos, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_Ops, acosOp, std::acos, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_Ops, acoshOp, std::acosh, 1.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_Ops, asinOp, std::asin, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_Ops, asinhOp, std::asinh, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_Ops, atanOp, std::atan, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_Ops, atanhOp, std::atanh, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_Ops, coshOp, std::cosh, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_Ops, logOp, std::log, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_Ops, log10Op, std::log10, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_Ops, sinhOp, std::sinh, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_Ops, tanOp, std::tan, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_Ops, tanhOp, std::tanh, 0.5)

AST_OP_TEST_CMPLX_MACRO(cmplx, Complex_UnaryFunc_Ast_Ops, sqrtOp,  std::sqrt, cmplx(4.0, 0.2))
AST_OP_TEST_CMPLX_MACRO(cmplx, Complex_UnaryFunc_Ast_Ops, expOp, std::exp, cmplx(.5, 0.2))
AST_OP_TEST_CMPLX_MACRO(cmplx, Complex_UnaryFunc_Ast_Ops, absOp, std::abs, cmplx(-.5, 0.2))
AST_OP_TEST_CMPLX_MACRO(cmplx, Complex_UnaryFunc_Ast_Ops, sinOp, std::sin, cmplx(.5, 0.2))
AST_OP_TEST_CMPLX_MACRO(cmplx, Complex_UnaryFunc_Ast_Ops, cosOp, std::cos, cmplx(.5, 0.2))
AST_OP_TEST_CMPLX_MACRO(cmplx, Complex_UnaryFunc_Ast_Ops, acosOp, std::acos, cmplx(.5, 0.2))
AST_OP_TEST_CMPLX_MACRO(cmplx, Complex_UnaryFunc_Ast_Ops, acoshOp, std::acosh, cmplx(1.5, 0.2))
AST_OP_TEST_CMPLX_MACRO(cmplx, Complex_UnaryFunc_Ast_Ops, asinOp, std::asin, cmplx(.5, 0.2))
AST_OP_TEST_CMPLX_MACRO(cmplx, Complex_UnaryFunc_Ast_Ops, asinhOp, std::asinh, cmplx(.5, 0.2))
AST_OP_TEST_CMPLX_MACRO(cmplx, Complex_UnaryFunc_Ast_Ops, atanOp, std::atan, cmplx(.5, 0.2))
AST_OP_TEST_CMPLX_MACRO(cmplx, Complex_UnaryFunc_Ast_Ops, atanhOp, std::atanh, cmplx(.5, 0.2))
AST_OP_TEST_CMPLX_MACRO(cmplx, Complex_UnaryFunc_Ast_Ops, coshOp, std::cosh, cmplx(.5, 0.2))
AST_OP_TEST_CMPLX_MACRO(cmplx, Complex_UnaryFunc_Ast_Ops, logOp, std::log, cmplx(.5, 0.2))
AST_OP_TEST_CMPLX_MACRO(cmplx, Complex_UnaryFunc_Ast_Ops, log10Op, std::log10, cmplx(.5, 0.2))
AST_OP_TEST_CMPLX_MACRO(cmplx, Complex_UnaryFunc_Ast_Ops, sinhOp, std::sinh, cmplx(.5, 0.2))
AST_OP_TEST_CMPLX_MACRO(cmplx, Complex_UnaryFunc_Ast_Ops, tanOp, std::tan, cmplx(.5, 0.2))
AST_OP_TEST_CMPLX_MACRO(cmplx, Complex_UnaryFunc_Ast_Ops, tanhOp, std::tanh, cmplx(.5, 0.2))

//-------------------------------------------------------------------------------
// constants
TEST ( Ast_Const_Test, Pi )
{  
  piConstOp<double> testPi;
  EXPECT_DOUBLE_EQ(testPi.val(), M_PI);  
}

//-------------------------------------------------------------------------------
// other functions
TEST ( Double_Ast_Func_Test, powOp )
{  
  RCP<astNode<double> > arg1 = rcp(new numval<double> (2.0));
  RCP<astNode<double> > arg2 = rcp(new numval<double> (3.0));
  RCP<astNode<double> > testPow = rcp(new powOp<double>(arg1,arg2)); 
  EXPECT_DOUBLE_EQ(testPow->val(), std::pow(2.0,3.0));  
}

TEST ( Double_Ast_Func_Test, pwrsOp )
{  
  RCP<astNode<double> > arg1 = rcp(new numval<double> (2.0));
  RCP<astNode<double> > arg2 = rcp(new numval<double> (3.0));
  RCP<astNode<double> > testPow = rcp(new pwrsOp<double>(arg1,arg2)); 
  EXPECT_DOUBLE_EQ(testPow->val(), std::pow(2.0,3.0));  
}

TEST ( Double_Ast_Func_Test, pwrsOp2 )
{  
  RCP<astNode<double> > arg1 = rcp(new numval<double> (0.0));
  RCP<astNode<double> > arg2 = rcp(new numval<double> (3.0));
  RCP<astNode<double> > testPow = rcp(new pwrsOp<double>(arg1,arg2)); 
  EXPECT_DOUBLE_EQ(testPow->val(), 0.0);  
}

TEST ( Double_Ast_Func_Test, pwrsOp3 )
{  
  RCP<astNode<double> > arg1 = rcp(new numval<double> (-2.0));
  RCP<astNode<double> > arg2 = rcp(new numval<double> (3.0));
  RCP<astNode<double> > testPow = rcp(new pwrsOp<double>(arg1,arg2)); 
  EXPECT_DOUBLE_EQ(testPow->val(), -std::pow(2.0,3.0) );
}

TEST ( Double_Ast_Func_Test, phaseOp )
{  
  double a1(-1.0);
  RCP<astNode<double> > arg1 = rcp(new numval<double> (a1));
  RCP<astNode<double> > testPhase = rcp(new phaseOp<double>(arg1)); 
  EXPECT_DOUBLE_EQ(testPhase->val(), std::arg(a1)*(180.0/M_PI));   // default behavior uses Degrrees, not radians
}

TEST ( Complex_Ast_Func_Test, powOp )
{  
  std::complex<double> a1(2.0,3.0);
  std::complex<double> a2(3.0,4.0);

  RCP<astNode<cmplx> > arg1 = rcp(new numval<cmplx> (a1));
  RCP<astNode<cmplx> > arg2 = rcp(new numval<cmplx> (a2));
  RCP<astNode<cmplx> > testPow = rcp(new powOp<cmplx>(arg1,arg2)); 

  EXPECT_DOUBLE_EQ(std::real(testPow->val()), std::real(std::pow(a1,a2)));
  EXPECT_DOUBLE_EQ(std::imag(testPow->val()), std::imag(std::pow(a1,a2)));
}

TEST ( Complex_Ast_Func_Test, pwrsOp )
{  
  std::complex<double> a1(2.0,3.0);
  std::complex<double> a2(3.0,4.0);

  RCP<astNode<cmplx> > arg1 = rcp(new numval<cmplx> (a1));
  RCP<astNode<cmplx> > arg2 = rcp(new numval<cmplx> (a2));
  RCP<astNode<cmplx> > testPow = rcp(new pwrsOp<cmplx>(arg1,arg2)); 

  EXPECT_DOUBLE_EQ(std::real(testPow->val()), std::real(std::pow(a1,a2)));
  EXPECT_DOUBLE_EQ(std::imag(testPow->val()), std::imag(std::pow(a1,a2)));
}

TEST ( Complex_Ast_Func_Test, pwrsOp2 )
{  
  std::complex<double> a1(0.0,0.0);
  std::complex<double> a2(3.0,4.0);

  RCP<astNode<cmplx> > arg1 = rcp(new numval<cmplx> (a1));
  RCP<astNode<cmplx> > arg2 = rcp(new numval<cmplx> (a2));
  RCP<astNode<cmplx> > testPow = rcp(new pwrsOp<cmplx>(arg1,arg2)); 

  EXPECT_DOUBLE_EQ(std::real(testPow->val()), 0.0);
  EXPECT_DOUBLE_EQ(std::imag(testPow->val()), 0.0);
}

TEST ( Complex_Ast_Func_Test, pwrsOp3 )
{
  std::complex<double> a1(-2.0,0.0);
  std::complex<double> a2(3.0,0.0);

  RCP<astNode<cmplx> > arg1 = rcp(new numval<cmplx> (a1));
  RCP<astNode<cmplx> > arg2 = rcp(new numval<cmplx> (a2));
  RCP<astNode<cmplx> > testPow = rcp(new pwrsOp<cmplx>(arg1,arg2)); 

  EXPECT_DOUBLE_EQ(std::real(testPow->val()),std::real( -std::pow( std::complex<double>(2.0,0.0),std::complex<double>(3.0,0.0) )));
  EXPECT_DOUBLE_EQ(std::imag(testPow->val()),std::imag( -std::pow( std::complex<double>(2.0,0.0),std::complex<double>(3.0,0.0) )));
}

TEST ( Complex_Ast_Func_Test, phaseOp )
{
  std::complex<double> a1(-1.0,0.0);
  RCP<astNode<cmplx> > arg1 = rcp(new numval<cmplx> (a1));
  RCP<astNode<cmplx> > testPhase = rcp(new phaseOp<cmplx>(arg1)); 
  EXPECT_DOUBLE_EQ(std::real(testPhase->val()), std::real(std::arg(a1)*(180.0/M_PI)));   // default behavior uses Degrrees, not radians
  EXPECT_DOUBLE_EQ(std::imag(testPhase->val()), std::imag(std::arg(a1)*(180.0/M_PI)));   // default behavior uses Degrrees, not radians
}

TEST ( Complex_Ast_Func_Test, realOp )
{  
  std::complex<double> a1(-1.0,0.3);
  RCP<astNode<cmplx> > arg1 = rcp(new numval<cmplx> (a1));
  RCP<astNode<cmplx> > testReal = rcp(new realOp<cmplx>(arg1)); 
  EXPECT_DOUBLE_EQ(std::real(testReal->val()), std::real(a1));  
}

TEST ( Complex_Ast_Func_Test, imagOp )
{  
  std::complex<double> a1(-1.0,0.3);
  RCP<astNode<cmplx> > arg1 = rcp(new numval<cmplx> (a1));
  RCP<astNode<cmplx> > testImag = rcp(new imagOp<cmplx>(arg1)); 
  EXPECT_DOUBLE_EQ(std::real(testImag->val()), std::imag(a1));   
}

TEST ( Double_Ast_Func_Test, unaryMinusOp)
{  
  double a1(1.0);
  RCP<astNode<double> > arg1 = rcp(new numval<double> (a1));
  RCP<astNode<double> > negarg1 = rcp(new unaryMinusOp<double>(arg1)); 

  EXPECT_DOUBLE_EQ(negarg1->val(), -1.0);
}

TEST ( Complex_Ast_Func_Test, unaryMinusOp)
{  
  std::complex<double> a1(1.0,1.0);
  RCP<astNode<std::complex<double> > > arg1 = rcp(new numval<std::complex<double> > (a1));
  RCP<astNode<std::complex<double> > > negarg1 = rcp(new unaryMinusOp<std::complex<double> >(arg1)); 

  EXPECT_DOUBLE_EQ(std::real(negarg1->val()), -1.0);
  EXPECT_DOUBLE_EQ(std::imag(negarg1->val()), -1.0);
}



//-------------------------------------------------------------------------------
// spice time-dependent source functions

TEST ( Double_Ast_Spice_Src_Test, spicePulseOp )
{  
  // these numbers are from the VPULSE regression test:
  //0V 1V 0S 10US 10US 0.1US 20.1US
  double v1=0.0, v2=1.0, td=0.0, tr=10.0e-6, tf=10.0e-6, pw=0.1e-6, per=20.1e-6, time = 0.0;

  RCP<astNode<double> > v1_op = rcp(new numval<double> (v1));
  RCP<astNode<double> > v2_op = rcp(new numval<double> (v2));
  RCP<astNode<double> > td_op = rcp(new numval<double> (td));
  RCP<astNode<double> > tr_op = rcp(new numval<double> (tr));
  RCP<astNode<double> > tf_op = rcp(new numval<double> (tf));
  RCP<astNode<double> > pw_op = rcp(new numval<double> (pw));
  RCP<astNode<double> > per_op = rcp(new numval<double> (per));
  RCP<astNode<double> > time_op = rcp(new specialsOp<double> (std::string("time")));
  RCP<specialsOp<double> > time_op2 = rcp_static_cast<specialsOp<double> > (time_op);

  std::vector< RCP<astNode<double> > > sourceArgs = { v1_op, v2_op, td_op, tr_op, tf_op, pw_op, per_op };

  // test initial value(v1)
  time_op2->setValue(time);
  spicePulseOp<double> pulse( sourceArgs, time_op);
  EXPECT_DOUBLE_EQ(pulse.val(), v1);

  // test post-rise value(v2)
  time=tr+0.5*pw;
  time_op2->setValue(time);
  EXPECT_DOUBLE_EQ(pulse.val(), v2);
}

TEST ( Double_Ast_Spice_Src_Test, spicePulseOp_breakPoints )
{  
  // these numbers are from the VPULSE regression test:
  //0V 1V 0S 10US 10US 0.1US 20.1US
  double v1=0.0, v2=1.0, td=0.0, tr=10.0e-6, tf=10.0e-6, pw=0.1e-6, per=20.1e-6, time = 0.0;

  std::vector<double> bpTestVec = { 0.0, tr, tr+pw, tr+pw+tf, per, per+tr, per+tr+pw, per+tr+pw+tf, per+per };

#if 1
  {
    // If this block of code is commented out, then the test will fail 
    std::sort ( bpTestVec.begin(), bpTestVec.end() );
    std::vector<double>::iterator it = std::unique ( bpTestVec.begin(), bpTestVec.end() );
    bpTestVec.resize( std::distance (bpTestVec.begin(), it ));
  }
#endif

  RCP<astNode<double> > v1_op = rcp(new numval<double> (v1));
  RCP<astNode<double> > v2_op = rcp(new numval<double> (v2));
  RCP<astNode<double> > td_op = rcp(new numval<double> (td));
  RCP<astNode<double> > tr_op = rcp(new numval<double> (tr));
  RCP<astNode<double> > tf_op = rcp(new numval<double> (tf));
  RCP<astNode<double> > pw_op = rcp(new numval<double> (pw));
  RCP<astNode<double> > per_op = rcp(new numval<double> (per));
  RCP<astNode<double> > time_op = rcp(new specialsOp<double> (std::string("time")));
  RCP<specialsOp<double> > time_op2 = rcp_static_cast<specialsOp<double> > (time_op);

  std::vector< RCP<astNode<double> > > sourceArgs = { v1_op, v2_op, td_op, tr_op, tf_op, pw_op, per_op };

  // test initial value(v1)
  time_op2->setValue(time);
  spicePulseOp<double> pulse( sourceArgs, time_op);

  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  bool ret = pulse.getBreakPoints(breakPointTimes);

  {
    double tol=1.0e-20;
    Xyce::Util::BreakPointLess breakPointLess(tol);
    Xyce::Util::BreakPointEqual breakPointEqual(tol);
    std::sort ( breakPointTimes.begin(), breakPointTimes.end(), breakPointLess );
    std::vector<Xyce::Util::BreakPoint>::iterator it = std::unique ( breakPointTimes.begin(), breakPointTimes.end(), breakPointEqual );
    breakPointTimes.resize( std::distance (breakPointTimes.begin(), it ));
  }

  std::vector<double> bpTimes(breakPointTimes.size());
  EXPECT_EQ(breakPointTimes.size(), bpTestVec.size());
  for (int ii=0;ii<breakPointTimes.size();ii++) 
  { 
    bpTimes[ii] = breakPointTimes[ii].value(); 
    EXPECT_DOUBLE_EQ(bpTimes[ii],bpTestVec[ii]);
  }
}

TEST ( Double_Ast_Spice_Src_Test, spiceSinOp )
{
  // these numbers are from the  VSIN/bug1679 regression test:
  // SIN ( 1.65 1.65 10000 0 0 -90 )
  //
  double v0(1.65), va(1.65), freq(10000), td(0.0), theta(0.0), phase(-90), time(0.0);

  RCP<astNode<double> > v0_op = rcp(new numval<double> (v0));
  RCP<astNode<double> > va_op = rcp(new numval<double> (va));
  RCP<astNode<double> > freq_op = rcp(new numval<double> (freq));
  RCP<astNode<double> > td_op = rcp(new numval<double> (td));
  RCP<astNode<double> > theta_op = rcp(new numval<double> (theta));
  RCP<astNode<double> > phase_op = rcp(new numval<double> (phase));

  std::vector< RCP<astNode<double> > > sourceArgs = { v0_op, va_op, freq_op, td_op, theta_op, phase_op };

  RCP<astNode<double> > time_op = rcp(new specialsOp<double> (std::string("time")));
  RCP<specialsOp<double> > time_op2 = rcp_static_cast<specialsOp<double> > (time_op);

  // test the DCOP value , which should be: DCOPValue = V0 + VA * sin (2.0*mpi*(PHASE/360));
  time_op2->setValue(time);
  spiceSinOp<double> sinOp (sourceArgs, time_op);

  double DCOPValue = v0 + va * std::sin (2.0*M_PI*(phase/360));
  EXPECT_DOUBLE_EQ(sinOp.val(), DCOPValue);

  // test a transient value, which should be 
  // TRANValue = (v0) + (va) * std::sin(2.0*mpi*((freq)*time + (phase)/360)) * std::exp( -(time*(theta)));
  time=1.0/10000.0;
  time_op2->setValue(time);
  double TRANValue = v0 + va * std::sin(2.0*M_PI*((freq)*time + (phase)/360)) * std::exp( -(time*(theta)));
  EXPECT_DOUBLE_EQ(sinOp.val(), TRANValue);
}

TEST ( Double_Ast_Spice_Src_Test, spiceExpOp )
{  
  // these numbers are from the   SOURCES/sources.cir
  //  V1=1.1 V2=2 TD1=2ns TAU1=15ns TD2=5ns TAU2=30ns
  //
  double v1(1.1), v2(2.0), td1(2e-9), tau1(15e-9), td2(5e-9), tau2(30e-9), time(0.0);

  RCP<astNode<double> > v1_op = rcp(new numval<double> (v1));
  RCP<astNode<double> > v2_op = rcp(new numval<double> (v2));
  RCP<astNode<double> > td1_op = rcp(new numval<double> (td1));
  RCP<astNode<double> > tau1_op = rcp(new numval<double> (tau1));
  RCP<astNode<double> > td2_op = rcp(new numval<double> (td2));
  RCP<astNode<double> > tau2_op = rcp(new numval<double> (tau2));

  std::vector< RCP<astNode<double> > > sourceArgs = { v1_op, v2_op, td1_op, tau1_op, td2_op, tau2_op };

  RCP<astNode<double> > time_op = rcp(new specialsOp<double> (std::string("time")));
  RCP<specialsOp<double> > time_op2 = rcp_static_cast<specialsOp<double> > (time_op);

  // test time <= td1, which should be v1
  time_op2->setValue(time);
  spiceExpOp<double> expOp ( sourceArgs, time_op );
  EXPECT_DOUBLE_EQ(expOp.val(), v1);

  // test td1 < time <= td2, which should be 
  // value = v1 + (v2-v1)*(1.0-std::exp(-(time-td1)/tau1));
  time=(td1+td2)*0.5;
  time_op2->setValue(time);
  double value = v1 + (v2-v1)*(1.0-std::exp(-(time-td1)/tau1));
  EXPECT_DOUBLE_EQ(expOp.val(), value);

  // test time > td2, which should be 
  // value = v1 + (v2-v1)*(1.0-std::exp(-(time-td1)/tau1)) + (v1-v2)*(1.0-std::exp(-(time-td2)/tau2)) ;
  time=td2*1.1;
  time_op2->setValue(time);
  value = v1 + (v2-v1)*(1.0-std::exp(-(time-td1)/tau1)) + (v1-v2)*(1.0-std::exp(-(time-td2)/tau2)) ;
  EXPECT_DOUBLE_EQ(expOp.val(), value);
}

TEST ( Double_Ast_Spice_Src_Test, spiceSffmOp )
{
  // these numbers are from the   SOURCES/sources.cir
  //  V0=-0.5 VA=2 FC=100meg MDI=0.3 FS=2.1meg
  //
  double v0(-0.5), va(2.0), fc(100e6), mdi(0.3), fs(2.1e6), time(0.0);
  //numval<double> v0_op(v0), va_op(va), fc_op(fc), mdi_op(mdi), fs_op(fs);
  //specialsOp<double> time_op(std::string("time"));

  RCP<astNode<double> > v0_op = rcp(new numval<double> (v0));
  RCP<astNode<double> > va_op = rcp(new numval<double> (va));
  RCP<astNode<double> > fc_op = rcp(new numval<double> (fc));
  RCP<astNode<double> > mdi_op = rcp(new numval<double> (mdi));
  RCP<astNode<double> > fs_op = rcp(new numval<double> (fs));
  
  std::vector< RCP<astNode<double> > > sourceArgs = { v0_op, va_op, fc_op, mdi_op, fs_op };

  RCP<astNode<double> > time_op = rcp(new specialsOp<double> (std::string("time")));
  RCP<specialsOp<double> > time_op2 = 
    rcp_static_cast<specialsOp<double> > (time_op);
  
  spiceSffmOp<double> sffmOp ( sourceArgs, time_op );

  // test, which should be    
  // value = v0 + va * sin((2 * mpi * fc * time) + mdi * sin (2 * mpi * fs * time));
  time=0.1;
  time_op2->setValue(time);
  double value = v0 + va * sin((2 * M_PI * fc * time) + mdi * sin (2 * M_PI * fs * time));
  EXPECT_NEAR(sffmOp.val(), value, 2e-8);
}

//-------------------------------------------------------------------------------
// derivatives of binary operators
TEST ( Double_Ast_Deriv_Test, test1)
{
  RCP<astNode<double> > val1 = rcp(new numval<double> (2.0));
  RCP<astNode<double> > val2 = rcp(new numval<double> (5.0));
  RCP<astNode<double> > arg1 = rcp(new paramOp<double> ( std::string("A"))); 
  RCP<astNode<double> > arg2 = rcp(new paramOp<double> ( std::string("B")));
  arg1->setNode(val1);
  arg2->setNode(val2);
  RCP<astNode<double> > binaryAddOp_1 = rcp(new binaryAddOp<double>( arg1, arg2));

  // value
  EXPECT_DOUBLE_EQ(binaryAddOp_1->val(),(2.0 + 5.0));

  // derivs
  arg1->setDerivIndex(0);  
  arg2->setDerivIndex(1);
  arg1->setIsVar();
  arg2->setIsVar();
  EXPECT_DOUBLE_EQ(binaryAddOp_1->dx(0),1.0);
  EXPECT_DOUBLE_EQ(binaryAddOp_1->dx(1),1.0);
}

#define AST_BINARY_DERIV_TEST_MACRO(NAME,OP,CPPOP, VAL1, VAL2, D1, D2) \
TEST ( NAME, OP ) \
{ \
  RCP<astNode<double> > val1 = rcp(new numval<double> (VAL1)); \
  RCP<astNode<double> > val2 = rcp(new numval<double> (VAL2)); \
  RCP<astNode<double> > arg1 = rcp(new paramOp<double> (std::string("A")));  \
  arg1->setNode(val1); \
  RCP<astNode<double> > arg2 = rcp(new paramOp<double> (std::string("B"))); \
  arg2->setNode(val2); \
  RCP<astNode<double> > OP_1 = rcp(new OP<double> (arg1,arg2)); \
  EXPECT_DOUBLE_EQ(OP_1->val(),(VAL1 CPPOP VAL2)); \
  arg1->setDerivIndex(0); \
  arg2->setDerivIndex(1); \
  arg1->setIsVar(); \
  arg2->setIsVar(); \
  EXPECT_DOUBLE_EQ(OP_1->dx(0),D1); \
  EXPECT_DOUBLE_EQ(OP_1->dx(1),D2); \
}

#define AST_BINARY_DERIV_TEST_MACRO_CMPLX(NAME,OP,CPPOP, VAL1, VAL2, D1, D2) \
TEST ( NAME, OP ) \
{ \
  RCP<astNode<cmplx> > val1 = rcp(new numval<cmplx> (VAL1)); \
  RCP<astNode<cmplx> > val2 = rcp(new numval<cmplx> (VAL2)); \
  RCP<astNode<cmplx> > arg1 = rcp(new paramOp<cmplx> (std::string("A")));  \
  arg1->setNode(val1); \
  RCP<astNode<cmplx> > arg2 = rcp(new paramOp<cmplx> (std::string("B"))); \
  arg2->setNode(val2); \
  RCP<astNode<cmplx> > OP_1 = rcp(new OP<cmplx> (arg1,arg2)); \
  EXPECT_DOUBLE_EQ(std::real(OP_1->val()),std::real((VAL1 CPPOP VAL2))); \
  EXPECT_DOUBLE_EQ(std::imag(OP_1->val()),std::imag((VAL1 CPPOP VAL2))); \
  arg1->setDerivIndex(0); \
  arg2->setDerivIndex(1); \
  arg1->setIsVar(); \
  arg2->setIsVar(); \
  EXPECT_DOUBLE_EQ(std::real(OP_1->dx(0)),std::real(D1)); \
  EXPECT_DOUBLE_EQ(std::imag(OP_1->dx(0)),std::imag(D1)); \
  EXPECT_DOUBLE_EQ(std::real(OP_1->dx(1)),std::real(D2)); \
  EXPECT_DOUBLE_EQ(std::imag(OP_1->dx(1)),std::imag(D2)); \
} 

AST_BINARY_DERIV_TEST_MACRO(Double_Ast_Deriv_Test,binaryAddOp,+,3.0,4.0,1.0,1.0)
AST_BINARY_DERIV_TEST_MACRO(Double_Ast_Deriv_Test,binaryMinusOp,-,3.0,4.0,1.0,-1.0)
AST_BINARY_DERIV_TEST_MACRO(Double_Ast_Deriv_Test,binaryMulOp,*,3.0,4.0,4.0,3.0)
AST_BINARY_DERIV_TEST_MACRO(Double_Ast_Deriv_Test,binaryDivOp,/,3.0,4.0,0.25,(-3.0/16.0))

AST_BINARY_DERIV_TEST_MACRO_CMPLX(Complex_Ast_Deriv_Test,binaryAddOp,+,cmplx(3.0,1.0),cmplx(4.0,2.0),cmplx(1.0,0.0),cmplx(1.0,0.0))
AST_BINARY_DERIV_TEST_MACRO_CMPLX(Complex_Ast_Deriv_Test,binaryMinusOp,-,cmplx(3.0,1.0),cmplx(4.0,2.0),cmplx(1.0,0.0),cmplx(-1.0,0.0))
AST_BINARY_DERIV_TEST_MACRO_CMPLX(Complex_Ast_Deriv_Test,binaryMulOp,*,cmplx(3.0,1.0),cmplx(4.0,2.0),cmplx(4.0,2.0),cmplx(3.0,1.0))
AST_BINARY_DERIV_TEST_MACRO_CMPLX(Complex_Ast_Deriv_Test,binaryDivOp,/,cmplx(3.0,1.0),cmplx(4.0,2.0), +cmplx(4.0,2.0)/(cmplx(4.0,2.0)*cmplx(4.0,2.0)), -cmplx(3.0,1.0)/(cmplx(4.0,2.0)*cmplx(4.0,2.0)))

//-------------------------------------------------------------------------------
// derivatives of unary std functions
#define AST_UNARY_DERIV_TEST_MACRO(TYPE,NAME,OP,CPPFUNC,VAL1,D1) \
TEST ( NAME, OP ) \
{ \
  RCP<astNode<TYPE> > val1 = rcp(new numval<TYPE> (VAL1)); \
  RCP<astNode<TYPE> > arg1 = rcp(new paramOp<TYPE> ( std::string("A")));  \
  arg1->setNode(val1); \
  RCP<astNode<TYPE> > OP_1 = rcp(new OP<TYPE> (arg1)); \
  EXPECT_DOUBLE_EQ(OP_1->val(),CPPFUNC(VAL1)); \
  arg1->setDerivIndex(0); \
  arg1->setIsVar(); \
  EXPECT_DOUBLE_EQ( (OP_1->dx(0)-D1), 0.0); \
}

#define AST_UNARY_DERIV_TEST_MACRO2(TYPE,NAME,OP,CPPFUNC,VAL1,D1) \
TEST ( NAME, OP ) \
{ \
  RCP<astNode<TYPE> > val1 = rcp(new numval<TYPE> (VAL1)); \
  RCP<astNode<TYPE> > arg1 = rcp(new paramOp<TYPE> ( std::string("A")));  \
  arg1->setNode(val1); \
  RCP<astNode<TYPE> > OP_1 = rcp(new OP<TYPE> (arg1)); \
  EXPECT_DOUBLE_EQ(std::real(OP_1->val()),std::real(CPPFUNC(VAL1))); \
  EXPECT_DOUBLE_EQ(std::imag(OP_1->val()),std::imag(CPPFUNC(VAL1))); \
  arg1->setDerivIndex(0); \
  arg1->setIsVar(); \
  EXPECT_DOUBLE_EQ( std::real(OP_1->dx(0)-D1), 0.0); \
  EXPECT_DOUBLE_EQ( std::imag(OP_1->dx(0)-D1), 0.0); \
}

#define AST_UNARY_DERIV_TEST_MACRO3(TYPE,NAME,OP,CPPFUNC,VAL1,D1,TOL) \
TEST ( NAME, OP ) \
{ \
  RCP<astNode<TYPE> > val1 = rcp(new numval<TYPE> (VAL1)); \
  RCP<astNode<TYPE> > arg1 = rcp(new paramOp<TYPE> ( std::string("A")));  \
  arg1->setNode(val1); \
  RCP<astNode<TYPE> > OP_1 = rcp(new OP<TYPE> (arg1)); \
  EXPECT_NEAR(std::real(OP_1->val()),std::real(CPPFUNC(VAL1)), (TOL)); \
  EXPECT_NEAR(std::imag(OP_1->val()),std::imag(CPPFUNC(VAL1)),(TOL)); \
  arg1->setDerivIndex(0); \
  arg1->setIsVar(); \
  EXPECT_NEAR( std::real(OP_1->dx(0)-D1), 0.0, (TOL)); \
  EXPECT_NEAR( std::imag(OP_1->dx(0)-D1), 0.0, (TOL)); \
}

// complex
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_Ops, sqrtOp, std::sqrt, 4.0, (0.5/std::sqrt(4.0)))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_Ops, expOp, std::exp, 0.5, std::exp(0.5))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_Ops, absOp, std::abs, -0.5, -1.0)
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_Ops, sinOp, std::sin, 0.5, std::cos(0.5))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_Ops, cosOp, std::cos, 0.5, -std::sin(0.5))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_Ops, tanOp, std::tan, 0.5, (1.0+(std::tan(0.5)*std::tan(0.5))))

AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_Ops, asinOp, std::asin, 0.5,(+1.0/std::sqrt(1.0-0.5*0.5)))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_Ops, acosOp, std::acos, 0.5,(-1.0/std::sqrt(1.0-0.5*0.5)))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_Ops, atanOp, std::atan, 0.5,(+1.0/(1+0.5*0.5)))

AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_Ops, sinhOp, std::sinh, 0.5, std::cosh(0.5))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_Ops, coshOp, std::cosh, 0.5, std::sinh(0.5))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_Ops, tanhOp, std::tanh, 0.1, 1.0/(std::cosh(0.1)*std::cosh(0.1)))

AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_Ops, asinhOp, std::asinh, 0.5, 1.0/(std::sqrt(1+0.5*0.5)))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_Ops, acoshOp, std::acosh, 1.5, 1.0/(std::sqrt( (1.5-1) * (1.5+1) )))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_Ops, atanhOp, std::atanh, 0.5, 1.0/(1.0-0.5*0.5))

AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_Ops, logOp, std::log, 0.5, 1.0/0.5)
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_Ops, log10Op, std::log10, 0.5, 1.0/(0.5*std::log(10)))

// complex
AST_UNARY_DERIV_TEST_MACRO2(cmplx, Complex_UnaryDeriv_Ast_Ops, sqrtOp, std::sqrt, cmplx(4.0,0.2), (0.5/std::sqrt(cmplx(4.0,0.2))))
AST_UNARY_DERIV_TEST_MACRO2(cmplx, Complex_UnaryDeriv_Ast_Ops, expOp, std::exp, cmplx(0.5,0.2), std::exp(cmplx(0.5,0.2)))
AST_UNARY_DERIV_TEST_MACRO2(cmplx, Complex_UnaryDeriv_Ast_Ops, absOp, std::abs, cmplx(-0.5,0.2), -1.0)
AST_UNARY_DERIV_TEST_MACRO2(cmplx, Complex_UnaryDeriv_Ast_Ops, sinOp, std::sin, cmplx(0.5,0.2), std::cos(cmplx(0.5,0.2)))
AST_UNARY_DERIV_TEST_MACRO2(cmplx, Complex_UnaryDeriv_Ast_Ops, cosOp, std::cos, cmplx(0.5,0.2), -std::sin(cmplx(0.5,0.2)))
AST_UNARY_DERIV_TEST_MACRO2(cmplx, Complex_UnaryDeriv_Ast_Ops, tanOp, std::tan, cmplx(0.5,0.2), (1.0+std::tan(cmplx(0.5,0.2))*std::tan(cmplx(0.5,0.2))))

AST_UNARY_DERIV_TEST_MACRO2(cmplx, Complex_UnaryDeriv_Ast_Ops, asinOp, std::asin, cmplx(0.5,0.2),(+1.0/std::sqrt(1.0-cmplx(0.5,0.2)*cmplx(0.5,0.2))))
AST_UNARY_DERIV_TEST_MACRO2(cmplx, Complex_UnaryDeriv_Ast_Ops, acosOp, std::acos, cmplx(0.5,0.2),(-1.0/std::sqrt(1.0-cmplx(0.5,0.2)*cmplx(0.5,0.2))))
AST_UNARY_DERIV_TEST_MACRO2(cmplx, Complex_UnaryDeriv_Ast_Ops, atanOp, std::atan, cmplx(0.5,0.2),(+1.0/(1.0+cmplx(0.5,0.2)*cmplx(0.5,0.2))))

AST_UNARY_DERIV_TEST_MACRO2(cmplx, Complex_UnaryDeriv_Ast_Ops, sinhOp, std::sinh, cmplx(0.5,0.2), std::cosh(cmplx(0.5,0.2)))
AST_UNARY_DERIV_TEST_MACRO2(cmplx, Complex_UnaryDeriv_Ast_Ops, coshOp, std::cosh, cmplx(0.5,0.2), std::sinh(cmplx(0.5,0.2)))
AST_UNARY_DERIV_TEST_MACRO2(cmplx, Complex_UnaryDeriv_Ast_Ops, tanhOp, std::tanh, cmplx(0.1,0.2), 1.0/(std::cosh(cmplx(0.1,0.2))*std::cosh(cmplx(0.1,0.2))))

AST_UNARY_DERIV_TEST_MACRO2(cmplx, Complex_UnaryDeriv_Ast_Ops, asinhOp, std::asinh, cmplx(0.5,0.2), 1.0/(std::sqrt(1.0+cmplx(0.5,0.2)*cmplx(0.5,0.2))))
AST_UNARY_DERIV_TEST_MACRO3(cmplx, Complex_UnaryDeriv_Ast_Ops, acoshOp, std::acosh, cmplx(1.5,0.2), (1.0/(std::sqrt( (cmplx(1.5,0.2)-1.0) * (cmplx(1.5,0.2)+1.0) ))), 1e-15)
AST_UNARY_DERIV_TEST_MACRO2(cmplx, Complex_UnaryDeriv_Ast_Ops, atanhOp, std::atanh, cmplx(0.5,0.2), 1.0/(1.0-cmplx(0.5,0.2)*cmplx(0.5,0.2)))

AST_UNARY_DERIV_TEST_MACRO3(cmplx, Complex_UnaryDeriv_Ast_Ops, logOp, std::log, cmplx(0.5,0.2), 1.0/cmplx(0.5,0.2), 1e-15)
AST_UNARY_DERIV_TEST_MACRO2(cmplx, Complex_UnaryDeriv_Ast_Ops, log10Op, std::log10, cmplx(0.5,0.2), 1.0/(cmplx(0.5,0.2)*std::log(10)))

// testing cosh and acosh for -10  (see the ascth.cir test)

TEST ( Double_UnaryDeriv_Ast_Ops , coshOp_ascth_test ) 
{
  RCP<astNode<double> > val1 = rcp(new numval<double> (-10.0));
  RCP<astNode<double> > arg1 = rcp(new paramOp<double> ( std::string("A")));
  arg1->setNode(val1);
  RCP<astNode<double> > cosh_1 = rcp(new coshOp<double> (arg1));
  EXPECT_DOUBLE_EQ(cosh_1->val(), std::cosh(-10.0)); 
  arg1->setDerivIndex(0);
  arg1->setIsVar();
  EXPECT_DOUBLE_EQ( (cosh_1->dx(0)-std::sinh(-10.0)), 0.0);
}

TEST ( Double_UnaryDeriv_Ast_Ops , acoshOp_ascth_test ) 
{
  double tmp = -10.0;
  double input = std::cosh(-10.0);
  RCP<astNode<double> > val1 = rcp(new numval<double> (input));
  RCP<astNode<double> > arg1 = rcp(new paramOp<double> ( std::string("A")));
  arg1->setNode(val1);
  RCP<astNode<double> > acosh_1 = rcp(new acoshOp<double> (arg1));
  EXPECT_DOUBLE_EQ(acosh_1->val(), std::acosh(input)); 
  arg1->setDerivIndex(0);
  arg1->setIsVar();
  double deriv = 1.0/(std::sqrt( (input-1) * (input+1) ));
  EXPECT_DOUBLE_EQ( (acosh_1->dx(0)-deriv), 0.0 );
}

//-------------------------------------------------------------------------------
// pow function
TEST ( Double_Ast_Deriv_Test, powOp )
{  
  double A=7.0;
  double B=4.0;
  RCP<astNode<double> > val1 = rcp(new numval<double> (A));
  RCP<astNode<double> > val2 = rcp(new numval<double> (B));
  RCP<astNode<double> > arg1 = rcp(new paramOp<double> (std::string("A")));
  arg1->setNode(val1);
  RCP<astNode<double> > arg2 = rcp(new paramOp<double> (std::string("B")));
  arg2->setNode(val2);

  // value
  RCP<astNode<double> > testPow = rcp(new powOp<double> (arg1,arg2));
  EXPECT_DOUBLE_EQ(testPow->val(), std::pow(A,B));  

  // derivative
  arg1->setDerivIndex(0);
  arg2->setDerivIndex(1);
  arg1->setIsVar();
  arg2->setIsVar();
  EXPECT_NEAR(testPow->dx(0)-((B/A)*std::pow(A,(B))), 0.0, 5.0e-13 );
  EXPECT_NEAR(testPow->dx(1)-(std::log(A))*std::pow(A,B), 0.0, 5.0e-13 );
}

TEST ( Double_Ast_Deriv_Test, pwrsOp )
{  
  double A=7.0;
  double B=4.0;
  RCP<astNode<double> > val1 = rcp(new numval<double> (A));
  RCP<astNode<double> > val2 = rcp(new numval<double> (B));
  RCP<astNode<double> > arg1 = rcp(new paramOp<double> (std::string("A")));
  arg1->setNode(val1);
  RCP<astNode<double> > arg2 = rcp(new paramOp<double> (std::string("B")));
  arg2->setNode(val2);

  // value
  RCP<astNode<double> > testPwrs = rcp(new pwrsOp<double> (arg1,arg2));
  EXPECT_DOUBLE_EQ(testPwrs->val(), std::pow(A,B));  

  // derivative
  arg1->setDerivIndex(0);
  arg2->setDerivIndex(1);
  arg1->setIsVar();
  arg2->setIsVar();
  EXPECT_NEAR(testPwrs->dx(0)-((B/A)*std::pow(A,(B))), 0.0, 5.0e-13 );
  EXPECT_NEAR(testPwrs->dx(1)-(std::log(A))*std::pow(A,B), 0.0, 5.0e-13 );
}

TEST ( Double_Ast_Deriv_Test, pwrsOp2 )
{
  double A=0.0;
  double B=4.0;
  RCP<astNode<double> > val1 = rcp(new numval<double> (A));
  RCP<astNode<double> > val2 = rcp(new numval<double> (B));
  RCP<astNode<double> > arg1 = rcp(new paramOp<double> (std::string("A")));
  arg1->setNode(val1);
  RCP<astNode<double> > arg2 = rcp(new paramOp<double> (std::string("B")));
  arg2->setNode(val2);

  // value
  RCP<astNode<double> > testPwrs = rcp(new pwrsOp<double> (arg1,arg2));
  EXPECT_DOUBLE_EQ(testPwrs->val(), std::pow(A,B));  

  // derivative
  arg1->setDerivIndex(0);
  arg2->setDerivIndex(1);
  arg1->setIsVar();
  arg2->setIsVar();
  EXPECT_DOUBLE_EQ(testPwrs->dx(0)-0.0, 0.0 );
  EXPECT_DOUBLE_EQ(testPwrs->dx(1)-0.0, 0.0 );
}

TEST ( Double_Ast_Deriv_Test, pwrsOp3 )
{  
  double A=-7.0;
  double B=4.0;
  RCP<astNode<double> > val1 = rcp(new numval<double> (A));
  RCP<astNode<double> > val2 = rcp(new numval<double> (B));
  RCP<astNode<double> > arg1 = rcp(new paramOp<double> (std::string("A")));
  arg1->setNode(val1);
  RCP<astNode<double> > arg2 = rcp(new paramOp<double> (std::string("B")));
  arg2->setNode(val2);

  RCP<astNode<double> > testPwrs = rcp(new pwrsOp<double> (arg1,arg2));

  // test expression, using powOp
  double negA=+7.0;
  RCP<astNode<double> > negval1 = rcp(new numval<double> (negA));
  RCP<astNode<double> > negarg1 = rcp(new paramOp<double> (std::string("A")));
  negarg1->setNode(negval1);
  RCP<astNode<double> > refPow = rcp(new powOp<double> (negarg1,arg2));
  RCP<astNode<double> > negPow = rcp(new unaryMinusOp<double> (refPow));

  // value
  EXPECT_DOUBLE_EQ(testPwrs->val(), negPow->val());

  // derivative
  negarg1->setDerivIndex(0);
  arg1->setDerivIndex(0);
  arg2->setDerivIndex(1);
  negarg1->setIsVar();
  arg1->setIsVar();
  arg2->setIsVar();
  EXPECT_DOUBLE_EQ(testPwrs->dx(0),negPow->dx(0));
  EXPECT_DOUBLE_EQ(testPwrs->dx(1),negPow->dx(1));
}

TEST ( Double_Ast_Deriv_Test, pwrsOp4 )
{  
  double A=-7.0;
  double B=4.0;
  RCP<astNode<double> > val1 = rcp(new numval<double> (A));
  RCP<astNode<double> > val2 = rcp(new numval<double> (B));
  RCP<astNode<double> > arg1 = rcp(new paramOp<double> (std::string("A")));
  arg1->setNode(val1);

  RCP<astNode<double> > testPwrs = rcp(new pwrsOp<double> (arg1,val2));

  // test expression, using powOp
  double negA=+7.0;
  RCP<astNode<double> > negval1 = rcp(new numval<double> (negA));
  RCP<astNode<double> > negarg1 = rcp(new paramOp<double> (std::string("A")));
  negarg1->setNode(negval1);
  RCP<astNode<double> > refPow = rcp(new powOp<double> (negarg1,val2));
  RCP<astNode<double> > negPow = rcp(new unaryMinusOp<double> (refPow));

  // value
  EXPECT_DOUBLE_EQ(testPwrs->val(), negPow->val());

  // derivative
  negarg1->setDerivIndex(0);
  arg1->setDerivIndex(0);
  negarg1->setIsVar();
  arg1->setIsVar();
  EXPECT_DOUBLE_EQ(testPwrs->dx(0),negPow->dx(0));
}

TEST ( Double_Ast_Deriv_Test, pwrsOp5 )
{
  double A=-7.0;
  double B=-4.0;
  RCP<astNode<double> > val1 = rcp(new numval<double> (A));
  RCP<astNode<double> > val2 = rcp(new numval<double> (B));
  RCP<astNode<double> > arg1 = rcp(new paramOp<double> (std::string("A")));
  arg1->setNode(val1);

  RCP<astNode<double> > testPwrs = rcp(new pwrsOp<double> (arg1,val2));

  // test expression, using powOp
  double negA=+7.0;
  RCP<astNode<double> > negval1 = rcp(new numval<double> (negA));
  RCP<astNode<double> > negarg1 = rcp(new paramOp<double> (std::string("A")));
  negarg1->setNode(negval1);
  RCP<astNode<double> > refPow = rcp(new powOp<double> (negarg1,val2));
  RCP<astNode<double> > negPow = rcp(new unaryMinusOp<double> (refPow));

  // value
  EXPECT_DOUBLE_EQ(testPwrs->val(), negPow->val());

  // derivative
  negarg1->setDerivIndex(0);
  arg1->setDerivIndex(0);
  negarg1->setIsVar();
  arg1->setIsVar();
  EXPECT_DOUBLE_EQ(testPwrs->dx(0),negPow->dx(0));
}

TEST ( Double_Ast_Deriv_Test, pwrsOp6 )
{  
  double A=-7.0;
  double B=-4.0;
  RCP<astNode<double> > val1 = rcp(new numval<double> (A));
  RCP<astNode<double> > val2 = rcp(new numval<double> (B));
  RCP<astNode<double> > arg2 = rcp(new paramOp<double> (std::string("B")));
  arg2->setNode(val2);

  RCP<astNode<double> > testPwrs = rcp(new pwrsOp<double> (val1,arg2));

  // test expression, using powOp
  double negA=+7.0;
  RCP<astNode<double> > negval1 = rcp(new numval<double> (negA));
  RCP<astNode<double> > refPow = rcp(new powOp<double> (negval1,arg2));
  RCP<astNode<double> > negPow = rcp(new unaryMinusOp<double> (refPow));

  // value
  EXPECT_DOUBLE_EQ(testPwrs->val(), negPow->val());

  // derivative
  arg2->setDerivIndex(0);
  arg2->setIsVar();
  EXPECT_DOUBLE_EQ(testPwrs->dx(0),negPow->dx(0));
}

TEST ( Complex_Ast_Deriv_Test, powOp )
{
  cmplx A=cmplx(7.0,2.0);
  cmplx B=cmplx(4.0,3.0);
  //numval <cmplx> val1(A), val2(B); 
  //paramOp <cmplx> arg1(std::string("A"),&val1), arg2(std::string("B"),&val2); 
  RCP<astNode<cmplx> > val1 = rcp(new numval<cmplx> (A));
  RCP<astNode<cmplx> > val2 = rcp(new numval<cmplx> (B));
  RCP<astNode<cmplx> > arg1 = rcp(new paramOp<cmplx> (std::string("A")));
  arg1->setNode(val1);
  RCP<astNode<cmplx> > arg2 = rcp(new paramOp<cmplx> (std::string("B")));
  arg2->setNode(val2);

  // value
  RCP<astNode<cmplx> > testPow = rcp(new powOp<cmplx> (arg1,arg2));
  EXPECT_DOUBLE_EQ(std::real(testPow->val()), std::real(std::pow(A,B)));
  EXPECT_DOUBLE_EQ(std::imag(testPow->val()), std::imag(std::pow(A,B)));

  // general:
  // (B.dx(i)*std::log(A.val())+B.val()*A.dx(i)/A.val())*std::pow(A.val(),B.val())
 
  // w.r.t. A:
  // (+B/A)*std::pow(A,B)
  
  // w.r.t. B:
  // (std::log(A))*std::pow(A,B)
  
  // derivative
  arg1->setDerivIndex(0);
  arg2->setDerivIndex(1);
  arg1->setIsVar();
  arg2->setIsVar();
  EXPECT_DOUBLE_EQ(std::real(testPow->dx(0)-((B/A)*std::pow(A,(B)))), 0.0 );
  EXPECT_DOUBLE_EQ(std::imag(testPow->dx(1)-(std::log(A))*std::pow(A,B)), 0.0 );
}

//-------------------------------------------------------------------------------
TEST ( Double_Ast_Param_Test, paramOp )
{
  double A=7.0;
  double B=4.0;
  RCP<astNode<double> > valA = rcp(new numval<double> (A));
  RCP<astNode<double> > valB = rcp(new numval<double> (B));

  RCP<astNode<double> > paramA = rcp(new paramOp<double> (std::string("A")));
  paramA->setNode(valA);
  RCP<astNode<double> > paramB = rcp(new paramOp<double> (std::string("B")));
  paramB->setNode(valB);

  EXPECT_DOUBLE_EQ(paramA->val()-A, 0.0);
  EXPECT_DOUBLE_EQ(paramB->val()-B, 0.0);
}

//-------------------------------------------------------------------------------
TEST ( Double_Ast_Param_Test, voltageOp1 )
{
  double A=7.0;
  double B=4.0;
  RCP<astNode<double> > valA = rcp(new numval<double> (A));
  RCP<astNode<double> > valB = rcp(new numval<double> (B));

  RCP<voltageOp<double> > voltageA = rcp(new voltageOp<double> ( std::string("A")) );
  RCP<voltageOp<double> > voltageB = rcp(new voltageOp<double> ( std::string("B")) );

  voltageA->setVal(A);
  voltageB->setVal(B);

  EXPECT_DOUBLE_EQ(voltageA->val()-A, 0.0);
  EXPECT_DOUBLE_EQ(voltageB->val()-B, 0.0);
}

//-------------------------------------------------------------------------------
// calling a .func
TEST ( Double_Ast_Func_Test, test1)
{
  // define the .func
  // This is what should result from a parser specification of .func function1(arg1,arg2) {arg1+arg2}
  RCP<astNode<double> > dummyArg1 = rcp(new paramOp<double> (std::string("arg1")));
  RCP<astNode<double> > dummyArg2 = rcp(new paramOp<double> (std::string("arg2")));
  std::vector<Teuchos::RCP<astNode<double> > > dummyBaseArgs = { dummyArg1, dummyArg2 };
  RCP<astNode<double> > function1 = rcp(new binaryAddOp<double>( dummyArg1, dummyArg2));
  std::string name = "function1";

  // define the expression that calls the .func
  // This is what should result from a parser specification of  .global_param function2 = {function1(5.0,6.0)}
  // Which should evaluate to function2 = 5.0+6.0=11.0
  RCP<astNode<double> > arg1 = rcp(new numval<double> (5.0));
  RCP<astNode<double> > arg2 = rcp(new numval<double> (6.0));
  std::vector<Teuchos::RCP<astNode<double> > > args = { arg1, arg2};
  RCP<funcOp<double> > function2 = rcp(new funcOp<double> (name, args));

  function2->setNode(function1);
  function2->setFuncArgs(dummyBaseArgs); // this does work
  EXPECT_DOUBLE_EQ(function2->val(), 11.0);
}

TEST ( Complex_Ast_Func_Test, test1)
{
  // define the .func
  // This is what should result from a parser specification of .func function1(arg1,arg2) {arg1+arg2}
  RCP<astNode<std::complex<double> > > dummyArg1 = rcp(new paramOp<std::complex<double> > (std::string("arg1")));
  RCP<astNode<std::complex<double> > > dummyArg2 = rcp(new paramOp<std::complex<double> > (std::string("arg2")));
  std::vector<Teuchos::RCP<astNode<std::complex<double> > > > dummyBaseArgs = { dummyArg1, dummyArg2};
  RCP<astNode<std::complex<double> > > function1 = rcp(new binaryAddOp<std::complex<double> >( dummyArg1, dummyArg2));

  std::string name = "function1";

  // define the expression that calls the .func
  // This is what should result from a parser specification of  .global_param function2 = {function1((5.0,2.0),(6.0,7.0))}
  // Which should evaluate to function2 = (5,2)+(6,7)=(11,9)
  RCP<astNode<std::complex<double> > > arg1 = rcp(new numval<std::complex<double> >(std::complex<double>(5.0,2.0)));
  RCP<astNode<std::complex<double> > > arg2 = rcp(new numval<std::complex<double> >(std::complex<double>(6.0,7.0)));
  std::vector<Teuchos::RCP<astNode<std::complex<double> > > > args = { arg1, arg2 };
  RCP<funcOp<std::complex<double> > > function2 = rcp(new funcOp<std::complex<double> > (name, args));

  function2->setNode(function1);
  function2->setFuncArgs(dummyBaseArgs);
  EXPECT_DOUBLE_EQ(std::real(function2->val()), 11.0);
  EXPECT_DOUBLE_EQ(std::imag(function2->val()), 9.0);
}

TEST ( Double_Ast_Func_Test, test2)
{
  // define the .func
  // This is what should result from a parser specification of .func f1(arg1,arg2) {arg1+arg2}
  RCP<astNode<double> > dummyArg1 = rcp(new paramOp<double> (std::string("arg1")));
  RCP<astNode<double> > dummyArg2 = rcp(new paramOp<double> (std::string("arg2")));
  std::vector<Teuchos::RCP<astNode<double> > > dummyBaseArgs = { dummyArg1, dummyArg2 };
  RCP<astNode<double> > f1 = rcp(new binaryAddOp<double>(dummyArg1, dummyArg2));
  std::string name = "f1";

  // define a new func that calls the first .func
  // This is what should result from a parser specification of  .func f2(arg1,arg2) = {arg1+arg2*f1(5.0,6.0)}
  // Which should evaluate to f2 = 5.0+6.0=11.0
  RCP<astNode<double> > arg1 = rcp(new numval<double> (5.0));
  RCP<astNode<double> > arg2 = rcp(new numval<double> (6.0));
  std::vector<Teuchos::RCP<astNode<double> > > args = { arg1,arg2 };
  RCP<astNode<double> > f1call = rcp(new funcOp<double> (name, args));
  f1call->setNode(f1);
  f1call->setFuncArgs(dummyBaseArgs);

  RCP<astNode<double> > d1 = rcp(new paramOp<double> (std::string("arg1")));
  RCP<astNode<double> > d2 = rcp(new paramOp<double> (std::string("arg2")));
  std::vector<Teuchos::RCP<astNode<double> > > dummyBaseArgs2 = { d1, d2};

  RCP<astNode<double> > prod = rcp(new binaryMulOp<double>(d2, f1call));
  RCP<astNode<double> > f2 = rcp(new binaryAddOp<double>(d1, prod));
  std::string name2 = "f2";

  // define expression that calls the second func
  // This is what should result from a parser specification of  .global_param final = {f2(3.0,2.0)}
  // Which should evaluate to function2 = 3.0+2.0*f1(5.0,6.0) = 3.0+2.0*11.0 = 25.0
  RCP<astNode<double> > argFinal1 = rcp(new numval<double> (3.0));
  RCP<astNode<double> > argFinal2 = rcp(new numval<double> (2.0));
  std::vector<Teuchos::RCP<astNode<double> > > argsFinal = { argFinal1, argFinal2 };
  RCP<funcOp<double> > function2 = rcp(new funcOp<double> (name2, argsFinal));

  function2->setNode(f2);
  function2->setFuncArgs(dummyBaseArgs2);
  EXPECT_DOUBLE_EQ(function2->val(), 25.0);
}

TEST ( Complex_Ast_Func_Test, test2)
{
  // define the .func
  // This is what should result from a parser specification of .func f1(arg1,arg2) {arg1+arg2}
  RCP<astNode<std::complex<double> > > dummyArg1 = rcp(new paramOp<std::complex<double> > (std::string("arg1")));
  RCP<astNode<std::complex<double> > > dummyArg2 = rcp(new paramOp<std::complex<double> > (std::string("arg2")));
  std::vector<Teuchos::RCP<astNode<std::complex<double> > > > dummyBaseArgs = { dummyArg1, dummyArg2 };
  RCP<astNode<std::complex<double> > > f1 = rcp(new binaryAddOp<std::complex<double> >(dummyArg1, dummyArg2));
  std::string name = "f1";

  // define a new func that calls the first .func
  // This is what should result from a parser specification of  .func f2(arg1,arg2) = {arg1+arg2*f1((5,2),(6,7))}
  // Which should evaluate to f2 = (5,2)+(6,7)=(11,9)
  RCP<astNode<std::complex<double> > > arg1 = rcp(new numval<std::complex<double> > (std::complex<double> (5.0,2.0)));
  RCP<astNode<std::complex<double> > > arg2 = rcp(new numval<std::complex<double> > (std::complex<double> (6.0,7.0)));
  std::vector<Teuchos::RCP<astNode<std::complex<double> > > > args = { arg1, arg2 };
  RCP<astNode<std::complex<double> > > f1call = rcp(new funcOp<std::complex<double> > (name, args));
  f1call->setNode(f1);
  f1call->setFuncArgs(dummyBaseArgs);

  RCP<astNode<std::complex<double> > > d1 = rcp(new paramOp<std::complex<double> > (std::string("arg1")));
  RCP<astNode<std::complex<double> > > d2 = rcp(new paramOp<std::complex<double> > (std::string("arg2")));
  std::vector<Teuchos::RCP<astNode<std::complex<double> > > > dummyBaseArgs2 = {d1,d2};

  RCP<astNode<std::complex<double> > > prod = rcp(new binaryMulOp<std::complex<double> >(d2, f1call));
  RCP<astNode<std::complex<double> > > f2 = rcp(new binaryAddOp<std::complex<double> >(d1, prod));
  std::string name2 = "f2";

  // define expression that calls the second func
  // This is what should result from a parser specification of  .global_param final = {f2((3,1),(2,4))}
  // Which should evaluate to function2 = (3,1)+(2,4)*f1((5,2),(6,7)) = (3,1)+(2,4)*(11,9) 
  RCP<astNode<std::complex<double> > > argFinal1 = rcp(new numval<std::complex<double> > (std::complex<double> (3.0,1.0)));
  RCP<astNode<std::complex<double> > > argFinal2 = rcp(new numval<std::complex<double> > (std::complex<double> (2.0,4.0)));
  std::vector<Teuchos::RCP<astNode<std::complex<double> > > > argsFinal = { argFinal1, argFinal2 };
  RCP<funcOp<std::complex<double> > > function2 = rcp(new funcOp<std::complex<double> > (name2, argsFinal));

  function2->setNode(f2);
  function2->setFuncArgs(dummyBaseArgs2);
  EXPECT_DOUBLE_EQ(std::real(function2->val()), 
      std::real (std::complex<double>(3,1)+std::complex<double>(2,4)*std::complex<double>(11,9)) );
  EXPECT_DOUBLE_EQ(std::imag(function2->val()), 
      std::imag (std::complex<double>(3,1)+std::complex<double>(2,4)*std::complex<double>(11,9)) );
}

//-------------------------------------------------------------------------------
TEST ( Double_Ast_Param_Test, test1)
{
  // define the .param
  double A=3.0;
  RCP<astNode<double> > paramValA = rcp(new numval<double> (A));
  RCP<astNode<double> > testParamA = rcp(new paramOp<double> (std::string("A")));
  testParamA->setNode(paramValA);

  // define the expression that uses the .param
  double B=5.0;
  RCP<astNode<double> > paramValB = rcp(new numval<double> (B));
  RCP<astNode<double> > finalExp = rcp(new binaryAddOp<double> (paramValB,testParamA));
  EXPECT_DOUBLE_EQ(finalExp->val(), A+B);
}

TEST ( Complex_Ast_Param_Test, test1)
{
  // define the .param
  std::complex<double> A(3.0,2.0);
  RCP<astNode<std::complex<double> > > paramValA = rcp(new numval<std::complex<double> > (A));
  RCP<astNode<std::complex<double> > > testParamA = rcp(new paramOp<std::complex<double> > (std::string("A")));
  testParamA->setNode(paramValA);

  // define the expression that uses the .param
  std::complex<double> B(5.0,4.0);
  RCP<astNode<std::complex<double> > > paramValB = rcp(new numval<std::complex<double> > (B));
  RCP<astNode<std::complex<double> > > finalExp = rcp(new binaryAddOp<std::complex<double> > (paramValB,testParamA));
  EXPECT_DOUBLE_EQ(std::real(finalExp->val()), std::real(A+B));
  EXPECT_DOUBLE_EQ(std::imag(finalExp->val()), std::imag(A+B));
}

TEST ( Double_Ast_Param_Test, test2)
{
  // define the .params
  double A=3.0;
  double B=4.0;
  RCP<astNode<double> > paramValA = rcp(new numval<double> (A));
  RCP<astNode<double> > paramValB = rcp(new numval<double> (B));

  RCP<astNode<double> > testParamA = rcp(new paramOp<double> (std::string("A")));
  RCP<astNode<double> > testParamB = rcp(new paramOp<double> (std::string("B")));

  testParamA->setNode(paramValA);
  testParamB->setNode(paramValB);

  // define the expression that uses the .params
  binaryMulOp <double> finalExp(testParamA, testParamB);

  EXPECT_DOUBLE_EQ(finalExp.val(), 12.0);
}

TEST ( Complex_Ast_Param_Test, test2)
{
  // define the .params
  std::complex<double> A(3.0,5.0);
  std::complex<double> B(4.0,6.0);
  RCP<astNode<std::complex<double> > > paramValA = rcp(new numval<std::complex<double> > (A));
  RCP<astNode<std::complex<double> > > paramValB = rcp(new numval<std::complex<double> > (B));
  RCP<astNode<std::complex<double> > > testParamA = rcp(new paramOp<std::complex<double> > (std::string("A")));
  RCP<astNode<std::complex<double> > > testParamB = rcp(new paramOp<std::complex<double> > (std::string("B")));
  testParamA->setNode(paramValA);
  testParamB->setNode(paramValB);

  // define the expression that uses the .params
  binaryMulOp <std::complex<double> > finalExp(testParamA, testParamB);

  EXPECT_DOUBLE_EQ(std::real(finalExp.val()), std::real(A*B));
  EXPECT_DOUBLE_EQ(std::imag(finalExp.val()), std::imag(A*B));
}

TEST ( Double_Ast_Param_Test, test3)
{
  // define the .params
  double A=3.0;
  double B=4.0;
  double C=A*B;
  double Result = (A+B)/C;  // this should be 7/12

  RCP<astNode<double> > paramValA = rcp(new numval<double> (A));
  RCP<astNode<double> > paramValB = rcp(new numval<double> (B));

  RCP<astNode<double> > testParamA = rcp(new paramOp<double> (std::string("A")));
  RCP<astNode<double> > testParamB = rcp(new paramOp<double> (std::string("B")));
  
  testParamA->setNode(paramValA);
  testParamB->setNode(paramValB);

  // define the expression that uses the .params
  RCP<astNode<double> > prodExp = rcp(new binaryMulOp<double> (testParamA,testParamB));

  // now create another .param and assign prodExp to it
  RCP<astNode<double> > testParamC = rcp(new paramOp<double> (std::string("C")));
  testParamC->setNode(prodExp);

  // now create another expression that uses all of them
  RCP<astNode<double> > numerator = rcp(new binaryAddOp<double> (testParamA, testParamB));
  RCP<astNode<double> > finalExp = rcp(new binaryDivOp<double> (numerator, testParamC));

  // final result should be = 7/12
  EXPECT_DOUBLE_EQ(finalExp->val(), Result);
}

TEST ( Complex_Ast_Param_Test, test3)
{
  // define the .params
  std::complex<double> A(3.0,5.0);
  std::complex<double> B(4.0,6.0);
  std::complex<double> C = A*B;
  std::complex<double> Result = (A+B)/C;

  RCP<astNode<std::complex<double> > > paramValA = rcp(new numval<std::complex<double> > (A));
  RCP<astNode<std::complex<double> > > paramValB = rcp(new numval<std::complex<double> > (B));
  RCP<astNode<std::complex<double> > > testParamA = rcp(new paramOp<std::complex<double> > (std::string("A")));
  RCP<astNode<std::complex<double> > > testParamB = rcp(new paramOp<std::complex<double> > (std::string("B")));
  testParamA->setNode(paramValA);
  testParamB->setNode(paramValB);

  // define the expression that uses the .params
  RCP<astNode<std::complex<double> > > prodExp = rcp(new binaryMulOp<std::complex<double> > (testParamA,testParamB));

  // now create another .param and assign prodExp to it
  RCP<astNode<std::complex<double> > > testParamC = rcp(new paramOp<std::complex<double> > (std::string("C")));
  testParamC->setNode(prodExp);

  // now create another expression that uses all of them
  RCP<astNode<std::complex<double> > > numerator = rcp(new binaryAddOp<std::complex<double> > (testParamA, testParamB));
  RCP<astNode<std::complex<double> > > finalExp = rcp(new binaryDivOp<std::complex<double> > (numerator, testParamC));

  // final result should be = (A+B)/C
  EXPECT_DOUBLE_EQ(std::real(finalExp->val()), std::real(Result) );
  EXPECT_DOUBLE_EQ(std::imag(finalExp->val()), std::imag(Result) );
}

//-------------------------------------------------------------------------------
// test values of conditional operators and if statements
//
#define AST_IF_OP_TEST_MACRO(NAME,SUBNAME,OP, C1, C2, VAL1, VAL2, RESULT) \
TEST ( NAME, SUBNAME ) \
{ \
  RCP<astNode<double> > c1 = rcp(new numval<double> (C1)); \
  RCP<astNode<double> > c2 = rcp(new numval<double> (C2)); \
  RCP<astNode<double> > cond1 = rcp(new OP<double> (c1,c2)); \
  RCP<astNode<double> > ifVal1 = rcp(new numval<double> (VAL1)); \
  RCP<astNode<double> > ifVal2 = rcp(new numval<double> (VAL2)); \
  RCP<astNode<double> > ifStmt = rcp(new ifStatementOp<double>(cond1,ifVal1,ifVal2)); \
  EXPECT_DOUBLE_EQ(ifStmt->val(), RESULT);\
}

AST_IF_OP_TEST_MACRO(Double_Ast_if_Test,eq1,eqOp,3.0,2.0,2.0,3.0,3.0) 
AST_IF_OP_TEST_MACRO(Double_Ast_if_Test,eq2,eqOp,3.0,3.0,2.0,3.0,2.0) 

AST_IF_OP_TEST_MACRO(Double_Ast_if_Test,ne1,neOp,3.0,2.0,2.0,3.0,2.0) 
AST_IF_OP_TEST_MACRO(Double_Ast_if_Test,ne2,neOp,3.0,3.0,2.0,3.0,3.0) 

AST_IF_OP_TEST_MACRO(Double_Ast_if_Test,gt1,gtOp,3.0,2.0,2.0,3.0,2.0) 
AST_IF_OP_TEST_MACRO(Double_Ast_if_Test,gt2,gtOp,2.0,3.0,2.0,3.0,3.0) 

AST_IF_OP_TEST_MACRO(Double_Ast_if_Test,lt1,ltOp,3.0,2.0,2.0,3.0,3.0) 
AST_IF_OP_TEST_MACRO(Double_Ast_if_Test,lt2,ltOp,2.0,3.0,2.0,3.0,2.0) 

AST_IF_OP_TEST_MACRO(Double_Ast_if_Test,ge1,geOp,3.0,2.0,2.0,3.0,2.0) 
AST_IF_OP_TEST_MACRO(Double_Ast_if_Test,ge2,geOp,3.0,3.0,2.0,3.0,2.0) 
AST_IF_OP_TEST_MACRO(Double_Ast_if_Test,ge3,geOp,2.0,3.0,2.0,3.0,3.0) 

AST_IF_OP_TEST_MACRO(Double_Ast_if_Test,le1,leOp,3.0,2.0,2.0,3.0,3.0) 
AST_IF_OP_TEST_MACRO(Double_Ast_if_Test,le2,leOp,3.0,3.0,2.0,3.0,2.0) 
AST_IF_OP_TEST_MACRO(Double_Ast_if_Test,le3,leOp,2.0,3.0,2.0,3.0,2.0) 

#define AST_LOGIC_OP_TEST_MACRO(NAME,SUBNAME,OP, C1, C2, RESULT) \
TEST ( NAME, SUBNAME ) \
{ \
  RCP<astNode<double> > c1 = rcp(new numval<double> (C1)); \
  RCP<astNode<double> > c2 = rcp(new numval<double> (C2)); \
  RCP<astNode<double> > cond1 = rcp(new OP<double> (c1,c2)); \
  EXPECT_DOUBLE_EQ(cond1->val(), RESULT);\
}

AST_LOGIC_OP_TEST_MACRO(Double_Ast_logical_Test,and1,andOp, 1.0, 1.0, 1.0) 
AST_LOGIC_OP_TEST_MACRO(Double_Ast_logical_Test,and2,andOp, 1.0, 0.0, 0.0) 
AST_LOGIC_OP_TEST_MACRO(Double_Ast_logical_Test,and3,andOp, 0.0, 1.0, 0.0) 
AST_LOGIC_OP_TEST_MACRO(Double_Ast_logical_Test,and4,andOp, 0.0, 0.0, 0.0) 

AST_LOGIC_OP_TEST_MACRO(Double_Ast_logical_Test,or1,orOp, 1.0, 1.0, 1.0) 
AST_LOGIC_OP_TEST_MACRO(Double_Ast_logical_Test,or2,orOp, 1.0, 0.0, 1.0) 
AST_LOGIC_OP_TEST_MACRO(Double_Ast_logical_Test,or3,orOp, 0.0, 1.0, 1.0) 
AST_LOGIC_OP_TEST_MACRO(Double_Ast_logical_Test,or4,orOp, 0.0, 0.0, 0.0) 

AST_LOGIC_OP_TEST_MACRO(Double_Ast_logical_Test,xor1,xorOp, 1.0, 1.0, 0.0) 
AST_LOGIC_OP_TEST_MACRO(Double_Ast_logical_Test,xor2,xorOp, 1.0, 0.0, 1.0) 
AST_LOGIC_OP_TEST_MACRO(Double_Ast_logical_Test,xor3,xorOp, 0.0, 1.0, 1.0) 
AST_LOGIC_OP_TEST_MACRO(Double_Ast_logical_Test,xor4,xorOp, 0.0, 0.0, 0.0) 

TEST ( Double_Ast_logical_Test, not1 ) 
{ 
  RCP<astNode<double> > c1 = rcp(new numval<double> (1.0)); \
  RCP<astNode<double> > op = rcp(new unaryNotOp<double> (c1)); \
  EXPECT_DOUBLE_EQ(op->val(), 0.0);
}

TEST ( Double_Ast_logical_Test, not2 ) 
{ 
  RCP<astNode<double> > c1 = rcp(new numval<double> (0.0)); \
  RCP<astNode<double> > op = rcp(new unaryNotOp<double> (c1)); \
  EXPECT_DOUBLE_EQ(op->val(), 1.0);
}

//-------------------------------------------------------------------------------
// table tests
//
// adapted from break.cir
TEST ( Double_Ast_table_Test, break1)
{ 
  RCP<astNode<double> > time_op = rcp(new specialsOp<double> (std::string("time")));

  std::vector<Teuchos::RCP<astNode<double> > > args;
  args.push_back( rcp(new numval<double> (0.0)) );
  args.push_back( rcp(new numval<double> (0.0)) );
  args.push_back( rcp(new numval<double>(0.3))); 
  args.push_back( rcp(new numval<double>(0))); 
  args.push_back( rcp(new numval<double>(0.301)));
  args.push_back( rcp(new numval<double>(2)));
  args.push_back( rcp(new numval<double>(0.302)));
  args.push_back( rcp(new numval<double>(2)));
  args.push_back( rcp(new numval<double>(0.6)));
  args.push_back( rcp(new numval<double>(1)));
  args.push_back( rcp(new numval<double>(1)));
  args.push_back( rcp(new numval<double>(1)));

  std::string keyword = std::string("TABLE");
  RCP<astNode<double> > table = rcp(new tableOp<double> (keyword, time_op, args));

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<double> refRes = { 0, 0, 2, 2, 1, 1 };
  int numpoints = times.size();
  std::vector<double> result(numpoints,0.0);
  for (int ii=0;ii<numpoints;ii++)
  {
    time_op->setValue(times[ii]);
    result[ii] = table->val();
    EXPECT_DOUBLE_EQ(refRes[ii],result[ii]);
  }
}

TEST ( Double_Ast_table_Test, array1)
{
  RCP<astNode<double> > time_op = rcp(new specialsOp<double> (std::string("time")));

  std::vector<double> xa = { 0.0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<double> ya = { 0.0, 0, 2, 2, 1, 1 }; 

  std::string keyword = std::string("TABLE");
  RCP<astNode<double> > table = rcp(new tableOp<double> (keyword, time_op, xa,ya));

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<double> refRes = { 0, 0, 2, 2, 1, 1 };
  int numpoints = times.size();
  std::vector<double> result(numpoints,0.0);
  for (int ii=0;ii<numpoints;ii++)
  {
    time_op->setValue(times[ii]);
    result[ii] = table->val();
    EXPECT_DOUBLE_EQ(refRes[ii],result[ii]);
  }
}

#define RCP_NV(VAL1) rcp(new numval<double>(VAL1))

// This test is adapted from Bsrc_C1.cir.  
TEST ( Double_Ast_table_Test, array2)
{ 
  RCP<astNode<double> > time_op = rcp(new specialsOp<double> (std::string("time")));

// v1 which is a sin source, a sinewave that goes between +20 and -20
// 0, 20, 1k, -.25e-3, 0, 0
  double v0(0), va(20), freq(1000), td1(-0.25e-3), theta(0.0), phase(0), time(0.0);
  RCP<astNode<double> > v0_op = rcp(new numval<double> (v0));
  RCP<astNode<double> > va_op = rcp(new numval<double> (va));
  RCP<astNode<double> > freq_op = rcp(new numval<double> (freq));
  RCP<astNode<double> > td1_op = rcp(new numval<double> (td1));
  RCP<astNode<double> > theta_op = rcp(new numval<double> (theta));
  RCP<astNode<double> > phase_op = rcp(new numval<double> (phase));
  std::vector< RCP<astNode<double> > > sourceArgs = { v0_op, va_op, freq_op, td1_op, theta_op, phase_op };
  time_op->setValue(time);
  RCP<astNode<double> > V1op = rcp(new spiceSinOp<double>( sourceArgs, time_op ));

// v2 which is a pulse source that goes between 0 and 1.  PW is short
// v1=0, v2=1, td=0, tr=0.5us, tf=0.5us, pw=2us, per=20us
  double v1=0.0, v2=1.0, td2=0, tr=0.5e-6, tf=0.5e-6, pw=2.0e-6, per=20.0e-6;
  RCP<astNode<double> > v1_op = rcp(new numval<double> (v1));
  RCP<astNode<double> > v2_op = rcp(new numval<double> (v2));
  RCP<astNode<double> > td2_op = rcp(new numval<double> (td2));
  RCP<astNode<double> > tr_op = rcp(new numval<double> (tr));
  RCP<astNode<double> > tf_op = rcp(new numval<double> (tf));
  RCP<astNode<double> > pw_op = rcp(new numval<double> (pw));
  RCP<astNode<double> > per_op = rcp(new numval<double> (per));

  sourceArgs = { v1_op, v2_op, td2_op, tr_op, tf_op, pw_op, per_op };

  RCP<astNode<double> > V2op = rcp(new spicePulseOp<double> ( sourceArgs, time_op));

  // Set up v2 * (v1+30)/60 = leftArgExpr
  // The expression, V(2) * (V(1) + 30) / 60  has roughly the scaled 
  // shape of V(1) but is spiked/digitized
  RCP<astNode<double> > thirty = rcp(new numval<double> (30));
  RCP<astNode<double> > denom = rcp(new numval<double> (60));
  RCP<astNode<double> > firstSum = rcp(new binaryAddOp<double> (V1op,thirty));
  RCP<astNode<double> > numer = rcp(new binaryMulOp<double> (firstSum,V2op));
  RCP<astNode<double> > leftArgExpr = rcp(new binaryDivOp<double> (numer,denom));

  // set up a table using the constructor that would be used by the parser
  // (i.e. invoked by ExpressionParser.yxx)
  std::vector<Teuchos::RCP<astNode<double> > > args = {
    RCP_NV(0.0000000), RCP_NV(0),
    RCP_NV(0.0312500), RCP_NV(0),
    RCP_NV(0.0312813), RCP_NV(1),
    RCP_NV(0.0625000), RCP_NV(1),
    RCP_NV(0.0625313), RCP_NV(2),
    RCP_NV(0.0937500), RCP_NV(2),
    RCP_NV(0.0937813), RCP_NV(3),
    RCP_NV(0.1250000), RCP_NV(3),
    RCP_NV(0.1250313), RCP_NV(4),
    RCP_NV(0.1562500), RCP_NV(4),
    RCP_NV(0.1562813), RCP_NV(5),
    RCP_NV(0.1875000), RCP_NV(5),
    RCP_NV(0.1875313), RCP_NV(6),
    RCP_NV(0.2187500), RCP_NV(6),
    RCP_NV(0.2187813), RCP_NV(7),
    RCP_NV(0.2500000), RCP_NV(7),
    RCP_NV(0.2500313), RCP_NV(8),
    RCP_NV(0.2812500), RCP_NV(8),
    RCP_NV(0.2812813), RCP_NV(9),
    RCP_NV(0.3125000), RCP_NV(9),
    RCP_NV(0.3125313), RCP_NV(10),
    RCP_NV(0.3437500), RCP_NV(10),
    RCP_NV(0.3437813), RCP_NV(11),
    RCP_NV(0.3750000), RCP_NV(11),
    RCP_NV(0.3750313), RCP_NV(12),
    RCP_NV(0.4062500), RCP_NV(12),
    RCP_NV(0.4062813), RCP_NV(13),
    RCP_NV(0.4375000), RCP_NV(13),
    RCP_NV(0.4375313), RCP_NV(14),
    RCP_NV(0.4687500), RCP_NV(14),
    RCP_NV(0.4687813), RCP_NV(15),
    RCP_NV(0.5000000), RCP_NV(15),
    RCP_NV(0.5000313), RCP_NV(16),
    RCP_NV(0.5312500), RCP_NV(16),
    RCP_NV(0.5312813), RCP_NV(17),
    RCP_NV(0.5625000), RCP_NV(17),
    RCP_NV(0.5625313), RCP_NV(18),
    RCP_NV(0.5937500), RCP_NV(18),
    RCP_NV(0.5937813), RCP_NV(19),
    RCP_NV(0.6250000), RCP_NV(19),
    RCP_NV(0.6250313), RCP_NV(20),
    RCP_NV(0.6562500), RCP_NV(20),
    RCP_NV(0.6562813), RCP_NV(21),
    RCP_NV(0.6875000), RCP_NV(21),
    RCP_NV(0.6875313), RCP_NV(22),
    RCP_NV(0.7187500), RCP_NV(22),
    RCP_NV(0.7187813), RCP_NV(23),
    RCP_NV(0.7500000), RCP_NV(23),
    RCP_NV(0.7500313), RCP_NV(24),
    RCP_NV(0.7812500), RCP_NV(24),
    RCP_NV(0.7812813), RCP_NV(25),
    RCP_NV(0.8125000), RCP_NV(25),
    RCP_NV(0.8125313), RCP_NV(26),
    RCP_NV(0.8437500), RCP_NV(26),
    RCP_NV(0.8437813), RCP_NV(27),
    RCP_NV(0.8750000), RCP_NV(27),
    RCP_NV(0.8750313), RCP_NV(28),
    RCP_NV(0.9062500), RCP_NV(28),
    RCP_NV(0.9062813), RCP_NV(29),
    RCP_NV(0.9375000), RCP_NV(29),
    RCP_NV(0.9375313), RCP_NV(30),
    RCP_NV(0.9687500), RCP_NV(30),
    RCP_NV(0.9687813), RCP_NV(31),
    RCP_NV(1.0000000), RCP_NV(31)
  };

  std::string keyword = std::string("TABLE");
  RCP<astNode<double> > tableFromParser = rcp(new tableOp<double> (keyword, leftArgExpr,  args));
  RCP<astNode<double> > tableFromParserCopy = tableFromParser;

  // set up the "pure array" table
    std::vector<double> xa = { 0.0000000, 0.0312500, 0.0312813, 0.0625000, 0.0625313, 0.0937500, 0.0937813, 0.1250000, 0.1250313, 0.1562500, 0.1562813, 0.1875000, 0.1875313, 0.2187500, 0.2187813, 0.2500000, 0.2500313, 0.2812500, 0.2812813, 0.3125000, 0.3125313, 0.3437500, 0.3437813, 0.3750000, 0.3750313, 0.4062500, 0.4062813, 0.4375000, 0.4375313, 0.4687500, 0.4687813, 0.5000000, 0.5000313, 0.5312500, 0.5312813, 0.5625000, 0.5625313, 0.5937500, 0.5937813, 0.6250000, 0.6250313, 0.6562500, 0.6562813, 0.6875000, 0.6875313, 0.7187500, 0.7187813, 0.7500000, 0.7500313, 0.7812500, 0.7812813, 0.8125000, 0.8125313, 0.8437500, 0.8437813, 0.8750000, 0.8750313, 0.9062500, 0.9062813, 0.9375000, 0.9375313, 0.9687500, 0.9687813, 1.0000000 };

    std::vector<double> ya = { 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31 };

  RCP<astNode<double> > tablePureArray = rcp(new tableOp<double> (keyword, leftArgExpr, xa,ya));
  RCP<astNode<double> > tablePureArrayCopy = tablePureArray;

  int numpoints=100; 
  double tfinal = 0.0005;
  double dt = tfinal/(numpoints-1);

  std::vector<double> result(numpoints,0.0);
  std::vector<double> refRes(numpoints,0.0);
  std::vector<double> resultCopy(numpoints,0.0);
  std::vector<double> refResCopy(numpoints,0.0);

  for (int ii=0;ii<numpoints;ii++,time+=dt) 
  {
    time_op->setValue(time);
    refRes[ii] = tablePureArray->val();
    result[ii] = tableFromParser->val();
    refResCopy[ii] = tablePureArrayCopy->val();
    resultCopy[ii] = tableFromParserCopy->val();

    EXPECT_DOUBLE_EQ(refRes[ii],result[ii]);
    EXPECT_DOUBLE_EQ(refRes[ii],resultCopy[ii]);
    EXPECT_DOUBLE_EQ(refResCopy[ii],result[ii]);
  }
}

//-------------------------------------------------------------------------------
TEST ( Double_Ast_calculus_Test, ddx1)
{ 
  // define the .param
  double A=3.0;
  RCP<astNode<double> > paramValA = rcp(new numval<double> (A));
  RCP<astNode<double> > testParamA = rcp(new paramOp<double> (std::string("A")));
  testParamA->setNode(paramValA);

  // define the expression that uses the .param
  double B=5.0;
  RCP<astNode<double> > paramValB = rcp(new numval<double> (B));
  RCP<astNode<double> > f_of_x = rcp(new binaryAddOp<double> (paramValB,testParamA));

  RCP<astNode<double> > finalExp1 = rcp(new ddxOp<double> (f_of_x,testParamA));
  EXPECT_DOUBLE_EQ(finalExp1->val(), 1.0);
}

TEST ( Double_Ast_calculus_Test, ddx2)
{ 
  // define the .params
  RCP<astNode<double> > paramValA = rcp(new numval<double> (3.0));
  RCP<astNode<double> > paramValB = rcp(new numval<double> (4.0));
  RCP<astNode<double> > testParamA = rcp(new paramOp<double> (std::string("A")));
  RCP<astNode<double> > testParamB = rcp(new paramOp<double> (std::string("B")));
  testParamA->setNode(paramValA);
  testParamB->setNode(paramValB);

  // define the expression that uses the .params
  RCP<astNode<double> > f_of_x = rcp(new binaryMulOp<double> (testParamA, testParamB));

  RCP<astNode<double> > finalExp1 = rcp(new ddxOp<double> (f_of_x,testParamA));
  EXPECT_DOUBLE_EQ(finalExp1->val(), 4.0);

  RCP<astNode<double> > finalExp2 = rcp(new ddxOp<double> (f_of_x,testParamB));
  EXPECT_DOUBLE_EQ(finalExp2->val(), 3.0);
}

TEST ( Double_Ast_calculus_Test, ddx3)
{ 
  RCP<astNode<double> > paramValA = rcp(new numval<double> (3.0));
  RCP<astNode<double> > testParamA = rcp(new paramOp<double> (std::string("A")));
  testParamA->setNode(paramValA);

  RCP<astNode<double> > f_of_x = rcp(new binaryMulOp<double> (testParamA, testParamA));

  ddxOp<double> finalExp1(f_of_x, testParamA);
  EXPECT_DOUBLE_EQ(finalExp1.val(), 6.0);
}

TEST ( Double_Ast_calculus_Test, ddx4)
{ 
  RCP<astNode<double> > paramValA = rcp(new numval<double> (3.0));
  RCP<astNode<double> > testParamA = rcp(new paramOp<double> (std::string("A")));
  testParamA->setNode(paramValA);

  RCP<astNode<double> > f_of_x = rcp(new sinOp<double> (testParamA));

  ddxOp<double> finalExp1(f_of_x, testParamA);
  EXPECT_DOUBLE_EQ(finalExp1.val(), std::cos(3.0));
}

TEST ( Double_Ast_floor_Test, test1)
{ 
  //Floor of 10.25 = 10
  RCP<astNode<double> > valA = rcp(new numval<double> (10.25));
  RCP<astNode<double> > floor1 = rcp(new floorOp<double> (valA));
  EXPECT_DOUBLE_EQ(floor1->val(), 10);
}

TEST ( Double_Ast_floor_Test, test2)
{ 
  //Floor of -34.251 = -35
  RCP<astNode<double> > valA = rcp(new numval<double> (-34.251));
  RCP<astNode<double> > floor1 = rcp(new floorOp<double> (valA));
  EXPECT_DOUBLE_EQ(floor1->val(), -35);
}

TEST ( Double_Ast_floor_Test, test3)
{ 
  //Floor of 0.71 = 0
  RCP<astNode<double> > valA = rcp(new numval<double> (0.71));
  RCP<astNode<double> > floor1 = rcp(new floorOp<double> (valA));
  EXPECT_DOUBLE_EQ(floor1->val(), 0);
}

TEST ( Double_Ast_floor_Test, test4)
{ 
  //Floor of 10.25 = 10
  RCP<astNode<double> > valA = rcp(new numval<double> (10.25));
  RCP<astNode<double> > testParamA = rcp(new paramOp<double> (std::string("A")));
  testParamA->setNode(valA);
  RCP<astNode<double> > floor1 = rcp(new floorOp<double> (testParamA));

  EXPECT_DOUBLE_EQ(floor1->val(), 10);

  RCP<astNode<double> > ddx1 = rcp(new ddxOp<double> (floor1,testParamA));
  EXPECT_DOUBLE_EQ(ddx1->val(), 0);
}

TEST ( Double_Ast_ceil_Test, test1)
{ 
  //Floor of 10.25 = 10
  RCP<astNode<double> > valA = rcp(new numval<double> (10.25));
  RCP<astNode<double> > ceil1 = rcp(new ceilOp<double> (valA));
  EXPECT_DOUBLE_EQ(ceil1->val(), 11);
}

TEST ( Double_Ast_ceil_Test, test2)
{ 
  //Floor of -34.251 = -35
  RCP<astNode<double> > valA = rcp(new numval<double> (-34.251));
  RCP<astNode<double> > ceil1 = rcp(new ceilOp<double> (valA));
  EXPECT_DOUBLE_EQ(ceil1->val(), -34);
}

TEST ( Double_Ast_ceil_Test, test3)
{ 
  //Floor of 0.71 = 0
  RCP<astNode<double> > valA = rcp(new numval<double> (0.71));
  RCP<astNode<double> > ceil1 = rcp(new ceilOp<double> (valA));
  EXPECT_DOUBLE_EQ(ceil1->val(), 1);
}

TEST ( Double_Ast_ceil_Test, test4)
{ 
  //Floor of 10.25 = 10
  RCP<astNode<double> > valA = rcp(new numval<double> (10.25));
  RCP<astNode<double> > testParamA = rcp(new paramOp<double> (std::string("A")));
  testParamA->setNode(valA);
  RCP<astNode<double> > ceil1 = rcp(new ceilOp<double> (testParamA));

  EXPECT_DOUBLE_EQ(ceil1->val(), 11);
  RCP<astNode<double> > ddx1 = rcp(new ddxOp<double> (ceil1,testParamA));
  EXPECT_DOUBLE_EQ(ddx1->val(), 0);
}

TEST ( Double_Ast_sgn_Test, test1)
{ 
  // +1 if x > 0 
  //  0 if x = 0 
  // -1 if x < 0
  RCP<astNode<double> > valA = rcp(new numval<double> (10.25));
  RCP<astNode<double> > sgn1 = rcp(new sgnOp<double> (valA));
  EXPECT_DOUBLE_EQ(sgn1->val(), 1);
}

TEST ( Double_Ast_sgn_Test, test2)
{ 
// +1 if x > 0 
//  0 if x = 0 
// -1 if x < 0
  RCP<astNode<double> > valA = rcp(new numval<double> (0.0));
  RCP<astNode<double> > sgn1 = rcp(new sgnOp<double> (valA));
  EXPECT_DOUBLE_EQ(sgn1->val(), 0);
}

TEST ( Double_Ast_sgn_Test, test3)
{ 
  // +1 if x > 0 
  //  0 if x = 0 
  // -1 if x < 0
  RCP<astNode<double> > valA = rcp(new numval<double> (-7.5));
  RCP<astNode<double> > sgn1 = rcp(new sgnOp<double> (valA));
  EXPECT_DOUBLE_EQ(sgn1->val(), -1);
}

TEST ( Double_Ast_sign_Test, test1)
{ 
  // sign(x,y) = sgn(y)|x|  sign of y times absolute value of x
  RCP<astNode<double> > valX = rcp(new numval<double> (-25));
  RCP<astNode<double> > valY = rcp(new numval<double> (10.25));
  RCP<astNode<double> > sign1 = rcp(new signOp<double> (valX,valY));
  EXPECT_DOUBLE_EQ(sign1->val(), 25);
} 

TEST ( Double_Ast_sign_Test, test2)
{ 
  // sign(x,y) = sgn(y)|x|  sign of y times absolute value of x
  RCP<astNode<double> > valX = rcp(new numval<double> (15));
  RCP<astNode<double> > valY = rcp(new numval<double> (-10.25));
  RCP<astNode<double> > sign1 = rcp(new signOp<double> (valX,valY));
  EXPECT_DOUBLE_EQ(sign1->val(), -15);
}

TEST ( Double_Ast_limit_Test, test1)
{ // limit(x,y,z)
  // x limited to range y to z
  // y if x < y
  // x if y < x < z
  // z if x > z
  RCP<astNode<double> > valX = rcp(new numval<double> (-25));   // test value
  RCP<astNode<double> > valY = rcp(new numval<double> (1.25));  // lower value
  RCP<astNode<double> > valZ = rcp(new numval<double> (11.25)); // upper value
  RCP<astNode<double> > limit1 = rcp(new limitOp<double> (valX,valY,valZ));

  EXPECT_DOUBLE_EQ(limit1->val(), 1.25);
}

TEST ( Double_Ast_limit_Test, test2)
{ // limit(x,y,z)
  // x limited to range y to z
  // y if x < y
  // x if y < x < z
  // z if x > z
  RCP<astNode<double> > valX = rcp(new numval<double> (+5));   // test value
  RCP<astNode<double> > valY = rcp(new numval<double> (1.25));  // lower value
  RCP<astNode<double> > valZ = rcp(new numval<double> (11.25)); // upper value
  RCP<astNode<double> > limit1 = rcp(new limitOp<double> (valX,valY,valZ));

  EXPECT_DOUBLE_EQ(limit1->val(), 5);
}

TEST ( Double_Ast_limit_Test, test3)
{ // limit(x,y,z)
  // x limited to range y to z
  // y if x < y
  // x if y < x < z
  // z if x > z
  RCP<astNode<double> > valX = rcp(new numval<double> (+17));   // test value
  RCP<astNode<double> > valY = rcp(new numval<double> (1.25));  // lower value
  RCP<astNode<double> > valZ = rcp(new numval<double> (11.25)); // upper value
  RCP<astNode<double> > limit1 = rcp(new limitOp<double> (valX,valY,valZ));
  EXPECT_DOUBLE_EQ(limit1->val(), 11.25);
}

TEST ( Double_Ast_int_Test, test1)
{ 
  // int(x)
  // integer part of the real variable x 
  RCP<astNode<double> > valX = rcp(new numval<double> (11.2423));   // test value
  RCP<astNode<double> > int1 = rcp(new intOp<double> (valX));
  EXPECT_DOUBLE_EQ(int1->val(), 11);
}

TEST ( Double_Ast_int_Test, test2)
{ 
  // int(x)
  // integer part of the real variable x 
  RCP<astNode<double> > valX = rcp(new numval<double> (-11.2423));   // test value
  RCP<astNode<double> > int1 = rcp(new intOp<double> (valX));
  EXPECT_DOUBLE_EQ(int1->val(),-11);
}

TEST ( Double_Ast_stp_Test, test1)
{
  RCP<astNode<double> > valA = rcp(new numval<double> (10.25));
  RCP<astNode<double> > stp1 = rcp(new stpOp<double> (valA));
  EXPECT_DOUBLE_EQ(stp1->val(), 1);
}

TEST ( Double_Ast_stp_Test, test2)
{
  RCP<astNode<double> > valA = rcp(new numval<double> (-2));
  RCP<astNode<double> > stp1 = rcp(new stpOp<double> (valA));
  EXPECT_DOUBLE_EQ(stp1->val(), 0);
}

int main (int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}


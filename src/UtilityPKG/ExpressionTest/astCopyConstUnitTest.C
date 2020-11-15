//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
#include <N_UTL_BreakPoint.h>
#include <N_UTL_ExtendedString.h>

using namespace Teuchos;

//ASSERT_TRUE(1 == 1);

//-------------------------------------------------------------------------------
// test values of binary operators 
//
#define AST_BINARY_OP_TEST_MACRO(TYPE,NAME,OP,VAL1, VAL2) \
TEST ( NAME, OP ) \
{  \
  RCP<astNode<TYPE> > arg1 = rcp(new numval<TYPE> (VAL1)); \
  RCP<astNode<TYPE> > arg2 = rcp(new numval<TYPE> (VAL2)); \
  RCP<astNode<TYPE> > arg1copy = (arg1); \
  RCP<astNode<TYPE> > arg2copy = (arg2); \
  EXPECT_EQ(arg1->val(),arg1copy->val()); \
  EXPECT_EQ(arg2->val(),arg2copy->val()); \
  RCP<OP<TYPE> > OP_1 = rcp(new OP<TYPE>(arg1,arg2));  \
  RCP<OP<TYPE> > OP_2 = OP_1;  \
  EXPECT_EQ(OP_1->val(),OP_2->val()); \
  RCP<astNode<TYPE> > OP_3 = OP_2; \
  EXPECT_EQ(OP_1->val(),OP_3->val()); \
}

AST_BINARY_OP_TEST_MACRO(double,Double_Binary_Ast_CloneOp,binaryAddOp,1.0,2.0) 
AST_BINARY_OP_TEST_MACRO(double,Double_Binary_Ast_CloneOp,binaryMinusOp,1.0,2.0) 
AST_BINARY_OP_TEST_MACRO(double,Double_Binary_Ast_CloneOp,binaryMulOp,1.0,2.0) 
AST_BINARY_OP_TEST_MACRO(double,Double_Binary_Ast_CloneOp,binaryDivOp,1.0,2.0) 
AST_BINARY_OP_TEST_MACRO(double,Double_Binary_Ast_CloneOp,binaryModOp,15,3)  // must be ints

typedef std::complex<double> cmplx;

AST_BINARY_OP_TEST_MACRO(cmplx,Complex_Binary_Ast_CloneOp,binaryAddOp,cmplx(1.0,0.5),cmplx(2.0,1.0)) 
AST_BINARY_OP_TEST_MACRO(cmplx,Complex_Binary_Ast_CloneOp,binaryMinusOp,cmplx(1.0,0.5),cmplx(2.0,1.0)) 
AST_BINARY_OP_TEST_MACRO(cmplx,Complex_Binary_Ast_CloneOp,binaryMulOp,cmplx(1.0,0.5),cmplx(2.0,1.0)) 
AST_BINARY_OP_TEST_MACRO(cmplx,Complex_Binary_Ast_CloneOp,binaryDivOp,cmplx(1.0,0.5),cmplx(2.0,1.0)) 

//-------------------------------------------------------------------------------
// test values of unary std functions 
//
#define AST_OP_TEST_MACRO(TYPE,NAME,OP,VAL) \
TEST ( NAME, OP )  \
{  \
  RCP<astNode<TYPE> > arg1 = rcp(new numval<TYPE> (VAL)); \
  RCP<OP<TYPE> > OP_1 = rcp(new OP<TYPE>(arg1));  \
  RCP<OP<TYPE> > OP_2 = OP_1;  \
  EXPECT_EQ(OP_1->val(), OP_2->val());  \
  RCP<astNode<TYPE> > OP_3 = OP_2; \
  EXPECT_EQ(OP_1->val(),OP_3->val()); \
}

AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_CloneOps, sqrtOp,  4.0)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_CloneOps, expOp, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_CloneOps, absOp, -0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_CloneOps, sinOp, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_CloneOps, cosOp, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_CloneOps, acosOp, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_CloneOps, acoshOp, 1.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_CloneOps, asinOp, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_CloneOps, asinhOp, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_CloneOps, atanOp, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_CloneOps, atanhOp, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_CloneOps, coshOp, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_CloneOps, logOp, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_CloneOps, log10Op, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_CloneOps, sinhOp, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_CloneOps, tanOp, 0.5)
AST_OP_TEST_MACRO(double, Double_UnaryFunc_Ast_CloneOps, tanhOp, 0.5)

AST_OP_TEST_MACRO(cmplx, Complex_UnaryFunc_Ast_CloneOps, sqrtOp,  cmplx(4.0, 0.2))
AST_OP_TEST_MACRO(cmplx, Complex_UnaryFunc_Ast_CloneOps, expOp, cmplx(.5, 0.2))
AST_OP_TEST_MACRO(cmplx, Complex_UnaryFunc_Ast_CloneOps, absOp, cmplx(-.5, 0.2))
AST_OP_TEST_MACRO(cmplx, Complex_UnaryFunc_Ast_CloneOps, sinOp, cmplx(.5, 0.2))
AST_OP_TEST_MACRO(cmplx, Complex_UnaryFunc_Ast_CloneOps, cosOp, cmplx(.5, 0.2))
AST_OP_TEST_MACRO(cmplx, Complex_UnaryFunc_Ast_CloneOps, acosOp, cmplx(.5, 0.2))
AST_OP_TEST_MACRO(cmplx, Complex_UnaryFunc_Ast_CloneOps, acoshOp, cmplx(1.5, 0.2))
AST_OP_TEST_MACRO(cmplx, Complex_UnaryFunc_Ast_CloneOps, asinOp, cmplx(.5, 0.2))
AST_OP_TEST_MACRO(cmplx, Complex_UnaryFunc_Ast_CloneOps, asinhOp, cmplx(.5, 0.2))
AST_OP_TEST_MACRO(cmplx, Complex_UnaryFunc_Ast_CloneOps, atanOp, cmplx(.5, 0.2))
AST_OP_TEST_MACRO(cmplx, Complex_UnaryFunc_Ast_CloneOps, atanhOp, cmplx(.5, 0.2))
AST_OP_TEST_MACRO(cmplx, Complex_UnaryFunc_Ast_CloneOps, coshOp, cmplx(.5, 0.2))
AST_OP_TEST_MACRO(cmplx, Complex_UnaryFunc_Ast_CloneOps, logOp, cmplx(.5, 0.2))
AST_OP_TEST_MACRO(cmplx, Complex_UnaryFunc_Ast_CloneOps, log10Op, cmplx(.5, 0.2))
AST_OP_TEST_MACRO(cmplx, Complex_UnaryFunc_Ast_CloneOps, sinhOp, cmplx(.5, 0.2))
AST_OP_TEST_MACRO(cmplx, Complex_UnaryFunc_Ast_CloneOps, tanOp, cmplx(.5, 0.2))
AST_OP_TEST_MACRO(cmplx, Complex_UnaryFunc_Ast_CloneOps, tanhOp, cmplx(.5, 0.2))

#if 0
//-------------------------------------------------------------------------------
// constants
TEST ( Ast_Const_CloneTest, Pi )
{  
  piConstOp<double> testPi1;
  piConstOp<double> testPi2(testPi1);
  EXPECT_EQ(testPi1.val(), testPi2.val());  
  astNode<double> * testPi3 = &testPi2;
  astNode<double> * testPi4 = testPi3->clone();
  EXPECT_EQ(testPi1.val(),testPi4->val());
  delete testPi4;
}

//-------------------------------------------------------------------------------
// other functions
TEST ( Double_Ast_Func_CloneTest, powOp )
{  
  numval <double> arg1(2.0), arg2(3.0); 
  powOp<double> testPow1(&arg1,&arg2);
  powOp<double> testPow2(testPow1);
  EXPECT_EQ(testPow1.val(), testPow2.val());

  astNode<double> * testPow3 = &testPow2;
  astNode<double> * testPow4 = testPow3->clone();
  EXPECT_EQ(testPow1.val(),testPow4->val());
  delete testPow4;
}

TEST ( Double_Ast_Func_CloneTest, phaseOp )
{  
  double a1(-1.0);
  numval <double> arg1(a1);
  phaseOp<double> testPhase1(&arg1);
  phaseOp<double> testPhase2(testPhase1);
  EXPECT_EQ(testPhase1.val(), testPhase2.val());

  astNode<double> * testPhase3 = &testPhase2;
  astNode<double> * testPhase4 = testPhase3->clone();
  EXPECT_EQ(testPhase1.val(),testPhase4->val());
  delete testPhase4;
}

TEST ( Complex_Ast_Func_CloneTest, powOp )
{  
  std::complex<double> a1(2.0,3.0);
  std::complex<double> a2(3.0,4.0);
  numval <std::complex<double> > arg1(a1), arg2(a2);
  powOp<std::complex<double> > testPow1(&arg1,&arg2);
  powOp<std::complex<double> > testPow2(testPow1);
  EXPECT_EQ(testPow1.val(), testPow2.val());

  astNode<std::complex<double> > * testPow3 = &testPow2;
  astNode<std::complex<double> > * testPow4 = testPow3->clone();
  EXPECT_EQ(testPow1.val(),testPow4->val());
  delete testPow4;
}

TEST ( Complex_Ast_Func_CloneTest, phaseOp )
{  
  std::complex<double> a1(-1.0,0.0);
  numval <std::complex<double> > arg1(a1);
  phaseOp<std::complex<double> > testPhase1(&arg1);
  phaseOp<std::complex<double> > testPhase2(testPhase1);
  EXPECT_EQ(testPhase1.val(), testPhase2.val());

  astNode<std::complex<double> > * testPhase3 = &testPhase2;
  astNode<std::complex<double> > * testPhase4 = testPhase3->clone();
  EXPECT_EQ(testPhase1.val(),testPhase4->val());
  delete testPhase4;
}

TEST ( Complex_Ast_Func_CloneTest, realOp )
{  
  std::complex<double> a1(-1.0,0.3);
  numval <std::complex<double> > arg1(a1);
  realOp<std::complex<double> > testReal1(&arg1);
  realOp<std::complex<double> > testReal2(testReal1);
  EXPECT_EQ(testReal1.val(), testReal2.val());

  astNode<std::complex<double> > * testReal3 = &testReal2;
  astNode<std::complex<double> > * testReal4 = testReal3->clone();
  EXPECT_EQ(testReal1.val(),testReal4->val());
  delete testReal4;
}

TEST ( Complex_Ast_Func_CloneTest, imagOp )
{  
  std::complex<double> a1(-1.0,0.3);
  numval <std::complex<double> > arg1(a1);
  imagOp<std::complex<double> > testImag1(&arg1);
  imagOp<std::complex<double> > testImag2(testImag1);
  EXPECT_EQ(testImag1.val(), testImag2.val());

  astNode<std::complex<double> > * testImag3 = &testImag2;
  astNode<std::complex<double> > * testImag4 = testImag3->clone();
  EXPECT_EQ(testImag1.val(),testImag4->val());
  delete testImag4;
}


// the copy ctor and clone tests of the spice independent source functions will pass
// if the spice source node classes are modified so that the "time" operator is
// copied as a raw pointer rather than cloned.
//
// I haven't yet decided how best to handled "specials".  Like .params and .funcs,
// they depend on stuff outside the AST tree.  However, unlike .params and .funcs, 
// they do not depend on external expressions.  Each of them just has a special
// singleton node that is owned by the newExpression class.
#if 0
//-------------------------------------------------------------------------------
// spice time-dependent source functions

TEST ( Double_Ast_Spice_Src_CloneTest, spicePulseOp )
{  
  // these numbers are from the VPULSE regression test:
  //0V 1V 0S 10US 10US 0.1US 20.1US
  double v1=0.0, v2=1.0, td=0.0, tr=10.0e-6, tf=10.0e-6, pw=0.1e-6, per=20.1e-6, time = 0.0;

  numval<double> v1_op(v1), v2_op(v2), td_op(td), tr_op(tr), tf_op(tf), pw_op(pw), per_op(per);
  specialsOp<double> time_op(std::string("time"));

  // test initial value(v1)
  time_op.setValue(time);
  spicePulseOp<double> pulse( &v1_op, &v2_op, &td_op, &tr_op, &tf_op, &pw_op, &per_op, &time_op);
  spicePulseOp<double> pulseCopy(pulse);
  EXPECT_EQ(pulseCopy.val(), v1);

  // test post-rise value(v2)
  time=tr+0.5*pw;
  time_op.setValue(time);
  EXPECT_EQ(pulseCopy.val(), v2);

  astNode<double> * testPulse3 = &pulseCopy;
  astNode<double> * testPulse4 = testPulse3->clone();
  EXPECT_EQ(pulse.val(),testPulse4->val());
  delete testPulse4;
}

TEST ( Double_Ast_Spice_Src_CloneTest, spicePulseOp_breakPoints )
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

  numval<double> v1_op(v1), v2_op(v2), td_op(td), tr_op(tr), tf_op(tf), pw_op(pw), per_op(per);
  specialsOp<double> time_op(std::string("time"));

  // test initial value(v1)
  time_op.setValue(time);
  spicePulseOp<double> pulse( &v1_op, &v2_op, &td_op, &tr_op, &tf_op, &pw_op, &per_op, &time_op);
  spicePulseOp<double> pulseCopy(pulse);
  astNode<double> * testPulse3 = &pulseCopy;
  astNode<double> * testPulse4 = testPulse3->clone();

  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  bool ret = testPulse4->getBreakPoints(breakPointTimes);
  delete testPulse4;

  {
    double tol=1.0e-20;
    Xyce::Util::BreakPointLess breakPointLess(tol);
    Xyce::Util::BreakPointEqual breakPointEqual(tol);
    std::sort ( breakPointTimes.begin(), breakPointTimes.end(), breakPointLess );
    std::vector<Xyce::Util::BreakPoint>::iterator it = std::unique ( breakPointTimes.begin(), breakPointTimes.end(), breakPointEqual );
    breakPointTimes.resize( std::distance (breakPointTimes.begin(), it ));
  }

  std::vector<double> bpTimes(breakPointTimes.size());
  for (int ii=0;ii<breakPointTimes.size();ii++) { bpTimes[ii] = breakPointTimes[ii].value(); }
  EXPECT_EQ(breakPointTimes.size(), bpTestVec.size());
  EXPECT_EQ(bpTimes,bpTestVec);
}

TEST ( Double_Ast_Spice_Src_CloneTest, spiceSinOp )
{
  // these numbers are from the  VSIN/bug1679 regression test:
  // SIN ( 1.65 1.65 10000 0 0 -90 )
  //
  double v0(1.65), va(1.65), freq(10000), td(0.0), theta(0.0), phase(-90), time(0.0);
  numval<double> v0_op(v0), va_op(va), freq_op(freq), td_op(td), theta_op(theta), phase_op(phase);

  specialsOp<double> time_op(std::string("time"));

  // test the DCOP value , which should be: DCOPValue = V0 + VA * sin (2.0*mpi*(PHASE/360));
  time_op.setValue(time);
  spiceSinOp<double> sinOp (&v0_op, &va_op, &freq_op, &td_op, &theta_op, &phase_op, &time_op);
  spiceSinOp<double> sinCopy(sinOp);
  astNode<double> * testSin3 = &sinCopy;
  astNode<double> * testSin4 = testSin3->clone();

  double DCOPValue = v0 + va * std::sin (2.0*M_PI*(phase/360));
  EXPECT_EQ(sinCopy.val(), DCOPValue);
  EXPECT_EQ(testSin4->val(), DCOPValue);

  // test a transient value, which should be 
  // TRANValue = (v0) + (va) * std::sin(2.0*mpi*((freq)*time + (phase)/360)) * std::exp( -(time*(theta)));
  time=1.0/10000.0;
  time_op.setValue(time);
  double TRANValue = v0 + va * std::sin(2.0*M_PI*((freq)*time + (phase)/360)) * std::exp( -(time*(theta)));
  EXPECT_EQ(sinCopy.val(), TRANValue);
  EXPECT_EQ(testSin4->val(), TRANValue);

  delete testSin4;
}

TEST ( Double_Ast_Spice_Src_CloneTest, spiceExpOp )
{
  // these numbers are from the   SOURCES/sources.cir
  //  V1=1.1 V2=2 TD1=2ns TAU1=15ns TD2=5ns TAU2=30ns
  //
  double v1(1.1), v2(2.0), td1(2e-9), tau1(15e-9), td2(5e-9), tau2(30e-9), time(0.0);
  numval<double> v1_op(v1), v2_op(v2), td1_op(td1), tau1_op(tau1), td2_op(td2), tau2_op(tau2);

  specialsOp<double> time_op(std::string("time"));

  // test time <= td1, which should be v1
  time_op.setValue(time);
  spiceExpOp<double> expOp ( &v1_op, &v2_op, &td1_op, &tau1_op, &td2_op, &tau2_op, &time_op);
  spiceExpOp<double> expCopy(expOp);
  astNode<double> * testExp3 = &expCopy;
  astNode<double> * testExp4 = testExp3->clone();

  EXPECT_EQ(expCopy.val(), v1);
  EXPECT_EQ(testExp4->val(), v1);

  // test td1 < time <= td2, which should be 
  // value = v1 + (v2-v1)*(1.0-std::exp(-(time-td1)/tau1));
  time=(td1+td2)*0.5;
  time_op.setValue(time);
  double value = v1 + (v2-v1)*(1.0-std::exp(-(time-td1)/tau1));
  EXPECT_EQ(expCopy.val(), value);
  EXPECT_EQ(testExp4->val(), value);

  // test time > td2, which should be 
  // value = v1 + (v2-v1)*(1.0-std::exp(-(time-td1)/tau1)) + (v1-v2)*(1.0-std::exp(-(time-td2)/tau2)) ;
  time=td2*1.1;
  time_op.setValue(time);
  value = v1 + (v2-v1)*(1.0-std::exp(-(time-td1)/tau1)) + (v1-v2)*(1.0-std::exp(-(time-td2)/tau2)) ;
  EXPECT_EQ(expCopy.val(), value);
  EXPECT_EQ(testExp4->val(), value);

  delete testExp4;
}

TEST ( Double_Ast_Spice_Src_CloneTest, spiceSffmOp )
{
  // these numbers are from the   SOURCES/sources.cir
  //  V0=-0.5 VA=2 FC=100meg MDI=0.3 FS=2.1meg
  //
  double v0(-0.5), va(2.0), fc(100e6), mdi(0.3), fs(2.1e6), time(0.0);
  numval<double> v0_op(v0), va_op(va), fc_op(fc), mdi_op(mdi), fs_op(fs);
  specialsOp<double> time_op(std::string("time"));
  spiceSffmOp<double> sffmOp ( &v0_op, &va_op, &fc_op, &mdi_op, &fs_op, &time_op );
  spiceSffmOp<double> sffmCopy (sffmOp);
  astNode<double> * testSffm3 = &sffmCopy;
  astNode<double> * testSffm4 = testSffm3->clone();

  // test, which should be    
  // value = v0 + va * sin((2 * mpi * fc * time) + mdi * sin (2 * mpi * fs * time));
  time=0.1;
  time_op.setValue(time);
  double value = v0 + va * sin((2 * M_PI * fc * time) + mdi * sin (2 * M_PI * fs * time));
  EXPECT_EQ(sffmCopy.val(), value);
  EXPECT_EQ(testSffm4->val(), value);

  delete testSffm4;
}
#endif

//-------------------------------------------------------------------------------
// derivatives of binary operators
TEST ( Double_Ast_Deriv_CloneTest, test1)
{
  numval <double> val1(2.0), val2(5.0);
  paramOp <double> arg1(std::string("A")), arg2(std::string("B"));
  arg1.setNode(&val1);
  arg2.setNode(&val2);
  binaryAddOp<double> binaryAddOp_1(&arg1,&arg2);

  binaryAddOp<double> binaryAddCopy_1(binaryAddOp_1);
  astNode<double> * testAdd3 = &binaryAddCopy_1;
  astNode<double> * testAdd4 = testAdd3->clone();

  // re-resolve the params on the copy and test
  std::vector<paramOp<double> *> paramVec; 
  binaryAddCopy_1.getParamOps (paramVec);
  for (int ii=0;ii<paramVec.size();ii++)
  {
    std::string name = paramVec[ii]->getName(); Xyce::Util::toLower(name);
    if (name == std::string("a")) { paramVec[ii]->setNode(&val1); }
    else if (name == std::string("b")) { paramVec[ii]->setNode(&val2); }
  }

  // value
  EXPECT_EQ(binaryAddCopy_1.val(),(2.0 + 5.0));

  // derivs
  paramVec[0]->setDerivIndex(0);
  paramVec[1]->setDerivIndex(1);
  EXPECT_EQ(binaryAddCopy_1.dx(0),1.0);
  EXPECT_EQ(binaryAddCopy_1.dx(1),1.0);


  // re-resolve the params on the clone and test
  paramVec.clear();
  testAdd4->getParamOps (paramVec);
  for (int ii=0;ii<paramVec.size();ii++)
  {
    std::string name = paramVec[ii]->getName(); Xyce::Util::toLower(name);
    if (name == std::string("a")) { paramVec[ii]->setNode(&val1); }
    else if (name == std::string("b")) { paramVec[ii]->setNode(&val2); }
  }

  // value
  EXPECT_EQ(testAdd4->val(),(2.0 + 5.0));

  // derivs
  paramVec[0]->setDerivIndex(0);
  paramVec[1]->setDerivIndex(1);
  EXPECT_EQ(testAdd4->dx(0),1.0);
  EXPECT_EQ(testAdd4->dx(1),1.0);

  delete testAdd4;
}

// macro tests both copy construction and clone
// As these tests rely on .params, they have to be re-resolved.
#define AST_BINARY_DERIV_TEST_MACRO(TYPE,NAME,OP,CPPOP, VAL1, VAL2, D1, D2) \
TEST ( NAME, OP ) \
{ \
  numval <TYPE> val1(VAL1), val2(VAL2); \
  paramOp <TYPE> arg1(std::string("A"),&val1), arg2(std::string("B"),&val2); \
  OP<TYPE> OP_1(&arg1,&arg2); \
  OP<TYPE> OP_copy(OP_1); \
  std::vector<paramOp<TYPE> *> paramVec;  \
  OP_copy.getParamOps (paramVec); \
  for (int ii=0;ii<paramVec.size();ii++) \
  { \
    std::string name = paramVec[ii]->getName(); Xyce::Util::toLower(name); \
    if (name == std::string("a")) { paramVec[ii]->setNode(&val1); } \
    else if (name == std::string("b")) { paramVec[ii]->setNode(&val2); } \
  } \
  EXPECT_EQ(OP_copy.val(),(VAL1 CPPOP VAL2)); \
  paramVec[0]->setDerivIndex(0); \
  paramVec[1]->setDerivIndex(1); \
  EXPECT_EQ(OP_copy.dx(0),D1); \
  EXPECT_EQ(OP_copy.dx(1),D2); \
\
  astNode<TYPE> * testOP3 = &OP_copy; \
  astNode<TYPE> * testOP4 = testOP3->clone(); \
  paramVec.clear(); \
  testOP4->getParamOps (paramVec); \
  for (int ii=0;ii<paramVec.size();ii++) \
  { \
    std::string name = paramVec[ii]->getName(); Xyce::Util::toLower(name); \
    if (name == std::string("a")) { paramVec[ii]->setNode(&val1); } \
    else if (name == std::string("b")) { paramVec[ii]->setNode(&val2); } \
  } \
  EXPECT_EQ(testOP4->val(),(VAL1 CPPOP VAL2)); \
  paramVec[0]->setDerivIndex(0); \
  paramVec[1]->setDerivIndex(1); \
  EXPECT_EQ(testOP4->dx(0),D1); \
  EXPECT_EQ(testOP4->dx(1),D2); \
  delete testOP4;\
}

AST_BINARY_DERIV_TEST_MACRO(double,Double_Ast_Deriv_CloneTest,binaryAddOp,+,3.0,4.0,1.0,1.0)
AST_BINARY_DERIV_TEST_MACRO(double,Double_Ast_Deriv_CloneTest,binaryMinusOp,-,3.0,4.0,1.0,-1.0)
AST_BINARY_DERIV_TEST_MACRO(double,Double_Ast_Deriv_CloneTest,binaryMulOp,*,3.0,4.0,4.0,3.0)
AST_BINARY_DERIV_TEST_MACRO(double,Double_Ast_Deriv_CloneTest,binaryDivOp,/,3.0,4.0,0.25,(-3.0/16.0))

AST_BINARY_DERIV_TEST_MACRO(cmplx,Complex_Ast_Deriv_CloneTest,binaryAddOp,+,cmplx(3.0,1.0),cmplx(4.0,2.0),cmplx(1.0,0.0),cmplx(1.0,0.0))
AST_BINARY_DERIV_TEST_MACRO(cmplx,Complex_Ast_Deriv_CloneTest,binaryMinusOp,-,cmplx(3.0,1.0),cmplx(4.0,2.0),cmplx(1.0,0.0),cmplx(-1.0,0.0))
AST_BINARY_DERIV_TEST_MACRO(cmplx,Complex_Ast_Deriv_CloneTest,binaryMulOp,*,cmplx(3.0,1.0),cmplx(4.0,2.0),cmplx(4.0,2.0),cmplx(3.0,1.0))
AST_BINARY_DERIV_TEST_MACRO(cmplx,Complex_Ast_Deriv_CloneTest,binaryDivOp,/,cmplx(3.0,1.0),cmplx(4.0,2.0), +cmplx(4.0,2.0)/(cmplx(4.0,2.0)*cmplx(4.0,2.0)), -cmplx(3.0,1.0)/(cmplx(4.0,2.0)*cmplx(4.0,2.0)))

//-------------------------------------------------------------------------------
// derivatives of unary std functions
// macro tests both copy construction and clone
// As these tests rely on .params, they have to be re-resolved.
#define AST_UNARY_DERIV_TEST_MACRO(TYPE,NAME,OP,CPPFUNC,VAL1,D1) \
TEST ( NAME, OP ) \
{ \
  numval <TYPE> val1(VAL1); \
  paramOp <TYPE> arg1(std::string("A"),&val1); \
  OP<TYPE> OP_1(&arg1); \
  \
  OP<TYPE> OP_copy(OP_1); \
  std::vector<paramOp<TYPE> *> paramVec;  \
  OP_copy.getParamOps (paramVec); \
  paramVec[0]->setNode(&val1);  \
  paramVec[0]->setDerivIndex(0); \
  EXPECT_EQ(OP_copy.val(),CPPFUNC(VAL1)); \
  EXPECT_EQ( (OP_copy.dx(0)-D1), 0.0); \
  \
  astNode<TYPE> * testOP3 = &OP_copy; \
  astNode<TYPE> * testOP4 = testOP3->clone(); \
  paramVec.clear(); \
  testOP4->getParamOps (paramVec); \
  paramVec[0]->setNode(&val1);  \
  paramVec[0]->setDerivIndex(0); \
  EXPECT_EQ(testOP4->val(),CPPFUNC(VAL1)); \
  EXPECT_EQ( (testOP4->dx(0)-D1), 0.0); \
  delete testOP4;\
} 

// complex
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_CloneOps, sqrtOp, std::sqrt, 4.0, (0.5/std::sqrt(4.0)))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_CloneOps, expOp, std::exp, 0.5, std::exp(0.5))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_CloneOps, absOp, std::abs, -0.5, -1.0)
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_CloneOps, sinOp, std::sin, 0.5, std::cos(0.5))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_CloneOps, cosOp, std::cos, 0.5, -std::sin(0.5))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_CloneOps, tanOp, std::tan, 0.5, (1.0+(std::tan(0.5)*std::tan(0.5))))

AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_CloneOps, asinOp, std::asin, 0.5,(+1.0/std::sqrt(1.0-0.5*0.5)))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_CloneOps, acosOp, std::acos, 0.5,(-1.0/std::sqrt(1.0-0.5*0.5)))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_CloneOps, atanOp, std::atan, 0.5,(+1.0/(1+0.5*0.5)))

AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_CloneOps, sinhOp, std::sinh, 0.5, std::cosh(0.5))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_CloneOps, coshOp, std::cosh, 0.5, std::sinh(0.5))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_CloneOps, tanhOp, std::tanh, 0.1, 1.0/(std::cosh(0.1)*std::cosh(0.1)))

AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_CloneOps, asinhOp, std::asinh, 0.5, 1.0/(std::sqrt(1+0.5*0.5)))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_CloneOps, acoshOp, std::acosh, 1.5, 1.0/(std::sqrt( (1.5-1) * (1.5+1) )))
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_CloneOps, atanhOp, std::atanh, 0.5, 1.0/(1.0-0.5*0.5))

AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_CloneOps, logOp, std::log, 0.5, 1.0/0.5)
AST_UNARY_DERIV_TEST_MACRO(double, Double_UnaryDeriv_Ast_CloneOps, log10Op, std::log10, 0.5, 1.0/(0.5*std::log(10)))

// complex
AST_UNARY_DERIV_TEST_MACRO(cmplx, Complex_UnaryDeriv_Ast_CloneOps, sqrtOp, std::sqrt, cmplx(4.0,0.2), (0.5/std::sqrt(cmplx(4.0,0.2))))
AST_UNARY_DERIV_TEST_MACRO(cmplx, Complex_UnaryDeriv_Ast_CloneOps, expOp, std::exp, cmplx(0.5,0.2), std::exp(cmplx(0.5,0.2)))
AST_UNARY_DERIV_TEST_MACRO(cmplx, Complex_UnaryDeriv_Ast_CloneOps, absOp, std::abs, cmplx(-0.5,0.2), -1.0)
AST_UNARY_DERIV_TEST_MACRO(cmplx, Complex_UnaryDeriv_Ast_CloneOps, sinOp, std::sin, cmplx(0.5,0.2), std::cos(cmplx(0.5,0.2)))
AST_UNARY_DERIV_TEST_MACRO(cmplx, Complex_UnaryDeriv_Ast_CloneOps, cosOp, std::cos, cmplx(0.5,0.2), -std::sin(cmplx(0.5,0.2)))
AST_UNARY_DERIV_TEST_MACRO(cmplx, Complex_UnaryDeriv_Ast_CloneOps, tanOp, std::tan, cmplx(0.5,0.2), (1.0+std::tan(cmplx(0.5,0.2))*std::tan(cmplx(0.5,0.2))))

AST_UNARY_DERIV_TEST_MACRO(cmplx, Complex_UnaryDeriv_Ast_CloneOps, asinOp, std::asin, cmplx(0.5,0.2),(+1.0/std::sqrt(1.0-cmplx(0.5,0.2)*cmplx(0.5,0.2))))
AST_UNARY_DERIV_TEST_MACRO(cmplx, Complex_UnaryDeriv_Ast_CloneOps, acosOp, std::acos, cmplx(0.5,0.2),(-1.0/std::sqrt(1.0-cmplx(0.5,0.2)*cmplx(0.5,0.2))))
AST_UNARY_DERIV_TEST_MACRO(cmplx, Complex_UnaryDeriv_Ast_CloneOps, atanOp, std::atan, cmplx(0.5,0.2),(+1.0/(1.0+cmplx(0.5,0.2)*cmplx(0.5,0.2))))

AST_UNARY_DERIV_TEST_MACRO(cmplx, Complex_UnaryDeriv_Ast_CloneOps, sinhOp, std::sinh, cmplx(0.5,0.2), std::cosh(cmplx(0.5,0.2)))
AST_UNARY_DERIV_TEST_MACRO(cmplx, Complex_UnaryDeriv_Ast_CloneOps, coshOp, std::cosh, cmplx(0.5,0.2), std::sinh(cmplx(0.5,0.2)))
AST_UNARY_DERIV_TEST_MACRO(cmplx, Complex_UnaryDeriv_Ast_CloneOps, tanhOp, std::tanh, cmplx(0.1,0.2), 1.0/(std::cosh(cmplx(0.1,0.2))*std::cosh(cmplx(0.1,0.2))))

AST_UNARY_DERIV_TEST_MACRO(cmplx, Complex_UnaryDeriv_Ast_CloneOps, asinhOp, std::asinh, cmplx(0.5,0.2), 1.0/(std::sqrt(1.0+cmplx(0.5,0.2)*cmplx(0.5,0.2))))
AST_UNARY_DERIV_TEST_MACRO(cmplx, Complex_UnaryDeriv_Ast_CloneOps, acoshOp, std::acosh, cmplx(1.5,0.2), (1.0/(std::sqrt( (cmplx(1.5,0.2)-1.0) * (cmplx(1.5,0.2)+1.0) ))))
AST_UNARY_DERIV_TEST_MACRO(cmplx, Complex_UnaryDeriv_Ast_CloneOps, atanhOp, std::atanh, cmplx(0.5,0.2), 1.0/(1.0-cmplx(0.5,0.2)*cmplx(0.5,0.2)))

AST_UNARY_DERIV_TEST_MACRO(cmplx, Complex_UnaryDeriv_Ast_CloneOps, logOp, std::log, cmplx(0.5,0.2), 1.0/cmplx(0.5,0.2))
AST_UNARY_DERIV_TEST_MACRO(cmplx, Complex_UnaryDeriv_Ast_CloneOps, log10Op, std::log10, cmplx(0.5,0.2), 1.0/(cmplx(0.5,0.2)*std::log(10)))

//-------------------------------------------------------------------------------
// pow function, test of clone and cctor
TEST ( Double_Ast_Deriv_Test, powOp )
{  
  double A=7.0;
  double B=4.0;
  numval <double> val1(A), val2(B); 
  paramOp <double> arg1(std::string("A"),&val1), arg2(std::string("B"),&val2); 

  // test copy ctor
  // value
  powOp<double> testPow(&arg1,&arg2);
  powOp<double> testPowCopy(testPow);
  std::vector<paramOp<double> *> paramVec;
  testPowCopy.getParamOps (paramVec);
  paramVec[0]->setNode(&val1);
  paramVec[1]->setNode(&val2);
  EXPECT_EQ(testPowCopy.val(), std::pow(A,B));  

  // derivative
  paramVec[0]->setDerivIndex(0);
  paramVec[1]->setDerivIndex(1);
  EXPECT_EQ(testPowCopy.dx(0)-((B/A)*std::pow(A,(B))), 0.0 );
  EXPECT_EQ(testPowCopy.dx(1)-(std::log(A))*std::pow(A,B), 0.0 );

  // test clone
  astNode<double> * testOP3 = &testPowCopy;
  astNode<double> * testOP4 = testOP3->clone();
  paramVec.clear();
  testOP4->getParamOps (paramVec);
  paramVec[0]->setNode(&val1); 
  paramVec[1]->setNode(&val2); 
  paramVec[0]->setDerivIndex(0);
  paramVec[1]->setDerivIndex(1);
  EXPECT_EQ(testOP4->val(), std::pow(A,B));  
  EXPECT_EQ(testOP4->dx(0)-((B/A)*std::pow(A,(B))), 0.0 );
  EXPECT_EQ(testOP4->dx(1)-(std::log(A))*std::pow(A,B), 0.0 );
  delete testOP4;
}

TEST ( Complex_Ast_Deriv_Test, powOp )
{
  cmplx A=cmplx(7.0,2.0);
  cmplx B=cmplx(4.0,3.0);
  numval <cmplx> val1(A), val2(B); 
  paramOp <cmplx> arg1(std::string("A"),&val1), arg2(std::string("B"),&val2); 

  // copy ctor tests
  // value
  powOp<cmplx> testPow(&arg1,&arg2);

  powOp<cmplx> testPowCopy(testPow);
  std::vector<paramOp<cmplx> *> paramVec;
  testPowCopy.getParamOps (paramVec);
  paramVec[0]->setNode(&val1);
  paramVec[1]->setNode(&val2);

  EXPECT_EQ(testPowCopy.val(), std::pow(A,B));  

  // general:
  // (B.dx(i)*std::log(A.val())+B.val()*A.dx(i)/A.val())*std::pow(A.val(),B.val())
 
  // w.r.t. A:
  // (+B/A)*std::pow(A,B)
  
  // w.r.t. B:
  // (std::log(A))*std::pow(A,B)
  
  // derivative
  paramVec[0]->setDerivIndex(0);
  paramVec[1]->setDerivIndex(1);
  EXPECT_EQ(testPowCopy.dx(0)-((B/A)*std::pow(A,(B))), 0.0 );
  EXPECT_EQ(testPowCopy.dx(1)-(std::log(A))*std::pow(A,B), 0.0 );

  // test clone
  astNode<cmplx> * testOP3 = &testPowCopy;
  astNode<cmplx> * testOP4 = testOP3->clone();
  paramVec.clear();
  testOP4->getParamOps (paramVec);
  paramVec[0]->setNode(&val1); 
  paramVec[1]->setNode(&val2); 
  paramVec[0]->setDerivIndex(0);
  paramVec[1]->setDerivIndex(1);
  EXPECT_EQ(testOP4->val(), std::pow(A,B));  
  EXPECT_EQ(testOP4->dx(0)-((B/A)*std::pow(A,(B))), 0.0 );
  EXPECT_EQ(testOP4->dx(1)-(std::log(A))*std::pow(A,B), 0.0 );
  delete testOP4;
}

#if 0
//-------------------------------------------------------------------------------
TEST ( Double_Ast_Param_Test, paramOp )
{
  double A=7.0;
  double B=4.0;
  numval <double> valA(A), valB(B); 
  paramOp <double> paramA(std::string("A"),&valA), paramB(std::string("B"),&valB); 

  EXPECT_EQ(paramA.val()-A, 0.0);
  EXPECT_EQ(paramB.val()-B, 0.0);
}

//-------------------------------------------------------------------------------
TEST ( Double_Ast_Param_Test, voltageOp1 )
{
  double A=7.0;
  double B=4.0;
  numval <double> valA(A), valB(B); 

  voltageOp <double> voltageA( std::vector<std::string> (1,std::string("A")) );
  voltageOp <double> voltageB( std::vector<std::string> (1,std::string("B")) );

  voltageA.setVals(std::vector<double>(1,A));
  voltageB.setVals(std::vector<double>(1,B));

  EXPECT_EQ(voltageA.val()-A, 0.0);
  EXPECT_EQ(voltageB.val()-B, 0.0);
}

TEST ( Double_Ast_Param_Test, voltageOp2 )
{
  double A=7.0, B=4.0, C=3.0, D=2.0;
  numval <double> valA(A), valB(B), valC(C), valD(D); 

  std::vector<std::string> names1 = { std::string("A"), std::string("B") };
  std::vector<std::string> names2 = { std::string("C"), std::string("D") };
  voltageOp <double> voltage1( names1 );
  voltageOp <double> voltage2( names2 );

  std::vector<double>vals1 = {A,B};
  std::vector<double>vals2 = {C,D};
  voltage1.setVals(vals1);
  voltage2.setVals(vals2);

  EXPECT_EQ(voltage1.val()-(A-B), 0.0);
  EXPECT_EQ(voltage2.val()-(C-D), 0.0);
}

//-------------------------------------------------------------------------------
// calling a .func
TEST ( Double_Ast_Func_Test, test1)
{

  // define the .func
  // This is what should result from a parser specification of .func function1(arg1,arg2) {arg1+arg2}
  std::vector<paramOp<double> *> dummyArgs;
  paramOp <double> dummyArg1(std::string("arg1")), dummyArg2(std::string("arg2")); 
  dummyArgs.push_back(&dummyArg1);
  dummyArgs.push_back(&dummyArg2);
  binaryAddOp <double> function1(&dummyArg1,&dummyArg2); 
  std::string name = "function1";

  // define the expression that calls the .func
  // This is what should result from a parser specification of  .global_param function2 = {function1(5.0,6.0)}
  // Which should evaluate to function2 = 5.0+6.0=11.0
  std::vector<astNode<double> * > args;
  numval <double> arg1(5.0), arg2(6.0); 
  args.push_back(&arg1);
  args.push_back(&arg2);
  funcOp<double> function2(name,&args);

  function2.setNode(&function1);
  function2.setFuncArgs(dummyArgs);
  EXPECT_EQ(function2.val(), 11.0);
}

TEST ( Complex_Ast_Func_Test, test1)
{

  // define the .func
  // This is what should result from a parser specification of .func function1(arg1,arg2) {arg1+arg2}
  std::vector<paramOp<std::complex<double> > *> dummyArgs;
  paramOp <std::complex<double> > dummyArg1(std::string("arg1")), dummyArg2(std::string("arg2")); 
  dummyArgs.push_back(&dummyArg1);
  dummyArgs.push_back(&dummyArg2);
  binaryAddOp <std::complex<double> > function1(&dummyArg1,&dummyArg2); 
  std::string name = "function1";

  // define the expression that calls the .func
  // This is what should result from a parser specification of  .global_param function2 = {function1((5.0,2.0),(6.0,7.0))}
  // Which should evaluate to function2 = (5,2)+(6,7)=(11,9)
  std::vector<astNode<std::complex<double> > * > args;
  numval <std::complex<double> > arg1(std::complex<double> (5.0,2.0)), arg2(std::complex<double> (6.0,7.0));
  args.push_back(&arg1);
  args.push_back(&arg2);
  funcOp<std::complex<double> > function2(name,&args);

  function2.setNode(&function1);
  function2.setFuncArgs(dummyArgs);
  EXPECT_EQ(function2.val(), std::complex<double> (11.0,9.0));
}

TEST ( Double_Ast_Func_Test, test2)
{
  // define the .func
  // This is what should result from a parser specification of .func f1(arg1,arg2) {arg1+arg2}
  paramOp <double> dummyArg1(std::string("arg1")), dummyArg2(std::string("arg2")); 
  std::vector<paramOp<double> *> dummyArgs = {&dummyArg1,&dummyArg2};
  binaryAddOp <double> f1(&dummyArg1,&dummyArg2); 
  std::string name = "f1";

  // define a new func that calls the first .func
  // This is what should result from a parser specification of  .func f2(arg1,arg2) = {arg1+arg2*f1(5.0,6.0)}
  // Which should evaluate to f2 = 5.0+6.0=11.0
  numval <double> arg1(5.0), arg2(6.0); 
  std::vector<astNode<double> * > args = {&arg1,&arg2};
  funcOp<double> f1call(name,&args);
  f1call.setNode(&f1);
  f1call.setFuncArgs(dummyArgs);

  paramOp <double> d1(std::string("arg1"));
  paramOp <double> d2(std::string("arg2"));
  std::vector<paramOp<double> *> dummyArgs2 = {&d1,&d2};
  binaryMulOp <double> prod(&d2,&f1call); 
  binaryAddOp <double> f2(&d1,&prod); 
  std::string name2 = "f2";

  // define expression that calls the second func
  // This is what should result from a parser specification of  .global_param final = {f2(3.0,2.0)}
  // Which should evaluate to function2 = 3.0+2.0*f1(5.0,6.0) = 3.0+2.0*11.0 = 25.0

  numval <double> argFinal1(3.0), argFinal2(2.0); 
  std::vector<astNode<double> * > argsFinal = { &argFinal1, &argFinal2 };
  funcOp<double> function2(name,&argsFinal);
  function2.setNode(&f2);
  function2.setFuncArgs(dummyArgs2);
  EXPECT_EQ(function2.val(), 25.0);
}

TEST ( Complex_Ast_Func_Test, test2)
{
  // define the .func
  // This is what should result from a parser specification of .func f1(arg1,arg2) {arg1+arg2}
  paramOp <std::complex<double> > dummyArg1(std::string("arg1")), dummyArg2(std::string("arg2")); 
  std::vector<paramOp<std::complex<double> > *> dummyArgs = {&dummyArg1,&dummyArg2};
  binaryAddOp <std::complex<double> > f1(&dummyArg1,&dummyArg2); 
  std::string name = "f1";

  // define a new func that calls the first .func
  // This is what should result from a parser specification of  .func f2(arg1,arg2) = {arg1+arg2*f1((5,2),(6,7))}
  // Which should evaluate to f2 = (5,2)+(6,7)=(11,9)
  numval <std::complex<double> > arg1(std::complex<double> (5.0,2.0)), arg2(std::complex<double> (6.0,7.0)); 
  std::vector<astNode<std::complex<double> > * > args = {&arg1,&arg2};
  funcOp<std::complex<double> > f1call(name,&args);
  f1call.setNode(&f1);
  f1call.setFuncArgs(dummyArgs);

  paramOp <std::complex<double> > d1(std::string("arg1"));
  paramOp <std::complex<double> > d2(std::string("arg2"));
  std::vector<paramOp<std::complex<double> > *> dummyArgs2 = {&d1,&d2};
  binaryMulOp <std::complex<double> > prod(&d2,&f1call); 
  binaryAddOp <std::complex<double> > f2(&d1,&prod); 
  std::string name2 = "f2";

  // define expression that calls the second func
  // This is what should result from a parser specification of  .global_param final = {f2((3,1),(2,4))}
  // Which should evaluate to function2 = (3,1)+(2,4)*f1((5,2),(6,7)) = (3,1)+(2,4)*(11,9) 

  numval <std::complex<double> > argFinal1(std::complex<double> (3.0,1.0)), argFinal2(std::complex<double> (2.0,4.0)); 
  std::vector<astNode<std::complex<double> > * > argsFinal = { &argFinal1, &argFinal2 };
  funcOp<std::complex<double> > function2(name,&argsFinal);
  function2.setNode(&f2);
  function2.setFuncArgs(dummyArgs2);
  EXPECT_EQ(function2.val(), std::complex<double>(3,1)+std::complex<double>(2,4)*std::complex<double>(11,9));
}

//-------------------------------------------------------------------------------
TEST ( Double_Ast_Param_Test, test1)
{
  // define the .param
  double A=3.0;
  numval <double> paramValA(A);
  paramOp <double> testParamA(std::string("A"));
  testParamA.setNode(&paramValA);

  // define the expression that uses the .param
  double B=5.0;
  numval <double> paramValB(B);
  binaryAddOp <double> finalExp(&paramValB,&testParamA);

  EXPECT_EQ(finalExp.val(), A+B);
}

TEST ( Complex_Ast_Param_Test, test1)
{
  // define the .param
  std::complex<double> A(3.0,2.0);
  numval <std::complex<double> > paramValA(A);
  paramOp <std::complex<double> > testParamA(std::string("A"));
  testParamA.setNode(&paramValA);

  // define the expression that uses the .param
  std::complex<double> B(5.0,4.0);
  numval <std::complex<double> > paramValB(B);
  binaryAddOp <std::complex<double> > finalExp(&paramValB,&testParamA);

  EXPECT_EQ(finalExp.val(), (A+B));
}

TEST ( Double_Ast_Param_Test, test2)
{
  // define the .params
  numval <double> paramValA(3.0), paramValB(4.0);
  paramOp <double> testParamA(std::string("A")), testParamB(std::string("B"));
  testParamA.setNode(&paramValA);
  testParamB.setNode(&paramValB);

  // define the expression that uses the .params
  binaryMulOp <double> finalExp(&testParamA, &testParamB);

  EXPECT_EQ(finalExp.val(), 12.0);
}

TEST ( Complex_Ast_Param_Test, test2)
{
  // define the .params
  std::complex<double> A(3.0,5.0);
  std::complex<double> B(4.0,6.0);
  numval <std::complex<double> > paramValA(A), paramValB(B);
  paramOp <std::complex<double> > testParamA(std::string("A")), testParamB(std::string("B"));
  testParamA.setNode(&paramValA);
  testParamB.setNode(&paramValB);

  // define the expression that uses the .params
  binaryMulOp <std::complex<double> > finalExp(&testParamA, &testParamB);

  EXPECT_EQ(finalExp.val(), A*B);
}

TEST ( Double_Ast_Param_Test, test3)
{
  // define the .params
  double A=3.0;
  double B=4.0;
  numval <double> paramValA(3.0), paramValB(4.0);
  paramOp <double> testParamA(std::string("A")), testParamB(std::string("B"));
  testParamA.setNode(&paramValA);
  testParamB.setNode(&paramValB);

  // define the expression that uses the .params
  binaryMulOp <double> prodExp(&testParamA, &testParamB);

  // now create another .param and assign prodExp to it
  paramOp <double> testParamC(std::string("C"));
  testParamC.setNode(&prodExp);
  double C = A*B;

  // now create another expression that uses all of them
  binaryAddOp <double> numerator(&testParamA, &testParamB);
  binaryDivOp <double> finalExp(&numerator, &testParamC);

  // final result should be = (A+B)/C
  EXPECT_EQ(finalExp.val(), (A+B)/C);
}

TEST ( Complex_Ast_Param_Test, test3)
{
  // define the .params
  std::complex<double> A(3.0,5.0);
  std::complex<double> B(4.0,6.0);
  numval <std::complex<double> > paramValA(A), paramValB(B);
  paramOp <std::complex<double> > testParamA(std::string("A")), testParamB(std::string("B"));
  testParamA.setNode(&paramValA);
  testParamB.setNode(&paramValB);

  // define the expression that uses the .params
  binaryMulOp <std::complex<double> > prodExp(&testParamA, &testParamB);

  // now create another .param and assign prodExp to it
  paramOp <std::complex<double> > testParamC(std::string("C"));
  testParamC.setNode(&prodExp);
  std::complex<double> C = A*B;

  // now create another expression that uses all of them
  binaryAddOp <std::complex<double> > numerator(&testParamA, &testParamB);
  binaryDivOp <std::complex<double> > finalExp(&numerator, &testParamC);

  // final result should be = (A+B)/C
  EXPECT_EQ(finalExp.val(), (A+B)/C);
}
#endif

//-------------------------------------------------------------------------------
// test values of conditional operators and if statements
// testing the copy ctor and clone
#define AST_IF_OP_TEST_MACRO(TYPE,NAME,SUBNAME,OP, C1, C2, VAL1, VAL2, RESULT) \
TEST ( NAME, SUBNAME ) \
{ \
  numval<double> c1(C1), c2(C2); \
  OP<double> cond1(&c1,&c2); \
  numval<double> ifVal1(VAL1), ifVal2(VAL2); \
  ifStatementOp <double> ifStmt(&cond1,&ifVal1,&ifVal2); \
  ifStatementOp <double> ifStmtCopy(ifStmt); \
  EXPECT_EQ(ifStmt.val(), RESULT);\
  astNode<TYPE> * testOP3 = &ifStmtCopy; \
  astNode<TYPE> * testOP4 = testOP3->clone(); \
  EXPECT_EQ(testOP4->val(), RESULT);\
  delete testOP4; \
}

AST_IF_OP_TEST_MACRO(double,Double_Ast_if_Test,eq1,eqOp,3.0,2.0,2.0,3.0,3.0) 
AST_IF_OP_TEST_MACRO(double,Double_Ast_if_Test,eq2,eqOp,3.0,3.0,2.0,3.0,2.0) 

AST_IF_OP_TEST_MACRO(double,Double_Ast_if_Test,ne1,neOp,3.0,2.0,2.0,3.0,2.0) 
AST_IF_OP_TEST_MACRO(double,Double_Ast_if_Test,ne2,neOp,3.0,3.0,2.0,3.0,3.0) 

AST_IF_OP_TEST_MACRO(double,Double_Ast_if_Test,gt1,gtOp,3.0,2.0,2.0,3.0,2.0) 
AST_IF_OP_TEST_MACRO(double,Double_Ast_if_Test,gt2,gtOp,2.0,3.0,2.0,3.0,3.0) 

AST_IF_OP_TEST_MACRO(double,Double_Ast_if_Test,lt1,ltOp,3.0,2.0,2.0,3.0,3.0) 
AST_IF_OP_TEST_MACRO(double,Double_Ast_if_Test,lt2,ltOp,2.0,3.0,2.0,3.0,2.0) 

AST_IF_OP_TEST_MACRO(double,Double_Ast_if_Test,ge1,geOp,3.0,2.0,2.0,3.0,2.0) 
AST_IF_OP_TEST_MACRO(double,Double_Ast_if_Test,ge2,geOp,3.0,3.0,2.0,3.0,2.0) 
AST_IF_OP_TEST_MACRO(double,Double_Ast_if_Test,ge3,geOp,2.0,3.0,2.0,3.0,3.0) 

AST_IF_OP_TEST_MACRO(double,Double_Ast_if_Test,le1,leOp,3.0,2.0,2.0,3.0,3.0) 
AST_IF_OP_TEST_MACRO(double,Double_Ast_if_Test,le2,leOp,3.0,3.0,2.0,3.0,2.0) 
AST_IF_OP_TEST_MACRO(double,Double_Ast_if_Test,le3,leOp,2.0,3.0,2.0,3.0,2.0) 


// testing the copy ctor and clone
#define AST_LOGIC_OP_TEST_MACRO(TYPE,NAME,SUBNAME,OP, C1, C2, RESULT) \
TEST ( NAME, SUBNAME ) \
{ \
  numval<double> c1(C1), c2(C2); \
  OP<double> op(&c1,&c2); \
  OP<double> opCopy(op); \
  EXPECT_EQ(opCopy.val(), RESULT);\
  astNode<TYPE> * testOP3 = &opCopy; \
  astNode<TYPE> * testOP4 = testOP3->clone(); \
  EXPECT_EQ(testOP4->val(), RESULT);\
  delete testOP4; \
}

AST_LOGIC_OP_TEST_MACRO(double,Double_Ast_logical_Test,and1,andOp, 1.0, 1.0, 1.0) 
AST_LOGIC_OP_TEST_MACRO(double,Double_Ast_logical_Test,and2,andOp, 1.0, 0.0, 0.0) 
AST_LOGIC_OP_TEST_MACRO(double,Double_Ast_logical_Test,and3,andOp, 0.0, 1.0, 0.0) 
AST_LOGIC_OP_TEST_MACRO(double,Double_Ast_logical_Test,and4,andOp, 0.0, 0.0, 0.0) 

AST_LOGIC_OP_TEST_MACRO(double,Double_Ast_logical_Test,or1,orOp, 1.0, 1.0, 1.0) 
AST_LOGIC_OP_TEST_MACRO(double,Double_Ast_logical_Test,or2,orOp, 1.0, 0.0, 1.0) 
AST_LOGIC_OP_TEST_MACRO(double,Double_Ast_logical_Test,or3,orOp, 0.0, 1.0, 1.0) 
AST_LOGIC_OP_TEST_MACRO(double,Double_Ast_logical_Test,or4,orOp, 0.0, 0.0, 0.0) 

AST_LOGIC_OP_TEST_MACRO(double,Double_Ast_logical_Test,xor1,xorOp, 1.0, 1.0, 0.0) 
AST_LOGIC_OP_TEST_MACRO(double,Double_Ast_logical_Test,xor2,xorOp, 1.0, 0.0, 1.0) 
AST_LOGIC_OP_TEST_MACRO(double,Double_Ast_logical_Test,xor3,xorOp, 0.0, 1.0, 1.0) 
AST_LOGIC_OP_TEST_MACRO(double,Double_Ast_logical_Test,xor4,xorOp, 0.0, 0.0, 0.0) 

TEST ( Double_Ast_logical_Test, not1 ) 
{ 
  numval<double> c1(1.0);
  unaryNotOp<double> op(&c1); 
  // test copy
  unaryNotOp<double> opCopy(op);
  EXPECT_EQ(opCopy.val(), 0.0);

  // test clone
  astNode<double> * testOP3 = &opCopy;
  astNode<double> * testOP4 = testOP3->clone();
  EXPECT_EQ(testOP4->val(), 0.0);
  delete testOP4;
}

TEST ( Double_Ast_logical_Test, not2 ) 
{ 
  numval<double> c1(0.0);
  unaryNotOp<double> op(&c1); 
  // test copy
  unaryNotOp<double> opCopy(op);
  EXPECT_EQ(opCopy.val(), 1.0);

  // test clone
  astNode<double> * testOP3 = &opCopy;
  astNode<double> * testOP4 = testOP3->clone();
  EXPECT_EQ(testOP4->val(), 1.0);
  delete testOP4;
}
#endif

#if 0
//-------------------------------------------------------------------------------
// table tests
//
// adapted from break.cir
//
// testing the copy ctor and the clone
TEST ( Double_Ast_table_Test, break1)
{
  specialsOp<double> time_op(std::string("time"));
  std::vector<astNode<double> * > args;
  args.push_back(new numval<double>(0));
  args.push_back(new numval<double>(0));
  args.push_back(new numval<double>(0.3)); 
  args.push_back(new numval<double>(0)); 
  args.push_back(new numval<double>(0.301));
  args.push_back(new numval<double>(2));
  args.push_back(new numval<double>(0.302));
  args.push_back(new numval<double>(2));
  args.push_back(new numval<double>(0.6));
  args.push_back(new numval<double>(1));
  args.push_back(new numval<double>(1));
  args.push_back(new numval<double>(1));

  tableOp<double> table(&time_op, &args);

  tableOp<double> tableCopy(table);

  std::vector<double> times = { 0, 0.3, 0.301, 0.302, 0.6, 1 };
  std::vector<double> refRes = { 0, 0, 2, 2, 1, 1 };
  int numpoints = times.size();
  std::vector<double> result(numpoints,0.0);
  for (int ii=0;ii<numpoints;ii++)
  {
    time_op.setValue(times[ii]);
    result[ii] = tableCopy.val();
  }

  for (int ii=0;ii<args.size();ii++) { delete args[ii]; }
 
  EXPECT_EQ(refRes,result);
}
#endif

#if 0
//-------------------------------------------------------------------------------
TEST ( Double_Ast_calculus_Test, ddx1)
{ 
  // define the .param
  double A=3.0;
  numval <double> paramValA(A);
  paramOp <double> testParamA(std::string("A"));
  testParamA.setNode(&paramValA);

  // define the expression that uses the .param
  double B=5.0;
  numval <double> paramValB(B);
  binaryAddOp <double> f_of_x(&paramValB,&testParamA);

  ddxOp<double> finalExp1(&f_of_x, &testParamA);
  EXPECT_EQ(finalExp1.val(), 1.0);
  ddxOp<double> finalExp2(&f_of_x, &testParamA);
  EXPECT_EQ(finalExp2.val(), 1.0);
}

TEST ( Double_Ast_calculus_Test, ddx2)
{ 
  // define the .params
  numval <double> paramValA(3.0), paramValB(4.0);
  paramOp <double> testParamA(std::string("A")), testParamB(std::string("B"));
  testParamA.setNode(&paramValA);
  testParamB.setNode(&paramValB);

  // define the expression that uses the .params
  binaryMulOp <double> f_of_x(&testParamA, &testParamB);

  ddxOp<double> finalExp1(&f_of_x, &testParamA);
  EXPECT_EQ(finalExp1.val(), 4.0);

  ddxOp<double> finalExp2(&f_of_x, &testParamB);
  EXPECT_EQ(finalExp2.val(), 3.0);
}


TEST ( Double_Ast_calculus_Test, ddx3)
{ 
  numval <double> paramValA(3.0);
  paramOp <double> testParamA(std::string("A"));
  testParamA.setNode(&paramValA);

  binaryMulOp <double> f_of_x(&testParamA, &testParamA);

  ddxOp<double> finalExp1(&f_of_x, &testParamA);
  EXPECT_EQ(finalExp1.val(), 6.0);
}

TEST ( Double_Ast_calculus_Test, ddx4)
{ 
  numval <double> paramValA(3.0);
  paramOp <double> testParamA(std::string("A"));
  testParamA.setNode(&paramValA);

  sinOp<double> f_of_x(&testParamA);

  ddxOp<double> finalExp1(&f_of_x, &testParamA);
  EXPECT_EQ(finalExp1.val(), std::cos(3.0));
}

TEST ( Double_Ast_modulus_Test, test1)
{ 
  numval <double> valA(15);
  numval <double> valB(4);
  binaryModOp<double> floor1(&valA,&valB);
  EXPECT_EQ(floor1.val(), 3);
}

TEST ( Double_Ast_modulus_Test, test2)
{ 
  numval <double> valA(30);
  numval <double> valB(7);
  binaryModOp<double> floor1(&valA,&valB);
  EXPECT_EQ(floor1.val(), 2);
}

TEST ( Double_Ast_modulus_Test, test3)
{ 
  numval <double> valA(20);
  numval <double> valB(7);
  binaryAddOp<double> plus1(&valA,&valB);

  numval <double> valC(2);
  numval <double> valD(8);
  binaryMulOp<double> mult1(&valC,&valD);

  binaryModOp<double> floor1(&plus1,&mult1);
  EXPECT_EQ(floor1.val(), 11);
}

TEST ( Double_Ast_floor_Test, test1)
{ 
//Floor of 10.25 = 10
  numval <double> valA(10.25);
  floorOp<double> floor1(&valA);
  EXPECT_EQ(floor1.val(), 10);
}

TEST ( Double_Ast_floor_Test, test2)
{ 
//Floor of -34.251 = -35
  numval <double> valA(-34.251);
  floorOp<double> floor1(&valA);
  EXPECT_EQ(floor1.val(), -35);
}

TEST ( Double_Ast_floor_Test, test3)
{ 
//Floor of 0.71 = 0
  numval <double> valA(0.71);
  floorOp<double> floor1(&valA);
  EXPECT_EQ(floor1.val(), 0);
}

TEST ( Double_Ast_floor_Test, test4)
{ 
//Floor of 10.25 = 10
  numval <double> valA(10.25);
  paramOp <double> testParamA(std::string("A"));
  testParamA.setNode(&valA);
  floorOp<double> floor1(&testParamA);

  EXPECT_EQ(floor1.val(), 10);
  ddxOp<double> ddx1(&floor1,&testParamA);
  EXPECT_EQ(ddx1.val(), 0);
}

TEST ( Double_Ast_ceil_Test, test1)
{ 
//Floor of 10.25 = 10
  numval <double> valA(10.25);
  ceilOp<double> ceil1(&valA);
  EXPECT_EQ(ceil1.val(), 11);
}

TEST ( Double_Ast_ceil_Test, test2)
{ 
//Floor of -34.251 = -35
  numval <double> valA(-34.251);
  ceilOp<double> ceil1(&valA);
  EXPECT_EQ(ceil1.val(), -34);
}

TEST ( Double_Ast_ceil_Test, test3)
{ 
//Floor of 0.71 = 0
  numval <double> valA(0.71);
  ceilOp<double> ceil1(&valA);
  EXPECT_EQ(ceil1.val(), 1);
}

TEST ( Double_Ast_ceil_Test, test4)
{ 
//Floor of 10.25 = 10
  numval <double> valA(10.25);
  paramOp <double> testParamA(std::string("A"));
  testParamA.setNode(&valA);
  ceilOp<double> ceil1(&testParamA);

  EXPECT_EQ(ceil1.val(), 11);
  ddxOp<double> ddx1(&ceil1,&testParamA);
  EXPECT_EQ(ddx1.val(), 0);
}


TEST ( Double_Ast_sgn_Test, test1)
{ 
// +1 if x > 0 
//  0 if x = 0 
// -1 if x < 0
  numval <double> valA(10.25);
  sgnOp<double> sgn1(&valA);
  EXPECT_EQ(sgn1.val(), 1);
}

TEST ( Double_Ast_sgn_Test, test2)
{ 
// +1 if x > 0 
//  0 if x = 0 
// -1 if x < 0
  numval <double> valA(0.0);
  sgnOp<double> sgn1(&valA);
  EXPECT_EQ(sgn1.val(), 0);
}

TEST ( Double_Ast_sgn_Test, test3)
{ 
// +1 if x > 0 
//  0 if x = 0 
// -1 if x < 0
  numval <double> valA(-7.5);
  sgnOp<double> sgn1(&valA);
  EXPECT_EQ(sgn1.val(), -1);
}

TEST ( Double_Ast_sign_Test, test1)
{ 
  // sign(x,y) = sgn(y)|x|  sign of y times absolute value of x
  numval <double> valX(-25);
  numval <double> valY(10.25);
  signOp<double> sign1(&valX,&valY);
  EXPECT_EQ(sign1.val(), 25);
}

TEST ( Double_Ast_sign_Test, test2)
{ 
  // sign(x,y) = sgn(y)|x|  sign of y times absolute value of x
  numval <double> valX(15);
  numval <double> valY(-10.25);
  signOp<double> sign1(&valX,&valY);
  EXPECT_EQ(sign1.val(), -15);
}

TEST ( Double_Ast_limit_Test, test1)
{ // limit(x,y,z)
  // x limited to range y to z
  // y if x < y
  // x if y < x < z
  // z if x > z
  numval <double> valX(-25); // test value
  numval <double> valY(1.25); // lower value
  numval <double> valZ(11.25); // upper value
  limitOp<double> limit1(&valX,&valY,&valZ);
  EXPECT_EQ(limit1.val(), 1.25);
}

TEST ( Double_Ast_limit_Test, test2)
{ // limit(x,y,z)
  // x limited to range y to z
  // y if x < y
  // x if y < x < z
  // z if x > z
  numval <double> valX(+5); // test value
  numval <double> valY(1.25); // lower value
  numval <double> valZ(11.25); // upper value
  limitOp<double> limit1(&valX,&valY,&valZ);
  EXPECT_EQ(limit1.val(), 5);
}

TEST ( Double_Ast_limit_Test, test3)
{ // limit(x,y,z)
  // x limited to range y to z
  // y if x < y
  // x if y < x < z
  // z if x > z
  numval <double> valX(+17); // test value
  numval <double> valY(1.25); // lower value
  numval <double> valZ(11.25); // upper value
  limitOp<double> limit1(&valX,&valY,&valZ);
  EXPECT_EQ(limit1.val(), 11.25);
}

TEST ( Double_Ast_int_Test, test1)
{ 
  // int(x)
  // integer part of the real variable x 
  numval <double> valX(11.2423); // test value
  intOp<double> int1(&valX);
  EXPECT_EQ(int1.val(), 11);
}

TEST ( Double_Ast_int_Test, test2)
{ 
  // int(x)
  // integer part of the real variable x 
  numval <double> valX(-11.2423); // test value
  intOp<double> int1(&valX);
  EXPECT_EQ(int1.val(),-11);
}
#endif

int main (int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}


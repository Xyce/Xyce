//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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

int main (int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}


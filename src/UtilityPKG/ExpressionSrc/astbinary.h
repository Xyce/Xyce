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

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 10/xx/2019
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef _astbinary_h_
#define _astbinary_h_

#define AST_BIN_OP_MACRO(NAME,FCTQUOTE,FCTCODE,VAL,DX)                                \
template <typename ScalarT>                                                            \
class NAME : public astNode<ScalarT>                                                   \
{                                                                                      \
  public:                                                                              \
    NAME (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):   \
        astNode<ScalarT>(left,right) {                                                 \
          rightConst_ = this->rightAst_->numvalType();                                 \
          leftConst_ = this->leftAst_->numvalType();                                   \
        };                                                                             \
                                                                                       \
    virtual ScalarT val(){                                                             \
      ScalarT leftVal=this->leftAst_->val();                                           \
      ScalarT rightVal=this->rightAst_->val();                                         \
      return VAL; }                                                                    \
                                                                                       \
    virtual ScalarT dx(int i) {                                                        \
      ScalarT leftVal=this->leftAst_->val();                                           \
      ScalarT rightVal=this->rightAst_->val();                                         \
      ScalarT leftDx, rightDx;                                                         \
      if (!leftConst_) { leftDx =this->leftAst_->dx(i); } else { leftDx=0.0; }         \
      if (!rightConst_) { rightDx =this->rightAst_->dx(i); } else {rightDx=0.0; }      \
      return DX; }                                                                     \
                                                                                       \
    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)                  \
    {                                                                                  \
      int numDerivs = derivs.size();                                                   \
      ScalarT leftVal, rightVal, leftDx=0.0, rightDx=0.0;                              \
      if (leftConst_) { leftVal = this->leftAst_->val(); }                             \
      else {                                                                           \
        if (lefDerivs_.empty()) { lefDerivs_.resize(numDerivs,0.0); }                  \
        this->leftAst_->dx2(leftVal,lefDerivs_); }                                     \
      if (rightConst_) { rightVal = this->rightAst_->val(); }                          \
      else {                                                                           \
        if (rigDerivs_.empty()) { rigDerivs_.resize(numDerivs,0.0); }                  \
        this->rightAst_->dx2(rightVal,rigDerivs_); }                                   \
      result=VAL;                                                                      \
      for (int i=0;i<numDerivs;i++) {                                                  \
        if (!leftConst_)  { leftDx = lefDerivs_[i]; }                                  \
        if (!rightConst_) { rightDx = rigDerivs_[i]; }                                \
        derivs[i] = DX; }                                                              \
    }                                                                                  \
                                                                                       \
    virtual void output(std::ostream & os, int indent=0)                               \
    {                                                                                  \
      os << std::setw(indent) << " ";                                                  \
      os << FCTQUOTE << " id = " << this->id_ << std::endl;                            \
      ++indent;                                                                        \
      this->leftAst_->output(os,indent+1);                                             \
      this->rightAst_->output(os,indent+1);                                            \
    }                                                                                  \
                                                                                       \
    virtual void compactOutput(std::ostream & os)                                      \
    { os << FCTQUOTE << " id = " << this->id_ << std::endl; }                          \
                                                                                       \
    virtual void codeGen (std::ostream & os )                                          \
    {                                                                                  \
      os << "(";                                                                       \
      this->leftAst_->codeGen(os);                                                     \
      os << FCTCODE;                                                                   \
      this->rightAst_->codeGen(os);                                                    \
      os << ")";                                                                       \
    }                                                                                  \
    bool rightConst_;                                                                  \
    bool leftConst_;                                                                   \
    std::vector<ScalarT> lefDerivs_;                                                   \
    std::vector<ScalarT> rigDerivs_;                                                   \
};

AST_BIN_OP_MACRO(
    binaryAddOp,
    "binary add ",
	  "+",
	  (leftVal + rightVal),
   (rightConst_)? (leftConst_?(0.0):(leftDx)) : (leftConst_?(rightDx):(leftDx+rightDx))
    )
AST_BIN_OP_MACRO(
    binaryMinusOp,
    "binary minus ",
	  "-",
	  (leftVal - rightVal),
    (rightConst_)?(leftConst_?(0.0):(leftDx)):(leftConst_?(-rightDx):(leftDx-rightDx))
    )
AST_BIN_OP_MACRO(
    binaryModOp,
    "modulus operator ",
	  "%",
	  (static_cast<int>(std::real(leftVal)) % static_cast<int>(std::real(rightVal))),
	  (0.0)
    )
AST_BIN_OP_MACRO(
    binaryMulOp,
    "binary multiply ",
	  "*",
	  (leftVal * rightVal),
    (rightConst_)?(leftConst_?(0.0):(leftDx * rightVal)):(leftConst_?(rightDx * leftVal):(leftDx*rightVal+rightDx*leftVal))
    )
AST_BIN_OP_MACRO(
    binaryDivOp,
    "binary division ",
	  "/",
	  (leftVal / rightVal),
    (rightConst_)?(leftConst_?(0.0):((leftDx*rightVal)/(rightVal*rightVal))):(leftConst_?((-rightDx*leftVal)/(rightVal*rightVal)):((leftDx*rightVal-rightDx*leftVal)/(rightVal*rightVal)))
    )

#endif


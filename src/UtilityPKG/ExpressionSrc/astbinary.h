//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
        astNode<ScalarT>(left,right), localDerivsSizeLef_(0), localDerivsSizeRig_(0)   \
  {                                                 \
          rightConst_ = this->childrenAstNodes_[1]->numvalType();                      \
          leftConst_ = this->childrenAstNodes_[0]->numvalType();                       \
        };                                                                             \
                                                                                       \
    virtual ScalarT val(){                                                             \
      ScalarT leftVal=this->childrenAstNodes_[0]->val();                               \
      ScalarT rightVal=this->childrenAstNodes_[1]->val();                              \
      return VAL; }                                                                    \
                                                                                       \
    virtual ScalarT dx(int i) {                                                        \
      ScalarT leftVal=this->childrenAstNodes_[0]->val();                               \
      ScalarT rightVal=this->childrenAstNodes_[1]->val();                              \
      ScalarT leftDx, rightDx;                                                         \
      if (!leftConst_) { leftDx =this->childrenAstNodes_[0]->dx(i); } else { leftDx=0.0; }         \
      if (!rightConst_) { rightDx =this->childrenAstNodes_[1]->dx(i); } else {rightDx=0.0; }      \
      return DX; }                                                                     \
                                                                                       \
    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs, int numDerivs)   \
    {                                                                                  \
      ScalarT leftVal, rightVal, leftDx=0.0, rightDx=0.0;                              \
      if (leftConst_) { leftVal = this->childrenAstNodes_[0]->val(); }                 \
      else {                                                                           \
        if (localDerivsSizeLef_ < numDerivs) { lefDerivs_.resize(numDerivs,0.0); localDerivsSizeLef_ = numDerivs; } \
        this->childrenAstNodes_[0]->dx2(leftVal,lefDerivs_,numDerivs); }               \
      if (rightConst_) { rightVal = this->childrenAstNodes_[1]->val(); }               \
      else {                                                                           \
        if (localDerivsSizeRig_ < numDerivs) { rigDerivs_.resize(numDerivs,0.0);  localDerivsSizeRig_ = numDerivs;}      \
        this->childrenAstNodes_[1]->dx2(rightVal,rigDerivs_,numDerivs); }              \
      result=VAL;                                                                      \
      for (int i=0;i<numDerivs;i++) {                                                  \
        if (!leftConst_)  { leftDx = lefDerivs_[i]; }                                  \
        if (!rightConst_) { rightDx = rigDerivs_[i]; }                                 \
        derivs[i] = DX; }                                                              \
    }                                                                                  \
                                                                                       \
    virtual bool getIsComplex ()                                                       \
    { return (this->childrenAstNodes_[1]->getIsComplex() || this->childrenAstNodes_[0]->getIsComplex()); }    \
                                                                                       \
    virtual void generateExpressionString (std::string & str)                          \
    {                                                                                  \
      std::string tmp1,tmp2;                                                           \
      this->childrenAstNodes_[0]->generateExpressionString(tmp1);                      \
      this->childrenAstNodes_[1]->generateExpressionString(tmp2);                      \
      str = "(" + tmp1 + std::string(FCTCODE) + tmp2 + ")";                            \
    }                                                                                  \
                                                                                       \
    virtual void output(std::ostream & os, int indent=0)                               \
    {                                                                                  \
      os << std::setw(indent) << " ";                                                  \
      os << FCTQUOTE << " id = " << this->id_ << std::endl;                            \
      ++indent;                                                                        \
      this->childrenAstNodes_[0]->output(os,indent+1);                                 \
      this->childrenAstNodes_[1]->output(os,indent+1);                                 \
    }                                                                                  \
                                                                                       \
    virtual void compactOutput(std::ostream & os)                                      \
    { os << FCTQUOTE << " id = " << this->id_ << std::endl; }                          \
                                                                                       \
    virtual void codeGen (std::ostream & os )                                          \
    {                                                                                  \
      os << "(";                                                                       \
      this->childrenAstNodes_[0]->codeGen(os);                                         \
      os << FCTCODE;                                                                   \
      this->childrenAstNodes_[1]->codeGen(os);                                         \
      os << ")";                                                                       \
    }                                                                                  \
    virtual bool getIsTreeConstant() { return                                          \
     (this->childrenAstNodes_[0]->getIsTreeConstant() && this->childrenAstNodes_[0]->getIsTreeConstant()); }   \
                                                                                       \
    virtual void accept                                                                \
          (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) \
    { Teuchos::RCP<NAME<ScalarT> > castToThis = Teuchos::rcp_static_cast<NAME<ScalarT> > (thisAst_); \
      visitor.visit( castToThis );  \
      this->childrenAstNodes_[0]->accept(visitor, this->childrenAstNodes_[0]);         \
      this->childrenAstNodes_[1]->accept(visitor, this->childrenAstNodes_[1]); }       \
                                                                                       \
    bool rightConst_;                                                                  \
    bool leftConst_;                                                                   \
    std::vector<ScalarT> lefDerivs_;                                                   \
    std::vector<ScalarT> rigDerivs_;                                                   \
    int localDerivsSizeLef_; \
    int localDerivsSizeRig_; \
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


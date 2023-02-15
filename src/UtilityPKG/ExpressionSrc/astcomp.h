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

#ifndef _astcomp_h_
#define _astcomp_h_

#define AST_CMP_OP_MACRO(NAME,FCTQUOTE,VAL,DX)                                              \
template <typename ScalarT>                                                                 \
class NAME : public astNode<ScalarT>                                                        \
{                                                                                           \
  public:                                                                                   \
    NAME (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):   \
        astNode<ScalarT>(left,right), bpTol_(0.0) {};                                       \
                                                                                            \
    virtual ScalarT val()                                                                   \
    {                                                                                       \
      bpTimes_.clear();                                                                     \
      computeBreakPoint ( this->leftAst_, this->rightAst_, timeOpVec_, bpTol_, bpTimes_);   \
      return VAL;                                                                           \
    }                                                                                       \
                                                                                            \
    virtual ScalarT dx(int i) { return DX; }                                                \
                                                                                            \
    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)  \
    { \
      result = val(); \
      std::fill(derivs.begin(),derivs.end(),0.0); \
    } \
                                                                                            \
    virtual bool getBreakPoints(std::vector<Xyce::Util::BreakPoint> & breakPointTimes)      \
    {                                                                                       \
      if(!(bpTimes_.empty()))                                                               \
      {                                                                                     \
        for (int ii=0;ii<bpTimes_.size();ii++)                                              \
        {                                                                                   \
          breakPointTimes.push_back(bpTimes_[ii]);                                          \
        }                                                                                   \
      }                                                                                     \
      return true;                                                                          \
    }                                                                                       \
    virtual void setBreakPointTol(double tol) { bpTol_ = tol; }                             \
    virtual bool getIsComplex () { return false; }                                          \
    virtual void output(std::ostream & os, int indent=0)                                    \
    {                                                                                       \
      os << std::setw(indent) << " ";                                                       \
      os << FCTQUOTE " operator id = " << this->id_ << std::endl;                           \
      ++indent;                                                                             \
      this->leftAst_->output(os,indent+1);                                                  \
      this->rightAst_->output(os,indent+1);                                                 \
    }                                                                                       \
    virtual void compactOutput(std::ostream & os)                                           \
    { os << FCTQUOTE " operator id = " << this->id_ << std::endl; }                         \
                                                                                            \
    virtual bool getIsTreeConstant() { return                                               \
     (this->leftAst_->getIsTreeConstant() && this->leftAst_->getIsTreeConstant()); }        \
    virtual bool compType() { return true; }                                                \
                                                                                            \
    virtual void codeGen (std::ostream & os )                                               \
    {                                                                                       \
      os << "(";                                                                            \
      this->leftAst_->codeGen(os);                                                          \
      os << FCTQUOTE;                                                                       \
      this->rightAst_->codeGen(os);                                                         \
      os << ")";                                                                            \
    }                                                                                       \
    std::vector<Teuchos::RCP<astNode<ScalarT> > > timeOpVec_;                               \
    double bpTol_;                                                                          \
    std::vector<Xyce::Util::BreakPoint> bpTimes_;                                           \
};       


#define AST_CMP_OP_MACRO(NAME,FCTQUOTE,VAL,DX)                                              \
template <typename ScalarT>                                                                 \
class NAME : public astNode<ScalarT>                                                        \
{                                                                                           \
  public:                                                                                   \
    NAME (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):   \
        astNode<ScalarT>(left,right), bpTol_(0.0) {};                                       \
                                                                                            \
    virtual ScalarT val()                                                                   \
    {                                                                                       \
      bpTimes_.clear();                                                                     \
      computeBreakPoint ( this->leftAst_, this->rightAst_, timeOpVec_, bpTol_, bpTimes_);   \
      return VAL;                                                                           \
    }                                                                                       \
                                                                                            \
    virtual ScalarT dx(int i) { return DX; }                                                \
                                                                                            \
    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)  \
    { \
      result = val(); \
      std::fill(derivs.begin(),derivs.end(),0.0); \
    } \
                                                                                            \
    virtual bool getBreakPoints(std::vector<Xyce::Util::BreakPoint> & breakPointTimes)      \
    {                                                                                       \
      if(!(bpTimes_.empty()))                                                               \
      {                                                                                     \
        for (int ii=0;ii<bpTimes_.size();ii++)                                              \
        {                                                                                   \
          breakPointTimes.push_back(bpTimes_[ii]);                                          \
        }                                                                                   \
      }                                                                                     \
      return true;                                                                          \
    }                                                                                       \
    virtual void setBreakPointTol(double tol) { bpTol_ = tol; }                             \
    virtual bool getIsComplex () { return false; }                                          \
    virtual void output(std::ostream & os, int indent=0)                                    \
    {                                                                                       \
      os << std::setw(indent) << " ";                                                       \
      os << FCTQUOTE " operator id = " << this->id_ << std::endl;                           \
      ++indent;                                                                             \
      this->leftAst_->output(os,indent+1);                                                  \
      this->rightAst_->output(os,indent+1);                                                 \
    }                                                                                       \
    virtual void compactOutput(std::ostream & os)                                           \
    { os << FCTQUOTE " operator id = " << this->id_ << std::endl; }                         \
                                                                                            \
    virtual bool getIsTreeConstant() { return                                               \
     (this->leftAst_->getIsTreeConstant() && this->leftAst_->getIsTreeConstant()); }        \
    virtual bool compType() { return true; }                                                \
                                                                                            \
    virtual void codeGen (std::ostream & os )                                               \
    {                                                                                       \
      os << "(";                                                                            \
      this->leftAst_->codeGen(os);                                                          \
      os << FCTQUOTE;                                                                       \
      this->rightAst_->codeGen(os);                                                         \
      os << ")";                                                                            \
    }                                                                                       \
    std::vector<Teuchos::RCP<astNode<ScalarT> > > timeOpVec_;                               \
    double bpTol_;                                                                          \
    std::vector<Xyce::Util::BreakPoint> bpTimes_;                                           \
};


#define AST_CMP_OP_MACRO2(NAME,FCTQUOTE,VAL,DX)                                             \
template <typename ScalarT>                                                                 \
class NAME : public astNode<ScalarT>                                                        \
{                                                                                           \
  public:                                                                                   \
    NAME (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):   \
        astNode<ScalarT>(left,right) {};                                                    \
                                                                                            \
    virtual ScalarT val(){return VAL; }                                                     \
                                                                                            \
    virtual ScalarT dx(int i) { return DX; }                                                \
                                                                                            \
    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)                       \
    {                                                                                       \
      result = val();                                                                       \
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0); }             \
    }                                                                                       \
                                                                                            \
    virtual bool getIsTreeConstant() { return                                               \
     (this->leftAst_->getIsTreeConstant() && this->leftAst_->getIsTreeConstant()); }        \
                                                                                            \
    virtual bool getIsComplex () { return false; }                                          \
                                                                                            \
    virtual void output(std::ostream & os, int indent=0)                                    \
    {                                                                                       \
      os << std::setw(indent) << " ";                                                       \
      os << FCTQUOTE " operator id = " << this->id_ << std::endl;                           \
      ++indent;                                                                             \
      this->leftAst_->output(os,indent+1);                                                  \
      this->rightAst_->output(os,indent+1);                                                 \
    }                                                                                       \
    virtual void compactOutput(std::ostream & os)                                           \
    { os << FCTQUOTE " operator id = " << this->id_ << std::endl; }                         \
                                                                                            \
    virtual void codeGen (std::ostream & os )                                               \
    {                                                                                       \
      os << "(";                                                                            \
      this->leftAst_->codeGen(os);                                                          \
      os << FCTQUOTE;                                                                       \
      this->rightAst_->codeGen(os);                                                         \
      os << ")";                                                                            \
    }                                                                                       \
}; 

#define FIXVAL(val)  Xyce::Util::fixNan(  Xyce::Util::fixInf( val ) )
#define FIXVALREAL(val)  std::real( FIXVAL(val) )

AST_CMP_OP_MACRO(  gtOp,  ">", ((FIXVALREAL(this->leftAst_->val()) > FIXVALREAL(this->rightAst_->val()))? 1 : 0), (0.0))
AST_CMP_OP_MACRO(  ltOp,  "<", ((FIXVALREAL(this->leftAst_->val()) < FIXVALREAL(this->rightAst_->val()))? 1 : 0), (0.0))
AST_CMP_OP_MACRO(  neOp, "!=", (((FIXVAL(this->leftAst_->val())) != (FIXVAL(this->rightAst_->val())))? 1 : 0), (0.0))
AST_CMP_OP_MACRO(  eqOp, "==", (((FIXVAL(this->leftAst_->val())) == (FIXVAL(this->rightAst_->val())))? 1 : 0), (0.0))
AST_CMP_OP_MACRO(  geOp, ">=", ((FIXVALREAL(this->leftAst_->val()) >= FIXVALREAL(this->rightAst_->val()))? 1 : 0), (0.0))
AST_CMP_OP_MACRO(  leOp, "<=", ((FIXVALREAL(this->leftAst_->val()) <= FIXVALREAL(this->rightAst_->val()))? 1 : 0), (0.0))

AST_CMP_OP_MACRO2(  orOp, "||", ((FIXVALREAL(this->leftAst_->val()) || FIXVALREAL(this->rightAst_->val()))? 1 : 0), (0.0))
AST_CMP_OP_MACRO2( andOp, "&&", ((FIXVALREAL(this->leftAst_->val()) && FIXVALREAL(this->rightAst_->val()))? 1 : 0), (0.0))
AST_CMP_OP_MACRO2( xorOp,  "^", (((FIXVALREAL(this->leftAst_->val()) > 0 && FIXVALREAL(this->rightAst_->val()) <= 0) || (FIXVALREAL(this->leftAst_->val()) <= 0 && FIXVALREAL(this->rightAst_->val()) > 0))?1.0:0.0), (0.0))
#endif


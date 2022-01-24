//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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

#ifndef _astfuncs_h_
#define _astfuncs_h_

#define AST_OP_MACRO(NAME,DX)                                                          \
template <typename ScalarT>                                                            \
class NAME ## Op : public astNode<ScalarT>                                             \
{                                                                                      \
  public:                                                                              \
    NAME ## Op (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {       \
          leftConst_ = this->leftAst_->numvalType(); };                                \
                                                                                       \
    virtual ScalarT val() { return std::NAME(this->leftAst_->val()); }                 \
                                                                                       \
    virtual ScalarT dx(int i)                                                          \
    {                                                                                  \
      if(leftConst_) { return 0.0; }                                                   \
      ScalarT leftVal=this->leftAst_->val();                                           \
      ScalarT leftDx =this->leftAst_->dx(i);                                           \
      return DX;                                                                       \
    }                                                                                  \
                                                                                       \
    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)                  \
    {                                                                                  \
      ScalarT leftVal, leftDx;                                                         \
      if (leftConst_)                                                                  \
      {                                                                                \
        leftVal=this->leftAst_->val();                                                 \
        result = std::NAME(leftVal);                                                   \
        std::fill(derivs.begin(),derivs.end(),0.0);                                    \
        return;                                                                        \
      }                                                                                \
      else                                                                             \
      {                                                                                \
        int numDerivs=derivs.size();                                                   \
        if (lefDerivs_.empty()) { lefDerivs_.resize(numDerivs,0.0); }                  \
        this->leftAst_->dx2(leftVal,lefDerivs_);                                       \
        result= std::NAME(leftVal);                                                    \
        for (int i=0;i<numDerivs;i++) { ScalarT leftDx=lefDerivs_[i]; derivs[i]=DX; }  \
      }                                                                                \
    }                                                                                  \
                                                                                       \
    virtual void output(std::ostream & os, int indent=0)                               \
    {                                                                                  \
      os << std::setw(indent) << " ";                                                  \
      os << #NAME << " operator id = " << this->id_ << std::endl;                      \
      ++indent;                                                                        \
      this->leftAst_->output(os,indent+1);                                             \
    }                                                                                  \
    virtual void compactOutput(std::ostream & os)                                      \
    {                                                                                  \
      os << #NAME << " operator id = " << this->id_ << std::endl;                      \
    }                                                                                  \
    virtual void codeGen (std::ostream & os )                                          \
    {                                                                                  \
      os << "std::" #NAME << "(";                                                      \
      this->leftAst_->codeGen(os);                                                     \
      os << ")";                                                                       \
    }                                                                                  \
    bool leftConst_;                                                                   \
    std::vector<ScalarT> lefDerivs_;                                                   \
};

AST_OP_MACRO( sqrt, (leftDx/(2.*std::sqrt(leftVal))))
AST_OP_MACRO( exp, (leftDx*std::exp(leftVal)))
AST_OP_MACRO( abs, (std::real(leftVal) >= 0 ? leftDx : ScalarT(-leftDx)))
AST_OP_MACRO( sin, (leftDx*std::cos(leftVal)))
AST_OP_MACRO( cos, (-leftDx*std::sin( leftVal)))
AST_OP_MACRO( acos, ( -leftDx/std::sqrt(1.-leftVal*leftVal) ))
AST_OP_MACRO( acosh, (leftDx/std::sqrt((leftVal-1.)*(leftVal+1.))))
AST_OP_MACRO( asin, (leftDx/std::sqrt(1.-leftVal*leftVal)))
AST_OP_MACRO( asinh, (leftDx/std::sqrt(1.+leftVal*leftVal)))
AST_OP_MACRO( cosh, (leftDx*std::sinh(leftVal)))
AST_OP_MACRO( log, (leftDx/leftVal))
AST_OP_MACRO( log10, (leftDx/(std::log(ScalarT(10))*leftVal)))
AST_OP_MACRO( sinh, (leftDx*std::cosh(leftVal)))
AST_OP_MACRO( tan, (leftDx*(1.+std::tan(leftVal)*std::tan(leftVal))))
AST_OP_MACRO( atan, (leftDx/(1.+leftVal*leftVal)))

template <typename ScalarT> 
class tanhOp : public astNode<ScalarT> 
{ 
  public: 
  tanhOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

  virtual ScalarT val() 
  { 
    ScalarT retval;
    ScalarT arg = this->leftAst_->val();
    if      (std::real(arg) > +20) { retval = +1.0; }
    else if (std::real(arg) < -20) { retval = -1.0; }
    else                           { retval = std::tanh(arg); }
    return retval;
  } 

  virtual ScalarT dx(int i) 
  { 
    ScalarT retdx=0.0;
    ScalarT arg = this->leftAst_->val();
    if (std::real(arg) <= 20 && std::real(arg) >= -20)
    {
      ScalarT cosh_arg = std::cosh(arg);
      retdx = (this->leftAst_->dx(i)/(cosh_arg*cosh_arg));
    }
    return retdx;
  } 

  virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs) 
  {
    int numDerivs = derivs.size();
    if (lefDerivs_.empty()) { lefDerivs_.resize(numDerivs,0.0); }

    ScalarT arg;
    this->leftAst_->dx2(arg,lefDerivs_);

    if      (std::real(arg) > +20) { result = +1.0; }
    else if (std::real(arg) < -20) { result = -1.0; }
    else                           { result = std::tanh(arg); }

    if (std::real(arg) <= 20 && std::real(arg) >= -20)
    {
      ScalarT cosh_arg = std::cosh(arg);
      for (int i=0;i<numDerivs;i++)
      {
        derivs[i] = (lefDerivs_[i]/(cosh_arg*cosh_arg));
      }
    }
    else
    {
      std::fill(derivs.begin(),derivs.end(),0.0);
    }
  }

  virtual void output(std::ostream & os, int indent=0) 
  { 
    os << std::setw(indent) << " ";
    os << "tanh" << " operator id = " << this->id_ << std::endl;
    ++indent;
    this->leftAst_->output(os,indent+1);
  } 

  virtual void compactOutput(std::ostream & os)
  { 
    os << "tanh" << " operator id = " << this->id_ << std::endl;
  }

  virtual void codeGen (std::ostream & os ) 
  { 
    os << "std::" "tanh" << "(";
    this->leftAst_->codeGen(os);
    os << ")";
  } 
  std::vector<ScalarT> lefDerivs_;
};

template <typename ScalarT> 
class atanhOp : public astNode<ScalarT> 
{
  public: 
  atanhOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

  virtual ScalarT val() 
  { 
    ScalarT Epsilon = 1.e-12;
    ScalarT arg = this->leftAst_->val();

    if      (std::real(arg) < std::real(Epsilon) - 1.0) { arg = Epsilon - 1.0; }
    else if (std::real(arg) > 1.0 - std::real(Epsilon)) { arg = 1.0 - Epsilon; }

    return std::atanh(arg); // old expression library returns (log((1.0 + arg) / (1.0 - arg)) / 2.0)
    //return (log((1.0 + arg) / (1.0 - arg)) / 2.0);  // ERK.  tried this, in hopes it would fixed bug 254 test. didn't help
  } 

  virtual ScalarT dx(int i) 
  { 
    ScalarT Epsilon = 1.e-12;
    ScalarT retdx=0.0;
    ScalarT arg = this->leftAst_->val();
    if (std::real(arg) >= (std::real(Epsilon) - 1.0) && std::real(arg) <= (1.0 - std::real(Epsilon)))
    {
      retdx = (this->leftAst_->dx(i)/(1.- arg*arg));
    }

    return retdx;
  } 

  virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs) 
  {
    ScalarT Epsilon = 1.e-12;


    // derivs
    int numDerivs = derivs.size();
    if (lefDerivs_.empty()) { lefDerivs_.resize(numDerivs,0.0); }

    ScalarT arg;
    this->leftAst_->dx2(arg,lefDerivs_);
    if (std::real(arg) >= (std::real(Epsilon) - 1.0) && std::real(arg) <= (1.0 - std::real(Epsilon)))
    {
      for (int i=0;i<numDerivs;i++) { derivs[i] = (lefDerivs_[i]/(1.- arg*arg)); }
    }
    else { std::fill(derivs.begin(),derivs.end(),0.0); }

    // result (this changes arg, so do second)
    if      (std::real(arg) < std::real(Epsilon) - 1.0) { arg = Epsilon - 1.0; }
    else if (std::real(arg) > 1.0 - std::real(Epsilon)) { arg = 1.0 - Epsilon; }

    result = std::atanh(arg);

  }

  virtual void output(std::ostream & os, int indent=0) 
  { 
   os << std::setw(indent) << " ";
   os << "atanh" << " operator id = " << this->id_ << std::endl;
   ++indent;
   this->leftAst_->output(os,indent+1);
  } 

  virtual void compactOutput(std::ostream & os)
  { 
   os << "atanh" << " operator id = " << this->id_ << std::endl;
  }

  virtual void codeGen (std::ostream & os ) 
  { 
    os << "std::" "atanh" << "(";
    this->leftAst_->codeGen(os);
    os << ")";
  } 
  std::vector<ScalarT> lefDerivs_;
};

#endif

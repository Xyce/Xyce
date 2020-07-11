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
    NAME ## Op (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};     \
                                                                                       \
    virtual ScalarT val() { return std::NAME(this->leftAst_->val()); }                 \
                                                                                       \
    virtual ScalarT dx(int i) { return DX; }                                           \
                                                                                       \
    virtual void output(std::ostream & os, int indent=0)                               \
    {                                                                                  \
      os << std::setw(indent) << " ";                                                  \
      os << #NAME << " operator " << std::endl;                                        \
      ++indent;                                                                        \
      this->leftAst_->output(os,indent+1);                                             \
    }                                                                                  \
                                                                                       \
    virtual void codeGen (std::ostream & os )                                          \
    {                                                                                  \
      os << "std::" #NAME << "(";                                                      \
      this->leftAst_->codeGen(os);                                                     \
      os << ")";                                                                       \
    }                                                                                  \
};

AST_OP_MACRO( sqrt, (this->leftAst_->dx(i)/(2.*std::sqrt(this->leftAst_->val()))))
AST_OP_MACRO( exp, (this->leftAst_->dx(i)*std::exp(this->leftAst_->val())))
AST_OP_MACRO( abs, (std::real(this->leftAst_->val()) >= 0 ? this->leftAst_->dx(i) : ScalarT(-this->leftAst_->dx(i))))
AST_OP_MACRO( sin, (this->leftAst_->dx(i)*std::cos(this->leftAst_->val())))
AST_OP_MACRO( cos, (-this->leftAst_->dx(i)*std::sin( this->leftAst_->val())))
AST_OP_MACRO( acos, ( -this->leftAst_->dx(i)/std::sqrt(1.-this->leftAst_->val()*this->leftAst_->val()) ))
AST_OP_MACRO( acosh, (this->leftAst_->dx(i)/std::sqrt((this->leftAst_->val()-1.)*(this->leftAst_->val()+1.))))
AST_OP_MACRO( asin, (this->leftAst_->dx(i)/std::sqrt(1.-this->leftAst_->val()*this->leftAst_->val())))
AST_OP_MACRO( asinh, (this->leftAst_->dx(i)/std::sqrt(1.+this->leftAst_->val()*this->leftAst_->val())))
AST_OP_MACRO( cosh, (this->leftAst_->dx(i)*std::sinh(this->leftAst_->val())))
AST_OP_MACRO( log, (this->leftAst_->dx(i)/this->leftAst_->val()))
AST_OP_MACRO( log10, (this->leftAst_->dx(i)/(std::log(ScalarT(10))*this->leftAst_->val())))
AST_OP_MACRO( sinh, (this->leftAst_->dx(i)*std::cosh(this->leftAst_->val())))
AST_OP_MACRO( tan, (this->leftAst_->dx(i)*(1.+std::tan(this->leftAst_->val())*std::tan(this->leftAst_->val()))))
AST_OP_MACRO( atan, (this->leftAst_->dx(i)/(1.+this->leftAst_->val()*this->leftAst_->val())))

//AST_OP_MACRO( tanh, (this->leftAst_->dx(i)/(std::cosh(this->leftAst_->val())*std::cosh(this->leftAst_->val()))))
//AST_OP_MACRO( atanh, (this->leftAst_->dx(i)/(1.-this->leftAst_->val()*this->leftAst_->val())))

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
      retdx = (this->leftAst_->dx(i)/(std::cosh(this->leftAst_->val())*std::cosh(this->leftAst_->val())));
    }
    return retdx;
  } 

  virtual void output(std::ostream & os, int indent=0) 
  { 
    os << std::setw(indent) << " ";
    os << "tanh" << " operator " << std::endl;
    ++indent;
    this->leftAst_->output(os,indent+1);
  } 

  virtual void codeGen (std::ostream & os ) 
  { 
    os << "std::" "tanh" << "(";
    this->leftAst_->codeGen(os);
    os << ")";
  } 
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

   // return std::atanh(arg); // old expression library returns (log((1.0 + arg) / (1.0 - arg)) / 2.0)
    return (log((1.0 + arg) / (1.0 - arg)) / 2.0);  // ERK.  tried this, in hopes it would fixed bug 254 test. didn't help
  } 

  virtual ScalarT dx(int i) 
  { 
    ScalarT Epsilon = 1.e-12;
    ScalarT retdx=0.0;
    ScalarT arg = this->leftAst_->val();

    //if (std::real(arg) >= (std::real(Epsilon) - 1.0) && std::real(arg) <= (1.0 - std::real(Epsilon)))
    {
      retdx = (this->leftAst_->dx(i)/(1.-this->leftAst_->val()*this->leftAst_->val()));
    }

    return retdx;
  } 

  virtual void output(std::ostream & os, int indent=0) 
  { 
   os << std::setw(indent) << " ";
   os << "atanh" << " operator " << std::endl;
   ++indent;
   this->leftAst_->output(os,indent+1);
  } 

  virtual void codeGen (std::ostream & os ) 
  { 
    os << "std::" "atanh" << "(";
    this->leftAst_->codeGen(os);
    os << ")";
  } 
};

#endif

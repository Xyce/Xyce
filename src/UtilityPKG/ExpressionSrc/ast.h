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

#ifndef ast_H
#define ast_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <list>
#include <unordered_map>
#include <string>
#include <complex>

#include <Teuchos_RCP.hpp>

#include <N_UTL_BreakPoint.h>
#include <N_UTL_Interpolators.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_CheckIfValidFile.h>
#include <N_ERH_Message.h>
#include <N_UTL_HspiceBools.h>
#include <N_UTL_Math.h>
#include <expressionParamTypes.h>

#define CONSTCtoK    (273.15)

namespace Xyce {
namespace Util {

template <typename ScalarA>
inline void fixNan(ScalarA & result) { if (std::isnan(std::real(result))) { result = 1.0e+50; } }

template <>
inline void fixNan(std::complex<double> & result)
{
  bool negReal = std::signbit(std::real(result));
  bool negImag = std::signbit(std::imag(result));
  if (std::isnan(std::real(result))) {std::complex<double> tmp = std::complex<double>((1.0e+50)*(negReal?-1.0:1.0),result.imag());result = tmp;}
  if (std::isnan(std::imag(result))) {std::complex<double> tmp = std::complex<double>(result.real(),(1.0e+50)*(negImag?-1.0:1.0));result = tmp;}
}

template <typename ScalarA>
inline ScalarA fixNan(const ScalarA & result)
{
  ScalarA tmp = result;
  if (std::isnan(std::real(result))) { tmp = 1.0e+50; }
  return tmp;
}

template <>
inline std::complex<double> fixNan(const std::complex<double> & result)
{
  bool negReal = std::signbit(std::real(result));
  bool negImag = std::signbit(std::imag(result));
  std::complex<double> tmp = result;
  if (std::isnan(std::real(result))) { tmp = std::complex<double>((1.0e+50)*(negReal?-1.0:1.0),result.imag()); }
  if (std::isnan(std::imag(result))) { tmp = std::complex<double>(result.real(),(1.0e+50)*(negImag?-1.0:1.0)); }
  return tmp;
}

template <typename ScalarA>
inline void fixInf(ScalarA & result) { if (std::isinf(std::real(result))) { bool neg = std::signbit(result); result = (1.0e+50)*(neg?-1.0:1.0); } }

template <>
inline void fixInf(std::complex<double> & result)
{
  bool negReal = std::signbit(std::real(result));
  bool negImag = std::signbit(std::imag(result));

  if (std::isinf(std::real(result))) {std::complex<double> tmp = std::complex<double>((1.0e+50)*(negReal?-1.0:1.0),result.imag());result = tmp;}
  if (std::isinf(std::imag(result))) {std::complex<double> tmp = std::complex<double>(result.real(),(1.0e+50)*(negImag?-1.0:1.0));result = tmp;}
}

template <typename ScalarA>
inline ScalarA fixInf(const ScalarA & result)
{
  ScalarA tmp = result;
  if (std::isinf(std::real(result))) { bool neg = std::signbit(result); tmp = (1.0e+50)*(neg?-1.0:1.0); }
  return tmp;
}

template <>
inline std::complex<double> fixInf(const std::complex<double> & result)
{
  bool negReal = std::signbit(std::real(result));
  bool negImag = std::signbit(std::imag(result));
  std::complex<double> tmp = result;
  if (std::isinf(std::real(result))) {tmp = std::complex<double>((1.0e+50)*(negReal?-1.0:1.0),result.imag());}
  if (std::isinf(std::imag(result))) {tmp = std::complex<double>(result.real(),(1.0e+50)*(negImag?-1.0:1.0));}
  return tmp;
}

}
}


inline void yyerror(std::vector<std::string> & s);

#include "ast_visitor.h"

//-------------------------------------------------------------------------------
// this is to make the call to "getStateOps" have a single
// function argument that never has to change.
template <typename ScalarT>
struct stateOpVectorContainers
{
public:
  stateOpVectorContainers(
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & sdt,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & ddt
      ):
    sdtOpVector(sdt),
    ddtOpVector(ddt)
  {};

  std::vector< Teuchos::RCP<astNode<ScalarT> > > & sdtOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & ddtOpVector;
};

struct staticsContainer
{
  static unsigned long int nextID;
  static unsigned long int stepNumber;
  static bool processSuccessfulStepFlag;
  static std::unordered_map<unsigned long int,int> processSuccessfulStepMap;
};

template <typename ScalarT>
struct sdtStateData : public staticsContainer
{
  sdtStateData(): id(0), val1(0.0), val2(0.0), integral_old(0.0), integral(0.0)
  { id = ++(this->nextID); };

  virtual void processSuccessfulTimeStep ()
  {
    integral_old = integral;
    val1 = val2;
  };

  unsigned long int id;
  ScalarT val1;
  ScalarT val2;
  ScalarT integral_old;
  ScalarT integral;
};

template <typename ScalarT>
struct ddtStateData : public staticsContainer
{
  ddtStateData(): id(0), val1(0.0), val2(0.0)
  { id = ++(this->nextID); };

  virtual void processSuccessfulTimeStep () { val1 = val2; };

  unsigned long int id;
  ScalarT val1;
  ScalarT val2;
};

//-------------------------------------------------------------------------------
// base node class
template <typename ScalarT>
class astNode : public staticsContainer
{
  public:
    astNode(): id_(0) { id_ = ++(this->nextID); };

    astNode( Teuchos::RCP<astNode<ScalarT> > &left ): leftAst_(left), id_(0) { id_ = ++(this->nextID); };

    astNode(Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):
      leftAst_(left),rightAst_(right) { id_ = ++(this->nextID); };


    virtual ~astNode() = default;
    astNode(const astNode&) = default;
    astNode& operator=(const astNode&) = default;
    astNode(astNode&&) = default;
    astNode& operator=(astNode&&) = default;

    virtual void processSuccessfulTimeStep ()
    {
      sdtState_.integral_old = sdtState_.integral;
      sdtState_.val1 = sdtState_.val2;
      ddtState_.val1 = ddtState_.val2;
    };

    virtual bool updateForStep() { return false; }

    virtual ScalarT val() = 0;
    virtual ScalarT dx(int i) = 0;

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs) = 0;

    virtual void output(std::ostream & os, int indent=0) = 0;
    virtual void compactOutput(std::ostream & os) = 0;
    virtual void codeGen (std::ostream & os )
    {
      os << "// This node has not implemented a code gen function yet" <<std::endl;
    }

    virtual void setNode(Teuchos::RCP<astNode<ScalarT> > & tmpNode) {};
    virtual void unsetNode() {};

    virtual void setFuncArgs(const std::vector< Teuchos::RCP<astNode<ScalarT> > > & tmpArgVec ) {};

    virtual bool getBreakPoints(std::vector<Xyce::Util::BreakPoint> & breakPointTimes) { return true;}
    virtual void setBreakPointTol(double tol){return;};
    virtual void setStartingTimeStep(double timeStep){return;};
    virtual void setFinalTime(double finalTime){return;};

    virtual void setDerivIndex(int i) {};
    virtual void unsetDerivIndex() {};

    virtual ScalarT getValue() { return 0.0; }
    virtual void setValue(ScalarT val) {}; // supports specialsOp, paramOp and globalParamLayerOp otherwise no-op
    virtual void unsetValue() {};          // supports specialsOp, paramOp and globalParamLayerOp otherwise no-op

    virtual bool getIsComplex () { return (typeid(ScalarT) == typeid(std::complex<double>)) ; }

    // base class no-ops.  Derived functions only in paramOp, base class version only called from ddx.
    virtual void setIsVar() {};
    virtual void unsetIsVar() {};
    virtual bool getIsVar() { return false; }

    // these 3 only apply to paramOp
    virtual void setIsConstant() {};
    virtual void unsetIsConstant() {};
    virtual bool getIsConstant() { return false; }

    // this applies to any node.  Use to determine if the AST (or part of it) is constant.
    virtual bool getIsTreeConstant() { return false; }

    // these 3 only apply to paramOp
    virtual void setIsAttached() {};
    virtual void unsetIsAttached() {};
    virtual bool getIsAttached() { return false; }

    // various "getType" functions.  There is a cleaner way to do this, but I haven't had the time
    virtual bool numvalType()      { return false; };
    virtual bool paramType()       { return false; };
    virtual bool funcType()        { return false; };
    virtual bool voltageType()     { return false; };
    virtual bool currentType()     { return false; };
    virtual bool leadCurrentType() { return false; };
    virtual bool bsrcCurrentType() { return false; };
    virtual bool powerType()       { return false; };
    virtual bool internalDeviceVarType()  { return false; };

    virtual bool dnoNoiseVarType() { return false; }
    virtual bool dniNoiseVarType() { return false; }
    virtual bool oNoiseType() { return false; }
    virtual bool iNoiseType() { return false; }

    virtual bool sdtType() { return false; }
    virtual bool ddtType() { return false; }
    virtual bool srcType() { return false; }
    virtual bool stpType() { return false; }
    virtual bool compType() { return false; }
    virtual bool limitType() { return false; }
    virtual bool phaseType()       { return false; };
    virtual bool sparamType()       { return false; };
    virtual bool yparamType()       { return false; };
    virtual bool zparamType()       { return false; };

// random op types
    virtual bool agaussType()      { return false; };
    virtual bool gaussType()       { return false; };
    virtual bool aunifType()       { return false; };
    virtual bool unifType()        { return false; };
    virtual bool randType()        { return false; };
    virtual bool twoArgLimitType() { return false; };

    virtual bool timeSpecialType() { return false; }
    virtual bool dtSpecialType()   { return false; }
    virtual bool tempSpecialType() { return false; }
    virtual bool vtSpecialType()   { return false; }
    virtual bool freqSpecialType() { return false; }
    virtual bool gminSpecialType() { return false; }

    virtual bool getFunctionArgType() { return false; };
    virtual void setFunctionArgType() {};
    virtual void unsetFunctionArgType() {};

    virtual bool scheduleType() { return false; }

    virtual std::string getName () { return std::string(""); };

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) {}

    virtual ddtStateData<ScalarT> & getDdtState() { return ddtState_; }
    virtual sdtStateData<ScalarT> & getSdtState() { return sdtState_; }

    virtual void setDdtState( const ddtStateData<ScalarT> & ddt) { ddtState_ = ddt; };
    virtual void setSdtState( const sdtStateData<ScalarT> & sdt) { sdtState_ = sdt; };

    unsigned long int getId () { return id_; }
    virtual unsigned long int getNodeId () { return id_; }

  protected:
    Teuchos::RCP<astNode<ScalarT> > leftAst_;
    Teuchos::RCP<astNode<ScalarT> > rightAst_;

    ddtStateData<ScalarT> ddtState_;
    sdtStateData<ScalarT> sdtState_;

    unsigned long int id_;
};

//-------------------------------------------------------------------------------
template <typename ScalarT>
class numval : public astNode<ScalarT>
{
  public:
    numval (ScalarT d): astNode<ScalarT>(),number(d) {};
    numval (std::complex<ScalarT> d): astNode<ScalarT>(),number(std::real(d)) {};

    virtual ScalarT val() { return number; }

    virtual ScalarT dx(int i) {return 0.0;}

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = number;
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0); }
      return;
    };

    virtual bool getIsComplex () { return false; }

    ScalarT number;

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "numval number = " << number
        << " id = " << this->id_ << std::endl;
    }

    virtual void compactOutput(std::ostream & os) { output(os,0); }

    virtual void codeGen (std::ostream & os )
    {
      // fix this later for formatting
      os << number;
    }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    { 
      Teuchos::RCP<numval<ScalarT> > castToThis = Teuchos::rcp_static_cast<numval<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
    } // 2nd dispatch

    virtual bool getIsTreeConstant() { return true; }
    virtual bool numvalType() { return true; };
};

//-------------------------------------------------------------------------------
template <>
class numval<std::complex<double>> : public astNode<std::complex<double>>
{
  public:

    numval (std::complex<double> d): astNode<std::complex<double> >(),number(d) {};

    virtual std::complex<double> val() {return number;}
    virtual std::complex<double> dx(int i) {return (std::complex<double>(0.0,0.0));}

    virtual void dx2(std::complex<double> & result, std::vector<std::complex<double> > & derivs)
    {
      result = number;
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),std::complex<double>(0.0,0.0)); }
      return;
    }

    virtual bool getIsComplex () { return (std::imag(number) != 0.0) ; }

    std::complex<double> number;

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "numval number = " << number
        << " id = " << id_ << std::endl;
    }

    virtual void compactOutput(std::ostream & os) { output(os,0); }

    virtual void codeGen (std::ostream & os )
    {
      // fix this later for formatting
      os << "std::complex<double>" << number;
    }

    virtual void accept (nodeVisitor<std::complex<double>> & visitor, Teuchos::RCP<astNode<std::complex<double>> > & thisAst_) 
    { 
      Teuchos::RCP<numval<std::complex<double>> > castToThis = Teuchos::rcp_static_cast<numval<std::complex<double>> > (thisAst_);
      visitor.visit( castToThis );
    } // 2nd dispatch

    virtual bool getIsTreeConstant() { return true; }
    virtual bool numvalType() { return true; };
};

#include "astbinary.h"


//-------------------------------------------------------------------------------
// This function is called by the various comparison operators (greater than,
// less than, etc) as well as the stpOp class to compute breakpoints.
// It is only intended to work for comparisons involving time.  It is not designed
// to catch other types of events, such as when a voltage exceeds a given value.
//
// This function uses Newton's method to compute the (hopefully) nearest
// breakpoint.  The tolerance is the breakpoint tolerance.  The max number
// of iterations is 20.
//
// This is typically called from the val() function of an operator.  The
// reason this is necessary is to handle the use case of the operator being
// inside of a .func.  If inside of a .func, then it has to be called as
// part of an AST traversal.  As a result, this function is called more than it
// needs to be  (every Newton step instead of every time step) and is a little
// bit wasteful in that regard.  Possibly, this could be revisited later to
// squeeze out more efficiency.
//-------------------------------------------------------------------------------
  template <typename ScalarT>
inline void computeBreakPoint(
    Teuchos::RCP<astNode<ScalarT> > & leftAst_,
    Teuchos::RCP<astNode<ScalarT> > & rightAst_,
    std::vector<Teuchos::RCP<astNode<ScalarT> > > & timeOpVec_,
    double bpTol_,
    std::vector<Xyce::Util::BreakPoint> & bpTimes_
    )
{
  timeOpVec_.clear();

  getTimeOpsVisitor<ScalarT> visitor(timeOpVec_);
  leftAst_->accept(visitor,leftAst_);
  rightAst_->accept(visitor,rightAst_);

  if (!(timeOpVec_.empty()))
  {
    // The following uses Newton's method to obtain the next breakpoint.
    // If it fails to converge, then it will not save one.
    // Possibly bpTol (which is the breakpoint tolerance used by the time integrator)
    // is too tight for this calculation.  Look into this later, perhaps.
    Teuchos::RCP<binaryMinusOp<ScalarT> > f_Ast_ = Teuchos::rcp(new binaryMinusOp<ScalarT>(leftAst_,rightAst_));

    int index = -99999;

    for (int ii=0; ii< timeOpVec_.size(); ii++) { timeOpVec_[ii]->setDerivIndex(index); }
    ScalarT f = f_Ast_->val();        Xyce::Util::fixNan(f);  Xyce::Util::fixInf(f);
    ScalarT dfdt = f_Ast_->dx(index); Xyce::Util::fixNan(dfdt);  Xyce::Util::fixInf(dfdt);

    // The Newton iterate is:  -(F(t)-A)/F'(t)
    double time = std::real(timeOpVec_[0]->val());
    double delta_bpTime =  0.0;

    // if initial dfdt==0, then algorithm has no hope of succeeding
    // if initial f==0, then we already at a BP, and we shouldn't set it again.
    if ( std::abs(std::real(f)) > bpTol_ && std::real(dfdt) != 0.0)
    {
      delta_bpTime =  -std::real(f)/std::real(dfdt);
      double bpTime = time+delta_bpTime;

      // test
      for (int ii=0; ii< timeOpVec_.size(); ii++) { timeOpVec_[ii]->setValue(bpTime); }
      f = f_Ast_->val();        Xyce::Util::fixNan(f);  Xyce::Util::fixInf(f);
      dfdt = f_Ast_->dx(index); Xyce::Util::fixNan(dfdt);  Xyce::Util::fixInf(dfdt);

      int iteration = 1;
      while (std::abs(std::real(f)) > bpTol_ && iteration < 20)  // iterate
      {
      delta_bpTime =  0.0;
      if (std::real(dfdt) != 0.0) { delta_bpTime =  -std::real(f)/std::real(dfdt); }
      bpTime +=delta_bpTime;

      for (int ii=0; ii< timeOpVec_.size(); ii++) { timeOpVec_[ii]->setValue(bpTime); }
      f = f_Ast_->val();        Xyce::Util::fixNan(f);  Xyce::Util::fixInf(f);
      dfdt = f_Ast_->dx(index); Xyce::Util::fixNan(dfdt);  Xyce::Util::fixInf(dfdt);

      ++iteration;
      }

      if (std::abs(std::real(f)) <= bpTol_) { bpTimes_.push_back( bpTime ); }// save this breakpoint if we converged
    }

    // unset the index, restore time
    for (int ii=0; ii< timeOpVec_.size(); ii++) { timeOpVec_[ii]->unsetDerivIndex(); }
    for (int ii=0; ii< timeOpVec_.size(); ii++) { timeOpVec_[ii]->setValue(time); }
  }
}

#include "astfuncs.h"
#include "astcomp.h"
#include "ast_spice_src.h"

//-------------------------------------------------------------------------------
// pow "^"  operator
template <typename ScalarT>
class powOp : public astNode<ScalarT>
{
  public:
    powOp (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):
      astNode<ScalarT>(left,right), rightConst_(true),leftConst_(false)
    {
      rightConst_ = this->rightAst_->numvalType();
      leftConst_ = this->leftAst_->numvalType();
    };

    virtual ScalarT val() { return std::pow(this->leftAst_->val(), this->rightAst_->val());}

    virtual ScalarT dx (int i)
    {
      Teuchos::RCP<astNode<ScalarT> > & lef = this->leftAst_;
      Teuchos::RCP<astNode<ScalarT> > & rig = this->rightAst_;
      ScalarT retVal = 0.0;
      ScalarT leftVal=lef->val();
      ScalarT righVal=rig->val();

      if (rightConst_ && !leftConst_) {if (leftVal != 0.0) { retVal = righVal*lef->dx(i)/leftVal*std::pow(leftVal,righVal) ; }}
      else if (!rightConst_ && leftConst_) {if (leftVal != 0.0) { retVal = std::log(leftVal)*std::pow(leftVal,righVal)*rig->dx(i); }}
      else {if (leftVal != 0.0) { retVal = (rig->dx(i)*std::log(leftVal)+righVal*lef->dx(i)/leftVal)*std::pow(leftVal,righVal);}}
      return  retVal;
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      Teuchos::RCP<astNode<ScalarT> > & lef = this->leftAst_;
      Teuchos::RCP<astNode<ScalarT> > & rig = this->rightAst_;

      std::vector<ScalarT> & retVal = derivs;

      int numDerivs = derivs.size();
      std::vector<ScalarT> lefDerivs_;
      std::vector<ScalarT> rigDerivs_;
      lefDerivs_.resize(numDerivs,0.0);
      rigDerivs_.resize(numDerivs,0.0);

      ScalarT leftVal, righVal;
      lef->dx2(leftVal,lefDerivs_);
      rig->dx2(righVal,rigDerivs_);
      result  = std::pow(leftVal, righVal);

      if (rightConst_ && !leftConst_)
      {
        if (leftVal != 0.0)
        {
          for (int ii=0;ii<numDerivs;ii++)
          {
            retVal[ii] = righVal*lefDerivs_[ii]/leftVal*std::pow(leftVal,righVal) ;
          }
        }
      }
      else if (!rightConst_ && leftConst_)
      {
        if (leftVal != 0.0)
        {
          for (int ii=0;ii<numDerivs;ii++)
          {
            retVal[ii] = std::log(leftVal)*std::pow(leftVal,righVal)*rigDerivs_[ii];
          }
        }
      }
      else
      {
        if (leftVal != 0.0)
        {
          for (int ii=0;ii<numDerivs;ii++)
          {
            retVal[ii] = (rigDerivs_[ii]*std::log(leftVal)+righVal*lefDerivs_[ii]/leftVal)*std::pow(leftVal,righVal);
          }
        }
      }
    }

    virtual bool getIsComplex () { return (this->rightAst_->getIsComplex() || this->leftAst_->getIsComplex()); }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "power operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "power operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::pow(";
      this->leftAst_->codeGen(os);
      os << ",";
      this->rightAst_->codeGen(os);
      os << ")";
    }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    {
      Teuchos::RCP<powOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<powOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
      this->leftAst_->accept(visitor, this->leftAst_); 
      this->rightAst_->accept(visitor, this->rightAst_); 
    }

    virtual bool getIsTreeConstant()
    { return (this->leftAst_->getIsTreeConstant() && this->rightAst_->getIsTreeConstant()); }

  private:
    bool rightConst_;
    bool leftConst_;
};

//-------------------------------------------------------------------------------
// atan2(y,x) operator
template <typename ScalarT>
class atan2Op : public astNode<ScalarT>
{
  public:
    atan2Op (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):
      astNode<ScalarT>(left,right), rightConst_(true),leftConst_(false)
    {
      rightConst_ = this->rightAst_->numvalType();
      leftConst_ = this->leftAst_->numvalType();
    };

    virtual ScalarT val() { return std::atan2(std::real(this->leftAst_->val()), std::real(this->rightAst_->val()));}

    virtual ScalarT dx (int i)
    {
      Teuchos::RCP<astNode<ScalarT> > & lef = this->leftAst_;
      Teuchos::RCP<astNode<ScalarT> > & rig = this->rightAst_;
      ScalarT retVal = 0.0;

      ScalarT leftVal=lef->val();
      ScalarT righVal=rig->val();

      if (rightConst_ && !leftConst_) { retVal = (righVal*lef->dx(i))/ (leftVal*leftVal + righVal*righVal); }
      else if (!rightConst_ && leftConst_) { retVal = (-leftVal*rig->dx(i)) / (leftVal*leftVal + righVal*righVal); }
      else { retVal = (righVal*lef->dx(i) - leftVal*rig->dx(i))/ (leftVal*leftVal + righVal*righVal) ; }
      return  retVal;
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      Teuchos::RCP<astNode<ScalarT> > & lef = this->leftAst_;
      Teuchos::RCP<astNode<ScalarT> > & rig = this->rightAst_;

      std::vector<ScalarT> & retVal = derivs;

      int numDerivs = derivs.size();
      std::vector<ScalarT> lefDerivs_;
      std::vector<ScalarT> rigDerivs_;
      lefDerivs_.resize(numDerivs,0.0);
      rigDerivs_.resize(numDerivs,0.0);

      ScalarT leftVal;
      ScalarT righVal;
      lef->dx2(leftVal,lefDerivs_);
      rig->dx2(righVal,rigDerivs_);

      result = std::atan2(std::real(leftVal), std::real(righVal));

      if (rightConst_ && !leftConst_)
      {
        for (int ii=0;ii<numDerivs;ii++)
        {
          retVal[ii] = (righVal* lefDerivs_[ii])/ (leftVal*leftVal + righVal*righVal);
        }
      }
      else if (!rightConst_ && leftConst_)
      {
        for (int ii=0;ii<numDerivs;ii++)
        {
          retVal[ii] = (-leftVal*rigDerivs_[ii]) / (leftVal*leftVal + righVal*righVal);
        }
      }
      else
      {
        for (int ii=0;ii<numDerivs;ii++)
        {
          retVal[ii] = (righVal*lefDerivs_[ii] - leftVal*rigDerivs_[ii])/ (leftVal*leftVal + righVal*righVal) ;
        }
      }
    }

    virtual bool getIsComplex () { return (this->rightAst_->getIsComplex() || this->leftAst_->getIsComplex()); }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "atan2 operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "atan2 operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::atan2(";
      this->leftAst_->codeGen(os);
      os << ",";
      this->rightAst_->codeGen(os);
      os << ")";
    }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    {
      Teuchos::RCP<atan2Op<ScalarT> > castToThis = Teuchos::rcp_static_cast<atan2Op<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
      this->leftAst_->accept(visitor, this->leftAst_); 
      this->rightAst_->accept(visitor, this->rightAst_); 
    }

    virtual bool getIsTreeConstant()
    { return (this->leftAst_->getIsTreeConstant() && this->rightAst_->getIsTreeConstant()); }

  private:
    bool rightConst_;
    bool leftConst_;
};

//-------------------------------------------------------------------------------
// phase operator
template <typename ScalarT>
class phaseOp : public astNode<ScalarT>
{
  public:
    phaseOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left), phaseOutputUsesRadians_(false)
    {};

    virtual ScalarT val()
    {
      return (std::arg(this->leftAst_->val())) * ((phaseOutputUsesRadians_)?1.0:(180.0/M_PI));
    }

    virtual ScalarT dx(int i)
    {
      std::vector<std::string> errStr(1,std::string("AST node (phase) is not differentiable"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();
      if ( !(derivs.empty() ) )
      {
        std::fill(derivs.begin(),derivs.end(),0.0);
        std::vector<std::string> errStr(1,std::string("AST node (phase) is not differentiable"));
        yyerror(errStr);
      }
    }

    virtual bool getIsComplex () { return false; }

    virtual bool phaseType()       { return true; };

    void setPhaseOutputUsesRadians( bool usesRadians )
    {
      phaseOutputUsesRadians_ = usesRadians;
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "phase operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "phase operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::arg(";
      this->leftAst_->codeGen(os);
      os << ")";
    }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    {
      Teuchos::RCP<phaseOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<phaseOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
      this->leftAst_->accept(visitor, this->leftAst_); 
    }

    virtual bool getIsTreeConstant() { return this->leftAst_->getIsTreeConstant() ; }

  private:
    bool phaseOutputUsesRadians_;
};

//-------------------------------------------------------------------------------
// real part  operator
template <typename ScalarT>
class realOp : public astNode<ScalarT>
{
  public:
    realOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val() { return std::real(this->leftAst_->val()); }

    virtual ScalarT dx(int i)
    {
      std::vector<std::string> errStr(1,std::string("AST node (real) is not differentiable"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();
      if ( !(derivs.empty() ) )
      {
        std::fill(derivs.begin(),derivs.end(),0.0);
        std::vector<std::string> errStr(1,std::string("AST node (real) is not differentiable"));
        yyerror(errStr);
      }
    }

    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "real operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "real operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::real(";
      this->leftAst_->codeGen(os);
      os << ")";
    }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    {
      Teuchos::RCP<realOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<realOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
      this->leftAst_->accept(visitor, this->leftAst_); 
    }

    virtual bool getIsTreeConstant() { return this->leftAst_->getIsTreeConstant() ; }
};

//-------------------------------------------------------------------------------
// imag part  operator
template <typename ScalarT>
class imagOp : public astNode<ScalarT>
{
  public:
    imagOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val() { return std::imag(this->leftAst_->val()); }

    virtual ScalarT dx(int i)
    {
      std::vector<std::string> errStr(1,std::string("AST node (imag) is not differentiable"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();
      if ( !(derivs.empty() ) )
      {
        std::fill(derivs.begin(),derivs.end(),0.0);
        std::vector<std::string> errStr(1,std::string("AST node (imag) is not differentiable"));
        yyerror(errStr);
      }
    }

    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "imag operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "imag operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::imag(";
      this->leftAst_->codeGen(os);
      os << ")";
    }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    {
      Teuchos::RCP<imagOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<imagOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
      this->leftAst_->accept(visitor, this->leftAst_); 
    }

    virtual bool getIsTreeConstant() { return this->leftAst_->getIsTreeConstant() ; }
};

//-------------------------------------------------------------------------------
// max operator
template <typename ScalarT>
class maxOp : public astNode<ScalarT>
{
  public:
    maxOp ( Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right ): astNode<ScalarT>(left,right) {};

    virtual ScalarT val()
    { return std::max( std::real(this->leftAst_->val()), std::real(this->rightAst_->val()) ); }

    virtual ScalarT dx(int i)
    {
      bool cmp = std::real(this->leftAst_->val()) < std::real(this->rightAst_->val());
      return cmp?(std::real(this->rightAst_->dx(i))):(std::real(this->leftAst_->dx(i)));
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      int numDerivs = derivs.size();
      std::vector<ScalarT> lefDerivs_;
      std::vector<ScalarT> rigDerivs_;
      lefDerivs_.resize(numDerivs,0.0);
      rigDerivs_.resize(numDerivs,0.0);

      ScalarT leftVal, rightVal;
      this->leftAst_->dx2(leftVal,lefDerivs_);
      this->rightAst_->dx2(rightVal,rigDerivs_);

      result = std::max( std::real(leftVal), std::real(rightVal));

      bool cmp = std::real(leftVal) < std::real(rightVal);
      for (int i=0;i<numDerivs;i++)
      {
        derivs[i] = cmp?(std::real(rigDerivs_[i])):(std::real(lefDerivs_[i]));
      }
    }

    virtual bool getIsComplex () { return false; } // this operator only uses the real part

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "max operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "max operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::max(";
      this->leftAst_->codeGen(os);
      os << ",";
      this->rightAst_->codeGen(os);
      os << ")";
    }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    {
      Teuchos::RCP<maxOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<maxOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
      this->leftAst_->accept(visitor, this->leftAst_); 
      this->rightAst_->accept(visitor, this->rightAst_); 
    }

    virtual bool getIsTreeConstant()
    { return (this->leftAst_->getIsTreeConstant() && this->rightAst_->getIsTreeConstant()); }
};

//-------------------------------------------------------------------------------
// min operator
template <typename ScalarT>
class minOp : public astNode<ScalarT>
{
  public:
    minOp ( Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right ): astNode<ScalarT>(left,right) {};

    virtual ScalarT val() { return std::min( std::real(this->leftAst_->val()), std::real(this->rightAst_->val()) ); }

    virtual ScalarT dx(int i)
    {
      bool cmp = std::real(this->rightAst_->val()) < std::real(this->leftAst_->val()) ;
      return (!cmp)?(std::real(this->leftAst_->dx(i))):(std::real(this->rightAst_->dx(i)));
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      int numDerivs = derivs.size();
      std::vector<ScalarT> lefDerivs_;
      std::vector<ScalarT> rigDerivs_;
      lefDerivs_.resize(numDerivs,0.0);
      rigDerivs_.resize(numDerivs,0.0);
      ScalarT leftVal, rightVal;
      this->leftAst_->dx2(leftVal,lefDerivs_);
      this->rightAst_->dx2(rightVal,rigDerivs_);

      result = std::min( std::real(leftVal), std::real(rightVal) );
      bool cmp = std::real(rightVal) < std::real(leftVal) ;
      for (int i=0;i<numDerivs;i++)
      {
        derivs[i] = (!cmp)?(std::real(lefDerivs_[i])):(std::real(rigDerivs_[i]));
      }
    }

    virtual bool getIsComplex () { return false; } // this operator only uses the real part

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "min operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "min operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::min(";
      this->leftAst_->codeGen(os);
      os << ",";
      this->rightAst_->codeGen(os);
      os << ")";
    }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    {
      Teuchos::RCP<minOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<minOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
      this->leftAst_->accept(visitor, this->leftAst_); 
      this->rightAst_->accept(visitor, this->rightAst_); 
    }

    virtual bool getIsTreeConstant()
    { return (this->leftAst_->getIsTreeConstant() && this->rightAst_->getIsTreeConstant()); }
};

//-------------------------------------------------------------------------------
// unary logical NOT sign
template <typename ScalarT>
class unaryNotOp : public astNode<ScalarT>
{
  public:
    unaryNotOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val()
    {
      return ((std::real(this->leftAst_->val())==0)?1:0);
    }

    virtual ScalarT dx(int i) {return 0.0;}

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = ((std::real(this->leftAst_->val())==0)?1:0);
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0); }
    }

    virtual bool getIsComplex () { return false; } // this operator only uses the real part

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "unary NOT operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "unary NOT operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "(!";
      this->leftAst_->codeGen(os);
      os << ")";
    }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    {
      Teuchos::RCP<unaryNotOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<unaryNotOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
      this->leftAst_->accept(visitor, this->leftAst_); 
    }

    virtual bool getIsTreeConstant() { return (this->leftAst_->getIsTreeConstant()); }
};

//-------------------------------------------------------------------------------
// unary minus sign
template <typename ScalarT>
class unaryMinusOp : public astNode<ScalarT>
{
  public:
    unaryMinusOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val() { return (-(this->leftAst_->val())); }
    virtual ScalarT dx(int i) { return (-(this->leftAst_->dx(i))); }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      int numDerivs = derivs.size();
      std::vector<ScalarT> lefDerivs_;
      lefDerivs_.resize(numDerivs,0.0);
      this->leftAst_->dx2(result,lefDerivs_);
      result *=  -1.0;
      for (int i=0;i<numDerivs;i++) { derivs[i] = (-(lefDerivs_[i])); }
    }

    virtual bool getIsComplex () { return this->leftAst_->getIsComplex(); }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "unary minus operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "unary minus operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "(-";
      this->leftAst_->codeGen(os);
      os << ")";
    }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    {
      Teuchos::RCP<unaryMinusOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<unaryMinusOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
      this->leftAst_->accept(visitor, this->leftAst_); 
    }

    virtual bool getIsTreeConstant() { return (this->leftAst_->getIsTreeConstant()); }
    virtual bool numvalType() { return (this->leftAst_->numvalType()); };
};

//-------------------------------------------------------------------------------
// unary plus sign
template <typename ScalarT>
class unaryPlusOp : public astNode<ScalarT>
{
  public:
    unaryPlusOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val() { return (+(this->leftAst_->val())); }

    virtual ScalarT dx(int i) { return (+(this->leftAst_->dx(i))); }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      this->leftAst_->dx2(result,derivs);
    }

    virtual bool getIsComplex () { return this->leftAst_->getIsComplex(); }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "unary plus operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "unary plus operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "(+";
      this->leftAst_->codeGen(os);
      os << ")";
    }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    {
      Teuchos::RCP<unaryPlusOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<unaryPlusOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
      this->leftAst_->accept(visitor, this->leftAst_); 
    }

    virtual bool getIsTreeConstant() { return (this->leftAst_->getIsTreeConstant()); }
    virtual bool numvalType() { return (this->leftAst_->numvalType()); };
};


//-------------------------------------------------------------------------------
// This class is designed to help process global parameters.
//
// It is part of a refactor to get rid of the "make_var" function in the
// newExpression class, which is used to set certain parameters as "variable".
// These were generally .global_params that were likely to be reset by a .STEP or
// .SAMPLING, or some other analysis that modifies params.
//
// This is used on expressions that are the RHS of a global_param statement.
// If an expression is a global_param, then it needs to optionally be
// replaced with a value.  This class makes it possible for that to happen.
//-------------------------------------------------------------------------------
template <typename ScalarT>
class globalParamLayerOp: public astNode<ScalarT>
{
  public:
    globalParamLayerOp ():
      astNode<ScalarT>()
    {
      numvalNode_ = Teuchos::rcp(new numval<ScalarT> (0.0));
      paramNode_ = numvalNode_;
      savedParamNode_ = numvalNode_;
    };

    virtual ScalarT val() { return paramNode_->val(); }
    virtual ScalarT dx(int i)
    {
      ScalarT retval=0.0;
      retval = paramNode_->dx(i);
      return retval;
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      paramNode_->dx2(result,derivs);
    }

    virtual bool getIsComplex () { return paramNode_->getIsComplex(); }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "globalParamLayer Op  val = " << val()
        << " id = " << this->id_
        << " node_id = " << paramNode_->getId()
        << std::endl;
        paramNode_->output(os,indent+2);
    }

    virtual void compactOutput(std::ostream & os)
    {
       os << "globalParamLayer Op  val = " << val() << " id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os ) { }

    virtual void setNode(Teuchos::RCP<astNode<ScalarT> > & tmpNode) { paramNode_ = tmpNode; savedParamNode_ = tmpNode; };
    virtual void unsetNode() { paramNode_ = numvalNode_; };

    virtual ScalarT getValue() { return numvalNode_->number; };
    virtual void setValue(ScalarT val) { numvalNode_->number = val; paramNode_ = numvalNode_; };
    virtual void unsetValue() { paramNode_ = savedParamNode_; };

    virtual void processSuccessfulTimeStep ()
    {
      paramNode_->processSuccessfulTimeStep ();
    };

    virtual bool getIsTreeConstant() { return (paramNode_->getIsTreeConstant()); }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    {
      Teuchos::RCP<globalParamLayerOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<globalParamLayerOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
      paramNode_->accept(visitor, paramNode_); 
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > paramNode_;
    Teuchos::RCP<astNode<ScalarT> > savedParamNode_;
    Teuchos::RCP<numval<ScalarT> > numvalNode_;
};

//-------------------------------------------------------------------------------
//
// This is the parameter Op class.
//
// This is still a work in progress.  I originally wrote it to primarily behave
// like a global_param.  ie, something that could dynamically change, and was
// attached to an external expression via the "paramNode" object, which points
// to the top of the ast tree of another expression.
//
// I have been modifying this slowly to make it so that it can also (under
// some circumstances) behave like a Xyce .param.  ie, a hardwired constant.
// This aspect is not 100% fleshed out yet.
//
// For a simple evaluation of the "val" function, this class will:
//
//    (1) call the "val" function of the underlying external syntax tree
//    (2) return the scalar quantity "number_"
//
// For derivative calculation, there is at least one use case that mixes these two
// modes of operation together.  So I still need to think about that.  Possibly
// there should be two different kinds of classes for this, to avoid confusion.
//
template <typename ScalarT>
class paramOp: public astNode<ScalarT>
{
  public:
    paramOp (const std::string & par):
      astNode<ScalarT>(),
      paramName_(par),
      thisIsAFunctionArgument_(false),
      isVar_(false),
      isConstant_(false),
      isAttached_(false),
      paramType_(DOT_GLOBAL_PARAM),
      derivIndex_(-1)
    {
      numvalNode_ = Teuchos::rcp(new numval<ScalarT> (0.0));
      paramNode_ = numvalNode_;
      savedParamNode_ = numvalNode_;
    };

    virtual ScalarT val() { return paramNode_->val(); }
    virtual ScalarT dx(int i)
    {
      ScalarT retval=0.0;
      if (isVar_) { retval = (derivIndex_==i)?1.0:0.0; }
      else        { retval = paramNode_->dx(i); }
      return retval;
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      if (isVar_)
      {
        result = paramNode_->val();
        if ( !(derivs.empty() ) )
        {
          std::fill(derivs.begin(),derivs.end(),0.0);
          if(derivIndex_>-1) { derivs[derivIndex_] = 1.0; }
        }
      }
      else { paramNode_->dx2(result,derivs); }
    }

    virtual bool getIsComplex () { return paramNode_->getIsComplex(); }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "parameter : " << paramName_ << " = " << val()
        << " id = " << this->id_
        << " node_id = " << paramNode_->getId()
        << std::endl;
        paramNode_->output(os,indent+2);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "parameter : " << paramName_ << " = " << val() << " id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << paramName_;
    }

    virtual void setNode(Teuchos::RCP<astNode<ScalarT> > & tmpNode) { paramNode_ = tmpNode; savedParamNode_ = tmpNode; };
    virtual void unsetNode() { paramNode_ = numvalNode_; };

    virtual ScalarT getValue() { return numvalNode_->number; };
    virtual void setValue(ScalarT val) { numvalNode_->number = val; paramNode_ = numvalNode_; };
    virtual void unsetValue() { paramNode_ = savedParamNode_; };

    virtual void setDerivIndex(int i) { derivIndex_=i; };
    virtual void unsetDerivIndex() {derivIndex_=-1;};

    void setName(const std::string & name) { paramName_ = name; }
    virtual std::string getName() { return paramName_; }

    virtual bool paramType() { return !thisIsAFunctionArgument_ ; };

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    {
      Teuchos::RCP<paramOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<paramOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
      paramNode_->accept(visitor, paramNode_); 
    }

    virtual bool getFunctionArgType() { return thisIsAFunctionArgument_; };
    virtual void setFunctionArgType() { thisIsAFunctionArgument_ = true;};
    virtual void unsetFunctionArgType() { thisIsAFunctionArgument_ = true;};

    // isAttached_, and isConstant_ are all checked by the
    // Expression::getUnresolvedParams function.
    //
    // isVar used to mean either:
    // (1) this is a global param of Util::DBLE type
    // or:
    // (2) we want the derivatives w.r.t. to this parameter.
    // Now it ONLY means (2).
    void setIsVar() { isVar_ = true; }
    void unsetIsVar() { isVar_ = false; }
    bool getIsVar() { return isVar_; }

    // this flag checks if this parameter is just a simple constant.
    // This is specific to the parameter operator and is different from "getIsTreeConstant"
    // Thus, it is not a derived function.  The class much be casted to paramOp for this to work.
    void setIsConstant() { isConstant_ = true; }
    void unsetIsConstant() { isConstant_ = false; }
    bool getIsConstant() { return isConstant_; }

    // getIsTreeConstant checks if the subordinate tree represents a constant expression.
    // If this parameter is a "global", meaning a parameter that is allowed to change via .STEP/.DC/.SAMPLING,
    // then it isn't constant and returns false.
    // If not a global, it might be constant, or it might not, depending on next node in tree.
    bool getIsTreeConstant()
    {
      if ( paramType_ == DOT_GLOBAL_PARAM ) { return false; }
      else { return paramNode_->getIsTreeConstant(); }
    }

    // this flag indicates if an external AST has been attached
    // to this class
    void setIsAttached() { isAttached_ = true; }
    void unsetIsAttached() { isAttached_ = false; }
    bool getIsAttached() { return isAttached_; }

    // the param type can be .param, .global_param or a subcircuit argument
    // The enum is defined as enum enumParamType {DOT_PARAM, DOT_GLOBAL_PARAM, SUBCKT_ARG_PARAM}
    void setParamType(enumParamType type) { paramType_ = type; }
    enumParamType getParamType() { return paramType_; }

    virtual void processSuccessfulTimeStep ()
    {
      paramNode_->processSuccessfulTimeStep ();
    };

    ddtStateData<ScalarT> & getDdtState() { return paramNode_->getDdtState(); }
    sdtStateData<ScalarT> & getSdtState() { return paramNode_->getSdtState(); }

    virtual unsigned long int getNodeId () { return paramNode_->getNodeId(); }

  private:
    // data:
    std::string paramName_;

    Teuchos::RCP<astNode<ScalarT> > paramNode_;
    Teuchos::RCP<astNode<ScalarT> > savedParamNode_;
    Teuchos::RCP<numval<ScalarT> > numvalNode_;

    bool thisIsAFunctionArgument_;

    bool isVar_;
    bool isConstant_;
    bool isAttached_;

    enumParamType paramType_;

    int derivIndex_;
};

//-------------------------------------------------------------------------------
template <typename ScalarT>
class voltageOp: public astNode<ScalarT>
{
  public:
    voltageOp (const std::string & voltageNode):
      astNode<ScalarT>(),
      voltageNode_(voltageNode),
      voltageVal_(0.0),
      derivIndex_(-1)
    {
      Xyce::Util::toUpper(voltageNode_);
    };

    virtual ScalarT val()
    {
      return voltageVal_;
    }

    virtual ScalarT dx(int i)
    {
      ScalarT retval = (derivIndex_==i)?1.0:0.0;
      return retval;
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = voltageVal_;
      if ( !(derivs.empty() ) )
      {
        std::fill(derivs.begin(),derivs.end(),0.0);
        if(derivIndex_>-1) { derivs[derivIndex_] = 1.0; }
      }
    }

    virtual bool getIsComplex () { return (typeid(ScalarT) == typeid(std::complex<double>)) ; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "Voltage node:" << " id = " << this->id_ << std::endl;

      os << std::setw(indent) << " ";
      os << "V(" << voltageNode_ << ") = " << voltageVal_ <<std::endl;

      os << std::setw(indent) << " ";
      os << "value = " << val() << std::endl;
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "Voltage node: id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "V_";
      os << voltageNode_;
    }

    virtual void setVal(const ScalarT & val)
    {
      voltageVal_ = val;
    }

    virtual void setDerivIndex(int i) { derivIndex_=i; };
    virtual void unsetDerivIndex() { derivIndex_=-1; };

    const std::string & getVoltageNode() { return voltageNode_; }
    void setVoltageNode(const std::string & name) { voltageNode_ = name; }
    ScalarT & getVoltageVal() { return voltageVal_; }

    virtual bool voltageType() { return true; };

    virtual std::string getName () { return voltageNode_; }

    virtual bool getIsTreeConstant() { return false; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<voltageOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<voltageOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
    } // 2nd dispatch

  private:
    std::string voltageNode_;
    ScalarT voltageVal_;
    int derivIndex_;
    std::string name_;
};

//-------------------------------------------------------------------------------
template <typename ScalarT>
class currentOp: public astNode<ScalarT>
{
  public:
    currentOp (const std::string & currentDevice):
      astNode<ScalarT>(),
      number_(0.0),
      currentDevice_(currentDevice),
      derivIndex_(-1),
      bsrcFlag_(false)
    {
      Xyce::Util::toUpper(currentDevice_);
    };

    virtual ScalarT val() {return number_;}

    virtual ScalarT dx(int i) { return (derivIndex_==i)?1.0:0.0; }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = number_;
      if ( !(derivs.empty() ) )
      {
        std::fill(derivs.begin(),derivs.end(),0.0);
        if (derivIndex_>-1) derivs[derivIndex_] = 1.0;
      }
    }

    virtual bool getIsComplex () { return (typeid(ScalarT) == typeid(std::complex<double>)) ; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "Isrc : device = " << currentDevice_ << " id = " << this->id_ << std::endl;
      os << std::setw(indent) << " " << "value = " << val() <<std::endl;
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "Isrc : device = " << currentDevice_ << " id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "I_";
      os << currentDevice_;
    }

    virtual void setDerivIndex(int i) {derivIndex_=i;};
    virtual void unsetDerivIndex() {derivIndex_=-1;};

    void setCurrentDevice(const std::string & devName) { currentDevice_ = devName; }
    const std::string & getCurrentDevice() { return currentDevice_; }

    ScalarT & getCurrentVal () { return number_; }
    void setCurrentVal (ScalarT n) { number_ = n; }

    virtual bool currentType() { return true; };
    virtual bool bsrcCurrentType() { return bsrcFlag_; }

    virtual std::string getName () { return currentDevice_; }

    bool getBsrcFlag   () { return bsrcFlag_; }
    void setBsrcFlag  () { bsrcFlag_ = true; }
    void unsetBsrcFlag  () { bsrcFlag_ = false; }

    virtual bool getIsTreeConstant() { return false; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<currentOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<currentOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
    } // 2nd dispatch

  private:
    ScalarT number_;
    std::string currentDevice_;
    int derivIndex_;
    bool bsrcFlag_;
};

//-------------------------------------------------------------------------------
template <typename ScalarT>
class sparamOp: public astNode<ScalarT>
{
  public:
    sparamOp (const std::vector<int> & args):
      astNode<ScalarT>(),
      number_(0.0),
      derivIndex_(-1),
      sparamArgs_(args)
    {
    };

    virtual ScalarT val() {return number_;}

    virtual ScalarT dx(int i) { return (derivIndex_==i)?1.0:0.0; }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = number_;
      if ( !(derivs.empty() ) )
      {
        std::fill(derivs.begin(),derivs.end(),0.0);
        if (derivIndex_>-1) derivs[derivIndex_] = 1.0;
      }
    }

    virtual bool getIsComplex () { return (typeid(ScalarT) == typeid(std::complex<double>)) ; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "SParam(";
      int size=sparamArgs_.size();
      for (int ii=0;ii<size;ii++)
      {
        os << sparamArgs_[ii];
        if (size>1 && ii < size-1) { os << ","; }
      }
      os << " id = " << this->id_ << std::endl;
      os << std::setw(indent) << " " << "value = " << val() <<std::endl;
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "SParam id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "S";
      int size=sparamArgs_.size();
      for (int ii=0;ii<size;ii++)
      {
        os << sparamArgs_[ii];
        if (size>1 && ii < size-1) { os << ","; }
      }
      os << " = " << val();
    }

    virtual void setDerivIndex(int i) {derivIndex_=i;};
    virtual void unsetDerivIndex() {derivIndex_=-1;};

    ScalarT & getSparamValue () { return number_; }
    virtual void setValue(ScalarT val) { number_ = val; };

    std::vector<int> & getSparamArgs () { return sparamArgs_; }

    virtual bool getIsTreeConstant() { return false; }
    virtual bool sparamType() { return true; };

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<sparamOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<sparamOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
    } // 2nd dispatch

  private:
// data:
    ScalarT number_;
    int derivIndex_;
    std::vector<int> sparamArgs_;
};

//-------------------------------------------------------------------------------
template <typename ScalarT>
class yparamOp: public astNode<ScalarT>
{
  public:
    yparamOp (const std::vector<int> & args):
      astNode<ScalarT>(),
      number_(0.0),
      derivIndex_(-1),
      yparamArgs_(args)
    {
    };

    virtual ScalarT val() {return number_;}

    virtual ScalarT dx(int i) { return (derivIndex_==i)?1.0:0.0; }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = number_;
      if ( !(derivs.empty() ) )
      {
        std::fill(derivs.begin(),derivs.end(),0.0);
        if (derivIndex_>-1) derivs[derivIndex_] = 1.0;
      }
    }

    virtual bool getIsComplex () { return (typeid(ScalarT) == typeid(std::complex<double>)) ; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "YParam(";
      int size=yparamArgs_.size();
      for (int ii=0;ii<size;ii++)
      {
        os << yparamArgs_[ii];
        if (size>1 && ii < size-1) { os << ","; }
      }
      os << " id = " << this->id_ << std::endl;
      os << std::setw(indent) << " " << "value = " << val() <<std::endl;
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "YParam id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "Y";
      int size=yparamArgs_.size();
      for (int ii=0;ii<size;ii++)
      {
        os << yparamArgs_[ii];
        if (size>1 && ii < size-1) { os << ","; }
      }
      os << " = " << val();
    }

    virtual void setDerivIndex(int i) {derivIndex_=i;};
    virtual void unsetDerivIndex() {derivIndex_=-1;};

    ScalarT & getYparamValue () { return number_; }
    virtual void setValue(ScalarT val) { number_ = val; };

    std::vector<int> & getYparamArgs () { return yparamArgs_; }

    virtual bool getIsTreeConstant() { return false; }
    virtual bool yparamType() { return true; };

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<yparamOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<yparamOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
    } // 2nd dispatch

  private:
// data:
    ScalarT number_;
    int derivIndex_;
    std::vector<int> yparamArgs_;
};

//-------------------------------------------------------------------------------
template <typename ScalarT>
class zparamOp: public astNode<ScalarT>
{
  public:
    zparamOp (const std::vector<int> & args):
      astNode<ScalarT>(),
      number_(0.0),
      derivIndex_(-1),
      zparamArgs_(args)
    {
    };

    virtual ScalarT val() {return number_;}

    virtual ScalarT dx(int i) { return (derivIndex_==i)?1.0:0.0; }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = number_;
      if ( !(derivs.empty() ) )
      {
        std::fill(derivs.begin(),derivs.end(),0.0);
        if (derivIndex_>-1) derivs[derivIndex_] = 1.0;
      }
    }

    virtual bool getIsComplex () { return (typeid(ScalarT) == typeid(std::complex<double>)) ; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "ZParam(";
      int size=zparamArgs_.size();
      for (int ii=0;ii<size;ii++)
      {
        os << zparamArgs_[ii];
        if (size>1 && ii < size-1) { os << ","; }
      }
      os << " id = " << this->id_ << std::endl;
      os << std::setw(indent) << " " << "value = " << val() <<std::endl;
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "ZParam id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "Z";
      int size=zparamArgs_.size();
      for (int ii=0;ii<size;ii++)
      {
        os << zparamArgs_[ii];
        if (size>1 && ii < size-1) { os << ","; }
      }
      os << " = " << val();
    }

    virtual void setDerivIndex(int i) {derivIndex_=i;};
    virtual void unsetDerivIndex() {derivIndex_=-1;};

    ScalarT & getZparamValue () { return number_; }
    virtual void setValue(ScalarT val) { number_ = val; };

    std::vector<int> & getZparamArgs () { return zparamArgs_; }

    virtual bool getIsTreeConstant() { return false; }
    virtual bool zparamType() { return true; };

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<zparamOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<zparamOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
    } // 2nd dispatch

  private:
// data:
    ScalarT number_;
    int derivIndex_;
    std::vector<int> zparamArgs_;
};

//-------------------------------------------------------------------------------
template <typename ScalarT>
class leadCurrentOp: public astNode<ScalarT>
{
  public:
    leadCurrentOp (const std::string & designator, const std::string & leadCurrentDevice):
      astNode<ScalarT>(),
      number_(0.0),
      leadCurrentDesignator_(designator),
      leadCurrentDevice_(leadCurrentDevice),
      derivIndex_(-1)
    {
      Xyce::Util::toUpper(leadCurrentDesignator_);
      Xyce::Util::toUpper(leadCurrentDevice_);
    };

    virtual ScalarT val() {return number_;}

    virtual ScalarT dx(int i) { return (derivIndex_==i)?1.0:0.0; }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = number_;
      if ( !(derivs.empty() ) )
      {
        std::fill(derivs.begin(),derivs.end(),0.0);
        if (derivIndex_>-1) derivs[derivIndex_] = 1.0;
      }
    }

    virtual bool getIsComplex () { return (typeid(ScalarT) == typeid(std::complex<double>)) ; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "Lead Current: device = " << leadCurrentDevice_ << " designator = " << leadCurrentDesignator_ << " id = " << this->id_ << std::endl;
      os << std::setw(indent) << " " << "value = " << val() <<std::endl;
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "Lead Current: device = " << leadCurrentDevice_ << " designator = " << leadCurrentDesignator_ << " id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "I_";
      os << leadCurrentDevice_;
    }

    virtual void setDerivIndex(int i) {derivIndex_=i;};
    virtual void unsetDerivIndex() {derivIndex_=-1;};

    void setLeadCurrentDevice(const std::string & devName) { leadCurrentDevice_ = devName; }
    const std::string & getLeadCurrentDevice() { return leadCurrentDevice_; }

    void setLeadCurrentDesignator(const std::string & desName) { leadCurrentDesignator_ = desName; }
    const std::string & getLeadCurrentDesignator() { return leadCurrentDesignator_; }

    ScalarT & getLeadCurrentVar () { return number_; }
    void setLeadCurrentVar (ScalarT n) { number_ = n; }

    virtual bool getIsTreeConstant() { return false; }
    virtual bool leadCurrentType() { return true; };

    virtual std::string getName () { return leadCurrentDevice_; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<leadCurrentOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<leadCurrentOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
    } // 2nd dispatch

  private:
// data:
    ScalarT number_;
    std::string leadCurrentDesignator_;
    std::string leadCurrentDevice_;
    int derivIndex_;
};

//-------------------------------------------------------------------------------
template <typename ScalarT>
class powerOp: public astNode<ScalarT>
{
  public:
    powerOp (const std::string & tag, const std::string & powerDevice):
      astNode<ScalarT>(),
      number_(0.0),
      tag_(tag),
      powerDevice_(powerDevice),
      derivIndex_(-1)
    {
      Xyce::Util::toUpper(powerDevice_);
    };

    virtual ScalarT val() {return number_;}

    virtual ScalarT dx(int i) { return (derivIndex_==i)?1.0:0.0; }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = number_;
      if ( !(derivs.empty() ) )
      {
        std::fill(derivs.begin(),derivs.end(),0.0);
        if (derivIndex_>-1) derivs[derivIndex_] = 1.0;
      }
    }

    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "Power : device = " << powerDevice_ << " id = " << this->id_ << std::endl;
      os << std::setw(indent) << " " << "value = " << val() <<std::endl;
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "Power : device = " << powerDevice_ << " id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "P_";
      os << powerDevice_;
    }

    virtual void setDerivIndex(int i) {derivIndex_=i;};
    virtual void unsetDerivIndex() {derivIndex_=-1;};

    void setPowerDevice(const std::string & devName) { powerDevice_ = devName; }
    const std::string & getPowerTag   () { return tag_; }
    const std::string & getPowerDevice() { return powerDevice_; }
    ScalarT & getPowerVal () { return number_; }
    void setPowerVal (ScalarT n) { number_ = n; }

    virtual bool getIsTreeConstant() { return false; }
    virtual bool powerType() { return true; };

    virtual std::string getName () { return powerDevice_; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<powerOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<powerOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
    } // 2nd dispatch

  private:
// data:
    ScalarT number_;
    std::string tag_; // either P or W
    std::string powerDevice_;
    int derivIndex_;
};

//-------------------------------------------------------------------------------
template <typename ScalarT>
class internalDevVarOp: public astNode<ScalarT>
{
  public:
    internalDevVarOp (const std::string & internalDevVarDevice):
      astNode<ScalarT>(),
      number_(0.0),
      internalDevVarDevice_(internalDevVarDevice),
      derivIndex_(-1)
    {
      Xyce::Util::toUpper(internalDevVarDevice_);
    };

    virtual ScalarT val() {return number_;}

    virtual ScalarT dx(int i) { return (derivIndex_==i)?1.0:0.0; }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = number_;
      if ( !(derivs.empty() ) )
      {
        std::fill(derivs.begin(),derivs.end(),0.0);
        if (derivIndex_>-1) derivs[derivIndex_] = 1.0;
      }
    }

    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "Internal device variable : device = " << internalDevVarDevice_ << " id = " << this->id_ << std::endl;
      os << std::setw(indent) << " " << "value = " << val() <<std::endl;
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "Internal device variable : device = " << internalDevVarDevice_ << " id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "I_";
      os << internalDevVarDevice_;
    }

    virtual void setDerivIndex(int i) {derivIndex_=i;};
    virtual void unsetDerivIndex() {derivIndex_=-1;};

    void setInternalVarDevice(const std::string & devName) { internalDevVarDevice_ = devName; }
    const std::string & getInternalVarDevice() { return internalDevVarDevice_; }
    ScalarT & getInternalDeviceVar () { return number_; }
    void setInternalDeviceVar (ScalarT n) { number_ = n; }

    virtual bool getIsTreeConstant() { return false; }
    virtual bool internalDeviceVarType()  { return true; };

    virtual std::string getName () { return internalDevVarDevice_; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<internalDevVarOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<internalDevVarOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
    } // 2nd dispatch

  private:
// data:
    ScalarT number_;
    std::string internalDevVarDevice_;
    int derivIndex_;
};

//-------------------------------------------------------------------------------
// noise ops.  Possibly combine these
//-------------------------------------------------------------------------------
template <typename ScalarT>
class dnoNoiseVarOp: public astNode<ScalarT>
{
  public:
    dnoNoiseVarOp (const std::vector<std::string> & noiseDevices):
      astNode<ScalarT>(),
      number_(0.0),
      noiseDevices_(noiseDevices),
      derivIndex_(-1)
    {
      for(int ii=0;ii<noiseDevices.size();ii++) { Xyce::Util::toUpper(noiseDevices_[ii]); }
    };

    virtual ScalarT val() {return number_;}

    virtual ScalarT dx(int i) { return (derivIndex_==i)?1.0:0.0; }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = number_;
      if ( !(derivs.empty() ) )
      {
        std::fill(derivs.begin(),derivs.end(),0.0);
        if (derivIndex_>-1) derivs[derivIndex_] = 1.0;
      }
    }

    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "DNO noise variable : devices = ";
      for (int ii=0;ii<noiseDevices_.size();ii++)
      {
        os << noiseDevices_[ii] << " ";
      }
      os << " id = " << this->id_ << std::endl;
      os << std::setw(indent) << " " << "value = " << val() <<std::endl;
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "DNO noise variable : devices = ";
      for (int ii=0;ii<noiseDevices_.size();ii++) { os << noiseDevices_[ii] << " "; }
      os << " id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "DNO_";
      for (int ii=0;ii<noiseDevices_.size();ii++)
      {
        os << "_" << noiseDevices_[ii];
      }
    }

    virtual void setDerivIndex(int i) {derivIndex_=i;};
    virtual void unsetDerivIndex() {derivIndex_=-1;};

    void setNoiseDevices (const std::vector<std::string> & devNames) { noiseDevices_ = devNames; }
    std::vector<std::string> getNoiseDevices () { return noiseDevices_; }

    ScalarT & getNoiseVar () { return number_; }
    void setNoiseVar (ScalarT n) { number_ = n; }

    virtual bool getIsTreeConstant() { return false; }
    virtual bool dnoNoiseVarType()  { return true; };

    //virtual std::string getName () { return noiseDevice_; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<dnoNoiseVarOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<dnoNoiseVarOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
    } // 2nd dispatch

  private:
    ScalarT number_;
    std::vector<std::string> noiseDevices_;
    int derivIndex_;
};

//-------------------------------------------------------------------------------
template <typename ScalarT>
class dniNoiseVarOp: public astNode<ScalarT>
{
  public:
    dniNoiseVarOp (const std::vector<std::string> & noiseDevices):
      astNode<ScalarT>(),
      number_(0.0),
      noiseDevices_(noiseDevices),
      derivIndex_(-1)
    {
      for(int ii=0;ii<noiseDevices.size();ii++) { Xyce::Util::toUpper(noiseDevices_[ii]); }
    };

    virtual ScalarT val() {return number_;}

    virtual ScalarT dx(int i) { return (derivIndex_==i)?1.0:0.0; }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = number_;
      if ( !(derivs.empty() ) )
      {
        std::fill(derivs.begin(),derivs.end(),0.0);
        if (derivIndex_>-1) derivs[derivIndex_] = 1.0;
      }
    }

    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "DNI noise variable : devices = ";
      for (int ii=0;ii<noiseDevices_.size();ii++)
      {
        os << noiseDevices_[ii] << " ";
      }
      os << " id = " << this->id_ << std::endl;
      os << std::setw(indent) << " " << "value = " << val() <<std::endl;
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "DNI noise variable : devices = ";
      for (int ii=0;ii<noiseDevices_.size();ii++) { os << noiseDevices_[ii] << " "; }
      os << " id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "DNI";
      for (int ii=0;ii<noiseDevices_.size();ii++)
      {
        os << "_" << noiseDevices_[ii];
      }
    }

    virtual void setDerivIndex(int i) {derivIndex_=i;};
    virtual void unsetDerivIndex() {derivIndex_=-1;};

    void setNoiseDevices (const std::vector<std::string> & devNames) { noiseDevices_ = devNames; }
    std::vector<std::string> getNoiseDevices () { return noiseDevices_; }
    ScalarT & getNoiseVar () { return number_; }
    void setNoiseVar (ScalarT n) { number_ = n; }

    virtual bool getIsTreeConstant() { return false; }
    virtual bool dniNoiseVarType()  { return true; };

    //virtual std::string getName () { return noiseDevice_; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<dniNoiseVarOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<dniNoiseVarOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
    } // 2nd dispatch

  private:
    ScalarT number_;
    std::vector<std::string> noiseDevices_;
    int derivIndex_;
};

//-------------------------------------------------------------------------------
template <typename ScalarT>
class oNoiseOp: public astNode<ScalarT>
{
  public:
    oNoiseOp (): astNode<ScalarT>(), number_(0.0), derivIndex_(-1) {};
    virtual ScalarT val() {return number_;}
    virtual ScalarT dx(int i) { return (derivIndex_==i)?1.0:0.0; }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = number_;
      if ( !(derivs.empty() ) )
      {
        std::fill(derivs.begin(),derivs.end(),0.0);
        if (derivIndex_>-1) derivs[derivIndex_] = 1.0;
      }
    }

    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "onoise variable id = " << this->id_ << std::endl;
      os << std::setw(indent) << " " << "value = " << val() <<std::endl;
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "onoise variable id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os ) { os << "ONOISE"; }

    virtual void setDerivIndex(int i) {derivIndex_=i;};
    virtual void unsetDerivIndex() {derivIndex_=-1;};
    ScalarT & getNoiseVar () { return number_; }
    void setNoiseVar (ScalarT n) { number_ = n; }
    virtual bool getIsTreeConstant() { return false; }
    virtual bool oNoiseType()  { return true; };

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<oNoiseOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<oNoiseOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
    } // 2nd dispatch

  private:
    ScalarT number_;
    int derivIndex_;
};

//-------------------------------------------------------------------------------
template <typename ScalarT>
class iNoiseOp: public astNode<ScalarT>
{
  public:
    iNoiseOp (): astNode<ScalarT>(), number_(0.0), derivIndex_(-1) {};
    virtual ScalarT val() {return number_;}
    virtual ScalarT dx(int i) { return (derivIndex_==i)?1.0:0.0; }
    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = number_;
      if ( !(derivs.empty() ) )
      {
        std::fill(derivs.begin(),derivs.end(),0.0);
        if (derivIndex_>-1) derivs[derivIndex_] = 1.0;
      }
    }

    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "inoise variable id = " << this->id_ << std::endl;
      os << std::setw(indent) << " " << "value = " << val() <<std::endl;
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "inoise variable id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os ) { os << "INOISE"; }

    virtual void setDerivIndex(int i) {derivIndex_=i;};
    virtual void unsetDerivIndex() {derivIndex_=-1;};
    ScalarT & getNoiseVar () { return number_; }
    void setNoiseVar (ScalarT n) { number_ = n; }
    virtual bool getIsTreeConstant() { return false; }
    virtual bool iNoiseType()  { return true; };

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<iNoiseOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<iNoiseOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
    } // 2nd dispatch

  private:
    ScalarT number_;
    int derivIndex_;
};

//-------------------------------------------------------------------------------
// This class represents the invocation of a .FUNC (function).  It represents a
// specific usage of a function that was declared elsewhere.
//
// This needs a variable number of arguments, and those args that are passed
// in are the unique ones for this usage.
//
// The functionNode_ pointer is the "top" node of an AST tree for the function
// that is being invoked by this class.  It has to be resolved externally, as
// it was defined somewhere else.  It is a separate expression.
//
// The "dummy" args are also defined externally.  They are the generic arguments
// that were part of the .FUNC declaration.
//
// So, for example:
//
// .func f(x,y) { x+2*y } <- function declaration, set up elsewhere.
// { 10*f(2,3) }          <- expression using the function
//
//external to this AST:
// {x+2*y} = pointed to by functionNode_
// x,y     = dummyFuncArgs_  - parameter nodes, used by functionNode_ AST. Their setNode function must be called to set them to specific values
//
//internal to this funcOp class:
// 2,3     = funcArgs_ - generic AST nodes.  Could be anything.  These are the nodes that are "set" onto the dummyArgs.
// f(2,3)  = function call handled by this class
//
//internal to this AST:
// 10      = handled elsewhere in this AST tree as a numval node
// *       = handled elsewhere in this AST tree as a binaryMultOp node
//
//-------------------------------------------------------------------------------
template <typename ScalarT>
class funcOp: public astNode<ScalarT>
{
  public:
    // functions:
    funcOp (const std::string & name, std::vector<Teuchos::RCP<astNode<ScalarT> > > & args):
      astNode<ScalarT>(),
      funcName_(name),
      funcArgs_(args),
      nodeResolved_(false),
      argsResolved_(false),
      sdtNodesResolved_(false),
      ddtNodesResolved_(false)
    {};

    //-------------------------------------------------------------------------------
    void setArgs()
    {
      for (int ii=0;ii<dummyFuncArgs_.size();++ii)
      {
        dummyFuncArgs_[ii]->setNode( funcArgs_[ii] );
      }

      if (!(sdtNodes_.empty()))
      {
        std::string stateKey;
        for (int ii=0;ii<dummyFuncArgs_.size();++ii)
        {
          stateKey += std::to_string(  funcArgs_[ii]->getNodeId() );
          if (ii<dummyFuncArgs_.size()-1) stateKey += std::string("_");
        }
        if( sdtStateVecMap_.find(stateKey) == sdtStateVecMap_.end() )
        {
          std::vector<sdtStateData<ScalarT> > sdtStateVec( sdtNodes_.size() );
          sdtStateVecMap_[stateKey] = sdtStateVec;
        }

        for (int ii=0;ii<sdtNodes_.size();ii++)
        {
          saved_sdtStateVec_[ii] = sdtNodes_[ii]->getSdtState();
          std::vector<sdtStateData<ScalarT> > & ssVec = sdtStateVecMap_[stateKey];
          sdtNodes_[ii]->setSdtState( ssVec[ii] );
        }
      }

      if (!(ddtNodes_.empty()))
      {
        std::string stateKey;
        for (int ii=0;ii<dummyFuncArgs_.size();++ii)
        {
          stateKey += std::to_string(  funcArgs_[ii]->getNodeId() );
          if (ii<dummyFuncArgs_.size()-1) stateKey += std::string("_");
        }
        if( ddtStateVecMap_.find(stateKey) == ddtStateVecMap_.end() )
        {
          std::vector<ddtStateData<ScalarT> > ddtStateVec( ddtNodes_.size() );
          ddtStateVecMap_[stateKey] = ddtStateVec;
        }

        for (int ii=0;ii<ddtNodes_.size();ii++)
        {
          saved_ddtStateVec_[ii] = ddtNodes_[ii]->getDdtState();
          std::vector<ddtStateData<ScalarT> > & ssVec = ddtStateVecMap_[stateKey];
          ddtNodes_[ii]->setDdtState( ssVec[ii] );
        }
      }
    }

    //-------------------------------------------------------------------------------
    void unsetArgs()
    {
      if (!(sdtNodes_.empty()))
      {
        std::string stateKey;
        for (int ii=0;ii<dummyFuncArgs_.size();++ii)
        {
          stateKey += std::to_string(  funcArgs_[ii]->getNodeId() );
          if (ii<dummyFuncArgs_.size()-1) stateKey += std::string("_");
        }
        for (int ii=0;ii<sdtNodes_.size();ii++)
        {
          std::vector<sdtStateData<ScalarT> > & ssVec = sdtStateVecMap_[stateKey];
          ssVec[ii] = sdtNodes_[ii]->getSdtState();
          sdtNodes_[ii]->setSdtState(saved_sdtStateVec_[ii]);
        }
      }

      if (!(ddtNodes_.empty()))
      {
        std::string stateKey;
        for (int ii=0;ii<dummyFuncArgs_.size();++ii)
        {
          stateKey += std::to_string(  funcArgs_[ii]->getNodeId() );
          if (ii<dummyFuncArgs_.size()-1) stateKey += std::string("_");
        }
        for (int ii=0;ii<ddtNodes_.size();ii++)
        {
          std::vector<ddtStateData<ScalarT> > & ssVec = ddtStateVecMap_[stateKey];
          ssVec[ii] = ddtNodes_[ii]->getDdtState();
          ddtNodes_[ii]->setDdtState(saved_ddtStateVec_[ii]);
        }
      }

      for (int ii=0;ii<dummyFuncArgs_.size();++ii)
      {
        dummyFuncArgs_[ii]->unsetNode();
      }
    }

    //-------------------------------------------------------------------------------
    virtual ScalarT val()
    {
      ScalarT val = 0.0;
      if (nodeResolved_ && argsResolved_)
      {
        if (funcArgs_.size() == dummyFuncArgs_.size())
        {
          setArgs();
          val = functionNode_->val();
          unsetArgs();
        }
      }
      else
      {
        std::string tmp = "Function " + funcName_ + " is not resolved";
        std::vector<std::string> errStr(1,tmp);
        yyerror(errStr);
      }
      return val;
    }

    //-------------------------------------------------------------------------------
    virtual ScalarT dx(int i)
    {
      ScalarT dfdx = 0.0;
      if (nodeResolved_ && argsResolved_)
      {
        if (funcArgs_.size() == dummyFuncArgs_.size())
        {
          // Two phases, do do a complete chain rule calculation:
          //
          // all the "d" symbols should be partials:
          // chain rule :  F′(x) = F'(x) + f′(g(x)) * g′(x)  ->  df/dx = df/dx + df/dp * dp/dx
          //
          // phase 1:  F'(x) = df/dx.
          //
          // For this phase, the funcArgs are in the "full" form -
          //   ie, if they represent an AST tree, we use the whole tree to evaluate.
          setArgs();
          dfdx = functionNode_->dx(i);
          unsetArgs();

          // phase 2:  f′(g(x)) * g′(x) = df/dp * dp/dx
          //
          // g(x) = funcArg->val().  This should be evaluated inside of dx call.
          // g'(x) = funcArg->dx(i)
          // f'(g(x)) = functionNode_->dx(ii);
          //
          // For this phase, the funcArgs are simple params.
          //   ie, they don't have an AST tree, just a number.
          for (int ii=0;ii<dummyFuncArgs_.size();++ii)
          {
            // the index is intentionally negative, starting at -1.  This is so it
            // doesn't conflict with derivative indices that were already set at
            // the top of the tree in the newExpression::evaluate function.
            int index=-ii-1;
            dummyFuncArgs_[ii]->setValue ( funcArgs_[ii]->val() );
            dummyFuncArgs_[ii]->setDerivIndex ( index );
          }

          // This can be a big slowdown, for nested function calls.  num1 (functionNode_->dx) is the problem.
          for (int ii=0;ii<dummyFuncArgs_.size();++ii)  // loop over args (p).  ii = p index, i = x index
          {
            int index=-ii-1;
            ScalarT delta = 0.0;
            ScalarT dpdx = funcArgs_[ii]->dx(i); // usually zero ...
            if (dpdx != 0.0) { delta = dpdx * functionNode_->dx(index); } // slow. do not evaluate if not needed.
            dfdx += delta;
          }

          for (int ii=0;ii<dummyFuncArgs_.size();++ii)
          {
            dummyFuncArgs_[ii]->unsetValue ();
            dummyFuncArgs_[ii]->unsetDerivIndex ();
          } // restore
        }
      }
      return dfdx;
    }

    //-------------------------------------------------------------------------------
    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      std::vector<ScalarT> & dfdx = derivs;
      int numDerivs = derivs.size();

      if (nodeResolved_ && argsResolved_)
      {
        if (funcArgs_.size() == dummyFuncArgs_.size())
        {
          // Two phases, do do a complete chain rule calculation:
          //
          // all the "d" symbols should be partials:
          // chain rule :  F′(x) = F'(x) + f′(g(x)) * g′(x)  ->  df/dx = df/dx + df/dp * dp/dx
          //
          // phase 1:  F'(x) = df/dx.
          //
          // For this phase, the funcArgs are in the "full" form -
          //   ie, if they represent an AST tree, we use the whole tree to evaluate.
          setArgs();
          functionNode_->dx2(result,dfdx);
          unsetArgs();

          // phase 2:  f′(g(x)) * g′(x) = df/dp * dp/dx
          //
          // g(x) = funcArg->val().  This should be evaluated inside of dx call.
          // g'(x) = funcArg->dx(i)
          // f'(g(x)) = functionNode_->dx(ii);
          //
          // For this phase, the funcArgs are simple params.
          //   ie, they don't have an AST tree, just a number.
          for (int ii=0;ii<dummyFuncArgs_.size();++ii)
          {
            // the index is intentionally negative, starting at -1.  This is so it
            // doesn't conflict with derivative indices that were already set at
            // the top of the tree in the newExpression::evaluate function.
            int index=-ii-1;
            dummyFuncArgs_[ii]->setValue ( funcArgs_[ii]->val() );
            dummyFuncArgs_[ii]->setDerivIndex ( index );
          }

          // This can be a big slowdown, for nested function calls.  num1 (functionNode_->dx) is the problem.
          for (int ii=0;ii<dummyFuncArgs_.size();++ii)  // loop over args (p).  ii = p index, i = x index
          {
            int index=-ii-1;

            ScalarT pval;
            std::vector<ScalarT> dpdx_;
            dpdx_.resize(numDerivs,0.0);
            funcArgs_[ii]->dx2(pval,dpdx_); // usually zero ...

            for (int jj=0;jj<numDerivs;jj++)
            {
              ScalarT delta = 0.0;
              if (dpdx_[jj] != 0.0) { delta = dpdx_[jj] * functionNode_->dx(index); } // slow. do not evaluate if not needed.
              dfdx[jj] += delta;
            }
          }

          for (int ii=0;ii<dummyFuncArgs_.size();++ii)
          {
            dummyFuncArgs_[ii]->unsetValue ();
            dummyFuncArgs_[ii]->unsetDerivIndex ();
          } // restore
        }
      }
    }

    //-------------------------------------------------------------------------------
    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " " << "function: " << funcName_ << ": id = " << this->id_ << std::endl;
      os << std::setw(indent) << " " << "function args: " << std::endl; indent++;
      for (int ii=0;ii<funcArgs_.size();++ii) { funcArgs_[ii]->output(os,indent+1); }

      if( !(Teuchos::is_null(functionNode_)) )
      {
        os << std::setw(indent) << " " << "functionNode_ ("<<funcName_<<") details: " << std::endl;

        if (funcArgs_.size() == dummyFuncArgs_.size())
        {
          for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->setNode( funcArgs_[ii] ); }
          functionNode_->output(os,indent+2);
          for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->unsetNode(); } // restore
        }

        os << std::setw(indent) << " " << "val = " << val() << std::endl;
      }
      else
      {
        os << std::setw(indent) << " " << "functionNode_ is not resolved " << std::endl;
      }
    }

    //-------------------------------------------------------------------------------
    virtual void compactOutput(std::ostream & os)
    {
      os << " " << "function: " << funcName_ << ": id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << funcName_;
      os << "(";
      int size = funcArgs_.size();
      for (int ii=0;ii<size;ii++)
      {
        funcArgs_[ii]->codeGen(os);
        if (ii < size-1) {os << ",";}
      }
      os << ")";
    }

    virtual void setNode(Teuchos::RCP<astNode<ScalarT> > & tmpNode) { functionNode_ = tmpNode; nodeResolved_ = true; };
    virtual void unsetNode()
    {
      nodeResolved_ = false;
    };

    virtual void setFuncArgs(const std::vector< Teuchos::RCP<paramOp<ScalarT> > > & tmpParamVec )
    {
      dummyFuncArgs_.clear(); dummyFuncArgs_.resize(tmpParamVec.size());
      for (int ii=0;ii<tmpParamVec.size();++ii)
      {
        dummyFuncArgs_[ii] = tmpParamVec[ii];
      }
      argsResolved_ = true;
    };

    virtual void setFuncArgs(const std::vector< Teuchos::RCP<astNode<ScalarT> > > & tmpArgVec )
    {
      dummyFuncArgs_.clear(); dummyFuncArgs_.resize(tmpArgVec.size());
      for (int ii=0;ii<tmpArgVec.size();++ii)
      {
        dummyFuncArgs_[ii] = tmpArgVec[ii];
      }
      argsResolved_ = true;
    };

    std::vector< Teuchos::RCP<astNode<ScalarT> > > & getFuncArgs()
    {
      return funcArgs_;
    };

    virtual void setSdtArgs(const std::vector< Teuchos::RCP<astNode<ScalarT> > > & tmpSdtVec )
    {
      sdtNodes_.clear(); sdtNodes_.resize(tmpSdtVec.size());
      sdtArgNodes_.clear(); sdtArgNodes_.resize(tmpSdtVec.size());
      saved_sdtStateVec_.clear(); saved_sdtStateVec_.resize(tmpSdtVec.size());
      for (int ii=0;ii<tmpSdtVec.size();++ii)
      {
        sdtNodes_[ii] = tmpSdtVec[ii];
        Teuchos::RCP<sdtOp<ScalarT> > castedSdtPtr = Teuchos::rcp_dynamic_cast<sdtOp<ScalarT> > (tmpSdtVec[ii]);
        sdtArgNodes_[ii] = castedSdtPtr->getArg();
      }
      sdtNodesResolved_ = true;
    };

    virtual void setDdtArgs(const std::vector< Teuchos::RCP<astNode<ScalarT> > > & tmpDdtVec )
    {
      ddtNodes_.clear();    ddtNodes_.resize(tmpDdtVec.size());
      ddtArgNodes_.clear(); ddtArgNodes_.resize(tmpDdtVec.size());
      saved_ddtStateVec_.clear(); saved_ddtStateVec_.resize(tmpDdtVec.size());
      for (int ii=0;ii<tmpDdtVec.size();++ii)
      {
        ddtNodes_[ii] = tmpDdtVec[ii];
        Teuchos::RCP<ddtOp<ScalarT> > castedDdtPtr = Teuchos::rcp_dynamic_cast<ddtOp<ScalarT> > (tmpDdtVec[ii]);
        ddtArgNodes_[ii] = castedDdtPtr->getArg();
      }
      ddtNodesResolved_ = true;
    };

    virtual bool funcType()    { return true; };

    virtual std::string getName() { return funcName_; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<funcOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<funcOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch

      if( !(Teuchos::is_null( functionNode_)))
      {
        if(dummyFuncArgs_.size() == funcArgs_.size())
          for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->setNode( funcArgs_[ii] ); }
        functionNode_->accept(visitor, functionNode_);
        if(dummyFuncArgs_.size() == funcArgs_.size())
          for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->unsetNode(); } // restore
      }
    } // 2nd dispatch

    bool getNodeResolved() { return nodeResolved_; }
    bool getArgsResolved() { return argsResolved_; }

    virtual void processSuccessfulTimeStep ()
    {
      functionNode_->processSuccessfulTimeStep ();
    };

    ddtStateData<ScalarT> & getDdtState() { return functionNode_->getDdtState(); }
    sdtStateData<ScalarT> & getSdtState() { return functionNode_->getSdtState(); }

    virtual unsigned long int getNodeId () { return functionNode_->getNodeId(); }

    virtual bool getIsComplex ()
    {
      bool isComplex = (typeid(ScalarT) == typeid(std::complex<double>));

      if( !(Teuchos::is_null( functionNode_)))
      {
        if(dummyFuncArgs_.size() == funcArgs_.size())
          for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->setNode( funcArgs_[ii] ); }

        isComplex = functionNode_->getIsComplex();

        if(dummyFuncArgs_.size() == funcArgs_.size())
          for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->unsetNode(); } // restore
      }

      return isComplex ;
    }

    virtual bool getIsTreeConstant()
    {
      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->setNode( funcArgs_[ii] ); }

      bool isConstant = functionNode_->getIsTreeConstant();

      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->unsetNode(); } // restore

      return isConstant;
    }

  private:
// data:
    std::string funcName_;
    std::vector<Teuchos::RCP<astNode<ScalarT> > > funcArgs_;  // the unique args that are passed in to this instance of a function
    std::vector<Teuchos::RCP<astNode<ScalarT> > > dummyFuncArgs_;  // generic args that the functionNode_ owns; ie, that are used to evaluate it.  They have to be temporarily replaced whenever the function is called.

    std::vector<Teuchos::RCP<astNode<ScalarT> > > sdtNodes_;
    std::vector<Teuchos::RCP<astNode<ScalarT> > > ddtNodes_;

    std::vector<Teuchos::RCP<astNode<ScalarT> > > sdtArgNodes_;
    std::vector<Teuchos::RCP<astNode<ScalarT> > > ddtArgNodes_;

    std::unordered_map<std::string, std::vector<sdtStateData<ScalarT> > > sdtStateVecMap_;
    std::unordered_map<std::string, std::vector<ddtStateData<ScalarT> > > ddtStateVecMap_;

    std::vector<sdtStateData<ScalarT> > saved_sdtStateVec_;
    std::vector<ddtStateData<ScalarT> > saved_ddtStateVec_;

    Teuchos::RCP<astNode<ScalarT> > functionNode_;
    bool nodeResolved_;
    bool argsResolved_;
    bool sdtNodesResolved_;
    bool ddtNodesResolved_;
};

//-------------------------------------------------------------------------------
// sign corrected x raised to y power
//
//  pwrs(x,y) = x^y     if x > 0
//  pwrs(x,y) = 0       if x = 0
//  pwrs(x,y) = -(-x)^y if x < 0
//
template <typename ScalarT>
class pwrsOp : public astNode<ScalarT>
{
  public:
    pwrsOp (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):
      astNode<ScalarT>(left,right), rightConst_(true),leftConst_(false)
    {
      rightConst_ = this->rightAst_->numvalType();
      leftConst_ = this->leftAst_->numvalType();
    }

    virtual ScalarT val()
    {
      ScalarT ret=0.0;
      ScalarT leftVal=this->leftAst_->val();
      if (std::real(leftVal) >= 0)
      {
        ret = std::pow(leftVal, this->rightAst_->val());
      }
      else if (std::real(leftVal) < 0)
      {
        ret = -std::pow(-(leftVal), this->rightAst_->val());
      }
      return ret;
    }

    virtual ScalarT dx (int i)
    {
      Teuchos::RCP<astNode<ScalarT> > & lef = this->leftAst_;
      Teuchos::RCP<astNode<ScalarT> > & rig = this->rightAst_;
      ScalarT retVal = 0.0;

      ScalarT leftVal=lef->val();
      ScalarT righVal=rig->val();

      if (leftVal != 0.0)
      {
        if (rightConst_ && !leftConst_)
        {
          if (std::real(leftVal) >= 0)
          {
            retVal = righVal*lef->dx(i)/leftVal*std::pow(leftVal,righVal) ;
          }
          else if (std::real(leftVal) < 0)
          {
            retVal = righVal*(-lef->dx(i))/(-leftVal)*std::pow((-leftVal),righVal) ;
          }
        }
        else if (!rightConst_ && leftConst_)
        {
          if (std::real(leftVal) >= 0)
          {
            retVal = std::log(leftVal)*std::pow(leftVal,righVal)*rig->dx(i);
          }
          else if (std::real(leftVal) < 0)
          {
            retVal = -std::log(-leftVal)*std::pow(-leftVal,righVal)*rig->dx(i);
          }
        }
        else
        {
          if (std::real(leftVal) >= 0)
          {
            retVal = (rig->dx(i)*std::log(leftVal)+righVal*lef->dx(i)/leftVal)*std::pow(leftVal,righVal);
          }
          else if (std::real(leftVal) < 0)
          {
            retVal = (-rig->dx(i)*std::log(-leftVal)+righVal*(-lef->dx(i))/(-leftVal))*std::pow((-leftVal),righVal);
          }
        }
      }
      return  retVal;
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      Teuchos::RCP<astNode<ScalarT> > & lef = this->leftAst_;
      Teuchos::RCP<astNode<ScalarT> > & rig = this->rightAst_;

      ScalarT leftVal;
      ScalarT righVal;
      std::vector<ScalarT> lefDerivs_;
      std::vector<ScalarT> rigDerivs_;

      int numDerivs=derivs.size();
      lefDerivs_.resize(numDerivs,0.0);
      rigDerivs_.resize(numDerivs,0.0);
      lef->dx2(leftVal,lefDerivs_);
      rig->dx2(righVal,rigDerivs_);

      if (std::real(leftVal) >= 0)     { result = std::pow(leftVal, righVal); }
      else if (std::real(leftVal) < 0) { result = -std::pow(-(leftVal), righVal); }

      if (leftVal != 0.0)
      {
        if (rightConst_ && !leftConst_)
        {
          lef->dx2(leftVal,lefDerivs_);
          righVal = rig->val();
          if (std::real(leftVal) >= 0)
          {
            for (int ii=0;ii<numDerivs;ii++)
            {
              derivs[ii] = righVal*lefDerivs_[ii]/leftVal*std::pow(leftVal,righVal) ;
            }
          }
          else if (std::real(leftVal) < 0)
          {
            for (int ii=0;ii<numDerivs;ii++)
            {
              derivs[ii] = righVal*(-lefDerivs_[ii])/(-leftVal)*std::pow((-leftVal),righVal) ;
            }
          }
        }
        else if (!rightConst_ && leftConst_)
        {
          if (std::real(leftVal) >= 0)
          {
            for (int ii=0;ii<numDerivs;ii++)
            {
              derivs[ii] = std::log(leftVal)*std::pow(leftVal,righVal)*rigDerivs_[ii];
            }
          }
          else if (std::real(leftVal) < 0)
          {
            for (int ii=0;ii<numDerivs;ii++)
            {
              derivs[ii] = -std::log(-leftVal)*std::pow(-leftVal,righVal)*rigDerivs_[ii];
            }
          }
        }
        else
        {

          if (std::real(leftVal) >= 0)
          {
            for (int ii=0;ii<numDerivs;ii++)
            {
              derivs[ii] = (rigDerivs_[ii]*std::log(leftVal)+righVal*lefDerivs_[ii]/leftVal)*std::pow(leftVal,righVal);
            }
          }
          else if (std::real(leftVal) < 0)
          {
            for (int ii=0;ii<numDerivs;ii++)
            {
              derivs[ii] = (-rigDerivs_[ii]*std::log(-leftVal)+righVal*(-lefDerivs_[ii])/(-leftVal))*std::pow((-leftVal),righVal);
            }
          }
        }
      }
    }

    virtual bool getIsComplex ()
    {
      bool isComplex = (this->rightAst_->getIsComplex() || this->leftAst_->getIsComplex());

      if (!isComplex)
      {
        if (  std::real(this->leftAst_->val()) < 0.0  && std::abs(this->rightAst_->val()) < 1.0 )
        {
          isComplex = true;
        }
      }
      return isComplex;
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "pwrs operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "pwrs operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      if (std::real(this->leftAst_->val()) < 0) { os << "-"; }
      os << "std::pow(";
      if (std::real(this->leftAst_->val()) < 0) { os << "-"; }
      this->leftAst_->codeGen(os);
      os << ",";
      this->rightAst_->codeGen(os);
      os << ")";
    }

    virtual bool getIsTreeConstant()
    { return (this->leftAst_->getIsTreeConstant() && this->rightAst_->getIsTreeConstant()); }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<pwrsOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<pwrsOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
      this->leftAst_->accept(visitor, this->leftAst_); 
      this->rightAst_->accept(visitor, this->rightAst_); 
    }

  private:
    bool rightConst_;
    bool leftConst_;
};

//-------------------------------------------------------------------------------
//
// +1 if x > 0
//  0 if x = 0
// -1 if x < 0
//
template <typename ScalarT>
class sgnOp : public astNode<ScalarT>
{
  public:
    sgnOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val()
    {
      ScalarT ret = 0.0;
      double leftVal = std::real(this->leftAst_->val());
      ret = ( (leftVal>0)?+1:ret );
      ret = ( (leftVal<0)?-1:ret );
      return ret;
    }

    virtual ScalarT dx(int i)
    {
      return ScalarT(0.0);
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0); }
    }

    virtual bool getIsComplex () { return false; } // val() only considers real part

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "sgn operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "sgn operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "signbit(";
      this->leftAst_->codeGen(os);
      os << ")";
    }

    virtual bool getIsTreeConstant() { return this->leftAst_->getIsTreeConstant(); }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<sgnOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<sgnOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
      this->leftAst_->accept(visitor, this->leftAst_); 
    }
};

//-------------------------------------------------------------------------------
//  sign(x,y) = sgn(y)|x|  sign of y times absolute value of x
template <typename ScalarT>
class signOp : public astNode<ScalarT>
{
  public:
    signOp (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right): astNode<ScalarT>(left,right) {};

    virtual ScalarT val()
    {
      ScalarT y = 0.0;
      double realRightVal = std::real(this->rightAst_->val());
      y = ( (realRightVal>0)?+1:y );
      y = ( (realRightVal<0)?-1:y );

      ScalarT x =  (std::abs(this->leftAst_->val()));
      return (y*x);
    }

    virtual ScalarT dx (int i)
    {
      ScalarT y = 0.0;
      double realRightVal = std::real(this->rightAst_->val());
      y = ( (realRightVal>0)?+1:y );
      y = ( (realRightVal<0)?-1:y );
      ScalarT dx = (std::real(this->leftAst_->val()) >= 0 ? this->leftAst_->dx(i) : ScalarT(-this->leftAst_->dx(i)));
      return (y*dx);
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      ScalarT y = 0.0;

      Teuchos::RCP<astNode<ScalarT> > & lef = this->leftAst_;
      Teuchos::RCP<astNode<ScalarT> > & rig = this->rightAst_;

      ScalarT righVal=rig->val();
      std::vector<ScalarT> lefDerivs_;
      int numDerivs = derivs.size();
      lefDerivs_.resize(numDerivs,0.0);
      ScalarT leftVal;
      lef->dx2(leftVal,lefDerivs_);

      double realRightVal = std::real(righVal);
      y = ( (realRightVal>0)?+1:y );
      y = ( (realRightVal<0)?-1:y );

      ScalarT x =  (std::abs(leftVal));
      result = (y*x);

      for (int i=0;i<numDerivs;i++)
      {
        ScalarT dx = (std::real(leftVal)) >= 0 ?  lefDerivs_[i] : ScalarT(-lefDerivs_[i]);
        derivs[i] = (y*dx);
      }

    }

    virtual bool getIsComplex () { return false; } // val() only considers real part of right, and mag(left)

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "sign(x,y) = (sgn(y)|x|) op id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "sign(x,y) = (sgn(y)|x|) op id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "signbit(";
      this->rightAst_->codeGen(os);
      os << ")*std::abs(";
      this->leftAst_->codeGen(os);
      os << ")";
    }

    virtual bool getIsTreeConstant()
    { return (this->leftAst_->getIsTreeConstant() && this->rightAst_->getIsTreeConstant()); }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<signOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<signOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
      this->leftAst_->accept(visitor, this->leftAst_); 
      this->rightAst_->accept(visitor, this->rightAst_);
    }
};

//-------------------------------------------------------------------------------
// fmod operator. returns the remainder of the division as a real number
//
template <typename ScalarT>
class fmodOp : public astNode<ScalarT>
{
  public:
    fmodOp (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):
      astNode<ScalarT>(left,right), rightConst_(true),leftConst_(false), bpTol_(0.0)
    {
      rightConst_ = this->rightAst_->numvalType();
      leftConst_ = this->leftAst_->numvalType();
    }

    virtual ScalarT val()
    {
#if 0
      // copied from the stpOp class:  rewrite this for fmod.
      // stpOp returns a 1 or a 0.
      Teuchos::RCP<astNode<ScalarT> > zeroAst_ = Teuchos::rcp(new numval<ScalarT>(0.0));
      bpTimes_.clear();
      computeBreakPoint ( this->leftAst_, zeroAst_, timeOpVec_, bpTol_, bpTimes_);
#endif

      return std::fmod ( std::real(this->leftAst_->val()) , std::real(this->rightAst_->val()));
    }

    virtual ScalarT dx(int i)
    {
      ScalarT leftVal=this->leftAst_->val();
      ScalarT rightVal=this->rightAst_->val();
      ScalarT leftDx = 0.0; ScalarT rightDx = 0.0;

      double res = fabs((std::real(leftVal))/(std::real(rightVal)));
      Xyce::Util::fixNan(res);
      double floorRes = ((std::real(leftVal))>0)?(std::floor(res)):(-std::floor(res));

      // 𝑎′(𝑥)−𝑏′(𝑥)[𝑎(𝑥)/𝑏(𝑥)]
      if (!leftConst_)
      {
        leftDx = this->leftAst_->dx(i);
      }
      if (!rightConst_)
      {
        rightDx = this->rightAst_->dx(i);
      }
      return leftDx-rightDx*floorRes;
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();

      int numDerivs = derivs.size();
      ScalarT leftVal, rightVal, leftDx=0.0, rightDx=0.0;
      std::vector<ScalarT> lefDerivs_;
      std::vector<ScalarT> rigDerivs_;
      if (leftConst_)
      {
        leftVal = this->leftAst_->val();
      }
      else
      {
        lefDerivs_.resize(numDerivs,0.0);
        this->leftAst_->dx2(leftVal,lefDerivs_);
      }
      if (rightConst_)
      {
        rightVal = this->rightAst_->val();
      }
      else
      {
        rigDerivs_.resize(numDerivs,0.0);
        this->rightAst_->dx2(rightVal,rigDerivs_);
      }

      double res = fabs((std::real(leftVal))/(std::real(rightVal)));
      Xyce::Util::fixNan(res);
      double floorRes = ((std::real(leftVal))>0)?(std::floor(res)):(-std::floor(res));

      // 𝑎′(𝑥)−𝑏′(𝑥)[𝑎(𝑥)/𝑏(𝑥)]
      for (int i=0;i<numDerivs;i++)
      {
        if (!leftConst_)
        {
          leftDx = lefDerivs_[i];
        }
        if (!rightConst_)
        {
          rightDx = rigDerivs_[i];
        }
        derivs[i] = leftDx-rightDx*floorRes;
      }
    }

    virtual bool getBreakPoints(std::vector<Xyce::Util::BreakPoint> & breakPointTimes)
    {
      if(!(bpTimes_.empty()))
      {
        for (int ii=0;ii<bpTimes_.size();ii++)
        {
          breakPointTimes.push_back(bpTimes_[ii]);
        }
      }

      return true;
    }

    virtual void setBreakPointTol(double tol) { bpTol_ = tol; }

    virtual bool getIsComplex () { return false; } // val() only considers real parts

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "fmod operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "fmod operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      if (std::real(this->leftAst_->val()) < 0) { os << "-"; }
      os << "std::fmod(";
      if (std::real(this->leftAst_->val()) < 0) { os << "-"; }
      this->leftAst_->codeGen(os);
      os << ",";
      this->rightAst_->codeGen(os);
      os << ")";
    }

    virtual bool getIsTreeConstant()
    { return (this->leftAst_->getIsTreeConstant() && this->rightAst_->getIsTreeConstant()); }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<fmodOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<fmodOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
      this->leftAst_->accept(visitor, this->leftAst_); 
      this->rightAst_->accept(visitor, this->rightAst_);
    }

  private:
    bool rightConst_;
    bool leftConst_;
    std::vector<Teuchos::RCP<astNode<ScalarT> > > timeOpVec_;
    double bpTol_;
    std::vector<Xyce::Util::BreakPoint> bpTimes_;
};

//-------------------------------------------------------------------------------
// nint(x) Rounds x up or down, to the nearest integer
template <typename ScalarT>
class roundOp : public astNode<ScalarT>
{
  public:
    roundOp(Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val() { return std::round(std::real(this->leftAst_->val())); }

    // derivative is undefined at integers and 0.0 elsewhere
    virtual ScalarT dx(int i) { return  0.0; }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0); }
    }

    virtual bool getIsComplex () { return false; } // val() only considers real part

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "round operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "round operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::round(";
      this->leftAst_->codeGen(os);
      os << ")";
    }

    virtual bool getIsTreeConstant()
    { return this->leftAst_->getIsTreeConstant(); }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<roundOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<roundOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
      this->leftAst_->accept(visitor, this->leftAst_); 
    }
};

//-------------------------------------------------------------------------------
// least integer greater or equal to variable x
template <typename ScalarT>
class ceilOp : public astNode<ScalarT>
{
  public:
    ceilOp(Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val() { return std::ceil(std::real(this->leftAst_->val())); }

    // derivative is undefined at integers and 0.0 elsewhere
    virtual ScalarT dx(int i) { return  0.0; }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0); }
    }

    virtual bool getIsComplex () { return false; } // val() only considers real part

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "ceil operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "ceil operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::ceil(";
      this->leftAst_->codeGen(os);
      os << ")";
    }
    virtual bool getIsTreeConstant()
    { return this->leftAst_->getIsTreeConstant(); }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<ceilOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<ceilOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
      this->leftAst_->accept(visitor, this->leftAst_); 
    }
};

//-------------------------------------------------------------------------------
// greatest integer less than or equal to variable x
template <typename ScalarT>
class floorOp : public astNode<ScalarT>
{
  public:
    floorOp(Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val() { return std::floor(std::real(this->leftAst_->val())); }

    // derivative is undefined at integers and 0.0 elsewhere
    virtual ScalarT dx(int i) { return  0.0; }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0); }
    }

    virtual bool getIsComplex () { return false; } // val() only considers real part

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "floor operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "floor operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::floor(";
      this->leftAst_->codeGen(os);
      os << ")";
    }
    virtual bool getIsTreeConstant()
    { return this->leftAst_->getIsTreeConstant(); }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<floorOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<floorOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
      this->leftAst_->accept(visitor, this->leftAst_); 
    }
};

//-------------------------------------------------------------------------------
// integer part of the real variable x
template <typename ScalarT>
class intOp : public astNode<ScalarT>
{
  public:
    intOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val()
    {
      int tmp = std::real(this->leftAst_->val()) ;
      ScalarT ret = tmp;
      return ret;
    }

    virtual ScalarT dx(int i) { return  0.0; }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0); }
    }

    virtual bool getIsComplex () { return false; } // val() only considers real part

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "int operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "int operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "static_cast<int>( std::real(";
      this->leftAst_->codeGen(os);
      os << "))";
    }

    virtual bool getIsTreeConstant()
    { return this->leftAst_->getIsTreeConstant(); }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<intOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<intOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
      this->leftAst_->accept(visitor, this->leftAst_); 
    }
};

//-------------------------------------------------------------------------------
// if statement op (or ternary op)
template <typename ScalarT>
class ifStatementOp : public astNode<ScalarT>
{
  public:
    ifStatementOp (Teuchos::RCP<astNode<ScalarT> > &xAst, Teuchos::RCP<astNode<ScalarT> > &yAst, Teuchos::RCP<astNode<ScalarT> > &zAst):
      astNode<ScalarT>(xAst,yAst),
      zAst_(zAst) {};

    virtual ScalarT val()
    {
      Teuchos::RCP<astNode<ScalarT> > & x = (this->leftAst_);
      Teuchos::RCP<astNode<ScalarT> > & y = (this->rightAst_);
      Teuchos::RCP<astNode<ScalarT> > & z = (zAst_);

      // not "fixing" x->val() b/c it is the result of a conditional, which is 1 or 0.
      // The correct place to fix this is in the comparison operators.  Fix later.
      ScalarT yFixed = y->val();  Xyce::Util::fixNan(yFixed);  Xyce::Util::fixInf(yFixed);
      ScalarT zFixed = z->val();  Xyce::Util::fixNan(zFixed);  Xyce::Util::fixInf(zFixed);
      return ((std::real(x->val()))?(yFixed):(zFixed));
    };

    virtual ScalarT dx (int i)
    {
      Teuchos::RCP<astNode<ScalarT> > & x = (this->leftAst_);
      Teuchos::RCP<astNode<ScalarT> > & y = (this->rightAst_);
      Teuchos::RCP<astNode<ScalarT> > & z = (zAst_);

      // not "fixing" x->val() b/c it is the result of a conditional, which is 1 or 0.
      // The correct place to fix this is in the comparison operators.  Fix later.
      ScalarT dyFixed = y->dx(i);  Xyce::Util::fixNan(dyFixed);  Xyce::Util::fixInf(dyFixed);
      ScalarT dzFixed = z->dx(i);  Xyce::Util::fixNan(dzFixed);  Xyce::Util::fixInf(dzFixed);
      return ((std::real(x->val()))?(dyFixed):(dzFixed));
    };

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      int numDerivs=derivs.size();
      std::vector<ScalarT> yDerivs_;
      std::vector<ScalarT> zDerivs_;
      yDerivs_.resize(numDerivs,0.0);
      zDerivs_.resize(numDerivs,0.0);

      ScalarT xVal = this->leftAst_->val();
      ScalarT yFixed, zFixed;
      this->rightAst_->dx2(yFixed,yDerivs_);
      zAst_->dx2(zFixed,zDerivs_);

      Xyce::Util::fixNan(yFixed);  Xyce::Util::fixInf(yFixed);
      Xyce::Util::fixNan(zFixed);  Xyce::Util::fixInf(zFixed);

      result = ((std::real(xVal))?(yFixed):(zFixed));

      // not "fixing" x->val() b/c it is the result of a conditional, which is 1 or 0.
      // The correct place to fix this is in the comparison operators.  Fix later.
      for (int ii=0;ii<numDerivs;ii++)
      {
        ScalarT dyFixed = yDerivs_[ii];  Xyce::Util::fixNan(dyFixed);  Xyce::Util::fixInf(dyFixed);
        ScalarT dzFixed = zDerivs_[ii];  Xyce::Util::fixNan(dzFixed);  Xyce::Util::fixInf(dzFixed);
        derivs[ii] = ((std::real(xVal))?(dyFixed):(dzFixed));
      }
    };

    // For the x-part (the conditional), only real part is used.  But y and z can be complex.
    virtual bool getIsComplex () { return ((this->rightAst_)->getIsComplex() || (zAst_)->getIsComplex()); }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "if statement operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
      zAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "if statement operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "((";
      this->leftAst_->codeGen(os);
      os << ")?(";
      this->rightAst_->codeGen(os);
      os << "):(";
      zAst_->codeGen(os);
      os << "))";
    }

    virtual bool getIsTreeConstant()
    { return
      (this->leftAst_->getIsTreeConstant() &&
       this->rightAst_->getIsTreeConstant() &&
       this->zAst_->getIsTreeConstant() ) ;
    }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<ifStatementOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<ifStatementOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
      this->leftAst_->accept(visitor, this->leftAst_); 
      this->rightAst_->accept(visitor, this->rightAst_); 
      zAst_->accept(visitor,zAst_);
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > zAst_;
};

//-------------------------------------------------------------------------------
// x limited to range y to z
// y if x < y
// x if y < x < z
// z if x > z
template <typename ScalarT>
class limitOp : public astNode<ScalarT>
{
  public:
    limitOp (Teuchos::RCP<astNode<ScalarT> > &xAst, Teuchos::RCP<astNode<ScalarT> > &yAst, Teuchos::RCP<astNode<ScalarT> > &zAst):
      astNode<ScalarT>(xAst,yAst),
      zAst_(zAst), bpTol_(0.0) {};

    virtual ScalarT val()
    {
      Teuchos::RCP<astNode<ScalarT> > & x = (this->leftAst_);
      Teuchos::RCP<astNode<ScalarT> > & y = (this->rightAst_);
      Teuchos::RCP<astNode<ScalarT> > & z = (zAst_);

      ScalarT xFixed = x->val();  Xyce::Util::fixNan(xFixed);  Xyce::Util::fixInf(xFixed);
      ScalarT yFixed = y->val();  Xyce::Util::fixNan(yFixed);  Xyce::Util::fixInf(yFixed);
      ScalarT zFixed = z->val();  Xyce::Util::fixNan(zFixed);  Xyce::Util::fixInf(zFixed);

      bpTimes_.clear();
      computeBreakPoint ( x, y, timeOpVec_, bpTol_, bpTimes_);
      computeBreakPoint ( x, z, timeOpVec_, bpTol_, bpTimes_);

      return ((std::real(xFixed)<std::real(yFixed))?
          std::real(yFixed):
          ((std::real(xFixed)>std::real(zFixed))?std::real(zFixed):std::real(xFixed)));
    };

    virtual ScalarT dx (int i)
    {
      Teuchos::RCP<astNode<ScalarT> > & x = (this->leftAst_);
      Teuchos::RCP<astNode<ScalarT> > & y = (this->rightAst_);
      Teuchos::RCP<astNode<ScalarT> > & z = (zAst_);

      ScalarT xFixed = x->val();  Xyce::Util::fixNan(xFixed);  Xyce::Util::fixInf(xFixed);
      ScalarT dxFixed = x->dx(i); Xyce::Util::fixNan(dxFixed); Xyce::Util::fixInf(dxFixed);
      ScalarT yFixed = y->val();  Xyce::Util::fixNan(yFixed);  Xyce::Util::fixInf(yFixed);
      ScalarT zFixed = z->val();  Xyce::Util::fixNan(zFixed);  Xyce::Util::fixInf(zFixed);

      return ((std::real(xFixed)<std::real(yFixed))?0.0:((std::real(xFixed)>std::real(zFixed))?0.0:(dxFixed)));
    };

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      Teuchos::RCP<astNode<ScalarT> > & x = (this->leftAst_);
      Teuchos::RCP<astNode<ScalarT> > & y = (this->rightAst_);
      Teuchos::RCP<astNode<ScalarT> > & z = (zAst_);

      int numDerivs = derivs.size();
      std::vector<ScalarT> xDerivs_;
      xDerivs_.resize(numDerivs,0.0);

      ScalarT xFixed;
      x->dx2(xFixed, xDerivs_);
      Xyce::Util::fixNan(xFixed);  Xyce::Util::fixInf(xFixed);

      ScalarT yFixed = y->val();  Xyce::Util::fixNan(yFixed);  Xyce::Util::fixInf(yFixed);
      ScalarT zFixed = z->val();  Xyce::Util::fixNan(zFixed);  Xyce::Util::fixInf(zFixed);

      bpTimes_.clear();
      computeBreakPoint ( x, y, timeOpVec_, bpTol_, bpTimes_);
      computeBreakPoint ( x, z, timeOpVec_, bpTol_, bpTimes_);

      result = ((std::real(xFixed)<std::real(yFixed))?
          std::real(yFixed):
          ((std::real(xFixed)>std::real(zFixed))?std::real(zFixed):std::real(xFixed)));

      for (int ii=0;ii<numDerivs;ii++)
      {
        ScalarT dxFixed = xDerivs_[ii]; Xyce::Util::fixNan(dxFixed); Xyce::Util::fixInf(dxFixed);
        derivs[ii] = ((std::real(xFixed)<std::real(yFixed))?0.0:((std::real(xFixed)>std::real(zFixed))?0.0:(dxFixed)));
      }
    };

    virtual bool getBreakPoints(std::vector<Xyce::Util::BreakPoint> & breakPointTimes)
    {
      if(!(bpTimes_.empty()))
      {
        for (int ii=0;ii<bpTimes_.size();ii++)
        {
          breakPointTimes.push_back(bpTimes_[ii]);
        }
      }

      return true;
    }

    virtual void setBreakPointTol(double tol) { bpTol_ = tol; }

    virtual bool getIsComplex () { return false; } // val() only uses real parts

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "limit operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
      zAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "limit operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      /// fix this
      os << "LIMIT";
    }

    virtual bool getIsTreeConstant()
    { return
      (this->leftAst_->getIsTreeConstant() &&
       this->rightAst_->getIsTreeConstant() &&
       this->zAst_->getIsTreeConstant() ) ;
    }

    virtual bool limitType() { return true; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<limitOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<limitOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
      this->leftAst_->accept(visitor, this->leftAst_); 
      this->rightAst_->accept(visitor, this->rightAst_); 
      zAst_->accept(visitor,zAst_);
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > zAst_;
    std::vector<Teuchos::RCP<astNode<ScalarT> > > timeOpVec_;
    double bpTol_;
    std::vector<Xyce::Util::BreakPoint> bpTimes_;
};

//-------------------------------------------------------------------------------
// step function :  1 if x > 0
template <typename ScalarT>
class stpOp : public astNode<ScalarT>
{
  public:
    stpOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left), bpTol_(0.0) { };

    virtual ScalarT val()
    {
      // stpOp returns a 1 or a 0.
      Teuchos::RCP<astNode<ScalarT> > zeroAst_ = Teuchos::rcp(new numval<ScalarT>(0.0));
      bpTimes_.clear();
      computeBreakPoint ( this->leftAst_, zeroAst_, timeOpVec_, bpTol_, bpTimes_);

      ScalarT xFixed = this->leftAst_->val();
      Xyce::Util::fixNan(xFixed);  Xyce::Util::fixInf(xFixed);
      return ((std::real(xFixed))>0)?1.0:0.0;
    }

    virtual ScalarT dx (int i) { return 0.0; }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0); }
    }

    virtual bool getBreakPoints(std::vector<Xyce::Util::BreakPoint> & breakPointTimes)
    {
      if(!(bpTimes_.empty()))
      {
        for (int ii=0;ii<bpTimes_.size();ii++)
        {
          breakPointTimes.push_back(bpTimes_[ii]);
        }
      }

      return true;
    }

    virtual void setBreakPointTol(double tol) { bpTol_ = tol; }

    virtual bool getIsComplex () { return false; } // val() only uses real part

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "step function operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "step function operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      /// fix this
      os << "STP";
    }

    virtual bool getIsTreeConstant()
    { return this->leftAst_->getIsTreeConstant(); }

    virtual bool stpType() { return true; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<stpOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<stpOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
      this->leftAst_->accept(visitor, this->leftAst_); 
    }

  private:
    std::vector<Teuchos::RCP<astNode<ScalarT> > > timeOpVec_;
    double bpTol_;
    std::vector<Xyce::Util::BreakPoint> bpTimes_;
};

//-------------------------------------------------------------------------------
// ramp function x if x > 0; 0 otherwise
template <typename ScalarT>
class urampOp : public astNode<ScalarT>
{
  public:
    urampOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val()
    {
      double leftVal = std::real(this->leftAst_->val());
      return ((leftVal)>0)?(leftVal):0.0;
    }

    // ERK check this.  Looks wrong
    virtual ScalarT dx (int i)
    {
      return ((std::real(this->leftAst_->val()))>0)?1.0:0.0;
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      double leftVal = std::real(this->leftAst_->val());
      result = ((leftVal)>0)?(leftVal):0.0;

      int numDerivs = derivs.size();
      for (int ii=0;ii<numDerivs;ii++)
      {
        derivs[ii] = ((leftVal)>0)?1.0:0.0;
      }
    }

    virtual bool getIsComplex () { return false; } // val() only uses real part

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "uramp operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "uramp operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "(((std::real(";
      this->leftAst_->codeGen(os);
      os << "))>0)?(std::real(";
      this->leftAst_->codeGen(os);
      os << ")):0.0)";
    }

    virtual bool getIsTreeConstant()
    { return this->leftAst_->getIsTreeConstant(); }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<urampOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<urampOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
      this->leftAst_->accept(visitor, this->leftAst_); 
    }
};

inline bool isLeftCurlyBrace(char c) { return (c=='{'); }

inline bool isLeftParen(char c) { return (c=='('); }

inline bool isQuoteSymbol(char c) { return (c=='"'); }


//-----------------------------------------------------------------------------
// Function      : checkAndFixPulse
//
// Purpose       : This function checks for and if possible, fixes various
//                 issues with input pulse.
//
//                 If it finds duplicate points, it uses only the first.
//
//                 If the times are not monotonically increasing, it exits with
//                 a fatal error.
//
//                 If the removeZeros boolean has been set to true, then it
//                 will exclude any zero points in the pulse.  This is mainly
//                 a problem for splined pulses.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 8/23/2017
//-----------------------------------------------------------------------------
template <typename ScalarT>
inline void checkAndFixPulse (
    bool removeZeros,
    std::vector<ScalarT> & timeVec,
    std::vector<ScalarT> & pulseVec)
{
  std::vector<ScalarT> tmpTime;
  std::vector<ScalarT> tmpPulse;

  int size = timeVec.size();

  bool badPointsRemoved=false;
  bool monotonic=true;
  if (size > 0)
  {
    tmpTime.reserve(timeVec.size());
    tmpPulse.reserve(timeVec.size());

    for (int i=1;i<size;++i)
    {
      if (std::real(timeVec[i]) < std::real(timeVec[i-1])) { monotonic=false; }
    }

    if (!monotonic)
    {
      std::vector<std::string> errStr(1,std::string("Time points in the pulse file are not monotonically increasing"));
      yyerror(errStr);
    }
    else
    {
      int firstNonzero=0;
      if (removeZeros)
      {
        for (int ii=0;ii<timeVec.size();ii++)
        {
          if (std::real(pulseVec[ii]) > 0.0)
          {
            tmpTime.push_back(timeVec[ii]);
            tmpPulse.push_back(pulseVec[ii]);
            firstNonzero=ii;
            break;
          }
          else { badPointsRemoved=true; }
        }
      }
      else
      {
        tmpTime.push_back(timeVec[0]);
        tmpPulse.push_back(pulseVec[0]);
      }

      for (int i=(firstNonzero+1);i<size;++i)
      {
        if (std::real(timeVec[i]) != std::real(timeVec[i-1])
           &&
          (std::real(pulseVec[i]) > 0.0)
            )
        {
          tmpTime.push_back(timeVec[i]);
          tmpPulse.push_back(pulseVec[i]);
        }
        else { badPointsRemoved=true; }
      }
    }
  }

  if (badPointsRemoved)
  {
    timeVec.clear();
    pulseVec.clear();
    timeVec = tmpTime;
    pulseVec = tmpPulse;
  }

  if (timeVec.size() < 1)
  {
    std::vector<std::string> errStr(1,std::string("After fixes, the specified pulse has size < 1, which is not valid."));
    yyerror(errStr);
  }
}

//-----------------------------------------------------------------------------
// Removes pulse points that are Nan or Inf (usually done after spline)
//-----------------------------------------------------------------------------
template <typename ScalarT>
inline void removePulseNaNs (
    std::vector<ScalarT> & timeVec,
    std::vector<ScalarT> & pulseVec)
{
  std::vector<ScalarT> tmpTime;
  std::vector<ScalarT> tmpPulse;
  int size = timeVec.size();

  bool badFound=false;
  tmpTime.clear();
  tmpPulse.clear();

  for (int i=0;i<size;++i)
  {
    if (  !(std::isnan(std::real( pulseVec[i])))
       && !(std::isinf(std::real( pulseVec[i])))
        )
    {
      badFound=true;
      tmpTime.push_back(timeVec[i]);
      tmpPulse.push_back(pulseVec[i]);
    }
  }

  if (badFound)
  {
    pulseVec.clear();
    timeVec = tmpTime;
    pulseVec = tmpPulse;
  }

  if (timeVec.size() < 1)
  {
    std::vector<std::string> errStr(1,std::string("After fixes, the specified pulse has size < 1, which is not valid."));
    yyerror(errStr);
  }
}

//-----------------------------------------------------------------------------
// Function      : trapezoidIntegral
// Purpose       :
// Special Notes : used by downSampleVector to normalize the pulse
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename ScalarT>
inline void trapezoidIntegral (
   const std::vector<ScalarT> & times,
   const std::vector<ScalarT> & values,
   std::vector<ScalarT> & testIntegral,
   ScalarT & integral)
{
  int cpSize = times.size();
  int midIndex = cpSize-1;
  integral=0.0;

  testIntegral.resize(cpSize,0.0);

  for (int is=0;is<cpSize-1;++is)
  {
    ScalarT deltaT = times[is+1]-times[is];
    ScalarT pulse1 = values[is];
    ScalarT pulse2 = values[is+1];
    ScalarT Tau1 = times[is];
    ScalarT Tau2 = times[is+1];
    ScalarT deltaI = 0.5*(pulse1+pulse2)*deltaT;
    integral += deltaI;
    testIntegral[is+1] = integral;
  }
}

//-------------------------------------------------------------------------------
// Function      : downSampleVector
// Purpose       : reduce size of a table to a specified length.
// Special Notes :
//
// Reducing the size of a table (or pulse file) can significantly reduce runtime,
// partly because the default behavior of tables is to create a breakpoint for
// every point in the table.
//
// This is only used for tables that have been specified via files, partly
// because that use case was much easier to parse, and partly because most
// really large tables are specified that way.
//
// This function does time sampling based on the ideas in these papers:
//
// Yiming Zhao and Panagiotis Tsiotras.
// "Mesh Refinement Using Density Function for Solving Optimal Control Problems",
// AIAA Infotech@Aerospace Conference, Infotech@Aerospace Conferences, ()
// http://dx.doi.org/10.2514/6.2009-2019
//
// Yiming Zhao and Panagiotis Tsiotras.
// "Density Functions for Mesh Refinement in Numerical Optimal Control",
// Journal of Guidance, Control, and Dynamics, Vol. 34, No. 1 (2011), pp. 271-277.
// http://dx.doi.org/10.2514/1.45852
//
// The general idea is very simple:
//
// (1) Using the spline functions obtain the derivative array of the pulse file.
//     Call it dfdt, where "f" is the values in the pulse file, and "t" is time.
// (2) Take absolute value of dfdt, so always >= 0.  Let dfdt be the "density function"
//     to be used in refining the mesh.
// (3) Create a new array that is the integral of dfdt.  Kind of like a CDF.
//     Call it F, where F(t_i) = \int_t_0^t_i dfdt dt
// (4) since dfdt is always >=0, the array/function F(t) is monotonically increasing.
// (5) Min of F is at F(t_0), Max of F is at F(t_N) where N=number of sample points.
// (6) Let the integral from F(t_i) to F(t_{i+1}) be a constant DF = (1/N)* F(t_N)
// (7) Set up a spline function, but with the purpose to provide F^{-1}.
//     In other words, let the F array be "x" and the time array be "y".
//     That way, for any given F, the spline will return time.
// (8) Loop i from 1 to N.   At each i, F = i*DF.  Use the F^{-1} spline to obtain time[i].
//     At the end of the loop, the gradient-based time points will be set up.
//
//  This can be done with either the original pulse, or with the log10
//  of the original pulse.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 1/27/2023
//-------------------------------------------------------------------------------
template <typename ScalarT>
inline void downSampleVector(
    const int numSamples,
    bool logSpacing,
    std::vector<ScalarT> & ta_,
    std::vector<ScalarT> & ya_
    )
{
  std::vector<ScalarT> pulseTimes = ta_;
  std::vector<ScalarT> pulseValues = ya_;

  std::vector<ScalarT> pulseValuesLog;
  std::vector<ScalarT> pulseIntegral;

  bool removeZeros=true;
  checkAndFixPulse (removeZeros, pulseTimes, pulseValues);

  // prepare to normalize the reduced data set by integrating original pulse.
  ScalarT originalIntegral=0.0;
  trapezoidIntegral( pulseTimes, pulseValues, pulseIntegral, originalIntegral);

  // redo the pulse times and the pulse values, by splining the
  // pulse data file.
  int pulseSize = pulseTimes.size();

  Xyce::Util::interpolator<ScalarT> * pulseInterpolator;
  pulseInterpolator = new Xyce::Util::akima<ScalarT>();
  pulseInterpolator->clear();
  pulseInterpolator->init(pulseTimes, pulseValues);

  ScalarT tstart = pulseTimes[0];
  ScalarT tstop = pulseTimes[pulseSize-1];


  Xyce::Util::interpolator<ScalarT> * pulseInterpolatorLog;
  if ( logSpacing )
  {
    pulseValuesLog.resize(pulseSize,0.0);
    for (int ip=0;ip<pulseSize;++ip)
    {
      pulseValuesLog[ip] = log10(pulseValues[ip]);
    }

    checkAndFixPulse (removeZeros, pulseTimes, pulseValuesLog);

    if (pulseValuesLog.size() != pulseValuesLog.size())
    {
      pulseSize = pulseValuesLog.size();
      pulseValues.clear();
      pulseValues.resize(pulseValuesLog.size());
      for (int ii=0;ii<pulseSize;ii++)
      {
        pulseValues[ii] = std::pow(10.0,pulseValuesLog[ii]);
      }
    }

    pulseInterpolatorLog = new Xyce::Util::akima<ScalarT>();
    pulseInterpolatorLog->clear();
    pulseInterpolatorLog->init(pulseTimes, pulseValuesLog);
  }

  std::vector<ScalarT> pulseValues_dfdtFabs(pulseSize,0.0);
  std::vector<ScalarT> pulseValues_F(pulseSize,0.0);

  for (int ic=0;ic<pulseSize;++ic)
  {
    ScalarT tmp=0.0;
    ScalarT time=pulseTimes[ic];

    if ( logSpacing )
    {
      pulseInterpolatorLog->evalDeriv( pulseTimes, pulseValuesLog, time, tmp);
    }
    else
    {
      pulseInterpolator->evalDeriv( pulseTimes, pulseValues, time, tmp);
    }

    pulseValues_dfdtFabs[ic] = fabs(std::real(tmp));
  }

  Xyce::Util::akima<ScalarT> newInterp;
  newInterp.clear();
  newInterp.init(pulseTimes, pulseValues_dfdtFabs);

  for (int ic=0;ic<pulseSize;++ic)
  {
    ScalarT tmp=0.0;
    ScalarT time0=tstart;
    ScalarT time=pulseTimes[ic];
    newInterp.evalInteg( pulseTimes, pulseValues_dfdtFabs, time0, time, tmp);
    pulseValues_F[ic] = tmp;
  }
  ScalarT F_N = pulseValues_F[pulseSize-1];
  ScalarT one_over_F_N = 1.0/pulseValues_F[pulseSize-1];

  // normalize F array:
  std::for_each(pulseValues_F.begin(), pulseValues_F.end(), [&one_over_F_N](ScalarT &el) { el *= one_over_F_N; } );

  Xyce::Util::akima<ScalarT> inverseFInterp;
  inverseFInterp.clear();
  inverseFInterp.init( pulseValues_F, pulseTimes);

  // Create a sparsified mesh based on gradients.
  ScalarT dF=1.0/static_cast<ScalarT>(numSamples-1);
  ScalarT F=0.0;

  ta_.clear();
  ya_.clear();
  ta_.resize(numSamples,0.0);
  ya_.resize(numSamples,0.0);

  for (int ic=0;ic<numSamples;++ic)
  {
    ScalarT t=0.0;
    inverseFInterp.eval(  pulseValues_F, pulseTimes, F, t);

    F += dF;
    ta_[ic]=t;

    ScalarT v=0.0;
    pulseInterpolator->eval( pulseTimes, pulseValues, t, v);
    ya_[ic]=v;
  }

  removePulseNaNs ( ta_, ya_ );

  int size = ya_.size();

  if ( logSpacing )
  {
    pulseInterpolatorLog->clear();
    delete pulseInterpolatorLog;
  }

  // complete the normalization:
  ScalarT newIntegral=0.0;
  std::vector<ScalarT> testIntegral;
  trapezoidIntegral(ta_,
    ya_,
    testIntegral,
    newIntegral);

  if (std::real(newIntegral)!=0.0)
  {
    ScalarT normalization=originalIntegral/newIntegral;

    int size = ya_.size();
    for (int ic=0;ic<size;ic++)
    {
      ya_[ic] *= normalization;
    }
  }

  pulseInterpolator->clear();
  delete pulseInterpolator;

  pulseTimes.clear();
  pulseValues.clear();
}

//-------------------------------------------------------------------------------
// This is an interpolation operator.
//
// Arbitrary number of (y,z) pairs can be specified
//
// This class supports the following:
//
// TABLE(x,y,z,*) = PWL interpolation
// SPLINE(x,y,z,*)or AKIMA(x,y,z,*)  = akima aplines
// CUBIC(x,y,z,*)  = traditional cubic spline
// WODICKA(x,y,z,*) = wodicka spline
// BLI(x,y,z,*) = Barycentric Lagrange Interpolation
//
// f(x) where f(y) = z
//-------------------------------------------------------------------------------
template <typename ScalarT>
class tableOp : public astNode<ScalarT>
{
  public:
    //-------------------------------------------------------------------------------
    // functions:
    tableOp (const std::string & kw, Teuchos::RCP<astNode<ScalarT> > &input, std::vector<Teuchos::RCP<astNode<ScalarT> > > & args):
      astNode<ScalarT>(), tableArgs_(args),
      allConst_(true),
      input_(input),
      useBreakPoints_(true),
      keyword_(kw)
      {
        allocateInterpolators();

        int size = tableArgs_.size();
        if (size % 2)
        {
          std::vector<std::string> errStr(1,std::string("AST node (table) needs an even number of arguments")); yyerror(errStr);
        }
        else
        {
          allConst_=true; ta_.resize(size/2); ya_.resize(size/2); dya_.resize(size/2,0.0);
          evaluatedAndConstant_.resize(size,0); // cannot know status yet, so initialize to 0.
          for (int ii=0,jj=0;ii<size;ii+=2,jj++)
          {
            ta_[jj] = (tableArgs_)[ii]->val();
            ya_[jj] = (tableArgs_)[ii+1]->val();
            if (!( (tableArgs_)[ii]->numvalType() && (tableArgs_)[ii+1]->numvalType() ) ) { allConst_ = false; }
          }
          yInterpolator_->init(ta_,ya_); // for linear, this isn't necessary, but for others it is

          if (ya_.size() > 2 && ( keyword_==std::string("TABLE") || keyword_==std::string("FASTTABLE") ) )
          {
            createOldStyleDerivativeTable ();
          }
        }
      };

    //-------------------------------------------------------------------------------
    // special constructor for values read in from a file, in which the file IO is
    // handled directly in this constructor.
    tableOp (
        const std::string & kw,
        Teuchos::RCP<astNode<ScalarT> > &input,
        const std::string & filename,
        Teuchos::RCP<astNode<ScalarT> > & numSamplesAst,
        Teuchos::RCP<astNode<ScalarT> > & logSpacingAst
        ):
      astNode<ScalarT>(), allConst_(true), input_(input),
      useBreakPoints_(true),
      keyword_(kw)
      {
        int numDownSamples =
          static_cast<int>( std::real(numSamplesAst->val()));
        bool logSpacing = ((static_cast<int>( std::real(logSpacingAst->val())))>0);

        allocateInterpolators();

        if ( !(Xyce::Util::checkIfValidFile(filename)) )
        {
          std::vector<std::string> errStr(1,std::string("Could not find file " + filename));
          yyerror(errStr);
        }
        else
        {
          std::ifstream dataIn;
          dataIn.open(filename.c_str(), std::ios::in);
          if ( !dataIn.good() )
          {
            std::vector<std::string> errStr(1,std::string("Could not open file " + filename));
            yyerror(errStr);
          }
          else
          {
            double time;
            double value;
            ta_.clear();
            ya_.clear();
            std::string lineOfFile;
            while(!std::getline(dataIn, lineOfFile).eof())
            {
              // eliminate comments (anything to the right of the # symbol, inclusive)
              std::size_t found = lineOfFile.find_first_of("#");
              if (found!=std::string::npos)
              {
                lineOfFile.erase(found);
              }

              // eliminate comments (anything to the right of the * symbol, inclusive)
              found = lineOfFile.find_first_of("*");
              if (found!=std::string::npos)
              {
                lineOfFile.erase(found);
              }

              // eliminate comments (anything to the right of the ; symbol, inclusive)
              found = lineOfFile.find_first_of(";");
              if (found!=std::string::npos)
              {
                lineOfFile.erase(found);
              }

              if (!lineOfFile.empty())
              {
                std::istringstream iss (lineOfFile,std::istringstream::in);
                std::vector<double> rowDoubles;

                while ( iss >> value )
                {
                  rowDoubles.push_back(value);
                }
                if (rowDoubles.size()==2)
                {
                  ta_.push_back(rowDoubles[0]);
                  ya_.push_back(rowDoubles[1]);
                  rowDoubles.clear();
                }
                else
                {
                  std::vector<std::string> errStr(1,std::string("Data file format error. " + filename));
                  yyerror(errStr);
                }
              }
            }
          }
          dataIn.close();

          // now downsample if necessary
          if (numDownSamples>0 )
          {
            downSampleVector(numDownSamples, logSpacing, ta_, ya_);
          } // downsample if-statement
        }

        evaluatedAndConstant_.resize(2*(ta_.size()),1); // tables from files are always pure numbers, so initialize to 1

        yInterpolator_->init(ta_,ya_); // for linear, this isn't necessary, but for others it is

        if (ya_.size() > 2 && ( keyword_==std::string("TABLE") || keyword_==std::string("FASTTABLE") ) )
        {
          createOldStyleDerivativeTable ();
        }
      };

    //-------------------------------------------------------------------------------
    // special constructor for values read in from a file, that are now stored in std::vector objects
    // ERK.  Currently, Xyce doesn't use this function, but it should, as it is more reliable than
    // the above numvalType test in the first constructor.
    tableOp (const std::string & kw, Teuchos::RCP<astNode<ScalarT> > & input, const std::vector<ScalarT> & xvals, const std::vector<ScalarT> & yvals):
      astNode<ScalarT>(), allConst_(true), input_(input),
      useBreakPoints_(true),
      keyword_(kw)
      {
        allocateInterpolators();

        allConst_=true; int size = xvals.size(); int size2=yvals.size();
        if (size != size2)
        {
          std::vector<std::string> errStr(1,std::string("AST node (table) needs x and y vectors to be the same size.")); yyerror(errStr);
        }
        ta_.resize(size); ya_.resize(size); dya_.resize(size,0.0);
        evaluatedAndConstant_.resize(2*size,1); // tables from files are always pure numbers, so initialize to 1

        for (int ii=0;ii<size;ii++)
        {
          ta_[ii] = xvals[ii];
          ya_[ii] = yvals[ii];
        }

        yInterpolator_->init(ta_,ya_); // for linear, this isn't necessary, but for others it is

        if (ya_.size() > 2 && ( keyword_==std::string("TABLE") || keyword_==std::string("FASTTABLE") ) )
        {
          createOldStyleDerivativeTable ();
        }
      };

    //-------------------------------------------------------------------------------
    // The interpolation functions use std::vector<double> objects.  However, the
    // table was first set up as an array of AST nodes, which could be pure numbers,
    // or could be expressions.  This function re-populates the std::vector<double>
    // objects by evaluating the AST node expressions.  The expectation is that these
    // expressions will NOT depend on solution variables (voltages or currents) and
    // also will not depend on time.  So, for table entries that aren't numbers,
    // they should just depend on .params and/or .global_params.
    //
    // This function should only be called when parameter values can change,
    // such as .STEP iterations.  For large tables with many expressions, it can
    // be expensive to set up.
    //
    // This function also attempts to be efficient and only update entries which
    // are non-constant.  So, numerical entries should not be re-evaluated, and
    // constant expressions should not be re-evaluated.  This aspect of the function
    // may be more trouble than it is worth.
    //-------------------------------------------------------------------------------
    virtual bool updateForStep ()
    {
      bool updated=false;

      // if the table was specified as pure numbers, then the allConst_ flag was "true"
      // in the constructor.  If the table contains expressions, then it might turn out to
      // be allConst_, or it might not.  But we can't know for sure at construction if
      // some expression-based entries are const or not.  If they depend on .params, then
      // that isn't known until later.
      //
      if (!allConst_)  // if not all constants, then might need to reinitialize the arrays
      {
        if(!(tableArgs_.empty()))
        {
          bool tmpAllConst=true;
          int size = tableArgs_.size();

          // If the table depends on a .param, then the construction of the table happened
          // before it was determined if the .param was constant or not.  So, the
          // "getIsTreeConstant" status of the tableArg that depends on that parameter
          // may have changed between construction and the first call to val().  If it
          // has changed, it must be re-evaluated at least once, since the function
          // that sets a parameter to be constant (make_constant) also sets the value.
          for (int ii=0,jj=0;ii<size;ii+=2,jj++)
          {
            if ( evaluatedAndConstant_[ii] == 0 )
            {
              bool tIsConstant = (tableArgs_)[ii]->getIsTreeConstant();
              if ( !tIsConstant ) { tmpAllConst = false; }

              evaluatedAndConstant_[ii] = tIsConstant?1:0;
              updated = true;
              ta_[jj] = (tableArgs_)[ii]->val();
            }

            if ( evaluatedAndConstant_[ii+1] == 0)
            {
              bool yIsConstant = (tableArgs_)[ii+1]->getIsTreeConstant();
              if ( !yIsConstant ) { tmpAllConst = false; }

              evaluatedAndConstant_[ii+1] = yIsConstant?1:0;
              updated = true;
              ya_[jj] = (tableArgs_)[ii+1]->val();
            }
          }
          allConst_ = tmpAllConst;

          if ( keyword_!=std::string("TABLE") && keyword_!=std::string("FASTTABLE") && updated )
          {
            yInterpolator_->init(ta_,ya_);
          }

          if (ya_.size() > 2 && ( keyword_==std::string("TABLE") || keyword_==std::string("FASTTABLE") ) && updated )
          {
            createOldStyleDerivativeTable ();
          }
        } // tableArgs_.empty()
      } // !allConst_

      return updated;
    }

    //-------------------------------------------------------------------------------
    void createOldStyleDerivativeTable ()
    {
      // create derivative table
      //
      // this code mimics the old expression library.   It uses finite differencing
      // to set up a new table of derivatives.  The new table is based on the midpoints of
      // the original table, so it has one extra entry.
      //
      // I initially tried to use the evalDeriv function in the yInterpolator object.
      // That method doen't use midpoints, it just differentiates the the linear
      // interpolation device.  That approach failed at least one regression test.
      //
      // This should only be called if using the linear option, i.e. TABLE or TABLEFILE
      int ya_size = ya_.size();
      ta2_.resize(ya_size+1);
      dya_.resize(ya_size+1);
      ta2_[0] = ta_[0]; ta2_[ya_size] = ta_[ya_size-1];
      dya_[0] = 0.0;    dya_[ya_size] = 0.0;
      for (int ii=1;ii<ya_size;++ii)
      {
        ta2_[ii] = 0.5* (ta_[ii-1]+ta_[ii]);
        ScalarT h = ( ta_[ii]- ta_[ii-1]);
        if (std::real(h) != 0.0)
        {
          dya_[ii] = ( ya_[ii]- ya_[ii-1])/ h;
        }
        else
        {
          dya_[ii] = 0.0;
        }
      }
      dyInterpolator_->init(ta2_,dya_); // for linear, this isn't necessary, but for others it is

      return;
    }

    //-------------------------------------------------------------------------------
    void allocateInterpolators()
    {
      // remove whitespace from the keyword
      keyword_.erase(std::remove_if(keyword_.begin(),keyword_.end(), ::isspace), keyword_.end());

      // remove left curly braces from the keyword
      keyword_.erase(std::remove_if(keyword_.begin(),keyword_.end(), isLeftCurlyBrace), keyword_.end());

      // remove left parens from the keyword
      keyword_.erase(std::remove_if(keyword_.begin(),keyword_.end(), isLeftParen), keyword_.end());

      // remove quote symbol from the keyword
      keyword_.erase(std::remove_if(keyword_.begin(),keyword_.end(), isQuoteSymbol), keyword_.end());

      // remove the word "FILE" from the keyword
      std::string fileString ("FILE");
      std::string::size_type n = fileString.length();
      for (std::string::size_type i = keyword_.find(fileString); i != std::string::npos; i = keyword_.find(fileString))
      {
        keyword_.erase(i, n);
      }

      if (keyword_==std::string("FASTTABLE"))
      {
        useBreakPoints_ = false;
        yInterpolator_  = Teuchos::RCP<Xyce::Util::interpolator<ScalarT> >(new Xyce::Util::linear<ScalarT>());
        dyInterpolator_ = Teuchos::RCP<Xyce::Util::interpolator<ScalarT> >(new Xyce::Util::linear<ScalarT>());
      }
      else if (keyword_==std::string("TABLE"))
      {
        useBreakPoints_ = true;
        yInterpolator_  = Teuchos::RCP<Xyce::Util::interpolator<ScalarT> >(new Xyce::Util::linear<ScalarT>());
        dyInterpolator_ = Teuchos::RCP<Xyce::Util::interpolator<ScalarT> >(new Xyce::Util::linear<ScalarT>());
      }
      else if (keyword_==std::string("CUBIC"))
      {
        useBreakPoints_ = false;
        yInterpolator_  = Teuchos::RCP<Xyce::Util::interpolator<ScalarT> >(new Xyce::Util::cubicSpline<ScalarT>());
        dyInterpolator_ = Teuchos::RCP<Xyce::Util::interpolator<ScalarT> >(new Xyce::Util::cubicSpline<ScalarT>());
      }
      else if (keyword_==std::string("SPLINE") || keyword_==std::string("AKIMA"))
      {
        useBreakPoints_ = false;
        yInterpolator_  = Teuchos::RCP<Xyce::Util::interpolator<ScalarT> >(new Xyce::Util::akima<ScalarT>());
        dyInterpolator_ = Teuchos::RCP<Xyce::Util::interpolator<ScalarT> >(new Xyce::Util::akima<ScalarT>());
      }
      else if (keyword_==std::string("WODICKA"))
      {
        useBreakPoints_ = false;
        yInterpolator_  = Teuchos::RCP<Xyce::Util::interpolator<ScalarT> >(new Xyce::Util::wodicka<ScalarT>());
        dyInterpolator_ = Teuchos::RCP<Xyce::Util::interpolator<ScalarT> >(new Xyce::Util::wodicka<ScalarT>());
      }
      else if (keyword_==std::string("BLI")) // Barycentric Lagrange Interpolation
      {
        useBreakPoints_ = false;
        yInterpolator_  = Teuchos::RCP<Xyce::Util::interpolator<ScalarT> >(new Xyce::Util::barycentricLagrange<ScalarT>());
        dyInterpolator_ = Teuchos::RCP<Xyce::Util::interpolator<ScalarT> >(new Xyce::Util::barycentricLagrange<ScalarT>());
      }
      else
      {
        std::vector<std::string> errStr(1,std::string("AST node (table) type not recognized.  Type = "));
        errStr[0] += keyword_;
        yyerror(errStr);
      }
    }

    //-------------------------------------------------------------------------------
    virtual ScalarT val()
    {
      ScalarT y = 0.0;

      if ( !(ta_.empty()) )
      {
        ScalarT input = std::real(this->input_->val());
        int arraySize=ta_.size();
        if (std::real(input) < std::real(ta_[0]))
        {
          y = ya_[0];
        }
        else if (std::real(input) > std::real(ta_[arraySize-1]))
        {
          y = ya_[arraySize-1];
        }
        else
        {
          yInterpolator_->eval(ta_,ya_, input, y);
        }
      }

      return y;
    };

    //-------------------------------------------------------------------------------
    virtual ScalarT dx(int i)
    {
      ScalarT dydx = 0.0;

      if ( ( keyword_==std::string("TABLE") || keyword_==std::string("FASTTABLE") ) )
      {
        dydx = dx_linear(i);
      }
      else
      {
        dydx = dx_splines(i);
      }

      return dydx;
    }

    //-------------------------------------------------------------------------------
    // ERK FIX THIS
    //-------------------------------------------------------------------------------
    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();
      for(int ii=0;ii<derivs.size();ii++)
      {
        derivs[ii] = dx(ii);
      }
    }

    // see the comments in the function "createOldStyleDerivativeTable"
    // for motivation for this func.
    ScalarT dx_linear(int i)
    {
      ScalarT dydx = 0.0;
      ScalarT dinput_dx = std::real(this->input_->dx(i));
      if (std::real(dinput_dx) != 0.0)
      {
        if (!allConst_)  // if not all pure numbers, then initialize the arrays again
        {
          if (!(tableArgs_.empty()))
          {
            int size = tableArgs_.size();
            for (int ii=0,jj=0;ii<size;ii+=2,jj++)
            {
              ta_[jj] = (tableArgs_)[ii]->val();
              ya_[jj] = (tableArgs_)[ii+1]->val();
            }
            yInterpolator_->init(ta_,ya_); // for linear, this isn't necessary, but for others it is

            if (ya_.size() > 2)
            {
              createOldStyleDerivativeTable ();
            }
          }
        }

        ScalarT input = std::real(this->input_->val());

        if ( !(ta2_.empty()) )
        {
          int arraySize=ta2_.size();

          if (std::real(input) <= std::real(ta2_[0]))
          {
            dydx = 0.0;
          }
          else if (std::real(input) >= std::real(ta2_[arraySize-1]))
          {
            dydx = 0.0;
          }
          else
          {
            dyInterpolator_->eval(ta2_,dya_, input, dydx);
            dydx *= dinput_dx;
          }
        }
        else
        {
          dydx=0.0;
          if  (ya_.size()==2)
          {
            if (std::real(input) <= std::real(ya_[1]) && std::real(input) >= std::real(ya_[0]))
            {
              ScalarT h = (ta_[1]-ta_[0]);
              if (std::real(h) != 0.0) { dydx = (ya_[1]-ya_[0])/h;}
              dydx *= dinput_dx;
            }
          }
        }
      }
#if 0
      // this code is slightly busted, due to the changes above. Fix later.
      // (dyInterpolator and dya are used a little differently)
      // This code is for derivatives w.r.t. table values, rather than the input.
      // This is a use case that the old library didn't handle.
      // So, for now, leaving it commented out.
      else
      {
        // derivative w.r.t. table y values
        //
        if (!allConst_)  // if not all pure numbers, then initialize the arrays again.
        {
          for (int ii=0,jj=0;ii<size;ii+=2,jj++)
          {
            ta_[jj] = (tableArgs_)[ii]->val();
            dya_[jj] = (tableArgs_)[ii+1]->dx(i);
          }
          dyInterpolator_->init(ta_,dya_); // for linear, this isn't necessary, but for others it is

          ScalarT input = std::real(this->input_->val());

          if ( !(ta_.empty()) )
          {
            int arraySize=ta_.size();
            if (std::real(input) < std::real(ta_[0]))
            {
              dydx = dya_[0];
            }
            else if (std::real(input) > std::real(ta_[arraySize-1]))
            {
              dydx = dya_[arraySize-1];
            }
            else
            {
              dyInterpolator_->eval(ta_,dya_, input, dydx);
            }
          }
        }
      }
#endif
      return dydx;
    }

    ScalarT dx_splines(int i)
    {
      ScalarT dydx = 0.0;

      ScalarT dinput_dx = std::real(this->input_->dx(i));

      if (std::real(dinput_dx) != 0.0)
      {
        // derivative w.r.t. input, using the "evalDeriv" function
        // The higher-order interpolators should produce smooth derivatives, so this
        // approach is much cleaner than what we do for the linear interpolator.
        if (!allConst_)  // if not all pure numbers, then initialize the arrays again
        {
          if (!(tableArgs_.empty()))
          {
            int size = tableArgs_.size();
            for (int ii=0,jj=0;ii<size;ii+=2,jj++)
            {
              ta_[jj] = (tableArgs_)[ii]->val();
              ya_[jj] = (tableArgs_)[ii+1]->val();
            }
            yInterpolator_->init(ta_,ya_); // for linear, this isn't necessary, but for others it is
          }
        }

        if (!(ta_.empty()))
        {
          ScalarT input = std::real(this->input_->val());
          yInterpolator_->evalDeriv(ta_,ya_, input, dydx);
          dydx *= dinput_dx;
        }

      }
      return dydx;
    }

    // in practice so far, tables are only used for real numbers, mostly data-based sources.
    // in theory, they could be used for complex numbers, but there is no use case for
    // this yet.   Update if that changes.
    //virtual bool getIsComplex () { return (typeid(ScalarT) == typeid(std::complex<double>)) ; }
    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0) // FIX THIS
    {
      os << std::setw(indent) << " ";
      os << "table operator id = " << this->id_ << std::endl;
      //++indent;
      //this->input_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "table operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this
      os << "TABLE";
    }

    virtual bool getBreakPoints(std::vector<Xyce::Util::BreakPoint> & breakPointTimes)
    {
      if (useBreakPoints_)
      {
        if (!(ta_.empty()))
        {
          if ( input_->timeSpecialType() )
          {
            ScalarT time = std::real(this->input_->val());
            size_t size = ta_.size();
            size_t index = yInterpolator_->binarySearch (ta_, time, 0, size - 1);

            if ( std::real(ta_[index]) < std::real(time))
            {
             int tmp=index;
             while(tmp < size && std::real(ta_[tmp]) < std::real(time) ) { tmp++; }
             index = tmp;
            }

            if (index < size)
            {
             size_t max = index+5;
             if ( max  > size ) { max = size; }
             int ii=index;
             for( ;ii<max;ii++)
             {
               breakPointTimes.push_back( std::real(ta_[ii]) );
             }
            }
          }
        }
      }
      return true;
    }

    virtual bool srcType() { return ( input_->timeSpecialType() ); }

    virtual bool getIsTreeConstant() { return allConst_; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<tableOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<tableOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch

      input_->accept(visitor, input_);

      if (!allConst_)
      {
        if (!(tableArgs_.empty()))
        {
          int size=tableArgs_.size();
          for(int ii=0;ii<size;ii++)
          {
            tableArgs_[ii]->accept(visitor, tableArgs_[ii]);
          }
        }
      }
    }

  private:
    std::vector<Teuchos::RCP<astNode<ScalarT> > > tableArgs_;
    bool allConst_;
    std::vector<ScalarT> ta_; // using ta for name instead of xa so as not to confuse meaning of dx function
    std::vector<ScalarT> ya_;
    std::vector<ScalarT> ta2_; // using ta for name instead of xa so as not to confuse meaning of dx function
    std::vector<ScalarT> dya_;
    std::vector<int> evaluatedAndConstant_;

    Teuchos::RCP<Xyce::Util::interpolator<ScalarT> > yInterpolator_;
    Teuchos::RCP<Xyce::Util::interpolator<ScalarT> > dyInterpolator_;

    Teuchos::RCP<astNode<ScalarT> > input_;
    bool useBreakPoints_;
    std::string keyword_;
};

//-------------------------------------------------------------------------------
// SCHEDULE (y,z,*)
// f(time) where f(y) = z
// kind of like a table, but with no interpolation
template <typename ScalarT>
class scheduleOp : public astNode<ScalarT>
{
  public:
    // functions:
    scheduleOp (
        std::vector<Teuchos::RCP<astNode<ScalarT> > > & args,
        Teuchos::RCP<astNode<ScalarT> > &time
        ): astNode<ScalarT>(), time_(time), tableArgs_(args), allNumVal_(true)
      {
        int size = tableArgs_.size();
        if (size % 2)
        {
          std::vector<std::string> errStr(1,std::string("AST node (schedule) needs an even number of arguments")); yyerror(errStr);
        }
        else
        {
          allNumVal_=true; ta_.resize(size/2); ya_.resize(size/2);
          //dya_.resize(size/2,0.0);
          for (int ii=0,jj=0;ii<size;ii+=2,jj++)
          {
            ta_[jj] = (tableArgs_)[ii]->val();
            ya_[jj] = (tableArgs_)[ii+1]->val();
            if (!( (tableArgs_)[ii]->numvalType() && (tableArgs_)[ii+1]->numvalType() ) ) { allNumVal_ = false; }
          }
        }
      };

     //If the schedule is (t0, dt0, t1, dt1, t2, dt2)
     // then the value is:
     // if time < t0      value = 0
     // if t0 < time < t1 value = dt0
     // if t1 < time < t2 value = dt1
     // if t2 < time      value = dt2
    virtual ScalarT val()
    {
      ScalarT y = 0.0;
      int size = tableArgs_.size();

      if (!allNumVal_)  // if not all pure numbers, then initialize the arrays again
      {
        for (int ii=0,jj=0;ii<size;ii+=2,jj++)
        {
          ta_[jj] = (tableArgs_)[ii]->val();
          ya_[jj] = (tableArgs_)[ii+1]->val();
        }
      }

      ScalarT time = std::real(this->time_->val());

      if ( !(ta_.empty()) )
      {
        int arraySize=ta_.size();
        if (std::real(time) < std::real(ta_[0]))
        {
          y = 0.0;
        }
        else if (std::real(time) > std::real(ta_[arraySize-1]))
        {
          y = ya_[arraySize-1];
        }
        else
        {
          for (int ii=0;ii<arraySize-1;ii++)
          {
            if (std::real(time) > std::real(ta_[ii]) && std::real(time) <= std::real(ta_[ii+1]) )
            {
              y = ya_[ii];
              break;
            }
          }
        }
      }

      return y;
    };

    virtual ScalarT dx(int i)
    {
      // ERK.  Not implemented.  Might not be needed
      ScalarT dydx = 0.0;
      return dydx;
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0); }
    }

    // in practice so far, schedule is only used for setting up time windows (i.e. real, not complex)
    // It also isn't generally used in outputs.
    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0) // FIX THIS
    {
      os << std::setw(indent) << " ";
      os << "schedule operator id = " << this->id_ << std::endl;
      //++indent;
      //this->time_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os) // FIX THIS
    {
      os << "schedule operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this
      os << "SCHEDULE ";
    }


    // ERK.  This is not efficient.
    virtual bool getBreakPoints(std::vector<Xyce::Util::BreakPoint> & breakPointTimes)
    {
     if ( time_->timeSpecialType() ) // this should always be true ....
     {
        double time = std::real(this->time_->val());
        int size = tableArgs_.size();
        for (int ii=0,jj=0;ii<size;ii+=2,jj++)
        {
          double bpTime =  std::real( (tableArgs_)[ii]->val() );
          breakPointTimes.push_back( bpTime );
        }
      }
      return true;
    }

    virtual bool scheduleType() { return true; }

    virtual bool getIsTreeConstant() { return allNumVal_; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<scheduleOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<scheduleOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch

      time_->accept(visitor, time_);

      if (!allNumVal_)
      {
        int size=tableArgs_.size();
        for(int ii=0;ii<size;ii++)
        {
          tableArgs_[ii]->accept(visitor, tableArgs_[ii]);
        }
      }
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > time_;
    std::vector<Teuchos::RCP<astNode<ScalarT> > > tableArgs_;
    bool allNumVal_;
    std::vector<ScalarT> ta_; // using ta for name instead of xa so as not to confuse meaning of dx function
    std::vector<ScalarT> ya_;
};

//-------------------------------------------------------------------------------
// time integral of x
template <typename ScalarT>
class sdtOp : public astNode<ScalarT>
{
  public:
    sdtOp(Teuchos::RCP<astNode<ScalarT> > &left,
        Teuchos::RCP<astNode<ScalarT> > &dt,
        Teuchos::RCP<astNode<ScalarT> > &time
        )
      : astNode<ScalarT>(left),
      dt_(dt),
      time_(time)
    {
    };

    virtual ScalarT val()
    {
      if (this->processSuccessfulStepFlag)
      {
        unsigned long int id = this->getSdtState().id;
        if ( this->processSuccessfulStepMap.find(id) == this->processSuccessfulStepMap.end() )
        {
          this->getSdtState().processSuccessfulTimeStep ();
          this->processSuccessfulStepMap[id] = 1;
        }
      }

      ScalarT time = 0.0;
      ScalarT deltaT = 0.0;

      if( !(Teuchos::is_null( time_ ))) { time = std::real(this->time_->val()); }
      else
      {
        std::vector<std::string> errStr(1,std::string("AST node (sdt) has a null time pointer"));
        yyerror(errStr);
      }

      if (time != 0.0) // at time point zero, treat dt as zero.
      {
        if( !(Teuchos::is_null( dt_ ))) { deltaT = std::real(this->dt_->val()); }
        else
        {
          std::vector<std::string> errStr(1,std::string("AST node (sdt) has a null dt pointer"));
          yyerror(errStr);
        }
      }

      sdtStateData<ScalarT> & state = this->getSdtState();
      state.val2 = this->leftAst_->val();
      ScalarT deltaI = 0.5*(state.val1+state.val2)*deltaT;
      state.integral = state.integral_old + deltaI;
      return state.integral;
    };

    virtual ScalarT dx(int i)
    {
      ScalarT time = 0.0;
      double deltaT = 0.0;

      if( !(Teuchos::is_null( time_ ))) { time = std::real(this->time_->val()); }
      else
      {
        std::vector<std::string> errStr(1,std::string("AST node (sdt) has a null time pointer"));
        yyerror(errStr);
      }

      if (time != 0.0) // at time point zero, treat dt as zero.
      {
        if( !(Teuchos::is_null( dt_ ))) { deltaT = std::real(this->dt_->val()); }
        else
        {
          std::vector<std::string> errStr(1,std::string("AST node (sdt) has a null dt pointer"));
          yyerror(errStr);
        }
      }

      ScalarT dVal2dx = this->leftAst_->dx(i);
      ScalarT dIdVal2 = 0.5*deltaT;
      ScalarT dIdx = dIdVal2*dVal2dx;
      return dIdx;
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      if (this->processSuccessfulStepFlag)
      {
        unsigned long int id = this->getSdtState().id;
        if ( this->processSuccessfulStepMap.find(id) == this->processSuccessfulStepMap.end() )
        {
          this->getSdtState().processSuccessfulTimeStep ();
          this->processSuccessfulStepMap[id] = 1;
        }
      }

      ScalarT time = 0.0;
      ScalarT deltaT = 0.0;

      if( !(Teuchos::is_null( time_ ))) { time = std::real(this->time_->val()); }
      else
      {
        std::vector<std::string> errStr(1,std::string("AST node (sdt) has a null time pointer"));
        yyerror(errStr);
      }

      if (time != 0.0) // at time point zero, treat dt as zero.
      {
        if( !(Teuchos::is_null( dt_ ))) { deltaT = std::real(this->dt_->val()); }
        else
        {
          std::vector<std::string> errStr(1,std::string("AST node (sdt) has a null dt pointer"));
          yyerror(errStr);
        }
      }

      sdtStateData<ScalarT> & state = this->getSdtState();
      int numDerivs = derivs.size();
      this->leftAst_->dx2(state.val2, derivs);

      ScalarT dIdVal2 = 0.5*deltaT;
      ScalarT deltaI = (state.val1+state.val2)*dIdVal2;
      state.integral = state.integral_old + deltaI;
      result = state.integral;

      for(int ii=0;ii<numDerivs;ii++) { derivs[ii] *= dIdVal2; }
    }

    // in practice this is (so far) only used in transient simulation, so real valued, not complex.
    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "sdt (time integral) operator " << " id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "sdt (time integral) operator " << " id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this
      os << "SDT";
    }

    virtual bool sdtType() { return true; }

    Teuchos::RCP<astNode<ScalarT> > & getArg() { return (this->leftAst_); }

    virtual bool getIsTreeConstant() { return false; }  // time dependent can't be constant

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<sdtOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<sdtOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
      this->leftAst_->accept(visitor, this->leftAst_);
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > dt_;
    Teuchos::RCP<astNode<ScalarT> > time_;
};

//-------------------------------------------------------------------------------
// time derivative of x
//
// This class can compute a Backward Euler time derivative locally in this class,
// or it can rely on an externally supplied time derivative.
template <typename ScalarT>
class ddtOp : public astNode<ScalarT>
{
  public:
    ddtOp(Teuchos::RCP<astNode<ScalarT> > &left,
        Teuchos::RCP<astNode<ScalarT> > &dt,
        Teuchos::RCP<astNode<ScalarT> > &time
        )
      : astNode<ScalarT>(left),
      dt_(dt),
      time_(time),
      timeDerivative_(0.0),
      useExternDeriv_(false)
    {};

    virtual ScalarT val()
    {
      if (this->processSuccessfulStepFlag)
      {
        unsigned long int id = this->getDdtState().id;
        if ( this->processSuccessfulStepMap.find(id) == this->processSuccessfulStepMap.end() )
        {
          this->getDdtState().processSuccessfulTimeStep ();
          this->processSuccessfulStepMap[id] = 1;
        }
      }

      ScalarT time = 0.0;
      ScalarT deltaT = 0.0;
      ddtStateData<ScalarT> & state = this->getDdtState();
      state.val2 = this->leftAst_->val();

      if (!useExternDeriv_ )
      {
        timeDerivative_ = 0.0;

        if( !(Teuchos::is_null( time_ ))) { time = std::real(this->time_->val()); }
        else
        { std::vector<std::string> errStr(1,std::string("AST node (ddt) has a null time pointer")); yyerror(errStr); }

        if (time != 0.0) // at time point zero, treat dt as zero.
        {
          if( !(Teuchos::is_null( dt_ ))) { deltaT = std::real(this->dt_->val()); }
          else
          { std::vector<std::string> errStr(1,std::string("AST node (ddt) has a null dt pointer")); yyerror(errStr); }

          // for now, hardwire to backward Euler
          timeDerivative_ = (state.val2-state.val1)/deltaT;
        }
      }
      return timeDerivative_;
    };

    virtual ScalarT dx(int i) // ERK.  this isn't right for the version that uses setDdtDeriv, as it is hardwired to BE.
    {
      ScalarT ddt_dx = 0.0;
      ScalarT time = 0.0;
      ScalarT deltaT = 0.0;

      if (!useExternDeriv_ )
      {
        if( !(Teuchos::is_null( time_ ))) { time = std::real(this->time_->val()); }
        else
        { std::vector<std::string> errStr(1,std::string("AST node (ddt) has a null time pointer")); yyerror(errStr); }

        if (time != 0.0) // at time point zero, treat dt as zero.
        {
          if( !(Teuchos::is_null( dt_ ))) { deltaT = std::real(this->dt_->val()); }
          else
          { std::vector<std::string> errStr(1,std::string("AST node (ddt) has a null dt pointer")); yyerror(errStr); }


          ScalarT dVal2dx = this->leftAst_->dx(i);
          ScalarT ddt_dVal2 = 1.0/deltaT;
          // for now, hardwire to backward Euler
          ddt_dx = ddt_dVal2 * dVal2dx;
        }
      }
      return ddt_dx;
    };

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      if (this->processSuccessfulStepFlag)
      {
        unsigned long int id = this->getDdtState().id;
        if ( this->processSuccessfulStepMap.find(id) == this->processSuccessfulStepMap.end() )
        {
          this->getDdtState().processSuccessfulTimeStep ();
          this->processSuccessfulStepMap[id] = 1;
        }
      }

      int numDerivs =  derivs.size();
      ScalarT time = 0.0;
      ScalarT deltaT = 0.0;
      ddtStateData<ScalarT> & state = this->getDdtState();

      if (!useExternDeriv_ )
      {
        timeDerivative_ = 0.0;

        if( !(Teuchos::is_null( time_ ))) { time = std::real(this->time_->val()); }
        else
        { std::vector<std::string> errStr(1,std::string("AST node (ddt) has a null time pointer")); yyerror(errStr); }

        if (time != 0.0) // at time point zero, treat dt as zero.
        {
          if( !(Teuchos::is_null( dt_ ))) { deltaT = std::real(this->dt_->val()); }
          else
          { std::vector<std::string> errStr(1,std::string("AST node (ddt) has a null dt pointer")); yyerror(errStr); }

          // for now, hardwire to backward Euler
          ScalarT ddt_dVal2 = 1.0/deltaT;
          this->leftAst_->dx2(state.val2, derivs);
          timeDerivative_ = (state.val2-state.val1)*ddt_dVal2;

          for(int ii=0;ii<numDerivs;ii++)
          {
            // for now, hardwire to backward Euler
            derivs[ii] *= ddt_dVal2;
          }
        }
      }
      else
      {
        state.val2 = this->leftAst_->val();
        if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0); }
      }
      result = timeDerivative_;
    };

    // in practice this is (so far) only used in transient simulation, so real valued, not complex.
    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "ddt (time derivative) operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "ddt (time derivative) operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this
      os << "DDT";
    }

    virtual bool ddtType() { return true; }

    ScalarT getDdtArg()
    {
      ddtStateData<ScalarT> & state = this->getDdtState();
      return state.val2;
    }

    void    setDdtDeriv(ScalarT deriv) { useExternDeriv_ = true; timeDerivative_ = deriv; };

    Teuchos::RCP<astNode<ScalarT> > & getArg() { return (this->leftAst_); }

    virtual bool getIsTreeConstant() { return false; }  // time dependent can't be constant

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<ddtOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<ddtOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
      this->leftAst_->accept(visitor, this->leftAst_);
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > dt_;
    Teuchos::RCP<astNode<ScalarT> > time_;
    ScalarT timeDerivative_;
    bool useExternDeriv_;
};

//-------------------------------------------------------------------------------
// ddx(f(x),x)
// partial derivative of f (x) with respect to x
template <typename ScalarT>
class ddxOp : public astNode<ScalarT>
{
  public:
    ddxOp (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):
      astNode<ScalarT>(left,right) , foundX_(false) { }

    void resolveArg_()
    {
      // check for params
      //
      // For most uses of "paramType" it is important to exclude function
      // arguments, which are also allocated as parameters, but are usually
      // not used in the same way as parameters.  This is the one use case where
      // function arguments must be included.
      //
      // ERK:
      // Note, as currently written, these are very simple linear searches
      // to match the differential variable with a parameter (or voltage or current)
      //
      // This is fine for expressions with very small numbers of parameters, etc.
      // Not fine for really large containers.  Fix later.
      //
      if ( this->rightAst_->paramType() || this->rightAst_->getFunctionArgType() )
      {
        std::vector<Teuchos::RCP<astNode<ScalarT> > > paramOpVector;
        getParamOpsVisitor<ScalarT> visitor(paramOpVector);
        this->leftAst_->accept(visitor,this->leftAst_);

        // now match the user-specified right hand argument with a parameter inside
        // the function specified by the left-hand argument.  This is a klunky,
        // inefficient search and should be replaced with something better.
        // (same is true below for voltages and currents)
        std::string tmp = this->rightAst_->getName();
        if (!(tmp.empty()))
        {
          Xyce::Util::toUpper(tmp);
          for (int ii=0;ii<paramOpVector.size();ii++)
          {
            std::string tmp2 = paramOpVector[ii]->getName();
            if (!(tmp2.empty()))
            {
              Xyce::Util::toUpper(tmp2);
              if (tmp==tmp2) { foundX_ = true; astNodeX_ = paramOpVector[ii]; break; }
            }
          }
        }
      }
      else if (this->rightAst_->voltageType())
      {
        std::vector<Teuchos::RCP<astNode<ScalarT> > > voltOpVector;
        getVoltageOpsVisitor<ScalarT> visitor(voltOpVector);
        this->leftAst_->accept(visitor,this->leftAst_);

        std::string tmp = this->rightAst_->getName();
        if (!(tmp.empty()))
        {
          Xyce::Util::toUpper(tmp);
          for (int ii=0;ii<voltOpVector.size();ii++)
          {
            std::string tmp2 = voltOpVector[ii]->getName();
            if (!(tmp2.empty()))
            {
              Xyce::Util::toUpper(tmp2);
              if (tmp==tmp2) { foundX_ = true; astNodeX_ = voltOpVector[ii]; break; }
            }
          }
        }
      }
      else if (this->rightAst_->currentType())
      {
        std::vector<Teuchos::RCP<astNode<ScalarT> > > currentOpVector;
        getCurrentOpsVisitor<ScalarT> visitor(currentOpVector);
        this->leftAst_->accept(visitor,this->leftAst_);

        std::string tmp = this->rightAst_->getName();
        if (!(tmp.empty()))
        {
          Xyce::Util::toUpper(tmp);
          for (int ii=0;ii<currentOpVector.size();ii++)
          {
            std::string tmp2 = currentOpVector[ii]->getName();
            if (!(tmp2.empty()))
            {
              Xyce::Util::toUpper(tmp2);
              if (tmp==tmp2) { foundX_ = true; astNodeX_ = currentOpVector[ii]; break; }
            }
          }
        }
      }
      else // unsupported type
      {
        std::vector<std::string> errStr(1,std::string("DDX unsupported type"));
        yyerror(errStr);
      }
    };

    virtual ScalarT val()
    {
      ScalarT ret = 0.0;

      if( !foundX_ ) { resolveArg_(); }

      if (!foundX_ || (Teuchos::is_null( astNodeX_)))
      {
        std::string msg = "DDX argument ";
        std::string tmp;
        // ERK: this block of code needs to be revised.  The set of if-statements, below often fail, so tmp=""
        if (this->rightAst_->paramType() || this->rightAst_->getFunctionArgType())
        {
          tmp = this->rightAst_->getName();
        }
        else if (this->rightAst_->currentType())
        {
          tmp = "I(" + this->rightAst_->getName() + ")";
        }
        else if (this->rightAst_->voltageType())
        {
          std::string name = this->rightAst_->getName();
          tmp = "V(";
          tmp += name;
          tmp += ")";
        }
        msg += tmp + " not resolved";

        std::vector<std::string> errStr(1,msg);
        yyerror(errStr);
      }
      else
      {
        astNodeX_->setDerivIndex(0);
        astNodeX_->setIsVar();
        ret = this->leftAst_->dx(0);
        astNodeX_->unsetDerivIndex();
        astNodeX_->unsetIsVar();
      }

      return ret;
    };

    virtual ScalarT dx(int i)
    {
      //std::vector<std::string> errStr(1,std::string("AST node (ddx) without a dx function"));
      //yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0); }
    }

    // revisit this.  Going with typeid is probably good enough for now.
    virtual bool getIsComplex () { return (typeid(ScalarT) == typeid(std::complex<double>)) ; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "ddx (derivative) operator " << " id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "ddx (derivative) operator " << " id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this
      os << "DDX";
    }

    virtual bool getIsTreeConstant() { return false; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<ddxOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<ddxOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
      this->leftAst_->accept(visitor, this->leftAst_);
      this->rightAst_->accept(visitor, this->rightAst_); 
    }

  private:
    bool foundX_;
    Teuchos::RCP<astNode<ScalarT> > astNodeX_;
};

#include <ast_random.h>

//-------------------------------------------------------------------------------
// specials
//  {"TIME", "TEMP", "VT", "FREQ"}
template <typename ScalarT>
class specialsOp : public astNode<ScalarT>
{
  public:
    specialsOp (std::string typeName) : astNode<ScalarT>(), type_(typeName), value_(0.0), derivIndex_(-1)
  {
    Xyce::Util::toUpper(type_);
  };

    virtual ScalarT val() { return value_; };
    virtual ScalarT dx (int i)
    {
      ScalarT retval = (derivIndex_==i)?1.0:0.0;
      return retval;
    };

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = value_;
      if ( !(derivs.empty() ) )
      {
        std::fill(derivs.begin(),derivs.end(),0.0);
        if(derivIndex_>-1) { derivs[derivIndex_] = 1.0; }
      }
    }

    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << type_ << " operator.  val = " << value_ << " id = " << this->id_ << std::endl;
    }

    virtual void compactOutput(std::ostream & os) { output(os,0); }

    virtual void codeGen (std::ostream & os ) { os << value_; }

    std::string getType() { return type_; }
    void getType(std::string t) { type_ = t; }
    ScalarT getValue() { return value_; }
    void setValue(ScalarT val) { value_ = val; }

    virtual void setDerivIndex(int i) { derivIndex_=i; };
    virtual void unsetDerivIndex() { derivIndex_=-1; };

    virtual bool timeSpecialType() { return (type_ == std::string("TIME")); }
    virtual bool dtSpecialType()   { return (type_ == std::string("DT")); }
    virtual bool tempSpecialType() { return (type_ == std::string("TEMP")); }
    virtual bool vtSpecialType()   { return (type_ == std::string("VT")); }
    virtual bool freqSpecialType() { return (type_ == std::string("FREQ")); }
    virtual bool gminSpecialType() { return (type_ == std::string("GMIN")); }

    virtual bool getIsTreeConstant() { return false; } // sometimes constant, sometimes not, so be conservative

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<specialsOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<specialsOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
    }

  private:
    std::string type_;
    ScalarT value_;
    int derivIndex_;
};

//-------------------------------------------------------------------------------
// constants.  These are obsolete.
template <typename ScalarT>
class piConstOp : public astNode<ScalarT>
{
  public:
    piConstOp (): astNode<ScalarT>() {};

    virtual ScalarT val() { return ScalarT(M_PI); };
    virtual ScalarT dx (int i) { return ScalarT(0.0); };
    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)  { result = ScalarT(M_PI);
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0); } }

    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "pi const operator.  val = " << ScalarT(M_PI) << " id = " << this->id_ << std::endl;
    }

    virtual void compactOutput(std::ostream & os) { output(os,0); }

    virtual void codeGen (std::ostream & os ) { os << ScalarT(M_PI); }

    virtual bool getIsTreeConstant() { return true; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<piConstOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<piConstOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
    }
  private:
};

//-------------------------------------------------------------------------------
// constants
template <typename ScalarT>
class CtoKConstOp : public astNode<ScalarT>
{
  public:
    CtoKConstOp (): astNode<ScalarT>() {};

    virtual ScalarT val() { return ScalarT(CONSTCtoK); };
    virtual ScalarT dx (int i) { return ScalarT(0.0); };
    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)  { result = ScalarT(CONSTCtoK);
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0); } }

    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "CtoK const operator.  val = " << ScalarT(CONSTCtoK) << " id = " << this->id_ << std::endl;
    }

    virtual void compactOutput(std::ostream & os) { output(os,0); }

    virtual void codeGen (std::ostream & os ) { os << ScalarT(CONSTCtoK); }

    virtual bool getIsTreeConstant() { return true; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_)
    { 
      Teuchos::RCP<CtoKConstOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<CtoKConstOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
    }

  private:
};

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
inline void yyerror(std::vector<std::string> & s)
{
  //Xyce::Report::UserError() << "ERROR!!!  in expression " << newExp.getExpressionString() << std::endl;
  for (int i=0;i<s.size();++i)
  {
    //std::cerr << "\t" << s[i] << std::endl;
    Xyce::Report::UserError() << s[i] ;
  }
}

#endif

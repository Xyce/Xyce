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

#ifndef ast_H
#define ast_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <list>
#include <unordered_map>
#include <string>
#include <cmath>
#include <complex>

#include <Teuchos_RCP.hpp>

#include <N_UTL_BreakPoint.h>
#include <N_UTL_Interpolators.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_CheckIfValidFile.h>
#include <N_ERH_Message.h>
#include <N_UTL_HspiceBools.h>
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



template <typename ScalarT>
class astNode;

template <typename ScalarT>
class funcOp;

template <typename ScalarT>
class paramOp;

template <typename ScalarT>
class voltageOp;

template <typename ScalarT>
class currentOp;

template <typename ScalarT>
class sdtOp;

template <typename ScalarT>
class ddtOp;

inline void yyerror(std::vector<std::string> & s);

#define AST_GET_INTERESTING_OPS(PTR) if( !(Teuchos::is_null(PTR)) ) {  \
  if (PTR->paramType()) { ovc.paramOpVector.push_back(PTR); }  \
  if (PTR->funcType())    { ovc.funcOpVector.push_back(PTR); } \
  if (PTR->voltageType()) { ovc.voltOpVector.push_back(PTR); } \
  if (PTR->currentType()) { ovc.currentOpVector.push_back(PTR); } \
  if (PTR->leadCurrentType()) { ovc.leadCurrentOpVector.push_back(PTR); } \
  if (PTR->bsrcCurrentType()) { ovc.bsrcCurrentOpVector.push_back(PTR); } \
  if (PTR->powerType()) { ovc.powerOpVector.push_back(PTR); } \
  if (PTR->internalDeviceVarType()) { ovc.internalDevVarOpVector.push_back(PTR); } \
  if (PTR->dnoNoiseVarType()) { ovc.dnoNoiseDevVarOpVector.push_back(PTR); } \
  if (PTR->dniNoiseVarType()) { ovc.dniNoiseDevVarOpVector.push_back(PTR); } \
  if (PTR->oNoiseType()) { ovc.oNoiseOpVector.push_back(PTR); } \
  if (PTR->iNoiseType()) { ovc.iNoiseOpVector.push_back(PTR); } \
  if (PTR->sdtType()) { ovc.sdtOpVector.push_back(PTR); } \
  if (PTR->ddtType()) { ovc.ddtOpVector.push_back(PTR); } \
  if (PTR->srcType()) { ovc.srcOpVector.push_back(PTR); } \
  if (PTR->stpType()) { ovc.stpOpVector.push_back(PTR); } \
  if (PTR->compType()) { ovc.compOpVector.push_back(PTR); } \
  if (PTR->limitType()) { ovc.limitOpVector.push_back(PTR); } \
  if (PTR->phaseType()) { ovc.phaseOpVector.push_back(PTR); } \
  if (PTR->sparamType()) { ovc.sparamOpVector.push_back(PTR); } \
  if (PTR->yparamType()) { ovc.yparamOpVector.push_back(PTR); } \
  if (PTR->zparamType()) { ovc.zparamOpVector.push_back(PTR); } \
  if (PTR->agaussType()) { ovc.agaussOpVector.push_back(PTR); } \
  if (PTR->gaussType()) { ovc.gaussOpVector.push_back(PTR); } \
  if (PTR->aunifType()) { ovc.aunifOpVector.push_back(PTR); } \
  if (PTR->unifType()) { ovc.unifOpVector.push_back(PTR); } \
  if (PTR->randType()) { ovc.randOpVector.push_back(PTR); } \
  if (PTR->twoArgLimitType()) { ovc.twoArgLimitOpVector.push_back(PTR); } \
  if (PTR->timeSpecialType() || PTR->dtSpecialType()) { ovc.isTimeDependent = true; } \
  if (PTR->tempSpecialType()) { ovc.isTempDependent = true; } \
  if (PTR->vtSpecialType()) { ovc.isVTDependent = true; } \
  if (PTR->freqSpecialType()) { ovc.isFreqDependent = true; } \
  if (PTR->gminSpecialType()) { ovc.isGminDependent = true; } \
  PTR->getInterestingOps(ovc); }

#define AST_GET_STATE_OPS(PTR) if( !(Teuchos::is_null(PTR)) ) {  \
  if (PTR->sdtType()) { ovc.sdtOpVector.push_back(PTR); } \
  if (PTR->ddtType()) { ovc.ddtOpVector.push_back(PTR); } \
  PTR->getStateOps(ovc); }

#define AST_GET_PARAM_OPS(PTR)  if( !(Teuchos::is_null(this->PTR)) ) { if (this->PTR->paramType()) { paramOpVector.push_back(this->PTR); } this->PTR->getParamOps(paramOpVector); }

#define AST_GET_FUNC_ARG_OPS(PTR)  if( !(Teuchos::is_null(this->PTR)) ) { if (this->PTR->getFunctionArgType()) { funcArgOpVector.push_back(this->PTR); } this->PTR->getFuncArgOps(funcArgOpVector); }

#define AST_GET_FUNC_OPS(PTR)  if( !(Teuchos::is_null(this->PTR)) ) { if (this->PTR->funcType()) { funcOpVector.push_back(this->PTR); } this->PTR->getFuncOps(funcOpVector); }

#define AST_GET_VOLT_OPS(PTR)  if( !(Teuchos::is_null(this->PTR)) ) { if (this->PTR->voltageType()) { voltOpVector.push_back(this->PTR); } this->PTR->getVoltageOps(voltOpVector); }

#define AST_GET_CURRENT_OPS(PTR)  if( !(Teuchos::is_null(this->PTR)) ) { if (this->PTR->currentType()) { currentOpVector.push_back(this->PTR); } this->PTR->getCurrentOps(currentOpVector); }

#define AST_GET_INTERNAL_DEV_VAR_OPS(PTR)  if( !(Teuchos::is_null(this->PTR)) ) { if (this->PTR->internalDeviceVarType()) { internalDevVarOpVector.push_back(this->PTR); } this->PTR->getInternalDevVarOps(internalDevVarOpVector); }

#define AST_GET_TIME_OPS(PTR)  if( !(Teuchos::is_null(this->PTR)) ) { if (this->PTR->timeSpecialType()) { timeOpVector.push_back(this->PTR); } this->PTR->getTimeOps(timeOpVector); }

// this one adds "this"
#define AST_GET_INTERESTING_OPS2(PTR) AST_GET_INTERESTING_OPS (this->PTR) 
#define AST_GET_STATE_OPS2(PTR) AST_GET_STATE_OPS (this->PTR) 


//-------------------------------------------------------------------------------
// this is to make the call to "getInterestingOps" have a single 
// function argument that never has to change.
template <typename ScalarT>
struct opVectorContainers
{
public:
  opVectorContainers(
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & param,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & func,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & volt,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & current,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & leadCurrent,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & bsrcCurrent,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & power,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & internalDevVar,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & dnoNoiseDevVar,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & dniNoiseDevVar,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & oNoise,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & iNoise,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & sdt,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & ddt,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & src,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & stp,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & comp,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & limit,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & phase,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & sparam,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & yparam,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & zparam,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & agauss,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & gauss,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & aunif,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & unif,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & rand,
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & twoArgLimit,
  bool timeDep,
  bool tempDep,
  bool vTDep,
  bool FreqDep,
  bool gminDep 
      ):
  paramOpVector(param),
    funcOpVector(func),
    voltOpVector(volt),
    currentOpVector(current),
    leadCurrentOpVector(leadCurrent),
    bsrcCurrentOpVector(bsrcCurrent),
    powerOpVector(power),
    internalDevVarOpVector(internalDevVar),
    dnoNoiseDevVarOpVector(dnoNoiseDevVar),
    dniNoiseDevVarOpVector(dniNoiseDevVar),
    oNoiseOpVector(oNoise),
    iNoiseOpVector(iNoise),
    sdtOpVector(sdt),
    ddtOpVector(ddt),
    srcOpVector(src),
    stpOpVector(stp),
    compOpVector(comp),
    limitOpVector(limit),
    phaseOpVector(phase),
    sparamOpVector(sparam),
    yparamOpVector(yparam),
    zparamOpVector(zparam),
    agaussOpVector(agauss),
    gaussOpVector(gauss),
    aunifOpVector(aunif),
    unifOpVector(unif),
    randOpVector(rand),
    twoArgLimitOpVector(twoArgLimit),
    isTimeDependent(timeDep),
    isTempDependent(tempDep),
    isVTDependent(vTDep),
    isFreqDependent(FreqDep),
    isGminDependent(gminDep)
  {};

  std::vector< Teuchos::RCP<astNode<ScalarT> > > & paramOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & funcOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & voltOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & currentOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & leadCurrentOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & bsrcCurrentOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & powerOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & internalDevVarOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & dnoNoiseDevVarOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & dniNoiseDevVarOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & oNoiseOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & iNoiseOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & sdtOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & ddtOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & srcOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & stpOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & compOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & limitOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & phaseOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & sparamOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & yparamOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & zparamOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & agaussOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & gaussOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & aunifOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & unifOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & randOpVector;
  std::vector< Teuchos::RCP<astNode<ScalarT> > > & twoArgLimitOpVector;

  bool isTimeDependent;
  bool isTempDependent;
  bool isVTDependent;
  bool isFreqDependent;
  bool isGminDependent;
};

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

    virtual ScalarT val() = 0;
    virtual ScalarT dx(int i) = 0;

    virtual std::vector<ScalarT> dx2(int numDerivs)
    {
      if (derivVec_.empty()) { derivVec_.resize(numDerivs,0.0); }
      else { std::fill(derivVec_.begin(),derivVec_.end(),0.0); }
      return derivVec_;
    };

    virtual void output(std::ostream & os, int indent=0) = 0;
    virtual void compactOutput(std::ostream & os) = 0;
    virtual void codeGen (std::ostream & os ) 
    {
      os << "// This node has not implemented a code gen function yet" <<std::endl;
    }

    virtual void setNode(Teuchos::RCP<astNode<ScalarT> > & tmpNode) {};
    virtual void unsetNode() {};

    virtual void setFuncArgs(const std::vector< Teuchos::RCP<astNode<ScalarT> > > & tmpArgVec ) {};

    virtual void setupBreakPoints() { return; }
    virtual bool getBreakPoints(std::vector<Xyce::Util::BreakPoint> & breakPointTimes) { return true;}
    virtual void setBreakPointTol(double tol){return;};
    virtual void setStartingTimeStep(double timeStep){return;};
    virtual void setFinalTime(double finalTime){return;};

    virtual void setDerivIndex(int i) {};
    virtual void unsetDerivIndex() {};

    virtual ScalarT getValue() { return 0.0; }
    virtual void setValue(ScalarT val) {}; // supports specialsOp, paramOp and globalParamLayerOp otherwise no-op
    virtual void unsetValue() {};          // supports specialsOp, paramOp and globalParamLayerOp otherwise no-op

    // base class no-ops.  Derived functions only in paramOp, base class version only called from ddx.
    virtual void setIsVar() {};
    virtual void unsetIsVar() {};
    virtual bool getIsVar() { return false; }

    virtual void setIsConstant() {};
    virtual void unsetIsConstant() {};
    virtual bool getIsConstant() { return false; }

    virtual void setIsAttached() {};
    virtual void unsetIsAttached() {};
    virtual bool getIsAttached() { return false; }

    virtual void setIsDotParam() {};
    virtual void unsetIsDotParam() {};
    virtual bool getIsDotParam() { return false; }

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

    virtual std::string getName () { return std::string(""); };
    //virtual std::vector<std::string> getNodeNames() { std::vector<std::string> tmp; return tmp; }

    virtual void getInterestingOps(opVectorContainers<ScalarT> & ovc)
    {
AST_GET_INTERESTING_OPS(leftAst_) AST_GET_INTERESTING_OPS(rightAst_) 
    }

    virtual void getStateOps(stateOpVectorContainers<ScalarT> & ovc)
    {
AST_GET_STATE_OPS(leftAst_) AST_GET_STATE_OPS(rightAst_) 
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
AST_GET_PARAM_OPS(leftAst_) AST_GET_PARAM_OPS(rightAst_) 
    }

    // func arg ops are of class paramOp, but have been identified as being
    // function arguments.  (thus being excluded from the getParamOps function, above)
    // This is needed by the ddx function
    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
AST_GET_FUNC_ARG_OPS(leftAst_) AST_GET_FUNC_ARG_OPS(rightAst_) 
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
AST_GET_FUNC_OPS(leftAst_) AST_GET_FUNC_OPS(rightAst_) 
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
AST_GET_VOLT_OPS(leftAst_) AST_GET_VOLT_OPS(rightAst_) 
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
AST_GET_CURRENT_OPS(leftAst_) AST_GET_CURRENT_OPS(rightAst_) 
    }
    virtual void getInternalDevVarOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & internalDevVarOpVector)
    {
AST_GET_INTERNAL_DEV_VAR_OPS (leftAst_) AST_GET_INTERNAL_DEV_VAR_OPS (rightAst_) 
    }

    virtual void getTimeOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & timeOpVector)
    {
AST_GET_TIME_OPS(leftAst_) AST_GET_TIME_OPS(rightAst_) 
    }

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

    std::vector<ScalarT> derivVec_;

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

    virtual std::vector<ScalarT> dx2(int numDerivs)
    {
      if (this->derivVec_.empty()) { this->derivVec_.resize(numDerivs,0.0); }
      else { std::fill(this->derivVec_.begin(),this->derivVec_.end(),0.0); }
      return this->derivVec_;
    };

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

    virtual std::vector<std::complex<double> > dx2(int numDerivs)
    {
      if (derivVec_.empty()) { derivVec_.resize(numDerivs,std::complex<double>(0.0,0.0)); }
      else { std::fill(derivVec_.begin(),derivVec_.end(),std::complex<double>(0.0,0.0)); }
      return derivVec_;
    }

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

    virtual bool numvalType() { return true; };

    std::vector<std::complex<double> > derivVec_;
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

  if (leftAst_->timeSpecialType()) { timeOpVec_.push_back(leftAst_); }
  if (rightAst_->timeSpecialType()) { timeOpVec_.push_back(rightAst_); }
  leftAst_->getTimeOps(timeOpVec_);
  rightAst_->getTimeOps(timeOpVec_);

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
    double delta_bpTime =  0.0;
    if (std::real(dfdt) != 0.0) { delta_bpTime =  -std::real(f)/std::real(dfdt); }
    double time = std::real(timeOpVec_[0]->val());
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

    // cleanup, restore
    for (int ii=0; ii< timeOpVec_.size(); ii++) { timeOpVec_[ii]->unsetDerivIndex(); }
    for (int ii=0; ii< timeOpVec_.size(); ii++) { timeOpVec_[ii]->setValue(time); }

    if (std::abs(std::real(f)) <= bpTol_) { bpTimes_.push_back( bpTime ); }// save this breakpoint if we converged
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

    virtual std::vector<ScalarT> dx2(int numDerivs)
    {
      Teuchos::RCP<astNode<ScalarT> > & lef = this->leftAst_;
      Teuchos::RCP<astNode<ScalarT> > & rig = this->rightAst_;

      if (this->derivVec_.empty()) { this->derivVec_.resize(numDerivs,0.0); }
      std::vector<ScalarT> & retVal = this->derivVec_;

      ScalarT leftVal=lef->val();
      ScalarT righVal=rig->val();

      if (lefDerivs_.empty()) { lefDerivs_.resize(numDerivs,0.0); }
      if (rigDerivs_.empty()) { rigDerivs_.resize(numDerivs,0.0); }

      lefDerivs_ = lef->dx2(numDerivs);
      rigDerivs_ = rig->dx2(numDerivs);

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
      return  retVal;
    }

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

  private:
    bool rightConst_;
    bool leftConst_;

    std::vector<ScalarT> lefDerivs_;
    std::vector<ScalarT> rigDerivs_;
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

    virtual std::vector<ScalarT> dx2(int numDerivs)
    {
      Teuchos::RCP<astNode<ScalarT> > & lef = this->leftAst_;
      Teuchos::RCP<astNode<ScalarT> > & rig = this->rightAst_;

      if (this->derivVec_.empty()) { this->derivVec_.resize(numDerivs,0.0); }
      std::vector<ScalarT> & retVal = this->derivVec_;

      ScalarT leftVal=lef->val();
      ScalarT righVal=rig->val();

      if (lefDerivs_.empty()) { lefDerivs_.resize(numDerivs,0.0); }
      if (rigDerivs_.empty()) { rigDerivs_.resize(numDerivs,0.0); }
      lefDerivs_ = lef->dx2(numDerivs);
      rigDerivs_ = rig->dx2(numDerivs);

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

      return  retVal;
    }

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

  private:
    bool rightConst_;
    bool leftConst_;

    std::vector<ScalarT> lefDerivs_;
    std::vector<ScalarT> rigDerivs_;
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

    virtual void getInterestingOps(opVectorContainers<ScalarT> & ovc)
    {
AST_GET_INTERESTING_OPS(paramNode_)
    }

    virtual void getStateOps(stateOpVectorContainers<ScalarT> & ovc)
    {
AST_GET_STATE_OPS(paramNode_)
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
AST_GET_PARAM_OPS(paramNode_)
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
AST_GET_FUNC_ARG_OPS(paramNode_)
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
AST_GET_FUNC_OPS(paramNode_)
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
AST_GET_VOLT_OPS(paramNode_)
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
AST_GET_CURRENT_OPS(paramNode_)
    }

    virtual void getTimeOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & timeOpVector)
    {
AST_GET_TIME_OPS(paramNode_)
    }

    virtual void processSuccessfulTimeStep () 
    {
      paramNode_->processSuccessfulTimeStep ();
    };

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

    virtual std::vector<ScalarT> dx2(int numDerivs)
    {
      if (this->derivVec_.empty()) { this->derivVec_.resize(numDerivs,0.0); }
      else { std::fill(this->derivVec_.begin(),this->derivVec_.end(),0.0); }
      std::vector<ScalarT> & retval = this->derivVec_;

      if (isVar_) { retval[derivIndex_] = 1.0; }
      else        { retval = paramNode_->dx2(numDerivs); }

      return retval;
    }

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

    virtual void getInterestingOps(opVectorContainers<ScalarT> & ovc)
    {
AST_GET_INTERESTING_OPS(paramNode_)
    }

    virtual void getStateOps(stateOpVectorContainers<ScalarT> & ovc)
    {
AST_GET_STATE_OPS(paramNode_)
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
AST_GET_PARAM_OPS(paramNode_)
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
AST_GET_FUNC_ARG_OPS(paramNode_)
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
AST_GET_FUNC_OPS(paramNode_)
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
AST_GET_VOLT_OPS(paramNode_)
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
AST_GET_CURRENT_OPS(paramNode_)
    }

    virtual void getTimeOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & timeOpVector)
    {
AST_GET_TIME_OPS(paramNode_)
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
    void setIsConstant() { isConstant_ = true; }
    void unsetIsConstant() { isConstant_ = false; }
    bool getIsConstant() { return isConstant_; }

    // this flag indicates if an external AST has been attached
    // to this class
    void setIsAttached() { isAttached_ = true; }
    void unsetIsAttached() { isAttached_ = false; }
    bool getIsAttached() { return isAttached_; }

    // the param type can be .param, .global_param or a subcircuit argument
    // The enum is defined as enum enumParamType {DOT_PARAM, DOT_GLOBAL_PARAM, SUBCKT_ARG_PARAM}
    void setParamType(enumParamType type) { paramType_ = type; }
    void unsetIsDotParam() { paramType_ = DOT_GLOBAL_PARAM; }
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

    virtual std::vector<ScalarT> dx2(int numDerivs)
    {
      if (this->derivVec_.empty()) { this->derivVec_.resize(numDerivs,0.0); }
      else { std::fill(this->derivVec_.begin(),this->derivVec_.end(),0.0); }
      std::vector<ScalarT> & retval = this->derivVec_;
      retval[derivIndex_] = 1.0;
      return retval;
    }

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

    std::string & getVoltageNode() { return voltageNode_; }
    ScalarT & getVoltageVal() { return voltageVal_; }

    virtual bool voltageType() { return true; };

    virtual std::string getName () { return voltageNode_; }

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

    virtual std::vector<ScalarT> dx2(int numDerivs)
    {
      if (this->derivVec_.empty()) { this->derivVec_.resize(numDerivs,0.0); }
      else { std::fill(this->derivVec_.begin(),this->derivVec_.end(),0.0); }
      std::vector<ScalarT> & retval = this->derivVec_;
      retval[derivIndex_] = 1.0;
      return retval;
    }

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
    std::string & getCurrentDevice() { return currentDevice_; }

    ScalarT & getCurrentVal () { return number_; }
    void setCurrentVal (ScalarT n) { number_ = n; }

    virtual bool currentType() { return true; };
    virtual bool bsrcCurrentType() { return bsrcFlag_; }

    virtual std::string getName () { return currentDevice_; }

    bool getBsrcFlag   () { return bsrcFlag_; }
    void setBsrcFlag  () { bsrcFlag_ = true; }
    void unsetBsrcFlag  () { bsrcFlag_ = false; }

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

    virtual void setValue(ScalarT val) { number_ = val; };

    std::vector<int> & getSparamArgs () { return sparamArgs_; }

    virtual bool sparamType() { return true; };

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

    virtual void setValue(ScalarT val) { number_ = val; };

    std::vector<int> & getYparamArgs () { return yparamArgs_; }

    virtual bool yparamType() { return true; };

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

    virtual void setValue(ScalarT val) { number_ = val; };

    std::vector<int> & getZparamArgs () { return zparamArgs_; }

    virtual bool zparamType() { return true; };

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
    std::string getLeadCurrentDevice() { return leadCurrentDevice_; }

    void setLeadCurrentDesignator(const std::string & desName) { leadCurrentDesignator_ = desName; }
    std::string getLeadCurrentDesignator() { return leadCurrentDesignator_; }

    ScalarT getLeadCurrentVar () { return number_; }
    void setLeadCurrentVar (ScalarT n) { number_ = n; }

    virtual bool leadCurrentType() { return true; };

    virtual std::string getName () { return leadCurrentDevice_; }

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
    std::string getPowerTag   () { return tag_; }
    std::string getPowerDevice() { return powerDevice_; }
    ScalarT getPowerVal () { return number_; }
    void setPowerVal (ScalarT n) { number_ = n; }

    virtual bool powerType() { return true; };

    virtual std::string getName () { return powerDevice_; }

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
    std::string getInternalVarDevice() { return internalDevVarDevice_; }
    ScalarT getInternalDeviceVar () { return number_; }
    void setInternalDeviceVar (ScalarT n) { number_ = n; }

    virtual bool internalDeviceVarType()  { return true; };

    virtual std::string getName () { return internalDevVarDevice_; }

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

    ScalarT getNoiseVar () { return number_; }
    void setNoiseVar (ScalarT n) { number_ = n; }

    virtual bool dnoNoiseVarType()  { return true; };

    //virtual std::string getName () { return noiseDevice_; }

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
    ScalarT getNoiseVar () { return number_; }
    void setNoiseVar (ScalarT n) { number_ = n; }

    virtual bool dniNoiseVarType()  { return true; };

    //virtual std::string getName () { return noiseDevice_; }

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
    ScalarT getNoiseVar () { return number_; }
    void setNoiseVar (ScalarT n) { number_ = n; }
    virtual bool oNoiseType()  { return true; };

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
    ScalarT getNoiseVar () { return number_; }
    void setNoiseVar (ScalarT n) { number_ = n; }
    virtual bool iNoiseType()  { return true; };

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
// *       = handled elsewhere in this AST tree as a binyarMultOp node
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
    virtual std::vector<ScalarT> dx2(int numDerivs)
    {
      if (this->derivVec_.empty()) { this->derivVec_.resize(numDerivs,0.0); }
      std::vector<ScalarT> & dfdx = this->derivVec_;

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
          dfdx = functionNode_->dx2(numDerivs);
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

            if (dpdx_.empty()) { dpdx_.resize(numDerivs,0.0); }
            dpdx_ = funcArgs_[ii]->dx2(numDerivs); // usually zero ...

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
      return dfdx;
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

    virtual void getInterestingOps(opVectorContainers<ScalarT> & ovc)
    {
      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->setNode( funcArgs_[ii] ); }
AST_GET_INTERESTING_OPS(functionNode_) 
      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->unsetNode(); } // restore
    }

    virtual void getStateOps(stateOpVectorContainers<ScalarT> & ovc)
    {
      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->setNode( funcArgs_[ii] ); }
AST_GET_STATE_OPS(functionNode_) 
      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->unsetNode(); } // restore
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->setNode( funcArgs_[ii] ); }
AST_GET_PARAM_OPS(functionNode_) 
      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->unsetNode(); } // restore
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->setNode( funcArgs_[ii] ); }
AST_GET_FUNC_ARG_OPS(functionNode_) 
      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->unsetNode(); } // restore
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->setNode( funcArgs_[ii] ); }
AST_GET_FUNC_OPS(functionNode_)
      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->unsetNode(); } // restore
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->setNode( funcArgs_[ii] ); }
AST_GET_VOLT_OPS(functionNode_)
      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->unsetNode(); } // restore
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->setNode( funcArgs_[ii] ); }
AST_GET_CURRENT_OPS(functionNode_) 
      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->unsetNode(); } // restore
    }

    virtual void getTimeOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & timeOpVector)
    {
      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->setNode( funcArgs_[ii] ); }
AST_GET_TIME_OPS(functionNode_) 
      if(dummyFuncArgs_.size() == funcArgs_.size())
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->unsetNode(); } // restore
    }

    bool getNodeResolved() { return nodeResolved_; }
    bool getArgsResolved() { return argsResolved_; }

    virtual void processSuccessfulTimeStep () 
    {
      functionNode_->processSuccessfulTimeStep ();
    };

    ddtStateData<ScalarT> & getDdtState() { return functionNode_->getDdtState(); }
    sdtStateData<ScalarT> & getSdtState() { return functionNode_->getSdtState(); }

    virtual unsigned long int getNodeId () { return functionNode_->getNodeId(); }

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

    std::vector<ScalarT> dpdx_;
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

    virtual std::vector<ScalarT> dx2(int numDerivs)
    {
      Teuchos::RCP<astNode<ScalarT> > & lef = this->leftAst_;
      Teuchos::RCP<astNode<ScalarT> > & rig = this->rightAst_;

      if (this->derivVec_.empty()) { this->derivVec_.resize(numDerivs,0.0); }
      std::vector<ScalarT> & retVal = this->derivVec_;

      ScalarT leftVal=lef->val();
      ScalarT righVal=rig->val();

      if (leftVal != 0.0) 
      {
        if (rightConst_ && !leftConst_) 
        {
          if (lefDerivs_.empty()) { lefDerivs_.resize(numDerivs,0.0); }
          lefDerivs_ = lef->dx2(numDerivs);
          if (std::real(leftVal) >= 0)
          {
            for (int ii=0;ii<numDerivs;ii++)
            {
              retVal[ii] = righVal*lefDerivs_[ii]/leftVal*std::pow(leftVal,righVal) ; 
            }
          }
          else if (std::real(leftVal) < 0)
          {
            for (int ii=0;ii<numDerivs;ii++)
            {
              retVal[ii] = righVal*(-lefDerivs_[ii])/(-leftVal)*std::pow((-leftVal),righVal) ; 
            }
          }
        }
        else if (!rightConst_ && leftConst_) 
        {
          if (rigDerivs_.empty()) { rigDerivs_.resize(numDerivs,0.0); }
          rigDerivs_ = rig->dx2(numDerivs);
          if (std::real(leftVal) >= 0)
          {
            for (int ii=0;ii<numDerivs;ii++)
            {
              retVal[ii] = std::log(leftVal)*std::pow(leftVal,righVal)*rigDerivs_[ii]; 
            }
          }
          else if (std::real(leftVal) < 0)
          {
            for (int ii=0;ii<numDerivs;ii++)
            {
              retVal[ii] = -std::log(-leftVal)*std::pow(-leftVal,righVal)*rigDerivs_[ii]; 
            }
          }
        }
        else 
        {
          if (lefDerivs_.empty()) { lefDerivs_.resize(numDerivs,0.0); }
          if (rigDerivs_.empty()) { rigDerivs_.resize(numDerivs,0.0); }
          lefDerivs_ = lef->dx2(numDerivs);
          rigDerivs_ = rig->dx2(numDerivs);
          if (std::real(leftVal) >= 0)
          {
            for (int ii=0;ii<numDerivs;ii++)
            {
              retVal[ii] = (rigDerivs_[ii]*std::log(leftVal)+righVal*lefDerivs_[ii]/leftVal)*std::pow(leftVal,righVal);
            }
          }
          else if (std::real(leftVal) < 0)
          {
            for (int ii=0;ii<numDerivs;ii++)
            {
              retVal[ii] = (-rigDerivs_[ii]*std::log(-leftVal)+righVal*(-lefDerivs_[ii])/(-leftVal))*std::pow((-leftVal),righVal);
            }
          }
        }
      }
      return  retVal;
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

  private:
    bool rightConst_;
    bool leftConst_;
    std::vector<ScalarT> lefDerivs_;
    std::vector<ScalarT> rigDerivs_;
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

};

//-------------------------------------------------------------------------------
// fmod operator. returns the remainder of the division as a real number
//
template <typename ScalarT>
class fmodOp : public astNode<ScalarT>
{
  public:
    fmodOp (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):
      astNode<ScalarT>(left,right), rightConst_(true),leftConst_(false)
    {
      rightConst_ = this->rightAst_->numvalType();
      leftConst_ = this->leftAst_->numvalType();
    }

    virtual ScalarT val()
    {
      return std::fmod ( std::real(this->leftAst_->val()) , std::real(this->rightAst_->val()));
    }

    virtual ScalarT dx (int i)
    {
      Teuchos::RCP<astNode<ScalarT> > & lef = this->leftAst_;
      Teuchos::RCP<astNode<ScalarT> > & rig = this->rightAst_;
      ScalarT retVal = 0.0;
      return  retVal;
    }

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

  private:
    bool rightConst_;
    bool leftConst_;
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

    virtual std::vector<ScalarT> dx2(int numDerivs)
    {
      Teuchos::RCP<astNode<ScalarT> > & x = (this->leftAst_);
      Teuchos::RCP<astNode<ScalarT> > & y = (this->rightAst_);
      Teuchos::RCP<astNode<ScalarT> > & z = (zAst_);

      if (yDerivs_.empty()) { yDerivs_.resize(numDerivs,0.0); }
      if (zDerivs_.empty()) { zDerivs_.resize(numDerivs,0.0); }
      if (this->derivVec_.empty()) { this->derivVec_.resize(numDerivs,0.0); }

      ScalarT xVal = x->val();
      yDerivs_ = y->dx2(numDerivs);
      zDerivs_ = z->dx2(numDerivs);

      // not "fixing" x->val() b/c it is the result of a conditional, which is 1 or 0.
      // The correct place to fix this is in the comparison operators.  Fix later.
      for (int ii=0;ii<numDerivs;ii++)
      {
        ScalarT dyFixed = yDerivs_[ii];  Xyce::Util::fixNan(dyFixed);  Xyce::Util::fixInf(dyFixed);
        ScalarT dzFixed = zDerivs_[ii];  Xyce::Util::fixNan(dzFixed);  Xyce::Util::fixInf(dzFixed);
        this->derivVec_[ii] = ((std::real(xVal))?(dyFixed):(dzFixed));
      }
      return this->derivVec_;
    };

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

    virtual void getInterestingOps(opVectorContainers<ScalarT> & ovc)
    {
AST_GET_INTERESTING_OPS2(leftAst_) AST_GET_INTERESTING_OPS2(rightAst_) AST_GET_INTERESTING_OPS(zAst_)
    }

    virtual void getStateOps(stateOpVectorContainers<ScalarT> & ovc)
    {
AST_GET_STATE_OPS2(leftAst_) AST_GET_STATE_OPS2(rightAst_) AST_GET_STATE_OPS(zAst_)
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
AST_GET_PARAM_OPS(leftAst_) AST_GET_PARAM_OPS(rightAst_) AST_GET_PARAM_OPS(zAst_) 
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
AST_GET_FUNC_ARG_OPS(leftAst_) AST_GET_FUNC_ARG_OPS(rightAst_) AST_GET_FUNC_ARG_OPS(zAst_) 
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
AST_GET_FUNC_OPS(leftAst_) AST_GET_FUNC_OPS(rightAst_) AST_GET_FUNC_OPS(zAst_) 
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
AST_GET_VOLT_OPS(leftAst_) AST_GET_VOLT_OPS(rightAst_) AST_GET_VOLT_OPS(zAst_) 
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
AST_GET_CURRENT_OPS(leftAst_) AST_GET_CURRENT_OPS(rightAst_) AST_GET_CURRENT_OPS(zAst_) 
    }

    virtual void getTimeOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & timeOpVector)
    {
AST_GET_TIME_OPS(leftAst_) AST_GET_TIME_OPS(rightAst_) AST_GET_TIME_OPS(zAst_)
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > zAst_;
    std::vector<ScalarT> yDerivs_;
    std::vector<ScalarT> zDerivs_;
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

    virtual std::vector<ScalarT> dx2 (int numDerivs)
    {
      Teuchos::RCP<astNode<ScalarT> > & x = (this->leftAst_);
      Teuchos::RCP<astNode<ScalarT> > & y = (this->rightAst_);
      Teuchos::RCP<astNode<ScalarT> > & z = (zAst_);

      if (xDerivs_.empty()) { xDerivs_.resize(numDerivs,0.0); }
      if (this->derivVec_.empty()) { this->derivVec_.resize(numDerivs,0.0); }

      ScalarT xFixed = x->val();  Xyce::Util::fixNan(xFixed);  Xyce::Util::fixInf(xFixed);
      ScalarT yFixed = y->val();  Xyce::Util::fixNan(yFixed);  Xyce::Util::fixInf(yFixed);
      ScalarT zFixed = z->val();  Xyce::Util::fixNan(zFixed);  Xyce::Util::fixInf(zFixed);

      xDerivs_ = x->dx2(numDerivs);

      std::vector<ScalarT> & retDerivs = this->derivVec_;
      for (int ii=0;ii<numDerivs;ii++)
      {
        ScalarT dxFixed = xDerivs_[ii]; Xyce::Util::fixNan(dxFixed); Xyce::Util::fixInf(dxFixed);
        retDerivs[ii] = ((std::real(xFixed)<std::real(yFixed))?0.0:((std::real(xFixed)>std::real(zFixed))?0.0:(dxFixed)));
      }

      return retDerivs;
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

    virtual void getInterestingOps(opVectorContainers<ScalarT> & ovc)
    {
AST_GET_INTERESTING_OPS2(leftAst_) AST_GET_INTERESTING_OPS2(rightAst_) AST_GET_INTERESTING_OPS(zAst_)
    }

    virtual void getStateOps(stateOpVectorContainers<ScalarT> & ovc)
    {
AST_GET_STATE_OPS2(leftAst_) AST_GET_STATE_OPS2(rightAst_) AST_GET_STATE_OPS(zAst_)
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
AST_GET_PARAM_OPS(leftAst_) AST_GET_PARAM_OPS(rightAst_) AST_GET_PARAM_OPS(zAst_) 
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
AST_GET_FUNC_ARG_OPS(leftAst_) AST_GET_FUNC_ARG_OPS(rightAst_) AST_GET_FUNC_ARG_OPS(zAst_) 
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
AST_GET_FUNC_OPS(leftAst_) AST_GET_FUNC_OPS(rightAst_) AST_GET_FUNC_OPS(zAst_)
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
AST_GET_VOLT_OPS(leftAst_) AST_GET_VOLT_OPS(rightAst_) AST_GET_VOLT_OPS(zAst_) 
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
AST_GET_CURRENT_OPS(leftAst_) AST_GET_CURRENT_OPS(rightAst_) AST_GET_CURRENT_OPS(zAst_) 
    }

    virtual void getTimeOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & timeOpVector)
    {
AST_GET_TIME_OPS(leftAst_) AST_GET_TIME_OPS(rightAst_) AST_GET_TIME_OPS(zAst_) 
    }

    virtual bool limitType() { return true; }

  private:
    Teuchos::RCP<astNode<ScalarT> > zAst_;
    std::vector<Teuchos::RCP<astNode<ScalarT> > > timeOpVec_;
    double bpTol_;
    std::vector<Xyce::Util::BreakPoint> bpTimes_;

    std::vector<ScalarT> xDerivs_;
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

    virtual bool stpType() { return true; }

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

    virtual ScalarT dx (int i)
    {
      return ((std::real(this->leftAst_->val()))>0)?1.0:0.0;
    }

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
};

inline bool isLeftCurlyBrace(char c) { return (c=='{'); }

inline bool isLeftParen(char c) { return (c=='('); }

inline bool isQuoteSymbol(char c) { return (c=='"'); }

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
      allNumVal_(true), 
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
          allNumVal_=true; ta_.resize(size/2); ya_.resize(size/2); dya_.resize(size/2,0.0);
          for (int ii=0,jj=0;ii<size;ii+=2,jj++)
          {
            ta_[jj] = (tableArgs_)[ii]->val();
            ya_[jj] = (tableArgs_)[ii+1]->val();
            if (!( (tableArgs_)[ii]->numvalType() && (tableArgs_)[ii+1]->numvalType() ) ) { allNumVal_ = false; }
          }
          yInterpolator_->init(ta_,ya_); // for linear, this isn't necessary, but for others it is

          if (ya_.size() > 2 && ( keyword_==std::string("TABLE") || keyword_==std::string("FASTTABLE") ) )
          {
            // create derivative table
            // this code mimics the old expression library.   It uses finite differencing
            // to set up a new table of derivatives.  The new table is based on the midpoints of
            // the original table, so it has one extra entry.
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
          }
        }
      };

    //-------------------------------------------------------------------------------
    // special constructor for values read in from a file, in which the file IO is 
    // handled directly in this constructor.
    tableOp (const std::string & kw, Teuchos::RCP<astNode<ScalarT> > &input, const std::string & filename):
      astNode<ScalarT>(), allNumVal_(true), input_(input),
      useBreakPoints_(true),
      keyword_(kw)
      {
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
            while ( dataIn >> time )
            {
              if ( dataIn >> value )
              {
                ta_.push_back(time);
                ya_.push_back(value);
              }
              else
              {
                std::vector<std::string> errStr(1,std::string("Reached end of file in " + filename + " while expecting another value"));
                yyerror(errStr);
              }
            }
          }
          dataIn.close();
        }

        yInterpolator_->init(ta_,ya_); // for linear, this isn't necessary, but for others it is

        if (ya_.size() > 2 && ( keyword_==std::string("TABLE") || keyword_==std::string("FASTTABLE") ) )
        {
          // create derivative table
          // this code mimics the old expression library.   It uses finite differencing
          // to set up a new table of derivatives.  The new table is based on the midpoints of
          // the original table, so it has one extra entry.
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
        }
      };

    //-------------------------------------------------------------------------------
    // special constructor for values read in from a file, that are now stored in std::vector objects
    // ERK.  Currently, Xyce doesn't use this function, but it should, as it is more reliable than 
    // the above numvalType test in the first constructor.
    tableOp (const std::string & kw, Teuchos::RCP<astNode<ScalarT> > & input, const std::vector<ScalarT> & xvals, const std::vector<ScalarT> & yvals):
      astNode<ScalarT>(), allNumVal_(true), input_(input),
      useBreakPoints_(true),
      keyword_(kw)
      {
        allocateInterpolators();

        allNumVal_=true; int size = xvals.size(); int size2=yvals.size();
        if (size != size2)
        {
          std::vector<std::string> errStr(1,std::string("AST node (table) needs x and y vectors to be the same size.")); yyerror(errStr);
        }
        ta_.resize(size); ya_.resize(size); dya_.resize(size,0.0);

        for (int ii=0;ii<size;ii++)
        {
          ta_[ii] = xvals[ii];
          ya_[ii] = yvals[ii];
        }

        yInterpolator_->init(ta_,ya_); // for linear, this isn't necessary, but for others it is

        if (ya_.size() > 2 && ( keyword_==std::string("TABLE") || keyword_==std::string("FASTTABLE") ) )
        {
          // create derivative table
          // this code mimics the old expression library.   It uses finite differencing
          // to set up a new table of derivatives.  The new table is based on the midpoints of
          // the original table, so it has one extra entry.
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
        }
      };


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

      if (!allNumVal_)  // if not all pure numbers, then initialize the arrays again
      {
        if(!(tableArgs_.empty()))
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

    ScalarT dx_linear(int i)
    {
      ScalarT dydx = 0.0;

      ScalarT dinput_dx = std::real(this->input_->dx(i));

      if (std::real(dinput_dx) != 0.0)
      {
        // derivative w.r.t. input
        //
        // this code mimics the old expression library.   It uses finite differencing
        // to set up a new table of derivatives.  The new table is based on the midpoints of
        // the original table, so it has one extra entry.
        //
        // I initially tried to use the evalDeriv function in the yInterpolator object.
        // That method doen't use midpoints, it just differentiates the the linear
        // interpolation device.  That approach failed at least one regression test.
        //
        if (!allNumVal_)  // if not all pure numbers, then initialize the arrays again
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

            int ya_size = ya_.size();
            if (ya_size > 2)
            {
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
        if (!allNumVal_)  // if not all pure numbers, then initialize the arrays again.  
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
        if (!allNumVal_)  // if not all pure numbers, then initialize the arrays again
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

    virtual void setupBreakPoints() {};
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
             while( std::real(ta_[tmp]) < std::real(time) && tmp <= size ) { tmp++; }
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

    virtual void getInterestingOps(opVectorContainers<ScalarT> & ovc)
    {

AST_GET_INTERESTING_OPS(input_)

      if (!allNumVal_)
      {
        if (!(tableArgs_.empty()))
        {
          int size=tableArgs_.size();
          for(int ii=0;ii<size;ii++)
          {
AST_GET_INTERESTING_OPS(tableArgs_[ii])
          }
        }
      }
    }

    virtual void getStateOps(stateOpVectorContainers<ScalarT> & ovc)
    {

AST_GET_STATE_OPS(input_)

      if (!allNumVal_)
      {
        if (!(tableArgs_.empty()))
        {
          int size=tableArgs_.size();
          for(int ii=0;ii<size;ii++)
          {
AST_GET_STATE_OPS(tableArgs_[ii])
          }
        }
      }
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
AST_GET_PARAM_OPS(input_) 

      if (!allNumVal_)
      {
        if (!(tableArgs_.empty()))
        {
          int size=tableArgs_.size();
          for(int ii=0;ii<size;ii++)
          {
AST_GET_PARAM_OPS(tableArgs_[ii]) 
          }
        }
      }
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
AST_GET_FUNC_ARG_OPS(input_) 

      if (!allNumVal_)
      {
        if (!(tableArgs_.empty()))
        {
          int size=tableArgs_.size();
          for(int ii=0;ii<size;ii++)
          {
AST_GET_FUNC_ARG_OPS(tableArgs_[ii]) 
          }
        }
      }
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
AST_GET_FUNC_OPS(input_) 

      if (!allNumVal_)
      {
        if (!(tableArgs_.empty()))
        {
          int size=tableArgs_.size();
          for(int ii=0;ii<size;ii++)
          {
AST_GET_FUNC_OPS(tableArgs_[ii]) 
          }
        }
      }
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
AST_GET_VOLT_OPS(input_) 

      if (!allNumVal_)
      {
        if (!(tableArgs_.empty()))
        {
          int size=tableArgs_.size();
          for(int ii=0;ii<size;ii++)
          {
AST_GET_VOLT_OPS(tableArgs_[ii] ) 
          }
        }
      }
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
AST_GET_CURRENT_OPS(input_) 

      if (!allNumVal_)
      {
        if (!(tableArgs_.empty()))
        {
          int size=tableArgs_.size();
          for(int ii=0;ii<size;ii++)
          {
AST_GET_CURRENT_OPS(tableArgs_[ii]) 
          }
        }
      }
    }

    virtual void getTimeOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & timeOpVector)
    {
AST_GET_TIME_OPS(input_) 

      if (!allNumVal_)
      {
        if (!(tableArgs_.empty()))
        {
          int size=tableArgs_.size();
          for(int ii=0;ii<size;ii++)
          {
AST_GET_TIME_OPS(tableArgs_[ii]) 
          }
        }
      }
    }

  private:
    std::vector<Teuchos::RCP<astNode<ScalarT> > > tableArgs_;
    bool allNumVal_;
    std::vector<ScalarT> ta_; // using ta for name instead of xa so as not to confuse meaning of dx function
    std::vector<ScalarT> ya_;
    std::vector<ScalarT> ta2_; // using ta for name instead of xa so as not to confuse meaning of dx function
    std::vector<ScalarT> dya_;

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

    virtual void getInterestingOps(opVectorContainers<ScalarT> & ovc)
    {

AST_GET_INTERESTING_OPS(time_)

      if (!allNumVal_)
      {
        int size=tableArgs_.size();
        for(int ii=0;ii<size;ii++)
        {
AST_GET_INTERESTING_OPS(tableArgs_[ii])
        }
      }
    }

    virtual void getStateOps(stateOpVectorContainers<ScalarT> & ovc)
    {

AST_GET_STATE_OPS(time_)

      if (!allNumVal_)
      {
        int size=tableArgs_.size();
        for(int ii=0;ii<size;ii++)
        {
AST_GET_STATE_OPS(tableArgs_[ii])
        }
      }
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
AST_GET_PARAM_OPS(time_) 

      if (!allNumVal_)
      {
        int size=tableArgs_.size();
        for(int ii=0;ii<size;ii++)
        {
AST_GET_PARAM_OPS(tableArgs_[ii]) 
        }
      }
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
AST_GET_FUNC_ARG_OPS(time_) 

      if (!allNumVal_)
      {
        int size=tableArgs_.size();
        for(int ii=0;ii<size;ii++)
        {
AST_GET_FUNC_ARG_OPS(tableArgs_[ii]) 
        }
      }
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
AST_GET_FUNC_OPS(time_) 

      if (!allNumVal_)
      {
        int size=tableArgs_.size();
        for(int ii=0;ii<size;ii++)
        {
AST_GET_FUNC_OPS(tableArgs_[ii]) 
        }
      }
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
AST_GET_VOLT_OPS(time_) 

      if (!allNumVal_)
      {
        int size=tableArgs_.size();
        for(int ii=0;ii<size;ii++)
        {
AST_GET_VOLT_OPS(tableArgs_[ii] ) 
        }
      }
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
AST_GET_CURRENT_OPS(time_) 

      if (!allNumVal_)
      {
        int size=tableArgs_.size();
        for(int ii=0;ii<size;ii++)
        {
AST_GET_CURRENT_OPS(tableArgs_[ii]) 
        }
      }
    }

    virtual void getTimeOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & timeOpVector)
    {
AST_GET_TIME_OPS(time_) 

      if (!allNumVal_)
      {
        int size=tableArgs_.size();
        for(int ii=0;ii<size;ii++)
        {
AST_GET_TIME_OPS(tableArgs_[ii]) 
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

    virtual std::vector<ScalarT> dx2 (int numDerivs)
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

      if (this->derivVec_.empty()) { this->derivVec_.resize(numDerivs,0.0); }
      std::vector<ScalarT> & dIdx = this->derivVec_;
      dIdx = this->leftAst_->dx2(numDerivs);
      ScalarT dIdVal2 = 0.5*deltaT;
      for(int ii=0;ii<numDerivs;ii++) { dIdx[ii] *= dIdVal2; }
      return dIdx;
    }

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

    virtual std::vector<ScalarT> dx2 (int numDerivs)
    {
      if (this->derivVec_.empty()) { this->derivVec_.resize(numDerivs,0.0); }
      std::vector<ScalarT> & ddt_dx = this->derivVec_;
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

          ddt_dx = this->leftAst_->dx2(numDerivs);
          for(int ii=0;ii<numDerivs;ii++)
          {
            ScalarT ddt_dVal2 = 1.0/deltaT;
            // for now, hardwire to backward Euler
            ddt_dx[ii] *= ddt_dVal2;
          }
        }
      }
      return ddt_dx;
    };

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
        if ( this->leftAst_->paramType() || this->leftAst_->getFunctionArgType() )
        { paramOpVector.push_back(this->leftAst_); }

        this->leftAst_->getParamOps(paramOpVector);
        this->leftAst_->getFuncArgOps(paramOpVector);

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
              if (tmp==tmp2) { foundX_ = true; astNodeX_ = paramOpVector[ii]; }
            }
          }
        }
      }
      else if (this->rightAst_->voltageType())
      {
        std::vector<Teuchos::RCP<astNode<ScalarT> > > voltOpVector;
        if (this->leftAst_->voltageType()) { voltOpVector.push_back(this->leftAst_); }

        this->leftAst_->getVoltageOps(voltOpVector);

        std::string tmp = this->rightAst_->getName();
        if (!(tmp.empty()))
        {
          Xyce::Util::toUpper(tmp);
          for (int ii=0;ii<voltOpVector.size();ii++)
          {
            std::string tmp2 = this->rightAst_->getName();
            if (!(tmp2.empty()))
            {
              Xyce::Util::toUpper(tmp2);
              if (tmp==tmp2) { foundX_ = true; astNodeX_ = voltOpVector[ii]; }
            }
          }
        }
      }
      else if (this->rightAst_->currentType())
      {
        std::vector<Teuchos::RCP<astNode<ScalarT> > > currentOpVector;

        if (this->leftAst_->currentType()) { currentOpVector.push_back(this->leftAst_); }

        this->leftAst_->getCurrentOps(currentOpVector);

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
              if (tmp==tmp2) { foundX_ = true; astNodeX_ = currentOpVector[ii]; }
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

  private:
    std::string type_;
    ScalarT value_;
    int derivIndex_;
};

//-------------------------------------------------------------------------------
// constants
template <typename ScalarT>
class piConstOp : public astNode<ScalarT>
{
  public:
    piConstOp (): astNode<ScalarT>() {};

    virtual ScalarT val() { return ScalarT(M_PI); };
    virtual ScalarT dx (int i) { return ScalarT(0.0); };

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "pi const operator.  val = " << ScalarT(M_PI) << " id = " << this->id_ << std::endl;
    }

    virtual void compactOutput(std::ostream & os) { output(os,0); }

    virtual void codeGen (std::ostream & os ) { os << ScalarT(M_PI); }

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

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "CtoK const operator.  val = " << ScalarT(CONSTCtoK) << " id = " << this->id_ << std::endl;
    }

    virtual void compactOutput(std::ostream & os) { output(os,0); }

    virtual void codeGen (std::ostream & os ) { os << ScalarT(CONSTCtoK); }

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

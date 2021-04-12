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

//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

#ifndef newExpression_H
#define newExpression_H

// std includes
#include<string>
#include<complex>
#include<utility>
#include<vector>
#include<list>

// trilinos includes
#include <Teuchos_RCP.hpp>

// Xyce includes.
#include <N_UTL_BreakPoint.h>
#include <expressionParamTypes.h>
#include <N_UTL_HspiceBools.h>
#include <N_UTL_Interface_Enum_Types.h>

// new code includes:
#include <ast.h>
#include <ExpressionType.h>
#include <expressionGroup.h>

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Class         : newExpression
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 10/28/2019
//-----------------------------------------------------------------------------
class newExpression
{
public:
  newExpression () :
    parsed_(false),
    derivsSetup_(false),
    astArraysSetup_(false),
    bpTol_(0.0),
    time_(0.0),
    timeStep_(0.0),
    timeStepAlpha_(0.0),
    timeStepPrefac_(0.0),
    stepNumber_(0),
    numDerivs_(0),
    traditionalParse_(true),
    externalDependencies_(false),
    isTimeDependent_(false),
    isTempDependent_(false),
    isVTDependent_(false),
    isFreqDependent_(false),
    isGminDependent_(false),

    isShallowTimeDependent_(false),
    isShallowTempDependent_(false),
    isShallowVTDependent_(false),
    isShallowFreqDependent_(false),
    isShallowGminDependent_(false),

    isVariableDependent_(false),
    isVoltageNodeDependent_(false),
    isDeviceCurrentDependent_(false),
    isLeadCurrentDependent_(false),
    isLeadCurrentDependentExcludeBsrc_(false),

    overrideGroupTemperature_(false),
    overrideTemp_(27.0),
    isConstant_(false),
    evaluateFunctionCalledBefore_(false),
    evaluateCalledBefore_(false),
    getTheSeedCalledBefore_(false),
    savedResult_(0.0),
    phaseOutputUsesRadians_(false),
    opVectors_(paramOpVec_,funcOpVec_, voltOpVec_, currentOpVec_, leadCurrentOpVec_, bsrcCurrentOpVec_, powerOpVec_, internalDevVarOpVec_, dnoNoiseDevVarOpVec_, dniNoiseDevVarOpVec_, oNoiseOpVec_, iNoiseOpVec_, sdtOpVec_, ddtOpVec_, srcAstNodeVec_, stpAstNodeVec_, compAstNodeVec_, limitAstNodeVec_, phaseOpVec_, sparamOpVec_, yparamOpVec_, zparamOpVec_, agaussOpVec_, gaussOpVec_, aunifOpVec_, unifOpVec_, randOpVec_, twoArgLimitOpVec_, isTimeDependent_, isTempDependent_, isVTDependent_, isFreqDependent_, isGminDependent_)
  {};

  // primary constructor
  newExpression ( std::string const & exp, Teuchos::RCP<baseExpressionGroup> & group ) :
    group_(group),
    expressionString_(exp),
    parsed_(false),
    derivsSetup_(false),
    astArraysSetup_(false),
    bpTol_(0.0),
    time_(0.0),
    timeStep_(0.0),
    timeStepAlpha_(0.0),
    timeStepPrefac_(0.0),
    stepNumber_(0),
    numDerivs_(0),
    traditionalParse_(true),
    externalDependencies_(false),
    isTimeDependent_(false),
    isTempDependent_(false),
    isVTDependent_(false),
    isFreqDependent_(false),
    isGminDependent_(false),

    isShallowTimeDependent_(false),
    isShallowTempDependent_(false),
    isShallowVTDependent_(false),
    isShallowFreqDependent_(false),
    isShallowGminDependent_(false),

    isVariableDependent_(false),
    isVoltageNodeDependent_(false),
    isDeviceCurrentDependent_(false),
    isLeadCurrentDependent_(false),
    isLeadCurrentDependentExcludeBsrc_(false),

    overrideGroupTemperature_(false),
    overrideTemp_(27.0),
    isConstant_(false),
    evaluateFunctionCalledBefore_(false),
    evaluateCalledBefore_(false),
    getTheSeedCalledBefore_(false),
    savedResult_(0.0),
    phaseOutputUsesRadians_(false),
    opVectors_(paramOpVec_,funcOpVec_, voltOpVec_, currentOpVec_, leadCurrentOpVec_, bsrcCurrentOpVec_, powerOpVec_, internalDevVarOpVec_, dnoNoiseDevVarOpVec_, dniNoiseDevVarOpVec_, oNoiseOpVec_, iNoiseOpVec_, sdtOpVec_, ddtOpVec_, srcAstNodeVec_, stpAstNodeVec_, compAstNodeVec_, limitAstNodeVec_, phaseOpVec_, sparamOpVec_, yparamOpVec_, zparamOpVec_, agaussOpVec_, gaussOpVec_, aunifOpVec_, unifOpVec_, randOpVec_, twoArgLimitOpVec_, isTimeDependent_, isTempDependent_, isVTDependent_, isFreqDependent_, isGminDependent_)
  {
    dtNodePtr_   = Teuchos::rcp(new specialsOp<usedType> (std::string("DT")));
    timeNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("TIME")));
    tempNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("TEMP")));
    vtNodePtr_   = Teuchos::rcp(new specialsOp<usedType> (std::string("VT")));
    freqNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("FREQ")));
    gminNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("GMIN")));
    piNodePtr_   = Teuchos::rcp(new piConstOp<usedType>  ());
    CtoKNodePtr_   = Teuchos::rcp(new CtoKConstOp<usedType>  ());
  };


  // special "big table" constructor
  newExpression (const std::vector<usedType> & xvals, const std::vector<usedType> & yvals,
      Teuchos::RCP<baseExpressionGroup> & group ) :
    group_(group),
    expressionString_("TIME"),
    parsed_(false),
    derivsSetup_(false),
    astArraysSetup_(false),
    bpTol_(0.0),
    time_(0.0),
    timeStep_(0.0),
    timeStepAlpha_(0.0),
    timeStepPrefac_(0.0),
    stepNumber_(0),
    numDerivs_(0),
    traditionalParse_(false),
    externalDependencies_(false),
    isTimeDependent_(false),
    isTempDependent_(false),
    isVTDependent_(false),
    isFreqDependent_(false),
    isGminDependent_(false),

    isShallowTimeDependent_(false),
    isShallowTempDependent_(false),
    isShallowVTDependent_(false),
    isShallowFreqDependent_(false),
    isShallowGminDependent_(false),

    isVariableDependent_(false),
    isVoltageNodeDependent_(false),
    isDeviceCurrentDependent_(false),
    isLeadCurrentDependent_(false),
    isLeadCurrentDependentExcludeBsrc_(false),

    overrideGroupTemperature_(false),
    overrideTemp_(27.0),
    isConstant_(false),
    evaluateFunctionCalledBefore_(false),
    evaluateCalledBefore_(false),
    getTheSeedCalledBefore_(false),
    savedResult_(0.0),
    phaseOutputUsesRadians_(false),
    opVectors_(paramOpVec_,funcOpVec_, voltOpVec_, currentOpVec_, leadCurrentOpVec_, bsrcCurrentOpVec_, powerOpVec_, internalDevVarOpVec_, dnoNoiseDevVarOpVec_, dniNoiseDevVarOpVec_, oNoiseOpVec_, iNoiseOpVec_, sdtOpVec_, ddtOpVec_, srcAstNodeVec_, stpAstNodeVec_, compAstNodeVec_, limitAstNodeVec_, phaseOpVec_, sparamOpVec_, yparamOpVec_, zparamOpVec_, agaussOpVec_, gaussOpVec_, aunifOpVec_, unifOpVec_, randOpVec_, twoArgLimitOpVec_, isTimeDependent_, isTempDependent_, isVTDependent_, isFreqDependent_, isGminDependent_)
  {
    dtNodePtr_   = Teuchos::rcp(new specialsOp<usedType> (std::string("DT")));
    timeNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("TIME")));
    tempNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("TEMP")));
    vtNodePtr_   = Teuchos::rcp(new specialsOp<usedType> (std::string("VT")));
    freqNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("FREQ")));
    gminNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("GMIN")));
    piNodePtr_   = Teuchos::rcp(new piConstOp<usedType>  ());
    CtoKNodePtr_   = Teuchos::rcp(new CtoKConstOp<usedType>  ());

    std::string keyword = std::string("TABLE");
    Teuchos::RCP<astNode<usedType> > time_base = timeNodePtr_;
    Teuchos::RCP<tableOp<usedType> > tableNodePtr_ = Teuchos::RCP<tableOp<usedType> >(new tableOp<usedType> (keyword, time_base, xvals, yvals));
    astNodePtr_ = Teuchos::RCP<astNode<usedType> >(tableNodePtr_);
  };

  // another "big table" constructor
  newExpression (Teuchos::RCP<astNode<usedType> > & left,
      const std::vector<usedType> & xvals, const std::vector<usedType> & yvals,
      Teuchos::RCP<baseExpressionGroup> & group ) :
    group_(group),
    expressionString_("TIME"),
    parsed_(false),
    derivsSetup_(false),
    astArraysSetup_(false),
    bpTol_(0.0),
    time_(0.0),
    timeStep_(0.0),
    timeStepAlpha_(0.0),
    timeStepPrefac_(0.0),
    stepNumber_(0),
    numDerivs_(0),
    traditionalParse_(false),
    externalDependencies_(false),
    isTimeDependent_(false),
    isTempDependent_(false),
    isVTDependent_(false),
    isFreqDependent_(false),
    isGminDependent_(false),

    isShallowTimeDependent_(false),
    isShallowTempDependent_(false),
    isShallowVTDependent_(false),
    isShallowFreqDependent_(false),
    isShallowGminDependent_(false),

    isVariableDependent_(false),
    isVoltageNodeDependent_(false),
    isDeviceCurrentDependent_(false),
    isLeadCurrentDependent_(false),
    isLeadCurrentDependentExcludeBsrc_(false),

    overrideGroupTemperature_(false),
    overrideTemp_(27.0),
    isConstant_(false),
    evaluateFunctionCalledBefore_(false),
    evaluateCalledBefore_(false),
    getTheSeedCalledBefore_(false),
    savedResult_(0.0),
    phaseOutputUsesRadians_(false),
    opVectors_(paramOpVec_,funcOpVec_, voltOpVec_, currentOpVec_, leadCurrentOpVec_, bsrcCurrentOpVec_, powerOpVec_, internalDevVarOpVec_, dnoNoiseDevVarOpVec_, dniNoiseDevVarOpVec_, oNoiseOpVec_, iNoiseOpVec_, sdtOpVec_, ddtOpVec_, srcAstNodeVec_, stpAstNodeVec_, compAstNodeVec_, limitAstNodeVec_, phaseOpVec_, sparamOpVec_, yparamOpVec_, zparamOpVec_, agaussOpVec_, gaussOpVec_, aunifOpVec_, unifOpVec_, randOpVec_, twoArgLimitOpVec_, isTimeDependent_, isTempDependent_, isVTDependent_, isFreqDependent_, isGminDependent_)
  {
    dtNodePtr_   = Teuchos::rcp(new specialsOp<usedType> (std::string("DT")));
    timeNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("TIME")));
    tempNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("TEMP")));
    vtNodePtr_   = Teuchos::rcp(new specialsOp<usedType> (std::string("VT")));
    freqNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("FREQ")));
    gminNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("GMIN")));
    piNodePtr_   = Teuchos::rcp(new piConstOp<usedType>  ());
    CtoKNodePtr_   = Teuchos::rcp(new CtoKConstOp<usedType>  ());

    std::string keyword = std::string("TABLE");
    Teuchos::RCP<tableOp<usedType> > tableNodePtr_ = Teuchos::RCP<tableOp<usedType> >(new tableOp<usedType> (keyword, left, xvals, yvals));
    astNodePtr_ = Teuchos::RCP<astNode<usedType> >(tableNodePtr_);

    if( !(Teuchos::is_null(left)) )
    {
      if (left->paramType())   { paramOpVec_.push_back(left); }
      if (left->funcType())    { funcOpVec_.push_back(left); }
      if (left->voltageType()) { voltOpVec_.push_back(left); }
      if (left->currentType()) { currentOpVec_.push_back(left); }
      if (left->leadCurrentType()) { leadCurrentOpVec_.push_back(left); }
      if (left->bsrcCurrentType()) { bsrcCurrentOpVec_.push_back(left); }
      if (left->powerType())   { powerOpVec_.push_back(left); }
      if (left->internalDeviceVarType()) { internalDevVarOpVec_.push_back(left); }

      if (left->dnoNoiseVarType()) { dnoNoiseDevVarOpVec_.push_back(left); }
      if (left->dniNoiseVarType()) { dniNoiseDevVarOpVec_.push_back(left); }
      if (left->oNoiseType()) { oNoiseOpVec_.push_back(left); }
      if (left->iNoiseType()) { iNoiseOpVec_.push_back(left); }
      if (left->sdtType()) { sdtOpVec_.push_back(left); }
      if (left->ddtType()) { ddtOpVec_.push_back(left); }
      if (left->srcType()) { srcAstNodeVec_.push_back(left); }
      if (left->stpType()) { stpAstNodeVec_.push_back(left); }
      if (left->compType()) { compAstNodeVec_.push_back(left); }
      if (left->limitType()) { limitAstNodeVec_.push_back(left); }
      if (left->phaseType()) { phaseOpVec_.push_back(left); }
      if (left->sparamType()) { sparamOpVec_.push_back(left); }
      if (left->yparamType()) { yparamOpVec_.push_back(left); }
      if (left->zparamType()) { zparamOpVec_.push_back(left); }
      if (left->agaussType()) { agaussOpVec_.push_back(left); }
      if (left->gaussType()) { gaussOpVec_.push_back(left); }
      if (left->aunifType()) { aunifOpVec_.push_back(left); }
      if (left->unifType()) { unifOpVec_.push_back(left); }
      if (left->randType()) { randOpVec_.push_back(left); }
      if (left->twoArgLimitType()) { twoArgLimitOpVec_.push_back(left); }

      left->getInterestingOps( opVectors_  );
    }
  };

  // copy constructor - this may need work
  newExpression (const newExpression & right) :
    group_(right.group_),
    expressionString_(right.expressionString_),
    parsed_(right.parsed_),
    derivsSetup_(right.derivsSetup_),
    astArraysSetup_(right.astArraysSetup_),
    functionArgStringVec_(right.functionArgStringVec_),
    functionArgOpVec_ (right.functionArgOpVec_),
    paramNameVec_(right.paramNameVec_),
    globalParamNameVec_(right.globalParamNameVec_),
    unresolvedParamNameVec_(right.unresolvedParamNameVec_),
    paramOpVec_(right.paramOpVec_),
    paramOpMap_(right.paramOpMap_),
    funcNameVec_(right.funcNameVec_),
    unresolvedFuncNameVec_(right.unresolvedFuncNameVec_),
    funcOpVec_(right.funcOpVec_),
    funcOpMap_(right.funcOpMap_),
    voltNameVec_(right.voltNameVec_),
    voltOpVec_(right.voltOpVec_),
    voltOpMap_(right.voltOpMap_),

    currentNameVec_(right.currentNameVec_),
    currentOpVec_(right.currentOpVec_),
    currentOpMap_(right.currentOpMap_),

    leadCurrentNameVec_(right.leadCurrentNameVec_),
    leadCurrentExcludeBsrcNameVec_(right.leadCurrentExcludeBsrcNameVec_),
    leadCurrentOpVec_(right.leadCurrentOpVec_),

    bsrcCurrentOpVec_(right.bsrcCurrentOpVec_),

    powerOpVec_(right.powerOpVec_),
    internalDevVarOpVec_(right.internalDevVarOpVec_),
    dnoNoiseDevVarOpVec_(right.dnoNoiseDevVarOpVec_),
    dniNoiseDevVarOpVec_(right.dniNoiseDevVarOpVec_),
    oNoiseOpVec_(right.oNoiseOpVec_),
    iNoiseOpVec_(right.iNoiseOpVec_),
    sdtOpVec_(right.sdtOpVec_),
    ddtOpVec_(right.ddtOpVec_),
    phaseOpVec_(right.phaseOpVec_),
    sparamOpVec_(right.sparamOpVec_),
    yparamOpVec_(right.yparamOpVec_),
    zparamOpVec_(right.zparamOpVec_),

    agaussOpVec_(right.agaussOpVec_),
    localAgaussOpVec_(right.localAgaussOpVec_),
    gaussOpVec_(right.gaussOpVec_),
    localGaussOpVec_(right.localGaussOpVec_),
    aunifOpVec_(right.aunifOpVec_),
    localAunifOpVec_(right.localAunifOpVec_),
    unifOpVec_(right.unifOpVec_),
    localUnifOpVec_(right.localUnifOpVec_),
    randOpVec_(right.randOpVec_),
    localRandOpVec_(right.localRandOpVec_),
    twoArgLimitOpVec_(right.twoArgLimitOpVec_),
    localTwoArgLimitOpVec_(right.localTwoArgLimitOpVec_),

    bpTol_(right.bpTol_),
    time_(right.time_),
    timeStep_(right.timeStep_),
    timeStepAlpha_(right.timeStepAlpha_),
    timeStepPrefac_(right.timeStepPrefac_),
    stepNumber_(right.stepNumber_),

    srcAstNodeVec_(right.srcAstNodeVec_),
    stpAstNodeVec_(right.stpAstNodeVec_),
    compAstNodeVec_(right.compAstNodeVec_),
    limitAstNodeVec_(right.limitAstNodeVec_),

    timeOpVec_(right.timeOpVec_),
    dtOpVec_(right.dtOpVec_),
    tempOpVec_(right.tempOpVec_),
    vtOpVec_(right.vtOpVec_),
    freqOpVec_(right.freqOpVec_),
    externalExpressions_(right.externalExpressions_),

    derivIndexVec_ (right.derivIndexVec_),
    derivNodeIndexMap_(right.derivNodeIndexMap_),
    numDerivs_(right.numDerivs_),
    traditionalParse_(right.traditionalParse_),
    externalDependencies_(right.externalDependencies_),
    isTimeDependent_(right.isTimeDependent_),
    isTempDependent_(right.isTempDependent_),
    isVTDependent_(right.isVTDependent_),
    isFreqDependent_(right.isFreqDependent_),
    isGminDependent_(right.isGminDependent_),

    isShallowTimeDependent_(right.isShallowTimeDependent_),
    isShallowTempDependent_(right.isShallowTempDependent_),
    isShallowVTDependent_(right.isShallowVTDependent_),
    isShallowFreqDependent_(right.isShallowFreqDependent_),
    isShallowGminDependent_(right.isShallowGminDependent_),

    isVariableDependent_(right.isVariableDependent_),
    isVoltageNodeDependent_(right.isVoltageNodeDependent_),
    isDeviceCurrentDependent_(right.isDeviceCurrentDependent_),
    isLeadCurrentDependent_(right.isLeadCurrentDependent_),
    isLeadCurrentDependentExcludeBsrc_(right.isLeadCurrentDependentExcludeBsrc_),

    overrideGroupTemperature_(right.overrideGroupTemperature_),
    overrideTemp_(right.overrideTemp_),
    isConstant_(right.isConstant_),
    evaluateFunctionCalledBefore_(right.evaluateFunctionCalledBefore_),
    evaluateCalledBefore_(right.evaluateCalledBefore_),
    getTheSeedCalledBefore_(right.getTheSeedCalledBefore_),
    savedResult_(right.savedResult_),
    phaseOutputUsesRadians_(right.phaseOutputUsesRadians_),
    opVectors_(paramOpVec_,funcOpVec_, voltOpVec_, currentOpVec_, leadCurrentOpVec_, bsrcCurrentOpVec_, powerOpVec_, internalDevVarOpVec_, dnoNoiseDevVarOpVec_, dniNoiseDevVarOpVec_, oNoiseOpVec_, iNoiseOpVec_, sdtOpVec_, ddtOpVec_, srcAstNodeVec_, stpAstNodeVec_, compAstNodeVec_, limitAstNodeVec_, phaseOpVec_, sparamOpVec_, yparamOpVec_, zparamOpVec_, agaussOpVec_, gaussOpVec_, aunifOpVec_, unifOpVec_, randOpVec_, twoArgLimitOpVec_, isTimeDependent_, isTempDependent_, isVTDependent_, isFreqDependent_, isGminDependent_)
  {
    dtNodePtr_   = right.dtNodePtr_;
    timeNodePtr_ = right.timeNodePtr_;
    tempNodePtr_ = right.tempNodePtr_;
    vtNodePtr_   = right.vtNodePtr_;
    freqNodePtr_ = right.freqNodePtr_;
    piNodePtr_   = right.piNodePtr_;
    CtoKNodePtr_ = right.CtoKNodePtr_;
    astNodePtr_ = right.astNodePtr_; // copy over the whole tree
  };

  // assignment operator
  newExpression & operator =(const newExpression & right)
  {
    group_ = right.group_;
    expressionString_ = right.expressionString_;
    parsed_ = right.parsed_;
    derivsSetup_ = right.derivsSetup_;
    astArraysSetup_ = right.astArraysSetup_;
    functionArgStringVec_ = right.functionArgStringVec_;
    functionArgOpVec_  = right.functionArgOpVec_;
    paramNameVec_ = right.paramNameVec_;
    globalParamNameVec_ = right.globalParamNameVec_;
    unresolvedParamNameVec_ = right.unresolvedParamNameVec_;
    paramOpVec_ = right.paramOpVec_;
    paramOpMap_ = right.paramOpMap_;
    unresolvedFuncNameVec_ = right.unresolvedFuncNameVec_;
    funcOpVec_ = right.funcOpVec_;
    funcOpMap_ = right.funcOpMap_;
    voltNameVec_ = right.voltNameVec_;
    voltOpVec_ = right.voltOpVec_;
    voltOpMap_ = right.voltOpMap_;
    currentNameVec_ = right.currentNameVec_;
    currentOpVec_ = right.currentOpVec_;
    currentOpMap_ = right.currentOpMap_;

    leadCurrentNameVec_ = right.leadCurrentNameVec_;
    leadCurrentExcludeBsrcNameVec_ = right.leadCurrentExcludeBsrcNameVec_;
    leadCurrentOpVec_ = right.leadCurrentOpVec_;

    bsrcCurrentOpVec_ = right.bsrcCurrentOpVec_;

    powerOpVec_ = right.powerOpVec_;
    internalDevVarOpVec_ = right.internalDevVarOpVec_;
    dnoNoiseDevVarOpVec_ = right.dnoNoiseDevVarOpVec_;
    dniNoiseDevVarOpVec_ = right.dniNoiseDevVarOpVec_;
    oNoiseOpVec_ = right.oNoiseOpVec_;
    iNoiseOpVec_ = right.iNoiseOpVec_;
    sdtOpVec_ = right.sdtOpVec_;
    ddtOpVec_ = right.ddtOpVec_;
    phaseOpVec_ = right.phaseOpVec_;
    sparamOpVec_ = right.sparamOpVec_;
    yparamOpVec_ = right.yparamOpVec_;
    zparamOpVec_ = right.zparamOpVec_;

    agaussOpVec_ = right.agaussOpVec_;
    localAgaussOpVec_ = right.localAgaussOpVec_;
    gaussOpVec_ = right.gaussOpVec_;
    localGaussOpVec_ = right.localGaussOpVec_;
    aunifOpVec_ = right.aunifOpVec_;
    localAunifOpVec_ = right.localAunifOpVec_;
    unifOpVec_ = right.unifOpVec_;
    localUnifOpVec_ = right.localUnifOpVec_;
    randOpVec_ = right.randOpVec_;
    localRandOpVec_ = right.localRandOpVec_;
    twoArgLimitOpVec_ = right.twoArgLimitOpVec_;
    localTwoArgLimitOpVec_ = right.localTwoArgLimitOpVec_;

    bpTol_ = right.bpTol_;
    time_ = right.time_;
    timeStep_ = right.timeStep_;
    timeStepAlpha_ = right.timeStepAlpha_;
    timeStepPrefac_ = right.timeStepPrefac_;
    stepNumber_ = right.stepNumber_;

    srcAstNodeVec_ = right.srcAstNodeVec_;
    stpAstNodeVec_ = right.stpAstNodeVec_;
    compAstNodeVec_ = right.compAstNodeVec_;
    limitAstNodeVec_ = right.limitAstNodeVec_;

    timeOpVec_ = right.timeOpVec_;
    dtOpVec_ = right.dtOpVec_;
    tempOpVec_ = right.tempOpVec_;
    vtOpVec_ = right.vtOpVec_;
    freqOpVec_ = right.freqOpVec_;
    externalExpressions_ = right.externalExpressions_;

    derivIndexVec_ = right.derivIndexVec_;
    derivNodeIndexMap_ = right.derivNodeIndexMap_;
    numDerivs_ = right.numDerivs_;
    traditionalParse_ = right.traditionalParse_;
    externalDependencies_ = right.externalDependencies_;

    isTimeDependent_ = right.isTimeDependent_;
    isTempDependent_ = right.isTempDependent_;
    isVTDependent_ = right.isVTDependent_;
    isFreqDependent_ = right.isFreqDependent_;
    isGminDependent_ = right.isGminDependent_;

    isShallowTimeDependent_ = right.isShallowTimeDependent_;
    isShallowTempDependent_ = right.isShallowTempDependent_;
    isShallowVTDependent_ = right.isShallowVTDependent_;
    isShallowFreqDependent_ = right.isShallowFreqDependent_;
    isShallowGminDependent_ = right.isShallowGminDependent_;

    isVariableDependent_ = right.isVariableDependent_;
    isVoltageNodeDependent_ = right.isVoltageNodeDependent_;
    isDeviceCurrentDependent_ = right.isDeviceCurrentDependent_;
    isLeadCurrentDependent_ = right.isLeadCurrentDependent_;
    isLeadCurrentDependentExcludeBsrc_ = right.isLeadCurrentDependentExcludeBsrc_;

    overrideGroupTemperature_ = right.overrideGroupTemperature_;
    overrideTemp_ = right.overrideTemp_;
    isConstant_ = right.isConstant_;
    evaluateFunctionCalledBefore_ = right.evaluateFunctionCalledBefore_;
    evaluateCalledBefore_ = right.evaluateCalledBefore_;
    getTheSeedCalledBefore_ = right.getTheSeedCalledBefore_;
    savedResult_ = right.savedResult_;
    phaseOutputUsesRadians_ = right.phaseOutputUsesRadians_;

    dtNodePtr_   = right.dtNodePtr_;
    timeNodePtr_ = right.timeNodePtr_;
    tempNodePtr_ = right.tempNodePtr_;
    vtNodePtr_   = right.vtNodePtr_;
    freqNodePtr_ = right.freqNodePtr_;
    piNodePtr_   = right.piNodePtr_;
    CtoKNodePtr_   = right.CtoKNodePtr_;
    astNodePtr_ = right.astNodePtr_; // copy over the whole tree

    return *this;
  };

  bool lexAndParseExpression();

  bool attachFunctionNode(const std::string & funcName, const Teuchos::RCP<Xyce::Util::newExpression> expPtr);
  bool attachParameterNode(const std::string & paramName, const Teuchos::RCP<Xyce::Util::newExpression> expPtr, enumParamType type=DOT_GLOBAL_PARAM);

  void clear(); // reset expression to the state it should be before lexAndParseExpression

  bool parsed() const { return parsed_; };
  bool derivsSetup () const { return derivsSetup_; };
  bool astArraysSetup () const { return astArraysSetup_; }

  bool make_constant (std::string const & var, usedType const & val, enumParamType type=DOT_GLOBAL_PARAM);
  bool make_var      (std::string const & var, usedType const & val, enumParamType type=DOT_GLOBAL_PARAM);

  void setAstPtr(Teuchos::RCP<astNode<usedType> > & astNodePtr) { astNodePtr_ = astNodePtr; };

  bool evaluate (usedType &result, std::vector< usedType > &derivs);
  bool evaluateFunction (usedType &result, bool efficiencyOn=false);

  // supporting "changed" boolean .... 
  void clearOldResult()  
  { 
    savedResult_ = 0.0; 
    evaluateFunctionCalledBefore_=false; 
  };

  void dumpParseTree(std::ostream & os) { if ( !(Teuchos::is_null(astNodePtr_)) ){astNodePtr_->output(os); }}

  void setupBreakPoints ();
  bool getBreakPoints (std::vector<Xyce::Util::BreakPoint> & breakPointTimes );

  Teuchos::RCP<astNode<usedType> > & getAst() {return astNodePtr_;}

  Teuchos::RCP<astNode<usedType> > getDtNode () { return dtNodePtr_; }
  Teuchos::RCP<astNode<usedType> > getTimeNode () { return timeNodePtr_; }
  Teuchos::RCP<astNode<usedType> > getTempNode () { return tempNodePtr_; }
  Teuchos::RCP<astNode<usedType> > getVtNode () { return vtNodePtr_; }
  Teuchos::RCP<astNode<usedType> > getFreqNode () { return freqNodePtr_; }
  Teuchos::RCP<astNode<usedType> > getGminNode () { return gminNodePtr_; }
  Teuchos::RCP<astNode<usedType> > getPiNode () { return piNodePtr_; }
  Teuchos::RCP<astNode<usedType> > getCtoKNode () { return CtoKNodePtr_; }

  // some of the parameter and function objects are stored in multiple containers.
  void setFunctionArgStringVec (const std::vector<std::string> & args);

  const std::vector<std::string> & getFunctionArgStringVec () { return functionArgStringVec_; };

  const std::vector< Teuchos::RCP<astNode<usedType> > > & getFunctionArgOpVec() { return functionArgOpVec_; };

  std::vector<Teuchos::RCP<astNode<usedType> > > & getParamOpVec () { return paramOpVec_; };
  std::vector<std::string> & getParamNameVec () { return paramNameVec_; };
  std::vector<std::string> & getGlobalParamNameVec () { return globalParamNameVec_; };
  std::vector<std::string> & getUnresolvedParamNameVec () { return unresolvedParamNameVec_; };

  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & getParamOpMap ()
  {
    setupVariousAstArrays ();
    return paramOpMap_;
  };

  std::vector<Teuchos::RCP<astNode<usedType> > > & getFuncOpVec () { return funcOpVec_; };
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & getFuncOpMap () { return funcOpMap_; };
  std::vector< std::string > & getFuncNameVec () { return funcNameVec_; };
  std::vector< std::string > & getUnresolvedFuncNameVec () { return unresolvedFuncNameVec_; };

  std::vector<Teuchos::RCP<astNode<usedType> > > & getVoltOpVec () { return voltOpVec_; };
  std::vector<std::string> & getVoltNameVec () { return voltNameVec_; };

  const std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & getVoltOpMap ()
  {
    setupVariousAstArrays ();
    return voltOpMap_;
  };

  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & getLocalVoltOpMap ()
  {
    return voltOpMap_;
  }

  std::vector<Teuchos::RCP<astNode<usedType> > > & getPowerOpVec() { return powerOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getInternalDevVarOpVec() { return internalDevVarOpVec_; }

  std::vector<Teuchos::RCP<astNode<usedType> > > & getDnoNoiseDevVarOpVec() { return dnoNoiseDevVarOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getDniNoiseDevVarOpVec() { return dniNoiseDevVarOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getONoiseOpVec() { return oNoiseOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getINoiseOpVec() { return iNoiseOpVec_; }

  std::vector<Teuchos::RCP<astNode<usedType> > > & getSdtOpVec() { return sdtOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getLocalSdtOpVec() { return localSdtOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getDdtOpVec() { return ddtOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getLocalDdtOpVec() { return localDdtOpVec_; }


  std::vector<Teuchos::RCP<astNode<usedType> > > & getPhaseOpVec() { return phaseOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getSparamOpVec() { return sparamOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getYparamOpVec() { return yparamOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getZparamOpVec() { return zparamOpVec_; }

  std::vector<Teuchos::RCP<astNode<usedType> > > & getAgaussOpVec() { return agaussOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getGaussOpVec() { return gaussOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getAunifOpVec() { return aunifOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getUnifOpVec() { return unifOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getRandOpVec() { return randOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getTwoArgLimitOpVec() { return twoArgLimitOpVec_; }

  std::vector<Teuchos::RCP<astNode<usedType> > > & getLocalAgaussOpVec() { return localAgaussOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getLocalGaussOpVec() { return localGaussOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getLocalAunifOpVec() { return localAunifOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getLocalUnifOpVec() { return localUnifOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getLocalRandOpVec() { return localRandOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getLocalTwoArgLimitOpVec() { return localTwoArgLimitOpVec_; }

  std::vector<Teuchos::RCP<astNode<usedType> > > & getTimeOpVec() { return timeOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getDtOpVec() { return dtOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getTempOpVec() { return tempOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getVtOpVec() { return vtOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getFreqOpVec() { return freqOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getGminOpVec() { return gminOpVec_; }

  std::vector<Teuchos::RCP<astNode<usedType> > > & getCurrentOpVec () { return currentOpVec_; };
  std::vector<std::string> & getCurrentNameVec () { return currentNameVec_; };
  const std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & getCurrentOpMap ()
  {
    setupVariousAstArrays ();
    return currentOpMap_;
  };

  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & getLocalCurrentOpMap ()
  {
    return currentOpMap_;
  };

  std::vector<std::string> & getLeadCurrentNameVec () { return leadCurrentNameVec_; }
  std::vector<std::string> & getLeadCurrentExcludeBsrcNameVec () { return leadCurrentExcludeBsrcNameVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getLeadCurrentOpVec () { return leadCurrentOpVec_; };

  std::vector<Teuchos::RCP<astNode<usedType> > > & getBsrcCurrentOpVec () { return bsrcCurrentOpVec_; };

  void codeGen( std::ostream & os )
  {
    if ( !(Teuchos::is_null(astNodePtr_)) )
    {
      astNodePtr_->codeGen(os); os << ";" <<std::endl;
    }
  };

  std::vector< Teuchos::RCP<astNode<usedType> > > & getSrcNodeVec() { return srcAstNodeVec_;}
  std::vector< Teuchos::RCP<astNode<usedType> > > & getStpNodeVec() { return stpAstNodeVec_;}
  std::vector< Teuchos::RCP<astNode<usedType> > > & getCompNodeVec() { return compAstNodeVec_;}
  std::vector< Teuchos::RCP<astNode<usedType> > > & getLimitNodeVec() { return limitAstNodeVec_;}

  const std::string & getExpressionString() { return expressionString_; };

  bool replaceName ( const std::string & old_name, const std::string & new_name);

  double getTime() { return std::real(timeNodePtr_->val()); };

  bool getTimeDependent() { return isTimeDependent_; }
  void setTimeDependent(bool val) { isTimeDependent_ = val; }

  bool getTempDependent() { return isTempDependent_; }
  void setTempDependent(bool val) { isTempDependent_ = val; }

  bool getVTDependent() { return isVTDependent_; }
  void setVTDependent(bool val) { isVTDependent_ = val; }

  bool getFreqDependent() { return isFreqDependent_; }
  void setFreqDependent(bool val) { isFreqDependent_ = val; }

  bool getGminDependent() { return isGminDependent_; }
  void setGminDependent(bool val) { isGminDependent_ = val; }

  bool getShallowTimeDependent() { return isShallowTimeDependent_; }
  void setShallowTimeDependent(bool val) { isShallowTimeDependent_ = val; }

  bool getShallowTempDependent() { return isShallowTempDependent_; }
  void setShallowTempDependent(bool val) { isShallowTempDependent_ = val; }

  bool getShallowVTDependent() { return isShallowVTDependent_; }
  void setShallowVTDependent(bool val) { isShallowVTDependent_ = val; }

  bool getShallowFreqDependent() { return isShallowFreqDependent_; }
  void setShallowFreqDependent(bool val) { isShallowFreqDependent_ = val; }

  bool getShallowGminDependent() { return isShallowGminDependent_; }
  void setShallowGminDependent(bool val) { isShallowGminDependent_ = val; }

  bool getVariableDependent() { return isVariableDependent_; }
  bool getVoltageNodeDependent() { return isVoltageNodeDependent_; }
  bool getDeviceCurrentDependent() { return isDeviceCurrentDependent_; }
  bool getLeadCurrentDependent() { return isLeadCurrentDependent_; }
  bool getLeadCurrentDependentExcludeBsrc() { return isLeadCurrentDependentExcludeBsrc_; }

  // this function is only used to determine function arguments of a function prototype
  // So if we have .func abc(x,y) {x+y+10}
  // At a certain point the prototype abc(x,y) will get parsed, and x,y are the prototype args.
  //
  // The old Xyce source code uses an expression to parse abc(x,y) so I'm attempting to
  // do the same here. When doing this, it is necessary to keep the args in order,
  // which is why I can't rely on the getParamOpVec function (as the paramOpVec is in
  // the same order as the paramOpMap)
  void getFuncPrototypeArgStrings ( std::vector< std::string > & funcArgStrings )
  {
    funcArgStrings.clear();
    if(!(funcOpVec_.empty()))
    {
      Teuchos::RCP<funcOp<usedType> > funcOperator = Teuchos::rcp_dynamic_cast<funcOp<usedType> > (funcOpVec_[0]) ;
      if ( !(Teuchos::is_null(funcOperator)) )
      {
        std::vector< Teuchos::RCP<astNode<usedType> > > & funcArgs = funcOperator->getFuncArgs();
        for(int ii=0;ii<funcArgs.size();++ii) { funcArgStrings.push_back(funcArgs[ii]->getName()); }
      }
    }
    return;
  }

  // when parsing the function prototype, the function name is needed as well.
  void getFuncPrototypeName ( std::string & prototypeName)
  {
    if(!(funcOpVec_.empty()))
    {
      prototypeName = funcOpVec_[0]->getName();
    }
  }

  void outputVariousAstArrays     ( std::ostream & os );
  void outputVariousAstArraySizes ( std::ostream & os );


  // "expression" traversal functions, as opposed to AST traversals.
  // Expression traverals make more sense for tracking specials (time, dt, temp, vt, freq, gmin), as each
  // expression object will have zero or one allocation of each special.
  // (which really should be 100% singletons, but for the time being are not)
  void getTimeNodes( std::vector<Teuchos::RCP<astNode<usedType> > > & timeVec)
  {
    if (!(timeOpVec_.empty())) { timeVec.push_back(timeOpVec_[0]); }
    for (int ii=0;ii<externalExpressions_.size();ii++) { externalExpressions_[ii]->getTimeNodes(timeVec); }
  }

  void getDtNodes( std::vector<Teuchos::RCP<astNode<usedType> > > & dtVec)
  {
    if (!(dtOpVec_.empty())) { dtVec.push_back(dtOpVec_[0]); }
    for (int ii=0;ii<externalExpressions_.size();ii++) { externalExpressions_[ii]->getDtNodes(dtVec); }
  }

  void getTempNodes( std::vector<Teuchos::RCP<astNode<usedType> > > & tempVec)
  {
    if (!(tempOpVec_.empty())) { tempVec.push_back(tempOpVec_[0]); }
    for (int ii=0;ii<externalExpressions_.size();ii++) { externalExpressions_[ii]->getTempNodes(tempVec); }
  }

  void getVtNodes( std::vector<Teuchos::RCP<astNode<usedType> > > & vtVec)
  {
    if (!(vtOpVec_.empty())) { vtVec.push_back(vtOpVec_[0]); }
    for (int ii=0;ii<externalExpressions_.size();ii++) { externalExpressions_[ii]->getVtNodes(vtVec); }
  }

  void getFreqNodes( std::vector<Teuchos::RCP<astNode<usedType> > > & freqVec)
  {
    if (!(freqOpVec_.empty())) { freqVec.push_back(freqOpVec_[0]); }
    for (int ii=0;ii<externalExpressions_.size();ii++) { externalExpressions_[ii]->getFreqNodes(freqVec); }
  }

  void getGminNodes( std::vector<Teuchos::RCP<astNode<usedType> > > & gminVec)
  {
    if (!(gminOpVec_.empty())) { gminVec.push_back(gminOpVec_[0]); }
    for (int ii=0;ii<externalExpressions_.size();ii++) { externalExpressions_[ii]->getGminNodes(gminVec); }
  }

  bool getIsConstant() 
  { 
    setupVariousAstArrays ();
    return isConstant_; 
  }

  void treatAsTempAndConvert();

  bool setTemperature (const double & temp);

  static void clearProcessSuccessfulTimeStepMap () { staticsContainer::processSuccessfulStepMap.clear(); }

  void processSuccessfulTimeStep ();

  int getNumDdt () { return ddtOpVec_.size(); }
  void getDdtVals   (std::vector<double> & vals);
  void getDdtVals   (std::vector<std::complex<double> > & vals);
  void setDdtDerivs (std::vector<double> & vals);
  void setDdtDerivs (std::vector<std::complex<double> > & vals);

  void setupVariousAstArrays ();

  void setGroup( Teuchos::RCP<baseExpressionGroup> & grp ) { group_ = grp; }
  Teuchos::RCP<baseExpressionGroup> getGroup() { return group_; }

private:
  void setupDerivatives_ ();
  void checkIsConstant_();
  bool getValuesFromGroup_();

  Teuchos::RCP<baseExpressionGroup> group_;
  std::string expressionString_;
  bool parsed_;
  bool derivsSetup_;
  bool astArraysSetup_;

  Teuchos::RCP<astNode<usedType> > astNodePtr_;

  // function argument vectors.  Both strings and operators
  std::vector< std::string > functionArgStringVec_;

  std::vector< Teuchos::RCP<astNode<usedType> > > functionArgOpVec_;

  std::vector<std::string> paramNameVec_;
  std::vector<std::string> globalParamNameVec_;
  std::vector<std::string> unresolvedParamNameVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > paramOpVec_; // mainly used for .func argument setup
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > paramOpMap_; // used by make_const, make_var and attachParamNode

  std::vector<std::string> funcNameVec_;
  std::vector<std::string> unresolvedFuncNameVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > funcOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > funcOpMap_; // used by attachFunctionNode

  std::vector<std::string> voltNameVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > voltOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > voltOpMap_; // used by replaceName

  std::vector<std::string> currentNameVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > currentOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > currentOpMap_; // used by replaceName

  std::vector<std::string> leadCurrentNameVec_;
  std::vector<std::string> leadCurrentExcludeBsrcNameVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > leadCurrentOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > bsrcCurrentOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > powerOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > internalDevVarOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > dnoNoiseDevVarOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > dniNoiseDevVarOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > oNoiseOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > iNoiseOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > sdtOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > localSdtOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > ddtOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > localDdtOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > phaseOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > sparamOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > yparamOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > zparamOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > agaussOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > localAgaussOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > gaussOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > localGaussOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > aunifOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > localAunifOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > unifOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > localUnifOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > randOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > localRandOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > twoArgLimitOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > localTwoArgLimitOpVec_;

  // time integration related variables
  double bpTol_;
  double startingTimeStep_;
  double finalTime_;

  double time_;
  double timeStep_;
  double timeStepAlpha_;
  double timeStepPrefac_;
  unsigned int stepNumber_;

  // vector of independent sources, but only those with breakpoints.  This
  // vector is ONLY used for obtaining breakpoints.
  std::vector< Teuchos::RCP<astNode<usedType> > > srcAstNodeVec_;

  // vector of STP objects.  Needed for breakpoints.  
  std::vector< Teuchos::RCP<astNode<usedType> > > stpAstNodeVec_;

  // vector of COMP objects.  Needed for breakpoints.  
  std::vector< Teuchos::RCP<astNode<usedType> > > compAstNodeVec_;

  // vector of LIMIT objects.  Needed for breakpoints.  
  std::vector< Teuchos::RCP<astNode<usedType> > > limitAstNodeVec_;

  // const and specials nodes:
  Teuchos::RCP<specialsOp<usedType> > dtNodePtr_;
  Teuchos::RCP<specialsOp<usedType> > timeNodePtr_;
  Teuchos::RCP<specialsOp<usedType> > tempNodePtr_;
  Teuchos::RCP<specialsOp<usedType> > vtNodePtr_;
  Teuchos::RCP<specialsOp<usedType> > freqNodePtr_;
  Teuchos::RCP<specialsOp<usedType> > gminNodePtr_;
  Teuchos::RCP<piConstOp<usedType> > piNodePtr_;
  Teuchos::RCP<CtoKConstOp<usedType> > CtoKNodePtr_;

  // to handle externally attached expressions, which have specials dependence, we need vectors of these things.
  std::vector<Teuchos::RCP<astNode<usedType> > > timeOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > dtOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > tempOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > vtOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > freqOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > gminOpVec_;

  std::vector<Teuchos::RCP<Xyce::Util::newExpression> > externalExpressions_;

  // derivative related stuff
  typedef  std::pair< Teuchos::RCP<astNode<usedType> > , int > derivIndexPair_;
  std::vector< derivIndexPair_  > derivIndexVec_;

  std::unordered_map<std::string,int> derivNodeIndexMap_;
  int numDerivs_;

  bool traditionalParse_;
  bool externalDependencies_; // true if expression includes a call to a .func, .param or .global_param

  bool isTimeDependent_;
  bool isTempDependent_;
  bool isVTDependent_;
  bool isFreqDependent_;
  bool isGminDependent_;

  bool isShallowTimeDependent_;
  bool isShallowTempDependent_;
  bool isShallowVTDependent_;
  bool isShallowFreqDependent_;
  bool isShallowGminDependent_;

  bool isVariableDependent_;
  bool isVoltageNodeDependent_;
  bool isDeviceCurrentDependent_;
  bool isLeadCurrentDependent_;
  bool isLeadCurrentDependentExcludeBsrc_;

  bool overrideGroupTemperature_;
  double overrideTemp_;

  bool isConstant_;
  bool evaluateFunctionCalledBefore_;
  bool evaluateCalledBefore_;
  bool getTheSeedCalledBefore_;
  usedType savedResult_;
  bool phaseOutputUsesRadians_;

  opVectorContainers<usedType> opVectors_;
  std::vector<usedType> oldSolVals_;
};

}
}

#endif

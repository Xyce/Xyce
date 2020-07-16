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
    astNodePtrPtr_(NULL),
    tableNodePtrPtr_(NULL),
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
    overrideGroupTemperature_(false),
    isConstant_(false),
    evaluateFunctionCalledBefore_(false),
    evaluateCalledBefore_(false),
    savedResult_(0.0),
    phaseOutputUsesRadians_(false),
    opVectors_(paramOpVec_,funcOpVec_, voltOpVec_, currentOpVec_, leadCurrentOpVec_, bsrcCurrentOpVec_, powerOpVec_, internalDevVarOpVec_, dnoNoiseDevVarOpVec_, dniNoiseDevVarOpVec_, oNoiseOpVec_, iNoiseOpVec_, sdtOpVec_, ddtOpVec_, srcAstNodeVec_, stpAstNodeVec_, compAstNodeVec_, phaseOpVec_, sparamOpVec_, yparamOpVec_, zparamOpVec_, isTimeDependent_, isTempDependent_, isVTDependent_, isFreqDependent_, isGminDependent_)
  {};

  // primary constructor
  newExpression ( std::string const & exp, Teuchos::RCP<baseExpressionGroup> & group ) :
    group_(group),
    expressionString_(exp),
    originalExpressionString_(exp),
    parsed_(false),
    derivsSetup_(false),
    astArraysSetup_(false),
    astNodePtrPtr_(NULL),
    tableNodePtrPtr_(NULL),
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
    overrideGroupTemperature_(false),
    isConstant_(false),
    evaluateFunctionCalledBefore_(false),
    evaluateCalledBefore_(false),
    savedResult_(0.0),
    phaseOutputUsesRadians_(false),
    opVectors_(paramOpVec_,funcOpVec_, voltOpVec_, currentOpVec_, leadCurrentOpVec_, bsrcCurrentOpVec_, powerOpVec_, internalDevVarOpVec_, dnoNoiseDevVarOpVec_, dniNoiseDevVarOpVec_, oNoiseOpVec_, iNoiseOpVec_, sdtOpVec_, ddtOpVec_, srcAstNodeVec_, stpAstNodeVec_, compAstNodeVec_, phaseOpVec_, sparamOpVec_, yparamOpVec_, zparamOpVec_, isTimeDependent_, isTempDependent_, isVTDependent_, isFreqDependent_, isGminDependent_)
  {
    // The bison file is officially case-insensitive.  So converting the
    // input string to all upper case is not necessary for it to work.
    //
    // However:
    //
    // Xyce mostly deals with netlist strings by converting to upper case.
    // So, when the code outside of bison code interacts with it, Bison will not
    // have converted it to upper or lower case.  It will be the original case,
    // whatever that was in the netlist.
    //
    // The simplest, easiest way to make all of this work is to simply
    // upcase the whole string.
    //
    Xyce::Util::toUpper(expressionString_);

    garbageParamOpPtr_ = Teuchos::rcp(new paramOp<usedType> (std::string("GARBAGE")));

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
    originalExpressionString_("TIME"),
    parsed_(false),
    derivsSetup_(false),
    astArraysSetup_(false),
    astNodePtrPtr_(NULL),
    tableNodePtrPtr_(NULL),
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
    overrideGroupTemperature_(false),
    isConstant_(false),
    evaluateFunctionCalledBefore_(false),
    evaluateCalledBefore_(false),
    savedResult_(0.0),
    phaseOutputUsesRadians_(false),
    opVectors_(paramOpVec_,funcOpVec_, voltOpVec_, currentOpVec_, leadCurrentOpVec_, bsrcCurrentOpVec_, powerOpVec_, internalDevVarOpVec_, dnoNoiseDevVarOpVec_, dniNoiseDevVarOpVec_, oNoiseOpVec_, iNoiseOpVec_, sdtOpVec_, ddtOpVec_, srcAstNodeVec_, stpAstNodeVec_, compAstNodeVec_, phaseOpVec_, sparamOpVec_, yparamOpVec_, zparamOpVec_, isTimeDependent_, isTempDependent_, isVTDependent_, isFreqDependent_, isGminDependent_)
  {
    garbageParamOpPtr_ = Teuchos::rcp(new paramOp<usedType> (std::string("GARBAGE")));

    dtNodePtr_   = Teuchos::rcp(new specialsOp<usedType> (std::string("DT")));
    timeNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("TIME")));
    tempNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("TEMP")));
    vtNodePtr_   = Teuchos::rcp(new specialsOp<usedType> (std::string("VT")));
    freqNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("FREQ")));
    gminNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("GMIN")));
    piNodePtr_   = Teuchos::rcp(new piConstOp<usedType>  ());
    CtoKNodePtr_   = Teuchos::rcp(new CtoKConstOp<usedType>  ());

    Teuchos::RCP<astNode<usedType> > time_base = timeNodePtr_;
    tableNodePtrPtr_ = new Teuchos::RCP<tableOp<usedType> >(new tableOp<usedType> (time_base, xvals, yvals));
    astNodePtrPtr_ = new Teuchos::RCP<astNode<usedType> >(*tableNodePtrPtr_);
    setAstPtr(*astNodePtrPtr_);
  };

  // another "big table" constructor
  newExpression (Teuchos::RCP<astNode<usedType> > & left,
      const std::vector<usedType> & xvals, const std::vector<usedType> & yvals,
      Teuchos::RCP<baseExpressionGroup> & group ) :
    group_(group),
    expressionString_("TIME"),
    originalExpressionString_("TIME"),
    parsed_(false),
    derivsSetup_(false),
    astArraysSetup_(false),
    astNodePtrPtr_(NULL),
    tableNodePtrPtr_(NULL),
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
    overrideGroupTemperature_(false),
    isConstant_(false),
    evaluateFunctionCalledBefore_(false),
    evaluateCalledBefore_(false),
    savedResult_(0.0),
    phaseOutputUsesRadians_(false),
    opVectors_(paramOpVec_,funcOpVec_, voltOpVec_, currentOpVec_, leadCurrentOpVec_, bsrcCurrentOpVec_, powerOpVec_, internalDevVarOpVec_, dnoNoiseDevVarOpVec_, dniNoiseDevVarOpVec_, oNoiseOpVec_, iNoiseOpVec_, sdtOpVec_, ddtOpVec_, srcAstNodeVec_, stpAstNodeVec_, compAstNodeVec_, phaseOpVec_, sparamOpVec_, yparamOpVec_, zparamOpVec_, isTimeDependent_, isTempDependent_, isVTDependent_, isFreqDependent_, isGminDependent_)
  {
    garbageParamOpPtr_ = Teuchos::rcp(new paramOp<usedType> (std::string("GARBAGE")));

    dtNodePtr_   = Teuchos::rcp(new specialsOp<usedType> (std::string("DT")));
    timeNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("TIME")));
    tempNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("TEMP")));
    vtNodePtr_   = Teuchos::rcp(new specialsOp<usedType> (std::string("VT")));
    freqNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("FREQ")));
    gminNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("GMIN")));
    piNodePtr_   = Teuchos::rcp(new piConstOp<usedType>  ());
    CtoKNodePtr_   = Teuchos::rcp(new CtoKConstOp<usedType>  ());

    tableNodePtrPtr_ = new Teuchos::RCP<tableOp<usedType> >(new tableOp<usedType> (left, xvals, yvals));
    astNodePtrPtr_ = new Teuchos::RCP<astNode<usedType> >(*tableNodePtrPtr_);
    setAstPtr(*astNodePtrPtr_);

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
      if (left->phaseType()) { phaseOpVec_.push_back(left); }
      if (left->sparamType()) { sparamOpVec_.push_back(left); }
      if (left->yparamType()) { yparamOpVec_.push_back(left); }
      if (left->zparamType()) { zparamOpVec_.push_back(left); }

      left->getInterestingOps( opVectors_  );
    }
  };

  // copy constructor - this may need work
  // This is necessary b/c of things like the masterAstNodeVec_ object, which uses raw pointers
  newExpression (const newExpression & right) :
    group_(right.group_),
    expressionString_(right.expressionString_),
    originalExpressionString_(right.originalExpressionString_),
    parsed_(right.parsed_),
    derivsSetup_(right.derivsSetup_),
    astArraysSetup_(right.astArraysSetup_),
    astNodePtrPtr_(NULL),
    tableNodePtrPtr_(NULL),
    functionArgStringVec_(right.functionArgStringVec_),
    functionArgOpVec_ (right.functionArgOpVec_),
    paramOpVec_(right.paramOpVec_),
    unresolvedParamOpVec_(right.unresolvedParamOpVec_),
    paramOpNames_(right.paramOpNames_),
    funcOpVec_(right.funcOpVec_),
    unresolvedFuncOpVec_(right.unresolvedFuncOpVec_),
    voltOpVec_(right.voltOpVec_),
    unresolvedVoltOpVec_(right.unresolvedVoltOpVec_),
    voltOpNames_(right.voltOpNames_),

    currentOpVec_(right.currentOpVec_),
    unresolvedCurrentOpVec_(right.unresolvedCurrentOpVec_),
    currentOpNames_(right.currentOpNames_),

    leadCurrentOpVec_(right.leadCurrentOpVec_),
    unresolvedLeadCurrentOpVec_(right.unresolvedLeadCurrentOpVec_),
    leadCurrentOpNames_(right.leadCurrentOpNames_),

    bsrcCurrentOpVec_(right.bsrcCurrentOpVec_),
    unresolvedBsrcCurrentOpVec_(right.unresolvedBsrcCurrentOpVec_),
    bsrcCurrentOpNames_(right.bsrcCurrentOpNames_),

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

    bpTol_(right.bpTol_),
    time_(right.time_),
    timeStep_(right.timeStep_),
    timeStepAlpha_(right.timeStepAlpha_),
    timeStepPrefac_(right.timeStepPrefac_),
    stepNumber_(right.stepNumber_),

    srcAstNodeVec_(right.srcAstNodeVec_),
    stpAstNodeVec_(right.stpAstNodeVec_),
    compAstNodeVec_(right.compAstNodeVec_),

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
    overrideGroupTemperature_(right.overrideGroupTemperature_),
    isConstant_(right.isConstant_),
    evaluateFunctionCalledBefore_(right.evaluateFunctionCalledBefore_),
    evaluateCalledBefore_(right.evaluateCalledBefore_),
    savedResult_(right.savedResult_),
    phaseOutputUsesRadians_(right.phaseOutputUsesRadians_),
    opVectors_(paramOpVec_,funcOpVec_, voltOpVec_, currentOpVec_, leadCurrentOpVec_, bsrcCurrentOpVec_, powerOpVec_, internalDevVarOpVec_, dnoNoiseDevVarOpVec_, dniNoiseDevVarOpVec_, oNoiseOpVec_, iNoiseOpVec_, sdtOpVec_, ddtOpVec_, srcAstNodeVec_, stpAstNodeVec_, compAstNodeVec_, phaseOpVec_, sparamOpVec_, yparamOpVec_, zparamOpVec_, isTimeDependent_, isTempDependent_, isVTDependent_, isFreqDependent_, isGminDependent_)
  {
    garbageParamOpPtr_ = right.garbageParamOpPtr_;
    dtNodePtr_   = right.dtNodePtr_;
    timeNodePtr_ = right.timeNodePtr_;
    tempNodePtr_ = right.tempNodePtr_;
    vtNodePtr_   = right.vtNodePtr_;
    freqNodePtr_ = right.freqNodePtr_;
    piNodePtr_   = right.piNodePtr_;
    CtoKNodePtr_ = right.CtoKNodePtr_;
    astNodePtr_ = right.astNodePtr_; // copy over the whole tree

    // pointers to RCP's. use constructors
    if (right.astNodePtrPtr_)
    {
      astNodePtrPtr_ = new Teuchos::RCP<astNode<usedType> >(*(right.astNodePtrPtr_));
    }
    if (right.tableNodePtrPtr_)
    {
      tableNodePtrPtr_ = new Teuchos::RCP< tableOp<usedType> >(*(right.tableNodePtrPtr_));
    }

    // vectors of pointers to RCPs of AST nodes; use constructors
    for (int ii=0;ii<right.masterAstNodeVec_.size();ii++)
    {
      masterAstNodeVec_.push_back( new Teuchos::RCP<astNode<usedType> >(*(right.masterAstNodeVec_[ii])) );
    }
  };

  // assignment operator
  // This is necessary b/c of things like the masterAstNodeVec_ object, which uses raw pointers
  newExpression & operator =(const newExpression & right)
  {
    group_ = right.group_;
    expressionString_ = right.expressionString_;
    originalExpressionString_ = right.originalExpressionString_;
    parsed_ = right.parsed_;
    derivsSetup_ = right.derivsSetup_;
    astArraysSetup_ = right.astArraysSetup_;
    astNodePtrPtr_ = NULL;
    tableNodePtrPtr_ = NULL;
    functionArgStringVec_ = right.functionArgStringVec_;
    functionArgOpVec_  = right.functionArgOpVec_;
    paramOpVec_ = right.paramOpVec_;
    unresolvedParamOpVec_ = right.unresolvedParamOpVec_;
    paramOpNames_ = right.paramOpNames_;
    funcOpVec_ = right.funcOpVec_;
    unresolvedFuncOpVec_ = right.unresolvedFuncOpVec_;
    voltOpVec_ = right.voltOpVec_;
    voltOpNames_ = right.voltOpNames_;
    unresolvedVoltOpVec_ = right.unresolvedVoltOpVec_;
    currentOpVec_ = right.currentOpVec_;
    currentOpNames_ = right.currentOpNames_;
    unresolvedCurrentOpVec_ = right.unresolvedCurrentOpVec_;

    leadCurrentOpVec_ = right.leadCurrentOpVec_;
    leadCurrentOpNames_ = right.leadCurrentOpNames_;
    unresolvedLeadCurrentOpVec_ = right.unresolvedLeadCurrentOpVec_;

    bsrcCurrentOpVec_ = right.bsrcCurrentOpVec_;
    bsrcCurrentOpNames_ = right.bsrcCurrentOpNames_;
    unresolvedBsrcCurrentOpVec_ = right.unresolvedBsrcCurrentOpVec_;

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

    bpTol_ = right.bpTol_;
    time_ = right.time_;
    timeStep_ = right.timeStep_;
    timeStepAlpha_ = right.timeStepAlpha_;
    timeStepPrefac_ = right.timeStepPrefac_;
    stepNumber_ = right.stepNumber_;

    srcAstNodeVec_ = right.srcAstNodeVec_;
    stpAstNodeVec_ = right.stpAstNodeVec_;
    compAstNodeVec_ = right.compAstNodeVec_;

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

    overrideGroupTemperature_ = right.overrideGroupTemperature_;
    isConstant_ = right.isConstant_;
    evaluateFunctionCalledBefore_ = right.evaluateFunctionCalledBefore_;
    evaluateCalledBefore_ = right.evaluateCalledBefore_;
    savedResult_ = right.savedResult_;
    phaseOutputUsesRadians_ = right.phaseOutputUsesRadians_;

    garbageParamOpPtr_ = right.garbageParamOpPtr_;
    dtNodePtr_   = right.dtNodePtr_;
    timeNodePtr_ = right.timeNodePtr_;
    tempNodePtr_ = right.tempNodePtr_;
    vtNodePtr_   = right.vtNodePtr_;
    freqNodePtr_ = right.freqNodePtr_;
    piNodePtr_   = right.piNodePtr_;
    CtoKNodePtr_   = right.CtoKNodePtr_;
    astNodePtr_ = right.astNodePtr_; // copy over the whole tree

    // pointers to RCP's. use constructors
    if (right.astNodePtrPtr_)
    {
      astNodePtrPtr_ = new Teuchos::RCP<astNode<usedType> >(*(right.astNodePtrPtr_));
    }
    if (right.tableNodePtrPtr_)
    {
      tableNodePtrPtr_ = new Teuchos::RCP< tableOp<usedType> >(*(right.tableNodePtrPtr_));
    }

    // vectors of pointers to RCPs of AST nodes; use constructors
    for (int ii=0;ii<right.masterAstNodeVec_.size();ii++)
    {
      masterAstNodeVec_.push_back( new Teuchos::RCP<astNode<usedType> >(*(right.masterAstNodeVec_[ii])) );
    }

    return *this;
  };

  // destructor
  ~newExpression ()
  {
    // pointers to RCP's. use constructors
    if (astNodePtrPtr_)
    {
      delete astNodePtrPtr_;
    }
    if ( tableNodePtrPtr_ )
    {
      delete tableNodePtrPtr_;
    }

    // vectors of pointers to RCPs of AST nodes; use constructors
    for (int ii=0;ii<masterAstNodeVec_.size();ii++)
    {
      delete masterAstNodeVec_[ii];
    }
    masterAstNodeVec_.clear();
  };

  bool lexAndParseExpression();

  bool attachFunctionNode(const std::string & funcName, const Teuchos::RCP<Xyce::Util::newExpression> expPtr);
  bool attachParameterNode(const std::string & paramName, const Teuchos::RCP<Xyce::Util::newExpression> expPtr, bool isDotParam=false);

  void clear(); // reset expression to the state it should be before lexAndParseExpression

  bool parsed() const { return parsed_; };
  bool derivsSetup () const { return derivsSetup_; };
  bool astArraysSetup () const { return astArraysSetup_; }

  bool make_constant (std::string const & var, usedType const & val, bool isDotParam=false);
  bool make_var (std::string const & var, bool isDotParam=false);

  void setAstPtr(Teuchos::RCP<astNode<usedType> > & astNodePtr) { astNodePtr_ = astNodePtr; };

  // these two functions return int error codes in the original expression library
  int evaluate (usedType &result, std::vector< usedType > &derivs);
  int evaluateFunction (usedType &result);

  void dumpParseTree(std::ostream & os) { if ( !(Teuchos::is_null(astNodePtr_)) ){astNodePtr_->output(os); }}

  bool getBreakPoints (std::vector<Xyce::Util::BreakPoint> & breakPointTimes );

  Teuchos::RCP<astNode<usedType> > & getAst() {return astNodePtr_;}

  Teuchos::RCP<paramOp<usedType> > & getGarbageParam() {return garbageParamOpPtr_;}

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
  std::vector<Teuchos::RCP<astNode<usedType> > > & getUnresolvedParamOpVector() {  return unresolvedParamOpVec_; };
  std::vector<std::string> & getParamNameVec () { return paramNameVec_; };

  const std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & getParamOpNames ()
  {
    if (!astArraysSetup_) { setupVariousAstArrays_ (); }
    return paramOpNames_;
  };

  std::vector<Teuchos::RCP<astNode<usedType> > > & getFuncOpVec () { return funcOpVec_; };
  std::vector<Teuchos::RCP<astNode<usedType> > > & getUnresolvedFuncOpVec() { return unresolvedFuncOpVec_; };
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & getFuncOpMap () { return funcOpMap_; };
  std::vector< std::string > & getFuncNameVec () { return funcNameVec_; };

  std::vector<Teuchos::RCP<astNode<usedType> > > & getVoltOpVec () { return voltOpVec_; };
  std::vector<Teuchos::RCP<astNode<usedType> > > & getUnresolvedVoltOpVec() { return unresolvedVoltOpVec_; };
  const std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & getVoltOpNames ()
  {
    if (!astArraysSetup_) { setupVariousAstArrays_ (); }
    return voltOpNames_;
  };

  std::vector<Teuchos::RCP<astNode<usedType> > > & getPowerOpVec() { return powerOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getInternalDevVarOpVec() { return internalDevVarOpVec_; }

  std::vector<Teuchos::RCP<astNode<usedType> > > & getDnoNoiseDevVarOpVec() { return dnoNoiseDevVarOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getDniNoiseDevVarOpVec() { return dniNoiseDevVarOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getONoiseOpVec() { return oNoiseOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getINoiseOpVec() { return iNoiseOpVec_; }

  std::vector<Teuchos::RCP<astNode<usedType> > > & getSdtOpVec() { return sdtOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getDdtOpVec() { return ddtOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getPhaseOpVec() { return phaseOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getSparamOpVec() { return sparamOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getYparamOpVec() { return yparamOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getZparamOpVec() { return zparamOpVec_; }

  std::vector<Teuchos::RCP<astNode<usedType> > > & getTimeOpVec() { return timeOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getDtOpVec() { return dtOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getTempOpVec() { return tempOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getVtOpVec() { return vtOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getFreqOpVec() { return freqOpVec_; }
  std::vector<Teuchos::RCP<astNode<usedType> > > & getGminOpVec() { return gminOpVec_; }

  std::vector<Teuchos::RCP<astNode<usedType> > > & getCurrentOpVec () { return currentOpVec_; };
  std::vector<Teuchos::RCP<astNode<usedType> > > & getUnresolvedCurrentOpVec() { return unresolvedCurrentOpVec_; };
  const std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & getCurrentOpNames ()
  {
    if (!astArraysSetup_) { setupVariousAstArrays_ (); }
    return currentOpNames_;
  };

  std::vector<Teuchos::RCP<astNode<usedType> > > & getLeadCurrentOpVec () { return leadCurrentOpVec_; };
  const std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & getLeadCurrentOpNames ()
  {
    if (!astArraysSetup_) { setupVariousAstArrays_ (); }
    return leadCurrentOpNames_;
  }

  std::vector<Teuchos::RCP<astNode<usedType> > > & getBsrcCurrentOpVec () { return bsrcCurrentOpVec_; };
  const std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & getBsrcCurrentOpNames ()
  {
    if (!astArraysSetup_) { setupVariousAstArrays_ (); }
    return bsrcCurrentOpNames_;
  }

  void codeGen( std::ostream & os ) { astNodePtr_->codeGen(os); os << ";" <<std::endl; };

  std::vector< Teuchos::RCP<astNode<usedType> > * > & getMasterNodeVec() { return masterAstNodeVec_; }
  std::vector< Teuchos::RCP<astNode<usedType> > > & getSrcNodeVec() { return srcAstNodeVec_;}
  std::vector< Teuchos::RCP<astNode<usedType> > > & getStpNodeVec() { return stpAstNodeVec_;}
  std::vector< Teuchos::RCP<astNode<usedType> > > & getCompNodeVec() { return compAstNodeVec_;}

  const std::string & getExpressionString() { return expressionString_; };
  const std::string & getOriginalExpressionString() { return originalExpressionString_; };

  bool replaceName ( const std::string & old_name, const std::string & new_name);
  bool replaceParamName ( const std::string & old_name, const std::string & new_name);

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

  // note: I don't particularly like these next 2 functions, but they are needed
  // to support the old expression API.
  //
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
  // Note; this was needed by the OLD API.  May not be needed now.
  void getFuncPrototypeName ( std::string & prototypeName) 
  {
    if(!(funcOpVec_.empty()))
    {
      prototypeName = funcOpVec_[0]->getName();
    }
  }

  void outputVariousAstArrays ( std::ostream & os );


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

  bool getIsConstant() { return isConstant_; }

  void treatAsTempAndConvert();

  bool setTemperature (const double & temp);

  void processSuccessfulTimeStep ();

  int getNumDdt () { return ddtOpVec_.size(); }
  void getDdtVals   (std::vector<double> & vals);
  void getDdtVals   (std::vector<std::complex<double> > & vals);
  void setDdtDerivs (std::vector<double> & vals);
  void setDdtDerivs (std::vector<std::complex<double> > & vals);

private:
  void setupDerivatives_ ();
  void setupVariousAstArrays_ ();
  void getValuesFromGroup_();

  Teuchos::RCP<baseExpressionGroup> group_;
  std::string expressionString_;
  std::string originalExpressionString_; // before toUpper, used for error messages
  bool parsed_;
  bool derivsSetup_;
  bool astArraysSetup_;

  Teuchos::RCP<astNode<usedType> > astNodePtr_;
  Teuchos::RCP<astNode<usedType> > * astNodePtrPtr_;
  Teuchos::RCP<tableOp<usedType> > * tableNodePtrPtr_;

  // function argument vectors.  Both strings and operators
  std::vector< std::string > functionArgStringVec_;

  Teuchos::RCP<paramOp<usedType> > garbageParamOpPtr_; // convenient for unused params/args

  std::vector< Teuchos::RCP<astNode<usedType> > > functionArgOpVec_;

  std::vector<std::string> paramNameVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > paramOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedParamOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > paramOpNames_;

  std::vector<std::string> funcNameVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > funcOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedFuncOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > funcOpMap_; 

  std::vector<Teuchos::RCP<astNode<usedType> > > voltOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedVoltOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > voltOpNames_;

  std::vector<Teuchos::RCP<astNode<usedType> > > currentOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedCurrentOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > currentOpNames_;

  std::vector<Teuchos::RCP<astNode<usedType> > > leadCurrentOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedLeadCurrentOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > leadCurrentOpNames_;

  std::vector<Teuchos::RCP<astNode<usedType> > > bsrcCurrentOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedBsrcCurrentOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > bsrcCurrentOpNames_;

  std::vector<Teuchos::RCP<astNode<usedType> > > powerOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedPowerOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > powerOpNames_;

  std::vector<Teuchos::RCP<astNode<usedType> > > internalDevVarOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedInternalDevVarOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > internalDevVarOpNames_;

  std::vector<Teuchos::RCP<astNode<usedType> > > dnoNoiseDevVarOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedDnoNoiseDevVarOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > dnoNoiseDevVarOpNames_;

  std::vector<Teuchos::RCP<astNode<usedType> > > dniNoiseDevVarOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedDniNoiseDevVarOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > dniNoiseDevVarOpNames_;

  std::vector<Teuchos::RCP<astNode<usedType> > > oNoiseOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedONoiseOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > oNoiseOpNames_;

  std::vector<Teuchos::RCP<astNode<usedType> > > iNoiseOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedINoiseOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > iNoiseOpNames_;

  std::vector<Teuchos::RCP<astNode<usedType> > > sdtOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedSdtOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > ddtOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedDdtOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > phaseOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedPhaseOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > sparamOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedSparamOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > yparamOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedYparamOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > zparamOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedZparamOpVec_;

  // master vector of nodes.  This is only used for deleting the ast tree in
  // the destructor.  The tree should be deleted by marching down the
  // branches of the tree, as some of the nodes use the same pointer.
  std::vector< Teuchos::RCP<astNode<usedType> > * > masterAstNodeVec_;

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

  // vector of STP objects.  Also needed for breakpoints.  STP devices will need to auto-detect breakpoints, if argument includes time
  std::vector< Teuchos::RCP<astNode<usedType> > > stpAstNodeVec_;

  // vector of COMP objects.  Also needed for breakpoints.  STP devices may need to auto-detect breakpoints, if comparing to a time
  std::vector< Teuchos::RCP<astNode<usedType> > > compAstNodeVec_;

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

  bool overrideGroupTemperature_;

  bool isConstant_;
  bool evaluateFunctionCalledBefore_;
  bool evaluateCalledBefore_;
  usedType savedResult_;
  bool phaseOutputUsesRadians_;

  opVectorContainers<usedType> opVectors_;

};

}
}

#endif

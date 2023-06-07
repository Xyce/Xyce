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
// Purpose        : ast_visitor
//
// Special Notes  : this header file contains functions and classes related to
//                  implementation of the "visitor pattern" for the
//                  expression library AST.
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 03/06/2023
//
//-----------------------------------------------------------------------------

#ifndef ast_visitor_H
#define ast_visitor_H

#include <checkGroundName.h>

// fwd declarations
template <typename ScalarT> struct opVectorContainers;

template <typename ScalarT> class astNode;
template <typename ScalarT> class numval;
template <typename ScalarT> class powOp;
template <typename ScalarT> class atan2Op;
template <typename ScalarT> class phaseOp;
template <typename ScalarT> class realOp;
template <typename ScalarT> class imagOp;
template <typename ScalarT> class maxOp;
template <typename ScalarT> class minOp;
template <typename ScalarT> class unaryNotOp;
template <typename ScalarT> class unaryMinusOp;
template <typename ScalarT> class unaryPlusOp;
template <typename ScalarT> class globalParamLayerOp;
template <typename ScalarT> class paramOp;
template <typename ScalarT> class voltageOp;
template <typename ScalarT> class currentOp;
template <typename ScalarT> class sparamOp;
template <typename ScalarT> class yparamOp;
template <typename ScalarT> class zparamOp;
template <typename ScalarT> class leadCurrentOp;
template <typename ScalarT> class powerOp;
template <typename ScalarT> class internalDevVarOp;
template <typename ScalarT> class dnoNoiseVarOp;
template <typename ScalarT> class dniNoiseVarOp;
template <typename ScalarT> class oNoiseOp;
template <typename ScalarT> class iNoiseOp;
template <typename ScalarT> class funcOp;
template <typename ScalarT> class pwrsOp;
template <typename ScalarT> class sgnOp;
template <typename ScalarT> class signOp;
template <typename ScalarT> class fmodOp;
template <typename ScalarT> class roundOp;
template <typename ScalarT> class ceilOp;
template <typename ScalarT> class floorOp;
template <typename ScalarT> class intOp;
template <typename ScalarT> class ifStatementOp;
template <typename ScalarT> class limitOp;
template <typename ScalarT> class stpOp;
template <typename ScalarT> class urampOp;
template <typename ScalarT> class tableOp;
template <typename ScalarT> class scheduleOp;
template <typename ScalarT> class sdtOp;
template <typename ScalarT> class ddtOp;
template <typename ScalarT> class ddxOp;
template <typename ScalarT> class specialsOp;
template <typename ScalarT> class piConstOp;
template <typename ScalarT> class CtoKConstOp;
template <typename ScalarT> class agaussOp;
template <typename ScalarT> class gaussOp;
template <typename ScalarT> class aunifOp;
template <typename ScalarT> class unifOp;
template <typename ScalarT> class randOp;
template <typename ScalarT> class twoArgLimitOp;
template <typename ScalarT> class spicePulseOp;
template <typename ScalarT> class spiceSinOp;
template <typename ScalarT> class spiceExpOp;
template <typename ScalarT> class spiceSffmOp;
template <typename ScalarT> class binaryAddOp;
template <typename ScalarT> class binaryMinusOp;
template <typename ScalarT> class binaryMulOp;
template <typename ScalarT> class binaryDivOp;
template <typename ScalarT> class gtOp;
template <typename ScalarT> class ltOp;
template <typename ScalarT> class neOp;
template <typename ScalarT> class eqOp;
template <typename ScalarT> class geOp;
template <typename ScalarT> class leOp;
template <typename ScalarT> class orOp;
template <typename ScalarT> class andOp;
template <typename ScalarT> class xorOp;
template <typename ScalarT> class sqrtOp;
template <typename ScalarT> class expOp;
template <typename ScalarT> class absOp;
template <typename ScalarT> class sinOp;
template <typename ScalarT> class cosOp;
template <typename ScalarT> class acosOp;
template <typename ScalarT> class acoshOp;
template <typename ScalarT> class asinOp;
template <typename ScalarT> class asinhOp;
template <typename ScalarT> class coshOp;
template <typename ScalarT> class logOp;
template <typename ScalarT> class log10Op;
template <typename ScalarT> class sinhOp;
template <typename ScalarT> class tanOp;
template <typename ScalarT> class atanOp;
template <typename ScalarT> class tanhOp;
template <typename ScalarT> class atanhOp;


//-------------------------------------------------------------------------------
// base visitor class
//-------------------------------------------------------------------------------
template <typename ScalarT>
class nodeVisitor
{
  public:
  virtual void visit(Teuchos::RCP<numval<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<powOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<atan2Op<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<phaseOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<realOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<imagOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<maxOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<minOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<unaryNotOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<unaryMinusOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<unaryPlusOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<globalParamLayerOp<ScalarT> > & astNode) {}

  virtual void visit(Teuchos::RCP<paramOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<voltageOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<currentOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<sparamOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<yparamOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<zparamOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<leadCurrentOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<powerOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<internalDevVarOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<dnoNoiseVarOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<dniNoiseVarOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<oNoiseOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<iNoiseOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<funcOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<pwrsOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<sgnOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<signOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<fmodOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<roundOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<ceilOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<floorOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<intOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<ifStatementOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<limitOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<stpOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<urampOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<tableOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<scheduleOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<sdtOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<ddtOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<ddxOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<specialsOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<piConstOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<CtoKConstOp<ScalarT> > & astNode) {}

  // random operators
  virtual void visit(Teuchos::RCP<agaussOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<gaussOp<ScalarT> > & astNode)  {}
  virtual void visit(Teuchos::RCP<aunifOp<ScalarT> > & astNode)  {}
  virtual void visit(Teuchos::RCP<unifOp<ScalarT> > & astNode)   {}
  virtual void visit(Teuchos::RCP<randOp<ScalarT> > & astNode)   {}
  virtual void visit(Teuchos::RCP<twoArgLimitOp<ScalarT> > & astNode) {}

  // spice pulse operators
  virtual void visit(Teuchos::RCP<spicePulseOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<spiceSinOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<spiceExpOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<spiceSffmOp<ScalarT> > & astNode) {}

  // binary operators
  virtual void visit(Teuchos::RCP<binaryAddOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<binaryMinusOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<binaryMulOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<binaryDivOp<ScalarT> > & astNode) {}

  // comparison operators
  virtual void visit(Teuchos::RCP<gtOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<ltOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<neOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<eqOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<geOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<leOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<orOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<andOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<xorOp<ScalarT> > & astNode) {}

  // math functions
  virtual void visit(Teuchos::RCP<sqrtOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<expOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<absOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<sinOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<cosOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<acosOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<acoshOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<asinOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<asinhOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<coshOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<logOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<log10Op<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<sinhOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<tanOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<atanOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<tanhOp<ScalarT> > & astNode) {}
  virtual void visit(Teuchos::RCP<atanhOp<ScalarT> > & astNode) {}
};

// each derived visitor class is for a different operation to be done as we "visit" all the nodes.

//-------------------------------------------------------------------------------
// getInterestingOps visitor.  
// This is the most important visitor.
//-------------------------------------------------------------------------------
template <typename ScalarT>
class getInterestingOpsVisitor : public nodeVisitor<ScalarT>
{
  public:
  getInterestingOpsVisitor (opVectorContainers<ScalarT> & tmp) : ovc(tmp) {}
  opVectorContainers<ScalarT> & ovc;

  void visit(Teuchos::RCP<paramOp<ScalarT> > & astNode)       
  { 
    if ( !(astNode->getFunctionArgType()) )  // parameters are occasionally function arguments.  Don't include those
    { ovc.paramOpVector.push_back(astNode); }
  }
  void visit(Teuchos::RCP<voltageOp<ScalarT> > & astNode)
  { // don't bother if this is ground
   if ( !Xyce::Util::checkGroundNodeName( astNode->getVoltageNode()) ) { ovc.voltOpVector.push_back(astNode); }
  }
  void visit(Teuchos::RCP<currentOp<ScalarT> > & astNode)     { ovc.currentOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<sparamOp<ScalarT> > & astNode)      { ovc.sparamOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<yparamOp<ScalarT> > & astNode)      { ovc.yparamOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<zparamOp<ScalarT> > & astNode)      { ovc.zparamOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<leadCurrentOp<ScalarT> > & astNode) { ovc.leadCurrentOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<powerOp<ScalarT> > & astNode)       { ovc.powerOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<internalDevVarOp<ScalarT> > & astNode) { ovc.internalDevVarOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<dnoNoiseVarOp<ScalarT> > & astNode)    { ovc.dnoNoiseDevVarOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<dniNoiseVarOp<ScalarT> > & astNode)    { ovc.dniNoiseDevVarOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<oNoiseOp<ScalarT> > & astNode)         { ovc.oNoiseOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<iNoiseOp<ScalarT> > & astNode)         { ovc.iNoiseOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<funcOp<ScalarT> > & astNode) { ovc.funcOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<limitOp<ScalarT> > & astNode) { ovc.limitOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<stpOp<ScalarT> > & astNode) { ovc.stpOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<tableOp<ScalarT> > & astNode) { if (astNode->srcType()) { ovc.srcOpVector.push_back(astNode); } }
  void visit(Teuchos::RCP<sdtOp<ScalarT> > & astNode) { ovc.sdtOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<ddtOp<ScalarT> > & astNode) { ovc.ddtOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<agaussOp<ScalarT> > & astNode) { ovc.agaussOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<gaussOp<ScalarT> > & astNode)  { ovc.gaussOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<aunifOp<ScalarT> > & astNode)  { ovc.aunifOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<unifOp<ScalarT> > & astNode)   { ovc.unifOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<randOp<ScalarT> > & astNode)   { ovc.randOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<twoArgLimitOp<ScalarT> > & astNode) { ovc.twoArgLimitOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<spicePulseOp<ScalarT> > & astNode) { ovc.srcOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<spiceSinOp<ScalarT> > & astNode) { ovc.srcOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<spiceExpOp<ScalarT> > & astNode) { ovc.srcOpVector.push_back(astNode); }
  void visit(Teuchos::RCP<spiceSffmOp<ScalarT> > & astNode) { ovc.srcOpVector.push_back(astNode); }

  virtual void visit(Teuchos::RCP<scheduleOp<ScalarT> > & astNode) { ovc.isScheduleDependent = true; }

  // comparison operators
  virtual void visit(Teuchos::RCP<gtOp<ScalarT> > & astNode) { ovc.compOpVector.push_back(astNode); }
  virtual void visit(Teuchos::RCP<ltOp<ScalarT> > & astNode) { ovc.compOpVector.push_back(astNode); }
  virtual void visit(Teuchos::RCP<neOp<ScalarT> > & astNode) { ovc.compOpVector.push_back(astNode); }
  virtual void visit(Teuchos::RCP<eqOp<ScalarT> > & astNode) { ovc.compOpVector.push_back(astNode); }
  virtual void visit(Teuchos::RCP<geOp<ScalarT> > & astNode) { ovc.compOpVector.push_back(astNode); }
  virtual void visit(Teuchos::RCP<leOp<ScalarT> > & astNode) { ovc.compOpVector.push_back(astNode); }
  virtual void visit(Teuchos::RCP<orOp<ScalarT> > & astNode) { ovc.compOpVector.push_back(astNode); }
  virtual void visit(Teuchos::RCP<andOp<ScalarT> > & astNode) { ovc.compOpVector.push_back(astNode); }
  virtual void visit(Teuchos::RCP<xorOp<ScalarT> > & astNode) { ovc.compOpVector.push_back(astNode); }

  virtual void visit(Teuchos::RCP<specialsOp<ScalarT> > & astNode) 
  { 
    if (astNode->timeSpecialType() || astNode->dtSpecialType()) { ovc.isTimeDependent = true; } 
    else if (astNode->tempSpecialType()) { ovc.isTempDependent = true; } 
    else if (astNode->vtSpecialType()) { ovc.isVTDependent = true; } 
    else if (astNode->freqSpecialType()) { ovc.isFreqDependent = true; } 
    else if (astNode->gminSpecialType()) { ovc.isGminDependent = true; } 
  }
};

template <typename ScalarT>
class getParamOpsVisitor : public nodeVisitor<ScalarT>
{
  public:
  getParamOpsVisitor (std::vector<Teuchos::RCP<astNode<ScalarT> > > & tmp) : parVec(tmp) {}
  std::vector<Teuchos::RCP<astNode<ScalarT> > > & parVec;
  // while paramOps *can* be either a function argument or not, for some use cases you need both.
  // In the ddxOp class, you need both.  As of this writing *only* the ddxOp class uses this visitor.
  // So this visitor does *not* check for functionArgType.
  void visit(Teuchos::RCP<paramOp<ScalarT> > & astNode) { parVec.push_back(astNode); }
};

template <typename ScalarT>
class getVoltageOpsVisitor : public nodeVisitor<ScalarT>
{
  public:
  getVoltageOpsVisitor (std::vector<Teuchos::RCP<astNode<ScalarT> > > & tmp) : voltVec(tmp) {}
  std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltVec;
  void visit(Teuchos::RCP<voltageOp<ScalarT> > & astNode)
  {  // don't bother if this is ground
    if ( !Xyce::Util::checkGroundNodeName( astNode->getVoltageNode()) ) { voltVec.push_back(astNode); }
  }
};

template <typename ScalarT>
class getCurrentOpsVisitor : public nodeVisitor<ScalarT>
{
  public:
  getCurrentOpsVisitor (std::vector<Teuchos::RCP<astNode<ScalarT> > > & tmp) :currVec(tmp) {}
  std::vector<Teuchos::RCP<astNode<ScalarT> > > & currVec;
  void visit(Teuchos::RCP<currentOp<ScalarT> > & astNode) { currVec.push_back(astNode); }
};

template <typename ScalarT>
class getTimeOpsVisitor : public nodeVisitor<ScalarT>
{
  public:
  getTimeOpsVisitor (std::vector<Teuchos::RCP<astNode<ScalarT> > > & tmp) : tov(tmp) {}
  std::vector<Teuchos::RCP<astNode<ScalarT> > > & tov;

  void visit(Teuchos::RCP<specialsOp<ScalarT> > & astNode) 
  { 
    if (astNode->timeSpecialType() || astNode->dtSpecialType()) { tov.push_back(astNode); }
  }
};

//-------------------------------------------------------------------------------
// this object is used by "getInterestingOps".  It is a catch-all, that contains 
// references to all the book-keeping STL containers needed by the expression library
//-------------------------------------------------------------------------------
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
  bool gminDep,
  bool scheduleDep
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
    isGminDependent(gminDep),
    isScheduleDependent(scheduleDep)
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
  bool isScheduleDependent;

  // Visitor functions.  
  void getInterestingOps(Teuchos::RCP<astNode<ScalarT> > & ast)
  {
    getInterestingOpsVisitor<ScalarT> visitor(*this);
    ast->accept(visitor,ast); // << 1st dispatch  ("accept" the visitor,ast)
  }
};

#endif

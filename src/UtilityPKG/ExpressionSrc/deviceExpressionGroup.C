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
// Creation Date  : 8/20/2020
//
//
//
//
//-----------------------------------------------------------------------------

#include <iostream>
#include <unordered_map>
#include <string>
#include <random>

#include <deviceExpressionGroup.h>
#include <ast.h>
#include <newExpression.h>

#include <N_TIA_DataStore.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TOP_Topology.h>
#include <N_LAS_Vector.h>
#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Manager.h>
#include <N_PDS_Serial.h>

#include <N_ANP_AnalysisManager.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_Op.h>
#include <N_DEV_Const.h>

#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Util {

//-------------------------------------------------------------------------------
// Function      : deviceExpressionGroup::deviceExpressionGroup 
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 8/20/2020
//-------------------------------------------------------------------------------
deviceExpressionGroup::deviceExpressionGroup ( 
   const Teuchos::RCP<Xyce::Util::mainXyceExpressionGroup> & mainGroup
    ):
  mainXyceExpressionGroup(
      mainGroup->comm_, mainGroup->top_,
      mainGroup->analysisManager_,
      mainGroup->deviceManager_,
      mainGroup->outputManager_
      )
{
}

//-------------------------------------------------------------------------------
// Function      : deviceExpressionGroup::deviceExpressionGroup 
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 8/20/2020
//-------------------------------------------------------------------------------
deviceExpressionGroup::~deviceExpressionGroup ()
{
}

//-------------------------------------------------------------------------------
// Function      : deviceExpressionGroup::setSolutionLIDs 
// Purpose       : Order of expVarNames:  nodes, leads, instances
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 8/22/2020
//-------------------------------------------------------------------------------
void deviceExpressionGroup::setSolutionLIDs( 
    const std::vector<std::string> & expVarNames, 
    const std::vector<int> & expVarLIDs, int lo, int hi)
{
  for (int ii=lo;ii<hi;ii++)
  {
    lidMap_[expVarNames[ii]] = expVarLIDs[ii-lo];
  }
}

//-------------------------------------------------------------------------------
// Function      : deviceExpressionGroup::setupGroup
// Purpose       : This group mostly sets up the order.
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 5/23/2021
//-------------------------------------------------------------------------------
bool deviceExpressionGroup::setupGroup(newExpression &expr)
{
  lidVec_.clear();
  bool success = true;

  if (lidMap_.empty() && !(expr.voltOpVec_.empty())) { success = false; }
  if (lidMap_.empty() && !(expr.currentOpVec_.empty())) { success = false; }

  if (!(lidMap_.empty()))
  {
    if ( !(expr.voltOpVec_.empty()) )
    {
      for (int ii=0;ii<expr.voltOpVec_.size();ii++)
      {
        Teuchos::RCP<voltageOp<usedType> > voltOp
          = Teuchos::rcp_static_cast<voltageOp<usedType> > (expr.voltOpVec_[ii]);
        const std::string & node = voltOp->getVoltageNode();

        if (lidMap_.find(node) != lidMap_.end()) { lidVec_.push_back( lidMap_[node] ); }
        else { success = false; }
      }
    }

    if ( !(expr.currentOpVec_.empty()) )
    {
      for (int ii=0;ii<expr.currentOpVec_.size();ii++)
      {
        Teuchos::RCP<currentOp<usedType> > currOp 
          = Teuchos::rcp_static_cast<currentOp<usedType> > (expr.currentOpVec_[ii]);
        const std::string & name = currOp->getCurrentDevice();
        if (lidMap_.find(name) != lidMap_.end()) { lidVec_.push_back( lidMap_[name] ); }
        else { success = false; }
      }
    }
  }

  return success;
}

//-------------------------------------------------------------------------------
// Function      : deviceExpressionGroup::putValues
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 5/23/2021
//-------------------------------------------------------------------------------
bool deviceExpressionGroup::putValues(newExpression & expr)
{
  bool noChange=true;
  int index=0;
  const Linear::Vector * nextSolVector = deviceManager_.getExternData().nextSolVectorPtr;

  if ( !(expr.voltOpVec_.empty()) )
  {
    for (int ii=0;ii<expr.voltOpVec_.size();ii++)
    {
      Teuchos::RCP<voltageOp<usedType> > voltOp
        = Teuchos::rcp_static_cast<voltageOp<usedType> > (expr.voltOpVec_[ii]);
      usedType & val = voltOp->getVoltageVal();
      usedType oldval = val;
      if (nextSolVector && index < lidVec_.size()) { val = (*nextSolVector)[lidVec_[index]]; index++; }
      if(val != oldval) noChange=false;
    }
  }

  if ( !(expr.currentOpVec_.empty()) )
  {
    for (int ii=0;ii<expr.currentOpVec_.size();ii++)
    {
      Teuchos::RCP<currentOp<usedType> > currOp 
        = Teuchos::rcp_static_cast<currentOp<usedType> > (expr.currentOpVec_[ii]);
      usedType & val = currOp->getCurrentVal();
      usedType oldval = val;
      if (nextSolVector && index < lidVec_.size()) { val = (*nextSolVector)[lidVec_[index]]; index++; }
      if (val != oldval) noChange=false;
    }
  }

  return noChange;
}

//-------------------------------------------------------------------------------
// Function      : deviceExpressionGroup::getSolutionVal
// Purpose       : 
// Special Notes : double precision version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 8/20/2020
//-------------------------------------------------------------------------------
bool deviceExpressionGroup::getSolutionVal(const std::string & nodeName, double & retval )
{
  retval = 0.0;
  if (lidMap_.find(nodeName) != lidMap_.end())
  {
    int lid = lidMap_[nodeName];
    const Linear::Vector * nextSolVector = deviceManager_.getExternData().nextSolVectorPtr;
    if (nextSolVector) { retval = (*nextSolVector)[lid]; }
    return true;
  }
  else { return false; }
}

//-------------------------------------------------------------------------------
// Function      : deviceExpressionGroup::getSolutionVal
// Purpose       : 
// Special Notes : std::complex<double> precision version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 8/20/2020
//-------------------------------------------------------------------------------
bool deviceExpressionGroup::getSolutionVal(const std::string & nodeName, std::complex<double> & retval)
{
  double real_val=0.0;
  double imag_val=0.0;
  if (lidMap_.find(nodeName) != lidMap_.end())
  {
    int lid = lidMap_[nodeName];
    const Linear::Vector * nextSolVector = deviceManager_.getExternData().nextSolVectorPtr;
    if (nextSolVector) { real_val = (*nextSolVector)[lid]; }
    retval = std::complex<double>(real_val,imag_val);
    return true;
  }
  else { return false; }
}

} // Util
} // Xyce

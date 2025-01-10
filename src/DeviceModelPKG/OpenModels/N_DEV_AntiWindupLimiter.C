//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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

//----------------------------------------------------------------------------
//
// Purpose        : Anti-Windup Limiter: See Section C.2 of "Power System
//                  Modeling and Scripting" by F. Milano
//
// Special Notes  : Experimental new device for an LDRD.
//
// Creator        : Pete Sholander, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 11/12/15
//
//
//
//
//----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_AntiWindupLimiter.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Message.h>
#include <N_DEV_SolverState.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_IndentStreamBuf.h>

namespace Xyce {
namespace Device {
namespace AntiWindupLimiter {

// enables debug output for just this device
//static const int DEBUG_DEVICE = 1;

void Traits::loadInstanceParameters(ParametricData<AntiWindupLimiter::Instance> &p)
{
  p.addPar("T", 1.0, &AntiWindupLimiter::Instance::T_)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_SECOND)
    .setDescription("Time Constant");

  p.addPar("UL", 1.0, &AntiWindupLimiter::Instance::upperLimit_)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_PERUNIT)
    .setDescription("Upper Limit");

  p.addPar("LL", -1.0, &AntiWindupLimiter::Instance::lowerLimit_)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_PERUNIT)
    .setDescription("Lower Limit");
}

void Traits::loadModelParameters(ParametricData<AntiWindupLimiter::Model> &p)
{}


std::vector< std::vector<int> > Instance::jacStamp;

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/12/15
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &               IB,
  Model & Miter,
  const FactoryBlock &  factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    T_(1.0),
    lowerLimit_(1.0),
    upperLimit_(1.0),
    holdAtUpper_(false),
    holdAtLower_(false),
    upperLimitGiven_(false),
    lowerLimitGiven_(false)
{
  // used by the Topology Checker.  See SON Bug 974 for more details.
  devConMap.resize(2);
  devConMap[0]=1;
  devConMap[1]=1;

  numIntVars   = 1;
  numExtVars   = 2;  
  numStateVars = 0; 

  if( jacStamp.empty() )
  {
    jacStamp.resize(3);
    jacStamp[0].resize(1);
    jacStamp[1].resize(1);
    jacStamp[2].resize(2);

    jacStamp[0][0] = 2;  
    jacStamp[1][0] = 2;
    jacStamp[2][0] = 0;
    jacStamp[2][1] = 1;
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/10/15
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/08/16
//-----------------------------------------------------------------------------
bool Instance::processParams()
{
  bool bsuccess=true;

  if (given("UL"))
  {
    upperLimitGiven_ = true;
  }
  if (given("LL"))
  {
    lowerLimitGiven_ = true;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << getName() << " (T,UL,LL) are: "
		 << "(" << T_ << " , " << upperLimit_ << " , " << lowerLimit_ << ")" 
                 << std::endl << std::endl;
  }

  // Sanity check on input parameters.
  // Upper limit must be greater than the lower limit. 
  // T must be strictly positive.
  if (lowerLimit_ >= upperLimit_)
  {
    UserError(*this) << "Upper limit (UL) must be greater than lower limit (LL)";
    bsuccess=false;
  }
  if (T_ <= 0)
  {
    UserError(*this) << "Time constant (T) must be positive";
    bsuccess=false;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/12/15
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef )
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  li_vin = extLIDVec[0];
  li_vout = extLIDVec[1];
  li_BranCurr = intLIDVec[0];
    
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_LIDS))
  {
    Xyce::dout() << getName() << " LIDs"
        << Util::push << std::endl
           << "li_vin = " << li_vin << std::endl
           << "li_vout = " << li_vout << std::endl
           << "li_BranCurr = " << li_BranCurr << std::endl
           << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadNodeSymbols
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
void Instance::loadNodeSymbols(Util::SymbolTable &symbol_table) const
{
  addInternalNode(symbol_table, li_BranCurr, getName(), "BranchCurr");
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       : Set up offsets so we can store quantities that need to be
//                 differentiated.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/12/15
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs(const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/20/01
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/27/01
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  VI_I_Offset = jacLIDVec[0][0];
  VO_I_Offset = jacLIDVec[1][0];
  I_VI_Offset = jacLIDVec[2][0];
  I_VO_Offset = jacLIDVec[2][1];

  if (DEBUG_DEVICE && (isActive(Diag::DEVICE_JACSTAMP) || isActive(Diag::DEVICE_LIDS)))
  {
    Xyce::dout() << getName() << ": In registerJacLIDs, the offsets are: (" 
                 << VI_I_Offset << " , " << VO_I_Offset << " ," 
                 << I_VI_Offset << " , " << I_VO_Offset << ")" 
                 << std::endl << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one anti-windup limiter instance
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/12/15
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  double * solVec = extData.nextSolVectorRawPtr;

  // add limits on output here
  VIN_ = solVec[li_vin];
  VOUT_ = solVec[li_vout];
  BranCurr_ = solVec[li_BranCurr];

  // QV_ is used as a proxy for the derivative of VOUT in the
  // absence of limiting.
  QV_ = VIN_ - VOUT_;

  // default is no limiting, and then test for limit conditions
  holdAtUpper_ = false;
  holdAtLower_ = false;
  if ( (VOUT_ >= upperLimit_) && upperLimitGiven_ && (QV_ >= 0) )
  {
    holdAtUpper_ = true;
  }
  else if  ( (VOUT_ <= lowerLimit_) && lowerLimitGiven_ && (QV_ <= 0))
  {
     holdAtLower_ = true;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/12/15
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  return updateIntermediateVars ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions 
//
// Special Notes :
//
//
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/12/15
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;
  
  // ensure that current into the limiter is equal to the current
  // out of the limitier.
  fVec[li_vin] += BranCurr_;
  fVec[li_vout]-= BranCurr_;

  // VOUT acts like a voltage source, with a value of UL or LL when
  // limiting is being enforced. In that case, it will have a B-vector
  // component. Otherwise, the device acts like a 1st order low-pass
  // filter, with time-constant T and gain K.  In the latter
  // case, it does not have a B-vector component.
  if (holdAtUpper_ || holdAtLower_)
  {
    fVec[li_BranCurr] += VOUT_;
  }
  else
  { 
    fVec[li_BranCurr] += (VOUT_ - VIN_)/T_;
  }

  if (getSolverState().dcopFlag)
  {
    // no op placeholder
  }
  else
  {
    // no op placeholder
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions.  
//
// Special Notes :
//
//
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/12/15
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  double * qVec = extData.daeQVectorRawPtr;

  // load contribution for derivative of VOUT in all cases.
  qVec[li_BranCurr] += VOUT_;
   
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEBVector
//
// Purpose       : Loads the B-vector contributions for a single
//                 anti-windup limiter instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool Instance::loadDAEBVector ()
{
  double * bVec = extData.daeBVectorRawPtr;

  // VOUT acts like a voltage source, with a value of UL or LL when
  // limiting is being enforced. In that case, it will have a B-vector
  // component. Otherwise, the device acts like a 1st order low-pass
  // filter, with time-constant T and gain K.  In the latter
  // case, it does not have a B-vector component.
  if (holdAtUpper_ )
  {
    bVec[li_BranCurr] += upperLimit_;
  }
  else if (holdAtLower_)
  {
    bVec[li_BranCurr] += lowerLimit_;
  }        

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the dFdx-matrix contributions 
//
// Special Notes :
//
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/12/15
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_vin][VI_I_Offset] += 1.0;
  dFdx[li_vout][VO_I_Offset] -= 1.0;

  if (holdAtUpper_ || holdAtLower_)
  {
    dFdx[li_BranCurr][I_VO_Offset] += 1.0;
  }
  else
  {
    dFdx[li_BranCurr][I_VI_Offset] -= 1.0/T_;
    dFdx[li_BranCurr][I_VO_Offset] += 1.0/T_;
  }

  if (getSolverState().dcopFlag)
  {
    // no op placeholder code;
  }
  else
  {
    // no op placeholder code;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx ()
//
// Purpose       : Loads the  dQdx-matrix contributions.  
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);  

  dQdx[li_BranCurr][I_VO_Offset] += 1;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : "Model block" constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------


// These are all just placeholder functions, as there is nothing to the model
// class.
// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : "Model block" constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block)
{
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
Model::~Model ()
{
  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    delete (*iter);
  }
}

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/12/15
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << std::endl;
  os << "    name     model name  Parameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "      ";
    os << getName();
    os << std::endl;
  }

  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Model::forEachInstance
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date : 2/4/2014
//-----------------------------------------------------------------------------
/// Apply a device instance "op" to all instances associated with this
/// model
/// 
/// @param[in] op Operator to apply to all instances.
/// 
/// 
void Model::forEachInstance(DeviceInstanceOp &op) const /* override */ 
{
  for (std::vector<Instance *>::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
    op(*it);
}


Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new DeviceMaster<Traits>(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() || (deviceMap.find("ANTIWINDUPLIMITER")!=deviceMap.end())
      || (deviceMap.find("AWL")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("AntiWindupLimiter", 1)
      .registerDevice("AWL", 1)
      .registerModelType("AWL", 1);
  }
}

} // namespace AntiWindupLimiter
} // namespace Device
} // namespace Xyce

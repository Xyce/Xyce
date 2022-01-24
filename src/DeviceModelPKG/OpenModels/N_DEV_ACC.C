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

//----------------------------------------------------------------------------
//
// Purpose        : Provide a solution to the simple second order
//                  initial value problem:
//                  d^2x/dt^2 = a(t); x(0) = x0, dx/dt(0) = v0
//
// Special Notes  : Intended for use in FY07 CoilGun LDRD work
//
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 10/25/07
//
//
//
//
//----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_ACC.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Message.h>
#include <N_DEV_SolverState.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_IndentStreamBuf.h>

namespace Xyce {
namespace Device {
namespace ACC {

void Traits::loadInstanceParameters(ParametricData<ACC::Instance> &p)
{
  p.addPar("V0", 0.0, &ACC::Instance::v0)
    .setUnit(U_MSM1)
    .setDescription("Initial Velocity");

  p.addPar("X0", 0.0, &ACC::Instance::x0)
    .setUnit(U_METER)
    .setDescription("Initial Position");
}

void Traits::loadModelParameters(ParametricData<ACC::Model> &p)
{}



std::vector< std::vector<int> > Instance::jacStamp;

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &               IB,
  Model & Miter,
  const FactoryBlock &  factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    v0(0.0),
    x0(0.0),
    li_Acc(-1),
    li_Velocity(-1),
    li_Position(-1),
    li_state_vel(-1),
    li_state_pos(-1),
    AVelEquAccNodeOffset(-1),
    AVelEquVelNodeOffset(-1),
    APosEquVelNodeOffset(-1),
    APosEquPosNodeOffset(-1)
{
  numIntVars   = 0;
  numExtVars   = 3;
  numStateVars = 2; // position and velocity saved in state so we can get
                    //derivs

  devConMap.resize(3);
  devConMap[0]=1;
  devConMap[1]=1;
  devConMap[2]=1;

  if( jacStamp.empty() )
  {
    jacStamp.resize(3);
    //jacStamp[0].resize(0); // the row for acceleration has nothing.
    jacStamp[1].resize(2); // velocity row
    jacStamp[1][0]=0;      //velocity-acceleration
    jacStamp[1][1]=1;      //velocity-velocity
    jacStamp[2].resize(2);
    jacStamp[2][0]=1;      // position-velocity
    jacStamp[2][1]=2;      // position-position
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:

  setParams (IB.params);

  // Set any non-constant parameter defaults:


}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef )
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  li_Acc = extLIDVec[0];
  li_Velocity = extLIDVec[1];
  li_Position = extLIDVec[2];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    dout() << getName() << " LIDs"
      << Util::push << std::endl
           << "li_Acc = " << li_Acc << std::endl
           << "li_Velocity = " << li_Velocity << std::endl
           << "li_Position = " << li_Position << std::endl
           << Util::pop << std::endl;
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
{}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       : Set up offsets so we can store quantities that need to be
//                 differentiated.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs(const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // Copy over the global ID lists:
  staLIDVec = staLIDVecRef;

  li_state_vel = staLIDVec[0];
  li_state_pos = staLIDVec[1];

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

  AVelEquAccNodeOffset = jacLIDVec[1][0];
  AVelEquVelNodeOffset = jacLIDVec[1][1];
  APosEquVelNodeOffset = jacLIDVec[2][0];
  APosEquPosNodeOffset = jacLIDVec[2][1];
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one diode instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;

  double * solVec = extData.nextSolVectorRawPtr;

  position = solVec[li_Position];
  velocity = solVec[li_Velocity];
  acceleration = solVec[li_Acc];

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::updatePrimaryState" <<std::endl;
  }

  bool bsuccess = updateIntermediateVars ();
  double * staVec = extData.nextStaVectorRawPtr;

  // We need to save the position and velocity so we can differentiate them
  // to give what we *think* the velocity and acceleration are.
  staVec[li_state_pos] = position;
  staVec[li_state_vel] = velocity;

  // if this is the first newton step of the first time step
  // of the transient simulation, we need to enforce that the
  // time derivatives w.r.t. charge are zero.  This is to maintain 3f5
  // compatibility.  ERK.

  if (!(getSolverState().dcopFlag) && (getSolverState().initTranFlag_) && getSolverState().newtonIter==0)
  {
    double * currStaVec = extData.currStaVectorRawPtr;
    currStaVec[li_state_pos] = position;
    currStaVec[li_state_vel] = velocity;
  }

  return  bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  double * staDerivVec = extData.nextStaDerivVectorRawPtr;

  xdot = staDerivVec[li_state_pos];
  vdot = staDerivVec[li_state_vel];

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 ACC  instance.
//
// Special Notes :
//
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;

  // Load DAE F-vector
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "  Instance::loadDAEFVector" << std::endl;
    Xyce::dout() << "  name = " << getName() <<std::endl;
    Xyce::dout() << "  velocity = " << velocity << std::endl;
    Xyce::dout() << "  position = " << position << std::endl;
  }

  if (getSolverState().dcopFlag)
  {

    fVec[li_Velocity] += velocity-v0;

    fVec[li_Position] += position-x0;
  }
  else
  {

    fVec[li_Velocity] += -acceleration;

    fVec[li_Position] += -velocity;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 ACC  instance.
//
// Special Notes :
//
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  double * qVec = extData.daeQVectorRawPtr;


  qVec[li_Velocity] += velocity;

  qVec[li_Position] += position;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 ACC  instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  if (getSolverState().dcopFlag)
  {

    dFdx[li_Velocity][AVelEquVelNodeOffset] += 1.0;

    dFdx[li_Position][APosEquPosNodeOffset] += 1.0;
  }
  else
  {

    dFdx[li_Velocity][AVelEquAccNodeOffset] += -1.0;

    dFdx[li_Position][APosEquVelNodeOffset] += -1.0;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx ()
//
// Purpose       : Loads the Q-vector contributions for a single
//                 ACC  instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);

  if (!getSolverState().dcopFlag)
  {

    dQdx[li_Velocity][AVelEquVelNodeOffset] += 1.0;

    dQdx[li_Position][APosEquPosNodeOffset] += 1.0;
  }

  return true;
}

// These are all just placeholder functions, as there is nothing to the model
// class.
// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : "Model block" constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
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
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
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
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
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

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() || (deviceMap.find("ACC")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("acc", 1);
  }
}

} // namespace ACC
} // namespace Device
} // namespace Xyce

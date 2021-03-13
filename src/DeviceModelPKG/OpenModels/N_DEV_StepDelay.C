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
// Purpose        :  Provide a step-delay element for ANN work
//
// Special Notes  :  Voltage across first two nodes should be the voltage
//                   at control nodes delayed by one timestep
//
// Creator        : Tom Russo
// Creation Date  : 14 Apr 2020
//
//
//
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_StepDelay.h>
#include <N_DEV_Bsrc.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {

namespace StepDelay {

void Traits::loadInstanceParameters(ParametricData<StepDelay::Instance> &p)
{
}

void Traits::loadModelParameters(ParametricData<StepDelay::Model> &p)
{
}

std::vector< std::vector<int> > Instance::jacStamp;

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 14 Apr 2020
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       Viter,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Viter),
    li_Pos(-1),
    li_Neg(-1),
    li_Bra(-1),
    li_ContPos(-1),
    li_ContNeg(-1),
    li_branch_data(0),
    li_storeVold(-1),
    ABraEquPosNodeOffset(-1),
    ABraEquNegNodeOffset(-1),
    ABraEquContPosNodeOffset(-1),
    ABraEquContNegNodeOffset(-1),
    APosEquBraVarOffset(-1),
    ANegEquBraVarOffset(-1),

    f_BraEquPosNodePtr(0),
    f_BraEquNegNodePtr(0),
    f_BraEquContPosNodePtr(0),
    f_BraEquContNegNodePtr(0),
    f_PosEquBraVarPtr(0),
    f_NegEquBraVarPtr(0)
{
  numIntVars   = 1;
  numExtVars   = 4;
  numStateVars = 0;
  setNumBranchDataVars(0);          // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1; // this is the space to allocate if lead current or power is needed.

  setNumStoreVars(1);
  
  if( jacStamp.empty() )
  {
    jacStamp.resize(5);
    jacStamp[0].resize(1);
    jacStamp[0][0]=4;
    jacStamp[1].resize(1);
    jacStamp[1][0]=4;
    jacStamp[4].resize(4);
    jacStamp[4][0]=0;
    jacStamp[4][1]=1;
    jacStamp[4][2]=2;
    jacStamp[4][3]=3;
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
//-----------------------------------------------------------------------------
Instance::~Instance ()
{

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
//-----------------------------------------------------------------------------
void Instance::registerLIDs ( const std::vector<int> & intLIDVecRef,
                                        const std::vector<int> & extLIDVecRef )
{
  std::string msg;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  StepDelayInstance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

  if (numInt != numIntVars)
  {
    DevelFatal(*this).in("Instance::registerLIDs")
      <<  "numInt != numIntVars";
  }

  if (numExt != numExtVars)
  {
    DevelFatal(*this).in("Instance::registerLIDs")
      <<  "numExt != numExtVars";
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];
  li_ContPos = extLIDVec[2];
  li_ContNeg = extLIDVec[3];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl;
    Xyce::dout() << "  li_Neg = " << li_Neg << std::endl;
    Xyce::dout() << "  li_ContPos = " << li_ContPos << std::endl;
    Xyce::dout() << "  li_ContNeg = " << li_ContNeg << std::endl;
  }

  li_Bra = intLIDVec[0];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    Xyce::dout() << "  li_Bra = " << li_Bra << std::endl;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    Xyce::dout() << section_divider << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs( const std::vector<int> & stoLIDVecRef )
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());

  // Copy over the global ID lists:
  stoLIDVec = stoLIDVecRef;

  int i=0;
  li_storeVold =  stoLIDVec[i++];

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerBranchDataLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
//-----------------------------------------------------------------------------
/// Register the local store IDs
///
/// In addition to state vector, Xyce maintains a separate datastructure
/// called a "branch data" vector.  As with other such vectors, the device
/// declares at construction time how many branch vector entries it needs,
/// and later Topology assigns locations for devices, returning LIDs.
///
/// These LIDs are stored in this method for later use.
///
/// The Resistor device uses exactly one "branch data vector" element, where
///
/// @param stoLIDVecRef Store variable local IDs
///
/// @author Richard Schiek, Electrical Systems Modeling
/// @date   12/18/2012

void Instance::registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef)
{
  AssertLIDs(branchLIDVecRef.size() == getNumBranchDataVars());

  if (loadLeadCurrent)
  {
    li_branch_data= branchLIDVecRef[0];
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadNodeSymbols
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
//-----------------------------------------------------------------------------
void Instance::loadNodeSymbols(Util::SymbolTable &symbol_table) const
{
  addInternalNode(symbol_table, li_Bra, getName(), "branch");

  if (loadLeadCurrent)
  {
    addBranchDataNode( symbol_table, li_branch_data, getName(), "BRANCH_D");
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
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
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  APosEquBraVarOffset = jacLIDVec[0][0];
  ANegEquBraVarOffset = jacLIDVec[1][0];
  ABraEquPosNodeOffset = jacLIDVec[4][0];
  ABraEquNegNodeOffset = jacLIDVec[4][1];
  ABraEquContPosNodeOffset = jacLIDVec[4][2];
  ABraEquContNegNodeOffset = jacLIDVec[4][3];
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  f_PosEquBraVarPtr = &(dFdx[li_Pos][ APosEquBraVarOffset ]);
  f_NegEquBraVarPtr = &(dFdx[li_Neg][ ANegEquBraVarOffset ]);
  f_BraEquPosNodePtr = &(dFdx[li_Bra][ ABraEquPosNodeOffset ]);
  f_BraEquNegNodePtr = &(dFdx[li_Bra][ ABraEquNegNodeOffset ]);
  f_BraEquContPosNodePtr = &(dFdx[li_Bra][ ABraEquContPosNodeOffset ]);
  f_BraEquContNegNodePtr = &(dFdx[li_Bra][ ABraEquContNegNodeOffset ]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::acceptStep
// Purpose       : This function saves the value of the control voltage drop
//                 at the current accepted time.  It is to be called ONLY at the
//                 point when the time integrator has determined we've got a
//                 converged, acceptable solution and is accepting it,
//                 but before it's updated its times and rotated vectors.
//
// Special Notes : In SPICE this same stuff was done in the "TRAaccept" function.
//
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
//-----------------------------------------------------------------------------
void Instance::acceptStep()
{
  double * solVec = extData.nextSolVectorRawPtr;
  double * stoVec = extData.nextStoVectorRawPtr;

  stoVec[li_storeVold]=solVec[li_ContPos]-solVec[li_ContNeg];

}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 StepDelay instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;
  double * solVec = extData.nextSolVectorRawPtr;

  double v_pos = solVec[li_Pos];
  double v_neg = solVec[li_Neg];

  double i_bra = solVec[li_Bra];

  double * stoVec = extData.currStoVectorRawPtr;
  double v_drop=0.0;
  if ((getSolverState().dcopFlag))
  {
    v_drop=solVec[li_ContPos]-solVec[li_ContNeg];
  }

  fVec[li_Pos] += i_bra;
  fVec[li_Neg] += -i_bra;

  double src = (v_pos - v_neg) -v_drop ;
  fVec[li_Bra] += src;

  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    leadF[li_branch_data] = i_bra;
    double * junctionV = extData.nextJunctionVCompRawPtr;
    junctionV[li_branch_data] = src; 
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEBVector
//
// Purpose       : Loads the B-vector contributions for a single
//                 StepDelay instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
//-----------------------------------------------------------------------------
bool Instance::loadDAEBVector ()
{
  if (!(getSolverState().dcopFlag))
  {
    double * stoVec = extData.currStoVectorRawPtr;
    double v_drop;
    double * bVec = extData.daeBVectorRawPtr;
    v_drop = stoVec[li_storeVold];
    bVec[li_Bra] += v_drop;
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 resistor  instance.
//
// Special Notes : The F-vector is an algebraic constaint.
//
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Pos][APosEquBraVarOffset] += 1.0;
  dFdx[li_Neg][ANegEquBraVarOffset] -= 1.0;

  dFdx[li_Bra][ABraEquPosNodeOffset] += 1.0;
  dFdx[li_Bra][ABraEquNegNodeOffset] -= 1.0;
  if ((getSolverState().dcopFlag))
  {
    dFdx[li_Bra][ABraEquContPosNodeOffset] -= 1.0;
    dFdx[li_Bra][ABraEquContNegNodeOffset] += 1.0;
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  varTypeVec.resize(1);
  varTypeVec[0] = 'I';
}

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : "Model block" constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
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
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
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
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
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
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
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
  if (deviceMap.empty() || (deviceMap.find("STEPDELAY")!=deviceMap.end()))
  {

    Config<Traits>::addConfiguration()
      .registerDevice("stepdelay", 1)
      .registerModelType("stepdelay", 1);
  }
}

} // namespace StepDelay
} // namespace Device
} // namespace Xyce

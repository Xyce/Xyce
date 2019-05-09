//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose        : An Ideal OpAmp Model
//
// Special Notes  : This model assumes infinite gain
//                   and unbounded output voltages
//
// Creator        : Brian Fett, SNL
//
// Creation Date  : 07/27/05
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_OpAmp.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {

namespace OpAmp {

void Traits::loadInstanceParameters(ParametricData<OpAmp::Instance> &p)
{
  p.addPar ("FAKEPARAM", 0.0, &OpAmp::Instance::FAKEPARAM)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("");
}

void Traits::loadModelParameters(ParametricData<OpAmp::Model> &p)
{
}

std::vector< std::vector<int> > Instance::jacStamp;

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Brian Fett, SNL
// Creation Date : 08/05/07
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       iter,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(iter),
    outCurrent(0.0),
    deltaVoltage(0.0),
    
    FAKEPARAM(0.0),
    
    v_pos(0.0),
    v_neg(0.0),
    v_out(0.0),
    i_bra(0.0),

    li_Pos(-1),
    li_Neg(-1),
    li_Out(-1),
    li_Bra(-1),
    
    ABraEquPosNodeOffset(-1),
    ABraEquNegNodeOffset(-1),
    AOutEquBraVarOffset(-1)
{
  numIntVars   = 1;
  numExtVars   = 3;
  numStateVars = 0;

  if( jacStamp.empty() )
  {
    jacStamp.resize(4);    //    V1 V2 V3 I
    jacStamp[2].resize(1); // N1
    jacStamp[2][0] = 3;    // N2
    jacStamp[3].resize(2); // N3          X
    jacStamp[3][0] = 0;    // Br X  X
    jacStamp[3][1] = 1;
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // calculate dependent (ie computed) params and check for errors:

}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Brian Fett
// Creation Date : 07/28/05
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Brian Fett, SNL
// Creation Date : 08/05/05
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{

  return true;
}


// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Brian Fett, SNL
// Creation Date : 08/05/05
//-----------------------------------------------------------------------------
void Instance::registerLIDs ( const std::vector<int> & intLIDVecRef,
	                                const std::vector<int> & extLIDVecRef)
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  OpAmpInstance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

  if (numInt != numIntVars)
  {
    DevelFatal(*this).in("Instance::registerLIDs")
      << "numInt != numIntVars";
  }

  if (numExt != numExtVars)
  {
    DevelFatal(*this).in("Instance::registerLIDs")
      << "numExt != numExtVars";
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];
  li_Out = extLIDVec[2];
  li_Bra = intLIDVec[0];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl;
    Xyce::dout() << "  li_Neg = " << li_Neg << std::endl;
    Xyce::dout() << "  li_Out = " << li_Out << std::endl;
    Xyce::dout() << "  li_Bra = " << li_Bra << std::endl;
    Xyce::dout() << section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
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
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/21/02
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
// Creator       : Brian Fett, SNL
// Creation Date : 08/05/05
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  AOutEquBraVarOffset  = jacLIDVec[2][0];
  ABraEquPosNodeOffset = jacLIDVec[3][0];
  ABraEquNegNodeOffset = jacLIDVec[3][1];
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
  addInternalNode(symbol_table, li_Bra, getName(), "branch");
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Brian Fett, SNL
// Creation Date : 08/05/05
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;

  // get the value for v_pos, v_neg, v_out, i_bra

  v_pos = (*solVectorPtr)[li_Pos];
  v_neg = (*solVectorPtr)[li_Neg];
  v_out = (*solVectorPtr)[li_Out];
  i_bra = (*solVectorPtr)[li_Bra];

  outCurrent = i_bra;
  deltaVoltage = v_pos - v_neg;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "  Instance::updateIntermediateVars" << std::endl;
    Xyce::dout() << "  name = " << getName() <<std::endl;
    Xyce::dout() << "  v_pos  = " << v_pos << std::endl;
    Xyce::dout() << "  v_neg  = " << v_neg << std::endl;
    Xyce::dout() << "  v_out  = " << v_out << std::endl;
    Xyce::dout() << "  i_bra  = " << i_bra << std::endl;
    Xyce::dout() << std::endl;
    Xyce::dout() << subsection_divider << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/29/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "  Instance::updatePrimaryState" << std::endl;
  }
  bool bsuccess = true;
  bsuccess = updateIntermediateVars ();
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 OpAmp instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Brian Fett, SNL
// Creation Date : 08/05/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  // bool bsuccess = true;
  // int i_pos_rhs, i_neg_rhs, i_bra_rhs;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "  Instance::loadDAEFVector" << std::endl;
    Xyce::dout() << "  name       = " << getName() <<std::endl;
    Xyce::dout() << "  Output Current = " << outCurrent << std::endl;
    Xyce::dout() << "  Delta  Voltage = " << deltaVoltage << std::endl;
  }

  // Using values determined during the loadRHS function call.
  // (outCurrent, deltaVoltage).

  //(*extData.daeFVectorPtr)[li_Pos] += 0;
  //(*extData.daeFVectorPtr)[li_Neg] += 0;

  (*extData.daeFVectorPtr)[li_Out] += outCurrent;

  (*extData.daeFVectorPtr)[li_Bra] += deltaVoltage;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 resistor  instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Brian Fett, SNL
// Creation Date : 08/02/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  // bool bsuccess = true;

  Linear::Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "name = " << getName() << std::endl;
    Xyce::dout() << "\nOPAMP dFdx LOADS\n";
    Xyce::dout() << "Pos,Bra: " << li_Out << ",";
    Xyce::dout() << AOutEquBraVarOffset << ": " << 1.0 << std::endl;
    Xyce::dout() << "Bra,Pos: " << li_Bra << ",";
    Xyce::dout() << ABraEquPosNodeOffset << ": " << 1.0 << std::endl;
    Xyce::dout() << "Bra,Neg: " << li_Bra << ",";
    Xyce::dout() << ABraEquNegNodeOffset << ": " << -1.0 << std::endl;
    Xyce::dout() << "DONE OPAMP dFdx LOAD\n";
  }

  (*dFdxMatPtr)[li_Out][AOutEquBraVarOffset] += 1.0;
  (*dFdxMatPtr)[li_Bra][ABraEquPosNodeOffset] += 1.0;
  (*dFdxMatPtr)[li_Bra][ABraEquNegNodeOffset] -= 1.0;

  return true;
}

// end of new-DAE functions

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Brian Fett, SNL
// Creation Date : 08/05/05
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
  FAKEPARAM(0.0)
{
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
Model::~Model ()
{
  std::vector<Instance*>::iterator iter = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for ( ; iter!=last; ++iter)
  {
    delete (*iter);
  }
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  os << std::endl;
  os << "    name     model name  Parameters" << std::endl;
  for (int i=0; iter!=last; ++iter, ++i)
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
  if (deviceMap.empty() || (deviceMap.find("OPAMP")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("opamp", 1)
      .registerModelType("opamp", 1);
  }
}

} // namespace OpAmp
} // namespace Device
} // namespace Xyce

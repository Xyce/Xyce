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
// Purpose        : Implementation of the a non-ideal voltage source
//
// Special Notes  :
//
// Creator        : Richard Schiek, Electrical Models & Simulations
//
// Creation Date  : 10/23/14
//
//
//
//
//----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_DEV_Battery.h>

#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>
#include <Sacado_No_Kokkos.hpp>

namespace Xyce {
namespace Device {
namespace Battery {

//
// Basic equations for the battery are:
//
// V1 - V2 = deltaV = Ifactor( I ) * Tfactor( Temp ) * Vfxn( CapUsedEff )
//
// CapUsedEff = Icomp( I ) * Tcomp( Temp ) * CapUsed 
//
// Ifactor( I ) = Ifactor0 + Ifactor1 * (I-If0)
//
// Tfactor( T ) = Tfactor0 + Tfactor1 * (T - Tnom)
//
// Icomp( I ) = Icomp0 + Icomp1 * (I-Ic0)
//
// Tcomp( T ) = Tcomp0 + Tcomp1 * (T - Tnom)
//
// 

//
// Temp Compensation on Capacity
//
template <typename ScalarT>
void Tcomp( const ScalarT & T, const ScalarT & Tnom, const ScalarT & Tcomp0, const ScalarT & Tcomp1, ScalarT & fval )
{
  fval = Tcomp0 + Tcomp1 * (T - Tnom);
}



//
// Current Compensation on Capacity
//
template <typename ScalarT>
void Icomp( const ScalarT & I, const ScalarT & Ic0, const ScalarT & Icomp0, const ScalarT & Icomp1, ScalarT & fval )
{
  fval = Icomp0 + Icomp1 * (I - Ic0);
}

//
// Temperature Prefix Factor
//
template <typename ScalarT>
void Tfactor( const ScalarT & T, const ScalarT & TCnom, const ScalarT & Tfactor0, const ScalarT & Tfactor1, ScalarT & fval )
{
  fval = Tfactor0 + Tfactor1 * (T - TCnom);
}

//
// Current Prefix Factor
//
template <typename ScalarT>
void Ifactor( const ScalarT & I, const ScalarT & If0, const ScalarT & Ifactor0, const ScalarT & Ifactor1, ScalarT & fval )
{
  fval = Ifactor0 + Ifactor1 * (I - If0);
}

//
// Effective Capacity Used
//
template <typename ScalarT>
void CapUsedEff( const ScalarT & I, const ScalarT & T, const ScalarT & CapUsed, 
  const ScalarT & Ic0, const ScalarT & Icomp0, const ScalarT & Icomp1,
  const ScalarT & Tnom, const ScalarT & Tcomp0, const ScalarT & Tcomp1, ScalarT & fval )
{
  ScalarT Icp;
  Icomp( I, Ic0, Icomp0, Icomp1, Icp );
  ScalarT Tcp;
  Tcomp( T, Tnom, Tcomp0, Tcomp1, Tcp);
  fval = Icp  * Tcp * CapUsed;
}

//
// Voltage Equation 
//
template <typename ScalarT>
void VoltageEqu( const ScalarT & T, const ScalarT & I, const ScalarT & CapUsed, 
  const ScalarT & C0, const ScalarT & V0, const ScalarT & V1, const ScalarT & V2, const ScalarT & V3,
  const ScalarT & Ic0, const ScalarT & Icomp0, const ScalarT & Icomp1,
  const ScalarT & If0, const ScalarT & Ifactor0, const ScalarT & Ifactor1,
  const ScalarT & Tnom, const ScalarT & Tcomp0, const ScalarT & Tcomp1,
  const ScalarT & TCnom, const ScalarT & Tfactor0, const ScalarT & Tfactor1, ScalarT & fval )
{
  ScalarT effCapUsed;
  CapUsedEff( I, T, CapUsed, Ic0, Icomp0, Icomp1, Tnom, Tcomp0, Tcomp1, effCapUsed);
  
  ScalarT Ifc;
  Ifactor( I, If0, Ifactor0, Ifactor1, Ifc );
  ScalarT Tfc;
  Tfactor( T, TCnom, Tfactor0, Tfactor1, Tfc );
  
  ScalarT dCap = effCapUsed - C0;
  
  fval = V0 + V1 * dCap + V2 * std::pow( dCap, 2 ) + V3 * std::pow( dCap, 3 ); 
}



// Common Jacobian Stamp for all Battery devices.
// Because all batteries have identical Jacobian stamps, this data is
// declared static and is shared by all battery instances.
// 
std::vector<std::vector<int> > Instance::jacStamp;

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Instance::initializeJacobianStamp
// Purpose       : 
// Special Notes : 
// Scope         : private
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 2/11/2014
//-----------------------------------------------------------------------------
//
// @brief Common Jacobian stamp initializer for all Battery devices.
//
// The Jacobian stamp is a sparse-matrix representation of the pattern
// of non-zero elements that a device will put into the Jacobian matrix.
//
// The Jacobian stamp is used by the Topology package to determine indices
// into the full Jacobian matrix for elements that correspond to this 
// device.
//
// There is one row of the Jacobian stamp for each equation associated with
// a device.  The number of elements in a row are the number of non-zero
// elements in that row of the device's contribution to the Jacobian.
// The values of the elements are numbers local to the device that
// represent the column in which the non-zero element appears.
//
// For this battery, there are two external nodes (the positive and negative
// terminals of the device).  The positive node is referred to as the 0th
// node of the device, and the negative node the 1st node.
//
void Instance::initializeJacobianStamp()
{
  if (jacStamp.empty())
  {
    //      V1   V2   CellTemp   I   CapUsed 
    // V1   a    b       c       d     e
    // V2  -a   -b      -c      -d    -e
    // CT                f
    // I    g    h 
    // CU                        i 
    //
    
    jacStamp.resize(5);
    jacStamp[0].resize(5);
    jacStamp[0][0] = 0;
    jacStamp[0][1] = 1;
    jacStamp[0][2] = 2;
    jacStamp[0][3] = 3;
    jacStamp[0][4] = 4;
    
    jacStamp[1].resize(5);
    jacStamp[1][0] = 0;
    jacStamp[1][1] = 1;
    jacStamp[1][2] = 2;
    jacStamp[1][3] = 3;
    jacStamp[1][4] = 4;
    
    jacStamp[2].resize(1);
    jacStamp[1][0] = 2;
    
    jacStamp[3].resize(2);
    jacStamp[3][0] = 0;
    jacStamp[3][1] = 1;
    
    jacStamp[4].resize(1);
    jacStamp[0][0] = 3;
    
  }
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Traits::loadInstanceParameters
// Purpose       : 
// Special Notes : 
// Scope         : private
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Loads the parameter definition into the instance parameter map.
//
// @param p     instance parameter map
//
// Each parameter supported by the battery device instance is
// defined via the addPar call.  The minimum information required is
// the name of the parameter, the default value, and a member pointer
// to a variable in the instance class to which the value will be
// assigned.
//
// Additional information such as documentation for the parameter, units
// (used only in output of tables for documentation), whether a special
// boolean should be set if the parameter is given, etc. may be specified
// using the various set* calls of the Xyce::Device::Descriptor class.
//
// It is important to note that since parameters defined by addPar are
// defined as metadata separate from any instance, it is not possible to
// establish interrelationships between parameter defaults here.  Parameter
// defaults specified in addPar are always constants.  If a device requires
// that parameter defaults depend on values of other parameters, this has to
// be done in the instance constructor.  Examples of such parameters in this 
// device arethe "DTEMP" and "W" parameters, which are set to special defaults
// at device instantiation time.  Defaults for those parameters in the addPar
// calls (and hence in the LaTeX tables in the reference guide produced from
// this data) are misleading.
//
void Traits::loadInstanceParameters(ParametricData<Battery::Instance> &p)
{
  p.addPar("RCELL", 1.0, &Battery::Instance::Rcell)
    .setUnit(U_OHM)
    .setDescription("Cell resistance.");


}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Traits::loadModelParameters
// Purpose       : 
// Special Notes : The addPar calls here were refactored and moved here
//                 from the model constructor.  Those addPars had been
//                 in place from 2/4/2005.
// Scope         : private
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Loads the parameter definition into the model parameter map.
//
// @param p     model parameter map
//
// @see Xyce::Device::Battery::Traits::loadInstanceParameters
//
//
void Traits::loadModelParameters(ParametricData<Battery::Model> &p)
{

  // Create parameter definitions for parameter member variables
  p.addPar("VOFFSET", 1.0, &Battery::Model::VCapOffset_)
    .setUnit(U_VOLT)
    .setDescription("Battery voltage-capacity offset.");
  p.addPar("V0", 1.0, &Battery::Model::VCap0_)
    .setUnit(U_VOLT)
    .setDescription("Battery voltage at zero capacity used .");
  p.addPar("V1", 0.0, &Battery::Model::VCap1_)
    .setUnit(U_VOLT)
    .setDescription("Battery voltage first derivative");
  p.addPar("V2", 0.0, &Battery::Model::VCap2_)
    .setUnit(U_VOLT)
    .setDescription("Battery voltage second derivative");  
  p.addPar("V3", 0.0, &Battery::Model::VCap3_)
    .setUnit(U_NONE)
    .setDescription("Battery voltage third derivative");
  p.addPar("IC0", 1.0, &Battery::Model::IC0_)
    .setUnit(U_AMP)
    .setDescription("Current compensation factor offset.");
  p.addPar("ICOMP0", 1.0, &Battery::Model::IComp0_)
    .setUnit(U_AMP)
    .setDescription("Current compensation factor intercept.");
  p.addPar("ICOMP1", 1.0, &Battery::Model::IComp1_)
    .setUnit(U_NONE)
    .setDescription("Current compensation rate slope "); 
  p.addPar("TNOM", 1.0, &Battery::Model::Tnom_)
    .setUnit(U_DEGC)
    .setDescription("Temperature compensation offset. ");
  p.addPar("TCOMP0", 1.0, &Battery::Model::TComp0_)
    .setUnit(U_DEGC)
    .setDescription("Temperature compensation intercept.");
  p.addPar("TCOMP1", 1.0, &Battery::Model::TComp1_)
    .setUnit(U_NONE)
    .setDescription("Temperature compensation slope.");
  p.addPar("IF0", 1.0, &Battery::Model::IF0_)
    .setUnit(U_AMP)
    .setDescription("Current-Capacity factor offset.");
  p.addPar("IFACT0", 1.0, &Battery::Model::IFact0_)
    .setUnit(U_AMP)
    .setDescription("Current-Capacity factor intercept.");
  p.addPar("IFACT1", 1.0, &Battery::Model::IFact1_)
    .setUnit(U_NONE)
    .setDescription("Current-Capacity rate slope "); 
  p.addPar("TFACTNOM", 1.0, &Battery::Model::TFactNom_)
    .setUnit(U_DEGC)
    .setDescription("Temperature-Capacity offset ");
  p.addPar("TFACT0", 1.0, &Battery::Model::TFact0_)
    .setUnit(U_DEGC)
    .setDescription("Temperature-Capacity intercept.");
  p.addPar("TFACT1", 1.0, &Battery::Model::TFact1_)
    .setUnit(U_NONE)
    .setDescription("Temperature-Capacity  slope.");
    
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Instance::Instance
// Purpose       : Instance constructor
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Construct a battery instance.
//

Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & instance_block,
  Model &               model,
  const FactoryBlock &  factory_block)
  : DeviceInstance(instance_block, configuration.getInstanceParameters(), factory_block),
    model_(model),
    li_Pos(-1),
    li_Neg(-1),
    li_CellTemp(-1),
    li_current(-1),
    li_capUsed(-1), 
    li_branch_data(0),
    APosEquPosNodeOffset(-1),
    APosEquNegNodeOffset(-1),
    APosEquCellTempOffset(-1),
    APosEquCapUsedNodeOffset(-1),
    ANegEquPosNodeOffset(-1),
    ANegEquNegNodeOffset(-1),
    ANegEquCellTempOffset(-1),
    ANegEquCapUsedNodeOffset(-1),
    CellTempEquCellTempOffset(-1),
    CurrentEquPosNodeOffset(-1),
    CurrentEquNegNodeOffset(-1),
    CurrentEquCurrentOffset(-1)
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    ,
    f_PosEquPosNodePtr(0),
    f_PosEquNegNodePtr(0),
    f_PosEquCellTempPtr(0),
    f_PosEquCurrentPtr(0),
    f_PosEquCapUsedNodePtr(0),
    f_NegEquPosNodePtr(0),
    f_NegEquNegNodePtr(0),
    f_NegEquCellTempPtr(0),
    f_NegEquCurrentPtr(0),
    f_NegEquCapUsedNodePtr(0),
    f_CellTempEquCellTempPtr(0),
    f_CurrentEquPosNodePtr(0),
    f_CurrentEquNegNodePtr(0),
    f_CurrentEquCurrentPtr(0),
    f_CapUsedEquCurrentPtr(0),
    q_CapUsedEquCurrentPtr(0)
  // ,
#endif
    //resSens(*this)
{
  // Initialize DeviceInstance values
  numIntVars   = 2;    // Initialize number if internal nodes in DeviceInstance
  numExtVars   = 3;    // Initialize number if external nodes in DeviceInstance
  numStateVars = 0;    // Initialize number if state variables in DeviceInstance
  setNumStoreVars(1);  // Initialize number if store variables in DeviceInstance 

  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1;    // this is the space to allocate if lead current or power is needed.

  initializeJacobianStamp();

  // Set params to constant default values from parameter definition
  setDefaultParams();

  // Set params according to instance line and constant defaults from metadata
  setParams(instance_block.params);

  // Calculate any parameters specified as expressions
  updateDependentParameters();

  // Process the parameters to complete initialization
  processParams();
  
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Instance::processParams
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
// Process parameters.
//
// @return true on success
//
// In general, the processParams method is intended as a place for complex
// manipulation of parameters that must happen if temperature or other
// external changes happen.
//
bool Instance::processParams()
{
  // now set the temperature related stuff.
  double temp=0;
  return updateTemperature(temp);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Instance::registerLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Register local IDs
//
// Register the local internal and external node IDs.
//
// @param intLIDVecRef internal local IDs from topology package
// @param extLIDVecRef external local IDs from topology package
// 
void Instance::registerLIDs(
  const std::vector<int> & intLIDVecRef,
  const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  // Copy the local ID lists.
  intLIDVec = intLIDVecRef;                           // Set the internal local IDs in DeviceInstance
  extLIDVec = extLIDVecRef;                           // Set the external local IDs in DeviceInstance

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];
  li_CellTemp = extLIDVec[2];

  // This lid is for the internal variable for layer thickess that determines resistence
  li_current = intLIDVec[0];
  li_capUsed = intLIDVec[1];

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Instance::registerStateLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Register the local state IDs
//
// @note The battery does not have any state vars, so this function
// does nothing.
//
// @param staLIDVecRef State variable local IDs
//
// In general, devices may declare at construction that they require storage
// locations in the "state vector."  Topology assigns locations in the 
// state vector and returns "Local IDs" (LIDs) for devices to use for their
// state vector entries.  If a device uses state vector entries, it
// uses the registerStateLIDs method to save the local IDs for later use.
// 
// @author Robert Hoekstra, SNL, Parallel Computational Sciences
// @date   06/12/02

void Instance::registerStateLIDs(const std::vector<int> & staLIDVecRef)
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Instance::registerStoreLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 12/18/2012
//-----------------------------------------------------------------------------
// Register the local store IDs
//
// In addition to state vector, Xyce maintains a separate datastructure
// called a "store" vector.  As with other such vectors, the device
// declares at construction time how many store vector entries it needs,
// and later Topology assigns locations for devices, returning LIDs.
//
//
//

void Instance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef)
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Instance::registerBranchDataLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 12/18/2012
//-----------------------------------------------------------------------------
// Register the local store IDs
//
// In addition to state vector, Xyce maintains a separate datastructure
// called a "branch data" vector.  As with other such vectors, the device
// declares at construction time how many branch vector entries it needs,
// and later Topology assigns locations for devices, returning LIDs.
//
// These LIDs are stored in this method for later use.
//
// The battery device uses exactly one "branch data vector" element, where
// it keeps the "lead current" that may be used on .PRINT lines as
// "I(ybattery)" for the current through ybattery. and a junction voltage.
//
// @param stoLIDVecRef Store variable local IDs
//
// @author Richard Schiek, Electrical Systems Modeling
// @date   12/18/2012

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
// Creator       : Richard Schiek, Electrical Sysetms Modeling 
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
void Instance::loadNodeSymbols(Util::SymbolTable &symbol_table) const
{
  addInternalNode(symbol_table, li_current, getName(), "i");

  if (loadLeadCurrent)
  {
    addBranchDataNode( symbol_table, li_branch_data, getName(), "BRANCH_D");
  }
}




//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Instance::registerJacLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Register the Jacobian local IDs
//
// @param jacLIDVec Jacobian local Ids
//
// @see Xyce::Device::Battery::Instance::initializeJacobianStamp
//
// Having established local IDs for the solution variables, Topology must
// also assign local IDs for the elements of the Jacobian matrix.
//
// For each non-zero element that was identified in the jacobianStamp,
// Topology will assign a Jacobian LID.  The jacLIDVec will have the 
// same structure as the JacStamp, but the values in the jacLIDVec will
// be offsets into the row of the sparse Jacobian matrix corresponding
// to the non-zero identified in the stamp.
// 
// These offsets are stored in this method for use later when we load
// the Jacobian.
//
// @author Robert Hoekstra, SNL, Parallel Computational Sciences
// @date   08/27/01

void Instance::registerJacLIDs(const std::vector< std::vector<int> > & jacLIDVec)
{
  // Let DeviceInstance do its work.
  DeviceInstance::registerJacLIDs(jacLIDVec);

  // Store offsets of the components of the Jacobian of this instance
  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNegNodeOffset = jacLIDVec[0][1];
  APosEquCellTempOffset = jacLIDVec[0][2];
  APosEquCurrentOffset  = jacLIDVec[1][3];
  APosEquCapUsedNodeOffset = jacLIDVec[0][4];
  
  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];
  ANegEquCellTempOffset = jacLIDVec[1][2];
  ANegEquCurrentOffset  = jacLIDVec[1][3];
  ANegEquCapUsedNodeOffset = jacLIDVec[1][4];
  
  CellTempEquCellTempOffset = jacLIDVec[2][0];
  
  CurrentEquPosNodeOffset = jacLIDVec[3][0];
  CurrentEquNegNodeOffset = jacLIDVec[3][1];
  
  // this term is not requested in the jacstamp
  // CurrentEquCurrentOffset = jacLIDVec[3][2];
  
  CapUsedEquCurrentOffset = jacLIDVec[4][0];

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Instance::setupPointers
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Setup direct access pointer to solution matrix and vectors.
//
// @see Xyce::Device::Battery::Instance::registerJacLIDs
//
// As an alternative to the row offsets defined in registerJacLIDs, it 
// is also possible to obtain direct pointers of the Jacobian elements.
//
// This method uses the offsets obtained in registerJacLIDs to retrieve
// the pointers.
//
// In this device the pointers to the matrix are only saved
// (and are only used for matrix access) if
// Xyce_NONPOINTER_MATRIX_LOAD is NOT defined at compile time.  For
// some devices the use of pointers instead of array indexing can be
// a performance enhancement.
//
// Use of pointers in this device is disabled by defining
// Xyce_NONPOINTER_MATRIX_LOAD at compile time.  When disabled, array
// indexing with the offsets from registerJacLIDs is used in
// the matrix load methods.
//
// @author Eric Keiter, SNL
// @date   11/30/08

void Instance::setupPointers()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);
  
  f_PosEquPosNodePtr = &(dFdx[li_Pos][APosEquPosNodeOffset]);
  f_PosEquNegNodePtr = &(dFdx[li_Pos][APosEquNegNodeOffset]);
  f_PosEquCellTempPtr = &(dFdx[li_Pos][APosEquCellTempOffset]);
  f_PosEquCurrentPtr = &(dFdx[li_Pos][APosEquCurrentOffset]);
  f_PosEquCapUsedNodePtr = &(dFdx[li_Pos][APosEquCapUsedNodeOffset]);
  
  f_NegEquPosNodePtr = &(dFdx[li_Neg][ANegEquPosNodeOffset]);
  f_NegEquNegNodePtr = &(dFdx[li_Neg][ANegEquNegNodeOffset]);
  f_NegEquCellTempPtr = &(dFdx[li_Neg][ANegEquCellTempOffset]);
  f_NegEquCurrentPtr = &(dFdx[li_Neg][ANegEquCurrentOffset]);
  f_NegEquCapUsedNodePtr = &(dFdx[li_Neg][ANegEquCapUsedNodeOffset]);

  f_CellTempEquCellTempPtr = &(dFdx[li_CellTemp][CellTempEquCellTempOffset]);
  
  f_CurrentEquPosNodePtr = &(dFdx[li_current][CurrentEquPosNodeOffset]);
  f_CurrentEquNegNodePtr = &(dFdx[li_current][CurrentEquNegNodeOffset]);
  //f_CurrentEquCurrentPtr = &(dFdx[li_current][CurrentEquCurrentOffset]);

  f_CapUsedEquCurrentPtr = &(dFdx[li_capUsed][CapUsedEquCurrentOffset]);
  q_CapUsedEquCurrentPtr = &(dQdx[li_capUsed][CapUsedEquCurrentOffset]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Instance::updateIntermediateVars
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Update the intermediate variables
//
// @return true on success
//
// The bulk of any device's computation is handled in the instance class'
// updateIntermediateVars method.  For the battery, this is
// merely the computation of the current through the device given the
// voltage difference between its terminals.
//
// Intermediate variables computed here are used in the methods that
// load data into the F, Q, dFdX and dQdX data structures.
//
// @note This method is called by the updatePrimaryState
// method. Since the Battery class reimplements the "Master"
// "loadState" function that loads the contributions from all
// such devices in a single loop, THIS FUNCTION IS NOT ACTUALLY
// USED!  
//
//
bool Instance::updateIntermediateVars()
{
  double * solVec = extData.nextSolVectorRawPtr;
  // external variables
  double v_pos    = solVec[li_Pos];
  double v_neg    = solVec[li_Neg];
  double cellTemp = solVec[li_CellTemp];
  // internal variables
  double i0       = solVec[li_current];
  double capUsed  = solVec[li_capUsed];
  

  // update values that will be loaded into 
  // the F & Q vectors as well as the jacobian terms 
  
  // Voltage equation (V1 and V2)
  {
    // extra scope to limit variable extent.  Otherwise one will have lots of 
    // similar variables 
    
    // variables for which we need derivatives 
    Sacado::Fad::SFad<double,3> varCellTemp( 3, 0, cellTemp );
    Sacado::Fad::SFad<double,3> varCurrent( 3, 1, i0 );
    Sacado::Fad::SFad<double,3> varCapUsed( 3, 2, capUsed );
    
    // parameters
    Sacado::Fad::SFad<double,3> pC0( model_.VCapOffset_ );
    Sacado::Fad::SFad<double,3> pV0( model_.VCap0_ );
    Sacado::Fad::SFad<double,3> pV1( model_.VCap1_ );
    Sacado::Fad::SFad<double,3> pV2( model_.VCap2_ );
    Sacado::Fad::SFad<double,3> pV3( model_.VCap3_ );
    
    Sacado::Fad::SFad<double,3> pIc0( model_.IC0_ );
    Sacado::Fad::SFad<double,3> pIcomp0( model_.IComp0_ );
    Sacado::Fad::SFad<double,3> pIcomp1( model_.IComp1_ );
    
    Sacado::Fad::SFad<double,3> pIf0( model_.IF0_ );
    Sacado::Fad::SFad<double,3> pIfactor0( model_.IFact0_ );
    Sacado::Fad::SFad<double,3> pIfactor1( model_.IFact1_ );
    
    Sacado::Fad::SFad<double,3> pTnom( model_.Tnom_ );
    Sacado::Fad::SFad<double,3> pTcomp0( model_.TComp0_ );
    Sacado::Fad::SFad<double,3> pTcomp1( model_.TComp1_ );
    
    Sacado::Fad::SFad<double,3> pTFnom( model_.TFactNom_ );
    Sacado::Fad::SFad<double,3> pTfactor0( model_.TFact0_ );
    Sacado::Fad::SFad<double,3> pTfactor1( model_.TFact1_ );
    
    Sacado::Fad::SFad<double,3> resultFad;  
    
    VoltageEqu( varCellTemp, varCurrent,  varCapUsed, 
      pC0, pV0, pV1, pV2, pV3,
      pIc0, pIcomp0,  pIcomp1,
      pIf0, pIfactor0, pIfactor1,
      pTnom, pTcomp0, pTcomp1,
      pTFnom, pTfactor0, pTfactor1, 
      resultFad );
     
    v1FEqu = resultFad.val();
    dv1dv1FEqu = 1.0;
    dv2dv1FEqu = -1.0;
    dCellTempdv1FEqu = resultFad.dx(0);
    dIdv1Equ = resultFad.dx(1); 
    dCapUseddv1Equ = resultFad.dx(2);
    
    v2Equ = -resultFad.val();
    dv1dv2FEqu = -1.0;
    dv2dv2FEqu = 1.0;
    dCellTempdv2FEqu = -resultFad.dx(0);
    dIdv2Equ = -resultFad.dx(1); 
    dCapUseddv2Equ = -resultFad.dx(2);
    
  }
  
  // Cell Temp equation
  {
    cellTempEqu = 0.0;
    dCellTempdCellTempEqu = 1.0;
  }
  
  // Current equation 
  {
    currentEqu = (v_pos - v_neg)/Rcell;
    dV1dCurrentEqu = 1.0/Rcell;
    dV2dCurrentEqu = -1.0/Rcell;
  }
  
  // Capacity Used equation
  {
    capUsedEqu = getSolverState().pdt_ * i0;
    dCurrentdCapUsedEqu = getSolverState().pdt_;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Instance::updatePrimaryState
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Update the state variables.
//
// @return true on success
//
// This function is the function that is called from the device manager
// each time a load of vectors and Jacobian is required.  Its first
// job is to call updateIntermediateVars.
//
// After calling updateIntermediateVars, the updatePrimaryState method
// may load state vector elements as needed.
// 
// The battery device has no state vector elements, so all this method does
// in the battery instance class is call updateIntermediateVars.
//
// There is no longer a "secondary" state.  The "primary" in
// this method's name is purely historical.
//
// @note This method is called by the default implementation of the
// loadState master function. Since the Battery class reimplements
// the "Master" "loadState" function that loads the contributions
// from all battery devices in a single loop, THIS FUNCTION IS NOT
// ACTUALLY USED, NOR is the updateIntermediateVars method it calls!
// The updatePrimaryState method is only called when a device class
// does not re-implement the master class.  This can be a source of
// confusion when attempting to modify the Battery device, or any
// other device That reimplements the Master classes.
//
// @see Xyce::Device::Battery::Master::updateState
//
// @author Eric Keiter, SNL, Parallel Computational Sciences
// @date   01/29/01
//
bool Instance::updatePrimaryState()
{
  return updateIntermediateVars();
}

bool Instance::loadDAEQVector()
{
  double * qVec = extData.daeQVectorRawPtr;
  double * solVec = extData.nextSolVectorRawPtr;
  

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Instance::loadDAEFVector
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Load the DAE F vector.
//
// @return true on success
//
// The Xyce DAE formulation solves the differential-algebraic
// equations \f$F(x)+dQ(x)/dt=0\f$ These are vector equations
// resulting from the modified nodal analysis of the circuit.
// 
// This method loads the F-vector contributions for a single battery
// instance.
//
// In this method, the offsets defined in registerLIDs are used to
// store the device's F contributions into the F vector.
//
// The Q vector is used for devices that store charge or magnetic
// energy, and since the battery does none of that, the F vector
// contributions are the whole of the battery's contribution to the
// full set of DAEs.
//
// @note This method is called by the default implementation of the
// loadDAEVectors master function. Since the Battery class
// reimplements the "Master" "loadDAEVectors" function that loads the
// contributions from all battery devices in a single loop, THIS
// FUNCTION IS NOT ACTUALLY USED.  The loadDAEFVector method is only
// called when a device class does not re-implement the master class.
// This can be a source of confusion when attempting to modify the Battery
// device, or any other device that reimplements the Master classes.
//
// @see Xyce::Device::Battery::Master::loadDAEVectors
//
// @author Eric Keiter, SNL, Parallel Computational Sciences
// @date   01/24/03
//
bool Instance::loadDAEFVector()
{
  double * fVec = extData.daeFVectorRawPtr;
  double * solVec = extData.nextSolVectorRawPtr;
  fVec[li_Pos] += i0;
  fVec[li_Neg] += -i0;

  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    double * junctionV = extData.nextJunctionVCompRawPtr;
    
    leadF[li_branch_data] = i0;
    junctionV[li_branch_data] = solVec[li_Pos] - solVec[li_Neg];
  }
  
  fVec[li_CellTemp] += cellTempEqu;
  fVec[li_current] += currentEqu;
  fVec[li_capUsed] += solVec[li_capUsed] - capUsedEqu;

  return true;
}

bool Instance::loadDAEdQdx()
{
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);
  
  return true;
}


//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Instance::loadDAEdFdx
// Purpose       : Loads the F-vector contributions for a single
//                 battery  instance.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Load the DAE the derivative of the F vector with respect to the
// solution vector x, dFdx
//
// Loads the contributions for a single battery instance to the 
// dFdx matrix (the F contribution to the Jacobian).
//
// This method uses the Jacobian LIDs (row offsets) that were stored by
// registerJacLIDs.
//
// @see Xyce::Device::Battery::Instance::registerJacLIDs
//
// @note This method is called by the default implementation of the
// loadDAEMatrices master function. Since the Battery class
// reimplements the "Master" "loadDAEMatrices" function that loads the
// contributions from all battery devices in a single loop, THIS
// FUNCTION IS NOT ACTUALLY USED.  The loadDAEdFdx method is only
// called when a device class does not re-implement the master class.
// This can be a source of confusion when attempting to modify the Battery
// device, or any other device that reimplements the Master classes.
//
// @see Xyce::Device::Battery::Master::loadDAEMatrices
//
// @return true on success
//
// @author Eric Keiter, SNL, Parallel Computational Sciences
// @date   03/05/04
bool Instance::loadDAEdFdx()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  dFdx[li_Pos][APosEquPosNodeOffset] += dv1dv1FEqu;
  dFdx[li_Pos][APosEquNegNodeOffset] += dv2dv1FEqu;
  dFdx[li_Pos][APosEquCellTempOffset] += dCellTempdv1FEqu;
  dFdx[li_Pos][APosEquCurrentOffset] += dIdv1Equ;
  dFdx[li_Pos][APosEquCapUsedNodeOffset]   += dCapUseddv1Equ;
  
  dFdx[li_Neg][ANegEquPosNodeOffset] += dv1dv2FEqu;
  dFdx[li_Neg][ANegEquNegNodeOffset] += dv2dv2FEqu;
  dFdx[li_Neg][ANegEquCellTempOffset] += dCellTempdv2FEqu;
  dFdx[li_Neg][ANegEquCurrentOffset] += dIdv2Equ;
  dFdx[li_Neg][ANegEquCapUsedNodeOffset]   += dCapUseddv2Equ;
  
  dFdx[li_CellTemp][CellTempEquCellTempOffset] += dCellTempdCellTempEqu;
  
  dFdx[li_current][CurrentEquPosNodeOffset] = dV1dCurrentEqu;
  dFdx[li_current][CurrentEquNegNodeOffset] = dV2dCurrentEqu;
  //dFdx[li_current][CurrentEquCurrentOffset] = 0.0;
  
  dFdx[li_capUsed][CapUsedEquCurrentOffset] += dCurrentdCapUsedEqu;
  
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Update the parameters that depend on the temperature of the device
//
// @param temp_tmp temperature
//
// Xyce has a number of mechanisms that allow temperature to be changed
// after a device has been instantiated.  These include .STEP loops over
// temperature.  When temperature is changed, any device that has parameters
// that depend on temperature must be updated.  That updating happens here.
//
// The Battery device supports temperature-dependent resistance through its
// TC1 (linear dependence) and TC2 (quadratic dependence) parameters.
// If these parameters are specified, the resistance must be updated.
//
// @return true on success
//
// @author Tom Russo, Component Information and Models
// @date   02/27/01
bool Instance::updateTemperature(const double & temp_tmp)
{
  bool bsuccess = true;
  
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Process model parameters
//
// @return true on success
//
// @author Eric Keiter, SNL, Parallel Computational Sciences
// @date   6/03/02
bool Model::processParams()
{
  return true;
}

//----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Model::processInstanceParams
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//----------------------------------------------------------------------------
//
// Process the instance parameters of instance owned by this model
//
// This method simply loops over all instances associated with this
// model and calls their processParams method.
//
// @return true
//
// @author Dave Shirely, PSSI
// @date   03/23/06

bool Model::processInstanceParams()
{
  for (InstanceVector::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
  {
    (*it)->processParams();
  }

  return true;
}
//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Construct a battery model from a "model block" that was created
// by the netlist parser.
//
// @param configuration
// @param model_block
// @param factory_block
//
// @author Eric Keiter, SNL, Parallel Computational Sciences
// @date   5/16/00
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    model_block,
  const FactoryBlock &  factory_block)
  : DeviceModel(model_block, configuration.getModelParameters(), factory_block),
  VCap0_(0.0),
  VCap1_(0.0),
  VCap2_(0.0),
  VCap3_(0.0),
  IComp0_(0.0),
  IComp1_(0.0),
  Tnom_(0.0),
  TComp0_(0.0),
  TComp1_(0.0)

{
  // Set params to constant default values.
  setDefaultParams();

  // Set params according to .model line and constant defaults from metadata.
  setModParams(model_block.params);

  // Set any non-constant parameter defaults.
  //if (!given("TNOM"))
  //  tnom = getDeviceOptions().tnom;

  // Calculate any parameters specified as expressions.
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors.
  processParams();


}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Destroy this model.
//
// Also destroys all instances that use this model.
//
// @author Eric Keiter, SNL, Parallel Computational Sciences
// @date   3/16/00
Model::~Model()
{
  // Destory all owned instances
  for (InstanceVector::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
  {
    delete (*it);
  }

}

//-----------------------------------------------------------------------------
// Function      : N_DEV_BatteryModel::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Print instances associated with this model.
//
// Used only for debugging
//
// @param os output stream
//
// @return reference to output stream
//
// @author Eric Keiter, SNL, Parallel Computational Sciences
// @date   4/03/00
//
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  os << std::endl;
  os << "Number of Battery Instances: " << instanceContainer.size() << std::endl;
  os << "    name     model name  Parameters" << std::endl;

  int i = 0;
  for (InstanceVector::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
  {
    os << "  " << i << ": " << (*it)->getName() << "\t";
    os << getName();
    //os << "\t\tR(Tnom) = " << (*it)->R;
    os << "\tG(T) = " << (*it)->G;
    os << std::endl;
    ++i;
  }

  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Model::forEachInstance
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
// Apply a device instance "op" to all instances associated with this
// model
// 
// @param[in] op Operator to apply to all instances.
// 
// 
void Model::forEachInstance(DeviceInstanceOp &op) const /* override */ 
{
  for (std::vector<Instance *>::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
    op(*it);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Master::updateState
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Update state for all battery instances, regardless of model.
//
// @param solVec solution vector
// @param staVec state vector
// @param stoVec store vector
//
// @return true on success
//
// @note Because the battery device re-implements the base-class
// Master::updateState, the Instance::updatePrimaryState method is never
// called, nor is the Instance::updateIntermediateVars method.  This method
// replaces those, and does the same work but inside a loop over all
// battery instances.
//
// @see Xyce::Device::Battery::Instance::updatePrimaryState
// @author Eric Keiter, SNL
// @date   11/26/08

bool Master::updateState(double * solVec, double * staVec, double * stoVec)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & bi = *(*it);
    
    // external variables
    double v_pos    = solVec[bi.li_Pos];
    double v_neg    = solVec[bi.li_Neg];
    double cellTemp = solVec[bi.li_CellTemp];
    // internal variables
    double capUsed  = solVec[bi.li_capUsed];
    bi.i0    = solVec[bi.li_current];

    // update values that will be loaded into 
    // the F & Q vectors as well as the jacobian terms 
  
    // Voltage equation (V1 and V2)
    {
      // extra scope to limit variable extent.  Otherwise one will have lots of 
      // similar variables 
    
      // variables for which we need derivatives 
      Sacado::Fad::SFad<double,3> varCellTemp( 3, 0, cellTemp );
      Sacado::Fad::SFad<double,3> varCurrent( 3, 1, bi.i0 );
      Sacado::Fad::SFad<double,3> varCapUsed( 3, 2, capUsed );
    
      // parameters
      Sacado::Fad::SFad<double,3> pC0( bi.model_.VCapOffset_ );
      Sacado::Fad::SFad<double,3> pV0( bi.model_.VCap0_ );
      Sacado::Fad::SFad<double,3> pV1( bi.model_.VCap1_ );
      Sacado::Fad::SFad<double,3> pV2( bi.model_.VCap2_ );
      Sacado::Fad::SFad<double,3> pV3( bi.model_.VCap3_ );
    
      Sacado::Fad::SFad<double,3> pIc0( bi.model_.IC0_ );
      Sacado::Fad::SFad<double,3> pIcomp0( bi.model_.IComp0_ );
      Sacado::Fad::SFad<double,3> pIcomp1( bi.model_.IComp1_ );
    
      Sacado::Fad::SFad<double,3> pIf0( bi.model_.IF0_ );
      Sacado::Fad::SFad<double,3> pIfactor0( bi.model_.IFact0_ );
      Sacado::Fad::SFad<double,3> pIfactor1( bi.model_.IFact1_ );
    
      Sacado::Fad::SFad<double,3> pTnom( bi.model_.Tnom_ );
      Sacado::Fad::SFad<double,3> pTcomp0( bi.model_.TComp0_ );
      Sacado::Fad::SFad<double,3> pTcomp1( bi.model_.TComp1_ );
    
      Sacado::Fad::SFad<double,3> pTFnom( bi.model_.TFactNom_ );
      Sacado::Fad::SFad<double,3> pTfactor0( bi.model_.TFact0_ );
      Sacado::Fad::SFad<double,3> pTfactor1( bi.model_.TFact1_ );
    
      Sacado::Fad::SFad<double,3> resultFad;  
    
      VoltageEqu( varCellTemp, varCurrent,  varCapUsed, 
        pC0, pV0, pV1, pV2, pV3,
        pIc0, pIcomp0,  pIcomp1,
        pIf0, pIfactor0, pIfactor1,
        pTnom, pTcomp0, pTcomp1,
        pTFnom, pTfactor0, pTfactor1, 
        resultFad );
     
      bi.v1FEqu = resultFad.val();
      bi.dv1dv1FEqu = 1.0;
      bi.dv2dv1FEqu = -1.0;
      bi.dCellTempdv1FEqu = resultFad.dx(0);
      bi.dIdv1Equ = resultFad.dx(1); 
      bi.dCapUseddv1Equ = resultFad.dx(2);
    
      bi.v2Equ = -resultFad.val();
      bi.dv1dv2FEqu = -1.0;
      bi.dv2dv2FEqu = 1.0;
      bi.dCellTempdv2FEqu = -resultFad.dx(0);
      bi.dIdv2Equ = -resultFad.dx(1); 
      bi.dCapUseddv2Equ = -resultFad.dx(2);
    
    }
  
    // Cell Temp equation
    {
      bi.cellTempEqu = 0.0;
      bi.dCellTempdCellTempEqu = 1.0;
    }
  
    // Current equation 
    {
      bi.currentEqu = (v_pos - v_neg)/bi.Rcell;
      bi.dV1dCurrentEqu = 1.0/bi.Rcell;
      bi.dV2dCurrentEqu = -1.0/bi.Rcell;
    }
  
    // Capacity Used equation
    {
      bi.capUsedEqu = getSolverState().pdt_ * bi.i0;
      bi.dCurrentdCapUsedEqu = getSolverState().pdt_;
    }
    
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Master::loadDAEVectors
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Load DAE vectors of all battery instances, regardless of model
//
// @param solVec solution vector
// @param fVec f vector
// @param qVec q vector
// @param leadF store lead current f vector
// @param leadQ store lead current q vector
//
// @return true on success
//
// @note Because the battery device re-implements the base-class
// Master::loadDAEVectors, the Instance::loadDAEFVector method is
// never called.  This method replaces those, and does the same work
// but inside a loop over all battery instances.
//
// @see Xyce::Device::Battery::Instance::loadDAEFVector
//
// @author Eric Keiter, SNL
// @date   11/26/08

bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, double * leadF, double * leadQ, double * junctionV)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & bi = *(*it);
    if( bi.loadLeadCurrent )
    {
      leadF[bi.li_branch_data] = bi.i0;
      junctionV[bi.li_branch_data] = solVec[bi.li_Pos] - solVec[bi.li_Neg];
    }
  
    fVec[bi.li_Pos] += bi.i0;
    fVec[bi.li_Neg] -= bi.i0;

    fVec[bi.li_CellTemp] += bi.cellTempEqu;
    fVec[bi.li_current] += bi.currentEqu;
    fVec[bi.li_capUsed] += solVec[bi.li_capUsed] - bi.capUsedEqu;
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Master::loadDAEMatrices
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Load DAE matrices for all battery instances, regardless of model
//
// @param dFdx matrix of derivatives of F vector with respect to solution
// @param dQdx matrix of derivatives of Q vector with respect to solution
//
// @return true on success
//
// @note Because the battery device re-implements the base-class
// Master::loadDAEMatrices, the Instance::loadDAEdFdx method is
// never called.  This method replaces those, and does the same work
// but inside a loop over all battery instances.
//
// @see Xyce::Device::Battery::Instance::loadDAEdFdx
//
// @author Eric Keiter, SNL
// @date   11/26/08

bool Master::loadDAEMatrices(Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & bi = *(*it);

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    *(bi.f_PosEquPosNodePtr) += bi.dv1dv1FEqu;
    *(bi.f_PosEquNegNodePtr) += bi.dv2dv1FEqu;
    *(bi.f_PosEquCellTempPtr) += bi.dCellTempdv1FEqu;
    *(bi.f_PosEquCurrentPtr)  += bi.dIdv1Equ;
    *(bi.f_PosEquCapUsedNodePtr) += bi.dCapUseddv1Equ;
    
    *(bi.f_NegEquPosNodePtr) -= bi.dv1dv2FEqu;
    *(bi.f_NegEquNegNodePtr) += bi.dv2dv2FEqu;
    *(bi.f_NegEquCellTempPtr ) += bi.dCellTempdv2FEqu;
    *(bi.f_NegEquCurrentPtr) += bi.dIdv2Equ;
    *(bi.f_NegEquCapUsedNodePtr) += bi.dCapUseddv2Equ;
    
    *(bi.f_CellTempEquCellTempPtr) += bi.dCellTempdCellTempEqu;
    
    *(bi.f_CurrentEquPosNodePtr) += bi.dV1dCurrentEqu;
    *(bi.f_CurrentEquNegNodePtr) += bi.dV2dCurrentEqu;
    //*(bi.f_CurrentEquCurrentPtr) += 0.0;
    
    *(bi.f_CapUsedEquCurrentPtr ) += bi.dCurrentdCapUsedEqu;
#else   
    dFdx[bi.li_Pos][bi.APosEquPosNodeOffset] += bi.dv1dv1FEqu;
    dFdx[bi.li_Pos][bi.APosEquNegNodeOffset] += bi.dv2dv1FEqu;
    dFdx[bi.li_Pos][bi.APosEquCellTempOffset] += bi.dCellTempdv1FEqu;
    dFdx[bi.li_Pos][bi.APosEquCurrentOffset] += bi.dIdv1Equ;
    dFdx[bi.li_Pos][bi.APosEquCapUsedNodeOffset]   += bi.dCapUseddv1Equ;
  
    dFdx[bi.li_Neg][bi.ANegEquPosNodeOffset] += bi.dv1dv2FEqu;
    dFdx[bi.li_Neg][bi.ANegEquNegNodeOffset] += bi.dv2dv2FEqu;
    dFdx[bi.li_Neg][bi.ANegEquCellTempOffset] += bi.dCellTempdv2FEqu;
    dFdx[bi.li_Neg][bi.ANegEquCurrentOffset] += bi.dIdv2Equ;
    dFdx[bi.li_Neg][bi.ANegEquCapUsedNodeOffset]   += bi.dCapUseddv2Equ;
  
    dFdx[bi.li_CellTemp][bi.CellTempEquCellTempOffset] += bi.dCellTempdCellTempEqu;
  
    dFdx[bi.li_current][bi.CurrentEquPosNodeOffset] = bi.dV1dCurrentEqu;
    dFdx[bi.li_current][bi.CurrentEquNegNodeOffset] = bi.dV2dCurrentEqu;
    // dFdx[bi.li_current][bi.CurrentEquCurrentOffset] = 0.0;
  
    dFdx[bi.li_capUsed][bi.CapUsedEquCurrentOffset] += bi.dCurrentdCapUsedEqu;

#endif
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::Traits::factory
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Create a new instance of the Battery device.
//
// @param configuration
// @param factory_block
//
Device *
Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{
  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Battery::registerDevice
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Define how to use the device in a netlist.
//
// This method is called from the Xyce::Device::registerOpenDevices
// function, which in turn is called by the device manager.
//
// The device is declared here to be an "battery" device, which must 
// take a model card of type "battery".  This device will correspond to model
// level 3 of battery models.
void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() || (deviceMap.find("BATTERY")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("battery", 1)
      .registerModelType("battery", 1);
  }
}

//-----------------------------------------------------------------------------
// Function      : BatterySensitivity::operator
// Purpose       : produces df/dp and dq/dp, where p=R.  
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
void BatterySensitivity::operator()(
    const ParameterBase &entity,
    const std::string & name,
    std::vector<double> & dfdp, 
    std::vector<double> & dqdp, 
    std::vector<double> & dbdp, 
    std::vector<int> & Findices,
    std::vector<int> & Qindices,
    std::vector<int> & Bindices
    ) const
{
  const ParameterBase * e1 = &entity;
  const Instance * in = dynamic_cast<const Instance *> (e1);

  double * solVec = in->extData.nextSolVectorRawPtr;
  double v_pos = solVec[in->li_Pos];
  double v_neg = solVec[in->li_Neg];

  double dfdpLoc = -(v_pos-v_neg)*in->G*in->G;

  dfdp.resize(2);
  dfdp[0] = +dfdpLoc;
  dfdp[1] = -dfdpLoc;

  Findices.resize(2);
  Findices[0] = in->li_Pos;
  Findices[1] = in->li_Neg;
}

} // namespace Battery
} // namespace Device
} // namespace Xyce

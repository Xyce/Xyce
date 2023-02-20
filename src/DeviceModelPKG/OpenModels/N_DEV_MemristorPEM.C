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

//----------------------------------------------------------------------------
//
// Purpose        : Implementation of the Sandia memristor model.  See.
//                  Sandia: Threshold Adaptive Memristor Model
//                  Shahar Kvatinsky, Eby G. Friedman, Uri Weiser,
//                  IEEE Transactions on Circuits and Systems-I Vol 60, No. 1, Jan. 2013
//                  model details.
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

#include <N_DEV_MemristorPEM.h>

#include <fstream>

#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>
#include <Sacado_No_Kokkos.hpp>

namespace Xyce {
namespace Device {
namespace MemristorPEM {


//
// voltage driven current function
//
template <typename ScalarT>
void h( const ScalarT & V1, const ScalarT & V2, const double paramV1, const double paramV2, const double paramG0, const double paramI1, const double paramI2, ScalarT & fval )
{
  fval = paramI1 * exp((V1-V2)/paramV1) - paramI2 * exp(-(V1-V2)/paramV2) + paramG0*(V1-V2) - (paramI1-paramI2);
}

//
// Threshold voltage function: G(V1, V2, Vp, Vn, d1, d2)
//
template <typename ScalarT>
void G( const ScalarT & V1, const ScalarT & V2, const double paramD1, const double paramD2, const double paramVp, const double paramVn, ScalarT & fval )
{
  if( (V1-V2) >= paramVp )
  {
    fval = (exp(paramD1 * (V1-V2-paramVp)) -1.0);
  }
  else if ( (V1-V2) <= paramVn )
  {
    fval = (exp(paramD2 * (V1-V2-paramVn)) -1.0);
  }
  else
  {
    fval = 0.0;
  }

}

//
// Device F Equation
//
//
// I = x h(V1,V2)
//
template <typename ScalarT>
void I_V_Fxn( const ScalarT & V1, const ScalarT & V2, const ScalarT & X, const double paramV1, const double paramV2, const double paramG0, const double paramI1, const double paramI2, ScalarT & fval )
{
  ScalarT hV;
  h( V1, V2, paramV1, paramV2, paramG0, paramI1, paramI2, hV );

  fval = X * hV;
}



// Common Jacobian Stamp for all MemristorPEM devices.
// Because all memristors have identical Jacobian stamps, this data is
// declared static and is shared by all memristor instances.
//
std::vector<std::vector<int> > Instance::jacStamp;

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Instance::initializeJacobianStamp
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 2/11/2014
//-----------------------------------------------------------------------------
///
/// @brief Common Jacobian stamp initializer for all MemristorPEM devices.
///
/// The Jacobian stamp is a sparse-matrix representation of the pattern
/// of non-zero elements that a device will put into the Jacobian matrix.
///
/// The Jacobian stamp is used by the Topology package to determine indices
/// into the full Jacobian matrix for elements that correspond to this
/// device.
///
/// There is one row of the Jacobian stamp for each equation associated with
/// a device.  The number of elements in a row are the number of non-zero
/// elements in that row of the device's contribution to the Jacobian.
/// The values of the elements are numbers local to the device that
/// represent the column in which the non-zero element appears.
///
/// For this memristor, there are two external nodes (the positive and negative
/// terminals of the device).  The positive node is referred to as the 0th
/// node of the device, and the negative node the 1st node.
///
void Instance::initializeJacobianStamp()
{
  if (jacStamp.empty())
  {
    jacStamp.resize(3);
    jacStamp[0].resize(3);
    jacStamp[0][0] = 0;
    jacStamp[0][1] = 1;
    jacStamp[0][2] = 2;
    jacStamp[1].resize(3);
    jacStamp[1][0] = 0;
    jacStamp[1][1] = 1;
    jacStamp[1][2] = 2;
    jacStamp[2].resize(3);
    jacStamp[2][0] = 0;
    jacStamp[2][1] = 1;
    jacStamp[2][2] = 2;
  }
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Traits::loadInstanceParameters
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Loads the parameter definition into the instance parameter map.
///
/// @param p     instance parameter map
///
/// Each parameter supported by the memristor device instance is
/// defined via the addPar call.  The minimum information required is
/// the name of the parameter, the default value, and a member pointer
/// to a variable in the instance class to which the value will be
/// assigned.
///
/// Additional information such as documentation for the parameter, units
/// (used only in output of tables for documentation), whether a special
/// boolean should be set if the parameter is given, etc. may be specified
/// using the various set* calls of the Xyce::Device::Descriptor class.
///
/// It is important to note that since parameters defined by addPar are
/// defined as metadata separate from any instance, it is not possible to
/// establish interrelationships between parameter defaults here.  Parameter
/// defaults specified in addPar are always constants.  If a device requires
/// that parameter defaults depend on values of other parameters, this has to
/// be done in the instance constructor.  Examples of such parameters in this
/// device arethe "DTEMP" and "W" parameters, which are set to special defaults
/// at device instantiation time.  Defaults for those parameters in the addPar
/// calls (and hence in the LaTeX tables in the reference guide produced from
/// this data) are misleading.
///
void Traits::loadInstanceParameters(ParametricData<MemristorPEM::Instance> &p)
{
  p.addPar("XO", 0.0, &MemristorPEM::Instance::XO_)
    .setGivenMember(&MemristorPEM::Instance::XOGiven_)
    .setUnit(U_NONE)
    .setDescription("Initial value for internal variable x");

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Traits::loadModelParameters
// Purpose       :
// Special Notes : The addPar calls here were refactored and moved here
//                 from the model constructor.  Those addPars had been
//                 in place from 2/4/2005.
// Scope         : private
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Loads the parameter definition into the model parameter map.
///
/// @param p     model parameter map
///
/// @see Xyce::Device::MemristorPEM::Traits::loadInstanceParameters
///
///
void Traits::loadModelParameters(ParametricData<MemristorPEM::Model> &p)
{

  // Create parameter definitions for parameter member variables
  p.addPar("V1", 1.0, &MemristorPEM::Model::V1_)
    .setUnit(U_VOLT)
    .setDescription("Voltage Scale factor.");
  p.addPar("V2", 1.0, &MemristorPEM::Model::V2_)
    .setUnit(U_VOLT)
    .setDescription("Voltage Scale factor.");
  p.addPar("I1", 1.0, &MemristorPEM::Model::I1_)
    .setUnit(U_AMP)
    .setDescription("Current Scale factor.");
  p.addPar("I2", 1.0, &MemristorPEM::Model::I2_)
    .setUnit(U_AMP)
    .setDescription("Current Scale factor.");
  p.addPar("G0", 1.0, &MemristorPEM::Model::G0_)
    .setUnit(U_NONE)
    .setDescription("Conductance factor.");
  p.addPar("VP", 1.0e-2, &MemristorPEM::Model::Vp_)
    .setUnit(U_VOLT)
    .setDescription("Positive Voltage Threshold");
  p.addPar("VN", -1.0e-2, &MemristorPEM::Model::Vn_)
    .setUnit(U_VOLT)
    .setDescription("Negative Voltage Threshold");
  p.addPar("D1", 1.0, &MemristorPEM::Model::D1_)
    .setUnit(U_NONE)
    .setDescription("Positive Voltage Threshold Magnitude Parameter");
  p.addPar("D2", 1.0, &MemristorPEM::Model::D2_)
    .setUnit(U_NONE)
    .setDescription("Negative Voltage Threshold Magnitude Parameter");
  p.addPar("C1", 1.0, &MemristorPEM::Model::C1_)
    .setUnit(U_NONE)
    .setDescription("State variable proportionality parameter for forward bias");
  p.addPar("C2", 1.0, &MemristorPEM::Model::C2_)
    .setUnit(U_NONE)
    .setDescription("State variable proportionality parameter for negative bias");
  p.addPar("FXPDATA", "filep.dat", &MemristorPEM::Model::fxPDataFileName_)
    .setUnit(U_NONE)
    .setDescription("File from which to read x,f+(x) data");
  p.addPar("FXMDATA", "filem.dat", &MemristorPEM::Model::fxMDataFileName_)
    .setUnit(U_NONE)
    .setDescription("File from which to read x,f-(x) data");

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Instance::Instance
// Purpose       : Instance constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Construct a memristor instance.
///

Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & instance_block,
  Model &               model,
  const FactoryBlock &  factory_block)
  : DeviceInstance(instance_block, configuration.getInstanceParameters(), factory_block),
    model_(model),
    XO_(0.0),
    G(0.0),
    i0(0.0),
    li_Pos(-1),
    li_Neg(-1),
    li_x(-1),
    li_store_tdt(-1),
    li_branch_data(0),
    APosEquPosNodeOffset(-1),
    APosEquNegNodeOffset(-1),
    APosEquXNodeOffset(-1),
    ANegEquPosNodeOffset(-1),
    ANegEquNegNodeOffset(-1),
    ANegEquXNodeOffset(-1),
    XEquVPosOffset(-1),
    XEquVNegOffset(-1),
    XEquXOffset(-1)
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    ,
    f_PosEquPosNodePtr(0),
    f_PosEquNegNodePtr(0),
    f_PosEquXNodePtr(0),
    f_NegEquPosNodePtr(0),
    f_NegEquNegNodePtr(0),
    f_NegEquXNodePtr(0),
    f_XEquPosNodePtr(0),
    f_XEquNegNodePtr(0),
    f_XEquXNodePtr(0),
    q_XEquXNodePtr(0)
  // ,
#endif
    //resSens(*this)
{
  // Initialize DeviceInstance values
  numIntVars   = 1;    // Initialize number if internal nodes in DeviceInstance
  numExtVars   = 2;    // Initialize number if external nodes in DeviceInstance
  numStateVars = 0;    // Initialize number if state variables in DeviceInstance
  setNumStoreVars(2);  // Initialize number if store variables in DeviceInstance

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
// Function      : Xyce::Device::MemristorPEM::Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
/// Process parameters.
///
/// @return true on success
///
/// In general, the processParams method is intended as a place for complex
/// manipulation of parameters that must happen if temperature or other
/// external changes happen.
///
bool Instance::processParams()
{
  // now set the temperature related stuff.
  double temp=0;
  return updateTemperature(temp);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Register local IDs
///
/// Register the local internal and external node IDs.
///
/// @param intLIDVecRef internal local IDs from topology package
/// @param extLIDVecRef external local IDs from topology package
///
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

  // This lid is for the internal variable for layer thickess that determines resistence
  li_x   = intLIDVec[0];

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Register the local state IDs
///
/// @note The memristor does not have any state vars, so this function
/// does nothing.
///
/// @param staLIDVecRef State variable local IDs
///
/// In general, devices may declare at construction that they require storage
/// locations in the "state vector."  Topology assigns locations in the
/// state vector and returns "Local IDs" (LIDs) for devices to use for their
/// state vector entries.  If a device uses state vector entries, it
/// uses the registerStateLIDs method to save the local IDs for later use.
///
/// @author Robert Hoekstra, SNL, Parallel Computational Sciences
/// @date   06/12/02

void Instance::registerStateLIDs(const std::vector<int> & staLIDVecRef)
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 12/18/2012
//-----------------------------------------------------------------------------
/// Register the local store IDs
///
/// In addition to state vector, Xyce maintains a separate datastructure
/// called a "store" vector.  As with other such vectors, the device
/// declares at construction time how many store vector entries it needs,
/// and later Topology assigns locations for devices, returning LIDs.
///
///
///

void Instance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef)
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());
  li_store_R = stoLIDVecRef[0];
  li_store_tdt = stoLIDVecRef[1];

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Instance::registerBranchDataLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 12/18/2012
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
/// The memristor device uses exactly one "branch data vector" element, where
/// it keeps the "lead current" that may be used on .PRINT lines as
/// "I(ymemristor)" for the current through ymemristor. and a junction voltage.
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
// Creator       : Richard Schiek, Electrical Sysetms Modeling
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
void Instance::loadNodeSymbols(Util::SymbolTable &symbol_table) const
{
  addInternalNode(symbol_table, li_x, getName(), "x");

  addStoreNode(symbol_table, li_store_R, getName(), "R");
  addStoreNode(symbol_table, li_store_tdt, getName(), "TDT");

  if (loadLeadCurrent)
  {
    addBranchDataNode( symbol_table, li_branch_data, getName(), "BRANCH_D");
  }
}




//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Register the Jacobian local IDs
///
/// @param jacLIDVec Jacobian local Ids
///
/// @see Xyce::Device::MemristorPEM::Instance::initializeJacobianStamp
///
/// Having established local IDs for the solution variables, Topology must
/// also assign local IDs for the elements of the Jacobian matrix.
///
/// For each non-zero element that was identified in the jacobianStamp,
/// Topology will assign a Jacobian LID.  The jacLIDVec will have the
/// same structure as the JacStamp, but the values in the jacLIDVec will
/// be offsets into the row of the sparse Jacobian matrix corresponding
/// to the non-zero identified in the stamp.
///
/// These offsets are stored in this method for use later when we load
/// the Jacobian.
///
/// @author Robert Hoekstra, SNL, Parallel Computational Sciences
/// @date   08/27/01

void Instance::registerJacLIDs(const std::vector< std::vector<int> > & jacLIDVec)
{
  // Let DeviceInstance do its work.
  DeviceInstance::registerJacLIDs(jacLIDVec);

  // Store offsets of the components of the Jacobian of this instance
  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNegNodeOffset = jacLIDVec[0][1];
  APosEquXNodeOffset   = jacLIDVec[0][2];
  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];
  ANegEquXNodeOffset   = jacLIDVec[1][2];
  XEquVPosOffset       = jacLIDVec[2][0];
  XEquVNegOffset       = jacLIDVec[2][1];
  XEquXOffset          = jacLIDVec[2][2];

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Setup direct access pointer to solution matrix and vectors.
///
/// @see Xyce::Device::MemristorPEM::Instance::registerJacLIDs
///
/// As an alternative to the row offsets defined in registerJacLIDs, it
/// is also possible to obtain direct pointers of the Jacobian elements.
///
/// This method uses the offsets obtained in registerJacLIDs to retrieve
/// the pointers.
///
/// In this device the pointers to the matrix are only saved
/// (and are only used for matrix access) if
/// Xyce_NONPOINTER_MATRIX_LOAD is NOT defined at compile time.  For
/// some devices the use of pointers instead of array indexing can be
/// a performance enhancement.
///
/// Use of pointers in this device is disabled by defining
/// Xyce_NONPOINTER_MATRIX_LOAD at compile time.  When disabled, array
/// indexing with the offsets from registerJacLIDs is used in
/// the matrix load methods.
///
/// @author Eric Keiter, SNL
/// @date   11/30/08

void Instance::setupPointers()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);
  f_PosEquPosNodePtr = &(dFdx[li_Pos][APosEquPosNodeOffset]);
  f_PosEquNegNodePtr = &(dFdx[li_Pos][APosEquNegNodeOffset]);
  f_PosEquXNodePtr   = &(dFdx[li_Pos][APosEquXNodeOffset]);
  f_NegEquPosNodePtr = &(dFdx[li_Neg][ANegEquPosNodeOffset]);
  f_NegEquNegNodePtr = &(dFdx[li_Neg][ANegEquNegNodeOffset]);
  f_NegEquXNodePtr   = &(dFdx[li_Neg][ANegEquXNodeOffset]);
  f_XEquPosNodePtr   = &(dFdx[li_x][XEquVPosOffset]);
  f_XEquNegNodePtr   = &(dFdx[li_x][XEquVNegOffset]);
  f_XEquXNodePtr     = &(dFdx[li_x][XEquXOffset]);
  q_XEquXNodePtr     = &(dQdx[li_x][XEquXOffset]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Update the intermediate variables
///
/// @return true on success
///
/// The bulk of any device's computation is handled in the instance class'
/// updateIntermediateVars method.  For the memristor, this is
/// merely the computation of the current through the device given the
/// voltage difference between its terminals.
///
/// Intermediate variables computed here are used in the methods that
/// load data into the F, Q, dFdX and dQdX data structures.
///
/// @note This method is called by the updatePrimaryState
/// method. Since the MemristorPEM class reimplements the "Master"
/// "loadState" function that loads the contributions from all
/// such devices in a single loop, THIS FUNCTION IS NOT ACTUALLY
/// USED!
///
///
bool Instance::updateIntermediateVars()
{
  double * solVec = extData.nextSolVectorRawPtr;
  double v_pos    = solVec[li_Pos];
  double v_neg    = solVec[li_Neg];
  double x        = solVec[li_x];

  G=1.0/Reff;
  i0= (v_pos-v_neg)*G;

  {
    // calculate the xVar F equation update and the
    // derivatives d(F(x))/dVpos, d(F(x))/dVneg and d(F(x))/dx
    /*
    Sacado::Fad::SFad<double,3> varVpos( 3, 0, v_pos);
    Sacado::Fad::SFad<double,3> varVneg( 3, 1, v_neg);
    Sacado::Fad::SFad<double,3> varX( 3, 2, x );
    Sacado::Fad::SFad<double,3> parami( i0 );
    Sacado::Fad::SFad<double,3> paramG( G );
    Sacado::Fad::SFad<double,3> paramiOff( model_.iOff_ );
    Sacado::Fad::SFad<double,3> paramiOn( model_.iOn_ );
    Sacado::Fad::SFad<double,3> paramkOff( model_.kOff_ );
    Sacado::Fad::SFad<double,3> paramkOn( model_.kOn_ );
    Sacado::Fad::SFad<double,3> paramAlphaOff( model_.alphaOff_ );
    Sacado::Fad::SFad<double,3> paramAlphaOn( model_.alphaOn_ );
    Sacado::Fad::SFad<double,3> resultFad;

    xVarFterm(varVpos, varVneg, varX, paramG, paramiOff, paramiOn, paramkOff, paramkOn,
      paramAlphaOff, paramAlphaOn, resultFad);

    // calculate the window function value
    Sacado::Fad::SFad<double,3> windowFunctionFad;

    if( model_.windowType_ ==  0)
    {
      // no window
      windowFunctionFad = 1.0;
    }
    else if( model_.windowType_ ==  1)
    {
      Sacado::Fad::SFad<double,3> paramD( model_.D_ );
      Sacado::Fad::SFad<double,3> paramP( model_.p_ );
      JogelkarWindowFunction( varX, paramD, paramP, windowFunctionFad );
    }
    else if( model_.windowType_ ==  2)
    {
      Sacado::Fad::SFad<double,3> paramD( model_.D_ );
      Sacado::Fad::SFad<double,3> paramP( model_.p_ );
      BiolekWindowFunction( varX, paramD, paramP, parami, windowFunctionFad );
    }
    else if( model_.windowType_ ==  3)
    {
      Sacado::Fad::SFad<double,3> paramD( model_.D_ );
      Sacado::Fad::SFad<double,3> paramP( model_.p_ );
      Sacado::Fad::SFad<double,3> paramJ( model_.j_ );
      ProdromakisWindowFunction( varX, paramD, paramP, paramJ, windowFunctionFad );
    }
    else if( model_.windowType_ ==  4)
    {
      Sacado::Fad::SFad<double,3> paramAOn( model_.aOn_ );
      Sacado::Fad::SFad<double,3> paramAOff( model_.aOff_ );
      Sacado::Fad::SFad<double,3> paramWc( model_.wc_ );
      TEAMWindowFunctionF( varX, parami, paramAOff, paramAOn, paramWc, windowFunctionFad );
    }

    xVarFContribution = resultFad.val();
    dxFEqdVpos = resultFad.dx(0);
    dxFEqdVneg = resultFad.dx(1);
    dxFEqdx = resultFad.dx(2);
    //xVarFContribution = 1.0;
    //dxFEqdVpos = 0.0;
    //dxFEqdVneg = 0.0;
    //dxFEqdx = 0.0;
    */
  }

  //i0 = (v_pos-v_neg)*G;

  // get derivative terms for jacobian load
  //
  // APosEquPosNodeOffset  -->  G
  // APosEquNegNodeOffset  --> -G
  // APosEquXNodeOffset    --> -Reff^2 dR/dx
  // ANegEquPosNodeOffset  --> -G
  // ANegEquNegNodeOffset  -->  G
  // ANegEquXNodeOffset    -->  Reff^2 dR/dx
  // XEquVPosOffset        --> d xVarFterm / dVpos
  // XEquVNegOffset        --> d xVarFterm / dVneg
  // XEquXOffset           --> d xVarFterm / dx

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Update the state variables.
///
/// @return true on success
///
/// This function is the function that is called from the device manager
/// each time a load of vectors and Jacobian is required.  Its first
/// job is to call updateIntermediateVars.
///
/// After calling updateIntermediateVars, the updatePrimaryState method
/// may load state vector elements as needed.
///
/// The memristor device has no state vector elements, so all this method does
/// in the memristor instance class is call updateIntermediateVars.
///
/// There is no longer a "secondary" state.  The "primary" in
/// this method's name is purely historical.
///
/// @note This method is called by the default implementation of the
/// loadState master function. Since the MemristorPEM class reimplements
/// the "Master" "loadState" function that loads the contributions
/// from all memristor devices in a single loop, THIS FUNCTION IS NOT
/// ACTUALLY USED, NOR is the updateIntermediateVars method it calls!
/// The updatePrimaryState method is only called when a device class
/// does not re-implement the master class.  This can be a source of
/// confusion when attempting to modify the MemristorPEM device, or any
/// other device That reimplements the Master classes.
///
/// @see Xyce::Device::MemristorPEM::Master::updateState
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   01/29/01
///
bool Instance::updatePrimaryState()
{
  return updateIntermediateVars();
}

bool Instance::loadDAEQVector()
{
  double * qVec = extData.daeQVectorRawPtr;
  double * solVec = extData.nextSolVectorRawPtr;
  double x     = solVec[li_x];

  qVec[li_x]   += x;

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Instance::loadDAEFVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Load the DAE F vector.
///
/// @return true on success
///
/// The Xyce DAE formulation solves the differential-algebraic
/// equations \f$F(x)+dQ(x)/dt=0\f$ These are vector equations
/// resulting from the modified nodal analysis of the circuit.
///
/// This method loads the F-vector contributions for a single memristor
/// instance.
///
/// In this method, the offsets defined in registerLIDs are used to
/// store the device's F contributions into the F vector.
///
/// The Q vector is used for devices that store charge or magnetic
/// energy, and since the memristor does none of that, the F vector
/// contributions are the whole of the memristor's contribution to the
/// full set of DAEs.
///
/// @note This method is called by the default implementation of the
/// loadDAEVectors master function. Since the MemristorPEM class
/// reimplements the "Master" "loadDAEVectors" function that loads the
/// contributions from all memristor devices in a single loop, THIS
/// FUNCTION IS NOT ACTUALLY USED.  The loadDAEFVector method is only
/// called when a device class does not re-implement the master class.
/// This can be a source of confusion when attempting to modify the MemristorPEM
/// device, or any other device that reimplements the Master classes.
///
/// @see Xyce::Device::MemristorPEM::Master::loadDAEVectors
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   01/24/03
///
bool Instance::loadDAEFVector()
{
  double * fVec = extData.daeFVectorRawPtr;
  fVec[li_Pos] += i0;
  fVec[li_Neg] += -i0;
  fVec[li_x]   += xVarFContribution;

  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    double * junctionV = extData.nextJunctionVCompRawPtr;
    double * solVec = extData.nextSolVectorRawPtr;
    leadF[li_branch_data] = i0;
    junctionV[li_branch_data] = solVec[li_Pos] - solVec[li_Neg];
  }


  return true;
}

bool Instance::loadDAEdQdx()
{
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);
  dQdx[li_x][XEquXOffset] = 1.0;
  return true;
}


//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Instance::loadDAEdFdx
// Purpose       : Loads the F-vector contributions for a single
//                 memristor  instance.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Load the DAE the derivative of the F vector with respect to the
/// solution vector x, dFdx
///
/// Loads the contributions for a single memristor instance to the
/// dFdx matrix (the F contribution to the Jacobian).
///
/// This method uses the Jacobian LIDs (row offsets) that were stored by
/// registerJacLIDs.
///
/// @see Xyce::Device::MemristorPEM::Instance::registerJacLIDs
///
/// @note This method is called by the default implementation of the
/// loadDAEMatrices master function. Since the MemristorPEM class
/// reimplements the "Master" "loadDAEMatrices" function that loads the
/// contributions from all memristor devices in a single loop, THIS
/// FUNCTION IS NOT ACTUALLY USED.  The loadDAEdFdx method is only
/// called when a device class does not re-implement the master class.
/// This can be a source of confusion when attempting to modify the MemristorPEM
/// device, or any other device that reimplements the Master classes.
///
/// @see Xyce::Device::MemristorPEM::Master::loadDAEMatrices
///
/// @return true on success
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   03/05/04
bool Instance::loadDAEdFdx()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  dFdx[li_Pos][APosEquPosNodeOffset] += G;
  dFdx[li_Pos][APosEquNegNodeOffset] -= G;
  dFdx[li_Pos][APosEquXNodeOffset]   += dReffdvpos;
  dFdx[li_Neg][ANegEquPosNodeOffset] -= G;
  dFdx[li_Neg][ANegEquNegNodeOffset] += G;
  dFdx[li_Neg][ANegEquXNodeOffset]   += dReffdvneg;
  dFdx[li_x][XEquVPosOffset]         += dxFEqdVpos;
  dFdx[li_x][XEquVNegOffset]         += dxFEqdVneg;
  dFdx[li_x][XEquXOffset]            += dxFEqdx;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Update the parameters that depend on the temperature of the device
///
/// @param temp_tmp temperature
///
/// Xyce has a number of mechanisms that allow temperature to be changed
/// after a device has been instantiated.  These include .STEP loops over
/// temperature.  When temperature is changed, any device that has parameters
/// that depend on temperature must be updated.  That updating happens here.
///
/// The MemristorPEM device supports temperature-dependent resistance through its
/// TC1 (linear dependence) and TC2 (quadratic dependence) parameters.
/// If these parameters are specified, the resistance must be updated.
///
/// @return true on success
///
/// @author Tom Russo, Component Information and Models
/// @date   02/27/01
bool Instance::updateTemperature(const double & temp_tmp)
{
  bool bsuccess = true;
  /*
  double difference, factor;

  if (temp_tmp != -999.0)
    temp = temp_tmp;
  difference = temp - model_.tnom;
  factor = model_.resistanceMultiplier*(1.0 + tempCoeff1*difference + tempCoeff2*difference*difference);

  if (R*factor != 0.0)
    G = 1.0/(R * factor);
  else
    G = 0.0;
  */
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Model::processParams
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
// Function      : Xyce::Device::MemristorPEM::Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//----------------------------------------------------------------------------
///
/// Process the instance parameters of instance owned by this model
///
/// This method simply loops over all instances associated with this
/// model and calls their processParams method.
///
/// @return true
///
/// @author Dave Shirely, PSSI
/// @date   03/23/06

bool Model::processInstanceParams()
{
  for (InstanceVector::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
  {
    (*it)->processParams();
  }

  return true;
}
//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Construct a memristor model from a "model block" that was created
/// by the netlist parser.
///
/// @param configuration
/// @param model_block
/// @param factory_block
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   5/16/00
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    model_block,
  const FactoryBlock &  factory_block)
  : DeviceModel(model_block, configuration.getModelParameters(), factory_block),
  V1_(0.0),
  V2_(0.0),
  I1_(0.0),
  I2_(0.0),
  G0_(0.0),
  Vp_(0.0),
  Vn_(0.0),
  D1_(0.0),
  D2_(0.0)
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

  // try to load the data arrays from the supplied file
  std::ifstream inputFile;
  inputFile.open( fxPDataFileName_.c_str(), std::ifstream::in );

  if( !inputFile.good() )
  {
    //couldn't open input file.  Signal fatal error.
    Report::UserFatal() << "Couldn't open data file \"" << fxPDataFileName_ << "\" for MemristorPEM device model " << getName() << std::endl;
  }
  const size_t maxLineLength=1024;
  char aLine[maxLineLength];
  char commentChar='#';
  while( inputFile.good() )
  {
    char c = inputFile.peek();
    if( c == commentChar )
    {
      // comment line so skip it
      inputFile.getline( aLine, maxLineLength );
    }
    else
    {
      // have a good line.  try to break it into two numbers
      double num1, num2;
      char separator;
      if( (inputFile >> num1 ) && (inputFile >> separator) && (inputFile >> num2 ) )
      {
        xpData.push_back(num1);
        fpData.push_back(num2);
      }
    }
  }

  Xyce::dout() << "MemristorPEM model " << getName() << " read " << xpData.size() << ", " << fpData.size() << " data points from " << fxPDataFileName_ << std::endl;
  // set up interpolator
  FpInterp.init( xpData, fpData );
  inputFile.close();

  inputFile.open( fxMDataFileName_.c_str(), std::ifstream::in );

  if( !inputFile.good() )
  {
    //couldn't open input file.  Signal fatal error.
    Report::UserFatal() << "Couldn't open data file \"" << fxMDataFileName_ << "\" for MemristorPEM device model " << getName() << std::endl;
  }
  while( inputFile.good() )
  {
    char c = inputFile.peek();
    if( c == commentChar )
    {
      // comment line so skip it
      inputFile.getline( aLine, maxLineLength );
    }
    else
    {
      // have a good line.  try to break it into two numbers
      double num1, num2;
      char separator;
      if( (inputFile >> num1 ) && (inputFile >> separator) && (inputFile >> num2 ) )
      {
        xmData.push_back(num1);
        fmData.push_back(num2);
      }
    }
  }

  Xyce::dout() << "MemristorPEM model " << getName() << " read " << xmData.size() << ", " << fmData.size() << " data points from " << fxMDataFileName_ << std::endl;
  xMax = xpData[ xpData.size() -1 ];
  xMin = xmData[ 0 ];
  // set up interpolator
  FmInterp.init( xmData, fmData );
  inputFile.close();

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Destroy this model.
///
/// Also destroys all instances that use this model.
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   3/16/00
Model::~Model()
{
  // Destory all owned instances
  for (InstanceVector::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
  {
    delete (*it);
  }

}

//-----------------------------------------------------------------------------
// Function      : N_DEV_MemristorPEMModel::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Print instances associated with this model.
///
/// Used only for debugging
///
/// @param os output stream
///
/// @return reference to output stream
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   4/03/00
///
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  os << std::endl;
  os << "Number of MemristorPEM Instances: " << instanceContainer.size() << std::endl;
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
/// Apply a device instance "op" to all instances associated with this
/// model
///
/// @param[in] op Operator to apply to all instances.
///
///
void Model::forEachInstance(DeviceInstanceOp &op) const
{
  for (std::vector<Instance *>::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
    op(*it);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Update state for all memristor instances, regardless of model.
///
/// @param solVec solution vector
/// @param staVec state vector
/// @param stoVec store vector
///
/// @return true on success
///
/// @note Because the memristor device re-implements the base-class
/// Master::updateState, the Instance::updatePrimaryState method is never
/// called, nor is the Instance::updateIntermediateVars method.  This method
/// replaces those, and does the same work but inside a loop over all
/// memristor instances.
///
/// @see Xyce::Device::MemristorPEM::Instance::updatePrimaryState
/// @author Eric Keiter, SNL
/// @date   11/26/08

bool Master::updateState(double * solVec, double * staVec, double * stoVec)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & ri = *(*it);

    // For this model I-V relationship is
    //
    // i(t) = x h(V)
    //
    // with h(V) = I1 exp( (v1-v2) / V1_ ) - I2_ exp( (v1-v2)/V2_ ) + G0 (v1-v2) - (I1_ - I2)
    //

    double v_pos    = solVec[ri.li_Pos];
    double v_neg    = solVec[ri.li_Neg];
    double x        = solVec[ri.li_x];


    {
      // extra scope to limit variable extent.  Otherwise one will have lots of
      // similar variables
      //Sacado::Fad::SFad<double,2> varDeltaV( 2, 0, (v_pos-v_neg) );
      Sacado::Fad::SFad<double,3> varV1( 3, 0, v_pos );
      Sacado::Fad::SFad<double,3> varV2( 3, 1, v_neg );
      Sacado::Fad::SFad<double,3> varX( 3, 2, x );

      Sacado::Fad::SFad<double,3> resultFad;
      I_V_Fxn( varV1, varV2, varX, ri.model_.V1_, ri.model_.V2_, ri.model_.G0_, ri.model_.I1_, ri.model_.I2_, resultFad );

      ri.i0 = resultFad.val();
      ri.G  = resultFad.dx(0);
      ri.Reff = 1.0/ri.G;
      ri.dIdx = resultFad.dx(2);
    }


    // state variable equation is
    //
    // dx/dt = c f(x) g(V)
    //
    // in DAE formulation is:
    //
    // -1 dx/dt + c f(x) g(V1,V2)
    //
    // f(x) is the piece-wise continuous function.
    //

    if( getSolverState().dcopFlag )
    {
      // during the dc op the internal variable should
      // follow the equation
      // x - x0 == 0
      // where x0 is the initial conditon if it's given
      // or x0 = 1 if Vp > Vn
      //    x0 = 0 if Vp <= Vn
      double xInf = ri.XO_;
      if(!ri.XOGiven_ && (v_pos > v_neg ) )
      {
        xInf = 1.0;
      }
      else if( !ri.XOGiven_ && (v_pos <= v_neg ) )
      {
        xInf = 0.0;
      }
      ri.xVarFContribution = x - xInf;
      ri.dxFEqdVpos = 0.0;
      ri.dxFEqdVneg = 0.0;
      ri.dxFEqdx = 1.0 ;

    }
    else
    {
      // evaluate the state variable equation
      //Sacado::Fad::SFad<double,2> varDeltaV( 2, 0, (v_pos-v_neg) );
      Sacado::Fad::SFad<double,3> varV1( 3, 0, v_pos );
      Sacado::Fad::SFad<double,3> varV2( 3, 1, v_neg );
      Sacado::Fad::SFad<double,3> varX( 3, 2, x );

      Sacado::Fad::SFad<double,3> g_function_resultFad;

      G( varV1, varV2, ri.model_.D1_, ri.model_.D2_, ri.model_.Vp_, ri.model_.Vn_, g_function_resultFad) ;

      // for now I've factored out the user supplied Fp(x) and Fm(x) so that they could
      // be PWL functions.  Interpolate the value of F_(x) and it's derivative
      double fx=0;
      double dfdx=0;
      if( v_pos >= v_neg )
      {
        ri.model_.FpInterp.eval (ri.model_.xpData, ri.model_.fpData, x, fx);
        ri.model_.FpInterp.evalDeriv (ri.model_.xpData, ri.model_.fpData, x, dfdx);
        if( x > ri.model_.xMax )
        {
          fx = 0.0;
          x = ri.model_.xMax;
          dfdx = 0.0;
        }

        ri.xVarFContribution = ri.model_.C1_*g_function_resultFad.val()*fx;

        ri.dxFEqdVpos = ri.model_.C1_*g_function_resultFad.dx(0)*fx;
        ri.dxFEqdVneg = ri.model_.C1_*g_function_resultFad.dx(1)*fx;
        ri.dxFEqdx = ri.model_.C1_*g_function_resultFad.val()*dfdx;

      }
      else
      {
        ri.model_.FmInterp.eval (ri.model_.xmData, ri.model_.fmData, x, fx);
        ri.model_.FmInterp.evalDeriv (ri.model_.xmData, ri.model_.fmData, x, dfdx);
        if( x < ri.model_.xMin )
        {
          fx = 0.0;
          x = ri.model_.xMin;
          dfdx = 0.0;
        }

        ri.xVarFContribution = ri.model_.C2_*g_function_resultFad.val()*fx;

        ri.dxFEqdVpos = ri.model_.C2_*g_function_resultFad.dx(0)*fx;
        ri.dxFEqdVneg = ri.model_.C2_*g_function_resultFad.dx(1)*fx;
        ri.dxFEqdx = ri.model_.C2_*g_function_resultFad.val()*dfdx;
      }
      //Xyce::dout() << " VP=" << v_pos << " VN=" << v_neg << "  x=" << x << " fx=" << fx << " dfdx=" << dfdx << " ri.xVarFContribution=" << ri.xVarFContribution << " ri.dxFEqdVneg=" << ri.dxFEqdVneg << " ri.dxFEqdx=" << ri.dxFEqdx << std::endl;


    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Load DAE vectors of all memristor instances, regardless of model
///
/// @param solVec solution vector
/// @param fVec f vector
/// @param qVec q vector
/// @param leadF store lead current f vector
/// @param leadQ store lead current q vector
///
/// @return true on success
///
/// @note Because the memristor device re-implements the base-class
/// Master::loadDAEVectors, the Instance::loadDAEFVector method is
/// never called.  This method replaces those, and does the same work
/// but inside a loop over all memristor instances.
///
/// @see Xyce::Device::MemristorPEM::Instance::loadDAEFVector
///
/// @author Eric Keiter, SNL
/// @date   11/26/08

bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, double * leadF, double * leadQ, double * junctionV)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & ri = *(*it);
    fVec[ri.li_Pos] += ri.i0;
    fVec[ri.li_Neg] += -ri.i0;
    fVec[ri.li_x]   += ri.xVarFContribution;
    qVec[ri.li_x]   -= solVec[ri.li_x];
    /*
    if( getSolverState().dcopFlag )
    {
      qVec[ri.li_x] -= ri.XO_;
    }
    */
    if( ri.G != 0 )
    {
      double * storeVec = ri.extData.nextStoVectorRawPtr;
      storeVec[ri.li_store_R] = 1.0/ri.G;
    }

    if( ri.loadLeadCurrent )
    {
      leadF[ri.li_branch_data] = ri.i0;
      junctionV[ri.li_branch_data] = solVec[ri.li_Pos] - solVec[ri.li_Neg];
    }

  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Load DAE matrices for all memristor instances, regardless of model
///
/// @param dFdx matrix of derivatives of F vector with respect to solution
/// @param dQdx matrix of derivatives of Q vector with respect to solution
///
/// @return true on success
///
/// @note Because the memristor device re-implements the base-class
/// Master::loadDAEMatrices, the Instance::loadDAEdFdx method is
/// never called.  This method replaces those, and does the same work
/// but inside a loop over all memristor instances.
///
/// @see Xyce::Device::MemristorPEM::Instance::loadDAEdFdx
///
/// @author Eric Keiter, SNL
/// @date   11/26/08

bool Master::loadDAEMatrices(Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & ri = *(*it);

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    *(ri.f_PosEquPosNodePtr) += ri.G;
    *(ri.f_PosEquNegNodePtr) -= ri.G;
    *(ri.f_PosEquXNodePtr)   += ri.dIdx;
    *(ri.f_NegEquPosNodePtr) -= ri.G;
    *(ri.f_NegEquNegNodePtr) += ri.G;
    *(ri.f_NegEquXNodePtr )  += ri.dIdx;
    *(ri.f_XEquPosNodePtr )  += ri.dxFEqdVpos;
    *(ri.f_XEquNegNodePtr )  += ri.dxFEqdVneg;
    *(ri.f_XEquXNodePtr )    += ri.dxFEqdx;
    *(ri.q_XEquXNodePtr )    = -1.0;
#else
    dFdx[ri.li_Pos][ri.APosEquPosNodeOffset] += ri.G;
    dFdx[ri.li_Pos][ri.APosEquNegNodeOffset] -= ri.G;
    dFdx[ri.li_Pos][ri.APosEquXNodeOffset]   += ri.dIdx;
    dFdx[ri.li_Neg][ri.ANegEquPosNodeOffset] -= ri.G;
    dFdx[ri.li_Neg][ri.ANegEquNegNodeOffset] += ri.G;
    dFdx[ri.li_Neg][ri.ANegEquXNodeOffset]   += ri.dIdx;
    dFdx[ri.li_x][ri.XEquVPosOffset]         += ri.dxFEqdVpos;
    dFdx[ri.li_x][ri.XEquVNegOffset]         += ri.dxFEqdVneg;
    dFdx[ri.li_x][ri.XEquXOffset]            += ri.dxFEqdx;
    dQdx[ri.li_x][ri.XEquXOffset] = -1.0;
#endif
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::Traits::factory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Create a new instance of the MemristorPEM device.
///
/// @param configuration
/// @param factory_block
///
Device *
Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{
  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorPEM::registerDevice
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
///
/// Define how to use the device in a netlist.
///
/// This method is called from the Xyce::Device::registerOpenDevices
/// function, which in turn is called by the device manager.
///
/// The device is declared here to be an "memristor" device, which must
/// take a model card of type "memristor".  This device will correspond to model
/// level 3 of memristor models.
void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() ||
      ((deviceMap.find("MEMRISTOR")!=deviceMap.end()) && (levelSet.find(4)!=levelSet.end())))
  {
    MemristorTEAM::registerDevice();

    Config<Traits>::addConfiguration()
      .registerDevice("memristor", 4)
      .registerModelType("memristor", 4);
  }
}

//-----------------------------------------------------------------------------
// Function      : memristorPEMSensitivity::operator
// Purpose       : produces df/dp and dq/dp, where p=R.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
void memristorPEMSensitivity::operator()(
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

} // namespace MemristorPEM
} // namespace Device
} // namespace Xyce

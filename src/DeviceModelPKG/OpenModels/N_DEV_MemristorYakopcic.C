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

//----------------------------------------------------------------------------
//
// Purpose        : Implementation of the Yakopcic memristor model.  See.
//                  Yakopcic: Threshold Adaptive Memristor Model
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

#include <N_DEV_MemristorYakopcic.h>

#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>
#include <Sacado.hpp>

namespace Xyce {
namespace Device {
namespace MemristorYakopcic {

// 
// Threshold voltage function: G(v, Ap, An, Vp, Vn)
//
template <typename ScalarT>
void G( const ScalarT & V1, const ScalarT & V2, const ScalarT & Ap, const ScalarT & An, ScalarT & Vp, ScalarT & Vn, ScalarT & fval )
{
  if( (V1-V2) <= Vp )
  {
    if( (V1-V2) >= -Vn )
    {
      fval = 0.0;
    }
    else
    {
      fval = -An * (exp(-(V1-V2)) - exp(Vn) );
    }
  }
  else
  {
    fval = Ap * (exp((V1-V2)) - exp(Vp) );
  }
}

//
// Equations for the state variable motion
// 

template <typename ScalarT>
void WP( const ScalarT & X, const ScalarT & Xp, ScalarT & fval )
{
  // equation 5
  fval = 1 + (Xp - X) / (1.0 - Xp);
}

template <typename ScalarT>
void WN( const ScalarT & X, const ScalarT & Xn, ScalarT & fval )
{
  // equation 6
  fval = X / (1.0 - Xn);
}


template <typename ScalarT>
void F_x_equation( const ScalarT & V1, const ScalarT & V2, const ScalarT & X, const ScalarT & Xp, const ScalarT & Xn, 
  const ScalarT & AlphaP, const ScalarT & AlphaN, const ScalarT & Eta, ScalarT & fval )
{
  if( Eta*(V1-V2) > 0 )
  {
    if(X > Xp) 
    {
      ScalarT WPval;
      WP( X, Xp, WPval );
      fval = WPval * exp(-AlphaP * (X - Xp ) );  // equation 3
    }
    else
    {
      fval = 1.0;
    }
  }
  else
  {
    if(X <= (1-Xn)) 
    {
      ScalarT WNval;
      WN( X, Xn, WNval );
      fval = WNval * exp(AlphaN * (X + Xn -1) );  // equation 4
    }
    else
    {
      fval = 1.0;
    }  
  }
}


template <typename ScalarT>
void X_var_F_equation( const ScalarT & V1, const ScalarT & V2, const ScalarT & X,
  const ScalarT & Ap, const ScalarT & An, ScalarT & Vp, ScalarT & Vn,
  const ScalarT & Xp, const ScalarT & Xn,  const ScalarT & AlphaP, const ScalarT & AlphaN, 
  const ScalarT & Eta, ScalarT & fval )
{
  ScalarT Gval;
  G( V1, V2, Ap, An, Vp, Vn, Gval );
  
  ScalarT FXval;
  F_x_equation( V1, V2, X, Xp, Xn, AlphaP, AlphaN, Eta, FXval );
  fval = Eta * Gval * FXval;
}


 // For this model I-V relationship is 
 //
 // i(t) = a1 x(t) sinh( b V(t) )  : iff V(t) >= 0
 //        a2 x(t) sinh( b V(t) )  : iff V(t) < 0 
 //
 // so need to update the x(t) equation and 
 // calculate I(t) based on V(t)
 //
 
template <typename ScalarT>
void I_V_Fxn( const ScalarT & V1, const ScalarT & V2, const ScalarT & X, const ScalarT & A1, const ScalarT & A2, const ScalarT & B, ScalarT & fval )
{
  // equation 1
  if( (V1-V2) >= 0.0 )
  {
    fval = A1 * X * sinh( B * (V1-V2) );
  }
  else
  {
    fval = A2 * X * sinh( B * (V1-V2) );
  }
  
}




// Common Jacobian Stamp for all MemristorYakopcic devices.
// Because all memristors have identical Jacobian stamps, this data is
// declared static and is shared by all memristor instances.
// 
std::vector<std::vector<int> > Instance::jacStamp;

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorYakopcic::Instance::initializeJacobianStamp
// Purpose       : 
// Special Notes : 
// Scope         : private
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 2/11/2014
//-----------------------------------------------------------------------------
//
// @brief Common Jacobian stamp initializer for all MemristorYakopcic devices.
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
// For this memristor, there are two external nodes (the positive and negative
// terminals of the device).  The positive node is referred to as the 0th
// node of the device, and the negative node the 1st node.
//
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
// Function      : Xyce::Device::MemristorYakopcic::Traits::loadInstanceParameters
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
// Each parameter supported by the memristor device instance is
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
void Traits::loadInstanceParameters(ParametricData<MemristorYakopcic::Instance> &p)
{
  p.addPar("XO", 0.0, &MemristorYakopcic::Instance::XO_)
    .setUnit(U_NONE)
    .setDescription("Initial value for internal variable x");

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorYakopcic::Traits::loadModelParameters
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
// @see Xyce::Device::MemristorYakopcic::Traits::loadInstanceParameters
//
//
void Traits::loadModelParameters(ParametricData<MemristorYakopcic::Model> &p)
{

  // Create parameter definitions for parameter member variables
  p.addPar("ETA", 1.0, &MemristorYakopcic::Model::Eta_)
    .setUnit(U_NONE)
    .setDescription("State variable motion relative to voltage.");
  p.addPar("VP", 1.0e-2, &MemristorYakopcic::Model::Vp_)
    .setUnit(U_VOLT)
    .setDescription("Positive Voltage Threshold");
  p.addPar("VN", -1.0e-2, &MemristorYakopcic::Model::Vn_)
    .setUnit(U_VOLT)
    .setDescription("Negative Voltage Threshold");  
  p.addPar("AP", 1.0, &MemristorYakopcic::Model::Ap_)
    .setUnit(U_NONE)
    .setDescription("Positive Voltage Threshold Magnitude Parameter"); 
  p.addPar("AN", 1.0, &MemristorYakopcic::Model::An_)
    .setUnit(U_NONE)
    .setDescription("Negative Voltage Threshold Magnitude Parameter");
  p.addPar("A1", 1.0, &MemristorYakopcic::Model::A1_)
    .setUnit(U_NONE)
    .setDescription("Dielectric layer thickness parameter [dimensionless] "); 
  p.addPar("A2", 1.0, &MemristorYakopcic::Model::A2_)
    .setUnit(U_NONE)
    .setDescription("Dielectric layer thickness parameter [dimensionless] ");
  p.addPar("B", 1.0, &MemristorYakopcic::Model::B_)
    .setUnit(U_NONE)
    .setDescription("Curvature in I-V relation.  Relates to how much conduction in the device is Ohmic and versus tunnel barrier.");
  p.addPar("ALPHAP", 1.0, &MemristorYakopcic::Model::AlphaP_)
    .setUnit(U_NONE)
    .setDescription("State variable motion.");
  p.addPar("ALPHAN", 1.0, &MemristorYakopcic::Model::AlphaN_)
    .setUnit(U_NONE)
    .setDescription("State variable motion.");
  p.addPar("XP", 1.0, &MemristorYakopcic::Model::XP_)
    .setUnit(U_NONE)
    .setDescription("State variable motion.");
  p.addPar("XN", 1.0, &MemristorYakopcic::Model::XN_)
    .setUnit(U_NONE)
    .setDescription("State variable motion.");

  p.addPar("XSCALING", 1.0, &MemristorYakopcic::Model::xScaling_)
    .setUnit(U_NONE)
    .setDescription("Scaling for x variable.  For example 1e9 if x will be in units of nanometers.");

  p.addPar("RESNOISE", false,&MemristorYakopcic::Model::randomResNoiseOn_)
    .setUnit(U_NONE)
    .setDescription("RTN model in resistance (on/off)" );

  p.addPar("RESSEED", 0,&MemristorYakopcic::Model::randomResNoiseSeed_)
    .setUnit(U_NONE)
    .setDescription("RTN model in resistance: seed" );

  p.addPar("RESLAMBDA", 0.0,&MemristorYakopcic::Model::randomResNoiseLambda_)
    .setUnit(U_NONE)
    .setDescription("RTN model: lambda" );

  p.addPar("RESTD", 0.0,&MemristorYakopcic::Model::randomResUpdateTime_)
    .setUnit(U_SECOND)
    .setDescription("RTN model in resistance: Update time" );
    
  p.addPar("RESEPTD", 1.0e-10,&MemristorYakopcic::Model::randomResEpsilonUpdateTime_)
    .setUnit(U_SECOND)
    .setDescription("RTN model in resistance: Minimum allowed update time" );
    
  p.addPar("RESDELTA", 0.0,&MemristorYakopcic::Model::randomResDelta_)
    .setUnit(U_OHM)
    .setDescription("RTN model in resistance: Base change in resistance for RTN" );  
    
  p.addPar("RESDELTAGRAD", 0.0,&MemristorYakopcic::Model::randomResDeltaGrad_)
    .setUnit(U_NONE)
    .setDescription("RTN model in resistance: Base change in resistance for RTN scaled by R" );

  p.addPar("XNOISE", false,&MemristorYakopcic::Model::randomXNoiseOn_)
    .setUnit(U_NONE)
    .setDescription("RTN model in growth (on/off)" );

  p.addPar("XSEED", 0,&MemristorYakopcic::Model::randomXNoiseSeed_)
    .setUnit(U_NONE)
    .setDescription("RTN model in growth: seed" );

  p.addPar("XLAMBDA", 0.0,&MemristorYakopcic::Model::randomXNoiseLambda_)
    .setUnit(U_NONE)
    .setDescription("RTN growth model: lambda" );

  p.addPar("XTD", 0.0,&MemristorYakopcic::Model::randomXUpdateTime_)
    .setUnit(U_SECOND)
    .setDescription("RTN model in growth: Update time" );
    
  p.addPar("XEPTD", 1.0e-10,&MemristorYakopcic::Model::randomXEpsilonUpdateTime_)
    .setUnit(U_SECOND)
    .setDescription("RTN model in growth: Minimum allowed update time" );
    
  p.addPar("XDELTA", 0.0,&MemristorYakopcic::Model::randomXDelta_)
    .setUnit(U_OHM)
    .setDescription("RTN model in growth: Base change in growth rate for RTN" );  
    
  p.addPar("XDELTAGRAD", 0.0,&MemristorYakopcic::Model::randomXDeltaGrad_)
    .setUnit(U_NONE)
    .setDescription("RTN model in growth: Base change in growth for RTN scaled by X" );

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorYakopcic::Instance::Instance
// Purpose       : Instance constructor
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Construct a memristor instance.
//

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
    resNoiseLastUpdateStep(-1),
    resNoiseLastUpdateTime(0.0),
    resNoiseNextUpdateTime(0.0),
    resNoiseHiLoState(0),
    resNoiseRfactor(1.0),
    xNoiseLastUpdateStep(-1),
    xNoiseLastUpdateTime(0.0),
    xNoiseNextUpdateTime(0.0),
    xNoiseHiLoState(0),
    xNoiseFactor(1.0),
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
  
  // if the random resistance noise model is on the update the "update" time.
  if( model_.randomResNoiseOn_ )
  {
    // get a poisson distributed number with mean lambda
    //int aNum = model_.randomNumberGen_->poissonRandom(model_.randomResNoiseLambda_);
    // the actual update interval is the poisson distributed number times the base update time.
    //resNoiseNextUpdateTime = aNum * model_.randomResUpdateTime_ + model_.randomResEpsilonUpdateTime_;
    resNoiseNextUpdateTime = -log(model_.randomNumberGen_->uniformRandom() ) * model_.randomResNoiseLambda_*model_.randomResUpdateTime_;
    // don't allow a zero time for the next update.
    //if( aNum == 0 ) 
    //  resNoiseNextUpdateTime = model_.randomResEpsilonUpdateTime_;
  }


}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorYakopcic::Instance::processParams
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
// Function      : Xyce::Device::MemristorYakopcic::Instance::registerLIDs
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

  // This lid is for the internal variable for layer thickess that determines resistence
  li_x   = intLIDVec[0];

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorYakopcic::Instance::registerStateLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Register the local state IDs
//
// @note The memristor does not have any state vars, so this function
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
// Function      : Xyce::Device::MemristorYakopcic::Instance::registerStoreLIDs
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
  li_store_R = stoLIDVecRef[0];
  li_store_tdt = stoLIDVecRef[1];

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorYakopcic::Instance::registerBranchDataLIDs
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
// The memristor device uses exactly one "branch data vector" element, where
// it keeps the "lead current" that may be used on .PRINT lines as
// "I(ymemristor)" for the current through ymemristor. and a junction voltage.
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
  addInternalNode(symbol_table, li_x, getName(), "x");

  addStoreNode(symbol_table, li_store_R, getName(), "R");
  addStoreNode(symbol_table, li_store_tdt, getName(), "TDT");

  if (loadLeadCurrent)
  {
    addBranchDataNode( symbol_table, li_branch_data, getName(), "BRANCH_D");
  }
}




//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorYakopcic::Instance::registerJacLIDs
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
// @see Xyce::Device::MemristorYakopcic::Instance::initializeJacobianStamp
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
  APosEquXNodeOffset   = jacLIDVec[0][2];
  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];
  ANegEquXNodeOffset   = jacLIDVec[1][2];
  XEquVPosOffset       = jacLIDVec[2][0]; 
  XEquVNegOffset       = jacLIDVec[2][1];
  XEquXOffset          = jacLIDVec[2][2];

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorYakopcic::Instance::setupPointers
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Setup direct access pointer to solution matrix and vectors.
//
// @see Xyce::Device::MemristorYakopcic::Instance::registerJacLIDs
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
// Function      : Xyce::Device::MemristorYakopcic::Instance::updateIntermediateVars
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
// updateIntermediateVars method.  For the memristor, this is
// merely the computation of the current through the device given the
// voltage difference between its terminals.
//
// Intermediate variables computed here are used in the methods that
// load data into the F, Q, dFdX and dQdX data structures.
//
// @note This method is called by the updatePrimaryState
// method. Since the MemristorYakopcic class reimplements the "Master"
// "loadState" function that loads the contributions from all
// such devices in a single loop, THIS FUNCTION IS NOT ACTUALLY
// USED!  
//
//
bool Instance::updateIntermediateVars()
{
  double * solVec = extData.nextSolVectorRawPtr;
  double v_pos    = solVec[li_Pos];
  double v_neg    = solVec[li_Neg];
  double x        = solVec[li_x];
  //double Reff     = model_.ROn_;
  //
  // Calclate Reff and deriative term dReff/dx 
  //
  {
    // extra scope to limit variable extent.  Otherwise one will have lots of 
    // similar variables 
    /*
    Sacado::Fad::SFad<double,1> varX( 1, 0, x );
    Sacado::Fad::SFad<double,1> paramRon( model_.ROn_ );
    Sacado::Fad::SFad<double,1> paramRoff( model_.ROff_ );
    Sacado::Fad::SFad<double,1> paramXon( model_.xOn_ );
    Sacado::Fad::SFad<double,1> paramXoff( model_.xOff_ );
    Sacado::Fad::SFad<double,1> resultFad;
    if(IVrelation==0)
    {
      // liner current voltage relationship
      ReffLin( varX, paramRon, paramRoff, paramXon, paramXoff, resultFad );

    }
    else if( IVrelation==1 )
    {
      ReffNonLin( varX, paramRon, paramRoff, paramXon, paramXoff, resultFad );
    }
    Reff = resultFad.val();
    dReffdx = resultFad.dx(0);
    dReffdvpos = -v_pos * pow( Reff, 2 ) * dReffdx; 
    dReffdvneg =  v_neg * pow( Reff, 2 ) * dReffdx; 
    */
  }

  // if the random noise in the resistance model is on
  // and this is the start of a new step, then get a new 
  // gaussian distributed variable
  double rfactor=1.0;
  if( model_.randomResNoiseOn_ && (getSolverState().timeStepNumber_ != resNoiseLastUpdateStep) )
  {
    /*
    resNoiseLastUpdateStep = getSolverState().timeStepNumber_; 
    rfactor = model_.randomNumberGen_->gaussianRandom(model_.randomResNoiseMean_, model_.randomResNoiseSD_);
    */
  }

  G=1.0/(rfactor*Reff);
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
// Function      : Xyce::Device::MemristorYakopcic::Instance::updatePrimaryState
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
// The memristor device has no state vector elements, so all this method does
// in the memristor instance class is call updateIntermediateVars.
//
// There is no longer a "secondary" state.  The "primary" in
// this method's name is purely historical.
//
// @note This method is called by the default implementation of the
// loadState master function. Since the MemristorYakopcic class reimplements
// the "Master" "loadState" function that loads the contributions
// from all memristor devices in a single loop, THIS FUNCTION IS NOT
// ACTUALLY USED, NOR is the updateIntermediateVars method it calls!
// The updatePrimaryState method is only called when a device class
// does not re-implement the master class.  This can be a source of
// confusion when attempting to modify the MemristorYakopcic device, or any
// other device That reimplements the Master classes.
//
// @see Xyce::Device::MemristorYakopcic::Master::updateState
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
  double x     = solVec[li_x];
  
  qVec[li_x]   += x;

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorYakopcic::Instance::loadDAEFVector
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
// This method loads the F-vector contributions for a single memristor
// instance.
//
// In this method, the offsets defined in registerLIDs are used to
// store the device's F contributions into the F vector.
//
// The Q vector is used for devices that store charge or magnetic
// energy, and since the memristor does none of that, the F vector
// contributions are the whole of the memristor's contribution to the
// full set of DAEs.
//
// @note This method is called by the default implementation of the
// loadDAEVectors master function. Since the MemristorYakopcic class
// reimplements the "Master" "loadDAEVectors" function that loads the
// contributions from all memristor devices in a single loop, THIS
// FUNCTION IS NOT ACTUALLY USED.  The loadDAEFVector method is only
// called when a device class does not re-implement the master class.
// This can be a source of confusion when attempting to modify the MemristorYakopcic
// device, or any other device that reimplements the Master classes.
//
// @see Xyce::Device::MemristorYakopcic::Master::loadDAEVectors
//
// @author Eric Keiter, SNL, Parallel Computational Sciences
// @date   01/24/03
//
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
// Function      : Xyce::Device::MemristorYakopcic::Instance::loadDAEdFdx
// Purpose       : Loads the F-vector contributions for a single
//                 memristor  instance.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Load the DAE the derivative of the F vector with respect to the
// solution vector x, dFdx
//
// Loads the contributions for a single memristor instance to the 
// dFdx matrix (the F contribution to the Jacobian).
//
// This method uses the Jacobian LIDs (row offsets) that were stored by
// registerJacLIDs.
//
// @see Xyce::Device::MemristorYakopcic::Instance::registerJacLIDs
//
// @note This method is called by the default implementation of the
// loadDAEMatrices master function. Since the MemristorYakopcic class
// reimplements the "Master" "loadDAEMatrices" function that loads the
// contributions from all memristor devices in a single loop, THIS
// FUNCTION IS NOT ACTUALLY USED.  The loadDAEdFdx method is only
// called when a device class does not re-implement the master class.
// This can be a source of confusion when attempting to modify the MemristorYakopcic
// device, or any other device that reimplements the Master classes.
//
// @see Xyce::Device::MemristorYakopcic::Master::loadDAEMatrices
//
// @return true on success
//
// @author Eric Keiter, SNL, Parallel Computational Sciences
// @date   03/05/04
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
// Function      : Xyce::Device::MemristorYakopcic::Instance::updateTemperature
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
// The MemristorYakopcic device supports temperature-dependent resistance through its
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
// Function      : Xyce::Device::MemristorYakopcic::Model::processParams
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
// Function      : Xyce::Device::MemristorYakopcic::Model::processInstanceParams
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
// Function      : Xyce::Device::MemristorYakopcic::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Construct a memristor model from a "model block" that was created
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
  Eta_(0.0),
  Vp_(0.0),
  Vn_(0.0),
  Ap_(0.0),
  An_(0.0),
  A1_(0.0),
  A2_(0.0),
  B_(0.0),
  AlphaP_(0.0),
  AlphaN_(0.0),
  XP_(0.0),
  XN_(0.0),
  randomResNoiseOn_(false),
  randomResNoiseSeed_(0),
  randomResNoiseLambda_(0.0),
  randomResNoiseMean_(0.0),
  randomResNoiseSD_(0.0),
  randomResUpdateTime_(0.0),
  randomResEpsilonUpdateTime_(0.0),
  randomResDelta_(0.0),
  randomResDeltaGrad_(0.0)

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

  // if the random noise in the resistance model is onj
  // seed the random number generator
  if( randomResNoiseOn_ )
  {
    randomNumberGen_ = new Xyce::Util::RandomNumbers( randomResNoiseSeed_ );
  }

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorYakopcic::Model::~Model
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
  // if the random noise in the resistance model is onj
  // delete the random number generator
  if( randomResNoiseOn_ )
  {
    delete randomNumberGen_;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_MemristorYakopcicModel::printOutInstances
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
  os << "Number of MemristorYakopcic Instances: " << instanceContainer.size() << std::endl;
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
// Function      : Xyce::Device::MemristorYakopcic::Master::updateState
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Update state for all memristor instances, regardless of model.
//
// @param solVec solution vector
// @param staVec state vector
// @param stoVec store vector
//
// @return true on success
//
// @note Because the memristor device re-implements the base-class
// Master::updateState, the Instance::updatePrimaryState method is never
// called, nor is the Instance::updateIntermediateVars method.  This method
// replaces those, and does the same work but inside a loop over all
// memristor instances.
//
// @see Xyce::Device::MemristorYakopcic::Instance::updatePrimaryState
// @author Eric Keiter, SNL
// @date   11/26/08

bool Master::updateState(double * solVec, double * staVec, double * stoVec)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & ri = *(*it);
    
    // For this model I-V relationship is 
    //
    // i(t) = a1 x(t) sinh( b V(t) )  : iff V(t) >= 0
    //        a2 x(t) sinh( b V(t) )  : iff V(t) < 0 
    //
    // so need to update the x(t) equation and 
    // calculate I(t) based on V(t)
    //
    
    double v_pos    = solVec[ri.li_Pos];
    double v_neg    = solVec[ri.li_Neg];
    double x        = solVec[ri.li_x];
    
    // If RTN noise is enabled, then check if it's time to 
    // update the effective a1 or a2 variable 
    if( ri.model_.randomResNoiseOn_ && 
        (getSolverState().timeStepNumber_ != ri.resNoiseLastUpdateStep) && 
        (fabs(ri.resNoiseLastUpdateTime - getSolverState().currTime_ ) > ri.resNoiseNextUpdateTime) )
    {
      
      ri.resNoiseLastUpdateStep = getSolverState().timeStepNumber_;
      ri.resNoiseLastUpdateTime = getSolverState().currTime_;

      ri.resNoiseNextUpdateTime = -log(ri.model_.randomNumberGen_->uniformRandom() ) * ri.model_.randomResNoiseLambda_*ri.model_.randomResUpdateTime_;
      
      // the above is only valid if the mean is small (i.e. less than 1).  Check if this is better 
      //int aNum = ri.model_.randomNumberGen_->poissonRandom(ri.model_.randomResNoiseLambda_);
      // the actual update interval is the poisson distributed number times the base update time.
      //ri.resNoiseNextUpdateTime = aNum * ri.model_.randomResUpdateTime_;
      

      // update the resistance factor 
      if( ri.resNoiseHiLoState == 0 )
      {
        // in low state so move to high state 
        ri.resNoiseHiLoState = 1;
        ri.resNoiseRfactor = 1 + (0.5 * ri.model_.randomResDelta_ * ri.model_.randomResDeltaGrad_);
      }
      else
      {
        // in high state so move to low state 
        ri.resNoiseHiLoState = 0;
        ri.resNoiseRfactor = 1 - (0.5 * ri.model_.randomResDelta_ * ri.model_.randomResDeltaGrad_);
      }

    }
    
    
    // If RTN noise for X's growth rate is enabled, then check if it's time to 
    // update the effective a1 or a2 variable 
    if( ri.model_.randomXNoiseOn_ && 
        (getSolverState().timeStepNumber_ != ri.xNoiseLastUpdateStep) && 
        (fabs(ri.xNoiseLastUpdateTime - getSolverState().currTime_ ) > ri.xNoiseNextUpdateTime) )
    {
      
      ri.xNoiseLastUpdateStep = getSolverState().timeStepNumber_;
      ri.xNoiseLastUpdateTime = getSolverState().currTime_;

      ri.xNoiseNextUpdateTime = -log(ri.model_.randomNumberGen_->uniformRandom() ) * ri.model_.randomXNoiseLambda_*ri.model_.randomResUpdateTime_;
      
      // the above is only valid if the mean is small (i.e. less than 1).  Check if this is better 
      //int aNum = ri.model_.randomNumberGen_->poissonRandom(ri.model_.randomXNoiseLambda_);
      // the actual update interval is the poisson distributed number times the base update time.
      //ri.XNoiseNextUpdateTime = aNum * ri.model_.randomResUpdateTime_;
      
      // unlike the resistance noise, this will be gaussian noise added to 
      // the growth rate
      ri.xNoiseFactor = 1.0 - 2.0 * (ri.model_.randomNumberGen_->uniformRandom()-0.5) * ri.model_.randomXDelta_;
      /*
      // update the growth rate factor 
      if( ri.xNoiseHiLoState == 0 )
      {
        // in low state so move to high state 
        ri.xNoiseHiLoState = 1;
        ri.xNoiseFactor = 1 + (0.5 * ri.model_.randomXDelta_ * ri.model_.randomXDeltaGrad_);
      }
      else
      {
        // in high state so move to low state 
        ri.xNoiseHiLoState = 0;
        ri.xNoiseFactor = 1 - (0.5 * ri.model_.randomXDelta_ * ri.model_.randomXDeltaGrad_);
      }
      */

    }
    
    
    
    
    if( ri.model_.randomResNoiseOn_ )
    {
      stoVec[ri.li_store_tdt] = ri.resNoiseNextUpdateTime;
    }
    
    {
      // extra scope to limit variable extent.  Otherwise one will have lots of 
      // similar variables 
      //Sacado::Fad::SFad<double,2> varDeltaV( 2, 0, (v_pos-v_neg) );
      Sacado::Fad::SFad<double,3> varV1( 3, 0, v_pos );
      Sacado::Fad::SFad<double,3> varV2( 3, 1, v_neg );
      Sacado::Fad::SFad<double,3> varX( 3, 2, x );
      Sacado::Fad::SFad<double,3> paramA1( ri.model_.A1_ / ri.resNoiseRfactor);
      Sacado::Fad::SFad<double,3> paramA2( ri.model_.A2_ / ri.resNoiseRfactor);
      Sacado::Fad::SFad<double,3> paramB( ri.model_.B_ );
      Sacado::Fad::SFad<double,3> resultFad;
      I_V_Fxn( varV1, varV2, varX, paramA1, paramA2, paramB, resultFad );  
          
      ri.i0 = resultFad.val();
      ri.G  = resultFad.dx(0);
      ri.dIdx = resultFad.dx(2);
    }
    
    {
      // evaluate the state variable equation 
      //Sacado::Fad::SFad<double,2> varDeltaV( 2, 0, (v_pos-v_neg) );
      Sacado::Fad::SFad<double,3> varV1( 3, 0, v_pos );
      Sacado::Fad::SFad<double,3> varV2( 3, 1, v_neg );
      Sacado::Fad::SFad<double,3> varX( 3, 2, x );
      Sacado::Fad::SFad<double,3> paramAp( ri.model_.Ap_ );
      Sacado::Fad::SFad<double,3> paramAn( ri.model_.An_ );
      Sacado::Fad::SFad<double,3> paramVp( ri.model_.Vp_ );
      Sacado::Fad::SFad<double,3> paramVn( ri.model_.Vn_ );
      Sacado::Fad::SFad<double,3> paramXp( ri.model_.XP_ );
      Sacado::Fad::SFad<double,3> paramXn( ri.model_.XN_ );
      Sacado::Fad::SFad<double,3> paramAlphaP( ri.model_.AlphaP_ * ri.xNoiseFactor);
      Sacado::Fad::SFad<double,3> paramAlphaN( ri.model_.AlphaN_ * ri.xNoiseFactor);
      Sacado::Fad::SFad<double,3> paramEta( ri.model_.Eta_ );
      Sacado::Fad::SFad<double,3> resultFad;
      
      X_var_F_equation( varV1, varV2, varX,  
        paramAp, paramAn, paramVp, paramVn, paramXp, paramXn, paramAlphaP, paramAlphaN, paramEta, resultFad );
      
      ri.xVarFContribution = resultFad.val();
      ri.dxFEqdVpos = resultFad.dx(0);
      ri.dxFEqdVneg = resultFad.dx(1);
      ri.dxFEqdx = resultFad.dx(2);
      
    }
    
    
    //double Reff     = ri.model_.ROn_;
    //
    // Calclate Reff and deriative term dReff/dx 
    //
    {
      // extra scope to limit variable extent.  Otherwise one will have lots of 
      // similar variables 
      /*
      Sacado::Fad::SFad<double,1> varX( 1, 0, x );
      Sacado::Fad::SFad<double,1> paramRon( ri.resNoiseRfactor*ri.model_.ROn_ );
      Sacado::Fad::SFad<double,1> paramRoff( ri.resNoiseRfactor*ri.model_.ROff_ );
      Sacado::Fad::SFad<double,1> paramXon( ri.model_.xOn_  * ri.model_.xScaling_  );
      Sacado::Fad::SFad<double,1> paramXoff( ri.model_.xOff_  * ri.model_.xScaling_ );
      Sacado::Fad::SFad<double,1> resultFad;
      if(ri.IVrelation==0)
      {
        // liner current voltage relationship
        //Xyce::dout() << ri.getName() << "About to call ReffLin..." << ri.model_.ROn_ << " , " << ri.model_.ROff_ << ", " << ri.model_.xOn_ << ", " << ri.model_.xOff_ << std::endl;

        ReffLin( varX, paramRon, paramRoff, paramXon, paramXoff, resultFad );
        //fval = Ron + (X - Xon) * (Roff-Ron)/(Xoff - Xon);
        //dReff/dV1 = 0
        //dReff/dV2 = 0;
        //dReff/dx = (Roff-Ron)/(Xoff - Xon);
      }
      else if( ri.IVrelation==1 )
      {
        //Xyce::dout() << ri.getName() << "About to call ReffNonLin..." << std::endl;
        ReffNonLin( varX, paramRon, paramRoff, paramXon, paramXoff, resultFad );
      }
      ri.Reff = resultFad.val();
      ri.dReffdx = resultFad.dx(0);
      ri.dReffdvpos =  (v_pos-v_neg)*(-pow( ri.Reff, 2 ) * ri.dReffdx) / (ri.model_.xScaling_);
      ri.dReffdvneg = -(v_pos-v_neg)*(-pow( ri.Reff, 2 ) * ri.dReffdx) / (ri.model_.xScaling_);
      //ri.dReffdvpos = 0.0;  
      //ri.dReffdvneg = 0.0; 
      */
    }

    //ri.G = 1.0/(ri.Reff);
    //ri.i0 = (v_pos-v_neg)*ri.G;

    // if the random noise in the resistance model is on, then
    // with the estimate for the current check if it's high enough to 
    // trigger a random change in ROn/ROff
    // and this is the start of a new step, then get a new 
    // gaussian distributed variable
    

    {
      // calculate the xVar F equation update and the 
      // derivatives d(F(x))/dVpos, d(F(x))/dVneg and d(F(x))/dx
      /*
      Sacado::Fad::SFad<double,3> varVpos( 3, 0, v_pos);
      Sacado::Fad::SFad<double,3> varVneg( 3, 1, v_neg);
      Sacado::Fad::SFad<double,3> varX( 3, 2, x );
      Sacado::Fad::SFad<double,3> parami( ri.i0 );
      Sacado::Fad::SFad<double,3> paramG( ri.G );
      Sacado::Fad::SFad<double,3> paramiOff( ri.model_.iOff_ );
      Sacado::Fad::SFad<double,3> paramiOn( ri.model_.iOn_ );

      //Sacado::Fad::SFad<double,3> paramkOff( ri.model_.kOff_ );
      //Sacado::Fad::SFad<double,3> paramkOn( ri.model_.kOn_ ); 
      Sacado::Fad::SFad<double,3> paramkOff( ri.model_.kOff_*ri.model_.xScaling_ );
      Sacado::Fad::SFad<double,3> paramkOn( ri.model_.kOn_ *ri.model_.xScaling_ );
      Sacado::Fad::SFad<double,3> paramAlphaOff( ri.model_.alphaOff_ );
      Sacado::Fad::SFad<double,3> paramAlphaOn( ri.model_.alphaOn_ );
      Sacado::Fad::SFad<double,3> resultFad;
  
      xVarFterm(varVpos, varVneg, varX, paramG, paramiOff, paramiOn, paramkOff, paramkOn,
        paramAlphaOff, paramAlphaOn, resultFad);

      // calculate the window function value
      Sacado::Fad::SFad<double,3> windowFunctionFad;

      if( ri.model_.windowType_ ==  0) 
      {
        // no window
        windowFunctionFad = 1.0;
      }
      else if( ri.model_.windowType_ ==  1) 
      {
        Sacado::Fad::SFad<double,3> paramD( ri.model_.D_ );
        Sacado::Fad::SFad<double,3> paramP( ri.model_.p_ );
        JogelkarWindowFunction( varX, paramD, paramP, windowFunctionFad );
      }
      else if( ri.model_.windowType_ ==  2) 
      {
        Sacado::Fad::SFad<double,3> paramD( ri.model_.D_ );
        Sacado::Fad::SFad<double,3> paramP( ri.model_.p_ );
        BiolekWindowFunction( varX, paramD, paramP, parami, windowFunctionFad );
      }
      else if( ri.model_.windowType_ ==  3) 
      {
        Sacado::Fad::SFad<double,3> paramD( ri.model_.D_ );
        Sacado::Fad::SFad<double,3> paramP( ri.model_.p_ );
        Sacado::Fad::SFad<double,3> paramJ( ri.model_.j_ );
        ProdromakisWindowFunction( varX, paramD, paramP, paramJ, windowFunctionFad );
      }
      else if( ri.model_.windowType_ ==  4) 
      {
        Sacado::Fad::SFad<double,3> paramAOn( ri.model_.aOn_ *ri.model_.xScaling_ );
        Sacado::Fad::SFad<double,3> paramAOff( ri.model_.aOff_ *ri.model_.xScaling_ );
        Sacado::Fad::SFad<double,3> paramWc( ri.model_.wc_ *ri.model_.xScaling_ );
        TEAMWindowFunctionF( varX, parami, paramAOff, paramAOn, paramWc, windowFunctionFad );
      }

      // xVarFterm * WindowFunction
      // dertivative terms then become: WindowFunction dxVarFterm/dv1  + xVarFterm dW/dv1
       
      ri.xVarFContribution = resultFad.val()*windowFunctionFad.val();
      ri.dxFEqdVpos = windowFunctionFad.val()*resultFad.dx(0);
      ri.dxFEqdVneg = windowFunctionFad.val()*resultFad.dx(1);
      ri.dxFEqdx = windowFunctionFad.val()*resultFad.dx(2) + resultFad.val()*windowFunctionFad.dx(2) / (ri.model_.xScaling_);
      */
    }
   
    //ri.G = 1.0/(rfactor*ri.Reff);
    //ri.i0 = (v_pos-v_neg)*ri.G;

  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorYakopcic::Master::loadDAEVectors
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Load DAE vectors of all memristor instances, regardless of model
//
// @param solVec solution vector
// @param fVec f vector
// @param qVec q vector
// @param leadF store lead current f vector
// @param leadQ store lead current q vector
//
// @return true on success
//
// @note Because the memristor device re-implements the base-class
// Master::loadDAEVectors, the Instance::loadDAEFVector method is
// never called.  This method replaces those, and does the same work
// but inside a loop over all memristor instances.
//
// @see Xyce::Device::MemristorYakopcic::Instance::loadDAEFVector
//
// @author Eric Keiter, SNL
// @date   11/26/08

bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, double * leadF, double * leadQ, double * junctionV)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & ri = *(*it);
    fVec[ri.li_Pos] += ri.i0;
    fVec[ri.li_Neg] += -ri.i0;
    fVec[ri.li_x]   += ri.xVarFContribution;
    qVec[ri.li_x]   -= solVec[ri.li_x];
    if( getSolverState().dcopFlag )
    {
      qVec[ri.li_x] -= ri.XO_;
    }
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
// Function      : Xyce::Device::MemristorYakopcic::Master::loadDAEMatrices
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Load DAE matrices for all memristor instances, regardless of model
//
// @param dFdx matrix of derivatives of F vector with respect to solution
// @param dQdx matrix of derivatives of Q vector with respect to solution
//
// @return true on success
//
// @note Because the memristor device re-implements the base-class
// Master::loadDAEMatrices, the Instance::loadDAEdFdx method is
// never called.  This method replaces those, and does the same work
// but inside a loop over all memristor instances.
//
// @see Xyce::Device::MemristorYakopcic::Instance::loadDAEdFdx
//
// @author Eric Keiter, SNL
// @date   11/26/08

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
// Function      : Xyce::Device::MemristorYakopcic::Traits::factory
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Create a new instance of the MemristorYakopcic device.
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
// Function      : Xyce::Device::MemristorYakopcic::registerDevice
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
// The device is declared here to be an "memristor" device, which must 
// take a model card of type "memristor".  This device will correspond to model
// level 3 of memristor models.
void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() ||
      ((deviceMap.find("MEMRISTOR")!=deviceMap.end()) && (levelSet.find(3)!=levelSet.end())))
  {
    MemristorTEAM::registerDevice();

    Config<Traits>::addConfiguration()
      .registerDevice("memristor", 3)
      .registerModelType("memristor", 3);
  }
}

//-----------------------------------------------------------------------------
// Function      : memristorYakopcicSensitivity::operator
// Purpose       : produces df/dp and dq/dp, where p=R.  
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
void memristorYakopcicSensitivity::operator()(
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

} // namespace MemristorYakopcic
} // namespace Device
} // namespace Xyce

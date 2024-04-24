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
// Purpose        :
//
// Special Notes  :  "Creator and Creation Date" are actually the dates
//                   this file was created.  When it was created it was
//                   just a placeholder.
//                   Actual implementation of the Xyce voltage controlled
//                   switch occurred on 5/22/2001, and was done by Tom Russo,
//                   SNL, Component Information and Models.
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <N_UTL_Math.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SW.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

#include <N_UTL_Expression.h>

namespace Xyce {
namespace Device {

namespace SW {

void Traits::loadInstanceParameters(ParametricData<SW::Instance> &p)
{
// Set up double precision variables:
    p.addPar ("CONTROL",0.0,&SW::Instance::CONTROL)
     .setExpressionAccess(ParameterType::SOLN_DEP)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("");

    // Set up non-double precision variables:
    p.addPar ("ON",false,&SW::Instance::ON);
    p.addPar ("OFF",false,&SW::Instance::OFF);
}

void Traits::loadModelParameters(ParametricData<SW::Model> &p)
{
    p.addPar ("RON",1.0,&SW::Model::RON)
     .setUnit(U_OHM)
     .setCategory(CAT_NONE)
     .setDescription("On resistance");

    p.addPar ("ROFF",1.0e6,&SW::Model::ROFF)
     .setUnit(U_OHM)
     .setCategory(CAT_NONE)
     .setDescription("Off resistance");

    p.addPar ("VON",1.0,&SW::Model::VON)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("On voltage");
     
    p.addPar ("VHON",1.0,&SW::Model::VHON)
     .setGivenMember(&SW::Model::VHONGiven)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("On voltage with hysteresis");
     
    p.addPar ("VOFF",0.0,&SW::Model::VOFF)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("Off voltage");
     
    p.addPar ("VHOFF",0.0,&SW::Model::VHOFF)
     .setGivenMember(&SW::Model::VHOFFGiven)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("Off voltage with hysteresis");
     
    p.addPar ("ION",0.001,&SW::Model::ION)
     .setUnit(U_AMP)
     .setCategory(CAT_NONE)
     .setDescription("On current");
     
    p.addPar ("IHON",0.001,&SW::Model::IHON)
     .setGivenMember(&SW::Model::IHONGiven)
     .setUnit(U_AMP)
     .setCategory(CAT_NONE)
     .setDescription("On current with hysteresis");

    p.addPar ("IOFF",0.0,&SW::Model::IOFF)
     .setUnit(U_AMP)
     .setCategory(CAT_NONE)
     .setDescription("Off current with hysteresis");
     
    p.addPar ("IHOFF",0.0,&SW::Model::IOFF)
     .setGivenMember(&SW::Model::IHOFFGiven)
     .setUnit(U_AMP)
     .setCategory(CAT_NONE)
     .setDescription("Off current");

    p.addPar ("ON",1.0,&SW::Model::ON)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("On control value");
     
    p.addPar ("ONH",1.0,&SW::Model::ONH)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("On control value with hysteresis");
     
    p.addPar ("OFFH",0.0,&SW::Model::OFFH)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Off control value with hysteresis");
    
    p.addPar ("OFF",0.0,&SW::Model::OFF)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Off control value");
}

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/16/05
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & IB,
  Model & SWiter,
  const FactoryBlock &  factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(SWiter),
    Exp_ptr(0),
    R(0.0),
    ON(false),
    G(0.0),
    SW_STATE(0.0),
    switch_state(0.0),
    li_switch_state(-1),
    li_Pos(-1),
    li_Neg(-1),
    li_branch_data(-1),
    APosEquPosNodeOffset(-1),
    APosEquNegNodeOffset(-1),
    ANegEquPosNodeOffset(-1),
    ANegEquNegNodeOffset(-1),
    fPosEquPosNodePtr(0),
    fPosEquNegNodePtr(0),
    fNegEquPosNodePtr(0),
    fNegEquNegNodePtr(0)
{
  numIntVars = 0;
  numExtVars = 2;
  numStateVars = 1;
  setNumStoreVars(1);                  // Initialize number if store variables in DeviceInstance
  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1;    // this is the space to allocate if lead current or power is needed.


  jacStamp.resize(2);
  jacStamp[0].resize(2);
  jacStamp[0][0]=0;
  jacStamp[0][1]=1;
  jacStamp[1].resize(2);
  jacStamp[1][0]=0;
  jacStamp[1][1]=1;

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (given("OFF"))
  {
    if (given("ON"))
    {
      UserError(*this) << "Cannot specify both 'on' and off' for switch";
    }
    ON = !OFF;
  }
  if (!given("CONTROL"))
  {
    UserError(*this) << "Must specify 'control' for switch";
  }

  std::vector<Depend>::const_iterator d;
  std::vector<Depend>::const_iterator begin = getDependentParams().begin();
  std::vector<Depend>::const_iterator end = getDependentParams().end();

  for  (d = begin ; d != end ; ++d)
  {
    if (d->name == "CONTROL")
    {
      expNumVars = d->numVars;
      expBaseVar = d->lowVarIndex;
      Exp_ptr = d->expr;

      jacStamp[0].resize(2+expNumVars);
      jacStamp[1].resize(2+expNumVars);
      for( int i = 0; i < expNumVars; ++i )
      {
        jacStamp[0][2+i] = 2+i;
        jacStamp[1][2+i] = 2+i;
      }
      expVarDerivs.resize(expNumVars);
      myVarVals.resize(expNumVars);

      dependentParamExcludeMap_[d->name] = 1;
    }
  }

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                                     const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       : Note that the SW does not have any state vars.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/20/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // copy over the global ID lists:
  staLIDVec = staLIDVecRef;
  li_switch_state = staLIDVec[0];
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       : One store var for device current.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/05/2013
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef )
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());

  // copy over the global ID lists.
  stoLIDVec = stoLIDVecRef;

  li_control = stoLIDVec[0];
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::registerBranchDataLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 12/21/2012
//-----------------------------------------------------------------------------

void Instance::registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef)
{
  AssertLIDs(branchLIDVecRef.size() == getNumBranchDataVars());

  if (loadLeadCurrent)
  {
    li_branch_data= branchLIDVecRef[0];
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::getDepSolnVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/06/01
//-----------------------------------------------------------------------------
const std::vector<std::string> & Instance::getDepSolnVars()
{
  return DeviceInstance::getDepSolnVars();
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/2/02
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
// Creation Date : 09/2/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNegNodeOffset = jacLIDVec[0][1];
  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];

  APosEquControlNodeOffset.resize( expNumVars );
  ANegEquControlNodeOffset.resize( expNumVars );
  for( int i = 0; i < expNumVars; ++i )
  {
    APosEquControlNodeOffset[i] = jacLIDVec[0][2+i];
    ANegEquControlNodeOffset[i] = jacLIDVec[1][2+i];
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/30/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  fPosEquPosNodePtr = &(dFdx[li_Pos][APosEquPosNodeOffset]);
  fPosEquNegNodePtr = &(dFdx[li_Pos][APosEquNegNodeOffset]);
  fNegEquPosNodePtr = &(dFdx[li_Neg][ANegEquPosNodeOffset]);
  fNegEquNegNodePtr = &(dFdx[li_Neg][ANegEquNegNodeOffset]);

  fPosEquControlNodePtr.resize( expNumVars );
  fNegEquControlNodePtr.resize( expNumVars );
  for( int i = 0; i < expNumVars; ++i )
  {
    fPosEquControlNodePtr[i] = &(dFdx[li_Pos][APosEquControlNodeOffset[i]]);
    fNegEquControlNodePtr[i] = &(dFdx[li_Neg][ANegEquControlNodeOffset[i]]);
  }
#endif
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
  double * staVec = extData.nextStaVectorRawPtr;
  bool bsuccess = updateIntermediateVars ();

  //  obtain the current value of the switch state
  switch_state = SW_STATE;

  staVec[li_switch_state] = switch_state;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/29/01
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  double current_state, last_state, current_stateHysOn, current_stateHysOff;
  double control;
  int i;

  double * solVec = extData.nextSolVectorRawPtr;
  double * stoVec = extData.nextStoVectorRawPtr;
  double * stoVecLast = extData.lastStoVectorRawPtr;

  // Evaluate Expression with corrected time derivative values

  Exp_ptr->evaluate( expVal, expVarDerivs );
  control = expVal;

  // This is not really correct, an interim hack.  This is supposed
  // to be where we deal with the specification of ON or OFF from the
  // netlist.
  if (getSolverState().initJctFlag_)
  {
    if (ON)
      current_state = 1;
    else
      current_state = 0;
  }
  else
  {
    // current state with NO hysteresis
    current_state = (control-model_.OFF)*model_.dInv;
    // current state with hysteresis on the ON state 
    current_stateHysOn = (control-model_.OFF)*model_.dInvHOn;
    // current state with hysteresis on the OFF state
    current_stateHysOff = (control-model_.OFFH)*model_.dInvHOff;
  }

  stoVec[li_control] = current_state;
  last_state = stoVecLast[li_control];
  v_pos = solVec[li_Pos];
  v_neg = solVec[li_Neg];

  if (current_state >= 1.0)
  {
    R = model_.RON;
    G = 1.0/R;
    for (i=0 ; i<expNumVars ; ++i)
      expVarDerivs[i] = 0;
  }
  else if ( current_state <= 0.0)
  {
    R = model_.ROFF;
    G = 1.0/R;
    for (i=0 ; i<expNumVars ; ++i)
      expVarDerivs[i] = 0;
  }
  else
  {
    current_state = 2*current_state - 1;
    G = exp(-model_.Lm - 0.75*model_.Lr*current_state +
		0.25*model_.Lr*current_state*current_state*current_state);
    R = 1.0/G;
    for (i=0 ; i<expNumVars ; ++i)
    {
      expVarDerivs[i] = G * (1.5 * (current_state*current_state-1) * model_.Lr *
                                model_.dInv * expVarDerivs[i]);
    }
  }

  return true;
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
  if (loadLeadCurrent)
    addBranchDataNode( symbol_table, li_branch_data, getName(), "BRANCH_D");
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one switch instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 1/10/01
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 resistor  instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/13/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * solVec = extData.nextSolVectorRawPtr;
  double * fVec = extData.daeFVectorRawPtr;

  // load RHS vector element for the positive circuit node KCL equ.
  double coef = (v_pos-v_neg)*G;

  fVec[li_Pos] += coef;
  fVec[li_Neg] += -coef;
  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    double * junctionV = extData.nextJunctionVCompRawPtr;
    leadF[li_branch_data] = coef;
    junctionV[li_branch_data] = solVec[li_Pos] - solVec[li_Neg];
  }

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
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/13/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Pos][APosEquPosNodeOffset] += G;
  dFdx[li_Pos][APosEquNegNodeOffset] -= G;
  dFdx[li_Neg][ANegEquPosNodeOffset] -= G;
  dFdx[li_Neg][ANegEquNegNodeOffset] += G;

  if( expNumVars )
  {
    for( int i = 0; i < expNumVars; ++i )
    {
      dFdx[li_Pos][APosEquControlNodeOffset[i]] += (v_pos-v_neg) * expVarDerivs[i];
      dFdx[li_Neg][ANegEquControlNodeOffset[i]] -= (v_pos-v_neg) * expVarDerivs[i];
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
  dtype(1),
  RON(0.0),
  ROFF(0.0),
  ON(0.0),
  OFF(0.0)
{
  if (getType() != "")
  {
    if (getType() == "SWITCH" ) {
      dtype = 1;
    }
    else if (getType() == "ISWITCH") {
      dtype = 2;
    }
    else if (getType() == "VSWITCH") {
      dtype = 3;
    }
    else
    {
      UserError(*this) << "Unrecognized model type " << getType();
    }
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

  if (dtype == 2)
  {
    if (!given("ON"))
      ON = ION;
    if (!given("OFF"))
      OFF = IOFF;
      
    if( IHONGiven )
      ONH = IHON;
    else
    {
      ONH = ON;
      IHON = ION;
    }
    if( IHOFFGiven )
      OFFH = IHOFF;
    else
    {
      OFFH = OFF;
      IHOFF = IOFF;
    }
  }
  else if (dtype == 3)
  {
    if (!given("ON"))
      ON = VON;
    if (!given("OFF"))
      OFF = VOFF;
    
    if( VHONGiven )
      ONH = VHON;
    else
    {
      ONH = ON;
      VHON = VON;
    }
    if( VHOFFGiven )
      OFFH = VHOFF;
    else
    {
      OFFH = OFF;
      VHOFF = VOFF;
    }
  }

  processParams();
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/16/05
//-----------------------------------------------------------------------------
bool Model::processParams ()
{
  double del, delOnH, delOffH;
  Lm = log (sqrt(RON*ROFF));
  Lr = log (RON/ROFF);

  del = ON-OFF;
  delOnH = ONH-OFF;
  delOffH = ON-OFFH;

  if (del < 0 && del > -1e-12)
    del = -1e-12;
  if (del >= 0 && del < 1e-12)
    del = 1e-12;
  dInv = 1.0/del;
  
  if (delOnH < 0 && delOnH > -1e-12)
    delOnH = -1e-12;
  if (delOnH >= 0 && delOnH < 1e-12)
    delOnH = 1e-12;
  dInvHOn = 1.0/delOnH;
  
  if (dInvHOff < 0 && dInvHOff > -1e-12)
    dInvHOff = -1e-12;
  if (dInvHOff >= 0 && dInvHOff < 1e-12)
    dInvHOff = 1e-12;
  dInvHOff = 1.0/delOffH;
  
  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirely, PSSI
// Creation Date : 03/23/06
//----------------------------------------------------------------------------
bool Model::processInstanceParams()
{

  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    (*iter)->processParams();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Model::~Model ()
{
  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
    delete (*iter);
}

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
    os << "    R = " << (*iter)->R;
    os << "  G = " << (*iter)->G;
    os << "  State = " << (*iter)->SW_STATE;
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

//-----------------------------------------------------------------------------
// SW Master functions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  bool bsuccess = true;

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & si = *(*it);

    bool btmp = si.updateIntermediateVars ();
    bsuccess = bsuccess && btmp;

    //  obtain the current value of the switch state
    si.switch_state = si.SW_STATE;
    staVec[si.li_switch_state] = si.switch_state;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateSecondaryState ( double * staDerivVec, double * stoVec )
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & si = *(*it);

    double current_state, last_state, current_stateHysOn, current_stateHysOff;;
    double control;
    double * solVec = si.extData.nextSolVectorRawPtr;
    double * stoVec = si.extData.nextStoVectorRawPtr;
    double * stoVecLast = si.extData.lastStoVectorRawPtr;

    // Evaluate Expression with corrected time derivative values

    si.Exp_ptr->evaluate( si.expVal, si.expVarDerivs );
    control = si.expVal;

    // This is not really correct, an interim hack.  This is supposed
    // to be where we deal with the specification of ON or OFF from the
    // netlist.
    if (getSolverState().initJctFlag_)
    {
      if (si.ON)
        current_state = 1;
      else
        current_state = 0;
    }
    else
    {
      // current state with NO hysteresis
      current_state = (control-si.getModel().OFF)*si.getModel().dInv;
      // current state with hysteresis on the ON state 
      current_stateHysOn = (control-si.getModel().OFF)*si.getModel().dInvHOn;
      // current state with hysteresis on the OFF state
      current_stateHysOff = (control-si.getModel().OFFH)*si.getModel().dInvHOff;
    }
    stoVec[si.li_control] = current_state;
    last_state = stoVecLast[si.li_control];
        
    si.v_pos = solVec[si.li_Pos];
    si.v_neg = solVec[si.li_Neg];

    if ((current_state >= 1.0) || ((last_state >= 1.0) && (current_stateHysOn >= 1.0)))
    {
      si.R = si.getModel().RON;
      si.G = 1.0/si.R;
      for (int i=0 ; i<si.expNumVars ; ++i)
        si.expVarDerivs[i] = 0;
    }
    else if (( current_state <= 0.0) || ((last_state <= 0.0) && (current_stateHysOff <= 0.0)))
    {
      si.R = si.getModel().ROFF;
      si.G = 1.0/si.R;
      for (int i=0 ; i<si.expNumVars ; ++i)
        si.expVarDerivs[i] = 0;
    }
    else
    {
      current_state = 2*current_state - 1;
      si.G = exp(-si.getModel().Lm - 0.75*si.getModel().Lr*current_state +
      0.25*si.getModel().Lr*current_state*current_state*current_state);
      si.R = 1.0/si.G;
      for (int i=0 ; i<si.expNumVars ; ++i)
      {
        si.expVarDerivs[i] = si.G * (1.5 * (current_state*current_state-1) * si.getModel().Lr *
                                  si.getModel().dInv * si.expVarDerivs[i]);
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, double * leadF, double * leadQ, double * junctionV)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
     Instance & si = *(*it);
    // F-vector:
    double coef = (si.v_pos-si.v_neg)*si.G;
    fVec[si.li_Pos] += coef;
    fVec[si.li_Neg] += -coef;
    if( si.loadLeadCurrent )
    {
      leadF[si.li_branch_data] = coef;
      junctionV[si.li_branch_data] = solVec[si.li_Pos] - solVec[si.li_Neg];
    }
    // Q-vector:
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
     Instance & si = *(*it);
#ifndef Xyce_NONPOINTER_MATRIX_LOAD

    *si.fPosEquPosNodePtr += si.G;
    *si.fPosEquNegNodePtr -= si.G;
    *si.fNegEquPosNodePtr -= si.G;
    *si.fNegEquNegNodePtr += si.G;

    if( si.expNumVars )
    {
      for( int j = 0; j < si.expNumVars; ++j )
      {
        *si.fPosEquControlNodePtr[j] += (si.v_pos-si.v_neg) * si.expVarDerivs[j];
        *si.fNegEquControlNodePtr[j] -= (si.v_pos-si.v_neg) * si.expVarDerivs[j];
      }
    }
#else

    dFdx[si.li_Pos][si.APosEquPosNodeOffset] += si.G;
    dFdx[si.li_Pos][si.APosEquNegNodeOffset] -= si.G;
    dFdx[si.li_Neg][si.ANegEquPosNodeOffset] -= si.G;
    dFdx[si.li_Neg][si.ANegEquNegNodeOffset] += si.G;

    if( si.expNumVars )
    {
      for( int i = 0; i < si.expNumVars; ++i )
      {
        dFdx[si.li_Pos][si.APosEquControlNodeOffset[i]] += (si.v_pos-si.v_neg) * si.expVarDerivs[i];
        dFdx[si.li_Neg][si.ANegEquControlNodeOffset[i]] -= (si.v_pos-si.v_neg) * si.expVarDerivs[i];
      }
    }
#endif
  }
  return true;
}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() || (deviceMap.find("S")!=deviceMap.end())
                        || (deviceMap.find("W")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("s", 1)
      .registerModelType("switch", 1)
      .registerModelType("iswitch", 1)
      .registerModelType("vswitch", 1)
      ;
  }
}

} // namespace SW
} // namespace Device
} // namespace Xyce

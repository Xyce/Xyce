//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Special Notes  :
//
// Creator        : Christina Warrender, SNL, Cognitive Modeling
//
// Creation Date  : 10/12/11
//
//
//
//
//----------------------------------------------------------------------------

#include <Xyce_config.h>
//#define Xyce_FullSynapseJac 1

// ---------- Standard Includes ----------
// used to get time in seconds to seed random number generator.
#include<time.h>

// ----------   Xyce Includes   ----------
//
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Synapse4.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_DEV_Synapse.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {

namespace Synapse4 {

void Traits::loadInstanceParameters(ParametricData<Synapse4::Instance> &p)
{
  p.addPar ("GMAX",0.01,&Synapse4::Instance::gMax)
   .setGivenMember(&Synapse4::Instance::gMaxGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("Maximal Synaptic Conductance");
}

void Traits::loadModelParameters(ParametricData<Synapse4::Model> &p)
{
  p.addPar ("VTHRESH",0.01,&Synapse4::Model::vThresh)
    .setUnit(U_VOLT)
    .setDescription("Presynaptic voltage spike threhsold");

  p.addPar ("DELAY",0.001,&Synapse4::Model::delay)
    .setUnit(U_SECOND)
    .setDescription("Time delay between presynaptic signal and postsynaptic response");

  p.addPar ("GMAX",0.01,&Synapse4::Model::gMax)
    .setUnit(U_OHMM1)
    .setDescription("Maximal Synaptic Conductance");

  p.addPar ("EREV",0.0,&Synapse4::Model::eRev)
    .setUnit(U_VOLT)
    .setDescription("Postsynaptic Reversal Potential");

  p.addPar ("TAU1",0.0001,&Synapse4::Model::tau1)
    .setUnit(U_SECM1)
    .setDescription("Rise time constant");

  p.addPar ("TAU2",0.01,&Synapse4::Model::tau2)
    .setUnit(U_SECM1)
    .setDescription("Decay time constant");
}

std::vector< std::vector<int> > Instance::jacStamp;

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  // initialization
  respondTime = std::numeric_limits<double>::max( );
  ready = true;
  //active = false;

  // Set any non-constant parameter defaults:

  if (!gMaxGiven )
  {
    gMax = model_.gMax;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & IB,
  Model & Riter,
  const FactoryBlock &  factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Riter),
    li_Prev(-1),
    li_Post(-1),
    li_A0_store(-1),
    li_B0_store(-1),
    li_t0_store(-1),
    li_branch_data(-1),
#ifdef Xyce_FullSynapseJac
    APostEquPostNodeOffset(-1),
    f_PostEquPostNodePtr(0),
#endif
    ipost(0),
    didVpost(0),
    randInitialized(false)
{
  numIntVars   = 0;   // A and B   2
  numExtVars   = 2;   // presynaptic V and postsynaptic V
  setNumStoreVars(3);   // A0, B0, t0
  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1;    // this is the space to allocate if lead current or power is needed.

  if( jacStamp.empty() )
  {
    jacStamp.resize(2);
    jacStamp[0].resize(0);    // presynaptic V not changed
#ifdef Xyce_FullSynapseJac
    jacStamp[1].resize(1);    // postsynaptic V depends on itself
#else
    jacStamp[1].resize(0);    // postsynaptic V depends on itself
#endif
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

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
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef )
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  SynapseInstance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  li_Prev = extLIDVec[0];
  li_Post = extLIDVec[1];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << "  li_Prev = " << li_Prev << std::endl;
    Xyce::dout() << "  li_Post = " << li_Post << std::endl;
  }

  /*
    li_AVar = intLIDVec[0];
    li_BVar = intLIDVec[1];
  */

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << section_divider << std::endl;
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
  addStoreNode(symbol_table, li_A0_store, getName().getEncodedName() + "_A0");
  addStoreNode(symbol_table, li_B0_store, getName().getEncodedName() + "_B0");
  addStoreNode(symbol_table, li_t0_store, getName().getEncodedName() + "_T0");

  if (loadLeadCurrent)
    addBranchDataNode(symbol_table, li_branch_data, getName(), "BRANCH_D");
}


//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef )
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());

  // copy over the global ID lists.
  stoLIDVec = stoLIDVecRef;

  li_A0_store = stoLIDVec[0];
  li_B0_store = stoLIDVec[1];
  li_t0_store = stoLIDVec[2];

}


//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::registerBranchDataLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 12/18/2012
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
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
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
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
#ifdef Xyce_FullSynapseJac
  APostEquPostNodeOffset = jacLIDVec[1][0];
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
#ifdef Xyce_FullSynapseJac
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  f_PostEquPostNodePtr = &(dFdx[li_Post][APostEquPostNodeOffset]);
#endif
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one synapse instance
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars()
{
  bool bsuccess = true;

  double * lastSolVecPtr = extData.lastSolVectorRawPtr;
  double * stoVector = extData.nextStoVectorRawPtr;

  // This check need to
  // Check for presynaptic spike start, set time to respond
  double vPre  = lastSolVecPtr[li_Prev];
  double vPost = lastSolVecPtr[li_Post];
  double time = getSolverState().currTime_;
  double vThresh = model_.vThresh;
  double delay = model_.delay;
  if (ready)
  {
    if (vPre > vThresh)
    {
      ready=false;
      respondTime = time + delay;
    }
  }
  else  // already had spike start, looking for end
  {
    if (vPre < vThresh)
    {
      ready=true;
    }
  }

  // only need to update state variables and synaptic current when synapse active
  // disabling testing for this because determining the point at which the synaptic current
  // becomes negligible is somewhat problematic, and there don't appear to be significant savings
  // from avoiding loads when synapse is inactive
  //if (active)
  {
    // handle decay of A and B
    double tau1 = model_.tau1;
    double tau2 = model_.tau2;

    double t0 = stoVector[li_t0_store];
    double Anow = stoVector[li_A0_store] * exp( - ( time - t0 ) / tau1 );
    double Bnow =  stoVector[li_B0_store] * exp( - ( time - t0 ) / tau2 );

    // set up variables for load methods
    // current equation is the same whether we're responding to spike or not,
    // assuming current A and B values set appropriately above
    // ipost = (B-A)*(V-Erev)

    double eRev = model_.eRev;
    double cond = Bnow-Anow;
    ipost = cond*(vPost-eRev);
    didVpost = cond;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << std::endl << section_divider << std::endl;
      Xyce::dout() << "  SynapseInstance::updateIntermediateVars" << std::endl;
      Xyce::dout() << "Anow:  " << Anow << std::endl;
      Xyce::dout() << "Bnow:  " << Bnow << std::endl;
      Xyce::dout() << "vPost:  " << vPost << std::endl;
      Xyce::dout() << "eRev:  " << eRev << std::endl;
      Xyce::dout() << "ipost:  " << ipost << std::endl;
      Xyce::dout() << "didVpost:  " << didVpost << std::endl;
    }

  }	// end if synapse active

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  bool bsuccess = updateIntermediateVars();
  return  bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState()
{
  return  true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 Synapse4
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
//                 The "Q" vector contains charges and fluxes, mostly.
//                 However, it is ordered like the solution vector, and as
//                 it is part of the KCL formulation, the terms in Q will
//                 actually be *sums* of charges, rather than single
//                 distinct charges.
//
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  return true;
}
//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 Synapse4  instance.
//
// Special Notes : This is an algebraic constaint, and as such the Synapse4
//                 does make a contribution to it.
//
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  Linear::Vector *  fVecPtr = extData.daeFVectorPtr;
  (*fVecPtr)[li_Prev] += 0.0;
  (*fVecPtr)[li_Post] += ipost;
  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    leadF[li_branch_data] = ipost;
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the dQdx-matrix contributions for a single
//                 Synapse4 instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
//                 The "Q" vector contains charges and fluxes, mostly.
//                 However, it is ordered like the solution vector, and as
//                 it is part of the KCL formulation, the terms in Q will
//                 actually be *sums* of charges, rather than single
//                 distinct charges.
//
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 Synapse4  instance.
//
// Special Notes : This is an algebraic constaint, and as such the Synapse4
//                 does make a contribution to it.
//
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
#ifdef Xyce_FullSynapseJac
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  dFdx[li_Post][APostEquPostNodeOffset] += didVpost;
#endif
  return true;
}

//-----------------------------------------------------------------------------
// Function      : InstanceInstance::outputPlotFiles
// Purpose       : If requested by the user output all the variables
//                 associated with the population
// Special Notes : We're actually using this method not for output, but because
//                 it's called after the system has converged.  In this case,
//                 that let us mark the end of handling a presynaptic event.
// Scope         : public
// Creator       : Christina Warrender, SNL, Cognitive Modeling
// Creation Date : 10/25/2011
//-----------------------------------------------------------------------------
bool Instance::outputPlotFiles(bool force_final_output)
{
  bool bsuccess = true;

  // cew 11/3/11:  changing this back to just storing A0 and B0 for updateIntermediateVars
  // to use in calculating A, B, and ipost
  // But when incrementing A0 (B0), current value of A (B) must be used.

  double time = getSolverState().currTime_;
  //Linear::Vector * staVector = extData.nextStaVectorPtr;
  double * stoVector = extData.nextStoVectorRawPtr;
  if (time >= respondTime)
  {
    //Xyce::dout() << "Instance::outputPlotFiles() adjusting A0, B0 and t0 = " << getSolverState().currTime_ << std::endl;
    // succesfully processed a step, so adjust the next
    // respondTime to distant future
    respondTime = std::numeric_limits<double>::max( );
    //active = true;

    // now we also update the A0, B0 and t0 in the state vector
    double factor = model_.factor;
    double deltaAB = factor*gMax;
    double tau1 = model_.tau1;
    double tau2 = model_.tau2;
    double t0 = stoVector[li_t0_store];
    double Anow = stoVector[li_A0_store] * exp( - ( time - t0 ) / tau1 );
    double Bnow =  stoVector[li_B0_store] * exp( - ( time - t0 ) / tau2 );
    stoVector[li_A0_store] = Anow + deltaAB;
    stoVector[li_B0_store] = Bnow + deltaAB;
    stoVector[li_t0_store] = getSolverState().currTime_;
    //Xyce::dout() << "A:  " << (*staVector)[li_A0_state] << " B: " << (*staVector)[li_B0_state] << std::endl;
  }

/* disabling this because determining the point at which the synaptic current becomes negligible
   is somewhat problematic, and there don't appear to be significant savings from avoiding loads
   when synapse is inactive
   if (active)
   {
   // determine when we can stop updating synaptic current, when it's negligible
   // can't use magnitude of current, because it starts off small initially
   // instead, determine whether time since last spike is greater than reasonable multiple of
   // larger time constant
   if ( (time-(*staVector)[li_t0_state])>3.0*model_.maxtau )
   {
   active = false;
   ipost = 0.0;
   didVpost = 0.0;
   }
   }
*/

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp_tmp)
{
  bool bsuccess = true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Model::processParams ()
{
  // initialize variables needed to calculate synaptic dynamics
  if (tau1/tau2 > .9999) {
    tau1 = .9999*tau2;
  }
  tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1);
  factor = -exp(-tp/tau1) + exp(-tp/tau2);
  factor = 1/factor;
  maxtau = (tau1>tau2)?tau1:tau2;

  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirely, PSSI
// Creation Date : 10/12/11
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
// Function      : Model::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    vThresh(0.0),
    gMax(0.0),
    delay(0.0),
    eRev(0.0),
    tau1(0.0),
    tau2(0.0)
{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams();
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
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
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i,isize;
  isize = instanceContainer.size();
  os << std::endl;
  os << "Number of Synapse4 Instances: " << isize << std::endl;
  os << "    name     model name  Parameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
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


//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL, Cognitive Modeling
// Creation Date : 07/18/12
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  if( getSolverState().dcopFlag )
  {
    // no synaptic activity during dcop
    return true;
  }

  bool bsuccess = true;

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    bool tmpBool = (*it)->updateIntermediateVars();	// skipping call to updatePrimaryState
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL, Cognitive Modeling
// Creation Date : 7/18/12
//-----------------------------------------------------------------------------
bool Master::updateSecondaryState (double * staDerivVec, double *stoVec)
{
  // not used for this device
  return true;
}


//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL, Cognitive Modeling
// Creation Date : 07/16/12
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors(double * solVec, double * fVec, double * qVec, double * bVec, double * leadF, double * leadQ, double * junctionV)
{
  if( getSolverState().dcopFlag )
  {
    // no synaptic activity during dcop
    return true;
  }

  bool bsuccess = true;

  // We don't need Q loads for this synapse device at all
  // Only need to do F loads if there's a nonzero synaptic current
  // disabling testing for the latter because determining the point at which the synaptic current
  // becomes negligible is somewhat problematic, and there don't appear to be significant savings
  // from avoiding loads when synapse is inactive

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    bool tmpBool = (*it)->loadDAEFVector();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL, Cognitive Modeling
// Creation Date : 07/16/12
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  if( getSolverState().dcopFlag )
  {
    // no synaptic activity during dcop
    return true;
  }

  bool bsuccess = true;

  // We don't need Q loads for this synapse device at all
  // Only need to do F loads if there's a nonzero synaptic current
  // disabling testing for the latter because determining the point at which the synaptic current
  // becomes negligible is somewhat problematic, and there don't appear to be significant savings
  // from avoiding loads when synapse is inactive

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    bool tmpBool = (*it)->loadDAEdFdx();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() ||
      ((deviceMap.find("SYNAPSE")!=deviceMap.end()) && (levelSet.find(4)!=levelSet.end())))
  { 
    Synapse::registerDevice();

    Config<Traits>::addConfiguration()
      .registerDevice("synapse", 4)
      .registerModelType("synapse", 4);
  }
}

} // namespace Synapse4
} // namespace Device
} // namespace Xyce

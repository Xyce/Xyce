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

//-------------------------------------------------------------------------
//
// Purpose        : Simple Clopath synapse
//
// Special Notes  :
//
// Creator        : Rich Schiek, SNL, Electrical Systems Modeling
//
// Creation Date : 01/25/2011
//
//
//
//
//-------------------------------------------------------------------------


#include <Xyce_config.h>
//#define Xyce_FullSynapseJac 1

// ----------   Xyce Includes   ----------
//
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Synapse3.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_DEV_Synapse.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {

namespace Synapse3 {

std::mt19937 * Instance::randomNumberGeneratorPtr_=0;

void Traits::loadInstanceParameters(ParametricData<Synapse3::Instance> &p)
{
  p.addPar ("GMAX",0.01,&Synapse3::Instance::gMax)
   .setGivenMember(&Synapse3::Instance::gMaxGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("Maximal Synaptic Conductance");

  p.addPar ("P",1.0,&Synapse3::Instance::transmissionProbability)
   .setGivenMember(&Synapse3::Instance::transmissionProbabilityValueGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Transmission Probability");

  p.addPar ("WINIT",0.01,&Synapse3::Instance::wInitialValue)
   .setGivenMember(&Synapse3::Instance::wInitialValueGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Synaptic weight,initial value");
}

void Traits::loadModelParameters(ParametricData<Synapse3::Model> &p)
{
  p.addPar ("VTHRESH",0.01,&Synapse3::Model::vThresh)
    .setUnit(U_VOLT)
    .setDescription("Presynaptic voltage spike threhsold");

  p.addPar ("DELAY",0.001,&Synapse3::Model::delay)
    .setUnit(U_SECOND)
    .setDescription("Time delay between presynaptic signal and postsynaptic response");

  p.addPar ("GMAX",0.01,&Synapse3::Model::gMax)
    .setUnit(U_OHMM1)
    .setDescription("Maximal Synaptic Conductance");

  p.addPar ("EREV",0.0,&Synapse3::Model::eRev)
    .setUnit(U_VOLT)
    .setDescription("Postsynaptic Reversal Potential");

  p.addPar ("TAU1",0.0001,&Synapse3::Model::tau1)
    .setUnit(U_SECM1)
    .setDescription("Rise time constant");

  p.addPar ("TAU2",0.01,&Synapse3::Model::tau2)
    .setUnit(U_SECM1)
    .setDescription("Decay time constant");

  p.addPar ("S",0.01,&Synapse3::Model::sParam)
    .setUnit(U_VOLT)
    .setDescription("Voltage threshold for a spike event");

  p.addPar ("R",0.01,&Synapse3::Model::rParam)
    .setUnit(U_VOLT)
    .setDescription("Resting voltage for resting event");

  p.addPar ("WMIN",0.01,&Synapse3::Model::wMin)
    .setDescription("Synaptic weight,minimum value");

  p.addPar ("WMAX",0.01,&Synapse3::Model::wMax)
    .setDescription("Synaptic weight,maximum value");

  p.addPar ("WINIT",0.01,&Synapse3::Model::wInitialValue)
    .setDescription("Synaptic weight,initial value");

  p.addPar ("L1TAU",0.01,&Synapse3::Model::vL1tau1)
    .setUnit(U_SECOND)
    .setDescription("Rate for Longterm potentiation factor (LPF) based on post-synaptic voltage (rate 1)");

  p.addPar ("L2TAU",0.01,&Synapse3::Model::vL2tau2)
    .setUnit(U_SECOND)
    .setDescription("Rate for Longterm potentiation factor (LPF) based on post-synaptic voltage (rate 2)");

  p.addPar ("L3TAU",0.01,&Synapse3::Model::vL3tau3)
    .setUnit(U_SECOND)
    .setDescription("Rate for Longterm potentiation factor (LPF) based on pre-synaptic voltage (rate 3)");

  p.addPar ("ALTD",0.01,&Synapse3::Model::aLTD)
    .setDescription("Long term depression coefficient");

  p.addPar ("ALTP",0.01,&Synapse3::Model::aLTP)
    .setDescription("Long term potentiation coefficient");

  p.addPar ("P",1.0,&Synapse3::Model::transmissionProbability)
    .setDescription("Transmission Probability");
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


  // Set any non-constant parameter defaults:

  if (!gMaxGiven )
  {
    gMax = model_.gMax;
  }

  if( !wInitialValueGiven )
  {
    // no value for initial weighting given on instance line.
    // set default from model value
    wInitialValue = model_.wInitialValue;
  }

  if( !transmissionProbabilityValueGiven )
  {
    transmissionProbability = model_.transmissionProbability;
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
    transmissionProbability(1.0),
    li_Prev(-1),
    li_Post(-1),
    li_A0_store(-1),
    li_B0_store(-1),
    li_t0_store(-1),
    li_weight_store(-1),
    li_VL1_store(-1),
    li_VL2_store(-1),
    li_VL3_store(-1),
    li_branch_data(-1),
#ifdef Xyce_FullSynapseJac
    APostEquPostNodeOffset(-1),
    f_PostEquPostNodePtr(0),
#endif
    ipost(0),
    didVpost(0),
    transmissionFactor(1),
    randInitialized(false)
{
  numIntVars   = 0;   // A and B   2
  numExtVars   = 2;   // presynaptic V and postsynaptic V
  setNumStoreVars(7);   // A0, B0, t0, weight, vl1, vl2, vl3
  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1;    // this is the space to allocate if lead current or power is needed.

  if( jacStamp.empty() )
  {
    jacStamp.resize(2);
    jacStamp[0].resize(0);    // presynaptic V not changed
#ifdef Xyce_FullSynapseJac
    jacStamp[1].resize(1);    // postsynaptic V depends on itself
#else
    jacStamp[1].resize(0);
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
  addStoreNode(symbol_table, li_weight_store, getName().getEncodedName() + "_W");
  addStoreNode(symbol_table, li_VL1_store, getName().getEncodedName() + "_VL1");
  addStoreNode(symbol_table, li_VL2_store, getName().getEncodedName() + "_VL2");
  addStoreNode(symbol_table, li_VL3_store, getName().getEncodedName() + "_VL3");

  if (loadLeadCurrent)
    addBranchDataNode(symbol_table, li_branch_data, getName(), "BRANCH_D");
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       : Note that the Synapse3 does not have any state vars.
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
  li_weight_store = stoLIDVec[3];
  li_VL1_store = stoLIDVec[4];
  li_VL2_store = stoLIDVec[5];
  li_VL3_store = stoLIDVec[6];

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
// Purpose       : update intermediate variables for one diode instance
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars()
{
  bool bsuccess = true;

  double * lastSolVecPtr = extData.lastSolVectorRawPtr;
  double * stoVec = extData.nextStoVectorRawPtr;
  double * currStoVec = extData.currStoVectorRawPtr;

  double vPre  = lastSolVecPtr[li_Prev];
  double vPost = lastSolVecPtr[li_Post];

  // initialized random number generator if needed
  if( !randInitialized )
  {
    if (randomNumberGeneratorPtr_ == 0)
    {
      if( getDeviceOptions().randomSeed != 0 )
      {
        randomNumberGeneratorPtr_=new std::mt19937(getDeviceOptions().randomSeed);
      }
      else
      {
        std::random_device rd;
        randomNumberGeneratorPtr_=new std::mt19937(rd());
      }
    }
    randInitialized=true;
  }

  // to do:  Need to adjust this conditional to do this same initialization
  // if the dcop is being skipped and this is the first iteration.
  if( getSolverState().dcopFlag )
  {
    // no firing during DCOP, so postsynaptic current is 0 and unchanging
    ipost = 0.0;
    didVpost = 0.0;
    synapticWeight = wInitialValue;
    synapticWeightUpdate = 0;
    stoVec[li_weight_store] = synapticWeight;
    stoVec[li_A0_store] = 0.0;
    stoVec[li_B0_store] = 0.0;
    stoVec[li_VL1_store] = vPost;
    stoVec[li_VL2_store] = vPost;
    stoVec[li_VL3_store] = 0.0;
  }
  else
  {
    // Check for presynaptic spikes, set time to respond
    double time = getSolverState().currTime_;
    double vThresh = model_.vThresh;
    double delay = model_.delay;
    if (ready)
    {
      if (vPre > vThresh)
      {
        ready=false;
        respondTime = time + delay;
        // check transmissionProbability to see if current is transmitted
        transmissionFactor=1;
        if( transmissionProbability < 1.0 )
        {
          std::uniform_real_distribution<double> uniformRandom(0.0,1.0);
          double arand = uniformRandom(*randomNumberGeneratorPtr_);
          if( arand > transmissionProbability )
          {
            transmissionFactor=0;
          }
        }
      }
    }
    else  // already had spike start, looking for end
    {
      // Need to improve the logic here so that no only is the pre snyapse going down (vPre < vThresh),
      // but synapse must be done with its work. (tau1 + tau2) after delay?
      if( (vPre < vThresh) && ((time - currStoVec[li_t0_store]) > ( model_.tau1 + model_.tau2 )))
      {
        ready=true;
      }
    }

    // handle decay of A and B
    double tau1 = model_.tau1;
    double tau2 = model_.tau2;

    double t0 = stoVec[li_t0_store];
    double Anow = stoVec[li_A0_store] * exp( - ( time - t0 ) / tau1 );
    double Bnow =  stoVec[li_B0_store] * exp( - ( time - t0 ) / tau2 );
    double vl1Now = stoVec[li_VL1_store];
    double vl2Now = stoVec[li_VL2_store];
    double vl3Now = stoVec[li_VL3_store];
    synapticWeight = stoVec[li_weight_store];

    // calculate update to synaptic weight
    synapticWeightUpdate = 0.0;
    // only if the synapse is transmitting should be update the synaptic weight.
    // transmissionFactor is only 0 or 1
    if( transmissionFactor == 1 )
    {
      if( (synapticWeight > model_.wMin) && (synapticWeight < model_.wMax) )
      {
        // weight in in range so update may be non zero.  Check LTD and LTP terms.
        if( (vPre > model_.sParam) && (vl1Now > model_.rParam) )
        {
          synapticWeightUpdate += -(model_.aLTD) * (vl1Now - model_.rParam);
        }
        if( (vPost > model_.sParam) && (vl2Now > model_.rParam) )
        {
          synapticWeightUpdate += (model_.aLTP) * vl3Now * (vPost - model_.sParam) * (vl2Now - model_.rParam);
        }
      }

      double wMin = model_.wMin;
      double wMax = model_.wMax;

      // enforce min/max bounds on weight
      if (synapticWeight+synapticWeightUpdate > wMax)
      {
        synapticWeightUpdate = wMax-synapticWeight;
      }
      else if (synapticWeight+synapticWeightUpdate < wMin)
      {
        synapticWeightUpdate = -(synapticWeight-wMin);
      }
    }


    vl1Update = getSolverState().currTimeStep_ * (vPost - vl1Now) / model_.vL1tau1;
    vl2Update = getSolverState().currTimeStep_ * (vPost - vl2Now) / model_.vL2tau2;
    if( vPre > model_.sParam )
    {
      vl3Update = getSolverState().currTimeStep_ * (1.0 - vl3Now) / model_.vL3tau3;
    }
    else
    {
      vl3Update = getSolverState().currTimeStep_ * (-vl3Now) / model_.vL3tau3;
    }

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << std::endl << section_divider << std::endl;
      Xyce::dout() << "  SynapseInstance::updateIntermediateVars" << std::endl;
      Xyce::dout() << "Anow:  " << Anow << std::endl;
      Xyce::dout() << "Bnow:  " << Bnow << std::endl;
    }

    // set up variables for load methods

    // current equation is the same whether we're responding to spike or not,
    // assuming current A and B values set appropriately above
    // ipost = (B-A)*(V-Erev)
    double eRev = model_.eRev;
    //ipost = (synapticWeight + synapticWeightUpdate) * (Bnow-Anow)*(vPost-eRev);
    ipost = transmissionFactor * synapticWeight * (Bnow-Anow)*(vPost-eRev);

    //didVpost = (synapticWeight + synapticWeightUpdate) * (Bnow - Anow);
    // update synaptic weight after time steps is successful.
    didVpost = synapticWeight * (Bnow - Anow);

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << std::endl << section_divider << std::endl;
      Xyce::dout() << "vPost:  " << vPost << std::endl;
      Xyce::dout() << "eRev:  " << eRev << std::endl;
      Xyce::dout() << "weight: " << synapticWeight << std::endl;
      Xyce::dout() << "weight update: " << synapticWeightUpdate << std::endl;
      Xyce::dout() << "ipost:  " << ipost << std::endl;
      Xyce::dout() << "didVpost:  " << didVpost << std::endl;
    }

  } // end else (!getSolverState().dcopFlag)

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
//                 Synapse3  instance.
//
// Special Notes : This is an algebraic constaint, and as such the Synapse3
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
//                 Synapse3 instance.
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
//                 Synapse3  instance.
//
// Special Notes : This is an algebraic constaint, and as such the Synapse3
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
  double * stoVec = extData.nextStoVectorRawPtr;

  if (time >= respondTime)
  {
    //Xyce::dout() << "Instance::outputPlotFiles() adjusting A0, B0 and t0 = " << getSolverState().currTime_ << ", " << respondTime << std::endl;
    // succesfully processed a step, so adjust the next
    // respondTime to distant future
    respondTime = std::numeric_limits<double>::max( );

    // now we also update the A0, B0 and t0 in the state vector
    double factor = model_.factor;
    double deltaAB = factor*gMax;
    double tau1 = model_.tau1;
    double tau2 = model_.tau2;

    double t0 = stoVec[li_t0_store];
    double Anow = stoVec[li_A0_store] * exp( - ( time - t0 ) / tau1 );
    double Bnow = stoVec[li_B0_store] * exp( - ( time - t0 ) / tau2 );
    stoVec[li_A0_store] = Anow + deltaAB;
    stoVec[li_B0_store] = Bnow + deltaAB;
    stoVec[li_t0_store] = getSolverState().currTime_;
  } // end if time >= respondTime

  stoVec[li_weight_store] += synapticWeightUpdate;
  stoVec[li_VL1_store] += vl1Update;
  stoVec[li_VL2_store] += vl2Update;
  stoVec[li_VL3_store] += vl3Update;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  SynapseInstance::outputPlotFiles" << std::endl;
    Xyce::dout() << "time:  " << getSolverState().currTime_ << std::endl;
    Xyce::dout() << "weight in store:  " << stoVec[li_weight_store] << std::endl;
  }


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
  os << "Number of Synapse3 Instances: " << isize << std::endl;
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


Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new DeviceMaster<Traits>(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() ||
      ((deviceMap.find("SYNAPSE")!=deviceMap.end()) && (levelSet.find(3)!=levelSet.end())))
  { 
    Synapse::registerDevice();

    Config<Traits>::addConfiguration()
      .registerDevice("synapse", 3)
      .registerModelType("synapse", 3);
  }
}

} // namespace Synapse3
} // namespace Device
} // namespace Xyce

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
// Purpose        :  Provide an ideal delay element for ANN work
//
// Special Notes  :  Voltage across first two nodes should be the voltage
//                   at control nodes delayed by given time delay
//                   This behaves as a voltage-controlled voltage source
//                   with unity gain, but whose output is simply delayed.
//                   Borrows the history technique from the
//                   lossless transmission line.
//
// Creator        : Tom Russo
// Creation Date  : 17 Nov 2020
//
//
//
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceState.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Delay.h>
#include <N_DEV_Bsrc.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Functors.h>
#include <N_UTL_MachDepParams.h>

namespace Xyce {
namespace Device {

namespace Delay {

void Traits::loadInstanceParameters(ParametricData<Delay::Instance> &p)
{
  p.addPar("TD", 0.0, &Delay::Instance::TD_)
    .setUnit(U_SECOND)
    .setDescription("Time delay");
  p.addPar("BPENABLED", true, &Delay::Instance::canSetBreakPoints_)
    .setDescription("Can this device set discontinuity breakpoints?");
  p.addPar("EXTRAPOLATION", true, &Delay::Instance::useExtrapolation_)
    .setDescription("Can this device use extrapolation on history?");
  p.addPar("LINEARINTERP", false, &Delay::Instance::useOnlyLinearInterpolation_)
    .setDescription("Should this device use only linear interpolation on history?");
}

void Traits::loadModelParameters(ParametricData<Delay::Model> &p)
{
}

std::vector< std::vector<int> > Instance::jacStamp;

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
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
    f_NegEquBraVarPtr(0),
    timeOld_(-1.0),
    newBreakPoint_(false),
    canSetBreakPoints_(true),
    useExtrapolation_(true),
    useOnlyLinearInterpolation_(false)
{
  numIntVars   = 1;
  numExtVars   = 4;
  numStateVars = 0;
  setNumBranchDataVars(0);          // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1; // this is the space to allocate if lead current or power is needed.

  setNumStoreVars(0);

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

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams();
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       : Check for errors in instance parameters
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  bool bsuccess = true;
  if (!given("TD"))
  {
    UserError(*this) << " Required time delay parameter TD not specified";
    bsuccess=false;
  }
  else
  {
    if (TD_ <= 0)
    {
      UserError(*this) << "Zero or negative time delay.";
      bsuccess = false;
    }
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
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
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------
void Instance::registerLIDs ( const std::vector<int> & intLIDVecRef,
                                        const std::vector<int> & extLIDVecRef )
{
  std::string msg;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  DelayInstance::registerLIDs" << std::endl;
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
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs( const std::vector<int> & stoLIDVecRef )
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerBranchDataLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
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
// Creation Date : 17 Nov 2020
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
// Creation Date : 17 Nov 2020
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
// Creation Date : 17 Nov 2020
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
// Creation Date : 17 Nov 2020
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
// Creation Date : 17 Nov 2020
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
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  return updateIntermediateVars ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  double * solVec = extData.nextSolVectorRawPtr;    // the current guess
  double Vpos, Vneg;
  Vpos = solVec[li_ContPos];
  Vneg = solVec[li_ContNeg];

  if ((getSolverState().dcopFlag))
  {
    v_drop_ = Vpos - Vneg;
  }
  else
  {
    double currentTime = getSolverState().currTime_;
    if (getSolverState().newtonIter == 0 && (currentTime != timeOld_))
    {
      timeOld_ = currentTime;
      if (getSolverState().initTranFlag_)
      {
        v_drop_ = Vpos - Vneg;
        history_.clear();
        history_.push_back(History(-2*TD_,v_drop_));
        history_.push_back(History(-TD_,v_drop_));
        history_.push_back(History(0,v_drop_));
      }
      else
      {
        double delayedTime = currentTime - TD_;
        lastInterpolationConverged_ =
          interpVoltageFromHistory_(delayedTime,v_drop_,
                                   currentTime,(Vpos-Vneg));
      }
    }
    else
    {
      // we're the second or later iteration of the second time step or
      // later.  If the previous iteration used unconverged data to get
      // historical value, we must recompute
      if (!lastInterpolationConverged_)
      {
        double delayedTime = currentTime - TD_;
        lastInterpolationConverged_ =
          interpVoltageFromHistory_(delayedTime,v_drop_,
                                   currentTime,(Vpos-Vneg));
      }
    }
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Instance:isConverged ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/22/08
//-----------------------------------------------------------------------------

inline bool Instance::isConverged()
{

  bool converged = true;

  if ((!getSolverState().dcopFlag) && !(getSolverState().initTranFlag_ &&  getSolverState().newtonIter == 0 ))
  {
    double currentTime = getSolverState().currTime_;
    double currentDeltaV;
    double d1, d2; // derivatives in history
    Linear::Vector *theSolVectorPtr = extData.nextSolVectorPtr;
    std::vector<History>::iterator last = history_.end();
    currentDeltaV = (*theSolVectorPtr)[li_ContPos]- (*theSolVectorPtr)[li_ContNeg];

    double t3=currentTime;
    double v3=currentDeltaV;
    last--;
    double t2=last->t_;
    double v2=last->v_;
    last--;
    double t1=last->t_;
    double v1=last->v_;
    d1 = (v3-v2)/(t3-t2);
    d2 = (v2-v1)/(t2-t1);

    if ((fabs(d1-d2) >= .99*std::max(fabs(d1),fabs(d2))+1))
    {

      if ( (currentTime - (t2 + TD_) ) > getSolverState().bpTol_ )
        converged = false;
    }
  }

  return converged;

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
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------
void Instance::acceptStep()
{
  if (!getSolverState().dcopFlag)
  {
    double currentTime = getSolverState().currTime_;
    double currentDeltaV;
    double d1, d2; // derivatives in history
    Linear::Vector *theSolVectorPtr = extData.nextSolVectorPtr;
    std::vector<History>::iterator last = history_.end();
    last--;
    currentDeltaV = (*theSolVectorPtr)[li_ContPos]- (*theSolVectorPtr)[li_ContNeg];
    history_.push_back(History(currentTime,currentDeltaV));

    // now check for discontinuity of derivatives:
    last = history_.end();
    last--;
    double t3=last->t_;
    double v3=last->v_;
    last--;
    double t2=last->t_;
    double v2=last->v_;
    last--;
    double t1=last->t_;
    double v1=last->v_;
    d1 = (v3-v2)/(t3-t2);
    d2 = (v2-v1)/(t2-t1);
    newBreakPoint_=false;
    if ((fabs(d1-d2) >= .99*std::max(fabs(d1),fabs(d2))+1))
    {
      // derivative changed dramatically, call it a discontinuity at t2
      // set a breakpoint if we have those enabled
      newBreakPointTime_=t2+TD_;

      if ( fabs(currentTime - newBreakPointTime_) > getSolverState().bpTol_ )
        newBreakPoint_ = true;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInstanceBreakPoints
// Purpose       : This function adds break points to a vector of breakpoints.
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------
bool Instance::getInstanceBreakPoints ( std::vector<Util::BreakPoint> & breakPointTimes )
{

  bool bsuccess = true;
  double currentTime = getSolverState().currTime_;
  int timeStep = getSolverState().timeStepNumber_;

  if (newBreakPoint_ && canSetBreakPoints_)
  {
    breakPointTimes.push_back(newBreakPointTime_);
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getMaxTimeStepSize
// Purpose       : Force a maximum time step
// Special Notes : If the last accepted step flagged a discontinuity
//                 we set the max time step to TD_, but otherwise return
//                 what is effectively no maximum.
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------
double Instance::getMaxTimeStepSize ()
{
  return 1e99;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 Delay instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
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
    v_drop=v_drop_;
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
//                 Delay instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------
bool Instance::loadDAEBVector ()
{
  if (!(getSolverState().dcopFlag))
  {
    double * bVec = extData.daeBVectorRawPtr;

    bVec[li_Bra] += v_drop_;
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
// Creation Date : 17 Nov 2020
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
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  varTypeVec.resize(1);
  varTypeVec[0] = 'I';
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInternalState
// Purpose       : Generates an DeviceState object and populates
//                 it with the contents of the history vector for use by
//                 restarts
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------

DeviceState * Instance::getInternalState()
{
  int hsize,i,j;
  // allocate object to return
  DeviceState * myState = new DeviceState;


  myState->ID=getName().getEncodedName();
  // We'll pack our history data into the single vector of doubles
  myState->data.resize(history_.size()*3);
  hsize=history_.size();
  for (i=0;i<hsize;++i)
  {
    j=i*2;
    myState->data[j]=history_[i].t_;
    myState->data[j+1]=history_[i].v_;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << Xyce::section_divider
         << std::endl;
    Xyce::dout() << " In Instance::getInternalState " << std::endl;
    Xyce::dout() << "   name=" << getName() << std::endl;
    Xyce::dout() << "   history size = " << hsize << std::endl;
    Xyce::dout() << "   history  data: " << std::endl;
    for (i = 0 ; i < hsize ; ++i)
    {
      Xyce::dout() << "   (" << history_[i].t_ << ", " << history_[i].v_  << ")"<< std::endl;
    }

    Xyce::dout() << "   DeviceState ID = " << myState->ID << std::endl;
    Xyce::dout() << "   DeviceState data size " << myState->data.size() << std::endl;
    Xyce::dout() << "   Device State data: " << std::endl;
    for (i = 0 ; i < myState->data.size() ; ++i)
    {
      Xyce::dout() << "    " << myState->data[i] << std::endl;
    }
    Xyce::dout() << Xyce::section_divider
         << std::endl;
  }

  return myState;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setInternalState
// Purpose       : Reload history data from restart
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------
bool Instance::setInternalState(const DeviceState &state)
{
  int dsize=state.data.size();
  int hsize,i,j;
  if (getName().getEncodedName() != state.ID)
  {
    DevelFatal(*this).in("Delay::Instance::setInternalState") << "ID(" << state.ID << ") from restart does not match my name (" << getName() << ")";
    return false;
  }

  if (dsize%2 != 0)
  {
    UserError(*this) << "Data size from restart (" << dsize << ") not a multiple of 2";
    return false;
  }

  hsize=dsize/2;
  history_.clear();
  history_.resize(hsize);
  for ( i=0; i<hsize; ++i)
  {
    j=i*2;
    history_[i].t_=state.data[j];
    history_[i].v_=state.data[j+1];
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << Xyce::section_divider
         << std::endl;
    Xyce::dout() << " In Instance::setInternalState " << std::endl;
    Xyce::dout() << "   name=" << getName() << std::endl;
    Xyce::dout() << "   history size = " << hsize << std::endl;
    Xyce::dout() << "   history  data: " << std::endl;
    for (i = 0 ; i < hsize ; ++i)
    {
      Xyce::dout() << "   (" << history_[i].t_ << ", " << history_[i].v_ << ")"<< std::endl;
    }

    Xyce::dout() << "   DeviceState ID = " << state.ID << std::endl;
    Xyce::dout() << "   DeviceState data size " << state.data.size() << std::endl;
    Xyce::dout() << "   Device State data: " << std::endl;
    for (i = 0 ; i < state.data.size() ; ++i)
    {
      Xyce::dout() << "    " << state.data[i] << std::endl;
    }
    Xyce::dout() << Xyce::section_divider
         << std::endl;
  }
  return true;
}
// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : "Model block" constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
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
// Creation Date : 17 Nov 2020
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
// Creation Date : 17 Nov 2020
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
// Creation Date : 17 Nov 2020
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
// Function      : Instance::InterpVoltageFromHistory_
// Purpose       : Use 3-point lagrange interpolation to determine
//                 the control voltage at a time in the past
// Special Notes : By default, if the history does not contain a point
//                 with time *after* the desired time, then we do extrapolation
//                 If the user has disabled extrapolation, then in this case
//                 we use the unconverged current value of time as the third
//                 interpolation point
//                 There is also an instance parameter that forces
//                 us to use only linear interpolation.
// Scope         : private
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------
bool Instance::interpVoltageFromHistory_(double t, double & voltage,
                                const double currentT, const double currentV)
{

  std::vector<History>::iterator first = history_.begin();
  std::vector<History>::iterator it1;
  std::vector<History>::iterator last = history_.end();
  double t1,t2,t3;
  double v1,v2,v3;
  bool usedOnlyConverged=true;
  if (history_.size() <= 0)
  {
    DevelFatal(*this).in("Instance::interpVoltageFromHistory_")
      << " interpVoltageFromHistory_ called but history list is"
      << " empty.";
  }

  if (t - first->t_ < -Util::MachineDependentParams::MachinePrecision())
  {
    UserError(*this) << "Cannot interpolate to a time (" << t << ") prior to oldest("
                     << first->t_ << ") in history";
    return false;
  }

  --last;

  // If we are within roundoff of the endpoints of the history, just use
  // the endpoints, otherwise interpolate to get it.
  if ( fabs(t-first->t_)<Util::MachineDependentParams::MachinePrecision())
  {
    voltage = first->v_;
  }
  else if ( fabs(t-last->t_)<Util::MachineDependentParams::MachinePrecision())
  {
    voltage = last->v_;
  }
  else
  {
        // If there are no elements of history with time later than t, there
    // is no point using lower_bound to search for one.

    it1=last;
    if (it1->t_ < t)
    {
      if (useExtrapolation_)
      {
        // just use this last point, which will cause us to extrapolate
        t3=it1->t_;
        v3=it1->v_;
        --it1;
      }
      else
      {
        // we don't extroplate, but instead use unconverged current value
        // to interpolate
        t3=currentT;
        v3=currentV;
        usedOnlyConverged=false;
      }
    }
    else
    {
      LessThan<History,double> lessFunct;
      it1 = lower_bound(history_.begin(),history_.end(),t,lessFunct);
      t3=it1->t_;
      v3=it1->v_;
      --it1;
    }
    t2=it1->t_;
    v2=it1->v_;
    --it1;
    t1=it1->t_;
    v1=it1->v_;
    //  now we have three values of each function to be interpolated, and three
    // times.  t3 is after the desired time, t1 and t2 are before (t2 might be
    // equal to the desired time)

    // Compute derivatives for discontinuity detection
    double d1 = (v3-v2)/(t3-t2);
    double d2 = (v2-v1)/(t2-t1);

    // If discontinuous, or user has asked for purely linear interpolation
    if (useOnlyLinearInterpolation_ || fabs(d1-d2) >= .99*std::max(fabs(d1),fabs(d2))+1)
    {
      // linear
      if (fabs(v3-v2)<Util::MachineDependentParams::MachinePrecision())
      {
        // this is a really pathological case where the history
        // after a breakpoint is totally flat.  Either extrapolation or
        // interpolation should just be the average of the two
        voltage = (v3+v2)/2.0;
      }
      else
      {
        voltage = v2+d1*(t-t2);
      }
    }
    else
    {
      // If we're not doing linear interpolation because of discontinuities,
      // and the user didn't specify linear-only,
      // then we're doing quadratic lagrange interpolation.

      // Set up the differences for lagrange interpolation:
      double dt12,dt13,dt23;
      double dt1,dt2,dt3;
      double f1,f2,f3;
      dt12 = t1-t2;
      dt13 = t1-t3;
      dt23 = t2-t3;
      dt1 = t-t1;
      dt2 = t-t2;
      dt3 = t-t3;
      // now we set up the lagrange interpolating functions
      // e.g. f1 = (t-t2)*(t-t3)/((t1-t2)*(t1-t3))
      // so that fi is 1 at ti and 0 at the other times.
      f1 = dt2*dt3;
      f2 = dt1*dt3;
      f3 = dt1*dt2;
      if (dt12 != 0)
      {
        f1 /= dt12;
        f2 /= -dt12;
      }
      else
      {
        f1 = f2 = 0.0;
      }
      if (dt13 != 0)
      {
        f1 /= dt13;
        f3 /= -dt13;
      }
      else
      {
        f1 = f2 = 0.0;
      }
      if (dt23 != 0)
      {
        f2 /= dt23;
        f3 /= -dt23;
      }
      else
      {
        f2 = f3 = 0.0;
      }
      // that's it, we have the interpolation functions evaluated at the time t,
      // and the values of v1 and v2 at the points, perform  the interpolation
      voltage = f1*v1+f2*v2+f3*v3;
    }
  }
  return usedOnlyConverged;
}

// History member (trivial) functions

//-----------------------------------------------------------------------------
// Function      : History::History
// Purpose       : default constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------
History::History()
  : t_(0),v_(0)
{
}

//-----------------------------------------------------------------------------
// Function      : History::History
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------
History::~History()
{
}
//-----------------------------------------------------------------------------
// Function      : History::History
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------
History::History(const History &right)
  : t_(right.t_),v_(right.v_)
{
}

//-----------------------------------------------------------------------------
// Function      : History::History
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------
History::History(double t, double v)
  : t_(t),v_(v)
{
}


Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new DeviceMaster<Traits>(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() || (deviceMap.find("DELAY")!=deviceMap.end()))
  {

    Config<Traits>::addConfiguration()
      .registerDevice("delay", 1)
      .registerModelType("delay", 1);
  }
}

} // namespace Delay
} // namespace Device
} // namespace Xyce

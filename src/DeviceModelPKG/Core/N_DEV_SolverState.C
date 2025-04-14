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

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 08/25/03
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_ANP_AnalysisManager.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_SolverState.h>
#include <N_NLS_NonLinInfo.h>
#include <N_TIA_WorkingIntegrationMethod.h>
#include <N_TIA_StepErrorControl.h>
#include <N_UTL_Expression.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Diagnostic.h>

#include <expressionGroup.h>
#include <N_DEV_ExpressionGroupWrapper.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : SolverState::SolverState
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/25/03
//-----------------------------------------------------------------------------
SolverState::SolverState ()
  : isPDESystem_(false),
    pdt_(0.0),
    currentOrder_(0),
    usedOrder_(0),
    //currTimeStep_(0.0),
    //lastTimeStep_(0.0),
    currTimeStep_(1.0e-10),
    lastTimeStep_(1.0e-10),
    oldeTimeStep_(1.0e-10),
    currTime_(0.0),
    finalTime_(0.0),
    startingTimeStep_(0.0),
    bpTol_(0.0),
    acceptedTime_(0.0),
    mpdeOnFlag_(false),
    currFastTime_(0.0),
    blockAnalysisFlag_(false),
    spAnalysisFlag_(false),
    earlyNoiseFlag_(false),
    doubleDCOPEnabled(false),
    doubleDCOPStep   (0),
    timeStepNumber_(0),
    ltraDevices_(false),
    ltraTimeIndex_(0),
    ltraTimeHistorySize_(0),
    ltraDoCompact_(false),
    ltraTimePoints_(),
    newtonIter       (0),
    continuationStepNumber (0),
    firstContinuationParam (true),
    firstSolveComplete (false),
    initTranFlag_(false),
    beginIntegrationFlag_(true),
    dcopFlag         (true),
    inputOPFlag      (false),
    transientFlag    (true),
    dcsweepFlag      (false),
    tranopFlag       (true),
    acopFlag         (false),
    noiseFlag         (false),
    locaEnabledFlag  (false),
    externalInitJctFlag_(false),
    externalStateFlag_(false),
    initJctFlag_(false),
    initFixFlag      (false),
    sweepSourceResetFlag(false),
    debugTimeFlag    (false),
    twoLevelNewtonCouplingMode (Nonlinear::FULL_PROBLEM),
    pdeAlpha_(1.0),
    PDEcontinuationFlag_(false),
    chargeHomotopy_(false),
    chargeAlpha_(1.0),
    artParameterFlag_(false),
    gainScale_(1.0),
    nltermScale_(1.0),
    sizeParameterFlag_(false),
    sizeScale_(1.0),
    currFreq_(0.0),
    groupWrapperPtr_(0)
{
  groupWrapperPtr_ = new Xyce::Device::expressionGroupWrapper();
}

SolverState::~SolverState ()
{
  delete groupWrapperPtr_;
}

//-----------------------------------------------------------------------------
// Function      : SolverState::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/05/2005
//-----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream & os, const SolverState & ss)
{

  os << section_divider << std::endl;
  os << "  Device Package Solver State:" << std::endl;

  os << "  pdt = " << ss.pdt_ << std::endl;
  os << "  currTimeStep = " << ss.currTimeStep_ << std::endl;
  os << "  lastTimeStep = " << ss.lastTimeStep_ << std::endl;
  os << "  oldeTimeStep = " << ss.oldeTimeStep_ << std::endl;
  os << "  currTime = " << ss.currTime_ << std::endl;
  os << "  finalTime = " << ss.finalTime_ << std::endl;
  os << "  startingTimeStep = " << ss.startingTimeStep_ << std::endl;
  os << "  bpTol = " << ss.bpTol_ << std::endl;

  os << "  acceptedTime = " << ss.acceptedTime_ << std::endl;
  os << "  currentOrder = " << ss.currentOrder_ << std::endl;
  os << "  usedOrder =  " << ss.usedOrder_ << std::endl;

  os << "  mpdeOnFlag = ";
  if (ss.mpdeOnFlag_)
  {
    os << "yes" << std::endl;
    os << "  currFastTime = " << ss.currFastTime_ << std::endl;
    os << "  blockAnalysisFlag  = " << ss.blockAnalysisFlag_ << std::endl;
  }
  else
  {
    os << "no" << std::endl;
  }

  os << "  timeStepNumber = " << ss.timeStepNumber_ << std::endl;
  os << "  ltraDevices = " << ss.ltraDevices_ << std::endl;
  os << "  ltraTimeIndex = " << ss.ltraTimeIndex_ << std::endl;
  os << "  ltraTimeStepHistorySize = " << ss.ltraTimeHistorySize_ << std::endl;
  os << "  ltraDoCompact = " << ss.ltraDoCompact_ << std::endl;
  os << "  newtonIter = " << ss.newtonIter << std::endl;
  os << "  continuationStepNumber = " << ss.continuationStepNumber << std::endl;
  os << "  firstContinuationParam = ";
  if (ss.firstContinuationParam) os << "yes" << std::endl;
  else                           os << "no" << std::endl;

  os << "  firstSolveComplete = ";
  if (ss.firstSolveComplete) os << "yes" << std::endl;
  else                       os << "no" << std::endl;

  os << "  currFreq = " << ss.currFreq_ << std::endl;

  os << "  initTranFlag = " << (ss.initTranFlag_ ? "yes" : "no") << std::endl;
  os << "  beginIntegrationFlag = " << (ss.beginIntegrationFlag_ ? "yes" : "no") << std::endl;

  os << "  dcopFlag = ";
  if (ss.dcopFlag) os << "yes" << std::endl;
  else             os << "no" << std::endl;

  os << "  inputOPFlag = ";
  if (ss.inputOPFlag) os << "yes" << std::endl;
  else             os << "no" << std::endl;

  os << "  transientFlag = ";
  if (ss.transientFlag) os << "yes" << std::endl;
  else                  os << "no" << std::endl;

  os << "  dcsweepFlag = ";
  if (ss.dcsweepFlag) os << "yes" << std::endl;
  else                os << "no" << std::endl;

  os << "  tranopFlag = ";
  if (ss.tranopFlag) os << "yes" << std::endl;
  else               os << "no" << std::endl;

  os << "  acopFlag = ";
  if (ss.acopFlag) os << "yes" << std::endl;
  else             os << "no" << std::endl;

  os << "  noiseFlag = ";
  if (ss.noiseFlag) os << "yes" << std::endl;
  else             os << "no" << std::endl;

  os << "  isPDESystem_ = ";
  if (ss.isPDESystem_) os << "yes" << std::endl;
  else                  os << "no" << std::endl;

  os << "  locaEnabledFlag = ";
  if (ss.locaEnabledFlag) os << "yes" << std::endl;
  else                    os << "no" << std::endl;

  os << "  initJctFlag = " << (ss.initJctFlag_ ? "yes" : "no") << std::endl;
  os << "  initFixFlag = ";
  if (ss.initFixFlag) os << "yes" << std::endl;
  else                os << "no" << std::endl;

  os << "  sweepSourceResetFlag = ";
  if (ss.sweepSourceResetFlag) os << "yes" << std::endl;
  else                         os << "no" << std::endl;

  os << "  debugTimeFlag = ";
  if (ss.debugTimeFlag) os << "yes" << std::endl;
  else                  os << "no" << std::endl;

  os << section_divider << std::endl;
  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : setupSolverInfo
//
// Purpose       :
//
// Special Notes : This function gets called a lot, and this can be kind of
//                 confusing.  Probably, this function is being used to handle
//                 too many different types of data.
//
//                 For example, it gets called at the beginning
//                 of the "updateSources" function.  Why?  Because the sources
//                 need to know which step we are at, and/or what the current
//                 time is, to do their update properly.
//
//                 However, at the moment, updateSources also provides
//                 information that setupSolverInfo needs.  Only after the
//                 sources have been updated do we know if a sweep source has
//                 been reset.  And, the sweepSourceResetFlag  is used by
//                 setupSolverInfo, to set up the initJctFlag boolean.
//                 So it will need to be called at least one more
//                 time, at the beginning of the RHS load.  (which it is).
//
//                 Anyway, any of the functions that are called from the
//                 outside, such as:  updateSources, loadRHSVector,
//                 loadJacobianMatrix, etc.... have no way of knowing when,
//                 w.r.t. the solvers they are being called.  The only way
//                 to do this properly is to have each of them request
//                 the current solver state, before they go to do their
//                 work.  Hence, this function is called a lot.
//
//                 Unfortunately, this has led to a somewhat sloppy and
//                 confusing interface between the solvers and the
//                 device package.  I wanted to avoid having a lot of
//                 function arguments being passed around for each of
//                 these functions, in part because the calling code
//                 (NLS) doesn't know everything. - NLS knows about the newton
//                 step, but it doesn't know the time step, for example.
//
//                 At some point I hope to refactor this.
//
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/30/01
//-----------------------------------------------------------------------------
bool setupSolverInfo(
  SolverState &                         solver_state,
  const Analysis::AnalysisManager &     analysis_manager,
  bool                                  all_devices_converged,
  const DeviceOptions &                 device_options,
  const Nonlinear::NonLinInfo &         nonlinear_info)
{
  bool bsuccess = true;

  // Time integrator information
  solver_state.pdt_                 = analysis_manager.getWorkingIntegrationMethod().partialTimeDeriv();        // system_state.pdt;
  solver_state.currentOrder_        = analysis_manager.getWorkingIntegrationMethod().getOrder();                // system_state.currentOrder;
  solver_state.usedOrder_           = analysis_manager.getWorkingIntegrationMethod().getUsedOrder();            // system_state.usedOrder;
  solver_state.integrationMethod_   = analysis_manager.getWorkingIntegrationMethod().getMethod(); 

  // Time step error control information
  solver_state.currTimeStep_        = analysis_manager.getStepErrorControl().currentTimeStep;                   // system_state.nextTimeStep;
  solver_state.lastTimeStep_        = analysis_manager.getStepErrorControl().lastTimeStep;                      // system_state.currTimeStep;
  solver_state.oldeTimeStep_        = analysis_manager.getStepErrorControl().oldeTimeStep;                      // 
  solver_state.currTime_            = analysis_manager.getStepErrorControl().nextTime;                          // system_state.nextTime;
  solver_state.finalTime_           = analysis_manager.getStepErrorControl().finalTime;                         // system_state.finalTime;
  solver_state.startingTimeStep_    = analysis_manager.getStepErrorControl().startingTimeStep;                  // system_state.startingTimeStep;
  solver_state.bpTol_               = analysis_manager.getStepErrorControl().getBreakPointLess().tolerance_;    // system_state.bpTol;

  solver_state.currFreq_        = analysis_manager.getCurrentFreq();

  if (solver_state.mpdeOnFlag_)
  {
    solver_state.initTranFlag_ =  true;
    solver_state.beginIntegrationFlag_ =  true;
  }
  else
  {
//    solver_state.dcopFlag             = analysis_manager.getAnalysisObject().getDCOPFlag();             // system_state.dcopFlag;
    solver_state.initTranFlag_         = analysis_manager.getInitTranFlag();                             // system_state.initTranFlag;
    solver_state.beginIntegrationFlag_ = analysis_manager.getBeginningIntegrationFlag();                 // system_state.beginIntegrationFlag;
  }

  solver_state.dcopFlag             = analysis_manager.getAnalysisObject().getDCOPFlag();             // system_state.dcopFlag;
  solver_state.inputOPFlag          = analysis_manager.getAnalysisObject().getInputOPFlag();
  solver_state.acopFlag             = analysis_manager.getACOPFlag();
  solver_state.noiseFlag             = analysis_manager.getNoiseFlag(); 
  solver_state.tranopFlag           = analysis_manager.getTranOPFlag();
  solver_state.transientFlag        = analysis_manager.getTransientFlag();
  solver_state.dcsweepFlag          = analysis_manager.getDCSweepFlag();
  solver_state.sweepSourceResetFlag = analysis_manager.getSweepSourceResetFlag();

  solver_state.timeStepNumber_      = analysis_manager.getStepNumber();

  solver_state.doubleDCOPStep       = analysis_manager.getDoubleDCOPStep();
  solver_state.doubleDCOPEnabled    = analysis_manager.getDoubleDCOPEnabled();

  // Nonlinear solver info:
  solver_state.newtonIter           = nonlinear_info.newtonIter;
  solver_state.twoLevelNewtonCouplingMode         = nonlinear_info.twoLevelNewtonCouplingMode;

  // Get LOCA-specific information.  Note - in general, LOCA is only used for
  // steady state calculations - there isn't much point in using it
  // for transient.  A common case is one where LOCA is used for the
  // tranOP, but not for the subsequent transient phase.  The
  // locaEnabledFlag should switch from true to false under
  // that scenario, once the transient phase starts.
  solver_state.locaEnabledFlag      = nonlinear_info.locaFlag;
  if (solver_state.locaEnabledFlag) // if no LOCA, these are 0, true, respectively.
  {
    solver_state.continuationStepNumber = nonlinear_info.continuationStep;
    solver_state.firstContinuationParam = nonlinear_info.firstContinuationParam;
  }
  else
  {
    solver_state.continuationStepNumber = 0;
    solver_state.firstContinuationParam = true;
  }
  solver_state.firstSolveComplete     = nonlinear_info.firstSolveComplete;
  // Done with LOCA information.

  // Setup the initialize junctions flag.
  // The initJct flag should only be true if we are at the first Newton step of
  // an initial point in the calculation.  Examples include:
  //     - 1st Newton step of the DCOP initialization for transient (tranOp)
  //     - 1st Newton step of first DC sweep step.
  //     - 1st Newton step of a sweep that has been reset.  That typically
  //         happens if the sweep is multi-dimensional, and the inner loop has
  //         cycled back to the beginning again.
  //
  bool resetFlag =  (solver_state.timeStepNumber_ == 0) || (solver_state.sweepSourceResetFlag);

  // Do this if using LOCA for a DC or tranop calculation.
  if (solver_state.dcopFlag && solver_state.locaEnabledFlag)
  {
    resetFlag = resetFlag && (solver_state.continuationStepNumber==0);
  }

  if (!device_options.disableInitJctFlag)
  {
    solver_state.initJctFlag_ = ((solver_state.dcopFlag) &&
                              (solver_state.newtonIter==0) &&
                               solver_state.firstContinuationParam &&
                               !solver_state.firstSolveComplete && resetFlag);

    // One final check.  See if the "external state" has been set.  If so,
    // check to see if it has the initJctFlag set.  If not, then we probably
    // shouldn't either.  The external state comes from a higher up level
    // in a multi-level newton solve.
    //
    // This should be made more detailed later.
    if (solver_state.externalStateFlag_)
    {
      if (solver_state.newtonIter==0 && solver_state.dcopFlag)
      {
        solver_state.initJctFlag_ = solver_state.externalInitJctFlag_;
      }
    }
  }
  else
  {
    solver_state.initJctFlag_ = false;
  }


  // initFixFlag: try to mimic "MODEINITFIX" of SPICE.  This is set if:
  //   DCOP or TranOP
  //   Not first iteration
  //   Any device not converged
  solver_state.initFixFlag =
    solver_state.dcopFlag
    && !all_devices_converged
    && (solver_state.newtonIter != 0)
    && solver_state.firstContinuationParam
    && !solver_state.firstSolveComplete
    && resetFlag;

  if (solver_state.ltraDevices_)
  {
    if (solver_state.dcopFlag || (solver_state.ltraTimeHistorySize_==0))
    {
      solver_state.ltraTimeIndex_ = 0;
      solver_state.ltraTimeHistorySize_ = 10;
      solver_state.ltraTimePoints_.resize(solver_state.ltraTimeHistorySize_);
    }
  }

  // The first DCOP step of a "double DCOP" simulation is a special case,
  // in which the nonlinear poisson is solved in place of drift-diffusion
  // equations for the PDE devices.  For this initialization problem, the
  // circuit subproblem should not be included.
  if (solver_state.doubleDCOPEnabled && solver_state.dcopFlag && (solver_state.doubleDCOPStep == 0))
  {
    solver_state.twoLevelNewtonCouplingMode       = Nonlinear::INNER_PROBLEM;
  }

  if (DEBUG_DEVICE) 
  {
    solver_state.debugTimeFlag =
      (solver_state.currTime_ >= device_options.debugMinTime
       && solver_state.currTime_ <= device_options.debugMaxTime)
      && (solver_state.timeStepNumber_ >= device_options.debugMinTimestep
          && solver_state.timeStepNumber_ <= device_options.debugMaxTimestep);

    if (isActive(Diag::DEVICE_SOLVER_STATE) && solver_state.debugTimeFlag)
    {
      dout() << solver_state;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : updateTimeInfo
//
// Purpose       : Updates the solver state variables related to time variables.
//
// Special Notes : This function is needed to ensure that breakpoints from expressions
//                 are correct at the beginning of each .STEP iteration.  This function
//                 must be called before any expressions are updated via the setParam
//                 or setParamRandomExpressionTerms functions.
//
//                 Normally, .STEP related updates are handled in the various notify
//                 functions, but that isn't possible for this  function, as the
//                 notify functions happen in the wrong order. The Transient::notify
//                 function is called after the DeviceMgr::notify function.
//
//                 So, this must be called after all the notifies are done, but before
//                 (or at the beginning) of the setParam calls.
//
//                 Also, while the setupSolverInfo function also updates these variables, it
//                 updates many, many things and is overkill for this .STEP/breakpoints issue.
//                 Also, it depends on pointers in the time integrator (particularly the
//                 method pointer) that are often NULL when the setParam functions are
//                 called.
//
// Creator       : Eric R. Keiter, SNL
// Creation Date : 1/25/2023
//-----------------------------------------------------------------------------
bool updateTimeInfo (SolverState & solver_state, const Analysis::AnalysisManager & analysis_manager)
{
  solver_state.currTimeStep_        = analysis_manager.getStepErrorControl().currentTimeStep;                   // system_state.nextTimeStep;
  solver_state.lastTimeStep_        = analysis_manager.getStepErrorControl().lastTimeStep;                      // system_state.currTimeStep;
  solver_state.oldeTimeStep_        = analysis_manager.getStepErrorControl().oldeTimeStep;                      // 
  solver_state.currTime_            = analysis_manager.getStepErrorControl().nextTime;                          // system_state.nextTime;
  solver_state.finalTime_           = analysis_manager.getStepErrorControl().finalTime;                         // system_state.finalTime;
  solver_state.startingTimeStep_    = analysis_manager.getStepErrorControl().startingTimeStep;                  // system_state.startingTimeStep;

  return true;
}

} // namespace Device
} // namespace Xyce


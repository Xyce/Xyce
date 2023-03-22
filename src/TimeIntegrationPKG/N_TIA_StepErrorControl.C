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


//-----------------------------------------------------------------------------
//
// Purpose       : This file contains the functions which define the
//		             time integration stepsize control algorithm.
//
// Special Notes :
//
// Creator       : Buddy Watts
//
// Creation Date : 6/1/00
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <iomanip>
#include <sstream>

#include <N_TIA_StepErrorControl.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ERH_ErrorMgr.h>
#include <N_LOA_Loader.h>
#include <N_PDS_Comm.h>
#include <N_PDS_Manager.h>
#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_TIAParams.h>
#include <N_TIA_WorkingIntegrationMethod.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Functors.h>
#include <N_UTL_MachDepParams.h>
#include <N_UTL_SaveIOSState.h>

namespace Xyce {
namespace TimeIntg {

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::StepErrorControl
// Purpose       : Non-argument constructor.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/06/01
//-----------------------------------------------------------------------------
StepErrorControl::StepErrorControl(
  const std::string &           netlist_filename,
  Analysis::AnalysisManager &   analysis_manager,
  WorkingIntegrationMethod &    working_integration_method,
  const TIAParams &             tia_params)
  : analysisManager_(analysis_manager),
    wimPtr_(working_integration_method),
    netlistFilename_(netlist_filename),
    startingTimeStep(1.0e-10),
    currentTimeStep(1.0e-10),
    lastAttemptedTimeStep(1.0e-10),
    lastTimeStep(1.0e-10),
    oldeTimeStep(1.0e-10),
    minTimeStep(0.0),
    maxTimeStep(0.0),
    maxTimeStepUser(1.0e+99),
    maxTimeStepBP(0.0),
    savedTimeStep(1.0e-10),
    lastTime(0.0),
    currentTime(0.0),
    nextTime(0.0),
    stopTime(0.0),
    initialTime(0.0),
    finalTime(0.0),
    currentTimeStepRatio(0.0),
    currentTimeStepSum(0.0),
    lastTimeStepRatio(0.0),
    lastTimeStepSum(0.0),
    newtonConvergenceStatus(-1),
    nIterations(0),
    numberSuccessiveFailures(0),
    stepAttemptStatus(true),
    previousCallStepSuccessful(false),
    estOverTol_(0.0),
    TimeStepLimitedbyBP(false),
    pauseTime(0.0),
    pauseSetAtZero(false),
    minStepPrecisionFac_(10.0),
    newtonStepReduction_(0.25),
    restartTimeStepScale_(0.005),
    tolAimFac_(0.5),
    breakPointLess_(Util::BreakPoint::defaultTolerance_),
    breakPointEqual_(Util::BreakPoint::defaultTolerance_),
    breakPoints_(),
    currentPauseBP(breakPoints_.end()),
    // define "heuristic" StepSize and Error Control parameters.

    // new-DAE variables:
    currentOrder_(1), // Current order of integration
    minOrder_(1),  // minimum order = max(1,user option minord) - see below.
    maxOrder_(2),  // maximum order = min(2,user option maxord) - see below.
    usedOrder_(1),  // order used in current step (used after currentOrder_ is updated)
    alphas_(-1.0),  // $\alpha_s$ fixed-leading coefficient of this BDF method
    alpha_(6,0.0),  // $\alpha_j(n)=h_n/\psi_j(n)$ coefficient used in local error test
                    // note:   $h_n$ = current step size, n = current time step
    alpha0_(0.0),   // $-\sum_{j=1}^k \alpha_j(n)$ coefficient used in local error test
    cj_ (0.0),      // $-\alpha_s/h_n$ coefficient used in local error test
    ck_ (0.0),      // local error coefficient gamma_[0] = 0; // $\gamma_j(n)=\sum_{l=1}^{j-1}1/\psi_l(n)$ coefficient used to
    sigma_(6,0.0),  // $\sigma_j(n) = \frac{h_n^j(j-1)!}{\psi_1(n)*\cdots *\psi_j(n)}$
    gamma_(6,0.0),  // calculate time derivative of history array for predictor
    beta_(6,0.0),   // coefficients used to evaluate predictor from history array
    psi_(6,0.0),    // $\psi_j(n) = t_n-t_{n-j}$ intermediary variable used to
                    // compute $\beta_j(n)$
    numberOfSteps_(0),   // number of total time integration steps taken
    nef_(0),
    usedStep_(0.0),
    nscsco_(0),
    Ek_(0.0),
    Ekm1_(0.0),
    Ekm2_(0.0),
    Ekp1_(0.0),
    Est_(0.0),
    Tk_(0.0),
    Tkm1_(0.0),
    Tkm2_(0.0),
    Tkp1_(0.0),
    newOrder_(1),
    initialPhase_(true),
    h0_safety_(2.0),
    h0_max_factor_(0.005),  // New value, to match old-DAE.
    //h0_max_factor_(0.0001), // this is the new-DAE equivalent of restartTimeStepScale_(0.005)
    h_phase0_incr_(2.0),
    h_max_inv_(0.0),
    Tkm1_Tk_safety_(2.0),
    Tkp1_Tk_safety_(0.5),
    r_factor_(1.0),
    r_safety_(2.0),
    r_fudge_(0.0001),
    r_min_(0.25),  // r_min_ is the same as the old-DAE variable, minFailStepFac_.
    r_max_(0.9),   // r_max_ is the same as the old-DAE variable, maxFailStepFac_.
    r_hincr_test_(2.0),
    r_hincr_(2.0),
    max_LET_fail_(10),
    maxNumfail_(15),
    reportedPauseBP(false)
{
  setFromTIAParams(tia_params);

  setBreakPoint(Util::BreakPoint(tia_params.finalTime, Xyce::Util::BreakPoint::PAUSE), tia_params.initialTime);
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::~StepErrorControl
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
StepErrorControl::~StepErrorControl()
{}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::setFromTIAParams
// Purpose       : This function copies stuff out of the tiaParams object into
//                 variables that are local to the step error control class.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 09/06/01
//-----------------------------------------------------------------------------
bool StepErrorControl::setFromTIAParams(
  const TIAParams &     tia_params)
{
  startingTimeStep = tia_params.initialTimeStep;
  currentTimeStep  = tia_params.initialTimeStep;
  initialTime      = tia_params.initialTime;
  finalTime        = tia_params.finalTime;
  currentTime      = tia_params.initialTime;
  nextTime         = tia_params.initialTime;
  lastTime         = tia_params.initialTime;

  // if initial time steps are baloney, set then to a default value.
  if (startingTimeStep <= 0.0) startingTimeStep = 1.0e-10;
  if (currentTimeStep <= 0.0)  currentTimeStep  = 1.0e-10;

  if (tia_params.maxTimeStepGiven)
  {
    maxTimeStepUser  = tia_params.maxTimeStep;
    maxTimeStep      = tia_params.maxTimeStep;
  }
  else
  {
    maxTimeStep  = 0.1*(tia_params.finalTime-tia_params.initialTime);
  }

  restartTimeStepScale_ = tia_params.restartTimeStepScale;
  h0_max_factor_ = tia_params.restartTimeStepScale;

  initializeBreakPoints(tia_params.initialOutputTime, tia_params.initialTime, tia_params.finalTime);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::resetAll
//
// Purpose       : This function resets everything so that a transient loop
//                 can be started from the beginning.
//
// Special Notes : This function was needed for the .STEP capability.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 11/04/03
//-----------------------------------------------------------------------------
bool StepErrorControl::resetAll(const TIAParams & tia_params)
{
  startingTimeStep = tia_params.initialTimeStep;
  currentTimeStep  = tia_params.initialTimeStep;
  initialTime      = tia_params.initialTime;
  finalTime        = tia_params.finalTime;
  currentTime      = tia_params.initialTime;
  lastTime         = tia_params.initialTime;
  nextTime         = tia_params.initialTime;

  breakPointLess_.tolerance_ = Util::BreakPoint::defaultTolerance_;
  breakPointEqual_.tolerance_ = Util::BreakPoint::defaultTolerance_;

  // if initial time steps are baloney, set then to a default value.
  if (startingTimeStep <= 0.0) startingTimeStep = 1.0e-10;
  if (currentTimeStep <= 0.0)  currentTimeStep  = 1.0e-10;

  if (tia_params.maxTimeStepGiven)
  {
    maxTimeStepUser  = tia_params.maxTimeStep;
    maxTimeStep      = tia_params.maxTimeStep;
  }
  else
  {
    maxTimeStep  = 0.1* (tia_params.finalTime - tia_params.initialTime);
  }

  restartTimeStepScale_ = tia_params.restartTimeStepScale;

  initializeBreakPoints(tia_params.initialOutputTime, tia_params.initialTime, tia_params.finalTime);

  pauseSetAtZero = false;
  pauseTime = 0.0;

  lastTimeStep     = tia_params.initialTimeStep;
  lastAttemptedTimeStep = tia_params.initialTimeStep;

  currentTimeStepRatio = 1.0;
  currentTimeStepSum   = 2.0*currentTimeStep;

  lastTimeStep      = currentTimeStep;
  lastTimeStepRatio = currentTimeStepRatio;
  lastTimeStepSum   = currentTimeStepSum;

  newtonConvergenceStatus = -1;
  numberSuccessiveFailures = 0;
  stepAttemptStatus        = true;

  minTimeStep = 0.0;
  estOverTol_ = 0.0;

  // need to set a pause breakpoint at the final time.
  // This is done inside of the "initializeBreakPoints" function, so this is redundant
  setBreakPoint(Util::BreakPoint(tia_params.finalTime ,Xyce::Util::BreakPoint::PAUSE), tia_params.initialTime);

  if (DEBUG_TIME && isActive(Diag::TIME_BREAKPOINTS))
  {
    Xyce::dout() << "  after resetAll:" << std::endl;
    printBreakPoints(Xyce::dout());
    if ( !(breakPoints_.empty()) && currentPauseBP != breakPoints_.end() ) // if not valid don't do this
    {
      Xyce::dout() <<"  currentPauseBP = " <<  currentPauseBP->value() << std::endl
                   << Xyce::section_divider << std::endl;
    }
    else
    {
      Xyce::dout() <<"  currentPauseBP not valid" << std::endl
                   << Xyce::section_divider << std::endl;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::getEstOverTol
// Purpose       : This function lets the controlling class (a transient
//                 analysis) get the estimated error over tol from the
//                 last step
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and MEMS Modeling
// Creation Date : 01/23/09
//-----------------------------------------------------------------------------
double StepErrorControl::getEstOverTol() const
{
  return estOverTol_;
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::setTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 08/09/09
//-----------------------------------------------------------------------------
void StepErrorControl::setTimeStep(double newTimeStep)
{
  newTimeStep = std::max(newTimeStep, minTimeStep);
  newTimeStep = std::min(newTimeStep, maxTimeStep);

  double nextTimePt = currentTime + newTimeStep;

  if (nextTimePt > stopTime)
  {
    nextTimePt  = stopTime;
    newTimeStep = stopTime - currentTime;
    TimeStepLimitedbyBP = true;
  }

  nextTime = nextTimePt;

  currentTimeStepRatio = newTimeStep/lastTimeStep;
  currentTimeStepSum   = newTimeStep + lastTimeStep;

  currentTimeStep = newTimeStep;
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::updateStopTime
// Purpose       : The "stop time" is either the next discontinuity point,
//                 a pause point, or the final time, whichever comes first.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel ComputationalSciences.
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
void StepErrorControl::updateStopTime(
  Parallel::Machine     comm,
  bool                  breakpoints_enabled,
  double                initial_time,
  bool                  min_time_steps_breakpoint_given,
  double                min_time_steps_breakpoint)
{
  double oldStopTime = stopTime;
  double diffStopTime = 0.0;

  if (breakpoints_enabled)
  {
    // Find the first breakpoint equal to or larger than the
    // current time.
    BreakPointVector::iterator itBP = std::upper_bound(breakPoints_.begin(), breakPoints_.end(), currentTime, breakPointLess_);

    stopTime =  std::min(finalTime, itBP->value());

    // The breakpoint could be a pause breakpoint, in which case we might
    // need to update the pauseTime:
    if (itBP->bptype() == Xyce::Util::BreakPoint::PAUSE)
    {
      updatePauseTime(*itBP, initial_time);
    }

    // this seems like a weird thing to do - we just:
    //
    // stopTIme = nearest breakpoint
    //
    //
    if ( !(breakPoints_.empty()) && currentPauseBP != breakPoints_.end() ) // if not valid don't do this
    {
      stopTime =  std::min(currentPauseBP->value(), stopTime);
    }

    // if this is a breakpoint step, make sure the new stop
    // time is for the next breakpoint, not the current one.
    // This check is neccessary because of roundoff error.

    diffStopTime = fabs(stopTime-oldStopTime);
    if (diffStopTime < breakPointLess_.tolerance_ &&
       analysisManager_.getBeginningIntegrationFlag() &&
       stopTime != pauseTime &&
       stopTime != finalTime )
    {
      ++itBP;
      stopTime = itBP->value();
    }

    Parallel::AllReduce(comm, MPI_MIN, &stopTime, 1);
  }
  else
  {
    stopTime = std::min(pauseTime, finalTime);
  }

  if ( analysisManager_.getBeginningIntegrationFlag())
  {
    double time_to_stop = stopTime - currentTime;
    if (min_time_steps_breakpoint_given && (min_time_steps_breakpoint > 0) )
    {
      maxTimeStepBP = time_to_stop/min_time_steps_breakpoint;
    }
  }


  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << std::endl
                 << "  stopTime    = " <<  stopTime << std::endl
                 << "  pauseTime    = " <<  pauseTime << std::endl
                 << "  currentTime = " <<  currentTime << std::endl
                 << "  oldStopTime = " <<  oldStopTime << std::endl
                 << "  finalTime   = " <<  finalTime << std::endl
                 << "  maxTimeStepBP = " <<  maxTimeStepBP << std::endl
                 << "  beginningIntegration = " <<  analysisManager_.getBeginningIntegrationFlag() << std::endl
                 << Xyce::section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::evaluateStepError
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel ComputationalSciences.
// Creation Date : 1/28/07
//-----------------------------------------------------------------------------
void StepErrorControl::evaluateStepError(
  const Loader::Loader &        loader,
  const TIAParams &             tia_params)
{
  bool step_attempt_status( newtonConvergenceStatus >= 0);
  bool sAStatus(false);
  bool errorOptionStatus(true);
  bool testTimeIntegrationError(false);

  // If we are running with constant step size, or are on the first pass
  // through the transient loop, only base success on the Newton loop.
  if (tia_params.newBPStepping)
  {
    if (currentTime == tia_params.initialTime)
    {
      testTimeIntegrationError = (analysisManager_.getStepNumber() >= 1 && !analysisManager_.getBeginningIntegrationFlag());
    }
    else
    {
      testTimeIntegrationError = (analysisManager_.getStepNumber() >= 1);
    }
  }
  else
  {
    testTimeIntegrationError = (analysisManager_.getStepNumber() >= 1 && !analysisManager_.getBeginningIntegrationFlag());
  }

  if (tia_params.testFirstStep)
  {
    testTimeIntegrationError = true;
  }

  // if the step status is already false, don't do any more work.
  if (!step_attempt_status)
  {
    testTimeIntegrationError = false;
  }


  if (testTimeIntegrationError)
  {
    // Needed for 2-level Solves:
    loader.getInnerLoopErrorSums(analysisManager_.getDataStore()->innerErrorInfoVec);

    estOverTol_ = wimPtr_.computeErrorEstimate();

    if (estOverTol_ <= tia_params.errTolAcceptance)
    {
      sAStatus = true;
    }
    else
    {
      sAStatus = false;
    }

    if (tia_params.timestepsReversal == true)
    {
      if (nIterations <= tia_params.NLmax)
        errorOptionStatus = true;
      else
        errorOptionStatus = false;
    }

    if (VERBOSE_TIME && tia_params.errorAnalysisOption == TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES)
    {
      Xyce::dout() << "ERROROPTION=1:  DOREJECTSTEP = ";
      if (tia_params.timestepsReversal == true)
      {
        Xyce::dout() << "1" << std::endl;
      }
      else
      {
        Xyce::dout() << "0" << std::endl;
      }
    }

    if ( tia_params.minTimeStepGiven && (currentTimeStep < tia_params.minTimeStep) )
    {
      // This step has dropped under the user specified min time step, so only
      // test if the solver converged to accept the step.
      // don't do the step_attempt_status && sAStatus;
      if (DEBUG_TIME && isActive(Diag::TIME_ERROR))
      {
        Xyce::dout() << "Trying to skip time integrator error checks: " << currentTimeStep
                     << " newton status " << step_attempt_status << std::endl;
      }
    }
//    else if (!tia_params.constantTimeStepFlag && ((tia_params.errorAnalysisOption == LOCAL_TRUNCATED_ESTIMATES)) )
//    {
//      step_attempt_status = step_attempt_status && sAStatus;
//    }
    else if (!tia_params.constantTimeStepFlag)
    {
      if (tia_params.errorAnalysisOption == TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES)
        step_attempt_status = step_attempt_status && errorOptionStatus;
      else
        step_attempt_status = step_attempt_status && sAStatus;
    }
  }

  if (DEBUG_TIME && isActive(Diag::TIME_ERROR))
  {
    integrationStepReport(Xyce::dout(), step_attempt_status, sAStatus, testTimeIntegrationError, tia_params);
  }
  else if (VERBOSE_TIME)
  {
    terseIntegrationStepReport(Xyce::dout(), step_attempt_status, sAStatus, testTimeIntegrationError, tia_params);
  }

  // Now that the status has been completely determined,
  // set the class variable for step attempt
  stepAttemptStatus = step_attempt_status;
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::terseIntegrationStepReport_
// Purpose       : This gives a one-line description of the step accept/reject.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/27/07
//-----------------------------------------------------------------------------
void StepErrorControl::terseIntegrationStepReport(
  std::ostream &        os,
  bool                  step_attempt_status,
  bool                  sAStatus,
  bool                  testedError,
  const TIAParams &     tia_params)
{
  os << (DEBUG_TIME ? netlistFilename_ : "")
     << "  STEP STATUS: " << (step_attempt_status ? " success" : " fail")
     << "  Newton: " << newtonConvergenceStatus
     << "   estOverTol: " << estOverTol_ << (testedError && !tia_params.constantTimeStepFlag ? "" : " (not used for this step)") << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::integrationStepReport_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void StepErrorControl::integrationStepReport(std::ostream &os, bool step_attempt_status, bool sAStatus, bool testedError, const TIAParams &tia_params)
{
  if (isActive(Diag::TIME_PARAMETERS))
  {
    os << "\n estOverTol      = " <<  estOverTol_ << std::endl
       << "  error tolerance = " << tia_params.errTolAcceptance << std::endl
       << std::endl
       << "\nSTEP ATTEMPT STATUS:" << std::endl
       << "NOTE:" << std::endl;

    if (!tia_params.constantTimeStepFlag &&
        analysisManager_.getStepNumber() >= 1 &&
        !analysisManager_.getBeginningIntegrationFlag())
    {
      os << "  We are running in variable stepsize mode " << std::endl
         << "  and we have NOT just passed a breakpoint.  As such " << std::endl
         << "  for an integration step to succeed the " << std::endl
         << "  nonlinear solver must succeed AND the predictor" << std::endl
         << "  and corrector need to be close within a tolerance." << std::endl;

      if (tia_params.errorAnalysisOption == TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES)
      {
        os << "ADDENDUM:  This is with erroption=1 so predictor-corrector is ignored for step error control." << std::endl;
      }
    }
    else
    {
      os << "  We are either running constant stepsize " << std::endl
         << "  or we just passed a breakpoint.  As such " << std::endl
         << "  the only criteria we use in accepting/rejecting" << std::endl
         << "  an integration step is the nonlinear solver" << std::endl
         << "  success/failure." << std::endl;
    }

    if (step_attempt_status)
    {
      os << "\n  This has been a successful step:" << std::endl;
    }
    else
    {
      os << "\n  This has NOT been a successful step:" << std::endl;
    }

    if ( newtonConvergenceStatus > 0)
    {
      os << "    - Newton solver succeded with return code " << newtonConvergenceStatus << std::endl << std::endl;
    }
    else
    {
      os << "    - Newton solver failed with return code " << newtonConvergenceStatus << std::endl;
    }

    if (testedError)
    {
      if (!tia_params.constantTimeStepFlag)
      {
        if (sAStatus)
        {
          os << "   - predictor vs. corrector analysis succeeded." << std::endl;
        }
        else
        {
          os << "   - predictor vs. corrector analysis failed." << std::endl;
        }

        os << "     (compare estOverTol with error tolerance above.)" << std::endl;
      }
      else
      {
        os << "If we had been using it <<  " << std::endl;

        if (sAStatus)
        {
          os << "   - predictor vs. corrector analysis would have succeeded." << std::endl;
        }
        else
        {
          os << "   - predictor vs. corrector analysis would have failed." << std::endl;
        }

        os << "     (compare estOverTol with error tolerance above.)" << std::endl;
      }
    }
    else
    {
      os << "  predictor vs. corrector was not tested" << std::endl;
    }

    os << Xyce::section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::initializeBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel ComputationalSciences.
// Creation Date : 06/11/01
//-----------------------------------------------------------------------------
bool StepErrorControl::initializeBreakPoints(
  double        initial_output_time,
  double        initial_time,
  double        final_time)
{
  bool bsuccess = true;

  breakPoints_.clear ();
  currentPauseBP = breakPoints_.end();

  // first breakpoint is the start time, last one is the final time.
  setBreakPoint(initialTime, initial_time);

  // if initial_output_time is nonzero, then make it a breakpoint.
  if (initial_output_time > initialTime && initial_output_time < finalTime)
  {
    setBreakPoint(initial_output_time, initial_time);
  }

  // The final time needs to be the very last breakpoint.
  setBreakPoint(Util::BreakPoint(final_time, Xyce::Util::BreakPoint::PAUSE), initial_time);

  if (analysisManager_.getBlockAnalysisFlag())
  {
    BreakPointVector::iterator lastBP = breakPoints_.end();
    --lastBP;
    updatePauseTime(*lastBP, initialTime);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::updateBreakPoints
// Purpose       : Requests dynamic breakpoint information from the
//                 loader.  Adds, subtracts from the breakpoints array.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel ComputationalSciences.
// Creation Date : 06/11/01
//-----------------------------------------------------------------------------
bool StepErrorControl::updateBreakPoints(
  const Loader::Loader &        loader, 
  double                        initial_time)
{
  bool bsuccess = true;

  if (DEBUG_TIME && isActive(Diag::TIME_BREAKPOINTS))
  {
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  StepErrorControl::updateBreakPoints.  time = " <<  currentTime << std::endl
                 << std::endl;
  }

  double bpTol = 2.0 * minTimeStep;
  //double bpTol = std::max( Xyce::Util::MachineDependentParams::MachinePrecision(), (2.0*minTimeStep) );
  breakPointLess_.tolerance_ = bpTol;
  breakPointEqual_.tolerance_ = bpTol;

  std::vector<Util::BreakPoint> tmpBP; tmpBP.clear ();
  std::vector<Util::BreakPoint> tmpPauseBP; tmpPauseBP.clear ();
  loader.getBreakPoints(tmpBP,tmpPauseBP);

  // save a copy of the currentPauseBP, as all the sorting, etc that is about to happen may break it.
  Util::BreakPoint savedCurrentPauseBP(finalTime, Xyce::Util::BreakPoint::PAUSE);
  if ( !(breakPoints_.empty()) && currentPauseBP != breakPoints_.end() ) // if not valid don't do this
  {
    savedCurrentPauseBP = *currentPauseBP;
  }

  bool BPchanged=false;

  // Add new breakpoints to the set:
  if ( !(tmpBP.empty()) )
  {
    // simple breakpoints:
    // sort,unique and eliminate too early points in the the tmpBP vector
    std::sort ( tmpBP.begin(), tmpBP.end(), breakPointLess_ );
    std::vector<Util::BreakPoint>::iterator it = std::unique ( tmpBP.begin(), tmpBP.end(), breakPointEqual_ );
    tmpBP.resize( std::distance (tmpBP.begin(), it ));

    std::vector<Util::BreakPoint>::iterator iter;
    std::vector<Util::BreakPoint>::iterator first = std::lower_bound(tmpBP.begin(), tmpBP.end(), lastTime, breakPointLess_);    
    std::vector<Util::BreakPoint>::iterator last  = std::upper_bound(tmpBP.begin(), tmpBP.end(), finalTime, breakPointLess_);    

    // add the contents of tmpBP to the master breakPoints_ container.
    // This will require a subsequent std::sort  and std::unique
    breakPoints_.insert(breakPoints_.end(), first,last);
    BPchanged=true;
  }

  if ( !(tmpPauseBP.empty()) )
  {
    reportedPauseBP = true;
    // pause breakpoints:
    // sort,unique and eliminate too early points in the tmpPauseBP vector
    std::sort ( tmpPauseBP.begin(), tmpPauseBP.end(), breakPointLess_ );
    std::vector<Util::BreakPoint>::iterator it = std::unique ( tmpPauseBP.begin(), tmpPauseBP.end(), breakPointEqual_ );
    tmpPauseBP.resize( std::distance (tmpPauseBP.begin(), it ));

    // add the contents of tmpPauseBP to the master breakPoints_ container.
    std::vector<Util::BreakPoint>::iterator iter;
    std::vector<Util::BreakPoint>::iterator first = std::lower_bound(tmpPauseBP.begin(), tmpPauseBP.end(), lastTime, breakPointLess_);    
    std::vector<Util::BreakPoint>::iterator last  = std::upper_bound(tmpPauseBP.begin(), tmpPauseBP.end(), finalTime, breakPointLess_);    
    for (iter=first; iter!=last; ++iter)
    {
      if (iter->value() < finalTime && iter->value() > lastTime)
      {
        setBreakPoint(*iter, initial_time);
      }
    }
    BPchanged=true;
  }

  if (BPchanged)
  {
    std::sort ( breakPoints_.begin(), breakPoints_.end(), breakPointLess_ );
    std::vector<Util::BreakPoint>::iterator it = std::unique ( breakPoints_.begin(), breakPoints_.end(), breakPointEqual_ );
    breakPoints_.resize( std::distance (breakPoints_.begin(), it ));
  }

  {
  // Remove breakpoints which are now obsolete (old):
  BreakPointVector::iterator itBP;
  itBP = std::lower_bound(breakPoints_.begin(), breakPoints_.end(), lastTime, breakPointLess_);

  if (itBP != breakPoints_.begin())
  {
    breakPoints_.erase(breakPoints_.begin(),itBP);
    BPchanged=true;
  }
  }

  // some of the iterators might have gotten mangled, so check and fix.
  if (BPchanged)
  {
    doubleCheckEndBreakPoint();
    currentPauseBP = std::find(breakPoints_.begin(), breakPoints_.end(), savedCurrentPauseBP);
  }

  if (DEBUG_TIME && isActive(Diag::TIME_BREAKPOINTS))
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() << "  breakPoints_ vector container after sort/unique/erase:" << std::endl;

    int i=0;
    BreakPointVector::iterator itBP;
    BreakPointVector::iterator itBP_2;
    for (itBP=breakPoints_.begin();itBP!=breakPoints_.end();++i,++itBP)
    {
      if (i==0)
      {
        Xyce::dout() << i << " " << itBP->value() << "  type=" << itBP->bptype() << std::endl;
      }
      else
      {
        Xyce::dout() << i << " " << itBP->value() << "  type=" << itBP->bptype() << " diff=" << (itBP->value()-itBP_2->value()) << std::endl;
      }

      itBP_2 = itBP;
    }

    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::updateMaxTimeStep
// Purpose       : Requests dynamic time step information from the loader.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
bool StepErrorControl::updateMaxTimeStep(
  Parallel::Machine     comm,
  Loader::Loader &      loader,
  const TIAParams &     tia_params,
  double                suggestedMaxTimeStep)
{
  bool bsuccess = true;

  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << Xyce::section_divider << std::endl
                 << "  StepErrorControl::updateMaxTimeStep" << std::endl;
  }

  double maxDevStep = 1.0e+99;
  if (tia_params.useDeviceTimeStepMaxFlag)
  {
    maxDevStep = loader.getMaxTimeStepSize();
  }

  if (tia_params.maxTimeStepGiven || tia_params.delmaxGiven)
  {
    maxTimeStep = std::min(tia_params.maxTimeStep, tia_params.delmax);
  }
  else
  {
    maxTimeStep  = 0.1*(tia_params.finalTime-tia_params.initialTime);
  }

  // if the default arg is not zero, then a suggested max time step was
  // passed in.  Test if it is feasible to use that at this time
  if( suggestedMaxTimeStep > 0.0 )
  {
    maxTimeStep = std::min( maxTimeStep, suggestedMaxTimeStep );
  }

  if ((maxTimeStepBP > 0.0) && (maxTimeStep > maxTimeStepBP))
  {
    maxTimeStep = maxTimeStepBP;
  }

  if (maxDevStep > 0.0)
  {
    maxTimeStep = std::min(maxTimeStep, maxDevStep);
  }

  if (tia_params.maxTimeStepGiven)
  {
    maxTimeStep = std::min(maxTimeStep, maxTimeStepUser);
  }

  Parallel::AllReduce(comm, MPI_MIN, &maxTimeStep, 1);

  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    if (!tia_params.maxTimeStepGiven)
    {
      Xyce::dout() << "  User did not specify a maximum time step." << std::endl;
    }
    else
    {
      Xyce::dout() << "  User specified a maximum time step. = " <<  maxTimeStepUser << std::endl;
    }

    Xyce::dout() << "  maxDevStep  = " <<  maxDevStep << std::endl
                 << "  maxTimeStep = " <<  maxTimeStep << std::endl
                 << Xyce::section_divider << std::endl;
  }

  if(maxTimeStep<=0.0)
  {
    Report::DevelFatal0() << "Maximum Time step is invalid!";
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::updateMinTimeStep
// Purpose       : Sets the minimum time step based on machine precision.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel ComputationalSciences.
// Creation Date : 08/02/01
//-----------------------------------------------------------------------------
bool StepErrorControl::updateMinTimeStep()
{
  bool bsuccess = true;

  minTimeStep = currentTime*minStepPrecisionFac_*Util::MachineDependentParams::MachinePrecision();
//  minTimeStep = minStepPrecisionFac_*Util::MachineDependentParams::MachinePrecision();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::setBreakPoint
// Purpose       : public method to set individual breakpoint
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/30/04
//-----------------------------------------------------------------------------
void StepErrorControl::setBreakPoint(const Util::BreakPoint & breakpoint, double initial_time)
{
  std::vector<Util::BreakPoint>::iterator found = std::find(breakPoints_.begin(), breakPoints_.end(), breakpoint);
  if (found != breakPoints_.end()) 
  {
    // if found, and the passed BP is a PAUSE BP, then replace (to change from SIMPLE to PAUSE)
    if (breakpoint.bptype() == Xyce::Util::BreakPoint::PAUSE)
    {
      *found = breakpoint;
      updatePauseTime(breakpoint, initial_time);
    }
  }
  else // if not found, then add it in
  {
    // save a copy of the currentPauseBP, as all the sorting, etc that is about to happen may break it.
    Util::BreakPoint savedCurrentPauseBP(finalTime, Xyce::Util::BreakPoint::PAUSE);
    if ( !(breakPoints_.empty()) && currentPauseBP != breakPoints_.end() ) // if not valid don't do this
    {
      savedCurrentPauseBP = *currentPauseBP;
    }

    breakPoints_.push_back(breakpoint);
    std::sort ( breakPoints_.begin(), breakPoints_.end(), breakPointLess_ );
    std::vector<Util::BreakPoint>::iterator it = std::unique ( breakPoints_.begin(), breakPoints_.end(), breakPointEqual_ );
    breakPoints_.resize( std::distance (breakPoints_.begin(), it ));

    // The iterator pointing to the currentPause time was mangled by 
    // all the above work, so ind it again and restore
    currentPauseBP = std::find(breakPoints_.begin(), breakPoints_.end(), savedCurrentPauseBP);

    // if the new breakpoint is a pause breakpoint, we need to update the pause time
    if (breakpoint.bptype() == Xyce::Util::BreakPoint::PAUSE)
    {
      updatePauseTime(breakpoint, initial_time);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::setBreakPoint
// Purpose       : public method to set individual SIMPLE breakpoint
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/30/04
//-----------------------------------------------------------------------------
void StepErrorControl::setBreakPoint(double time)
{
  std::vector<Util::BreakPoint>::iterator found = std::find(breakPoints_.begin(), breakPoints_.end(), time);

  if (found == breakPoints_.end()) // if found, do nothing, we already have this one
  {
    // save a copy of the currentPauseBP, as all the sorting, etc that is about to happen may break it.
    Util::BreakPoint savedCurrentPauseBP(finalTime, Xyce::Util::BreakPoint::PAUSE);
    if ( !(breakPoints_.empty()) && currentPauseBP != breakPoints_.end() ) // if not valid don't do this
    {
      savedCurrentPauseBP = *currentPauseBP;
    }

    breakPoints_.push_back(time);
    std::sort ( breakPoints_.begin(), breakPoints_.end(), breakPointLess_ );
    std::vector<Util::BreakPoint>::iterator it = std::unique ( breakPoints_.begin(), breakPoints_.end(), breakPointEqual_ );
    breakPoints_.resize( std::distance (breakPoints_.begin(), it ));

    // now, if the iterator pointing to the currentPause time was mangled by 
    // all the above work, find it again and restore
    currentPauseBP = std::find(breakPoints_.begin(), breakPoints_.end(), savedCurrentPauseBP);
  }
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::doubleCheckEndBreakPoint
// Purpose       : check that last iterator is at the final time and is a PAUSE BP
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/18/2019
//-----------------------------------------------------------------------------
void StepErrorControl::doubleCheckEndBreakPoint()
{
  BreakPointVector::iterator itBP;
  Util::BreakPoint candidateBP(finalTime, Xyce::Util::BreakPoint::PAUSE);

  if (!breakPoints_.empty())
  {
    itBP = breakPoints_.end()-1;
    if ( breakPointEqual_( *itBP , candidateBP ) )
    {
      itBP->setType(Xyce::Util::BreakPoint::PAUSE);
    }
    else if ( breakPointLess_ (*itBP, candidateBP) )
    {
      breakPoints_.push_back( Util::BreakPoint(finalTime, Xyce::Util::BreakPoint::PAUSE) );
    }
    else 
    {
      // if this happens, then the lower_bound failed
    }
  }
  else
  {
    breakPoints_.push_back( Util::BreakPoint(finalTime, Xyce::Util::BreakPoint::PAUSE) );
    currentPauseBP = breakPoints_.end()-1;
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::updatePauseTime
// Purpose       : private method to recalculate pause time
//
// Special Notes : If the breakpoint being passed in is a pause breakpoint, 
//                 and it is earlier than the current value for the puase 
//                 time, then the pauseTime variable will be reset.
//
// Scope         : private
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/30/04
//-----------------------------------------------------------------------------
void StepErrorControl::updatePauseTime(
  Util::BreakPoint      breakpoint,
  double                initial_time)
{
  // gotta handle case where pauseTime is still at its initial value of
  // 0.0, or we already passed the last pause time!
  // But we mustn't reset it if it's equal to the current time, because
  // that means we need to stop NOW and would overwrite that.
  //
  // If a pause break point is specifically set at 0, then
  // we shouldn't ignore that here.  So set the pauseSetAtZero flag
  // here if needed.
  //
  // ERK: type>0 means breakpoint is "not simple", ie PAUSE breakpoint.
  if ((breakpoint.bptype() == Xyce::Util::BreakPoint::PAUSE) && (breakpoint.value() == 0.0))
  {
    pauseSetAtZero = true;
  }

  if (pauseTime < currentTime || ((pauseTime == initial_time) && !pauseSetAtZero))
  {
    pauseTime = breakpoint.value();
  }
  else
  {
    pauseTime = std::min(pauseTime, breakpoint.value());
  }

  // If we used this breakpoint for the pause time, save the iterator to this bp in
  // the list so we can use it later.
  if (pauseTime == breakpoint.value())
  {
    currentPauseBP = std::find(breakPoints_.begin(), breakPoints_.end(), breakpoint);

    if (DEBUG_TIME && isActive(Diag::TIME_BREAKPOINTS))
    {
      Xyce::dout() << "\n" << netlistFilename_;
      if(currentPauseBP != breakPoints_.end()) 
      {
        Xyce::dout() << "  UPDATING PAUSE TIME TO " << currentPauseBP->value();
      }
      else
      {
        Xyce::dout() << "  UPDATING PAUSE TIME FAILED TO UPDATE ";
      }

      Xyce::dout() << " encountered breakpoint " << breakpoint.value() << " current time  is "
                   << currentTime << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::simulationPaused
// Purpose       : public method to clear out breakpoint list and reset
//                 pause time when pause time is reached
// Special Notes : Need this method because the prior values get in the way
//                 when we resume.
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/30/04
//-----------------------------------------------------------------------------
void StepErrorControl::simulationPaused(double initial_time)
{
  BreakPointVector::iterator lb = std::lower_bound(breakPoints_.begin(), breakPoints_.end(), currentTime, breakPointLess_);
  if (lb != breakPoints_.end() )
  {
    breakPoints_.erase(breakPoints_.begin(), lb); 
  }
  currentPauseBP = breakPoints_.end(); // make this invalid
  pauseTime = initial_time;          // unset this
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::printBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/17/05
//-----------------------------------------------------------------------------
void StepErrorControl::printBreakPoints (std::ostream & os) const
{
  BreakPointVector::const_iterator itBP;
  BreakPointVector::const_iterator itBP2;
  BreakPointVector::const_iterator firstBP = breakPoints_.begin();
  BreakPointVector::const_iterator lastBP  = breakPoints_.end();

  char tmp[128];

  int i;
  for (i=0, itBP=firstBP;itBP!=lastBP;++i,++itBP)
  {
    if (i==0)
      sprintf(tmp,"%4d %16.8e  type=%d",i,itBP->value(),itBP->bptype());
    else
      sprintf(tmp,"%4d %16.8e type=%d diff=%16.8e", i, itBP->value(),
	      itBP->bptype(),(itBP->value()-itBP2->value()));

    os << std::string(tmp);
    itBP2 = itBP;
  }
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::restartDataSize
// Purpose       :
// Special Notes : This gives the *total* size:  both the base
//                 StepErrorControl and the derived
//                 StepErrorControlDAE, summed together.
//
//                 Don't sum them 2x!
//
//                 ERK:  6/20/2010:
//                 There are no longer 2 distinct classes.  The DAE version is
//                 now part of the original, as one single class.  The two
//                 blocks of code in this function correspond to the original
//                 class (first block), and the DAE class (second block).
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 07/27/06
//-----------------------------------------------------------------------------
int StepErrorControl::getRestartDataSize( bool pack )
{
  int totalSize = 0;

  // original set of vars:
  int numdoubles =  22;
  int numints = 9;
  int count = sizeof(double)*numdoubles;
  count += sizeof(int)*numints;

  // Must include the bp type now, not just the value
  count += sizeof(Util::BreakPoint)*breakPoints_.size();

  //overestimate buffer size for unpacked data
  // assumes there are fewer than 100 possible breakpoints types (i.e.
  // because the type can be represented by two ascii characters)
  if( !pack )
  {                  // 34
    count = 24*(numdoubles + numints + breakPoints_.size()) + 2*breakPoints_.size();
  }

  int baseClassSize = count;

  // another set of vars: (newDAE)
  numdoubles = 57;
  numints = 10;

  totalSize = baseClassSize;
  totalSize += sizeof(double) * numdoubles;
  totalSize += sizeof(int) * numints;

  //overestimate buffer size for unpacked data
  if ( !pack )
  {
    totalSize = baseClassSize + 24*(numdoubles+numints);
  }

  return totalSize;
}
//-----------------------------------------------------------------------------
// Function      : StepErrorControl::dumpRestartData
// Purpose       :
// Special Notes : ERK:  6/20/2010:
//                 There are no longer 2 distinct classes.  The DAE version is
//                 now part of the original, as one single class.  The two
//                 blocks of code in this function correspond to the original
//                 class (first block), and the DAE class (second block).
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 7/28/06
//-----------------------------------------------------------------------------
bool StepErrorControl::dumpRestartData(
  char *                    buf,
  int                       bsize,
  int &                     pos,
  Parallel::Communicator *  comm,
  bool                      pack)
{

  // Set this variable up for later.  Note that pos means different things
  // for packed vs. unpacked.  For unpacked, it is the current index into
  // the buf array.
  int newPos = pos + getRestartDataSize(false);

  // if unpacked, initialize the buf array before calling the base class
  // function.  At this point pos is zero, probably.  The derived DAE
  // dataSize will be the size of the entire buf array, including the
  // base class size.
  if ( !pack )
  {
    for( int i = pos; i < (newPos); ++i) buf[i] = ' ';
  }

  if (DEBUG_RESTART)
  {
    Xyce::dout() << "TIA Restart Data DUMP!  " << netlistFilename_ << "\n"
                 << Xyce::section_divider << std::endl
                 << "startingTimeStep: " << startingTimeStep << std::endl
                 << "currentTimeStep: " << currentTimeStep << std::endl
                 << "lastAttemptedTimeStep: " << lastAttemptedTimeStep << std::endl
                 << "lastTimeStep: " << lastTimeStep << std::endl
                 << "minTimeStep: " << minTimeStep << std::endl
                 << "maxTimeStep: " << maxTimeStep << std::endl
                 << "maxTimeStepUser: " << maxTimeStepUser << std::endl
                 << "lastTime: " << lastTime << std::endl
                 << "currentTime: " << currentTime << std::endl
                 << "nextTime: " << nextTime << std::endl
                 << "initialTime: " << initialTime << std::endl
                 << "estOverTol_: " << estOverTol_ << std::endl
                 << "breakpts: ";

    for (BreakPointVector::iterator iterSD = breakPoints_.begin(); iterSD != breakPoints_.end(); ++iterSD )
      Xyce::dout() << iterSD->value() << " ";

    Xyce::dout() << std::endl
                 << "integMethod: " << analysisManager_.getIntegrationMethod() << std::endl
                 << "stepNumber: " << analysisManager_.getStepNumber() << std::endl
                 << "transStepNumber: " << analysisManager_.getTranStepNumber() << std::endl
                 << "breakPointRestartNumber: " << analysisManager_.breakPointRestartStep << std::endl
                 << Xyce::section_divider << std::endl << std::endl;
  }

  if( pack )
  {
    comm->pack( &startingTimeStep, 1, buf, bsize, pos );
    comm->pack( &currentTimeStep, 1, buf, bsize, pos );
    comm->pack( &lastAttemptedTimeStep, 1, buf, bsize, pos );
    comm->pack( &lastTimeStep, 1, buf, bsize, pos );
    comm->pack( &minTimeStep, 1, buf, bsize, pos );
    comm->pack( &maxTimeStep, 1, buf, bsize, pos );
    comm->pack( &maxTimeStepUser, 1, buf, bsize, pos );
    comm->pack( &lastTime, 1, buf, bsize, pos );
    comm->pack( &currentTime, 1, buf, bsize, pos );
    comm->pack( &nextTime, 1, buf, bsize, pos );
    comm->pack( &initialTime, 1, buf, bsize, pos );
    comm->pack( &currentTimeStepRatio, 1, buf, bsize, pos );
    comm->pack( &currentTimeStepSum, 1, buf, bsize, pos );
    comm->pack( &lastTimeStepRatio, 1, buf, bsize, pos );
    comm->pack( &lastTimeStepSum, 1, buf, bsize, pos );
    comm->pack( &newtonConvergenceStatus, 1, buf, bsize, pos );
    comm->pack( &numberSuccessiveFailures, 1, buf, bsize, pos );
    int flag = stepAttemptStatus;
    comm->pack( &flag, 1, buf, bsize, pos );
    comm->pack( &minStepPrecisionFac_, 1, buf, bsize, pos );
    comm->pack( &newtonStepReduction_, 1, buf, bsize, pos );
    comm->pack( &tolAimFac_, 1, buf, bsize, pos );
    comm->pack( &estOverTol_, 1, buf, bsize, pos );
    // Subtract one, because we won't write out the pause breakpoint at the
    // final time
    int size = breakPoints_.size() -1 ;
    BreakPointVector::iterator bpStart = breakPoints_.begin();
    BreakPointVector::iterator bpEnd = breakPoints_.end();
    comm->pack( &size, 1, buf, bsize, pos );

    double val;
    int bptype;
    {
    int i=0;
    for(BreakPointVector::iterator iterSD = bpStart; iterSD != bpEnd; ++iterSD, ++i)
    {
      if (i>=size) break;

      val=iterSD->value();
      bptype=iterSD->bptype();
      if (!(bptype == Xyce::Util::BreakPoint::PAUSE && val == finalTime))
      {
        comm->pack( &(val), 1, buf, bsize, pos );
        comm->pack( &(bptype), 1, buf, bsize, pos );
      }
    }
    }

    comm->pack( &savedTimeStep, 1, buf, bsize, pos);

    int sN = analysisManager_.getStepNumber();
    comm->pack( &sN, 1, buf, bsize, pos );
    int tSN = analysisManager_.getTranStepNumber();
    comm->pack( &tSN, 1, buf, bsize, pos );
    int bPRS = analysisManager_.breakPointRestartStep;
    comm->pack( &bPRS, 1, buf, bsize, pos );
    int beginFlag = (analysisManager_.getBeginningIntegrationFlag())?1:0;
    comm->pack( &beginFlag, 1, buf, bsize, pos );
  }
  else
  {
    // count here will be the size for the base StepErrorControl
    // class *only*.
    int count = getRestartDataSize( false );
    int startIndex = pos;

    // Clobber any data in buf lest we leave garbage
    for( int i = startIndex; i < (startIndex+count); ++i) buf[i] = ' ';

    std::ostringstream ost;
    ost.width(24);ost.precision(16);ost.setf(std::ios::scientific);
    ost << startingTimeStep << " ";
    ost << currentTimeStep << " ";
    ost << lastAttemptedTimeStep << " ";
    ost << lastTimeStep << " ";
    ost << minTimeStep << " ";
    ost << maxTimeStep << " ";
    ost << maxTimeStepUser << " ";
    ost << lastTime << " ";
    ost << currentTime << " ";
    ost << nextTime << " ";
    ost << initialTime << " ";
    ost << currentTimeStepRatio << " ";
    ost << currentTimeStepSum << " ";
    ost << lastTimeStepRatio << " ";
    ost << lastTimeStepSum << " ";
    ost << newtonConvergenceStatus << " ";
    ost << numberSuccessiveFailures << " ";
    int flag = (stepAttemptStatus)?1:0;
    ost << flag << " ";
    ost << minStepPrecisionFac_ << " ";
    ost << newtonStepReduction_ << " ";
    ost << tolAimFac_ << " ";
    ost << estOverTol_ << " ";
    // Subtract one because we won't write out the pause breakpoint at the
    // final time
    int size = breakPoints_.size() - 1;
    ost << size << " ";

    BreakPointVector::iterator bpStart = breakPoints_.begin();
    BreakPointVector::iterator bpEnd = breakPoints_.end();
   
    { 
    int i=0;
    for( BreakPointVector::iterator iterSD = bpStart; iterSD != bpEnd; ++iterSD, ++i)
    {
      if (i>=size) break;

      if (!(iterSD->bptype() == Xyce::Util::BreakPoint::PAUSE && iterSD->value() == finalTime))
      {
        ost << iterSD->value() << " ";
        ost << iterSD->bptype() << " ";
      }
    }
    }

    ost << savedTimeStep << " ";
    int sN = analysisManager_.getStepNumber();
    ost << sN << " ";
    int tSN = analysisManager_.getTranStepNumber();
    ost << tSN << " ";
    int bPRS = analysisManager_.breakPointRestartStep;
    ost << bPRS << " ";
    int beginFlag = (analysisManager_.getBeginningIntegrationFlag())?1:0;
    ost << beginFlag << " ";
    std::string data( ost.str() );

    for( unsigned int i = 0; i < data.length(); ++i ) buf[startIndex+i] = data[i];
    // The line above copies the characters of the data string into buf,
    // but doesn't null-terminate buf.
    // it is essential to terminate the buffer with a null, or attempts
    // to construct a string object from it will get memory access problems.
    buf[startIndex+data.length()] = '\0';
    pos += data.length();

    if (DEBUG_RESTART)
    {
      std::string outputString(buf);

      Xyce::dout() << "StepErrorControl UNPACKED output buffer:" << std::endl
                   << outputString << std::endl;
    }
  }

  if (DEBUG_RESTART)
  {
    Xyce::dout() << "TIA Restart Data DUMP (DAE)!  " << netlistFilename_ << "\n"
                 << Xyce::section_divider << std::endl<<std::endl
                 << "alphas_ = " <<alphas_<<std::endl
                 << "alpha0_ = " <<alpha0_<<std::endl
                 << "cj_ = " <<cj_<<std::endl
                 << "ck_ = " <<ck_<<std::endl
                 << "usedStep_ = " <<usedStep_<<std::endl
                 << "Ek_ = " <<Ek_<<std::endl
                 << "Ekm1_ = " <<Ekm1_<<std::endl
                 << "Ekm2_ = " <<Ekm2_<<std::endl
                 << "Ekp1_ = " <<Ekp1_<<std::endl
                 << "Est_ = " <<Est_<<std::endl
                 << "Tk_ = " <<Tk_<<std::endl
                 << "Tkm1_ = " <<Tkm1_<<std::endl
                 << "Tkm2_ = " <<Tkm2_<<std::endl
                 << "Tkp1_ = " <<Tkp1_<<std::endl
                 << "h0_safety_ = " <<h0_safety_<<std::endl
                 << "h0_max_factor_ = " <<h0_max_factor_<<std::endl
                 << "h_phase0_incr_ = " <<h_phase0_incr_<<std::endl
                 << "h_max_inv_ = " <<h_max_inv_<<std::endl
                 << "Tkm1_Tk_safety_ = " <<Tkm1_Tk_safety_<<std::endl
                 << "Tkp1_Tk_safety_ = " <<Tkp1_Tk_safety_<<std::endl
                 << "r_factor_ = " <<r_factor_<<std::endl
                 << "r_safety_ = " <<r_safety_<<std::endl
                 << "r_fudge_ = " <<r_fudge_<<std::endl
                 << "r_min_ = " <<r_min_<<std::endl
                 << "r_max_ = " <<r_max_<<std::endl
                 << "r_hincr_test_ = " <<r_hincr_test_<<std::endl
                 << "r_hincr_ = " <<r_hincr_<<std::endl;

    for (int i=0;i<6;++i)
    {
      Xyce::dout() << "  alpha_["<<i<<"] = " <<  alpha_[i]<<std::endl
                   << "  sigma_["<<i<<"] = " <<  sigma_[i]<<std::endl
                   << "  gamma_["<<i<<"] = " <<  gamma_[i]<<std::endl
                   << "  beta_["<<i<<"] = " <<  beta_[i]<<std::endl
                   << "  psi_["<<i<<"] = " <<  psi_[i]<<std::endl;
    }

    Xyce::dout() << Xyce::section_divider << std::endl << std::endl << std::endl;
  }

  if( pack )
  {
    // doubles:
    comm->pack( &alphas_  , 1, buf, bsize, pos );
    comm->pack( &alpha0_  , 1, buf, bsize, pos );
    comm->pack( &cj_  , 1, buf, bsize, pos );
    comm->pack( &ck_  , 1, buf, bsize, pos );
    comm->pack( &usedStep_  , 1, buf, bsize, pos );
    comm->pack( &Ek_  , 1, buf, bsize, pos );
    comm->pack( &Ekm1_  , 1, buf, bsize, pos );
    comm->pack( &Ekm2_  , 1, buf, bsize, pos );
    comm->pack( &Ekp1_  , 1, buf, bsize, pos );
    comm->pack( &Est_ , 1, buf, bsize, pos );
    comm->pack( &Tk_  , 1, buf, bsize, pos );
    comm->pack( &Tkm1_  , 1, buf, bsize, pos );
    comm->pack( &Tkm2_  , 1, buf, bsize, pos );
    comm->pack( &Tkp1_  , 1, buf, bsize, pos );
    comm->pack( &h0_safety_ , 1, buf, bsize, pos );
    comm->pack( &h0_max_factor_ , 1, buf, bsize, pos );
    comm->pack( &h_phase0_incr_ , 1, buf, bsize, pos );
    comm->pack( &h_max_inv_ , 1, buf, bsize, pos );
    comm->pack( &Tkm1_Tk_safety_  , 1, buf, bsize, pos );
    comm->pack( &Tkp1_Tk_safety_  , 1, buf, bsize, pos );
    comm->pack( &r_factor_  , 1, buf, bsize, pos );
    comm->pack( &r_safety_  , 1, buf, bsize, pos );
    comm->pack( &r_fudge_ , 1, buf, bsize, pos );
    comm->pack( &r_min_ , 1, buf, bsize, pos );
    comm->pack( &r_max_ , 1, buf, bsize, pos );
    comm->pack( &r_hincr_test_  , 1, buf, bsize, pos );
    comm->pack( &r_hincr_ , 1, buf, bsize, pos );

    // vectors of doubles:
    for (int i=0;i<6;++i)
    {
      comm->pack( &alpha_[i], 1, buf, bsize, pos );
      comm->pack( &sigma_[i], 1, buf, bsize, pos );
      comm->pack( &gamma_[i], 1, buf, bsize, pos );
      comm->pack( &beta_[i], 1, buf, bsize, pos );
      comm->pack( &psi_[i], 1, buf, bsize, pos );
    }

    // ints:
    comm->pack ( &currentOrder_, 1, buf, bsize, pos);
    comm->pack ( &maxOrder_, 1, buf, bsize, pos);
    comm->pack ( &minOrder_, 1, buf, bsize, pos);
    comm->pack ( &usedOrder_, 1, buf, bsize, pos);
    comm->pack ( &numberOfSteps_, 1, buf, bsize, pos);
    comm->pack ( &nef_, 1, buf, bsize, pos);
    comm->pack ( &nscsco_, 1, buf, bsize, pos);
    comm->pack ( &newOrder_, 1, buf, bsize, pos);
    comm->pack ( &maxNumfail_, 1, buf, bsize, pos);

    // bools:
    int iP = (initialPhase_)?1:0;
    comm->pack ( &iP, 1, buf, bsize, pos);
  }
  else
  {
    std::ostringstream ost;
    ost.width(24);ost.precision(16);ost.setf(std::ios::scientific);

    // doubles:
    ost << alphas_  << " ";
    ost << alpha0_  << " ";
    ost << cj_  << " ";
    ost << ck_  << " ";
    ost << usedStep_  << " ";
    ost << Ek_  << " ";
    ost << Ekm1_  << " ";
    ost << Ekm2_  << " ";
    ost << Ekp1_  << " ";
    ost << Est_ << " ";
    ost << Tk_  << " ";
    ost << Tkm1_  << " ";
    ost << Tkm2_  << " ";
    ost << Tkp1_  << " ";
    ost << h0_safety_ << " ";
    ost << h0_max_factor_ << " ";
    ost << h_phase0_incr_ << " ";
    ost << h_max_inv_ << " ";
    ost << Tkm1_Tk_safety_  << " ";
    ost << Tkp1_Tk_safety_  << " ";
    ost << r_factor_  << " ";
    ost << r_safety_  << " ";
    ost << r_fudge_ << " ";
    ost << r_min_ << " ";
    ost << r_max_ << " ";
    ost << r_hincr_test_  << " ";
    ost << r_hincr_ << " ";

    // vectors of doubles:
    for (int i=0;i<6;++i)
    {
      ost << alpha_[i]<< " " ;
      ost << sigma_[i]<< " " ;
      ost << gamma_[i]<< " " ;
      ost << beta_[i]<< " " ;
      ost << psi_[i]<< " " ;
    }

    // ints:
    ost << currentOrder_ << " ";
    ost << maxOrder_ << " ";
    ost << minOrder_ << " ";
    ost << usedOrder_ << " ";
    ost << numberOfSteps_ << " ";
    ost << nef_ << " ";
    ost << nscsco_ << " ";
    ost << newOrder_ << " ";
    ost << maxNumfail_ << " ";

    // bools:
    int iP = (initialPhase_)?1:0;
    ost << iP << " ";

    std::string data( ost.str() );
    for( unsigned int i = 0; i < data.length(); ++i ) buf[pos+i] = data[i];
    // null terminate buf, essential if buf is used as a string anywhere,
    // including in the string constructor below:
    buf[pos+data.length()] = '\0';

    if (DEBUG_RESTART) {
      std::string outputString(buf);

      Xyce::dout() << "StepErrorControlDAE  UNPACKED output buffer:" << std::endl
                   << outputString << std::endl;
    }

    pos = newPos;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::restoreRestartData
// Purpose       :
// Special Notes : ERK:  6/20/2010:
//                 There are no longer 2 distinct classes.  The DAE version is
//                 now part of the original, as one single class.  The two
//                 blocks of code in this function correspond to the original
//                 class (first block), and the DAE class (second block).
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 7/28/06
//-----------------------------------------------------------------------------
bool StepErrorControl::restoreRestartData(
  char *                    buf,
  int                       bsize,
  int &                     pos,
  Parallel::Communicator *  comm,
  bool                      pack,
  double &                  initial_time)
{
  // original class variables:
  if( pack )
  {
    comm->unpack(buf, bsize, pos, &startingTimeStep, 1);
    comm->unpack(buf, bsize, pos, &currentTimeStep, 1);
    comm->unpack(buf, bsize, pos, &lastAttemptedTimeStep, 1);
    comm->unpack(buf, bsize, pos, &lastTimeStep, 1);
    comm->unpack(buf, bsize, pos, &minTimeStep, 1);
    comm->unpack(buf, bsize, pos, &maxTimeStep, 1);
    comm->unpack(buf, bsize, pos, &maxTimeStepUser, 1);
    comm->unpack(buf, bsize, pos, &lastTime, 1);
    comm->unpack(buf, bsize, pos, &currentTime, 1);
    comm->unpack(buf, bsize, pos, &nextTime, 1);
    comm->unpack(buf, bsize, pos, &initialTime, 1);
    comm->unpack(buf, bsize, pos, &currentTimeStepRatio, 1);
    comm->unpack(buf, bsize, pos, &currentTimeStepSum, 1);
    comm->unpack(buf, bsize, pos, &lastTimeStepRatio, 1);
    comm->unpack(buf, bsize, pos, &lastTimeStepSum, 1);
    comm->unpack(buf, bsize, pos, &newtonConvergenceStatus, 1);
    comm->unpack(buf, bsize, pos, &numberSuccessiveFailures, 1);
    int flag;
    comm->unpack(buf, bsize, pos, &flag, 1);
    stepAttemptStatus = flag;
    comm->unpack(buf, bsize, pos, &minStepPrecisionFac_, 1);
    comm->unpack(buf, bsize, pos, &newtonStepReduction_, 1);
    comm->unpack(buf, bsize, pos, &tolAimFac_, 1);
    comm->unpack(buf, bsize, pos, &estOverTol_, 1);

    double bpTol = 2.0 * minTimeStep;
    breakPointLess_.tolerance_ = bpTol;
    initial_time = initialTime;
    pauseTime = initialTime;

    std::vector<Util::BreakPoint>::iterator lastToDelete = std::lower_bound(breakPoints_.begin(), breakPoints_.end(), currentTime, breakPointLess_);    
    breakPoints_.erase( breakPoints_.begin(), lastToDelete);

    int size;
    double val;
    int bptype;
    comm->unpack(buf, bsize, pos, &size, 1);

    for (int i = 0; i < size; ++i)
    {
      comm->unpack(buf, bsize, pos, &val, 1);
      comm->unpack(buf, bsize, pos, &bptype, 1);
      if (val > currentTime)
      {
        Util::BreakPoint TmpBP(val, bptype);
        breakPoints_.push_back(TmpBP);
      }
    }

    std::sort ( breakPoints_.begin(), breakPoints_.end(), breakPointLess_ );
    std::vector<Util::BreakPoint>::iterator it = std::unique ( breakPoints_.begin(), breakPoints_.end(), breakPointEqual_ );
    breakPoints_.resize( std::distance (breakPoints_.begin(), it )); 
    for (int i=0;i<breakPoints_.size();++i)
    {
      Util::BreakPoint & TmpBP = breakPoints_[i];
      if (TmpBP.bptype() == Xyce::Util::BreakPoint::PAUSE)
      {
        updatePauseTime(TmpBP, initialTime);
      }
    }
   
    comm->unpack(buf, bsize, pos, &savedTimeStep, 1);

    int im;
    comm->unpack(buf, bsize, pos, &im, 1);
    analysisManager_.setStepNumber(im);
    comm->unpack(buf, bsize, pos, &im, 1);
    analysisManager_.setTranStepNumber(im);
    comm->unpack(buf, bsize, pos, &im, 1);
    analysisManager_.breakPointRestartStep = im;
    comm->unpack(buf, bsize, pos, &im, 1);
    analysisManager_.setBeginningIntegrationFlag(im==1);
  }
  else
  {
    std::string str1(buf);
    int length = str1.size() - pos;
    std::string str2(str1,pos,length);

    std::istringstream ist( str2 );

    // count here will be the size for the base StepErrorControl
    // class *only*.

    // THIS IS INCORRECT!  restartDataSize returns the MAXIMUM the thing
    // is allowed to be, which allows for 24 characters for each integer
    // value in the output.  This is almost always an overestimate, and
    // using it as a way of pointing to the next character after our
    // data is a sure-fire way to break downstream processing!
    // Instead, we'll use the tellg function from the stringstream to
    // return the offset after we read everything out.
    //    int count = StepErrorControl::restartDataSize( false );
    //    pos += count;

    ist >> startingTimeStep;
    ist >> currentTimeStep;
    ist >> lastAttemptedTimeStep;
    ist >> lastTimeStep;
    ist >> minTimeStep;
    ist >> maxTimeStep;
    ist >> maxTimeStepUser;
    ist >> lastTime;
    ist >> currentTime;
    ist >> nextTime;
    ist >> initialTime;
    ist >> currentTimeStepRatio;
    ist >> currentTimeStepSum;
    ist >> lastTimeStepRatio;
    ist >> lastTimeStepSum;
    ist >> newtonConvergenceStatus;
    ist >> numberSuccessiveFailures;
    int flag;
    ist >> flag;
    stepAttemptStatus = flag;
    ist >> minStepPrecisionFac_;
    ist >> newtonStepReduction_;
    ist >> tolAimFac_;
    ist >> estOverTol_;

    double bpTol = 2.0 * minTimeStep;
    breakPointLess_.tolerance_ = bpTol;
    initial_time = initialTime;
    pauseTime=initialTime;

    std::vector<Util::BreakPoint>::iterator lastToDelete = std::lower_bound(breakPoints_.begin(), breakPoints_.end(), currentTime, breakPointLess_);    
    breakPoints_.erase( breakPoints_.begin(), lastToDelete);

    int size;
    double val;
    int bptype;
    ist >> size;

    for( int i = 0; i < size; ++i )
    {
      ist >> val;
      ist >> bptype;
      if (val > currentTime)
      {
        Util::BreakPoint TmpBP(val, bptype);
        breakPoints_.push_back(TmpBP);
      }
    }

    std::sort ( breakPoints_.begin(), breakPoints_.end(), breakPointLess_ );
    std::vector<Util::BreakPoint>::iterator it = std::unique ( breakPoints_.begin(), breakPoints_.end(), breakPointEqual_ );
    breakPoints_.resize( std::distance (breakPoints_.begin(), it )); 
    for (int i=0;i<breakPoints_.size();++i)
    {
      Util::BreakPoint & TmpBP = breakPoints_[i];
      if (TmpBP.bptype() == Xyce::Util::BreakPoint::PAUSE)
      {
        updatePauseTime(TmpBP, initialTime);
      }
    }

    ist >> savedTimeStep;

    int tmpInt;
    ist >> tmpInt;
    analysisManager_.setStepNumber(tmpInt);
    ist >> tmpInt;
    analysisManager_.setTranStepNumber(tmpInt);
    ist >> tmpInt;
    analysisManager_.breakPointRestartStep = tmpInt;
    ist >> tmpInt;
    analysisManager_.setBeginningIntegrationFlag(tmpInt==1);

    pos += ist.tellg();
  }

  if (DEBUG_RESTART)
  {
    Xyce::dout() << "TIA Restart Data RESTORE!  " << netlistFilename_ << "\n"
                 << Xyce::section_divider << std::endl
                 << "startingTimeStep: " << startingTimeStep <<  std::endl
                 << "currentTimeStep: " << currentTimeStep << std::endl
                 << "lastAttemptedTimeStep: " << lastAttemptedTimeStep << std::endl
                 << "lastTimeStep: " << lastTimeStep << std::endl
                 << "minTimeStep: " << minTimeStep << std::endl
                 << "maxTimeStep: " << maxTimeStep << std::endl
                 << "maxTimeStepUser: " << maxTimeStepUser << std::endl
                 << "lastTime: " << lastTime << std::endl
                 << "currentTime: " << currentTime << std::endl
                 << "nextTime: " << nextTime << std::endl
                 << "initialTime: " << initialTime << std::endl
                 << "estOverTol_: " << estOverTol_ << std::endl
                 << "breakpts: ";
    for (BreakPointVector::iterator iterSD = breakPoints_.begin();
         iterSD != breakPoints_.end(); ++iterSD)
    {
      Xyce::dout() << iterSD->value() << " " << iterSD->bptype() << std::endl;
    }
    Xyce::dout() << std::endl
                 << "integMethod: " << analysisManager_.getIntegrationMethod() << std::endl
                 << "stepNumber: " << analysisManager_.getStepNumber() << std::endl
                 << "tranStepNumber: " << analysisManager_.getTranStepNumber() << std::endl
                 << "breakPointRestartStep: " << analysisManager_.breakPointRestartStep << std::endl
                 << Xyce::section_divider << std::endl << std::endl;
  }

  // DAE class variables:
  if( pack )
  {
    // doubles:
    comm->unpack(buf, bsize, pos, &alphas_  , 1);
    comm->unpack(buf, bsize, pos, &alpha0_  , 1);
    comm->unpack(buf, bsize, pos, &cj_  , 1);
    comm->unpack(buf, bsize, pos, &ck_  , 1);
    comm->unpack(buf, bsize, pos, &usedStep_  , 1);
    comm->unpack(buf, bsize, pos, &Ek_  , 1);
    comm->unpack(buf, bsize, pos, &Ekm1_  , 1);
    comm->unpack(buf, bsize, pos, &Ekm2_  , 1);
    comm->unpack(buf, bsize, pos, &Ekp1_  , 1);
    comm->unpack(buf, bsize, pos, &Est_ , 1);
    comm->unpack(buf, bsize, pos, &Tk_  , 1);
    comm->unpack(buf, bsize, pos, &Tkm1_  , 1);
    comm->unpack(buf, bsize, pos, &Tkm2_  , 1);
    comm->unpack(buf, bsize, pos, &Tkp1_  , 1);
    comm->unpack(buf, bsize, pos, &h0_safety_ , 1);
    comm->unpack(buf, bsize, pos, &h0_max_factor_ , 1);
    comm->unpack(buf, bsize, pos, &h_phase0_incr_ , 1);
    comm->unpack(buf, bsize, pos, &h_max_inv_ , 1);
    comm->unpack(buf, bsize, pos, &Tkm1_Tk_safety_  , 1);
    comm->unpack(buf, bsize, pos, &Tkp1_Tk_safety_  , 1);
    comm->unpack(buf, bsize, pos, &r_factor_  , 1);
    comm->unpack(buf, bsize, pos, &r_safety_  , 1);
    comm->unpack(buf, bsize, pos, &r_fudge_ , 1);
    comm->unpack(buf, bsize, pos, &r_min_ , 1);
    comm->unpack(buf, bsize, pos, &r_max_ , 1);
    comm->unpack(buf, bsize, pos, &r_hincr_test_  , 1);
    comm->unpack(buf, bsize, pos, &r_hincr_ , 1);

    // vectors of doubles:
    for (int i=0;i<6;++i)
    {
      comm->unpack(buf, bsize, pos, &alpha_[i], 1);
      comm->unpack(buf, bsize, pos, &sigma_[i], 1);
      comm->unpack(buf, bsize, pos, &gamma_[i], 1);
      comm->unpack(buf, bsize, pos, &beta_[i], 1);
      comm->unpack(buf, bsize, pos, &psi_[i], 1);
    }

    // ints:
    comm->unpack(buf, bsize, pos, &currentOrder_, 1);
    comm->unpack(buf, bsize, pos, &maxOrder_, 1);
    comm->unpack(buf, bsize, pos, &minOrder_, 1);
    comm->unpack(buf, bsize, pos, &usedOrder_, 1);
    comm->unpack(buf, bsize, pos, &numberOfSteps_, 1);
    comm->unpack(buf, bsize, pos, &nef_, 1);
    comm->unpack(buf, bsize, pos, &nscsco_, 1);
    comm->unpack(buf, bsize, pos, &newOrder_, 1);
    comm->unpack(buf, bsize, pos, &maxNumfail_, 1);

    // bools:
    int iP;
    comm->unpack (buf, bsize, pos, &iP, 1);
    if (iP == 0) initialPhase_ = false;
    else         initialPhase_ = true;
  }
  else
  {
    // want the string stream to only represent the new-DAE part of
    // the buf array.
    std::string str1(buf);
    int length = str1.size() - pos;
    std::string str2(str1,pos,length);

    std::istringstream ist( str2 );

    int count = getRestartDataSize(false);
    pos = count;  // can update here, as pos isn't used after this.

    // doubles:
    ist >> alphas_;
    ist >> alpha0_;
    ist >> cj_;
    ist >> ck_;
    ist >> usedStep_;
    ist >> Ek_;
    ist >> Ekm1_;
    ist >> Ekm2_;
    ist >> Ekp1_;
    ist >> Est_;
    ist >> Tk_;
    ist >> Tkm1_;
    ist >> Tkm2_;
    ist >> Tkp1_;
    ist >> h0_safety_;
    ist >> h0_max_factor_;
    ist >> h_phase0_incr_;
    ist >> h_max_inv_;
    ist >> Tkm1_Tk_safety_;
    ist >> Tkp1_Tk_safety_;
    ist >> r_factor_;
    ist >> r_safety_;
    ist >> r_fudge_;
    ist >> r_min_;
    ist >> r_max_;
    ist >> r_hincr_test_;
    ist >> r_hincr_;

    // vectors of doubles:
    for (int i=0;i<6;++i)
    {
      ist >> alpha_[i];
      ist >> sigma_[i];
      ist >> gamma_[i];
      ist >> beta_[i];
      ist >> psi_[i];
    }

    // ints:
    ist >> currentOrder_;
    ist >> maxOrder_;
    ist >> minOrder_;
    ist >> usedOrder_;
    ist >> numberOfSteps_;
    ist >> nef_;
    ist >> nscsco_;
    ist >> newOrder_;
    ist >> maxNumfail_;

    // bools:
    int iP;
    ist >> iP;
    if (iP == 0) initialPhase_ = false;
    else         initialPhase_ = true;
  }

  if (DEBUG_RESTART)
  {
    
    Xyce::dout() << "TIA Restart Data RESTORE (DAE)!  " << netlistFilename_ << "\n"
                 << Xyce::section_divider << std::endl<<std::endl
                 << "alphas_ = " <<alphas_<<std::endl
                 << "alpha0_ = " <<alpha0_<<std::endl
                 << "cj_ = " <<cj_<<std::endl
                 << "ck_ = " <<ck_<<std::endl
                 << "usedStep_ = " <<usedStep_<<std::endl
                 << "Ek_ = " <<Ek_<<std::endl
                 << "Ekm1_ = " <<Ekm1_<<std::endl
                 << "Ekm2_ = " <<Ekm2_<<std::endl
                 << "Ekp1_ = " <<Ekp1_<<std::endl
                 << "Est_ = " <<Est_<<std::endl
                 << "Tk_ = " <<Tk_<<std::endl
                 << "Tkm1_ = " <<Tkm1_<<std::endl
                 << "Tkm2_ = " <<Tkm2_<<std::endl
                 << "Tkp1_ = " <<Tkp1_<<std::endl
                 << "h0_safety_ = " <<h0_safety_<<std::endl
                 << "h0_max_factor_ = " <<h0_max_factor_<<std::endl
                 << "h_phase0_incr_ = " <<h_phase0_incr_<<std::endl
                 << "h_max_inv_ = " <<h_max_inv_<<std::endl
                 << "Tkm1_Tk_safety_ = " <<Tkm1_Tk_safety_<<std::endl
                 << "Tkp1_Tk_safety_ = " <<Tkp1_Tk_safety_<<std::endl
                 << "r_factor_ = " <<r_factor_<<std::endl
                 << "r_safety_ = " <<r_safety_<<std::endl
                 << "r_fudge_ = " <<r_fudge_<<std::endl
                 << "r_min_ = " <<r_min_<<std::endl
                 << "r_max_ = " <<r_max_<<std::endl
                 << "r_hincr_test_ = " <<r_hincr_test_<<std::endl
                 << "r_hincr_ = " <<r_hincr_<<std::endl;

    for (int i=0;i<6;++i)
    {
      Xyce::dout() << "  alpha_["<<i<<"] = " <<  alpha_[i]<<std::endl
                   << "  sigma_["<<i<<"] = " <<  sigma_[i]<<std::endl
                   << "  gamma_["<<i<<"] = " <<  gamma_[i]<<std::endl
                   << "  beta_["<<i<<"] = " <<  beta_[i]<<std::endl
                   << "  psi_["<<i<<"] = " <<  psi_[i]<<std::endl;
    }

    Xyce::dout() << Xyce::section_divider << std::endl << std::endl << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::updateTwoLevelTimeInfo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel ComputationalSciences.
// Creation Date : 3/14/06
//-----------------------------------------------------------------------------
void StepErrorControl::updateTwoLevelTimeInfo(
  Parallel::Machine     comm,
  double                tiInfonextTimeStep,
  double                tiInfonextTime,
  int                   tiInfocurrentOrder,
  bool                  breakpoints_enabled,
  double                initial_time,
  bool                  min_time_steps_breakpoint_given,
  double                min_time_steps_breakpoint)
{
  updateStopTime(comm, breakpoints_enabled, initial_time, min_time_steps_breakpoint_given, min_time_steps_breakpoint);

  // Need to get this right.
  if (previousCallStepSuccessful == true)
  {
    lastTimeStep      = currentTimeStep;
    lastTimeStepRatio = currentTimeStepRatio;
    lastTimeStepSum   = currentTimeStepSum;
    previousCallStepSuccessful = false;
  }

  // If the step needs to be adjusted:
  lastAttemptedTimeStep = currentTimeStep;
  double newTimeStep = tiInfonextTimeStep;
  nextTime = tiInfonextTime;
  currentTimeStepRatio = newTimeStep/lastTimeStep;
  currentTimeStepSum   = newTimeStep + lastTimeStep;
  currentTimeStep = newTimeStep;
  currentOrder_ = tiInfocurrentOrder;
}

//-----------------------------------------------------------------------------
// Function      : StepErrorControl::outputTimeInfo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 03/04/08
//-----------------------------------------------------------------------------
void StepErrorControl::outputTimeInfo(std::ostream &os)
{
  {
    Xyce::ios_flags_saver flagsave(os);
    
    os << " " 
      //<< (DEBUG_TIME ? netlistFilename_ : "   ") << " "
      << netlistFilename_  << " "
       << "Current,Next,Step: "
       //<< std::setw(DEBUG_TIME ? 14 : 16) << std::setprecision(Xyce::DEBUG_TIME ? 7 : 9)
       << std::scientific
       << std::setw(14) << std::setprecision(7)
       << currentTime << ", " << nextTime << ", " << currentTimeStep;
  }

  os << "    ("<<numberOfSteps_<<") Used, Next Order:  " << usedOrder_ << ", " << currentOrder_ << std::endl;
}


//-----------------------------------------------------------------------------
// Function      : StepErrorControl::isPauseTime
// Purpose       : public method to flag whether the current time is a pause
//                 time
// Special Notes : Need this method because the prior values get in the way 
//                 when we resume.
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/30/04
//-----------------------------------------------------------------------------
bool StepErrorControl::isPauseTime()
{
  // Check whether pauseTime == currentTime (to within bpTol).
  // If the pause time is equal to the final time, we are at the end, not
  // at a pause time.  
  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << "StepErrorControl::isPauseTime\n";
    if ( !(breakPoints_.empty()) && currentPauseBP != breakPoints_.end() ) // if not valid don't do this
    {
      Xyce::dout() << "currentPauseBP.value = " << currentPauseBP->value () << std::endl;
    }
    else
    {
      Xyce::dout() << "currentPauseBP not valid" <<std::endl;
    }
    Xyce::dout() << "final Time           = " << finalTime  << std::endl;
    Xyce::dout() << "current Time         = " << currentTime  << std::endl;
  }

  if ( !(breakPoints_.empty()) && currentPauseBP != breakPoints_.end() ) // if not valid don't do this
  {
    return (breakPointLess_(*currentPauseBP, finalTime) || breakPointLess_(finalTime, *currentPauseBP))
      && (!breakPointLess_(*currentPauseBP, currentTime) && !breakPointLess_(currentTime, *currentPauseBP));
  }
  else
  {
    return false;
  }
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : "<<" operator for step error control class.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/17/05
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const StepErrorControl & sec)
{
  os << "\n\n-----------------------------------------" << std::endl
     << "\tStepErrorControl:\n"
     << "\t\tstartingTimeStep      = " << sec.startingTimeStep << std::endl 
     << "\t\tcurrentTimeStep       = " << sec.currentTimeStep << std::endl 
     << "\t\tlastAttemptedTimeStep = " << sec.lastAttemptedTimeStep << std::endl 
     << "\t\tlastTimeStep          = " << sec.lastTimeStep << std::endl 
     << "\t\tminTimeStep           = " << sec.minTimeStep << std::endl 
     << "\t\tmaxTimeStep           = " << sec.maxTimeStep << std::endl 
     << "\t\tmaxTimeStepUser       = " << sec.maxTimeStepUser << std::endl 
     << "\t\tlastTime              = " << sec.lastTime << std::endl 
     << "\t\tcurrentTime           = " << sec.currentTime << std::endl 
     << "\t\tnextTime              = " << sec.nextTime << std::endl 
     << "\t\tstopTime              = " << sec.stopTime << std::endl 
     << "\t\tinitialTime           = " << sec.initialTime << std::endl 
     << "\t\tfinalTime             = " << sec.finalTime << std::endl 
     << std::endl
     << "\t\tBreak Points:" << std::endl;

  sec.printBreakPoints (os);


  os << Xyce::section_divider << std::endl
     << std::endl;

  return os;
}

} // namespace TimeIntg
} // namespace Xyce

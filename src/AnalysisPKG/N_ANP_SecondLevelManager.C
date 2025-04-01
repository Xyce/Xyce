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
// Purpose       : This file contains the functions which define the time
//                 domain & integration algorithm classes.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_SecondLevelManager.h>
#include <N_ANP_Transient.h>
#include <N_ERH_Message.h>
#include <N_IO_ActiveOutput.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_LOA_CktLoader.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_NLS_Manager.h>
#include <N_PDS_Comm.h>
#include <N_PDS_Manager.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_WorkingIntegrationMethod.h>
#include <N_TIA_fwd.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Diagnostic.h>
#include <N_ANP_OutputConductanceFile.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setExternalSolverState 
// Purpose       : 
// Special Notes : Used for multi-level Newton solves, for levels other
//                 than the top level.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/22/2014
//-----------------------------------------------------------------------------
void SecondLevelManager::setExternalSolverState(
  Loader::CktLoader &           loader,
  bool                          external_initJctFlag)
{
  loader.setExternalSolverState(external_initJctFlag);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::runSecondLevelStep
//
// Purpose       : This function is similar to "run" except that only a single
//                 integration (or DC sweep) step will be executed.
//
// Special Notes : Used for multi-level Newton solves, for levels other
//                 than the top level.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/06/2006
//-----------------------------------------------------------------------------
bool SecondLevelManager::runSecondLevelStep(TimeIntg::TwoLevelError & tlError)
{
  bool status = twoLevelAnalysisObject_->twoLevelStep();

  getWorkingIntegrationMethod().getTwoLevelError(tlError);

  return status;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::startupSolvers
// Purpose       :
// Special Notes : Used only for 2-level solves.
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 3/10/2006
//-----------------------------------------------------------------------------
bool SecondLevelManager::startupSecondLevelSolvers(
  Linear::System &              linear_system,
  Nonlinear::Manager &          nonlinear_manager)
{
  bool bsuccess = true;

  // Hardwire the erroption parameter to 1, which will force the inner
  // solve (initiated by this function) to only use the Newton success/failure
  // as step criteria.  Predictor-corrector information is handled in the
  // upper level solver.
  getTIAParams().errorAnalysisOption = TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES;


  // If running Xyce-to-Xyce coupling, there should already be a primary 
  // analysis object allocated at this stage.  So grab that one.  
  twoLevelAnalysisObject_ = 0;
  twoLevelAnalysisObject_ = getAnalysisObjectPtr();

  if (twoLevelAnalysisObject_==0)
  {
    Report::UserError() << "Primary Analysis Object not allocated";
    return false;
  }

  twoLevelAnalysisObject_->init();

  activeOutput_ = new IO::ActiveOutput(getOutputManagerAdapter().getOutputManager());
  activeOutput_->add(getPDSManager()->getPDSComm()->comm(), getAnalysisMode());

  // Reset the solvers timer.
//  xyceTranTimerPtr_.resetStartTime();

 return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::finishSolvers
// Purpose       :
// Special Notes : Used only for 2-level solves.
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 3/10/2006
//-----------------------------------------------------------------------------
bool SecondLevelManager::finishSecondLevelSolvers()
{
  bool bsuccess = true;

  twoLevelAnalysisObject_->finish();

  delete activeOutput_;
  activeOutput_ = 0;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::homotopyStepSuccess
// Purpose       : Lower-level processing of a successful homotopy step,
//                 which was controlled from the upper level of a 2-level solve.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void SecondLevelManager::homotopyStepSuccess(
  const std::vector<std::string> &      paramNames,
  const std::vector<double> &           paramVals)
{
  if (DEBUG_ANALYSIS)
    Xyce::dout() << "\n " << getNetlistFilename()
                 << " AnalysisManager::homotopyStepSuccess " << std::endl;

  // output:
  getOutputManagerAdapter().outputHomotopy( paramNames, paramVals, *getDataStore()->nextSolutionPtr );

  // update the data arrays:
  getDataStore()->updateSolDataArrays();

  // pass info to the next level down, if it exists.
  getNonlinearEquationLoader().homotopyStepSuccess(paramNames,paramVals);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::homotopyStepFailure
//
// Purpose       : Lower-level processing of a failed homotopy step,
//                 which was controlled from the upper level of a
//                 2-level solve.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/30/06
//-----------------------------------------------------------------------------
void SecondLevelManager::homotopyStepFailure()
{
  if (DEBUG_ANALYSIS)
    Xyce::dout() << "\n " << getNetlistFilename()
                 << " AnalysisManager::homotopyStepFailure " << std::endl;

  // The solutions currently in place represent failure.  Get rid of them.
  getDataStore()->usePreviousSolAsPredictor ();

  // pass info to the next level down, if it exists.
  getNonlinearEquationLoader().homotopyStepFailure ();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::stepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void SecondLevelManager::stepSecondLevelSuccess(TwoLevelMode analysisUpper)
{
  if (DEBUG_ANALYSIS)
    Xyce::dout() << "\n " << getNetlistFilename()
                 << " AnalysisManager::stepSuccess " << std::endl;

  setTwoLevelMode(analysisUpper);
  getStepErrorControl().stepAttemptStatus = true;
  switch (analysisUpper)
  {
    case Analysis::TWO_LEVEL_MODE_TRANSIENT_DCOP:
      {
        Transient * twoLevelTransientAnalysisObject = dynamic_cast<Transient *>(twoLevelAnalysisObject_);
        if (twoLevelTransientAnalysisObject)
        {
          twoLevelTransientAnalysisObject->processSuccessfulDCOP();
        }
        else
        {
          Report::DevelFatal().in("AnalysisManager::stepSuccess") 
            << "Failed dynamic_cast of twoLevelAnalysisObject to Transient.";
        }
      }
      break;

    case Analysis::TWO_LEVEL_MODE_TRANSIENT:
      twoLevelAnalysisObject_->processSuccessfulStep();
      break;

    case Analysis::TWO_LEVEL_MODE_DC_SWEEP:
      twoLevelAnalysisObject_->processSuccessfulStep();
      break;

    default:
      Report::DevelFatal().in("AnalysisManager::stepSecondLevelSuccess") 
        << "TwoLevelMode " << analysisUpper << " is not available";
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::stepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void SecondLevelManager::stepSecondLevelFailure(TwoLevelMode analysisUpper)
{
  setTwoLevelMode(analysisUpper);
  getStepErrorControl().stepAttemptStatus = false;
  switch (analysisUpper)
  {
    case Analysis::TWO_LEVEL_MODE_TRANSIENT_DCOP:
      {
        Transient * twoLevelTransientAnalysisObject = dynamic_cast<Transient *>(twoLevelAnalysisObject_);
        if (twoLevelTransientAnalysisObject)
        {
          twoLevelTransientAnalysisObject->processFailedDCOP();
        }
        else
        {
          Report::DevelFatal().in("AnalysisManager::stepFailure") << "Failed dynamic_cast of twoLevelAnalysisObject to Transient.";
        }
      }
      break;
    case Analysis::TWO_LEVEL_MODE_TRANSIENT:
      twoLevelAnalysisObject_->processFailedStep();
      break;
    case Analysis::TWO_LEVEL_MODE_DC_SWEEP:
      twoLevelAnalysisObject_->processFailedStep();
      break;
    default:
      Report::DevelFatal().in("AnalysisManager::stepSecondLevelFailure") << "TwoLevelMode " << analysisUpper << " is not available";
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getInitialQnorm
// Purpose       : Used for 2-level solves.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
bool SecondLevelManager::getSecondLevelInitialQnorm
  (TimeIntg::TwoLevelError & tle) const
{
  getWorkingIntegrationMethod().getInitialQnorm(tle);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getBreakPoints
// Purpose       : Used for 2-level solves.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
bool SecondLevelManager::getSecondLevelBreakPoints(
  Loader::CktLoader &                   loader,
  std::vector<Util::BreakPoint> &       breakPointTimes,
  std::vector<Util::BreakPoint> &       pauseBreakPointTimes
  ) //const
{
  if (!breakPointsRequestedBefore_)
  {
    loader.resetBreakPoints();
  }

  loader.getBreakPoints(breakPointTimes, pauseBreakPointTimes);
  breakPointsRequestedBefore_ = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::startSecondLevelTimeStep
// Purpose       : used by 2-level solves.
// Special Notes : One of the primary purposes for this function is to impose
//                 a lot of upper level information from the top level time
//                 integrator on the inner level time integrator.
//
//                 In general, in a 2-level solve, an inner solver doesn't
//                 have enough information to correctly determine the
//                 step size, order, etc.  This is in part because the
//                 inner solver, while it knows about the top level solver,
//                 it cannot know about any OTHER inner solvers.
//
//                 The top level solver, however, does have enough information.
//                 It gathers break point, error analysis, and other info
//                 from all the inner solves.  So, the top level solver
//                 makes all the decisions and imposes them on the inner
//                 solves.  This function is where it does that impose.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/06/2006
//-----------------------------------------------------------------------------
bool SecondLevelManager::startSecondLevelTimeStep(
  const TimeIntg::TIAParams &   tia_params,
  Nonlinear::Manager &          nonlinear_manager,
  bool                          beginIntegrationFlag,
  double                        nextTimeStep,
  double                        nextTime,
  int                           currentOrder)
{
  twoLevelAnalysisObject_->setBeginningIntegrationFlag(beginIntegrationFlag);

  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    Xyce::dout() << "AnalysisManager::startTimeStep:" << std::endl;
  }

  if (VERBOSE_TIME && twoLevelAnalysisObject_)
  {
    twoLevelAnalysisObject_->printStepHeader(Xyce::lout());
  }

  if (getSwitchIntegrator())
    createTimeIntegratorMethod(tia_params, twoLevelAnalysisObject_->getIntegrationMethod());

  // ------------------------------------------------------------------------
  // Set the step size, current time and next time.
#ifndef Xyce_CHARON
  if ( twoLevelAnalysisObject_->getIntegrationMethod() != TimeIntg::methodsEnum::NO_TIME_INTEGRATION)
#endif
  {
    getStepErrorControl().updateTwoLevelTimeInfo(
      getPDSManager()->getPDSComm()->comm(),
      nextTimeStep,
      nextTime,
      currentOrder,
      tia_params.bpEnable,
      tia_params.initialTime,
      tia_params.minTimeStepsBPGiven,
      tia_params.minTimeStepsBP);
  }

  if (twoLevelAnalysisObject_->getBeginningIntegrationFlag() && getStepErrorControl().stepAttemptStatus)
  {
    getWorkingIntegrationMethod().setTwoLevelTimeInfo(); // new-dae only
  }

  // ------------------------------------------------------------------------
  // If we've switched the integration method, we need to obtain the
  // corrector derivative only after we've updated the TimeInfo.
  if (getSwitchIntegrator())
  {
    setSwitchIntegrator(false);
    getWorkingIntegrationMethod().obtainCorrectorDeriv();
  }

  bool dcopFlag = true;
  Transient *twoLevelTransientAnalysisObject = dynamic_cast<Transient *>(twoLevelAnalysisObject_);
  if (twoLevelTransientAnalysisObject) {
    dcopFlag = twoLevelTransientAnalysisObject->getDCOPFlag();
  }

  if (VERBOSE_TIME && !dcopFlag)
    getStepErrorControl().outputTimeInfo(lout());

  // ------------------------------------------------------------------------
  // Set the nonlinear solver parameters to those appropriate for the
  // transient solution, if neccessary.
  if (!dcopFlag)
  {
    nonlinear_manager.setAnalysisMode(nonlinearAnalysisMode(ANP_MODE_TRANSIENT));
  }

  // Ask the method to update its coefficients
  getWorkingIntegrationMethod().updateCoeffs(); 
  twoLevelAnalysisObject_->handlePredictor();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : SecondLevelManager::setTwoLevelParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/15/2022
//-----------------------------------------------------------------------------
bool SecondLevelManager::setTwoLevelParams(
  const Util::OptionBlock &     option_block)
{
  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; ++it)
  {
    const Util::Param &param = (*it);

    if (param.uTag() == "OUTPUT_DAE_VECTORS")
    {
      outputDAEvectors_ = param.getImmutableValue<bool>();
    }
    else if (param.uTag() == "OUTPUT_DAE_VECTORS_NOPORT")
    {
      outputDAEvectors_noport_ = param.getImmutableValue<bool>();
    }
    else if (param.uTag() == "OUTPUT_DAE_MATRICES")
    {
      outputDAEmatrices_ = param.getImmutableValue<bool>();
    }
    else if (param.uTag() == "OUTPUT_REDUCED_CONDUCTANCES")
    {
      condOutputFlag_ = static_cast<bool> (param.getImmutableValue<int>());
    }
    else if (param.uTag() == "OUTPUT_PORT_CURRENTS")
    {
      portCurrentOutputFlag_ = static_cast<bool> (param.getImmutableValue<int>());
    }
    else
    {
      Report::UserError() << param.uTag() << " is not a recognized two-level analysis option";
      return false;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : registerTwoLevelPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/2022
//-----------------------------------------------------------------------------
bool registerTwoLevelPkgOptionsMgr(
    SecondLevelManager &second_level_manager,
    IO::PkgOptionsMgr &options_manager)
{
  Xyce::Analysis::SecondLevelManager::populateMetadata(options_manager);

  options_manager.addCommandProcessor("TWOLEVEL", 
    IO::createRegistrationOptions(second_level_manager, &SecondLevelManager::setTwoLevelParams));

  return true;
}

//-----------------------------------------------------------------------------
// Function      : SecondLevelManager::populateMetadata
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/2022
//-----------------------------------------------------------------------------
void SecondLevelManager::populateMetadata(IO::PkgOptionsMgr &options_manager)
{
  Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("TWOLEVEL");
  parameters.insert(Util::ParamMap::value_type("OUTPUT_DAE_VECTORS", Util::Param("OUTPUT_DAE_VECTORS", 0)));
  parameters.insert(Util::ParamMap::value_type("OUTPUT_DAE_VECTORS_NOPORT", Util::Param("OUTPUT_DAE_VECTORS_NOPORT", 0)));
  parameters.insert(Util::ParamMap::value_type("OUTPUT_DAE_MATRICES", Util::Param("OUTPUT_DAE_MATRICES", 0)));
  parameters.insert(Util::ParamMap::value_type("OUTPUT_REDUCED_CONDUCTANCES", Util::Param("OUTPUT_REDUCED_CONDUCTANCES", 0)));
  parameters.insert(Util::ParamMap::value_type("OUTPUT_PORT_CURRENTS", Util::Param("OUTPUT_PORT_CURRENTS", 0)));
}

} // namespace Analysis
} // namespace Xyce

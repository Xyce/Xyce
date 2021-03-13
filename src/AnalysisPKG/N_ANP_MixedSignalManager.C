//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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

#include <N_ANP_MixedSignalManager.h>
#include <N_ANP_Transient.h>
#include <N_ERH_Message.h>
#include <N_LOA_Loader.h>
#include <N_NLS_Manager.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_WorkingIntegrationMethod.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : MixedSignalManager::provisionalStep
// Purpose       : Used by mixed-signal.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/04/2009
//-----------------------------------------------------------------------------
bool
MixedSignalManager::provisionalMixedSignalStep(
  const TimeIntg::TIAParams &   tia_params,
  Linear::System &              linear_system,
  Nonlinear::Manager &          nonlinear_manager,
  double                        maxTimeStep,
  double &                      timeStep)
{
  if (!mixedSignalAnalysisObject_)
  {
    if (!getCreatorVector().empty()) {
      mixedSignalAnalysisObject_ = dynamic_cast<Transient *>((*getCreatorVector().front()).create());
      if (!mixedSignalAnalysisObject_)
      {
        Report::UserError() << "Mixed signal requires transient analysis";
        return false;
      }
    }
    else
    {
      Report::DevelFatal().in("MixedSignalManager::provisionalStep") << "unknown type of analysis";
      return false;
    }
    
    mixedSignalAnalysisObject_->init();

    // Start the solvers timer.
    getXyceTranTimer().resetStartTime();
    setPrimaryAnalysisObject(mixedSignalAnalysisObject_);
  }

  bool bsuccess = true;
  bool dcopFlag = true;

  if (mixedSignalAnalysisObject_)
  {
    dcopFlag = mixedSignalAnalysisObject_->getDCOPFlag();
  }

  // Now save time step info, in case this step gets rejected.

  if (!getStepErrorControl().isFinished())
  {
    bool stepSuccess = false;

    if (dcopFlag) // if dcop step, make one attempt.
    {
      mixedSignalAnalysisObject_->mixedSignalStep(maxTimeStep);

      // only call finalize step here if we have failed.
      if (!getStepErrorControl().stepAttemptStatus)
      {
        mixedSignalAnalysisObject_->finalizeMixedSignalStep();
      }
      stepSuccess = getStepErrorControl().stepAttemptStatus;
    }
    else // else, if transient step, keep re-taking the step
         // until it succeeds, or gets an unrecoverable failure,
         // such as time-step-too-small.
    {
      bool recoverableFailureFlag=true;
      while (!stepSuccess && recoverableFailureFlag)
      {
        mixedSignalAnalysisObject_->mixedSignalStep(maxTimeStep);

        // Only call finalize step here if step has failed.
        // If we succeed, we want to give Habanero the opportunity
        // to reject the step, after this function (provisionalStep)
        // exits.
        if (!getStepErrorControl().stepAttemptStatus)
        {
          recoverableFailureFlag = mixedSignalAnalysisObject_->finalizeMixedSignalStep();
        }
        else
        {
          stepSuccess = true;
        }
      }
    }
    bsuccess = stepSuccess;
  }

  // get the step information.
  //
  if (dcopFlag)
  {
    timeStep = 0.0;
  }
  else
  {
    timeStep = getStepErrorControl().currentTimeStep;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MixedSignalManager::acceptMixedSignalProvisionalStep
// Purpose       : Used by mixed-signal.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/2009
//-----------------------------------------------------------------------------
void MixedSignalManager::acceptMixedSignalProvisionalStep()
{
  mixedSignalAnalysisObject_->finalizeMixedSignalStep();
}

//-----------------------------------------------------------------------------
// Function      : MixedSignalManager::rejectMixedSignalProvisionalStep
// Purpose       : Used by mixed-signal.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/2009
//-----------------------------------------------------------------------------
void MixedSignalManager::rejectMixedSignalProvisionalStep(
  Loader::Loader &              loader,
  const TimeIntg::TIAParams &   tia_params)
{
  getStepErrorControl().stepAttemptStatus = false;
  getStepErrorControl().updateBreakPoints(loader, tia_params.initialTime);

  bool dcopFlag = false;
  if (mixedSignalAnalysisObject_)
  {
    dcopFlag = mixedSignalAnalysisObject_->getDCOPFlag();
  }

  if (dcopFlag)
  {
    mixedSignalAnalysisObject_->finalizeMixedSignalStep();
  }

  else // Transient
  {
    loader.stepFailure(getTwoLevelMode());
    getWorkingIntegrationMethod().rejectStepForHabanero();

    mixedSignalAnalysisObject_->stats_.failedStepsAttempted_  += 1;
    getStepErrorControl().numberSuccessiveFailures += 1;
  }
}

} // namespace Analysis
} // namespace Xyce

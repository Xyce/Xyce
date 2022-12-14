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

//-----------------------------------------------------------------------------
// Purpose       : DC Sweep class analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL
// Creation Date : 01/24/08
//
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_ANP_AnalysisManager.h>
#include <N_ANP_DCSweep.h>
#include <N_ANP_HB.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_SweepParam.h>
#include <N_ANP_SweepParamFreeFunctions.h>
#include <N_ANP_OutputConductanceFile.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>
#include <N_IO_InitialConditions.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_LOA_Loader.h>
#include <N_NLS_Manager.h>
#include <N_NLS_fwd.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_fwd.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_WorkingIntegrationMethod.h>
#include <N_TOP_Topology.h>
#include <N_LAS_System.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_Factory.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_NoCase.h>
#include <N_UTL_OptionBlock.h>
#include <N_ERH_Message.h>

#include <Teuchos_RCP.hpp>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : DCSweep::DCSweep
// Purpose       :
// Special Notes : 
// Scope         : 
// Creator       : Eric R. Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
DCSweep::DCSweep(
  AnalysisManager &                     analysis_manager,
  //Linear::System &                      linear_system,
  Linear::System *                      linear_system_ptr,
  Nonlinear::Manager &                  nonlinear_manager,
  Loader::Loader &                      loader,
  Topo::Topology &                      topology,
  IO::InitialConditionsManager &        initial_conditions_manager,
  HB *                                  hb_analysis)
  : AnalysisBase(analysis_manager, "DC Sweep"),
    analysisManager_(analysis_manager),
    loader_(loader),
    linearSystemPtr_(linear_system_ptr),
    nonlinearManager_(nonlinear_manager),
    topology_(topology),
    initialConditionsManager_(initial_conditions_manager),
    outputManagerAdapter_(analysis_manager.getOutputManagerAdapter()),
    hbAnalysis_(hb_analysis),
    tiaParams_(),
    sensFlag_(analysis_manager.getSensFlag()),
    dcLoopInitialized_(false),
    condTestFlag_(false),
    condTestDeviceNames_(),
    numSensParams_(0),
    dcLoopSize_(0)
{}

//-----------------------------------------------------------------------------
// Function      : DCSweep::~DCSweep
// Purpose       :
// Special Notes : 
// Scope         : 
// Creator       : Eric R. Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
DCSweep::~DCSweep()
{}

//-----------------------------------------------------------------------------
// Function      : DCSweep::setTranOptions
// Purpose       : Process '.options timeint' parameters
//
// Special Notes : Most time integration options are not relevant to the DC sweep class.  
//
//                 But, if one has a transient circuit that uses some of these, and then 
//                 changes the analysis in that circuit from .TRAN to .DC, but leaves the 
//                 options line in place, it makes sense that Xyce should not throw a 
//                 fatal error.  So, all the valid options should be included in 
//                 the if-then-else block below.
//                 
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool DCSweep::setTimeIntegratorOptions(const Util::OptionBlock & option_block)
{
  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; ++it)
  {
    const Util::Param &param = (*it);

    if (param.uTag() == "METHOD") { ; }
    else if (param.uTag()=="EXITTIME" ) { ; }
    else if (param.uTag()=="EXITSTEP" ) { ; }
    else if (param.uTag() == "HISTORYTRACKINGDEPTH" ) { ; }
    else if (param.uTag() == "PASSNLSTALL") { ; }
    else if (param.uTag() == "CONDTEST")
    {
      condTestFlag_ = static_cast<bool> (param.getImmutableValue<int>());
    }
    else if (std::string( param.uTag() ,0,18) == "CONDTESTDEVICENAME")
    {
      condTestDeviceNames_.push_back(param.stringValue() );
    }
    else if (param.uTag() == "DAESTATEDERIV")
    {
      analysisManager_.setDAEStateDerivFlag(static_cast<bool> (param.getImmutableValue<int>()));
    }
    else if (param.uTag() == "DEBUGLEVEL")
    {
      IO::setTimeIntegratorDebugLevel(analysisManager_.getCommandLine(), param.getImmutableValue<int>());
    }
    else if ( std::string( param.uTag() ,0,11) == "BREAKPOINTS") { ; }
    else if (nonlinearManager_.setReturnCodeOption(param)) { ; }
    else if (tiaParams_.setTimeIntegratorOption(param)) { ; }
    else if (setDCOPOption(param)) { ; }
    else
    {
      Report::UserError() << param.uTag() << " is not a recognized time integration option";
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::setAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 6/22/10
//-----------------------------------------------------------------------------
bool DCSweep::setAnalysisParams(const Util::OptionBlock & paramsBlock)
{
  if (isDataSpecified(paramsBlock))
  {
    // This handle the case of having multiple .DC lines in the netlist, of
    // which only some might use DATA=<tableName>
    dataSpecification_ = true;
  }
  dcSweepVector_.push_back(parseSweepParams(paramsBlock.begin(), paramsBlock.end()));
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DCPSweep::setDataStatements
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/5/18
//-----------------------------------------------------------------------------
bool DCSweep::setDataStatements(const Util::OptionBlock & paramsBlock)
{
  return processDataStatements(paramsBlock, dataNamesMap_, dataTablesMap_);
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::convertDataToSweepParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/5/18
//-----------------------------------------------------------------------------
bool  DCSweep::convertDataToSweepParams()
{
  return convertData( dcSweepVector_, dataNamesMap_, dataTablesMap_);
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::outputFailureStats
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/21/2010
//-----------------------------------------------------------------------------
bool DCSweep::outputFailureStats(std::ostream &os)
{
  if (!(dcSweepFailures_.empty()))
  {
    os << "\tFailed DC sweep steps:\t\t" << std::endl;

    for (std::vector<int>::iterator iter = dcSweepFailures_.begin(); iter != dcSweepFailures_.end(); ++iter)
    {
      os << "\t\tDC Step # " << *iter << std::endl;
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::finalExpressionBasedSetup()
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/4/2021
//-----------------------------------------------------------------------------
void DCSweep::finalExpressionBasedSetup()
{
  if (sensFlag_)
  {
    Stats::StatTop _sensitivityStat("Sensitivity");

    nonlinearManager_.enableSensitivity(
        *analysisManager_.getDataStore(), 
        analysisManager_.getStepErrorControl(),
        *analysisManager_.getPDSManager(), 
        topology_, 
        outputManagerAdapter_.getOutputManager(),
        numSensParams_);
  }
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::doRun()
// Purpose       : This is the main controlling loop for DC sweep analysis.
//                 This loop calls a series of operating point calculations
//                 (ie calculations in which there is no time integration,
//                 and the time integration method is always set to "none").
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool DCSweep::doRun()
{
  return doInit() && doLoopProcess() && doFinish();
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::doInit()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool DCSweep::doInit()
{
  bool bsuccess = true;

  // check if the "DATA" specification was used.  If so, create a new vector of 
  // SweepParams, in the "TABLE" style.
  if (dataSpecification_ && !convertDataToSweepParams())
  {
    Report::UserFatal() << "Invalid data=<name> parameter on .DC line.";
    return false;
  }

  // set up the various sweep variables
  // the following should only happen once, but must be isolated to
  // prevent a step loop outside of this dc sweep from causing to
  // occur more often
  if( !dcLoopInitialized_ )
  {
    if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
    {
      Xyce::dout() << std::endl << std::endl
                   << Xyce::subsection_divider << std::endl
                   << "DCSweep::run()" << std::endl;
    }

    dcLoopSize_ = setupSweepLoop(analysisManager_.getComm(), loader_, dcSweepVector_.begin(), dcSweepVector_.end());

    outputManagerAdapter_.setDCAnalysisMaxSteps( dcLoopSize_ );
    outputManagerAdapter_.setDCSweepVector(dcSweepVector_);

    dcLoopInitialized_ = true;
  }

  //setup for operating pt calculation
  baseIntegrationMethod_ = TimeIntg::methodsEnum::NO_TIME_INTEGRATION;
  analysisManager_.createTimeIntegratorMethod(tiaParams_, baseIntegrationMethod_);

  stepNumber = 0;
  setDoubleDCOPEnabled(loader_.isPDESystem());
  if (getDoubleDCOPEnabled() && getDoubleDCOPStep() == 0)
  {
    nonlinearManager_.setAnalysisMode(nonlinearAnalysisMode(ANP_MODE_DC_NLPOISSON));
  }

  initializeSolution_();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::getDCSweepVec
// Purpose       : used to provide the DC Sweep Vector to remeasure. 
//                 
// Special Notes : For "normal" (not remeasured) simulation runs, it is more 
//                 typical to use the getDCSweepVector() function from the 
//                 OutputMgrAdapter class to get the current information in 
//                 the DC Sweep Vector.  So, this function may only be for
//                 use with remeasure.
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystem Modeling
// Creation Date : 6/27/2017
//-----------------------------------------------------------------------------
SweepVector DCSweep::getDCSweepVec()
{
  return dcSweepVector_;
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::initializeSolution()
// Purpose       : Move solution-initialization steps to separate routine so
//                 the process may be repeated cleanly when necessary.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystem Modeling
// Creation Date : 2/17/2010
//-----------------------------------------------------------------------------
void DCSweep::initializeSolution_()
{
  // set initial guess, if there is one to be set.
  // this setInitialGuess call is to up an initial guess in the
  // devices that have them (usually PDE devices).  This is different than
  // the "intializeProblem" call, which sets IC's.  (initial conditions are
  // different than initial guesses.
  loader_.setInitialGuess(analysisManager_.getDataStore()->nextSolutionPtr);
  if ( !hbAnalysis_|| hbAnalysis_->getHBtranFlags())

  // If available, set initial solution (.IC, .NODESET, etc).
    setInputOPFlag(
      initialConditionsManager_.setupInitialConditions(outputManagerAdapter_.getComm(),
                                                     topology_.getSolutionNodeNameMap(),
						     outputManagerAdapter_.getAliasNodeMap(),
                                                     *analysisManager_.getDataStore()->nextSolutionPtr,
                                                     *linearSystemPtr_));

  // Set a constant history for operating point calculation
  analysisManager_.getDataStore()->setConstantHistory();
  analysisManager_.getWorkingIntegrationMethod().obtainCorrectorDeriv();
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::loopProcess()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool DCSweep::doLoopProcess()
{
  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::INITIALIZE, AnalysisEvent::DC));

  int currentStep = 0;
  int finalStep = dcLoopSize_;
  while (currentStep < finalStep)
  {
    outputManagerAdapter_.setDCAnalysisStepNumber(currentStep);

    if (VERBOSE_TIME)
    {
      printStepHeader(Xyce::dout());
      analysisManager_.getStepErrorControl().outputTimeInfo(Xyce::dout());
    }

    bool reset = updateSweepParams(loader_, currentStep, dcSweepVector_.begin(), dcSweepVector_.end(), true);

    // Tell the manager if any of our sweeps are being reset in this loop iteration.
    analysisManager_.setSweepSourceResetFlag(reset);

    outputManagerAdapter_.setDCSweepVector(dcSweepVector_);

    if (currentStep != 0 && reset)
    {
      analysisManager_.getDataStore()->setZeroHistory();
      initializeSolution_();
    }

    // Perform the step:
    static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_STARTED, AnalysisEvent::DC, 0.0, currentStep));

    takeStep_();

    // Set things up for the next time step, based on if this one was
    // successful.
    if (analysisManager_.getStepErrorControl().stepAttemptStatus)
    {
      static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_SUCCESSFUL, AnalysisEvent::DC, 0.0, currentStep));
      doProcessSuccessfulStep();
    }
    else // stepAttemptStatus  (ie do this if the step FAILED)
    {
      static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_FAILED, AnalysisEvent::DC, 0.0, currentStep));
      doProcessFailedStep();
    }

    // we don't really control the loop counter here.  processSuccessfulStep
    // and processFailedStep do that work.  This needs to be cleaned up. RLS
    currentStep = stepNumber;
  } // end of sweep loop

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::FINISH, AnalysisEvent::DC));

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::processSuccessfulStep()
//
// Purpose       : Used by both function dcSweepLoop and 2-level solve
//                 function calls to process successful DC steps.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool DCSweep::doProcessSuccessfulStep()
{
  Stats::StatTop _processSuccessfulStepStat("Successful Step");
  Stats::TimeBlock _processSuccessfulStepTimer(_processSuccessfulStepStat);

  loader_.stepSuccess(TWO_LEVEL_MODE_DC_SWEEP); // analysisManager_.getTwoLevelMode());

  // This output call is for device-specific output, such as .OP,
  // or internal plot output from PDE(TCAD) devices.
  loader_.outputPlotFiles();

  if (sensFlag_ && !firstDoubleDCOPStep() )
  {
    nonlinearManager_.calcSensitivity(objectiveVec_, dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);
  }

  // Do some statistics, as long as this isn't the first "double"
  // DCOP step. (that one doesn't count)
  if ( !firstDoubleDCOPStep() )
  {
    stepNumber += 1;
    stats_.successStepsThisParameter_ += 1;
    stats_.successfulStepsTaken_ += 1;
  }

  // update the data arrays, output:
  analysisManager_.getDataStore()->updateSolDataArrays();

  dcSweepOutput();

  // now that output has been called, update the doubleDCOP step
  // if neccessary. (pde-only)
  nextDCOPStep();
  nonlinearManager_.setAnalysisMode(nonlinearAnalysisMode(ANP_MODE_DC_SWEEP));

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::processFailedStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool DCSweep::doProcessFailedStep()
{
  Stats::StatTop _processFailedStat("Failed Steps");
  Stats::TimeBlock _processFailedTimer(_processFailedStat);

  loader_.stepFailure(analysisManager_.getTwoLevelMode());

  stepNumber += 1;
  dcSweepFailures_.push_back(stepNumber);
  stats_.failedStepsAttempted_  += 1;
  analysisManager_.getStepErrorControl().numberSuccessiveFailures += 1;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::doFinish()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool DCSweep::doFinish()
{
  bool bsuccess = true;

  if (DEBUG_ANALYSIS)
    Xyce::dout() << "Calling DCSweep::doFinish() outputs!" << std::endl;

  outputManagerAdapter_.finishOutput();
  if (!(dcSweepFailures_.empty()))
  {
    bsuccess = false;
  }
  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : DCSweep::doHandlePredictor
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 06/24/2013
//-----------------------------------------------------------------------------
bool DCSweep::doHandlePredictor()
{
  analysisManager_.getDataStore()->setErrorWtVector(tiaParams_, topology_.getVarTypes());
  analysisManager_.getWorkingIntegrationMethod().obtainPredictor();
  analysisManager_.getWorkingIntegrationMethod().obtainPredictorDeriv();

  // In case this is the upper level of a 2-level sim, tell the
  // inner solve to do its prediction:
  bool        beginIntegrationFlag = analysisManager_.getBeginningIntegrationFlag();
  double      nextTimeStep = analysisManager_.getStepErrorControl().currentTimeStep;
  double      nextTime = analysisManager_.getStepErrorControl().nextTime;          
  int         currentOrder = analysisManager_.getWorkingIntegrationMethod().getOrder();

  loader_.startTimeStep(beginIntegrationFlag, nextTimeStep, nextTime, currentOrder);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::takeStep_
// Purpose       : Take a DC Sweep integration step.
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, SNL
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
void DCSweep::takeStep_()
{
  // Integration step predictor
  doHandlePredictor();

  // Load B/V source devices with time data
  loader_.updateSources();

  // Nonlinear solve
  {
    Stats::StatTop _nonlinearSolveStat("Solve");
    Stats::TimeBlock _nonlinearSolveTimer(_nonlinearSolveStat);

    analysisManager_.getStepErrorControl().newtonConvergenceStatus = nonlinearManager_.solve();
  }

  // Add change to solution
  analysisManager_.getWorkingIntegrationMethod().stepLinearCombo();

  gatherStepStatistics(stats_, nonlinearManager_.getNonlinearSolver(), analysisManager_.getStepErrorControl().newtonConvergenceStatus);

  analysisManager_.getStepErrorControl().evaluateStepError(loader_, tiaParams_);
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::twoLevelStep
//
// Purpose       : Used by 2-level Newton solves to execute a single DC sweep
//                 step.
//
// Special Notes : This is mostly what happens on the inner loop of
//                 dcSweepLoop, except that DC parameters are not updated,
//                 and success/failure of the step is not determined.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool DCSweep::twoLevelStep()
{
  loader_.updateSources();
  analysisManager_.getStepErrorControl().newtonConvergenceStatus = nonlinearManager_.solve();
  analysisManager_.getWorkingIntegrationMethod().stepLinearCombo();
  gatherStepStatistics(stats_, nonlinearManager_.getNonlinearSolver(), analysisManager_.getStepErrorControl().newtonConvergenceStatus);
  analysisManager_.getStepErrorControl().evaluateStepError(loader_, tiaParams_);

  return analysisManager_.getStepErrorControl().stepAttemptStatus;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void DCSweep::outputDAEvectors()
{
  nonlinearManager_.outputDAEvectors();
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void DCSweep::outputDAEmatrices()
{
  nonlinearManager_.outputDAEmatrices();
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::dcSweepOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 1/31/08
//-----------------------------------------------------------------------------
void DCSweep::dcSweepOutput()
{
  if (hbAnalysis_)
  {
    std::vector<double> timePoints, freqPoints;
    Teuchos::RCP<Linear::BlockVector> timeDomainSolnVec;
    Teuchos::RCP<Linear::BlockVector> freqDomainSolnVecReal;
    Teuchos::RCP<Linear::BlockVector> freqDomainSolnVecImaginary;
    Teuchos::RCP<Linear::BlockVector> timeDomainLeadCurrentVec;
    Teuchos::RCP<Linear::BlockVector> freqDomainLeadCurrentVecReal;
    Teuchos::RCP<Linear::BlockVector> freqDomainLeadCurrentVecImaginary;
    Teuchos::RCP<Linear::BlockVector> timeDomainJunctionVoltageVec;
    Teuchos::RCP<Linear::BlockVector> freqDomainJunctionVoltageVecReal;
    Teuchos::RCP<Linear::BlockVector> freqDomainJunctionVoltageVecImaginary;

    // RLS Todo:  Need to update HB call here 
    // for leadcurrent and junctionvoltage vector 


    hbAnalysis_->prepareHBOutput(
      *(analysisManager_.getDataStore()->currSolutionPtr),
      timePoints,
      freqPoints,
      timeDomainSolnVec,
      freqDomainSolnVecReal,
      freqDomainSolnVecImaginary,
      timeDomainLeadCurrentVec,
      freqDomainLeadCurrentVecReal,
      freqDomainLeadCurrentVecImaginary,
      timeDomainJunctionVoltageVec,
      freqDomainJunctionVoltageVecReal,
      freqDomainJunctionVoltageVecImaginary);

    // output associated with HB frequency output, such as .PRINT HB_FD lines
    outputManagerAdapter_.outputHB_FD(
      freqPoints,
      *freqDomainSolnVecReal,
      *freqDomainSolnVecImaginary,
      *freqDomainLeadCurrentVecReal,
      *freqDomainLeadCurrentVecImaginary,
      *freqDomainJunctionVoltageVecReal,
      *freqDomainJunctionVoltageVecImaginary);

    // output associated with HB time-domain output, such as .PRINT HB_TD lines.
    // This function is not used by .PRINT HB_STARTUP or .PRINT HB_IC lines though.
    outputManagerAdapter_.outputHB_TD(
      timePoints,
      *timeDomainSolnVec,
      *timeDomainLeadCurrentVec,
      *timeDomainJunctionVoltageVec);
  }

  // Make sure this isn't the NL_POISSON step of a PDE DCOP.
  else if (!firstDoubleDCOPStep())
  {
    // conventional .PRINT output
    outputManagerAdapter_.dcOutput(
      stepNumber,
      *(analysisManager_.getDataStore()->currSolutionPtr),
      *(analysisManager_.getDataStore()->currStatePtr),
      *(analysisManager_.getDataStore()->currStorePtr),
      *(analysisManager_.getDataStore()->currLeadCurrentPtr),
      *(analysisManager_.getDataStore()->currLeadDeltaVPtr),
      *(analysisManager_.getDataStore()->currLeadCurrentQDerivPtr),
      objectiveVec_,
      dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);

    // output for SAVE
    // outputManagerAdapter_.outputDCOP(*analysisManager_.getDataStore()->currSolutionPtr);
    initialConditionsManager_.outputDCOP(outputManagerAdapter_.getComm(), topology_.getSolutionNodeNameMap(), *analysisManager_.getDataStore()->currSolutionPtr);
  
    // ERK. embedded sampling call back.  This is (probably) a hack until a better method for
    // handling outputs for embedded sampling is devised.
    if ( !(parentAnalysisPtrVec_.empty()) )
    {
      std::vector < AnalysisBase * >::iterator iter = parentAnalysisPtrVec_.begin();
      std::vector < AnalysisBase * >::iterator end = parentAnalysisPtrVec_.end();
      for ( ; iter!=end; ++iter)
      {
        (*iter)->stepCallBack();
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::printStepHeader()
// Purpose       : Prints out time step information.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/26/00
//-----------------------------------------------------------------------------
void DCSweep::printStepHeader(std::ostream &os)
{
  if (VERBOSE_TIME)
    dout() << std::endl << std::endl
           << "***** "<< (DEBUG_ANALYSIS ? analysisManager_.getNetlistFilename() : "")
           << "  Start of DCOP STEP                        # " << stepNumber+1
           << std::endl << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::printLoopInfo
// Purpose       :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 12/14/2010
//-----------------------------------------------------------------------------
bool DCSweep::printLoopInfo(int start, int finish)
{
  bool bsuccess = AnalysisBase::printLoopInfo(start, finish);
  if (start == 0 && finish == 0)
  {
    outputFailureStats(Xyce::lout());
  }
  return bsuccess;
}

namespace {

typedef Util::Factory<AnalysisBase, DCSweep>  DCSweepFactoryBase;

//-----------------------------------------------------------------------------
// Class         : DCSweepFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jan 29 12:53:02 2015
//-----------------------------------------------------------------------------
///
/// Factory for parsing DCSweep parameters from the netlist and creating DCSweep analysis.
///
class DCSweepFactory : public DCSweepFactoryBase
{
public:
  //-----------------------------------------------------------------------------
  // Function      : DCSweepFactory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:54:09 2015
  //-----------------------------------------------------------------------------
  ///
  /// Constructs the DCSweep analysis factory
  ///
  /// @invariant Stores the results of parsing.  Multiple DCSweep analysis options may be
  /// applied and each generates and additional dcsweep.
  ///
  /// @invariant The existence of the parameters specified in the constructor cannot
  /// change.
  ///
  /// @param analysis_manager 
  /// @param linear_system 
  /// @param nonlinear_manager 
  /// @param topology 
  ///
  DCSweepFactory(
    Analysis::AnalysisManager &         analysis_manager,
    Linear::System *                    linear_system_ptr,
    Nonlinear::Manager &                nonlinear_manager,
    Loader::Loader &                    loader,
    Topo::Topology &                    topology,
    IO::InitialConditionsManager &      initial_conditions_manager)
    : DCSweepFactoryBase(),
      analysisManager_(analysis_manager),
      linearSystemPtr_(linear_system_ptr),
      nonlinearManager_(nonlinear_manager),
      loader_(loader),
      topology_(topology),
      initialConditionsManager_(initial_conditions_manager),
      dcSweepAnalysisOptionBlock_(),
      timeIntegratorOptionBlock_()
  {}

  virtual ~DCSweepFactory()
  {}

  //-----------------------------------------------------------------------------
  // Function      : create
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:59:00 2015
  //-----------------------------------------------------------------------------
  ///
  /// Create a new DCSweep analysis and applies the analysis and time integrator option blocks.
  ///
  /// @return new DCSweep analysis object
  ///
  ///
  DCSweep *create() const
  {
    analysisManager_.setAnalysisMode(ANP_MODE_DC_SWEEP);
    //DCSweep *dc_sweep = new DCSweep(analysisManager_, linearSystem_, nonlinearManager_, loader_, topology_, initialConditionsManager_);
    DCSweep *dc_sweep = new DCSweep(analysisManager_, linearSystemPtr_, nonlinearManager_, loader_, topology_, initialConditionsManager_);
    for (std::vector<Util::OptionBlock>::const_iterator it = dcSweepAnalysisOptionBlock_.begin(), end = dcSweepAnalysisOptionBlock_.end(); it != end; ++it)
      dc_sweep->setAnalysisParams(*it);
    dc_sweep->setTimeIntegratorOptions(timeIntegratorOptionBlock_);

    for (std::vector<Util::OptionBlock>::const_iterator it = dataOptionBlockVec_.begin(), end = dataOptionBlockVec_.end(); it != end; ++it)
    {
      dc_sweep->setDataStatements(*it);
    }



    return dc_sweep;
  }

  //-----------------------------------------------------------------------------
  // Function      : setDCSweepAnalysisOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 13:00:14 2015
  //-----------------------------------------------------------------------------
  ///
  /// Saves the analysis parsed options block in the factory.
  ///
  /// @invariant Appends to any previously specified analysis option block.
  ///
  /// @param option_block parsed option block
  ///
  void setDCSweepAnalysisOptionBlock(const Util::OptionBlock &option_block)
  {
    for (std::vector<Util::OptionBlock>::iterator it = dcSweepAnalysisOptionBlock_.begin(), end = dcSweepAnalysisOptionBlock_.end(); it != end; ++it)
    {
      if (Util::compareParamLists(option_block, *it))
      {
        (*it) = option_block;
        return;
      }
    }

    // save the new one.
    dcSweepAnalysisOptionBlock_.push_back(option_block); // save a copy for later.
  }

  //-----------------------------------------------------------------------------
  // Function      : setTimeIntegratorOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 13:01:27 2015
  //-----------------------------------------------------------------------------
  ///
  /// Saves the time integrator parsed option block.
  ///
  /// @invariant Overwrites any previously specified time integrator option block.
  ///
  /// @param option_block parsed option block
  ///
  bool setTimeIntegratorOptionBlock(const Util::OptionBlock &option_block)
  {
    timeIntegratorOptionBlock_ = option_block;

    return true;
  }

  //-----------------------------------------------------------------------------
  bool setDotDataBlock(const Util::OptionBlock &option_block)
  {
    dataOptionBlockVec_.push_back(option_block);
    return true;
  }

public:
  AnalysisManager &                     analysisManager_;
  Linear::System *                      linearSystemPtr_;
  Nonlinear::Manager &                  nonlinearManager_;
  Loader::Loader &                      loader_;
  Topo::Topology &                      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;

private:
  std::vector<Util::OptionBlock>        dcSweepAnalysisOptionBlock_;
  Util::OptionBlock                     timeIntegratorOptionBlock_;
  std::vector<Util::OptionBlock>        dataOptionBlockVec_;
};


// .DC
struct DCSweepAnalysisReg : public IO::PkgOptionsReg
{
  DCSweepAnalysisReg(
    DCSweepFactory &          factory)
    : factory_(factory)
  {}

  bool operator()(const Util::OptionBlock &option_block)
  {
    factory_.setDCSweepAnalysisOptionBlock(option_block);

    factory_.analysisManager_.addAnalysis(&factory_);

    return true;
  }

  DCSweepFactory &            factory_;
};

//-----------------------------------------------------------------------------
// Function      : extractDCData
// Purpose       : Determine number of sweep variables on this DC line
//                : and create option blocks for each, storing the blocks
//                : in the referenced vector.
//-----------------------------------------------------------------------------
bool extractDCData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  std::vector<Util::OptionBlock> option_block_vec = 
    extractDCDataInternals("DC", options_manager, netlist_filename, parsed_line); 

  if (option_block_vec.size())
  {
    std::vector<Util::OptionBlock>::iterator it = option_block_vec.begin();
    std::vector<Util::OptionBlock>::const_iterator it_end = option_block_vec.end();
    for ( ; it != it_end; ++it ) 
      circuit_block.addOptions(*it);
  }
  else
   return false;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : extractOPData
// Purpose       : Extract the parameters from a netlist .OP line held in
//                 parsed_line.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/05/2001
//-----------------------------------------------------------------------------
bool
extractOPData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("OP", Util::OptionBlock::NO_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

  int numFields = parsed_line.size();

  // Check that the minimum required number of fields are on the line.
  if ( numFields > 1 )
  {
    Report::UserWarning0().at(netlist_filename, parsed_line[0].lineNumber_) << "Ignoring extra fields on .OP line";
  }

  circuit_block.addOptions(option_block);

  return true;
}

} // namespace <unnamed>


namespace {
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
struct isTokenString 
{
  isTokenString (std::string & test) : testString(test) {};
  bool operator() (const IO::StringToken & t1)
  {
    return compare_nocase(t1.string_.c_str(), testString.c_str()) == 0;
  }
  std::string & testString;
};
}

//-----------------------------------------------------------------------------
// Function      : processDCblock
// Purpose       : helper function for extractDCDataInternals
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 12/14/2022
//-----------------------------------------------------------------------------
void processDCblock(
    int & linePosition,
    int numFields,
    Util::OptionBlock & option_block,
    const std::string & type,
    const std::string & nextStr,
    const std::string & currStr,
    const std::string & netlist_filename,
    const IO::TokenVector & parsed_line)
{
  option_block.addParam(Util::Param("TYPE", type));
  option_block.addParam(Util::Param("PARAM", currStr == type ? nextStr : currStr)); linePosition+=2;
  if( numFields <= linePosition + 2 )
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_) << 
      ".DC line not formatted correctly, found unexpected number of fields";
    linePosition = numFields;
  }
  else 
  {
    option_block.addParam(Util::Param("START", parsed_line[linePosition++].string_));
    option_block.addParam(Util::Param("STOP", parsed_line[linePosition++].string_));
    option_block.addParam(Util::Param("NUMSTEPS", parsed_line[linePosition++].string_));
  }
}

//-----------------------------------------------------------------------------
// Function      : extractDCDataInternals
// Purpose       : Determine number of sweep variables on this DC line
//                : and create option blocks for each, storing the blocks
//                : in the referenced vector.
//-----------------------------------------------------------------------------
std::vector<Util::OptionBlock>
extractDCDataInternals(
  const std::string &           name,
  IO::PkgOptionsMgr &           options_manager,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  std::vector<Util::OptionBlock> option_block_vec;
 
  // length of the original .DC line
  int numFields = parsed_line.size();

  if (numFields < 4)
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << ".DC line has an unexpected number of fields";
    return option_block_vec;
  }

  // number of sweep sources on this line, and index to current source
  int sourcesFound = 0;

  // start of parameters (skip over the ".DC")
  int linePosition = 1;

  // check for "DATA" first.  If DATA is found, then skip everything else ----
  std::string tmp = std::string("DATA");
  IO::TokenVector::const_iterator startPL = parsed_line.begin();  startPL++;
  IO::TokenVector::const_iterator endPL = parsed_line.end();
  IO::TokenVector::const_iterator iter = std::find_if(startPL, endPL, isTokenString(tmp)); 
  if (iter != parsed_line.end())
  {
    int dataPos = std::distance(parsed_line.begin(),iter);
    if (numFields != 4)
    {
      Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
        << ".DC line not formatted correctly.  numFields = " << numFields;
      return option_block_vec;
    }

    Util::OptionBlock option_block(name, Util::OptionBlock::ALLOW_EXPRESSIONS, netlist_filename, parsed_line[linePosition].lineNumber_);
    Util::Param parameter("", "");
    parameter.setTag( "TYPE" );
    parameter.setVal( "DATA" );
    option_block.addParam( parameter );

    parameter.setTag( "DATASET" );
    parameter.setVal( parsed_line[ dataPos+2 ].string_ );
    option_block.addParam( parameter );

    option_block_vec.push_back( option_block );

    return option_block_vec;
  }
  // End of the DATA block ----

  // line can be variable in length, with muliple sources to a line
  linePosition = 1;
  while( linePosition < numFields )
  {
    Util::OptionBlock option_block(name, Util::OptionBlock::ALLOW_EXPRESSIONS, netlist_filename, parsed_line[linePosition].lineNumber_);

    if (linePosition + 1 == numFields) {
      Report::UserError0().at(netlist_filename, parsed_line[linePosition].lineNumber_) << "Extraneous values on .DC line";
      break;
    }

    std::string currStr = parsed_line[linePosition].string_;
    Util::toUpper(currStr);

    std::string nextStr = parsed_line[linePosition + 1].string_;
    Util::toUpper(nextStr);

    if ( (nextStr == "LIST" || currStr == "LIST") )
    {
      // sweep type is LIST, get sweep variable name and move to just before beginning of list
      option_block.addParam(Util::Param("TYPE", "LIST"));
      option_block.addParam(Util::Param("PARAM", currStr == "LIST" ? nextStr : currStr));
      ++linePosition;

      // collect values in this list until next sweep variable name is found
      while (++linePosition < numFields
          && (Util::isValue(parsed_line[linePosition].string_) || parsed_line[linePosition].string_[0]=='{') )
      {
        option_block.addParam(Util::Param("VAL", parsed_line[linePosition].string_));
      }
    }
    else if ( (nextStr == "DEC" || currStr == "DEC") )
    {
      processDCblock( linePosition, numFields, option_block,"DEC",nextStr,currStr, netlist_filename, parsed_line);
    }
    else if ( (nextStr == "OCT" || currStr == "OCT") )
    {
      processDCblock( linePosition, numFields, option_block,"OCT",nextStr,currStr, netlist_filename, parsed_line);
    }
    else if ( (nextStr == "LIN" || currStr == "LIN") )
    {
      processDCblock( linePosition, numFields, option_block,"LIN",nextStr,currStr, netlist_filename, parsed_line);
    }
    else // no explicit type given, meaning LIN:
    {
      option_block.addParam(Util::Param("TYPE", "LIN"));
      if( numFields <= linePosition + 3 )
      {
        Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_) << 
          ".DC line not formatted correctly, found unexpected number of fields";
        linePosition = numFields;
      }
      else 
      {
        option_block.addParam(Util::Param("PARAM", parsed_line[linePosition++].string_));
        option_block.addParam(Util::Param("START", parsed_line[linePosition++].string_));
        option_block.addParam(Util::Param("STOP", parsed_line[linePosition++].string_));
      }
      option_block.addParam(Util::Param("STEP", parsed_line[linePosition++].string_));
    }

    option_block_vec.push_back( option_block );

    // record this source (and move on to the next)
    ++sourcesFound;
  }

  return option_block_vec;
}

bool registerDCSweepFactory(FactoryBlock & factory_block)
{
  //DCSweepFactory *factory = new DCSweepFactory(factory_block.analysisManager_, factory_block.linearSystem_, factory_block.nonlinearManager_, factory_block.loader_, factory_block.topology_, factory_block.initialConditionsManager_);
  DCSweepFactory *factory = new DCSweepFactory(factory_block.analysisManager_, &(factory_block.linearSystem_), factory_block.nonlinearManager_, factory_block.loader_, factory_block.topology_, factory_block.initialConditionsManager_);

  addAnalysisFactory(factory_block, factory);

  factory_block.optionsManager_.addCommandParser(".DC", extractDCData);
  factory_block.optionsManager_.addCommandParser(".OP", extractOPData);

  factory_block.optionsManager_.addCommandProcessor("DC", new DCSweepAnalysisReg(*factory));

  factory_block.optionsManager_.addOptionsProcessor("TIMEINT", IO::createRegistrationOptions(*factory, &DCSweepFactory::setTimeIntegratorOptionBlock));

  factory_block.optionsManager_.addCommandProcessor("DATA", 
      IO::createRegistrationOptions(*factory, &DCSweepFactory::setDotDataBlock) );

  return true;
}

} // namespace Analysis
} // namespace Xyce

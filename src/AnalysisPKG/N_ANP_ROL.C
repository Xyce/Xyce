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
// Purpose       : .ROL class analysis functions.
// Special Notes :
// Creator       : 
// Creation Date : 
//
//-----------------------------------------------------------------------------
#include <Xyce_config.h>
 
#include <N_ANP_AnalysisManager.h>
#include <N_ANP_ROL.h>

#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_SweepParam.h>
#include <N_ANP_SweepParamFreeFunctions.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>
#include <N_IO_InitialConditions.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_LOA_Loader.h>
#include <N_NLS_Manager.h> // TT: was not included
#include <N_NLS_fwd.h> // TT: was not included
#include <N_TIA_DataStore.h>
#include <N_TIA_fwd.h> // TT: was not included
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_WorkingIntegrationMethod.h> // TT: was not included

#include <N_UTL_Diagnostic.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_Factory.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_ERH_Message.h>

#include <N_TOP_Topology.h>
#include <N_LAS_System.h>
#include <N_LAS_Builder.h>


#include <N_ANP_StepEvent.h>
//#include <N_ANP_ROLEvent.h>
#include <N_PDS_Manager.h>
#include <N_UTL_Timer.h>

#include <N_LOA_NonlinearEquationLoader.h>

// ROL includes
// Trilinos ROL includes
#include <N_ANP_DCSweep.h>

#ifdef Xyce_ROL

#ifndef OBJTYPE
#define OBJTYPE 1 // 0 - L2Norm, 1 - amplifier circuit
#endif

#ifndef UQ
#define UQ 0
#endif

#include "ROL_StdVector.hpp"
#include "ROL_Vector.hpp"

#include "ROL_Algorithm.hpp"
#include "ROL_CompositeStep.hpp"
#include "ROL_LineSearchStep.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "ROL_BoundConstraint.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Types.hpp"
#include "ROL_BundleStep.hpp"
#include "ROL_BundleStatusTest.hpp"

#include "ROL_XyceVector.hpp"

#include <N_ANP_ROL_DC_Optimization.h>
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_ParametrizedObjective_SimOpt.hpp"
#include "ROL_Reduced_ParametrizedObjective_SimOpt.hpp"

#if UQ==1
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "ROL_RiskNeutralObjective.hpp"
#include "ROL_RiskAverseObjective.hpp"
// ROL sample generators
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_QuasiMonteCarloGenerator.hpp"
#include "ROL_StdTeuchosBatchManager.hpp"
#include "ROL_UserInputGenerator.hpp"
// ROL CVaR definitions
#include "ROL_PlusFunction.hpp"
#include "ROL_CVaR.hpp"
#include "ROL_CVaRVector.hpp"
#include "ROL_CVaRBoundConstraint.hpp"

#endif

#endif

// std includes
#include <assert.h>
#include <fstream>

namespace Xyce {
namespace Analysis {
    
//-----------------------------------------------------------------------------
// Function      : ROL::ROL
// Purpose       :
//-----------------------------------------------------------------------------
ROL::ROL(
   AnalysisManager &analysis_manager, 
   Nonlinear::Manager &nonlinear_manager, // TT
   Loader::Loader &loader, 
   Linear::System & linear_system,
   Topo::Topology & topology,
   IO::InitialConditionsManager & initial_conditions_manager)
  : AnalysisBase(analysis_manager, "ROL"),
    solutionPtrVector_(0),
    statePtrVector_(0),
    constraintPtrVector_(0),
    jvecPtrVector_(0),
    testPtrVector_(0),
    mydfdpPtrVector_(0),
    mydqdpPtrVector_(0),
    mydbdpPtrVector_(0),
    mysensRHSPtrVector_(0),
    analysisManager_(analysis_manager),
    nonlinearManager_(nonlinear_manager), // TT
    loader_(loader),
    topology_(topology),
    initialConditionsManager_(initial_conditions_manager),
    linearSystem_(linear_system),
    outputManagerAdapter_(analysis_manager.getOutputManagerAdapter()),
    stepLoopSize_(0),
    rolLoopInitialized_(false),
    numParams_(0),
    numSensParams_(0)
{
  // For testing purposes get parameters from txt file
  std::ifstream param_file;
  std::string deviceName,paramName,dummy;
  RealT temp;
  numParams_ = 0;
  param_file.open("parameters.txt");
  getline(param_file,dummy); // skip the first line
  while (true)
  {
    if(!(param_file >> deviceName)) break;
    if(!(param_file >> paramName)) break;
    paramName = deviceName+":"+paramName;
    paramNameVec_.push_back(paramName);
    if(!(param_file >> temp)) break; // init guess
    if(!(param_file >> temp)) break; // lo bound
    if(!(param_file >> temp)) break; // up bound
    numParams_ += 1;
  }
  param_file.close();
}

//-----------------------------------------------------------------------------
// Function      : ROL::~ROL
// Purpose       :
//-----------------------------------------------------------------------------
ROL::~ROL(){}
  
//TT : copied from DCSweep;
//-----------------------------------------------------------------------------
// Function      : ROL::setTranOptions
// Purpose       :
// Special Notes : These are from '.options timeint'
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool ROL::setTimeIntegratorOptions(
   const Util::OptionBlock &     option_block)
{
  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; ++it)
  {
    const Util::Param &param = (*it);

    if (param.uTag() == "DAESTATEDERIV")
      analysisManager_.setDAEStateDerivFlag(static_cast<bool> (param.getImmutableValue<int>()));
    else if (param.uTag() == "DEBUGLEVEL")
      IO::setTimeIntegratorDebugLevel(analysisManager_.getCommandLine(), param.getImmutableValue<int>());
    else if (nonlinearManager_.setReturnCodeOption(param))
      ;
    else if (tiaParams_.setTimeIntegratorOption(param))
      ;
    else if (setDCOPOption(param))
      ;
    else if (param.uTag() == "METHOD")
      ;
    else
      Report::UserError() << param.uTag() << " is not a recognized time integration option";
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ROL::setAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 6/22/10
//-----------------------------------------------------------------------------
bool ROL::setAnalysisParams(const Util::OptionBlock & paramsBlock)
{
  stepSweepVector_.push_back(parseSweepParams(paramsBlock.begin(), paramsBlock.end()));
  outputManagerAdapter_.setStepSweepVector(stepSweepVector_);

  return true;
}

const TimeIntg::TIAParams &
ROL::getTIAParams() const
{
  return tiaParams_; // TT
}

TimeIntg::TIAParams &
ROL::getTIAParams()
{
  return tiaParams_; // TT
}

//-----------------------------------------------------------------------------
// Function      : ROL::run()
// Purpose       : This is the main controlling loop for ROL analysis.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/04/00
//-----------------------------------------------------------------------------
// Revised: Timur Takhtaganov, 06/03/2015
//         
bool ROL::doRun()
{
  // Instead of running DCSweep type analysis, we will run ROL analysis; the following will be only used in solve() function in the EqualityConstraint class: doInit() && doLoopProcess() && doFinish();
  return doInit() && runROLAnalysis() && doFinish();
  // return doInit() && doLoopProcess() && doFinish();
}

//-----------------------------------------------------------------------------
// Function      : ROL::init()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 06/02/2015
//-----------------------------------------------------------------------------
bool ROL::doInit()
{
  bool status;
  
  nonlinearManager_.enableSensitivity(
      *analysisManager_.getDataStore(), 
      analysisManager_.getStepErrorControl(),
      *analysisManager_.getPDSManager(), topology_, 
      outputManagerAdapter_.getOutputManager(),
      numSensParams_); // this only needs to be called once

  if ( !rolLoopInitialized_ )
  {
    if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
    {
      Xyce::dout() << std::endl << std::endl;
      Xyce::dout() << section_divider << std::endl;
      Xyce::dout() << "ROL Sweep::init" << std::endl;
    }

    // TT: processes sweep parameters; location N_ANP_SweepParam.C
    stepLoopSize_ = setupSweepLoop(analysisManager_.getComm(), loader_, stepSweepVector_.begin(), stepSweepVector_.end());

    // status = doAllocations(stepLoopSize_,numParams_); // TT: call here if runROLAnalysis is not called

    // TT: created ROLEvent, but this still doesn't work, probably needs be included in a bunch of places
    //Util::publish<ROLEvent>(analysisManager_, ROLEvent(ROLEvent::INITIALIZE, stepSweepVector_, stepLoopSize_));

    // TT: two lines from DCSweep, changed to Step instead of DC
    //outputManagerAdapter_.setStepAnalysisMaxSteps( stepLoopSize_ );
    outputManagerAdapter_.setStepSweepVector(stepSweepVector_);
    
    rolLoopInitialized_ = true;

    Util::publish<StepEvent>(analysisManager_, StepEvent(StepEvent::INITIALIZE, stepSweepVector_, stepLoopSize_)); // TT
    
  }

  // analysisManager_.setStepLoopInitialized(true);
  
  // TT: copied from DCSweep
  //setup for operating pt calculation
  baseIntegrationMethod_ = TimeIntg::methodsEnum::NO_TIME_INTEGRATION;
  analysisManager_.createTimeIntegratorMethod(tiaParams_, baseIntegrationMethod_);

  // analysisManager_.setAnalysisMode(ANP_MODE_ROL); // TT: sets analysis mode; added ANP_MODE_ROL to N_ANP_AnalysisManager.C and N_ANP_fwd.C 

  // TT: copied from DCSweep;
  stepNumber = 0;
  setDoubleDCOPEnabled(loader_.isPDESystem());
  if (getDoubleDCOPEnabled() && getDoubleDCOPStep() == 0)
  {
    nonlinearManager_.setAnalysisMode(nonlinearAnalysisMode(ANP_MODE_DC_NLPOISSON));
  }


  initializeSolution_(); // TT: this function is copied from DCSweep

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ROL::initializeSolution()
// Purpose       : Imitate the behavior of DCSweep
//                 
// Special Notes : Copied from DCSweep
// Scope         : public
// Creator       : Timur Takhtaganov
// Creation Date : 06/02/2015
//-----------------------------------------------------------------------------
void ROL::initializeSolution_()
{
  // set initial guess, if there is one to be set.
  // this setInitialGuess call is to up an initial guess in the
  // devices that have them (usually PDE devices).  This is different than
  // the "intializeProblem" call, which sets IC's.  (initial conditions are
  // different than initial guesses.
  loader_.setInitialGuess(analysisManager_.getDataStore()->nextSolutionPtr); //TT: this function calls a function by the same name in DeviceManager; sets initial guess for devices that have one (PDE?)

  // If available, set initial solution (.IC, .NODESET, etc).
  setInputOPFlag(
     initialConditionsManager_.setupInitialConditions(outputManagerAdapter_.getComm(),
                                                      topology_.getSolutionNodeNameMap(),
                                                      outputManagerAdapter_.getAliasNodeMap(),
                                                      *analysisManager_.getDataStore()->nextSolutionPtr,
                                                      linearSystem_));

  // Set a constant history for operating point calculation
  analysisManager_.getDataStore()->setConstantHistory();
  analysisManager_.getWorkingIntegrationMethod().obtainCorrectorDeriv();
}

//-----------------------------------------------------------------------------
// TT: this function gets called by the value() function in EqualityConstraint class; it updates sweep parameter so that correct source term is loaded into RHS vector
//-----------------------------------------------------------------------------
void ROL::setSweepValue(int step)
{
  bool reset = updateSweepParams(loader_, step, stepSweepVector_.begin(), stepSweepVector_.end(), false);// TT: update sweep parameters; location N_ANP_SweepParam.C
}

//-----------------------------------------------------------------------------
// Function      : ROL::loopProcess()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool ROL::doLoopProcess()
{
  bool integration_status = true;
  TimeIntg::DataStore & ds = *(analysisManager_.getDataStore()); // TT
  
  //Util::publish<StepEvent>(analysisManager_, StepEvent(StepEvent::INITIALIZE, stepSweepVector_, stepLoopSize_)); // TT: moved here from doInit();
  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::INITIALIZE, AnalysisEvent::DC));
  
  // StepEvent step_event(StepEvent::STEP_STARTED, stepSweepVector_, currentStep);

  int currentStep = 0;
  int finalStep = stepLoopSize_;
  while (currentStep < finalStep)
  {
    outputManagerAdapter_.setDCAnalysisStepNumber(currentStep);
      
    // Tell the manager if any of our sweeps are being reset in this loop iteration.
    bool reset = updateSweepParams(loader_, currentStep, stepSweepVector_.begin(), stepSweepVector_.end(), false);// TT: update sweep parameters; location N_ANP_SweepParam.C
      
    analysisManager_.setSweepSourceResetFlag(reset);// TT: sets sweep source flag to reset in AnalysisManager
          
    outputManagerAdapter_.setStepSweepVector(stepSweepVector_);
      
    if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
      for (SweepVector::const_iterator it = stepSweepVector_.begin(), end = stepSweepVector_.end(); it != end; ++it)
        Xyce::dout() << "ROL DC Sweep # " << currentStep <<"\t" << (*it);
      
    // TT
    if (currentStep != 0 && reset)
    {
      analysisManager_.getDataStore()->setZeroHistory();
      initializeSolution_();
    }
    
          
    // TT: the following line causes currSol and nextSol reset to zero
    // Util::publish<StepEvent>(analysisManager_, StepEvent(StepEvent::STEP_STARTED, stepSweepVector_, currentStep)); 

    static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_STARTED, AnalysisEvent::DC, 0.0, currentStep));
    //step_event.state = StepEvent::STEP_STARTED;
    //Util::publish<StepEvent>(analysisManager_, step_event);
  
    takeStep_(); // TT: take integration step

    // Set things up for the next time step, based on if this one was
    // successful.
    if (analysisManager_.getStepErrorControl().stepAttemptStatus)
    {
      static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_SUCCESSFUL, AnalysisEvent::DC, 0.0, currentStep));
      doProcessSuccessfulStep();
      // collect solutions here
      *(solutionPtrVector_[currentStep]) = *(ds.currSolutionPtr); 
      // solutionPtrVector_[currentStep]->print(Xyce::dout());
    }
    else // stepAttemptStatus  (ie do this if the step FAILED)
    {
      static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_FAILED, AnalysisEvent::DC, 0.0, currentStep));
      doProcessFailedStep();
    }
      
    currentStep = stepNumber;
  } // end of sweep loop

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::FINISH, AnalysisEvent::DC));

  // Util::publish<StepEvent>(analysisManager_, StepEvent(StepEvent::STEP_COMPLETED, stepSweepVector_, currentStep));
  
  // step_event.state_ = StepEvent::STEP_COMPLETED;
  // step_event.finalSimTime_ = getTIAParams().finalTime;
  // Util::publish<StepEvent>(analysisManager_, step_event);

  return integration_status;
}

//-----------------------------------------------------------------------------
// TT: copied from DCSweep
//-----------------------------------------------------------------------------
bool ROL::doHandlePredictor()
{
  analysisManager_.getDataStore()->setErrorWtVector(tiaParams_, topology_.getVarTypes());
  analysisManager_.getWorkingIntegrationMethod().obtainPredictor();
  analysisManager_.getWorkingIntegrationMethod().obtainPredictorDeriv();
  
  // In case this is the upper level of a 2-level sim, tell the
  // inner solve to do its prediction:
  bool        beginIntegrationFlag = analysisManager_.getBeginningIntegrationFlag();           // system_state.beginIntegrationFlag;
  double      nextTimeStep = analysisManager_.getStepErrorControl().currentTimeStep;           // system_state.nextTimeStep;
  double      nextTime = analysisManager_.getStepErrorControl().nextTime;                      // system_state.nextTime;
  int         currentOrder = analysisManager_.getWorkingIntegrationMethod().getOrder();        // system_state.currentOrder;

  loader_.startTimeStep(beginIntegrationFlag, nextTimeStep, nextTime, currentOrder);

  return true;
}

//-----------------------------------------------------------------------------
// TT: copied from DCSweep
//-----------------------------------------------------------------------------
void ROL::takeStep_()
{
  { // Integration step predictor
    Stats::StatTop _predictorStat("Predictor");
    Stats::TimeBlock _predictorTimer(_predictorStat);

    doHandlePredictor();
  }
  
  { // Load B/V source devices with time data
    Stats::StatTop _updateDeviceSourceStat("Update Device Sources");
    Stats::TimeBlock _updateDeviceSourceTimer(_updateDeviceSourceStat);

    loader_.updateSources();
  }
  
  { // Nonlinear solve
    Stats::StatTop _nonlinearSolveStat("Solve");
    Stats::TimeBlock _nonlinearSolveTimer(_nonlinearSolveStat);

    analysisManager_.getStepErrorControl().newtonConvergenceStatus = nonlinearManager_.solve();
  }

  { // Add change to solution
    Stats::StatTop _errorStat("Error Estimation");
    Stats::TimeBlock _errorTimer(_errorStat);

    analysisManager_.getWorkingIntegrationMethod().stepLinearCombo();

    gatherStepStatistics(stats_, nonlinearManager_.getNonlinearSolver(), analysisManager_.getStepErrorControl().newtonConvergenceStatus);

    analysisManager_.getStepErrorControl().evaluateStepError(loader_, tiaParams_);
  }
}
  
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
bool ROL::doAllocations(int nc, int nz)
{
  // Allocate space for solution and simulation space vectors
  solutionPtrVector_.resize(nc);
  statePtrVector_.resize(nc);
  constraintPtrVector_.resize(nc);
  jvecPtrVector_.resize(nc);
  testPtrVector_.resize(nc);

  for (int i=0;i<nc;i++)
  {
    solutionPtrVector_[i]   = linearSystem_.builder().createVector();
    statePtrVector_[i]      = linearSystem_.builder().createVector();
    constraintPtrVector_[i] = linearSystem_.builder().createVector();
    jvecPtrVector_[i]       = linearSystem_.builder().createVector();
    testPtrVector_[i]       = linearSystem_.builder().createVector();
  }
  // Allocate space for sensitivity vectors
  mydfdpPtrVector_.resize(nz);
  mydqdpPtrVector_.resize(nz);
  mydbdpPtrVector_.resize(nz);
  mysensRHSPtrVector_.resize(nz);
  for (int i=0;i<nz;i++)
  {
    mydfdpPtrVector_[i] = linearSystem_.builder().createVector();
    mydqdpPtrVector_[i] = linearSystem_.builder().createVector();
    mydbdpPtrVector_[i] = linearSystem_.builder().createVector();
    mysensRHSPtrVector_[i] = linearSystem_.builder().createVector();
  }

  return true;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
bool ROL::doFree()
{
  paramNameVec_.clear();
  for (int i=0;i<stepLoopSize_;i++)
  {
    delete solutionPtrVector_[i];
    solutionPtrVector_[i] = 0;
    delete statePtrVector_[i];
    statePtrVector_[i] = 0;
    delete constraintPtrVector_[i];
    constraintPtrVector_[i] = 0;
    delete jvecPtrVector_[i];
    jvecPtrVector_[i] = 0;
    delete testPtrVector_[i];
    testPtrVector_[i] = 0;
  }
  solutionPtrVector_.clear();
  statePtrVector_.clear();
  constraintPtrVector_.clear();
  jvecPtrVector_.clear();
  testPtrVector_.clear();

  for (int i=0;i<numParams_;i++)
  {
    delete mydfdpPtrVector_[i];
    mydfdpPtrVector_[i] = 0;
    delete mydqdpPtrVector_[i];
    mydqdpPtrVector_[i] = 0;
    delete mydbdpPtrVector_[i];
    mydbdpPtrVector_[i] = 0;
    delete mysensRHSPtrVector_[i];
    mysensRHSPtrVector_[i] = 0;
  }
  mydfdpPtrVector_.clear();
  mydqdpPtrVector_.clear();
  mydbdpPtrVector_.clear();
  mysensRHSPtrVector_.clear();

  return true;
}

//-----------------------------------------------------------------------------
// TT: copied from DCSweep
//-----------------------------------------------------------------------------
bool ROL::doProcessSuccessfulStep()
{
  TimeIntg::DataStore & ds = *(analysisManager_.getDataStore());

  Stats::StatTop _processSuccessfulStepStat("Successful Step");
  Stats::TimeBlock _processSuccessfulStepTimer(_processSuccessfulStepStat);

  loader_.stepSuccess(TWO_LEVEL_MODE_DC_SWEEP); // analysisManager_.getTwoLevelMode());

  // This output call is for device-specific output, such as .OP,
  // or internal plot output from PDE(TCAD) devices.
  loader_.outputPlotFiles();

  // TT : not using
  // if (sensFlag_ && !firstDoubleDCOPStep() )
  // {
  //   nonlinearManager_.calcSensitivity(objectiveVec_, dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);
  // }

  // Do some statistics, as long as this isn't the first "double"
  // DCOP step. (that one doesn't count)
  if ( !firstDoubleDCOPStep() )
  {
    stepNumber += 1; // TT: this is inherited from AnalysisBase
    stats_.successStepsThisParameter_ += 1;
    stats_.successfulStepsTaken_ += 1;
      
  }
    
  // update the data arrays, output:
  analysisManager_.getDataStore()->updateSolDataArrays();

  //TT: print out current solution and rhs
  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    if ( !firstDoubleDCOPStep() )
    {
      std::cout << "Current solution" << std::endl;
      for (unsigned int k = 0; k < ds.solutionSize; k++)
      {
        std::cout << (*(ds.currSolutionPtr))[k] << std::endl;
      }

      for (int i=0;i<3;i++)
      {
        std::cout << (*(ds.dFdxdVpVectorPtr))[i] << std::endl;
      }

      std::cout << "Current f vector" << std::endl;
      for (unsigned int k = 0; k < ds.daeFVectorPtr->globalLength(); k++)
      {
        std::cout << (*ds.daeFVectorPtr)[k] << std::endl;
      }
      std::cout << "Current q vector" << std::endl;
      for (unsigned int k = 0; k < ds.daeQVectorPtr->globalLength(); k++)
      {
        std::cout << (*ds.daeQVectorPtr)[k] << std::endl;
      }
      std::cout << "Current b vector" << std::endl;
      for (unsigned int k = 0; k < ds.daeBVectorPtr->globalLength(); k++)
      {
        std::cout << (*ds.daeBVectorPtr)[k] << std::endl;
      }
    }
  }
    
  //dcSweepOutput(); // TT: need something else here
  
  // now that output has been called, update the doubleDCOP step
  // if neccessary. (pde-only)
  nextDCOPStep();
  // std::cout << "Setting analysis mode to DC_SWEEP" << std::endl;
  nonlinearManager_.setAnalysisMode(nonlinearAnalysisMode(ANP_MODE_DC_SWEEP));

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ROL::processFailedStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool ROL::doProcessFailedStep()
{
  Stats::StatTop _processFailedStat("Failed Steps");
  Stats::TimeBlock _processFailedTimer(_processFailedStat);

  loader_.stepFailure(analysisManager_.getTwoLevelMode());

  stepNumber += 1;
  rolSweepFailures_.push_back(stepNumber);
  stats_.failedStepsAttempted_  += 1;
  analysisManager_.getStepErrorControl().numberSuccessiveFailures += 1;

  return true;
}


//-----------------------------------------------------------------------------
// Function      : ROL::doFinish()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool ROL::doFinish()
{
  Util::publish<StepEvent>(analysisManager_, StepEvent(StepEvent::FINISH, stepSweepVector_, stepLoopSize_));

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ROL::twoLevelStep
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
// TT: moved from DCSweep; I think its needed for pde devices
bool ROL::twoLevelStep()
{
  loader_.updateSources();
  analysisManager_.getStepErrorControl().newtonConvergenceStatus = nonlinearManager_.solve();
  analysisManager_.getDataStore()->stepLinearCombo();
  gatherStepStatistics(stats_, nonlinearManager_.getNonlinearSolver(), analysisManager_.getStepErrorControl().newtonConvergenceStatus);
  analysisManager_.getStepErrorControl().evaluateStepError(loader_, tiaParams_);

  return analysisManager_.getStepErrorControl().stepAttemptStatus;
}

#if Xyce_ROL
typedef double RealT;
//-----------------------------------------------------------------------------
// Function      : ROL::runROLAnalysis
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Timur Takhtaganov
// Creation Date : 06/02/2015
//-----------------------------------------------------------------------------
bool ROL::runROLAnalysis()
{

  // Util::publish<StepEvent>(analysisManager_, StepEvent(StepEvent::STEP_STARTED, stepSweepVector_, currentStep)); 
  
#if UQ==1
  //Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
#endif
  
  std::string msg;
  bool status = true;
  int errorFlag = 0;

  try
  {    

    int nu = analysisManager_.getDataStore()->nextSolutionPtr->globalLength(); // number of solution variables
    int nc = stepLoopSize_; // number of constraint equations
    int nz = numParams_; // number of optimization parameters

    status = doAllocations(nc,nz); // allocate solution and sensitivity arrays

    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, Teuchos::Ptr<Teuchos::ParameterList>(&*parlist) );

    std::string outputname     = parlist->get("Output File","rol_output.txt");
    std::ofstream out(outputname.c_str());
    Teuchos::RCP<std::ostream> outStream = Teuchos::rcp(&out,false);

    out << "nu = " << nu << " nc = " << nc << " nz = " << nz << std::endl;    
    out << "Entering ROL loop" << std::endl;

    bool do_checks      = parlist->get("Do Checks",true);
    bool use_scale      = parlist->get("Use Scaling For Epsilon-Active Sets",true);
    bool use_sqp        = parlist->get("Use SQP", true);
    bool use_lsearch    = parlist->get("Use Line Search", true);    
    bool use_TR         = parlist->get("Use Trust Region", true);
    bool use_bundle     = parlist->get("Use Bundle Method", true);
    bool use_bcon       = parlist->get("Use Bound Constraints", true);
    RealT alpha         = parlist->get("Penalty Parameter", 1.e-4);
    // amplifier problem parameters
    RealT ampl          = parlist->get("Amplifier Gain", 4.0);
    int ptype           = parlist->get("Penalty Type", 1);

#if UQ==1
    int nSamp           = parlist->get("Number of Samples (UQ)", 10);
    int samplertype     = parlist->get("Sampler Type", 1); // 1 - MC, 2 - QMC
    RealT gamma         = parlist->get("CVaR: gamma", 1.e-4);
    RealT prob          = parlist->get("CVaR: prob", 0.99);
    RealT coeff         = parlist->get("CVaR: coeff", 1.0);
#endif

    //************ Equality constraint ***************//
    //Teuchos::RCP< ::ROL::EqualityConstraint_SimOpt<RealT> > con;
    Teuchos::RCP< ::ROL::ParametrizedEqualityConstraint_SimOpt<RealT> > con;
    con = Teuchos::rcp(new EqualityConstraint_ROL_DC<RealT>(nu,nc,nz,analysisManager_,analysisManager_.getNonlinearEquationLoader(),nonlinearManager_.getNonlinearSolver(),linearSystem_,*this));
    
    //************ Optimization-space vectors ***************//
    Teuchos::RCP<std::vector<RealT> > z_rcp     = Teuchos::rcp( new std::vector<RealT> (nz, 0.0) );
    Teuchos::RCP<std::vector<RealT> > z1_rcp     = Teuchos::rcp( new std::vector<RealT> (nz, 0.0) );
    Teuchos::RCP<std::vector<RealT> > z2_rcp     = Teuchos::rcp( new std::vector<RealT> (nz, 0.0) );
    Teuchos::RCP<std::vector<RealT> > z3_rcp     = Teuchos::rcp( new std::vector<RealT> (nz, 0.0) );
    Teuchos::RCP<std::vector<RealT> > yz_rcp    = Teuchos::rcp( new std::vector<RealT> (nz, 0.0) );
    Teuchos::RCP<std::vector<RealT> > ajvz_rcp  = Teuchos::rcp( new std::vector<RealT> (nz, 0.0) );
    ::ROL::StdVector<RealT> z(z_rcp);
    ::ROL::StdVector<RealT> z1(z1_rcp);
    ::ROL::StdVector<RealT> z2(z2_rcp);
    ::ROL::StdVector<RealT> z3(z3_rcp);
    ::ROL::StdVector<RealT> yz(yz_rcp);
    ::ROL::StdVector<RealT> ajvz(ajvz_rcp);
    
    // Simulation and constraint space vectors

    // Using default constructor
    // Linear::ROL_XyceVector<RealT> u(nc,*analysisManager_.getDataStore()->nextSolutionPtr);
    // Linear::ROL_XyceVector<RealT> yu(nc,*analysisManager_.getDataStore()->nextSolutionPtr);
    // Linear::ROL_XyceVector<RealT> c(nc,*analysisManager_.getDataStore()->nextSolutionPtr);
    // Linear::ROL_XyceVector<RealT> jv(nc,*analysisManager_.getDataStore()->nextSolutionPtr);

    // Using copy constructor
    Linear::ROL_XyceVector<RealT> u(statePtrVector_);
    Linear::ROL_XyceVector<RealT> yu(testPtrVector_);
    Linear::ROL_XyceVector<RealT> c(constraintPtrVector_);
    Linear::ROL_XyceVector<RealT> jv(jvecPtrVector_);

    u.randomize();
    yu.randomize();
    c.randomize();
    jv.randomize();

    // // Testing cloning
    // Teuchos::RCP< ::ROL::Vector<RealT> > copy = u.clone();
    // // assert (!is_null(copy));
    // Linear::ROL_XyceVector<RealT> cop = dynamic_cast< Linear::ROL_XyceVector<RealT> &>(*copy);
    // cop.randomize();
    // cop.print(*outStream);

    // *outStream << "Checking linear algebra" << std::endl;
    // std::vector<RealT> consistency = jv.checkVector(u,c,true,*outStream);
    
    // Get Initial Guess and bounds from file
    std::ifstream param_file;
    param_file.open("parameters.txt");
    RealT temp;
    std::string parName, dummy;
    std::vector<RealT> up_bounds(nz,0.0);
    std::vector<RealT> lo_bounds(nz,0.0);
    getline(param_file,dummy); // skip the first line
    for(int i=0;i<nz;i++)
    {
      param_file >> parName;
      param_file >> parName;
      param_file >> temp; // init guess
      (*z_rcp)[i] = temp;
      (*yz_rcp)[i] = temp;
      param_file >> temp; // lo bound
      lo_bounds[i] = temp;
      param_file >> temp; // up bound
      up_bounds[i] = temp;
    }
    param_file.close();
    
#if UQ==1
    // Get uncertain parameters and their bounds
    param_file.open("uncertain_parameters.txt");
    std::vector<std::vector<RealT> > ubounds;
    std::vector<RealT> tempv(2,0);
    std::string deviceName,paramName;
    getline(param_file,dummy);
    while(true)
    {
      if(!(param_file >> deviceName))
        break;
      if(deviceName[0]=='*')
      {
        getline(param_file,dummy);
        continue;
      }
      if(!(param_file >> paramName))
        break;
      paramName = deviceName+":"+paramName;
      std::cout << "Uncertain parameter " << paramName << std::endl;
      uncertainParams_.push_back(paramName);
      if(!(param_file >> temp)) break;
      tempv[0] = temp;
      if(!(param_file >> temp)) break;
      tempv[1] = temp;
      ubounds.push_back(tempv);
    }
    param_file.close();
#endif

    RealT tol = 1.e-12;

    Teuchos::RCP< ::ROL::Vector<RealT> > up = Teuchos::rcp(&u,false);
    Teuchos::RCP< ::ROL::Vector<RealT> > yup = Teuchos::rcp(&yu,false);
    Teuchos::RCP< ::ROL::Vector<RealT> > cp = Teuchos::rcp(&c,false);
    Teuchos::RCP< ::ROL::Vector<RealT> > jvp = Teuchos::rcp(&jv,false);

    Teuchos::RCP< ::ROL::Vector<RealT> > zp = Teuchos::rcp(&z,false);
    Teuchos::RCP< ::ROL::Vector<RealT> > z3p = Teuchos::rcp(&z3,false);
    Teuchos::RCP< ::ROL::Vector<RealT> > yzp = Teuchos::rcp(&yz,false); 
    Teuchos::RCP< ::ROL::Vector<RealT> > ajvzp = Teuchos::rcp(&ajvz,false); 

    // SimOpt vectors
    ::ROL::Vector_SimOpt<RealT> x(up,zp);
    ::ROL::Vector_SimOpt<RealT> y(yup,yzp);
    ::ROL::Vector_SimOpt<RealT> ajv(jvp,ajvzp);

#if UQ==1
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build samplers
    Teuchos::RCP< ::ROL::BatchManager<double> > bman = Teuchos::rcp(new ::ROL::StdTeuchosBatchManager<double,int>(comm));
    Teuchos::RCP< ::ROL::SampleGenerator<double> > sampler;
    std::ifstream dirNums("new-joe-kuo-6.21201");
    if(samplertype==1)
      sampler = Teuchos::rcp(new ::ROL::MonteCarloGenerator<double>(50*nSamp,ubounds,bman,false,false,100));
    else if(samplertype==2)
    {
      if(!dirNums)
        out << "DirNums file not found!" << std::endl;
      int skip = 0;
      sampler = Teuchos::rcp(new ::ROL::QuasiMonteCarloGenerator<double>(nSamp,ubounds,bman,dirNums,skip,false,false,100));
    }
    con->setParameter(sampler->getMyPoint(0));
#endif
   
    con->solve(u,z,tol);
#if UQ==1
    std::ofstream allsamples("all_samples.txt");
    std::ofstream samples("good_samples.txt");
    int nGSamp = 0; // number of "good" samples
    out << "Generating good samples..." << std::endl;
    for (int k=0;k<50*nSamp;k++)
    {
      con->setParameter(sampler->getMyPoint(k));
#endif
      //======================= Separate testing of Jacobian_1 and 2 ==============================//
      
      out << "Checking Jacobian_1" << std::endl;
      int numSteps = 13;
      std::vector<RealT> steps(numSteps);
      for(int i=0;i<numSteps;++i)
      {
        steps[i] = pow(10,-i);
      }
      
      using ::ROL::Finite_Difference_Arrays::shifts;
      using ::ROL::Finite_Difference_Arrays::weights;
      
      int numVals = 4;
      std::vector<RealT> tmp(numVals);
      std::vector<std::vector<RealT> > jvCheck(numSteps, tmp);
      
      // Compute constraint value at x.
      Teuchos::RCP< ::ROL::Vector<RealT> > cc = c.clone();
      con->value(*cc,*(x.get_1()),*(x.get_2()),tol);
      
      // Compute (Jacobian at x) times (vector v).
      Teuchos::RCP< ::ROL::Vector<RealT> > Jv = jv.clone();
      RealT normJv;
      con->applyJacobian_1(*Jv, *(y.get_1()), *(x.get_1()), *(x.get_2()), tol);
      normJv = Jv->norm();
      
      // Temporary vectors.
      Teuchos::RCP< ::ROL::Vector<RealT> > cdif = c.clone();
      Teuchos::RCP< ::ROL::Vector<RealT> > cnew = c.clone();
      Teuchos::RCP< ::ROL::Vector<RealT> > xnew = x.clone();
      
      // FD order
      int order = 1;
      RealT minGradError = 1.e+10; // for "good" samples
      
      for (int i=0; i<numSteps; i++)
      {
        RealT eta = steps[i];
        
        xnew->set(x);
        
        cdif->set(*cc);
        cdif->scale(weights[order-1][0]);
        
        for(int j=0; j<order; ++j)
              {
          
          xnew->axpy(eta*shifts[order-1][j], y);
          
          if( weights[order-1][j+1] != 0 )
                {
            //this->update(*xnew);
            ::ROL::Vector_SimOpt<RealT> &xnews = Teuchos::dyn_cast< ::ROL::Vector_SimOpt<RealT> >(*xnew);
            con->value(*cnew,*(xnews.get_1()),z,tol);
            cdif->axpy(weights[order-1][j+1],*cnew);    
          }
          
        }
        
        cdif->scale(1.0/eta);    
        
        // Compute norms of Jacobian-vector products, finite-difference approximations, and error.
        jvCheck[i][0] = eta;
        jvCheck[i][1] = normJv;
        jvCheck[i][2] = cdif->norm();
        cdif->axpy(-1.0, *Jv);
        jvCheck[i][3] = cdif->norm();
        
        //if (printToStream)
              //{
        std::stringstream hist;
        if (i==0)
        {
          hist << std::right
               << std::setw(20) << "Step size"
               << std::setw(20) << "norm(Jac*vec)"
               << std::setw(20) << "norm(FD approx)"
               << std::setw(20) << "norm(abs error)"
               << "\n"
               << std::setw(20) << "---------"
               << std::setw(20) << "-------------"
               << std::setw(20) << "---------------"
               << std::setw(20) << "---------------"
               << "\n";
        }
        hist << std::scientific << std::setprecision(11) << std::right
             << std::setw(20) << jvCheck[i][0]
             << std::setw(20) << jvCheck[i][1]
             << std::setw(20) << jvCheck[i][2]
             << std::setw(20) << jvCheck[i][3]
             << "\n";
        out << hist.str();
        
        minGradError = std::min(jvCheck[i][3],minGradError);
      }
      out << "\n" << std::endl;
      
#if UQ==1
      std::cout << "minGradError = " << minGradError << std::endl;
      for(int j=0;j<ubounds.size();j++)
        allsamples << std::scientific << std::setprecision(8) << (sampler->getMyPoint(k))[j] << ' ';
      allsamples << '\n';
      if( minGradError < 2.e-6 )
      {
        for(int j=0;j<ubounds.size();j++)
          samples << std::scientific << std::setprecision(8) << (sampler->getMyPoint(k))[j] << ' ';
        samples << '\n';
        nGSamp++;
      }
#endif
      
      out << "Checking Jacobian_2" << std::endl;
        
      // Compute constraint value at x.
      con->value(*cc,*(x.get_1()),*(x.get_2()),tol);
      
      // Compute (Jacobian at x) times (vector v).
      con->applyJacobian_2(*Jv, *(y.get_2()), *(x.get_1()), *(x.get_2()), tol);
      normJv = Jv->norm();
      
      for (int i=0; i<numSteps; i++)
      {
        RealT eta = steps[i];

        xnew->set(x);

        cdif->set(*cc);
        cdif->scale(weights[order-1][0]);

        for(int j=0; j<order; ++j)
        {
          xnew->axpy(eta*shifts[order-1][j], y);
          
          if( weights[order-1][j+1] != 0 )
          {
            //this->update(*xnew);
            ::ROL::Vector_SimOpt<RealT> &xnews = Teuchos::dyn_cast< ::ROL::Vector_SimOpt<RealT> >(*xnew);
            con->value(*cnew,u,*(xnews.get_2()),tol);
            cdif->axpy(weights[order-1][j+1],*cnew);    
          }
        }
        cdif->scale(1.0/eta);    
          
        // Compute norms of Jacobian-vector products, finite-difference approximations, and error.
        jvCheck[i][0] = eta;
        jvCheck[i][1] = normJv;
        jvCheck[i][2] = cdif->norm();
        cdif->axpy(-1.0, *Jv);
        jvCheck[i][3] = cdif->norm();

        //if (printToStream)
              //{
        std::stringstream hist;
        if (i==0)
        {
          hist << std::right
               << std::setw(20) << "Step size"
               << std::setw(20) << "norm(Jac*vec)"
               << std::setw(20) << "norm(FD approx)"
               << std::setw(20) << "norm(abs error)"
               << "\n"
               << std::setw(20) << "---------"
               << std::setw(20) << "-------------"
               << std::setw(20) << "---------------"
               << std::setw(20) << "---------------"
               << "\n";
        }
        hist << std::scientific << std::setprecision(11) << std::right
             << std::setw(20) << jvCheck[i][0]
             << std::setw(20) << jvCheck[i][1]
             << std::setw(20) << jvCheck[i][2]
             << std::setw(20) << jvCheck[i][3]
             << "\n";
        out << hist.str();
      }
      out << "\n" << std::endl;
            
            
      //===================== end of checking Jacobian ========================//
#if UQ==1
    }// end (loop over random params)
    out << "Generated " << nGSamp << " good samples." << std::endl;
    samples.close();
    allsamples.close();
#endif
    
    if(do_checks)
    {
      out << "Checking full Jacobian" << std::endl;
      con->checkApplyJacobian(x,y,jv,true,*outStream);
      
      // // TT: the following requires implementing basis for ROL_XyceVector
      // std::cout << "Checking full adjoint Jacobian" << std::endl;
      // con->applyAdjointJacobian_1(const_cast< ::ROL::Vector<RealT> &>(*(ajv.get_1())),yu,u,z,tol);
      // std::cout << "norm ajvu = " << ajv.get_1()->norm() << std::endl;
      // con->applyAdjointJacobian_2(const_cast< ::ROL::Vector<RealT> &>(*(ajv.get_2())),yu,u,z,tol);
      // std::cout << "norm ajvz = " << ajv.get_2()->norm() << std::endl;
      // std::cout << "norm ajv = " << std::sqrt(ajv.get_1()->norm()*ajv.get_1()->norm() + ajv.get_2()->norm()*ajv.get_2()->norm()) << std::endl;
      // con->checkApplyAdjointJacobian(x,yu,c,ajv);
      
      out << "Checking Jacobian consistency" << std::endl;
      con->checkAdjointConsistencyJacobian_1(jv,yu,u,z,true,*outStream);
      con->checkAdjointConsistencyJacobian_2(jv,yz,u,z,true,*outStream);
      
      out << "Checking consistency of solves" << std::endl;
      con->checkSolve(u,z,c,true,*outStream);
      con->checkInverseJacobian_1(jv,yu,u,z,true,*outStream);
      con->checkInverseAdjointJacobian_1(yu,jv,u,z,true,*outStream);
      
    } // end(do_checks)
    
    //*********************** Initialize objective and penalty functions ***************************//
    Teuchos::RCP< ::ROL::ParametrizedObjective_SimOpt<RealT> > obj;
    Teuchos::RCP< ::ROL::ParametrizedObjective_SimOpt<RealT> > pen;
    
#if OBJTYPE==0
    // L2-norm objective
    obj = Teuchos::rcp(new Objective_DC_L2Norm<RealT>(1.e-4,nc,nz));
#elif OBJTYPE==1
    // AMPlifier circuit objective
    std::cout << "Amplifier circuit problem" << std::endl;
    obj = Teuchos::rcp(new Objective_DC_AMP<RealT>(nc,nz));
    pen = Teuchos::rcp(new Penalty_DC_AMP<RealT>(ptype,alpha,ampl,nc,nz));
#endif
    std::vector< Teuchos::RCP< ::ROL::Objective<RealT> > > objVec;
    std::vector<bool> types(1,true);
    
    con->checkSolve(u,z,c,true,*outStream);
    if(do_checks)
    {
      out << "Checking objective" << std::endl;
      obj->checkGradient(x,x,y,true,*outStream);
      obj->checkHessVec(x,x,y,true,*outStream);
#if OBJTYPE==1
      out << "Checking penalty" << std::endl;
      pen->checkGradient(x,x,y,true,*outStream);
      pen->checkHessVec(x,x,y,true,*outStream);
#endif
    }
    
    //********************* Initialize reduced objective and penalty functions **********************//
    Teuchos::RCP< ::ROL::ParametrizedObjective<RealT> > robj = Teuchos::rcp(new ::ROL::Reduced_ParametrizedObjective_SimOpt<RealT>(obj,con,up,cp));
#if OBJTYPE==1
    Teuchos::RCP< ::ROL::ParametrizedObjective<RealT> > rpen = Teuchos::rcp(new ::ROL::Reduced_ParametrizedObjective_SimOpt<RealT>(pen,con,up,cp));
#endif

    if(do_checks)
    {
      *outStream << "Derivatives of reduced objective" << std::endl; 
      robj->checkGradient(z,z,yz,true,*outStream);
      robj->checkHessVec(z,z,yz,true,*outStream);
#if OBJTYPE==1
      *outStream << "Derivatives of reduced penalty" << std::endl; 
      rpen->checkGradient(z,z,yz,true,*outStream);
      rpen->checkHessVec(z,z,yz,true,*outStream);
#endif
    }

    //****************************** Bound constraints ****************************************//
    Teuchos::RCP<std::vector< RealT > > g0_rcp = Teuchos::rcp( new std::vector<RealT> (nz, 0.0) );
    ::ROL::StdVector<RealT> g0p(g0_rcp);
    robj->gradient(g0p,z,tol);
    *outStream << std::scientific << "Norm of initial gradient = " << g0p.norm() << "\n";
    // Define scaling for epsilon-active sets (used in inequality constraints)
    RealT scale;
    if(use_scale) 
      scale = 1.0e-2/g0p.norm();
    else
      scale = 1.0;
    *outStream << std::scientific << "Scaling: " << scale << "\n";
    Teuchos::RCP<::ROL::BoundConstraint<RealT> > bcon = Teuchos::rcp(new BoundConstraint_ROL_DC<RealT>(scale,lo_bounds,up_bounds));
    if(!use_bcon)
      bcon->deactivate();
    
    /**********************************************************************************************/    
    //***************************** Optimization ************************************************//
    /**********************************************************************************************/
    *outStream << "\n Initial guess " << std::endl;
    for (int i=0;i<nz;i++)
      *outStream << paramNameVec_[i] << " = " << (*z_rcp)[i] << std::endl;
    RealT gtol = parlist->get("Gradient Tolerance",1.e-14); // norm of gradient tolerance
    RealT stol = parlist->get("Step Tolerance",1.e-16); // norm of step tolerance
    int maxit  = parlist->get("Maximum Number of Iterations",100);// maximum number of iterations

    con->solve(u,z,tol); // start feasible

#if UQ==1
    //*************** Initial objective values and initial states (UQ) ***************//
    std::ofstream sampleWeights("weights.txt");
    for(int i=0;i<nSamp;i++)
    {
      sampleWeights << 1./(RealT)nSamp << '\n';
    }
    sampleWeights.close();
    std::string input_samples = "good_samples.txt";
    std::string input_weights = "weights.txt";
    // This sampler used to solve problem (with nSamp samples)
    Teuchos::RCP< ::ROL::SampleGenerator<double> > samplerSolve = Teuchos::rcp(new ::ROL::UserInputGenerator<RealT>(input_samples,input_weights,nSamp,ubounds.size(),bman));
    std::vector<RealT> sampleMean(ubounds.size(),0);
    // This sampler used to test solution (with nGSamp samples)
    Teuchos::RCP< ::ROL::SampleGenerator<double> > samplerTest  = Teuchos::rcp(new ::ROL::UserInputGenerator<RealT>(input_samples,input_weights,nGSamp,ubounds.size(),bman));
    std::ofstream data_IG_UQ("data_IG_UQ.txt");
    std::ofstream plot_IG_UQ("plot_IG_UQ.txt");
    RealT sumlin = 0, sumlin2 = 0, sumgain = 0, sumgain2 = 0, gain, oval, param;
    Teuchos::RCP<std::vector<Teuchos::RCP<Linear::Vector> > > uup = Teuchos::rcp_const_cast<std::vector<Teuchos::RCP<Linear::Vector> > >((Teuchos::dyn_cast<Linear::ROL_XyceVector<RealT> >(u)).getVector());
    for(int i=0;i<nSamp;i++)
    {
      con->setParameter(samplerSolve->getMyPoint(i));
      con->solve(u,z,tol);
      oval = obj.value(u,z,tol);
      for(int j=0;j<nc;j++)
      {
        plot_IG_UQ << std::scientific << std::setprecision(8) << (*(*uup)[j])[3] << ' ';
      }
      plot_IG_UQ << '\n';
      gain = ((*(*uup)[nc-1])[3]-(*(*uup)[0])[3])/2.0;
      sumgain += gain;
      sumgain2 += gain*gain;
      for(int j=0;j<uncertainParams_.size();j++)
      {
        param = (samplerSolve->getMyPoint(i))[j];
        data_IG_UQ << std::scientific << std::setprecision(8) << param << ' ';
        sampleMean[j] += param/nSamp;
      }
      data_IG_UQ << oval << ' ' << gain << '\n';
      sumlin += oval;
      sumlin2 += oval*oval;
    }
    plot_IG_UQ.close();
    RealT meanlin = sumlin/nSamp;
    RealT meangain = sumgain/nSamp;
    RealT varlin = (sumlin2 - nSamp*meanlin*meanlin)/(RealT)(nSamp-1);
    RealT vargain = (sumgain2 - nSamp*meangain*meangain)/(RealT)(nSamp-1);
    data_IG_UQ << '\n' << "------------------------------------------" << '\n';
    data_IG_UQ << "mean(lin) = " << meanlin << " " << '\n';
    data_IG_UQ << "sqrt(sample variance)(lin) = " << std::sqrt(varlin) << '\n';
    data_IG_UQ << "mean(gain) = " << meangain << " " << '\n';
    data_IG_UQ << "sqrt(sample variance)(gain) = " << std::sqrt(vargain) << '\n';
    data_IG_UQ.close();

    /**********************************************************************************************/
    //************************* Build risk-averse objective function *****************************//
    /**********************************************************************************************/
    Teuchos::RCP< ::ROL::Distribution<double> > dist;
    Teuchos::RCP< ::ROL::PlusFunction<double> > pf;
    Teuchos::RCP< ::ROL::RiskMeasure<RealT> > rmobj;
    bool storage = true;
    Teuchos::RCP< ::ROL::RiskNeutralObjective<RealT> > neutobj;
    Teuchos::RCP< ::ROL::RiskAverseObjective<RealT> > riskobj;
    Teuchos::RCP< ::ROL::RiskNeutralObjective<RealT> > neutpen;
    Teuchos::RCP< ::ROL::RiskAverseObjective<RealT> > riskpen;    
    Teuchos::RCP< ::ROL::BoundConstraint<RealT> >bconcvar = Teuchos::rcp(new ::ROL::CVaRBoundConstraint<RealT>(bcon));
    if(!use_bcon)
      bconcvar->deactivate();    

    //****************************** Trust Region ********************************************//
    if (use_TR)
    {
      Teuchos::RCP< ::ROL::StatusTest<double> > status_tr;
      Teuchos::RCP< ::ROL::Step<double> > step_tr;
      Teuchos::RCP< ::ROL::Algorithm<double> > algo_tr;
      std::clock_t timer_tr = std::clock();
      
      *outStream << "\nSOLVE DETERMINISTIC MEAN-VALUE PROBLEM\n";
      z1.set(z);
      con->setParameter(sampleMean);
      objVec.clear();
      objVec.push_back(robj);
#if OBJTYPE==1
      objVec.push_back(rpen);
      types.push_back(1);
#endif
      SumObjective<RealT> fullobjDET(objVec,types);
      con->checkSolve(u,z,c,true,*outStream);
      if(do_checks)
      {
        *outStream << "Derivatives of reduced objective" << std::endl; 
        robj->checkGradient(z1,z1,yz,true,*outStream);
        robj->checkHessVec(z1,z1,yz,true,*outStream);
#if OBJTYPE==1
        *outStream << "Derivatives of reduced penalty" << std::endl; 
        rpen->checkGradient(z1,z1,yz,true,*outStream);
        rpen->checkHessVec(z1,z1,yz,true,*outStream);
#endif
      }
      status_tr = Teuchos::rcp(new ::ROL::StatusTest<RealT>(gtol, stol, maxit));    
      step_tr = Teuchos::rcp(new ::ROL::TrustRegionStep<RealT>(*parlist));
      algo_tr = Teuchos::rcp(new ::ROL::Algorithm<RealT>(step_tr,status_tr,false));
      algo_tr->run(z1,fullobjDET,*bcon,true,*outStream);
      *outStream << "\n Solution " << std::endl;
      for (int i=0;i<nz;i++)
      {
        *outStream << paramNameVec_[i] << " = " << (*z1_rcp)[i] << std::endl;
      }

      *outStream << "\nSOLVE EXPECTED VALUE WITH TRUST REGION\n";
      out << "Using " << nSamp << " samples" << std::endl;
      z2.set(z1);
      *outStream << "\n Initial guess " << std::endl;
      for (int i=0;i<nz;i++)
        *outStream << paramNameVec_[i] << " = " << (*z2_rcp)[i] << std::endl;
      neutobj = Teuchos::rcp(new ::ROL::RiskNeutralObjective<RealT>(robj,samplerSolve,storage));
      objVec.clear();
      objVec.push_back(neutobj);
      neutpen = Teuchos::rcp(new ::ROL::RiskNeutralObjective<RealT>(rpen,samplerSolve,storage));
      objVec.push_back(neutpen);
      SumObjective<RealT> fullobjEXP(objVec,types);
      if(do_checks)
      {
        out << "Derivatives of risk neutral objective" << std::endl;
        neutobj->checkGradient(z2,z2,yz,true,*outStream);
        neutobj->checkHessVec(z2,z2,yz,true,*outStream);
#if OBJTYPE==1
        out << "Derivatives of risk neutral penalty" << std::endl;
        neutpen->checkGradient(z2,z2,yz,true,*outStream);
        neutpen->checkHessVec(z2,z2,yz,true,*outStream);
#endif
      }
      status_tr = Teuchos::rcp(new ::ROL::StatusTest<RealT>(gtol, stol, maxit));    
      step_tr = Teuchos::rcp(new ::ROL::TrustRegionStep<RealT>(*parlist));
      algo_tr = Teuchos::rcp(new ::ROL::Algorithm<RealT>(step_tr,status_tr,false));
      algo_tr->run(z2,fullobjEXP,*bcon,true,*outStream);
      *outStream << "\n Solution " << std::endl;
      for (int i=0;i<nz;i++)
      {
        *outStream << paramNameVec_[i] << " = " << (*z2_rcp)[i] << std::endl;
      }

      *outStream << "\nSOLVE SMOOTHED CONDITIONAL VALUE AT RISK WITH TRUST REGION\n";
      *outStream << "prob = " << prob << ", coeff = " << coeff << ", gamma = " << gamma << std::endl;  
      out << "Using " << nSamp << " samples" << std::endl;
      z3.set(z2);
      *outStream << "\n Initial guess " << std::endl;
      for (int i=0;i<nz;i++)
        *outStream << paramNameVec_[i] << " = " << (*z3_rcp)[i] << std::endl;
      // Build CVaR objective function
      std::vector<RealT> data(2,0.0);
      data[0] = -0.5; data[1] = 0.5;
      dist = Teuchos::rcp( new ::ROL::Distribution<RealT>(::ROL::DISTRIBUTION_PARABOLIC,data) );
      pf   = Teuchos::rcp( new ::ROL::PlusFunction<RealT>(dist,gamma) );
      rmobj  = Teuchos::rcp( new ::ROL::CVaR<RealT>(prob,coeff,pf) );
#if 1
      neutobj = Teuchos::rcp(new ::ROL::RiskNeutralObjective<RealT>(robj,samplerSolve,storage));
      objVec.clear();
      objVec.push_back(neutobj);
      if(coeff > 0)
        types[0] = 0; // if riskpen
#endif
#if 0
      riskobj = Teuchos::rcp(new ::ROL::RiskAverseObjective<RealT>(robj,rmobj,samplerSolve,storage));
      objVec.clear();
      objVec.push_back(riskobj);
#endif
#if OBJTYPE==1
      if(coeff == 0)
      {
        neutpen = Teuchos::rcp(new ::ROL::RiskNeutralObjective<RealT>(rpen,samplerSolve,storage));
        objVec.push_back(neutpen);
        //types.push_back(0); // if riskobj
        types.push_back(1); // if neutobj
      }
      else
      {
        riskpen = Teuchos::rcp(new ::ROL::RiskAverseObjective<RealT>(rpen,rmobj,samplerSolve,storage));
        objVec.push_back(riskpen);
        //types.push_back(1);
        types[1] = 1;
      }
#endif
      SumObjective<RealT> fullobj(objVec,types);
      status_tr = Teuchos::rcp(new ::ROL::StatusTest<RealT>(gtol, stol, maxit));    
      step_tr = Teuchos::rcp(new ::ROL::TrustRegionStep<RealT>(*parlist));
      algo_tr = Teuchos::rcp(new ::ROL::Algorithm<RealT>(step_tr,status_tr,false));
      // Build CVaR vectors
      RealT z3v = (RealT)rand()/(RealT)RAND_MAX;
      RealT yzv = (RealT)rand()/(RealT)RAND_MAX;
      ::ROL::CVaRVector<RealT> z3c(z3v,z3p);
      ::ROL::CVaRVector<RealT> yzc(yzv,yzp);
      if(do_checks)
      {
        out << "Derivatives of risk neutral objective" << std::endl;
        neutobj->checkGradient(z3,z3,yz,true,*outStream);
        neutobj->checkHessVec(z3,z3,yz,true,*outStream);
        // out << "Derivatives of risk averse objective" << std::endl;
        // riskobj->checkGradient(zc,zc,yzc,true,*outStream);
        // riskobj->checkHessVec(zc,zc,yzc,true,*outStream);
#if OBJTYPE==1
        if(coeff == 0)
        {
          out << "Derivatives of risk neutral penalty" << std::endl;
          neutpen->checkGradient(z3,z3,yz,true,*outStream);
          neutpen->checkHessVec(z3,z3,yz,true,*outStream);
        }
        else
        {
          out << "Derivatives of risk averse penalty" << std::endl;
          riskpen->checkGradient(z3c,z3c,yzc,true,*outStream);
          riskpen->checkHessVec(z3c,z3c,yzc,true,*outStream);
        }
#endif
      }
      // Run ROL algorithm
      if(coeff > 0)
        algo_tr->run(z3c,fullobj,*bconcvar,true,*outStream);
      else
        algo_tr->run(z3,fullobj,*bcon,true,*outStream);
      *outStream << "\n Solution " << std::endl;
      for (int i=0;i<nz;i++)
      {
        *outStream << paramNameVec_[i] << " = " << (*z3_rcp)[i] << std::endl;
      }
      *outStream << "var = " << z3v << std::endl;
      *outStream << "Trust-Region required " << (std::clock()-timer_tr)/(RealT)CLOCKS_PER_SEC
                 << " seconds.\n";
    }

    //***************************** Bundle Method *********************************//
    if(use_bundle)
    {
      *outStream << "\nSOLVE NONSMOOTH CVAR PROBLEM WITH BUNDLE TRUST REGION\n";
      // Build CVaR objective function
      *outStream << "prob = " << prob << ", coeff = " << coeff << std::endl;
      dist = Teuchos::rcp( new ::ROL::Distribution<RealT>(::ROL::DISTRIBUTION_DIRAC) );
      pf   = Teuchos::rcp( new ::ROL::PlusFunction<RealT>(dist,1.0) );
      rmobj  = Teuchos::rcp( new ::ROL::CVaR<RealT>(prob,coeff,pf) );
      neutobj = Teuchos::rcp(new ::ROL::RiskNeutralObjective<RealT>(robj,samplerSolve,storage));
      objVec.clear();
      objVec.push_back(neutobj);
      if(coeff > 0)
        types[0] = 0; // if riskpen
#if OBJTYPE==1
      if(coeff == 0)
      {
        neutpen = Teuchos::rcp(new ::ROL::RiskNeutralObjective<RealT>(rpen,samplerSolve,storage));
        objVec.push_back(neutpen);
        //types.push_back(0); // if riskobj
        types.push_back(1); // if neutobj
      }
      else
      {
        riskpen = Teuchos::rcp(new ::ROL::RiskAverseObjective<RealT>(rpen,rmobj,samplerSolve,storage));
        objVec.push_back(riskpen);
        types.push_back(1);
      }
#endif
      SumObjective<RealT> fullobj(objVec,types);
      // Build CVaR vectors
      RealT zv = (RealT)rand()/(RealT)RAND_MAX;
      ::ROL::CVaRVector<RealT> zc(zv,zp);
      // Run ROL algorithm
      parlist->set("Bundle Step: Epsilon Solution Tolerance",gtol);
      ::ROL::BundleStatusTest<RealT> status_bm(gtol, maxit);    
      ::ROL::BundleStep<RealT> step_bm(*parlist);
      ::ROL::Algorithm<RealT> algo_bm(Teuchos::rcp(&step_bm,false),Teuchos::rcp(&status_bm,false),false);
      std::clock_t timer_bm = std::clock();
      if(coeff > 0)
        algo_bm.run(zc,fullobj,*bconcvar,true,*outStream);
      else
        algo_bm.run(z,fullobj,*bcon,true,*outStream);
      *outStream << "\n Solution " << std::endl;
      for (int i=0;i<nz;i++)
      {
        *outStream << paramNameVec_[i] << " = " << (*z_rcp)[i] << std::endl;
      }
      *outStream << "var = " << zv << std::endl;
      *outStream << "Bundle Method required " << (std::clock()-timer_bm)/(RealT)CLOCKS_PER_SEC
                 << " seconds.\n";
    }

    // Post-processing (UQ)
    //****************************** Final objective values and final states (UQ) *******************//
    std::vector<std::string> output_data_UQ;
    output_data_UQ.push_back("data_UQ_DET.txt");
    output_data_UQ.push_back("data_UQ_EV.txt");    
    output_data_UQ.push_back("data_UQ_CVaR.txt");
    std::vector<std::string> output_plot_UQ;
    output_plot_UQ.push_back("plot_UQ_DET.txt");
    output_plot_UQ.push_back("plot_UQ_EV.txt");    
    output_plot_UQ.push_back("plot_UQ_CVaR.txt");
    std::vector< ::ROL::StdVector<RealT> > solutions;
    solutions.push_back(z1);
    solutions.push_back(z2);
    solutions.push_back(z3);
    // Test with the rest of "good" samples
    int nTestSamp = nGSamp - nSamp;
    out << "Using " << nTestSamp << " to approximate densities" << std::endl;
    for(int k=0;k<solutions.size();k++)
    {
      std::ofstream data_UQ(output_data_UQ[k].c_str());
      std::ofstream plot_UQ(output_plot_UQ[k].c_str());
      sumlin = 0; sumlin2 = 0;
      sumgain = 0; sumgain2 = 0;
      for(int i=nSamp;i<nGSamp;i++)
      {
        con->setParameter(samplerTest->getMyPoint(i));
        con->solve(u,solutions[k],tol);
        oval = obj.value(u,solutions[k],tol);
        for(int j=0;j<nc;j++)
          plot_UQ << std::scientific << std::setprecision(8) << (*(*uup)[j])[3] << ' ';
        plot_UQ << '\n';
        gain = ((*(*uup)[nc-1])[3]-(*(*uup)[0])[3])/2.0;
        sumgain += gain;
        sumgain2 += gain*gain;
        for(int j=0;j<uncertainParams_.size();j++)
        {
          data_UQ << std::scientific << std::setprecision(8) << (samplerTest->getMyPoint(i))[j] << ' ';
        }
        data_UQ << oval << ' ' << gain << '\n';
        sumlin += oval;
        sumlin2 += oval*oval;
      }
      plot_UQ.close();
      meanlin  = sumlin/nTestSamp;
      meangain = sumgain/nTestSamp;
      varlin  = (sumlin2  - nTestSamp*meanlin*meanlin)/(RealT)(nTestSamp-1);
      vargain = (sumgain2 - nTestSamp*meangain*meangain)/(RealT)(nTestSamp-1);
      data_UQ << "\n ------------------------------------------ \n";
      data_UQ << "mean(lin) = " << meanlin << " " << '\n';
      data_UQ << "sqrt(sample variance)(lin) = " << std::sqrt(varlin) << '\n';
      data_UQ << "mean(gain) = " << meangain << '\n';
      data_UQ << "sqrt(sample variance)(gain) = " << std::sqrt(vargain) << '\n';
      data_UQ.close();
    }

    //********************************************************************************************// 
    
#else 
    //**********************************************************************************************// 
    //******************************** Deterministic problem **************************************//
    //********************************************************************************************// 
    if(use_bundle)
    {
      parlist->set("Bundle Step: Epsilon Solution Tolerance",gtol);
      ::ROL::BundleStatusTest<RealT> status_bm(gtol, maxit);    
      ::ROL::BundleStep<RealT> step_bm(*parlist);
      ::ROL::Algorithm<RealT> algo_bm(Teuchos::rcp(&step_bm,false),Teuchos::rcp(&status_bm,false),false);
      std::clock_t timer_bm = std::clock();
      objVec.clear();
      objVec.push_back(robj);
#if OBJTYPE==1
      objVec.push_back(rpen);
      types.push_back(1);
#endif
      SumObjective<RealT> fullobj(objVec,types);
      algo_bm.run(z,fullobj,*bcon,true,*outStream);
      *outStream << "\n Solution " << std::endl;
      for (int i=0;i<nz;i++)
      {
        *outStream << paramNameVec_[i] << " = " << (*z_rcp)[i] << std::endl;
      }
      *outStream << "Bundle Method required " << (std::clock()-timer_bm)/(RealT)CLOCKS_PER_SEC
                 << " seconds.\n";
    }
    if (use_TR)
    {
      // Trust Region
      ::ROL::StatusTest<RealT> status_tr(gtol, stol, maxit);    
      ::ROL::TrustRegionStep<RealT> step_tr(*parlist);
      ::ROL::Algorithm<RealT> algo_tr(Teuchos::rcp(&step_tr,false),Teuchos::rcp(&status_tr,false),false);
      std::clock_t timer_tr = std::clock();
      objVec.clear();
      objVec.push_back(robj);
#if OBJTYPE==1
      objVec.push_back(rpen);
      types.push_back(1);
#endif
      SumObjective<RealT> fullobj(objVec,types);
      algo_tr.run(z,fullobj,*bcon,true,*outStream);
      *outStream << "\n Solution " << std::endl;
      for (int i=0;i<nz;i++)
      {
        *outStream << paramNameVec_[i] << " = " << (*z_rcp)[i] << std::endl;
      }
      *outStream << "Trust-Region required " << (std::clock()-timer_tr)/(RealT)CLOCKS_PER_SEC
                 << " seconds.\n";
    }
    if (use_lsearch)
    {
      // Line Search
      ::ROL::StatusTest<RealT> status_ls(gtol, stol, maxit);    
      ::ROL::LineSearchStep<RealT> step_ls(*parlist);
      ::ROL::Algorithm<RealT> algo_ls(Teuchos::rcp(&step_ls,false),Teuchos::rcp(&status_ls,false),false);
      std::clock_t timer_ls = std::clock();
      objVec.clear();
      objVec.push_back(robj);
#if OBJTYPE==1
      objVec.push_back(rpen);
      types.push_back(1);
#endif
      SumObjective<RealT> fullobj(objVec,types);
      algo_ls.run(z,fullobj,*bcon,true,*outStream);
      *outStream << "\n Solution " << std::endl;
      for (int i=0;i<nz;i++)
      {
        *outStream << paramNameVec_[i] << " = " << (*z_rcp)[i] << std::endl;
      }
      *outStream << "Line-Search required " << (std::clock()-timer_ls)/(RealT)CLOCKS_PER_SEC
                 << " seconds.\n";
    }
    if (use_sqp)
    {
      // SQP.
      Teuchos::RCP<std::vector<RealT> > gz_rcp = Teuchos::rcp( new std::vector<RealT> (nz, 0.0) );
      ::ROL::StdVector<RealT> gz(gz_rcp);
      Teuchos::RCP< ::ROL::Vector<RealT> > gzp = Teuchos::rcp(&gz,false);
      Linear::ROL_XyceVector<RealT> gu(nc,*analysisManager_.getDataStore()->nextSolutionPtr);
      Teuchos::RCP< ::ROL::Vector<RealT> > gup = Teuchos::rcp(&gu,false);
      ::ROL::Vector_SimOpt<RealT> g(gup,gzp);
      Linear::ROL_XyceVector<RealT> cc(nc,*analysisManager_.getDataStore()->nextSolutionPtr);
      Linear::ROL_XyceVector<RealT> l(nc,*analysisManager_.getDataStore()->nextSolutionPtr);
      RealT ctol = 1.e-16;
      ::ROL::ConstraintStatusTest<RealT> status_sqp(gtol,ctol,stol,maxit);
      ::ROL::CompositeStep<RealT> step_sqp(*parlist);
      ::ROL::Algorithm<RealT> algo_sqp(Teuchos::rcp(&step_sqp,false),Teuchos::rcp(&status_sqp, false),false);
      std::clock_t timer_sqp = std::clock();
      objVec.clear();
      objVec.push_back(obj);
#if OBJTYPE==1
      objVec.push_back(pen);
      types.push_back(1);
#endif
      SumObjective<RealT> fullobj(objVec,types);
      algo_sqp.run(x,g,l,cc,fullobj,*con,true,*outStream);
      out << "\n Solution " << std::endl;
      for (int i=0;i<nz;i++)
      {
        *outStream << paramNameVec_[i] << " = " << (*z_rcp)[i] << std::endl;
      }
      out << "Composite-Step SQP required " << (std::clock()-timer_sqp)/(RealT)CLOCKS_PER_SEC
          << " seconds.\n";
    }
#endif

  }  
  catch (std::logic_error err) 
  {
    //out << err.what() << "\n";
    errorFlag = -1000;
  }; // end try       
  
  // if (errorFlag != 0)
  //   out << "End Result: TEST FAILED\n";
  // else
  //   out << "End Result: TEST PASSED\n";
  
  doFree(); // deallocate solution and sensitivity arrays
  
  return status;
}
#else
//-----------------------------------------------------------------------------
//TT: do the same as doLoopProcess
//-----------------------------------------------------------------------------
bool ROL::runROLAnalysis()
{
  bool integration_status = true;
  TimeIntg::DataStore & ds = *(analysisManager_.getDataStore()); // TT
  
  //Util::publish<StepEvent>(analysisManager_, StepEvent(StepEvent::INITIALIZE, stepSweepVector_, stepLoopSize_)); // TT: moved here from doInit();
  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::INITIALIZE, AnalysisEvent::DC));
  
  // StepEvent step_event(StepEvent::STEP_STARTED, stepSweepVector_, currentStep);

  //for (int currentStep = 0; currentStep < stepLoopSize_; ++i)
  int currentStep = 0;
  int finalStep = stepLoopSize_;
  while (currentStep < finalStep) // TT: I suspect that this might be necessary
  {
    outputManagerAdapter_.setDCAnalysisStepNumber(currentStep);
      
    // Tell the manager if any of our sweeps are being reset in this loop iteration.
    bool reset = updateSweepParams(loader_, currentStep, stepSweepVector_.begin(), stepSweepVector_.end(), false);// TT: update sweep parameters; location N_ANP_SweepParam.C
      
    analysisManager_.setSweepSourceResetFlag(reset);// TT: sets sweep source flag to reset in AnalysisManager
          
    outputManagerAdapter_.setStepSweepVector(stepSweepVector_);// TT: not sure what this is for
      
    // if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
    for (SweepVector::const_iterator it = stepSweepVector_.begin(), end = stepSweepVector_.end(); it != end; ++it)
      Xyce::dout() << "ROL DC Sweep # " << currentStep <<"\t" << (*it);

    // TT
    if (currentStep != 0 && reset)
    {
      analysisManager_.getDataStore()->setZeroHistory();
      initializeSolution_();
    }
    
          
    // TT: the following line causes currSol and nextSol reset to zero
    // Util::publish<StepEvent>(analysisManager_, StepEvent(StepEvent::STEP_STARTED, stepSweepVector_, currentStep)); 

    static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_STARTED, AnalysisEvent::DC, 0.0, currentStep));
    //step_event.state = StepEvent::STEP_STARTED;
    //Util::publish<StepEvent>(analysisManager_, step_event);
  
    takeStep_(); // TT: take integration step

    // Set things up for the next time step, based on if this one was
    // successful.
    if (analysisManager_.getStepErrorControl().stepAttemptStatus)
    {
      static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_SUCCESSFUL, AnalysisEvent::DC, 0.0, currentStep));
      doProcessSuccessfulStep();
      // collect solutions here
      *(solutionPtrVector_[currentStep]) = *(ds.currSolutionPtr); 
    }
    else // stepAttemptStatus  (ie do this if the step FAILED)
    {
      static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_FAILED, AnalysisEvent::DC, 0.0, currentStep));
      doProcessFailedStep();
    }
      
    currentStep = stepNumber;
  } // end of sweep loop

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::FINISH, AnalysisEvent::DC));

  // Util::publish<StepEvent>(analysisManager_, StepEvent(StepEvent::STEP_COMPLETED, stepSweepVector_, currentStep));
  
  // step_event.state_ = StepEvent::STEP_COMPLETED;
  // step_event.finalSimTime_ = getTIAParams().finalTime;
  // Util::publish<StepEvent>(analysisManager_, step_event);

  return integration_status;
}

#endif


namespace {

//-----------------------------------------------------------------------------
// Class         : ROLFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jan 29 12:53:02 2015
//-----------------------------------------------------------------------------
///
/// Factory for parsing ROL parameters from the netlist and creating ROL analysis.
///
class ROLFactory : public Util::Factory<AnalysisBase, ROL>
{
public:
  //-----------------------------------------------------------------------------
  // Function      : ROLFactory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:54:09 2015
  //-----------------------------------------------------------------------------
  ///
  /// Constructs the ROL analysis factory
  ///
  /// @invariant Stores the results of parsing.  Multiple ROL analysis options may be
  /// applied and each generates and additional step.
  ///
  /// @invariant The existence of the parameters specified in the constructor cannot
  /// change.
  ///
  /// @param analysis_manager 
  /// @param linear_system 
  /// @param nonlinear_manager 
  ///
  ROLFactory(
     Analysis::AnalysisManager &         analysis_manager,
     Linear::System &                    linear_system,
     Nonlinear::Manager &                nonlinear_manager,
     Loader::Loader &                    loader,
     Topo::Topology &                    topology,
     IO::InitialConditionsManager &      initial_conditions_manager)
    : Util::Factory<AnalysisBase, ROL>(),
    analysisManager_(analysis_manager),
    linearSystem_(linear_system),
    nonlinearManager_(nonlinear_manager),
    loader_(loader),
    topology_(topology),
    initialConditionsManager_(initial_conditions_manager),
    stepSweepAnalysisOptionBlock_(), // TT
    timeIntegratorOptionBlock_() // TT
  {}

  virtual ~ROLFactory()
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
  /// Create a new ROL analysis and applies the analysis and time integrator option blocks.
  ///
  /// @return new ROL analysis object
  ///
  ROL *create() const
  {
    // analysisManager_.setAnalysisMode(ANP_MODE_ROL); // TT: 
    analysisManager_.setAnalysisMode(ANP_MODE_DC_SWEEP); // TT: 
    ROL *step = new ROL(analysisManager_, nonlinearManager_, loader_, linearSystem_, topology_, initialConditionsManager_);
    for (std::vector<Util::OptionBlock>::const_iterator it = stepSweepAnalysisOptionBlock_.begin(), end = stepSweepAnalysisOptionBlock_.end(); it != end; ++it)
      step->setAnalysisParams(*it);
    step->setTimeIntegratorOptions(timeIntegratorOptionBlock_);// TT

    return step;
  }

  //-----------------------------------------------------------------------------
  // Function      : setROLAnalysisOptionBlock
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
  void setROLAnalysisOptionBlock(const Util::OptionBlock &option_block)
  {
    for (std::vector<Util::OptionBlock>::iterator it = stepSweepAnalysisOptionBlock_.begin(), end = stepSweepAnalysisOptionBlock_.end(); it != end; ++it)
    {
      if (Util::compareParamLists(option_block, *it))
      {
        (*it) = option_block;
        return;
      }
    }

    // save the new one.
    stepSweepAnalysisOptionBlock_.push_back(option_block); // save a copy for later.

    //TT
    std::string param_name;
    for (std::vector<Util::OptionBlock>::iterator it = stepSweepAnalysisOptionBlock_.begin(), end = stepSweepAnalysisOptionBlock_.end(); it != end; ++it)
    {
      param_name = (*it).getName();
      std::cout << "parameter = " << param_name << std::endl; 
    }
  }
  

  // TT: blindly copied from DCSweep
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
    
public:
  AnalysisManager &             analysisManager_;
  Linear::System &              linearSystem_;
  Nonlinear::Manager &          nonlinearManager_;
  Loader::Loader &              loader_;
  Topo::Topology &              topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;

private:
  std::vector<Util::OptionBlock>        stepSweepAnalysisOptionBlock_;
  Util::OptionBlock                     timeIntegratorOptionBlock_; // TT
};
  
// .ROL
struct ROLAnalysisReg : public IO::PkgOptionsReg
{
  ROLAnalysisReg(
     ROLFactory &             factory)
    : factory_(factory)
  {}

  bool operator()(const Util::OptionBlock &option_block)
  {
    factory_.setROLAnalysisOptionBlock(option_block);

    factory_.analysisManager_.addAnalysis(&factory_);

    return true;
  }

  ROLFactory &               factory_;
};

//-----------------------------------------------------------------------------
// Function      : extractROLData
// Purpose       : Extract the parameters from a netlist .ROL line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 10/30/2003
//-----------------------------------------------------------------------------
// TT: this function is copied from STEP, the parameters need to be modified
// TT: I will want here a way to extract optimization parameters
bool extractROLData(
   IO::PkgOptionsMgr &           options_manager,
   IO::CircuitBlock &            circuit_block,
   const std::string &           netlist_filename,
   const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("ROL", Util::OptionBlock::ALLOW_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

  int numFields = parsed_line.size();

  // First check if the type has been explicitly set.
  // If not, set it to the default, LIN.
  int pos1=1;

  bool typeExplicitSetLinDecOct = false;
  bool typeExplicitSetList = false;
  std::string type("LIN");
  while ( pos1 < numFields )
  {
    ExtendedString stringVal ( parsed_line[pos1].string_ );
    stringVal.toUpper ();
    if (stringVal == "LIN" ||
        stringVal == "DEC" ||
        stringVal == "OCT")
    {
      typeExplicitSetLinDecOct = true;
      type = stringVal;
    }
    else if (stringVal == "LIST")
    {
      typeExplicitSetList = true;
      type = stringVal;
    }

    ++pos1;
  }

  // Check that the minimum required number of fields are on the line.
  int offset = 1;
  if (typeExplicitSetLinDecOct)
  {
    offset = 2;
  }

  if (!typeExplicitSetList)// if this is a list, number of fields is arbitrary.
  {
    if ( (numFields-offset)%4 != 0 )
    {
      Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
        << ".ROL line not formatted correctly.";
      return false;
    }
  }

  int linePosition = 1;   // Start of parameters on .param line.
  Util::Param parameter("", "");

  // Add the type (which was determined above) to the parameter list.
  parameter.setTag( "TYPE" );
  parameter.setVal( type );
  option_block.addParam( parameter );

  if (type=="LIN")
  {
    if (typeExplicitSetLinDecOct) linePosition=2;
    while ( linePosition < numFields )
    {
      parameter.setTag( "PARAM" );
      parameter.setVal(std::string(ExtendedString(parsed_line[linePosition].string_).toUpper()));
      option_block.addParam( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "START" );
      parameter.setVal( parsed_line[linePosition].string_ );
      option_block.addParam( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "STOP" );
      parameter.setVal( parsed_line[linePosition].string_ );
      option_block.addParam( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "ROL" );
      parameter.setVal( parsed_line[linePosition].string_ );
      option_block.addParam( parameter );
      ++linePosition;     // Advance to next parameter.
    }
  }
  else if (type=="DEC")
  {
    if (typeExplicitSetLinDecOct) linePosition=2;

    while ( linePosition < numFields )
    {
      parameter.setTag( "PARAM" );
      parameter.setVal(std::string(ExtendedString(parsed_line[linePosition].string_).toUpper()));
      option_block.addParam( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "START" );
      parameter.setVal( parsed_line[linePosition].string_ );
      option_block.addParam( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "STOP" );
      parameter.setVal( parsed_line[linePosition].string_ );
      option_block.addParam( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "NUMROLS" );
      parameter.setVal( parsed_line[linePosition].string_ );
      option_block.addParam( parameter );
      ++linePosition;     // Advance to next parameter.
    }
  }
  else if (type=="OCT")
  {
    if (typeExplicitSetLinDecOct) linePosition=2;

    while ( linePosition < numFields )
    {
      parameter.setTag( "PARAM" );
      parameter.setVal(std::string(ExtendedString(parsed_line[linePosition].string_).toUpper()));
      option_block.addParam( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "START" );
      parameter.setVal( parsed_line[linePosition].string_ );
      option_block.addParam( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "STOP" );
      parameter.setVal( parsed_line[linePosition].string_ );
      option_block.addParam( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "NUMROLS" );
      parameter.setVal( parsed_line[linePosition].string_ );
      option_block.addParam( parameter );
      ++linePosition;     // Advance to next parameter.
    }

  }
  else if (type=="LIST")
  {
    parameter.setTag( "PARAM" );
    parameter.setVal(std::string(ExtendedString(parsed_line[1].string_).toUpper()));
    option_block.addParam( parameter );

    int linePosition=3;
    while (linePosition<numFields)
    {
      parameter.setTag( "VAL" );
      parameter.setVal( parsed_line[linePosition].string_ );
      option_block.addParam( parameter );
      ++linePosition;
    }
  }
  else
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << ".ROL line contains an unrecognized type";
  }

  circuit_block.addOptions(option_block);

  return true;
}

} // namespace <unnamed>


bool registerROLFactory(
   FactoryBlock &        factory_block)
{
  ROLFactory *factory = new ROLFactory(factory_block.analysisManager_, factory_block.linearSystem_, factory_block.nonlinearManager_, factory_block.loader_, factory_block.topology_, factory_block.initialConditionsManager_);

  addAnalysisFactory(factory_block, factory);

  factory_block.optionsManager_.addCommandParser(".ROL", extractROLData);

  factory_block.optionsManager_.addCommandProcessor("ROL", new ROLAnalysisReg(*factory));

  factory_block.optionsManager_.addOptionsProcessor("TIMEINT", IO::createRegistrationOptions(*factory, &ROLFactory::setTimeIntegratorOptionBlock));

  return true;
}

} // namespace Analysis
} // namespace Xyce
 

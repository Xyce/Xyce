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
#define OBJTYPE 0 // 0 - L2Norm, 1 - amplifier circuit
#endif

#include "ROL_BoundConstraint.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Problem.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Solver.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Vector.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"

#include "ROL_XyceVector.hpp"
#include <N_ANP_ROL_DC_Optimization.h>

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

#ifdef Xyce_ROL
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
    
  std::string msg;
  bool status = true;
  int errorFlag = 0;

  try
  { 

    int nz = numParams_;     // Number of optimization variables                                           
    // Number of simulation variables
    int nu = analysisManager_.getDataStore()->nextSolutionPtr->globalLength();
    int nc = stepLoopSize_;  // Number of constraints.

    RealT tol = 1.e-12;
   
    // STEP 1A: Initialize vectors.  //////////////////////////////////////////

    // Configure the optimization vector.
    auto p = ::ROL::makePtr<  std::vector   <RealT>>(nz, 0.0);
    auto z = ::ROL::makePtr<::ROL::StdVector<RealT>>(p);
    // p := pointer to a standard vector representing z

    // Read in z's initial guess and upper and lower bounds.
    std::ifstream file;  file.open("parameters.txt");
    std::vector<RealT> zUpperBoundVector(nz, 0.0);
    std::vector<RealT> zLowerBoundVector(nz ,0.0);
    std::string strStreamIn;
    RealT       numStreamIn;
    getline(file, strStreamIn);  // Skips the first line
    for(int i = 0; i < nz; i++)
    {
      file                 >> strStreamIn;
      file                 >> strStreamIn;
      file                 >> numStreamIn;  // Initial guess
      (*p)[i]              =  numStreamIn;
      file                 >> numStreamIn;  // Lower bound
      zLowerBoundVector[i] =  numStreamIn;
      file                 >> numStreamIn;  // Upper bound
      zUpperBoundVector[i] =  numStreamIn;
    }
    file.close();

    // Configure the state and constraint vectors.
    status = doAllocations(nc, nz);
    auto u = ::ROL::makePtr<Linear::ROL_XyceVector<RealT>>(statePtrVector_);
    auto l = ::ROL::makePtr<Linear::ROL_XyceVector<RealT>>(constraintPtrVector_);

    // STEP 1B: Initialize parameters.  ///////////////////////////////////////

    // Load parameters.
    ::ROL::Ptr<::ROL::ParameterList> parlist 
      = ::ROL::getParametersFromXmlFile("input.xml");
    std::string outName     = parlist->get("Output File", "rol_output.txt");
    bool useBoundConstraint = parlist->get("Use Bound Constraints", true);
    bool useScale           = parlist->get("Use Scaling For Epsilon-Active Sets", true);
    // TODO (asjavee): Remove hardwired code. Eventually we shouldn't need the 
    //                 subsequent lines in this code block.
    // bool doChecks           = parlist->get("Do Checks",true);
    // bool useSQP             = parlist->get("Use SQP", true);
    // bool useLineSearch      = parlist->get("Use Line Search", true);    
    // bool useTrustRegion     = parlist->get("Use Trust Region", true);
    // bool useBundleMethod    = parlist->get("Use Bundle Method", true);
    // Parameters for the amplifier problem. 
    int ptype               = parlist->get("Penalty Type", 1);
    RealT alpha             = parlist->get("Penalty Parameter", 1.0e-4);
    RealT ampl              = parlist->get("Amplifier Gain", 4.0);

    // Configure the output stream.
    std::ofstream out(outName.c_str());
    // TODO (asjavee): I don't think we need the line below; ofstream inherits 
    //                 from ostream.
    // Teuchos::RCP<std::ostream> outStream = Teuchos::rcp(&out,false);

    // STEP 2: Create the SimOpt equality constraint.  ////////////////////////

    auto con = ::ROL::makePtr<EqualityConstraint_ROL_DC<RealT>>
               (
                 nu, nc, nz, analysisManager_, 
                             analysisManager_.getNonlinearEquationLoader(),
                             nonlinearManager_.getNonlinearSolver(),
                             linearSystem_,
                             *this
               );
    // Here, con is the constraint c(u, z) = 0. Its solve member function 
    // specifies u as a function of z (u = S(z)); see slide 26 of
    //   <https://trilinos.github.io/pdfs/ROL.pdf>.

    con->solve(*l, *u, *z, tol);

    // STEP 3: Build our objective.  /////////////////////////////////////////

    // Create a full space objective (as well as a penalty when solving the
    // amplifier problem).
#if OBJTYPE==0
    auto obj = ::ROL::makePtr<Objective_DC_L2Norm<RealT>>(1.0e-4, nc, nz);
#elif OBJTYPE==1
    auto obj = ::ROL::makePtr<Objective_DC_AMP<RealT>>(nc, nz);
    auto pen = ::ROL::makePtr<Penalty_DC_AMP  <RealT>>(ptype, alpha, ampl, nc, nz);
#endif
   
    // Create a reduced objective (and penalty) from the full space counterpart.
    auto robj = ::ROL::makePtr<::ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, u, z, l);
#if OBJTYPE==1
    auto rpen = ::ROL::makePtr<::ROL::Reduced_Objective_SimOpt<RealT>>(pen, con, u, z, l);
#endif

    // Create the (possibly penalized) objective that we wish to optimize.
    std::vector<Teuchos::RCP<::ROL::Objective<RealT>>> objVec(1, robj);
#if OBJTYPE==1
    objVec.push_back(rpen);
#endif
    std::vector<bool> types(objVec.size(), true);
    auto pobj = ::ROL::makePtr<SumObjective<RealT>>(objVec, types);

    // STEP 4: Create z's BoundConstraint.  ///////////////////////////////////
    
    RealT scale = 1.0;
    if(useScale)
    {
      // Evaluate the reduced objective gradient at the initial value of z.
      auto g0P = ::ROL::makePtr<  std::vector   <RealT>>(nz, 0.0);
      auto g0  = ::ROL::makePtr<::ROL::StdVector<RealT>>(g0P);
      robj->gradient(*g0, *z, tol);

      // std::cout << std::scientific << "Norm of initial gradient = " << g0->norm() << "\n";
    
      scale = 1.0e-2/(g0->norm());
    }

    // std::cout << std::scientific << "Scaling: " << scale << "\n";
    
    auto bnd = ::ROL::makePtr<BoundConstraint_ROL_DC<RealT>>(scale, zLowerBoundVector, zUpperBoundVector);
    
    // STEP 5A: Define the optimization problem.  /////////////////////////////

    auto problem = ::ROL::makePtr<::ROL::Problem<RealT>>(robj, z);
    if (useBoundConstraint)
      problem->addBoundConstraint(bnd);
    // problem->check(true);

    // STEP 5B: Define the optimization solver.  //////////////////////////////
    
    parlist->sublist("General").set("Output Level", 1);    
    parlist->sublist("Status Test").set("Gradient Tolerance", 1.e-10);
    parlist->sublist("Status Test").set("Step Tolerance", 1.e-30); 
 
    ::ROL::Solver<RealT> solver(problem, *parlist);

    // STEP 6: Solve.  ////////////////////////////////////////////////////////
    
    // Print initial guess.
    out << "\nInitial guess " << std::endl;
    for (int i = 0; i < nz; i++)
      out << paramNameVec_[i] << " = " << (*p)[i] << std::endl;

    // TODO (asjavee): Manage the default below by setting them in parlist.
    // int maxit  = parlist->get("Maximum Number of Iterations", 100);
    
    std::clock_t timer_bm = std::clock();
    
    solver.solve(std::cout);
    con->solve(*l, *u, *z, tol);

    // Print results.
    out << "\nSolution " << std::endl;
    for (int i = 0; i < nz; i++)
    {
      out << paramNameVec_[i] << " = " << (*p)[i] << std::endl;
    }
    out << "Solve required " << (std::clock()-timer_bm)/(RealT)CLOCKS_PER_SEC
         << " seconds.\n";

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
 

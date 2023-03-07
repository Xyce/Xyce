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
// Purpose       : This file contains the functions which define the time
//                 domain & integration algorithm classes.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <ctime>
#include <iostream>
#include <sstream>

#include <N_ANP_AnalysisManager.h>

#include <N_ANP_Transient.h>
#include <N_ANP_DCSweep.h>
#include <N_ANP_AC.h>
#include <N_ANP_HB.h>
#include <N_ANP_MPDE.h>
#include <N_ANP_Step.h>
#include <N_ANP_Sampling.h>
#include <N_ANP_EmbeddedSampling.h>
#ifdef Xyce_STOKHOS_ENABLE
#include <N_ANP_PCE.h>
#endif
#include <N_ANP_ROL.h>

#include <N_ANP_OutputMgrAdapter.h>
#include <N_DEV_DeviceMgr.h>
#include <N_IO_ActiveOutput.h>
#include <N_IO_CmdParse.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_Vector.h>
#include <N_LOA_Loader.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_NLS_Manager.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_fwd.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TIAParams.h>
#include <N_TIA_WorkingIntegrationMethod.h>
#include <N_TIA_TwoLevelError.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Stats.h>
#include <N_UTL_Timer.h>

#include <N_PDS_Comm.h>
#include <N_PDS_Manager.h>

#include <N_TOP_Topology.h>
#include <N_LAS_Builder.h>
#include <N_PDS_ParMap.h>

#include <expressionGroup.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : analysisModeName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Oct 21 11:22:31 2014
//-----------------------------------------------------------------------------
///
/// Returns a const char pointer to a name representing the mode.
///
/// @param mode Analysis mode to get name of
///
/// @return Name corresponding to the mode
///
//-----------------------------------------------------------------------------
const char *
analysisModeName(
  Mode          mode) 
{
  static const char * const mode_names[] = {"Invalid", "DC OP", "DC Sweep", "DC NL Poisson", "Transient", "MPDE", "HB", "AC", "NOISE", "MOR", "ROL"};

  if (mode < sizeof(mode_names)/sizeof(mode_names[0]))
    return mode_names[mode];
  else
    return mode_names[0];
}

//-----------------------------------------------------------------------------
// Function      : nonlinearAnalysisMode
// Purpose       : Convert between Nonlinear::Manager AnalysisMode enum and
//               : expanded AnalysisManager::Manager::Mode enum.
// Special Notes :
// Scope         :
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
///
/// Returns the nonlinear analysis mode given the analysis mode
///
/// @param mode Analysis mode
///
/// @return Nonlinear analysis mode
///
//-----------------------------------------------------------------------------
Nonlinear::AnalysisMode nonlinearAnalysisMode(Mode mode)
{
  Nonlinear::AnalysisMode outMode;
  if (mode == ANP_MODE_TRANSIENT)
  {
    outMode = Nonlinear::TRANSIENT;
  }
  else if (mode == ANP_MODE_DC_OP)
  {
    outMode = Nonlinear::DC_OP;
  }
  else if (mode == ANP_MODE_DC_SWEEP)
  {
    outMode = Nonlinear::DC_SWEEP;
  }
  else if (mode == ANP_MODE_DC_NLPOISSON)
  {
    outMode = Nonlinear::DC_NLPOISSON;
  }
  else if (mode == ANP_MODE_HB)
  {
    outMode = Nonlinear::HB_MODE;
  }
  else
  {
    outMode = Nonlinear::NUM_MODES; // Should be this be TRANSIENT?
  }
  return outMode;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::AnalysisManager
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
///
/// Constructs the analysis manager.  The analysis manager is responsible for 
/// the construction and control of the top level analysis.
///
/// @invariant
///
/// @param command_line              Command line that started Xyce
/// @param restart_manager              
/// @param output_manager_adapter 
/// @param analysis_stat             Base node of the analysis stats collection
///
//-----------------------------------------------------------------------------
AnalysisManager::AnalysisManager(
  const IO::CmdParse &  command_line,
  OutputMgrAdapter &    output_manager_adapter,
  Stats::Stat           analysis_stat)
  : StepEventNotifier(),
    AnalysisEventNotifier(),
    StepEventListener(this),
    AnalysisEventListener(this),
    commandLine_(command_line),
    netlistFilename_(command_line.getArgumentValue("netlist")),
    outputManagerAdapter_(output_manager_adapter),
    workingIntgMethod_(0),
    stepErrorControl_(0),
    nonlinearEquationLoader_(0),
    parallelManager_(0),
    dataStore_(0),
    analysisMode_(ANP_MODE_TRANSIENT),
    twoLevelMode_(TWO_LEVEL_MODE_TRANSIENT_DCOP),
    resumeSimulation_(false),
    blockAnalysisFlag_(false),
    daeStateDerivFlag_(true),
    // dcopRestartFlag_(false),
    // saveFlag_(false),
    dotOpSpecified_(false),
    gui_(command_line.argExists("-gui")),
    progressFlag_(true),
    saveTimeGiven_(false),
    savedAlready_(false),
    sensFlag_(false),
    sweepSourceResetFlag_(true),
    switchIntegrator_(false),
    diagnosticMode_(false),
    diagnosticFileNameGiven_(false),
    diagnosticExtremaLimitGiven_(false),
    diagnosticVoltageLimitGiven_(false),
    diagnosticCurrentLimitGiven_(false),
    diagnosticDiscontinuityLimitGiven_(false),
    diagnosticExtremaLimit_(0.0),
    diagnosticVoltageLimit_(0.0),
    diagnosticCurrentLimit_(0.0),
    diagnosticDiscontinuityLimit_(0.0),
    diagnosticFileName_("XyceDiag.out"),
    diagnosticOutputStreamPtr_(NULL),
    xyceTranTimerPtr_(),
    elapsedTimerPtr_(0),
    solverStartTime_(0.0),
    // saveTime_(0.0),
    nextOutputTime_(0.0),
    // nextRestartSaveTime_(0.0),
    analysisObject_(0),
    primaryAnalysisObject_(0),
    analysisStat_(analysis_stat),
    breakPointRestartStep(0)
{
  IO::setTimeIntegratorDebugLevel(command_line, 1);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::~AnalysisManager
//
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
AnalysisManager::~AnalysisManager()
{
  delete nonlinearEquationLoader_;
  delete workingIntgMethod_;
  delete dataStore_;
  delete stepErrorControl_;

  for (std::vector<ProcessorBase *>::iterator it = analysisVector_.begin(), end = analysisVector_.end(); it != end; ++it)
  {
    delete (*it);
  }
  if( diagnosticOutputStreamPtr_ != NULL)
  {
    diagnosticOutputStreamPtr_->close();
    delete diagnosticOutputStreamPtr_;
    diagnosticOutputStreamPtr_ = NULL;
  }
}

Parallel::Machine
AnalysisManager::getComm() const
{
  return parallelManager_->getPDSComm()->comm();
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::notify
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Jul 14 12:15:30 2014
//-----------------------------------------------------------------------------
///
/// Notification that there is a StepEvent.
///
/// @param step_event   information about the event
///
//-----------------------------------------------------------------------------
void AnalysisManager::notify(
  const StepEvent &     step_event)
{
  if (step_event.state_ == StepEvent::STEP_STARTED) {
    dataStore_->setZeroHistory();
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::notify
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Jul 14 12:15:30 2014
//-----------------------------------------------------------------------------
///
/// Notification that there is a AnalysisEvent.
///
/// @param time_integrator_event   information about the event
///
//-----------------------------------------------------------------------------
void AnalysisManager::notify(
  const AnalysisEvent &     analysis_event)
{  
  if( diagnosticMode_ )
  {
    this->OutputDiagnosticInfo(analysis_event);
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::createTimeIntegratorMethod
//
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
//-----------------------------------------------------------------------------
void AnalysisManager::createTimeIntegratorMethod(
  const TimeIntg::TIAParams &   tia_params,
  const unsigned int            integration_method)
{
  workingIntgMethod_->createTimeIntegMethod(integration_method, tia_params, *stepErrorControl_, *dataStore_);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::resetSolverSystem
//
// Purpose       : just like a destructor without the death
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek, Ting Mei
// Creation Date : 7/24/08
//-----------------------------------------------------------------------------
void AnalysisManager::resetSolverSystem()
{
  // dataStore_ is created in initializeAll
  delete dataStore_;
  dataStore_ = 0;

  // stepErrorControl_ is created in initializeAll
  delete stepErrorControl_;
  stepErrorControl_ = 0;

  // workingIntgMethod_ is created in initializeAll
  delete workingIntgMethod_;
  workingIntgMethod_ = 0;

  // Reset step statistics to zero.
  primaryAnalysisObject_->resetAll();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::isPaused
//
// Purpose       : return the true if the simulation is currently paused
//
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical and Microsystems Simulation
// Creation Date : 02/18/2008
//-----------------------------------------------------------------------------
bool AnalysisManager::isPaused() const
{
  return stepErrorControl_->isPauseTime();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::initializeSolverSystem
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Mar 16 14:06:00 2015
//-----------------------------------------------------------------------------
///
/// Initializes the solver system.  Called on initial analysis assembly and on change of time integrator.
///
/// @invariant There be gremlins here.  The datas store, step error control, time integrator (managed by the working
/// integrator method) and nonlinear equation loader are destroyed and recreated.  Hopefully none of the subordinate
/// users of those objects are stll hanging on.  The linear system will have much of its information replaced by data
/// from the newly created datastore and vice versa.
///
///
///
/// @param tia_params           Time integrator parameters to construct new time integrator
/// @param loader               Loader to loader linear system
/// @param linear_system        Linear system
/// @param nonlinear_manager    Nonlinear system
/// @param device_manager       Device Manager
///
/// @return true if initialization was successful
///
//-----------------------------------------------------------------------------
bool AnalysisManager::initializeSolverSystem(
  const TimeIntg::TIAParams &   tia_params,
  Loader::Loader &              loader,
  Linear::System &              linear_system,
  Nonlinear::Manager &          nonlinear_manager,
  Device::DeviceMgr &           device_manager)
{
  Stats::Stat stat(std::string(analysisObject_->getName()), analysisStat_);

  // allocate data store class, which will allocate all the vectors.
  delete dataStore_;
  delete stepErrorControl_;
  delete workingIntgMethod_;
  delete nonlinearEquationLoader_;

  dataStore_ = new TimeIntg::DataStore(tia_params.maxOrder, linear_system.builder());

  workingIntgMethod_ = new TimeIntg::WorkingIntegrationMethod(stat);
  stepErrorControl_ = new TimeIntg::StepErrorControl(netlistFilename_, *this, *workingIntgMethod_, tia_params);
  nonlinearEquationLoader_ = new Loader::NonlinearEquationLoader(*dataStore_, loader, device_manager, *workingIntgMethod_, daeStateDerivFlag_);

  dataStore_->resetAll(tia_params.absErrorTol, tia_params.relErrorTol);

  nextOutputTime_ = stepErrorControl_->initialTime;
  // nextRestartSaveTime_ = stepErrorControl_->initialTime;

  // register the device mask
  linear_system.registerDeviceMaskVector(dataStore_->deviceErrorWeightMask_);

  // Get the RHS and the Jacobian
  dataStore_->JMatrixPtr    = linear_system.getJacobianMatrix();
  dataStore_->RHSVectorPtr  = linear_system.getRHSVector();

  // Get the limiter vectors
  dataStore_->dFdxdVpVectorPtr = linear_system.getdFdxdVpVector();
  dataStore_->dQdxdVpVectorPtr = linear_system.getdQdxdVpVector();

  dataStore_->limiterFlag = loader.getLimiterFlag ();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getBlockAnalysisFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Mar 16 14:12:17 2015
//-----------------------------------------------------------------------------
///
/// Return true if primary analysis is HB or MPDE.
///
///
/// @return true if primary analysis is HB or MPDE.
///
//-----------------------------------------------------------------------------
bool AnalysisManager::getBlockAnalysisFlag() const
{
#ifdef Xyce_STOKHOS_ENABLE
  // ERK: check this
  return 
    dynamic_cast<Analysis::PCE *>(primaryAnalysisObject_) || 
    dynamic_cast<Analysis::EmbeddedSampling *>(primaryAnalysisObject_) || 
    dynamic_cast<Analysis::HB *>(primaryAnalysisObject_) || 
    dynamic_cast<Analysis::MPDE *>(primaryAnalysisObject_);
#else
  return 
    dynamic_cast<Analysis::EmbeddedSampling *>(primaryAnalysisObject_) || 
    dynamic_cast<Analysis::HB *>(primaryAnalysisObject_) || 
    dynamic_cast<Analysis::MPDE *>(primaryAnalysisObject_);
#endif
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::finalExpressionBasedSetup
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 5/4/2021
//-----------------------------------------------------------------------------
void AnalysisManager::finalExpressionBasedSetup()
{
  // This performs the early setup for the outputters for this analysis mode.
  (outputManagerAdapter_.getOutputManager()).earlyPrepareOutput
      (parallelManager_->getPDSComm()->comm(), analysisMode_);

  // now inform the analysis to complete expression setup.  After this, the 
  // maps will be deleted.
  analysisObject_->finalExpressionBasedSetup();

  // ERK "processors" are basically .RESULT statements AFAIK.  Set them up here
  for(int ii=0;ii<processorVector_.size();ii++)
  {
    processorVector_[ii]->finalExpressionBasedSetup();
  }

  Report::safeBarrier(parallelManager_->getPDSComm()->comm());
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::run
// Purpose       : Execute the control loop for the set analysis type.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
///
/// Runs the top level analysis.
///
/// @invariant
///
///
/// @return 
///
//-----------------------------------------------------------------------------
bool AnalysisManager::run()
{
  Stats::StatTop _analysisStat("Analysis");
  Stats::TimeBlock _analysisTimer(_analysisStat);

  if (!primaryAnalysisObject_)
  {
    Report::UserError0() << "No analysis statement in the netlist";
    return false;
  }

  bool runStatus = false;
  {
    // This prepares the outputters for this analysis mode
    IO::ActiveOutput active(outputManagerAdapter_.getOutputManager());

    active.setStepSweepVector(outputManagerAdapter_.getStepSweepVector());
    active.setDCSweepVector(outputManagerAdapter_.getDCSweepVector());
    active.add(parallelManager_->getPDSComm()->comm(), analysisMode_);

    if (VERBOSE_TIME)
      getTIAParams().printParams(Xyce::lout(), nonlinearAnalysisMode(analysisMode_));

    nonlinearEquationLoader_->getLoader().loadDeviceErrorWeightMask((dataStore_->deviceErrorWeightMask_));

    Report::safeBarrier(parallelManager_->getPDSComm()->comm());

    solverStartTime_ = elapsedTimerPtr_->elapsedTime();

    // Start the solvers timers.
    xyceTranTimerPtr_.resetStartTime();

    runStatus = analysisObject_->run();
  }

  return runStatus;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::allocateAnalysisObject
// Purpose       : Allocate analysis objects, and also setup params.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 06/24/10
//-----------------------------------------------------------------------------
///
/// Creates the primary analysis and driving analysis (.STEP, .SAMPLING, dakota).
///
/// If no analysis was specified and .OP was specified, create a DC analysis as the primary analysis
///
/// @invariant
///
//-----------------------------------------------------------------------------
void AnalysisManager::allocateAnalysisObject(AnalysisCreatorRegistry & analysis_registry)
{
  // .OP, but no analysis creates a dcsweep analysis.
  if (dotOpSpecified_ && 
       (analysisCreatorVector_.empty() || (analysisCreatorVector_.size() == 1 && 
             (analysisCreatorVector_.front()->isType<Step>() || 
              analysisCreatorVector_.front()->isType<Sampling>() || 
              analysisCreatorVector_.front()->isType<EmbeddedSampling>()
#ifdef Xyce_STOKHOS_ENABLE
              ||
              analysisCreatorVector_.front()->isType<PCE>()
#endif
              ))))
  {
    CreatorVector::const_iterator it = analysis_registry.begin(); 
    CreatorVector::const_iterator end = analysis_registry.end();
    for ( ; it != end; ++it) 
    {
      if ((*it)->isType<DCSweep>()) 
      {
        analysisMode_ = ANP_MODE_DC_SWEEP;
        primaryAnalysisObject_ = (*it)->create();
        analysisVector_.push_back(primaryAnalysisObject_);
      }
    }
  }

  // CRAZY HACK UNTIL PRIORITY QUEUE OR SOMETHING CLEVER ON ANALYSIS TYPES
  CreatorVector::const_iterator it = analysisCreatorVector_.begin(); 
  CreatorVector::const_iterator end = analysisCreatorVector_.end();
  for ( ; it != end; ++it) 
  {
    if (!(*it)->isType<Step>() 
        && !(*it)->isType<Sampling>()
        && !(*it)->isType<EmbeddedSampling>()
#ifdef Xyce_STOKHOS_ENABLE
        && !(*it)->isType<PCE>()
#endif
        ) 
    {
      primaryAnalysisObject_ = (*it)->create();
      analysisVector_.push_back(primaryAnalysisObject_);
    }
  }

  it = analysisCreatorVector_.begin(); 
  end = analysisCreatorVector_.end();
  for ( ; it != end; ++it) 
  {
    if ((*it)->isType<Step>()) 
    {
      analysisObject_ = (*it)->create();
      analysisVector_.push_back(analysisObject_);
      pushActiveAnalysis(analysisObject_);
    }

    if ((*it)->isType<Sampling>()) 
    {
      analysisObject_ = (*it)->create();
      analysisVector_.push_back(analysisObject_);
      pushActiveAnalysis(analysisObject_);
    }

    if ((*it)->isType<EmbeddedSampling>()) 
    {
      analysisObject_ = (*it)->create();
      analysisVector_.push_back(analysisObject_);
      pushActiveAnalysis(analysisObject_);
    }

#ifdef Xyce_STOKHOS_ENABLE
    if ((*it)->isType<PCE>()) 
    {
      analysisObject_ = (*it)->create();
      analysisVector_.push_back(analysisObject_);
      pushActiveAnalysis(analysisObject_);
    }
#endif

#ifdef Xyce_ROL
    if ((*it)->isType<ROL>()) 
    {
      primaryAnalysisObject_ = (*it)->create(); // TT: We want ROL analysis to by primary analysis (unlike STEP), hence call it primaryAnalysisObject
      analysisVector_.push_back(primaryAnalysisObject_);
      //pushActiveAnalysis(primaryAnalysisObject_);
    }
#endif
  }

  // ERK "processors" are basically .RESULT statements AFAIK.
  // I have no idea why someone thought it was a good idea to add
  // them to the analysis vector, or to have processors and
  // analyses be derived from the same base class, but that is
  // what happens.
  for(CreatorSet::const_iterator it = processorCreatorSet_.begin(),
      end = processorCreatorSet_.end(); it != end; ++it)
  {
    ProcessorBase * procPtr  = (*it)->create();
    processorVector_.push_back(procPtr);
    analysisVector_.push_back(procPtr);
  }

  if (!primaryAnalysisObject_)
  {
    Report::UserError() << "Analysis mode " << analysisMode_ << " is not available";
    return;
  }

  pushActiveAnalysis(primaryAnalysisObject_);

  if (!analysisObject_)
  {
    analysisObject_ = primaryAnalysisObject_;
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::printLoopInfo
// Purpose       : Prints out time loop information.
// Special Notes : Prints stats from save point start to save point finish.
//                 Special case 0,0 is entire run to this point
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/26/00
//-----------------------------------------------------------------------------
bool AnalysisManager::printLoopInfo(int start, int finish)
{
  return primaryAnalysisObject_->printLoopInfo(start, finish);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setOPAnalysisParams
// Purpose       : Handle OP statement. (.OP)
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/15/07
//-----------------------------------------------------------------------------
bool AnalysisManager::setOPAnalysisParams(
  const Util::OptionBlock &     paramsBlock)
{
  dotOpSpecified_ = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setSensOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, Sandia
// Creation Date : 6/5/13
//-----------------------------------------------------------------------------
bool AnalysisManager::setSensOptions(const Util::OptionBlock & OB)
{
  sensFlag_ = true;
  return true;
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setDiagnosticMode
// Purpose       : Set the Analysis Manager to gather diagnostic info for 
//                 user circuit debugging.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Sandia
// Creation Date : 11/17/22
//-----------------------------------------------------------------------------
bool AnalysisManager::setDiagnosticMode(const Util::OptionBlock & OB)
{
  diagnosticMode_ = true;
  bool result = false;
  bool subResult = false;
  for (Util::ParamList::const_iterator it = OB.begin(), end = OB.end(); it != end; ++it)
  {
    const Util::Param &param = *it;
    // reset flag on if a parameter was set each time through the loop
    subResult = false; 

    if (param.uTag() == "DIAGFILENAME")
    {
      diagnosticFileName_ = it->stringValue();
      diagnosticFileNameGiven_=true;
      subResult = true;
    }

    if (!subResult)
    {
      bool subResult = Util::setValue( param, "EXTREMALIMIT", diagnosticExtremaLimit_, diagnosticExtremaLimitGiven_)
        || Util::setValue( param, "VOLTAGELIMIT", diagnosticVoltageLimit_, diagnosticVoltageLimitGiven_)
        || Util::setValue( param, "CURRENTLIMIT", diagnosticCurrentLimit_, diagnosticCurrentLimitGiven_)
        || Util::setValue( param, "DISCLIMIT", diagnosticDiscontinuityLimit_, diagnosticDiscontinuityLimitGiven_);
    }

    // need to or in successive results of parameter setting 
    result = result || subResult;
  }
  
  // if a diagnostic file name wasn't given then make a default name from the netlist file name.
  if( !diagnosticFileNameGiven_)
  {
    diagnosticFileName_ = netlistFilename_ + ".dia";
  }
  return result;  
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::completeOPStartStep
// Purpose       : Call to rotate next state to current state following a
//               : constrained DCOP solve when using a previous operating point
// Special Notes : Called from the Xyce-LOCA interface.
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 06/27/2006
//-----------------------------------------------------------------------------
bool AnalysisManager::completeOPStartStep( )
{
  bool bsuccess = true;

  bsuccess = dataStore_->updateStateDataArrays ();

  dataStore_->setConstantHistory();
  dataStore_->equateTmpVectors();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::completeHomotopyStep
// Purpose       : Call to do final cleanup after a single,
//                 successful homotopy step.
// Special Notes : Called from the Xyce-LOCA interface.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/20/2006
//-----------------------------------------------------------------------------
bool AnalysisManager::completeHomotopyStep(
  Loader::NonlinearEquationLoader &     loader,
  const std::vector<std::string> &      paramNames,
  const std::vector<double> &           paramVals,
  Linear::Vector *                      solnVecPtr )
{
  bool bsuccess = true;

  if (DEBUG_ANALYSIS)
    Xyce::dout() << "\n " << netlistFilename_
                 << " AnalysisManager::completeHomotopyStep " << std::endl;

  // Rotate the data vectors:
  bool bs1 = dataStore_->updateStateDataArrays ();
  bsuccess = bsuccess && bs1;
  dataStore_->setConstantHistory();
  dataStore_->equateTmpVectors();

  // Pass info in to the lower level solver:
  loader.homotopyStepSuccess(paramNames,paramVals);

  // Call output
  outputManagerAdapter_.outputHomotopy( paramNames, paramVals, *solnVecPtr );

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::failHomotopyStep
// Purpose       : Call to do final cleanup after a single,
//                 failed homotopy step.
// Special Notes : Called from the Xyce-LOCA interface.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/30/2006
//-----------------------------------------------------------------------------
bool AnalysisManager::failHomotopyStep(Loader::NonlinearEquationLoader &loader)
{
  if (DEBUG_ANALYSIS)
    Xyce::dout() << "\n " << netlistFilename_
                 << " AnalysisManager::failHomotopyStep " << std::endl;

// Pass info in to the lower level solver:
  loader.homotopyStepFailure();

  return true;
}

// ***** Accessor methods *****
//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setBeginningIntegrationFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void AnalysisManager::setBeginningIntegrationFlag(bool bif)
{
  currentAnalysisStack_.back()->setBeginningIntegrationFlag(bif);
  primaryAnalysisObject_->setBeginningIntegrationFlag(bif);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getBeginningIntegrationFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AnalysisManager::getBeginningIntegrationFlag() const
{
  return currentAnalysisStack_.back()->getBeginningIntegrationFlag();
//  return primaryAnalysisObject_->getBeginningIntegrationFlag();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setIntegrationMethod
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void AnalysisManager::setIntegrationMethod(int im)
{
  currentAnalysisStack_.back()->setIntegrationMethod(im);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getIntegrationMethod
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
int AnalysisManager::getIntegrationMethod ()
{
  return currentAnalysisStack_.back()->getIntegrationMethod();
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTotalLinearSolutionTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
double AnalysisManager::getTotalLinearSolutionTime() const
{
  return primaryAnalysisObject_->getTotalLinearSolutionTime();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTotalResidualLoadTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
double AnalysisManager::getTotalResidualLoadTime() const
{
  return primaryAnalysisObject_->getTotalResidualLoadTime();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTotalJacobianLoadTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
double AnalysisManager::getTotalJacobianLoadTime() const
{
  return primaryAnalysisObject_->getTotalJacobianLoadTime();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getDoubleDCOPEnabled ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool AnalysisManager::getDoubleDCOPEnabled() const
{
  return primaryAnalysisObject_->getDoubleDCOPEnabled ();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::isSimulationComplete
// Purpose       : return boolean signifying whether simulation complete or
//                 not.
//
// Special Notes :THIS VERSION IS ONLY VALID FOR TRANSIENT RUNS, where
//                 completion of the simulation means integration to final
//                 time.
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 05/06/2004
//-----------------------------------------------------------------------------
bool AnalysisManager::isSimulationComplete()
{
  if (analysisMode_ == ANP_MODE_TRANSIENT)
  {
    return stepErrorControl_->isFinished();
  }
  else
  {
    Report::DevelFatal0().in("AnalysisManager::simulationComplete") 
      << "Called for non-transient run, not currently valid";
    return false;
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getPauseTime
//
// Purpose       : return the time at which the simulation will pause
//
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/29/2004
//-----------------------------------------------------------------------------
double AnalysisManager::getPauseTime() const
{
  return stepErrorControl_->pauseTime;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getStepNumber
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
unsigned int AnalysisManager::getStepNumber() const
{
  if (!currentAnalysisStack_.empty()){
    return currentAnalysisStack_.back()->getStepNumber();
  }

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setStepNumber
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
void AnalysisManager::setStepNumber(int step)
{
  if (primaryAnalysisObject_){
       primaryAnalysisObject_->setStepNumber(step);
  }

}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTranStepNumber
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
int AnalysisManager::getTranStepNumber()
{
  int number = 0;

  if (analysisMode_ == ANP_MODE_TRANSIENT && primaryAnalysisObject_)
    number = primaryAnalysisObject_->getTranStepNumber();

  return number;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setTranStepNumber
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
void AnalysisManager::setTranStepNumber (int step)
{
  if (primaryAnalysisObject_)
    primaryAnalysisObject_->setTranStepNumber(step);
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getInitTranFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
bool AnalysisManager::getInitTranFlag() const
{
  return getStepNumber() <= 0;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTime
// Purpose       : Gets the next time value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
double AnalysisManager::getTime() const
{
  return stepErrorControl_->nextTime;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getFinalTime
// Purpose       : Gets the final time-step value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
double AnalysisManager::getFinalTime() const
{
  return stepErrorControl_->finalTime;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getFinalTimeForRemeasure
// Purpose       : Gets the final time-step value.
// Special Notes : This version should ONLY be used for code related to -remeasure
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 8/13/2020
//-----------------------------------------------------------------------------
double AnalysisManager::getFinalTimeForRemeasure() const
{
  return getTIAParams().finalTime;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getInitialTime
// Purpose       : Gets the initial time-step value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
double AnalysisManager::getInitialTime() const
{
  return stepErrorControl_->initialTime;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getDCOPFlag
// Purpose       : Gets a flag indicating we are in steady state.
//                  (steady=true)
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/02/01
//-----------------------------------------------------------------------------
bool AnalysisManager::getDCOPFlag() const
{
  return primaryAnalysisObject_ ? primaryAnalysisObject_->getDCOPFlag() : false;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTranOPFlag
// Purpose       : Gets a flag indicating we are in a DCOP calculation that is
//                 the initialization for a transient simulation.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/28/05
//-----------------------------------------------------------------------------
bool AnalysisManager::getTranOPFlag() const
{
  return ((analysisMode_ == ANP_MODE_TRANSIENT || primaryAnalysisObject_->isAnalysis(ANP_MODE_TRANSIENT))
          && (primaryAnalysisObject_->getIntegrationMethod()) == TimeIntg::methodsEnum::NO_TIME_INTEGRATION);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getACOPFlag
// Purpose       : Gets a flag indicating we are in a DCOP calculation that is
//                 the initialization for a AC simulation.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/16/12
//-----------------------------------------------------------------------------
bool AnalysisManager::getACOPFlag() const
{
  bool return1 = (analysisMode_ == ANP_MODE_AC || primaryAnalysisObject_->isAnalysis(ANP_MODE_AC))
      && (primaryAnalysisObject_->getIntegrationMethod()) == TimeIntg::methodsEnum::NO_TIME_INTEGRATION;

  // if the analysis type is noise, this flag must also be set, as it is based on an 
  // AC-style calculation.
  bool return2 = (analysisMode_ == ANP_MODE_NOISE || primaryAnalysisObject_->isAnalysis(ANP_MODE_NOISE))
      && (primaryAnalysisObject_->getIntegrationMethod()) == TimeIntg::methodsEnum::NO_TIME_INTEGRATION;

  return (return1 || return2);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getACFlag
// Purpose       : Gets a flag indicating we are in an AC calculation
// Special Notes : This function needs to check whether the primary analysis
//                 object exists, in order to balance the competing needs of
//                 making I() invalid, in some cases, for .AC and .NOISE (see SON
//                 Bug 855) and allowing for -remeasure of lead currents to work
//                 for .TRAN and .DC.  The primary analysis object is not made
//                 during -remeasure.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 04/08/2019
//-----------------------------------------------------------------------------
bool AnalysisManager::getACFlag() const
{
  if (primaryAnalysisObject_ != 0)
    return(analysisMode_ == ANP_MODE_AC || primaryAnalysisObject_->isAnalysis(ANP_MODE_AC));
  else
    return(analysisMode_ == ANP_MODE_AC);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getNoiseFlag
// Purpose       : Gets a flag indicating we are in a noise calculation
// Special Notes : This function needs to check whether the primary analysis
//                 object exists, in order to balance the competing needs of
//                 making I() invalid, in some cases, for .AC and .NOISE (see SON
//                 Bug 855) and allowing for -remeasure of lead currents to work
//                 for .TRAN and .DC.  The primary analysis object is not made
//                 during -remeasure.
// Scope         : public
// Creator       : Tom Russo, SNL, Parallel Computational Sciences
// Creation Date : 02/19/2015
//-----------------------------------------------------------------------------
bool AnalysisManager::getNoiseFlag() const
{
  if (primaryAnalysisObject_ != 0)
    return(analysisMode_ == ANP_MODE_NOISE || primaryAnalysisObject_->isAnalysis(ANP_MODE_NOISE));
  else
    return(analysisMode_ == ANP_MODE_NOISE);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getDCSweepFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/28/05
//-----------------------------------------------------------------------------
bool AnalysisManager::getDCSweepFlag() const
{
  return ((analysisMode_ == ANP_MODE_DC_SWEEP || primaryAnalysisObject_->isAnalysis(ANP_MODE_DC_SWEEP)));
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTransientFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/28/05
//-----------------------------------------------------------------------------
bool AnalysisManager::getTransientFlag () const
{
  return (((analysisMode_ == ANP_MODE_TRANSIENT || primaryAnalysisObject_->isAnalysis(ANP_MODE_TRANSIENT))
     && (primaryAnalysisObject_->getIntegrationMethod()) != TimeIntg::methodsEnum::NO_TIME_INTEGRATION));
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getDoubleDCOPStep
// Purpose       : Gets the double DC Operating Point step value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/25/01
//-----------------------------------------------------------------------------
int AnalysisManager::getDoubleDCOPStep() const
{
  return primaryAnalysisObject_->getDoubleDCOPStep();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::updateDerivs
// Purpose       : Calls the time  int. method to update the corrector
//                 derivatives.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
bool AnalysisManager::updateDerivs()
{
  workingIntgMethod_->obtainCorrectorDeriv();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerElapsedTimer
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 02/24/06
//-----------------------------------------------------------------------------
bool AnalysisManager::registerElapsedTimer(Util::Timer * et)
{
  elapsedTimerPtr_ = et;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerExpressionGroup
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/19/2020
//-----------------------------------------------------------------------------
bool AnalysisManager::registerExpressionGroup(Teuchos::RCP<Xyce::Util::baseExpressionGroup> & group)
{
  expressionGroup_ = group;
  return (!(Teuchos::is_null(expressionGroup_)));
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getRestartDataSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
int AnalysisManager::getRestartDataSize( bool pack ) const
{
  return stepErrorControl_->getRestartDataSize( pack );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::dumpRestartData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::dumpRestartData(
  char *                        buf,
  int                           bsize,
  int &                         pos,
  Parallel::Communicator *      comm,
  bool                          pack)
{
  return stepErrorControl_->dumpRestartData(buf, bsize, pos, comm, pack);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::restoreRestartData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::restoreRestartData(
  char *                        buf,
  int                           bsize,
  int &                         pos,
  Parallel::Communicator *      comm,
  bool                          pack)
{
  bool ret = stepErrorControl_->restoreRestartData(
      buf, bsize, pos, comm, pack, getTIAParams().initialTime);

  // Create the time integration method before continuting to restore restart data.
  // NOTE: It is important that the time integrator exists before we restore the data.
  if (ret)
  {
    createTimeIntegratorMethod(getTIAParams(), getIntegrationMethod()); 
  }

  return ret;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getSolnVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::getSolnVarData( const int & gid,
                                             std::vector<double> & varData ) const
{
  return workingIntgMethod_->getSolnVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getStateVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::getStateVarData(
  const int &           gid,
  std::vector<double> & varData ) const
{
  return workingIntgMethod_->getStateVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getStoreVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool AnalysisManager::getStoreVarData( const int & gid,
                                              std::vector<double> & varData ) const
{
  return workingIntgMethod_->getStoreVarData( gid, varData );
}
//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setSolnVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::setSolnVarData( const int & gid,
                                             const std::vector<double> & varData )
{
  return workingIntgMethod_->setSolnVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setStateVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::setStateVarData( const int & gid,
                                              const std::vector<double> & varData )
{
  return workingIntgMethod_->setStateVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setStoreVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool AnalysisManager::setStoreVarData( const int & gid,
                                              const std::vector<double> & varData )
{
  return workingIntgMethod_->setStoreVarData( gid, varData );
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTIAParams const
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
const TimeIntg::TIAParams & AnalysisManager::getTIAParams() const
{
  return primaryAnalysisObject_->getTIAParams();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTIAParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
TimeIntg::TIAParams & AnalysisManager::getTIAParams()
{
  return primaryAnalysisObject_->getTIAParams();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerParallelServices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/19/02
//-----------------------------------------------------------------------------
bool AnalysisManager::registerParallelServices(Parallel::Manager * pds)
{
  parallelManager_ = pds;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setNextSolVectorPtr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,SNL, Computational Sciences
// Creation Date : 04/22/2003
//-----------------------------------------------------------------------------
bool AnalysisManager::setNextSolVectorPtr (Linear::Vector * solVecPtr)
{
  return dataStore_->setNextSolVectorPtr (solVecPtr);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setPauseTime
//
// Purpose       : Set the time at which to pause the simulation
//
// Special Notes : This is a temporary implementation that I'm using to
//                 begin the process of hiding time integrator internals from
//                 the "N_CIR_Xyce::simulateUntil" method, so that I can
//                 ultimately change those internals without the simulateUntil
//                 method
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/29/2004
//-----------------------------------------------------------------------------
void AnalysisManager::setPauseTime(
  double        pause_time,
  double        initial_time)
{
  stepErrorControl_->setBreakPoint(
      Util::BreakPoint(pause_time, Util::BreakPoint::PAUSE), initial_time);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getDCOPSolve
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
DCOPType AnalysisManager::getDCOPSolve() const
{
  if (!getDoubleDCOPEnabled())
    return OFF;
  else if (getAnalysisObject().getIntegrationMethod() == 0 && getDoubleDCOPStep() == 0)
    return NL_POISSON;
  else
    return DRIFT_DIFFUSION;
  
  // if (!getDoubleDCOPEnabled())
  //     return OFF;
  //   else if (getDoubleDCOPStep() == 0)
  //     return NL_POISSON;
  //   else
  //     return DRIFT_DIFFUSION;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getCurrentFreq
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/7/2018
//-----------------------------------------------------------------------------
double AnalysisManager::getCurrentFreq() const
{
  return analysisObject_->getCurrentFreq();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setRFParamsRequested
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/31/2019
//-----------------------------------------------------------------------------
void AnalysisManager::setRFParamsRequested(const std::string & type)
{
  analysisObject_->setRFParamsRequested(type);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getNodeNameFromIndex
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 11/21/2022
//-----------------------------------------------------------------------------
std::string AnalysisManager::getNodeNameFromIndex( const int varIndex ) const
{
  std::string node_name = "N/A";
  // This is a bit messy as all of the primary analysis classes 
  // derive off of AnalysisBase.  But only the derived classes have a reference 
  // to Topology which is needed to look up a node name.  We can't just make a 
  // Topology Class object to hold the reference we find in the primaryAnalysis 
  // object because the Topology class doesn't have a simple default constructor.
  // It might make sense to move Topology to the AnalysisBase class, but some 
  // classes derived form AnalysisBase, like STEP, don't have or need Topology.
  
  const Topo::Topology * topologyPtr = NULL;
  
  if(getTransientFlag())
  {
    Transient * transientAnalysisObj = dynamic_cast<Transient *>(primaryAnalysisObject_);
    if ( transientAnalysisObj != NULL)
    {
      topologyPtr = &(transientAnalysisObj->getTopology());
    }
  }
  
  if(getDCSweepFlag())
  {
    DCSweep * dcAnalysisObj = dynamic_cast<DCSweep *>(primaryAnalysisObject_);
    if ( dcAnalysisObj != NULL)
    {
      topologyPtr = &(dcAnalysisObj->getTopology());
    }
  }
  
  if(getACFlag())
  {
    AC * acAnalysisObj = dynamic_cast<AC *>(primaryAnalysisObject_);
    if ( acAnalysisObj != NULL)
    {
      topologyPtr = &(acAnalysisObj->getTopology());
    }
  }
  
  // if we have a topology object then try to look up the name
  if ( topologyPtr != NULL)
  {
    const std::vector<const std::string *> name_vec = topologyPtr->getSolutionNodeNames();
    
    Xyce::Parallel::Communicator& pdsComm = *(this->getPDSManager()->getPDSComm());

    // Get the node name, communicate it in parallel.
    if ( pdsComm.isSerial() )
    {
      if ((varIndex > -1) && (varIndex < name_vec.size()))
        node_name = *name_vec[varIndex];
    }
    else
    {
      // Convert outIndex from GID to LID.  If this processor does not own that solution id, it will be -1.
      // Also, if this failure is being printed for a block analysis type, there isn't a
      // clear way to associate the LID with the node name (HB, PCE, ES, will map this differently)
      //Teuchos::RCP<Parallel::ParMap> solnMap = (dataStore_->builder_).getSolutionMap();
      int varIndex_LID = ((dataStore_->builder_).getSolutionMap())->globalToLocalIndex( varIndex ); 

      // Now determine which processor owns this solution node
      int proc, tmp_proc = -1;
      if ((varIndex_LID > -1) && (varIndex_LID < name_vec.size()))
      {
        tmp_proc = Parallel::rank(parallelManager_->getPDSComm()->comm());
        node_name = *name_vec[varIndex_LID];
      }
      pdsComm.maxAll( &tmp_proc, &proc, 1 );

      // Make sure someone owns this GID, otherwise proc will be -1.
      if (proc != -1)
      {
        // Communicate the node name length and character string from the processor that owns it
        int length = node_name.length();
        pdsComm.bcast( &length, 1, proc ); 
        char *buffer = (char *)malloc( length+1 );

        if (varIndex_LID > -1)
          strcpy( buffer, node_name.c_str() );

        pdsComm.bcast( buffer, length+1, proc );
        node_name = std::string( buffer );
        free( buffer );
      }
    }
  }
  return node_name; 
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getNodeNameFromLocalIndex
// Purpose       :
// Special Notes : No parallel communication used in this routine.
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 11/21/2022
//-----------------------------------------------------------------------------
std::string AnalysisManager::getNodeNameFromLocalIndex( const int varIndex ) const
{
  std::string node_name = "N/A";
  // This is a bit messy as all of the primary analysis classes 
  // derive off of AnalysisBase.  But only the derived classes have a reference 
  // to Topology which is needed to look up a node name.  We can't just make a 
  // Topology Class object to hold the reference we find in the primaryAnalysis 
  // object because the Topology class doesn't have a simple default constructor.
  // It might make sense to move Topology to the AnalysisBase class, but some 
  // classes derived form AnalysisBase, like STEP, don't have or need Topology.
  
  const Topo::Topology * topologyPtr = NULL;
  
  if(getTransientFlag())
  {
    Transient * transientAnalysisObj = dynamic_cast<Transient *>(primaryAnalysisObject_);
    if ( transientAnalysisObj != NULL)
    {
      topologyPtr = &(transientAnalysisObj->getTopology());
    }
  }
  
  if(getDCSweepFlag())
  {
    DCSweep * dcAnalysisObj = dynamic_cast<DCSweep *>(primaryAnalysisObject_);
    if ( dcAnalysisObj != NULL)
    {
      topologyPtr = &(dcAnalysisObj->getTopology());
    }
  }
  
  if(getACFlag())
  {
    AC * acAnalysisObj = dynamic_cast<AC *>(primaryAnalysisObject_);
    if ( acAnalysisObj != NULL)
    {
      topologyPtr = &(acAnalysisObj->getTopology());
    }
  }
  
  // if we have a topology object then try to look up the name
  if(topologyPtr != NULL)
  {
    const std::vector<const std::string *> name_vec = topologyPtr->getSolutionNodeNames();
    if ((varIndex > -1) && (varIndex < name_vec.size()))
    {
        node_name = *name_vec[varIndex];
    }
  }
  
  return node_name; 
}



//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getNodeTypeFromIndex
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 11/21/2022
//-----------------------------------------------------------------------------
char AnalysisManager::getNodeTypeFromIndex( const int varIndex ) const
{
  char node_type = '\0';
  char node_type_found = '\0';
  int nodeAsInt = 0;
  int nodeTypeAsInt = 0;

  // This is a bit messy as all of the primary analysis classes 
  // derive off of AnalysisBase.  But only the derived classes have a reference 
  // to Topology which is needed to look up a node name.  We can't just make a 
  // Topology Class object to hold the reference we find in the primaryAnalysis 
  // object because the Topology class doesn't have a simple default constructor.
  // It might make sense to move Topology to the AnalysisBase class, but some 
  // classes derived form AnalysisBase, like STEP, don't have or need Topology.
  
  const Topo::Topology * topologyPtr = NULL;
  
  if(getTransientFlag())
  {
    Transient * transientAnalysisObj = dynamic_cast<Transient *>(primaryAnalysisObject_);
    if ( transientAnalysisObj != NULL)
    {
      topologyPtr = &(transientAnalysisObj->getTopology());
    }
  }
  
  if(getDCSweepFlag())
  {
    DCSweep * dcAnalysisObj = dynamic_cast<DCSweep *>(primaryAnalysisObject_);
    if ( dcAnalysisObj != NULL)
    {
      topologyPtr = &(dcAnalysisObj->getTopology());
    }
  }
  
  if(getACFlag())
  {
    AC * acAnalysisObj = dynamic_cast<AC *>(primaryAnalysisObject_);
    if ( acAnalysisObj != NULL)
    {
      topologyPtr = &(acAnalysisObj->getTopology());
    }
  }
  
  // if we have a topology object then try to look up the name
  if ( topologyPtr != NULL)
  {
    const std::vector<char> type_vec = topologyPtr->getVarTypes(); 
    const std::vector<const std::string *> name_vec = topologyPtr->getSolutionNodeNames();
    
    Xyce::Parallel::Communicator& pdsComm = *(this->getPDSManager()->getPDSComm());

    // Get the node name, communicate it in parallel.
    if ( pdsComm.isSerial() )
    {
      if ((varIndex > -1) && (varIndex < type_vec.size()))
        node_type = type_vec[varIndex];
    }
    else
    {
      // Convert outIndex from GID to LID.  If this processor does not own that solution id, it will be -1.
      // Also, if this failure is being printed for a block analysis type, there isn't a
      // clear way to associate the LID with the node name (HB, PCE, ES, will map this differently)
      //Teuchos::RCP<Parallel::ParMap> solnMap = (dataStore_->builder_).getSolutionMap();
      int varIndex_LID = ((dataStore_->builder_).getSolutionMap())->globalToLocalIndex( varIndex ); 

      // Now determine which processor owns this solution node
      int proc = 0;
      int tmp_proc = -1;
      tmp_proc = Parallel::rank(parallelManager_->getPDSComm()->comm());
      if ((varIndex_LID > -1) && (varIndex_LID < type_vec.size()))
      {
        node_type = type_vec[varIndex_LID];
      }
      //pdsComm.maxAll( &tmp_proc, &proc, 1 );
      nodeTypeAsInt = (int) node_type;
      pdsComm.maxAll( &nodeTypeAsInt, &nodeAsInt, 1 );
      node_type_found = (char) nodeAsInt;
      node_type = node_type_found;
    }
  }
  return node_type; 
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getNodeTypeFromLocalIndex
// Purpose       :
// Special Notes : No parallel communication if index is not local
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 11/21/2022
//-----------------------------------------------------------------------------
char AnalysisManager::getNodeTypeFromLocalIndex( const int varIndex ) const
{
  char node_type = '\0';
  // This is a bit messy as all of the primary analysis classes 
  // derive off of AnalysisBase.  But only the derived classes have a reference 
  // to Topology which is needed to look up a node name.  We can't just make a 
  // Topology Class object to hold the reference we find in the primaryAnalysis 
  // object because the Topology class doesn't have a simple default constructor.
  // It might make sense to move Topology to the AnalysisBase class, but some 
  // classes derived form AnalysisBase, like STEP, don't have or need Topology.
  
  const Topo::Topology * topologyPtr = NULL;
  
  if(getTransientFlag())
  {
    Transient * transientAnalysisObj = dynamic_cast<Transient *>(primaryAnalysisObject_);
    if ( transientAnalysisObj != NULL)
    {
      topologyPtr = &(transientAnalysisObj->getTopology());
    }
  }
  
  if(getDCSweepFlag())
  {
    DCSweep * dcAnalysisObj = dynamic_cast<DCSweep *>(primaryAnalysisObject_);
    if ( dcAnalysisObj != NULL)
    {
      topologyPtr = &(dcAnalysisObj->getTopology());
    }
  }
  
  if(getACFlag())
  {
    AC * acAnalysisObj = dynamic_cast<AC *>(primaryAnalysisObject_);
    if ( acAnalysisObj != NULL)
    {
      topologyPtr = &(acAnalysisObj->getTopology());
    }
  }
  Transient * transientAnalysisObj = dynamic_cast<Transient *>(primaryAnalysisObject_);
  if ( topologyPtr != NULL)
  {
    const std::vector<char> type_vec =topologyPtr->getVarTypes(); 
    if ((varIndex > -1) && (varIndex < type_vec.size()))
    {
        node_type = type_vec[varIndex];
    }
  }
  return node_type; 
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::OutputDiagnosticInfo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 11/21/2022
//-----------------------------------------------------------------------------

void AnalysisManager::OutputDiagnosticInfo(const AnalysisEvent & analysis_event)
{
  bool outputHeaderLine = true;
  Xyce::Parallel::Communicator& pdsComm = *(this->getPDSManager()->getPDSComm());
  
  if((diagnosticOutputStreamPtr_ == NULL) && (pdsComm.isSerial() || (pdsComm.procID() == 0)))
  {
    //open the output file
    diagnosticOutputStreamPtr_ = new std::ofstream();
    diagnosticOutputStreamPtr_->open(diagnosticFileName_);
    (*diagnosticOutputStreamPtr_) << "Xyce Diagnostic Output" << std::endl;
    if(diagnosticExtremaLimitGiven_)
    {
      (*diagnosticOutputStreamPtr_) 
        << " Solution Extrema output requested with extremalimit = " << diagnosticExtremaLimit_ << std::endl;
    }
    if(diagnosticVoltageLimitGiven_)
    {
      (*diagnosticOutputStreamPtr_) 
        << " Node voltage output requested with voltagelimit = " << diagnosticVoltageLimit_ << std::endl;
    }
    if(diagnosticCurrentLimitGiven_)
    {
      (*diagnosticOutputStreamPtr_) 
        << " Lead current output requested with currentlimit = " << diagnosticCurrentLimit_ << std::endl;
    }
    if(diagnosticDiscontinuityLimitGiven_)
    {
      (*diagnosticOutputStreamPtr_) 
        << " Discontinuity output requested with disclimit = " << diagnosticDiscontinuityLimit_ << std::endl;
    }
  }
  
  // on step==0 for DC or TRAN there can be useful, general information about
  // solver status and solver method used like gmin stepping.  Output this general 
  // info only when step==0
  //
  // Look for DC_OP_GMIN_STEPPING_FAILED and DC_OP_SOURCE_STEPPING_FAILED 
  // as these are useful messages to pass on to the user. 
  //
  // only output if the stream is open
  if(diagnosticOutputStreamPtr_ != NULL)
  {
    if(analysis_event.step_ == 0) 
    {
      (*diagnosticOutputStreamPtr_) << "Analysis event " << analysis_event.state_ << " "
        << analysis_event.outputType_ <<  std::endl;
    }
    else if(analysis_event.state_ == AnalysisEvent::DC_OP_GMIN_STEPPING_FAILED)
    {
      (*diagnosticOutputStreamPtr_) << "Analysis event  " << analysis_event.state_ << " "
        << analysis_event.outputType_ << " smallest gmin = " << std::scientific << analysis_event.step_ << std::defaultfloat <<  std::endl;
    }
    else if(analysis_event.state_ == AnalysisEvent::DC_OP_SOURCE_STEPPING_FAILED)
    {
      (*diagnosticOutputStreamPtr_) << "Analysis event  " << analysis_event.state_ << " "
        << analysis_event.outputType_ << " source value = " << analysis_event.step_ <<  std::endl;
    }
  }
  
  
  if( diagnosticExtremaLimitGiven_)
  {
    // get largest absolute value from solution vector
    double value = 0.0;
    int localId = 0;
    (this->getDataStore())->currSolutionPtr->infNorm( &value, &localId );      
    if( fabs(value) > diagnosticExtremaLimit_)
    {
      if( (diagnosticOutputStreamPtr_ != NULL) &&  outputHeaderLine )
      {
        (*diagnosticOutputStreamPtr_) << " Extreme value found in " 
          << analysis_event.outputType_ 
          << " analysis at " 
          << analysis_event.state_;
        if (this->getAnalysisMode() == ANP_MODE_TRANSIENT)
        {
          (*diagnosticOutputStreamPtr_) << " time=" << this->getTime() << std::endl;
        }
        else
        {
          (*diagnosticOutputStreamPtr_) << " Step=" << analysis_event.step_ << std::endl;
        }
        outputHeaderLine = false;
      }
      std::string nodeName = this->getNodeNameFromIndex( localId ); 
      char varType = this->getNodeTypeFromIndex( localId ); 
    
      if(diagnosticOutputStreamPtr_ != NULL) 
      {
        (*diagnosticOutputStreamPtr_) 
          << "     " << varType << "(" << nodeName << ")=" << value << std::endl;
      }
    }
    
  }
  
  
  
  if( diagnosticVoltageLimitGiven_ || diagnosticCurrentLimitGiven_)
  {
    // stream buffers for the output to keep it organized.
    std::stringstream voltageOutput;
    std::stringstream currentOutput;
    voltageOutput.flush();
    currentOutput.flush();
    bool outputVoltageHeaderLine = true;
    bool outputCurrentHeaderLine = true;
    
    // loop over the full solution vector 
    double value = 0.0;
    //auto numSolVars = (this->getDataStore())->solutionSize;
    auto numSolVars = (this->getDataStore())->currSolutionPtr->localLength();
    for( auto localID=0 ; localID < numSolVars; localID++)
    {
      value =(*(this->getDataStore()->currSolutionPtr))[localID];
      auto globalID = ((dataStore_->builder_).getSolutionMap())->localToGlobalIndex( localID );
      char varType = this->getNodeTypeFromLocalIndex( localID ); 
      
      if( diagnosticVoltageLimitGiven_ && (varType=='V') && (fabs(value) > diagnosticVoltageLimit_))
      {
        if( outputVoltageHeaderLine )
        {
          voltageOutput << " Voltage over limit value found in " 
            << analysis_event.outputType_ 
            << " analysis at " 
            << analysis_event.state_;
          if (this->getAnalysisMode() == ANP_MODE_TRANSIENT)
          {
            voltageOutput << " time=" << this->getTime() << std::endl;
          }
          else
          {
            voltageOutput << " Step=" << analysis_event.step_ << std::endl;
          }
          outputVoltageHeaderLine = false;
        }
        std::string nodeName = this->getNodeNameFromLocalIndex( localID ); 
        voltageOutput
          << "     V" << "(" << nodeName << ")=" << value << std::endl;
      }
      if( diagnosticCurrentLimitGiven_ && (varType=='I') && (fabs(value) > diagnosticCurrentLimit_))
      {
        if( outputCurrentHeaderLine )
        {
          currentOutput << " Current over limit value found in " 
            << analysis_event.outputType_ 
            << " analysis at " 
            << analysis_event.state_;
          if (this->getAnalysisMode() == ANP_MODE_TRANSIENT)
          {
            currentOutput << " time=" << this->getTime() << std::endl;
          }
          else
          {
            currentOutput << " Step=" << analysis_event.step_ << std::endl;
          }
          outputCurrentHeaderLine = false;
        }
        std::string nodeName = this->getNodeNameFromLocalIndex( localID ); 
        currentOutput
          << "     solution I" << "(" << nodeName << ")=" << value << std::endl;
      } 
    }
    
    // if current output is requested then also loop over the lead-current vector
    if( diagnosticCurrentLimitGiven_ )
    {
      auto numLeadVars = (this->getDataStore())->leadCurrentSize;
      for( auto localID=0 ; localID < numLeadVars; localID++)
      {
        value = (*(this->getDataStore()->currLeadCurrentPtr))[localID];
        std::string nodeName("NF");
        if( fabs(value) > diagnosticCurrentLimit_)
        {
          NodeNameMap branchVarMap = outputManagerAdapter_.getBranchVarsNodeMap();
          auto mapItr = branchVarMap.begin();
          auto endItr = branchVarMap.end();
          bool nameFound = false;
          while( !nameFound && (mapItr != endItr))
          {
            if( mapItr->second == localID)
            {
              nodeName = mapItr->first;
              nameFound = true;
            }
            mapItr++;
          }
          if( outputCurrentHeaderLine )
          {
            currentOutput << " Current over limit value found in " 
              << analysis_event.outputType_ 
              << " analysis at " 
              << analysis_event.state_;
            if (this->getAnalysisMode() == ANP_MODE_TRANSIENT)
            {
              currentOutput << " time=" << this->getTime() << std::endl;
            }
            else
            {
              currentOutput << " Step=" << analysis_event.step_ << std::endl;
            }
            outputCurrentHeaderLine = false;
          }
          currentOutput
            << "     lead I" << "(" << nodeName << ")=" << value << std::endl;
        } 
      }
    
    }
    // send any accumulated output to the diagnostic file
    std::string vOutput(voltageOutput.str() );
    std::string iOutput(currentOutput.str() );
    std::string viOutput = vOutput + iOutput;
    std::stringstream combinedOutput;
    combinedOutput.flush();
    //
    // each process will have a different value for viOutput.  For clear output 
    // to the user we need to collect these strings on the lead processor 
    //
    if( pdsComm.isSerial() )
    {
      combinedOutput << viOutput;
    }
    else
    {
      int numProc = pdsComm.numProc();
      int thisProc = pdsComm.procID();
      int stringLenPerProcPreCom[ numProc ];
      int stringLenPerProc[ numProc ];
      for( int p=0; p<numProc; p++) 
      {
        stringLenPerProcPreCom[p] = 0;
        stringLenPerProc[p] = 0;
      }
      stringLenPerProcPreCom[thisProc] = viOutput.length();
      pdsComm.maxAll(stringLenPerProcPreCom, stringLenPerProc, numProc);
    
      // now we know the length of a diagnostic message on each processor.  So, accumulate them   
      for( int p=0; p<numProc; p++) 
      {
        if( stringLenPerProc[p] > 0 )
        {
          char *buffer = (char *)malloc( stringLenPerProc[p]+1 );
          if( stringLenPerProcPreCom[p] > 0 )
          {
            strcpy( buffer, viOutput.c_str() );
          }
          pdsComm.bcast( buffer, stringLenPerProc[p]+1, p );
          combinedOutput << std::string(buffer);
          free(buffer);
        }
      }
    }
    
    // only output if the stream is open
    if(diagnosticOutputStreamPtr_ != NULL)
    {
      (*diagnosticOutputStreamPtr_)  << combinedOutput.str(); 
    }
  }
  
  if(diagnosticDiscontinuityLimitGiven_)
  {
    // loop over the predictor vector compared to the solution vector and output any
    // values with an absolute difference greater than diagnosticDiscontinuityLimit_
    
    // stream buffers for the output to keep it organized.
    std::stringstream discOutput;
    discOutput.flush();
    bool outputDiscHeaderLine = true;

    // loop over the full solution vector 
    double value = 0.0;
    double xPredictor = 0.0;
    double qPredictor = 0.0;
    //auto numSolVars = (this->getDataStore())->solutionSize;
    auto numSolVars = (this->getDataStore())->currSolutionPtr->localLength();
    for( auto localID=0 ; localID < numSolVars; localID++)
    {
      value =(*(this->getDataStore()->currSolutionPtr))[localID];
      xPredictor =(*(this->getDataStore()->xn0Ptr))[localID];
      qPredictor =(*(this->getDataStore()->qn0Ptr))[localID];
      // Can't just look at the difference between the solution value and the 
      // xPredictor and qPredictor.  Most values will have all of their contribution
      // from x or q.  Thus if the difference between value and x is significant 
      // then q is likely zero and will be a false positive for a discontinuity.
      // So, get the absolute maximum of the x and q predictors and then use that
      // to compare with the solution value.
      double maxPredictor = std::max( fabs(xPredictor), fabs(qPredictor));
      auto globalID = ((dataStore_->builder_).getSolutionMap())->localToGlobalIndex( localID );
      char varType = this->getNodeTypeFromLocalIndex( localID ); 

      if( fabs(fabs(value)-maxPredictor) > diagnosticDiscontinuityLimit_) 
      {
        if( outputDiscHeaderLine )
        {
          discOutput << " Discontinuity over limit value found in " 
            << analysis_event.outputType_ 
            << " analysis at " 
            << analysis_event.state_;
          if (this->getAnalysisMode() == ANP_MODE_TRANSIENT)
          {
            discOutput << " time=" << this->getTime() << std::endl;
          }
          else
          {
            discOutput << " Step=" << analysis_event.step_ << std::endl;
          }
          outputDiscHeaderLine = false;
        }
        std::string nodeName = this->getNodeNameFromLocalIndex( localID ); 
        discOutput
          << "     " << varType << "(" << nodeName << ")=" << value
          << " xPredictor  = " << xPredictor
          << " qPredictor  = " << qPredictor 
          << " |diff| = " <<  fabs(fabs(value)-maxPredictor)<< std::endl;
      }
    }

    std::string discStringOutput = discOutput.str();
    std::stringstream combinedOutput;
    combinedOutput.flush();
    //
    // each process will have a different value for discStringOutput.  For clear output 
    // to the user we need to collect these strings on the lead processor 
    //
    if( pdsComm.isSerial() )
    {
      combinedOutput << discStringOutput;
    }
    else
    {
      int numProc = pdsComm.numProc();
      int thisProc = pdsComm.procID();
      int stringLenPerProcPreCom[ numProc ];
      int stringLenPerProc[ numProc ];
      for( int p=0; p<numProc; p++) 
      {
        stringLenPerProcPreCom[p] = 0;
        stringLenPerProc[p] = 0;
      }
      stringLenPerProcPreCom[thisProc] = discStringOutput.length();
      pdsComm.maxAll(stringLenPerProcPreCom, stringLenPerProc, numProc);
    
      // now we know the length of a diagnostic message on each processor.  So, accumulate them   
      for( int p=0; p<numProc; p++) 
      {
        if( stringLenPerProc[p] > 0 )
        {
          char *buffer = (char *)malloc( stringLenPerProc[p]+1 );
          if( stringLenPerProcPreCom[p] > 0 )
          {
            strcpy( buffer, discStringOutput.c_str() );
          }
          pdsComm.bcast( buffer, stringLenPerProc[p]+1, p );
          combinedOutput << std::string(buffer);
          free(buffer);
        }
      }
    }
    // only output if the stream is open
    if(diagnosticOutputStreamPtr_ != NULL)
    {
      (*diagnosticOutputStreamPtr_)  << combinedOutput.str(); 
      diagnosticOutputStreamPtr_->flush();
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool registerPkgOptionsMgr(
    AnalysisManager &analysis_manager, 
    IO::PkgOptionsMgr &options_manager)
{
  Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("DIAGNOSTIC");
  
  parameters.insert(Util::ParamMap::value_type("EXTREMALIMIT", Util::Param("EXTREMALIMIT", 0.0)));
  parameters.insert(Util::ParamMap::value_type("VOLTAGELIMIT", Util::Param("VOLTAGELIMIT", 0.0)));
  parameters.insert(Util::ParamMap::value_type("CURRENTLIMIT", Util::Param("CURRENTLIMIT", 0.0)));
  parameters.insert(Util::ParamMap::value_type("DISCLIMIT", Util::Param("DISCLIMIT", 0.0)));
  parameters.insert(Util::ParamMap::value_type("DIAGFILENAME", Util::Param("DIAGFILENAME", "XyceDiag.out")));
  
  options_manager.addCommandProcessor("OP", 
    IO::createRegistrationOptions(analysis_manager, &AnalysisManager::setOPAnalysisParams));

  options_manager.addCommandProcessor("SENS", 
    IO::createRegistrationOptions(analysis_manager, &AnalysisManager::setSensOptions));
    
  options_manager.addOptionsProcessor("DIAGNOSTIC",
    IO::createRegistrationOptions(analysis_manager, &AnalysisManager::setDiagnosticMode));

  return true;
}

} // namespace Analysis
} // namespace Xyce

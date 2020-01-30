//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
// Purpose       :
//
// Special Notes :
//
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date : 01/24/08
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_AnalysisManager_h
#define Xyce_N_ANP_AnalysisManager_h

#include <list>
#include <set>
#include <vector>

#include <N_ANP_fwd.h>
#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_LOA_fwd.h>
#include <N_NLS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TIA_fwd.h>
#include <N_TOP_fwd.h>
#include <N_UTL_fwd.h>

#include <N_ANP_AnalysisBase.h>
#include <N_ANP_AnalysisEvent.h>
#include <N_ANP_RegisterAnalysis.h>
#include <N_ANP_StepEvent.h>
#include <N_TIA_TIAParams.h>
#include <N_UTL_Factory.h>
#include <N_UTL_Listener.h>
#include <N_UTL_Stats.h>
#include <N_UTL_Timer.h>

namespace Xyce {
namespace Analysis {

typedef std::vector<Util::Factory<AnalysisBase, void> *> CreatorVector;
typedef std::set<Util::Factory<ProcessorBase, void> *> CreatorSet;

typedef Util::Notifier<StepEvent> StepEventNotifier;
typedef Util::Notifier<AnalysisEvent> AnalysisEventNotifier;
typedef Util::ListenerAutoSubscribe<StepEvent> StepEventListener;
typedef Util::ListenerAutoSubscribe<AnalysisEvent> AnalysisEventListener;

//-----------------------------------------------------------------------------
// Class         :
// Purpose       : This function converts between Nonlinear::AnalysisMode
//               : and AnalysisManager.h ANP_Analysis_Mode
// Special Notes :
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
Nonlinear::AnalysisMode nonlinearAnalysisMode(Mode mode);

//-----------------------------------------------------------------------------
// Function      : analysisModeName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Mar 16 13:56:09 2015
//-----------------------------------------------------------------------------
///
/// Returns the name of the analysis mode given by mode.
///
/// @param mode         Analysis mode
///
/// @return name of the mode
///
///
const char *analysisModeName(Mode mode);

//-----------------------------------------------------------------------------
// Class         : AnalysisManager
//
// Purpose       : This class manages, allocates, and sets up the
//                 various analysis types, such as DC, Tran, HB, etc.
//
// Special Notes : Some of this class was once in the TimeIntg::ControlAlgorithm
//                 class, which was set up back when Xyce only did transient
//                 simulations.  As we added analysis types it became necessary
//                 to refactor the code so that each analysis type had its own
//                 set of classes.
//
// Creator       : Eric Keiter, SNL
// Creation Date : 6/01/00 (TimeIntg::ControlAlgorithm, now deprecated)
// Creation Date : 1/24/08 (for this version of the class. Date is approximate)
//-----------------------------------------------------------------------------
class AnalysisManager : public StepEventNotifier,
                        public AnalysisEventNotifier,
                        public StepEventListener,
                        public AnalysisEventListener
{
public:
  AnalysisManager(
    const IO::CmdParse &        command_line,
    OutputMgrAdapter &          output_manager_adapter,
    Stats::Stat                 analysis_stat);
  virtual ~AnalysisManager();

private:
  AnalysisManager(const AnalysisManager &);
  AnalysisManager &operator=(const AnalysisManager &);

public:
  void notify(const StepEvent &step_event);
  void notify(const AnalysisEvent &analysis_event);

  void allocateAnalysisObject(AnalysisCreatorRegistry &analysis_registry);

  bool initializeSolverSystem(
    const TimeIntg::TIAParams & tia_params,
    Loader::Loader &            loader,
    Linear::System &            linear_system,
    Nonlinear::Manager &        nonlinear_manager,
    Device::DeviceMgr &         device_manager);

  void resetSolverSystem();

  // Execute the control loop for the set analysis type.
  bool run();

  // Gets the next time-step value.
  double getTime() const;

  // Gets the final time-step value.
  double getFinalTime() const;

  // Gets the initial time-step value.
  double getInitialTime() const;

  // Calls the time int. method to update the corrector derivatives.
  bool updateDerivs();

  // updates state vectors following initial solve with previous operating point constraints
  bool completeOPStartStep();

  // updates the State vectors. This function is called from the LOCA interface.
  bool completeHomotopyStep(
    Loader::NonlinearEquationLoader &   loader,
    const std::vector<std::string> &    paramNames,
    const std::vector<double> &         paramVals,
    Linear::Vector *                    solnVecPtr );

  bool failHomotopyStep(Loader::NonlinearEquationLoader &loader);

  // Prints out time loop information.
  bool printLoopInfo(int start, int finish);

  DCOPType getDCOPSolve() const;

  // Get the steady-state flag (true if time int mode is none)
  bool getDCOPFlag() const;

  // Get the dcop flag for transient.(true if doing DCOP for transient initialization)
  bool getTranOPFlag() const;

  // Get the dcop flag for AC.(true if doing DCOP for AC initialization)
  bool getACOPFlag() const;

  bool getACFlag() const;

  bool getACLinFlag() const;

  bool getNoiseFlag() const;

  bool getDCSweepFlag() const;

  bool getDotOpSpecified()
  {
    return dotOpSpecified_;
  }

  bool getSweepSourceResetFlag() const
  {
    return sweepSourceResetFlag_;
  }

  void setSweepSourceResetFlag (bool reset)
  {
    sweepSourceResetFlag_ = reset;
  }

  bool getTransientFlag () const;

  // Is the doubleDCOP algorithm enabled?
  bool getDoubleDCOPEnabled() const;

  void setTwoLevelMode(TwoLevelMode two_level_mode)
  {
    twoLevelMode_ = two_level_mode;
  }

  TwoLevelMode getTwoLevelMode() const
  {
    return twoLevelMode_;
  }

  // Get block analysis information for HB
  bool getBlockAnalysisFlag() const;

  // gets the index of the DCOP step.
  // 0 = nonlinear poisson, 1=full TCAD
  int getDoubleDCOPStep() const;

  // Gets/sets the step number.
  int getStepNumber() const;
  int getTranStepNumber();

  void setStepNumber(int step);
  void setTranStepNumber(int step);

    // This is true only at the beginning of integration, not at a breakpoint.
  bool getInitTranFlag() const;

  const IO::CmdParse &getCommandLine() const
  {
    return commandLine_;
  }

  const std::string &getNetlistFilename() const
  {
    return netlistFilename_;
  }

  // returns whether transient analysis is completed
  bool isSimulationComplete();

  // Gets the size of the restart data (bytes?).
  int getRestartDataSize( bool pack ) const;

  // Sets the DC sweep calculation parameters.
  bool setDCAnalysisParams(const Util::OptionBlock & paramsBlock);

  // Method to handle OP statements.
  bool setOPAnalysisParams(const Util::OptionBlock & paramsBlock);

  // // Sets the DCOP restart parameters.
  // bool setDCOPRestartParams(const Util::OptionBlock & OB);

  // Sets the AC Analysis options.
  bool setACAnalysisParams(const Util::OptionBlock & OB);

  // Sets the NOISE Analysis options.
  bool setNOISEAnalysisParams(const Util::OptionBlock & OB);

  // sets a time at which to pause the simulation
  void setPauseTime(double pauseTime, double initial_time);

  // returns time at which to pause the simulation
  double getPauseTime() const;

  // returns true if the simulation is currently paused
  bool isPaused() const;

  // Sets the sensitivity options.
  bool setSensOptions(const Util::OptionBlock & OB);

  // Registers the parallel services manager pointer.
  bool registerParallelServices(Parallel::Manager * pds_tmp);

  // Registers the elapsed time timer
  bool registerElapsedTimer(Util::Timer *);

  // Writes-out the restart data.
  bool dumpRestartData(char * buf, int bsize, int & pos, Parallel::Communicator * comm, bool pack);

  // Restores the restart data.
  bool restoreRestartData(char * buf, int bsize, int & pos, Parallel::Communicator * comm, bool pack );

  // Gets the solution variable data.
  bool getSolnVarData(const int & gid, std::vector< double > & varData) const;

  // Gets the state variable data.
  bool getStateVarData(const int & gid, std::vector< double > & varData) const;

  // Gets the store variable data.
  bool getStoreVarData(const int & gid, std::vector< double > & varData) const;

  // Sets the solution variable data.
  bool setSolnVarData(const int & gid, const std::vector< double > & varData);

  // Sets the state variable data.
  bool setStateVarData(const int & gid, const std::vector< double > & varData);

  // Sets the store variable data.
  bool setStoreVarData(const int & gid, const std::vector< double > & varData);

  // set/get beginning integration flag
  void setBeginningIntegrationFlag(bool bif);
  bool getBeginningIntegrationFlag() const;

  // set/get integration method
  void setIntegrationMethod(int im);
  int getIntegrationMethod();

  // Gets the total time spent in the linear solvers.
  double getTotalLinearSolutionTime() const;

  // Gets the total time spent in the residual load calculations.
  double getTotalResidualLoadTime() const;

  // Gets the total time spent in the Jacobian load calculations.
  double getTotalJacobianLoadTime() const;

  // set the next solution vector pointer.  Needed for NOX...
  bool setNextSolVectorPtr (Linear::Vector * solVecPtr);

  void createTimeIntegratorMethod(const TimeIntg::TIAParams &tia_params, const unsigned int integration_method);

  //-----------------------------------------------------------------------------
  // Function      : AnalysisManager::getTIAParams
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric R. Keiter, SNL, Computational Sciences
  // Creation Date : 11/05/04
  //-----------------------------------------------------------------------------
  const TimeIntg::TIAParams &getTIAParams() const;
  TimeIntg::TIAParams &getTIAParams();

  bool getSensFlag() const
  {
    return sensFlag_;
  }

  void addAnalysis(Util::Factory<AnalysisBase, void> *factory)
  {
    analysisCreatorVector_.push_back(factory);
  }

  void addProcessor(Util::Factory<ProcessorBase, void> *factory)
  {
    processorCreatorSet_.insert(factory);
  }

  void setDAEStateDerivFlag(bool state)
  {
    daeStateDerivFlag_ = state;
  }

public:
  TimeIntg::DataStore *getDataStore()
  {
    return dataStore_;
  }

  const AnalysisBase &getAnalysisObject() const
  {
    return *primaryAnalysisObject_;
  }

  AnalysisBase &getAnalysisObject()
  {
    return *primaryAnalysisObject_;
  }

  void setPrimaryAnalysisObject(AnalysisBase *primary)
  {
    primaryAnalysisObject_ = primary;
  }

  CreatorVector &getCreatorVector()
  {
    return analysisCreatorVector_;
  }

  Parallel::Manager *getPDSManager() const
  {
    return parallelManager_;
  }

  Parallel::Machine getComm() const;

  bool getSwitchIntegrator() const
  {
    return switchIntegrator_;
  }

  void setSwitchIntegrator(bool switch_itegrator)
  {
    switchIntegrator_ = switch_itegrator;
  }

  void setNextOutputTime(double next_output_time)
  {
    nextOutputTime_ = next_output_time;
  }

  double getNextOutputTime() const
  {
    return nextOutputTime_;
  }

  Util::Timer &getXyceTranTimer()
  {
    return xyceTranTimerPtr_;
  }

  OutputMgrAdapter &getOutputManagerAdapter() const
  {
    return outputManagerAdapter_;
  }

  TimeIntg::WorkingIntegrationMethod &getWorkingIntegrationMethod()
  {
    return *workingIntgMethod_;
  }

  const TimeIntg::WorkingIntegrationMethod &getWorkingIntegrationMethod() const
  {
    return *workingIntgMethod_;
  }

  TimeIntg::StepErrorControl &getStepErrorControl()
  {
    return *stepErrorControl_;
  }

  const TimeIntg::StepErrorControl &getStepErrorControl() const
  {
    return *stepErrorControl_;
  }

  Loader::NonlinearEquationLoader &getNonlinearEquationLoader()
  {
    return *nonlinearEquationLoader_;
  }

  void setAnalysisMode(Mode analysis_mode)
  {
    analysisMode_ = analysis_mode;
  }

  Mode getAnalysisMode() const
  {
    return analysisMode_;
  }

  double getSolverStartTime() const
  {
    return solverStartTime_;
  }

  void silenceProgress()
  {
    progressFlag_ = false;
  }

  bool getProgressFlag() const
  {
    return progressFlag_;
  }

  void pushActiveAnalysis(AnalysisBase *analysis)
  {
    currentAnalysisStack_.push_back(analysis);
  }

  void popActiveAnalysis()
  {
    currentAnalysisStack_.pop_back();
  }

  const AnalysisBase *getActiveAnalysis() const
  {
    return currentAnalysisStack_.front();
  }

  const AnalysisBase *getCurrentAnalysis() const
  {
    return currentAnalysisStack_.back();
  }
  void setResumeSimulation(bool resume)
  {
    resumeSimulation_ = resume;
  }

  bool getResumingSimulation() const
  {
    return resumeSimulation_;
  }

  bool getSavedAlready() const
  {
    return savedAlready_;
  }

  void setSavedAlready(bool saved_already)
  {
    savedAlready_ = saved_already;
  }

  double getCurrentFreq() const;

  bool getSparcalc() const;

  void setRFParamsRequested(const std::string & type);


private:
  const IO::CmdParse &                  commandLine_;                   ///< Command line object
  const std::string                     netlistFilename_;               ///< Netlist file name

  OutputMgrAdapter &                    outputManagerAdapter_;          ///< Output manager adapter

  TimeIntg::WorkingIntegrationMethod *  workingIntgMethod_;             ///< Working integration method
  TimeIntg::StepErrorControl *          stepErrorControl_;              ///< Pointer to the TIA step-error control object
  Loader::NonlinearEquationLoader *     nonlinearEquationLoader_;    ///< Pointer to the nonlinear equation loader

  Parallel::Manager *                   parallelManager_;               ///< Pointer to the parallel services manager
  TimeIntg::DataStore *                 dataStore_;                     ///< Data store object
  IO::ActiveOutput *                    activeOutput_;


  Mode                  analysisMode_;
  TwoLevelMode          twoLevelMode_;

  bool                  resumeSimulation_;              ///< Resume simulation from a paused transient
  bool                  blockAnalysisFlag_;             ///< HB Analysis (maybe something with MPDE too)
  bool                  daeStateDerivFlag_;             ///< .OPTIONS TIMEINT DAESTATEDERIV=
  bool                  dotOpSpecified_;                ///< Set if .OP
  bool                  gui_;                           ///< Set if -giu appears on command line
  bool                  progressFlag_;
  bool                  saveTimeGiven_;
  bool                  savedAlready_;
  bool                  sensFlag_;
  bool                  sweepSourceResetFlag_;
  bool                  switchIntegrator_;              ///< Set to true when Transient::integrationMethod_ is changed

  Util::Timer           xyceTranTimerPtr_;              /// Xyce timing utility for timing the transient simulation CPU time.
  Util::Timer *         elapsedTimerPtr_;               /// Xyce timing utility for timing elapsed run time

  double                solverStartTime_;

  double                nextOutputTime_;
  // double                nextRestartSaveTime_;

  AnalysisBase *        analysisObject_;                ///< .STEP, Dakota
  AnalysisBase *        primaryAnalysisObject_;         ///< .TRAN, .AC, .HB, ...

  std::vector<ProcessorBase *> analysisVector_;
  std::vector<AnalysisBase *> currentAnalysisStack_;

  CreatorVector                 analysisCreatorVector_;
  CreatorSet                    processorCreatorSet_;

  Stats::Stat           analysisStat_;

public:
    unsigned int breakPointRestartStep;
};

bool registerPkgOptionsMgr(AnalysisManager &analysis_manager, IO::PkgOptionsMgr &options_manager);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_AnalysisManager_h

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
//
// Purpose        : Transient analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/24/08
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_Transient_h
#define Xyce_N_ANP_Transient_h

#include <N_ANP_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_TOP_fwd.h>
#include <N_UTL_fwd.h>

#include <fstream>

#include <N_ANP_AnalysisBase.h>
#include <N_ANP_RegisterAnalysis.h>
#include <N_ANP_StepEvent.h>
#include <N_TIA_TIAParams.h>
#include <N_UTL_FixedQueue.h>
#include <N_UTL_Listener.h>
#include <N_UTL_OptionBlock.h>

class N_MPDE_Manager;

namespace Xyce {
namespace Analysis {

typedef Util::ListenerAutoSubscribe<StepEvent> StepEventListener;

//-------------------------------------------------------------------------
// Class         : Transient
// Purpose       : Transient analysis class
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class Transient : public AnalysisBase, public StepEventListener
{
public:
  Transient(
    AnalysisManager &                   analysis_manager,
    //Linear::System &                    linear_system,
    Linear::System *                    linear_system_ptr,
    Nonlinear::Manager &                nonlinear_manager,
    Loader::Loader &                    loader,
    Topo::Topology &                    topology,
    IO::InitialConditionsManager &      initial_conditions_manager,
    IO::RestartMgr &                    restart_manager,
    OutputAdapter *                     output_adapter = 0,
    HB *                                hb_analysis = 0,        // HACK TO GET US MOVING FORWARD, should pass output adapter or something
    N_MPDE_Manager *                    mpde_manager = 0);      // HACK TO GET US MOVING FORWARD, should pass output adapter or something

  virtual ~Transient()
  {}

  void notify(const StepEvent &event);

  void setTIAParams(const TimeIntg::TIAParams &tia_params)
  {
    tiaParams_ = tia_params;
  }

  const TimeIntg::TIAParams &getTIAParams() const
  {
    return tiaParams_;
  }

  TimeIntg::TIAParams &getTIAParams()
  {
    return tiaParams_;
  }

  bool setAnalysisParams(const Util::OptionBlock & paramsBlock);
  bool setTimeIntegratorOptions(const Util::OptionBlock &option_block);
  bool setOutputOptions(const Util::OptionBlock &option_block);
  bool setSensitivityOptions(const Util::OptionBlock &option_block);
  bool setSensAnalysisParams(const Util::OptionBlock &option_block);


  void registerLinearSystem(Linear::System * linear_system_ptr)
  {
    linearSystemPtr_ = linear_system_ptr;
  }

  void stepCallback();
  void registerParentAnalysis(AnalysisBase * parentPtr)
  {
    parentAnalysisPtrVec_.push_back(parentPtr);
  };

  void finalExpressionBasedSetup();

protected:
  bool doRun();
  bool doInit();
  bool resuming();
  bool doLoopProcess();
  bool doTranOP ();
  bool doProcessSuccessfulStep();
  bool doProcessFailedStep();
  bool doFinish();
  bool doHandlePredictor();

  void allocateDODP();
  bool doTransientAdjointSensitivity ();
  bool saveTransientAdjointSensitivityInfo ();
  bool saveTransientAdjointSensitivityInfoDCOP ();

  void transientLambdaOutputHeader ();
  void transientLambdaOutputZone (int itGlobal);
  void transientLambdaOutput (int it);
  void transientLambdaOutputFooter ();

  void transientAdjointSensOutput (int itGlobal);
  void transientAdjointSensOutputFooter ();

public:
  bool processSuccessfulDCOP();
  bool processFailedDCOP();

  bool resetForHB();
  bool finalVerboseOutput();

  void printStepHeader(std::ostream &os);
  void printProgress(std::ostream &os);

  // mixed-signal specific
  bool mixedSignalStep(double maxTimeStepFromHabanero);
  bool finalizeMixedSignalStep();

  // Two Level specific
  bool twoLevelStep();

  bool getDCOPFlag() const
  {
    return dcopFlag_;
  }

  int getDCStats() { return dcStats; }
  int getTranStats() { return tranStats; }

  void setSaveTimeSteps(bool save_time_steps)
  {
    saveTimeStepsFlag = save_time_steps;
  }

private:
  void preMixedSignalStepDetails(double maxTimeStepFromHabanero);
  void noopOutputs ();
  void tranopOutputs ();
  void tranStepOutputs ();

  void takeAnIntegrationStep_();

  bool retakeAndAcceptTimeStep( double aTimeStep );

  void logQueuedData();

  void outputFailedStepData();

private:
  Parallel::Machine                     comm_;
  AnalysisManager &                     analysisManager_;
  Loader::Loader &                      loader_;
  //Linear::System &                      linearSystem_;
  Linear::System *                      linearSystemPtr_;
  Nonlinear::Manager &                  nonlinearManager_;
  Topo::Topology &                      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;
  IO::RestartMgr &                      restartManager_;
  OutputMgrAdapter &                    outputManagerAdapter_;
  OutputAdapter *                       outputAdapter_;
  TimeIntg::TIAParams                   tiaParams_;
  bool                                  sensFlag_;
  bool                                  solveAdjointSensitivityFlag_;
  bool                                  solveDirectSensitivityFlag_;
  bool                                  outputTransientLambda_;
  bool                                  outputAdjointSensitivity_;
  unsigned int                          initialIntegrationMethod_;
  bool                                  firstTranOutput_;
  bool                                  fullAdjointTimeRange_;

  // a flag to indicate of the simulation is paused
  bool isPaused;
  bool dcopFlag_;               // true if this is a DCOP calculation.

  double startDCOPtime;
  double startTRANtime_;
  double endTRANtime_;
  bool quiet_;                  // command line arg -quiet is present

  bool historyTrackingOn_;      // bool to indicate if history tracking is on.

  // These are used to track the minimum estimated error over tol. of failed
  // time steps so that if we're going to exit with a time step too small error,
  // we can have the option of accepting whatever step had the minimum error.
  double minEstErrorOverTol;
  int    stepNumberAtMinEstErrorOverTol;
  double timeStepAtMinEstErrorOverTol;

  std::string                         maxTimeStepExpressionString_;
  Util::ExpressionData *              maxTimeStepExpression_;

  // here we store stats on the last few time steps
  // to report if Xyce fails. A user can use this info
  // to figure how how to make the simulation work.
  // The number of items saved is set in the constructor
  Util::FixedQueue<double> timeQueue_;
  Util::FixedQueue<double> timeStepQueue_;
  Util::FixedQueue<int> stepStatusQueue_;
  Util::FixedQueue<double> estErrorOverTolQueue_;
  Util::FixedQueue<int> nonlinearSolverStatusQueue_;
  Util::FixedQueue<int> nonlinearSolverNumIterationsQueue_;
  Util::FixedQueue<double> nonlinearSolverMaxNormQueue_;
  Util::FixedQueue<double> nonlinearSolverMaxNormIndexQueue_;
  // Util::FixedQueue<double> nonlinearSolverNormQueue_;

  bool          firstTime;
  double        oldPercentComplete;
  double        startSimTime;
  double        currRestartSaveTime_;
  double        nextRestartSaveTime_;
  int           dcStats;
  int           tranStats;
  double        exitTime;               ///< Exit when it exceeds this time.
  int           exitStep;               ///< Exit after taking this many steps.
  unsigned int  integrationMethod;      ///< Time-integration method


  // Xyce tracks how the last few timesteps have gone for error reporting
  // if Xyce has to exit.  This var specifies how many time steps (both
  // passing and failing) to keep in its history.  If it's zero, then
    // history tracking is off.
  int historyTrackingDepth;

  bool          passNLStall;            ///< option to pass some non-linear solver failures
  bool          saveTimeStepsFlag;      ///< flag to save timestpes in data store for later use
  bool          condTestFlag;           ///< flag for conductance test
  std::vector<std::string> condTestDeviceNames; ///< names for conductance test

  HB *                  hbAnalysis_;
  N_MPDE_Manager *      mpdeManager_;

  // Sensitivity data.  
  int numSensParams_;
  std::vector<double>   objectiveVec_;
  std::vector<double>   dOdpVec_;
  std::vector<double>   dOdpAdjVec_;
  std::vector<double>   scaled_dOdpVec_;
  std::vector<double>   scaled_dOdpAdjVec_;

  //int maxParamStringSize_;
  std::vector<std::string> paramNameVec_;

  // finite difference variables
  int difference_;
  double sqrtEta_;
  bool sqrtEtaGiven_;
  bool forceFD_;
  bool forceDeviceFD_;
  bool forceAnalytic_;
  bool newLowMem_;
  bool sparseAdjointStorage_;

  double adjointBeginTime_;
  bool adjointBeginTimeGiven_;
  double adjointFinalTime_ ;
  bool adjointFinalTimeGiven_ ;

  std::vector<double> adjointTimePoints_;
  bool adjointTimePointsGiven_;

  std::vector<double> outputTimePoints_;
  bool outputTimePointsGiven_;

  std::vector<double> userBreakPoints_;
  bool userBreakPointsGiven_;

  // temporary debug output files:
  std::ofstream         lambdaFile;
  std::ofstream         sensitivityFile;

  // step management
  bool resetForStepCalledBefore_;


  std::vector < AnalysisBase * > parentAnalysisPtrVec_;
};

std::vector<double> computeOutputInterpolationTimes(
  double                        current_time,
  double                        next_output_time, 
  double                        final_output_time, 
  double                        initial_output_interval,
  const IO::IntervalVector &    output_intervals);

double updateOutputTime(
  double                        current_time,
  double                        next_output_time,
  double                        final_output_time, 
  double                        initial_output_interval,
  const IO::IntervalVector &    output_intervals,
  const std::vector<double> &   outputTimePoints,
  bool                          outputTimePointsGiven);

bool testOutputTime(
  double                        current_time,
  double                        next_output_time,
  double                        start_time);

bool testSaveOutputTime(
  Analysis::AnalysisManager &           analysis_manager,
  IO::InitialConditionsManager &        initial_conditions_manager);

// Queries about output and restart times
bool testRestartSaveTime(
  AnalysisManager &             analysis_manager,
  IO::RestartMgr &              restart_manager,
  double                        current_time,
  double &                      curr_restart_save_time,
  double &                      next_restart_save_time);

bool registerTransientFactory(FactoryBlock &factory_block);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_Transient_h

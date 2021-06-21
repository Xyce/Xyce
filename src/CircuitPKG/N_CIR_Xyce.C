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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/27/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <stdexcept>
#include <ctime>
#include <numeric>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_DIRECT_H
#include <direct.h>
#endif

// ----------   Xyce Includes   ----------

#include <N_CIR_Xyce.h>

#include <N_DEV_fwd.h>
#include <N_DEV_ADC.h>
#include <N_DEV_Algorithm.h>
#include <N_DEV_Configuration.h>
#include <N_DEV_DAC.h>
#include <N_DEV_Device.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_ExternalSimulationData.h>
#include <N_DEV_OpBuilders.h>
#include <N_DEV_Print.h>
#include <N_DEV_RegisterDevices.h>
#include <N_DEV_SolverState.h>

#include <N_ERH_ErrorMgr.h>
#include <N_ERH_Messenger.h>

#include <N_IO_NetlistImportTool.h>

#include <N_IO_OutputMgr.h>
#include <N_IO_ParsingMgr.h>
#include <N_IO_OpBuilders.h>
#include <N_IO_MeasureManager.h>
#include <N_IO_FourierMgr.h>
#include <N_IO_FFTMgr.h>
#include <N_IO_InitialConditions.h>
#include <N_IO_LoadManager.h>
#include <N_IO_OutputMacroResults.h>
#include <N_IO_OutputResponse.h>
#include <N_IO_OutputResults.h>
#include <N_IO_PrintDeviceCount.h>
#include <N_IO_RestartMgr.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_ParsingHelpers.h>

#include <N_LAS_fwd.h>
#include <N_LAS_System.h>
#include <N_LAS_Builder.h>
#include <N_LAS_QueryUtil.h>

#include <N_LOA_CktLoader.h>

#include <N_NLS_Manager.h>

#include <N_PDS_Manager.h>
#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_DCSweep.h>
#include <N_ANP_NOISE.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_OpBuilders.h>
#include <N_ANP_RegisterAnalysis.h>

#include <N_TIA_StepErrorControl.h>
#include <N_TIA_WorkingIntegrationMethod.h>

#include <N_TOP_Topology.h>
#include <N_TOP_CktNode.h>

#include <N_UTL_Algorithm.h>
#include <N_UTL_CheckIfValidFile.h>
#include <N_UTL_Misc.h>
#include <N_UTL_Timer.h>
#include <N_UTL_Expression.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_PrintStats.h>
#include <N_UTL_Platform.h>
#include <N_UTL_JSON.h>
#include <N_UTL_SendCURL.h>

#include <N_UTL_Version.h>
#include <N_UTL_BreakPoint.h>

#include <N_DEV_DeviceSupport.h>
#include <mainXyceExpressionGroup.h>
#include <N_DEV_DeviceSupport.h>

namespace Xyce {
namespace Circuit {

#ifdef Xyce_TRACKING_URL
static const char *trackingURL = Xyce_TRACKING_URL;
#else
static const char *trackingURL = 0;
#endif

//--------  Global Declarations ------------
void report_handler(const char *message, unsigned report_mask);

namespace {

int s_errorWrap = 78;

struct ADCDeviceInstanceParameterOp: public Device::DeviceInstanceOp
{
  ADCDeviceInstanceParameterOp(std::map<std::string, std::map<std::string, double> > &adc_device_parameter_map)
      : ADCDeviceParameterMap_(adc_device_parameter_map)
  {}

  virtual bool operator()(Device::DeviceInstance *instance) {
    Device::ADC::Instance *adc_instance = dynamic_cast<Device::ADC::Instance *>(instance);
      if (adc_instance) {
        adc_instance->getParam("WIDTH", ADCDeviceParameterMap_[instance->getName().getEncodedName()]["width"]);
        adc_instance->getParam("R", ADCDeviceParameterMap_[instance->getName().getEncodedName()]["r"]);
        adc_instance->getModel().getParam("LOWERVOLTAGELIMIT", ADCDeviceParameterMap_[instance->getName().getEncodedName()]["lowerVoltageLimit"]);
        adc_instance->getModel().getParam("UPPERVOLTAGELIMIT", ADCDeviceParameterMap_[instance->getName().getEncodedName()]["upperVoltageLimit"]);
        adc_instance->getModel().getParam("SETTLINGTIME", ADCDeviceParameterMap_[instance->getName().getEncodedName()]["settlingTime"]);
      }

      return true;
    }

  std::map<std::string, std::map<std::string, double> > &ADCDeviceParameterMap_;
};

struct TimeVoltagePairsOp: public Device::DeviceInstanceOp
{
  TimeVoltagePairsOp(std::map<std::string, std::vector< std::pair<double,double> > >&time_voltage_map)
    : TimeVoltageMap_(time_voltage_map)
  {}

  virtual bool operator()(Device::DeviceInstance *instance) {
    Device::ADC::Instance &adc_instance = static_cast<Device::ADC::Instance &>(*instance);

    std::vector<std::pair <double,double> > TmpVec;
    adc_instance.getTVVEC(TmpVec);
    double current_time = adc_instance.getSolverState().currTime_;
    double * solVector = adc_instance.getExternData().nextSolVectorRawPtr;

    double vPos = solVector[adc_instance.getLIPos()];
    double vNeg = solVector[adc_instance.getLINeg()];
    TmpVec.push_back(std::pair<double,double>(current_time, vPos - vNeg));

    TimeVoltageMap_[instance->getName().getEncodedName()] = TmpVec;

    return true;
  }

  std::map<std::string, std::vector< std::pair<double,double> > > &      TimeVoltageMap_;
};

struct TimeStatePairsOp: public Device::DeviceInstanceOp
{
  TimeStatePairsOp(std::map<std::string, std::vector< std::pair<double,int> > >&time_state_map)
    : TimeStateMap_(time_state_map)
  {}

  virtual bool operator()(Device::DeviceInstance *instance) {
    Device::ADC::Instance &adc_instance = static_cast<Device::ADC::Instance &>(*instance);

    std::vector<std::pair <double,double> > TmpVec;
    std::vector<std::pair <double,int> > stateVec;
    adc_instance.getTVVEC(TmpVec);
    for (std::vector<std::pair <double,double> >::iterator it=TmpVec.begin(); it!=TmpVec.end();it++)
    {
      stateVec.push_back(std::pair<double,int>(it->first,adc_instance.deltaVToStateVal(it->second)));
    }

    double current_time = adc_instance.getSolverState().currTime_;
    double * solVector = adc_instance.getExternData().nextSolVectorRawPtr;

    double vPos = solVector[adc_instance.getLIPos()];
    double vNeg = solVector[adc_instance.getLINeg()];
    stateVec.push_back(std::pair<double,int>(current_time, adc_instance.deltaVToStateVal(vPos - vNeg)));
    
    TimeStateMap_[instance->getName().getEncodedName()] = stateVec;

    return true;
  }

  std::map<std::string, std::vector< std::pair<double,int> > > &      TimeStateMap_;
};

} // namespace <unnamed>


//-----------------------------------------------------------------------------
// Function      : Simulator::Simulator
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
Simulator::Simulator(Parallel::Machine comm)
  : runState_(START),
    comm_(comm),
    parsingManager_(0),
    deviceManager_(0),
    analysisRegistry_(0),
    processorRegistry_(0),
    topology_(0),
    linearSystem_(0),
    builder_(0),
    analysisManager_(0),
    circuitLoader_(0),
    outputManagerAdapter_(0),
    nonlinearManager_(0),
    parallelManager_(0),
    opBuilderManager_(0),
    outputManager_(0),
    measureManager_(0),
    fourierManager_(0),
    fftManager_(0),
    initialConditionsManager_(0),
    outputResponse_(0),
    restartManager_(0),
    loadManager_(0),
    optionsManager_(0),
    rootStat_(Stats::createRootStat("Xyce", Stats::StatSet(Stats::STAT_ALL))),
    analysisStat_("Analysis", rootStat_),
    auditJSON_(),
    XyceTimerPtr_(0),
    ElapsedTimerPtr_(0),
    commandLine_(),
    hangingResistor_(),
    externalNetlistParams_(),
    dacDeviceMap_(),
    adcDeviceMap_()
{
  previousReportHandler_ = set_report_handler(report_handler);
  Xyce::Report::reset_message_counts();

  rootStat_.start();

  TimeIntg::registerTimeIntegrationMethods();
}

//-----------------------------------------------------------------------------
// Function      : ~Xyce
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
Simulator::~Simulator()
{
  delete analysisRegistry_;
  delete processorRegistry_;
  delete parsingManager_;
  delete deviceManager_;
  delete nonlinearManager_;
  delete outputManagerAdapter_;
  delete circuitLoader_;
  delete analysisManager_;
  delete outputManager_;
  delete measureManager_;
  delete fourierManager_;
  delete fftManager_;
  delete loadManager_;
  delete initialConditionsManager_;
  delete outputResponse_;
  delete opBuilderManager_;
  delete linearSystem_;
  delete builder_;
  delete XyceTimerPtr_;
  delete ElapsedTimerPtr_;
  delete parallelManager_;
  delete topology_;
  delete restartManager_;
  delete optionsManager_;

  set_report_handler(previousReportHandler_);

  Stats::deleteRootStat(rootStat_);
}

Loader::CktLoader &
Simulator::getCircuitLoader() {
  return *circuitLoader_;
}


Analysis::AnalysisManager *
Simulator::newAnalysisManager(
  const IO::CmdParse &                command_line,
  IO::RestartMgr &                    restart_manager,
  Analysis::OutputMgrAdapter &        output_manager_adapter,
  Stats::Stat                         analysis_stat)
{
  return new Analysis::AnalysisManager(command_line, output_manager_adapter, analysis_stat);
}


//---------------------------------------------------------------------------
// Function      : Simulator::setNetlistParameters
// Purpose       : This passes a vector of pairs "key" "value" that will
//                 be substituted during the processing of the netlist.  This
//                 more easily allows Dakota to change any netlist parameter
//                 during netlist setup.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and MEMS Modeling
// Creation Date : 10/9/2008
//---------------------------------------------------------------------------
void Simulator::setNetlistParameters( const std::vector< std::pair< std::string, std::string > > & externalParams )
{
  externalNetlistParams_ = externalParams;
}

//---------------------------------------------------------------------------
// Function      : Simulator::setNetlistParameters
// Purpose       : Call through to the output manager to set the suffix to
//                 be used on the output file, as in circuit + suffix + prn
//                 This is useful in Dakota controlled runs to keep each
//                 simulation from overwritting the last one.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and MEMS Modeling
// Creation Date : 10/9/2008
//---------------------------------------------------------------------------
void Simulator::setOutputFileSuffix( const std::string newSuffix )
{
  if( outputManager_ )
  {
    outputManager_->setOutputFilenameSuffix( newSuffix );
  }
}

//-----------------------------------------------------------------------------
// Function      : Simulator::doAllocations_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool Simulator::doAllocations_()
{
  const std::string &netlist_filename = commandLine_.getArgumentValue("netlist");

  analysisRegistry_         = new Analysis::AnalysisCreatorRegistry();
  processorRegistry_        = new Analysis::ProcessorCreatorRegistry();
  parsingManager_           = new IO::ParsingMgr(commandLine_);
  opBuilderManager_         = new Util::Op::BuilderManager();
  optionsManager_           = new IO::PkgOptionsMgr();
  nonlinearManager_         = new Nonlinear::Manager(commandLine_);
  topology_                 = new Topo::Topology(commandLine_, hangingResistor_, *parallelManager_);
  deviceManager_            = new Device::DeviceMgr(comm_, *topology_, *opBuilderManager_, commandLine_);
  outputManager_            = new IO::OutputMgr(commandLine_, *opBuilderManager_, *topology_);
  outputResponse_           = new IO::OutputResponse();
  measureManager_           = new IO::Measure::Manager(commandLine_);
  fourierManager_           = new IO::FourierMgr(commandLine_);
  fftManager_               = new IO::FFTMgr(commandLine_);
  loadManager_              = new IO::LoadManager();
  initialConditionsManager_ = new IO::InitialConditionsManager(netlist_filename);
  restartManager_           = new IO::RestartMgr();
  builder_                  = new Linear::Builder();
  linearSystem_             = new Linear::System();
  outputManagerAdapter_     = new Analysis::OutputMgrAdapter(comm_, *outputManager_, *measureManager_, *fourierManager_,
                                                             *fftManager_, *deviceManager_);

  analysisManager_          = newAnalysisManager(commandLine_, *restartManager_, *outputManagerAdapter_, analysisStat_);
  circuitLoader_            = new Loader::CktLoader(*deviceManager_, *builder_);

  Util::subscribe<Analysis::StepEvent>(*analysisManager_, *fourierManager_);
  Util::subscribe<Analysis::StepEvent>(*analysisManager_, *fftManager_);
  Util::subscribe<Analysis::StepEvent>(*analysisManager_, *measureManager_);
  Util::subscribe<Analysis::StepEvent>(*analysisManager_, *outputManager_);
  Util::subscribe<Analysis::StepEvent>(*analysisManager_,  *outputManagerAdapter_);
  Util::subscribe<Analysis::StepEvent>(*analysisManager_, *deviceManager_);

  Device::registerOpBuilders(*opBuilderManager_, comm_, *deviceManager_);
  IO::registerOpBuilders(*opBuilderManager_, comm_, *outputManager_, *analysisManager_);
  IO::registerOpBuilders(*opBuilderManager_, comm_, *measureManager_);
  Analysis::registerOpBuilders(*opBuilderManager_, comm_, *analysisManager_);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::doRegistrations_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool Simulator::doRegistrations_()
{
  const std::string &netlist_filename = commandLine_.getArgumentValue("netlist");

  bool bsuccess = true;
  bool bs1 = true;

  Analysis::FactoryBlock factory_block(*analysisRegistry_, *processorRegistry_,
                                       *optionsManager_, *analysisManager_, *outputManager_, *linearSystem_,
                                       *nonlinearManager_, *circuitLoader_, *deviceManager_, *builder_, *topology_,
                                       *initialConditionsManager_, *restartManager_);
  Analysis::registerAnalysisFactory(factory_block);

  IO::registerOutputResultsFactory(factory_block, comm_);

  // Register options processors
  bs1 = IO::registerPkgOptionsMgr(*parsingManager_, *optionsManager_ ); bsuccess = bsuccess && bs1;
  bs1 = Device::registerPkgOptionsMgr(*deviceManager_, *optionsManager_ ); bsuccess = bsuccess && bs1;
  bs1 = Topo::registerPkgOptionsMgr(*topology_, *optionsManager_ ); bsuccess = bsuccess && bs1;
  bs1 = IO::registerPkgOptionsMgr(*restartManager_, *optionsManager_, Parallel::size(comm_), Parallel::rank(comm_)); bsuccess = bsuccess && bs1;
  bs1 = IO::registerPkgOptionsMgr(*outputManager_, *optionsManager_); bsuccess = bsuccess && bs1;
  bs1 = IO::registerPkgOptionsMgr(*loadManager_, *optionsManager_); bsuccess = bsuccess && bs1;
  bs1 = IO::Measure::registerPkgOptionsMgr(*measureManager_, *optionsManager_); bsuccess = bsuccess && bs1;
  bs1 = IO::registerPkgOptionsMgr(*fourierManager_, *optionsManager_); bsuccess = bsuccess && bs1;
  bs1 = IO::registerPkgOptionsMgr(*fftManager_, *optionsManager_); bsuccess = bsuccess && bs1;
  bs1 = IO::registerPkgOptionsMgr(*initialConditionsManager_, *optionsManager_); bsuccess = bsuccess && bs1;
  bs1 = Analysis::registerPkgOptionsMgr(*analysisManager_, *optionsManager_ ); bsuccess = bsuccess && bs1;
  bs1 = Nonlinear::registerPkgOptionsMgr(*nonlinearManager_, *optionsManager_ ); bsuccess = bsuccess && bs1;

  // Device Manager registrations
  bs1 = deviceManager_->registerNonlinearSolver(nonlinearManager_); bsuccess = bsuccess && bs1;
  bs1 = deviceManager_->registerAnalysisManager(analysisManager_);  bsuccess = bsuccess && bs1;

  // Analysis manager registrations:
  bs1 = analysisManager_->registerParallelServices(parallelManager_); bsuccess = bsuccess && bs1;
  bs1 = analysisManager_->registerElapsedTimer(ElapsedTimerPtr_); bsuccess = bsuccess && bs1;

  // Linear Solver registrations:
  bs1 = linearSystem_->registerPDSManager(parallelManager_);    bsuccess = bsuccess && bs1;
  bs1 = linearSystem_->registerBuilder( builder_ ); bsuccess = bsuccess && bs1;

  bs1 = builder_->registerPDSManager(parallelManager_); bsuccess = bsuccess && bs1;
  bs1 = builder_->registerQueryUtil(topology_->getLinearSolverUtility()); bsuccess = bsuccess && bs1;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::setupTopology
// Purpose       : This function handles a lot of the initial setup.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/08/00
//-----------------------------------------------------------------------------
Simulator::RunStatus
Simulator::setupTopology()
{
  // topology query's device manager to see if any devices are bad (i.e. a resistor with zero resistance)
  // if so, a list of nodes to be supernoded is created
  topology_->verifyNodesAndDevices(*deviceManager_);

  // create a union of the supernode list on all processors
  topology_->mergeOffProcTaggedNodesAndDevices();

  // combine nodes into supernodes and remove now redundant devices (i.e. those only connected to 1 processor )
  topology_->removeTaggedNodesAndDevices();

  return SUCCESS;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::setUpMatrixStructure_
// Purpose       : This function needs to set up the various linear algebra
//                 entities which are owned by the LAS system class.  This
//                 includes, most importantly, the Jacobian matrix and the
//                 right hand side vector.  It should also set the solution
//                 vector size and the state vector size.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/08/00
//-----------------------------------------------------------------------------
bool Simulator::setUpMatrixStructure_()
{
  Stats::Stat matrixStat("Setup Matrix Structure", rootStat_);
  Stats::TimeBlock mat(matrixStat);

  builder_->generateParMaps();
  builder_->generateGraphs();

  linearSystem_->initializeSystem();
 
  topology_->registerLIDswithDevs();

  deviceManager_->setupExternalDevices(*parallelManager_->getPDSComm());

  int lasSize = linearSystem_->numGlobalRows();
  Xyce::lout() << "***** Number of Unknowns = " << lasSize << std::endl;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::doInitializations_
// Purpose       : This function calls "initializeAll" functions in all the
//                 packages which have them.  Packages that have the LAS system
//                 class registered with them will have one of these functions.
//
// Special Notes : This is called once the size of the linear system is known,
//                 as various packages will need to allocate vectors, matrices,
//                 etc.
//
//                 These probably can be done in any order, but to be safe,
//                 make sure for now that the time integrator's initializeAll
//                 function is called first.  The other two are essentially
//                 secondary registrations, while the TIA one includes a lot of
//                 allocations.
//
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
bool Simulator::doInitializations_()
{
  bool bsuccess = true;
  bool bs1 = true;

  analysisManager_->allocateAnalysisObject(*analysisRegistry_);
  bsuccess = bsuccess && bs1;

  analysisManager_->initializeSolverSystem(analysisManager_->getTIAParams(), *circuitLoader_, *linearSystem_, *nonlinearManager_, *deviceManager_);

  bs1 = deviceManager_->initializeAll(*linearSystem_);
  bsuccess = bsuccess && bs1;

  bs1 = nonlinearManager_->initializeAll(*analysisManager_, analysisManager_->getNonlinearEquationLoader(), *linearSystem_, *analysisManager_->getDataStore(), *parallelManager_, *initialConditionsManager_, *outputManager_, *topology_);
  bsuccess = bsuccess && bs1;

  if (restartManager_->isRestarting())
    restartManager_->restoreRestartData(*parallelManager_->getPDSComm(), *topology_, *analysisManager_, *deviceManager_, restartManager_->path() );

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::runSolvers_
// Purpose       : This function runs the solvers.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
bool Simulator::runSolvers_()
{
  return analysisManager_->run();
}


//-----------------------------------------------------------------------------
// Function      : Simulator::run
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/19/01
//-----------------------------------------------------------------------------
Simulator::RunStatus
Simulator::run(
  int           argc,
  char **       argv)
{
  RunStatus run_status = SUCCESS;

  try
  {
    run_status = initialize(argc, argv);
  }
  catch (std::exception &x)
  {
    Xyce::lout() << "Exception:\n" << x.what() << std::endl;
    run_status = ERROR;
  }

  if (run_status == ERROR)
  {
    if (runState_ > PARSE_COMMAND_LINE) {
      reportTotalElapsedTime();
      Xyce::lout() << "Xyce Initialization Phase failed" << std::endl;
    }

    return ERROR;
  }

  if (run_status == SUCCESS)
  {
    try
    {
      run_status = runSimulation();
    }
    catch (std::exception &x)
    {
      Xyce::lout() << "Exception " << x.what() << std::endl;
      run_status = ERROR;
    }
  }

  if (run_status != DONE)
  {
    try
    {
      finalize();
    }
    catch (std::exception &x)
    {
      Xyce::lout() << "Exception " << x.what() << std::endl;
      run_status = ERROR;
    }
  }

  return run_status;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::runSimulation
// Purpose       : Main simulation driver.
// Special Notes : Not private as this is also called from N_DAK_DakotaInterface
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
Simulator::RunStatus Simulator::runSimulation()
{
  return runSolvers_() ? SUCCESS : ERROR;
}

// ---------------------------------------------------------------------------
// API METHODS NEEDED FOR MIXED-SIGNAL and other external applications
//
//-----------------------------------------------------------------------------
// Function      : Simulator::initialize
// Purpose       : capture all "initialization-type" activities in one
//                 method
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 02/17/2004
//-----------------------------------------------------------------------------
Simulator::RunStatus Simulator::initialize(
  int           argc,
  char **       argv)
{
  RunStatus run_status = SUCCESS;
  
  run_status = initializeEarly(argc,argv);
  if (run_status == SUCCESS)
    run_status = initializeLate();

  return run_status;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::initializeEarly
// Purpose       : capture early "initialization-type" activities in one
//                 method
// Special Notes : This is basically the first half of the what used to be
//                 the initialize function.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 02/28/2017
//-----------------------------------------------------------------------------
///
/// Parse command line, read netlist, check for errors, do early reporting.
/// STOP before Topology needs to query devices for number of internal
/// variables and jacobian stamp.
///
/// @param[in] argc Number of command line arguments
/// @param[in] argv Command line arguments
/// @return simulator status from type RunStatus
///
Simulator::RunStatus Simulator::initializeEarly(
  int           argc,
  char **       argv)
{
  runState_ = PARALLEL_INIT;

  // Setup the Parallel Mgr. with a default load-balance based on the numProc value.
  parallelManager_ = new Parallel::Manager(argc, argv, comm_);
  if (comm_ == MPI_COMM_NULL)
    comm_ = parallelManager_->getPDSComm()->comm();

  Report::registerComm(comm_);

  if (DEBUG_ALL_PROCS_SAME_WD)
  {
    char path[255];

    if (Parallel::rank(comm_) == 0)
#ifdef HAVE_WIN_DIRCOMMANDS
      _getcwd(path, sizeof(path));
#else
      getcwd(path, sizeof(path));
#endif
      
    Parallel::Broadcast(comm_, path, sizeof(path), 0);

    if (Parallel::rank(comm_) != 0 )
#ifdef HAVE_WIN_DIRCOMMANDS
      _chdir(path);
#else
      chdir(path);
#endif
  }


  runState_ = PARSE_COMMAND_LINE;
  // Parse command line arguments
  IO::CmdParse::ParseState parse_state = commandLine_.parseCommandLine(comm_, argc, argv);
  if (parse_state == IO::CmdParse::DONE)
    return DONE;
  else if (parse_state == IO::CmdParse::ERROR) {
    return ERROR;
  }

  s_errorWrap = commandLine_.argExists("-error-test") ? 2024 : 78;

  // Load any device plugins requested
  // This has to be done very early, because -param, -doc, and -doc_cat
  // all need the plugins loaded in order to work.
  // Sadly, it also needs to be silent on the console so as not to
  // pollute -param output.
  if (commandLine_.argExists(std::string("-plugin")))
  {
    const std::string plugin = commandLine_.getArgumentValue("-plugin");
    
    for (std::string::size_type i = 0, j = plugin.find_first_of(", "); i != std::string::npos; i = (j == std::string::npos ? j : j + 1), j = plugin.find_first_of(", ", i))
    {
      Device::registerPlugin(plugin.substr(i, j).c_str());
    }
  }
  Report::safeBarrier(comm_);

  // Handle "-param", "-doc" and "-doc_cat"
  bool gotParamOpt=false;
  for (int paramArgSelect=0; paramArgSelect < 3; paramArgSelect++)
  {
    std::string paramArg=(paramArgSelect==0)?"-param":((paramArgSelect==1)?"-doc":"-doc_cat");
    if (commandLine_.argExists(paramArg))
    {
      if (Parallel::rank(comm_) == 0)
      {
        // Register the devices before printing out documentation
        Device::registerDevices();
        std::string deviceName = commandLine_.getArgumentValue(paramArg);
        // argExists returns false if the value is "", so we had to force
        // " " in there, now have to undo that.  
        if (deviceName == " ") deviceName = "";
        int modelLevel = commandLine_.getArgumentIntValue(paramArg+std::string("_level"),-1);
        int printFlags = (commandLine_.getArgumentIntValue(paramArg+std::string("_flags"), 3));
        bool printModel = ((printFlags & 2)==2);
        bool printInstance = ((printFlags & 1)==1);
        processParamOrDoc_(paramArg,deviceName,modelLevel,printModel,printInstance);
      }
      gotParamOpt=true;
    }
  }
  if (gotParamOpt)
  {
    return DONE;     // we do nothing further if we got one of those args.
  }
  
  Report::safeBarrier(comm_);
  runState_ = CHECK_NETLIST;
  // Don't bother doing anything else if we weren't given a netlist!
  const std::string &netlist_filename = commandLine_.getArgumentValue("netlist");
  if (netlist_filename.empty())
  {
    if (Parallel::rank(comm_) == 0)
      Xyce::lout() << "Netlist not found on command line" << std::endl
                   << "Usage: " << argv[0] << " [arguments] netlist" << std::endl
                   << "Use -h argument to get additional help info" << std::endl;
    return DONE;
  }

  if (Parallel::rank(comm_) == 0) 
  {
    // Error out if the user-specified netlist file does not exist, cannot be opened,
    // or is a directory name rather than a file name.  See SON Bugs 730 
    // and 785 for more details.  
    if ( !(Util::checkIfValidFile(netlist_filename)) )
    {
      Report::UserError() << "Could not open netlist file " << netlist_filename << " for reading.";
    }

    if ( commandLine_.argExists("-o") && !Util::checkIfValidDashoFileName(commandLine_.getArgumentValue("-o")) )
    {
      Report::UserError() << "Invalid basename " << commandLine_.getArgumentValue("-o") << " specified with -o";
    }
  }

  Report::safeBarrier(comm_);

  runState_ = OPEN_LOGSTREAM;
  initializeLogStream(Parallel::rank(comm_), Parallel::size(comm_));

  // Set the output stream of the "-l" flag exists
  if (commandLine_.argExists("-l"))
  {
    openLogFile(commandLine_.getArgumentValue("-l"), commandLine_.argExists("-per-processor"));
  }

  // Set the output stream of the "-l" flag exists
  if (commandLine_.argExists("-verbose"))
  {
    openDiagnosticFile(commandLine_.getArgumentValue("-verbose"), commandLine_.argExists("-per-processor"));
  }

  // Construct the date and time, in a way that all systems can
  // do right (Windows has Issues with more than 2 specifiers
  time_t now = time(NULL);
  char dayOfWeek[6],month[6],dayOfMonth[4],timeNow[20],timeZone[6], year[6];
  strftime(dayOfWeek,6,"%a",localtime(&now));
  strftime(month,6,"%b",localtime(&now));
  strftime(dayOfMonth,4,"%d",localtime(&now));
  strftime(timeNow,20,"%X",localtime(&now));
  strftime(timeZone,6,"%Z",localtime(&now));
  strftime(year,6,"%Y",localtime(&now));

  
  Xyce::lout() << "\n"
               << "*****\n"
               << "***** Welcome to the Xyce(TM) Parallel Electronic Simulator\n"
               << "*****\n"
               << "***** This is version " << Util::Version::getFullVersionString() << "\n"
               << "***** Date: " 
               << dayOfWeek << " " << month << " " << dayOfMonth << " "
               << timeNow << " " << timeZone<< " "<< year << "\n\n\n"
               << "***** Executing netlist " << netlist_filename << "\n\n";

  // check if a parameter file was specified for param substitution during parsing
  if (commandLine_.argExists(std::string("-prf")))
  {
    IO::readExternalParamsFromFile(*parallelManager_->getPDSComm(),
                                   commandLine_.getArgumentValue("-prf"), externalNetlistParams_);
  }

  if (commandLine_.argExists("-randseed"))
  {
    unsigned long theSeed;
    std::stringstream iss(commandLine_.getArgumentValue("-randseed"));
    iss >> theSeed;
    Xyce::Util::Expression::seedRandom((long)theSeed);
    Xyce::Device::DeviceSupport theDeviceSupport;
    theDeviceSupport.SetSeed((long)theSeed);
  }

  Report::safeBarrier(comm_);

  // Start the global timer.
  ElapsedTimerPtr_ = new Util::Timer();

  runState_ = ALLOCATE_SUBSYSTEMS;
  // Allocate all the various packages:
  doAllocations_();

  if (DEBUG_CIRCUIT)
    dout() << "Allocation was successful" << std::endl;

  // Register the external package pointers:
  if (!doRegistrations_())
    Report::DevelFatal() << "Registration was NOT successful" << std::endl;

  if (DEBUG_CIRCUIT)
    dout() << "Registration was successful" << std::endl;

  // ERK.  Allocate up the main expression group
  Teuchos::RCP<Util::mainXyceExpressionGroup>  mainExprGroup_ =  Teuchos::rcp(new Xyce::Util::mainXyceExpressionGroup (
    *parallelManager_->getPDSComm(),
    *topology_, *analysisManager_, *deviceManager_, *outputManager_));
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> baseGroupCast = mainExprGroup_;

  IO::NetlistImportTool netlist_import_tool(*opBuilderManager_, *parsingManager_, baseGroupCast);
  IO::registerPkgOptionsMgr(netlist_import_tool, *optionsManager_);

  // ERK. these could be in "doRegistrations" but the mainXyceGroup must be allocated after netlist_import_tool, which happens after ...
  bool bs1=deviceManager_->registerExpressionGroup(baseGroupCast); 
  bool bs2=analysisManager_->registerExpressionGroup(baseGroupCast);
  bool bs3=measureManager_->registerExpressionGroup(baseGroupCast);
  bool bs4=nonlinearManager_->registerExpressionGroup(baseGroupCast);

  runState_ = PARSE_NETLIST;
  Xyce::lout() << "***** Reading and parsing netlist..." << std::endl;

  {
    Stats::StatTop _netlistImportStat("Netlist Import");
    Stats::TimeBlock _netlistImportTimer(_netlistImportStat);

    netlist_import_tool.constructCircuitFromNetlist(
      commandLine_,
      hangingResistor_,
      netlist_filename,
      externalNetlistParams_,
      *topology_,
      *parallelManager_->getPDSComm(),
      *optionsManager_,
      *outputManager_,
      *deviceManager_,
      *measureManager_,
      *fourierManager_,
      *fftManager_);

    if (netlist_import_tool.getUseMOR())
      return DONE;

    Report::safeBarrier(comm_);

    if ( commandLine_.argExists("-syntax") || commandLine_.argExists("-count") )
    {
      if (commandLine_.argExists("-syntax"))
      {
        Xyce::lout() << "***** Netlist syntax OK\n";
      }

      Xyce::lout() << std::endl;

      Xyce::lout() << "***** Device Type Counts ...\n" << std::endl;

      IO::printDeviceCount(Xyce::lout(), deviceManager_->getDeviceCountMap());

      Xyce::lout() << std::endl;

      reportTotalElapsedTime ();

      return DONE;
    }

    mainExprGroup_->setAliasNodeMap(netlist_import_tool.getAliasNodeMap());
    outputManager_->setAliasNodeMap(netlist_import_tool.getAliasNodeMap());
    outputManager_->setMainContextFunctionMap(netlist_import_tool.getMainContextFunctions());
    outputManager_->setMainContextParamMap(netlist_import_tool.getMainContextParams().begin(), netlist_import_tool.getMainContextParams().end());
    outputManager_->setMainContextGlobalParamMap(netlist_import_tool.getMainContextGlobalParams().begin(), netlist_import_tool.getMainContextGlobalParams().end());

    runState_ = SETUP_TOPOLOGY;
    Xyce::lout() << "***** Setting up topology...\n" << std::endl;

    {
      Stats::StatTop _verifyStat("Verify Devices");
      Stats::TimeBlock _verifyTimer(_verifyStat);
      RunStatus run_status = setupTopology();
      if (run_status != SUCCESS)
        return run_status;
    }

    // if "-remeasure" was on the command line, then we don't need to
    // instantiate the devices.  We do need to partially allocate the
    // the analysisObject though for a DC remeasure.  That will be done 
    // later in the constructor for the RemeasureDC object.
    // 
    if (commandLine_.argExists("-remeasure"))
    {
      char analysisName = netlist_import_tool.getAnalysisName()[0];
      Parallel::Communicator &pds_comm = *parallelManager_->getPDSComm();
      measureManager_->remeasure(pds_comm, commandLine_.getArgumentValue("netlist"), commandLine_.getArgumentValue("-remeasure"), 
                                 analysisName, *opBuilderManager_, *outputManager_, *fftManager_, *analysisManager_,
                                 *analysisRegistry_, topology_->getNodeSymbols());

      Xyce::lout() << "***** Remeasure analysis complete\n" << std::endl;
      return DONE;
    }

    runState_ = INSTANTIATE_DEVICES;
    {
      Stats::StatTop _instantiateStat("Instantiate");
      Stats::TimeBlock _instantiateTimer(_instantiateStat);
      topology_->instantiateDevices();
    }

    IO::deferredParameterDiagnostics(comm_, netlist_import_tool.getDeferredParameterCheck(), *deviceManager_);

    Report::safeBarrier(comm_);

    deviceManager_->setGlobalFlags();
  }

  Report::safeBarrier(comm_);

  return SUCCESS;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::initializeLate
// Purpose       : capture late "initialization-type" activities in one
//                 method
// Special Notes : This is basically the second half of what used to be
//                 the initialize function.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 02/28/2017
//-----------------------------------------------------------------------------
///
/// Perform final initialization steps, after the external simulator has
/// had the opportunity to modify the number of internal variables of the
/// interface devices
///
Simulator::RunStatus Simulator::initializeLate()
{
  Stats::StatTop _lateInitStat("Late Initialization");
  Stats::TimeBlock _lateInitTimer(_lateInitStat);
  {
    // Setup of indices including global reordering.
    {
      Stats::StatTop _globalIndexStat("Global Indices");
      Stats::TimeBlock _globalIndexTimer(_globalIndexStat);
      topology_->setupGlobalIndices();
    }
  }

  Report::safeBarrier(comm_);

  const time_t now = time(NULL);
  char timeDate[40];

//For whatever reason, Windows only allows two format specifiers, otherwise it won't return the time and date
#ifdef HAVE_WINDOWS_H
    strftime(timeDate,40,"%x %X",localtime(&now));
#else
    strftime(timeDate,40,"%x %X %Z",localtime(&now));
#endif

  Xyce::lout() << "***** Device Count Summary ..." << std::endl;
  {
    IO::DeviceCountMap global_device_count_map;

    IO::gatherGlobalDeviceCount(comm_, global_device_count_map, deviceManager_->getDeviceCountMap());

    auditJSON_ << Util::nameValuePair("StartTime", timeDate) << Util::JSON::sep
               << Util::nameValuePair("DeviceCount", global_device_count_map) << Util::JSON::sep
               << Util::nameValuePair("DeviceTotalCount", std::accumulate(global_device_count_map.begin(), global_device_count_map.end(), 0, IO::DeviceCountMapSum()));

    IO::printDeviceCount(Xyce::lout(), global_device_count_map);
  }

  if (commandLine_.argExists("-norun") || commandLine_.argExists("-namesfile") || 
      commandLine_.argExists("-noise_names_file") )
  {
    Xyce::lout() << "\n***** Syntax and topology analysis complete" << std::endl;

    // Want to allow the user to use both -namesfile and -noise_names_file on a 
    // Xyce command line.  However, can only call setUpMatrixStructure_() once.  
    bool matrixSetup = false;
 
    if (commandLine_.argExists("-namesfile"))
    {
      const std::string &namesfile = commandLine_.getArgumentValue("-namesfile");

      matrixSetup = setUpMatrixStructure_();  // matrixSetup is set to true if this succeeds
      topology_->outputNameFile(comm_, namesfile, true);
    }

    if (commandLine_.argExists("-noise_names_file"))
    {
      const std::string &noisenamesfile = commandLine_.getArgumentValue("-noise_names_file");

      if (!matrixSetup) setUpMatrixStructure_();
      Analysis::outputNoiseNameFile(comm_, noisenamesfile, *circuitLoader_);
    }
   
    _lateInitTimer.stop();
    rootStat_.stop();

    if (Parallel::rank(comm_) == 0 && Parallel::size(comm_) > 1) {
      pout() << std::endl
             << "Timing summary of processor " << Parallel::rank(comm_) << std::endl;
      Stats::printStatsTable(pout(), rootStat_, Stats::METRICS_ALL, false);
    }

    Xyce::lout() << std::endl
                 << "Timing summary of " << Parallel::size(comm_) << " processor" << (Parallel::size(comm_) == 1 ? "" : "s") << std::endl;
    Stats::printStatsTable(Xyce::lout(), rootStat_, Stats::METRICS_ALL, false, comm_);

    reportTotalElapsedTime ();

    return DONE;
  }

  runState_ = SETUP_MATRIX_STRUCTURE;
  Xyce::lout() << "\n***** Setting up matrix structure..." << std::endl;
  if (!setUpMatrixStructure_())
    return ERROR;

  runState_ = INITIALIZE_SYSTEM;
  Xyce::lout() << "***** Initializing...\n" << std::endl;
  if (!doInitializations_())
    return ERROR;

  // optional diagnostic output file:
  const std::string &netlist_filename = commandLine_.getArgumentValue("netlist");
  topology_->outputNameFile(comm_, "namesMap_" + netlist_filename + ".txt");

  outputManager_->checkPrintParameters(comm_, *opBuilderManager_);
  fourierManager_->fixupFourierParameters(comm_, *opBuilderManager_);
  fftManager_->enableFFTAnalysis(analysisManager_->getAnalysisMode());
  fftManager_->fixupFFTParameters(comm_, *outputManager_, *opBuilderManager_, analysisManager_->getFinalTime(),
                                  analysisManager_->getStepErrorControl());

  // Make the measure operators, and then check for agreement between the 
  // measures' requested modes (e.g, TRAN) and the analysis type (e.g, .TRAN).
  // It is a fatal error, if they disagree.  At this point, we can also set 
  // the suffix (e.g., .mt versus .ms versus .ma) for the .MEASURE output files.
  // Note that a different code path, contained in the measureManager code, is 
  // used for these three steps for when remeasure is used.
  measureManager_->makeMeasureOps(comm_, *opBuilderManager_);
  if ( !measureManager_->checkMeasureModes(analysisManager_->getAnalysisMode()) )
  {
    return ERROR;
  }

  measureManager_->fixupFFTMeasures(comm_, *fftManager_);
  measureManager_->setMeasureOutputFileSuffix(analysisManager_->getAnalysisMode());

  // if we loaded a parameter file from the command line, then the output manager
  // should scan the results as well as such files can specify response functions
  // than need to be reported by the output manager.
  if (commandLine_.argExists(std::string("-prf")))
  {
    outputResponse_->setExternalNetlistParams( externalNetlistParams_ );
  }

  if (commandLine_.argExists("-rsf"))
  {
    std::string responseFilename = commandLine_.getArgumentValue("-rsf");

    // Error out, before starting the simulation, if the user-specified response file
    // can not be opened.  See SON Bug 785 for more details.  
    if (Parallel::rank(comm_) == 0)
    {
      std::ofstream os(responseFilename.c_str());
      if ( !os.good() )
      {
        Report::UserError() << "Could not open response file " << responseFilename << " for output";
      } 
      os.close();
    }
    Report::safeBarrier(comm_);
    
    // response file is valid, so continue onwards
    outputResponse_->setResponseFilename(responseFilename);
  }

  // initialize all the outputters and other objects (like sensitivity objective 
  // functions) which depend on expressions.   Expressions can depend on .global_params,  
  // .params as well as .funcs, and this is the last chance to resolve them before the 
  // maps get deleted.
  //
  // Also, the outputters must be set up after "allocateAnalysisObject", 
  // which is called from "doInitializations".  
  analysisManager_->finalExpressionBasedSetup();

  // at this point, since setup is complete, the maps for 
  // .global_param, .param and .func can be deleted.
  outputManager_->deleteMainContextFunctionMap();
  outputManager_->deleteMainContextParamMap();
  outputManager_->deleteMainContextGlobalParamMap();

  // Start the Xyce solver timer, now that the setup is complete.
  XyceTimerPtr_ = new Util::Timer();

  return SUCCESS;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::getDeviceNames
// Purpose       : get all names of devices of specified type (modelGroupName) 
//                 in the netlist
// Special Notes : "modelGroupName" takes a string of the form the devices would
//                 have when instantiated in a netlist, e.g. "R" for resistors or 
//                 "YGENEXT" for GenExt devices.  For U devices, the modelGroupName would
//                 be "BUF" and not UBUF".  The returned device name(s) will be the 
//                 fully qualified name(s), including any subcircuit hierarchy.
//                 For a YADC device named ADC1, it would be YADC!ADC1.
//                 This function will return false if the requested modelGroupName
//                 is invalid.  It will also return false if the netlist does not
//                 contain any devices for a valid modelGroupName.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/25/08
//-----------------------------------------------------------------------------
bool Simulator::getDeviceNames(const std::string &modelGroupName, std::vector<std::string> &deviceNames)
{
  // this call will succeed for any valid non-Y device, including U devices
  Device::EntityTypeId model_group = deviceManager_->getModelGroup(modelGroupName);

  // try again, assuming it's a Y device.  All Y device names must have at least two characters        
  if ( !model_group.defined() && (modelGroupName.length() > 1) && (modelGroupName[0] == 'Y') )
  {
    model_group = deviceManager_->getModelGroup(Device::InstanceName(modelGroupName).getDeviceType());
  }

  if (!model_group.defined())
  {
    Report::UserWarning0() << "No devices from model group " << modelGroupName << " found in netlist";
    return false;
  }

  Device::Device *device = deviceManager_->getDevice(model_group);

  if (!device)
  {
    Report::UserWarning0() << "No devices from model group " << modelGroupName << " found in netlist";
    return false;
  }

  Device::getDeviceInstanceNames(*device, std::back_inserter(deviceNames));
  
  return true;
}

//----------------------------------------------------------------------------
// Function       : Simulator::getDACDeviceNames
// Purpose        : Gets the fully qualified names of the DAC devices, including
//                  any subcircuit hierarchy.  For a YDAC device named DAC1, it 
//                  would be YDAC!DAC1.
// Special Notes  :
// Scope          :
// Creator        : Lisa Maynes
// Creation Date  : 06/13/2003
//----------------------------------------------------------------------------
bool Simulator::getDACDeviceNames(std::vector< std::string >& dacNames)
{
  dacNames.clear();

  Device::Device *device = deviceManager_->getDevice(Device::DAC::Traits::modelGroup());
  if (!device)
  {
    Report::UserWarning0() << "No DAC devices found in netlist";
    return false;
  }

  Device::getDeviceInstanceNames(*device, std::back_inserter(dacNames));

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::getAllDeviceNames
// Purpose       : get the names of all devices in the netlist.
// Special Notes : The returned device name(s) will be the fully qualified
//                 name(s), including any subcircuit hierarchy.  This function
//                 will return false if the netlist does not have any devices.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 12/11/2019
//-----------------------------------------------------------------------------
bool Simulator::getAllDeviceNames(std::vector<std::string> &deviceNames)
{
  // This map will only be populated with device types that exist in the netlist
  Device::InstanceVector instance_ptr_vec = deviceManager_->getInstancePtrVec();

  if (instance_ptr_vec.size() == 0)
  {
    Report::UserWarning0() << "No devices found in netlist";
    return false;
  }

  Device::InstanceVector::const_iterator it = instance_ptr_vec.begin();
  Device::InstanceVector::const_iterator end = instance_ptr_vec.end();
  for ( ; it != end; ++it)
  {
    deviceNames.push_back((*it)->getName().getEncodedName());
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::checkDeviceParamName
// Purpose       : check if the specified parameter name (e.g., X1:R1:R) is
//                 a valid parameter for a device that exists in the netlist.
// Special Notes : This function will return false if the device or parameter
//                 does not exist.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 12/18/2019
//-----------------------------------------------------------------------------
bool Simulator::checkDeviceParamName(const std::string full_param_name) const
{
  Device::DeviceEntity *device_entity = deviceManager_->getDeviceEntity(full_param_name);
  if (!device_entity)
  {
    Report::UserWarning0() << "Device entity not found for " << full_param_name << std::endl;
    return false;
  }

  if ( !device_entity->findParam(Util::paramNameFromFullParamName(full_param_name)) )
  {
    Report::UserWarning0() << "Device parameter not found for " << full_param_name << std::endl;
    return false;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::getDeviceParamVal
// Purpose       : get the value of a specified device parameter.
// Special Notes : This function will return false if the device or parameter
//                 does not exist.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 12/11/2019
//-----------------------------------------------------------------------------
bool Simulator::getDeviceParamVal(const std::string full_param_name, double& val) const
{
  Device::DeviceEntity *device_entity = deviceManager_->getDeviceEntity(full_param_name);
  if (!device_entity)
  {
    Report::UserWarning0() << "Device entity not found for " << full_param_name << std::endl;
    return false;
  }

  if ( !device_entity->getParam(Util::paramNameFromFullParamName(full_param_name), val) )
  {
    Report::UserWarning0() << "Device parameter not found for " << full_param_name << std::endl;
    return false;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::getNumAdjNodesForDevice
// Purpose       : get the number of the nodes that are adjacent to a specified
//                 device, including the ground node
// Special Notes : This function will return false if the device does not exist.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/7/2020
//-----------------------------------------------------------------------------
bool Simulator::getNumAdjNodesForDevice(const std::string deviceName, int& numAdjNodes) const
{
  bool retVal = true;
  ExtendedString deviceNameUC(deviceName);
  deviceNameUC.toUpper();
  const Topo::CktNode * cNodePtr = topology_->findCktNode(NodeID(deviceNameUC, Xyce::_DNODE));

  if (cNodePtr == 0)
  {
    Report::UserWarning0() << "Device " << deviceName << " not found" << std::endl;
    numAdjNodes = 0;
    retVal = false;
  }
  else
  {
    numAdjNodes = topology_->numAdjNodesWithGround(cNodePtr->get_gID());
  }

  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::getAdjGIDsForDevice
// Purpose       : get the GIDs of the nodes that are adjacent to a specified
//                 device, including the ground node (GID of -1).  The other
//                 GIDs are non-negative integers.
// Special Notes : This function will return false if the device does not exist.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/7/2020
//-----------------------------------------------------------------------------
bool Simulator::getAdjGIDsForDevice(const std::string deviceName, std::vector<int> & adj_GIDs) const
{
  bool retVal = true;
  ExtendedString deviceNameUC(deviceName);
  deviceNameUC.toUpper();
  const Topo::CktNode * cNodePtr = topology_->findCktNode(NodeID(deviceNameUC, Xyce::_DNODE));

  if (cNodePtr == 0)
  {
    Report::UserWarning0() << "Device " << deviceName << " not found" << std::endl;
    retVal = false;
  }
  else
  {
    topology_->returnAdjGIDsWithGround(cNodePtr->get_gID(), adj_GIDs);
  }

  return retVal;
}

//----------------------------------------------------------------------------
// Function       : getADCMap
// Purpose        : Gets the names of the ADC devices in the circuit (as the key 
//                  of the map) with an inner map of their instance parameters 
//                  (keyed by parameter name) for each device.  For a YADC 
//                  device named ADC1, the returned name would be YADC!ADC1.
// Special Notes  :
// Scope          :
// Creator        : Tom Russo
// Creation Date  : 05/07/2004
//----------------------------------------------------------------------------
bool Simulator::getADCMap(std::map<std::string, std::map<std::string, double> >&ADCMap)
{
  ADCDeviceInstanceParameterOp op(ADCMap);

  Device::Device *device = deviceManager_->getDevice(Device::ADC::Traits::modelGroup());
  if (device)
  {
    device->forEachInstance(op);
    return true;
  }
  else
  {
    return false;
  }
}

//----------------------------------------------------------------------------
// Function       : updateTimeVoltagePairs
// Purpose        : Update the DAC devices in a circuit by adding the set
//                  of time and voltage pairs built up on the "digital side"
//                  since the last update and by removing the time-voltage
//                  pairs for times that pre-date the given simulation time.
// Special Notes  : The current method for locating DAC instances
//                  works for the serial case only. The parallel
//                  case will be added in later modifications.
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 06/09/2003
//----------------------------------------------------------------------------
bool Simulator::updateTimeVoltagePairs(const std::map< std::string, std::vector<std::pair<double,double> > *> & timeVoltageUpdateMap)
{
  // assume success, but then return false if the update fails at any DAC
  bool success = true;

  for (std::map<std::string, std::vector< std::pair<double,double> >* >::const_iterator it = timeVoltageUpdateMap.begin(), end = timeVoltageUpdateMap.end(); it != end; ++it) 
  {
    const std::string &dacName = (*it).first;
    const std::vector<std::pair<double,double> > &tv_pair_vector = *(*it).second;

    Device::DAC::Instance *dacInstancePtr = getDACInstance_(dacName);

    if (dacInstancePtr)
    {
      // Update the time-voltage pairs for the given DAC instance.
      if (!dacInstancePtr->updateTVVEC(tv_pair_vector))
      {
        Report::UserWarning0() << "Failed to update the time-voltage pairs for the DAC " << dacName;
        success = false;
      }
    }
    else
    {
      Report::UserWarning0() << "Failed to update the time-voltage pairs for the DAC " << dacName;
      success = false;
    }
  }

  return success;
}

//----------------------------------------------------------------------------
// Function       : getTimeVoltagePairs
// Purpose        : get a map of all time-voltage pairs from all ADC instances
//
// Special Notes  : The current method for locating ADC instances
//                  works for the serial case only. The parallel
//                  case will be added in later modifications.
// Scope          : public
// Creator        : Tom Russo
// Creation Date  : 05/10/2004
//----------------------------------------------------------------------------
bool Simulator::getTimeVoltagePairs(std::map< std::string, std::vector< std::pair<double,double> > > & timeVoltageUpdateMap)
{
  bool success = false;

  Device::Device *device = deviceManager_->getDevice(Device::ADC::Traits::modelGroup());
  if (device) 
  {
    TimeVoltagePairsOp op(timeVoltageUpdateMap);
    timeVoltageUpdateMap.clear();

    device->forEachInstance(op);
    success = true;
  }

  return success;
}

//----------------------------------------------------------------------------
// Function       : getTimeStatePairs
// Purpose        : get a map of all time-state pairs from all ADC instances
//
// Special Notes  : The current method for locating ADC instances
//                  works for the serial case only. The parallel
//                  case will be added in later modifications.
// Scope          : public
// Creator        : Pete Sholander
// Creation Date  : 11/13/2018
//----------------------------------------------------------------------------
bool Simulator::getTimeStatePairs(std::map< std::string, std::vector< std::pair<double,int> > > & timeStateUpdateMap)
{
  bool success = false;

  Device::Device *device = deviceManager_->getDevice(Device::ADC::Traits::modelGroup());
  if (device) 
  {
    TimeStatePairsOp op(timeStateUpdateMap);
    timeStateUpdateMap.clear();

    device->forEachInstance(op);
    success = true;
  }

  return success;
}

//----------------------------------------------------------------------------
// Function       : setADCWidths
// Purpose        : Update the ADC devices in a circuit by informing them
//                  of the width of their bitvector output on the
//                  "digital side"
// Special Notes  :
// Scope          :
// Creator        : Tom Russo
// Creation Date  : 05/07/2004
//----------------------------------------------------------------------------
bool Simulator::setADCWidths(const std::map<std::string, int> &ADCWidthMap)
{
  // handle pathological case
  if (ADCWidthMap.size() == 0) 
  {
    Report::UserWarning0() << "setADCWidths() called with empty list of ADC names";
    return false;
  }

  // assume success, but then return false if the width is not set for any ADC
  bool success = true;

  for (std::map<std::string, int>::const_iterator it = ADCWidthMap.begin(), end = ADCWidthMap.end(); it != end; ++it) 
  {
    const std::string &adcName = (*it).first;
    const int width = (*it).second;

    Device::ADC::Instance *adcInstancePtr = getADCInstance_(adcName);

    if (adcInstancePtr) 
    {
      // Update output bit-vector width for this ADC.
      if (!adcInstancePtr->setBitVectorWidth(width))
      {
        Report::UserWarning0() << "Failed to update the width for ADC " << adcName;
        success = false;
      }     
    }
    else
    {
      Report::UserWarning0() << "Failed to update the width for ADC " << adcName;
      success = false;
    }
  }

  return success;
}

//----------------------------------------------------------------------------
// Function       : getADCWidths
// Purpose        : get the width of the specified ADC devices bitvector output 
//                  on the "digital side"
// Special Notes  :
// Scope          :
// Creator        : Pete Sholander
// Creation Date  : 11/16/2018
//----------------------------------------------------------------------------
bool Simulator::getADCWidths(std::map<std::string, int> &ADCWidthMap)
{
  // handle pathological case
  if (ADCWidthMap.size() == 0) 
  {
    Report::UserWarning0() << "getADCWidths() called with empty list of ADC names";
    return false;
  }

  // assume success, but then return false if the width is not obtained for any ADC
  bool success = true;
  for (std::map<std::string, int>::iterator it = ADCWidthMap.begin(), end = ADCWidthMap.end(); it != end; ++it)
  {
    const std::string &adcName = (*it).first;
    Device::ADC::Instance *adcInstancePtr = getADCInstance_(adcName);

    if (adcInstancePtr) 
    {
      // get the bit-width for the given ADC instance.
      if (adcInstancePtr->getBitVectorWidth())
      {
        (*it).second = adcInstancePtr->getBitVectorWidth();
      }
      else
      {
        Report::UserWarning0() << "Failed to get the width for ADC " << adcName;
        success = false;
      }     
    }
    else
    {
      Report::UserWarning0() << "Failed to get the width for ADC " << adcName;
      success = false;
    }
  }

  return success;
}

//---------------------------------------------------------------------------
// Function      : Simulator::simulateUntil
// Purpose       : To continue the existing analog circuit simulation
//                 until either the given <requestedUntilTime> is reached
//                 or the simulation termination criterion is met.
//                 Return a Boolean indicating whether the simulation
//                 run was successful. (Note that the run is successful
//                 even when the given <requestedUntilTime> is not reached,
//                 so long as the run completed normally.)
// Special Notes : The time variables are in units of seconds.
// Scope         : public
// Creator       : Lon Waters, SNL
//               : Tom Russo, SNL, Component Information and Models
// Creation Date : 06/03/2003
//---------------------------------------------------------------------------
bool Simulator::simulateUntil(
  double        requestedUntilTime,
  double &      completedUntilTime)
{
  bool bsuccess = false;
  double currentTimeBeforeSim = analysisManager_->getTime();
  double finalTime = analysisManager_->getFinalTime();
  double initialTime = analysisManager_->getInitialTime();

  if (requestedUntilTime <= currentTimeBeforeSim)
    Report::UserError0() << "requestedUntilTime <= current simulation time in simulateUntil() call.  Simulation will abort.";

  // Silence the "Percent Complete" noise, we don't want it when using this
  // interface
  analysisManager_->silenceProgress();

  if (DEBUG_CIRCUIT)
    dout() << "simulateUntil: finalTime = " << finalTime<< ", currentTimeBeforeSim = " << currentTimeBeforeSim << std::endl;

  if (currentTimeBeforeSim >= finalTime)
  {
    // We have already simulated as far as the netlist requested.
    bsuccess = true;
    completedUntilTime = currentTimeBeforeSim;

    if (DEBUG_CIRCUIT)
      dout() << "Case1: completedUntilTime = " << completedUntilTime;
  }
  else
  {
    // We are not already past the end of the netlist time
    analysisManager_->setPauseTime(std::min(requestedUntilTime, finalTime), analysisManager_->getTIAParams().initialTime);

    if (DEBUG_CIRCUIT)
      dout() << "simulateUntil currentTimeBeforeSim = " << currentTimeBeforeSim << "  initialTime = " << initialTime << std::endl;

    if (currentTimeBeforeSim > initialTime)
    {
      analysisManager_->setResumeSimulation(true);
    }

    if (DEBUG_CIRCUIT)
      dout() << "simulateUntil: Case2: requestedUntilTime = " << requestedUntilTime
             << ", pauseTime = " << analysisManager_->getPauseTime() << std::endl;

    bsuccess = runSolvers_();
    completedUntilTime = analysisManager_->getTime();
    if (DEBUG_CIRCUIT)
      dout() << "simulateUntil: Case2: completedUntilTime = " << completedUntilTime << std::endl;
  }

  if (DEBUG_CIRCUIT)
    dout() << std::endl;

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::finalize
// Purpose       : To clean up after driving Xyce with the SIMBUS
//                 simulation backplane. This includes the following:
//                    Free any dynamically allocated memory...
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
//               : Lisa Maynes, CoMeT Solutions, Inc.
// Creation Date : 06/03/2003
//---------------------------------------------------------------------------
Simulator::RunStatus Simulator::finalize()
{
  Xyce::lout() << "\n***** Solution Summary *****"  << std::endl;

  analysisManager_->printLoopInfo(0, 0);

  // The Stats will replace this bit of ugliness.  But for now I'll write a function to just get the value.
  Analysis::StatCounts analysis_stat_counts = analysisManager_->getAnalysisObject().getStatCounts() - analysisManager_->getAnalysisObject().getStatCounts(0);

  IO::outputMacroResults(comm_,
                         *measureManager_,
                         *fourierManager_,
                         *fftManager_,
                         outputManager_->getNetlistFilename(),
                         outputResponse_->getResponseFunctions(),
                         outputResponse_->getResponseFilename(),
                         outputManager_->getStepLoopNumber(),
                         analysisManager_->getFinalTime());

  rootStat_.stop();

  Xyce::lout() << std::endl
               << "***** Total Simulation Solvers Run Time: " << XyceTimerPtr_->elapsedTime() << " seconds" << std::endl
               << "***** Total Elapsed Run Time:            " << ElapsedTimerPtr_->elapsedTime() << " seconds" << std::endl
               << "*****" << std::endl
               << "***** End of Xyce(TM) Simulation" << std::endl
               << "*****" << std::endl;

  const char *xyce_no_tracking = ::getenv("XYCE_NO_TRACKING");
  if (Parallel::rank(comm_) == 0 && trackingURL && !xyce_no_tracking)
  {
    const time_t now=time(NULL);
    char timeDate[40];

//For whatever reason, Windows only allows two format specifiers, otherwise it won't return the time and date
#ifdef HAVE_WINDOWS_H
    strftime(timeDate,40,"%x %X",localtime(&now));
#else
    strftime(timeDate,40,"%x %X %Z",localtime(&now));
#endif

    auditJSON_ << Util::JSON::sep
               << Util::nameValuePair("Hostname", hostname()) << Util::JSON::sep
               << Util::nameValuePair("Domainname", domainname()) << Util::JSON::sep
               << Util::nameValuePair("Username", username()) << Util::JSON::sep
               << Util::nameValuePair("Hardware", hardware()) << Util::JSON::sep
               << Util::nameValuePair("OSname", osname()) << Util::JSON::sep
               << Util::nameValuePair("OSversion", osversion()) << Util::JSON::sep
               << Util::nameValuePair("Version", Util::Version::getFullVersionString()) << Util::JSON::sep
               << Util::nameValuePair("Processors", Parallel::size(comm_)) << Util::JSON::sep
               << Util::nameValuePair("PrimaryAnalysis", Analysis::analysisModeName(analysisManager_->getAnalysisMode())) << Util::JSON::sep
               << Util::nameValuePair("EndTime", timeDate) << Util::JSON::sep
               << Util::nameValuePair("RuntimeStats", analysis_stat_counts);
    Util::JSON json;
    json << Util::nameValuePair("audit", auditJSON_);

    const std::string &json_string = json.str();

    Util::sendTrackingData(trackingURL, 0, json_string);
  }

  if (Parallel::rank(comm_) == 0 && Parallel::size(comm_) > 1) {
    pout() << std::endl
           << "Timing summary of processor " << Parallel::rank(comm_) << std::endl;
    Stats::printStatsTable(pout(), rootStat_, Stats::METRICS_ALL, false);
  }

  Xyce::lout() << std::endl
               << "Timing summary of " << Parallel::size(comm_) << " processor" << (Parallel::size(comm_) == 1 ? "" : "s") << std::endl;
  Stats::printStatsTable(Xyce::lout(), rootStat_, Stats::METRICS_ALL, false, comm_);

  // Close the output stream:
  closeLogFile();

  return SUCCESS;
}

//---------------------------------------------------------------------------
// Function      : Simulator::reportTotalElapsedTime ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/01/2009
//---------------------------------------------------------------------------
void Simulator::reportTotalElapsedTime()
{
  Xyce::lout() <<  "\n***** Total Elapsed Run Time: " << ElapsedTimerPtr_->elapsedTime() << " seconds" << std::endl;
}

//---------------------------------------------------------------------------
// Function      : Simulator::simulationComplete
// Purpose       : Simply report whether we've reached the end of the
//                 simulation
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 05/06/2004
//---------------------------------------------------------------------------
bool Simulator::simulationComplete()
{
  return analysisManager_->isSimulationComplete();
}

//---------------------------------------------------------------------------
// Function      : Simulator::checkResponseVar
// Purpose       :
// Special Notes : Used when Dakota or an external program calls Xyce to tell
//                 Xyce what .measure lines are to be used as response functions
//                 to pass back to Dakota.  This call checks the measure manager
//                 in the I/O package has set up measure objects for each label.
// Scope         : public
// Creator       : Heidi Thornquist, Sandia National Labs
// Creation Date : 02/24/2014
//---------------------------------------------------------------------------
bool Simulator::checkResponseVar(
  const std::string &   variable_name) const
{
  if (!measureManager_)
    Report::DevelFatal0().in("Simulator::checkResponseVar") << "measureManager_ is null";

  return measureManager_->find(variable_name);
}

//---------------------------------------------------------------------------
// Function      : Simulator::obtainResponse
// Purpose       :
// Special Notes : Used when Dakota or an external program calls Xyce to tell
//                 Xyce what .measure lines are to be used as response functions
//                 to pass back to Dakota.  This call obtains the responses for
//                 each labelled measure.
// Scope         : public
// Creator       : Heidi Thornquist, Sandia National Labs
// Creation Date : 02/24/2014
//---------------------------------------------------------------------------
bool Simulator::obtainResponse(
  const std::string &   variable_name,
  double &              result) const
{
  if (!measureManager_)
    Report::DevelFatal0().in("Simulator::obtainResponse") << "measureManager_ is null";

  return measureManager_->getMeasureValue(variable_name, result);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getDACInstancePtr_
// Purpose       : Returns the pointer to a named DAC device instance
// Special Notes :
// Scope         : PRIVATE
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
Device::DAC::Instance *
Simulator::getDACInstance_(const std::string &deviceName)
{
  // See if we've looked this up before.
  if (dacDeviceMap_.empty())
  {
    Device::Device *device = deviceManager_->getDevice(Device::DAC::Traits::modelGroup());
    if (device)
      Device::mapDeviceInstances(*device, dacDeviceMap_);
  }

  std::map<std::string, Device::DAC::Instance *>::iterator mapIter = dacDeviceMap_.find(deviceName);
  if (mapIter == dacDeviceMap_.end())
    return 0;

  return (*mapIter).second;
}


//-----------------------------------------------------------------------------
// Function      : Simulator::getADCInstance_
// Purpose       : Returns the pointer to a named ADC device instance
// Special Notes :
// Scope         : PRIVATE
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
Device::ADC::Instance *
Simulator::getADCInstance_(const std::string &deviceName)
{
  // See if we've looked this up before.
  if (adcDeviceMap_.empty())
  {
    Device::Device *device = deviceManager_->getDevice(Device::ADC::Traits::modelGroup());
    if (device)
      Device::mapDeviceInstances(*device, adcDeviceMap_);
  }

  std::map<std::string, Device::ADC::Instance *>::iterator mapIter = adcDeviceMap_.find(deviceName);
  if (mapIter == adcDeviceMap_.end())
    return 0;

  return (*mapIter).second;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::processParamOrDoc_
// Purpose       : Handle "-param", "-doc" or "-doc_cat" options
// Special Notes :
// Scope         : PRIVATE
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/09/16
//-----------------------------------------------------------------------------
void Simulator::processParamOrDoc_(std::string & optionName, std::string & deviceName,
                                  int modelLevel, bool printModel, bool printInstance)
{
  Device::OutputMode::Mode format = Device::OutputMode::DEFAULT;
  if (optionName == "-param")
    format = Device::OutputMode::PARAM;
  else if (optionName == "-doc")
    format = Device::OutputMode::DOC;
  else if (optionName == "-doc_cat")
    format = Device::OutputMode::DOC_CAT;

  Device::handleParameterOutputs(format, deviceName, modelLevel, printModel,
                                 printInstance);


}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void report_handler(
  const char *  message,
  unsigned      report_mask)
{

  std::ostringstream oss;
  Util::word_wrap(oss, message, s_errorWrap, " ", "");
  bool asymmetric;
  
  // If symmetric then all processors are getting the same message,
  // only write to p0 and not to ~p0 backlog.
  // If asymetric then one processor is getting the message, write to
  // per processor stream which writes to per processor log file and
  // to backlog.
  if (report_mask & Report::MSG_SYMMETRIC)
  {
    Xyce::lout() << oss.str();
    asymmetric=false;
  }
  else
  {
    pout() << oss.str();
    asymmetric=true;
  }

  // If fatal error also send the message to the standard error file, but
  // make sure we only do this from proc0 if a symmetric error in parallel.
  if (report_mask & Report::MSG_TERMINATE)
  {
    bool outputToErr=true;

#ifdef Xyce_PARALLEL_MPI
    // If we're a symmetric error, only output to stderr on proc0
    // If we're asymmetric, output to stderr from wherever we are.
    if (!asymmetric)
    {
      int me;
      MPI_Comm_rank(MPI_COMM_WORLD,&me);
      if (me != 0)
        outputToErr=false;
    }
#endif

    // don't send the message itself to lout, because we've already done that
    // above.  
    Xyce::lout() << "*** Xyce Abort ***" << std::endl;
    if (outputToErr)
      std::cerr << oss.str() << std::endl << std::endl
                << "*** Xyce Abort ***" << std::endl;
      
    Xyce_exit(1,asymmetric);
  }
}

//-----------------------------------------------------------------------------
// Function      : Xyce_exit
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 07/26/04
//-----------------------------------------------------------------------------
void Xyce_exit(int code,bool asymmetric)
{
#ifdef Xyce_PARALLEL_MPI
  int me;
  if (!asymmetric)
    MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  if ((me == 0 || asymmetric) && code != 0)
    MPI_Abort(MPI_COMM_WORLD,code);
  MPI_Finalize();
#else
  exit(code);
#endif
}

} // namespace Circuit
} // namespace Xyce

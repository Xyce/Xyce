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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DeviceMgr_h
#define Xyce_N_DEV_DeviceMgr_h

#if defined(HAVE_UNORDERED_MAP)
#include <unordered_map>
using std::unordered_map;
#elif defined(HAVE_TR1_UNORDERED_MAP)
#include <tr1/unordered_map>
using std::tr1::unordered_map;
#else
#error neither unordered_map or tr1/unordered_map found
#endif
#include <map>
#include <set>
#include <string>
#include <vector>

#include <N_ANP_fwd.h>
#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_NLS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TIA_fwd.h>
#include <N_TOP_fwd.h>
#include <N_UTL_fwd.h>

#include <N_ANP_StepEvent.h>
#include <N_LOA_Loader.h>
#include <N_DEV_ArtificialParameters.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceSensitivities.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_UTL_Listener.h>
#include <N_UTL_Op.h>

class expressionGroup;

namespace Xyce {
namespace Device {

bool
getParamAndReduce(
  Parallel::Machine     comm,
  const DeviceMgr &     device_manager,
  const std::string &   name,
  double &              val);

//-----------------------------------------------------------------------------
// Class         : DeviceMgr
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class DeviceMgr : public Util::Listener<Analysis::StepEvent>
{
  friend struct ArtificialParameters::ArtificialParameter;
  friend bool getParamAndReduce(Parallel::Machine comm, const DeviceMgr &device_manager, const std::string &name, double &value);

public:
  typedef std::vector<Device *> DeviceVector;
  typedef std::vector<DeviceEntity *> EntityVector;
  typedef std::vector<DeviceInstance *> InstanceVector;
  typedef std::vector<DeviceModel *> ModelVector;
  typedef std::map<ModelTypeId, ModelVector> ModelTypeModelVectorMap;
  typedef std::map<ModelTypeId, InstanceVector> ModelTypeInstanceVectorMap;
  typedef unordered_map<std::string, DeviceEntity *, HashNoCase, EqualNoCase> DeviceEntityMap;
  typedef unordered_map<std::string, ModelTypeId, HashNoCase, EqualNoCase> ModelTypeNameModelTypeIdMap;

  typedef unordered_map<std::string, SourceInstance *, HashNoCase, EqualNoCase> IndependentSourceMap;
  typedef std::vector<SourceInstance *> IndependentSourceVector;

  typedef unordered_map<std::string, Util::Op::Operator *> OpMap;
  
  DeviceMgr(Parallel::Machine comm, Topo::Topology &topology, Util::Op::BuilderManager &op_builder_manager, const IO::CmdParse &command_line);
  ~DeviceMgr();

private:
  DeviceMgr(const DeviceMgr &);                       ///< No copying
  DeviceMgr &operator=(const DeviceMgr &);            ///< No assignment

public:
  void notify(const Analysis::StepEvent &event);

  bool registerAnalysisManager(Analysis::AnalysisManager *analysis_manager);

  bool registerExpressionGroup(Teuchos::RCP<Xyce::Util::baseExpressionGroup> & group);

  bool registerNonlinearSolver (Nonlinear::Manager * tmp_nlsMgrPtr)
  {
    nlsMgrPtr_              = tmp_nlsMgrPtr;

    return nlsMgrPtr_ != 0;
  }

  // this function is called from the output manager (through the
  // device interface) to inform the device package of the devices for
  // which lead currents have been requested. The device manager will
  // take care of doing isolated F and Q loads for these devices so
  // the lead currents can be calculated
  bool setLeadCurrentRequests(const std::set<std::string> & deviceNames );

  // MPDE related registrations:
  std::vector<double> getFastSourcePeriod(Parallel::Machine comm, const std::vector<std::string> &sourceNames);
  std::vector<double> registerFastSources(Parallel::Machine comm, const std::vector<std::string> &sourceNames);
  void deRegisterFastSources(const std::vector<std::string> &sourceNames);
  void deactivateSlowSources();
  void activateSlowSources();

  void setMPDEFlag( bool flagVal );
  void setBlockAnalysisFlag( bool flagVal );
  void setFastTime( double timeVal );

  void setSPAnalysisFlag(bool flagVal);

  // Initialization function, to be called after all registrations are
  // finished, and the linear system class is completely set up.
  bool initializeAll(Linear::System &linear_system);

  // Device accessor functions:

  bool addDeviceModel(const ModelBlock & MB);

  std::pair<ModelTypeId, ModelTypeId> getModelType(const InstanceBlock &instance_block);

  bool verifyDeviceInstance(const InstanceBlock & IB);

  DeviceInstance * addDeviceInstance(const InstanceBlock & IB);

  bool deleteDeviceInstance (const std::string & name);

  int getHomotopyBlockSize() const;

  bool outputPlotFiles(bool force_final_output);
  bool finishOutput ();

  void  dotOpOutput ();

  // Load functions:
  bool setInitialGuess (Linear::Vector * solVectorPtr);
  bool loadErrorWeightMask(Linear::Vector * deviceMaskPtr);

  void debugOutput1();
  void debugOutput2();

  void getAnalyticSensitivities(
      const std::string & name,
      std::vector<double> & dfdpVec,
      std::vector<double> & dqdpVec,
      std::vector<double> & dbdpVec,
      std::vector<int> & FindicesVec,
      std::vector<int> & QindicesVec,
      std::vector<int> & BindicesVec);

  void getNumericalSensitivities(
      const std::string & name,
      std::vector<double> & dfdpVec,
      std::vector<double> & dqdpVec,
      std::vector<double> & dbdpVec,
      std::vector<int> & FindicesVec,
      std::vector<int> & QindicesVec,
      std::vector<int> & BindicesVec);

  bool analyticSensitivitiesAvailable(const std::string & name);
  bool numericalSensitivitiesAvailable(const std::string & name);
//

  void getAnalyticalBSensVectorsforAC (const std::string & name,
          std::vector< std::complex<double> > &     dbdp,
          std::vector<int> &        BindicesVec);

  void getNumericalBSensVectorsforAC (const std::string & name,
          std::vector< std::complex<double> > &     dbdp,
          std::vector<int> &        BindicesVec);

  void getAnalyticMatrixSensitivities(
      const std::string & name,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs);

  void getNumericalMatrixSensitivities(
      const std::string & name,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs);

  bool analyticMatrixSensitivitiesAvailable(const std::string & name);
  bool numericalMatrixSensitivitiesAvailable(const std::string & name);

  bool analyticBVecSensAvailable(const std::string & name);
  bool numericalBVecSensAvailable(const std::string & name);
//
  bool setParam(const std::string & name, double val, bool overrideOriginal = false);

  void setSeparateLoadFlag (bool flag) { devOptions_.separateLoad = flag; }
  bool getSeparateLoadFlag ()          { return devOptions_.separateLoad; }

  bool updateTemperature(double val);

  bool updateSources();

  bool resetRHSLoadFlags (int index);

  bool isLinearSystem()
  {
    return isLinearSystem_;
  }

  bool isPDESystem()
  {
    return solState_.isPDESystem_;
  }

  const DeviceOptions &getDeviceOptions() const
  {
    return devOptions_;
  }

  const SolverState &getSolverState() const
  {
    return solState_;
  }

  const ExternData &getExternData() const
  {
    return externData_;
  }

  const ArtificialParameterMap &getArtificialParameterMap() const
  {
    return artificialParameterMap_;
  }

  const PassthroughParameterSet &getPassthroughParameterSet() const
  {
    return passthroughParameterSet_;
  }

  bool getVoltageLimiterFlag()
  {
    return devOptions_.voltageLimiterFlag;
  }

  double getGmin()
  {
    return devOptions_.gmin;
  }

  // setup initial conditions on devices
  bool setICs (Linear::Vector * tmpSolVectorPtr,
               Linear::Vector * tmpCurrSolVectorPtr,
               Linear::Vector * tmpLastSolVectorPtr,
               Linear::Vector * tmpStaVectorPtr,
               Linear::Vector * tmpCurrStaVectorPtr,
               Linear::Vector * tmpLasStaVectorPtr,
               Linear::Vector * tmpStaDerivVectorPtr,
               Linear::Vector * tmpStoVectorPtr,
               Linear::Vector * tmpCurrStoVectorPtr,
               Linear::Vector * tmpQVectorPtr,
               Linear::Vector * tmpFVectorPtr,
               Linear::Vector * tmpBVectorPtr,
               Linear::Vector * tmpdFdxdVpVectorPtr,
               Linear::Vector * tmpdQdxdVpVectorPtr);

  // time integration stuff:
  void   resetBreakPoints();

  bool   getBreakPoints(
      std::vector<Util::BreakPoint> & breakPointTimes,
      std::vector<Util::BreakPoint> & pauseBreakPointTimes
      );
  double getMaxTimeStepSize ();

  // two-level newton and pde-continuation
  int  enablePDEContinuation();
  bool disablePDEContinuation ();
  void getNumInterfaceNodes (std::vector<int> & numInterfaceNodes);
  bool loadCouplingRHS (int iPDEDevice, int iElectrode, Linear::Vector * dfdvPtr);
  bool calcCouplingTerms (int iSubProblem, int iElectrode, const Linear::Vector * dxdvPtr);

  bool calcPDESubProblemInfo ();

  // load functions:
  bool loadDAEMatrices  (Linear::Vector * tmpSolVectorPtr,
                         Linear::Vector * tmpStaVectorPtr,
                         Linear::Vector * tmpStaDerivVectorPtr,
                         Linear::Vector * tmpStoVectorPtr,
                         Linear::Matrix * tmpdQdxMatrixPtr,
                         Linear::Matrix * tmpdFdxMatrixPtr,
                         int loadType = Xyce::Device::ALL);

  bool loadDAEVectors   (Linear::Vector * tmpNextSolVectorPtr,
                         Linear::Vector * tmpCurrSolVectorPtr,
                         Linear::Vector * tmpLastSolVectorPtr,
                         Linear::Vector * tmpNextStaVectorPtr,
                         Linear::Vector * tmpCurrStaVectorPtr,
                         Linear::Vector * tmpLastStaVectorPtr,
                         Linear::Vector * tmpStaDerivVectorPtr,
                         Linear::Vector * tmpNextStoVectorPtr,
                         Linear::Vector * tmpCurrStoVectorPtr,
                         Linear::Vector * tmpLeadFCompVectorPtr,
                         Linear::Vector * tmpLeadQCompVectorPtr,
                         Linear::Vector * tmpJunctionVCompVectorPtr,
                         Linear::Vector * tmpQVectorPtr,
                         Linear::Vector * tmpFVectorPtr,
                         Linear::Vector * tmpBVectorPtr,
                         Linear::Vector * tmpdFdxdVpVectorPtr,
                         Linear::Vector * tmpdQdxdVpVectorPtr,
                         int loadType = Xyce::Device::ALL);

  // For frequency-based analysis techniques. Specifically, allows all devices
  // to set internal variables needed by the load functions in advance.
  bool updateFDIntermediateVars(double frequency,
                                std::complex<double>* freqSolVec);

  // For frequency-based analysis techniques. Specifically, for devices that
  // have an analytic frequency reponse.  Method called to load the DAE
  // contributions to the Jacobian matrix.
  bool loadFreqDAEMatrices(double frequency,
                           std::complex<double>* freqSolVec,
                           std::vector<Util::FreqMatEntry>& dFdxMatrix);

  // For frequency-based analysis techniques. Specifically, for devices that
  // have an analytic frequency reponse.  Method called to load the DAE vectors
  // that contribute to the residual (RHS) vector.
  bool loadFreqDAEVectors(double frequency, std::complex<double>* freqSolVec,
                          std::vector<Util::FreqVecEntry>& freqFVector,
                          std::vector<Util::FreqVecEntry>& freqBVector);

  bool updateState(
    Linear::Vector * nextSolVectorPtr,
    Linear::Vector * currSolVectorPtr,
    Linear::Vector * lastSolVectorPtr,
    Linear::Vector * nextStaVectorPtr,
    Linear::Vector * currStaVectorPtr,
    Linear::Vector * lastStaVectorPtr,
    Linear::Vector * nextStoVectorPtr,
    Linear::Vector * currStoVectorPtr,
    int loadType = Xyce::Device::ALL);

  bool loadBVectorsforAC (Linear::Vector * bVecRealPtr,
                          Linear::Vector * bVecImagPtr);

  bool loadBVectorsforSources ();

  bool getBMatrixEntries(std::vector<int>& bMatEntriesVec, std::vector<int>& portsVec, std::vector<double> * Z0sVec  = NULL); 

  int getNumNoiseDevices ();
  int getNumNoiseSources ();
  void setupNoiseSources(std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec);
  void getNoiseSources(std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec);

  // voltlim doesn't work for all analysis types, so neccessary to
  // be able to turn it off/on.
  void unsetVoltageLimiterFlag ();
  void setVoltageLimiterFlag ( bool flagVal );

  void addGlobalPar(const Util::Param &param);
  const double *findGlobalPar(const std::string &name) const;

  // functions related to options registration
  bool setDeviceOptions(const Util::OptionBlock & option_block);
  bool registerSensParams(const Util::OptionBlock & option_block);
  bool setOPAnalysisParams(const Util::OptionBlock & option_block);
  bool setHBAnalysisParams(const Util::OptionBlock & option_block);
  bool setACAnalysisParams(const Util::OptionBlock & option_block);
  bool setNOISEAnalysisParams(const Util::OptionBlock & option_block);
  bool setSensitivityOptions (const Util::OptionBlock & option_block);

  bool getHBSpecified() const
  {
    return HBSpecified_;
  }

  bool getACSpecified() const
  {
    return ACSpecified_;
  }

  bool getIStarRequested() const
  {
    return iStarRequested_;
  }

  void setIStarRequested(bool flagVal);

  // convergence:  allow devices to signal back to the solvers that
  // they've played some game that invalidates normal convergence tests,
  // and so the solution should be considered unconverged no matter how
  // small the various norms are.
  bool allDevicesConverged(Parallel::Machine comm) const;

  // Functions needed for power node (2-level) algorithm):

  // for the parallel case, we need to give all the processors a copy
  // of the device so all the parallel synchronized calls such are
  // called by all processors together.
  bool setupExternalDevices(N_PDS_Comm &parallel_comm);

  const DeviceCountMap &getDeviceCountMap()
  {
    return localDeviceCountMap_;
  }

  void addDeviceToCount(const std::string & device_name, int num_devs = 1)
  {
    localDeviceCountMap_[device_name] += num_devs;
  }

  void addDevicesToCount(const DeviceCountMap &device_map);

  DeviceEntity *getDeviceEntity(const std::string &full_param_name) const;

  void acceptStep();

  bool getInitialQnorm (std::vector<TimeIntg::TwoLevelError> & tleVec );
  bool getInnerLoopErrorSums (std::vector<TimeIntg::TwoLevelError> & tleVec) const;
 
  const std::set<std::string>& getDevicesNeedingLeadCurrentLoads() const
  {
    return devicesNeedingLeadCurrentLoads_;
  }

  bool updateStateArrays();
  bool startTimeStep(
    bool                          beginIntegrationFlag,
    double                        nextTimeStep,
    double                        nextTime,
    int                           currentOrder);
  void setExternalSolverState(bool external_initJctFlag);

  int restartDataSize(bool pack) const;

  // Output restart data.
  bool dumpRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack ) const;

  // Load restart data.
  bool restoreRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack );

  // needed for parallel only:
  void setGlobalFlags();

  Device *getDevice(EntityTypeId model_type_id)
  {
    EntityTypeIdDeviceMap::iterator it = deviceMap_.find(model_type_id);
    return it == deviceMap_.end() ? 0 : (*it).second;
  }

  EntityTypeId getModelGroup(const std::string &model_or_device_type_name);

  const InstanceVector &getInstancePtrVec() const
  {
    return instancePtrVec_;
  }

  void addArtificialParameter(const std::string &name, ArtificialParameters::ArtificialParameter *artificial_parameter) 
  {
    artificialParameterMap_[name] = artificial_parameter;
    passthroughParameterSet_.insert(name);
  }

  const InstanceVector &getDevices(ModelTypeId model_type_id) const;

  // voltage limiter toggle functions
  bool getVoltageLimiterStatus();
  void setVoltageLimiterStatus(bool voltageLimterStatus);

  Util::Op::Operator *getOp(Parallel::Machine comm, const std::string &name) const;

  bool parameterExists(Parallel::Machine comm, const std::string & name) const;

// private:
//   bool getParamNoReduce(const std::string &name, double &value) const;

private:
  Device &getDeviceByModelType(const EntityTypeId model_type);

  bool setupRawVectorPointers_ ();
  bool setupRawMatrixPointers_ ();

  bool updateIntermediateVars_();
  bool updatePrimaryState_();
  bool updateSecondaryState_();

  bool updateDependentParameters_();

  // Do the actual solve/calculation for the external devices
  void updateExternalDevices_();

  // add external devices for processors that don't own it
  ExternDevice::Instance * addExtDeviceInstance_(const InstanceBlock & IB);

private:
  const IO::CmdParse &          commandLine_;           ///< Command line
  Topo::Topology &              topology_;              ///< Topology
  DeviceOptions                 devOptions_;            ///< user-defined options:
  EntityTypeIdDeviceMap         deviceMap_;

  bool                          sensFlag_;              ///< .SENS present in netlist
  bool                          isLinearSystem_;        ///< True if all devices in netlist have isLinearDevice() true
  bool                          firstDependent_;        ///< True until updateDependentParameters_ is called.
  bool                          parameterChanged_;      ///< Only used locally in updateDependentParameters_, don't know if stateful of just a member for fun
  bool                          breakPointInstancesInitialized;
  double                        timeParamsProcessed_;   ///< Time updateDependentParameters was called

  ExternData                    externData_;
  MatrixLoadData                matrixLoadData_;        ///< temporary jacobian load structures:
  SolverState                   solState_;              ///< real time solver data:
  Globals &                     globals_;               ///< global variables
  bool                          externalInitJctFlag_;
  bool                          externalStateFlag_;

  Analysis::AnalysisManager *   analysisManager_;       ///< To search for non-device parameters.  This needs to be removed
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> expressionGroup_; ///< required for setting up expressions

  // IO::Measure::Manager *        measureManager_;        ///< To search for non-device parameters.  This needs to be removed
  ArtificialParameterMap        artificialParameterMap_; ///< Specially named parameters.  This needs to be removed
  PassthroughParameterSet       passthroughParameterSet_;  ///< Parameters to pass through to external devices when set

  Parallel::Machine             comm_;                  ///< Communicator (should be passed in when needed, not a member)
  Nonlinear::Manager *          nlsMgrPtr_;             ///< To get Nonlinear solver information.  This needs to be removed

  ModelTypeNameModelTypeIdMap   modelTypeMap_;          ///< Model type name to model
  ModelTypeNameModelTypeIdMap   modelGroupMap_;         ///< Model type name to model group
  ModelTypeInstanceVectorMap    modelGroupInstanceVector_;
  ModelTypeInstanceVectorMap    modelTypeInstanceVector_;

  ModelVector                   modelVector_;
  ModelTypeModelVectorMap       modelTypeModelVector_;

  DeviceCountMap                localDeviceCountMap_;

  DeviceVector                  devicePtrVec_;
  DeviceVector                  pdeDevicePtrVec_;

  InstanceVector                instancePtrVec_;
  InstanceVector                devicesWithMaxTimeStepFuncsPtrVec_;

  //InstanceVector                bpInstancePtrVec_;    ///< instances with breakpoints functions
  InstanceVector                pauseBpInstancePtrVec_;    ///< instances with pause breakpoints functions
  InstanceVector                simpleBpInstancePtrVec_;    ///< instances with simple breakpoints functions
  InstanceVector                pdeInstancePtrVec_;
  InstanceVector                nonPdeInstancePtrVec_;
  InstanceVector                plotFileInstancePtrVec_;

  IndependentSourceVector       independentSourceVector_;
  IndependentSourceMap          independentSourceMap_;

  std::vector<DeviceInstance *> pauseBPDeviceVector_; /// < devices with pause breakpoints
  std::vector<DeviceInstance *> nonpauseBPDeviceVector_; /// < devices w/o pause BP(no BP, or only simple BP)

  // this is used to store the contents of the independentSourceVector_
  // during an mpde initialization where we'll remove slow sources from
  // the that vector so that they don't get updated
  IndependentSourceVector       indepSourceInstanceBackupPtrVec_;

  InstanceVector                testJacDevicePtrVec_;   ///< Devices under jacobian test
  EntityVector                  dependentPtrVec_;

  std::set<std::string>         devicesNeedingLeadCurrentLoads_;

  Util::Op::BuilderManager &    opBuilderManager_;
  mutable OpMap                 opMap_;

  mutable DeviceEntityMap       parameterDeviceCache_;                  ///< Full parameter name to device entity cache

  std::vector<int>              numInterfaceNodes_;
  bool                          calledBeforeCSPI;

  // sensitivities:
  //DeviceSensitivities *         devSensPtr_;
  DeviceSensitivities         devSens_;

  // .OP output flags.
  bool                          dotOpOutputRequested_;
  bool                          dotOpOutputFlag_;

  // analysis options  (for now bools that say what type of "." line was present in netlist)
  bool                          ACSpecified_;
  bool                          HBSpecified_;

  // used to enable lead-current calcuations for all devices.  This is set, during
  // netlist importation, when I(*), P(*) or W(*) appears on any .PRINT line.
  bool                          iStarRequested_;
};

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::unsetVoltageLimiterFlag ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/04
//-----------------------------------------------------------------------------
inline void DeviceMgr::unsetVoltageLimiterFlag ()
{
  devOptions_.voltageLimiterFlag = false;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setVoltageLimiterFlag ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 02/18/14
//-----------------------------------------------------------------------------
inline void DeviceMgr::setVoltageLimiterFlag ( bool flagVal )
{
  devOptions_.voltageLimiterFlag = flagVal;
}

bool registerPkgOptionsMgr(DeviceMgr &device_manager, IO::PkgOptionsMgr &options_manager);
void registerOpBuilders(Util::Op::BuilderManager &builder_manager, Parallel::Machine comm, DeviceMgr &device_manager);

double
getParamAndReduce(
  Parallel::Machine     comm,
  const DeviceMgr &     device_manager,
  const std::string &   name);

void addGlobalParameter(SolverState &solver_state, double temp, Globals &global, const Util::Param &param);
const double *findGlobalParameter(const GlobalParameterMap &global_map, const std::string &name);

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_DeviceMgr_h

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

//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_SweepParam.h>
#include <N_DEV_Algorithm.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_ExternDevice.h>
#include <N_DEV_MutIndLin.h>
#include <N_DEV_InstanceName.h>
#include <N_DEV_Op.h>
#include <N_DEV_Print.h>
#include <N_DEV_RegisterDevices.h>
#include <N_DEV_Source.h>
#include <N_ERH_ErrorMgr.h>
#include <N_ERH_Message.h>
#include <N_IO_CmdParse.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_LAS_Builder.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_NLS_Manager.h>
#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>
#include <N_PDS_Manager.h>
#include <N_TIA_WorkingIntegrationMethod.h>
#include <N_TOP_Topology.h>
#include <N_UTL_Algorithm.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Op.h>
#include <N_UTL_OpBuilder.h>
#include <N_UTL_Expression.h>

#include <expressionGroup.h>

// These are slowly being removed as DeviceMgr should not depend on Devices.  Devices should plug in to DeviceMgr.
#include <N_DEV_Bsrc.h>
#include <N_DEV_ISRC.h>
#include <N_DEV_MOSFET_B3.h>
#include <N_DEV_MOSFET_B3SOI.h>
#include <N_DEV_MOSFET_B4.h>
#include <N_DEV_Resistor.h>
#include <N_DEV_Resistor3.h>
#include <N_DEV_Vsrc.h>
#include <N_DEV_LTRA.h>
#include <N_DEV_TRA.h>
#include <N_DEV_ADC.h>

namespace Xyce {
namespace Device {

namespace {

// bool
// getParameter(
//   const ArtificialParameterMap &        artificial_parameter_map,
//   const GlobalParameterMap &            global_parameter_map,
//   const DeviceMgr &                     device_manager,
//   const IO::Measure::Manager &          measure_manager,
//   const std::string &                   name,
//   double &                              value);

//-----------------------------------------------------------------------------
bool setParameter(
  Parallel::Machine             comm,
  ArtificialParameterMap &      artificial_parameter_map,
  PassthroughParameterSet &     passthrough_parameter_map,
  Globals &                     globals,
  DeviceMgr &                   device_manager,
  EntityVector &                dependent_entity_vector,
  const InstanceVector &        extern_device_vector,
  const std::string &           name,
  double                        value,
  bool                          override_original);

//-----------------------------------------------------------------------------
//                 AGAUSS, GAUSS, AUNIF, UNIF, RAND and LIMIT
bool setParameterRandomExpressionTerms(
  Parallel::Machine             comm,
  ArtificialParameterMap &      artificial_parameter_map,
  PassthroughParameterSet &     passthrough_parameter_map,
  Globals &                     globals,
  DeviceMgr &                   device_manager,
  EntityVector &                dependent_entity_vector,
  const InstanceVector &        extern_device_vector,
  const std::string &           name,
  const std::string &           opName,
  int opIndex,
  //enum Util::astRandTypes astType,
  int astType,
  double                        value,
  bool                          override_original);

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::DeviceMgr
// Purpose       : constructor
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
DeviceMgr::DeviceMgr(
  Parallel::Machine             comm,
  Topo::Topology &              topology,
  Util::Op::BuilderManager &    op_builder_manager,
  const IO::CmdParse &          command_line)
  : commandLine_(command_line),
    topology_(topology),
    devOptions_(),
    deviceMap_(),
    sensFlag_(false),
    isLinearSystem_(true),
    firstDependent_(true),
    breakPointInstancesInitialized(false),
    timeParamsProcessed_(0.0),
    freqParamsProcessed_(0.0),
    externData_(),
    matrixLoadData_(),
    solState_(),
    globals_(solState_.globals_),
    analysisManager_(0),
    // measureManager_(0),
    artificialParameterMap_(),
    passthroughParameterSet_(),
    comm_(comm),
    nlsMgrPtr_(0),
    modelTypeMap_(),
    modelGroupMap_(),
    modelGroupInstanceVector_(),
    modelTypeInstanceVector_(),
    modelVector_(),
    modelTypeModelVector_(),
    localDeviceCountMap_(),
    devicePtrVec_(),
    pdeDevicePtrVec_(),
    instancePtrVec_(),
    //bpInstancePtrVec_(),
    pauseBpInstancePtrVec_(),
    simpleBpInstancePtrVec_(),
    pdeInstancePtrVec_(),
    nonPdeInstancePtrVec_(),
    plotFileInstancePtrVec_(),
    independentSourceVector_(),
    independentSourceMap_(),
    pauseBPDeviceVector_(),
    nonpauseBPDeviceVector_(),
    indepSourceInstanceBackupPtrVec_(),
    testJacDevicePtrVec_(),
    dependentPtrVec_(),
    devicesNeedingLeadCurrentLoads_(),
    opBuilderManager_(op_builder_manager),
    parameterDeviceCache_(),
    numInterfaceNodes_(),
    calledBeforeCSPI (false),
    devSens_(*this, devOptions_),
    dotOpOutputRequested_(false),
    dotOpOutputFlag_(false),
    ACSpecified_(false),
    HBSpecified_(false),
    iStarRequested_(false),
    expressionBasedSamplingEnabled_(false)
{
  addArtificialParameter("MOSFET:GAINSCALE", new ArtificialParameters::MOSFETGainScaleParam());
  addArtificialParameter("MOSFET:GAIN", new ArtificialParameters::MOSFETGainScaleParam());
  addArtificialParameter("MOSFET:NLTERMSCALE", new ArtificialParameters::MOSFETNLTermScaleParam());
  addArtificialParameter("MOSFET:NLTERM", new ArtificialParameters::MOSFETNLTermScaleParam());
  addArtificialParameter("MOSFET:L", new ArtificialParameters::MOSFETLParam());
  addArtificialParameter("MOSFET:W", new ArtificialParameters::MOSFETWParam());
  addArtificialParameter("MOSFET:SIZESCALE", new ArtificialParameters::MOSFETSizeScaleParam());
  addArtificialParameter("MOSFET:TOX", new ArtificialParameters::MOSFETTOXParam());
  addArtificialParameter("BJT:BF", new ArtificialParameters::BJTBFParam());
  addArtificialParameter("BJT:NF", new ArtificialParameters::BJTNFParam());
  addArtificialParameter("BJT:NR", new ArtificialParameters::BJTNRParam());
  addArtificialParameter("DIODE:N", new ArtificialParameters::DiodeNParam());
  addArtificialParameter("VSRCSCALE", new ArtificialParameters::VsrcScaleParam());
  addArtificialParameter("PDEALPHA", new ArtificialParameters::PDEAlphaParam());
  addArtificialParameter("PDEBETA", new ArtificialParameters::PDEBetaParam());
  addArtificialParameter("PDECHARGEALPHA", new ArtificialParameters::PDEChargeAlphaParam());
  addArtificialParameter("GSTEPPING", new ArtificialParameters::GSteppingParam());
  addArtificialParameter("GMIN", new ArtificialParameters::GMinParam());
  addArtificialParameter("VT", new ArtificialParameters::VtParam());
  addArtificialParameter("TEMP", new ArtificialParameters::TempParam());

  passthroughParameterSet_.insert("MOSFET_ALL:GAINSCALE");
  passthroughParameterSet_.insert("MOSFET_ALL:NLTERMSCALE");
  passthroughParameterSet_.insert("MOSFET1:GAINSCALE");
  passthroughParameterSet_.insert("MOSFET1:NLTERMSCALE");

  if (DEBUG_DEVICE) 
  {
    IO::setSensitivityDebugLevel(command_line, 0);
    IO::setDeviceDebugLevel(command_line, 1);
  }

  if (command_line.argExists("-jacobian_test"))
  {
    devOptions_.testJacobianFlag = true;
  }

  solState_.debugTimeFlag = true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerAnalysisManager
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
bool DeviceMgr::registerAnalysisManager(Analysis::AnalysisManager * analysis_manager)
{
  analysisManager_ = analysis_manager;

  return analysisManager_ != 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerExpressionGroup
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/19/2020
//-----------------------------------------------------------------------------
bool DeviceMgr::registerExpressionGroup(Teuchos::RCP<Xyce::Util::baseExpressionGroup> & group)
{
  expressionGroup_ = group;
  solState_.registerExpressionGroup(expressionGroup_);
  return (!(Teuchos::is_null(expressionGroup_)));
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::~DeviceMgr
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
DeviceMgr::~DeviceMgr()
{
  delete externData_.tmpdIdXPtr;
  delete externData_.tmpdQdXPtr;

  for (EntityTypeIdDeviceMap::iterator it = deviceMap_.begin(), end = deviceMap_.end(); it != end; ++it)
  {
    delete (*it).second;
  }

  for (OpMap::iterator it = opMap_.begin(), end = opMap_.end(); it != end; ++it)
  {
    delete (*it).second;
  }

  for (ArtificialParameterMap::iterator it = artificialParameterMap_.begin(), end = artificialParameterMap_.end(); it != end; ++it)
  {
    delete (*it).second;
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerSensParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool DeviceMgr::registerSensParams (const Util::OptionBlock & OB)
{
  sensFlag_ = true;
  return devSens_.registerSensParams (OB);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setSensitivityOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool DeviceMgr::setSensitivityOptions (const Util::OptionBlock & OB)
{
  return devSens_.setSensitivityOptions (OB);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setDeviceOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool DeviceMgr::setDeviceOptions (const Util::OptionBlock & OB)
{
  return devOptions_.setOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerLeadCurrentRequests
// Purpose       : this function is called from the output manager (through the
//                 device interface) to inform the device package of the devices
//                 for which lead currents have been requested.  The device
//                 manager will take care of doing isolated F and Q loads for
//                 these devices so the lead currents can be calculated
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical Systems Modeling
// Creation Date : 03/20/13
//-----------------------------------------------------------------------------
bool DeviceMgr::setLeadCurrentRequests(const std::set<std::string> & deviceNames)
{
  // this is called prior to fully constructing the devices.  This function is called 
  // multiple times, once for each "manager" (output, measure, fourier, ...) that can
  // use lead currents or power.
  for (std::set<std::string>::iterator it = deviceNames.begin(), end = deviceNames.end(); it != end; ++it)
  {
    devicesNeedingLeadCurrentLoads_.insert(*it);
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    if (! devicesNeedingLeadCurrentLoads_.empty())
    {
      dout() << "DeviceMgr::registerLeadCurrentRequests Devices for which lead currents were requested: ";
      for (std::set<std::string>::iterator it = devicesNeedingLeadCurrentLoads_.begin(), 
          end = devicesNeedingLeadCurrentLoads_.end(); it != end; ++it)
      {
        dout() << (*it) << "  ";
      }
      dout() << std::endl;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setOPAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/14/11
//-----------------------------------------------------------------------------
bool DeviceMgr::setOPAnalysisParams (const Util::OptionBlock & OB)
{
  dotOpOutputRequested_ = true;
  return true;
}
//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setHBAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/14/11
//-----------------------------------------------------------------------------
bool DeviceMgr::setHBAnalysisParams (const Util::OptionBlock & OB)
{
  HBSpecified_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setACAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/14/11
//-----------------------------------------------------------------------------
bool DeviceMgr::setACAnalysisParams (const Util::OptionBlock & OB)
{
  ACSpecified_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setNOISEAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/14/2014
//-----------------------------------------------------------------------------
bool DeviceMgr::setNOISEAnalysisParams (const Util::OptionBlock & OB)
{
  ACSpecified_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getDevices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
const InstanceVector & DeviceMgr::getDevices(ModelTypeId model_type_id) const
{
  static InstanceVector s_emptyVector;

  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(model_type_id);
  if (model_type_it != modelTypeInstanceVector_.end())
    return (*model_type_it).second;
  else
    return s_emptyVector;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getVoltageLimiterStatus
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 6/22/2015
//---------------------------------------------------------------------------
bool DeviceMgr::getVoltageLimiterStatus()
{
  return devOptions_.voltageLimiterFlag;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setVoltageLimiterStatus
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 6/22/2015 
//---------------------------------------------------------------------------
void DeviceMgr::setVoltageLimiterStatus(bool voltageLimterStatus)
{
  devOptions_.voltageLimiterFlag = voltageLimterStatus;
  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setIStarRequested
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 9/23/2018
//---------------------------------------------------------------------------
void DeviceMgr::setIStarRequested(bool iStarRequested)
{
  iStarRequested_ = iStarRequested;
  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getFastSourcePeriod
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/27/04
//-----------------------------------------------------------------------------
std::vector<double> DeviceMgr::getFastSourcePeriod(
  Parallel::Machine                     comm,
  const std::vector<std::string> &      sourceNames)
{
  int numFastSrcs = sourceNames.size();

  // Setup return of source periods
  std::vector<double> srcPeriods(numFastSrcs);

  // Now loop over them, and mark them.
  for (int i = 0; i < numFastSrcs; ++i)
  {
    IndependentSourceMap::const_iterator it = independentSourceMap_.find(sourceNames[i]);
    if (it != independentSourceMap_.end())
    {
      const SourceInstance &source_instance = *(*it).second;
      srcPeriods[i] = source_instance.period();
    }
    else
    {
#ifndef Xyce_PARALLEL_MPI
      Report::UserError message;
      message << "Unable to find source: " <<  sourceNames[i] << "\n"
              << "Potential names are: ";
      for (IndependentSourceMap::const_iterator it = independentSourceMap_.begin(); 
          it != independentSourceMap_.end(); ++it)
      {
        message << (*it).first << " ";
      }
#endif
    }
  }

  Parallel::AllReduce(comm, MPI_MAX, srcPeriods);

  return srcPeriods;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerFastSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/27/04
//-----------------------------------------------------------------------------
std::vector<double> DeviceMgr::registerFastSources(
  Parallel::Machine                     comm,
  const std::vector<std::string> &      sourceNames)
{
  int numFastSrcs = sourceNames.size();

  std::vector<double> srcPeriods;

  if (numFastSrcs > 0)
  {
    srcPeriods.resize(numFastSrcs, 0.0);

    // Now loop over them, and mark them.
    for (int i = 0; i < numFastSrcs; ++i)
    {
      IndependentSourceMap::iterator it = independentSourceMap_.find(sourceNames[i]);
      if (it != independentSourceMap_.end())
      {
        SourceInstance &source_instance = *(*it).second;
        source_instance.setFastSourceFlag(true);
        srcPeriods[i] = source_instance.period();
      }
      else
      {
#ifndef Xyce_PARALLEL_MPI
      Report::UserError message;
      message << "Unable to find source: " << sourceNames[i] << "\n"
              << "Potential names are: ";
      for (IndependentSourceMap::const_iterator it = independentSourceMap_.begin(); 
          it != independentSourceMap_.end(); ++it)
      {
        message << (*it).first << " ";
      }
#endif
      }
    }

  }
  else
  {
    // tscoffe/tmei 09/16/08
    // Special case:  Use all sources
    // Compute the total number of fast sources for all processors.
    // NOTE:  In parallel, this will not work correctly if more than one processor has fast sources.

    int myNumFastSrcs = independentSourceVector_.size();
    Parallel::AllReduce(comm, MPI_SUM, &myNumFastSrcs, &numFastSrcs, 1);

    if (myNumFastSrcs > independentSourceVector_.size())
      throw std::runtime_error("registerFastSources() does not handle parallel");

    srcPeriods.resize(numFastSrcs, -1.0);
    for (int i = 0; i < myNumFastSrcs; ++i)
    {
      independentSourceVector_[i]->setFastSourceFlag(true);
      srcPeriods[i] = independentSourceVector_[i]->period();
    }
  }

  Parallel::AllReduce(comm, MPI_MAX, srcPeriods);

  return srcPeriods;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::deRegisterFastSources
// Purpose       : reverses the effect of registerFastSources
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/12/2013
//-----------------------------------------------------------------------------
void DeviceMgr::deRegisterFastSources(const std::vector<std::string> &sourceNames)
{
  int numFastSrcs = sourceNames.size();

  if (numFastSrcs > 0)
  {
    // Now loop over them, and mark them.
    for (int i = 0; i < numFastSrcs; ++i)
    {
      IndependentSourceMap::iterator it = independentSourceMap_.find(sourceNames[i]);
      if (it != independentSourceMap_.end())
      {
        SourceInstance &source_instance = *(*it).second;
        source_instance.setFastSourceFlag(false);
      }
      else
      {
#ifndef Xyce_PARALLEL_MPI
        Report::DevelFatal message;
        message.in("DeviceMgr::deRegisterFastSources");
        message << "Unable to find source: " <<  sourceNames[i] << std::endl
                << "Potential names are: ";
        for (IndependentSourceMap::const_iterator it = independentSourceMap_.begin(); 
            it != independentSourceMap_.end(); ++it)
        {
          message << (*it).first << " ";
        }
#endif
      }
    }
  }
  else
  {
    // Special case:  Use all sources
    numFastSrcs = independentSourceVector_.size();
    for (int i = 0; i < numFastSrcs; ++i)
    {
      independentSourceVector_[i]->setFastSourceFlag(false);
    }
  }
  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::deactivateSlowSources
// Purpose       : traverse fast source list and remove any slow sources from
//                 the deviceArray
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 3/22/07
//-----------------------------------------------------------------------------
void DeviceMgr::deactivateSlowSources()
{
  // first back-up a copy of the deviceArray so we can edit out the slow sources
  indepSourceInstanceBackupPtrVec_ = independentSourceVector_;

  // erase the existing list of sources
  independentSourceVector_.clear();

  // now copy back only those that are fast sources
  for (IndependentSourceVector::iterator it = indepSourceInstanceBackupPtrVec_.begin(), 
      end = indepSourceInstanceBackupPtrVec_.end(); it != end; ++it)
  {
    if ((*it)->getFastSourceFlag())
    {
      independentSourceVector_.push_back(*it);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::activateSlowSources
// Purpose       : restore any slow sources to the device array.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 3/22/07
//-----------------------------------------------------------------------------
void DeviceMgr::activateSlowSources()
{
  // restore the independent source list from backup
  independentSourceVector_ = indepSourceInstanceBackupPtrVec_;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setSPAnalysisFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void DeviceMgr::setSPAnalysisFlag(bool flagVal)
{
  solState_.spAnalysisFlag_ = flagVal;

}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setMPDEFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 07/21/08
//-----------------------------------------------------------------------------
void DeviceMgr::setMPDEFlag(bool flagVal)
{
  solState_.mpdeOnFlag_  = flagVal;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setBlockAnalysisFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 07/21/08
//-----------------------------------------------------------------------------
void DeviceMgr::setBlockAnalysisFlag(bool flagVal)
{
  solState_.blockAnalysisFlag_ = flagVal;

  // tscoffe/tmei 07/30/08:  Note, if we call this with "false", then the
  // newExcessPhase and defaultNewExcessPhase flags will not go back to their
  // "user-set" values, they will be set to false.
  devOptions_.newExcessPhase = flagVal;
  devOptions_.defaultNewExcessPhase = flagVal;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setFastTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 07/21/08
//-----------------------------------------------------------------------------
void DeviceMgr::setFastTime(double timeVal)
{
  solState_.currFastTime_ = timeVal;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::initializeAll
// Purpose       : This function, via the LAS system class, sets up
//                 the pointers to the various linear algebra entities.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/12/00
//-----------------------------------------------------------------------------
bool DeviceMgr::initializeAll(Linear::System &linear_system)
{
  bool bsuccess = true;

  externData_.lasSysPtr = &linear_system;

  // nullify ptrs that are passed in at each step  (see the loadDAEVectors function args)
  externData_.nextSolVectorPtr = 0;
  externData_.currSolVectorPtr = 0;
  externData_.lastSolVectorPtr = 0;
  externData_.daeQVectorPtr    = 0;
  externData_.daeFVectorPtr    = 0;
  externData_.daeBVectorPtr    = 0;
  externData_.dFdxdVpVectorPtr = 0;
  externData_.dQdxdVpVectorPtr = 0;
  externData_.nextStaVectorPtr = 0;
  externData_.currStaVectorPtr = 0;
  externData_.lastStaVectorPtr = 0;
  externData_.nextStaDerivVectorPtr = 0;
  externData_.nextStoVectorPtr = 0;
  externData_.currStoVectorPtr = 0;
  externData_.nextLeadCurrFCompPtr = 0;
  externData_.nextLeadCurrQCompPtr = 0;
  externData_.nextJunctionVCompPtr = 0;

  if (DEBUG_DEVICE)
  {
    // create Jdxp vector pointer:
    externData_.JdxpVectorPtr = linear_system.getJDXPVector();
    bsuccess = bsuccess && (externData_.JdxpVectorPtr != 0);
  }

  // get flag solution pointer-pointer:
  externData_.flagSolVectorPtr = linear_system.getFlagSolVector();
  bsuccess = bsuccess && (externData_.flagSolVectorPtr != 0);

  // get device mask pointer.
  externData_.deviceErrorWeightMask_ = linear_system.getDeviceMaskVector ();
  bsuccess = bsuccess && (externData_.deviceErrorWeightMask_ != 0);

  externData_.tmpdIdXPtr = linear_system.builder().createVector();
  externData_.tmpdQdXPtr = linear_system.builder().createVector();

  externData_.initializeAllFlag = true;

  // For Homotopy on block gainscale
  solState_.initializeHomotopyBlockSize(devOptions_.numGainScaleBlocks);

#ifdef Xyce_SIZEOF
  int size = sizeof(*this);
  dout() << "Size of device package after initializeAll  = " << size << std::endl;
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::notify
// Purpose       : Receive "StepEvent" notifications from publishers of such.
// Special Notes : This method is part of the Xyce::Util::Listener interface,
//                 defined in the N_UTL_Listener.h file.  DeviceMgr derives
//                 from the Listener class.  The device manager is subscribed
//                 to StepEvent notifications from the analysis manager by
//                 a call to Util::subscribe<Analysis::StepEvent> up in the
//                 Xyce::Simulator::doAllocations() method.  StepEvent
//                 events are published by the
//                 Xyce::Analysis::Step::doLoopProcess and doInit methods of
//                 the step analysis class.
//
//                 This method replaces the "resetForStepAnalysis" method
//                 that used to be called more explicitly and less flexibly.
//                 The notify method was introduced in revision 1.626 by
//                 Dave Baur on 11 July 2014.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/12/2013
//-----------------------------------------------------------------------------
void DeviceMgr::notify(const Analysis::StepEvent &event)
{
  if (event.state_ == Analysis::StepEvent::STEP_STARTED)
  {
    delete externData_.tmpdIdXPtr;  externData_.tmpdIdXPtr = 0;
    delete externData_.tmpdQdXPtr;  externData_.tmpdQdXPtr = 0;

    if ( solState_.ltraDevices_ )
    {
      solState_.ltraTimeIndex_ = 0;
      solState_.ltraTimeHistorySize_ = 10;
      solState_.ltraTimePoints_.resize(solState_.ltraTimeHistorySize_);
    }
  }

  // Breakpoint management
  InstanceVector::iterator iterI;
  ModelVector::iterator iterM;
  for (iterM = modelVector_.begin() ; iterM != modelVector_.end() ; ++iterM)
  {
    if (!(*iterM)->getDependentParams().empty()) { (*iterM)->setupParamBreakpoints(); }
  }

  for (iterI = instancePtrVec_.begin() ; iterI != instancePtrVec_.end() ; ++iterI)
  {
    if (!(*iterI)->getDependentParams().empty()) { (*iterI)->setupParamBreakpoints(); }
  }

  std::vector<Util::Expression>::iterator globalExp_i = globals_.global_expressions.begin();
  std::vector<Util::Expression>::iterator globalExp_end = globals_.global_expressions.end();
  for (; globalExp_i != globalExp_end; ++globalExp_i)
  {
    globalExp_i->setupBreakPoints();
  }

  InstanceVector::iterator iter = instancePtrVec_.begin();
  InstanceVector::iterator end = instancePtrVec_.end();
  for (; iter!= end; iter++)
  {
    (*iter)->setupBreakPoints();
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getDeviceByModelType
// Purpose       : This function creates a single device based on the passed
//                 index.
// Special Notes :
// Scope         : protected
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Device &DeviceMgr::getDeviceByModelType(const ModelTypeId model_type_id)
{
  Device *device = 0;

  if (model_type_id.defined())
  {
    EntityTypeIdDeviceMap::const_iterator it = deviceMap_.find(model_type_id);
    if (it == deviceMap_.end())
    {
      FactoryBlock factory_block(*this, devOptions_, solState_, matrixLoadData_, externData_, commandLine_);

      const Configuration *configuration = Configuration::findConfiguration(model_type_id);
      device = configuration->createDevice(factory_block);
      deviceMap_[model_type_id] = device;
      devicePtrVec_.push_back(device);

      if (device->isPDEDevice())
      {
        pdeDevicePtrVec_.push_back(device);
      }
    }
    else
    {
      device = (*it).second;
    }
  }

  return *device;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getModelGroup
// Purpose       : This function returns the device type index for a given
//                 named string.  This assumes that the device names used
//                 in the simulation will obey the spice3f5 netlist language
//                 convention.
//
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
/**
 * Return the ModelGroup of the device associated with the model type name or device type name.
 *
 * @param model_type_name model type name or device type name
 *
 * @return
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355
 * @date   Mon Sep 23 07:53:04 2013
 */
EntityTypeId DeviceMgr::getModelGroup(const std::string &model_or_device_type_name)
{
  return Configuration::getModelGroup(model_or_device_type_name);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addDeviceModel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
bool DeviceMgr::addDeviceModel(const ModelBlock & model_block)
{
  ModelTypeId model_type;
  ModelTypeId model_group = Configuration::getModelGroup(model_block.getType());

  if (!model_block.getName().empty())
  {
    model_type = Configuration::getModelType(model_block.getType(), model_block.getLevel());

    if (!model_type.defined())
    {
      Report::UserError() << "There is no device " << model_block.getType() 
        << " of level " << model_block.getLevel() << " to define model " << model_block.getName();
    }
  }

  if (!model_type.defined())
    model_type = model_group;

  if (!model_type.defined())
    return false;

  FactoryBlock factory_block(*this, devOptions_, solState_, matrixLoadData_, externData_, commandLine_);
  Device &device = getDeviceByModelType(model_type);
  DeviceModel *device_model = device.addModel(model_block, factory_block);

  modelTypeMap_[model_block.getName()] = model_type;
  modelGroupMap_[model_block.getName()] = model_group;

  // add the various model vectors:
  if (device_model != 0)
  {
    modelVector_.push_back(device_model);
    modelTypeModelVector_[model_type].push_back(device_model);
  }

  return device_model != 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getModelType
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Baur
// Creation Date :
//-----------------------------------------------------------------------------
std::pair<ModelTypeId, ModelTypeId> DeviceMgr::getModelType(const InstanceBlock & instance_block)
{
  // If the ModelName string is not a null string, use it to get the
  // device type index.  Otherwise, use the name of the device itself to determine the device.
  ModelTypeId model_type;
  ModelTypeId model_group;
  if (instance_block.getModelName().empty())
  {
    model_type = getModelGroup(instance_block.getInstanceName().getDeviceType());
    model_group = model_type;
  }
  else
  {
    model_type = modelTypeMap_[instance_block.getModelName()];
    model_group = modelGroupMap_[instance_block.getModelName()];
  }

  if (!model_type.defined())
  {
    Report::UserError message;
    message << "Unable to determine type of device for instance name " << instance_block.getInstanceName();
    if (!instance_block.getModelName().empty())
    {
      message << " with model name " << instance_block.getModelName();
    }
    return std::pair<ModelTypeId, ModelTypeId>(ModelTypeId(), ModelTypeId());
  }
  else if (instance_block.bsourceFlag)
  {
    // This is an E, F, G, or H source that is to be treated as a B source,
    // set its type now to BSRC.
    model_type = Bsrc::Traits::modelType();
  }

  // Check if this is a simple resistor, but with a resistance of zero.  If so,
  // then change the type to RESISTOR3.  This variant of the resistor acts like
  // a voltage souce with zero voltage difference.
  //
  // ERK.  7/30/2020.  This block of code is pretty fragile, and broke my new 
  // expression library multiple times.  The reason is that it indirectly 
  // modifies parameters that are of type Util::EXPR, as long as they 
  // don't have any dependencies.    
  //
  // It isn't obvious from staring at the code below that this happens.  However, 
  // the call to currentParam->getImmutableValue<double> changes the type from 
  // Util::EXPR to Util::DBLE.  It performs this change even if the result of 
  // the if-statement using this call is "false".  In a few cases, I didn't 
  // want this change to happen and it messed things up.  Note, the Util::Param object 
  // is supposed to check this stuff as well, and if an expression is non-constant,
  // then it emits a fatal error when getImmutableValue is called.  I've updated
  // that code as well.
  if (devOptions_.checkForZeroResistance && (model_type == Resistor::Traits::modelType()))
  {
    const double zeroResistanceValue = devOptions_.zeroResistanceTol;
    // loop over the parameters
    std::vector<Param>::const_iterator currentParam = instance_block.params.begin();
    std::vector<Param>::const_iterator endParam = instance_block.params.end();
    while(currentParam != endParam)
    {
      if ((currentParam->uTag() == "R"))
      {
        if (currentParam->given())
        {
          std::vector<std::string> variables, specials;
          bool isRandomDependent=false;

          // check if this is a time-dependent, or variable-dependent expression.
          // If it is, then skip.
          const Param * devPar = &(*(currentParam));
          const Util::Param * tmpPar = (dynamic_cast<const Util::Param*> (devPar));
          // only check if this is an expression-type parameter.
          if (tmpPar->getType() == Util::EXPR)
          {
            Util::Expression tmpExp = tmpPar->getValue<Util::Expression>();
            tmpExp.getVariables(variables); 
            tmpExp.getSpecials(specials);      
            isRandomDependent = tmpExp.isRandomDependent();
          }

          if (specials.empty() && variables.empty() && !isRandomDependent)
          {
            if (fabs(currentParam->getImmutableValue<double>()) < devOptions_.zeroResistanceTol) // call here will change param type from EXPR to DBLE
            {
              // change device type to be the level 3 resistor
              model_type = Resistor3::Traits::modelType();
            }
          }
        }
        break;
      }
      currentParam++;
    }
  }

  return std::pair<ModelTypeId, ModelTypeId>(model_type, model_group);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::verifyDeviceInstance
// Purpose       : This function verifies a device instance prior to
//                 instantiating.
//
//                 Theoretically, we could do this in addDeviceInstance() and
//                 not make a device if it fails some verification criteria (a
//                 resistor with zero resistance is the primary case here).
//                 However, later unlinking of redundant devices that are
//                 connected to just one node is difficult after
//                 addDeviceInstance() is called because the device instance
//                 pointer can be placed in many other containers
//                 It could be done, but this is a simpler first step to having the
//                 device manager be in charge of device verification -- rather
//                 than have it in topology or IO.
//
// Special Notes : return true if this device is ok to instantiate, false otherwise
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/18/2010
//-----------------------------------------------------------------------------
bool DeviceMgr::verifyDeviceInstance(const InstanceBlock & instance_block)
{
  ModelTypeId model_type;
  ModelTypeId model_group;
  std::pair<ModelTypeId, ModelTypeId> x = getModelType(instance_block);
  model_type = x.first;
  model_group = x.second;

  if (!model_type.defined())
  {
    Report::UserError message;
    message << "Unable to determine model type of device for instance name " << instance_block.getInstanceName();
    if (!instance_block.getModelName().empty())
    {
      message << " with model name" << instance_block.getModelName();
    }

    return false;
  }

  // Return false if we found a supernode.
  if (model_type == Resistor3::Traits::modelType())
  {
    return false;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addDeviceInstance
// Purpose       : addDeviceInstance will create a new instance of the
//                 designated device type.  This version of the function
//                 accepts a parameter list as one of the arguments,
//                 so it is assumed that a parameter instance will
//                 also have to be allocated for it.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
DeviceInstance * DeviceMgr::addDeviceInstance(
  const InstanceBlock &         instance_block)
{
  if (sensFlag_ && ACSpecified_) 
  { 
    devOptions_.matrixSensitivityFlag = true;
  }

  // experiment
  //bool all_devices_converged = allDevicesConverged(comm_);
  //setupSolverInfo(solState_, *analysisManager_, all_devices_converged, devOptions_, nlsMgrPtr_->getNonLinInfo());

  // If the ModelName string is not a null string, use it to get the
  // device type index.  Otherwise, use the name of the device itself to determine the device.
  ModelTypeId model_type;
  ModelTypeId model_group;
  std::pair<ModelTypeId, ModelTypeId> x = getModelType(instance_block);
  model_type = x.first;
  model_group = x.second;

  if (!model_type.defined())
  {
    Report::UserError message;
    message << "Unable to determine type of device for instance name " << instance_block.getInstanceName();
    if (!instance_block.getModelName().empty())
    {
      message << " with model name " << instance_block.getModelName();
    }

    DeviceInstance *instance = 0;
    return instance;
  }

  // Add an instance of this type.
  Device &device = getDeviceByModelType(model_type);
  DeviceInstance *instance = device.addInstance(
      instance_block,
      FactoryBlock(*this, devOptions_, solState_, matrixLoadData_, externData_, commandLine_));

  std::string outputName = (instance->getName()).getEncodedName();
  std::string deviceName = (instance->getName()).getDeviceName();
  // Special handling for mutual inductors.  We need to know the names of the
  // component inductors and then check if those component inductor names are
  // in the set devicesNeedingLeadCurrentLoads_
  bool mutualInductorNeedsLeadCurrents = false;
  if ( deviceName[0] == 'K' )
  {
    std::vector<std::string> inductorNames = instance->getInductorNames();
    for(int i=0; i < inductorNames.size(); ++i)
    {
      if (devicesNeedingLeadCurrentLoads_.find(inductorNames[i]) != devicesNeedingLeadCurrentLoads_.end())
      {
        mutualInductorNeedsLeadCurrents = true;
        break;  // break after finding the first match
      }
    }
  }
  
  if ( mutualInductorNeedsLeadCurrents || devOptions_.calculateAllLeadCurrents || 
       (devicesNeedingLeadCurrentLoads_.find(outputName) != devicesNeedingLeadCurrentLoads_.end()) ||
       iStarRequested_ )
  {
    instance->enableLeadCurrentCalc();

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      dout() << "DeviceMgr::addDeviceInstance Enabling lead current load for device \""
        << instance->getName()
        << "\" ->  \""
        << outputName
        << "\"" << std::endl;
    }
  }
  else if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    dout() << "DeviceMgr::addDeviceInstance Cannot enable lead current load for device \""
           << instance->getName()
           << "\" ->  \""
           << outputName
           << "\""
           << std::endl;
  }

  localDeviceCountMap_[device.getDefaultModelName()]++;

  isLinearSystem_ = isLinearSystem_ && instance->isLinearDevice();

  solState_.isPDESystem_ = solState_.isPDESystem_ || device.isPDEDevice();

  // Set up the instance vectors.  These are the main containers used in the load procedures.
  instancePtrVec_.push_back(instance);

  // set up the list of pde device instances
  // and the list of non-pde devices instances.
  if (device.isPDEDevice())
  {
    pdeInstancePtrVec_.push_back(instance);
  }
  else
  {
    nonPdeInstancePtrVec_.push_back(instance);
  }

  modelGroupInstanceVector_[model_group].push_back(instance);
  modelTypeInstanceVector_[model_type].push_back(instance);

  // set up the independent source map.
  if (model_type == Vsrc::Traits::modelType() || model_type == ISRC::Traits::modelType())
  {
    independentSourceMap_[instance_block.getInstanceName().getEncodedName()] = dynamic_cast<SourceInstance *>(instance);
    independentSourceVector_.push_back(dynamic_cast<SourceInstance *>(instance));
  }

  // add to the list of devices that have PAUSE breakpoints, and those that do not
  if ( model_type == ADC::Traits::modelType() )
  {
    pauseBPDeviceVector_.push_back(instance);
  }
  else
  {
    nonpauseBPDeviceVector_.push_back(instance);
  }

  // if this is an LTRA device, let the solver state know.
  if (model_type == LTRA::Traits::modelType())
  {
    solState_.ltraDevices_ = true;
  }

  if (instance->maxTimeStepSupported())
  {
    devicesWithMaxTimeStepFuncsPtrVec_.push_back(instance);
  }

  if (instance->plotfileFlag ())
  {
    plotFileInstancePtrVec_.push_back(instance);
  }

  // Set up the vector of devices subject to the jacobian test.
  if (instance->getName() == devOptions_.testJacDeviceName)
  {
    testJacDevicePtrVec_.push_back(instance);
  }

  return instance;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::deleteDeviceInstance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/08/01
//-----------------------------------------------------------------------------
bool DeviceMgr::deleteDeviceInstance(const std::string & name)
{
  Report::DevelFatal().in("DeviceMgr::deleteDeviceInstance")
    << "Not ready with the new boilerplate-free device package";

  bool bsuccess = true;

  // bool tmpBool = true;
  // EntityTypeId type = getModelGroup(name);

  // DeviceMap::iterator it = deviceMap_.find(type);
  // if (it != deviceMap_.end())
  // {
  //   // bsuccess &= (*it).second->deleteInstance(name);
  //   DeviceEntity *entity = (*it).second->findEntity(name);
  //   DeviceInstance *instance = dynamic_cast<DeviceInstance *>(entity);
  //   if (instance)
  //     bsuccess = (*it).second->deleteInstance(instance);
  // }

  // note as this is written it ignores lots of other containers
  // this may be used for clean up at the end of a run, but
  // it is not sufficient to use during simulation setup.
  //
  // need to remove pointer to the instance "name" from other arrays
  // candidate lists are:
  //     InstanceVector instancePtrVec_;
  //     InstanceVector bpInstancePtrVec_; // instances with breakpoints functions
  //     InstanceVector pdeInstancePtrVec_;
  //     InstanceVector nonPdeInstancePtrVec_;
  //     InstanceVector mosfetInstancePtrVec_;
  //     InstanceVector vsrcInstancePtrVec_;
  //     InstanceVector bjtInstancePtrVec_;
  //     std::map<std::string, VsrcInstance*> vsrcInstancePtrMap_;
  //
  //     InstanceVector plotFileInstancePtrVec_;
  //
  //     std::map<std::string,SourceInstance*> independentSourceMap_;
  //     IndependentSourceVector independentSourceVector_;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::debugOutput1
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 1/29/07
//-----------------------------------------------------------------------------
void DeviceMgr::debugOutput1()
{
  // dump the fvector and the jdxp vector to files.
  if (isActive(Diag::DEVICE_DUMP_VECTORS) && solState_.debugTimeFlag)
  {
    // To match the numbering scheme of the NLS debug output files,
    // it is neccessary to add +1 to the newton iterartion number.
    int newton_iteration_count = solState_.newtonIter + 1;

    // jdxp-vector
    // char fn_jdxp[256]; for (int ich = 0; ich < 256; ++ich) fn_jdxp[ich] = 0;
    // sprintf(fn_jdxp, "Jdxp.%03d.txt", newton_iteration_count);
    {
      std::ostringstream oss;
      oss << "Jdxp." << std::setw(3) << std::setfill('0') << newton_iteration_count << ".txt";
      (externData_.JdxpVectorPtr)->writeToFile(oss.str().c_str());
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::debugOutput2
// Purpose       : new-dae version
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/19/08
//-----------------------------------------------------------------------------
void DeviceMgr::debugOutput2()
{
  // dump the fvector and the jdxp vector to files.
  if (isActive(Diag::DEVICE_DUMP_VECTORS) && solState_.debugTimeFlag)
  {
    // To match the numbering scheme of the NLS debug output files,
    // it is neccessary to add +1 to the newton iterartion number.
    int newton_iteration_count = solState_.newtonIter + 1;
    int outputStepNumber = 0;

    if (solState_.tranopFlag)
    {
      outputStepNumber = 0;
    }
    else if (solState_.initTranFlag_)
    {
      outputStepNumber = solState_.timeStepNumber_ + 1;
    }
    else
    {
      outputStepNumber = solState_.timeStepNumber_ + 1;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setInitialGuess
// Purpose       : This is a function call that sets the initial guess for
// devices that have initial guesses.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/02/03
//-----------------------------------------------------------------------------
bool DeviceMgr::setInitialGuess (Linear::Vector * solVectorPtr)
{
  bool bsuccess = true;

  if (solVectorPtr != 0)
  {
    externData_.nextSolVectorPtr = solVectorPtr;

    // if two-level, and just the inner problem, only load the PDE devices.
    for (InstanceVector::iterator it = instancePtrVec_.begin(), 
        end = instancePtrVec_.end (); it != end; ++it)
    {
      bool tmpBool = (*it)->setInitialGuess();
      bsuccess = bsuccess && tmpBool;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::analyticSensitivitiesAvailable
// Purpose       :
// Special Notes : Includes a reduce, so works in parallel.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/30/2014
//-----------------------------------------------------------------------------
bool DeviceMgr::analyticSensitivitiesAvailable(const std::string & name)
{
  // This assumes that the only available analytic sensitivities are for
  // non-artificial parameters.
  DeviceEntity * device_entity = getDeviceEntity(name);

  int available = 0;
  if (device_entity)
  {
    std::string paramName = Util::paramNameFromFullParamName(name);

    if (paramName == "")
    {
      available = device_entity->analyticSensitivityAvailableDefaultParam();
    }
    else
    {
      available = device_entity->analyticSensitivityAvailable(paramName);
    }
  }

  Parallel::AllReduce(comm_, MPI_LOR, &available, 1);

  return available != 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::numericalSensitivitiesAvailable
// Purpose       :
// Special Notes : Includes a reduce, so works in parallel.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool DeviceMgr::numericalSensitivitiesAvailable(const std::string & name)
{
  // This assumes that the only available analytic sensitivities are for
  // non-artificial parameters.
  DeviceEntity * device_entity = getDeviceEntity(name);

  // if we find a valid device entity then device-level numerical derivatives will 
  // work.  Don't need to check anything else
  int available = 0;
  if (device_entity)
  {
    available = 1;
  }
  Parallel::AllReduce(comm_, MPI_LOR, &available, 1);

  return available != 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getAnalyticSensitivities
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/30/2014
//-----------------------------------------------------------------------------
void DeviceMgr::getAnalyticSensitivities
  (const std::string & name,
   std::vector<double> & dfdpVec,
   std::vector<double> & dqdpVec,
   std::vector<double> & dbdpVec,
   std::vector<int> & FindicesVec,
   std::vector<int> & QindicesVec,
   std::vector<int> & BindicesVec)
{
  // If not artificial, then search for the appropriate natural param(s).
  DeviceEntity * device_entity = getDeviceEntity(name);

  bool found=false;
  if (device_entity)
  {
    std::string paramName = Util::paramNameFromFullParamName(name);
    if (paramName == "")
    {
      found = device_entity->getAnalyticSensitivityDefaultParam(
                                                    dfdpVec, dqdpVec, dbdpVec,
                                                    FindicesVec, QindicesVec, BindicesVec);
    }
    else
    {
      found = device_entity->getAnalyticSensitivity(paramName,
                                                    dfdpVec, dqdpVec, dbdpVec,
                                                    FindicesVec, QindicesVec, BindicesVec);
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getNumericalSensitivities
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 
//-----------------------------------------------------------------------------
void DeviceMgr::getNumericalSensitivities
  (const std::string & name,
   std::vector<double> & dfdpVec,
   std::vector<double> & dqdpVec,
   std::vector<double> & dbdpVec,
   std::vector<int> & FindicesVec,
   std::vector<int> & QindicesVec,
   std::vector<int> & BindicesVec)
{

  // what kind of parameter is this:  instance? model?  global_param?
  //
  // instance only needs the same function args as getAnalyticSensitivity
  //
  // model, however, needs a list of instances.
  //
  // global_param needs even more work.  Currently not supported here.
  //

  DeviceEntity * device_entity = getDeviceEntity(name);

  bool found=false;
  if (device_entity)
  {
    std::string paramName = Util::paramNameFromFullParamName(name);

    if (paramName == "")
    {
      found = device_entity->getNumericalSensitivityDefaultParam(
                                                    dfdpVec, dqdpVec, dbdpVec,
                                                    FindicesVec, QindicesVec, BindicesVec);
    }
    else
    {
      found = device_entity->getNumericalSensitivity(paramName,
                                                    dfdpVec, dqdpVec, dbdpVec,
                                                    FindicesVec, QindicesVec, BindicesVec);
    }
  }

  return;
}

//

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::analyticMatrixSensitivitiesAvailable
// Purpose       :
// Special Notes : Includes a reduce, so works in parallel.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool DeviceMgr::analyticMatrixSensitivitiesAvailable(const std::string & name)
{
  DeviceEntity * device_entity = getDeviceEntity(name);

  int available = 0;
  if (device_entity)
  {
    std::string paramName = Util::paramNameFromFullParamName(name);

    if (paramName == "")
    {
      available = device_entity->analyticMatrixSensitivityAvailableDefaultParam();
    }
    else
    {
      available = device_entity->analyticMatrixSensitivityAvailable(paramName);
    }
  }

  Parallel::AllReduce(comm_, MPI_LOR, &available, 1);

  return available != 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::numericalMatrixSensitivitiesAvailable
// Purpose       :
// Special Notes : Includes a reduce, so works in parallel.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool DeviceMgr::numericalMatrixSensitivitiesAvailable(const std::string & name)
{
  DeviceEntity * device_entity = getDeviceEntity(name);

  // if we find a valid device entity then device-level numerical derivatives will 
  // work.  Don't need to check anything else
  int available = 0;
  if (device_entity)
  {
    available = 1;
  }
  Parallel::AllReduce(comm_, MPI_LOR, &available, 1);

  return available != 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::analyticBVecSensAvailable
// Purpose       :
// Special Notes : Includes a reduce, so works in parallel.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool DeviceMgr::analyticBVecSensAvailable(const std::string & name)
{
  // non-artificial parameters.
  DeviceEntity * device_entity = getDeviceEntity(name);

  int available = 0;
  if (device_entity)
  {
    std::string paramName = Util::paramNameFromFullParamName(name);

    if (paramName == "")
    {
      available = device_entity->getAnalyticACSensitivityAvailableDefaultParam();
    }
    else
    {
      available = device_entity->getAnalyticACSensitivityAvailable(paramName);
    }
  }

  Parallel::AllReduce(comm_, MPI_LOR, &available, 1);

  return available != 0;

}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::numericalBVecSensAvailable
// Purpose       :
// Special Notes : Includes a reduce, so works in parallel.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool DeviceMgr::numericalBVecSensAvailable(const std::string & name)
{
#if 0
  // This assumes that the only available analytic sensitivities are for
  // non-artificial parameters.
  DeviceEntity * device_entity = getDeviceEntity(name);

  // if we find a valid device entity then device-level numerical derivatives will 
  // work.  Don't need to check anything else
  int available = 0;
  if (device_entity)
  {
    available = 1;
  }
  Parallel::AllReduce(comm_, MPI_LOR, &available, 1);

  return available != 0;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getAnalyticalBSensVectorsforAC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void DeviceMgr::getAnalyticalBSensVectorsforAC (const std::string & name,
          std::vector< std::complex<double> > &     dbdp,
          std::vector<int> &        BindicesVec)
{
  // If not artificial, then search for the appropriate natural param(s).
  DeviceEntity * device_entity = getDeviceEntity(name);

  bool found=false;
  if (device_entity)
  {
    std::string paramName = Util::paramNameFromFullParamName(name);
    if (paramName == "")
    {
      found = device_entity->getAnalyticBSensVectorsforACDefaultParam( 
          dbdp, BindicesVec);
    }
    else
    {
      found = device_entity->getAnalyticBSensVectorsforAC(
          paramName, dbdp, BindicesVec);
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getNumericalBSensVectorsforAC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void DeviceMgr::getNumericalBSensVectorsforAC (const std::string & name,
          std::vector< std::complex<double> > &     dbdp,
          std::vector<int> &        BindicesVec)
{
  DeviceEntity * device_entity = getDeviceEntity(name);

  bool found=false;
  if (device_entity)
  {
    std::string paramName = Util::paramNameFromFullParamName(name);
    if (paramName == "")
    {
      found = device_entity->getNumericalBSensVectorsforACDefaultParam(
          dbdp, BindicesVec);
    }
    else
    {
      found = device_entity->getNumericalBSensVectorsforAC(
          paramName, dbdp, BindicesVec);
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getAnalyticMatrixSensitivities
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void DeviceMgr::getAnalyticMatrixSensitivities
  (const std::string & name,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs)
{
  // If not artificial, then search for the appropriate natural param(s).
  DeviceEntity * device_entity = getDeviceEntity(name);

  bool found=false;
  if (device_entity)
  {
    std::string paramName = Util::paramNameFromFullParamName(name);
    if (paramName == "")
    {
      found = device_entity->getAnalyticMatrixSensitivityDefaultParam(
                d_dfdx_dp, d_dqdx_dp, F_lids, Q_lids, F_jacLIDs, Q_jacLIDs);
    }
    else
    {
      found = device_entity->getAnalyticMatrixSensitivity(paramName,
                d_dfdx_dp, d_dqdx_dp, F_lids, Q_lids, F_jacLIDs, Q_jacLIDs);
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getNumericalMatrixSensitivities
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 
//-----------------------------------------------------------------------------
void DeviceMgr::getNumericalMatrixSensitivities
  (const std::string & name,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs)
{

  // what kind of parameter is this:  instance? model?  global_param?
  //
  // instance only needs the same function args as getAnalyticMatrixSensitivity
  //
  // model, however, needs a list of instances.
  //
  // global_param needs even more work.  Currently not supported here.
  //

  DeviceEntity * device_entity = getDeviceEntity(name);

  bool found=false;
  if (device_entity)
  {
    std::string paramName = Util::paramNameFromFullParamName(name);

    if (paramName == "")
    {
      found = device_entity->getNumericalMatrixSensitivityDefaultParam(
                d_dfdx_dp, d_dqdx_dp, F_lids, Q_lids, F_jacLIDs, Q_jacLIDs);
    }
    else
    {
      found = device_entity->getNumericalMatrixSensitivity(paramName,
                d_dfdx_dp, d_dqdx_dp, F_lids, Q_lids, F_jacLIDs, Q_jacLIDs);
    }
  }

  return;
}
//

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setParam
//
// Purpose       : This function sets named parameters (name) to a
//                 specified value (val).
//
// Special Notes : Used for continuation calculations, as well as possibly
//                 intrusive sensitivity/optimization calculations.  It is
//                 assumed that this is called after everything (devices,
//                 solvers, etc.) is set up.
//
//                 The specified parameter can be either a natural or
//                 artificial parameter.
//
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/02/03
//-----------------------------------------------------------------------------
bool DeviceMgr::setParam(
  const std::string &   name,
  double                val,
  bool                  overrideOriginal)
{
  return setParameter(comm_, artificialParameterMap_, passthroughParameterSet_, globals_, *this,
                      dependentPtrVec_, getDevices(ExternDevice::Traits::modelType()), name, val, overrideOriginal);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setParamRandomExpressionTerms
//
// Purpose       : 
// 
// Special Notes : AGAUSS, GAUSS, AUNIF, UNIF, RAND and LIMIT
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 7/30/2020
//-----------------------------------------------------------------------------
bool DeviceMgr::setParamRandomExpressionTerms(
  const std::string &   name, const std::string &   opName, int opIndex,
  //enum Util::astRandTypes astType,
  int astType,
  double                val,
  bool                  overrideOriginal)
{
  return setParameterRandomExpressionTerms(comm_, artificialParameterMap_, passthroughParameterSet_, globals_, *this,
      dependentPtrVec_, getDevices(ExternDevice::Traits::modelType()), name, opName, opIndex, astType, val, overrideOriginal);
  return true;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::findParam
// Purpose       : Returns the current value of a named parameter.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DeviceMgr::parameterExists(
  Parallel::Machine     comm,
  const std::string &   name) const
{
  Util::ParamList param_list;
  param_list.push_back(Param(name, ""));
  Util::ParamList::const_iterator it = param_list.begin();
  
  Util::Op::Operator *op = opBuilderManager_.createOp(it);

  delete op;

  int found = op ? 1 : 0;

  Parallel::AllReduce(comm, MPI_SUM, &found, 1);

  return found > 0;
  
//  return getOp(comm, name);
  
  // double value = 0.0;
  // return getParameter(artificialParameterMap_, globals_.global_params, *this, *measureManager_, name, value);
}


// this macro is to support the populateSweepParam function.
// The OP argument is only used for the debug output part of this.  The functional part only needs FUNC
#define GET_RANDOM_OP_DATA(FUNC,OP) \
  { \
    std::vector<Xyce::Analysis::SweepParam> tmpParams; expr.FUNC(tmpParams);  \
    if ( !(tmpParams.empty()) )  {\
      if (DEBUG_DEVICE) std::cout << "Parameter " << paramName << " contains "<< #OP << std::endl;  \
      for(int jj=0;jj<tmpParams.size();jj++) {  \
        tmpParams[jj].baseName = paramName; \
        tmpParams[jj].name = paramName ; } \
      SamplingParams.insert( SamplingParams.end(), tmpParams.begin(), tmpParams.end() );\
    } \
    else  \
    {  \
      if (DEBUG_DEVICE) std::cout << "Parameter " << paramName << " does not contain " << #OP << std::endl; \
    }\
  }

//-----------------------------------------------------------------------------
// Function      : populateSweepParam
//
// Purpose       : populates a vector of "sweep param" objects for a single 
//                 expression/parameter
//
// Special Notes : Currently, each type of random operator has a separate 
//                 unique "get" function in the expression library.  
//
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 7/29/2020
//-----------------------------------------------------------------------------
void populateSweepParam (  
    Util::Expression & expr, 
    const std::string & paramName, 
    std::vector<Xyce::Analysis::SweepParam> & SamplingParams)
{
  GET_RANDOM_OP_DATA(getAgaussData,Agauss)
  GET_RANDOM_OP_DATA(getGaussData,Gauss)
  GET_RANDOM_OP_DATA(getAunifData,Aunif)
  GET_RANDOM_OP_DATA(getUnifData,Unif)
  GET_RANDOM_OP_DATA(getRandData,Rand)
  GET_RANDOM_OP_DATA(getLimitData,Limit)
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getRandomParams
//
// Purpose       : To conduct Hspice-style sampling analysis, it is necessary 
//                 to gather all the parameters (global and/or device) which 
//                 depend upon expressions containing operators such as 
//                 AGAUSS, GAUSS, AUNIF, UNIF, RAND and LIMIT
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/28/2020
//-----------------------------------------------------------------------------
void DeviceMgr::getRandomParams(std::vector<Xyce::Analysis::SweepParam> & SamplingParams,
   Parallel::Communicator & parallel_comm)
{
  std::vector<Util::Expression> & global_expressions  = globals_.global_expressions;
  std::vector<std::string> & global_exp_names = globals_.global_exp_names;

  // get random expressions from global parameters first.
  // ERK Q: in parallel, is the order of global params preserved?  This global param 
  // part works in parallel when there is a single parameter.
  // If the order is preserved, then nothing extra has to be done for parallel.
  for (int ii=0;ii<global_expressions.size();ii++)
  {
    populateSweepParam(global_expressions[ii], global_exp_names[ii], SamplingParams );
  }

  // now do device params.  When this function is called, 
  // the "dependentPtrVec_" hasn't been set up yet, so it can't be used.
  //
  // This part currently doensn't work in parallel.   It is just considering local
  // params, which means that each proc will have a different list.
  {
    std::vector<Xyce::Analysis::SweepParam>  deviceSamplingParams;

    // models
    for (int ii=0;ii<modelVector_.size();ii++)
    {
      if (!(*(modelVector_[ii])).getDependentParams().empty())
      {
        const std::vector<Depend> & depVec = (*(modelVector_[ii])).getDependentParams();
        std::string entityName = (*(modelVector_[ii])).getName();
        for (int jj=0;jj<depVec.size();jj++) 
        {
          std::string fullName = entityName + ":" + depVec[jj].name;
          populateSweepParam( *(depVec[jj].expr), fullName, deviceSamplingParams ); 
        }
      }
    }

    // instances
    for (int ii=0;ii<instancePtrVec_.size();ii++)
    {
      if (!(*(instancePtrVec_[ii])).getDependentParams().empty())
      {
        const std::vector<Depend> & depVec = (*(instancePtrVec_[ii])).getDependentParams();
        const InstanceName & instName = (*(instancePtrVec_[ii])).getName();
        const std::string entityName = instName.getEncodedName(); // ERK. Check this!
        for (int jj=0;jj<depVec.size();jj++) 
        { 
          std::string fullName = entityName + ":" + depVec[jj].name;
          populateSweepParam( *(depVec[jj].expr), fullName, deviceSamplingParams ); 
        }
      }
    }

    // if parallel, get a combined vector of random params from all processors.
    if (Parallel::is_parallel_run(parallel_comm.comm()))
    {
      int procID = parallel_comm.procID();
      int numProc = parallel_comm.numProc();
      int numDevParam = deviceSamplingParams.size();
      int numDevParamTotal = 0;
      parallel_comm.sumAll(&numDevParam, &numDevParamTotal, 1);

      if (numDevParamTotal > 0)
      {
        std::vector<Xyce::Analysis::SweepParam> orig_deviceSamplingParams = deviceSamplingParams;

        deviceSamplingParams.resize(numDevParamTotal);

        int loc=0;
        for (int proc=0; proc<numProc; ++proc)
        {
          int cnt = 0;
          if (proc == procID) { cnt = numDevParam; }
          parallel_comm.bcast(&cnt, 1, proc);

          for (int i = 0; i < cnt; ++i)
          {
            if (proc == procID)
            {
              Xyce::Analysis::SweepParam & tmpSweepParam = deviceSamplingParams[loc];
              int size = Xyce::packedByteCount(tmpSweepParam);
              int bufSize = size + 100;
              char *buf = new char[bufSize];
              parallel_comm.bcast(&size, 1, proc);
              int pos = 0;
              Xyce::pack(tmpSweepParam, buf, bufSize, pos, &parallel_comm);
              parallel_comm.bcast(buf, size, proc);
              deviceSamplingParams[loc] = orig_deviceSamplingParams[i];
              delete[] buf;
            }
            else
            {
              int size = 0;
              parallel_comm.bcast(&size, 1, proc);
              int bufSize = size + 100;
              char *buf = new char[bufSize];
              parallel_comm.bcast(buf, size, proc);
              int pos = 0;
              Xyce::Analysis::SweepParam tmpSweepParam;
              Xyce::unpack(tmpSweepParam, buf, bufSize, pos, &parallel_comm);
              deviceSamplingParams[loc] = tmpSweepParam;
              delete[] buf;
            }
            ++loc;
          }
        }
      }
    }

    // now insert the model/instance parameters into the vector.
    SamplingParams.insert
      (SamplingParams.end(), deviceSamplingParams.begin(), deviceSamplingParams.end());
  }

  if ( !(SamplingParams.empty()) ) expressionBasedSamplingEnabled_ = true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateDependentParams
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/11/2020
//-----------------------------------------------------------------------------
void DeviceMgr::updateDependentParams()
{
  updateDependentParameters_();
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::resetScaledParams()
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/11/2020
//-----------------------------------------------------------------------------
void DeviceMgr::resetScaledParams()
{
  ModelVector::iterator iterM;
  ModelVector::iterator beginM =modelVector_.begin();
  ModelVector::iterator endM =modelVector_.end();
  for (iterM=beginM; iterM!=endM;++iterM)
  {
    if (!(*iterM)->getDependentParams().empty())
    {
      (*iterM)->resetScaledParams();
    }
  }

  // do the instances
  InstanceVector::iterator iter;
  InstanceVector::iterator begin =instancePtrVec_.begin();
  InstanceVector::iterator end =instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    if (!(*iter)->getDependentParams().empty())
    {
      (*iter)->resetScaledParams();
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateState
// Purpose       : This should be called prior to loadDAEVectors.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/17/05
//-----------------------------------------------------------------------------
bool DeviceMgr::updateState(
  Linear::Vector * nextSolVectorPtr,
  Linear::Vector * currSolVectorPtr,
  Linear::Vector * lastSolVectorPtr,
  Linear::Vector * nextStaVectorPtr,
  Linear::Vector * currStaVectorPtr,
  Linear::Vector * lastStaVectorPtr,
  Linear::Vector * nextStoVectorPtr,
  Linear::Vector * currStoVectorPtr,
  int loadType)
{
  bool bsuccess = true;
  bool tmpBool = true;

  bool all_devices_converged = allDevicesConverged(comm_);

  tmpBool = setupSolverInfo(solState_, *analysisManager_, all_devices_converged, devOptions_, nlsMgrPtr_->getNonLinInfo());
  bsuccess = bsuccess && tmpBool;

  // copy over the passed pointers:
  externData_.nextSolVectorPtr = nextSolVectorPtr;
  externData_.currSolVectorPtr = currSolVectorPtr;
  externData_.lastSolVectorPtr = lastSolVectorPtr;
  externData_.nextStaVectorPtr = nextStaVectorPtr;
  externData_.currStaVectorPtr = currStaVectorPtr;
  externData_.lastStaVectorPtr = lastStaVectorPtr;
  externData_.nextStoVectorPtr = nextStoVectorPtr;
  externData_.currStoVectorPtr = currStoVectorPtr;

#ifdef Xyce_PARALLEL_MPI
  externData_.nextSolVectorPtr->importOverlap();
#endif

  // Now reset the relevant RAW pointers:
  externData_.nextSolVectorRawPtr = &((*externData_.nextSolVectorPtr)[0]);
  externData_.currSolVectorRawPtr = &((*externData_.currSolVectorPtr)[0]);
  externData_.lastSolVectorRawPtr = &((*externData_.lastSolVectorPtr)[0]);
  externData_.nextStaVectorRawPtr = &((*externData_.nextStaVectorPtr)[0]);
  externData_.currStaVectorRawPtr = &((*externData_.currStaVectorPtr)[0]);
  externData_.lastStaVectorRawPtr = &((*externData_.lastStaVectorPtr)[0]);
  externData_.nextStoVectorRawPtr = &((*externData_.nextStoVectorPtr)[0]);
  externData_.currStoVectorRawPtr = &((*externData_.currStoVectorPtr)[0]);

  updateDependentParameters_();

  std::vector<Device*>::iterator iter;
  std::vector<Device*>::iterator begin;
  std::vector<Device*>::iterator end;

  // if inner problem, only do the PDE updates.
  if (loadType == Xyce::Device::PDE)
  {
    begin = pdeDevicePtrVec_.begin ();
    end   = pdeDevicePtrVec_.end ();
  }
  else
  {
    begin = devicePtrVec_.begin ();
    end   = devicePtrVec_.end ();
  }

  for (iter=begin; iter!=end;++iter)
  {
    tmpBool = (*iter)->updateState (externData_.nextSolVectorRawPtr, externData_.nextStaVectorRawPtr, externData_.nextStoVectorRawPtr, loadType);
    bsuccess = bsuccess && tmpBool;
  }

  updateExternalDevices_();

#ifdef Xyce_PARALLEL_MPI
  externData_.nextStaVectorPtr->importOverlap();
  externData_.nextStoVectorPtr->importOverlap();
#endif

  Report::safeBarrier(comm_);

  return true;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadDAEMatrices
// Purpose       : This function loads the various DAE related matrices to
//                 set up the following expression:
//
//                 residual:  f(x) = dQ/dt + F(x) - B(t) = 0
//
//                 jacobian:  J(x) = d(dQ/dt)dx + dFdx
//                                 = d(dQdx)dt + dFdx
//
// Special Notes : ERK.  This function is only called if using the new
//                 (new-DAE) integrator.  As such, it is also used for
//                 MPDE.  It is not used for the old (ODE) integrator.
//                 The corresponding function for the old integrator
//                 is loadJacobianMatrix
//
//                 Note that this function, unlike the loadJacobianMatrix
//                 function *requires* that vector/matrix pointers
//                 be passed in to the function.  When
//                 running in new-DAE or MPDE, we can't rely (much) on
//                 registered vectors. (set up in initializeAll).  In
//                 particular, MPDE cycles through different vector and
//                 matrix blocks, which each correspond to different
//                 "fast" time points, and for each block the state
//                 information is different.
//
//                 NOTE:  This function assumes that all the load matrices
//                 are zeroed out.  That is, dFdx, dQdx are all zeros.
//
//                 If that isn't true, then this function will not produce
//                 the correct answer.  The reason for doing this
//                 is MPDE.  For the warped case, the entire linear system
//                 needs to be zeroed out, but the full system includes
//                 phase equations that never get passed in here.
//                 That zeroing now happens upstream from here, in either
//                 the MPDE loader or the time integrator.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
bool DeviceMgr::loadDAEMatrices(
  Linear::Vector * tmpSolVectorPtr,
  Linear::Vector * tmpStateVectorPtr,
  Linear::Vector * tmpStateDerivVectorPtr,
  Linear::Vector * tmpStoreVectorPtr,
  Linear::Matrix * tmpdQdxMatrixPtr,
  Linear::Matrix * tmpdFdxMatrixPtr,
  int loadType)
{
  bool bsuccess = true;

  // copy over the passed pointers:
  externData_.nextSolVectorPtr = tmpSolVectorPtr;

#ifdef Xyce_PARALLEL_MPI
  externData_.nextSolVectorPtr->importOverlap();
#endif

  bool resetRawMatrixPointers = true;
  if ((externData_.dQdxMatrixPtr == tmpdQdxMatrixPtr) && (externData_.dFdxMatrixPtr == tmpdFdxMatrixPtr))
  {
    resetRawMatrixPointers = false;
  }
  if (resetRawMatrixPointers)
  {
    externData_.dQdxMatrixPtr = tmpdQdxMatrixPtr;
    externData_.dFdxMatrixPtr = tmpdFdxMatrixPtr;
  }

  externData_.nextStaVectorPtr = tmpStateVectorPtr;
  externData_.nextStaDerivVectorPtr = tmpStateDerivVectorPtr;
  externData_.nextStoVectorPtr = tmpStoreVectorPtr;

  // setup the relevant RAW vector pointers:
  externData_.nextSolVectorRawPtr = &((*externData_.nextSolVectorPtr)[0]);
  externData_.nextStaVectorRawPtr = &((*externData_.nextStaVectorPtr)[0]);
  externData_.nextStaDerivVectorRawPtr = &((*externData_.nextStaDerivVectorPtr)[0]);
  externData_.nextStoVectorRawPtr = &((*externData_.nextStoVectorPtr)[0]);

  // setup the relevant RAW matrix pointers (down in the devices that need them):
  if (resetRawMatrixPointers || solState_.blockAnalysisFlag_)
  {
    setupRawMatrixPointers_();
  }

  // if this is an "inner problem" phase of a Two-Level Newton
  // simulation, then only load the PDE devices.  Everything else just
  // gets "1" on the diagonal.
  //
  // Note, it is possible to just load 1's using a petra call, so I may
  // get rid of the "trivial" matrix stamp stuff soon.

  if (loadType == Xyce::Device::PDE)
  {
    InstanceVector::iterator iter = nonPdeInstancePtrVec_.begin(); 
    InstanceVector::iterator end = nonPdeInstancePtrVec_.end();
    for ( ; iter!=end; ++iter)
    {
      (*iter)->loadTrivialDAE_FMatrixStamp ();
    }

    iter = pdeInstancePtrVec_.begin(); 
    end = pdeInstancePtrVec_.end();
    for ( ; iter!=end; ++iter)
    {
      bsuccess = bsuccess && (*iter)->loadDAEdQdx();
      bsuccess = bsuccess && (*iter)->loadDAEdFdx();
    }
  }
  // Else, do a normal analytical matrix load.
  else
  {
    for (DeviceVector::iterator it = devicePtrVec_.begin(), end = devicePtrVec_.end(); it != end; ++it)
    {
      bsuccess = bsuccess && (*it)->loadDAEMatrices(*externData_.dFdxMatrixPtr , *externData_.dQdxMatrixPtr, loadType);
    }
  }

  // Run jacobian diagnostic.
  if (devOptions_.testJacobianFlag &&
      (solState_.timeStepNumber_ >= devOptions_.testJacStartStep &&
       solState_.timeStepNumber_ <= devOptions_.testJacStopStep))
  {
    // Test just the specified device(s), if the user specified any.
    if ((devOptions_.testJacDeviceNameGiven))
    {
      InstanceVector::iterator iter = testJacDevicePtrVec_.begin(), end = testJacDevicePtrVec_.end();
      for ( ; iter != end; ++iter)
      {
        if ( (( loadType == LINEAR ) && (*iter)->isLinearDevice())
             || (( loadType == NONLINEAR ) && !(*iter)->isLinearDevice()) 
             || (( loadType == PDE ) && (*iter)->isPDEDevice())
             || (loadType == ALL) )
        {
          (*iter)->testDAEMatrices(topology_.getSolutionNodeNames());
        }
      }
    }
    else // Test all the devices:
    {
      InstanceVector::iterator iter = instancePtrVec_.begin(), end = instancePtrVec_.end();
      for ( ; iter != end; ++iter)
      {
        if ( (( loadType == LINEAR ) && (*iter)->isLinearDevice())
             || (( loadType == NONLINEAR ) && !(*iter)->isLinearDevice())
             || (( loadType == PDE ) && (*iter)->isPDEDevice())
             || ( loadType == ALL ) )
        {
          (*iter)->testDAEMatrices(topology_.getSolutionNodeNames());
        }
      }
    }
  }

  //Tell Jacobian, fill is complete allowing accumulation if necessary
  //externData_.dQdxMatrixPtr->fillComplete();
  //externData_.dFdxMatrixPtr->fillComplete();

  Report::safeBarrier(comm_);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PRINT_VECTORS) && solState_.debugTimeFlag)
  {
    int newtonIter = solState_.newtonIter;
    dout() << section_divider << std::endl;
    dout() <<  "Q-matrix: nonlinear iteration = " << newtonIter << "\n";
    externData_.dQdxMatrixPtr->printPetraObject(dout());
    dout() << std::endl;
    dout() << section_divider << std::endl;
    dout() <<  "F-matrix: nonlinear iteration = " << newtonIter << "\n";
    externData_.dFdxMatrixPtr->printPetraObject(dout());
    dout() << std::endl;
    dout() << section_divider << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadDAEVectors
// Purpose       : This function loads the various DAE related vectors to
//                 set up the following expression:
//
//                 f(x) = dQ/dt + F(x) - B(t) = 0
//
// Special Notes : ERK.  This function is only called if using the new
//                 (new-DAE) integrator.  As such, it is also used for
//                 MPDE.  It is not used for the old (ODE) integrator.
//                 The corresponding function for the old integrator
//                 is loadRHSVector.
//
//                 Note that this function, unlike the loadRHSVector
//                 function *requires* that vectors be passed in.  When
//                 running in new-DAE or MPDE, we can't rely (much) on
//                 registered vectors. (set up in initializeAll).  In
//                 particular, MPDE cycles through different vector and
//                 matrix blocks, which each correspond to different
//                 "fast" time points, and for each block the state
//                 information is different.
//
//                 NOTE:  This function assumes that all the load vectors
//                 are zeroed out.  That is, F, Q, B, dFdxdVp and dQdxVp
//                 are all zeros.
//
//                 If that isn't true, then this function will not produce
//                 the correct answer.  The reason for doing this
//                 is MPDE.  For the warped case, the entire linear system
//                 needs to be zeroed out, but the full system includes
//                 phase equations that never get passed in here.
//
//                 That zeroing now happens upstream from here, in either
//                 the MPDE loader or the time integrator.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/26/03
//-----------------------------------------------------------------------------
bool DeviceMgr::loadDAEVectors(
  Linear::Vector * tmpSolVectorPtr,
  Linear::Vector * tmpCurrSolVectorPtr,
  Linear::Vector * tmpLastSolVectorPtr,
  Linear::Vector * tmpStaVectorPtr,
  Linear::Vector * tmpCurrStaVectorPtr,
  Linear::Vector * tmpLastStaVectorPtr,
  Linear::Vector * tmpStaDerivVectorPtr,
  Linear::Vector * tmpStoVectorPtr,
  Linear::Vector * tmpCurrStoVectorPtr,
  Linear::Vector * tmpLeadFCompVectorPtr,
  Linear::Vector * tmpLeadQCompVectorPtr,
  Linear::Vector * tmpJunctionVCompVectorPtr,
  Linear::Vector * tmpQVectorPtr,
  Linear::Vector * tmpFVectorPtr,
  Linear::Vector * tmpBVectorPtr,
  Linear::Vector * tmpdFdxdVpVectorPtr,
  Linear::Vector * tmpdQdxdVpVectorPtr,
  int loadType)
{
  bool bsuccess = true;
  bool tmpBool = true;

  // copy over the passed pointers:
  externData_.nextSolVectorPtr = tmpSolVectorPtr;
  externData_.currSolVectorPtr = tmpCurrSolVectorPtr;
  externData_.lastSolVectorPtr = tmpLastSolVectorPtr;
  externData_.daeQVectorPtr    = tmpQVectorPtr;
  externData_.daeFVectorPtr    = tmpFVectorPtr;
  externData_.daeBVectorPtr    = tmpBVectorPtr;
  externData_.dFdxdVpVectorPtr = tmpdFdxdVpVectorPtr;
  externData_.dQdxdVpVectorPtr = tmpdQdxdVpVectorPtr;
  externData_.nextStaVectorPtr = tmpStaVectorPtr;
  externData_.currStaVectorPtr = tmpCurrStaVectorPtr;
  externData_.lastStaVectorPtr = tmpLastStaVectorPtr;
  externData_.nextStaDerivVectorPtr = tmpStaDerivVectorPtr;
  externData_.nextStoVectorPtr = tmpStoVectorPtr;
  externData_.currStoVectorPtr = tmpCurrStoVectorPtr;
  externData_.nextLeadCurrFCompPtr = tmpLeadFCompVectorPtr;
  externData_.nextLeadCurrQCompPtr = tmpLeadQCompVectorPtr;
  externData_.nextJunctionVCompPtr = tmpJunctionVCompVectorPtr;

  // Make sure all boundary data is valid in the solution vector
#ifdef Xyce_PARALLEL_MPI
  externData_.nextSolVectorPtr->importOverlap();
  externData_.nextStaDerivVectorPtr->importOverlap();
#endif

  // Set up the relevant RAW Pointers:
  setupRawVectorPointers_ ();

  // call all the intermediate vars loads:
  std::vector<Device*>::iterator iter;
  std::vector<Device*>::iterator begin;
  std::vector<Device*>::iterator end;

  // if inner problem, only do the PDE loads.
  if (loadType == Xyce::Device::PDE)
  {
    begin = pdeDevicePtrVec_.begin ();
    end   = pdeDevicePtrVec_.end ();
  }
  else
  {
    begin = devicePtrVec_.begin ();
    end   = devicePtrVec_.end ();
  }

  for (iter=begin; iter!=end;++iter)
  {
    bsuccess &= (*iter)->updateSecondaryState
                (externData_.nextStaDerivVectorRawPtr, externData_.nextStoVectorRawPtr);
  }

  // I'M NOT SURE ABOUT THIS ONE.
#ifdef Xyce_PARALLEL_MPI
  externData_.nextStaVectorPtr->importOverlap();
  externData_.nextStoVectorPtr->importOverlap();
#endif

  for (iter=begin; iter!=end;++iter)
  {
    bsuccess=(*iter)->loadDAEVectors(externData_.nextSolVectorRawPtr,
                                     externData_.daeFVectorRawPtr,
                                     externData_.daeQVectorRawPtr,
                                     externData_.daeBVectorRawPtr,
                                     externData_.nextLeadCurrFCompRawPtr,
                                     externData_.nextLeadCurrQCompRawPtr,
                                     externData_.nextJunctionVCompRawPtr,
                                     loadType);
  }

  // dump to the screen:
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PRINT_VECTORS) && solState_.debugTimeFlag)
  {
    int newtonIter = solState_.newtonIter;
    dout() <<  "Q-vector: nonlinear iteration = " << newtonIter << "\n";
    externData_.daeQVectorPtr->printPetraObject(std::cout);
    dout() << std::endl;
    dout() <<  "F-vector: nonlinear iteration = " << newtonIter << "\n";
    externData_.daeFVectorPtr->printPetraObject(std::cout);
    dout() << std::endl;

    if (devOptions_.voltageLimiterFlag)
    {
      dout() << "\n\n  dFdxdVp vector: nonlinear iteration = " << newtonIter << "\n";
      externData_.dFdxdVpVectorPtr->printPetraObject(std::cout);
      dout() << std::endl;
      dout() << "\n\n  dQdxdVp vector: nonlinear iteration = " << newtonIter << "\n";
      externData_.dQdxdVpVectorPtr->printPetraObject(std::cout);
      dout() << std::endl;
    }

    debugOutput2();
  }

  // Update parallel if necessary
  //externData_.daeQVectorPtr->fillComplete();
  //externData_.daeFVectorPtr->fillComplete();
  //externData_.daeBVectorPtr->fillComplete();
  //externData_.dFdxdVpVectorPtr->fillComplete();
  //externData_.dQdxdVpVectorPtr->fillComplete();

  Report::safeBarrier(comm_);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateFDIntermediateVars
// Purpose       : Let devices update internal frequency domain quantities
//
// Scope         : public
// Creator       : Tom Russo and Ting Mei
// Creation Date : 19 March 2018
//-----------------------------------------------------------------------------
/// Update device internal frequency domain variables
///
/// \author Tom Russo and Ting Mei
/// \date 19 March 2018
///
bool DeviceMgr::updateFDIntermediateVars(double frequency,
                                         std::complex<double>* freqSolVec)
{
  bool bsuccess(true);

  std::vector<Device*>::iterator iter;
  std::vector<Device*>::iterator begin;
  std::vector<Device*>::iterator end;

  begin = devicePtrVec_.begin ();
  end   = devicePtrVec_.end ();

  for (iter=begin; iter!=end; ++iter)
  {
    bsuccess &= (*iter)->updateFDIntermediateVars(frequency, freqSolVec);
  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadFreqDAEMatrices
// Purpose       : Load frequency domain contributions into the Jacobian
//
// Scope         : public
// Creator       : Jason C. Verley
// Creation Date : 07/18/17
//-----------------------------------------------------------------------------
/// Load frequency domain matrix contributions
///
/// Loads the contributions to the Jacobian matrix.  To be used by
/// frequency-based analyses techniques.  
///
/// \author Jason C. Verley
/// \date 07/18/17
///
bool DeviceMgr::loadFreqDAEMatrices(double frequency, std::complex<double>* freqSolVec,
                                    std::vector<Util::FreqMatEntry>& dFdxMatrix)
{ 
  bool bsuccess(true);

  std::vector<Device*>::iterator iter;
  std::vector<Device*>::iterator begin;
  std::vector<Device*>::iterator end;

  begin = devicePtrVec_.begin ();
  end   = devicePtrVec_.end ();

  for (iter=begin; iter!=end; ++iter)
  {
    bsuccess &= (*iter)->loadFreqDAEMatrices(frequency, freqSolVec, dFdxMatrix);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadFreqDAEVectors
// Purpose       : Load frequency domain contributions into the vectors
//
// Scope         : public
// Creator       : Jason C. Verley
// Creation Date : 07/18/17
//-----------------------------------------------------------------------------
/// Load frequency domain F and B contributions
///
/// Loads the contributions to the `F` vector and `B` vector.  To be used by 
/// frequency-based analyses techniques.  Only devices that have a frequency
/// response will add to the freqFVector and freqBVector.
///
/// \author Jason C. Verley
/// \date 07/18/17
///
bool DeviceMgr::loadFreqDAEVectors(double frequency, std::complex<double>* freqSolVec, 
                                   std::vector<Util::FreqVecEntry>& freqFVector,
                                   std::vector<Util::FreqVecEntry>& freqBVector)
{
  bool bsuccess(true);

  std::vector<Device*>::iterator iter;
  std::vector<Device*>::iterator begin;
  std::vector<Device*>::iterator end;

  begin = devicePtrVec_.begin ();
  end   = devicePtrVec_.end ();

  for (iter=begin; iter!=end; ++iter)
  {
    bsuccess &= (*iter)->loadFreqDAEVectors(frequency, freqSolVec, freqFVector, freqBVector);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadErrorWeightMask ()
// Purpose       : let devices set elements of a mask telling the time
//                 integrator what equations should be ignored in taking
//                 weighted norms for error control purposes.
// Special Notes : Devices should *only* zero internal variables, and then
//                 only those that absolutely should never be used to
//                 control step size (e.g. excess phase variables in BJTs)
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 01/18/07
//-----------------------------------------------------------------------------
bool DeviceMgr::loadErrorWeightMask(Linear::Vector * deviceMask)
{
  externData_.deviceErrorWeightMask_ = deviceMask;

  InstanceVector::const_iterator it = instancePtrVec_.begin(), end = instancePtrVec_.end();
  for ( ; it != end; ++it)
  {
    (*it)->loadErrorWeightMask();
  }

  externData_.deviceErrorWeightMask_->fillComplete();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addGlobalPar()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/18/05
//-----------------------------------------------------------------------------
void DeviceMgr::addGlobalPar(const Util::Param & param)
{
  // get the device manager's current temp
  double temp = getDeviceOptions().temp.getImmutableValue<double>();

  addGlobalParameter(solState_, temp, globals_, param);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::findGlobalPar
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
const double * DeviceMgr::findGlobalPar(
  const std::string & parName) const
{
  return findGlobalParameter(globals_.global_params, parName);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateIntermediateVars_
// Purpose       : This function calls updateIntermediateVars
//                 for the current devices.
// Special Notes :
// Scope         : private
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 12/14/2006
//-----------------------------------------------------------------------------
bool DeviceMgr::updateIntermediateVars_()
{
  bool bsuccess = true;

  InstanceVector::iterator iter;
  InstanceVector::iterator begin =instancePtrVec_.begin();
  InstanceVector::iterator end =instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->updateIntermediateVars ();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updatePrimaryState_
// Purpose       : This function updates primary states
//                 for the present time step.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
bool DeviceMgr::updatePrimaryState_()
{
  bool bsuccess = true;

  InstanceVector::iterator iter;
  InstanceVector::iterator begin =instancePtrVec_.begin();
  InstanceVector::iterator end =instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->updatePrimaryState ();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateSecondaryState_
// Purpose       : This function function updates secondary states
//                 for the present time step.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
bool DeviceMgr::updateSecondaryState_()
{
  bool bsuccess = true;

  InstanceVector::iterator iter;
  InstanceVector::iterator begin =instancePtrVec_.begin();
  InstanceVector::iterator end =instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->updateSecondaryState ();
  }

  return bsuccess;
}

//----------------------------------------------------------------------------
// Function       : DeviceMgr::updateDependentParameters_
// Purpose        : This function updates all dependent parameters for
//                  the current time step.
// Special Notes  : This was evolved from updateTimeDependentParameters_
// Scope          : private
// Creator        : Dave Shirley
// Creation Date  : 08/17/06
//----------------------------------------------------------------------------
void DeviceMgr::updateDependentParameters_()
{
  GlobalParameterMap & globalParamMap = globals_.global_params;
  std::vector<Util::Expression> & globalExpressionsVec = globals_.global_expressions;

  bool timeChanged = false;
  bool freqChanged = false;
  bool globalParamChanged = false;

  if (timeParamsProcessed_ != solState_.currTime_)
    timeChanged = true;

  if (freqParamsProcessed_ != solState_.currFreq_)
    freqChanged = true;

  // Update global params for new time and other global params
  int pos = 0;
  std::vector<Util::Expression>::iterator globalExprIter = globalExpressionsVec.begin(); 
  std::vector<Util::Expression>::iterator globalExprEnd  = globalExpressionsVec.end();
  for ( ; globalExprIter != globalExprEnd; ++globalExprIter)
  {
    // ERK.  4/24/2020. This (the changed bool) was conditionally true/false 
    // depending on various calls such as set_sim_time, which let each expression 
    // report back if its internal time variable (or temp, or freq, etc) had changed.
    //
    // Those function calls don't exist anymore due to the new expression refactor.
    //
    // This logic should maybe be updated to use the same logic that devices do w.r.t things like
    // limiting.  If on a new time step, time changes, otherwise not.  Etc.  
    // If this isn't changed, then a lot of "processParam" function calls are going 
    // to be called unneccessarily.
    //
    // But that will have to come later.

    double val;
    if (globalExprIter->evaluateFunction(val))
    {
      globalParamChanged = true;
      globalParamMap[globals_.global_exp_names[pos]] = val;
    }
    ++pos;
  }

  // do the models:
  if (firstDependent_)
  {
    firstDependent_ = false;

    dependentPtrVec_.clear();

    ModelVector::iterator iterM;
    ModelVector::iterator beginM =modelVector_.begin();
    ModelVector::iterator endM =modelVector_.end();
    for (iterM=beginM; iterM!=endM;++iterM)
    {
      if (!(*iterM)->getDependentParams().empty())
      {
        dependentPtrVec_.push_back(static_cast<DeviceEntity *>(*iterM));
        bool globalParamChangedLocal=true, timeChangedLocal=true, freqChangedLocal=true; // local to this if-stement
        bool tmpBool = (*iterM)->updateGlobalAndDependentParameters(globalParamChangedLocal,timeChangedLocal,freqChangedLocal);
        (*iterM)->processParams();
        (*iterM)->processInstanceParams();
      }
    }

    // do the instances
    InstanceVector::iterator iter;
    InstanceVector::iterator begin =instancePtrVec_.begin();
    InstanceVector::iterator end =instancePtrVec_.end();
    for (iter=begin; iter!=end;++iter)
    {
      if (!(*iter)->getDependentParams().empty())
      {
        dependentPtrVec_.push_back(static_cast<DeviceEntity *>(*iter));
        bool globalParamChangedLocal=true, timeChangedLocal=true, freqChangedLocal=true; // local to this if-stement
        bool tmpBool = (*iter)->updateGlobalAndDependentParameters(globalParamChangedLocal,timeChangedLocal,freqChangedLocal);
        (*iter)->processParams();
      }
    }
  }
  else
  {
    EntityVector::iterator iter;
    EntityVector::iterator begin = dependentPtrVec_.begin();
    EntityVector::iterator end = dependentPtrVec_.end();
    for (iter=begin; iter!=end;++iter)
    {
      bool changed = false;
      changed = (*iter)->updateGlobalAndDependentParameters(globalParamChanged, timeChanged, freqChanged);

      if (changed)
      {
        (*iter)->processParams();
        (*iter)->processInstanceParams();
      }
    }
  }
  timeParamsProcessed_ = solState_.currTime_;

  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadBVectorsforAC
// Purpose       : This function loads the B-vector contributions for sources.
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool DeviceMgr::loadBVectorsforAC(Linear::Vector * bVecRealPtr, Linear::Vector * bVecImagPtr)
{
  bool bsuccess = true;

  // make sure that the correct value of "currFreq" is set in solver state, 
  // which will then get propagated to expressions.  This is to allow expressions to be
  // dependent on freqency.
  bool all_devices_converged = allDevicesConverged(comm_);
  bool tmpBool = setupSolverInfo(solState_, *analysisManager_, all_devices_converged, devOptions_, nlsMgrPtr_->getNonLinInfo());
  // to make sure global parameters are updated properly in AC analysis.
  updateDependentParameters_();

#ifdef Xyce_PARALLEL_MPI
  externData_.nextSolVectorPtr->importOverlap();
#endif
  
  // Now reset the relevant RAW pointers:
  double * bVecRealRawPtr = &((*bVecRealPtr)[0]);
  double * bVecImagRawPtr = &((*bVecImagPtr)[0]);

  IndependentSourceVector::iterator it = independentSourceVector_.begin();
  IndependentSourceVector::iterator end = independentSourceVector_.end();
  for ( ; it != end; ++it)
  {
    (*it)->loadBVectorsforAC(bVecRealRawPtr, bVecImagRawPtr);
  }

  bVecRealPtr->fillComplete();

  bVecImagPtr->fillComplete();

  return bsuccess;  

}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadBVectorsforSources
// Purpose       : This function loads the B-vector contributions for sources.
// Special Notes : This is only done for linear sources.
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool DeviceMgr::loadBVectorsforSources()
{
  bool bsuccess = true;

  IndependentSourceVector::iterator it = independentSourceVector_.begin();
  IndependentSourceVector::iterator end = independentSourceVector_.end();
  for ( ; it != end; ++it)
  {
    if ( (*it)->isLinearDevice() )
    {
      (*it)->loadDAEBVector();
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getNumNoiseDevices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
int DeviceMgr::getNumNoiseDevices()
{
  int totalNum=0;

  if (!(instancePtrVec_.empty()))
  {
    for (InstanceVector::iterator iter = instancePtrVec_.begin();
        iter != instancePtrVec_.end(); ++iter)
    {
      int numNoiseSources=(*iter)->getNumNoiseSources();
      if (numNoiseSources > 0)
      {
        totalNum++;
      }
    }
  }

  return totalNum;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getNumNoiseSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
int DeviceMgr::getNumNoiseSources()
{
  int totalNum=0;
  for (InstanceVector::iterator iter = instancePtrVec_.begin();
      iter != instancePtrVec_.end(); ++iter)
  {
    int numNoiseSources=(*iter)->getNumNoiseSources();
    totalNum+=numNoiseSources;
  }

  return totalNum;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setupNoiseSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
void DeviceMgr::setupNoiseSources
  (std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec)
{
  int i=0;
  for (InstanceVector::iterator iter = instancePtrVec_.begin();
      iter != instancePtrVec_.end(); ++iter)
  {
    int numNoiseSources=(*iter)->getNumNoiseSources();
    if (numNoiseSources>0)
    {
      Xyce::Analysis::NoiseData * ndPtr = noiseDataVec[i];
      Xyce::Analysis::NoiseData & nd = *ndPtr;
      (*iter)->setupNoiseSources ( nd );
      i++;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getNoiseSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
void DeviceMgr::getNoiseSources
  (std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec)
{
  int i=0;
  for (InstanceVector::iterator iter = instancePtrVec_.begin();
      iter != instancePtrVec_.end(); ++iter)
  {
    int numNoiseSources=(*iter)->getNumNoiseSources();
    if (numNoiseSources>0)
    {
      Xyce::Analysis::NoiseData * ndPtr = noiseDataVec[i];
      Xyce::Analysis::NoiseData & nd = *ndPtr;
      (*iter)->getNoiseSources( nd );
      i++;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getBMatrixEntries()
// Purpose       : This function obtains the indices for the B-vector contributions for sources.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool DeviceMgr::getBMatrixEntries(
  std::vector<int> &    bMatEntriesVec,
  std::vector<int> &    portVec,
  std::vector<double> * Z0sVec)

{
  bool bsuccess = true;

  IndependentSourceVector::iterator it = independentSourceVector_.begin();
  IndependentSourceVector::iterator end = independentSourceVector_.end();
  for ( ; it != end; ++it)
  {
    Vsrc::Instance * vsrc = dynamic_cast<Vsrc::Instance *>(*it);
     if (vsrc != 0)
     {
       int lpos, lneg, lbra;

       vsrc->getLIDs(lpos, lneg, lbra);

       if (solState_.spAnalysisFlag_ )
       {
         if (vsrc->isPortSpecified())
         {
           portVec.push_back(vsrc->getPort());        
           (*Z0sVec).push_back(vsrc->getZ0());
           bMatEntriesVec.push_back(lbra);
         }
       }
       else
       {
         portVec.push_back(lpos);
         bMatEntriesVec.push_back(lbra);
       }
     }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateSources
// Purpose       : This function function updates sources for the present
//                 time step.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/26/00
//-----------------------------------------------------------------------------
bool DeviceMgr::updateSources()
{
  bool bsuccess = true;

  bool all_devices_converged = allDevicesConverged(comm_);

  setupSolverInfo(solState_, *analysisManager_, all_devices_converged, devOptions_, nlsMgrPtr_->getNonLinInfo());

  IndependentSourceVector::iterator it = independentSourceVector_.begin();
  IndependentSourceVector::iterator end = independentSourceVector_.end();
  for ( ; it != end; ++it)
  {
    (*it)->updateSource();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setICs
// Purpose       : This function function sets initial conditions for devices
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/13/00
//-----------------------------------------------------------------------------
bool DeviceMgr::setICs(
  Linear::Vector * tmpSolVectorPtr,
  Linear::Vector * tmpCurrSolVectorPtr,
  Linear::Vector * tmpLastSolVectorPtr,
  Linear::Vector * tmpStaVectorPtr,
  Linear::Vector * tmpCurrStaVectorPtr,
  Linear::Vector * tmpLastStaVectorPtr,
  Linear::Vector * tmpStaDerivVectorPtr,
  Linear::Vector * tmpStoVectorPtr,
  Linear::Vector * tmpCurrStoVectorPtr,
  Linear::Vector * tmpQVectorPtr,
  Linear::Vector * tmpFVectorPtr,
  Linear::Vector * tmpBVectorPtr,
  Linear::Vector * tmpdFdxdVpVectorPtr,
  Linear::Vector * tmpdQdxdVpVectorPtr)
{
  bool bsuccess = true;

  // copy over the passed pointers:
  externData_.nextSolVectorPtr = tmpSolVectorPtr;
  externData_.currSolVectorPtr = tmpCurrSolVectorPtr;
  externData_.lastSolVectorPtr = tmpLastSolVectorPtr;
  externData_.daeQVectorPtr    = tmpQVectorPtr;
  externData_.daeFVectorPtr    = tmpFVectorPtr;
  externData_.daeBVectorPtr    = tmpBVectorPtr;
  externData_.dFdxdVpVectorPtr = tmpdFdxdVpVectorPtr;
  externData_.dQdxdVpVectorPtr = tmpdQdxdVpVectorPtr;
  externData_.nextStaVectorPtr = tmpStaVectorPtr;
  externData_.currStaVectorPtr = tmpCurrStaVectorPtr;
  externData_.lastStaVectorPtr = tmpLastStaVectorPtr;
  externData_.nextStaDerivVectorPtr = tmpStaDerivVectorPtr;
  externData_.nextStoVectorPtr = tmpStoVectorPtr;
  externData_.currStoVectorPtr = tmpCurrStoVectorPtr;

  // Make sure all boundary data is valid in the solution vector
#ifdef Xyce_PARALLEL_MPI
  externData_.nextSolVectorPtr->importOverlap();
  externData_.nextStaDerivVectorPtr->importOverlap();
#endif

  // if IC's on devices are set, we need to ensure that the
  // raw pointers are up to date first.
  setupRawVectorPointers_();

  for (InstanceVector::iterator iter = instancePtrVec_.begin();
       iter != instancePtrVec_.end(); ++iter)
  {
    (*iter)->setIC();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::output
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool DeviceMgr::outputPlotFiles(bool force_final_output)
{
  bool bsuccess = true;

  InstanceVector::iterator it = plotFileInstancePtrVec_.begin();
  InstanceVector::iterator end = plotFileInstancePtrVec_.end ();
  for ( ; it != end; ++it)
  {
    bool tmpBool = (*it)->outputPlotFiles(force_final_output);
    bsuccess = bsuccess && tmpBool;
  }

  dotOpOutput();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::finishOutput
// Purpose       : Same as output, only this one forces the output.
//
// Special Notes : This function was neccessary with the implementation of
//                 outputInterval.  The final output, which needs to
//                 happen after the transient is over, won't neccessarily
//                 happen with outputInterval enabled.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/19/04
//-----------------------------------------------------------------------------
bool DeviceMgr::finishOutput()
{
  bool bsuccess = true;
  bool tmpBool = true;

  tmpBool = outputPlotFiles(true);
  bsuccess = bsuccess && tmpBool;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::dotOpOutput
// Purpose       :
// Special Notes : This is a quick-and-dirty implementation, to get something
//                 working quickly.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/3/12
//-----------------------------------------------------------------------------
void DeviceMgr::dotOpOutput()
{
  // only output .OP information once.
  if (dotOpOutputRequested_ && !dotOpOutputFlag_)
  {
    dotOpOutputFlag_ = true;

    std::map<std::string, Device *> device_map;
    for (EntityTypeIdDeviceMap::const_iterator it = deviceMap_.begin();
         it != deviceMap_.end(); ++it)
    {
      device_map[(*it).second->getName()] = (*it).second;
    }

    lout() << section_divider << "\n"
           << "Operating point information:";

    for (std::map<std::string, Device *>::const_iterator it = device_map.begin();
         it != device_map.end(); ++it)
    {
      Xyce::Device::printDotOpOutput(lout(), *(*it).second);
    }

    lout() << section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setGlobalFlags
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 11/2/12
//-----------------------------------------------------------------------------
void DeviceMgr::setGlobalFlags()
{
  int pde_flag = solState_.isPDESystem_ ? 1 : 0;
  Parallel::AllReduce(comm_, MPI_LOR, &pde_flag, 1);
  solState_.isPDESystem_ = pde_flag != 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::resetBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/17/2019
//-----------------------------------------------------------------------------
void DeviceMgr::resetBreakPoints()
{
  std::vector<Util::Expression>::iterator globalExp_i = globals_.global_expressions.begin();
  std::vector<Util::Expression>::iterator globalExp_end = globals_.global_expressions.end();
  for (; globalExp_i != globalExp_end; ++globalExp_i)
  {
    globalExp_i->setupBreakPoints();
  }

  InstanceVector::iterator iter = instancePtrVec_.begin();
  InstanceVector::iterator end = instancePtrVec_.end();
  for (; iter!= end; iter++)
  {
    (*iter)->setupBreakPoints();
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/08/01
//-----------------------------------------------------------------------------
bool DeviceMgr::getBreakPoints (std::vector<Util::BreakPoint> & breakPointTimes,
      std::vector<Util::BreakPoint> & pauseBreakPointTimes)
{
  InstanceVector::iterator iterI;
  ModelVector::iterator iterM;
  bool bsuccess = true;
  bool tmpBool = true;

  bool all_devices_converged = allDevicesConverged(comm_);

  tmpBool = setupSolverInfo(solState_, *analysisManager_, all_devices_converged, devOptions_, nlsMgrPtr_->getNonLinInfo());
  bsuccess = bsuccess && tmpBool;
  setupRawVectorPointers_ ();

  // For some devices we only need to set breakpoints caused by discontinuities
  // in their parameters:

  for (iterM = modelVector_.begin() ; iterM != modelVector_.end() ; ++iterM)
  {
    if (!(*iterM)->getDependentParams().empty())
    {
      tmpBool = (*iterM)->getParamBreakpoints(breakPointTimes);
      bsuccess = bsuccess && tmpBool;
    }
  }

  for (iterI = instancePtrVec_.begin() ; iterI != instancePtrVec_.end() ; ++iterI)
  {
    if (!(*iterI)->getDependentParams().empty())
    {
      tmpBool = (*iterI)->getParamBreakpoints(breakPointTimes);
      bsuccess = bsuccess && tmpBool;
    }
  }

  // Breakpoints for global params:
  std::vector<Util::Expression>::iterator globalExp_i = globals_.global_expressions.begin();
  std::vector<Util::Expression>::iterator globalExp_end = globals_.global_expressions.end();
  for (; globalExp_i != globalExp_end; ++globalExp_i)
  {
    globalExp_i->getBreakPoints(breakPointTimes);
  }

  if (!breakPointInstancesInitialized)
  {
    InstanceVector::iterator beginI = nonpauseBPDeviceVector_.begin();
    InstanceVector::iterator endI = nonpauseBPDeviceVector_.end();

    for (iterI=beginI;iterI!=endI;++iterI)
    {
      // this function returns false if it is the base class, and true otherwise.
      bool functionSetup = (*iterI)->getInstanceBreakPoints (breakPointTimes);
      if (functionSetup)
      {
        simpleBpInstancePtrVec_.push_back(*iterI);
      }
    }
    beginI = pauseBPDeviceVector_.begin();
    endI = pauseBPDeviceVector_.end();
    for (iterI=beginI;iterI!=endI;++iterI)
    {
      // this function returns false if it is the base class, and true otherwise.
      bool functionSetup = (*iterI)->getInstanceBreakPoints (breakPointTimes);
      if (functionSetup)
      {
        pauseBpInstancePtrVec_.push_back(*iterI);
      }
    }
    breakPointInstancesInitialized = true;
  }
  else
  {
    InstanceVector::iterator beginI = simpleBpInstancePtrVec_.begin();
    InstanceVector::iterator endI = simpleBpInstancePtrVec_.end();
    for (iterI=beginI;iterI!=endI;++iterI)
    {
      bool functionSetup = (*iterI)->getInstanceBreakPoints (breakPointTimes);
    }
    beginI = pauseBpInstancePtrVec_.begin();
    endI = pauseBpInstancePtrVec_.end();
    for (iterI=beginI;iterI!=endI;++iterI)
    {
      bool functionSetup = (*iterI)->getInstanceBreakPoints (pauseBreakPointTimes);
    }
  }

  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end())
  {
    for (InstanceVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      extern_device.getBreakPoints(breakPointTimes,pauseBreakPointTimes);
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setupRawVectorPointers_
// Purpose       : set up raw pointers
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
bool DeviceMgr::setupRawVectorPointers_ ()
{
  if (externData_.daeQVectorPtr != 0)
  {
    externData_.daeQVectorRawPtr    = &((*externData_.daeQVectorPtr)[0]);
  }

  if (externData_.daeFVectorPtr != 0)
  {
    externData_.daeFVectorRawPtr    = &((*externData_.daeFVectorPtr)[0]);
  }

  if (externData_.daeBVectorPtr != 0)
  {
    externData_.daeBVectorRawPtr    = &((*externData_.daeBVectorPtr)[0]);
  }

  if (externData_.dFdxdVpVectorPtr != 0)
  {
    externData_.dFdxdVpVectorRawPtr = &((*externData_.dFdxdVpVectorPtr)[0]);
  }

  if (externData_.dQdxdVpVectorPtr != 0)
  {
    externData_.dQdxdVpVectorRawPtr = &((*externData_.dQdxdVpVectorPtr)[0]);
  }

  if (externData_.nextSolVectorPtr != 0)
  {
    externData_.nextSolVectorRawPtr = &((*externData_.nextSolVectorPtr)[0]);
  }

  if (externData_.currSolVectorPtr != 0)
  {
    externData_.currSolVectorRawPtr = &((*externData_.currSolVectorPtr)[0]);
  }

  if (externData_.lastSolVectorPtr != 0)
  {
    externData_.lastSolVectorRawPtr = &((*externData_.lastSolVectorPtr)[0]);
  }

  if (externData_.nextStaVectorPtr != 0)
  {
    externData_.nextStaVectorRawPtr = &((*externData_.nextStaVectorPtr)[0]);
  }

  if (externData_.currStaVectorPtr != 0)
  {
    externData_.currStaVectorRawPtr = &((*externData_.currStaVectorPtr)[0]);
  }

  if (externData_.lastStaVectorPtr != 0)
  {
    externData_.lastStaVectorRawPtr = &((*externData_.lastStaVectorPtr)[0]);
  }

  if (externData_.nextStaDerivVectorPtr != 0)
  {
    externData_.nextStaDerivVectorRawPtr = &((*externData_.nextStaDerivVectorPtr)[0]);
  }

  if (externData_.nextStoVectorPtr != 0)
  {
    externData_.nextStoVectorRawPtr = &((*externData_.nextStoVectorPtr)[0]);
  }

  if (externData_.currStoVectorPtr != 0)
  {
    externData_.currStoVectorRawPtr = &((*externData_.currStoVectorPtr)[0]);
  }

  // lead current and junction voltage vectors
  if (externData_.nextLeadCurrFCompPtr != 0)
  {
    externData_.nextLeadCurrFCompRawPtr = &((*externData_.nextLeadCurrFCompPtr)[0]);
  }

  if (externData_.nextLeadCurrQCompPtr != 0)
  {
    externData_.nextLeadCurrQCompRawPtr = &((*externData_.nextLeadCurrQCompPtr)[0]);
  }

  if (externData_.nextJunctionVCompPtr != 0)
  {
    externData_.nextJunctionVCompRawPtr = &((*externData_.nextJunctionVCompPtr)[0]);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setupRawMatrixPointers_
// Purpose       : set up raw pointers for matrices
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 06/03/09
//-----------------------------------------------------------------------------
bool DeviceMgr::setupRawMatrixPointers_()
{
  InstanceVector::iterator it = instancePtrVec_.begin();
  InstanceVector::iterator end = instancePtrVec_.end();
  for ( ; it != end; ++it)
  {
    (*it)->setupPointers();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/31/01
//-----------------------------------------------------------------------------
double DeviceMgr::getMaxTimeStepSize()
{
  double maxStep = devOptions_.defaultMaxTimeStep;

  InstanceVector::iterator it = devicesWithMaxTimeStepFuncsPtrVec_.begin();
  InstanceVector::iterator end = devicesWithMaxTimeStepFuncsPtrVec_.end();
  for ( ; it != end; ++it)
  {
    double step = (*it)->getMaxTimeStepSize ();
    if (!((*it)->getFastSourceFlag()))
    {
      maxStep = std::min(step, maxStep);
    }
  }

  return maxStep;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::enablePDEContinuation
// Purpose       : This function turns on the Continuation flag, which lets
//                 the devices which are Continuation enabled know that they
//                 need to set up their variable parameters.
//
// Special Notes : Currently, only PDE devices can take advantage of this
//                 capability, and it is only used in the context of two-level
//                 Newton.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
int DeviceMgr::enablePDEContinuation()
{
  if (DEBUG_DEVICE)
  {
    dout() << "DeviceMgr::enablePDEContinuation" << std::endl;
  }

  bool bsuccess = true;
  solState_.PDEcontinuationFlag_  = true;
  int target_max_PDE_continuation_steps = 0;

  int max_PDE_continuation_steps = 1;
  while (max_PDE_continuation_steps != target_max_PDE_continuation_steps) {
    target_max_PDE_continuation_steps = max_PDE_continuation_steps;
    for (InstanceVector::iterator it = instancePtrVec_.begin(), end = instancePtrVec_.end(); it != end; ++it)
    {
      bool tmpSuccess = (*it)->enablePDEContinuation(max_PDE_continuation_steps);
      bsuccess = bsuccess && tmpSuccess;
    }
  }

  return bsuccess ? target_max_PDE_continuation_steps : -1;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::disablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
bool DeviceMgr::disablePDEContinuation()
{
  bool bsuccess = true;
  solState_.PDEcontinuationFlag_ = false;

  InstanceVector::iterator it = instancePtrVec_.begin();
  InstanceVector::iterator end = instancePtrVec_.end();
  for ( ; it != end; ++it)
  {
    bool tmpSuccess = (*it)->disablePDEContinuation();
    bsuccess = bsuccess && tmpSuccess;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::calcPDESubProblemInfo
//
// Purpose       : Determines the number of PDE sub-problems.,
//
//                 This is mainly used/needed for 2-level problems.
//
//                 Also determines the number of interface nodes per
//                 sub-problem (ie number of electrodes on each device).
//
// Special Notes : Need to modify to work correctly in parallel, probably.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
bool DeviceMgr::calcPDESubProblemInfo()
{
  // now set up numInterfaceNodes_;
  numInterfaceNodes_.reserve(pdeInstancePtrVec_.size());

  InstanceVector::const_iterator it = pdeInstancePtrVec_.begin();
  InstanceVector::const_iterator end = pdeInstancePtrVec_.end();
  for ( ; it != end; ++it)
  {
    numInterfaceNodes_.push_back((*it)->getNumExtVars());
  }

  calledBeforeCSPI = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getNumInterfaceNodes
// Purpose       : returns the vector calculaed in calcPDESubProblemInfo.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
void DeviceMgr::getNumInterfaceNodes (std::vector<int> & numINodes)
{
  if (!calledBeforeCSPI)
  {
    calcPDESubProblemInfo ();
  }

  int size = numINodes.size ();
  int size2 = numInterfaceNodes_.size();

  if (size < size2) numINodes.resize(size2);

  for (int i=0;i<size2;++i)
  {
    numINodes[i] = numInterfaceNodes_[i];
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadCouplingRHS
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
bool DeviceMgr::loadCouplingRHS (int iPDEDevice, int iElectrode, Linear::Vector * dfdvPtr)
{
  return pdeInstancePtrVec_[iPDEDevice]->loadDFDV(iElectrode,dfdvPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::calcCouplingTerms
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
bool DeviceMgr::calcCouplingTerms (int iPDEDevice, int iElectrode, const Linear::Vector * dxdvPtr)
{
  return pdeInstancePtrVec_[iPDEDevice]->calcConductance(iElectrode, dxdvPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getHomotopyBlockSize
// Purpose       : Returns the number of mosfet gainscale blocks (a value
//                 in the device options).
// Special Notes : Needed for block homotopy on gainscale.
// Scope         : public
// Creator       : Roger P. Pawlowski, SNL, Computational Sciences
// Creation Date : 01/26/2005
//-----------------------------------------------------------------------------
int DeviceMgr::getHomotopyBlockSize() const
{
  return devOptions_.numGainScaleBlocks;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/26/04
//-----------------------------------------------------------------------------
bool DeviceMgr::updateTemperature (double val)
{
  // convert to kelvin:
  double Ctemp = val;
  double Ktemp = val + CONSTCtoK;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
  {
    dout() << "In DeviceMgr::updateTemperature.  new C temp = " << Ctemp << " K temp = " << Ktemp << std::endl;
  }

  // First set the global temp.  This is used in each device if the "tempGiven"
  // variable is false.  This should be in Kelvin.
  devOptions_.temp.setVal(Ktemp);

  {
    // loop over the bsim3 models and delete the size dep params.
    ModelTypeModelVectorMap::const_iterator model_type_it = modelTypeModelVector_.find(MOSFET_B3::Traits::modelType());
    if (model_type_it != modelTypeModelVector_.end())
    {
      for (ModelVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
      {
        (*it)->clearTemperatureData ();
      }
    }
  }

  {
    // loop over the bsim4 models and delete the size dep params.
    ModelTypeModelVectorMap::const_iterator model_type_it = modelTypeModelVector_.find(MOSFET_B4::Traits::modelType());
    if (model_type_it != modelTypeModelVector_.end())
    {
      for (ModelVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
      {
        (*it)->clearTemperatureData ();
      }
    }
  }

  {
    // loop over the b3soi models and delete the size dep params.
    ModelTypeModelVectorMap::const_iterator model_type_it = modelTypeModelVector_.find(MOSFET_B3SOI::Traits::modelType());
    if (model_type_it != modelTypeModelVector_.end())
    {
      for (ModelVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
      {
        (*it)->clearTemperatureData ();
      }
    }
  }

  // Loop over all models, call processParams, with CTemp.
  // "XYCEADMS*TEMP is there to force Verilog devices, which might have
  // temperature dependence through "$temperature" instead of a "TEMP"
  // parameter, to work properly.  If so, they need the temperature set in
  // Kelvin.
  std::string tname("TEMP");
  std::string tname2("XYCEADMSMODTEMP");
  std::string tname3("XYCEADMSINSTTEMP");
  for (ModelVector::const_iterator it = modelVector_.begin(), end = modelVector_.end(); it != end; ++it)
  {
    bool success = (*it)->setParam(tname, Ctemp);
    success = (*it)->setParam(tname2, Ktemp) || success;
    success = (*it)->updateDependentParameters(Ktemp) || success;
    success = success && (*it)->processParams();
  }

  // Loop over device instances, and set the temperature.  This should be
  // in C, if going through processParams, and K if going through
  // the updateTemperature function.
  for (InstanceVector::const_iterator it = instancePtrVec_.begin(), end = instancePtrVec_.end(); it != end; ++it)
  {
    bool success = (*it)->setParam(tname, Ctemp);
    success = (*it)->setParam(tname3, Ktemp) || success;
    success = (*it)->updateDependentParameters(Ktemp) || success;
    success = success && (*it)->processParams();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::allDevicesConverged
// Purpose       : Check whether any device has taken an action that renders
//                  normal convergence checks invalid (i.e. that the current
//                  step must be assumed unconverged).
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 03/22/05
//-----------------------------------------------------------------------------
bool DeviceMgr::allDevicesConverged(
  Parallel::Machine             comm) const
{
  int allDevsConv = true;

  // if two-level, and just the inner problem, only check the PDE devices,
  // coz those are all we loaded.

  if (solState_.twoLevelNewtonCouplingMode == Nonlinear::INNER_PROBLEM)
  {
    for (InstanceVector::const_iterator it = pdeInstancePtrVec_.begin(), end = pdeInstancePtrVec_.end(); it != end; ++it)
    {
      bool tmpBool = (*it)->isConverged();
      allDevsConv = allDevsConv && tmpBool;
    }
  }
  else
  {
    for (DeviceVector::const_iterator iter = devicePtrVec_.begin(), end = devicePtrVec_.end(); iter != end; ++iter)
    {
      bool tmpBool = (*iter)->isConverged();
      allDevsConv = allDevsConv && tmpBool;
    }
  }

  Parallel::AllReduce(comm, MPI_LAND, &allDevsConv, 1);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    if (allDevsConv)
    {
      dout() << "All devices converged!" << std::endl;
    }
    else
    {
      dout() << "At least one device NOT converged!" << std::endl;
    }
  }

  return allDevsConv != 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setupExternalDevices
// Purpose       : In parallel, we need to setup all external devices
//                 and appropriately setup the list of instances
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Elec. & MicroSystems Modeling
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool DeviceMgr::setupExternalDevices(
  Parallel::Communicator &              parallel_comm)
{
  if (Parallel::is_parallel_run(parallel_comm.comm()))
  {
    InstanceVector &extern_device_vector = modelTypeInstanceVector_[ExternDevice::Traits::modelType()];

    int procID = parallel_comm.procID();
    int numProc = parallel_comm.numProc();

    int numExt = extern_device_vector.size();
    int numExtTotal = 0;
    parallel_comm.sumAll(&numExt, &numExtTotal, 1);

    InstanceVector orig_extern_device_vector = extern_device_vector;

    //resetup the instance vector to have a size totally the global
    //number of ext devices
    if (numExtTotal > 0)
    {
      extern_device_vector.resize(numExtTotal);

      int loc = 0;
      for (int proc = 0; proc < numProc; ++proc)
      {
        int cnt = 0;
        if (proc == procID)
          cnt = numExt;
        parallel_comm.bcast(&cnt, 1, proc);

        for (int i = 0; i < cnt; ++i)
        {
          if (proc == procID)
          {
            ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*extern_device_vector[loc]);
            int size = Xyce::packedByteCount(extern_device.getInstanceBlock());
            int bufSize = size + 100;
            char *buf = new char[bufSize];
            parallel_comm.bcast(&size, 1, proc);
            int pos = 0;
            Xyce::pack(extern_device.getInstanceBlock(), buf, bufSize, pos, &parallel_comm);
            parallel_comm.bcast(buf, size, proc);
            extern_device_vector[loc] = orig_extern_device_vector[i];
            delete[] buf;
          }
          else
          {
            int size = 0;
            parallel_comm.bcast(&size, 1, proc);
            int bufSize = size + 100;
            char *buf = new char[bufSize];
            parallel_comm.bcast(buf, size, proc);
            int pos = 0;
            InstanceBlock instance_block;
            Xyce::unpack(instance_block, buf, bufSize, pos, &parallel_comm);
            extern_device_vector[loc] = addExtDeviceInstance_(instance_block);
            delete[] buf;
          }
          ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*extern_device_vector[loc]);
          extern_device.setOwningProc(proc);
          extern_device.setComm(&parallel_comm);
          ++loc;
        }
      }
    }
  }
  else // serial run
  {
    ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
    if (model_type_it != modelTypeInstanceVector_.end())
    {
      for (InstanceVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
      {
        ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

        extern_device.setComm(&parallel_comm);
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateExternalDevices_
// Purpose       : Do the actual solve of the external devices
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Elec. & MicroSystems Modeling
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
void DeviceMgr::updateExternalDevices_()
{
  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end())
  {
    for (InstanceVector::const_iterator it = (*model_type_it).second.begin(); it != (*model_type_it).second.end(); ++it)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      extern_device.runExternalDevice();
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addExtDeviceInstance_
// Purpose       : adds an external device instance on the processors
//               : that don't actually own it.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/16/06
//-----------------------------------------------------------------------------
ExternDevice::Instance *
DeviceMgr::addExtDeviceInstance_(
  const InstanceBlock &         instance_block)
{
  ModelTypeId model_type;

  if (instance_block.getModelName().empty())
  {
    model_type = getModelGroup(instance_block.getInstanceName().getDeviceType());
  }
  else
  {
    model_type = modelTypeMap_[instance_block.getModelName()];
  }

  if (!model_type.defined())
  {
    Report::UserError message;
    message << "Unable to determine type of device for instance name " << instance_block.getInstanceName();
    if (!instance_block.getModelName().empty())
    {
      message << " with model name " << instance_block.getModelName();
    }
  }

  // Add an instance of this type.
  Device &device = getDeviceByModelType(model_type);
  DeviceInstance *instance = device.addInstance(instance_block, FactoryBlock(*this, devOptions_, solState_, matrixLoadData_, externData_, commandLine_));

  return static_cast<ExternDevice::Instance *>(instance);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::acceptStep
// Purpose       : Communicate to devices that the current step is accepted
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL
// Creation Date : 01/23/07
//-----------------------------------------------------------------------------
void DeviceMgr::acceptStep()
{
  // The time history for the LTRA device(s) has to be tracked
  // separately because it can be truncated if the user specifies that
  // option. This has to be called before
  // LTRAInstance::acceptStep(). Note that the DCOP is stored at
  // index zero and the first time step is stored at index 1 for the
  // time history vectors associated with the LTRA
  if (solState_.ltraDevices_)
  {
    if (solState_.dcopFlag)
    {
      solState_.ltraTimeIndex_ = 0;
      solState_.ltraTimeHistorySize_ = 10;
      solState_.ltraTimePoints_.resize(solState_.ltraTimeHistorySize_);
    }
    else
    {
      solState_.ltraTimeIndex_++;
      if (solState_.ltraTimeIndex_ >= solState_.ltraTimeHistorySize_)
      {
        solState_.ltraTimeHistorySize_ += 10;
        solState_.ltraTimePoints_.resize(solState_.ltraTimeHistorySize_);
      }
      solState_.ltraTimePoints_[solState_.ltraTimeIndex_] = solState_.currTime_;
    }
  }

  // make sure we have the right solver state.  Putting this here
  // makes it essential that acceptStep be called BEFORE any changes have
  // been made to the solution or state vectors (e.g. rotating them) and
  // while "nextTime" in the solver still references the time for which the
  // solution vector is valid.

  bool all_devices_converged = allDevicesConverged(comm_);

  bool tmpBool = setupSolverInfo(solState_, *analysisManager_, all_devices_converged, devOptions_, nlsMgrPtr_->getNonLinInfo());

  solState_.acceptedTime_ = solState_.currTime_;

  for (InstanceVector::iterator iter = instancePtrVec_.begin(); iter != instancePtrVec_.end(); ++iter)
  {
    (*iter)->acceptStep();
  }

  // If the TRYCOMPACT option is set then the LTRA model will try to
  // compact the amount of data stored and speed up the convolutions. If
  // any of the LTRA instances request this, done in their acceptStep()
  // member function, then perform the time-step compaction here.
  if (solState_.ltraDevices_)
  {
    if (solState_.ltraDoCompact_)
    {
      solState_.ltraTimePoints_[solState_.ltraTimeIndex_-1] = solState_.ltraTimePoints_[solState_.ltraTimeIndex_];

      solState_.ltraTimeIndex_--;

      // reset the flag for the next time step
      solState_.ltraDoCompact_ = false;
    }
  }


  // tell the various dependent parameters (ie ones with expressions) 
  // to advance their time steps as needed.

  Xyce::Util::Expression::clearProcessSuccessfulTimeStepMap(); // kludge but neccessary for now

  std::vector<Util::Expression> & globalExpressionsVec = globals_.global_expressions;

  // Update global params for new time and other global params
  std::vector<Util::Expression>::iterator globalExprIter = globalExpressionsVec.begin(); 
  std::vector<Util::Expression>::iterator globalExprEnd  = globalExpressionsVec.end();
  for ( ; globalExprIter != globalExprEnd; ++globalExprIter)
  {
    globalExprIter->processSuccessfulTimeStep();
  }

  EntityVector::iterator iter;
  EntityVector::iterator begin = dependentPtrVec_.begin();
  EntityVector::iterator end = dependentPtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->processSuccessfulTimeStep();
  }

}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getInitialQnorm
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/07
//-----------------------------------------------------------------------------
bool DeviceMgr::getInitialQnorm (std::vector<TimeIntg::TwoLevelError> & tleVec)
{
  bool bsuccess = true;

  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end())
  {
    int numExt = (*model_type_it).second.size();

    tleVec.resize(numExt);
    int iext = 0;
    InstanceVector::const_iterator it = (*model_type_it).second.begin();
    for ( ; it != (*model_type_it).second.end(); ++it, ++iext)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      bool bs1 = extern_device.getInitialQnorm(tleVec[iext]);
      bsuccess = bsuccess && bs1;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getInnerLoopErrorSums
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
bool DeviceMgr::getInnerLoopErrorSums(
  std::vector<TimeIntg::TwoLevelError> & tleVec) const
{
  bool bsuccess = true;

  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end())
  {
    int numExt = (*model_type_it).second.size();

    tleVec.resize(numExt);
    int iext = 0;
    InstanceVector::const_iterator it = (*model_type_it).second.begin();
    for ( ; it != (*model_type_it).second.end(); ++it, ++iext)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      bool bs1 = extern_device.getInnerLoopErrorSum(tleVec[iext]);
      bsuccess = bsuccess && bs1;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateStateArrays()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
bool DeviceMgr::updateStateArrays()
{
  bool bsuccess = true;

  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end())
  {
    InstanceVector::const_iterator it = (*model_type_it).second.begin();
    for ( ; it != (*model_type_it).second.end(); ++it)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      bool bs1 = extern_device.updateStateArrays();
      bsuccess = bsuccess && bs1;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::startTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
bool DeviceMgr::startTimeStep(
  bool                          beginIntegrationFlag,
  double                        nextTimeStep,
  double                        nextTime,
  int                           currentOrder)
{
  bool bsuccess = true;

  ModelTypeInstanceVectorMap::const_iterator model_type_it = modelTypeInstanceVector_.find(ExternDevice::Traits::modelType());
  if (model_type_it != modelTypeInstanceVector_.end())
  {
    InstanceVector::const_iterator it = (*model_type_it).second.begin();
    for ( ; it != (*model_type_it).second.end(); ++it)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      bool bs1 = extern_device.startTimeStep(
        beginIntegrationFlag,
        nextTimeStep,
        nextTime,
        currentOrder);

      bsuccess = bsuccess && bs1;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setExternalSolverState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void DeviceMgr::setExternalSolverState(
  bool  external_initJctFlag)
{
  solState_.externalStateFlag_ = true;
  solState_.externalInitJctFlag_ = external_initJctFlag;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::restartDataSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
int DeviceMgr::restartDataSize(bool pack) const
{
  int numdoubles = solState_.ltraTimePoints_.size();
  int numints = 3;
  int count = sizeof(double) * (numdoubles);
  count += sizeof(int) * numints;

  // bump up size for unpacked data.  This is an empirical multiplier.
  if (!pack)
  {
    count *= 3;
  }

  return count;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::dumpRestartData
// Purpose       : Output restart data.
// Special Notes : This function is called by the restart manager to output
//                 persistent data for the device package.  It should NOT
//                 include any data from individual devices, as that restart
//                 data is collected elsewhere.
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
bool DeviceMgr::dumpRestartData(
  char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack) const
{
  bool retval=true;

  if (pack)
  {
    int size=solState_.ltraTimePoints_.size();
    comm->pack(&(solState_.ltraTimeIndex_), 1, buf, bsize, pos);
    comm->pack(&(solState_.ltraTimeHistorySize_), 1, buf, bsize, pos);
    comm->pack(&(size), 1, buf, bsize, pos);
    comm->pack(&(solState_.ltraTimePoints_[0]), size, buf, bsize, pos);
  }
  else
  {
    int count = restartDataSize(false);
    int startIndex = pos;
    for(int i = startIndex; i < (startIndex+count); ++i) buf[i] = ' ';

    int size=solState_.ltraTimePoints_.size();
    std::ostringstream ost;
    ost.width(24);ost.precision(16);ost.setf(std::ios::scientific);
    ost << solState_.ltraTimeIndex_ << " ";
    ost << solState_.ltraTimeHistorySize_ << " ";

    if (DEBUG_RESTART)
    {
      dout() << "DeviceMgr::getRestartData:  ltraTimeIndex = " 
             << solState_.ltraTimeIndex_ << std::endl
             << "DeviceMgr::getRestartData:  ltraTimeHistorySize = " 
             << solState_.ltraTimeHistorySize_ << std::endl;
    }

    ost << size << " ";
    for (int i=0;i<size;i++)
    {
      ost << solState_.ltraTimePoints_[i] << " ";
      if (DEBUG_RESTART)
      {
        dout() << "DeviceMgr::dumpRestartData:  ltraTimePoints["<<i<<"] =" 
          << solState_.ltraTimePoints_[i]<<std::endl;
      }
    }

    std::string data(ost.str());
    for(unsigned int i = 0; i < data.length(); ++i) buf[startIndex+i] = data[i];

    // The line above copies the characters of the data string into buf,
    // but doesn't null-terminate buf.
    // it is essential to terminate the buffer with a null, or attempts
    // to construct a string object from it will get memory access problems.
    buf[startIndex+data.length()] = '\0';
    pos += data.length();
  }

  return retval;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::restoreRestartData
// Purpose       : Load restart data.
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
bool DeviceMgr::restoreRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack)
{
  bool retval=true;

  if (pack)
  {
    comm->unpack(buf, bsize, pos, &(solState_.ltraTimeIndex_), 1);
    comm->unpack(buf, bsize, pos, &(solState_.ltraTimeHistorySize_), 1);
    int size=0;
    comm->unpack(buf, bsize, pos, &(size), 1);
    solState_.ltraTimePoints_.resize(size);
    comm->unpack(buf, bsize, pos, &(solState_.ltraTimePoints_[0]), size);
  }
  else
  {
    std::string str1(buf);
    int length = str1.size() - pos;
    std::string str2(str1,pos,length);

    std::istringstream ist(str2);

    ist >> solState_.ltraTimeIndex_;
    ist >> solState_.ltraTimeHistorySize_;
    if (DEBUG_RESTART)
    {
      dout() << "DeviceMgr::restoreRestartData:  ltraTimeIndex = " 
             << solState_.ltraTimeIndex_ <<std::endl
             << "DeviceMgr::restoreRestartData:  ltraTimeHistorySize = " 
             << solState_.ltraTimeHistorySize_ <<std::endl;
    }

    int size=0;
    ist >> size;
    solState_.ltraTimePoints_.resize(size);
    for (int i=0;i<size;++i)
    {
      ist >> solState_.ltraTimePoints_[i];
      if (DEBUG_RESTART)
      {
        dout() << "DeviceMgr::restoreRestartData:  ltraTimePoints["<<i<<"] = " 
          << solState_.ltraTimePoints_[i] << std::endl;
      }
    }

    pos += ist.tellg();
  }

  return retval;
}

// ----------------------------------------------------------------------------
// Function      : findDeviceEntity
// Scope         : public
// Creator       : Dave Baur
// ----------------------------------------------------------------------------
DeviceEntity *findDeviceEntity(
    EntityTypeIdDeviceMap::const_iterator begin,
    EntityTypeIdDeviceMap::const_iterator end,
    const std::string entity_name)
{
  for (EntityTypeIdDeviceMap::const_iterator it = begin; it != end; ++it)
  {
    DeviceEntity *device_entity = (*it).second->findInstance(InstanceName(entity_name));
    if (device_entity)
      return device_entity;

    device_entity = (*it).second->findModel(ModelName(entity_name));
    if (device_entity)
      return device_entity;
  }

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getDeviceEntity
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/02/03
//-----------------------------------------------------------------------------
DeviceEntity * DeviceMgr::getDeviceEntity(
  const std::string &   full_param_name) const
{
  std::string entity_name = Xyce::Util::entityNameFromFullParamName(full_param_name).getEncodedName();
  
  DeviceEntityMap::iterator it = parameterDeviceCache_.find(full_param_name);
  if (it == parameterDeviceCache_.end())
  {
    DeviceEntity *device_entity = findDeviceEntity(deviceMap_.begin(), deviceMap_.end(), entity_name);
    parameterDeviceCache_[full_param_name] = device_entity;
    return device_entity;
  }
  else if (!(*it).second)
  {
    DeviceEntity *device_entity = findDeviceEntity(deviceMap_.begin(), deviceMap_.end(), entity_name);
    if (device_entity!=0)
    {
      (*it).second = device_entity;
      return device_entity;
    }
    else
    {
      // The less-common case of a device in a subcircuit, such as the R device, that has
      // a default instance parameter.  So, this handles X1:R1 (instead of the more verbose X1:R1:R)
      std::string entity_name2 = Xyce::Util::entityNameFromDefaultParamName(full_param_name).getEncodedName();
      DeviceEntity *device_entity2 = findDeviceEntity(deviceMap_.begin(), deviceMap_.end(), entity_name2);
      (*it).second = device_entity2;
      return device_entity2;
    }
  }
  else
    return (*it).second;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addDevicesToCount
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 8/06/14
//-----------------------------------------------------------------------------
void DeviceMgr::addDevicesToCount(const DeviceCountMap & device_map)
{
  // Add in devices to map for counting.
  for (DeviceCountMap::const_iterator it = device_map.begin(), end = device_map.end(); it != end; ++it)
  {
    if ( localDeviceCountMap_[(*it).first] )
    {
       localDeviceCountMap_[(*it).first] += (*it).second;
    }
    else
    {
       localDeviceCountMap_[(*it).first] = (*it).second;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Mems Modeling
// Creation Date : 10/20/2008
//-----------------------------------------------------------------------------
bool registerPkgOptionsMgr(DeviceMgr &device_manager, IO::PkgOptionsMgr &options_manager)
{
  DeviceOptions::populateMetadata(options_manager);

  options_manager.addCommandProcessor("SENS", IO::createRegistrationOptions(device_manager, &DeviceMgr::registerSensParams));
  options_manager.addCommandProcessor("OP", IO::createRegistrationOptions(device_manager, &DeviceMgr::setOPAnalysisParams));
  options_manager.addCommandProcessor("HB", IO::createRegistrationOptions(device_manager, &DeviceMgr::setHBAnalysisParams));
  options_manager.addCommandProcessor("AC", IO::createRegistrationOptions(device_manager, &DeviceMgr::setACAnalysisParams));
  options_manager.addCommandProcessor("NOISE", IO::createRegistrationOptions(device_manager, &DeviceMgr::setNOISEAnalysisParams));
  options_manager.addOptionsProcessor("DEVICE", IO::createRegistrationOptions(device_manager, &DeviceMgr::setDeviceOptions));
  options_manager.addOptionsProcessor("SENSITIVITY", IO::createRegistrationOptions(device_manager, &DeviceMgr::setSensitivityOptions));

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getOp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Baur
// Creation Date : 
//-----------------------------------------------------------------------------
Util::Op::Operator * DeviceMgr::getOp(
  Parallel::Machine     comm,
  const std::string &   name) const
{
  OpMap::const_iterator it = opMap_.find(name);
  if (it != opMap_.end())
  {
    return (*it).second;
  }
  else 
  {
    Util::ParamList param_list;
    param_list.push_back(Param(name, ""));
    Util::ParamList::const_iterator it = param_list.begin();

    return opMap_[name] = makeOp(comm, opBuilderManager_, it);
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getParamAndReduce
// Purpose       : Returns the current value of a named parameter.
//
// Special Notes : This works in parallel.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/26/03
//-----------------------------------------------------------------------------
bool getParamAndReduce(
  Parallel::Machine     comm,
  const DeviceMgr &     device_manager,
  const std::string &   name,
  double &              value)
{
  Util::Op::Operator *op = device_manager.getOp(comm, name);
  value = op ? (*op)(comm, Util::Op::OpData()).real() : 0.0;

  return op;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getParamAndReduce
// Purpose       : Returns the current value of a named parameter.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/26/03
//-----------------------------------------------------------------------------
double getParamAndReduce(
  Parallel::Machine     comm,
  const DeviceMgr &     device_manager,
  const std::string &   name)
{
  double value = 0.0;
  bool found = getParamAndReduce(comm, device_manager, name, value);

  if (!found)
  {
    if (DEBUG_DEVICE)
    {
      Report::DevelWarning() <<
        "DeviceMgr::getParamAndReduce.  Unable to find parameter " << name;
    }
    else
    {
      Report::UserError() << "Unable to find parameter " << name;
    }
  }

  return value;
}

namespace {

//-----------------------------------------------------------------------------
// Function      : setParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool setParameter(
  Parallel::Machine             comm,
  ArtificialParameterMap &      artificial_parameter_map,
  PassthroughParameterSet &     passthrough_parameter_map,
  Globals &                     globals,               ///< global variables
  DeviceMgr &                   device_manager,
  EntityVector &                dependent_entity_vector, // dependentPtrVec_, if called from DeviceMgr::setParam
  const InstanceVector &        extern_device_vector,
  const std::string &           name,
  double                        value,
  bool                          override_original)
{
  bool bsuccess = true, success = true;
  GlobalParameterMap &  global_parameter_map = globals.global_params;

  ArtificialParameterMap::iterator artificial_param_it = artificial_parameter_map.find(name);
  if (artificial_param_it != artificial_parameter_map.end())
  {
    (*artificial_param_it).second->setValue(device_manager, value);
  }
  else
  {
    GlobalParameterMap::iterator global_param_it = global_parameter_map.find(name);
    if (global_param_it != global_parameter_map.end())
    {
      if ((*global_param_it).second != value)
      {
        (*global_param_it).second = value;
        EntityVector::iterator it = dependent_entity_vector.begin();
        EntityVector::iterator end = dependent_entity_vector.end();
        for ( ; it != end; ++it)
        {
          bool globalParamChangedLocal=true;
          bool timeChangedLocal=false;
          bool freqChangedLocal=false;
          if ((*it)->updateGlobalAndDependentParameters(globalParamChangedLocal,timeChangedLocal,freqChangedLocal))
          {
            (*it)->processParams();
            (*it)->processInstanceParams();
          }
        }
      }

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
        Xyce::dout() << " in DeviceMgr setParam, setting global parameter "<< name << " to " << value << std::endl;
        if (override_original)
        {
          Xyce::dout()  << " overriding original";
        }
        Xyce::dout() << std::endl;
      }

    }
    else
    {
      // If not artificial, then search for the appropriate natural param(s).
      DeviceEntity * device_entity = device_manager.getDeviceEntity(name);

      int entity_found = (device_entity != 0);
      if (entity_found)
      {
        bool found;
        std::string paramName = Util::paramNameFromFullParamName(name);
        if (paramName == "")
        {
          //if (DEBUG_DEVICE)
          if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
          {
            Xyce::dout() << " in DeviceMgr setParam, setting default parameter to " << value ;
            if (override_original)
            {
              Xyce::dout()  << " overriding original";
            }
            Xyce::dout() << std::endl;
          }
          found = device_entity->setDefaultParam(value, override_original);
        }
        else
        {
          //if (DEBUG_DEVICE)
          if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
          {
            Xyce::dout() << " in DeviceMgr setParam, setting parameter "<< paramName 
              << " to " << value << std::endl;
          }
          found = device_entity->setParam(paramName, value);
        }

        if (found)
        {
          device_entity->processParams (); // if this "entity" is a model, then need to
          // also do a  "processParams" on the related instances.
          device_entity->processInstanceParams();
        }

        entity_found = found;
      }

      Parallel::AllReduce(comm, MPI_LOR, &entity_found, 1);

      if (entity_found == 0)
      {
        //if (DEBUG_DEVICE)
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
        {
          Report::DevelWarning() << "DeviceMgr::setParam.  Unable to find parameter " << name;
        }
        else
        {
          Report::UserError() << "Unable to find parameter " << name;
        }
      }
    }
  }

  // Certain parameters should be passed through to the inner solve,
  // if there is an inner circuit problem.  The names of parameters
  // that should be passed through are stored in the map.
  if (passthrough_parameter_map.find(name) != passthrough_parameter_map.end())
  {
    InstanceVector::const_iterator it = extern_device_vector.begin();
    InstanceVector::const_iterator end = extern_device_vector.end();
    for ( ; it != end; ++it)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      bool bs1 = extern_device.setInternalParam(name, value);
    }
  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : setParameterRandomExpressionTerms
// Purpose       : Update params that depend on random ops such as
// Special Notes : AGAUSS, GAUSS, AUNIF, UNIF, RAND and LIMIT
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 7/30/2020
//-----------------------------------------------------------------------------
bool setParameterRandomExpressionTerms(
  Parallel::Machine             comm,
  ArtificialParameterMap &      artificial_parameter_map,
  PassthroughParameterSet &     passthrough_parameter_map,
  Globals &                     globals,               ///< global variables
  DeviceMgr &                   device_manager,
  EntityVector &                dependent_entity_vector,   // dependentPtrVec_, if called from DeviceMgr::setParam
  const InstanceVector &        extern_device_vector,
  const std::string &           name,
  const std::string &           opName,
  int opIndex,
  //enum Util::astRandTypes astType,
  int astType,
  double                        value,
  bool                          override_original)
{
  bool bsuccess = true, success = true;

  GlobalParameterMap &          global_parameter_map = globals.global_params;
  std::vector<Util::Expression> & global_expressions = globals.global_expressions;
  std::vector<std::string> & global_exp_names = globals.global_exp_names;

#if 0
  // artificial params probabaly can't work for AGAUSS, etc.
  ArtificialParameterMap::iterator artificial_param_it = artificial_parameter_map.find(name);
  if (artificial_param_it != artificial_parameter_map.end())
  {
    (*artificial_param_it).second->setValue(device_manager, value);
  }
  else
#endif
  {
    GlobalParameterMap::iterator global_param_it = global_parameter_map.find(name);

    if (global_param_it != global_parameter_map.end())
    {
      // find the expression:
      std::vector<std::string>::iterator name_it = std::find(global_exp_names.begin(), global_exp_names.end(), name);
      int globalIndex = std::distance (global_exp_names.begin(), name_it );

      Util::Expression &expression = global_expressions[globalIndex];

      // set value .  The variable "value" was determined in one of the UQ classes.  
      // It represents the value of a particular operator within the expression.  
      // Use the opIndex to identify exactly which operator receives this value.
      switch(astType)
      {
        case Util::AST_AGAUSS:
          expression.setAgaussValue(opIndex,value);
          break;
        case Util::AST_GAUSS:
          expression.setGaussValue(opIndex,value);
          break;
        case Util::AST_AUNIF:
          expression.setAunifValue(opIndex,value);
          break;
        case Util::AST_UNIF:
          expression.setUnifValue(opIndex,value);
          break;
        case Util::AST_RAND:
          expression.setRandValue(opIndex,value);
          break;
        case Util::AST_LIMIT:
          expression.setLimitValue(opIndex,value);
          break;
      }
      double paramValue=0.0;
      expression.evaluateFunction(paramValue);

      if ((*global_param_it).second != value)
      {
        (*global_param_it).second = value;
        EntityVector::iterator it = dependent_entity_vector.begin();
        EntityVector::iterator end = dependent_entity_vector.end();
        for ( ; it != end; ++it)
        {
          bool globalParamChangedLocal=true;
          bool timeChangedLocal=false;
          bool freqChangedLocal=false;
          if ((*it)->updateGlobalAndDependentParameters(globalParamChangedLocal,timeChangedLocal,freqChangedLocal))
          {
            (*it)->processParams();
            (*it)->processInstanceParams();
          }
        }
      }

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
        Xyce::dout() << " in DeviceMgr setParam, setting global parameter "<< name << " to " << value << std::endl;
        if (override_original)
        {
          Xyce::dout()  << " overriding original";
        }
        Xyce::dout() << std::endl;
      }

    }
    else
    {

      // ERK.  This needs to be update/fixed.  Using a setDefaultParam or setParam, below is incorrect.
      //
      // The correct thing to do is pull out the expression, and then do a setAgauss, setAunif, or whichever
      // is appropriate.

      // If not artificial, then search for the appropriate natural param(s).
      DeviceEntity * device_entity = device_manager.getDeviceEntity(name);

      int entity_found = (device_entity != 0);
      if (entity_found)
      {
        bool found;
        std::string paramName = Util::paramNameFromFullParamName(name);
        if (paramName == "")
        {
          //if (DEBUG_DEVICE)
          if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
          {
            Xyce::dout() << " in DeviceMgr setParam, setting default parameter to " << value ;
            if (override_original)
            {
              Xyce::dout()  << " overriding original";
            }
            Xyce::dout() << std::endl;
          }
          found = device_entity->setDefaultParam(value, override_original);
        }
        else
        {
          //if (DEBUG_DEVICE)
          if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
          {
            Xyce::dout() << " in DeviceMgr setParam, setting parameter "<< paramName 
              << " to " << value << std::endl;
          }
          found = device_entity->setParam(paramName, value);
        }

        if (found)
        {
          device_entity->processParams (); // if this "entity" is a model, then need to
          // also do a  "processParams" on the related instances.
          device_entity->processInstanceParams();
        }

        entity_found = found;
      }

      Parallel::AllReduce(comm, MPI_LOR, &entity_found, 1);

      if (entity_found == 0)
      {
        //if (DEBUG_DEVICE)
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
        {
          Report::DevelWarning() << "DeviceMgr.C:  setParameterRandomExpressionTerms.  Unable to find parameter " << name;
        }
        else
        {
          Report::UserError() << "Unable to find parameter " << name;
        }
      }
    }
  }

#if 0
  // ERK. Check this later.

  // Certain parameters should be passed through to the inner solve,
  // if there is an inner circuit problem.  The names of parameters
  // that should be passed through are stored in the map.
  if (passthrough_parameter_map.find(name) != passthrough_parameter_map.end())
  {
    InstanceVector::const_iterator it = extern_device_vector.begin();
    InstanceVector::const_iterator end = extern_device_vector.end();
    for ( ; it != end; ++it)
    {
      ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

      bool bs1 = extern_device.setInternalParam(name, value);
    }
  }
#endif

  return bsuccess;
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addGlobalParameter()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/18/05
//-----------------------------------------------------------------------------
void addGlobalParameter(
  SolverState &         solver_state,
  double                temp,
  Globals &             globals,
  const Util::Param &   param)
{
  if (param.getType() == Util::EXPR)
  {
    globals.global_expressions.push_back(param.getValue<Util::Expression>());
    globals.global_exp_names.push_back(param.uTag());

    Util::Expression &expression = globals.global_expressions.back();
    double val;
    expression.evaluateFunction(val);
    expression.clearOldResult();
    globals.global_params[param.uTag()] = val;
  }
  else
  {
    globals.global_params[param.uTag()] = param.getImmutableValue<double>();
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::findGlobalPar
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
const double * findGlobalParameter(
  const GlobalParameterMap &    global_map,
  const std::string &           name)
{
  GlobalParameterMap::const_iterator it = global_map.find(name);
  if (it != global_map.end())
    return &(*it).second;

  return 0;
}

} // namespace Device
} // namespace Xyce

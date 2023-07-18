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
#include <numeric>
#include <sstream>
#include <stdexcept>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_SweepParam.h>
#include <N_DEV_Algorithm.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_ExternDevice.h>
#include <N_DEV_MutIndLin.h>
#include <N_DEV_MutIndNonLin.h>
#include <N_DEV_MutIndNonLin2.h>
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
#include <N_TIA_WorkingIntegrationMethod.h>
#include <N_TOP_Topology.h>
#include <N_UTL_Algorithm.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Op.h>
#include <N_UTL_OpBuilder.h>
#include <N_UTL_Expression.h>
#include <N_UTL_HspiceBools.h>

#include <Teuchos_RCP.hpp>
#include <expressionGroup.h>
#include <N_DEV_ExpressionGroupWrapper.h>

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
  UserDefinedParams &           globals,
  DeviceMgr &                   device_manager,
  EntityVector &                dependent_entity_vector,
  const InstanceVector &        extern_device_vector,
  const std::string &           name,
  double                        value,
  bool                          override_original);

//-----------------------------------------------------------------------------
// Function to set: AGAUSS, GAUSS, AUNIF, UNIF, RAND and LIMIT inside of 
// expressions/params/etc.  This version of the function sets all the parameters
// in a single function call, and is much faster than the original version of
// this function (above).
//-----------------------------------------------------------------------------
bool setParameterRandomExpressionTerms2(
  Parallel::Machine             comm,
  ArtificialParameterMap &      artificial_parameter_map,
  PassthroughParameterSet &     passthrough_parameter_map,
  UserDefinedParams &           globals,
  DeviceMgr &                   device_manager,
  EntityVector &                dependent_entity_vector,
  const InstanceVector &        extern_device_vector,
  const std::vector<Xyce::Analysis::SweepParam> & SamplingParams,
  std::unordered_map<DeviceEntity *, std::vector<Depend> > & masterGlobalParamDependencies,
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
    subcktGlobals_(solState_.subcktGlobals_),
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
    SAMPLINGSpecified_(false),
    ESSpecified_(false),
    PCESpecified_(false),
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
  solState_.getGroupWrapper()->expressionGroup_ = group;
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
// Function      : DeviceMgr::setParserOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/29/2020
//-----------------------------------------------------------------------------
bool DeviceMgr::setParserOptions (const Util::OptionBlock & OB)
{
  return devOptions_.setParserOptions (OB);
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
// Function      : DeviceMgr::setDiagnosticOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool DeviceMgr::setDiagnosticOptions (const Util::OptionBlock & OB)
{
  //return devOptions_.setOptions(OB);
  bool retval = false;
  for (Util::ParamList::const_iterator it = OB.begin(), end = OB.end(); it != end; ++it)
  {
    const Util::Param &param = *it;
    if( param.tag() == "CURRENTLIMIT")
    {
      devOptions_.calculateAllLeadCurrents = true;
      retval=true;
    }
  }
  return retval;
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
// Function      : DeviceMgr::setSamplingParams 
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/11/2021
//-----------------------------------------------------------------------------
bool DeviceMgr::setSamplingParams (const Util::OptionBlock & option_block)
{
  SAMPLINGSpecified_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setEmbeddedSamplingParams 
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/11/2021
//-----------------------------------------------------------------------------
bool DeviceMgr::setEmbeddedSamplingParams (const Util::OptionBlock & option_block)
{
  ESSpecified_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setPCEParams 
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/11/2021
//-----------------------------------------------------------------------------
bool DeviceMgr::setPCEParams (const Util::OptionBlock & option_block)
{
  PCESpecified_ = true;
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
// Function      : DeviceMgr::setEarlyNoiseFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void DeviceMgr::setEarlyNoiseFlag(std::string analysisName)
{
  if (analysisName == "NOISE")
    solState_.earlyNoiseFlag_ = true;

}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setSPAnalysisFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/20/2021
//-----------------------------------------------------------------------------
void DeviceMgr::setDisableInitJctFlags(bool flag)
{
  devOptions_.disableInitJctFlag = flag;
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
  externData_.lastStoVectorPtr = 0;
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

  InstanceVector::iterator iter = instancePtrVec_.begin();
  InstanceVector::iterator end = instancePtrVec_.end();
  for (; iter!= end; iter++)
  {
    (*iter)->setupBreakPoints();
  }
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::notify
// Purpose       : Receive "AnalysisEvent" notifications from publishers of such.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/10/2023
//-----------------------------------------------------------------------------
void DeviceMgr::notify(const Analysis::AnalysisEvent &analysis_event)
{
  if (analysis_event.state_ == Analysis::AnalysisEvent::STEP_STARTED)
  {
    GlobalParameterMap & globalParamMap = globals_.paramMap;
    std::vector<Util::Expression> & globalExpressionsVec = globals_.expressionVec;

    bool globalParamChanged = false, timeChanged = false, freqChanged = false;
    std::vector<Util::Expression>::iterator globalExprIter = globalExpressionsVec.begin(); 
    std::vector<Util::Expression>::iterator globalExprEnd  = globalExpressionsVec.end();

    if ( analysis_event.outputType_ == Analysis::AnalysisEvent::AC ||
         analysis_event.outputType_ == Analysis::AnalysisEvent::NOISE)
    {
      // this section is mainly called when frequency is changed.  However, AC and NOISE can update other things
      // on their sweep, including device and/or global parameters.  So, we can't strictly update only frequency 
      // dependent things.
      // Note also, that it isn't necessary to explicitly update the global parameters to support device parameters.
      // Device parameter updates will automatically account for them anyway, since they share the same AST.
      // The only reason to explicitly update global params here is to support printing them out on the .PRINT line, 
      // when they are requested outside an expression.  (ie. ".print tran param" rather than ".print tran {param}"; these two use cases follow different paths thru the code).
      // If they are specified inside an expression, then the AST is set up automatically.
      // So, if the .PRINT line can be refactored a bit, this loop over global params can go away.
      freqChanged = true;
      bool all_devices_converged = true; // this should be true at this stage
      setupSolverInfo(solState_, *analysisManager_, all_devices_converged, devOptions_, nlsMgrPtr_->getNonLinInfo());

      // Update global params as needed.  This is mainly to support updated frequency.
      if ( !(globalExpressionsVec.empty()) )
      {
        int pos = 0;
        for ( ; globalExprIter != globalExprEnd; ++globalExprIter)
        {
          //if (globalExprIter->isFreqDependent())
          {
            double val;
            if (globalExprIter->evaluateFunction(val))
            {
              globalParamChanged = true;
              globalParamMap[globals_.expNameVec[pos]] = val;
            }
          }
          ++pos;
        }
      }

      if (firstDependent_)
      {
        firstDependent_ = false;
        setupDependentEntities();
      }

      if ( ! (dependentPtrVec_.empty()) ) // update all device parameters, not just freq-dependent ones
      {
        EntityVector::iterator iter;
        EntityVector::iterator begin = dependentPtrVec_.begin();
        EntityVector::iterator end = dependentPtrVec_.end();
        for (iter=begin; iter!=end;++iter)
        {
          bool changed = (*iter)->updateGlobalAndDependentParametersForStep(true,true,true);
          (*iter)->processParams();
          (*iter)->processInstanceParams();
        }
      }
    }
    else if ( analysis_event.outputType_ == Analysis::AnalysisEvent::TRAN )
    {
      timeChanged = true;
      bool all_devices_converged = true; // this should be true at this stage
      setupSolverInfo(solState_, *analysisManager_, all_devices_converged, devOptions_, nlsMgrPtr_->getNonLinInfo());

      // Update global params as needed.  This is mainly to support updated time, so only check time-dependent params
      // Note also, that it isn't necessary to explicitly update the global parameters to support device parameters.
      // Device parameter updates will automatically account for them anyway, since they share the same AST.
      // The only reason to explicitly update global params here is to support printing them out on the .PRINT line, 
      // when they are requested outside an expression.  (ie. ".print tran param" rather than ".print tran {param}"; these two use cases follow different paths thru the code).
      // If they are specified inside an expression, then the AST is set up automatically.
      // So, if the .PRINT line can be refactored a bit, this loop over global params can go away.
      if ( !(globalExpressionsVec.empty()) )
        { 
        int pos = 0;
        for ( ; globalExprIter != globalExprEnd; ++globalExprIter)
        {
          if (globalExprIter->isTimeDependent())
          {
            double val;
            if (globalExprIter->evaluateFunction(val))
            {
              globalParamChanged = true;
              globalParamMap[globals_.expNameVec[pos]] = val;
            }
          }
          ++pos;
        }
      }

      if (firstDependent_)
      {
        firstDependent_ = false;
        setupDependentEntities();
      }

      if ( ! (timeDepEntityPtrVec_.empty()) )
      {
        EntityVector::iterator iter;
        EntityVector::iterator begin = timeDepEntityPtrVec_.begin();
        EntityVector::iterator end = timeDepEntityPtrVec_.end();

        for (iter=begin; iter!=end;++iter)
        {
          bool changed = (*iter)->updateTimeDependentParameters();
          if (changed)
          {
            (*iter)->processParams();
            (*iter)->processInstanceParams();
          }
        }
      }
    }

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
  if (devOptions_.checkForZeroResistance && ( model_type == Resistor::Traits::modelType()))
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
          bool isVariablesDependent=false; // this means global parameters
          bool isSpecialsDependent=false;
          bool isRandomDependent=false;
          bool isVoltageNodeDependent=false;
          bool isDeviceCurrentDependent=false;

          // check if this is a time-dependent, or variable-dependent expression.
          // If it is, then skip.
          const Param * devPar = &(*(currentParam));
          const Util::Param * tmpPar = (dynamic_cast<const Util::Param*> (devPar));
          // only check if this is an expression-type parameter.
          if (tmpPar->getType() == Util::EXPR)
          {
            Util::Expression tmpExp = tmpPar->getValue<Util::Expression>();
            isVariablesDependent=tmpExp.getVariableDependent(); 
            isSpecialsDependent=tmpExp.getSpecialsDependent();      
            isRandomDependent = tmpExp.isRandomDependent();
            isVoltageNodeDependent = tmpExp.getVoltageNodeDependent();
            isDeviceCurrentDependent = tmpExp.getDeviceCurrentDependent();
          }

          if (!isSpecialsDependent && !isVariablesDependent && !isRandomDependent &&
            !isVoltageNodeDependent  && !isDeviceCurrentDependent )
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
  DeviceInstance *instance=0;
  if (model_type == Resistor3::Traits::modelType() && !(instance_block.getModelName().empty()) )
  {
    // This is designed to handle a zero-valued resistor that refers to a .model.
    // To avoid "can't find model" errors, the model name is deleted.
    // The problem is that if a model was logged as a Resistor::Traits::modelType, then
    // Resistor3::Traits::modelType device won't be able to find it.  And, it doesn't
    // use it anyway.  Deleting it here is a kludgey way to solve the problem.
    InstanceBlock instance_block_copy = instance_block;
    instance_block_copy.setModelName(std::string(""));
    instance = device.addInstance(
        instance_block_copy,
        FactoryBlock(*this, devOptions_, solState_, matrixLoadData_, externData_, commandLine_));
  }
  else
  {
    instance = device.addInstance(
        instance_block,
        FactoryBlock(*this, devOptions_, solState_, matrixLoadData_, externData_, commandLine_));
  }

  // if the addInstance function fails, it returns instance=0.   
  // This will (correctly) trigger an error trap upstream from this function
  if(instance != 0) 
  {
    localDeviceCountMap_[device.getDefaultModelName()]++;

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
// Function      : DeviceMgr::getSourceDeviceNamesDCVal
// Purpose       : Return the fully expanded names for the source devices in
//                 independentSourceMap_ and the DC value.  These will be used 
//                 to perform source stepping for the DCOP calculation.
// Special Notes : This list is distributed in parallel and needs to be broadcast
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 5/2/23
//-----------------------------------------------------------------------------
std::map<std::string, std::pair<double,int> > DeviceMgr::getSourceDeviceNamesDCVal(Parallel::Machine comm) const
{
  std::map<std::string, std::pair<double,int> > globalSources, localSources;
  IndependentSourceMap::const_iterator it = independentSourceMap_.begin();
  IndependentSourceMap::const_iterator it_end = independentSourceMap_.end();
  for ( ; it != it_end; ++it )
  {
    // Determine how many devices are adjacent to this one.
    NodeID sourceNode( it->first, _DNODE );
    std::vector<NodeID> adjNodes;
    topology_.returnAdjIDs( sourceNode, adjNodes );
    int numAdjDevices = 0;
    for (std::vector<NodeID>::iterator adj_it = adjNodes.begin(); adj_it!= adjNodes.end(); ++adj_it)
    {
      std::vector<NodeID> adjDevs;
      topology_.returnAdjIDs( *adj_it, adjDevs);
      numAdjDevices += adjDevs.size();
    }
    it->second->processParams();
    it->second->updateSource();
    localSources[ it->first ] = std::make_pair( it->second->getDefaultParam(), numAdjDevices );
  }

  // Collect the global list of source devices
  if (Parallel::is_parallel_run(comm))
  {
    int numSources = localSources.size();
    std::vector<int> procSources( Parallel::size(comm) );
    procSources[Parallel::rank(comm)] = numSources;
    Parallel::AllReduce(comm, MPI_SUM, &procSources[0], Parallel::size(comm));

    std::map<std::string, std::pair<double,int> >::iterator map_it = localSources.begin();
    std::map<std::string, std::pair<double,int> >::iterator map_end = localSources.end();

    for (int proc=0; proc < Parallel::size(comm); ++proc)
    {
      for (int pNumSources=0; pNumSources<procSources[proc]; ++pNumSources)
      {
        int numAdj = 0;
        double dcVal = 0.0;
        std::string sourceName;
        if ((proc == Parallel::rank(comm)) && (map_it != map_end))
        {
          // Get the source name and value, then increment the iterator
          sourceName = map_it->first;
          dcVal = (map_it->second).first;
          numAdj = (map_it->second).second;
          map_it++;
        }
        // Broadcast the source name and value
        Parallel::Broadcast( comm, sourceName, proc ); 
        Parallel::Broadcast( comm, &dcVal, 1, proc );
        Parallel::Broadcast( comm, &numAdj, 1, proc );

        // Add it to the global map
        globalSources[ sourceName ] = std::make_pair( dcVal, numAdj ); 
      }
    }
  }
  else
  {
    return localSources;
  }
 
  // Return the sources on all processors, if run in parallel
  return globalSources;  
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
// Function      : DeviceMgr::finalizeLeadCurrentRequests
// Purpose       : Do the actual calls to enable lead currents in devices
//                 that need them.
// Special Notes : This must be called after devices have been constructed,
//                 but before the matrix is set up, when the lead current
//                 data structures are actually allocated.
//                 This function used to be performed in addDeviceInstance,
//                 but that is too early to capture any lead current
//                 requests that have been made by external outputters (all
//                 of which happen after initial netlist parsing and
//                 construction of devices)
//
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 7/22/2021
//-----------------------------------------------------------------------------
void DeviceMgr::finalizeLeadCurrentRequests()
{
  for (InstanceVector::iterator it = instancePtrVec_.begin(),
         end = instancePtrVec_.end (); it != end; ++it)
  {
    DeviceInstance *instance=*it;

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
        if (devicesNeedingLeadCurrentLoads_.find(inductorNames[i])
            != devicesNeedingLeadCurrentLoads_.end())
        {
          mutualInductorNeedsLeadCurrents = true;
          break;  // break after finding the first match
        }
      }
    }

    if ( mutualInductorNeedsLeadCurrents
         || devOptions_.calculateAllLeadCurrents
         || (devicesNeedingLeadCurrentLoads_.find(outputName)
             != devicesNeedingLeadCurrentLoads_.end())
         || iStarRequested_ )
    {
      instance->enableLeadCurrentCalc();

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
        dout() << "DeviceMgr::finalizeLeadCurrentRequests: Enabling lead current load for device \""
          << instance->getName()
          << "\" ->  \""
          << outputName
          << "\"" << std::endl;
      }
    }
    else if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      dout() << "DeviceMgr::finalizeLeadCurrentRequests: Cannot enable lead current load for device \""
             << instance->getName()
             << "\" ->  \""
             << outputName
             << "\""
             << std::endl;
    }

    // We can't determine this until after lead currents finalized
    isLinearSystem_ = isLinearSystem_ && instance->isLinearDevice();

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

  if (!available)
  {
    // handle the special case of inductor instances contained inside mutual inductors.
    int inductorIndex=-1;
    DeviceInstance *device_instance=getMutualInductorDeviceInstance(name,inductorIndex);
    available = (device_instance != 0);
    Parallel::AllReduce(comm_, MPI_LOR, &available, 1);
  }

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
      found = device_entity->getNumericalSensitivity(name,
                                                    dfdpVec, dqdpVec, dbdpVec,
                                                    FindicesVec, QindicesVec, BindicesVec);
    }
  }
  else
  {
    // handle the special case of inductor instances contained inside mutual inductors.
    int inductorIndex=-1;
    DeviceInstance *device_instance=getMutualInductorDeviceInstance(name,inductorIndex);

    if (device_instance)
    {
      found = device_instance->getNumericalSensitivity(name,
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
// Function      : DeviceMgr::setupDependentEntities
//
// Purpose       : This function sets up vectors of pointers to entities which
//                 contain dependent params
//
// Special Notes : The "dependentPtrVec_" contains entities which depend 
//                 on any kind of active expression.  So, they can depend on 
//                 time, frequency, global parameters and/or solution variables.
//
//                 The "solnDepEntityPtrVec_" contains ONLY entities with 
//                 solution dependent parameters.  This should be a much smaller 
//                 vector than dependentPtrVec_.  There are only a few devices 
//                 that allow solution dependent params, including the Bsrc,
//                 capacitor, resistor and mutual inductor.
//
//                 Solution-dependent parameters need to be updated at every
//                 Newton step.  Other dependencies, like time, freq or global 
//                 params can be updated much less frequently.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/13/2023
//-----------------------------------------------------------------------------
void DeviceMgr::setupDependentEntities ()
{
  dependentPtrVec_.clear();
  solnDepEntityPtrVec_.clear();
  timeDepEntityPtrVec_.clear();
  freqDepEntityPtrVec_.clear();

  // do the models:
  ModelVector::iterator iterM;
  ModelVector::iterator beginM =modelVector_.begin();
  ModelVector::iterator endM =modelVector_.end();
  for (iterM=beginM; iterM!=endM;++iterM)
  {
    const std::vector<Depend> & depParams = (*iterM)->getDependentParams();
    if (!(depParams.empty()))
    {
      dependentPtrVec_.push_back(static_cast<DeviceEntity *>(*iterM));
    }

    // find solution, time, freq dependent params
    // Note, as of this writing, I don't think any .MODELs are allowed 
    // to have solution-dependent parameters (that may change in the future).  
    // So this model loop is usually a no-op.
    bool solnDep=false, timeDep=false, freqDep=false; 
    for (int ii=0;ii<depParams.size();ii++) 
    {
      if (depParams[ii].expr->isSolutionDependent()) { solnDep=true; }
      if (depParams[ii].expr->isTimeDependent() && !(  depParams[ii].expr->isSolutionDependent())) { timeDep=true; }
      if (depParams[ii].expr->isFreqDependent()  && !(  depParams[ii].expr->isSolutionDependent())) { freqDep=true; }
    }
    if (solnDep) { solnDepEntityPtrVec_.push_back(static_cast<DeviceEntity *>(*iterM)); }
    if (timeDep) { timeDepEntityPtrVec_.push_back(static_cast<DeviceEntity *>(*iterM)); }
    if (freqDep) { freqDepEntityPtrVec_.push_back(static_cast<DeviceEntity *>(*iterM)); }
  }

  // do the instances
  InstanceVector::iterator iter;
  InstanceVector::iterator begin =instancePtrVec_.begin();
  InstanceVector::iterator end =instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    const std::vector<Depend> & depParams =
      (*iter)->getDependentParams();

    const std::unordered_map <std::string, int> & dependentParamExcludeMap =
      (*iter)->getDependentParamExcludeMap();

    if (!(depParams.empty()))
    {
      dependentPtrVec_.push_back(static_cast<DeviceEntity *>(*iter));
    }

    // find solution, time, freq dependent params
    bool solnDep=false, timeDep=false, freqDep=false; 
    for (int ii=0;ii<depParams.size();ii++)
    {
      // Only set the solnDep flag if :
      // (1) the expression is solution dependent and
      // (2) this parameter has not been tagged for exclusion from this list.
      // Most solution-dependent parameters (but not all) are just handled directly in
      // the device, so they don't need to be updated from the deviceMgr. Hence,
      // they are kept out of this list.
      if (depParams[ii].expr->isSolutionDependent())
      {
        if ( !(dependentParamExcludeMap.empty()) )
        {
          if ( dependentParamExcludeMap.find( depParams[ii].name ) == dependentParamExcludeMap.end() )
          {
            solnDep=true;
          }
        }
        else
        {
          solnDep=true;
        }
      }
      if (depParams[ii].expr->isTimeDependent() && !(  depParams[ii].expr->isSolutionDependent())) { timeDep=true; }
      if (depParams[ii].expr->isFreqDependent()  && !(  depParams[ii].expr->isSolutionDependent())) { freqDep=true; }
    }
    if (solnDep) { solnDepEntityPtrVec_.push_back(static_cast<DeviceEntity *>(*iter)); }
    if (timeDep) { timeDepEntityPtrVec_.push_back(static_cast<DeviceEntity *>(*iter)); }
    if (freqDep) { freqDepEntityPtrVec_.push_back(static_cast<DeviceEntity *>(*iter)); }
  }

  return;
}

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
  if (firstDependent_)
  {
    firstDependent_ = false;
    setupDependentEntities();
  }

  // The call to updateTimeInfo is necessary for expresssion breaktpoints.
  // It cannot be called from the DeviceMgr::notify
  // function, as the notify functions happen in the wrong order.
  // (device gets "notified" before the time integrator)
  updateTimeInfo (solState_, *analysisManager_); 
  return setParameter(comm_, artificialParameterMap_, passthroughParameterSet_, globals_, *this,
                      dependentPtrVec_, getDevices(ExternDevice::Traits::modelType()), name, val, overrideOriginal);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setParamRandomExpressionTerms2
//
// Purpose       : 
// 
// Special Notes : Handles AGAUSS, GAUSS, AUNIF, UNIF, RAND and LIMIT
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 7/30/2020
//-----------------------------------------------------------------------------
bool DeviceMgr::setParamRandomExpressionTerms2(
      const std::vector<Xyce::Analysis::SweepParam> & SamplingParams,
      bool overrideOriginal)
{
  // The call to updateTimeInfo is necessary for expresssion breaktpoints.
  // It cannot be called from the DeviceMgr::notify
  // function, as the notify functions happen in the wrong order.
  // (device gets "notified" before the time integrator)
  updateTimeInfo (solState_, *analysisManager_);
  return setParameterRandomExpressionTerms2(comm_, artificialParameterMap_, passthroughParameterSet_, globals_, *this,
      dependentPtrVec_, getDevices(ExternDevice::Traits::modelType()), 
      SamplingParams,
      masterGlobalParamDependencies,
      overrideOriginal);
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
  // return getParameter(artificialParameterMap_, globals_.paramMap, *this, *measureManager_, name, value);
}


// this macro is to support the populateSweepParam function.
// The OP argument is only used for the debug output part of this.  The functional part only needs FUNC
      // std::cout << "Parameter " << paramName << " = " << expr.get_expression() << " contains "<< #OP << std::endl;  \

#define GET_RANDOM_OP_DATA(FUNC,OP,INDEX,GPITER) \
  { \
    std::vector<Xyce::Analysis::SweepParam> tmpParams; expr.FUNC(tmpParams);  \
    if ( !(tmpParams.empty()) )  {\
      if (DEBUG_DEVICE) std::cout << "Parameter " << paramName << " contains "<< #OP << std::endl;  \
      for(int jj=0;jj<tmpParams.size();jj++) {  \
        tmpParams[jj].globalIndex = INDEX; \
        tmpParams[jj].gpIter = GPITER; \
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
    std::vector<Xyce::Analysis::SweepParam> & SamplingParams,
    int globalIndex,
    GlobalParameterMap::iterator global_param_it
    )
{
  GET_RANDOM_OP_DATA(getAgaussData,Agauss,globalIndex,global_param_it)
  GET_RANDOM_OP_DATA(getGaussData,Gauss,globalIndex,global_param_it)
  GET_RANDOM_OP_DATA(getAunifData,Aunif,globalIndex,global_param_it)
  GET_RANDOM_OP_DATA(getUnifData,Unif,globalIndex,global_param_it)
  GET_RANDOM_OP_DATA(getRandData,Rand,globalIndex,global_param_it)
  GET_RANDOM_OP_DATA(getLimitData,Limit,globalIndex,global_param_it)
}

struct SweepParam_lesser
{
  bool operator ()(Xyce::Analysis::SweepParam const& a, Xyce::Analysis::SweepParam const& b) const 
  {
    return (a.name < b.name);
  }
};

struct SweepParam_greater
{
  bool operator ()(Xyce::Analysis::SweepParam const& a, Xyce::Analysis::SweepParam const& b) const 
  {
    return (a.name > b.name);
  }
};

struct SweepParam_equal
{
  bool operator ()(Xyce::Analysis::SweepParam const& a, Xyce::Analysis::SweepParam const& b) const 
  {
    return (a.name == b.name);
  }
};

//-----------------------------------------------------------------------------
// Diagnostic output for global params, serial version.
//-----------------------------------------------------------------------------
void printOutGlobalParamsInfoSerial(
    std::string & extra,
  std::vector<Util::Expression> & expressionVec,
  std::vector<std::string> & expNameVec,
  GlobalParameterMap & global_parameter_map,
  std::vector< std::vector<entityDepend> > & deviceEntityDependVec)
{
  {
  for (int ii=0;ii< expNameVec.size(); ++ii)
  {
    std::cout << extra << "expNameVec["<<ii<<"] = " << expNameVec[ii] ;
    std::cout << " = " << expressionVec[ii].get_expression() << std::endl;

    if ( !(deviceEntityDependVec[ii].empty()))
    {
    for (int jj=0;jj< deviceEntityDependVec[ii].size();jj++)
    {
      std::string entityName;
      DeviceInstance * inst = dynamic_cast<DeviceInstance *>(deviceEntityDependVec[ii][jj].entityPtr);
      if (inst != 0)
      {
        entityName = inst->getName().getEncodedName();
      }
      else
      {
        DeviceModel * model = dynamic_cast<DeviceModel *>(deviceEntityDependVec[ii][jj].entityPtr);
        entityName = model->getName();
      }
      std::cout << extra << "  entity["<<jj<<"] = " << entityName;
      std::cout <<std::endl;
      for (int kk=0;kk< deviceEntityDependVec[ii][jj].parameterVec.size(); kk++)
      {
        std::cout << extra << "    Depend["<<kk<<"].name = " << deviceEntityDependVec[ii][jj].parameterVec[kk].name;
        std::cout <<std::endl;
      }
    }
    std::cout <<std::endl;
    }
  }
  }
}

//-----------------------------------------------------------------------------
// Diagnostic output for global params, parallel version.
//-----------------------------------------------------------------------------
void printOutGlobalParamsInfo(
    std::string & extra,
  std::vector<Util::Expression> & expressionVec,
  std::vector<std::string> & expNameVec,
  GlobalParameterMap & global_parameter_map,
  std::vector< std::vector<entityDepend> > & deviceEntityDependVec,
   Parallel::Communicator & parallel_comm)
{
  // diagnostic: print out the expNames on each processor
  //if (Parallel::is_parallel_run(parallel_comm.comm())) 
  {
    int procID = parallel_comm.procID();
    int numProc = parallel_comm.numProc();

    for (int proc=0; proc<numProc; ++proc)
    {
      if (proc == procID)
      {
        std::cout << extra << "expNameVec for proc = " << procID << std::endl;

        std::string extra2 = extra + std::string("proc ") + std::to_string(procID) + std::string(" :");

        printOutGlobalParamsInfoSerial(extra2,
            expressionVec, expNameVec, global_parameter_map, deviceEntityDependVec);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : mergeGlobals
// Purpose       : Merge the subcircuit globals into the main container of globals.
// Special Notes : Helper function for DeviceMgr::getRandomParams.  
//                 This is probably where parallelism needs to be added for 
//                 subcircuit global parameters.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/12/2023
//-----------------------------------------------------------------------------
void mergeGlobals(
  UserDefinedParams & globals_, 
  UserDefinedParams & subcktGlobals_,
   Parallel::Communicator & parallel_comm
    )
{
  //if (DEBUG_DEVICE)
  if (false)
  {
    std::vector<Util::Expression> & subcktExpressionVec = subcktGlobals_.expressionVec;
    std::vector<std::string> & subcktExpNameVec = subcktGlobals_.expNameVec;
    std::vector< std::vector<entityDepend> > & subcktDeviceEntityDependVec = subcktGlobals_.deviceEntityDependVec;
    GlobalParameterMap & global_parameter_map = subcktGlobals_.paramMap;

    std::string extra("mergeGlobal: ");
    printOutGlobalParamsInfo(extra,
        subcktExpressionVec, subcktExpNameVec, global_parameter_map, subcktDeviceEntityDependVec,parallel_comm);

  }

  if (!(subcktGlobals_.paramMap.empty()))
  {
    GlobalParameterMap & paramMap = globals_.paramMap;
    GlobalParameterMap & subcktParamMap = subcktGlobals_.paramMap;
    paramMap.insert(subcktParamMap.begin(), subcktParamMap.end());

    subcktParamMap.clear();
  }

  if (!(subcktGlobals_.expressionVec.empty()))
  {
    std::vector<Util::Expression> & expressionVec = globals_.expressionVec;
    std::vector<Util::Expression> & subcktExpressionVec = subcktGlobals_.expressionVec;
    expressionVec.insert(expressionVec.end(), subcktExpressionVec.begin(), subcktExpressionVec.end());

    subcktExpressionVec.clear();
  }

  if (!(subcktGlobals_.expNameVec.empty()))
  {
    std::vector<std::string> & expNameVec = globals_.expNameVec;
    std::vector<std::string> & subcktExpNameVec = subcktGlobals_.expNameVec;
    expNameVec.insert(expNameVec.end(), subcktExpNameVec.begin(), subcktExpNameVec.end());

    subcktExpNameVec.clear();
  }

  if (!(subcktGlobals_.deviceEntityDependVec.empty()))
  {
    std::vector< std::vector<entityDepend> > & deviceEntityDependVec = globals_.deviceEntityDependVec;
    std::vector< std::vector<entityDepend> > & subcktDeviceEntityDependVec = subcktGlobals_.deviceEntityDependVec;
    deviceEntityDependVec.insert(deviceEntityDependVec.end(), subcktDeviceEntityDependVec.begin(), subcktDeviceEntityDependVec.end());

    subcktDeviceEntityDependVec.clear();
  }

}

//-----------------------------------------------------------------------------
// Function      : sortGlobals
//
// Purpose       : Perform sort operations on the std::vectors 
//                 associated with global parameters to ensure the same order 
//                 on each processor.   Ideally these 3 objects would be 
//                 in a struct together and thus sorted once, rather than 
//                 3 separate sorts.  
//
// Special Notes : Helper function for DeviceMgr::getRandomParams.  
//                 Originally this code was inside that function.
//
//                 The three objects that are sorted in this class are:
//
//                 1. expNameVec, which is a std::vector of strings, where 
//                 each string is the name of an expression-based global 
//                 parameter.
//
//                 2. expressionVec, which is a std::vector of expressions,
//                 corresponding to expression-based global parameters.
//
//                 3. deviceEntityDependVec, which is a std::vector of std::vectors.
//                 The subordinate std::vector contains dependencies.  This vector was 
//                 set during the various DeviceEntity setup functions (specifically setParams).
//
//                 All 3 of these objects need to be in the same order.
//
//                 Other relevant objects, NOT dealt with in this function include:
//
//                 1.  globals_.paramMap, which is a std::unordered_map.  Unlike 
//                 the 3 objects above, this map contains *all* the global parameters, 
//                 including non-expression params. (i.e. params that are just a raw number).
//                 This object's size will be >= the size of the above 3 vectors.
//                 The key for this object is a string (the param name) and the value 
//                 is a double-precision number.   This "value" mostly isn't used in 
//                 the device package anymore and is a legacy of the old expression 
//                 library, which needed to have global param values explicitly pushed into it.
//                 The new expression library doesn't need this.  So, in practice, 
//                 these values are only used by things like the .PRINT line, via 
//                 "global expression Op" objects.
//
//                 2. masterGlobalParamDependencies, which is a std::unordered_map.  
//                 This contains the same information (mostly) as the deviceEntityDependVec.  
//                 But, the key is a deviceEntityPtr.  For sampling-based parameter 
//                 updates, this is the only object that gets used after setup.  But it 
//                 only contains parameters that are actually used.  The other objects, 
//                 above, don't get pruned in that way as of this writing (5/2023).
//
//                 Once the masterGlobalParamDepdendencies object is set up,
//                 the deviceEntityDependVec probably can be deleted.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/12/2023
//-----------------------------------------------------------------------------
void sortGlobals(
  std::vector<std::string> & expNameVec,
  std::vector<Util::Expression> & expressionVec,
  std::vector< std::vector<entityDepend> > & deviceEntityDependVec
    )
{
  // sort indices based on name vector
  std::vector<std::size_t> indices(expNameVec.size(),0);
  std::iota(indices.begin(), indices.end(),0);
  std::sort(indices.begin(), indices.end(), 
      [&](const int& a, const int& b) { return (expNameVec[a] < expNameVec[b]); } ) ;

  // sort the the global expression vector, using permutation based on name vector
  {
    std::vector<Util::Expression> sortedExpressionVec = expressionVec;
    std::transform(indices.begin(), indices.end(), sortedExpressionVec.begin(),
        [&](std::size_t i){ return expressionVec[i]; });
    expressionVec = sortedExpressionVec;
  }

  // sort the the global depend vector, using permutation based on name vector
  {
    std::vector< std::vector<entityDepend> > sortedDeviceEntityDependVec(deviceEntityDependVec.size());
    std::transform(indices.begin(), indices.end(), sortedDeviceEntityDependVec.begin(),
        [&](std::size_t i){ return deviceEntityDependVec[i]; });
    deviceEntityDependVec = sortedDeviceEntityDependVec;
  }

  // sort the global name vector.  Do this last.
  std::sort(expNameVec.begin(), expNameVec.end(), 
      [&](const std::string & a, const std::string & b) { return (a < b); } ) ;
}

//-----------------------------------------------------------------------------
// Function      : pruneGlobals
//
// Purpose       : Perform sort operations on the std::vectors 
//                 associated with global parameters to ensure Xyce isn't 
//                 storing and processing any unused parameters.
//
//                 Issue #306 work led to adding this function.  Some of 
//                 the params being sampled had zero dependencies in large 
//                 PDK test circuit.  
//
// Special Notes : Helper function for DeviceMgr::getRandomParams.  
//                 Originally this code was inside that function.
//
//                 Prior to this call, the size of expressionVec,expNameVec & 
//                 deviceEntityDependVec must be the same size on all procs.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/12/2023
//-----------------------------------------------------------------------------
void pruneGlobals(
    std::vector<std::string> & expNameVec,
    std::vector<Util::Expression> & expressionVec,
    std::vector< std::vector<entityDepend> > & deviceEntityDependVec,
    Parallel::Communicator & parallel_comm
    )
{
  std::vector<int> pruneIndices(expressionVec.size(),0);
  int count=0;
  for (int ii=0;ii<expressionVec.size();ii++)
  {
    bool empty = deviceEntityDependVec[ii].empty();
    if (Parallel::is_parallel_run(parallel_comm.comm()))
    {
      int tmpEmpty = empty?1:0; int globalEmpty = 0;
      parallel_comm.minAll(&tmpEmpty, &globalEmpty, 1);
      empty = (globalEmpty==0)?false:true;
    }
    if ( empty ) { pruneIndices[ii] = 1; }
    else { count++; }
  }

  std::vector<std::string> tmpExpNameVec;
  std::vector<Util::Expression> tmpExpressionVec;
  std::vector< std::vector<entityDepend> > tmpDeviceEntityDependVec;
  tmpExpNameVec.reserve(count);
  tmpExpressionVec.reserve(count);
  tmpDeviceEntityDependVec.reserve(count);

  for (int ii=0; ii<pruneIndices.size(); ii++)
  {
    if ( pruneIndices[ii] == 0 )
    {
      tmpExpNameVec.push_back(expNameVec[ii]);
      tmpExpressionVec.push_back(expressionVec[ii]);
      tmpDeviceEntityDependVec.push_back(deviceEntityDependVec[ii]);
    }
  }

  expNameVec = tmpExpNameVec;
  expressionVec = tmpExpressionVec;
  deviceEntityDependVec = tmpDeviceEntityDependVec;
}

//-----------------------------------------------------------------------------
// Function      : resolveParam
//
// Purpose       : Resolve (global) parameter that depencies on other global
//                 parameters.
//
// Special Notes : This is only needed in parllel, for "subcircuit global"
//                 parameters.   These parameters (as of this writing) only
//                 happen when applying "local variation" with expression-based
//                 random operators.  This process turns subcircuit .params
//                 (which normally are constant) into .global_params
//                 (which are allowed to change).
//
//                 This is a helper function for broadcastSubcktGlobals.
//
//                 Ordinarily, parameters are fully resolved in the parser,
//                 specifically in functions like CircuitContext::resolve
//                 and its subordinate functions.
//
//                 However, in Xyce's current design it is necessary for the full
//                 set of global parameters to be known on every processor, at
//                 least for UQ methods.  At the end of parsing, only proc 0
//                 has a full set, and they must be broadcast to the other
//                 processors.
//
//                 When they are broadcast, the expressions are sent as strings.
//                 The strings that are sent are generated from fully resolved
//                 expressions.  So, *most* external dependencies just become
//                 part of that string.  Any .param will be converted into a raw
//                 number.  But, for global_params (i.e. variable) it is necessary
//                 for them to be attached as an external AST.  These external
//                 global param dependencies are not included in the newly
//                 generated string, so they must be re-resolved after they've
//                 been received by the non-0  processor.
//
//                 If everything has been done properly, they should only
//                 depend on globals, and they should be able to find those
//                 globals in the device manager globals_ containers.  Also,
//                 to get to this point they must have been successfully
//                 resolved on proc 0 during parsing.  So, if they don't resolve
//                 successfully here, it won't be because of a netlist problem.
//                 It will instead be a mistake in the code design.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/16/2023
//-----------------------------------------------------------------------------
bool resolveParam(
  UserDefinedParams &  globals_,
  const std::string & paramName,
  Util::Expression & expression,
  const std::vector<std::string> & unresolvedParams
    )
{
  std::vector<Util::Expression> & expressionVec = globals_.expressionVec;
  std::vector<std::string> & expNameVec = globals_.expNameVec;
  std::vector< std::vector<entityDepend> > & deviceEntityDependVec = globals_.deviceEntityDependVec;
  GlobalParameterMap & global_parameter_map = globals_.paramMap;

  for (int ii=0;ii<unresolvedParams.size();ii++)
  {
    GlobalParameterMap::iterator global_param_it = global_parameter_map.find(unresolvedParams[ii]);
    if (global_param_it != global_parameter_map.end())
    {
      std::vector<std::string>::iterator name_it =
        std::find(expNameVec.begin(), expNameVec.end(), unresolvedParams[ii]);
      if (name_it != expNameVec.end())
      {
        int index = std::distance( expNameVec.begin(), name_it );
        Util::Expression & expToBeAttached = expressionVec[index];
        expression.attachParameterNode(unresolvedParams[ii], expToBeAttached);
      }
      else
      {
        //enumParamType paramType=DOT_PARAM; // possibly this should be type GLOBAL
        //expression.make_constant(unresolvedParams[ii], global_param_it->second, paramType);
        expression.make_constant(unresolvedParams[ii], global_param_it->second);
      }
    }
    else
    {
      Report::UserError() << "Could not find " << unresolvedParams[ii] << " in " << paramName;
      return false;
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::broadcastSubcktGlobals
//
// Purpose       : Send the full list of subcircuit global parameters from
//                 proc 0 to the other processors.
//
// Special Notes : Helper function for the DeviceMgr::getRandomParams function.
//
//                 Historically, Xyce only allowed global parameters to come
//                 from the top level of the netlist.  global_params from
//                 subcircuits was not allowed.
//
//                 However, in order to support local variation for UQ
//                 methods, it is necessary to support "global" parameters
//                 from subcircuits.  A better name would be "variable"
//                 parameters.  Subcircuit global parameters are parameters
//                 that were originally declared using a (non-global)
//                 .param statement, but that got converted to "global"
//                 during parsing, due to dependence on local random
//                 operators such as AGAUSS.
//
//                 In the current parser design, the top level of the netlist
//                 is (mostly) processed on proc 0 during the first parser pass.
//                 These top-level parameters are broadcast to the other
//                 processors at the end of pass 1.
//
//                 Subcircuits, in contrast, are expanded during pass 2.
//                 With respect to parameter resolution (handled in the
//                 CircuitContext::resolve function), *all* subcircuits are
//                 resolved on proc 0.  On other processors, expanded
//                 subcirrcuits are resolved as needed.  The result is that
//                 proc 0 has a complete vector of all the subcircuit global
//                 params, and other processors have incompletely vectors.
//                 Currently, the UQ methods in Xyce assume that every
//                 processor has a full vector.  If the vectors are
//                 inconsistent from processor to processor, the current
//                 code will fail.
//
//                 This function performs a broadcast to ensure that all
//                 processors have the same vectors of subcircuit params.
//
//                 After to this call, the size of expressionVec,expNameVec &
//                 deviceEntityDependVec must be the same size on all procs.
//
//                 NOTE: possible fragile aspects of this function include:
//
//                 1. It assumes that the ONLY use case that will produce
//                 subcircuit globals is local variation via random operators.
//                 As such, it assumes that the only subcircuit params we need
//                 to broadcast are of EXPR type.  So, it is oriented towards
//                 the expNameVec and expressionVec.  It adds params to the
//                 paramMap as well (which can also contain DBL type params),
//                 but it doesn't broadcast from it.
//
//                 2. Expressions are broadcast as strings, and then reparsed
//                 via the expression library. They may need a "resolve" step.
//                 For simple cases (like expressions that are just agauss()),
//                 this isn't necessary.  If they need a resolve step, I think
//                 that should work from here.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/16/2023
//-----------------------------------------------------------------------------
void DeviceMgr::broadcastSubcktGlobals(Parallel::Communicator & parallel_comm)
{
  if (Parallel::is_parallel_run(parallel_comm.comm()))
  {

    std::vector<Util::Expression> & subcktExpressionVec = subcktGlobals_.expressionVec;
    std::vector<std::string> & subcktExpNameVec = subcktGlobals_.expNameVec;
    std::vector< std::vector<entityDepend> > & subcktDeviceEntityDependVec = subcktGlobals_.deviceEntityDependVec;

    int procID = parallel_comm.procID();

    std::vector<std::string> tmpExpNames;
    std::vector<std::string> tmpExpressionStringVec;

    //  broadcast subckt expNameVec
    {
      int byteCount = 0;
      int bsize=0;

      if (procID == 0)
      {
        int size = subcktExpNameVec.size();
        byteCount += sizeof( int );
        std::vector<std::string>::const_iterator it = subcktExpNameVec.begin();
        std::vector<std::string>::const_iterator end = subcktExpNameVec.end();
        for (; it != end; ++it)
        {
          byteCount += sizeof( int );
          byteCount += it->length();
        }
        bsize = byteCount + 100;
      }

      parallel_comm.bcast( &bsize, 1, 0 );
      char * nameBuffer = new char[bsize];

      int pos =0;
      if (procID == 0)
      {
        int numExpNameVecSize = subcktExpNameVec.size();
        parallel_comm.pack( &numExpNameVecSize, 1, nameBuffer, bsize, pos );

        for (int ii=0;ii<subcktExpNameVec.size();ii++)
        {
          int length = subcktExpNameVec[ii].length();
          parallel_comm.pack( &length, 1, nameBuffer, bsize, pos );
          parallel_comm.pack( subcktExpNameVec[ii].c_str(), length, nameBuffer, bsize, pos );
        }
      }
      parallel_comm.bcast( nameBuffer, bsize, 0);

      if (procID != 0)
      {
        int numExpNameVecSize = 0;
        parallel_comm.unpack( nameBuffer, bsize, pos, &numExpNameVecSize, 1 );
        for (int ii=0; ii<numExpNameVecSize; ++ii)
        {
          int length=0;
          parallel_comm.unpack( nameBuffer, bsize, pos, &length, 1 );
          std::string paramName(std::string( ( nameBuffer + pos ), length ));
          pos += length;
          tmpExpNames.push_back(paramName);
        }
      }
      delete [] nameBuffer;
    }

    //  broadcast subckt expressionVec.
    {
      int byteCount = 0;
      int bsize=0;

      if (procID == 0)
      {
        int expressionVecSize = subcktExpressionVec.size();
        tmpExpressionStringVec.resize(expressionVecSize);
        {
          for(int ii=0;ii<expressionVecSize;ii++)
          {
            std::string exprString;
            subcktExpressionVec[ii].generateExpressionString(exprString); // regenerates the string, gets rid of .param dependencies
            tmpExpressionStringVec[ii] = exprString;
            byteCount += sizeof( int );
            byteCount += exprString.length();
          }
        }
        bsize = byteCount + 100;
      }

      parallel_comm.bcast( &bsize, 1, 0 );
      char * exprBuffer = new char[bsize];

      int pos =0;
      if (procID == 0)
      {
        int expressionVecSize = subcktExpressionVec.size();
        parallel_comm.pack( &expressionVecSize, 1, exprBuffer, bsize, pos );

        for (int ii=0;ii<expressionVecSize;ii++)
        {
          int length = tmpExpressionStringVec[ii].length();
          parallel_comm.pack( &length, 1, exprBuffer, bsize, pos );
          parallel_comm.pack( tmpExpressionStringVec[ii].c_str(), length, exprBuffer, bsize, pos );
        }
      }
      parallel_comm.bcast( exprBuffer, bsize, 0);

      if (procID != 0)
      {
        int expressionVecSize = 0;
        parallel_comm.unpack( exprBuffer, bsize, pos, &expressionVecSize, 1 );
        for (int ii=0; ii<expressionVecSize; ++ii)
        {
          int length=0;
          parallel_comm.unpack( exprBuffer, bsize, pos, &length, 1 );
          std::string exprString(std::string( ( exprBuffer + pos ), length ));
          pos += length;
          tmpExpressionStringVec.push_back(exprString);
        }
      }

      delete [] exprBuffer;
    }

    // now figure it all out
    if (procID != 0)
    {
      for (int ii=0;ii<tmpExpNames.size();ii++)
      {
        std::vector<std::string>::iterator name_it =
          std::find(subcktExpNameVec.begin(), subcktExpNameVec.end(), tmpExpNames[ii]);

        if (name_it == subcktExpNameVec.end())
        {
          subcktExpNameVec.push_back(tmpExpNames[ii]);
          std::string & expressionString = tmpExpressionStringVec[ii];
          Util::Expression expression(expressionGroup_, expressionString);
          subcktExpressionVec.push_back(expression);
          subcktDeviceEntityDependVec.push_back(std::vector<entityDepend>(0));
        }
      }

      for (int ii=0;ii<subcktExpNameVec.size();ii++)
      {
        Util::Expression & expression = subcktExpressionVec[ii];

        // this code is to check if things are resolved or not
        const std::vector<std::string> unresolvedParams = expression.getUnresolvedParams();
        if (!(unresolvedParams.empty()))
        {
          bool resolved = false;
          std::string & nameOfParamBeingResolved = subcktExpNameVec[ii];
          if (!(resolveParam( globals_,nameOfParamBeingResolved, expression, unresolvedParams)))
          {
            resolved = resolveParam( subcktGlobals_,nameOfParamBeingResolved, expression, unresolvedParams);
          }
          else
          {
            resolved = true;
          }
          if (!resolved)
          {
            Report::UserError() << "Failed to resolve " << subcktExpNameVec[ii];
          }
        }

        double val;
        expression.evaluateFunction(val);
        subcktGlobals_.paramMap[tmpExpNames[ii]] = val;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getRandomParams
//
// Purpose       : To conduct Hspice-style sampling analysis, it is necessary
//                 to gather all the parameters (global and/or device) which
//                 depend upon expressions containing operators such as
//                 AGAUSS, GAUSS, AUNIF, UNIF, RAND and LIMIT
//
// Special Notes : Setup function, called once.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/28/2020
//-----------------------------------------------------------------------------
void DeviceMgr::getRandomParams(std::vector<Xyce::Analysis::SweepParam> & SamplingParams,
   Parallel::Communicator & parallel_comm)
{
  Stats::StatTop _samplingStat("Setup Expression Sampling Params");
  Stats::TimeBlock _samplingTimer(_samplingStat);

  std::vector<Util::Expression> & expressionVec  = globals_.expressionVec;
  std::vector<std::string> & expNameVec = globals_.expNameVec;
  GlobalParameterMap & global_parameter_map = globals_.paramMap;
  std::vector< std::vector<entityDepend> > & deviceEntityDependVec = globals_.deviceEntityDependVec;

  if (false) // useful debug
  {
    std::string extra("DeviceMgr::getRandomParams: ");
    printOutGlobalParamsInfo(extra,expressionVec, expNameVec, global_parameter_map, deviceEntityDependVec,parallel_comm);
  }

  if (false)
  {
    std::vector<Util::Expression> & subcktExpressionVec = subcktGlobals_.expressionVec;
    std::vector<std::string> & subcktExpNameVec = subcktGlobals_.expNameVec;
    std::vector< std::vector<entityDepend> > & subcktDeviceEntityDependVec = subcktGlobals_.deviceEntityDependVec;
    pruneGlobals(subcktExpNameVec, subcktExpressionVec, subcktDeviceEntityDependVec,parallel_comm);
  }

  broadcastSubcktGlobals(parallel_comm);
  mergeGlobals( globals_, subcktGlobals_,parallel_comm);
  sortGlobals( expNameVec, expressionVec, deviceEntityDependVec);

  //  this minAll should not be necessary, but just being extra-safe.
  bool emptyExprVec = expressionVec.empty();
  if (Parallel::is_parallel_run(parallel_comm.comm()))
  {
    int tmpEmpty = emptyExprVec?1:0; int globalEmpty = 0;
    parallel_comm.minAll(&tmpEmpty, &globalEmpty, 1);
    emptyExprVec = (globalEmpty==0)?false:true;
  }

  if (! ( emptyExprVec) )
  {
    pruneGlobals(expNameVec, expressionVec, deviceEntityDependVec,parallel_comm); // is this necessary?

    for (int ii=0;ii<expressionVec.size();ii++)
    {
      {
        GlobalParameterMap::iterator global_param_it = global_parameter_map.find(expNameVec[ii]);
        populateSweepParam(expressionVec[ii], expNameVec[ii], SamplingParams, ii, global_param_it);

        {
          // new code to set up the "master" global parameter dependency map.
          // May need to separate this code when refactoring the above code for parallel.
          //
          // The goal for the masterGlobalParamDependencies container is
          // to only contain unique objects.
          //
          // So in this container, each entity is unique, ie only appears once.
          // This is enforced by having the container be a std::unordered_map,
          // and using the entityPtr as the key, or "first" part of the pair.
          //
          // The parameter vector is the value, or "second" part of the pair.
          // So each parameterVec is associated with a specific entity (model or instance).
          //
          // The elements of this parameter vector are also unique.  So,
          // each parameter in the vector will appear at most once. This uniqueness
          // is enforced by the code below.
          //
          // The params in the parameterVec ONLY include parameters that
          // have one or more .param/.global_param dependencies.  Any params
          // that lack these dependencies have been excluded.  This screening
          // happened upstream, in the device entity, during the DeviceEntity::setParams
          // and subordinate  DeviceEntity::setDependentParameter functions,
          // when the globals_.deviceEntityDependVec object was set up.
          std::vector<entityDepend> & edVec = globals_.deviceEntityDependVec[ii];
          for (int ie=0;ie< edVec.size(); ie++)
          {
            entityDepend & ed = edVec[ie];

            std::unordered_map<DeviceEntity *, std::vector<Depend> >::iterator mgpDepIter =
              masterGlobalParamDependencies.find(ed.entityPtr);

            if(mgpDepIter!=masterGlobalParamDependencies.end())
            {
              std::vector<Depend> & oldParamVec = mgpDepIter->second;
              oldParamVec.insert(oldParamVec.end(), ed.parameterVec.begin(), ed.parameterVec.end());

              // sort and remove duplicates, if any:
              std::sort(oldParamVec.begin(), oldParamVec.end(), Depend_greater());
              std::vector<Depend>::iterator it = std::unique(oldParamVec.begin(), oldParamVec.end(), Depend_equal());
              oldParamVec.resize( std::distance (oldParamVec.begin(), it ));
            }
            else
            {
              masterGlobalParamDependencies[ed.entityPtr] = ed.parameterVec;
            }
          }
        }
      }
    }

    // in parallel, the order is not guaranteed from processor to processor
    // do in serial as well, for testing purposes.
    //if (Parallel::is_parallel_run(parallel_comm.comm()))
    {
      if ( !(SamplingParams.empty()) )
      {
        std::sort(SamplingParams.begin(), SamplingParams.end(), SweepParam_lesser());
      }
    }
  }

  // now do device params
  {
    std::vector<Xyce::Analysis::SweepParam>  deviceSamplingParams;
    std::vector<Xyce::Analysis::SweepParam>  deviceModelSamplingParams;

    // models
    for (int ii=0;ii<modelVector_.size();ii++)
    {
      if (!(*(modelVector_[ii])).getDependentParams().empty())
      {
        const std::vector<Depend> & depVec = (*(modelVector_[ii])).getDependentParams();
        std::string entityName = (*(modelVector_[ii])).getName();
        for (int jj=0;jj<depVec.size();jj++)
        {
          std::string fullName = entityName + Xyce::Util::separator + depVec[jj].name;
          GlobalParameterMap::iterator global_param_it = global_parameter_map.end();
          populateSweepParam( *(depVec[jj].expr), fullName, deviceModelSamplingParams, -1, global_param_it);
        }
      }
    }

    std::sort(deviceModelSamplingParams.begin(), deviceModelSamplingParams.end(), SweepParam_lesser());

    // if parallel, get a combined vector of random model params from all processors.
    // ERK Note: this parallel reduction may be unneccessary.  As of this writing,
    // I am not 100% sure what the distribution strategy is for device models.
    //
    // So, I am writing this assuming that
    //   (1) models are not global,
    //   (2) the same model may appear on 2 or more processors
    //
    if (Parallel::is_parallel_run(parallel_comm.comm()))
    {
      int procID = parallel_comm.procID();
      int numProc = parallel_comm.numProc();
      int numDevParam = deviceModelSamplingParams.size();
      int numDevParamTotal = 0;
      parallel_comm.sumAll(&numDevParam, &numDevParamTotal, 1);

      std::vector<Xyce::Analysis::SweepParam> externalDeviceSamplingParams;
      if (numDevParamTotal > 0)
      {
        // Count the bytes for packing
        int byteCount = 0;

        // First count the size of the local vector
        byteCount += sizeof(int);

        // Now count all the SweepParams in the vector
        for (int ii=0;ii<deviceModelSamplingParams.size();ii++)
        {
          byteCount += Xyce::packedByteCount(deviceModelSamplingParams[ii]);
        }

        for (int proc=0; proc<numProc; ++proc)
        {
          parallel_comm.barrier();

          // Broadcast the buffer size for this processor.
          int bsize=0;
          if (proc == procID) { bsize = byteCount+100; }
          parallel_comm.bcast( &bsize, 1, proc );

          // Create buffer.
          int pos = 0;
          char * paramBuffer = new char[bsize];

          if (proc == procID)
          {
            parallel_comm.pack( &numDevParam, 1, paramBuffer, bsize, pos );
            for (int ii=0;ii<deviceModelSamplingParams.size();ii++)
            {
              Xyce::pack(deviceModelSamplingParams[ii], paramBuffer, bsize, pos, &parallel_comm);
            }
            parallel_comm.bcast( paramBuffer, bsize, proc );
          }
          else
          {
            parallel_comm.bcast( paramBuffer, bsize, proc );
            int devParams = 0;
            parallel_comm.unpack( paramBuffer, bsize, pos, &devParams, 1 );
            Xyce::Analysis::SweepParam tmpSweepParam;
            for (int i=0; i<devParams; ++i)
            {
              Xyce::unpack(tmpSweepParam, paramBuffer, bsize, pos, &parallel_comm);
              externalDeviceSamplingParams.push_back( tmpSweepParam );
            }
          }
          delete [] paramBuffer;
        }
      }

      deviceModelSamplingParams.insert
        (deviceModelSamplingParams.end(), externalDeviceSamplingParams.begin(), externalDeviceSamplingParams.end());

      // sort the new vector and then eliminate duplicates.
      std::sort(deviceModelSamplingParams.begin(), deviceModelSamplingParams.end(), SweepParam_lesser());
      std::vector<Xyce::Analysis::SweepParam>::iterator it = std::unique(deviceModelSamplingParams.begin(), deviceModelSamplingParams.end(), SweepParam_equal());
      deviceModelSamplingParams.resize( std::distance (deviceModelSamplingParams.begin(), it ));
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
          std::string fullName = entityName + Xyce::Util::separator + depVec[jj].name;
          GlobalParameterMap::iterator global_param_it = global_parameter_map.end();
          populateSweepParam( *(depVec[jj].expr), fullName, deviceSamplingParams, -2, global_param_it);
        }
      }
    }

    std::sort(deviceSamplingParams.begin(), deviceSamplingParams.end(), SweepParam_lesser());

    // if parallel, get a combined vector of random instance params from all processors.
    if (Parallel::is_parallel_run(parallel_comm.comm()))
    {
      int procID = parallel_comm.procID();
      int numProc = parallel_comm.numProc();
      int numDevParam = deviceSamplingParams.size();
      int numDevParamTotal = 0;
      parallel_comm.sumAll(&numDevParam, &numDevParamTotal, 1);

      std::vector<Xyce::Analysis::SweepParam> externalDeviceSamplingParams;
      if (numDevParamTotal > 0)
      {
        // Count the bytes for packing
        int byteCount = 0;

        // First count the size of the local vector
        byteCount += sizeof(int);

        // Now count all the SweepParams in the vector
        for (int ii=0;ii<deviceSamplingParams.size();ii++)
        {
          byteCount += Xyce::packedByteCount(deviceSamplingParams[ii]);
        }

        for (int proc=0; proc<numProc; ++proc)
        {
          parallel_comm.barrier();

          // Broadcast the buffer size for this processor.
          int bsize=0;
          if (proc == procID) { bsize = byteCount+100; }
          parallel_comm.bcast( &bsize, 1, proc );

          // Create buffer.
          int pos = 0;
          char * paramBuffer = new char[bsize];

          if (proc == procID)
          {
            parallel_comm.pack( &numDevParam, 1, paramBuffer, bsize, pos );
            for (int ii=0;ii<deviceSamplingParams.size();ii++)
            {
              Xyce::pack(deviceSamplingParams[ii], paramBuffer, bsize, pos, &parallel_comm);
            }
            parallel_comm.bcast( paramBuffer, bsize, proc );
          }
          else
          {
            parallel_comm.bcast( paramBuffer, bsize, proc );
            int devParams = 0;
            parallel_comm.unpack( paramBuffer, bsize, pos, &devParams, 1 );
            Xyce::Analysis::SweepParam tmpSweepParam;
            for (int i=0; i<devParams; ++i)
            {
              Xyce::unpack(tmpSweepParam, paramBuffer, bsize, pos, &parallel_comm);
              externalDeviceSamplingParams.push_back( tmpSweepParam );
            }
          }
          delete [] paramBuffer;
        }
      }

      deviceSamplingParams.insert
        (deviceSamplingParams.end(), externalDeviceSamplingParams.begin(), externalDeviceSamplingParams.end());

      // sort the new vector, so they match on each proc.
      std::sort(deviceSamplingParams.begin(), deviceSamplingParams.end(), SweepParam_lesser());
    }

    // now insert the model/instance parameters into the master vector.
    SamplingParams.insert
      (SamplingParams.end(), deviceModelSamplingParams.begin(), deviceModelSamplingParams.end());

    SamplingParams.insert
      (SamplingParams.end(), deviceSamplingParams.begin(), deviceSamplingParams.end());
  }

  if ( !(SamplingParams.empty()) ) expressionBasedSamplingEnabled_ = true;
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
  Linear::Vector * lastStoVectorPtr,
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
  externData_.lastStoVectorPtr = lastStoVectorPtr;

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
  externData_.lastStoVectorRawPtr = &((*externData_.lastStoVectorPtr)[0]);

  updateSolutionDependentParameters_();

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
    externData_.dQdxMatrixPtr->print(dout());
    dout() << std::endl;
    dout() << section_divider << std::endl;
    dout() <<  "F-matrix: nonlinear iteration = " << newtonIter << "\n";
    externData_.dFdxMatrixPtr->print(dout());
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
  Linear::Vector * tmpLastStoVectorPtr,
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
  externData_.lastStoVectorPtr = tmpLastStoVectorPtr;
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
    externData_.daeQVectorPtr->print(std::cout);
    dout() << std::endl;
    dout() <<  "F-vector: nonlinear iteration = " << newtonIter << "\n";
    externData_.daeFVectorPtr->print(std::cout);
    dout() << std::endl;

    if (devOptions_.voltageLimiterFlag)
    {
      dout() << "\n\n  dFdxdVp vector: nonlinear iteration = " << newtonIter << "\n";
      externData_.dFdxdVpVectorPtr->print(std::cout);
      dout() << std::endl;
      dout() << "\n\n  dQdxdVp vector: nonlinear iteration = " << newtonIter << "\n";
      externData_.dQdxdVpVectorPtr->print(std::cout);
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
  addGlobalParameter(globals_, param, expressionGroup_);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addGlobalPars()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/12/2023
//-----------------------------------------------------------------------------
void DeviceMgr::addGlobalPars(const Util::UParamList & ioGlobalParams)
{
  addGlobalParameters(globals_, ioGlobalParams, expressionGroup_);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addSubcktGlobalPars()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/12/2023
//-----------------------------------------------------------------------------
void DeviceMgr::addSubcktGlobalPars(const Util::UParamList & ioGlobalParams)
{
  addGlobalParameters(subcktGlobals_, ioGlobalParams, expressionGroup_);
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
  return findGlobalParameter(globals_.paramMap, parName);
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
// Function       : DeviceMgr::updateSolutionDependentParameters_
//
// Purpose        : This function updates all solution dependent parameters for
//                  the current time step.
//
//                  ERK. This function is a little misleading.  Most of 
//                  the devices which can have solution dependent parameters 
//                  (Bsrc, resistor, capacitor) handle them entirely inside 
//                  functions such as updatePrimaryState or 
//                  updateIntermediateVars.  So, those devices are excluded
//                  from the loop below.
//
//                  The main device in the "open models" directory, which 
//                  can have a solution dependent variable, but that does not 
//                  handle it for itself is the linear mutual inductor device.
//
//                  As a result, in most circuits this function doesn't
//                  do anything.
//
// Scope          : private
// Creator        : Eric Keiter
// Creation Date  : 05/11/2023
//----------------------------------------------------------------------------
void DeviceMgr::updateSolutionDependentParameters_()
{
  if (firstDependent_)
  {
    firstDependent_ = false;
    setupDependentEntities();
  }

  if (!(solnDepEntityPtrVec_.empty()))
  {
    EntityVector::iterator iter;
    EntityVector::iterator begin = solnDepEntityPtrVec_.begin();
    EntityVector::iterator end = solnDepEntityPtrVec_.end();
    for (iter=begin; iter!=end;++iter)
    {
      bool changed = (*iter)->updateSolutionDependentParameters ();
      if (changed)
      {
        (*iter)->processParams();
        (*iter)->processInstanceParams();
      }
    }
  }

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
// Function      : DeviceMgr::loadFreqBVectorsforSources
// Purpose       : This function loads the B-vector contributions for sources.
// Special Notes : This is only done for linear sources.
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool DeviceMgr::loadFreqBVectorsforSources(double frequency,
                                           std::vector<Util::FreqVecEntry>& BVecEntries)

{
  bool bsuccess = true;

  IndependentSourceVector::iterator it = independentSourceVector_.begin();
  IndependentSourceVector::iterator end = independentSourceVector_.end();
  for ( ; it != end; ++it)
  {
    if ( (*it)->isLinearDevice() )
    {
      (*it)->loadFreqBVector(frequency, BVecEntries);
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadTwoLevelVsrcs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2022
//-----------------------------------------------------------------------------
bool DeviceMgr::loadTwoLevelVsrcs (
       const std::vector<std::string> & names, 
       Xyce::Linear::Vector * fptr,
       Xyce::Linear::Vector * bptr,
       Xyce::Linear::Vector * sol)
{
  bool bsuccess = true;
  int loadType = Xyce::Device::ALL;

  Xyce::Linear::Vector * saveF = externData_.daeFVectorPtr;
  Xyce::Linear::Vector * saveB = externData_.daeBVectorPtr;
  Xyce::Linear::Vector * saveSol = externData_.nextSolVectorPtr;

  externData_.daeFVectorPtr    = fptr;
  externData_.daeBVectorPtr    = bptr;
  externData_.nextSolVectorPtr = sol;
  externData_.daeFVectorRawPtr    = &((*externData_.daeFVectorPtr)[0]);
  externData_.daeBVectorRawPtr    = &((*externData_.daeBVectorPtr)[0]);
  externData_.nextSolVectorRawPtr = &((*externData_.nextSolVectorPtr)[0]);

  for (int ii=0;ii<names.size();ii++)
  {
    DeviceEntity * device_entity = getDeviceEntity(names[ii]);
    if (device_entity != 0)
    {
      Vsrc::Instance * vsrc = dynamic_cast<Vsrc::Instance *>(device_entity);
      if (vsrc != 0)
      {
        vsrc->loadDAEFVector();
        vsrc->loadDAEBVector();
      }
    }
  }

  externData_.daeFVectorPtr = saveF;
  externData_.daeBVectorPtr = saveB;
  externData_.nextSolVectorPtr = saveSol;
  externData_.daeFVectorRawPtr    = &((*externData_.daeFVectorPtr)[0]);
  externData_.daeBVectorRawPtr    = &((*externData_.daeBVectorPtr)[0]);
  externData_.nextSolVectorRawPtr = &((*externData_.nextSolVectorPtr)[0]);

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
  Linear::Vector * tmpLastStoVectorPtr,
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
  externData_.lastStoVectorPtr = tmpLastStoVectorPtr;

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
  std::vector<Util::Expression>::iterator globalExp_i = globals_.expressionVec.begin();
  std::vector<Util::Expression>::iterator globalExp_end = globals_.expressionVec.end();
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

  if (externData_.lastStoVectorPtr != 0)
  {
    externData_.lastStoVectorRawPtr = &((*externData_.lastStoVectorPtr)[0]);
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

  std::vector<Util::Expression> & globalExpressionsVec = globals_.expressionVec;

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
  char * buf, int bsize, int & pos, Parallel::Communicator * comm, bool pack) const
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
bool DeviceMgr::restoreRestartData(char * buf, int bsize, int & pos, Parallel::Communicator * comm, bool pack)
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
// Function      : getMutualInductor
// Purpose       : helper function for DeviceMgr::getMutualInductorDeviceInstance
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/12/2021
//-----------------------------------------------------------------------------
DeviceInstance * getMutualInductor (
    const std::string & deviceName,
      const InstanceVector &MIdevices, int & index
      ) 
{
  DeviceInstance * device_instance = 0;
  index=-1;

  for (InstanceVector::const_iterator instance_it=MIdevices.begin(); 
      instance_it!=MIdevices.end(); instance_it++)
  {
    std::vector<std::string> inductorNames = (*instance_it)->getInductorNames();
    std::vector<std::string>::iterator it = std::find(inductorNames.begin(), inductorNames.end(), deviceName);
    if (it != inductorNames.end())
    {
      index = std::distance( inductorNames.begin(), it );
      device_instance = const_cast<Xyce::Device::DeviceInstance *>(*instance_it);
    }
  }

  return device_instance;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getMutualInductorDeviceInstance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/30/2021
//-----------------------------------------------------------------------------
DeviceInstance * DeviceMgr::getMutualInductorDeviceInstance (
      const std::string &full_param_name, int & index
      ) const
{
  std::string entity_name = Xyce::Util::entityNameFromFullParamName(full_param_name).getEncodedName();
  std::string param_name  = Xyce::Util::paramNameFromFullParamName(full_param_name);
  InstanceName foo(entity_name);
  DeviceInstance * device_instance = 0;
  index=-1;

  if ( foo.getDeviceName()[0] == 'L')
  {
    bool found=false;
    const InstanceVector &MILdevices = getDevices(MutIndLin::Traits::modelType());
    if (MILdevices.size() > 0)
    {
      device_instance = getMutualInductor(foo.getDeviceName(),MILdevices,index);
      found = (index!=-1);
    }

    if (!found)
    {
      const InstanceVector &MINdevices = getDevices(MutIndNonLin::Traits::modelType());
      if (MINdevices.size() > 0)
      {
        device_instance = getMutualInductor(foo.getDeviceName(),MINdevices,index);
        found = (index!=-1);
      }
    }

    if (!found)
    {
      const InstanceVector &MINdevices = getDevices(MutIndNonLin2::Traits::modelType());
      if (MINdevices.size() > 0)
      {
        device_instance = getMutualInductor(foo.getDeviceName(),MINdevices,index);
        found = (index!=-1);
      }
    }
  }

  return device_instance;
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
  options_manager.addOptionsProcessor("PARSER", IO::createRegistrationOptions(device_manager, &DeviceMgr::setParserOptions));

  options_manager.addCommandProcessor("SAMPLING", IO::createRegistrationOptions(device_manager, &DeviceMgr::setSamplingParams));
  options_manager.addCommandProcessor("EMBEDDEDSAMPLING", IO::createRegistrationOptions(device_manager, &DeviceMgr::setEmbeddedSamplingParams));
  options_manager.addCommandProcessor("PCE", IO::createRegistrationOptions(device_manager, &DeviceMgr::setPCEParams));
  options_manager.addOptionsProcessor("DIAGNOSTIC", IO::createRegistrationOptions(device_manager, &DeviceMgr::setDiagnosticOptions));

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

  if (!op)
  {
    if (DEBUG_DEVICE) { Report::DevelWarning() << "Xyce::Device::getParamAndReduce.  Unable to find parameter " << name; }
    else { Report::UserError() << "Xyce::Device::getParamAndReduce.  Unable to find parameter " << name; }
  }

  return op;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getParamAndReduce
// Purpose       : Returns the current value of a named parameter.
//
// Special Notes : This works in parallel.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/26/23
//-----------------------------------------------------------------------------
bool getParamAndReduce(
  Parallel::Machine     comm,
  const DeviceMgr &     device_manager,
  const std::string &   name,
  std::complex<double> &  value)
{
  Util::Op::Operator *op = device_manager.getOp(comm, name);
  value = op ? (*op)(comm, Util::Op::OpData()) : 0.0;

  if (!op)
  {
    if (DEBUG_DEVICE) { Report::DevelWarning() << "Xyce::Device::getParamAndReduce.  Unable to find parameter " << name; }
    else { Report::UserError() << "Xyce::Device::getParamAndReduce.  Unable to find parameter " << name; }
  }

  return op;
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
  UserDefinedParams &           globals,               ///< global variables
  DeviceMgr &                   device_manager,
  EntityVector &                dependent_entity_vector, // dependentPtrVec_, if called from DeviceMgr::setParam
  const InstanceVector &        extern_device_vector,
  const std::string &           name,
  double                        value,
  bool                          override_original)
{
  bool bsuccess = true, success = true;
  GlobalParameterMap &  global_parameter_map = globals.paramMap;
  std::vector<Util::Expression> & global_expressions = globals.expressionVec;
  std::vector<std::string> & global_exp_names = globals.expNameVec;

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
      // find the expression.  ERK.  This should be the "real" place to update.
      std::string tmpName = name; Util::toUpper(tmpName);
      std::vector<std::string>::iterator name_it = std::find(global_exp_names.begin(), global_exp_names.end(), tmpName);

      if (name_it == global_exp_names.end())
      {
        Report::UserError() << "Xyce::Device::setParameter. Unable to find parameter " << name;
      }
      else
      {
        int globalIndex = std::distance (global_exp_names.begin(), name_it );
        Util::Expression &expression = global_expressions[globalIndex];
        double exprValue(0.0);
        expression.evaluateFunction(exprValue);

        if ((*global_param_it).second != value || exprValue != value)
        {
          (*global_param_it).second = value;

          expression.setValue(value);

          EntityVector::iterator it = dependent_entity_vector.begin();
          EntityVector::iterator end = dependent_entity_vector.end();
          for ( ; it != end; ++it)
          {
            bool globalParamChangedLocal=true;
            bool timeChangedLocal=false;
            bool freqChangedLocal=false;
            if ((*it)->updateGlobalAndDependentParametersForStep(globalParamChangedLocal,timeChangedLocal,freqChangedLocal))
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
        // handle the special case of inductor instances contained inside mutual inductors.
        int inductorIndex=-1;
        DeviceInstance *device_instance=device_manager.getMutualInductorDeviceInstance(name,inductorIndex);
        int instance_found = (device_instance != 0);
        if (instance_found)
        {
          std::vector< double > inductorInductances = device_instance->getInductorInductances();
          inductorInductances[inductorIndex] = value ;
          device_instance->setInductorInductances( inductorInductances );
          device_instance->processParams ();
        }
        Parallel::AllReduce(comm, MPI_LOR, &instance_found, 1);

        if (!instance_found)
        {
          //if (DEBUG_DEVICE)
          if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
          {
            Report::DevelWarning() << "Xyce::Device::setParameter.  Unable to find parameter " << name;
          }
          else
          {
            Report::UserError() << "Xyce::Device::setParameter.  Unable to find parameter " << name;
          }
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
// Function      : setParameterRandomExpressionTerms2
// Purpose       : Update params that depend on random ops such as
//
// Special Notes : AGAUSS, GAUSS, AUNIF, UNIF, RAND and LIMIT.
//                 This version does all the params in one call.
//
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 3/3/2022
//-----------------------------------------------------------------------------
bool setParameterRandomExpressionTerms2(
  Parallel::Machine             comm,
  ArtificialParameterMap &      artificial_parameter_map,
  PassthroughParameterSet &     passthrough_parameter_map,
  UserDefinedParams &           globals,               ///< global variables
  DeviceMgr &                   device_manager,
  EntityVector &                dependent_entity_vector,   // dependentPtrVec_, if called from DeviceMgr::setParam
  const InstanceVector &        extern_device_vector,
  const std::vector<Xyce::Analysis::SweepParam> & SamplingParams,
  std::unordered_map<DeviceEntity *, std::vector<Depend> > & masterGlobalParamDependencies, // should be const
  bool                          override_original)
{
  bool bsuccess = true, success = true;

  GlobalParameterMap &          global_parameter_map = globals.paramMap;
  std::vector<Util::Expression> & expressionVec = globals.expressionVec;
  std::vector<std::string> & expNameVec = globals.expNameVec;

  std::set<int> uniqueGlobalIndices;
  std::set <DeviceEntity *> uniqueDeviceEntities;

  for(int isamp=0;isamp<SamplingParams.size();isamp++)
  {
    std::string name   = SamplingParams[isamp].setParamName;
    std::string opName = SamplingParams[isamp].opName;
    int opIndex        = SamplingParams[isamp].astOpIndex;
    int astType        = SamplingParams[isamp].astType;
    int globalIndex    = SamplingParams[isamp].globalIndex;
    double value       = SamplingParams[isamp].currentVal;
    GlobalParameterMap::iterator global_param_it = SamplingParams[isamp].gpIter;

    if (globalIndex >= 0)
    {
      Util::Expression &expression = expressionVec[globalIndex];

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

      uniqueGlobalIndices.insert(globalIndex);
    }
    else
    {
      // If not .param/.global_param, then search for the appropriate device param(s).
      DeviceEntity * device_entity = device_manager.getDeviceEntity(name);

      int entity_found = (device_entity != 0);
      if (entity_found)
      {
        bool found;
        std::string paramName = Util::paramNameFromFullParamName(name);
        found = device_entity->setParameterRandomExpressionTerms(paramName,opIndex,astType,value,override_original);

        if (found)
        {
          uniqueDeviceEntities.insert(device_entity);
        }

        entity_found = found;
      }

      Parallel::AllReduce(comm, MPI_LOR, &entity_found, 1);

      if (entity_found == 0)
      {
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
        {
          Report::DevelWarning() << "DeviceMgr.C:  setParameterRandomExpressionTerms.  Unable to find parameter " << name;
        }
        else
        {
          Report::UserError() << "Xyce::Device::setParameterRandomExpressionTerms.  Unable to find parameter " << name;
        }
      }
    }
  }


  if ( !(uniqueGlobalIndices.empty()) )
  {
    // ERK.  5/11/2023.
    // Note that the ONLY reason for running these "evaluateFunction" calls here is
    // to (if necessary) provide a correct value for a global parameter being printed on the .PRINT line
    // (and maybe the *res file as well).  It isn't needed for updating the device parameter 
    // values, below, since they share the same ASTs.  And, w.r.t. the .PRINT line, 
    // it is only needed if it is requested as a raw global parameter, rather than an expression.
    // So, for example, ".PRINT TRAN paramName" instead of ".PRINT TRAN {paramName}".  The first
    // case will request the value directly from the device manager, so the newValue, below is needed
    // for that to be correct.  The second case is evaluated as an expression, so no need to make such a request.
    // If the Op classes can be modified so that "globalParameterOps" are always converted to expression Ops, 
    // then this block of code can go.
    //
    // See also: the values in the *res file.
    std::set<int>::iterator iter = uniqueGlobalIndices.begin();
    std::set<int>::iterator end = uniqueGlobalIndices.end();
    for (;iter!=end;iter++)
    {
      int globalIndex = *iter;
      Util::Expression &expression = expressionVec[globalIndex];
      double newValue=0.0;
      expression.evaluateFunction(newValue);

      GlobalParameterMap::iterator global_param_it = global_parameter_map.find(expNameVec[globalIndex]); // this "find" should only have happened 1x, in the getRandomParams function.  That is what the gpIter was for.  Fix later.
      if (global_param_it != global_parameter_map.end())
      {
        (*global_param_it).second = newValue; 
      }
      else
      {
        Report::DevelFatal().in("Xyce::Device::setParameterRandomExpressionTerms2") // we should never get here.
          << "Cannot find global parameter";
      }
    }
  }

  if ( !(uniqueDeviceEntities.empty()) )
  {
    std::set<DeviceEntity *>::iterator iter = uniqueDeviceEntities.begin();
    std::set<DeviceEntity *>::iterator end = uniqueDeviceEntities.end();
    for (;iter!=end;iter++)
    {
      DeviceEntity * device_entity = *iter;
      device_entity->processParams (); 
      device_entity->processInstanceParams();
    }
  }

  // loop over the "master" list of global param dependencies 
  // and perform updates.  In this list, each entity is unique, ie only appears once.
  // The parameter vector for each entity is also unique.  So, 
  // within each unique entity, each param only appears, at most, once.
  //
  // The params in the parameterVec ONLY include parameters that 
  // have one or more .param/.global_param dependencies.  Any params 
  // that lack these dependencies have been excluded.
  {
  std::unordered_map<DeviceEntity *, std::vector<Depend> >::iterator mgpIter = 
              masterGlobalParamDependencies.begin();
  std::unordered_map<DeviceEntity *, std::vector<Depend> >::iterator mgpEnd = 
              masterGlobalParamDependencies.end();

  for ( ; mgpIter != mgpEnd; mgpIter++)
  {
    DeviceEntity * entityPtr = mgpIter->first;
    std::vector<Depend> & parameterVec = mgpIter->second;

    if (!(parameterVec.empty()))
    {
      bool globalParamChangedLocal=true;
      bool timeChangedLocal=false;
      bool freqChangedLocal=false;
      if (entityPtr->updateGlobalAndDependentParameters(
              globalParamChangedLocal,
              timeChangedLocal,
              freqChangedLocal,
              parameterVec))
      {
        entityPtr->processParams();
        entityPtr->processInstanceParams();
      }
    }
  }
  }

#if 0
  // ERK. This will need to be set up later.  Currently, the random operators 
  // cannot be applied thru the extern device interface.  This code below is just
  // copied from the original setParameter function, and won't work correctly
  // for this use case.  Leaving in place as a reminder.

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
  UserDefinedParams &   globals,
  const Util::Param &   param,
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> & expressionGroup
  )
{
  if (param.getType() == Util::EXPR)
  {
    globals.expressionVec.push_back(param.getValue<Util::Expression>());
    Util::Expression &expression = globals.expressionVec.back();

    globals.expNameVec.push_back(param.uTag());

    globals.deviceEntityDependVec.push_back(std::vector<entityDepend>(0));

    double val;
    expression.evaluateFunction(val);
    globals.paramMap[param.uTag()] = val;

    // the group may get updatd later 
    expression.setGroup(expressionGroup);
  }
  else
  {
    globals.paramMap[param.uTag()] = param.getImmutableValue<double>();
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addGlobalParameters()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/12/2023
//-----------------------------------------------------------------------------
void addGlobalParameters(
  UserDefinedParams &   globals,
  const Util::UParamList & ioGlobalParams,
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> & expressionGroup
  )
{
  if (!(ioGlobalParams.empty()))
  {
    int size = ioGlobalParams.size();

    if (globals.expressionVec.empty()) { globals.expressionVec.reserve(size); }
    if (globals.expNameVec.empty()) { globals.expNameVec.reserve(size); }
    if (globals.deviceEntityDependVec.empty()) { globals.deviceEntityDependVec.reserve(size); }

    Util::UParamList::const_iterator begin = ioGlobalParams.begin();
    Util::UParamList::const_iterator end= ioGlobalParams.end();
    for (; begin != end; ++begin)
    {
      const Util::Param & param = *begin;

      if (param.getType() == Util::EXPR)
      {
        globals.expressionVec.push_back(param.getValue<Util::Expression>());
        Util::Expression &expression = globals.expressionVec.back();

        globals.expNameVec.push_back(param.uTag());

        globals.deviceEntityDependVec.push_back(std::vector<entityDepend>(0));

        double val;
        expression.evaluateFunction(val);
        globals.paramMap[param.uTag()] = val;

        // the group may get updated later 
        expression.setGroup(expressionGroup);
      }
      else
      {
        globals.paramMap[param.uTag()] = param.getImmutableValue<double>();
      }
    }
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

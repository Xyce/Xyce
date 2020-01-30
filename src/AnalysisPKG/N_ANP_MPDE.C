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
// Purpose       : MPDE analysis functions.
// Special Notes :
// Creator       : Todd Coffey, 1414, Ting Mei, 1437
// Creation Date : 07/23/08
//
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_MPDE.h>
#include <N_DEV_DeviceMgr.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_LOA_Loader.h>
#include <N_MPDE_Manager.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_NoTimeIntegration.h>
#include <N_UTL_Factory.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : MPDE::MPDE( AnalysisManager * )
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
MPDE::MPDE(
  AnalysisManager &                     analysis_manager,
  Linear::System &                      linear_system,
  Nonlinear::Manager &                  nonlinear_manager,
  Loader::Loader &                      loader,
  Device::DeviceMgr &                   device_manager,
  Linear::Builder &                     builder,
  Topo::Topology &                      topology,
  IO::InitialConditionsManager &        initial_conditions_manager,
  IO::RestartMgr &                      restart_manager)
  : AnalysisBase(analysis_manager, "MPDE"),
    StepEventListener(&analysis_manager),
    analysisManager_(analysis_manager),
    loader_(loader),
    linearSystem_(linear_system),
    nonlinearManager_(nonlinear_manager),
    topology_(topology),
    mpdeManager_(new N_MPDE_Manager(analysisManager_, loader, device_manager, builder, topology, initial_conditions_manager, restart_manager, analysisManager_.getCommandLine()))
{}

void MPDE::notify(const StepEvent &event)
{
  if (event.state_ == StepEvent::STEP_STARTED)
  {
    AnalysisBase::resetForStepAnalysis();
  }
}

//-----------------------------------------------------------------------------
// Function      : MPDE::~MPDE
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
MPDE::~MPDE()
{
  delete mpdeManager_;
}

const TimeIntg::TIAParams &
MPDE::getTIAParams() const
{
  return mpdeManager_->getTIAParams();
}

TimeIntg::TIAParams &
MPDE::getTIAParams()
{
  return mpdeManager_->getTIAParams();
}

//-----------------------------------------------------------------------------
// Function      : MPDE::getDCOPFlag()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/24/2014
//-----------------------------------------------------------------------------
bool MPDE::getDCOPFlag() const
{
  return (analysisManager_.getCurrentAnalysis() )->getDCOPFlag();
}

//-----------------------------------------------------------------------------
// Function      : MPDE::run()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::doRun()
{
  // mpdeManager_->registerApplicationLoader(&loader_);
  mpdeManager_->run(linearSystem_, nonlinearManager_, topology_);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : MPDE::init()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::doInit()
{
//  analysisManager_.getMPDEManager()->initializeAll();

  return false;
}

//-----------------------------------------------------------------------------
// Function      : MPDE::loopProcess()
// Purpose       : Conduct the time stepping loop.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::doLoopProcess()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : MPDE::processSuccessfulDCOP()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::processSuccessfulDCOP()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : MPDE::processSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::doProcessSuccessfulStep()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : MPDE::processFailedStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::doProcessFailedStep()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : MPDE::processFailedDCOP
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::processFailedDCOP()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : MPDE::finish
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::doFinish()
{
  return false;
}

bool MPDE::doHandlePredictor()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : MPDE::finalVerboseOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::finalVerboseOutput()
{
  return false;
}

namespace {

typedef Util::Factory<AnalysisBase, MPDE>  MPDEFactoryBase;

//-----------------------------------------------------------------------------
// Class         : MPDEFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jan 29 12:53:02 2015
//-----------------------------------------------------------------------------
///
/// Factory for parsing MPDE parameters from the netlist and creating MPDE analysis.
///
class MPDEFactory : public MPDEFactoryBase
{
public:
  //-----------------------------------------------------------------------------
  // Function      : MPDEFactory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:54:09 2015
  //-----------------------------------------------------------------------------
  ///
  /// Constructs the MPDE analysis factory
  ///
  /// @invariant Stores the results of parsing, so if more than one of the analysis and
  /// associated lines are parsed, the second options simply overwrite the previously parsed
  /// values.
  ///
  /// @invariant The existence of the parameters specified in the constructor cannot
  /// change.
  ///
  /// @param analysis_manager 
  /// @param linear_system 
  /// @param nonlinear_manager 
  /// @param device_manager 
  /// @param builder 
  /// @param topology 
  ///
  MPDEFactory(
    Analysis::AnalysisManager &         analysis_manager,
    Linear::System &                    linear_system,
    Nonlinear::Manager &                nonlinear_manager,
    Loader::Loader &                    loader,
    Device::DeviceMgr &                 device_manager,
    Linear::Builder &                   builder,
    Topo::Topology &                    topology,
    IO::InitialConditionsManager &      initial_conditions_manager,
    IO::RestartMgr &                    restart_manager)
    : MPDEFactoryBase(),
      analysisManager_(analysis_manager),
      linearSystem_(linear_system),
      nonlinearManager_(nonlinear_manager),
      loader_(loader),
      deviceManager_(device_manager),
      builder_(builder),
      topology_(topology),
      initialConditionsManager_(initial_conditions_manager),
      restartManager_(restart_manager)
  {}

  virtual ~MPDEFactory()
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
  /// Create a new MPDE analysis and applies the analysis and time integrator option blocks.
  ///
  /// @return new MPDE analysis object
  ///
  MPDE *create() const
  {
    analysisManager_.setAnalysisMode(ANP_MODE_MPDE);

    MPDE *mpde = new MPDE(analysisManager_, linearSystem_, nonlinearManager_, loader_, deviceManager_, builder_, topology_, initialConditionsManager_, restartManager_);

    mpde->getMPDEManager().setMPDEAnalysisParams(mpdeAnalysisOptionBlock_);
    mpde->getMPDEManager().setMPDEOptions(mpdeIntOptionBlock_);
    mpde->getMPDEManager().setTransientOptions(mpdeTimeIntegratorOptionBlock_);
    mpde->getMPDEManager().setLinSolOptions(linSolOptionBlock_);

    return mpde;
  }

  //-----------------------------------------------------------------------------
  // Function      : setMPDEAnalysisOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 13:00:14 2015
  //-----------------------------------------------------------------------------
  ///
  /// Saves the analysis parsed options block in the factory.
  ///
  /// @invariant Overwrites any previously specified analysis option block.
  ///
  /// @param option_block parsed option block
  ///
  void setMPDEAnalysisOptionBlock(const Util::OptionBlock &option_block)
  {
    mpdeAnalysisOptionBlock_ = option_block;
  }

  //-----------------------------------------------------------------------------
  // Function      : setMPDETimeIntegratorOptionBlock
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
  bool setMPDETimeIntegratorOptionBlock(const Util::OptionBlock &option_block)
  {
    mpdeTimeIntegratorOptionBlock_ = option_block;

    return true;
  }

  //-----------------------------------------------------------------------------
  // Function      : setMPDEIntOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 13:01:27 2015
  //-----------------------------------------------------------------------------
  ///
  /// Saves the MPDEINT parsed option block.
  ///
  /// @invariant Overwrites any previously specified time integrator option block.
  ///
  /// @param option_block parsed option block
  ///
  bool setMPDEIntOptionBlock(const Util::OptionBlock &option_block)
  {
    mpdeIntOptionBlock_ = option_block;

    return true;
  }

  //-----------------------------------------------------------------------------
  // Function      : setLinSolOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 13:01:27 2015
  //-----------------------------------------------------------------------------
  ///
  /// Saves the LINSOL parsed option block.
  ///
  /// @invariant Overwrites any previously specified linear solver option block.
  ///
  /// @param option_block parsed option block
  ///

  bool setLinSolOptionBlock(const Util::OptionBlock &option_block)
  {
    linSolOptionBlock_ = option_block;

    return true;
  }

public:
  AnalysisManager &                     analysisManager_;
  Linear::System &                      linearSystem_;
  Nonlinear::Manager &                  nonlinearManager_;
  Loader::Loader &                      loader_;
  Device::DeviceMgr &                   deviceManager_;
  Linear::Builder &                     builder_;
  Topo::Topology &                      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;
  IO::RestartMgr &                      restartManager_;

private:
  Util::OptionBlock     mpdeAnalysisOptionBlock_;
  Util::OptionBlock     mpdeIntOptionBlock_;
  Util::OptionBlock     mpdeTimeIntegratorOptionBlock_;
  Util::OptionBlock     timeIntegratorOptionBlock_;
  Util::OptionBlock     linSolOptionBlock_;
};

// .MPDE
struct MPDEAnalysisReg : public IO::PkgOptionsReg
{
  MPDEAnalysisReg(
    MPDEFactory &   factory )
    : factory_(factory)
  {}

  bool operator()(const Util::OptionBlock &option_block)
  {
    factory_.setMPDEAnalysisOptionBlock(option_block);
    factory_.deviceManager_.setBlockAnalysisFlag(true);

    factory_.analysisManager_.addAnalysis(&factory_);

    return true;
  }

  MPDEFactory &               factory_;
};

//-----------------------------------------------------------------------------
// Function      : extractMPDEData
// Purpose       : Extract the parameters from a netlist .DC line held in
//                 parsed_line.
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool
extractMPDEData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("MPDE", Util::OptionBlock::NO_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

  int numFields = parsed_line.size();

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 3 || numFields > 6 )
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << ".MPDE line has an unexpected number of fields";
  }

  int linePosition = 1;   // Start of parameters on .param line.
  int endPosition = numFields - 1;

  Util::Param parameter("", "");

  // TSTEP and TSTOP are required, get them now.
  parameter.setTag( "TSTEP" );
  parameter.setVal( parsed_line[linePosition].string_ );
  option_block.addParam( parameter );
  ++linePosition;     // Advance to next parameter.

  parameter.setTag( "TSTOP" );
  parameter.setVal( parsed_line[linePosition].string_ );
  option_block.addParam( parameter );
  ++linePosition;     // Advance to next parameter.

  // Next check last field to see if it is UIC or NOOP.
  parameter.setTag( parsed_line[endPosition].string_ );
  if ( parameter.uTag() == "NOOP" || parameter.uTag() == "UIC" )
  {
    parameter.setVal( "1" );
    option_block.addParam( parameter );
    --endPosition;
  }
  else if ( numFields == 6)
  {
    Report::UserError0().at(netlist_filename, parsed_line[endPosition].lineNumber_) << "expected NOOP/UIC field on .MPDE line but found" << parameter.usVal();
  }

  if ( linePosition <= endPosition )
  {
    parameter.setTag( "TSTART" );
    parameter.setVal( parsed_line[linePosition].string_ );
    option_block.addParam( parameter );
    ++linePosition;     // Advance to next parameter.
  }

  if ( linePosition <= endPosition )
  {
    parameter.setTag( "DTMAX" );
    parameter.setVal( parsed_line[linePosition].string_ );
    option_block.addParam( parameter );
    ++linePosition;     // Advance to next parameter.
  }

  circuit_block.addOptions(option_block);

  return true;
}

void
populateMetadata(
  IO::PkgOptionsMgr &   options_manager)
{
  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("TIMEINT-MPDE");

    parameters.insert(Util::ParamMap::value_type("METHOD", Util::Param("METHOD", 1)));
    if (DEBUG_ANALYSIS)
      parameters.insert(Util::ParamMap::value_type("CONSTSTEP", Util::Param("CONSTSTEP", 0)));
    parameters.insert(Util::ParamMap::value_type("USEDEVICEMAX", Util::Param("USEDEVICEMAX", 1)));
    parameters.insert(Util::ParamMap::value_type("RELTOL", Util::Param("RELTOL", 1.0E-2)));
    parameters.insert(Util::ParamMap::value_type("ABSTOL", Util::Param("ABSTOL", 1.0E-6)));
    parameters.insert(Util::ParamMap::value_type("RESTARTSTEPSCALE", Util::Param("RESTARTSTEPSCALE", .005)));
    parameters.insert(Util::ParamMap::value_type("NLNEARCONV", Util::Param("NLNEARCONV", 0)));
    parameters.insert(Util::ParamMap::value_type("NLSMALLUPDATE", Util::Param("NLSMALLUPDATE", 1)));
    parameters.insert(Util::ParamMap::value_type("DOUBLEDCOPSTEP", Util::Param("DOUBLEDCOPSTEP", 0)));
    parameters.insert(Util::ParamMap::value_type("FIRSTDCOPSTEP", Util::Param("FIRSTDCOPSTEP", 0)));
    parameters.insert(Util::ParamMap::value_type("LASTDCOPSTEP", Util::Param("LASTDCOPSTEP", 1)));
    parameters.insert(Util::ParamMap::value_type("RESETTRANNLS", Util::Param("RESETTRANNLS", 1)));
    parameters.insert(Util::ParamMap::value_type("BPENABLE", Util::Param("BPENABLE", 1)));
    parameters.insert(Util::ParamMap::value_type("EXITTIME", Util::Param("EXITTIME", 0.0)));
    parameters.insert(Util::ParamMap::value_type("EXITSTEP", Util::Param("EXITSTEP", 0)));
    parameters.insert(Util::ParamMap::value_type("ERROPTION", Util::Param("ERROPTION", 0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGLEVEL", Util::Param("DEBUGLEVEL", 0)));
    parameters.insert(Util::ParamMap::value_type("JACLIMITFLAG", Util::Param("JACLIMITFLAG", 0)));
    parameters.insert(Util::ParamMap::value_type("JACLIMIT", Util::Param("JACLIMIT", 1.0e17)));
    parameters.insert(Util::ParamMap::value_type("DAESTATEDERIV", Util::Param("DAESTATEDERIV", 0)));
    parameters.insert(Util::ParamMap::value_type("TESTFIRSTSTEP", Util::Param("TESTFIRSTSTEP", 0)));
    parameters.insert(Util::ParamMap::value_type("DTMIN", Util::Param("DTMIN", 0.0)));
    parameters.insert(Util::ParamMap::value_type("MAXORD", Util::Param("MAXORD", 5)));
    parameters.insert(Util::ParamMap::value_type("MINORD", Util::Param("MINORD", 1)));
    parameters.insert(Util::ParamMap::value_type("OUTPUTINTERPMPDE", Util::Param("OUTPUTINTERPMPDE", 1)));
    parameters.insert(Util::ParamMap::value_type("INTERPOUTPUT", Util::Param("INTERPOUTPUT", 1)));
    parameters.insert(Util::ParamMap::value_type("CONDTEST", Util::Param("CONDTEST", 0)));
    parameters.insert(Util::ParamMap::value_type("CONDTESTDEVICENAME", Util::Param("CONDTESTDEVICENAME", "dev_name")));
    parameters.insert(Util::ParamMap::value_type("ISOCONDTEST", Util::Param("ISOCONDTEST", 0)));
    parameters.insert(Util::ParamMap::value_type("ISOCONDTESTDEVICENAME", Util::Param("ISOCONDTESTDEVICENAME", "dev_name")));
    parameters.insert(Util::ParamMap::value_type("MINTIMESTEPSBP", Util::Param("MINTIMESTEPSBP", 10)));
    parameters.insert(Util::ParamMap::value_type("NLMIN", Util::Param("NLMIN", 3)));
    parameters.insert(Util::ParamMap::value_type("NLMAX", Util::Param("NLMAX", 8)));
    parameters.insert(Util::ParamMap::value_type("DELMAX", Util::Param("DELMAX", 1.0e+99)));
    parameters.insert(Util::ParamMap::value_type("TIMESTEPSREVERSAL", Util::Param("TIMESTEPSREVERSAL", false)));
  }

  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("MPDEINT");

    parameters.insert(Util::ParamMap::value_type("AUTON2", Util::Param("AUTON2", false)));
    parameters.insert(Util::ParamMap::value_type("AUTON2MAX", Util::Param("AUTON2MAX", 100)));
    parameters.insert(Util::ParamMap::value_type("DCOPEXIT", Util::Param("DCOPEXIT", 0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGLEVEL", Util::Param("DEBUGLEVEL", 0)));
    parameters.insert(Util::ParamMap::value_type("DIFF", Util::Param("DIFF", 0)));
    parameters.insert(Util::ParamMap::value_type("DIFFORDER", Util::Param("DIFFORDER", 1)));
    parameters.insert(Util::ParamMap::value_type("EXITSAWTOOTHSTEP", Util::Param("EXITSAWTOOTHSTEP", -1)));
    parameters.insert(Util::ParamMap::value_type("FREQDOMAIN", Util::Param("FREQDOMAIN", 0)));
    parameters.insert(Util::ParamMap::value_type("IC", Util::Param("IC", 0)));
    parameters.insert(Util::ParamMap::value_type("ICEXIT", Util::Param("ICEXIT", 0)));
    parameters.insert(Util::ParamMap::value_type("ICPER", Util::Param("ICPER", 10)));
    parameters.insert(Util::ParamMap::value_type("N2", Util::Param("N2", 10)));
    parameters.insert(Util::ParamMap::value_type("NONLTESTEPS", Util::Param("NONLTESTEPS", 10)));
    parameters.insert(Util::ParamMap::value_type("OSCOUT", Util::Param("OSCOUT", "")));
    parameters.insert(Util::ParamMap::value_type("OSCSRC", Util::Param("OSCSRC", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("PHASE", Util::Param("PHASE", 0)));
    parameters.insert(Util::ParamMap::value_type("PHASECOEFF", Util::Param("PHASECOEFF", 0)));
    parameters.insert(Util::ParamMap::value_type("SAVEICDATA", Util::Param("SAVEICDATA", false)));
    parameters.insert(Util::ParamMap::value_type("STARTUPPERIODS", Util::Param("STARTUPPERIODS", 0)));
    parameters.insert(Util::ParamMap::value_type("T2", Util::Param("T2", 0.0)));
    parameters.insert(Util::ParamMap::value_type("TEST", Util::Param("TEST", 0)));
    parameters.insert(Util::ParamMap::value_type("WAMPDE", Util::Param("WAMPDE", 0)));
  }
}

} // namespace <unnamed>

bool
registerMPDEFactory(
  FactoryBlock &        factory_block)
{
  MPDEFactory *factory = new MPDEFactory(factory_block.analysisManager_, factory_block.linearSystem_, factory_block.nonlinearManager_, factory_block.loader_, factory_block.deviceManager_, factory_block.builder_, factory_block.topology_, factory_block.initialConditionsManager_, factory_block.restartManager_);

  addAnalysisFactory(factory_block, factory);

  populateMetadata(factory_block.optionsManager_);

  factory_block.optionsManager_.addCommandParser(".MPDE", extractMPDEData);

  factory_block.optionsManager_.addCommandProcessor("MPDE", new MPDEAnalysisReg(*factory));

  factory_block.optionsManager_.addOptionsProcessor("MPDEINT", IO::createRegistrationOptions(*factory, &MPDEFactory::setMPDEIntOptionBlock));
  factory_block.optionsManager_.addOptionsProcessor("TIMEINT-MPDE", IO::createRegistrationOptions(*factory, &MPDEFactory::setMPDETimeIntegratorOptionBlock));

  factory_block.optionsManager_.addOptionsProcessor("LINSOL", IO::createRegistrationOptions(*factory, &MPDEFactory::setLinSolOptionBlock));

  return true;
}

} // namespace Analysis
} // namespace Xyce

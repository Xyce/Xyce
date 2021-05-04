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
// Purpose       : .STEP Sweep class analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_SweepParam.h>
#include <N_ANP_SweepParamFreeFunctions.h>
#include <N_ANP_Step.h>
#include <N_ANP_StepEvent.h>
#include <N_ERH_Message.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_DataStore.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_Factory.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Timer.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : Step::setAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 6/22/10
//-----------------------------------------------------------------------------
bool Step::setAnalysisParams(const Util::OptionBlock & paramsBlock)
{
  if (isDataSpecified(paramsBlock))
  {
    // This handle the case of having multiple .STEP lines in the netlist, of
    // which only some might use DATA=<tableName>
    dataSpecification_ = true;
  }
  stepSweepVector_.push_back(parseSweepParams(paramsBlock.begin(), paramsBlock.end()));
  outputManagerAdapter_.setStepSweepVector(stepSweepVector_);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Step::setDataStatements
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/5/18
//-----------------------------------------------------------------------------
bool Step::setDataStatements(const Util::OptionBlock & paramsBlock)
{
  return processDataStatements(paramsBlock, dataNamesMap_, dataTablesMap_);
}

//-----------------------------------------------------------------------------
// Function      : Step::convertDataToSweepParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/5/18
//-----------------------------------------------------------------------------
bool  Step::convertDataToSweepParams()
{
  return convertData( stepSweepVector_, dataNamesMap_, dataTablesMap_);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
const TimeIntg::TIAParams & Step::getTIAParams() const
{
  return childAnalysis_.getTIAParams();
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
TimeIntg::TIAParams & Step::getTIAParams()
{
  return childAnalysis_.getTIAParams();
}

//-----------------------------------------------------------------------------
// Function      : Step::getDCOPFlag()
// Purpose       :
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/24/2014
//-----------------------------------------------------------------------------
bool Step::getDCOPFlag() const
{
  return childAnalysis_.getDCOPFlag();
}

//-----------------------------------------------------------------------------
// Function      : Step::finalExpressionBasedSetup()
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/4/2021
//-----------------------------------------------------------------------------
void Step::finalExpressionBasedSetup()
{
  childAnalysis_.finalExpressionBasedSetup();
}

//-----------------------------------------------------------------------------
// Function      : Step::run()
// Purpose       : This is the main controlling loop for Step analysis.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/04/00
//-----------------------------------------------------------------------------
bool Step::doRun()
{
  return doInit() && doLoopProcess() && doFinish();
}

//-----------------------------------------------------------------------------
// Function      : Step::init()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool Step::doInit()
{
  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    Xyce::dout() << std::endl << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "Step::init" << std::endl;
  }

  // check if the "DATA" specification was used.  If so, create a new vector of 
  // SweepParams, in the "TABLE" style.
  if (dataSpecification_ && !convertDataToSweepParams())
  {
    Report::UserFatal() << "Invalid data=<name> parameter on .STEP line.";
    return false;
  }

  stepLoopSize_ = setupSweepLoop(
      analysisManager_.getComm(), 
      loader_, 
      stepSweepVector_.begin(), 
      stepSweepVector_.end());

  Util::publish<StepEvent>(
      analysisManager_, 
      StepEvent(StepEvent::INITIALIZE, stepSweepVector_, stepLoopSize_));

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Step::loopProcess()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool Step::doLoopProcess()
{
  bool integration_status = true;

  for (int i = 0; i < stepLoopSize_; ++i)
  {
    // Tell the manager if any of our sweeps are being reset in this loop iteration.
    bool reset = updateSweepParams(loader_, i, stepSweepVector_.begin(), stepSweepVector_.end(), false);

    analysisManager_.setSweepSourceResetFlag(reset);

    outputManagerAdapter_.setStepSweepVector(stepSweepVector_);

    if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
    {
      for (SweepVector::const_iterator it = stepSweepVector_.begin(), end = stepSweepVector_.end(); it != end; ++it)
      {
        Xyce::dout() << "Step Analysis # " << i<<"\t";
        Xyce::dout() << (*it);
      }
    }

    StepEvent step_event(StepEvent::STEP_STARTED, stepSweepVector_, i);
    Util::publish<StepEvent>(analysisManager_, step_event);

    // solve the loop.
    integration_status = childAnalysis_.run();

    step_event.state_ = StepEvent::STEP_COMPLETED;
    step_event.finalSimTime_ = getTIAParams().finalTime;
    Util::publish<StepEvent>(analysisManager_, step_event);
  }

  return integration_status;
}

//-----------------------------------------------------------------------------
// Function      : Step::processSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool Step::doProcessSuccessfulStep()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Step::processFailedStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool Step::doProcessFailedStep()
{
  return true;
}


//-----------------------------------------------------------------------------
// Function      : Step::doFinish()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool Step::doFinish()
{
  Util::publish<StepEvent>(analysisManager_, StepEvent(StepEvent::FINISH, stepSweepVector_, stepLoopSize_));

  return true;
}

namespace {

//-----------------------------------------------------------------------------
// Class         : StepFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jan 29 12:53:02 2015
//-----------------------------------------------------------------------------
///
/// Factory for parsing Step parameters from the netlist and creating Step analysis.
///
class StepFactory : public Util::Factory<AnalysisBase, Step>
{
public:
  //-----------------------------------------------------------------------------
  // Function      : StepFactory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:54:09 2015
  //-----------------------------------------------------------------------------
  ///
  /// Constructs the Step analysis factory
  ///
  /// @invariant Stores the results of parsing.  Multiple Step analysis options may be
  /// applied and each generates and additional step.
  ///
  /// @invariant The existence of the parameters specified in the constructor cannot
  /// change.
  ///
  /// @param analysis_manager 
  /// @param linear_system 
  /// @param nonlinear_manager 
  ///
  StepFactory(
    Analysis::AnalysisManager & analysis_manager,
    Linear::System &            linear_system,
    Nonlinear::Manager &        nonlinear_manager,
    Loader::Loader &            loader)
    : Util::Factory<AnalysisBase, Step>(),
      analysisManager_(analysis_manager),
      linearSystem_(linear_system),
      nonlinearManager_(nonlinear_manager),
      loader_(loader)
  {}

  virtual ~StepFactory()
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
  /// Create a new Step analysis and applies the analysis and time integrator option blocks.
  ///
  /// @return new Step analysis object
  ///
  Step *create() const
  {
    Step *step = new Step(analysisManager_, loader_, analysisManager_.getAnalysisObject());

    for (std::vector<Util::OptionBlock>::const_iterator it = stepSweepAnalysisOptionBlock_.begin(), end = stepSweepAnalysisOptionBlock_.end(); it != end; ++it)
    {
      step->setAnalysisParams(*it);
    }

    for (std::vector<Util::OptionBlock>::const_iterator it = dataOptionBlockVec_.begin(), end = dataOptionBlockVec_.end(); it != end; ++it)
    {
      step->setDataStatements(*it);
    }


    return step;
  }

  //-----------------------------------------------------------------------------
  // Function      : setStepAnalysisOptionBlock
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
  void setStepAnalysisOptionBlock(const Util::OptionBlock &option_block)
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
  }

  //-----------------------------------------------------------------------------
  bool setDotDataBlock(const Util::OptionBlock &option_block)
  {
    dataOptionBlockVec_.push_back(option_block);
    return true;
  }

public:
  AnalysisManager &             analysisManager_;
  Linear::System &              linearSystem_;
  Nonlinear::Manager &          nonlinearManager_;
  Loader::Loader &              loader_;

private:
  std::vector<Util::OptionBlock>        stepSweepAnalysisOptionBlock_;
  std::vector<Util::OptionBlock>        dataOptionBlockVec_;
};

// .STEP
struct StepAnalysisReg : public IO::PkgOptionsReg
{
  StepAnalysisReg(
    StepFactory &             factory)
    : factory_(factory)
  {}

  bool operator()(const Util::OptionBlock &option_block)
  {
    factory_.setStepAnalysisOptionBlock(option_block);

    factory_.analysisManager_.addAnalysis(&factory_);

    return true;
  }

  StepFactory &               factory_;
};

//-----------------------------------------------------------------------------
// Function      : extractSTEPData
// Purpose       : Extract the parameters from a netlist .STEP line held in
//                 parsed_line.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 10/30/2003
//-----------------------------------------------------------------------------
bool extractSTEPData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("STEP", Util::OptionBlock::NO_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

  int numFields = parsed_line.size();

  if (numFields < 4)
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << ".STEP line has an unexpected number of fields";
    return false;
  }

  // First check if the type has been explicitly set.
  // If not, set it to the default, LIN.
  int pos1=1;
  int tablePos=1;
  int dataPos=1;
  bool typeExplicitSetLinDecOct = false;
  bool typeExplicitSetList = false;
  std::string type("LIN");

  // check for "DATA" first, as data set names could be "TABLE".
  bool dataFound=false;
  while ( pos1 < numFields )
  {
    ExtendedString stringVal ( parsed_line[pos1].string_ );
    stringVal.toUpper ();

    if (stringVal == "DATA")
    {
      type = stringVal;
      typeExplicitSetList = true;
      dataPos=pos1;
      dataFound=true;
      break;
    }
    ++pos1;
  }

  // check for everything else: LIN, DEC, OCT, LIST, TABLE
  if (!dataFound)
  {
    pos1=1;
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
      else if (stringVal == "LIST" || stringVal == "TABLE")
      {
        typeExplicitSetList = true;
        type = stringVal;

        if (stringVal == "TABLE") { tablePos = pos1; }
      }
      ++pos1;
    }
  }

  // Check that the minimum required number of fields are on the line.
  int offset = 1;
  if (typeExplicitSetLinDecOct)
  {
    offset = 2;
  }

  if (dataFound)
  {
    if (numFields != 4)
    {
      Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
        << ".STEP line not formatted correctly.  numFields = " << numFields;
      return false;
    }
  }
  else
  {
    if (!typeExplicitSetList)// if this is a list, number of fields is arbitrary.
    {
      if ( (numFields-offset)%4 != 0 )
      {
        Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
          << ".STEP line not formatted correctly.";
        return false;
      }
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

      parameter.setTag( "STEP" );
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

      parameter.setTag( "NUMSTEPS" );
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

      parameter.setTag( "NUMSTEPS" );
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
  // experiment with this for now.  Change later?
  else if (type=="TABLE")
  {
    int linePosition=1; // if the parameter is first, and the word TABLE is second
    if (tablePos==1) { linePosition=2; } // if the word TABLE is first, and the parameter name is second

    parameter.setTag( "PARAM" );
    parameter.setVal(std::string(ExtendedString(parsed_line[linePosition].string_).toUpper()));
    option_block.addParam( parameter );

    linePosition=3;
    while (linePosition<numFields)
    {
      parameter.setTag( "VAL" );
      parameter.setVal( parsed_line[linePosition].string_ );
      option_block.addParam( parameter );
      ++linePosition;
    }
  }
  else if (type=="DATA")
  {
    parameter.setTag( "DATASET" );
    parameter.setVal( parsed_line[ dataPos+2 ].string_ );
    option_block.addParam( parameter );
  }
  else
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << ".STEP line contains an unrecognized type";
    return false;
  }

  circuit_block.addOptions(option_block);

  return true;
}

} // namespace <unnamed>


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
bool registerStepFactory(FactoryBlock & factory_block)
{
  StepFactory *factory = new StepFactory(
      factory_block.analysisManager_, 
      factory_block.linearSystem_, 
      factory_block.nonlinearManager_, 
      factory_block.loader_);

  addAnalysisFactory(factory_block, factory);

  factory_block.optionsManager_.addCommandProcessor("DATA", 
      IO::createRegistrationOptions(*factory, &StepFactory::setDotDataBlock) );

  factory_block.optionsManager_.addCommandParser(".STEP", extractSTEPData);
  factory_block.optionsManager_.addCommandProcessor("STEP", new StepAnalysisReg(*factory));

  return true;
}

} // namespace Analysis
} // namespace Xyce

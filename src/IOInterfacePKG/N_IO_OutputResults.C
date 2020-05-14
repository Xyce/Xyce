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
// Creator        : Dave Baur
//
// Creation Date  :
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_SweepParam.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ERH_Message.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_OutputResults.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>
#include <N_TIA_DataStore.h>
#include <N_UTL_ExpressionData.h>
#include <N_UTL_Factory.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : OutputResults::OutputResults
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
OutputResults::OutputResults(
  Parallel::Machine             comm,
  Analysis::AnalysisManager &   analysis_manager,
  OutputMgr &                   output_manager)
  : StepEventListener(&analysis_manager),
    comm_(comm),
    analysisManager_(analysis_manager),
    outputManager_(output_manager),
    os_(0),
    noIndexResult_(false)
{}

//-----------------------------------------------------------------------------
// Function      : OutputResults::~OutputResults
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
OutputResults::~OutputResults()
{
  delete os_;

  std::vector<Util::ExpressionData *>::iterator iter = resultVector_.begin();
  std::vector<Util::ExpressionData *>::iterator end = resultVector_.end();
  for (; iter!=end; ++iter)
  {
    if ( (*iter) )
    {
      delete (*iter);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputResults::addResultParams
// Purpose       : Sets the RESULT calculation parameters.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/29/04
//-----------------------------------------------------------------------------
bool OutputResults::addResultParams(
  const Teuchos::RCP<Xyce::Util::baseExpressionGroup> & exprGroup_,
  const Util::OptionBlock &     option_block)
{
// first check to see that there is only 1 PARAM set.  They need to be in
// separate lines for this to work.
  int countPar = 0;
  for (Util::ParamList::const_iterator it = option_block.begin(), 
      end = option_block.end(); it != end; ++it)
  {
    if ((*it).uTag() == "EXPRESSION")
    {
      ++countPar;
    }
  }

  if (countPar > 1)
  {
    Report::UserFatal0() 
      << "Only one expression per .RESULT command.  Each parameter needs its own .RESULT line.";
  }

  Util::ExpressionData * expDataPtr;

  for (Util::ParamList::const_iterator it = option_block.begin(), 
      end = option_block.end(); it != end; ++it)
  {
    const Util::Param &expParam = *it;

    if (!expParam.hasExpressionValue())
    {
      Report::DevelFatal0() << "Parameter must be an expression in .RESULT command";
    }
    else
    {
      expDataPtr = new Util::ExpressionData(exprGroup_,expParam.stringValue());
      resultVector_.push_back(expDataPtr);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputResults::setup
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void OutputResults::setup(
  Parallel::Machine             comm,
  OutputMgr &                   output_manager)
{
  for (ResultVector::const_iterator it = resultVector_.begin(), 
      end = resultVector_.end(); it != end; ++it)
  {
    (*it)->setup(comm,
                 output_manager.getOpBuilderManager(),
                 output_manager.getMainContextFunctionMap(),
                 output_manager.getMainContextParamMap(),
                 output_manager.getMainContextGlobalParamMap());
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputResults::outputRESULT
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/30/04
//-----------------------------------------------------------------------------
void OutputResults::output(
  Parallel::Machine             comm,
  double                        circuit_time,
  double                        circuit_dt,
  const Analysis::SweepVector & step_sweep_vector,
  int                           step_loop_number,
  const Linear::Vector &        solution_vector,
  const Linear::Vector &        state_vector,
  const Linear::Vector &        store_vector,
  const Linear::Vector &        lead_current_vector,
  const Linear::Vector &        junction_voltage_vector)
{
  std::string delim = " ";
  int width = 20;
  int width2 = 10;  // used to control spacing between STEP index and first sweep column
  int precision = 8;

  if (Parallel::rank(comm) == 0)
  {
    if (!os_)
    {
      std::string resultfilename = analysisManager_.getNetlistFilename() + ".res";

      os_ = new std::ofstream(resultfilename.c_str());

      os_->setf(std::ios::scientific);
      os_->precision(precision);

      if (!noIndexResult_)
      {
        (*os_) << "STEP";
      }

      for (Analysis::SweepVector::const_iterator it = step_sweep_vector.begin(), 
          end = step_sweep_vector.end(); it != end; ++it)
      {
        (*os_) << delim << std::setw(width) << (*it).name;
      }

      for (ResultVector::const_iterator it = resultVector_.begin(), 
          end = resultVector_.end(); it != end; ++it)
      {
        const Util::ExpressionData *expDataPtr = (*it);
        (*os_) << delim << std::setw(width) << expDataPtr->getExpression();
      }
      (*os_) << std::endl;
    }
  }

  if (Parallel::rank(comm) == 0)
  {
    os_->setf(std::ios::left, std::ios::adjustfield);
    if (!noIndexResult_)
    {
      (*os_) << std::setw(width2) << step_loop_number;
    }

    for (Analysis::SweepVector::const_iterator it = step_sweep_vector.begin(), 
        end = step_sweep_vector.end(); it != end; ++it)
    {
      (*os_) << delim << std::setw(width) << (*it).currentVal;
    }
  }

  for (ResultVector::const_iterator it = resultVector_.begin(), 
      end = resultVector_.end(); it != end; ++it)
  {
    Util::ExpressionData &expression_data = *(*it);

    double result = 0.0;
    if (expression_data.parsed()) 
    {
      Util::Op::OpData opDataTmp(0, &solution_vector, 0, &state_vector, &store_vector, 0);
      result = expression_data.evaluate(comm, circuit_time, circuit_dt, opDataTmp);
    }

    if (Parallel::rank(comm) == 0)
    {
      (*os_) << delim << std::setw(width) << result;
    }
  }

  if (Parallel::rank(comm) == 0)
    (*os_) << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : OutputResults::notify
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void OutputResults::notify(
  const Analysis::StepEvent &   step_event)
{
  switch (step_event.state_) 
  {
    case Analysis::StepEvent::INITIALIZE:
      setup(comm_, outputManager_);
      break;

    case Analysis::StepEvent::STEP_STARTED:
      break;

    case Analysis::StepEvent::STEP_COMPLETED:
      output(comm_, 
          step_event.finalSimTime_, 
          step_event.finalSimDt_, 
          step_event.stepSweepVector_, step_event.count_, 
        *analysisManager_.getDataStore()->currSolutionPtr,
        *analysisManager_.getDataStore()->currStatePtr,
        *analysisManager_.getDataStore()->currStorePtr,
        *analysisManager_.getDataStore()->currLeadCurrentPtr,
        *analysisManager_.getDataStore()->currLeadDeltaVPtr);
      break;

    case Analysis::StepEvent::FINISH:
      steppingComplete();
      break;
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputResults::steppingComplete
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void OutputResults::steppingComplete()
{
  // Deal with the result file:
  if (os_)
  {
    (*os_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
  }

  delete os_;
  os_ = 0;
}

namespace {

//-----------------------------------------------------------------------------
// Class         : ResultsFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jan 29 12:53:02 2015
//-----------------------------------------------------------------------------
///
/// Factory for parsing Results parameters from the netlist and creating Results analysis.
///
class ResultsFactory : public Util::Factory<Analysis::ProcessorBase, OutputResults>
{
public:
  //-----------------------------------------------------------------------------
  // Function      : ResultsFactory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:54:09 2015
  //-----------------------------------------------------------------------------
  ///
  /// Constructs the Results output factory
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
  /// @param topology
  ///
  ResultsFactory(
    Parallel::Machine           comm,
    Analysis::AnalysisManager & analysis_manager,
    OutputMgr &                 output_manager)
    : comm_(comm),
      analysisManager_(analysis_manager),
      outputManager_(output_manager)
  {}

  virtual ~ResultsFactory()
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
  /// Create a new Results analysis and applies the analysis and time integrator option blocks.
  ///
  /// @return new Results analysis object
  ///
  OutputResults *create() const
  {
    OutputResults *results = new OutputResults(comm_, analysisManager_, outputManager_);
    for (std::vector<Util::OptionBlock>::const_iterator it = resultsOptionBlocks_.begin(), 
        end = resultsOptionBlocks_.end(); it != end; ++it)
    {
      results->addResultParams( analysisManager_.getExpressionGroup(), *it);
    }

    return results;
  }

public:
  Parallel::Machine                     comm_;
  Analysis::AnalysisManager &           analysisManager_;
  OutputMgr &                           outputManager_;

  std::vector<Util::OptionBlock>        resultsOptionBlocks_;
};

// .RESULTS
struct ResultOptionsReg : public IO::PkgOptionsReg
{
  ResultOptionsReg(ResultsFactory &factory)
    : factory_(factory)
  {}

  bool operator()(const Util::OptionBlock & option_block)
  {
    factory_.resultsOptionBlocks_.push_back(option_block);

    factory_.analysisManager_.addProcessor(&factory_);

    return true;
  }

private:
  ResultsFactory &              factory_;
};

// .STEP
struct StepAnalysisReg : public IO::PkgOptionsReg
{
  StepAnalysisReg(ResultsFactory &factory)
    : factory_(factory)
  {}

  bool operator()(const Util::OptionBlock & option_block)
  {
    factory_.analysisManager_.addProcessor(&factory_);

    return true;
  }

private:
  ResultsFactory &              factory_;
};

//-----------------------------------------------------------------------------
// Function      : extractRESULTData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 08/29/2004
//-----------------------------------------------------------------------------
bool extractRESULTData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("RESULT", 
      Util::OptionBlock::ALLOW_EXPRESSIONS, 
      netlist_filename, 
      parsed_line[0].lineNumber_);

  // print out the parsed line
  if (DEBUG_IO)
  {
    for (int ieric=0;ieric<parsed_line.size();++ieric)
    {
      Xyce::dout() << "parsed_line["<<ieric<<"] = " << parsed_line[ieric].string_ << std::endl;
    }
  }

  int linePosition = 1;   // Start of parameters on .param line.
  Util::Param parameter("", "");

  if (parsed_line.size () <= 1)
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << "Too few fields in the .RESULT line";
    return false;
  }

  parameter.setTag( "EXPRESSION" );
  parameter.setVal( parsed_line[linePosition].string_ );
  option_block.addParam( parameter );

  if (linePosition != ((int)parsed_line.size () - 1))
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << "Too many fields in the .RESULT line";
    return false;
  }

  circuit_block.addOptions(option_block);

  return true;
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : registerOutputResultsFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool registerOutputResultsFactory(
  Analysis::FactoryBlock &      factory_block,
  Parallel::Machine             comm)
{
  ResultsFactory *factory = new ResultsFactory(comm, factory_block.analysisManager_, factory_block.outputManager_);

  addProcessorFactory(factory_block, factory);

  factory_block.optionsManager_.addCommandParser(".RESULT", extractRESULTData);

  factory_block.optionsManager_.addCommandProcessor("RESULT", new ResultOptionsReg(*factory));
  factory_block.optionsManager_.addCommandProcessor("STEP", new StepAnalysisReg(*factory));
  factory_block.optionsManager_.addCommandProcessor("SAMPLING", new StepAnalysisReg(*factory));

  return true;
}

} // namespace IO
} // namespace Xyce

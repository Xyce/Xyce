//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Purpose       : This file contains the functions for a manager class for
//                 the handling analysis objects generated by .fft statements
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_DEV_DeviceMgr.h>
#include <N_ERH_ErrorMgr.h>

#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>
#include <N_IO_FFTMgr.h>
#include <N_IO_FFTAnalysis.h>
#include <N_IO_NetlistImportTool.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_PkgOptionsMgr.h>

#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Interpolators.h>
#include <N_UTL_Math.h>
#include <N_UTL_Op.h>
#include <N_UTL_OpBuilder.h>
#include <N_UTL_OptionBlock.h>

#include <Teuchos_ScalarTraits.hpp>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : FFTMgr::FFTMgr
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
FFTMgr::FFTMgr(const CmdParse &cp)
  : commandLine_(cp),
    fftAnalysisEnabled_(false),
    fft_accurate_(true),
    fftout_(false),
    fft_mode_(0)
{}

//-----------------------------------------------------------------------------
// Function      : FFTMgr::FFTMgr
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
FFTMgr::~FFTMgr()
{
  for (FFTAnalysisVector::iterator it = FFTAnalysisList_.begin(); it != FFTAnalysisList_.end(); ++it)
    delete (*it);
}

//-----------------------------------------------------------------------------
// Function      : FFTMgr::notify
// Purpose       : Reset each FFTAnalysis object at the start of a STEP iteration,
//                 and output the FFT results at the end of each STEP iteration.
//                 Output for the non-step case is currently handled by
//                 outputMacroResults().
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/27/2021
//-----------------------------------------------------------------------------
void FFTMgr::notify( const Analysis::StepEvent & step_event)
{
  switch (step_event.state_)
  {
    case Analysis::StepEvent::INITIALIZE:
      break;

    case Analysis::StepEvent::STEP_STARTED:
      resetFFTAnalyses();
      break;

    case Analysis::StepEvent::STEP_COMPLETED:
      outputResultsToFFTfile(step_event.count_);
      outputVerboseResults(Xyce::lout());
      break;

    case Analysis::StepEvent::FINISH:
      break;
  }
}

//-----------------------------------------------------------------------------
// Function      : FFTMgr::resetFFTAnalyses
// Purpose       : Resets all of the FFTAnalysis objects.  Typically called at
//                 the start of a .STEP loop iteration.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 2/15/2021
//-----------------------------------------------------------------------------
void FFTMgr::resetFFTAnalyses()
{
  for (FFTAnalysisVector::iterator it = FFTAnalysisList_.begin(); it != FFTAnalysisList_.end(); ++it)
    (*it)->reset();
}

//-----------------------------------------------------------------------------
// Function      : FFTMgr::enableFFTAnalysis
// Purpose       : Determines whether the FFT analysis should actually be done,
//                 based on the primary analysis mode
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 2/28/2021
//-----------------------------------------------------------------------------
void FFTMgr::enableFFTAnalysis(const Analysis::Mode analysisMode)
{
  if ( (analysisMode == Xyce::Analysis::ANP_MODE_TRANSIENT) && !FFTAnalysisList_.empty() )
    fftAnalysisEnabled_ = true;
}

//-----------------------------------------------------------------------------
// Function      : FFTMgr::fixupFFTParameters
// Purpose       : This sets parameters in the FFTAnalysis objects, that could
//                 not be determined when those objects were constructed.  It
//                 may also modify the fft_accurate_ setting if it is incompatible
//                 with the .OPTIONS OUTPUT settings.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
void FFTMgr::fixupFFTParameters(
  Parallel::Machine comm,
  const IO::OutputMgr &output_manager,
  const Util::Op::BuilderManager &op_builder_manager,
  double endSimTime,
  TimeIntg::StepErrorControl & sec)
{
  if (fftAnalysisEnabled_)
  {
    // coordinate settings between OutputMgr and FFTMgr
    if ( (fft_accurate_ == 1) && (output_manager.getInitialOutputInterval() > 0.0) )
    {
      fft_accurate_=0;
      Report::UserWarning0() << "FFT_ACCURATE reset to 0, because .OPTIONS OUTPUT INITIAL_INTERVAL used";
    }

    // now fixup the individual FFTAnalysis objects
    for (FFTAnalysisVector::iterator it = FFTAnalysisList_.begin(); it != FFTAnalysisList_.end(); ++it)
      (*it)->fixupFFTParameters(comm, op_builder_manager, endSimTime, sec, fft_accurate_, fftout_, fft_mode_);
  }
}

//-----------------------------------------------------------------------------
// Function      : FFTMgr::fixupFFTParametersForRemeasure
// Purpose       : This sets parameters in the FFTAnalysis objects, that could
//                 not be determined when those objects were constructed.
// Special Notes : fft_accurate_ is set to false in all of the FFTAnalysis objects
//                 for remeasure, because StepErrorControl is not instantiated
//                 during remeasure.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
void FFTMgr::fixupFFTParametersForRemeasure(
  Parallel::Machine comm,
  const Util::Op::BuilderManager &op_builder_manager,
  double endSimTime,
  TimeIntg::StepErrorControl & sec)
{
  if (fftAnalysisEnabled_)
  {
    for (FFTAnalysisVector::iterator it = FFTAnalysisList_.begin(); it != FFTAnalysisList_.end(); ++it)
      (*it)->fixupFFTParameters(comm, op_builder_manager, endSimTime, sec, false, fftout_, fft_mode_);
  }
}

//-----------------------------------------------------------------------------
// Function      : FFTMgr::addFFTAnalysis
// Purpose       : Entry point when .fft lines are parsed in the netlist.  An
//                 individual FFTAnalysis object is constructed for each .FFT
//                 line in the netlist.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
bool FFTMgr::addFFTAnalysis(const Util::OptionBlock & fftBlock)
{
  // based on what's in the option block passed in, we create the needed FFTAnalysis instance
  if (DEBUG_IO)
    Xyce::dout() << "In FFTMgr::addFFTAnalysis" << std::endl
                 << "dot FFT line passed in was: " << std::endl << fftBlock << std::endl;

  FFTAnalysis  * theFFTObject = 0;
  theFFTObject = new FFTAnalysis(fftBlock);

  if (theFFTObject)
  {
    FFTAnalysisList_.push_back( theFFTObject );
    // Used to help register lead current requests with device manager.
    // FFT manager keeps the combined list, based on parsing of
    // dependent solution variable vector generated by each .FFT line.
    getLeadCurrentDevices(theFFTObject->getDepSolVarIterVector(), devicesNeedingLeadCurrents_);
  }
  else
  {
    delete theFFTObject;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : FFTMgr::updateFFTData
// Purpose       : Called during the simulation to update the vectors (in each
//                 FFTAnalysis object) with the values of the time-points and
//                 their respective output variable.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
void FFTMgr::updateFFTData(Parallel::Machine comm, double circuitTime, const Linear::Vector *solnVec,
  const Linear::Vector *stateVec, const Linear::Vector * storeVec,
  const Linear::Vector *lead_current_vector, const Linear::Vector *junction_voltage_vector,
  const Linear::Vector *lead_current_dqdt_vector)
{
  if (fftAnalysisEnabled_)
  {
    // loop over FFTAnalysis objects and get them to update themselves.
    for (FFTAnalysisVector::iterator it = FFTAnalysisList_.begin(); it != FFTAnalysisList_.end(); ++it)
    {
      (*it)->updateFFTData(comm, circuitTime, solnVec, stateVec, storeVec,
                        lead_current_vector, junction_voltage_vector, lead_current_dqdt_vector);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : FFTMgr::outputResultsToFFTFile
// Purpose       : Output all of the FFT results at end of simulation for the
//                 .STEP and -remeasure cases
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 2/15/2021
//-----------------------------------------------------------------------------
void FFTMgr::outputResultsToFFTfile(int stepNumber)
{
  if (isFFTActive())
  {
    std::string filename = IO::makeOutputFileNameWithStepNum(commandLine_, ".fft", stepNumber);
    std::ofstream outputFileStream;
    outputFileStream.open( filename.c_str() );
    outputResults(outputFileStream);
    outputFileStream.close();
  }
}

//-----------------------------------------------------------------------------
// Function      : FFTMgr::outputResults
// Purpose       : Output all of the FFT results at end of simulation
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
void FFTMgr::outputResults(std::ostream& outputStream)
{
  for (FFTAnalysisVector::iterator it = FFTAnalysisList_.begin(); it != FFTAnalysisList_.end(); ++it)
  {
    (*it)->outputResults(outputStream);
  }
}

//-----------------------------------------------------------------------------
// Function      : FFTMgr::outputVerboseResults
// Purpose       : Output metrics values and sorted harmonic list to stdout,
//                 at end of simulation, for .OPTIONS FFT FFTOUT=1
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 8/29/2021
//-----------------------------------------------------------------------------
void FFTMgr::outputVerboseResults(std::ostream& outputStream)
{
  if (isFFTActive() && fftout_)
  {
    outputStream << std::endl
                 << " ***** FFT Analyses ***** " << std::endl
                 << std::endl;

    for (FFTAnalysisVector::iterator it = FFTAnalysisList_.begin(); it != FFTAnalysisList_.end(); ++it)
    {
      (*it)->outputVerboseResults(outputStream);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : FFTMgr::registerFFTOptions
// Purpose       : registers set of variables to set for .OPTIONS FFT
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
bool FFTMgr::registerFFTOptions(const Util::OptionBlock &option_block)
{
  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; )
  {
    if ((*it).tag() == "FFT_ACCURATE")
    {
      int fft_accurate_entered = (*it).getImmutableValue<int>();
      // need to point at next parameter
      ++it;
      if ((fft_accurate_entered < 0) || (fft_accurate_entered > 1))
      {
        fft_accurate_ = true;
	Report::UserWarning0() << "FFT_ACCURATE values of 0 or 1 are supported.  Setting to default of 1";
      }
      else
      {
        fft_accurate_ = fft_accurate_entered;
      }
    }
    else if ((*it).tag() == "FFTOUT")
    {
      fftout_ = (*it).getImmutableValue<int>();
      ++it;
    }
    else if ((*it).tag() == "FFT_MODE")
    {
      int fft_mode_entered = (*it).getImmutableValue<int>();
      // need to point at next parameter
      ++it;
      if ((fft_mode_entered < 0) || (fft_mode_entered > 1))
      {
        fft_mode_ = 0;
	Report::UserWarning0() << "FFT_MODE values of 0 or 1 are supported.  Setting to default of 0";
      }
      else
      {
        fft_mode_ = fft_mode_entered;
      }
    }
    else
    {
      // silently ignore?
      ++it;
    }
  }

  return true;
}

namespace {

//-----------------------------------------------------------------------------
// Class         : FFTOptionsReg
// Purpose       : functor for registering FFT options
// Special Notes : Used by package manager addOptionsProcessor method
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
struct FFTOptionsReg : public PkgOptionsReg
{
  FFTOptionsReg( FFTMgr &mgr )
    : fft_manager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  {
    return fft_manager_.addFFTAnalysis( options );
  }

  FFTMgr &fft_manager_;
};

//-----------------------------------------------------------------------------
// Function      : extractFFTData
// Purpose       : Convert a .fft line to an options block
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 1/4/2020
//-----------------------------------------------------------------------------
bool extractFFTData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  // This is the set of allowed qualifiers on the .FFT lines
  std::set<std::string> keywords;
  keywords.insert( std::string("START") );
  keywords.insert( std::string("FROM") );
  keywords.insert( std::string("STOP") );
  keywords.insert( std::string("TO") );
  keywords.insert( std::string("NP") );
  keywords.insert( std::string("FORMAT") );
  keywords.insert( std::string("WINDOW") );
  keywords.insert( std::string("ALFA") );
  keywords.insert( std::string("FREQ") );
  keywords.insert( std::string("FMIN") );
  keywords.insert( std::string("FMAX") );

  Util::OptionBlock option_block("DOT_FFT_LINE", Util::OptionBlock::ALLOW_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

  const int numFields = parsed_line.size();

  // used to ensure that the output variable is not preceded by any of the allowed qualifiers
  bool outputVarFound=false;

  // Create an option block to temporarily store the default options.
  Util::OptionBlock defaultOptions;

  // Get the default options from metadata.
  //addDefaultOptionsParameters(options_manager, defaultOptions, "FFT" );

  Util::Param parameter;
  ExtendedString nextWord("");

  if(numFields < 2)
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << "Error: .FFT line needs at least 2 arguments '.FFT ov1 '";
    return false;
  }

  int position = 1;
  int endPosition = numFields;
  while (position < endPosition)
  {

    nextWord = parsed_line[position].string_;
    nextWord.toUpper();
    // nextWord[0] is used for I because branch currents for transistors can have two
    // characters.  An example is IS for the M-Device.
    if( (nextWord[0] == 'I') || (nextWord == "V") || (nextWord == "P") || (nextWord == "W") || (nextWord == "N") )
    {
      // need to do a bit of look ahead here to see if this is a V(a) or V(a,b)
      int numNodes = 1;
      if( ((position+3) < endPosition) && (parsed_line[position+3].string_ == ",") )
      {
        numNodes = 2;
      }
      parameter.set(nextWord, numNodes);
      option_block.addParam( parameter );

      if( (position+2) < endPosition )
      {
        position+=2;
        nextWord = parsed_line[position].string_;
        nextWord.toUpper();
        parameter.set( nextWord, 0.0 );
        option_block.addParam( parameter );

        if( ((position+2) < endPosition) && (parsed_line[position+1].string_ == ",") )
        {
          // may have a voltage difference request
          position+=2;
          nextWord = parsed_line[position].string_;
          nextWord.toUpper();
          parameter.set( nextWord, 0.0 );
          option_block.addParam( parameter );
        }
      }
      else
      {
         Report::UserError0().at(netlist_filename, parsed_line[position].lineNumber_) << "Could not parse .FFT variable";
         return false;
      }
      outputVarFound=true;
      position+=2;
    }
    else if ( (nextWord[0]=='{') && (nextWord[nextWord.size()-1]=='}') )
    {
      parameter.set(nextWord, std::string(nextWord));
      option_block.addParam(parameter);
      outputVarFound=true;
      position++;
    }
    else if (!outputVarFound)
    {
      Report::UserError().at(netlist_filename, parsed_line[position].lineNumber_) << "Output variable must be first field on .FFT line";
      return false;
    }
    else if ( (keywords.find(nextWord) != keywords.end()) && (position+2 <= endPosition) &&
              (parsed_line[position+1].string_[0] == '=') )
    {
      parameter.setTag(nextWord);
      parameter.setVal(parsed_line[position+2].string_ );
      option_block.addParam(parameter);
      position +=3;
    }
    else
    {
      Report::UserError().at(netlist_filename, parsed_line[position].lineNumber_) << "Could not parse .FFT variable or keyword " << nextWord;
      return false;
    }
  }

  circuit_block.addOptions(option_block);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : populateMetadata
// Purpose       :
// Special Notes : option blocks with the name FFT will come from .OPTIONS FFT
//                 lines.  option blocks from .FFT lines will have the name
//                 DOT_FFT_LINE.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 01/4/2021
//-----------------------------------------------------------------------------
void populateMetadata(IO::PkgOptionsMgr &   options_manager)
{
  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("FFT");

    parameters.insert(Util::ParamMap::value_type("FFT_ACCURATE", Util::Param("FFT_ACCURATE", 1)));
    parameters.insert(Util::ParamMap::value_type("FFTOUT", Util::Param("FFTOUT", 1)));
    parameters.insert(Util::ParamMap::value_type("FFT_MODE", Util::Param("FFT_MODE", 1)));
  }
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
bool registerPkgOptionsMgr(FFTMgr &fft_manager, PkgOptionsMgr &options_manager)
{
  populateMetadata(options_manager);

  options_manager.addCommandParser(".FFT", extractFFTData);

  // This handles .FFT lines
  options_manager.addCommandProcessor("DOT_FFT_LINE", new FFTOptionsReg(fft_manager));

  // This handles .OPTIONS FFT lines
  options_manager.addOptionsProcessor("FFT",
      IO::createRegistrationOptions(fft_manager, &FFTMgr::registerFFTOptions));

  return true;
}

} // namespace IO
} // namespace Xyce

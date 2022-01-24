//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Purpose        : Used to send input options to registered pkg's
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL
//
// Creation Date  : 1/28/2003
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ERH_Message.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Param.h>

#include <N_LAS_QueryUtil.h>            // Hacked in while fixing the default option/metadata data stuff
#include <N_TIA_TIAParams.h>            // Hacked in while fixing the default option/metadata data stuff

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : PkgOptionsMgr::PkgOptionsMgr
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 10/21/2008
//-----------------------------------------------------------------------------
PkgOptionsMgr::PkgOptionsMgr()
{
  addCommandParser(".OPTION", extractOptionsData);
  addCommandParser(".OPTIONS", extractOptionsData);
  addCommandParser(".DATA", extractDotDataStatement);

  Linear::QueryUtil::populateMetadata(*this);
  TimeIntg::TIAParams::populateMetadata(*this);
}

//-----------------------------------------------------------------------------
// Function      : PkgOptionsMgr::~PkgOptionsMgr
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL
// Creation Date : 1/28/2003
//-----------------------------------------------------------------------------
PkgOptionsMgr::~PkgOptionsMgr()
{
  for (ProcessorMap::iterator it = processorMap_.begin(), end = processorMap_.end(); it != end; ++it)
    delete (*it).second;
}

//-----------------------------------------------------------------------------
// Function      : PkgOptionsMgr::addOptionsProcessor
// Purpose       : Register an object for recieving options
//
// Special Notes : The design of this function and the other 'submit' function
//                 is intended to be robust enough that it doesn't matter
//                 which call is done first.  If the options are submitted
//                 first, then the addOptionsProcessor function will send
//                 them out to the various registration functors.  If the
//                 registration functors are submitted first, then when the
//                 options are submitted, they'll be sent out.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL
// Creation Date : 1/28/2003
//-----------------------------------------------------------------------------
bool
PkgOptionsMgr::addOptionsProcessor(
  const std::string &   option_block_name,
  PkgOptionsReg *       processor)
{
  processorMap_.insert(ProcessorMap::value_type(option_block_name, processor));

  // Call the processor functions if they are interested in this circuit
  std::pair<OptionsMap::iterator, OptionsMap::iterator> range = optionsMap_.equal_range(option_block_name);
  for (OptionsMap::iterator it = range.first; it != range.second; ++it)
    (*processor)((*it).second);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PkgOptionsMgr::submitOptions
// Purpose       : Register an option block to be given to pkg's
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL
// Creation Date : 1/28/2003
//-----------------------------------------------------------------------------
bool PkgOptionsMgr::submitOptions(
  const Util::OptionBlock &     option_block)
{
  const std::string &option_block_name = option_block.getName();

  Util::OptionBlock dispatch_option_block;
  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; ++it)
  {
    if ((*it).tag() == "SINGLETON" && (*it).getImmutableValue<int>()
        && optionsMap_.count(option_block_name))
      Report::UserError().at(option_block.getNetlistLocation()) << "Only 1 option line of type " << option_block_name << " allowed";
    else
      dispatch_option_block.addParam(*it);
  }

  optionsMap_.insert(OptionsMap::value_type(option_block_name, dispatch_option_block));

  std::pair<ProcessorMap::iterator, ProcessorMap::iterator> range = processorMap_.equal_range(option_block_name);
  for (ProcessorMap::iterator it = range.first; it != range.second; ++it)
    (*(*it).second)(dispatch_option_block);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : extractOptionsData
// Purpose       : Extract the parameters from a netlist .OPTIONS line held in
//                 parsed_line.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/05/2001
//-----------------------------------------------------------------------------
bool
extractOptionsData(
  PkgOptionsMgr &       options_manager,
  CircuitBlock &        circuit_block,
  const std::string &   netlist_filename,
  const TokenVector &   parsed_line)
{
  const int numFields = parsed_line.size();
  // catch case where the input line is .OPTIONS, without any further info
  if (numFields == 1)
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
          << ".OPTIONS line missing package name";
    return false;
  }

  // The type of options is given by the second field on the .options line
  // UNLESS they are SPICE style options. Check the second field and set the
  // name attribute appropriately.
  ExtendedString optionName ( parsed_line[1].string_ );
  optionName.toUpper();

  int parameterStartPos = 2;

  if (!options_manager.findOptionsMetadata(optionName))
  {
    // Check to see if this is a convertable SPICE style option line or if an unrecognized
    // package name was given. Do this by checking for an "=" in either the
    // 3rd or 4th position on the line.
    if (numFields > 2 && parsed_line[2].string_ == "=")
    {
      // look for TEMP and/or TNOM  If they are the only tags present on the
      // rest of the line, they try setting the package name to DEVICE
      //
      // assume we only found TEMP and TNOM for now.
      bool foundOnlyTempAndTnom=true;

      // format of parsed line at this point is <tag> = <value>
      // with the first tag starting at position 1 (hence why
      // parsed_line[2].string_ == "="  was true to get into this if block.)
      // Check every 3rd item starting from i=1 to only check the tags.
      for (int i = 1; i < numFields; i = i + 3)
      {
        ExtendedString tag ( parsed_line[i].string_ );
        tag.toUpper();
        if ((tag != "TNOM") && (tag != "TEMP"))
        {
          foundOnlyTempAndTnom=false;
        }
      }

      // Now do a similar procedure for SCALE.
      // If SCALE is here by itself, try setting the optionName to PARSER
      bool foundOnlyScale=true;
      for (int i = 1; i < numFields; i = i + 3)
      {
        ExtendedString tag ( parsed_line[i].string_ );
        tag.toUpper();
        if ((tag != "SCALE"))
        {
          foundOnlyScale=false;
        }
      }

      if (foundOnlyTempAndTnom)
      {
        optionName = "DEVICE";
        parameterStartPos = 1;
      }
      else if (foundOnlyScale)
      {
        optionName = "PARSER";
        parameterStartPos = 1;
      }
      else
      {
        optionName = "PSPICE";  // labeled as "PSPICE" for historical reasons, but this covers many varieties of SPICE.
        parameterStartPos = 1;

        Report::UserWarning0().at(netlist_filename, parsed_line[0].lineNumber_)
          << "Option name missing, or SPICE-style option not supported in Xyce";

        return false;
      }
    }
    else
    {
      Report::UserWarning0().at(netlist_filename, parsed_line[0].lineNumber_)
        << "Unrecognized .OPTIONS package or parameter " << optionName << " will be ignored";
      return true;
    }
  }

  Util::OptionBlock option_block(
      optionName,
    (  equal_nocase(optionName, "SENSITIVITY") || 
       equal_nocase(optionName, "SAMPLES") ||
       equal_nocase(optionName, "EMBEDDEDSAMPLES") ||
       equal_nocase(optionName, "PCES" )
       ) ? Util::OptionBlock::ALLOW_EXPRESSIONS : Util::OptionBlock::NO_EXPRESSIONS,
      netlist_filename, parsed_line[0].lineNumber_);

  // Create an option block to temporarily store the default options.
  Util::OptionBlock defaultOptions;

  // Get the default options from metadata.
  addDefaultOptionsParameters(options_manager, defaultOptions, optionName);

  // generate a warning message for lines like .OPTIONS LINSOL TYPE
  if ( parameterStartPos > (numFields - 2 ) )
  {
    Report::UserWarning0().at(netlist_filename, parsed_line[0].lineNumber_)
        << ".OPTIONS line is missing one or more required parameters, and was ignored.";
  }
  // Extract the parameters from parsed_line.
  Util::ParamList inputParameters;
  int intervalParameterStart = -1;
  for (int i = parameterStartPos, parameterEndPos = numFields - 1; i <= parameterEndPos - 1; )
  {
    // Check for equal sign.  
    if ( parsed_line[i + 1].string_ != "=" )
    {
      if ( optionName != "OUTPUT" && optionName != "RESTART" )
      {
        // this error occurs if <name>=<val> syntax is mis-formatted.  Examples include:
        //   .OPTIONS LINSOL=KSPARSE
        //   .OPTIONS LINSOL =KSPARSE
        //   .OPTIONS LINSOL type KSPARE
        Report::UserError0().at(netlist_filename, parsed_line[i].lineNumber_) 
          << "Misformatted or incorrect <name>=<val> syntax in .OPTIONS " << optionName;
        return false;
      }
      else
      {
        // Stop after the tagged parameters have been extracted
        // from a .OPTIONS RESTART or .OPTIONS OUTPUT line, they
        // will be handled later.
        intervalParameterStart = i;
        break;
      }
    }

    // Extract parameter name and value from parsedLine and add to
    // parameter list. Check to see if the parameter is "VECTOR"
    // valued and treat accordingly.
    const std::string &parameter_name = parsed_line[i].string_;
    Util::Param *parameterPtr = Util::findParameter(defaultOptions.begin(), defaultOptions.end(), parameter_name);
    if (parameterPtr == NULL)
    {
      Report::UserWarning0().at(netlist_filename, parsed_line[0].lineNumber_)
        << "No options parameter " << parameter_name << " found, parameter will be ignored.";
      i += 3;
    }
    else if (parameterPtr->stringValue() != "VECTOR")
    {
      // error out if the next push_back will cause a segfault.  The .OPTIONS line is too
      // short, and is missing a parameter
      if ( (i+2) > (numFields-1) )
      {
        Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
          << ".OPTIONS line is missing one or more required parameters";
        return false;
      }
      inputParameters.push_back(Util::Param(parsed_line[i].string_, parsed_line[i + 2].string_));
      i += 3;

      if (i < parameterEndPos - 1 && parsed_line[i].string_ == ",")
      {
        Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
          << "Options parameter " << parameter_name << " is not a vector, but has comma in value";
        return false;
      }
    }
    else
    {
      // We have a vector valued parameter.
      // Name the jth component of the parameter of the vector by appending
      // "j" to the parameter name.
      std::ostringstream paramName;
      std::string paramBaseName = ExtendedString(parameter_name).toUpper();
      int j = 1;

      // used to help stop reading off the end of parsed_line
      int testSize= parsed_line.size() - 1;

      paramName << paramBaseName << j;
      i += 2;
      if ( i > testSize)
      {
        Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
          << ".OPTIONS line is missing one or more elements in a vector parameter";
        return false;
      }
      else if (parsed_line[i].string_ == ",")
      {
        // catch the invalid case of option=,value,value
        Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
            << ".OPTIONS line is mis-formatted. Vector parameter begins with a comma";
        return false;
      }

      option_block.addParam(Util::Param(paramName.str(), parsed_line[i].string_));

      // The nominal format is option=value,value,value   
      // However, this loop should error out on invalid lines that have:
      // errors like:
      //
      //        option=value,value,value,
      //        option=value,value,,value
      //
      // or basically any case where the "parameter value" is a comma,
      // which indicates that something went wrong during parsing.
      while ((i < testSize) && (parsed_line[i + 1].string_ == ",") )
      {
        paramName.str("");
        ++j;
        paramName << paramBaseName << j;
        i += 2;
        if ( (i <= testSize) && (parsed_line[i].string_ != ",") )
        {         
          option_block.addParam(Util::Param(paramName.str(), parsed_line[i].string_));
        }
        else
        {
          Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
            << ".OPTIONS line is mis-formatted. Error parsing vector parameter";
          return false;
        }
      }

      ++i;
      
      // Also need to guard against these cases:
      //
      //     option1=value,value value
      //     option1=value,value value option2=value 
      if ( (i==testSize) || (( i+1 <= testSize) && (parsed_line[i+1].string_ != "=")) )
      {
        Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
          << ".OPTIONS line is mis-formatted. Possibly a missing comma in vector parameter";
        return false;
      }
    }
  }

  // For each input parameter, check that it is in the default
  // set and if so, set its value in "parameters" to the input
  // value, otherwise flag it as an unknown parameter.
  for (Util::ParamList::const_iterator it = inputParameters.begin(), end = inputParameters.end(); it != end; ++it)
  {
    Util::Param *parameterPtr = Util::findParameter(defaultOptions.begin(), defaultOptions.end(), (*it).tag());
    if ( parameterPtr != NULL )
    {
      parameterPtr->setVal(*it);
      option_block.addParam(*parameterPtr);
    }
    else
    {
      Report::UserWarning0().at(netlist_filename, parsed_line[0].lineNumber_)
        << "No options parameter " << (*it).tag() << " found, parameter will be ignored.";
    }
  }

  // If this is a ".OPTIONS OUTPUT" line or ".OPTIONS RESTART" line
  // then get the time and interval pairs.
  if (intervalParameterStart != -1 && (optionName == "OUTPUT" || optionName == "RESTART"))
  {
    for (int i = intervalParameterStart, parameterEndPos = numFields - 1;  i < parameterEndPos; )
    {
      option_block.addParam(Util::Param("TIME", parsed_line[i].string_));
      ++i;

      option_block.addParam(Util::Param("INTERVAL", parsed_line[i].string_));
      ++i;
    }
  }

  circuit_block.addOptions(option_block);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : extractDotDataStatement
// Purpose       : read in a table of data. 
// Special Notes : this is used by multiple packages, so putting here
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/5/18
//-----------------------------------------------------------------------------
bool extractDotDataStatement(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  // length of the original .DATA line
  const int numFields = parsed_line.size();

  // .DATA line must have a NAME field and at least one param and value
  if (numFields < 4)
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
        << ".DATA line lacks enough fields";
    return false;
  }

  // start of parameters (skip over the ".DATA")
  int linePosition = 1;

  // get the data name
  std::string dataName = parsed_line[linePosition++].string_;
  Util::toUpper(dataName);

  Util::OptionBlock option_block("DATA", Util::OptionBlock::ALLOW_EXPRESSIONS, netlist_filename, parsed_line[linePosition].lineNumber_);

  option_block.addParam(Util::Param("NAME",dataName));

  // get the parameter list
  while (linePosition < numFields && !(Util::isValue(parsed_line[linePosition].string_)))
  {
    std::string paramName = parsed_line[linePosition].string_;
    Util::toUpper(paramName);
    option_block.addParam(Util::Param("PARAM", paramName));
    linePosition++;
  }

  while (linePosition < numFields && (Util::isValue(parsed_line[linePosition].string_)))
  {
    option_block.addParam(Util::Param("VAL", parsed_line[linePosition].string_));
    linePosition++;
  }

  circuit_block.addOptions(option_block);

  return true;
}

ParamWarning::ParamWarning(const Util::OptionBlock &option_block, const Util::Param &param)
  : Report::UserWarning()
{
  at(option_block.getNetlistLocation());

  os() << param.tag() << ": ";
}

ParamError::ParamError(const Util::OptionBlock &option_block, const Util::Param &param)
  : Report::UserError()
{
  at(option_block.getNetlistLocation());

  os() << param.tag() << ": ";
}

} // namespace IO
} // namespace Xyce

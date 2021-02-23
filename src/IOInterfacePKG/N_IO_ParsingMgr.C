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
// Purpose        : Parsing Manager
//
// Special Notes  :
//
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 04/09/19
//
//
//
//
//-------------------------------------------------------------------------

#include <string>

#include <Xyce_config.h>

#include <N_IO_CmdParse.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_ParsingMgr.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_HspiceBools.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : ParsingMgr::ParsingMgr
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 04/09/19
//-----------------------------------------------------------------------------
ParsingMgr::ParsingMgr(
  const CmdParse &              command_line)
  : hspiceExtFlag_(command_line.argExists("-hspice-ext")),
    useHspiceUnits_(false),
    useHspiceMath_(false),
    enableRandomExpression_(true),
    modelBinningFlag_(false),
    lengthScale_(1.0)
{
  if (hspiceExtFlag_)
  {
    std::string hspiceExtArgStr = command_line.getArgumentValue("-hspice-ext");
    Util::toLower(hspiceExtArgStr);

    // hspiceExtArgStr may contain a comma-separated list of options,
    // so we need to pull it apart.
    std::vector<std::string> hspiceExtVec;
    int tokStart = 0;
    for (int i=0; i< hspiceExtArgStr.length(); ++i)
    {
      if ( hspiceExtArgStr[i] == ',' )
      {
        hspiceExtVec.push_back(hspiceExtArgStr.substr(tokStart,i-tokStart));
        tokStart = i+1;
      }
      else if ( i == hspiceExtArgStr.length()-1 )
      {
        hspiceExtVec.push_back(hspiceExtArgStr.substr(tokStart,hspiceExtArgStr.length()-tokStart));
      }
    }

    // now parse the strings, and set the appropriate flags to true.
    // Error out, if any invalid or blank entries are found.
    for (std::vector<std::string>::const_iterator it = hspiceExtVec.begin(); it != hspiceExtVec.end(); ++it)
    {
      if (*it == "all")
      {
        useHspiceUnits_ = true;
        useHspiceMath_ = true;
        enableRandomExpression_ = false;
      }
      else if (*it == "units")
        useHspiceUnits_ = true;
      else if (*it == "math")
        useHspiceMath_ = true;
      else if (*it == "random")
        enableRandomExpression_ = false;
      else
        Report::UserFatal0() << "Invalid value " << *it << " for -hspice-ext command line option";
    }
  }

  // These variables are used, in lieu of passing const references to the
  // ParsingMgr into various Util functions like isValue() and Value()
  // and ExpressionInternals::tokenize_.  If enableRandomExpression_ is
  // set to false, via -hspice-ext random, then AGAUSS() and GAUSS() will
  // just return the mean rather than a random number.
  Xyce::Util::useHspiceUnits = useHspiceUnits_;
  Xyce::Util::useHspiceMath = useHspiceMath_;
  Xyce::Util::enableRandomExpression = enableRandomExpression_;
}

//-----------------------------------------------------------------------------
// Function      : ParsingMgr::~ParsingMgr
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL
// Creation Date : 04/09/19
//-----------------------------------------------------------------------------
ParsingMgr::~ParsingMgr()
{}

//-----------------------------------------------------------------------------
// Function      : ParsingMgr::setParserOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool ParsingMgr::setParserOptions(const Util::OptionBlock & OB)
{
  for (Util::ParamList::const_iterator it = OB.begin(), end = OB.end(); it != end; ++it)
  {
    const std::string &tag = (*it).uTag();

    if (tag == "MODEL_BINNING")
    {
      modelBinningFlag_ = static_cast<bool>((*it).getImmutableValue<int>());
    }
    else if (tag == "SCALE")
    {
      lengthScale_ = ((*it).getImmutableValue<double>());
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Namespace     : Unnamed
// Purpose       : file-local scoped methods and data
// Special Notes :
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
namespace {

//-----------------------------------------------------------------------------
// Function      : populateMetadata
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Baur
// Creation Date : May 12, 2015
//-----------------------------------------------------------------------------
void populateMetadata(
  IO::PkgOptionsMgr &   options_manager)
{
   Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("PARSER");
   parameters.insert(Util::ParamMap::value_type("MODEL_BINNING", Util::Param("MODEL_BINNING", 0)));
   parameters.insert(Util::ParamMap::value_type("SCALE", Util::Param("SCALE", 1.0)));
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : registerPkgOptionsMgr
// Purpose       : Add an options parser for .OPTIONS PARSER
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 4/9/2019
//-----------------------------------------------------------------------------
///
/// Register various OutputManager functions with the options manager.
///
/// @param output_manager Reference to an OutputMgr object
/// @param options_manager Reference to a PkgOptionsMgr object
///
/// Registers options block processors for PARSING options block
/// types, not all of them associated with output.
///
/// @note Many of the functions registered here as command processors (or
///  rather, options block processors) have "command parsers" registered
///  in completely different files.  Command parsers take input netlist
///  lines and create options blocks, command processors take options blocks
///  and do something with them.  The organization of the code is somewhat
///  haphazard, and for some reason (probably an abandoned refactor)
///  the command processors and command parsers are not consistently in the
///  same files.  For example, the "registerSens" function is registered here
///  to process "SENS" options blocks, but the associated command parser
///  to create "SENS" options blocks from .SENS lines is actually in the
///  nonlinear solver package's N_NLS_Sensitivity.C file.
///  It does make this code design rather hard to follow.
///  This may be a good reason to go through a reorganization of
///  code at a later date.
bool registerPkgOptionsMgr(ParsingMgr & parsing_manager, PkgOptionsMgr &options_manager)
{
  populateMetadata(options_manager);

  options_manager.addOptionsProcessor("PARSER",
      IO::createRegistrationOptions(parsing_manager, &ParsingMgr::setParserOptions));

  return true;
}

} // namespace IO
} // namespace Xyce

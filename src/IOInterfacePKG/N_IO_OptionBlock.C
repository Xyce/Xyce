//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose        : Define the N_IO_OptionBlock class an instantiation of
//                  which is associated with a netlist .PARAM or .OPTIONS
//                  line.
//
// Special Notes  :
//
// Creator        : Lon Waters, SNL
//
// Creation Date  : 09/10/2001
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iterator>
#include <iostream>
#include <algorithm>
#include <set>

#include <sstream>

// ----------   Xyce Includes   ----------
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_UTL_Expression.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_MachDepParams.h>

namespace Xyce {
namespace IO {


// //-----------------------------------------------------------------------------
// // Function      : OptionBlock::printDiagnostic
// // Purpose       : Output the details of a device block to standard out.
// // Special Notes :
// // Scope         : public
// // Creator       : Lon Waters, SNL
// // Creation Date : 09/10/2001
// //-----------------------------------------------------------------------------
// void OptionBlock::printDiagnostic() const
// {
//   Xyce::dout() << std::endl;
//   Xyce::dout() << "Option Information" << std::endl;
//   Xyce::dout() << "------------------" << std::endl;

//   if (getName() != "")
//   {
//     Xyce::dout() << "input line:" << std::endl;
//     unsigned int size = parsedLine_.size();
//     for (unsigned int i = 0; i < size; ++i)
//     {
//       Xyce::dout() << "  " << parsedLine_[i].string_;
//     }
//     Xyce::dout() << std::endl;
//     Xyce::dout() << "  name: " << getName() << std::endl;
//   }

//   Xyce::dout() << "  parameters: " << std::endl;
//   Util::ParamList::const_iterator paramIter, paramEnd;
//   paramEnd = optionData_.end();
//   paramIter = optionData_.begin();
//   for ( ; paramIter != paramEnd ; ++paramIter)
//   {
//     Xyce::dout() << "  " << paramIter->tag() << "  ";
//     Xyce::dout() << paramIter->stringValue() << std::endl;
//   }
//   Xyce::dout() << std::endl << std::endl;
// }

//-----------------------------------------------------------------------------
// Function      : extractData
// Purpose       : Determine option type and extract the data appropriately.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 05/22/2002
// Completely revamped by Dave Baur, 5/12/2015
//-----------------------------------------------------------------------------
///
/// Attempt to process a netlist line that may or may not be an options line
///
/// @param options_manager Reference to options manager
/// @param circuit_block reference to current circuit block being processed
/// @param[in] netlist_filename File containing line being processed
/// @param[in] parsed_line vector of strings ("tokens") representing the fields of the input line being processed.
///
/// We look at the first word of the input line, and try to find a registered
/// command parser associated with it.  If we find such a parser, run it
/// over the line.
///
/// The commmand parser will create an options_block and add it to the
/// circuit_block for later processing by a registered "command processor"
/// @author Dave Baur
/// @date 5/12/2015
bool extractData(
  PkgOptionsMgr &       options_manager,
  CircuitBlock &        circuit_block,
  const std::string &   netlist_filename,
  const TokenVector &   parsed_line)
{
  ParseFunction f = options_manager.getCommandParser(parsed_line[0].string_);
  if (f)
    return (*f)(options_manager, circuit_block, netlist_filename, parsed_line);

  return false;
}

//-----------------------------------------------------------------------------
// Function      : addDefaultOptionsParameters
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void addDefaultOptionsParameters(
  const PkgOptionsMgr & options_manager,
  Util::OptionBlock &   option_block,
  const std::string &   default_option_name)
{
  const Util::OptionsMetadataMap::mapped_type *option_map = options_manager.findOptionsMetadata(default_option_name);
  if (!option_map)
  {
    throw std::runtime_error("Missing option");
  }

  for (Util::ParamMap::const_iterator it = option_map->begin(), end = option_map->end(); it != end; ++it)
  {
    option_block.addParam((*it).second);
  }
}

} // namespace IO
} // namespace Xyce

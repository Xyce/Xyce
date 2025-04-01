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
//
// Purpose        : Declare the N_IO_OptionBlock class an instantiation of
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

#ifndef N_IO_OPTIONBLOCK_H
#define N_IO_OPTIONBLOCK_H

#include <string>
#include <vector>
#include <iosfwd>

#include <N_IO_fwd.h>

#include <N_ERH_Message.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Param.h>

namespace Xyce {
namespace IO {

void addDefaultOptionsParameters(const PkgOptionsMgr &options_manager, Util::OptionBlock &option_block, const std::string &default_option_name);

// Determine the line type and extract the data appropriately.
bool extractData(PkgOptionsMgr &options_manager, CircuitBlock &circuit_block, const std::string &netlist_filename, const TokenVector &parsed_line);

} // namespace IO
} // namespace Xyce

#endif

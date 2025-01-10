//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : Linear Solver Factory
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 08/12/03
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_DistToolFactory_h
#define Xyce_N_IO_DistToolFactory_h

// ---------- Standard Includes ----------

#include <string>

#include <N_IO_fwd.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_DistributionTool.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_PDS_Comm.h>

#include <N_IO_DistributionTool.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : DistToolFactory
// Purpose       :
// Special Notes :
// Creator       : Heidi Thornquist, SNL
// Creation Date : 08/12/03
//-----------------------------------------------------------------------------
struct DistToolFactory
{
  // Creates a new IO distribution tool
  static DistributionTool * create(
    Parallel::Communicator *                 pdsCommPtr,    
    Util::OptionBlock &                      distOptions,
    CircuitBlock &                           circuitBlock,
    std::map<std::string,FileSSFPair>      & ssfMap,
    std::map<std::string, IncludeFileInfo> & iflMap,
    const std::vector< std::pair< std::string, std::string> > & externalNetlistParams,
    const ParsingMgr                       & parsing_manager);
};

} // namespace IO 
} // namespace Xyce

#endif // Xyce_N_IO_DistToolFactory_h


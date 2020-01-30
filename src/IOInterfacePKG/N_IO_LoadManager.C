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
// Purpose        : Output Manager
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_LoadManager.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Param.h>
#include <N_ERH_Message.h>

namespace Xyce {
namespace IO {

namespace {

//-----------------------------------------------------------------------------
// Class         : LoadOptionsRegx
// Purpose       : functor for registering Load options
// Special Notes : Used by package manager addOptionsProcessor method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct LoadOptionsReg : public PkgOptionsReg
{
  LoadOptionsReg()
  {}

  bool operator()( const Util::OptionBlock & options )
  {
    Report::UserWarning0() << ".LOAD not supported yet.  Use .INCLUDE instead";
    return true;
  }
};

//-----------------------------------------------------------------------------
// Function      : extractLoadData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/11/2007
//-----------------------------------------------------------------------------
bool
extractLoadData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("LOAD", Util::OptionBlock::NO_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

  if (DEBUG_IO)
  {
    for (int ieric=0;ieric<parsed_line.size();++ieric)
    {
      Xyce::dout() << "parsed_line["<<ieric<<"] = " << parsed_line[ieric].string_ << std::endl;
    }
  }

  circuit_block.addOptions(option_block);

  return true;
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool registerPkgOptionsMgr(LoadManager &load, PkgOptionsMgr &options_manager)
{
  options_manager.addCommandParser(".LOAD", extractLoadData);

  options_manager.addCommandProcessor("LOAD", new LoadOptionsReg());

  return true;
}

} // namespace IO
} // namespace Xyce

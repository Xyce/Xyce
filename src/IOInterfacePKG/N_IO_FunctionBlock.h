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
// Purpose        : Declare an N_IO_FunctionBlock instantiations of which are
//                  associated with netlist .FUNC lines.
//
// Special Notes  :
//
// Creator        : Lon Waters, SNL
//
// Creation Date  : 12/26/2001
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef N_IO_FUNCTIONBLOCK_H
#define N_IO_FUNCTIONBLOCK_H

#include <string>
#include <vector>

#include <N_IO_fwd.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_UTL_NetlistLocation.h>
#include <N_UTL_Pack.h> 

namespace Xyce {
namespace IO {

class FunctionBlock
{
  friend class Pack<FunctionBlock>;

public:
  FunctionBlock()
  {}

  FunctionBlock(
    const std::string & fileName,
    const TokenVector & parsedInputLine);

  FunctionBlock(
    const std::string & fileName,
    const int lineNumber);

  FunctionBlock(FunctionBlock const& rhsFB);

  ~FunctionBlock()
  {}

  // Public data.
  std::string functionName;

  std::string functionNameAndArgs;

  std::vector<std::string> functionArgs;

  std::string functionBody;

  // Public methods.

  void print();
  // Prints the details of an ParameterBlock to standard out.

private:
  bool extractData(const TokenVector & parsedInputLine);
  //- Extract model data from parsedInputLine using model metadata.

  NetlistLocation netlistLocation_;
};

} // namespace IO
} // namespace Xyce

#endif


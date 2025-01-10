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

//-------------------------------------------------------------------------
//
// Purpose        : Forward declarations
//
// Special Notes  : Forward declaring everything as a class breaks if the implementation of the type changes (like during
//                  templatization)
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_NodeSymbols_h
#define Xyce_N_UTL_NodeSymbols_h

#include <string>
#include <vector>

#include <N_UTL_fwd.h>

namespace Xyce {
namespace Util {

enum SymbolType {SOLUTION_SYMBOL, STATE_SYMBOL, STORE_SYMBOL, EXTERN_SYMBOL, VSRC_SYMBOL, BRANCH_SYMBOL, 
                 NOISE_DEVICE_SYMBOL, NOISE_TYPE_SYMBOL, SYMBOL_TYPE_END};

typedef std::vector<NodeNameMap> SymbolTable;

void addSymbol(SymbolTable &symbol_table, SymbolType symbol_type, int index, const std::string &name);

} // namespace Util
} // namespace Xyce

#endif // Xyce_N_UTL_NodeSymbols_h

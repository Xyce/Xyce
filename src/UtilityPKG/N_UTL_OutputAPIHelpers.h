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
// Purpose        : Provide some string and ParamList manipulating helper
//                  functions
//
// Special Notes  : These are intended for use with the external coupling
//                  output interface, but may have more general utility,
//                  so I'm putting them here.
//
// Creator        : Tom Russo
// Creation Date  : 02/07/2018
//-----------------------------------------------------------------------------

#ifndef N_UTL_OutputAPIHelpers_H
#define N_UTL_OutputAPIHelpers_H

#include <Xyce_config.h>
#include <string>
#include <vector>
#include <N_UTL_Param.h>

namespace Xyce {
namespace Util {

std::string stripWhiteSpace(std::string s);
bool parseFunctionString(const std::string &s, std::string &funName,
                         std::vector<std::string> &funArgs);
bool stringToParamList(const std::string &s, ParamList &paramList);
bool stringsToParamList(const std::vector<std::string> & stringVec,
                        ParamList & paramList,
                        std::vector<bool> & stringValidities);
}
}
#endif


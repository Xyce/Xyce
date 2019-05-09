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
// Purpose        : set version string
// Special Notes  :
//
// Creator        : Eric Rankin
//
// Creation Date  :
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef N_UTL_VERSION_H
#define N_UTL_VERSION_H

#include <string>

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Class         : Version
// Purpose       : Wrapper for retrieving Xyce version string info.
// Special Notes :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
class Version
{

public:

  // get full banner string for Xyce version
  static std::string getFullVersionString();
  // get the maj-min-rev number for Xyce version
  static std::string getShortVersionString();

  static std::string getBuildVariant();

  static std::string getXyceName();
  
  static std::string getCapabilities();

  static std::string getLicense();

private:

  // no need to create objects
  Version();

};

} // namespace Util
} // namespace Xyce

#endif

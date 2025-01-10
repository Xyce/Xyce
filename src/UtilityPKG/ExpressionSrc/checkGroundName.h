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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 10/xx/2023
//-------------------------------------------------------------------------

#ifndef checkGroundName_H
#define checkGroundName_H

// std includes
#include<string>

// Xyce includes
#include <N_UTL_HspiceBools.h>
#include <N_UTL_ExtendedString.h>
#include <N_IO_fwd.h>

namespace Xyce {
namespace Util {

//-------------------------------------------------------------------------
// Function      : checkGroundName
// Purpose       : checks if the name of a voltage node is a synonym for ground
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 9/3/2020
//-------------------------------------------------------------------------
inline bool checkGroundNodeName(const std::string & name)
{
  std::string tmpName = name;
  Xyce::Util::toUpper(tmpName);

  if (tmpName == std::string("0")) { return true; }

  if (Xyce::Util::preprocessFilter.size() >= Xyce::IO::PreprocessType::NUM_PREPROCESS)
  {
  if (Xyce::Util::preprocessFilter[Xyce::IO::PreprocessType::REPLACE_GROUND])
  {
    if (tmpName == std::string("GND")) { return true; }
    if (tmpName == std::string("GND!")) { return true; }
    if (tmpName == std::string("GROUND")) { return true; }
  }
  }

  bool thisIsGround=false;
  if (tmpName.size() > 1)
  {
    std::size_t lastColonIndex = tmpName.find_last_of(":");
    std::size_t last = tmpName.size()-1;
    if (lastColonIndex < last)
    {
      std::string endOfName = tmpName.substr(lastColonIndex,last);
      if (endOfName == ":0") { thisIsGround=true; }
      else if (endOfName == ":GND") { thisIsGround=true; }
      else if (endOfName == ":GND!") { thisIsGround=true; }
      else if (endOfName == ":GROUND") { thisIsGround=true; }
    }
  }

  return thisIsGround;
}

}
}

#endif

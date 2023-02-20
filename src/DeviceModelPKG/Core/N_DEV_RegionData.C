//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Creation Date  : 07/19/06
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <iostream>
#include <N_UTL_Math.h>

// ----------    Xyce Includes  ----------
#include <N_DEV_Const.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_RegionData.h>
#include <N_UTL_ExtendedString.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<RegionData>::ParametricData()
{
  // Set up map for double precision variables:
  addPar("AREA", 1.0e+15, &RegionData::area);
  addPar("XLOC", 0.0, &RegionData::xloc);

  // Set up map for non-double precision variables:
  addPar("NAME", std::string("none"), &RegionData::name);
  addPar("TYPE", std::string("none"), &RegionData::type);
  addPar("FILE", std::string(""), &RegionData::reactionFile);
}

ParametricData<RegionData> &RegionData::getParametricData() {
  static ParametricData<RegionData> parMap;

  return parMap;
}

//-----------------------------------------------------------------------------
// RegionData functions:
//
//-----------------------------------------------------------------------------
// Function      : RegionData::RegionData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/24/06
//-----------------------------------------------------------------------------
RegionData::RegionData ():
  CompositeParam (getParametricData()),
  name("RXNREGION"),
  outName("RXNREGION"),
  type("JUNCTION"),
  reactionFile("reaction_spec_full"),
  area(1.0e-4),
  xloc(1.83e-4),
  doNothing(false)
{}

//-----------------------------------------------------------------------------
// Function      : RegionData::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/24/06
//-----------------------------------------------------------------------------
void RegionData::processParams()
{
  ParameterMap::const_iterator p_i = getParameterMap().find("TYPE");
  const Descriptor &p = *(*p_i).second;
  ExtendedString tmp = p.value<std::string>(*this);
  p.value<std::string>(*this) = tmp.toLower();
}

//-----------------------------------------------------------------------------
// Function      : RegionData::operator<<
// Purpose       : "<<" operator
// Special Notes : doesn't print everything.
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/23/05
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const RegionData & rd)
{
  os << " Region Data: name = " << rd.name <<
    " x=" << rd.xloc <<
    " reaction file = " << rd.reactionFile <<
    " type = " << rd.type <<
    std::endl;

  return os;
}

} // namespace Device
} // namespace Xyce

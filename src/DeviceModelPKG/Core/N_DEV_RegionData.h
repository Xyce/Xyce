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

//-----------------------------------------------------------------------------
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
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_RegionData_h
#define Xyce_N_DEV_RegionData_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_CompositeParam.h>

// ---------- Forward Declarations -------

//-----------------------------------------------------------------------------
// Class         : N_DEV_RegionData
// Purpose       :
// Special Notes : This class is intended to be a class for passing data into
//                 the Rxn Region constructor.
//
// Creator       : Eric Keiter, SNL
// Creation Date : 4/17/11
//-----------------------------------------------------------------------------

namespace Xyce {
namespace Device {

class RegionData : public CompositeParam
{
  friend class ParametricData<RegionData>;
  friend std::ostream & operator<<(std::ostream & os, const RegionData & rd);

public:
  static ParametricData<RegionData> &getParametricData();

  RegionData ();

  void processParams ();

public:
  std::string name;
  std::string outName;
  std::string type;
  std::string reactionFile;
  double area;
  double xloc;

  bool doNothing; // set to true if reaction set should be ignored.
};

} // namespace Device
} // namespace Xyce

#endif


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
// Purpose        : This file contains the PDE device instance base class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/15/01
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_SpecieSource_h
#define Xyce_N_DEV_SpecieSource_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_Const.h>
#include <N_DEV_CompositeParam.h>

namespace Xyce {
namespace Device {


//-----------------------------------------------------------------------------
// Class         : SpecieSource
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
class SpecieSource : public CompositeParam
{
  friend class ParametricData<SpecieSource>;

public:
  static ParametricData<SpecieSource> &getParametricData();

  SpecieSource ();
  bool processParam (Param & ndParam, std::string & param, DevicePDEInstance & di);
  void processParams ();

public:
  std::string name;         // this is also the index into the map.
};

// inline functions
//-----------------------------------------------------------------------------
// Function      : SpecieSource::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
inline std::ostream & operator<<(std::ostream & os, const SpecieSource & ds)
{
  os << ds.name << ":\n";
  os << std::endl;
  return os;
}

} // namespace Device
} // namespace Xyce

#endif


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

//----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//----------------------------------------------------------------------------

#include <Xyce_config.h>
#include <N_DEV_PlaceHolder.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : PlaceHolder::factory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Device * PlaceHolder::factory(SolverState & solver_state, DeviceOptions & device_options)
{
  static PlaceHolder placeHolder(solver_state, device_options);

  return &placeHolder;
}

//-----------------------------------------------------------------------------
// Function      : PlaceHolder::PlaceHolder
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
PlaceHolder::PlaceHolder(SolverState & ss1, DeviceOptions & do1)
  : Device(ss1,do1)
{}

//-----------------------------------------------------------------------------
// Function      : PlaceHolder::~PlaceHolder
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
PlaceHolder::~PlaceHolder()
{}

} // namespace Device
} // namespace Xyce

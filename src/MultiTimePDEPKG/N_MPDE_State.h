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
// Purpose       : State info for MPDE Solver
//
// Special Notes : Useful for outside objects to lookup info
//
// Creator       : Robert Hoekstra, 9233, Computational Sciences
//
// Creation Date : 3/23/04
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_MPDE_STATE_H
#define Xyce_MPDE_STATE_H

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

// ---------- Forward Declarations ----------

// ---------- Enum Definitions ----------

//-----------------------------------------------------------------------------
// Class         : N_MPDE_State
// Purpose       : 
// Special Notes :
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 3/23/04
//-----------------------------------------------------------------------------
struct N_MPDE_State
{
  // Default constructor
  N_MPDE_State()
  : fastTime(0.0)
  {}

  //Fast Time for MPDE driving source
  double fastTime;
};

#endif //Xyce_MPDE_STATE_H

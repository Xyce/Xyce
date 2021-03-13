//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Purpose        : This file contains some enumerated types which are needed 
//                  both inside and outside the expression interface.
//
// Special Notes  :
//
// Creator        : 
//
// Creation Date  : 
//
//-----------------------------------------------------------------------------

#ifndef astRandEnum_H
#define astRandEnum_H

namespace Xyce {
namespace Util {
// ---------- Enumerations ----------

enum astRandTypes
{
  AST_AGAUSS,          //  0
  AST_GAUSS,           //  1
  AST_AUNIF,           //  2
  AST_UNIF,            //  3
  AST_RAND,            //  4
  AST_LIMIT            //  5
};

} // Util
} // Xyce

#endif // astRandEnum_H

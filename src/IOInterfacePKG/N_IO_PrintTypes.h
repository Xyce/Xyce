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

//-------------------------------------------------------------------------
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : Eric Keiter, SNL
//
// Creation Date  : 2016/04/7 18:01:27
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_IO_PrintTypes_h
#define Xyce_N_IO_PrintTypes_h


namespace Xyce {
namespace IO {

namespace PrintType {
  enum PrintType {NONE, DC, TRAN, AC, AC_IC, HB, HB_TD, HB_FD, HB_IC, HB_STARTUP, HOMOTOPY, MPDE, MPDE_IC, MPDE_STARTUP, RAW_OVERRIDE, SENS, TRANADJOINT, NOISE, SPARAM, ES, PCE};
}

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_PrintTypes_h


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
// Purpose        : Class for re-reading Xyce file output, of simulation results,
//                  that is in HDF5 format.
//
// Special Notes  :
//
// Creator        : Richard Schiek, Electrical Systems Modeling, Sandia National Laboratories
//
// Creation Date  : 12/06/12
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputHDF5_h
#define Xyce_N_IO_OutputHDF5_h

// ----------   Standard Includes   ----------

// ----------   Xyce Includes   ----------
#include<N_IO_OutputFileBase.h>

// ---------- Forward Declarations ----------


namespace Xyce {
namespace IO {

class OutputHDF5 : public OutputFileBase
{
  public:
  OutputHDF5();
  ~OutputHDF5();

};

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputHDF5_h

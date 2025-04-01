//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Purpose        : This is the base class for mesh processing/ownership.
//
// Special Notes  : Classes derived off of this one will be developed
//                  as needed.
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/21/02
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_PDE_Mesh_Container_h
#define Xyce_N_DEV_PDE_Mesh_Container_h

#include <string>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : PDEMeshContainer
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/21/02
//-----------------------------------------------------------------------------
class PDEMeshContainer
{
public:
  PDEMeshContainer ();
  ~PDEMeshContainer ();

  virtual bool initializeMesh (const std::string & meshFileName);
};

} // namespace Device
} // namespace Xyce

#endif


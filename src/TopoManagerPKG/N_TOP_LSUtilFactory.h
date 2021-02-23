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
// Purpose        : Specification file for the interface for creating linear
//                  algebra objects using the GoF Abstract Factory design
//                  pattern.  It must be used with a compatible linear algebra
//                  package (e.g., Trilinos/Petra to provide the concrete
//                  implementations.
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/20/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TOP_LSUtilFactory_h
#define Xyce_N_TOP_LSUtilFactory_h

// ---------- Standard Includes ----------

#include <N_IO_fwd.h>
#include <N_TOP_fwd.h>
#include <N_PDS_fwd.h>
#include <N_LAS_fwd.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : LSUtilFactory
// Purpose       : Provides an interface for creating linear system query objects
// Special Notes :
// Creator       : Heidi K. Thornquist, SNL, Parallel Compuational Sciences
// Creation Date : 02/23/15
//-----------------------------------------------------------------------------
class LSUtilFactory
{

 public:

  // Creates a new LAS vector
  static Linear::QueryUtil * newLSUtil( Topology &topology, Parallel::Manager &pds_manager );

 private:

  // The Topo Algebra Services default factory
  LSUtilFactory();

};

} // namespace Topo
} // namespace Xyce

#endif // Xyce_N_TOP_LSUtilFactory_h

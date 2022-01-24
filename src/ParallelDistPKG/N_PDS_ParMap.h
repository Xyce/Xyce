//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Purpose        : Specification file for abstract base class for the parallel
//                  map data and functions.
//
// Special Notes  : Part of a GoF Abstract Factory.
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/08/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_PDS_ParMap_h
#define Xyce_N_PDS_ParMap_h

// ---------- Standard Includes ----------

#include <iostream>

// ----------   Xyce Includes   ----------

#include <N_PDS_fwd.h>

namespace Xyce {
namespace Parallel {

//-----------------------------------------------------------------------------
// Class         : ParMap
// Purpose       : Abstract base class for the parallel map data and functions.
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 03/08/00
//-----------------------------------------------------------------------------
class ParMap
{

public:

  ParMap(const Xyce::Parallel::Communicator& aComm)
  : pdsComm_(aComm)
  {}

  // Destructor
  virtual ~ParMap() {}

  // Clone method
  virtual ParMap* clone() const = 0;

  // Number of global "entities" represented as vertices in the graph. These
  // may be, for example, equations for the linear algebra quantities or
  // devices/nodes for the circuit graph.
  virtual int numGlobalEntities() const = 0;

  // Number of local (on processor) "entities".
  virtual int numLocalEntities() const = 0;

  // Indexing base (0 or 1) used for the maps.
  virtual int indexBase() const = 0;

  // Minimum globally-numbered identifier on this processor.
  virtual int minMyGlobalEntity() const = 0;

  // Maximum globally-numbered identifier on this processor.
  virtual int maxMyGlobalEntity() const = 0;

  // Maximum globally-numbered identifier.
  virtual int maxGlobalEntity() const = 0;

  // Comm object
  const Xyce::Parallel::Communicator& pdsComm() const { return pdsComm_; }

  // dereference global index to get local index
  virtual int globalToLocalIndex(int global_index) const = 0;

  // dereference local index to get global index
  virtual int localToGlobalIndex(int local_index) const = 0;

  // Output the map to a file
  virtual void writeToFile(const char * filename) const = 0;

  virtual void print(std::ostream &os) const {}

protected:

  // Comm object.
  const Xyce::Parallel::Communicator &          pdsComm_;
};

} // namespace Parallel
} // namespace Xyce

#endif

//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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

#ifndef Xyce_N_PDS_EpetraParMap_h
#define Xyce_N_PDS_EpetraParMap_h

// ---------- Standard Includes ----------

#include <ostream>

// ----------   Xyce Includes   ----------

#include <N_PDS_fwd.h>
#include <N_PDS_ParMap.h>

class Epetra_Map;

namespace Xyce {
namespace Parallel {

//-----------------------------------------------------------------------------
// Class         : EpetraParMap
// Purpose       : Epetra implementation of the parallel map data and functions.
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 03/08/00
//-----------------------------------------------------------------------------
class EpetraParMap : public ParMap
{

public:
  // Constructor which takes a Epetra map.
  EpetraParMap( const Epetra_Map * pMap,
                const Communicator & aComm,
                bool mapOwned = false );

  // Destructor
  virtual ~EpetraParMap();

  // Clone method
  EpetraParMap* clone() const;

private:

  // Copy constructor (private).
  EpetraParMap(const ParMap & right);

  // Assignment operator (private).
  EpetraParMap & operator=(const ParMap & right);

  // Equality operator (private).
  bool operator==(const EpetraParMap & right) const;

  // Non-equality operator (private).
  bool operator!=(const EpetraParMap & right) const;

public:

  // Number of global "entities" represented as vertices in the graph. These
  // may be, for example, equations for the linear algebra quantities or
  // devices/nodes for the circuit graph.
  int numGlobalEntities() const;

  // Number of local (on processor) "entities".
  int numLocalEntities() const;

  // Indexing base (0 or 1) used for the maps.
  int indexBase() const;

  // Minimum globally-numbered identifier on this processor.
  int minMyGlobalEntity() const;

  // Maximum globally-numbered identifier on this processor.
  int maxMyGlobalEntity() const;

  // Maximum globally-numbered identifier.
  int maxGlobalEntity() const;

  // dereference global index to get local index
  int globalToLocalIndex(int global_index) const;

  // dereference local index to get global index
  int localToGlobalIndex(int local_index) const;

  // Accessor functions (overridden in derived classes) for the pointer to the
  // library map object.
  const Epetra_Map * petraMap() const { return petraMap_; }

  void writeToFile(const char * filename) const;

  void print(std::ostream &os) const;

private:

  // Pointer to Petra map object.
  const Epetra_Map * petraMap_;
  bool mapOwned_;
};

} // namespace Parallel
} // namespace Xyce

#endif

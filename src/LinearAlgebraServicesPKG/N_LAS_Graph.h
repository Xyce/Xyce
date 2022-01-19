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
// Purpose        : Abstract interface to linear solver type.
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL
//
// Creation Date  : 06/19/20
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_Graph_h
#define Xyce_N_LAS_Graph_h

#include <iostream>

#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : Graph
// Purpose       : Abstract interface to a graph object
// Special Notes : This is necessary to define the non-zero pattern of the matrix.
// Creator       : Heidi Thornquist, SNL
// Creation Date : 06/19/20
//-----------------------------------------------------------------------------
class Graph
{

public:

  // Default constructor
  Graph() {}

  // Destructor
  virtual ~Graph() {}

  // Clone this graph
  virtual Graph* cloneCopy() const = 0;

  // Create new graph exporting values from this one.
  virtual Graph* exportGraph( const Parallel::ParMap& map ) const = 0;

  // Get the base index for this graph
  virtual int indexBase() const = 0;

  // Get the maximum number of indices for any row on this processor.
  virtual int maxNumIndices() const = 0;

  // Get the number of rows on this processor.
  virtual int numLocalEntities() const = 0;

  // Get the number of nonzero entries on this processor.
  virtual int numLocalNonzeros() const = 0;

  // Convert a local ID to the global ID using the row map
  virtual int localToGlobalRowIndex(int localIndex) const = 0;

  // Convert a local ID to the global ID using the column map
  virtual int localToGlobalColIndex(int localIndex) const = 0;

  // Convert a global ID to the local ID using the row map
  virtual int globalToLocalRowIndex(int globalIndex) const = 0;

  // Convert a global ID to the local ID using the column map
  virtual int globalToLocalColIndex(int globalIndex) const = 0;

  // Get a copy of the row on this processor using local IDs
  virtual void extractLocalRowCopy(int localRow, int length, int& numIndices, int* indices) const = 0;

  // Get a copy of the row on this processor using global IDs
  virtual void extractGlobalRowCopy(int globalRow, int length, int& numIndices, int* indices) const = 0;

  // Get a pointer to the local row with local IDs
  virtual void extractLocalRowView(int localRow, int& numIndices, int*& indices) const = 0;

  // Insert indices into the graph with global IDs
  virtual void insertGlobalIndices(int globalRow, int numIndices, int* indices) = 0;

  // Accumulate off processor fill contributions if necessary
  virtual void fillComplete() = 0;

  // Accumulate off processor fill contributions using rowMap/colMap objects
  virtual void fillComplete( Parallel::ParMap& rowMap, Parallel::ParMap& colMap ) = 0;

  // Print the underlying object.
  virtual void print(std::ostream &os) const = 0;

private:

  // Copy constructor
  Graph( const Graph& graph ) {}
};

} // namespace Linear
} // namespace Xyce

#endif

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

#ifndef Xyce_N_LAS_EpetraGraph_h
#define Xyce_N_LAS_EpetraGraph_h

#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>

#include <N_LAS_Graph.h>

#include <Teuchos_RCP.hpp>
#include <Epetra_CrsGraph.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : EpetraGraph
// Purpose       : Epetra interface to a graph object
// Special Notes : This is necessary to define the non-zero pattern of the matrix.
// Creator       : Heidi Thornquist, SNL
// Creation Date : 12/23/20
//-----------------------------------------------------------------------------
class EpetraGraph : public Graph
{

public:

  EpetraGraph(const Parallel::ParMap & solution_overlap,
              const Parallel::ParMap & solution_overlap_ground,
              const std::vector<int>& numIndicesPerRow,
              const std::vector<std::vector<int> >& rcData);

  // Basic constructor with map and number of entries per row
  EpetraGraph( const Parallel::ParMap & map, const std::vector<int>& numIndicesPerRow );

  // Basic constructor with map and maximum number of entries per row
  EpetraGraph( const Parallel::ParMap & map, int maxIndicesPerRow );

  // Simple constructor using Epetra_CrsGraph
  EpetraGraph( const Teuchos::RCP<Epetra_CrsGraph>& graph );

  // Destructor
  virtual ~EpetraGraph() {}

  // Clone this graph
  EpetraGraph* cloneCopy() const;

  // Create new graph exporting values from this one.
  EpetraGraph* exportGraph( const Parallel::ParMap& map ) const;

  // Get the base index for this graph
  int indexBase() const
  { return epetraGraph_->IndexBase(); }

  // Get the maximum number of indices for any row on this processor.
  int maxNumIndices() const
  { return epetraGraph_->MaxNumIndices(); }

  // Get the number of rows on this processor.
  int numLocalEntities() const 
  { return epetraGraph_->NumMyRows(); } 

  // Get the number of nonzero entries on this processor.
  int numLocalNonzeros() const
  { return epetraGraph_->NumMyNonzeros(); }

  int localToGlobalRowIndex(int localIndex) const 
  { return epetraGraph_->GRID( localIndex ); }

  int localToGlobalColIndex(int localIndex) const 
  { return epetraGraph_->GCID( localIndex ); }

  int globalToLocalRowIndex(int globalIndex) const
  { return epetraGraph_->LRID( globalIndex ); }

  int globalToLocalColIndex(int globalIndex) const
  { return epetraGraph_->LCID( globalIndex ); }

  void extractLocalRowCopy(int localRow, int length, int& numIndices, int* indices) const
  { epetraGraph_->ExtractMyRowCopy( localRow, length, numIndices, indices ); }

  void extractGlobalRowCopy(int globalRow, int length, int& numIndices, int* indices) const
  { epetraGraph_->ExtractGlobalRowCopy( globalRow, length, numIndices, indices ); }

  void extractLocalRowView(int localRow, int& numIndices, int*& indices) const
  { epetraGraph_->ExtractMyRowView( localRow, numIndices, indices ); }

  //Accumulate off processor fill contributions if necessary
  void fillComplete()
  {
    epetraGraph_->FillComplete();
    epetraGraph_->OptimizeStorage(); 
  }

  void fillComplete( Parallel::ParMap& rowMap, Parallel::ParMap& colMap );

  const Teuchos::RCP<Epetra_CrsGraph>& epetraObj() const { return epetraGraph_; }

  // Print the underlying object.
  void print(std::ostream &os) const
  {  epetraGraph_->Print( os ); }

private:

  // Copy constructor
  EpetraGraph( const EpetraGraph& graph );

  Teuchos::RCP<Epetra_CrsGraph> epetraGraph_;

};

} // namespace Linear
} // namespace Xyce

#endif

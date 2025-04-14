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
// Purpose        : Implementation file for the Abstract interface to the
//                  vector types (RDP, RSP, CDP or CSP).
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Computational Sciences
//
// Creation Date  : 10/13/00
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_FeatureTest.h>
#include <N_LAS_EpetraGraph.h>
#include <N_PDS_EpetraParMap.h>

// ---------  Other Includes  -----------

#include <Epetra_Map.h>
#include <Epetra_Export.h>

namespace Xyce {
namespace Linear {

EpetraGraph::EpetraGraph(  const Parallel::ParMap & solution_overlap,
              const Parallel::ParMap & solution_overlap_ground,
              const std::vector<int>& numIndicesPerRow,
              const std::vector<std::vector<int> >& rcData)
  {
    const Epetra_Map* epetraMap = dynamic_cast<const Parallel::EpetraParMap&>(solution_overlap).petraMap();
    epetraGraph_ = Teuchos::rcp( new Epetra_CrsGraph( Copy, *epetraMap, numIndicesPerRow.data() ) );

    int numLocalRows_Overlap = rcData.size();
    for(int i = 0; i < numLocalRows_Overlap; ++i)
    {
      std::vector<int>& rcData_i = const_cast<std::vector<int> &>(rcData[i]);
      if( solution_overlap_ground.localToGlobalIndex(i) != -1 && numIndicesPerRow[i] )
      {
        if( rcData_i[0] == -1 )
          epetraGraph_->InsertGlobalIndices( solution_overlap.localToGlobalIndex(i), numIndicesPerRow[i]-1, &rcData_i[1] );
        else
          epetraGraph_->InsertGlobalIndices( solution_overlap.localToGlobalIndex(i), numIndicesPerRow[i], &rcData_i[0] );
      }
    }
    
    epetraGraph_->FillComplete();
    epetraGraph_->OptimizeStorage();
  }

  EpetraGraph::EpetraGraph( const Parallel::ParMap & map, const std::vector<int>& numIndicesPerRow )
  {
    const Epetra_Map* epetraMap = dynamic_cast<const Parallel::EpetraParMap&>(map).petraMap();
    epetraGraph_ = Teuchos::rcp( new Epetra_CrsGraph( Copy, *epetraMap, &numIndicesPerRow[0] ) );
  }

  // Basic constructor with map and maximum number of entries per row
  EpetraGraph::EpetraGraph( const Parallel::ParMap & map, int maxIndicesPerRow )
  {
    const Epetra_Map* epetraMap = dynamic_cast<const Parallel::EpetraParMap&>(map).petraMap();
    epetraGraph_ = Teuchos::rcp( new Epetra_CrsGraph( Copy, *epetraMap, maxIndicesPerRow ) );
  }

  EpetraGraph::EpetraGraph( const Teuchos::RCP<Epetra_CrsGraph>& graph )
  : epetraGraph_(graph)
  {}

  EpetraGraph::EpetraGraph( const EpetraGraph& graph )
  : epetraGraph_(graph.epetraGraph_)
  {}

  EpetraGraph* EpetraGraph::cloneCopy() const
  {
    return( new EpetraGraph( *this ) ); 
  }

  EpetraGraph* EpetraGraph::exportGraph( const Parallel::ParMap& map ) const
  {
    const Epetra_Map* exportMap = dynamic_cast<const Parallel::EpetraParMap&>(map).petraMap();

    Epetra_Export exporter( epetraGraph_->Map(), *exportMap );
    Epetra_CrsGraph * newGraph = new Epetra_CrsGraph( Copy, *exportMap, 0 );
    newGraph->Export( *epetraGraph_, exporter, Add );
    newGraph->FillComplete();
    newGraph->OptimizeStorage();
    EpetraGraph* retGraph = new EpetraGraph( Teuchos::rcp( newGraph ) );

    return retGraph;
  }

  void EpetraGraph::fillComplete( Parallel::ParMap& rowMap, Parallel::ParMap& colMap )
  {
    const Epetra_Map* rMap = dynamic_cast<Parallel::EpetraParMap&>(rowMap).petraMap();
    const Epetra_Map* cMap = dynamic_cast<Parallel::EpetraParMap&>(colMap).petraMap();

    epetraGraph_->FillComplete( *rMap, *cMap );
    epetraGraph_->OptimizeStorage();
  }

} // namespace Linear
} // namespace Xyce

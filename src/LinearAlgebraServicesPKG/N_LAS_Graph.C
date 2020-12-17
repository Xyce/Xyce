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

#include <N_LAS_Graph.h>
#include <N_UTL_FeatureTest.h>
#include <N_PDS_EpetraParMap.h>

// ---------  Other Includes  -----------

#include <Epetra_Map.h>
#include <Epetra_Export.h>

namespace Xyce {
namespace Linear {

  Graph::Graph( Parallel::ParMap & map, const std::vector<int>& numIndicesPerRow )
  {
    Epetra_Map* epetraMap = dynamic_cast<Parallel::EpetraParMap&>(map).petraMap();
    epetraGraph_ = Teuchos::rcp( new Epetra_CrsGraph( Copy, *epetraMap, &numIndicesPerRow[0] ) );
  }

  // Basic constructor with map and maximum number of entries per row
  Graph::Graph( Parallel::ParMap & map, int maxIndicesPerRow )
  {
    Epetra_Map* epetraMap = dynamic_cast<Parallel::EpetraParMap&>(map).petraMap();
    epetraGraph_ = Teuchos::rcp( new Epetra_CrsGraph( Copy, *epetraMap, maxIndicesPerRow ) );
  }

  Graph::Graph( const Teuchos::RCP<Epetra_CrsGraph>& graph )
  : epetraGraph_(graph)
  {}

  Graph::Graph( const Graph& graph )
  : epetraGraph_(graph.epetraGraph_)
  {}

  Graph* Graph::cloneCopy() const
  {
    return( new Graph( *this ) ); 
  }

  Graph* Graph::exportGraph( Parallel::ParMap& map ) const
  {
    Epetra_Map* exportMap = dynamic_cast<Parallel::EpetraParMap&>(map).petraMap();

    Epetra_Export exporter( epetraGraph_->Map(), *exportMap );
    Epetra_CrsGraph * newGraph = new Epetra_CrsGraph( Copy, *exportMap, 0 );
    newGraph->Export( *epetraGraph_, exporter, Add );
    newGraph->FillComplete();
    newGraph->OptimizeStorage();
    Graph* retGraph = new Graph( Teuchos::rcp( newGraph ) );

    return retGraph;
  }

  void Graph::fillComplete( Parallel::ParMap& rowMap, Parallel::ParMap& colMap )
  {
    Epetra_Map* rMap = dynamic_cast<Parallel::EpetraParMap&>(rowMap).petraMap();
    Epetra_Map* cMap = dynamic_cast<Parallel::EpetraParMap&>(colMap).petraMap();

    epetraGraph_->FillComplete( *rMap, *cMap );
    epetraGraph_->OptimizeStorage();
  }

} // namespace Linear
} // namespace Xyce

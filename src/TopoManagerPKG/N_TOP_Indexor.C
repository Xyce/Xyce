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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/12/02
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <map>

#include <N_TOP_Indexor.h>

#include <N_PDS_Manager.h>
#include <N_PDS_Comm.h>
#include <N_PDS_ParMap.h>
#include <N_LAS_Graph.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : Indexor::globalToLocal
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/12/02
//-----------------------------------------------------------------------------
bool
Indexor::globalToLocal(
  int                   graph_id,
  std::vector<int> &    ids )
{
  Parallel::ParMap *map = pdsMgr_.getParallelMap(graph_id);

  for (unsigned int i = 0; i < ids.size(); ++i)
    ids[i] = map->globalToLocalIndex( ids[i] );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Indexor::localToGlobal
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/28/10
//-----------------------------------------------------------------------------
bool
Indexor::localToGlobal(
  int                   graph_id,
  std::vector<int> &    ids )
{
  Parallel::ParMap *map = pdsMgr_.getParallelMap(graph_id);

  for( unsigned int i = 0; i < ids.size(); ++i )
    ids[i] = map->localToGlobalIndex( ids[i] );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Indexor::setupAcceleratedMatrixIndexing
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/23/02
//-----------------------------------------------------------------------------
bool
Indexor::setupAcceleratedMatrixIndexing(
  int           graph_id)
{
  const Linear::Graph * graph = pdsMgr_.getMatrixGraph( graph_id );

  int NumRows = graph->numLocalEntities();
  matrixIndexMap_.clear();
  matrixIndexMap_.resize( NumRows );

  int NumElements;
  int * Elements;
  for( int i = 0; i < NumRows; ++i )
  {
    graph->extractLocalRowView( i, NumElements, Elements );
    for( int j = 0; j < NumElements; ++j ) matrixIndexMap_[i][ Elements[j] ] = j;
  }

  accelMatrixIndex_ = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Indexor::deleteAcceleratedMatrixIndexing
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/23/02
//-----------------------------------------------------------------------------
bool Indexor::deleteAcceleratedMatrixIndexing()
{
  matrixIndexMap_.clear();
  accelMatrixIndex_ = false;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Indexor::matrixGlobalToLocal
// Purpose       : This method converts the GIDs in device stamps to LIDs
// Special Notes : If the GID is ground it will also be -1 as the LID
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/23/02
//-----------------------------------------------------------------------------
bool
Indexor::matrixGlobalToLocal(
  int                                   graph_id,
  const std::vector<int> &              gids,
  std::vector< std::vector<int> > &     stamp )
{
  const Linear::Graph* graph = pdsMgr_.getMatrixGraph( graph_id );

  int numRows = stamp.size();

  int numElements;
  int * elements;

  if( accelMatrixIndex_ )
  {
    for( int i = 0; i < numRows; ++i )
    {
      int numCols = stamp[i].size();
      if (gids[i] != -1)
      {
        int rowLID = graph->globalToLocalRowIndex(gids[i]);
        for( int j = 0; j < numCols; ++j )
        {
          int lid = graph->globalToLocalColIndex(stamp[i][j]);
          if (stamp[i][j] != -1)
            stamp[i][j] = matrixIndexMap_[rowLID][lid];
          else
            stamp[i][j] = -1;
        }
      }
      else
      {
        for( int j = 0; j < numCols; ++j )
          stamp[i][j] = -1;
      }
    }
  }
  else
  {
    for( int i = 0; i < numRows; ++i )
    {
      int numCols = stamp[i].size();
      if (gids[i] != -1)
      { 
        graph->extractLocalRowView( graph->globalToLocalRowIndex(gids[i]), numElements, elements );

        std::map<int,int> indexToOffsetMap;
        for( int j = 0; j < numElements; ++j ) 
          indexToOffsetMap[ elements[j] ] = j;

        for( int j = 0; j < numCols; ++j )
        {
          int lid = graph->globalToLocalColIndex(stamp[i][j]);
          stamp[i][j] = indexToOffsetMap[lid];
        }
      }
      else
      {
        for( int j = 0; j < numCols; ++j )
          stamp[i][j] = -1;
      } 
    }
  }

  return true;
}

} // namespace Topo
} // namespace Xyce

//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose        : Implementation file for abstract base class for the
//                  parallel map data and functions.
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
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>
#include <N_UTL_FeatureTest.h>
#include <N_ERH_Message.h>

using Xyce::DEBUG_PARALLEL;

// ----------   Other Includes   ----------

#include <Epetra_Map.h>

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::N_PDS_ParMap
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/2/00
//-----------------------------------------------------------------------------
N_PDS_ParMap::N_PDS_ParMap(
  int &                         numGlobalEntities,
  int                           numLocalEntities,
  const std::vector<int> &      lbMap,
  const int                     index_base,
  N_PDS_Comm &                  aComm)
  : petraMap_(0),
    mapOwned_(true),
    pdsComm_(aComm)
{
  const int * mArray = lbMap.empty() ? 0 : static_cast<const int *>(&lbMap[0]);

  // fix for empty maps
  int nGE = std::max( -1, numGlobalEntities );
  int nLE = std::max( 0, numLocalEntities );
  // Call the Petra constructor for the true Petra map.
  petraMap_ = new Epetra_Map( nGE,
                              nLE,
                              mArray,
                              index_base,
                              *pdsComm_.petraComm());

  if (DEBUG_PARALLEL)
    std::cout << "New Petra Map: " << numGlobalEntities << " " << numLocalEntities << std::endl
              << "  " << petraMap_->NumMyElements() << " " << petraMap_->NumGlobalElements() << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::N_PDS_ParMap
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/2/00
//-----------------------------------------------------------------------------
N_PDS_ParMap::N_PDS_ParMap(
  int &         numGlobalEntities,
  int           numLocalEntities,
  const int     index_base,
  N_PDS_Comm &  aComm)
  : petraMap_(0),
    mapOwned_(true),
    pdsComm_(aComm)
{
  // Call the Petra constructor for the true Petra map.
  petraMap_ = new Epetra_Map( numGlobalEntities,
                              numLocalEntities,
                              index_base,
                              *pdsComm_.petraComm());

  if (DEBUG_PARALLEL)
    std::cout << "New Petra Map: " << numGlobalEntities << " " << numLocalEntities << std::endl
              << "  " << petraMap_->NumMyElements() << " " << petraMap_->NumGlobalElements() << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::N_PDS_ParMap
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/2/00
//-----------------------------------------------------------------------------
N_PDS_ParMap::N_PDS_ParMap(
  Epetra_Map *          map,
  N_PDS_Comm &          aComm )
  : petraMap_(map),
    mapOwned_(false),
    pdsComm_(aComm)
{}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::~N_PDS_ParMap
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
N_PDS_ParMap::~N_PDS_ParMap()
{
  if (mapOwned_)
    delete petraMap_;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::numGlobalEntities
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int N_PDS_ParMap::numGlobalEntities() const
{
  return petraMap_->NumGlobalElements();
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::numLocalEntities
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int N_PDS_ParMap::numLocalEntities() const
{
  return petraMap_->NumMyElements();
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::indexBase
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int N_PDS_ParMap::indexBase() const
{
  return petraMap_->IndexBase();
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::maxGlobalEntity
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int N_PDS_ParMap::maxGlobalEntity() const
{
  return petraMap_->MaxAllGID();
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::petraBlockMap
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
Epetra_BlockMap * N_PDS_ParMap::petraBlockMap()
{
  return dynamic_cast<Epetra_BlockMap*>(petraMap_);
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::globalToLocalIndex
// Purpose       : dereference Global to Local Index
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/6/00
//-----------------------------------------------------------------------------
int N_PDS_ParMap::globalToLocalIndex(int global_index) const
{
  return petraMap_->LID(global_index);
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::localToGlobalIndex
// Purpose       : dereference Local to Global Index
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/10/06
//-----------------------------------------------------------------------------
int N_PDS_ParMap::localToGlobalIndex(int local_index) const
{
  return petraMap_->GID(local_index);
}

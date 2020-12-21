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

#include <N_PDS_Comm.h>
#include <N_PDS_EpetraParMap.h>

// ----------   Other Includes   ----------

#include <Epetra_Map.h>
#include <EpetraExt_BlockMapOut.h>

namespace Xyce {
namespace Parallel {

//-----------------------------------------------------------------------------
// Function      : EpetraParMap::EpetraParMap
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/2/00
//-----------------------------------------------------------------------------
EpetraParMap::EpetraParMap(
  Epetra_Map *          map,
  Communicator &        aComm,
  bool                  mapOwned )
  : ParMap(aComm),
    petraMap_(map),
    mapOwned_(mapOwned)
{}

//-----------------------------------------------------------------------------
// Function      : EpetraParMap::~EpetraParMap
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
EpetraParMap::~EpetraParMap()
{
  if (mapOwned_)
  {
    delete petraMap_;
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraParMap::clone
// Purpose       : Create a copy of the map 
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 12/21/21
//-----------------------------------------------------------------------------
EpetraParMap* EpetraParMap::clone() const
{
  return new EpetraParMap( new Epetra_Map( *petraMap_ ), pdsComm_, true );
}

//-----------------------------------------------------------------------------
// Function      : EpetraParMap::numGlobalEntities
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int EpetraParMap::numGlobalEntities() const
{
  return petraMap_->NumGlobalElements();
}

//-----------------------------------------------------------------------------
// Function      : EpetraParMap::numLocalEntities
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int EpetraParMap::numLocalEntities() const
{
  return petraMap_->NumMyElements();
}

//-----------------------------------------------------------------------------
// Function      : EpetraParMap::indexBase
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int EpetraParMap::indexBase() const
{
  return petraMap_->IndexBase();
}


//-----------------------------------------------------------------------------
// Function      : EpetraParMap::minMyGlobalEntity
// Purpose       : Minimum globally-numbered identifier on this processor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int EpetraParMap::minMyGlobalEntity() const
{
  return petraMap_->MinMyGID();
}
 
//-----------------------------------------------------------------------------
// Function      : EpetraParMap::maxMyGlobalEntity
// Purpose       : Maximum globally-numbered identifier on this processor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int EpetraParMap::maxMyGlobalEntity() const
{
  return petraMap_->MaxMyGID();
}

//-----------------------------------------------------------------------------
// Function      : EpetraParMap::maxGlobalEntity
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int EpetraParMap::maxGlobalEntity() const
{
  return petraMap_->MaxAllGID();
}

//-----------------------------------------------------------------------------
// Function      : EpetraParMap::globalToLocalIndex
// Purpose       : dereference Global to Local Index
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/6/00
//-----------------------------------------------------------------------------
int EpetraParMap::globalToLocalIndex(int global_index) const
{
  return petraMap_->LID(global_index);
}

//-----------------------------------------------------------------------------
// Function      : EpetraParMap::localToGlobalIndex
// Purpose       : dereference Local to Global Index
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/10/06
//-----------------------------------------------------------------------------
int EpetraParMap::localToGlobalIndex(int local_index) const
{
  return petraMap_->GID(local_index);
}

//-----------------------------------------------------------------------------
// Function      : EpetraParMap::writeToFile
// Purpose       : write out map
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist 
// Creation Date : 12/14/20
//-----------------------------------------------------------------------------
void EpetraParMap::writeToFile(const char * filename) const
{
  EpetraExt::BlockMapToMatrixMarketFile( filename, *petraMap_ );
}

//-----------------------------------------------------------------------------
// Function      : EpetraParMap::print
// Purpose       : print map
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist 
// Creation Date : 09/30/20
//-----------------------------------------------------------------------------
void EpetraParMap::print(std::ostream &os) const
{
  petraMap_->Print(os);
}

} // namespace Parallel
} // namespace Xyce


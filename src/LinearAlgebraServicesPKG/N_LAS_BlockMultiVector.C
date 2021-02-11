//----------------------------------------------------------------------
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

//-------------------------------------------------------------------------
//
// Purpose        : Implementation file for Block MultiVector
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Computational Sciences
//
// Creation Date  : 3/13/04
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------


#include <N_LAS_BlockMultiVector.h>
#include <N_PDS_ParMap.h>
#include <N_PDS_EpetraParMap.h>
#include <N_PDS_Comm.h>

#include <N_LAS_BlockSystemHelpers.h>

// ---------  Other Includes  -----------

#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Epetra_MultiVector.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : BlockMultiVector::BlockMultiVector
// Purpose       : constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
BlockMultiVector::BlockMultiVector( int numBlocks, int numVectors,
                                    const Teuchos::RCP<const Parallel::ParMap> & globalMap,
                                    const Teuchos::RCP<const Parallel::ParMap> & subBlockMap
                                  )
: MultiVector( *globalMap, numVectors ),
  globalBlockSize_(subBlockMap->numGlobalEntities()),
  localBlockSize_(subBlockMap->numLocalEntities()),
  numBlocks_(numBlocks),
  startBlock_(0),
  endBlock_(numBlocks),
  newBlockMap_(subBlockMap),
  blocks_(numBlocks)
{

//  Using this Epetra constructor to view each block of multivectors:
//  Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map,
//         double **ArrayOfPointers, int NumVectors);

  //Setup Views of blocks using Block Map
  double ** Ptrs, ** Loc;
  Loc = (double**)malloc(sizeof(double*) * numVectors);

  epetraObj().ExtractView( &Ptrs );

  const Parallel::EpetraParMap& e_map = dynamic_cast<const Parallel::EpetraParMap&>(*newBlockMap_);

  for( int i = 0; i < numBlocks; ++i )
  {
    for( int j = 0; j < numVectors; ++j )
    {
      Loc[j] = Ptrs[j] + localBlockSize_*i;
    }
    blocks_[i] =  Teuchos::rcp( new MultiVector( new Epetra_MultiVector( View, *(e_map.petraMap()), Loc, numVectors ), true ) );
  }

  free(Loc);
}

//-----------------------------------------------------------------------------
// Function      : BlockMultiVector:::print
// Purpose       : Output
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
void BlockMultiVector::print(std::ostream &os) const
{
  os << "BlockMultiVector Object (Number of Blocks =" << numBlocks_ << ", Number of Vectors =" << numVectors() << std::endl;

  os << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  for( int i = 0; i < numBlocks_; ++i )
  {
    if (i >= startBlock_ && i < endBlock_)
    {
      os << "Block[" << i << "]\n";
    }
    blocks_[i]->print( os );
  }
  os << "Base Object\n";
  os << epetraObj();
  os << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

} // namespace Linear
} // namespace Xyce

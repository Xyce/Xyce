//----------------------------------------------------------------------
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
                                    const Teuchos::RCP<N_PDS_ParMap> & globalMap,
                                    const Teuchos::RCP<N_PDS_ParMap> & subBlockMap
                                  )
: MultiVector( *globalMap, numVectors ),
  blocksViewGlobalVec_(true),
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

  aMultiVector_->ExtractView( &Ptrs );

  for( int i = 0; i < numBlocks; ++i )
  {
    for( int j = 0; j < numVectors; ++j )
    {
      Loc[j] = Ptrs[j] + localBlockSize_*i;
    }
    blocks_[i] =  Teuchos::rcp( new MultiVector( new Epetra_MultiVector( View, dynamic_cast<const Epetra_BlockMap&>(*newBlockMap_->petraMap()), Loc, numVectors ), true ) );
  }

  free(Loc);
}

//-----------------------------------------------------------------------------
// Function      : BlockMultiVector::BlockMultiVector
// Purpose       : copy constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
BlockMultiVector::BlockMultiVector( const BlockMultiVector & rhs )
: MultiVector( dynamic_cast<const MultiVector&>(rhs) ),
  blocksViewGlobalVec_( rhs.blocksViewGlobalVec_ ),
  globalBlockSize_( rhs.globalBlockSize_ ),
  localBlockSize_( rhs.localBlockSize_ ),
  numBlocks_( rhs.numBlocks_ ),
  startBlock_( rhs.startBlock_ ),
  endBlock_( rhs.endBlock_ ),
  newBlockMap_( rhs.newBlockMap_ ),
  blocks_( rhs.blocks_.size() )
{
  if (blocksViewGlobalVec_)
  {
    // If the startBlock_ and endBlock_ cover every block in this vector than this is a time-domain representation
    // or serial simulation, in which case a frequency-domain distinction need not be made.
    if ((startBlock_ == 0) && (endBlock_ == numBlocks_))
    {
      int numBlocks = blocks_.size();

      // Setup Views of blocks using Block Map
      double ** Ptrs, ** Loc;
      Loc = (double**)malloc(sizeof(double*) * numVectors());

      aMultiVector_->ExtractView( &Ptrs );

      for( int i = 0; i < numBlocks; ++i )
      {
        for ( int j = 0; j < numVectors(); ++j )
        {
          Loc[j] = Ptrs[j] + localBlockSize_*i;
        }
        blocks_[i] =  Teuchos::rcp( new MultiVector( new Epetra_MultiVector( View, dynamic_cast<const Epetra_BlockMap&>(*(newBlockMap_->petraMap())), Loc, numVectors() ), true ) );
      }

      free(Loc);
    }
    else
    {
      // This is a frequency-domain representation of the block vector, so create views accordingly.
      int blockSize = globalBlockSize_;

      // Setup Views of blocks using Block Map
      double ** Ptrs, **Loc;
      Loc = (double**)malloc(sizeof(double*) * numVectors());

      aMultiVector_->ExtractView( &Ptrs );
      for ( int j = 0; j < numVectors(); ++j )
      {
        Loc[j] = Ptrs[j];
      } 

      for( int i = 0; i < numBlocks_; ++i )
      {
        // Create a MultiVector that views all the block data that is local.
        blocks_[i] =  Teuchos::rcp( new MultiVector( new Epetra_MultiVector( View, ((rhs.blocks_[i])->epetraObj()).Map(), Loc, numVectors() ), true ) );

        if ( (i >= startBlock_) && (i < endBlock_) )
        {
          for ( int j = 0; j < numVectors(); ++j )
          {
            // Advance the pointer for the local data.
            Loc[j] += blockSize;
          }
        }
      }

      free(Loc);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : BlockMultiVector::assembleGlobalVector
// Purpose       : Fills global MultiVector with the values in each block.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void BlockMultiVector::assembleGlobalVector()
{
  if (!blocksViewGlobalVec_)
  {

  }
}

//-----------------------------------------------------------------------------
// Function      : BlockMultiVector:::printPetraObject
// Purpose       : Output
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
void BlockMultiVector::printPetraObject(std::ostream &os) const
{
  os << "BlockMultiVector Object (Number of Blocks =" << numBlocks_ << ", Number of Vectors =" << numVectors() << ", View =" << blocksViewGlobalVec_ << std::endl;

  os << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  for( int i = 0; i < numBlocks_; ++i )
  {
    if (i >= startBlock_ && i < endBlock_)
    {
      os << "Block[" << i << "]\n";
    }
    blocks_[i]->printPetraObject( os );
  }
  os << "Base Object\n";
  os << *aMultiVector_;
  os << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

} // namespace Linear
} // namespace Xyce

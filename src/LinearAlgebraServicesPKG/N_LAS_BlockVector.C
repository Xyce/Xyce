//-------------------------------------------------------------------------
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
// Purpose        : Implementation file for Block Vector
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


#include <N_LAS_BlockVector.h>
#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>
#include <N_PDS_ParHelpers.h>
#include <N_PDS_EpetraParMap.h>

#include <N_LAS_BlockSystemHelpers.h>

// ---------  Other Includes  -----------

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : BlockVector::BlockVector
// Purpose       : constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
BlockVector::BlockVector( int numBlocks,
                          const Teuchos::RCP<const Parallel::ParMap> & globalMap,
                          const Teuchos::RCP<const Parallel::ParMap> & subBlockMap,
                          int augmentRows )
: Vector( *globalMap ),
  blocksViewGlobalVec_(true),
  globalBlockSize_(subBlockMap->numGlobalEntities()),
  localBlockSize_(subBlockMap->numLocalEntities()),
  overlapBlockSize_(subBlockMap->numLocalEntities()),
  numBlocks_(numBlocks),
  augmentCount_(augmentRows),
  startBlock_(0),
  endBlock_(numBlocks),
  newBlockMap_(subBlockMap),
  blocks_(numBlocks)
{
  //Setup Views of blocks using Block Map
  double ** Ptrs;
  epetraObj().ExtractView( &Ptrs );
  double * Loc;

  const Parallel::EpetraParMap& e_map = dynamic_cast<const Parallel::EpetraParMap&>(*newBlockMap_);

  for( int i = 0; i < numBlocks; ++i )
  {
    Loc = Ptrs[0] + overlapBlockSize_*i;
    blocks_[i] =  Teuchos::rcp( new Vector( new Epetra_Vector( View, *e_map.petraMap(), Loc ), true ) );
  }
}

//-----------------------------------------------------------------------------
// Function      : BlockVector::BlockVector
// Purpose       : constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
BlockVector::BlockVector( int blockSize,
                          const Teuchos::RCP<const Parallel::ParMap> & globalMap,
                          int augmentRows )
: Vector( *globalMap ),
  blocksViewGlobalVec_( true ),
  globalBlockSize_( blockSize ),
  localBlockSize_( blockSize ),
  overlapBlockSize_( blockSize ),
  numBlocks_( (globalMap->numGlobalEntities()-augmentRows) / blockSize ),
  augmentCount_( augmentRows ),
  startBlock_( 0 ),
  endBlock_( (globalMap->numGlobalEntities()-augmentRows) / blockSize ),
  blocks_( (globalMap->numGlobalEntities()-augmentRows) / blockSize )
{
  newBlockMap_ = Teuchos::rcp( Parallel::createPDSParMap( blockSize, blockSize, 
                               globalMap->indexBase(), globalMap->pdsComm() ) );

  // Determine where these blocks start and end in the grand scheme of things.
  startBlock_ = (int) std::floor( (double)(globalMap->minMyGlobalEntity() + 1) / (double)blockSize );
  endBlock_ = (int) std::floor( (double)(globalMap->maxMyGlobalEntity() + 1) / (double)blockSize );

  // Check for the augmented rows
  // Assume they are being placed on one processor.
  if (augmentRows && (globalMap->numLocalEntities() % blockSize))
  {
    endBlock_ = (int) std::floor( (double)(globalMap->maxMyGlobalEntity()-augmentRows + 1) / (double)blockSize );
  }

  //Setup Views of blocks using Block Map
  double ** Ptrs;
  epetraObj().ExtractView( &Ptrs );
  double * Loc = 0;
  if (globalMap->numLocalEntities() > 0)
  {
    Loc = Ptrs[0];
  }
  
  for( int i = 0; i < numBlocks_; ++i )
  {
    int myBlockSize = 0;
 
    // Generate maps where all the entries of the block are owned by one processor.
    if ( (i >= startBlock_) && (i < endBlock_) )
      myBlockSize = blockSize;

    Teuchos::RCP<Parallel::ParMap> currBlockMap = Teuchos::rcp( Parallel::createPDSParMap( blockSize, myBlockSize, 
                                                            globalMap->indexBase(), globalMap->pdsComm() ) );
    Teuchos::RCP<Parallel::EpetraParMap> e_currBlockMap = Teuchos::rcp_dynamic_cast<Parallel::EpetraParMap>(currBlockMap);

    // Create a Vector that views all the block data that is local.
    blocks_[i] =  Teuchos::rcp( new Vector( new Epetra_Vector( View, *e_currBlockMap->petraMap(), Loc ), true ) );

    if ( (i >= startBlock_) && (i < endBlock_) )
    {
      // Advance the pointer for the local data.
      Loc += blockSize;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : BlockVector::operator=
// Purpose       : assignment
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
BlockVector & BlockVector::operator=( const BlockVector & right )
{
  MultiVector::operator=( right ); 

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : BlockVector::BlockVector
// Purpose       : copy constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
BlockVector::BlockVector( const BlockVector & rhs )
: Vector( dynamic_cast<const Vector&>(rhs) ),
  blocksViewGlobalVec_( rhs.blocksViewGlobalVec_ ),
  globalBlockSize_( rhs.globalBlockSize_ ),
  localBlockSize_( rhs.localBlockSize_ ),
  overlapBlockSize_( rhs.overlapBlockSize_ ),
  numBlocks_( rhs.numBlocks_ ),
  augmentCount_( rhs.augmentCount_ ),
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
      double ** Ptrs;
      epetraObj().ExtractView( &Ptrs );
      double * Loc;

      const Parallel::EpetraParMap& e_map = dynamic_cast<const Parallel::EpetraParMap&>(*newBlockMap_);

      for( int i = 0; i < numBlocks; ++i )
      {
        Loc = Ptrs[0] + overlapBlockSize_*i;
        blocks_[i] =  Teuchos::rcp( new Vector( new Epetra_Vector( View, *e_map.petraMap(), Loc ), true ) );
      }
    }
    else
    {
      // This is a frequency-domain representation of the block vector, so create views accordingly.
      int blockSize = globalBlockSize_;

      // Setup Views of blocks using Block Map
      double ** Ptrs;
      epetraObj().ExtractView( &Ptrs );
      double * Loc = Ptrs[0];

      for( int i = 0; i < numBlocks_; ++i )
      {
        // Create a Vector that views all the block data that is local.
        blocks_[i] =  Teuchos::rcp( new Vector( new Epetra_Vector( View, ((rhs.blocks_[i])->epetraObj()).Map(), Loc ), true ) );

        if ( (i >= startBlock_) && (i < endBlock_) )
        {
          // Advance the pointer for the local data.
          Loc += blockSize;
        }
      }
    }
  }
  else
  {
    for( int i = 0; i < numBlocks_; ++i )
    {
      blocks_[i] =  Teuchos::rcp( rhs.blocks_[i]->cloneCopyVector() );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : BlockVector:::BlockVector
// Purpose       : view constructor
// Special Notes : Memory management is assumed to be outside this constructor
// Scope         : Public
// Creator       : Heidi Thornquist
// Creation Date : 11/12/20
//-----------------------------------------------------------------------------
BlockVector::BlockVector( const Vector * right, int blockSize )
: Vector( (const_cast<Vector*>(right)->epetraObj())(0), false ),
  blocksViewGlobalVec_( true ), 
  globalBlockSize_( blockSize ),
  localBlockSize_( blockSize ),
  overlapBlockSize_( blockSize ),
  numBlocks_( right->globalLength() / blockSize ),
  augmentCount_( right->globalLength() % blockSize ),
  startBlock_( 0 ),
  endBlock_( right->globalLength() / blockSize ),
  blocks_( right->globalLength() / blockSize )
{
  // If the oscillating HB algorithm is being used then augmentCount_ is probably not zero.
  int localAugmentCount = right->localLength() % blockSize;
  if (augmentCount_)
  {
    endBlock_ = (right->globalLength() - augmentCount_) / blockSize;
    blocks_.resize( endBlock_ );
    numBlocks_ = endBlock_; 
  }

  // Create the new maps for each block that places all the entries of the block on one processor.
  newBlockMap_ = Teuchos::rcp( Parallel::createPDSParMap( blockSize, blockSize, 
                               epetraObj().Map().IndexBase(), *right->pdsComm() ) );

  // Determine where these blocks start and end in the grand scheme of things.
  int minMyGID = (epetraObj().Map()).MinMyGID();
  int maxMyGID = (epetraObj().Map()).MaxMyGID();
  if ( localAugmentCount )
    maxMyGID -= localAugmentCount;

  startBlock_ = (int) std::floor( (double)(minMyGID + 1) / (double)blockSize );
  endBlock_ = (int) std::floor( (double)(maxMyGID + 1) / (double)blockSize );

  //Setup Views of blocks using Block Map
  double ** Ptrs;
  epetraObj().ExtractView( &Ptrs );
  double * Loc = Ptrs[0];

  for( int i = 0; i < numBlocks_; ++i )
  {
    int myBlockSize = 0;
 
    // Generate maps where all the entries of the block are owned by one processor.
    if ( (i >= startBlock_) && (i < endBlock_) )
      myBlockSize = blockSize;

    Teuchos::RCP<Parallel::ParMap> currBlockMap = Teuchos::rcp( Parallel::createPDSParMap( blockSize, myBlockSize,
                                                                epetraObj().Map().IndexBase(), *right->pdsComm() ) );
 
    Teuchos::RCP<Parallel::EpetraParMap> e_currBlockMap = Teuchos::rcp_dynamic_cast<Parallel::EpetraParMap>(currBlockMap);

    // Create a Vector that views all the block data that is local.
    blocks_[i] =  Teuchos::rcp( new Vector( new Epetra_Vector( View, *e_currBlockMap->petraMap(), Loc ), true ) );

    if ( (i >= startBlock_) && (i < endBlock_) )
    {
      // Advance the pointer for the local data.
      Loc += blockSize;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : BlockVector:::print
// Purpose       : Output
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
void BlockVector::print(std::ostream &os) const
{
  os << "BlockVector Object (Number of Blocks =" << numBlocks_ << ", View =" << blocksViewGlobalVec_ << ", Augmented Rows=" << augmentCount_ << ")" << std::endl;

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

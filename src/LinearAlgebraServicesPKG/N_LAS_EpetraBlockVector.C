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


#include <N_LAS_EpetraBlockVector.h>
#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>
#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>
#include <N_PDS_ParHelpers.h>
#include <N_PDS_EpetraParMap.h>

#include <N_LAS_BlockSystemHelpers.h>

#include <N_ERH_ErrorMgr.h>

// ---------  Other Includes  -----------

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <EpetraExt_MultiVectorOut.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : EpetraBlockVector::EpetraBlockVector
// Purpose       : constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
EpetraBlockVector::EpetraBlockVector( int numBlocks,
                          const Teuchos::RCP<const Parallel::ParMap> & globalMap,
                          const Teuchos::RCP<const Parallel::ParMap> & subBlockMap,
                          int augmentRows )
: BlockVector( *globalMap ),
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
// Function      : EpetraBlockVector::EpetraBlockVector
// Purpose       : constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
EpetraBlockVector::EpetraBlockVector( int blockSize,
                          const Teuchos::RCP<const Parallel::ParMap> & globalMap,
                          int augmentRows )
: BlockVector( *globalMap ),
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
// Function      : EpetraBlockVector::operator=
// Purpose       : assignment
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
EpetraBlockVector & EpetraBlockVector::operator=( const EpetraBlockVector & right )
{
  BlockVector::operator=( right ); 

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockVector::EpetraBlockVector
// Purpose       : view constructor
// Special Notes : Memory management is assumed to be outside this constructor
// Scope         : Public
// Creator       : Heidi Thornquist
// Creation Date : 11/12/20
//-----------------------------------------------------------------------------
EpetraBlockVector::EpetraBlockVector( const Vector * right, int blockSize )
: BlockVector( right ),
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
// Function      : cloneVector
// Purpose       : vector clone function
// Special Notes : clones shape, not values
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
Vector* EpetraBlockVector::cloneVector() const
{
  BlockVector* new_vec = new EpetraBlockVector( numBlocks_, Teuchos::rcp( pmap(), false ),
                                                newBlockMap_, augmentCount_ );
  return new_vec;
}
 
//-----------------------------------------------------------------------------
// Function      : cloneCopyVector
// Purpose       : vector clone function
// Special Notes : clones shape and values
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
Vector* EpetraBlockVector::cloneCopyVector() const
{
  BlockVector* new_vec = new EpetraBlockVector( numBlocks_, Teuchos::rcp( pmap(), false ),
                                                newBlockMap_, augmentCount_ );
  *new_vec = *this;  

  return new_vec;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockVector::dotProduct
// Purpose       : dot product
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
double EpetraBlockVector::dotProduct( const Vector & y ) const
{
  double result = 0.0;
  epetraObj().Dot(y.epetraObj(), &result);

  return result;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockVector::dotProduct
// Purpose       : dot product 
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
void EpetraBlockVector::dotProduct(const MultiVector & y, std::vector<double>& d) const
{
  int ynum = y.numVectors();
  for (int j=0; j<ynum; ++j)
  {
    epetraObj().Dot(*(y.epetraObj()(j)), &d[j]);
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockVector::lpNorm
// Purpose       : Returns lp norms of each vector in MultiVector
// Special Notes : Only p=1 and p=2 implemented now since this is all Petra
//                 supports.
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
int EpetraBlockVector::lpNorm(const int p, double * result) const
{ 
  if (p == 1)
    epetraObj().Norm1(result);
  else if (p == 2)
    epetraObj().Norm2(result);
  else
    Xyce::Report::DevelFatal0().in("EpetraBlockVector::lpNorm") << "Requested norm is not supported";

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockVector::wRMSNorm
// Purpose       : Returns weighted root-mean-square of each vector in
//                 MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/12/00
//-----------------------------------------------------------------------------
int EpetraBlockVector::wRMSNorm(const MultiVector & weights, double * result) const
{
  epetraObj().NormWeighted( weights.epetraObj(), result );
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockVector::infNorm
// Purpose       : Returns infinity norm of each blockVector
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/16/01
//-----------------------------------------------------------------------------
int EpetraBlockVector::infNorm(double * result, int * index) const
{
  if (!index)
  {
    epetraObj().NormInf(result);
  }
  else
  {
    int numProcs = pdsComm()->numProc();
    int numVecs = numVectors();
    int myLength = localLength();
    std::vector<int> indexTemp( numVecs, 0 ), indexTempAll( numVecs*numProcs, 0 );
    std::vector<double> doubleTemp( numVecs, 0.0 ), doubleTempAll( numVecs*numProcs, 0.0 );
    double ** pointers = epetraObj().Pointers();

    for (int i=0; i < numVecs; i++)
    {
      indexTemp[i] = -1;
      doubleTemp[i] = 0.0;
      for (int j=0; j < myLength; j++)
      {
        double tmp = std::abs(pointers[i][j]);
        if ( tmp > doubleTemp[i] )
        {
          doubleTemp[i] = tmp;
          indexTemp[i] = j;
        }
      }
      // Convert indexTemp from local to global ID
      if (indexTemp[i] > -1)
        indexTemp[i] = epetraObj().Map().GID(indexTemp[i]);
    }

    if (numProcs > 1)
    {
      // Use the communicator to gather all the local maximum values and indices
      Parallel::AllGather( pdsComm()->comm(), indexTemp, indexTempAll );
      Parallel::AllGather( pdsComm()->comm(), doubleTemp, doubleTempAll );

      // Compute the global infNorm and index
      for (int i=0; i < numVecs; i++)
      {
        result[i] = doubleTempAll[i];
        index[i] = indexTempAll[i];
        for (int j=1; j < numProcs; j++)
        {
          if ( doubleTempAll[j*numVecs + i] > result[i] )
          {
            result[i] = doubleTempAll[j*numVecs + i];
            index[i] = indexTempAll[j*numVecs + i];
          }
        }
      }
    }
    else
    {
      for (int i=0; i < numVecs; i++)
      {
        result[i] = doubleTemp[i];
        index[i] = indexTemp[i];
      }     
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockVector::wMaxNorm
// Purpose       : Returns the weighted inf-norm
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 03/19/01
//-----------------------------------------------------------------------------
int EpetraBlockVector::wMaxNorm(const MultiVector & weights, double * result) const
{
  int length  = localLength();
  int numVecs = numVectors();
  double tmpVal = 0.0;

  for (int i = 0;  i < numVecs; ++i)
  {
    double localMax = 0.0;
    if (length)
    {
      localMax = fabs(*(*this)(0,i)) / (*weights(0,i));
      for (int j = 1; j < length; ++j)
      {
        tmpVal = fabs(*(*this)(j,i)) / (*weights(j,i));
        if (tmpVal > localMax)
          localMax = tmpVal;
      }
    }
    // Determine global maximum.
    pdsComm()->maxAll( &localMax, &(result[i]), 1 );
  }

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockVector::getElementByGlobalIndex
// Purpose       : Get element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/6/00
//-----------------------------------------------------------------------------
const double & EpetraBlockVector::getElementByGlobalIndex(
  const int & global_index, const int & vec_index) const
{
  if( pmap() == NULL )
    return epetraObj()[vec_index][ epetraObj().Map().LID(global_index) ];
  else
  {
    int i = pmap()->globalToLocalIndex(global_index);

    if (i != -1)
      return epetraObj()[vec_index][i];
    else 
    {
      Xyce::Report::DevelFatal().in("getElementByGlobalIndex") 
        << "Failed to find BlockVector global index: " << global_index;
      return groundNode_;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockVector::setElementByGlobalIndex
// Purpose       : Set element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/6/00
//-----------------------------------------------------------------------------
bool EpetraBlockVector::setElementByGlobalIndex(const int & global_index,
                                                const double & val,
                                                const int & vec_index)
{
  if( pmap() == NULL )
    epetraObj()[vec_index][ epetraObj().Map().LID(global_index) ] = val;
  else
  {
    if (global_index != -1)
    {
      int i = pmap()->globalToLocalIndex(global_index);

      if (i != -1)
        epetraObj()[vec_index][i] = val;
      else
      {
        Xyce::Report::DevelFatal().in("setElementByGlobalIndex") 
          << "Failed to find BlockVector global index: " << global_index;
        return false;
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockVector::sumElementByGlobalIndex
// Purpose       : Set element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/7/00
//-----------------------------------------------------------------------------
bool EpetraBlockVector::sumElementByGlobalIndex(const int & global_index,
                                                const double & val,
                                                const int & vec_index)
{ 
  if( pmap() == NULL )
    epetraObj()[vec_index][ epetraObj().Map().LID(global_index) ] += val;
  else
  { 
    if (global_index != -1 )
    { 
      int i = pmap()->globalToLocalIndex(global_index);
      
      if (i != -1)
        epetraObj()[vec_index][i] += val;
      else
      { 
        Report::DevelFatal() 
          << " sumElementByGlobalIndex: failed to find BlockVector global index ";
        return false;
      }
    }
  }
  
  return true;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockVector::writeToFile
// Purpose       : Dumps out the multivector entries to a file.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/19/00
//-----------------------------------------------------------------------------
void EpetraBlockVector::writeToFile( const char * filename, bool useLIDs, bool mmFormat ) const
{
  int numProcs = pdsComm()->numProc();
  int localRank = pdsComm()->procID();
  int masterRank = 0;

  if (!mmFormat)
  {
    for( int p = 0; p < numProcs; ++p )
    {
      //A barrier inside the loop so each processor waits its turn.
      pdsComm()->barrier();

      if(p == localRank)
      {
        FILE *file = NULL;

        if(masterRank == localRank)
        {
          //This is the master processor, open a new file.
          file = fopen(filename,"w");

          //Write the dimension n into the file.
          fprintf(file,"%d\n",globalLength());
        }
        else
        {
          //This is not the master proc, open file for appending
          file = fopen(filename,"a");
        }

        //Now loop over the local portion of the block vector.
        int length  = localLength();
        int numVecs = numVectors();

        for (int i = 0; i < numVecs; ++i)
          for (int j = 0; j < length; ++j)
          {
            int loc = epetraObj().Map().GID(j);
            if( useLIDs ) loc = j;
            fprintf(file,"%d %d %20.13e\n",i,loc,epetraObj()[i][j]);
          } 
        fclose(file);
      } 
    } 
  } 
  else
  {
    EpetraExt::MultiVectorToMatrixMarketFile( filename, epetraObj() );
  } 
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockVector:::print
// Purpose       : Output
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
void EpetraBlockVector::print(std::ostream &os) const
{
  os << "EpetraBlockVector Object (Number of Blocks =" << numBlocks_ << ", Augmented Rows=" << augmentCount_ << ")" << std::endl;

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

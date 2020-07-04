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
// Purpose        : Implementation  for Block Matrix
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/13/04
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_fwd.h>

#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockVector.h>

// ----------   Other Includes -----------

#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Teuchos_TestForException.hpp>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : BlockMatrix::BlockMatrix
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
BlockMatrix::BlockMatrix( int size,
                          int offset,
                          const std::vector< std::vector<int> > & blockColumns,
                          const Graph* globalGraph,
                          const Graph* subBlockGraph,
                          int augmentCount )

: Matrix( new Epetra_CrsMatrix( Copy, *(globalGraph->epetraObj()) ), true ),
  blocksViewGlobalMat_(true),
  blockSize_( subBlockGraph->epetraObj()->NumMyRows() ),
  offset_( offset ),
  numBlockRows_( size ),
  augmentCount_( augmentCount ),
  cols_( blockColumns ),
  blocks_( size )
{
  // Individual blocks cannot be a view of the global matrix because of graph ordering and communication
  if ( globalGraph->epetraObj()->Comm().NumProc() > 1 )
  {
    blocksViewGlobalMat_ = false;

    for( int i = 0; i < numBlockRows_; ++i )
    {
      int numCols = cols_[i].size();
      blocks_[i].resize( numCols );
      for( int j = 0; j < numCols; ++j )
      {
        Epetra_CrsMatrix * bMat = new Epetra_CrsMatrix( Copy, *(subBlockGraph->epetraObj()) );
        blocks_[i][j] = Teuchos::rcp( new Matrix( bMat ) );
      }
    }

    // Get the local indices for the sub block so assembling is easier
    baseNumCols_.resize( blockSize_ );
    baseIndices_.resize( subBlockGraph->epetraObj()->NumMyNonzeros() );
    int ptr = 0;
    for( int i = 0; i < blockSize_; ++i )
    {
      subBlockGraph->epetraObj()->ExtractMyRowCopy( i, subBlockGraph->epetraObj()->NumMyNonzeros()-ptr, baseNumCols_[i], &baseIndices_[ptr] );
      ptr += baseNumCols_[i];
    }
  }
  else
  {
    std::vector<int> baseNumCols( blockSize_ );
    std::vector< int* > baseIndices( blockSize_ );
    for( int i = 0; i < blockSize_; ++i )
      subBlockGraph->epetraObj()->ExtractMyRowView( i, baseNumCols[i], baseIndices[i] );

    std::vector< double* > Values( blockSize_ );
    int NumEntries;

    for( int i = 0; i < numBlockRows_; ++i )
    {
      int numCols = cols_[i].size();
      blocks_[i].resize( numCols );

      for( int j = 0; j < blockSize_; ++j )
        aDCRSMatrix_->ExtractMyRowView( j+blockSize_*i, NumEntries, Values[j] );

      for( int j = 0; j < numCols; ++j )
      {
        Epetra_CrsMatrix * bMat = new Epetra_CrsMatrix( View, *(subBlockGraph->epetraObj()) );

        for( int k = 0; k < blockSize_; ++k )
          bMat->InsertMyValues( k, baseNumCols[k], Values[k]+j*baseNumCols[k], baseIndices[k] );

        blocks_[i][j] = Teuchos::rcp( new Matrix( bMat ) );
      }
    }
  }

  // Generate the augmented GIDs list.
  if( augmentCount_ )
  {
    augmentGIDs_.resize(augmentCount_);
    int augStart = blockSize_ * size;
    for( int i = 0; i < augmentCount_; ++i )
    {
      augmentGIDs_[i] = globalGraph->epetraObj()->RowMap().GID(augStart+i);
    }
  }

  // Communicate the augmented GIDs to all processors.
  // All other processors other than the one that owns the augmented GID will have -1.
  std::vector<int> tmpAugmentGIDs = augmentGIDs_;
  globalGraph->epetraObj()->Comm().MaxAll( &tmpAugmentGIDs[0], &augmentGIDs_[0], augmentCount_ );
}

//-----------------------------------------------------------------------------
// Function      : BlockMatrix::block
// Purpose       : Block Accessor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
Matrix & BlockMatrix::block( int row, int col )
{
  for( int i = 0; i < (int)(cols_[row].size()); ++i )
    if( cols_[row][i] == col )
      return *blocks_[row][i];

  TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
      "Error!  BlockMatrix::block("<<row<<","<<col<<"):  This block does not exist!"
      );
}

//-----------------------------------------------------------------------------
// Function      : BlockMatrix::put
// Purpose       : Put function for the sparse-matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void BlockMatrix::put( double s )
{
  aDCRSMatrix_->PutScalar(s);

  if (!blocksViewGlobalMat_)
  {
    for( int i = 0; i < numBlockRows_; ++i )
    {
      int numCols = cols_[i].size();
      for ( int j = 0; j < numCols; ++j )
      {
        blocks_[i][j]->put( s );
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : BlockMatrix::fillComplete
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/29/03
//-----------------------------------------------------------------------------
void BlockMatrix::fillComplete()
{
  // Call fillComplete on all the individual Matrix blocks, 
  // then assemble the global matrix.
  for( int i = 0; i < numBlockRows_; ++i )
  {
    int numCols = cols_[i].size();
    for ( int j = 0; j < numCols; ++j )
    {
      blocks_[i][j]->fillComplete();
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : BlockMatrix::assembleGlobalMatrix
// Purpose       : Fills global Matrix with the values in each block.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 10/03/13
//-----------------------------------------------------------------------------
void BlockMatrix::assembleGlobalMatrix()
{
  if (!blocksViewGlobalMat_)
  {
    for( int i = 0; i < numBlockRows_; ++i )
    {
      int numCols = cols_[i].size();

      for( int k = 0; k < blockSize_; ++k )
      {
        // Create memory for all the entries for one whole row.
        int length = numCols*baseNumCols_[k];
        std::vector<int> Indices( length );
        std::vector<double> Values( length );

        // For each block column extract the current values from the subblock and correct the column indices.
        // NOTE:  All extractions and insertions are done using global ids, it seems easier that way.
        int ptr = 0;
        for( int j = 0; j < numCols; ++j )
        {
          int numIndices = 0;
          int col = cols_[i][j];
          int globalRow = (blocks_[i][j])->epetraObj().Graph().GRID(k);
          (*blocks_[i][j]).getRowCopy(globalRow, length-ptr, numIndices, &Values[ptr], &Indices[ptr]);
        
          // Correct the column indices for this subblock. 
          for (int idx = 0; idx < numIndices; idx++)
          {
            Indices[ptr+idx] += offset_*col;
          }
          ptr += numIndices;
        }
       
        // Insert the values for all the global columns at the same time.
        aDCRSMatrix_->ReplaceGlobalValues(aDCRSMatrix_->Graph().GRID(k + i*blockSize_), length, &Values[0], &Indices[0]); 
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : BlockMatrix::printPetraObject
// Purpose       : Output
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
void BlockMatrix::printPetraObject(std::ostream &os) const
{
  os << "BlockMatrix Object (Size=" << numBlockRows_ << ", View =" << blocksViewGlobalMat_ << ")" << std::endl;
  os << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  for( int i = 0; i < numBlockRows_; ++i )
  {
    int numCols = cols_[i].size();
    for( int j = 0; j < numCols; ++j )
    {
      os << "Block[" << i << "][" << cols_[i][j] << "]\n";
      blocks_[i][j]->printPetraObject(os);
    }
  }
  os << "Base Object\n";
  os << *aDCRSMatrix_;
  os << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

//-----------------------------------------------------------------------------
// Function      : BlockMatrix::replaceAugmentedColumn
// Purpose       : Replace augmented column of matrix with block vector
// Special Notes : This places values directly in the base object.
//               : This should NOT be used to replace values that are within blocks.
// Scope         : Public
// Creator       : Todd Coffey, SNL, 1414
// Creation Date : 08/12/05
//-----------------------------------------------------------------------------
void BlockMatrix::replaceAugmentedColumn(int colGID, const BlockVector & vec)
{
  const Epetra_BlockMap & RowMap = const_cast<BlockVector&>(vec).epetraObj().Map();
  int NumRows = RowMap.NumMyElements();

  int lCol = aDCRSMatrix_->Graph().LCID(colGID);
  for( int i = 0; i < NumRows; ++i )
  {
    int Row = i;
    double Val = vec[i];
    aDCRSMatrix_->ReplaceMyValues( Row, 1, &Val, &lCol );
  }
}

//-----------------------------------------------------------------------------
// Function      : BlockMatrix::replaceAugmentedRow
// Purpose       : Replace augmented row values
// Special Notes : This places values directly in the base object.
//               : This should NOT be used to replace values that are within blocks.
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 11/12/13
//-----------------------------------------------------------------------------
void BlockMatrix::replaceAugmentedRow(int rowGID, int length, double * coeffs, int * colIndices)
{
  if ( aDCRSMatrix_->Graph().LRID( rowGID ) >= 0 )
  {
    aDCRSMatrix_->ReplaceGlobalValues( rowGID, length, coeffs, colIndices );
  }
}

} // namespace Linear
} // namespace Xyce

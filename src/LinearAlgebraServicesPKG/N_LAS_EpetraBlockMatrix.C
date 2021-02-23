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
// Purpose        : Epetra implementation for BlockMatrix
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
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_FeatureTest.h>

#include <N_LAS_EpetraBlockMatrix.h>
#include <N_LAS_EpetraMatrix.h>
#include <N_LAS_EpetraGraph.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>
#include <N_PDS_EpetraParMap.h>

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <EpetraExt_RowMatrixOut.h>
#include <Epetra_OffsetIndex.h>
// ----------   Other Includes -----------

#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_MultiVector.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::EpetraBlockMatrix
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
EpetraBlockMatrix::EpetraBlockMatrix( int size,
                                      int offset,
                                      const std::vector< std::vector<int> > & blockColumns,
                                      const Graph* globalGraph,
                                      const Graph* subBlockGraph,
                                      int augmentCount )

: aDCRSMatrix_(0),
  aColMap_(0),
  baseGraph_(0),
  isOwned_(true),
  blocksViewGlobalMat_(true),
  blockSize_( subBlockGraph->numLocalEntities() ),
  offset_( offset ),
  numBlockRows_( size ),
  augmentCount_( augmentCount ),
  cols_( blockColumns ),
  blocks_( size )
{
  aDCRSMatrix_ = new Epetra_CrsMatrix( Copy, *(dynamic_cast<const EpetraGraph *>(globalGraph)->epetraObj()) );

  baseGraph_ = new EpetraGraph( Teuchos::rcp( const_cast<Epetra_CrsGraph*>(&(aDCRSMatrix_->Graph())), false ) );

  // Individual blocks cannot be a view of the global matrix because of graph ordering and communication
  if ( aDCRSMatrix_->Comm().NumProc() > 1 )
  {
    blocksViewGlobalMat_ = false;

    for( int i = 0; i < numBlockRows_; ++i )
    {
      int numCols = cols_[i].size();
      blocks_[i].resize( numCols );
      for( int j = 0; j < numCols; ++j )
      {
        Epetra_CrsMatrix * bMat = 
          new Epetra_CrsMatrix( Copy, *(dynamic_cast<const EpetraGraph *>(subBlockGraph)->epetraObj()) );
        blocks_[i][j] = Teuchos::rcp( new EpetraMatrix( bMat ) );
      }
    }

    // Get the local indices for the sub block so assembling is easier
    baseNumCols_.resize( blockSize_ );
    baseIndices_.resize( subBlockGraph->numLocalNonzeros() );
    int ptr = 0;
    for( int i = 0; i < blockSize_; ++i )
    {
      subBlockGraph->extractLocalRowCopy( i, subBlockGraph->numLocalNonzeros()-ptr, baseNumCols_[i], &baseIndices_[ptr] );
      ptr += baseNumCols_[i];
    }
  }
  else
  {
    std::vector<int> baseNumCols( blockSize_ );
    std::vector< int* > baseIndices( blockSize_ );
    for( int i = 0; i < blockSize_; ++i )
      subBlockGraph->extractLocalRowView( i, baseNumCols[i], baseIndices[i] );

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
        Epetra_CrsMatrix * bMat = 
          new Epetra_CrsMatrix( View, *(dynamic_cast<const EpetraGraph *>(subBlockGraph)->epetraObj()) );

        for( int k = 0; k < blockSize_; ++k )
          bMat->InsertMyValues( k, baseNumCols[k], Values[k]+j*baseNumCols[k], baseIndices[k] );

        blocks_[i][j] = Teuchos::rcp( new EpetraMatrix( bMat ) );
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
      augmentGIDs_[i] = globalGraph->localToGlobalRowIndex( augStart+i );
    }
  }

  // Communicate the augmented GIDs to all processors.
  // All other processors other than the one that owns the augmented GID will have -1.
  std::vector<int> tmpAugmentGIDs = augmentGIDs_;
  aDCRSMatrix_->Comm().MaxAll( &tmpAugmentGIDs[0], &augmentGIDs_[0], augmentCount_ );
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::~EpetraBlockMatrix
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
EpetraBlockMatrix::~EpetraBlockMatrix()
{ 
  if ( isOwned_ )
  { 
    delete aDCRSMatrix_;
  }

  delete aColMap_;
  delete baseGraph_;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::block
// Purpose       : Block Accessor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
Matrix & EpetraBlockMatrix::block( int row, int col )
{
  for( int i = 0; i < (int)(cols_[row].size()); ++i )
  {
    if( cols_[row][i] == col )
      return *blocks_[row][i];
  }

  Report::UserFatal0() << " EpetraBlockMatrix::block( "<<row<<", "<<col<<" ):  This block does not exist!";
  return *blocks_[row][0];  // We shouldn't get here.
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::put
// Purpose       : Put function for the sparse-matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void EpetraBlockMatrix::put( double s )
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
// Function      : EpetraBlockMatrix::assembleGlobalMatrix
// Purpose       : Fills global Matrix with the values in each block.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 10/03/13
//-----------------------------------------------------------------------------
void EpetraBlockMatrix::assembleGlobalMatrix()
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
// Function      : EpetraBlockMatrix::print
// Purpose       : Output
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
void EpetraBlockMatrix::print(std::ostream &os) const
{
  os << "BlockMatrix Object (Size=" << numBlockRows_ << ", View =" << blocksViewGlobalMat_ << ")" << std::endl;
  os << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  for( int i = 0; i < numBlockRows_; ++i )
  {
    int numCols = cols_[i].size();
    for( int j = 0; j < numCols; ++j )
    {
      os << "Block[" << i << "][" << cols_[i][j] << "]\n";
      blocks_[i][j]->print(os);
    }
  }
  os << "Base Object\n";
  os << *aDCRSMatrix_;
  os << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::replaceAugmentedColumn
// Purpose       : Replace augmented column of matrix with block vector
// Special Notes : This places values directly in the base object.
//               : This should NOT be used to replace values that are within blocks.
// Scope         : Public
// Creator       : Todd Coffey, SNL, 1414
// Creation Date : 08/12/05
//-----------------------------------------------------------------------------
void EpetraBlockMatrix::replaceAugmentedColumn(int colGID, const BlockVector & vec)
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
// Function      : EpetraBlockMatrix::replaceAugmentedRow
// Purpose       : Replace augmented row values
// Special Notes : This places values directly in the base object.
//               : This should NOT be used to replace values that are within blocks.
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 11/12/13
//-----------------------------------------------------------------------------
void EpetraBlockMatrix::replaceAugmentedRow(int rowGID, int length, double * coeffs, int * colIndices)
{
  if ( aDCRSMatrix_->Graph().LRID( rowGID ) >= 0 )
  {
    aDCRSMatrix_->ReplaceGlobalValues( rowGID, length, coeffs, colIndices );
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::add
// Purpose       : Sums in a matrix contribution
// Special Notes : WARNING: only works if graphs match, no checking
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/11/04
//-----------------------------------------------------------------------------
void EpetraBlockMatrix::add( const Matrix & A )
{
  int NumRows = A.getLocalNumRows();
  int* Indices;
  double* Values;
  int NumIndices;

  for( int i = 0; i < NumRows; ++i )
  {
    A.getLocalRowView( i, NumIndices, Values, Indices );
    aDCRSMatrix_->SumIntoMyValues( i, NumIndices, Values, Indices );
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::linearCombo
// Purpose       : Sums in a matrix contribution
// Special Notes : WARNING: only works if graphs EXACTLY match no checking
//
//                 this = a*A + b*B
//
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 2/13/07
//-----------------------------------------------------------------------------
void EpetraBlockMatrix::linearCombo( const double a, const Matrix & A,
                                     const double b, const Matrix & B)
{
  int NumRows = aDCRSMatrix_->NumMyRows();

  int *aIndices, *bIndices;
  int aNumIndices, bNumIndices;
  double *aValues, *bValues;

  for( int i = 0; i < NumRows; ++i )
  {
    // Get a view of the i-th row for A and B.
    A.getLocalRowView( i, aNumIndices, aValues, aIndices );
    B.getLocalRowView( i, bNumIndices, bValues, bIndices );

    // Add in the entries from each matrix.
    for ( int j = 0; j < aNumIndices; ++j )
      (*aDCRSMatrix_)[i][j] = a*aValues[j] + b*bValues[j];
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::matvec
// Purpose       : Sparse-matrix vector multiply - multivector version.  This
//                 function forms the product y = Ax where x and y are
//                 multivectors.  If transA is true, multiply by the transpose
//                 of matrix, otherwise just use matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void EpetraBlockMatrix::matvec(bool transA, const MultiVector &x,
                          MultiVector &y)
{
  int PetraError = aDCRSMatrix_->Multiply(transA, x.epetraObj(), y.epetraObj());

  if (DEBUG_LINEAR)
    processError( "EpetraBlockMatrix::matvec - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::getDiagonal
// Purpose       : Return the diagonal entries of the sparse matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void EpetraBlockMatrix::getDiagonal( Vector & diagonal ) const
{
  int PetraError = aDCRSMatrix_->ExtractDiagonalCopy( diagonal.epetraObj() );

  if (DEBUG_LINEAR)
    processError( "EpetraBlockMatrix::getDiagonal - ", PetraError );
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::replaceDiagonal
// Purpose       : Replace values of diagonal elements
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/05/03
//-----------------------------------------------------------------------------
bool EpetraBlockMatrix::replaceDiagonal( const Vector & vec )
{
  const Epetra_Vector * eVec = vec.epetraObj()(0);
  int PetraError = aDCRSMatrix_->ReplaceDiagonalValues( *eVec );

  if (DEBUG_LINEAR)
    processError( "EpetraBlockMatrix::replaceDiagonal - ", PetraError );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::getLocalRowView
// Purpose       : Returns row coefficients and associated column indices.
// Special Notes : Uses Petra's ExtractRowView which does not require user
//               : to setup space.
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
int EpetraBlockMatrix::getLocalRowView(int lidRow, int& numEntries, double*& values, int*& indices) const
{
  return aDCRSMatrix_->ExtractMyRowView(lidRow, numEntries, values, indices);
}


//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::getRowCopy
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/05/06
//-----------------------------------------------------------------------------
void EpetraBlockMatrix::getRowCopy
  (int row, int length, int & numEntries, double *coeffs, int *colIndices) const
{
  int PetraError = aDCRSMatrix_->ExtractGlobalRowCopy
       (row, length, numEntries, coeffs, colIndices);

  if (DEBUG_LINEAR)
    processError( "EpetraBlockMatrix::getRowCopy - ", PetraError );
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::getLocalRowCopy
// Purpose       :
// Special Notes :
//               :
//               :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/05/06
//-----------------------------------------------------------------------------
void EpetraBlockMatrix::getLocalRowCopy
  (int row, int length, int & numEntries, double *coeffs, int *colIndices) const
{
  int PetraError = aDCRSMatrix_->ExtractMyRowCopy
       (row, length, numEntries, coeffs, colIndices);

  if (DEBUG_LINEAR)
    processError( "EpetraBlockMatrix::getLocalRowCopy - ", PetraError );
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::getColMap
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 9/6/17
//-----------------------------------------------------------------------------
const Parallel::ParMap* EpetraBlockMatrix::getColMap( const Parallel::Communicator& comm ) const
{
  if (!aColMap_)
    aColMap_ = new Parallel::EpetraParMap( &aDCRSMatrix_->ColMap(), comm );
 
  return aColMap_;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::sumIntoLocalRow
// Purpose       : Sum values into a row into the sparse matrix, using local indices.
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/30/10
//-----------------------------------------------------------------------------
bool EpetraBlockMatrix::addIntoLocalRow(int row, int length, const double * coeffs,
                                        const int * colIndices)
{
  int PetraError = aDCRSMatrix_->SumIntoMyValues(row, length, coeffs, colIndices);

  if (DEBUG_LINEAR | DEBUG_DEVICE)
    processError( "EpetraBlockMatrix::addIntoLocalRow - ", PetraError );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::putLocalRow
// Purpose       : Put values into a row into the sparse matrix, using local indices.
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/29/12
//-----------------------------------------------------------------------------
bool EpetraBlockMatrix::putLocalRow(int row, int length, const double * coeffs,
                                    const int * colIndices)
{
  int PetraError = aDCRSMatrix_->ReplaceMyValues(row, length, coeffs, colIndices);

  if (DEBUG_LINEAR | DEBUG_DEVICE)
    processError( "EpetraBlockMatrix::putLocalRow - ", PetraError );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::putRow
// Purpose       : Put a row into the sparse matrix.
// Special Notes : Replace already allocated values
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
bool EpetraBlockMatrix::putRow(int row, int length, const double *coeffs, const int *colIndices)
{
  int PetraError = aDCRSMatrix_->ReplaceGlobalValues(row, length, coeffs, colIndices);

  if (DEBUG_LINEAR)
    processError( "EpetraBlockMatrix::putRow - ", PetraError );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMatrix::processError
// Purpose       : Concrete implementation which processes Petra (in this case)
//                 error codes taken from the Petra member function returns.
// Special Notes : Petra specific.  NOTE ALSO - this function is currently
//                 within the "Xyce_DEBUG_LINEAR" ifdef and so any calls to
//                 this should also be so bracketed.
// Scope         : Private
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void EpetraBlockMatrix::processError(std::string methodMsg, int error) const
{
  const std::string PetraError("Function returned with an error.\n");

  // Process the error
  if( error < 0 )
    Report::DevelFatal0() << methodMsg + PetraError;
}

} // namespace Linear
} // namespace Xyce

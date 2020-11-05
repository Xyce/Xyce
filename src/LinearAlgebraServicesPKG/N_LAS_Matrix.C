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
// Purpose        : Implemenation file for the Abstract interface to sparse
//                  matrix type.
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/20/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <algorithm>

// ----------   Xyce Includes   ----------

#include <N_UTL_fwd.h>

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Graph.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>
#include <N_PDS_EpetraParMap.h>
#include <N_UTL_FeatureTest.h>

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <EpetraExt_RowMatrixOut.h>
#include <Epetra_OffsetIndex.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : Matrix::~Matrix
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
Matrix::~Matrix()
{
  if ( isOwned_ ) 
  {
    if( oDCRSMatrix_ != aDCRSMatrix_ )
    {
      delete aDCRSMatrix_;
    }

    if( oDCRSMatrix_ ) 
      delete oDCRSMatrix_;
  }

  delete exporter_;
  delete offsetIndex_;
  delete aColMap_;
  delete oColMap_;

  if (overlapGraph_ != baseGraph_)
  {
    delete baseGraph_;
  }
  if (overlapGraph_)
    delete overlapGraph_;
}


//-----------------------------------------------------------------------------
// Function      : Matrix::Matrix
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
Matrix::Matrix( N_PDS_ParMap & map, std::vector<int> & diagArray )
: aDCRSMatrix_(0),
  oDCRSMatrix_(0),
  exporter_(0),
  offsetIndex_(0),
  aColMap_(0),
  oColMap_(0),
  overlapGraph_(0),
  baseGraph_(0),
  proxy_( 0, *this ),
  groundLID_(-1),
  groundNode_(0.0),
  isOwned_(true)
{
  N_PDS_EpetraParMap& e_map = dynamic_cast<N_PDS_EpetraParMap&>( map );
  aDCRSMatrix_ = new Epetra_CrsMatrix( Copy, *e_map.petraMap() , &(diagArray[0]) );
  oDCRSMatrix_ = aDCRSMatrix_;

  baseGraph_ = new Graph( Teuchos::rcp( &(aDCRSMatrix_->Graph()), false ) );
  overlapGraph_ = baseGraph_;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::Matrix
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
Matrix::Matrix( Epetra_CrsMatrix * origMatrix, bool isOwned )
: aDCRSMatrix_( origMatrix ),
  exporter_(0),
  offsetIndex_(0),
  aColMap_(0),
  oColMap_(0),
  overlapGraph_(0),
  baseGraph_(0),
  proxy_( 0, *this ),
  groundLID_(-1),
  groundNode_(0.0),
  isOwned_(isOwned)
{
  oDCRSMatrix_ = aDCRSMatrix_;

  baseGraph_ = new Graph( Teuchos::rcp( &(aDCRSMatrix_->Graph()), false ) );
  overlapGraph_ = baseGraph_;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::Matrix
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/21/02
//-----------------------------------------------------------------------------
Matrix::Matrix( const Graph* overlapGraph,
                const Graph* baseGraph )
: aDCRSMatrix_(0),
  oDCRSMatrix_(0),
  exporter_(0),
  offsetIndex_(0),
  aColMap_(0),
  oColMap_(0),
  overlapGraph_(0),
  baseGraph_(0),
  proxy_( 0, *this ),
  groundLID_(-1),
  groundNode_(0.0),
  isOwned_(true)
{
  if ( baseGraph!= overlapGraph )
  {
    oDCRSMatrix_ = new Epetra_CrsMatrix( Copy, *(overlapGraph->epetraObj()) );

    // Get ground node, if there is one.
    groundLID_ = overlapGraph->epetraObj()->LRID( -1 );

    aDCRSMatrix_ = new Epetra_CrsMatrix( Copy, *(baseGraph->epetraObj()) );
    exporter_ = new Epetra_Export( overlapGraph->epetraObj()->RowMap(), baseGraph->epetraObj()->RowMap() );
    offsetIndex_ = new Epetra_OffsetIndex( *(overlapGraph->epetraObj()), *(baseGraph->epetraObj()), *exporter_ );
  }
  else
  {
    aDCRSMatrix_ = new Epetra_CrsMatrix( Copy, *(baseGraph->epetraObj()) );
    oDCRSMatrix_ = aDCRSMatrix_;
  }

  overlapGraph_ = new Graph( *overlapGraph );
  baseGraph_ = new Graph( *baseGraph );
}

//-----------------------------------------------------------------------------
// Function      : Matrix::fillComplete
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/29/03
//-----------------------------------------------------------------------------
void Matrix::fillComplete()
{
  if( exporter_ )
  {
    aDCRSMatrix_->Export( *oDCRSMatrix_, *exporter_, Add, offsetIndex_ );
  }
}

//-----------------------------------------------------------------------------
// Function      : Matrix::matvec
// Purpose       : Sparse-matrix vector multiply - multivector version.  This
//                 function forms the product y = Ax where x and y are
//                 multivectors.  If transA is true, multiply by the transpose
//                 of matrix, otherwise just use matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void Matrix::matvec(bool transA, const MultiVector &x,
                          MultiVector &y)
{
  int PetraError = aDCRSMatrix_->Multiply(transA, *(x.aMultiVector_),
					  *(y.aMultiVector_));

  if (DEBUG_LINEAR)
    processError( "Matrix::matvec - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : Matrix::put
// Purpose       : Put function for the sparse-matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void Matrix::put( double s )
{
  if ( exporter_ )
  {
    aDCRSMatrix_->PutScalar(s);
  }
  oDCRSMatrix_->PutScalar(s);
  groundNode_ = s;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::scale
// Purpose       : Scale function for the sparse-matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void Matrix::scale(double scaleFactor)
{
  if ( exporter_ )
  {
    aDCRSMatrix_->Scale(scaleFactor);
  }
  oDCRSMatrix_->Scale(scaleFactor);
}

//-----------------------------------------------------------------------------
// Function      : Matrix::getRowLength
// Purpose       : Returns the number of nonzeroes in the row.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
int Matrix::getRowLength(int row) const
{
  return aDCRSMatrix_->NumGlobalEntries(row);
}

//-----------------------------------------------------------------------------
// Function      : Matrix::getLocalRowView
// Purpose       : Returns row coefficients and associated column indices.
// Special Notes : Uses Petra's ExtractRowView which does not require user
//               : to setup space.
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void Matrix::getLocalRowView(int row, int &length, double *coeffs, int *colIndices) const
{
  int PetraError = aDCRSMatrix_->ExtractMyRowView(row, length, coeffs, colIndices);

  if (DEBUG_LINEAR)
    processError( "Matrix::getLocalRowView - ", PetraError );
}


//-----------------------------------------------------------------------------
// Function      : Matrix::getRowCopy
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/05/06
//-----------------------------------------------------------------------------
void Matrix::getRowCopy
  (int row, int length, int & numEntries, double *coeffs, int *colIndices) const
{
  int PetraError = aDCRSMatrix_->ExtractGlobalRowCopy
       (row, length, numEntries, coeffs, colIndices);

  if (DEBUG_LINEAR)
    processError( "Matrix::getRowCopy - ", PetraError );
}

//-----------------------------------------------------------------------------
// Function      : Matrix::getLocalRowCopy
// Purpose       :
// Special Notes :
//               :
//               :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/05/06
//-----------------------------------------------------------------------------
void Matrix::getLocalRowCopy
  (int row, int length, int & numEntries, double *coeffs, int *colIndices) const
{
  int PetraError = aDCRSMatrix_->ExtractMyRowCopy
       (row, length, numEntries, coeffs, colIndices);

  if (DEBUG_LINEAR)
    processError( "Matrix::getLocalRowCopy - ", PetraError );
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : Matrix::putRow
// Purpose       : Put a row into the sparse matrix.
// Special Notes : Replace already allocated values
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
bool Matrix::putRow(int row, int length, const double *coeffs, const int *colIndices)
{
  int PetraError = oDCRSMatrix_->ReplaceGlobalValues(row, length, coeffs, colIndices);

  if (DEBUG_LINEAR)
    processError( "Matrix::putRow - ", PetraError );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::getLocalNumRows
// Purpose       : Returns the number of nonzeroes in the row.
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/24/06
//-----------------------------------------------------------------------------
int Matrix::getLocalNumRows() const
{
  return aDCRSMatrix_->NumMyRows();
}

//-----------------------------------------------------------------------------
// Function      : Matrix::getLocalRowLength
// Purpose       : Returns the number of nonzeroes in the row.
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/05/06
//-----------------------------------------------------------------------------
int Matrix::getLocalRowLength(int row) const
{
  return aDCRSMatrix_->NumMyEntries(row);
}

//-----------------------------------------------------------------------------
// Function      : Matrix::putGlobalRow
// Purpose       : Put a row into the sparse matrix.
// Special Notes : Replace already allocated values.
//                 erkeite: note; unlike putRow, this function uses the
//                 assembled matrix and global ids.
//
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/24/06
//-----------------------------------------------------------------------------
bool Matrix::putGlobalRow(int row, int length, double *coeffs, int *colIndices)
{
  int PetraError = aDCRSMatrix_->ReplaceGlobalValues(row, length, coeffs, colIndices);

  if (DEBUG_LINEAR)
    processError( "Matrix::putRow - ", PetraError );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::putLocalRow
// Purpose       : Put values into a row into the sparse matrix, using local indices.
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/29/12
//-----------------------------------------------------------------------------
bool Matrix::putLocalRow(int row, int length, const double * coeffs,
                                                   const int * colIndices)
{
  double * tmp_c = const_cast<double *>(coeffs);
  int * tmp_i = const_cast<int *>(colIndices);
  int PetraError = aDCRSMatrix_->ReplaceMyValues(row, length, tmp_c, tmp_i);

  if (DEBUG_LINEAR | DEBUG_DEVICE)
    processError( "Matrix::putLocalRow - ", PetraError );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::replaceLocalRow
// Purpose       : Replace a row in the sparse matrix.
// Special Notes : replace allocated locations
// Scope         : Public
// Creator       : Todd Coffey, 1414, Heidi Thornquist, 1437
// Creation Date : 01/31/07
//-----------------------------------------------------------------------------
void Matrix::replaceLocalRow(int row, int length, double *coeffs, int *colIndices)
{
  int PetraError = oDCRSMatrix_->ReplaceMyValues(row, length, coeffs, colIndices);

  if (DEBUG_LINEAR | DEBUG_DEVICE)
    processError( "Matrix::replaceLocalRow - ", PetraError );
}

//-----------------------------------------------------------------------------
// Function      : Matrix::getDiagonal
// Purpose       : Return the diagonal entries of the sparse matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void Matrix::getDiagonal( Vector & diagonal ) const
{
  int PetraError = aDCRSMatrix_->ExtractDiagonalCopy( *((*(diagonal.aMultiVector_))(0)) );

  if (DEBUG_LINEAR)
    processError( "Matrix::getDiagonal - ", PetraError );
}

//-----------------------------------------------------------------------------
// Function      : Matrix::replaceDiagonal
// Purpose       : Replace values of diagonal elements
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/05/03
//-----------------------------------------------------------------------------
bool Matrix::replaceDiagonal( const Vector & vec )
{
  Epetra_Vector * eVec = vec.epetraVector();
  int PetraError = aDCRSMatrix_->ReplaceDiagonalValues( *eVec );
  delete eVec;

  if (DEBUG_LINEAR)
    processError( "Matrix::replaceDiagonal - ", PetraError );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::sumIntoRow
// Purpose       : Sum values into a row into the sparse matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
bool Matrix::sumIntoRow(int row, int length, const double * coeffs,
                                                   const int * colIndices)
{
  double * tmp_c = const_cast<double *>(coeffs);
  int * tmp_i = const_cast<int *>(colIndices);
  int PetraError = oDCRSMatrix_->SumIntoGlobalValues(row, length, tmp_c, tmp_i);

  if (DEBUG_LINEAR | DEBUG_DEVICE)
    processError( "Matrix::sumIntoRow - ", PetraError );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::sumIntoLocalRow
// Purpose       : Sum values into a row into the sparse matrix, using local indices.
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/30/10
//-----------------------------------------------------------------------------
bool Matrix::sumIntoLocalRow(int row, int length, const double * coeffs,
                                                   const int * colIndices)
{
  double * tmp_c = const_cast<double *>(coeffs);
  int * tmp_i = const_cast<int *>(colIndices);
  int PetraError = oDCRSMatrix_->SumIntoMyValues(row, length, tmp_c, tmp_i);

  if (DEBUG_LINEAR | DEBUG_DEVICE)
    processError( "Matrix::sumIntoLocalRow - ", PetraError );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::returnRawEntryPointer
//
// Purpose       : This function returns a raw double* pointer to a single
//                 matrix element, specified by the local row,col indices.
//
// Special Notes : This function is much more convenient for developers of the
//                 device package than dealing with local compressed rows and
//                 the offsets required to use the bracket operators.
//                 THIS METHOD IS NOT TO BE USED BY DEVELOPERS AFTER 5/24/18
//                 THE COLUMN LID HAS BEEN PASSED IN USING LIDS FROM THE ROW MAP.
//                 THIS METHOD WILL CONVERT THE COL LID USING THE COLUMN MAP
//                 UNTIL ALL DEVICES HAVE BEEN CORRECTED.
//
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/18/2010
//-----------------------------------------------------------------------------
double * Matrix::returnRawEntryPointer (int lidRow, int lidCol)
{
  double * retPtr = &groundNode_;

  if (lidRow >= 0 && lidCol >= 0)
  {
    int num_entries = 0;
    int * indices = 0;
    double * values = 0;

    // Convert the lidCol, which is based on the row map, to be based on the column map.
    int newColLID = oDCRSMatrix_->ColMap().LID( oDCRSMatrix_->Graph().RowMap().GID(lidCol) );
    if (newColLID >= 0)
    {
      oDCRSMatrix_->ExtractMyRowView( lidRow, num_entries, values, indices );

      for( int j = 0; j < num_entries; ++j )
      {
         if (indices[j] == newColLID)
         {
           retPtr = &(values[j]);
           break;
         }
      }
    }
  }
  return retPtr;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::extractLocalRowView
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/18/2010
//-----------------------------------------------------------------------------
int Matrix::extractLocalRowView(int lidRow, int& numEntries, double*& values, int*& indices) const
{
  return oDCRSMatrix_->ExtractMyRowView( lidRow, numEntries, values, indices );
}

//-----------------------------------------------------------------------------
// Function      : Matrix::extractLocalRowView
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/18/2010
//-----------------------------------------------------------------------------
int Matrix::extractLocalRowView(int lidRow, int& numEntries, double*& values) const
{
  return oDCRSMatrix_->ExtractMyRowView( lidRow, numEntries, values);
}

//-----------------------------------------------------------------------------
// Function      : Matrix::add
// Purpose       : Sums in a matrix contribution
// Special Notes : WARNING: only works if graphs match, no checking
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/11/04
//-----------------------------------------------------------------------------
void Matrix::add( const Matrix & A )
{
  int NumRows = A.aDCRSMatrix_->NumMyRows();
  int* Indices;
  double* Values;
  int NumIndices;

  for( int i = 0; i < NumRows; ++i )
  {
    A.aDCRSMatrix_->ExtractMyRowView( i, NumIndices, Values, Indices );
    aDCRSMatrix_->SumIntoMyValues( i, NumIndices, Values, Indices );
  }
}

//-----------------------------------------------------------------------------
// Function      : Matrix::addOverlap
// Purpose       : Sums in a matrix contribution
// Special Notes : WARNING: only works if graphs match, no checking
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/11/04
//-----------------------------------------------------------------------------
void Matrix::addOverlap( const Matrix & A )
{
  int NumRows = A.oDCRSMatrix_->NumMyRows();
  int* Indices;
  double* Values;
  int NumIndices;

  for( int i = 0; i < NumRows; ++i )
  {
    A.oDCRSMatrix_->ExtractMyRowView( i, NumIndices, Values, Indices );
    oDCRSMatrix_->SumIntoMyValues( i, NumIndices, Values, Indices );
  }
}

//-----------------------------------------------------------------------------
// Function      : Matrix::linearCombo
// Purpose       : Sums in a matrix contribution
// Special Notes : WARNING: only works if graphs EXACTLY match no checking
//
//                 this = a*A + b*B
//
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 2/13/07
//-----------------------------------------------------------------------------
void Matrix::linearCombo( const double a, const Matrix & A,
                                const double b, const Matrix & B)
{
  int NumRows = (*aDCRSMatrix_).NumMyRows();

  int *aIndices, *bIndices;
  int aNumIndices, bNumIndices;
  double *aValues, *bValues;

  for( int i = 0; i < NumRows; ++i ) {
    // Get a view of the i-th row for A and B.
    A.aDCRSMatrix_->ExtractMyRowView( i, aNumIndices, aValues, aIndices );
    B.aDCRSMatrix_->ExtractMyRowView( i, bNumIndices, bValues, bIndices );

    // Add in the entries from each matrix.
    for ( int j = 0; j < aNumIndices; ++j )
      (*aDCRSMatrix_)[i][j] = a*aValues[j] + b*bValues[j];
  }
}

//-----------------------------------------------------------------------------
// Function      : Matrix::operator[]
// Purpose       : Direct access into matrix rows using local indexing
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/5/02
//-----------------------------------------------------------------------------
Matrix::bracketProxy& Matrix::operator[]( int row )
{
  proxy_.rowLID_ = row;
  return proxy_;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::operator[] const
// Purpose       : Direct access into matrix rows using local indexing
// Special Notes : const version
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/5/02
//-----------------------------------------------------------------------------
const Matrix::bracketProxy& Matrix::operator[]( int row ) const
{
  proxy_.rowLID_ = row;
  return proxy_;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::operator()
// Purpose       : Direct access into matrix values using local indexing
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/5/02
//-----------------------------------------------------------------------------
double * Matrix::operator()( int row, int col_offset )
{
  if ( row != groundLID_ && col_offset >= 0)
    return (*oDCRSMatrix_)[row]+col_offset;
  else
    return &groundNode_;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::operator() const
// Purpose       : Direct access into matrix values using local indexing
// Special Notes : const version
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/5/02
//-----------------------------------------------------------------------------
const double * Matrix::operator()( int row, int col_offset ) const
{
  if ( row != groundLID_ && col_offset >= 0)
    return (*oDCRSMatrix_)[row]+col_offset;
  else
    return &groundNode_;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::writeToFile
// Purpose       : Dumps out the sparse matrix to a file.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/19/00
//-----------------------------------------------------------------------------
void Matrix::writeToFile(const char *filename, bool useLIDs, bool mmFormat ) const
{
  if (!mmFormat)
  {
    int numProcs = aDCRSMatrix_->Comm().NumProc();
    int thisProc = aDCRSMatrix_->Comm().MyPID();
    int masterRank = 0;

    int MaxNumEntries = aDCRSMatrix_->MaxNumEntries();
    std::vector<int> Indices( MaxNumEntries );
    std::vector<double> Values( MaxNumEntries );
    int NumEntries;
    int NumMyRows = aDCRSMatrix_->NumMyRows();

    if( !aDCRSMatrix_->Filled() )
    {
      std::cerr << "Matrix: can't writeToFile unless Filled!" << std::endl;
      return;
    }

    for( int p = 0; p < numProcs; ++p )
    {
      aDCRSMatrix_->Comm().Barrier();

      if( p == thisProc )
      {
        FILE *file = NULL;

        if( masterRank == thisProc )
        {
          file = fopen( filename, "w" );
          fprintf( file, "%d\n", aDCRSMatrix_->NumGlobalNonzeros() );
        }
        else
          file = fopen( filename, "a" );

        for( int i = 0; i < NumMyRows; ++i )
        {

          if( useLIDs )
          {
            int num_entries;
            int * indices;
            double * values;
            aDCRSMatrix_->ExtractMyRowView( i, num_entries, values, indices );
            for( int j = 0; j < num_entries; ++j )
             fprintf( file, "%d %d %26.18e\n", i, indices[j], values[j] );
          }
          else
          {
            int Row = aDCRSMatrix_->Graph().RowMap().GID(i);

            aDCRSMatrix_->ExtractGlobalRowCopy( Row, MaxNumEntries, NumEntries, &Values[0], &Indices[0] );

            for( int j = 0; j < NumEntries; ++j )
             fprintf( file, "%d %d %26.18e\n", Row, Indices[j], Values[j] );
          }
        }

        fclose( file );
      }
    }
  }
  else
  {
    std::string sandiaReq = "Sandia National Laboratories is a multimission laboratory managed and operated by National Technology and\n%";
    sandiaReq += " Engineering Solutions of Sandia LLC, a wholly owned subsidiary of Honeywell International Inc. for the\n%";
    sandiaReq += " U.S. Department of Energy’s National Nuclear Security Administration under contract DE-NA0003525.\n%\n% Xyce circuit matrix.\n%%";

    EpetraExt::RowMatrixToMatrixMarketFile( filename, *aDCRSMatrix_, sandiaReq.c_str() );
  }
}

//-----------------------------------------------------------------------------
// Function      : Matrix::processError
// Purpose       : Concrete implementation which processes Petra (in this case)
//                 error codes taken from the Petra member function returns.
// Special Notes : Petra specific.  NOTE ALSO - this function is currently
//                 within the "Xyce_DEBUG_LINEAR" ifdef and so any calls to
//                 this should also be so bracketed.
// Scope         : Private
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void Matrix::processError(std::string methodMsg, int error) const
{

  const std::string PetraError("Function returned with an error.\n");

  // Process the error
  if( error < 0 )
    Report::DevelFatal0() << methodMsg + PetraError;

}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : output stream pipe operator for Matrix
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/14/00
//-----------------------------------------------------------------------------
void Matrix::printPetraObject(std::ostream &os) const
{
  if (oDCRSMatrix_ != aDCRSMatrix_)
  {
    os << *oDCRSMatrix_;
  }
  os << *aDCRSMatrix_;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::setUseTranspose
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/20/02
//-----------------------------------------------------------------------------
int Matrix::setUseTranspose (bool useTranspose)
{
  return aDCRSMatrix_->SetUseTranspose(useTranspose);
}

//-----------------------------------------------------------------------------
// Function      : Matrix::useTranspose
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/20/02
//-----------------------------------------------------------------------------
bool Matrix::useTranspose ()
{
  return aDCRSMatrix_->UseTranspose();
}

//-----------------------------------------------------------------------------
// Function      : Matrix::getOverlapColMap
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 9/6/17
//-----------------------------------------------------------------------------
N_PDS_ParMap* Matrix::getOverlapColMap( N_PDS_Comm& comm )
{
  if (!oColMap_)
    oColMap_ = new N_PDS_EpetraParMap( const_cast<Epetra_Map *>(&oDCRSMatrix_->ColMap()), comm );
  
  return oColMap_;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::getColMap
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 9/6/17
//-----------------------------------------------------------------------------
N_PDS_ParMap* Matrix::getColMap( N_PDS_Comm& comm )
{
  if (!aColMap_)
    aColMap_ = new N_PDS_EpetraParMap( const_cast<Epetra_Map *>(&aDCRSMatrix_->ColMap()), comm );
  
  return aColMap_;
}

} // namespace Linear
} // namespace Xyce

//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
#include <vector>

// ----------   Xyce Includes   ----------

#include <N_UTL_fwd.h>

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_EpetraMatrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>
#include <N_LAS_EpetraGraph.h>
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
// Function      : EpetraMatrix::~EpetraMatrix
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
EpetraMatrix::~EpetraMatrix()
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
// Function      : EpetraMatrix::EpetraMatrix
// Purpose       : Default Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
EpetraMatrix::EpetraMatrix()
: aDCRSMatrix_(0),
  oDCRSMatrix_(0),
  exporter_(0),
  offsetIndex_(0),
  aColMap_(0),
  oColMap_(0),
  overlapGraph_(0),
  baseGraph_(0),
  isOwned_(false)
{}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::EpetraMatrix
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
EpetraMatrix::EpetraMatrix( Epetra_CrsMatrix * origMatrix, bool isOwned )
: aDCRSMatrix_( origMatrix ),
  exporter_(0),
  offsetIndex_(0),
  aColMap_(0),
  oColMap_(0),
  overlapGraph_(0),
  baseGraph_(0),
  isOwned_(isOwned)
{
  oDCRSMatrix_ = aDCRSMatrix_;

  baseGraph_ = new EpetraGraph( Teuchos::rcp( const_cast<Epetra_CrsGraph*>(&(aDCRSMatrix_->Graph())), false ) );
  overlapGraph_ = baseGraph_;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::EpetraMatrix
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/21/02
//-----------------------------------------------------------------------------
EpetraMatrix::EpetraMatrix( const Graph* overlapGraph,
                            const Graph* baseGraph )
: aDCRSMatrix_(0),
  oDCRSMatrix_(0),
  exporter_(0),
  offsetIndex_(0),
  aColMap_(0),
  oColMap_(0),
  overlapGraph_(0),
  baseGraph_(0),
  isOwned_(true)
{
  const EpetraGraph* e_overlapGraph = dynamic_cast<const EpetraGraph *>( overlapGraph );
  const EpetraGraph* e_baseGraph = dynamic_cast<const EpetraGraph *>( baseGraph );

  if ( baseGraph!= overlapGraph )
  {
    oDCRSMatrix_ = new Epetra_CrsMatrix( Copy, *(e_overlapGraph->epetraObj()) );

    // Get ground node, if there is one.
    groundLID_ = overlapGraph->globalToLocalRowIndex( -1 );

    aDCRSMatrix_ = new Epetra_CrsMatrix( Copy, *(e_baseGraph->epetraObj()) );
    exporter_ = new Epetra_Export( e_overlapGraph->epetraObj()->RowMap(), e_baseGraph->epetraObj()->RowMap() );
    offsetIndex_ = new Epetra_OffsetIndex( *(e_overlapGraph->epetraObj()), *(e_baseGraph->epetraObj()), *exporter_ );
  }
  else
  {
    aDCRSMatrix_ = new Epetra_CrsMatrix( Copy, *(e_baseGraph->epetraObj()) );
    oDCRSMatrix_ = aDCRSMatrix_;
  }

  overlapGraph_ = overlapGraph->cloneCopy();
  baseGraph_ = baseGraph->cloneCopy();
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::fillComplete
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/29/03
//-----------------------------------------------------------------------------
void EpetraMatrix::fillComplete()
{
  if( exporter_ )
  {
    aDCRSMatrix_->Export( *oDCRSMatrix_, *exporter_, Add, offsetIndex_ );
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::matvec
// Purpose       : Sparse-matrix vector multiply - multivector version.  This
//                 function forms the product y = Ax where x and y are
//                 multivectors.  If transA is true, multiply by the transpose
//                 of matrix, otherwise just use matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void EpetraMatrix::matvec(bool transA, const MultiVector &x,
                          MultiVector &y)
{
  const EpetraVectorAccess* e_x = dynamic_cast<const EpetraVectorAccess *>( &x );
  EpetraVectorAccess* e_y = dynamic_cast<EpetraVectorAccess *>( &y );
  int PetraError = aDCRSMatrix_->Multiply(transA, e_x->epetraObj(), e_y->epetraObj());

  if (DEBUG_LINEAR)
    processError( "EpetraMatrix::matvec - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::put
// Purpose       : Put function for the sparse-matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void EpetraMatrix::put( double s )
{
  if ( exporter_ )
  {
    aDCRSMatrix_->PutScalar(s);
  }
  oDCRSMatrix_->PutScalar(s);
  groundNode_ = s;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::scale
// Purpose       : Scale function for the sparse-matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void EpetraMatrix::scale(double scaleFactor)
{
  if ( exporter_ )
  {
    aDCRSMatrix_->Scale(scaleFactor);
  }
  oDCRSMatrix_->Scale(scaleFactor);
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::getRowLength
// Purpose       : Returns the number of nonzeroes in the row.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
int EpetraMatrix::getRowLength(int row) const
{
  return aDCRSMatrix_->NumGlobalEntries(row);
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::getLocalRowView
// Purpose       : Returns row coefficients and associated column indices.
// Special Notes : Uses Petra's ExtractRowView which does not require user
//               : to setup space.
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
int EpetraMatrix::getLocalRowView(int lidRow, int& numEntries, double*& values, int*& indices) const
{
  return aDCRSMatrix_->ExtractMyRowView(lidRow, numEntries, values, indices);
}


//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::getRowCopy
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/05/06
//-----------------------------------------------------------------------------
void EpetraMatrix::getRowCopy
  (int row, int length, int & numEntries, double *coeffs, int *colIndices) const
{
  int PetraError = aDCRSMatrix_->ExtractGlobalRowCopy
       (row, length, numEntries, coeffs, colIndices);

  if (DEBUG_LINEAR)
    processError( "EpetraMatrix::getRowCopy - ", PetraError );
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::getLocalRowCopy
// Purpose       :
// Special Notes :
//               :
//               :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/05/06
//-----------------------------------------------------------------------------
void EpetraMatrix::getLocalRowCopy
  (int row, int length, int & numEntries, double *coeffs, int *colIndices) const
{
  int PetraError = aDCRSMatrix_->ExtractMyRowCopy
       (row, length, numEntries, coeffs, colIndices);

  if (DEBUG_LINEAR)
    processError( "EpetraMatrix::getLocalRowCopy - ", PetraError );
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::putRow
// Purpose       : Put a row into the sparse matrix.
// Special Notes : Replace already allocated values
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
bool EpetraMatrix::putRow(int row, int length, const double *coeffs, const int *colIndices)
{
  int PetraError = oDCRSMatrix_->ReplaceGlobalValues(row, length, coeffs, colIndices);

  if (DEBUG_LINEAR)
    processError( "EpetraMatrix::putRow - ", PetraError );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::getLocalNumRows
// Purpose       : Returns the number of rows on this processor.
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/24/06
//-----------------------------------------------------------------------------
int EpetraMatrix::getLocalNumRows() const
{
  return aDCRSMatrix_->NumMyRows();
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::getLocalNumRowsOverlap
// Purpose       : Returns the number of rows on this processor.
// Special Notes : This is for the overlapped matrix
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/24/06
//-----------------------------------------------------------------------------
int EpetraMatrix::getLocalNumRowsOverlap() const
{
  return oDCRSMatrix_->NumMyRows();
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::getNumRows
// Purpose       : Returns the number of rows on all processors.
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/24/06
//-----------------------------------------------------------------------------
int EpetraMatrix::getNumRows() const
{
  return aDCRSMatrix_->NumGlobalRows();
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::getLocalRowLength
// Purpose       : Returns the number of nonzeroes in the row.
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/05/06
//-----------------------------------------------------------------------------
int EpetraMatrix::getLocalRowLength(int row) const
{
  return aDCRSMatrix_->NumMyEntries(row);
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::sumIntoLocalRow
// Purpose       : Sum values into a row into the sparse matrix, using local indices.
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/30/10
//-----------------------------------------------------------------------------
bool EpetraMatrix::addIntoLocalRow(int row, int length, const double * coeffs,
                                                   const int * colIndices)
{
  int PetraError = aDCRSMatrix_->SumIntoMyValues(row, length, coeffs, colIndices);

  if (DEBUG_LINEAR | DEBUG_DEVICE)
    processError( "EpetraMatrix::addIntoLocalRow - ", PetraError );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::putLocalRow
// Purpose       : Put values into a row into the sparse matrix, using local indices.
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/29/12
//-----------------------------------------------------------------------------
bool EpetraMatrix::putLocalRow(int row, int length, const double * coeffs,
                                                   const int * colIndices)
{
  int PetraError = aDCRSMatrix_->ReplaceMyValues(row, length, coeffs, colIndices);

  if (DEBUG_LINEAR | DEBUG_DEVICE)
    processError( "EpetraMatrix::putLocalRow - ", PetraError );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::replaceLocalRow
// Purpose       : Replace a row in the sparse matrix.
// Special Notes : replace allocated locations
// Scope         : Public
// Creator       : Todd Coffey, 1414, Heidi Thornquist, 1437
// Creation Date : 01/31/07
//-----------------------------------------------------------------------------
void EpetraMatrix::replaceLocalRow(int row, int length, double *coeffs, int *colIndices)
{
  int PetraError = oDCRSMatrix_->ReplaceMyValues(row, length, coeffs, colIndices);

  if (DEBUG_LINEAR | DEBUG_DEVICE)
    processError( "EpetraMatrix::replaceLocalRow - ", PetraError );
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::getDiagonal
// Purpose       : Return the diagonal entries of the sparse matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void EpetraMatrix::getDiagonal( Vector & diagonal ) const
{
  EpetraVectorAccess* e_diag = dynamic_cast<EpetraVectorAccess *>( &diagonal );
  Epetra_Vector * ediag = e_diag->epetraObj()(0);
  int PetraError = aDCRSMatrix_->ExtractDiagonalCopy( *ediag );

  if (DEBUG_LINEAR)
    processError( "EpetraMatrix::getDiagonal - ", PetraError );
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::replaceDiagonal
// Purpose       : Replace values of diagonal elements
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/05/03
//-----------------------------------------------------------------------------
bool EpetraMatrix::replaceDiagonal( const Vector & vec )
{
  const EpetraVectorAccess* e_diag = dynamic_cast<const EpetraVectorAccess *>( &vec );
  const Epetra_Vector * eVec = e_diag->epetraObj()(0);
  int PetraError = aDCRSMatrix_->ReplaceDiagonalValues( *eVec );

  if (DEBUG_LINEAR)
    processError( "EpetraMatrix::replaceDiagonal - ", PetraError );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::sumIntoLocalRow
// Purpose       : Sum values into a row into the sparse matrix, using local indices.
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/30/10
//-----------------------------------------------------------------------------
bool EpetraMatrix::sumIntoLocalRow(int row, int length, const double * coeffs,
                                                   const int * colIndices)
{
  int PetraError = oDCRSMatrix_->SumIntoMyValues(row, length, coeffs, colIndices);

  if (DEBUG_LINEAR | DEBUG_DEVICE)
    processError( "EpetraMatrix::sumIntoLocalRow - ", PetraError );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::returnRawEntryPointer
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
double * EpetraMatrix::returnRawEntryPointer (int lidRow, int lidCol)
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
// Function      : EpetraMatrix::extractLocalRowView
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/18/2010
//-----------------------------------------------------------------------------
int EpetraMatrix::extractLocalRowView(int lidRow, int& numEntries, double*& values, int*& indices) const
{
  return oDCRSMatrix_->ExtractMyRowView( lidRow, numEntries, values, indices );
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::add
// Purpose       : Sums in a matrix contribution
// Special Notes : WARNING: only works if graphs match, no checking
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/11/04
//-----------------------------------------------------------------------------
void EpetraMatrix::add( const Matrix & A )
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
// Function      : EpetraMatrix::addOverlap
// Purpose       : Sums in a matrix contribution
// Special Notes : WARNING: only works if graphs match, no checking
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/11/04
//-----------------------------------------------------------------------------
void EpetraMatrix::addOverlap( const Matrix & A )
{
  int NumRows = A.getLocalNumRowsOverlap();
  int* Indices;
  double* Values;
  int NumIndices;

  for( int i = 0; i < NumRows; ++i )
  {
    A.extractLocalRowView( i, NumIndices, Values, Indices );
    oDCRSMatrix_->SumIntoMyValues( i, NumIndices, Values, Indices );
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::linearCombo
// Purpose       : Sums in a matrix contribution
// Special Notes : WARNING: only works if graphs EXACTLY match no checking
//
//                 this = a*A + b*B
//
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 2/13/07
//-----------------------------------------------------------------------------
void EpetraMatrix::linearCombo( const double a, const Matrix & A,
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
// Function      : EpetraMatrix::operator()
// Purpose       : Direct access into matrix values using local indexing
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/5/02
//-----------------------------------------------------------------------------
double * EpetraMatrix::operator()( int row, int col_offset )
{
  if ( row != groundLID_ && col_offset >= 0)
    return (*oDCRSMatrix_)[row]+col_offset;
  else
    return &groundNode_;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::operator() const
// Purpose       : Direct access into matrix values using local indexing
// Special Notes : const version
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/5/02
//-----------------------------------------------------------------------------
const double * EpetraMatrix::operator()( int row, int col_offset ) const
{
  if ( row != groundLID_ && col_offset >= 0)
    return (*oDCRSMatrix_)[row]+col_offset;
  else
    return &groundNode_;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::writeToFile
// Purpose       : Dumps out the sparse matrix to a file.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/19/00
//-----------------------------------------------------------------------------
void EpetraMatrix::writeToFile(const char *filename, bool useLIDs, bool mmFormat ) const
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
// Function      : EpetraMatrix::processError
// Purpose       : Concrete implementation which processes Petra (in this case)
//                 error codes taken from the Petra member function returns.
// Special Notes : Petra specific.  NOTE ALSO - this function is currently
//                 within the "Xyce_DEBUG_LINEAR" ifdef and so any calls to
//                 this should also be so bracketed.
// Scope         : Private
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void EpetraMatrix::processError(std::string methodMsg, int error) const
{

  const std::string PetraError("Function returned with an error.\n");

  // Process the error
  if( error < 0 )
    Report::DevelFatal0() << methodMsg + PetraError;

}

//-----------------------------------------------------------------------------
// Function      : print
// Purpose       : print method for Matrix
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/14/00
//-----------------------------------------------------------------------------
void EpetraMatrix::print(std::ostream &os) const
{
  if (oDCRSMatrix_ != aDCRSMatrix_)
  {
    os << *oDCRSMatrix_;
  }
  os << *aDCRSMatrix_;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::setUseTranspose
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/20/02
//-----------------------------------------------------------------------------
int EpetraMatrix::setUseTranspose (bool useTranspose)
{
  return aDCRSMatrix_->SetUseTranspose(useTranspose);
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::useTranspose
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/20/02
//-----------------------------------------------------------------------------
bool EpetraMatrix::useTranspose () const
{
  return aDCRSMatrix_->UseTranspose();
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::getOverlapColMap
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 9/6/17
//-----------------------------------------------------------------------------
const Parallel::ParMap* EpetraMatrix::getOverlapColMap( const Parallel::Communicator& comm )
{
  if (!oColMap_)
    oColMap_ = new Parallel::EpetraParMap( &oDCRSMatrix_->ColMap(), comm );
  
  return oColMap_;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMatrix::getColMap
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 9/6/17
//-----------------------------------------------------------------------------
const Parallel::ParMap* EpetraMatrix::getColMap( const Parallel::Communicator& comm ) const
{
  if (!aColMap_)
    aColMap_ = new Parallel::EpetraParMap( &aDCRSMatrix_->ColMap(), comm );
  
  return aColMap_;
}

} // namespace Linear
} // namespace Xyce

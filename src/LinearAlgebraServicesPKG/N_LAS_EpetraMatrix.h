//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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

//-----------------------------------------------------------------------------
//
// Purpose        : Epetra implementation of the abstract interface to sparse
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
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_EpetraMatrix_h
#define Xyce_N_LAS_EpetraMatrix_h

// ---------- Standard Includes ----------
#include <string>

// ----------   Xyce Includes   ----------
#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_EpetraHelpers.h>

class Epetra_CrsMatrix;

class Epetra_Export;
class Epetra_OffsetIndex;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : EpetraMatrix
// Purpose       : Abstract interface to sparse matrix type.
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Compuational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
class EpetraMatrix : public Matrix, public EpetraMatrixAccess
{

public:

  // Default constructor 
  EpetraMatrix();

  //Constructors
  EpetraMatrix( const Graph* overlapGraph,
                const Graph* baseGraph );

  //Constructor from an existing Epetra_CrsMatrix (makes copy of origMatrix)
  EpetraMatrix( Epetra_CrsMatrix * origMatrix, bool isOwned = true );

  //Destructor
  virtual ~EpetraMatrix();

  // This function needs to be invoked for a transpose solve.
  int setUseTranspose (bool useTranspose);
  bool useTranspose () const;

  //Accumulate off processor fill contributions if necessary
  void fillComplete();

  // --------------------------------------------------------------------------------
  // Overlapped matrix methods
  // - Underlying this class is both an overlapped matrix and assembled matrix.
  // - The overlapped matrix is used during the assembly phase for device model loading.
  // - So these methods should only be used before fillComplete is called.
  // --------------------------------------------------------------------------------

  // Get number of rows in the overlapped matrix
  int getLocalNumRowsOverlap() const;

  // Add in a matrix contribution
  void addOverlap( const Matrix& A );

  // Put a set of values into a row
  bool putRow(int row, int length, const double * coeffs, const int * colIndices);

  // Replace a set of values into a row
  void replaceLocalRow(int row, int length, double * coeffs, int * colIndices);

  // Sum values into a row into the sparse matrix, using local indices
  bool sumIntoLocalRow(int row, int length, const double * coeffs, const int * colIndices);

  // get a pointer to the compressed local row.
  int extractLocalRowView(int lidRow, int& numEntries, double*& values, int*& indices) const;

  // Return a pointer to a single row, col element.
  double * returnRawEntryPointer (int lidRow, int lidCol);

  // Direct access into matrix rows and columns using local indexing, with column offset.
  double * operator()(int row, int col_offset);
  const double * operator()(int row, int col_offset) const;

  // Get column map for overlapped matrix
  const Parallel::ParMap* getOverlapColMap( const Parallel::Communicator& comm );

  // Get graphs for overlapped matrix
  const Graph* getOverlapGraph() const { return overlapGraph_; }

  // --------------------------------------------------------------------------------
  // Assembled matrix methods
  // - Underlying this class is both an overlapped matrix and assembled matrix.
  // - The assembled matrix is used after device model loading by the time integrator, nonlinear
  //   solver, and linear solver.
  // - Any matrix used after the loadDAEMatrices will be assembled, so use these 
  //   computational methods.
  // --------------------------------------------------------------------------------

  // Add in a matrix contribution
  void add( const Matrix & A );

  // Sparse-matrix vector multiply - multivector version.  If transA is true,
  // multiply by the transpose of matrix, otherwise just use matrix.
  void matvec(bool transA, const MultiVector & x, MultiVector & y) const;

  // Performs the operation this <- a*A + b*B
  void linearCombo ( const double a, const Matrix & A, 
                     const double b, const Matrix & B);

  // Get the matrix diagonal (stored in a Vector)
  void getDiagonal(Vector & diagonal) const;

  // Replace a vector's values into the matrix diagonal
  bool replaceDiagonal(const Vector & vec);

  // Get a row's length (nonzero) using global and local row ids.
  int getRowLength(int row) const;
  int getLocalRowLength(int row) const;

  // Get the non-zero values in a row, using local indices
  // NOTE:  The global indices version of this method is not available after fillComplete is called.
  int getLocalRowView(int lidRow, int& numEntries, double*& values, int*& indices) const;

  // Get the non-zero values in a row
  void getRowCopy(int row, int length, int & numEntries, double * coeffs, int * colIndices) const;
  void getLocalRowCopy(int row, int length, int & numEntries, double * coeffs, int * colIndices) const;

  int getNumRows() const;
  int getLocalNumRows() const;

  // Sum values into a row into the sparse matrix, using local indices, without overlap contributions
  bool addIntoLocalRow(int row, int length, const double * coeffs, const int * colIndices);

  // Put a set of values into a row, using local indices
  bool putLocalRow(int row, int length, const double * coeffs, const int * colIndices);

  // Output the matrix to a file
  void writeToFile(const char * filename, bool useLIDs = false, bool mmFormat=false ) const;

  // Get column map for assembled matrix
  const Parallel::ParMap* getColMap( const Parallel::Communicator& comm ) const;

  // Get graph for assembled matrix
  const Graph* getGraph() const { return baseGraph_; }

  Epetra_CrsMatrix & epetraObj() { return *aDCRSMatrix_; }
  const Epetra_CrsMatrix & epetraObj() const { return *aDCRSMatrix_; }

  // --------------------------------------------------------------------------------
  // Overlapped/assembled matrix methods
  // - Underlying this class is both an overlapped matrix and assembled matrix.
  // - These methods perform the same operation on both matrices, so fillComplete 
  //   is not necessary to migrate the data from overlapped to assembled in parallel.
  // --------------------------------------------------------------------------------

  // Initialize matrix values to s
  void put(double s);

  // Scale the matrix
  void scale(double scaleFactor);

  // Print the underlying objects
  void print(std::ostream &os) const;

protected:

  // Pointer the Petra multi-vector object.
  Epetra_CrsMatrix * aDCRSMatrix_;

  // Overlapped version of matrix
  Epetra_CrsMatrix * oDCRSMatrix_;

  // Importing Tools
  Epetra_Export * exporter_;
  Epetra_OffsetIndex * offsetIndex_;

  // Column maps, assembled and overlapped.
  mutable const Parallel::ParMap *aColMap_, *oColMap_;

  // Graphs, assembled and overlapped.
  const Graph *overlapGraph_, *baseGraph_;

  // isOwned flag
  bool isOwned_;

private:

  // Copy constructor (private)
  EpetraMatrix(const Matrix & right);
  // Assignment operator (private)
  EpetraMatrix & operator = (const Matrix & right);

  bool operator == (const EpetraMatrix & right) const;
  bool operator != (const EpetraMatrix & right) const;

  // Process library error codes.
  void processError(std::string methodMsg, int error) const;
};

} // namespace Linear
} // namespace Xyce

#endif

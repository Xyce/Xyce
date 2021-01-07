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

//-----------------------------------------------------------------------------
//
// Purpose        : Specification file for the Abstract interface to sparse
//                  block matrix type.
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/12/04
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_EpetraBlockMatrix_h
#define Xyce_N_LAS_EpetraBlockMatrix_h

// ---------- Standard Includes ----------
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_LAS_fwd.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_EpetraMatrix.h>
#include <N_LAS_EpetraHelpers.h>

// ----------  Other Includes   ----------
#include <Teuchos_RCP.hpp>
using Teuchos::RCP;

#include <Epetra_CrsMatrix.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : BlockMatrix
// Purpose       : Abstract interface to sparse block matrix type.
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
class EpetraBlockMatrix : public BlockMatrix, public EpetraMatrixAccess
{
 public:

  EpetraBlockMatrix( int size,
                     int offset,
                     const std::vector< std::vector<int> > & blockColumns,
                     const Graph* globalGraph,
                     const Graph* subBlockGraph,
                     int augmentCount = 0 );

  //Destructor
  virtual ~EpetraBlockMatrix();

  // --------------------------------------------------------------------------------
  // Linear::Matrix assembled matrix methods
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
  void matvec(bool transA, const MultiVector & x, MultiVector & y);

  // Performs the operation this <- a*A + b*B
  void linearCombo ( const double a, const Matrix & A,
                     const double b, const Matrix & B);

  // Get the matrix diagonal (stored in an Vector)
  void getDiagonal(Vector & diagonal) const;

  // Replace a vector's values into the matrix diagonal
  bool replaceDiagonal(const Vector & vec);

  // Get a row's length (nonzero) using global and local row ids.
  int getRowLength(int row) const
  { return aDCRSMatrix_->NumGlobalEntries(row); }
  int getLocalRowLength(int row) const
  { return aDCRSMatrix_->NumMyEntries(row); }

  // Get the non-zero values in a row, using local indices
  // NOTE:  The global indices version of this method is not available after fillComplete is called.
  int getLocalRowView(int lidRow, int& numEntries, double*& values, int*& indices) const;

  // Get the non-zero values in a row
  void getRowCopy(int row, int length, int & numEntries, double * coeffs, int * colIndices) const;
  void getLocalRowCopy(int row, int length, int & numEntries, double * coeffs, int * colIndices) const;

  int getNumRows() const
  { return aDCRSMatrix_->NumGlobalRows(); }
  int getLocalNumRows() const
  { return aDCRSMatrix_->NumMyRows(); }

  // Sum values into a row into the sparse matrix, using local indices, without overlap contributions
  bool addIntoLocalRow(int row, int length, const double * coeffs, const int * colIndices);

  // Put a set of values into a row, using local indices
  bool putLocalRow(int row, int length, const double * coeffs, const int * colIndices);

  // Output the matrix to a file
  void writeToFile(const char * filename, bool useLIDs = false, bool mmFormat=false ) const {}

  // Get column map for assembled matrix
  const Parallel::ParMap* getColMap( const Parallel::Communicator& comm ) const;

  // Get graph for assembled matrix
  const Graph* getGraph() const
  { return baseGraph_; }

  // This function needs to be invoked for a transpose solve.
  int setUseTranspose (bool useTranspose)
  { return aDCRSMatrix_->SetUseTranspose(useTranspose); }
  bool useTranspose () const
  { return aDCRSMatrix_->UseTranspose(); }

  // Put function for the block sparse-matrix.
  void put(double s);

  // Scale the matrix
  void scale(double scaleFactor)
  { aDCRSMatrix_->Scale(scaleFactor); }

  void print(std::ostream &os) const;

  // --------------------------------------------------------------------------------
  // Linear::Matrix overlapped matrix methods
  // --------------------------------------------------------------------------------

  // Put a set of values into a row
  bool putRow(int row, int length, const double * coeffs, const int * colIndices);

  // Return a pointer to a single row, col element.
  double * returnRawEntryPointer (int lidRow, int lidCol) { return 0; }

  //Accumulate off processor fill contributions if necessary
  //NOTE:  This method sums the contributions from the overlapped matrix and places it in
  //       the assembled matrix.
  void fillComplete() {}

  // --------------------------------------------------------------------------------
  // Linear::BlockMatrix methods
  // --------------------------------------------------------------------------------

  // Block Access
  Matrix & block( int row, int col );

  int blockSize() const
  { return blockSize_; }
  
  int numBlockRows() const
  { return numBlockRows_; }

  // Replace the entries of an augmented row using the row GID.
  void replaceAugmentedRow(int rowGID, int length, double * coeffs, int * colIndices); 

  void replaceAugmentedColumn(int augmentedColumn, const BlockVector & vec);

  // Assemble global matrix with blocks
  // NOTE:  The global matrix is not always a view of the local matrix, so this function ensures
  // that the values are sync'ed up.  Call this before using the global matrix for computations.
  void assembleGlobalMatrix();
 
  // Underlying object access
  Epetra_CrsMatrix & epetraObj() { return *aDCRSMatrix_; }
  const Epetra_CrsMatrix & epetraObj() const { return *aDCRSMatrix_; }

 private:

  // Process library error codes.
  void processError(std::string methodMsg, int error) const;

  // Pointer the Petra multi-vector object.
  Epetra_CrsMatrix * aDCRSMatrix_;

  // Column maps, assembled and overlapped.
  mutable const Parallel::ParMap *aColMap_;

  // Graphs, assembled and overlapped.
  const Graph *baseGraph_;

  // isOwned flag
  bool isOwned_;

  bool blocksViewGlobalMat_;
  const int blockSize_;
  const int offset_;
  const int numBlockRows_;
  const int augmentCount_;

  std::vector<int> augmentGIDs_, baseNumCols_, baseIndices_;
  const std::vector< std::vector<int> > cols_;
  std::vector< std::vector<Teuchos::RCP<EpetraMatrix> > > blocks_;
};

} // namespace Linear
} // namespace Xyce

#endif // Xyce_N_LAS_EpetraBlockMatrix_h

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

//-----------------------------------------------------------------------------
//
// Purpose        : Specification file for the Abstract interface to sparse
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

#ifndef Xyce_N_LAS_Matrix_h
#define Xyce_N_LAS_Matrix_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : Matrix
// Purpose       : Abstract interface to sparse matrix type.
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Compuational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
class Matrix
{

public:

  // Default constructor 
  Matrix()
  : proxy_( 0, *this ),
    groundLID_(-1),
    groundNode_(0.0)
  {}

  struct bracketProxy {

    bracketProxy( int row, Matrix& thisMatrix )
    : rowLID_( row ),
      matrix_( thisMatrix ) 
    {}

    ~bracketProxy() {}

    int rowLID_;
    Matrix& matrix_;

    double& operator[] (int col_offset);

    const double& operator[] (int col_offset) const;
  };

  //Destructor
  virtual ~Matrix() {}

  // This function needs to be invoked for a transpose solve.
  virtual int setUseTranspose (bool useTranspose) = 0;
  virtual bool useTranspose () const = 0;

  //Accumulate off processor fill contributions if necessary
  //NOTE:  This method sums the contributions from the overlapped matrix and places it in
  //       the assembled matrix.
  virtual void fillComplete() {};

  // --------------------------------------------------------------------------------
  // Overlapped matrix methods
  // - Underlying this class is possibly both an overlapped matrix and assembled matrix.
  // - The overlapped matrix is used during the assembly phase for device model loading.
  // - So these methods should only be used before fillComplete is called.
  // - These methods can often default to assembled matrix methods.
  // --------------------------------------------------------------------------------
 
  // Get number of rows in the overlapped matrix  
  virtual int getLocalNumRowsOverlap() const
  { return this->getLocalNumRows(); } 

  // Add in a matrix contribution
  virtual void addOverlap( const Matrix & A )
  { this->add( A ); }

  // Put a set of values into a row
  virtual bool putRow(int row, int length, const double * coeffs, const int * colIndices) = 0;

  // Replace a set of values into a row
  virtual void replaceLocalRow(int row, int length, double * coeffs, int * colIndices)
  { this->putLocalRow( row, length, coeffs, colIndices ); }

  // Sum values into a row into the sparse matrix, using local indices
  virtual bool sumIntoLocalRow(int row, int length, const double * coeffs, const int * colIndices)
  { return this->addIntoLocalRow( row, length, coeffs, colIndices ); }

  // get a pointer to the compressed local row.
  virtual int extractLocalRowView(int lidRow, int& numEntries, double*& values, int*& indices) const
  { return this->getLocalRowView( lidRow, numEntries, values, indices ); }

  // Return a pointer to a single row, col element.
  virtual double * returnRawEntryPointer (int lidRow, int lidCol) = 0;

  // Direct access into matrix rows using local indexing.
  bracketProxy& operator[]( int row );

  const bracketProxy& operator[]( int row ) const;

  // Direct access into matrix rows and columns using local indexing, with column offset.
  virtual double * operator()(int row, int col_offset) { return &groundNode_; }
  virtual const double * operator()(int row, int col_offset) const { return &groundNode_; }

  // Get column map for overlapped matrix
  virtual const Parallel::ParMap* getOverlapColMap( const Parallel::Communicator& comm )
  { return this->getColMap( comm ); }

  // Get graphs for overlapped matrix
  virtual const Graph* getOverlapGraph() const
  { return this->getGraph(); }

  // --------------------------------------------------------------------------------
  // Assembled matrix methods
  // - Underlying this class is both an overlapped matrix and assembled matrix.
  // - The assembled matrix is used after device model loading by the time integrator, nonlinear
  //   solver, and linear solver.
  // - Any matrix used after the loadDAEMatrices will be assembled, so use these 
  //   computational methods.
  // --------------------------------------------------------------------------------

  // Add in a matrix contribution
  virtual void add( const Matrix & A ) = 0;

  // Sparse-matrix vector multiply - multivector version.  If transA is true,
  // multiply by the transpose of matrix, otherwise just use matrix.
  virtual void matvec(bool transA, const MultiVector & x, MultiVector & y) const = 0;

  // Performs the operation this <- a*A + b*B
  virtual void linearCombo ( const double a, const Matrix & A, 
                             const double b, const Matrix & B) = 0;

  // Get the matrix diagonal (stored in an Vector)
  virtual void getDiagonal(Vector & diagonal) const = 0;

  // Replace a vector's values into the matrix diagonal
  virtual bool replaceDiagonal(const Vector & vec) = 0;

  // Get a row's length (nonzero) using global and local row ids.
  virtual int getRowLength(int row) const = 0;
  virtual int getLocalRowLength(int row) const = 0;

  // Get the non-zero values in a row, using local indices
  // NOTE:  The global indices version of this method is not available after fillComplete is called.
  virtual int getLocalRowView(int lidRow, int& numEntries, double*& values, int*& indices) const = 0;

  // Get the non-zero values in a row
  virtual void getRowCopy(int row, int length, int & numEntries, double * coeffs, int * colIndices) const = 0;
  virtual void getLocalRowCopy(int row, int length, int & numEntries, double * coeffs, int * colIndices) const = 0;

  virtual int getNumRows() const = 0;
  virtual int getLocalNumRows() const = 0;

  // Sum values into a row into the sparse matrix, using local indices, without overlap contributions
  virtual bool addIntoLocalRow(int row, int length, const double * coeffs, const int * colIndices) = 0;

  // Put a set of values into a row, using local indices
  virtual bool putLocalRow(int row, int length, const double * coeffs, const int * colIndices) = 0;

  // Output the matrix to a file
  virtual void writeToFile(const char * filename, bool useLIDs = false, bool mmFormat=false ) const = 0;

  // Get column map for assembled matrix
  virtual const Parallel::ParMap* getColMap( const Parallel::Communicator& comm ) const = 0;

  // Get graph for assembled matrix
  virtual const Graph* getGraph() const = 0;

  // --------------------------------------------------------------------------------
  // Overlapped/assembled matrix methods
  // - Underlying this class is both an overlapped matrix and assembled matrix.
  // - These methods perform the same operation on both matrices, so fillComplete 
  //   is not necessary to migrate the data from overlapped to assembled in parallel.
  // --------------------------------------------------------------------------------

  // Initialize matrix values to s
  virtual void put(double s) = 0;

  // Scale the matrix
  virtual void scale(double scaleFactor) = 0;

  // Print the underlying objects
  virtual void print(std::ostream &os) const = 0;

protected:

  // Dummy variable for loading ground node contributions.
  mutable bracketProxy proxy_;
  int groundLID_;
  double groundNode_;

private:

  // Copy constructor (private)
  Matrix(const Matrix & right);
  // Assignment operator (private)
  Matrix & operator = (const Matrix & right);

  bool operator == (const Matrix & right) const;
  bool operator != (const Matrix & right) const;

};

//-----------------------------------------------------------------------------
// Function      : Matrix::operator[]
// Purpose       : Direct access into matrix rows using local indexing
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 9/5/02
//-----------------------------------------------------------------------------
inline Matrix::bracketProxy& Matrix::operator[]( int row )
{
  proxy_.rowLID_ = row;
  return proxy_;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::operator[] const
// Purpose       : Direct access into matrix rows using local indexing
// Special Notes : const version
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 9/5/02
//-----------------------------------------------------------------------------
inline const Matrix::bracketProxy& Matrix::operator[]( int row ) const
{
  proxy_.rowLID_ = row;
  return proxy_;
}

//-----------------------------------------------------------------------------
// Function      : Matrix::bracketProxy::operator[] 
// Purpose       : Direct access into matrix rows using local indexing
// Special Notes : 
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 9/5/02
//-----------------------------------------------------------------------------
inline double& Matrix::bracketProxy::operator[] (int col_offset)
{
  return *(matrix_( rowLID_, col_offset ));
}

//-----------------------------------------------------------------------------
// Function      : Matrix::bracketProxy::operator[] const
// Purpose       : Direct access into matrix rows using local indexing
// Special Notes : const version
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 9/5/02
//-----------------------------------------------------------------------------
inline const double& Matrix::bracketProxy::operator[] (int col_offset) const
{
  return *(matrix_( rowLID_, col_offset ));
}

} // namespace Linear
} // namespace Xyce

#endif

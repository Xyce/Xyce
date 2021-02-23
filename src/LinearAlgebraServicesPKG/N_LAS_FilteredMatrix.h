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

//-----------------------------------------------------------------------------
//
// Purpose        : Specification file for the Abstract interface to sparse
//                  matrix type.
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, Sandia National Labs
//
// Creation Date  : 
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_FilteredMatrix_h
#define Xyce_N_LAS_FilteredMatrix_h

// ---------- Standard Includes ----------
#include <string>
#include <vector>
#include <utility>

// ----------   Xyce Includes   ----------
#include <N_LAS_fwd.h>
#include <N_NLS_fwd.h>
#include <N_LAS_Matrix.h>

#include <Teuchos_RCP.hpp>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : FilteredMatrix
// Purpose       : Abstract interface to sparse matrix type.
// Special Notes :
// Creator       : Heidi Thornquist, Sandia National Labs
// Creation Date : 
//-----------------------------------------------------------------------------
class FilteredMatrix
{
public:

  //Constructors
  FilteredMatrix( const Matrix* matrix, const Parallel::ParMap* map, 
                  bool filterOverlap = true );

  FilteredMatrix( const std::vector<int>& ptr, const std::vector<int>& indices,
                  const std::vector<double>& values, bool isCSR = true );

  //Destructor
  virtual ~FilteredMatrix() {}

  // Filter current matrix.
  bool filterMatrix( const Matrix* matrix, const Parallel::ParMap* map, bool reset = false );

  bool isEmpty() const { return (totalNZRows_ == 0); }

  void clearFilter();

  // This is used during assembly, A += alpha*this.
  void addToMatrix( Matrix& A, double alpha=1.0 ); 

  // Sparse-matrix vector multiply - multivector version, y = A*x + y.  
  void axpy(const MultiVector & x, MultiVector & y);

  // Sparse-matrix vector multiply - multivector version, y = A*x.  
  void matvec(const MultiVector & x, MultiVector & y);

  // Return the local row IDs for the rows with nonzero values.
  const std::vector<int>& getNZRows() const { return nzRows_; }

  // Return the local col IDs for the cols with nonzero values.
  const std::vector<int>& getNZCols() const { return nzCols_; }

  // Return the row pointer for the nonzero values.
  const std::vector<int>& getRowPtr() const { return rowPtr_; }

  // Return the local col IDs for the nonzero values.
  const std::vector<int>& getIndices() const { return colIndices_; }

  // Return the nonzero values.
  const std::vector<double>& getValues() const { return values_; }

  // Print the filtered matrix (non-zero) entries from the original matrix.
  void printFilteredMatrix(std::ostream &os);

private:

  // Default constructor (private)
  FilteredMatrix();

  // Copy constructor (private)
  FilteredMatrix(const FilteredMatrix & right);

  // Assignment operator (private)
  FilteredMatrix & operator = (const FilteredMatrix & right);

  bool operator == (const FilteredMatrix & right) const;
  bool operator != (const FilteredMatrix & right) const;

  // Whether the overlapped matrix is filtered.
  // NOTE:  This is used during device loading, before assembly occurs.
  bool filterOverlap_;

  // Import object, if needed.
  Teuchos::RCP<Importer> importer_;

  // PDS_ParMap object, if needed.
  Teuchos::RCP<Parallel::ParMap> targetMap_;

  // Local x vector for matvec/axpy
  Teuchos::RCP<MultiVector> targetX_;

  // Filtered matrix in CRS format.
  int totalNZRows_;
  double minValue_, maxValue_;
  std::vector<int> colIndices_, vecIndices_, rowPtr_, nzRows_, nzCols_;
  std::vector<double> values_;

};

} // namespace Linear
} // namespace Xyce

#endif

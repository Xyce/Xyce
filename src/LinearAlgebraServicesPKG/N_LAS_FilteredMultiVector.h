//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
//                  multivector type.
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

#ifndef Xyce_N_LAS_FilteredMultiVector_h
#define Xyce_N_LAS_FilteredMultiVector_h

// ---------- Standard Includes ----------
#include <string>
#include <vector>
#include <utility>

// ----------   Xyce Includes   ----------
#include <N_LAS_fwd.h>
#include <N_NLS_fwd.h>
#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : FilteredMultiVector
// Purpose       : Abstract interface to sparse multivector type.
// Special Notes :
// Creator       : Heidi Thornquist, Sandia National Labs
// Creation Date : 
//-----------------------------------------------------------------------------
class FilteredMultiVector
{
public:

  //Constructors

  // Construct with pointers, indices, and values.
  FilteredMultiVector( const std::vector<int>& ptr, const std::vector<int>& indices,
                  const std::vector<double>& values );

  // Construct with pointers and indices, but no values.
  FilteredMultiVector( const std::vector<int>& ptr, const std::vector<int>& indices );

  // Construct with number of columns, empty multivector.
  FilteredMultiVector( const int numCols );

  //Destructor
  virtual ~FilteredMultiVector() {}

  // Copy constructor 
  FilteredMultiVector(const FilteredMultiVector & right);

  // Assignment operator 
  FilteredMultiVector & operator = (const FilteredMultiVector & right);

  // Insert column.
  bool insertColumn( const std::vector<int>& indices, const std::vector<double>& values,
                     const int col );

  // Filtered input vector and place in given column.
  bool filterVector( const Vector* vector, const int col );

  // Extract nonzero values from input MultiVector.
  bool replaceValues( const std::vector<std::vector<int> >& indices, const MultiVector& V );

  bool isEmpty() const { return (colPtr_.empty()); }

  // This is used during assembly, V += alpha*this.
  void addToMultiVector( MultiVector& V, double alpha=1.0 ) const; 

  // Insert values from this sparse multivector into a dense multivector.
  void convertToMultiVector( MultiVector& V ) const;

  // Dot product each column with dense vector.
  void dotProduct( const Vector& V, std::vector<double>& dot ) const;

  // Return the col pointer for the nonzero values.
  const std::vector<int>& getColPtr() const { return colPtr_; }

  // Return the local row IDs for the nonzero values.
  const std::vector<int>& getIndices() const { return rowIndices_; }

  // Return the nonzero values.
  const std::vector<double>& getValues() const { return values_; }

private:

  // Default constructor (private)
  FilteredMultiVector();

  bool operator == (const FilteredMultiVector & right) const;
  bool operator != (const FilteredMultiVector & right) const;

  // Filtered vector in CSC format.
  std::vector<int> rowIndices_, colPtr_;
  std::vector<double> values_;

};

} // namespace Linear
} // namespace Xyce

#endif

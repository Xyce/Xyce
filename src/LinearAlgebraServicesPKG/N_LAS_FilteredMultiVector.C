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
// Creator        : Heidi Thornquist, Sandia National Labs
//
// Creation Date  : 
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
#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>
#include <N_LAS_FilteredMultiVector.h>
#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : FilteredMultiVector::FilteredMultiVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, Sandia National Labs
// Creation Date : 02/01/17
//-----------------------------------------------------------------------------
FilteredMultiVector::FilteredMultiVector(  const std::vector<int>& ptr, const std::vector<int>& indices,
                                           const std::vector<double>& values )
: rowIndices_( indices ),
  colPtr_( ptr ),
  values_( values )
{
}

//-----------------------------------------------------------------------------
// Function      : FilteredMultiVector::FilteredMultiVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, Sandia National Labs
// Creation Date : 02/01/17
//-----------------------------------------------------------------------------
FilteredMultiVector::FilteredMultiVector( const std::vector<int>& ptr, const std::vector<int>& indices )
: rowIndices_( indices ),
  colPtr_( ptr )
{
  values_.resize( rowIndices_.size() );
}

//-----------------------------------------------------------------------------
// Function      : FilteredMultiVector::FilteredMultiVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, Sandia National Labs
// Creation Date : 02/01/17
//-----------------------------------------------------------------------------
FilteredMultiVector::FilteredMultiVector( const int numCols )
{
  colPtr_.resize( numCols+1, 0 );
}

//-----------------------------------------------------------------------------
// Function      : FilteredMultiVector::FilteredMultiVector
// Purpose       : Copy Constructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, Sandia National Labs
// Creation Date : 02/01/17 
//-----------------------------------------------------------------------------
FilteredMultiVector::FilteredMultiVector(const FilteredMultiVector & right)
{
  rowIndices_ = right.rowIndices_;
  colPtr_ = right.colPtr_;
  values_ = right.values_;
}

//-----------------------------------------------------------------------------
// Function      : FilteredMultiVector::operator=
// Purpose       : Assignment Operator
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, Sandia National Labs
// Creation Date : 02/01/17
//-----------------------------------------------------------------------------
FilteredMultiVector & FilteredMultiVector::operator = (const FilteredMultiVector & right)
{
  rowIndices_ = right.rowIndices_;
  colPtr_ = right.colPtr_;
  values_ = right.values_;

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : FilteredMultiVector::insertColumn
// Purpose       : Insertion Operator
// Special Notes : Sparse column is replaced by this operation
// Scope         : Public
// Creator       : Heidi Thornquist, Sandia National Labs
// Creation Date : 02/01/17
//-----------------------------------------------------------------------------
bool FilteredMultiVector::insertColumn( const std::vector<int>& indices, 
                                        const std::vector<double>& values,
                                        const int col )
{
  int numCols = colPtr_.size()-1;

  if (col > numCols)
  {
    Report::DevelFatal().in("FilteredMultiVector::insertColumn")
      << "Filtered multivector has too few columns.";
  }
  if ( indices.size() != values.size() )
  {
    Report::DevelFatal().in("FilteredMultiVector::insertColumn")
      << "Column being inserted has a different number of indices and values.";
  }

  int numEntries = colPtr_[col+1] - colPtr_[col];
  int newEntries = indices.size();
  int diffEntries = newEntries - numEntries;  // This can be negative.
  int ptr = colPtr_[col];

  // We are replacing this column, so erase the current indices and values,
  // if there are any.
  if (numEntries) 
  {
    std::vector<int>::iterator r_it = rowIndices_.begin();
    rowIndices_.erase( r_it+ptr, r_it+(ptr+numEntries) );

    std::vector<double>::iterator v_it = values_.begin();
    values_.erase( v_it+ptr, v_it+(ptr+numEntries) );
  }

  // Now insert the new entries for this column.
  // Insert indices.
  std::vector<int>::iterator r_it = rowIndices_.begin();
  rowIndices_.insert( r_it+ptr, indices.begin(), indices.end() );
 
  // Insert values.
  std::vector<double>::iterator v_it = values_.begin();
  values_.insert( v_it+ptr, values.begin(), values.end() );

  // Update column pointers by adding the difference between the new column size and old column size.
  for (int j=col+1; j<numCols; j++)
  {
    colPtr_[j] += diffEntries;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : FilteredMultiVector::filterVector
// Purpose       : Insertion Operator
// Special Notes : Sparse column is replaced by this operation
// Scope         : Public
// Creator       : Heidi Thornquist, Sandia National Labs
// Creation Date : 02/01/17
//-----------------------------------------------------------------------------
bool FilteredMultiVector::filterVector( const Vector* vector, const int col )
{
  int numCols = colPtr_.size()-1;

  if (col > numCols)
  {
    Report::DevelFatal().in("FilteredMultiVector::filterVector")
      << "Filtered multivector has too few columns.";
  }
  
  int length = vector->localLength();

  std::vector<int> newIndices;
  std::vector<double> newValues;
  for (int i=0; i<length; i++)
  {
    if ( (*vector)[i] != 0.0 )
    {
      newIndices.push_back( i );
      newValues.push_back( (*vector)[i] );
    }
  } 

  // Now insert the column.
  return( insertColumn( newIndices, newValues, col ) );
}

//-----------------------------------------------------------------------------
// Function      : FilteredMultiVector::insertValues
// Purpose       : Replaces current filtered matrix.
// Special Notes : 
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 2/01/17
//-----------------------------------------------------------------------------
bool FilteredMultiVector::replaceValues( const std::vector<std::vector<int> >& indices, const MultiVector& V )
{
  // Clear current filtered matrix.
  colPtr_.clear();
  rowIndices_.clear();
  values_.clear();

  // Start placing indices and values as given by indices and V.
  int numCols = indices.size();
  colPtr_.resize( numCols + 1 );
  colPtr_.push_back( 0 );

  for (int j=0; j<numCols; ++j)
  {
    int count=0;
    std::vector<int>::const_iterator it = indices[j].begin();
    std::vector<int>::const_iterator it_end = indices[j].end();
    for ( ; it != it_end; ++it, ++count)
    {
      rowIndices_.push_back( *it );
      values_.push_back( *V(*it,j) );
    }
    colPtr_[j+1] = colPtr_[j] + count;   
  } 

  return true;
}

//-----------------------------------------------------------------------------
// Function      : FilteredMultiVector::addToMultiVector
// Purpose       : Sums in a filtered multivector contribution
// Special Notes : 
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 2/01/17
//-----------------------------------------------------------------------------
void FilteredMultiVector::addToMultiVector( MultiVector& V, double alpha ) const
{
  if (!isEmpty())
  {
    int numCols = colPtr_.size()-1;
   
    // Check if the multivectors are compatible. 
    if ( V.numVectors() != numCols )  
    {
      Report::DevelFatal().in("FilteredMultiVector::addToMultiVector")
        << "Filtered multivector does not have same number of columns as input MultiVector.";
    }

    for (int j=0; j<numCols; j++)
    {
      for (int ptr=colPtr_[j]; ptr<colPtr_[j+1]; ptr++)
      {
        *V(rowIndices_[ptr],j) += alpha*values_[ptr];
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : FilteredMultiVector::convertToMultiVector
// Purpose       : Insert values from sparse multivector into a dense multivector.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 2/01/17
//-----------------------------------------------------------------------------
void FilteredMultiVector::convertToMultiVector( MultiVector& V ) const
{
  // Initialize vector to zero, then add sparse multivector to it.
  V.putScalar( 0.0 );
  addToMultiVector( V );
}

//-----------------------------------------------------------------------------
// Function      : FilteredMultiVector::dotProduct
// Purpose       : Dot product each column with dense vector.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 2/01/17
//-----------------------------------------------------------------------------
void FilteredMultiVector::dotProduct( const Vector& V, std::vector<double>& dot ) const
{
  int numCols = colPtr_.size()-1;

  // Make sure the resulting vector is the right size.
  dot.resize( numCols );
  std::vector<double> tmpDot( numCols, 0.0 );

  // Perform sparse dot product with dense multivector.
  for (int j=0; j<numCols; j++)
  {
    double tmp=0.0;
    for (int ptr=colPtr_[j]; ptr<colPtr_[j+1]; ptr++)
    {
      tmp += V[rowIndices_[ptr]]*values_[ptr];
    }
    tmpDot[j] = tmp;
  }

  // Perform parallel summation.
  (V.pdsComm())->sumAll( &tmpDot[0], &dot[0], numCols ); 
}

} // namespace Linear
} // namespace Xyce

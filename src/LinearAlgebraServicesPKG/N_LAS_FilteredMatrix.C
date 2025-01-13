//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
#include <N_UTL_FeatureTest.h>

#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Graph.h>
#include <N_LAS_FilteredMatrix.h>
#include <N_LAS_SystemHelpers.h>
#include <N_LAS_Importer.h>
#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>

namespace Xyce {
namespace Linear {

namespace {

struct valComp
{
  using result_type = bool;
  using first_argument_type = std::pair<double,int>;
  using second_argument_type = std::pair<double,int>;

  bool operator()(const std::pair<double,int>& a, const std::pair<double,int>& b) const
  {
    return std::abs(a.first) > std::abs(b.first); 
  }
  bool operator()(const std::pair<double,int>* a, const std::pair<double,int>* b) const
  {
    return std::abs(a->first) > std::abs(b->first); 
  }
};

}

//-----------------------------------------------------------------------------
// Function      : FilteredMatrix::FilteredMatrix
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, Sandia National Labs
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
FilteredMatrix::FilteredMatrix( const Matrix* matrix, const Parallel::ParMap* map, 
                                bool filterOverlap )
: filterOverlap_(filterOverlap),
  totalNZRows_( 0 ),
  minValue_(0.0),
  maxValue_(0.0)
{
  filterMatrix( matrix, map, true );
}

FilteredMatrix::FilteredMatrix( const std::vector<int>& ptr, const std::vector<int>& indices,
                                const std::vector<double>& values, bool isCSR )
: filterOverlap_(false),
  totalNZRows_( 0 ),
  minValue_(0.0),
  maxValue_(0.0)
{
  if (isCSR)
  {
    rowPtr_ = ptr;
    colIndices_ = indices;
    vecIndices_ = indices;
    values_ = values;
  
    // Compute minValue_ and maxValue_
    minValue_ = std::abs(values[0]);
    maxValue_ = std::abs(values[0]);
    for (std::vector<double>::iterator it=values_.begin(); it!=values_.end(); it++)
    {
      double absVal = std::abs(*it);
      if (absVal < minValue_)
        minValue_ = absVal;
      if (absVal > maxValue_)
        maxValue_ = absVal;
    }

    // Compute nonzero columns. 
    nzCols_ = vecIndices_;
    std::sort( nzCols_.begin(), nzCols_.end() );
    nzCols_.erase( std::unique( nzCols_.begin(), nzCols_.end() ), nzCols_.end() );
   
    // Compute nonzero rows.
    for (unsigned int i=1; i<rowPtr_.size(); i++)
    {
      if ((rowPtr_[i]-rowPtr_[i-1]) > 0)
      {
        nzRows_.push_back( i-1 );
      }
    }
    totalNZRows_ = nzRows_.size(); 
  }
  else
  {
    // Convert CSC to CSR.
    // Compute nonzero rows.
    nzRows_ = indices;
    std::sort( nzRows_.begin(), nzRows_.end() );
    nzRows_.erase( std::unique( nzRows_.begin(), nzRows_.end() ), nzRows_.end() );
    totalNZRows_ = nzRows_.size(); 

    rowPtr_.resize( ptr.size(), 0 );

    // Compute nonzero columns.
    for (unsigned int i=1; i<ptr.size(); i++)
    {
      if ((ptr[i]-ptr[i-1]) > 0)
      {
        nzCols_.push_back( i-1 );
      }
    }

    // Loop over each nonzero row and collect nonzero column entries.
    int lastNZRow = nzRows_[ totalNZRows_-1 ];
    for (std::vector<int>::const_iterator row_it = nzRows_.begin(); row_it != nzRows_.end(); row_it++)
    {
      int row = *row_it;
      int newEntries = 0;
      for (std::vector<int>::const_iterator col_it = nzCols_.begin(); col_it != nzCols_.end(); col_it++)
      {
        int col = *col_it;
        for (int k=ptr[ col ]; k<ptr[ col+1 ]; k++)
        {
          if (indices[k] == row)
          {
            colIndices_.push_back( col );
            values_.push_back( values[k] );
          
            double absVal = std::abs(values[k]);
            if (absVal < minValue_ || values_.size() == 1)
              minValue_ = absVal;
            if (absVal > maxValue_)
              maxValue_ = absVal;

            newEntries++;
          }
        }
      }

      // Update pointer for this row and all empty rows following it.
      int newPtr = rowPtr_[ row ] + newEntries;
      if ( row != lastNZRow )
      {
        
        for (int k=row; k<*(row_it+1); k++)
        {
          rowPtr_[ k+1 ] = newPtr;
        }
      }
      else
      {
        rowPtr_[ row+1 ] = newPtr;
      }
    }

    // Set row pointers for the remainder of the empty rows, if there are any.
    for (unsigned int k=lastNZRow+1; k<rowPtr_.size()-1; k++)
    {  
      rowPtr_[ k+1 ] = rowPtr_[ k ];
    }
 
    // Set the vector indices to the column indices.
    vecIndices_ = colIndices_;
  }

}

//-----------------------------------------------------------------------------
// Function      : FilteredMatrix::filterMatrix
// Purpose       : Filter out zero values from input matrix
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, Sandia National Labs
// Creation Date : 3/4/01
//-----------------------------------------------------------------------------
bool FilteredMatrix::filterMatrix( const Matrix* matrix, const Parallel::ParMap* map, bool reset )
{
  bool isReset = reset;
  bool bSuccess = true;

  const Graph* matrixGraph = 0;

  if (filterOverlap_)
    matrixGraph = matrix->getOverlapGraph();
  else
    matrixGraph = matrix->getGraph();

  int maxEntries = matrixGraph->numLocalNonzeros();
  int numMyRows = matrixGraph->numLocalEntities();
  int indexBase = matrixGraph->indexBase();

  // Delete/clear all internal objects and recreate them
  if (reset)
  {
    clearFilter();

    rowPtr_.resize( numMyRows+1 );
    rowPtr_[0] = 0;
  }

  // Get the local ID for the ground node, we don't need to keep these values.
  int row_groundID = -1, col_groundID = -1;
  if (indexBase == -1)
  {
    row_groundID = matrixGraph->globalToLocalRowIndex( -1 );
    col_groundID = matrixGraph->globalToLocalColIndex( -1 );
  }

  int numEntries, currEntries=0;
  int * indices;
  double * values;
  std::vector<std::pair<double,int> > valColPair;

  // Collect information about whether there are off processor columns.
  std::set<int> colIdxs;
  int colsOffProc = 0;

  // If there isn't a reset, use as much information about the previous filtered matrix as possible.
  // If this is a reset and the previous matrix was empty, just make sure that this one is empty too.
  if (!reset && nzRows_.size())
  {
    numMyRows = nzRows_.size();
  }

  for( int i = 0; i < numMyRows; ++i )
  {
    valColPair.clear();

    int rowIdx = i;
    if (!reset && nzRows_.size())
      rowIdx = nzRows_[i];

    if (rowIdx != row_groundID)
    {
      if (filterOverlap_)
      {
        matrix->extractLocalRowView( rowIdx, numEntries, values, indices );
      }
      else
      {
        matrix->getLocalRowView( rowIdx, numEntries, values, indices );
      }
      for( int j = 0; j < numEntries; ++j )
      {
        // Filter nonzero entries that aren't the ground node.
        if( (values[j] != 0.0) && (indices[j] != col_groundID) )
        {
          valColPair.push_back( std::pair<double,int>(values[j], indices[j]) );

          double absVal = std::abs(values[j]);
          if (absVal < minValue_ || values_.size() == 0)
            minValue_ = absVal;
          if (absVal > maxValue_)
            maxValue_ = absVal;
        }
      }
    }

    // Now update pointers for this row, copy over sorted entries. 
    int new_nz = valColPair.size();
   
    // Check that this is the same number of nonzeros as before, if this is not a reset.
    if (!reset && (new_nz != (rowPtr_[rowIdx+1]-rowPtr_[rowIdx])))
    {
      bSuccess = false;
      break;
    }
 
    if (bSuccess && new_nz)
    {
      // Sort the entries largest to smallest.
      // NOTE:  In general smallest to largest is the best practice, but large
      //        entries are, most likely, the result of small resistors which will generate
      //        cancellation.  So, this will let the cancellation happen before any additional
      //        computations.
      if (new_nz > 1)
      {
        std::sort( valColPair.begin(), valColPair.end(), valComp() ); 
      }

      for (int j = 0; j < new_nz; ++j )
      {
        if (reset)
        {
          values_.push_back( valColPair[j].first );
          colIndices_.push_back( valColPair[j].second );
          if (filterOverlap_)
          {
            vecIndices_.push_back( map->globalToLocalIndex(matrixGraph->localToGlobalColIndex(valColPair[j].second)) );
          }
          else
          {
            int idx = matrixGraph->localToGlobalColIndex(valColPair[j].second);
            vecIndices_.push_back( idx );
            colIdxs.insert( idx );
            if ( map->globalToLocalIndex( idx ) < 0 )
              colsOffProc++;
          }
        }
        else
        {
          values_[currEntries] = valColPair[j].first;
          colIndices_[currEntries] = valColPair[j].second;
          if (filterOverlap_)
          {
            vecIndices_[currEntries] = map->globalToLocalIndex(matrixGraph->localToGlobalColIndex(valColPair[j].second));
          }
          else
          {
            // Reuse the current targetMap_, if it exists.
            if (!Teuchos::is_null(targetMap_))
            {
              vecIndices_[currEntries] = targetMap_->globalToLocalIndex(matrixGraph->localToGlobalColIndex(valColPair[j].second));
            }
            else
            {
              vecIndices_[currEntries] = map->globalToLocalIndex(matrixGraph->localToGlobalColIndex(valColPair[j].second));
            } 
            if (vecIndices_[currEntries] < 0)
            {
              bSuccess = false;
            }
          }
          currEntries++;
        }
      }

      if (reset)
      {
        rowPtr_[i+1] = rowPtr_[i] + new_nz;
        nzRows_.push_back( i );
      }
    }
    else
    {
      if (reset)
        rowPtr_[i+1] = rowPtr_[i];
    }
  }

  if (reset)
  {
    // At this point, all of the sanity checks have passed, so create a new targetMap_/importer_ if necessary.
    // Get the total number of nonzero rows in the global filtered matrix.
    int localNZRows = nzRows_.size();
    (map->pdsComm()).sumAll( &localNZRows, &totalNZRows_, 1 );

    if (!filterOverlap_)
    {
      // Determine if we need to create an exporter to perform matvec/axpy.
      int globalColsOffProc = 0;
      (map->pdsComm()).sumAll( &colsOffProc, &globalColsOffProc, 1 );

      if (globalColsOffProc)
      {
        const Parallel::ParMap* colMap = matrix->getColMap( map->pdsComm() );
        importer_ = Teuchos::rcp( Xyce::Linear::createImporter( *colMap, *map ) ); 
        targetMap_ = Teuchos::rcp( colMap->clone() );

        // Compute the local indicies for the target map.
        for (unsigned int i=0; i<vecIndices_.size(); i++)
        {
          vecIndices_[i] = targetMap_->globalToLocalIndex( vecIndices_[i] );
        }
      }
      else
      {
        // There is no need to import information to perform a matvec/axpy, so use local IDs.
        for (unsigned int i=0; i<vecIndices_.size(); i++)
        {
          vecIndices_[i] = map->globalToLocalIndex( vecIndices_[i] );
        }
      }
    }

    // Compute the nonzero columns from the indices.
    nzCols_ = vecIndices_;
    std::sort( nzCols_.begin(), nzCols_.end() );
    nzCols_.erase( std::unique( nzCols_.begin(), nzCols_.end() ), nzCols_.end() );
  }
  else
  {
    // Check bSuccess across all processors.  If unsuccessful, call filterMatrix with reset=true.
    int global_bSuccess = 0;
    int local_bSuccess = bSuccess ? 1 : 0;
    (map->pdsComm()).minAll( &local_bSuccess, &global_bSuccess, 1 );
    if (global_bSuccess == 0)
    {
      filterMatrix( matrix, map, true );
      isReset = true;
    }
  }

  return isReset;
}

//-----------------------------------------------------------------------------
// Function      : FilteredMatrix::clearFilter
// Purpose       : Clear filtered matrix objects.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, Sandia National Labs
// Creation Date : 3/4/01
//-----------------------------------------------------------------------------
void FilteredMatrix::clearFilter()
{
  nzRows_.clear();
  nzCols_.clear();
  colIndices_.clear();
  vecIndices_.clear();
  values_.clear();
}

//-----------------------------------------------------------------------------
// Function      : FilteredMatrix::axpy
// Purpose       : Perform y = A*x + y using the filtered matrix 
// Special Notes : Using assembled filtered matrix here.
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/11/04
//-----------------------------------------------------------------------------
void FilteredMatrix::axpy(const MultiVector & x, MultiVector & y) 
{ 
  if (!isEmpty())
  {
    // Untransposed filtered matrix axpy 
    int numCols_x = x.numVectors();
    int numCols_y = y.numVectors();
    int numRows_x = x.localLength();
    int numRows_y = y.localLength();
    if (numCols_x != numCols_y)
    {
      Report::DevelFatal()
          << "FilteredMatrix matrix vector product, x and y have different numbers of columns.";
    }
    if (numRows_y != numRows_x)
    {
      Report::DevelFatal()
          << "FilteredMatrix matrix vector product, x and y have different numbers of rows.";
    }

    if (!Teuchos::is_null(importer_) && !Teuchos::is_null(targetMap_))
    {
      if (Teuchos::is_null(targetX_))
      {
        targetX_ = Teuchos::rcp( Linear::createMultiVector( *targetMap_, x.numVectors() ) );
      }
      else
      {
        if (targetX_->numVectors() != x.numVectors())
        {
          targetX_ = Teuchos::rcp( Linear::createMultiVector( *targetMap_, x.numVectors() ) );
        }
      }
      targetX_->vectorImport( &x, &*importer_ );
    } 

    for (int j=0; j<numCols_x; j++)
    {
      for( std::vector<int>::iterator it = nzRows_.begin(); it != nzRows_.end(); it++ )
      {
        int row = *it;
        int rowPtr = rowPtr_[row];
        int rowPtr_next = rowPtr_[row+1];
        double sum=0.0;
        //double ysave = y[j][row];
        for (int j2=rowPtr; j2<rowPtr_next; j2++)
        {
          if (filterOverlap_)
          {
            sum += values_[j2]*(*x(vecIndices_[j2],j));
          }
          else
          {
            if (!Teuchos::is_null(targetX_))
            {
              sum += values_[j2]*(*(*targetX_)(vecIndices_[j2],j));
            }
            else
            {
              sum += values_[j2]*(*x(vecIndices_[j2],j));
            }
          }
        }
        (*y(row,j)) += sum;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : FilteredMatrix::matvec
// Purpose       : Perform y = A*x using the filtered matrix
// Special Notes : Using assembled filtered matrix here.
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/11/04
//-----------------------------------------------------------------------------
void FilteredMatrix::matvec(const MultiVector & x, MultiVector & y)
{
  if (!isEmpty())
  {
    // Untransposed filtered matrix axpy
    int numCols_x = x.numVectors();
    int numCols_y = y.numVectors();
    int numRows_x = x.localLength();
    int numRows_y = y.localLength();
    if (numCols_x != numCols_y)
    {
      Report::DevelFatal()
          << "FilteredMatrix matrix vector product, x and y have different numbers of columns.";
    }
    if (numRows_y != numRows_x)
    {
      Report::DevelFatal()
          << "FilteredMatrix matrix vector product, x and y have different numbers of rows.";
    }

    if (!Teuchos::is_null(importer_) && !Teuchos::is_null(targetMap_))
    {
      if (Teuchos::is_null(targetX_))
      {
        targetX_ = Teuchos::rcp( Linear::createMultiVector( *targetMap_, x.numVectors() ) );
      }
      else
      {
        if (targetX_->numVectors() != x.numVectors())
        {
          targetX_ = Teuchos::rcp( Linear::createMultiVector( *targetMap_, x.numVectors() ) );
        }
      }
      targetX_->vectorImport( &x, &*importer_ );
    } 

    for (int j=0; j<numCols_x; j++)
    {
      for( std::vector<int>::iterator it = nzRows_.begin(); it != nzRows_.end(); it++ )
      {
        int row = *it;
        int rowPtr = rowPtr_[row];
        int rowPtr_next = rowPtr_[row+1];
        double sum=0.0;
        //double ysave = y[j][row];
        for (int j2=rowPtr; j2<rowPtr_next; j2++)
        {
          if (filterOverlap_)
          {
            sum += values_[j2]*(*x(vecIndices_[j2],j));
          }
          else
          {
            if (!Teuchos::is_null( targetX_ ))
            {
              sum += values_[j2]*(*(*targetX_)(vecIndices_[j2],j));
 /*
              std::cout << "values_[" << j2 << "] = " << values_[j2] << ", vecIndices_[" << j2 << "] = "
                        << vecIndices_[j2] << ", targetX_[" << j << "][" << vecIndices_[j2] << " = "
                        << *(*targetX_)(vecIndices_[j2],j) << std::endl;
 */
            } 
            else
            {
              sum += values_[j2]*(*x(vecIndices_[j2],j));
            }
          }
        }
        (*y(row,j)) = sum;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : FilteredMatrix::addToMatrix
// Purpose       : Sums in a matrix contribution
// Special Notes : WARNING: only works if graphs match, no checking
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/11/04
//-----------------------------------------------------------------------------
void FilteredMatrix::addToMatrix( Matrix & A, double alpha )
{
  if (!isEmpty())
  {
    double* scaledValues = values_.data();
    if (alpha != 1.0)
    {
      int i=0;
      scaledValues = new double[values_.size()];
      for ( std::vector<double>::iterator it = values_.begin(); it != values_.end(); it++, i++ )
      {
        scaledValues[i] = alpha*(*it);
      }
    }

    for( std::vector<int>::iterator it = nzRows_.begin(); it != nzRows_.end(); it++ )
    {
      int row = *it;
      int rowPtr = rowPtr_[row];
      if (filterOverlap_)
      {
        A.sumIntoLocalRow( row, rowPtr_[row+1]-rowPtr, 
                           scaledValues+rowPtr, &colIndices_[rowPtr] );
      }
      else
      {
        A.addIntoLocalRow( row, rowPtr_[row+1]-rowPtr, 
                           scaledValues+rowPtr, &colIndices_[rowPtr] );
      }
    }

    if (alpha != 1.0)
      delete [] scaledValues;
  }
}

//-----------------------------------------------------------------------------
// Function      : FilteredMatrix::printFilteredMatrix
// Purpose       : 
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/11/04
//-----------------------------------------------------------------------------
void FilteredMatrix::printFilteredMatrix(std::ostream& os)
{
  os << "Filtered Matrix:  nzRows = " << nzRows_.size() << ", total_nzRows = " 
     << totalNZRows_ << ", minValue = " << minValue_
     << ", maxValue = " << maxValue_ << ")" << std::endl;
  for( std::vector<int>::iterator it = nzRows_.begin(); it != nzRows_.end(); it++ )
  {
    int row = *it;
    int rowPtr = rowPtr_[row];
    int rowPtr_next = rowPtr_[row+1];
    if (filterOverlap_)
      os << "Row " << row << " : ";
    else
      os << "Row " << row << " : ";
    for (int j2=rowPtr; j2<rowPtr_next; j2++)
    { 
      if (filterOverlap_)
        os << values_[j2] << " [" << colIndices_[j2] << "], ";
      else
        os << values_[j2] << " [" << colIndices_[j2] << "], ";
    }
    os << std::endl;
  }
}

} // namespace Linear
} // namespace Xyce

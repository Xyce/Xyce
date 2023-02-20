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
// Purpose        : HB Block Matrix
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Parallel Computational Sciences
//
// Creation Date  : 12/7/16
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_BlockMatrixEntry_h
#define Xyce_N_LAS_BlockMatrixEntry_h

#include <Xyce_config.h>
#include <N_ERH_ErrorMgr.h>

#include <vector>
#include <iostream>
#include <numeric>

#ifdef Xyce_AMESOS2_BASKER

#include "Amesos2_config.h"
#include "Amesos2_Basker_TypeMap.hpp"
#include "Amesos2_Basker.hpp"

#endif

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"

namespace Xyce 
{
// Create a block matrix for each Xyce matrix entry.
// This matrix can be dense (SerialDenseMatrix) or diagonal (std::vector).
template<class T>
struct genericBlockMatrixEntry 
{
  typedef T val_type;
  typedef typename Teuchos::ScalarTraits<val_type>::magnitudeType mag_type;
  typedef Teuchos::SerialDenseMatrix<int, val_type> mtx_type_nonlin;
  typedef Teuchos::SerialDenseSolver<int, val_type> nonlin_solver;
  typedef std::vector<val_type> mtx_type_lin;

  // Default constructor
  genericBlockMatrixEntry()
  : rows( 0 ),
    cols( 0 )
  {}

  // Null constructor
  genericBlockMatrixEntry( int should_be_null )
  : rows( 0 ),
    cols( 0 )
  {
    if ( should_be_null )
      std::cout << "Null constructor received value: " << should_be_null << std::endl;
  }

  // Basic constructor
  genericBlockMatrixEntry( int numRows, int numCols, bool isDense )
  : rows( numRows ),
    cols( numCols )
  {
    if (isDense)
      denseMtx.reshape( rows, cols );
    else
      diagVector.resize( rows );
  }

  // Copy contructor
  genericBlockMatrixEntry( const genericBlockMatrixEntry& Source )
  : rows( Source.rows ),
    cols( Source.cols )
  {
    if (Source.isDense())
    {
      denseMtx.reshape( rows, cols );
      denseMtx.assign( Source.denseMtx );
    }
    else
    {
      diagVector = Source.diagVector;
    }
  }

  // Is this entry related to a linear block (diagonal)
  bool isDiag() const { return (diagVector.size() > 0); }

  // Is this entry related to a nonlinear block (dense)
  bool isDense() const { return ( !denseMtx.empty() ); }   

  // Overload operators for this matrix.
  template<class R>genericBlockMatrixEntry& operator=  (const R value); //R = return type
  genericBlockMatrixEntry& operator=  (const genericBlockMatrixEntry& Source);

  const genericBlockMatrixEntry operator+  (const genericBlockMatrixEntry& Source) const;
  const genericBlockMatrixEntry operator-  (const genericBlockMatrixEntry& Source) const;
  const genericBlockMatrixEntry operator*  (const genericBlockMatrixEntry& Source) const;
  const genericBlockMatrixEntry operator/  (const genericBlockMatrixEntry& Source) const;

  genericBlockMatrixEntry& operator+= (const genericBlockMatrixEntry& Source);
  genericBlockMatrixEntry& operator-= (const genericBlockMatrixEntry& Source);
  genericBlockMatrixEntry& operator*= (const genericBlockMatrixEntry& Source);
  genericBlockMatrixEntry& operator/= (const genericBlockMatrixEntry& Source);

  bool operator==(const genericBlockMatrixEntry& Source) const;
  bool operator!=(const genericBlockMatrixEntry& Source) const;
  bool operator> (const genericBlockMatrixEntry& Source) const;
  bool operator< (const genericBlockMatrixEntry& Source) const;

  // Initialize values of either diagonal or dense block with value.
  void putScalar(const val_type value);

  // This function will be used to expand any 1x1s to be compatible with block sizes
  void expandDiag( int nrows );

  // Add a value to the diagonal entry of either dense or diagonal matrix.
  void addToDiag( int index, val_type val ); 

  // Frobenius norm
  double normFrobenius() const;

  void print(std::ostream& os) const;
 
  int rows, cols;
  mtx_type_nonlin denseMtx;
  mtx_type_lin diagVector;
};

template<class T>
template<class R>
inline genericBlockMatrixEntry<T>& genericBlockMatrixEntry<T>::operator= (const R value)
{
  if ( isDense() )
  {
    denseMtx.putScalar( value );
  }
  else
  {
    // As a workaround for Basker, if this is an empty matrix, create a 1x1.
    // The inconsistent sizing will be resolved in a later operation if it is used.
    if (rows == 0)
    {
      rows = 1;
      cols = 1;
      diagVector.resize( rows );
    }

    for (int i=0; i<rows; i++)
    {
      diagVector[i] = value;
    }
  }

  return *this;
}

// ERK. specialization for assigning double to a std::complex.  This might not be necessary.
template<>
template<>
inline genericBlockMatrixEntry< std::complex<double> >& genericBlockMatrixEntry< std::complex<double>  >::operator= (const double value)
{
  if ( isDense() )
  {
    denseMtx.putScalar( std::complex<double>( value, 0.0 ) );
  }
  else
  {
    // As a workaround for Basker, if this is an empty matrix, create a 1x1.
    // The inconsistent sizing will be resolved in a later operation if it is used.
    if (rows == 0)
    {
      rows = 1;
      cols = 1;
      diagVector.resize( rows );
    }

    for (int i=0; i<rows; i++)
    {
      diagVector[i] = std::complex<double>( value, 0.0 );
    }
  }

  return *this;
}

// Assign one block matrix to another.
template<class T>
inline genericBlockMatrixEntry<T>& genericBlockMatrixEntry<T>::operator=  (const genericBlockMatrixEntry& Source)
{
  rows = Source.rows;
  cols = Source.cols;
  if ( Source.isDense() )
  {
    denseMtx.reshape( Source.rows, Source.cols );
    denseMtx.assign( Source.denseMtx );
    diagVector.clear();
  }
  else
  {
    diagVector = Source.diagVector;
    denseMtx.reshape( 0, 0 );
  }

  return *this;
}

// Return this + b
template<class T>
inline const genericBlockMatrixEntry<T> genericBlockMatrixEntry<T>::operator+  (const genericBlockMatrixEntry& Source) const
{
  genericBlockMatrixEntry result( *this );

  if ( rows != Source.rows )
  {
    if ( (rows == 1) && isDiag() )
    {
      result.expandDiag( Source.rows );
      result += Source;
    }

    if ( (Source.rows == 1) && Source.isDiag() )
    {
      genericBlockMatrixEntry newSource( Source );
      newSource.expandDiag( rows );
      result += newSource;
    }
  }
  else
  {
    result += Source;
  }

  return result;
}

// Return this - b
template<class T>
inline const genericBlockMatrixEntry<T> genericBlockMatrixEntry<T>::operator-  (const genericBlockMatrixEntry& Source) const
{
  genericBlockMatrixEntry result( *this );

  if ( rows != Source.rows )
  {
    if ( (rows == 1) && isDiag() )
    {
      result.expandDiag( Source.rows );
      result -= Source;
    }
  
    if ( (Source.rows == 1) && Source.isDiag() )
    {
      genericBlockMatrixEntry newSource( Source );
      newSource.expandDiag( rows );
      result -= newSource;
    }
  }
  else
  {
    result -= Source;
  
  } 

  return result;
}

// Return this * b
template<class T>
inline const genericBlockMatrixEntry<T> genericBlockMatrixEntry<T>::operator*  (const genericBlockMatrixEntry& Source) const
{
  genericBlockMatrixEntry result;
  bool isMatrix = ((Source.rows == Source.cols) && (rows == cols));
 
  if ( isMatrix )
  {
    result = *this;
  
    if ( rows != Source.rows )
    {
      if ( (rows == 1) && isDiag() )
      {
        if ( diagVector[0] != Teuchos::ScalarTraits<val_type>::one() )
        {
          result.expandDiag( Source.rows );
          result *= Source;
        }
        else
          result = Source;
      }

      if ( (Source.rows == 1) && Source.isDiag() )
      {
        if ( Source.diagVector[0] != Teuchos::ScalarTraits<val_type>::one() )
        {
          genericBlockMatrixEntry newSource( Source );
          newSource.expandDiag( rows );
          result *= newSource;
        }
      }
    }
    else
    {
      result *= Source;
    }
  }
  else
  {
    // This is a mat-vec.
    result.rows = Source.rows;
    result.cols = Source.cols;
    result.denseMtx.reshape( result.rows, result.cols );

    if ( isDense() )
    {
      // Create a dense matrix to collect the result of the GEMM.
      result.denseMtx.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, denseMtx, Source.denseMtx, 0.0 );
    }
    else
    {
      // Copy over dense matrix, scale rows.
      result.denseMtx.assign( Source.denseMtx );

      val_type alpha = diagVector[0];
      for (int i=0; i<result.rows; i++)
      {
        if (rows > 1)
          alpha = diagVector[i];
      
        if ( alpha != Teuchos::ScalarTraits<val_type>::one() )
        {
          for (int j=0; j<result.cols; j++)
          {
            result.denseMtx(i,j) *= alpha;
          }
        }
      }
    }
  }

  return result;
}

// Return this / b
template<class T>
inline const genericBlockMatrixEntry<T> genericBlockMatrixEntry<T>::operator/  (const genericBlockMatrixEntry& Source) const
{
  genericBlockMatrixEntry result;
  bool isMatrix = ((Source.rows == Source.cols) && (rows == cols));

  if ( isMatrix )
  {
    result = *this;
  
    if ( rows != Source.rows )
    {
      if ( (rows == 1) && isDiag() )
      {
        result.expandDiag( Source.rows );
        result /= Source;
      }

      if ( (Source.rows == 1) && Source.isDiag() )
      {
        if ( Source.diagVector[0] != Teuchos::ScalarTraits<val_type>::one() )
        {
          genericBlockMatrixEntry newSource( Source );
          newSource.expandDiag( rows );
          result /= newSource;
        }
      }
    }
    else
    {
      result /= Source;
    }
  }
  else
  {
    // This is a back solve with a vector.
    result.rows = rows;
    result.cols = cols;
    result.denseMtx.reshape( result.rows, result.cols );

    if ( Source.isDense() )
    {
      // Compute inverse of source matrix.
      mtx_type_nonlin A( Teuchos::Copy, Source.denseMtx );
      mtx_type_nonlin b( Teuchos::Copy, denseMtx );
      nonlin_solver denseSolver;
      denseSolver.factorWithEquilibration( true );
      denseSolver.setMatrix( Teuchos::rcp( &A, false ) );
      denseSolver.setVectors( Teuchos::rcp( &(result.denseMtx), false ), Teuchos::rcp( &b, false ) );
      denseSolver.factor();
      denseSolver.solve();
    }
    else
    {
      // This is a set of vectors scaled by the diagonal of Source.
      result.denseMtx.assign( denseMtx );

      // Now perform row scaling of the inverted matrix.
      val_type alpha = Source.diagVector[0];
      for (int i=0; i<result.rows; i++)
      {
        if (Source.rows > 1)
          alpha = Source.diagVector[i];

        if ( alpha != Teuchos::ScalarTraits<val_type>::one() )
        {
          for (int j=0; j<result.cols; j++)
          {
            result.denseMtx(i,j) /= alpha;
          }
        }
      }
    }
  }

  return result;
}

// Add one block matrix to another.
template<class T>
inline genericBlockMatrixEntry<T>& genericBlockMatrixEntry<T>::operator+= (const genericBlockMatrixEntry& Source)
{
  if ( (Source.rows != rows) || Source.cols != cols )
    Report::DevelFatal0() << "genericBlockMatrixEntry::operator+= : matrices are not compatible!";

  // Check all four possible scenarios of addition.
  if ( Source.isDense() )
  {
    if ( isDense() )
    {
      // Add into current dense matrix.
      denseMtx += Source.denseMtx;
    }
    else
    {
      // Copy over dense matrix and then add diagonals.
      denseMtx.reshape( Source.rows, Source.cols );
      denseMtx.assign( Source.denseMtx );
      for (int i=0; i<Source.rows; i++)
      {
        denseMtx(i,i) += diagVector[i];
      }
      diagVector.clear();
    }
  }
  else
  {
    // Add diagonals into current matrix.
    for (int i=0; i<rows; i++)
    {
      if ( isDense() )
        denseMtx(i,i) += Source.diagVector[i];
      else
        diagVector[i] += Source.diagVector[i];
    }
  }
 
  return *this;
}

template<class T>
inline genericBlockMatrixEntry<T>& genericBlockMatrixEntry<T>::operator-= (const genericBlockMatrixEntry& Source)
{
  if ( (Source.rows != rows) || Source.cols != cols )
    Report::DevelFatal0() << "genericBlockMatrixEntry::operator-= : matrices are not compatible!";

  // Check all four possible scenarios of subtraction.
  if ( Source.isDense() )
  {
    if ( isDense() )
    {
      // Add into current dense matrix.
      denseMtx -= Source.denseMtx;
    }
    else
    {
      // Copy over dense matrix and then add diagonals.
      denseMtx.reshape( Source.rows, Source.cols );
      denseMtx.putScalar( Teuchos::ScalarTraits<val_type>::zero() );
      for (int i=0; i<Source.rows; i++)
      {
        denseMtx(i,i) = diagVector[i];
      }
      denseMtx -= Source.denseMtx;
      diagVector.clear();
    }
  }
  else
  {
    // Add diagonals into current matrix.
    for (int i=0; i<rows; i++)
    {
      if ( isDense() )
        denseMtx(i,i) -= Source.diagVector[i];
      else
        diagVector[i] -= Source.diagVector[i];
    }
  }

  return *this;
}

template<class T>
inline genericBlockMatrixEntry<T>& genericBlockMatrixEntry<T>::operator*= (const genericBlockMatrixEntry& Source)
{
  if ( (Source.rows != rows) || Source.cols != cols )
    Report::DevelFatal0() << "genericBlockMatrixEntry::operator*= : matrices are not compatible!";

  // Check all four possible scenarios of multiplication.
  if ( Source.isDense() )
  {
    if ( isDense() )  
    {
      // Create a dense matrix to collect the result of the GEMM.
      mtx_type_nonlin result( rows, cols );
      result.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, denseMtx, Source.denseMtx, 0.0 );
      denseMtx.assign( result );
    }
    else
    {
      // Copy over dense matrix, scale rows.
      denseMtx.reshape( Source.rows, Source.cols );
      denseMtx.assign( Source.denseMtx ); 
      for (int i=0; i<rows; i++)
      {
        val_type alpha = diagVector[i];
        if ( alpha != Teuchos::ScalarTraits<val_type>::one() )
        {
          mtx_type_nonlin row_i( Teuchos::View, denseMtx, 1, cols, i, 0 );
          row_i.scale( alpha );
        }
      }
      diagVector.clear();
    }
  }
  else
  {
    if ( isDense() )
    {
      // Scale cols.
      for (int i=0; i<cols; i++)
      { 
        val_type alpha = Source.diagVector[i];
        if ( alpha != Teuchos::ScalarTraits<val_type>::one() )
        {
          mtx_type_nonlin col_i( Teuchos::View, denseMtx, rows, 1, 0, i );
	  col_i.scale( alpha );
        }
      }
    }
    else
    {
      // Multiply diagonals together.
      for (int i=0; i<rows; i++)
      {
        diagVector[i] *= Source.diagVector[i];
      }
    }
  }
 
  return *this;
}

template<class T>
inline genericBlockMatrixEntry<T>& genericBlockMatrixEntry<T>::operator/= (const genericBlockMatrixEntry& Source)
{
  if ( (Source.rows != rows) || Source.cols != cols )
    Report::DevelFatal0() << "genericBlockMatrixEntry::operator/= : matrices are not compatible!";

  // Check all four possible scenarios of division.
  if ( Source.isDense() )
  {
    if ( isDense() )
    {
      // Compute inverse of source matrix.
      // this = this / Source = this * inv(Source) = ( Source' \ this' )'
      mtx_type_nonlin A( Source.denseMtx, Teuchos::TRANS );
      mtx_type_nonlin b( denseMtx, Teuchos::TRANS );
      nonlin_solver denseSolver;
      denseSolver.factorWithEquilibration( true ); 
      denseSolver.setMatrix( Teuchos::rcp( &A, false ) );
      denseSolver.setVectors( Teuchos::rcp( &denseMtx, false ), Teuchos::rcp( &b, false ) );
      denseSolver.factor();
      denseSolver.solve();

      // Now transpose the solution back.
      for (int i=0; i<rows; i++)
      {
        for (int j=0; j<i; j++)
        {
          val_type tmp = denseMtx(i,j);
          denseMtx(i,j) = denseMtx(j,i);
          denseMtx(j,i) = tmp;
        }
      }
    } 
    else
    {
      // This is now a dense matrix that is the row scaling of the inv(Source).
      denseMtx.reshape( Source.rows, Source.cols );
      denseMtx.assign( Source.denseMtx );
      nonlin_solver denseSolver;
      denseSolver.setMatrix( Teuchos::rcp( &denseMtx, false ) );
      denseSolver.factor();
      denseSolver.invert();

      // Now perform row scaling of the inverted matrix.
      for (int i=0; i<rows; i++)
      {
        val_type alpha = diagVector[i];
        if ( alpha != Teuchos::ScalarTraits<val_type>::one() )
        {
          mtx_type_nonlin row_i( Teuchos::View, denseMtx, 1, cols, i, 0 );
          row_i.scale( alpha );
        }
      }
      diagVector.clear();
    }
  }
  else
  {
    if ( isDense() )
    {
      // Scale cols.
      for (int i=0; i<cols; i++)
      {
        val_type alpha = Source.diagVector[i];
        if ( alpha != Teuchos::ScalarTraits<val_type>::one() )
        {
          mtx_type_nonlin col_i( Teuchos::View, denseMtx, rows, 1, 0, i );
          col_i.scale( Teuchos::ScalarTraits<val_type>::one() / alpha );
        }
      }
    }
    else
    {
      // Multiply diagonals together.
      for (int i=0; i<rows; i++)
      {
        diagVector[i] /= Source.diagVector[i];
      }
    }
  }  

  return *this;
}

template<class T>
inline bool genericBlockMatrixEntry<T>::operator==(const genericBlockMatrixEntry& Source) const
{
  // Check for special case of null matrix.
  bool ret = false;

  if (Source.rows==0 && Source.cols==0 && Source.denseMtx.empty() && Source.diagVector.empty())
  {
    if (rows && cols)
    {
      // This doesn't make a lot of sense, but the abstractions in Basker require that
      // this matrix can be the same as an empty one if it is zero.
      val_type suma = std::accumulate( diagVector.begin(), diagVector.end(), Teuchos::ScalarTraits<val_type>::zero() );
      mag_type sumb = denseMtx.normFrobenius();
      if ( (suma == Teuchos::ScalarTraits<val_type>::zero()) && (sumb == 0.0) )
        ret = true; 
    }
    else
    {
      // This matrix is also empty, so they are the same.
      ret = true;
    }
  }
  else
  {
    if ( isDense() )
    {
      ret = ( denseMtx == Source.denseMtx );
    }
    else
    {
      ret = ( diagVector == Source.diagVector );
    }
  }
     
  return ret; 
}

template<class T>
inline bool genericBlockMatrixEntry<T>::operator!=(const genericBlockMatrixEntry& Source) const
{
  return( !(*this == Source) );
}

template<class T>
inline bool genericBlockMatrixEntry<T>::operator> (const genericBlockMatrixEntry& Source) const
{
  bool ret = false;

  if ( isDiag() )
  {
    // Need to check if this diagonal is singular.
    if ( Source.isDense() )
    {
      bool singA = false;
      for (int i=0; i<rows; i++)
      {
        if ( Teuchos::ScalarTraits<val_type>::magnitude( diagVector[i] ) == 0.0 )
        {
          singA = true;
          break;
        }
      } 
      if (!singA)
        ret = true;
    }
    else
    {
      mag_type mina = Teuchos::ScalarTraits<val_type>::magnitude( diagVector[0] );
      mag_type maxa = Teuchos::ScalarTraits<val_type>::magnitude( diagVector[0] );
      for (int i=1; i<rows; i++)
      {
        mag_type val = Teuchos::ScalarTraits<val_type>::magnitude( diagVector[i] );
        if (val > maxa)
          maxa = val;
        if (val < mina)
          mina = val;
      }

      mag_type minb = Teuchos::ScalarTraits<val_type>::magnitude( Source.diagVector[0] );
      mag_type maxb = Teuchos::ScalarTraits<val_type>::magnitude( Source.diagVector[0] );
      for (int i=1; i<Source.rows; i++)
      {
        mag_type val = Teuchos::ScalarTraits<val_type>::magnitude( Source.diagVector[i] );
        if (val > maxb)
          maxb = val;
        if (val < minb)
          minb = val;
      }

      if (mina == 0.0)
        ret = false;
      else if (minb == 0.0) 
        ret = true;
      else
        ret = ( maxa/mina ) < ( maxb/minb );
    }
  }
  else
  {
    if ( Source.isDense() )
    {
      mag_type condEst, condEst2;
      mtx_type_nonlin inverse( Teuchos::Copy, denseMtx ), inverse2( Teuchos::Copy, Source.denseMtx );
      nonlin_solver denseSolver;

      denseSolver.setMatrix( Teuchos::rcp( &inverse, false ) );
      denseSolver.factor();
      denseSolver.reciprocalConditionEstimate( condEst );

      denseSolver.setMatrix( Teuchos::rcp( &inverse2, false ) );
      denseSolver.factor();
      denseSolver.reciprocalConditionEstimate( condEst2 );

      if ( condEst == 0.0 )
        ret = false;
      else if ( condEst2 == 0.0 )
        ret = true;
      else
        ret = ( 1.0/condEst ) < ( 1.0/condEst2 );
    }
    else
    {
      // Need to check if the Source a zero diagonal.
      bool singB = false;
      for (int i=0; i<rows; i++)
      {
        if ( Teuchos::ScalarTraits<val_type>::magnitude( Source.diagVector[i] ) == 0.0 )
        {
          singB = true;
          break;
        }
      } 
      if (singB)
        ret = true;
    }
  }

  return ret; 
}

template<class T>
inline bool genericBlockMatrixEntry<T>::operator< (const genericBlockMatrixEntry& Source) const
{
  return ( !(*this > Source) );
}

template<class T>
inline void genericBlockMatrixEntry<T>::putScalar(const val_type value)
{
  if ( isDiag() )
  {
    for (unsigned int i=0; i<diagVector.size(); i++)
    {
      diagVector[i] = value;
    }
  }
  if ( isDense() )
  {
    denseMtx.putScalar( value );
  }
}

template<class T>
inline void genericBlockMatrixEntry<T>::print(std::ostream& os) const
{
  if ( isDense() )
  {
    os << "genericBlockMatrixEntry Dense: " << std::endl;
    denseMtx.print( os );
  }
  else
  {
    os << "genericBlockMatrixEntry Diagonal: " << std::endl;
    os << "Rows : " << rows << std::endl;
    os << "Columns : " << cols << std::endl;
    os << "Values : ";
    for (unsigned int i=0; i<diagVector.size(); i++)
    {
      os << diagVector[i] << " ";
    }
    os << std::endl;
  }
}

template<class T>
inline void genericBlockMatrixEntry<T>::expandDiag( int nrows )
{
  val_type value = diagVector[0];

  rows = nrows;
  cols = nrows;
  diagVector.resize( nrows, value );
}

template<class T>
inline void genericBlockMatrixEntry<T>::addToDiag( int index, val_type val )
{
  if (rows && cols)
  {
    if (isDiag())
    {
      diagVector[index] += val;
    }
    else
    {
      denseMtx(index, index) += val;
    }
  }
}

template<class T>
inline double genericBlockMatrixEntry<T>::normFrobenius() const
{
  mag_type norm = 0.0;
  if (isDiag())
  {
    for (int i=0; i<rows; i++)
    {
      norm += Teuchos::ScalarTraits<val_type>::magnitude( diagVector[i]*diagVector[i] );
    }
    norm = Teuchos::ScalarTraits<val_type>::magnitude(
             Teuchos::ScalarTraits<val_type>::squareroot(norm));
  }
  else
  {
    norm = denseMtx.normFrobenius();
  }

  return norm;
}

typedef genericBlockMatrixEntry< std::complex<double> > HBBlockMatrixEntry;
typedef genericBlockMatrixEntry< double > ESBlockMatrixEntry;
typedef ESBlockMatrixEntry PCEBlockMatrixEntry;

template<class T>
inline int packGenericBlockMatrix(const genericBlockMatrixEntry<T>& Source, std::vector<double>& vec)
{
  int sizeVec = vec.size();
  int reqStorage = 2 * Source.rows;
  if ( Source.isDense() )
  {
    reqStorage *= Source.cols;
  }

  // Resize if necessary.
  if ( sizeVec < reqStorage )
  {
    vec.resize( reqStorage );
  } 

  // Now copy values from Source into vec.
  if ( Source.isDense() )
  {
    for (int i=0; i<Source.rows; i++)
    {
      for (int j=0; j<Source.cols; j++)
      {
        vec[ j*Source.rows + i ] = Source.denseMtx( i, j ).real();
        vec[ Source.rows*Source.cols + j*Source.rows + i ] = Source.denseMtx( i, j ).imag();
      }
    }
  }
  else
  {
    for (int i=0; i<Source.rows; i++)
    {
      vec[ i ] = Source.diagVector[ i ].real();
      vec[ Source.rows + i ] = Source.diagVector[ i ].imag();
    } 
  }

  // Return the length of the vector that has been packed.
  return reqStorage;
}

template<class T>
inline void unpackGenericBlockMatrixUpdate(const std::vector<double>& vec, bool isDense, genericBlockMatrixEntry<T>& Source)
{
  // Now copy values from vec into Source.
  if (isDense)
  {
    for (int i=0; i<Source.rows; i++)
    {
      for (int j=0; j<Source.cols; j++)
      {
        //Source.denseMtx( i, j ) += genericBlockMatrixEntry<T>::val_type( vec[ j*Source.rows + i ],
        Source.denseMtx( i, j ) += T( vec[ j*Source.rows + i ],
                                                                 vec[ Source.rows*Source.cols + j*Source.rows + i ] );
      }
    }
  }
  else
  {
    for (int i=0; i<Source.rows; i++)
    {
      if ( Source.isDense() )
      {
        //Source.denseMtx( i, i ) += genericBlockMatrixEntry<T>::val_type( vec[ i ], vec[ Source.rows + i ] );
        Source.denseMtx( i, i ) += T( vec[ i ], vec[ Source.rows + i ] );
      }
      else
      {
        //Source.diagVector[ i ] += genericBlockMatrixEntry<T>::val_type( vec[ i ], vec[ Source.rows + i ] );
        Source.diagVector[ i ] += T( vec[ i ], vec[ Source.rows + i ] );
      }
    }
  }
}

inline int packHBBlockMatrix(const HBBlockMatrixEntry& Source, std::vector<double>& vec)
{
  return packGenericBlockMatrix(Source, vec);
}

inline void unpackHBBlockMatrixUpdate(const std::vector<double>& vec, bool isDense, HBBlockMatrixEntry& Source)
{
  unpackGenericBlockMatrixUpdate(vec, isDense, Source);
  return;
}

// Overload print operator.
template<class T>
inline std::ostream& operator<< (std::ostream& os, const Xyce::genericBlockMatrixEntry<T>& obj)
{
  obj.print( os );
  return os;
}
}

#ifdef Xyce_AMESOS2_BASKER

// Specialization of BASKER_ScalarTraits for block vectors
template <>
struct BASKER_ScalarTraits< Xyce::HBBlockMatrixEntry > {
  typedef Xyce::HBBlockMatrixEntry valType;
  typedef Xyce::HBBlockMatrixEntry magnitudeType;
  // Fix this one.
  static inline valType reciprocal(valType c){ return c; }
  static inline valType divide(valType a, valType b){ return a/b; }
  static inline magnitudeType approxABS(valType a) { return a; }
  static inline magnitudeType abs(valType a) { return a; }
  // Fix this one.
  static inline bool gt (valType a, valType b){ return a>b; }
};


// Specialization of BASKER_ScalarTraits for block vectors
template <>
struct BASKER_ScalarTraits< Xyce::ESBlockMatrixEntry > {
  typedef Xyce::ESBlockMatrixEntry valType;
  typedef Xyce::ESBlockMatrixEntry magnitudeType;
  // Fix this one.
  static inline valType reciprocal(valType c){ return c; }
  static inline valType divide(valType a, valType b){ return a/b; }
  static inline magnitudeType approxABS(valType a) { return a; }
  static inline magnitudeType abs(valType a) { return a; }
  // Fix this one.
  static inline bool gt (valType a, valType b){ return a>b; }
};

namespace Amesos2 {

  // Enable HBBlockMatrixEntry as a valid Scalar type for Basker
  template <>
  struct TypeMap< Basker, Xyce::HBBlockMatrixEntry > {
    static Xyce::HBBlockMatrixEntry dtype;
    typedef Xyce::HBBlockMatrixEntry type;
    typedef Xyce::HBBlockMatrixEntry magnitude_type;
  };

  // Enable ESBlockMatrixEntry as a valid Scalar type for Basker
  template <>
  struct TypeMap< Basker, Xyce::ESBlockMatrixEntry > {
    static Xyce::ESBlockMatrixEntry dtype;
    typedef Xyce::ESBlockMatrixEntry type;
    typedef Xyce::ESBlockMatrixEntry magnitude_type;
  };
}

#endif

#endif // Xyce_N_LAS_BlockMatrixEntry_h

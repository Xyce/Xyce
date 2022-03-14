//----------------------------------------------------------------------
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
// Purpose        : Epetra implementation for Block MultiVector
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Computational Sciences
//
// Creation Date  : 3/13/04
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_LAS_EpetraBlockMultiVector.h>
#include <N_LAS_BlockSystemHelpers.h>
#include <N_LAS_EpetraImporter.h>
#include <N_LAS_EpetraVector.h>
#include <N_LAS_EpetraMultiVector.h>

#include <N_PDS_EpetraParMap.h>
#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>

#include <N_ERH_ErrorMgr.h>
#include <N_UTL_FeatureTest.h>

// ---------  Other Includes  -----------

#include <Epetra_Map.h>
#include <Epetra_LocalMap.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>

#include <Teuchos_BLAS.hpp>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::EpetraBlockMultiVector
// Purpose       : constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
EpetraBlockMultiVector::EpetraBlockMultiVector( int numBlocks, int numVectors,
                                                const Teuchos::RCP<const Parallel::ParMap> & globalMap,
                                                const Teuchos::RCP<const Parallel::ParMap> & subBlockMap
                                              )
: parallelMap_(globalMap.get()),
  vecOwned_(true),
  mapOwned_(false),
  groundNode_(0.0),
  blocksViewGlobalVec_(true),
  globalBlockSize_(subBlockMap->numGlobalEntities()),
  localBlockSize_(subBlockMap->numLocalEntities()),
  numBlocks_(numBlocks),
  startBlock_(0),
  endBlock_(numBlocks),
  subBlockMap_(subBlockMap),
  blocks_(numBlocks)
{
  pdsComm_ = rcp( &(globalMap->pdsComm()),false );

  if (parallelMap_->numGlobalEntities() < 0)
  {
    Report::DevelFatal().in("EpetraBlockMultiVector::EpetraBlockMultiVector")
      << "vector length too short. Vectors must be > 0 in length.";
  }
  else if (numVectors < 1)
  {
     Report::DevelFatal().in("EpetraBlockMultiVector::EpetraBlockMultiVector")
      << "numVectors < 1";
  }

  // Create a new Petra MultiVector and set the pointer.
  const Parallel::EpetraParMap& e_map = dynamic_cast<const Parallel::EpetraParMap&>( *parallelMap_ );
  aMultiVector_ = new Epetra_MultiVector( *e_map.petraMap(), numVectors );

  //Setup Views of blocks using Block Map
  double ** Ptrs, ** Loc;
  Loc = (double**)malloc(sizeof(double*) * numVectors);

  aMultiVector_->ExtractView( &Ptrs );

  const Parallel::EpetraParMap& e_submap = dynamic_cast<const Parallel::EpetraParMap&>(*subBlockMap_);

  for( int i = 0; i < numBlocks; ++i )
  {
    for( int j = 0; j < numVectors; ++j )
    {
      Loc[j] = Ptrs[j] + localBlockSize_*i;
    }
    blocks_[i] =  Teuchos::rcp( new EpetraMultiVector( new Epetra_MultiVector( View, *(e_submap.petraMap()), Loc, numVectors ), true ) );
  }

  free(Loc);
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::~EpetraBlockMultiVector
// Purpose       : Default destructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 01/06/21
//-----------------------------------------------------------------------------
EpetraBlockMultiVector::~EpetraBlockMultiVector()
{
  if (vecOwned_)
  {
    delete aMultiVector_; aMultiVector_=0;
  }
  if (mapOwned_)
  {
    delete parallelMap_; parallelMap_=0;
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::EpetraBlockMultiVector
// Purpose       : Copy Constructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
EpetraBlockMultiVector::EpetraBlockMultiVector( const EpetraBlockMultiVector & right )
: parallelMap_( right.parallelMap_ ),
  pdsComm_( right.pdsComm_ ),
  vecOwned_(true),
  mapOwned_(false),
  groundNode_(0.0),
  blocksViewGlobalVec_(true),
  globalBlockSize_(right.globalBlockSize_),
  localBlockSize_(right.localBlockSize_),
  numBlocks_(right.numBlocks_),
  startBlock_(right.startBlock_),
  endBlock_(right.endBlock_),
  subBlockMap_(right.subBlockMap_),
  blocks_(right.numBlocks_)
{
  const Parallel::EpetraParMap& e_map = dynamic_cast<const Parallel::EpetraParMap&>( *parallelMap_ );
  aMultiVector_ = new Epetra_MultiVector( *e_map.petraMap(), numVectors() );

  //Setup Views of blocks using Block Map
  double ** Ptrs, ** Loc;
  Loc = (double**)malloc(sizeof(double*) * numVectors());

  aMultiVector_->ExtractView( &Ptrs );

  const Parallel::EpetraParMap& e_submap = dynamic_cast<const Parallel::EpetraParMap&>(*subBlockMap_);

  for( int i = 0; i < numBlocks_; ++i )
  { 
    for( int j = 0; j < numVectors(); ++j )
    { 
      Loc[j] = Ptrs[j] + localBlockSize_*i;
    }
    blocks_[i] =  Teuchos::rcp( new EpetraMultiVector( new Epetra_MultiVector( View, *(e_submap.petraMap()), Loc, numVectors() ), true ) );
  }

  free(Loc);
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::operator=
// Purpose       : assignment
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
BlockMultiVector& EpetraBlockMultiVector::operator=( const BlockMultiVector & right )
{
  if ((this != &right) && globalLength())
  {
    const EpetraVectorAccess* e_right = dynamic_cast<const EpetraVectorAccess *>( &right );
    if( (globalLength() == right.globalLength()) && (localLength() == right.localLength()) )
    { 
      *aMultiVector_ = e_right->epetraObj();
    }
    else
    {
      if (VERBOSE_LINEAR)
        Report::DevelFatal0() <<"BlockMultiVector being assigned with different map";
    }
  }

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::operator=
// Purpose       : assignment
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
BlockMultiVector& EpetraBlockMultiVector::operator=( const MultiVector & right )
{
  if ((this != &right) && globalLength())
  {
    const EpetraVectorAccess* e_right = dynamic_cast<const EpetraVectorAccess *>( &right );
    if( (globalLength() == right.globalLength()) && (localLength() == right.localLength()) )
    {
      *aMultiVector_ = e_right->epetraObj();
    }
    else
    {
      if (VERBOSE_LINEAR)
        Report::DevelFatal0() <<"BlockMultiVector being assigned with different map";
    }
  }

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::clone
// Purpose       : clone multivector
// Special Notes : clone shape, not values
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 11/18/20
//-----------------------------------------------------------------------------
MultiVector* EpetraBlockMultiVector::clone() const
{
  EpetraBlockMultiVector* new_vec = 0;
  if ( parallelMap_ )
  {
    new_vec = new EpetraBlockMultiVector( numBlocks_, numVectors(),
                                          Teuchos::rcp( parallelMap_, false ),
                                          subBlockMap_ );
  }
  else
  {
    // We don't have a map, so perform a cloneCopy
    new_vec = new EpetraBlockMultiVector( *this );
  }
  return new_vec;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::cloneCopy
// Purpose       : clone multivector
// Special Notes : clone shape and values
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 11/18/20
//-----------------------------------------------------------------------------
MultiVector* EpetraBlockMultiVector::cloneCopy() const
{
  return new EpetraBlockMultiVector( *this );
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::dotProduct
// Purpose       : Returns the dot product of "this" vector and another.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
void EpetraBlockMultiVector::dotProduct(const MultiVector & y, std::vector<double>& d) const
{
  const EpetraVectorAccess* e_y = dynamic_cast<const EpetraVectorAccess *>( &y );

  int xnum = numVectors();
  int ynum = y.numVectors();

  if (xnum == ynum)
  {
    // Let Epetra handle this.
    int PetraError = aMultiVector_->Dot(e_y->epetraObj(), &d[0]);

    if (DEBUG_LINEAR)
      processError( "EpetraBlockMultiVector::dotProduct - ", PetraError );
  } 
  else if (xnum == 1 || ynum == 1)
  {
    int maxDim = ( xnum > ynum ) ? xnum : ynum;

    Epetra_LocalMap LocalMap(maxDim, 0, (e_y->epetraObj()).Map().Comm());
    Epetra_MultiVector tmpVec(View, LocalMap, &d[0], maxDim, 1 );
    if (maxDim == xnum)
      tmpVec.Multiply('T', 'N', 1.0, *aMultiVector_, e_y->epetraObj(), 0.0 );
    else
      tmpVec.Multiply('T', 'N', 1.0, e_y->epetraObj(), *aMultiVector_, 0.0 );
  }
  else
  {
    Xyce::Report::DevelFatal().in("dotProduct") 
      << "Cannot perform dot product with vectors of dimension " << xnum << " and " << ynum;
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::scale
// Purpose       : Scales a MultiVector by a constant value.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void EpetraBlockMultiVector::scale(const double a)
{
  if ( globalLength() )
  {
    int PetraError = aMultiVector_->Scale(a);

    if (DEBUG_LINEAR)
      processError( "EpetraBlockMultiVector::scale - ", PetraError);
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::reciprocal
// Purpose       : Reciprocal elements of MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/07/01
//-----------------------------------------------------------------------------
void EpetraBlockMultiVector::reciprocal(const MultiVector & A)
{
  const EpetraVectorAccess* e_A = dynamic_cast<const EpetraVectorAccess *>( &A );
  int PetraError = aMultiVector_->Reciprocal( e_A->epetraObj() );

  if (DEBUG_LINEAR)
    processError( "EpetraBlockMultiVector::reciprocal - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::getVectorView
// Purpose       : Const view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
const Vector* EpetraBlockMultiVector::getVectorView(int index) const
{
  const Vector* vec = new EpetraVector((*aMultiVector_)(index),
                                       aMultiVector_->Map(),false);
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::getNonConstVectorView
// Purpose       : NonConst view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
Vector* EpetraBlockMultiVector::getNonConstVectorView(int index)
{
  Vector* vec = new EpetraVector((*aMultiVector_)(index),
                                 aMultiVector_->Map(),false);
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::getVectorView
// Purpose       : Const view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
const Vector* EpetraBlockMultiVector::getVectorViewAssembled(int index) const
{
  const Vector* vec = new EpetraVector( new
                          Epetra_Vector( View, *aMultiVector_, index ), true );
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::getNonConstVectorView
// Purpose       : NonConst view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
Vector* EpetraBlockMultiVector::getNonConstVectorViewAssembled(int index)
{
  Vector* vec = new EpetraVector( new
                    Epetra_Vector( View, *aMultiVector_, index ), true );
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::multiply
// Purpose       : Performs element-wise multiplication of two vectors
//                 this = this @ x
//                 where @ represents element-wise multiplication
// Special Notes :
// Scope         : Public
// Creator       : Roger P. Pawlowski, SNL, Parallel Computational Sciences
// Creation Date : 3/24/03
//-----------------------------------------------------------------------------
void EpetraBlockMultiVector::multiply(const MultiVector &x)
{
  const EpetraVectorAccess* e_x = dynamic_cast<const EpetraVectorAccess *>( &x );
  int PetraError = aMultiVector_->Multiply(1.0, *aMultiVector_,
                                           e_x->epetraObj(), 0.0);

  if (DEBUG_LINEAR)
    processError( "EpetraBlockMultiVector::multiply - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::update
// Purpose       :
// Special Notes : ERK. From the epetra documentation:
//
//                 this = s*this + a*A
//
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 02/04/02
//-----------------------------------------------------------------------------
void EpetraBlockMultiVector::update( double a, const MultiVector & A,
                                     double s )
{
  if ( globalLength() )
  {
    const EpetraVectorAccess* e_A = dynamic_cast<const EpetraVectorAccess *>( &A );
    aMultiVector_->Update( a, e_A->epetraObj(), s );
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::update
// Purpose       :
// Special Notes : ERK.  From the epetra documentation:
//
//                 this = s*this + a*A + b*B
//
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 02/04/02
//-----------------------------------------------------------------------------
void EpetraBlockMultiVector::update( double a, const MultiVector & A,
                                     double b, const MultiVector & B,
                                     double s )
{
  if ( globalLength() )
  {
    const EpetraVectorAccess* e_A = dynamic_cast<const EpetraVectorAccess *>( &A );
    const EpetraVectorAccess* e_B = dynamic_cast<const EpetraVectorAccess *>( &B );
    aMultiVector_->Update( a, e_A->epetraObj(),
                           b, e_B->epetraObj(),
                           s );
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::lpNorm
// Purpose       : Returns lp norms of each vector in MultiVector
// Special Notes : Only p=1 and p=2 implemented now since this is all Petra
//                 supports.
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
int EpetraBlockMultiVector::lpNorm(const int p, double * result) const
{
  int PetraError = -1;
  static const char *methodMsg = "EpetraBlockMultiVector::lpNorm - ";

  if (p == 1)
    PetraError = aMultiVector_->Norm1(result);
  else if (p == 2)
    PetraError = aMultiVector_->Norm2(result);
  else
    Xyce::Report::DevelFatal0().in(methodMsg) << "Requested norm is not supported";

  if (DEBUG_LINEAR)
    processError(methodMsg, PetraError);

  return PetraError;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::infNorm
// Purpose       : Returns infinity norm of each vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/16/01
//-----------------------------------------------------------------------------
int EpetraBlockMultiVector::infNorm(double * result, int * index) const
{
  static const char *methodMsg = "MultiVector::infNorm - ";

  int PetraError = 0;

  if (!index)
  {
    PetraError = aMultiVector_->NormInf(result);

    if (DEBUG_LINEAR)
      processError(methodMsg, PetraError);
  }
  else
  {
    int numProcs = pdsComm_->numProc();
    int numVectors = aMultiVector_->NumVectors();
    int myLength = aMultiVector_->MyLength();
    std::vector<int> indexTemp( numVectors, 0 ), indexTempAll( numVectors*numProcs, 0 );
    std::vector<double> doubleTemp( numVectors, 0.0 ), doubleTempAll( numVectors*numProcs, 0.0 );
    double ** pointers = aMultiVector_->Pointers();

    for (int i=0; i < numVectors; i++)
    {
      indexTemp[i] = -1;
      doubleTemp[i] = 0.0;
      for (int j=0; j < myLength; j++)
      {
        double tmp = fabs(pointers[i][j]);
        if ( tmp > doubleTemp[i] )
        {
          doubleTemp[i] = tmp;
          indexTemp[i] = j;
        }
      }
      // Convert indexTemp from local to global ID
      if (indexTemp[i] > -1)
        indexTemp[i] = aMultiVector_->Map().GID(indexTemp[i]);
    }

    if (numProcs > 1)
    {
      // Use the communicator to gather all the local maximum values and indices
      Parallel::AllGather( pdsComm_->comm(), indexTemp, indexTempAll );
      Parallel::AllGather( pdsComm_->comm(), doubleTemp, doubleTempAll );

      // Compute the global infNorm and index
      for (int i=0; i < numVectors; i++)
      {
        result[i] = doubleTempAll[i];
        index[i] = indexTempAll[i];
        for (int j=1; j < numProcs; j++)
        {
          if ( doubleTempAll[j*numVectors + i] > result[i] )
          {
            result[i] = doubleTempAll[j*numVectors + i];
            index[i] = indexTempAll[j*numVectors + i];
          }
        }
      }
    }
    else
    {
      for (int i=0; i < numVectors; i++)
      {
        result[i] = doubleTemp[i];
        index[i] = indexTemp[i];
      }     
    }
  }

  return PetraError;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::wRMSNorm
// Purpose       : Returns weighted root-mean-square of each vector in
//                 MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/12/00
//-----------------------------------------------------------------------------
int EpetraBlockMultiVector::wRMSNorm(const MultiVector & weights, double * result) const
{
  const EpetraVectorAccess* e_weights = dynamic_cast<const EpetraVectorAccess *>( &weights );
  int PetraError = aMultiVector_->NormWeighted( e_weights->epetraObj(), result );

  if (DEBUG_LINEAR)
    processError( "EpetraBlockMultiVector::wRMSNorm - ", PetraError);

  return PetraError;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::wMaxNorm
// Purpose       : Returns the weighted
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 03/19/01
//-----------------------------------------------------------------------------
int EpetraBlockMultiVector::wMaxNorm(const MultiVector & weights, double * result, int * index) const
{
  int length  = localLength();
  int numVecs = numVectors();
  int numProcs = pdsComm_->numProc();
  double tmpVal = 0.0;

  std::vector<int> indexTemp( numVecs, 0 ), indexTempAll( numVecs*numProcs, 0 );
  std::vector<double> doubleTemp( numVecs, 0.0 ), doubleTempAll( numVecs*numProcs, 0.0 );

  for (int i = 0;  i < numVecs; ++i)
  { 
    indexTemp[i] = -1;
    doubleTemp[i] = 0.0;
    if (length)
    { 
      indexTemp[i] = 0;
      doubleTemp[i] = fabs(*(*this)(0,i)) / (*weights(0,i));
      for (int j = 1; j < length; ++j)
      { 
        tmpVal = fabs(*(*this)(j,i)) / (*weights(j,i));
        if (tmpVal > doubleTemp[i])
        { 
          doubleTemp[i] = tmpVal;
          indexTemp[i] = j;
        }
      }
    }
  }

  if (numProcs > 1)
  {
    // Use the communicator to gather all the local maximum values and indices
    Parallel::AllGather( pdsComm_->comm(), indexTemp, indexTempAll );
    Parallel::AllGather( pdsComm_->comm(), doubleTemp, doubleTempAll );

    // Compute the global infNorm and index
    for (int i=0; i < numVecs; i++)
    {
      result[i] = doubleTempAll[i];
      if (index)
        index[i] = indexTempAll[i];
      for (int j=1; j < numProcs; j++)
      {
        if ( doubleTempAll[j*numVecs + i] > result[i] )
        {
          result[i] = doubleTempAll[j*numVecs + i];
          if (index)
            index[i] = indexTempAll[j*numVecs + i];
        }
      }
    }
  }
  else
  {
    for (int i=0; i < numVecs; i++)
    {
      result[i] = doubleTemp[i];
      if (index)
        index[i] = indexTemp[i];
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::putScalar
// Purpose       : Fills MultiVector with the value "scalar".
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void EpetraBlockMultiVector::putScalar(const double scalar)
{
  if ( globalLength() )
  {
    int PetraError = aMultiVector_->PutScalar(scalar);

    groundNode_ = scalar;

    if (DEBUG_LINEAR)
      processError( "EpetraBlockMultiVector::putScalar - ", PetraError);
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::addScalar
// Purpose       : Adds to MultiVector with the value "scalar".
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 01/18/01
//-----------------------------------------------------------------------------
void EpetraBlockMultiVector::addScalar(const double scalar)
{
  int length  = aMultiVector_->MyLength();
  int numVecs = numVectors();

  for (int i = 0; i < numVecs; ++i)
    for (int j = 0; j < length; ++j)
      (*aMultiVector_)[i][j] += scalar;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::absValue
// Purpose       : Abs value of elements of MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 01/18/01
//-----------------------------------------------------------------------------
void EpetraBlockMultiVector::absValue(const MultiVector & A)
{
  const EpetraVectorAccess* e_A = dynamic_cast<const EpetraVectorAccess *>( &A );
  int PetraError = aMultiVector_->Abs( e_A->epetraObj() );

  if (DEBUG_LINEAR)
    processError( "EpetraBlockMultiVector::absValue - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::vectorImport
// Purpose       : Import using Petra_Import object
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/09/01
//-----------------------------------------------------------------------------
bool EpetraBlockMultiVector::vectorImport(const MultiVector * vec,
                                          Importer * importer)
{
  EpetraImporter * e_importer = dynamic_cast<EpetraImporter *>( importer );
  const EpetraVectorAccess* e_vec = dynamic_cast<const EpetraVectorAccess *>( vec );
  aMultiVector_->Import(e_vec->epetraObj(), e_importer->epetraObj(), Insert);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::getElementByGlobalIndex
// Purpose       : Get element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/6/00
//-----------------------------------------------------------------------------
const double & EpetraBlockMultiVector::getElementByGlobalIndex(
  const int & global_index, const int & vec_index) const
{ 
  if( parallelMap_ == NULL )
    return (*aMultiVector_)[vec_index][ aMultiVector_->Map().LID(global_index) ];
  else
  { 
    int i = parallelMap_->globalToLocalIndex(global_index);
    if (i != -1)
      return ((*aMultiVector_)[vec_index])[i];
    else
      return groundNode_;
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::setElementByGlobalIndex
// Purpose       : Set element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/6/00
//-----------------------------------------------------------------------------
bool EpetraBlockMultiVector::setElementByGlobalIndex(const int & global_index,
                                                     const double & val,
                                                     const int & vec_index)
{ 
  if( parallelMap_ == NULL )
    (*aMultiVector_)[vec_index][ aMultiVector_->Map().LID(global_index) ] = val;
  else
  { 
    if (global_index != -1)
    { 
      int i = parallelMap_->globalToLocalIndex(global_index);
      if (i != -1)
      { 
        ( (*aMultiVector_)[vec_index] )[i] = val;
        return true;
      }
      else
      { 
        Xyce::Report::DevelFatal().in("setElementByGlobalIndex") 
          << "Failed to find EpetraBlockMultiVector global index: " << global_index;
        return false;
      }
    }
  }
  
  return true;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::sumElementByGlobalIndex
// Purpose       : Set element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/7/00
//-----------------------------------------------------------------------------
bool EpetraBlockMultiVector::sumElementByGlobalIndex(const int & global_index,
                                                     const double & val,
                                                     const int & vec_index)
{ 
  if( parallelMap_ == NULL )
    (*aMultiVector_)[vec_index][ aMultiVector_->Map().LID(global_index) ] += val;
  else
  { 
    if (global_index != -1 )
    { 
      int i = parallelMap_->globalToLocalIndex(global_index);
      if (i != -1)
      { 
        ( (*aMultiVector_)[vec_index] )[i] += val;
        return true;
      }
      else
      { 
        Report::DevelFatal()
          << " sumElementByGlobalIndex: failed to find EpetraBlockMultiVector global index ";
        return false;
      }
    }
  }
  
  return true;
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector::processError
// Purpose       : Concrete implementation which processes Petra (in this case)
//                 error codes taken from the Petra member function returns.
// Special Notes : Petra specific.  NOTE ALSO - this function is currently
//                 within the "Xyce_DEBUG_LINEAR" ifdef and so any calls to
//                 this should also be so bracketed.
// Scope         : Private
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
void EpetraBlockMultiVector::processError(const char *methodMsg, int error) const
{
  // Process the error
  switch (error)
  {
  case 0:
    Xyce::dout() << methodMsg << ": Function returned without warnings or errors." << std::endl;
    break;

  default:
    Xyce::Report::DevelFatal0().in(methodMsg) << "Function returned with an error.";
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraBlockMultiVector:::print
// Purpose       : Output
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
void EpetraBlockMultiVector::print(std::ostream &os) const
{
  os << "EpetraBlockMultiVector Object (Number of Blocks =" << numBlocks_ << ", Number of Vectors =" << numVectors() << ", View =" << blocksViewGlobalVec_ << std::endl;

  os << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  for( int i = 0; i < numBlocks_; ++i )
  {
    if (i >= startBlock_ && i < endBlock_)
    {
      os << "Block[" << i << "]\n";
    }
    blocks_[i]->print( os );
  }
  os << "Base Object\n";
  os << *aMultiVector_;
  os << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

} // namespace Linear
} // namespace Xyce

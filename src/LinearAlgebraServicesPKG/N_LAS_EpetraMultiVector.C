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

//-------------------------------------------------------------------------
//
// Purpose        : Implementation file for the Abstract interface to the
//                  multi-vector types (RDP, RSP, CDP or CSP).
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
#include <cstdio>

// ----------   Xyce Includes   ----------
#include <N_UTL_fwd.h>

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_EpetraMultiVector.h>
#include <N_LAS_EpetraVector.h>
#include <N_LAS_EpetraImporter.h>
#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>
#include <N_PDS_EpetraHelpers.h>
#include <N_PDS_EpetraParMap.h>
#include <N_UTL_FeatureTest.h>

// ---------  Other Includes  -----------

#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include <Epetra_Map.h>

#include <EpetraExt_View_MultiVector.h>
#include <EpetraExt_MultiVectorOut.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::EpetraMultiVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
EpetraMultiVector::EpetraMultiVector(const Parallel::ParMap & map, int numVectors)
:  parallelMap_(&map),
   overlapMap_(&map),
   importer_(0),
   exporter_(0),
   viewTransform_(0),
   vecOwned_(true),
   mapOwned_(false),
   groundNode_(0.0)
{
  pdsComm_ = rcp( &map.pdsComm(),false ); 

  if (map.numGlobalEntities() < 0)
  {
    Report::DevelFatal().in("EpetraMultiVector::EpetraMultiVector")
      << "vector length too short. Vectors must be > 0 in length.";
  }
  else if (numVectors < 1)
  {
     Report::DevelFatal().in("EpetraMultiVector::EpetraMultiVector")
      << "numVectors < 1";
  }

  // Create a new Petra EpetraMultiVector and set the pointer.
  const Parallel::EpetraParMap& e_map = dynamic_cast<const Parallel::EpetraParMap&>( map );
  aMultiVector_ = new Epetra_MultiVector( *e_map.petraMap(), numVectors );

  oMultiVector_ = aMultiVector_;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::EpetraMultiVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/06/02
//-----------------------------------------------------------------------------
EpetraMultiVector::EpetraMultiVector(
  const Parallel::ParMap &        map,
  const Parallel::ParMap &        ol_map,
  int                   numVectors )
  : parallelMap_(&map),
    overlapMap_(&ol_map),
    importer_(0),
    exporter_(0),
    viewTransform_(0),
    vecOwned_(true),
    mapOwned_(false),
    groundNode_(0.0)
{
  pdsComm_ = rcp( &map.pdsComm(),false ); 

  if (map.numGlobalEntities() < 0)
    Report::DevelFatal().in("EpetraMultiVector::EpetraMultiVector")
      << "vector length too short. Vectors must be > 0 in length.";

  // Create a new Petra MultiVector and set the pointer.
  const Parallel::EpetraParMap& e_map = dynamic_cast<const Parallel::EpetraParMap&>( map );
  const Parallel::EpetraParMap& e_ol_map = dynamic_cast<const Parallel::EpetraParMap&>( ol_map );
  oMultiVector_ = new Epetra_MultiVector( *e_ol_map.petraMap(), numVectors);

  viewTransform_ = new EpetraExt::MultiVector_View( *e_ol_map.petraMap(), *e_map.petraMap() );
  aMultiVector_ = &((*viewTransform_)(*oMultiVector_));
  if (map.pdsComm().numProc() > 1)
  {
    exporter_ = new Epetra_Export( *e_ol_map.petraMap(), *e_map.petraMap() );
  }

  importer_ = new Epetra_Import( *e_ol_map.petraMap(), *e_map.petraMap() );
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::EpetraMultiVector
// Purpose       : Copy Constructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
EpetraMultiVector::EpetraMultiVector( const EpetraMultiVector & right )
: parallelMap_( right.pmap() ),
  overlapMap_( right.pmap() ),
  oMultiVector_( new Epetra_MultiVector( *(right.oMultiVector_) ) ),
  importer_(0),
  exporter_(0),
  viewTransform_(0),
  pdsComm_( Teuchos::rcp( right.pdsComm(), false ) ),
  vecOwned_(true),
  mapOwned_(false),
  groundNode_(0.0)
{
  if (right.aMultiVector_ == right.oMultiVector_)
    aMultiVector_ = oMultiVector_;
  else
  {
    const Parallel::EpetraParMap* e_map = dynamic_cast<const Parallel::EpetraParMap*>( parallelMap_ );
    const Parallel::EpetraParMap* e_ol_map = dynamic_cast<const Parallel::EpetraParMap*>( overlapMap_ );

    viewTransform_ = new EpetraExt::MultiVector_View( *e_ol_map->petraMap(), *e_map->petraMap() );
    aMultiVector_ = &((*viewTransform_)( *oMultiVector_ ));
  }

  // Generate new exporter instead of using copy constructor, there is an issue with Epetra_MpiDistributor
  if( right.exporter_ ) 
  {
    const Parallel::EpetraParMap* e_map = dynamic_cast<const Parallel::EpetraParMap*>( parallelMap_ );
    const Parallel::EpetraParMap* e_ol_map = dynamic_cast<const Parallel::EpetraParMap*>( overlapMap_ );

    exporter_ = new Epetra_Export( *e_ol_map->petraMap(), *e_map->petraMap() );
  }
  if( right.importer_ ) importer_ = new Epetra_Import( *right.importer_ );
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::EpetraMultiVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
EpetraMultiVector::EpetraMultiVector( Epetra_MultiVector * overlapMV, const Epetra_BlockMap& parMap, bool isOwned )
: parallelMap_(0),
  overlapMap_(0),
  oMultiVector_( overlapMV ),
  importer_(0),
  exporter_(0),
  viewTransform_(0),
  vecOwned_(isOwned),
  mapOwned_(false),
  groundNode_(0.0) 
{
  pdsComm_ = Teuchos::rcp( Xyce::Parallel::createPDSComm( &overlapMV->Comm() ) );  

  // Make sure there is anything to communicate before creating a transform, importer, or exporter
  if (overlapMV->MyLength() == parMap.NumMyElements())
  {
    aMultiVector_ = oMultiVector_;
  }
  else
  {
    viewTransform_ = new EpetraExt::MultiVector_View( overlapMV->Map(), parMap );
    aMultiVector_ = &((*viewTransform_)(*oMultiVector_));
    if( pdsComm_->numProc() > 1 )
      exporter_ = new Epetra_Export( overlapMV->Map(), parMap );

    importer_ = new Epetra_Import( overlapMV->Map(), parMap );
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::EpetraMultiVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
EpetraMultiVector::EpetraMultiVector( Epetra_MultiVector * origMV, bool isOwned )
: parallelMap_(0),
  overlapMap_(0),
  aMultiVector_( origMV ),
  oMultiVector_( origMV ),
  importer_(0),
  exporter_(0),
  viewTransform_(0),
  vecOwned_(isOwned),
  mapOwned_(false),
  groundNode_(0.0)
{
  pdsComm_ = Teuchos::rcp( Xyce::Parallel::createPDSComm( &origMV->Comm() ) );
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::operator=
// Purpose       : assignment
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
MultiVector & EpetraMultiVector::operator=( const MultiVector & right )
{
  if (this != &right && globalLength())
  {
    const EpetraVectorAccess* e_right = dynamic_cast<const EpetraVectorAccess *>( &right );
    const Epetra_MultiVector & e_aMV = e_right->epetraObj();
    const Epetra_MultiVector & e_oMV = e_right->epetraOverlapObj();
    if( (oMultiVector_->Map().NumGlobalElements() == e_oMV.Map().NumGlobalElements())
        && (oMultiVector_->Map().NumMyElements() == e_oMV.Map().NumMyElements()) )
    {
      *oMultiVector_ = e_oMV;
    }

    if( (globalLength() == right.globalLength()) && (localLength() == right.localLength()) )
    {
      *aMultiVector_ = e_aMV;
    }
    else
    {
      if (VERBOSE_LINEAR)
        Report::DevelFatal0() <<"EpetraMultiVector being assigned with different mapping";
    }
  }

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::~EpetraMultiVector
// Purpose       : Default destructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
EpetraMultiVector::~EpetraMultiVector()
{
  delete importer_;
  delete exporter_;
  delete viewTransform_; //destroys of aMultiVector_ as well
  if (vecOwned_)
  {
    delete oMultiVector_;
  }
  if (mapOwned_)
  {
    delete parallelMap_; parallelMap_=0;
    if (overlapMap_)
      delete overlapMap_;
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::clone
// Purpose       : clone multivector
// Special Notes : clone shape, not values
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 11/18/20
//-----------------------------------------------------------------------------
MultiVector* EpetraMultiVector::clone() const
{
  MultiVector* new_vec = 0;
  if ( parallelMap_ )
  {
    if ( parallelMap_ == overlapMap_ )
      new_vec = new EpetraMultiVector( *parallelMap_, this->numVectors() );
    else
      new_vec = new EpetraMultiVector( *parallelMap_, *overlapMap_, this->numVectors() );
  }
  else
  {
    // We don't have a map, so perform a cloneCopy
    new_vec = new EpetraMultiVector( *this );
  }
  return new_vec;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::cloneCopy
// Purpose       : clone multivector
// Special Notes : clone shape and values
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 11/18/20
//-----------------------------------------------------------------------------
MultiVector* EpetraMultiVector::cloneCopy() const
{
  return new EpetraMultiVector( *this );
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::dotProduct
// Purpose       : Returns the dot product of "this" vector and another.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
void EpetraMultiVector::dotProduct(const MultiVector & y, std::vector<double>& d) const
{
  const EpetraVectorAccess* e_y = dynamic_cast<const EpetraVectorAccess *>( &y );

  int xnum = numVectors();
  int ynum = y.numVectors();
  if (xnum == 1 || ynum == 1)
  {
    for (int i=0; i<xnum; ++i)
    {
      for (int j=0; j<ynum; ++j)
      {
        (*aMultiVector_)(i)->Dot(*(e_y->epetraObj())(j), &d[i*ynum + j]);
      }
    }
  }
  else
  {    
    // Let Epetra handle this.
    int PetraError = aMultiVector_->Dot(e_y->epetraObj(), &d[0]);
  
    if (DEBUG_LINEAR)
      processError( "EpetraMultiVector::dotProduct - ", PetraError );
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::scale
// Purpose       : Scales a MultiVector by a constant value.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void EpetraMultiVector::scale(const double a)
{
  int PetraError = aMultiVector_->Scale(a);

  if (DEBUG_LINEAR)
    processError( "EpetraMultiVector::scale - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::multiply
// Purpose       : Performs element-wise multiplication of two vectors
//                 this = this @ x
//                 where @ represents element-wise multiplication
// Special Notes :
// Scope         : Public
// Creator       : Roger P. Pawlowski, SNL, Parallel Computational Sciences
// Creation Date : 3/24/03
//-----------------------------------------------------------------------------
void EpetraMultiVector::multiply(const MultiVector &x)
{
  const EpetraVectorAccess* e_x = dynamic_cast<const EpetraVectorAccess *>( &x );
  int PetraError = aMultiVector_->Multiply(1.0, *aMultiVector_,
					   e_x->epetraObj(), 0.0);

  if (DEBUG_LINEAR)
    processError( "EpetraMultiVector::multiply - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::update
// Purpose       :
// Special Notes : ERK. From the epetra documentation:
//
//                 this = s*this + a*A
//
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 02/04/02
//-----------------------------------------------------------------------------
void EpetraMultiVector::update( double a, const MultiVector & A,
                          double s )
{
  const EpetraVectorAccess* e_A = dynamic_cast<const EpetraVectorAccess *>( &A );
  aMultiVector_->Update( a, e_A->epetraObj(), s );
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::update
// Purpose       :
// Special Notes : ERK.  From the epetra documentation:
//
//                 this = s*this + a*A + b*B
//
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 02/04/02
//-----------------------------------------------------------------------------
void EpetraMultiVector::update( double a, const MultiVector & A,
                          double b, const MultiVector & B,
                          double s )
{
  const EpetraVectorAccess* e_A = dynamic_cast<const EpetraVectorAccess *>( &A );
  const EpetraVectorAccess* e_B = dynamic_cast<const EpetraVectorAccess *>( &B );
  aMultiVector_->Update( a, e_A->epetraObj(),
                         b, e_B->epetraObj(),
                         s );
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::lpNorm
// Purpose       : Returns lp norms of each vector in MultiVector
// Special Notes : Only p=1 and p=2 implemented now since this is all Petra
//                 supports.
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
int EpetraMultiVector::lpNorm(const int p, double * result) const
{
  int PetraError = -1;
  static const char *methodMsg = "EpetraMultiVector::lpNorm - ";

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
// Function      : EpetraMultiVector::infNorm
// Purpose       : Returns infinity norm of each vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/16/01
//-----------------------------------------------------------------------------
int EpetraMultiVector::infNorm(double * result, int * index) const
{
  static const char *methodMsg = "EpetraMultiVector::infNorm - ";

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
        double tmp = std::abs(pointers[i][j]);
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
// Function      : EpetraMultiVector::wRMSNorm
// Purpose       : Returns weighted root-mean-square of each vector in
//                 MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/12/00
//-----------------------------------------------------------------------------
int EpetraMultiVector::wRMSNorm(const MultiVector & weights, double * result) const
{
  const EpetraVectorAccess* e_weights = dynamic_cast<const EpetraVectorAccess *>( &weights );
  int PetraError = aMultiVector_->NormWeighted( e_weights->epetraObj(), result );

  if (DEBUG_LINEAR)
    processError( "EpetraMultiVector::wRMSNorm - ", PetraError);

  return PetraError;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::wMaxNorm
// Purpose       : Returns the weighted inf-norm
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 03/19/01
//-----------------------------------------------------------------------------
int EpetraMultiVector::wMaxNorm(const MultiVector & weights, double * result) const
{
  int length  = aMultiVector_->MyLength();
  int numVecs = numVectors();
  double tmpVal = 0.0;

  for (int i = 0;  i < numVecs; ++i)
  {
    double localMax = 0.0;
    if (length)
    { 
      localMax = fabs(*(*this)(0,i)) / (*weights(0,i)); 
      for (int j = 1; j < length; ++j)
      {
        tmpVal = fabs(*(*this)(j,i)) / (*weights(j,i));
        if (tmpVal > localMax)
          localMax = tmpVal;
      }
    } 
    // Determine global maximum.
    pdsComm_->maxAll( &localMax, &(result[i]), 1 );
  }

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::random
// Purpose       : Generates random numbers drawn from a uniform distribution
//                 on the interval (-1,1) using a multiplicative congruential
//                 generator with modulus 2^31 - 1.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void EpetraMultiVector::random()
{
  int PetraError = aMultiVector_->Random();

  if (DEBUG_LINEAR)
    processError( "EpetraMultiVector::random - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::putScalar
// Purpose       : Fills MultiVector with the value "scalar".
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void EpetraMultiVector::putScalar(const double scalar)
{
  int PetraError = oMultiVector_->PutScalar(scalar);

  groundNode_ = scalar;

  if (DEBUG_LINEAR)
    processError( "EpetraMultiVector::putScalar - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::addScalar
// Purpose       : Adds to MultiVector with the value "scalar".
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 01/18/01
//-----------------------------------------------------------------------------
void EpetraMultiVector::addScalar(const double scalar)
{
  int length  = aMultiVector_->MyLength();
  int numVecs = numVectors();

  for (int i = 0; i < numVecs; ++i)
    for (int j = 0; j < length; ++j)
      (*aMultiVector_)[i][j] += scalar;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::absValue
// Purpose       : Abs value of elements of MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 01/18/01
//-----------------------------------------------------------------------------
void EpetraMultiVector::absValue(const MultiVector & A)
{
  const EpetraVectorAccess* e_A = dynamic_cast<const EpetraVectorAccess *>( &A );
  int PetraError = oMultiVector_->Abs(e_A->epetraOverlapObj());

  if (DEBUG_LINEAR)
    processError( "EpetraMultiVector::absValue - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::reciprocal
// Purpose       : Reciprocal elements of MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/07/01
//-----------------------------------------------------------------------------
void EpetraMultiVector::reciprocal(const MultiVector & A)
{
  const EpetraVectorAccess* e_A = dynamic_cast<const EpetraVectorAccess *>( &A );
  int PetraError = oMultiVector_->Reciprocal(e_A->epetraOverlapObj());

  if (DEBUG_LINEAR)
    processError( "EpetraMultiVector::reciprocal - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::getVectorView
// Purpose       : Const view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
const Vector* EpetraMultiVector::getVectorView(int index) const
{
  const Vector* vec = new EpetraVector((*oMultiVector_)(index),
                                  aMultiVector_->Map(),false);
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::getNonConstVectorView
// Purpose       : NonConst view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
Vector* EpetraMultiVector::getNonConstVectorView(int index)
{
  Vector* vec = new EpetraVector((*oMultiVector_)(index),
                            aMultiVector_->Map(),false);
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::getVectorView
// Purpose       : Const view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
const Vector* EpetraMultiVector::getVectorViewAssembled(int index) const
{
  const Vector* vec = new EpetraVector( new 
                      Epetra_Vector( View, *aMultiVector_, index ), true );
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::getNonConstVectorView
// Purpose       : NonConst view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
Vector* EpetraMultiVector::getNonConstVectorViewAssembled(int index)
{
  Vector* vec = new EpetraVector( new 
                Epetra_Vector( View, *aMultiVector_, index ), true );
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::fillComplete
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/29/03
//-----------------------------------------------------------------------------
void EpetraMultiVector::fillComplete()
{
  if ( exporter_ )
  {
    aMultiVector_->Export( *oMultiVector_, *exporter_, Add );
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::vectorImport
// Purpose       : Import using Petra_Import object
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/09/01
//-----------------------------------------------------------------------------
bool EpetraMultiVector::vectorImport(const MultiVector * vec,
                               Importer * importer)
{
  const EpetraVectorAccess* e_vec = dynamic_cast<const EpetraVectorAccess *>( vec );
  EpetraImporter * e_importer = dynamic_cast<EpetraImporter *>( importer );
  aMultiVector_->Import(e_vec->epetraObj(), e_importer->epetraObj(), Insert);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::importOverlap
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/06/02
//-----------------------------------------------------------------------------
bool EpetraMultiVector::importOverlap()
{
  bool flag = false;

  if( importer_ )
    flag = oMultiVector_->Import( *aMultiVector_, *importer_, Insert );

  return flag;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::getElementByGlobalIndex
// Purpose       : Get element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/6/00
//-----------------------------------------------------------------------------
const double & EpetraMultiVector::getElementByGlobalIndex(
  const int & global_index, const int & vec_index) const
{
  if( aMultiVector_ != oMultiVector_ )
    return (*oMultiVector_)[vec_index][overlapMap_->globalToLocalIndex(global_index)];
  else if( parallelMap_ == NULL )
    return (*aMultiVector_)[vec_index][ aMultiVector_->Map().LID(global_index) ];
  else
  {
    int i = parallelMap_->globalToLocalIndex(global_index);

    if (i != -1)
      return ((*aMultiVector_)[vec_index])[i];
    else {
      std::map<int,double>::const_iterator it = externVectorMap_.find(global_index);
      if (it != externVectorMap_.end())
        return (*it).second;
      else
      {
        char message[128];
        sprintf(message, "getElementByGlobalIndex: failed to find MultiVector "
                "global index. global_index = %d", global_index);
        std::string msg(message);

        Report::DevelFatal() << msg;
        return (*externVectorMap_.find(-1)).second;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::setElementByGlobalIndex
// Purpose       : Set element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/6/00
//-----------------------------------------------------------------------------
bool EpetraMultiVector::setElementByGlobalIndex(const int & global_index,
                                                const double & val,
                                                const int & vec_index)
{
  if( aMultiVector_ != oMultiVector_ )
    (*oMultiVector_)[vec_index][overlapMap_->globalToLocalIndex(global_index)] = val;
  else if( parallelMap_ == NULL )
    (*oMultiVector_)[vec_index][ oMultiVector_->Map().LID(global_index) ] = val;
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
        Xyce::Report::DevelFatal().in("setElementByGlobalIndex") << "Failed to find MultiVector global index: " << global_index;
        return false;
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : EpetraMultiVector::sumElementByGlobalIndex
// Purpose       : Set element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/7/00
//-----------------------------------------------------------------------------
bool EpetraMultiVector::sumElementByGlobalIndex(const int & global_index,
                                                const double & val,
                                                const int & vec_index)
{
  if( aMultiVector_ != oMultiVector_ )
    (*oMultiVector_)[vec_index][overlapMap_->globalToLocalIndex(global_index)] += val;
  else if( parallelMap_ == NULL )
    (*oMultiVector_)[vec_index][ oMultiVector_->Map().LID(global_index) ] += val;
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
          << " sumElementByGlobalIndex: failed to find MultiVector global index ";
        return false;
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : print
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/14/00
//-----------------------------------------------------------------------------
void EpetraMultiVector::print(std::ostream &os) const
{
  if (aMultiVector_ != oMultiVector_)
  {
    os << *aMultiVector_;
  }
  os << *oMultiVector_;

  if (VERBOSE_LINEAR)
  {
    std::map<int, double>::const_iterator it_idM = externVectorMap_.begin();
    std::map<int, double>::const_iterator end_idM = externVectorMap_.end();
    if (it_idM != end_idM) os << "<Extern Vector Map>" << std::endl;
    for (; it_idM != end_idM; ++it_idM)
    {
      os << "  " << it_idM->first << "\t" << it_idM->second << std::endl;
    }
    os << std::endl;
  }
}

} // namespace Linear

} // namespace Xyce

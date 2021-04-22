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

//-------------------------------------------------------------------------
//
// Purpose        : Implementation file for the Abstract interface to the
//                  vector types (RDP, RSP, CDP or CSP).
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Computational Sciences
//
// Creation Date  : 10/13/00
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ----------   Xyce Includes   ----------
#include <N_UTL_fwd.h>

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_EpetraVector.h>
#include <N_LAS_EpetraImporter.h>
#include <N_LAS_EpetraVector.h>
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

// ---------- Standard Includes ----------

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : EpetraVector::EpetraVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
EpetraVector::EpetraVector(const Parallel::ParMap & map)
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
    Report::DevelFatal().in("EpetraVector::EpetraVector")
      << "vector length too short. Vectors must be > 0 in length.";
  }

  // Create a new Petra MultiVector and set the pointer.
  const Parallel::EpetraParMap& e_map = dynamic_cast<const Parallel::EpetraParMap&>( map );
  aMultiVector_ = new Epetra_MultiVector( *e_map.petraMap(), 1 );

  oMultiVector_ = aMultiVector_;
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::EpetraVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/06/02
//-----------------------------------------------------------------------------
EpetraVector::EpetraVector( const Parallel::ParMap &        map,
                            const Parallel::ParMap &        ol_map )
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
    Report::DevelFatal().in("EpetraVector::EpetraVector")
      << "vector length too short. Vectors must be > 0 in length.";

  // Create a new Petra MultiVector and set the pointer.
  const Parallel::EpetraParMap& e_map = dynamic_cast<const Parallel::EpetraParMap&>( map );
  const Parallel::EpetraParMap& e_ol_map = dynamic_cast<const Parallel::EpetraParMap&>( ol_map );
  oMultiVector_ = new Epetra_MultiVector( *e_ol_map.petraMap(), 1 );

  viewTransform_ = new EpetraExt::MultiVector_View( *e_ol_map.petraMap(), *e_map.petraMap() );
  aMultiVector_ = &((*viewTransform_)(*oMultiVector_));
  if (map.pdsComm().numProc() > 1)
  {
    exporter_ = new Epetra_Export( *e_ol_map.petraMap(), *e_map.petraMap() );
  }

  importer_ = new Epetra_Import( *e_ol_map.petraMap(), *e_map.petraMap() );
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::EpetraVector
// Purpose       : Copy Constructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
EpetraVector::EpetraVector( const EpetraVector& right )
: parallelMap_( right.parallelMap_ ),
  overlapMap_( right.overlapMap_ ),
  aMultiVector_(0),
  oMultiVector_(0),
  importer_(0),
  exporter_(0),
  viewTransform_(0),
  pdsComm_( right.pdsComm_ ),
  vecOwned_(true),
  mapOwned_(false),
  groundNode_(0.0)
{
  oMultiVector_ = new Epetra_MultiVector( *right.oMultiVector_ );

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
// Function      : EpetraVector::EpetraVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
EpetraVector::EpetraVector( Epetra_Vector * overlapMV, const Epetra_BlockMap& parMap, bool isOwned )
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
// Function      : EpetraVector::EpetraVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
EpetraVector::EpetraVector( Epetra_Vector * origMV, bool isOwned )
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
// Function      : cloneVector
// Purpose       : vector clone function 
// Special Notes : clones shape, not values
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
Vector* EpetraVector::cloneVector() const
{
  EpetraVector* new_vec = 0;
  if ( parallelMap_ )
  {
    if ( parallelMap_ == overlapMap_ )
      new_vec = new EpetraVector( *parallelMap_ );
    else
      new_vec = new EpetraVector( *parallelMap_, *overlapMap_ );
  }
  else
  {
    // We don't have a map, so perform a cloneCopy
    new_vec = new EpetraVector( *this );
  }
  return new_vec;
}
  
//-----------------------------------------------------------------------------
// Function      : cloneCopyVector
// Purpose       : vector clone function 
// Special Notes : clones shape and values
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
Vector* EpetraVector::cloneCopyVector() const
{
  return new EpetraVector( *this );
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::operator=
// Purpose       : assignment
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
Vector& EpetraVector::operator=( const Vector& right )
{
  if (this != &right && globalLength())
  {
    const EpetraVectorAccess* e_right = dynamic_cast<const EpetraVectorAccess *>( &right );
    const Epetra_MultiVector & e_aMV = e_right->epetraObj();
    const Epetra_MultiVector & e_oMV = e_right->epetraOverlapObj();

    if( (oMultiVector_->Map().NumGlobalElements() == e_aMV.Map().NumGlobalElements())
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
        Report::DevelFatal0() <<"EpetraVector being assigned with different map";
    }
  }

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::~EpetraVector
// Purpose       : Default destructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
EpetraVector::~EpetraVector()
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
// Function      : EpetraVector::dotProduct
// Purpose       : Returns the dot product of "this" vector and another.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 02/08/17
//-----------------------------------------------------------------------------
double EpetraVector::dotProduct( const Vector & y ) const
{
  double result = 0.0;
  const EpetraVectorAccess* e_y = dynamic_cast<const EpetraVectorAccess *>( &y );
  int PetraError = aMultiVector_->Dot(e_y->epetraObj(), &result);

  if (DEBUG_LINEAR)
    processError( "EpetraVector::dotProduct - ", PetraError );

  return result;
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::dotProduct
// Purpose       : Returns the dot product of "this" vector and another.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
void EpetraVector::dotProduct(const MultiVector & y, std::vector<double>& d) const
{
  const EpetraVectorAccess* e_y = dynamic_cast<const EpetraVectorAccess *>( &y );
  int ynum = y.numVectors();
  for (int j=0; j<ynum; ++j)
  {
    aMultiVector_->Dot(*(e_y->epetraObj()(j)), &d[j]);
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::scale
// Purpose       : Scales a MultiVector by a constant value.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void EpetraVector::scale(const double a)
{
  if (globalLength()) 
  {
    int PetraError = aMultiVector_->Scale(a);

    if (DEBUG_LINEAR)
      processError( "EpetraVector::scale - ", PetraError);
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::multiply
// Purpose       : Performs element-wise multiplication of two vectors
//                 this = this @ x
//                 where @ represents element-wise multiplication
// Special Notes :
// Scope         : Public
// Creator       : Roger P. Pawlowski, SNL, Parallel Computational Sciences
// Creation Date : 3/24/03
//-----------------------------------------------------------------------------
void EpetraVector::multiply(const MultiVector &x)
{
  const EpetraVectorAccess* e_x = dynamic_cast<const EpetraVectorAccess *>( &x );
  int PetraError = aMultiVector_->Multiply(1.0, *aMultiVector_,
					   e_x->epetraObj(), 0.0);

  if (DEBUG_LINEAR)
    processError( "EpetraVector::multiply - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::update
// Purpose       :
// Special Notes : ERK. From the epetra documentation:
//
//                 this = s*this + a*A
//
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 02/04/02
//-----------------------------------------------------------------------------
void EpetraVector::update( double a, const MultiVector & A,
                          double s )
{
  if (globalLength())
  {
    const EpetraVectorAccess* e_A = dynamic_cast<const EpetraVectorAccess *>( &A );
    aMultiVector_->Update( a, e_A->epetraObj(), s );
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::update
// Purpose       :
// Special Notes : ERK.  From the epetra documentation:
//
//                 this = s*this + a*A + b*B
//
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 02/04/02
//-----------------------------------------------------------------------------
void EpetraVector::update( double a, const MultiVector & A,
                          double b, const MultiVector & B,
                          double s )
{
  if (globalLength()) 
  {
    const EpetraVectorAccess* e_A = dynamic_cast<const EpetraVectorAccess *>( &A );
    const EpetraVectorAccess* e_B = dynamic_cast<const EpetraVectorAccess *>( &B );
    aMultiVector_->Update( a, e_A->epetraObj(),
                           b, e_B->epetraObj(), 
                           s );
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::lpNorm
// Purpose       : Returns lp norms of each vector in MultiVector
// Special Notes : Only p=1 and p=2 implemented now since this is all Petra
//                 supports.
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
int EpetraVector::lpNorm(const int p, double * result) const
{
  int PetraError = -1;
  static const char *methodMsg = "EpetraVector::lpNorm - ";

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
// Function      : EpetraVector::infNorm
// Purpose       : Returns infinity norm of each vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/16/01
//-----------------------------------------------------------------------------
int EpetraVector::infNorm(double * result, int * index) const
{
  static const char *methodMsg = "EpetraVector::infNorm - ";

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
// Function      : EpetraVector::wRMSNorm
// Purpose       : Returns weighted root-mean-square of each vector in
//                 MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/12/00
//-----------------------------------------------------------------------------
int EpetraVector::wRMSNorm(const MultiVector & weights, double * result) const
{
  const EpetraVectorAccess* e_weights = dynamic_cast<const EpetraVectorAccess *>( &weights );
  int PetraError = aMultiVector_->NormWeighted( e_weights->epetraObj(), result );

  if (DEBUG_LINEAR)
    processError( "EpetraVector::wRMSNorm - ", PetraError);

  return PetraError;
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::wMaxNorm
// Purpose       : Returns the weighted inf-norm
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 03/19/01
//-----------------------------------------------------------------------------
int EpetraVector::wMaxNorm(const MultiVector & weights, double * result) const
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
// Function      : EpetraVector::random
// Purpose       : Generates random numbers drawn from a uniform distribution
//                 on the interval (-1,1) using a multiplicative congruential
//                 generator with modulus 2^31 - 1.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void EpetraVector::random()
{
  int PetraError = aMultiVector_->Random();

  if (DEBUG_LINEAR)
    processError( "EpetraVector::random - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::putScalar
// Purpose       : Fills MultiVector with the value "scalar".
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void EpetraVector::putScalar(const double scalar)
{
  if (globalLength()) 
  {
    int PetraError = oMultiVector_->PutScalar(scalar);

    groundNode_ = scalar;

    if (DEBUG_LINEAR)
      processError( "EpetraVector::putScalar - ", PetraError);
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::addScalar
// Purpose       : Adds to MultiVector with the value "scalar".
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 01/18/01
//-----------------------------------------------------------------------------
void EpetraVector::addScalar(const double scalar)
{
  int length  = aMultiVector_->MyLength();
  int numVecs = numVectors();

  for (int i = 0; i < numVecs; ++i)
    for (int j = 0; j < length; ++j)
      (*aMultiVector_)[i][j] += scalar;
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::absValue
// Purpose       : Abs value of elements of MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 01/18/01
//-----------------------------------------------------------------------------
void EpetraVector::absValue(const MultiVector & A)
{
  const EpetraVectorAccess* e_A = dynamic_cast<const EpetraVectorAccess *>( &A );
  int PetraError = oMultiVector_->Abs(e_A->epetraOverlapObj());

  if (DEBUG_LINEAR)
    processError( "EpetraVector::absValue - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::reciprocal
// Purpose       : Reciprocal elements of MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/07/01
//-----------------------------------------------------------------------------
void EpetraVector::reciprocal(const MultiVector & A)
{
  const EpetraVectorAccess* e_A = dynamic_cast<const EpetraVectorAccess *>( &A );
  int PetraError = oMultiVector_->Reciprocal(e_A->epetraOverlapObj());

  if (DEBUG_LINEAR)
    processError( "EpetraVector::reciprocal - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::getVectorView
// Purpose       : Const view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
const Vector* EpetraVector::getVectorView(int index) const
{
  const Vector* vec = new EpetraVector((*oMultiVector_)(index),
                                       aMultiVector_->Map(),false);
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::getNonConstVectorView
// Purpose       : NonConst view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
Vector* EpetraVector::getNonConstVectorView(int index)
{
  Vector* vec = new EpetraVector((*oMultiVector_)(index),
                                 aMultiVector_->Map(),false);
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::getVectorView
// Purpose       : Const view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
const Vector* EpetraVector::getVectorViewAssembled(int index) const
{
  const Vector* vec = new EpetraVector( new 
                      Epetra_Vector( View, *aMultiVector_, index ), true );
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::getNonConstVectorView
// Purpose       : NonConst view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
Vector* EpetraVector::getNonConstVectorViewAssembled(int index)
{
  Vector* vec = new EpetraVector( new 
                Epetra_Vector( View, *aMultiVector_, index ), true );
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::fillComplete
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/29/03
//-----------------------------------------------------------------------------
void EpetraVector::fillComplete()
{
  if ( exporter_ )
  {
    aMultiVector_->Export( *oMultiVector_, *exporter_, Add );
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::vectorImport
// Purpose       : Import using Petra_Import object
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/09/01
//-----------------------------------------------------------------------------
bool EpetraVector::vectorImport(const MultiVector * vec,
                               Importer * importer)
{
  EpetraImporter * e_importer = dynamic_cast<EpetraImporter *>( importer );
  const EpetraVectorAccess* e_vec = dynamic_cast<const EpetraVectorAccess *>( vec );
  aMultiVector_->Import(e_vec->epetraObj(), e_importer->epetraObj(), Insert);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::importOverlap
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/06/02
//-----------------------------------------------------------------------------
bool EpetraVector::importOverlap()
{
  bool flag = false;

  if( importer_ )
    flag = oMultiVector_->Import( *aMultiVector_, *importer_, Insert );

  return flag;
}

//-----------------------------------------------------------------------------
// Function      : EpetraVector::getElementByGlobalIndex
// Purpose       : Get element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/6/00
//-----------------------------------------------------------------------------
const double & EpetraVector::getElementByGlobalIndex(
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
// Function      : EpetraVector::setElementByGlobalIndex
// Purpose       : Set element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/6/00
//-----------------------------------------------------------------------------
bool EpetraVector::setElementByGlobalIndex(const int & global_index,
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
// Function      : EpetraVector::sumElementByGlobalIndex
// Purpose       : Set element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/7/00
//-----------------------------------------------------------------------------
bool EpetraVector::sumElementByGlobalIndex(const int & global_index,
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
// Function      : EpetraVector::processError
// Purpose       : Concrete implementation which processes Petra (in this case)
//                 error codes taken from the Petra member function returns.
// Special Notes : Petra specific.  NOTE ALSO - this function is currently
//                 within the "Xyce_DEBUG_LINEAR" ifdef and so any calls to
//                 this should also be so bracketed.
// Scope         : Private
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
void EpetraVector::processError(const char *methodMsg, int error) const
{
  // Process the error
  switch (error)
  {
  case 0:
    //Xyce::dout() << methodMsg << ": Function returned without warnings or errors." << std::endl;
    break;

  default:
    Xyce::Report::DevelFatal0().in(methodMsg) << "Function returned with an error.";
  }
}

//-----------------------------------------------------------------------------
// Function      : print
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/14/00
//-----------------------------------------------------------------------------
void EpetraVector::print(std::ostream &os) const
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

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
#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>
#include <N_PDS_Comm.h>
#include <N_PDS_EpetraHelpers.h>
#include <N_PDS_EpetraParMap.h>
#include <N_UTL_FeatureTest.h>

// ---------  Other Includes  -----------

#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

#ifdef Xyce_PARALLEL_MPI
#include <Epetra_MpiComm.h>
#endif

#include <EpetraExt_View_MultiVector.h>
#include <EpetraExt_MultiVectorOut.h>
#include <Teuchos_BLAS.hpp>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : MultiVector::MultiVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
MultiVector::MultiVector(N_PDS_ParMap & map, int numVectors)
:  parallelMap_(&map),
   overlapMap_(&map),
   importer_(0),
   exporter_(0),
   viewTransform_(0),
   pdsComm_(rcp(&map.pdsComm(),false)), 
   isOwned_(true),
   groundNode_(0.0)
{
  if (map.numGlobalEntities() < 0)
  {
    Report::DevelFatal().in("MultiVector::MultiVector")
      << "vector length too short. Vectors must be > 0 in length.";
  }
  else if (numVectors < 1)
  {
     Report::DevelFatal().in("MultiVector::MultiVector")
      << "numVectors < 1";
  }

  // Create a new Petra MultiVector and set the pointer.
  N_PDS_EpetraParMap& e_map = dynamic_cast<N_PDS_EpetraParMap&>( map );
  aMultiVector_ = new Epetra_MultiVector( *dynamic_cast<Epetra_BlockMap*>(e_map.petraMap()), numVectors );

  oMultiVector_ = aMultiVector_;
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::MultiVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/06/02
//-----------------------------------------------------------------------------
MultiVector::MultiVector(
  N_PDS_ParMap &        map,
  N_PDS_ParMap &        ol_map,
  int                   numVectors )
  : parallelMap_(&map),
    overlapMap_(&ol_map),
    importer_(0),
    exporter_(0),
    viewTransform_(0),
    pdsComm_(rcp(&map.pdsComm(),false)), 
    isOwned_(true),
    groundNode_(0.0)
{
  if (map.numGlobalEntities() < 0)
    Report::DevelFatal().in("MultiVector::MultiVector")
      << "vector length too short. Vectors must be > 0 in length.";

  // Create a new Petra MultiVector and set the pointer.
  N_PDS_EpetraParMap& e_map = dynamic_cast<N_PDS_EpetraParMap&>( map );
  N_PDS_EpetraParMap& e_ol_map = dynamic_cast<N_PDS_EpetraParMap&>( ol_map );
  oMultiVector_ = new Epetra_MultiVector(*dynamic_cast<Epetra_BlockMap*>(e_ol_map.petraMap()), numVectors);

  viewTransform_ = new EpetraExt::MultiVector_View(
    *dynamic_cast<Epetra_BlockMap*>(e_ol_map.petraMap()), *dynamic_cast<Epetra_BlockMap*>(e_map.petraMap()));
  aMultiVector_ = &((*viewTransform_)(*oMultiVector_));
  if (map.pdsComm().numProc() > 1)
  {
    exporter_ = new Epetra_Export( *dynamic_cast<Epetra_BlockMap*>(e_ol_map.petraMap()),
                                   *dynamic_cast<Epetra_BlockMap*>(e_map.petraMap()) );
  }

  importer_ = new Epetra_Import( *dynamic_cast<Epetra_BlockMap*>(e_ol_map.petraMap()), 
                                 *dynamic_cast<Epetra_BlockMap*>(e_map.petraMap()) );
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::MultiVector
// Purpose       : Copy Constructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
MultiVector::MultiVector( const MultiVector & right )
: parallelMap_( right.parallelMap_ ),
  overlapMap_( right.overlapMap_ ),
  oMultiVector_( new Epetra_MultiVector( *(right.oMultiVector_) ) ),
  importer_(0),
  exporter_(0),
  viewTransform_(0),
  pdsComm_( right.pdsComm_ ),
  isOwned_(true),
  groundNode_(0.0)
{
  if (right.aMultiVector_ == right.oMultiVector_)
    aMultiVector_ = oMultiVector_;
  else
  {
    N_PDS_EpetraParMap* e_map = dynamic_cast<N_PDS_EpetraParMap*>( parallelMap_ );
    N_PDS_EpetraParMap* e_ol_map = dynamic_cast<N_PDS_EpetraParMap*>( overlapMap_ );

    viewTransform_ = new EpetraExt::MultiVector_View( *dynamic_cast<Epetra_BlockMap*>(e_ol_map->petraMap()),
                                                      *dynamic_cast<Epetra_BlockMap*>(e_map->petraMap()) );
    aMultiVector_ = &((*viewTransform_)( *oMultiVector_ ));
  }

  // Generate new exporter instead of using copy constructor, there is an issue with Epetra_MpiDistributor
  if( right.exporter_ ) 
  {
    N_PDS_EpetraParMap* e_map = dynamic_cast<N_PDS_EpetraParMap*>( parallelMap_ );
    N_PDS_EpetraParMap* e_ol_map = dynamic_cast<N_PDS_EpetraParMap*>( overlapMap_ );

    exporter_ = new Epetra_Export( *dynamic_cast<Epetra_BlockMap*>(e_ol_map->petraMap()),
                                   *dynamic_cast<Epetra_BlockMap*>(e_map->petraMap()) );
  }
  if( right.importer_ ) importer_ = new Epetra_Import( *right.importer_ );
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::MultiVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
MultiVector::MultiVector( Epetra_MultiVector * overlapMV, const Epetra_BlockMap& parMap, bool isOwned )
: parallelMap_(0),
  overlapMap_(0),
  oMultiVector_( overlapMV ),
  importer_(0),
  exporter_(0),
  viewTransform_(0),
  isOwned_(isOwned),
  groundNode_(0.0) 
{
  // Make sure there is anything to communicate before creating a transform, importer, or exporter
  if (overlapMV->MyLength() == parMap.NumMyElements())
    aMultiVector_ = oMultiVector_;
  else
  {
    viewTransform_ = new EpetraExt::MultiVector_View( overlapMV->Map(), parMap );
    aMultiVector_ = &((*viewTransform_)(*oMultiVector_));
    if( parMap.Comm().NumProc() > 1 )
      exporter_ = new Epetra_Export( overlapMV->Map(), parMap );

    importer_ = new Epetra_Import( overlapMV->Map(), parMap );
  }

  Epetra_Comm& ecomm = const_cast<Epetra_Comm &>( overlapMV->Comm() );
  pdsComm_ = Teuchos::rcp( Xyce::Parallel::createPDSComm( &ecomm ) );  
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::MultiVector
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
MultiVector::MultiVector( Epetra_MultiVector * origMV, bool isOwned )
: parallelMap_(0),
  overlapMap_(0),
  aMultiVector_( origMV ),
  oMultiVector_( origMV ),
  importer_(0),
  exporter_(0),
  viewTransform_(0),
  isOwned_(isOwned),
  groundNode_(0.0)
{
  Epetra_Comm& ecomm = const_cast<Epetra_Comm &>( origMV->Comm() );
  pdsComm_ = Teuchos::rcp( Xyce::Parallel::createPDSComm( &ecomm ) );
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::operator=
// Purpose       : assignment
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
MultiVector & MultiVector::operator=( const MultiVector & right )
{
  if (this != &right)
  {
    if( (oMultiVector_->Map().NumGlobalElements() == right.oMultiVector_->Map().NumGlobalElements())
        && (oMultiVector_->Map().NumMyElements() == right.oMultiVector_->Map().NumMyElements()) )
    {
      *oMultiVector_ = *right.oMultiVector_;
    }

    if( (aMultiVector_->Map().NumGlobalElements() == right.aMultiVector_->Map().NumGlobalElements())
        && (aMultiVector_->Map().NumMyElements() == right.aMultiVector_->Map().NumMyElements()) )
    {
      *aMultiVector_ = *right.aMultiVector_;
    }
    else
    {
      if (VERBOSE_LINEAR)
        Report::DevelFatal0() <<"MultiVector being assigned with different Mapping";
    }
  }

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::~MultiVector
// Purpose       : Default destructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
MultiVector::~MultiVector()
{
  delete importer_;
  delete exporter_;
  delete viewTransform_; //destroys of aMultiVector_ as well
  if (isOwned_)
  {
    delete oMultiVector_;
  }
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::globalLength
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
int MultiVector::globalLength() const
{
  return aMultiVector_->GlobalLength();
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::localLength
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
int MultiVector::localLength() const
{
  return aMultiVector_->MyLength();
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::numVectors
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
int MultiVector::numVectors() const
{
  return aMultiVector_->NumVectors();
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::dotProduct
// Purpose       : Returns the dot product of "this" vector and another.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
void MultiVector::dotProduct(const MultiVector & y, std::vector<double>& d) const
{
  int xnum = aMultiVector_->NumVectors();
  int ynum = y.aMultiVector_->NumVectors();
  if (xnum == 1 || ynum == 1)
  {
    for (int i=0; i<xnum; ++i)
    {
      for (int j=0; j<ynum; ++j)
      {
        (*aMultiVector_)(i)->Dot(*(*y.aMultiVector_)(j), &d[i*ynum + j]);
      }
    }
  }
  else
  {    
    // Let Epetra handle this.
    int PetraError = aMultiVector_->Dot(*(y.aMultiVector_), &d[0]);
  
    if (DEBUG_LINEAR)
      processError( "MultiVector::dotProduct - ", PetraError );
  }
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::scale
// Purpose       : Scales a MultiVector by a constant value.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void MultiVector::scale(const double a)
{
  int PetraError = aMultiVector_->Scale(a);

  if (DEBUG_LINEAR)
    processError( "MultiVector::scale - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::scale
// Purpose       : Scales a MultiVector by a constant value, but
//                 puts it into this.
//                 this = a*x
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/16/00
//-----------------------------------------------------------------------------
void MultiVector::scale(const double a, const MultiVector &x)
{
  int PetraError = aMultiVector_->Scale(a, *(x.aMultiVector_));

  if (DEBUG_LINEAR)
    processError( "MultiVector::scale - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::multiply
// Purpose       : Performs element-wise multiplication of two vectors
//                 this = this @ x
//                 where @ represents element-wise multiplication
// Special Notes :
// Scope         : Public
// Creator       : Roger P. Pawlowski, SNL, Parallel Computational Sciences
// Creation Date : 3/24/03
//-----------------------------------------------------------------------------
void MultiVector::multiply(const MultiVector &x)
{
  int PetraError = aMultiVector_->Multiply(1.0, *aMultiVector_,
					   *(x.aMultiVector_), 0.0);

  if (DEBUG_LINEAR)
    processError( "MultiVector::scale - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::daxpy
// Purpose       : Linear combination of two MultiVectors:
//                 this = y + a*x
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void MultiVector::daxpy(const MultiVector &y, const double a,
                              const MultiVector &x)
{
  int PetraError = aMultiVector_->Update(1.0, *(y.aMultiVector_), a,
                                         *(x.aMultiVector_), 0.0);

  if (DEBUG_LINEAR)
    processError( "MultiVector::daxpy - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::linearCombo
// Purpose       : Linear combination of two MultiVectors:
//                 this = a*x + b*y
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/08/01
//-----------------------------------------------------------------------------
void MultiVector::linearCombo(const double a, const MultiVector &x,
                                    const double b, const MultiVector &y)
{
  int PetraError = aMultiVector_->Update(a, *(x.aMultiVector_), b,
                                         *(y.aMultiVector_), 0.0);

  if (DEBUG_LINEAR)
    processError( "MultiVector::linearCombo - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::update
// Purpose       :
// Special Notes : ERK. From the epetra documentation:
//
//                 this = s*this + a*A
//
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 02/04/02
//-----------------------------------------------------------------------------
void MultiVector::update( double a, const MultiVector & A,
                                double s )
{
  aMultiVector_->Update( a, *(A.aMultiVector_), s );
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::update
// Purpose       :
// Special Notes : ERK.  From the epetra documentation:
//
//                 this = s*this + a*A + b*B
//
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 02/04/02
//-----------------------------------------------------------------------------
void MultiVector::update( double a, const MultiVector & A,
                                double b, const MultiVector & B,
                                double s )
{
  aMultiVector_->Update( a, *(A.aMultiVector_),
                         b, *(B.aMultiVector_),
                         s );
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::addVec
// Purpose       : Add multiple of a MultiVector:  this = this + a*y
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void MultiVector::addVec(const double a, const MultiVector &y)
{
  int PetraError = aMultiVector_->Update(a, *(y.aMultiVector_), 1.0);

  if (DEBUG_LINEAR)
    processError( "MultiVector::addVec - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::lpNorm
// Purpose       : Returns lp norms of each vector in MultiVector
// Special Notes : Only p=1 and p=2 implemented now since this is all Petra
//                 supports.
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
int MultiVector::lpNorm(const int p, double * result) const
{
  int PetraError = -1;
  static const char *methodMsg = "MultiVector::lpNorm - ";

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
// Function      : MultiVector::infNorm
// Purpose       : Returns infinity norm of each vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/16/01
//-----------------------------------------------------------------------------
int MultiVector::infNorm(double * result) const
{
  static const char *methodMsg = "MultiVector::infNorm - ";
  int PetraError = aMultiVector_->NormInf(result);

  if (DEBUG_LINEAR) 
    processError(methodMsg, PetraError);

  return PetraError;
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::infNormIndex
// Purpose       : Returns index of the maximum absolute entry in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Heidi K. Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 11/21/11
//-----------------------------------------------------------------------------
int MultiVector::infNormIndex(int * index) const
{
  Teuchos::BLAS<int,double> blas;

  int numProcs = aMultiVector_->Comm().NumProc();
  int numVectors = aMultiVector_->NumVectors();
  int myLength = aMultiVector_->MyLength();
  std::vector<int> indexTemp( numVectors, 0 ), indexTempAll( numVectors*numProcs, 0 );
  std::vector<double> doubleTemp( numVectors, 0.0 ), doubleTempAll( numVectors*numProcs, 0.0 );
  double ** pointers = aMultiVector_->Pointers();

  for (int i=0; i < numVectors; i++)
  {
    // Remember that IAMAX returns 1-based indexing, so subtract 1 to get the actual index.
    int jj = blas.IAMAX(myLength, pointers[i], 1) - 1;
    if (jj>-1)
    {
      indexTemp[i] = aMultiVector_->Map().GID(jj);
      doubleTemp[i] = std::abs(pointers[i][jj]);
    }
  }

  // Use the Epetra communicator to gather all the local maximum values and indices
  int result = aMultiVector_->Comm().GatherAll(&indexTemp[0], &indexTempAll[0], numVectors);
  result += aMultiVector_->Comm().GatherAll(&doubleTemp[0], &doubleTempAll[0], numVectors);

  // Compute the global infNorm and index
  for (int i=0; i < numVectors; i++)
  {
    // Remember that IAMAX returns 1-based indexing, so subtract 1 to get the actual index.
    int ii = blas.IAMAX( numProcs, &doubleTempAll[i], numVectors ) - 1;
    index[i] = indexTempAll[ii*numVectors + i];
  }

  return result;
}


//-----------------------------------------------------------------------------
// Function      : MultiVector::wRMSNorm
// Purpose       : Returns weighted root-mean-square of each vector in
//                 MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/12/00
//-----------------------------------------------------------------------------
int MultiVector::wRMSNorm(const MultiVector & weights, double * result) const
{
  int PetraError = aMultiVector_->NormWeighted( *(weights.aMultiVector_), result );

  if (DEBUG_LINEAR)
    processError( "MultiVector::wRMSNorm - ", PetraError);

  return PetraError;
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::wMaxNorm
// Purpose       : Returns the weighted
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 03/19/01
//-----------------------------------------------------------------------------
int MultiVector::wMaxNorm(const MultiVector & weights, double * result) const
{
  int length  = aMultiVector_->MyLength();
  int numVecs = numVectors();
  double tmpVal = 0.0;

  for (int i = 0;  i < numVecs; ++i)
  {
    double localMax = 0.0;
    if (length)
    { 
      localMax = fabs((*(this))[i][0]) / weights[i][0]; 
      for (int j = 1; j < length; ++j)
      {
        tmpVal = fabs((*(this))[i][j]) / weights[i][j];
        if (tmpVal > localMax)
          localMax = tmpVal;
      }
    } 
    // Determine global maximum.
    aMultiVector_->Comm().MaxAll( &localMax, &(result[i]), 1 );
  }

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::minValue
// Purpose       : Return the minimum value for the multivector
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/16/00
//-----------------------------------------------------------------------------
int MultiVector::minValue(double * result) const
{
  int PetraError = aMultiVector_->MinValue(result);

  return PetraError;
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::maxValue
// Purpose       : Return the maximum value for the multivector
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/16/00
//-----------------------------------------------------------------------------
int MultiVector::maxValue(double * result) const
{
  int PetraError = aMultiVector_->MaxValue(result);

  return PetraError;
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::random
// Purpose       : Generates random numbers drawn from a uniform distribution
//                 on the interval (-1,1) using a multiplicative congruential
//                 generator with modulus 2^31 - 1.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void MultiVector::random()
{
  int PetraError = aMultiVector_->Random();

  if (DEBUG_LINEAR)
    processError( "MultiVector::random - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::putScalar
// Purpose       : Fills MultiVector with the value "scalar".
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void MultiVector::putScalar(const double scalar)
{
  int PetraError = oMultiVector_->PutScalar(scalar);

  groundNode_ = scalar;

  if (DEBUG_LINEAR)
    processError( "MultiVector::putScalar - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::addScalar
// Purpose       : Adds to MultiVector with the value "scalar".
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 01/18/01
//-----------------------------------------------------------------------------
void MultiVector::addScalar(const double scalar)
{
  int length  = aMultiVector_->MyLength();
  int numVecs = numVectors();

  for (int i = 0; i < numVecs; ++i)
    for (int j = 0; j < length; ++j)
      (*aMultiVector_)[i][j] += scalar;
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::absValue
// Purpose       : Abs value of elements of MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 01/18/01
//-----------------------------------------------------------------------------
void MultiVector::absValue(const MultiVector & A)
{
  int PetraError = oMultiVector_->Abs(*(A.oMultiVector_));

  if (DEBUG_LINEAR)
    processError( "MultiVector::absValue - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::reciprocal
// Purpose       : Reciprocal elements of MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/07/01
//-----------------------------------------------------------------------------
void MultiVector::reciprocal(const MultiVector & A)
{
  int PetraError = oMultiVector_->Reciprocal(*(A.oMultiVector_));

  if (DEBUG_LINEAR)
    processError( "MultiVector::reciprocal - ", PetraError);
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::getVectorView
// Purpose       : Const view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
const Vector* MultiVector::getVectorView(int index) const
{
  const Vector* vec = new Vector((*oMultiVector_)(index),
                                  aMultiVector_->Map(),false);
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::getNonConstVectorView
// Purpose       : NonConst view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
Vector* MultiVector::getNonConstVectorView(int index)
{
  Vector* vec = new Vector((*oMultiVector_)(index),
                            aMultiVector_->Map(),false);
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::getVectorView
// Purpose       : Const view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
const Vector* MultiVector::getVectorViewAssembled(int index) const
{
  const Vector* vec = new Vector( new 
                      Epetra_Vector( View, *aMultiVector_, index ), true );
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::getNonConstVectorView
// Purpose       : NonConst view of individual vector in MultiVector
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
Vector* MultiVector::getNonConstVectorViewAssembled(int index)
{
  Vector* vec = new Vector( new 
                Epetra_Vector( View, *aMultiVector_, index ), true );
  return vec;
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::fillComplete
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/29/03
//-----------------------------------------------------------------------------
void MultiVector::fillComplete()
{
  if ( exporter_ )
  {
    aMultiVector_->Export( *oMultiVector_, *exporter_, Add );
  }
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::vectorImport
// Purpose       : Import using Petra_Import object
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/09/01
//-----------------------------------------------------------------------------
bool MultiVector::vectorImport(const MultiVector * vec,
                               Epetra_Import * importer)
{
  aMultiVector_->Import(*(vec->aMultiVector_), *importer, Insert);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::importOverlap
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/06/02
//-----------------------------------------------------------------------------
bool MultiVector::importOverlap()
{
  bool flag = false;

  if( importer_ )
    flag = oMultiVector_->Import( *aMultiVector_, *importer_, Insert );

  return flag;
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::writeToFile
// Purpose       : Dumps out the multivector entries to a file.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/19/00
//-----------------------------------------------------------------------------
void MultiVector::writeToFile( const char * filename, bool useLIDs, bool mmFormat ) const
{
  int numProcs = aMultiVector_->Comm().NumProc();
  int localRank = aMultiVector_->Comm().MyPID();
  int masterRank = 0;

  if (!mmFormat)
  {
    for( int p = 0; p < numProcs; ++p )
    {
      //A barrier inside the loop so each processor waits its turn.
      aMultiVector_->Comm().Barrier();

      if(p == localRank)
      {
        FILE *file = NULL;

        if(masterRank == localRank)
        {
          //This is the master processor, open a new file.
          file = fopen(filename,"w");

          //Write the RDP_MultiVector dimension n into the file.
          fprintf(file,"%d\n",globalLength());
        }
        else
        {
          //This is not the master proc, open file for appending
          file = fopen(filename,"a");
        }

        //Now loop over the local portion of the RDP_MultiVector.
        int length  = localLength();
        int numVecs = numVectors();

        for (int i = 0; i < numVecs; ++i)
          for (int j = 0; j < length; ++j)
          {
            int loc = aMultiVector_->Map().GID(j);
            if( useLIDs ) loc = j;
            fprintf(file,"%d %d %20.13e\n",i,loc,(*aMultiVector_)[i][j]);
          }
        fclose(file);
      }
    }
  }
  else 
  {
    EpetraExt::MultiVectorToMatrixMarketFile( filename, *aMultiVector_ );
  }
}

//-----------------------------------------------------------------------------
// Function      : MultiVector::processError
// Purpose       : Concrete implementation which processes Petra (in this case)
//                 error codes taken from the Petra member function returns.
// Special Notes : Petra specific.  NOTE ALSO - this function is currently
//                 within the "Xyce_DEBUG_LINEAR" ifdef and so any calls to
//                 this should also be so bracketed.
// Scope         : Private
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
void MultiVector::processError(const char *methodMsg, int error) const
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
// Function      : MultiVector::getElementByGlobalIndex
// Purpose       : Get element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/6/00
//-----------------------------------------------------------------------------
const double & MultiVector::getElementByGlobalIndex(
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
// Function      : MultiVector::setElementByGlobalIndex
// Purpose       : Set element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/6/00
//-----------------------------------------------------------------------------
bool MultiVector::setElementByGlobalIndex(const int & global_index,
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
// Function      : MultiVector::sumElementByGlobalIndex
// Purpose       : Set element from vector using global index.
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/7/00
//-----------------------------------------------------------------------------
bool MultiVector::sumElementByGlobalIndex(const int & global_index,
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
void MultiVector::print(std::ostream &os) const
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

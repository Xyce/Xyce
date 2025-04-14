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

//-------------------------------------------------------------------------
//
// Purpose        : Simple migrator utility using Zoltan utilities
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/06/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// --------- Standard Includes ----------

// --------- Xyce Includes --------------

#include <N_ERH_ErrorMgr.h>

#include <N_LAS_MultiVector.h>

#include <N_PDS_GlobalAccessor.h>
#include <N_PDS_Comm.h>
#include <N_UTL_FeatureTest.h>

#include <N_PDS_EpetraHelpers.h> 

#ifdef Xyce_PARALLEL_MPI
 #include <Epetra_MpiComm.h>
 #include <Epetra_MpiDistributor.h>
#endif

#include <Epetra_SerialComm.h>
#include <Epetra_SerialDistributor.h>

using Xyce::DEBUG_PARALLEL;

namespace Xyce {
namespace Parallel {

//-----------------------------------------------------------------------------
// Function      : GlobalAccessor::GlobalAccessor
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 02/06/01
//-----------------------------------------------------------------------------
GlobalAccessor::GlobalAccessor( const Communicator &comm )
  : pdsComm_(comm),
    numReceiveObjs_(0),
    arrayReceiveGIDs_(0),
    arrayReceiveProcs_(0),
    recvBuf_(0),
    recvBufSize_(0),
    numSendObjs_(0),
    arraySendGIDs_(0),
    arraySendProcs_(0),
    sendBuf_(0),
    sendBufSize_(0),
    distributor_(0),
    petraComm_(Xyce::Parallel::getEpetraComm(&comm))
{
}

//-----------------------------------------------------------------------------
// Function      : GlobalAccessor::~GlobalAccessor
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 06/07/00
//-----------------------------------------------------------------------------
GlobalAccessor::~GlobalAccessor()
{
  delete[] arrayReceiveGIDs_;
  delete[] arrayReceiveProcs_;
  delete[] recvBuf_;
  delete[] arraySendGIDs_;
  delete[] arraySendProcs_;
  delete[] sendBuf_;
  delete distributor_;
}

//-----------------------------------------------------------------------------
// Function      : GlobalAccessor::generateMigrationPlan
// Purpose       : generates migration plan (Comm_Obj) using Zoltan
// 		   utilities
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/07/00
//-----------------------------------------------------------------------------
void GlobalAccessor::generateMigrationPlan()
{
  numReceiveObjs_ = externGIDVector_.size();

  std::map<int,int> sortMap;

  // sort externGIDVector_ numerically by processor
  for( int i = 0; i < numReceiveObjs_; ++i )
    sortMap[ externGIDVector_[i].second ]++;

  int count = 0;
  int tmp;
  for( std::map<int,int>::iterator it_iiM = sortMap.begin();
       it_iiM != sortMap.end(); ++it_iiM )
  {
    tmp = count;
    count += it_iiM->second;
    it_iiM->second = tmp;
  }

  if( !arrayReceiveGIDs_ )
  {
    arrayReceiveGIDs_ = new int[ numReceiveObjs_ ];
    arrayReceiveProcs_ = new int[ numReceiveObjs_ ];
  }
  else
  {
    delete [] arraySendGIDs_;
    delete [] arraySendProcs_;

    delete [] sendBuf_;
    delete [] recvBuf_;
  }

  int loc, val1, val2;
  for( int i = 0; i < numReceiveObjs_; ++i )
  {
    val1 = externGIDVector_[i].first;
    val2 = externGIDVector_[i].second;
    loc = sortMap[ val2 ];
    arrayReceiveGIDs_[ loc ] = val1;
    arrayReceiveProcs_[ loc ] = val2;
    sortMap[ val2 ]++;
  }

  if (DEBUG_PARALLEL)
  {
    std::cout << "GlobalAccessor::generateMigrationPlan:" << std::endl;
    std::cout << " setup numRecvObjs: " << numReceiveObjs_ << std::endl;
    std::cout << " setup numSendObjs: " << numSendObjs_ << std::endl;
    
    for( int i = 0; i < numReceiveObjs_; ++i )
      std::cout << "  " << arrayReceiveGIDs_[i] << " " << arrayReceiveProcs_[i]
                << std::endl;
  }

#ifdef Xyce_PARALLEL_MPI

  if( !distributor_ ) 
  {
    if( pdsComm_.isSerial() )
    {
      distributor_ = new Epetra_SerialDistributor(
                    *(dynamic_cast<const Epetra_SerialComm*>(petraComm_)));
    }
    else
    {
      distributor_ = new Epetra_MpiDistributor(
                    *(dynamic_cast<const Epetra_MpiComm*>(petraComm_)));
      // can only call CreateFromRecvs if we're parallel with more than one proc
      distributor_->CreateFromRecvs( numReceiveObjs_, arrayReceiveGIDs_, arrayReceiveProcs_,
	                         true, numSendObjs_, arraySendGIDs_, arraySendProcs_ );
    }
  }

  sendGIDVector_.resize( numSendObjs_ );

  for( int i = 0; i < numSendObjs_; ++i )
    sendGIDVector_[i] = std::pair<int,int>( arraySendGIDs_[i],
			arraySendProcs_[i] );

  sendBufSize_ = numSendObjs_ * (sizeof(int)+sizeof(double));
  recvBufSize_ = numReceiveObjs_ * (sizeof(int)+sizeof(double));
  sendBuf_ = new char[ sendBufSize_ ];
  recvBuf_ = new char[ recvBufSize_ ];

  if (DEBUG_PARALLEL)
  {
    std::cout << "Created Migration Plan: " << std::endl;
    std::cout << " numRecvObjs: " << numReceiveObjs_ << std::endl;
    for( int i = 0; i < numReceiveObjs_; ++i )
      std::cout << "  " << arrayReceiveGIDs_[i] << " " << arrayReceiveProcs_[i]
                << std::endl;
    std::cout << " numSendObjs: " << numSendObjs_ << std::endl;
    for( int i = 0; i < numSendObjs_; ++i )
      std::cout << "  " << arraySendGIDs_[i] << " " << arraySendProcs_[i]
                << std::endl;
  }
#endif /* Xyce_PARALLEL_MPI */
}

//-----------------------------------------------------------------------------
// Function      : GlobalAccessor::migrateMultiVector
// Purpose       : migrates nonlocal parts of multivector based on
//                 migration plan (theZoltanCommObjPtr_)
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/07/00
//-----------------------------------------------------------------------------
void GlobalAccessor::migrateMultiVector( Xyce::Linear::MultiVector * mVector )
{

  if (DEBUG_PARALLEL)
    std::cout << "migrating vector<<<<<<<<<<<<<<" << std::endl
              << " " << numSendObjs_ << " " << numReceiveObjs_ << std::endl;

#ifdef Xyce_PARALLEL_MPI
  int pos = 0;
  int idx;
  double val;
  for( int i = 0; i < numSendObjs_; ++i )
  {
    val = mVector->getElementByGlobalIndex( sendGIDVector_[i].first );

    idx = sendGIDVector_[i].first;
    pdsComm_.pack( &idx, 1, sendBuf_, sendBufSize_, pos );
    pdsComm_.pack( &val, 1, sendBuf_, sendBufSize_, pos );
  }

  distributor_->Do( sendBuf_, sizeof(int)+sizeof(double), recvBufSize_, recvBuf_ );

  mVector->clearExternVectorMap();

  pos = 0;
  for( int i = 0; i < numReceiveObjs_; ++i )
  {
    pdsComm_.unpack( recvBuf_, numReceiveObjs_*(sizeof(int)+sizeof(double)), pos, &idx, 1 );
    pdsComm_.unpack( recvBuf_, numReceiveObjs_*(sizeof(int)+sizeof(double)), pos, &val, 1 );

    mVector->addElementToExternVectorMap( idx, val );
  }
#endif /* Xyce_PARALLEL_MPI */
}

//-----------------------------------------------------------------------------
// Function      : GlobalAccessor::migrateIntArray
// Purpose       : migrates nonlocal parts of integer array based on
//                 migration plan (theZoltanCommObjPtr_)
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/07/00
//-----------------------------------------------------------------------------
void GlobalAccessor::migrateIntArray( std::map<int,int> & sendMap,
	std::map<int,int> & recvMap )
{

#ifdef Xyce_PARALLEL_MPI

  std::map<int,int>::iterator it_i2M, end_i2M;

  if (DEBUG_PARALLEL)
  {
    std::cout << "Send Map" << std::endl;
    it_i2M = sendMap.begin();
    end_i2M = sendMap.end();
    for( ; it_i2M != end_i2M; ++it_i2M )
      std::cout << "  " << it_i2M->first << " " << it_i2M->second << std::endl;
    std::cout << std::endl;
  }

  int pos = 0;
  int idx;
  int val;
  for( int i = 0; i < numSendObjs_; ++i )
  {
    idx = sendGIDVector_[i].first;
    val = sendMap[ sendGIDVector_[i].first ];

    if (DEBUG_PARALLEL)
      std::cout << "send values: " << sendGIDVector_[i].first
                << " " << val << std::endl;

    pdsComm_.pack( &idx, 1, sendBuf_, 2 * numSendObjs_ * sizeof( int ), pos );
    pdsComm_.pack( &val, 1, sendBuf_, 2 * numSendObjs_ * sizeof( int ), pos );
  }

  distributor_->Do( sendBuf_, 2*sizeof(int), recvBufSize_, recvBuf_ );

  recvMap.clear();

  pos = 0;
  for( int i = 0; i < numReceiveObjs_; ++i )
  {
    pdsComm_.unpack( recvBuf_, 2 * numReceiveObjs_ * sizeof( int ), pos, &idx, 1 );
    pdsComm_.unpack( recvBuf_, 2 * numReceiveObjs_ * sizeof( int ), pos, &val, 1 );

    recvMap[ idx ] = val;
  }

  if (DEBUG_PARALLEL)
  {
    std::cout << "Recv Map" << std::endl;
    it_i2M = recvMap.begin();
    end_i2M = recvMap.end();
    for( ; it_i2M != end_i2M; ++it_i2M )
      std::cout << "  " << it_i2M->first << " " << it_i2M->second << std::endl;
    std::cout << std::endl;
  }
#endif /* Xyce_PARALLEL_MPI */
}

//-----------------------------------------------------------------------------
// Function      : GlobalAccessor::migrateIntVecs
// Purpose       : migrates nonlocal parts of integer vectors
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
void GlobalAccessor::migrateIntVecs( std::map< int,std::vector<int> > & sendMap,
	std::map< int,std::vector<int> > & recvMap )
{

#ifdef Xyce_PARALLEL_MPI

  std::map< int,std::vector<int> >::iterator it_i2M, end_i2M;

  int maxSize = 0;
  int maxGlobalSize;

  double tmpVar1, tmpVar2;

  it_i2M = sendMap.begin();
  end_i2M = sendMap.end();
  for( ; it_i2M != end_i2M; ++it_i2M )
    if( maxSize < it_i2M->second.size() ) maxSize = it_i2M->second.size();

  tmpVar1 = maxSize;
  pdsComm_.maxAll( &tmpVar1, &tmpVar2, 1 );
  maxGlobalSize = tmpVar2;

  int pos = 0;
  int idx, val;

  std::map<int,int> recvSizeMap;

  for( int i = 0; i < numSendObjs_; ++i )
  {
    idx = sendGIDVector_[i].first;
    val = sendMap[ idx ].size();

    pdsComm_.pack( &idx, 1, sendBuf_, 2 * numSendObjs_ * sizeof( int ), pos );
    pdsComm_.pack( &val, 1, sendBuf_, 2 * numSendObjs_ * sizeof( int ), pos );
  }

  distributor_->Do( sendBuf_, 2*sizeof(int), recvBufSize_, recvBuf_ );

  recvSizeMap.clear();

  pos = 0;
  for( int i = 0; i < numReceiveObjs_; ++i )
  {
    pdsComm_.unpack( recvBuf_, 2 * numReceiveObjs_ * sizeof( int ), pos, &idx, 1 );
    pdsComm_.unpack( recvBuf_, 2 * numReceiveObjs_ * sizeof( int ), pos, &val, 1 );

    recvSizeMap[ idx ] = val;
    recvMap[ idx ] = std::vector<int>(val,0);
  }

  for( int loc = 0; loc < maxGlobalSize; ++loc )
  {
    pos = 0;

    for( int i = 0; i < numSendObjs_; ++i )
    {
      idx = sendGIDVector_[i].first;
      val = -1;
      if( loc < sendMap[ idx ].size() )
        val = (sendMap[ idx ])[loc];

      pdsComm_.pack( &idx, 1, sendBuf_, 2 * numSendObjs_ * sizeof( int ), pos );
      pdsComm_.pack( &val, 1, sendBuf_, 2 * numSendObjs_ * sizeof( int ), pos );
    }

    distributor_->Do( sendBuf_, 2*sizeof(int), recvBufSize_, recvBuf_ );

    pos = 0;
    for( int i = 0; i < numReceiveObjs_; ++i )
    {
      pdsComm_.unpack( recvBuf_, 2 * numReceiveObjs_ * sizeof( int ), pos, &idx, 1 );
      pdsComm_.unpack( recvBuf_, 2 * numReceiveObjs_ * sizeof( int ), pos, &val, 1 );

      if( loc < recvMap[ idx ].size() ) (recvMap[ idx ])[ loc ] = val;
    }

  }

#endif /* Xyce_PARALLEL_MPI */

} // namespace Parallel
} // namespace Xyce

}


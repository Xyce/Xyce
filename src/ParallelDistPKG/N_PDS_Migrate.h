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

//-----------------------------------------------------------------------------
//
// Purpose        : Migrate tool using Zoltan utilities
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 08/11/03
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_PDS_Migrate_h
#define Xyce_N_PDS_Migrate_h

// ---------- Standard Includes ----------

#include <string>
#include <vector>
#include <map>

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;
// ---------- Forward Declarations -------

// ----------   Xyce Includes   ----------

#include <N_PDS_Comm.h>

#include <N_PDS_PackTraits.h>
#include <N_ERH_ErrorMgr.h>

// ----------   Other Includes   ----------

#ifdef Xyce_PARALLEL_MPI
#include <Epetra_MpiComm.h>
#include <Epetra_MpiDistributor.h>
#endif


namespace Xyce {
namespace Parallel {

//-----------------------------------------------------------------------------
// Class         : Xyce::Parallel::Migrate
// Purpose       : Migrate functor object
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 08/11/03
//-----------------------------------------------------------------------------
template <typename KT, typename DT>
class Migrate
{
 public:

  typedef typename std::map< KT, Teuchos::RCP<DT> >         DataMap;
  typedef typename DataMap::iterator        DataMapIter;
  typedef typename DataMap::const_iterator  DataMapCIter;

  typedef typename DataMap::value_type        DataPair;

  typedef typename std::vector<KT>          KeyList;
  typedef typename KeyList::iterator        KeyListIter;
  typedef typename KeyList::const_iterator  KeyListCIter;

  typedef typename std::vector<int>   ProcList;
  typedef typename ProcList::iterator ProcListIter;

  typedef typename std::vector<char> Buffer;

  // Constructor
  Migrate( Communicator & comm )
  : comm_(comm),
    imports_(0),
    importSize_(0)
  {}

  // Destructor
  ~Migrate()
  { 
    if( importSize_ ) delete [] imports_; 
  }

 private:

  // No public copy construction, assignment, or equality operators.
  Migrate();

  bool operator==( Migrate const & right ) const;
  bool operator!=( Migrate const & right ) const;

 public:

  void operator()( std::vector<int> const & pList,
                   std::vector<KT> const & iKeys,
                   std::vector<KT> & oKeys );
  
  void operator()( std::vector<int> const & pList,
                   std::map< KT, RCP<DT> > const & iData,
                   std::multimap< KT, RCP<DT> > & oData );
  
  void rvs( std::vector<int> const & pList,
            std::vector<KT> const & keys,
            std::map< KT, RCP<DT> > & iData,
            std::map< KT, RCP<DT> > & oData );
  
 protected:

  Communicator & comm_;

  char * imports_;
  int    importSize_;

  Buffer exports_;
};

template <typename DT>
class Migrate1
{
 public:

  typedef typename Teuchos::RCP<DT> DataPtr;
  typedef typename std::vector<DataPtr>   DataContainer;
  typedef typename DataContainer::iterator        DataContainerIter;
  typedef typename DataContainer::const_iterator  DataContainerCIter;

  typedef typename std::vector<int>   ProcList;
  typedef typename ProcList::iterator ProcListIter;

  typedef typename std::vector<char>  Buffer;

  // Constructor
  Migrate1( Communicator & comm )
  : comm_(comm),
    imports_(0),
    importSize_(0)
  {}

  // Destructor
  ~Migrate1()
  { if( importSize_ ) delete [] imports_; }

 private:

  // No public copy construction, assignment, or equality operators.
  Migrate1();

  bool operator==( Migrate1 const & right ) const;
  bool operator!=( Migrate1 const & right ) const;

 public:

  void operator()( std::vector<int> const & pList,
                   std::vector< RCP<DT> > const & iData,
                   std::vector< RCP<DT> > & oData );
  
  void rvs( std::vector<int> const & pList,
            std::vector< RCP<DT> > const & iData,
            std::vector< RCP<DT> > & oData );
  
 protected:

  Communicator & comm_;

  char * imports_;
  int    importSize_;

  Buffer exports_;

};

//-----------------------------------------------------------------------------
// Function      : Xyce::Parallel::Migrate::operator()
// Purpose       : Fwd Migrate Objects 
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/11/03
//-----------------------------------------------------------------------------
template <typename KT, typename DT>
void
Migrate<KT,DT>::
operator()( std::vector<int> const & pList,
            std::vector<KT> const & iKeys,
            std::vector<KT> & oKeys )
{
#ifdef Xyce_PARALLEL_MPI
  if( !comm_.isSerial() )
  {
    Epetra_MpiComm petraComm( comm_.comm() );
    Epetra_MpiDistributor distributor( petraComm );

    int exportCnt = pList.size();

    int max_size = 0;
    KeyListCIter citKL = iKeys.begin();
    KeyListCIter cendKL = iKeys.end();
    for( ; citKL != cendKL; ++citKL )
      max_size = std::max( max_size, PackTraits<KT>::size( *citKL ) );

    int importCnt;
    distributor.CreateFromSends( exportCnt, &(pList[0]), true, importCnt ); 

    double d_max_size = static_cast<double>(max_size);
    double d_max_all;
    comm_.maxAll( &d_max_size, &d_max_all, 1 );
    int max_all = static_cast<int>(d_max_all);

    exports_.resize( max_all * exportCnt );

    if( importSize_ < (max_all*importCnt) )
    {
      if( importSize_ ) delete [] imports_;
      importSize_ = (max_all*importCnt);
      imports_ = new char[importSize_];
    }

    int pos = 0;
    citKL = iKeys.begin();
    for( int i = 0; citKL != cendKL; ++citKL, ++i )
    {
      pos = max_all * i;
      PackTraits<KT>::pack( *citKL, &(exports_[0]), (max_all*exportCnt ), pos, comm_ );
    }

    distributor.Do( &(exports_[0]), max_all, importSize_, imports_ );

    oKeys.resize( importCnt );
    for( int i = 0; i < importCnt; ++i )
    {
      pos = max_all * i;
      PackTraits<KT>::unpack( oKeys[i], &(imports_[0]), (max_all*importCnt), pos, comm_ );
    }
  }
  else
  {
    // Parallel binary running in serial so just copy data
    oKeys = iKeys;
  }
#else
  //Just Copy Data
  oKeys = iKeys;
#endif
}
  
//-----------------------------------------------------------------------------
// Function      : Xyce::Parallel::Migrate::operator()
// Purpose       : Fwd Migrate Objects 
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/11/03
//-----------------------------------------------------------------------------
template <typename KT, typename DT>
void
Migrate<KT,DT>::
operator()( std::vector<int> const & pList,
            std::map< KT, RCP<DT> > const & iData,
            std::multimap< KT, RCP<DT> > & oData )
{
#ifdef Xyce_PARALLEL_MPI
  if( !comm_.isSerial() )
  {
    Epetra_MpiComm petraComm( comm_.comm() );
    Epetra_MpiDistributor distributor( petraComm );

    int exportCnt = pList.size();

    int max_size = 0;
    DataMapCIter citDM  = iData.begin();
    DataMapCIter cendDM = iData.end();
    for( ; citDM != cendDM; ++citDM )
      max_size = std::max( max_size, PackTraits<KT>::size( citDM->first )
                 + PackTraits<DT>::size( *(citDM->second) ) );

    int importCnt;
    distributor.CreateFromSends( exportCnt, &(pList[0]), true, importCnt ); 

    double d_max_size = static_cast<double>(max_size);
    double d_max_all;
    comm_.maxAll( &d_max_size, &d_max_all, 1 );
    int max_all = static_cast<int>(d_max_all);

    exports_.resize( max_all * exportCnt );

    if( importSize_ < (max_all*importCnt) )
    {
      if( importSize_ ) delete [] imports_;
      importSize_ = (max_all*importCnt);
      imports_ = new char[importSize_];
    }

    int pos = 0;
    citDM = iData.begin();
    for( int i = 0; citDM != cendDM; ++citDM, ++i )
    {
      pos = max_all * i;
      PackTraits<KT>::pack( citDM->first, &(exports_[0]), (max_all*exportCnt ), pos, comm_ );
      PackTraits<DT>::pack( *(citDM->second), &(exports_[0]), (max_all*exportCnt ), pos, comm_ );
    }

    distributor.Do( &(exports_[0]), max_all, importSize_, imports_ );

    oData.clear();
    KT key;
    for( int i = 0; i < importCnt; ++i )
    {
      pos = max_all * i;
      PackTraits<KT>::unpack( key, &(imports_[0]), (max_all*importCnt), pos, comm_ );
      RCP<DT> data(new DT);
      PackTraits<DT>::unpack( *data, &(imports_[0]), (max_all*importCnt), pos, comm_ );
      oData.insert( DataPair( key, data ) );
    }
  }
  else
  {
    // paralle binary running in serial so just copy data
    DataMapCIter citDM  = iData.begin();
    DataMapCIter cendDM = iData.end();
    for( ; citDM != cendDM; ++citDM )
      oData.insert( *citDM );
  }
#else
  //Just Copy Data
  DataMapCIter citDM  = iData.begin();
  DataMapCIter cendDM = iData.end();
  for( ; citDM != cendDM; ++citDM )
    oData.insert( *citDM );
#endif
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Parallel::Migrate::rvs
// Purpose       : Rvs Migrate Objects 
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/11/03
//-----------------------------------------------------------------------------
template <typename KT, typename DT>
void
Migrate<KT,DT>::
rvs( std::vector<int> const & pList,
     std::vector<KT> const & keys,
     std::map< KT, RCP<DT> > & iData,
     std::map< KT, RCP<DT> > & oData )
{
#ifdef Xyce_PARALLEL_MPI
  if( !comm_.isSerial() )
  {
    Epetra_MpiComm petraComm( comm_.comm() );
    Epetra_MpiDistributor distributor( petraComm );

    int importCnt = pList.size();
    int exportCnt;

    distributor.CreateFromSends( importCnt, &(pList[0]), true, exportCnt );

    if( exportCnt != (int)keys.size() )
      Xyce::Report::DevelFatal().in("Xyce::Parallel::Migrate::rvs")
        << "Failed Size Match!";

    int max_size = 0;
    KeyListCIter citKL  = keys.begin();
    KeyListCIter cendKL = keys.end();
    for( ; citKL != cendKL; ++citKL )
      max_size = std::max( max_size, PackTraits<KT>::size( *citKL )
                                   + PackTraits<DT>::size( *(iData[*citKL]) ) );

    double d_max_size = static_cast<double>(max_size);
    double d_max_all;
    comm_.maxAll( &d_max_size, &d_max_all, 1 );
    int max_all = static_cast<int>(d_max_all);

    exports_.resize( max_all * exportCnt );

    if( importSize_ < (max_all*importCnt) )
    {
      if( importSize_ ) delete [] imports_;
      importSize_ = (max_all*importCnt);
      imports_ = new char[importSize_];
    }

    int pos = 0;
    int i = 0;
    citKL  = keys.begin();
    for( ; citKL != cendKL; ++citKL, ++i )
    {
      pos = max_all * i;
      PackTraits<KT>::pack( *citKL, &(exports_[0]), (max_all*exportCnt ), pos, comm_ );
      PackTraits<DT>::pack( *(iData[*citKL]), &(exports_[0]), (max_all*exportCnt ), pos, comm_ );
    }

    distributor.DoReverse( &(exports_[0]), max_all, importSize_, imports_ );

    oData.clear();
    KT key;
    for( int i = 0; i < importCnt; ++i )
    {
      pos = max_all * i;
      PackTraits<KT>::unpack( key, &(imports_[0]), (max_all*importCnt), pos, comm_ );
      RCP<DT> data(new DT);
      PackTraits<DT>::unpack( *data, &(imports_[0]), (max_all*importCnt), pos, comm_ );
      oData[key] = data;
    }
  }
  else
  {
    // parallel binary running in serial so just copy data
    oData = iData;
  }
#else
  oData = iData;
#endif
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Parallel::Migrate1::operator()
// Purpose       : Fwd Migrate Objects 
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/11/03
//-----------------------------------------------------------------------------
template <typename DT>
void
Migrate1<DT>::
operator()( std::vector<int> const & pList,
            std::vector< RCP<DT> > const & iData,
            std::vector< RCP<DT> > & oData )
{
#ifdef Xyce_PARALLEL_MPI
  if( !comm_.isSerial() )
  {
    Epetra_MpiComm petraComm( comm_.comm() );
    Epetra_MpiDistributor distributor( petraComm );

    int exportCnt = pList.size();

    int max_size = 0;
    DataContainerCIter citDC = iData.begin();
    DataContainerCIter cendDC = iData.end();
    for( ; citDC != cendDC; ++citDC )
      max_size = std::max( max_size, PackTraits<DT>::size( **citDC ) );

    int importCnt;
    distributor.CreateFromSends( exportCnt, &(pList[0]), true, importCnt ); 

    double d_max_size = static_cast<double>(max_size);
    double d_max_all;
    comm_.maxAll( &d_max_size, &d_max_all, 1 );
    int max_all = static_cast<int>(d_max_all);

    exports_.resize( max_all * exportCnt );

    if( importSize_ < (max_all*importCnt) )
    {
      if( importSize_ ) delete [] imports_;
      importSize_ = (max_all*importCnt);
      imports_ = new char[importSize_];
    }


    int pos = 0;
    citDC = iData.begin();
    for( int i = 0; citDC != cendDC; ++citDC, ++i )
    {
      pos = max_all * i;
      PackTraits<DT>::pack( **citDC, &(exports_[0]), (max_all*exportCnt ), pos, comm_ );
    }

    distributor.Do( &(exports_[0]), max_all, importSize_, imports_ );

    oData.clear();
    for( int i = 0; i < importCnt; ++i )
    {
      pos = max_all * i;
      RCP<DT> data(rcp(new DT));
      PackTraits<DT>::unpack( *data, &(imports_[0]), (max_all*importCnt), pos, comm_ );
      oData.push_back( data );
    }
  }
  else
  {
    // parallel binary running in serial so just copy data
    oData = iData;
  }
#else
  //Just Copy Data
  oData = iData;
#endif
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Parallel::Migrate1::rvs
// Purpose       : Rvs Migrate Objects 
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/11/03
//-----------------------------------------------------------------------------
template <typename DT>
void
Migrate1<DT>::
rvs( std::vector<int> const & pList,
     std::vector< RCP<DT> > const & iData,
     std::vector< RCP<DT> > & oData )
{
#ifdef Xyce_PARALLEL_MPI
  if( !comm_.isSerial() )
  {
    Epetra_MpiComm petraComm( comm_.comm() );
    Epetra_MpiDistributor distributor( petraComm );

    int importCnt = pList.size();
    int exportCnt;

    distributor.CreateFromSends( importCnt, &(pList[0]), true, exportCnt );

    int max_size = 0;
    DataContainerCIter citDC  = iData.begin();
    DataContainerCIter cendDC = iData.end();
    for( ; citDC != cendDC; ++citDC )
      max_size = std::max( max_size, PackTraits<DT>::size( **citDC ) );

    double d_max_size = static_cast<double>(max_size);
    double d_max_all;
    comm_.maxAll( &d_max_size, &d_max_all, 1 );
    int max_all = static_cast<int>(d_max_all);

    exports_.resize( max_all * exportCnt );

    if( importSize_ < (max_all*importCnt) )
    {
      if( importSize_ ) delete [] imports_;
      importSize_ = (max_all*importCnt);
      imports_ = new char[importSize_];
    }

    int i = 0;
    int pos = 0;
    citDC  = iData.begin();
    for( ; citDC != cendDC; ++citDC, ++i )
    {
      pos = max_all * i;
      PackTraits<DT>::pack( **citDC, &(exports_[0]), (max_all*exportCnt ), pos, comm_ );
    }

    distributor.DoReverse( &(exports_[0]), max_all, importSize_, imports_ );

    oData.clear();
    for( int i = 0; i < importCnt; ++i )
    {
      pos = max_all * i;
      RCP<DT> data( rcp(new DT) );
      PackTraits<DT>::unpack( *data, &(imports_[0]), (max_all*importCnt), pos, comm_ );
      oData.push_back( data );
    }
  }
  else
  {
    // parallel binary running in serial so just copy data
    oData = iData;
  }
#else
  oData = iData;
#endif
}

} //namespace Parallel
} //namespace Xyce

#endif // Xyce_N_PDS_Migrate_h

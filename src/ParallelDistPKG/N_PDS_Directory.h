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

//-----------------------------------------------------------------------------
//
// Purpose        : Generic templated distributed directory
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 08/07/03
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_PDS_Directory_h
#define Xyce_N_PDS_Directory_h

// ---------- Standard Includes ----------
#include <map>
#include <iostream>
#include <vector>

#include <N_UTL_Misc.h> // To define NodeID for the topology package

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;

// ---------- Forward Declarations -------

// ----------   Xyce Includes   ----------

#include <N_ERH_ErrorMgr.h>

#include <N_UTL_Functors.h>

// ----------   Other Includes   ----------

namespace Xyce {
namespace Parallel {

//-----------------------------------------------------------------------------
// Class         : N_PDS_Directory
// Purpose       : Distributed directory for object lookup
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/18/01
//-----------------------------------------------------------------------------
template <typename KT, typename DT, class DH, class AC, class MG>
class Directory
{

public:

  typedef typename std::map< KT, RCP<DT> >         DataMap;
  typedef typename DataMap::iterator        DataMapIter;
  typedef typename DataMap::const_iterator  DataMapCIter;

  typedef typename std::multimap< KT, RCP<DT> >    DataRecvMap;
  typedef typename DataRecvMap::iterator        DataRecvMapIter;
  typedef typename DataRecvMap::const_iterator  DataRecvMapCIter;

  typedef typename std::vector<KT>          KeyList;
  typedef typename KeyList::iterator        KeyListIter;
  typedef typename KeyList::const_iterator  KeyListCIter;

  typedef typename std::vector<int>   ProcList;
  typedef typename ProcList::iterator ProcListIter;

  typedef typename AC::iterator       ContainerIter;
  typedef typename AC::const_iterator ContainerCIter;

  // Constructors
  Directory( MG migrate,
	     DH distHash )
  : migrate_(migrate),
    distHash_(distHash)
  {}

  // Destructor
  ~Directory() {}

private:
  // No public copy construction, assignment, or equality operators
  Directory( const Directory & );

  Directory & operator=( const Directory & );

  bool operator==( const Directory & ) const;
  bool operator!=( const Directory & ) const;

public:

  // Add objects from directory.
  void addEntries( DataMap const & entries );

  // Remove objects from directory.
  void deleteEntries( KeyList & keys );

  // Get the items in the directory.
  void getEntries( KeyList & keys,
                   DataMap & entries );

  AC & container() { return container_; }
  ContainerIter & begin() { return container_.begin(); }
  ContainerIter & end() { return container_.end(); }

protected:

  void pushKeys_( KeyList &, KeyList &, ProcList & );
  void pushData_( DataMap const &, DataRecvMap &, ProcList & );

  MG migrate_;
  DH distHash_;
  AC container_;

};

//-----------------------------------------------------------------------------
// Function      : Xyce::Parallel::Hash
// Purpose       : Basic hash class with impl. for std::string
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/20/03
//-----------------------------------------------------------------------------
template <typename T>
class Hash
{
 public:
  int operator()( const T & in ) { return 0; }
  int size();
};

template <>
class Hash<std::string>
{
  int size_;

 public:

  Hash( int size )
  : size_( size )
  {}

  int operator()( const std::string & in )
  {
    int slen = in.length();
    int sum = 0;
    for( int i = 0; i < slen; ++i )
      sum += static_cast<int>( in[i] ); 

    return static_cast<int>( fmod( static_cast<double>( sum ), static_cast<double>(size_) ) );
  }

  int size() {return size_;}
};

template <>
class Hash<NodeID>
{
  int size_;

 public:

  Hash( int size )
  : size_( size )
  {}

  // Perform a hash on a pair<string,int>
  int operator()( const NodeID & in )
  {
    int slen = in.first.length();
    int sum = 0;
    for( int i = 0; i < slen; ++i )
      sum += static_cast<int>( (in.first)[i] ); 
    sum += in.second;

    return static_cast<int>( fmod( static_cast<double>( sum ), static_cast<double>(size_) ) );
  }

  int size() {return size_;}
};


//-----------------------------------------------------------------------------
// Function      : Xyce::Parallel::Directory::addEntries
// Purpose       : Add entries to directory
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/7/03
//-----------------------------------------------------------------------------
template <typename KT, typename DT, class DH, class AC, class MG>
void
Directory<KT,DT,DH,AC,MG>::
addEntries( DataMap const & entries )
{

  DataRecvMap newEntries;
  ProcList procs;
  pushData_( entries, newEntries, procs );

  DataRecvMapCIter citDM = newEntries.begin();
  DataRecvMapCIter cendDM = newEntries.end();

  for( ; citDM != cendDM; ++citDM )
      container_.insert( *citDM );
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Parallel::Directory::deleteEntries
// Purpose       : Delete entries from directory
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/7/03
//-----------------------------------------------------------------------------
template <typename KT, typename DT, class DH, class AC, class MG>
void
Directory<KT,DT,DH,AC,MG>::
deleteEntries( KeyList & keys )
{
  KeyList newKeys;
  ProcList procs;
  pushKeys_( keys, newKeys, procs );

  KeyListCIter citKL = newKeys.begin();
  KeyListCIter cendKL = newKeys.end();

  for( ; citKL != cendKL; ++citKL )
    container_.erase( *citKL );
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Parallel::Directory::getEntries
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/7/03
//-----------------------------------------------------------------------------
template <typename KT, typename DT, class DH, class AC, class MG>
void
Directory<KT,DT,DH,AC,MG>::
getEntries( KeyList & keys,
            DataMap & entries )
{
  //Push Keys to owning processors
  KeyList newKeys;
  ProcList procs;
  pushKeys_( keys, newKeys, procs );
  //int numProcs = procs.size();

  KeyListCIter citKL  = newKeys.begin();
  KeyListCIter cendKL = newKeys.end();

  //Rvs migrate to move data from directory back to requesting procs
  DataMap newEntries;
  for( ; citKL != cendKL; ++citKL )
  {
    ContainerCIter cit = container_.find( *citKL );
    if( cit == container_.end() )
    {
      Xyce::Report::DevelFatal().in("Xyce::Parallel::Director::getEntries")
        << "Data not in directory: " << *citKL;
    }
    newEntries[*citKL] = cit->second;
  }

  migrate_.rvs( procs, newKeys, newEntries, entries );
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Parallel::Directory::pushKeys_
// Purpose       : Pushes keys to owning processors based on hash function
// Special Notes : Forces keys to be ordered by proc number
// Scope         : Protected
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/7/03
//-----------------------------------------------------------------------------
template <typename KT, typename DT, class DH, class AC, class MG>
void
Directory<KT,DT,DH,AC,MG>::
pushKeys_( KeyList & sKeys,
           KeyList & rKeys,
	   ProcList & procs )
{

  KeyListCIter itKL  = sKeys.begin();
  KeyListCIter endKL = sKeys.end();

  procs.clear();
  for( ; itKL != endKL; ++itKL )
    procs.push_back( distHash_(*itKL) );

  if( !IsSorted( procs ) ) SortContainer2( procs, sKeys );
  
  migrate_( procs, sKeys, rKeys );
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Parallel::Directory::pushData_
// Purpose       : Pushes data to owning processors based on hash function
// Special Notes :
// Scope         : Protected
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/7/03
//-----------------------------------------------------------------------------
template <typename KT, typename DT, class DH, class AC, class MG>
void
Directory<KT,DT,DH,AC,MG>::
pushData_( DataMap const & sData,
           DataRecvMap & rData,
	   ProcList & procs )
{
  DataMapCIter itDM  = sData.begin();
  DataMapCIter endDM = sData.end();

  procs.clear();
  for( ; itDM != endDM; ++itDM )
    procs.push_back( distHash_(itDM->first) );

  migrate_( procs, sData, rData );
}

} //namespace Parallel
} //namespace Xyce

#endif

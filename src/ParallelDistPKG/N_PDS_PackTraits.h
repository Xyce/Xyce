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
// Purpose        : 
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

#ifndef Xyce_N_PDS_PackTraits_h
#define Xyce_N_PDS_PackTraits_h

#include <string>
#include <vector>

#include <N_PDS_Comm.h>
#include <N_UTL_Misc.h>
#include <N_UTL_Pack.h>

#include <N_TOP_Node.h>
#include <N_TOP_ParNode.h>

namespace Xyce {
namespace Parallel {

//-----------------------------------------------------------------------------
// Class         : Xyce::Parallel::PackTraits
// Purpose       : 
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 08/11/03
//-----------------------------------------------------------------------------
template <typename T>
struct PackTraits
{
  static int size( T const & object )
  { return Pack<T>::packedByteCount(object); }

  static void pack( T const & object, char * buf, int size, int & pos, Communicator & comm )
  { Pack<T>::pack(object, buf, size, pos, &comm ); }

  static void unpack( T & object, char * buf, int size, int & pos, Communicator & comm )
  { Pack<T>::unpack(object, buf, size, pos, &comm ); }
};

template <>
struct PackTraits<std::string>
{
  static int size( std::string const & object )
  { return object.length() + sizeof(int); }

  static void pack( std::string const & object, char * buf, int size, int & pos, Communicator & comm )
  {
    int len = object.length();
    comm.pack( &len, 1, buf, size, pos );
    comm.pack( object.c_str(), len, buf, size, pos );
  }

  static void unpack( std::string & object, char * buf, int size, int & pos, Communicator & comm )
  {
    int len = 0;
    comm.unpack( buf, size, pos, &len, 1 );
    object = std::string( buf+pos, len );
    pos += len;
  }
};

template <>
struct PackTraits< std::vector<int> >
{
  static int size( std::vector<int> const & object )
  { return (object.size()+1) * sizeof(int); }

  static void pack( std::vector<int> const & object, char * buf, int size, int & pos, Communicator & comm )
  {
    int len = object.size();
    comm.pack( &len, 1, buf, size, pos );
    for( int i = 0; i < len; ++i )
      comm.pack( &object[i], 1, buf, size, pos );
  }

  static void unpack( std::vector<int> & object, char * buf, int size, int & pos, Communicator & comm )
  {
    int len = 0;
    comm.unpack( buf, size, pos, &len, 1 );
    object.resize(len);
    for( int i = 0; i < len; ++i )
      comm.unpack( buf, size, pos, &object[i], 1 );
  }
};

template <>
struct PackTraits< NodeID >
{
  static int size( NodeID const & object )
  { return (object.first.length() + 2*sizeof(int)); }

  static void pack( NodeID const & object, char * buf, int size, int & pos, Communicator & comm )
  {
    int len = object.first.length();
    comm.pack( &len, 1, buf, size, pos );
    comm.pack( object.first.c_str(), len, buf, size, pos );
    comm.pack( &object.second, 1, buf, size, pos );
  }

  static void unpack( NodeID & object, char * buf, int size, int & pos, Communicator & comm )
  {
    int len = 0;
    comm.unpack( buf, size, pos, &len, 1 );
    object.first = std::string( buf+pos, len );
    pos += len;
    comm.unpack( buf, size, pos, &object.second, 1 );
  }
};

} //namespace Parallel
} //namespace Xyce

#endif

//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose        : Node data for parallel applications
//
// Special Notes  : 
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 08/20/03
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_PDS_Node_h
#define Xyce_N_PDS_Node_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <Xyce_config.h>
#include <N_PDS_PackTraits.h>

// ----------  Other Includes   ----------

// ----------  Fwd Declarations ----------

namespace Xyce {
namespace Parallel {

//-----------------------------------------------------------------------------
// Class         : Xyce::Parallel::IndexNode
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 08/20/03
//-----------------------------------------------------------------------------
struct IndexNode
{
  IndexNode()
  : gid(-99),
    pid(-99)
  {}

  IndexNode( int g, int p )
  : gid(g),
    pid(p)
  {}

  int gid;
  int pid;

};

std::ostream & operator<<( std::ostream & os, IndexNode const & in );

template <>
struct PackTraits<IndexNode>
{
  static int size( IndexNode const & object )
  { return 2*sizeof(int); }

  static void pack( IndexNode const & object, char * buf, int size, int & pos, Comm & comm )
  {
    comm.pack( &object.gid, 1, buf, size, pos );
    comm.pack( &object.pid, 1, buf, size, pos );
  }

  static void unpack( IndexNode & object, char * buf, int size, int & pos, Comm & comm )
  {
    comm.unpack( buf, size, pos, &object.gid, 1 );
    comm.unpack( buf, size, pos, &object.pid, 1 );
  }
};

} //namespace Parallel
} //namespace Xyce 

#endif

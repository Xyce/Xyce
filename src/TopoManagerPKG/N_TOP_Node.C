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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/11/01
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_TOP_fwd.h>
#include <N_TOP_Node.h>
#include <N_PDS_Comm.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/11/01
//-----------------------------------------------------------------------------
std::ostream & operator<<( std::ostream & os, Node & node )
{
  return node.put(os);
}


//-----------------------------------------------------------------------------
// Function      : Node::put
// Purpose       :
// Special Notes :
// Scope         : protected
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/11/01
//-----------------------------------------------------------------------------
std::ostream & Node::put(std::ostream & os) const
{
  os << "NodeID:\t" << Xyce::get_node_id(nodeID_) << "\t" << Xyce::get_node_type(nodeID_);
  if (owned_) os << "\tOWNED";
  return os << std::endl;
}

} // namespace Topo

//-----------------------------------------------------------------------------
// Function      : Node::packedByteCount
// Purpose       : Constructor
// Scope         : public
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/01
//-----------------------------------------------------------------------------
template<>
int Pack<Topo::Node>::packedByteCount(
  const Topo::Node &    node)
{
  int count = 0;

  // Pack NodeID
  count += Xyce::packedByteCount( node.nodeID_ );

  count += sizeof(int); //owned

  return count;
}

//-----------------------------------------------------------------------------
// Function      : Node::pack
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/11/01
//-----------------------------------------------------------------------------
template<>
void
Pack<Topo::Node>::pack(
  const Topo::Node &       node,
  char *                   buf,
  int                      bsize,
  int &                    pos,
  Parallel::Communicator * comm ) 
{
  // Pack NodeID
  Xyce::pack( node.nodeID_, buf, bsize, pos, comm );

  int owned = (node.owned_?1:0);
  comm->pack( &owned, 1, buf, bsize, pos );
}

//-----------------------------------------------------------------------------
// Function      : Node::unpack
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/11/01
//-----------------------------------------------------------------------------
template<>
void
Pack<Topo::Node>::unpack(
  Topo::Node &              node,
  char *                    buf,
  int                       bsize,
  int &                     pos,
  Parallel::Communicator *  comm )
{
  // Unpack NodeID
  Xyce::unpack( node.nodeID_, buf, bsize, pos, comm );

  int owned;
  comm->unpack( buf, bsize, pos, &owned, 1 );
  node.owned_ = (owned!=0);
}

} // namespace Xyce

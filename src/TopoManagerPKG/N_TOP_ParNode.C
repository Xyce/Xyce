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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/17/01
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_TOP_fwd.h>
#include <N_TOP_ParNode.h>
#include <N_PDS_Comm.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/17/01
//-----------------------------------------------------------------------------
std::ostream & operator<<( std::ostream & os, const ParNode & pn )
{
  pn.put(os);
  return os;
}

//-----------------------------------------------------------------------------
// Function      : ParNode::put
// Purpose       : Constructor
// Scope         : public
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/01
//-----------------------------------------------------------------------------
std::ostream & ParNode::put(std::ostream & os) const
{
  os << Xyce::subsection_divider << std::endl;
  os << "PARALLEL Node" << std::endl;
  Node::put(os);
  os << "Proc Owner:\t" << proc_ << std::endl;
  os << Xyce::subsection_divider << std::endl << std::endl;
  return os;
}

} // namespace Topo

//-----------------------------------------------------------------------------
// Function      : ParNode::packedByteCount
// Purpose       : Constructor
// Scope         : public
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/01
//-----------------------------------------------------------------------------
template<>
int
Pack<Topo::ParNode>::packedByteCount(
  const Topo::ParNode &         par_node)
{
  int count = Xyce::packedByteCount(static_cast<const Topo::Node &>(par_node));

  count += sizeof(int); //proc

  return count;
}

//-----------------------------------------------------------------------------
// Function      : ParNode::pack
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/17/01
//-----------------------------------------------------------------------------
template<>
void
Pack<Topo::ParNode>::pack(
  const Topo::ParNode &         par_node,
  char *                        buf,
  int                           bsize,
  int &                         pos,
  Parallel::Communicator *      comm ) 
{
  Xyce::pack(static_cast<const Topo::Node &>(par_node), buf, bsize, pos, comm );

  comm->pack( &par_node.proc_, 1, buf, bsize, pos );
}

//-----------------------------------------------------------------------------
// Function      : ParNode::unpack
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/17/01
//-----------------------------------------------------------------------------
template<>
void
Pack<Topo::ParNode>::unpack(
  Topo::ParNode &          par_node,
  char *                   buf,
  int                      bsize,
  int &                    pos,
  Parallel::Communicator * comm )
{
  Xyce::unpack(static_cast<Topo::Node &>(par_node), buf, bsize, pos, comm );

  comm->unpack( buf, bsize, pos, &par_node.proc_, 1 );
}

} // namespace Xyce

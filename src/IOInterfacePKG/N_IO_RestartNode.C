//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Creation Date  : 7/19/01
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <iomanip>
#include <iterator>

#include <N_IO_RestartNode.h>

#include <N_DEV_DeviceState.h>

#include <N_PDS_Comm.h>

#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : RestartNode::RestartNode
// Purpose       : Copy
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/22/01
//-----------------------------------------------------------------------------
RestartNode::RestartNode(const RestartNode & right)
: ID(right.ID),
  type(right.type),
  solnVarData(right.solnVarData),
  stateVarData(right.stateVarData),
  storeVarData(right.storeVarData),
  devState(0)
{
  if (right.devState)
    devState = new Device::DeviceState(*right.devState);
}

//-----------------------------------------------------------------------------
// Function      : RestartNode::operator=
// Purpose       : Assign
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/22/01
//-----------------------------------------------------------------------------
RestartNode & RestartNode::operator=(const RestartNode & right)
{
  if (this != &right)
  {
    ID = right.ID;
    type = right.type;

    solnVarData = right.solnVarData;
    stateVarData = right.stateVarData;
    storeVarData = right.storeVarData;
    delete devState;
    if (right.devState)
      devState = new Device::DeviceState(*right.devState);
  }

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : RestartNode::~RestartNode
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/22/01
//-----------------------------------------------------------------------------
RestartNode::~RestartNode()
{
  delete devState;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 8/22/01
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const RestartNode & rn)
{
  os << Xyce::subsection_divider << std::endl
     << "Restart Node: " << rn.ID << " ( type = " << rn.type << " )" << std::endl;
  std::ostream_iterator<double> out( os, " " );
  if( !rn.solnVarData.empty() )
  {
    os << " SolnVarData: " << std::endl;
    for( unsigned int i = 0; i < rn.solnVarData.size(); ++i )
    {
      os << " " << i << " ";
      copy( rn.solnVarData[i].begin(), rn.solnVarData[i].end(), out );
      os << std::endl;
    }
    os << std::endl;
  }
  if( !rn.stateVarData.empty() )
  {
    os << " StateVarData: " << std::endl;
    for( unsigned int i = 0; i < rn.stateVarData.size(); ++i )
    {
      os << " " << i << " ";
      copy( rn.stateVarData[i].begin(), rn.stateVarData[i].end(), out );
      os << std::endl;
    }
    os << std::endl;
  }

  if( !rn.storeVarData.empty() )
  {
    os << " StoreVarData: " << std::endl;
    for( unsigned int i = 0; i < rn.storeVarData.size(); ++i )
    {
      os << " " << i << " ";
      copy( rn.storeVarData[i].begin(), rn.storeVarData[i].end(), out );
      os << std::endl;
    }
    os << std::endl;
  }

  if( rn.devState ) os << *(rn.devState) << std::endl;

  os << Xyce::subsection_divider << std::endl;
  return os;
}

//-----------------------------------------------------------------------------
// Function      : RestartNode::dump
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 7/23/03
//-----------------------------------------------------------------------------
void RestartNode::dump( std::ostream & os ) const
{
  os << ID << " ";
  os << type << " ";

  int dsize = solnVarData.size();
  os << dsize << " ";
  for( int i = 0; i < dsize; ++i )
  {
    int nsize = solnVarData[i].size();
    os << nsize << " ";
    for( int j = 0; j < nsize; ++j )
      os << std::scientific << std::setw(24) << std::setprecision(17) << solnVarData[i][j] << " ";
  }

  dsize = stateVarData.size();
  os << dsize << " ";
  for( int i = 0; i < dsize; ++i )
  {
    int nsize = stateVarData[i].size();
    os << nsize << " ";
    for( int j = 0; j < nsize; ++j )
      os << std::scientific << std::setw(24) << std::setprecision(17) << stateVarData[i][j] << " ";
  }

  dsize = storeVarData.size();
  os << dsize << " ";
  for( int i = 0; i < dsize; ++i )
  {
    int nsize = storeVarData[i].size();
    os << nsize << " ";
    for( int j = 0; j < nsize; ++j )
      os << std::scientific << std::setw(24) << std::setprecision(17) << storeVarData[i][j] << " ";
  }

  int flag = 1;
  if( devState )
  {
    os << flag << " ";
    devState->dump( os );
  }
  else
  {
    flag = 0;
    os << flag << " ";
  }
}

//-----------------------------------------------------------------------------
// Function      : RestartNode::restore
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 7/23/03
//-----------------------------------------------------------------------------
void RestartNode::restore( std::istream & is )
{
  is >> ID;
  is >> type;

  int dsize;
  is >> dsize;
  solnVarData.resize(dsize);
  for( int i = 0; i < dsize; ++i )
  {
    int nsize;
    is >> nsize;
    solnVarData[i].resize(nsize);
    for( int j = 0; j < nsize; ++j )
      is >> solnVarData[i][j];
  }

  is >> dsize;
  stateVarData.resize(dsize);
  for( int i = 0; i < dsize; ++i )
  {
    int nsize;
    is >> nsize;
    stateVarData[i].resize(nsize);
    for( int j = 0; j < nsize; ++j )
      is >> stateVarData[i][j];
  }

  is >> dsize;
  storeVarData.resize(dsize);
  for( int i = 0; i < dsize; ++i )
  {
    int nsize;
    is >> nsize;
    storeVarData[i].resize(nsize);
    for( int j = 0; j < nsize; ++j )
      is >> storeVarData[i][j];
  }

  int flag;
  is >> flag;
  if( flag == 1 )
  {
    devState = new Device::DeviceState();
    devState->restore(is);
  }
}

} // namespace IO

//-----------------------------------------------------------------------------
// Function      : RestartNode::packedByteCount
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/22/01
//-----------------------------------------------------------------------------
template<>
int
Pack<IO::RestartNode>::packedByteCount(
  const IO::RestartNode &       restart_node)
{
  int byteCount = restart_node.ID.length() + 2*sizeof(int); // ID + length + type

  int svdSize = restart_node.solnVarData.size();

  byteCount += sizeof(int) * (1 + svdSize);
  for (int i = 0; i < svdSize; ++i)
    byteCount += restart_node.solnVarData[i].size() * sizeof(double);

  svdSize = restart_node.stateVarData.size();
  byteCount += sizeof(int) * (1 + svdSize);
  for (int i = 0; i < svdSize; ++i)
    byteCount += restart_node.stateVarData[i].size() * sizeof(double);

  svdSize = restart_node.storeVarData.size();
  byteCount += sizeof(int) * (1 + svdSize);
  for (int i = 0; i < svdSize; ++i)
    byteCount += restart_node.storeVarData[i].size() * sizeof(double);

  byteCount += sizeof(int); //devState flag
  if (restart_node.devState) byteCount += Xyce::packedByteCount(*restart_node.devState);

  return byteCount;
}

//-----------------------------------------------------------------------------
// Function      : RestartNode::pack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 8/22/01
//-----------------------------------------------------------------------------
template<>
void Pack<IO::RestartNode>::pack(
  const IO::RestartNode &       restart_node,
  char *                        buf,
  int                           bsize,
  int &                         pos,
  Parallel::Communicator *      comm)
{
  int size, size2;
  int predictedPos = pos+packedByteCount(restart_node);

  //pack tag
  size = restart_node.ID.length();
  comm->pack( &size, 1, buf, bsize, pos );
  comm->pack( restart_node.ID.c_str(), size, buf, bsize, pos );
  comm->pack( &restart_node.type, 1, buf, bsize, pos );

  size = restart_node.solnVarData.size();
  comm->pack( &size, 1, buf, bsize, pos );

  for( int i = 0; i < size; ++i )
  {
    size2 = restart_node.solnVarData[i].size();
    comm->pack( &size2, 1, buf, bsize, pos );
    for( int ii = 0; ii < size2; ++ii )
      comm->pack( &(restart_node.solnVarData[i][ii]), 1, buf, bsize, pos );
  }

  size = restart_node.stateVarData.size();
  comm->pack( &size, 1, buf, bsize, pos );

  for( int i = 0; i < size; ++i )
  {
    size2 = restart_node.stateVarData[i].size();
    comm->pack( &size2, 1, buf, bsize, pos );
    for( int ii = 0; ii < size2; ++ii )
      comm->pack( &(restart_node.stateVarData[i][ii]), 1, buf, bsize, pos );
  }

  size = restart_node.storeVarData.size();
  comm->pack( &size, 1, buf, bsize, pos );

  for( int i = 0; i < size; ++i )
  {
    size2 = restart_node.storeVarData[i].size();
    comm->pack( &size2, 1, buf, bsize, pos );
    for( int ii = 0; ii < size2; ++ii )
      comm->pack( &(restart_node.storeVarData[i][ii]), 1, buf, bsize, pos );
  }

  int flag = restart_node.devState?1:0;
  comm->pack( &flag, 1, buf, bsize, pos );
  if( restart_node.devState ) Xyce::pack(*restart_node.devState, buf, bsize, pos, comm );
  if (pos != predictedPos)
  {
    std::string msg("Predicted pos does not match actual pos in RestartNode::pack");
    Report::UserFatal() << msg;
  }
}

//-----------------------------------------------------------------------------
// Function      : RestartNode::unpack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 8/22/01
//-----------------------------------------------------------------------------
template<>
void
Pack<IO::RestartNode>::unpack(
  IO::RestartNode &        restart_node,
  char *                   buf,
  int                      bsize,
  int &                    pos,
  Parallel::Communicator * comm)
{
  int size, size2;
  comm->unpack( buf, bsize, pos, &size, 1 );
  restart_node.ID = std::string( (buf+pos), size);
  pos += size;
  comm->unpack( buf, bsize, pos, &restart_node.type, 1 );

  comm->unpack( buf, bsize, pos, &size, 1 );
  restart_node.solnVarData.resize(size);
  for( int i = 0; i < size; ++i )
  {
    comm->unpack( buf, bsize, pos, &size2, 1 );
    restart_node.solnVarData[i].resize(size2);
    for( int ii = 0; ii < size2; ++ii )
      comm->unpack( buf, bsize, pos, &(restart_node.solnVarData[i][ii]), 1 );
  }

  comm->unpack( buf,  bsize, pos, &size, 1 );
  restart_node.stateVarData.resize(size);
  for( int i = 0; i < size; ++i )
  {
    comm->unpack( buf, bsize, pos, &size2, 1 );
    restart_node.stateVarData[i].resize(size2);
    for( int ii = 0; ii < size2; ++ii )
      comm->unpack( buf, bsize, pos, &(restart_node.stateVarData[i][ii]), 1 );
  }

  comm->unpack( buf,  bsize, pos, &size, 1 );
  restart_node.storeVarData.resize(size);
  for( int i = 0; i < size; ++i )
  {
    comm->unpack( buf, bsize, pos, &size2, 1 );
    restart_node.storeVarData[i].resize(size2);
    for( int ii = 0; ii < size2; ++ii )
      comm->unpack( buf, bsize, pos, &(restart_node.storeVarData[i][ii]), 1 );
  }

  int flag;
  comm->unpack( buf, bsize, pos, &flag, 1 );
  if( flag )
  {
    restart_node.devState = new Device::DeviceState();
    Xyce::unpack(*restart_node.devState, buf, bsize, pos, comm );
  }
}

} // namespace Xyce

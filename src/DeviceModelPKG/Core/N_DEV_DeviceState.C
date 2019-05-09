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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 09/02/01
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


#include <iostream>
#include <iomanip>

#include <N_DEV_DeviceState.h>

#include <N_PDS_Comm.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : DeviceState::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
std::ostream & operator<<( std::ostream & os, const DeviceState & ds )
{
  os << "Device State: " << ds.ID << std::endl;
  os << " -------------" << std::endl;
  for( int i = 0; i < ds.data.size(); ++i )
    os << " " << i << ": " << ds.data[i] << std::endl;
  os << " -------------" << std::endl;
  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : DeviceState::dump
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 11/03/04
//-----------------------------------------------------------------------------
void DeviceState::dump( std::ostream & os )
{
  os << ID << " ";

  int size = data.size();
  os << size << " ";
  for( int i = 0; i < size; ++i )
    os << std::scientific << std::setw(24) << std::setprecision(17) << data[i] << " ";


  size = dataInt.size();
  os << size << " ";
  for( int i = 0; i < size; ++i )
    os << dataInt[i] << " ";

}

//-----------------------------------------------------------------------------
// Function      : DeviceState::restore
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 11/03/04
//-----------------------------------------------------------------------------
void DeviceState::restore( std::istream & is )
{
  is >> ID;

  int size;
  is >> size;
  data.resize(size);
  for( int i = 0; i < size; ++i )
    is >> data[i];

  is >> size;
  dataInt.resize(size);
  for( int i = 0; i < size; ++i )
    is >> dataInt[i];

}

} // namespace Device

//-----------------------------------------------------------------------------
// Function      : DeviceState::packedByteCount
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoektra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
template<>
int
Pack<Device::DeviceState>::packedByteCount(const Device::DeviceState &device_state)
{
  int bCnt = sizeof(int);  //ID length
  bCnt += device_state.ID.length();

  bCnt += sizeof(int);  //data double length
  bCnt += device_state.data.size() * sizeof(double);

  bCnt += sizeof(int);  //data int length
  bCnt += device_state.dataInt.size() * sizeof(int);

  return bCnt;
}

//-----------------------------------------------------------------------------
// Function      : DeviceState::pack
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
template<>
void
Pack<Device::DeviceState>::pack(const Device::DeviceState &device_state, char * buf, int bsize, int & pos, N_PDS_Comm * comm)
{
  int length;

  //----- pack ID
  length = device_state.ID.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( device_state.ID.c_str(), length, buf, bsize, pos );

  //----- pack double data
  length = device_state.data.size();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( &(device_state.data[0]), length, buf, bsize, pos );

  //----- pack int data
  length = device_state.dataInt.size();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( &(device_state.dataInt[0]), length, buf, bsize, pos );
}

//-----------------------------------------------------------------------------
// Function      : DeviceState::unpack
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
template<>
void
Pack<Device::DeviceState>::unpack(Device::DeviceState &device_state, char * buf, int bsize, int & pos, N_PDS_Comm * comm)
{
  int length;

  //----- unpack ID
  comm->unpack( buf, bsize, pos, &length, 1 );
  device_state.ID = std::string( (buf+pos), length);
  pos += length;

  //----- unpack data
  comm->unpack( buf, bsize, pos, &length, 1 );
  device_state.data.resize(length);
  comm->unpack( buf, bsize, pos, &(device_state.data[0]), length );

  //----- unpack int data
  comm->unpack( buf, bsize, pos, &length, 1 );
  device_state.dataInt.resize(length);
  comm->unpack( buf, bsize, pos, &(device_state.dataInt[0]), length );
}

} // namespace Xyce

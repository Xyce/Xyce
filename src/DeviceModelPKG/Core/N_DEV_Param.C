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
// Purpose       :
//
// Special Notes :
//
// Creator       : Robert Hoekstra, SNL
//
// Creation Date : 5/15/01
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_DEV_Param.h>
#include <N_PDS_Comm.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {

//-----------------------------------------------------------------------------
// Function      : packedByteCount
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
template<>
int
Pack<Device::Param>::packedByteCount(
  const Device::Param &         param)
{
  // Util::Param info
  int byteCount = Xyce::packedByteCount(static_cast<const Util::Param &>(param));

  // given & default
  byteCount += sizeof(int);

  return byteCount;

}

//-----------------------------------------------------------------------------
// Function      : pack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
template<>
void
Pack<Device::Param>::pack(
  const Device::Param &    param,
  char *                   buf,
  const int                bsize,
  int &                    pos,
  Parallel::Communicator * comm )
{
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos + Xyce::packedByteCount(param);
#endif

  Xyce::pack(static_cast<const Util::Param &>(param), buf, bsize, pos, comm );

  //pack given_
  int dg = (param.isGiven_?1:0) + 2*(param.isDefault_ ? 1 : 0);
  comm->pack( &dg, 1, buf, bsize, pos );

#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    DevelFatal(*this, "Param::pack") << "Predicted pos does not match actual pos";
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : unpack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
template <>
void
Pack<Device::Param>::unpack(Device::Param &param, char * pB, int bsize, int & pos, Parallel::Communicator* comm )
{
  Xyce::unpack(static_cast<Util::Param &>(param), pB, bsize, pos, comm );

  //unpack given_
  int dg;
  comm->unpack( pB, bsize, pos, &dg, 1 );
  param.isGiven_ = ( dg%2 != 0 );
  param.isDefault_ = ( dg >= 2 );
}

} // namespace Xyce

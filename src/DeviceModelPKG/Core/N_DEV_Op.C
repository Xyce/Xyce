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
// Purpose        : Provide tools for accessing output data in parallel or
//                  serial
//
// Special Notes  :
//
// Creator        : David Baur
//
// Creation Date  : 11/15/2013
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_Op.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_DeviceMgr.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : DeviceMgrGlobalParameterOp::get
// Purpose       : get the current value of a global param from the device
//                 package
// Special Notes : It is inappropriate for the Device package to be in charge
//                 of global params, but that's where they are right now.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
DeviceMgrGlobalParameterOp::get(const DeviceMgrGlobalParameterOp &op, const Util::Op::OpData &op_data)
{
  return op.globalParameterValue_;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntityParameterOp::get
// Purpose       : get the current value of a device parameter from a device
//                 entity
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
DeviceEntityParameterOp::get(const DeviceEntityParameterOp &op, const Util::Op::OpData &op_data)
{
  double result;

  const_cast<DeviceEntity &>(op.deviceEntity_).getParam(op.deviceParameterName_, result);

  return result;
}

//-----------------------------------------------------------------------------
// Function      : ArtificialParameterOp::get
// Purpose       : get the current value of a device parameter from a device
//                 entity
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
ArtificialParameterOp::get(const ArtificialParameterOp &op, const Util::Op::OpData &op_data)
{
  return op.artificialParameter_.getValue(op.deviceManager_);
}

// //-----------------------------------------------------------------------------
// // Function      : DeviceMgrParameterOp::get
// // Purpose       : get the current value of a device parameter from the device
// //                 package
// // Special Notes :
// // Scope         : public
// // Creator       : David Baur, Raytheon
// // Creation Date : 11/15/2013
// //-----------------------------------------------------------------------------
// complex
// DeviceMgrParameterOp::get(const DeviceMgrParameterOp &op, const Util::Op::OpData &op_data)
// {
//   return op.deviceManager_.getParamNoReduce(op.deviceParameterName_);
// }

//-----------------------------------------------------------------------------
// Function      : DeviceMgrParameterOp::get
// Purpose       : get the current value of a device parameter from the device
//                 package
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
DeviceOptionsOp::get(const DeviceOptionsOp &op, const Util::Op::OpData &op_data)
{
  return op.deviceOptions_.gmin;
}

} // namespace Device
} // namespace Xyce

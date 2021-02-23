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
//
//
// Special Notes  :
//
//
// Creator        : David Baur
//
// Creation Date  :
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_ADC.h>
#include <N_DEV_Algorithm.h>
#include <N_DEV_Device.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternDevice.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : getName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Apr 22 12:22:31 2014
//-----------------------------------------------------------------------------
///
/// @param device_instance device instance to return name of
///
/// @return name of device instance
///
///
template<>
const std::string &getName(const DeviceInstance *device_instance) 
{
  return device_instance->getName().getEncodedName();
}


//-----------------------------------------------------------------------------
// Function      : getName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Apr 22 12:22:31 2014
//-----------------------------------------------------------------------------
///
/// @param device_model device model to return name of
///
/// @return name of device model
///
///
template<>
const std::string &getName(const DeviceModel *device_model) 
{
  return device_model->getName();
}


bool
devicesConverged(
  Parallel::Machine             comm,
  const InstanceVector &        extern_devices)
{
  int converged = true;
  for (InstanceVector::const_iterator it = extern_devices.begin(); it != extern_devices.end() && converged; ++it)
  {
    ExternDevice::Instance &extern_device = static_cast<ExternDevice::Instance &>(*(*it));

    converged = extern_device.isInnerSolveConverged();
  }

  Parallel::AllReduce(comm, MPI_LAND, &converged, 1);

  return converged;
}

} // namespace Device
} // namespace Xyce


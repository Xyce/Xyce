//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/15/02
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DeviceSensitivities_h
#define Xyce_N_DEV_DeviceSensitivities_h

#include <N_DEV_fwd.h>
#include <N_UTL_fwd.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : DeviceSensitivities
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/15/02
//-----------------------------------------------------------------------------
class DeviceSensitivities
{
public:
  DeviceSensitivities(
     DeviceMgr &               device_manager,
     const DeviceOptions &     device_options);
  ~DeviceSensitivities();

private:
  DeviceSensitivities(const DeviceSensitivities &);
  DeviceSensitivities &operator=(const DeviceSensitivities &);

public:
  bool registerSensParams(const Util::OptionBlock & OB);
  bool setSensitivityOptions (const Util::OptionBlock & OB);

  bool forceFD;
  bool forceFDgiven;

private:
  DeviceMgr &                 deviceManager_;
  const DeviceOptions &       deviceOptions_;
};

} // namespace Device
} // namespace Xyce

#endif


//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <algorithm>
#include <vector>
#include <list>

#include <N_DEV_DeviceSensitivities.h>

#include <N_DEV_DeviceMgr.h>
#include <N_DEV_Algorithm.h>

#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Algorithm.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_Param.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : DeviceSensitivities::DeviceSensitivities
// Purpose       : constructor
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/15/02
//-----------------------------------------------------------------------------
DeviceSensitivities::DeviceSensitivities(
  DeviceMgr &           device_manager,
  const DeviceOptions & device_options)
  : deviceManager_(device_manager),
    deviceOptions_(device_options)
{}

//-----------------------------------------------------------------------------
// Function      : DeviceSensitivities::~DeviceSensitivities
// Purpose       : destructor
// Special Notes : De-allocates all the devices pointed  to by deviceArray
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/15/02
//-----------------------------------------------------------------------------
DeviceSensitivities::~DeviceSensitivities()
{}

//-----------------------------------------------------------------------------
// Function      : DeviceSensitivites::registerSensParams
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/15/02
//-----------------------------------------------------------------------------
bool DeviceSensitivities::registerSensParams(const Util::OptionBlock & OB)
{
  bool bsuccess = true;

  if (DEBUG_DEVICE && isActive(Diag::SENS_PARAMETERS))
  {
    Xyce::dout() << "DeviceSensitivites::registerSensParams called!" <<std::endl;
  }
  int numSensParams = 0;

  for (Util::ParamList::const_iterator iter = OB.begin(); iter != OB.end(); ++iter)
  {
    if ( std::string(iter->uTag(), 0, 5) == "PARAM") // this is a vector
    {
      const std::string &tag = iter->stringValue();

      if (DEBUG_DEVICE)
      {
        Xyce::dout() << "name = " << iter->uTag() << "  tag = " << tag << std::endl;
      ++numSensParams;
      }
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::SENS_PARAMETERS))
  {
    Xyce::dout() << "number of sensitivity parameters = "<< numSensParams << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceSensitivites::setSensitivityOptions
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool DeviceSensitivities::setSensitivityOptions(const Util::OptionBlock & OB)
{
  if (DEBUG_DEVICE && isActive(Diag::SENS_PARAMETERS))
  {
    Xyce::dout() << "DeviceSensitivites::setSensitivityOptions called!" <<std::endl;
  }

  bool bsuccess = true;
  Util::ParamList::const_iterator it  = OB.begin();
  Util::ParamList::const_iterator end = OB.end();
  for ( ; it != end; ++ it)
  {
    if ((*it).uTag() == "FORCEFD")
    {
      forceFD = static_cast<bool>((*it).getImmutableValue<bool>());
      forceFDgiven = true;
    }
  }

  return bsuccess;
}

} // namespace Device
} // namespace Xyce

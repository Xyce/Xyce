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

#ifndef Xyce_N_UTL_Algorithm_h
#define Xyce_N_UTL_Algorithm_h

#include <string>

#include <N_DEV_fwd.h>
#include <N_DEV_InstanceName.h>
#include <N_UTL_HspiceBools.h>

namespace Xyce {
namespace Util {

inline
Device::InstanceName entityNameFromFullParamName(const std::string &full_param_name)
{
  std::string::size_type pos = full_param_name.find_last_of(Xyce::Util::separator);

  if (pos == std::string::npos)
    return Device::InstanceName(full_param_name);
  else
    return Device::InstanceName(std::string(full_param_name, 0, pos));
}

inline
Device::InstanceName entityNameFromDefaultParamName(const std::string &def_param_name)
{
  return Device::InstanceName(def_param_name);
}

inline
std::string paramNameFromFullParamName(const std::string &full_param_name)
{
  std::string::size_type pos = full_param_name.find_last_of(Xyce::Util::separator);

  if (pos == std::string::npos)
    return std::string();
  else
    return std::string(full_param_name, pos + 1);
}

} // namespace Util
} // namespace Xyce

#endif // Xyce_N_UTL_Algorithm_h

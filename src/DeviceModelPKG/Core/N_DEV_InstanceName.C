//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Creation Date  : 02/28/00
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_DeviceBlock.h>
#include <N_UTL_DeviceNameConverters.h>
#include <N_UTL_HspiceBools.h>

namespace Xyce {
namespace Device {

void
InstanceName::decode()
{
  deviceType_ = decodeDeviceType(*this);
  deviceName_ = decodeDeviceName(*this);
}

std::string decodeDeviceName(const InstanceName& instance_name) 
{
  // Skip subcircuit
  const std::string& name_ = instance_name.getEncodedName();
  std::string::size_type i = name_.find_last_of(Xyce::Util::separator);
  i = (i == std::string::npos ? 0 : i + 1);

  // For Y return the string following the the type terminating
  // exclamation point (!).  For U return the string between the type
  // terminating exclamation point (!) and the device name terminating
  // exclamation point (!).
  if (i < name_.size()) {
    if (name_[i] == 'Y')
      return name_.substr(name_.find("!") + 1);
    else if (name_[i] == 'U') {
      i = name_.find('!', i + 1);
      if (i < name_.size()) {
        std::string::size_type j = name_.find('!', i + 1);
        return name_.substr(i, j);
      }
    }
    else
      return name_.substr(i);
  }
  return "";
}

std::string decodeDeviceType(const InstanceName& instance_name)
{
  // Skip subcircuit
  const std::string& name_ = instance_name.getEncodedName();
  std::string::size_type i = name_.find_last_of(Xyce::Util::separator);
  i = (i == std::string::npos ? 0 : i + 1);

  // For multiletter (Y or U type), return the device type without the Y
  // or U.  For single letter, just return it.
  if (i < name_.size()) {
    if (name_[i] == 'Y' || name_[i] == 'U')
      return name_.substr(i + 1, name_.find("!", i) - i - 1);
    else
      return name_.substr(i, 1);
  }
  else
    return "";
}

std::string getSubcircuitName(const InstanceName& instance_name)
{
  const std::string& name_ = instance_name.getEncodedName();
  std::string::size_type i = name_.find_last_of(Xyce::Util::separator);

  if (i != std::string::npos)
    return name_.substr(0, i);
  else
    return "";
}

/// Return the first letter of the device specification after the subcircuit.
char getDeviceLetter(const InstanceName& instance_name)
{
  const std::string& name_ = instance_name.getEncodedName();

  // Skip subcircuit
  std::string::size_type i = name_.find_last_of(Xyce::Util::separator);
  i = (i == std::string::npos ? 0 : i + 1);

  // Return first letter after subcircuit.
  return i < name_.size() ? name_[i] : '\0';
}

/// Return the number of inputs which have been encoded into the device name.
int getNumInputs(const InstanceName& instance_name)
{
  const std::string& name_ = instance_name.getEncodedName();

  if (getDeviceLetter(instance_name) == 'U') {
    std::string::size_type i = name_.find_last_of("!");
    if (i == std::string::npos)
      return 0;

    std::istringstream iss(name_.substr(i + 1));
    int num_inputs = 0;

    iss >> num_inputs;

    return num_inputs;
  }

  return 0;
}

// ----------------------------------------------------------------------------
// Function      : setupOutputName
//
// Purpose       : This function takes the device instance name and creates
//                 an appropriate "outputName" to be used for file outputs.
//
//                 At this point PDE devices are all specified as "Y" devices,
//                 meaning that the device instance name will almost always
//                 start with "YPDE!".  Left unchanged, this has been
//                 resulting in (for example) tecplot files named thing like,
//                 "YPDEDIODE1.dat".  Obviously, the YPDE prefix is
//                 not needed, so this function removes it, if it exists.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/03/05
// ----------------------------------------------------------------------------
std::string setupOutputName(const InstanceName &name)
{
  std::string outputName;

  std::string s = name.getEncodedName();

  std::string pdeString("YPDE!");
  std::string neutString("YNEUTRON!");
  std::string::size_type pos1 = s.find(pdeString);
  std::string::size_type pos2 = s.find(neutString);

  if (pos1 != std::string::npos)
  {
    std::string tmp1 = "";
    if (pos1 > 0)
      tmp1 = s.substr(0,pos1);
    std::string tmp2 = s.substr(pos1+5, s.length()-1);
    outputName = tmp1 + tmp2;
  }
  else if (pos2 != std::string::npos)
  {
    std::string tmp1 = "";
    if (pos2 > 0) tmp1 = s.substr(0,pos2);
    std::string tmp2 = s.substr(pos2+9, s.length()-1);
    outputName = tmp1 + tmp2;
  }
  else
  {
    outputName = s;
  }

  // Tecplot doesn't like file names with the character, ":", so change all
  // colons to underscores.  I personally don't like "%" characters in
  // filenames, so remove those as well.
  for (int i=0;i<outputName.size();++i)
  {
    if (outputName[i]==':') outputName[i] = '_';
    if (outputName[i]=='%') outputName[i] = '_';
  }

  return outputName;
}

//-----------------------------------------------------------------------------
// Function      : spiceInternalName
// Purpose       : Converts a conventional "Xyce style" internal variable
//                 name to a spice style name.
// Special Notes :
//                 Example:
//                 chilespice:    d:subckt1:test1_internal
//                 xyce:          xsubckt1:dtest1_internal
//
//                 This function will not remove subcircuit "x" characters,
//                 but does move the first character of the local name
//                 to the beginning of the full std::string.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/12/01
//-----------------------------------------------------------------------------
std::string spiceInternalName(const InstanceName &entity_name, const std::string &lead)
{
  std::string s = entity_name.getEncodedName();

  s = s + "_" + lead;

  s = Util::xyceDeviceNameToSpiceName(s);
  
  return s;
}

std::string spiceStoreName(const InstanceName &entity_name, const std::string &lead)
{
  std::string s = entity_name.getEncodedName();

  s = Util::xyceDeviceNameToSpiceName(s);
  return s + Xyce::Util::separator + lead;
}

std::ostream &operator<<(std::ostream &os, const InstanceName &entity_name) {
  return os << entity_name.getEncodedName();
}

} // namespace Device
} // namespace Xyce

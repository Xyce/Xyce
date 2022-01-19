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

//-------------------------------------------------------------------------
//
// Purpose        : Forward declarations
//
// Special Notes  : Forward declaring everything as a class breaks if the implementation of the type changes (like during
//                  templatization)
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355
//
// Creation Date  : 2013/04/18 18:01:27
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_InstanceName_h
#define Xyce_N_DEV_InstanceName_h

#include <string>

#include <N_DEV_fwd.h>
#include <N_UTL_NoCase.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : InstanceName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon May 12 15:53:34 2014
//-----------------------------------------------------------------------------
///
/// Devices and models are each named.  Models are not encoded, so
/// simple string representation is sufficient.  However Devices are a
/// different lot.  They are encoded as
///
/// [s:]*xname
/// [s:]*Ytype!name
/// [s:]*Utype!name!count
///
/// where s is a subcircuit name, x is a single letter device type, type
/// is a multiletter device type (no Y or U prefix) and count is a
/// special input count for the U device.
///
/// Currently encoded names are accepted and then decoded using the
/// getter's.  In the future, these will be stored in component form and
/// then encoded onyl as needed.
///
class InstanceName
{
public:
  //-----------------------------------------------------------------------------
  // Function      : InstanceName
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon May 12 15:57:34 2014
  //-----------------------------------------------------------------------------
  ///
  /// Creates an empty entity name.
  ///
  ///
  InstanceName()
    : deviceType_(),
      deviceName_(),
      name_()
  {}

  //-----------------------------------------------------------------------------
  // Function      : InstanceName
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon May 12 15:58:19 2014
  //-----------------------------------------------------------------------------
  ///
  /// Creates an entity
  ///
  /// @invariant
  ///
  /// @param name model or encoded device entity name
  ///
  /// @return 
  ///
  ///
  explicit InstanceName(const std::string &name)
    : deviceType_(),
      deviceName_(),
      name_(name)
  {
    decode();
  }

  //-----------------------------------------------------------------------------
  // Function      : InstanceName
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon May 12 15:58:19 2014
  //-----------------------------------------------------------------------------
  ///
  /// Copies an instance name
  ///
  /// @invariant
  ///
  /// @param name model or encoded device entity name
  ///
  /// @return 
  ///
  ///
  InstanceName(const InstanceName &entity_name)
    : deviceType_(entity_name.deviceType_),
      deviceName_(entity_name.deviceName_),
      name_(entity_name.name_)
  {}

  //-----------------------------------------------------------------------------
  // Function      : InstanceName::operator=
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon May 12 15:58:19 2014
  //-----------------------------------------------------------------------------
  ///
  /// Assigns an instance name
  ///
  /// @invariant
  ///
  /// @param name model or encoded device entity name
  ///
  /// @return 
  ///
  ///
  InstanceName &operator=(const InstanceName &instance_name)
  {
    deviceType_ = instance_name.deviceType_;
    deviceName_ = instance_name.deviceName_;
    name_ = instance_name.name_;

    return *this;
  }

  //-----------------------------------------------------------------------------
  // Function      : getDeviceName
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon May 12 16:04:37 2014
  //-----------------------------------------------------------------------------
  ///
  /// Decodes the device name. 
  ///
  /// The device name is a string containing the first letter for a
  /// single letter device name, or the string starting after the Y or
  /// U up to but noe including the exclamation point (!).
  ///
  /// Subcircuit prefixes are not included.
  ///
  /// @return string representing the device name
  ///
  ///
  const std::string &getDeviceName() const
  {
    return deviceName_;
  }

  //-----------------------------------------------------------------------------
  // Function      : getDeviceType
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon May 12 16:02:24 2014
  //-----------------------------------------------------------------------------
  ///
  /// Decodes the device type.
  ///
  /// The device type is a string containing the first letter for a
  /// single letter device name, or the string starting after the Y or
  /// U up to but noe including the exclamation point (!).
  ///
  /// Subcircuit prefixes are not included.
  ///
  /// @return string representing the device type
  ///
  ///
  const std::string &getDeviceType() const
  {
    return deviceType_;
  }

  //-----------------------------------------------------------------------------
  // Function      : getEncodedName
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon May 19 07:40:55 2014
  //-----------------------------------------------------------------------------
  ///
  /// Return the instance name encoded as:
  ///   [s:]*xname
  ///   [s:]*Ytype!name
  ///   [s:]*Utype!name!count
  ///
  /// @return encoded instance name
  ///
  ///
  const std::string &getEncodedName() const
  {
    return name_;
  }

private:
  void decode();

private:
  std::string         deviceType_;            ///< Device type from netlist (does NOT include Y or U prefix)
  std::string         deviceName_;            ///< Device name from netlist (subcircuit NOT included)
  std::string         name_;                  ///< Complete encoded name
};

inline bool operator==(const InstanceName &instance_name, const std::string &name)
{
  return equal_nocase(instance_name.getEncodedName(), name);
}

inline bool operator<(const InstanceName &instance_name, const std::string &name)
{
  return less_nocase(instance_name.getEncodedName(), name);
}

inline bool operator!=(const InstanceName &instance_name, const std::string &name)
{
  return !equal_nocase(instance_name.getEncodedName(),  name);
}

std::string setupOutputName(const InstanceName &name);

std::ostream &operator<<(std::ostream &os, const InstanceName &instance_name);

/// Return the first letter of the device specification after the subcircuit.
char getDeviceLetter(const InstanceName& instance_name);

/// Decodes the subcircuit prefix to the device.
std::string getSubcircuitName(const InstanceName& instance_name);

/// Decodes the device type, returns character
std::string decodeDeviceType(const InstanceName& instance_name);

/// Decodes the device name, returns name string
std::string decodeDeviceName(const InstanceName& instance_name);

/// Return the number of inputs which have been encoded into the device name.
int getNumInputs(const InstanceName& instance_name);

std::string spiceInternalName(const InstanceName &instance_name, const std::string &lead);
std::string spiceStoreName(const InstanceName &instance_name, const std::string &lead);

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_InstanceName_h


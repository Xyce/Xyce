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
// Purpose        : Provide utility functions to convert device names
//                  with circuit context from Xyce to SPICE style, and back
//
// Special Notes  : A version of the Xyce-to-SPICE converter here used to
//                  live in the Device package, but was removed and replaced
//                  by something very device-instance specific.
//                  There should not be any reason for a SPICE-to-Xyce
//                  converter, except for a bad design choice that
//                  saves internal variable names (branch currents, etc.) using
//                  the SPICE name, and one specific spot where this breaks
//                  everything (c.f. bugs 715 and 982 on the SON bugzilla).
//                  When bug 982 is addressed, it is likely that the
//                  SPICE-to-Xyce function can be discarded
//
// Creator        : Tom Russo
//
// Creation Date  : 14 March 2018
//
//-------------------------------------------------------------------------

#include <N_UTL_DeviceNameConverters.h>
#include <N_UTL_HspiceBools.h>

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Function      : xyceDeviceNameToSpiceName
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
// Scope         : free function
// Creator       : Tom Russo
// Creation Date : 14 March 2018
//-----------------------------------------------------------------------------
/// given a Xyce-style device internal variable, create the SPICE-style equivalent
///
/// SPICE makes sure to keep the device type character the first character of
/// a device's name, even if it's inside a subcircuit.  Xyce uses a different
/// style of name, where the circuit context (path through subcircuits) comes
/// first, and the full device name (including device type character) is last.
///
/// For reasons related primarily to output of full-solution rawfiles and
/// direct comparison with SPICE, it is convenient to be able to get at the
/// SPICE-style name given a Xyce-style name.  Ideally, this SPICE name
/// should only be used for output purposes, but there are some cases where
/// the SPICE name is stored in lookup maps.  This function is intended to be
/// used everywhere that the Xyce name is given but the SPICE name is needed.
///
/// Example Xyce name:     xsub1:xsub2:dtest1_internal
/// SPICE name:            d:sub1:sub2:test1_internal
/// This function doesn't remove the "x" from the subcircuit context, but
/// does move the device type character.
///
/// @param[in] xdName    Xyce device variable name
/// @returns SPICE-style name string.
std::string xyceDeviceNameToSpiceName(std::string & xdName)
{
  std::string modifiedName="";
  std::string::size_type lastColonInName = xdName.find_last_of(Xyce::Util::separator);
  if ((lastColonInName != std::string::npos) && (lastColonInName + 1 < xdName.length()))
  {
    std::string::iterator deviceVarName = xdName.begin() + lastColonInName+1;
    std::string::iterator namePrefixEnd = xdName.begin() + lastColonInName;
    modifiedName.append(deviceVarName, deviceVarName + 1);
    modifiedName.append(":");
    modifiedName.append(xdName.begin(), namePrefixEnd + 1);
    modifiedName.append(deviceVarName + 1, xdName.end());
  }
  else
  {
    modifiedName = xdName;
  }

  if (Xyce::Util::useHspiceSeparator)
  {
    std::replace(modifiedName.begin(),modifiedName.end(),Xyce::Util::separator,':');
  }

  return modifiedName;
}

//-----------------------------------------------------------------------------
// Function      : spiceDeviceNameToXyceName
// Purpose       : Converts a "SPICE style" internal variable to a Xyce-style
//                 name.
// Special Notes :
//                 Example:
//                 SPICE:         d:xsubckt1:test1_internal
//                 xyce:          xsubckt1:dtest1_internal
//                 Note that real SPICE names don't have the "x", but
//                 those we actually need to work with here WILL.
//
// Scope         : free function
// Creator       : Tom Russo
// Creation Date : 14 March 2018
//-----------------------------------------------------------------------------
/// given a SPICE-style device internal variable, create the Xyce-style equivalent
///
/// There is no good reason for this function to exist, as it is here
/// only to work around a bad design decision documented in bug 982 on
/// the SON.  A decision was made to store the SPICE-style name of
/// internal device variables in SPICE form in a lookup map.  Later,
/// when print lines are handled, it is therefore necessary that
/// Xyce-style device names given on the print line have to be
/// converted into SPICE style in order to look them up to obtain data to print.
///
/// An unfortunate effect of  this is that in exactly one place in the code,
/// where print lines that contain "I(*)" are processed, the I(*) is replaced
/// with an explicit list of all branch currents, taken from the very map that
/// contains SPICE names.  Unfortunately, that breaks the bit of output code
/// that tries to rewrite Xyce names into SPICE names (because they're
/// already in SPICE style).  So to fake all of this out, we have to use
/// this function in that one place of the output manager, to recover
/// the Xyce name from the stored SPICE name.
///  
/// Once bug 982 is addressed, this function can go away.
///
/// @param[in] sdName  device variable name in SPICE style.
/// @returns device variable name in Xyce style.
///
std::string spiceDeviceNameToXyceName(std::string & sdName)
{
  std::string modifiedName="";
  std::string::size_type firstColonInName = sdName.find_first_of(":");
  std::string::size_type lastColonInName = sdName.find_last_of(":");

  if (sdName[0] != 'X' && (firstColonInName != std::string::npos) && (lastColonInName != std::string::npos) && (firstColonInName != lastColonInName) )
  {
    std::string::iterator deviceType = sdName.begin();
    std::string::iterator circuitContextBegin = sdName.begin()+firstColonInName+1;
    std::string::iterator circuitContextEnd = sdName.begin()+lastColonInName;
    std::string::iterator deviceName = sdName.begin()+lastColonInName+1;

    modifiedName.append(circuitContextBegin,circuitContextEnd);
    modifiedName.append(":");
    modifiedName.append(deviceType,deviceType+1);
    modifiedName.append(deviceName,sdName.end());
  }
  else
  {
    modifiedName=sdName;
  }

  if (Xyce::Util::useHspiceSeparator)
  {
    std::replace(modifiedName.begin(),modifiedName.end(), ':', Xyce::Util::separator);
  }

  return modifiedName;
}
}
}

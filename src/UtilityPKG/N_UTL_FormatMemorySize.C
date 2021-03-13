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

#include <sstream>
#include <iomanip>
#include <N_UTL_Math.h>

#include <N_UTL_FormatMemorySize.h>

#define KILO_ONLY

namespace Xyce {

//-----------------------------------------------------------------------------
// Function      : formatMemorySize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Sep  2 17:06:44 2014
//-----------------------------------------------------------------------------
std::string
formatMemorySize(
  double                size)
{
  std::ostringstream os;

  static const double kb = 1024.0;
  // static const double mb = kb * kb;
  // static const double gb = kb * kb * kb;

  if (size < 0.0) {
    os << "-";
    size = -size;
  }

#ifdef KILO_ONLY
  // output size in kilo bytes
  os << static_cast<unsigned long>(size/kb) << " KB";

#else
  if (size < kb) {
    // output size in bytes
    result += boost::lexical_cast<std::string>(static_cast<unsigned long>(size));
    result += " B";
  }
  else if (size < mb) {
    // output size in kilo bytes
    result += boost::lexical_cast<std::string>(static_cast<unsigned long>(size / kb));
    result += " KB";
  }
  else if (size < gb) {
    // output size in mega bytes
    result += boost::lexical_cast<std::string>(static_cast<unsigned long>(size / mb));
    result += " MB";
  }
  else {
    // everything else output in giga bytes
    result += boost::lexical_cast<std::string>(static_cast<unsigned long>(size / gb));
    result += " GB";
  }
#endif

  return os.str();
}


std::string
formatMemorySize(
  MemorySize            size)
{
  std::ostringstream os;

  static const MemorySize kb = 1024;
  // static const MemorySize mb = kb * kb;
  // static const MemorySize gb = kb * kb * kb;

#ifdef KILO_ONLY
  // output size in kilo bytes
  os << size/kb << " KB";

#else
  if (size < kb) {
    // output size in bytes
    result = boost::lexical_cast<std::string>(size);
    result += " B";
  }
  else if (size < mb) {
    // output size in kilo bytes
    result = boost::lexical_cast<std::string>(size / kb);
    result += " KB";
  }
  else if (size < gb) {
    // output size in mega bytes
    result = boost::lexical_cast<std::string>(size / mb);
    result += " MB";
  }
  else {
    // everything else output in giga bytes
    result = boost::lexical_cast<std::string>(size / gb);
    result += " GB";
  }
#endif

  return os.str();
}

} // namespace Xyce

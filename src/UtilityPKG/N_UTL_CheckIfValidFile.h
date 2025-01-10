//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : Verify that a user-specified file (e.g, netlist file on
//                  the Xyce command line) is valid and can be opened.  This
//                  function will return false if the file does not exist, cannot 
//                  be opened, or if the user accidentally specified a directory 
//                  name rather than a file name.  See SON Bugs 730 and 785 
//                  for more details.
// Special Notes  : 
// Creator        : Pete Sholander Sandia National Laboratories, 1355 
// Creation Date  : 2017/07/7 18:01:27
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_CheckIfValidFile_h
#define Xyce_N_UTL_CheckIfValidFile_h

#include <string>

namespace Xyce {
namespace Util {

  bool checkIfValidFile(std::string netlist_filename);
  bool checkIfValidDashoFileName(std::string dashoFilename);

} // namespace Util
} // namespace Xyce

#endif // Xyce_N_UTL_CheckIfValidFile_h

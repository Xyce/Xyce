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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 5/22/00
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_NetlistLocation_h
#define Xyce_N_UTL_NetlistLocation_h

#include <iosfwd>
#include <string>

#include <N_UTL_Marshal.h>

namespace Xyce {

struct NetlistLocation
{
  NetlistLocation();

  NetlistLocation(const std::string &filename, const int line_number);

  NetlistLocation(const NetlistLocation &netlist_location)
    : fileNumber_(netlist_location.fileNumber_),
      lineNumber_(netlist_location.lineNumber_)
  {}

  NetlistLocation &operator=(const NetlistLocation &netlist_location) 
  {
    fileNumber_ = netlist_location.fileNumber_;
    lineNumber_ = netlist_location.lineNumber_;

    return *this;
  }

  NetlistLocation &setFilename(const std::string &filename); 

  const std::string &getFilename() const;

  NetlistLocation &setFileNumber(int file_number)
  {
    fileNumber_ = file_number;
    return *this;
  }

  int getFileNumber() const
  {
    return fileNumber_;
  }

  NetlistLocation &setLineNumber(int line_number) 
  {
    lineNumber_ = line_number;
    return *this;
  }

  int getLineNumber() const 
  {
    return lineNumber_;
  }

private:
  int           fileNumber_;
  int           lineNumber_;
};

bool operator<(const NetlistLocation &s0, const NetlistLocation &s1);

std::ostream &operator<<(std::ostream &os, const NetlistLocation &netlist_locationc);

Util::Marshal &operator<<(Util::Marshal &mout, const NetlistLocation &netlist_location);

Util::Marshal &operator>>(Util::Marshal &min, NetlistLocation &netlist_location);

class Filename
{
public:
  typedef std::vector<std::string> FilenameVector;

  /// Returns the filename vector
  ///
  /// The filename vector contains all the netlist filenames.
  ///
  /// @return const reference to the filename vector.
  static const FilenameVector &getFilenameVector();

  /// Returns the filename associated with the netlist location
  ///
  /// @param loc               netlist location
  ///
  /// @return const reference to the filename
  ///
  static const std::string& getFilename(const NetlistLocation& loc);

  /// Returns the file number associated with the filename
  ///
  /// @param file_name         file name
  ///
  /// @return integer of the line number in the filename
  static int getFileNumber(const std::string& filename);
};

} // namespace Xyce

#endif // Xyce_N_UTL_NetlistLocation_h


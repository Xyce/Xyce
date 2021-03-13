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
// Purpose        : This file contains specializations for the APFT
//                  for various vector types.
//
// Special Notes  :
//
// Creator        : David G Baur
//
// Creation Date  : 4/9/14
//
//-------------------------------------------------------------------------

#ifndef N_UTL_Dump_h
#define N_UTL_Dump_h

namespace Xyce {

template<class T>
class Dump
{
public:
  Dump(const T &t)
    : t_(t)
  {}

  std::ostream &dump(std::ostream &os) const;

  const T &     t_;
};

template<class T>
std::ostream &operator<<(std::ostream &os, const Dump<T> &dumper) {
  return dumper.dump(os);
}

} // namespace Xyce

#endif // N_UTL_Dump_h

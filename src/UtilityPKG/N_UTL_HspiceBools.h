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


#ifndef N_UTL_HspiceBools_H
#define N_UTL_HspiceBools_H

namespace Xyce {
namespace Util {

// This value is derived from the -hspice-ext command line option.  It is
// set, based on that command line option, in the constructor for the
// IO::ParsingMgr class.  If set to true then logical AND is &&, logical
// OR is || and ^ is a synonym for exponentiation.
extern bool useHspiceMath;
extern bool useHspiceSeparator;
extern char separator;
extern std::vector<bool> preprocessFilter;
}
}

#endif


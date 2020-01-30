//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
// Purpose        : Demangle C++ type_info names
//
// Special Notes  : 
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_Demangle_hpp
#define Xyce_N_UTL_Demangle_hpp

#include <string>

#if __GNUC__ == 3 || __GNUC__ == 4 || defined __xlC__
#define Xyce__USE_PLATFORM_DEMANGLER
#endif

namespace Xyce {

/**
 * @brief Function <b>demangle</b> returns the demangled C++ symbol from the mangled
 * C++ symbol.  The mangled named is obtained from the <b>type_info</b>
 * <b>name()</b> function.  From some compilers, the name is already demangled.
 *
 * @param symbol	a <b>char</b> const pointer to the symbol.
 *
 * @return		a <b>std::string</b> value of the demangled name.
 */
std::string demangle(const char *symbol);

} // namespace Xyce

#endif // Xyce_N_UTL_Demangle_hpp

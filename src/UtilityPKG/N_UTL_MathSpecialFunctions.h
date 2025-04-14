//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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

#ifndef Xyce_N_UTL_MathSpecialFunctions_h
#define Xyce_N_UTL_MathSpecialFunctions_h

namespace Xyce {
namespace Util {

  // Modified bessel functions of the first kind
  double besselI0(double x);
  double besselI1(double x);
  double besselI1xOverX(double x);
  
  // compute erfcx(z) = exp(z^2) erfc(z)
  double erfcx(double x); // special case for real x

  // compute erf(z), the error function of complex arguments
  double erf(double x); // special case for real x

  // compute erfc(z) = 1 - erf(z), the complementary error function
  double erfc(double x); // special case for real x


} // Util namespace
} // Xyce namespace



#endif // Xyce_N_UTL_Math_h

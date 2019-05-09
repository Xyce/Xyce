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

//-----------------------------------------------------------------------------
//
// Purpose       : This file defines the machine dependent parameters in terms
//                 of those returned by "numeric_limits"
//
// Special Notes : This class tacitly assumes that "numeric_limits" is 
//                 available, which is the case in all current C++ compilers
//                 in 2012.  It replaces the old class 
//                 "N_TIA_MachineDependentParams" that was inappropriately
//                 under the TimeIntegrationPKG, and which was so old that 
//                 it labored under the burden of having to support compilers
//                 that could not be counted on to provide numeric_limits *OR*
//                 reliable template support.  In 2012, Trilinos is heavily
//                 dependent on the presence of numeric_limits, and doesn't 
//                 even bother checking if it is present during its 
//                 configuration by cmake.  We are highly dependent on 
//                 Trilinos.  Thus, it is OK to drop all pretense that Xyce 
//                 can be built with a compiler that does not provide 
//                 numeric_limits.
//
// Creator       : Tom Russo
//
// Creation Date : 8/1/2012
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_UTL_MACH_DEP_PARAMS_H_
#define Xyce_N_UTL_MACH_DEP_PARAMS_H_

#include <N_UTL_Math.h>
#include <limits>

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Class         : N_UTL_MachineDependentParams
// Purpose       :
// Special Notes :
// Creator       : Tom Russo
// Creation Date : 8/01/2012
//-----------------------------------------------------------------------------
class MachineDependentParams
{
public:

  // 10 * Minimum Floating Point Value
  inline static double MachineZero()
  { return 10.0 * std::numeric_limits<double>::min(); }

  // SquareRoot(Maximum Floating Point Value)
  inline static double MachineBig()
  { return sqrt(std::numeric_limits<double>::max()); }

  // Machine precision
  // This is NOT 4*epsilon as N_TIA_MachineDependentParams used to call it,
  // because the N_TIA_NumericalLimits class actually defined epsilon 
  // differently than the C++ std::numeric_limits does.  We do here precisely what 
  // the older class used to, because it turns out to matter in a few
  // cases where this MachinePrecision is used.  
  // 4*std::numeric_limits<double>::epsilon turns out to be even smaller than this.
  inline static double MachinePrecision()
  { return 2 * pow(10.0, -(std::numeric_limits<double>::digits10)); }

  // Machine Unit Roundoff
  inline static double MachineEpsilon()
  { return std::numeric_limits<double>::epsilon(); }

  inline static double DoubleMin()
  { return std::numeric_limits<double>::min(); }
  inline static double DoubleMax()
  { return std::numeric_limits<double>::max(); }
  inline static int IntMin()
  { return std::numeric_limits<int>::min(); }
  inline static int IntMax()
  { return std::numeric_limits<int>::max(); }
};

} // namespace Util
} // namespace Xyce

#endif     // Xyce_N_UTL_MACH_DEP_PARAMS_H

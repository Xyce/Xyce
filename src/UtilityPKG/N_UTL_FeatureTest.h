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
// Purpose        : This file is similar to Petra_Petra.h.  It  contains
//                  a number of Xyce-wide definitions and includes.
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

///
/// @file   N_UTL_FeatureTest.h
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
/// @date   Wed Jan 28 11:18:31 2015
///
/// @brief  Feature test translation from c pre-processor directive to Xyce namespace
/// static constants.
///
/// @see N_UTL_Diagnostic.h
///
/// Since using #if/#ifdef statements does not run the code through the compiler, the
/// maintenance of that code becomes difficult.  To have the code compile and not be
/// included in optimized, utilizes these static const int definition in a simple if
/// statement where shortcutting results in the optimizer trimming the code from the
/// object file.
///
/// Recommended use is:
///
///   if (DEBUG_DEVICE && Diag::isActive(Diag::DEVICE_PARAMETER))
///   {
///     .
///     .
///     .
///   }
///
/// Note that for a particular module, you can simply define the same state const int in a
/// nearer scope than namespace Xyce.  Please be sure to remove such uses prior to pushing
/// back to the repository.
///

#ifndef Xyce_N_UTL_FeatureTest_h
#define Xyce_N_UTL_FeatureTest_h

namespace Xyce {

#ifdef Xyce_DEBUG_DEVICE
static const int DEBUG_DEVICE = 1;
#else
static const int DEBUG_DEVICE = 0;
#endif

#ifdef Xyce_DEBUG_ANALYSIS
static const int DEBUG_ANALYSIS = 1;
#else
static const int DEBUG_ANALYSIS = 0;
#endif

#ifdef Xyce_DEBUG_ES
static const int DEBUG_ES = 1;
#else
static const int DEBUG_ES = 0;
#endif

#ifdef Xyce_DEBUG_PCE
static const int DEBUG_PCE = 1;
#else
static const int DEBUG_PCE = 0;
#endif

#ifdef Xyce_DEBUG_SAMPLING
static const int DEBUG_SAMPLING = 1;
#else
static const int DEBUG_SAMPLING = 0;
#endif

#ifdef Xyce_DEBUG_HB
static const int DEBUG_HB = 1;
#else
static const int DEBUG_HB = 0;
#endif

#ifdef Xyce_DEBUG_MPDE
static const int DEBUG_MPDE = 1;
#else
static const int DEBUG_MPDE = 0;
#endif

#ifdef Xyce_DEBUG_MOR
static const int DEBUG_MOR = 1;
#else
static const int DEBUG_MOR = 0;
#endif

#ifdef Xyce_DEBUG_IO
static const int DEBUG_IO = 1;
#else
static const int DEBUG_IO = 0;
#endif

#ifdef Xyce_DEBUG_EXPRESSION
static const int DEBUG_EXPRESSION = 1;
#else
static const int DEBUG_EXPRESSION = 0;
#endif

#ifdef Xyce_DEBUG_CONDUCTANCE
static const int DEBUG_CONDUCTANCE = 1;
#else
static const int DEBUG_CONDUCTANCE = 0;
#endif

#ifdef Xyce_VERBOSE_CONDUCTANCE
static const int VERBOSE_CONDUCTANCE = 1;
#else
static const int VERBOSE_CONDUCTANCE = 0;
#endif

#ifdef Xyce_DEBUG_RESTART
static const int DEBUG_RESTART = 1;
#else
static const int DEBUG_RESTART = 0;
#endif

#ifdef Xyce_DEBUG_TIME
static const int DEBUG_TIME = 1;
#else
static const int DEBUG_TIME = 0;
#endif

#ifdef Xyce_VERBOSE_TIME
static const int VERBOSE_TIME = 1;
#else
static const int VERBOSE_TIME = 0;
#endif

#ifdef Xyce_DEBUG_CIRCUIT
static const int DEBUG_CIRCUIT = 1;
#else
static const int DEBUG_CIRCUIT = 0;
#endif

#ifdef Xyce_DEBUG_NONLINEAR
static const int DEBUG_NONLINEAR = 1;
#else
static const int DEBUG_NONLINEAR = 0;
#endif

#ifdef Xyce_VERBOSE_NONLINEAR
static const int VERBOSE_NONLINEAR = 1;
#else
static const int VERBOSE_NONLINEAR = 0;
#endif

#ifdef Xyce_DEBUG_VOLTLIM
static const int DEBUG_VOLTLIM = 1;
#else
static const int DEBUG_VOLTLIM = 0;
#endif

#ifdef Xyce_DEBUG_LINEAR
static const int DEBUG_LINEAR = 1;
#else
static const int DEBUG_LINEAR = 0;
#endif

#ifdef Xyce_VERBOSE_LINEAR
static const int VERBOSE_LINEAR = 1;
#else
static const int VERBOSE_LINEAR = 0;
#endif

#ifdef Xyce_DEBUG_TOPOLOGY
static const int DEBUG_TOPOLOGY = 1;
#else
static const int DEBUG_TOPOLOGY = 0;
#endif

#ifdef Xyce_DEBUG_PARALLEL
static const int DEBUG_PARALLEL = 1;
#else
static const int DEBUG_PARALLEL = 0;
#endif

#ifdef Xyce_DEBUG_DISTRIBUTION
static const int DEBUG_DISTRIBUTION = 1;
#else
static const int DEBUG_DISTRIBUTION = 0;
#endif

#ifdef Xyce_Dakota_Debug
static const int DEBUG_DAKOTA = 1;
#else
static const int DEBUG_DAKOTA = 0;
#endif

#ifdef Xyce_Dakota
static const int DAKOTA = 1;
#else
static const int DAKOTA = 0;
#endif

#ifdef Xyce_DEBUG_ALL_PROCS_SAME_WD
static const int DEBUG_ALL_PROCS_SAME_WD = 1;
#else
static const int DEBUG_ALL_PROCS_SAME_WD = 0;
#endif

} // namespace Xyce

#endif // Xyce_N_UTL_FeatureTest_h

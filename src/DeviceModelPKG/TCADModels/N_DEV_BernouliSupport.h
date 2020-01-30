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

//-----------------------------------------------------------------------------
//
// Purpose        :  This class contains the range limit constants and
//                   supporting functions for the
//                   Bernoulli, Aux1 and Aux2 functions as well as their
//                   derivatives.
//
//                   This is derived (somewhat) from the SG Framework
//                   program, brkpnts.c, by Kevin M. Kramer.
//
//                   For more information consult "Analysis
//                   and Simulation of Electronic Devices" by Siegfried
//                   Selberherr pages 158 and 159.
//
// Special Notes  :  From the brkpnts.c program:
//
// This program determines how to evaluate the following functions to machine
// precision.  The breakpoints of these functions are dependent upon the host
// computer's floating point architecture and floating point library.
//
//
//                 x
//       B(x) = -------
//              e^x - 1
//
//
//       d      (1-x)*e^x - 1
//       --B(x) = -------------
//       dx       (e^x - 1)^2
//
//
//                     x
//       Aux1(x) =  -------
//                  sinh(x)
//
//       d           sinh(x) - x*cosh(x)
//       --Aux1(x) = -------------------
//       dx             (sinh(x))^2
//
//                    1
//       Aux2(x) = -------
//                 1 + e^x
//
//       d             - e^x
//       --Aux2(x) = -----------
//       dx          (1 + e^x)^2
//
//
// To achieve machine precision, the above functions should  be evaluated as
// follows:
//
//
//         /
//        | -x                                                x <=  +0.0e+00
//        |  x / (e^x - 1)                        +0.0e+00 <  x <   +0.0e+00
// B(x) =<   1 - x/2*(1 - x/6*(1 - x*x/60))       +0.0e+00 <= x <= +2.7e-314
//        |  x*e^-x / (1 - e^-x)                 +2.7e-314 <  x <   +0.0e+00
//        |  x*e^-x                               +0.0e+00 <= x <   +0.0e+00
//        |  0                                    +0.0e+00 <= x
//         \
//
//          /
//         | -1                                               x <= +3.7e-314
//         |  {(1-x)*e^x - 1}                    +3.7e-314 <  x <=  -2.0e+00
// d       |  {(1-x)*e^x - 1} / (e^x - 1)^2       -2.0e+00 <  x <  -3.1e-231
// --B(x)=<  -1/2 + x/6*(1 - x*x/30)             -3.1e-231 <= x <= -4.3e-232
// dx      |  {(1-x)*e^-x - e^-2x}/(1 - e^-x)^2  -4.3e-232 <  x <   +0.0e+00
//         |  {(1-x)*e^-x - e^-2x}                +0.0e+00 <= x <   -2.0e+00
//         |  0                                   -2.0e+00 <= x
//          \
//
//          /
//         |  x / sinh(x)                                     x <=  -8.0e-03
// Aux1(x)=<  1 - x*x/6*(1 - 7*x*x/60)            -8.0e-03 <  x <   +8.0e-03
//         |  x / sinh(x)                         +8.0e-03 <= x
//          \
//
//            /
// d          |  {sinh(x) - x*cosh(x)}/{sinh(x)}^2              x <=  -4.8e-03
// --Aux1(x)=<  -x/3*(1 - 7*x*x/30)                 -4.8e-03 <  x <   +4.8e-03
// dx         |  {sinh(x) - x*cosh(x)}/{sinh(x)}^2  +4.8e-03 <= x
//            \
//
//            /
//            |  1                                              x <=  -3.7e+01
//  Aux2(x) =<   1 / (1 + e^x)                      -3.7e+01 <  x <   +3.7e+01
//            |  e^-x                               +3.7e+01 <  x <   +7.5e+02
//            |  0                                  +7.5e+02 <= x
//            \
//
//            /
//            |  0                                              x <=  -7.5e+02
// d          |  - e^x                              -7.5e+02 <  x <   -3.7e+01
// --Aux2(x)=<  - e^x / {1 + e^x}^2                 -3.7e+01 <= x <=  +3.7e+01
// dx         |  - e^-x                             +3.7e+01 <  x <   +7.5e+02
//            |  0                                  +7.5e+02 <= x
//            \
//
//  Additional comment, ERK.  The auxilliary functions (aux1 and aux2) are
//  the only functions that are actually required.  The Bernouli functions
//  appear early in Selberherr's derivation of the electron and hole current
//  densities, but are not used in the final expressions.
//
//  Calculating the numerical limits for the Bernouli functions has proved
//  to be problematic on some platforms, so in this class only the aux1
//  and aux2 limits will be calculated.
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 08/01/04
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_BernouliSupport_h
#define Xyce_N_DEV_BernouliSupport_h

// ---------- Standard Includes ----------
#include <string>

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_UTL_Param.h>
#include <N_UTL_MachDepParams.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {

typedef double (* FUNC)(double);

//-----------------------------------------------------------------------------
// Class         : BernouliSupport
//
// Purpose       : This class contains support for functions required by the 
//                 Scharfetter-Gummel discretization.  That discretization, as 
//                 described in Selberherr's book, "Analysis and Simulation of 
//                 Semiconductor Devices", describes the discretization in terms 
//                 of Bernouli functions.   Hence, the name of this file and 
//                 class.
//
//                 However, the final form of the discretization doesn't 
//                 directly have Bernouli functions in it.  Instead, the 
//                 functions that are used are referred to as "aux1" and 
//                 "aux2".  For this reason, breakpoints associated with Bernouli 
//                 functions are not computed.  Only breakpoints associated with 
//                 aux1 and aux2 are computed.
//
//                 Aux1 and Aux2 are defined above.  For more details see 
//                 Selberherr's book.
//
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class BernouliSupport
{

public:
  BernouliSupport (bool regenerate = false);
  BernouliSupport (const BernouliSupport & right);
  virtual ~BernouliSupport ();

private:
  int sign(double x);
  double Bisection(FUNC func1, FUNC func2, double Xpos, double Xneg);
  double Secant(FUNC func1, FUNC func2, double x1);
  double Asymptotic(FUNC func1, FUNC func2, double x, double dx);

public:

  // aux1 function breakpoints
  double bp0_AUX1;
  double bp1_AUX1;
  double bp0_DAUX1;
  double bp1_DAUX1;

  // aux2 function breakpoints
  double bp0_AUX2;
  double bp1_AUX2;
  double bp2_AUX2;
  double bp0_DAUX2;
  double bp1_DAUX2;
  double bp2_DAUX2;
  double bp3_DAUX2;

  // Miscellaneous breakpoints
  double bp0_MISC;
};

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_BernouliSupport_h

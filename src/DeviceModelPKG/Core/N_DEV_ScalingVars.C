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
// Purpose        : This file contains a lot of the
//                  implementation of the model class for the two
//                  dimensional PDE based semiconductor device.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/05/03
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ----------  Standard Includes ----------
#include <iostream>

#include <N_UTL_fwd.h>

#include <N_DEV_ScalingVars.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : ScalingVars::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/15/09
//-----------------------------------------------------------------------------
std::ostream &operator<<(std::ostream & os, const ScalingVars & scaleVars)
{
  os << "\n\n-----------------------------------------" << std::endl;
  os << "\tPDE Scaling Vars:\n";
  //os << "\t\tdefad                 = " << scaleVars. <<"\n";

  os << "	x0   = " << scaleVars.x0  << "\n";  // distance scaling (cm)
  os << "	a0   = " << scaleVars.a0  << "\n";  // area scaling (cm^2)
  os << "	T0   = " << scaleVars.T0  << "\n";  // temperature scaling (K)
  os << "	V0   = " << scaleVars.V0  << "\n";  // electrostatic potential scaling (V)
  os << "	rV0  = " << scaleVars.rV0 << "\n"; // reciprocal of V0
  os << "	C0   = " << scaleVars.C0  << "\n";  // concentration scaling (cm^-3)
  os << "	D0   = " << scaleVars.D0  << "\n";  // diffusion coefficient scaling (cm^2/s)
  os << "	u0   = " << scaleVars.u0  << "\n";  // mobility coefficient scaling (cm^2/V/s)
  os << "	R0   = " << scaleVars.R0  << "\n";  // recombination rate scaling (cm^-3/s)
  os << "	rR0  = " << scaleVars.rR0 << "\n"; // reciprocal of R0
  os << "	t0   = " << scaleVars.t0  << "\n";  // time scaling (s)
  os << "	E0   = " << scaleVars.E0  << "\n";  // electric field scaling (V/s)
  os << "	F0   = " << scaleVars.F0  << "\n";  // particle flux scaling (cm^-2/s)
  os << "	J0   = " << scaleVars.J0  << "\n";  // current density scaling (A/cm^2)
  os << "	L0   = " << scaleVars.L0  << "\n";  // Laplacian scaling constant

  os << "	k0   = " << scaleVars.k0  << "\n";
  os << "	rt0  = " << scaleVars.rt0 << "\n";
  os << "	rk0  = " << scaleVars.rk0 << "\n";

  os << Xyce::section_divider << std::endl;
  os << std::endl;

  return os;
}

} // namespace Device
} // namespace Xyce

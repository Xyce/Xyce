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

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/04/08
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_ScalingVars_h
#define Xyce_ScalingVars_h

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : ScalingVars
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/04/08
//-----------------------------------------------------------------------------
class ScalingVars
{
public:
  ScalingVars()
    : x0(1.0), a0(1.0), T0(1.0), V0(1.0),
      rV0(1.0), C0(1.0), D0(1.0), u0(1.0),
      R0(1.0), rR0(1.0), t0(1.0), E0(1.0),
      F0(1.0), J0(1.0), L0(1.0), k0(1.0),
      rt0(1.0), rk0(1.0)
  {}

public:
  double x0;  // distance scaling (cm)
  double a0;  // area scaling (cm^2)
  double T0;  // temperature scaling (K)
  double V0;  // electrostatic potential scaling (V)
  double rV0; // reciprocal of V0
  double C0;  // concentration scaling (cm^-3);
  double D0;  // diffusion coefficient scaling (cm^2/s)
  double u0;  // mobility coefficient scaling (cm^2/V/s)
  double R0;  // recombination rate scaling (cm^-3/s)
  double rR0; // reciprocal of R0
  double t0;  // time scaling (s)
  double E0;  // electric field scaling (V/s)
  double F0;  // particle flux scaling (cm^-2/s)
  double J0;  // current density scaling (A/cm^2)
  double L0;  // Laplacian scaling constant

  double k0;
  double rt0;
  double rk0;
};

//-----------------------------------------------------------------------------
// Function      : ScalingVars::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/15/09
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const ScalingVars & scaleVars);

} // namespace Device
} // namespace Xyce

#endif


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
// Purpose        : This is a simple container class for different convergence
//                  status "return codes".
//
// Special Notes  : This point of having a class like this is to make the
//                  code more flexible in terms of what is considered to be
//                  a successful Newton solve.
//
//                  In the class N_NLS_DampedNewton, the function
//                  "converged_" performs a series of tests to determine
//                  convergence status, and depending on the results of
//                  these tests, it returns a value.  If the value is
//                  positive, then the solve is considered to be converged,
//                  and if it is <=0, it isn't.  
//
//                  This convention is used both inside the damped newton
//                  class, and outside.  It is used inside when the solver
//                  is trying to determine for itself whether or not to
//                  continue with the solve, or exit.  It is used outside
//                  by the time integrator, and also by continuation loops
//                  in the 2-level solver, in determining if the step was
//                  successful.
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/27/03
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_NLS_ReturnCodes.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Function      : ReturnCodes::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/02/03
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const ReturnCodes & rc)
{
  os << "\n\n-----------------------------------------" << std::endl
     << "\tNonlinear Solver Return Codes:\n"
     << "\t\tnormTooSmall      = " << rc.normTooSmall << "\n"
     << "\t\tnormalConvergence = " << rc.normalConvergence << "\n"
     << "\t\tnearConvergence   = " << rc.nearConvergence << "\n"
     << "\t\tsmallUpdate       = " << rc.smallUpdate << "\n"
     << "\t\tnanFail           = " << rc.nanFail << "\n"
     << "\t\ttooManySteps      = " << rc.tooManySteps << "\n"
     << "\t\ttooManyTranSteps  = " << rc.tooManyTranSteps << "\n"
     << "\t\tupdateTooBig      = " << rc.updateTooBig << "\n"
     << "\t\tstalled           = " << rc.stalled << "\n"
     << "\t\twrmsExactZero     = " << rc.wrmsExactZero << "\n"
     << section_divider << std::endl
     << std::endl;

  return os;
}

} // namespace Nonlinear
} // namespace Xyce

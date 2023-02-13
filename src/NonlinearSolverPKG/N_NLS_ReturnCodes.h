//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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

#ifndef Xyce_N_NLS_ReturnCodes_h
#define Xyce_N_NLS_ReturnCodes_h

#include <iosfwd>

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Class         : ReturnCodes
// Purpose       : Container class for solver success/failure return codes.
//
// Special Notes : Any result that you want the code to perceive as a
//                 "failure" should be <= 0.  Any result you want to code
//                 to perceive as a "success" should be >0.
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/27/03
//-----------------------------------------------------------------------------
class ReturnCodes
{
public:
  ReturnCodes ()
    : normTooSmall      (1),
      normalConvergence (2),
      nearConvergence   (-3),   // (near convergence, but assume failure)
      smallUpdate       (4),
      nanFail           (-6),
      tooManySteps      (-1),
      updateTooBig      (-2),
      stalled           (-3),  // (near convergence, but fail anyway)
      innerSolveFailed  (-5),
      linearSolverFailed (-9)
  {}

public:
  int normTooSmall;        // default = 1
  int normalConvergence;   // default = 2
  int nearConvergence;     // default = 3
  int smallUpdate;         // default = 4
  int nanFail;             // default = -6
  int tooManySteps;        // default = -1 
  int updateTooBig;        // default = -2
  int stalled;             // default = -3;
  int innerSolveFailed;    // default = -5;
  int linearSolverFailed;  // default = -9;
};

std::ostream & operator<<(std::ostream & os, const ReturnCodes & rc);

} // namespace Nonlinear
} // namespace Xyce

#endif // Xyce_N_NLS_ReturnCodes_h


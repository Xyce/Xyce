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
// Purpose        : Simple Direct Linear Solver Interface
//
// Special Notes  : This direct solver is used for trivial 1x1 linear systems
//
// Creator        : Heidi Thornquist, SNL
//
// Creation Date  : 03/07/2013
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_SimpleSolver_h
#define Xyce_N_LAS_SimpleSolver_h

#include <N_LAS_Solver.h>
#include <N_UTL_fwd.h>
#include <N_LAS_fwd.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : SimpleSolver
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
class SimpleSolver : public Solver
{

public:
  // Constructor
  SimpleSolver(
    Problem &     problem,
    Util::OptionBlock & options);

  // Destructor
  ~SimpleSolver();

  // Set the solver options
  bool setOptions(const Util::OptionBlock & OB);

  // Solve function: x = A^(-1) b.
  // This class is only used when A is a 1x1 matrix so x = b / A(1,1)
  int doSolve( bool reuse_factors, bool transpose = false );

private:

  //Options
  Util::OptionBlock * options_;

  //Timer
  Util::Timer * timer_;

};

} // namespace Linear
} // namespace Xyce

#endif // Xyce_N_LAS_SimpleSolver_h

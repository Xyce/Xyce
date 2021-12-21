//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Purpose        : Abstract interface to linear solver type.
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/17/04
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_Solver_h
#define Xyce_N_LAS_Solver_h

#include <N_LAS_fwd.h>
#include <N_UTL_fwd.h>

#include <N_UTL_Stats.h>

#include <Teuchos_RCP.hpp>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : Solver
// Purpose       : Abstract interface to linear solver type.
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/17/04
//-----------------------------------------------------------------------------
class Solver
{

public:
  //Constructors
  Solver(Problem& prob, bool isIterative)
  : solutionTime_(0.0),
    lasProblem_(prob),
    isIterative_(isIterative)
  {}

  //Destructor
  virtual ~Solver() {}

  // Set the solver options
  virtual bool setOptions(const Util::OptionBlock & OB) = 0;
  virtual bool setDefaultOptions() { return true; }
  virtual bool setPreconditioner( const Teuchos::RCP<Preconditioner>& precond ) { return true; }

  // Set individual options
  virtual bool setParam( const Util::Param & param ) { return true; }

  // Get info such as Num Iterations, Residual, etc.
  virtual bool getInfo( Util::Param & info ) { return true; }

  //Residual
  virtual double residual() { return 0.0; }
  virtual void setTolerance( const double & tol ) {}

  // Solve function: x = A^(-1) b.
  // input parameter 'ReuseFactors': If 'true', do not factor A, rather reuse
  // factors from previous solve.  Useful for inexact nonlinear techniques and
  // multiple RHS solves.
  virtual int doSolve( bool reuse_factors, bool transpose = false ) = 0;

  int solve( bool reuse_factors = false )
  {
    Stats::StatTop _linearSolveStat("Linear Solve");
    Xyce::Stats::TimeBlock _linearSolveTimer(_linearSolveStat);

    return doSolve(reuse_factors, false);
  }

  int solveTranspose( bool reuse_factors = false )
  {
    Stats::StatTop _linearSolveStat("Linear Solve Transpose");
    Xyce::Stats::TimeBlock _linearSolveTimer(_linearSolveStat);

    return doSolve(reuse_factors, true);
  }
  
  const Problem& getProblem() { return lasProblem_; }

  double solutionTime() { return solutionTime_; }

  bool isIterative() { return isIterative_; }

protected:
  double solutionTime_;

  // Linear Problem
  Problem & lasProblem_;

private:
  const bool isIterative_;

  //No copying
  Solver(const Solver & right);
  Solver & operator=(const Solver & right);

  //No comparison
  bool operator==(const Solver & right) const;
  bool operator!=(const Solver & right) const;

};

} // namespace Linear
} // namespace Xyce

#endif

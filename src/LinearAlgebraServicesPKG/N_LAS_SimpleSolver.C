//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Purpose        : Simple direct solver when the matrix is 1x1
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL
//
// Creation Date  : 03/07/2013
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ---------- Xyce Includes ----------

#include <N_LAS_SimpleSolver.h>

#include <N_LAS_Problem.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Graph.h>

#include <N_UTL_fwd.h>
#include <N_UTL_Timer.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>

#include <N_PDS_Comm.h>

#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : SimpleSolver::SimpleSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
SimpleSolver::SimpleSolver(
  Problem &       prob,
  Util::OptionBlock &   options)
  : Solver(prob, false),
   options_( new Util::OptionBlock( options ) ),
   timer_( new Util::Timer())
{
  setOptions( options );
}

//-----------------------------------------------------------------------------
// Function      : SimpleSolver::~SimpleSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
SimpleSolver::~SimpleSolver()
{
  delete timer_;
  delete options_;
}

//-----------------------------------------------------------------------------
// Function      : SimpleSolver::setOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool SimpleSolver::setOptions( const Util::OptionBlock & OB )
{
  if( options_ ) delete options_;
  options_ = new Util::OptionBlock( OB );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : SimpleSolver::doSolve
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
int SimpleSolver::doSolve( bool reuse_factors, bool transpose )
{
  // Start the timer...
  timer_->resetStartTime();

  // Perform linear solve using factorization
  double begSolveTime = timer_->elapsedTime();

  // Make sure this is a trivial linear system.
  int NumGlobalRows = lasProblem_.getMatrix()->getNumRows();
  if (NumGlobalRows > 1) {
    // Inform user that a nontrivial matrix was found and linear solve has failed.
    Report::UserError0()
      << "Nontrivial matrix has been found, this cannot be handled by this linear solver!";
  }

  int NumIndices = 0;
  int Length = lasProblem_.getMatrix()->getGraph()->maxNumIndices();
  std::vector<int> Indices(Length);
  std::vector<double> Values(Length);
  int NumMyRows = lasProblem_.getMatrix()->getLocalNumRows();
  int localSingularMat = 0, singularMat = 0;

  for (int i=0; i<NumMyRows; ++i) {

    // Get ith row
    lasProblem_.getMatrix()->getLocalRowCopy(i, Length, NumIndices, &Values[0], &Indices[0]);
    if (NumIndices != 1) {
      // Inform user that an empty matrix was found and linear solve has failed.
    Report::UserError0()
      << "Empty matrix has been found, this linear solve has failed!";
    }

    double pivot = Values[0];
    if (pivot != 0.0)
      lasProblem_.getLHS()->update( 1.0/pivot, *(lasProblem_.getRHS()), 0.0 );
    else
      localSingularMat = 1;
  }

#ifdef Xyce_PARALLEL_MPI
  // Communicate singular matrix
  lasProblem_.getRHS()->pdsComm()->maxAll(&localSingularMat, &singularMat, 1);
#else
  singularMat = localSingularMat;
#endif

  if (singularMat) 
  {
    // Inform user that singular matrix was found and linear solve has failed.
    Report::UserWarning0()
      << "Numerically singular matrix found, returning zero solution to nonlinear solver!";

    return 1;  // see bug 414
  }

  if (VERBOSE_LINEAR)
  {
    double endSolveTime = timer_->elapsedTime();
    Xyce::lout() << "  Simple (1x1 Matrix) Solve Time: " << (endSolveTime-begSolveTime) << std::endl;
  }

  // Update the total solution time
  solutionTime_ = timer_->elapsedTime();

  if (VERBOSE_LINEAR)
    Xyce::lout() << "Total Linear Solution Time (Simple): " << solutionTime_ << std::endl;

  return 0;
}

} // namespace Linear
} // namespace Xyce

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

#include <Epetra_LinearProblem.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Export.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

// ---------- Xyce Includes ----------

#include <N_LAS_SimpleSolver.h>

#include <N_LAS_Problem.h>

#include <N_LAS_TransformTool.h>

#include <N_UTL_fwd.h>
#include <N_UTL_Timer.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>

#include <N_ERH_ErrorMgr.h>

#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_BlockMapOut.h>

#include <Teuchos_Utils.hpp>

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
  : Solver(false),
   lasProblem_(prob),
   problem_(prob.epetraObj()),
   outputLS_(0),
   outputBaseLS_(0),
   outputFailedLS_(0),
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
  for( Util::ParamList::const_iterator it_tpL = OB.begin();
         it_tpL != OB.end(); ++it_tpL )
  {
    std::string tag = it_tpL->uTag();

    if( tag == "OUTPUT_LS" ) outputLS_ = it_tpL->getImmutableValue<int>();

    if( tag == "OUTPUT_BASE_LS" ) outputBaseLS_ = it_tpL->getImmutableValue<int>();

    if( tag == "OUTPUT_FAILED_LS" ) outputFailedLS_ = it_tpL->getImmutableValue<int>();
  }

  if( options_ ) delete options_;
  options_ = new Util::OptionBlock( OB );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : SimpleSolver::setDefaultOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool SimpleSolver::setDefaultOptions()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : SimpleSolver::setParam
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool SimpleSolver::setParam( const Util::Param & param )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : SimpleSolver::getInfo
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool SimpleSolver::getInfo( Util::Param & info )
{
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

  Epetra_LinearProblem * prob = &problem_;

  // Output the linear system to a Matrix Market file every outputLS_ calls if outputLS_ > 0 (or outputBaseLS_)
  static int failure_number = 0, file_number = 1, base_file_number = 1;
  if (outputBaseLS_ || outputLS_) {
    if (!(base_file_number % outputBaseLS_) || !(file_number % outputLS_)) {
      char file_name[40];
      if (!reuse_factors) {
        if (base_file_number == 1) {
          EpetraExt::BlockMapToMatrixMarketFile( "Simple_BlockMap.mm", (problem_.GetMatrix())->Map() );
        }
        sprintf( file_name, "Simple_Matrix%d.mm", base_file_number );

        std::string sandiaReq = "Sandia National Laboratories is a multimission laboratory managed and operated by National Technology and\n%";
        sandiaReq += " Engineering Solutions of Sandia LLC, a wholly owned subsidiary of Honeywell International Inc. for the\n%";
        sandiaReq += " U.S. Department of Energy’s National Nuclear Security Administration under contract DE-NA0003525.\n%\n% Xyce circuit matrix.\n%%";

        EpetraExt::RowMatrixToMatrixMarketFile( file_name, *(problem_.GetMatrix()), sandiaReq.c_str() );
        sprintf( file_name, "Simple_RHS%d.mm", base_file_number );
        EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(problem_.GetRHS()) );
      }
    }
    // base_file_number++;  This will be incremented after the solution vector is written to file.
  }

  // Perform linear solve using factorization
  double begSolveTime = timer_->elapsedTime();

  // Make sure this is a trivial linear system.
  int NumGlobalRows = problem_.GetMatrix()->NumGlobalRows();
  if (NumGlobalRows > 1) {
    // Inform user that a nontrivial matrix was found and linear solve has failed.
    Report::UserError0()
      << "Nontrivial matrix has been found, this cannot be handled by this linear solver!";
  }

  int NumIndices = 0;
  int Length = problem_.GetMatrix()->MaxNumEntries();
  std::vector<int> Indices(Length);
  std::vector<double> Values(Length);
  int NumMyRows = problem_.GetMatrix()->NumMyRows();
  int localSingularMat = 0, singularMat = 0;

  for (int i=0; i<NumMyRows; ++i) {

    // Get ith row
    EPETRA_CHK_ERR(problem_.GetMatrix()->ExtractMyRowCopy(i, Length, NumIndices, &Values[0], &Indices[0]));
    if (NumIndices != 1) {
      // Inform user that an empty matrix was found and linear solve has failed.
    Report::UserError0()
      << "Empty matrix has been found, this linear solve has failed!";
    }

    double pivot = Values[0];
    if (pivot != 0.0)
      problem_.GetLHS()->Scale( 1.0/pivot, *(problem_.GetRHS()) );
    else
      localSingularMat = true;
  }

#ifdef Xyce_PARALLEL_MPI
  // Communicate singular matrix
  problem_.GetRHS()->Comm().GatherAll(&localSingularMat, &singularMat, 1);
#else
  singularMat = localSingularMat;
#endif

  if (singularMat) {

    // Inform user that singular matrix was found and linear solve has failed.
    Report::UserWarning0()
      << "Numerically singular matrix found, returning zero solution to nonlinear solver!";

    // Put zeros in the solution since Amesos was not able to solve this problem
    problem_.GetLHS()->PutScalar( 0.0 );
    // Output the singular linear system to a Matrix Market file if outputFailedLS_ > 0
    if (outputFailedLS_) {
      failure_number++;
      char file_name[40];
      if (failure_number== 1) {
        EpetraExt::BlockMapToMatrixMarketFile( "Failed_BlockMap.mm", (prob->GetMatrix())->Map() );
      }
      sprintf( file_name, "Failed_Matrix%d.mm", failure_number );
      std::string sandiaReq = "Sandia National Laboratories is a multimission laboratory managed and operated by National Technology and\n%";
      sandiaReq += " Engineering Solutions of Sandia LLC, a wholly owned subsidiary of Honeywell International Inc. for the\n%";
      sandiaReq += " U.S. Department of Energy’s National Nuclear Security Administration under contract DE-NA0003525.\n%\n% Xyce circuit matrix.\n%%";

      EpetraExt::RowMatrixToMatrixMarketFile( file_name, *(prob->GetMatrix()), sandiaReq.c_str() );
      sprintf( file_name, "Failed_RHS%d.mm", failure_number );
      EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(prob->GetRHS()) );
    }
    return 1;  // see bug 414
  }

  if (VERBOSE_LINEAR)
  {
    double endSolveTime = timer_->elapsedTime();
    Xyce::lout() << "  Simple (1x1 Matrix) Solve Time: " << (endSolveTime-begSolveTime) << std::endl;
  }

  // Output computed solution vectors, if requested.
  if (outputBaseLS_ || outputLS_) {
    if (!(base_file_number % outputBaseLS_) || !(file_number % outputLS_)) {
      char file_name[40];
      sprintf( file_name, "Simple_Soln%d.mm", base_file_number );
      EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(problem_.GetLHS()) );
    }
    base_file_number++; file_number++;
  }

  // Update the total solution time
  solutionTime_ = timer_->elapsedTime();

  if (VERBOSE_LINEAR)
    Xyce::lout() << "Total Linear Solution Time (Simple): " << solutionTime_ << std::endl;

  return 0;
}

} // namespace Linear
} // namespace Xyce

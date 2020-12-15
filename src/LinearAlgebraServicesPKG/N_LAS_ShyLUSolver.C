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
// Purpose        : Implementation file for the ShyLU linear solver interface.
//
// Special Notes  :
//
// Creator        : Heidi K. Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 09/25/2007
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------



// ----------   Xyce Includes   ----------


#include <N_LAS_ShyLUSolver.h>

#include <N_UTL_OptionBlock.h>
#include <N_UTL_FeatureTest.h>

#include <N_PDS_ParMap.h>

#include <N_UTL_Timer.h>

#include <N_ERH_ErrorMgr.h>

#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>

#include <N_LAS_TransformTool.h>
#include <N_LAS_Problem.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Utils.hpp>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : ShyLUSolver::~ShyLUSolver
// Purpose       : Default destructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
ShyLUSolver::~ShyLUSolver()
{
}

//-----------------------------------------------------------------------------
// Function      : ShyLUSolver::ShyLUSolver
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
ShyLUSolver::ShyLUSolver( Problem & problem,
                          Util::OptionBlock & options )
: Solver(true),
  symmetryDefault_(1),
  innerMaxIterDefault_(30),
  innerTolDefault_(1.0e-12),
  relThreshDefault_(1e-3),
  diagFactorDefault_(0.05),
  outerSolverDefault_("Belos"),
  separatorTypeDefault_("Wide"),
  schurSolverDefault_("AztecOO-Exact"),  // "AztecOO-Inexact", "Amesos", "Guided Probing"
  schurApproxTypeDefault_("Threshold"),
  outputLS_(0),
  outputBaseLS_(0),
  lasProblem_(problem),
  problem_(&(problem.epetraObj())),
  updatedParams_(false),
  linearResidual_(1.0),
  tProblem_(0)
{
  options_ = Teuchos::rcp( new Util::OptionBlock( options ) );
  timer_ = Teuchos::rcp( new Util::Timer() );

  shyluParams_ = Teuchos::rcp( new Teuchos::ParameterList() );
  setDefaultOptions();
  setOptions( *options_ );
}

//-----------------------------------------------------------------------------
// Function      : ShyLUSolver::setDefaultOptions
// Purpose       : resets Aztec options
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 9/25/07
//-----------------------------------------------------------------------------
bool ShyLUSolver::setDefaultOptions()
{
  symmetry_ = symmetryDefault_;
  innerMaxIter_ = innerMaxIterDefault_;
  innerTol_ = innerTolDefault_;
  relThresh_ = relThreshDefault_;
  diagFactor_ = diagFactorDefault_;
  outerSolver_ = outerSolverDefault_;
  separatorType_ = separatorTypeDefault_;
  schurSolver_ = schurSolverDefault_;
  schurApproxType_ = schurApproxTypeDefault_;
  updatedParams_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : ShyLUSolver::setOptions
// Purpose       : sets Aztec options and params from modelblock
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool ShyLUSolver::setOptions(const Util::OptionBlock& OB)
{
  shyluParams_->set("Symmetry", symmetry_); 
  shyluParams_->set("Inner Solver MaxIters", innerMaxIter_); 
  shyluParams_->set("Inner Solver Tolerance", innerTol_);
  shyluParams_->set("Relative Threshold", relThresh_);
  shyluParams_->set("Diagonal Factor", diagFactor_); 
  shyluParams_->set("Outer Solver Library", outerSolver_); 
  shyluParams_->set("Separator Type", separatorType_); 
  shyluParams_->set("Schur Approximation Method", schurApproxType_); 
  shyluParams_->set("Schur Complement Solver", schurSolver_); 
  //shyluParams_->set("Schur Recompute Iteration", 5);

  std::list<Util::Param>::const_iterator it_tpL = OB.begin();
  std::list<Util::Param>::const_iterator end_tpL = OB.end();
  for (; it_tpL != end_tpL; ++it_tpL)
  {
    setParam( *it_tpL );
  }
  
  // store for restart of solver_
  if( &OB != &*options_ )
  {
    options_ = Teuchos::rcp( new Util::OptionBlock(OB) );
  }

  if (!lasProblem_.matrixFree()) 
  {
    // create the transformation object if needed.
    if( Teuchos::is_null(transform_) ) transform_ = TransformTool()( OB );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ShyLUSolver::setParam
// Purpose       : sets Aztec option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool ShyLUSolver::setParam( const Util::Param & param )
{
  std::string tag = param.tag();
  std::string uTag = param.uTag();

  // Set our copies of these parameters.
  if( tag == "AZ_max_iter" )
    setInnerMaxIter(param.getImmutableValue<int>());
  else if( tag == "AZ_tol" )
    setInnerTol(param.getImmutableValue<double>());
  else if( tag == "ShyLU_rthresh" )
    setRelThresh(param.getImmutableValue<double>());
  else if( uTag == "OUTPUT_LS" )
    outputLS_ = param.getImmutableValue<int>();
  else if( uTag == "OUTPUT_BASE_LS" )
    outputBaseLS_ = param.getImmutableValue<int>();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ShyLUSolver::getInfo
// Purpose       : sets Aztec option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool ShyLUSolver::getInfo( Util::Param & param )
{
  if( param.tag() == "AZ_max_iter" )
    param.setVal( innerMaxIter_ );
  else if( param.tag() == "Iterations" )
    param.setVal( (int)numLinearIters_ );
  else if( param.tag() == "AZ_tol" )
    param.setVal( innerTol_ );
  else 
    return false;
  
  return true;
}

//-----------------------------------------------------------------------------
// Function      : ShyLUSolver::setShyLUOption
// Purpose       : sets ShyLU option
// Special Notes : Takes a string as the option identifier
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Computational Sciences
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool ShyLUSolver::setShyLUOption_(const char * paramName,
                                            const int val)
{
  return setShyLUCntl_( Util::Param(paramName, val) );
}

//-----------------------------------------------------------------------------
// Function      : ShyLUSolver::setShyLUParam
// Purpose       : sets ShyLU parameter
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool ShyLUSolver::setShyLUParam_(const char * paramName,
                                         const double val)
{
  return setShyLUCntl_( Util::Param(paramName, val) );
}

//-----------------------------------------------------------------------------
// Function      : ShyLUSolver::doSolve
// Purpose       : Calls the actual solver to solve Ax=b.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
int ShyLUSolver::doSolve( bool reuse_factors, bool transpose )
{
  if ( transpose ) {
    // Inform user that a nontrivial matrix was found and linear solve has failed.
    Report::UserError0() << 
      "ShyLU linear solver does not support transpose solves at this time.";
  }

  // Start the timer...
  timer_->resetStartTime();
  double transTime=0.0;
  double time1 = 0.0, time2 = 0.0;

  if (VERBOSE_LINEAR)
    time1 = timer_->wallTime();

  // The Epetra_LinearProblem, prob, is the linear system being solved.
  // It will point to either the original linear system or transformed system.
  Epetra_LinearProblem * prob = problem_;

  if( !Teuchos::is_null(transform_) )
  {
    if( !tProblem_ )
      tProblem_ = &((*transform_)( *problem_ ));
    prob = tProblem_;
    transform_->fwd();
  }

  if (VERBOSE_LINEAR)
  {
    time2 = timer_->wallTime();
    transTime += (time2-time1);
  }

  if( Teuchos::is_null(solver_) ) {

    if (lasProblem_.matrixFree()) {
      Xyce::Report::DevelFatal0().in("ShyLUSolver::solve()") << "cannot work on matrix-free linear systems!";
    }

    Epetra_CrsMatrix * epetraA = dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix());
    solver_ = Teuchos::rcp( new Ifpack_ShyLU( epetraA ) );

    // Set the parameters
    IFPACK_CHK_ERR(solver_->SetParameters(*shyluParams_));

    // Compute symbolic factorization stage of preconditioner.
    IFPACK_CHK_ERR(solver_->Initialize());

    if (VERBOSE_LINEAR)
    {
      time1 = timer_->wallTime();
      Xyce::dout() << "ShyLU Initialize Time: " << (time1 - time2) << std::endl;
    }
    
  }
      
  // Output the linear system to a Matrix Market file every outputLS_ calls if outputLS_ > 0
  static int file_number = 1, base_file_number = 1;
  if (outputLS_ && !lasProblem_.matrixFree()) {
    if (!(file_number % outputLS_)) {
      if (!reuse_factors) {
        Xyce::Linear::writeToFile( *prob, "Transformed", file_number, (file_number == 1) );
      }
    }
    file_number++;
  }
  if (outputBaseLS_ && !lasProblem_.matrixFree()) {
    if (!(base_file_number % outputBaseLS_)) {
      if (!reuse_factors) {
        Xyce::Linear::writeToFile( *problem_, "Base", base_file_number, (base_file_number == 1) );
      }
    }
    base_file_number++;
  }

  if (VERBOSE_LINEAR)
    time1 = timer_->wallTime();

  // Compute the preconditioner
  IFPACK_CHK_ERR(solver_->Compute());
 
  if (VERBOSE_LINEAR)
  { 
    time2 = timer_->wallTime();
    Xyce::dout() << "ShyLU Compute Time: " << (time2 - time1) << std::endl;
  }

  // Solve the linear system
  int solveRet = solver_->ApplyInverse( *prob->GetRHS(), *problem_->GetLHS() );

  if (solveRet)
    Xyce::Report::DevelFatal0().in("ShyLUSolver::solve()") << "ShyLU solver could not be applied!";
  
  //numLinearIters_ = solver_->getNumIters();

  if (VERBOSE_LINEAR)
  {
    time1 = timer_->wallTime();
    Xyce::dout() << "ShyLU Solve Time: " << (time1 - time2) << std::endl;
  }

  if (DEBUG_LINEAR)
  {
    Epetra_MultiVector* b = prob->GetRHS();
    int numrhs = b->NumVectors();
    std::vector<double> actual_resids( numrhs ), rhs_norm( numrhs );
    Epetra_MultiVector resid( b->Map(), numrhs );
    prob->GetOperator()->Apply( *(prob->GetLHS()), resid );
    resid.Update( -1.0, *b, 1.0 );
    resid.Norm2( &actual_resids[0] );
    b->Norm2( &rhs_norm[0] );
    for (int i=0; i<numrhs; i++ ) {
      Xyce::lout() << "Problem " << i << " : \t" <<(actual_resids[i]/rhs_norm[i]) << std::endl;
    }
  }

  if( !Teuchos::is_null(transform_) ) transform_->rvs();

  if (VERBOSE_LINEAR)
  {
    time2 = timer_->wallTime();
    transTime += (time2-time1);
  }

  // Update the total solution time
  solutionTime_ = timer_->elapsedTime();
  if (VERBOSE_LINEAR)
  {
    Xyce::dout() << "Transform Time: " << transTime << std::endl;
    Xyce::dout() << "Total Solve Time: " << solutionTime_ << std::endl;
  }

  linearResidual_ = innerTol_/10;
  return 0;
}

} // namespace Linear
} // namespace Xyce

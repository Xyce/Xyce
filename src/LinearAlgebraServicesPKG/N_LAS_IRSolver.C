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

//-------------------------------------------------------------------------
//
// Purpose        : Iterative refinement solver wrapper
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL
//
// Creation Date  : 05/20/04
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------

#include <Amesos.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Export.h>
#include <Epetra_Map.h>

// ---------- Xyce Includes ----------

#include <N_UTL_fwd.h>

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_IRSolver.h>
#include <N_LAS_Problem.h>
#include <N_LAS_EpetraProblem.h>
#include <N_LAS_EpetraHelpers.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_TransformTool.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Timer.h>
#include <N_UTL_MachDepParams.h>

#include <Teuchos_Utils.hpp>

namespace Xyce {
namespace Linear {

//Solver defaults.
const std::string IRSolver::type_default_ = "KLU";
const double IRSolver::tol_default_ = 1e-9;

//-----------------------------------------------------------------------------
// Function      : IRSolver::IRSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
IRSolver::IRSolver(
  Problem &       problem,
  Util::OptionBlock &   options)
  : Solver(problem, false),
    type_(type_default_),
    ir_tol_(tol_default_),
    solver_(0),
    repivot_(true),
    outputLS_(0),
    outputBaseLS_(0),
    outputFailedLS_(0),
    tProblem_(0),
    options_( new Util::OptionBlock( options ) ),
    timer_( new Util::Timer() )
{
  EpetraProblem& eprob = dynamic_cast<EpetraProblem&>(lasProblem_);
  problem_ = &(eprob.epetraObj()); 

  setOptions( options );
}

//-----------------------------------------------------------------------------
// Function      : IRSolver::~IRSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
IRSolver::~IRSolver()
{
  delete solver_;
  delete timer_;
  delete options_;
}

//-----------------------------------------------------------------------------
// Function      : IRSolver::setOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool IRSolver::setOptions( const Util::OptionBlock & OB )
{
  bool foundAMD = false, foundPartition = false, foundSingleton = false;

  for( Util::ParamList::const_iterator it_tpL = OB.begin();
         it_tpL != OB.end(); ++it_tpL )
  {
    std::string tag = it_tpL->uTag();
    if( tag == "IR_SOLVER_TYPE" ) type_ = it_tpL->usVal();

    if( tag == "IR_SOLVER_TOL" ) ir_tol_ = it_tpL->getImmutableValue<double>();

    if( tag == "KLU_REPIVOT" ) repivot_ = static_cast<bool>(it_tpL->getImmutableValue<int>());
    
    if( tag == "OUTPUT_LS" ) outputLS_ = it_tpL->getImmutableValue<int>();

    if( tag == "OUTPUT_BASE_LS" ) outputBaseLS_ = it_tpL->getImmutableValue<int>();

    if( tag == "OUTPUT_FAILED_LS" ) outputFailedLS_ = it_tpL->getImmutableValue<int>();

    if( tag == "TR_AMD" ) foundAMD = true;

    if( tag == "TR_PARTITION" ) foundPartition = true;

    if( tag == "TR_SINGLETON_FILTER" ) foundSingleton = true;
  }

  if( options_ ) delete options_;
  options_ = new Util::OptionBlock( OB );

  // Unless the user set these options, they should not be used in the transforms for serial
  // or parallel runs.
  if ( !foundAMD )
    options_->addParam(Util::Param("TR_amd", 0));
  if ( !foundPartition )
    options_->addParam(Util::Param("TR_partition", 0));
  if ( !foundSingleton )
    options_->addParam(Util::Param("TR_singleton_filter", 0));

#ifdef Xyce_PARALLEL_MPI
  options_->addParam(Util::Param("TR_reindex", 1));
#endif

  if( Teuchos::is_null(transform_) ) transform_ = TransformTool()( *options_ );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : IRSolver::setDefaultOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool IRSolver::setDefaultOptions()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : IRSolver::setParam
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool IRSolver::setParam( const Util::Param & param )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : IRSolver::getInfo
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool IRSolver::getInfo( Util::Param & info )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : IRSolver::doSolve
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
int IRSolver::doSolve( bool reuse_factors, bool transpose )
{
  // Start the timer...
  timer_->resetStartTime();

  int linearStatus = 0;

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

  // Output the linear system to a Matrix Market file every outputLS_ calls if outputLS_ > 0
  static int failure_number = 0, file_number = 1, base_file_number = 1;
  if (outputLS_) {
    if (!(file_number % outputLS_)) {
      if (!reuse_factors) {
        Xyce::Linear::writeToFile( *prob, "Transformed", file_number, (file_number == 1) );
      }
    }
    // file_number++;  This will be incremented after the solution vector is written to file.
  }
  if (outputBaseLS_) {
    if (!(base_file_number % outputBaseLS_)) {
      if (!reuse_factors) {
        Xyce::Linear::writeToFile( *problem_, "Base", base_file_number, (base_file_number == 1) );
      }
    }
    // base_file_number++;  This will be incremented after the solution vector is written to file.
  }

  // Set the traceback mode in Epetra so it prints out warnings
  if (DEBUG_LINEAR)
    dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix())->SetTracebackMode( 2 );
  else
    dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix())->SetTracebackMode( 0 );


  Amesos localAmesosObject;
  if( !solver_ )
  {
    // the Query() function expects a string
    // in lower case with the first letter in upper case
    // So, our "KLU" must become "Klu"
    std::string solverType( type_ );
    if( type_ == "KLU" )
    {
      solverType = "Amesos_Klu";
    }
    else if( type_ == "SUPERLU" )
    {
      solverType = "Amesos_Superlu";
    }
    else if( type_ == "SUPERLUDIST" )
    {
      solverType = "Amesos_Superludist";
    }
    else if( type_ == "PARAKLETE" )
    {
      solverType = "Amesos_Paraklete";
    }
    else if( type_ == "PARDISO" )
    {
      solverType = "Amesos_Pardiso";
    }
    else if( type_ == "LAPACK" )
    {
      solverType = "Amesos_Lapack";
    }
    else if( type_ == "SCALAPACK" )
    {
      solverType = "Amesos_Scalapack";
    }
    else if( type_ == "MUMPS" )
    {
      solverType = "Amesos_Mumps";
    }

    if( !localAmesosObject.Query( solverType ) )
      Report::DevelFatal0() 
        << "Unknown or Unavailable Linear Solver: " << type_;

    solver_ = localAmesosObject.Create( solverType, *prob );

    Teuchos::ParameterList params;

#ifndef Xyce_PARALLEL_MPI
    // Inform solver not to check inputs to reduce overhead.
    params.set( "TrustMe", true );
    // If repivot == true (default), recompute the pivot order each numeric factorization,
    // else try to re-use pivot order to expedite numeric factorization.
    params.set( "Refactorize", !repivot_ );
#else
    if (type_ == "SUPERLUDIST") {
      Teuchos::ParameterList& sludistParams = params.sublist("Superludist");
      sludistParams.set("ReuseSymbolic", true );
    }
#endif

    if (VERBOSE_LINEAR)
      Xyce::dout() << "IRSolver::solve() setting solver : " << type_ << "\n"
                   << "IRSolver::solve() setting parameters : " << params << std::endl;

    solver_->SetParameters( params );

    double begSymTime = timer_->elapsedTime();

    // Perform symbolic factorization and check return value for failure
    linearStatus = solver_->SymbolicFactorization();
    if (linearStatus != 0)
    {
      // Update the total solution time
      solutionTime_ = timer_->elapsedTime();

      return linearStatus;
    }

    if (VERBOSE_LINEAR)
    {
      double endSymTime = timer_->elapsedTime();
      Xyce::dout() << "  Amesos (" << type_ << ") Symbolic Factorization Time: "
                   << (endSymTime - begSymTime) << std::endl;
    }
  }

  // Set the transpose flag only if that has changed since the last solve.
  if ( solver_->UseTranspose() != transpose )
  {
    solver_->SetUseTranspose( transpose );
  }

  // Try to solve with previous factors
  solver_->Solve();

  int numrhs = prob->GetLHS()->NumVectors();
  std::vector<double> resNorm(numrhs,0.0), bNorm(numrhs,0.0);
  Epetra_MultiVector res( prob->GetLHS()->Map(), prob->GetLHS()->NumVectors() );
  prob->GetOperator()->Apply( *(prob->GetLHS()), res );
  res.Update( 1.0, *(prob->GetRHS()), -1.0 );
  res.Norm2( &resNorm[0] );
  prob->GetRHS()->Norm2( &bNorm[0] );
  bool reuse_success = true;
  double max_residual = 0.0;
  double min_tol_ = 1e-2;
  int max_iters = 10;
  for (int i=0; i<numrhs; i++)
  {
    //Xyce::dout() << "IRSolver()::bNorm[" << i << "] = " << bNorm[i] << std::endl;
    double residual = resNorm[i];
    if (bNorm[i] > Util::MachineDependentParams::MachineEpsilon())
      residual /= bNorm[i];
    if (residual > ir_tol_)
    {
      reuse_success = false;
      //Xyce::dout() << "IRSolver()::Solve() numeric factors did not meet tolerance: " << residual << std::endl;
    }
    if (residual > max_residual)
      max_residual = residual;
  }

  if (reuse_success)
  {
    // Update the total solution time
    solutionTime_ = timer_->elapsedTime();

    if (VERBOSE_LINEAR)
      Xyce::dout() << "Total Linear Solution Time (Reused, Amesos " << type_ << "): "
                   << solutionTime_ << std::endl;

    return linearStatus;
  }
  else if (max_residual < min_tol_) 
  {
    // Save the original RHS
    Epetra_MultiVector x( prob->GetLHS()->Map(), numrhs );
    Epetra_MultiVector bsave( prob->GetLHS()->Map(), numrhs );
    bsave = *(prob->GetRHS());
    x = *(prob->GetLHS());

    // Try refinement with current factorization
    int iter = 0;
    bool diverging = false;
    //Xyce::dout() << "Trying iterative refinement " << max_residual << " and min_tol_ = " << min_tol_ << " !" << std::endl;
    while ((max_residual > ir_tol_)&&(iter < max_iters)&&(!diverging))
    {
      // Copy current residual to the right-hand side of the problem
      *(prob->GetRHS()) = res;

      // Try to solve with previous factors
      solver_->Solve();
 
      // Now update the soultion and check the residual 
      x.Update( 1.0, *(prob->GetLHS()), 1.0 );
      prob->GetOperator()->Apply( x, res );
      res.Update( 1.0, bsave, -1.0 );
      res.Norm2( &resNorm[0] );
 
      double new_max_residual = 0.0; 
      for (int i=0; i<numrhs; i++)
      {
        double residual = resNorm[i];
        if (bNorm[i] > Util::MachineDependentParams::MachineEpsilon())
          residual /= bNorm[i];
    
        if (residual > new_max_residual)
          new_max_residual = residual;
      }

      //Xyce::dout() << "IRSolver: Refinement step:  old residual " << max_residual << ", new residual " << new_max_residual << std::endl; 
      if (new_max_residual > max_residual)
        diverging = true;

      // Set the max residual to the residual from this refinement step
      max_residual = new_max_residual;

      iter++;
    }

    // Set the RHS to the original vector.
    *(prob->GetRHS()) = bsave;

    if (max_residual < ir_tol_)
    {
      // Set the LHS to the refined solution.
      *(prob->GetLHS()) = x;

      // Update the total solution time
      solutionTime_ = timer_->elapsedTime();

      prob->GetOperator()->Apply( *(prob->GetLHS()), res );
      res.Update( 1.0, *(prob->GetRHS()), -1.0 );
      res.Norm2( &resNorm[0] );
      prob->GetRHS()->Norm2( &bNorm[0] );
      Xyce::lout() << "Linear System Residual (AMESOS_" << type_ << "): " << std::endl;
      for (int i=0; i<numrhs; i++)
      {
        //Xyce::dout() << "IRSolver()::bNorm[" << i << "] = " << bNorm[i] << std::endl;
        if (bNorm[i] > Util::MachineDependentParams::MachineEpsilon())
          std::cout << "  Problem " << i << " : " << (resNorm[i]/bNorm[i]) << std::endl;
        else
          std::cout << "  Problem " << i << " : " << resNorm[i] << std::endl;
      } 
      if (VERBOSE_LINEAR)
        Xyce::dout() << "Total Linear Solution Time (Reused, Amesos " << type_ << "): "
                     << solutionTime_ << std::endl;
 

      return linearStatus;
    }
  }

  // Perform numeric factorization and check return value for failure
  if( !reuse_factors ) {

    double begNumTime = timer_->elapsedTime();

    linearStatus = solver_->NumericFactorization();
    if (VERBOSE_LINEAR)
    {
      double endNumTime = timer_->elapsedTime();
      Xyce::dout() << "  Amesos (" << type_ << ") Numeric Factorization Time: "
                   << (endNumTime - begNumTime) << std::endl;
    }
    
    if (linearStatus != 0) {

      // Inform user that singular matrix was found and linear solve has failed.
      Report::UserWarning0() 
        << "Numerically singular matrix found by Amesos, returning zero solution to nonlinear solver!";

      // Put zeros in the solution since Amesos was not able to solve this problem
      prob->GetLHS()->PutScalar( 0.0 );

      // Output the singular linear system to a Matrix Market file if outputFailedLS_ > 0
      if (outputFailedLS_) {
        failure_number++;
        Xyce::Linear::writeToFile( *prob, "Failed", failure_number, (failure_number == 1) );
      }

      // Update the total solution time
      solutionTime_ = timer_->elapsedTime();

      return linearStatus;  // return the actual status (see bug 414 SON)
    }
  }

  // Perform linear solve using factorization
  double begSolveTime = timer_->elapsedTime();

  solver_->Solve();

  if (VERBOSE_LINEAR)
  {
    double endSolveTime = timer_->elapsedTime();
    Xyce::dout() << "  Amesos (" << type_ << ") Solve Time: "
                 << (endSolveTime - begSolveTime) << std::endl;
  }

  if (DEBUG_LINEAR) 
  {
    int numrhs = prob->GetLHS()->NumVectors();
    std::vector<double> resNorm(numrhs,0.0), bNorm(numrhs,0.0);
    Epetra_MultiVector res( prob->GetLHS()->Map(), prob->GetLHS()->NumVectors() );
    prob->GetOperator()->Apply( *(prob->GetLHS()), res );
    res.Update( 1.0, *(prob->GetRHS()), -1.0 );
    res.Norm2( &resNorm[0] );
    prob->GetRHS()->Norm2( &bNorm[0] );
    Xyce::lout() << "Linear System Residual (AMESOS_" << type_ << "): " << std::endl;
    for (int i=0; i<numrhs; i++)
      std::cout << "  Problem " << i << " : " << (resNorm[i]/bNorm[i]) << std::endl;
  }

  if( !Teuchos::is_null(transform_) ) transform_->rvs();

  // Output computed solution vectors, if requested.
  if (outputLS_) {
    if (!(file_number % outputLS_)) {
      Teuchos::RCP<Problem> las_prob = Teuchos::rcp( new EpetraProblem( Teuchos::rcp( prob, false ) ) );
      char file_name[40];
      sprintf( file_name, "Transformed_Soln%d.mm", file_number );
      las_prob->getLHS()->writeToFile( file_name, false, true );
    }
    file_number++;
  }
  if (outputBaseLS_) {
    if (!(base_file_number % outputBaseLS_)) {
      char file_name[40];
      sprintf( file_name, "Base_Soln%d.mm", base_file_number );
      lasProblem_.getLHS()->writeToFile( file_name, false, true );
    }
    base_file_number++;
  }

  // Update the total solution time
  solutionTime_ = timer_->elapsedTime();

  if (VERBOSE_LINEAR)
    Xyce::dout() << "Total Linear Solution Time (Amesos " << type_ << "): "
                 << solutionTime_ << std::endl;

  return 0;
}

} // namespace Linear
} // namespace Xyce

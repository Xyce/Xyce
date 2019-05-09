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
// Purpose        : Implementation file for the Iterative linear solver
//                  interface.
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/20/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_LAS_AztecOOSolver.h>

#include <N_UTL_fwd.h>

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Preconditioner.h>
#include <N_LAS_Problem.h>
#include <N_LAS_TransformTool.h>
#include <N_LAS_TrilinosPrecondFactory.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Timer.h>

#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_BlockMapOut.h>

// ---------- Other Includes  -----------

#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>
#include <AztecOO.h>
#include <Teuchos_OrdinalTraits.hpp>

namespace Xyce {
namespace Linear {

int AztecOOSolver::num_AZ_options_ = 27;
const char * AztecOOSolver::AZ_options_[] =
        { "AZ_solver", "AZ_scaling", "AZ_precond", "AZ_conv", "AZ_output",
          "AZ_pre_calc", "AZ_max_iter", "AZ_poly_ord", "AZ_overlap",
          "AZ_type_overlap", "AZ_kspace", "AZ_orthog", "AZ_aux_vec",
          "AZ_reorder", "AZ_keep_info", "AZ_recursion_level", "AZ_print_freq",
          "AZ_graph_fill", "AZ_subdomain_solve", "AZ_init_guess",
          "AZ_keep_kvecs", "AZ_apply_kvecs", "AZ_orth_kvecs",
          "AZ_ignore_scaling", "AZ_check_update_size", "AZ_extreme",
          "AZ_diagnostics" };

int AztecOOSolver::num_AZ_params_ = 7;
const char * AztecOOSolver::AZ_params_[] =
        { "AZ_tol", "AZ_drop", "AZ_ilut_fill", "AZ_omega", "AZ_rthresh",
          "AZ_athresh", "AZ_weights" };


//-----------------------------------------------------------------------------
// Function      : AztecOOSolver::~AztecOOSolver
// Purpose       : Default destructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
AztecOOSolver::~AztecOOSolver()
{
  delete solver_;
  delete options_;
  delete timer_;
}

//-----------------------------------------------------------------------------
// Function      : AztecOOSolver::AztecOOSolver
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
AztecOOSolver::AztecOOSolver(
  Problem &       problem,
  Util::OptionBlock &   options)
  : Solver(true),
  reduceKSpace_(false),
  maxKSpace_(50),
  probDiff_(0),
  preCondDefault_(14),
  solverDefault_(1),
  scalingDefault_(0),
  subdomainSolveDefault_(9),
  convergenceDefault_(0),
  outputDefault_(VERBOSE_LINEAR ? 50 : AZ_none),
  diagnosticsDefault_(AZ_none),
  precalcDefault_(1),
  maxIterDefault_(200),
  KSpaceDefault_(50),
  reorderDefault_(0),
  keepInfoDefault_(1),
  orthogDefault_(1),
  overlapDefault_(0),
  toleranceDefault_(1.0e-9),
  dropDefault_(1.0e-03),
  ilutFillDefault_(2.0),
  rThreshDefault_(1.0001),
  aThreshDefault_(0.0001),
  linearResidual_(1.0),
  numLinearIters_(0),
  outputLS_(0),
  outputBaseLS_(0),
  lasProblem_(problem),
  problem_(&(problem.epetraObj())),
  options_( new Util::OptionBlock( options ) ),
  useAztecPrecond_(false),
  isPrecSet_(false),
  solver_(0),
  tProblem_(0),
  timer_( new Util::Timer())
{
  problem_->SetPDL((ProblemDifficultyLevel) probDiff_);

  setDefaultOptions();

  setOptions( *options_ );

  // Defaults for parameters above - see comment above for descriptions.
  // !!!NOTE: these values are based on the Aztec include files on or about
  // June, 2001.  We will eventually need a more consistent way of setting
  // these - SAH

}

//-----------------------------------------------------------------------------
// Function      : AztecOOSolver::setDefaultOptions
// Purpose       : resets Aztec options
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/04
//-----------------------------------------------------------------------------
bool AztecOOSolver::setDefaultOptions()
{
  // Set defaults
  resetPreCond();
  resetSolver();
  resetScaling();
  resetSubdomainSolve();
  resetConvergence();
  resetOutput();
  resetDiagnostics();
  resetPrecalc();
  resetMaxIter();
  resetOverlap();
  resetKSpace();
  resetReorder();
  resetKeepInfo();
  resetOrthog();
  resetTolerance();
  resetDrop();
  resetILUTFill();
  resetRThresh();
  resetAThresh();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AztecOOSolver::setOptions
// Purpose       : sets Aztec options and params from modelblock
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
bool AztecOOSolver::setOptions(const Util::OptionBlock & OB)
{
  if( solver_ )
  {
    setAztecOption_( "AZ_output", output_ );
    setAztecOption_( "AZ_diagnostics", diagnostics_ );
    setAztecOption_( "AZ_precond", preCond_ );
    setAztecOption_( "AZ_subdomain_solve", subdomainSolve_ );
    setAztecOption_( "AZ_overlap", overlap_ );
    setAztecParam_( "AZ_ilut_fill", ilutFill_ );
    setAztecParam_( "AZ_drop", drop_ );
    setAztecOption_( "AZ_kspace", KSpace_ );
    setAztecParam_( "AZ_athresh", aThresh_ );
    setAztecParam_( "AZ_rthresh", rThresh_ );
    setAztecParam_( "AZ_tol", tolerance_ );
    setAztecOption_( "AZ_max_iter", maxIter_ );

    Util::ParamList::const_iterator it_tpL = OB.begin();
    Util::ParamList::const_iterator end_tpL = OB.end();
    for (; it_tpL != end_tpL; ++it_tpL)
    {
      setParam( *it_tpL );
    }
  }

  // store for restart of solver_
  if( &OB != options_ )
  {
    if( options_ ) delete options_;
    options_ = new Util::OptionBlock(OB);
  }

  // if default Krylov subspace size is too large for AztecOO,
  // and user did not specify a smaller one, reset KSpace_ to maxKSpace_.
  if (reduceKSpace_)
  {
    if (KSpace_ > maxKSpace_)
    {
      Report::UserWarning0() 
        << "AztecOOSolver::solve():  Krylov subspace memory requirements too large, resizing Krylov subspace to " << maxKSpace_ << "!";
      setAztecOption_( "AZ_kspace", maxKSpace_ );
    }
  }

  // set singleton filtering as default for iterative solvers
  const Util::Param *param = Util::findParameter(OB.begin(), OB.end(), "TR_singleton_filter");
  if (!param)
  {
    options_->addParam(Util::Param("TR_singleton_filter", 1));
  }

  if (!lasProblem_.matrixFree())
  {
    if( Teuchos::is_null(transform_) ) transform_ = TransformTool()( *options_ );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AztecOOSolver::setParam
// Purpose       : sets Aztec option
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
bool AztecOOSolver::setParam( const Util::Param & param )
{
  std::string tag = param.tag();
  std::string uTag = param.uTag();

  // Set our copies of these parameters that get passed to the solver in the
  // "iterate" command
  if( tag == "AZ_max_iter" )
    setMaxIter(param.getImmutableValue<int>());
  else if( tag == "AZ_tol" )
    setTolerance(param.getImmutableValue<double>());
  else if( uTag == "USE_AZTEC_PRECOND" )
    useAztecPrecond_ = param.getImmutableValue<int>();
  else if( uTag == "OUTPUT_LS" )
    outputLS_ = param.getImmutableValue<int>();
  else if( uTag == "OUTPUT_BASE_LS" )
    outputBaseLS_ = param.getImmutableValue<int>();
  else
    setAztecCntl_( param );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AztecOOSolver::getInfo
// Purpose       : sets Aztec option
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
bool AztecOOSolver::getInfo( Util::Param & param )
{
  if     ( param.tag() == "AZ_precond" )
    param.setVal( preCond_ );
  else if( param.tag() == "AZ_subdomain_solve" )
    param.setVal( subdomainSolve_ );
  else if( param.tag() == "AZ_kspace" )
    param.setVal( KSpace_ );
  else if( param.tag() == "AZ_athresh" )
    param.setVal( aThresh_ );
  else if( param.tag() == "AZ_rthresh" )
    param.setVal( rThresh_ );
  else if( param.tag() == "AZ_ilut_fill" )
    param.setVal( ilutFill_ );
  else if( param.tag() == "AZ_drop" )
    param.setVal( drop_ );
  else if( param.tag() == "AZ_overlap" )
    param.setVal( overlap_ );
  else if( param.tag() == "AZ_output" )
    param.setVal( output_ );
  else if( param.tag() == "AZ_diagnostics" )
    param.setVal( diagnostics_ );
  else if( param.tag() == "AZ_max_iter" )
    param.setVal( maxIter_ );
  else if( param.tag() == "Iterations" )
    param.setVal( (int)numLinearIters_ );
  else if( param.tag() == "AZ_tol" )
    param.setVal( tolerance_ );
  else if( param.tag() == "use_aztec_precond" )
    param.setVal( useAztecPrecond_ );
  else
    return false;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AztecOOSolver::setAztecOption
// Purpose       : sets Aztec option
// Special Notes : Takes a string as the option identifier
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
bool AztecOOSolver::setAztecOption_(const char * paramName,
                                            const int val)
{
  return setAztecCntl_( Util::Param(paramName, val) );
}

//-----------------------------------------------------------------------------
// Function      : AztecOOSolver::setAztecParam
// Purpose       : sets Aztec parameter
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
bool AztecOOSolver::setAztecParam_(const char * paramName,
                                         const double val)
{
  return setAztecCntl_( Util::Param(paramName, val) );
}

//-----------------------------------------------------------------------------
// Function      : AztecOO::printParams_
// Purpose       : Print out the linear solver parameter values.
// Special Notes :
// Scope         : Private
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
void AztecOOSolver::printParams_() const
{
  Xyce::lout() << "\n" << std::endl
               << Xyce::section_divider << std::endl
               << "\n***** Linear solver options:\n" << std::endl
               << "\tPreconditioner:\t\t" << preCond_ << std::endl
               << "\tTolerance:\t\t" << tolerance_ << std::endl
               << "\tMax Iterations:\t\t" << maxIter_ << std::endl
               << Xyce::section_divider << std::endl
               << "\n" << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : AztecOOSolver::doSolve
// Purpose       : Calls the actual solver to solve Ax=b.
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
int AztecOOSolver::doSolve( bool reuse_factors, bool transpose )
{
  int linearStatus = 0;

  double transTime=0.0, precTime=0.0, solveTime=0.0;
  double time1=0.0, time2=0.0;

  // Start the timer...
  timer_->resetStartTime();

  if (VERBOSE_LINEAR)
    time1 = timer_->wallTime();

  if( transform_.get() )
  {
    if( !tProblem_ )
    {
      tProblem_ = &((*transform_)( *problem_ ));
      tProblem_->SetPDL(unsure);
      if( solver_ ) delete solver_;
    }
    std::swap( tProblem_, problem_ );
    transform_->fwd();
  }

  if (VERBOSE_LINEAR)
  {
    time2 = timer_->wallTime();
    transTime += (time2-time1);
  }

  if( !solver_ )
  {
    // Set the traceback mode in Epetra so it prints out warnings
    if (DEBUG_LINEAR)
      dynamic_cast<Epetra_CrsMatrix*>(problem_->GetMatrix())->SetTracebackMode( 2 );

    // Check to see if AztecOO can generate a Krylov subspace large enough for the default settings.
    // AztecOO's Krylov subspace allocation line:  vblock = AZ_manage_memory((kspace+1)*aligned_N_total*sizeof(double), ...
    // Resize KSpace_ if necessary.
    int aligned_dim = problem_->GetRHS()->MyLength() + 2;
    maxKSpace_ = Teuchos::OrdinalTraits<int>::max() / (aligned_dim * sizeof(double)) - 1;
    reduceKSpace_ = true;

    solver_ = new AztecOO( *problem_ );
    setDefaultOptions();
    if( options_ ) setOptions( *options_ );
  }

  // Output the linear system to a Matrix Market file every outputLS_ calls if outputLS_ > 0
  static int file_number = 1, base_file_number = 1;
  if (outputLS_ && !lasProblem_.matrixFree()) {
    if (!(file_number % outputLS_)) {
      char file_name[40];
      if (file_number == 1) {
        EpetraExt::BlockMapToMatrixMarketFile( "Transformed_BlockMap.mm", (problem_->GetMatrix())->Map() );
      }
      sprintf( file_name, "Transformed_Matrix%d.mm", file_number );
      std::string sandiaReq = "Sandia National Laboratories is a multimission laboratory managed and operated by National Technology and\n%";
      sandiaReq += " Engineering Solutions of Sandia LLC, a wholly owned subsidiary of Honeywell International Inc. for the\n%";
      sandiaReq += " U.S. Department of Energy’s National Nuclear Security Administration under contract DE-NA0003525.\n%\n% Xyce circuit matrix.\n%%";
      EpetraExt::RowMatrixToMatrixMarketFile( file_name, *(problem_->GetMatrix()), sandiaReq.c_str() );

      sprintf( file_name, "Transformed_RHS%d.mm", file_number );
      EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(problem_->GetRHS()) );
    }
    // file_number++;  This will be incremented after the solution vector is written to file.
  }
  if (outputBaseLS_ && !lasProblem_.matrixFree()) {
    if (!(base_file_number % outputBaseLS_)) {
      char file_name[40];
      if (base_file_number == 1) {
        EpetraExt::BlockMapToMatrixMarketFile( "Base_BlockMap.mm", (tProblem_->GetMatrix())->Map() );
      }
      sprintf( file_name, "Base_Matrix%d.mm", base_file_number );
      std::string sandiaReq = "Sandia National Laboratories is a multimission laboratory managed and operated by National Technology and\n%";
      sandiaReq += " Engineering Solutions of Sandia LLC, a wholly owned subsidiary of Honeywell International Inc. for the\n%";
      sandiaReq += " U.S. Department of Energy’s National Nuclear Security Administration under contract DE-NA0003525.\n%\n% Xyce circuit matrix.\n%%";
      EpetraExt::RowMatrixToMatrixMarketFile( file_name, *(tProblem_->GetMatrix()), sandiaReq.c_str() );

      sprintf( file_name, "Base_RHS%d.mm", base_file_number );
      EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(tProblem_->GetRHS()) );
    }
    // base_file_number++;  This will be incremented after the solution vector is written to file.
  }

  if (VERBOSE_LINEAR)
    time1 = timer_->wallTime();

  // Change the operator's transpose state if necessary before creating preconditioner.
  bool currTrans = problem_->GetOperator()->UseTranspose();
  if ( transpose != currTrans )
  {
    if (VERBOSE_LINEAR)
      Xyce::dout() << "AztecOO: Operator with transpose state " << currTrans << " is being set to transpose state " << transpose << std::endl;

    // Set the system to solve with the transposed matrix.  Make sure the new preconditioner is set with the solver.
    problem_->GetOperator()->SetUseTranspose( transpose );
    isPrecSet_ = false;
  }

  if( !useAztecPrecond_ )
  {
    Teuchos::RCP<Problem> tmpProblem = Teuchos::rcp( new Problem( Teuchos::rcp(problem_,false) ) );

    // Create the preconditioner if we don't have one.
    if ( Teuchos::is_null( precond_ ) ) {
      TrilinosPrecondFactory factory( *options_ );
      precond_ = factory.create( tmpProblem );
      isPrecSet_ = false;
    }

    // Initialize the values compute the preconditioner.
    bool initRet = precond_->initValues( tmpProblem );
    if (!initRet)
      Report::UserError0() 
        << "Preconditioner::initValues() preconditioner could not be initialized!\n";

    // Compute the preconditioner
    bool compRet = precond_->compute();
    if (!compRet)
      Report::UserError0() 
        <<"Preconditioner::compute() preconditioner could not be computed!";
    if (VERBOSE_LINEAR)
    {
      time2 = timer_->wallTime();
      precTime = time2-time1;
    }

    // Set the preconditioner as an operator if it was just constructed
    if( !isPrecSet_ && precond_->epetraObj()!=Teuchos::null ) {
      solver_->SetPrecOperator( &*precond_->epetraObj() );
      isPrecSet_ = true;
    }

    // If the preconditioner object is null, inform the Aztec solver to use no preconditioning.
    if ( Teuchos::is_null( precond_->epetraObj() ) )
      setAztecOption_( "AZ_precond", 0 );
  }

  linearStatus = solver_->Iterate(maxIter_, tolerance_);

  if (VERBOSE_LINEAR)
  {
    time1 = timer_->wallTime();
    solveTime = time1-time2;
  }

/*
  Epetra_MultiVector* b = problem_->GetRHS();
  int numrhs = b->NumVectors();
  std::vector<double> actual_resids( numrhs ), rhs_norm( numrhs );
  Epetra_MultiVector resid( b->Map(), numrhs );
  problem_->GetOperator()->Apply( *(problem_->GetLHS()), resid );
  resid.Update( -1.0, *b, 1.0 );
  resid.Norm2( &actual_resids[0] );
  b->Norm2( &rhs_norm[0] );
  for (int i=0; i<numrhs; i++ ) {
  Xyce::dout() << "Problem (before transform) " << i << " : \t" << actual_resids[i]/rhs_norm[i] << std::endl;
  }
*/
  if( transform_.get() )
  {
    transform_->rvs();
    std::swap( tProblem_, problem_ );
  }
  
  if (VERBOSE_LINEAR)
  {
    time2 = timer_->wallTime();
    transTime += (time2-time1);
  }

  // Output computed solution vectors, if requested.
  if (outputLS_ && !lasProblem_.matrixFree()) {
    if (!(file_number % outputLS_)) {
      char file_name[40];
      sprintf( file_name, "Transformed_Soln%d.mm", file_number );
      EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(tProblem_->GetLHS()) );
    }
    file_number++;
  }
  if (outputBaseLS_ && !lasProblem_.matrixFree()) {
    if (!(base_file_number % outputBaseLS_)) {
      char file_name[40];
      sprintf( file_name, "Base_Soln%d.mm", base_file_number );
      EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(problem_->GetLHS()) );
    }
    base_file_number++;
  }

/*
  b = problem_->GetRHS();
  Epetra_MultiVector resid2( b->Map(), numrhs );
  problem_->GetOperator()->Apply( *(problem_->GetLHS()), resid2 );
  resid2.Update( -1.0, *b, 1.0 );
  resid2.Norm2( &actual_resids[0] );
  b->Norm2( &rhs_norm[0] );
  for (int i=0; i<numrhs; i++ ) {
  Xyce::dout() << "Problem (after transform)" << i << " : \t" << actual_resids[i]/rhs_norm[i] << std::endl;
  }
*/

  // Get information back from the solver
  linearResidual_ = solver_->TrueResidual();
  numLinearIters_ = solver_->NumIters();

  // Update the total solution time
  solutionTime_ = timer_->elapsedTime();
  if (VERBOSE_LINEAR)
  {
    Xyce::dout() << "Transform Time: " << transTime << std::endl;
    Xyce::dout() << "Preconditioner Assembly Time: " << precTime << std::endl;
    Xyce::dout() << "AztecOO Solve Time: " << solveTime << std::endl;
    Xyce::dout() << "Total Solve Time: " << solutionTime_ << std::endl;
  }
  linearStatus = 0;

  return linearStatus;

}

//-----------------------------------------------------------------------------
// Function      : AztecOOSolver::setAztecCntl_
// Purpose       : Sets aztec controls
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
bool AztecOOSolver::setAztecCntl_( const Util::Param & param )
{
  bool found = false;

  if( param.tag() == "AZ_precond" )         preCond_ = param.getImmutableValue<int>();
  if( param.tag() == "AZ_subdomain_solve" ) subdomainSolve_ = param.getImmutableValue<int>();
  if( param.tag() == "AZ_kspace" )          KSpace_ = param.getImmutableValue<int>();
  if( param.tag() == "AZ_athresh" )         aThresh_ = param.getImmutableValue<double>();
  if( param.tag() == "AZ_rthresh" )         rThresh_ = param.getImmutableValue<double>();
  if( param.tag() == "AZ_ilut_fill" )       ilutFill_ = param.getImmutableValue<double>();
  if( param.tag() == "AZ_drop" )            drop_ = param.getImmutableValue<double>();
  if( param.tag() == "AZ_overlap" )         overlap_ = param.getImmutableValue<int>();
  if( param.tag() == "AZ_output" )          output_ = param.getImmutableValue<int>();
  if( param.tag() == "AZ_diagnostics" )     diagnostics_ = param.getImmutableValue<int>();
  if( param.tag() == "AZ_max_iter" )        maxIter_ = param.getImmutableValue<int>();
  if( param.tag() == "AZ_tol" )             tolerance_ = param.getImmutableValue<double>();

  if( solver_ )
  {
    for( int i = 0; (i < num_AZ_options_) && !found; ++i )
      if( param.tag() == AZ_options_[i] )
      {
        solver_->SetAztecOption( i, param.getImmutableValue<int>() );
        found = true;
        continue;
      }

    if( !found )
      for( int i = 0; (i < num_AZ_params_) && !found; ++i )
        if( param.tag() == AZ_params_[i] )
        {
          solver_->SetAztecParam( i, param.getImmutableValue<double>() );
          found = true;
          continue;
        }
  }

  return found;
}

} // namespace Linear
} // namespace Xyce

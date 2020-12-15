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
// Purpose        : Implementation file for the Belos linear solver interface.
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
#include <sstream>

// ----------   Xyce Includes   ----------
#include <N_UTL_fwd.h>

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_BelosSolver.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Preconditioner.h>
#include <N_LAS_EpetraProblem.h>
#include <N_LAS_EpetraHelpers.h>
#include <N_LAS_TransformTool.h>
#include <N_LAS_TrilinosPrecondFactory.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Timer.h>

#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_InvOperator.h>

#include <BelosBlockGmresSolMgr.hpp>
#include <BelosGCRODRSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosLinearProblem.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Utils.hpp>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : BelosSolver::BelosSolver
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
BelosSolver::BelosSolver(
  Problem &       problem,
  Util::OptionBlock &   options)
  : Solver(true),
  trivialLS_(false),
  outputDefault_(VERBOSE_LINEAR ? Belos::Errors + Belos::Warnings + Belos::StatusTestDetails : 0),
  maxIterDefault_(200),
  KSpaceDefault_(50),
  toleranceDefault_(1.0e-9),
  recycleDefault_(10),
  belosSolverDefault_("Block GMRES"),
  outputLS_(0),
  outputBaseLS_(0),
  lasProblem_(problem),
  updatedParams_(false),
  linearResidual_(1.0),
  tProblem_(0),
  isPrecSet_(false)
{
  EpetraProblem& eprob = dynamic_cast<EpetraProblem&>(lasProblem_);
  problem_ = &(eprob.epetraObj());  

  options_ = Teuchos::rcp( new Util::OptionBlock( options ) );
  timer_ = Teuchos::rcp( new Util::Timer( ) );

  belosParams_ = Teuchos::rcp( new Teuchos::ParameterList() );
  setDefaultOptions();
  setOptions( *options_ );
}

//-----------------------------------------------------------------------------
// Function      : BelosSolver::setDefaultOptions
// Purpose       : resets Aztec options
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 9/25/07
//-----------------------------------------------------------------------------
bool BelosSolver::setDefaultOptions()
{
  output_ = outputDefault_;
  maxIter_ = maxIterDefault_;
  numLinearIters_ = 0;
  KSpace_ = KSpaceDefault_;
  recycle_ = recycleDefault_;
  tolerance_ = toleranceDefault_;
  belosSolver_ = belosSolverDefault_;
  updatedParams_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : BelosSolver::setOptions
// Purpose       : sets Aztec options and params from modelblock
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool BelosSolver::setOptions(const Util::OptionBlock& OB)
{
  belosParams_->set("Verbosity",output_);
  if (VERBOSE_LINEAR)
  {
    belosParams_->set("Output Frequency", 50);
    belosParams_->set("Output Style", (int)(Belos::Brief));
  }
  belosParams_->set("Maximum Iterations",maxIter_);
  belosParams_->set("Num Blocks", KSpace_);
  belosParams_->set("Block Size", 1);
  belosParams_->set("Convergence Tolerance", tolerance_);
  belosParams_->set("Orthogonalization", "ICGS");

  Util::ParamList::const_iterator it_tpL = OB.begin();
  Util::ParamList::const_iterator end_tpL = OB.end();
  for (; it_tpL != end_tpL; ++it_tpL)
  {
    setParam( *it_tpL );
  }

  // Number of vectors in recycle space for GCRODR
  if (belosSolver_ == "GCRODR")
  {
    belosParams_->set("Num Recycled Blocks", recycle_);
  }

  // Request flexible GMRES
  if (belosSolver_ == "FGMRES")
  {
    belosParams_->set("Flexible Gmres", true);
  }

  // store for restart of solver_
  if( &OB != &*options_ )
  {
    options_ = Teuchos::rcp( new Util::OptionBlock(OB) );
  }

  // set singleton filtering as default for iterative solvers
  const Util::Param *param = Util::findParameter(OB.begin(), OB.end(), "TR_singleton_filter");
  if (!param)
  {
    options_->addParam(Util::Param("TR_singleton_filter", 1));
  }

  if (!lasProblem_.matrixFree())
  {
    // create the transformation object if needed.
    if( Teuchos::is_null(transform_) ) transform_ = TransformTool()( OB );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : BelosSolver::setParam
// Purpose       : sets Aztec option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool BelosSolver::setParam( const Util::Param & param )
{
  std::string tag = param.tag();
  std::string uTag = param.uTag();

  // Set our copies of these parameters.
  if( tag == "AZ_max_iter" )
    setMaxIter(param.getImmutableValue<int>());
  if( tag == "AZ_kspace" )
    setKSpace(param.getImmutableValue<int>());
  else if( tag == "AZ_tol" )
    setTolerance(param.getImmutableValue<double>());
  else if( uTag == "OUTPUT_LS" )
    outputLS_ = param.getImmutableValue<int>();
  else if( uTag == "OUTPUT_BASE_LS" )
    outputBaseLS_ = param.getImmutableValue<int>();
  else if( tag == "BELOS_SOLVER_TYPE" )
   belosSolver_ = param.usVal();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : BelosSolver::getInfo
// Purpose       : sets Aztec option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool BelosSolver::getInfo( Util::Param & param )
{
  if( param.tag() == "AZ_max_iter" )
    param.setVal( maxIter_ );
  else if( param.tag() == "Iterations" )
    param.setVal( (int)numLinearIters_ );
  else if( param.tag() == "AZ_tol" )
    param.setVal( tolerance_ );
  else
    return false;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : BelosSolver::setBelosOption
// Purpose       : sets Belos option
// Special Notes : Takes a string as the option identifier
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Computational Sciences
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool BelosSolver::setBelosOption_(const char * paramName,
                                            const int val)
{
  return setBelosCntl_( Util::Param(paramName, val) );
}

//-----------------------------------------------------------------------------
// Function      : BelosSolver::setBelosParam
// Purpose       : sets Belos parameter
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool BelosSolver::setBelosParam_(const char * paramName,
                                         const double val)
{
  return setBelosCntl_( Util::Param(paramName, val) );
}

//-----------------------------------------------------------------------------
// Function      : BelosSolver::doSolve
// Purpose       : Calls the actual solver to solve Ax=b.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
int BelosSolver::doSolve( bool reuse_factors, bool transpose )
{
  double transTime=0.0, precTime=0.0, solveTime=0.0;
  double time1=0.0, time2=0.0;

  // Start the timer...
  timer_->resetStartTime();

  if (VERBOSE_LINEAR)
    time1 = timer_->wallTime();

  // The Epetra_LinearProblem, prob, is the linear system being solved.
  // It will point to either the original linear system or transformed system.
  Epetra_LinearProblem * prob = problem_;

  if( !Teuchos::is_null(transform_) )
  {
    if( !tProblem_ )
    {
      tProblem_ = &((*transform_)( *problem_ ));
      solver_ = Teuchos::null;
    }
    prob = tProblem_;
    transform_->fwd();
  }

  if (VERBOSE_LINEAR)
  {
    time2 = timer_->wallTime();
    transTime += (time2-time1);
  }

  if( Teuchos::is_null(solver_) && !trivialLS_ ) {

    Teuchos::RCP< Epetra_Operator > A;
    if (lasProblem_.matrixFree()) {
      A = Teuchos::rcp( prob->GetOperator(), false );
    }
    else {
      A = Teuchos::rcp( prob->GetMatrix(), false );
    }
    Teuchos::RCP< Epetra_MultiVector> X = Teuchos::rcp( prob->GetLHS(), false );
    Teuchos::RCP< Epetra_MultiVector> B = Teuchos::rcp( prob->GetRHS(), false );
    belosProblem_ =
      Teuchos::rcp( new Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator>( A, X, B ) );

    if ( X->GlobalLength() != 0 )
    {
      // Create Belos solver (Block GMRES is default)
      if (belosSolver_ == "GCRODR") {
        solver_ = Teuchos::rcp( new Belos::GCRODRSolMgr<double,Epetra_MultiVector,Epetra_Operator>() );
      }
      else {
        solver_ = Teuchos::rcp( new Belos::BlockGmresSolMgr<double,Epetra_MultiVector,Epetra_Operator>() );
      }

      // Reduce the size of the Krylov space if the linear system is smaller than the default size
      // to avoid warning messages.
      int lsDim = X->GlobalLength();
      if (lsDim < KSpace_) {
        KSpace_ = lsDim;
        belosParams_->set("Num Blocks", KSpace_);
      }
    }
    else
    { 
      trivialLS_ = true;
    }
  }

  // If it is a trivial system, the transform has reduced the linear system to zero unknowns.
  // Just perform the reverse transform and return.
  if( trivialLS_ && !Teuchos::is_null(transform_) )
  {
    transform_->rvs();
    return 0;
  }

  // Output the linear system to a Matrix Market file every outputLS_ calls if outputLS_ > 0
  static int file_number = 1, base_file_number = 1;
  if (outputLS_ && !lasProblem_.matrixFree()) {
    if (!(file_number % outputLS_)) {
      if (!reuse_factors) {
        Xyce::Linear::writeToFile( *prob, "Transformed", file_number, (file_number == 1) );
      }
    }
    // file_number++;  This is incremented below after the solution is written to file.
  }
  if (outputBaseLS_ && !lasProblem_.matrixFree()) {
    if (!(base_file_number % outputBaseLS_)) {
      if (!reuse_factors) {
        Xyce::Linear::writeToFile( *problem_, "Base", base_file_number, (base_file_number == 1) );
      }
    }
    // base_file_number++;  This is incremented below after the solution is written to file.
  }

  if (VERBOSE_LINEAR)
    time1 = timer_->wallTime();

  Teuchos::RCP<Problem> tmpProblem = Teuchos::rcp( new EpetraProblem( Teuchos::rcp(prob,false) ) );

  // Create the preconditioner if we don't have one.
  if ( Teuchos::is_null( precond_ ) ) {
    TrilinosPrecondFactory factory( *options_ );
    precond_ = factory.create( tmpProblem );
    isPrecSet_ = false;
  }

  // Change the operator's transpose state if necessary before creating preconditioner.
  bool currTrans = prob->GetOperator()->UseTranspose();
  if ( transpose != currTrans )
  {
    if (VERBOSE_LINEAR)
      Xyce::dout() << "Belos: Operator with transpose state " << currTrans << " is being set to transpose state " << transpose << std::endl;

    // Set the system to solve with the transposed matrix.  Make sure the new preconditioner is set with the solver.
    prob->GetOperator()->SetUseTranspose( transpose );
  }

  // Initialize the values compute the preconditioner.
  bool initRet = precond_->initValues( tmpProblem );
  if (!initRet)
    Xyce::Report::DevelFatal0().in("Preconditioner::initValues()") << "preconditioner could not be initialized!";

  // Compute the preconditioner
  bool compRet = precond_->compute();
  if (!compRet)
    Xyce::Report::DevelFatal0().in("Preconditioner::compute()") << "preconditioner could not be computed!";

  // Set the preconditioner as an operator if it was just constructed, else inform the Belos
  // solver to use no preconditioning.
  if ( !isPrecSet_ ) {
    if ( !Teuchos::is_null( precond_->epetraObj() ) ) {
      belosPrecond_ = Teuchos::rcp( new Epetra_InvOperator( &*precond_->epetraObj() ) );
    }
      else {
      belosPrecond_ = Teuchos::null;
    }
    belosProblem_->setRightPrec( belosPrecond_ );
    try {
      belosProblem_->setProblem();
    }
    catch ( std::exception& e )
    {
      Xyce::dout() << "BELOS ERROR: " << e.what() << std::endl;
      return 1;
    }
    solver_->setProblem( belosProblem_ );
    isPrecSet_ = true;
  }

  // Make sure the operators are wrapped if we are solving the adjoint system.
  if ( transpose != currTrans )
  {
    if ( transpose )
    {
      // Wrap the operator.
      Teuchos::RCP<EpetraTransOp> transA = Teuchos::rcp( new EpetraTransOp( Teuchos::rcp( prob->GetOperator(), false ) ) );
      belosProblem_->setOperator( transA );   
    
      // Wrap the preconditioner, if we have one.
      if (belosPrecond_ != Teuchos::null)
      {
        Teuchos::RCP<EpetraTransOp> transP = Teuchos::rcp( new EpetraTransOp( belosPrecond_ ) );
        belosProblem_->setRightPrec( transP );
      }
    }
    else
    {
      // Change the linear problem from the transposed problem, back to the original forward problem.
      belosProblem_->setOperator( Teuchos::rcp( prob->GetOperator(), false ) );
      belosProblem_->setRightPrec( belosPrecond_ );   
    } 
  }

  if (VERBOSE_LINEAR)
  { 
    time2 = timer_->wallTime();
    precTime = time2-time1;
  }

  // Reset the problem for the next solve.
  try {
    belosProblem_->setProblem();
  }
  catch ( std::exception& e )
  {
    Xyce::dout() << "BELOS ERROR: " << e.what() << std::endl;
    return 1;
  }
 
  solver_->setParameters( belosParams_ );

  // Solve the problem.
  Belos::ReturnType linearStatus;
  try {
    linearStatus = solver_->solve();
  }
  catch ( std::exception& e )
  {
    linearStatus = Belos::Unconverged;
    Xyce::dout() << "BELOS ERROR: " << e.what() << std::endl;
    return 1;
  }
  numLinearIters_ = solver_->getNumIters();

  if (VERBOSE_LINEAR)
  {
    time1 = timer_->wallTime();
    solveTime = time1 - time2;
  }
/*
  Epetra_MultiVector* b = prob->GetRHS();
  int numrhs = b->NumVectors();
  std::vector<double> actual_resids( numrhs ), rhs_norm( numrhs );
  Epetra_MultiVector resid( b->Map(), numrhs );
  prob->GetOperator()->Apply( *(prob->GetLHS()), resid );
  resid.Update( -1.0, *b, 1.0 );
  resid.Norm2( &actual_resids[0] );
  b->Norm2( &rhs_norm[0] );
  for (int i=0; i<numrhs; i++ ) {
    Xyce::dout() << "Problem " << i << " : \t" << actual_resids[i]/rhs_norm[i] << std::endl;
  }
*/

 if( !Teuchos::is_null(transform_) ) transform_->rvs();

 if (VERBOSE_LINEAR)
 {
   time2 = timer_->wallTime();
   transTime += (time2-time1);
 }

/*
  Epetra_MultiVector* b2 = problem_->GetRHS();
  Epetra_MultiVector resid2( b2->Map(), numrhs );
  problem_->GetOperator()->Apply( *(problem_->GetLHS()), resid2 );
  resid2.Update( -1.0, *b2, 1.0 );
  resid2.Norm2( &actual_resids[0] );
  b2->Norm2( &rhs_norm[0] );
  for (int i=0; i<numrhs; i++ ) {
    Xyce::dout() << "Problem " << i << " : \t" << actual_resids[i]/rhs_norm[i] << std::endl;
  }
*/

  // Output computed solution vectors, if requested.
  if (outputLS_ && !lasProblem_.matrixFree()) {
    if (!(file_number % outputLS_)) {
      Teuchos::RCP<Problem> las_prob = Teuchos::rcp( new EpetraProblem( Teuchos::rcp( prob, false ) ) );
      char file_name[40];
      sprintf( file_name, "Transformed_Soln%d.mm", file_number );
      las_prob->getLHS()->writeToFile( file_name, false, true );
    }
    file_number++;
  }
  if (outputBaseLS_ && !lasProblem_.matrixFree()) {
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
  {
    Xyce::dout() << "Transform Time: " << transTime << std::endl;
    Xyce::dout() << "Preconditioner Assembly Time: " << precTime << std::endl;
    Xyce::dout() << "Belos Solve Time: " << solveTime << std::endl;
    Xyce::dout() << "Total Solve Time: " << solutionTime_ << std::endl;
  }

  // Belos does not return the residual, so use a fake residual that indicates convergence
  // if the solver converged.
  if (linearStatus == Belos::Unconverged)
  {
    if (VERBOSE_LINEAR && solver_->isLOADetected())
    {
      Report::UserWarning0() 
        << "Belos::BlockGmresSolMgr has detected a loss of accuracy!";
    }
  }

  linearResidual_ = tolerance_/10;
  return 0;
}

} // namespace Linear
} // namespace Xyce

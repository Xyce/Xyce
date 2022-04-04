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
// Filename       : $RCSfile: N_LAS_Amesos2Solver.C,v $
//
// Purpose        : Amesos2 direct solver wrapper
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/20/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.3 $
//
// Revision Date  : $Date: 2016/03/16 22:14:27 $
//
// Current Owner  : $Author: hkthorn $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------

#include <Epetra_LinearProblem.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Export.h>
#include <Epetra_Map.h>

// ---------- Xyce Includes ----------

#include <N_UTL_fwd.h>

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Amesos2Solver.h>
#include <N_LAS_Problem.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_EpetraProblem.h>
#include <N_LAS_EpetraHelpers.h>
#include <N_LAS_TransformTool.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Timer.h>

#include <Teuchos_Utils.hpp>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : Amesos2Solver::Amesos2Solver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
Amesos2Solver::Amesos2Solver(
  const std::string &   type,
  Problem &       problem,
  Util::OptionBlock &   options)
  : Solver(problem, false),
    type_(type),
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
// Function      : Amesos2Solver::~Amesos2Solver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
Amesos2Solver::~Amesos2Solver()
{
  delete timer_;
  delete options_;
}

//-----------------------------------------------------------------------------
// Function      : Amesos2Solver::setOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool Amesos2Solver::setOptions( const Util::OptionBlock & OB )
{
  bool foundAMD = false, foundPartition = false, foundSingleton = false;
  
  for( Util::ParamList::const_iterator it_tpL = OB.begin();
         it_tpL != OB.end(); ++it_tpL )
  {
    std::string tag = it_tpL->uTag();

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
// Function      : Amesos2Solver::doSolve
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
int Amesos2Solver::doSolve( bool reuse_factors, bool transpose )
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

  if( Teuchos::is_null( solver_ ) )
  {
    Teuchos::ParameterList amesos2_params("Amesos2");

    // Set some parameters on the solver, which are not used at this time.
    if (type_ == "SHYLU_BASKER") {
      //amesos2_params.sublist("ShyLUBasker").set("pivot", true);
      //amesos2_params.sublist("ShyLUBasker").set("realloc", true);
      //amesos2_params.sublist("ShyLUBasker").set("verbose", false );
    } 

#ifdef Xyce_AMESOS2_SHYLUBASKER
    if (type_ == "SHYLU_BASKER") {

      solver_ = Amesos2::create<Epetra_CrsMatrix,Epetra_MultiVector>("ShyLUBasker", 
                                                                     Teuchos::rcp(dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix()),false));

    }
#endif
#ifdef Xyce_AMESOS2_BASKER
    else if (type_ == "BASKER") {

      solver_ = Amesos2::create<Epetra_CrsMatrix,Epetra_MultiVector>("Basker", 
                                                                     Teuchos::rcp(dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix()),false));

    }
#endif
#ifdef Xyce_AMESOS2_KLU2
    else if (type_ == "KLU2") {
      
      solver_ = Amesos2::create<Epetra_CrsMatrix,Epetra_MultiVector>("klu2", 
                                                                     Teuchos::rcp(dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix()),false));
    }
#endif
    else {
      Report::DevelFatal0()
        << "Unknown or Unavailable Linear Solver: " << type_;
 
    }

    solver_->setParameters( rcpFromRef(amesos2_params) );

    double begSymTime = timer_->elapsedTime();

    try 
    {
      // Perform symbolic factorization and check return value for failure
      solver_->symbolicFactorization();
    }
    catch (std::runtime_error& e)
    {
      if (VERBOSE_LINEAR)
      {
        Xyce::dout() << "  Amesos2 (" << type_ << ") error: " << e.what() << std::endl;
      }
      linearStatus = -1;
      solver_ = Teuchos::null;

      // Inform user that singular matrix was found and linear solve has failed.
      Report::UserWarning0() 
        << "Symbolically singular matrix found by Amesos2, returning zero solution to nonlinear solver!";

      // Put zeros in the solution since Amesos2 was not able to solve this problem
      prob->GetLHS()->PutScalar( 0.0 );
      // Output the singular linear system to a Matrix Market file if outputFailedLS_ > 0
      if (outputFailedLS_) {
        failure_number++;
        Xyce::Linear::writeToFile( *prob, "Failed", failure_number, (failure_number == 1) );
      }

      return linearStatus;  // return the actual status (see bug 414 SON)
    }

    if (VERBOSE_LINEAR)
    {
      double endSymTime = timer_->elapsedTime();
      Xyce::dout() << "  Amesos2 (" << type_ << ") Symbolic Factorization Time: "
                   << (endSymTime - begSymTime) << std::endl;
    }
  }

  // Set transpose flag
  if (transpose) {
    Teuchos::ParameterList amesos2_params("Amesos2");
    amesos2_params.set("Transpose", true);
    if (type_ == "SHYLU_BASKER") {
      amesos2_params.sublist("ShyLUBasker").set("transpose", transpose);
    }
    else if (type_ == "KLU2") {
      amesos2_params.sublist("KLU2").set("Trans","TRANS","Solve with transpose");
    }
    solver_->setParameters( rcpFromRef(amesos2_params) );
  }


  if( !reuse_factors ) 
  {
    double begNumTime = timer_->elapsedTime();

    try 
    {
      solver_->numericFactorization();
      if (VERBOSE_LINEAR)
      {
        double endNumTime = timer_->elapsedTime();
        Xyce::dout() << "  Amesos2 (" << type_ << ") Numeric Factorization Time: "
                     << (endNumTime - begNumTime) << std::endl;
      }  
    }
    catch (std::runtime_error& e)
    {
      if (VERBOSE_LINEAR)
      {
        Xyce::dout() << "  Amesos2 (" << type_ << ") error: " << e.what() << std::endl;
      }
      linearStatus = -1;
    }
 
    if (linearStatus != 0) {

      // Inform user that singular matrix was found and linear solve has failed.
      Report::UserWarning0() 
        << "Numerically singular matrix found by Amesos2, returning zero solution to nonlinear solver!";

      // Put zeros in the solution since Amesos2 was not able to solve this problem
      prob->GetLHS()->PutScalar( 0.0 );
      // Output the singular linear system to a Matrix Market file if outputFailedLS_ > 0
      if (outputFailedLS_) {
        failure_number++;
        Xyce::Linear::writeToFile( *prob, "Failed", failure_number, (failure_number == 1) );
      }

      return linearStatus;  // return the actual status (see bug 414 SON)
    }
  }

  // Perform linear solve using factorization
  double begSolveTime = timer_->elapsedTime();

  solver_->solve( prob->GetLHS(), prob->GetRHS() );

  double endSolveTime = timer_->elapsedTime();

  if (VERBOSE_LINEAR)
  {
    Xyce::dout() << "  Amesos2 (" << type_ << ") Solve Time: "
                 << (endSolveTime - begSolveTime) << std::endl;
  }

  if (DEBUG_LINEAR)
  {
    int numrhs = prob->GetLHS()->NumVectors();
    std::vector<double> resNorm(numrhs,0.0), xNorm(numrhs,0.0), bNorm(numrhs,0.0);
    Epetra_MultiVector res( prob->GetLHS()->Map(), prob->GetLHS()->NumVectors() );
    bool oldTrans = prob->GetOperator()->UseTranspose();
    prob->GetOperator()->SetUseTranspose( transpose );
    prob->GetOperator()->Apply( *(prob->GetLHS()), res );
    prob->GetOperator()->SetUseTranspose( false );
    res.Update( 1.0, *(prob->GetRHS()), -1.0 );
    res.Norm2( &resNorm[0] );
    prob->GetRHS()->Norm2( &bNorm[0] );
    prob->GetLHS()->Norm2( &xNorm[0] );
    Xyce::lout() << "Linear System Residual (AMESOS_" << type_ << "): " << std::endl;
    for (int i=0; i<numrhs; i++)
    {
      if (bNorm[i] != 0.0)
        std::cout << "  Problem " << i << " : " << (resNorm[i]/bNorm[i]) << ", soln : " << xNorm[i] << std::endl;
      else
        std::cout << "  Problem " << i << " : " << resNorm[i] << ", soln : " << xNorm[i] << std::endl;
    }
  }

  // Unset transpose flag
  if (transpose) {
    Teuchos::ParameterList amesos2_params("Amesos2");
    amesos2_params.set("Transpose", false);
    if (type_ == "SHYLU_BASKER") {
      amesos2_params.sublist("ShyLUBasker").set("transpose", false);
    }
    else if (type_ == "KLU2") {
      amesos2_params.sublist("KLU2").set("Trans","NOTRANS","Solve with transpose");
    }
    solver_->setParameters( rcpFromRef(amesos2_params) );
  }

  if( transform_.get() ) transform_->rvs();

  // Output computed solution vectors, if requested.
  if (outputLS_) {
    if (!(file_number % outputLS_)) {
      Teuchos::RCP<Problem> las_prob = Teuchos::rcp( new EpetraProblem( Teuchos::rcp( problem_, false ) ) );
      char file_name[40];
      sprintf( file_name, "Transformed_Soln%d.mm", file_number );
      las_prob->getLHS()->writeToFile( file_name, false, true );
    }
  }
  file_number++;
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
    Xyce::dout() << "Total Linear Solution Time (Amesos2 " << type_ << "): "
                 << solutionTime_ << std::endl;

  return 0;
}

} // namespace Linear
} // namespace Xyce

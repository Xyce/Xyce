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
#include <N_LAS_TransformTool.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Timer.h>

#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_BlockMapOut.h>

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
  : Solver(false),
    type_(type),
    lasProblem_(problem),
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

  if( !transform_.get() ) transform_ = TransformTool()( *options_ );

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

  Epetra_LinearProblem * prob = problem_;

  if( transform_.get() )
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
      char file_name[40];
      if (!reuse_factors) {
        if (file_number == 1) {
          EpetraExt::BlockMapToMatrixMarketFile( "Transformed_BlockMap.mm", (prob->GetMatrix())->Map() );
        }
        sprintf( file_name, "Transformed_Matrix%d.mm", file_number );

        std::string sandiaReq = "Sandia National Laboratories is a multimission laboratory managed and operated by National Technology and\n%";
        sandiaReq += " Engineering Solutions of Sandia LLC, a wholly owned subsidiary of Honeywell International Inc. for the\n%";
        sandiaReq += " U.S. Department of Energy’s National Nuclear Security Administration under contract DE-NA0003525.\n%\n% Xyce circuit matrix.\n%%";

        EpetraExt::RowMatrixToMatrixMarketFile( file_name, *(prob->GetMatrix()), sandiaReq.c_str() );
        sprintf( file_name, "Transformed_RHS%d.mm", file_number );
        EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(prob->GetRHS()) );
      }
    }
    // file_number++;  This will be incremented after the solution vector is written to file.
  }
  if (outputBaseLS_) {
    if (!(base_file_number % outputBaseLS_)) {
      char file_name[40];
      if (!reuse_factors) {
        if (base_file_number == 1) {
          EpetraExt::BlockMapToMatrixMarketFile( "Base_BlockMap.mm", (problem_->GetMatrix())->Map() );
        }
        sprintf( file_name, "Base_Matrix%d.mm", base_file_number );
        lasProblem_.getMatrix()->writeToFile( file_name, false, true );
        sprintf( file_name, "Base_RHS%d.mm", base_file_number );
        lasProblem_.getRHS()->writeToFile( file_name, false, true );
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

    if (type_ == "BASKER") {
      //amesos2_params.sublist("Basker").set("num_threads", 12);
      amesos2_params.sublist("Basker").set("matching", true);
      amesos2_params.sublist("Basker").set("matching_type", 0);
      amesos2_params.sublist("Basker").set("btf",true);
      amesos2_params.sublist("Basker").set("amd_btf", true);
      amesos2_params.sublist("Basker").set("amd_dom", true);
      amesos2_params.sublist("Basker").set("transpose", false);
      amesos2_params.sublist("Basker").set("symmetric", false);
      amesos2_params.sublist("Basker").set("pivot", false);
      amesos2_params.sublist("Basker").set("pivot_tol", .001);
      amesos2_params.sublist("Basker").set("realloc", false);
      //amesos2_params.sublist("Basker").set("verbose", (bool)VERBOSE_LINEAR);

      solver_ = Amesos2::create<Epetra_CrsMatrix,Epetra_MultiVector>("Basker", 
                                                                     Teuchos::rcp(dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix()),false));

      solver_->setParameters( rcpFromRef(amesos2_params) );

    }
    else if (type_ == "KLU2") {
      
      solver_ = Amesos2::create<Epetra_CrsMatrix,Epetra_MultiVector>("klu2", 
                                                                     Teuchos::rcp(dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix()),false));
    }
      
    double begSymTime = timer_->elapsedTime();

    // Perform symbolic factorization and check return value for failure
    solver_->symbolicFactorization();

    if (VERBOSE_LINEAR)
    {
      double endSymTime = timer_->elapsedTime();
      Xyce::dout() << "  Amesos2 (" << type_ << ") Symbolic Factorization Time: "
                   << (endSymTime - begSymTime) << std::endl;
    }
  }

  if( !reuse_factors ) {

    double begNumTime = timer_->elapsedTime();

    solver_->numericFactorization();
    if (VERBOSE_LINEAR)
    {
      double endNumTime = timer_->elapsedTime();
      Xyce::dout() << "  Amesos2 (" << type_ << ") Numeric Factorization Time: "
                   << (endNumTime - begNumTime) << std::endl;
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

  if (DEBUG_LINEAR) {
    double resNorm = 0.0, bNorm = 0.0;
    Epetra_MultiVector res( prob->GetLHS()->Map(), prob->GetLHS()->NumVectors() );
    prob->GetOperator()->Apply( *(prob->GetLHS()), res );
    res.Update( 1.0, *(prob->GetRHS()), -1.0 );
    res.Norm2( &resNorm );
    prob->GetRHS()->Norm2( &bNorm );
    Xyce::lout() << "Linear System Residual (AMESOS2_" << type_ << "): "
                 << (resNorm/bNorm) << std::endl;
  }

  if( transform_.get() ) transform_->rvs();

  // Output computed solution vectors, if requested.
  if (outputLS_) {
    if (!(file_number % outputLS_)) {
      char file_name[40];
      sprintf( file_name, "Transformed_Soln%d.mm", file_number );
      EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(problem_->GetLHS()) );
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

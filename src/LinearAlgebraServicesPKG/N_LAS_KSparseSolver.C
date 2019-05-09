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
// Purpose        : KSparse direct solver wrapper
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/20/04
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
#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

// ---------- Xyce Includes ----------

#include <N_LAS_KSparseSolver.h>

#include <N_LAS_Problem.h>

#include <N_LAS_TransformTool.h>

#include <N_UTL_Timer.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_FeatureTest.h>

#include <N_ERH_ErrorMgr.h>

#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_BlockMapOut.h>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Utils.hpp>

#include <Epetra_CrsKundertSparse.h>
#include <Amesos.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : KSparseSolver::KSparseSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
KSparseSolver::KSparseSolver(Problem & prob, Util::OptionBlock & options)
 : Solver(false),
   lasProblem_(prob),
   problem_(prob.epetraObj()),
   outputLS_(0),
   outputBaseLS_(0),
   outputFailedLS_(0),
   tProblem_(0),
   options_( new Util::OptionBlock( options ) )
{
  timer_ = Teuchos::rcp( new Util::Timer( ) );
  setOptions( options );
}

//-----------------------------------------------------------------------------
// Function      : KSparseSolver::~KSparseSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
KSparseSolver::~KSparseSolver()
{
  if( options_ )   delete options_;
}

//-----------------------------------------------------------------------------
// Function      : KSparseSolver::setOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool KSparseSolver::setOptions( const Util::OptionBlock & OB )
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

#ifdef Xyce_PARALLEL_MPI
  options_->addParam( Util::Param( "TR_reindex", 1 ) );

  // Turn off partitioning and AMD if we're doing a parallel load serial solve
  options_->addParam( Util::Param( "TR_partition", 0 ) );
  options_->addParam( Util::Param( "TR_amd", 0 ) );
#endif

  if( !transform_.get() ) transform_ = TransformTool()( *options_ );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : KSparseSolver::setDefaultOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool KSparseSolver::setDefaultOptions()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : KSparseSolver::setParam
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool KSparseSolver::setParam( const Util::Param & param )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : KSparseSolver::getInfo
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool KSparseSolver::getInfo( Util::Param & info )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : KSparseSolver::doSolve
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
int KSparseSolver::doSolve( bool reuse_factors, bool transpose )
{
  if ( transpose ) {
    // Inform user that a nontrivial matrix was found and linear solve has failed.
    Report::UserError0()
      <<"KSparse linear solver does not support transpose solves at this time.";
  }

  // Start the timer...
  timer_->resetStartTime();

  int linearStatus = 0;

  Epetra_LinearProblem * prob = &problem_;

  if( transform_.get() )
  {
    if( !tProblem_ )
      tProblem_ = &((*transform_)( problem_ ));
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
    file_number++;
  }
  if (outputBaseLS_) {
    if (!(base_file_number % outputBaseLS_)) {
      char file_name[40];
      if (!reuse_factors) {
        if (base_file_number == 1) {      
          EpetraExt::BlockMapToMatrixMarketFile( "Base_BlockMap.mm", (problem_.GetMatrix())->Map() );
        }
        sprintf( file_name, "Base_Matrix%d.mm", base_file_number );
        std::string sandiaReq = "Sandia National Laboratories is a multimission laboratory managed and operated by National Technology and\n%";
        sandiaReq += " Engineering Solutions of Sandia LLC, a wholly owned subsidiary of Honeywell International Inc. for the\n%";
        sandiaReq += " U.S. Department of Energy’s National Nuclear Security Administration under contract DE-NA0003525.\n%\n% Xyce circuit matrix.\n%%";
        EpetraExt::RowMatrixToMatrixMarketFile( file_name, *(problem_.GetMatrix()), sandiaReq.c_str() );
        sprintf( file_name, "Base_RHS%d.mm", base_file_number );
        EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(problem_.GetRHS()) );
      }
    } 
    base_file_number++;
  }
 
  // Set the traceback mode in Epetra so it prints out warnings
  if (DEBUG_LINEAR)
    dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix())->SetTracebackMode( 2 );
  else
    dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix())->SetTracebackMode( 0 );

  // Import the parallel matrix to a serial one, if necessary.
  prob = importToSerial();

  // Create solver if one doesn't exist.
  if (solver_ == Teuchos::null)
    solver_ = Teuchos::rcp( new Epetra_CrsKundertSparse( prob ) );

  // Perform linear solve using factorization
  double begSolveTime = timer_->elapsedTime();

  linearStatus = solver_->Solve( !reuse_factors );

  int ksparseStatus = linearStatus;

  // Ksparse is not as robust as KLU at this time, so use KLU if a numerical failure occurs.
  if (linearStatus != 0)
  {
    // Create a KLU solver if we don't have one.
    if (kluSolver_ == Teuchos::null)
    {
      Amesos amesosFactory;
      kluSolver_ = Teuchos::rcp( amesosFactory.Create( "Amesos_Klu", *prob ) );

      // Perform symbolic factorization.
      linearStatus = kluSolver_->SymbolicFactorization();
    }

    // Perform numeric factorization with KLU
    linearStatus = kluSolver_->NumericFactorization();

    // Perform solve with KLU
    if (linearStatus == 0)
      kluSolver_->Solve();
  }
  

  // Export solution back to global system, if necessary.
  exportToGlobal();

  if (VERBOSE_LINEAR)
  {
    double endSolveTime = timer_->elapsedTime();
    Xyce::lout() << "  KSparse Solve Time: " << (endSolveTime - begSolveTime) << std::endl;
  }

  if (linearStatus != 0) {

    // Inform user that singular matrix was found and linear solve has failed.
    Report::UserWarning0() 
      << "Numerically singular matrix found by KSparse (err = " << ksparseStatus << "), returning zero solution to nonlinear solver!";

    // Put zeros in the solution since KSparse was not able to solve this problem
    prob->GetLHS()->PutScalar( 0.0 );
    // Output the singular linear system to a Matrix Market file if outputFailedLS_ > 0
    if (outputFailedLS_) {
      failure_number++;
      char file_name[40];
      if (failure_number == 1) {
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
  }

  if (DEBUG_LINEAR && (linearStatus==0))
  {
    double resNorm = 0.0, bNorm = 0.0;
    Epetra_MultiVector res( prob->GetLHS()->Map(), prob->GetLHS()->NumVectors() );
    prob->GetOperator()->Apply( *(prob->GetLHS()), res );
    res.Update( 1.0, *(prob->GetRHS()), -1.0 );
    res.Norm2( &resNorm );
    prob->GetRHS()->Norm2( &bNorm );
    Xyce::lout() << "Linear System Residual (KSparse) : " << (resNorm/bNorm) << std::endl;
  }

  if( transform_.get() ) transform_->rvs();
  
  // Update the total solution time
  solutionTime_ = timer_->elapsedTime();

  if (VERBOSE_LINEAR)
    Xyce::lout() << "Total Linear Solution Time (KSparse): " << solutionTime_ << std::endl;

  return linearStatus; // see bug 414 SON
}

//-----------------------------------------------------------------------------
// Function      : KSparseSolver::importToSerial
// Purpose       : Import a distributed matrix to a serial one for the solver
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystems Modeling
// Creation Date : 07/01/10
//-----------------------------------------------------------------------------
Epetra_LinearProblem * KSparseSolver::importToSerial()
{
#ifdef Xyce_PARALLEL_MPI
    Epetra_CrsMatrix * origMat = dynamic_cast<Epetra_CrsMatrix *>(problem_.GetOperator());

    // If we haven't set up the maps for serial matrix storage on proc 0, then do it now.
    if (serialMap_ == Teuchos::null) {
      const Epetra_Map& origMap = origMat->RowMap();
      int MyPID = origMap.Comm().MyPID(); 
      int NumGlobalElements = origMap.NumGlobalElements();
      int NumMyElements = NumGlobalElements;
      if (MyPID != 0)
        NumMyElements = 0;
      serialMap_ = Teuchos::rcp(new Epetra_Map(-1, NumMyElements, 0, origMap.Comm()));
      int NumVectors = problem_.GetRHS()->NumVectors() ;
    
      serialLHS_ = Teuchos::rcp( new Epetra_MultiVector( *serialMap_, NumVectors ));
      serialRHS_ = Teuchos::rcp (new Epetra_MultiVector( *serialMap_, NumVectors ));
      serialImporter_ = Teuchos::rcp(new Epetra_Import( *serialMap_, origMap ));
      serialMat_ = Teuchos::rcp( new Epetra_CrsMatrix( Copy, *serialMap_, 0 )) ;
      serialProblem_ = Teuchos::rcp( new Epetra_LinearProblem( &*serialMat_, &*serialLHS_, &*serialRHS_ ) );
    }

    // Import linear system from problem
    serialMat_->Import( *origMat, *serialImporter_, Insert );
    serialMat_->FillComplete();
    serialMat_->OptimizeStorage(); 
    serialLHS_->Import( *problem_.GetLHS(), *serialImporter_, Insert );
    serialRHS_->Import( *problem_.GetRHS(), *serialImporter_, Insert );   

    return &*serialProblem_;
#else
    return &problem_;
#endif
} 

//-----------------------------------------------------------------------------
// Function      : KSparseSolver::exportToGlobal
// Purpose       : Export the serial solution to the global one
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystems Modeling
// Creation Date : 07/01/10
//-----------------------------------------------------------------------------
int KSparseSolver::exportToGlobal()
{
#ifdef Xyce_PARALLEL_MPI
  // Return solution back to global problem
  problem_.GetLHS()->Export( *serialLHS_, *serialImporter_, Insert ) ;
#endif
  return 0;
}

} // namespace Linear
} // namespace Xyce

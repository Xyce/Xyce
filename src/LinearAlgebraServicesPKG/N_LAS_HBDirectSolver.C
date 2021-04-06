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
// Purpose        :  direct solver wrapper
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

// ---------- Xyce Includes ----------

#include <N_UTL_fwd.h>
#include <N_UTL_Math.h>

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_HBDirectSolver.h>
#include <N_LAS_HBBuilder.h>
#include <N_LOA_HBLoader.h>
#include <N_LAS_Builder.h>
#include <N_LAS_Solver.h>
#include <N_LAS_Problem.h>
#include <N_LAS_FilteredMatrix.h>
#include <N_LAS_TransformTool.h>
#include <N_LAS_EpetraVector.h>
#include <N_LAS_EpetraHelpers.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Timer.h>
#include <N_UTL_AssemblyTypes.h>
#include <N_PDS_Comm.h>
#include <N_PDS_EpetraParMap.h>

#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Import.h>
#include <Epetra_Util.h>
#include <Epetra_Comm.h>

#include <Teuchos_Utils.hpp>
#include <Teuchos_BLAS.hpp>

#include <set>
#include <utility>
#include <numeric>
#include <iostream>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : HBDirectSolver::HBDirectSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
HBDirectSolver::HBDirectSolver(
  Builder &       builder,
  Problem &       problem,
  Util::OptionBlock &   options)
  : Solver(false),
    builder_(builder),
    lasProblem_(problem),
    isInit_(false),
    hbOsc_(false),
    N_(0),
    n_(0),
    M_(0),
    numAugRows_(0),
    outputLS_(0),
    solver_(""),
    solverDefault_("LAPACK"),
    options_( new Util::OptionBlock( options ) ),
    timer_( new Util::Timer() )
{
  setDefaultOptions();

  setOptions( options );
}

//-----------------------------------------------------------------------------
// Function      : HBDirectSolver::~HBDirectSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
HBDirectSolver::~HBDirectSolver()
{
  delete timer_;
  delete options_;
}

//-----------------------------------------------------------------------------
// Function      : HBDirectSolver::setOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool HBDirectSolver::setOptions( const Util::OptionBlock & OB )
{
  Util::ParamList::const_iterator it_tpL = OB.begin();
  Util::ParamList::const_iterator end_tpL = OB.end();
  for (; it_tpL != end_tpL; ++it_tpL)
  {
    setParam( *it_tpL );
  }

  if ( solver_ == "DEFAULT" )
  {
    solver_ = solverDefault_;
  }

#ifdef Xyce_AMESOS2_BASKER
  if ( solver_ != "LAPACK" && solver_ != "BASKER" && solver_ != "BLOCK_BASKER" )
#else
  if ( solver_ != "LAPACK" )
#endif
  {
    Report::UserWarning0()
        << "HBDirectSolver does not recognize solver type " << solver_ << " setting to LAPACK";
    solver_ = "LAPACK";
  }

  if( options_ ) delete options_;
  options_ = new Util::OptionBlock( OB );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBDirectSolver::setDefaultOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool HBDirectSolver::setDefaultOptions()
{
  solver_ = solverDefault_;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBDirectSolver::setParam
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool HBDirectSolver::setParam( const Util::Param & param )
{
  std::string tag = param.tag();
  std::string uTag = param.uTag();

  if( uTag == "DIRECT_SOLVER" )
    solver_ = param.usVal();

  if( uTag == "OUTPUT_LS" ) 
    outputLS_ = param.getImmutableValue<int>();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBDirectSolver::doSolve
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
int HBDirectSolver::doSolve( bool reuse_factors, bool transpose )
{
  // Start the timer...
  timer_->resetStartTime();

  int linearStatus = 0;

  if (!isInit_)
  {
    // Get the number of frequencies and the number of time-domain unknowns.
    N_ = freqs_.size();
    M_ = (int)((N_-1)/2);

    if (hbOsc_)
    {
      numAugRows_ = (hbBuilderPtr_->getAugmentedLIDs()).size();
    }
    n_ = ((lasProblem_.getRHS())->globalLength()-numAugRows_) / (2*N_);

    // Create the block structure of the HB Jacobian.
    createBlockStructures();

    isInit_ = true;
  }

  double begAssembleTime = timer_->elapsedTime();

  // Generate the block HB Jacobian.
  formHBJacobian();

  if (VERBOSE_LINEAR)
  {
    double endAssembleTime = timer_->elapsedTime();
    Xyce::dout() << "Assembly Time: "
                 << (endAssembleTime - begAssembleTime) << std::endl;
  }

  // Output linear system, if requested.
  static int file_number = 1;
  if (outputLS_) 
  {
    if (!(file_number % outputLS_)) 
    {
      char file_name[40];
      sprintf( file_name, "Base_HB_Matrix%d.mm", file_number );
      printHBJacobian( std::string( file_name ) ); 
      sprintf( file_name, "Base_HB_RHS%d.mm", file_number );
      printHBResidual( std::string( file_name ) );
    }
  }

  double begNumTime = timer_->elapsedTime();

  // Factor the block HB Jacobian.
  linearStatus = numericFactorization();

  if (VERBOSE_LINEAR)
  {
    double endNumTime = timer_->elapsedTime();
    Xyce::dout() << "Numeric Factorization Time: "
                 << (endNumTime - begNumTime) << std::endl;
  }
   
  if (linearStatus != 0) {

    // Inform user that singular matrix was found and linear solve has failed.
    Report::UserWarning0() 
      << "Numerically singular matrix found by " << solver_ << ", returning zero solution to nonlinear solver!";

    // Put zeros in the solution since  was not able to solve this problem
    (lasProblem_.getLHS())->putScalar( 0.0 );

    return linearStatus;  // return the actual status (see bug 414 SON)
  }

  // Perform linear solve using factorization
  double begSolveTime = timer_->elapsedTime();

  linearStatus = solve();

  if (VERBOSE_LINEAR)
  {
    double endSolveTime = timer_->elapsedTime();
    Xyce::dout() << "Solve Time: "
                 << (endSolveTime - begSolveTime) << std::endl;
  }

  // Update the total solution time
  solutionTime_ = timer_->elapsedTime();

  if (VERBOSE_LINEAR)
    Xyce::dout() << "Total Linear Solution Time: "
                 << solutionTime_ << std::endl;

  // Write out solution, if requested.
  if (outputLS_) 
  {
    if (!(file_number % outputLS_)) 
    {
      char file_name[40];
      sprintf( file_name, "Base_HB_Soln%d.mm", file_number );
      printHBSolution( std::string( file_name ) );

    }
    file_number++;
  }

  return 0;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void HBDirectSolver::createBlockStructures()
{
  int numProcs = (builder_.getPDSComm())->numProc();
  int myProc = (builder_.getPDSComm())->procID();

  // Allocate space for the solution and RHS vectors. 
  if (solver_ == "LAPACK" || solver_ == "BASKER")
  {
    X_.reshape(N_*n_, 1);
    B_.reshape(N_*n_, 1);
  }
  else if (solver_ == "BLOCK_BASKER")
  {
    bX_.resize(n_, HBBlockMatrixEntry( N_, 1, true ));
    bB_.resize(n_, HBBlockMatrixEntry( N_, 1, true ));
  }

  if (solver_ == "LAPACK")
  {
    // Create a dense matrix that is the right size for the HB jacobian
    denseHBJacobian_.rows = N_*n_;
    denseHBJacobian_.cols = N_*n_;
    denseHBJacobian_.denseMtx.reshape(N_*n_, N_*n_);

    if (myProc == 0)
    {
      lapackSolver_ = Teuchos::rcp( new Teuchos::SerialDenseSolver<int,std::complex<double> >() );
    }
  }
  else if (solver_ == "BASKER" || solver_ == "BLOCK_BASKER")
  {

    Teuchos::RCP<Matrix> parMatrix;
    Teuchos::RCP<Vector> parVector;
    const Parallel::ParMap * columnMapPtr, * rowMapPtr;
    if (numProcs > 1)
    {
      parMatrix = Teuchos::rcp( builder_.createMatrix() );
      parVector = Teuchos::rcp( builder_.createVector() );

      columnMapPtr = parMatrix->getColMap( *builder_.getPDSComm() );
      rowMapPtr = parVector->pmap();
    }

    // Get the separated stored Jacobian matrices from the HB loader.
    Teuchos::RCP<FilteredMatrix>& linAppdQdx = hbLoaderPtr_->getStoreLindQdx();
    Teuchos::RCP<FilteredMatrix>& linAppdFdx = hbLoaderPtr_->getStoreLindFdx();
    std::vector<Teuchos::RCP<Linear::FilteredMatrix> >& nonlinAppdQdx = hbLoaderPtr_->getStoreNLdQdx();
    std::vector<Teuchos::RCP<Linear::FilteredMatrix> >& nonlinAppdFdx = hbLoaderPtr_->getStoreNLdFdx();
 
    // Determine the number of unique unknowns for each row.
    // The linear and nonlinear entries are kept separate and will be collected later.
    std::vector< std::vector<int> > nnzCol( n_ ), nnzCol_nl( n_ );

    // Collect linear nonzeros.
    const std::vector<int>& lindQdxNZs = linAppdQdx->getNZRows();
    const std::vector<int>& lindQdxRowPtr = linAppdQdx->getRowPtr();
    const std::vector<int>& lindQdxIndices = linAppdQdx->getIndices();
    for (std::vector<int>::const_iterator iter = lindQdxNZs.begin(); iter != lindQdxNZs.end(); iter++)
    {
      int row = *iter;
      int parRow = row;
      if (numProcs > 1)
      {
        parRow = rowMapPtr->localToGlobalIndex( row );
      }

      int numCols = lindQdxRowPtr[row+1] - lindQdxRowPtr[row];
      for (int j = 0; j < numCols; j++)
      {
        int col = lindQdxIndices[ lindQdxRowPtr[row] + j ];
        if (numProcs > 1)
        {
          col = columnMapPtr->localToGlobalIndex( col );
        }
        nnzCol[col].push_back( parRow );
      }
    }

    const std::vector<int>& lindFdxNZs = linAppdFdx->getNZRows();
    const std::vector<int>& lindFdxRowPtr = linAppdFdx->getRowPtr();
    const std::vector<int>& lindFdxIndices = linAppdFdx->getIndices();
    for (std::vector<int>::const_iterator iter = lindFdxNZs.begin(); iter != lindFdxNZs.end(); iter++)
    {
      int row = *iter;
      int parRow = row;
      if (numProcs > 1)
      {
        parRow = rowMapPtr->localToGlobalIndex( row );
      }

      int numCols = lindFdxRowPtr[row+1] - lindFdxRowPtr[row];
      for (int j = 0; j < numCols; j++)
      {
        int col = lindFdxIndices[ lindFdxRowPtr[row] + j ];
        if (numProcs > 1)
        {
          col = columnMapPtr->localToGlobalIndex( col );
        }
        nnzCol[col].push_back( parRow );
      }
    }

    // Get the frequency-domain Jacobian matrix entries from the HB loader.
    const std::vector< std::vector< Util::FreqMatEntry > >& linFreqdFdx = 
      hbLoaderPtr_->getFreqDFDXMatrix();
    if (linFreqdFdx.size())
    {
      for (std::vector< Util::FreqMatEntry >::const_iterator iter = linFreqdFdx[0].begin();
           iter != linFreqdFdx[0].end(); iter++)
      {
        int row = (*iter).row_lid;
        int col = (*iter).col_lid;
        int parRow = row;
        if (numProcs > 1)
        {
          parRow = rowMapPtr->localToGlobalIndex( row );
          col = columnMapPtr->localToGlobalIndex( col );
        }
        if (col != -1)
        {
          nnzCol[col].push_back( parRow );
        }
      }
    }

    int size_nldFdx = nonlinAppdFdx.size(); 
    std::vector<int> lin_nldFdx_ptr, lin_nldFdx_ind;
    std::vector<double> lin_nldFdx_val;
    if (size_nldFdx)
    {
      lin_nldFdx_ptr = nonlinAppdFdx[0]->getRowPtr();
      lin_nldFdx_ind = nonlinAppdFdx[0]->getIndices(); 
      lin_nldFdx_val = nonlinAppdFdx[0]->getValues();
    }

    for (int i=0; i < size_nldFdx; i++)
    {
      const std::vector<int>& nonlindFdxNZs = nonlinAppdFdx[i]->getNZRows();
      const std::vector<int>& nonlindFdxRowPtr = nonlinAppdFdx[i]->getRowPtr();
      const std::vector<int>& nonlindFdxIndices = nonlinAppdFdx[i]->getIndices();
      const std::vector<double>& nonlindFdxValues = nonlinAppdFdx[i]->getValues();

      for (std::vector<int>::const_iterator iter = nonlindFdxNZs.begin(); iter != nonlindFdxNZs.end(); iter++)
      {
        int row = *iter;
        int parRow = row;
        if (numProcs > 1)
        {
          parRow = rowMapPtr->localToGlobalIndex( row );
        }

        int numCols = nonlindFdxRowPtr[row+1] - nonlindFdxRowPtr[row];
        for (int j = 0; j < numCols; j++)
        {
          int col = nonlindFdxIndices[ nonlindFdxRowPtr[row] + j ];
          double val = nonlindFdxValues[ nonlindFdxRowPtr[row] + j ];

          for (int ptr=lin_nldFdx_ptr[row]; ptr<lin_nldFdx_ptr[row+1]; ptr++)
          {
            if ( col == lin_nldFdx_ind[ptr] )
            {
              // If the values changed, then change index to something invalid.
              if ( val != lin_nldFdx_val[ptr] )
              {
                lin_nldFdx_ind[ptr] = -1;
              }
            }
          } 
          if (numProcs > 1)
          { 
            col = columnMapPtr->localToGlobalIndex( col );     
          }
          nnzCol_nl[col].push_back( parRow );
        }
      }
    }

    for (std::vector<int>::const_iterator iter = nonlinAppdFdx[0]->getNZRows().begin();
         iter != nonlinAppdFdx[0]->getNZRows().end(); iter++)
    {
      for (int ptr=lin_nldFdx_ptr[*iter]; ptr<lin_nldFdx_ptr[*iter+1]; ptr++)
      { 
        if ( lin_nldFdx_ind[ptr] != -1 )
        {
          lin_nldFdx_.insert( std::make_pair( *iter, lin_nldFdx_ind[ptr] ) ); 
        }
      }
    }

    int size_nldQdx = nonlinAppdQdx.size(); 
    std::vector<int> lin_nldQdx_ptr, lin_nldQdx_ind;
    std::vector<double> lin_nldQdx_val;
    if (size_nldQdx)
    {
      lin_nldQdx_ptr = nonlinAppdQdx[0]->getRowPtr();
      lin_nldQdx_ind = nonlinAppdQdx[0]->getIndices(); 
      lin_nldQdx_val = nonlinAppdQdx[0]->getValues();
    }

    for (int i=0; i < size_nldQdx; i++)
    {
      const std::vector<int>& nonlindQdxNZs = nonlinAppdQdx[i]->getNZRows();
      const std::vector<int>& nonlindQdxRowPtr = nonlinAppdQdx[i]->getRowPtr();
      const std::vector<int>& nonlindQdxIndices = nonlinAppdQdx[i]->getIndices();
      const std::vector<double>& nonlindQdxValues = nonlinAppdQdx[i]->getValues();

      for (std::vector<int>::const_iterator iter = nonlindQdxNZs.begin(); iter != nonlindQdxNZs.end(); iter++)
      {
        int row = *iter;
        int parRow = row;
        if (numProcs > 1)
        {
          parRow = rowMapPtr->localToGlobalIndex( row );
        }

        int numCols = nonlindQdxRowPtr[row+1] - nonlindQdxRowPtr[row];
        for (int j = 0; j < numCols; j++)
        {
          int col = nonlindQdxIndices[ nonlindQdxRowPtr[row] + j ];
          double val = nonlindQdxValues[ nonlindQdxRowPtr[row] + j ];
          
          for (int ptr=lin_nldQdx_ptr[row]; ptr<lin_nldQdx_ptr[row+1]; ptr++)
          {
            if ( col == lin_nldQdx_ind[ptr] )
            {
              // If the values changed, then change index to something invalid.
              if ( val != lin_nldQdx_val[ptr] )
              {
                lin_nldQdx_ind[ptr] = -1;
              }
            }
          }
          if (numProcs > 1)
          { 
            col = columnMapPtr->localToGlobalIndex( col );        
          }
          nnzCol_nl[col].push_back( parRow );
        }
      }
    }

    for (std::vector<int>::const_iterator iter = nonlinAppdQdx[0]->getNZRows().begin();
         iter != nonlinAppdQdx[0]->getNZRows().end(); iter++)
    {
      for (int ptr=lin_nldQdx_ptr[*iter]; ptr<lin_nldQdx_ptr[*iter+1]; ptr++)
      {   
        if ( lin_nldQdx_ind[ptr] != -1 )
        {
          lin_nldQdx_.insert( std::make_pair( *iter, lin_nldQdx_ind[ptr] ) ); 
        }
      }
    }

    // Now compute the unique column ids.
    for (int i=0; i < n_; i++)
    {
      // First remove duplicates from nnzCol_nl.
      std::sort( nnzCol_nl[i].begin(), nnzCol_nl[i].end() );
      nnzCol_nl[i].erase(std::unique(nnzCol_nl[i].begin(), nnzCol_nl[i].end()), nnzCol_nl[i].end());

      // Then remove duplicates from nnzCol.
      std::sort( nnzCol[i].begin(), nnzCol[i].end() );
      nnzCol[i].erase(std::unique(nnzCol[i].begin(), nnzCol[i].end()), nnzCol[i].end());
    }
    
    // Communicate the nonzero entries to processor zero.
    if (numProcs > 1)
    {
      for (int i=0; i < n_; i++)
      {
        // Each processor should send nonzero values to proc 0 for column i.
        for (int proc = 1; proc < numProcs; proc++)
        {
          // Each processor, other than zero, should send their nonzero row indices.
          if ( myProc == proc )
          {
            // Send linear row indices.
            int nnzs = nnzCol[i].size();
            (builder_.getPDSComm())->send( &nnzs, 1, 0 );
            if (nnzs)
            {
              (builder_.getPDSComm())->send( &(nnzCol[i][0]), nnzs, 0 );
            }  
            // Send nonlinear row indices. 
            nnzs = nnzCol_nl[i].size();
            (builder_.getPDSComm())->send( &nnzs, 1, 0 );
            if (nnzs)
            {
              (builder_.getPDSComm())->send( &(nnzCol_nl[i][0]), nnzs, 0 );
            }  
          }
          if ( myProc == 0 )
          {
            // Receive linear row indices.
            int nnzs = 0;
            (builder_.getPDSComm())->recv( &nnzs, 1, proc );
            if (nnzs)
            {
              std::vector<int> rowIndices( nnzs );
              (builder_.getPDSComm())->recv( &rowIndices[0], nnzs, proc );
              nnzCol[i].insert( nnzCol[i].end(), rowIndices.begin(), rowIndices.end() );
            }
            // Receive nonlinear row indices.
            (builder_.getPDSComm())->recv( &nnzs, 1, proc );
            if (nnzs)
            {
              std::vector<int> rowIndices( nnzs );
              (builder_.getPDSComm())->recv( &rowIndices[0], nnzs, proc );
              nnzCol_nl[i].insert( nnzCol_nl[i].end(), rowIndices.begin(), rowIndices.end() );
            }
          }
        }

        // Compute unique entries again on processor 0.
        if ( myProc == 0 )
        {
          std::sort( nnzCol[i].begin(), nnzCol[i].end() );
          nnzCol[i].erase(std::unique(nnzCol[i].begin(), nnzCol[i].end()), nnzCol[i].end());

          std::sort( nnzCol_nl[i].begin(), nnzCol_nl[i].end() );
          nnzCol_nl[i].erase(std::unique(nnzCol_nl[i].begin(), nnzCol_nl[i].end()), nnzCol_nl[i].end()); 
        }
      }
    } 
        
    Acol_ptr_.clear();
    Arow_idx_.clear();
    Aval_.clear();

    Acol_ptr_.push_back( 0 );
    for (int i=0; i < n_; i++)
    {
      // First insert nonlinear non-zero entries, then sort and compute unique entries.
      nnzCol[i].insert( nnzCol[i].end(), nnzCol_nl[i].begin(), nnzCol_nl[i].end() );
      std::sort( nnzCol[i].begin(), nnzCol[i].end() );
      nnzCol[i].erase(std::unique(nnzCol[i].begin(), nnzCol[i].end()), nnzCol[i].end());

      // Now fill out the block structure of Acol_ptr, Arow_idx, and Aval.
      Acol_ptr_.push_back( Acol_ptr_[i] + nnzCol[i].size() );
      for (unsigned int j=0; j < nnzCol[i].size(); j++)
      {
        Arow_idx_.push_back( nnzCol[i][j] );
      }
    }

    // Allocate block matrix entries of the HB Jacobian.
    // These are initialized to diagonal matrices and then corrected to dense matrices for nonlinear entries.
    Aval_.resize( Acol_ptr_[n_], HBBlockMatrixEntry( N_, N_, false ) );
    for (int i=0; i < n_; i++)
    {
      for (std::vector<int>::const_iterator iter=nnzCol_nl[i].begin(); iter != nnzCol_nl[i].end(); iter++)
      {
        for (int ptr=Acol_ptr_[i]; ptr<Acol_ptr_[i+1]; ptr++)
        {
          // Place a dense matrix anywhere there is a nonlinear device entry.
          if (Arow_idx_[ptr] == *iter && Aval_[ptr].isDiag())
          {
            Aval_[ptr].denseMtx.reshape( N_, N_ );
            Aval_[ptr].diagVector.clear();
          }
        }
      }  
    }
  }

  // Create serial objects for parallel
  if (numProcs > 1)
  {
    Teuchos::RCP<Parallel::EpetraParMap> e_solnMap = Teuchos::rcp_dynamic_cast<Parallel::EpetraParMap>(hbBuilderPtr_->getSolutionMap());
    serialMap_ = Teuchos::rcp( new Epetra_Map( Epetra_Util::Create_Root_Map( *(e_solnMap->petraMap()), 0 ) ) );
    serialImporter_ = Teuchos::rcp( new Epetra_Import( *serialMap_, *(e_solnMap->petraMap()) ) );
    serialX_ = Teuchos::rcp( new Epetra_Vector( *serialMap_ ) );
    serialB_ = Teuchos::rcp( new Epetra_Vector( *serialMap_ ) );
  }

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void HBDirectSolver::initializeBlockCRS( std::complex<double> val )
{
  // Initialize the dense or diagonal blocks to the input value.
  for (unsigned int i=0; i < Aval_.size(); i++)
  {
    Aval_[i].putScalar( val );
  }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void HBDirectSolver::formHBJacobian()
{
  int myProc = (builder_.getPDSComm())->procID();
  int numProcs = (builder_.getPDSComm())->numProc();

  // Get the separated stored Jacobian matrices from the HB loader.
  Teuchos::RCP<FilteredMatrix>& linAppdQdx = hbLoaderPtr_->getStoreLindQdx();
  Teuchos::RCP<FilteredMatrix>& linAppdFdx = hbLoaderPtr_->getStoreLindFdx();
  std::vector<Teuchos::RCP<Linear::FilteredMatrix> >& nonlinAppdQdx = hbLoaderPtr_->getStoreNLdQdx();
  std::vector<Teuchos::RCP<Linear::FilteredMatrix> >& nonlinAppdFdx = hbLoaderPtr_->getStoreNLdFdx();

  int posFreq = (N_-1)/2;
  std::complex<double> IMAG(0.0,1.0);

  if ( solver_ == "LAPACK" )
  { 
    // Initialize values of HB jacobian to zero.
    denseHBJacobian_.putScalar( std::complex<double>( 0.0, 0.0 ) );
  }
  else if ( solver_ == "BASKER" || solver_ == "BLOCK_BASKER" )
  {
    initializeBlockCRS( std::complex<double>( 0.0, 0.0 ) );
  }

  Teuchos::RCP<Matrix> parMatrix;
  Teuchos::RCP<Vector> parVector;
  const Parallel::ParMap * columnMapPtr = 0, * rowMapPtr = 0;
  if (numProcs > 1)
  {
    parMatrix = Teuchos::rcp( builder_.createMatrix() );
    parVector = Teuchos::rcp( builder_.createVector() );

    columnMapPtr = parMatrix->getColMap( *builder_.getPDSComm() );
    rowMapPtr = parVector->pmap();
  }

  // load nonlin dFdx
  int size_nldFdx = nonlinAppdFdx.size(); 
  Teuchos::BLAS<int,std::complex<double> > blas;
  Teuchos::RCP< Teuchos::SerialDenseMatrix<int, std::complex<double> > > blockView;
  std::set<std::pair<int,int> > dFdx_set;

  // first check if any of the nonlinear dFdx entries is linear.
  std::set<std::pair<int, int> > curr_lin_nldFdx( lin_nldFdx_ );

  if ( !curr_lin_nldFdx.empty() )
  {
    std::set<std::pair<int, int> >::iterator it = curr_lin_nldFdx.begin();
    while ( it != curr_lin_nldFdx.end() )
    {
      int row = it->first;
      int col = it->second;
      double val = 0.0;

      bool removePair = true;
      for (int i=0; i < size_nldFdx; i++)
      {
        const std::vector<int>& nonlindFdxRowPtr = nonlinAppdFdx[i]->getRowPtr();
        const std::vector<int>& nonlindFdxIndices = nonlinAppdFdx[i]->getIndices();
        const std::vector<double>& nonlindFdxValues = nonlinAppdFdx[i]->getValues();

        bool isFound = false, isSame = false;
        for (int ptr = nonlindFdxRowPtr[row]; ptr < nonlindFdxRowPtr[row+1]; ptr++)
        {
          if (nonlindFdxIndices[ptr] == col)
          {
            if (i==0)
            {
              val = nonlindFdxValues[ptr];
              isSame = true;
            }
            else if (nonlindFdxValues[ptr] == val)
            {
              isSame = true;
            }
            else
            {
              isSame = false;
              break;
            }
            isFound = true;
          }
        }

        if (isFound && isSame)
        {
          removePair = false;
        }
        else
        {
          removePair = true;
          break;
        }
      }

      if (removePair)
      {
        std::set<std::pair<int, int> >::iterator next_it = it;
        next_it++;
        curr_lin_nldFdx.erase( it );
        it = next_it;
      }
      else
      {
        it++;
      }
    }
  } 

  for (int i=0; i < size_nldFdx; i++)
  {
    const std::vector<int>& nonlindFdxNZs = nonlinAppdFdx[i]->getNZRows();
    const std::vector<int>& nonlindFdxRowPtr = nonlinAppdFdx[i]->getRowPtr();
    const std::vector<int>& nonlindFdxIndices = nonlinAppdFdx[i]->getIndices();
    const std::vector<double>& nonlindFdxValues = nonlinAppdFdx[i]->getValues();

    Teuchos::SerialDenseMatrix<int, std::complex<double> > W_i( N_, 1 );
    if ( (nonlindFdxIndices.size()-curr_lin_nldFdx.size()) )
    {
      for (int j=0; j<=M_; j++)
      {
        std::complex<double> x( 0, -2.0 * M_PI * freqs_[posFreq+j] * times_[i] );
        W_i(j,0) = std::exp( x );
 
        if (j > 0)
        {
          W_i(N_-j,0) = std::exp( -x );
        } 
      }

      for (std::vector<int>::const_iterator iter = nonlindFdxNZs.begin(); iter != nonlindFdxNZs.end(); iter++)
      {
        int row = *iter;
        int parRow = row;
        if (numProcs > 1)
        {
          parRow = rowMapPtr->localToGlobalIndex( row );
        }
 
        int numCols = nonlindFdxRowPtr[row+1] - nonlindFdxRowPtr[row];
        for (int j = 0; j < numCols; j++)
        {
          int col = nonlindFdxIndices[ nonlindFdxRowPtr[row] + j ];
          
          // See if this entry is actually linear and will be loaded later.
          if ( curr_lin_nldFdx.find( std::make_pair( row, col ) ) != curr_lin_nldFdx.end() )
          {
            continue;
          }

          if (numProcs > 1)
          {
            col = columnMapPtr->localToGlobalIndex( col );
          }
          double val = nonlindFdxValues[ nonlindFdxRowPtr[row] + j ];

          // Collect this row and col.  
          dFdx_set.insert( std::make_pair( parRow, col ) );

          if ( solver_ == "LAPACK" )
          {
            blockView = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int, std::complex<double> >( Teuchos::View, 
                                      denseHBJacobian_.denseMtx, N_, N_, parRow*N_, col*N_ ) );
          }
          else if ( solver_ == "BASKER" || solver_ == "BLOCK_BASKER" )
          {
            for (int ptr = Acol_ptr_[col]; ptr < Acol_ptr_[col+1]; ptr++)
            {
              if ( Arow_idx_[ptr] == parRow )
              {
                blockView = Teuchos::rcp( &(Aval_[ptr].denseMtx), false );
                break;
              }
            }
          }
  
          // Add in contributions for time point i.
          blas.HERK( Teuchos::UPPER_TRI, Teuchos::NO_TRANS, N_, 1, val/N_, W_i.values(), W_i.stride(), 
                     std::complex<double>(1.0,0.0), blockView->values(), blockView->stride() );
        }
      }

      // Now copy the conjugate values over to the other side of the dense matrices.
      for (std::set<std::pair<int,int> >::const_iterator iter = dFdx_set.begin(); iter != dFdx_set.end(); iter++)
      {
        if ( solver_ == "LAPACK" )
        {
          blockView = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int, std::complex<double> >( Teuchos::View, 
                                        denseHBJacobian_.denseMtx, N_, N_, iter->first*N_, iter->second*N_ ) );
        }
        else if ( solver_ == "BASKER" || solver_ == "BLOCK_BASKER" )
        {
          for (int ptr = Acol_ptr_[iter->second]; ptr < Acol_ptr_[iter->second+1]; ptr++)
          {
            if ( Arow_idx_[ptr] == iter->first )
            {
              blockView = Teuchos::rcp( &(Aval_[ptr].denseMtx), false );
              break;
            }
          }
        }

        for (int i=0; i<N_; i++)
        {
          for (int j=0; j<i; j++)
          {
            (*blockView)(i,j) = Teuchos::ScalarTraits<std::complex<double> >::conjugate( (*blockView)(j,i) );
          }
        }
      }
    }
  }

  // load lin_nldFdx entries!
  if ( !curr_lin_nldFdx.empty() )
  {
    std::set<std::pair<int, int> >::iterator it = curr_lin_nldFdx.begin();
    while ( it != curr_lin_nldFdx.end() )
    {
      int row = it->first;
      int col = it->second;
      double val = 0.0;

      // Get the value.
      const std::vector<int>& nonlindFdxRowPtr = nonlinAppdFdx[0]->getRowPtr();
      const std::vector<int>& nonlindFdxIndices = nonlinAppdFdx[0]->getIndices();
      const std::vector<double>& nonlindFdxValues = nonlinAppdFdx[0]->getValues();
      for (int ptr = nonlindFdxRowPtr[row]; ptr < nonlindFdxRowPtr[row+1]; ptr++)
      {
        if (nonlindFdxIndices[ptr] == col)
        {
          val = nonlindFdxValues[ptr];
          break;
        }
      }

      // Enter the values into the appropriate block entry.
      int parRow = row;
      if (numProcs > 1)
      {
        parRow = rowMapPtr->localToGlobalIndex( row );
        col = columnMapPtr->localToGlobalIndex( col );
      }

      if ( solver_ == "LAPACK" )
      {
        blockView = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int, std::complex<double> >( Teuchos::View, denseHBJacobian_.denseMtx,
                                                                          N_, N_, parRow*N_, col*N_) );

        for (int ii = 0; ii < N_; ii++)
        {
          (*blockView)(ii, ii) += val;
        }
      }
      else if ( solver_ == "BASKER" || solver_ == "BLOCK_BASKER" )
      {
        for (int ptr = Acol_ptr_[col]; ptr < Acol_ptr_[col+1]; ptr++)
        {
          if ( Arow_idx_[ptr] == parRow )
          {
            for (int ii = 0; ii < N_; ii++)
            {
              Aval_[ptr].addToDiag( ii, val );
            }
          }
        }
      }

      it++;
    }
  }

  // load nonlin dQdx
  int size_nldQdx = nonlinAppdQdx.size(); 
  std::set<std::pair<int,int> > dQdx_set;
  Teuchos::SerialDenseMatrix<int, std::complex<double> > nonlindQdx_lapack;
  std::vector<Xyce::HBBlockMatrixEntry> nonlindQdx_basker;

  // first check if any of the nonlinear dQdx entries is linear.
  std::set<std::pair<int, int> > curr_lin_nldQdx( lin_nldQdx_ );

  if ( !curr_lin_nldQdx.empty() )
  {
    std::set<std::pair<int, int> >::iterator it = curr_lin_nldQdx.begin();
    while ( it != curr_lin_nldQdx.end() )
    {
      int row = it->first;
      int col = it->second;
      double val = 0.0;

      bool removePair = true;
      for (int i=0; i < size_nldQdx; i++)
      {
        const std::vector<int>& nonlindQdxRowPtr = nonlinAppdFdx[i]->getRowPtr();
        const std::vector<int>& nonlindQdxIndices = nonlinAppdFdx[i]->getIndices();
        const std::vector<double>& nonlindQdxValues = nonlinAppdFdx[i]->getValues();

        bool isFound = false, isSame = false;
        for (int ptr = nonlindQdxRowPtr[row]; ptr < nonlindQdxRowPtr[row+1]; ptr++)
        {
          if (nonlindQdxIndices[ptr] == col)
          {
            if (i==0)
            {
              val = nonlindQdxValues[ptr];
              isSame = true;
            }
            else if (nonlindQdxValues[ptr] == val)
            {
              isSame = true;
            }
            else
            {
              isSame = false;
              break;
            }
            isFound = true;
          }
        }

        if (isFound && isSame)
        {
          removePair = false;
        }
        else
        {
          removePair = true;
          break;
        }
      }

      if (removePair)
      {
        std::set<std::pair<int, int> >::iterator next_it = it;
        next_it++;
        curr_lin_nldQdx.erase( it );
        it = next_it;
      }
      else
      {
        it++;
      }
    }
  }

  if ( solver_ == "LAPACK" )
  {
    nonlindQdx_lapack.reshape( denseHBJacobian_.rows, denseHBJacobian_.cols );
  }
  else if ( solver_ == "BASKER" || solver_ == "BLOCK_BASKER" )
  {
    nonlindQdx_basker.resize( Acol_ptr_[n_], HBBlockMatrixEntry( N_, N_, true ) );
  }

  for (int i=0; i < size_nldQdx; i++)
  {
    const std::vector<int>& nonlindQdxNZs = nonlinAppdQdx[i]->getNZRows();
    const std::vector<int>& nonlindQdxRowPtr = nonlinAppdQdx[i]->getRowPtr();
    const std::vector<int>& nonlindQdxIndices = nonlinAppdQdx[i]->getIndices();
    const std::vector<double>& nonlindQdxValues = nonlinAppdQdx[i]->getValues();

    Teuchos::SerialDenseMatrix<int, std::complex<double> > W_i( N_, 1 );
    if ( (nonlindQdxIndices.size()-curr_lin_nldQdx.size()) > 0 )
    {
      for (int j=0; j<=M_; j++)
      {
        std::complex<double> x( 0, -2.0 * M_PI * freqs_[posFreq+j] * times_[i] );
        W_i(j,0) = std::exp( x );
  
        if (j > 0)
        {
          W_i(N_-j,0) = std::exp( -x );
        }
      }

      for (std::vector<int>::const_iterator iter = nonlindQdxNZs.begin(); iter != nonlindQdxNZs.end(); iter++)
      {
        int row = *iter;
        int parRow = row;
        if (numProcs > 1)
        {
          parRow = rowMapPtr->localToGlobalIndex( row );
        }

        int numCols = nonlindQdxRowPtr[row+1] - nonlindQdxRowPtr[row];
        for (int j = 0; j < numCols; j++)
        {
          int col = nonlindQdxIndices[ nonlindQdxRowPtr[row] + j ];

          if ( curr_lin_nldQdx.find( std::make_pair( row, col ) ) != curr_lin_nldQdx.end() )
          {
            continue;
          }

          if (numProcs > 1)
          {
            col = columnMapPtr->localToGlobalIndex( col );
          }
          double val = nonlindQdxValues[ nonlindQdxRowPtr[row] + j ];

          // Insert row and col for post processing.
          dQdx_set.insert( std::make_pair( parRow, col ) );

          if ( solver_ == "LAPACK" )
          {
            blockView = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int, std::complex<double> >( Teuchos::View, nonlindQdx_lapack,
                                                                            N_, N_, parRow*N_, col*N_ ) );
          }
          else if ( solver_ == "BASKER" || solver_ == "BLOCK_BASKER" )
          {
            for (int ptr = Acol_ptr_[col]; ptr < Acol_ptr_[col+1]; ptr++)
            {
              if ( Arow_idx_[ptr] == parRow )
              {
                blockView = Teuchos::rcp( &(nonlindQdx_basker[ptr].denseMtx), false );
                break;
              }
            }
          }

          // Add in contributions for time point i.
          blas.HERK( Teuchos::UPPER_TRI, Teuchos::NO_TRANS, N_, 1, val/N_, W_i.values(), W_i.stride(),
                     std::complex<double>(1.0,0.0), blockView->values(), blockView->stride() );
        }
      }
    }  
  }

  // Now copy the conjugate values over to the other side of the dense matrices.
  for (std::set<std::pair<int,int> >::const_iterator iter = dQdx_set.begin(); iter != dQdx_set.end(); iter++)
  {
    if ( solver_ == "LAPACK" )
    {
      blockView = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int, std::complex<double> >( Teuchos::View, nonlindQdx_lapack,
                                                                    N_, N_, iter->first*N_, iter->second*N_ ) );
    }
    else if ( solver_ == "BASKER" || solver_ == "BLOCK_BASKER" )
    {
      for (int ptr = Acol_ptr_[iter->second]; ptr < Acol_ptr_[iter->second+1]; ptr++)
      {
        if ( Arow_idx_[ptr] == iter->first )
        {
          blockView = Teuchos::rcp( &(nonlindQdx_basker[ptr].denseMtx), false );
          break;
        }
      }
    }

    // Symmetrize.
    for (int i=0; i<N_; i++)
    {
      // First symmetrize.
      for (int j=0; j<i; j++)
      {
        (*blockView)(i,j) = Teuchos::ScalarTraits<std::complex<double> >::conjugate( (*blockView)(j,i) );
      }
    }

    // Now scale.
    for (int i=0; i<=M_; i++)
    {
      std::complex<double> cplxCoeff = IMAG * 2.0 * M_PI * freqs_[posFreq+i]; 
      for (int j=0; j<N_; j++)
      {
        (*blockView)(i,j) *= cplxCoeff;

        if (i > 0)
        {
          (*blockView)(N_-i, j) *= -cplxCoeff;
        }
      }
    }

    // Now update dense Jacobian with the contribution from dQdx
    Teuchos::RCP< Teuchos::SerialDenseMatrix<int, std::complex<double> > > finalBlockView;

    if ( solver_ == "LAPACK" )
    {
      finalBlockView = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int, std::complex<double> >( Teuchos::View, 
                                         denseHBJacobian_.denseMtx, N_, N_, iter->first*N_, iter->second*N_ ) );
    }
    else if ( solver_ == "BASKER" || solver_ == "BLOCK_BASKER" )
    {
      for (int ptr = Acol_ptr_[iter->second]; ptr < Acol_ptr_[iter->second+1]; ptr++)
      {
        if ( Arow_idx_[ptr] == iter->first )
        {
          finalBlockView = Teuchos::rcp( &(Aval_[ptr].denseMtx), false );
          break;
        }
      }
    }  

    // Sum the dQdx portion of the Jacobian into the dFdx portion.
    (*finalBlockView) += (*blockView);
  }

   // load lin_nldQdx entries!
  if ( !curr_lin_nldQdx.empty() )
  {
    std::set<std::pair<int, int> >::iterator it = curr_lin_nldQdx.begin();
    while ( it != curr_lin_nldQdx.end() )
    {
      int row = it->first;
      int col = it->second;
      double val = 0.0;

      // Get the value.
      const std::vector<int>& nonlindQdxRowPtr = nonlinAppdQdx[0]->getRowPtr();
      const std::vector<int>& nonlindQdxIndices = nonlinAppdQdx[0]->getIndices();
      const std::vector<double>& nonlindQdxValues = nonlinAppdQdx[0]->getValues();
      for (int ptr = nonlindQdxRowPtr[row]; ptr < nonlindQdxRowPtr[row+1]; ptr++)
      {
        if (nonlindQdxIndices[ptr] == col)
        {
          val = nonlindQdxValues[ptr];
          break;
        }
      }

      // Enter the values into the appropriate block entry.
      int parRow = row;
      if (numProcs > 1)
      {
        parRow = rowMapPtr->localToGlobalIndex( row );
        col = columnMapPtr->localToGlobalIndex( col );
      }

      if ( solver_ == "LAPACK" )
      {
        blockView = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int, std::complex<double> >( Teuchos::View, denseHBJacobian_.denseMtx,
                                                                          N_, N_, parRow*N_, col*N_ ) );

        for (int ii = 0; ii <= M_; ii++)
        {
          std::complex<double> cplxCoeff = IMAG * 2.0 * M_PI * freqs_[posFreq+ii];

          (*blockView)(ii, ii) += cplxCoeff * val;

          if (ii > 0)
          {
            (*blockView)(N_-ii, N_-ii) -= cplxCoeff * val;
          }
        }
      }
      else if ( solver_ == "BASKER" || solver_ == "BLOCK_BASKER" )
      {
        for (int ptr = Acol_ptr_[col]; ptr < Acol_ptr_[col+1]; ptr++)
        {
          if ( Arow_idx_[ptr] == parRow )
          {
            for (int ii = 0; ii <= M_; ii++)
            {
              std::complex<double> cplxCoeff = IMAG * 2.0 * M_PI * freqs_[posFreq+ii];

              Aval_[ptr].addToDiag( ii, cplxCoeff * val );

              if (ii > 0)
              {
                Aval_[ptr].addToDiag( N_-ii, -cplxCoeff * val );
              }
            }
          }
        }
      }

      it++;
    }
  }
 
/*
  std::cout << "Nonlinear Dense HB Jacobian Matrix: " << std::endl;
  denseHBJacobian_.denseMtx.print(std::cout);
*/

  // load linear dQdx
  const std::vector<int>& lindQdxNZs = linAppdQdx->getNZRows();
  const std::vector<int>& lindQdxRowPtr = linAppdQdx->getRowPtr();
  const std::vector<int>& lindQdxIndices = linAppdQdx->getIndices();
  const std::vector<double>& lindQdxValues = linAppdQdx->getValues();

  for (std::vector<int>::const_iterator iter = lindQdxNZs.begin(); iter != lindQdxNZs.end(); iter++)
  {
    int row = *iter;
    int parRow = row;
    if (numProcs > 1)
    {
      parRow = rowMapPtr->localToGlobalIndex( row );
    }

    int numCols = lindQdxRowPtr[row+1] - lindQdxRowPtr[row];
    for (int j = 0; j < numCols; j++)
    {
      int col = lindQdxIndices[ lindQdxRowPtr[row] + j ];
      if (numProcs > 1)
      {
        col = columnMapPtr->localToGlobalIndex( col );
      }
      double val = lindQdxValues[ lindQdxRowPtr[row] + j ];

      if ( solver_ == "LAPACK" )
      {
        blockView = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int, std::complex<double> >( Teuchos::View, denseHBJacobian_.denseMtx,
                                                                          N_, N_, parRow*N_, col*N_ ) );

        for (int ii = 0; ii <= M_; ii++)
        {
          std::complex<double> cplxCoeff = IMAG * 2.0 * M_PI * freqs_[posFreq+ii];

          (*blockView)(ii, ii) += cplxCoeff * val;

          if (ii > 0)
          {
            (*blockView)(N_-ii, N_-ii) -= cplxCoeff * val;
          }
        }
      }
      else if ( solver_ == "BASKER" || solver_ == "BLOCK_BASKER" )
      {
        for (int ptr = Acol_ptr_[col]; ptr < Acol_ptr_[col+1]; ptr++)
        {
          if ( Arow_idx_[ptr] == parRow )
          {
            for (int ii = 0; ii <= M_; ii++)
            {
              std::complex<double> cplxCoeff = IMAG * 2.0 * M_PI * freqs_[posFreq+ii];

              Aval_[ptr].addToDiag( ii, cplxCoeff * val );

              if (ii > 0)
              {
                Aval_[ptr].addToDiag( N_-ii, -cplxCoeff * val );
              }
            }
          }
        }
      }
    }
  }

  // load dFdx
  const std::vector<int>& lindFdxNZs = linAppdFdx->getNZRows();
  const std::vector<int>& lindFdxRowPtr = linAppdFdx->getRowPtr();
  const std::vector<int>& lindFdxIndices = linAppdFdx->getIndices();
  const std::vector<double>& lindFdxValues = linAppdFdx->getValues();

  for (std::vector<int>::const_iterator iter = lindFdxNZs.begin(); iter != lindFdxNZs.end(); iter++)
  {
    int row = *iter;
    int parRow = row;
    if (numProcs > 1)
    {
      parRow = rowMapPtr->localToGlobalIndex( row );
    }

    int numCols = lindFdxRowPtr[row+1] - lindFdxRowPtr[row];
    for (int j = 0; j < numCols; j++)
    {
      int col = lindFdxIndices[ lindFdxRowPtr[row] + j ];
      if (numProcs > 1)
      {
        col = columnMapPtr->localToGlobalIndex( col );
      }
      double val = lindFdxValues[ lindFdxRowPtr[row] + j ];

      if ( solver_ == "LAPACK" )
      {
        blockView = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int, std::complex<double> >( Teuchos::View, denseHBJacobian_.denseMtx,
                                                                          N_, N_, parRow*N_, col*N_) );

        for (int ii = 0; ii < N_; ii++)
        {
          (*blockView)(ii, ii) += val;
        }
      }
      else if ( solver_ == "BASKER" || solver_ == "BLOCK_BASKER" )
      {
        for (int ptr = Acol_ptr_[col]; ptr < Acol_ptr_[col+1]; ptr++)
        {
          if ( Arow_idx_[ptr] == parRow )
          {
            for (int ii = 0; ii < N_; ii++)
            {
              Aval_[ptr].addToDiag( ii, val );
            }
          }
        }
      }
    }
  }


  // load frequency-domain dFdx, which also contains dQdx contributions.
  const std::vector< std::vector< Util::FreqMatEntry > >& linFreqdFdx = 
    hbLoaderPtr_->getFreqDFDXMatrix();
  if ( linFreqdFdx.size() )
  {
    std::vector< Util::FreqMatEntry >::const_iterator iter = linFreqdFdx[0].begin();
    for (int entry=0; iter != linFreqdFdx[0].end(); iter++, entry++)
    {
      int row = (*iter).row_lid;
      int col = (*iter).col_lid;
      int parRow = row;
      if (numProcs > 1)
      {
        parRow = rowMapPtr->localToGlobalIndex( row );
        col = columnMapPtr->localToGlobalIndex( col );
      }

      if ( solver_ == "LAPACK" )
      {
        blockView = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int, std::complex<double> >( Teuchos::View, denseHBJacobian_.denseMtx,
                                                                          N_, N_, parRow*N_, col*N_ ) );
        for (int ii = 0; ii <= M_; ii++)
        {
          std::complex<double> val = linFreqdFdx[ii][entry].val;

          (*blockView)(ii, ii) += val;

          if (ii > 0)
          {
            (*blockView)(N_-ii, N_-ii) += Teuchos::ScalarTraits<std::complex<double> >::conjugate( val );
          }
        }
      }
      else if ( solver_ == "BASKER" || solver_ == "BLOCK_BASKER" )
      {
        for (int ptr = Acol_ptr_[col]; ptr < Acol_ptr_[col+1]; ptr++)
        {
          if ( Arow_idx_[ptr] == parRow )
          {
            for (int ii = 0; ii <= M_; ii++)
            {
              std::complex<double> val = linFreqdFdx[ii][entry].val;
  
              Aval_[ptr].addToDiag( ii, val );
  
              if (ii > 0)
              {
                Aval_[ptr].addToDiag( N_-ii, Teuchos::ScalarTraits<std::complex<double> >::conjugate( val ) );
              }
            }
          }
        }
      }
    }
  }

  // Now collect the global matrix on processor 0 for solving.
  if (numProcs > 1)
  {
    if ( solver_ == "LAPACK" )
    {
      int numRows = denseHBJacobian_.rows;
      int numCols = denseHBJacobian_.cols;
      std::vector< double > sendRecvVec( 2 * numRows * numCols );
      for (int proc = 1; proc < numProcs; proc++)
      {
        if ( myProc == proc )
        {
          Xyce::packHBBlockMatrix( denseHBJacobian_, sendRecvVec );
          (builder_.getPDSComm())->send( &sendRecvVec[0], 2*numRows*numCols, 0 );
        }

        if ( myProc == 0 )
        {
          (builder_.getPDSComm())->recv( &sendRecvVec[0], 2*numRows*numCols, proc );
          Xyce::unpackHBBlockMatrixUpdate( sendRecvVec, 1, denseHBJacobian_ );
        }
      }
    }
    else
    {
      // Create a vector that is twice the dimension of a dense block.
      std::vector<int> rcd( 3, 0 );
      std::vector< double > sendRecvVec( 2 * N_ * N_ );
      for (int proc = 1; proc < numProcs; proc++ )
      {
        int localNZs = 0;
        if ( myProc == proc )
        {
          // Communicate how many nonzero entries you are sending to processor 0.
          int localNZs = Acol_ptr_[n_];
          (builder_.getPDSComm())->send( &localNZs, 1, 0 );
          
          // Now pack each block and send them one by one.
          for ( int j=0; j < n_; j++ )
          {
            for (int ptr = Acol_ptr_[j]; ptr < Acol_ptr_[j+1]; ptr++)
            {
              // Send row, column, and density.
              rcd[0] = Arow_idx_[ptr]; rcd[1] = j; rcd[2] = Aval_[ptr].isDense();
              (builder_.getPDSComm())->send( &rcd[0], 3, 0 );

              // Send values.
              int len = Xyce::packHBBlockMatrix( Aval_[ptr], sendRecvVec );
              (builder_.getPDSComm())->send( &sendRecvVec[0], len, 0 );
            } 
          }
        }
        if ( myProc == 0 )
        {
          // Get the number of nonzero entries to be received from processor 'proc'.
          (builder_.getPDSComm())->recv( &localNZs, 1, proc );
          for (int blocks = 0; blocks < localNZs; blocks++)
          {
            // Get row, column, and dense.
            (builder_.getPDSComm())->recv( &rcd[0], 3, proc );

            // Compute length of value vector.
            int len = 2 * N_;
            if ( rcd[2] )
            {
              len *= N_;
            }

            // Receive values.
            (builder_.getPDSComm())->recv( &sendRecvVec[0], len, proc );

            // Get block entry from local Aval_ vector
            for (int ptr = Acol_ptr_[ rcd[1] ]; ptr < Acol_ptr_[ rcd[1]+1 ]; ptr++)
            {
              if ( rcd[0] == Arow_idx_[ ptr ] )
              {
                Xyce::unpackHBBlockMatrixUpdate( sendRecvVec, rcd[2], Aval_[ptr] );
              }
            }
          }
        }
      }
    }      
  }
  
  if ( solver_ == "LAPACK" )
  {
/*
    if ( (builder_.getPDSComm())->procID() == 0 )
    {
      std::cout << "Final Dense HB Jacobian Matrix: " << std::endl;
      denseHBJacobian_.print(std::cout);
    }
*/
  } 
  else if ( solver_ == "BASKER" || outputLS_ )
  {
    // Copy over matrix into non-block data structures for Basker.
    if (myProc == 0)
    {
      // Copy matrix from Acol_ptr_, Arow_idx_, and Aval_.
      Anewcol_ptr_.clear();
      Anewrow_idx_.clear();
      Anewval_.clear();
      Anewcol_ptr_.push_back( 0 );
      for (int i=0; i<n_; i++)
      {
        for (int idx=0; idx<N_; idx++)
        {
          for (int j=Acol_ptr_[i]; j<Acol_ptr_[i+1]; j++)
          {
            int blockRow = Arow_idx_[j];

            if ( Aval_[j].isDiag() )
            {
              Anewval_.push_back( Aval_[j].diagVector[ idx ] );
              Anewrow_idx_.push_back( blockRow*N_ + idx );
            }
            else
            {
              for (int row=0; row<N_; row++)
              {
                std::complex<double> val = Aval_[j].denseMtx(row, idx);
                if (val != Teuchos::ScalarTraits<std::complex<double> >::zero())
                {
                  Anewval_.push_back( val );
                  Anewrow_idx_.push_back( blockRow*N_ + row );
                }
              }
            }
          }
          Anewcol_ptr_.push_back( Anewrow_idx_.size() );
        }
      }
    }
  }

  // form frequency-domain RHS vector
  MultiVector* B = lasProblem_.getRHS();
  int numVectors = B->numVectors();
  EpetraVectorAccess* e_B = dynamic_cast<EpetraVectorAccess *>( B );

  for (int j=0; j<numVectors; j++)
  {
    Teuchos::RCP<Linear::Vector> B_j;

    if (numProcs > 1)
    {
      serialB_->Import( *((e_B->epetraObj())( j )), *serialImporter_, Insert );
      B_j = Teuchos::rcp( new EpetraVector( &*serialB_, *serialMap_, false ) );
    }
    else
    {
      B_j = Teuchos::rcp( B->getNonConstVectorView( j ) );
    }

    if ( myProc == 0 )
    {

      for (int nB = 0; nB < n_; nB++)
      {
        // Copy values from B_j to B_ vector for solver.
        for (int i = 0; i <= M_; i++)
        {
          std::complex<double> val( (*B_j)[nB*2*N_ + 2*i], (*B_j)[nB*2*N_ + 2*i+1] );

          if (solver_ == "LAPACK" || solver_ == "BASKER")
            B_( nB*N_ + i, 0 ) = val;
          else
            bB_[nB].denseMtx(i,0) = val;

          if (i > 0)
          {
            if (solver_ == "LAPACK" || solver_ == "BASKER")
              B_( (nB+1)*N_ - i, 0 ) = Teuchos::ScalarTraits<std::complex<double> >::conjugate( val );
            else
              bB_[nB].denseMtx(N_-i,0) = Teuchos::ScalarTraits<std::complex<double> >::conjugate( val );
          }
        }
      }
    }
  }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int HBDirectSolver::numericFactorization()
{
  int linearStatus = 0;

  if ( (builder_.getPDSComm())->procID() == 0 )
  {
    if ( solver_ == "LAPACK" )
    {
      // Compute numeric factorization of dense Jacobian matrix.
      if (DEBUG_LINEAR)
      {
        A_ = denseHBJacobian_.denseMtx;
      }
      lapackSolver_->setMatrix( Teuchos::rcp( &(denseHBJacobian_.denseMtx), false ) );
      lapackSolver_->setVectors( Teuchos::rcp( &X_, false ), Teuchos::rcp( &B_, false ) );
      lapackSolver_->factorWithEquilibration(true);
      linearStatus = lapackSolver_->factor();
    }
#ifdef Xyce_AMESOS2_BASKER
    else if ( solver_ == "BASKER" )
    {
      // Create Basker solver and factor block diagonal matrix.
      basker_.factor(N_*n_, N_*n_, Anewcol_ptr_[N_*n_], &Anewcol_ptr_[0], &Anewrow_idx_[0], &Anewval_[0]);
    }
    else if ( solver_ == "BLOCK_BASKER" )
    {
      // Create Basker solver and factor block diagonal matrix.
      blockBasker_.factor(n_, n_, Acol_ptr_[n_], &Acol_ptr_[0], &Arow_idx_[0], &(Aval_[0]));
    }
#endif
  }

  // Get the return code to all processors.
  (builder_.getPDSComm())->bcast( &linearStatus, 1, 0 );

  return linearStatus;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int HBDirectSolver::solve()
{
  int linearStatus = 0;

  // Determine number of time-domain variables.
  MultiVector* X = lasProblem_.getLHS();
  MultiVector* B = lasProblem_.getRHS();
  EpetraVectorAccess* e_X = dynamic_cast<EpetraVectorAccess *>( X );

  // Initialize solution vector X.
  X->putScalar( 0.0 );

  int numVectors = X->numVectors();
  int numProcs = (builder_.getPDSComm())->numProc();
  int myProc = (builder_.getPDSComm())->procID();

  for (int j=0; j<numVectors; j++)
  {
    Teuchos::RCP<Linear::Vector> X_j;

    if (numProcs > 1)
    {
      X_j = Teuchos::rcp( new EpetraVector( &*serialX_, *serialMap_, false ) ); 
    }
    else
    {
      X_j = Teuchos::rcp( X->getNonConstVectorView( j ) );
    } 

    if ( myProc == 0 )
    {
      if ( solver_ == "LAPACK" )
      {
        // Solve the dense Jacobian in the frequency domain, complex-valued.
        double bnorm, rnorm;
        if (DEBUG_LINEAR)
        {
          R_ = B_;
          bnorm = R_.normFrobenius();
        }

        linearStatus = lapackSolver_->solve();

        if (DEBUG_LINEAR)
        {
          R_.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, std::complex<double>(-1.0,0.0), A_, X_, std::complex<double>(1.0,0.0) );
          rnorm = R_.normFrobenius();
          Xyce::dout() << "Linear System Residual (LAPACK) : " << rnorm/bnorm << std::endl; 
        }
      }
#ifdef Xyce_AMESOS2_BASKER
      else if ( solver_ == "BASKER" )
      {
        double bnorm, rnorm;
        if (DEBUG_LINEAR)
        {
          bnorm = B_.normFrobenius();
        }

        basker_.solve(B_.values(), X_.values());

        if (DEBUG_LINEAR)
        {
          for (int j=0; j<n_*N_; j++)
          {
            for (int ptr = Anewcol_ptr_[j]; ptr < Anewcol_ptr_[j+1]; ptr++)
            {
              B_(Anewrow_idx_[ptr],0) -= Anewval_[ptr]*X_(j,0);
            }
          }

          rnorm = B_.normFrobenius();
          if (bnorm > 0.0)
            Xyce::dout() << "Linear System Residual (BASKER) : " << rnorm/bnorm << std::endl; 
          else
            Xyce::dout() << "Linear System Residual (BASKER) : " << rnorm << std::endl;
        }
      }
      else if ( solver_ == "BLOCK_BASKER" )
      {
        double bnorm = 0.0, rnorm = 0.0;
        if (DEBUG_LINEAR)
        {
          for (int j=0; j<n_; j++)
          {
            double bjnorm = bB_[j].normFrobenius();
            bnorm += ( bjnorm*bjnorm );
          }
          bnorm = Teuchos::ScalarTraits<double>::magnitude(
                    Teuchos::ScalarTraits<double>::squareroot( bnorm ) );
        }

        blockBasker_.solve(&bB_[0], &bX_[0]);

        if (DEBUG_LINEAR)
        {
          for (int j=0; j<n_; j++)
          {
            for (int ptr = Acol_ptr_[j]; ptr < Acol_ptr_[j+1]; ptr++)
            {
              bB_[Arow_idx_[ptr]] -= Aval_[ptr]*bX_[j];
            }
          }

          for (int j=0; j<n_; j++)
          {
            double bjnorm = bB_[j].normFrobenius();
            Xyce::dout() << "Residual norm of block " << j << " is " << bjnorm << std::endl;
            rnorm += ( bjnorm*bjnorm );
          }
          rnorm = Teuchos::ScalarTraits<double>::magnitude(
                    Teuchos::ScalarTraits<double>::squareroot( rnorm ) );
          if ( bnorm > 0.0 )
            Xyce::dout() << "Linear System Residual (BLOCK BASKER) : " << rnorm/bnorm << std::endl; 
          else
            Xyce::dout() << "Linear System Residual (BLOCK BASKER) : " << rnorm << std::endl;
        }
      }
#endif

      for (int nB = 0; nB < n_; nB++)
      {
        // Copy values from X_ (from solver) to X_j.
        for (int i = 0; i <= M_; i++)
        {
          if (solver_ == "LAPACK" || solver_ == "BASKER")
          {
            (*X_j)[nB*2*N_ + 2*i] = X_(nB*N_ + i, 0).real();
            (*X_j)[nB*2*N_ + 2*i+1] = X_(nB*N_ + i, 0).imag();
          }
          else
          {
            (*X_j)[nB*2*N_ + 2*i] = bX_[nB].denseMtx(i,0).real();
            (*X_j)[nB*2*N_ + 2*i+1] = bX_[nB].denseMtx(i,0).imag();
          }

          if (i > 0)
          {
            if (solver_ == "LAPACK" || solver_ == "BASKER")
            {
              (*X_j)[(nB+1)*2*N_ - 2*i] = X_(nB*N_ + i, 0).real();
              (*X_j)[(nB+1)*2*N_ - 2*i+1] = - X_(nB*N_ + i, 0).imag();
            }
            else
            {
              (*X_j)[(nB+1)*2*N_ - 2*i] = bX_[nB].denseMtx(i,0).real();
              (*X_j)[(nB+1)*2*N_ - 2*i+1] = - bX_[nB].denseMtx(i,0).imag();
            }
          }
        }
      }
    }

    if (numProcs > 1)
    {
      (e_X->epetraObj())( j )->Export( *serialX_, *serialImporter_, Add );
    }
  }

  // Get the return code to all processors.
  (builder_.getPDSComm())->bcast( &linearStatus, 1, 0 );

  return linearStatus;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void HBDirectSolver::printHBJacobian( const std::string& fileName )
{
  int myProc = (builder_.getPDSComm())->procID();

  if (myProc == 0)
  {
    std::ofstream out;
    out.open( fileName.c_str() );

    // Write out banner.
    out << "%%MatrixMarket matrix ";

    if (solver_ == "LAPACK")
    {
      // Output dense format.
      out << "array complex general" << std::endl;

      // Write the dimensions of the sparse matrix: (# rows, #
      // columns, # matrix entries (counting duplicates as
      // separate entries)).
      out << denseHBJacobian_.rows << " " << denseHBJacobian_.cols << std::endl;

      // Set precision.
      out.precision( 16 );
      out << std::scientific;

      for (int j=0; j<denseHBJacobian_.cols; j++)
      {
        for (int i=0; i<denseHBJacobian_.rows; i++)
        {
          out << denseHBJacobian_.denseMtx( i, j ).real() << " " << denseHBJacobian_.denseMtx( i, j ).imag() << std::endl;
        }
      }
    }
  
    if (solver_ == "BASKER" || solver_ == "BLOCK_BASKER")
    {
      // Output sparse format.
      out << "coordinate complex general" << std::endl;

      // Write the dimensions of the sparse matrix: (# rows, #
      // columns, # matrix entries (counting duplicates as
      // separate entries)).
      out << n_*N_ << " " << n_*N_ << " " << Anewcol_ptr_[n_*N_] << std::endl;

      // Set precision.
      out.precision( 16 );
 
      for (int j=0; j<n_*N_; j++)
      {
        for (int ptr=Anewcol_ptr_[j]; ptr<Anewcol_ptr_[j+1]; ptr++)
        {
          out << Anewrow_idx_[ptr]+1 << " " << j+1 << " " 
              << std::scientific << Anewval_[ptr].real() << " " << Anewval_[ptr].imag() 
              << std::resetiosflags(std::ios_base::floatfield) << std::endl;
        }
      }
    }
    // Close the file.
    out.close(); 
  }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void HBDirectSolver::printHBResidual( const std::string& fileName )
{
  int myProc = (builder_.getPDSComm())->procID();

  // Determine number of time-domain variables.
  MultiVector* B = lasProblem_.getRHS();
  int numVectors = B->numVectors();

  if ( myProc == 0 )
  {
    std::ofstream out;
    out.open( fileName.c_str() );

    // Write out banner.
    out << "%%MatrixMarket matrix array complex general" << std::endl;

    // Write the dimensions of the sparse matrix: (# rows, #
    // columns, # matrix entries (counting duplicates as
    // separate entries)).
    out << n_*N_ << " " << numVectors << std::endl;

    // Set precision.
    out.precision( 16 );
    out << std::scientific;

    for (int j=0; j<numVectors; j++)
    {
      for (int nB = 0; nB < n_; nB++)
      {
        for (int i = 0; i < N_; i++)
        {
          if (solver_ == "LAPACK" || solver_ == "BASKER")
            out << B_( nB*N_ + i, 0 ).real() << " " << B_( nB*N_ + i, 0 ).imag() << std::endl;
          else
            out << bB_[nB].denseMtx(i,0).real() << " " << bB_[nB].denseMtx(i,0).imag() << std::endl;
        }
      }
    }
    // Close the file.
    out.close();
  }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void HBDirectSolver::printHBSolution( const std::string& fileName )
{
  int myProc = (builder_.getPDSComm())->procID();

  // Determine number of time-domain variables.
  MultiVector* X = lasProblem_.getLHS();
  int numVectors = X->numVectors();

  if ( myProc == 0 )
  {
    std::ofstream out;
    out.open( fileName.c_str() );

    // Write out banner.
    out << "%%MatrixMarket matrix array complex general" << std::endl;

    // Write the dimensions of the sparse matrix: (# rows, #
    // columns, # matrix entries (counting duplicates as
    // separate entries)).
    out << n_*N_ << " " << numVectors << std::endl;

    // Set precision.
    out.precision( 16 );
    out << std::scientific;

    for (int j=0; j<numVectors; j++)
    {
      for (int nB = 0; nB < n_; nB++)
      {
        for (int i = 0; i < N_; i++)
        {
          if (solver_ == "LAPACK" || solver_ == "BASKER")
          {
            out << X_(nB*N_ + i, 0).real() << " " << X_(nB*N_ + i, 0).imag() << std::endl;
          }
          else
          {
            out << bX_[nB].denseMtx(i,0).real() << " " << bX_[nB].denseMtx(i,0).imag() << std::endl;
          }
        }
      }
    }
    // Close the file.
    out.close();
  }
}

} // namespace Linear
} // namespace Xyce

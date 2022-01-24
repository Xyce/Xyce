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
// Purpose        : ES direct solver wrapper
// Special Notes  :
// Creator        : Eric Keiter, SNL
// Creation Date  : 6/1/2018
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Xyce Includes ----------

#include <N_UTL_fwd.h>

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_ESDirectSolver.h>
#include <N_LAS_ESBuilder.h>
#include <N_LOA_ESLoader.h>
#include <N_LAS_Builder.h>
#include <N_LAS_Solver.h>
#include <N_LAS_Problem.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_EpetraVector.h>
#include <N_LAS_EpetraHelpers.h>
#include <N_LAS_FilteredMatrix.h>
#include <N_LAS_TransformTool.h>
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
// Function      : ESDirectSolver::ESDirectSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//-----------------------------------------------------------------------------
ESDirectSolver::ESDirectSolver(
  Builder &       builder,
  Problem &       problem,
  Util::OptionBlock &   options)
  : Solver(problem, false),
    builder_(builder),
    isInit_(false),
    N_(0),
    n_(0),
    outputLS_(0),
    solver_(""),
    solverDefault_("LAPACK"),
    options_( new Util::OptionBlock( options ) ),
    timer_( new Util::Timer() ),
    numSamples_(1),
    paramsOuterLoop_(true)
{
  setDefaultOptions();

  setOptions( options );
}

//-----------------------------------------------------------------------------
// Function      : ESDirectSolver::~ESDirectSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//-----------------------------------------------------------------------------
ESDirectSolver::~ESDirectSolver()
{
  delete timer_;
  delete options_;
}

//-----------------------------------------------------------------------------
// Function      : ESDirectSolver::setOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//-----------------------------------------------------------------------------
bool ESDirectSolver::setOptions( const Util::OptionBlock & OB )
{
  Util::ParamList::const_iterator it_tpL = OB.begin();
  Util::ParamList::const_iterator end_tpL = OB.end();
  for (; it_tpL != end_tpL; ++it_tpL)
  {
    setParam( *it_tpL );
  }
  
  if (solver_ == "")
    solver_ = solverDefault_;

#ifdef Xyce_AMESOS2_BASKER
  if ( solver_ != "LAPACK" && solver_ != "BLOCK_BASKER" )
#else
  if ( solver_ != "LAPACK" )
#endif
  {
    Report::UserWarning0()
        << "ESDirectSolver does not recognize solver type " << solver_ << " setting to LAPACK";
    solver_ = "LAPACK";
  }

  if( options_ ) delete options_;
  options_ = new Util::OptionBlock( OB );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ESDirectSolver::setDefaultOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//-----------------------------------------------------------------------------
bool ESDirectSolver::setDefaultOptions()
{
  solver_ = solverDefault_;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ESDirectSolver::setParam
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//-----------------------------------------------------------------------------
bool ESDirectSolver::setParam( const Util::Param & param )
{
  std::string tag = param.tag();
  std::string uTag = param.uTag();

  if( uTag == "TYPE" ) 
    solver_ = param.usVal();

  if( uTag == "OUTPUT_LS" ) 
    outputLS_ = param.getImmutableValue<int>();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ESDirectSolver::doSolve
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//-----------------------------------------------------------------------------
int ESDirectSolver::doSolve( bool reuse_factors, bool transpose )
{
  // Start the timer...
  timer_->resetStartTime();
  int linearStatus = 0;
  if (!isInit_)
  {
    // Get the number of samples and the number of unknowns.
    N_ = numSamples_;
    n_ = (lasProblem_.getRHS())->globalLength() / (N_);

    // Create the block structure of the ES Jacobian.
    createBlockStructures();

    isInit_ = true;
  }

#if 0
  {
  {
    std::cout << "Acol_ptr_ array:" << std::endl;
    int size=Acol_ptr_.size();
    for (int ii=0;ii<size;++ii)
    {
      std::cout << "Acol_ptr_["<<ii<<"] = " << Acol_ptr_[ii] << std::endl;
    }
  }

  {
    std::cout << "Arow_idx_ array:" << std::endl;
    int size=Arow_idx_.size();
    for (int ii=0;ii<size;++ii)
    {
      std::cout << "Arow_idx_["<<ii<<"] = " << Arow_idx_[ii] << std::endl;
    }
  }

  {
    std::cout << "Aval_ array:" << std::endl;
    int size=Aval_.size();
    for (int ii=0;ii<size;++ii)
    {
      std::cout << "Aval_["<<ii<<"] = " << Aval_[ii] << std::endl;
    }
  }
  }
#endif


  double begAssembleTime = timer_->elapsedTime();

  // Generate the block ES Jacobian.
  formESJacobian();

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
      sprintf( file_name, "Base_ES_Matrix%d.mm", file_number );
      printESJacobian( std::string( file_name ) ); 
      sprintf( file_name, "Base_ES_RHS%d.mm", file_number );
      printESResidual( std::string( file_name ) );
    }
  }

  double begNumTime = timer_->elapsedTime();

  // Factor the block ES Jacobian.
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
      sprintf( file_name, "Base_ES_Soln%d.mm", file_number );
      printESSolution( std::string( file_name ) );
    }
    file_number++;
  }

  return 0;
}

//---------------------------------------------------------------------------
// Function      : ESDirectSolver::createBlockStructures
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//---------------------------------------------------------------------------
void ESDirectSolver::createBlockStructures()
{
  int numProcs = (builder_.getPDSComm())->numProc();
  int myProc = (builder_.getPDSComm())->procID();

  // Allocate space for the solution and RHS vectors. 
  if (solver_ == "LAPACK")
  {
    X_.reshape(N_*n_, 1);
    B_.reshape(N_*n_, 1);
  }
  else if (solver_ == "BLOCK_BASKER")
  {
    bool isDense=true; // check this
    if (paramsOuterLoop_) // not implemented yet for block_basker
    {
      bX_.resize(N_, ESBlockMatrixEntry( n_, 1, isDense ));
      bB_.resize(N_, ESBlockMatrixEntry( n_, 1, isDense ));
    }
    else
    {
      bX_.resize(n_, ESBlockMatrixEntry( N_, 1, isDense ));
      bB_.resize(n_, ESBlockMatrixEntry( N_, 1, isDense ));
    }
  }

  if (solver_ == "LAPACK") // structure is same for either ordering
  {
    // Create a dense matrix that is the right size for the ES jacobian
    denseESJacobian_.rows = N_*n_;
    denseESJacobian_.cols = N_*n_;
    denseESJacobian_.denseMtx.reshape(N_*n_, N_*n_);

    if (myProc == 0)
    {
      lapackSolver_ = Teuchos::rcp( new Teuchos::SerialDenseSolver<int, double>() );
    }
  }
  else if (solver_ == "BLOCK_BASKER") 
  {
    if (paramsOuterLoop_) // not implemented yet for block_basker
    {
      Report::UserFatal0() << "This ordering is not supported for the specialized BLOCK_BASKER solver" <<std::endl;
    }
    else
    {
      Teuchos::RCP<Matrix> parMatrix = Teuchos::rcp( builder_.createMatrix() );
      Teuchos::RCP<Vector> parVector = Teuchos::rcp( builder_.createVector() );

      const Parallel::ParMap * columnMapPtr = parMatrix->getColMap( *builder_.getPDSComm() );
      const Parallel::ParMap * rowMapPtr = parVector->pmap();

      // Determine the number of unique unknowns for each row.
      // This code is inspired by the similar code in the HBDirectSolver that Heidi 
      // wrote.  However, it is much simpler for several reasons.  
      // (1) I don't have to treat the linear and nonlinear parts separately (at least not so far; I may need to consider this later) 
      // (2) as such I am not using filtered matrices
      // (3) I also don't need to treat F and Q separately
      // (4) so the total amount of code is much smaller
      std::vector< std::vector<int> > nnzCol( n_ );
      Matrix * Jac = lasProblem_.getMatrix();
      BlockMatrix * bJac =  dynamic_cast<BlockMatrix*>(Jac); 
      Matrix & subMat = bJac->block(0,0); 

      for (int row=0;row<n_;++row) // loop over ckt unknowns
      {
        int parRow = row;
        if (numProcs > 1)
        {
          parRow = rowMapPtr->localToGlobalIndex( row );
        }

        int lengthRef = subMat.getLocalRowLength(row);
        int length=0; 
        double * coeffs; 
        int * colIndices;
        subMat.extractLocalRowView(row, length, coeffs, colIndices);

        for (int icol=0;icol<length;++icol)
        {
          int col = colIndices[icol];
          if (numProcs > 1)
          {
            col = columnMapPtr->localToGlobalIndex( col );
          }
          nnzCol[col].push_back( parRow );
        }
      }

      // Now compute the unique column ids.
      // Note, ERK: for ES this step probably isn't needed.
      for (int i=0; i < n_; i++)
      {
        // remove duplicates from nnzCol.
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
              // Send row indices.
              int nnzs = nnzCol[i].size();
              (builder_.getPDSComm())->send( &nnzs, 1, 0 );
              if (nnzs)
              {
                (builder_.getPDSComm())->send( &(nnzCol[i][0]), nnzs, 0 );
              }  
            }
            if ( myProc == 0 )
            {
              // Receive row indices.
              int nnzs = 0;
              (builder_.getPDSComm())->recv( &nnzs, 1, proc );
              if (nnzs)
              {
                std::vector<int> rowIndices( nnzs );
                (builder_.getPDSComm())->recv( &rowIndices[0], nnzs, proc );
                nnzCol[i].insert( nnzCol[i].end(), rowIndices.begin(), rowIndices.end() );
              }
            }
          }

          // Compute unique entries again on processor 0.
          if ( myProc == 0 )
          {
            std::sort( nnzCol[i].begin(), nnzCol[i].end() );
            nnzCol[i].erase(std::unique(nnzCol[i].begin(), nnzCol[i].end()), nnzCol[i].end());
          }
        }
      }

      Acol_ptr_.clear();
      Arow_idx_.clear();
      Aval_.clear();

      Acol_ptr_.push_back( 0 );
      for (int i=0; i < n_; i++)
      {
        // sort and compute unique NZ entries.
        // ERK; this may not be neccessary
        std::sort( nnzCol[i].begin(), nnzCol[i].end() );
        nnzCol[i].erase(std::unique(nnzCol[i].begin(), nnzCol[i].end()), nnzCol[i].end());

        // Now fill out the block structure of Acol_ptr, Arow_idx, and Aval.
        Acol_ptr_.push_back( Acol_ptr_[i] + nnzCol[i].size() );
        for (unsigned int j=0; j < nnzCol[i].size(); j++)
        {
          Arow_idx_.push_back( nnzCol[i][j] );
        }
      }

      // Allocate block matrix entries of the ES Jacobian
      // Unlike in the HB case, all of these are all diagonals
      //bool isDense=true;
      bool isDense=false; // experiment with dense
      Aval_.resize( Acol_ptr_[n_], ESBlockMatrixEntry( N_, N_, isDense ) );

#if 0
      // debug outputs
      {
        std::cout << "nnzCol 2D array:" << std::endl;
        int size=nnzCol.size();
        for (int ii=0;ii<size;++ii)
        {
          int sizeJJ=nnzCol[ii].size();
          for (int jj=0;jj<sizeJJ;++jj)
          {
            std::cout << "nnzCol["<<ii<<"]["<<jj<<"] = " << nnzCol[ii][jj] <<std::endl;
          }
        }
      }
#endif
    }
  }
  else
  {
    Report::UserWarning0() << "Solver type not recognized.  Using LAPACK" <<std::endl;
    // Create a dense matrix that is the right size for the ES jacobian
    denseESJacobian_.rows = N_*n_;
    denseESJacobian_.cols = N_*n_;
    denseESJacobian_.denseMtx.reshape(N_*n_, N_*n_);
    if (myProc == 0)
    {
      lapackSolver_ = Teuchos::rcp( new Teuchos::SerialDenseSolver<int, double>() );
    }
  }

  // Create serial objects for parallel
  // ERK. Not ready for this yet.
  if (numProcs > 1)
  {
    // ERK note, currently the builder that is passed into this class is an ES builder.  Not sure if this is 
    // correct builder for these function calls.  Check this.
    Teuchos::RCP<Parallel::EpetraParMap> e_solnMap = Teuchos::rcp_dynamic_cast<Parallel::EpetraParMap>(esBuilderPtr_->getSolutionMap());
    serialMap_ = Teuchos::rcp( new Epetra_Map( Epetra_Util::Create_Root_Map( *(e_solnMap->petraMap()), 0 ) ) );
    serialImporter_ = Teuchos::rcp( new Epetra_Import( *serialMap_, *(e_solnMap->petraMap()) ) );
    serialX_ = Teuchos::rcp( new Epetra_Vector( *serialMap_ ) );
    serialB_ = Teuchos::rcp( new Epetra_Vector( *serialMap_ ) );
  }
}

//---------------------------------------------------------------------------
// Function      : ESDirectSolver::initializeBlockCRS
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//---------------------------------------------------------------------------
void ESDirectSolver::initializeBlockCRS( double val )
{
  // Initialize the dense or diagonal blocks to the input value.
  for (unsigned int i=0; i < Aval_.size(); i++)
  {
    Aval_[i].putScalar( val );
  }
}

//---------------------------------------------------------------------------
// Function      : ESDirectSolver::formESJacobian
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//---------------------------------------------------------------------------
void ESDirectSolver::formESJacobian()
{
  int myProc = (builder_.getPDSComm())->procID();
  int numProcs = (builder_.getPDSComm())->numProc();

  Matrix * Jac = lasProblem_.getMatrix();
  BlockMatrix * bJac =  dynamic_cast<BlockMatrix*>(Jac); 
  Matrix & subMatRef = bJac->block(0,0); // use this to get the structure

  int numLocalRowsRef = subMatRef.getLocalNumRows(); // num ckt vars = n_
  int numBlockRows = bJac->numBlockRows(); // = num params = N_

  if (n_ != numLocalRowsRef)
  {
    Report::UserFatal0() << "n_ != numLocalRows. numLocalRows = " << numLocalRowsRef << "  n_ = " << n_ << " N_ = " << N_ <<std::endl;
  }

  if (N_ != numBlockRows)
  {
    Report::UserFatal0() << "N_ != numBlockRows" <<std::endl;
  }

  if ( solver_ == "LAPACK" )
  { 
    // Initialize values of ES jacobian to zero.
    denseESJacobian_.putScalar( 0.0 );
  }
  else if ( solver_ == "BLOCK_BASKER" )
  {
    initializeBlockCRS( double(0.0) );
  }

  int blockSize = bJac->blockSize();
  int numLocalRows = Jac->getLocalNumRows();

  // copy the matrix that was formed in the loader into the dense matrix used by LAPACK here.
  if (numProcs > 1)
  {
    Report::UserFatal0() << "Specialized ES solver not set up for parallel yet" <<std::endl;
  }

  if (paramsOuterLoop_)
  {
    // use the same ordering as the ES loader
    if ( solver_ == "LAPACK" )
    {
      for(int ii=0;ii<numLocalRows; ++ii)
      {
        int length=0;
        double * coeffs;
        int * colIndices;
        Jac->extractLocalRowView(ii, length, coeffs, colIndices);
        for (int icol=0;icol<length;++icol)
        {
          int jj=colIndices[icol];
          denseESJacobian_.denseMtx(ii,jj) = coeffs[icol];
        }
      }
    }
    else if ( solver_ == "BLOCK_BASKER" )
    {
      Report::UserFatal0() << "Specialized ES solver for conventional ordering not set up for BLOCK_BASKER yet" <<std::endl;
    }
  }
  else
  {
    // use inverted ordering (params inner loop, ckt unknowns outer loop)
    for (int ii=0;ii<n_;++ii) // loop over ckt unknowns
    {
      int lengthRef = subMatRef.getLocalRowLength(ii);

      for (int ipar=0;ipar<N_;++ipar) // loop over the paramters.  
      {
        Matrix & subMat = bJac->block(ipar,ipar);
        int length=0; 
        double * coeffs; 
        int * colIndices;
        subMat.extractLocalRowView(ii, length, coeffs, colIndices);
        if (length!=lengthRef) 
        { 
          // this may get triggered in parallel
          Report::UserFatal0() << "Ack!  lengths in block matrix don't match!!!" <<std::endl; 
        }

        for (int icol=0;icol<length;++icol)
        {
          if ( solver_ == "LAPACK" )
          {
            // The block row and block col indices in the new matrix are the
            // same as the row,col indices of the original ckt matrix = (ii,jj)
            // Within the current block, there is a diagonal of length N_.
            int jj=colIndices[icol];
            int Row= ii*N_ + ipar;
            int Col= jj*N_ + ipar;
            denseESJacobian_.denseMtx(Row,Col) = coeffs[icol];
          }
          else if ( solver_ == "BLOCK_BASKER" )
          {
            int jj=colIndices[icol];
            int AvalIndex = Acol_ptr_[ii] + icol;

            if (AvalIndex >= Aval_.size() || AvalIndex < 0)
            {
              Report::UserFatal0() << "AvalIndex is too big = " << AvalIndex <<std::endl;
            }

            if (Aval_[AvalIndex].isDense())
            {
              Aval_[AvalIndex].denseMtx(ipar,ipar) = coeffs[icol]; // experiment
            }
            else
            { // this should be used for embedded sampling but both can work
              Aval_[AvalIndex].diagVector[ipar] = coeffs[icol]; // experiment
            }
          }
        }
      }
    }
  }

  if ( solver_ == "BLOCK_BASKER"  && outputLS_ )
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
                double val = Aval_[j].denseMtx(row, idx);
                if (val != 0.0)
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

#if 0
  if ( solver_ == "LAPACK" )
  {
    std::cout.precision(4); std::cout.setf(std::ios::scientific);
    std::cout << "Dense Matrix:" <<std::endl;
    std::cout << denseESJacobian_.denseMtx;
    std::cout << std::endl;
  }

  if ( solver_ == "BLOCK_BASKER" )
  {
    std::cout << "Aval_ array:" << std::endl;
    int size=Aval_.size();
    for (int ii=0;ii<size;++ii)
    {
      std::cout << "Aval_["<<ii<<"] = " << Aval_[ii] << std::endl;
    }
    //exit(0);
  }
#endif

  MultiVector* B = lasProblem_.getRHS();
  EpetraVectorAccess* e_B = dynamic_cast<EpetraVectorAccess *>( B );

  int numVectors = B->numVectors();
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
      int size = B_j->globalLength();

      if (paramsOuterLoop_)
      {
        for (int ii=0;ii<size;++ii)
        {
          // Copy values from B_j to B_ vector for solver.
          if (solver_ == "LAPACK")
          {
              B_( ii, 0 ) =  (*B_j)[ii];
          }
          else
          {
            //bB_[nB].denseMtx(i,0) = val;
          }
        }
      }
      else // invert the ordering
      {
        for (int ii=0;ii<n_;++ii)
        {
          for (int ipar=0;ipar<N_;++ipar) // loop over the paramters.  
          {
            // Copy values from B_j to B_ vector for solver.
            if (solver_ == "BASKER")
            {
              int B_Row   = ii*N_ + ipar;
              int B_j_Row = ipar*n_ + ii;
              B_( B_Row, 0 ) =  (*B_j)[B_j_Row];
            }
            else
            {
              int B_j_Row = ipar*n_ + ii;
              bB_[ii].denseMtx(ipar,0) =  (*B_j)[B_j_Row];;
            }
          }
        }
      }
#if 0
      {
        std::cout << "B vector" << std::endl;
        for (int ii=0;ii<size;++ii)
        {
          std::cout << "B_("<<ii<<",0) = " << B_(ii,0) <<std::endl;
        }
      }
#endif
    }
  }
}

//---------------------------------------------------------------------------
// Function      : ESDirectSolver::numericFactorization
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//---------------------------------------------------------------------------
int ESDirectSolver::numericFactorization()
{
  int linearStatus = 0;

  if ( (builder_.getPDSComm())->procID() == 0 )
  {
    if ( solver_ == "LAPACK" )
    {
      lapackSolver_->setMatrix( Teuchos::rcp( &(denseESJacobian_.denseMtx), false ) );
      lapackSolver_->setVectors( Teuchos::rcp( &X_, false ), Teuchos::rcp( &B_, false ) );
      lapackSolver_->factorWithEquilibration(true);
      linearStatus = lapackSolver_->factor();
    }
#ifdef Xyce_AMESOS2_BASKER
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
// Function      : ESDirectSolver::solve
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//---------------------------------------------------------------------------
int ESDirectSolver::solve()
{
  int linearStatus = 0;

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
        linearStatus = lapackSolver_->solve();
      }
#ifdef Xyce_AMESOS2_BASKER
      else if ( solver_ == "BLOCK_BASKER" )
      {
        blockBasker_.solve(&bB_[0], &bX_[0]);
      }
#endif

      int size = lasProblem_.getRHS()->globalLength();

      if (paramsOuterLoop_)
      {
        for (int ii=0;ii<size;++ii)
        {
          if (solver_ == "LAPACK")
          {
            (*X_j)[ii] = X_(ii,0);
          }
          else
          {
            //(*X_j)[nB*2*N_ + 2*i] = bX_[nB].denseMtx(i,0).real();
          }
        }
      }
      else // invert the ordering
      {
        for (int ii=0;ii<n_;++ii)
        {
          for (int ipar=0;ipar<N_;++ipar) // loop over the paramters.  
          {
            // Copy values from X_j to X_ vector for solver.
            if (solver_ == "LAPACK")
            {
              int X_Row   = ii*N_ + ipar;
              int X_j_Row = ipar*n_ + ii;
              (*X_j)[X_j_Row] =  X_( X_Row, 0 ) ;
            }
            else
            {
              int X_Row   = ii*N_ + ipar;
              int X_j_Row = ipar*n_ + ii;
              (*X_j)[X_j_Row] = bX_[ii].denseMtx(ipar,0);
            }
          }
        }
      }

#if 0
      {
        std::cout << "X vector" << std::endl;
        for (int ii=0;ii<size;++ii)
        {
          std::cout << "X_("<<ii<<",0) = " << X_(ii,0) <<std::endl;
        }
      }
#endif
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
// Function      : ESDirectSolver::printESJacobian
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//---------------------------------------------------------------------------
void ESDirectSolver::printESJacobian( const std::string& fileName )
{
  int myProc = (builder_.getPDSComm())->procID();

  std::ofstream out;
  out.open( fileName.c_str() );

  // Write out banner.
  out << "%%MatrixMarket matrix ";

  if (myProc == 0)
  {
    if (solver_ == "LAPACK")
    {
      // Output dense format.
      out << "array real general" << std::endl;

      // Write the dimensions of the sparse matrix: (# rows, #
      // columns, # matrix entries (counting duplicates as
      // separate entries)).
      out << denseESJacobian_.rows << " " << denseESJacobian_.cols << std::endl;

      // Set precision.
      out.precision( 16 );
      out << std::scientific;

      for (int j=0; j<denseESJacobian_.cols; j++)
      {
        for (int i=0; i<denseESJacobian_.rows; i++)
        {
          out << i+1 << " " << j+1 << " " ;
          out << denseESJacobian_.denseMtx( i, j ) << std::endl;
        }
      }
    }
  
    if (solver_ == "BLOCK_BASKER")
    {
      // Output sparse format.
      out << "coordinate real general" << std::endl;

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
              << std::scientific << Anewval_[ptr]
              << std::resetiosflags(std::ios_base::floatfield) << std::endl;
        }
      }
    }
  }

  // Close the file.
  out.close(); 
}

//---------------------------------------------------------------------------
// Function      : ESDirectSolver::printESResidual
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//---------------------------------------------------------------------------
void ESDirectSolver::printESResidual( const std::string& fileName )
{
  int numProcs = (builder_.getPDSComm())->numProc();
  int myProc = (builder_.getPDSComm())->procID();

  // Determine number of time-domain variables.
  MultiVector* B = lasProblem_.getRHS();
  EpetraVectorAccess* e_B = dynamic_cast<EpetraVectorAccess *>( B );
  int numVectors = B->numVectors();

  std::ofstream out;
  out.open( fileName.c_str() );

  // Write out banner.
  out << "%%MatrixMarket matrix array real general" << std::endl;

  // Write the dimensions of the sparse matrix: (# rows, #
  // columns, # matrix entries (counting duplicates as
  // separate entries)).
  out << n_*N_ << " " << numVectors << std::endl;

  // Set precision.
  out.precision( 16 );
  out << std::scientific;

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
        for (int i = 0; i < N_; i++)
        {
          int B_j_Row = i*n_ + nB;
          out << (*B_j)[B_j_Row] << std::endl;
        }
      }
    }
  }

  // Close the file.
  out.close();
}

//---------------------------------------------------------------------------
// Function      : ESDirectSolver::printESSolution
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//---------------------------------------------------------------------------
void ESDirectSolver::printESSolution( const std::string& fileName )
{
  int numProcs = (builder_.getPDSComm())->numProc();
  int myProc = (builder_.getPDSComm())->procID();

  // Determine number of time-domain variables.
  MultiVector* X = lasProblem_.getLHS();
  int numVectors = X->numVectors();

  std::ofstream out;
  out.open( fileName.c_str() );

  // Write out banner.
  out << "%%MatrixMarket matrix array real general" << std::endl;

  // Write the dimensions of the sparse matrix: (# rows, #
  // columns, # matrix entries (counting duplicates as
  // separate entries)).
  out << n_*N_ << " " << numVectors << std::endl;

  // Set precision.
  out.precision( 16 );
  out << std::scientific;

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
      for (int nB = 0; nB < n_; nB++)
      {
        // Copy values from B_j to B_ vector for solver.
        for (int i = 0; i < N_; i++)
        {
          int X_j_Row = i*n_ + nB;
          out << (*X_j)[X_j_Row] << std::endl;
        }
      }
    }
  }

  // Close the file.
  out.close();
}

} // namespace Linear
} // namespace Xyce

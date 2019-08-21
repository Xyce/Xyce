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

//-----------------------------------------------------------------------------
//
// Purpose        : PCE Direct Linear Solver Interface
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL 
//
// Creation Date  : 06/27/2019
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_PCEDirectSolver_h
#define Xyce_N_LAS_PCEDirectSolver_h

#include <string>

#include <N_LAS_fwd.h>
#include <N_UTL_fwd.h>
#include <N_LOA_fwd.h>

#include <N_LAS_Solver.h>
#include <N_LAS_TransformTool.h>
#include <N_LAS_HBBlockMatrixEntry.h>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_RCP.hpp>

#if defined(Xyce_AMESOS2) && !defined(SHYLUBASKER)
#include "Amesos2_Basker_TypeMap.hpp"
#include "Amesos2_Basker.hpp"
#endif

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : PCEDirectSolver
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL 
// Creation Date : 06/27/2019
//-----------------------------------------------------------------------------
class PCEDirectSolver : public Solver
{

public:
  // Constructor
  PCEDirectSolver(
    Builder &                   builder,
    Problem &                   problem,
    Util::OptionBlock &         options);

  // Destructor
  virtual ~PCEDirectSolver();

  // Set the solver options
  bool setOptions(const Util::OptionBlock & OB);
  bool setDefaultOptions();

  // Set individual options
  bool setParam( const Util::Param & param );

  // Get info such as Num Iterations, Residual, etc.
  bool getInfo( Util::Param & info );

  // Register the PCE builder
  void registerPCEBuilder( const Teuchos::RCP<PCEBuilder> & pceBuilderPtr )
    { pceBuilderPtr_ = pceBuilderPtr; }

  // Register the PCE loader
  void registerPCELoader( const Teuchos::RCP<Loader::PCELoader> & pceLoaderPtr )
    { pceLoaderPtr_ = pceLoaderPtr; }

  void setNumSamples(int numS)
    { numSamples_ = numS; }

  void setParameterOuterLoop (bool paramsOL)
    { paramsOuterLoop_ = paramsOL; }

  // Solve function: x = A^(-1) b.
  // input parameter 'ReuseFactors': If 'true', do not factor A, rather reuse
  // factors from previous solve.  Useful for inexact nonlinear techniques and
  // multiple RHS solves.
  int doSolve( bool reuse_factors, bool transpose = false );

private:

  // Initialize the nonzero entries of the block CRS structure with a value.
  void initializeBlockCRS( double val );

  // This function will allocate the space needed to store the PCE Jacobian.
  // For "LAPACK", this function will create a dense matrix of the entire PCE Jacobian size (very large).
  // For "BASKER", this function will use the time domain graph to create a structure for the block
  // Jacobian, filling out Ap_, Ai_, and Av_.
  void createBlockStructures();

  // This function will use the matrices from the PCE loader to form the PCE Jacobian.
  void formPCEJacobian();

  // Compute numeric factorization with updated PCE Jacobian.
  int numericFactorization();

  // Solve linear systems with direct factors.
  int solve();

  // Print methods.
  void printPCEJacobian( const std::string& fileName );
  void printPCEResidual( const std::string& fileName );
  void printPCESolution( const std::string& fileName );

  // Time-domain builder
  Builder & builder_;

  //Primary problem access
  Problem & lasProblem_;

  bool isInit_;

  // N_ is the number of samples
  // n_ is the original problem size
  int N_, n_;

  // How often the linear system should be written to file, if at all.
  int outputLS_;

  // Solver type.
  std::string solver_, solverDefault_;

  // Embedded Sampling loader.
  Teuchos::RCP<Loader::PCELoader> pceLoaderPtr_;

  // Embedded Sampling builder.
  Teuchos::RCP<PCEBuilder> pceBuilderPtr_;

  // Solution and RHS vectors
  Teuchos::SerialDenseMatrix<int,double> X_, B_, R_, A_;
  std::vector<Xyce::PCEBlockMatrixEntry> bX_, bB_;

  // Dense matrix for LAPACK implementation of direct solver.
  Xyce::PCEBlockMatrixEntry densePCEJacobian_;
  Teuchos::RCP< Teuchos::SerialDenseSolver<int,double> > lapackSolver_;

  // Block CCS matrix for Basker.
  std::vector<int> Acol_ptr_, Arow_idx_;
  std::vector<Xyce::PCEBlockMatrixEntry> Aval_;

  // <double> matrix for single-value Basker.
  std::vector<int> Anewcol_ptr_, Anewrow_idx_;
  std::vector<double> Anewval_;

#if defined(Xyce_AMESOS2) && !defined(SHYLUBASKER)
  Basker::Basker<int, double> basker_;
  Basker::Basker<int, Xyce::PCEBlockMatrixEntry> blockBasker_;
#endif

  // Serialized objects for parallel loading.
  Teuchos::RCP<Epetra_Map> serialMap_;
  Teuchos::RCP<Epetra_Import> serialImporter_;
  Teuchos::RCP<Epetra_Vector> serialX_, serialB_;

  //Options
  Util::OptionBlock * options_;

  //Timer
  Util::Timer * timer_;

  int numSamples_;
  bool paramsOuterLoop_;
};

} // namespace Linear
} // namespace Xyce

#endif // Xyce_N_LAS_PCEDirectSolver_h


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

//-----------------------------------------------------------------------------
//
// Purpose        : ES Direct Linear Solver Interface
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL 
//
// Creation Date  : 05/24/04
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_ESDirectSolver_h
#define Xyce_N_LAS_ESDirectSolver_h

#include <string>

#include <N_LAS_fwd.h>
#include <N_UTL_fwd.h>
#include <N_LOA_fwd.h>

#include <N_LAS_Solver.h>
#include <N_LAS_TransformTool.h>
#include <N_LAS_BlockMatrixEntry.h>

#include <N_TIA_DataStore.h>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_RCP.hpp>

#ifdef Xyce_AMESOS2_BASKER
#include "Amesos2_Basker_TypeMap.hpp"
#include "Amesos2_Basker.hpp"
#endif

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : ESDirectSolver
// Purpose       :
// Special Notes :
// Creator       : Heidi Thornquist, SNL 
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
class ESDirectSolver : public Solver
{

public:
  // Constructor
  ESDirectSolver(
    Builder &                   builder,
    Problem &                   problem,
    Util::OptionBlock &         options);

  // Destructor
  virtual ~ESDirectSolver();

  // Set the solver options
  bool setOptions(const Util::OptionBlock & OB);
  bool setDefaultOptions();

  // Set individual options
  bool setParam( const Util::Param & param );

  // Register the ES builder
  void registerESBuilder( const Teuchos::RCP<ESBuilder> & esBuilderPtr )
    { esBuilderPtr_ = esBuilderPtr; }

  // Register the ES loader
  void registerESLoader( const Teuchos::RCP<Loader::ESLoader> & esLoaderPtr )
    { esLoaderPtr_ = esLoaderPtr; }

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

  // This function will allocate the space needed to store the ES Jacobian.
  // For "LAPACK", this function will create a dense matrix of the entire ES Jacobian size (very large).
  // For "BASKER", this function will use the time domain graph to create a structure for the block
  // Jacobian, filling out Ap_, Ai_, and Av_.
  void createBlockStructures();

  // This function will use the matrices from the ES loader to form the ES Jacobian.
  void formESJacobian();

  // Compute numeric factorization with updated ES Jacobian.
  int numericFactorization();

  // Solve linear systems with direct factors.
  int solve();

  // Print methods.
  void printESJacobian( const std::string& fileName );
  void printESResidual( const std::string& fileName );
  void printESSolution( const std::string& fileName );

  // Time-domain builder
  Builder & builder_;

  bool isInit_;

  // N_ is the number of samples
  // n_ is the original problem size
  int N_, n_;

  // How often the linear system should be written to file, if at all.
  int outputLS_;

  // Solver type.
  std::string solver_, solverDefault_;

  // Embedded Sampling loader.
  Teuchos::RCP<Loader::ESLoader> esLoaderPtr_;

  // Embedded Sampling builder.
  Teuchos::RCP<ESBuilder> esBuilderPtr_;

  // Solution and RHS vectors 
  Teuchos::SerialDenseMatrix<int,double> X_, B_, R_, A_;
  std::vector<Xyce::ESBlockMatrixEntry> bX_, bB_;

  // Dense matrix for LAPACK implementation of direct solver.
  Xyce::ESBlockMatrixEntry denseESJacobian_;
  Teuchos::RCP< Teuchos::SerialDenseSolver<int,double> > lapackSolver_;

  // Block CCS matrix for Basker.
  std::vector<int> Acol_ptr_, Arow_idx_;
  std::vector<Xyce::ESBlockMatrixEntry> Aval_;

  // <double> matrix for single-value Basker.
  std::vector<int> Anewcol_ptr_, Anewrow_idx_;
  std::vector<double> Anewval_;

#ifdef Xyce_AMESOS2_BASKER

#ifdef Xyce_NEW_BASKER
  BaskerClassicNS::BaskerClassic<int, Xyce::ESBlockMatrixEntry > blockBasker_;
#else
  Basker::Basker<int, Xyce::ESBlockMatrixEntry> blockBasker_;
#endif

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

#endif // Xyce_N_LAS_ESDirectSolver_h


//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Purpose        : Amesos Direct Linear Solver Interface
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/24/04
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_IRSolver_h
#define Xyce_N_LAS_IRSolver_h

#include <string>

#include <N_LAS_fwd.h>
#include <N_UTL_fwd.h>

#include <N_LAS_Solver.h>
#include <N_LAS_TransformTool.h>
#include <Teuchos_RCP.hpp>

#ifdef Xyce_AMESOS2
#include <Amesos2.hpp>
#endif

class Amesos_BaseSolver;
class Epetra_LinearProblem;
class Epetra_CrsMatrix;
class Epetra_Export;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : IRSolver
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
class IRSolver : public Solver
{

public:
  // Constructor
  IRSolver(
    Problem &                   problem,
    Util::OptionBlock &         options);

  // Destructor
  ~IRSolver();

  // Set the solver options
  bool setOptions(const Util::OptionBlock & OB);
  bool setDefaultOptions();
  bool setNewtonIter( int nIter ) { nIter_ = nIter; return true; }

  // Set individual options
  bool setParam( const Util::Param & param );

  // Get info such as Num Iterations, Residual, etc.
  bool getInfo( Util::Param & info );

  // Solve function: x = A^(-1) b.
  // input parameter 'ReuseFactors': If 'true', do not factor A, rather reuse
  // factors from previous solve.  Useful for inexact nonlinear techniques and
  // multiple RHS solves.
  int doSolve( bool reuse_factors, bool transpose = false );

private:

  // Perform a standard solve without trying to reuse the solver
  int doStandardSolve( Epetra_LinearProblem * prob );

  //Solver Type
  std::string type_;

  //Solver tolerance.
  double ir_min_tol_, ir_tol_;

  //Solver defaults.
  static const std::string type_default_;
  static const double tol_default_;
  static const double min_tol_default_;

  //Primary problem access
  Epetra_LinearProblem * problem_;

  //Wrapped solver object
  Amesos_BaseSolver * asolver_;

#ifdef Xyce_AMESOS2
  Teuchos::RCP<Amesos2::Solver<Epetra_CrsMatrix,Epetra_MultiVector> > a2solver_;
#endif

  //Repivot every time or use static pivoting
  bool repivot_;
  
  //Output linear system every outputLS_ calls
  int outputLS_;
  int outputBaseLS_;
  int outputFailedLS_;

  // Newton Iteration
  int nIter_;

  // Transform Support
  Teuchos::RCP<Transform> transform_;
  Epetra_LinearProblem * tProblem_;

  //Options
  Util::OptionBlock * options_;

  //Timer
  Util::Timer * timer_;

};

} // namespace Linear
} // namespace Xyce

#endif // Xyce_N_LAS_IRSolver_h


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

//-----------------------------------------------------------------------------
//
// Purpose        : KSparse Direct Linear Solver Interface
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

#ifndef Xyce_N_LAS_KSparseSolver_h
#define Xyce_N_LAS_KSparseSolver_h

#include <string>

#include <N_UTL_fwd.h>
#include <N_LAS_fwd.h>

class Epetra_LinearProblem;
class Epetra_CrsMatrix;
class Epetra_Export;
class Epetra_CrsKundertSparse;
class Epetra_Import;
class Epetra_Map;
class Epetra_MultiVector;
class Amesos_BaseSolver;

// ----------   Xyce Includes   ----------

#include <N_LAS_Solver.h>
#include <Teuchos_RCP.hpp>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : KSparseSolver
// Purpose       : 
// Special Notes : 
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
class KSparseSolver : public Solver
{

public:
  // Constructor
  KSparseSolver(Problem & prob, Util::OptionBlock & options);

  // Destructor
  ~KSparseSolver();

  // Set the solver options
  bool setOptions(const Util::OptionBlock & OB);
   
  // Solve function: x = A^(-1) b. 
  // input parameter 'ReuseFactors': If 'true', do not factor A, rather reuse 
  // factors from previous solve.  Useful for inexact nonlinear techniques and
  // multiple RHS solves.  
  int doSolve( bool reuse_factors, bool transpose = false );
 
private:

  //Primary problem access
  Problem & lasProblem_;
  Epetra_LinearProblem * problem_;

  //Repivot every time or use static pivoting
  bool repivot_;

  //Output linear system every outputLS_ calls
  int outputLS_;
  int outputBaseLS_;
  int outputFailedLS_;

  // Transform Support
  Teuchos::RCP<Transform> transform_;
  Epetra_LinearProblem * tProblem_;

  // Serialized Matrix (if using parallel load serial solve scenario) 
  Teuchos::RCP<Epetra_Map> serialMap_;
  Teuchos::RCP<Epetra_LinearProblem> serialProblem_;
  Teuchos::RCP<Epetra_CrsMatrix> serialMat_;
  Teuchos::RCP<Epetra_MultiVector> serialLHS_, serialRHS_;
  Teuchos::RCP<Epetra_Import> serialImporter_;

  // Import and export matrices in parallel load serial solve scenario
  Epetra_LinearProblem * importToSerial();
  int exportToGlobal();

  // KSparse Solver
  Teuchos::RCP<Epetra_CrsKundertSparse> solver_;

  // KLU Solver (only for numerical failures)
  Teuchos::RCP<Amesos_BaseSolver> kluSolver_;

  //Options
  Util::OptionBlock * options_;

  //Timer
  Teuchos::RCP<Util::Timer> timer_;

};

} // namespace Linear
} // namespace Xyce

#endif // Xyce_N_LAS_KSparseSolver_h

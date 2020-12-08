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
// Purpose        : Amesos2 Direct Linear Solver Interface
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

#ifndef Xyce_N_LAS_Amesos2Solver_h
#define Xyce_N_LAS_Amesos2Solver_h

#include <string>

#include <N_LAS_fwd.h>
#include <N_UTL_fwd.h>

#include <N_LAS_Solver.h>
#include <N_LAS_TransformTool.h>
#include <Teuchos_RCP.hpp>

#include <Amesos2.hpp>

class Epetra_LinearProblem;
class Epetra_CrsMatrix;
class Epetra_Export;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : Amesos2Solver
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
class Amesos2Solver : public Solver
{

public:
  // Constructor
  Amesos2Solver(
    const std::string &         type,
    Problem &                   problem,
    Util::OptionBlock &         options);

  // Destructor
  ~Amesos2Solver();

  // Set the solver options
  bool setOptions(const Util::OptionBlock & OB);

  // Solve function: x = A^(-1) b.
  // input parameter 'ReuseFactors': If 'true', do not factor A, rather reuse
  // factors from previous solve.  Useful for inexact nonlinear techniques and
  // multiple RHS solves.
  int doSolve( bool reuse_factors, bool transpose = false );

private:

  //Solver Type
  const std::string type_;

  //Primary problem access
  Problem & lasProblem_;
  Epetra_LinearProblem * problem_;

  //Wrapped solver object
  Teuchos::RCP<Amesos2::Solver<Epetra_CrsMatrix,Epetra_MultiVector> > solver_;

  //Output linear system every outputLS_ calls
  int outputLS_;
  int outputBaseLS_;
  int outputFailedLS_;

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

#endif // Xyce_N_LAS_Amesos2Solver_h


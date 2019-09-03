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
// Purpose        : interface to linear problem
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/17/04
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_Problem_h
#define Xyce_N_LAS_Problem_h

#include <N_PDS_fwd.h>
#include <N_PDS_ParMap.h>
#include <N_LAS_MultiVector.h>

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;

class Epetra_LinearProblem;
class Epetra_Operator;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : Problem
// Purpose       : interface to linear problem
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/17/04
//-----------------------------------------------------------------------------
class Problem
{

public:
  //Constructors
  Problem( const RCP<Matrix> & A, const RCP<MultiVector> & x, const RCP<MultiVector> & b );
  Problem( Matrix* A, MultiVector* x, MultiVector* b );

  //Epetra constructors
  Problem( const RCP<Epetra_Operator> & Op, const RCP<MultiVector> & x, const RCP<MultiVector> & b );
  Problem( const RCP<Epetra_LinearProblem> & epetraProblem );
  Epetra_LinearProblem & epetraObj() { return *epetraProblem_; }

  //Destructor
  ~Problem();

  // Access solution and right-hand side vectors
  RCP<MultiVector>& getRHS() { return b_; }
  RCP<MultiVector>& getLHS() { return x_; }

  RCP<Matrix>& getJac () { return A_; }

  bool matrixFree() const { return(matrixFreeFlag_); }

private:

  RCP<Matrix> A_;
  RCP<Epetra_Operator> Op_;
  RCP<MultiVector> x_;
  RCP<MultiVector> b_;

  RCP<Epetra_LinearProblem> epetraProblem_;

  bool matrixFreeFlag_;
};

} // namespace Linear
} // namespace Xyce

#endif

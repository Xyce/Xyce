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
// Purpose        : 
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
#include <Epetra_CrsMatrix.h>

// ----------   Xyce Includes   ----------

#include <N_LAS_Problem.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : Problem::Problem
// Purpose       : 
// Special Notes : 
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
Problem::Problem( const RCP<Matrix> & A, const RCP<MultiVector> & x, const RCP<MultiVector> & b )
 : A_(A),
   x_(x),
   b_(b),
   epetraProblem_( rcp( new Epetra_LinearProblem( dynamic_cast<Epetra_RowMatrix*>(&(A_->epetraObj())),
                                             &(x_->epetraObj()),
                                             &(b_->epetraObj()) ) ) )
{
  matrixFreeFlag_ = false;
}

//-----------------------------------------------------------------------------
// Function      : Problem::Problem
// Purpose       : 
// Special Notes : 
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
Problem::Problem( Matrix* A, MultiVector* x, MultiVector* b )
 : A_(rcp(A,false)),
   x_(rcp(x,false)),
   b_(rcp(b,false)),
   epetraProblem_( rcp( new Epetra_LinearProblem( dynamic_cast<Epetra_RowMatrix*>(&(A_->epetraObj())),
                                             &(x_->epetraObj()),
                                             &(b_->epetraObj()) ) ) )
{
  matrixFreeFlag_ = false;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : Problem::Problem
// Purpose       : 
// Special Notes : 
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/04/08
//-----------------------------------------------------------------------------
Problem::Problem( const RCP<Epetra_Operator> & Op, const RCP<MultiVector> & x, const RCP<MultiVector> & b )
 : Op_(Op),
   x_(x),
   b_(b),
   epetraProblem_( rcp( new Epetra_LinearProblem( &*Op,
                                             &(x_->epetraObj()),
                                             &(b_->epetraObj()) ) ) )
{
  matrixFreeFlag_ = true;
}

//-----------------------------------------------------------------------------
// Function      : Problem::Problem
// Purpose       : 
// Special Notes : 
// Scope         : Public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 10/08/08
//-----------------------------------------------------------------------------
Problem::Problem( const RCP<Epetra_LinearProblem> & epetraProblem )
 : x_(Teuchos::rcp(new MultiVector(epetraProblem->GetLHS(), false))),
   b_(Teuchos::rcp(new MultiVector(epetraProblem->GetRHS(), false))),
   epetraProblem_(epetraProblem)
{
  if (epetraProblem_->GetMatrix()) {
    matrixFreeFlag_ = false;
    A_ = Teuchos::rcp(new Matrix(dynamic_cast<Epetra_CrsMatrix *>(epetraProblem_->GetMatrix()), false));
  }
  else {
    matrixFreeFlag_ = true;
    Op_ = Teuchos::rcp(epetraProblem_->GetOperator(), false);
  }
}

//-----------------------------------------------------------------------------
// Function      : Problem::~Problem
// Purpose       : 
// Special Notes : 
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
Problem::~Problem()
{
}

} // namespace Linear
} // namespace Xyce

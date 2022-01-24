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
// Purpose       : Implementation of the Epetra_Operator interface to define
//                 inv(A)*B, where inv(A) is computed by Amesos.
// Special Notes :
//
// Creator       : Heidi Thornquist, SNL
//
// Creation Date : 06/04/12
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <algorithm>

// ----------   Xyce Includes   ----------

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_MOROperators.h>

// ----------  Other Includes   ----------

#include "Epetra_MultiVector.h"

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : AmesosGenOp::AmesosGenOp
// Purpose       : Constructor 
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL 
// Creation Date : 06/04/12
//-----------------------------------------------------------------------------
AmesosGenOp::AmesosGenOp( const Teuchos::RCP<Amesos_BaseSolver>& solver,
                                      const Teuchos::RCP<Epetra_Operator>& B,
                                      bool useTranspose )
  : useTranspose_(useTranspose),
    solver_(solver),
    B_(B)
{
  problem_ = Teuchos::rcp( const_cast<Epetra_LinearProblem*>( solver->GetProblem() ), false);

  if ( solver_->UseTranspose() )
    solver_->SetUseTranspose(!useTranspose);
  else
    solver_->SetUseTranspose(useTranspose);

  if ( B_->UseTranspose() )
    B_->SetUseTranspose(!useTranspose);
  else
    B_->SetUseTranspose(useTranspose);
}

//-----------------------------------------------------------------------------
// Function      : AmesosGenOp::Apply()
// Purpose       : Applies the operator inv(A)*B*X = Y
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL 
// Creation Date : 06/04/12
//-----------------------------------------------------------------------------
int AmesosGenOp::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
{ 
  if (!useTranspose_) {

    // Storage for B*X 
    Epetra_MultiVector BX(X.Map(),X.NumVectors());

    // Apply B*X
    B_->Apply(X, BX);
    Y.PutScalar(0.0);

    // Set the LHS and RHS
    problem_->SetRHS(&BX);
    problem_->SetLHS(&Y);

    // Solve the linear system A*Y = BX
    solver_->Solve();
  }
  else {
    // Storage for A^{-T}*X
    Epetra_MultiVector ATX(X.Map(),X.NumVectors());
    Epetra_MultiVector tmpX = const_cast<Epetra_MultiVector&>(X);

    // Set the LHS and RHS
    problem_->SetRHS(&tmpX);
    problem_->SetLHS(&ATX);

    // Solve the linear system A^T*Y = X 
    solver_->Solve();

    // Apply B*ATX
    B_->Apply(ATX, Y);
  }
  
  return 0;
}

} // namespace Linear
} // namespace Xyce

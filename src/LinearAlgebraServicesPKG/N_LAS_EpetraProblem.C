//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
#include <Epetra_MultiVector.h>

// ----------   Xyce Includes   ----------

#include <N_LAS_EpetraProblem.h>

#include <N_LAS_Operator.h>
#include <N_LAS_EpetraMatrix.h>
#include <N_LAS_EpetraBlockMatrix.h>
#include <N_LAS_EpetraMultiVector.h>
#include <N_LAS_MatrixFreeEpetraOperator.h>

#include <Teuchos_RCP.hpp>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : EpetraProblem::EpetraProblem
// Purpose       : 
// Special Notes : 
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
EpetraProblem::EpetraProblem( Matrix* A, MultiVector* x, MultiVector* b )
 : Problem(A,x,b),
   isOwned_(false)
{
  EpetraMatrixAccess* e_A = dynamic_cast<EpetraMatrixAccess *>( A_ );
  EpetraVectorAccess* e_x = dynamic_cast<EpetraVectorAccess *>( x_ );
  EpetraVectorAccess* e_b = dynamic_cast<EpetraVectorAccess *>( b_ );
  if (e_A)
  {
    epetraProblem_ = Teuchos::rcp( new Epetra_LinearProblem(
                                   dynamic_cast<Epetra_RowMatrix*>(&(e_A->epetraObj())),
                                   &(e_x->epetraObj()),
                                   &(e_b->epetraObj()) ) );
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraProblem::EpetraProblem
// Purpose       : 
// Special Notes : 
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/04/08
//-----------------------------------------------------------------------------
EpetraProblem::EpetraProblem( Operator* Op, MultiVector* x, MultiVector* b )
 : Problem(Op,x,b),
   isOwned_(false)
{
  EpetraVectorAccess* e_x = dynamic_cast<EpetraVectorAccess *>( x_ );
  EpetraVectorAccess* e_b = dynamic_cast<EpetraVectorAccess *>( b_ );
  epetraOp_ = matrixFreeEpetraOperator( Teuchos::rcp( Op, false ), 
                                        Teuchos::rcp( x->pmap(), false ) );
  epetraProblem_ =  Teuchos::rcp( new Epetra_LinearProblem( &*epetraOp_,
                                             &(e_x->epetraObj()),
                                             &(e_b->epetraObj()) ) ); 
}

//-----------------------------------------------------------------------------
// Function      : EpetraProblem::EpetraProblem
// Purpose       : 
// Special Notes : 
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 10/08/08
//-----------------------------------------------------------------------------
EpetraProblem::EpetraProblem( const Teuchos::RCP<Epetra_LinearProblem> & epetraProblem )
 : Problem(),
   isOwned_(true),
   epetraProblem_(epetraProblem)
{
  x_ = new EpetraMultiVector(epetraProblem->GetLHS(), false);
  b_ = new EpetraMultiVector(epetraProblem->GetRHS(), false);

  if (epetraProblem_->GetMatrix()) {
    matrixFreeFlag_ = false;
    A_ = new EpetraMatrix(dynamic_cast<Epetra_CrsMatrix *>(epetraProblem_->GetMatrix()), false);
  }
  else {
    matrixFreeFlag_ = true;
    epetraOp_ = Teuchos::rcp(epetraProblem_->GetOperator(),false);
  }
}

//-----------------------------------------------------------------------------
// Function      : EpetraProblem::~EpetraProblem
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL 
// Creation Date : 12/16/20
//-----------------------------------------------------------------------------
EpetraProblem::~EpetraProblem()
{
  if (isOwned_)
  {
    delete A_;
    delete x_;
    delete b_;
  }
}

} // namespace Linear
} // namespace Xyce

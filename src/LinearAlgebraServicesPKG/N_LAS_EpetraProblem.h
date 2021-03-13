//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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

#ifndef Xyce_N_LAS_EpetraProblem_h
#define Xyce_N_LAS_EpetraProblem_h

#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_PDS_ParMap.h>

#include <N_LAS_Problem.h>

#include <Teuchos_RCP.hpp>

class Epetra_LinearProblem;
class Epetra_Operator;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : EpetraProblem
// Purpose       : interface to Epetra linear problem
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/17/04
//-----------------------------------------------------------------------------
class EpetraProblem : public Problem
{
public:

  //Constructors
  EpetraProblem( Operator* Op, MultiVector* x, MultiVector* b );
  EpetraProblem( Matrix* A, MultiVector* x, MultiVector* b );

  //Epetra constructors
  EpetraProblem( const Teuchos::RCP<Epetra_LinearProblem> & epetraProblem );

  //Destructor
  virtual ~EpetraProblem();

  Epetra_LinearProblem & epetraObj() { return *epetraProblem_; }

private:

  bool isOwned_;

  Teuchos::RCP<Epetra_LinearProblem> epetraProblem_;
  Teuchos::RCP<Epetra_Operator> epetraOp_;
};

} // namespace Linear
} // namespace Xyce

#endif

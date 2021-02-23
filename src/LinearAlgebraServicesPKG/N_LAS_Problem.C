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

// ----------   Xyce Includes   ----------

#include <N_LAS_Problem.h>
#include <N_LAS_Operator.h>
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
Problem::Problem( Matrix* A, MultiVector* x, MultiVector* b )
 : A_(A),
   x_(x),
   b_(b),
   matrixFreeFlag_(false)
{}

//-----------------------------------------------------------------------------
// Function      : Problem::Problem
// Purpose       : 
// Special Notes : 
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/04/08
//-----------------------------------------------------------------------------
Problem::Problem( Operator* Op, MultiVector* x, MultiVector* b )
 : Op_(Op),
   x_(x),
   b_(b),
   matrixFreeFlag_(true)
{}

//-----------------------------------------------------------------------------
// Function      : Problem::Problem
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/04/08
//-----------------------------------------------------------------------------
Problem::Problem( const Problem& prob ) 
 : A_(prob.A_),
   Op_(prob.Op_),
   x_(prob.x_),
   b_(prob.b_)
{}

} // namespace Linear
} // namespace Xyce

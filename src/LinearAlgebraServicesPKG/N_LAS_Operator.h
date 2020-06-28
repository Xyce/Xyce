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
// Purpose        : Abstract interface to linear solver type.
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL
//
// Creation Date  : 06/19/20
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_Operator_h
#define Xyce_N_LAS_Operator_h

#include <N_LAS_fwd.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_BLAS_types.hpp>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : Operator
// Purpose       : Abstract interface to an operator.
// Special Notes :
// Creator       : Heidi Thornquist, SNL
// Creation Date : 06/19/20
//-----------------------------------------------------------------------------
class Operator
{

public:

  //Apply method
  virtual int apply(const Linear::MultiVector& X, Linear::MultiVector& Y, 
                    Teuchos::ETransp mode = Teuchos::NO_TRANS) const = 0;

  //Whether or not the operator supports Teuchos::TRANS or Teuchos::CONJ_TRANS
  virtual bool hasTransposeApply() const { return false; }

};

} // namespace Linear
} // namespace Xyce

#endif

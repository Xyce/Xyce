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
// Purpose        : This is a container class for solver information.
//                  It may occasionally contain stuff that isn't strictly
//                  pertaining to the solver state, but that is its primary
//                  intention.
//
//                  In general, stuff that goes into this class should
//                  be stuff needed by more than one device instance type.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/02/21
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_ExpressionGroupWrapper_h
#define Xyce_N_DEV_ExpressionGroupWrapper_h

#include <Teuchos_RCP.hpp>
#include <expressionGroup.h>

namespace Xyce {
namespace Device {

struct expressionGroupWrapper
{
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> expressionGroup_; ///< required for setting up expressions
};

}
}

#endif


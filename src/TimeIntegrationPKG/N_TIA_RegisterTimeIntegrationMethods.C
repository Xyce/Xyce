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


//-----------------------------------------------------------------------------
//
// Purpose       : This file contains the functions which define the
//		             time integration methods classes.
//
// Special Notes :
//
// Creator       : Buddy Watts
//
// Creation Date : 6/1/00
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_TIA_WorkingIntegrationMethod.h>

#include <N_TIA_DataStore.h>
#include <N_TIA_Gear12.h>
#include <N_TIA_NoTimeIntegration.h>
#include <N_TIA_OneStep.h>
#include <N_TIA_StepErrorControl.h>

namespace Xyce {
namespace TimeIntg {

void
registerTimeIntegrationMethods()
{
  registerTimeIntegrationMethod<NoTimeIntegration>();
  registerTimeIntegrationMethod<Gear12>();
  registerTimeIntegrationMethod<OneStep>();
}

} // namespace TimeIntg
} // namespace Xyce

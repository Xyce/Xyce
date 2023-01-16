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

//-----------------------------------------------------------------------------
//
// Purpose        : ROL Transient analysis classes
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Heidi Thornquist, SNL
//
// Creation Date  : 01/24/08
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ANP_ROL_TRAN_Optimization.h>
#include <N_ANP_ROL.h>

namespace Xyce {
namespace Analysis {

#ifdef Xyce_ROL

bool ROL_TRAN::doInit()
{ 
  bool ret = Transient::doInit();

  return ret;
}

bool ROL_TRAN::doRun()
{
  bool ret = doInit() && doTranOP() && doLoopProcess();

  return ret;
}

bool ROL_TRAN::doProcessSuccessfulStep()
{ 
  bool ret = Transient::doProcessSuccessfulStep();

  return ret;
}

#endif

} // namespace Analysis
} // namespace Xyce

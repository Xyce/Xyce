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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_LAS_HBSolverFactory.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL
//
// Creation Date  : 08/12/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.35 $
//
// Revision Date  : $Date: 2016/03/16 17:25:05 $
//
// Current Owner  : $Author: hkthorn $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ----------   Xyce Includes   ----------

#include <N_LAS_HBSolverFactory.h>
#include <N_LAS_Problem.h>
#include <N_LAS_MultiVector.h>

#include <N_ERH_ErrorMgr.h>

#include <N_LAS_AztecOOSolver.h>
#include <N_LAS_BelosSolver.h>
#include <N_LAS_HBDirectSolver.h>

#include <N_IO_CmdParse.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_ExtendedString.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : HBSolverFactory::HBSolverFactory
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 10/01/07
//-----------------------------------------------------------------------------
HBSolverFactory::HBSolverFactory(
  Linear::Builder &             builder)
  : hbOsc_(false),
    builder_(builder)
{
}

//-----------------------------------------------------------------------------
// Function      : HBSolverFactory::create
// Purpose       :
// Special Notes : Static
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
Solver *
HBSolverFactory::create(
  Util::OptionBlock &   options,
  Problem &             problem,
  const IO::CmdParse &  command_line) const
{
  std::string type = "AZTECOO";

  Util::ParamList::const_iterator itPI = options.begin();
  Util::ParamList::const_iterator endPI = options.end();
  for( ; itPI != endPI; ++itPI )
  {
    if( itPI->uTag() == "TYPE" && itPI->usVal() != "DEFAULT" )
    {
      type = itPI->usVal();
    }
  }

  // If the linear problem is matrix free, make sure an iterative method is being used.
  if (problem.matrixFree())
  {
    if ((type != "AZTECOO") && (type != "BELOS") && (type != "DIRECT"))
    {
      std::string msg = "The linear solver option that was specified is not compatible with a matrix free analysis type, changing to AZTECOO";
      Report::UserWarning0() << msg;
      type = "AZTECOO";
    }
  }

  if( type == "AZTECOO" )
  {
    return new AztecOOSolver( problem, options);
  }
  else if( type == "BELOS" )
  {
    return new BelosSolver( problem, options);
  }
  else if( type == "DIRECT" )
  {
    HBDirectSolver* newSolver = new HBDirectSolver(builder_, problem, options);

    newSolver->registerHBLoader( hbLoaderPtr_ );
    newSolver->registerHBBuilder( hbBuilderPtr_ );
    newSolver->setHBFreqs( freqs_ );
    newSolver->setFastTimes( times_ );
    newSolver->setHBOsc( hbOsc_ );

    return newSolver;
  }
  else
  {
    return new AztecOOSolver( problem, options);
  }

  return 0;
}

} // namespace Linear
} // namespace Xyce

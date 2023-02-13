//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Filename       : $RCSfile: N_LAS_ESSolverFactory.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL
//
// Creation Date  : 06/08/2018
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

#include <N_LAS_ESSolverFactory.h>
#include <N_LAS_Problem.h>
#include <N_LAS_MultiVector.h>

#include <N_ERH_ErrorMgr.h>

#include <N_LAS_SimpleSolver.h>
#include <N_LAS_AmesosSolver.h>
#include <N_LAS_AztecOOSolver.h>
#include <N_LAS_KSparseSolver.h>
#include <N_LAS_BelosSolver.h>
#ifdef Xyce_SHYLU
#include <N_LAS_ShyLUSolver.h>
#endif
#ifdef Xyce_AMESOS2
#include <N_LAS_Amesos2Solver.h>
#endif

#include <N_LAS_ESDirectSolver.h>

#include <N_IO_CmdParse.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_ExtendedString.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : ESSolverFactory::ESSolverFactory
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 10/01/07
//-----------------------------------------------------------------------------
ESSolverFactory::ESSolverFactory(
  Linear::Builder &             builder)
  : builder_(builder),
  numSamples_(1),
  paramsOuterLoop_(true)
{
}

//-----------------------------------------------------------------------------
// Function      : ESSolverFactory::create
// Purpose       :
// Special Notes : Static
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
Solver *
ESSolverFactory::create(
  Util::OptionBlock &   options,
  Problem &             problem,
  const IO::CmdParse &  command_line) const
{

  int lsDim = (problem.getRHS())->globalLength();

  //If the linear system is trivial, i.e. the matrix is 1x1, then just create a simple solver
  if (lsDim == 1)
  {
    return new SimpleSolver(problem, options);
  }

#ifdef Xyce_PARALLEL_MPI
  std::string type = "AZTECOO";
  if (!problem.matrixFree() && lsDim < Xyce::Linear::iterativeMin)
  {
    type = "KLU";
  }
#else
  std::string type = "KLU";
  if (problem.matrixFree())
  {
    type = "AZTECOO";
  }
#endif  
 
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
    if ((type != "AZTECOO") && (type != "BELOS"))
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
  else if( type == "KSPARSE" )
  {
    return new KSparseSolver( problem, options);
  }
#ifdef Xyce_SHYLU
  else if( type == "SHYLU" )
  {
    return new ShyLUSolver( problem, options);
  }
#endif
#ifdef Xyce_AMESOS2
  else if( type == "SHYLU_BASKER" || type == "BASKER" || type == "KLU2" ) 
  {
    return new Amesos2Solver( type, problem, options );
  }
#ifdef Xyce_AMESOS2_BASKER
  else if( type == "BLOCK_BASKER" || type == "LAPACK" ) 
  {
    ESDirectSolver* newSolver = new ESDirectSolver(builder_, problem, options);

    newSolver->registerESLoader( esLoaderPtr_ );
    newSolver->registerESBuilder( esBuilderPtr_ );
    newSolver->setNumSamples( numSamples_ );
    newSolver->setParameterOuterLoop (paramsOuterLoop_);

    return newSolver;
  } 
#endif
#endif
  else
  {
    return new AmesosSolver( type, problem, options);
  }

  return 0;
}

} // namespace Linear
} // namespace Xyce

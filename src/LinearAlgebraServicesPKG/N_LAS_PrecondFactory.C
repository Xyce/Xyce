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
//
// Purpose        : 
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 10/01/07
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_LAS_PrecondFactory.h>
#include <N_LAS_IfpackPrecond.h>
#include <N_LAS_NoPrecond.h>

#include <N_LAS_Problem.h>

#include <N_ERH_ErrorMgr.h>

#include <N_UTL_OptionBlock.h>

// ---------- Trilinos Includes ----------

//-----------------------------------------------------------------------------
// Function      : N_LAS_PrecondFactory::create
// Purpose       :
// Special Notes : Static
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 10/01/07
//-----------------------------------------------------------------------------
Teuchos::RCP<N_LAS_Preconditioner> N_LAS_PrecondFactory::create( Util::OptionBlock & options,
								 const Teuchos::RCP<N_LAS_Problem> & problem)
{
  string type = "IFPACK";
  Teuchos::RCP<N_LAS_Preconditioner> precond;

  Util::OptionBlock::ParamIter itPI = options.begin();
  Util::OptionBlock::ParamIter endPI = options.end();
  for( ; itPI != endPI; ++itPI )
  {
    if( itPI->uTag() == "PREC_TYPE" && itPI->usVal() != "DEFAULT" )
      type = itPI->usVal();
  }
  if (problem->matrixFree())
  {
    // matrix free!   (duh!)
    type = "NONE";
  }
   
  if (type == "IFPACK") 
  {
    precond = Teuchos::rcp( new N_LAS_IfpackPrecond() );

    // First set the options for the preconditioner.
    precond->setOptions( options );
 
    // Now the preconditioner can create the graph for the preconditioner
    // Note:  The fill and overlap must be set before this point.
    precond->initGraph( problem );
  }
  else if (type == "NONE") {
    // Create an empty preconditioner, which does nothing.
    precond = Teuchos::rcp( new N_LAS_NoPrecond() );
  }
  else
  {
    Xyce::Report::DevelFatal0().in("N_LAS_PrecondFactory::create()") << "preconditioning type " << type << " unrecognized!";
  }                                                                                    
  return precond;
}

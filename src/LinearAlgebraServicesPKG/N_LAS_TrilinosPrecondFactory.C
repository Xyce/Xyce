//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose        : Implementation of Trilinos Preconditioning Factory
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

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <Xyce_config.h>

#include <N_LAS_TrilinosPrecondFactory.h>
#include <N_LAS_IfpackPrecond.h>
#include <N_LAS_NoPrecond.h>

#include <N_LAS_Problem.h>

#include <N_ERH_ErrorMgr.h>

#include <N_UTL_OptionBlock.h>

// ---------- Trilinos Includes ----------

#include <Epetra_Operator.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : TrilinosPrecondFactory::TrilinosPrecondFactory
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 10/01/07
//-----------------------------------------------------------------------------
TrilinosPrecondFactory::TrilinosPrecondFactory( const Util::OptionBlock & OB )
{
  OB_ = Teuchos::rcp( &OB, false );
  precType_ = "IFPACK";

  Util::ParamList::const_iterator itPI = OB.begin();
  Util::ParamList::const_iterator endPI = OB.end();

  for( ; itPI != endPI; ++itPI )
  {
    if( itPI->uTag() == "PREC_TYPE" && itPI->usVal() != "DEFAULT" )
      precType_ = itPI->usVal();
  }
}
//-----------------------------------------------------------------------------
// Function      : TrilinosPrecondFactory::create
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 10/01/07
//-----------------------------------------------------------------------------
Teuchos::RCP<Preconditioner>
TrilinosPrecondFactory::create( const Teuchos::RCP<Problem> & problem ) const
{
  Teuchos::RCP<Preconditioner> precond;

  std::string prec_type = precType_;

  // If the problem is matrix free that we can't use any of the available preconditioners.
  if (problem->matrixFree())
  {
    prec_type = "NONE";
  }

  if (prec_type == "IFPACK")
  {
    precond = Teuchos::rcp( new IfpackPrecond() );

    // First set the options for the preconditioner.
    precond->setOptions( *OB_ );

    // Now the preconditioner can create the graph for the preconditioner
    // Note:  The fill and overlap must be set before this point.
    precond->initGraph( problem );
  }
  else if (prec_type == "NONE") {
    // Create an empty preconditioner, which does nothing.
    precond = Teuchos::rcp( new NoPrecond() );
  }
  else
  {
    Xyce::Report::DevelFatal0().in("TrilinosPrecondFactory::create()") << "preconditioning type " << prec_type << " unrecognized!";
  }
  return precond;
}

} // namespace Linear
} // namespace Xyce

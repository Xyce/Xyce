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

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_LAS_Preconditioner.h>
#include <N_LAS_HBPrecondFactory.h>
#include <N_LAS_HBBlockJacobiPrecond.h>

#include <N_LAS_Problem.h>
#include <N_LAS_System.h>
#include <N_LAS_NoPrecond.h>

#include <N_ERH_ErrorMgr.h>

// ---------- Trilinos Includes ----------

#include <Epetra_Operator.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : HBPrecondFactory::HBPrecondFactory
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 10/01/07
//-----------------------------------------------------------------------------
HBPrecondFactory::HBPrecondFactory(
  const Util::OptionBlock &     option_block,
  Linear::Builder &             builder)
  : hbOsc_(false),
    precType_("BLOCK_JACOBI"),
    builder_(builder)
{
  optionBlock_ = Teuchos::rcp( new Util::OptionBlock( option_block ) );

  Util::ParamList::const_iterator it_tpL = option_block.begin();
  Util::ParamList::const_iterator end_tpL = option_block.end();
  for( ; it_tpL != end_tpL; ++it_tpL )
  {
    if( it_tpL->uTag() == "PREC_TYPE" && it_tpL->usVal() != "DEFAULT" )
      precType_ = it_tpL->usVal();
  }
}

//-----------------------------------------------------------------------------
// Function      : HBPrecondFactory::create
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 10/01/07
//-----------------------------------------------------------------------------
Teuchos::RCP<Preconditioner>
HBPrecondFactory::create( const Teuchos::RCP<System> & lasSysPtr ) const
{
  Teuchos::RCP<Preconditioner> precond;

  if (precType_ == "NONE") {
    // Create an empty preconditioner, which does nothing.
    precond = Teuchos::rcp( new NoPrecond() );
  }
  else if (precType_ == "BLOCK_JACOBI") {
    precond = Teuchos::rcp( new HBBlockJacobiPrecond(builder_) );
    precond->setOptions( *optionBlock_ );

    // Register necessary classes for block Jacobi preconditioner.
    Teuchos::RCP<HBBlockJacobiPrecond> tmpPrecond = Teuchos::rcp_dynamic_cast<HBBlockJacobiPrecond>( precond );

    tmpPrecond->registerLinearSystem( lasSysPtr );
    tmpPrecond->registerHBLoader( hbLoaderPtr_ );
    tmpPrecond->registerHBBuilder( hbBuilderPtr_ );
    tmpPrecond->setHBFreqs( freqs_ );
    tmpPrecond->setFastTimes( times_ );
    tmpPrecond->setHBOsc( hbOsc_ );

    // Initialize the graph for the preconditioner
    Teuchos::RCP<Problem> tmpProblem;
    tmpPrecond->initGraph( tmpProblem );
  }
  else {
    Xyce::Report::DevelFatal0().in("HBPrecondFactory::create()") << "preconditioning type " << precType_ << " unrecognized!";
  }

  return precond;
}

} // namespace Linear
} // namespace Xyce

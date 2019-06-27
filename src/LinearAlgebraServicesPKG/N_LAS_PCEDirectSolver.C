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
// Purpose        : PCE direct solver wrapper
// Special Notes  :
// Creator        : Eric Keiter, SNL
// Creation Date  : 6/27/2019
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Xyce Includes ----------

#include <N_UTL_fwd.h>

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_PCEDirectSolver.h>
#include <N_LAS_PCEBuilder.h>
#include <N_LOA_PCELoader.h>
#include <N_LAS_Builder.h>
#include <N_LAS_Solver.h>
#include <N_LAS_Problem.h>
#include <N_LAS_Vector.h>
#include <N_LAS_FilteredMatrix.h>
#include <N_LAS_TransformTool.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Timer.h>
#include <N_UTL_AssemblyTypes.h>
#include <N_PDS_Comm.h>

#include <Epetra_MultiVector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Import.h>
#include <Epetra_Util.h>
#include <Epetra_Comm.h>

#include <Teuchos_Utils.hpp>
#include <Teuchos_BLAS.hpp>

#include <set>
#include <utility>
#include <numeric>
#include <iostream>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : PCEDirectSolver::PCEDirectSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//-----------------------------------------------------------------------------
PCEDirectSolver::PCEDirectSolver(
  Builder &       builder,
  Problem &       problem,
  Util::OptionBlock &   options)
  : Solver(false),
    builder_(builder),
    lasProblem_(problem),
    isInit_(false),
    N_(0),
    n_(0),
    M_(0),
    numAugRows_(0),
    outputLS_(0),
    solver_(""),
    solverDefault_("LAPACK"),
    options_( new Util::OptionBlock( options ) ),
    timer_( new Util::Timer() )
{
  setDefaultOptions();

  setOptions( options );
}

//-----------------------------------------------------------------------------
// Function      : PCEDirectSolver::~PCEDirectSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//-----------------------------------------------------------------------------
PCEDirectSolver::~PCEDirectSolver()
{
  delete timer_;
  delete options_;
}

//-----------------------------------------------------------------------------
// Function      : PCEDirectSolver::setOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//-----------------------------------------------------------------------------
bool PCEDirectSolver::setOptions( const Util::OptionBlock & OB )
{
  Util::ParamList::const_iterator it_tpL = OB.begin();
  Util::ParamList::const_iterator end_tpL = OB.end();
  for (; it_tpL != end_tpL; ++it_tpL)
  {
    setParam( *it_tpL );
  }

  if ( solver_ == "DEFAULT" )
  {
    solver_ = solverDefault_;
  }

#if defined(Xyce_AMESOS2) && !defined(SHYLUBASKER)
  if ( solver_ != "LAPACK" && solver_ != "BASKER" && solver_ != "BLOCK_BASKER" )
#else
  if ( solver_ != "LAPACK" )
#endif
  {
    Report::UserWarning0()
        << "PCEDirectSolver does not recognize solver type " << solver_ << " setting to LAPACK";
    solver_ = "LAPACK";
  }

  if( options_ ) delete options_;
  options_ = new Util::OptionBlock( OB );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PCEDirectSolver::setDefaultOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//-----------------------------------------------------------------------------
bool PCEDirectSolver::setDefaultOptions()
{
  solver_ = solverDefault_;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PCEDirectSolver::setParam
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//-----------------------------------------------------------------------------
bool PCEDirectSolver::setParam( const Util::Param & param )
{
  std::string tag = param.tag();
  std::string uTag = param.uTag();

  if( uTag == "DIRECT_SOLVER" )
    solver_ = param.usVal();

  if( uTag == "OUTPUT_LS" ) 
    outputLS_ = param.getImmutableValue<int>();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PCEDirectSolver::getInfo
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//-----------------------------------------------------------------------------
bool PCEDirectSolver::getInfo( Util::Param & info )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : PCEDirectSolver::doSolve
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//-----------------------------------------------------------------------------
int PCEDirectSolver::doSolve( bool reuse_factors, bool transpose )
{
  // Start the timer...
  timer_->resetStartTime();
  int linearStatus = 0;

  return 0;
}

//---------------------------------------------------------------------------
// Function      : PCEDirectSolver::createBlockStructures
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//---------------------------------------------------------------------------
void PCEDirectSolver::createBlockStructures()
{
  int numProcs = (builder_.getPDSComm())->numProc();
  int myProc = (builder_.getPDSComm())->procID();

}

//---------------------------------------------------------------------------
// Function      : PCEDirectSolver::initializeBlockCRS
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//---------------------------------------------------------------------------
void PCEDirectSolver::initializeBlockCRS( double val )
{
  // Initialize the dense or diagonal blocks to the input value.
  for (unsigned int i=0; i < Aval_.size(); i++)
  {
    Aval_[i].putScalar( val );
  }
}

//---------------------------------------------------------------------------
// Function      : PCEDirectSolver::formPCEJacobian
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//---------------------------------------------------------------------------
void PCEDirectSolver::formPCEJacobian()
{
}

//---------------------------------------------------------------------------
// Function      : PCEDirectSolver::numericFactorization
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//---------------------------------------------------------------------------
int PCEDirectSolver::numericFactorization()
{
  int linearStatus = 0;
  return linearStatus;
}

//---------------------------------------------------------------------------
// Function      : PCEDirectSolver::solve
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//---------------------------------------------------------------------------
int PCEDirectSolver::solve()
{
  int linearStatus = 0;
  return linearStatus;
}

//---------------------------------------------------------------------------
// Function      : PCEDirectSolver::printPCEJacobian
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//---------------------------------------------------------------------------
void PCEDirectSolver::printPCEJacobian( const std::string& fileName )
{
}

//---------------------------------------------------------------------------
// Function      : PCEDirectSolver::printPCEResidual
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//---------------------------------------------------------------------------
void PCEDirectSolver::printPCEResidual( const std::string& fileName )
{
}

//---------------------------------------------------------------------------
// Function      : PCEDirectSolver::printPCESolution
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2018
//---------------------------------------------------------------------------
void PCEDirectSolver::printPCESolution( const std::string& fileName )
{
}

} // namespace Linear
} // namespace Xyce

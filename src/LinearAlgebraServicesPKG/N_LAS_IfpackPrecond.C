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
// Purpose        : Implementation file for the Iterative linear solver
//                  interface.
//
// Special Notes  :
//
// Creator        : Heidi K. Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 09/27/07
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <sstream>

// ----------   Xyce Includes   ----------
#include <N_ERH_ErrorMgr.h>
#include <N_LAS_IfpackPrecond.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Problem.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Timer.h>

#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>
#include <Ifpack_IlukGraph.h>
#include <Ifpack_CrsRiluk.h>
#include <Ifpack.h>
#include <Ifpack_Preconditioner.h>

namespace Xyce {
namespace Linear {

// static class member initializations
// Default preconditioner values
const bool   IfpackPrecond::useFactory_default_ = false;
const std::string IfpackPrecond::ifpackType_default_ = "Amesos";
const double IfpackPrecond::diagPerturb_default_= 0.0;
const int    IfpackPrecond::overlap_default_    = 0;
const double IfpackPrecond::dropTol_default_    = 1.0e-03;
const double IfpackPrecond::ilutFill_default_   = 2.0;
const double IfpackPrecond::rThresh_default_    = 1.0001;
const double IfpackPrecond::aThresh_default_    = 0.0001;

using Xyce::VERBOSE_LINEAR;

//-----------------------------------------------------------------------------
// Function      : IfpackPrecond::IfpackPrecond
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
IfpackPrecond::IfpackPrecond()
  : Preconditioner()
{
  setDefaultOptions();
}

//-----------------------------------------------------------------------------
// Function      : IfpackPrecond::setDefaultOptions
// Purpose       : resets Ifpack options
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool IfpackPrecond::setDefaultOptions()
{
  // Set defaults
  useFactory_ = useFactory_default_;
  ifpackType_ = ifpackType_default_;
  diagPerturb_ = diagPerturb_default_;
  dropTol_ = dropTol_default_;
  ilutFill_ = ilutFill_default_;
  rThresh_ = rThresh_default_;
  aThresh_ = aThresh_default_;
  overlap_ = overlap_default_;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : IfpackPrecond::setOptions
// Purpose       : sets Ifpack options and params from modelblock
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool IfpackPrecond::setOptions( const Util::OptionBlock & OB )
{
  // Set the parameters from the list
  Util::ParamList::const_iterator it_tpL = OB.begin();
  Util::ParamList::const_iterator end_tpL = OB.end();
  for (; it_tpL != end_tpL; ++it_tpL)
    {
      setParam( *it_tpL );
    }

  // store for restart of solver_
  if( &OB != options_.get() )
    {
      options_ = Teuchos::rcp( new Util::OptionBlock(OB) );
    }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : IfpackPrecond::setParam
// Purpose       : sets Ifpack option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool IfpackPrecond::setParam( const Util::Param & param )
{
  std::string tag = param.tag();
  std::string uTag = param.uTag();

  // Set our copies of these parameters that get passed to the solver in the
  // "iterate" command
  if( tag == "AZ_overlap" )
    overlap_ = param.getImmutableValue<int>();
  else if( tag == "AZ_athresh")
    aThresh_ = param.getImmutableValue<double>();
  else if( tag == "AZ_rthresh")
    rThresh_ = param.getImmutableValue<double>();
  else if( tag == "AZ_drop")
    dropTol_ = param.getImmutableValue<double>();
  else if( tag == "AZ_ilut_fill")
    ilutFill_ = param.getImmutableValue<double>();
  else if( tag == "use_ifpack_factory")
    useFactory_ = param.getImmutableValue<int>();
  else if( tag == "ifpack_type")
    ifpackType_ = param.usVal();
  else if( tag == "diag_perturb")
    diagPerturb_ = param.getImmutableValue<double>();
  else
    return false;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : IfpackPrecond::initGraph
// Purpose       : Set up the graph pattern for the preconditioner.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool IfpackPrecond::initGraph( const Teuchos::RCP<Problem> & problem )
{
  bool precStatus = true;

  if (useFactory_)
  {
    // Because a new graph is being used, recreate the preconditioner.

    // Create Ifpack factory.
    Ifpack Factory;

    // Create the preconditioner.
    Epetra_CrsMatrix * epetraA = dynamic_cast<Epetra_CrsMatrix*>(problem->epetraObj().GetMatrix());
    ifpackPrecond_ = Teuchos::rcp( Factory.Create(ifpackType_, epetraA, overlap_) );

    if (ifpackPrecond_ == Teuchos::null) {
      Xyce::Report::DevelFatal0().in("IfpackPrecond::initGraph()") << "preconditioning type " << ifpackType_ << " unrecognized!";
    }

    // Pass solver type to Amesos.
    Teuchos::ParameterList List;
    if (ifpackType_ == "Amesos") 
    {
      List.set("amesos: solver type", "Amesos_Klu");
      if (diagPerturb_ != 0.0)
        List.set("AddToDiag", diagPerturb_);
    }
    else if (ifpackType_ == "ILUT")
    {
      List.set("fact: absolute threshold", aThresh_);
      List.set("fact: relative threshold", rThresh_);
      List.set("fact: ilut level-of-fill", (double)ilutFill_);
      List.set("fact: drop tolerance", dropTol_);
    }
    else if (ifpackType_ == "ILU")
    {
      List.set("fact: absolute threshold", aThresh_);
      List.set("fact: relative threshold", rThresh_);
      List.set("fact: level-of-fill", (int)ilutFill_);
      List.set("fact: drop tolerance", dropTol_);
    }

    // Set the parameters
    IFPACK_CHK_ERR(ifpackPrecond_->SetParameters(List));

    // Compute symbolic factorization stage of preconditioner.
    IFPACK_CHK_ERR(ifpackPrecond_->Initialize());

    // Set the Epetra object to point to this preconditioner.
    epetraPrec_ = ifpackPrecond_;
  }
  else
  {
    // Create the graph.
    Epetra_CrsMatrix * epetraA = dynamic_cast<Epetra_CrsMatrix*>(problem->epetraObj().GetMatrix());
    const Epetra_CrsGraph & Graph = epetraA->Graph();
    ilukGraph_ = Teuchos::rcp( new Ifpack_IlukGraph( Graph, static_cast<int>(ilutFill_), overlap_ ) );
    int graphRet = ilukGraph_->ConstructFilledGraph();
    if ( graphRet != 0 )
      precStatus = false;

    // Because a new graph was created, destroy the preconditioner to make sure it's recreated.
    rILUK_ = Teuchos::null;
    epetraPrec_ = Teuchos::null;
  }

  return precStatus;
}

//-----------------------------------------------------------------------------
// Function      : IfpackPrecond::initValues
// Purpose       : Set the values for the preconditioner.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool IfpackPrecond::initValues( const Teuchos::RCP<Problem> & problem )
{
  bool precStatus = true;
  problem_ = problem;

  if (useFactory_)
  {
    // It is currently assumed that the same Matrix object you are passing in
    // is the one you used to set the graph.  This method takes a pointer to
    // it and will use it during the preconditioner computation.
    if ( Teuchos::is_null( ifpackPrecond_ ) )
      initGraph( problem_ );
  }
  else
  {
    // Initialize the graph if we need to.
    if ( Teuchos::is_null( ilukGraph_ ) )
      initGraph( problem_ );

    // Create the preconditioner if one doesn't exist.
    if ( Teuchos::is_null( rILUK_ ) )
    {
      rILUK_ = Teuchos::rcp( new Ifpack_CrsRiluk( *ilukGraph_ ) );
      rILUK_->SetAbsoluteThreshold( aThresh_ );
      rILUK_->SetRelativeThreshold( rThresh_ );
    }

    Epetra_CrsMatrix * epetraA = dynamic_cast<Epetra_CrsMatrix*>(problem_->epetraObj().GetMatrix());
    int initErr = rILUK_->InitValues( *epetraA );
    if (initErr < 0)
      precStatus = false;

    epetraPrec_ = rILUK_;
  }

  return precStatus;
}

//-----------------------------------------------------------------------------
// Function      : IfpackPrecond::compute
// Purpose       : Compute a preconditioner M such that M ~= A^{-1}.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool IfpackPrecond::compute()
{
  bool precStatus = true;
  if ( Teuchos::is_null( epetraPrec_ ) )
    return false;

  if (useFactory_)
  {
    // Build the preconditioner using the values in problem_; numeric factorization.
    // This may fail if "Amesos" used to compute a exact factorization of each subdomain and
    // any of those subdomains are numerically singular.
    int ierr = ifpackPrecond_->Compute();
    if (ierr < 0 && ifpackType_ == "Amesos")
    {
      Report::UserWarning0()
        << "IfpackPrecond::compute():  Subdomain solve failed, numerically singular, changing to ILU.";
      ifpackPrecond_ = Teuchos::null;
      ifpackType_ = "ILU";
      initValues( problem_ );
      ierr = ifpackPrecond_->Compute();
    }
    IFPACK_CHK_ERR(ierr)
  }
  else
  {
    Epetra_CrsMatrix * epetraA = dynamic_cast<Epetra_CrsMatrix*>(problem_->epetraObj().GetMatrix());
    bool transpose = epetraA->UseTranspose();

    int factErr = rILUK_->Factor();
    if (factErr < 0)
      return false;

    double condest = 0.0;
    if (VERBOSE_LINEAR)
    {
      rILUK_->Condest( transpose, condest );
    }

    // Define label for printing out during the solve phase
    std::ostringstream ost;
    ost << "Ifpack_CrsRiluk Preconditioner: LevelFill = " << ilutFill_ << std::endl
        << "                                Overlap = " << overlap_ << std::endl
        << "                                Athresh = " << aThresh_ << std::endl
        << "                                Rthresh = " << rThresh_ << std::endl;

    if (VERBOSE_LINEAR)
      ost << "                                CondEst = " << condest  << std::endl;

    ost << "                                ErrCode = " << factErr  << std::endl;

    std::string label = ost.str();
    rILUK_->SetLabel(label.c_str());

    rILUK_->SetUseTranspose( transpose );
  }

  return precStatus;
}

//-----------------------------------------------------------------------------
// Function      : IfpackPrecond::apply
// Purpose       : Calls the actual preconditioner to apply y = M*x.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
int IfpackPrecond::apply( MultiVector & x, MultiVector & y )
{
  int precStatus = 0;

  // If there is no preconditioner to apply return a nonzero code
  if( Teuchos::is_null(epetraPrec_) )
    precStatus = -1;
  else
    precStatus = epetraPrec_->ApplyInverse( x.epetraObj(), y.epetraObj() );

  return precStatus;
}

} // namespace Linear
} // namespace Xyce

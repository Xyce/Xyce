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
// Purpose        : Interface to Xyce vectors for NOX.
//
// Special Notes  :
//
// Creator        : Tammy Kolda, NLS, 8950
//
// Creation Date  : 01/31/02
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <list>

// ----------   Xyce Includes   ----------
#include <N_NLS_fwd.h>
#include <N_NLS_NOX_Interface.h>
#include <N_NLS_NOX_XyceTests.h>
#include <N_ERH_ErrorMgr.h>
#include <N_ERH_Message.h>
#include <N_UTL_Param.h>
#include <N_PDS_fwd.h>
#include <N_UTL_OptionBlock.h>
#include <LOCA.H>
#include <N_NLS_ReturnCodes.h>
#include <N_LOA_Loader.h>

// The following are needed to build AugmentLinSys strategies
#include <N_NLS_NOX_AugmentLinSys.h>
#include <N_NLS_NOX_AugmentLinSys_PseudoTransient.h>
#include <N_NLS_NOX_AugmentLinSys_GStepping.h>
#include <N_NLS_NOX_AugmentLinSys_IC.h>
#include <N_NLS_NOX_AugmentLinSys_IC_Gmin.h>
#include <N_LAS_Builder.h>
#include <N_LAS_System.h>

#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_MachDepParams.h>

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

//-----------------------------------------------------------------------------
// Function      : ParameterSet::ParameterSet
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
ParameterSet::ParameterSet(Nonlinear::AnalysisMode mode) :
  allParams_(Teuchos::rcp(new Teuchos::ParameterList)),
  noxParams_(allParams_->sublist("NOX")),
  locaParams_(allParams_->sublist("LOCA")),
  debugParams_(allParams_->sublist("DEBUG")),
  comboPtr_(Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR))),
  isParamsSet_(false),
  isStatusTestsSet_(false),
  continuationSpecified_(false),
  mode_(mode),
  noxSolver(0),
  voltageListType_(VLT_None),
  voltageScaleFactor_(1.0),
  gstepping_minimum_conductance_(0.0),
  savedLocaOptions_(false),
  debugLevel_(0),
  debugMinTimeStep_(0),
  debugMaxTimeStep_(Util::MachineDependentParams::IntMax()),
  debugMinTime_(0.0),
  debugMaxTime_(Util::MachineDependentParams::DoubleMax()),
  screenOutputFlag_(false),
  maskingFlag_(false)
{
  // Add the main status test to the list of tests to delete in dtor
  tests_.push_back(comboPtr_);

  // Default printing options
#ifdef Xyce_VERBOSE_NOX
  noxParams_.sublist("Printing")
    .set("Output Information",
		  NOX::Utils::Error +
		  NOX::Utils::Warning +
		  NOX::Utils::OuterIteration +
		  NOX::Utils::OuterIterationStatusTest +
		  NOX::Utils::InnerIteration +
                  NOX::Utils::Details +
		  NOX::Utils::StepperIteration +
		  NOX::Utils::StepperDetails  +
		  NOX::Utils::Parameters
		  );

#else
  if (VERBOSE_NONLINEAR)
    noxParams_.sublist("Printing")
      .set("Output Information",
           NOX::Utils::Error+
           NOX::Utils::Warning +
           NOX::Utils::OuterIteration +
           NOX::Utils::StepperIteration
           );
  else
    noxParams_.sublist("Printing")
      .set("Output Information", NOX::Utils::Error);
#endif
  noxParams_.sublist("Printing")
    .set("Output Stream", Teuchos::RCP<std::ostream>(&dout(), false));

  // Defaults that are mode dependent
  switch (mode_)
  {
    case Nonlinear::TRANSIENT:
    // These values correspond to N_NLS_DampedNewton.C; see the constructor
    statusTestParams_.set("ABSTOL", 1.0e-6);
    statusTestParams_.set("RELTOL", 1.0e-2);
    statusTestParams_.set("DELTAXTOL", 0.33);
    statusTestParams_.set("RHSTOL", 1.0e-2);
    statusTestParams_.set("MAXSTEP", 20);
    noxParams_.set("Nonlinear Solver", "Line Search Based");
    noxParams_.sublist("Line Search").set("Method", "Full Step");
    break;
  case Nonlinear::HB_MODE:
    // These values correspond to N_NLS_DampedNewton.C; see the constructor
    statusTestParams_.set("ABSTOL", 1.0e-9);
    statusTestParams_.set("RELTOL", 1.0e-3);
    statusTestParams_.set("DELTAXTOL", 1.0);
    statusTestParams_.set("RHSTOL", 1.0e-4);
    statusTestParams_.set("MAXSTEP", 200);
    noxParams_.set("Nonlinear Solver", "Line Search Based");
    noxParams_.sublist("Line Search").set("Method", "Full Step");
    break;
  default:
    statusTestParams_.set("ABSTOL", 1.0e-12);
    statusTestParams_.set("RELTOL", 1.0e-3);
    statusTestParams_.set("DELTAXTOL", 1.0);
    statusTestParams_.set("RHSTOL", 1.0e-6);
    statusTestParams_.set("MAXSTEP", 200);
    noxParams_.set("Nonlinear Solver", "Line Search Based");
    noxParams_.sublist("Line Search").set("Method", "Full Step");
    noxParams_.sublist("Direction").sublist("Newton")
      .sublist("Linear Solver").set("Tolerance", 1.0e-12);
    break;
  }

  // Parameters that should always be set
  noxParams_.sublist("Line Search").sublist("Polynomial")
    .set("Recovery Step Type", "Last Computed Step");


  // Set default loca options in case this is a loca run.
  Teuchos::ParameterList& stepperList = locaParams_.sublist("Stepper");
  Teuchos::ParameterList& predictorList = locaParams_.sublist("Predictor");
  Teuchos::ParameterList& stepSizeList = locaParams_.sublist("Step Size");

  stepperList.set("Continuation Method", "Natural");
  stepperList.set("Skip df/dp", true);
  predictorList.set("Method", "Tangent");
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::~ParameterSet
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
ParameterSet::~ParameterSet()
{
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::setOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool ParameterSet::setOptions(const Util::OptionBlock& OB)
{
  // Parse the option block
  bool parseok = parseOptionBlock_(OB);
  if (!parseok)
  {
    return false;
  }

  isParamsSet_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::setOutputOptions
// Purpose       : Set output for parallel runs
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool ParameterSet::setOutputOptions(int myPID, int outputProcess)
{
  noxParams_.sublist("Printing").set("MyPID", myPID);
  noxParams_.sublist("Printing").set("Output Processor",
					      outputProcess);
  locaParams_.sublist("Utilities").set("MyPID", myPID);
  locaParams_.sublist("Utilities").set("Output Processor",
						outputProcess);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::createStatusTests
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool ParameterSet::createStatusTests(
  Parallel::Machine                     comm,
  TimeIntg::DataStore*                  data_store,
  Loader::NonlinearEquationLoader &     loader, 
  Linear::Solver &                      lsolver, 
  Linear::Vector *                      maskVectorPtr)          
{
  // Tests All - replaces tests 1-7
  bool isTransient = false;
  if (mode_ == Nonlinear::TRANSIENT)
  {
    isTransient = true;
  }

  Teuchos::RCP<XyceTests> allTests;

  // here we make the default XyceTests object
  allTests = Teuchos::rcp(new XyceTests(comm,
                                        isTransient,
                                        statusTestParams_.get("RHSTOL", 1.0e-6),
                                        Util::MachineDependentParams::MachineEpsilon(),
                                        data_store,
                                        statusTestParams_.get("ABSTOL", 1.0e-12),
                                        statusTestParams_.get("RELTOL", 1.0e-3),
                                        statusTestParams_.get("DELTAXTOL", 1.0),
                                        statusTestParams_.get("MAXSTEP", 200),
                                        1.0,
                                        0.9,
                                        0.5*Util::MachineDependentParams::DoubleMax(),
                                        1.0e-3,
                                        5,
                                        statusTestParams_.get("ENFORCEDEVICECONV", 1),
                                        statusTestParams_.get("SMALLUPDATETOL", 1.0e-6),
                                        &loader,
                                        &lsolver, 
                                        maskingFlag_,
                                        maskVectorPtr));

  tests_.push_back(allTests);
  comboPtr_->addStatusTest(allTests);
  isStatusTestsSet_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getStatusTests
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Teuchos::RCP<NOX::StatusTest::Generic> ParameterSet::getStatusTests()
{
  if (!isStatusTestsSet_)
  {
    Report::DevelFatal0().in("N_NLS::NOX::ParameterSet::getStatusTests")
      << "Status tests are not set!";
  }

  return comboPtr_;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getAllParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Teuchos::RCP<Teuchos::ParameterList> ParameterSet::getAllParams()
{
  return allParams_;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getNoxParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Teuchos::RCP<Teuchos::ParameterList> ParameterSet::getNoxParams()
{
  return Teuchos::rcp(&noxParams_, false);
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getLocaParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Teuchos::RCP<Teuchos::ParameterList> ParameterSet::getLocaParams()
{
  return Teuchos::rcp(&locaParams_, false);
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getDebugParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Teuchos::RCP<Teuchos::ParameterList> ParameterSet::getDebugParams()
{
  return Teuchos::rcp(&debugParams_, false);
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::unsupportedOption_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void ParameterSet::unsupportedOption_(const std::string& tag)
{
  Report::UserWarning0()
    << "Tag \"" << tag << "\" is unsupported by the NOX interface at this time.\n";
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::parseOptionBlock_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool ParameterSet::parseOptionBlock_(const Util::OptionBlock& OB)
{
  // RPP: Some parameters can't be immediately set in the nox
  // parameter list because they are in a sublist dependent upon
  // another parameter being set first.  We have to store these
  // parameters until the entire option block is parsed and then set
  // them in the correct sublist.  These parameters are listed below.
  int maxSearchStep = 2;
  int in_Forcing = 0;
  double AZ_tol = 1.0e-12;
  int recoveryStepType = 0;
  double recoveryStep = 1.0;
  int memory = 400;

  // Loop over all parameters in the option block
  for (Util::ParamList::const_iterator it_tpL = OB.begin();
       it_tpL != OB.end(); ++ it_tpL)
  {
    const std::string tag = it_tpL->uTag();


    // Parameters for nonlinear convergence tests
    if (tag == "ABSTOL")
    {
      statusTestParams_.set("ABSTOL", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "RELTOL")
    {
      statusTestParams_.set("RELTOL", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "DELTAXTOL")
    {
      statusTestParams_.set("DELTAXTOL", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "RHSTOL")
    {
      statusTestParams_.set("RHSTOL", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "MAXSTEP")
    {
      statusTestParams_.set("MAXSTEP", it_tpL->getImmutableValue<int>());
    }
    else if (tag == "SMALLUPDATETOL")
    {
      statusTestParams_.set("SMALLUPDATETOL", it_tpL->getImmutableValue<double>());
    }

    // Check devices for convergence (by calls to cktloader)
    else if (tag == "ENFORCEDEVICECONV")
    {
      statusTestParams_.set("ENFORCEDEVICECONV", it_tpL->getImmutableValue<int>());
    }
    else if (tag == "DEBUGLEVEL")
    {
      debugLevel_ = it_tpL->getImmutableValue<int>();
    }
    else if (tag == "SCREENOUTPUT")
    {
      screenOutputFlag_ = it_tpL->getImmutableValue<bool>();
    }
    else if (tag == "USEMASKING")
    {
      maskingFlag_ = it_tpL->getImmutableValue<bool>();
    }
    else if (tag == "DEBUGMINTIMESTEP")
    {
      debugMinTimeStep_ = it_tpL->getImmutableValue<int>();
    }
    else if (tag == "DEBUGMAXTIMESTEP")
    {
      debugMaxTimeStep_ = it_tpL->getImmutableValue<int>();
    }
    else if (tag == "DEBUGMINTIME")
    {
      debugMinTime_ = it_tpL->getImmutableValue<double>();
    }
    else if (tag == "DEBUGMAXTIME")
    {
      debugMaxTime_ = it_tpL->getImmutableValue<double>();
    }

    // Nonlinear Strategy
    else if (tag == "NLSTRATEGY")
    {
      int val = it_tpL->getImmutableValue<int>();
      if (val == 0)		// Newton
      {
        noxParams_.set("Nonlinear Solver", "Line Search Based");
        noxParams_.sublist("Direction").set("Method", "Newton");
      }
      else if (val == 1) 	// Steepest descent
      {
        noxParams_.set("Nonlinear Solver", "Line Search Based");
        noxParams_.sublist("Direction")
          .set("Method", "Steepest Descent");
      }
      else if (val == 2)	// Trust Region
      {
        noxParams_.set("Nonlinear Solver", "Trust Region Based");
      }
      else if (val == 3) 	// Modified Newton
      {
        noxParams_.set("Nonlinear Solver", "Line Search Based");
        noxParams_.sublist("Direction")
          .set("Method", "Modified-Newton");
      }
      else if (val == 4)	// BFGS
      {
        noxParams_.set("Nonlinear Solver", "Line Search Based");
        noxParams_.sublist("Direction").set("Method", "Quasi-Newton");
      }
      else if (val == 5)	// Broyden
      {
        noxParams_.set("Nonlinear Solver", "Line Search Based");
        noxParams_.sublist("Direction").set("Method", "Broyden");
      }
      else if (val == 6)	// Tensor
      {
        noxParams_.set("Nonlinear Solver", "Tensor Based");
        noxParams_.sublist("Direction").set("Method", "Tensor");
        noxParams_.sublist("Direction").sublist("Tensor")
          .sublist("Linear Solver").set("Compute Step", "Newton");
        noxParams_.sublist("Direction").sublist("Tensor")
          .sublist("Linear Solver").set("Reorthogonalize", "Always");
        noxParams_.sublist("Line Search").set("Method", "Tensor");
        noxParams_.sublist("Line Search").sublist("Tensor")
          .set("Submethod", "Full Step");
      }
      else if (val == 7)	// Fast Newton Direction
      {
        // RPP: No longer supported
        Report::UserWarning0() 
          << "NLStrategy = 7 is no longer supported.";
      }
      else if (val == 8)	// Newton/Steepest Descent Combo Direction
      {
        // RPP: No longer supported
        Report::UserWarning0() 
          << "NLStrategy = 8 is no longer supported.";
      }
      else
      {
        Report::UserWarning0() 
          << "NLStrategy is not found!";
      }
    }

    // Line search method
    else if (tag == "SEARCHMETHOD")
    {
      int val = it_tpL->getImmutableValue<int>();
      if (val == 0)
      {
        noxParams_.sublist("Line Search").set("Method", "Full Step");
      }
      else if (val == 1)
      {
        noxParams_.sublist("Line Search").set("Method", "Backtrack");
      }
      else if (val == 2)
      {
        noxParams_.sublist("Line Search").set("Method", "Polynomial");
        noxParams_.sublist("Line Search").sublist("Polynomial")
          .set("Interpolation Type", "Quadratic");
      }
      else if (val == 3)
      {
        noxParams_.sublist("Line Search").set("Method", "Polynomial");
        noxParams_.sublist("Line Search").sublist("Polynomial")
          .set("Interpolation Type", "Cubic");
      }
      else if (val == 4)
      {
        noxParams_.sublist("Line Search")
          .set("Method", "More'-Thuente");
      }
      else
      {
        Report::UserWarning0() 
          << "ParameterSet::parseOptionBlock_ - "
          << "SEARCHMETHOD = " << val
          << " not supported by Xyce at this time.";
      }
    }

    // Trust Region auxiliary parameters
    else if (tag == "TRMINRADIUS")
    {
      noxParams_.sublist("Trust Region")
        .set("Minimum Trust Region Radius", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "TRMAXRADIUS")
    {
      noxParams_.sublist("Trust Region")
        .set("Maximum Trust Region Radius", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "TRMINIMPROVEMENTRATIO")
    {
      noxParams_.sublist("Trust Region")
        .set("Minimum Improvement Ratio", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "TRCONTRACTIONRATIO")
    {
      noxParams_.sublist("Trust Region")
        .set("Contraction Trigger Ratio", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "TRCONTRACTIONFACTOR")
    {
      noxParams_.sublist("Trust Region")
        .set("Contraction Factor", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "TREXPANSIONRATIO")
    {
      noxParams_.sublist("Trust Region")
        .set("Expansion Trigger Ratio", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "TREXPANSIONFACTOR")
    {
      noxParams_.sublist("Trust Region")
        .set("Expansion Factor", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "TRRECOVERYSTEP")
    {
      noxParams_.sublist("Trust Region")
        .set("Recovery Step", it_tpL->getImmutableValue<double>());
    }

    // RPP: Why this is here???
    else if (tag == "NOX")
    {
      // do nothing (this option is handled in the manager)
    }

    // LOCA Continuation control
    // 0 = Nox solve (no continuation)
    // 1 = Natural Parameter Continuation
    // 2 = Mosfet Specific Dual Parameter Continuation
    // 3 = gmin stepping.
    // 33 = Artificial Parameter Continuation
    else if (tag == "CONTINUATION")
    {
      continuationSpecified_=true;
      if (it_tpL->isNumeric())
      {
        noxSolver = it_tpL->getImmutableValue<int>();
      }
      else
      {
        ExtendedString p(it_tpL->stringValue());
        p.toUpper();
        if (p.substr(0,4) == "STAN")
        {
          noxSolver = 0;
        }
        else if (p.substr(0,3) == "NAT")
        {
          noxSolver = 1;
        }
        else if (p.substr(0,3) == "MOS")
        {
          noxSolver = 2;
        }
        else if (p.substr(0,4) == "GMIN")
        {
          noxSolver = 3;
        }
        else if (p.substr(0,6) == "NEWMOS")
        {
          noxSolver = 4;
        }
        else if (p.substr(0,9) == "BSIM3INV1")
        {
          noxSolver = 5;
        }
        else if (p.substr(0,9) == "BSIM3INV2")
        {
          noxSolver = 6;
        }
        else if (p.substr(0,9) == "BLOCKGAIN")
        {
          noxSolver = 7;
        }
        else if (p.substr(0,4) == "TEST")
        {
          noxSolver = 8;
        }
        else if (p.substr(0,6) == "PSEUDO")
        {
          noxSolver = 9;
        }
        else if (p.substr(0,5) == "POWER")
        {
          noxSolver = 10;
        }
        else if (p.substr(0,3) == "ART")
        {
          noxSolver = 33;
        }
        else if (p.substr(0,10) == "SOURCESTEP")
        {
          noxSolver = 34;
        }
        else
        {
          Report::DevelFatal0() 
            << "Unknown specification in .options for 'continuation': "
            << it_tpL->stringValue();
        }

	// Bail out if option 33 (artif. homotopy) is chosen but
	// required support is not built in.
#ifndef Xyce_NOX_LOCA_ARTIFICIAL_HOMOTOPY_SUPPORT
        if (noxSolver == 33)
        {
          Report::UserFatal0() << "Nonlinear Solver (NOX::Interface) Artificial parameter continuation requires "
                                     << "building xyce with the define: -DXyce_NOX_LOCA_ARTIFICIAL_HOMOTOPY_SUPPORT to "
                                     << "allow LOCA to augment the diagonal of Jacobian! Either rebuild Xyce or do not "
                                     << "run Xyce with \"continuation=33\"";
        }
#endif
      }
    }

    // Parameters that can't be set in the list until all options
    // have been parsed
    else if (tag == "MAXSEARCHSTEP")
    {
      maxSearchStep = it_tpL->getImmutableValue<int>();
    }
    else if (tag == "IN_FORCING")
    {
      in_Forcing = it_tpL->getImmutableValue<int>();
    }
    else if (tag == "AZ_TOL")
    {
      AZ_tol = it_tpL->getImmutableValue<double>();
    }
    else if (tag == "RECOVERYSTEPTYPE")
    {
      recoveryStepType = it_tpL->getImmutableValue<int>();
    }
    else if (tag == "RECOVERYSTEP")
    {
      recoveryStep = it_tpL->getImmutableValue<double>();
    }

    // Parameters that can't be set in the list until all options
    // have been parsed
    else if (tag == "MAXSEARCHSTEP")
    {
      maxSearchStep = it_tpL->getImmutableValue<int>();
    }
    else if (tag == "IN_FORCING")
    {
      in_Forcing = it_tpL->getImmutableValue<int>();
    }
    else if (tag == "AZ_TOL")
    {
      AZ_tol = it_tpL->getImmutableValue<double>();
    }
    else if (tag == "RECOVERYSTEPTYPE")
    {
      recoveryStepType = it_tpL->getImmutableValue<int>();
    }
    else if (tag == "RECOVERYSTEP")
    {
      recoveryStep = it_tpL->getImmutableValue<double>();
    }
    // Warn user about unrecognized solver option
    else
    {
      Report::UserWarning() << tag << " is not a recognized nonlinear solver option.\n" << std::endl;
    }

  } // end loop over all options in block


  /*
     Set parameters that are dependent upon other parameters

     RPP: Some parameters can't be immediately set in the nox
     parameter list because they are in a sublist dependent upon
     another parameter being set first.  We have to store these
     parameters until the entire option block is parsed and then set
     them in the correct sublist.
  */
  std::string directionMethod =
    noxParams_.sublist("Direction").get("Method", "Newton");

  std::string lineSearchMethod =
    noxParams_.sublist("Line Search").get("Method", "Full Step");

  // MAXSEARCHSTEP
  noxParams_.sublist("Line Search").sublist(lineSearchMethod)
    .set("Max Iters", maxSearchStep);

  // In_FORCING
  if (in_Forcing == 0)
  {
    noxParams_.sublist("Direction").sublist(directionMethod)
      .set("Forcing Term Method", "Constant");
  }
  else if (in_Forcing == 1)
  {
    noxParams_.sublist("Direction").sublist(directionMethod)
      .set("Forcing Term Method", "Type 1");
    noxParams_.sublist("Direction").sublist("Newton")
      .set("Forcing Term Minimum Tolerance", AZ_tol);
  }
  else
  {
    noxParams_.sublist("Direction").sublist(directionMethod)
      .set("Forcing Term Method", "Type 2");
    noxParams_.sublist("Direction").sublist("Newton")
      .set("Forcing Term Minimum Tolerance", AZ_tol);
  }

  // RECOVERYSTEPTYPE
  if (recoveryStepType == 1)
  {
    // Recovery the NOX way
    noxParams_.sublist("Line Search").sublist(lineSearchMethod)
      .set("Recovery Step Type", "Constant");
  }
  else
  {
    // Recovery the Xyce way
    noxParams_.sublist("Line Search").sublist(lineSearchMethod)
      .set("Recovery Step Type", "Last Computed Step");
  }

  // RECOVERYSTEP
  noxParams_.sublist("Line Search").sublist(lineSearchMethod)
    .set("Recovery Step", recoveryStep);

  // MEMORY
  if ((directionMethod == "Quasi-Newton") ||
      (directionMethod == "Broyden"))
    noxParams_.sublist("Direction").sublist(directionMethod)
      .set("Memory", memory);


  return true;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::setLocaOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool ParameterSet::setLocaOptions(const Util::OptionBlock& OB, bool saveCopy)
{
  Teuchos::ParameterList& stepperList = locaParams_.sublist("Stepper");
  Teuchos::ParameterList& predictorList = locaParams_.sublist("Predictor");
  Teuchos::ParameterList& stepSizeList = locaParams_.sublist("Step Size");

  bool stepperGiven=false;
  bool predictorGiven=false;

  if (saveCopy)
  {
    savedLocaOptions_ = true;
    savedLocaOB_ = OB; // save a copy to re-assert defaults later, if needed.
  }

  for (Util::ParamList::const_iterator it_tpL = OB.begin();
       it_tpL != OB.end(); ++ it_tpL)
  {
    const std::string tag = it_tpL->uTag();
    std::string baseTag=tag.substr(0,8);
    bool isVectorParam=false;
    if (baseTag == "CONPARAM" ||
        baseTag == "MINVALUE" ||
        baseTag == "MAXVALUE" ||
        baseTag == "INITIALV" ||
        baseTag == "INITIALS" ||
        baseTag == "MINSTEPS" ||
        (baseTag == "MAXSTEPS" && tag.substr(0,11)=="MAXSTEPSIZE")||
        baseTag == "AGGRESSI")
    {
      if (baseTag == "INITIALV")
      {
        baseTag = "INITIALVALUE";
      }
      else if (baseTag == "INITIALS")
      {
        baseTag = "INITIALSTEPSIZE";
      }
      else if (baseTag == "MINSTEPS")
      {
        baseTag = "MINSTEPSIZE";
      }
      else if (baseTag == "MAXSTEPS")
      {
        baseTag = "MAXSTEPSIZE";
      }
      else if (baseTag == "AGGRESSI")
      {
        baseTag = "AGGRESSIVENESS";
      }

      std::string num=tag.substr(baseTag.size(),tag.size()-baseTag.size());
      int index=ExtendedString(num).Ival();
      vectorParams[baseTag].push_back(*it_tpL);
      isVectorParam=true;
    }
    if (tag == "STEPPER")
    {
      stepperGiven=true;
      if (it_tpL->isNumeric())
      {
        int iType = it_tpL->getImmutableValue<int>();
        if (iType == 1)
        {
          stepperList.set("Continuation Method", "Arc Length");
          stepperList.set("Skip df/dp", false);
        }
        else
        {
          stepperList.set("Continuation Method", "Natural");
          stepperList.set("Skip df/dp", true);
        }
      }
      else
      {
        ExtendedString p(it_tpL->stringValue());
        p.toUpper();
        if (p.substr(0,3) == "ARC")
        {
          stepperList.set("Continuation Method", "Arc Length");
          stepperList.set("Skip df/dp", false);
        }
        else if (p.substr(0,3) == "NAT")
        {
          stepperList.set("Continuation Method", "Natural");
          stepperList.set("Skip df/dp", true);
        }
        else
        {
          Report::DevelFatal0()
            << "Unknown specification in .options for 'stepper': "
            << it_tpL->stringValue()
            <<  ".  Legal choices are ARCLENGTH or NATURAL, which may be abbreviated to three characters.";
        }
      }
    }
    else if (tag == "PREDICTOR")
    {
      predictorGiven=true;
      if (it_tpL->isNumeric())
      {
        int iType = it_tpL->getImmutableValue<int>();
        if (iType == 1)
        {
          predictorList.set("Method", "Tangent");
        }
        else if (iType == 2)
        {
          predictorList.set("Method", "Secant");
        }
        else if (iType == 3)
        {
          predictorList.set("Method", "Random");
        }
        else
        {
         predictorList.set("Method", "Constant");
        }
      }
      else
      {
        ExtendedString p(it_tpL->stringValue());
        p.toUpper();
        if (p.substr(0,3) == "TAN")
        {
          predictorList.set("Method", "Tangent");
        }
        else if (p.substr(0,3) == "SEC")
        {
          predictorList.set("Method", "Secant");
        }
        else if (p.substr(0,3) == "RAN")
        {
          predictorList.set("Method", "Random");
        }
        else if (p.substr(0,3) == "CON")
        {
         predictorList.set("Method", "Constant");
        }
        else
        {
          Report::DevelFatal0()
            << "Unknown specification in .options for 'predictor': "
            <<  it_tpL->stringValue()
            <<  ".  Legal choices are TANGENT, SECANT, RANDOM, CONSTANT, which may be abbreviated to three characters.";
        }
      }
    }
    // if we're the very first instance of one of these vector params,
    // then set our stepperList (this is the old behavior from the
    // "alternate" vector handling
    else if (tag == "CONPARAM1") // continuation parameter
    {
      stepperList.set("Continuation Parameter", it_tpL->stringValue());
    }
    else if (tag == "INITIALVALUE1")
    {
      stepperList.set("Initial Value", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "MINVALUE1")
    {
      gstepping_min_value_ = it_tpL->getImmutableValue<double>();
      stepperList.set("Min Value", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "RESIDUALCONDUCTANCE")
    {
      gstepping_minimum_conductance_ = it_tpL->getImmutableValue<double>();
      if (gstepping_minimum_conductance_ > 0)
      {
        Report::UserWarning0() << "A non-zero value for the GMIN Stepping residual conductance has been specified (RESIDUALCONDUCTANCE= " << it_tpL->stringValue() << ")." << std::endl
                               << "This option should never be used unless absolutely necessary to obtain an initial condition for transient runs with ill-posed DC operating points.  The operating point obtained by GMIN Stepping will not be a valid steady state condition for the circuit as defined in the netlist, but might possibly produce a reasonable initial condition for transient runs.";
      }
    }
    else if (tag == "INITIALSTEPSIZE1")
    {
      stepSizeList.set("Initial Step Size", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "MINSTEPSIZE1")
    {
      stepSizeList.set("Min Step Size", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "MAXSTEPSIZE1")
    {
      stepSizeList.set("Max Step Size", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "AGGRESSIVENESS1")
    {
      stepSizeList.set("Aggressiveness", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "MAXVALUE1")
    {
      stepperList.set("Max Value", it_tpL->getImmutableValue<double>());
    }
    else if (tag == "BIFPARAM") // bifurcation parameter
    {
      stepperList.set("Bifurcation Parameter", it_tpL->stringValue());
    }
    else if (tag == "MAXSTEPS")
    {
      stepperList.set("Max Steps", it_tpL->getImmutableValue<int>());
    }
    else if (tag == "MAXNLITERS")
    {
      stepperList.set("Max Nonlinear Iterations", it_tpL->getImmutableValue<int> ());
    }
    else if (tag == "STEPCONTROL")
    {
      if (it_tpL->isNumeric())
      {
        int iType = it_tpL->getImmutableValue<int>();
        if (iType == 1)
          stepSizeList.set("Method", "Adaptive");
        else
          stepSizeList.set("Method", "Constant");
      }
      else
      {
        ExtendedString p(it_tpL->stringValue());
        p.toUpper();
        if (p.substr(0,3) == "ADA")
        {
          stepSizeList.set("Method", "Adaptive");
        }
        else if (p.substr(0,3) == "CON")
        {
          stepSizeList.set("Method", "Constant");
        }
        else
        {
          Report::DevelFatal0()
            << "Unknown specification in .options for 'stepcontrol': "
            << it_tpL->stringValue()
            << ".  Legal choices are ADAPTIVE or CONSTANT, which may be abbreviated to three characters.";
        }
      }
    }
    else if (tag == "POWERNODE") // continuation parameter
    {
      stepperList.set("Power Node", it_tpL->stringValue());
    }
    else if (tag == "VOLTAGELIST")
    {
      ExtendedString p(it_tpL->stringValue());
      p.toUpper();
      if (p.substr(0,3) == "DOF") // DOFS - degrees of freedom
      {
        voltageListType_ = VLT_DOFS;
      }
      else if (p.substr(0,3) == "NOD") // NODES - voltage nodes
      {
        voltageListType_ = VLT_Node;
      }
      else
      {
          Report::DevelFatal0()
            << "Unknown specification in .options loca for 'voltagelist': "
            << it_tpL->stringValue()
            << ".  Legal choices are DOFS or NODES, which may be abbreviated to three characters.";
      }
    }
    else if(tag == "VOLTAGESCALEFACTOR")
    {
      voltageScaleFactor_ = it_tpL->getImmutableValue<double>();
    }
    // Start of the new parameter section
    else if (std::string(tag,0,10) == "PARAMLIST")
    {
      // don't know what to do yet.
      dout() << "tag = " << tag << std::endl;
    }
    else
    {
      // if "isVectorParam" we've already handled this at the beginning
      // of the loop.  Otherwise it's an unrecognized parameter.
      if (!isVectorParam)
      {
        Report::UserWarning0()
          << tag << " is not a recognized loca option.\n";
      }
    }
  }

  // Insure that the correct defaults are set.  Sometimes the LOCA
  // defaults change in the LOCA library from release to release.  However
  // Xyce's defaults should not necessarily change in that case.
  if (!stepperGiven)
  {
    stepperList.set("Continuation Method", "Natural");
    stepperList.set("Skip df/dp", true);
  }

  if (!predictorGiven)
  {
    predictorList.set("Method", "Tangent");
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getVectorParam
// Purpose       : Obtain a parameter specified in the option line in vector syntax
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 02/02/06
//-----------------------------------------------------------------------------
bool ParameterSet::getVectorParam (const std::string & tag, int index, double & value)
{
  if (vectorParams.find(tag) != vectorParams.end() &&
      vectorParams[tag].size() > index)
  {
    value = vectorParams[tag][index].getImmutableValue<double>();
    return true;
  }
  else
  {
    return false;
  }
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getVectorParam
// Purpose       : Obtain a parameter specified in the option line in vector syntax
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 02/02/06
//-----------------------------------------------------------------------------
bool ParameterSet::getVectorParam (const std::string & tag, int index, std::string & value)
{
  if (vectorParams.find(tag) != vectorParams.end() &&
      vectorParams[tag].size() > index)
  {
    value = vectorParams[tag][index].stringValue();
    return true;
  }
  else
  {
    return false;
  }
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getVectorParamSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int ParameterSet::getVectorParamSize(const std::string& tag)
{
  if (vectorParams.find(tag) != vectorParams.end())
  {
    return static_cast<int>(vectorParams[tag].size());
  }
  else
  {
    Report::DevelFatal0().in("ParameterSet::getVectorParam")
      <<  "the parameter \""
      <<  tag
      << "\" is required for parameter continuation!";
  }

  return -1;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getStatusTestReturnCode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int ParameterSet::getStatusTestReturnCode() const
{
  // Get the main Xyce Test
  Teuchos::RCP<XyceTests> testPtr =
    Teuchos::rcp_dynamic_cast<XyceTests>(tests_[1]);
  if (Teuchos::is_null(testPtr))
  {
    Report::DevelFatal0().in("ParameterSet::getStatusTestReturnCode")
      << "Dynamic cast on Xyce Tests check failed.";
  }

  return testPtr->getXyceReturnCode();
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::setStatusTestReturnCode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void ParameterSet::setStatusTestReturnCodes
  (const Nonlinear::ReturnCodes & retCodesTmp)
{
  // Get the main Xyce Test
  Teuchos::RCP<XyceTests> testPtr =
    Teuchos::rcp_dynamic_cast<XyceTests>(tests_[1]);
  if (Teuchos::is_null(testPtr))
  {
    Report::DevelFatal0().in("ParameterSet::setStatusTestReturnCode")
      << "Dynamic cast on Xyce Tests check failed.";
  }

  return testPtr->setReturnCodes (retCodesTmp);
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getMaxNormF
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
double ParameterSet::getMaxNormF() const
{
  // Get the main Xyce Test
  Teuchos::RCP<XyceTests> testPtr =
    Teuchos::rcp_dynamic_cast<XyceTests>(tests_[1]);
  if (Teuchos::is_null(testPtr))
  {
    Report::DevelFatal0().in("ParameterSet::getMaxNormF")
      << "Dynamic cast on Xyce Tests check failed.";
  }

  return testPtr->getMaxNormF();
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getMaxNormFindex
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int ParameterSet::getMaxNormFindex () const
{
  // Get the main Xyce Test
  Teuchos::RCP<XyceTests> testPtr =
    Teuchos::rcp_dynamic_cast<XyceTests>(tests_[1]);
  if (Teuchos::is_null(testPtr))
  {
    Report::DevelFatal0().in("ParameterSet::getMaxNormFindex")
      << "Dynamic cast on Xyce Tests check failed.";
  }

  return testPtr->getMaxNormFindex ();
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getNoxSolverType
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int ParameterSet::getNoxSolverType() const
{
  return noxSolver;
}

void ParameterSet::setNoxSolverType(int type) 
{
  noxSolver = type;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::createAugmentLinearSystem
// Purpose       : creates an AugmentLinSys strategy object.
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 1416
// Creation Date : 03/08/06
//-----------------------------------------------------------------------------
Teuchos::RCP<AugmentLinSys>
ParameterSet::createAugmentLinearSystem(Linear::System* ls) const
{
  Teuchos::RCP<AugmentLinSys> als;

  if (noxSolver == 9)
  {
    if (voltageScaleFactor_ == 1.0)
    {
      als = Teuchos::rcp( new 
			  AugmentLinSysPseudoTransient(ls->builder().createSolnColoring(),
						       ls->getRHSVector()) );
    }
    else
    {
      als = Teuchos::rcp( new 
			  AugmentLinSysPseudoTransient(ls->builder().createSolnColoring(),
						       ls->getRHSVector(),
						       true,
						       voltageScaleFactor_) );
    }

  }
  else if(noxSolver == 1 || noxSolver == 3)
  {
    if (voltageListType_ == VLT_DOFS)
    {
      als = Teuchos::rcp( new GStepping(GStepping::NLT_AllVoltageUnknowns,
                                        ls->builder().createSolnColoring(),
					ls->getRHSVector(),
					gstepping_min_value_,
                                        gstepping_minimum_conductance_) );
    }
    else
    {
      als = Teuchos::rcp( new GStepping(GStepping::NLT_VoltageNodes,
                                        ls->builder().vnodeGIDVec(),
				        ls->getRHSVector(),
				        gstepping_min_value_,
                                        gstepping_minimum_conductance_) );
    }
  }
  else
  {
    Report::DevelFatal0().in("ParameterSet::createAugmentLinearSystem")
      << "- The \'continuation\' "
      << "parameter in the .options nox list must be set to PSEUDO or NATURAL for "
      <<  "this function to be called!";
  }

  return als;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::createAugmentLinearSystem
// Purpose       : creates an AugmentLinSys strategy object for .IC with
//                 or without gmin stepping
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/08/07
//-----------------------------------------------------------------------------
Teuchos::RCP<AugmentLinSys>
ParameterSet::createAugmentLinearSystem(Linear::System* ls, IO::InitialConditionsData::NodeLidValueMap & op,
   bool gminStepping) const
{
  Teuchos::RCP<AugmentLinSys> als;

  if (gminStepping==false)
  {
    als = Teuchos::rcp( new AugmentLinSysIC(op, 
                                            ls->builder().createInitialConditionColoring(),
                                            ls->getRHSVector() ) );
  }
  else
  {
    if (voltageListType_ == VLT_DOFS)
    {
      als = Teuchos::rcp( new AugmentLinSysIC_Gmin(
                AugmentLinSysIC_Gmin::NLT_AllVoltageUnknowns,
                op,
                ls->builder().createInitialConditionColoring(),
                ls->builder().createSolnColoring(),
                ls->getRHSVector(),
                gstepping_min_value_,
                gstepping_minimum_conductance_) );
    }
    else
    {
      als = Teuchos::rcp( new AugmentLinSysIC_Gmin(
                AugmentLinSysIC_Gmin::NLT_VoltageNodes,
                op,
                ls->builder().createInitialConditionColoring(),
                ls->builder().vnodeGIDVec(),
                ls->getRHSVector(),
	        gstepping_min_value_,
                gstepping_minimum_conductance_) );
    }
  }

  return als;
}

}}} // namespace N_NLS_NOX

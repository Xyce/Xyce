//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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

#include <sstream>

#include <LOCA_GlobalData.H>
#include <LOCA_StatusTest_Wrapper.H>
#include <NOX_Solver_Factory.H>
#include <N_ANP_AnalysisManager.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CmdParse.h>
#include <N_IO_OutputMgr.h>
#include <N_LAS_Builder.h>
#include <N_LAS_QueryUtil.h>
#include <N_LAS_Solver.h>
#include <N_LAS_System.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_NLS_LOCA_Group.h>
#include <N_NLS_NOX_AugmentLinSys.h>
#include <N_NLS_NOX_Group.h>
#include <N_NLS_NOX_Interface.h>
#include <N_NLS_NOX_PseudoTransientSolver.h>
#include <N_NLS_NOX_PseudoTransientTest.h>
#include <N_NLS_NOX_SharedSystem.h>
#include <N_NLS_NOX_XyceTests.h>
#include <N_NLS_SensitivityResiduals.h>
#include <N_UTL_FeatureTest.h>
#include <N_PDS_Comm.h>


// ----------   NOX Includes   ----------
#include <LOCA.H>

// -----------  Forward declarations  -------
namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

//-----------------------------------------------------------------------------
// Function      : Interface::Interface
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Interface::Interface(const IO::CmdParse & cp)
  : Nonlinear::NonLinearSolver(cp),
    dcParams_(Nonlinear::DC_OP),
    ICspecified_(false),
    NODESETspecified_(false),
    transientParams_(Nonlinear::TRANSIENT),
    hbParams_(Nonlinear::HB_MODE),
    nlpParams_(Nonlinear::DC_NLPOISSON),
    sharedSystemPtr_(0),
    mode_(Nonlinear::DC_OP),
    usemode_(true),
    lastParametersMode_(Nonlinear::DC_OP),
    parametersMode_(Nonlinear::DC_OP),
    copiedGroupFlag_(false),
    isFirstContinuationParam_(true),
    firstSolveComplete_(false),
    iParam_(0)
{
}

//-----------------------------------------------------------------------------
// Function      : Interface::~Interface
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Interface::~Interface()
{
  delete sharedSystemPtr_;
  if (!Teuchos::is_null(globalDataPtr_))
  {
    LOCA::destroyGlobalData(globalDataPtr_);
  }
}

//-----------------------------------------------------------------------------
// Function      : Interface::setOptions
// Purpose       : Passes option block corresponding to "NONLIN" onto
//                 nonlinear solver. These parameters set convergence
//                 tolerances, the type of method used, and so on.
// Special Notes :
// Return Type   : boolean
//
// See also      : setTranOptions, setAnalysisMode
//
// - Input Arguments -
//
//    OB         : Option block containing options corresponding to
//                 "NONLIN" in the netlist.
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::setOptions(const Util::OptionBlock& OB)
{
  return dcParams_.setOptions(OB);
}


//-----------------------------------------------------------------------------
// Function      : Interface::setTranOptions
//
// Purpose       : Passes option block corresponding to "NONLIN-TRAN" onto
//                 nonlinear solver. This affects the settings used when
//                 the mode is Transient.
//
// Special Notes :
// Return Type   : boolean
//
// See also      : setOptions, setAnalysisMode
//
// - Input Arguments -
//    OB         : Option block containing options corresponding to
//                 "NONLIN-TRAN" in the netlist.
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::setTranOptions(const Util::OptionBlock& OB)
{
  return transientParams_.setOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : Interface::setHBOptions
// Purpose       :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::setHBOptions(const Util::OptionBlock& OB)
{
  return hbParams_.setOptions(OB);
}


//-----------------------------------------------------------------------------
// Function      : Interface::setNLPOptions
// Purpose       :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::setNLPOptions(const Util::OptionBlock& OB)
{
  return nlpParams_.setOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : Interface::setLocaOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::setLocaOptions(const Util::OptionBlock& OB)
{
  hbParams_.setLocaOptions(OB);
  return dcParams_.setLocaOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : Interface::setICOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::setICOptions(const Util::OptionBlock& OB)
{
  ICspecified_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Interface::setNodeSetOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::setNodeSetOptions(const Util::OptionBlock& OB)
{
  NODESETspecified_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Interface::initializeAll
//
// Purpose       : Called after all register and set functions.
//                 Once the various registrations have taken place,
//                 this function sets the remaining pointers.
//
// Special Notes :
// Derived Notes : This function also calls the base object initializeAll.
//
// Special Notes:  This function obtains the solution, temporary solution and
//                 f vectors from the LAS system class.
//
// Return Type   : boolean
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::initializeAll()
{
  // Use the base class initialization
  bool initok = Nonlinear::NonLinearSolver::initializeAll();
  if (!initok)
  {
    return false;
  }

  // Set the processor ID for printing in the nonlinear solver
  // For now, processor 0 will output information to the screen
  int myPID = pdsMgrPtr_->getPDSComm()->procID();
  dcParams_.setOutputOptions(myPID, 0);
  transientParams_.setOutputOptions(myPID, 0);

  hbParams_.setOutputOptions(myPID, 0);

  // Set up the status tests
  bool testsok = false;
  testsok =
    dcParams_.createStatusTests(pdsMgrPtr_->getPDSComm()->comm(), dsPtr_, getLoader(), lasSysPtr_->getDeviceMaskVector()) &&
    transientParams_.createStatusTests(pdsMgrPtr_->getPDSComm()->comm(), dsPtr_, getLoader(), lasSysPtr_->getDeviceMaskVector()) &&
    hbParams_.createStatusTests(pdsMgrPtr_->getPDSComm()->comm(), dsPtr_, getLoader(), lasSysPtr_->getDeviceMaskVector());
    nlpParams_.createStatusTests(pdsMgrPtr_->getPDSComm()->comm(), dsPtr_, getLoader(), lasSysPtr_->getDeviceMaskVector());

  if (!testsok)
    return false;

  // Set up any linear solver options
  // We only set the tolerance if adaptive forcing is being used.
  setAZ_Tol_DC = false;
  setAZ_Tol_Transient = false;
  if (!(dcParams_.getNoxParams()->sublist("Direction").sublist("Newton")
        .get("Forcing Term Method","Constant") == std::string("Constant")))
  {
    setAZ_Tol_DC = true;
  }
  if (!(transientParams_.getNoxParams()->sublist("Direction").sublist("Newton")
        .get("Forcing Term Method", "Constant") == std::string("Constant")))
  {
    setAZ_Tol_Transient = true;
  }
  if (!(hbParams_.getNoxParams()->sublist("Direction").sublist("Newton")
        .get("Forcing Term Method", "Constant") == std::string("Constant")))
  {
    setAZ_Tol_DC = true;
  }
  if (!(nlpParams_.getNoxParams()->sublist("Direction").sublist("Newton")
        .get("Forcing Term Method", "Constant") == std::string("Constant")))
  {
    setAZ_Tol_Transient = true;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Interface::spiceStrategy
//
// Purpose       : This function attempts a stdNewton solve first, and if that
//                 fails then does gmin stepping.  This is only done for DCOP
//                 calculations.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 5/7/2014
//-----------------------------------------------------------------------------
int Interface::spiceStrategy ( ParameterSet* paramsPtr )
{
  // Save a copy of nextSolutionPtr before stdNewtonSolve, just in case
  // of failure, so we can start gmin and source stepping with the same vector.
  // NOTE:  This is necessary if .IC statements are being enforced.
  Linear::Vector * initVec = dsPtr_->nextSolutionPtr->cloneCopyVector();

  groupPtr_->setNonContinuationFlag (true);
  int isuccess = stdNewtonSolve (paramsPtr);
 
  if (isuccess < 0) // attempt gmin stepping.
  {
    int saveSolverType=paramsPtr->getNoxSolverType();
    paramsPtr->setNoxSolverType(3);
    groupPtr_->setNonContinuationFlag (false);

    rhsVectorPtr_->putScalar(0.0);
    NewtonVectorPtr_->putScalar(0.0);
    gradVectorPtr_->putScalar(0.0);

    dsPtr_->setZeroHistory();

    // Copy saved initial solution vector.
    (*dsPtr_->nextSolutionPtr) = *initVec;
    Vector tmpVec(*dsPtr_->nextSolutionPtr, *lasSysPtr_);
    groupPtr_->setX(tmpVec);

    sharedSystemPtr_->reset(*dsPtr_->nextSolutionPtr,
        *rhsVectorPtr_,
        *jacobianMatrixPtr_,
        *NewtonVectorPtr_,
        *gradVectorPtr_,
        *lasSysPtr_,
        *this);

    isuccess=gminSteppingSolve ( paramsPtr );

    if (isuccess < 0)
    {
      paramsPtr->setNoxSolverType(34);
      groupPtr_->setNonContinuationFlag (false);
        
      rhsVectorPtr_->putScalar(0.0);
      NewtonVectorPtr_->putScalar(0.0);
      gradVectorPtr_->putScalar(0.0);

      dsPtr_->setZeroHistory();

      // Copy saved initial solution vector.
      (*dsPtr_->nextSolutionPtr) = *initVec;
      Vector tmpVec(*dsPtr_->nextSolutionPtr, *lasSysPtr_);
      groupPtr_->setX(tmpVec);
      
      sharedSystemPtr_->reset(*dsPtr_->nextSolutionPtr,
                              *rhsVectorPtr_,
                              *jacobianMatrixPtr_,
                              *NewtonVectorPtr_,
                              *gradVectorPtr_,
                              *lasSysPtr_,
                              *this);
      
      isuccess=sourceSteppingSolve ( paramsPtr );
      paramsPtr->setNoxSolverType(saveSolverType);

      nonlinearEquationLoader_->resetScaledParams();
    }
  }

  delete initVec;

  return isuccess;

}

//-----------------------------------------------------------------------------
// Function      : Interface::stdNewtonSolve
// Purpose       : 
// Special Notes : This corresponds to 
//                    .options nonlin continuation=stan or 
//                    .options nonlin continuation=0
// Scope         : public
// Creator       : Roger Pawlowski, Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::stdNewtonSolve ( ParameterSet* paramsPtr )
{
  bool usedIC=false;
  bool usedNODESET=false;

  // RPP: needed for tensor method.  DO NOT UNCOMMENT!
  // This somehow breaks the one-shot test circuit.
  // ERK:  it breaks it (and probably many other circuits) b/c voltage
  //       limiting is fragile and cannot handle extra (unexpected)
  //       F-loads.
  //groupPtr_->computeF();

  // Set up nox nonlinear solver.  solverPtr is only needed for
  // non-LOCA (i.e. just Newton's method, no continuation) solves.
  if ( Teuchos::is_null(solverPtr_) || ((usemode_) && (lastParametersMode_ != parametersMode_)) )
  {
    if (DEBUG_NONLINEAR)
    {
      if ( Teuchos::is_null(solverPtr_) )
        dout() << "Creating new NLS solver b/c it is 0." <<std::endl;
      if ((usemode_) && (lastParametersMode_ != parametersMode_))
        dout() << "Creating new NLS solver b/c starting next phase, post-DC." <<std::endl;
    }

    solverPtr_ = NOX::Solver::buildSolver(groupPtr_,
        paramsPtr->getStatusTests(),
        paramsPtr->getNoxParams());
  }
  else  // solverPtr is not null, and the mode didn't change since last call.
  {
    if (DEBUG_NONLINEAR)
      dout() << "NOT Creating new NLS solver, just resetting." <<std::endl;

    solverPtr_->reset(groupPtr_->getX());
  }

  // If .IC, .NODESET or .OP have been specified, then set up the
  // augmented linear systems
  if ((usemode_) && (mode_ != Nonlinear::TRANSIENT) && (mode_ != Nonlinear::HB_MODE) ) 
  {
    if (ICspecified_)
    {
      usedIC=icCont (paramsPtr);
    }
    else if (NODESETspecified_)
    {
      usedNODESET=nodesetCont0 (paramsPtr);
    }
  }

  NOX::StatusTest::StatusType status = solverPtr_->solve();
  firstSolveComplete_ = true;

  if (usedIC)
  {
    groupPtr_->setAugmentLinearSystem(false, Teuchos::null);
  }
  // Send back the correct return code
  if (DEBUG_NONLINEAR)
    dout() << "return code for transient params: " << transientParams_.getStatusTestReturnCode() <<std::endl
                 << "return code for hb params: " << hbParams_.getStatusTestReturnCode() << std::endl
                 << "return code for nlp params: " << nlpParams_.getStatusTestReturnCode() << std::endl
                 << "return code for dc params: " << dcParams_.getStatusTestReturnCode () << std::endl;

  return paramsPtr->getStatusTestReturnCode();
}

//-----------------------------------------------------------------------------
// Function      : Interface::naturalParameterContinuationSolve 
// Purpose       : 
// Special Notes : This corresponds to 
//                    .options nonlin continuation=nat or 
//                    .options nonlin continuation=1
// Scope         : public
// Creator       : Roger Pawlowski, Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::naturalParameterContinuationSolve ( ParameterSet* paramsPtr )
{
  std::vector<std::string> pars;
  std::string con;
  int j, numParam = 0;
  double value;
  bool usedOP=false;
  bool usedNODESET=false;
  bool usedIC=false;

  if ((usemode_) && (mode_ != Nonlinear::TRANSIENT))
  {
    if (ICspecified_)
    {
      usedIC=icCont (paramsPtr);
    }
    else if (NODESETspecified_)
    {
      usedNODESET=nodesetCont1 (paramsPtr);
    }
  }

  Teuchos::RCP<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();

  // Create Parameter Vector and get the stepper parameter list.
  LOCA::ParameterVector locaPVec;
  Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");
  Teuchos::ParameterList& stepSizeList = locaList->sublist("Step Size");

  while (paramsPtr->getVectorParam("CONPARAM", numParam, con))
  {
    pars.push_back(con);
    numParam++;
  }

  // Fail out if no parameters have been specified.
  if ( numParam == 0 ) {
    Report::UserFatal0() << "Using \"continuation=1\" requires a parameter to be set with the conparam keyword in the loca option block!";
  }

  // check if "skip df/dp" is true or false.  If false, then this is probably an arclength 
  // continuation solve, and we need to make sure sensitivity arrays are allocated
  const std::string strArg("Skip df/dp");
  bool skipDFDP = stepperList.get(strArg, false);
  if (!skipDFDP)
  {
    dsPtr_->allocateSensitivityArrays(numParam,true, false); // last arg is true only for transient adjoints, and 2nd arg is for direct
  }

  Teuchos::RCP<AugmentLinSys> als;

  std::vector<double> minValue(numParam, 0.);
  std::vector<double> maxValue(numParam, 0.);
  std::vector<double> initialValue(numParam, 0.);
  std::vector<double> initialStepSize(numParam, 0.);
  std::vector<double> minStepSize(numParam, 0.);
  std::vector<double> maxStepSize(numParam, 0.);
  std::vector<double> aggressiveness(numParam, 0.);

  // Check the size of params to make sure users provided all required data
  std::vector<std::string> paramNames(0);
  paramNames.push_back("MINVALUE");
  paramNames.push_back("MAXVALUE");
  paramNames.push_back("INITIALVALUE");
  paramNames.push_back("INITIALSTEPSIZE");
  paramNames.push_back("MINSTEPSIZE");
  paramNames.push_back("MAXSTEPSIZE");
  paramNames.push_back("AGGRESSIVENESS");

  for (std::size_t p = 0 ; p < paramNames.size() ; ++p) {
    if ( paramsPtr->getVectorParamSize(paramNames[p]) != numParam) {
      Report::UserFatal0() << "The parameter \"" 
                           << paramNames[p]  
                           << "\" must have a list of values with size equal to the numParamber of parameters specified in \"conparam\".";
    }
  }

  for (iParam_=0 ; iParam_<numParam ; ++iParam_)
  {
    paramsPtr->getVectorParam("MINVALUE", iParam_, value);
    minValue[iParam_] = value;

    paramsPtr->getVectorParam("MAXVALUE", iParam_, value);
    maxValue[iParam_] = value;

    paramsPtr->getVectorParam("INITIALVALUE", iParam_, value);
    initialValue[iParam_] = value;

    paramsPtr->getVectorParam("INITIALSTEPSIZE", iParam_, value);
    initialStepSize[iParam_] = value;

    paramsPtr->getVectorParam("MINSTEPSIZE", iParam_, value);
    minStepSize[iParam_] = value;

    paramsPtr->getVectorParam("MAXSTEPSIZE", iParam_, value);
    maxStepSize[iParam_] = value;

    paramsPtr->getVectorParam("AGGRESSIVENESS", iParam_, value);
    aggressiveness[iParam_] = value;

    locaPVec.addParameter (pars[iParam_], initialValue[iParam_]);

  }

  if (usedOP || usedNODESET)
  {
    const N_NLS_LOCA::Group & conLocaGrp =
      dynamic_cast<const N_NLS_LOCA::Group&>(solverPtr_->getSolutionGroup());
    *groupPtr_ = const_cast<N_NLS_LOCA::Group&>(conLocaGrp);
    solverPtr_ = Teuchos::null;
  }

  groupPtr_->setParams(locaPVec);

  LOCA::Abstract::Iterator::IteratorStatus locaStatus;

  for (iParam_=0 ; iParam_<numParam ; ++iParam_)
  {
    // Copy out the solution and use it in the next run
    if (iParam_ > 0)
    {
      groupPtr_->copy(*(stepperPtr_->getSolutionGroup()));
    }

    stepperList.set("Continuation Parameter", pars[iParam_]);
    stepperList.set("Initial Value", initialValue[iParam_]);
    stepperList.set("Max Value", maxValue[iParam_]);
    stepperList.set("Min Value", minValue[iParam_]);
    stepSizeList.set("Initial Step Size", initialStepSize[iParam_]);
    stepSizeList.set("Min Step Size", minStepSize[iParam_]);
    stepSizeList.set("Max Step Size", maxStepSize[iParam_]);
    stepSizeList.set("Aggressiveness", aggressiveness[iParam_]);

    for (j=0 ; j<iParam_ ; ++j)
    {
      locaPVec.setValue(pars[j], maxValue[j]);
    }
    for (j=iParam_ ; j<numParam ; ++j)
    {
      locaPVec.setValue(pars[j], initialValue[j]);
    }
    groupPtr_->setParams(locaPVec);

    if (iParam_==0)
    {
      if (!usedOP && !usedNODESET) // (usedOP and usedNODESET have already loaded F)
      {
        groupPtr_->computeF();
      }
    }

    // RPP 03/08/2006 If this is a special parameter used in voltage
    // node resistance (Spice's GMIN stepping) then we need to
    // intercept and handle the parameter in the solver (in the LOCA
    // Group) - the device package does nothing for this parameter.
    if (pars[iParam_] != "GSTEPPING")
    {
      // Do the continuation run
      resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
      locaStatus = stepperPtr_->run();
      if (usedIC || usedNODESET)
      {
        groupPtr_->setAugmentLinearSystem(false, Teuchos::null);
      }
    }
    else
    {
      if ((usemode_) && (mode_ != Nonlinear::TRANSIENT))
      {
        if (ICspecified_)
        {
          usedIC=icCont3 (paramsPtr);
        }
      }

      // Initialize parameters in xyce
      if (!usedOP && !usedNODESET) // (usedOP and usedNODESET have already loaded F)
      {
        groupPtr_->computeF();
      }

      if (!usedIC)
      {
        Teuchos::RCP<AugmentLinSys> als =
          paramsPtr->createAugmentLinearSystem(lasSysPtr_);
        groupPtr_->setAugmentLinearSystem(true, als);
      }

      // Do the continuation run
      resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
      locaStatus = stepperPtr_->run();

      groupPtr_->setAugmentLinearSystem(false, Teuchos::null);
    }

    // Kick out if continuation failed
    if (locaStatus != LOCA::Abstract::Iterator::Finished)
      return (-1);

    // Increment Param Number Tracking
    isFirstContinuationParam_ = false;
    firstSolveComplete_ = true;
  }
  return (paramsPtr->getStatusTestReturnCode());
}

//-----------------------------------------------------------------------------
// Function      : Interface::mosfetContinuationSolve 
// Purpose       : 
// Special Notes : This corresponds to 
//                    .options nonlin continuation=mos or 
//                    .options nonlin continuation=2
// Scope         : public
// Creator       : Roger Pawlowski, Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::mosfetContinuationSolve ( ParameterSet* paramsPtr )
{
  bool usedOP=false;
  bool usedNODESET=false;
  bool usedIC=false;

  if ((usemode_) && (mode_ != Nonlinear::TRANSIENT))
  {
    if (ICspecified_)
    {
      usedIC=icCont (paramsPtr);
    }
    else if (NODESETspecified_)
    {
      usedNODESET=nodesetCont1 (paramsPtr);
    }
  }

  Teuchos::RCP<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();

  // Create the continuation parameter names
  std::string gain = "mosfet:gainscale";
  std::string nonlinear = "mosfet:nltermscale";

  // Create Parameter Vector and get the stepper parameter list.
  LOCA::ParameterVector locaPVec;
  Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");
  Teuchos::ParameterList& predictorList = locaList->sublist("Predictor");
  Teuchos::ParameterList& stepSizeList = locaList->sublist("Step Size");

  // Continuation solve from alpha2 (gain) = 0.0 -> 1.0
  //                    with alpha1 (nlterm) = 0.0 (constant)
  locaPVec.addParameter(gain, 0.0);
  locaPVec.addParameter(nonlinear, 0.0);
  groupPtr_->setParams(locaPVec);
  stepperList.set("Continuation Parameter", gain);

  stepperList.set("Initial Value", 0.0);
  stepperList.set("Max Value", 1.0);
  stepperList.set("Min Value",-1.0);

  stepSizeList.set("Initial Step Size", 0.2);
  stepSizeList.set("Min Step Size", 1.0e-4);
  stepSizeList.set("Max Step Size", 1.0);
  stepSizeList.set("Aggressiveness", 1.0);

  // assert the user-specified defaults, if any.
  dcParams_.applySavedLocaOptions();

  // Initialize parameters in xyce
  if (!usedOP && !usedNODESET) // (usedOP and usedNODESET have already loaded F)
  {
    groupPtr_->computeF();
  }

  // Do the continuation run
  iParam_ = 0;
  resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
  LOCA::Abstract::Iterator::IteratorStatus locaStatus = stepperPtr_->run();

  // Kick out if continuation failed
  if (locaStatus != LOCA::Abstract::Iterator::Finished)
    return (-1);

  // Increment Param Number Tracking
  isFirstContinuationParam_ = false;
  firstSolveComplete_ = true;

  // Copy out the solution and use it in the next run
  groupPtr_->copy(*(stepperPtr_->getSolutionGroup()));

  // Continuation solve from alpha1 (nlterm) = 0.0 -> 1.0
  //                    with alpha2 (gain) = 1.0 (constant)
  stepperList.set("Continuation Parameter", nonlinear);

  stepperList.set("Initial Value", 0.0);
  stepperList.set("Max Value", 1.0);
  stepperList.set("Min Value",-1.0);

  stepSizeList.set("Initial Step Size", 0.2);
  stepSizeList.set("Min Step Size", 1.0e-4);
  stepSizeList.set("Max Step Size", 1.0);
  stepSizeList.set("Aggressiveness", 1.0);

  // assert the user-specified defaults, if any.
  dcParams_.applySavedLocaOptions();

  locaPVec.setValue(gain, 1.0);
  locaPVec.setValue(nonlinear, 0.0);
  groupPtr_->setParams(locaPVec);

  // Do the continuation run
  iParam_ = 1;
  resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
  locaStatus = stepperPtr_->run();

  if (usedIC)
  {
    groupPtr_->setAugmentLinearSystem(false, Teuchos::null);
  }

  // Return the solution status
  if (locaStatus != LOCA::Abstract::Iterator::Finished)
    return (-1);
  else
    return (paramsPtr->getStatusTestReturnCode());
}

//-----------------------------------------------------------------------------
// Function      : Interface::mosfetContinuationSolve2 
// Purpose       : 
// Return Type   : 
// Special Notes : This corresponds to 
//                    .options nonlin continuation=newmos or 
//                    .options nonlin continuation=4
// Scope         : public
// Creator       : Roger Pawlowski, Eric Keiter
// Creation Date : 2/21/2004
//-----------------------------------------------------------------------------
int Interface::mosfetContinuationSolve2 ( ParameterSet* paramsPtr )
{
  Teuchos::RCP<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();

  // Create the continuation parameter names
  std::string gain = "mosfet:gainscale";
  std::string nonlinear = "mosfet:nltermscale";
  std::string size = "mosfet:sizescale";

  // Create Parameter Vector and get the stepper parameter list.
  LOCA::ParameterVector locaPVec;
  Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");

  // Continuation solve from alpha2 (gain) = 0.0 -> 1.0
  //                    with alpha1 (nlterm) = 0.0 (constant)
  //                    sizeScale = 0.0 (constant)
  locaPVec.addParameter(gain, 0.0);
  locaPVec.addParameter(nonlinear, 0.0);
  locaPVec.addParameter(size, 0.0);
  groupPtr_->setParams(locaPVec);
  stepperList.set("Continuation Parameter", gain);
  stepperList.set("Initial Value", 0.0);
  stepperList.set("Max Value", 1.0);

  // Initialize parameters in xyce
  groupPtr_->computeF();

  // Do the continuation run
  resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
  LOCA::Abstract::Iterator::IteratorStatus locaStatus = stepperPtr_->run();

  // Kick out if continuation failed
  if (locaStatus != LOCA::Abstract::Iterator::Finished)
    return (-1);

  // Increment Param Number Tracking
  isFirstContinuationParam_ = false;
  firstSolveComplete_ = true;

  // Copy out the solution and use it in the next run
  groupPtr_->copy(*(stepperPtr_->getSolutionGroup()));

  // Continuation solve from alpha1 (nlterm) = 0.0 -> 1.0
  //                    with alpha2 (gain) = 1.0 (constant)
  //                    with size   (scale) = 0.0 (constant)
  stepperList.set("Continuation Parameter", nonlinear);
  stepperList.set("Initial Value", 0.0);
  stepperList.set("Max Value", 1.0);
  locaPVec.setValue(gain, 1.0);
  locaPVec.setValue(nonlinear, 0.0);
  locaPVec.setValue(size, 0.0);
  groupPtr_->setParams(locaPVec);

  // Do the continuation run
  resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
  locaStatus = stepperPtr_->run();

  // Kick out if continuation failed
  if (locaStatus != LOCA::Abstract::Iterator::Finished)
    return (-1);

  // Copy out the solution and use it in the next run
  groupPtr_->copy(*(stepperPtr_->getSolutionGroup()));

  // Continuation solve from alpha1 (nlterm) = 1.0 -> 1.0
  //                    with alpha2 (gain) = 1.0 (constant)
  //                    with size   (scale) = 0.0 -> 1.0
  stepperList.set("Continuation Parameter", size);
  stepperList.set("Initial Value", 0.0);
  stepperList.set("Max Value", 1.0);
  locaPVec.setValue(gain, 1.0);
  locaPVec.setValue(nonlinear, 1.0);
  locaPVec.setValue(size, 0.0);
  groupPtr_->setParams(locaPVec);

  // Do the continuation run
  resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
  locaStatus = stepperPtr_->run();

  // Return the solution status
  if (locaStatus != LOCA::Abstract::Iterator::Finished)
    return (-1);
  else
    return (paramsPtr->getStatusTestReturnCode());
}

//-----------------------------------------------------------------------------
// Function      : Interface::mosfetContinuationSolve3 
// Purpose       : (5) Mosfet:BSIM3:Inverter specific continuation 
// Special Notes : This corresponds to 
//                    .options nonlin continuation=bsim3inv1 or 
//                    .options nonlin continuation=5
// Scope         : public
// Creator       : Roger Pawlowski, Eric Keiter
// Creation Date : 2/25/2004
//-----------------------------------------------------------------------------
int Interface::mosfetContinuationSolve3 ( ParameterSet* paramsPtr )
{
  Teuchos::RCP<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();

  // Create the continuation parameter names
  //std::string gain = "mosfet:gainscale";
  std::string nonlinear = "mosfet:nltermscale";
  std::string size = "mosfet:sizescale";

  // Create Parameter Vector and get the stepper parameter list.
  LOCA::ParameterVector locaPVec;
  Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");

  // Continuation solve from alpha2 (gain) = 0.0 -> 1.0
  //                    with alpha1 (nlterm) = 0.0 (constant)
  //                    sizeScale = 0.0 (constant)
  //locaPVec.addParameter(gain, 0.0);
  locaPVec.addParameter(nonlinear, 0.0);
  locaPVec.addParameter(size, 0.0);
  groupPtr_->setParams(locaPVec);
  stepperList.set("Continuation Parameter", nonlinear);
  stepperList.set("Initial Value", 0.0);
  stepperList.set("Max Value", 1.0);

  // Initialize parameters in xyce
  groupPtr_->computeF();

  // Do the continuation run
  resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
  LOCA::Abstract::Iterator::IteratorStatus locaStatus = stepperPtr_->run();

  // Kick out if continuation failed
  if (locaStatus != LOCA::Abstract::Iterator::Finished)
    return (-1);

  // Increment Param Number Tracking
  isFirstContinuationParam_ = false;
  firstSolveComplete_ = true;

  // Copy out the solution and use it in the next run
  groupPtr_->copy(*(stepperPtr_->getSolutionGroup()));

  // Continuation solve from alpha1 (nlterm) = 0.0 -> 1.0
  //                    with alpha2 (gain) = 1.0 (constant)
  //                    with size   (scale) = 0.0 (constant)
  stepperList.set("Continuation Parameter", size);
  stepperList.set("Initial Value", 0.0);
  stepperList.set("Max Value", 1.0);
  //locaPVec.setValue(gain, 1.0);
  locaPVec.setValue(nonlinear, 1.0);
  locaPVec.setValue(size, 0.0);
  groupPtr_->setParams(locaPVec);

  // Do the continuation run
  resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
  locaStatus = stepperPtr_->run();

  // Return the solution status
  if (locaStatus != LOCA::Abstract::Iterator::Finished)
    return (-1);
  else
    return (paramsPtr->getStatusTestReturnCode());
}

//-----------------------------------------------------------------------------
// Function      : Interface::mosfetContinuationSolve4 
// Purpose       : (6) Mosfet:BSIM3:Inverter specific continuation 
// Special Notes : This corresponds to 
//                    .options nonlin continuation=bsim3inv2 or 
//                    .options nonlin continuation=6
// Scope         : public
// Creator       : Roger Pawlowski, Eric Keiter
// Creation Date : 2/25/2004
//-----------------------------------------------------------------------------
int Interface::mosfetContinuationSolve4 ( ParameterSet* paramsPtr )
{
  Teuchos::RCP<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();

  // Create the continuation parameter names
  std::string gain = "mosfet:gainscale";
  std::string nonlinear = "mosfet:nltermscale";
  std::string size = "mosfet:sizescale";

  // Create Parameter Vector and get the stepper parameter list.
  LOCA::ParameterVector locaPVec;
  Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");

  // Continuation solve
  locaPVec.addParameter(gain, 1.0);
  locaPVec.addParameter(nonlinear, 1.0);
  locaPVec.addParameter(size, 0.0);
  groupPtr_->setParams(locaPVec);
  stepperList.set("Continuation Parameter", size);
  stepperList.set("Initial Value", 0.0);
  stepperList.set("Max Value", 1.0);

  // Initialize parameters in xyce
  groupPtr_->computeF();

  // Do the continuation run
  resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
  LOCA::Abstract::Iterator::IteratorStatus locaStatus = stepperPtr_->run();

  // Return the solution status
  if (locaStatus != LOCA::Abstract::Iterator::Finished)
    return (-1);
  else
  {
    firstSolveComplete_ = true;
    return (paramsPtr->getStatusTestReturnCode());
  }
}

//-----------------------------------------------------------------------------
// Function      : Interface::mosfetContinuationSolve5 
// Purpose       : 
// Return Type   : 
// Special Notes : This corresponds to 
//                    .options nonlin continuation=test or 
//                    .options nonlin continuation=8
// Scope         : public
// Creator       : Roger Pawlowski, Eric Keiter
// Creation Date : 2/25/2004
//-----------------------------------------------------------------------------
int Interface::mosfetContinuationSolve5 ( ParameterSet* paramsPtr )
{
  // Get some initial objects
  Teuchos::RCP<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();
  Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");
  LOCA::ParameterVector locaPVec;

  // Create storage for continuation objects
  int numHomotopyContinuationRuns = 1;
  std::vector<std::string> names(numHomotopyContinuationRuns);
  std::vector<double> initialVal(numHomotopyContinuationRuns);
  std::vector<double> finalVal(numHomotopyContinuationRuns);
  std::vector<double> minVal(numHomotopyContinuationRuns);
  std::vector<double> maxVal(numHomotopyContinuationRuns);

  // Set up continuation steps
  // ***************************************************
  // Changes vals below for continuation
  // ***************************************************
  names[0] = "mosfet:gainscale";
  //names[1] = "mosfet:nltermscale";
  //names[2] = "mosfet:sizescale";
  //names[2] = "mosfet:l";
  initialVal[0] = 0.0;
  //initialVal[1] = 0.0;
  //initialVal[2] = 0.0;
  finalVal[0] = 1.0;
  //finalVal[1] = 1.0;
  //finalVal[2] = 1.0;
  // ***************************************************
  std::string n1 = "mosfet:nltermscale";
  locaPVec.addParameter(n1, 0.0);
  //std::string n2 = "mosfet:sizescale";
  //locaPVec.addParameter(n2, 0.0);
  // ***************************************************
  // ***************************************************

  // Ste up max/min bounds
  for (int i = 0; i < names.size(); ++i) {
    if (finalVal[i] > initialVal[i]) {
      minVal[i] = initialVal[i];
      maxVal[i] = finalVal[i];
    }
    else {
      minVal[i] = finalVal[i];
      maxVal[i] = initialVal[i];
    }
  }

  // Initialize loca parameter vector
  for (int i = 0; i < names.size(); ++i)
    locaPVec.addParameter(names[i], initialVal[i]);

  LOCA::Abstract::Iterator::IteratorStatus locaStatus;

  LOCA::StatusTest::Wrapper test(paramsPtr->getStatusTests());

  // Loop over the number of homotopy steps
  for (int hs = 0; hs < names.size(); ++hs) {
    for (int i = 0; i < names.size(); ++i) {
      if (i >= hs)
        locaPVec.setValue(names[i], initialVal[i]);
      else
        locaPVec.setValue(names[i], finalVal[i]);
    }
    groupPtr_->setParams(locaPVec);
    stepperList.set("Continuation Parameter", names[hs]);
    stepperList.set("Initial Value", initialVal[hs]);
    stepperList.set("Min Value", minVal[hs]);
    stepperList.set("Max Value", maxVal[hs]);

    // Initialize parameters in xyce
    groupPtr_->computeF();

    // Do the continuation run
    resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
    locaStatus = stepperPtr_->run();

    // Kick out if continuation failed
    if (locaStatus != LOCA::Abstract::Iterator::Finished)
      return (-1);

    // Increment Param Number Tracking
    isFirstContinuationParam_ = false;
    firstSolveComplete_ = true;

    // Copy out the solution and use it in the next run
    groupPtr_->copy(*(stepperPtr_->getSolutionGroup()));
  }

  // Return converged solver code
  return (paramsPtr->getStatusTestReturnCode());
}

//-----------------------------------------------------------------------------
// Function      : Interface::mosfetContinuationSolve6 
// Purpose       : 
// Special Notes : This corresponds to 
//                    .options nonlin continuation=power or 
//                    .options nonlin continuation=10
// Scope         : public
// Creator       : Roger Pawlowski, Eric Keiter
// Creation Date : 2/25/2004
//-----------------------------------------------------------------------------
int Interface::mosfetContinuationSolve6 ( ParameterSet* paramsPtr )
{
  // Get some initial objects
  Teuchos::RCP<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();
  Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");
  LOCA::ParameterVector locaPVec;

  // Create storage for continuation objects
  int numHomotopyContinuationRuns = 4;
  std::vector<std::string> names(numHomotopyContinuationRuns);
  std::vector<double> initialVal(numHomotopyContinuationRuns);
  std::vector<double> finalVal(numHomotopyContinuationRuns);
  std::vector<double> minVal(numHomotopyContinuationRuns);
  std::vector<double> maxVal(numHomotopyContinuationRuns);

  // Set up continuation steps
  // ***************************************************
  // Changes vals below for continuation
  // ***************************************************
  names[0] = stepperList.get("Power Node", "VA:V0");
  names[1] = "mosfet:gainscale";
  names[2] = "mosfet:nltermscale";
  names[3] = "mosfet:sizescale";
  initialVal[0] = 0.0;
  initialVal[1] = 0.0;
  initialVal[2] = 0.0;
  initialVal[3] = 0.0;
  finalVal[0] = 1.0;
  finalVal[1] = 1.0;
  finalVal[2] = 1.0;
  finalVal[3] = 1.0;
  // ***************************************************
  //std::string n1 = "mosfet:nltermscale";
  //locaPVec.addParameter(n1, 0.0);
  //std::string n2 = "mosfet:sizescale";
  //locaPVec.addParameter(n2, 0.0);
  // ***************************************************
  // ***************************************************

  // Set up max/min bounds
  for (int i = 0; i < names.size(); ++i) {
    if (finalVal[i] > initialVal[i]) {
      minVal[i] = initialVal[i];
      maxVal[i] = finalVal[i];
    }
    else {
      minVal[i] = finalVal[i];
      maxVal[i] = initialVal[i];
    }
  }

  // Initialize loca parameter vector
  for (int i = 0; i < names.size(); ++i)
    locaPVec.addParameter(names[i], initialVal[i]);

  LOCA::Abstract::Iterator::IteratorStatus locaStatus;

  // Loop over the number of homotopy steps
  for (int hs = 0; hs < names.size(); ++hs) {
    for (int i = 0; i < names.size(); ++i) {
      if (i >= hs)
        locaPVec.setValue(names[i], initialVal[i]);
      else
        locaPVec.setValue(names[i], finalVal[i]);
    }
    groupPtr_->setParams(locaPVec);
    stepperList.set("Continuation Parameter", names[hs]);
    stepperList.set("Initial Value", initialVal[hs]);
    stepperList.set("Min Value", minVal[hs]);
    stepperList.set("Max Value", maxVal[hs]);

    // Initialize parameters in xyce
    groupPtr_->computeF();

    // Do the continuation run
    resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
    locaStatus = stepperPtr_->run();

    // Kick out if continuation failed
    if (locaStatus != LOCA::Abstract::Iterator::Finished)
      return (-1);

    firstSolveComplete_ = true;

    // Copy out the solution and use it in the next run
    groupPtr_->copy(*(stepperPtr_->getSolutionGroup()));
  }

  // Return converged solver code
  return (paramsPtr->getStatusTestReturnCode());
}

//-----------------------------------------------------------------------------
// Function      : Interface::blockGainscaleMosfetSolve 
// Purpose       : 
// Special Notes : This corresponds to 
//                    .options nonlin continuation=blockgain or 
//                    .options nonlin continuation=7
// Scope         : public
// Creator       : Roger Pawlowski, Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::blockGainscaleMosfetSolve ( ParameterSet* paramsPtr )
{
  // Get some initial objects
  Teuchos::RCP<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();
  Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");
  LOCA::ParameterVector locaPVec;

  // Create storage for continuation objects
  int numGainBlocks = nonlinearEquationLoader_->getHomotopyBlockSize();
  int numHomotopyContinuationRuns = 1 + numGainBlocks;
  std::vector<std::string> names(numHomotopyContinuationRuns);
  std::vector<double> initialVal(numHomotopyContinuationRuns);
  std::vector<double> finalVal(numHomotopyContinuationRuns);
  std::vector<double> minVal(numHomotopyContinuationRuns);
  std::vector<double> maxVal(numHomotopyContinuationRuns);

  // Set up continuation steps
  // ***************************************************
  // Changes vals below for continuation
  // ***************************************************
  for (int i = 0; i < numGainBlocks; ++i) {
    std::stringstream s;
    s << i;
    names[i] = "mosfet:gainscale_block_" + s.str();
    initialVal[i] = 0.0;
    finalVal[i] = 1.0;
  }
  names[numHomotopyContinuationRuns - 1] = "mosfet:nltermscale";
  initialVal[numHomotopyContinuationRuns - 1] = 0.0;
  finalVal[numHomotopyContinuationRuns - 1] = 1.0;
  // ***************************************************
  // ***************************************************

  // Ste up max/min bounds
  for (int i = 0; i < names.size(); ++i) {
    if (finalVal[i] > initialVal[i]) {
      minVal[i] = initialVal[i];
      maxVal[i] = finalVal[i];
    }
    else {
      minVal[i] = finalVal[i];
      maxVal[i] = initialVal[i];
    }
  }

  // Initialize loca parameter vector
  for (int i = 0; i < names.size(); ++i)
    locaPVec.addParameter(names[i], initialVal[i]);

  LOCA::Abstract::Iterator::IteratorStatus locaStatus;

  // Loop over the number of homotopy steps
  for (int hs = 0; hs < names.size(); ++hs) {
    for (int i = 0; i < names.size(); ++i) {
      if (i >= hs)
        locaPVec.setValue(names[i], initialVal[i]);
      else
        locaPVec.setValue(names[i], finalVal[i]);
    }
    groupPtr_->setParams(locaPVec);
    stepperList.set("Continuation Parameter", names[hs]);
    stepperList.set("Initial Value", initialVal[hs]);
    stepperList.set("Min Value", minVal[hs]);
    stepperList.set("Max Value", maxVal[hs]);

    // Initialize parameters in xyce
    groupPtr_->computeF();

    // Do the continuation run
    resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
    locaStatus = stepperPtr_->run();

    // Kick out if continuation failed
    if (locaStatus != LOCA::Abstract::Iterator::Finished)
      return (-1);

    // Increment Param Number Tracking
    isFirstContinuationParam_ = false;
    firstSolveComplete_ = true;

    // Copy out the solution and use it in the next run
    groupPtr_->copy(*(stepperPtr_->getSolutionGroup()));
  }

  // Return converged solver code
  return (paramsPtr->getStatusTestReturnCode());
}

//-----------------------------------------------------------------------------
// Function      : Interface::gminSteppingSolve
// Purpose       : 
// Special Notes : This corresponds to 
//                    .options nonlin continuation=gmin or 
//                    .options nonlin continuation=3
// Scope         : public
// Creator       : Roger Pawlowski, Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::gminSteppingSolve ( ParameterSet* paramsPtr )
{
  Teuchos::RCP<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();

  // Create the continuation parameter names
  std::string gmin = "GSTEPPING";

  // Create Parameter Vector and get the stepper parameter list.
  LOCA::ParameterVector locaPVec;
  Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");
  Teuchos::ParameterList& predictorList = locaList->sublist("Predictor");
  Teuchos::ParameterList& stepSizeList = locaList->sublist("Step Size");

  // Continuation solve using Gmin stepping.
  locaPVec.addParameter(gmin, 0.0);
  groupPtr_->setParams(locaPVec);
  stepperList.set("Continuation Parameter", gmin);
  stepperList.set("Continuation Method", "Natural");

  stepSizeList.set("Method", "Adaptive");
  predictorList.set("Method", "Constant");

  stepperList.set("Initial Value", 4.0);
  stepperList.set("Min Value", -4.0);
  paramsPtr->set_gstepping_min_value (-4.0);
  stepperList.set("Max Value", 4.0);

  stepSizeList.set("Initial Step Size", -2.0);
  stepSizeList.set("Min Step Size", 1.0e-6);
  stepSizeList.set("Max Step Size", 1.0e+12);
  stepSizeList.set("Aggressiveness", 0.01);

  stepperList.set("Max Steps", 400);
  stepperList.set("Max Nonlinear Iterations", 20);

  // the following set of if-statemtents call functions that
  // allocate augmented linear systems for various scenarios.
  // It is important that the augmented systems get allocated
  // after the paramter (above) have been set.
  bool usedOP=false;
  bool usedNODESET=false;
  bool usedIC=false;

  if ((usemode_) && (mode_ != Nonlinear::TRANSIENT))
  {
    if (ICspecified_)
    {
      usedIC=icCont3 (paramsPtr);
    }
    else if (NODESETspecified_)
    {
      usedNODESET=nodesetCont1 (paramsPtr);
    }
  }

  // Initialize parameters in xyce
  if (!usedOP && !usedNODESET) // (usedOP and usedNODESET have already loaded F)
  {
    groupPtr_->computeF();
  }

  // Do the continuation run
  iParam_ = 0;

  if (!usedIC)
  {
    Teuchos::RCP<AugmentLinSys> als =
      paramsPtr->createAugmentLinearSystem(lasSysPtr_);
    groupPtr_->setAugmentLinearSystem(true, als);
  }

  // Do the continuation run
  resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
  LOCA::Abstract::Iterator::IteratorStatus locaStatus = stepperPtr_->run();

  groupPtr_->setAugmentLinearSystem(false, Teuchos::null);

  // Kick out if continuation failed
  if (locaStatus != LOCA::Abstract::Iterator::Finished)
    return (-1);
  else
    return (paramsPtr->getStatusTestReturnCode());
}

//-----------------------------------------------------------------------------
// Function      : Interface::pseudoTransientSolve
// Purpose       : 
// Special Notes : This corresponds to 
//                    .options nonlin continuation=pseudo or 
//                    .options nonlin continuation=9
// Scope         : 
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::pseudoTransientSolve ( ParameterSet* paramsPtr )
{
  Teuchos::RCP<NOX::StatusTest::Combo> ctest =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<Teuchos::ParameterList> locaList =
    paramsPtr->getLocaParams();
  Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");
  Teuchos::ParameterList& stepSizeList = locaList->sublist("Step Size");
  double initialStepSize =
    stepSizeList.get("Initial Step Size", 1.0e-3);
  double minStepSize = stepSizeList.get("Min Step Size", 1.0e-12);
  double maxStepSize = stepSizeList.get("Max Step Size", 1.0e4);
  Teuchos::RCP<Teuchos::ParameterList> noxList =
    paramsPtr->getNoxParams();

  // Create Pseudo Transient status tests.
  //paramsPtr->getStatusTests()
  Teuchos::RCP<NOX::StatusTest::MaxIters> mi =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(stepperList.get("Max Steps", 200)));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<PseudoTransientTest> pt =
    Teuchos::rcp(new PseudoTransientTest(maxStepSize, 1.0e-8));

  ctest->addStatusTest(mi);
  ctest->addStatusTest(fv);
  ctest->addStatusTest(pt);

  // First solve - pseudo transient

  Teuchos::RCP<AugmentLinSys> als =
    paramsPtr->createAugmentLinearSystem(lasSysPtr_);

  solverPtr_ = Teuchos::rcp(new PseudoTransientBased(als,
        groupPtr_,
        ctest,
        noxList,
        initialStepSize,
        minStepSize,
        maxStepSize));

  NOX::StatusTest::StatusType status = solverPtr_->solve();
  firstSolveComplete_ = true;

  // RPP 3/7/2006: We don't care if pseudo transient solve fails the
  // solve() call above.  This is just to get an inital guess for
  // the corrector step in the next solve.  So we don't check the
  // status at this point.

  // Turn off pseudo transient in groups. (this is also done in the
  // pseudo transient solver, but just being safe - groups could be
  // different).
  groupPtr_->setAugmentLinearSystem(false, Teuchos::null);

  // Second solve is the correct steady state solve
  solverPtr_ = NOX::Solver::buildSolver(groupPtr_,
      paramsPtr->getStatusTests(),
      paramsPtr->getNoxParams());
  status = solverPtr_->solve();

  // Send back the correct return code
  return (paramsPtr->getStatusTestReturnCode());
}

//-----------------------------------------------------------------------------
// Function      : Interface::artificialParameterHomotopy 
// Purpose       : 
// Special Notes : This corresponds to 
//                    .options nonlin continuation=art or 
//                    .options nonlin continuation=33
// Scope         : 
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::artificialParameterHomotopy ( ParameterSet* paramsPtr )
{
#ifdef Xyce_NOX_LOCA_ARTIFICIAL_HOMOTOPY_SUPPORT
  Teuchos::RCP<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();
  Teuchos::ParameterList& locaUtilsList = locaList->sublist("Utilities");

  Teuchos::RCP<LOCA::Homotopy::Group> hGrp =
    Teuchos::rcp(new LOCA::Homotopy::Group(*locaList, globalDataPtr_, groupPtr_, 1.0, 0.0));

  hGrp->computeF();

  locaList->sublist("Predictor").set("Secant", 0.999);
  locaList->sublist("Stepper").set("Max Value", 0.999);

  resetStepper(globalDataPtr_, hGrp, locaStatusTestPtr_, paramsPtr->getAllParams());

  LOCA::Abstract::Iterator::IteratorStatus locaStatus = stepperPtr_->run();
  firstSolveComplete_ = true;

  Teuchos::RCP<LOCA::Homotopy::Group> hGrp2 =
    Teuchos::rcp(new LOCA::Homotopy::Group(*locaList, globalDataPtr_, groupPtr_, 0.1, 1.0));

  locaList->sublist("Predictor").set("Secant", 0.999);
  locaList->sublist("Stepper").set("Initial Value", 0.999);
  locaList->sublist("Stepper").set("Max Value", 1.0);
  locaList->sublist("Step Size").set("Method", "Constant");
  locaList->sublist("Step Size").set("Initial Step Size", 0.0001);
  locaList->sublist("Step Size").set("Min Step Size", 0.0001);

  resetStepper(globalDataPtr_, hGrp2, locaStatusTestPtr_, paramsPtr->getAllParams());

  locaStatus = stepperPtr_->run();

  // Return the solution status
  if (locaStatus != LOCA::Abstract::Iterator::Finished)
    return (-1);
  else
    return (paramsPtr->getStatusTestReturnCode());

#else
  Report::UserFatal0() << "Nonlinear Solver (NOX::Interface) Artificial parameter continuation requires "
    << "building xyce with the define: -DXyce_NOX_LOCA_ARTIFICIAL_HOMOTOPY_SUPPORT to "
    << "allow LOCA to augment the diagonal of Jacobian! Either rebuild Xyce or do not "
    << "run Xyce with \"continuation=33\"";
  return -1;
#endif
}

//-----------------------------------------------------------------------------
// Function      : Interface::sourceSteppingSolve
// Purpose       : 
// Special Notes : This corresponds to 
//                    .options nonlin continuation=sourcestep or  
//                    .options nonlin continuation=34
// Scope         : public
// Creator       : Tom Russo
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::sourceSteppingSolve ( ParameterSet* paramsPtr )
{
  Teuchos::RCP<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();

  // Create the continuation parameter names
  std::string vsrcscale = "VSRCSCALE";

  // Create Parameter Vector and get the stepper parameter list.
  LOCA::ParameterVector locaPVec;
  Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");
  Teuchos::ParameterList& predictorList = locaList->sublist("Predictor");
  Teuchos::ParameterList& stepSizeList = locaList->sublist("Step Size");

  // Continuation solve using source stepping.
  locaPVec.addParameter(vsrcscale, 0.0);
  groupPtr_->setParams(locaPVec);
  stepperList.set("Continuation Parameter", vsrcscale);
  stepperList.set("Continuation Method", "Natural");

  stepSizeList.set("Method", "Adaptive");
  predictorList.set("Method", "Constant");

  stepperList.set("Initial Value", 0.0);
  stepperList.set("Min Value", -1.0);
  stepperList.set("Max Value", 1.0);

  stepSizeList.set("Initial Step Size", 0.2);
  stepSizeList.set("Min Step Size", 1.0e-4);
  stepSizeList.set("Max Step Size", 0.2);
  stepSizeList.set("Aggressiveness", 1.0);

  stepperList.set("Max Steps", 400);
  stepperList.set("Max Nonlinear Iterations", 20);

  // the following set of if-statemtents call functions that
  // allocate augmented linear systems for various scenarios.
  // It is important that the augmented systems get allocated
  // after the paramter (above) have been set.
  bool usedOP=false;
  bool usedNODESET=false;
  bool usedIC=false;

  if ((usemode_) && (mode_ != Nonlinear::TRANSIENT))
  {
    if (ICspecified_)
    {
      usedIC=icCont3 (paramsPtr);
    }
    else if (NODESETspecified_)
    {
      usedNODESET=nodesetCont1 (paramsPtr);
    }
  }

  // Initialize parameters in xyce
  if (!usedOP && !usedNODESET) // (usedOP and usedNODESET have already loaded F)
  {
    groupPtr_->computeF();
  }
  nonlinearEquationLoader_->resetScaledParams();

  // Do the continuation run
  resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
  LOCA::Abstract::Iterator::IteratorStatus locaStatus = stepperPtr_->run();

  groupPtr_->setAugmentLinearSystem(false, Teuchos::null);

  nonlinearEquationLoader_->resetScaledParams();

  // Kick out if continuation failed
  if (locaStatus != LOCA::Abstract::Iterator::Finished)
    return (-1);
  else
    return (paramsPtr->getStatusTestReturnCode());
}

//-----------------------------------------------------------------------------
// Function      : Interface::solve
//
// Purpose       : Reset all the counters and parameters and solve the
//                 nonlinear problem for this time step. The solution is
//                 stored in nextSolVector (obtained from the N_LAS_System
//                 registered above by registerLinearSystem).
//
// Special Notes : Should not be called until *after* initializeAll() has
//                 been called.
//
// Return Type   : Integer - postive for sucess, negative for failure.
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::solve (Nonlinear::NonLinearSolver * nlsTmpPtr)
{

  try
  {
    // For base object
    Nonlinear::NonLinearSolver::resetCountersAndTimers_();

    if (DEBUG_NONLINEAR)
      Nonlinear::NonLinearSolver::setDebugFlags(getAnalysisManager().getStepNumber() + 1, getAnalysisManager().getTime());

    // Setup the status tests
    if (Teuchos::is_null(locaStatusTestPtr_))
    {
      locaDCOpStatusTestPtr_ =
        Teuchos::rcp(new LOCA::StatusTest::Wrapper(dcParams_.getStatusTests()));
      locaTransientStatusTestPtr_ =
        Teuchos::rcp(new LOCA::StatusTest::Wrapper(transientParams_.getStatusTests()));
      locaHBStatusTestPtr_ =
        Teuchos::rcp(new LOCA::StatusTest::Wrapper(hbParams_.getStatusTests()));
      locaDC_NLPStatusTestPtr_ = 
        Teuchos::rcp(new LOCA::StatusTest::Wrapper(nlpParams_.getStatusTests()));
    }

    // Pick the parameter set to use.
    ParameterSet* paramsPtr;
    if ((usemode_) && (mode_ == Nonlinear::TRANSIENT))
    {
      paramsPtr = &transientParams_;
      locaStatusTestPtr_ = locaTransientStatusTestPtr_;
      lastParametersMode_ = parametersMode_;
      parametersMode_ = Nonlinear::TRANSIENT;
    }
    else if ((usemode_) && (mode_ == Nonlinear::HB_MODE))
    {
      paramsPtr = &hbParams_;
      locaStatusTestPtr_ = locaHBStatusTestPtr_;
      lastParametersMode_ = parametersMode_;
      parametersMode_ = Nonlinear::HB_MODE;
    }
    else if ((usemode_) && (mode_ == Nonlinear::DC_NLPOISSON))
    {
      paramsPtr = &nlpParams_;
      locaStatusTestPtr_ = locaDC_NLPStatusTestPtr_;
      lastParametersMode_ = parametersMode_;
      parametersMode_ = Nonlinear::DC_NLPOISSON;
    }
    else
    {
      paramsPtr = &dcParams_;
      locaStatusTestPtr_ = locaDCOpStatusTestPtr_;
      lastParametersMode_ = parametersMode_;
      parametersMode_ = Nonlinear::DC_OP;
    }

    if (Teuchos::is_null(globalDataPtr_))
    {
      globalDataPtr_ = LOCA::createGlobalData(paramsPtr->getAllParams());
    }

    // set the xyce return codes:
    paramsPtr->setStatusTestReturnCodes(retCodes_);

    // Set up the shared system (we have to redo this every time because
    // the object pointed to by nextSolVectorPtrPtr may have changed.
    //delete sharedSystemPtr_;
    if (sharedSystemPtr_ == 0)
    {
      sharedSystemPtr_ = new SharedSystem(*dsPtr_->nextSolutionPtr,
                                          *rhsVectorPtr_,
                                          *jacobianMatrixPtr_,
                                          *NewtonVectorPtr_,
                                          *gradVectorPtr_,
                                          *lasSysPtr_,
                                          *this);
    }
    else
    {
      sharedSystemPtr_->reset(*dsPtr_->nextSolutionPtr,
                              *rhsVectorPtr_,
                              *jacobianMatrixPtr_,
                              *NewtonVectorPtr_,
                              *gradVectorPtr_,
                              *lasSysPtr_,
                              *this);
    }

    //////////////////////////////////////////////////////////////////////////
    // erkeite: Group handling required by 2-level Newton:
    // Reset up the corresponding group as well
    if (nlsTmpPtr==0)
    {
      if (Teuchos::is_null(groupPtr_))
      {
        groupPtr_ = Teuchos::rcp(new N_NLS_LOCA::Group(globalDataPtr_,
                                                       *sharedSystemPtr_,
                                                       getLoader(),
                                                       *outMgrPtr_,
                                                       getAnalysisManager())
                                 );
      }
      else
      {
        Vector tmpVec(*dsPtr_->nextSolutionPtr, *lasSysPtr_);
        groupPtr_->setX(tmpVec);
      }
    }
    else
    {
      copiedGroupFlag_ = true;
      Interface * nlsOtherPtr = dynamic_cast<Interface*>(nlsTmpPtr);
      groupPtr_ = nlsOtherPtr->getSolutionGroup();
    }
    // End of block needed by 2-level Newton.
    //////////////////////////////////////////////////////////////////////////

    int solverType = paramsPtr->getNoxSolverType();
    bool continuationSpecified = paramsPtr->getContinuationSpecifiedFlag();

    if ((mode_ == Nonlinear::DC_OP || mode_ == Nonlinear::DC_SWEEP) && !continuationSpecified)
    {
      return spiceStrategy(paramsPtr);
    }
    else
    {

      // Setting the nonContinuation flag is required to prevent incorrect
      // Device::DeviceMgr::setParam calls.
      if (solverType==0)
      {
        groupPtr_->setNonContinuationFlag (true);
      }
      else
      {
        groupPtr_->setNonContinuationFlag (false);
      }

      if (DEBUG_NONLINEAR)
        dout() << "solverType is " << solverType << std::endl;

      // (0) Standard Newton Method Solve (includes line search and
      // trust region based methods
      if (solverType == 0)
      {
        return stdNewtonSolve(paramsPtr);
      }
      // (1) Natural Parameter Continuation
      else if (solverType == 1)
      {
        return naturalParameterContinuationSolve ( paramsPtr );
      }
      // (2) Mosfet specific continuation
      else if (solverType == 2)
      {
        return mosfetContinuationSolve ( paramsPtr );
      }
      else if (solverType == 3)  // GMIN stepping, simple specification
      {
        return gminSteppingSolve ( paramsPtr );
      }
      // (4) Mosfet specific continuation (ERK, new 2/21/2004)
      else if (solverType == 4)
      {
        return mosfetContinuationSolve2 (paramsPtr);
      }
      // (5) Mosfet:BSIM3:Inverter specific continuation (RPP, new 2/25/2004)
      else if (solverType == 5)
      {
        return mosfetContinuationSolve3 (paramsPtr);
      }
      // (6) Mosfet:BSIM3:Inverter specific continuation (RPP, new 2/25/2004)
      else if (solverType == 6)
      {
        return mosfetContinuationSolve4 (paramsPtr);
      }  // Block gainscale
      else if (solverType == 7)
      {
        return blockGainscaleMosfetSolve ( paramsPtr );
      }  // Test suite
      else if (solverType == 8)
      {
        return mosfetContinuationSolve5 (paramsPtr);
      } // Pseudo Transient
      else if (solverType == 9)
      {
        return pseudoTransientSolve( paramsPtr );
      }  // continuation = 4 + power node
      else if (solverType == 10)
      {
        return mosfetContinuationSolve6 (paramsPtr);
      }
      else if (solverType == 33)  // artificial parameter
      {
        return artificialParameterHomotopy (paramsPtr);
      } 
      else if (solverType == 34)  // source stepping, simple specification
      {
        return sourceSteppingSolve ( paramsPtr );
      }// End of if (solverType == )
    }
  } // try
  catch (const char* error_msg) 
  {
    std::string nox_error = "NOX Error";
    std::string err_msg = std::string(error_msg);
    if (err_msg == nox_error) 
    {
      Report::DevelFatal()
        << "Caught a NOX Exception in Interface::solve()!";
    }
    else // Otherwise, rethrow...
    {
      throw;
    }
  }
#ifndef Xyce_CHARON
  catch (const std::exception& e) 
  {
    dout() << e.what() << std::endl;
    Report::DevelFatal() 
      << "Caught std::exception in Interface::solve()!";
  }
  catch (...) 
  {
    Report::DevelFatal()
      << "Caught Unknown Exception in Interface::solve()!";
  }
#endif

  // Should never get this far
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : Interface::icCont
// Purpose       :
// Special Notes : returns true if IC is being used.
//
//                 The "found" variable indicates if any of the nodes specified
//                 in the dcop start file were found in this circuit.  If not,
//                 then don't bother with this.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystem Modeling
// Creation Date : 09/15/07
//-----------------------------------------------------------------------------
bool Interface::icCont (ParameterSet* paramsPtr)
{
  bool usedIC(false);

#ifdef Xyce_DEBUG_IC
  dout() << "NOX_Interface:  Inside continuation=0 .IC code." << std::endl;
#endif
  int found = 0;
  int icType;
  IO::InitialConditionsData::NodeLidValueMap & op = initialConditionsManager_->getICData(found, icType);

  // The builder may need to update op due to analysis type (ex. embedded sampling)
  (lasSysPtr_->builder()).createInitialConditionOp( op );

  usedIC = (icType==IO::InitialConditionsData::IC_TYPE_IC && found > 0);
  if (usedIC)
  {
    bool useGminStepping=false;
    Teuchos::RCP<AugmentLinSys> als =
      paramsPtr->createAugmentLinearSystem(lasSysPtr_, op, useGminStepping);
    groupPtr_->setAugmentLinearSystem(true, als);
  }
  return usedIC;
}

//-----------------------------------------------------------------------------
// Function      : Interface::icCont3
// Purpose       : IC with gmin stepping
// Special Notes : returns true if IC is being used.
//
//                 The "found" variable indicates if any of the nodes specified
//                 in the dcop start file were found in this circuit.  If not,
//                 then don't bother with this.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/29/2012
//-----------------------------------------------------------------------------
bool Interface::icCont3 (ParameterSet* paramsPtr)
{
  bool usedIC(false);

#ifdef Xyce_DEBUG_IC
  dout() << "NOX_Interface:  Inside continuation=3 .IC code." << std::endl;
#endif
  int found = 0;
  int icType;
  IO::InitialConditionsData::NodeLidValueMap & op = initialConditionsManager_->getICData(found, icType);

  // The builder may need to update op due to analysis type (ex. embedded sampling)
  (lasSysPtr_->builder()).createInitialConditionOp( op );

  usedIC = (icType==IO::InitialConditionsData::IC_TYPE_IC && found > 0);
  if (usedIC)
  {
    bool useGminStepping=true;
    Teuchos::RCP<AugmentLinSys> als =
      paramsPtr->createAugmentLinearSystem(lasSysPtr_, op, useGminStepping);
    groupPtr_->setAugmentLinearSystem(true, als);
  }
  return usedIC;
}

//-----------------------------------------------------------------------------
// Function      : Interface::nodesetCont0
// Purpose       :
// Special Notes : returns true if is being used.
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystem Modeling
// Creation Date : 09/15/07
//-----------------------------------------------------------------------------
bool Interface::nodesetCont0 (ParameterSet* paramsPtr)
{
  bool usedNODESET(false);

#ifdef Xyce_DEBUG_IC
  dout() << "NOX_Interface:  Inside continuation=0 .NODESET code (case 1)" << std::endl;
#endif
  int found = 0;
  int icType;
  IO::InitialConditionsData::NodeLidValueMap & op = initialConditionsManager_->getICData(found, icType);

  // The builder may need to update op due to analysis type (ex. embedded sampling)
  (lasSysPtr_->builder()).createInitialConditionOp( op );

  usedNODESET = (icType==IO::InitialConditionsData::IC_TYPE_NODESET && found > 0);
  if (usedNODESET)
  {
    Teuchos::RCP<AugmentLinSys> als =
      paramsPtr->createAugmentLinearSystem(lasSysPtr_, op);

    groupPtr_->setAugmentLinearSystem(true, als);
    NOX::StatusTest::StatusType status = solverPtr_->solve();

    // Create a reset solver after performing the initial nodeset solve.
    groupPtr_->setAugmentLinearSystem(false, Teuchos::null);
    solverPtr_->reset(groupPtr_->getX());
    getAnalysisManager().completeOPStartStep();
    firstSolveComplete_ = true;
  }
  return usedNODESET;
}


//-----------------------------------------------------------------------------
// Function      : Interface::nodesetCont1
// Purpose       :
// Special Notes : returns true if .NODESET is being used.
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystem Modeling
// Creation Date : 09/15/07
//-----------------------------------------------------------------------------
bool Interface::nodesetCont1 (ParameterSet* paramsPtr)
{
  bool usedNODESET(false);

#ifdef Xyce_DEBUG_IC
  dout() << "NOX_Interface:  Inside continuation=1 .NODESET code (case 2)" << std::endl;
#endif

  int found = 0;
  int icType;
  IO::InitialConditionsData::NodeLidValueMap & op = initialConditionsManager_->getICData(found, icType);

  // The builder may need to update op due to analysis type (ex. embedded sampling)
  (lasSysPtr_->builder()).createInitialConditionOp( op );

  usedNODESET = (icType==IO::InitialConditionsData::IC_TYPE_NODESET && found > 0);
  if (usedNODESET)
  {
    // Set up nox nonlinear solver
    if (Teuchos::is_null(solverPtr_))
    {
      solverPtr_ = NOX::Solver::buildSolver(groupPtr_,
          paramsPtr->getStatusTests(),
          paramsPtr->getNoxParams());
    }

    Teuchos::RCP<AugmentLinSys> als =
      paramsPtr->createAugmentLinearSystem(lasSysPtr_, op);

    groupPtr_->setAugmentLinearSystem(true, als);
    NOX::StatusTest::StatusType status = solverPtr_->solve();

    firstSolveComplete_ = true;
    groupPtr_->setAugmentLinearSystem(false, Teuchos::null);
    getAnalysisManager().completeOPStartStep();
  }
  return usedNODESET;
}

//-----------------------------------------------------------------------------
// Function      : Interface::takeFirstSolveStep
// Purpose       : same as Interface::solve, except that solverPtr_->iterate is
//                 called instead of solverPtr_->solve.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::takeFirstSolveStep (Nonlinear::NonLinearSolver * nlsTmpPtr)
{
  // For base object
  Nonlinear::NonLinearSolver::resetCountersAndTimers_();

  // Pick the parameter set to use.
  ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == Nonlinear::TRANSIENT))
    paramsPtr = &transientParams_;
  else if ((usemode_) && (mode_ == Nonlinear::HB_MODE))
    paramsPtr = &hbParams_;
  else if ((usemode_) && (mode_ == Nonlinear::DC_NLPOISSON))
    paramsPtr = &nlpParams_;
  else
    paramsPtr = &dcParams_;

  // set the xyce return codes:
  paramsPtr->setStatusTestReturnCodes(retCodes_);

  if (Teuchos::is_null(globalDataPtr_))
    globalDataPtr_ = LOCA::createGlobalData(paramsPtr->getAllParams());

  // Set up the shared system (we have to redo this every time because
  // the object pointed to by nextSolVectorPtrPtr may have changed.
  delete sharedSystemPtr_;
  sharedSystemPtr_ = new SharedSystem(*dsPtr_->nextSolutionPtr,
                                      *rhsVectorPtr_,
                                      *jacobianMatrixPtr_,
                                      *NewtonVectorPtr_,
                                      *gradVectorPtr_,
                                      *lasSysPtr_,
                                      *this);

  // Reset up the corresponding group as well
  //delete groupPtr_;
  if (nlsTmpPtr==0)
  {
    if (Teuchos::is_null(groupPtr_))
    {
      dout() << "takeFirstSolveStep: allocating a new group!" << std::endl;
      groupPtr_ = Teuchos::rcp(new N_NLS_LOCA::Group(globalDataPtr_,
                                                     *sharedSystemPtr_,
                                                     getLoader(),
                                                     *outMgrPtr_,
                                                     getAnalysisManager()));
    }
    else
    {
      dout() << "takeFirstSolveStep: using the old group!" << std::endl;
      Vector tmpVec(*dsPtr_->nextSolutionPtr, *lasSysPtr_);
      groupPtr_->setX(tmpVec);
    }
  }
  else
  {
    dout() << "takeFirstSolveStep: copying over the passed group!" << std::endl;
    copiedGroupFlag_ = true;
    Interface * nlsOtherPtr = dynamic_cast<Interface*>(nlsTmpPtr);
    groupPtr_ = nlsOtherPtr->getSolutionGroup();
  }

  // Set up solver
  if (Teuchos::is_null(solverPtr_))
    solverPtr_ = NOX::Solver::buildSolver(groupPtr_,
                                          paramsPtr->getStatusTests(),
                                          paramsPtr->getNoxParams());
  else
    solverPtr_->reset(groupPtr_->getX());

  // Solve
  NOX::StatusTest::StatusType status = solverPtr_->step();

  // Return the solution status
  return (status == NOX::StatusTest::Converged) ? 1 : -1;
}

//-----------------------------------------------------------------------------
// Function      : Interface::takeOneSolveStep
// Purpose       : same as Interface::takeFirstSolveStep, except that none of the
//                 set up stuff (like allocating the solverPtr) is done here.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::takeOneSolveStep ()
{
  // Solve
  NOX::StatusTest::StatusType status = solverPtr_->step();

  // Return the solution status
  return (status == NOX::StatusTest::Converged) ? 1 : -1;
}

//-----------------------------------------------------------------------------
// Function      : Interface::getNumIterations
// Purpose       :
// Special Notes :
// Return Type   : Integer (current number of nonlinear iterations)
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::getNumIterations() const
{
  // Pick the parameter set to use.
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == Nonlinear::TRANSIENT))
    paramsPtr = &transientParams_;
  else if ((usemode_) && (mode_ == Nonlinear::HB_MODE))
    paramsPtr = &hbParams_;
  else if ((usemode_) && (mode_ == Nonlinear::DC_NLPOISSON))
    paramsPtr = &nlpParams_;
  else
    paramsPtr = &dcParams_;

  int solverType = paramsPtr->getNoxSolverType();

  if ((!Teuchos::is_null(solverPtr_)) && (solverType == 0))
    return solverPtr_->getNumIterations();
  else if ((!Teuchos::is_null(solverPtr_)) && (solverType == 1))
    return solverPtr_->getNumIterations();
  else if ((!Teuchos::is_null(solverPtr_)) && (solverType == 9))
    return solverPtr_->getNumIterations();
  else if ((!Teuchos::is_null(stepperPtr_)) && (solverType != 0))
  {
    return stepperPtr_->getSolver()->getNumIterations();
  }

  // Sometimes this is called before solve() itself, in which calse
  // the solverPtr_ has not yet been initialized, so we just return 0.
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : Interface::getMaxNormF() const
// Purpose       :
// Special Notes :
// Return Type   : double (norm of F)
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
double Interface::getMaxNormF() const
{
  // Pick the parameter set to use.
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == Nonlinear::TRANSIENT))
  {
    paramsPtr = &transientParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::HB_MODE))
  {
    paramsPtr = &hbParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::DC_NLPOISSON))
  {
    paramsPtr = &nlpParams_;
  }
  else
  {
    paramsPtr = &dcParams_;
  }

  double maxNormF = paramsPtr->getMaxNormF();
  return maxNormF;
}


//-----------------------------------------------------------------------------
// Function      : Interface::getMaxNormFindex() const
// Purpose       :
// Special Notes :
// Return Type   : int (vector index norm of F)
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::getMaxNormFindex() const
{
  // Pick the parameter set to use.
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == Nonlinear::TRANSIENT))
  {
    paramsPtr = &transientParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::HB_MODE))
  {
    paramsPtr = &hbParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::DC_NLPOISSON))
  {
    paramsPtr = &nlpParams_;
  }
  else
  {
    paramsPtr = &dcParams_;
  }

  int maxNormFindex = paramsPtr->getMaxNormFindex();
  return maxNormFindex;
}

//-----------------------------------------------------------------------------
// Function      : Interface::getDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2007
//-----------------------------------------------------------------------------
int Interface::getDebugLevel() const
{
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == Nonlinear::TRANSIENT))
  {
    paramsPtr = &transientParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::HB_MODE))
  {
    paramsPtr = &hbParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::DC_NLPOISSON))
  {
    paramsPtr = &nlpParams_;
  }
  else
  {
    paramsPtr = &dcParams_;
  }

  return (paramsPtr->getDebugLevel());
}

//-----------------------------------------------------------------------------
// Function      : Interface::getScreenOutputFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2007
//-----------------------------------------------------------------------------
bool Interface::getScreenOutputFlag () const
{
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == Nonlinear::TRANSIENT))
  {
    paramsPtr = &transientParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::HB_MODE))
  {
    paramsPtr = &hbParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::DC_NLPOISSON))
  {
    paramsPtr = &nlpParams_;
  }
  else
  {
    paramsPtr = &dcParams_;
  }

  return (paramsPtr->getScreenOutputFlag());
}

//-----------------------------------------------------------------------------
// Function      : Interface::getDebugMinTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2007
//-----------------------------------------------------------------------------
double Interface::getDebugMinTime() const
{
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == Nonlinear::TRANSIENT))
  {
    paramsPtr = &transientParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::HB_MODE))
  {
    paramsPtr = &hbParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::DC_NLPOISSON))
  {
    paramsPtr = &nlpParams_;
  }
  else
  {
    paramsPtr = &dcParams_;
  }

  return (paramsPtr->getDebugMinTime());
}

//-----------------------------------------------------------------------------
// Function      : Interface::getDebugMaxTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2007
//-----------------------------------------------------------------------------
double Interface::getDebugMaxTime() const
{
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == Nonlinear::TRANSIENT))
  {
    paramsPtr = &transientParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::HB_MODE))
  {
    paramsPtr = &hbParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::DC_NLPOISSON))
  {
    paramsPtr = &nlpParams_;
  }
  else
  {
    paramsPtr = &dcParams_;
  }

  return (paramsPtr->getDebugMaxTime());
}

//-----------------------------------------------------------------------------
// Function      : Interface::getDebugMinTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2007
//-----------------------------------------------------------------------------
int Interface::getDebugMinTimeStep() const
{
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == Nonlinear::TRANSIENT))
  {
    paramsPtr = &transientParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::HB_MODE))
  {
    paramsPtr = &hbParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::DC_NLPOISSON))
  {
    paramsPtr = &nlpParams_;
  }
  else
  {
    paramsPtr = &dcParams_;
  }

  return (paramsPtr->getDebugMinTimeStep());
}

//-----------------------------------------------------------------------------
// Function      : Interface::getDebugMaxTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2007
//-----------------------------------------------------------------------------
int Interface::getDebugMaxTimeStep() const
{
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == Nonlinear::TRANSIENT))
  {
    paramsPtr = &transientParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::HB_MODE))
  {
    paramsPtr = &hbParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::DC_NLPOISSON))
  {
    paramsPtr = &nlpParams_;
  }
  else
  {
    paramsPtr = &dcParams_;
  }

  return (paramsPtr->getDebugMaxTimeStep());
}

//-----------------------------------------------------------------------------
// Function      : Interface::getMMFormat
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/2/2011
//-----------------------------------------------------------------------------
bool Interface::getMMFormat () const
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : Interface::isFirstContinuationParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::isFirstContinuationParam() const
{
  return isFirstContinuationParam_;
}

//-----------------------------------------------------------------------------
// Function      : Interface::isFirstSolveComplete
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::isFirstSolveComplete() const
{
  return firstSolveComplete_;
}

//-----------------------------------------------------------------------------
// Function      : Interface::getContinuationStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::getContinuationStep() const
{
  if (!Teuchos::is_null(stepperPtr_))
  {
    return stepperPtr_->getStepNumber();
  }
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : Interface::getContinuationStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::getParameterNumber() const
{
  return iParam_;
}

//-----------------------------------------------------------------------------
// Function      : Interface::setAnalysisMode
//
// Purpose       : Specify the analysis mode to be used by the nonlinear
//                 solver in the next call to solve(). This *may* affect
//                 the parameters used by the solver.
//
// See Also      : setOptions, setTranOptions
//
// - Input Arguments -
//
//    mode       : Mode to be used in the next nonlinear solve.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void Interface::setAnalysisMode(Nonlinear::AnalysisMode mode)
{
  mode_ = mode;
}

//-----------------------------------------------------------------------------
// Function      : Interface::resetAll
// Purpose       : This is used when Xyce is doing a STEP loop, and
//                 needs to act like it is at the beginning of a transient
//                 simulation again, for the next parameter in the STEP loop.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void Interface::resetAll (Nonlinear::AnalysisMode mode)
{
  setAnalysisMode(mode);

  firstSolveComplete_ = false;
  isFirstContinuationParam_ = true;

  rhsVectorPtr_->putScalar(0.0);
  NewtonVectorPtr_->putScalar(0.0);
  gradVectorPtr_->putScalar(0.0);
  dsPtr_->setZeroHistory();

  stepperPtr_ = Teuchos::null;
  groupPtr_   = Teuchos::null;
}

//-----------------------------------------------------------------------------
// Function      : Interface::computeF()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::computeF()
{
  return Nonlinear::NonLinearSolver::rhs_();
}

//-----------------------------------------------------------------------------
// Function      : Interface::computeNewton
// Purpose       : Set up the parameters for the linear solver and then
//                 call newton_()
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::computeNewton(Teuchos::ParameterList& params)
{
  if (mode_ == Nonlinear::DC_OP && setAZ_Tol_DC)
  {
    lasSolverRCPtr_->setTolerance(params.get("Tolerance", 1.0e-12));
  }
  else if (mode_ == Nonlinear::TRANSIENT && setAZ_Tol_Transient)
  {
    lasSolverRCPtr_->setTolerance(params.get("Tolerance", 1.0e-12));
  }
  else if (mode_ == Nonlinear::HB_MODE && setAZ_Tol_DC)
  {
    lasSolverRCPtr_->setTolerance(params.get("Tolerance", 1.0e-12));
  }

  return Nonlinear::NonLinearSolver::newton_();
}

//-----------------------------------------------------------------------------
// Function      : Interface::computeJacobian()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::computeJacobian()
{
  bool status = Nonlinear::NonLinearSolver::jacobian_();
  return status;
}

//-----------------------------------------------------------------------------
// Function      : Interface::applyJacobian(const N_LAS_Vector& input, N_LAS_Vector& result)
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::applyJacobian(const Linear::Vector& input, Linear::Vector& result)
{
  return Nonlinear::NonLinearSolver::applyJacobian(input,result);
}

//-----------------------------------------------------------------------------
// Function      : Interface::computeDfDpMulti	
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::computeDfDpMulti	
  (const std::vector< int > & paramIDs, 
   NOX::Abstract::MultiVector & dfdp, 
   bool isValidF)
{
  // this could be cleaner:
  ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == Nonlinear::TRANSIENT))
  {
    paramsPtr = &transientParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::HB_MODE))
  {
    paramsPtr = &hbParams_;
  }
  else if ((usemode_) && (mode_ == Nonlinear::DC_NLPOISSON))
  {
    paramsPtr = &nlpParams_;
  }
  else
  {
    paramsPtr = &dcParams_;
  }

  Teuchos::RCP<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();
  Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");
  const std::string strArg("Skip df/dp");
  bool skipDFDP = stepperList.get(strArg, false);

  if (!skipDFDP)
  {
    // get the separate df/dp, dq/dp and db/dp
    int difference = SENS_FWD;
    bool forceFD = false;
    bool forceDeviceFD = false;
    bool forceAnalytic = false;

    double sqrtEta = 1.0e-8;

    // populate the paramNameVec.
    LOCA::ParameterVector locaPVec = groupPtr_->getParams();

    int size = locaPVec.length();
    std::vector<std::string> paramNameVec(size);

    for (int i = 0; i < locaPVec.length(); ++i) 
    {
      int index = paramIDs[i];
      paramNameVec[i] = locaPVec.getLabel(index);
    }

    loadSensitivityResiduals (difference, 
        forceFD, forceDeviceFD, forceAnalytic, 
        sqrtEta, netlistFilename_, 
        *dsPtr_, *nonlinearEquationLoader_, paramNameVec, getAnalysisManager());

    // get the complete residual (i.e. assemble the 3 vectors)
    nonlinearEquationLoader_->loadSensitivityResiduals();

    // now copy residuals back out into the NOX multivector, dfdp.
    Linear::MultiVector * sensRHSPtrVector = dsPtr_->sensRHSPtrVector;
    for (int i = 0; i < locaPVec.length(); ++i) 
    {
      // NOTE: LOCA stores f in dfdp[0], so the indexing for derivatives in dfdp starts at 1.
      int index = paramIDs[i];
      NOX::Abstract::Vector *DFDP = &dfdp[index+1]; 
      DFDP->init(0.0);

      Teuchos::RCP<Linear::Vector> tmp = Teuchos::rcp( sensRHSPtrVector->getNonConstVectorView(index) );
      Vector tmpNox(*tmp, *lasSysPtr_);
      (*DFDP) = tmpNox;
      DFDP->scale(-1.0);

      // debug output:
      //DFDP->print(std::cout);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Interface::getSolutionGroup
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Teuchos::RCP<N_NLS_LOCA::Group> Interface::getSolutionGroup ()
{
  return groupPtr_;
}

//-----------------------------------------------------------------------------
// Function      : Interface::getLoader
//
// Purpose       : LOCA needs access to loader to set parameters
//
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Loader::NonlinearEquationLoader& Interface::getLoader() const
{
  return *nonlinearEquationLoader_;
}

//-----------------------------------------------------------------------------
// Function      : Interface::resetStepper
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void Interface::resetStepper(const Teuchos::RCP<LOCA::GlobalData>& gd,
    const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& initialGuess,
    const Teuchos::RCP<NOX::StatusTest::Generic>& test,
    const Teuchos::RCP<Teuchos::ParameterList>& p)
{
  stepperPtr_ =
    Teuchos::rcp(new LOCA::Stepper(gd, initialGuess, test, p));
}

//-----------------------------------------------------------------------------
// Function      : Interface::getLocaFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::getLocaFlag ()
{
  // Pick the parameter set to use.
  const ParameterSet* paramsPtr;
  bool retCode;

  if ((usemode_) && (mode_ == Nonlinear::TRANSIENT))
  {
    firstSolveComplete_ = false;
    paramsPtr = &transientParams_;
    int solverType = paramsPtr->getNoxSolverType();
    retCode = false;
    if (solverType != 0) retCode = true;
  }
  else
  {
    if ((usemode_) && (mode_ == Nonlinear::HB_MODE))
    {
      paramsPtr = &hbParams_;
    }
    else if ((usemode_) && (mode_ == Nonlinear::DC_NLPOISSON))
    {
      paramsPtr = &nlpParams_;
    }
    else
    {
      paramsPtr = &dcParams_;
    }

    {
    int solverType = paramsPtr->getNoxSolverType();
    retCode=false;
    if (solverType != 0) retCode = true;
    }
  }

  return retCode;
}

}}} // namespace N_NLS_NOX

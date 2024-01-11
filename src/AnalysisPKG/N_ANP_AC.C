//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose       : AC analysis functions.
// Special Notes :
// Creator       : Ting Mei
// Creation Date :  7/11
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iosfwd>
#include <iomanip>
#include <fstream>
#include <complex>

#include <N_ANP_AC.h>

#include<N_ANP_SweepParamFreeFunctions.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ERH_Message.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>
#include <N_IO_InitialConditions.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockSystemHelpers.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_Builder.h>
#include <N_LAS_Graph.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_System.h>
#include <N_LAS_SystemHelpers.h>
#include <N_LAS_Solver.h>
#include <N_LAS_Problem.h>
#include <N_LAS_TranSolverFactory.h>
#include <N_LOA_Loader.h>
#include <N_NLS_Manager.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_WorkingIntegrationMethod.h>
#include <N_TIA_fwd.h>
#include <N_UTL_Diagnostic.h>
#include <N_TOP_Topology.h>
#include <N_UTL_Factory.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_Math.h>
#include <N_UTL_RFparams.h>
#include <N_UTL_SaveIOSState.h>
#include <N_UTL_Timer.h>
#include <N_UTL_MachDepParams.h>

#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>
#include <N_PDS_ParMap.h>

#include <N_DEV_DeviceMgr.h>
#include <N_UTL_ExtendedString.h>
#include <N_NLS_ReturnCodes.h>
#include <N_NLS_SensitivityResiduals.h>
#include <N_NLS_ObjectiveFunctions.h>

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;

namespace Xyce {
namespace Analysis {


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
bool ACExpressionGroup::getSolutionVal
    (const std::string & nodeName, std::complex<double> & retval)
{
  double real_val=0.0;
  double imag_val=0.0;
  int tmpGID = -1;

  tmpGID = getSolutionGID_(nodeName);
  if (tmpGID >= 0)
  {
    Linear::Vector & Xreal = X_.block( 0 );
    Linear::Vector & Ximag = X_.block( 1 );
    real_val = Xreal.getElementByGlobalIndex(tmpGID, 0); 
    imag_val = Ximag.getElementByGlobalIndex(tmpGID, 0); 
  }

  Xyce::Parallel::AllReduce(getComm().comm(), MPI_SUM, &real_val, 1);
  Xyce::Parallel::AllReduce(getComm().comm(), MPI_SUM, &imag_val, 1);

  retval = std::complex<double>(real_val,imag_val);
  return (tmpGID>=0);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setTimeIntegratorOptions
// Purpose       :
// Special Notes : These are from '.options timeint'
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool AC::setTimeIntegratorOptions(
  const Util::OptionBlock &     option_block)
{
  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; ++it)
  {
    const Util::Param &param = *it;

    if (param.uTag() == "DEBUGLEVEL" )
      IO::setTimeIntegratorDebugLevel(analysisManager_.getCommandLine(), param.getImmutableValue<int>());
    else if (nonlinearManager_.setReturnCodeOption(param))
      ;
    else if (tiaParams_.setTimeIntegratorOption(param))
      ;
    else if (setDCOPOption(param))
      ;
    else
    {
      Report::UserError() << param.uTag() << " is not a recognized time integration option";
      return false;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : setACLinSolOptionBlock
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/23/2018
//-----------------------------------------------------------------------------
///
/// Saves the LINSOL-AC parsed options block in the factory.
///
/// @invariant Overwrites any previously specified LINSOL-AC option block.
///
/// @param option_block parsed option block
///
bool AC::setACLinSolOptions(const Util::OptionBlock &option_block)
{
  acLinSolOptionBlock_ = option_block;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : setDCLinSolOptionBlock
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 12/10/2020
//-----------------------------------------------------------------------------
///
/// Saves the LINSOL parsed options block in the factory.
///
/// @invariant Overwrites any previously specified LINSOL option block.
///
/// @param option_block parsed option block
///
bool AC::setDCLinSolOptions(const Util::OptionBlock &option_block)
{
  dcLinSolOptionBlock_ = option_block;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::setDataStatements
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/5/18
//-----------------------------------------------------------------------------
bool AC::setDataStatements(const Util::OptionBlock & paramsBlock)
{
  return processDataStatements(paramsBlock, dataNamesMap_, dataTablesMap_);
}

//-----------------------------------------------------------------------------
// Function      : AC::convertDataToSweepParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/5/18
//-----------------------------------------------------------------------------
bool  AC::convertDataToSweepParams()
{
  return convertData(acSweepVector_, dataNamesMap_, dataTablesMap_);
}

//-----------------------------------------------------------------------------
// Function      : AC::setSensAnalysisParams
// Purpose       : This function processes the .SENS line
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AC::setSensAnalysisParams(const Util::OptionBlock & OB)
{
  Util::ParamList::const_iterator iter = OB.begin();
  Util::ParamList::const_iterator end   = OB.end();
  bool isObjVar = false, isObjFunc = false;
  numSensParams_ = 0;

  for ( ; iter != end; ++ iter)
  {
    isObjFunc = std::string( iter->uTag() ,0,9) == "ACOBJFUNC";
    isObjVar = std::string(iter->uTag(), 0, 7) == "OBJVARS";
    if(isObjFunc || isObjVar)
    {
      Xyce::Nonlinear::objectiveFunctionData<std::complex<double> > * ofDataPtr = new Xyce::Nonlinear::objectiveFunctionData<std::complex<double> >();
      ofDataPtr->objFuncString = isObjFunc ? iter->stringValue() : "{V(" + iter->stringValue() + ")}";
      objFuncDataVec_.push_back(ofDataPtr);
      objFuncStrings_.push_back(ofDataPtr->objFuncString);
    }
    else if ( std::string( iter->uTag() ,0,5) == "PARAM") // this is a vector
    {
      ExtendedString tag = iter->stringValue();
      // set up the initial skeleton of the maps:
      ++numSensParams_;
      paramNameVec_.push_back(tag);

      int sz = tag.size();
      if (sz > maxParamStringSize_)
      {
        maxParamStringSize_ = sz;
      }
    }
    else
    {
      Xyce::Report::UserWarning() << iter->uTag()
        << " is not a recognized sensitivity solver option.\n" << std::endl;
    }
  }

  objFuncGiven_ = !objFuncDataVec_.empty();
  if (objFuncGiven_)
  {
    IO::OutputMgr & output_manager = outputManagerAdapter_.getOutputManager();

    Xyce::Nonlinear::setupObjectiveFunctions(
      analysisManager_.getExpressionGroup(),
         objFuncDataVec_, output_manager, linearSystem_, analysisManager_.getCommandLine() );
  }

  // check that the .sens line exists and actually has all of the components it needs (objective functions and parameters)
  if(sensFlag_ && !objFuncGiven_)
  {
    Report::UserError0() << "No objective functions specified for .SENS";
    return false;
  }

  if(sensFlag_ && !(numSensParams_ > 0))
  {
    Report::UserError0() << "No objective function parameters specified for .SENS";
    return false;    
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::setSensitivityOptions
// Purpose       :
// Special Notes : These are from '.options sensitivity'
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AC::setSensitivityOptions(const Util::OptionBlock & option_block)
{
  Util::ParamList::const_iterator it  = option_block.begin();
  Util::ParamList::const_iterator end = option_block.end();
  for ( ; it != end; ++ it)
  {
    const Util::Param &param = (*it);

    if (param.uTag() == "ADJOINT")
    {
      solveAdjointSensitivityFlag_ =
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if (param.uTag() == "DIRECT")
    {
      solveDirectSensitivityFlag_ =
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "OUTPUTSCALED")
    {
      outputScaledFlag_ =
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "OUTPUTUNSCALED")
    {
      outputUnscaledFlag_ =
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "STDOUTPUT")
    {
      stdOutputFlag_ =
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "FORCEFD") // deprecated, not applicable
    {
      forceFD_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "FORCEDEVICEFD") // use the device numerical sens
    {
      forceDeviceFD_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "FORCEANALYTIC")
    {
      forceAnalytic_= 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "REUSEFACTORS")
    {
      reuseFactors_ = (*it).getImmutableValue<double>();
    }
    else
    {
      // do nothing.
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::AC
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 5/11
//-----------------------------------------------------------------------------
AC::AC(
  AnalysisManager &                     analysis_manager,
  Linear::System &                      linear_system,
  Nonlinear::Manager &                  nonlinear_manager,
  Device::DeviceMgr &                   device_manager,
  Loader::Loader &                      loader,
  Topo::Topology &                      topology,
  IO::InitialConditionsManager &        initial_conditions_manager)
  : AnalysisBase(analysis_manager, "AC"),
    StepEventListener(&analysis_manager),
    analysisManager_(analysis_manager),
    loader_(loader),
    linearSystem_(linear_system),
    nonlinearManager_(nonlinear_manager),
    topology_(topology),
    deviceManager_(device_manager),
    initialConditionsManager_(initial_conditions_manager),
    outputMOR_(analysisManager_.getNetlistFilename()),
    outputManagerAdapter_(analysis_manager.getOutputManagerAdapter()),
    bVecRealPtr(linearSystem_.builder().createVector()),
    bVecImagPtr(linearSystem_.builder().createVector()),
    acLoopSize_(0),
    stepFlag_(false),
    type_("DEC"),
    np_(10.0),
    fStart_(1.0),
    fStop_(1.0),
    stepMult_(0.0),
    fstep_(0.0),
    currentFreq_(0.0),
    sparcalc_(false),
    numPorts_(0),
    hParamsRequested_(false),
    sParamsRequested_(false),
    zParamsRequested_(false),
    ACMatrix_(0),
    B_(0),
    X_(0),

    dbdpVecRealPtr(0),
    dbdpVecImagPtr(0),
    dOdxVecRealPtr(0),
    dOdxVecImagPtr(0),
    dCdp_(0),
    dGdp_(0),
    origC_(0),
    origG_(0),
    dBdp_(0),
    dXdp_(0),
    lambda_(0),
    dOdXreal_(0),
    dOdXimag_(0),
    savedX_(0),
    sensRhs_(0),

    blockSolver_(0),
    blockProblem_(0),

    sensFlag_(analysis_manager.getSensFlag()),
    solveAdjointSensitivityFlag_(true),
    solveDirectSensitivityFlag_(false),
    outputScaledFlag_(false),
    outputUnscaledFlag_(false),
    maxParamStringSize_(0),
    stdOutputFlag_(false),
    numSensParams_(0),

    objFuncGiven_(false),
    objFuncGIDsetup_(false),

    forceFD_(false),
    forceDeviceFD_(false),
    forceAnalytic_(false),
    newLowMem_(false),
    sparseAdjointStorage_(true),
    reuseFactors_(true)

{
  bVecRealPtr->putScalar(0.0);
  bVecImagPtr->putScalar(0.0);

  RFparams_.insert(std::make_pair("Y",&Yparams_));
  RFparams_.insert(std::make_pair("S",&Sparams_));
  RFparams_.insert(std::make_pair("Z",&Zparams_));

  outputManagerAdapter_.setDotACSpecified(true);
}

//-----------------------------------------------------------------------------
// Function      : AC::~AC()
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 5/11
//-----------------------------------------------------------------------------
AC::~AC()
{
  delete bVecRealPtr;
  delete bVecImagPtr;
  delete ACMatrix_;
  delete B_;
  delete X_;
  delete blockSolver_;
  delete blockProblem_;

  if (sensFlag_)
  {
    delete origB_;
    delete dbdpVecRealPtr;
    delete dbdpVecImagPtr;
    delete dOdxVecRealPtr;
    delete dOdxVecImagPtr;
    delete dBdp_;
    delete dXdp_;
    delete sensRhs_;
    delete lambda_;
    delete dOdXreal_;
    delete dOdXimag_;
    delete savedX_;

    for (int iobj=0;iobj<objFuncDataVec_.size();++iobj)
    {
      delete objFuncDataVec_[iobj]->dOdXVectorRealPtr;
      objFuncDataVec_[iobj]->dOdXVectorRealPtr = 0;

      delete objFuncDataVec_[iobj]->dOdXVectorImagPtr;
      objFuncDataVec_[iobj]->dOdXVectorImagPtr = 0;

      delete objFuncDataVec_[iobj]->expPtr;
      objFuncDataVec_[iobj]->expPtr = 0;

      delete objFuncDataVec_[iobj];
      objFuncDataVec_[iobj] = 0;
    } 

    for (int ii=0;ii<numSensParams_;++ii) { delete dJdpVector_[ii]; }
    for (int ii=0;ii<numSensParams_;++ii) { delete dBdpVector_[ii]; }
  }
}

//-----------------------------------------------------------------------------
// Function      : AC::notify
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void AC::notify(const StepEvent &event)
{
  if (event.state_ == StepEvent::STEP_STARTED)
  {
    AnalysisBase::resetForStepAnalysis();

    stepFlag_ = true;

    analysisManager_.getStepErrorControl().resetAll(tiaParams_);

    bVecRealPtr->putScalar(0.0);
    bVecImagPtr->putScalar(0.0);
  }
}

//-----------------------------------------------------------------------------
// Function      : AC::setAnalysisParams
// Purpose       :
// Special Notes : These are from the .AC statement.
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 6/11
//-----------------------------------------------------------------------------
bool AC::setAnalysisParams(
  const Util::OptionBlock & paramsBlock)
{
  // Check for DATA first.  If DATA is present, then use the sweep functions,
  // rather than the AC specific built-in ones.
  if (isDataSpecified(paramsBlock))
  {
    dataSpecification_ = true;
    type_="TYPE";
    acSweepVector_.push_back(parseSweepParams(paramsBlock.begin(), paramsBlock.end()));
    return true;
  }

  bool retval = true;
  // -- original code ---
  for (Util::ParamList::const_iterator it = paramsBlock.begin(), end = paramsBlock.end(); it != end; ++it)
  {
    if ((*it).uTag() == "TYPE")
    {
      type_ = (*it).stringValue();
    }
    else if ((*it).uTag() == "NP")
    {
      np_ = (*it).getImmutableValue<double>();
      ExtendedString npStr((*it).stringValue());
      if ( !npStr.isInt() )
      {
        Report::UserError0() << "Points Value parameter on .AC line must be an integer";
        retval = false;
      }
    }
    else if ((*it).uTag() == "FSTART")
    {
      fStart_ = (*it).getImmutableValue<double>();
    }
    else if ((*it).uTag() == "FSTOP")
    {
      fStop_ = (*it).getImmutableValue<double>();
    }
  }

  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    dout() << section_divider << std::endl
           << "AC simulation parameters"
           << std::endl
           << "number of points  = " << np_ << std::endl
           << "starting frequency = " << fStart_ << std::endl
           << "stop frequency = " << fStop_
           << std::endl;
  }

  // error checking of parameters
  if ( np_ < 1 )
  {
    Report::UserError0() << "Points Value parameter on .AC line must be >= 1";
    retval = false;
  }
  if ( (fStart_ <=0 || fStop_ <= 0)  && (type_ == "DEC" || type_ == "OCT") )
  {
    Report::UserError0() << "Illegal values for start or end frequencies on .AC line. " <<
       "Both values must be > 0";
    retval = false;
  }
  if ( fStop_ < fStart_ )
  {
    Report::UserError0() << "End frequency must not be less than start frequency on .AC line";
    retval = false;
  }

  return retval;
}

//-----------------------------------------------------------------------------
// Function      : AC::setACLinOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 02/28/19
//-----------------------------------------------------------------------------
bool AC::setACLinOptions(const Util::OptionBlock & OB)
{
  for(Util::ParamList::const_iterator iterPL = OB.begin(), endPL = OB.end(); iterPL != endPL; ++iterPL )
  {
    ExtendedString tag = iterPL->tag();
    tag.toUpper();

    if ( tag == "SPARCALC"   ) 
    {
      // a .LIN analysis will happen if any of the .LIN lines have SPARCALC=1
      sparcalc_ = sparcalc_ | static_cast<bool> (iterPL->getImmutableValue<int>());

      // The output manager needs to know, for -r output, whether a
      // .LIN analysis is being done.
      outputManagerAdapter_.setEnableSparCalcFlag(sparcalc_);
    }
    else if ( tag == "LINTYPE")
    {
      const std::string linType = iterPL->getImmutableValue<std::string>();
      setRFParamsRequested(linType);
    }
    else
    {
      Report::UserError() << "Unrecognized option for .LIN line" << tag;
      return false;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::setRFParamsRequested()
// Purpose       : Determine which RF parameter types (S, Y or Z) have
//                 been requested by a .LIN or .PRINT AC line.  This
//                 helps with memory management and performance since the
//                 Sparams_ and Zparams_ matrices don't have to be resized,
//                 or converted from YParams_, if the netlist doesn't need
//                 S or Z parameters.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/31/2019
//-----------------------------------------------------------------------------
void AC::setRFParamsRequested(const std::string & type)
{
  if (type == "S")
    sParamsRequested_ = true;
  else if (type == "Z")
    zParamsRequested_ = true;
}

//-----------------------------------------------------------------------------
// Function      : AC::getDCOPFlag()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/24/2014
//-----------------------------------------------------------------------------
bool AC::getDCOPFlag() const
{
  return getIntegrationMethod() == TimeIntg::methodsEnum::NO_TIME_INTEGRATION; 
}

//-----------------------------------------------------------------------------
// Function      : AC::finalExpressionBasedSetup()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/4/2021
//-----------------------------------------------------------------------------
void AC::finalExpressionBasedSetup()
{
  return;
}

//-----------------------------------------------------------------------------
// Function      : AC::doRun()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 5/11
//-----------------------------------------------------------------------------
bool AC::doRun()
{
  return doInit() && doLoopProcess() && doFinish();
}

//-----------------------------------------------------------------------------
// Function      : AC::doInit()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 5/11
//-----------------------------------------------------------------------------
bool AC::doInit()
{
  bool bsuccess = true;

  // check if the "DATA" specification was used.  If so, create a new vector of
  // SweepParams, in the "TABLE" style.
  if (dataSpecification_)
  {
    if (!convertDataToSweepParams())
    {
      Report::UserFatal() << "Invalid data=<name> parameter on .AC line.";
      return false;
    }

    std::vector<SweepParam>::iterator begin = acSweepVector_.begin();
    std::vector<SweepParam>::iterator end = acSweepVector_.end();
    std::vector<SweepParam>::iterator it = begin;
    for ( ; it != end; ++it)
    {
      SweepParam &sweep_param = (*it);
      std::string name = (*it).name; Util::toUpper(name);
      if (name == "FREQ" || name == "HERTZ")
      {
        // frequency values for .AC must be > 0
        for (int i=0; i<(*it).valList.size(); ++i)
        {
          if ( (*it).valList[i] <= 0 )
          {
            Report::UserFatal() << "Frequency values in .DATA for .AC analysis must be > 0";
            return false;
          }
        }
      }
      else
      {
        loader_.getParamAndReduce(analysisManager_.getComm(), sweep_param.name);
      }
    }

    // now set up the looping, etc
    acLoopSize_ = setSweepLoopVals(begin, end);
  }
  else
  {
    acLoopSize_ = setupSweepParam_();
  }

  // Get set to do the operating point.
  baseIntegrationMethod_ = TimeIntg::methodsEnum::NO_TIME_INTEGRATION;
  analysisManager_.createTimeIntegratorMethod(tiaParams_, baseIntegrationMethod_);

  stepNumber = 0;
  setDoubleDCOPEnabled(loader_.isPDESystem());

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::INITIALIZE, AnalysisEvent::AC_IC));

  // set initial guess, if there is one to be set.
  // this setInitialGuess call is to up an initial guess in the
  // devices that have them (usually PDE devices).  This is different than
  // the "intializeProblem" call, which sets IC's.  (initial conditions are
  // different than initial guesses.
  loader_.setInitialGuess(analysisManager_.getDataStore()->nextSolutionPtr);

  // If available, set initial solution (.IC, .NODESET, etc).
  setInputOPFlag(
    initialConditionsManager_.setupInitialConditions(outputManagerAdapter_.getComm(),
                                                     topology_.getSolutionNodeNameMap(),
                                                     outputManagerAdapter_.getAliasNodeMap(),
                                                     *analysisManager_.getDataStore()->nextSolutionPtr,
                                                     linearSystem_));

  // Set a constant history for operating point calculation
  analysisManager_.getDataStore()->setConstantHistory();
  analysisManager_.getWorkingIntegrationMethod().obtainCorrectorDeriv();

  {
    Stats::StatTop _nonlinearStat("DC Nonlinear Solve");
    Stats::TimeBlock _nonlinearTimer(_nonlinearStat);

    // solving for DC op
    doHandlePredictor();
    loader_.updateSources();
    analysisManager_.getStepErrorControl().newtonConvergenceStatus = nonlinearManager_.solve();
    analysisManager_.getWorkingIntegrationMethod().stepLinearCombo ();
    gatherStepStatistics(stats_, nonlinearManager_.getNonlinearSolver(), analysisManager_.getStepErrorControl().newtonConvergenceStatus);
    analysisManager_.getStepErrorControl().evaluateStepError(loader_, tiaParams_);
  }

  if ( analysisManager_.getStepErrorControl().newtonConvergenceStatus <= 0)
  {
    Report::UserError() << "Solving for DC operating point failed! Cannot continue AC analysis";
    return false;
  }

  // only output DC op if the .op was specified in the netlist
  // or if this AC analysis is called from a .step loop
  if ( stepFlag_ || analysisManager_.getDotOpSpecified() )
  {
    outputManagerAdapter_.dcOutput(
        stepNumber,
        *analysisManager_.getDataStore()->nextSolutionPtr,
        *analysisManager_.getDataStore()->nextStatePtr,
        *analysisManager_.getDataStore()->nextStorePtr,
        *analysisManager_.getDataStore()->nextLeadCurrentPtr,
        *analysisManager_.getDataStore()->nextLeadDeltaVPtr,
        *analysisManager_.getDataStore()->nextLeadCurrentQDerivPtr,
        objectiveVec_,
          dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);
    outputManagerAdapter_.finishOutput();
  }

  // This for saving the data from the DC op.  different from the above where we are
  // concerned with generating normal output.
  // outputManagerAdapter_.outputDCOP(*analysisManager_.getDataStore()->nextSolutionPtr);
  initialConditionsManager_.outputDCOP(outputManagerAdapter_.getComm(), topology_.getSolutionNodeNameMap(), *analysisManager_.getDataStore()->nextSolutionPtr);

  loader_.loadBVectorsforAC (bVecRealPtr, bVecImagPtr);

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::FINISH, AnalysisEvent::AC_IC));

  // Create B matrix stamp

  deviceManager_.setSPAnalysisFlag( sparcalc_ );

  if  (sparcalc_ )
  {
    std::vector<int> tempVec;
    std::vector<double> tempZ0s;

    loader_.getBMatrixEntries(tempVec, portNumVec_, &tempZ0s);

    Parallel::Manager * pdsMgrPtr = analysisManager_.getPDSManager();

    int myPID = pdsMgrPtr->getPDSComm()->procID();
    Parallel::Machine comm =  pdsMgrPtr->getPDSComm()->comm();

  // Determine how many global ports there are by performing a sumAll()

    int  len = tempVec.size();

    Parallel::AllReduce(comm, MPI_SUM, &len, &numPorts_, 1); 
  
    if (myPID == 0)
    {
      if ( numPorts_ ==  0 )
        Report::UserFatal() << "No port device is found for S parameter analysis";
    }

    portPID_.assign(numPorts_, -1);

    Z0sVec_.assign(numPorts_, 0);

    for (int i=0; i<len; ++i)
    {
      portPID_[ portNumVec_[i] - 1] = myPID;
      Z0sVec_[ portNumVec_[i] - 1] = tempZ0s[i];
    }

    Xyce::Parallel::AllReduce(comm, MPI_MAX, &portPID_[0], numPorts_);

    Xyce::Parallel::AllReduce(comm, MPI_SUM, &Z0sVec_[0], numPorts_);       

    // error checking
    if (myPID == 0)
    {
      for (int i=0; i< numPorts_;  ++i)
      {
        if ( portPID_[i] == -1 )
        {
          Report::UserFatal() << "Did not find port " << i + 1 << " for .LIN analysis";
        }

        if ( Z0sVec_[i] < 0.0 )
        {
          Report::UserFatal() << " The negative impedance " << Z0sVec_[i] << " for port " << i + 1
                              <<  " is not supported for .LIN analysis";
        }
      }
    }

    for (int i=0; i<len; ++i)
    {
      portMap_[portNumVec_[i]] = std::pair<int, double> ( tempVec[i], tempZ0s[i]);
    }

    Yparams_.shape(numPorts_, numPorts_);
    // only reshape these matrices if S-, Z- or H-parameter output was requested
    // on the .LIN or .PRINT AC lines
    if (sParamsRequested_) { Sparams_.shape(numPorts_, numPorts_); }
    if (zParamsRequested_) { Zparams_.shape(numPorts_, numPorts_); }
    if (hParamsRequested_) { Hparams_.shape(numPorts_, numPorts_); }
  }

  if (sensFlag_ && !objFuncGIDsetup_)
  {
      IO::OutputMgr & output_manager = outputManagerAdapter_.getOutputManager();
      Parallel::Machine comm =  analysisManager_.getPDSManager()->getPDSComm()->comm();
      Xyce::Nonlinear::setupObjectiveFuncGIDs ( objFuncDataVec_, comm, topology_, output_manager );

      objFuncGIDsetup_ = true;
  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : AC::doLoopProcess()
// Purpose       : Conduct the stepping loop.
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 5/11
//-----------------------------------------------------------------------------
bool AC::doLoopProcess()
{
  updateLinearSystem_C_and_G_();

  createLinearSystem_();

  if(sensFlag_)
  {
    precomputeDCsensitivities_ ();
  }

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::INITIALIZE, AnalysisEvent::AC));

  for (int currentStep = 0; currentStep < acLoopSize_; ++currentStep)
  {
    if (dataSpecification_)
    {
      updateDataParams_(currentStep);
    }
    else
    {
      updateCurrentFreq_(currentStep);
    }

    static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_STARTED, AnalysisEvent::AC, currentFreq_, currentStep));

    updateLinearSystem_C_and_G_();
    updateLinearSystemFreq_();
    updateLinearSystemMagAndPhase_();

    bool stepAttemptStatus;
    {
      Stats::StatTop _ACsolveStat("AC Linear Solve");
      Stats::TimeBlock _AC_Timer(_ACsolveStat);

      stepAttemptStatus = solveLinearSystem_();
    }

    if (stepAttemptStatus)
    {
      if(sensFlag_)
      {
        bool sensSuccess = solveSensitivity_();
      }

      static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_SUCCESSFUL, AnalysisEvent::AC, currentFreq_, currentStep));
      doProcessSuccessfulStep();
    }
    else // stepAttemptStatus  (ie do this if the step FAILED)
    {
      static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_FAILED, AnalysisEvent::AC, currentFreq_, currentStep));
      doProcessFailedStep();
    }
  }

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::FINISH, AnalysisEvent::AC));



  return true;
}


//-----------------------------------------------------------------------------
// Function      : AC::createLinearSystem_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/20/2011
//-----------------------------------------------------------------------------
bool AC::createLinearSystem_()
{
  Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();

  RCP<Parallel::ParMap> baseMap = rcp(pds_manager.getParallelMap( Parallel::SOLUTION ), false);
  const Linear::Graph* baseFullGraph = pds_manager.getMatrixGraph(Parallel::JACOBIAN);

  int numBlocks = 2;
  int offset = baseMap->maxGlobalEntity() + 1;  // Use this offset to create a contiguous gid map for direct solvers.

  RCP<Parallel::ParMap> blockMap = Linear::createBlockParMap(numBlocks, *baseMap, 0, 0, offset);

  // Create a block vector
  delete B_;
  B_ = Xyce::Linear::createBlockVector(numBlocks, blockMap, baseMap);

  // -----------------------------------------------------
  // Now test block graphs.
  // -----------------------------------------------------

  std::vector<std::vector<int> > blockPattern(2);
  blockPattern[0].resize(2);
  blockPattern[0][0] = 0; blockPattern[0][1] = 1;
  blockPattern[1].resize(2);
  blockPattern[1][0] = 0; blockPattern[1][1] = 1;

  RCP<Linear::Graph> blockGraph = Linear::createBlockGraph( offset, blockPattern, *blockMap, *baseFullGraph );

  delete ACMatrix_;
  ACMatrix_ = Xyce::Linear::createBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), baseFullGraph );

  ACMatrix_->put( 0.0 ); // Zero out whole matrix.
  // Matrix will be loaded with nonzero (C,G) sub-matrices later.

  // Copy the values loaded into the blocks into the global matrix for the solve.
  ACMatrix_->assembleGlobalMatrix();

  B_->block( 0 ).update( 1.0, *bVecRealPtr, 0.0 );
  B_->block( 1 ).update( 1.0, *bVecImagPtr, 0.0 );

  delete X_;
  X_ = Xyce::Linear::createBlockVector (numBlocks, blockMap, baseMap);
  X_->putScalar( 0.0 );

  delete blockProblem_;
  blockProblem_ = Xyce::Linear::createProblem( ACMatrix_, X_, B_ );

  delete blockSolver_;
  Linear::TranSolverFactory factory;
  blockSolver_ = factory.create( acLinSolOptionBlock_, *blockProblem_, analysisManager_.getCommandLine() );

  if (sensFlag_)
  {
    dbdpVecRealPtr = linearSystem_.builder().createVector();
    dbdpVecImagPtr = linearSystem_.builder().createVector();
    dOdxVecRealPtr = linearSystem_.builder().createVector();
    dOdxVecImagPtr = linearSystem_.builder().createVector();

    dCdp_ = linearSystem_.builder().createMatrix ();
    dGdp_ = linearSystem_.builder().createMatrix ();

    origC_ = linearSystem_.builder().createMatrix ();
    origG_ = linearSystem_.builder().createMatrix ();

    origB_    = Xyce::Linear::createBlockVector(numBlocks, blockMap, baseMap);
    dBdp_     = Xyce::Linear::createBlockVector(numBlocks, blockMap, baseMap);
    dXdp_     = Xyce::Linear::createBlockVector(numBlocks, blockMap, baseMap);
    sensRhs_  = Xyce::Linear::createBlockVector(numBlocks, blockMap, baseMap);
    lambda_   = Xyce::Linear::createBlockVector(numBlocks, blockMap, baseMap);
    dOdXreal_ = Xyce::Linear::createBlockVector(numBlocks, blockMap, baseMap);
    dOdXimag_ = Xyce::Linear::createBlockVector(numBlocks, blockMap, baseMap);
    savedX_   = Xyce::Linear::createBlockVector(numBlocks, blockMap, baseMap);

    if (objFuncGiven_)
    {
      const Teuchos::RCP<Xyce::Util::mainXyceExpressionGroup> group = 
        Teuchos::rcp_dynamic_cast<Xyce::Util::mainXyceExpressionGroup>(analysisManager_.getExpressionGroup());

      Teuchos::RCP<ACExpressionGroup>  acGroup =
         Teuchos::rcp(new ACExpressionGroup(group, *X_));

      Teuchos::RCP<Xyce::Util::baseExpressionGroup> newGroup = acGroup;
 
      for (int iobj=0;iobj<objFuncDataVec_.size();++iobj)
      {
        objFuncDataVec_[iobj]->expPtr->setGroup( newGroup );
      }
    }

    dBdpVector_.resize(numSensParams_);
    for (int ii=0;ii<numSensParams_;++ii) { dBdpVector_[ii] = Xyce::Linear::createBlockVector(numBlocks, blockMap, baseMap); }

    dJdpVector_.resize(numSensParams_);
    for (int ii=0;ii<numSensParams_;++ii) { dJdpVector_[ii] = Xyce::Linear::createBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), baseFullGraph); }

  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::updateLinearSystem_C_and_G_()
// Purpose       :
// Special Notes : This is only needed if device parameters are being swept.
//
//                 For sweeps, this only works for linear problems!
//
//                 For sensitivity analysis, this works for nonlinear 
//                 problems and devices.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/10/2018
//-----------------------------------------------------------------------------
bool AC::updateLinearSystem_C_and_G_()
{
  analysisManager_.getDataStore()->daeQVectorPtr->putScalar(0.0);
  analysisManager_.getDataStore()->daeFVectorPtr->putScalar(0.0);
  analysisManager_.getDataStore()->daeBVectorPtr->putScalar(0.0);

  analysisManager_.getDataStore()->dFdxdVpVectorPtr->putScalar(0.0);
  analysisManager_.getDataStore()->dQdxdVpVectorPtr->putScalar(0.0);
  analysisManager_.getDataStore()->dQdxMatrixPtr->put(0.0);
  analysisManager_.getDataStore()->dFdxMatrixPtr->put(0.0);

  loader_.updateState
                ((analysisManager_.getDataStore()->nextSolutionPtr),
               (analysisManager_.getDataStore()->currSolutionPtr),
               (analysisManager_.getDataStore()->lastSolutionPtr),
               (analysisManager_.getDataStore()->nextStatePtr),
               (analysisManager_.getDataStore()->currStatePtr),
               (analysisManager_.getDataStore()->lastStatePtr),
               (analysisManager_.getDataStore()->nextStorePtr),
               (analysisManager_.getDataStore()->currStorePtr),
               (analysisManager_.getDataStore()->lastStorePtr)
               );

  loader_.loadDAEVectors
              ((analysisManager_.getDataStore()->nextSolutionPtr),
               (analysisManager_.getDataStore()->currSolutionPtr),
               (analysisManager_.getDataStore()->lastSolutionPtr),
               (analysisManager_.getDataStore()->nextStatePtr),
               (analysisManager_.getDataStore()->currStatePtr),
               (analysisManager_.getDataStore()->lastStatePtr),
               (analysisManager_.getDataStore()->nextStateDerivPtr),
               (analysisManager_.getDataStore()->nextStorePtr),
               (analysisManager_.getDataStore()->currStorePtr),
               (analysisManager_.getDataStore()->lastStorePtr),
               (analysisManager_.getDataStore()->nextLeadCurrentPtr),
               (analysisManager_.getDataStore()->nextLeadCurrentQPtr),
               (analysisManager_.getDataStore()->nextLeadDeltaVPtr),
               (analysisManager_.getDataStore()->daeQVectorPtr),
               (analysisManager_.getDataStore()->daeFVectorPtr),
               (analysisManager_.getDataStore()->daeBVectorPtr),
               (analysisManager_.getDataStore()->dFdxdVpVectorPtr),
               (analysisManager_.getDataStore()->dQdxdVpVectorPtr) );

  loader_.loadDAEMatrices(analysisManager_.getDataStore()->nextSolutionPtr,
      analysisManager_.getDataStore()->nextStatePtr, analysisManager_.getDataStore()->nextStateDerivPtr,
      analysisManager_.getDataStore()->nextStorePtr,
      analysisManager_.getDataStore()->dQdxMatrixPtr,  analysisManager_.getDataStore()->dFdxMatrixPtr);

  C_ = analysisManager_.getDataStore()->dQdxMatrixPtr;
  G_ = analysisManager_.getDataStore()->dFdxMatrixPtr;

  if (false)
  {
  if (DEBUG_TIME && isActive(Diag::TIME_PARAMETERS))
  {
    Xyce::dout() << "dQdxMatrixPtr:" << std::endl;
    analysisManager_.getDataStore()->dQdxMatrixPtr->print( Xyce::dout() );

    Xyce::dout() << "dFdxMatrixPtr:" << std::endl;
    analysisManager_.getDataStore()->dFdxMatrixPtr->print( Xyce::dout() );

    Xyce::dout() << std::endl;
  }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::updateLinearSystemFreq_()
// Purpose       :
// Special Notes : ERK.  This has been modified to also reload G, to accomodate
//                 model parameter sweeps on linear circuits.
//
// Scope         : public
// Creator       : Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/20/2011
//-----------------------------------------------------------------------------
bool AC::updateLinearSystemFreq_()
{
  // First diagonal block
  ACMatrix_->put( 0.0 ); // Zero out whole matrix
  ACMatrix_->block( 0, 0 ).add(*G_);

  // Second diagonal block
  ACMatrix_->block( 1, 1 ).add(*G_);

  double omega =  2.0 * M_PI * currentFreq_;

  ACMatrix_->block( 0, 1).put( 0.0);
  ACMatrix_->block( 0, 1).add(*C_);
  ACMatrix_->block( 0, 1).scale(-omega);

  ACMatrix_->block(1, 0).put( 0.0);
  ACMatrix_->block(1, 0).add(*C_);
  ACMatrix_->block(1, 0).scale(omega);

  // Copy the values loaded into the blocks into the global matrix for the solve.
  ACMatrix_->assembleGlobalMatrix();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::updateLinearSystemMagAndPhase_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/7/2018
//-----------------------------------------------------------------------------
bool AC::updateLinearSystemMagAndPhase_()
{
  // ACMAG and ACPHASE have already been updated, so all that
  // remains is to load the B-vector and apply to the block system
  bVecRealPtr->putScalar(0.0);
  bVecImagPtr->putScalar(0.0);

  // re-load the B-vectors
  loader_.loadBVectorsforAC (bVecRealPtr, bVecImagPtr);

  B_->block( 0 ).update( 1.0, *bVecRealPtr, 0.0 );
  B_->block( 1 ).update( 1.0, *bVecImagPtr, 0.0 );

  return true;
}


//-----------------------------------------------------------------------------
// Function      : AC::solveLinearSystem_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :  Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/2011
//-----------------------------------------------------------------------------
bool AC::solveLinearSystem_()
{
  bool bsuccess = true;
 // Solve the block problem

  if (sparcalc_   == 0 )
  {
    int linearStatus = blockSolver_->solve();

    if (linearStatus != 0)
    {
      bsuccess = false;
    }
  }
  else
  // Loop over number of I/O ports here
  {

    Parallel::Manager * pdsMgrPtr = analysisManager_.getPDSManager();

    int myPID = pdsMgrPtr->getPDSComm()->procID();
    Parallel::Machine comm =  pdsMgrPtr->getPDSComm()->comm();

    Yparams_.putScalar(std::complex<double>(0.0, 0.0));

    for (unsigned int j=0; j < numPorts_; ++j)
    {

      B_->putScalar( 0.0 );

      if ( portPID_[j] == myPID )
        (B_->block( 0 ))[portMap_[j+1].first] = 1.0;

      int linearStatus = blockSolver_->solve();


      if (linearStatus != 0)
      {
  //      Xyce::dout() << "Amesos solve exited with error: " << linearStatus;
        bsuccess = false;
      }

      // Compute Y params entries for all I/O
      for (unsigned int i=0; i < numPorts_; ++i)
      {
         // Populate Y for all ports in L
         // L is the same as B and also a set of canonical basis vectors (e_i), so
         // we can pick off the appropriate entries of REFXPtr to place into X.
         if ( portPID_[i] == myPID )
           Yparams_(i,j) = std::complex<double>(-( X_->block( 0 ))[ portMap_[i+1].first],
                                            -( X_->block( 1 ))[portMap_[i+1].first]);  
      }
    }

    Xyce::Parallel::AllReduce(comm, MPI_SUM, Yparams_.values(), numPorts_*numPorts_);

  }

  return bsuccess;
}

// sensitivity functions
//-----------------------------------------------------------------------------
// Function      : AC::precomputeDCsensitivities_
//
// Purpose       : The goal of this function is to perform computations in 
//                 support of computing the Jacobian derivative, dJdp.
//
//                 Most of the operations only have to be performed once.
//
//                 To do this requires solving DCOP sensitivities for every p, 
//                 and using these to compute dGdp and dCdp for every p.  These 
//                 are independent of things like frequency, so for many 
//                 problems they only need to be computed once.
//
//
// (1) compute dx/dp direct sensitivity at the DCOP.  Using the direct method, 
//     you get the entire dx/dp vector for each param.  We need entire X vector, 
//     meaning adjoints aren't a realistic option for this.  
//
// (2) then, once you have this, perturb the DC x-vector using 
//          x_pert = x_orig + dx/dp * dP, where dx/dp is the DC sensitivity for p.
//
// (3) perturb p in the device package using : p_new = p + dP
//
// (3) re-load dFdx and dQdx this is the accurately perturbed dFdx and dQdx
//
// (4) then compute dCdp and dGdp via finite difference.
//
// (5) assemble dJdp using dCdp and dGdp.  However, this assembled dJdp will 
//     intentionally lack the omega factor in the C blocks.  As the 
//     calculation progresses thru a range of frequency, those blocks will be 
//     scaled and rescaled, to get the correct dJdp for that freq.  As long 
//     as we are sweeping frequency (and not model parameters), then this 
//     update to omega is all that is needed, and that is much cheaper than 
//     re-assembling everything.
//
//     If we ever want sensitivities to work with model parameter sweeps, 
//     then the full dJdp matrix has to be re-assembled each time those
//     parameter values are updated.  As of this writing the code cannot do that.
//     But sweeping parameters in an AC calculation isn't done very much.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 02/15/2022
//-----------------------------------------------------------------------------
bool AC::precomputeDCsensitivities_ ()
{
  Stats::StatTop _ACsensStat("AC sensitivities pre-compute");
  Stats::TimeBlock _ACsens_Timer(_ACsensStat);

  TimeIntg::DataStore & ds = *(analysisManager_.getDataStore());

  // do initial parameter setup
  bool origParamsDone =  Nonlinear::setupOriginalParams (
      ds, analysisManager_.getNonlinearEquationLoader(),
      paramNameVec_, analysisManager_);

  double epsilon = fabs(Util::MachineDependentParams::MachineEpsilon());
  double sqrtEta= std::sqrt(epsilon);

  param_dP_.clear();
  numericalDiff_.clear();
  numericalDiff_.resize(ds.paramOrigVals_.size(),0);

  for (int ii=0;ii< ds.paramOrigVals_.size();ii++)
  {
    const std::string name = paramNameVec_[ii];
    double origParamValue = ds.paramOrigVals_[ii];
    param_dP_.push_back(sqrtEta * fabs( origParamValue ));
  }

  {
  analysisManager_.getDataStore()->allocateSensitivityArrays(paramNameVec_.size(), true, false);
  int difference_ = 0;
  std::string netlist = analysisManager_.getNetlistFilename();
  Xyce::Nonlinear::loadSensitivityResiduals (difference_, 
      false, false, false,
      //forceFD_, forceDeviceFD_, forceAnalytic_, 
      sqrtEta, netlist, ds,
      analysisManager_.getNonlinearEquationLoader(),
      paramNameVec_,
      analysisManager_);

  analysisManager_.getNonlinearEquationLoader().loadSensitivityResiduals();
  }

  // compute dJdp for each p, but leave out the omega scalar from the C-blocks.
  //
  // If the analytical matrix sensitivities are available (unlikely), then use those.
  // otherwise compute a numerical dJdp.

  // Attempt to compute J' analytically first.   
  // Check if parameter is a B-term.  
  // If so we don't have to do anything.
  bool numerical_dJdp_Needed=false;
  for (int ipar=0;ipar<numSensParams_;++ipar)
  {
    const std::string name = paramNameVec_[ipar];
    bool analyticBVecSensAvailable = loader_.analyticBVecSensAvailable (name);
    bool deviceLevelBVecSensNumericalAvailable = loader_.numericalBVecSensAvailable (name);

    // check for B' sensitivities:
    std::vector< std::complex<double> > dbdp;
    if (analyticBVecSensAvailable && !forceDeviceFD_) { numericalDiff_[ipar] = 0; dBdpVector_[ipar]->putScalar(0.0); }
    else if (deviceLevelBVecSensNumericalAvailable && !forceAnalytic_) { numericalDiff_[ipar] = 1;  dBdpVector_[ipar]->putScalar(0.0); }

    if (!analyticBVecSensAvailable && !deviceLevelBVecSensNumericalAvailable)
    {
      bool analyticMatrixAvailable = loader_.analyticMatrixSensitivitiesAvailable (name);
      bool deviceLevelMatrixNumericalAvailable = loader_.numericalMatrixSensitivitiesAvailable (name);

      std::vector< std::vector<double> > d_dfdx_dp, d_dqdx_dp;
      std::vector< std::vector<int> > F_jacLIDs, Q_jacLIDs;
      std::vector<int> F_lids, Q_lids;

      if (analyticMatrixAvailable && !forceDeviceFD_)
      {
        numericalDiff_[ipar] = 0;

        loader_.getAnalyticMatrixSensitivities(name,
                      d_dfdx_dp, d_dqdx_dp, F_lids, Q_lids, F_jacLIDs, Q_jacLIDs);

        dGdp_->put(0.0);
        for (int ii=0;ii<F_lids.size();++ii)
        {
          for (int jj=0;jj<F_jacLIDs[ii].size();++jj)
          {
            // NOTE: this used to be a sum, ie, "+=".  
            // Now it is not, as we are doing one param at a time, so there 
            // shouldn't be more than one legitimate contribution to the 
            // dGdp matrix.  However, the d_dfdx_dp object may have some
            // duplicates in it.  So, should not sum 2x.
            (*dGdp_)[F_lids[ii]][F_jacLIDs[ii][jj]] = d_dfdx_dp[ii][jj];
          }
        }
        dGdp_->fillComplete();

        dCdp_->put(0.0);
        for (int ii=0;ii<Q_lids.size();++ii)
        {
          for (int jj=0;jj<Q_jacLIDs[ii].size();++jj)
          {
            // NOTE: this used to be a sum, ie, "+=".  
            // Now it is not, as we are doing one param at a time, so there 
            // shouldn't be more than one legitimate contribution to the 
            // dCdp matrix.  However, the d_dqdx_dp object may have some
            // duplicates in it.  So, should not sum 2x.
            (*dCdp_)[Q_lids[ii]][Q_jacLIDs[ii][jj]] = d_dqdx_dp[ii][jj];
          }
        }
        dCdp_->fillComplete();

        dJdpVector_[ipar]->put(0.0);// Zero out whole matrix
        dJdpVector_[ipar]->block(0, 0).add(*dGdp_);
        dJdpVector_[ipar]->block(1, 1).add(*dGdp_);
        dJdpVector_[ipar]->block(0, 1).add(*dCdp_);
        dJdpVector_[ipar]->block(1, 0).add(*dCdp_);
        dJdpVector_[ipar]->assembleGlobalMatrix();

        numericalDiff_[ipar] = 0;
      }
      else
      {
        numericalDiff_[ipar] = 1;
        numerical_dJdp_Needed=true;
      }
    }
  }

  if (numerical_dJdp_Needed)
  {
    // save copies of the original C and G DCOP matrices.  
    origC_->put(0.0);
    origC_->addOverlap(*C_);
    origC_->fillComplete();
    origG_->put(0.0);
    origG_->addOverlap(*G_);
    origG_->fillComplete();

    // save original B vector:
    *origB_ = *B_;

    // compute dxdp for each p, and scale by dP to get dx
    for (int ipar=0;ipar<numSensParams_;++ipar)
    {
      if (numericalDiff_[ipar] == 0) continue;

      Linear::Vector *dcDXdp = linearSystem_.getNewtonVector();
      Linear::Vector* rhsVectorPtr_ = linearSystem_.getRHSVector();
      Linear::MultiVector * sensRHSPtrVector = ds.sensRHSPtrVector;

      *rhsVectorPtr_ = *Teuchos::rcp(sensRHSPtrVector->getNonConstVectorView(ipar));

      Teuchos::RCP<Linear::Solver> sensSolver = nonlinearManager_.getNonlinearSolver().getLinearSolver();
      bool reuseFactors=true;
      sensSolver->solve(reuseFactors);

      // Use dxdp to get dx = dxdp*dP 
      dcDXdp->scale(param_dP_[ipar]);

      // save next solution
      analysisManager_.getDataStore()->tmpSolVectorPtr->update(1.0, *(analysisManager_.getDataStore()->nextSolutionPtr), 0.0);
      // update next Solution to the perturbed value of x+dx
      analysisManager_.getDataStore()->nextSolutionPtr->update(1.0, *(dcDXdp), 1.0);

      // perturb p to p+dP
      const std::string name = paramNameVec_[ipar];
      double origParamValue = ds.paramOrigVals_[ipar];
      double dP = param_dP_[ipar];
      double newParamValue = origParamValue + dP;
      std::string tmpName = name; // remove const
      loader_.setParam(tmpName, newParamValue, true);

      // reload C and G, with updated x and p
      updateLinearSystem_C_and_G_();

      // compute numerical B' = dBdp = (perturbB - origB)/dP
      // reload B vector, sett up perturbed B
      bVecRealPtr->putScalar(0.0);
      bVecImagPtr->putScalar(0.0);
      loader_.loadBVectorsforAC (bVecRealPtr, bVecImagPtr);
      B_->block( 0 ).update( 1.0, *bVecRealPtr, 0.0 );
      B_->block( 1 ).update( 1.0, *bVecImagPtr, 0.0 );

      dBdpVector_[ipar]->update(1.0, *B_, -1.0, *origB_, 0.0 );
      dBdpVector_[ipar]->scale(1.0/dP);

      // compute numerical dGdp 
      // dGdp = (perturbG - origG)/dP
      dGdp_->put(0.0);
      dGdp_->addOverlap(*origG_);
      dGdp_->scale(-1.0);
      dGdp_->addOverlap(*G_);
      dGdp_->scale(1.0/dP);
      dGdp_->fillComplete();

      // compute numerical dCdp 
      // dCdp = (perturbC - origC)/dP
      dCdp_->put(0.0);
      dCdp_->addOverlap(*origC_);
      dCdp_->scale(-1.0);
      dCdp_->addOverlap(*C_);
      dCdp_->scale(1.0/dP);
      dCdp_->fillComplete();

      // put together the whole thing and save.
      // This is an incomplete dJdp.  It (intentionally) lacks the omega scalar on the C blocks.
      dJdpVector_[ipar]->put(0.0);// Zero out whole matrix
      dJdpVector_[ipar]->block(0, 0).add(*dGdp_);
      dJdpVector_[ipar]->block(1, 1).add(*dGdp_);
      dJdpVector_[ipar]->block(0, 1).add(*dCdp_);
      dJdpVector_[ipar]->block(1, 0).add(*dCdp_);
      dJdpVector_[ipar]->assembleGlobalMatrix();

      // restore original p and x.
      loader_.setParam(tmpName, origParamValue, true);
      analysisManager_.getDataStore()->nextSolutionPtr->update(1.0, *(analysisManager_.getDataStore()->tmpSolVectorPtr), 0.0);
    }

    // restore the original pre-sensitivity C and G
    C_->put(0.0);
    C_->addOverlap(*origC_);
    C_->fillComplete();
    G_->put(0.0);
    G_->addOverlap(*origG_);
    G_->fillComplete();

    // restore original B vector:
    *B_ = *origB_;

  }

  // at this point everything should be restored.   Not sure if this call is necessary.  check later
  updateLinearSystem_C_and_G_();

  // these are not needed anymore.  But if we ever upate AC .sens to work with 
  // model parameter sweeps, then this function could be called multiple times, 
  // and then they would be needed throughout the calculation.
  delete dCdp_;
  delete dGdp_;
  delete origC_;
  delete origG_;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::solveSensitivity_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AC::solveSensitivity_()
{
  TimeIntg::DataStore & ds = *(analysisManager_.getDataStore());

  ds.dOdpVec_.clear();
  ds.dOdpAdjVec_.clear();
  ds.scaled_dOdpVec_.clear();
  ds.scaled_dOdpAdjVec_.clear();

  Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();

  if ( objFuncGiven_ )
  {
    std::string netlist = analysisManager_.getNetlistFilename();

    Parallel::Machine comm =  analysisManager_.getPDSManager()->getPDSComm()->comm();
    Xyce::Nonlinear::evaluateObjFuncs (
        objFuncDataVec_, comm,
        analysisManager_.getNonlinearEquationLoader(), 
        netlist); 

    objectiveVec_.clear();
    for (int iobj=0;iobj<objFuncDataVec_.size();++iobj)
    {
      double xr = std::real(objFuncDataVec_[iobj]->expVal);
      double xi = std::imag(objFuncDataVec_[iobj]->expVal);
      double sumOfSquares = (xr*xr + xi*xi);
      double xm = sqrt(sumOfSquares);
      double xp = (std::arg(objFuncDataVec_[iobj]->expVal)); 
      if (!outputManagerAdapter_.getPhaseOutputUsesRadians())
        xp *= 180.0/M_PI;

      objectiveVec_.push_back(xr);
      objectiveVec_.push_back(xi);
      objectiveVec_.push_back(xm);
      objectiveVec_.push_back(xp);
    }
  }

  if(solveDirectSensitivityFlag_)
  {
    Stats::StatTop _ACsensStat("AC Direct Sensitivity");
    Stats::TimeBlock _ACsens_Timer(_ACsensStat);

    bool directSuccess = solveDirectSensitivity_();
  }

  if(solveAdjointSensitivityFlag_)
  {
    Stats::StatTop _ACsensStat("AC Adjoint Sensitivity");
    Stats::TimeBlock _ACsens_Timer(_ACsensStat);

    bool adjointSuccess = solveAdjointSensitivity_();
  }


  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::solve_mag_phase_Sensitivities
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter, SNL
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
void AC::solve_mag_phase_Sensitivities_(
      const double dxdpReal,
      const double dxdpImag,
      const double xr,
      const double xi,

      double & dxdp_mag,
      double & dxdp_phase
      )
{
  // --------
  // magnitude:  vm = sqrt(vr^2 + vi^2)
  //  dvm/dvr = 1/2 * (vr^2+vi^2)^(-1/2)  *  (2*vr)
  //  dvm/dvi = 1/2 * (vr^2+vi^2)^(-1/2)  *  (2*vi)
  double sumOfSquares = (xr*xr + xi*xi);
  double mag = sqrt(sumOfSquares);
  double dxmag_dxr = 0.0;
  double dxmag_dxi = 0.0;

  if (mag != 0.0)
  {
    dxmag_dxr = xr/mag;
    dxmag_dxi = xi/mag;
  }

  dxdp_mag   = (dxmag_dxr   * dxdpReal + dxmag_dxi   * dxdpImag);

  // --------
  // phase = atan(vi/vr)
  // derivative of atan(vi/vr) w.r.t. vr = - vi/(vi^2 + vr^2)
  // derivative of atan(vi/vr) w.r.t. vi =   vr/(vi^2 + vr^2)
  //double phase = atan(xi/xr);
  double dxphase_dxr = 0.0;
  double dxphase_dxi = 0.0;
  if (sumOfSquares != 0.0)
  {
    dxphase_dxr = - xi / sumOfSquares;
    dxphase_dxi =   xr / sumOfSquares;
  }

  dxdp_phase = (dxphase_dxr * dxdpReal + dxphase_dxi * dxdpImag);
  // account for whether the phase output is in radians or degrees
  if (!outputManagerAdapter_.getPhaseOutputUsesRadians())
    dxdp_phase *= 180.0/M_PI;
}

//-----------------------------------------------------------------------------
// Function      : sensStdOutput
// Purpose       : output sensitivity stuff to the screen
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/28/2019
//-----------------------------------------------------------------------------
std::ostream& sensStdOutput (
       const std::string idString,
       const std::vector<double> & paramVals,
       const std::vector<double> & sensitivities,
       const std::vector<double> & scaled_sensitivities,
       const std::vector<std::string> & paramNameVec,
       const std::vector<double> & param_dP,
       const std::vector<int> & numericalDiff,
       const std::vector<Xyce::Nonlinear::objectiveFunctionData<std::complex<double> > *> & objFuncDataVec,
       const std::vector<double> & objectiveVec,
       OutputMgrAdapter & outputManagerAdapter,
       std::ostream& os
       )
{
  // save current stream state, and then set the stream to use scientific notation.
  // Otherwise the info for the stepped parameters may not be output in
  // scientific notation.
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os.setf(std::ios::scientific);

  if ( !(outputManagerAdapter.getStepSweepVector().empty()) )
  {
    // Output for .STEP.  Step count should be output as an integer.
    // The other values should be output in scientific notation.
    os << "\nFor Step " << outputManagerAdapter.getStepAnalysisStepNumber() << ":" << std::endl;
    for (std::vector<Analysis::SweepParam>::const_iterator it = outputManagerAdapter.getStepSweepVector().begin();
         it != outputManagerAdapter.getStepSweepVector().end(); ++it)
    {
      os << it->name << " = " << it->currentVal << std::endl;
    }
  }

  int maxSize=5;
  for (int ip=0;ip<paramVals.size();++ip){if(paramNameVec[ip].size() > maxSize) maxSize = paramNameVec[ip].size();}

  int numWV=11;
  int numW=14;
  int numWP=18;
  int ii=0;
  int sensIndex=0;
  if (!(objFuncDataVec.empty()))
  {
    int ii=0;
    int sensIndex=0;
    for (int iobj=0;iobj<objFuncDataVec.size();++iobj)
    {
      double xr = std::real(objFuncDataVec[iobj]->expVal);
      double xi = std::imag(objFuncDataVec[iobj]->expVal);
      double sumOfSquares = (xr*xr + xi*xi);
      double xm = sqrt(sumOfSquares);
      double xp = (std::arg(objFuncDataVec[iobj]->expVal)); 
      if (!outputManagerAdapter.getPhaseOutputUsesRadians())
        xp *= 180.0/M_PI;

      os << "\n"<<idString << " Sensitivities for "<< objFuncDataVec[iobj]->objFuncString <<std::endl;

      os << " Re(" << objFuncDataVec[iobj]->objFuncString << ") = " 
        << std::setw(numW)<< std::scientific<< std::setprecision(4) 
        << xr << "  ";

      os << "Img(" << objFuncDataVec[iobj]->objFuncString << ") = " 
        << std::setw(numW)<< std::scientific<< std::setprecision(4) 
        << xi <<std::endl;

      if (!outputManagerAdapter.getPhaseOutputUsesRadians()){ xp *= 180.0/M_PI; }

      os << "  M(" << objFuncDataVec[iobj]->objFuncString << ") = " 
        << std::setw(numW)<< std::scientific<< std::setprecision(4) 
        << xm << "  ";

      os << " Ph(" << objFuncDataVec[iobj]->objFuncString << ") = " 
        << std::setw(numW)<< std::scientific<< std::setprecision(4) 
        << xp << std::endl;

      os << std::setw(maxSize)<<std::left<< "Name";
      os << "\t"<<std::setw(numWV)<<"Value";
      os << "\t"<<std::setw(numW)<<"Sensitivity_Re";
      os << "\t"<<std::setw(numW)<<"Sensitivity_Im";
      os << "\t"<<std::setw(numW)<<"Sensitivity_Mag";
      os << "\t"<<std::setw(numWP)<<"Sensitivity_Phase";
      os << "\t"<<std::setw(numWV)<<"Delta P";
      os << std::endl;

      for (int iparam=0; iparam< paramVals.size(); ++iparam)
      {
        os << std::setw(maxSize)<<std::left<<paramNameVec[iparam];
        os <<"\t"<< std::setw(numWV)<< std::scientific<< std::setprecision(4) << paramVals[iparam];
        os <<"\t"<< std::setw(numW)<< std::scientific<< std::setprecision(4)<< sensitivities[sensIndex++];
        os <<"\t"<< std::setw(numW)<< std::scientific<< std::setprecision(4)<< sensitivities[sensIndex++];
        os <<"\t"<< std::setw(numW)<< std::scientific<< std::setprecision(4)<< sensitivities[sensIndex++];
        os <<"\t"<< std::setw(numWP)<< std::scientific<< std::setprecision(4)<< sensitivities[sensIndex++];
        os <<"\t"<< std::setw(numWV)<< std::scientific<< std::setprecision(4)<< param_dP[iparam];
        if (numericalDiff[iparam]==1) { os << "\t" << "FD used"; }
        else { os << "\t" << "FD not used";}
        os << std::endl;
      }
    }
  }

  if ( !(outputManagerAdapter.getStepSweepVector().empty()) )
  {
    if ( (outputManagerAdapter.getStepAnalysisStepNumber()+1) <
          outputManagerAdapter.getStepAnalysisMaxSteps() )
    {
      // add a blank line after each block of .STEP output, except for the last one
      os << std::endl;
    }
  }

  // previous state of os is restored when the function returns
  return os;
}

//-----------------------------------------------------------------------------
// Function      : AC::solveDirectSensitivity_()
// Purpose       : Solves parameter sensitivities using the direct method
//
//    The AC analysis calculation is a linear calculation which solves:
//
//    J * x = b, where J is the Jacobian matrix from the DCOP, but
//
//    converted to frequency domain.  b is the source vector, and x
//    is the solution vector.
//
//    We want x'  where x' = \partial x/ \partial p
//
//    Take derivative of both sides:
//
//    J' * x + J * x' = b'
//
//    giving:
//
//    J * x' = b' - J' * x
//
//    So we can obtain x' by solving the same matrix, but with a different RHS.
//
//    Unfortunately, for the direct method this must be done for each p.
//
// Special Notes : 
//
// Each solution output has a real and imaginary part and
// each of those parts is treated as a separate objective in part of this analysis.
// For the direct method that means that for every dO/dp = dO/dx * dx/dp dot 
// product, two of them are performed.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/21/2019
//-----------------------------------------------------------------------------
bool AC::solveDirectSensitivity_()
{
  TimeIntg::DataStore & ds = *(analysisManager_.getDataStore());
  Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();
  Parallel::Communicator & comm = *(pds_manager.getPDSComm());
  int myPID = comm.procID();

  // list of objective functions specified as expressions.
  {
    int ii=0;
    int numOutFuncs = objFuncDataVec_.size();
    for (int iobj=0;iobj<numOutFuncs;++iobj)
    {
      double xr=std::real(objFuncDataVec_[iobj]->expVal);
      double xi=std::imag(objFuncDataVec_[iobj]->expVal);

      for (int ipar=0;ipar<numSensParams_;++ipar)
      {
        bool rhsSuccess = loadSensitivityRHS_(ipar);

        savedX_->update( 1.0, *X_, 0.0 );

        B_->update( 1.0, *sensRhs_, 0.0 );

        int linearStatus = blockSolver_->solve(reuseFactors_);

        dXdp_->update( 1.0, *X_, 0.0 );
        X_->update( 1.0, *savedX_, 0.0 ); // restore X

        double dOdp_real = objFuncDataVec_[iobj]->dOdXVectorRealPtr->dotProduct ( dXdp_->block(0) );
        double dOdp_imag = objFuncDataVec_[iobj]->dOdXVectorImagPtr->dotProduct ( dXdp_->block(1) );
        std::complex<double> dOdp = std::complex<double>(dOdp_real,dOdp_imag);

        ds.dOdpVec_.push_back(std::real(dOdp));
        ds.dOdpVec_.push_back(std::imag(dOdp));
        
        // do polar
        double dOdp_mag, dOdp_phase;
        solve_mag_phase_Sensitivities_(std::real(dOdp), std::imag(dOdp), xr, xi, dOdp_mag, dOdp_phase);
        ds.dOdpVec_.push_back(dOdp_mag);
        ds.dOdpVec_.push_back(dOdp_phase);
      }
    }
  }

  if (stdOutputFlag_)
  {
    Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();

    Analysis::sensStdOutput(std::string("Direct"), ds.paramOrigVals_, ds.dOdpVec_, ds.scaled_dOdpVec_,
       paramNameVec_, 
       param_dP_,
       numericalDiff_,
       objFuncDataVec_,
       objectiveVec_,
       outputManagerAdapter_, Xyce::lout ());
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::solveAdjointSensitivity_()
// Purpose       : solves sensitivities using the adjoint method
// Special Notes :
//
// Each solution output has a real and imaginary part and
// each of those parts is treated as a separate objective in part of this analysis.
// For the direct method that means that for every J^T * lambda = dOdX linear solve,
// product, two of them are performed, one for dOdX_real and another for dOdX_imag.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/21/2019
//-----------------------------------------------------------------------------
bool AC::solveAdjointSensitivity_()
{
  TimeIntg::DataStore & ds = *(analysisManager_.getDataStore());
  Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();
  Parallel::Communicator & comm = *(pds_manager.getPDSComm());
  int myPID = comm.procID();

  // list of objective functions specified as expressions.
  {
    int ii=0;
    int numOutFuncs = objFuncDataVec_.size();
    for (int iobj=0;iobj<numOutFuncs;++iobj)
    {
      dOdXreal_->block( 0 ).update( 1.0, *(objFuncDataVec_[iobj]->dOdXVectorRealPtr), 0.0 );
      dOdXreal_->block( 1 ).putScalar( 0.0 );

      dOdXimag_->block( 0 ).putScalar( 0.0 );
      dOdXimag_->block( 1 ).update( 1.0, *(objFuncDataVec_[iobj]->dOdXVectorImagPtr), 0.0 );

      double xr=std::real(objFuncDataVec_[iobj]->expVal);
      double xi=std::imag(objFuncDataVec_[iobj]->expVal);

      savedX_->update( 1.0, *X_, 0.0 );

      B_->update( 1.0, *dOdXreal_, 0.0 );

      int linearStatus= blockSolver_->solveTranspose(reuseFactors_); 

      lambda_->update( 1.0, *X_, 0.0 );

      X_->update( 1.0, *savedX_, 0.0 ); // restore X
      
      std::vector<double> dOdpReal(numSensParams_,0.0);
      std::vector<double> dOdpImag(numSensParams_,0.0);
      std::vector<double> dOdpMag (numSensParams_,0.0);
      std::vector<double> dOdpPhase(numSensParams_,0.0);

      for (int ipar=0;ipar<numSensParams_;++ipar)
      {
        bool rhsSuccess = loadSensitivityRHS_(ipar);
        // compute dot products for this param
        dOdpReal[ipar] = sensRhs_->dotProduct( *lambda_ );
      }

      savedX_->update( 1.0, *X_, 0.0 );

      B_->update( 1.0, *dOdXimag_, 0.0 );

      linearStatus= blockSolver_->solveTranspose(reuseFactors_); 

      lambda_->update( 1.0, *X_, 0.0 );

      X_->update( 1.0, *savedX_, 0.0 ); // restore X
      for (int ipar=0;ipar<numSensParams_;++ipar)
      {
        bool rhsSuccess = loadSensitivityRHS_(ipar);
        // compute dot products for this param
        dOdpImag[ipar] = sensRhs_->dotProduct( *lambda_ );
      }

      // get mag, phase sensitivities:
      for (int ipar=0;ipar<numSensParams_;++ipar)
      {
        solve_mag_phase_Sensitivities_(dOdpReal[ipar], dOdpImag[ipar], xr, xi, dOdpMag[ipar], dOdpPhase[ipar]);
      }

      for (int ipar=0;ipar<numSensParams_;++ipar)
      {
        ds.dOdpAdjVec_.push_back(dOdpReal[ipar]);
        ds.dOdpAdjVec_.push_back(dOdpImag[ipar]);
        ds.dOdpAdjVec_.push_back(dOdpMag[ipar]);
        ds.dOdpAdjVec_.push_back(dOdpPhase[ipar]);
      }
    }
  }

  if (stdOutputFlag_)
  {
    Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();

    Analysis::sensStdOutput (std::string("Adjoint"), ds.paramOrigVals_, ds.dOdpAdjVec_, ds.scaled_dOdpAdjVec_,
       paramNameVec_, 
       param_dP_,
       numericalDiff_,
       objFuncDataVec_,
       objectiveVec_,
       outputManagerAdapter_, Xyce::lout ());
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::applyOmega_dJdp
//
// Purpose       : scales the C-blocks by omega=2*pi*freq
//
// Special Notes : applies matvec without matrix assembly
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 3/02/2022
//-----------------------------------------------------------------------------
void AC::applyOmega_dJdp(bool transA, const Linear::BlockMatrix * A, 
                         const Linear::BlockVector& x, Linear::BlockVector& y)
{
  double omega =  2.0 * M_PI * currentFreq_;
  Linear::Vector* tmpVec = y.block(0).cloneVector(); 

  // A is scaled by omega in the off-diagonal blocks
  // block(0, 1).scale(-omega);
  // block(1, 0).scale(omega);

  if (transA)
  {
    // y[0] = A[0,0]'*x[0] + omega*A[1,0]'*x[1];
    (*A).block(0, 0).matvec( transA, x.block(0), *tmpVec );
    (*A).block(1, 0).matvec( transA, x.block(1), y.block(0) );
    y.block(0).update( 1.0, *tmpVec, omega );

    // y[1] = A[1,1]'*x[1] - omega*A[0,1]'*x[0]
    (*A).block(1, 1).matvec( transA, x.block(1), *tmpVec );
    (*A).block(0, 1).matvec( transA, x.block(0), y.block(1) );
    y.block(1).update( 1.0, *tmpVec, -omega );
  }
  else
  {
    // y[0] = A[0,0]*x[0] - omega*A[0,1]*x[1];
    (*A).block(0, 0).matvec( transA, x.block(0), *tmpVec );
    (*A).block(0, 1).matvec( transA, x.block(1), y.block(0) );
    y.block(0).update( 1.0, *tmpVec, -omega );

    // y[1] = A[1,1]*x[1] + omega*A[1,0]*x[0];
    (*A).block(1, 1).matvec( transA, x.block(1), *tmpVec );
    (*A).block(1, 0).matvec( transA, x.block(0), y.block(1));
    y.block(1).update( 1.0, *tmpVec, omega );
  }

  delete tmpVec;
}

//-----------------------------------------------------------------------------
// Function      : AC::loadSensitivityRHS_()
// Purpose       : puts values into the rhs block vector.  The is the RHS for 
//                 the linear system in the direct method.  It is also used in 
//                 the adjoint method for the dot product that must be computed
//                 after the transpose linear solve.
//
//            RHS = B' - J' * X
//
// Special Notes :
//
// for most device parameters, there will either be a B' contribution, or 
// there will be a J'*X contribution.   Most parameters will not produce both, 
// unless the parameter is a global parameter, which applies to multiple 
// devices.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/15/2019
//-----------------------------------------------------------------------------
bool AC::loadSensitivityRHS_(int ipar)
{
  const std::string name = paramNameVec_[ipar];

  bool analyticBVecSensAvailable = loader_.analyticBVecSensAvailable (name);
  bool deviceLevelBVecSensNumericalAvailable = loader_.numericalBVecSensAvailable (name);

  // B' sensitivities:
  std::vector< std::complex<double> > dbdp;
  std::vector<int> BindicesVec;
  if (analyticBVecSensAvailable && !forceDeviceFD_)
  {
    loader_.getAnalyticBSensVectorsforAC (name, dbdp, BindicesVec);
    numericalDiff_[ipar] = 0;
  }
  else if (deviceLevelBVecSensNumericalAvailable && !forceAnalytic_)
  {
    loader_.getNumericalBSensVectorsforAC (name, dbdp, BindicesVec);
    numericalDiff_[ipar] = 1;
  }

  if (analyticBVecSensAvailable || deviceLevelBVecSensNumericalAvailable)
  {
    dbdpVecRealPtr->putScalar(0.0);
    dbdpVecImagPtr->putScalar(0.0);

    for (int ib=0;ib<BindicesVec.size();++ib)
    {
      (*dbdpVecRealPtr)[BindicesVec[ib]] += dbdp[ib].real();
      (*dbdpVecImagPtr)[BindicesVec[ib]] += dbdp[ib].imag();
    }

#ifdef Xyce_PARALLEL_MPI
    dbdpVecRealPtr->importOverlap();
    dbdpVecImagPtr->importOverlap();
#endif

    dBdp_->block( 0 ).update( 1.0, *dbdpVecRealPtr, 0.0 );
    dBdp_->block( 1 ).update( 1.0, *dbdpVecImagPtr, 0.0 );

    sensRhs_->update( 1.0, *dBdp_, 0.0 );  
  }
  else 
  // If couldn't obtain B', then try J'.   
  // Most stuff related to J' was pre-computed in the 
  // "AC::precomputeDCsensitivities_ ()" function. so not much work to be done here.
  // Just scale the C-block of the dJdp matrix by omega and then do a matvec.
  {
    // compute the matvec and then sum into the rhs vector.
    bool Transpose = false;
    applyOmega_dJdp( Transpose, dJdpVector_[ipar], *X_, *sensRhs_ );

    sensRhs_->update( 1.0, *dBdpVector_[ipar], -1.0 );  
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::outputSensitivity_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AC::outputSensitivity_()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::processSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AC::doProcessSuccessfulStep()
{
  if  (sparcalc_ == 0)
  {
    // Output x.
    outputManagerAdapter_.outputAC (currentFreq_, fStart_,fStop_,
	    X_->block(0), X_-> block(1), RFparams_);

    if (sensFlag_ && !(objFuncDataVec_.empty()))
    {
        // this outputs both the direct and adjoint sensitivity information
        const TimeIntg::DataStore & ds = *(analysisManager_.getDataStore());
        outputManagerAdapter_.outputSensitivityAC(currentFreq_, X_->block(0), X_-> block(1),
          ds.paramOrigVals_, paramNameVec_, objFuncStrings_, objectiveVec_,
          ds.dOdpVec_, ds.dOdpAdjVec_, ds.scaled_dOdpVec_, ds.scaled_dOdpAdjVec_);
    }
  }
  else
  {
    // only do these conversions if S- or Z-parameter output was requested on the
    // .LIN or .AC line
    if (sParamsRequested_) { Util::ytos(Yparams_, Sparams_, Z0sVec_ );}
    if (zParamsRequested_) { Util::ytoz(Yparams_, Zparams_); }

    // Outputter for Touchstone1 and/or Touchstone2 formatted files.
    // acLoopSize_ is the total number of frequency points in the analyses.
    outputManagerAdapter_.outputSParams(currentFreq_, acLoopSize_, Z0sVec_, RFparams_);

    outputManagerAdapter_.outputAC (currentFreq_, fStart_,fStop_,
            X_->block(0), X_-> block(1), RFparams_);

//    outputMOR_.output(outputManagerAdapter_.getComm(), 1, currentFreq_,  Yparams_ );
  }

  if ( !firstDoubleDCOPStep() )
  {
      stepNumber += 1;
      stats_.successStepsThisParameter_ += 1;
      stats_.successfulStepsTaken_ += 1;
  }

  // This output call is for device-specific output (like from a PDE device,
  // outputting mesh-based tecplot files).  It will only work in parallel if on
  // a machine where all processors have I/O capability.
  loader_.outputPlotFiles();

  return true;
}


//-----------------------------------------------------------------------------
// Function      : AC::processFailedStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AC::doProcessFailedStep()
{

  stepNumber += 1;

  acSweepFailures_.push_back(stepNumber);
  stats_.failedStepsAttempted_  += 1;
  analysisManager_.getStepErrorControl().numberSuccessiveFailures += 1;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::finish
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AC::doFinish()
{
  bool bsuccess = true;

  if (DEBUG_ANALYSIS)
  {
    Xyce::dout() << "Calling AC::doFinish() outputs!" << std::endl;
  }

  outputManagerAdapter_.finishOutput();
  if (!(acSweepFailures_.empty()))
  {
    bsuccess = false;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AC::doHandlePredictor
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 06/24/2013
//-----------------------------------------------------------------------------
bool AC::doHandlePredictor()
{
  analysisManager_.getDataStore()->setErrorWtVector(tiaParams_, topology_.getVarTypes());
  analysisManager_.getWorkingIntegrationMethod().obtainPredictor();
  analysisManager_.getWorkingIntegrationMethod().obtainPredictorDeriv();

  // In case this is the upper level of a 2-level sim, tell the
  // inner solve to do its prediction:
  bool   beginIntegrationFlag = analysisManager_.getBeginningIntegrationFlag();
  double nextTimeStep = analysisManager_.getStepErrorControl().currentTimeStep;
  double nextTime = analysisManager_.getStepErrorControl().nextTime;
  int    currentOrder = analysisManager_.getWorkingIntegrationMethod().getOrder();

  loader_.startTimeStep(beginIntegrationFlag, nextTimeStep, nextTime, currentOrder);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::updateDataParams_
// Purpose       :
// Special Notes : Used for AC analysis classes, when .DATA is used
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/10/2018
//-----------------------------------------------------------------------------
bool AC::updateDataParams_ (int stepNumber)
{
  bool reset = updateSweepParams(stepNumber, acSweepVector_.begin(), acSweepVector_.end());

  bool nonFreqPresent=false;
  for (int iac=0;iac<acSweepVector_.size();++iac)
  {
    std::string name = acSweepVector_[iac].name; Util::toUpper(name);
    double val = acSweepVector_[iac].currentVal;
    if (name == "FREQ" || name == "HERTZ")
    {
      currentFreq_ = val;
    }
    else
    {
      nonFreqPresent=true;
      loader_.setParam(name, val, true);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::updateCurrentFreq_(int stepNumber )
// Purpose       :
// Special Notes : Used for AC analysis classes.
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AC::updateCurrentFreq_(int stepNumber)
{
  if (type_ == "LIN")
  {

    currentFreq_  = fStart_ + static_cast<double>(stepNumber)*fstep_;
  }
  else if(type_ == "DEC" || type_ == "OCT")
  {

    currentFreq_ = fStart_*pow(stepMult_, static_cast<double>(stepNumber) );
  }
  else
  {
    Report::DevelFatal().in("AC::updateCurrentFreq_") << "AC::updateCurrentFreq_: unsupported STEP type";
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::setupSweepParam_
// Purpose       : Processes sweep parameters.
// Special Notes : Used for AC analysis classes.
// Scope         : public
// Creator       : Ting Mei, SNL.
// Creation Date :
//-----------------------------------------------------------------------------
int AC::setupSweepParam_()
{
  double fstart, fstop;
  double fcount = 0.0;

  fstart = fStart_;
  fstop = fStop_;

  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    Xyce::dout() << std::endl << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "AC::setupSweepParam_" << std::endl;
  }

    if (type_ == "LIN")
    {
      int np = static_cast<int>(np_);

      if ( np == 1)
        fstep_ = 0;
      else
        fstep_  = (fstop - fstart)/(np_ - 1.0);

      fcount = np_;
      if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
      {
        Xyce::dout() << "fstep   = " << fstep_  << std::endl;
      }
    }
    else if (type_ == "DEC")
    {
      stepMult_ = pow(10.0, 1.0/np_);
      fcount   = floor(fabs(log10(fstart) - log10(fstop)) * np_ + 1.0);
      if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
      {
        Xyce::dout() << "stepMult_ = " << stepMult_  << std::endl;
      }
    }
    else if (type_ == "OCT")
    {
      stepMult_ = pow(2.0, 1.0/np_);

      // changed to remove dependence on "log2" function, which apparently
      // doesn't exist in the math libraries of FreeBSD or the mingw
      // cross-compilation suite.   Log_2(x)=log_e(x)/log_e(2.0)
      double ln2 = log(2.0);
      fcount   = floor(fabs(log(fstart) - log(fstop))/ln2 * np_ + 1.0);
      if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
      {
        Xyce::dout() << "stepMult_ = " << stepMult_  << std::endl;
      }
    }
    else
    {
      Report::UserFatal0() << "Unsupported AC sweep type: " << type_;
    }

  // At this point, pinterval equals the total number of steps
  // for the step loop.
  return static_cast<int> (fcount);
}

namespace {

typedef Util::Factory<AnalysisBase, AC>  ACFactoryBase;

//-----------------------------------------------------------------------------
// Class         : ACFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jan 29 12:53:02 2015
//-----------------------------------------------------------------------------
///
/// Factory for parsing AC parameters from the netlist and creating AC analysis.
///
class ACFactory : public ACFactoryBase
{
public:
  //-----------------------------------------------------------------------------
  // Function      : ACFactory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:54:09 2015
  //-----------------------------------------------------------------------------
  ///
  /// Constructs the AC analysis factory
  ///
  /// @invariant Stores the results of parsing, so if more than one of the analysis and
  /// associated lines are parsed, the second options simply overwrite the previously parsed
  /// values.
  ///
  /// @invariant The existence of the parameters specified in the constructor cannot
  /// change.
  ///
  /// @param analysis_manager
  /// @param linear_system
  /// @param nonlinear_manager
  /// @param topology
  ///
  ACFactory(
    Analysis::AnalysisManager & analysis_manager,
    Linear::System &            linear_system,
    Nonlinear::Manager &        nonlinear_manager,
    Loader::Loader &            loader,
    Topo::Topology &            topology,
    Device::DeviceMgr &                 device_manager,
    IO::InitialConditionsManager &      initial_conditions_manager)
    : ACFactoryBase(),
      analysisManager_(analysis_manager),
      linearSystem_(linear_system),
      nonlinearManager_(nonlinear_manager),
      loader_(loader),
      topology_(topology),
      deviceManager_(device_manager),
      initialConditionsManager_(initial_conditions_manager)
  {}

  virtual ~ACFactory()
  {}

  //-----------------------------------------------------------------------------
  // Function      : create
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:59:00 2015
  //-----------------------------------------------------------------------------
  ///
  /// Create a new AC analysis and applies the analysis and time integrator option blocks.
  ///
  /// @return new AC analysis object
  ///
  AC *create() const
  {
    analysisManager_.setAnalysisMode(ANP_MODE_AC);
    AC *ac = new AC(analysisManager_, linearSystem_, nonlinearManager_, deviceManager_, loader_, topology_, initialConditionsManager_);
    ac->setAnalysisParams(acAnalysisOptionBlock_);
    ac->setTimeIntegratorOptions(timeIntegratorOptionBlock_);
    ac->setACLinSolOptions(acLinSolOptionBlock_);
    ac->setDCLinSolOptions(dcLinSolOptionBlock_);

    for (std::vector<Util::OptionBlock>::const_iterator it = ACLinOptionBlockVec_.begin(), end = ACLinOptionBlockVec_.end(); it != end; ++it)
    {
      ac->setACLinOptions(*it);
    }

    for (std::vector<Util::OptionBlock>::const_iterator it = dataOptionBlockVec_.begin(), end = dataOptionBlockVec_.end(); it != end; ++it)
    {
      ac->setDataStatements(*it);
    }

    ac->setSensAnalysisParams(sensAnalysisOptionBlock_);
    ac->setSensitivityOptions(sensitivityOptionBlock_);

    return ac;
  }

  //-----------------------------------------------------------------------------
  // Function      : setACAnalysisOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 13:00:14 2015
  //-----------------------------------------------------------------------------
  ///
  /// Saves the analysis parsed options block in the factory.
  ///
  /// @invariant Overwrites any previously specified analysis option block.
  ///
  /// @param option_block parsed option block
  ///
  void setACAnalysisOptionBlock(const Util::OptionBlock &option_block)
  {
    acAnalysisOptionBlock_ = option_block;
  }

  //-----------------------------------------------------------------------------
  // Function      : setTimeIntegratorOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 13:01:27 2015
  //-----------------------------------------------------------------------------
  ///
  /// Saves the time integrator parsed option block.
  ///
  /// @invariant Overwrites any previously specified time integrator option block.
  ///
  /// @param option_block parsed option block
  ///
  bool setTimeIntegratorOptionBlock(const Util::OptionBlock &option_block)
  {
    timeIntegratorOptionBlock_ = option_block;

    return true;
  }

  bool setACLinSolOptionBlock(const Util::OptionBlock &option_block)
  {
    acLinSolOptionBlock_ = option_block;
    return true;
  }

  bool setDCLinSolOptionBlock(const Util::OptionBlock &option_block)
  {
    dcLinSolOptionBlock_ = option_block;
    return true;
  }

  bool setDotDataBlock(const Util::OptionBlock &option_block)
  {
    dataOptionBlockVec_.push_back(option_block);
    return true;
  }

//-----------------------------------------------------------------------------
// Function      : setACLinOptionBlock
// Purpose       : Saves the parsed options blocks, from .LIN lines, that are
//                 relevant to the AC object in the factory.  There may be
//                 multiple .LIN lines in the netlist, since those lines
//                 also function as print lines for Touchstone output.
// Special Notes :
// Scope         : public
// Creator       : Ting Mei  
// Creation Date : 2/27/2019
//-----------------------------------------------------------------------------
  bool setACLinOptionBlock(const Util::OptionBlock &option_block)
  {
    ACLinOptionBlockVec_.push_back(option_block);

    return true;
  }

  bool setSensitivityOptionBlock(const Util::OptionBlock &option_block)
  {
    sensitivityOptionBlock_ = option_block;
    return true;
  }

  bool setSensAnalysisOptionBlock(const Util::OptionBlock &option_block)
  {
    sensAnalysisOptionBlock_ = option_block;
    return true;
  }

public:
  AnalysisManager &                     analysisManager_;
  Linear::System &                      linearSystem_;
  Nonlinear::Manager &                  nonlinearManager_;
  Loader::Loader &                      loader_;
  Topo::Topology &                      topology_;
  Device::DeviceMgr &                   deviceManager_; 
  IO::InitialConditionsManager &        initialConditionsManager_;

private:
  Util::OptionBlock     acAnalysisOptionBlock_;
  Util::OptionBlock     timeIntegratorOptionBlock_;
  Util::OptionBlock     acLinSolOptionBlock_;
  Util::OptionBlock     dcLinSolOptionBlock_;
  std::vector<Util::OptionBlock>        dataOptionBlockVec_;

  std::vector<Util::OptionBlock>        ACLinOptionBlockVec_;
  Util::OptionBlock     sensitivityOptionBlock_;
  Util::OptionBlock     sensAnalysisOptionBlock_;

};

// .AC
struct ACAnalysisReg : public IO::PkgOptionsReg
{
  ACAnalysisReg(
    ACFactory &   factory )
    : factory_(factory)
  {}

  bool operator()(const Util::OptionBlock &option_block)
  {
    factory_.setACAnalysisOptionBlock(option_block);

    factory_.analysisManager_.addAnalysis(&factory_);

    return true;
  }

  ACFactory &         factory_;
};

namespace {
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
struct isTokenString 
{
  isTokenString (std::string & test) : testString(test) {};
  bool operator() (const IO::StringToken & t1)
  {
    return compare_nocase(t1.string_.c_str(), testString.c_str()) == 0;
  }
  std::string & testString;
};
}

//-----------------------------------------------------------------------------
// Function      : extractACData
// Purpose       : Extract the parameters from a netlist .AC line held in
//                 parsed_line.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool extractACData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("AC", Util::OptionBlock::ALLOW_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);
  bool ret = extractACDataInternals(option_block, options_manager, netlist_filename, parsed_line);

  if (ret)
    circuit_block.addOptions(option_block);

  return ret;
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : extractACDataInternals
// Purpose       : Extract the parameters from a netlist .AC line
//-----------------------------------------------------------------------------
bool
extractACDataInternals(
  Util::OptionBlock &           option_block,
  IO::PkgOptionsMgr &           options_manager,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  int numFields = parsed_line.size();
  int linePosition = 1;   // Start of parameters on .param line.
  Util::Param parameter("", "");

  // check for "DATA" first.  If DATA is found, then skip everything else ----
  std::string tmp = std::string("DATA");
  IO::TokenVector::const_iterator startPL = parsed_line.begin();  startPL++;
  IO::TokenVector::const_iterator endPL = parsed_line.end();
  IO::TokenVector::const_iterator iter = std::find_if(startPL, endPL, isTokenString(tmp)); 
  if (iter != parsed_line.end())
  {
    int dataPos = std::distance(parsed_line.begin(),iter);
    if (numFields != 4)
    {
      Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
        << ".AC line not formatted correctly.  numFields = " << numFields;
      return false;
    }

    option_block.addParam( Util::Param("TYPE", "DATA"));
    option_block.addParam( Util::Param( "DATASET", parsed_line[ dataPos+2 ].string_ ));

    return true;
  }
  // End of the DATA block

  // Check that the minimum required number of fields are on the line.
  if (numFields != 5)
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << ".AC line has an unexpected number of fields";
    return false;
  }

  // type is required
  parameter.setTag( "TYPE" );
  ExtendedString stringVal ( parsed_line[linePosition].string_ );
  stringVal.toUpper();
  parameter.setVal(std::string(stringVal));
  option_block.addParam( parameter );

  ++linePosition;     // Advance to next parameter.

  // np is required
  parameter.setTag( "NP" );
  parameter.setVal( parsed_line[linePosition].string_ );
  option_block.addParam( parameter );
  ++linePosition;     // Advance to next parameter.

  // fstart is required
  parameter.setTag( "FSTART" );
  parameter.setVal( parsed_line[linePosition].string_ );
  option_block.addParam( parameter );
  ++linePosition;     // Advance to next parameter.


  // fstop is required
  parameter.setTag( "FSTOP" );
  parameter.setVal( parsed_line[linePosition].string_ );
  option_block.addParam( parameter );

  return true;
}


void
populateMetadata(
  IO::PkgOptionsMgr &   options_manager)
{
  {


    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("ACLIN");

//    parameters.insert(Util::ParamMap::value_type("AZ_max_iter", Util::Param("AZ_max_iter", 200)));

    parameters.insert(Util::ParamMap::value_type("sparcalc", Util::Param("sparcalc", 1)));
    parameters.insert(Util::ParamMap::value_type("lintype", Util::Param("lintype", "S")));
  }

  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("LINSOL-AC");

    parameters.insert(Util::ParamMap::value_type("AZ_max_iter", Util::Param("AZ_max_iter", 200)));
    parameters.insert(Util::ParamMap::value_type("AZ_precond", Util::Param("AZ_precond", 14)));
    parameters.insert(Util::ParamMap::value_type("AZ_solver", Util::Param("AZ_solver", 1)));
    parameters.insert(Util::ParamMap::value_type("AZ_conv", Util::Param("AZ_conv", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_pre_calc", Util::Param("AZ_pre_calc", 1)));
    parameters.insert(Util::ParamMap::value_type("AZ_keep_info", Util::Param("AZ_keep_info", 1)));
    parameters.insert(Util::ParamMap::value_type("AZ_orthog", Util::Param("AZ_orthog", 1)));
    parameters.insert(Util::ParamMap::value_type("AZ_subdomain_solve", Util::Param("AZ_subdomain_solve", 9)));
    parameters.insert(Util::ParamMap::value_type("AZ_ilut_fill", Util::Param("AZ_ilut_fill", 3.0)));
    parameters.insert(Util::ParamMap::value_type("AZ_drop", Util::Param("AZ_drop", 1.0E-3)));
    parameters.insert(Util::ParamMap::value_type("AZ_reorder", Util::Param("AZ_reorder", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_scaling", Util::Param("AZ_scaling", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_kspace", Util::Param("AZ_kspace", 50)));
    parameters.insert(Util::ParamMap::value_type("AZ_tol", Util::Param("AZ_tol", 1.0E-9)));
    parameters.insert(Util::ParamMap::value_type("AZ_output", Util::Param("AZ_output", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_diagnostics", Util::Param("AZ_diagnostics", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_overlap", Util::Param("AZ_overlap", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_rthresh", Util::Param("AZ_rthresh", 1.0001)));
    parameters.insert(Util::ParamMap::value_type("AZ_athresh", Util::Param("AZ_athresh", 1.0E-4)));
    parameters.insert(Util::ParamMap::value_type("AZ_filter", Util::Param("AZ_filter", 0.0)));
    parameters.insert(Util::ParamMap::value_type("TR_partition", Util::Param("TR_partition", 1)));
    parameters.insert(Util::ParamMap::value_type("TR_partition_type", Util::Param("TR_partition_type", "HYPERGRAPH")));
#ifdef Xyce_SHYLU
    parameters.insert(Util::ParamMap::value_type("ShyLU_rthresh", Util::Param("ShyLU_rthresh", 1.0E-3)));
#endif
    parameters.insert(Util::ParamMap::value_type("TR_reindex", Util::Param("TR_reindex", 1)));
    parameters.insert(Util::ParamMap::value_type("TR_solvermap", Util::Param("TR_solvermap", 1)));
    parameters.insert(Util::ParamMap::value_type("TR_amd", Util::Param("TR_amd", 1)));
    parameters.insert(Util::ParamMap::value_type("TR_btf", Util::Param("TR_btf", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_global_btf", Util::Param("TR_global_btf", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_global_btf_droptol", Util::Param("TR_global_btf_droptol", 1.0E-16)));
    parameters.insert(Util::ParamMap::value_type("TR_global_btf_verbose", Util::Param("TR_global_btf_verbose", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_global_amd", Util::Param("TR_global_amd", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_global_amd_verbose", Util::Param("TR_global_amd_verbose", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_singleton_filter", Util::Param("TR_singleton_filter", 0)));
    parameters.insert(Util::ParamMap::value_type("adaptive_solve", Util::Param("adaptive_solve", 0)));
    parameters.insert(Util::ParamMap::value_type("use_aztec_precond", Util::Param("use_aztec_precond", 1)));
    parameters.insert(Util::ParamMap::value_type("use_ifpack_factory", Util::Param("use_ifpack_factory", 0)));
    parameters.insert(Util::ParamMap::value_type("ifpack_type", Util::Param("ifpack_type", "Amesos")));
    parameters.insert(Util::ParamMap::value_type("diag_perturb", Util::Param("diag_perturb", 0.0)));
    parameters.insert(Util::ParamMap::value_type("TR_rcm", Util::Param("TR_rcm", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_scale", Util::Param("TR_scale", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_scale_left", Util::Param("TR_scale_left", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_scale_right", Util::Param("TR_scale_right", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_scale_exp", Util::Param("TR_scale_exp", 1.0)));
    parameters.insert(Util::ParamMap::value_type("TR_scale_iter", Util::Param("TR_scale_iter", 0)));
    parameters.insert(Util::ParamMap::value_type("TYPE", Util::Param("TYPE", "DEFAULT")));
    parameters.insert(Util::ParamMap::value_type("PREC_TYPE", Util::Param("PREC_TYPE", "DEFAULT")));
    parameters.insert(Util::ParamMap::value_type("IR_SOLVER_TYPE", Util::Param("IR_SOLVER_TYPE", "DEFAULT")));
    parameters.insert(Util::ParamMap::value_type("IR_SOLVER_TOL", Util::Param("IR_SOLVER_TOL", "DEFAULT")));
    parameters.insert(Util::ParamMap::value_type("BELOS_SOLVER_TYPE", Util::Param("BELOS_SOLVER_TYPE", "Block GMRES")));
    parameters.insert(Util::ParamMap::value_type("KLU_REPIVOT", Util::Param("KLU_REPIVOT", 1)));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_LS", Util::Param("OUTPUT_LS", 1)));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_BASE_LS", Util::Param("OUTPUT_BASE_LS", 1)));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_FAILED_LS", Util::Param("OUTPUT_FAILED_LS", 1)));
  }
}


bool
registerACFactory(
  FactoryBlock &        factory_block)
{
  ACFactory *factory = new ACFactory(factory_block.analysisManager_, factory_block.linearSystem_, factory_block.nonlinearManager_, factory_block.loader_, factory_block.topology_, factory_block.deviceManager_, factory_block.initialConditionsManager_);

  addAnalysisFactory(factory_block, factory);

  populateMetadata(factory_block.optionsManager_);

  factory_block.optionsManager_.addCommandParser(".AC", extractACData);

  factory_block.optionsManager_.addCommandProcessor("AC", new ACAnalysisReg(*factory));

  factory_block.optionsManager_.addOptionsProcessor("TIMEINT", IO::createRegistrationOptions(*factory, &ACFactory::setTimeIntegratorOptionBlock));

  factory_block.optionsManager_.addOptionsProcessor("LINSOL-AC", IO::createRegistrationOptions(*factory, &ACFactory::setACLinSolOptionBlock));

  factory_block.optionsManager_.addOptionsProcessor("LINSOL", IO::createRegistrationOptions(*factory, &ACFactory::setDCLinSolOptionBlock));

  factory_block.optionsManager_.addCommandProcessor("DATA",
    IO::createRegistrationOptions(*factory, &ACFactory::setDotDataBlock) );

  factory_block.optionsManager_.addOptionsProcessor("ACLIN", IO::createRegistrationOptions(*factory, &ACFactory::setACLinOptionBlock));

  factory_block.optionsManager_.addOptionsProcessor("SENS",
      IO::createRegistrationOptions(*factory, &ACFactory::setSensAnalysisOptionBlock));

  factory_block.optionsManager_.addOptionsProcessor("SENSITIVITY",
      IO::createRegistrationOptions(*factory, &ACFactory::setSensitivityOptionBlock));

  return true;
}

} // namespace Analysis
} // namespace Xyce

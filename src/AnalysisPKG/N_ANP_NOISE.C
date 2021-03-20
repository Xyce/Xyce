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


//-----------------------------------------------------------------------------
// Purpose       : NOISE analysis functions.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date :  12/8/2014
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iomanip>

#include <fstream>

#include <N_ANP_NOISE.h>
#include <N_ANP_NoiseData.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_SweepParam.h>
#include <N_ANP_SweepParamFreeFunctions.h>
#include <N_ERH_Message.h>
#include <N_ERH_ErrorMgr.h>
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
#include <N_TIA_WorkingIntegrationMethod.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_NoTimeIntegration.h>

#include <N_UTL_Diagnostic.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_Math.h>
#include <N_UTL_NodeSymbols.h>
#include <N_UTL_SaveIOSState.h>
#include <N_UTL_Timer.h>
#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Manager.h>
#include <N_PDS_Serial.h>

#include <N_PDS_ParMap.h>

#include <N_TOP_Topology.h>

#include<N_UTL_ExtendedString.h>
#include<N_NLS_ReturnCodes.h>

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;

// putting these here for now (from spice3's noisedefs.h file)
#define N_MINLOG          1E-38       /* the smallest number we can take the log of */
#define N_MINGAIN         1E-20  // the smallest input-output gain we can tolerate
                                 // (to calculate input-referred noise we divide
                                 // the output noise by the gain)

#define N_INTFTHRESH   1E-10       // the largest slope (of a log-log noise spectral
                                   //    density vs. freq plot) at which the noise
                                   //    spectum is still considered flat. (no need for
                                   //    log curve fitting)
#define N_INTUSELOG      1E-10       // decides which expression to use for the integral of
                                     //    x**k.  If k is -1, then we must use a 'ln' form.
                                     //    Otherwise, we use a 'power' form.  This
                                     //    parameter is the region around (k=) -1 for which
                                     //    use the 'ln' form.

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : NOISE::setTimeIntegratorOptions
// Purpose       :
// Special Notes : These are from '.options timeint'
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool NOISE::setTimeIntegratorOptions(
  const Util::OptionBlock &     option_block)
{
  for (Util::ParamList::const_iterator it = option_block.begin(),
      end = option_block.end(); it != end; ++it)
  {
    const Util::Param &param = *it;

    if (param.uTag() == "DEBUGLEVEL" )
    {
      IO::setTimeIntegratorDebugLevel(analysisManager_.getCommandLine(), param.getImmutableValue<int>());
    }
    else if (nonlinearManager_.setReturnCodeOption(param))
    {
      ;
    }
    else if (tiaParams_.setTimeIntegratorOption(param))
    {
      ;
    }
    else if (setDCOPOption(param))
    {
      ;
    }
    else
    {
      Report::UserError() << param.uTag()
        << " is not a recognized time integration option";
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
bool NOISE::setACLinSolOptions(const Util::OptionBlock &option_block)
{
  acLinSolOptionBlock_ = option_block;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : NOISE::setDataStatements
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/5/18
//-----------------------------------------------------------------------------
bool NOISE::setDataStatements(const Util::OptionBlock & paramsBlock)
{
  return processDataStatements(paramsBlock, dataNamesMap_, dataTablesMap_);
}

//-----------------------------------------------------------------------------
// Function      : NOISE::convertDataToSweepParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/5/18
//-----------------------------------------------------------------------------
bool NOISE::convertDataToSweepParams()
{
  return convertData(noiseSweepVector_, dataNamesMap_, dataTablesMap_);
}

//-----------------------------------------------------------------------------
// Function      : NOISE::NOISE
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
NOISE::NOISE(
  AnalysisManager &                     analysis_manager,
  Linear::System &                      linear_system,
  Nonlinear::Manager &                  nonlinear_manager,
  Loader::Loader &                      loader,
  Topo::Topology &                      topology,
  IO::InitialConditionsManager &        initial_conditions_manager)
  : AnalysisBase(analysis_manager, "NOISE"),
    StepEventListener(&analysis_manager),
    analysisManager_(analysis_manager),
    loader_(loader),
    linearSystem_(linear_system),
    nonlinearManager_(nonlinear_manager),
    topology_(topology),
    initialConditionsManager_(initial_conditions_manager),
    outputManagerAdapter_(analysis_manager.getOutputManagerAdapter()),
    bVecRealPtr(linearSystem_.builder().createVector()),
    bVecImagPtr(linearSystem_.builder().createVector()),
    bNoiseVecRealPtr(linearSystem_.builder().createVector()),
    bNoiseVecImagPtr(linearSystem_.builder().createVector()),
    hackOutputCalledBefore_(false),
    outputNodeSingle_(true),
    outputNode1_(""),
    outputNode2_(""),
    specifiedSource_(""),
    noiseLoopSize_(0),
    stepFlag_(false),
    type_("DEC"),
    np_(10.0),
    fStart_(1.0),
    fStop_(1.0),
    stepMult_(0.0),
    fstep_(0.0),
    pts_per_summary_(0),
    pts_per_summary_Given(false),
    calcNoiseIntegrals_(true),
    delFreq_(0.0),
    lastFreq_(0.0),
    currentFreq_(0.0),
    lnFreq_(0.0),
    lnLastFreq_(0.0),
    delLnFreq_(0.0),
    GainSqInv_(0.0),
    lnGainInv_(0.0),
    totalOutputNoise_(0.0),
    totalInputNoise_(0.0),
    totalOutputNoiseDens_(0.0),
    totalInputNoiseDens_(0.0),
    ACMatrix_(0),
    B_(0),
    X_(0),
    saved_AC_X_(0),
    blockSolver_(0),
    blockProblem_(0)
{
  bVecRealPtr->putScalar(0.0);
  bVecImagPtr->putScalar(0.0);
  bNoiseVecRealPtr->putScalar(0.0);
  bNoiseVecImagPtr->putScalar(0.0);

  outputManagerAdapter_.setDotNoiseSpecified(true);

  // noiseDataVec holds the total noise results and the noise
  // results for each device
  int numNoiseDevices = loader_.getNumNoiseDevices();
  noiseDataVec_.resize(numNoiseDevices);
  for (int i=0;i<numNoiseDevices;++i)
  {
    noiseDataVec_[i] = new NoiseData();
  }

  // Note: setting up the noise sources and putting the relevant entries
  // into the symbol table must be done after the devices are created,
  // but before the creation of the DNI() and DNO() operators.

  // set up the noise sources in each device
  loader_.setupNoiseSources(noiseDataVec_);

  // Put an entry for each device (with noise source(s)) into the
  // symbol table owned by Topology.  This will be used during
  // operator creation for DNI() and DNO().  DNO or DNI operators come
  // in two forms, DNO(deviceName) or DNO(deviceName,noiseSource)
  Util::SymbolTable& symbol_table = topology_.getNodeSymbols();
  std::multimap<std::string,int> noiseNamesMap;
  for (int i=0;i<numNoiseDevices;++i)
  {
    // noiseDataVec_[i]->deviceName is the individual device name (e.g., Q1)
    addSymbol(symbol_table, Util::NOISE_DEVICE_SYMBOL, i, noiseDataVec_[i]->deviceName +"_ND");

    // Account for duplicate names in noiseDataVec_[i]->noiseNames, which happens with
    // some of the ADMS device, by placing them into a multimap first before adding them
    // to the symbol table.
    for (int j=0;j<noiseDataVec_[i]->noiseNames.size();++j)
    {
      std::string prefix = "noise_" + noiseDataVec_[i]->deviceName;
      if (prefix == noiseDataVec_[i]->noiseNames[j])
      {
        // noiseDataVec_[i]->noiseNames[j] are the names of the noise types (e.g., rc, rb
        // re, ic, ib and fn for a Q device).  Don't add entries if noiseName[j] is equal
        // to the string "noise_ + deviceName" (e.g., for R devices).  Those entries are
        // superfluous since (for example) DNO(R1,R1) doesn't work by design.  Only DNO(R1)
        // works.  This block should be changed from a "no op" if that design decision
        // changes.
      }
      else if (prefix.length() < noiseDataVec_[i]->noiseNames[j].length())
      {
        // For more complex devices, noiseNames[j] will be (for example) noise_Q1_RC .
        // For some ADMS device, the noise type (RC in this example) may have inconvenient
        // characters like ( or ) in it.  The ExtendedString method removeBadChars()
        // will remove them before insertion into noiseNamesMap.
        ExtendedString noiseType(noiseDataVec_[i]->noiseNames[j].substr(prefix.length()+1));
	std::string noiseName = prefix + "_" + noiseType.removeBadChars();
        noiseNamesMap.insert(std::pair<std::string,int>(noiseName, j));
      }
    }

    std::string prevName="";
    int suffix = 0;
    for (std::multimap<std::string,int>::iterator it=noiseNamesMap.begin(); it!=noiseNamesMap.end(); ++it)
    {
      // For ADMS devices, that may have duplicate entries for a given noise type, the entries
      // are "suffixed" with _0, _1, _2, etc.  If there are no duplicate entries (e.g., for
      // the Q device) then just the _0 suffix is used.
      (*it).first != prevName ? suffix=0 : ++suffix;
      std::ostringstream s;
      s << suffix;
      addSymbol(symbol_table, Util::NOISE_TYPE_SYMBOL, (*it).second, (*it).first + "_" + s.str());
      prevName = (*it).first;
    }

    noiseNamesMap.clear();
  }
}

//-----------------------------------------------------------------------------
// Function      : NOISE::~NOISE()
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
NOISE::~NOISE()
{
  delete bVecRealPtr;
  delete bVecImagPtr;
  delete bNoiseVecRealPtr;
  delete bNoiseVecImagPtr;
  delete ACMatrix_;
  delete B_;
  delete X_;
  delete saved_AC_X_;
  delete blockSolver_;
  delete blockProblem_;

  int size = noiseDataVec_.size();
  for (int i=0;i<size;++i)
  {
    delete (noiseDataVec_[i]);
  }
  noiseDataVec_.clear();
}

//-----------------------------------------------------------------------------
// Function      : NOISE::notify
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void NOISE::notify(const StepEvent &event)
{
  if (event.state_ == StepEvent::STEP_STARTED)
  {
    AnalysisBase::resetForStepAnalysis();

    stepFlag_ = true;
    analysisManager_.getStepErrorControl().resetAll(tiaParams_);

    bVecRealPtr->putScalar(0.0);
    bVecImagPtr->putScalar(0.0);
    bNoiseVecRealPtr->putScalar(0.0);
    bNoiseVecImagPtr->putScalar(0.0);

    // reset accumulators for total input and output noise
    totalInputNoise_ = 0.0;
    totalOutputNoise_ = 0.0;
  }
}

//-----------------------------------------------------------------------------
// Function      : NOISE::setAnalysisParams
// Purpose       :
// Special Notes : These are from the .NOISE statement.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
bool NOISE::setAnalysisParams(const Util::OptionBlock & paramsBlock)
{
  bool retval=true;

  // Check for DATA first.  If DATA is present, then use the sweep functions,
  // rather than the NOISE specific built-in ones.  This also supports the case
  // of having multiple .NOISE lines in the netlist, wherein only the last .NOISE
  // line is used.
  if (isDataSpecified(paramsBlock))
  {
    dataSpecification_ = true;
    type_="TYPE";
    noiseSweepVector_.push_back(parseSweepParams(paramsBlock.begin(), paramsBlock.end()));
  }

  for (Util::ParamList::const_iterator it = paramsBlock.begin(),
      end = paramsBlock.end(); it != end; ++it)
  {
    if ((*it).uTag() == "V")
    {
      if ((*it).getImmutableValue<double>()==1.0)
      {
        outputNodeSingle_ = true;
        Util::ParamList::const_iterator itNode = it;
        itNode++;
        outputNode1_ = (*itNode).uTag();
      }
      else if ((*it).getImmutableValue<double>()==2.0)
      {
        outputNodeSingle_ = false;
        Util::ParamList::const_iterator itNode = it;
        itNode++;
        outputNode1_ = (*itNode).uTag();
        itNode++;
        outputNode2_ = (*itNode).uTag();
      }
    }
    else if ((*it).uTag() == "SOURCE")
    {
      specifiedSource_ = (*it).stringValue();
    }
    else if ((*it).uTag() == "TYPE" && !dataSpecification_)
    {
      type_ = (*it).stringValue();
    }
    else if ((*it).uTag() == "NP")
    {
      np_ = (*it).getImmutableValue<double>();
      ExtendedString npStr((*it).stringValue());
      if ( !npStr.isInt() )
      {
        Report::UserError0() << "Points Value parameter on .NOISE line must be an integer";
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
    else if ((*it).uTag() == "PTS_PER_SUMMARY")
    {
      pts_per_summary_ = (*it).getImmutableValue<int>();
    }
  }

  // exit from here if DATA=<name> is used on the .NOISE line
  if (dataSpecification_) return retval;

  // debug output, when DATA=<name> is not used
  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    dout() << section_divider << std::endl
           << "NOISE simulation parameters"
           << std::endl;

    if (outputNodeSingle_)
    {
      dout() << "Output Node: V(" << outputNode1_ << ")" <<std::endl;
    }
    else
    {
      dout() << "Output Node: V(" << outputNode1_ << ","<<outputNode2_<<")" <<std::endl;
    }

    dout() << "specified source = " << specifiedSource_ << std::endl
           << "number of points  = " << np_ << std::endl
           << "starting frequency = " << fStart_ << std::endl
           << "stop frequency = " << fStop_ << std::endl
           << "pts_per_summary = " << pts_per_summary_
             << std::endl;
  }

  // error checking of parameters, when DATA=<name> is not used
  if ( np_ < 1 )
  {
    Report::UserError0() << "Points Value parameter on .NOISE line must be >= 1";
    retval = false;
  }
  if ( (fStart_ <=0) || (fStop_ <= 0) )
  {
    Report::UserError0() << "Illegal values for start or end frequencies on .NOISE line. " <<
       "Both values must be > 0";
    retval = false;
  }
  if ( fStop_ < fStart_ )
  {
    Report::UserError0() << "End frequency must not be less than start frequency on .NOISE line";
    retval = false;
  }

  return retval;
}

//-----------------------------------------------------------------------------
// Function      : NOISE::getDCOPFlag()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
bool NOISE::getDCOPFlag() const
{
  return getIntegrationMethod() == TimeIntg::NoTimeIntegration::type;
}

//-----------------------------------------------------------------------------
// Function      : NOISE::run()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
bool NOISE::doRun()
{
  return doInit() && doLoopProcess() && doFinish();
}

//-----------------------------------------------------------------------------
// Function      : NOISE::doInit()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
bool NOISE::doInit()
{
  bool bsuccess = true;

  // check if the "DATA" specification was used.  If so, create a new vector of
  // SweepParams, in the "TABLE" style.
  if (dataSpecification_)
  {
    if (!convertDataToSweepParams())
    {
      Report::UserFatal() << "Invalid data=<name> parameter on .NOISE line.";
      return false;
    }

    std::vector<SweepParam>::iterator begin = noiseSweepVector_.begin();
    std::vector<SweepParam>::iterator end = noiseSweepVector_.end();
    std::vector<SweepParam>::iterator it = begin;
    for ( ; it != end; ++it)
    {
      SweepParam &sweep_param = (*it);
      std::string name = (*it).name; Util::toUpper(name);
      if (name == "FREQ" || name == "HERTZ")
      {
        // used to check that the specified frequencies are monotonically
        // increasing, to determine whether the noise integrals can be
        // calculated when DATA=<name> is used on the .NOISE line
        double prevFreq=-1.0;

        // frequency values for .NOISE must be > 0
        for (int i=0; i<(*it).valList.size(); ++i)
	{
          if ( (*it).valList[i] <= 0 )
	  {
            Report::UserFatal() << "Frequency values in .DATA for .NOISE analysis must be > 0";
            return false;
          }
          if ( calcNoiseIntegrals_ && ((*it).valList[i] <= prevFreq) )
	  {
            calcNoiseIntegrals_ = false;
	    Report::UserWarning0() << "Total Noise Integrals will not be calculated, "
		<< "since frequencies in .DATA table are not monotonically increasing";
          }
          prevFreq = (*it).valList[i];
        }
      }
      else
      {
        loader_.getParamAndReduce(analysisManager_.getComm(), sweep_param.name);
      }
    }

    // now set up the looping, etc
    noiseLoopSize_ = setSweepLoopVals(begin, end);
  }
  else
  {
    noiseLoopSize_ = setupSweepParam_();
  }

  // Get set to do the operating point.
  baseIntegrationMethod_ = TimeIntg::NoTimeIntegration::type;
  analysisManager_.createTimeIntegratorMethod(tiaParams_, baseIntegrationMethod_);

  stepNumber = 0;
  setDoubleDCOPEnabled(loader_.isPDESystem());

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish
    (AnalysisEvent(AnalysisEvent::INITIALIZE, AnalysisEvent::NOISE_IC));

  // set initial guess, if there is one to be set.
  // this setInitialGuess call is to up an initial guess in the
  // devices that have them (usually PDE devices).  This is different than
  // the "intializeProblem" call, which sets IC's.  (initial conditions are
  // different than initial guesses.
  loader_.setInitialGuess (analysisManager_.getDataStore()->nextSolutionPtr);

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

  // solving for DC op
  doHandlePredictor();
  loader_.updateSources();
  analysisManager_.getStepErrorControl().newtonConvergenceStatus = nonlinearManager_.solve();
  analysisManager_.getWorkingIntegrationMethod().stepLinearCombo ();
  gatherStepStatistics(stats_, nonlinearManager_.getNonlinearSolver(), analysisManager_.getStepErrorControl().newtonConvergenceStatus);
  analysisManager_.getStepErrorControl().evaluateStepError(loader_, tiaParams_);

  if ( analysisManager_.getStepErrorControl().newtonConvergenceStatus <= 0)
  {
    Report::UserError() << "Solving for DC operating point failed! Cannot continue NOISE analysis";
    return false;
  }

  // only output DC op if the .op was specified in the netlist
  // or if this NOISE analysis is called from a .step loop
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

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish
    (AnalysisEvent(AnalysisEvent::FINISH, AnalysisEvent::NOISE_IC));

  processOutputNodes ();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : NOISE::loopProcess()
// Purpose       : Conduct the stepping loop.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
bool NOISE::doLoopProcess()
{
  updateACLinearSystem_C_and_G_();
  createACLinearSystem_();

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish
    (AnalysisEvent(AnalysisEvent::INITIALIZE, AnalysisEvent::NOISE));

  int numNoiseDevices = loader_.getNumNoiseDevices();
  // clear out the integral arrays
  for (int i=0;i<numNoiseDevices;++i)
  {
    int numNoiseThisDevice = noiseDataVec_[i]->numSources;
    for (int j=0;j<numNoiseThisDevice;++j)
    {
      noiseDataVec_[i]->inputNoiseTotal[j] = 0.0;
      noiseDataVec_[i]->outputNoiseTotal[j] = 0.0;
    }
  }

  setupAdjointRHS_();

  Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();
  Parallel::Communicator &comm = *(pds_manager.getPDSComm());
  int myPID = comm.procID();

  ///////////////////////////////////////////////////////////////////////////
  // frequency loop
  // loop over all the specified frequencies, do an AC solve and a NOISE solve at each.
  for (int currentStep = 0; currentStep < noiseLoopSize_; ++currentStep)
  {
    // solve the AC system, to get up-to-date currents and voltages
    if (dataSpecification_)
    {
      updateDataParams_(currentStep);
    }
    else
    {
      updateCurrentFreq_(currentStep);
    }

    updateACLinearSystem_C_and_G_();
    updateACLinearSystemFreq_();
    updateACLinearSystemMagAndPhase_();

    static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish
      (AnalysisEvent(AnalysisEvent::STEP_STARTED, AnalysisEvent::NOISE, currentFreq_, currentStep));

    bool stepAttemptStatus;
    {
      Stats::StatTop _nonlinearStat("Nonlinear Solve");
      Stats::TimeBlock _nonlinearTimer(_nonlinearStat);

      stepAttemptStatus = solveACLinearSystem_();
    }

    // save a copy of X_ (the AC solution, already computed), for output purposes, etc.
    *saved_AC_X_ = *X_;

    // Compute AC gain.
    Linear::Vector & Xreal = X_->block( 0 );
    Linear::Vector & Ximag = X_->block( 1 );

    double v1r = 0.0;
    double v1i = 0.0;
    double v2r = 0.0;
    double v2i = 0.0;

    //comm.barrier();
    int root=-1;
    if (outputVarGIDs_.size()>0)
    {
      if (outputVarGIDs_[0] > -1)
      {
        v1r = Xreal.getElementByGlobalIndex(outputVarGIDs_[0]);
        v1i = Ximag.getElementByGlobalIndex(outputVarGIDs_[0]);
        root = myPID;
      }
      Xyce::Parallel::AllReduce(comm.comm(), MPI_MAX, &root, 1);
      comm.bcast( &v1r, 1, root );
      comm.bcast( &v1i, 1, root );
    }

    root=-1;
    if (outputVarGIDs_.size()>1)
    {
      if (outputVarGIDs_[1] > -1)
      {
        v2r = Xreal.getElementByGlobalIndex(outputVarGIDs_[1]);
        v2i = Ximag.getElementByGlobalIndex(outputVarGIDs_[1]);
        root = myPID;
      }
      Xyce::Parallel::AllReduce(comm.comm(), MPI_MAX, &root, 1);
      comm.bcast( &v2r, 1, root );
      comm.bcast( &v2i, 1, root );
    }

    double realVal = v1r-v2r;
    double imagVal = v1i-v2i;
    GainSqInv_ = 1.0 / std::max(((realVal*realVal) + (imagVal*imagVal)),N_MINGAIN);
    lnGainInv_ = std::log(GainSqInv_);

    // save previous (last) lnNoise densities
    for (int i=0;i<noiseDataVec_.size();++i)
    {
      int numNoiseThisDevice = noiseDataVec_[i]->numSources;
      for (int j=0;j<numNoiseThisDevice;++j)
      {
        noiseDataVec_[i]->lastLnNoiseDens[j] = noiseDataVec_[i]->lnNoiseDens[j];
      }
    }

    // do NOISE analysis for this frequency.
    resetAdjointNOISELinearSystem_();
    solveAdjointNOISE_();

    // Perform total noise integrals, if the specified frequency values are
    // monotonically increasing.  This is always true if DATA=<name> is NOT
    // used on the .NOISE line.
    if (currentStep != 0 && calcNoiseIntegrals_)
    {
      for (int i=0;i<noiseDataVec_.size();++i)
      {
        int numNoiseThisDevice = noiseDataVec_[i]->numSources;
        for (int j=0;j<numNoiseThisDevice;++j)
        {
          double noizDens = noiseDataVec_[i]->outputNoiseDens[j];
          double lnDens = noiseDataVec_[i]->lnNoiseDens[j];
          double lnlastDens = noiseDataVec_[i]->lastLnNoiseDens[j];

          double tempOutNoise = noiseIntegral( noizDens, lnDens, lnlastDens,
                 delLnFreq_, delFreq_, lnFreq_, lnLastFreq_);

          double tempInNoise = noiseIntegral(
                 noizDens * GainSqInv_,
                 lnDens + lnGainInv_,
                 lnlastDens + lnGainInv_,
                 delLnFreq_, delFreq_, lnFreq_, lnLastFreq_);


          noiseDataVec_[i]->outputNoiseTotal[j] += tempOutNoise;
          noiseDataVec_[i]->inputNoiseTotal[j] += tempInNoise;

          totalOutputNoise_+=tempOutNoise;
          totalInputNoise_+=tempInNoise;
        }
      }
    }

    // process success/failure
    if (stepAttemptStatus)
    {
      static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish
        (AnalysisEvent(AnalysisEvent::STEP_SUCCESSFUL, AnalysisEvent::NOISE, currentFreq_, currentStep));
      doProcessSuccessfulStep();
    }
    else // stepAttemptStatus  (ie do this if the step FAILED)
    {
      static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish
        (AnalysisEvent(AnalysisEvent::STEP_FAILED, AnalysisEvent::NOISE, currentFreq_, currentStep));
      doProcessFailedStep();
    }
  }

  Xyce::Parallel::AllReduce(comm.comm(), MPI_SUM, &totalOutputNoise_, 1);
  Xyce::Parallel::AllReduce(comm.comm(), MPI_SUM, &totalInputNoise_, 1);

  if (calcNoiseIntegrals_)
  {
    // Outputs to the screen
    noiseOutputToScreen_( Xyce::lout() );
  }

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish
    (AnalysisEvent(AnalysisEvent::FINISH, AnalysisEvent::NOISE));

  return true;
}

//-----------------------------------------------------------------------------
// Function      : NOISE::noiseOutputToScreen_()
// Purpose       : sends .NOISE output to the "screen" which is typically
//               : the stream Xyce::lout().  This allows that output to
//               : use scientific notation, without changing the existing
//               : state of the stream (os).
// Special Notes :
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 10/12/2017
//-----------------------------------------------------------------------------
std::ostream& NOISE::noiseOutputToScreen_(std::ostream& os)
{
  // save current stream state, and then set the stream to use scientific notation.
  // Otherwise the info for the stepped parameters may not be output in
  // scientific notation.
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os.setf(std::ios::scientific);

  if ( outputManagerAdapter_.getStepSweepVector().empty() )
  {
    // output for non-step case
    os << "Total Output Noise = " << totalOutputNoise_ <<std::endl;
    os << "Total Input Noise = " << totalInputNoise_ <<std::endl;
  }
  else
  {
    // Output for .STEP.  Step count should be output as an integer.
    // The other values should be output in scientific notation.
    os << "For Step " << outputManagerAdapter_.getStepAnalysisStepNumber() << ":" << std::endl;
    for (std::vector<Analysis::SweepParam>::const_iterator it = outputManagerAdapter_.getStepSweepVector().begin();
         it != outputManagerAdapter_.getStepSweepVector().end(); ++it)
    {
      os << it->name << " = " << it->currentVal << std::endl;
    }

    os << "Total Output Noise = " << totalOutputNoise_ <<std::endl;
    os << "Total Input Noise = " << totalInputNoise_ <<std::endl;

    if ( (outputManagerAdapter_.getStepAnalysisStepNumber()+1) <
          outputManagerAdapter_.getStepAnalysisMaxSteps() )
    {
      // add a blank line after each block of .STEP output, except for the last one
      os << std::endl;
    }
  }

  // previous state of os is restored when the function returns
  return os;
}

//-----------------------------------------------------------------------------
// Function      : NOISE::createACLinearSystem_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
bool NOISE::createACLinearSystem_()
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

  RCP<Linear::Graph> blockGraph = Linear::createBlockGraph( offset, blockPattern, *blockMap, *baseFullGraph);

  delete ACMatrix_;
  ACMatrix_ = Xyce::Linear::createBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), baseFullGraph);

  ACMatrix_->put( 0.0 ); // Zero out whole matrix.
  // Matrix will be loaded with nonzero (C,G) sub-matrices later.

  // Copy the values loaded into the blocks into the global matrix for the solve.
  ACMatrix_->assembleGlobalMatrix();

  B_->putScalar( 0.0 );
  B_->block( 0 ).update( 1.0, *bVecRealPtr);
  B_->block( 1 ).update( 1.0, *bVecImagPtr);

  delete X_;
  X_ = Xyce::Linear::createBlockVector (numBlocks, blockMap, baseMap);
  X_->putScalar( 0.0 );

  delete saved_AC_X_;
  saved_AC_X_ = Xyce::Linear::createBlockVector (numBlocks, blockMap, baseMap);
  saved_AC_X_->putScalar( 0.0 );

  delete blockProblem_;
  blockProblem_ = Xyce::Linear::createProblem( ACMatrix_, X_, B_ );

  delete blockSolver_;
  Linear::TranSolverFactory factory;
  blockSolver_ = factory.create( acLinSolOptionBlock_, *blockProblem_, analysisManager_.getCommandLine() );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : NOISE:updateLinearSystem_C_and_G_()
// Purpose       :
// Special Notes : This is only needed if device parameters are being swept.
//                 This only works for linear problems!
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/10/2018
//-----------------------------------------------------------------------------
bool NOISE::updateACLinearSystem_C_and_G_()
{
  analysisManager_.getDataStore()->daeQVectorPtr->putScalar(0.0);
  analysisManager_.getDataStore()->daeFVectorPtr->putScalar(0.0);

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
               (analysisManager_.getDataStore()->currStorePtr)
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

  if (DEBUG_TIME && isActive(Diag::TIME_PARAMETERS))
  {
    Xyce::dout() << "dQdxMatrixPtr:" << std::endl;
    analysisManager_.getDataStore()->dQdxMatrixPtr->print( Xyce::dout() );

    Xyce::dout() << "dFdxMatrixPtr:" << std::endl;
    analysisManager_.getDataStore()->dFdxMatrixPtr->print( Xyce::dout() );

    Xyce::dout() << std::endl;
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : NOISE::updateACLinearSystemFreq_()
// Purpose       :
// Special Notes : ERK.  This has been modified to also reload G, to accomodate
//                 model parameter sweeps on linear circuits.
//
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool NOISE::updateACLinearSystemFreq_()
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
// Function      : NOISE::updateACLinearSystemMagAndPhase_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/7/2018
//-----------------------------------------------------------------------------
bool NOISE::updateACLinearSystemMagAndPhase_()
{
  // ACMAG and ACPHASE have already been updated, so all that
  // remains is to load the B-vector and apply to the block system
  bVecRealPtr->putScalar(0.0);
  bVecImagPtr->putScalar(0.0);

  // re-load the B-vectors
  loader_.loadBVectorsforAC (bVecRealPtr, bVecImagPtr);

  B_->putScalar( 0.0 );
  B_->block( 0 ).update( 1.0, *bVecRealPtr);
  B_->block( 1 ).update( 1.0, *bVecImagPtr);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : NOISE::solveACLinearSystem_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
bool NOISE::solveACLinearSystem_()
{
  bool bsuccess = true;
  int linearStatus = blockSolver_->solve();
  if (linearStatus != 0)
  {
    Xyce::dout() << "Linear solve exited with error: " << linearStatus;
    bsuccess = false;
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : NOISE::processOutputNodes
// Purpose       : determines the GIDs for the nodes specified in the first argument
//                 of the .NOISE line.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/15/2014
//-----------------------------------------------------------------------------
void NOISE::processOutputNodes ()
{
  // setup the names:
  outputVarNames_.clear();
  outputVarNames_.push_back(outputNode1_);
  if (!outputNodeSingle_)
  {
    outputVarNames_.push_back(outputNode2_);
  }
  int numOutVars = outputVarNames_.size();

  // set up the gid's:
  int found(0);
  int found2(0);
  bool foundLocal(false);
  bool foundLocal2(false);

  Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();
  Parallel::Communicator &comm = *(pds_manager.getPDSComm());

  outputVarGIDs_.resize( numOutVars, -1 );
  for (int iout = 0; iout < numOutVars; ++iout)
  {
    std::vector<int> svGIDList1, dummyList;
    char type1;
    foundLocal = topology_.getNodeSVarGIDs(NodeID(outputVarNames_[iout], Xyce::_VNODE),
        svGIDList1, dummyList, type1);

    found = static_cast<int>(foundLocal);
    Xyce::Parallel::AllReduce(comm.comm(), MPI_LOR, &found, 1);

    foundLocal2 = false;
    if (!found)// if looking for this as a voltage node failed, try a "device" (i.e. current) node.
    {
      foundLocal2 = topology_.getNodeSVarGIDs(NodeID(outputVarNames_[iout], Xyce::_DNODE),
          svGIDList1, dummyList, type1);
    }
    found2 = static_cast<int>(foundLocal2);
    Xyce::Parallel::AllReduce(comm.comm(), MPI_LOR, &found2, 1);

    if (!found && !found2)
    {
      Report::UserError() << "Output function variable " << outputVarNames_[iout] << " not found";
    }

    if (found || found2)
    {
      int tmpGID=-1;
      if(svGIDList1.size()==1)
      {
        tmpGID = svGIDList1.front();
      }
      outputVarGIDs_[iout] = tmpGID;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : NOISE::setupAdjointRHS_
// Purpose       : Sets up stuff that only needs to be set up once,
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/15/2014
//-----------------------------------------------------------------------------
void NOISE::setupAdjointRHS_()
{
  // set up the RHS vector of 1's and 0's vector to represent dOutput/dX.
  // If thinking of this as the input to the adjoint circuit problem
  // (or adjoint operator) this can also be thought of as current sources to
  // the adjoint circuit.
  //
  // For the adjoint problem
  // assume that the phase is zero (as this is essentially the input to the
  // adjoint operator) and just set 1's on the real part of the block vector.
  // For this purpose re-use the B_ vector from the AC problem.

  // replace RHS in linear system.  Assume that the phase is zero, so no
  // imaginary contribution.

  Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();
  Parallel::Communicator &comm = *(pds_manager.getPDSComm());

  bNoiseVecRealPtr->putScalar(0.0);
  bNoiseVecImagPtr->putScalar(0.0);

  int numOutVars = outputVarNames_.size();
  for (int iout=0;iout<numOutVars;++iout)
  {
    int tmpGID=outputVarGIDs_[iout];
    if (tmpGID > -1)
    {
      double val=1.0;
      if (iout>0) val=-1.0;
      bNoiseVecRealPtr->setElementByGlobalIndex( tmpGID, val, 0);
    }
  }
  bNoiseVecRealPtr->fillComplete();
}

//-----------------------------------------------------------------------------
// Function      : NOISE::resetAdjointNOISELinearSystem_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/15/2014
//-----------------------------------------------------------------------------
void NOISE::resetAdjointNOISELinearSystem_()
{
  // clear out the X_ vector used by the linear solver.
  X_->putScalar(0.0);

  // setup the B_ vector RHS for the adjoint solve:
  B_->putScalar( 0.0 );
  B_->block( 0 ).update( 1.0, *bNoiseVecRealPtr);
  B_->block( 1 ).update( 1.0, *bNoiseVecImagPtr);

  if (DEBUG_ANALYSIS)
  {
    Xyce::dout()<<"adjoint noise B vector:"<<std::endl;
    B_->print(Xyce::dout());
  }

  if (DEBUG_ANALYSIS)
  {
    Xyce::dout()<<"matrix:"<<std::endl;
    ACMatrix_->print(Xyce::dout());
  }
}

//-----------------------------------------------------------------------------
// Function      : NOISE::solveAdjointNOISE_
// Purpose       : Solves for NOISE using the adjoint method.
// Special Notes : This is the production solve.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/15/2014
//-----------------------------------------------------------------------------
bool NOISE::solveAdjointNOISE_()
{
  bool bsuccess = true;
  int linearStatus= blockSolver_->solveTranspose();
  if (linearStatus != 0)
  {
    Xyce::dout() << "Linear solve exited with error: " << linearStatus;
    bsuccess = false;
  }

  double omega =  2.0 * M_PI * currentFreq_;

  int numNoiseDevices = noiseDataVec_.size();
  for (int i=0;i<numNoiseDevices;++i)
  {
    (noiseDataVec_[i])->omega = omega;
    (noiseDataVec_[i])->freq = currentFreq_;
  }

  loader_.getNoiseSources(noiseDataVec_);

  std::vector< Teuchos::RCP<Linear::Vector> > outputVectors;
  outputVectors.push_back( rcp(linearSystem_.builder().createVector()) );
  outputVectors.push_back( rcp(linearSystem_.builder().createVector()) );

  copyFromBlockVector( *X_, outputVectors );
  Linear::Vector & Xreal = *(outputVectors[0]);
  Linear::Vector & Ximag = *(outputVectors[1]);

  if (DEBUG_ANALYSIS)
  {
    std::cout << "Xreal: ------------------------------------"<<std::endl;
    Xreal.print( Xyce::dout() );
    std::cout << "Ximag: ------------------------------------"<<std::endl;
    Ximag.print( Xyce::dout() );
  }

  Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();
  Parallel::Communicator &comm = *(pds_manager.getPDSComm());

  totalOutputNoiseDens_ = 0.0;
  totalInputNoiseDens_ = 0.0;
  for (int i=0;i<numNoiseDevices;++i)
  {
    int numNoiseThisDevice = noiseDataVec_[i]->numSources;

    noiseDataVec_[i]->totalNoise = 0.0;
    noiseDataVec_[i]->totalOutputNoise = 0.0;
    for (int j=0;j<numNoiseThisDevice;++j)
    {
      int li_Pos = noiseDataVec_[i]->li_Pos[j];
      int li_Neg = noiseDataVec_[i]->li_Neg[j];

      double realVal = ((li_Pos!=-1)?Xreal[li_Pos]:0) - ((li_Neg!=-1)?Xreal[li_Neg]:0);
      double imagVal = ((li_Pos!=-1)?Ximag[li_Pos]:0) - ((li_Neg!=-1)?Ximag[li_Neg]:0);
      double gain = (realVal*realVal) + (imagVal*imagVal);

      noiseDataVec_[i]->totalNoise += noiseDataVec_[i]->noiseDens[j];
      noiseDataVec_[i]->outputNoiseDens[j] = gain * noiseDataVec_[i]->noiseDens[j];
      noiseDataVec_[i]->lnNoiseDens[j] = std::log(std::max( noiseDataVec_[i]->outputNoiseDens[j],N_MINLOG) );
      noiseDataVec_[i]->inputNoiseDens[j] = noiseDataVec_[i]->outputNoiseDens[j] * GainSqInv_;
      noiseDataVec_[i]->totalOutputNoise += noiseDataVec_[i]->outputNoiseDens[j];
    }
    noiseDataVec_[i]->totalInputNoise = noiseDataVec_[i]->totalOutputNoise * GainSqInv_;
    totalOutputNoiseDens_ += noiseDataVec_[i]->totalOutputNoise;
  }
  Xyce::Parallel::AllReduce(comm.comm(), MPI_SUM, &totalOutputNoiseDens_, 1);
  totalInputNoiseDens_ += totalOutputNoiseDens_ * GainSqInv_;

  if (comm.isSerial() )
  {
    // FIX:  replace this output call!
    hackTecplotOutput();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : NOISE::hackTecplotOutput
// Purpose       : provides an output function until I setup the output
//                 manager properly for NOISE.
// Author        : Eric Keiter
// Date          : 1/5/2015
//-----------------------------------------------------------------------------
void NOISE::hackTecplotOutput()
{
  int numNoiseDevices = noiseDataVec_.size();
  std::string fileName = analysisManager_.getNetlistFilename() + "_noise.dat";

  std::ofstream output_stream(fileName.c_str(),
            !hackOutputCalledBefore_ ? std::ios_base::out : std::ios_base::app);

  if (!hackOutputCalledBefore_)
  {
    // header output
    output_stream << "TITLE=\"noise output\"\tVARIABLES=\"frequency\" "<<std::endl;
    output_stream << "\t\"Re(V(" << outputNode1_ << "))\""<<std::endl;
    output_stream << "\t\"Im(V(" << outputNode1_ << "))\""<<std::endl;

    if(!outputNodeSingle_)
    {
      output_stream << "\t\"Re(V(" << outputNode2_ << "))\""<<std::endl;
      output_stream << "\t\"Im(V(" << outputNode2_ << "))\""<<std::endl;
    }
    output_stream << "\t\"onoise_spectrum \""<<std::endl;
    output_stream << "\t\"inoise_spectrum \""<<std::endl;

    if(pts_per_summary_>0)
    {
      // output noise variable names
      for (int i=0;i<numNoiseDevices;++i)
      {
        int numNoiseThisDevice = noiseDataVec_[i]->numSources;

        for (int j=0;j<numNoiseThisDevice;++j)
        {
          std::string namestring = "o" + noiseDataVec_[i]->noiseNames[j];
          std::transform(namestring.begin(), namestring.end(), namestring.begin(), ::tolower);
          output_stream << "\t\"" << namestring << "\"";
        }

        if (numNoiseThisDevice > 1)
        {
          std::string totalOutputNoiseName = "onoise_" + noiseDataVec_[i]->deviceName;
          std::transform(totalOutputNoiseName.begin(), totalOutputNoiseName.end(), totalOutputNoiseName.begin(), ::tolower);
          output_stream << "\t\""<<totalOutputNoiseName<<"\"";
        }
      }

      // input noise variable names
      for (int i=0;i<numNoiseDevices;++i)
      {
        int numNoiseThisDevice = noiseDataVec_[i]->numSources;

        for (int j=0;j<numNoiseThisDevice;++j)
        {
          std::string namestring = "i" + noiseDataVec_[i]->noiseNames[j];
          std::transform(namestring.begin(), namestring.end(), namestring.begin(), ::tolower);
          output_stream << "\t\"" << namestring << "\"";
        }

        // spice doesn't output input noise for individual devices, so exclude.
        if (numNoiseThisDevice > 1)
        {
          std::string totalOutputNoiseName = "inoise_" + noiseDataVec_[i]->deviceName;
          std::transform(totalOutputNoiseName.begin(), totalOutputNoiseName.end(), totalOutputNoiseName.begin(), ::tolower);
          output_stream << "\t\""<<totalOutputNoiseName<<"\"";
        }
      }
    }

    output_stream << std::endl;
    output_stream << "ZONE F=POINT  T=\"Xyce data\""<<std::endl;
    hackOutputCalledBefore_ = true;
  }

  // data output
  output_stream.setf(std::ios::scientific);

  // output the AC solution for the output node(s) (for double checking)
  Linear::Vector & Xreal = saved_AC_X_->block( 0 );
  Linear::Vector & Ximag = saved_AC_X_->block( 1 );
  double v1r = Xreal.getElementByGlobalIndex(outputVarGIDs_[0]);
  double v1i = Ximag.getElementByGlobalIndex(outputVarGIDs_[0]);

  output_stream << currentFreq_ << "\t" << v1r << "\t" << v1i;

  if(!outputNodeSingle_)
  {
    double v2r = 0.0;
    double v2i = 0.0;
    v2r = Xreal.getElementByGlobalIndex(outputVarGIDs_[1]);
    v2i = Ximag.getElementByGlobalIndex(outputVarGIDs_[1]);
    output_stream << "\t" << v2r << "\t" << v2i;
  }

  output_stream << "\t" << totalOutputNoiseDens_ << "\t" << totalInputNoiseDens_;

  if(pts_per_summary_>0)
  {
    // output noise values
    for (int i=0;i<numNoiseDevices;++i)
    {
      int numNoiseThisDevice = noiseDataVec_[i]->numSources;

      for (int j=0;j<numNoiseThisDevice;++j)
      {
        output_stream << "\t" << noiseDataVec_[i]->outputNoiseDens[j];
      }

      if (numNoiseThisDevice > 1)
      {
        output_stream << "\t" << noiseDataVec_[i]->totalOutputNoise;
      }
      //output_stream << std::endl;
    }

    // input noise values
    for (int i=0;i<numNoiseDevices;++i)
    {
      int numNoiseThisDevice = noiseDataVec_[i]->numSources;

      for (int j=0;j<numNoiseThisDevice;++j)
      {
        output_stream << "\t" << noiseDataVec_[i]->inputNoiseDens[j];
      }

      if (numNoiseThisDevice > 1)
      {
        output_stream << "\t" << noiseDataVec_[i]->totalInputNoise;
      }
    }
  }
  output_stream << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : NOISE::noiseIntegral
// Purpose       : based on the spice3 function: Nintegrate.
//  This subroutine evaluates the integral of the function
//
//                                             EXPONENT
//                      NOISE = a * (FREQUENCY)
//
//   given two points from the curve.  If EXPONENT is relatively close
//   to 0, the noise is simply multiplied by the change in frequency.
//   If it isn't, a more complicated expression must be used.  Note that
//   EXPONENT = -1 gives a different equation than EXPONENT <> -1.
//   Hence, the reason for the constant 'N_INTUSELOG'.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
double NOISE::noiseIntegral(double noizDens, double lnNdens, double lnNlstDens,
      double delLnFreq, double delFreq, double lnFreq, double lnLastFreq)
{
  double exponent;
  double a;

  exponent = (lnNdens - lnNlstDens) / delLnFreq;
  if ( fabs(exponent) < N_INTFTHRESH )
  {
    return (noizDens * delFreq);
  } else
  {
    a = std::exp(lnNdens - exponent*lnFreq);
    exponent += 1.0;
    if (fabs(exponent) < N_INTUSELOG)
    {
      return (a * (lnFreq - lnLastFreq));
    }
    else
    {
      return (a * ((std::exp(exponent * lnFreq) - std::exp(exponent * lnLastFreq)) /
            exponent));
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : NOISE::doProcessSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
bool NOISE::doProcessSuccessfulStep()
{
  bool bsuccess = true;

  // The raw output (via -r and -a) still uses .outputAC.  So, this possibly still needs
  // to be fixed or changed.  The outputACwoMeasureUpdates() function does not try to
  // update any AC mode measures.
  outputManagerAdapter_.outputACwoMeasureUpdates (currentFreq_, fStart_,fStop_,
	    saved_AC_X_->block(0), saved_AC_X_-> block(1), RFparams_);

  outputManagerAdapter_.outputNoise (currentFreq_, fStart_, fStop_, saved_AC_X_->block(0), saved_AC_X_-> block(1),
     totalOutputNoiseDens_, totalInputNoiseDens_, noiseDataVec_);

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

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : NOISE::processFailedStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
bool NOISE::doProcessFailedStep()
{
  bool bsuccess = true;

  stepNumber += 1;
  noiseSweepFailures_.push_back(stepNumber);
  stats_.failedStepsAttempted_  += 1;
  analysisManager_.getStepErrorControl().numberSuccessiveFailures += 1;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : NOISE::finish
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
bool NOISE::doFinish()
{
  bool bsuccess = true;

  if (DEBUG_ANALYSIS)
  {
    Xyce::dout() << "Calling NOISE::doFinish() outputs!" << std::endl;
  }

  // FIX OUTPUT
  outputManagerAdapter_.finishOutput();

  if (!(noiseSweepFailures_.empty()))
  {
    bsuccess = false;
  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : NOISE::doHandlePredictor
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
bool NOISE::doHandlePredictor()
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
// Function      : NOISE::updateDataParams_
// Purpose       :
// Special Notes : Used for AC analysis classes, when .DATA is used
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/10/2018
//-----------------------------------------------------------------------------
bool NOISE::updateDataParams_ (int stepNumber)
{
  bool reset = updateSweepParams(stepNumber, noiseSweepVector_.begin(), noiseSweepVector_.end());

  bool nonFreqPresent=false;
  lastFreq_ = currentFreq_;
  for (int iac=0;iac<noiseSweepVector_.size();++iac)
  {
    std::string name = noiseSweepVector_[iac].name; Util::toUpper(name);
    double val = noiseSweepVector_[iac].currentVal;
    if (name == "FREQ" || name == "HERTZ")
    {
      currentFreq_ = val;

      delFreq_ = currentFreq_ - lastFreq_;
      lnFreq_     = std::log(std::max(currentFreq_,N_MINLOG));
      lnLastFreq_ = std::log(std::max(lastFreq_,N_MINLOG));
      delLnFreq_  = lnFreq_ - lnLastFreq_;
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
// Function      : NOISE::updateCurrentFreq_(int stepNumber )
// Purpose       :
// Special Notes : Used for NOISE analysis classes.
// Scope         : public
// Creator       : Eric Keiter, SNL.
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
bool NOISE::updateCurrentFreq_(int stepNumber)
{
  lastFreq_ = currentFreq_;
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
    Report::DevelFatal().in("NOISE::updateCurrentFreq_")
      << "NOISE::updateCurrentFreq_: unsupported STEP type";
  }

  delFreq_ = currentFreq_ - lastFreq_;

  lnFreq_     = std::log(std::max(currentFreq_,N_MINLOG));
  lnLastFreq_ = std::log(std::max(lastFreq_,N_MINLOG));
  delLnFreq_  = lnFreq_ - lnLastFreq_;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : NOISE::setupSweepParam_
// Purpose       : Processes sweep parameters.
// Special Notes : Used for NOISE analysis classes.
// Scope         : public
// Creator       : Eric Keiter, SNL.
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
int NOISE::setupSweepParam_()
{
  double fstart, fstop;
  double fcount = 0.0;

  fstart = fStart_;
  fstop = fStop_;

  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    Xyce::dout() << std::endl << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "NOISE::setupSweepParam_" << std::endl;
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
      stepMult_ = std::pow(10.0, 1.0/np_);
      fcount   = floor(fabs(std::log10(fstart) - std::log10(fstop)) * np_ + 1.0);
      if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
      {
        Xyce::dout() << "stepMult_ = " << stepMult_  << std::endl;
      }
    }
    else if (type_ == "OCT")
    {
      stepMult_ = std::pow(2.0, 1.0/np_);

      // changed to remove dependence on "log2" function, which apparently
      // doesn't exist in the math libraries of FreeBSD or the mingw
      // cross-compilation suite.   Log_2(x)=log_e(x)/log_e(2.0)
      double ln2 = std::log(2.0);
      fcount   = floor(fabs(std::log(fstart) - std::log(fstop))/ln2 * np_ + 1.0);
      if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
      {
        Xyce::dout() << "stepMult_ = " << stepMult_  << std::endl;
      }
    }
    else
    {
      Report::UserFatal0() << "Unsupported NOISE sweep type: " << type_;
    }

  // At this point, pinterval equals the total number of steps
  // for the step loop.
  return static_cast<int> (fcount);
}

namespace {

typedef Util::Factory<AnalysisBase, NOISE>  NOISEFactoryBase;

//-----------------------------------------------------------------------------
// Class         : NOISEFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jan 29 12:53:02 2015
//-----------------------------------------------------------------------------
///
/// Factory for parsing NOISE parameters from the netlist and creating NOISE analysis.
///
class NOISEFactory : public NOISEFactoryBase
{
public:
  //-----------------------------------------------------------------------------
  // Function      : NOISEFactory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:54:09 2015
  //-----------------------------------------------------------------------------
  ///
  /// Constructs the NOISE analysis factory
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
  NOISEFactory(
    Analysis::AnalysisManager & analysis_manager,
    Linear::System &            linear_system,
    Nonlinear::Manager &        nonlinear_manager,
    Loader::Loader &            loader,
    Topo::Topology &            topology,
    IO::InitialConditionsManager &      initial_conditions_manager)
    : NOISEFactoryBase(),
      analysisManager_(analysis_manager),
      linearSystem_(linear_system),
      nonlinearManager_(nonlinear_manager),
      loader_(loader),
      topology_(topology),
      initialConditionsManager_(initial_conditions_manager)
  {}

  virtual ~NOISEFactory()
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
  /// Create a new NOISE analysis and applies the analysis and time integrator option blocks.
  ///
  /// @return new NOISE analysis object
  ///
  NOISE *create() const
  {
    analysisManager_.setAnalysisMode(ANP_MODE_NOISE);
    NOISE *noise = new NOISE(analysisManager_, linearSystem_, nonlinearManager_, loader_, topology_, initialConditionsManager_);
    noise->setAnalysisParams(noiseAnalysisOptionBlock_);
    noise->setTimeIntegratorOptions(timeIntegratorOptionBlock_);
    noise->setACLinSolOptions(acLinSolOptionBlock_);

    for (std::vector<Util::OptionBlock>::const_iterator it = dataOptionBlockVec_.begin(), end = dataOptionBlockVec_.end(); it != end; ++it)
    {
      noise->setDataStatements(*it);
    }

    return noise;
  }

  //-----------------------------------------------------------------------------
  // Function      : setNOISEAnalysisOptionBlock
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
  void setNOISEAnalysisOptionBlock(const Util::OptionBlock &option_block)
  {
    noiseAnalysisOptionBlock_ = option_block;
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

  bool setDotDataBlock(const Util::OptionBlock &option_block)
  {
    dataOptionBlockVec_.push_back(option_block);
    return true;
  }

public:
  AnalysisManager &                     analysisManager_;
  Linear::System &                      linearSystem_;
  Nonlinear::Manager &                  nonlinearManager_;
  Loader::Loader &                      loader_;
  Topo::Topology &                      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;

private:
  Util::OptionBlock     noiseAnalysisOptionBlock_;
  Util::OptionBlock     timeIntegratorOptionBlock_;
  Util::OptionBlock     acLinSolOptionBlock_;
  std::vector<Util::OptionBlock>        dataOptionBlockVec_;

};

// .NOISE
struct NOISEAnalysisReg : public IO::PkgOptionsReg
{
  NOISEAnalysisReg(
    NOISEFactory &   factory )
    : factory_(factory)
  {}

  bool operator()(const Util::OptionBlock &option_block)
  {
    factory_.setNOISEAnalysisOptionBlock(option_block);

    factory_.analysisManager_.addAnalysis(&factory_);

    return true;
  }

  NOISEFactory &              factory_;
};

//-----------------------------------------------------------------------------
// Function      : extractNOISEData
// Purpose       : Extract the parameters from a netlist .NOISE line held in
//                 parsed_line.
// Special Notes : noise specification:
//
//  .noise v(output <,ref>) src ( dec | lin | oct ) pts fstart fstop + <pts_per_summary>
//
// for now ignoring the optional points per summary
//
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
bool
extractNOISEData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("NOISE", Util::OptionBlock::NO_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

  int numFields = parsed_line.size();

  // check for "DATA" first.
  bool dataFound=false;
  int pos1=1;
  while ( pos1 < numFields )
  {
    ExtendedString stringVal ( parsed_line[pos1].string_ );
    stringVal.toUpper ();

    if (stringVal == "DATA")
    {
      dataFound=true;
      break;
    }
    ++pos1;
  }


  // Check that the minimum required number of fields are on the line.
  if ((!dataFound && numFields < 10) || (dataFound && numFields < 9))
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << ".NOISE line has an unexpected number of fields.  NumFields = " << numFields;
    return false;
  }

  int linePosition = 1;   // Start of parameters on .param line.

  Util::Param parameter("", "");

  // output node(s):
  ExtendedString stringVal("");
  std::ostringstream msg;
  int p_err=0;
  if( parsed_line[linePosition].string_ == "V" || parsed_line[linePosition].string_ == "v")
  {
    if( parsed_line[linePosition+3].string_ == ")" )
    {
      stringVal = parsed_line[linePosition].string_;
      stringVal.toUpper();
      parameter.setTag(stringVal);
      parameter.setVal( 1.0 );
      option_block.addParam( parameter );

      stringVal = parsed_line[linePosition+2].string_;
      stringVal.toUpper();
      parameter.setTag( stringVal );
      parameter.setVal( 0.0 );
      option_block.addParam( parameter );

      linePosition += 4;
    }
    else if( parsed_line[linePosition+5].string_ == ")" )
    {
      stringVal = parsed_line[linePosition].string_;
      stringVal.toUpper();
      parameter.setTag(stringVal);
      parameter.setVal( 2.0 );
      option_block.addParam( parameter );

      stringVal = parsed_line[linePosition+2].string_;
      stringVal.toUpper();
      parameter.setTag( stringVal );
      parameter.setVal( 0.0 );
      option_block.addParam( parameter );

      stringVal = parsed_line[linePosition+4].string_;
      stringVal.toUpper();
      parameter.setTag( stringVal );
      parameter.setVal( 0.0 );
      option_block.addParam( parameter );

      linePosition += 6;
    }
    else
    {
      msg << "Unrecognized parenthetical specification for NOISE output ";
      p_err = linePosition;
    }
  }
  else
  {
    msg << "Incorrect format for NOISE output.";
    p_err = linePosition;
  }

  // source is required
  stringVal = parsed_line[linePosition].string_;
  stringVal.toUpper();
  parameter.setTag( "SOURCE" );
  parameter.setVal(std::string(stringVal));
  option_block.addParam( parameter );
  ++linePosition;     // Advance to next parameter.

  // type is required
  stringVal = parsed_line[linePosition].string_;
  stringVal.toUpper();
  parameter.setTag( "TYPE" );
  parameter.setVal(std::string(stringVal));
  option_block.addParam( parameter );
  ++linePosition;     // Advance to next parameter.

  if (dataFound)
  {
    // handle DATA=<name> format
    ++linePosition;  // skip over the = sign
    parameter.setTag( "DATASET" );
    parameter.setVal( parsed_line[ linePosition ].string_ );
    option_block.addParam( parameter );
  }
  else
  {
    // handle format of <points value> <start frequency value> <end frequency value>
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
  }

  // pts_per_summary is optional.  If value is negative, assume it wasn't set.
  parameter.setTag( "PTS_PER_SUMMARY" );
  if ( (!dataFound && numFields >= 11) || (dataFound && numFields >= 10) )
  {
    ++linePosition;     // Advance to next parameter.
    parameter.setVal( parsed_line[linePosition].string_ );
  }
  else
  {
    parameter.setVal( std::string("-1") );
  }

  option_block.addParam( parameter );
  circuit_block.addOptions(option_block);

  return true;
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : outputNoiseNameFile
// Purpose       : Output the device names (and their noise sources) that can be
//                 used with the DNI() and DNO() operators for a specified netlist
// Special Notes : This is a free function so it can be called from N_CIR_Xyce.C
//                 without fully instantiating the Analysis Manager.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 2/12/18
//-----------------------------------------------------------------------------
bool outputNoiseNameFile(
  Parallel::Machine     comm,
  const std::string &   path,
  Loader::Loader &      circuitLoader)
{
  // make the noiseDataVec, in the same fashion as in the
  // NOISE constructor.  This is not the private variable for
  // the NOISE object though.
  int numNoiseDevices = circuitLoader.getNumNoiseDevices();
  std::vector<Xyce::Analysis::NoiseData*> noiseDataVec;
  noiseDataVec.resize(numNoiseDevices);
  for (int i=0;i<numNoiseDevices;++i)
  {
    noiseDataVec[i] = new Analysis::NoiseData();
  }
  circuitLoader.setupNoiseSources(noiseDataVec);

  // get names from noiseDataVec and pass them into an ostringstream
  std::ostringstream oss;
  for (int i=0;i<numNoiseDevices;++i)
  {
    // noiseDataVec_[i]->deviceName is the individual device name (e.g., Q1)
    oss << noiseDataVec[i]->deviceName << std::endl;

    // need to remove the duplicate noise source names, if any.  The sort-unique-erase
    // sequence used is "parallel safe" because each entry (e.g., noiseDataVec[i])
    // exists on only one processor.
    std::sort(noiseDataVec[i]->noiseNames.begin(),noiseDataVec[i]->noiseNames.end());
    std::vector<std::string>::iterator end_unique = std::unique(noiseDataVec[i]->noiseNames.begin(),
                                                          noiseDataVec[i]->noiseNames.end());
    noiseDataVec[i]->noiseNames.erase(end_unique,noiseDataVec[i]->noiseNames.end());

    // output the noiseNames vector, that has had any duplicate entries removed
    for (int j=0;j<noiseDataVec[i]->noiseNames.size();++j)
    {
      // noiseDataVec_[i]->noiseNames[j] are the names of the noise sources, prefaced
      // by noise_<device_name>_
      ExtendedString outStr =  noiseDataVec[i]->noiseNames[j];
      std::string prefixStr = "noise_" + noiseDataVec[i]->deviceName + "_";
      if (outStr.find(prefixStr) != std::string::npos)
      {
        // Remove the prefixStr from outStr, and also make it lower case. We also need
        // to remove any "bad characters" like ( or ) that would interfere with the
        // creation of the DNO and DNI operators.  That "removal" is also done in the
        // NOISE:NOISE() constructor before insertion of the entries into the
        // symbol table.
        outStr.erase(0,prefixStr.length());
        oss << noiseDataVec[i]->deviceName << "\t" << outStr.toLower().removeBadChars() << std::endl;
      }
    }
  }

  // Output entries, in a way that works in parallel.
  // Same approach is used for -namesfile output.
  Util::Marshal mout;
  mout << oss.str();

  std::vector<std::string> dest;
  Parallel::GatherV(comm, 0, mout.str(), dest);

  if (Parallel::rank(comm) == 0)
  {
    std::ofstream output_stream(path.c_str(), std::ios_base::out);

    if (output_stream.fail())
    {
      Report::UserWarning() << "Unable to open noise names file" <<std::endl;
    }
    else
    {
      output_stream << "Valid DNO() and DNI() operators for this netlist." << std::endl
                    << "Format is DNO(deviceName) for the total output noise for a device" << std::endl
                    << "or DNO(deviceName, noiseSource) for an individual output noise contribution." << std::endl
                    << "If a device (e.g., R1) only has one noise source then use DNO(R1) to" << std::endl
                    << "get the total output noise." << std::endl
                    << std::endl
                    << "deviceName \t noiseType" << std::endl;

     for (int p = 0; p < Parallel::size(comm); ++p)
     {
       Util::Marshal min(dest[p]);

        std::string s;
        min >> s;
        output_stream << s;
      }
    }
  }

  // clean up noiseDataVec
  for (int i=0;i<numNoiseDevices;++i)
  {
    delete (noiseDataVec[i]);
  }
  noiseDataVec.clear();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : NOISEFactory::registerFactory
// Purpose       :
//-----------------------------------------------------------------------------
bool
registerNOISEFactory(
  FactoryBlock &        factory_block)
{
  NOISEFactory *factory = new NOISEFactory(factory_block.analysisManager_, factory_block.linearSystem_, factory_block.nonlinearManager_, factory_block.loader_, factory_block.topology_, factory_block.initialConditionsManager_);

  addAnalysisFactory(factory_block, factory);

  factory_block.optionsManager_.addCommandParser(".NOISE", extractNOISEData);

  factory_block.optionsManager_.addCommandProcessor("NOISE", new NOISEAnalysisReg(*factory));

  factory_block.optionsManager_.addOptionsProcessor("TIMEINT", IO::createRegistrationOptions<NOISEFactory>(*factory, &NOISEFactory::setTimeIntegratorOptionBlock));

  factory_block.optionsManager_.addOptionsProcessor("LINSOL-AC", IO::createRegistrationOptions(*factory, &NOISEFactory::setACLinSolOptionBlock));

  factory_block.optionsManager_.addCommandProcessor("DATA",
    IO::createRegistrationOptions(*factory, &NOISEFactory::setDotDataBlock) );

  return true;
}

} // namespace Analysis
} // namespace Xyce

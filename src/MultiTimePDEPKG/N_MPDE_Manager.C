//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Purpose       :
// Special Notes :
// Creator       : Eric R. Keiter, 9233, Computational Sciences
// Creation Date : 03/12/04
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_DCSweep.h>
#include <N_ANP_Transient.h>
#include <N_DEV_DeviceMgr.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_ActiveOutput.h>
#include <N_IO_CmdParse.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_PrintTypes.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_Builder.h>
#include <N_LAS_PrecondFactory.h>
#include <N_LAS_System.h>
#include <N_LAS_Vector.h>
#include <N_LOA_Loader.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_MPDE_Builder.h>
#include <N_MPDE_Discretization.h>
#include <N_MPDE_Loader.h>
#include <N_MPDE_Manager.h>
#include <N_MPDE_SawtoothLoader.h>
#include <N_NLS_Manager.h>
#include <N_PDS_Manager.h>
#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TIAParams.h>
#include <N_TIA_WorkingIntegrationMethod.h>
#include <N_TOP_Topology.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Param.h>
#include <N_UTL_Timer.h>

#include <Teuchos_RCP.hpp>
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

using Xyce::DEBUG_ANALYSIS;
using Xyce::DEBUG_MPDE;
using Xyce::VERBOSE_TIME;


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::N_MPDE_Manager
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
N_MPDE_Manager::N_MPDE_Manager(
  Xyce::Analysis::AnalysisManager &     analysis_manager,
  Xyce::Loader::Loader &                loader,
  Xyce::Device::DeviceMgr &             device_manager,
  Xyce::Linear::Builder &               builder,
  Xyce::Topo::Topology &                topology,
  Xyce::IO::InitialConditionsManager &  initial_conditions_manager,
  Xyce::IO::RestartMgr &                restart_manager,
  const Xyce::IO::CmdParse &            command_line)
 : commandLine_(command_line),
   analysisManager_(analysis_manager),
   nonlinearManager_(commandLine_),
   mpdeLinearSystem_(),
   deviceManager_(device_manager),
   pdsManager_(*analysis_manager.getPDSManager()),
   appBuilder_(builder),
   topology_(topology),
   initialConditionsManager_(initial_conditions_manager),
   restartManager_(restart_manager),
   appLoader_(loader),
   tiaMPDEParams_(),
   timeIntegratorOptionBlock_(),
   outputAdapter_(analysis_manager.getOutputManagerAdapter(), analysis_manager, *this),
   mpdeState_(),
   mpdeLoaderPtr_(0),
   mpdeBuilderPtr_(0),
   mpdeDiscPtr_(0),
   dcOpSolVecPtr_(0),
   dcOpStateVecPtr_(0),
   dcOpQVecPtr_(0),
   dcOpStoreVecPtr_(0),
   mpdeICVectorPtr_(0),
   mpdeICStateVectorPtr_(0),
   mpdeICQVectorPtr_(0),
   mpdeICStoreVectorPtr_(0),
   test_(false),
   size_(21),
   tranRunForSize_(false),
   maxCalcSize_(20),
   maxCalcSizeGiven_(false),
   NOOP_(false), 
   fastSrcGiven_(false),
   oscOut_(""),
   oscOutGiven_(false),
   nonLteSteps_(0),
   nonLteStepsGiven_(false),
   period_(1.0),
   periodGiven_(false),
   startUpPeriods_(0),
   startUpPeriodsGiven_(false),
   saveIcData_(false),
   transientNeedsToLoadInitialConditionsAndInitializeProblem_(true),
   transientNowCanOutputTrueMPDEResults_(false),
   fastTimes_(),
   freqPoints_(),
   fastTimeDisc_(0),
   fastTimeDiscOrder_(1),
   warpMPDE_(false),
   warpMPDEPhasePtr_(0),
   warpMPDEOSCOUT_(-1),
   warpPhase_(0),
   warpPhaseGiven_(false),
   warpPhaseCoeff_(0.0),
   warpPhaseCoeffGiven_(false),
   fftFlag_(false),
   icPer_(10),
   initialCondition_(0),
   outputInterpMPDE_(true),
   warpMPDEICFlag_(false),
   dcopExitFlag_(false),
   icExitFlag_(false),
   exitSawtoothStep_(-1)
{
  tiaMPDEParams_.setMaxOrder(command_line.getArgumentIntValue("-maxord", tiaMPDEParams_.maxOrder));
}

N_MPDE_Manager::~N_MPDE_Manager()
{
  delete mpdeLoaderPtr_;
  delete dcOpSolVecPtr_;
  delete dcOpStateVecPtr_;
  delete dcOpQVecPtr_;
  delete dcOpStoreVecPtr_;
  delete mpdeICVectorPtr_;
  delete mpdeICStateVectorPtr_;
  delete mpdeICQVectorPtr_;
  delete mpdeICStoreVectorPtr_;
  delete mpdeBuilderPtr_;
  delete mpdeDiscPtr_;
  delete warpMPDEPhasePtr_;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::finalExpressionBasedSetup
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/4/2021
//-----------------------------------------------------------------------------
void N_MPDE_Manager::finalExpressionBasedSetup()
{
  // This performs the early setup for the outputters for this analysis mode.
  (analysisManager_.getOutputManagerAdapter().getOutputManager()).earlyPrepareOutput
      (pdsManager_.getPDSComm()->comm(), Xyce::Analysis::ANP_MODE_TRANSIENT);

  (analysisManager_.getOutputManagerAdapter().getOutputManager()).earlyPrepareOutput
      (pdsManager_.getPDSComm()->comm(), Xyce::Analysis::ANP_MODE_MPDE);

}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::run
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool
N_MPDE_Manager::run(
  Xyce::Linear::System &        linear_system,
  Xyce::Nonlinear::Manager &    nonlinear_manager,
  Xyce::Topo::Topology &        topology)
{
  bool returnValue = true;

  initializeOscOut(topology);

  transientNeedsToLoadInitialConditionsAndInitializeProblem_ = true;
  transientNowCanOutputTrueMPDEResults_ = false;

  returnValue = initializeAll(linear_system, nonlinear_manager);

  deviceManager_.setMPDEFlag(false);
  bool ret1 = runInitialCondition(linear_system, nonlinear_manager);
  returnValue = returnValue && ret1;

  deviceManager_.setMPDEFlag(true);

  //Turn off voltage limiting when in MPDE phase. This is a temporary
  //code and should be removed once voltage limiting for MPDE is developed.
  deviceManager_.setVoltageLimiterFlag (false);

  bool ret2 = setupMPDEProblem_();
  returnValue = returnValue && ret2;

  bool ret3 = runMPDEProblem_();
  returnValue = returnValue && ret3;

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::setMPDEAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool
N_MPDE_Manager::setMPDEAnalysisParams(
  const Xyce::Util::OptionBlock &       tiaParamsBlock)
{
  for (Xyce::Util::ParamList::const_iterator it = tiaParamsBlock.begin(), end = tiaParamsBlock.end(); it != end; ++it)
  {
    const Xyce::Util::Param &param = *it;

    if (tiaMPDEParams_.setAnalysisOption(param))
      ;
    else if (param.uTag() == "NOOP" ||
             param.uTag() == "UIC")
    {
      NOOP_ = true;
    }
    else
      Xyce::Report::UserError() << param.uTag() << " is not a recognized analysis option";
  }

  if (tiaMPDEParams_.finalTime <= tiaMPDEParams_.initialOutputTime || tiaMPDEParams_.finalTime <= 0 || tiaMPDEParams_.initialOutputTime < 0)
  {
    Xyce::Report::UserFatal0()
      << "Final time of " << tiaMPDEParams_.finalTime
      << " is earlier or same as start time of "
      << tiaMPDEParams_.initialOutputTime
      << " Check netlist for invalid .MPDE specification ";
  }

  // if starting time steps is baloney, set then to a default value.
  // ERK.  This trap is redudant with an identical one in step error
  // control.
  if ( tiaMPDEParams_.initialTimeStep <= 0.0 )
  {
    tiaMPDEParams_.initialTimeStep = 1.0e-10;
  }

  if (DEBUG_ANALYSIS && Xyce::isActive(Xyce::Diag::TIME_PARAMETERS))
  {
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  Transient simulation parameters" << std::endl
                 << "  initial time = " << tiaMPDEParams_.initialTime << std::endl
                 << "  final time   = " << tiaMPDEParams_.finalTime << std::endl
                 << "  starting time step = " << tiaMPDEParams_.initialTimeStep << std::endl
                 << "  initial output time = " << tiaMPDEParams_.initialOutputTime << std::endl
                 << "  NOOP/UIC is " << (NOOP_ ? "" : "NOT") << " set" << std::endl
                 << Xyce::section_divider << std::endl;
  }


  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::setMPDEOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool
N_MPDE_Manager::setMPDEOptions(
  const Xyce::Util::OptionBlock & option_block)
{
  //MOST OF THIS SHOULD BE PUSHED TO INITCOND AND PROBLEM CLASSES
  for (Xyce::Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; ++it)
  {
    bool value_set = false;

    if (it->uTag() == "OSCOUT")
    {
      oscOut_ = it->stringValue();
      oscOutGiven_ = true;
      value_set = true;
    }

    if (!value_set)
    {
      value_set = Xyce::Util::setValue(*it, "N2", size_)
      || Xyce::Util::setValue(*it, "AUTON2", tranRunForSize_)
      || Xyce::Util::setValue(*it, "AUTON2MAX", maxCalcSize_, maxCalcSizeGiven_)
      || Xyce::Util::setValue(*it, "NONLTESTEPS", nonLteSteps_, nonLteStepsGiven_)
      || Xyce::Util::setValue(*it, "STARTUPPERIODS", startUpPeriods_, startUpPeriodsGiven_)
      || Xyce::Util::setValue(*it, "SAVEICDATA", saveIcData_)
      || Xyce::Util::setValue(*it, "PHASE", warpPhase_, warpPhaseGiven_)
      || Xyce::Util::setValue(*it, "PHASECOEFF", warpPhaseCoeff_, warpPhaseCoeffGiven_)
      || Xyce::Util::setValue(*it, "TEST", test_)
      || Xyce::Util::setValue(*it, "T2", period_, periodGiven_)
      || Xyce::Util::setValue(*it, "WAMPDE", warpMPDE_)
      || Xyce::Util::setValue(*it, "FREQDOMAIN", fftFlag_)
      || Xyce::Util::setValue(*it, "ICPER", icPer_)
      || Xyce::Util::setValue(*it, "IC", initialCondition_)
      || Xyce::Util::setValue(*it, "DIFF", fastTimeDisc_)
      || Xyce::Util::setValue(*it, "DIFFORDER", fastTimeDiscOrder_)
      || Xyce::Util::setValue(*it, "DCOPEXIT", dcopExitFlag_)
      || Xyce::Util::setValue(*it, "ICEXIT", icExitFlag_)
      || Xyce::Util::setValue(*it, "EXITSAWTOOTHSTEP", exitSawtoothStep_)
      || Xyce::Util::setValue(*it, "OSCKICKOFF", warpMPDEICFlag_)
      || Xyce::Util::setValue(*it, "DEBUGLEVEL", Xyce::setMPDEDebugLevel)
      || Xyce::Util::setValue(*it, "OSCSRC", srcVec_, fastSrcGiven_);
    }

    if (value_set)
      ;
    else
    {
      Xyce::Report::UserWarning() << " Unrecognized MPDEINT option: " << (*it).tag();
    }
  }

  if( (fastTimeDisc_ == 1) && (fastTimeDiscOrder_ < 2 ) )
  {
    // if user asked for centered differences, then default to
    // second order if none was specified otherwise the Discretization
    // class will return backwards differences which may not be what
    // the user expects
    fastTimeDiscOrder_ = 2;
  }

  printParams_ ();

  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerTranMPDEOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 08/20/07
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::setTransientOptions(const Xyce::Util::OptionBlock &option_block)
{
  timeIntegratorOptionBlock_ = option_block;

  for (Xyce::Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; ++it)
  {
    const Xyce::Util::Param &param = (*it);

    if (param.uTag() == "OUTPUTINTERPMPDE")
      outputInterpMPDE_ = static_cast<bool> (param.getImmutableValue<int>());
  }
  
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::setLinSolOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 03/13/2017
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::setLinSolOptions(const Xyce::Util::OptionBlock &option_block)
{
  linSolOptionBlock_ = option_block;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::printParams_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 03/30/04
//-----------------------------------------------------------------------------
void N_MPDE_Manager::printParams_ ()
{
  Xyce::dout() << "\n" << std::endl;
  Xyce::dout() << Xyce::section_divider << std::endl;
  Xyce::dout() <<
                         "\n***** MPDE options:\n" << std::endl;

  Xyce::dout() <<
       "\tN2:\t\t\t" <<  size_ << std::endl;

  std::string msg;
  if (fastSrcGiven_)
  {
    for (unsigned int i1=0;i1<srcVec_.size();++i1)
    {
      msg = "\toscsrc:\t\t\t"+srcVec_[i1];
      Xyce::dout() <<msg << std::endl;
    }
  }
  else
  {
    msg = "\toscsrc:\t\t\tNot Specified";
    Xyce::dout() <<msg << std::endl;
  }

  if (oscOutGiven_)
  {
    msg = "\toscout:\t\t\t"+oscOut_;
  }
  else
  {
    msg = "\toscout:\t\t\tNot Specified";
  }
  Xyce::dout() <<msg << std::endl;

  Xyce::dout() <<
       "\tT2:\t\t\t" <<  period_ << std::endl;

  if( fastTimeDisc_ == 0 )
    msg = "\tdiscretization:\tBackward";
  else if( fastTimeDisc_ == 1 )
    msg = "\tdiscretization:\tCentered";
  else
    msg = "\tdiscretization:\tUNKNOWN";
  Xyce::dout() <<msg << std::endl;

  Xyce::dout() <<
       "\tDisc Order:\t\t" <<  fastTimeDiscOrder_ << std::endl;

  if (test_)
  {
    Xyce::dout() <<
      "\tTest Mode:" << std::endl;
  }
  else
  {
    Xyce::dout() <<
      "\tFull MPDE Mode(not test mode)" << std::endl;
  }

  if (warpMPDE_)
  {
    Xyce::dout() <<
      "\tWarpedMPDE:\t\tON" << std::endl;
  }
  else
  {
    Xyce::dout() <<
      "\tWarpedMPDE:\t\tOFF" << std::endl;
  }

  if (fftFlag_)
  {
    Xyce::dout() <<
      "\tFrequency domain:\tON" << std::endl;
  }
  else
  {
    Xyce::dout() <<
      "\tFrequency domain:\tOFF" << std::endl;
  }

  Xyce::dout() <<
       "\tICPER:\t\t\t" <<  icPer_ << std::endl;

  if (initialCondition_ == MPDE_IC_DCOP)
  {
    Xyce::dout() <<
       "\tInitial Condition:\tDCOP" << std::endl;
  }
  else if (initialCondition_ == MPDE_IC_SAWTOOTH)
  {
    Xyce::dout() <<
       "\tInitial Condition:\tSAWTOOTH" << std::endl;
  }

  Xyce::dout() << "\n" << std::endl;
  Xyce::dout() << Xyce::section_divider << std::endl;
  Xyce::dout() << "\n" << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::initializeAll
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool
N_MPDE_Manager::initializeAll(
  Xyce::Linear::System &        linear_system,
  Xyce::Nonlinear::Manager &    nonlinear_manager)
{
  bool returnValue = true;

  if (DEBUG_MPDE)
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  N_MPDE_Manager::initializeAll" << std::endl;

  if (DEBUG_MPDE)
    Xyce::dout() << "N_MPDE_Manager::initializeAll[Construct Builder]\n";

  std::vector<double> srcPeriods;
  srcPeriods = deviceManager_.getFastSourcePeriod(pdsManager_.getPDSComm()->comm(), srcVec_ );

  // If running straight MPDE, the only option is to get the period
  // from the "oscsrc", or fast source.  If running warped, it is still
  // possible to get it from a fast source, but the T2 (period_) parameter
  // takes precedence.
  if (periodGiven_ && fastSrcGiven_)
  {
    Xyce::dout() << "Note:  oscsrc is being ignored, as the user has set T2" << std::endl;
  }

  if (fastSrcGiven_  && !periodGiven_)
  {
    // don't assume the periods are all the same:
    period_ = srcPeriods[0];
    int numSrc = srcPeriods.size();
    for( int i=1; i<numSrc; i++ )
    {
      if( srcPeriods[i] > period_ )
      {
        period_ = srcPeriods[i];
      }
    }
  }

  if (warpMPDE_ && !periodGiven_ && !fastSrcGiven_)
  {
    Xyce::Report::UserError() <<"The Warped MPDE needs for the user to set oscsrc or T2 (preferably T2)";
    return false;
  }

  if (!warpMPDE_ && !fastSrcGiven_)
  {
    Xyce::Report::UserError() << "The MPDE algorithm needs the user to set oscsrc";
    return false;
  }

  Xyce::TimeIntg::TIAParams tia_params;

  // set DAE initial time = 0.0
  tia_params.initialTime = 0.0;
  tia_params.finalTime = startUpPeriods_ * period_;
  analysisManager_.getStepErrorControl().pauseTime = tia_params.finalTime;

  // If it was requested, advance the solution a fixed number of startup periods
  if (startUpPeriodsGiven_)
  {
    bool startupPeriodsSuccess = runStartupPeriods(tia_params, linear_system, nonlinear_manager);
    if (!startupPeriodsSuccess)
    {
      Xyce::Report::UserError() << "Failed to calculate the startup periods";
      return false;
    }
    returnValue = returnValue && startupPeriodsSuccess;

    // tell mpde to start after this startup period
    tia_params.initialTime = tia_params.finalTime;
  }

  tiaMPDEParams_.initialTime = tia_params.initialTime;
  
  if (!tranRunForSize_)
  {
    // Simple discrete time points for fast time scale
    fastTimes_.resize(size_+1);
    double fastStep = period_/static_cast<double>(size_);
    int i=0;
    for( i = 0; i < size_; ++i )
    {
      fastTimes_[i] = static_cast<double>(i) * fastStep;
    }
    fastTimes_[size_] = period_;
    //-----------------------------------------
    if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS) )
    {
      Xyce::dout() << "MPDE Fast Step: " << fastStep << std::endl;
      for( i = 0; i < size_; ++i )
        Xyce::dout() << "fastTimes_["<<i<<"] = " << fastTimes_[i] << std::endl;
      Xyce::dout() << Xyce::section_divider << std::endl;
    }
  }
  else
  {
    tia_params.finalTime = tia_params.initialTime + period_;

    if (initialCondition_ == MPDE_IC_TWO_PERIOD)
      tia_params.finalTime += period_;

    analysisManager_.getStepErrorControl().pauseTime = tia_params.finalTime;

    returnValue = runTransientIC(tia_params, linear_system, nonlinear_manager);
    if (!returnValue)
    {
      Xyce::Report::UserError() << "Failed to compute the transient initial condition";
      return false;
    }

    filterFastTimePoints(tia_params);
    tia_params.initialTime = tia_params.finalTime;
  }

  // set the fast source flag in the devices
  srcPeriods = deviceManager_.registerFastSources(pdsManager_.getPDSComm()->comm(), srcVec_);

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::initializeMPDEAll
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Genie Hsieh, 00434
// Creation Date : 11/14/2013
//-----------------------------------------------------------------------------
bool
N_MPDE_Manager::initializeMPDEAll()
{
  bool returnValue = true;

  mpdeDiscPtr_ = new N_MPDE_Discretization(static_cast<N_MPDE_Discretization::Type>(fastTimeDisc_), fastTimeDiscOrder_);

  mpdeBuilderPtr_ = new N_MPDE_Builder(*this, size_, *mpdeDiscPtr_, warpMPDE_);

  if (DEBUG_MPDE)
    std::cout << "N_MPDE_Manager::initializeMPDEAll[Generate Maps]\n";

  mpdeBuilderPtr_->generateMaps( rcp( pdsManager_.getParallelMap( Xyce::Parallel::SOLUTION ), false) );
  mpdeBuilderPtr_->generateStateMaps( rcp( pdsManager_.getParallelMap( Xyce::Parallel::STATE ), false) );
  mpdeBuilderPtr_->generateStoreMaps( rcp( pdsManager_.getParallelMap( Xyce::Parallel::STORE ), false) );
  mpdeBuilderPtr_->generateLeadCurrentMaps( rcp( pdsManager_.getParallelMap( Xyce::Parallel::LEADCURRENT ), false) );


  // Setup WaMPDE Phase Condition Object
  if (warpMPDE_)
  {

    // These three values are initialized during generateMaps in the MPDE Builder.
    int offset = mpdeBuilderPtr_->getMPDEOffset();
    int omegaGID = mpdeBuilderPtr_->getMPDEomegaGID();
    int phiGID = mpdeBuilderPtr_->getMPDEphiGID();
    int augProcID = mpdeBuilderPtr_->getMPDEaugProcID();

    if (DEBUG_MPDE)
      std::cout << "N_MPDE_Manager::initializeMPDEAll[ offset = " << offset << ", omegaGID = " << omegaGID
                << ", phiGID = " << phiGID << ", augProcID = " << augProcID << " ] " << std::endl;

    warpMPDEPhasePtr_ = new N_MPDE_WarpedPhaseCondition(warpPhase_, warpPhaseCoeff_, warpMPDEOSCOUT_, omegaGID, offset, size_, phiGID, augProcID);
    mpdeBuilderPtr_->setWarpedPhaseCondition(warpMPDEPhasePtr_);
  }

  // Now that the builder is set up, allocate petra vectors, etc.
  mpdeICVectorPtr_ = dynamic_cast<Xyce::Linear::BlockVector *>(mpdeBuilderPtr_->createVector());
  mpdeICStateVectorPtr_ = dynamic_cast<Xyce::Linear::BlockVector *>(mpdeBuilderPtr_->createStateVector());
  mpdeICQVectorPtr_ = dynamic_cast<Xyce::Linear::BlockVector *>(mpdeBuilderPtr_->createVector());
  mpdeICStoreVectorPtr_ = dynamic_cast<Xyce::Linear::BlockVector *>(mpdeBuilderPtr_->createStoreVector());

  if (DEBUG_MPDE)
    std::cout << "N_MPDE_Manager::initializeMPDEAll[Finished]\n";

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::initializeOscOut
// Purpose       : Sets up the warpMPDEOSCOUT_  global id (GID) index.
//
// Special Notes : This is called from initializeAll.
//
//   This has to be done here  (rather than setMPDEOptions) because
//   it can only be done after the non-MPDE topology has been set up.
//
// Scope         : private
// Creator       : Eric Keiter, 1437, Computational Sciences
// Creation Date : 12/16/07
//-----------------------------------------------------------------------------
bool
N_MPDE_Manager::initializeOscOut(
  Xyce::Topo::Topology &        topology)
{
  // Set up the gid index for the osc out variable

  if (oscOutGiven_)
  {
    bool oscOutError = false;

    std::string varName = oscOut_;
    Xyce::Util::toUpper(varName);

    std::vector<int> oscoutList, dummyList;
    char type1;

    int pos1 = varName.find("I(");
    int pos2 = varName.find("V(");
    int pos3 = varName.find(")");

    if (pos1!=(int)std::string::npos && pos3!=(int)std::string::npos) // this is a current variable
    {
      pos1+=2;
      int len = pos3-pos1;
      std::string tmpVar (varName,pos1,len);
      if (DEBUG_MPDE)
        Xyce::dout() << "tmpVar (for I-oscout) = " << tmpVar << std::endl;

      topology.getNodeSVarGIDs(Xyce::NodeID(tmpVar, Xyce::_DNODE), oscoutList, dummyList, type1);
    }
    else if (pos2 !=(int)std::string::npos && pos3!=(int)std::string::npos) // this is a voltage variable
    {
      pos2+=2;
      int len = pos3-pos2;

      std::string tmpVar (varName,pos2,len);
      if (DEBUG_MPDE)
        Xyce::dout() << "tmpVar (for V-oscout) = " << tmpVar << std::endl;

      topology.getNodeSVarGIDs(Xyce::NodeID(tmpVar, Xyce::_VNODE), oscoutList, dummyList, type1);
    }
    else
    {
      // do nothing, assume that the varName from the netlist requires no
      // modification.
      topology.getNodeSVarGIDs(Xyce::NodeID(varName,-1), oscoutList, dummyList, type1);
    }

    if (oscoutList.size()==1)
    {
      warpMPDEOSCOUT_=oscoutList.front();  // This is the GID of the OSCOUT.
      //oscOutGiven_ = true;

      if (warpMPDEOSCOUT_>=0)
        Xyce::dout() << "warpMPDEOSCOUT = " << warpMPDEOSCOUT_ << std::endl;

    }
#ifndef Xyce_PARALLEL_MPI
    else  // This is not an error in parallel.
    {
      oscOutError = true;
    }
#else
    int tmpOscOut = warpMPDEOSCOUT_;
    pdsManager_.getPDSComm()->maxAll( &tmpOscOut, &warpMPDEOSCOUT_, 1 );
    if ( warpMPDEOSCOUT_ < 0 )
    {
      oscOutError = true;
    }
#endif

    if (oscOutError)
    {
      Xyce::Report::UserWarning() << "Unrecognized value for MPDE option oscOut:  "<< oscOut_;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::runInitialCondition
// Purpose       : Execution of initial condition phase
// Special Notes : Produces an "MPDE" size initial value
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
bool
N_MPDE_Manager::runInitialCondition(
  Xyce::Linear::System &        linear_system,
  Xyce::Nonlinear::Manager &    nonlinear_manager)
{
  bool returnValue=true;
  int goodStartPoint = 0;
  int goodEndPoint = 0;
  bool foundStart = false;
  bool foundEnd = false;
  double shiftTimeBy = 0;

  Xyce::lout() << " ***** Running MPDE initial conditions....\n" << std::endl;
  if (DEBUG_MPDE)
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  N_MPDE_Manager::runInitialCondition" << std::endl;

  int n2 = size_;

  Xyce::TimeIntg::DataStore * dsPtr = analysisManager_.getDataStore();

  if (!warpMPDE_)
  {
    if (!dcOpSolVecPtr_)
    {
      // haven't done a dc op yet, so do it now
      runDCOP(linear_system, nonlinear_manager);
    }

    initializeMPDEAll();

    // place dc op values in mpde block vectors
    mpdeICVectorPtr_->block(0) = *dcOpSolVecPtr_;
    mpdeICStateVectorPtr_->block(0) = *dcOpStateVecPtr_;
    mpdeICQVectorPtr_->block(0) = *dcOpQVecPtr_;
    mpdeICStoreVectorPtr_->block(0) = *dcOpStoreVecPtr_;
  }
  else
  {
    if (initialCondition_ == WAMPDE_IC_TRANSIENT)
    {
      std::cout << "OSCOUT's GID = " << warpMPDEOSCOUT_ << std::endl;
      std::cout << "dsPtr->fastTimeSolutionVec[0] = " << std::endl;
      (dsPtr->fastTimeSolutionVec[0])->print(Xyce::dout());
      Xyce::Linear::Vector &firstSolVec = *dsPtr->fastTimeSolutionVec[0];
      Xyce::Linear::Vector &nextSolVec = *dsPtr->fastTimeSolutionVec[1];
      std::cout << "time=0 OSCOUT solution = " << firstSolVec[warpMPDEOSCOUT_] << std::endl;
      std::cout << "time=1 OSCOUT solution = " << nextSolVec[warpMPDEOSCOUT_] << std::endl;


      // Originally here saves the solution at the end of startupperiod to the 1st location in ICvector
      // Improve the IC so only when currSolution = 0 the solution is saved as the start of the IC.
      // Besides, improve the WaMPDE IC by checking the last solution stored (goodEndPoint) in MPDE IC Vector
      // has the same value (i.e. solution=0) and slope as the first solution (goodStartPoint).
      // This is to ensure we have a perfect WaMPDE IC.
      for (int i=0 ; i<n2-1 ; ++i)
      {
        firstSolVec = *dsPtr->fastTimeSolutionVec[i];
        nextSolVec = *dsPtr->fastTimeSolutionVec[i+1];

        std::cout << "time=" << i << ", OSCOUT solution = " << firstSolVec[warpMPDEOSCOUT_] << std::endl;

        if (!foundStart && firstSolVec[warpMPDEOSCOUT_] <=warpPhaseCoeff_ && nextSolVec[warpMPDEOSCOUT_] > warpPhaseCoeff_)
        {
          goodStartPoint = i;
          foundStart = true;
        }
        else if (foundStart && firstSolVec[warpMPDEOSCOUT_] <=warpPhaseCoeff_ && nextSolVec[warpMPDEOSCOUT_] > warpPhaseCoeff_)
        {
          goodEndPoint = i;
          foundEnd = true;
          break;
        }
      }


      if (!foundEnd)
      {
        Xyce::Report::UserError() << "Not enough fast time points to find a good WaMPDE initial condition\n"
                            <<"Try to icrease the fast time period, T2, to fix this error.";
        return false;
      }

      std::cout << "At time = " << goodStartPoint << ", oscout solution is approaching " << warpPhaseCoeff_ << std::endl;
      std::cout << "At time = " << goodEndPoint << ", oscout solution is again approaching " << warpPhaseCoeff_ << std::endl;

      size_ = goodEndPoint - goodStartPoint;
      initializeMPDEAll(); // initialize MPDE IC vectors to have
      // the updated size of the good initial condition


      mpdeICVectorPtr_->block(0) = *(dsPtr->fastTimeSolutionVec[indicesUsed_[goodStartPoint]]);
      mpdeICStateVectorPtr_->block(0) = *(dsPtr->fastTimeStateVec[indicesUsed_[goodStartPoint]]);
      mpdeICQVectorPtr_->block(0) = *(dsPtr->fastTimeQVec[indicesUsed_[goodStartPoint]]);
      mpdeICStoreVectorPtr_->block(0) = *(dsPtr->fastTimeStoreVec[indicesUsed_[goodStartPoint]]);

      std::cout << "mpdeIC vectors blockCount = " << mpdeICVectorPtr_->blockCount() << std::endl;
    } //end if initialCondition_ = WAMPDE_IC_TRANSIENT
    else
    {
      initializeMPDEAll();

      mpdeICVectorPtr_->block(0) = *(dsPtr->currSolutionPtr);
      mpdeICStateVectorPtr_->block(0) = *(dsPtr->currStatePtr);
      mpdeICQVectorPtr_->block(0) = *(dsPtr->daeQVectorPtr);
      mpdeICStoreVectorPtr_->block(0) = *(dsPtr->currStorePtr);
    }
  }

  if (DEBUG_MPDE && dcopExitFlag_)
    Xyce::Report::DevelFatal() << "Exit on DCOPEXIT";

  switch (initialCondition_)
  {
    case MPDE_IC_DCOP:
      // 03/23/04 tscoffe:  Here's the meta algorithm for the
      // constant DCOP IC strategy:
      for (int i=1 ; i<n2 ; ++i)
      {
        //std::cout << "Computing the initial condition for IC block " << i << std::endl;
        mpdeState_.fastTime = fastTimes_[i];
        if (!NOOP_)
        {
          // as long as the user didn't request "noop"
          // try and do the operating point calculation
          // try to run the problem
          {
            //Xyce::Analysis::DCSweep dc_sweep(analysisManager_, linear_system, nonlinear_manager, appLoader_, topology_, initialConditionsManager_, 0);
            Xyce::Analysis::DCSweep dc_sweep(analysisManager_, &linear_system, nonlinear_manager, appLoader_, topology_, initialConditionsManager_, 0);
            analysisManager_.pushActiveAnalysis(&dc_sweep);
            dc_sweep.run();

            // // Add in simulation times
            // accumulateStatistics_(dc_sweep);

            // print out analysis info
            Xyce::lout() << " ***** Harmonic Balance Computation Summary *****" << std::endl;
            dc_sweep.printLoopInfo( 0, 0 );

            analysisManager_.popActiveAnalysis();
          }
        }

        mpdeICVectorPtr_->block(i) = *analysisManager_.getDataStore()->currSolutionPtr;
        mpdeICStateVectorPtr_->block(i) = *analysisManager_.getDataStore()->currStatePtr;
        mpdeICQVectorPtr_->block(i) = *analysisManager_.getDataStore()->daeQVectorPtr;
        mpdeICStoreVectorPtr_->block(i) = *analysisManager_.getDataStore()->currStorePtr;
      }

      //std::cout << "Done computing the initial condition" << std::endl;

      break;

    case MPDE_IC_SAWTOOTH:
    {

      // 03/23/04 tscoffe:  Here's the meta algorithm for the
      // "sawtooth" IC strategy:

      //Setting up specialized Loader
      N_MPDE_SawtoothLoader stLoader( mpdeState_, analysisManager_, appLoader_ );
      stLoader.setTimeShift( fastTimes_[0] );

      // get the time integrator params, and save a copy.
      Xyce::TimeIntg::TIAParams tia_params;

      for (int i=1 ; i<n2 ; ++i)
      {
        // set DAE initial time = 0.0
        tia_params.initialTime = 0.0;
        // set DAE final time = h2
        double h2 = fastTimes_[i]-fastTimes_[i-1];
        tia_params.finalTime = h2;

        //tia_params.finalTime = fastTimes_[i];
        analysisManager_.getStepErrorControl().pauseTime = tia_params.finalTime;

        // We need to do an initialization every time...
        analysisManager_.setBeginningIntegrationFlag(true);
        analysisManager_.getStepErrorControl().resetAll(tia_params);
        analysisManager_.getWorkingIntegrationMethod().initialize(tia_params);

        if (DEBUG_MPDE)
          Xyce::dout() << "Beginning SAWTOOTH IC: i = " << i
                       << " initialTime = " << tia_params.initialTime
                       << " finalTime = " << tia_params.finalTime << std::endl;

        // Set t = t + (i-1.0)*h2 in MPDE source in bhat
        stLoader.setTimeShift( fastTimes_[i-1] );

        // set DAE initial condition to mpdeICVectorPtr_[i-1]
        *(analysisManager_.getDataStore()->nextSolutionPtr) = mpdeICVectorPtr_->block(i-1);
        *(analysisManager_.getDataStore()->nextStatePtr) = mpdeICStateVectorPtr_->block(i-1);
        *(analysisManager_.getDataStore()->daeQVectorPtr) = mpdeICQVectorPtr_->block(i-1);
        *(analysisManager_.getDataStore()->nextStorePtr) = mpdeICStoreVectorPtr_->block(i-1);

        // solve the DAE to time h2
        {
          // This prepares the outputters for this analysis mode
          Xyce::IO::ActiveOutput active(analysisManager_.getOutputManagerAdapter().getOutputManager());

          active.setStepSweepVector(analysisManager_.getOutputManagerAdapter().getStepSweepVector());
          active.setDCSweepVector(analysisManager_.getOutputManagerAdapter().getDCSweepVector());
          active.add(pdsManager_.getPDSComm()->comm(), Xyce::Analysis::ANP_MODE_TRANSIENT);

          if (VERBOSE_TIME)
            tia_params.printParams(Xyce::lout(), nonlinearAnalysisMode(Xyce::Analysis::ANP_MODE_TRANSIENT));

          analysisManager_.setAnalysisMode(Xyce::Analysis::ANP_MODE_TRANSIENT);

          //Xyce::Analysis::Transient transient(analysisManager_, linear_system, nonlinear_manager, stLoader, topology_, initialConditionsManager_, restartManager_, &outputAdapter_, 0, 0);
          Xyce::Analysis::Transient transient(analysisManager_, &linear_system, nonlinear_manager, stLoader, topology_, initialConditionsManager_, restartManager_, &outputAdapter_, 0, 0);
          analysisManager_.pushActiveAnalysis(&transient);
          transient.setTIAParams(tia_params);
          transient.setAnalysisParams(Xyce::Util::OptionBlock());
          transient.setNOOP(true);
          analysisManager_.getStepErrorControl().resetAll(tia_params);

          returnValue = transient.run();

          transientNeedsToLoadInitialConditionsAndInitializeProblem_ = false;

          analysisManager_.popActiveAnalysis();
        }

        // put the dsPtr->currentSolutionPtr into mpdeICVectorPtr_[i]
        mpdeICVectorPtr_->block(i) = *analysisManager_.getDataStore()->currSolutionPtr;
        mpdeICStateVectorPtr_->block(i) = *analysisManager_.getDataStore()->currStatePtr;
        mpdeICStoreVectorPtr_->block(i) = *analysisManager_.getDataStore()->currStorePtr;
        if (DEBUG_MPDE)
        {
          Xyce::dout() << "End SAWTOOTH IC: i = " << i <<std::endl;
          if (i == exitSawtoothStep_)
            Xyce::Report::DevelFatal() << "Exiting on exitSawtoothStep_";
        }
      }
    }
    break;

    case WAMPDE_IC_TRANSIENT:
    {
      if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS) )
      {
        for( unsigned int i=0; i<indicesUsed_.size(); i++ )
          Xyce::dout() << "indicesUsed_[ " << i << " ] = " << indicesUsed_[i] << std::endl;
      }

      // we did an initial transient run so pick up ic data from that
      // Need to seperate cases for MPDE & WaMPDE
      int numPoints = 0;
      if (warpMPDE_)
      {
        // Also update the time in ICVector/fastTimeSol here.
        // genie 120913. Update the fastTimes_ here to reflect the times of the "good" WaMPDE initial condition
        // genie 121113. Note: It is important to keep the fastStep the same as before.

        int i = 0;

        std::cout << "121813 debug: goodStart = " << goodStartPoint << std::endl;
        for( i = 0; i < (int)fastTimes_.size(); ++i )
          std::cout << "fastTimes_["<<i<<"] = " << fastTimes_[i] << std::endl;


        shiftTimeBy = fastTimes_[goodStartPoint];
        for( i = 0; i < (goodEndPoint - goodStartPoint); ++i )
        {
          fastTimes_[i] = fastTimes_[goodStartPoint + i] - shiftTimeBy;
        }
        fastTimes_[i] = fastTimes_[goodStartPoint + i] - shiftTimeBy; // i = goodEndPoint-goodStartPoint
        period_ = fastTimes_[i];

        std::cout << "shiftTimeBy = " << shiftTimeBy << ", new period_ (T2) is " << period_ << std::endl;

        numPoints = goodEndPoint - goodStartPoint;
        fastTimes_.resize(numPoints);

        //-----------------------------------------
        for( i = 0; i < (int)fastTimes_.size(); ++i )
          std::cout << "fastTimes_["<<i<<"] = " << fastTimes_[i] << std::endl;

      }
      else
        numPoints = n2;


      for (int i=0 ; i < numPoints ; i++)
      {

        Xyce::TimeIntg::DataStore *dsPtr = analysisManager_.getDataStore();

        if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS) )
        {
          Xyce::dout() << "Loading initial condition data from time: fastTimes_["
                       << i << "] = " << fastTimes_[i] << std::endl;
        }
        if ( Xyce::isActive(Xyce::Diag::MPDE_PRINT_VECTORS) )
        {
          Xyce::dout() << "mpdeICVectorPtr_->block(" << i << ") = dsPtr->fastTimeSolutionVec[" << indicesUsed_[i] << "] = " << std::endl;
          (dsPtr->fastTimeSolutionVec[indicesUsed_[i]])->print(Xyce::dout());
          Xyce::dout() << "mpdeICStateVectorPtr_->block(" << i << ") = dsPtr->fastTimeStateVec[" << indicesUsed_[i] << "] = " << std::endl;
          (dsPtr->fastTimeStateVec[indicesUsed_[i]])->print(Xyce::dout());
          Xyce::dout() << "mpdeICQVectorPtr_->block(" << i << ") = dsPtr->fastTimeQVec[" << indicesUsed_[i] << "] = " << std::endl;
          (dsPtr->fastTimeQVec[indicesUsed_[i]])->print(Xyce::dout());
          Xyce::dout() << "mpdeICStoreVectorPtr_->block(" << i << ") = dsPtr->fastTimeStoreVec[" << indicesUsed_[i] << "] = " << std::endl;
          (dsPtr->fastTimeStoreVec[indicesUsed_[i]])->print(Xyce::dout());
        }

        if (warpMPDE_)
        {
          // Store the "good" initial condition in the mpdeIC vectors. The fastTime vectors
          // is a temporary storage and doesn't have to be updated.
          // A good initial condition is solutions between sol_A and sol_B where sol_A = sol_B =0
          // and has the same slope (e.g. increasing)
          // Also, print the "mpde_ic" file here using values in mpdeIC vectors.

          mpdeICVectorPtr_->block(i) = *(dsPtr->fastTimeSolutionVec[indicesUsed_[i+goodStartPoint]]);
          mpdeICStateVectorPtr_->block(i) = *(dsPtr->fastTimeStateVec[indicesUsed_[i+goodStartPoint]]);
          mpdeICQVectorPtr_->block(i) = *(dsPtr->fastTimeQVec[indicesUsed_[i+goodStartPoint]]);
          mpdeICStoreVectorPtr_->block(i) = *(dsPtr->fastTimeStoreVec[indicesUsed_[i+goodStartPoint]]);

          std::cout << "Loading initial condition data from updated time: fastTimes_["
                    << i << "] = " << fastTimes_[i] << std::endl;
          std::cout << "mpdeICVectorPtr_->block(" << i << ") = dsPtr->fastTimeSolutionVec[" << indicesUsed_[i+goodStartPoint] << "] = " << std::endl;
          ////(dsPtr->fastTimeSolutionVec[indicesUsed_[i]])->print(Xyce::dout());
          (mpdeICVectorPtr_->block(i)).print(Xyce::dout());
          std::cout << "mpdeICStateVectorPtr_->block(" << i << ") = dsPtr->fastTimeStateVec[" << indicesUsed_[i+goodStartPoint] << "] = " << std::endl;
          ////(dsPtr->fastTimeStateVec[indicesUsed_[i]])->print(Xyce::dout());
          std::cout << "mpdeICQVectorPtr_->block(" << i << ") = dsPtr->fastTimeQVec[" << indicesUsed_[i+goodStartPoint] << "] = " << std::endl;
          ////(dsPtr->fastTimeQVec[indicesUsed_[i]])->print(Xyce::dout());
          std::cout << "mpdeICStoreVectorPtr_->block(" << i << ") = dsPtr->fastTimeStoreVec[" << indicesUsed_[i+goodStartPoint] << "] = " << std::endl;
          ////(dsPtr->fastTimeStoreVec[indicesUsed_[i]])->print(Xyce::dout());
        }
        else
        {
          mpdeICVectorPtr_->block(i) = *(dsPtr->fastTimeSolutionVec[indicesUsed_[i]]);
          mpdeICStateVectorPtr_->block(i) = *(dsPtr->fastTimeStateVec[indicesUsed_[i]]);
          mpdeICQVectorPtr_->block(i) = *(dsPtr->fastTimeQVec[indicesUsed_[i]]);
          mpdeICStoreVectorPtr_->block(i) = *(dsPtr->fastTimeStoreVec[indicesUsed_[i]]);

          std::cout << "Loading initial condition data from time: fastTimes_["
                    << i << "] = " << fastTimes_[i] << std::endl;
          std::cout << "mpdeICVectorPtr_->block(" << i << ") = dsPtr->fastTimeSolutionVec[" << indicesUsed_[i] << "] = " << std::endl;
          ////(dsPtr->fastTimeSolutionVec[indicesUsed_[i]])->print();
          std::cout << "mpdeICStateVectorPtr_->block(" << i << ") = dsPtr->fastTimeStateVec[" << indicesUsed_[i] << "] = " << std::endl;
          ////(dsPtr->fastTimeStateVec[indicesUsed_[i]])->print();
          std::cout << "mpdeICQVectorPtr_->block(" << i << ") = dsPtr->fastTimeQVec[" << indicesUsed_[i] << "] = " << std::endl;
          ////(dsPtr->fastTimeQVec[indicesUsed_[i]])->print();
          std::cout << "mpdeICStoreVectorPtr_->block(" << i << ") = dsPtr->fastTimeStoreVec[" << indicesUsed_[i] << "] = " << std::endl;
          ////(dsPtr->fastTimeStoreVec[indicesUsed_[i]])->print();
        }

      }


      if (warpMPDE_)
      {
        std::cout << "mpdeIC vectors NEW blockCount = " << mpdeICVectorPtr_->blockCount() << std::endl;

        // We've basically skipped some time for better WaMPDE initial condition.  So
        // have MPDE think it's starting a bit later.

        // get the time integrator params, and save a copy.
        Xyce::TimeIntg::TIAParams tia_params;

        // set DAE initial time = previous initial time + skipped time for better MPDE IC
        tiaMPDEParams_.initialTime = tiaMPDEParams_.initialTime + shiftTimeBy;


        std::cout << "Beginning WAMPDE_IC_TRANSIENT: MPDE_initialTime = " << tiaMPDEParams_.initialTime
                  << " MPDE_finalTime = " << tiaMPDEParams_.finalTime << std::endl;

        std::cout << "Beginning WAMPDE_IC_TRANSIENT: initialTime = " << tia_params.initialTime
                  << " finalTime = " << tia_params.finalTime << std::endl;
      }
    }
    break;

    case MPDE_IC_TRAN_MAP:
    {
      // at this point fastTimes_ is either setup
      // with fixed spaced points or a variable mesh.

      if( !tranRunForSize_ )
      {
        // now try the real thing
        Xyce::TimeIntg::TIAParams tia_params = tiaMPDEParams_;

        tia_params.initialTime = tiaMPDEParams_.initialTime;

        if( initialCondition_ == MPDE_IC_TWO_PERIOD )
        {
          // in this case we need to integrate forward 2 periods
          tia_params.finalTime = tiaMPDEParams_.initialTime + 2.0 * period_;
        }
        else
        {
          tia_params.finalTime = tiaMPDEParams_.initialTime + period_;
        }

        analysisManager_.getStepErrorControl().pauseTime = tia_params.finalTime;
        returnValue = runTransientIC(tia_params, linear_system, nonlinear_manager);
        // normally after runTransientIC() we would filter the resulting
        // points down.  We only get this this part of the code
        // if one called this initial condition type, but wanted fixed
        // spaced points on the fast time scale.  So, we don't call
        // filterFastTimePoints() here.
      }

      // use the last point for the endICvector
      Xyce::TimeIntg::DataStore *dsPtr = analysisManager_.getDataStore();

      const Xyce::Linear::Vector &endIcSolVecPtr_ = *dsPtr->fastTimeSolutionVec[indicesUsed_[n2-1]];
      const Xyce::Linear::Vector &endIcStateVecPtr_ = *dsPtr->fastTimeStateVec[indicesUsed_[n2-1]];
      const Xyce::Linear::Vector &endIcQVecPtr_ = *dsPtr->fastTimeQVec[indicesUsed_[n2-1]];
      const Xyce::Linear::Vector &endIcStoreVecPtr_ = *dsPtr->fastTimeStoreVec[indicesUsed_[n2-1]];

      // now we linearlly interpolate between the data in
      // dcOpSolVecPtr_    and  endIcSolVecPtr_
      // dcOpStateVecPtr_  and  endIcStateVecPtr_
      // dcOpQVecPtr_      and  endIcQVecPtr_
      //
      // and, if enabled:
      // dcOpStoreVecPtr_  and  endIcStoreVecPtr_
      //
      // for each point in fastTimes_ and use that value
      // as an initial condition to go from t=fastTimes_[i] to period_
      //

      mpdeICVectorPtr_->block(0) = endIcSolVecPtr_;
      mpdeICStateVectorPtr_->block(0) = endIcStateVecPtr_;
      mpdeICQVectorPtr_->block(0) = endIcQVecPtr_;
      mpdeICStoreVectorPtr_->block(0) = endIcStoreVecPtr_;

      //Setting up specialized Loader
      N_MPDE_SawtoothLoader stLoader( mpdeState_, analysisManager_, appLoader_ );
      stLoader.setTimeShift( fastTimes_[0] );

      // get the time integrator params, and save a copy.
      Xyce::TimeIntg::TIAParams tia_params;

      Xyce::Linear::Vector * interpIcSolVecPtr = endIcSolVecPtr_.cloneVector();
      Xyce::Linear::Vector * interpIcStateVecPtr = endIcStateVecPtr_.cloneVector();
      Xyce::Linear::Vector * interpIcQVecPtr = endIcQVecPtr_.cloneVector();
      Xyce::Linear::Vector * interpIcStoreVecPtr = endIcStoreVecPtr_.cloneVector();

      for (int i=1 ; i<n2 ; ++i)
      {
        // set DAE initial time = 0.0
        tia_params.initialTime = tiaMPDEParams_.initialTime + fastTimes_[n2-i];
        tia_params.finalTime = tiaMPDEParams_.initialTime + period_;
        tia_params.restartTimeStepScale = 1.0;  // try to take big steps.

        //tia_params.finalTime = fastTimes_[i];
        analysisManager_.getStepErrorControl().pauseTime = tia_params.finalTime;

        // We need to do an initialization every time...
        analysisManager_.setBeginningIntegrationFlag(true);
        analysisManager_.getStepErrorControl().resetAll(tia_params);
        analysisManager_.getWorkingIntegrationMethod().initialize(tia_params);

        // here's were we interpolate to find the initial condition
        double fraction = fastTimes_[n2-i] / period_;

        if (DEBUG_MPDE)
          Xyce::dout() << "Beginning MPDE_IC_TRAN_MAP IC: i = " << i
                       << " fraction = " << fraction
                       << " initialTime = " << tia_params.initialTime
                       << " finalTime = " << tia_params.finalTime << std::endl;

        // Set t = t + (i-1.0)*h2 in MPDE source in bhat
        stLoader.setTimeShift( -fastTimes_[n2-i] );

        interpIcSolVecPtr->update( (1.0-fraction), *dcOpSolVecPtr_, fraction, endIcSolVecPtr_, 0.0 );
        interpIcStateVecPtr->update( (1.0-fraction), *dcOpStateVecPtr_, fraction, endIcStateVecPtr_, 0.0 );
        interpIcQVecPtr->update( (1.0-fraction), *dcOpQVecPtr_, fraction, endIcQVecPtr_, 0.0 );
        interpIcStoreVecPtr->update( (1.0-fraction), *dcOpStoreVecPtr_, fraction, endIcStoreVecPtr_, 0.0 );

        // set DAE initial condition to mpdeICVectorPtr_[i-1]
        *(analysisManager_.getDataStore()->nextSolutionPtr) = *interpIcSolVecPtr;
        *(analysisManager_.getDataStore()->nextStatePtr) = *interpIcStateVecPtr;
        *(analysisManager_.getDataStore()->daeQVectorPtr) = *interpIcQVecPtr;
        *(analysisManager_.getDataStore()->nextStorePtr) = *interpIcStoreVecPtr;

        // throw std::runtime_error("This has not been tested using Transient locally");

        {
          // This prepares the outputters for this analysis mode
          Xyce::IO::ActiveOutput active(analysisManager_.getOutputManagerAdapter().getOutputManager());

          active.setStepSweepVector(analysisManager_.getOutputManagerAdapter().getStepSweepVector());
          active.setDCSweepVector(analysisManager_.getOutputManagerAdapter().getDCSweepVector());
          active.add(pdsManager_.getPDSComm()->comm(), Xyce::Analysis::ANP_MODE_TRANSIENT);

          if (VERBOSE_TIME)
            tia_params.printParams(Xyce::lout(), nonlinearAnalysisMode(Xyce::Analysis::ANP_MODE_TRANSIENT));

          analysisManager_.setAnalysisMode(Xyce::Analysis::ANP_MODE_TRANSIENT);

          //Xyce::Analysis::Transient transient(analysisManager_, linear_system, nonlinear_manager, stLoader, topology_, initialConditionsManager_, restartManager_, &outputAdapter_, 0, 0);
          Xyce::Analysis::Transient transient(analysisManager_, &linear_system, nonlinear_manager, stLoader, topology_, initialConditionsManager_, restartManager_, &outputAdapter_, 0, 0);
          analysisManager_.pushActiveAnalysis(&transient);
          transient.setTIAParams(tia_params);
          transient.setAnalysisParams(Xyce::Util::OptionBlock());
          transient.setNOOP(true);
          analysisManager_.getStepErrorControl().resetAll(tia_params);

          returnValue = transient.run();

          transientNeedsToLoadInitialConditionsAndInitializeProblem_ = false;

          analysisManager_.popActiveAnalysis();
        }

        // put the dsPtr->currentSolutionPtr into mpdeICVectorPtr_[i]
        mpdeICVectorPtr_->block(i) = *analysisManager_.getDataStore()->currSolutionPtr;
        mpdeICStateVectorPtr_->block(i) = *analysisManager_.getDataStore()->currStatePtr;
        mpdeICQVectorPtr_->block(i) = *analysisManager_.getDataStore()->daeQVectorPtr;
        mpdeICStoreVectorPtr_->block(i) = *analysisManager_.getDataStore()->currStorePtr;
        if (DEBUG_MPDE)
          Xyce::dout() << "End MPDE_IC_TRAN_MAP IC: i = " << i <<std::endl;
      }

      delete interpIcSolVecPtr;
      delete interpIcStateVecPtr;
      delete interpIcQVecPtr; 
      delete interpIcStoreVecPtr;

      tiaMPDEParams_.initialTime = tia_params.finalTime;
    }
    break;

    case MPDE_IC_TWO_PERIOD:
    {

      // The runTransientIC function would have done two periods of
      // the fast time source.  Now we need to look at the points in
      // the second period and match them up with the ones in the first
      // period and then do the inerpolation.

      // fastTimes_[] will already have been offset to run between 0
      // and one period_.  We can use indicesUsed_[] and the timeStep[]
      // array still in the dataStore to find solutions in the first
      // cycle that lie close to those in the fastTimes_[] cycle.

      Xyce::TimeIntg::DataStore *dsPtr = analysisManager_.getDataStore();
      int numPoints = dsPtr->timeSteps.size();
      double startUpOffset = tiaMPDEParams_.initialTime;
      int lastJused = 0;
      for( int i=0; i < size_; i++ )
      {
        for( int j=lastJused; j < numPoints; j++ )
        {
          if( (fastTimes_[i] >= (dsPtr->timeSteps[j] - startUpOffset)) &&
              (fastTimes_[i] <  (dsPtr->timeSteps[j+1] - startUpOffset)) )
          {
            // found time points in the first period that bracket this one
            if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS) )
            {
              Xyce::dout() << "For fastTime_[" << i << "] = " << fastTimes_[i] <<
                " Found bracketing times : " << (dsPtr->timeSteps[j] - startUpOffset) << ", "
                           << (dsPtr->timeSteps[j+1] - startUpOffset)
                           << " at j = " << j;
            }

            // assume we're closest to the first point, but check if we're
            // really closest to the second one
            int interpolationPoint = j;
            if( fabs( fastTimes_[i] - (dsPtr->timeSteps[j] - startUpOffset)) >
                fabs( fastTimes_[i] - (dsPtr->timeSteps[j+1] - startUpOffset)) )
            {
              interpolationPoint = j+1;
            }
            if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS) )
              Xyce::dout() << " closest to point " << interpolationPoint;

            // now do the interpolation
            Xyce::Linear::Vector &firstPeriodSolVecPtr = *dsPtr->fastTimeSolutionVec[ interpolationPoint ];
            Xyce::Linear::Vector &firstPeriodStateVecPtr = *dsPtr->fastTimeStateVec[ interpolationPoint ];
            Xyce::Linear::Vector &firstPeriodQVecPtr = *dsPtr->fastTimeQVec[ interpolationPoint ];
            Xyce::Linear::Vector &firstPeriodStoreVecPtr = *dsPtr->fastTimeStoreVec[ interpolationPoint ];

            Xyce::Linear::Vector &secondPeriodSolVecPtr = *dsPtr->fastTimeSolutionVec[indicesUsed_[i]];
            Xyce::Linear::Vector &secondPeriodStateVecPtr = *dsPtr->fastTimeStateVec[indicesUsed_[i]];
            Xyce::Linear::Vector &secondPeriodQVecPtr = *dsPtr->fastTimeQVec[indicesUsed_[i]];
            Xyce::Linear::Vector &secondPeriodStoreVecPtr = *dsPtr->fastTimeStoreVec[indicesUsed_[i]];

            Xyce::Linear::Vector * interpIcSolVecPtr = secondPeriodSolVecPtr.cloneVector();
            Xyce::Linear::Vector * interpIcStateVecPtr = secondPeriodStateVecPtr.cloneVector();
            Xyce::Linear::Vector * interpIcQVecPtr = secondPeriodQVecPtr.cloneVector();
            Xyce::Linear::Vector * interpIcStoreVecPtr = secondPeriodStoreVecPtr.cloneVector();

            double fraction = fastTimes_[i] / period_;
            if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS) )
              Xyce::dout() << " fraction = " << fraction << std::endl;

            interpIcSolVecPtr->update( fraction, firstPeriodSolVecPtr, (1.0-fraction), secondPeriodSolVecPtr, 0.0 );
            interpIcStateVecPtr->update( fraction, firstPeriodStateVecPtr, (1.0-fraction), secondPeriodStateVecPtr, 0.0 );
            interpIcQVecPtr->update( fraction, firstPeriodQVecPtr, (1.0-fraction), secondPeriodQVecPtr, 0.0 );
            interpIcStoreVecPtr->update( fraction, firstPeriodStoreVecPtr, (1.0-fraction), secondPeriodStoreVecPtr, 0.0 );

            mpdeICVectorPtr_->block(i) = *interpIcSolVecPtr;
            mpdeICStateVectorPtr_->block(i) = *interpIcStateVecPtr;
            mpdeICQVectorPtr_->block(i) = *interpIcQVecPtr;
            mpdeICStoreVectorPtr_->block(i) = *interpIcStoreVecPtr;

            lastJused = interpolationPoint;

            delete interpIcSolVecPtr;
            delete interpIcStateVecPtr;
            delete interpIcQVecPtr; 
            delete interpIcStoreVecPtr;
            
            break;  // break out of for(j..) loop.  We don't need to cycle more
          }
        }
      }

      // We've basically skipped a period with this initial condition.  So
      // have MPDE think it's starting a bit later.
      tiaMPDEParams_.initialTime += period_;

      // stil under development RLS
      //double error = checkPeriodicity_();
      //Xyce::dout() << "N_MPDE_Manager::runInitialConditon.  Based on two period anaysis, non-periodic error in this problem is " << error << std::endl;
    }
    break;
    case MPDE_IC_TRANSIENT:
    {
      if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS) )
      {
        for( unsigned int i=0; i<indicesUsed_.size(); i++ )
          Xyce::dout() << "indicesUsed_[ " << i << " ] = " << indicesUsed_[i] << std::endl;
      }

      // we did an initial transient run so pick up ic data from that
      for (int i=0 ; i<n2 ; ++i)
      {
        Xyce::TimeIntg::DataStore *dsPtr = analysisManager_.getDataStore();
        if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS) )
        {
          Xyce::dout() << "Loading initial condition data from time: fastTimes_["
                       << i << "] = " << fastTimes_[i] << std::endl;
        }
        if ( Xyce::isActive(Xyce::Diag::MPDE_PRINT_VECTORS) )
        {
          Xyce::dout() << "mpdeICVectorPtr_->block(" << i << ") = dsPtr->fastTimeSolutionVec[" << indicesUsed_[i] << "] = " << std::endl;
          (dsPtr->fastTimeSolutionVec[indicesUsed_[i]])->print(Xyce::dout());
          Xyce::dout() << "mpdeICStateVectorPtr_->block(" << i << ") = dsPtr->fastTimeStateVec[" << indicesUsed_[i] << "] = " << std::endl;
          (dsPtr->fastTimeStateVec[indicesUsed_[i]])->print(Xyce::dout());
          Xyce::dout() << "mpdeICQVectorPtr_->block(" << i << ") = dsPtr->fastTimeQVec[" << indicesUsed_[i] << "] = " << std::endl;
          (dsPtr->fastTimeQVec[indicesUsed_[i]])->print(Xyce::dout());
          Xyce::dout() << "mpdeICStoreVectorPtr_->block(" << i << ") = dsPtr->fastTimeStoreVec[" << indicesUsed_[i] << "] = " << std::endl;
          (dsPtr->fastTimeStoreVec[indicesUsed_[i]])->print(Xyce::dout());
        }

        mpdeICVectorPtr_->block(i) = *(dsPtr->fastTimeSolutionVec[indicesUsed_[i]]);
        mpdeICStateVectorPtr_->block(i) = *(dsPtr->fastTimeStateVec[indicesUsed_[i]]);
        mpdeICQVectorPtr_->block(i) = *(dsPtr->fastTimeQVec[indicesUsed_[i]]);
        mpdeICStoreVectorPtr_->block(i) = *(dsPtr->fastTimeStoreVec[indicesUsed_[i]]);
      }

    }
    break;

    default:
      Xyce::Report::UserError() << "Invalid IC option specified";
      break;
  }

  if (warpMPDE_)
  {
    //std::cout << "Setting omega = 1 into MPDE solution vector." << std::endl;
    // tscoffe/tmei 08/10/05:  put omega = 1 into MPDE solution vector
    // This is true for MPDE and WaMPDE because in WaMPDE we scale omega by T2.
    // tscoffe 12/14/06:  put phi(0) = omega(0) into MPDE solution vector.
    int omegaGID = warpMPDEPhasePtr_->getOmegaGID();
    int phiGID = warpMPDEPhasePtr_->getPhiGID();
    int omegaLID = mpdeICVectorPtr_->pmap()->globalToLocalIndex(omegaGID);
    int phiLID = mpdeICVectorPtr_->pmap()->globalToLocalIndex(phiGID);
    if (omegaLID >= 0)
      (*mpdeICVectorPtr_)[omegaLID] = 1.0;
    if (phiLID >= 0)
      (*mpdeICVectorPtr_)[phiLID] = 0.0;
  }

  if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PRINT_VECTORS) )
  {
    Xyce::dout() << "MPDE Initial Condition Solution!\n";
    mpdeICVectorPtr_->print(std::cout);
    Xyce::dout() << "MPDE Initial Condition State Vector!\n";
    mpdeICStateVectorPtr_->print(std::cout);
    Xyce::dout() << "MPDE Initial Condition Store Vector!\n";
    mpdeICStoreVectorPtr_->print(std::cout);
  }

  if (DEBUG_MPDE)
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
    if (icExitFlag_)
      Xyce::Report::DevelFatal() << "Exiting on icExitFlag_";
  }

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::runDCop()
// Purpose       : Runs the dc op problem
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, 1437, Electrical and Micro Modeling
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool
N_MPDE_Manager::runDCOP(
  Xyce::Linear::System &    linear_system,
  Xyce::Nonlinear::Manager &    nonlinear_manager)
{
  bool success = true;

  // Perform DCOP
  if (!NOOP_)
  {
    // as long as the user didn't request "noop"
    // try and do the operating point calculation
    {
      //Xyce::Analysis::DCSweep dc_sweep(analysisManager_, linear_system, nonlinear_manager, appLoader_, topology_, initialConditionsManager_, 0);
      Xyce::Analysis::DCSweep dc_sweep(analysisManager_, &linear_system, nonlinear_manager, appLoader_, topology_, initialConditionsManager_, 0);
      analysisManager_.pushActiveAnalysis(&dc_sweep);
      dc_sweep.run();

      analysisManager_.popActiveAnalysis();
    }
  }
  else
    success = false;

  // store the dc op results in case we need them later
  dcOpSolVecPtr_ = analysisManager_.getDataStore()->currSolutionPtr->cloneCopyVector();
  dcOpStateVecPtr_ = analysisManager_.getDataStore()->currStatePtr->cloneCopyVector();
  dcOpQVecPtr_ = analysisManager_.getDataStore()->daeQVectorPtr->cloneCopyVector();
  dcOpStoreVecPtr_ = analysisManager_.getDataStore()->currStorePtr->cloneCopyVector();

  return success;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::runStartupPeriods()
// Purpose       : Runs normal transient problem through the requested
//                 number of startup periods
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, 1437, Electrical and Micro Modeling
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool
N_MPDE_Manager::runStartupPeriods(
  const Xyce::TimeIntg::TIAParams &     tia_params,
  Xyce::Linear::System &                linear_system,
  Xyce::Nonlinear::Manager &            nonlinear_manager)
{
  bool returnValue = true;

  // startUpPeriodsFlag_ = true;

  Xyce::dout() << "Advancing time through " << startUpPeriods_ << " startup periods"
               << " initialTime = " << tia_params.initialTime
               << " initialTimeStep = " << tia_params.initialTimeStep
               << " finalTime = " << tia_params.finalTime << std::endl;

  {
    Xyce::IO::ActiveOutput x(analysisManager_.getOutputManagerAdapter().getOutputManager());
    x.add(Xyce::IO::PrintType::MPDE_STARTUP, Xyce::Analysis::ANP_MODE_MPDE);

    if (VERBOSE_TIME)
      tia_params.printParams(Xyce::lout(), nonlinearAnalysisMode(Xyce::Analysis::ANP_MODE_TRANSIENT));

    analysisManager_.setAnalysisMode(Xyce::Analysis::ANP_MODE_TRANSIENT);

    //Xyce::Analysis::Transient transient(analysisManager_, linear_system, nonlinear_manager, appLoader_, topology_, initialConditionsManager_, restartManager_, 0, 0, this);
    Xyce::Analysis::Transient transient(analysisManager_, &linear_system, nonlinear_manager, appLoader_, topology_, initialConditionsManager_, restartManager_, 0, 0, this);
    analysisManager_.pushActiveAnalysis(&transient);
    transient.setTIAParams(tia_params);
    transient.setAnalysisParams(Xyce::Util::OptionBlock());
    transient.setNOOP(!warpMPDEICFlag_);
    analysisManager_.getStepErrorControl().resetAll(tia_params);

    returnValue = transient.run();

    transientNeedsToLoadInitialConditionsAndInitializeProblem_ = false;

    analysisManager_.popActiveAnalysis();

    analysisManager_.getOutputManagerAdapter().getOutputManager().finishOutput();
  }

  // put the dsPtr->currentSolutionPtr into dcOpSol and State Vec so that it
  // is used as our initial condition for the pending fast time scale runs
  dcOpSolVecPtr_ = analysisManager_.getDataStore()->currSolutionPtr->cloneCopyVector();
  dcOpStateVecPtr_ =  analysisManager_.getDataStore()->currStatePtr->cloneCopyVector();
  dcOpQVecPtr_ = analysisManager_.getDataStore()->daeQVectorPtr->cloneCopyVector();
  dcOpStoreVecPtr_ = analysisManager_.getDataStore()->currStorePtr->cloneCopyVector();

  // startUpPeriodsFlag_ = false;

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::runTransientIC
// Purpose       : Conducts a regular transient run for MPDE initial conditions
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, 1437, Electrical and Micro Modeling
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool
N_MPDE_Manager::runTransientIC(
  const Xyce::TimeIntg::TIAParams &     tia_params,
  Xyce::Linear::System &                linear_system,
  Xyce::Nonlinear::Manager &            nonlinear_manager)
{
  bool returnValue = false;

  if (DEBUG_MPDE)
    Xyce::dout() << "N_MPDE_Manager::runTransientIC() Doing initial transient run." << std::endl;

  // flag for Xyce::TimeIntg::ControlgAlgorithm that we're calculating an IC
  // this prevents extra DC op data from being printed.
  // mpdeICFlag_ = true;

  if (!warpMPDE_)
  {
    if (!dcOpSolVecPtr_)
    {
      // haven't done a dc op yet, so do it now
      if (!NOOP_)
      {
        // as long as the user didn't request "noop"
        // try and do the operating point calculation
        {
          Xyce::Analysis::DCSweep dc_sweep(analysisManager_, &linear_system, nonlinear_manager, appLoader_, topology_, initialConditionsManager_, 0);
          analysisManager_.pushActiveAnalysis(&dc_sweep);
          dc_sweep.run();

          // // Add in simulation times
          // accumulateStatistics_(dc_sweep);

          // print out analysis info
          Xyce::lout() << " ***** Harmonic Balance Computation Summary *****" << std::endl;
          dc_sweep.printLoopInfo( 0, 0 );

          analysisManager_.popActiveAnalysis();
        }
      }
    }
  }

  if (dcOpSolVecPtr_)
  {
    *(analysisManager_.getDataStore()->nextSolutionPtr) = *dcOpSolVecPtr_;
    *(analysisManager_.getDataStore()->nextStatePtr) = *dcOpStateVecPtr_;
    *(analysisManager_.getDataStore()->daeQVectorPtr) = *dcOpQVecPtr_;
    *(analysisManager_.getDataStore()->nextStorePtr) = *dcOpStoreVecPtr_;
  }

  Xyce::IO::ActiveOutput active(analysisManager_.getOutputManagerAdapter().getOutputManager());

  if (saveIcData_)
  {
    active.add(Xyce::IO::PrintType::MPDE_IC, Xyce::Analysis::ANP_MODE_TRANSIENT);
  }

  if (VERBOSE_TIME)
    tia_params.printParams(Xyce::lout(), nonlinearAnalysisMode(Xyce::Analysis::ANP_MODE_TRANSIENT));

  analysisManager_.setAnalysisMode(Xyce::Analysis::ANP_MODE_TRANSIENT);

  Xyce::Analysis::Transient transient(analysisManager_, &linear_system, nonlinear_manager, appLoader_, topology_, initialConditionsManager_, restartManager_, &outputAdapter_, 0, this);
  analysisManager_.pushActiveAnalysis(&transient);
  transient.setTIAParams(tia_params);
  transient.setAnalysisParams(Xyce::Util::OptionBlock());
  transient.setNOOP(true);
  transient.setSaveTimeSteps(true);
  analysisManager_.getStepErrorControl().resetAll(tia_params);

  returnValue = transient.run();

  transientNeedsToLoadInitialConditionsAndInitializeProblem_ = false;

  analysisManager_.popActiveAnalysis();

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::filterFastTimePoints()
// Purpose       : Tries to filter the fast time points from a transient run
//                 so that points are not too close together
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, 1437, Electrical and Micro Modeling
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool
N_MPDE_Manager::filterFastTimePoints(
  const Xyce::TimeIntg::TIAParams &     tia_params)
{
  Xyce::TimeIntg::DataStore *dsPtr = analysisManager_.getDataStore();
  int numPoints = dsPtr->timeSteps.size();
  int startIndex = 0;
  double extraPeriodShift = 0.0;

  // Note, if we integrated more then one period in the runTransientIC()
  // function, the we'll try to pull in points only from the second period.
  if( initialCondition_ == MPDE_IC_TWO_PERIOD )
  {
    extraPeriodShift = period_;
    for(int i=0; i<numPoints; i++ )
    {
      if(dsPtr->timeSteps[i] > (tia_params.initialTime + period_))
      {
        startIndex = i;
        break;
      }
    }
  }

  if (DEBUG_MPDE)
    Xyce::dout() << "Initial transient run produced " << numPoints << " points." << std::endl;

  std::list<int> goodIndicies;
  int numGoodPoints = 0;
  int breakpoints = 0;          // need to keep track of how many breakpoints there are

  if( !maxCalcSizeGiven_  || (numPoints < maxCalcSize_ ) )
  {
    // keep all of the points from the transient run

    for( int i=startIndex; i < numPoints; i++ )
    {
      goodIndicies.push_back( i );
    }
    numGoodPoints = goodIndicies.size();

    if (DEBUG_MPDE)
      Xyce::dout() << " keeping all points for MPDE calculations." << std::endl;
  }
  else
  {
    // locate the set of points and indicies that are "good", i.e. those that
    // are not too close together
    const double relTimeTol = 0.005;
    const double absTimeDiffTol = 5.0e-9;

    // always keep first point
    goodIndicies.push_back(startIndex);
    int lastGoodIndex = startIndex;

    for( int i=startIndex + 1; i < numPoints; i++ )
    {
      // count up breakpoints
      if( dsPtr->timeStepsBreakpointFlag[i] == true )
      {
        breakpoints++;
      }

      if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS) )
      {
        Xyce::dout() << "\t\t timeStep[ " << i << " ] = " << dsPtr->timeSteps[i];
        if( dsPtr->timeStepsBreakpointFlag[i] == true )
        {
          Xyce::dout() << "  Breakpoint";
        }
        Xyce::dout() << std::endl;
      }

      double delta_t = (dsPtr->timeSteps[i] - dsPtr->timeSteps[lastGoodIndex]);
      if( ((delta_t > relTimeTol * period_) && (delta_t > absTimeDiffTol)) ||  dsPtr->timeStepsBreakpointFlag[i] == true )
      {
        // found a good point so save the index
        goodIndicies.push_back( i );
        lastGoodIndex = i;
      }
    }
    numGoodPoints = goodIndicies.size();

    if (DEBUG_MPDE)
      Xyce::dout() << " of which " << numGoodPoints << " fit tolerance requirements of abstol = "
                   << absTimeDiffTol << " and reltol*period = " << relTimeTol * period_ <<std::endl;
  }

  if( (maxCalcSizeGiven_)  && (maxCalcSize_ < numGoodPoints ) )
  {
    // can't just copy all of the fast times,
    // filter down to maxCalcSizeGiven
    //vector<double> tmpTimePoints = dsPtr->timeSteps;
    fastTimes_.resize( maxCalcSize_ + 1);
    indicesUsed_.resize( maxCalcSize_ );
    // always grab the first and last point
    fastTimes_[0] = dsPtr->timeSteps[startIndex] - tia_params.initialTime - extraPeriodShift;
    indicesUsed_[0] = startIndex;
    fastTimes_[maxCalcSize_] = period_;

    // this is a little confusing but here is what we're trying to do
    // we have a list of indicies that are considered good points.
    // we need to sample this list a regular intervals to get it
    // into the max size requested by the users.  We also have
    // points that are breakpoints.  We must keep all the breakpoints.

    double sampleRate = (1.0*(numGoodPoints-1-breakpoints)) / (1.0*(maxCalcSize_-1));

    std::list<int>::iterator currentIndex = goodIndicies.begin();
    std::list<int>::iterator endIndex = goodIndicies.end();

    int indexCounter=0, k=0;
    while( currentIndex != endIndex )
    {
      int targetIndex = static_cast<int>( k * sampleRate + 0.5 );
      if( ((indexCounter == targetIndex) && (k < maxCalcSize_)) ||
        (dsPtr->timeStepsBreakpointFlag[*currentIndex] == true) )
      {
        indicesUsed_[k] = *currentIndex;
        if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS) )
        {
          Xyce::dout() << " indicesUsed_[" << k << "] = " << indicesUsed_[k] << std::endl;
        }

        fastTimes_[k] = dsPtr->timeSteps[*currentIndex] - tia_params.initialTime - extraPeriodShift;
        k++;
      }
      currentIndex++;
      indexCounter++;
    }
  }
  else
  {
    // number of good points is less than the requested number
    // of points so just copy them all over
    //
    fastTimes_.resize( numGoodPoints + 1 );

    indicesUsed_.resize( numGoodPoints );
    fastTimes_[0] = 0.0;
    indicesUsed_[0] = startIndex;

    fastTimes_[numGoodPoints] = period_;

    std::list<int>::iterator currentIndex = goodIndicies.begin();
    std::list<int>::iterator endIndex = goodIndicies.end();

    int k=1;
    currentIndex++;
    while( currentIndex != endIndex )
    {
      indicesUsed_[k] = *currentIndex;
      if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS) )
      {
        Xyce::dout() << " indicesUsed_[" << k << "] = " << indicesUsed_[k] << std::endl;
      }

      fastTimes_[k] = dsPtr->timeSteps[*currentIndex] - tia_params.initialTime - extraPeriodShift;
      k++;
      currentIndex++;
    }
  }
  size_ = fastTimes_.size() - 1;

  Xyce::dout() << "MPDE: " << size_ << " fast time points added to the problem." << std::endl
    << "fast time point range: " << fastTimes_[0] << ", " << fastTimes_[size_] << std::endl;
  if( Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS) )
  {
    Xyce::dout() << "MPDE: new fast times are:" << std::endl;
    int i=0;
    for( i = 0; i < size_; ++i )
      Xyce::dout() << "fastTimes_["<<i<<"] = " << fastTimes_[i]
        << " forward difference is " << (fastTimes_[i+1] - fastTimes_[i]) << std::endl;
    Xyce::dout() << "period = " << period_
      << " fastTimes_[numGoodPoints] = " << fastTimes_[numGoodPoints]
      << " extraPeriodShift = " << extraPeriodShift << std::endl;
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::setupMPDEProblem_
// Purpose       : Configures solvers, etc. for MPDE run
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::setupMPDEProblem_()
{
  if (DEBUG_MPDE)
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  N_MPDE_Manager::setupMPDEProblem" << std::endl;

  Xyce::lout() << " ***** Setting up full MPDE problem....\n" << std::endl;

  // Destroy Solvers, etc. from IC phase
  analysisManager_.resetSolverSystem();

  // Finish setup of MPDE Builder
  mpdeBuilderPtr_->generateGraphs( *(pdsManager_.getMatrixGraph( Xyce::Parallel::JACOBIAN )) );

  // Setup MPDE Loader
  mpdeLoaderPtr_ = new N_MPDE_Loader(mpdeState_, deviceManager_, appLoader_, *mpdeDiscPtr_, *this, warpMPDEPhasePtr_);
  mpdeLoaderPtr_->setFastTimes( fastTimes_ );
  mpdeLoaderPtr_->setPeriodFlags( nonPeriodicFlags );

  // Include these for the loader to use with devices
  mpdeLoaderPtr_->registerAppNextVec ( rcp(appBuilder_.createVector()) );
  mpdeLoaderPtr_->registerAppCurrVec ( rcp(appBuilder_.createVector()) );
  mpdeLoaderPtr_->registerAppLastVec ( rcp(appBuilder_.createVector()) );

  mpdeLoaderPtr_->registerAppNextStaVec ( rcp(appBuilder_.createStateVector()) );
  mpdeLoaderPtr_->registerAppCurrStaVec ( rcp(appBuilder_.createStateVector()) );
  mpdeLoaderPtr_->registerAppLastStaVec ( rcp(appBuilder_.createStateVector()) );
  mpdeLoaderPtr_->registerOmegadQdt2( rcp_dynamic_cast<Xyce::Linear::BlockVector>(rcp(mpdeBuilderPtr_->createVector())) );

  mpdeLoaderPtr_->registerAppNextStoVec ( rcp(appBuilder_.createStoreVector()) );
  mpdeLoaderPtr_->registerAppCurrStoVec ( rcp(appBuilder_.createStoreVector()) );
  mpdeLoaderPtr_->registerAppLastStoVec ( rcp(appBuilder_.createStoreVector()) );
  
  mpdeLoaderPtr_->registerAppNextLeadCurrentVec ( rcp(appBuilder_.createLeadCurrentVector()) );
  mpdeLoaderPtr_->registerAppLeadCurrentQVec ( rcp(appBuilder_.createLeadCurrentVector()) );
  mpdeLoaderPtr_->registerAppNextJunctionVVec ( rcp(appBuilder_.createLeadCurrentVector()) );
  
  mpdeLoaderPtr_->registerAppdQdx( rcp(appBuilder_.createMatrix()) );
  mpdeLoaderPtr_->registerAppdFdx( rcp(appBuilder_.createMatrix()) );

  mpdeLoaderPtr_->registerMPDEdQdx( rcp_dynamic_cast<Xyce::Linear::BlockMatrix>(rcp(mpdeBuilderPtr_->createMatrix())) );

  mpdeLoaderPtr_->registerMPDEdFdx( rcp_dynamic_cast<Xyce::Linear::BlockMatrix>(rcp(mpdeBuilderPtr_->createMatrix())) );

  // Construct Solvers, etc. for MPDE Phase
  //-----------------------------------------

  // modify time integrator parameters for MPDE phase
  Xyce::dout() << "in setupMPDE: new period = " << period_ << std::endl;

  tiaMPDEParams_.initialTimeStep = period_;

  // if the user requested that a few steps be taken without LTE control,
  // then assert that option here.
  if( nonLteStepsGiven_ )
  {
    tiaMPDEParams_.errorAnalysisOption = Xyce::TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES;
    tiaMPDEParams_.errorAnalysisOptionResetCount = nonLteSteps_;
  }

  // Device package registrations.
  deviceManager_.registerNonlinearSolver( &nonlinearManager_ );
  deviceManager_.registerAnalysisManager( &analysisManager_ );

  nonlinearManager_.setLinSolOptions(linSolOptionBlock_);

  //hack needed by TIA initialization currently
  mpdeBuilderPtr_->registerPDSManager( &pdsManager_ );

  mpdeLinearSystem_.registerPDSManager( &pdsManager_ );
  mpdeLinearSystem_.registerBuilder(dynamic_cast<Xyce::Linear::Builder*>(&*mpdeBuilderPtr_));

  mpdeLinearSystem_.initializeSystem();

  // Initialization of Solvers, etc. for MPDE Phase
  analysisManager_.initializeSolverSystem(tiaMPDEParams_, *mpdeLoaderPtr_, mpdeLinearSystem_, nonlinearManager_, deviceManager_);

  nonlinearManager_.initializeAll(analysisManager_, analysisManager_.getNonlinearEquationLoader(), mpdeLinearSystem_, *analysisManager_.getDataStore(), pdsManager_, initialConditionsManager_, analysisManager_.getOutputManagerAdapter().getOutputManager(), topology_);

  if (DEBUG_MPDE)
    Xyce::dout() << Xyce::section_divider << std::endl;

  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::runMPDEProblem_
// Purpose       : Actual execution of MPDE problem
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::runMPDEProblem_()
{
  bool returnValue = true;

  transientNowCanOutputTrueMPDEResults_ = true;
  
  if (DEBUG_MPDE)
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  N_MPDE_Manager::runMPDEProblem_" << std::endl;

  *(analysisManager_.getDataStore()->nextSolutionPtr) = *mpdeICVectorPtr_;
  *(analysisManager_.getDataStore()->nextStatePtr) = *mpdeICStateVectorPtr_;
  *(analysisManager_.getDataStore()->daeQVectorPtr) = *mpdeICQVectorPtr_;
  *(analysisManager_.getDataStore()->nextStorePtr) = *mpdeICStoreVectorPtr_;

  Xyce::lout() << " ***** Beginning full MPDE simulation....\n" << std::endl;
  {
    Xyce::IO::ActiveOutput x(analysisManager_.getOutputManagerAdapter().getOutputManager());
    x.add(Xyce::IO::PrintType::MPDE, Xyce::Analysis::ANP_MODE_MPDE);
    x.add(pdsManager_.getPDSComm()->comm(), Xyce::Analysis::ANP_MODE_TRANSIENT);

    if (VERBOSE_TIME)
      tiaMPDEParams_.printParams(Xyce::lout(), nonlinearAnalysisMode(Xyce::Analysis::ANP_MODE_TRANSIENT));

    analysisManager_.setAnalysisMode(Xyce::Analysis::ANP_MODE_TRANSIENT);

    //Xyce::Analysis::Transient transient(analysisManager_, mpdeLinearSystem_, nonlinearManager_, *mpdeLoaderPtr_, topology_, initialConditionsManager_, restartManager_, &outputAdapter_, 0, this);
    Xyce::Analysis::Transient transient(analysisManager_, &mpdeLinearSystem_, nonlinearManager_, *mpdeLoaderPtr_, topology_, initialConditionsManager_, restartManager_, &outputAdapter_, 0, this);
    analysisManager_.pushActiveAnalysis(&transient);
    transient.setTIAParams(tiaMPDEParams_);
    transient.setAnalysisParams(Xyce::Util::OptionBlock());
    transient.setTimeIntegratorOptions(timeIntegratorOptionBlock_);
    transient.setNOOP(true);
    analysisManager_.getStepErrorControl().resetAll(tiaMPDEParams_);

    returnValue = transient.run();

    transientNeedsToLoadInitialConditionsAndInitializeProblem_ = false;

    analysisManager_.popActiveAnalysis();
  }

  if (DEBUG_MPDE)
    Xyce::dout() << Xyce::section_divider << std::endl;

  return returnValue;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::runTests_
// Purpose       : Testing MPDE objects
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::runTests_()
{
  Xyce::dout() << "N_MPDE_Manager::runTests_\n";

  //Finish setup of MPDE Builder
  //-----------------------------------------
  Xyce::dout() << "N_MPDE_Manager::runTests_[Generate Graphs]\n";
  mpdeBuilderPtr_->generateGraphs( *(pdsManager_.getMatrixGraph( Xyce::Parallel::JACOBIAN )) );

  Xyce::dout() << "N_MPDE_Manager::runTests_[Finished Graphs]\n";
  //-----------------------------------------

  //Setup MPDE Loader
  //-----------------------------------------
  Xyce::dout() << "N_MPDE_Manager::runTests_[Construct Loader]\n";
  mpdeLoaderPtr_ = new N_MPDE_Loader( mpdeState_, deviceManager_, appLoader_, *mpdeDiscPtr_, *this, warpMPDEPhasePtr_ );
  mpdeLoaderPtr_->setFastTimes( fastTimes_ );

  //Include these for the loader to use with devices
  mpdeLoaderPtr_->registerAppNextVec( rcp(appBuilder_.createVector()) );
  mpdeLoaderPtr_->registerAppCurrVec( rcp(appBuilder_.createVector()) );
  mpdeLoaderPtr_->registerAppLastVec( rcp(appBuilder_.createVector()) );

  mpdeLoaderPtr_->registerAppNextStaVec( rcp(appBuilder_.createStateVector()) );
  mpdeLoaderPtr_->registerAppCurrStaVec( rcp(appBuilder_.createStateVector()) );
  mpdeLoaderPtr_->registerAppLastStaVec( rcp(appBuilder_.createStateVector()) );
  mpdeLoaderPtr_->registerAppNextStoVec( rcp(appBuilder_.createStoreVector()) );
  mpdeLoaderPtr_->registerAppCurrStoVec( rcp(appBuilder_.createStoreVector()) );
  mpdeLoaderPtr_->registerAppdQdx( rcp(appBuilder_.createMatrix()) );
  mpdeLoaderPtr_->registerAppdFdx( rcp(appBuilder_.createMatrix()) );
  Xyce::dout() << "N_MPDE_Manager::runTests_[Finished Loader]\n";
  //-----------------------------------------

  //Test Construction of "MPDE Size" Objects
  //-----------------------------------------
  Xyce::dout() << "N_MPDE_Manager::runTests_[Construct Vectors and Matrices]\n";

  Xyce::Linear::Vector *Q(mpdeBuilderPtr_->createVector());
  Xyce::Linear::Vector *F(mpdeBuilderPtr_->createVector());
  Xyce::Linear::Vector *B(mpdeBuilderPtr_->createVector());
  Xyce::Linear::Vector *Res(mpdeBuilderPtr_->createVector());
  Xyce::Linear::Vector *nextX(mpdeBuilderPtr_->createVector());
  Xyce::Linear::Vector *currX(mpdeBuilderPtr_->createVector());
  Xyce::Linear::Vector *lastX(mpdeBuilderPtr_->createVector());

  Xyce::Linear::Vector *nextS(mpdeBuilderPtr_->createStateVector());
  Xyce::Linear::Vector *currS(mpdeBuilderPtr_->createStateVector());
  Xyce::Linear::Vector *lastS(mpdeBuilderPtr_->createStateVector());
  Xyce::Linear::Vector *dSdt(mpdeBuilderPtr_->createStateVector());
  Xyce::Linear::Vector *nextStore(mpdeBuilderPtr_->createStoreVector());
  Xyce::Linear::Vector *currStore(mpdeBuilderPtr_->createStoreVector());
  Xyce::Linear::Vector *lastStore(mpdeBuilderPtr_->createStoreVector());

  Xyce::Linear::Vector *nextLeadCurrent(mpdeBuilderPtr_->createLeadCurrentVector());
  Xyce::Linear::Vector *nextLeadCurrentQ(mpdeBuilderPtr_->createLeadCurrentVector());
  Xyce::Linear::Vector *nextJunctionV(mpdeBuilderPtr_->createLeadCurrentVector());
  
  Xyce::Linear::Matrix *dQdx(mpdeBuilderPtr_->createMatrix());
  Xyce::Linear::Matrix *dFdx(mpdeBuilderPtr_->createMatrix());
  Xyce::Linear::Matrix *Jac(mpdeBuilderPtr_->createMatrix());

  Xyce::Linear::BlockVector *bQ = dynamic_cast<Xyce::Linear::BlockVector *>(Q);
  Xyce::Linear::BlockVector *bF = dynamic_cast<Xyce::Linear::BlockVector *>(F);
  Xyce::Linear::BlockVector *bRes = dynamic_cast<Xyce::Linear::BlockVector *>(Res);
  Xyce::Linear::BlockVector *bNextX = dynamic_cast<Xyce::Linear::BlockVector *>(nextX);

  Xyce::Linear::BlockMatrix *bdQdx = dynamic_cast<Xyce::Linear::BlockMatrix *>(dQdx);
  Xyce::Linear::BlockMatrix *bdFdx = dynamic_cast<Xyce::Linear::BlockMatrix *>(dFdx);
  Xyce::Linear::BlockMatrix *bJac = dynamic_cast<Xyce::Linear::BlockMatrix *>(Jac);

  bNextX->print(Xyce::dout());
  bJac->print(Xyce::dout());
  
  Xyce::dout() << "N_MPDE_Manager::runTests_[Finished Vectors and Matrices]\n";

  //-----------------------------------------

  //Test Load of MPDE Objects
  //-----------------------------------------
  Xyce::dout() << "N_MPDE_Manager::runTests_[Load Vectors]\n";
  mpdeLoaderPtr_->loadDAEVectors(
      nextX, currX, lastX,
      nextS, currS, lastS, dSdt,
      nextStore, currStore, lastStore,
      nextLeadCurrent, nextLeadCurrentQ, 
      nextJunctionV, 
      Q, F, B,
      0, 0 );
  Xyce::dout() << "N_MPDE_Manager::runTests_[Load Matrices]\n";

  mpdeLoaderPtr_->loadDAEMatrices( &*nextX, &*nextS, &*dSdt, &*nextStore, &*dQdx, &*dFdx );
  Xyce::dout() << "N_MPDE_Manager::runTests_[Finished Loads]\n";

  bQ->print(Xyce::dout());
  bF->print(Xyce::dout());

  bdQdx->print(Xyce::dout());
  bdFdx->print(Xyce::dout());
  //-----------------------------------------

  //Pretend to solve DCOP problems
  //-----------------------------------------
  Res->update( +1.0, *F );
  bRes->print(Xyce::dout());

  Jac->add( *dFdx );
  bJac->print(Xyce::dout());

  delete Q;
  delete F;
  delete B;
  delete Res;
  delete nextX;
  delete currX;
  delete lastX;
  delete nextS;
  delete currS;
  delete lastS;
  delete dSdt;
  delete nextStore;
  delete currStore;
  delete lastStore;
  delete nextLeadCurrent;
  delete nextLeadCurrentQ;
  delete nextJunctionV;

  delete dQdx;
  delete dFdx;
  delete Jac;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::checkPeriodicity
// Purpose       : Return's an L2 norm of the difference of the initial
//                 condition signal's over two periods (if 2 were
//                 calculated) or just the endpoints if only one
//                 period was calculated.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, 1437, Electrical and Microsystem Simulation
// Creation Date : 06/13/2008
//-----------------------------------------------------------------------------
double N_MPDE_Manager::checkPeriodicity_()
{
  double returnValue = 0.0;

  // get access to time history of simulation
  Xyce::TimeIntg::DataStore *dsPtr = analysisManager_.getDataStore();
  int numPoints = dsPtr->timeSteps.size();

  // Note, if we integrated more then one period in the runTransientIC()
  // function, the we'll try to pull in points only from the second period.
  if (initialCondition_ == MPDE_IC_TWO_PERIOD)
  {
    double startUpOffset = tiaMPDEParams_.initialTime - period_;
    int lastJused = 0;
    for( int i=0; i < size_; i++ )
    {
      for( int j=lastJused; j < numPoints; j++ )
      {
        if( (fastTimes_[i] >= (dsPtr->timeSteps[j] - startUpOffset)) &&
            (fastTimes_[i] <  (dsPtr->timeSteps[j+1] - startUpOffset)) )
        {
          // found time points in the first period that bracket this one
          if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS) )
          {
            Xyce::dout() << "For fastTime_[" << i << "] = " << fastTimes_[i] <<
              " Found bracketing times : " << (dsPtr->timeSteps[j] - startUpOffset) << ", "
              << (dsPtr->timeSteps[j+1] - startUpOffset)
              << " at j = " << j;
          }

          // assume we're closest to the first point, but check if we're
          // really closest to the second one
          int interpolationPoint = j;
          if( fabs( fastTimes_[i] - (dsPtr->timeSteps[j] - startUpOffset)) >
              fabs( fastTimes_[i] - (dsPtr->timeSteps[j+1] - startUpOffset)) )
          {
            interpolationPoint = j+1;
          }
          if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS) )
            Xyce::dout() << " closest to point " << interpolationPoint << " index used = " << indicesUsed_[i] << std::endl;

          Xyce::Linear::Vector *thisPeriod = dsPtr->fastTimeSolutionVec[indicesUsed_[i]];
          Xyce::Linear::Vector *lastPeriod = dsPtr->fastTimeSolutionVec[interpolationPoint];
          Xyce::Linear::Vector *scratchVec = thisPeriod->cloneCopyVector();
          scratchVec->update( -1.0, *lastPeriod, 1.0 );
          scratchVec->infNorm(&returnValue);
          delete scratchVec;
          Xyce::dout() << i << " returnValue = " << returnValue << std::endl;
        }
      }
    }
  }

  return returnValue;
}

void
MPDEOutputAdapter::outputMPDE(
  double                        time,
  const Xyce::Linear::Vector *  solution_vector)
{
  if (dynamic_cast<const Xyce::Linear::BlockVector *>(solution_vector))
    outputManagerAdapter_.outputMPDE(time, mpdeManager_.getFastTimePoints(), dynamic_cast<const Xyce::Linear::BlockVector &>(*solution_vector));
}

bool
MPDEOutputAdapter::outputFunkyMPDE()
{
  // If we are actually in the MPDE phase, rather than the initial
  // condition, then output here.
  if (mpdeManager_.getOutputInterpolateMPDE() && mpdeManager_.getTransientNowCanOutputTrueMPDEResults())
  {
    if (!mpdeManager_.getWaMPDEFlag())
    {
      const std::vector<double> &fastTimes = mpdeManager_.getFastTimePoints();
      analysisManager_.getWorkingIntegrationMethod().printMPDEOutputSolution(
        outputManagerAdapter_, analysisManager_.getStepErrorControl().currentTime,
        analysisManager_.getDataStore()->currSolutionPtr,
        fastTimes );
    }
    else
    {
      const std::vector<double> &fastTimes = mpdeManager_.getFastTimePoints();
      int phiGID = mpdeManager_.getPhiGID();
      analysisManager_.getWorkingIntegrationMethod().printWaMPDEOutputSolution(
        outputManagerAdapter_, analysisManager_.getStepErrorControl().currentTime,
        analysisManager_.getDataStore()->currSolutionPtr,
        fastTimes, phiGID );
    }
    return true;
    
  }
  else
  {
    return false;
  }
}


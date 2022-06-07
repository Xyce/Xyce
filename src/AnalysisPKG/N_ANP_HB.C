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

//-----------------------------------------------------------------------------
// Purpose       : HB analysis functions.
// Special Notes :
// Creator       : Todd Coffey, 1414, Ting Mei, 1437
// Creation Date : 07/23/08
//
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_DCSweep.h>
#include <N_ANP_HB.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_Report.h>
#include <N_ANP_Transient.h>
#include <N_DEV_DeviceMgr.h>
#include <N_IO_ActiveOutput.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_IO_PrintTypes.h>
#include <N_LAS_BlockSystemHelpers.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_HBBuilder.h>
#include <N_LAS_HBPrecondFactory.h>
#include <N_LAS_HBSolverFactory.h>
#include <N_LAS_PrecondFactory.h>
#include <N_LAS_System.h>
#include <N_LOA_HBLoader.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_NLS_Manager.h>
#include <N_PDS_ParMap.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_fwd.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_WorkingIntegrationMethod.h>
#include <N_UTL_APFT.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FFTInterface.hpp>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_MachDepParams.h>
#include <N_UTL_Math.h>
#include <N_UTL_Timer.h>

#include <Teuchos_BLAS.hpp>
#include <Teuchos_Utils.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseHelpers.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

#include <N_TOP_Topology.h>
#include <N_PDS_Comm.h>

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : HB::HB( AnalysisManager * )
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
HB::HB(
  AnalysisManager &                     analysis_manager,
  Linear::System &                      linear_system,
  Nonlinear::Manager &                  nonlinear_manager,
  Loader::Loader &                      loader,
  Device::DeviceMgr &                   device_manager,
  Linear::Builder &                     builder,
  Topo::Topology &                      topology,
  IO::InitialConditionsManager &        initial_conditions_manager,
  IO::RestartMgr &                      restart_manager)
  : AnalysisBase(analysis_manager, "HB"),
    StepEventListener(&analysis_manager),
    analysisManager_(analysis_manager),
    loader_(loader),
    linearSystem_(linear_system),
    nonlinearManager_(nonlinear_manager), 
    deviceManager_(device_manager),
    builder_(builder),
    topology_(topology),
    initialConditionsManager_(initial_conditions_manager),
    restartManager_(restart_manager),
    hbLoaderPtr_(0),
    hbBuilderPtr_(),
    hbLinearSystem_(0),
    isPaused(false),
    startDCOPtime(0.0),
    endTRANtime(0.0),
    isTransient_(false),
    test_(false),
    size_(21),
    freqs_(),
    freqsGiven_(false),
    period_(1.0),
    relErrorTol_(0.0),
    startUpPeriods_(0),
    startUpPeriodsGiven_(false),
    saveIcData_(false),
    useStartupICs_(false),
    taHB_(1),
    hbOsc_(false),
    refID_(-1),
    refNode_(""),
    numTimePts_(0),
    voltLimFlag_(1),
    intmodMax_(0),
    loadTimeB_(1),
    method_("APFT"),
    intmodMaxGiven_(false),
    selectHarm_("BOX"),
    selectHarmGiven_(false),
    solverFactory_(0),
    precFactory_(0),
    resetForStepCalledBefore_(false)
{
  pdsMgrPtr_ = analysisManager_.getPDSManager();

  if (taHB_ == 1)
    analysisManager_.getOutputManagerAdapter().setTaHBSpecified(true);
}

//-----------------------------------------------------------------------------
// Function      : HB::~HB
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
HB::~HB()
{
  delete hbLoaderPtr_;
  delete hbLinearSystem_;
  delete solverFactory_;
  delete precFactory_;
}

//-----------------------------------------------------------------------------
// Function      : HB::notify
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Baur
// Creation Date : 
//-----------------------------------------------------------------------------
void HB::notify(const StepEvent &event) 
{
  if (event.state_ == StepEvent::STEP_STARTED)
  {
    AnalysisBase::resetForStepAnalysis();

    if (resetForStepCalledBefore_)
    {
      goodSolutionVec_.clear();
      goodStateVec_.clear();
      goodQVec_.clear();
      goodStoreVec_.clear();

      analysisManager_.setNextOutputTime(0.0);

      nonlinearManager_.resetAll(Nonlinear::DC_OP);
      nonlinearManager_.setMatrixFreeFlag( false );
      nonlinearManager_.setLinSolOptions( saved_lsOB_ );
      nonlinearManager_.registerSolverFactory( NULL );

      // un-set the fast source flag in the devices
      std::vector<std::string> srcVec;
      deviceManager_.deRegisterFastSources( srcVec );

      analysisManager_.initializeSolverSystem(TimeIntg::TIAParams(), loader_, linearSystem_, nonlinearManager_, deviceManager_);

      deviceManager_.initializeAll(linearSystem_);
      deviceManager_.setMPDEFlag( false );

      nonlinearManager_.initializeAll(
        analysisManager_,
        analysisManager_.getNonlinearEquationLoader(),
        linearSystem_,
        *analysisManager_.getDataStore(),
        *analysisManager_.getPDSManager(),
        initialConditionsManager_,
        analysisManager_.getOutputManagerAdapter().getOutputManager(),
        topology_);

      TimeIntg::DataStore * dsPtr = analysisManager_.getDataStore();
      hbLoaderPtr_->loadDeviceErrorWeightMask(dsPtr->deviceErrorWeightMask_);

      analysisManager_.getXyceTranTimer().resetStartTime();
    }

    resetForStepCalledBefore_ = true;
  }
}

//-----------------------------------------------------------------------------
// Function      : HB::getDoubleDCOPStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 3/24/2014
//-----------------------------------------------------------------------------
int HB::getDoubleDCOPStep() const
{
  if (currentAnalysisObject_)
    return currentAnalysisObject_->getDoubleDCOPStep();
  else
    return 0;
}

//-----------------------------------------------------------------------------
// Function      : HB::getDCOPFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 3/24/2014
//-----------------------------------------------------------------------------
bool HB::getDCOPFlag() const
{
  if (currentAnalysisObject_ && isTransient_ )
      return currentAnalysisObject_->getDCOPFlag ();

  // DCSweep is a special case in the HB analysis type.  The HB calculation is
  // performed through a DC analysis object.  However, this function (getDCOPFlag)
  // is called (ultimately) from the device package to determine a bunch of
  // state-dependent load decisions.  For an HB calculation, it should NOT
  // do a DCOP load.  It needs to do a full transient load.
  else 
  {
    if (taHB_ == 0 )
      return true;
    else
      return false;

  }
}

//-----------------------------------------------------------------------------
// Function      : HB::finalExpressionBasedSetup()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/4/2021
//-----------------------------------------------------------------------------
void HB::finalExpressionBasedSetup()
{
  return;
}

//-----------------------------------------------------------------------------
// Function      : HB::doRun()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::doRun()
{
  return doInit() && doLoopProcess() && doFinish();
}

//-----------------------------------------------------------------------------
// Function      : HB::init()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::doInit()
{
  Xyce::lout() << " ***** Running HB initial conditions....\n" << std::endl;

  Stats::StatTop _initStepStat("Initialization");
  Stats::TimeBlock _initStepTimer(_initStepStat);

  if (DEBUG_HB)
    Xyce::dout() << std::endl
                 << section_divider << std::endl
                 << "  HB::init()" << std::endl;

  setDoubleDCOPEnabled(loader_.isPDESystem());

  initializeOscOut( );
  setFreqPoints_();



  if (method_ == "AFM")
     mapFreqs_();


//  period_ = 1.0/freqPoints_[(size_ - 1)/2 + 1];

  int posFreqSize = (size_ - 1)/2;
   
  std::vector<double> freqPos;

  freqPos.resize(posFreqSize);

  for (int i=0; i < posFreqSize; i++)
  {
    freqPos[i] = fabs( freqPoints_[i + posFreqSize + 1]);

//    Xyce::dout() << "freqPos =" <<  freqPos[i] << std::endl;

  }

   double minFreq = *std::min_element( freqPos.begin(), freqPos.end());

   period_ = 1.0/minFreq;

  if (DEBUG_HB)
    Xyce::dout() << "HB period =" <<  period_ << std::endl;

  fastTimes_.resize(size_);

  goodTimePoints_.resize(size_);

  if ((freqs_.size() == 1) || ((taHB_ != 1) && (method_ == "AFM")) )
  {
    double TimeStep = period_/size_;
    timeSteps_.push_back( TimeStep );

    for( int i = 0; i < size_; ++i )
    {
      fastTimes_[i] = i*TimeStep;
    }
  }
  else
  {
    setTimePoints_();
  }

  goodTimePoints_ = fastTimes_;

  {
    Stats::StatTop _guessStepStat("Initial Guess");
    Stats::TimeBlock _guessStepTimer(_guessStepStat);
  
    setInitialGuess();
  }

  goodTimePoints_ = fastTimes_;

  // now that we have size_, continue with the initialization of objects for HB

  hbBuilderPtr_ = rcp(new Linear::HBBuilder(size_, hbOsc_));

  if (DEBUG_HB)
    Xyce::dout() << "HB::init():  Generate Maps\n";

  {
    Stats::StatTop _setupStepStat("Setup Maps/Graphs");
    Stats::TimeBlock _setupStepTimer(_setupStepStat);

    hbBuilderPtr_->generateMaps( rcp(pdsMgrPtr_->getParallelMap( Parallel::SOLUTION ), false),
                             rcp(pdsMgrPtr_->getParallelMap( Parallel::SOLUTION_OVERLAP_GND ), false) );
    hbBuilderPtr_->generateStateMaps( rcp(pdsMgrPtr_->getParallelMap( Parallel::STATE ),false) );
    hbBuilderPtr_->generateStoreMaps( rcp(pdsMgrPtr_->getParallelMap( Parallel::STORE ),false) );
    hbBuilderPtr_->generateLeadCurrentMaps( rcp(pdsMgrPtr_->getParallelMap( Parallel::LEADCURRENT ),false) );
    hbBuilderPtr_->generateGraphs( *pdsMgrPtr_->getMatrixGraph( Parallel::JACOBIAN ));
  }

  HBICVectorPtr_ = hbBuilderPtr_->createTimeDomainBlockVector();
  HBICStateVectorPtr_ = hbBuilderPtr_->createTimeDomainStateBlockVector();
  HBICStoreVectorPtr_ = hbBuilderPtr_->createTimeDomainStoreBlockVector();
  HBICQVectorPtr_ = hbBuilderPtr_->createTimeDomainBlockVector();

  HBICVectorFreqPtr_ = hbBuilderPtr_->createExpandedRealFormTransposeBlockVector();

  // set the fast source flag in the devices
  std::vector<std::string> srcVec;
  deviceManager_.registerFastSources(pdsMgrPtr_->getPDSComm()->comm(), srcVec );

  //nonlinearManager_.getHBOptions(saved_nlHBOB_);

  // Create HB Loader.
  delete hbLoaderPtr_;

  hbLoaderPtr_ = new Loader::HBLoader(deviceManager_, builder_, refID_, hbOsc_);
  hbLoaderPtr_->registerHBBuilder(hbBuilderPtr_);
  hbLoaderPtr_->registerAppLoader( rcp(&loader_, false) );

  // Create DFT for HB Loader
  // NOTE:  For single-tone HB the DFT will probably be a FFT, for multi-tone a specialized 
  //        implementation of the Util::DFTInterfaceDecl will need to be made and registered with
  //        the HB loader.

  ftInData_.resize( size_ );
  ftOutData_.resize( size_ +1 );
  iftInData_.resize( size_  +1 );
  iftOutData_.resize( size_ );
  if ((freqs_.size() == 1 ) || ((taHB_ != 1) && (method_ == "AFM")) )
  {
    if (ftInterface_ == Teuchos::null)
    {
      ftInterface_ = Teuchos::rcp( new N_UTL_FFTInterface<std::vector<double> >( size_ ) );
      ftInterface_->registerVectors( ftInData_, &ftOutData_, iftInData_, &iftOutData_ );
    } 
    else if (ftInterface_->getFFTInterface()->getSignalLength() != size_)
    {
      ftInterface_ = Teuchos::rcp( new N_UTL_FFTInterface<std::vector<double> >( size_ ) );
      ftInterface_->registerVectors( ftInData_, &ftOutData_, iftInData_, &iftOutData_ );
    }
    hbLoaderPtr_->registerDFTInterface( ftInterface_->getFFTInterface() );
  }
  else
  {
    createFT_();

    dftInterface_ = Teuchos::rcp( new N_UTL_APFT<std::vector<double> >( idftMatrix_, dftMatrix_ ) );
    dftInterface_->registerVectors( Teuchos::rcp( &ftInData_, false ), Teuchos::rcp( &ftOutData_, false ), Teuchos::rcp( &iftInData_, false ), Teuchos::rcp( &iftOutData_, false ) );
    hbLoaderPtr_->registerDFTInterface( dftInterface_ );
  }

  if (taHB_==1)
  {
    // Pick up IC data from the initial transient.
    for (int i=0 ; i<size_ ; ++i)
    {
      if (DEBUG_HB && isActive(Diag::HB_FAST_TIMES))
        Xyce::dout() << "HB::init():  Loading initial condition data from time: fastTimes_["
                     << i << "] = " << fastTimes_[i] << std::endl;

      HBICVectorPtr_->block(i) = *(goodSolutionVec_[i]);
      HBICStateVectorPtr_->block(i) = *(goodStateVec_[i]);
      HBICQVectorPtr_->block(i) = *(goodQVec_[i]);
      HBICStoreVectorPtr_->block(i) = *(goodStoreVec_[i]);
    }

    hbLoaderPtr_->permutedFFT(*HBICVectorPtr_, &*HBICVectorFreqPtr_);

    if (hbOsc_)
    {
      const std::vector<int> & augmentedLIDs = hbBuilderPtr_->getAugmentedLIDs();

      if ( augmentedLIDs.size() )
        (*HBICVectorFreqPtr_)[(augmentedLIDs)[0 ]] = 1.0;
    }

  }
  else if (taHB_== 2)

  {
    for (int i=0 ; i<size_ ; ++i)
    {

      HBICVectorPtr_->block(i) = *dcOpSolVecPtr_;
      HBICStateVectorPtr_->block(i) = *dcOpStateVecPtr_;
      HBICQVectorPtr_->block(i) = *dcOpQVecPtr_;
      HBICStoreVectorPtr_->block(i) =  *dcOpStoreVecPtr_;
    }

    hbLoaderPtr_->permutedFFT(*HBICVectorPtr_, &*HBICVectorFreqPtr_);

    if (hbOsc_)
    {
      const std::vector<int> & augmentedLIDs = hbBuilderPtr_->getAugmentedLIDs();

      if ( augmentedLIDs.size() )
        (*HBICVectorFreqPtr_)[(augmentedLIDs)[0 ]] = 1.0;
    }
  }

  if (DEBUG_HB && isActive(Diag::HB_PRINT_VECTORS))
  {
    Xyce::dout() << "HB Initial Condition Solution!\n";
    HBICVectorPtr_->print(std::cout);
    Xyce::dout() << "HB Initial Condition State Vector!\n";
    HBICStateVectorPtr_->print(std::cout);
    Xyce::dout() << "HB Initial Condition Store Vector!\n";
    HBICStoreVectorPtr_->print(std::cout);
  }

  if (method_ == "AFM")
  {

    if (ftInterface_ == Teuchos::null)
    {
      ftInterface_ = Teuchos::rcp( new N_UTL_FFTInterface<std::vector<double> >( size_ ) );
      ftInterface_->registerVectors( ftInData_, &ftOutData_, iftInData_, &iftOutData_ );
    } 
    else if (ftInterface_->getFFTInterface()->getSignalLength() != size_)
    {
      ftInterface_ = Teuchos::rcp( new N_UTL_FFTInterface<std::vector<double> >( size_ ) );
      ftInterface_->registerVectors( ftInData_, &ftOutData_, iftInData_, &iftOutData_ );
    }
    hbLoaderPtr_->registerDFTInterface( ftInterface_->getFFTInterface() );
  }

  //Destroy Solvers, etc. from IC phase and prepare for HB
  //-----------------------------------------

  deviceManager_.setVoltageLimiterFlag(voltLimFlag_);

  deviceManager_.setMPDEFlag( true );
  analysisManager_.resetSolverSystem();

  //-----------------------------------------

  //Finish setup of HB Loader
  //-----------------------------------------
  goodTimePoints_.resize(size_+1);
  goodTimePoints_[size_] = period_;

  {
    Stats::StatTop _nlsStepStat("Setup Solvers");
    Stats::TimeBlock _nlsStepTimer(_nlsStepStat);

    hbLoaderPtr_->setHBFreqs(freqPoints_);
    hbLoaderPtr_->setFastTimes(goodTimePoints_);

    //-----------------------------------------
    //Construct Solvers, etc. for HB Phase
    //-----------------------------------------
    delete hbLinearSystem_;
    hbLinearSystem_ = new Linear::System();
    //-----------------------------------------

    //hack needed by TIA initialization currently
    hbBuilderPtr_->registerPDSManager( pdsMgrPtr_ );
    hbLinearSystem_->registerPDSManager( pdsMgrPtr_ );
    hbLinearSystem_->registerBuilder( &*hbBuilderPtr_ );

    hbLinearSystem_->initializeSystem();

    nonlinearManager_.setLinSolOptions( saved_lsHBOB_ );
    nonlinearManager_.setMatrixFreeFlag( true );

    // Let the HB loader know that the application of the operator is matrix free
    hbLoaderPtr_->setMatrixFreeFlag( true );

    if ( (method_ == "AFM") || !loadTimeB_ )  
      hbLoaderPtr_->setLoadTimeBFlag( false );

    if (!solverFactory_)
    {
      // Generate the HB solver factory.
      solverFactory_ = new Linear::HBSolverFactory( builder_ );
    }

    // Register application loader with solver factory
    solverFactory_->registerHBLoader( rcp(hbLoaderPtr_, false) );
    solverFactory_->registerHBBuilder( hbBuilderPtr_ );
    solverFactory_->setFastTimes( goodTimePoints_ );
    solverFactory_->setHBFreqs( freqPoints_ );

    if (!precFactory_)
    {
      // Generate the HB preconditioner factory.
      precFactory_ = new Linear::HBPrecondFactory(saved_lsHBOB_, builder_);
    }

    // Register application loader with preconditioner factory
    precFactory_->registerHBLoader( rcp(hbLoaderPtr_, false) );
    precFactory_->registerHBBuilder( hbBuilderPtr_ );
    precFactory_->setFastTimes( goodTimePoints_ );
    precFactory_->setTimeSteps( timeSteps_ );
    precFactory_->setHBFreqs( freqPoints_ );

    nonlinearManager_.registerSolverFactory( solverFactory_ );
    nonlinearManager_.registerPrecondFactory( precFactory_ );

    precFactory_->setHBOsc(hbOsc_);

  //Initialization of Solvers, etc. for HB Phase
    analysisManager_.initializeSolverSystem(TimeIntg::TIAParams(), *hbLoaderPtr_, *hbLinearSystem_, nonlinearManager_, deviceManager_);

    nonlinearManager_.initializeAll(analysisManager_, analysisManager_.getNonlinearEquationLoader(), *hbLinearSystem_, *analysisManager_.getDataStore(), *analysisManager_.getPDSManager(), initialConditionsManager_, analysisManager_.getOutputManagerAdapter().getOutputManager(), topology_);
  
    nonlinearManager_.setAnalysisMode(nonlinearAnalysisMode(ANP_MODE_HB));
  }

  TimeIntg::DataStore * dsPtr = analysisManager_.getDataStore();
  hbLoaderPtr_->loadDeviceErrorWeightMask(dsPtr->deviceErrorWeightMask_);
  dsPtr->allocateHBVectors();

  if (DEBUG_HB)
  {
    Xyce::dout() << section_divider << std::endl;
  }

  return true;

}

//-----------------------------------------------------------------------------
// Function      : HB::loopProcess()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::doLoopProcess()
{
  bool returnValue = true;

  Xyce::lout() << " ***** Beginning full HB simulation....\n" << std::endl;

  if (DEBUG_HB)
  {
    Xyce::dout() << std::endl
                 << section_divider << std::endl
                 << "  HB::loopProcess" << std::endl;
  }

  TimeIntg::DataStore * dsPtr = analysisManager_.getDataStore();
  *(dsPtr->nextSolutionPtr) = *(HBICVectorFreqPtr_);
  *(dsPtr->nextStatePtr) = *(HBICStateVectorPtr_);
  *(dsPtr->nextStorePtr) = *(HBICStoreVectorPtr_);

  // try to run the problem
  {
    // add active outputters for both HB_FD and HB_TD output, even if only one of
    // those types of output has been requested in the netlist.
    Xyce::IO::ActiveOutput x(analysisManager_.getOutputManagerAdapter().getOutputManager());
    x.add(Xyce::IO::PrintType::HB_FD, ANP_MODE_HB);
    x.add(Xyce::IO::PrintType::HB_TD, ANP_MODE_HB);

    //DCSweep dc_sweep(analysisManager_, linearSystem_, nonlinearManager_, loader_, topology_, initialConditionsManager_, this);
    DCSweep dc_sweep(analysisManager_, &linearSystem_, nonlinearManager_, loader_, topology_, initialConditionsManager_, this);
    currentAnalysisObject_ = &dc_sweep;
    analysisManager_.pushActiveAnalysis(&dc_sweep);
    returnValue = dc_sweep.run();

    // Add in simulation times
    accumulateStatistics_(dc_sweep);

    // print out analysis info
    Xyce::lout() << " ***** Harmonic Balance Computation Summary *****" << std::endl;
    dc_sweep.printLoopInfo( 0, 0 );

    analysisManager_.popActiveAnalysis();

    currentAnalysisObject_ = 0;
  }

  if (DEBUG_HB)
  {
    dout() << section_divider << std::endl;
  }

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : HB::processSuccessfulDCOP()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::processSuccessfulDCOP()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : HB::processSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::doProcessSuccessfulStep()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : HB::processFailedStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::doProcessFailedStep()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : HB::processFailedDCOP
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::processFailedDCOP()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : HB::finish
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::doFinish()
{
  stats_ = hbStatCounts_;

  return true;
}

bool HB::doHandlePredictor()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::finalVerboseOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::finalVerboseOutput()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setHBAnalysisParams
// Purpose       : Sets the HB sweep calculation parameters (from .HB)
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 7/30/08
//-----------------------------------------------------------------------------
bool
HB::setAnalysisParams(
  const Util::OptionBlock &     option_block)
{
  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; ++it)
  {
    if ((*it).uTag() == "FREQ")
    {
      freqs_ = (*it).getValue<std::vector<double> >();
      freqsGiven_ = true;
    }
  }

  if (freqs_[0] <= 0.0 )
  {
    Report::UserError() << "Frequency of oscillation " << freqs_[0] << " is less than or equal to zero, invalid .HB specification";
  }

  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    dout() << section_divider << std::endl
           << "HB transient simulation parameters" 
           //<< Util::push << std::endl
           << std::endl
           << "HB frequency = " << freqs_[0] << std::endl
           //<< Util::pop << std::endl;
           << std::endl;
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : HB::setHBIntParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 05/13/13
//-----------------------------------------------------------------------------
bool HB::setHBIntParams(const Util::OptionBlock & OB)
{
  for(Util::ParamList::const_iterator iterPL = OB.begin(), endPL = OB.end(); iterPL != endPL; ++iterPL )
  {
    ExtendedString tag = iterPL->tag();
    tag.toUpper();

    if (std::string(tag,0,7) == "NUMFREQ" ) 
    {

      size_ = iterPL->getImmutableValue<int>() *2 + 1;

      numFreqs_.push_back(size_);
    }
    else if ( tag == "STARTUPPERIODS" )
    {
      startUpPeriods_ = iterPL->getImmutableValue<int>();

      if (startUpPeriods_ > 0)
        startUpPeriodsGiven_ = true;
    }
    else if( tag == "SAVEICDATA" )
    {
      saveIcData_ = true;
    }
    else if( tag == "TEST" )
    {
      test_     = static_cast<bool> (iterPL->getImmutableValue<int>());
    }
    else if (tag == "DEBUGLEVEL" )
    {
      setHBDebugLevel(iterPL->getImmutableValue<int>());
    }
    else if ( tag == "TAHB" )
    {
      taHB_ = iterPL->getImmutableValue<int>();
       if (taHB_ != 1)
         analysisManager_.getOutputManagerAdapter().setTaHBSpecified(false);
    }
    else if ( tag == "VOLTLIM" )
    {
      voltLimFlag_ = static_cast<bool> (iterPL->getImmutableValue<int>());
    }
    else if ( tag == "INTMODMAX" )
    {
      intmodMax_ = iterPL->getImmutableValue<int>();

      if ( intmodMax_  > 0)
        intmodMaxGiven_ = true;
    }
    else if ( tag == "LOADTIMESOURCES") 
    {
      loadTimeB_ = static_cast<bool> (iterPL->getImmutableValue<int>());
    }
    else if ( tag == "METHOD" )
    {
      ExtendedString stringVal ( iterPL->stringValue() );
      stringVal.toUpper();
      method_ = stringVal;
    }
    else if ( tag == "NUMTPTS") 
    {
      numTimePts_ = iterPL->getImmutableValue<int>();
    }
    else if ( tag == "HBOSC") 
    {
      hbOsc_ = static_cast<bool> (iterPL->getImmutableValue<int>());
    }
    else if ( tag == "SELECTHARMS" )
    {
      ExtendedString stringVal ( iterPL->stringValue() );
      stringVal.toUpper();
      selectHarm_ = stringVal;

      if ( selectHarm_ != "" )
        selectHarmGiven_ = true;
    }
    else if ( tag == "REFNODE") 
    {
      ExtendedString stringVal ( iterPL->stringValue() );
      stringVal.toUpper();
      refNode_ = stringVal;
    }
    else
    {
      UserWarning(*this) << "Unrecognized HBINT option " << tag;
    }
  }

  if (numFreqs_.size() != 0 && numFreqs_.size() != freqs_.size() ) 
  {
    Report::UserError() << "The size of numFreq does not match the number of tones in .hb!";
  }


  if (numFreqs_.size() == 0 )
  {

    numFreqs_.resize( freqs_.size() ) ;
    for (int i=0; i < freqs_.size(); i++ )
    {
      numFreqs_[i] = size_;
    }

  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::setLinSol
// Purpose       : this is needed for .STEP to work with HB
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 7/12/2013
//-----------------------------------------------------------------------------
bool HB::setLinSol(const Util::OptionBlock & OB)
{
  // Save the non-HB linear solver option block
  saved_lsOB_ = OB;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::setTimeInt
// Purpose       : this is needed for .STEP to work with HB
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 7/12/2013
//-----------------------------------------------------------------------------
bool HB::setTimeInt(const Util::OptionBlock & OB)
{
  // Save the time integrator option block
  saved_timeIntOB_ = OB;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::setHBLinSol
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/13/13
//-----------------------------------------------------------------------------
bool HB::setHBLinSol(const Util::OptionBlock & OB, Linear::Builder &builder)
{
  // Save the HB linear solver option block
  saved_lsHBOB_ = OB;

  // Generate the HB preconditioner factory.
  precFactory_ = new Linear::HBPrecondFactory(OB, builder);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::isAnalysis
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/13/13
// Notes         : Alternatively, we could try to cast the analysis object
//               : However, this method is called a lot.
//-----------------------------------------------------------------------------
bool HB::isAnalysis( int analysis_type ) const
{
  return analysis_type == ANP_MODE_TRANSIENT && isTransient_;
}

//-----------------------------------------------------------------------------
// Function      : HB::prepareHBOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 08/20/07
//-----------------------------------------------------------------------------
void HB::prepareHBOutput(
  Linear::Vector & solnVecPtr,
  std::vector<double> & timePoints,
  std::vector<double> & freqPoints,
  RCP<Linear::BlockVector> & timeDomainSolnVec,
  RCP<Linear::BlockVector> & freqDomainSolnVecReal,
  RCP<Linear::BlockVector> & freqDomainSolnVecImag,
  RCP<Linear::BlockVector> & timeDomainLeadCurrentVec,
  RCP<Linear::BlockVector> & freqDomainLeadCurrentVecReal,
  RCP<Linear::BlockVector> & freqDomainLeadCurrentVecImaginary,
  RCP<Linear::BlockVector> & timeDomainJunctionVoltageVec,
  RCP<Linear::BlockVector> & freqDomainJunctionVoltageVecReal,
  RCP<Linear::BlockVector> & freqDomainJunctionVoltageVecImaginary
                         ) 
{
  Linear::BlockVector & blockSolVecPtr = dynamic_cast<Linear::BlockVector &>(solnVecPtr);
//  int tdsampleRate = tdsampleRate_;
  
  if ( numTimePts_ < size_ )
    numTimePts_ = size_ ;
  if ( numTimePts_ % 2 == 0)
    numTimePts_ = numTimePts_ + 1;
  
  int blockCounttd = numTimePts_;

 //  int blockCounttd = (size_-1)/2*tdsampleRate*2 + 1; 

  int N = blockSolVecPtr.blockCount(); 
  Teuchos::RCP<Parallel::ParMap> baseMap = Teuchos::rcp_const_cast<Parallel::ParMap>( hbBuilderPtr_->getBaseSolutionMap() );
  Teuchos::RCP<Parallel::ParMap> globalMap = Linear::createBlockParMap( blockCounttd, *baseMap );
  
  timeDomainSolnVec = Teuchos::rcp( Xyce::Linear::createBlockVector(blockCounttd, globalMap, baseMap ) );

  Teuchos::RCP<Linear::BlockVector> bLeadCurrentVecFreqPtr_ = hbLoaderPtr_->getLeadCurrentVecFreqPtr();
  freqPoints = freqPoints_;
  
  timePoints.resize(blockCounttd);

  if (numTimePts_== size_ )
  {
    timePoints = fastTimes_;

    if ( ( method_ == "AFM" ) && (freqs_.size() > 1 ) )
    {
      updateIFT_( timePoints);
//      iftInData_.resize( size_+1 );

      dftInterface_ = Teuchos::rcp( new N_UTL_APFT<std::vector<double> >( idftMatrix_, dftMatrix_ ) );
      dftInterface_->registerVectors( Teuchos::rcp( &ftInData_, false ), Teuchos::rcp( &ftOutData_, false ), Teuchos::rcp( &iftInData_, false ), Teuchos::rcp( &iftOutData_, false ) );
      hbLoaderPtr_->registerDFTInterface( dftInterface_ ); 
    }
  }
  else
  {

    double TimeStep = period_/blockCounttd;

    for( int i = 0; i < blockCounttd; ++i )
      timePoints[i] = i*TimeStep;

    ftInData_.resize( blockCounttd);
    ftOutData_.resize( blockCounttd +1 );
    iftInData_.resize( blockCounttd +1 );
    iftOutData_.resize( blockCounttd);


    if (freqs_.size() == 1 )
    {
      if  (ftInterface_->getFFTInterface()->getSignalLength() !=  blockCounttd)
      {
        ftInterface_ = Teuchos::rcp( new N_UTL_FFTInterface<std::vector<double> >(blockCounttd) );
        ftInterface_->registerVectors( ftInData_, &ftOutData_, iftInData_, &iftOutData_ );

        hbLoaderPtr_->registerDFTInterface( ftInterface_->getFFTInterface() );
      }

    }
    else
    {
      updateIFT_( timePoints);
      iftInData_.resize( size_+1 );

      dftInterface_ = Teuchos::rcp( new N_UTL_APFT<std::vector<double> >( idftMatrix_, dftMatrix_ ) );
      dftInterface_->registerVectors( Teuchos::rcp( &ftInData_, false ), Teuchos::rcp( &ftOutData_, false ), Teuchos::rcp( &iftInData_, false ), Teuchos::rcp( &iftOutData_, false ) );
      hbLoaderPtr_->registerDFTInterface( dftInterface_ );

    }
  }

  int blockCount = size_;
  Teuchos::RCP<Parallel::ParMap> globalMapfreq = Linear::createBlockParMap( blockCount, *baseMap );
  freqDomainSolnVecReal = Teuchos::rcp( Xyce::Linear::createBlockVector( blockCount, globalMapfreq, baseMap ) );
  freqDomainSolnVecImag = Teuchos::rcp( Xyce::Linear::createBlockVector( blockCount, globalMapfreq, baseMap ) ); 

  hbLoaderPtr_->permutedIFT(blockSolVecPtr, &*timeDomainSolnVec, numTimePts_);

  // Now copy over the frequency domain solution, real and imaginary parts separately, into the output vectors.

  for (int j=0; j<N; j++)
  {
    // See if this time-domain solution variable is owned by the local processor.
    // If so, this processor owns the entire j-th block of the blockSolVecPtr vector,
    // and the j-th entry of every block in the freqDomainSo
    // lnVec[Real/Imag] vector.
    int lid = baseMap->globalToLocalIndex( j );
    Linear::Vector& solBlock = blockSolVecPtr.block( j );

    std::vector<std::pair<double, double>> realList, imagList;
    int sizePos = (size_ - 1)/2;

    if ( ( selectHarm_ == "BOX") && ( lid >= 0 ) )
    { 

//      std::vector<std::pair<double, double>> realList, imagList;

      realList.resize(size_);
      imagList.resize(size_);

      for (int m=0; m<size_; m++)
      {

        if  ( m <=  sizePos )
        {
          realList[m] = std::make_pair( freqPoints_[m+sizePos] , solBlock[2*m]); 
          imagList[m] = std::make_pair( freqPoints_[m+sizePos] , solBlock[2*m+1]);
        }
        else
        {
          realList[m] = std::make_pair( freqPoints_[m-sizePos-1] , solBlock[2*m]); 
          imagList[m] = std::make_pair( freqPoints_[m-sizePos-1] , solBlock[2*m+1]);
        }
      }

      std::sort(realList.begin(),realList.end());

      std::sort(imagList.begin(),imagList.end());

      if (j == 0)
      {
        for (int i=0; i< freqPoints_.size(); i++)
        {
          freqPoints[i]  = realList[i].first;
        }
      } 
    }


    Linear::Vector& realVecRef =  freqDomainSolnVecReal->block((blockCount-1)/2);
    Linear::Vector& imagVecRef =  freqDomainSolnVecImag->block((blockCount-1)/2);

    if (lid >= 0)
    { 
      realVecRef[lid] = solBlock[0];  
      imagVecRef[lid] = solBlock[1];
    }

    for (int i=1; i <= (size_-1)/2; ++i)
//  for (int i=1; i <= (blockCount-1)/2; ++i)
    {
      Linear::Vector& realVecRef_neg =  freqDomainSolnVecReal->block((blockCount-1)/2 - i);
      Linear::Vector& imagVecRef_neg =  freqDomainSolnVecImag->block((blockCount-1)/2 - i);
      Linear::Vector& realVecRef_pos =  freqDomainSolnVecReal->block((blockCount-1)/2 + i);
      Linear::Vector& imagVecRef_pos =  freqDomainSolnVecImag->block((blockCount-1)/2 + i);

      if (lid >= 0)
      {

        if ( selectHarm_ ==  "BOX" )
        {
        realVecRef_neg[lid] = realList[ sizePos - i].second;
        imagVecRef_neg[lid] = imagList[ sizePos - i].second;
        realVecRef_pos[lid] = realList[ sizePos + i].second; 
        imagVecRef_pos[lid] = imagList[ sizePos + i].second;
           
        }
        else
        {
        realVecRef_neg[lid] = solBlock[ 2*(size_-i) ];
        imagVecRef_neg[lid] = solBlock[ 2*(size_-i) + 1 ];
        realVecRef_pos[lid] = solBlock[ 2*i ];
        imagVecRef_pos[lid] = solBlock[ 2*i+1 ];
       
        }
      }
    }
  }

    // lead current vector
  Teuchos::RCP<Parallel::ParMap> baseLeadCurrentMap = Teuchos::rcp_const_cast<Parallel::ParMap>( hbBuilderPtr_->getBaseLeadCurrentMap() );
  Teuchos::RCP<Parallel::ParMap> globalLeadCurrentMap = Linear::createBlockParMap( blockCounttd, *baseLeadCurrentMap );
  Teuchos::RCP<Parallel::ParMap> globalLeadCurrentMapfq = Linear::createBlockParMap( blockCount, *baseLeadCurrentMap ); 

  freqDomainLeadCurrentVecReal  = Teuchos::rcp( Xyce::Linear::createBlockVector( blockCount, globalLeadCurrentMapfq, baseLeadCurrentMap ) );
  freqDomainLeadCurrentVecImaginary = Teuchos::rcp( Xyce::Linear::createBlockVector( blockCount, globalLeadCurrentMapfq, baseLeadCurrentMap ) );
  freqDomainJunctionVoltageVecReal = Teuchos::rcp( Xyce::Linear::createBlockVector( blockCount, globalLeadCurrentMapfq, baseLeadCurrentMap ) );
  freqDomainJunctionVoltageVecImaginary = Teuchos::rcp( Xyce::Linear::createBlockVector( blockCount, globalLeadCurrentMapfq, baseLeadCurrentMap ) );

  timeDomainLeadCurrentVec = Teuchos::rcp( Xyce::Linear::createBlockVector( blockCounttd, globalLeadCurrentMap, baseLeadCurrentMap ));
  timeDomainJunctionVoltageVec =  Teuchos::rcp( Xyce::Linear::createBlockVector( blockCounttd, globalLeadCurrentMap, baseLeadCurrentMap ) );

  N = timeDomainLeadCurrentVec->block(0).globalLength(); 

  for (int j=0; j<N; j++) 
  {
    // See if this time-domain solution variable is owned by the local processor.
    // If so, this processor owns the entire j-th block of the blockSolVecPtr vector,
    // and the j-th entry of every block in the freqDomainSolnVec[Real/Imag] vector.
    int lid = baseLeadCurrentMap->globalToLocalIndex( j );
    Linear::Vector& leadCurrentBlock =  bLeadCurrentVecFreqPtr_->block( j );

    Linear::Vector& realVecRef1 = freqDomainLeadCurrentVecReal->block((blockCount-1)/2);
    Linear::Vector& imagVecRef1 = freqDomainLeadCurrentVecImaginary->block((blockCount-1)/2);
    Linear::Vector& realVecRef2 = freqDomainJunctionVoltageVecReal->block((blockCount-1)/2);
    Linear::Vector& imagVecRef2 = freqDomainJunctionVoltageVecImaginary->block((blockCount-1)/2);

    if (lid >= 0)
    { 
      realVecRef1[lid] = leadCurrentBlock[0];  
      imagVecRef1[lid] = leadCurrentBlock[1];
      realVecRef2[lid] = leadCurrentBlock[0];  
      imagVecRef2[lid] = leadCurrentBlock[1];
    }

    for (int i=1; i <= (size_ -1)/2; ++i)
    {
      Linear::Vector& realVecRef_neg1 =  freqDomainLeadCurrentVecReal->block((blockCount-1)/2 - i);
      Linear::Vector& imagVecRef_neg1 =  freqDomainLeadCurrentVecImaginary->block((blockCount-1)/2 - i);
      Linear::Vector& realVecRef_pos2 =  freqDomainLeadCurrentVecReal->block((blockCount-1)/2 + i);
      Linear::Vector& imagVecRef_pos2 =  freqDomainLeadCurrentVecImaginary->block((blockCount-1)/2 + i);

      Linear::Vector& realVecRef_neg3 =  freqDomainJunctionVoltageVecReal->block((blockCount-1)/2 - i);
      Linear::Vector& imagVecRef_neg3 =  freqDomainJunctionVoltageVecImaginary->block((blockCount-1)/2 - i);
      Linear::Vector& realVecRef_pos4 =  freqDomainJunctionVoltageVecReal->block((blockCount-1)/2 + i);
      Linear::Vector& imagVecRef_pos4 =  freqDomainJunctionVoltageVecImaginary->block((blockCount-1)/2 + i);

      if (lid >= 0)
      {
        realVecRef_neg1[lid] = leadCurrentBlock[ 2*(size_ -i) ];
        imagVecRef_neg1[lid] = leadCurrentBlock[ 2*(size_ -i) + 1 ];
        realVecRef_pos2[lid] = leadCurrentBlock[ 2*i ];
        imagVecRef_pos2[lid] = leadCurrentBlock[ 2*i+1 ];
        realVecRef_neg3[lid] = leadCurrentBlock[ 2*(size_ -i) ];
        imagVecRef_neg3[lid] = leadCurrentBlock[ 2*(size_ -i) + 1 ]; 
        realVecRef_pos4[lid] = leadCurrentBlock[ 2*i ];
        imagVecRef_pos4[lid] = leadCurrentBlock[ 2*i+1 ];
      }
    }
  }

  if (bLeadCurrentVecFreqPtr_->blockCount() > 0 )
  {
    hbLoaderPtr_->permutedIFT(*bLeadCurrentVecFreqPtr_, &*timeDomainLeadCurrentVec, numTimePts_);
  }
 
  if (DEBUG_HB)
  {
    freqDomainSolnVecReal->print(std::cout);
    freqDomainSolnVecImag->print(std::cout);
  }
}


//-----------------------------------------------------------------------------
// Function      : HB::accumulateStatistics()
// Purpose       : Add in the statistics from the current analysis object
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist, 1355, Electrical Models & Simulation
// Creation Date : 05/29/13
//-----------------------------------------------------------------------------
void HB::accumulateStatistics_(AnalysisBase &analysis)
{
  hbStatCounts_ += analysis.stats_;
}


 
//-----------------------------------------------------------------------------
// Function      : HB::mapFreqs_
// Purpose       : Map frequency spectrum for HB analysis.
// Special Notes :
// Scope         : private
// Creator       : Ting Mei, SNL
// Creation Date : 04/09/2021
//-----------------------------------------------------------------------------
bool HB::mapFreqs_()
{
  int numAnalysisFreqs = freqs_.size();
  
  mappedFreqs_.resize(numAnalysisFreqs);

  mappedFreqs_[0] = 1.0;

  for (int i=1; i < numAnalysisFreqs; i++)
  {
    mappedFreqs_[i] = numFreqs_[i-1] * mappedFreqs_[i-1];

//    dout() << " mapped frequency point " << mappedFreqs_[i] << std::endl;
    
  }


  return true;
}


//-----------------------------------------------------------------------------
// Function      : HB::setFreqPoints_
// Purpose       : Set frequency spectrum for HB analysis.
// Special Notes :
// Scope         : private
// Creator       : Ting Mei, SNL
// Creation Date : 03/03/2014
//-----------------------------------------------------------------------------
bool HB::setFreqPointsAPFT_()
{
  if ( !intmodMaxGiven_)
  {
    int maxValue = 0;
    if (numFreqs_.size() != 0 )
      // find the max of numFreqs
    {
      maxValue = (numFreqs_[0] - 1)/2;
      for (int i=1; i<numFreqs_.size(); ++i)
      {
        if ((numFreqs_[i] - 1)/2 > maxValue)
          maxValue = (numFreqs_[i] - 1)/2 ;
      }
    } 
    else
    {
      maxValue = (size_ - 1)/2;
    }
    intmodMax_ = maxValue;
  }

  std::vector<int> k;
 
  int numAnalysisFreqs = freqs_.size();

  numPosFreqs.resize(numAnalysisFreqs);

  k.resize(numAnalysisFreqs);

  k[0] = 1;

  if (numFreqs_.size() != 0 )
    numPosFreqs[0] = (numFreqs_[0] - 1)/2;     
  else
    numPosFreqs[0]= (size_ - 1)/2;

  int numTotalFrequencies;

  int numExtraFreqs = 0;

  if (numPosFreqs[0] > intmodMax_) 
  {
    numExtraFreqs += ( numPosFreqs[0] - intmodMax_ );
    numFreqs_[0] = (intmodMax_*2 + 1);
  }

  numTotalFrequencies = numFreqs_[0];

  for (int i=1; i < numAnalysisFreqs; i++)
  {
    numPosFreqs[i] = (numFreqs_[i] - 1)/2;
   
    if (numPosFreqs[i] > intmodMax_)
    {
      numExtraFreqs += ( numPosFreqs[i] - intmodMax_ );
      numFreqs_[i] = (intmodMax_*2 + 1);
    }

    k[i] = k[i-1] * numFreqs_[i -1 ];

    numTotalFrequencies *= numFreqs_[i];  
  }

  if (DEBUG_HB) 
  {
    for (int i=0; i< numAnalysisFreqs; i++)
    {
      Xyce::dout() << "HB index " << i << std::endl;
      Xyce::dout() << "HB numPosFreqs =" << numPosFreqs[i] << std::endl;
      Xyce::dout() << "HB k =" << k[i] << std::endl;
    }
    Xyce::dout() << "HB numTotalFrequencies =" << numTotalFrequencies<< std::endl;
    Xyce::dout() << "HB numextrafreqs =" << numExtraFreqs << std::endl;
  }

  int numIndex = numTotalFrequencies;

  Teuchos::SerialDenseMatrix<int,double> indexMatrix(numAnalysisFreqs, numTotalFrequencies);

  if (DEBUG_HB)
  {
    Xyce::dout() << "HB intmodMax =" <<  intmodMax_ << std::endl;  
  }

  int nextIndex; 

  int idxMod,  idxValues;
  int sumIndex;

  std::vector<int> goodIndex;
  for (int i=0; i < numIndex; i++)      //  column
  {
    nextIndex = i;
    sumIndex = 0;

    for (int j= (numAnalysisFreqs - 1); j >= 0;  j-- )       // row 
    {
      idxMod = nextIndex%k[j];
      idxValues =  (nextIndex - idxMod)/k[j];

      indexMatrix (j, i) = static_cast<double>(idxValues - (numFreqs_[j] - 1)/2 );
      nextIndex = idxMod;
      sumIndex += abs(idxValues - (numFreqs_[j] - 1)/2 );

    }

    if( sumIndex <= intmodMax_) 
      goodIndex.push_back(i);
  }

  int diaindexSize = goodIndex.size();
 
  Teuchos::SerialDenseMatrix<int,double> diaindexMatrix( numAnalysisFreqs, (diaindexSize + numExtraFreqs) );
  diaindexMatrix.putScalar(0.0); 

  for (int i=0; i < diaindexSize; i++)       
  {
    for (int j= (numAnalysisFreqs - 1); j >= 0;  j-- )
      diaindexMatrix (j, i) = indexMatrix (j, goodIndex[i]);
  }

  if (DEBUG_HB)
  {
    for (int i=0; i< diaindexSize; i++)
    {
      dout() << "good index i = " << i << goodIndex[i]  << std::endl;
    }
    dout() <<  " checking diamond indexMatrix" << std::endl;
    diaindexMatrix.print(dout());

    dout() <<  " checking indexMatrix" << std::endl;
    indexMatrix.print(dout());
  }

  
  int extraIndexPos = diaindexSize;
  for (int i=0; i < numAnalysisFreqs ; i++)      //  column
  {

    if (numPosFreqs[i] > intmodMax_)
    {

      for (int j=0; j < (numPosFreqs[i] - intmodMax_ ) ; j++)
        diaindexMatrix (i, extraIndexPos + j ) = static_cast<double>(intmodMax_ + j + 1 );

      extraIndexPos += (numPosFreqs[i] - intmodMax_ );
    }  
  }


  if (DEBUG_HB)
  {
    dout() <<  " checking diamond indexMatrix after axis" << std::endl;
    diaindexMatrix.print(dout());
  }

// get the positive frequencies

  std::vector<double> posfreqPoints_;

  int posindexSize = (diaindexSize - 1)/2;
  posfreqPoints_.resize(posindexSize + numExtraFreqs); 

  Teuchos::SerialDenseMatrix<int,double> currindexMatrix( Teuchos::View, diaindexMatrix, numAnalysisFreqs, (posindexSize + numExtraFreqs), 0, posindexSize+1 );
  Teuchos::SerialDenseVector<int,double> currfreqPoints( Teuchos::View, &posfreqPoints_[0], (posindexSize + numExtraFreqs ) );

  Teuchos::SerialDenseVector<int,double> hbFreqs( Teuchos::View, &freqs_[0], numAnalysisFreqs);
//    Teuchos::SerialDenseVector<int,double> currWeightVector( Teuchos::View, &weightVector[i], oversampleRate*size_-(i+1) );
  currfreqPoints.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, currindexMatrix, hbFreqs, 0.0 );


  if (DEBUG_HB)
  {
    dout() << "checking positive frequencies" << std::endl;
    currindexMatrix.print(dout());
    hbFreqs.print(dout()); 
    currfreqPoints.print(dout());
  }

  for (int i=0; i < posindexSize; i++)       
  {
    if (posfreqPoints_[i] < 0.0)
      posfreqPoints_[i] = fabs( posfreqPoints_[i]);
  }

  std::sort(posfreqPoints_.begin(), posfreqPoints_.end() );


  if (DEBUG_HB) 
  {
    for (int i=0; i< posfreqPoints_.size(); i++)
      dout() << "pos frequency point " <<  posfreqPoints_[i] << std::endl;
  }

  posfreqPoints_.erase(std::unique(posfreqPoints_.begin(), posfreqPoints_.end() ), posfreqPoints_.end() );


  if (abs( posfreqPoints_[0]) < 2.0*Util::MachineDependentParams::MachinePrecision() )
    posfreqPoints_.erase( posfreqPoints_.begin()); 

  size_ = ( posfreqPoints_.size() ) *2 + 1;

  if (DEBUG_HB)
  {
    for (int i=0; i<posfreqPoints_.size(); i++)
    {
      dout() << "pos frequency point after " <<  posfreqPoints_[i] << std::endl;
    }

    Xyce::dout() << "HB size =" << size_ << std::endl;
  }

  int i=0;
  freqPoints_.resize(size_);

  for( i = 0; i < size_; ++i )
  {
    if (i < (size_-1)/2)
      freqPoints_[i] = - posfreqPoints_[ (size_-1)/2 - i - 1 ];
    else if (i > (size_-1)/2)
      freqPoints_[i] =  posfreqPoints_[ i - (size_-1)/2 - 1 ]; 
    else
      freqPoints_[i] = 0.0;
  } 

  if (DEBUG_HB)
  {
    for (int i=0; i<freqPoints_.size(); i++)
    {
      dout() << " frequency point " << freqPoints_[i] << std::endl;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::setFreqPoints_
// Purpose       : Set frequency spectrum for HB analysis.
// Special Notes :
// Scope         : private
// Creator       : Ting Mei, SNL
// Creation Date : 04/02/2021
//-----------------------------------------------------------------------------
bool HB::setFreqPoints_()
{

  if (  !selectHarmGiven_ )
  {
    if ( method_ == "APFT" )
    {

      selectHarm_ =  "HYBRID";
    }
    else if ( method_ == "AFM" )
    {
      selectHarm_ = "BOX";
    }
    else
    {
      Report::UserError() << "Unsupported method for HB";
      return false;
    }
  }


  if (  selectHarm_ == "HYBRID" )
  {

    setFreqPointsAPFT_();

    if ( method_ == "AFM" )
    {
      Report::UserError() << "Unsupported frequency truncation method for FM based HB";
      return false;
    }

  }
  else if (  selectHarm_  == "BOX" )
  {
    setFreqPointsFM_();
  }
  else if (  selectHarm_  == "DIAMOND" )
  {
    setFreqPointsDia_();
  }
  else
  {
    Report::UserError() << "Unsupported frequency truncation method for HB";
    return false;
  }

  return true;
}



//-----------------------------------------------------------------------------
// Function      : HB::setFreqPoints_
// Purpose       : Set frequency spectrum for HB analysis.
// Special Notes :
// Scope         : private
// Creator       : Ting Mei, SNL
// Creation Date : 04/02/2021
//-----------------------------------------------------------------------------
bool HB::setFreqPointsFM_()
{

  std::vector<int> k;

  int numAnalysisFreqs = freqs_.size();

//  numPosFreqs.resize(numAnalysisFreqs);

  k.resize(numAnalysisFreqs);

  k[0] = 1;

  int numTotalFrequencies;

  numTotalFrequencies = numFreqs_[0];

  for (int i=1; i < numAnalysisFreqs; i++)
  {
    k[i] = k[i-1] * numFreqs_[i -1 ];

    numTotalFrequencies *= numFreqs_[i];
  }

  if (DEBUG_HB)
  {
    for (int i=0; i< numAnalysisFreqs; i++)
    {
      Xyce::dout() << "HB index " << i << std::endl;
      Xyce::dout() << "HB k =" << k[i] << std::endl;
    }
    Xyce::dout() << "HB numTotalFrequencies =" << numTotalFrequencies<< std::endl;
  }

  int numIndex = numTotalFrequencies;

  Teuchos::SerialDenseMatrix<int,double> indexMatrix(numAnalysisFreqs, numTotalFrequencies);

  int nextIndex;

  int idxMod,  idxValues;

  for (int i=0; i < numIndex; i++)      //  column
  {
    nextIndex = i;

    for (int j= (numAnalysisFreqs - 1); j >= 0;  j-- )       // row
    {
      idxMod = nextIndex%k[j];
      idxValues =  (nextIndex - idxMod)/k[j];

      indexMatrix (j, i) = static_cast<double>(idxValues - (numFreqs_[j] - 1)/2 );
      nextIndex = idxMod;

    }

  }


  freqPoints_.resize(numTotalFrequencies);

  Teuchos::SerialDenseVector<int,double> currfreqPoints( Teuchos::View, &freqPoints_[0], numTotalFrequencies );

  Teuchos::SerialDenseVector<int,double> hbFreqs( Teuchos::View, &freqs_[0], numAnalysisFreqs);
//    Teuchos::SerialDenseVector<int,double> currWeightVector( Teuchos::View, &weightVector[i], oversampleRate*size_-(i+1) );
  currfreqPoints.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, indexMatrix, hbFreqs, 0.0 );


//  if (DEBUG_HB)
  {
    dout() << "checking frequencies" << std::endl;
    indexMatrix.print(dout());
    hbFreqs.print(dout());
    currfreqPoints.print(dout());
  }

  size_ = freqPoints_.size();

//  Xyce::dout() << "size = " << size_ << std::endl;

  return true;
}


//-----------------------------------------------------------------------------
// Function      : HB::setFreqPoints_
// Purpose       : Set frequency spectrum for HB analysis.
// Special Notes :
// Scope         : private
// Creator       : Ting Mei, SNL
// Creation Date : 06/02/2022
//-----------------------------------------------------------------------------
bool HB::setFreqPointsDia_()
{

  std::vector<int> k;

  int numAnalysisFreqs = freqs_.size();

//  numPosFreqs.resize(numAnalysisFreqs);

  k.resize(numAnalysisFreqs);

  numFreqs_[0] = (intmodMax_*2 + 1);

  k[0] = 1;

  int numTotalFrequencies;

  numTotalFrequencies = numFreqs_[0];

  for (int i=1; i < numAnalysisFreqs; i++)
  {
    k[i] = k[i-1] * numFreqs_[i -1 ];

    numFreqs_[i] = (intmodMax_*2 + 1);

    numTotalFrequencies *= numFreqs_[i];
  }

  if (DEBUG_HB)
  {
    for (int i=0; i< numAnalysisFreqs; i++)
    {
      Xyce::dout() << "HB index " << i << std::endl;
      Xyce::dout() << "HB k =" << k[i] << std::endl;
    }
    Xyce::dout() << "HB numTotalFrequencies =" << numTotalFrequencies<< std::endl;
  }

  int numIndex = numTotalFrequencies;

  Teuchos::SerialDenseMatrix<int,double> indexMatrix(numAnalysisFreqs, numTotalFrequencies);
//  if (DEBUG_HB)
  {
    Xyce::dout() << "HB intmodMax =" <<  intmodMax_ << std::endl;
  }

  int nextIndex;

  int idxMod,  idxValues, sumIndex ;   

  std::vector<int> goodIndex;

  for (int i=0; i < numIndex; i++)      //  column
  {
    nextIndex = i;
    sumIndex = 0;

    for (int j= (numAnalysisFreqs - 1); j >= 0;  j-- )       // row
    {
      idxMod = nextIndex%k[j];
      idxValues =  (nextIndex - idxMod)/k[j];

      indexMatrix (j, i) = static_cast<double>(idxValues - (numFreqs_[j] - 1)/2 );
      nextIndex = idxMod;

      sumIndex += abs(idxValues - (numFreqs_[j] - 1)/2 );
    }


    if( sumIndex <= intmodMax_)
      goodIndex.push_back(i);
  }

  int diaindexSize = goodIndex.size();

  Teuchos::SerialDenseMatrix<int,double> diaindexMatrix( numAnalysisFreqs, diaindexSize );

  diaindexMatrix.putScalar(0.0);

  for (int i=0; i < diaindexSize; i++)
  {
    for (int j= (numAnalysisFreqs - 1); j >= 0;  j-- )
      diaindexMatrix (j, i) = indexMatrix (j, goodIndex[i]);
  }

//  if (DEBUG_HB)
  {
    for (int i=0; i< diaindexSize; i++)
    {
      dout() << "good index i = " << i << goodIndex[i]  << std::endl;
    }
    dout() <<  " checking diamond indexMatrix" << std::endl;
    diaindexMatrix.print(dout());

    dout() <<  " checking indexMatrix" << std::endl;
    indexMatrix.print(dout());
  }

  freqPoints_.resize( diaindexSize );

  Teuchos::SerialDenseVector<int,double> currfreqPoints( Teuchos::View, &freqPoints_[0], diaindexSize );

  Teuchos::SerialDenseVector<int,double> hbFreqs( Teuchos::View, &freqs_[0], numAnalysisFreqs);
//    Teuchos::SerialDenseVector<int,double> currWeightVector( Teuchos::View, &weightVector[i], oversampleRate*size_-(i+1) );
  currfreqPoints.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, diaindexMatrix, hbFreqs, 0.0 );


//  if (DEBUG_HB)
  {
    dout() << "checking frequencies" << std::endl;
    diaindexMatrix.print(dout());
    hbFreqs.print(dout());
    currfreqPoints.print(dout());
  }

  size_ = freqPoints_.size();

  Xyce::dout() << "size = " << size_ << std::endl;

  return true;
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
HB::initializeOscOut( )
{
  // Set up the gid index for the osc out variable

  if ( hbOsc_)
  {

    bool oscOutError = false;

    std::string varName = refNode_;
//    Xyce::Util::toUpper(varName);

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

      topology_.getNodeSVarGIDs(NodeID(tmpVar, Xyce::_DNODE), oscoutList, dummyList, type1);
    }
    else if (pos2 !=(int)std::string::npos && pos3!=(int)std::string::npos) // this is a voltage variable
    {
      pos2+=2;
      int len = pos3-pos2;

      std::string tmpVar (varName,pos2,len);
      if (DEBUG_MPDE)
        Xyce::dout() << "tmpVar (for V-oscout) = " << tmpVar << std::endl;

      topology_.getNodeSVarGIDs(NodeID(tmpVar, Xyce::_VNODE), oscoutList, dummyList, type1);
    }
    else
    {
      // do nothing, assume that the varName from the netlist requires no
      // modification.
      topology_.getNodeSVarGIDs(NodeID(varName,-1), oscoutList, dummyList, type1);
    }

    if (oscoutList.size()==1)
    {
      refID_ =oscoutList.front();  // This is the GID of the OSCOUT.

//      if (refID_ >=0)
//        Xyce::dout() << " refNode ID = " << refID_ << std::endl;

    }
#ifndef Xyce_PARALLEL_MPI
    else  // This is not an error in parallel.
    {
      oscOutError = true;
    }
#else
    int tmpOscOut = refID_;
    pdsMgrPtr_->getPDSComm()->maxAll( &tmpOscOut, &refID_, 1 );

    if (refID_ < 0 )
    {
      oscOutError = true;
    }
#endif

    if (oscOutError)
    {
      Xyce::Report::UserWarning() << "Unrecognized value for HB option REFNode: "<< refNode_ ;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::setTimePoints_
// Purpose       : Set time points for multi-tone HB analysis.
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 03/05/2014
//-----------------------------------------------------------------------------
bool HB::setTimePoints_()
{
  // NOTE:  Need to make this parallel safe.

  // Resize fastTimes_ vector to be ready to receive from root processor. 
  fastTimes_.resize (size_); 

  int posFreq = (size_-1)/2;
  int oversampleRate = 5;
//  int periodSampleMultiplier = 3;
  int periodSampleMultiplier = 1;

  Teuchos::BLAS<int,double> blas;
  std::vector<double> testPoints(oversampleRate*size_);

  int myPID = pdsMgrPtr_->getPDSComm()->procID();

  if (myPID == 0)
  {

    for (int i=0; i<oversampleRate*size_; ++i)
    {
      testPoints[i] = periodSampleMultiplier*period_*((Teuchos::ScalarTraits<double>::random()+1)/2);
    }

    Teuchos::SerialDenseMatrix<int,double> testMatrix(size_,oversampleRate*size_);
    // NOTE: i represents frequency, j represents sample time
    // Set DC values first.
    for (int j=0; j<oversampleRate*size_; ++j)
    {
      testMatrix(0,j) = 1.0;
    }
    // Set rest of frequency values
    for (int i=1; i<=posFreq; i++)
    {
//    for (int j=0; j<oversampleRate*size_; j+2)
      for (int j=0; j<oversampleRate*size_; j++)
      {
        testMatrix(2*i-1,j) = cos(2*M_PI*freqPoints_[posFreq+i]*testPoints[j]);
        testMatrix(2*i,j) = sin(2*M_PI*freqPoints_[posFreq+i]*testPoints[j]);
      }
    }

    // Now for orthogonalization (yipeee)
    std::vector<double> weightVector(oversampleRate*size_);
    for (int i=0; i<size_; ++i)
    {
      // Find column with largest norm, choose as next vector for orthogonalization
      int maxIndex = 0;
      double maxValue = 0.0;
      for (int j=i; j<oversampleRate*size_; ++j)
      {
        Teuchos::SerialDenseMatrix<int,double> tempVector( Teuchos::View, testMatrix, size_, 1, 0, j );
        weightVector[j] = tempVector.normFrobenius();
    
        if (weightVector[j] > maxValue)
        {
          maxIndex = j;
          maxValue = weightVector[j];
        }
      }   
   

      // Swap time and vector.
      std::swap( testPoints[i], testPoints[maxIndex] );
      Teuchos::SerialDenseVector<int,double> newSwapVector2 = Teuchos::getCol<int,double>( Teuchos::Copy, testMatrix, maxIndex );
      Teuchos::SerialDenseVector<int,double> newSwapVector = Teuchos::getCol<int,double>( Teuchos::View, testMatrix, i );
      Teuchos::setCol<int,double>( newSwapVector, maxIndex, testMatrix );
      Teuchos::setCol<int,double>( newSwapVector2, i, testMatrix );

//    dout() << "Checking testMatrix for i= " << i << std::endl;
//    testMatrix.print(dout());   
 
      // Compute inner product with vector from last time point with rest of time points
      Teuchos::SerialDenseMatrix<int,double> currTestMatrix( Teuchos::View, testMatrix, size_, oversampleRate*size_-(i+1), 0, i+1 );
      Teuchos::SerialDenseVector<int,double> currWeightVector( Teuchos::View, &weightVector[i+1], oversampleRate*size_-(i+1) );
//    Teuchos::SerialDenseVector<int,double> currWeightVector( Teuchos::View, &weightVector[i], oversampleRate*size_-(i+1) );
      currWeightVector.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, currTestMatrix, newSwapVector, 0.0 );

//    dout() << "The current norm is" << std::endl;
//    currWeightVector.print(dout());

      // Subtract off scaled vector from rest of time points.
      // NOTE:  maxValue is the norm of the last time point.
      for (int j=i+1; j<oversampleRate*size_; ++j)
      {
         Teuchos::SerialDenseMatrix<int,double> currVector( Teuchos::View, testMatrix, size_, 1, 0, j );
         blas.AXPY( size_, -(currWeightVector[j-(i+1)]/(maxValue*maxValue) ), newSwapVector.values(), 1, currVector.values(), 1 );
      }


    }
//  dout() << "Checking testMatrix after orthogonalization." << std::endl;
//  testMatrix.print(dout());

    // Sort the chosen test points and then copy them into the fastTimes_ vector.
    std::sort( testPoints.begin(), testPoints.begin()+size_ );
  
    for (int i=0; i<size_; ++i)
    {
      fastTimes_[i] = testPoints[i];
    }
  }
  
  // Communicate testPoints to all processors (broadcast)
  pdsMgrPtr_->getPDSComm()->bcast( &fastTimes_[0], size_, 0 );


  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::createFT_
// Purpose       : Create the DFT and IFT matrices using the time points for multi-tone HB analysis.
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 03/05/2014
//-----------------------------------------------------------------------------
bool HB::createFT_()
{
  int posFreq = (size_-1)/2;
  idftMatrix_.reshape(size_,size_);

  // NOTE: i represents frequency, j represents sample time
  // Set DC values first.
  for (int i=0; i<size_; ++i)
  {
    idftMatrix_(i,0) = 1.0;
  }

  // Set rest of frequency values
  for (int i=0; i<size_; i++)
  {
    for (int j=1; j<=posFreq; j++)
    {
      idftMatrix_(i,2*j-1) = cos(2*M_PI*freqPoints_[posFreq+j]*fastTimes_[i]);
      idftMatrix_(i,2*j) = sin(2*M_PI*freqPoints_[posFreq+j]*fastTimes_[i]);
    }

  }

  // Compute DFT matrix.
  dftMatrix_ = idftMatrix_;
  Teuchos::SerialDenseSolver<int,double> ftSolver;
  ftSolver.setMatrix( Teuchos::rcp( &dftMatrix_, false ) );
  ftSolver.invert();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::setInitialGuess_
// Purpose       : Set initial guess for HB analysis.
// Special Notes :
// Scope         : private
// Creator       : Ting Mei, SNL
// Creation Date : 03/03/2014
//-----------------------------------------------------------------------------
bool HB::setInitialGuess()
{
  bool success = true;

  if (taHB_ == 1)
  {
    bool retTol1 = runTol();
    success = success && retTol1;

    // Start up periods need to be run before the initial condition is computed, otherwise
    // just used the solution from the tolerance calculation.
    if (startUpPeriodsGiven_)
    {
      bool startupPeriodsSuccess = runStartupPeriods();
      if (!startupPeriodsSuccess)
      {
        Report::UserError() << "Failed to calculate the startup periods";
        return false;
      }

      success = success && startupPeriodsSuccess;

      useStartupICs_ = true;
      bool icSuccess = runTransientIC();
      if (!icSuccess)
      {
        Report::UserError() << "Initial HB Transient failed";
        return false;
      }
      success = success && icSuccess;

      deviceManager_.setMPDEFlag( false );
    }

    interpolateIC(startUpPeriods_/freqs_[0]);
  }
  else if (taHB_ == 2 )
  {
    success = runDCOP();

  }

  return success;
}


 //-----------------------------------------------------------------------------
// Function      : HB::updateIFT_
// Purpose       : Create the IFT matrices using the time points for multi-tone HB analysis.
// Special Notes :
// Scope         : private
// Creator       : Ting Mei, SNL
// Creation Date : 07/2/2015
//-----------------------------------------------------------------------------
bool HB::updateIFT_( std::vector<double> & tPoints)
{  
  int posFreq = (size_-1)/2;
 //  idftMatrix_.reshape(size_,size_);

  int ntpt = tPoints.size();
  idftMatrix_.reshape( ntpt, size_);

  // NOTE: i represents frequency, j represents sample time
  // Set DC values first.
  for (int i=0; i<ntpt; ++i)
  {
    idftMatrix_(i,0) = 1.0;
  }

  // Set rest of frequency values
  for (int i=0; i< ntpt; i++)
  {
    for (int j=1; j<=posFreq; j++)
    {
      idftMatrix_(i,2*j-1) = cos(2*M_PI*freqPoints_[posFreq+j]* tPoints[i]);
      idftMatrix_(i,2*j) = sin(2*M_PI*freqPoints_[posFreq+j]*tPoints[i]);
    }
  } 


  return true; 
}

//-----------------------------------------------------------------------------
// Function      : HB::runTol
// Purpose       : Conducts transient run to determine right tolerance
//                 parameters for IC calculation
// Special Notes :
// Scope         : private
// Creator       : T. Mei, 1437, Electrical and Micro Modeling
// Creation Date : 02/23/09
//-----------------------------------------------------------------------------
bool
HB::runTol()
{
  Xyce::lout() << " ***** Computing tolerance parameters for HB IC calculation....\n" << std::endl;

  // Create a transient analysis object for this section.
  int numPoints = 0;
  {
    Xyce::IO::ActiveOutput x(analysisManager_.getOutputManagerAdapter().getOutputManager());
    if (!startUpPeriodsGiven_)
      x.add(Xyce::IO::PrintType::HB_IC, ANP_MODE_HB);

    isTransient_ = true;
    //Transient transient(analysisManager_, linearSystem_, nonlinearManager_, loader_, topology_, initialConditionsManager_, restartManager_, 0, this, 0);
    Transient transient(analysisManager_, &linearSystem_, nonlinearManager_, loader_, topology_, initialConditionsManager_, restartManager_, 0, this, 0);
    currentAnalysisObject_ = &transient;
    analysisManager_.pushActiveAnalysis(&transient);

    // Set options on transient object using options from .options timeint line
    transient.setTimeIntegratorOptions(saved_timeIntOB_); 

    // Modify TIAParams based on the running of one period.
    TimeIntg::TIAParams& new_tia_params = transient.getTIAParams();
    new_tia_params.initialTime = 0;
    new_tia_params.finalTime = 1.0/freqs_[0];
    analysisManager_.getStepErrorControl().pauseTime = new_tia_params.finalTime;
    new_tia_params.maxOrder = 1;
    transient.setAnalysisParams(Util::OptionBlock());
    transient.setSaveTimeSteps(!startUpPeriodsGiven_);
    transient.resetForHB();


    if (hbOsc_)
      transient.setNOOP(true);

    nonlinearManager_.resetAll(Nonlinear::DC_OP);
    analysisManager_.getStepErrorControl().resetAll(new_tia_params);
    analysisManager_.getDataStore()->resetAll(new_tia_params.absErrorTol, new_tia_params.relErrorTol);
    relErrorTol_ = new_tia_params.relErrorTol;
    analysisManager_.setNextOutputTime(0.0);

    {
      if (!transient.run())
      {
        Report::UserError() << "Calculation of tolerance parameters failed for relErrorTol = " << new_tia_params.relErrorTol;
        return false;
      }
    }

    numPoints = transient.getStepNumber();

    // Add in simulation times
    accumulateStatistics_(transient);

    analysisManager_.popActiveAnalysis();
    currentAnalysisObject_ = 0;
    isTransient_ = false;

    // If transient-assisted HB is used with .STEP then the data from each transient
    // run is copied to a tmp file.  This block copies that tmp data to the actual
    // HB_IC output file, if the transient run was acceptable.  See SON Bug 928
    // for more details.
    if (!startUpPeriodsGiven_ && acceptICTranRun(numPoints))
      x.copyTmpFileToOutputFile();
  }

  while ( !acceptICTranRun(numPoints) )
  {
    {
      Report::UserWarning() << "Tolerance parameters refined, re-running with relErrorTol = " << relErrorTol_/10;
    }

    if (!startUpPeriodsGiven_)
    {
      // Clear the fast time data storage before performing the next transient
      analysisManager_.getDataStore()->resetFastTimeData();
    }

    relErrorTol_ = relErrorTol_/10;

    // Create a transient analysis object for this section.
    {
      Xyce::IO::ActiveOutput x(analysisManager_.getOutputManagerAdapter().getOutputManager());
      if (!startUpPeriodsGiven_) {
        x.add(Xyce::IO::PrintType::HB_IC, ANP_MODE_HB);
        x.resetIndex();
        x.reopenTmpFile();  // tmp file used to handle HB_IC data from each transient run
      }

      isTransient_ = true;
      //Transient transient(analysisManager_, linearSystem_, nonlinearManager_, loader_, topology_, initialConditionsManager_, restartManager_, 0, this, 0);
      Transient transient(analysisManager_, &linearSystem_, nonlinearManager_, loader_, topology_, initialConditionsManager_, restartManager_, 0, this, 0);

      if (hbOsc_)
        transient.setNOOP(true);

      currentAnalysisObject_ = &transient;
      analysisManager_.pushActiveAnalysis(&transient);
      transient.setSaveTimeSteps(!startUpPeriodsGiven_);

      // Set options on transient object using options from .options timeint line
      transient.setTimeIntegratorOptions(saved_timeIntOB_);

      // Modify TIAParams based on the running of one period.
      TimeIntg::TIAParams& new_tia_params = transient.getTIAParams();
      new_tia_params.initialTime = 0;
      new_tia_params.finalTime = 1.0/freqs_[0];
      analysisManager_.getStepErrorControl().pauseTime = new_tia_params.finalTime;
      new_tia_params.maxOrder = 1;
      new_tia_params.relErrorTol = relErrorTol_;
      transient.setAnalysisParams(Util::OptionBlock());
      transient.resetForHB();
      nonlinearManager_.resetAll(Nonlinear::DC_OP);
      analysisManager_.getStepErrorControl().resetAll(new_tia_params);
      analysisManager_.getDataStore()->resetAll(new_tia_params.absErrorTol, new_tia_params.relErrorTol);
      analysisManager_.setNextOutputTime(0.0);

      {
        if (!transient.run())
        {
          Report::UserError() << "Calculation of tolerance parameters failed for relErrorTol = " << new_tia_params.relErrorTol;
          return false;
        }
      }

      numPoints = transient.getStepNumber();

      // Add in simulation times
      accumulateStatistics_(transient);

      analysisManager_.popActiveAnalysis();
      currentAnalysisObject_ = 0;
      isTransient_ = false;
  
      // copies tmp data for this transient run to the actual HB_IC output file,
      // if the run was acceptable.
      if ( !startUpPeriodsGiven_ && acceptICTranRun(numPoints) )
        x.copyTmpFileToOutputFile();
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::acceptICTranRun()
// Purpose       : Determine if a transient run, done for the HB_IC calculation,
//                 was acceptable.
// Special Notes : This function is also used to trigger the copying of the tmp
//                 file to the actual HB_IC output file.
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 11/19/2019
//-----------------------------------------------------------------------------
bool HB::acceptICTranRun(int numPoints)
{
  return !(numPoints < (1.2*size_)) || !(relErrorTol_ >= 1e-6);
}

//-----------------------------------------------------------------------------
// Function      : HB::runStartupPeriods_()
// Purpose       : Runs normal transient problem through the requested
//                 number of startup periods
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, 1437, Electrical and Micro Modeling
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool
HB::runStartupPeriods()
{
  bool returnValue = true;

  Xyce::lout() << "  ***** Computing " << startUpPeriods_ << " start up periods for HB IC calculation...." << std::endl;

  // << "HB::runStartupPeriods_():  Double DCOP tia_params:"
  // << " doubleDCOPStep = " << tia_params.doubleDCOPStep
  // << " firstDCOPStep = " << tia_params.firstDCOPStep
  // << " lastDCOPStep = " << tia_params.lastDCOPStep << std::endl

  // << "HB::runStartupPeriods_():  Double DCOP tiaParamsSave:"
  // << " doubleDCOPStep = " << tiaParamsSave.doubleDCOPStep
  // << " firstDCOPStep = " << tiaParamsSave.firstDCOPStep
  // << " lastDCOPStep = " << tiaParamsSave.lastDCOPStep << std::endl;

  {
    Xyce::IO::ActiveOutput x(analysisManager_.getOutputManagerAdapter().getOutputManager());
    x.add(Xyce::IO::PrintType::HB_STARTUP, ANP_MODE_HB);

    // Create a transient analysis object for this section.
    {
      isTransient_ = true;
      //Transient transient(analysisManager_, linearSystem_, nonlinearManager_, loader_, topology_, initialConditionsManager_, restartManager_, 0, this, 0);
      Transient transient(analysisManager_, &linearSystem_, nonlinearManager_, loader_, topology_, initialConditionsManager_, restartManager_, 0, this, 0);

      if (hbOsc_)
        transient.setNOOP(true);

      analysisManager_.pushActiveAnalysis(&transient);
      currentAnalysisObject_ = &transient;

      // Set options on transient object using options from .options timeint line
      transient.setTimeIntegratorOptions(saved_timeIntOB_);

      // Modify TIAParams based on the running of one period.
      TimeIntg::TIAParams& new_tia_params = transient.getTIAParams();
      new_tia_params.initialTime = 0.0;
      new_tia_params.finalTime = startUpPeriods_/freqs_[0];
      new_tia_params.relErrorTol = relErrorTol_;
      analysisManager_.getStepErrorControl().pauseTime = new_tia_params.finalTime;
      transient.setAnalysisParams(Util::OptionBlock());
      transient.resetForHB();
      nonlinearManager_.resetAll(Nonlinear::DC_OP);
      analysisManager_.getStepErrorControl().resetAll(new_tia_params);
      analysisManager_.getDataStore()->resetAll(new_tia_params.absErrorTol, new_tia_params.relErrorTol);
      analysisManager_.setNextOutputTime(0.0);

      if (DEBUG_HB)
        Xyce::dout() << "HB::runStartupPeriods_():  Advancing time through "
                     << startUpPeriods_ << " startup periods"
                     << " initialTime = " << new_tia_params.initialTime
                     << " finalTime = " << new_tia_params.finalTime << std::endl;
  
      {
        returnValue = transient.run();
      }

      isTransient_ = false;

      // Add in simulation times
      accumulateStatistics_(transient);

      analysisManager_.getOutputManagerAdapter().finishOutput();

      analysisManager_.popActiveAnalysis();
      currentAnalysisObject_ = 0;
    }
  }

  // reset the output filename suffix
  // analysisManager_.outMgrPtr->setOutputFilenameSuffix("");

  // put the dsPtr->currentSolutionPtr into dcOpSol and State Vec so that it
  // is used as our initial condition for the pending fast time scale runs
  TimeIntg::DataStore * dsPtr = analysisManager_.getDataStore();
  dcOpSolVecPtr_ = rcp( dsPtr->currSolutionPtr->cloneCopyVector() );
  dcOpStateVecPtr_ = rcp( dsPtr->currStatePtr->cloneCopyVector() ); 
  dcOpQVecPtr_ = rcp( dsPtr->daeQVectorPtr->cloneCopyVector() );
  dcOpStoreVecPtr_ = rcp( dsPtr->currStorePtr->cloneCopyVector() );

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : HB::runTransientIC_
// Purpose       : Conducts a regular transient run for HB initial conditions
// Special Notes :
// Scope         : private
// Creator       : Ting Mei, SNL
// Creation Date : 10/03/2008
//-----------------------------------------------------------------------------
bool
HB::runTransientIC()
{
  bool returnValue = true;

  Xyce::lout() << " ***** Running transient to compute HB initial condition....\n" << std::endl;

  // this prevents extra DC op data from being printed.
  deviceManager_.setMPDEFlag( true );

  // Initial conditions will be set if startup periods were run.
  TimeIntg::DataStore * dsPtr = analysisManager_.getDataStore();
  *(dsPtr->nextSolutionPtr) = *(dcOpSolVecPtr_);
  *(dsPtr->nextStatePtr) = *(dcOpStateVecPtr_);
  *(dsPtr->daeQVectorPtr) = *(dcOpQVecPtr_);
  *(dsPtr->nextStorePtr) = *(dcOpStoreVecPtr_);

  // Create a transient analysis object for this section.
  {
    Xyce::IO::ActiveOutput x(analysisManager_.getOutputManagerAdapter().getOutputManager());
    if (saveIcData_ && startUpPeriodsGiven_)
      x.add(Xyce::IO::PrintType::HB_IC, ANP_MODE_HB);

    isTransient_ = true;
    //Transient transient(analysisManager_, linearSystem_, nonlinearManager_, loader_, topology_, initialConditionsManager_, restartManager_, 0, this, 0);
    Transient transient(analysisManager_, &linearSystem_, nonlinearManager_, loader_, topology_, initialConditionsManager_, restartManager_, 0, this, 0);
    transient.setNOOP(true);
    transient.setSaveTimeSteps(true);
    currentAnalysisObject_ = &transient;
    analysisManager_.pushActiveAnalysis(&transient);

    // Set options on transient object using options from .options timeint line
    transient.setTimeIntegratorOptions(saved_timeIntOB_);

    // Modify TIAParams based on the running of one period.
    TimeIntg::TIAParams& new_tia_params = transient.getTIAParams();
    new_tia_params.initialTime = startUpPeriods_/freqs_[0]; ;
    new_tia_params.finalTime = new_tia_params.initialTime + 1.0/freqs_[0];
    new_tia_params.relErrorTol = relErrorTol_;
    analysisManager_.getStepErrorControl().pauseTime = new_tia_params.finalTime;
    new_tia_params.maxOrder = 1;
    transient.setAnalysisParams(Util::OptionBlock());
    transient.resetForHB();
    nonlinearManager_.resetAll(Nonlinear::DC_OP);
    analysisManager_.getStepErrorControl().resetAll(new_tia_params);
    analysisManager_.getDataStore()->resetAll(new_tia_params.absErrorTol, new_tia_params.relErrorTol);
    analysisManager_.setNextOutputTime(0.0);

    if (DEBUG_HB)
      Xyce::dout() << "HB::runTransientIC_():  Advancing time from"
                   << " initialTime = " << new_tia_params.initialTime
                   << " finalTime = " << new_tia_params.finalTime << std::endl;

    {
      returnValue = transient.run();
    }

    isTransient_ = false;

    // Add in simulation times
    accumulateStatistics_(transient);

    analysisManager_.popActiveAnalysis();
    currentAnalysisObject_ = 0;
  }

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : HB::runDCop()
// Purpose       : Runs the dc op problem
// Special Notes :
// Scope         : private
// Creator       : 
// Creation Date : 04/18
//-----------------------------------------------------------------------------
bool
HB::runDCOP()
{
  bool returnValue = true;

    DCSweep dc_sweep(analysisManager_, &linearSystem_, nonlinearManager_, loader_, topology_, initialConditionsManager_, 0);
    currentAnalysisObject_ = &dc_sweep;
    analysisManager_.pushActiveAnalysis(&dc_sweep);
    returnValue = dc_sweep.run();

    analysisManager_.popActiveAnalysis();

    currentAnalysisObject_ = 0;

  dcOpSolVecPtr_ =  rcp( analysisManager_.getDataStore()->currSolutionPtr->cloneCopyVector() );
  dcOpStateVecPtr_ = rcp( analysisManager_.getDataStore()->currStatePtr->cloneCopyVector() );
  dcOpQVecPtr_ =  rcp( analysisManager_.getDataStore()->daeQVectorPtr->cloneCopyVector() );
  dcOpStoreVecPtr_ = rcp( analysisManager_.getDataStore()->currStorePtr->cloneCopyVector() );

  return returnValue;


}
         
//-----------------------------------------------------------------------------
// Function      : HB::interpolateIC()
// Purpose       : Tries to filter the fast time points from a transient run
//                 so that points are not too close together
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, 1437, Electrical and Micro Modeling
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool
HB::interpolateIC(
  double                initial_time)
{
  Xyce::lout() << " ***** Interpolating transient solution for IC calculation....\n" << std::endl;

  TimeIntg::DataStore *dsPtr = analysisManager_.getDataStore();
  int numPoints = dsPtr->timeSteps.size();

  if  (DEBUG_HB)
    Xyce::dout() << "HB::interpolateIC_(): Initial transient run produced " << numPoints << " points." << std::endl;

  std::vector<int> goodIndicies;
  for( int i = 0; i < size_; ++i )
  {
    goodTimePoints_[i] = initial_time + goodTimePoints_[i];
  }

  bool sortedTimesFlag = true;

  if (freqs_.size() > 1)
    sortedTimesFlag = false;

  int breakpoints = 0;          // need to keep track of how many breakpoints there are
  int startIndex = 0;

  if (sortedTimesFlag)
  {
    goodIndicies.push_back(startIndex);
    int GoodTimePointIndex = startIndex + 1;

    for( int i=startIndex; i < numPoints - 1 ; i++ )
    {
    // count up breakpoints

      if( dsPtr->timeStepsBreakpointFlag[i] == true )
      {
        breakpoints++;
      }

      if (DEBUG_HB && isActive(Diag::HB_TIMESTEP))
      {
        Xyce::dout() << "\t\t timeStep[ " << i << " ] = " << dsPtr->timeSteps[i];
        if( dsPtr->timeStepsBreakpointFlag[i] == true )
        {
          Xyce::dout() << "  Breakpoint";
        }
        Xyce::dout() << std::endl;
      }

      while( ( GoodTimePointIndex < size_ )  && (dsPtr->timeSteps[i] <= goodTimePoints_[GoodTimePointIndex]) && (goodTimePoints_[GoodTimePointIndex] < dsPtr->timeSteps[i+1]))
      {
      // found a good point so save the index
        goodIndicies.push_back( i );
        GoodTimePointIndex =  GoodTimePointIndex+1;
      }
    }
  }
  else 
  {
    double tmpTimePoints;
    bool found;
    for(int i=0; i<size_; i++ )
    {
      tmpTimePoints = fmod(fastTimes_[i],  1.0/freqs_[0]) + initial_time;

      goodTimePoints_[i] = tmpTimePoints;
      found = false;

      for(int j=0; (j< (numPoints - 1) && !found ); j++)
      {
        if((dsPtr->timeSteps[j] <= tmpTimePoints ) && ( tmpTimePoints < dsPtr->timeSteps[j+1]))
        {
          // found a good point so save the index
          goodIndicies.push_back( j );
          found=true;
        }
      }
    }
  }

  for(int i=0; i<size_; i++ )
  {
    int currentIndex = goodIndicies[i];
    Linear::Vector * firstSolVecPtr = dsPtr->fastTimeSolutionVec[currentIndex];
    Linear::Vector * secondSolVecPtr = dsPtr->fastTimeSolutionVec[currentIndex+1];

    Linear::Vector * firstStateVecPtr = dsPtr->fastTimeStateVec[currentIndex];
    Linear::Vector * secondStateVecPtr = dsPtr->fastTimeStateVec[currentIndex+1];

    Linear::Vector * firstQVecPtr = dsPtr->fastTimeQVec[currentIndex];
    Linear::Vector * secondQVecPtr = dsPtr->fastTimeQVec[currentIndex+1];

    Linear::Vector * firstStoreVecPtr = dsPtr->fastTimeStoreVec[currentIndex];
    Linear::Vector * secondStoreVecPtr = dsPtr->fastTimeStoreVec[currentIndex+1];

    double fraction = (goodTimePoints_[i] -  dsPtr->timeSteps[currentIndex])/(dsPtr->timeSteps[currentIndex+1] -  dsPtr->timeSteps[currentIndex]);

    RCP<Linear::Vector> InterpICSolVecPtr = rcp( secondSolVecPtr->cloneCopyVector() );
    RCP<Linear::Vector> InterpICStateVecPtr = rcp( secondStateVecPtr->cloneCopyVector() );
    RCP<Linear::Vector> InterpICQVecPtr = rcp( secondQVecPtr->cloneCopyVector() );
    RCP<Linear::Vector> InterpICStoreVecPtr = rcp( secondStoreVecPtr->cloneCopyVector() );

    InterpICSolVecPtr->update(-1.0, *firstSolVecPtr, 1.0);
    InterpICSolVecPtr->update(1.0, *firstSolVecPtr, fraction);

    InterpICStateVecPtr->update(-1.0, *firstStateVecPtr, 1.0);
    InterpICStateVecPtr->update(1.0, *firstStateVecPtr, fraction);

    InterpICQVecPtr->update(-1.0, *firstQVecPtr, 1.0);
    InterpICQVecPtr->update(1.0, *firstQVecPtr, fraction);

    InterpICStoreVecPtr->update(-1.0, *firstStoreVecPtr, 1.0);
    InterpICStoreVecPtr->update(1.0, *firstStoreVecPtr, fraction);

    goodSolutionVec_.push_back(InterpICSolVecPtr);
    goodStateVec_.push_back(InterpICStateVecPtr);
    goodQVec_.push_back(InterpICQVecPtr);
    goodStoreVec_.push_back(InterpICStoreVecPtr);
  }

  // Clean up the fast time data since we are finished computing the initial condition.
  // The fast time data can take a considerable amount of memory for large problems.
  dsPtr->resetFastTimeData();

  return true;
}

namespace {

typedef Util::Factory<AnalysisBase, HB>  HBFactoryBase;


//-----------------------------------------------------------------------------
// Class         : HBFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jan 29 12:53:02 2015
//-----------------------------------------------------------------------------
///
/// Factory for parsing HB parameters from the netlist and creating HB analysis.
///
class HBFactory : public HBFactoryBase
{
public:
  //-----------------------------------------------------------------------------
  // Function      : HBFactory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:54:09 2015
  //-----------------------------------------------------------------------------
  ///
  /// Constructs the HB analysis factory
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
  /// @param device_manager 
  /// @param builder 
  /// @param topology 
  ///
  HBFactory(
    Analysis::AnalysisManager & analysis_manager,
    Linear::System &                    linear_system,
    Nonlinear::Manager &                nonlinear_manager,
    Loader::Loader &                    loader,
    Device::DeviceMgr &                 device_manager,
    Linear::Builder &                   builder,
    Topo::Topology &                    topology,
    IO::InitialConditionsManager &      initial_conditions_manager,
    IO::RestartMgr &                    restart_manager)
    : HBFactoryBase(),
      analysisManager_(analysis_manager),
      linearSystem_(linear_system),
      nonlinearManager_(nonlinear_manager),
      loader_(loader),
      deviceManager_(device_manager),
      builder_(builder),
      topology_(topology),
      initialConditionsManager_(initial_conditions_manager),
      restartManager_(restart_manager)
  {}

  virtual ~HBFactory()
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
  /// Create a new HB analysis and applies the analysis and time integrator option blocks.
  ///
  /// @return new HB analysis object
  ///
  HB *create() const
  {
    analysisManager_.setAnalysisMode(ANP_MODE_HB);

    HB *hb = new HB(analysisManager_, linearSystem_, nonlinearManager_, loader_, deviceManager_, builder_, topology_, initialConditionsManager_, restartManager_);
    hb->setAnalysisParams(hbAnalysisOptionBlock_);
    hb->setHBIntParams(hbIntOptionBlock_);
    hb->setHBLinSol(hbLinSolOptionBlock_, builder_);
    hb->setLinSol(linSolOptionBlock_);
    hb->setTimeInt(timeIntOptionBlock_);

    return hb;
  }

  //-----------------------------------------------------------------------------
  // Function      : setHBAnalysisOptionBlock
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
  bool setHBAnalysisOptionBlock(const Util::OptionBlock &option_block)
  {
    hbAnalysisOptionBlock_ = option_block;

    return true;
  }

  //-----------------------------------------------------------------------------
  // Function      : setHBIntOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 13:00:14 2015
  //-----------------------------------------------------------------------------
  ///
  /// Saves the HBINT parsed options block in the factory.
  ///
  /// @invariant Overwrites any previously specified HBINT option block.
  ///
  /// @param option_block parsed option block
  ///
  bool setHBIntOptionBlock(const Util::OptionBlock &option_block)
  {
    hbIntOptionBlock_ = option_block;

    return true;
  }

  //-----------------------------------------------------------------------------
  // Function      : setHBIntOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 13:00:14 2015
  //-----------------------------------------------------------------------------
  ///
  /// Saves the LINSOL-HB parsed options block in the factory.
  ///
  /// @invariant Overwrites any previously specified LINSOL-HB option block.
  ///
  /// @param option_block parsed option block
  ///
  bool setHBLinSolOptionBlock(const Util::OptionBlock &option_block)
  {
    hbLinSolOptionBlock_ = option_block;

    return true;
  }

  //-----------------------------------------------------------------------------
  // Function      : setLinSolOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 13:00:14 2015
  //-----------------------------------------------------------------------------
  ///
  /// Saves the LINSOL parsed options block in the factory.
  ///
  /// @invariant Overwrites any previously specified LINSOL option block.
  ///
  /// @param option_block parsed option block
  ///
  bool setLinSolOptionBlock(const Util::OptionBlock &option_block)
  {
    linSolOptionBlock_ = option_block;

    return true;
  }

  //-----------------------------------------------------------------------------
  // Function      : setTimeIntegratorOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 13:00:14 2015
  //-----------------------------------------------------------------------------
  ///
  /// Saves the TIMEINT parsed options block in the factory.
  ///
  /// @invariant Overwrites any previously specified TIMEINT option block.
  ///
  /// @param option_block parsed option block
  ///
  bool setTimeIntegratorOptionBlock(const Util::OptionBlock &option_block)
  {
    timeIntOptionBlock_ = option_block;

    return true;
  }

public:
  AnalysisManager &                     analysisManager_;
  Linear::System &                      linearSystem_;
  Nonlinear::Manager &                  nonlinearManager_;
  Loader::Loader &                      loader_;
  Device::DeviceMgr &                   deviceManager_;
  Linear::Builder &                     builder_;
  Topo::Topology &                      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;
  IO::RestartMgr &                      restartManager_;

private:
  Util::OptionBlock     hbAnalysisOptionBlock_;
  Util::OptionBlock     timeIntegratorOptionBlock_;
  Util::OptionBlock     hbIntOptionBlock_;
  Util::OptionBlock     hbLinSolOptionBlock_;
  Util::OptionBlock     linSolOptionBlock_;
  Util::OptionBlock     timeIntOptionBlock_;
};

// .HB
struct HBAnalysisReg : public IO::PkgOptionsReg
{
  HBAnalysisReg(
    HBFactory &   factory )
    : factory_(factory)
  {}

  bool operator()(const Util::OptionBlock &option_block)
  {
    factory_.setHBAnalysisOptionBlock(option_block);
    factory_.deviceManager_.setBlockAnalysisFlag(true);

    factory_.analysisManager_.addAnalysis(&factory_);

    return true;
  }

  HBFactory &         factory_;
};

//-----------------------------------------------------------------------------
// Function      : extractHBData
// Purpose       : Extract the parameters from a netlist .DC line held in
//                 parsed_line.
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool
extractHBData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("HB", Util::OptionBlock::NO_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

  int numFields = parsed_line.size();

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 2 )
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << ".HB line has an unexpected number of fields";
  }

  int linePosition = 1;   // Start of parameters on .param line.
  int endPosition = numFields;

  Util::Param parameter("", "");

// frequency of oscillation is required
  std::vector<double> freqs(numFields - 1);

  int i = 0;
  while( linePosition < endPosition )
  {
    const std::string & value = parsed_line[linePosition].string_;
    if (Util::isValue(value))
    {
      freqs[i] = Util::Value(value);
    }
    else
    {
      Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
        << "Attempt to assign value for FREQ from " << value;
    }
    ++linePosition;
    ++i;
  }

  parameter.setTag( "FREQ" );
  parameter.setVal( freqs );
  option_block.addParam( parameter );
//  ++linePosition;

  circuit_block.addOptions(option_block);

  return true;
}

void
populateMetadata(
  IO::PkgOptionsMgr &   options_manager)
{
  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("LINSOL-HB");

    parameters.insert(Util::ParamMap::value_type("AZ_max_iter", Util::Param("AZ_max_iter", 200)));
    parameters.insert(Util::ParamMap::value_type("AZ_solver", Util::Param("AZ_solver", 1)));
    parameters.insert(Util::ParamMap::value_type("AZ_conv", Util::Param("AZ_conv", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_pre_calc", Util::Param("AZ_pre_calc", 1)));
    parameters.insert(Util::ParamMap::value_type("AZ_keep_info", Util::Param("AZ_keep_info", 1)));
    parameters.insert(Util::ParamMap::value_type("AZ_orthog", Util::Param("AZ_orthog", 1)));
    parameters.insert(Util::ParamMap::value_type("AZ_reorder", Util::Param("AZ_reorder", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_scaling", Util::Param("AZ_scaling", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_kspace", Util::Param("AZ_kspace", 50)));
    parameters.insert(Util::ParamMap::value_type("AZ_tol", Util::Param("AZ_tol", 1.0E-9)));
    parameters.insert(Util::ParamMap::value_type("AZ_output", Util::Param("AZ_output", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_diagnostics", Util::Param("AZ_diagnostics", 0)));
    parameters.insert(Util::ParamMap::value_type("TYPE", Util::Param("TYPE", "DEFAULT")));
    parameters.insert(Util::ParamMap::value_type("DIRECT_SOLVER", Util::Param("DIRECT_SOLVER", "DEFAULT")));
    parameters.insert(Util::ParamMap::value_type("PREC_TYPE", Util::Param("PREC_TYPE", "BLOCK_JACOBI")));
    parameters.insert(Util::ParamMap::value_type("BELOS_SOLVER_TYPE", Util::Param("BELOS_SOLVER_TYPE", "Block GMRES")));
    parameters.insert(Util::ParamMap::value_type("BLOCK_JACOBI_CORRECTED", Util::Param("BLOCK_JACOBI_CORRECTED", false)));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_LS", Util::Param("OUTPUT_LS", 1)));
  }

  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("HBINT");

    parameters.insert(Util::ParamMap::value_type("TEST", Util::Param("TEST", false)));
    parameters.insert(Util::ParamMap::value_type("NUMFREQ", Util::Param("NUMFREQ", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("STARTUPPERIODS", Util::Param("STARTUPPERIODS", 0)));
    parameters.insert(Util::ParamMap::value_type("SAVEICDATA", Util::Param("SAVEICDATA", false)));
    parameters.insert(Util::ParamMap::value_type("DEBUGLEVEL", Util::Param("DEBUGLEVEL", 0)));
    parameters.insert(Util::ParamMap::value_type("TAHB", Util::Param("TAHB", 1)));
    parameters.insert(Util::ParamMap::value_type("VOLTLIM", Util::Param("VOLTLIM", 1)));
    parameters.insert(Util::ParamMap::value_type("INTMODMAX", Util::Param("INTMODMAX",0)));
    parameters.insert(Util::ParamMap::value_type("METHOD", Util::Param("METHOD", "APFT")));
    parameters.insert(Util::ParamMap::value_type("NUMTPTS", Util::Param("NUMTPTS", 1)));

    parameters.insert(Util::ParamMap::value_type("HBOSC", Util::Param("HBOSC", false)));
    parameters.insert(Util::ParamMap::value_type("LOADTIMESOURCES", Util::Param("LOADTIMESOURCES", 1)));
    parameters.insert(Util::ParamMap::value_type("SELECTHARMS", Util::Param("SELECTHARMS", "BOX")));
    parameters.insert(Util::ParamMap::value_type("REFNODE", Util::Param("REFNODE",  "")));
  }
}

} // namespace <unnamed>


bool
registerHBFactory(
  FactoryBlock &        factory_block)
{
  HBFactory *factory = new HBFactory(factory_block.analysisManager_, factory_block.linearSystem_, factory_block.nonlinearManager_,
                                     factory_block.loader_, factory_block.deviceManager_, factory_block.builder_, factory_block.topology_,
                                     factory_block.initialConditionsManager_, factory_block.restartManager_);

  addAnalysisFactory(factory_block, factory);

  populateMetadata(factory_block.optionsManager_);

  factory_block.optionsManager_.addCommandParser(".HB", extractHBData);

  factory_block.optionsManager_.addCommandProcessor("HB", new HBAnalysisReg(*factory));

  factory_block.optionsManager_.addOptionsProcessor("HBINT", IO::createRegistrationOptions(*factory, &HBFactory::setHBIntOptionBlock));
  factory_block.optionsManager_.addOptionsProcessor("LINSOL-HB", IO::createRegistrationOptions(*factory, &HBFactory::setHBLinSolOptionBlock));
  factory_block.optionsManager_.addOptionsProcessor("LINSOL", IO::createRegistrationOptions(*factory, &HBFactory::setLinSolOptionBlock));
  factory_block.optionsManager_.addOptionsProcessor("TIMEINT", IO::createRegistrationOptions(*factory, &HBFactory::setTimeIntegratorOptionBlock));
 
  return true;
}

} // namespace Analysis
} // namespace Xyce

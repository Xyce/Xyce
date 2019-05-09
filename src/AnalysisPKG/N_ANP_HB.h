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

//-----------------------------------------------------------------------------
//
// Purpose        : HB analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Todd Coffey, 1414, Ting Mei 1437
//
// Creation Date  : 07/23/08
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_HB_h
#define Xyce_N_ANP_HB_h

#include <vector>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_RCP.hpp>

#include <N_ANP_fwd.h>
#include <N_LOA_fwd.h>
#include <N_LAS_fwd.h>
#include <N_TOP_fwd.h>

#include <N_ANP_AnalysisBase.h>
#include <N_ANP_StepEvent.h>
#include <N_ANP_RegisterAnalysis.h>
#include <N_UTL_DFTInterfaceDecl.hpp>
#include <N_UTL_FFTInterface.hpp>
#include <N_UTL_Listener.h>
#include <N_UTL_OptionBlock.h>

namespace Xyce {
namespace Analysis {

//-------------------------------------------------------------------------
// Class         : HB
// Purpose       : HB analysis class
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class HB : public AnalysisBase, public Util::ListenerAutoSubscribe<StepEvent>
{
public:
  HB(
    AnalysisManager &                   analysis_manager,
    Linear::System &                    linear_system,
    Nonlinear::Manager &                nonlinear_manager,
    Loader::Loader &                    loader,
    Device::DeviceMgr &                 device_manager,
    Linear::Builder &                   builder,
    Topo::Topology &                    topology,
    IO::InitialConditionsManager &      initial_conditions_manager,
    IO::RestartMgr &                    restart_manager);

  virtual ~HB();

  void notify(const StepEvent &event);

  const TimeIntg::TIAParams &getTIAParams() const
  {
    static TimeIntg::TIAParams s_tiaParams;

    return s_tiaParams;
  }

  TimeIntg::TIAParams &getTIAParams()
  {
    static TimeIntg::TIAParams s_tiaParams;

    return s_tiaParams;
  }

  // Method to set HB options
  bool setAnalysisParams(const Util::OptionBlock & option_block);
  bool setHBIntParams(const Util::OptionBlock & option_block);

    // Method to set HB linear solver / preconditioning options
  bool setHBLinSol(const Util::OptionBlock & option_block, Linear::Builder &builder);

  // Method to set non-HB linear solver / preconditioning options (needed for .STEP)
  bool setLinSol(const Util::OptionBlock & option_block);

  // Method to set time integrator options (needed for initial condition / startup periods)
  bool setTimeInt(const Util::OptionBlock & option_block);

  int getDoubleDCOPStep() const;

  bool getDCOPFlag() const;

  bool getHBtranFlags() const 
  {  return isTransient_; }


  bool useStartupICs() const
  { return useStartupICs_ ; }

protected:
  bool doRun(); 
  bool doInit(); 
  bool doLoopProcess(); 
  bool doProcessSuccessfulStep();
  bool doProcessFailedStep();
  bool doFinish();
  bool doHandlePredictor();
  bool finalVerboseOutput();

private:
  bool processSuccessfulDCOP();
  bool processFailedDCOP();

public:
  bool isAnalysis( int analysis_type ) const;

  // Transform the current solution vector for time domain and frequency domain output
  void prepareHBOutput(Linear::Vector & solnVecPtr,
                       std::vector<double> & timePoints,
                       std::vector<double> & freqPoints,
                       Teuchos::RCP<Linear::BlockVector> & timeDomainSolnVec,
                       Teuchos::RCP<Linear::BlockVector> & freqDomainSolnVecReal,
                       Teuchos::RCP<Linear::BlockVector> & freqDomainSolnVecImaginary,
                       Teuchos::RCP<Linear::BlockVector> & timeDomainLeadCurrentVec,
                       Teuchos::RCP<Linear::BlockVector> & freqDomainLeadCurrentVecReal,
                       Teuchos::RCP<Linear::BlockVector> & freqDomainLeadCurrentVecImaginary,
                       Teuchos::RCP<Linear::BlockVector> & timeDomainJunctionVoltageVec,
                       Teuchos::RCP<Linear::BlockVector> & freqDomainJunctionVoltageVecReal,
                       Teuchos::RCP<Linear::BlockVector> & freqDomainJunctionVoltageVecImaginary ) ;

private:

  // Add in solver info and timing info from current analysisObject_
  void accumulateStatistics_(AnalysisBase &analysis);

  bool runTol();
  bool runStartupPeriods();
  bool runTransientIC();
  bool interpolateIC(double initial_time);

  bool setFreqPoints_();

  bool runDCOP();
  bool setInitialGuess();


  bool setTimePoints_();

  bool createFT_(); 


  bool updateIFT_( std::vector<double>& tPoints);
  
  bool initializeOscOut( );
private:
  AnalysisManager &                     analysisManager_;
  Loader::Loader &                      loader_;
  Linear::System &                      linearSystem_;
  Nonlinear::Manager &                  nonlinearManager_;
  Device::DeviceMgr &                   deviceManager_;
  Linear::Builder &                     builder_;
  Topo::Topology &                      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;
  IO::RestartMgr &                      restartManager_;
  N_PDS_Manager *                       pdsMgrPtr_;
  AnalysisBase *                        currentAnalysisObject_;
  Loader::HBLoader *                    hbLoaderPtr_; /// HB loader, builder, system, and DFT
  Teuchos::RCP<Linear::HBBuilder>       hbBuilderPtr_;
  Linear::System *                      hbLinearSystem_;

  bool                  isPaused;               ///< Flag to indicate of the simulation is paused

  double                startDCOPtime;          /// Timing/loop count info
  double                endTRANtime;

  bool                  isTransient_;           ///< Current analysis state flags

  bool                  test_;                  ///< Testing Flag

  int                   size_;                  /// Problem Size

  std::vector<double>   freqs_;
  bool                  freqsGiven_;

  std::vector<int>      numPosFreqs;
  std::vector<int>      numFreqs_;

  double                period_;                /// Periodicity Information
  double                relErrorTol_;

  // Number of fast time periods to integrate over and IGNORE before
  // getting initial conditions for HB.  Default is zero.
  int startUpPeriods_;
  bool startUpPeriodsGiven_;

  bool saveIcData_;


  bool useStartupICs_;
// Transient assisted HB.
  int taHB_;

  bool hbOsc_;

 
  int refID_;
  std::string refNode_;

  int numTimePts_;
  bool          voltLimFlag_;
  int           intmodMax_;

  std::string   method_;

  bool intmodMaxGiven_;

  Teuchos::RCP<N_UTL_FFTInterface<std::vector<double> > > ftInterface_;
  std::vector<double> ftInData_, ftOutData_, iftInData_, iftOutData_; 

  std::vector<double>                   fastTimes_;
  std::vector<double>                   timeSteps_;
  std::vector<double>                   freqPoints_;

  // Fourier matrices
  Teuchos::RCP<N_UTL_DFTInterfaceDecl<std::vector<double> > > dftInterface_; 
  Teuchos::SerialDenseMatrix<int,double> idftMatrix_, dftMatrix_;

  // Linear solver and nonlinear solver options
  Util::OptionBlock                     saved_lsHBOB_;
  Util::OptionBlock                     saved_lsOB_;

  // Time integrator options
  Util::OptionBlock                     saved_timeIntOB_;

  // An analysis-dependent linear solver factory.
  Linear::HBSolverFactory *            solverFactory_;

  // An analysis-dependent preconditioner factory.
  Linear::HBPrecondFactory *            precFactory_;

  // Local storage vectors
  Teuchos::RCP<Linear::Vector> dcOpSolVecPtr_;
  Teuchos::RCP<Linear::Vector> dcOpStateVecPtr_;
  Teuchos::RCP<Linear::Vector> dcOpQVecPtr_;
  Teuchos::RCP<Linear::Vector> dcOpStoreVecPtr_;

  std::vector<double> goodTimePoints_;
  std::vector<Teuchos::RCP<Linear::Vector> > goodSolutionVec_;
  std::vector<Teuchos::RCP<Linear::Vector> > goodStateVec_;
  std::vector<Teuchos::RCP<Linear::Vector> > goodQVec_;
  std::vector<Teuchos::RCP<Linear::Vector> > goodStoreVec_;

  // HB initial condition
  Teuchos::RCP<Linear::BlockVector> HBICVectorPtr_;
  Teuchos::RCP<Linear::BlockVector> HBICVectorFreqPtr_;
  Teuchos::RCP<Linear::BlockVector> HBICStateVectorPtr_;  ///< HB initial state condition
  Teuchos::RCP<Linear::BlockVector> HBICQVectorPtr_;      ///< HB initial Q condition
  Teuchos::RCP<Linear::BlockVector> HBICStoreVectorPtr_;  ///< HB initial store condition

  // HB statistics
  StatCounts         hbStatCounts_;

  bool resetForStepCalledBefore_;
};

bool registerHBFactory(FactoryBlock &factory_block);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_HB_h


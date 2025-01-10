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
//
// Purpose       : This file defines the manager class for the MPDE package
//
// Special Notes :
//
// Creator       : Robert Hoekstra, 9233, Computational Sciences
//
// Creation Date : 3/11/04
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_MPDE_MANAGER_H
#define Xyce_MPDE_MANAGER_H

#include <string>
#include <map>

#include <N_ANP_fwd.h>
#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_LOA_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TIA_fwd.h>
#include <N_TOP_fwd.h>
#include <N_UTL_fwd.h>

#include <N_ANP_OutputMgrAdapter.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_System.h>
#include <N_MPDE_State.h>
#include <N_MPDE_WarpedPhaseCondition.h>
#include <N_NLS_Manager.h>
#include <N_TIA_TIAParams.h>

#include <Teuchos_RCP.hpp>

class N_MPDE_Loader;
class N_MPDE_Builder;
class N_MPDE_Discretization;

class MPDEOutputAdapter : public Xyce::Analysis::OutputAdapter
{
public:
  MPDEOutputAdapter(Xyce::Analysis::OutputMgrAdapter &adapter, Xyce::Analysis::AnalysisManager &analysis_manager, const N_MPDE_Manager &mpde_manager)
    : Xyce::Analysis::OutputAdapter(adapter),
      analysisManager_(analysis_manager),
      mpdeManager_(mpde_manager)
  {}

  void outputMPDE(double time, const Xyce::Linear::Vector *solution_vector); // override
  bool outputFunkyMPDE();
  
private:
  Xyce::Analysis::AnalysisManager &     analysisManager_;
  const N_MPDE_Manager &                mpdeManager_;
};

//-----------------------------------------------------------------------------
// Class         : N_MPDE_Manager
// Purpose       : MPDE Manager Class
// Special Notes :
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 3/11/04
//-----------------------------------------------------------------------------
class N_MPDE_Manager
{
  enum MPDE_IC
  {
    MPDE_IC_DCOP,               // 0  -- use DC sweep to approx IC
    MPDE_IC_SAWTOOTH,           // 1  -- use many small transients to approx IC.
    WAMPDE_IC_TRANSIENT,        // 2  -- use a transient sim as an approx IC.

// New algorithms are added to get better IC for complex WaMPDE circuits. genie 042914.
    MPDE_IC_TRAN_MAP,           // 3  -- use many transient sims map out an IC.
    MPDE_IC_TWO_PERIOD,         // 4  -- use two periods to interpolate IC.
    MPDE_IC_TRANSIENT           // 5  -- use a transient sim as an approx IC. This is the old ic=2. genie 042914
  };

public:
  N_MPDE_Manager(
    Xyce::Analysis::AnalysisManager &           analysis_manager,
    Xyce::Loader::Loader &                      loader,
    Xyce::Device::DeviceMgr &                   device_manager,
    Xyce::Linear::Builder &                     builder,
    Xyce::Topo::Topology &                      topology,
    Xyce::IO::InitialConditionsManager &        initial_conditions_manager,
    Xyce::IO::RestartMgr &                      restart_manager,
    const Xyce::IO::CmdParse &                  command_line);
  ~N_MPDE_Manager();

  void finalExpressionBasedSetup();

  //Runs MPDE analysis
  bool run(Xyce::Linear::System &linear_system, Xyce::Nonlinear::Manager &nonlinear_manager, Xyce::Topo::Topology &topology);

  void registerRestartManager( Xyce::IO::RestartMgr *resMgrPtr );

  // Method to register the utility options.
  bool setMPDEAnalysisParams(const Xyce::Util::OptionBlock & OB);

  // Method to register the utility options.
  bool setMPDEOptions(const Xyce::Util::OptionBlock & OB);

  // Method to register the utility options.
  bool setTransientOptions(const Xyce::Util::OptionBlock &option_block);

  // Method to register the linear solver options.
  bool setLinSolOptions(const Xyce::Util::OptionBlock & option_block);

  bool getTransientNeedsToLoadInitialConditionsAndInitializeProblem() const {
    return transientNeedsToLoadInitialConditionsAndInitializeProblem_;
  }

  bool getTransientNowCanOutputTrueMPDEResults() const {
    return transientNowCanOutputTrueMPDEResults_;
  }

  const Xyce::TimeIntg::TIAParams &getTIAParams() const {
    return tiaMPDEParams_;
  }

  Xyce::TimeIntg::TIAParams &getTIAParams() {
    return tiaMPDEParams_;
  }

  // "get" function for WaMPDE flag. (true if not IC)
  bool getWaMPDEFlag() const
  {
    return warpMPDE_;
  }

  const std::vector<double> & getFastTimePoints () const;

    const std::vector<double> & getFreqPoints () const;

  // "get" function for GID of phi variable in Warped MPDE case
  int getPhiGID () const
  {
    return warpMPDEPhasePtr_->getPhiGID();
  }

  bool initializeAll(Xyce::Linear::System &linear_system, Xyce::Nonlinear::Manager &nonlinear_manager);

  bool getOutputInterpolateMPDE() const {
    return outputInterpMPDE_;
  }

private :
  bool initializeOscOut(Xyce::Topo::Topology &topology);
  void printParams_ ();

  //Run Stages
  bool initializeMPDEAll();
  bool runInitialCondition(Xyce::Linear::System &linear_system, Xyce::Nonlinear::Manager &nonlinear_manager);
  bool runDCOP(Xyce::Linear::System &linear_system, Xyce::Nonlinear::Manager &nonlinear_manager);
  bool runStartupPeriods(const Xyce::TimeIntg::TIAParams &tia_params, Xyce::Linear::System &linear_system, Xyce::Nonlinear::Manager &nonlinear_manager);
  bool runTransientIC(const Xyce::TimeIntg::TIAParams &tia_params, Xyce::Linear::System &linear_system, Xyce::Nonlinear::Manager &nonlinear_manager);
  bool filterFastTimePoints(const Xyce::TimeIntg::TIAParams &tia_params);
  bool setupMPDEProblem_();
  bool runMPDEProblem_();
  bool runTests_();

  // diagnostic function
  double checkPeriodicity_();

private:
  const Xyce::IO::CmdParse &            commandLine_;
  Xyce::Analysis::AnalysisManager &     analysisManager_;
  Xyce::Nonlinear::Manager              nonlinearManager_;
  Xyce::Linear::System                  mpdeLinearSystem_;
  Xyce::Device::DeviceMgr &             deviceManager_;
  Xyce::Parallel::Manager &             pdsManager_;
  Xyce::Linear::Builder &               appBuilder_;
  Xyce::Topo::Topology &                topology_;
  Xyce::IO::InitialConditionsManager &  initialConditionsManager_;
  Xyce::IO::RestartMgr &                restartManager_;

  Xyce::Loader::Loader &                appLoader_;
  
  Xyce::TimeIntg::TIAParams             tiaMPDEParams_;
  Xyce::Util::OptionBlock               timeIntegratorOptionBlock_;
  Xyce::Util::OptionBlock               linSolOptionBlock_;
  MPDEOutputAdapter                     outputAdapter_;

  N_MPDE_State                          mpdeState_;
  N_MPDE_Loader *                       mpdeLoaderPtr_;
  N_MPDE_Builder *                      mpdeBuilderPtr_;
  N_MPDE_Discretization *               mpdeDiscPtr_;

  Xyce::Linear::Vector *                dcOpSolVecPtr_;
  Xyce::Linear::Vector *                dcOpStateVecPtr_;
  Xyce::Linear::Vector *                dcOpQVecPtr_;
  Xyce::Linear::Vector *                dcOpStoreVecPtr_;

  Xyce::Linear::BlockVector *           mpdeICVectorPtr_;                       /// MPDE initial condition
  Xyce::Linear::BlockVector *           mpdeICStateVectorPtr_;                  /// MPDE initial state condition
  Xyce::Linear::BlockVector *           mpdeICQVectorPtr_;                      /// MPDE initial Q condition
  Xyce::Linear::BlockVector *           mpdeICStoreVectorPtr_;                  /// MPDE initial store condition

  bool                          test_;                  ///< Testing Flag
  int                           size_;                  ///< MPDE Problem Size Factor
  bool                          tranRunForSize_;        ///< use an initial transient run to calculate size_
  int                           maxCalcSize_;           ///< max size to use from transient run.
  bool                          maxCalcSizeGiven_;

  bool                          NOOP_;                  ///< No operating point calculation

  bool                          fastSrcGiven_;
  std::vector<std::string>      srcVec_;

  std::string                   oscOut_;                ///< Independent variable for warped MPDE.
  bool                          oscOutGiven_;

  int                           nonLteSteps_;           ///< Number of steps during the start of the full MPDE
  bool                          nonLteStepsGiven_;      /// calculation before turning on truncation error control

  double                        period_;                ///< MPDE Fast Time Scale Period
  bool                          periodGiven_;

  // MPDE number of fast time periods to integrate over and IGNORE before
  // getting initial conditions for MPDE.  Default is zero.
  int                           startUpPeriods_;
  bool                          startUpPeriodsGiven_;
  bool                          saveIcData_;

  bool                          transientNeedsToLoadInitialConditionsAndInitializeProblem_;
  bool                          transientNowCanOutputTrueMPDEResults_;

  //MPDE fast time points
  // 12/5/06 tscoffe:  Note, the period T2 is the last element in fastTimes.
  // This means that the number of fast time points is fastTimes_.size()-1
  std::vector<double>           fastTimes_;

  std::vector<double>           freqPoints_;

  // Fast time discretization
  int                           fastTimeDisc_;
  int                           fastTimeDiscOrder_;

  //if we pull data directly from an initial transient run, keep
  // a list of the indices we used so that we can pull out the solution
  // and state data too.
  std::vector<int>              indicesUsed_;
  std::vector<bool>             nonPeriodicFlags;
  bool                          warpMPDE_;                      ///< warped MPDE setting in netlist

  // WaMPDE phase condition class
  N_MPDE_WarpedPhaseCondition *         warpMPDEPhasePtr_;

  // WaMPDE OSCOUT GID:
  int warpMPDEOSCOUT_;

  // WaMPDE Phase equation:
  int warpPhase_;
  bool warpPhaseGiven_;

  // WaMPDE Phase equation constant:
  double warpPhaseCoeff_;
  bool warpPhaseCoeffGiven_;

  // frequency domain flag
  bool fftFlag_;

  // Number of fast periods used for initial condition.
  int icPer_;

  // Initial condition strategy.
  int initialCondition_;

  bool          outputInterpMPDE_;       ///< flag for interpolating the MPDE output. Sometimes, the interpolation is the hard part.

  // myhsieh 081213
  // WaMPDE initial transient IC flag; if true then the initial condition for the 
  // initial fast transient run is given (e.g. .IC is pecified in the netlist to
  // kickoff the oscillation). This parameter is specifically for WaMPDE since
  // MPDE circuits have oscillations coming from OSCSRC.
  bool warpMPDEICFlag_;

  // debug flags:
  bool dcopExitFlag_;
  bool icExitFlag_;
  int  exitSawtoothStep_;

  // An analysis-dependent preconditioner factory.
  Teuchos::RCP<Xyce::Linear::PrecondFactory> precFactory_;
};

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::getFreqTimePoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 09/02/08
//-----------------------------------------------------------------------------
inline const std::vector<double> & N_MPDE_Manager::getFreqPoints() const
{
  return freqPoints_;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::getFastTimePoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 9233, Computational Sciences
// Creation Date : 04/06/04
//-----------------------------------------------------------------------------
inline const std::vector<double> & N_MPDE_Manager::getFastTimePoints() const
{
  return fastTimes_;
}

#endif //Xyce_MPDE_MANAGER_H

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
// Purpose        : Defines Manager class.
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/28/00
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_Manager_h
#define Xyce_N_NLS_Manager_h

// ---------- Standard Includes ----------
#include <map>
#include <vector>

#include <N_ANP_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_LOA_fwd.h>
#include <N_NLS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TIA_fwd.h>
#include <N_TOP_fwd.h>

#include <N_UTL_OptionBlock.h>
#include <N_UTL_Stats.h>
#include <N_NLS_ReturnCodes.h>
#include <N_NLS_NonLinearSolver.h>

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Class         : Manager
// Purpose       : Interface to the nonlinear solvers (NLS). All communication
//                 between Xyce and the nonlinear solver should be via this
//                 interface.
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------

class Manager
{
public:
  typedef unordered_map<int, Util::OptionBlock> OptionBlockMap;

  Manager(const IO::CmdParse & command_line);
  ~Manager();


  bool setOptions(const Util::OptionBlock& option_block );
  bool setTranOptions(const Util::OptionBlock& option_block );
  bool setHBOptions(const Util::OptionBlock& option_block );
  bool setNLPOptions(const Util::OptionBlock& option_block );
  bool setTwoLevelOptions (const Util::OptionBlock & option_block);
  bool setTwoLevelTranOptions (const Util::OptionBlock & option_block);
  bool setSensOptions (const Util::OptionBlock & option_block);
  bool setSensitivityOptions (const Util::OptionBlock & option_block);
  bool setLinSolOptions(const Util::OptionBlock& option_block );
  bool setLocaOptions(const Util::OptionBlock& option_block );
  bool setTwoLevelLocaOptions(const Util::OptionBlock& option_block );
  bool setDCOPRestartOptions (const Util::OptionBlock& option_block );
  bool setICOptions (const Util::OptionBlock& option_block );
  bool setNodeSetOptions (const Util::OptionBlock& option_block );
  bool setReturnCodeOption(const Util::Param &param);

  Util::OptionBlock& getHBOptions();

  bool registerSolverFactory(const Linear::SolverFactory *solver_factory)
  {
    lasSolverPtr_ = solver_factory;
    return lasSolverPtr_;
  }

  bool registerPrecondFactory(const Linear::PrecondFactory *preconditioner_factory)
  {
    lasPrecPtr_ = preconditioner_factory;
    return lasPrecPtr_;
  }

  void setReturnCodes (const ReturnCodes & retCodeTmp);

  NonLinearSolver &getNonlinearSolver() {
    return *nonlinearSolver_;
  }

  ConductanceExtractor &getConductanceExtractor() {
    return *conductanceExtractorPtr_;
  }

  const ReturnCodes &getReturnCodes() const;

  bool initializeAll(
    Analysis::AnalysisManager &         analysis_manager,
    Loader::NonlinearEquationLoader &   nonlinear_equation_loader, 
    Linear::System &                    linear_system,
    TimeIntg::DataStore &               data_store,
    Parallel::Manager &                 parallel_manager,
    IO::InitialConditionsManager &      initial_conditions_manager,
    IO::OutputMgr &                     output_manager,
    Topo::Topology &                    topology);

  int solve();

  void setAnalysisMode(AnalysisMode mode);

  // This is a more extensive version of setAnalysisMode.
  // It makes it like we are starting over.
  void resetAll(AnalysisMode mode);

  NonLinInfo getNonLinInfo() const;

  bool enableSensitivity(
      TimeIntg::DataStore & data_store, 
      TimeIntg::StepErrorControl & sec,
      Parallel::Manager &parallel_manager, 
      Topo::Topology &topology, 
      IO::OutputMgr & output_manager,
      int & numSensParams);

  bool icSensitivity(
      std::vector<double> & objectiveVec,
      std::vector<double> & dOdpVec, std::vector<double> & dOdpAdjVec,
      std::vector<double> & scaled_dOdpVec, std::vector<double> & scaled_dOdpAdjVec);
  bool calcSensitivity(
      std::vector<double> & objectiveVec,
      std::vector<double> & dOdpVec, std::vector<double> & dOdpAdjVec,
      std::vector<double> & scaled_dOdpVec, std::vector<double> & scaled_dOdpAdjVec);

  bool calcTransientAdjoint(bool timePoint,
      std::vector<double> & objectiveVec,
      std::vector<double> & dOdpVec, std::vector<double> & dOdpAdjVec,
      std::vector<double> & scaled_dOdpVec, std::vector<double> & scaled_dOdpAdjVec);

  void setMatrixFreeFlag(bool matrixFreeFlag)
  {
    matrixFreeFlag_ = matrixFreeFlag;
  }

  void allocateTranSolver(
    Analysis::AnalysisManager &         analysis_manager,
    Loader::NonlinearEquationLoader &   nonlinear_equation_loader, 
    Linear::System &                    linear_system,
    TimeIntg::DataStore &               data_store,
    Parallel::Manager &                 parallel_manager,
    IO::OutputMgr &                     output_manager,
    Topo::Topology &                    topology);

private:
  bool allocateSolver(
    Analysis::AnalysisManager &         analysis_manager,
    Loader::NonlinearEquationLoader &   nonlinear_equation_loader, 
    Linear::System &                    linear_system,
    TimeIntg::DataStore &               data_store,
    Parallel::Manager &                 parallel_manager,
    IO::InitialConditionsManager &      initial_conditions_manager,
    IO::OutputMgr &                     output_manager);

  void usingNox();
  bool setupSensitivity(
      TimeIntg::DataStore & data_store, 
      TimeIntg::StepErrorControl & sec,
      Parallel::Manager &parallel_manager, 
      Topo::Topology &topology, 
      IO::OutputMgr & output_manager,
      int & numSensParams);

private:
  const IO::CmdParse &                  commandLine_;
  NonLinearSolver *                     nonlinearSolver_;
  ConductanceExtractor *                conductanceExtractorPtr_;
  Sensitivity *                         nlsSensitivityPtr_;
  const Linear::SolverFactory *         lasSolverPtr_;
  const Linear::PrecondFactory *        lasPrecPtr_;

  bool                                  matrixFreeFlag_;
  bool                                  twoLevelNewtonFlag_;            /// Flag to determine if we are doing 2-level newton or not.
  bool                                  noxFlag_;                       /// Flag to determine if NOX is the solver in use.
  bool                                  noxFlagInner_;                  /// For 2-level newton, option for inner loop to use nox.
  bool                                  noxFlagTransient_;              /// Use nox in transient phase of calculation.
  OptionBlockMap                        optionBlockMap_;                /// netlist option blocks until we know which
  bool                                  initializeAllFlag_;
  ReturnCodes                           retCodes_;                      /// Return Codes.
  Util::Expression *                    exprPtr_;
};

bool registerPkgOptionsMgr(Manager &manager, IO::PkgOptionsMgr &options_manager);

} // namespace Nonlinear
} // namespace Xyce

#endif // Xyce_N_NLS_Manager_h

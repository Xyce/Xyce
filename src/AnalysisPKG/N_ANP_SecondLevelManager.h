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
//
// Purpose       :
//
// Special Notes :
//
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date : 01/24/08
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_SecondLevelManager_h
#define Xyce_N_ANP_SecondLevelManager_h

#include <N_ANP_AnalysisManager.h>

namespace Xyce {
namespace Analysis {

class SecondLevelManager : public AnalysisManager
{
public:
  SecondLevelManager(
    const IO::CmdParse &        command_line,
    OutputMgrAdapter &          output_manager_adapter,
    Stats::Stat                 analysis_stat)
    : AnalysisManager(command_line, output_manager_adapter, analysis_stat),
      twoLevelAnalysisObject_(0),
      activeOutput_(0),
      breakPointsRequestedBefore_(false)
  {}

  virtual ~SecondLevelManager()
  {}

private:
  SecondLevelManager(const SecondLevelManager &);
  SecondLevelManager &operator=(const SecondLevelManager &);

public:
  // Two-level Newton API functions:
  // Execute the control loop for the set analysis type,
  // for a set number of steps.
  void setExternalSolverState(Loader::CktLoader &loader, bool external_initJctFlag);
  bool runSecondLevelStep(TimeIntg::TwoLevelError & tlError);

  bool startupSecondLevelSolvers(Linear::System &linear_system, Nonlinear::Manager &nonlinear_manager);
  bool finishSecondLevelSolvers();

  void homotopyStepSuccess(const std::vector<std::string> & paramNames, const std::vector<double> & paramVals);
  void homotopyStepFailure();

  void stepSecondLevelSuccess(TwoLevelMode analysisUpper);
  void stepSecondLevelFailure(TwoLevelMode analysisUpper);
  bool getSecondLevelInitialQnorm (TimeIntg::TwoLevelError & tle) const;
  bool getSecondLevelBreakPoints(
      Loader::CktLoader &loader, 
      std::vector<Util::BreakPoint> &breakPointTimes,
      std::vector<Util::BreakPoint> &pauseBreakPointTimes
      );// const;
  bool startSecondLevelTimeStep(
    const TimeIntg::TIAParams & tia_params,
    Nonlinear::Manager &        nonlinear_manager,
    bool                        beginIntegrationFlag,
    double                      nextTimeStep,
    double                      nextTime,
    int                         currentOrder);

private:
  AnalysisBase *        twoLevelAnalysisObject_;
  IO::ActiveOutput *    activeOutput_;


  bool breakPointsRequestedBefore_;
};

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_SecondLevelManager_h

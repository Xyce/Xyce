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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/27/00
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_CIR_SecondLevelSimulator_h
#define Xyce_N_CIR_SecondLevelSimulator_h

#include <N_CIR_Xyce.h>
#include <N_ANP_SecondLevelManager.h>

namespace Xyce {
namespace Circuit {

//-----------------------------------------------------------------------------
// Class         : Xyce
// Purpose       : This is the main "top level" class for Xyce.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
class SecondLevelSimulator : public Simulator
{
 public:
  SecondLevelSimulator(Parallel::Machine comm = MPI_COMM_NULL)
    : Simulator(comm),
      secondLevelManager_(0)
  {}

  virtual ~SecondLevelSimulator()
  {}

  virtual Analysis::AnalysisManager *newAnalysisManager(
    const IO::CmdParse &                command_line,
    IO::RestartMgr &                    restart_manager,
    Analysis::OutputMgrAdapter &        output_manager_adapter,
    Stats::Stat                         analysis_stat);

  bool simulateStep(
    bool                                        external_initJctFlag,
    const std::map<std::string,double> &        inputMap,
    std::vector<double> &                       outputVector,
    std::vector< std::vector<double> > &        jacobian,
    TimeIntg::TwoLevelError &                   tlError);

  bool simulateStep(
    const Device::ExternalSimulationData &      ext_data,
    const std::map<std::string,double> &        inputMap,
    std::vector<double> &                       outputVector,
    std::vector< std::vector<double> > &        jacobian,
    TimeIntg::TwoLevelError &                   tlError);

  // 2-level, power node functions:
  bool startupSolvers();
  bool finishSolvers();

  void homotopyStepSuccess(
    const std::vector<std::string> & paramNames,
    const std::vector<double> & paramVals);

  void homotopyStepFailure();

  void stepSuccess(Analysis::TwoLevelMode analysis);
  void stepFailure(Analysis::TwoLevelMode analysis);
  bool getBreakPoints (
      std::vector<Util::BreakPoint> &breakPointTimes,
      std::vector<Util::BreakPoint> &pauseBreakPointTimes
      );
  bool updateStateArrays ();
  bool startTimeStep(
    bool                          beginIntegrationFlag,
    double                        nextTimeStep,
    double                        nextTime,
    int                           currentOrder);
  
  bool startTimeStep(const Device::ExternalSimulationData & ext_data);
  bool endTimeStep (Device::ExternalSimulationData & ext_data);

  bool setInternalParam (const std::string & name, double val);

  bool getInitialQnorm (TimeIntg::TwoLevelError & tle);

private:
  Analysis::SecondLevelManager *        secondLevelManager_;
};

} // namespace Circuit
} // namespace Xyce

#endif // Xyce_N_CIR_SecondLevelSimulator_h

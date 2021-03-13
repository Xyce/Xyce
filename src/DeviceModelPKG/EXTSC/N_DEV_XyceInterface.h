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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/15/05
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_XyceInterface_h
#define Xyce_N_DEV_XyceInterface_h

#include <string>
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_ExternCodeInterface.h>
#include <N_UTL_fwd.h>
#include <N_ANP_fwd.h>
#include <N_CIR_fwd.h>
#include <N_TIA_fwd.h>
#include <N_IO_CmdParse.h>

#include <N_CIR_SecondLevelSimulator.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : XyceInterface
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, Parallel Computational Sciences
// Creation Date : 04/15/05
//-----------------------------------------------------------------------------
class XyceInterface : public ExternCodeInterface
{
  public:
    XyceInterface(
      const DeviceOptions & do1,
      const IO::CmdParse & cp,
      const std::string & netlist);

    virtual ~XyceInterface();

  private:
    XyceInterface (const XyceInterface &right);

  public:
    bool initialize(Parallel::Communicator* comm = 0);

  bool simulateStep(
    const SolverState &                         solState,
    const std::map<std::string,double> &        inputMap,
    std::vector<double> &                       outputVector,
    std::vector< std::vector<double> > &        jacobian,
    TimeIntg::TwoLevelError &                   tlError);

    bool finalize ();
    bool run ();

    void homotopyStepSuccess(
      const std::vector<std::string> &  paramNames,
      const std::vector<double> &       paramVals);

    void homotopyStepFailure ();

    void stepSuccess(Analysis::TwoLevelMode analysis);
    void stepFailure(Analysis::TwoLevelMode analysis);

    bool getBreakPoints (
        std::vector<Util::BreakPoint> &breakPointTimes,
        std::vector<Util::BreakPoint> &pauseBreakPointTimes);

    bool updateStateArrays ();
  bool startTimeStep(
    bool        beginIntegrationFlag,
    double      nextTimeStep,
    double      nextTime,
    int         currentOrder);

    bool setInternalParam (const std::string & name, double val);

    bool getInitialQnorm (TimeIntg::TwoLevelError & tle);

  private:
    std::string                         netlistFilename_;
    Circuit::SecondLevelSimulator *     simulator_;
    IO::CmdParse                        commandLine_;
};

} // namespace Device
} // namespace Xyce

#endif


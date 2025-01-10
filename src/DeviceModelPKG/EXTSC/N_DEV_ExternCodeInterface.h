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

#ifndef Xyce_N_DEV_ExternCodeInterface_h
#define Xyce_N_DEV_ExternCodeInterface_h

#include <vector>
#include <string>
#include <map>

#include <N_DEV_fwd.h>
#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>
#include <N_ANP_fwd.h>
#include <N_TIA_fwd.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : N_DEV_ExternCodeInterface
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, Parallel Computational Sciences
// Creation Date : 04/15/05
//-----------------------------------------------------------------------------
class ExternCodeInterface
{
  public:
    ExternCodeInterface ();

    ExternCodeInterface (const ExternCodeInterface &right);
    virtual ~ExternCodeInterface();

    virtual bool initialize(Parallel::Machine comm = MPI_COMM_NULL) = 0;

    virtual bool simulateStep(
      const SolverState & solState,
      const std::map< std::string, double > & inputMap,
      std::vector<double> & outputVector,
      std::vector< std::vector<double> > & jacobian,
      TimeIntg::TwoLevelError & tlError) = 0;

    virtual bool finalize () = 0;
    virtual bool run () = 0;

    virtual void homotopyStepSuccess(
      const std::vector<std::string> & paramNames,
      const std::vector<double> & paramVals) = 0;

    virtual void homotopyStepFailure () = 0;

    virtual void stepSuccess(Analysis::TwoLevelMode analysis) = 0;
    virtual void stepFailure(Analysis::TwoLevelMode analysis) = 0;
    virtual bool getBreakPoints (
        std::vector<Util::BreakPoint> &breakPointTimes,
        std::vector<Util::BreakPoint> &pauseBreakPointTimes
        ) = 0;
    virtual bool updateStateArrays () = 0;

  virtual bool startTimeStep(
    bool        beginIntegrationFlag,
    double      nextTimeStep,
    double      nextTime,
    int         currentOrder) = 0;

  virtual bool setInternalParam (const std::string & name, double val) = 0;

    virtual bool getInitialQnorm (TimeIntg::TwoLevelError & tle) = 0;
};

} // namespace Device
} // namespace Xyce

#endif


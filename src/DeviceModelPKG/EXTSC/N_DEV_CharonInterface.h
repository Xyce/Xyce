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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/13/06
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_CharonInterface_h
#define Xyce_N_DEV_CharonInterface_h

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include <N_DEV_fwd.h>
#include <N_CIR_fwd.h>
#include <N_DEV_ExternCodeInterface.h>
#include <N_UTL_fwd.h>
#include <N_TIA_fwd.h>

namespace Teuchos {
  class ParamList;
}

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : CharonInterface
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, Parallel Computational Sciences
// Creation Date : 03/13/06
//-----------------------------------------------------------------------------
class CharonInterface : public ExternCodeInterface
{
  public:
    CharonInterface (
      const DeviceOptions &     do1,
      const std::string &       netlist,
      const SolverState &       ss1);

    CharonInterface (const CharonInterface &right);
    virtual ~CharonInterface();

    bool initialize(N_PDS_Comm * comm = 0);

  bool simulateStep(
    const SolverState &                         solState,
    const std::map<std::string,double> &        inputMap,
    std::vector<double> &                       outputVector,
    std::vector< std::vector<double> > &        jacobian,
    TimeIntg::TwoLevelError &                   tlError);

    bool finalize ();
    bool run ();

    void homotopyStepSuccess
      (const std::vector<std::string> & paramNames,
       const std::vector<double> & paramVals);

    void homotopyStepFailure ();

    void stepSuccess(Analysis::TwoLevelMode analysis);
    void stepFailure(Analysis::TwoLevelMode analysis);
    bool getBreakPoints (
        std::vector<Util::BreakPoint> &breakPointTimes,
        std::vector<Util::BreakPoint> &pauseBreakPointTimes) { return true; }
    bool updateStateArrays () {return true;}

  bool startTimeStep(
    bool        beginIntegrationFlag,
    double      nextTimeStep,
    double      nextTime,
    int         currentOrder)
  {
    return true;
  }

  bool setInternalParam (const std::string & name, double val) {return true;}

    bool getInitialQnorm (TimeIntg::TwoLevelError & tle);

  private:
    std::string inputFileName_;
    const DeviceOptions& devOptions_;
    const SolverState & solState_;

    // The "command line" arguments
    std::vector<char*> cargs_;

    //! Input list for Charon.
    Teuchos::RCP<Teuchos::ParameterList> input_list_;

    //! Output list from Charon.
    Teuchos::RCP<Teuchos::ParameterList> output_list_;
};

} // namespace Device
} // namespace Xyce

#endif


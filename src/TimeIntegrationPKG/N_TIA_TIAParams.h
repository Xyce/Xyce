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
// Purpose       : This file identifies the class associated with all user
//                 specified parameters which relate to the time integration
//                 algorithms and problem definition.
//
// Special Notes :
//
// Creator       : Buddy Watts
//
// Creation Date : 6/1/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_PARAMS_H
#define Xyce_N_TIA_PARAMS_H

#include <list>
#include <string>
#include <vector>

#include <N_IO_fwd.h>
#include <N_UTL_Param.h>

namespace Xyce {
namespace TimeIntg {

int maxOrder(const IO::CmdParse &command_line);

enum ErrorAnalysisOption {LOCAL_TRUNCATED_ESTIMATES = 0, NO_LOCAL_TRUNCATED_ESTIMATES = 1};

//-----------------------------------------------------------------------------
// Class         : TIAParams::TIAParams
// Purpose       : This is a class that sets up time integration information,
//                 associated with user related input only, that defines the
//                 problem setup.
// Special Notes :
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
class TIAParams
{
public:
  static void populateMetadata(IO::PkgOptionsMgr &options_manager);

  TIAParams();
  TIAParams(const TIAParams &);

  TIAParams &operator=(const TIAParams & right);
  ~TIAParams();

  void setMaxOrder(int max_order);

  bool setTimeIntegratorOption(const Util::Param &param);
  bool setAnalysisOption(const Util::Param &param);
  bool updateAnalysisOptions();

  // Print out time-integration parameters.
  void printParams(std::ostream &os, int analysis) const;

public:
  double        initialTime;                    ///< Beginning time for the time integrator (StepErrorControl, integrators access from StepErrorControl)
  double        finalTime;                      ///< End time for simulation
  double        initialTimeStep;                ///< User specified initial time step
  double        minTimeStep;                    ///< User specified minimum time step
  bool          minTimeStepGiven;
  int           minTimeStepsBP;                 ///< User specified mininum number of steps between breakpoints
  bool          minTimeStepsBPGiven;
  double        maxTimeStep;                    ///< User specified maximum time step
  bool          maxTimeStepGiven;
  bool          constantTimeStepFlag;           ///< Constant time step integration flag
  double        restartTimeStepScale;           ///< Time step is scaled coming out of a breakpoint (StepErrorControl)

  double        initialOutputTime;              ///< Time at which output starts (StepErrorControl)
  bool          initialOutputTimeGiven;

  bool          useDeviceTimeStepMaxFlag;       ///< True if devices impose a time step maximum

  int           errorAnalysisOption;            ///< Error analysis option
  int           errorAnalysisOptionResetCount;  ///< Iteration count down to reset errorAnalysisOption to LOCAL_TRUNCATED_ESTIMATES

  bool          bpEnable;                       ///< Enable breakpoints flag

  // Iteration count algorithm parameters
  int           NLmin;
  int           NLmax;
  double        delmax;
  bool          delmaxGiven;
  bool          timestepsReversal;
  bool          testFirstStep;

  bool          newBPStepping;

  bool maskIVars;
  int newLte;

  // Error Tolerances:

    // Relative error tolerance.  This value should be selected to be 10^(-(m+1))
    // where m is the desired number of significant digits in the solution.
    double relErrorTol;
    bool relErrorTolGiven;

    // Absolute error tolerance.  In general, this should be selected to be small
    // enough such that any solution value smaller than this tolerance will be
    // considered "insignificant".
    double absErrorTol;

    // Error acceptance tolerance.
    double errTolAcceptance;

    // Jacobian limit - to prevent capacitive spiral of death!
  bool jacLimitFlag;
  double jacLimit;

    // Maximum/Minimum order desired (BackwardDifferentiation15 specific)
    int maxOrder;
    int minOrder;

    // flag for interpolated output
    bool interpOutputFlag;

    // if we have rejected several time steps and are about to fail due to a
    // time step too smal error, we'll go back and accept the step that had
    // the minimum estimated error over tol if this flag is true.
    // (set by the user in the netlist via .options timeint mintimesteprecovery=<int>)
    int minTimeStepRecoveryCounter;

    bool bpPrune;

    std::vector< std::pair< Util::Param, double * > > dependentOptions;
};

// To be moved??

// errorAnalysisOptionResetCount
// minTimeStepRecoveryCounter

// Used by TimeIntegrationMethods and StepErrorControl
// NLmax
// constantTimeStepFlag
// errorAnalysisOption

// Used by WorkingIntegrationMethod
// jacLimitFlag
// jacLimit

// Used by StepErrorControl
// NLmax
// constantTimeStepFlag
// delmax
// delmaxGiven
// errTolAcceptance
// errorAnalysisOption
// finalTime
// initialOutputTime
// initialTime
// initialTimeStep
// maxTimeStep
// maxTimeStepGiven
// minTimeStep
// minTimeStepGiven
// newBPStepping
// restartTimeStepScale
// testFirstStep
// timestepsReversal
// useDeviceTimeStepMaxFlag
// userSpecBreakPoints

// Used by Gear12
// NLmax
// NLmin
// constantTimeStepFlag
// errorAnalysisOption
// interpOutputFlag
// maxOrder
// minOrder

// Used by BackwardDifferentiation15
// constantTimeStepFlag
// errorAnalysisOption
// interpOutputFlag
// maxOrder
// minOrder

// Used by OneStep
// NLmax
// NLmin
// constantTimeStepFlag
// errorAnalysisOption
// interpOutputFlag
// maxOrder
// minOrder

// Used by DataStore
// newLte
// maxOrder (indirectly)
// absErrorTol (indirectly)
// relErrorTol (indirectly)

} // namespace TimeIntg
} // namespace Xyce

#endif // Xyce_N_TIA_TIAParams_H

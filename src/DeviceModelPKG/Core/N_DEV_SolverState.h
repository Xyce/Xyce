//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
// Purpose        : This is a container class for solver information.
//                  It may occasionally contain stuff that isn't strictly
//                  pertaining to the solver state, but that is its primary
//                  intention.
//
//                  In general, stuff that goes into this class should
//                  be stuff needed by more than one device instance type.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_SolverState_h
#define Xyce_N_DEV_SolverState_h

#include <map>
#include <string>
#include <vector>

#include <N_ANP_fwd.h>
#include <N_UTL_fwd.h>
#include <N_DEV_fwd.h>
#include <N_PDS_fwd.h>

#include <N_NLS_TwoLevelEnum.h>
#include <N_NLS_fwd.h>


namespace Xyce {
namespace Device {

struct Globals 
{
  GlobalParameterMap global_params;
  std::vector<Util::Expression> global_expressions;
  std::vector<std::string>global_exp_names;
};

//-----------------------------------------------------------------------------
// Class         : SolverState
// Purpose       : Container for current solver data.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/25/02
//-----------------------------------------------------------------------------
class SolverState
{
  friend std::ostream& operator<<(std::ostream& os, const SolverState & ss);

public:
  SolverState();

  void initializeHomotopyBlockSize(int numBlocks);

  Globals &getGlobals() const
  {
    return const_cast<Globals &>(globals_);
  }

public:

  // Startup information
  bool                  isPDESystem_;           ///< true if circuit includes a PDE device

  // Time integrator information
  double                pdt_;                   ///< Previous delta time alpha/dt (Many devices)
  int                   currentOrder_;          ///< ROM
  int                   usedOrder_;             ///< ROM

  // Time step error control information
  double                currTimeStep_;          ///< Region, BJT, Digital, ThermalResistor, ROM, Charon, Others
  double                lastTimeStep_;          ///< BJT, Others
  double                currTime_;              ///< DeviceEntity for expression time, breakpoints
                                                ///< DeviceMgr for dependent parameters, breakpoints, extern device
                                                ///< SourceData devices, ADC, DAC LTRA, TRA, Region, NumericalJacobian, RxnSet, Xygra, Digital
                                                ///< 2D PDE, Diode PDE, Charon, Synapse, Neuron, Others
  double                finalTime_;             ///< Analysis final time, SourceData devices
  double                startingTimeStep_;      ///< SourceData devices
  double                bpTol_;                 ///< Break point tolerance, SourceData devices, Neuron devices

  // Mixed signal
  double                acceptedTime_;          ///< DeviceMgr::acceptStep(), DAC (for habanero)

  // MPDE stuff
  bool                  mpdeOnFlag_;            ///< MPDE phase of MPDE problem (ie not initial condition)
  double                currFastTime_;          ///< Source devices, 

  bool                  blockAnalysisFlag_;     ///< Source devices, BJTDW, This indicates an MPDE/HB run.  This is true during both IC and MPDE/HB phase.

  bool spAnalysisFlag_;

  // output flag:
  bool   doubleDCOPEnabled;   // true if taking 2 DCOP steps for PDE sim.
  int    doubleDCOPStep;      // 0 or 1.  (first or second "double" DCOP).

  int                   timeStepNumber_;        ///< Memristor, LTRA, TRA, testing if debug or jacobian for testing

  // The following "ltra*" data structures are used to track time history for
  // the LTRA device. It requires an independent time history because
  // it can be compacted if the user specifies that option. With no
  // compaction ltraTimeStepNumber will be equal to timeStepNumber+1.
  bool                  ltraDevices_;           ///< LTRA, true when LTRA devices are in the circuit
  int ltraTimeIndex_;                           ///< LTRA, DeviceMgr::acceptStep()
  int ltraTimeHistorySize_;                     ///< LTRA, this looks like c code array sizing
  mutable bool          ltraDoCompact_;         ///< LTRA
  std::vector<double>   ltraTimePoints_;        ///< LTRA

  int                   newtonIter;
  int                   continuationStepNumber;
  bool                  firstContinuationParam;
  bool                  firstSolveComplete;

  bool                  initTranFlag_;          ///< RxnSet, TRA, LTRA, ACC, MOSFET, BJT, true only on very first(t=0) time step.
  bool                  beginIntegrationFlag_;  ///< BJT, true if 1st time step out of breakpoint (incl. t=0)

  bool dcopFlag;         // true if we are in a DCOP calculation (sweep, tranop or OP)
  bool inputOPFlag;       // true if starting from a previous OP calculation

  bool transientFlag;    // true if transient analysis(even during tranop)
  bool dcsweepFlag;      // true if DC Sweep or OP calculation.
  bool tranopFlag;       // true if in dcop phase of transient sim.
  bool acopFlag;         // true if in acop phase of ac sim.
  bool noiseFlag;         // true if we are in a noise calculation 

  bool locaEnabledFlag;  // true if LOCA is enabled for DCOP.

  bool                  externalInitJctFlag_;
  bool                  externalStateFlag_;
  bool                  initJctFlag_;           ///< true if on the first newton step of the first dcop solve of the first .STEP iteration. BJT, JFET, Diode, MOSFET, SW, Extern

  bool                  initFixFlag;            // true if DCOP solve, not first iteration *AND* any device not converged.  Allows "OFF" to be applied.

  bool sweepSourceResetFlag;
  bool debugTimeFlag;

  Nonlinear::TwoLevelNewtonMode twoLevelNewtonCouplingMode;

  // Homotopy

  // PDE device BC homotopy/two-level newton
  double                pdeAlpha_;              ///< PDEAlphaParam of ArtificialParameters

  bool                  PDEcontinuationFlag_;   ///< PDE enable/disablePDSContinuation(), VsrcScaleParam, PDEBetaParam, PDEAlphaParam, Diode PDE Device, 2d PDE devices

  bool                  chargeHomotopy_;        ///< 2d PDE Devices, ArtificialParameters
  double                chargeAlpha_;           ///< 2d PDE Devices, ArtificialParameters

  // MOSFET homotopy variables
  bool                  artParameterFlag_;      ///< MOSFET Devices, ArtificialParameters
  std::vector<double>   gainScale_;             ///< MOSFET Devices, ArtificialParameters
  double                nltermScale_;           ///< MOSFET Devices, ArtificialParameters

  bool                  sizeParameterFlag_;     ///< ArtificialParameters, not sure these are really used
  double                sizeScale_;             ///< ArtificialParameters

  Globals               globals_;

  double                currFreq_;              ///< current frequency
};

bool setupSolverInfo(
  SolverState &                         solver_state,
  const Analysis::AnalysisManager &     analysis_manager,
  bool                                  all_devices_converged,
  const DeviceOptions &                 device_options,
  const Nonlinear::NonLinInfo &         nonlinear_info);

} // namespace Device
} // namespace Xyce

#endif


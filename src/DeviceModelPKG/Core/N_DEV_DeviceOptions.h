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
// Purpose        :
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


#ifndef Xyce_N_DEV_DeviceOptions_h
#define Xyce_N_DEV_DeviceOptions_h

#include <string>

#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_UTL_fwd.h>
#include <N_UTL_Param.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : DeviceOptions
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class DeviceOptions
{
  friend std::ostream & operator<<(std::ostream & os, const DeviceOptions & devOp);

public:
  static void populateMetadata(IO::PkgOptionsMgr &options_manager);

  DeviceOptions();
  ~DeviceOptions();

  bool setOptions(const Util::OptionBlock &option_block);
  bool setParserOptions(const Util::OptionBlock &option_block);

private:
  DeviceOptions(DeviceOptions const &);                 ///< No copying
  DeviceOptions & operator=(DeviceOptions const &);     ///< No assignment

public:
  // some general MOS
  double        defad;                                ///< MOS drain diffusion area.
  double        defas;                                ///< MOS source diffusion area.
  double        defl;                                 ///< MOS channel length.
  double        defw;                                 ///< MOS channel width.
  bool          modelBinningFlag;
  double        lengthScale;

  double        abstol;                               ///< absolute current error tolerance.
  double        reltol;                               ///< relative current error tolerance.
  double        chgtol;                               ///< absolute charge error tolerance.

  double        gmin;                                 ///< minimum allowed conductance.
  double        gmin_orig;                            ///< this is needed for gmin-homotopy.
  double        gmin_init;                            ///< this is needed for gmin-homotopy.
  double        gmin_scalar;                          ///< this is needed for gmin-homotopy.

  double        gmax;                                 ///< maximum allowed conductance.

  double        testJac_relTol;                       ///< reltol for num. jacobian diagnostic
  double        testJac_absTol;                       ///< abstol for num. jacobian diagnostic.
  double        testJac_SqrtEta;                      ///< dx = numJacSqrtEta * (1.0 + std::fabs(soln[i]));
  double        deviceSens_dp;                        ///< similar to eta, but for numerical device sensitivities

  double        tnom;                                 ///< nominal temperature for device params.
  Util::Param   temp;                                 ///< operating temperature of ckt.

  bool          matrixSensitivityFlag;
  bool          testJacobianFlag;
  int           testJacStartStep;
  int           testJacStopStep;
  bool          testJacWarn;
  bool          testJacDeviceNameGiven;
  std::string   testJacDeviceName;
  bool          voltageLimiterFlag;
  int           lambertWFlag;

  bool          newMeyerFlag;

  double        icMultiplier;

  double        defaultMaxTimeStep;

  ///< mosfet homotopy:
  double        vdsScaleMin;
  double        vgstConst;
  int           numGainScaleBlocks;
  bool          staggerGainScale;
  bool          randomizeVgstConst;
  double        length0;
  double        width0;
  double        tox0;
  double        minRes;
  double        minCap;
  double        exp_order;

  ///< tolerance on resistance below which it will be treated as zero
  double        zeroResistanceTol;
  bool          checkForZeroResistance;

  int           debugMinTimestep;
  int           debugMaxTimestep;
  double        debugMinTime;
  double        debugMaxTime;

  int           verboseLevel;

  double        rc_const;
  bool          newABM;
  bool          newExcessPhase;
  bool          defaultNewExcessPhase;                ///< default is true for MPDE, false for non-MPDE.

  double        excessPhaseScalar1;
  double        excessPhaseScalar2;

  long  randomSeed;                           ///< seed for random number generator used by some devices.
  ///< note: each device gets its own random number generator so
  ///< it must initialize thing correctly. (See N_DEV_Synapse3 for an
  ///< example)

  bool          tryToCompact;                         ///< Try to compact past history for LTRA device(s).

  bool          calculateAllLeadCurrents;             ///< configure all devices to load lead currents

  int           digInitState;   ///< used to initialize all of the digital gates in a circuit

  bool          separateLoad;   ///< used to enable separated device loading

  bool          pwl_BP_off;     ///< if true, then PWL sources have no breakpoints
};

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_DeviceOptions_h

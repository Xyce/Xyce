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
// Purpose       : 
// Special Notes :
// Creator       : Dave Baur
// Creation Date : 
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iomanip>

#include <N_UTL_Diagnostic.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_SaveIOSState.h>

namespace Xyce {
namespace Diag {

template<>
int &getMask<TimeIntegrator>()
{
  static int s_timeIntegratorMask = 0;

  return s_timeIntegratorMask;
}

template<>
int &getMask<Device>()
{
  static int s_deviceMask = 0;

  return s_deviceMask;
}

template<>
int &getMask<HB>()
{
  static int s_hbMask = 0;

  return s_hbMask;
}

template<>
int &getMask<Sensitivity>()
{
  static int s_sensitivitiyMask = 0;

  return s_sensitivitiyMask;
}

template<>
int &getMask<MPDE>()
{
  static int s_mpdeMask = 0;

  return s_mpdeMask;
}

template<>
int &getMask<Nonlinear>()
{
  static int s_nonlinearMask = 0;

  return s_nonlinearMask;
}

template<>
int &getMask<IO>()
{
  static int s_ioMask = 0;

  return s_ioMask;
}

} // namespace Diag


//-----------------------------------------------------------------------------
// Function      : setTimeIntegratorDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:59:47 2015
//-----------------------------------------------------------------------------
///
/// Sets the time intergrator selection mask based on the specified level.
///
/// @param level Selection level
///
void
setTimeIntegratorDebugLevel(
  int           level)
{
  Diag::getMask<Diag::TimeIntegrator>() = 0;

  if (DEBUG_TIME) {
    if (level == 1)
      Diag::getMask<Diag::TimeIntegrator>() |= Diag::TIME_DEFAULT;
    if (level == 2)
      Diag::getMask<Diag::TimeIntegrator>() |= Diag::TIME_STEP | Diag::TIME_OUTPUT;
    if (level == 3 || level == 4)
      Diag::getMask<Diag::TimeIntegrator>() |= Diag::TIME_COEFFICIENTS | Diag::TIME_PREDICTOR | Diag::TIME_RESIDUAL | Diag::TIME_JACOBIAN | Diag::TIME_HISTORY | Diag::TIME_BREAKPOINTS;
    if (level == 5)
      Diag::getMask<Diag::TimeIntegrator>() |= Diag::TIME_DUMP_SOLUTION_ARRAYS | Diag::TIME_ERROR;

    basic_ios_all_saver<std::ostream::char_type> save(Xyce::lout());
    Xyce::lout() << "Time integrator debug is set to " << level << " (0x" << std::setw(8) << std::setfill('0') << std::hex << Diag::getMask<Diag::TimeIntegrator>() << ")" << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : setDeviceDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:59:47 2015
//-----------------------------------------------------------------------------
///
/// Sets the device selection mask based on the specified level.
///
/// @param level Selection level
///
void
setDeviceDebugLevel(
  int           level)
{
  Diag::getMask<Diag::Device>() = 0;

  if (DEBUG_DEVICE) {
    if (level > 0)
      Diag::getMask<Diag::Device>() |= Diag::DEVICE_DEFAULT;
    if (level > 1)
      Diag::getMask<Diag::Device>() |= Diag::DEVICE_JACSTAMP | Diag::DEVICE_LIDS | Diag::DEVICE_LOAD_VECTOR | Diag::DEVICE_PRINT_VECTORS | Diag::DEVICE_SOLVER_STATE;
    if (level > 2)
      Diag::getMask<Diag::Device>() |= Diag::DEVICE_DUMP_VECTORS;

    basic_ios_all_saver<std::ostream::char_type> save(Xyce::lout());

    Xyce::lout() << "Device debug is set to " << level << " (0x" << std::setw(8) << std::setfill('0') << std::hex << Diag::getMask<Diag::Device>() << ")" << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : setHBDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:59:47 2015
//-----------------------------------------------------------------------------
///
/// Sets the harmonic balance selection mask based on the specified level.
///
/// @param level Selection level
///
void
setHBDebugLevel(
  int           level)
{
  Diag::getMask<Diag::HB>() = 0;

  if (DEBUG_HB) {
    if (level > 0)
      Diag::getMask<Diag::HB>() |= Diag::HB_DEFAULT;
    if (level > 1)
      Diag::getMask<Diag::HB>() |= Diag::HB_FAST_TIMES | Diag::HB_TIMESTEP;
    if (level > 2)
      Diag::getMask<Diag::HB>() |= Diag::HB_PRINT_VECTORS;

    basic_ios_all_saver<std::ostream::char_type> save(Xyce::lout());

    Xyce::lout() << "Harmonic balance is set to " << level << " (0x" << std::setw(8) << std::setfill('0') << std::hex << Diag::getMask<Diag::HB>() << ")" << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : setMPDEDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:59:47 2015
//-----------------------------------------------------------------------------
///
/// Sets the multi-time PDE selection mask based on the specified level.
///
/// @param level Selection level
///
void
setMPDEDebugLevel(
  int           level)
{
  Diag::getMask<Diag::MPDE>() = 0;

  if (DEBUG_MPDE) {
    if (level > 0)
      Diag::getMask<Diag::MPDE>() |= Diag::MPDE_DEFAULT;
    if (level > 1)
      Diag::getMask<Diag::MPDE>() |= Diag::MPDE_PRINT_VECTORS;

    basic_ios_all_saver<std::ostream::char_type> save(Xyce::lout());

    Xyce::lout() << "Multi-time PDE debug level is set to " << level << " (0x" << std::setw(8) << std::setfill('0') << std::hex << Diag::getMask<Diag::MPDE>() << ")" << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : setDeviceSensitivityDebugLevel
// Purpose       : This applies to a few debug statements in the sensitivity 
//                 code in the device package
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:59:47 2015
//-----------------------------------------------------------------------------
///
/// Sets the sensitivity selection mask based on the specified level.
///
/// @param level Selection level
///
void
setDeviceSensitivityDebugLevel(
  int           level)
{
  Diag::getMask<Diag::Sensitivity>() &= ~Diag::SENS_PARAMETERS;

  if (DEBUG_NONLINEAR) {
    if (level > 0)
      Diag::getMask<Diag::Sensitivity>() |= Diag::SENS_PARAMETERS;

    basic_ios_all_saver<std::ostream::char_type> save(Xyce::lout());

    Xyce::lout() << "Sensitivity debug level is set to " << level << " (0x" << std::setw(8) << std::setfill('0') << std::hex << Diag::getMask<Diag::Sensitivity>() << ")" << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : setSensitivityDebugLevel
// Purpose       : This applies to a few debug statements in the sensitivity 
//                 code in the nonlinear solver (ie, most of the .SENS capability)
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:59:47 2015
//-----------------------------------------------------------------------------
///
/// Sets the nonlinear solver mask based on the specified sensitivity level.
///
/// @param level Selection level
///
void
setSensitivityDebugLevel(
  int           level)
{
  Diag::getMask<Diag::Sensitivity>() &= ~Diag::SENS_SOLVER;

  if (DEBUG_NONLINEAR) {
    if (level > 1)
      Diag::getMask<Diag::Sensitivity>() |= Diag::SENS_SOLVER;

    basic_ios_all_saver<std::ostream::char_type> save(Xyce::lout());

    Xyce::lout() << "Sensitivity debug level is set to " << level << " (0x" << std::setw(8) << std::setfill('0') << std::hex << Diag::getMask<Diag::Sensitivity>() << ")" << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : setNonlinearDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:59:47 2015
//-----------------------------------------------------------------------------
///
/// Sets the nonlinear solver selection mask based on the specified level.
///
/// @param level Selection level
///
void
setNonlinearDebugLevel(
  int           level)
{
  Diag::getMask<Diag::Nonlinear>() &= ~Diag::NONLINEAR_PARAMETERS;

  if (DEBUG_NONLINEAR) {
    if (level > 1)
      Diag::getMask<Diag::Nonlinear>() |= Diag::NONLINEAR_PARAMETERS;

    basic_ios_all_saver<std::ostream::char_type> save(Xyce::lout());

    Xyce::lout() << "Nonlinear debug level is set to " << level << " (0x" << std::setw(8) << std::setfill('0') << std::hex << Diag::getMask<Diag::Nonlinear>() << ")" << std::endl;
  }
}


//-----------------------------------------------------------------------------
// Function      : setNonlinearConductanceDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:59:47 2015
//-----------------------------------------------------------------------------
///
/// Sets the nonlinear selection mask based on the specified conductance level.
///
/// @param level Selection level
///
void
setNonlinearConductanceDebugLevel(
  int           level)
{
  Diag::getMask<Diag::Nonlinear>() &= ~Diag::NONLINEAR_CONDUCTANCE;

  if (DEBUG_NONLINEAR) {
    if (level > 1)
      Diag::getMask<Diag::Nonlinear>() |= Diag::NONLINEAR_CONDUCTANCE;

    basic_ios_all_saver<std::ostream::char_type> save(Xyce::lout());

    Xyce::lout() << "Nonlinear conductance debug level is set to " << level << " (0x" << std::setw(8) << std::setfill('0') << std::hex << Diag::getMask<Diag::Nonlinear>() << ")" << std::endl;
  }
}


//-----------------------------------------------------------------------------
// Function      : setNonlinearDumpDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:59:47 2015
//-----------------------------------------------------------------------------
///
/// Sets the nonlinear solver selection mask based on the specified dump level.
///
/// @param level Selection level
///
void
setNonlinearDumpDebugLevel(
  int           level)
{
  Diag::getMask<Diag::Nonlinear>() &= ~Diag::NONLINEAR_DUMP_MASK;

  if (DEBUG_NONLINEAR) {
    if (level > 1)
      Diag::getMask<Diag::Nonlinear>() |= Diag::NONLINEAR_DUMP_PARAM_NUMBER;
    else if (level > 2)
      Diag::getMask<Diag::Nonlinear>() |= Diag::NONLINEAR_DUMP_STEP;

    basic_ios_all_saver<std::ostream::char_type> save(Xyce::lout());

    //Xyce::lout() << "Nonlinear dump debug level is set to " << level << " (0x" << std::setw(8) << std::setfill('0') << std::hex << Diag::getMask<Diag::Nonlinear>() << ")" << std::endl;
  }
}


} // namespace Xyce

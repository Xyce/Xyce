//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
//
// Creator        : Dave Baur
//
// Creation Date  : 01/11
//
//
//-----------------------------------------------------------------------------

///
/// @file   N_UTL_Diagnostic.h
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
/// @date   Thu Oct  3 06:38:04 2013
///
/// @brief  Diagnostic output selection
///
/// @see N_UTL_FeatureTest.h
/// @see N_IO_CmdParse.h
///
/// Diagnostic output is selected to be written to the log file based on the combination
/// of the FeatureTest and Diag::isActive<type>().  The Diag::isActive(mask) returns true if any
/// of the masked bits are set for the corresponding mask of the mask's type.
///
/// The setXxxDebugLevel() and in the future setXxxDebugMask() functions are used to set
/// the selection mask for the corresponding type.  The setXxxDebugLevel() functions set
/// the bit mask based on the level integer value.  Note that each selection type has a
/// Xxx_DEFAULT mask.  The intent is that when the feature test is active, this mask
/// should be applied as the default output for merely enabling that diagnostic output
/// selection.
///
/// The IO::setXxxDebugLevel() functions are used to implement the command line override
/// of the netlist selected values.
///
/// Recommended use is:
///
///   if (DEBUG_DEVICE && Diag::isActive(Diag::DEVICE_PARAMETER))
///   {
///     .
///     .
///     .
///   }

#ifndef Xyce_N_UTL_Diagnostic_h
#define Xyce_N_UTL_Diagnostic_h


namespace Xyce {
namespace Diag {

//-----------------------------------------------------------------------------
// Class         : TimeIntegrator
// Purpose       : Diagnostic mask type for TimeIntegrator
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:18:31 2015
//-----------------------------------------------------------------------------
///
/// Defines the Time Integrator mask type and values
///
enum TimeIntegrator {
  TIME_PARAMETERS               = 0x001,                        ///< Log parameter settings
  TIME_COEFFICIENTS             = 0x002,                        ///< Log coefficient
  TIME_PREDICTOR                = 0x004,                        ///< Log preditor results
  TIME_RESIDUAL                 = 0x008,                        ///< Log residual calculation results
  TIME_JACOBIAN                 = 0x010,                        ///< Log jacobian calculation results
  TIME_HISTORY                  = 0x020,                        ///< Log time integration history
  TIME_ERROR                    = 0x040,                        ///< Log error weights
  TIME_STEP                     = 0x080,                        ///< Log stepping calculation
  TIME_BREAKPOINTS              = 0x100,                        ///< Log breakpoint values
  TIME_DUMP_SOLUTION_ARRAYS     = 0x200,                        ///< Log solution vectors (should probably write to separate file)
  TIME_OUTPUT                   = 0x400,                        ///< Log ?
  TIME_DEFAULT                  = TIME_PARAMETERS               ///< Default mask when feature test enabled
};

//-----------------------------------------------------------------------------
// Class         : Device
// Purpose       : Diagnostic mask type for devices
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:18:31 2015
//-----------------------------------------------------------------------------
///
/// Defines the device mask type and values
///
enum Device {
  DEVICE_PARAMETERS             = 0x01,
  DEVICE_JACSTAMP               = 0x02,
  DEVICE_LIDS                   = 0x04,
  DEVICE_LOAD_VECTOR            = 0x08,
  DEVICE_PRINT_VECTORS          = 0x10,
  DEVICE_DUMP_VECTORS           = 0x20,
  DEVICE_SOLVER_STATE           = 0x40,
  DEVICE_DEFAULT                = DEVICE_PARAMETERS
};

//-----------------------------------------------------------------------------
// Class         : HB
// Purpose       : Diagnostic mask type for harmonic balance
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:18:31 2015
//-----------------------------------------------------------------------------
///
/// Defines the harmonic balance mask type and values
///
enum HB {
  HB_FAST_TIMES                 = 0x01,
  HB_PRINT_VECTORS              = 0x02,
  HB_TIMESTEP                   = 0x04,
  HB_DEFAULT                    = 0x00
};

//-----------------------------------------------------------------------------
// Class         : Sensitivity
// Purpose       : Diagnostic mask type for sensitivity calculation
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:18:31 2015
//-----------------------------------------------------------------------------
///
/// Defines the sensitivity mask type and values
///
enum Sensitivity {
  SENS_PARAMETERS               = 0x01,
  SENS_SOLVER                   = 0x02,
  SENS_DEFAULT                  = SENS_PARAMETERS
};

//-----------------------------------------------------------------------------
// Class         : MPDE
// Purpose       : Diagnostic mask type for multi-time PDE
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:18:31 2015
//-----------------------------------------------------------------------------
///
/// Defines the multi-time PDE mask type and values
///
enum MPDE {
  MPDE_PARAMETERS               = 0x01,
  MPDE_PRINT_VECTORS            = 0x02,
  MPDE_DEFAULT                  = MPDE_PARAMETERS
};

//-----------------------------------------------------------------------------
// Class         : Nonlinear
// Purpose       : Diagnostic mask type for the nonlinear solver
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:18:31 2015
//-----------------------------------------------------------------------------
///
/// Defines the nonlinear solver mask type and values
///
enum Nonlinear {
  NONLINEAR_PARAMETERS          = 0x01,
  NONLINEAR_CONDUCTANCE         = 0x02,
  NONLINEAR_DUMP                = 0x04,
  NONLINEAR_DUMP_STEP           = 0x08,
  NONLINEAR_DUMP_PARAM_NUMBER   = 0x10,
  NONLINEAR_DUMP_MASK           = NONLINEAR_DUMP | NONLINEAR_DUMP_STEP | NONLINEAR_DUMP_PARAM_NUMBER,
  NONLINEAR_DEFAULT             = NONLINEAR_PARAMETERS
};

//-----------------------------------------------------------------------------
// Class         : IO
// Purpose       : Diagnostic mask type for I/O
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:18:31 2015
//-----------------------------------------------------------------------------
///
/// Defines the I/O mask type and values
///
enum IO {
  IO_PARSE                      = 0x01,
  IO_DEVICE_PARAMETERS          = 0x02,
  IO_DEFAULT                    = IO_PARSE
};

//-----------------------------------------------------------------------------
// Function      : getMask
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:55:06 2015
//-----------------------------------------------------------------------------
///
/// @return Current output selection mask for type specified type.
///
/// @param T    Selection mask type
///
template<class T>
int &getMask();

} // namespace Diag

//-----------------------------------------------------------------------------
// Function      : isActive
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jan 28 11:56:10 2015
//-----------------------------------------------------------------------------
///
/// Logically ands the specified bit mask with the current selection bit mask for the
/// specifies type.
///
/// @param T    Selection mask type
/// @param mask Mask
///
/// @return true if any masked bits are set
///
///
template<class T>
inline bool isActive(const T &mask) {
  return Diag::getMask<T>() & mask;
}

// template<class T>
// inline void setActive(const T &mask) {
//   Diag::getMask<T>() = mask;
// }

void setTimeIntegratorDebugLevel(int level);
void setDeviceDebugLevel(int level);
void setHBDebugLevel(int level);
void setDeviceSensitivityDebugLevel(int level);
void setMPDEDebugLevel(int level);
void setNonlinearDebugLevel(int level);
void setNonlinearConductanceDebugLevel(int level);
void setSensitivityDebugLevel(int level);
void setNonlinearDumpDebugLevel(int level);

} // namespace Xyce

#endif // Xyce_N_UTL_Diagnostic_h

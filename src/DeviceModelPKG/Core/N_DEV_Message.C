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
// Purpose       : Contains the class definition for the N_ERH_ErrorMgr
//                 class.  This class is the error handler for Xyce.
//
// Special Notes : 
//
// Creator       : Eric Keiter, SNL,  Parallel Computational Sciences
//
// Creation Date : 3/15/00
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_Message.h>
#include <N_DEV_Device.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_UTL_Demangle.h>

namespace Xyce {
namespace Device {

// For DeviceBlock: netlistFilename_, parsedLine_[getNumberOfNodes()].lineNumber_
// For OptionBlock: netlistFilename_, parsedLine[0].lineNumber_
// For ParameterBlock: netlistFilename_, parsedLine[0].lineNumber_

// Currently the model and instance headers are the same, so just use DeviceEntity.  If the output wants to be
// different, then copy the DeviceEntity code to make DeviceModel and DeviceInstance.

struct deviceEntityHeader {
  deviceEntityHeader(const DeviceEntity &device_entity)
    : deviceEntity_(device_entity)
  {}

  const DeviceEntity &        deviceEntity_;
};


std::ostream &operator<<(std::ostream &os, const deviceEntityHeader &x)
{
  os << "Device ";
  x.deviceEntity_.printName(os);

  return os;
}

// *****************************************************************
// IMPORTANT DEVELOPER INFORMATION
//
// The warning, error, and fatal error messages here have very
// different meanings downstream, and even more important, the "Fatal0"
// and "Fatal" variants have a special meaning that must be understood
//
// Warnings are messages that are non-fatal.
// Errors are fatal, but the messaging function does not do the exit.
//    When calling the Error messagers, you are responsible for making sure
//    the code exits cleanly.
// Fatal errors cause an immediate abort.  Due to the nature of running in
//    parallel, it is dangerous to use these, and their use should be
//    avoided where possible.  It is better to use "Error" messages and to
//    make sure the code is structured so that after the error message is
//    emitted you can flag an error condition and return.  Unfortunately
//    many parts of Xyce are not set up this way, and Fatal errors are
//    often necessary.
//
// There are two variants for each messager below, a symmetric version
// (with a "0" suffix) and an asymmetric version (without any suffix).
//
// Symmetric messages are those that are emitted identically by all processors
// in a parallel run, generally after having been synchronized at a barrier.
//
// Asymmetric messages are those that are emitted by only one processor.
// ALMOST ALL DEVICE ERRORS ARE ASYMMETRIC.
//
// Because it is necessary to follow proper procedure on exit in parallel,
// and because FATAL errors signal that an immediate exit is necessary,
// When emitting a FATAL error, it is essential to call only the correct
// variant for the type of message  being emitted.  Using DevelFatal0
// (the symmetric variant) to emit  a fatal error from a single processor
// (e.g. asymmetrically) will cause Xyce to hang on exit unless it just happens
// to be the case that the processor emitting the error is proc0.
//
// Similarly, using DevelFatal (the asymmetric variant) to throw an error
// from all processors (i.e. symmetrically) will randomly hang on exit,
// because all processors will be calling MPI_Abort at the same time, tripping
// a race condition in mpiron.
//
// The entire fatal error handling in Xyce is ripe for a refactor.  In the
// meantime, it is incumbent on all developers to use the error messaging
// system properly.


ParamWarning::ParamWarning(const DeviceEntity &device_entity)
  : Report::UserWarning()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

ParamError::ParamError(const DeviceEntity &device_entity)
  : Report::UserError()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

UserInfo::UserInfo(const DeviceEntity &device_entity)
  : Report::UserInfo()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

UserInfo0::UserInfo0(const DeviceEntity &device_entity)
  : Report::UserInfo0()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

UserWarning::UserWarning(const DeviceEntity &device_entity)
  : Report::UserWarning()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

UserWarning0::UserWarning0(const DeviceEntity &device_entity)
  : Report::UserWarning0()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

UserError::UserError(const DeviceEntity &device_entity)
  : Report::UserError()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

UserError0::UserError0(const DeviceEntity &device_entity)
  : Report::UserError0()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

UserFatal::UserFatal(const DeviceEntity &device_entity)
  : Report::UserFatal()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

UserFatal0::UserFatal0(const DeviceEntity &device_entity)
  : Report::UserFatal0()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

DevelFatal::DevelFatal(const DeviceEntity &device_entity, const char *function_name)
  : Report::DevelFatal()
{
  at(device_entity.netlistLocation());
  in(function_name);

  os() << deviceEntityHeader(device_entity) << ": ";
}

DevelFatal0::DevelFatal0(const DeviceEntity &device_entity, const char *function_name)
  : Report::DevelFatal0()
{
  at(device_entity.netlistLocation());
  in(function_name);

  os() << deviceEntityHeader(device_entity) << ": ";
}

struct deviceHeader {
  deviceHeader(const Device &device_)
    : device_(device_)
  {}

  const Device &        device_;
};

std::ostream &operator<<(std::ostream &os, const deviceHeader &x)
{
  os << "Device " << x.device_.getName();

  return os;
}

UserWarning::UserWarning(const Device &device)
  : Report::UserWarning()
{
  os() << deviceHeader(device) << ": ";
}

UserWarning0::UserWarning0(const Device &device)
  : Report::UserWarning0()
{
  os() << deviceHeader(device) << ": ";
}

UserError::UserError(const Device &device)
  : Report::UserError()
{
  os() << deviceHeader(device) << ": ";
}

UserError0::UserError0(const Device &device)
  : Report::UserError0()
{
  os() << deviceHeader(device) << ": ";
}

UserFatal::UserFatal(const Device &device)
  : Report::UserFatal()
{
  os() << deviceHeader(device) << ": ";
}

UserFatal0::UserFatal0(const Device &device)
  : Report::UserFatal0()
{
  os() << deviceHeader(device) << ": ";
}

DevelFatal::DevelFatal(const Device &device, const char *function_name)
  : Report::DevelFatal()
{
  in(function_name);

  os() << deviceHeader(device) << ": ";
}

DevelFatal0::DevelFatal0(const Device &device, const char *function_name)
  : Report::DevelFatal0()
{
  in(function_name);

  os() << deviceHeader(device) << ": ";
}

void device_assertion_error(const DeviceEntity &device_entity, const std::type_info &type, const char *label) 
{
  DevelFatal0(device_entity).in(demangle(type.name()).c_str()) << "Assertion " << label << " failed";
}

} // namespace Device
} // namespace Xyce

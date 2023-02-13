//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Purpose        : Provides templated versions of some boilerplate functions
//                  that are device-specific (so they can't easily be included
//                  in the base device, instance, or model classes).
//
// Special Notes  : Much of the functionality of the device classes, like
//                  N_DEV_Capacitor, is simply to manage STL containers
//                  of model and instance pointers.  That management is pretty
//                  generic, so templating that functionality makes sense.
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 01/31/06
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_DeviceMaster.h>
#include <N_DEV_Message.h>
#include <N_DEV_DeviceInstance.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : instance_must_reference_model_error
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Feb  4 10:16:23 2014
//-----------------------------------------------------------------------------
///
/// reports that the type of instance requires that a model be specified
///
/// @param device                const reference to the device
/// @param model_name            const reference to the model name
/// @param netlist_filename          const reference to the netlist path
/// @param line_number           line number in the netlist path
///
void instance_must_reference_model_error(const Device &device, const std::string &model_name, const NetlistLocation &netlist_location)
{
  UserError(device).at(netlist_location) << model_name << " instance must reference a model";
}

//-----------------------------------------------------------------------------
// Function      : could_not_find_model_error
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Feb  4 10:18:05 2014
//-----------------------------------------------------------------------------
///
/// reports that the model name is note defined
///
/// @param device                const reference to the device
/// @param model_name            const reference to the model name
/// @param instance_name         const reference to the instance name
/// @param netlist_filename          const reference to the netlist path
/// @param line_number           line number in the netlist path
///
void could_not_find_model_error(const Device &device, const std::string &model_name, const std::string &instance_name, const NetlistLocation &netlist_location)
{
  UserError(device).at(netlist_location) << "Could not find model " << model_name << " which is referenced by instance " << instance_name;
}

//-----------------------------------------------------------------------------
// Function      : duplicate_model_warning
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Feb  4 10:18:41 2014
//-----------------------------------------------------------------------------
///
/// duplicate_model_warning reports that the model name is duplicated.
///
///
/// @param device                const reference to the device
/// @param model_name            const reference to the model name
/// @param netlist_filename          const reference to the netlist path
/// @param line_number           line number in the netlist path
///
void duplicate_model_warning(const Device &device, const DeviceModel &model, const NetlistLocation &netlist_location)
{
  UserWarning message(device);
  message.at(netlist_location) << "Attempted to add model ";
  model.printName(message.os());
  message << " that already exists, ignoring all but the first definition";
}


//-----------------------------------------------------------------------------
// Function      : duplicate_instance_warninga
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Feb  4 10:12:14 2014
//-----------------------------------------------------------------------------
///
/// duplicate_instance_warning reports a duplication of instance names.
///
/// Currently models and devices can share a name and the current implementation of instanceMap_ results in lost
/// information.
///
/// @param instance        const reference to the instance that is being added
///
void duplicate_instance_warning(const Device &device, const DeviceInstance &instance, const NetlistLocation &netlist_location)
{
  UserWarning message(device);
  message.at(netlist_location) << "Attempted to add instance ";
  instance.printName(message.os());
  message << " that already exists, ignoring all but the first definition";
}


//-----------------------------------------------------------------------------
// Function      : duplicate_entity_warninga
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Feb  4 10:12:14 2014
//-----------------------------------------------------------------------------
///
/// duplicate_entity_warning reports a duplication of entity names.
///
/// Currently models and devices can share a name and the current implementation of entityMap_ results in lost
/// information.
///
/// @param entity        const reference to the entity that is being added
///
void duplicate_entity_warning(const Device &device, const DeviceEntity &entity, const NetlistLocation &netlist_location)
{
  UserWarning message(device);
  message.at(netlist_location) << "Duplicated model and device name ";
  entity.printName(message.os());
}

} // namespace Device
} // namespace Xyce

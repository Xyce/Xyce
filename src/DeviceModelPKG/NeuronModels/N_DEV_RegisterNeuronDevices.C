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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : David Baur
//
// Creation Date  : 3/15/2013
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_RegisterDevices.h>

#include <N_DEV_Neuron.h>
#include <N_DEV_Neuron2.h>
#include <N_DEV_Neuron3.h>
#include <N_DEV_Neuron4.h>
#include <N_DEV_Neuron5.h>
#include <N_DEV_Neuron6.h>
#include <N_DEV_Neuron7.h>
#include <N_DEV_Neuron8.h>
#include <N_DEV_Neuron9.h>
#include <N_DEV_NeuronPop1.h>
#include <N_DEV_Synapse.h>
#include <N_DEV_Synapse2.h>
#include <N_DEV_Synapse3.h>
#include <N_DEV_Synapse4.h>

namespace Xyce {
namespace Device {

void 
registerNeuronDevices(
  const DeviceCountMap& deviceMap, 
  const std::set<int>& levelSet)
{
  Neuron::registerDevice(deviceMap, levelSet);
  Neuron2::registerDevice(deviceMap, levelSet);
  Neuron3::registerDevice(deviceMap, levelSet);
  Neuron4::registerDevice(deviceMap, levelSet);
  Neuron5::registerDevice(deviceMap, levelSet);
  Neuron6::registerDevice(deviceMap, levelSet);
  Neuron7::registerDevice(deviceMap, levelSet);
  Neuron8::registerDevice(deviceMap, levelSet);
  Neuron9::registerDevice(deviceMap, levelSet);
  NeuronPop1::registerDevice(deviceMap, levelSet);
  Synapse::registerDevice(deviceMap, levelSet);
  Synapse2::registerDevice(deviceMap, levelSet);
  Synapse3::registerDevice(deviceMap, levelSet);
  Synapse4::registerDevice(deviceMap, levelSet);
}

} // namespace Device
} // namespace Xyce

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

#ifdef HAVE_DLFCN_H
#include <dlfcn.h>
#endif

#include <N_DEV_RegisterDevices.h>
#include <N_DEV_RegisterOpenDevices.h>

#include <N_DEV_ACC.h>
#include <N_DEV_ADC.h>
#include <N_DEV_BJT.h>
#include <N_DEV_Bsrc.h>
#include <N_DEV_Battery.h>
#include <N_DEV_Capacitor.h>
#include <N_DEV_DAC.h>
#include <N_DEV_Digital.h>
#include <N_DEV_Diode.h>
#include <N_DEV_GeneralExternal.h>
#include <N_DEV_Inductor.h>
#include <N_DEV_ISRC.h>
#include <N_DEV_JFET.h>
#include <N_DEV_LTRA.h>
#include <N_DEV_MemristorPEM.h>
#include <N_DEV_MemristorTEAM.h>
#include <N_DEV_MemristorYakopcic.h>
#include <N_DEV_MESFET.h>
#include <N_DEV_MOSFET1.h>
#include <N_DEV_MOSFET2.h>
#include <N_DEV_MOSFET3.h>
#include <N_DEV_MOSFET6.h>
#include <N_DEV_MOSFET_B3.h>
#include <N_DEV_MOSFET_B3SOI.h>
#include <N_DEV_MOSFET_B4.h>
#include <N_DEV_MutIndLin.h>
#include <N_DEV_MutIndNonLin2.h>
#include <N_DEV_MutIndNonLin.h>
#include <N_DEV_OpAmp.h>
#include <N_DEV_Resistor3.h>
#include <N_DEV_Resistor.h>
#include <N_DEV_ThermalResistor.h>
#include <N_DEV_ROM.h>
#include <N_DEV_RxnSet.h>
#include <N_DEV_SW.h>
#include <N_DEV_TRA.h>
#include <N_DEV_VCCS.h>
#include <N_DEV_Vcvs.h>
#include <N_DEV_VDMOS.h>
#include <N_DEV_Vsrc.h>
#include <N_DEV_Xygra.h>
#include <N_DEV_TransLine.h>
#include <N_DEV_PowerGrid.h>
#include <N_DEV_PowerGridBranch.h>
#include <N_DEV_PowerGridTransformer.h>
#include <N_DEV_PowerGridBusShunt.h>
#include <N_DEV_PowerGridGenBus.h>
#include <N_DEV_AntiWindupLimiter.h>

namespace Xyce {
namespace Device {

void 
registerOpenDevices(const DeviceCountMap& deviceMap,
                    const std::set<int>& levelSet,
                    bool includeMI)
{
  Resistor::registerDevice(deviceMap, levelSet);
  ThermalResistor::registerDevice(deviceMap, levelSet);
  Resistor3::registerDevice(deviceMap, levelSet);
  Capacitor::registerDevice(deviceMap, levelSet);
  Diode::registerDevice(deviceMap, levelSet);
  BJT::registerDevice(deviceMap, levelSet); 
  JFET::registerDevice(deviceMap, levelSet);
  MemristorPEM::registerDevice(deviceMap, levelSet);
  MemristorTEAM::registerDevice(deviceMap, levelSet);
  MemristorYakopcic::registerDevice(deviceMap, levelSet);
  MESFET::registerDevice(deviceMap, levelSet);
  MOSFET1::registerDevice(deviceMap, levelSet);  
  MOSFET2::registerDevice(deviceMap, levelSet);
  MOSFET3::registerDevice(deviceMap, levelSet);
  MOSFET6::registerDevice(deviceMap, levelSet);
  MOSFET_B3::registerDevice(deviceMap, levelSet);
  MOSFET_B3SOI::registerDevice(deviceMap, levelSet);
  MOSFET_B4::registerDevice(deviceMap, levelSet);
  ROM::registerDevice(deviceMap, levelSet);
  VDMOS::registerDevice(deviceMap, levelSet);
  ISRC::registerDevice(deviceMap, levelSet);
  Vcvs::registerDevice(deviceMap, levelSet);
  Bsrc::registerDevice(deviceMap, levelSet);
  VCCS::registerDevice(deviceMap, levelSet);
  Vsrc::registerDevice(deviceMap, levelSet);
  LTRA::registerDevice(deviceMap, levelSet);
  TRA::registerDevice(deviceMap, levelSet);
  SW::registerDevice(deviceMap, levelSet);
  ADC::registerDevice(deviceMap, levelSet);
  DAC::registerDevice(deviceMap, levelSet);
  OpAmp::registerDevice(deviceMap, levelSet);
  Digital::registerDevice(deviceMap, levelSet);
  ACC::registerDevice(deviceMap, levelSet);
  RxnSet::registerDevice(deviceMap, levelSet);
  TransLine::registerDevice(deviceMap, levelSet);
  PowerGrid::registerDevice(deviceMap, levelSet);
  PowerGridBranch::registerDevice(deviceMap, levelSet);
  PowerGridTransformer::registerDevice(deviceMap, levelSet);
  PowerGridBusShunt::registerDevice(deviceMap, levelSet);
  PowerGridGenBus::registerDevice(deviceMap, levelSet);
  AntiWindupLimiter::registerDevice(deviceMap, levelSet);
  Battery::registerDevice(deviceMap, levelSet);
  Xygra::registerDevice(deviceMap, levelSet);
  GeneralExternal::registerDevice(deviceMap, levelSet);

  if (includeMI)
    registerMutualInductors();
}

void
registerMutualInductors()
{
  MutIndLin::registerDevice();
  MutIndNonLin::registerDevice();
  MutIndNonLin2::registerDevice();
  Inductor::registerDevice();
}

} // namespace Device
} // namespace Xyce

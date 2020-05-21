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

#include <N_DEV_ADMSHBT_X.h>
#include <N_DEV_ADMSPSP103VA.h>
#include <N_DEV_ADMSPSP103TVA.h>
#include <N_DEV_ADMSPSP102VA.h>
#include <N_DEV_ADMSJUNCAP200.h>
#include <N_DEV_ADMSvbic13.h>
#include <N_DEV_ADMSvbic13_4t.h>
#include <N_DEV_ADMSbsimcmg.h>
#include <N_DEV_ADMSbsimcmg_108.h>
#include <N_DEV_ADMSbsimcmg_110.h>
#include <N_DEV_ADMSbsim6.h>
#include <N_DEV_ADMSbsimsoi.h>
#include <N_DEV_ADMSbsimsoi450.h>
#include <N_DEV_ADMSbjt504va.h>
#include <N_DEV_ADMSbjt504tva.h>
#include <N_DEV_ADMSmvs_2_0_0_etsoi.h>
#include <N_DEV_ADMSmvs_2_0_0_hemt.h>
#include <N_DEV_ADMSmvsg_cmc.h>
#include <N_DEV_ADMShicumL2va.h>
#include <N_DEV_ADMShic0_full.h>
#include <N_DEV_MOSFET1.h>
#include <N_DEV_BJT.h>

namespace Xyce {
namespace Device {

void 
registerADMSDevices(const DeviceCountMap& deviceMap,
                    const std::set<int>& levelSet)
{
  BJT::registerDevice();
  MOSFET1::registerDevice();
  Diode::registerDevice();
  ADMSvbic13::registerDevice(deviceMap, levelSet);
  ADMSvbic13_4t::registerDevice(deviceMap, levelSet);
  ADMSHBT_X::registerDevice(deviceMap, levelSet);
  ADMSPSP103VA::registerDevice(deviceMap, levelSet);
  ADMSPSP103TVA::registerDevice(deviceMap, levelSet);
  ADMSPSP102VA::registerDevice(deviceMap, levelSet);
  ADMSJUNCAP200::registerDevice(deviceMap, levelSet);
  ADMSbsimcmg::registerDevice(deviceMap, levelSet);
  ADMSbsimcmg_108::registerDevice(deviceMap,levelSet);
  ADMSbsimcmg_110::registerDevice(deviceMap,levelSet);
  ADMSbsim6::registerDevice(deviceMap, levelSet);
  ADMSbsimsoi::registerDevice(deviceMap, levelSet);
  ADMSbsimsoi450::registerDevice(deviceMap, levelSet);
  ADMSbjt504va::registerDevice(deviceMap, levelSet);
  ADMSbjt504tva::registerDevice(deviceMap, levelSet);
  ADMSmvs_2_0_0_etsoi::registerDevice(deviceMap, levelSet);
  ADMSmvs_2_0_0_hemt::registerDevice(deviceMap, levelSet);
  ADMSmvsg_cmc::registerDevice(deviceMap, levelSet);
  ADMShicumL2va::registerDevice(deviceMap, levelSet);
  ADMShic0_full::registerDevice(deviceMap, levelSet);
}

} // namespace Device
} // namespace Xyce

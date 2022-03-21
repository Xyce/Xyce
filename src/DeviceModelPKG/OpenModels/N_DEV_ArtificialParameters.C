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

//-------------------------------------------------------------------------
//
// Purpose        : 
//
// Special Notes  :
//
// Creator        : Dave Baur
//
// Creation Date  : 
//
//
//
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_DEV_ArtificialParameters.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceMgr.h>
#include <N_ERH_Message.h>

#include <N_DEV_MOSFET1.h>
#include <N_DEV_BJT.h>
#include <N_DEV_Diode.h>
#include <N_DEV_Vsrc.h>
#include <N_DEV_ISRC.h>

namespace Xyce {
namespace Device {
namespace ArtificialParameters {

SolverState &
ArtificialParameter::getSolverState(DeviceMgr &device_manager)
{
  return device_manager.solState_;
}

const SolverState &
ArtificialParameter::getSolverState(const DeviceMgr &device_manager) const
{
  return device_manager.solState_;
}

DeviceOptions &
ArtificialParameter::getDeviceOptions(DeviceMgr &device_manager)
{
  return device_manager.devOptions_;
}

const DeviceOptions &
ArtificialParameter::getDeviceOptions(const DeviceMgr &device_manager) const
{
  return device_manager.devOptions_;
}


ModelTypeInstanceVectorMap &
ArtificialParameter::getModelTypeInstanceVectorMap(DeviceMgr &device_manager)
{
  return device_manager.modelGroupInstanceVector_;
}

const ModelTypeInstanceVectorMap &
ArtificialParameter::getModelTypeInstanceVectorMap(const DeviceMgr &device_manager) const
{
  return device_manager.modelGroupInstanceVector_;
}

InstanceVector &
ArtificialParameter::getInstanceVector(DeviceMgr &device_manager)
{
  return device_manager.instancePtrVec_;
}

const InstanceVector &
ArtificialParameter::getInstanceVector(const DeviceMgr &device_manager) const
{
  return device_manager.instancePtrVec_;
}

bool MOSFETGainScaleParam::setValue(DeviceMgr &device_manager, double value)
{
  getSolverState(device_manager).gainScale_ = value;
  getSolverState(device_manager).artParameterFlag_ = true;
  return true;
}

double MOSFETGainScaleParam::getValue(const DeviceMgr &device_manager) const
{
  return getSolverState(device_manager).gainScale_;
}

bool MOSFETNLTermScaleParam::setValue(DeviceMgr &device_manager, double value)
{
  getSolverState(device_manager).nltermScale_ = value;
  getSolverState(device_manager).artParameterFlag_ = true;
  return true;
}

double MOSFETNLTermScaleParam::getValue(const DeviceMgr &device_manager) const
{
  return getSolverState(device_manager).nltermScale_;
}

bool MOSFETLParam::setValue(DeviceMgr &device_manager, double value)
{
  bool success = true;

  getSolverState(device_manager).sizeParameterFlag_ = true;

  // now loop over all the mosfet instances, and change the l.
  DeviceMgr::ModelTypeInstanceVectorMap::const_iterator model_group_it = getModelTypeInstanceVectorMap(device_manager).find(MOSFET1::Traits::modelGroup());
  if (model_group_it != getModelTypeInstanceVectorMap(device_manager).end()) {
    for (InstanceVector::const_iterator it = (*model_group_it).second.begin(); it != (*model_group_it).second.end(); ++it)
    {
      success = (*it)->setParam("l", value);
      success = success && (*it)->processParams ();
    }
  }

  // getDeviceOptions(device_manager).defl = value; // Shouldn't this be set?

  return success;
}

double MOSFETLParam::getValue(const DeviceMgr &device_manager) const
{
  return getDeviceOptions(device_manager).defl;
}

bool MOSFETWParam::setValue(DeviceMgr &device_manager, double value)
{
  bool success = true;

  getSolverState(device_manager).sizeParameterFlag_ = true;

  // now loop over all the mosfet instances, and change the w.
  DeviceMgr::ModelTypeInstanceVectorMap::const_iterator model_group_it = getModelTypeInstanceVectorMap(device_manager).find(MOSFET1::Traits::modelGroup());
  if (model_group_it != getModelTypeInstanceVectorMap(device_manager).end()) {
    for (InstanceVector::const_iterator it = (*model_group_it).second.begin(); it != (*model_group_it).second.end(); ++it)
    {
      success = (*it)->setParam("w", value);
      success = success && (*it)->processParams ();
    }
  }

  // getDeviceOptions(device_manager).defw = value; // Shouldn't this be set?

  return success;
}

double MOSFETWParam::getValue(const DeviceMgr &device_manager) const
{
  return getDeviceOptions(device_manager).defw;
}

bool MOSFETSizeScaleParam::setValue(DeviceMgr &device_manager, double value)
{
  bool success = true;

  getSolverState(device_manager).sizeParameterFlag_ = true;
  getSolverState(device_manager).sizeScale_ = value;

  // What we want at val = 0.0.  For all the lengths and widths to have a
  // ratio of:  L=5u W=175u  or L=5u W=270u.  Compromise.  L/W = 5/200,
  // with L actually being 5u.

  double length0 = getDeviceOptions(device_manager).length0;
  double width0  = getDeviceOptions(device_manager).width0;

  DeviceMgr::ModelTypeInstanceVectorMap::const_iterator model_group_it = getModelTypeInstanceVectorMap(device_manager).find(MOSFET1::Traits::modelGroup());
  if (model_group_it != getModelTypeInstanceVectorMap(device_manager).end()) {
    for (InstanceVector::const_iterator it = (*model_group_it).second.begin(); it != (*model_group_it).second.end(); ++it)
    {
      success = (*it)->scaleParam("l", value, length0);
      success = success || (*it)->scaleParam("w", value, width0);
      success = success && (*it)->processParams();
    }
  }

  return success;
}

double MOSFETSizeScaleParam::getValue(const DeviceMgr &device_manager) const
{
  return 0.0;
}

bool MOSFETTOXParam::setValue(DeviceMgr &device_manager, double value)
{
  bool success = true;

  getSolverState(device_manager).sizeParameterFlag_ = true;
  getSolverState(device_manager).sizeScale_ = value;

  // What we want at val = 0.0.  For all the lengths and widths to have a
  // ratio of:  L=5u W=175u  or L=5u W=270u.  Compromise.  L/W = 5/200,
  // with L actually being 5u.

  double tox0    = getDeviceOptions(device_manager).tox0;

  DeviceMgr::ModelTypeInstanceVectorMap::const_iterator model_group_it = getModelTypeInstanceVectorMap(device_manager).find(MOSFET1::Traits::modelGroup());
  if (model_group_it != getModelTypeInstanceVectorMap(device_manager).end()) {
    for (InstanceVector::const_iterator it = (*model_group_it).second.begin(); it != (*model_group_it).second.end(); ++it)
    {
      success = (*it)->scaleParam("tox", value, tox0);
      success = success && (*it)->processParams();
      success = success && (*it)->processInstanceParams();
    }
  }

  return success;
}

double MOSFETTOXParam::getValue(const DeviceMgr &device_manager) const
{
  return 0.0;
}

bool BJTBFParam::setValue(DeviceMgr &device_manager, double value)
{
  bool success = true;

  DeviceMgr::ModelTypeInstanceVectorMap::const_iterator model_group_it = getModelTypeInstanceVectorMap(device_manager).find(BJT::Traits::modelGroup());
  if (model_group_it != getModelTypeInstanceVectorMap(device_manager).end()) {
    for (InstanceVector::const_iterator it = (*model_group_it).second.begin(); it != (*model_group_it).second.end(); ++it)
    {
      success = (*it)->scaleParam("bf", value, 0.0);
      success = success && (*it)->processParams();
      success = success && (*it)->processInstanceParams();
    }
  }

  return success;
}

double BJTBFParam::getValue(const DeviceMgr &device_manager) const
{
  return 0.0;
}

bool BJTNFParam::setValue(DeviceMgr &device_manager, double value)
{
  bool success = true;

  DeviceMgr::ModelTypeInstanceVectorMap::const_iterator model_group_it = getModelTypeInstanceVectorMap(device_manager).find(BJT::Traits::modelGroup());
  if (model_group_it != getModelTypeInstanceVectorMap(device_manager).end()) {
    for (InstanceVector::const_iterator it = (*model_group_it).second.begin(); it != (*model_group_it).second.end(); ++it)
    {
      success = (*it)->scaleParam("nf", value, 0.0);
      success = success && (*it)->processParams();
      success = success && (*it)->processInstanceParams();
    }
  }

  return success;
}

double BJTNFParam::getValue(const DeviceMgr &device_manager) const
{
  return 0.0;
}

bool BJTNRParam::setValue(DeviceMgr &device_manager, double value)
{
  bool success = true;

  DeviceMgr::ModelTypeInstanceVectorMap::const_iterator model_group_it = getModelTypeInstanceVectorMap(device_manager).find(BJT::Traits::modelGroup());
  if (model_group_it != getModelTypeInstanceVectorMap(device_manager).end()) {
    for (InstanceVector::const_iterator it = (*model_group_it).second.begin(); it != (*model_group_it).second.end(); ++it)
    {
      success = (*it)->scaleParam("nr", value, 0.0);
      success = success && (*it)->processParams();
      success = success && (*it)->processInstanceParams();
    }
  }

  return success;
}

double BJTNRParam::getValue(const DeviceMgr &device_manager) const
{
  return 0.0;
}

bool DiodeNParam::setValue(DeviceMgr &device_manager, double value)
{
  bool success = true;

  DeviceMgr::ModelTypeInstanceVectorMap::const_iterator model_group_it = getModelTypeInstanceVectorMap(device_manager).find(Diode::Traits::modelGroup());
  if (model_group_it != getModelTypeInstanceVectorMap(device_manager).end()) {
    for (InstanceVector::const_iterator it = (*model_group_it).second.begin(); it != (*model_group_it).second.end(); ++it)
    {
      success = (*it)->scaleParam("n", value, 10.0);
      success = success && (*it)->processParams();
      success = success && (*it)->processInstanceParams();
    }
  }

  return success;
}

double DiodeNParam::getValue(const DeviceMgr &device_manager) const
{
  return getDeviceOptions(device_manager).defw;
}

bool VsrcScaleParam::setValue(DeviceMgr &device_manager, double value)
{
  bool success = true;

  {
    DeviceMgr::ModelTypeInstanceVectorMap::const_iterator model_group_it = getModelTypeInstanceVectorMap(device_manager).find(Vsrc::Traits::modelGroup());
    if (model_group_it != getModelTypeInstanceVectorMap(device_manager).end()) {
      for (InstanceVector::const_iterator it = (*model_group_it).second.begin(); it != (*model_group_it).second.end(); ++it)
      {
        success = (*it)->scaleDefaultParam(value);
        success = (*it)->processParams();
      }
    }
  }

  {
    DeviceMgr::ModelTypeInstanceVectorMap::const_iterator model_group_it = getModelTypeInstanceVectorMap(device_manager).find(ISRC::Traits::modelGroup());
    if (model_group_it != getModelTypeInstanceVectorMap(device_manager).end()) {
      for (InstanceVector::const_iterator it = (*model_group_it).second.begin(); it != (*model_group_it).second.end(); ++it)
      {
        success = (*it)->scaleDefaultParam(value);
        success = (*it)->processParams();
      }
    }
  }

  return success;
}

double VsrcScaleParam::getValue(const DeviceMgr &device_manager) const
{
  return getDeviceOptions(device_manager).defw;
}

bool PDEAlphaParam::setValue(DeviceMgr &device_manager, double value)
{
  getSolverState(device_manager).pdeAlpha_ = value; // not important - part of planned refactor.

  if (!getSolverState(device_manager).PDEcontinuationFlag_)
  {
    Report::DevelFatal().in("DeviceMgr::setParam") << "Tried to set pdeAlpha without first calling enablePDEContinuation";
  }

  for (InstanceVector::const_iterator it = getInstanceVector(device_manager).begin(), end = getInstanceVector(device_manager).end(); it != end; ++it)
  {
    (*it)->setPDEContinuationAlpha(value);
  }

  return true;
}

double PDEAlphaParam::getValue(const DeviceMgr &device_manager) const
{
  return getSolverState(device_manager).pdeAlpha_;
}

bool PDEBetaParam::setValue(DeviceMgr &device_manager, double value)
{
  getSolverState(device_manager).PDEcontinuationFlag_ = true;

  for (InstanceVector::const_iterator it = getInstanceVector(device_manager).begin(), end = getInstanceVector(device_manager).end(); it != end; ++it)
  {
    (*it)->setPDEContinuationBeta(value);
  }

  return true;
}

double PDEBetaParam::getValue(const DeviceMgr &device_manager) const
{
  return 0.0;
}

bool PDEChargeAlphaParam::setValue(DeviceMgr &device_manager, double value)
{
  getSolverState(device_manager).chargeAlpha_ = value;
  getSolverState(device_manager).chargeHomotopy_ = true;

  return true;
}

double PDEChargeAlphaParam::getValue(const DeviceMgr &device_manager) const
{
  return getSolverState(device_manager).chargeAlpha_;
}

bool GSteppingParam::setValue(DeviceMgr &device_manager, double value)
{
  return true;
}

double GSteppingParam::getValue(const DeviceMgr &device_manager) const
{
  return getSolverState(device_manager).chargeAlpha_;
}

bool GMinParam::setValue(DeviceMgr &device_manager, double value)
{
  getDeviceOptions(device_manager).gmin = getDeviceOptions(device_manager).gmin_orig*value + getDeviceOptions(device_manager).gmin_init*(1.0 - value);

  return true;
}

double GMinParam::getValue(const DeviceMgr &device_manager) const
{
  return getDeviceOptions(device_manager).gmin;
}

bool VtParam::setValue(DeviceMgr &device_manager, double value)
{
  return true;
}

double VtParam::getValue(const DeviceMgr &device_manager) const
{
  return getDeviceOptions(device_manager).temp.getImmutableValue<double>()*CONSTKoverQ;
}

bool TempParam::setValue(DeviceMgr &device_manager, double value)
{
  device_manager.updateTemperature(value);

  return true;
}

double TempParam::getValue(const DeviceMgr &device_manager) const
{
  return getDeviceOptions(device_manager).temp.getImmutableValue<double>() - CONSTCtoK;
}

} // namespace ArtificialParameters
} // namespace Device
} // namespace Xyce

//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Filename      : N_IO_CircuitMetadata.C
//
// Purpose       :
//
// Special Notes :
//
// Creator       :
//
// Creation Date :
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <N_IO_CircuitMetadata.h>

#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_SourceData.h>
#include <N_ERH_Message.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_MachDepParams.h>

namespace Xyce {
namespace IO {

//----------------------------------------------------------------------------
// Function       : CircuitMetadata::CircuitMetadata
// Purpose        : constructor
// Special Notes  :
// Scope          : public
// Creator        : 
// Creation Date  : 
//----------------------------------------------------------------------------
CircuitMetadata::CircuitMetadata()
  : deviceMetadata_(),
    deviceMetadataIndex()
{}

//----------------------------------------------------------------------------
// Function       : CircuitMetadata::getDeviceMetadata
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 09/28/2003
//----------------------------------------------------------------------------
const DeviceMetadata &
CircuitMetadata::getDeviceMetadata(
  const std::string &   deviceTypeIn,
  int                   level) const
{
  if (level == -1)
  {
    level = 1;
  }

  std::string deviceType = deviceTypeIn;
  if (deviceTypeIn == "K")
  {
    deviceType = "L";
  }

  DeviceMetadataIndexMap::const_iterator it_device_meta_index = 
    deviceMetadataIndex.find(NameLevelKey(deviceType, level));

  if (it_device_meta_index != deviceMetadataIndex.end())
  {
    return deviceMetadata_[it_device_meta_index->second];
  }

  // Handle default model:
  DeviceMetadata &device_metadata = deviceMetadata_[NameLevelKey(deviceType, level)];

  const Device::Configuration *configuration = Device::Configuration::findConfiguration(deviceType, level);
  if (configuration) {
    device_metadata.configuration = configuration;

    device_metadata.levelValid = true;
    device_metadata.numOptionalNodes = configuration->getNumOptionalNodes();
    device_metadata.numFillNodes = configuration->getNumFillNodes();
    device_metadata.numNodes = configuration->getNumNodes();
    device_metadata.modelRequired = configuration->getModelRequired();
    device_metadata.primaryParameter = configuration->getPrimaryParameter();
    device_metadata.modelTypes = configuration->getModelTypeNames();

    Xyce::Device::populateParams(configuration->getModelParameters().getMap(),
        device_metadata.modelParameters, device_metadata.modelCompositeParameterMap);
    Xyce::Device::populateParams(configuration->getInstanceParameters().getMap(),
        device_metadata.instanceParameters, device_metadata.instanceCompositeParameterMap);
  }

  for (std::vector<std::string>::iterator it_device_model_type = device_metadata.modelTypes.begin(); 
      it_device_model_type != device_metadata.modelTypes.end(); ++it_device_model_type)
  {
    deviceMetadataIndex.insert(DeviceMetadataIndexMap::value_type
        (NameLevelKey(*it_device_model_type, level), NameLevelKey(deviceType, level)));
  }
  deviceMetadataIndex[NameLevelKey(deviceType, level)] = NameLevelKey(deviceType, level);
  if (deviceType == "L")
    deviceMetadataIndex[NameLevelKey(std::string("K"), level)] = NameLevelKey(deviceType, level);

  return device_metadata;
}

//-----------------------------------------------------------------------------
// Function      : CircuitMetadata::getInstanceParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 09/30/2003
//-----------------------------------------------------------------------------
const std::vector<Device::Param> &
CircuitMetadata::getInstanceParameters(
  const std::string &   deviceType,
  int                   modelLevel) const
{
  if (modelLevel == -1 || deviceMetadataIndex.find(NameLevelKey(deviceType, modelLevel)) == deviceMetadataIndex.end())
  {
    return getDeviceMetadata(deviceType, modelLevel).instanceParameters;
  }
  else
  {
    return deviceMetadata_[deviceMetadataIndex[NameLevelKey(deviceType, modelLevel)]].instanceParameters;
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitMetadata::getInstanceCompositeComponents
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 10/07/2003
//----------------------------------------------------------------------------
void
CircuitMetadata::getInstanceCompositeComponents(
  const std::string &           deviceType,
  const std::string &           parameterName,
  int                           modelLevel,
  std::vector<Device::Param> &  components) const
{
  const DeviceMetadata &device_metadata = getDeviceMetadata(deviceType, modelLevel);
  const DeviceParamMap &icpMap = device_metadata.instanceCompositeParameterMap;
  const DeviceParamMap::const_iterator iterIcp = icpMap.find(parameterName);

  if ( iterIcp != icpMap.end() )
  {
    components = iterIcp->second;
  }
  else
  {
    Report::UserError() 
      << "There are no component parameters in metadata for the VECTOR-COMPOSITE parameter: " 
      << parameterName;
  }
}

//-----------------------------------------------------------------------------
// Function      : CircuitMetadata::getModelParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 09/30/2003
//-----------------------------------------------------------------------------
const std::vector<Device::Param> &CircuitMetadata::getModelParameters(
  const std::string &   modelType,
  int                   modelLevel) const
{
  if (modelLevel == -1 || deviceMetadataIndex.find(NameLevelKey(modelType, modelLevel)) == deviceMetadataIndex.end())
  {
    return getDeviceMetadata(modelType, modelLevel).modelParameters;
  }
  else
  {
    return deviceMetadata_[deviceMetadataIndex[NameLevelKey(modelType, modelLevel)]].modelParameters;
  }
}


//----------------------------------------------------------------------------
// Function       : CircuitMetadata::getModelCompositeComponents
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Eric Keiter
// Creation Date  : 5/07/2008
//----------------------------------------------------------------------------
void CircuitMetadata::getModelCompositeComponents(
  const std::string &                modelType,
  const std::string &                parameterName, int modelLevel,
  std::vector<Device::Param> &       components) const
{
  const DeviceMetadata &device_metadata = getDeviceMetadata(modelType, modelLevel);
  const DeviceParamMap &mcpMap = device_metadata.modelCompositeParameterMap;

  DeviceParamMap::const_iterator iterIcp = mcpMap.find(parameterName);

  if ( iterIcp != mcpMap.end() )
  {
    components = iterIcp->second;
  }
  else
  {
    Report::UserError() 
      << "There are no component parameters in metadata for the VECTOR-COMPOSITE parameter " 
      << parameterName;
  }
}

} // namespace IO
} // namespace Xyce

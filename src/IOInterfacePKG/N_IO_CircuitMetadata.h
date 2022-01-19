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
// Purpose       :
//
// Special Notes :
//
// Creator       : Lon Waters
//
// Creation Date : 12/11/00
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_IO_CircuitMetadata_h
#define Xyce_N_IO_CircuitMetadata_h

#include <map>
#include <string>
#include <vector>

#include <N_DEV_fwd.h>
#include <N_DEV_Param.h>
#include <N_UTL_NoCase.h>
#include <N_UTL_OpBuilder.h>

namespace Xyce {
namespace IO {

typedef std::map<std::string, std::vector<Device::Param>, LessNoCase> DeviceParamMap;

//-----------------------------------------------------------------------------
// Class         : DeviceMetadata
// Purpose       :
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 12/11/00
//-----------------------------------------------------------------------------

struct DeviceMetadata
{
  DeviceMetadata()
    : configuration(0),
      levelValid(false),
      deviceType(),
      level(0),
      numNodes(0),
      numOptionalNodes(0),
      numFillNodes(0),
      modelRequired(0),
      primaryParameter(),
      modelTypes(),
      instanceParameters(),
      modelParameters(),
      instanceCompositeParameterMap(),
      modelCompositeParameterMap()
  {}

  bool isModelTypeValid(const std::string & modelType) const 
  {
    return std::find_if(modelTypes.begin(), modelTypes.end(), EqualNoCasePred(modelType)) != modelTypes.end();
  }

  bool isModelLevelValid() const 
  {
    return levelValid;
  }

    const Device::Configuration *       configuration;

  bool                        levelValid;
  std::string                 deviceType;
  int                         level;
  int                         numNodes;
  int                         numOptionalNodes;
  int                         numFillNodes;
  int                         modelRequired;
  std::string                 primaryParameter;
  std::vector<std::string>    modelTypes;
  std::vector<Device::Param>  instanceParameters;
  std::vector<Device::Param>  modelParameters;
  DeviceParamMap              instanceCompositeParameterMap;
  DeviceParamMap              modelCompositeParameterMap;
};

//-----------------------------------------------------------------------------
// Class         : CircuitMetadata
// Purpose       :
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 12/11/00
//-----------------------------------------------------------------------------
class CircuitMetadata
{
public:
  typedef std::map<NameLevelKey, NameLevelKey, NameLevelKeyLess> DeviceMetadataIndexMap;
  typedef std::map<NameLevelKey, DeviceMetadata, NameLevelKeyLess> DeviceMetadataMap;

  CircuitMetadata();

  ~CircuitMetadata()
  {}

  CircuitMetadata(const CircuitMetadata &);
  CircuitMetadata &operator=(const CircuitMetadata &);

  // Find the metadata for a given device. If not found, return NULL.
  const DeviceMetadata &getDeviceMetadata(const std::string & deviceType, int level) const;

  // Determine if the given model type is valid for the given device.
  bool isModelLevelValid(const std::string & deviceType, int modelLevel) const 
  {
    return getDeviceMetadata(deviceType, modelLevel).isModelLevelValid();
  }

  // Determine if the given model type is valid for the given device.
  bool isModelTypeValid(const std::string & deviceType, const std::string & modelType, int modelLevel) const 
  {
    return getDeviceMetadata(deviceType, modelLevel).isModelTypeValid(modelType);
  }

  // Determine if the given parameter name is a valid parameter for the given device.
  bool isDeviceParameter(const std::string & deviceType, int modelLevel, const std::string & parameterName) const 
  {
    const DeviceMetadata &device_metadata = getDeviceMetadata(deviceType, modelLevel);

    std::vector<Device::Param>::const_iterator it = std::find_if(device_metadata.instanceParameters.begin(),
                                                                 device_metadata.instanceParameters.end(),
                                                                 Util::EqualParam(parameterName));

    return it != device_metadata.instanceParameters.end();
  }


  const std::string &getPrimaryParameter(const std::string & deviceType, int modelLevel) const 
  {
    return getDeviceMetadata(deviceType, modelLevel).primaryParameter;
  }

  int getNumberOfNodes(const std::string & deviceType, int modelLevel) const 
  {
    return getDeviceMetadata(deviceType, modelLevel).numNodes;
  }

  int getNumberOfOptionalNodes(const std::string & deviceType, int modelLevel) const 
  {
    return getDeviceMetadata(deviceType, modelLevel).numOptionalNodes;
  }

  int getNumberOfFillNodes(const std::string & deviceType, int modelLevel) const 
  {
    return getDeviceMetadata(deviceType, modelLevel).numFillNodes;
  }

  bool isModelRequired(const std::string & deviceType, int modelLevel) const 
  {
    return getDeviceMetadata(deviceType, modelLevel).modelRequired;
  }

  const std::vector<Device::Param> &getInstanceParameters(const std::string & deviceType, int modelLevel) const;

  void getInstanceCompositeComponents(const std::string & deviceType, const std::string & parameterName, int modelLevel, std::vector<Device::Param> & components) const;

  void getModelCompositeComponents(const std::string & modelType, const std::string & parameterName, int modelLevel, std::vector<Device::Param> & components) const;

  const std::vector<Device::Param> &getModelParameters(const std::string & modelType, int modelLevel) const;

private:
  // void optionsMetadata();

private:
  mutable DeviceMetadataMap             deviceMetadata_;                ///< Recognizable devices in a netlist.
  mutable DeviceMetadataIndexMap        deviceMetadataIndex;            ///< Translate from model line indexes to deviceMetadata indexes (NPN_1 to Q_1 for BJT)
};

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_CircuitMetadata_h

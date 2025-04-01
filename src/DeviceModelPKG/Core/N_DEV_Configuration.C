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
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : David Baur
//
// Creation Date  : 04/18/2013
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <algorithm>
#include <iomanip>
#include <ostream>
#include <stdexcept>

#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMgr.h>
#include <N_ERH_Message.h>
#include <N_UTL_Demangle.h>

namespace Xyce {
namespace Device {

namespace {

typedef Configuration::ConfigurationMap ConfigurationMap;
typedef std::map<EntityTypeId, Configuration *> EntityTypeIdConfigurationMap;

typedef unordered_map<std::string, EntityTypeId, HashNoCase, EqualNoCase> NameEntityTypeIdMap;
typedef unordered_map<NameLevelKey, EntityTypeId> NameLevelKeyEntityTypeIdMap;

//-----------------------------------------------------------------------------
// Class         : 
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Fri Mar 14 13:38:33 2014
//-----------------------------------------------------------------------------
///
/// Container for the name, model and level forward and backward maps.
///
struct Data
{
  ConfigurationMap              configurationMap_;                      ///< Device registration: Maps (name, level) -> configuration
  EntityTypeIdConfigurationMap  modelTypeConfigurationMap_;             ///< Device registration: Maps model_type_id -> configuration
  NameEntityTypeIdMap           modelTypeNameModelGroupMap_;            ///< Model registration: Maps model_group_name -> model_group_id
  NameLevelKeyEntityTypeIdMap   modelTypeNameLevelModelTypeMap_;        ///< Model registration: Maps (name, level) -> model_type_id


  //-----------------------------------------------------------------------------
  // Function      : ~Data
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Fri Mar 14 13:39:41 2014
  //-----------------------------------------------------------------------------
  ///
  /// Free all the created device configurations.
  ///
  ~Data()
  {
    std::vector<Configuration *> c;
    for (ConfigurationMap::const_iterator it = configurationMap_.begin(); it != configurationMap_.end(); ++it)
      c.push_back((*it).second);
    std::sort(c.begin(), c.end());
    c.erase(std::unique(c.begin(), c.end()), c.end());
    for (std::vector<Configuration *>::iterator it = c.begin(); it != c.end(); ++it)
      delete *it;
  }
};

//-----------------------------------------------------------------------------
// Function      : getData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Feb  5 11:49:42 2014
//-----------------------------------------------------------------------------
/// returns the configuration data singleton.
///
/// @return reference to the configuration data singleton.
///
Data &getData()
{
  static Data data_;

  return data_;
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : Configuration::getConfigurationMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Fri Mar 14 13:23:31 2014
//-----------------------------------------------------------------------------
const Configuration::ConfigurationMap &
Configuration::getConfigurationMap() {
  return getData().configurationMap_;
}

//-----------------------------------------------------------------------------
// Function      : Configuration::findConfiguration
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Fri Mar 14 13:23:52 2014
//-----------------------------------------------------------------------------
const Configuration *
Configuration::findConfiguration(
  ModelTypeId           model_type_id)
{
  EntityTypeIdConfigurationMap::const_iterator it = getData().modelTypeConfigurationMap_.find(model_type_id);
  return it == getData().modelTypeConfigurationMap_.end() ? 0 : (*it).second;
}

//-----------------------------------------------------------------------------
// Function      : Configuration::findConfiguration
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Fri Mar 14 13:23:52 2014
//-----------------------------------------------------------------------------
const Configuration *
Configuration::findConfiguration(
  const std::string &   device_name,
  const int             level)
{
  ConfigurationMap::const_iterator it = getData().configurationMap_.find(ConfigurationMap::key_type(device_name, level));
  return it == getData().configurationMap_.end() ? 0 : (*it).second;
}

//-----------------------------------------------------------------------------
// Function      : Configuration::createDevice
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Fri Mar 14 13:24:07 2014
//-----------------------------------------------------------------------------
Device *
Configuration::createDevice(
  const FactoryBlock &  factory_block) const
{
  return factory(factory_block);
}

//-----------------------------------------------------------------------------
// Function      : Configuration::getModelType
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Fri Mar 14 13:24:19 2014
//-----------------------------------------------------------------------------
EntityTypeId
Configuration::getModelType(const std::string &model_type_name, const int level)
{
  NameLevelKeyEntityTypeIdMap::const_iterator it = getData().modelTypeNameLevelModelTypeMap_.find(NameLevelKey(model_type_name, level));
  if (it == getData().modelTypeNameLevelModelTypeMap_.end())
    return EntityTypeId();
  else
    return (*it).second;
}

//-----------------------------------------------------------------------------
// Function      : Configuration::getModelGroup
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Fri Mar 14 13:37:53 2014
//-----------------------------------------------------------------------------
EntityTypeId
Configuration::getModelGroup(const std::string &device_name)
{
  NameEntityTypeIdMap::const_iterator it = getData().modelTypeNameModelGroupMap_.find(device_name);
  if (it == getData().modelTypeNameModelGroupMap_.end())
    return EntityTypeId();
  else
    return (*it).second;
}

//-----------------------------------------------------------------------------
// Function      : Configuration::addDevice
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Fri Mar 14 13:38:03 2014
//-----------------------------------------------------------------------------
void
Configuration::addDevice(
  const char *          model_name,
  const int             model_level,
  ModelTypeId           model_type_id,
  ModelTypeId           model_group_id,
  int                   model_type_nodes,
  int                   model_group_nodes)
{
  if (model_type_id == model_group_id) {
    std::pair<NameEntityTypeIdMap::iterator, bool> result = getData().modelTypeNameModelGroupMap_.insert(NameEntityTypeIdMap::value_type(model_name, model_type_id));
    if (!result.second && (*result.first).second != model_type_id)
      Report::DevelFatal0().in("Configuration::addDevice")
        << "Attempt to register more than one device model group to the name " << model_name;
  }
  else {
    if (model_type_nodes < model_group_nodes) {
      // Report::DevelWarning0().in("Configuration::addDevice")
      //   << "Registering " << model_name << " level " << model_level << " with " << model_type_nodes
      //   << " nodes which is less than the model group of " << model_group_nodes
      //   << " nodes, effects model name search which skips the model group value before searching for model name";
    }
  }

  {
    std::pair<ConfigurationMap::iterator, bool> result
      = getData().configurationMap_.insert(ConfigurationMap::value_type(NameLevelKey(model_name, model_level), this));
//    if (!result.second)
//      Report::DevelFatal0().in("Configuration::addDevice")
//        << "Device with name " << model_name << " level " << model_level
//        << " already registered as " << (*result.first).second->getName();
  }

  {
    std::pair<EntityTypeIdConfigurationMap::iterator, bool> result
      = getData().modelTypeConfigurationMap_.insert(EntityTypeIdConfigurationMap::value_type(model_type_id, this));
//    if (!result.second && (*result.first).second != this)
 //     Report::DevelFatal().in("Configuration::addDevice")
//        << "Model " << demangle(model_type_id.type().name()) << " already registered to device " << (*result.first).second->getName()
//        << " while trying to register with device " << model_name << " level " << model_level;
  }
}


//-----------------------------------------------------------------------------
// Function      : Configuration::addModel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Fri Mar 14 13:38:16 2014
//-----------------------------------------------------------------------------
void
Configuration::addModel(
  const char *          model_name,
  const int             level,
  ModelTypeId           model_type_id,
  ModelTypeId           model_group_id)
{
  if (model_type_id == model_group_id) {
    std::pair<NameEntityTypeIdMap::iterator, bool> result = getData().modelTypeNameModelGroupMap_.insert(NameEntityTypeIdMap::value_type(model_name, model_type_id));
    if (!result.second && (*result.first).second != model_type_id)
      Report::DevelFatal0().in("Configuration::addDevice")
        << "Attempt to register more than one device model group to the name " << model_name;
  }

  std::pair<NameLevelKeyEntityTypeIdMap::iterator, bool> result
    = getData().modelTypeNameLevelModelTypeMap_.insert(NameLevelKeyEntityTypeIdMap::value_type(NameLevelKey(model_name, level), model_type_id));
  if (!result.second && (*result.first).second != model_type_id)
    Report::DevelFatal0() << "Attempt to register more than one model type to the device " << model_name << " level " << level;

  if (std::find_if(modelTypeNames_.begin(), modelTypeNames_.end(), EqualNoCasePred(model_name)) == modelTypeNames_.end())
    modelTypeNames_.push_back(model_name);
}

} // namespace Device
} // namespace Xyce

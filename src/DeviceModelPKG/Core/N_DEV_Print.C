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
//
//
// Special Notes  :
//
//
// Creator        : David Baur
//
// Creation Date  :
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <map>
#include <utility>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <N_DEV_Algorithm.h>
#include <N_DEV_CompositeParam.h>
#include <N_DEV_Const.h>
#include <N_DEV_Device.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_Dump.h>
#include <N_UTL_IndentStreamBuf.h>
#include <N_DEV_LaTexDoc.h>
#include <N_DEV_Print.h>

namespace Xyce {
namespace Device {

namespace {

struct DeviceInstanceCmp
{
  using result_type = bool;
  using first_argument_type = DeviceInstance;
  using second_argument_type = DeviceInstance;

  bool operator()(const DeviceInstance &entity_0, const DeviceInstance &entity_1) const 
  {
    return less_nocase(entity_0.getName().getEncodedName(), entity_1.getName().getEncodedName());
  }
  bool operator()(const DeviceInstance *entity_0, const DeviceInstance *entity_1) const 
  {
    return less_nocase(entity_0->getName().getEncodedName(), entity_1->getName().getEncodedName());
  }
};

}
//-----------------------------------------------------------------------------
// Function      : printOutModels
// Purpose       : debugging tool
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, factory_block, Parallel Computational Sciences
// Creation Date : 2/03/06
//-----------------------------------------------------------------------------
std::ostream &print(std::ostream &os, const Device &device)
{
  std::vector<DeviceModel *> models;
  getDeviceModels(device, std::back_inserter(models));

  os << std::endl
     << std::endl << section_divider << std::endl
     << "Number of " << device.getName() << " models: " << models.size() << std::endl;

  int i = 0;
  for (std::vector<DeviceModel *>::const_iterator it = models.begin(); it != models.end(); ++it, ++i)
  {
    os << i << ": name = " << (*it)->getName() << " type = " << (*it)->getType() << " level = " << (*it)->getLevel() << std::endl;

    (*it)->printOutInstances(os);
  }
  os << section_divider << std::endl;

  return os;
}

typedef std::map<std::string, std::pair<DeviceModel *, std::vector<DeviceInstance *> > > UglyMap;

struct UglyDeviceModelOp: public DeviceModelOp
{
  UglyDeviceModelOp(const Device &device, UglyMap &model_map)
    : device_(device),
      modelMap_(model_map),
      deviceCount_(0)
    {}

  virtual bool operator()(DeviceModel *model) {
    std::pair<UglyMap::iterator, bool> result = modelMap_.insert(UglyMap::value_type(model->getName(), std::make_pair(model, std::vector<DeviceInstance *>())));

    if (result.second) {
      getDeviceInstances(*model, std::back_inserter((*result.first).second.second));
      if ((*result.first).second.second.empty())
        modelMap_.erase(result.first);
      else {
        deviceCount_ += (*result.first).second.second.size();
        std::sort((*result.first).second.second.begin(), (*result.first).second.second.end(), DeviceInstanceCmp());
      }
    }
    return true;
  }

  const Device &        device_;
  UglyMap &             modelMap_;
  int                   deviceCount_;
};



//-----------------------------------------------------------------------------
// Function      : DeviceMaster:printDotOpOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, factory_block, Parallel Computational Sciences
// Creation Date : 5/4/12
//-----------------------------------------------------------------------------
std::ostream &printDotOpOutput (std::ostream &os, const Device &device)
{
  static const int stdWidth = 15;

// Output models
  UglyMap model_map;

  UglyDeviceModelOp op(device, model_map);
  device.forEachModel(op);

  os << std::endl
     << "Number of " << device.getName() << " models: " << model_map.size() << std::endl
     << std::setw(stdWidth) << "name ";

  if (!model_map.empty())
  {
    int column_count = 1;
    for (UglyMap::const_iterator it_model = model_map.begin(); it_model != model_map.end(); ++it_model, ++column_count)
    {
      if (column_count%8 == 0)
        os << std::endl;
      os << std::setw(stdWidth) << (*it_model).first ;
    }
    os << std::endl;

    os << std::setw(stdWidth) << "type ";
    column_count = 1;
    for (UglyMap::const_iterator it_model = model_map.begin(); it_model != model_map.end(); ++it_model, ++column_count)
    {
      if (column_count%8 == 0)
        os << std::endl;
      os << std::setw(stdWidth) << (*it_model).second.first->getType() ;
    }
    os << std::endl;

    os << std::setw(stdWidth) << "level ";
    column_count = 1;
    for (UglyMap::const_iterator it_model = model_map.begin(); it_model != model_map.end(); ++it_model, ++column_count)
    {
      if (column_count%8 == 0)
        os << std::endl;
      os << std::setw(stdWidth) << (*it_model).second.first->getLevel() ;
    }
    os << std::endl;

    const ParameterMap &parMap = (*model_map.begin()).second.first->getParameterMap();

    column_count = 1;
    for (ParameterMap::const_iterator it_parameter = parMap.begin(); it_parameter != parMap.end(); ++it_parameter)
    {
      for (UglyMap::const_iterator it_model = model_map.begin(); it_model != model_map.end(); ++it_model, ++column_count)
      {
        if (it_model == model_map.begin())
        {
          if (column_count%8 == 0)
            os << std::endl;
          os << std::setw(stdWidth) << (*it_parameter).first;
        }

        os << std::setw(stdWidth);
        printParameter(os, *(*it_model).second.first, (*it_parameter).first, *(*it_parameter).second);
      }
      os << std::endl;
    }
    os << std::endl;

// Output instances
    os << "Number of " << device.getName() << " instances: " << op.deviceCount_ << std::endl
       << std::setw(stdWidth) << "name ";

    column_count = 1;
    for (UglyMap::const_iterator it_model = model_map.begin(); it_model != model_map.end(); ++it_model)
    {
      const std::vector<DeviceInstance *> &instance_list = (*it_model).second.second;

      for (std::vector<DeviceInstance *>::const_iterator it_instance = instance_list.begin(); it_instance != instance_list.end(); ++it_instance, ++column_count)
      {
        if (column_count%8 == 0)
          os << std::endl;
        os << std::setw(stdWidth) << (*it_instance)->getName();
      }
    }
    os << std::endl
       << std::setw(stdWidth) << "model ";


    const DeviceInstance &first_instance = *(*model_map.begin()).second.second.front();

    column_count = 1;
    for (UglyMap::const_iterator it_model = model_map.begin(); it_model != model_map.end(); ++it_model)
    {
      const DeviceModel &model = *(*it_model).second.first;
      const std::vector<DeviceInstance *> &instance_list = (*it_model).second.second;

      if ( !instance_list.empty () ) {
        for (std::vector<DeviceInstance *>::const_iterator it_instance = instance_list.begin(); it_instance != instance_list.end(); ++it_instance, ++column_count)
        {
          if (column_count%8 == 0)
            os << std::endl;
          os << std::setw(stdWidth) << model.getName();
        }
      }
    }

    os << std::endl;

    // output instance parameters:
    const ParameterMap &parameter_map = first_instance.getParameterMap();
    for (ParameterMap::const_iterator it_parameter = parameter_map.begin() ; it_parameter != parameter_map.end(); ++it_parameter)
    {
      os << std::setw(stdWidth) << it_parameter->first;

      column_count = 1;
      for (UglyMap::const_iterator it_model = model_map.begin(); it_model != model_map.end(); ++it_model)
      {
        const DeviceModel &model = *(*it_model).second.first;
        const std::vector<DeviceInstance *> &instance_list = (*it_model).second.second;

        if ( !instance_list.empty () ) {
          for (std::vector<DeviceInstance *>::const_iterator it_instance = instance_list.begin(); it_instance != instance_list.end(); ++it_instance, ++column_count)
          {
            if (column_count%8 == 0)
              os << std::endl;
            os << std::setw(stdWidth);
            printParameter(os, *(*it_instance), (*it_parameter).first, *(*it_parameter).second);
          }
        }
      }
      os << std::endl;
    }
    os << std::endl;
  }

  return os;
}


//-----------------------------------------------------------------------------
// Function      : printParameter
//
// Purpose       : This function finds a parameter in the par table, and then
//                 formats its value in a string.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/4/12
//-----------------------------------------------------------------------------
std::ostream &printParameter(std::ostream &os, const ParameterBase &entity, const std::string &name, const Descriptor &param)
{
  if (param.isType<double>())
  {
    if (isTempParam(name) && param.getAutoConvertTemperature())
    {
      os << param.value<double>(entity) - CONSTCtoK;
    }
    else
    {
      os << param.value<double>(entity);
    }
  }
  else if (param.isType<bool>())
  {
    os << param.value<bool>(entity);
  }
  else if (param.isType<int>())
  {
    os << param.value<int>(entity);
  }
  else if (param.isType<long>())
  {
    os << param.value<long>(entity);
  }
  else if (param.isType<std::string>())
  {
    os << param.value<std::string>(entity);
  }
  else if (param.isType<std::vector<std::string> >())
  {
    os << "string[]";
    // const std::vector<std::string> &string_vector = param.value<std::vector<std::string> >(entity);
    // os << " string[" << string_vector.size() << "]" << std::endl;
    // {
    //   for (std::vector<std::string>::const_iterator it = string_vector.begin(), end = string_vector.end(); it != end; ++it)
    //   {
    //     os << "  " << *it;
    //   }
    // }
  }
  else if (param.isType<std::vector<int> >())
  {
    os << "int[]";
    // const std::vector<int> &int_vector = param.value<std::vector<int> >(entity);
    // os << " int[" << int_vector.size() << "]" << std::endl;
    // {
    //   for (std::vector<int>::const_iterator it = int_vector.begin(), end = int_vector.end(); it != end; ++it)
    //   {
    //     os << "  " << *it;
    //   }
    // }
  }
  else if (param.isType<std::vector<double> >())
  {
    os << "double[]";
    // const std::vector<double> &double_vector = param.value<std::vector<double> >(entity);
    // os << " double[" << double_vector.size() << "]" << std::endl;
    // {
    //   for (std::vector<double>::const_iterator it = double_vector.begin(), end = double_vector.end(); it != end; ++it)
    //   {
    //     os << "  " << *it;
    //   }
    // }
  }
  else if (param.hasCompositeData())
  {

    os << "composite";
    // if (param.isType<CompositeMap>()) {
    //   const CompositeMap &composite_map = param.value<CompositeMap>(entity);
    //   for (CompositeMap::const_iterator it = composite_map.begin(), end = composite_map.end(); it != end; ++it)
    //   {
    //     os << "[" << (*it).first << "] ";
    //     printCompositeParameters(os, *(*it).second);
    //   }
    // }
    // else if (param.isType<CompositeVector>()) {
    //   const CompositeVector &composite_vector = param.value<CompositeVector>(entity);
    //   for (CompositeVector::const_iterator it = composite_vector.begin(), end = composite_vector.end(); it != end; ++it)
    //   {
    //     os << "[" << std::distance(composite_vector.begin(), it) << "] ";
    //     printCompositeParameters(os, *(*it));
    //   }
    // }
   }

  return os;
}

std::ostream &printCompositeParameters(std::ostream &os, const CompositeParam &composite)
{
  const ParameterMap &parameter_map = composite.getParameterMap();
  for (ParameterMap::const_iterator it = composite.getParameterMap().begin(), end = composite.getParameterMap().end(); it != end ;++it)
  {
    printParameter(os, composite, (*it).first, *(*it).second);
    os << " ";
  }

  return os;
}


//-----------------------------------------------------------------------------
// Function      : handleParameterOutputs_
//
// Purpose       : Handles command line arguments like -doc, to output model
//                 parameters either to stdout or to files, depending on which
//                 option is chosen.
//
// Special Notes : erkeite: This function is called from parseCommandLine.  I
//                 moved this code out of that function because that function had
//                 gotten too long and confusing.
//
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 10/31/07
//-----------------------------------------------------------------------------
void
handleParameterOutputs(
  OutputMode::Mode      format,
  std::string           option_device_name,
  int                   option_device_level,
  bool                  print_model,
  bool                  print_instance)
{
  typedef std::map<std::pair<std::string, int>, Device *> DeviceMap;

  if (format == OutputMode::PARAM) {
    typedef std::map<NameLevelKey, Configuration *, NameLevelKeyLess> OrderedConfigurationMap;
    const OrderedConfigurationMap configuration_map(Configuration::getConfigurationMap().begin(), Configuration::getConfigurationMap().end());
    for (OrderedConfigurationMap::const_iterator it = configuration_map.begin(), end = configuration_map.end(); it != end; ++it)
    {
      Xyce::dout() << (*it).first << Xyce::Util::push << std::endl
                   << *(*it).second << Xyce::Util::pop << std::endl;
    }
  }
  else {
    const Configuration::ConfigurationMap &configuration_map = Configuration::getConfigurationMap();
    for (Configuration::ConfigurationMap::const_iterator it = configuration_map.begin(), end = configuration_map.end(); it != end; ++it)
    {
      const NameLevelKey &device_key = (*it).first;

      std::string device_name = (*it).first.first;
      const int device_level = (*it).first.second;

      device_name[0] = toupper(device_name[0]);

      if ((option_device_name.empty() || Xyce::equal_nocase(option_device_name, device_name)) && (option_device_level == -1 || option_device_level == device_level)) {
        Configuration &configuration = *(*it).second;

        std::string device_description = configuration.getName();
        ParametricData<void> &instance_parameters = configuration.getInstanceParameters();
        ParametricData<void> &model_parameters = configuration.getModelParameters();

        if (configuration.getName() == "Behavioral Digital")
          device_name = "Digital";

        if (print_instance && !instance_parameters.getMap().empty()) {
          std::ostringstream path;
          path << device_name << "_" << device_level << "_Device_Instance" << "_Params.tex";

          std::ofstream os(path.str().c_str(), std::ios_base::out);
          os << "% This table was generated by Xyce:" << std::endl
             << "%   Xyce " << (format == OutputMode::DOC_CAT ? "-doc_cat " : "-doc ")
             << device_name << " " << device_level << std::endl
             << "%" << std::endl;

          laTexDevice(os, device_name, device_level, 0, device_description, instance_parameters, format);
        }

        if (print_model && !model_parameters.getMap().empty()) {
          std::ostringstream path;
          path << device_name << "_" << device_level << "_Device_Model" << "_Params.tex";

          std::ofstream os(path.str().c_str(), std::ios_base::out);

          os << "% This table was generated by Xyce:" << std::endl
             << "%   Xyce " << (format == OutputMode::DOC_CAT ? "-doc_cat " : "-doc ")
             << device_name << " " << device_level << std::endl
             << "%" << std::endl;

            laTexDevice(os, device_name, device_level, 1, device_description, model_parameters, format);
        }
      }
    }
  }
}

} // namespace Device
} // namespace Xyce


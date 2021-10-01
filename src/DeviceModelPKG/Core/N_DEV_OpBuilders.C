//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Purpose        : Output Manager
//
// Special Notes  :
//
// Creator        : David G. Baur, Raytheon
//
// Creation Date  : 11/11/2014
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <string>
#include <vector>

#include <N_DEV_fwd.h>
#include <N_UTL_fwd.h>

#include <N_DEV_Op.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_DeviceEntity.h>
#include <N_UTL_Algorithm.h>
#include <N_UTL_OpBuilder.h>
#include <N_UTL_Param.h>
#include <N_UTL_HspiceBools.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Namespace     : Unnamed
// Purpose       : file-local scoped methods and data
// Special Notes :
// Creator       : Dave Baur, Raytheon
// Creation Date : 11/11/2014
//-----------------------------------------------------------------------------
namespace {

//-----------------------------------------------------------------------------
// Function      : parameterNameAndArgs
// Purpose       : given a parameter list iterator, construct the string
//                 name of the function with arguments
// Special Notes : This function is never used here.  Similar functions exist
//                 in other OpBuilder files such as the IOInterfacePKG one
// Scope         : file-local
// Creator       : Dave Baur, Raytheon
// Creation Date : 11/11/2014
//-----------------------------------------------------------------------------
///
/// Construct a parameter name with arguments from a ParamList iterator
///
/// @param[out] name the string name of the function
/// @param[out] args a vector of strings for the arguments
/// @param[in] it iterator pointing into a ParamList
///
/// @note This function is currently unused in the device package.  It
/// is able to reconstruct accessors such as V(A,B), I(INSTANCE), and N(var)
/// from a parameter list.
///
/// @author Dave Baur, Raytheon
/// @date 11/11/2014
///
void parameterNameAndArgs(std::string &name, std::vector<std::string> &args, Util::ParamList::const_iterator &it)
{
  const std::string &param_tag = (*it).tag();

  if (param_tag[0] == 'V' || param_tag[0] == 'I' || param_tag[0] == 'N')
  {
    std::ostringstream oss;
    oss << param_tag << "(";
    int arg_count = (*it).getImmutableValue<int>();
    for (int i = 0; i < arg_count; ++i)
    {
      ++it;
      if (i != 0)
        oss << ",";
      oss << (*it).tag();
      args.push_back((*it).tag());
    }
    oss << ")";
    name = oss.str();
  }
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Class         : DeviceGlobalParameterOpBuilder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon
// Creation Date : 11/11/14
//-----------------------------------------------------------------------------
///
/// Provide operator access to global parameters
///
/// @note Global parameters are currently owned by the device package,
/// which is probably inappropriate
///
/// @author Dave Baur, Raytheon
/// @date 11/11/2014
///
struct DeviceGlobalParameterOpBuilder : public Util::Op::Builder
{
  /// Constructor
  DeviceGlobalParameterOpBuilder(const DeviceMgr &device_manager)
    : deviceManager_(device_manager)
  {}

  /// Destructor
  virtual ~DeviceGlobalParameterOpBuilder()
  {}

//-----------------------------------------------------------------------------
// Function      : DeviceGlobalParameterOpBuilder::registerCreateFunctions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : 11/11/14
//----------------------------------------------------------------------------
/// Register ops that this builder can create
  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<DeviceMgrGlobalParameterOp>();
  }

//-----------------------------------------------------------------------------
// Function      : DeviceGlobalParameterOpBuilder::makeOp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : 11/11/14
//----------------------------------------------------------------------------
///
/// Attempt to create an op for the global parameter provided
///
/// @param it ParamList iterator pointing to the parameter
/// @result a new operator if succeeded, null if failed
///
/// @author Dave Baur, Raytheon
/// @date 11/11/14
///
  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();
    const std::string &param_string = (*it).stringValue();

    if (param_tag == "GLOBAL_PARAMETER")
    {
      // new_op  = new DeviceMgrGlobalParameterOp(param_string, deviceManager_, param_string);
      // new_op->addArg(param_string);
      const double *result = deviceManager_.findGlobalPar(param_string);
      if (result)
      {
        // Refactor: [DGB] Should use the address.
        new_op = new DeviceMgrGlobalParameterOp(param_string, deviceManager_, *result);
        new_op->addArg(param_string);
      }
    }
    else
    {
      // Don't even look for this parameter name if we are actually
      // a function call with arguments!  (all of these things are
      // the first characters of valid access operators (voltage, current
      // internal vars, power, noise, S params, Y params, Z params)
      if ( !((*it).getType() == Util::INT
             && (param_tag[0] == 'V'
                                  || param_tag[0] == 'I' || param_tag[0] == 'N'
                                  || param_tag[0] == 'P' || param_tag[0] == 'W'
                                  || param_tag[0] == 'D' || param_tag[0] == 'S'
                 || param_tag[0] == 'Y' || param_tag[0] == 'Z')
             && (*it).getImmutableValue<int>() >0 )
           )
      {
        const double *result = deviceManager_.findGlobalPar(param_tag);
        if (result)
        {
          // Refactor: [DGB] Should use the address.
          new_op = new DeviceMgrGlobalParameterOp(param_tag, deviceManager_, *result);
        }
      }
    }

    return new_op;
  }

private:
  const DeviceMgr &     deviceManager_;
};

//-----------------------------------------------------------------------------
// Class         : DeviceEntityOpBuilder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon
// Creation Date : 11/11/14
//-----------------------------------------------------------------------------
///
/// Provide operator access to parameters in device entities
///
/// @author Dave Baur, Raytheon
/// @date 11/11/2014
///
struct DeviceEntityOpBuilder : public Util::Op::Builder
{
  /// Constructor
  DeviceEntityOpBuilder(const DeviceMgr &device_manager)
    : deviceManager_(device_manager)
  {}

  /// Destructor
  virtual ~DeviceEntityOpBuilder()
  {}

//-----------------------------------------------------------------------------
// Function      : DeviceEntityOpBuilder::registerCreateFunctions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : 11/11/14
//----------------------------------------------------------------------------
/// Register ops that this builder can create
  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<DeviceEntityParameterOp>();
  }

//-----------------------------------------------------------------------------
// Function      : DeviceEntityOpBuilder::makeOp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : 11/11/14
//----------------------------------------------------------------------------
///
/// Attempt to create an op for the device entity parameter provided
///
/// @param it ParamList iterator pointing to the parameter
/// @result a new operator if succeeded, null if failed
///
/// @author Dave Baur, Raytheon
/// @date 11/11/14
///
  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    int arg_count = 0;
    if ((*it).getType() == Util::INT)
      arg_count = (*it).getImmutableValue<int>();
    if (arg_count == 0)
    {
      const DeviceEntity *device_entity = deviceManager_.getDeviceEntity(param_tag);
      if (device_entity)
      {
        std::string param_name = Util::paramNameFromFullParamName(param_tag);
        if (device_entity->findParam(param_name))
        {
          // the typical case of a fully-specified <deviceName:paramName> pair like R1:R
          new_op = new DeviceEntityParameterOp(param_tag, *device_entity, param_name);
        }
        else
        {
          // The less-common case of a device, such as the R device, that has
          // a default instance parameter.  If the device has a default parameter
          // then try to use it to make the op.
          std::string default_param_name = device_entity->getDefaultParamName();
          if ( !(default_param_name.empty()) && deviceManager_.getDeviceEntity(param_tag + Xyce::Util::separator + default_param_name))
          {
            new_op = new DeviceEntityParameterOp(param_tag, *device_entity, default_param_name);
          }
        }
      }
    }

    return new_op;
  }

private:
  const DeviceMgr &     deviceManager_;
};

//-----------------------------------------------------------------------------
// Class         : MutualInductorInstancesOpBuilder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/30/2021
//-----------------------------------------------------------------------------
///
/// Special Op for inductors that have been merged into mutual inductor devices.
///
/// @author Eric Keiter, SNL
/// @date 09/30/2021
///
struct MutualInductorInstancesOpBuilder : public Util::Op::Builder
{
  /// Constructor
  MutualInductorInstancesOpBuilder(const DeviceMgr &device_manager)
    : deviceManager_(device_manager)
  {}

  /// Destructor
  virtual ~MutualInductorInstancesOpBuilder()
  {}

//-----------------------------------------------------------------------------
// Function      : MutualInductorInstancesOpBuilder::registerCreateFunctions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/30/2021
//----------------------------------------------------------------------------
/// Register ops that this builder can create
  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<DeviceEntityParameterOp>();
  }

//-----------------------------------------------------------------------------
// Function      : MutualInductorInstancesOpBuilder::makeOp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/30/2021
//----------------------------------------------------------------------------
///
/// @author Eric Keiter, SNL
/// @date 09/30/2021
///
  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const 
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    const std::string full_inductor_name = Xyce::Util::entityNameFromFullParamName(param_tag).getEncodedName();
    Xyce::Device::InstanceName foo(full_inductor_name);
    std::string inductor_name = foo.getDeviceName();

    int arg_count = 0;
    if ((*it).getType() == Util::INT) // double check; does this part make sense here?
      arg_count = (*it).getImmutableValue<int>();
    if (arg_count == 0)
    {
      int inductorIndex=-1;
      DeviceInstance *device_instance=deviceManager_.getMutualInductorDeviceInstance(inductor_name,inductorIndex);

      if (device_instance)
      {
        new_op = new MutualInductorInstancesOp(param_tag, inductor_name, *device_instance, inductorIndex);
      }
    }

    return new_op;
  }

private:
  const DeviceMgr &     deviceManager_;
};

//-----------------------------------------------------------------------------
// Class         : ArtificialParameterOpBuilder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon
// Creation Date : 11/11/14
//-----------------------------------------------------------------------------
///
/// Provide operator access to artificial parameters (continuation params)
///
/// @author Dave Baur, Raytheon
/// @date 11/11/2014
///
struct ArtificialParameterOpBuilder : public Util::Op::Builder
{
  /// Constructor
  ArtificialParameterOpBuilder(const DeviceMgr &device_manager, const ArtificialParameterMap &artificial_parameter_map, const PassthroughParameterSet &passthrough_parameter_set)
    : deviceManager_(device_manager),
      artificialParameterMap_(artificial_parameter_map),
      passthroughParameterSet_(passthrough_parameter_set)
  {}

  /// Destructor
  virtual ~ArtificialParameterOpBuilder()
  {}

//-----------------------------------------------------------------------------
// Function      : ArtificialParameterOpBuilder::registerCreateFunctions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : 11/11/14
//----------------------------------------------------------------------------
/// Register ops that this builder can create
  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<ArtificialParameterOp>();
  }

//-----------------------------------------------------------------------------
// Function      : ArtificialParameterOpBuilder::makeOp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : 11/11/14
//----------------------------------------------------------------------------
///
/// Attempt to create an op for the artificial continuation parameter provided
///
/// @param it ParamList iterator pointing to the parameter
/// @result a new operator if succeeded, null if failed
///
/// @author Dave Baur, Raytheon
/// @date 11/11/14
///
  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();
    const std::string &param_string = (*it).stringValue();

    const ArtificialParameterMap::const_iterator artificial_parameter_it = artificialParameterMap_.find(param_tag);
    if (artificial_parameter_it != artificialParameterMap_.end())
    {
      std::string param_name = Util::paramNameFromFullParamName(param_tag);
      new_op = new ArtificialParameterOp(param_tag, deviceManager_, *(*artificial_parameter_it).second, param_name);
    }

    return new_op;
  }

private:
  const DeviceMgr &                     deviceManager_;
  const ArtificialParameterMap &        artificialParameterMap_;
  const PassthroughParameterSet &       passthroughParameterSet_;
};


//-----------------------------------------------------------------------------
// Class         : DeviceOptionsOpBuilder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon
// Creation Date : 11/11/14
//-----------------------------------------------------------------------------
///
/// Provide operator access to device package options
///
/// @author Dave Baur, Raytheon
/// @date 11/11/2014
///
struct DeviceOptionsOpBuilder : public Util::Op::Builder
{
  /// Constructor
  DeviceOptionsOpBuilder(const DeviceOptions &device_options)
    : deviceOptions_(device_options)
  {}

  /// Destructor
  virtual ~DeviceOptionsOpBuilder()
  {}

//-----------------------------------------------------------------------------
// Function      : DeviceOptionsOpBuilder::registerCreateFunctions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : 11/11/14
//----------------------------------------------------------------------------
/// Register ops that this builder can create
  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<DeviceOptionsOp>();
  }

//-----------------------------------------------------------------------------
// Function      : DeviceOptionsOpBuilder::makeOp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : 11/11/14
//----------------------------------------------------------------------------
///
/// Attempt to create an op for the artificial continuation parameter provided
///
/// @param it ParamList iterator pointing to the parameter
/// @result a new operator if succeeded, null if failed
///
/// @note Currently, the only recognized parameter is "gmin"
///
/// @author Dave Baur, Raytheon
/// @date 11/11/14
///
  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();
    const std::string &param_string = (*it).stringValue();

    if (compare_nocase(param_tag.c_str(), "gmin") == 0)
    {
      new_op = new DeviceOptionsOp(param_tag, deviceOptions_, param_tag);
    }
    return new_op;
  }

private:
  const DeviceOptions &     deviceOptions_;
};

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::registerOpBuilders
// Purpose       : Register all Device package op builders with builder manager
// Special Notes :
// Scope         : public
// Creator       : Dave Baur, Raytheon
// Creation Date : 11/11/2014
//-----------------------------------------------------------------------------
///
/// Register all Device package op builders with builder manager
///
/// param[in] builder_manager  Reference to the builder manager
/// param[in] comm  parallel communicator
/// param[in] device_manager reference to the device manager
///
/// This function registers the DeviceGlobalParameter, DeviceEntity,
/// DeviceOptions, and ArtificialParameter op builders with the builder
/// manager.
///
/// @author Dave Baur, Raytheon
/// @date 11/11/2014
///
void registerOpBuilders(Util::Op::BuilderManager &builder_manager, Parallel::Machine comm, DeviceMgr &device_manager)
{
  builder_manager.addBuilder(new DeviceGlobalParameterOpBuilder(device_manager));
  builder_manager.addBuilder(new DeviceEntityOpBuilder(device_manager));
  builder_manager.addBuilder(new MutualInductorInstancesOpBuilder(device_manager));
  builder_manager.addBuilder(new DeviceOptionsOpBuilder(device_manager.getDeviceOptions()));
  builder_manager.addBuilder(new ArtificialParameterOpBuilder(device_manager, device_manager.getArtificialParameterMap(), device_manager.getPassthroughParameterSet()));
//  builder_manager.addBuilder(new DeviceMgrOpBuilder(device_manager));
}

} // namespace Device
} // namespace Xyce

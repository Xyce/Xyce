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
// Purpose        : Holds device configuration information
//
// Special Notes  :
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355
//
// Creation Date  : 2013/04/18 18:01:27
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Configuration_h
#define Xyce_N_DEV_Configuration_h

#include <unordered_map>
using std::unordered_map;
#include <string>
#include <vector>

#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_DEV_Pars.h>

namespace Xyce {
namespace Device {

struct ModelGroupLeadTraits
{};

template<class M, class G>
struct ModelGroupType_
{
  typedef typename G::ModelType GroupType_;
};

template<class M>
struct ModelGroupType_<M, ModelGroupLeadTraits>
{
  typedef M GroupType_;
};

//-----------------------------------------------------------------------------
// Class         : DeviceTraits
// Purpose       : Defines the device configuration data and grouping
// Special Notes :
// Scope         : public
// Creator       : Dave Baur, Raytheon
// Creation Date : 1/27/2014
//-----------------------------------------------------------------------------
///
/// The DeviceTraits template describes the configuration of a device.
///
/// Generally a device uses a typedef to define its traits:
///
///   typedef DeviceTraits<Instance, Model> Traits;                      ///< Device traits for a group device
///
/// or
///
///   typedef DeviceTraits<Instance, Model, Resistor::Traits> Traits     ///< Device is part of the Resistor device group
///
/// Every device has an Instance class, a Model class and a ModelGroup class.  An Instance is created for every device in
/// the netlist.  A Model class is created for each unique model definition in the net list.  The ModelGroup serves as a
/// device type grouping mechanism and is generally the level 1 device of a model with the other models having higher
/// level values.
///
/// The EntityTypeId, ModelTypeId and InstanceTypeId are all descriptive typedefs for Util::type_index.  These model
/// type identifiers and instance type identifiers are uniquely defined for each instance and model.  The model type
/// identifier also serves as the model group identifier.  The model groups are identified when the model traits is
/// ModeGroupLeadTraits which is the case for the two parameter DeviceTraits trait.
///
/// Every device also has or may have nodes, optional nodes, fill nodes, a model required flags, a primary parameter, a
/// linear device flags and PDE device flag.  These values are all provided via static traits functions, most of which
/// have default implementations.
///
/// Every device must specify not only the device relationships above as template parameters, but must also define the
/// numNodes() and isLinearDevice() functions as providing a default value was too error prone.
///
///
/// DEVICE REGISTRATION:
///   Each device calls registerDevice() to register the device with the configuration.  Each device must also call
///   registerModelType() for each model type to be defined.
///
///   registerDevice()
///     - adds to the configurationMap_: (device_name, device_level) -> configuration
///     - adds to the modelTypeConfigurationMap_: model_type_id -> configuration
///
///   registerModel()
///     - if model_type_id == model_group_id, adds to the modelTypeNameModelGroupMap_: model_name -> model_group_id
///     - adds to the modelTypeNameLevelModelTypeMap_: (model_name, level) -> model_type_id
///     - adds the model_name to the modelTypeNames list of device configuration
///
/// THE CONFIGURATION IS CURRENTLY STORED IN A SINGLETON.  This may not be as issue in that the configuration is
/// constructed exactly once via device registration and plugin loading at Xyce startup.
///
/// The default parameter for a device is the parameter that will be used
/// if a single untagged parameter value is given on an instance line.  For
/// example, the resistor's default parameter is R, and on the
/// line
/// @code
/// R1 A B 5k
/// @endcode
/// the untagged value (5k) will be assigned to the resistor's
/// R parameter.
///
/// A primary parameter is the parameter that will be changed by a
/// .STEP loop or natural parameter homotopy when only the device name
/// is given.  In the resistor, the primary parameter is also "R", and
/// so
/// @code
/// .STEP R1 1 10 1
/// @endcode
/// will step R1's R parameter from 1 to 10.  In
/// most existing devices the primary parameter is the same as the
/// default parameter, but it is not required that this be the
/// case.
///
/// Optional nodes are nodes that may be specified on the instance line, and
/// if not specified are internal nodes instead of external.  Fill nodes are
/// like optional nodes, but if not specified they are connected to ground.
/// Examples of optional nodes include the internal and external body contacts
/// in the BSIM3 SOI model.  Examples of fill nodes include the substrate
/// node of the level 1 BJT model.
template <class M, class I, class G>
class DeviceTraits
{
public:
  typedef I InstanceType;                                     ///< Make instance template parameter available
  typedef M ModelType;                                        ///< Make model template parameter available
  typedef G ModelGroupTraits;                                 ///< Make model group traits template parameter available
  typedef typename ModelGroupType_<M, G>::GroupType_ ModelGroupType;         ///< Make model group template parameter available

  //-----------------------------------------------------------------------------
  // Function      : instanceType
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Fri Mar 14 13:00:58 2014
  //-----------------------------------------------------------------------------
  ///
  /// Returns the instance type identifier
  ///
  /// @return the typeid wrapper for the instance type identifier
  ///
  static const InstanceTypeId instanceType()
  {
    return EntityTypeId(typeid(InstanceType));
  }

  //-----------------------------------------------------------------------------
  // Function      : modelType
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Feb  5 10:41:46 2014
  //-----------------------------------------------------------------------------
  ///
  /// Returns the model type identifier
  ///
  /// @return the typeid wrapper for the model type identifier
  ///
  static const ModelTypeId modelType()
  {
    return EntityTypeId(typeid(ModelType));
  }

  //-----------------------------------------------------------------------------
  // Function      : modelGroupType
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Feb  5 10:42:39 2014
  //-----------------------------------------------------------------------------
  ///
  /// Returns the model group identifier
  ///
  /// @return the typeid wrapper for the model group identifier
  ///
  static const ModelTypeId modelGroup()
  {
    return EntityTypeId(typeid(ModelGroupType));
  }

  //-----------------------------------------------------------------------------
  // Function      : factory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Feb  6 10:20:07 2014
  //-----------------------------------------------------------------------------
  ///
  /// Creates a device via the device template
  ///
  /// Generally implemented as:
  ///
  /// @code
  /// static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block)
  /// {
  ///  return new DeviceMaster<Traits>(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
  /// }
  /// @endcode
  ///
  /// or
  ///
  /// @code
  /// static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block)
  /// {
  ///  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
  /// }
  /// @endcode
  ///
  /// when Master extends DeviceMaster.
  ///
  /// @param configuration     const reference to the device configuration
  /// @param factory_block     const reference to the factory parameters
  ///
  /// @return pointer to the newly created device
  ///
  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);

  static int numNodes();                                ///< Number of nodes must be provided in the deriving class
  static int numOptionalNodes() {return 0;}             ///< Default number of optional nodes
  static int numFillNodes() {return 0;}                 ///< Default number of fill nodes
  static bool modelRequired() {return false;}           ///< By default, model is not required
  static const char *primaryParameter() {return "";}          ///< By default, there is no primary parameter
  static const char *instanceDefaultParameter() {return "";}  ///< By default, there is no instance default parameter
  static bool isLinearDevice();                         ///< Linear device flag must be provided in the deriving class
  static bool isPDEDevice() {return false;}             ///< By default, device is not a PDE device
};

//-----------------------------------------------------------------------------
// Class         : FactoryBlock
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Jan 27 14:49:12 2014
//-----------------------------------------------------------------------------
///
/// The FactoryBlock contains parameters needed by the device,
/// instance and model creation functions.  This allows additional
/// parameter to be added without the need to change the interface.
///
/// The DeviceMgr class generally calls the factory functions and owns
/// these objects, however this is by no means a requirement.
///
struct FactoryBlock
{
    //-----------------------------------------------------------------------------
    // Function      : FactoryBlock
    // Purpose       :
    // Special Notes :
    // Scope         : public
    // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
    // Creation Date : Fri Mar 14 13:05:27 2014
    //-----------------------------------------------------------------------------
    ///
    /// The FactoryBlock constructs serves to pass data to the device factory functions
    ///
    /// @invariant These references must exist through the execution of Xyce
    ///
    /// @param device_options 
    /// @param solver_state 
    /// @param matrix_load_data 
    /// @param extern_data 
    /// @param command_line 
    ///
    FactoryBlock(
      const DeviceMgr &         device_manager,
      const DeviceOptions &     device_options,
      const SolverState &       solver_state,
      MatrixLoadData &          matrix_load_data,
      const ExternData &        extern_data,
      const IO::CmdParse &      command_line)
      : deviceManager_(device_manager),
        deviceOptions_(device_options),
        solverState_(solver_state),
        externData_(extern_data),
        matrixLoadData_(matrix_load_data),
        commandLine_(command_line)
    {}

    const DeviceMgr &           deviceManager_;
    const DeviceOptions &       deviceOptions_;
    const SolverState &         solverState_;
    const ExternData &          externData_;
    MatrixLoadData &            matrixLoadData_;
    const IO::CmdParse &        commandLine_;
};

/// Class Configuration contains device configuration data.
///
/// The DeviceTraits template describes the device and configuration is an object that holds the data obtained via the
/// traits.
///
/// The configuration also maintains several containers that allow device, instance and model data and creation functions
/// to be located by name and level.
///
class Configuration
{
public:
  typedef unordered_map<NameLevelKey, Configuration *> ConfigurationMap;

  //-----------------------------------------------------------------------------
  // Function      : getConfigurationMap
  // Purpose       : Return the configuration map
  // Special Notes :
  // Scope         : public
  // Creator       : Dave Baur, Raytheon
  // Creation Date : 2/5/2014
  //-----------------------------------------------------------------------------
  ///
  /// Returns the configuration map of the system
  ///
  /// The configuration map maps the device name and level pair to its configuration.
  ///
  /// @return const reference to the device configuration map.
  static const ConfigurationMap &getConfigurationMap();

  //-----------------------------------------------------------------------------
  // Function      : findConfiguration
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Feb  5 10:54:24 2014
  //-----------------------------------------------------------------------------
  ///
  /// Returns the configuration associated with the device name and level or 0 if not found
  ///
  /// @param device_name       const reference to the device name
  /// @param level             device level
  ///
  /// @return const pointer to the specified device's configuration
  ///
  static const Configuration *findConfiguration(ModelTypeId model_type_id);

  //-----------------------------------------------------------------------------
  // Function      : findConfiguration
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Feb  5 10:54:24 2014
  //-----------------------------------------------------------------------------
  ///
  /// Returns the configuration associated with the device name and level or 0 if not found
  ///
  /// @param device_name       const reference to the device name
  /// @param level             device level
  ///
  /// @return const pointer to the specified device's configuration
  ///
  static const Configuration *findConfiguration(const std::string &device_name, const int level);

  //-----------------------------------------------------------------------------
  // Function      : getModelType
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Feb  5 10:54:24 2014
  //-----------------------------------------------------------------------------
  ///
  /// Returns the model type identifier for the specified device, if
  /// the model type and level is not found, return an invalid type.
  ///
  /// @param model_name       const reference to the model name
  /// @param level             device level
  ///
  /// @return the model type identifier for the specified device model
  ///
  static ModelTypeId getModelType(const std::string &model_name, const int level);

  //-----------------------------------------------------------------------------
  // Function      : getModelGroup
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Feb  5 10:54:24 2014
  //-----------------------------------------------------------------------------
  ///
  /// Returns the model group of the device name
  ///
  /// A device is registered as the model group when its model type identifier equals its model group identifier.
  ///
  /// @param device_name       const reference to the device name
  ///
  /// @return model group identifier for the specified device
  ///
  static ModelTypeId getModelGroup(const std::string &device_name);

protected:
  //-----------------------------------------------------------------------------
  // Function      : Configuration
  // Purpose       :
  // Special Notes :
  // Scope         : protected
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon Jan 27 15:06:50 2014
  //-----------------------------------------------------------------------------
  ///
  /// Populates the device configuration object.
  ///
  /// Note that the instanceParameters_ and modelParameters_ each reference the instanceParameters__ and
  /// modelParameters__ object of the derived Config<Traits> class.  This really just keeps that casting down and could
  /// be implemented as ParametricData<void> objects here.
  ///
  /// @param instance_parameters
  /// @param model_parameters
  /// @param name
  /// @param device_type_name
  /// @param num_nodes
  /// @param num_optional_nodes
  /// @param num_fill_nodes
  /// @param model_required
  /// @param primary_parameter
  ///
  Configuration(
    ParametricData<void> &    instance_parameters,
    ParametricData<void> &    model_parameters,
    const char *              name,
    const char *              device_type_name,
    const char *              instance_default_parameter_name,
    int                       num_nodes,
    int                       num_optional_nodes,
    int                       num_fill_nodes,
    bool                      model_required,
    bool                      linear_device,
    bool                      pde_device,
    const char *              primary_parameter)
  : instanceParameters_(instance_parameters),
    modelParameters_(model_parameters),
    name_(name),
    deviceTypeName_(device_type_name),
    instanceDefaultParameterName_(instance_default_parameter_name),
    numNodes_(num_nodes),
    numOptionalNodes_(num_optional_nodes),
    numFillNodes_(num_fill_nodes),
    modelRequired_(model_required),
    linearDevice_(linear_device),
    pdeDevice_(pde_device),
    primaryParameter_(primary_parameter),
    modelTypeNames_()
  {}

public:
  virtual ~Configuration()
  {}

private:
  Configuration(const Configuration &);               ///< No copies allowed
  Configuration &operator=(const Configuration &);    ///< No assignment allowed

public:
  ///  Returns the instance parameter descriptions
  ///
  ///  @return reference to the instance parameter descriptions
  ///
  ParametricData<void> &getInstanceParameters() const
  {
    return instanceParameters_;
  }

  ///  Returns the model parameter descriptions
  ///
  ///  @return reference to the model parameter descriptions
  ///
  ParametricData<void> &getModelParameters() const
  {
    return modelParameters_;
  }

  ///  Returns the device name
  ///
  ///  @return const reference to the device name
  ///
  const std::string &getName() const
  {
    return name_;
  }

  ///  Returns the device type name
  ///
  ///  @return const reference to the device type name
  ///
  const std::string &getDeviceTypeName() const
  {
    return deviceTypeName_;
  }

  ///  Returns the instance default parameter name
  ///
  ///  @return const reference to the instance default parameter name
  ///
  const std::string &getInstanceDefaultParameterName() const
  {
    return instanceDefaultParameterName_;
  }

  ///
  ///  Returns the number of nodes of this device
  ///
  ///  @return number of nodes of this device
  ///
  int getNumNodes() const
  {
    return numNodes_;
  }

  ///
  ///  Returns the number of optional nodes of this device
  ///
  ///  @return number of optional nodes of this device
  ///
  int getNumOptionalNodes() const
  {
    return numOptionalNodes_;
  }

  ///
  ///  Returns the number of fill nodes of this device
  ///
  ///  @return number of fill nodes of this device
  ///
  int getNumFillNodes() const
  {
    return numFillNodes_;
  }

  ///  Returns true of the model must be specified when defining an instance of this device
  ///
  ///  @return true if the model must be specified
  ///
  bool getModelRequired() const
  {
    return modelRequired_;
  }

  ///  Returns true of the model is a linear device
  ///
  ///  @return true if the model is a linear deive
  ///
  bool getLinearDevice() const
  {
    return linearDevice_;
  }

  ///  Returns true of the model is a PDE device
  ///
  ///  @return true if the model is a PDE device
  ///
  bool getPDEDevice() const
  {
    return pdeDevice_;
  }

  ///  Returns the default primary parameter of this device or the empty string if there is no primary parameter
  ///
  ///  @return primary parameter for this device.
  ///
  const std::string &getPrimaryParameter() const
  {
    return primaryParameter_;
  }

  ///  Returns a vector of strings that name all the model types defined of this model
  ///
  ///  @return const reference to the vector of string names of defined model types
  ///
  const std::vector<std::string> &getModelTypeNames() const
  {
    return modelTypeNames_;
  }

  //-----------------------------------------------------------------------------
  // Function      : findConfiguration
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Feb  5 10:54:24 2014
  //-----------------------------------------------------------------------------
  ///
  /// Creates the specified device
  ///
  /// @param model_type        model type of the device to create
  /// @param factory_block     parameters provided to the factory function
  ///
  /// @return
  ///
  Device *createDevice(const FactoryBlock &factory_block) const;

  //-----------------------------------------------------------------------------
  // Function      : factory
  // Purpose       :
  // Special Notes :
  // Scope         : protected
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon Jan 27 15:06:50 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Overriding function creates an instance this device
  ///
  ///  Ownership of this object is passed to the caller which is responsible to destory the created object.
  ///
  ///  @param factory_block reference to the parameter block
  ///
  ///  @return pointer to the newly created device.
  ///
  virtual Device *factory(const FactoryBlock &factory_block) const = 0;

  //-----------------------------------------------------------------------------
  // Function      : instancetype
  // Purpose       :
  // Special Notes :
  // Scope         : protected
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon Jan 27 15:06:50 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Returns the instance type identifier
  ///
  ///  @return the instance type identifier
  ///
  virtual InstanceTypeId instanceType() const = 0;

  //-----------------------------------------------------------------------------
  // Function      : modelType
  // Purpose       :
  // Special Notes :
  // Scope         : protected
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon Jan 27 15:06:50 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Returns the model type identifier
  ///
  ///  @return the model type identifier
  ///
  ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  ///  @date   Wed Feb  5 10:08:52 2014
  ///
  virtual ModelTypeId modelType() const = 0;

  //-----------------------------------------------------------------------------
  // Function      : modelGroup
  // Purpose       :
  // Special Notes :
  // Scope         : protected
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon Jan 27 15:06:50 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Returns the model group identifier
  ///
  ///  @return the model group identifier
  ///
  virtual ModelTypeId modelGroup() const = 0;

protected:
  //-----------------------------------------------------------------------------
  // Function      : addDevice
  // Purpose       :
  // Special Notes :
  // Scope         : protected
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon Jan 27 15:06:50 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Adds the device to the configuration
  ///
  ///  @invariant only one device can designate itself as the group lead
  ///  @invariant device and level must be unique
  ///
  ///  @param device_name       const character pointer to the device name
  ///  @param device_level      device level
  ///  @param model_type_id     model type identifier of the device
  ///
  void addDevice(const char *model_name, const int model_level, ModelTypeId model_type_id, ModelTypeId model_group_id, int model_type_nodes, int model_group_nodes);

  //-----------------------------------------------------------------------------
  // Function      : addModel
  // Purpose       :
  // Special Notes :
  // Scope         : protected
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon Jan 27 15:06:50 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Adds a model name for the device to the configuration
  ///
  ///  @invariant only one device can designate itself as the group lead
  ///  @invariant model names must be unique within a device configuration
  ///
  ///  @param model_name        const character pointer to the model name
  ///  @param level             model level of this model of the device
  ///  @param model_type_id     model type identifier
  ///  @param model_group_id    mode group identifier
  ///
  void addModel(const char *model_name, const int level, ModelTypeId model_type_id, ModelTypeId model_group_id);

private:
  ParametricData<void> &      instanceParameters_;            ///< Reference to the Instance specific parameter map
  ParametricData<void> &      modelParameters_;               ///< Reference to the Model specific parameter map
  const std::string           name_;                          ///< Name of the devive
  const std::string           deviceTypeName_;                ///< Type name of the device
  const std::string           instanceDefaultParameterName_;  ///< Default parameter name for the device instance
  int                         numNodes_;                      ///< Number of nodes of this device
  int                         numOptionalNodes_;              ///< Number of optional nodes of this device
  int                         numFillNodes_;                  ///< Number of fill nodes of this device
  bool                        modelRequired_;                 ///< True if model is required
  bool                        linearDevice_;                  ///< True if model is a linear device
  bool                        pdeDevice_;                     ///< True if model is a PDE device
  std::string                 primaryParameter_;              ///< Primary parameter name or the empty string if none
  std::vector<std::string>    modelTypeNames_;                ///< Vector of defined model type of this device's model
};

template<class T, class G>
struct ModelGroupTraits_
{
    typedef G GroupTraits_;
};

template<class T>
struct ModelGroupTraits_<T, ModelGroupLeadTraits>
{
    typedef T GroupTraits_;
};

//-----------------------------------------------------------------------------
// Class         : Config
// Purpose       :
// Special Notes :
// Scope         : protected
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Jan 27 15:06:50 2014
//-----------------------------------------------------------------------------
///
///  Config template derives from the Configuration class and provides the instance and model parameter decriptions.
///
template <class T>
class Config : public Configuration
{
  public:
    typedef T ModelTraits;
    typedef typename ModelGroupTraits_<T, typename T::ModelGroupTraits>::GroupTraits_ ModelGroupTraits;

    ///  Adds the device to the Xyce device configuration
    ///
    ///  @return reference to newly created configuration
    ///
    ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355
    ///  @date   Thu Feb  6 07:15:11 2014
    ///
    static Config<T> &addConfiguration()
    {
      return *(new Config<T>());
    }

  private:
    ///  Creates the configuration object, populating the instance and model parameter definitions from the Traits
    ///
    ///  The registerDevice() and registerModelType() functions are used by the device registration function to name the
    ///  device and device levels.
    ///
    ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355
    ///  @date   Mon Jan 27 15:28:59 2014
    ///
    Config()
      : Configuration(instanceParameters__, modelParameters__, ModelTraits::name(), ModelTraits::deviceTypeName(), ModelTraits::instanceDefaultParameter(),
                      ModelTraits::numNodes(), ModelTraits::numOptionalNodes(), ModelTraits::numFillNodes(), ModelTraits::modelRequired(),
                      ModelTraits::isLinearDevice(), ModelTraits::isPDEDevice(), ModelTraits::primaryParameter()),
        instanceParameters__(),
        modelParameters__()
    {
      ModelTraits::loadInstanceParameters(instanceParameters__);        /// loadInstanceParameters() is called populate the instance parameter desciptions
      ModelTraits::loadModelParameters(modelParameters__);              /// loadModelParameters() is called to populate the model parameter desciptions
    }

    //-----------------------------------------------------------------------------
    // Function      : ~Config
    // Purpose       :
    // Special Notes :
    // Scope         : public
    // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
    // Creation Date : Fri Mar 14 13:19:39 2014
    //-----------------------------------------------------------------------------
    ///
    /// Destroys the Config
    ///
    virtual ~Config()
    {}

    Config(const Config<T>  &);                 ///< No copying
    Config &operator=(const Config<T> &);       ///< No assignment

  public:
    //-----------------------------------------------------------------------------
    // Function      : registerDevice
    // Purpose       :
    // Special Notes :
    // Scope         : public
    // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
    // Creation Date : Mon Jan 27 15:31:43 2014
    //-----------------------------------------------------------------------------
    ///
    ///  Adds the device into the device registry giving it the specified device name and level
    ///
    ///  Registering a device maps the model type to the device name level key.
    ///
    ///  @param device_name       const character pointer name of the device
    ///  @param level             level of the device
    ///
    Config<T> &registerDevice(const char *device_name, const int level)
    {
      addDevice(device_name, level, ModelTraits::modelType(), ModelTraits::modelGroup(), ModelTraits::numNodes(), ModelGroupTraits::numNodes());

      return *this;
    }

    //-----------------------------------------------------------------------------
    // Function      : registerModelType
    // Purpose       :
    // Special Notes :
    // Scope         : public
    // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
    // Creation Date : Mon Jan 27 15:31:43 2014
    //-----------------------------------------------------------------------------
    ///
    ///  Adds a model type name and level to the device model type list
    ///
    ///  @param model_name        const character pointer model name
    ///  @param level             model level of the model
    ///
    ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355
    ///  @date   Mon Jan 27 15:32:50 2014
    ///
    Config<T> &registerModelType(const char *model_name, const int level)
    {
      addModel(model_name, level, ModelTraits::modelType(), ModelTraits::modelGroup());

      return *this;
    }

    //-----------------------------------------------------------------------------
    // Function      : instanceType
    // Purpose       :
    // Special Notes :
    // Scope         : public
    // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
    // Creation Date : Mon Jan 27 15:31:43 2014
    //-----------------------------------------------------------------------------
    ///
    ///  Returns the instance type identifier
    ///
    ///  @return the instance type identifier
    ///
    ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355
    ///  @date   Wed Feb  5 10:08:52 2014
    ///
    virtual InstanceTypeId instanceType() const /* override */
    {
      return ModelTraits::instanceType();
    }

    //-----------------------------------------------------------------------------
    // Function      : modelType
    // Purpose       :
    // Special Notes :
    // Scope         : public
    // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
    // Creation Date : Mon Jan 27 15:31:43 2014
    //-----------------------------------------------------------------------------
    ///
    ///  Returns the model type identifier
    ///
    ///  @return the model type identifier
    ///
    ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355
    ///  @date   Wed Feb  5 10:08:52 2014
    ///
    virtual ModelTypeId modelType() const /* override */
    {
      return ModelTraits::modelType();
    }

    //-----------------------------------------------------------------------------
    // Function      : modelGroup
    // Purpose       :
    // Special Notes :
    // Scope         : public
    // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
    // Creation Date : Mon Jan 27 15:31:43 2014
    //-----------------------------------------------------------------------------
    ///
    ///  Returns the model group identifier
    ///
    ///  @return the model group identifier
    ///
    virtual ModelTypeId modelGroup() const /* override */
    {
      return ModelTraits::modelGroup();
    }

    //-----------------------------------------------------------------------------
    // Function      : factory
    // Purpose       :
    // Special Notes :
    // Scope         : public
    // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
    // Creation Date : Mon Jan 27 15:31:43 2014
    //-----------------------------------------------------------------------------
    ///
    ///  Calls the device creation function
    ///
    ///  @param factory_block factory parameters
    ///
    ///  @return pointer to the new device
    ///
    ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355
    ///  @date   Mon Jan 27 15:35:19 2014
    ///
    Device *factory(const FactoryBlock &factory_block) const /* override */
    {
      return ModelTraits::factory(*this, factory_block);
    }

  private:
    ParametricData<typename ModelTraits::InstanceType>    instanceParameters__;           ///< Instance parameter description
    ParametricData<typename ModelTraits::ModelType>       modelParameters__;              ///< Model parameter description
};

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_Configuration_h

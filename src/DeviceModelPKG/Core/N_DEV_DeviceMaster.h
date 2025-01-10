//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : Provides templated versions of some boilerplate functions
//                  that are device-specific (so they can't easily be included
//                  in the base device, instance, or model classes).
//
// Special Notes  : Much of the functionality of the device classes, like
//                  N_DEV_Capacitor, is simply to manage STL containers
//                  of model and instance pointers.  That management is pretty
//                  generic, so templating that functionality makes sense.
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 01/31/06
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Device_Template_h
#define Xyce_N_DEV_Device_Template_h

#include <unordered_map>
using std::unordered_map;
#include <string>
#include <vector>

#include <N_DEV_Device.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_Configuration.h>

namespace Xyce {
namespace Device {

void instance_must_reference_model_error(const Device &device, const std::string &name, const NetlistLocation &netlist_location);
void could_not_find_model_error(const Device &device, const std::string &model_name, const std::string &instance_name, const NetlistLocation &netlist_location);
void duplicate_model_warning(const Device &device, const DeviceModel &model, const NetlistLocation &netlist_location);
void duplicate_instance_warning(const Device &device, const DeviceInstance &instance, const NetlistLocation &netlist_location);
void duplicate_entity_warning(const Device &device, const DeviceEntity &entity, const NetlistLocation &netlist_location);

//-----------------------------------------------------------------------------
// Class         : DeviceMaster
//
// Purpose       : This class contains a lot of boilerplate functions, such
//                 as addInstance, which are needed by every device.  However,
//                 as functions of this nature are device-specific, they
//                 didn't naturally fit as base class Device functions.
//
//                 This class is a template, which takes 2 typenames;
//                 model, instance and model-ptr STL vector.  A specific
//                 example of these would be CapacitorModel and
//                 CapacitorInstance
//
//                 Each derived device class will need to set up its own
//                 version of this template.
//
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/06
//-----------------------------------------------------------------------------

/**
 * DeviceMaster instantiates a device as described by the device traits T.
 *
 * The Device Traits are described in in the device configuration.
 *
 * @see Xyce::Device::Configuration
 *
 * @param T     device traits class
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Tue Feb  4 10:33:44 2014
 */
template<class T>
class DeviceMaster : public Device
{
public:
  typedef typename T::ModelType ModelType;            ///< Make the model begin defined available
  typedef typename T::InstanceType InstanceType;      ///< Make the instance being define available

protected:
  typedef std::vector<InstanceType *> InstanceVector;
  typedef unordered_map<std::string, ModelType *, HashNoCase, EqualNoCase> ModelMap;
  typedef unordered_map<std::string, InstanceType *, HashNoCase, EqualNoCase> InstanceMap;

public:
  /**
   *  Constructs a device
   *
   * When a device is constructed, a default model is created.
   *
   * @param configuration             const reference to device configuration
   * @param factory_block             const reference to the factory provided parameters
   * @param solver_state              const reference to the solver state to use for the device
   * @param device_options            const reference to the device options to use for the device
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:20:40 2014
   */
  DeviceMaster(
     const Configuration &     configuration,
     const FactoryBlock &      factory_block,
     const SolverState &       solver_state,
     const DeviceOptions &     device_options)
    : Device(),
      deviceName_(T::name()),
      defaultModelName_(std::string(T::deviceTypeName()) + " (" + T::name() + ")"),
      configuration_(configuration),
      solverState_(solver_state),
      deviceOptions_(device_options),
      modelMap_(),
      instanceVector_(),
      instanceMap_()
  {}

  /**
   * Destroys the device
   *
   * Delete the models created by this device.  The model will handle deleting the instances.
   *
   * @return
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:23:31 2014
   */
  virtual ~DeviceMaster()
  {
    for (typename ModelMap::iterator model_it = modelMap_.begin(); model_it != modelMap_.end(); ++model_it)
    {
      delete (*model_it).second;
    }
  }

private:
  DeviceMaster(const DeviceMaster &right);                ///< No copying
  DeviceMaster(const Device &);                             ///< Eh?
  DeviceMaster &operator=(const DeviceMaster &right);     ///< No assignments

public:
  /**
   * Returns the name of this device
   *
   * @return const reference to the device's name
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:25:37 2014
   */
  virtual const std::string &getName() const /* override */ 
  {
    return deviceName_;
  }

  /**
   * Returns the default model name to use if the instance being created does not specify one
   *
   * @return const reference to the default model name
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:26:10 2014
   */
  virtual const std::string &getDefaultModelName() const /* override */ 
  {
    return defaultModelName_;
  }

  /**
   * Returns a pointer to the model or model model with the specified name
   *
   * @param model_name       const reference to the model name
   *
   * @return pointer to the model or 0 is not found
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:36:44 2014
   */
  virtual DeviceModel *findModel(const ModelName &model_name) /* override */ 
  {
    typename ModelMap::iterator it = modelMap_.find(model_name);
    if (it != modelMap_.end())
      return (*it).second;

    return 0;
  }

  /**
   * Returns a pointer to the model or model model with the specified name
   *
   * @param model_name       const reference to the model name
   *
   * @return const pointer to the model or 0 is not found
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:36:44 2014
   */
  virtual const DeviceModel *findModel(const ModelName &model_name) const /* override */ 
  {
    typename ModelMap::const_iterator it = modelMap_.find(model_name);
    if (it != modelMap_.end())
      return (*it).second;

    return 0;
  }

  /**
   * Returns a pointer to the model or instance entity with the specified name
   *
   * @param instance_name       const reference to the entity name
   *
   * @return pointer to the entity or 0 is not found
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:36:44 2014
   */
  virtual DeviceEntity *findInstance(const InstanceName &instance_name) /* override */ 
  {
    typename InstanceMap::iterator it = instanceMap_.find(instance_name.getEncodedName());
    if (it != instanceMap_.end())
      return (*it).second;

    return 0;
  }

  /**
   * Returns a pointer to the model or instance entity with the specified name
   *
   * @param instance_name       const reference to the entity name
   *
   * @return const pointer to the entity or 0 is not found
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:36:44 2014
   */
  virtual const DeviceEntity *findInstance(const InstanceName &instance_name) const /* override */ 
  {
    typename InstanceMap::const_iterator it = instanceMap_.find(instance_name.getEncodedName());
    if (it != instanceMap_.end())
      return (*it).second;

    return 0;
  }

  /**
   * Returns true if this device is a linear device
   *
   * @return the device is linear value provided via the device trait.
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:51:52 2014
   */
  virtual bool isLinearDevice() const /* override */ 
  {
    return T::isLinearDevice();
  }

  /**
   * Returns true if this device is a PDE device
   *
   * @return the device is PDE value provided via the device trait.
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:53:10 2014
   */
  virtual bool isPDEDevice() const /* override */ 
  {
    return T::isPDEDevice();
  }

  /**
   * Executes operator op, passing its DeviceModel pointer, for each device model
   *
   * Since the model information is managed in a map of ModelType pointers, it is not possible to return the map as an
   * association to DeviceModel pointers.  This function iterates over the models and passes it the DeviceModel
   * pointer for each.
   *
   * @param op        functor taking a DeviceModel pointer as an argument
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:53:45 2014
   */
  virtual void forEachModel(DeviceModelOp &op) const /* override */ 
  {
    for (typename ModelMap::const_iterator it = modelMap_.begin(); it != modelMap_.end(); ++it)
      op((*it).second);
  }

  /**
   * Executes operator op, passing its DeviceInstance pointer, for each device instance
   *
   * Since the instance information is managed in a map of InstanceType pointers, it is not possible to return the map as an
   * association to DeviceInstance pointers.  This function iterates over the instances and passes it the DeviceInstance
   * pointer for each.
   *
   * @param op        functor taking a DeviceInstance pointer as an argument
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:53:45 2014
   */
  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */ 
  {
    for (typename InstanceVector::const_iterator it = instanceVector_.begin(); it != instanceVector_.end(); ++it)
      op(*it);
  }

  // virtual bool deleteInstance(DeviceInstance *instance); /* override */

  virtual ModelType *addModel(const ModelBlock & MB, const FactoryBlock &factory_block) /* override */;
  virtual InstanceType *addInstance(const InstanceBlock &instance_block, const FactoryBlock &factory_block); /* override */
  virtual void storeInstance(const FactoryBlock& factory_block, InstanceType* instance); /* override */

  virtual bool updateSources() /* override */;
  virtual bool updateState (double * solVec, double * staVec, double * stoVec)/* override */;
  virtual bool updateState (double * solVec, double * staVec, double * stoVec, int loadType) /* override */;
  virtual bool updateSecondaryState (double * staDerivVec, double * stoVec) /* override */;
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, 
                               double * leadF, double * leadQ, double * junctionV ) /* override */;
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, 
                               double * leadF, double * leadQ, double * junctionV, int loadType ) /* override */;
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx) /* override */;
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx, int loadType) /* override */;

  virtual bool updateFDIntermediateVars(double frequency,
                                        std::complex<double>* solVec)
  {
    return true;
  }

  virtual bool loadFreqDAEVectors(double frequency, std::complex<double>* solVec, 
                                  std::vector<Util::FreqVecEntry>& fVec,
                                  std::vector<Util::FreqVecEntry>& bVec) /* override */
  { 
    return true;
  }
  virtual bool loadFreqDAEMatrices(double frequency, std::complex<double>* solVec, 
                                   std::vector<Util::FreqMatEntry>& dFdx)  /* override */
  { 
    return true; 
  }

  virtual bool isConverged() const /* override */;

protected:
  /**
   * Returns the solver state given during device construction
   *
   * @return const reference to the solver state
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:27:14 2014
   */
  const SolverState &getSolverState() const 
  {
    return solverState_;
  }

  /**
   * Returns the device options given during device construction
   *
   * @return const reference to the device options
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:27:48 2014
   */
  const DeviceOptions &getDeviceOptions() const 
  {
    return deviceOptions_;
  }

  /**
   * Returns an iterator to the beginning of the vector of all instances created for this device
   *
   * While a device instance is created, the device model owns the pointer to the device.  The instanceVector_ gets a
   * copy so that all instances of this device may be iterated over without needing to go through the model.
   *
   * @return iterator to the beginning of the device instance vector
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:28:29 2014
   */
  typename InstanceVector::const_iterator getInstanceBegin() const 
  {
    return instanceVector_.begin();
  }

  /**
   * Returns an iterator to the ending of the vector of all instances created for this device
   *
   * While a device instance is created, the device model owns the pointer to the device.  The instanceVector_ gets a
   * copy so that all instances of this device may be iterated over without needing to go through the model.
   *
   * @return iterator to the ending of the device instance vector
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:28:29 2014
   */
  typename InstanceVector::const_iterator getInstanceEnd() const 
  {
    return instanceVector_.end();
  }

  /**
   * Returns true if the model name must be specified for each instance
   *
   * @return the device model required value provided via the device trait.
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:32:10 2014
   */
  bool isModelRequired() const 
  {
    return T::modelRequired();
  }

  // Separates the instance vector into linear and nonlinear instances.
  void separateInstanceTypes( InstanceVector& linearInstances, InstanceVector& nonlinearInstances ) const;


private:
    virtual bool getBreakPoints(std::vector<Util::BreakPoint> & breakPointTimes);

private:
  const std::string             deviceName_;
  const std::string             defaultModelName_;
  const Configuration &         configuration_;
  const SolverState &           solverState_;
  const DeviceOptions &         deviceOptions_;
  ModelMap                      modelMap_;
  InstanceVector                instanceVector_;
  InstanceMap                   instanceMap_;
};

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::addInstance
// Purpose       : This function adds an instance to the list of instances.
//
// Special Notes : All device instances are contained in one of many lists.
//                 Each model owns exactly one list of instances.  If
//                 an instance is added which does not reference a model, or
//                 that references a model that does not exist, then a new
//                 model has to be added.
//
//                 For now, this function will only add a model for instances
//                 that do not specifically reference one.  If an instance
//                 references a non-existant model, this function will return a
//                 fatal error.  (This function should only be called from the
//                 device manager function, AddDeviceInstance, and that
//                 function will not call this one if the instance refers to a
//                 non-existant model.)
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/06
//-----------------------------------------------------------------------------
template<class T>
typename DeviceMaster<T>::InstanceType *
DeviceMaster<T>::addInstance(
   const InstanceBlock &         instance_block,
   const FactoryBlock &          factory_block)
{
  std::string model_name = instance_block.getModelName();

  if (model_name.empty())  // If no model name, use the default model.
  {
    if (isModelRequired())
    {
      instance_must_reference_model_error(*this, model_name, instance_block.getNetlistLocation());
      return 0;
    }
    else
    {
      typename ModelMap::iterator it = modelMap_.find(defaultModelName_);
      if (it == modelMap_.end()) 
      {
        addModel(ModelBlock(defaultModelName_), factory_block);
      }
      model_name = defaultModelName_;
    }
  }

  typename ModelMap::iterator it = modelMap_.find(model_name);
  if (it == modelMap_.end()) 
  {
    could_not_find_model_error(*this, model_name, instance_block.getInstanceName().getDeviceName(), instance_block.getNetlistLocation());
    return 0;
  }

  ModelType &model = *(*it).second;

  std::pair<typename InstanceMap::iterator, bool> result = instanceMap_.insert(typename InstanceMap::value_type(instance_block.getInstanceName().getEncodedName(), typename InstanceMap::mapped_type(0)));
  if (result.second)
  {
    InstanceType *instance = new InstanceType(configuration_, instance_block, model, factory_block);
    instance->setDefaultParamName(T::instanceDefaultParameter());

    (*result.first).second = instance;

    model.addInstance(instance);

    storeInstance(factory_block, instance);

    if (modelMap_.find(instance_block.getInstanceName().getDeviceName()) != modelMap_.end())
      duplicate_entity_warning(*this, *instance, instance_block.getNetlistLocation());
  }
  else
  {
    duplicate_instance_warning(*this, *(*result.first).second, instance_block.getNetlistLocation());
  }

  return (*result.first).second;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::addModel
// Purpose       : This function adds a model to the list of device models.
//
// Special Notes : This is the  "public" version of the function, meaning that
//                 it returns a bool to indicate success or failure (rather
//                 than a pointer to the model).  Also, it checks to see if the
//                 model being added is already part of the list.  If the model
//                 already exists, it generates an error.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/06
//-----------------------------------------------------------------------------
template<class T>
typename DeviceMaster<T>::ModelType *
DeviceMaster<T>::addModel(
   const ModelBlock &    model_block,
   const FactoryBlock &  factory_block)
{
  // Check to make sure that there doesn't already exist a model of this name.
  std::pair<typename ModelMap::iterator, bool> result = modelMap_.insert(typename ModelMap::value_type(model_block.getName(), typename ModelMap::mapped_type(0)));
  if (result.second) 
  {
    ModelType *model = new ModelType(configuration_, model_block, factory_block);

    (*result.first).second = model;

    if (instanceMap_.find(model_block.getName()) != instanceMap_.end())
      duplicate_entity_warning(*this, *model, model_block.getNetlistLocation());
  }
  else
  {
    duplicate_model_warning(*this, *(*result.first).second, model_block.getNetlistLocation());
  }

  return (*result.first).second;
}

// //-----------------------------------------------------------------------------
// // Function      : DeviceMaster::deleteInstance
// // Purpose       :
// // Special Notes :
// // Scope         : public
// // Creator       : Eric Keiter, SNL
// // Creation Date : 2/03/06
// //-----------------------------------------------------------------------------
// template<class T>
// bool DeviceMaster<T>::deleteInstance(const std::string & instance_name)
// {
//   for (typename ModelMap::iterator it_model = modelMap_.begin(); it_model != modelMap_.end(); ++it_model)
//   {
//     InstanceVector &instance_list = (*it_model).second->getInstanceVector();

//     for (typename InstanceVector::iterator it_instance = instance_list.begin(); it_instance != instance_list.end(); ++it_instance)
//     {
//       if (equal_nocase(instance_name, (*it_instance)->getName()))
//       {
//         delete *it_instance;
//         instance_list.erase(it_instance);

//         return true;
//       }
//     }
//   }

//   return false;
// }

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::getBreakPoints
// Purpose       : This function adds break points to a vector of breakpoints.
//
// Special Notes : For most devices, this function won't get called.  The
//                 device manager only calls getBreakPoints for certain
//                 pre-selected devices.  Ultimately, this function probably
//                 isn't neccessary, as the device manager usually accesses
//                 device instances directly, rather than going through a
//                 device class middleman.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/05/06
//-----------------------------------------------------------------------------
template<class T>
bool DeviceMaster<T>::getBreakPoints( std::vector<Util::BreakPoint> & breakPointTimes )
{
  bool bsuccess = true;

  for (typename InstanceVector::iterator it = instanceVector_.begin(); it != instanceVector_.end(); ++it)
    bsuccess = (*it)->getInstanceBreakPoints(breakPointTimes) && bsuccess;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::updateSources
//
// Purpose       :
//
// Special Notes : Just like for the getBreakPoints function, for most devices, 
//                 this function won't get called.  The
//                 device manager only calls updateSources for certain
//                 pre-selected devices.  Ultimately, this function probably
//                 isn't neccessary, as the device manager usually accesses
//                 device instances directly, rather than going through a
//                 device class middleman.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 02/05/06
//-----------------------------------------------------------------------------
template<class T>
bool DeviceMaster<T>::updateSources()
{
  bool bsuccess = true;

  for (typename InstanceVector::iterator it = instanceVector_.begin(); it != instanceVector_.end(); ++it)
    bsuccess = (*it)->updateSource() && bsuccess;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
template<class T>
bool DeviceMaster<T>::updateState (double * solVec, double * staVec, double * stoVec)
{
  bool bsuccess = true;
  for (typename InstanceVector::const_iterator it = instanceVector_.begin(); it != instanceVector_.end(); ++it)
  {
    bool tmpBool = (*it)->updatePrimaryState();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
template<class T>
bool DeviceMaster<T>::updateState (double * solVec, double * staVec, double * stoVec, int loadType)
{
  // By default we call the original updateState.
  // Those devices that distinguish between load types can override this method.
  if ( (( loadType == LINEAR ) && isLinearDevice())
      || (( loadType == NONLINEAR ) && !isLinearDevice())
      || (( loadType == NONLINEAR_FREQ ) && !isLinearDevice())
      || (( loadType == PDE ) && isPDEDevice())
      || ( loadType == ALL ) )
  {
    return this->updateState( solVec, staVec, stoVec );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
template<class T>
bool DeviceMaster<T>::updateSecondaryState (double * staDerivVec, double * stoVec)
{
  bool bsuccess = true;
  for (typename InstanceVector::const_iterator it = instanceVector_.begin(); it != instanceVector_.end(); ++it)
  {
    bool tmpBool = (*it)->updateSecondaryState ();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
template<class T>
bool DeviceMaster<T>::loadDAEVectors(double * solVec, double * fVec, double * qVec, double * bVec, double * leadF, double * leadQ, double * junctionV)
{
  bool bsuccess = true;
  for (typename InstanceVector::const_iterator it = instanceVector_.begin(); it != instanceVector_.end(); ++it)
  {
    bool tmpBool = (*it)->loadDAEFVector();
    bsuccess = bsuccess && tmpBool;

    tmpBool = (*it)->loadDAEQVector();
    bsuccess = bsuccess && tmpBool;

    tmpBool = (*it)->loadDAEBVector();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
template<class T>
bool DeviceMaster<T>::loadDAEVectors(double * solVec, double * fVec, double * qVec, double * bVec, double * leadF, double * leadQ, double * junctionV, int loadType)
{
  // By default we call the original loadDAEVectors.
  // Those devices that distinguish between load types can override this method.
  if ( (( loadType == LINEAR ) && isLinearDevice()) 
      || (( loadType == NONLINEAR ) && !isLinearDevice())
      || (( loadType == NONLINEAR_FREQ ) && !isLinearDevice())
      || (( loadType == PDE ) && isPDEDevice())
      || ( loadType == ALL ) )
  {
    return this->loadDAEVectors( solVec, fVec, qVec, bVec, leadF, leadQ, junctionV );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
template<class T>
bool DeviceMaster<T>::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  bool bsuccess = true;
  for (typename InstanceVector::const_iterator it = instanceVector_.begin(); it != instanceVector_.end(); ++it)
  {
    bool tmpBool = (*it)->loadDAEdFdx();
    bsuccess = bsuccess && tmpBool;
    tmpBool = (*it)->loadDAEdQdx();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
template<class T>
bool DeviceMaster<T>::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx, int loadType)
{
  // By default we call the original loadDAEMatrices.
  // Those devices that distinguish between load types can override this method.
  if ( (( loadType == LINEAR ) && isLinearDevice()) 
      || (( loadType == NONLINEAR ) && !isLinearDevice())
      || (( loadType == NONLINEAR_FREQ ) && !isLinearDevice())
      || (( loadType == PDE ) && isPDEDevice())
      || ( loadType == ALL ) )
  {
    return this->loadDAEMatrices( dFdx, dQdx );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::isConverged
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 7/14/15
//-----------------------------------------------------------------------------
template<class T>
bool DeviceMaster<T>::isConverged() const
{
  bool bsuccess = true;
  
  if ( !isLinearDevice() )
  {
    for (typename InstanceVector::const_iterator it = instanceVector_.begin(); it != instanceVector_.end(); ++it)
    {
      bool tmpBool = (*it)->isConverged();
      bsuccess = bsuccess && tmpBool;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::separateInstanceTypes()
// Purpose       : loop over instance vector and separate linear and nonlinear instances.
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 7/18/2014
//-----------------------------------------------------------------------------
template<class T>
void DeviceMaster<T>::separateInstanceTypes( InstanceVector& linearInstances,
                                             InstanceVector& nonlinearInstances ) const
{
  for (typename InstanceVector::const_iterator it = instanceVector_.begin(); it != instanceVector_.end(); ++it)
  {
    if ((*it)->isLinearDevice())
    {
      linearInstances.push_back( *it );
    }
    else
    {
      nonlinearInstances.push_back( *it );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::storeInstance()
// Purpose       : store new instance in local vector storage, this is used for 
//               : instance type separation (i.e. linear vs. nonlinear)
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 8/1/2017
//-----------------------------------------------------------------------------
template<class T>
void DeviceMaster<T>::storeInstance( const FactoryBlock& factory_block, InstanceType* instance )
{
  instanceVector_.push_back(instance);
}

} // namespace Device
} // namespace Xyce

#endif

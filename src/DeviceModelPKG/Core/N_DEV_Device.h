//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Device_h
#define Xyce_N_DEV_Device_h

#include <iosfwd>
#include <map>
#include <string>
#include <vector>
#include <functional>

#include <N_DEV_fwd.h>
#include <N_LAS_fwd.h>
#include <N_UTL_fwd.h>
#include <N_DEV_InstanceName.h>

namespace Xyce {
namespace Device {

struct DeviceModelOp
{
  using result_type = bool;
  using first_argument_type = DeviceModel*;

  virtual ~DeviceModelOp()
  {}

  virtual bool operator()(DeviceModel *model) = 0;
};

struct DeviceInstanceOp
{
  using result_type = bool;
  using first_argument_type = DeviceInstance*;

  virtual ~DeviceInstanceOp()
  {}

  virtual bool operator()(DeviceInstance *instance) = 0;
};

//-----------------------------------------------------------------------------
// Class         : Device
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Apr 22 12:47:22 2014
//-----------------------------------------------------------------------------
///
///  The Device class is an interface for device implementations.
///
///  In general, DeviceMaster is the only class that actually inherits
///  from Device.  Most devices either use the DeviceMaster class, but
///  some derive from DeviceMaster in what is known as a Master class.
///
///  The interfaces are unfortunately used for but derived
///  implementation as well as usage.  This should be changed at some
///  point, but this is the currently basic design of most interface
///  classes.
///
///  @author Eric Keiter, SNL, Parallel Computational Sciences
///  @date   3/16/00
class Device
{
public:
  Device()
  {}

  virtual ~Device()
  {}

private:
  Device(const Device &);                     ///< No copying
  Device &operator=(const Device &);          ///< No assignment

public:
  //-----------------------------------------------------------------------------
  // Function      : isLinearDevice
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 16:36:27 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Returns true if the device is linear
  ///
  ///  @return  true if the device is linear.
  ///
  virtual bool isLinearDevice() const = 0;

  //-----------------------------------------------------------------------------
  // Function      : isPDEDevice
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 16:37:02 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Returns true is the device is a PDE device
  ///
  ///  @return true is the device is a PDE device.
  ///
  virtual bool isPDEDevice() const = 0;

  //-----------------------------------------------------------------------------
  // Function      : getName
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 16:37:34 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Returns the name given to the device
  ///
  ///  @return const reference to the device name.
  ///
  virtual const std::string &getName() const = 0;

  //-----------------------------------------------------------------------------
  // Function      : getDefaultModelName
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 16:38:02 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Returns the name of the default model that would to used for this device
  ///
  ///  @return const reference to the name of the default model
  ///
  virtual const std::string &getDefaultModelName() const = 0;

  //-----------------------------------------------------------------------------
  // Function      : findModel
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 16:39:23 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Returns the device entity with the specified name
  ///
  ///  @param model_name       const reference to the name of the entity
  ///
  ///  @return pointer to the device entity
  ///
  virtual DeviceModel *findModel(const ModelName &model_name) = 0;

  //-----------------------------------------------------------------------------
  // Function      : findModel
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 16:39:23 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Returns the device entity with the specified name
  ///
  ///  @param model_name       const reference name of the entity
  ///
  ///  @return const pointer to the device entity
  ///
  virtual const DeviceModel *findModel(const ModelName &model_name) const = 0;

  //-----------------------------------------------------------------------------
  // Function      : findInstance
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 16:39:23 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Returns the device entity with the specified name
  ///
  ///  @param instance_name       const reference to the name of the entity
  ///
  ///  @return pointer to the device entity
  ///
  virtual DeviceEntity *findInstance(const InstanceName &instance_name) = 0;

  //-----------------------------------------------------------------------------
  // Function      : findInstance
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 16:39:23 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Returns the device entity with the specified name
  ///
  ///  @param instance_name       const reference name of the entity
  ///
  ///  @return const pointer to the device entity
  ///
  virtual const DeviceEntity *findInstance(const InstanceName &instance_name) const = 0;

  //-----------------------------------------------------------------------------
  // Function      : addModel
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 16:40:30 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Creates a device model and adds it to the device's list of models
  ///
  ///  @param model_block       const reference to the model block describing the model to create
  ///  @param factory_block     const reference to the factory data needed to create the model
  ///
  ///  @return pointer to the newly created device model
  ///
  virtual DeviceModel *addModel(const ModelBlock &model_block, const FactoryBlock &factory_block) = 0;

  //-----------------------------------------------------------------------------
  // Function      : addInstance
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 16:43:08 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Creates a device instance and adds to the device model's instance list
  ///
  ///  @param instance_block    const reference to the model block describing the instance to create
  ///  @param factory_block     const reference to the factory data needed to create the instance
  ///
  ///  @return pointer to the newly creates device instance
  ///
  virtual DeviceInstance *addInstance(const InstanceBlock &instance_block, const FactoryBlock &factory_block) = 0;

  //-----------------------------------------------------------------------------
  // Function      : updateSources
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 16:44:34 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Updates the devices source information
  ///
  ///  This function is called by the analysis subsystem when it is time to update the device source information.
  ///
  ///  @return true if the update was successful
  ///
  virtual bool updateSources()
  {
    return true;
  }

  //-----------------------------------------------------------------------------
  // Function      : updateState
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 16:44:34 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Updates the devices state information
  ///
  ///  This function is called by the analysis subsystem when it is time to update the device state information.
  ///
  ///  @return true if the update was successful
  ///
  virtual bool updateState(double * solVec, double * staVec, double * stoVec)
  {
    return true;
  }

  //-----------------------------------------------------------------------------
  // Function      : updateState
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 16:44:34 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Updates the devices state information
  ///
  ///  This function is called by the analysis subsystem when it is time to update the device state information.
  ///
  ///  @return true if the update was successful
  ///
  virtual bool updateState(double * solVec, double * staVec, double * stoVec, int loadType)
  {
    return true;
  }

  //-----------------------------------------------------------------------------
  // Function      : updateSecondaryState
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 16:44:34 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Updates the devices secondary state information
  ///
  ///  This function is called by the analysis subsystem when it is time to update the device secondary state information.
  ///
  ///  @return true if the update was successful
  ///
  virtual bool updateSecondaryState(double * staDerivVec, double * stoVec)
  {
    return true;
  }

  //-----------------------------------------------------------------------------
  // Function      : loadDAEVectors
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 16:47:20 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Populates the device's ExternData object with these pointers
  ///
  ///  THIS FUNCTION MUST BE CALLED PRIOR TO CALLING loadDAEMatrices.
  ///
  ///  @param solVec            pointer to the analysis solution vector for this device
  ///  @param fVec              pointer to the analysis f vector for this device
  ///  @param qVec              pointer to the analysis q vector for this device
  ///  @param leadF             pointer to the analysis lead f vector for this device
  ///  @param leadQ             poitner to the analysis load q vector for this device
  ///
  ///  @return true if the update was successful
  ///
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double *bVec, double * leadF, double * leadQ, double * junctionV)
  {
    return true;
  }

  //-----------------------------------------------------------------------------
  // Function      : loadDAEVectors
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 16:47:20 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Populates the device's ExternData object with these pointers
  ///
  ///  THIS FUNCTION MUST BE CALLED PRIOR TO CALLING loadDAEMatrices.
  ///
  ///  @param solVec            pointer to the analysis solution vector for this device
  ///  @param fVec              pointer to the analysis f vector for this device
  ///  @param qVec              pointer to the analysis q vector for this device
  ///  @param leadF             pointer to the analysis lead f vector for this device
  ///  @param leadQ             poitner to the analysis load q vector for this device
  ///
  ///  @return true if the update was successful
  ///
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, double * leadF, double * leadQ, double * junctionV, int loadType ) 
  {
    return true;
  }

  //-----------------------------------------------------------------------------
  // Function      : loadDAEMatrices
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 17:06:52 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Populates the device's Jacobian object with these pointers
  ///
  ///  THIS FUNCTION MUST BE CALLED AFTER CALLING loadDAEVectors.
  ///
  ///  @param dFdx      pointer to the analysis dFdx matrix
  ///  @param dQdx      pointer to the analysis dQdx matrix
  ///
  ///  @return
  ///
  virtual bool loadDAEMatrices(Linear::Matrix & dFdx, Linear::Matrix & dQdx)
  {
    return true;
  }

  //-----------------------------------------------------------------------------
  // Function      : loadDAEMatrices
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 17:06:52 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Populates the device's Jacobian object with these pointers
  ///
  ///  THIS FUNCTION MUST BE CALLED AFTER CALLING loadDAEVectors.
  ///
  ///  @param dFdx      pointer to the analysis dFdx matrix
  ///  @param dQdx      pointer to the analysis dQdx matrix
  ///
  ///  @return
  ///
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx, int loadType)
  {
    return this->loadDAEMatrices( dFdx, dQdx );
  }

  //----------------------------------------------------------------------------
  // Function      : updateFDIntermediateVars
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Tom Russo and Ting Mei
  // Creation Date : 19 March 2018
  //----------------------------------------------------------------------------  /// Update device internal variables in frequency domain
  ///
  /// \author Tom Russo and Ting Mei
  /// \date 19 March 2018
  ///
  virtual bool updateFDIntermediateVars(double frequency,
                                        std::complex<double>* solVec)
  {
    return true;
  }

  //----------------------------------------------------------------------------
  // Function      : loadFreqDAEMatrices
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Jason C. Verley
  // Creation Date : 07/13/17
  //----------------------------------------------------------------------------  /// Load frequency domain jacobian matrices
  ///
  /// Loads the DAE contributions to the Jacobian matrix.  To be used by
  /// frequency-based analysis techniques.  Overridden by devices that have a
  /// frequency reponse, so this does nothing.
  ///
  /// \author Jason C. Verley
  /// \date 07/13/17
  ///
  virtual bool loadFreqDAEVectors(double frequency, std::complex<double>* solVec, 
                                  std::vector<Util::FreqVecEntry>& fVec,
                                  std::vector<Util::FreqVecEntry>& bVec)
  { 
    return true; 
  }

  //----------------------------------------------------------------------------
  // Function      : loadFreqDAEMatrices
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Jason C. Verley
  // Creation Date : 07/13/17
  //----------------------------------------------------------------------------
  /// Load frequency domain F and B vectors
  ///
  /// Loads the contributions to the `F` and `B` vectors.  To be used by
  /// frequency-based analyses techniques.  Overridden by devices that have a
  /// frequency reponse, so this does nothing.
  ///
  /// \author Jason C. Verley
  /// \date 07/13/17
  ///
  ///////////////////////////////////////////////////////////////////////////////
  virtual bool loadFreqDAEMatrices(double frequency, std::complex<double>* solVec, 
                                   std::vector<Util::FreqMatEntry>& dFdx)
  { 
    return true; 
  }

  //-----------------------------------------------------------------------------
  // Function      : isConverged 
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 17:06:52 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Checks that the device is converged
  //
  virtual bool isConverged() const
  {
    return true;
  }

  //-----------------------------------------------------------------------------
  // Function      : deleteInstance
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Wed Jan 29 17:08:58 2014
  //-----------------------------------------------------------------------------
  // ///
  // ///  Delete the specified instance from the device
  // ///
  // ///  TODO: Check this implementation in DeviceMaster...does it really work?
  // ///
  // ///  @param instance_name     const reference to the instance name of the device instance to be deleted
  // ///
  // ///  @return true if the delete was successful
  // ///
  // virtual bool deleteInstance(DeviceInstance *instance) = 0;

  //-----------------------------------------------------------------------------
  // Function      : forEachModel
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon Feb  3 10:36:04 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Executes op on each DeviceModel pointer of the device
  ///
  ///  To use this function, create a class which derives from DeviceModelOp.  Then implement the
  ///  operator()(DeviceModel/// ).  Sample classes are defined in N_DEV_Algorithm.C.
  ///
  ///  @param op        reference to operator functor
  ///
  virtual void forEachModel(DeviceModelOp &op) const = 0;

  //-----------------------------------------------------------------------------
  // Function      : forEachInstance
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon Feb  3 10:36:04 2014
  //-----------------------------------------------------------------------------
  ///
  ///  Executes op on each DeviceInstance pointer of the device
  ///
  ///  To use this function, create a class which derives from DeviceInstanceOp.  Then implement the
  ///  operator()(DeviceInstance *).  Sample classes are defined in N_DEV_Algorithm.C.
  ///
  ///  @param op        reference to operator functor
  ///
  virtual void forEachInstance(DeviceInstanceOp &op) const = 0;
};

} // namespace Device
} // namespace Xyce

#endif

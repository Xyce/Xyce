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

//-----------------------------------------------------------------------------
//
// Purpose        : This file contains the device entity base class.
//
// Special Notes  :
//
// Creator        : Dave Shirley, PSSI
//
// Creation Date  : 05/05/05
//
//
//
//
//-----------------------------------------------------------------------------
#ifndef Xyce_N_DEV_CompositeParam_h
#define Xyce_N_DEV_CompositeParam_h

#include <string>

#include <N_DEV_Pars.h>

namespace Xyce {
namespace Device {

///
/// CompositeParam is the base class for classes that wish to only manage the processing of parameter data.
///
/// The DeviceEntity class is vary similar, except that it manages a device as well as the device's parameter data.
/// During DeviceEntity's processing of parameters, it may create several object of classes derived from CompositeParam
/// which hold the processes parametric data.  The parametricData_ member holds the parameter descriptions the
/// Device::setParameters() function populates the values in the object of the derived class while the processParams()
/// virtual function can handle any additional processing of the parameters after they have been set.
/// parametricData_ object.
///
/// See Device::populateParams() and Device::setParameters() in the DeviceEntity implementation file.
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
/// @date   Wed Jan 29 17:26:35 2014
///
class CompositeParam : public ParameterBase
{
public:
  ///
  /// CompositeParam sets the parametric data description.
  ///
  ///
  /// @param parametric_data   reference to the parametric data description
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  /// @date   Wed Jan 29 17:33:22 2014
  /// 
  CompositeParam(ParametricData<void> &parametric_data)
    : parametricData_(parametric_data)
  {}

  virtual ~CompositeParam()
  {}

private:
  CompositeParam(const CompositeParam &);                     ///< No copying
  CompositeParam &operator=(const CompositeParam &);          ///< No assignment

public:
  ///
  /// processParams post processes the parameters that have been set in the object of the derived class.
  ///
  /// See Device::populateParams() and Device::setParameters() in the DeviceEntity implementation file.
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  /// @date   Wed Jan 29 17:34:53 2014
  ///
  virtual void processParams() = 0;

  bool given(const std::string &parameter_name) const;

  
  /// getParameterMap returns the parameter map which describes the parameters.
  /// 
  /// @return reference to the parameter map
  ///   
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  /// @date   Wed Jan 29 17:49:10 2014
  ///
  const ParameterMap &getParameterMap() const 
  {
    return parametricData_.getMap();
  }

private:
  ParametricData<void> &      parametricData_;                ///< Parameter data desciptions
};

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_CompositeParam_h

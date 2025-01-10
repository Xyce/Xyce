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
// Purpose        :
//
// Special Notes  :
//
// Creator        : David Baur
//
// Creation Date  : 3/28/2013
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <ostream>
#include <string>

#include <N_DEV_Pars.h>

#include <N_DEV_DeviceOptions.h>
#include <N_ERH_Message.h>
#include <N_UTL_Demangle.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : setDefaultParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Mar 10 11:14:01 2014
//-----------------------------------------------------------------------------
///
/// Iterates over parameters of parameter_base object, setting parameter_base members variables to the default value
/// provided when the addPar() function created the parameter descrption
///
/// @param parameter_base       DeviceEntity or CompositeParam
/// @param begin                Begin iterator of parameter descriptions
/// @param end                  End iterator of parameters descriptions
/// @param device_options       Device options
///
void setDefaultParameters(ParameterBase &parameter_base, ParameterMap::const_iterator begin, ParameterMap::const_iterator end, const DeviceOptions &device_options)
{
// First, allocate and zero out original and given vals
  for (ParameterMap::const_iterator it = begin ; it != end ; ++it)
  {
    Descriptor &param = *(*it).second;
    Xyce::Device::setValueGiven(parameter_base, param.getSerialNumber(), false);
    if (param.hasGivenMember())
      param.setGiven(parameter_base, false);

    if (param.isType<double>()) 
    {
      if (param.getExpressionAccess() & MIN_RES)
      {
        setDefaultValue<double>(param, device_options.minRes);
      }
      else if (param.getExpressionAccess() & MIN_CAP)
      {
        setDefaultValue<double>(param, device_options.minCap);
      }
      param.value<double>(parameter_base) = getDefaultValue<double>(param);
    }
    else if (param.isType<bool>())
      param.value<bool>(parameter_base) = getDefaultValue<bool>(param);
    else if (param.isType<int>())
      param.value<int>(parameter_base) = getDefaultValue<int>(param);
    else if (param.isType<long>())
      param.value<long>(parameter_base) = getDefaultValue<long>(param);
    else if (param.isType<std::string>())
      param.value<std::string>(parameter_base) = getDefaultValue<std::string>(param);
    else if (param.isType<std::vector<int> >())
      (param.value<std::vector<int> >(parameter_base)).clear();
    else if (param.isType<std::vector<double> >())
      (param.value<std::vector<double> >(parameter_base)).clear();
    else if (param.isType<std::vector<std::string> >())
      (param.value<std::vector<std::string> >(parameter_base)).clear();
  }
}

//-----------------------------------------------------------------------------
// Function      : nonexistentParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Mar 10 11:28:31 2014
//-----------------------------------------------------------------------------
///
/// Report casting error when attempting to cast from from_type to to_type
///
/// @param name          parameter name
/// @param entity_type   entity
///
void nonexistentParameter(const std::string &name, const std::type_info &entity_type)
{
  Report::DevelFatal0() << "Parameter " << name << " does not exist in " << demangle(entity_type.name());
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::typeMismatch
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : David Baur
// Creation Date : 8/6/2014
//-----------------------------------------------------------------------------
///
/// Report casting error when attempting to cast from from_type to to_type
///
/// @param from_type Typeinfo casting from
/// @param to_type Typeinfo casting to
///
void typeMismatch(const std::type_info &from_type, const std::type_info &to_type)
{
  Report::DevelFatal0() << "Attempting to cast parameter of type " << demangle(from_type.name()) << " to type " << demangle(to_type.name());
}

//-----------------------------------------------------------------------------
// Function      : checkExprAccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Mar 10 11:28:53 2014
//-----------------------------------------------------------------------------
///
/// Report error if both MIN_CAP and MIN_RES have been specified.
///
/// @param name Parameter name
/// @param expr_access Parameter expr access to verify
/// @param parameter_data_class Typeinfo to display name of class on error
///
void checkExprAccess(const std::string &name, ParameterType::ExprAccess &expr_access, const std::type_info &parameter_data_class)
{
  if ((expr_access & ParameterType::MIN_CAP) && (expr_access & ParameterType::MIN_RES))
    Report::DevelFatal0() << "Attempt to set MIN_CAP and MIN_RES on ParameterType::ExprAccess for parameter " << name << " in class " << parameter_data_class.name();
}


//-----------------------------------------------------------------------------
// Function      : ParametricData<void>::addDescriptor
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Mar 10 11:22:43 2014
//-----------------------------------------------------------------------------
///
/// Adds an entry to descriptor map.
///
/// The serial number of the parameter is set in the descriptor.
///
/// @invariant descriptor serial number is unique for each descriptor
/// @invariant parameter descriptor has unique name
///
/// @param name                 Name or parameter to declare
/// @param descriptor           Type information of parameter
/// @param parameter_data_class Typeinfo to display name of class on error
///
void ParametricData<void>::addDescriptor(const std::string &name, Descriptor *descriptor, const std::type_info &parameter_data_class)
{
  descriptor->setSerialNumber(map_.size());

  std::pair<ParameterMap::iterator, bool> result = map_.insert(ParameterMap::value_type(name, descriptor));

  if (!result.second)
    Report::DevelFatal0() << "Parameter " << name << " already added to class " << demangle(parameter_data_class.name());
}

} // namespace Device
} // namespace Xyce

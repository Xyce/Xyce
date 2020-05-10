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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/03/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <string>
#include <iostream>
#if defined(HAVE_UNORDERED_MAP)
#include <unordered_map>
using std::unordered_map;
#elif defined(HAVE_TR1_UNORDERED_MAP)
#include <tr1/unordered_map>
using std::tr1::unordered_map;
#else
#error neither unordered_map or tr1/unordered_map found
#endif
#include <N_DEV_fwd.h>
#include <N_DEV_CompositeParam.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_Message.h>
#include <N_DEV_Param.h>
#include <N_DEV_SolverState.h>
#include <N_LAS_Vector.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_Expression.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::DeviceEntity
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 12/13/04
//-----------------------------------------------------------------------------
DeviceEntity::DeviceEntity(
  ParametricData<void> &        parametric_data,
  const SolverState &           solver_state,
  const DeviceOptions &         device_options,
  const NetlistLocation &       netlist_location)
  : defaultParamName_(),
    parametricData_(parametric_data),
    netlistLocation_(netlist_location),
    solState_(solver_state),
    globals_(solver_state.getGlobals()),
    devOptions_(device_options)
{}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::~DeviceEntity
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/30/00
//-----------------------------------------------------------------------------
DeviceEntity::~DeviceEntity()
{
  for (std::vector<Depend>::iterator d = dependentParams_.begin(), end = dependentParams_.end(); d != end; ++d)
  {
    delete d->expr;
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::scaleParam
//
// Purpose : Scales the original value of the specified parameter by the specified value.
// The parameter is never specified by the user so errors are developer caused.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/21/04
//-----------------------------------------------------------------------------
bool DeviceEntity::scaleParam( const std::string & paramName, double val, double val0)
{
  ParameterMap::const_iterator p_i = getParameterMap().find(paramName);
  if (p_i == getParameterMap().end())
  {
    DevelFatal(*this).in("DeviceEntity::scaleParam") << "Unrecognized parameter " << paramName;
    return false;
  }

  const Descriptor &param = *(*p_i).second;
  if (!param.hasOriginalValueStored())
  {
    DevelFatal(*this).in("DeviceEntity::scaleParam") << "Original value not available for parameter " << paramName;
    return false;
  }

  if (!param.isType<double>())
  {
    DevelFatal(*this).in("DeviceEntity::scaleParam") << "Can scale only double parameters, parameter " << paramName << " is not double";
    return false;
  }

  // Scale the parameter
  setValue<double, DeviceEntity>(*this, param, Xyce::Device::getOriginalValue(*this, param.getSerialNumber())*val + val0*(1.0-val));

  if (param.hasGivenMember())
    param.setGiven(*this, true);

  Xyce::Device::setValueGiven(*this, param.getSerialNumber(), true);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::scaleParam
//
// Purpose       : Scales the specified parameter by a specified value.
// The parameter is never specified by the user so errors are developer caused.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/21/04
//-----------------------------------------------------------------------------
bool DeviceEntity::scaleParam( const std::string & paramName, double val)
{
  ParameterMap::const_iterator p_i = getParameterMap().find(paramName);
  if (p_i == getParameterMap().end())
  {
    DevelFatal(*this).in("DeviceEntity::scaleParam") << "Unrecognized parameter " << paramName;
    return false;
  }

  const Descriptor &param = *(*p_i).second;
  if (!param.hasOriginalValueStored())
  {
    DevelFatal(*this).in("DeviceEntity::scaleParam") << "Original value not available for parameter " << paramName;
    return false;
  }

  if (!param.isType<double>())
  {
    DevelFatal(*this).in("DeviceEntity::scaleParam") << "Can scale only double parameters, parameter " << paramName << " is not double";
    return false;
  }

  // Scale the parameter
  param.value<double>(*this) = Xyce::Device::getOriginalValue(*this, param.getSerialNumber())*val;

  if (param.hasGivenMember())
    param.setGiven(*this, true);

  Xyce::Device::setValueGiven(*this, param.getSerialNumber(), true);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::scaleDefaultParam
// Purpose       :
// The parameter is never specified by the user so errors are developer caused.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/19/05
//-----------------------------------------------------------------------------
bool DeviceEntity::scaleDefaultParam(double val)
{
  if (defaultParamName_.empty())
  {
    DevelFatal(*this).in("DeviceEntity::scaleDefaultParam") << "Device does not have a default parameter";
    return false;
  }

  return scaleParam(defaultParamName_, val);
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::analyticSensitivityAvailable
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 7/17/2014
//-----------------------------------------------------------------------------
bool DeviceEntity::analyticSensitivityAvailable (const std::string & paramName)
{
  ParameterMap::const_iterator p_i = getParameterMap().find(paramName);
  if (p_i == getParameterMap().end())
  {
    DevelFatal(*this).in("DeviceEntity::analyticSensitivityAvailable") << "Unrecognized parameter " << paramName;
    return false;
  }

  const Descriptor &param = *(*p_i).second;
  return param.getAnalyticSensitivityAvailable();
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::analyticSensitivityAvailableDefaultParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/16/2015
//-----------------------------------------------------------------------------
bool DeviceEntity::analyticSensitivityAvailableDefaultParam ()
{
  if (defaultParamName_.empty())
  {
    DevelFatal(*this).in("DeviceEntity::analyticSensitivityAvailableDefaultParam") 
      << "Device does not have a default parameter";
    return false;
  }

  return analyticSensitivityAvailable(defaultParamName_);
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::getAnalyticSensitivity
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 7/17/2014
//-----------------------------------------------------------------------------
bool DeviceEntity::getAnalyticSensitivity ( const std::string & paramName,
                                            std::vector<double> & dfdpVec,
                                            std::vector<double> & dqdpVec,
                                            std::vector<double> & dbdpVec,
                                            std::vector<int> & FindicesVec,
                                            std::vector<int> & QindicesVec,
                                            std::vector<int> & BindicesVec)
{
  ParameterMap::const_iterator p_i = getParameterMap().find(paramName);
  if (p_i == getParameterMap().end())
  {
    DevelFatal(*this).in("DeviceEntity::analyticSensitivityAvailable") << "Unrecognized parameter " << paramName;
    return false;
  }

  const Descriptor &param = *(*p_i).second;

  return param.getAnalyticSensitivity (*this, paramName, dfdpVec, dqdpVec, dbdpVec,
      FindicesVec, QindicesVec, BindicesVec);
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::getAnalyticSensitivityDefaultParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 7/17/2014
//-----------------------------------------------------------------------------
bool DeviceEntity::getAnalyticSensitivityDefaultParam (
                                            std::vector<double> & dfdpVec,
                                            std::vector<double> & dqdpVec,
                                            std::vector<double> & dbdpVec,
                                            std::vector<int> & FindicesVec,
                                            std::vector<int> & QindicesVec,
                                            std::vector<int> & BindicesVec)
{
  if (defaultParamName_.empty())
  {
    DevelFatal(*this).in("DeviceEntity::getAnalyticSensitivityDefaultParam") 
      << "Device does not have a default parameter";
    return false;
  }

  return getAnalyticSensitivity(defaultParamName_,
      dfdpVec, dqdpVec, dbdpVec, FindicesVec, QindicesVec, BindicesVec);
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::getNumericalSensitivityDefaultParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool DeviceEntity::getNumericalSensitivityDefaultParam (
                                            std::vector<double> & dfdpVec,
                                            std::vector<double> & dqdpVec,
                                            std::vector<double> & dbdpVec,
                                            std::vector<int> & FindicesVec,
                                            std::vector<int> & QindicesVec,
                                            std::vector<int> & BindicesVec)
{
  if (defaultParamName_.empty())
  {
    DevelFatal(*this).in("DeviceEntity::getNumericalSensitivityDefaultParam") 
      << "Device does not have a default parameter";
    return false;
  }

  return getNumericalSensitivity(defaultParamName_,
      dfdpVec, dqdpVec, dbdpVec, FindicesVec, QindicesVec, BindicesVec);
}

//

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::analyticMatrixSensitivityAvailable
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 7/17/2014
//-----------------------------------------------------------------------------
bool DeviceEntity::analyticMatrixSensitivityAvailable (const std::string & paramName)
{
  ParameterMap::const_iterator p_i = getParameterMap().find(paramName);
  if (p_i == getParameterMap().end())
  {
    DevelFatal(*this).in("DeviceEntity::analyticMatrixSensitivityAvailable") << "Unrecognized parameter " << paramName;
    return false;
  }

  const Descriptor &param = *(*p_i).second;
  return param.getAnalyticMatrixSensitivityAvailable();
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::analyticMatrixSensitivityAvailableDefaultParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/16/2015
//-----------------------------------------------------------------------------
bool DeviceEntity::analyticMatrixSensitivityAvailableDefaultParam ()
{
  if (defaultParamName_.empty())
  {
    DevelFatal(*this).in("DeviceEntity::analyticMatrixSensitivityAvailableDefaultParam") 
      << "Device does not have a default parameter";
    return false;
  }

  return analyticMatrixSensitivityAvailable(defaultParamName_);
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::getAnalyticMatrixSensitivity
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 7/17/2014
//-----------------------------------------------------------------------------
bool DeviceEntity::getAnalyticMatrixSensitivity ( const std::string & paramName,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs )
{
  ParameterMap::const_iterator p_i = getParameterMap().find(paramName);
  if (p_i == getParameterMap().end())
  {
    DevelFatal(*this).in("DeviceEntity::analyticMatrixSensitivityAvailable") << "Unrecognized parameter " << paramName;
    return false;
  }

  const Descriptor &param = *(*p_i).second;

  return param.getAnalyticMatrixSensitivity (*this, paramName, 
    d_dfdx_dp, d_dqdx_dp, F_lids, Q_lids, F_jacLIDs, Q_jacLIDs);
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::getAnalyticMatrixSensitivityDefaultParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 7/17/2014
//-----------------------------------------------------------------------------
bool DeviceEntity::getAnalyticMatrixSensitivityDefaultParam (
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs )
{
  if (defaultParamName_.empty())
  {
    DevelFatal(*this).in("DeviceEntity::getAnalyticMatrixSensitivityDefaultParam") 
      << "Device does not have a default parameter";
    return false;
  }

  return getAnalyticMatrixSensitivity(defaultParamName_,
    d_dfdx_dp, d_dqdx_dp, F_lids, Q_lids, F_jacLIDs, Q_jacLIDs);
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::getNumericalMatrixSensitivityDefaultParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool DeviceEntity::getNumericalMatrixSensitivityDefaultParam (
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs )
{
  if (defaultParamName_.empty())
  {
    DevelFatal(*this).in("DeviceEntity::getNumericalMatrixSensitivityDefaultParam") 
      << "Device does not have a default parameter";
    return false;
  }

  DevelFatal(*this).in("DeviceEntity::getNumericalMatrixSensitivity") 
      << "Function not implemented yet ";
  return false;

  return getNumericalMatrixSensitivity(defaultParamName_,
    d_dfdx_dp, d_dqdx_dp, F_lids, Q_lids, F_jacLIDs, Q_jacLIDs);
}
//
//-----------------------------------------------------------------------------
// Function      : DeviceEntity::getAnalyticACSensitivityAvailable 
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/2/2019
//-----------------------------------------------------------------------------
bool DeviceEntity::getAnalyticACSensitivityAvailable (const std::string & paramName)
{
  ParameterMap::const_iterator p_i = getParameterMap().find(paramName);
  if (p_i == getParameterMap().end())
  {
    DevelFatal(*this).in("DeviceEntity::getAnalyticACSensitivityAvailable") << "Unrecognized parameter " << paramName;
    return false;
  }

  const Descriptor &param = *(*p_i).second;
  return param.getAnalyticACSensitivityAvailable ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::getAnalyticACSensitivityAvailableDefaultParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/2/2019
//-----------------------------------------------------------------------------
bool DeviceEntity::getAnalyticACSensitivityAvailableDefaultParam ()
{
  if (defaultParamName_.empty())
  {
    DevelFatal(*this).in("DeviceEntity::getAnalyticACSensitivityAvailableDefaultParam") 
      << "Device does not have a default parameter";
    return false;
  }

  return getAnalyticACSensitivityAvailable(defaultParamName_);
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::getAnalyticBSensVectorsforAC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/2/2019
//-----------------------------------------------------------------------------
bool DeviceEntity::getAnalyticBSensVectorsforAC (const std::string & paramName,
          std::vector< std::complex<double> > & dbdp,
          std::vector<int> &        BindicesVec)
{
  ParameterMap::const_iterator p_i = getParameterMap().find(paramName);
  if (p_i == getParameterMap().end())
  {
    DevelFatal(*this).in("DeviceEntity::analyticSensitivityAvailable") << "Unrecognized parameter " << paramName;
    return false;
  }

  const Descriptor &param = *(*p_i).second;
  return param.getAnalyticBSensVectorsforAC (*this, paramName, dbdp, BindicesVec);
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::getAnalyticBSensVectorsforACDefaultParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/2/2019
//-----------------------------------------------------------------------------
bool DeviceEntity::getAnalyticBSensVectorsforACDefaultParam (
          std::vector< std::complex<double> > & dbdp,
          std::vector<int> &        BindicesVec)
{
  if (defaultParamName_.empty())
  {
    DevelFatal(*this).in("DeviceEntity::getAnalyticBSensVectorsforACDefaultParam") 
      << "Device does not have a default parameter";
    return false;
  }

  return getAnalyticBSensVectorsforAC(defaultParamName_, dbdp, BindicesVec);
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::getNumericalBSensVectorsforACDefaultParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/2/2019
//-----------------------------------------------------------------------------
bool DeviceEntity::getNumericalBSensVectorsforACDefaultParam (
          std::vector< std::complex<double> > & dbdp,
          std::vector<int> &        BindicesVec)
{
  return true;
}

//

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::setParam
//
// Purpose       : This function loops over the vector of parameters, and
//                 sets the specified one (if found) to a specified value.
//
// Special Notes : This is kind of tricky, b/c some parameters are actually
//                 deep inside other classes (like a source class, for
//                 example)
//
//                 This function always returns a true, b/c there are many
//                 instances, (for example running in parallel), where one
//                 could set a param that didn't exist locally on
//                 processor.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/25/03
//-----------------------------------------------------------------------------
bool DeviceEntity::setParam(const std::string & paramName, double val, bool overrideOriginal)
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "DeviceEntity::setParam  with paramname = " << paramName << " value = " << val << " overrideOriginal = " << overrideOriginal << std::endl;
  }
  
  ParameterMap::const_iterator p_i = getParameterMap().find(paramName);
  if (p_i == getParameterMap().end())
    return false;


  const Descriptor &param = *(*p_i).second;

  if (isTempParam(paramName) && param.getAutoConvertTemperature())
    val += CONSTCtoK;

  if (param.isType<double>())
    param.value<double>(*this) = val;
  else if (param.isType<int>())
    param.value<int>(*this) = static_cast <int> (val);
  else if (param.isType<long>())
    param.value<long>(*this) =  static_cast <long> (val);
  else if (param.isType<bool>())
    param.value<bool>(*this) = (val != 0);
  else
    DevelFatal(*this) << "Illegal type for parameter " << paramName;

  if (param.hasGivenMember())
    param.setGiven(*this, true);

  Xyce::Device::setValueGiven(*this, param.getSerialNumber(), true);

  if (param.hasOriginalValueStored() && overrideOriginal)
  {
    //if (DEBUG_DEVICE)
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << " Overriding original value for parameter " << paramName << "  with value " << val << std::endl;
    }
    Xyce::Device::setOriginalValue(*this,param.getSerialNumber(), val);
  }
  
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::findParam
//
// Purpose       : returns true if parameter exists for device
//
// Scope         : public
// Creator       : David G. Baur, Raytheon
// Creation Date : 07/01/2016
//-----------------------------------------------------------------------------
bool DeviceEntity::findParam(const std::string & param_name) const
{
  return getParameterMap().find(param_name) != getParameterMap().end();
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::getParam
//
// Purpose       : returns the value of the requested param.
//
// Special Notes : This  function currently assumes that the requested
//                 param is a double-precision number.
//
//                 Parameters are not case-dependent.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/25/03
//-----------------------------------------------------------------------------
bool DeviceEntity::getParam(const std::string & name, double & result) const
{
  double val = 0.0;
  bool found = false;

  ParameterMap::const_iterator p_i = getParameterMap().find(name);
  if (p_i != getParameterMap().end())
  {
    found = true;
    const Descriptor &param = *(*p_i).second;
    if (param.isType<double>())
      val = param.value<double>(*this);
    else if (param.isType<int>())
      val = static_cast <double> (param.value<int>(*this));
    else if (param.isType<long>())
      val = static_cast <double> (param.value<long>(*this));
    else if (param.isType<bool>())
    {
      if (param.value<bool>(*this))
        val = 1;
      else
        val = 0;
    }
    else
    {
      DevelFatal(*this).in("DeviceEntity::getParam") << "Illegal type for parameter " << name;
    }
    if (isTempParam(name)  && param.getAutoConvertTemperature())
      val -= CONSTCtoK;
  }
  else
  {
    // If not recognized, just do nothing
  }
  result = val;

  return found;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::setDefaultParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/06/03
//-----------------------------------------------------------------------------
bool DeviceEntity::setDefaultParam (double val, bool overrideOriginal)
{
  if (defaultParamName_.empty())
  {
    DevelFatal(*this).in("DeviceEntity::setDefaultParam") << "Device does not have a default parameter";
  }

  return setParam(defaultParamName_, val, overrideOriginal);
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::getDefaultParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/06/03
//-----------------------------------------------------------------------------
double DeviceEntity::getDefaultParam() const
{
  if (defaultParamName_.empty())
  {
    return 0.0;
  }

  double result = 0.0;

  getParam(defaultParamName_, result);

  return result;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::setDependentParameter
// Purpose       : Add expression, param pairs for future updates
// Special Notes :  This is an overloaded method, used instead of the old
//                  monolithic one.
// Scope         : protected
// Creator       : Tom Russo
// Creation Date : 6 Nov 07
//-----------------------------------------------------------------------------
///
/// Create a Depend object for a parameter, fill it in, and add to the dependParam_ vector
///
/// @param[in] par Util::Param object with name/value pair
/// @param[in] res A pointer to the member variable of the entity into which parameter values are stored
/// @param[in] depend expression access permissions for the parameter
///
/// This function is used when it is determined that a single device
/// parameter (as opposed to a parameter inside a vector composite)
/// depends (via an expression) on external quantities.
///
/// We create a Depend object and add it to the list, then call the
/// utility version of setDependentParameter to populate it
///
/// This version of setDependentParameter was split out of an older
/// "monolithic" method that was festooned with if statements to
/// perform multiple functions.  This function is now single purpose.
/// What remains of the old monolithic function is now in the
/// overloaded utility function below, also named
/// "setDependentParameter," but with different prototype.
///
/// @note The Depend objects contain an expression pointer, and are
/// copied using shallow copy.  Thus, when the local variable
/// "dependentParam" is added to the dependentParams_ vector, the
/// expression pointer is copied, not the underlying expression
/// object.  This fact is used after the push to the back of
/// dependentParams_, when the expression pointer in dependentParam is
/// used to access the expression for evaluation.  Because the pointer
/// is shared between dependentParam and the copy that has been pushed
/// back, this works.
///
/// @author Tom Russo
/// @date 11/6/2007
///
double DeviceEntity::setDependentParameter (Util::Param & par,
                                            double *res,
                                            ParameterType::ExprAccess depend)

{
  Depend dependentParam;
  setDependentParameter(par, dependentParam, depend);

  dependentParam.resultU.result = res;
  dependentParam.vectorIndex = -1;
  dependentParams_.push_back(dependentParam);

  double rval;
  dependentParam.expr->evaluateFunction (rval);

  return rval;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::setDependentParameter
// Purpose       : Add expression, param pairs for future updates
// Special Notes :  This is an overloaded method, used instead of the old
//                  monolithic one, and is specifically to set an element
//                  of a double vector.
// Scope         : protected
// Creator       : Tom Russo
// Creation Date : 6 Nov 07
//-----------------------------------------------------------------------------
///
/// Create a Depend object for a parameter, fill it in, and add to the dependParam_ vector
///
/// @param[in] par Util::Param object with name/value pair
/// @param[in] res A pointer to a vector of doubles, which should be a member variable of the entity, and in which the values of the parameter should be stored
/// @param[in] ind  integer index into the res vector, denoting the location in the vector where this parameter should be stored
/// @param[in] depend expression access permissions for the parameter
///
/// This function is used when it is determined that a parameter
/// inside a vector composite depends (via an expression) on external
/// quantities.
///
/// We create a Depend object and add it to the list, then call the
/// utility version of setDependentParameter to populate it
///
/// This version of setDependentParameter was split out of an older
/// "monolithic" method that was festooned with if statements to
/// perform multiple functions.  This function is now single purpose.
/// What remains of the old monolithic function is now in the
/// overloaded utility function below, also named
/// "setDependentParameter," but with different prototype.
///
/// @note The Depend objects contain an expression pointer, and are
/// copied using shallow copy.  Thus, when the local variable
/// "dependentParam" is added to the dependentParams_ vector, the
/// expression pointer is copied, not the underlying expression
/// object.  This fact is used after the push to the back of
/// dependentParams_, when the expression pointer in dependentParam is
/// used to access the expression for evaluation.  Because the pointer
/// is shared between dependentParam and the copy that has been pushed
/// back, this works.
///
/// @author Tom Russo
/// @date 11/6/2007
///
double DeviceEntity::setDependentParameter (Util::Param & par,
                                            std::vector<double> *res,
                                            int ind,
                                            ParameterType::ExprAccess depend)

{
  Depend dependentParam;
  setDependentParameter(par,dependentParam, depend);

  dependentParam.resultU.resVec = res;
  dependentParam.vectorIndex = ind;
  dependentParams_.push_back(dependentParam);

  double rval;
  dependentParam.expr->evaluateFunction (rval);

  return rval;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::setDependentParameter
// Purpose       : Add expression, param pairs for future updates
// Special Notes :  This is a utility version, used by overloaded methods.
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/18/04
//-----------------------------------------------------------------------------
///
/// Fill in a Depend structure that establishes that a given parameter
/// is dependent on external quantities (solution variables, global parameters,
/// etc.), associates an expression with the parameter, and sets up necessary
/// information to evaluate the parameter as needed.
///
/// The DeviceEntity maintains a list of Depend structures in vector
/// dependParams_.  During simulation, this list is consulted whenever a
/// device load needs to happen, and it is used to update parameters whose
/// dependent quantities have changed.
///
/// @param[in] par  The Util::Param object with the name/value pair
/// @param dependentParam A pre-constructed Depend object to be filled in
/// @param[in] depend expression access permissions for the parameter
///
/// This function is called as a utility function by the two other
/// overloaded setDependentParameter methods.
///
/// @author Dave Shirley, PSSI
/// @date 11/18/04
///
/// ERK. 4/24/2020.  A lot of this function isn't needed anymore with the 
/// new expression library.
///
void DeviceEntity::setDependentParameter (Util::Param & par,
                                          Depend & dependentParam,
                                          ParameterType::ExprAccess depend)
{
  std::vector<std::string> instances, leads, names, variables;

  dependentParam.name = par.tag();
  ParameterMap::const_iterator p_i = getParameterMap().find(dependentParam.name);
  const Descriptor &param = *(*p_i).second;
  if (isTempParam(par.tag()) && param.getAutoConvertTemperature())
  {
    //User input for temperature is in degrees C.  Devices use degrees K internally.
    dependentParam.expr = new Util::Expression (
        solState_.expressionGroup_, "(" + par.stringValue() + ")+CONSTCtoK"); 

    // Need to get the string names of the variables that are in par's expression,
    // and then use that vector of string names to make them variables in the
    // dependentParam's expression also.  This is needed to support global parameters
    // in expressions for the TEMP instance paramter.
    std::vector<std::string> parVars;
    par.getValue<Util::Expression>().getVariables(parVars); 
    for (std::vector<std::string>::const_iterator pv_i=parVars.begin(); pv_i != parVars.end(); ++pv_i)
    {
      dependentParam.expr->make_var(*pv_i);
    }

    dependentParam.expr->make_constant (std::string("CONSTCTOK"), CONSTCtoK);
  }
  else
  {
    // ERK.  Why do we need to copy construct this?
    dependentParam.expr = new Util::Expression (par.getValue<Util::Expression>());
  }

  if (param.hasOriginalValueStored())
  {
    dependentParam.storeOriginal=true;
    dependentParam.serialNumber=param.getSerialNumber();
  }
  else
  {
    dependentParam.storeOriginal=false;
  }

  names.clear();
  leads.clear();
  instances.clear();
  variables.clear();

  dependentParam.expr->getVoltageNodes(names);
  dependentParam.expr->getLeadCurrents(leads);
  dependentParam.expr->getDeviceCurrents(instances);
  dependentParam.expr->getVariables(variables); 

  //std::vector<std::string>::iterator s;
  std::vector<std::string>::iterator iterS;

  if (!(depend & ParameterType::SOLN_DEP))
  {
    if (names.size() > 0 || instances.size() > 0)
    {
      UserError(*this) << "Parameter " << par.tag() << " is not allowed to depend on voltage/current values";
      return;
    }
    if (depend & ParameterType::NO_DEP)
    {
      if (dependentParam.expr->get_num(XEXP_SPECIAL) > 0)
      {
        UserError(*this) << "Parameter " << par.tag() << " is not allowed to depend on time";
        return;
      }
    }
  }

  if (leads.size() > 0)
  {
    char type;
    int index;
    for (std::vector<std::string>::const_iterator n_i=leads.begin(); n_i != leads.end(); ++n_i)
    {
      index = n_i->find_last_of(":");
      if (index == std::string::npos )
        type = (*n_i)[0];
      else
        type = (*n_i)[index+1];

      if (type != 'B' && type != 'E' && type != 'H')
      {
        UserError(*this) << "Illegal use of lead current specification in expression '" << dependentParam.expr->get_expression()
                         << "' in parameter " << par.tag();
      }
    }
    names.insert( names.end(), leads.begin(), leads.end() );
  }

  names.insert( names.end(), instances.begin(), instances.end() );

  dependentParam.lo_var = expVarNames.size();
  dependentParam.n_vars = names.size();
  dependentParam.vals.resize(dependentParam.n_vars);
  int expVarLen = dependentParam.lo_var+dependentParam.n_vars;
  expVarGIDs.resize(expVarLen);
  expVarLIDs.resize(expVarLen);
  expVarVals.resize(expVarLen);

  if (!variables.empty())
    names.insert( names.end(), variables.begin(), variables.end() );

  if ( !names.empty() )
  {
    // Order the names in the expression so that it agrees with the order
    // in names.
    dependentParam.expr->order_names( names );
  }

  for (int i=0 ; i<dependentParam.n_vars ; ++i)
    expVarNames.push_back(names[i]);

#if 0
  // ERK.  FIX THIS!   commenting out so this will compile
  // Is this needed?
  if (dependentParam.n_vars > 0)
  {
    std::vector<double> zeros;
    zeros.resize(dependentParam.n_vars);
    for (int i=0 ; i<dependentParam.n_vars ; ++i)
      zeros[i] = 0;

    dependentParam.expr->set_vars(zeros);
  }
#endif

  dependentParam.global_params.clear();
  if (!variables.empty())
  {
    for (iterS=variables.begin() ; iterS!=variables.end() ; ++iterS)
    {
      GlobalParameterMap::iterator global_param_it = globals_.global_params.find(*iterS);
      if (global_param_it == globals_.global_params.end())
      {
        UserError(*this) << "Global parameter " << *iterS << " not found";
      }
      else 
      {
#if 0
    // ERK.  FIX THIS!   commenting out so this will compile
        dependentParam.expr->set_var(*iterS, (*global_param_it).second);
#endif
        dependentParam.global_params.push_back(*iterS);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::updateDependentParameters
// Purpose       : Update values of parameters defined as expressions
// Special Notes :
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/15/05
//-----------------------------------------------------------------------------
bool DeviceEntity::updateDependentParameters(const Linear::Vector & vars, bool changed)
{
  std::vector<Depend>::iterator dpIter = dependentParams_.begin();
  std::vector<Depend>::iterator end = dependentParams_.end();
  double rval(0.0);

  for ( ; dpIter != end ; ++dpIter)
  {
    // ERK.  4/24/2020. This (the changed bool) was conditionally true/false 
    // depending on various calls such as set_sim_time, which let each expression 
    // report back if its internal time variable (or temp, or freq, etc) had changed.
    //
    // Those function calls don't exist anymore due to the new expression refactor.
    //
    // This logic should maybe be updated to use the same logic that devices do w.r.t things like
    // limiting.  If on a new time step, time changes, otherwise not.  Etc.  
    // If this isn't changed, then a lot of "processParam" function calls are going 
    // to be called unneccessarily.
    //
    // But that will have to come later.

    if ( !(dpIter->expr->getIsConstant()) ) // ERK.  5/3/2020.  This works, but refine later.
    {
      changed = true;
    }

    if (changed)
    {
      dpIter->expr->evaluateFunction (rval);
      if (dpIter->vectorIndex==-1)
        *(dpIter->resultU.result) = rval;
      else
        (*(dpIter->resultU.resVec))[dpIter->vectorIndex] = rval;

      if (dpIter->storeOriginal)
        Xyce::Device::setOriginalValue(*this,dpIter->serialNumber,rval);
    }
  }

  return changed;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::updateGlobalParameters
// Purpose       : Update values of global parameters in expressions
// Special Notes :
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/17/05
//-----------------------------------------------------------------------------
bool DeviceEntity::updateGlobalParameters(GlobalParameterMap & global_map)
{
  std::vector<Depend>::iterator dpIter = dependentParams_.begin();
  std::vector<Depend>::iterator end = dependentParams_.end();
  double rval;
  int i, hi;
  bool changed = false;

  for ( ; dpIter != end ; ++dpIter)
  {
    if (!dpIter->global_params.empty())
    {
      std::vector<std::string>::iterator gp=dpIter->global_params.begin();
      std::vector<std::string>::iterator gend=dpIter->global_params.end();
      for ( ; gp != gend; ++gp)
      {
        if (global_map.find(*gp) == global_map.end())
        {
          DevelFatal(*this).in("DeviceEntity::updateGlobalParameters") << "Failed to find global parameter " << *gp;
        }
#if 0
    // ERK.  FIX THIS!   commenting out so this will compile
        if (dpIter->expr->set_var(*gp, global_map[*gp]))
          changed = true;
#endif
      }
    }
  }

  return changed;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::updateDependentParameters
// Purpose       : Update values of parameters defined as expressions
// Special Notes :
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/18/04
//-----------------------------------------------------------------------------
bool DeviceEntity::updateDependentParameters()
{
  double rval;
  bool changed = false;

  std::vector<Depend>::iterator dpIter = dependentParams_.begin();
  std::vector<Depend>::iterator end = dependentParams_.end();
  for ( ; dpIter != end; ++dpIter)
  {
    // ERK.  4/24/2020. This (the changed bool) was conditionally true/false 
    // depending on various calls such as set_sim_time, which let each expression 
    // report back if its internal time variable (or temp, or freq, etc) had changed.
    //
    // Those function calls don't exist anymore due to the new expression refactor.
    //
    // This logic should maybe be updated to use the same logic that devices do w.r.t things like
    // limiting.  If on a new time step, time changes, otherwise not.  Etc.  
    // If this isn't changed, then a lot of "processParam" function calls are going 
    // to be called unneccessarily.
    //
    // But that will have to come later.
    changed = true;

    dpIter->expr->evaluateFunction (rval);
    if (dpIter->vectorIndex == -1)
      *(dpIter->resultU.result) = rval;
    else
      (*(dpIter->resultU.resVec))[dpIter->vectorIndex] = rval;
  }

  return changed;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::updateDependentParameters
// Purpose       : Update values of parameters defined as expressions with a
//                 specified temperature
// Special Notes :
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/11/06
//-----------------------------------------------------------------------------
bool DeviceEntity::updateDependentParameters(double tempIn)
{
  double rval;
  bool changed = false;

  std::vector<Depend>::iterator dpIter = dependentParams_.begin();
  std::vector<Depend>::iterator end = dependentParams_.end();
  for ( ; dpIter != end; ++dpIter)
  {
    // ERK.  4/24/2020. This (the changed bool) was conditionally true/false 
    // depending on various calls such as set_sim_time, which let each expression 
    // report back if its internal time variable (or temp, or freq, etc) had changed.
    //
    // Those function calls don't exist anymore due to the new expression refactor.
    //
    // This logic should maybe be updated to use the same logic that devices do w.r.t things like
    // limiting.  If on a new time step, time changes, otherwise not.  Etc.  
    // If this isn't changed, then a lot of "processParam" function calls are going 
    // to be called unneccessarily.
    //
    // But that will have to come later.
    changed = true;

    dpIter->expr->evaluateFunction (rval);
    if (dpIter->vectorIndex == -1)
      *(dpIter->resultU.result) = rval;
    else
      (*(dpIter->resultU.resVec))[dpIter->vectorIndex] = rval;

  }

  return changed;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::getParamBreakpoints
// Purpose       : Add breakpoints caused by discontinuities in computed params
// Special Notes :
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/18/04
//-----------------------------------------------------------------------------
bool DeviceEntity::getParamBreakpoints( std::vector<Util::BreakPoint> & breakPointTimes )
{
  double bTime;

  std::vector<Depend>::iterator dpIter = dependentParams_.begin();
  std::vector<Depend>::iterator end = dependentParams_.end();
  for ( ; dpIter != end; ++dpIter)
  {
    bTime = dpIter->expr->getBreakPoints(breakPointTimes);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::given
// Purpose       : Return whether param was given
// Special Notes :
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/23/04
//-----------------------------------------------------------------------------
bool DeviceEntity::given( const std::string & parameter_name ) const
{
  ParameterMap::const_iterator it = getParameterMap().find(parameter_name);
  if (it == getParameterMap().end())
    DevelFatal(*this).in("DeviceEntity::given") << "Unrecognized parameter " << parameter_name;

  return Xyce::Device::wasValueGiven(*this, (*it).second->getSerialNumber());
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::setParams
// Purpose       : Set parameters according to a vector of params.  Used to
//                 set instance or model parameter sets to netlist values
//
//                 This function also takes care of setting up necessary
//                 infrastructure when it is detected that a parameter
//                 has a value that is an expression dependent on time or
//                 other variables.  These will be added to a dependent
//                 variables vector for recomputation when the things they
//                 depend on change.
// Special Notes : Converted to member function from free function setParameters
//                 by Tom Russo.  The history of this function is somewhat
//                 bizarre, like its implementation, but at one point it was
//                 a protected member function, then was converted to a
//                 free function and renamed by Dave Baur, but only ever
//                 called by a wrapper member function called setParams.
//                 Since it was unclear why it was ever necessary for this
//                 two-layered calling sequence, I've converted it back
//                 to a single member function on 9/6/2016.
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/20/04
//-----------------------------------------------------------------------------
void DeviceEntity::setParams(const std::vector<Param> &params)
{
  std::vector<Param>::const_iterator begin=params.begin();
  std::vector<Param>::const_iterator end=params.end();
  
  std::vector<std::string> composite_name_list;
  unordered_map<std::string, std::vector<CompositeParam *>, HashNoCase, EqualNoCase> composite_parameter_map;
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl << "In DeviceEntity::setParams, for ";
    printName(Xyce::dout());
    Xyce::dout() << " parameters are:" << std::endl;
    for (std::vector<Param>::const_iterator it = begin; it != end; ++it)
    {
      const Param &param = *it;
      Xyce::dout() << "Param = " << param.tag() << ", Type = ";
      int tmpType = param.getType();
      switch (tmpType)
      {
        case Util::STR:
          Xyce::dout() << "STR";
          break;
        case Util::DBLE:
          Xyce::dout() << "DBLE";
          break;
        case Util::INT:
          Xyce::dout() << "INT";
          break;
        case Util::LNG:
          Xyce::dout() << "LNG";
          break;
        case Util::EXPR:
          Xyce::dout() << "EXPR";
          break;
        case Util::BOOL:
          Xyce::dout() << "BOOL";
          break;
        case Util::STR_VEC:
          Xyce::dout() << "STR_VEC";
          break;
        case Util::INT_VEC:
          Xyce::dout() << "INT_VEC";
          break;
        case Util::DBLE_VEC:
          Xyce::dout() << "DBLE_VEC";
          break;
        case Util::DBLE_VEC_IND:
          Xyce::dout() << "DBLE_VEC_IND";
          break;
        case Util::COMPOSITE:
          Xyce::dout() << "COMPOSITE";
          break;
        default:
          Xyce::dout() << "Unknown";
      }
      Xyce::dout() << ", Value = " << param.stringValue();

      if (param.given())
      {
        Xyce::dout() << "  given=TRUE";
      }
      else
      {
        Xyce::dout() << "  given=FALSE";
      }

      if (param.default_val())
      {
        Xyce::dout() << "  default=TRUE" << std::endl;
      }
      else
      {
        Xyce::dout() << "  default=FALSE" << std::endl;
      }
    }
    Xyce::dout() << std::endl;
  }

  for (std::vector<Param>::const_iterator param_it = begin; param_it != end; ++param_it)
  {
    Param &param = const_cast<Param &>(*param_it);

    const std::string &tag = param.tag();

    // Is this parameter in the Entity?
    ParameterMap::const_iterator entity_parameter_it = getParameterMap().find(tag);
    if (entity_parameter_it != getParameterMap().end())
    {
      const Descriptor &descriptor = *(*entity_parameter_it).second;
      if (descriptor.hasGivenMember())
      {
        if (param.given())
        {
          descriptor.setGiven(*this, true);
        }
        else if (descriptor.getGiven(*this))
        {
          continue;
        }
      }

      Xyce::Device::setValueGiven(*this, descriptor.getSerialNumber(), param.given());
      if (param.given() || param.default_val())
      {
        if ( param.getType() == Util::EXPR )
        {
          if (descriptor.isType<double>())
          {
            double val = setDependentParameter(param, &(descriptor.value<double>(*this)), descriptor.getExpressionAccess());
            param.setVal(val);
            if (descriptor.hasOriginalValueStored())
              Xyce::Device::setOriginalValue(*this,descriptor.getSerialNumber(), val);
          }
          else if (descriptor.isType<std::vector<double> >())
          {
            int ind = (descriptor.value<std::vector<double> >(*this)).size();
            double val = setDependentParameter (param, &(descriptor.value<std::vector<double> >(*this)), ind, descriptor.getExpressionAccess());
            (descriptor.value<std::vector<double> >(*this)).push_back(val);
          }
          else
          {
            DevelFatal(*this).in("DeviceEntity::setParams") << "Non double param " <<  tag << " cannot be set to expression";
          }
        }
        else
        {
          if (descriptor.isType<double>())
          {
            if (!param.isNumeric())
            {
              UserError(*this) << "Cannot convert parameter " << tag <<  " to a numeric value from " << param.stringValue();
              continue;
            }

            descriptor.value<double>(*this) = param.getImmutableValue<double>();
            if (isTempParam(tag) && descriptor.getAutoConvertTemperature())
            {
              descriptor.value<double>(*this) += CONSTCtoK;
            }
            if (descriptor.hasOriginalValueStored())
            {
              Xyce::Device::setOriginalValue(*this, descriptor.getSerialNumber(), descriptor.value<double>(*this));
            }
          }
          else if (descriptor.isType<std::string>())
          {
            descriptor.value<std::string>(*this) = param.stringValue();
          }
          else if (descriptor.isType<int>())
          {
            if (!param.isInteger())
            {
              UserError(*this) << "Cannot convert parameter " << tag << " to an integer value from " << param.stringValue();
              continue;
            }
            descriptor.value<int>(*this) = param.getImmutableValue<int>();
            if (descriptor.hasOriginalValueStored())
            {
              Xyce::Device::setOriginalValue(*this, descriptor.getSerialNumber(), static_cast<double> (descriptor.value<int>(*this)));
            }
          }
          else if (descriptor.isType<long>())
          {
            if (!param.isInteger())
            {
              UserError(*this) << "Cannot convert parameter " << tag << " to an integer value from " << param.stringValue();
              continue;
            }
            descriptor.value<long>(*this) = param.getImmutableValue<long>();
            if (descriptor.hasOriginalValueStored())
            {
              Xyce::Device::setOriginalValue(*this, descriptor.getSerialNumber(), static_cast<double> (descriptor.value<long>(*this)));
            }
          }
          else if (descriptor.isType<bool>())
          {
            if (!param.isBool())
            {
              UserError(*this) << "Cannot convert parameter " << tag << " to a logical value from " << param.stringValue();
              continue;
            }
            descriptor.value<bool>(*this) = param.getImmutableValue<bool>();
            if (descriptor.hasOriginalValueStored())
            {
              if (descriptor.value<bool>(*this))
              {
                Xyce::Device::setOriginalValue(*this, descriptor.getSerialNumber(), 1.0);
              }
              else
              {
                Xyce::Device::setOriginalValue(*this, descriptor.getSerialNumber(), 0.0);
              }
            }
          }
          else if (descriptor.isType<std::vector<int> >())
          {
            if (param.getType() == Util::INT_VEC)
            {
              (descriptor.value<std::vector<int> >(*this)) = param.getValue<std::vector<int> >();
            }
            else if (param.getType() == Util::INT)
            {
              (descriptor.value<std::vector<int> >(*this)).push_back(param.getImmutableValue<int>());
            }
          }
          else if (descriptor.isType<std::vector<double> >())
          {
            if (param.getType() == Util::DBLE_VEC)
            {
              (descriptor.value<std::vector<double> >(*this)) = param.getValue<std::vector<double> >();
            }
            else if (param.getType() == Util::DBLE)
            {
              (descriptor.value<std::vector<double> >(*this)).push_back(param.getImmutableValue<double>());
            }
          }
          else if (descriptor.isType<std::vector<std::string> >())
          {
            if (param.getType() == Util::STR_VEC)
            {
              (descriptor.value<std::vector<std::string> >(*this)) = param.getValue<std::vector<std::string> >();
            }
            else if (param.getType() == Util::STR)
            {
              (descriptor.value<std::vector<std::string> >(*this)).push_back(param.stringValue());
            }
          }
          else if (descriptor.hasCompositeData())
          {
            composite_parameter_map[tag].clear();
            composite_name_list.push_back(tag);

            // Note: ERK.  This push-back is done to process the base-param tag of a vector composite.
            // For example, if the composite parameters are things like REGION0.XWIDTH, where REGION
            // is the base parameter tag, 0 is the index, and XWIDTH is the subcomponent, the
            // tag that should be pushed back is tagES=REGION.
            //
            // Note: ERK:  This function seems to implicitly rely on the base parameter always preceeding
            // subcomponent parameters.  So REGION (alone) should preceed REGION0.XWIDTH in  the
            // STL vector params that is passed into this function.  If it doesn't then vc_stat will
            // stay false and a fatal error will get thrown.
            if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
            {
              Xyce::dout() << "pushing back composite " << tag << std::endl;
            }
          }
          else
          {
            DevelFatal(*this).in("DeviceEntity::setParams") << "Unknown type";
          }
        }
      }
    }

    // Is it a vector (why do nothing?)
    else if (param.stringValue() == "VECTOR")
    {
    }

    // Must be a composite
    else
    {
      bool vc_stat = false;

      std::string::size_type dot = tag.find_first_of('.');
      if (dot != std::string::npos)
      {
        for (std::vector<std::string>::const_iterator it = composite_name_list.begin(); it != composite_name_list.end(); ++it)
        {
          const std::string &composite_name = *it;

          if (tag.find(composite_name) == 0) // Tag starts with the composite parameter name
          {
            std::string param_name(tag.begin() + dot + 1, tag.end());

            vc_stat = true;

            int n = 0;
            {
              std::istringstream is(std::string(tag.begin() + composite_name.size(), tag.begin() + dot));
              is >> n;
            }

            if (param_name == "NAME")
            {
              if (n != composite_parameter_map[composite_name].size())
              {
                DevelFatal(*this).in("DeviceEntity::setParams") << "Error filling 'NAME' vector param " <<  composite_name;
              }

              std::string name = param.stringValue();
#if 0
              // ERK
              std::cout << "Must be a composite: composite_name = " << composite_name << " name " << name << std::endl;

              if ( param.getType() == Util::EXPR ) 
              { 
                std::cout << "param " << tag << " is Util::EXPR type" << std::endl; 
                Util::Expression & expToBeDumped = param.getValue<Util::Expression>();
                expToBeDumped.dumpParseTree();
              }
              else { std::cout << "param " << tag << " is NOT Util::EXPR type" << std::endl; }
#endif
              CompositeParam *composite = constructComposite(composite_name, name);
              composite_parameter_map[composite_name].push_back(composite);
              setDefaultParameters(*composite, composite->getParameterMap().begin(), composite->getParameterMap().end(), devOptions_);
#if 0
              // ERK
              ParameterMap::const_iterator pIter = composite->getParameterMap().begin();
              for (;pIter!= composite->getParameterMap().end(); pIter++)
              {
                const std::string &name = (*pIter).first;
                const Descriptor &param = *(*pIter).second;
                //std::cout << "name = " << name << "paramDescriptor = " << param << std::endl;
              }

              //std::cout << "
#endif
            }
            else
            {
#if 0
              // ERK
              std::cout << "Must be a composite: composite_name = " << composite_name << " name " << param.stringValue() << std::endl;

              if ( param.getType() == Util::EXPR ) 
              { 
                std::cout << "param " << tag << " is Util::EXPR type" << std::endl; 
                Util::Expression & expToBeDumped = param.getValue<Util::Expression>();
                expToBeDumped.dumpParseTree();
              }
              else { std::cout << "param " << tag << " is NOT Util::EXPR type" << std::endl; }
#endif

              if (n >= composite_parameter_map[composite_name].size())
              {
                UserFatal(*this) << "Error in definition of vector-composite parameter. "
                                 << "NAME parameter must be first parameter given in definition of " << composite_name;
              }
            }
            setParamFromVCParam(*composite_parameter_map[composite_name][n], composite_name, param_name, param);
          }
        }
      }
      if (!vc_stat)
      {
        UserFatal(*this) << "Unknown error while parsing possible vector-composite parameter";
      }
    }
  }

  if (!composite_name_list.empty())
  {
    for (std::vector<std::string>::const_iterator name_it = composite_name_list.begin(); name_it != composite_name_list.end(); ++name_it)
    {
      const std::string &vcs = *name_it;

      for (std::vector<CompositeParam *>::iterator it =  composite_parameter_map[vcs].begin(); it != composite_parameter_map[vcs].end(); ++it)
      {
        (*it)->processParams();
      }
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : DeviceEntity::setParamFromVCParam
// Purpose       : Set values for a parameter contained within a 
//               : vector-composite parameter. 
// Special Notes : This was changed from a "free function" to a member function
//               : by Pete Sholander in December 2017.  That allowed the error 
//               : message to be improved (e.g., provide info on the device name).
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/06/05
//-----------------------------------------------------------------------------
void DeviceEntity::setParamFromVCParam(CompositeParam &composite_param, 
				       const std::string & composite_name,
                                       const std::string & pName, 
                                       const Param & ndParam)
{
  ParameterMap::const_iterator p_i = composite_param.getParameterMap().find(pName);
  if (p_i != composite_param.getParameterMap().end())
  {
    const Descriptor &p = *(*p_i).second;
    if (p.hasGivenMember())
    {
      if (ndParam.given())
        p.setGiven(composite_param, true);
      else if (p.getGiven(composite_param))
        return;
    }
    Xyce::Device::setValueGiven(composite_param, p.getSerialNumber(), ndParam.given());
    if (ndParam.given() || ndParam.default_val())
    {
      if (p.isType<double>())
      {
        p.value<double>(composite_param) = ndParam.getImmutableValue<double>();
        if (isTempParam(pName) && p.getAutoConvertTemperature())
          p.value<double>(composite_param) += CONSTCtoK;
        if (p.hasOriginalValueStored())
          Xyce::Device::setOriginalValue(composite_param, p.getSerialNumber(), p.value<double>(composite_param));
      }
      else if (p.isType<std::string>())
      {
        p.value<std::string>(composite_param) = ndParam.stringValue();
      }
      else if (p.isType<int>())
      {
        p.value<int>(composite_param) = ndParam.getImmutableValue<int>();
        if (p.hasOriginalValueStored())
          Xyce::Device::setOriginalValue(composite_param, p.getSerialNumber(), static_cast<double>(p.value<int>(composite_param)));
      }
      else if (p.isType<long>())
      {
        p.value<long>(composite_param) = ndParam.getImmutableValue<long>();
        if (p.hasOriginalValueStored())
          Xyce::Device::setOriginalValue(composite_param, p.getSerialNumber(), static_cast<double>(p.value<long>(composite_param)));
      }
      else if (p.isType<bool>())
      {
        p.value<bool>(composite_param) = (ndParam.getImmutableValue<double>() != 0.0);
        if (p.hasOriginalValueStored())
        {
          if (p.value<bool>(composite_param))
            Xyce::Device::setOriginalValue(composite_param, p.getSerialNumber(), 1.0);
          else
            Xyce::Device::setOriginalValue(composite_param, p.getSerialNumber(), 0.0);
        }
      }
      else if (p.isType<std::vector<double> >())
      {
       (p.value<std::vector<double> >(composite_param)).push_back(ndParam.getImmutableValue<double>());
      }
      else if (p.isType<std::vector<std::string> >())
      {
        p.value<std::vector<std::string> >(composite_param).push_back(ndParam.stringValue());
      }
      else
      {
        Report::DevelFatal().in("DeviceEntity::setParamFromVCParam") << "Unknown parameter type for " << pName;
      }
    }
  }
  else
  {
    UserFatal(*this) << "Either vector-composite parameter " << composite_name 
                     << " is improperly formatted or parameter " <<  pName  
                     << " is not defined for this device or model type";
  }
}


// this function is called from N_IO_CircuitMetadata.C in the function CircuitMetadata::getDeviceMetadata.
void populateParams(const ParameterMap &parameter_map, std::vector<Param> &param_list, CompositeParamMap &composite_param_map)
{
  for (ParameterMap::const_iterator it = parameter_map.begin(); it != parameter_map.end(); ++it)
  {
    const Descriptor &param = *(*it).second;

    if (param.isType<double>())
    {
      if (param.getVec() == 0)
      {
        double val;
        if (isTempParam((*it).first) && param.getAutoConvertTemperature())
          val = getDefaultValue<double>(param) - CONSTCtoK;
        else
          val = getDefaultValue<double>(param);
        param_list.push_back(Param((*it).first, val));
      }
      else if (param.getVec() > 0)
      {
        if (param.getVec() == 1)
        {
          // This converts parameters link IC1, IC2 to just IC and type vector
          std::string vPar((*it).first.substr(0, (*it).first.size()-1));
          param_list.push_back(Param(vPar, "VECTOR"));
        }
        // We will also output IC1, IC2 as type double so they can
        // be specified as individual elements
        // This allows TC=a, b to also be specified as TC1=a TC2=b
        double val = getDefaultValue<double>(param);
        param_list.push_back(Param((*it).first, val));
      }
    }
    else if (param.isType<bool>())
    {
      if (param.getVec() == 0)
        param_list.push_back(Param((*it).first, getDefaultValue<bool>(param)));
      else if (param.getVec() > 0)
      {
        if (param.getVec() == 1)
        {
          std::string vPar((*it).first.substr(0, (*it).first.size()-1));
          param_list.push_back(Param(vPar, "VECTOR"));
        }
        // We will also output IC1, IC2 as type double so they can
        // be specified as individual elements
        // This allows TC=a, b to also be specified as TC1=a TC2=b
        bool val = getDefaultValue<bool>(param);
        param_list.push_back(Param((*it).first, val));
      }
    }
    else if (param.isType<int>())
    {
      if (param.getVec() == 0)
        param_list.push_back(Param((*it).first, getDefaultValue<int>(param)));
      else if (param.getVec() > 0)
      {
        if (param.getVec() == 1)
        {
          std::string vPar((*it).first.substr(0, (*it).first.size()-1));
          param_list.push_back(Param(vPar, "VECTOR"));
        }
        int val = getDefaultValue<int>(param);
        param_list.push_back(Param((*it).first, val));
      }
    }
    else if (param.isType<std::string>())
    {
      if (param.getVec() == 0)
        param_list.push_back(Param((*it).first, getDefaultValue<std::string>(param)));
      else if (param.getVec() > 0)
      {
        if (param.getVec() == 1)
        {
          std::string vPar((*it).first.substr(0, (*it).first.size()-1));
          param_list.push_back(Param(vPar, "VECTOR"));
        }
        std::string val = getDefaultValue<std::string>(param);
        param_list.push_back(Param((*it).first, val));
      }
    }
    else if (param.isType<std::vector<std::string> >())
    {
      Param vc((*it).first, std::vector<std::string>()); 
      param_list.push_back(vc);
    }
    else if (param.isType<std::vector<double> >())
    {
      Param vc((*it).first, std::vector<double>()); 
      param_list.push_back(vc);
    }
    else if (param.hasCompositeData())
    {
      Param vc2((*it).first, "VECTOR-COMPOSITE");
      vc2.setDefault(true);
      param_list.push_back(vc2);

      std::vector<Param> compositeParams;
      const ParametricData<CompositeParam> *c = param.getCompositeParametricData<CompositeParam>();

      if (c == 0)
      {
        Report::DevelFatal().in("populateParams") << "Vector-composite map for device type entity empty.";
      }

      // TODO: [DGB] I think when the Descriptor is refactored this will be clearer.  But this basically adds the
      //   type to the composite list with 'NAME' first.
      const ParametricData<CompositeParam> &d = *c;
      const ParameterMap &e = d.getMap();

      for (ParameterMap::const_iterator it4 = e.find("NAME"); it4 != e.end();) {
        const Descriptor &p = *(*it4).second;
        if (p.isType<double>())
          compositeParams.push_back(Param((*it4).first, getDefaultValue<double>(p)));
        else if (p.isType<bool>())
          compositeParams.push_back(Param((*it4).first, getDefaultValue<bool>(p)));
        else if (p.isType<int>())
          compositeParams.push_back(Param((*it4).first, getDefaultValue<int>(p)));
        else if (p.isType<std::string>())
          compositeParams.push_back(Param((*it4).first, getDefaultValue<std::string>(p)));
        if ((*it4).first == "NAME")
          it4 = e.begin();
        else
          it4++;
        if (it4 != e.end() && (*it4).first == "NAME")
          it4++;
      }

      composite_param_map[(*it).first] = compositeParams;
    }
    else
    {
      // Just skip these, like list of coupled inductors because not needed for metadata
      Xyce::dout() << "In final else clause of DeviceEntity::setParams().";
      if( param.isType<std::vector<std::string> >() )
        Xyce::dout() << " type is STR_VEC ";
      if( param.isType<std::vector<double> >() )
        Xyce::dout() << " type is DBLE_VEC ";
      Xyce::dout() << it->first << " this item is NOT being added to default parameter list." << std::endl;
    }
  }
}

} // namespace Device
} // namespace Xyce

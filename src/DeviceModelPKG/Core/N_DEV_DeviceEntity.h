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
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/03/00
//
//
//
//
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_DeviceEntity_h
#define Xyce_N_DEV_DeviceEntity_h

#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_UTL_fwd.h>

#include <N_UTL_NetlistLocation.h>
#include <N_DEV_Pars.h>
#include <N_DEV_InstanceName.h>

namespace Xyce {
namespace Device {

typedef std::map<std::string, std::vector<Param>, LessNoCase> CompositeParamMap;

void populateParams(const ParameterMap &parameter_map, std::vector<Param> &param_list, CompositeParamMap &composite_param_map);

//-----------------------------------------------------------------------------
// Class         : Depend
// Purpose       : Used to record information about dependent parameters
// Special Notes :
// Creator       : Dave Shirley
// Creation Date : 
//-----------------------------------------------------------------------------
///
///  The Depend struct is used to keep track of dependent parameters
///
struct Depend
{
  std::string                 name;      ///< parameter name
  Util::Expression *          expr;      ///< expression used comput value
  union resUnion
  {
    double *                result;
    std::vector<double> *   resVec;
  } resultU;                            ///< Holds a pointer to where the
                                        ///< parameter is stored.
  int                         vectorIndex; ///< Used if parameter is in a vector
  std::vector<double>         vals;     ///< values on which expression depends
  std::vector<std::string>    global_params; ///< global params on which
                                             ///< expression depends
  int                         n_vars, lo_var; 
  bool                     storeOriginal;    ///< true if original value stored
  int                      serialNumber;     ///< used if original value stored

  // Constructor
  Depend()
    : vectorIndex(-1), n_vars(0), lo_var(0)
  {};

};

// ERK.  this could be replaced by a lambda
struct MatchDependName
{
  MatchDependName(const std::string& name) : matchName_(name) {}
  bool operator()(const Depend & dep) const
  {
    return dep.name == matchName_;
  }
  private:
    const std::string& matchName_;
};

//-----------------------------------------------------------------------------
// Class         : DeviceEntity
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/11/02
//-----------------------------------------------------------------------------
class DeviceEntity : public ParameterBase
{
public:
  DeviceEntity(
     ParametricData<void> &    parametric_data,
     const SolverState &       solver_state,
     const DeviceOptions &     device_options,
     const NetlistLocation &         netlist_location);

private:
  DeviceEntity(const DeviceEntity &);                   ///< No copying
  DeviceEntity &operator=(const DeviceEntity &);        ///< No assignment

public:
  virtual ~DeviceEntity();

  virtual bool processParams() = 0;
  virtual bool processInstanceParams() = 0;
  virtual void processSuccessfulTimeStep();

  virtual CompositeParam *constructComposite(const std::string &composite_name, const std::string &param_name) 
  {
    return NULL;
  }

  bool setDefaultParam(double val, bool overrideOriginal=false);
  double getDefaultParam() const;

  bool scaleParam(const std::string & paramName, double val, double val0);
  bool scaleParam(const std::string & paramName, double val);
  bool scaleDefaultParam(double val);

  bool analyticSensitivityAvailable (const std::string & paramName);
  bool analyticSensitivityAvailableDefaultParam ();

  bool getAnalyticSensitivity ( const std::string & paramName,
                                std::vector<double> & dfdpVec,
                                std::vector<double> & dqdpVec,
                                std::vector<double> & dbdpVec,
                                std::vector<int> & FindicesVec,
                                std::vector<int> & QindicesVec,
                                std::vector<int> & BindicesVec );

  bool getAnalyticSensitivityDefaultParam ( std::vector<double> & dfdpVec,
                                std::vector<double> & dqdpVec,
                                std::vector<double> & dbdpVec,
                                std::vector<int> & FindicesVec,
                                std::vector<int> & QindicesVec,
                                std::vector<int> & BindicesVec );

  virtual bool getNumericalSensitivity ( const std::string & paramName,
                                std::vector<double> & dfdpVec,
                                std::vector<double> & dqdpVec,
                                std::vector<double> & dbdpVec,
                                std::vector<int> & FindicesVec,
                                std::vector<int> & QindicesVec,
                                std::vector<int> & BindicesVec ) = 0;

  bool getNumericalSensitivityDefaultParam ( std::vector<double> & dfdpVec,
                                std::vector<double> & dqdpVec,
                                std::vector<double> & dbdpVec,
                                std::vector<int> & FindicesVec,
                                std::vector<int> & QindicesVec,
                                std::vector<int> & BindicesVec );



//
  bool analyticMatrixSensitivityAvailable (const std::string & paramName);
  bool analyticMatrixSensitivityAvailableDefaultParam ();

  bool getAnalyticMatrixSensitivity ( const std::string & paramName,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs );

  bool getAnalyticMatrixSensitivityDefaultParam (
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs );

  virtual bool getNumericalMatrixSensitivity ( const std::string & paramName,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs ) = 0;

  bool getNumericalMatrixSensitivityDefaultParam (
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs );

//
  bool getAnalyticACSensitivityAvailable (const std::string & paramName);
  bool getAnalyticACSensitivityAvailableDefaultParam ();

  bool getAnalyticBSensVectorsforAC ( const std::string & paramName,
          std::vector< std::complex<double> > & dbdp,
          std::vector<int> &        BindicesVec);

  bool getAnalyticBSensVectorsforACDefaultParam (
          std::vector< std::complex<double> > & dbdp,
          std::vector<int> &        BindicesVec);

  virtual bool getNumericalBSensVectorsforAC ( const std::string & paramName,
          std::vector< std::complex<double> > & dbdp,
          std::vector<int> &        BindicesVec) = 0;

  bool getNumericalBSensVectorsforACDefaultParam (
          std::vector< std::complex<double> > & dbdp,
          std::vector<int> &        BindicesVec);
//
  bool setParam(const std::string & paramName, double val, bool overrideOriginal=false);
  bool getParam(const std::string & paramName, double & result) const;
  bool findParam(const std::string &param_name) const;
  void setupParamBreakpoints();
  bool getParamBreakpoints( std::vector<Util::BreakPoint> & );

  bool setParameterRandomExpressionTerms(const std::string & paramName, int opIndex, int astType, double value, bool override_original);

  bool updateGlobalAndDependentParameters (
      bool globalParameterChanged,
      bool timeChanged, 
      bool freqChanged);

  bool updateDependentParameters(double temp_tmp);
  bool updateDependentParameters();
  void applyDepSolnLIDs();

  double setDependentParameter(Util::Param &, double *, ParameterType::ExprAccess);
  double setDependentParameter(Util::Param &, std::vector<double> *, int , ParameterType::ExprAccess);
  void setDependentParameter(Util::Param & par, Depend & dependentParam, ParameterType::ExprAccess depend);

  void setDefaultParams() 
  {
    setDefaultParameters(*this, getParameterMap().begin(), getParameterMap().end(), devOptions_);
  }

  void setParams(const std::vector<Param> & params);

  // special handling for parameters contained within a vector-composite parameter
  void setParamFromVCParam(CompositeParam &composite_param,
                           const std::string &composite_name, 
                           const std::string &pName, 
                           const Param &ndParam);

public:
  bool given(const std::string & parameter_name) const;

  virtual std::ostream &printName(std::ostream &os) const = 0;

  void setDefaultParamName(const std::string &default_param_name) 
  {
    defaultParamName_ = default_param_name;
  }

  const std::string& getDefaultParamName() const { return defaultParamName_;}

  const std::vector<Depend> &getDependentParams() 
  {
    return dependentParams_;
  }

  void addDependentParameter(const Depend &param)
  {
    dependentParams_.push_back(param);
  }

  const DeviceOptions &getDeviceOptions() const 
  {
    return devOptions_;
  }

  const SolverState &getSolverState() const 
  {
    return solState_;
  }

  const NetlistLocation &netlistLocation() const 
  {
    return netlistLocation_;
  }

  const ParameterMap &getParameterMap() const 
  {
    return parametricData_.getMap();
  }

  void resetScaledParams()
  {
    dependentScaleParamExcludeMap_.clear();
  }

private:
  void escape(std::string &) const;
  void checkDepend(ParameterType::ExprAccess &);

private:
  std::string                 defaultParamName_;
  ParametricData<void> &      parametricData_;
  NetlistLocation             netlistLocation_;

  const SolverState &         solState_;
  Globals &                   globals_;

  const DeviceOptions &       devOptions_;
  std::vector<Depend>         dependentParams_;

  std::unordered_map <std::string, int> dependentParamExcludeMap_;
  std::unordered_map <std::string, int> dependentScaleParamExcludeMap_;

protected:
  std::vector<int>            expVarGIDs;
  std::vector<int>            expVarLIDs;
  std::vector<int>            expVarTypes;
  std::vector<std::string>    expVarNames;
  std::vector<double>         expVarVals;
  std::vector<double>         eVarVals;
};

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_DeviceEntity_h

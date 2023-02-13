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
// Purpose        : This file contains the device model base class.
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

#ifndef Xyce_N_DEV_DeviceModel_h
#define Xyce_N_DEV_DeviceModel_h

#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include <N_DEV_fwd.h>

#include <N_DEV_Device.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_DeviceInstance.h>

namespace Xyce {
namespace Device {

/**
 * @class DeviceModel N_DEV_DeviceModel.h
 *
 * @author Eric Keiter, SNL, Parallel Computational Sciences
 * @date   4/03/00
 */
class DeviceModel : public DeviceEntity
{
  enum mType {TEMP, DOSE};
  enum iType {LIN, QUAD, PWL};
  enum fitType {LINEAR_FIT, LOG_FIT};

public:
  /**
   * Add the parameter "TEMPMODEL" to the parametric_data.
   *
   * @param parametric_data
   */
  template<class T>
  static void initThermalModel(ParametricData<T> &parametric_data)
  {
    parametric_data.addPar("TEMPMODEL", "NONE", &DeviceModel::temperatureModel)
      .setCategory(CAT_CONTROL)
      .setDescription("Specifies the type of parameter interpolation over temperature");
  }

  /**
   * Add the parameter "DOSEMODEL" to the parametric_data.
   *
   * @param parametric_data
   */
  template<class T>
  static void initDoseModel(ParametricData<T> &parametric_data)
  {
    parametric_data.addPar("DOSEMODEL", "NONE", &DeviceModel::doseModel);
  }

  DeviceModel(
     const ModelBlock &        model_block,
     ParametricData<void> &    parametric_data,
     const FactoryBlock &      factory_block);

  virtual ~DeviceModel();

private:
  DeviceModel();
  DeviceModel(const DeviceModel &);
  DeviceModel &operator=(const DeviceModel &);

public:
  const std::string &getName() const
  {
    return name_;
  }

  void setModParams(const std::vector<Param> &params);

  virtual void forEachInstance(DeviceInstanceOp &op) const = 0;

  virtual std::ostream &printName(std::ostream &os) const;

  virtual std::ostream &printOutInstances(std::ostream &os) const = 0;

  /**
   * processParams
   *
   * @return true if parameter processing was successful
   */
  virtual bool processParams() = 0;

  /**
   * processInstanceParams
   *
   * @return true if parameter processing was successful
   */
  virtual bool processInstanceParams() = 0;

  void addBaseInstance(Xyce::Device::DeviceInstance *instance) 
  {
    baseInstanceContainer.push_back(instance);
  }

  virtual void setupBaseInstanceContainer();

  virtual bool clearTemperatureData ()
  {
    return true;
  }

  // ERK.  these 4 functions are used for interpolating model parameters w.r.t. 
  // temperature and/or dose.
  void saveParams ();
  bool interpolateTNOM (double);
  bool interpolateDOSE (double);
  void restoreParams ();

  virtual bool getBinPrefixFlag ()
  {
    return false;
  }


  bool getNumericalSensitivity ( const std::string & paramName,
                                std::vector<double> & dfdpVec,
                                std::vector<double> & dqdpVec,
                                std::vector<double> & dbdpVec,
                                std::vector<int> & FindicesVec,
                                std::vector<int> & QindicesVec,
                                std::vector<int> & BindicesVec );

  bool getNumericalMatrixSensitivity ( const std::string & paramName,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs );

  bool getNumericalBSensVectorsforAC ( const std::string & paramName,
    std::vector< std::complex<double> > & dbdp,
    std::vector<int> &        BindicesVec) {return true;}

private:

  // these two functions are used for interpolating model parameters w.r.t. 
  // temperature and/or dose
  bool interpolated ();
  bool interpolate (double);

public:
  int getLevel() const
  {
    return level_;
  }

  void setLevel(int level)
  {
    level_ = level;
  }

  const std::string &getType() const
  {
    return type_;
  }

public:
  std::vector<Xyce::Device::DeviceInstance *> baseInstanceContainer;

private:
  std::string                           name_;
  std::string                           type_;
  int                                   level_;
  std::string                           temperatureModel;
  std::string                           doseModel;
  mType                                 iModel;
  iType                                 iMethod;
  double                                base_temp;
  std::map<std::string, int>            fitMap;
  std::vector<double DeviceEntity::*>   fitParams;
  std::vector<double>                   oldParams;
  std::vector<double>                   base;
  std::vector< std::vector<double> >    fit;
  std::vector<double>                   min_par;
  std::vector<double>                   max_par;
  std::vector<fitType>                  parType;
};

} // namespace Device
} // namespace Xyce

#endif


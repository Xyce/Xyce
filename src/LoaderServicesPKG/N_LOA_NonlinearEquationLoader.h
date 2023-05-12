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
// Purpose        : This file contains the interface class that sits between
//                  the nonlinear solver and (ultimately) the time integrator.
//                  For transient calculations only the time integrator knows
//                  what to sum into the Jacobian and Residual vector.
//
// Special Notes  :
//
// Creator        : Todd Coffey, SNL
//
// Creation Date  : 07/29/08
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_LOA_NonlinearEquationLoader_H
#define Xyce_LOA_NonlinearEquationLoader_H

#include <N_DEV_fwd.h>
#include <N_LAS_fwd.h>
#include <N_LOA_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TIA_fwd.h>
#include <N_UTL_fwd.h>

#include <N_UTL_Stats.h> 

namespace Xyce {
namespace Loader {

//-----------------------------------------------------------------------------
// Class         : NonlinearEquationLoader
//
// Purpose        : This class contains the interface class that sits between
//                  the nonlinear solver and (ultimately) the time integrator.
//                  For transient calculations only the time integrator knows
//                  what to sum into the Jacobian and Residual vector.
//
// Special Notes :
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
class NonlinearEquationLoader
{
public:

  // Default constructor
  NonlinearEquationLoader(
    TimeIntg::DataStore &                       ds,
    Loader &                                    loader,
    Device::DeviceMgr &                         device_manager,
    TimeIntg::WorkingIntegrationMethod &        wim,
    bool                                        daeStateDerivFlag);

  // Destructor
  virtual ~NonlinearEquationLoader();

  // // Method which is called to load the nonlinear Jacobian matrix.
  bool loadJacobian();

  // // Method which is called to apply the nonlinear Jacobian matrix.
  bool applyJacobian (const Linear::Vector& input, Linear::Vector& result);

  // // Method which is called to load the nonlinear residual (RHS) vector.
  bool loadRHS();

  void setSeparateLoadFlag (bool flag);
  bool getSeparateLoadFlag ();

  bool loadSensitivityResiduals ();

  bool updateSensitivityHistoryAdjoint();
  bool updateSensitivityHistoryAdjoint2();
  bool loadFunctionDerivativesForTranAdjoint ();
  bool loadSparseFunctionDerivativesForTranAdjoint ();
  bool loadAdjointSensitivityResidual ();

  bool isLinearSystem() const;

  // two-level newton functions:
  int  enablePDEContinuation();
  bool disablePDEContinuation ();

  void getNumInterfaceNodes (std::vector<int> & numINodes);
  bool loadCouplingRHS   (int iSubProblem, int iCouple, Linear::Vector * dfdvPtr);
  bool calcCouplingTerms (int iSubProblem, int iCouple, const Linear::Vector * dxdvPtr);

  // Gets the nonlinear residual load time.
  double getResidualTime() { return residualTime_; }

  // Gets the nonlinear Jacobian load time.
  double getJacobianTime() { return jacobianTime_; }

  // Get convergence info from devices
  bool allDevicesConverged(Parallel::Machine comm);

  // Get convergence info from inner-solves
  bool innerDevicesConverged(Parallel::Machine comm);

  // toggle if the initial junction voltages are applied to devices or not
  void setDisableInitJctFlags(bool flag);

  // Function for determining if an analytic sensitivity (df/dp, dq/dp, db/dp) is available.
  bool analyticSensitivitiesAvailable (const std::string & name);
  bool numericalSensitivitiesAvailable (const std::string & name);

  void getAnalyticSensitivities(
      std::string &             name, 
      std::vector<double> &     dfdpVec, 
      std::vector<double> &     dqdpVec,
      std::vector<double> &     dbdpVec,
      std::vector<int> &        FindicesVec,
      std::vector<int> &        QindicesVec,
      std::vector<int> &        BindicesVec) const;

  void getNumericalSensitivities(
      std::string &             name, 
      std::vector<double> &     dfdpVec, 
      std::vector<double> &     dqdpVec,
      std::vector<double> &     dbdpVec,
      std::vector<int> &        FindicesVec,
      std::vector<int> &        QindicesVec,
      std::vector<int> &        BindicesVec) const;

  bool analyticBVecSensAvailable(const std::string & name);
  bool numericalBVecSensAvailable(const std::string & name);

  bool analyticMatrixSensitivitiesAvailable(const std::string & name);
  bool numericalMatrixSensitivitiesAvailable(const std::string & name);

  void getAnalyticBSensVectorsforAC (const std::string & name,
          std::vector< std::complex<double> > &     dbdp,
          std::vector<int> &        BindicesVec) const;

  void getNumericalBSensVectorsforAC (const std::string & name,
          std::vector< std::complex<double> > &     dbdp,
          std::vector<int> &        BindicesVec) const;

  void getAnalyticMatrixSensitivities(
      std::string & name,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs) const;

  void getNumericalMatrixSensitivities(
      std::string & name,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs) const;


  bool setParam(std::string & name, double val, bool overrideOriginal=false);

  bool getParamAndReduce(Parallel::Machine comm, const std::string & name, double & val) const;

  const Loader &getLoader() const {
    return loader_;
  }

  // Functions needed by the NEW (power node) 2-level algorithm:
  void homotopyStepSuccess(const std::vector<std::string> & paramNames, const std::vector<double> & paramVals);
  void homotopyStepFailure();

  // voltage limiter toggle functions
  bool getVoltageLimiterStatus();
  void setVoltageLimiterStatus(bool voltageLimterStatus);

#if 0
  virtual void updateDependentParams ();
#endif
  virtual void resetScaledParams();

private:
  Util::Timer * residualTimerPtr_;
  Util::Timer * jacobianTimerPtr_;
  double                                residualTime_;
  double                                jacobianTime_;
  bool                                  daeStateDerivFlag_;
  TimeIntg::DataStore &                 ds_;
  Loader &                              loader_;
  TimeIntg::WorkingIntegrationMethod &  wim_;
  Device::DeviceMgr &                   deviceManager_; ///< Device manager
};

} // namespace Loader
} // namespace Xyce

#endif


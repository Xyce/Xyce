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
// Purpose       : This file defines the classes for the "no time integration"
//                 method.
//
// Special Notes :
//
// Creator       : Eric Keiter
//
// Creation Date : 7/21/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_NO_TIME_INTEGRATION_H
#define Xyce_N_TIA_NO_TIME_INTEGRATION_H

// ---------- Standard Includes ----------
#include <N_UTL_Math.h>

// ----------   Xyce Includes   ----------
#include <N_TIA_fwd.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_TIA_StepErrorControl.h>

#define MATRIX_FAILSAFE 1
#define NUM_LIMIT  1.0e-20

namespace Xyce {
namespace TimeIntg {

//-----------------------------------------------------------------------------
// Class         : NoTimeIntegration
// Purpose       : Class objects for use during DC analysis (applies only when
//                 all time derivatives are set to 0) (derived from
//                 TimeIntegrationMethod)
// Special Notes :
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
class NoTimeIntegration : public TimeIntegrationMethod
{
public:
  static const int type = Xyce::TimeIntg::methodsEnum::NO_TIME_INTEGRATION;
  static const char *name;

  static TimeIntegrationMethod *factory(const TIAParams &tia_params, StepErrorControl &step_error_control, DataStore &data_store);

  NoTimeIntegration(
    const TIAParams & tiaP,
    StepErrorControl & secTmp,
    DataStore & dsTmp);

  ~NoTimeIntegration();

  const char *getName() const {
    return "None";
  }

  // Predict solution at next time point (No Integration)
  void obtainPredictor()
  {
    ds.usePreviousSolAsPredictor ();
  }

  // Evaluate the predictor derivative formula (No Integration)
  void obtainPredictorDeriv()
  {}

  // Evaluate the corrector derivative formula (No Integration)
  void obtainCorrectorDeriv();

  // Compute an estimate of the error in the integration step (No Integration)
  double computeErrorEstimate() const { return 0.0; }

  // Interpolate solution approximation at prescribed time point (No
  // Integration)
  bool interpolateSolution(double timepoint, Linear::Vector * tmpSolVectorPtr, std::vector<Linear::Vector*> & historyVec )
  {
    return false;
  }

  bool printOutputSolution(
    Analysis::OutputMgrAdapter &        outputManagerAdapter,
    const TIAParams &                   tia_params,
    const double                        time,
    Linear::Vector *                      solnVecPtr,
    const bool                          doNotInterpolate,
    const std::vector<double> &         outputInterpolationTimes,
    bool                                skipPrintLineOutput);

  bool saveOutputSolution(
    Parallel::Machine                   comm,
    IO::InitialConditionsManager &      initial_conditions_manager,
    const NodeNameMap &                 node_name_map,
    const TIAParams &                   tia_params,
    Linear::Vector *                    solnVecPtr,
    const double                        saveTime,
    const bool                          doNotInterpolate);

  // Gets the time-integration order (No Integration)
  int getOrder() const
  {
    return 0;
  }
  
  int getUsedOrder() const
  {
    return 0;
  }

  int getMethod() const
  {
    return type;
  }

  int getNumberOfSteps() const
  {
    return 0;
  }

  int getNscsco() const
  {
    return 0;
  }

  int getMaxOrder() const
  {
    return -2;
  }

  void getInitialQnorm(TwoLevelError & tle) const;
  void getTwoLevelError(TwoLevelError & tle) const;

  double partialTimeDeriv() const
  {
#ifdef MATRIX_FAILSAFE
    return NUM_LIMIT;
#else
    return 0.0;
#endif
  }

  // Gets the leading coefficient for the specified time-integration method.
  double getLeadingCoeff() const { return leadingCoeff; }

  // sets the leading coefficient for the specified time-integration method.
  void setLeadingCoeff(double & LC) { leadingCoeff = LC; }

  // 03/08/04 erkeite:  New functions necessary new-DAE:
  // Evaluate residual for nonlinear solver
  void obtainResidual();

  // obtain residuals for transient sensitivity analysis
  void obtainSensitivityResiduals();

  // obtain function derivatives for transient adjoint sensitivity analysis
  void obtainFunctionDerivativesForTranAdjoint ();

  // obtain function derivatives for transient adjoint sensitivity analysis, sparse version
  void obtainSparseFunctionDerivativesForTranAdjoint ();

  // obtain residuals for transient adjoint sensitivity analysis
  void obtainAdjointSensitivityResidual ();

  // Evaluate Jacobian for nonlinear solver
  void obtainJacobian();

  // Apply Jacobian for nonlinear solver
  void applyJacobian(const Linear::Vector& input, Linear::Vector& result);

  void initialize(const TIAParams &tia_params)
  {}

  void initializeAdjoint(int index)
  {}

  void setTwoLevelTimeInfo()
  {}
  
  void rejectStep(const TIAParams &tia_params);
  void completeStep(const TIAParams &tia_params);
  void completeAdjointStep(const TIAParams &tia_params) {}

  // Compute "delta" vectors in data store.
  void stepLinearCombo() { ds.stepLinearCombo(); }

  // Restart methods.
  bool getSolnVarData( const int & gid, std::vector<double> & varData )
    { return ds.getSolnVarData( gid, varData ); }
  bool setSolnVarData( const int & gid, const std::vector<double> & varData )
    { return ds.setSolnVarData( gid, varData ); }
  bool getStateVarData( const int & gid, std::vector<double> & varData )
    { return ds.getStateVarData( gid, varData ); }
  bool setStateVarData( const int & gid, const std::vector<double> & varData )
    { return ds.setStateVarData( gid, varData ); }
  bool getStoreVarData( const int & gid, std::vector<double> & varData )
    { return ds.getStoreVarData( gid, varData ); }
  bool setStoreVarData( const int & gid, const std::vector<double> & varData )
    { return ds.setStoreVarData( gid, varData ); }

private:
  double                alphas;         ///< $\alpha_s$ fixed-leading coefficient of this BDF method
  DataStore &           ds;             ///< Reference to the TIA data-store object.
  StepErrorControl &    sec;            ///< Reference to step-error control object
  double                leadingCoeff;   ///< Time-integration method leading coefficient value.
};

} // namespace TimeIntg
} // namespace Xyce

#endif // Xyce_N_TIA_NO_TIME_INTEGRATION_H

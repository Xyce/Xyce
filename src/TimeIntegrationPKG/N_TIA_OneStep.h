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


//-----------------------------------------------------------------------------
//
// Purpose       : This file defines the classes for trapezoidal method 
//                 order 1-2 method.
//
// Special Notes :
//
// Creator       : Ting Mei
//
// Creation Date : 10/31/07
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_ONE_STEP_H
#define Xyce_N_TIA_ONE_STEP_H


// ---------- Standard Includes ----------
// ----------   Xyce Includes   ----------
#include <N_TIA_DataStore.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_PDS_ParMap.h>
#include <N_UTL_Math.h>

namespace Xyce {
namespace TimeIntg {

//-----------------------------------------------------------------------------
// Class         : OneStep
// Purpose       : OneStep Formula Integration Class
//                (derived from TimeIntegrationMethod)
// Special Notes :
// Creator       : Ting Mei, SNL
// Creation Date : 10/31/07
//-----------------------------------------------------------------------------
class OneStep : public TimeIntegrationMethod
{
public:
  static const int type = 7;
  static const char *name;

  static TimeIntegrationMethod *factory(const TIAParams &tia_params, StepErrorControl &step_error_control, DataStore &data_store);

  OneStep(
    const TIAParams & tiaP,
    StepErrorControl & secTmp,
    DataStore & dsTmp);

  ~OneStep() 
  {}

  const char *getName() const {
    return "Onestep: Trapezoidal";
  }

  // Predict solution at next time point.
  void obtainPredictor() ;

  // Evaluate the predictor derivative formula.
  void obtainPredictorDeriv() {return;}

  // Evaluate the corrector derivative formula.
  void obtainCorrectorDeriv() {return;}

  // Compute an estimate of the error in the integration step.
  double computeErrorEstimate() const { return sec.ck_*ds.WRMS_errorNorm(); }

  // Interpolate solution approximation at prescribed time point.
  bool interpolateSolution(double timepoint,
                  Linear::Vector * tmpSolVectorPtr, std::vector<Linear::Vector*> & historyVec);

  // Interpolate MPDE solution approximation at prescribed time point.
  bool interpolateMPDESolution(std::vector<double>& timepoint,
                  Linear::Vector * tmpSolVectorPtr) ;
  
  // Print output using interpolation when order is high
  bool printOutputSolution(
    Analysis::OutputMgrAdapter &        outputManagerAdapter,
    const TIAParams &                   tia_params,
    const double                        time,
    Linear::Vector *                      solnVecPtr,
    const bool                          doNotInterpolate,
    const std::vector<double> &         outputInterpolationTimes,
    bool                                skipPrintLineOutput );

  // Print transient output from MPDE simulation
  bool printMPDEOutputSolution(
                  Analysis::OutputMgrAdapter & outputManagerAdapter,
                  const double time,
	                Linear::Vector * solnVecPtr,
                  const std::vector<double> & fastTimes );

  // Print transient output from WaMPDE simulation
  bool printWaMPDEOutputSolution(
                  Analysis::OutputMgrAdapter & outputManagerAdapter,
                  const double time,
	                Linear::Vector * solnVecPtr,
                  const std::vector<double> & fastTimes,
                  const int phiGID );

  // .SAVE output using interpolation when order is high
  bool saveOutputSolution(
    Parallel::Machine                   comm,
    IO::InitialConditionsManager &      initial_conditions_manager,
    const NodeNameMap &                 node_name_map,
    const TIAParams &                   tia_params,
    Linear::Vector *                    solnVecPtr,
    const double                        saveTime,
    const bool                          doNotInterpolate); 

  // Gets the time-integration order.
  int getOrder() const
  {
    return sec.currentOrder_;
  }

  int getNumberOfSteps() const
  {
    return sec.numberOfSteps_;
  }

  int getUsedOrder() const
  {
    return sec.usedOrder_;
  }

  int getNscsco() const
  {
    return sec.nscsco_;
  }

  int getMaxOrder() const
  {
    return sec.maxOrder_;
  }

  void getInitialQnorm(TwoLevelError & tle) const;
  void getTwoLevelError(TwoLevelError & tle) const;

  // This is for new-DAE.
  void updateStateDeriv();
  
  // calculates dQ/dt component of lead current Q vector and adds it to the lead current vector
  void updateLeadCurrentVec();

  double partialTimeDeriv() const;

  // Gets the leading coefficient for the specified time-integration method.
  double getLeadingCoeff() const
  {
    return leadingCoeff;
  }

  // sets the leading coefficient for the specified time-integration method.
  void setLeadingCoeff(double & LC)
  {
    leadingCoeff = LC;
  }

  // Evaluate corrector residual for nonlinear solver
  void obtainResidual();

  // obtain residuals for transient sensitivity analysis
  void obtainSensitivityResiduals();

  // obtain function derivatives for transient adjoint sensitivity analysis
  void obtainFunctionDerivativesForTranAdjoint ();

  // obtain function derivatives for transient adjoint sensitivity analysis, sparse version
  void obtainSparseFunctionDerivativesForTranAdjoint ();

  // obtain residuals for transient adjoint sensitivity analysis
  void obtainAdjointSensitivityResidual ();

  // Evaluate corrector Jacobian for nonlinear solver
  void obtainJacobian();

  // Update history array after a successful step 
  void updateHistory();

  // Restore history array after a failed step
  void restoreHistory();

  // Update method coefficients
  void updateCoeffs();

  // Update method coefficients
  void updateAdjointCoeffs();

  // Initialize method with initial solution & step-size
  void initialize(const TIAParams &tia_params);

  // Initialize method for adjoint sensitivities
  void initializeAdjoint(int index);

  // setup 2-level data.
  void setTwoLevelTimeInfo();

  // Reject a step(this turns back history and chooses new order & step-size)
  void rejectStep(const TIAParams &tia_params);

  // Reject a step, but only restore the history, as the new step will be
  // imposed by an outside program. 
  void rejectStepForHabanero();

  // Complete a step(this updates history and chooses new order & step-size)
  void completeStep(const TIAParams &tia_params);

  // Complete an adjoint step(this updates history and chooses new order & step-size)
  void completeAdjointStep(const TIAParams &tia_params);

  // Compute "delta" vectors in data store.
  void stepLinearCombo() { ds.stepLinearCombo(); }

  // Restart methods.
  // Save data local to this integration method.
  bool getSolnVarData( const int & gid, std::vector<double> & varData ); 
  bool setSolnVarData( const int & gid, const std::vector<double> & varData ); 
  bool getStateVarData( const int & gid, std::vector<double> & varData ); 
  bool setStateVarData( const int & gid, const std::vector<double> & varData ); 
  bool getStoreVarData( const int & gid, std::vector<double> & varData ); 
  bool setStoreVarData( const int & gid, const std::vector<double> & varData ); 
  
private:
  // Initialize arrays related to sensitivities
  void initializeSensitivities();

  // Predict solution at next time point.
  void obtainSensitivityPredictors();

  // Update sensitivity history arrays after a successful step 
  void updateSensitivityHistory();

  // Update adjoint sensitivity history arrays after a successful step 
  void updateAdjointSensitivityHistory();

  // Check whether to reduce order independent of local error test
  void checkReduceOrder();

private:
  // used for second order interpolation where the last two times steps are needed.
  // lastTimeStep is stored in StepErrorControl and the last, last time step is stored
  // here as it's only needed by trap.  Called timeStepForHistory2 as it's the 
  // time step associated with the step in xHistory[2]
  double                timeStepForHistory2_;
  double                timept_;        ///< Keep track of last interpolation point in printMPDEOutputSolution
  DataStore &           ds;             ///< Reference to the TIA data-store object.
  StepErrorControl &    sec;            ///< Reference to step-error control object
  double                leadingCoeff;   ///< Time-integration method leading coefficient value.
};

} // namespace TimeIntg
} // namespace Xyce

#endif     //Xyce_N_TIA_ONE_STEP_H


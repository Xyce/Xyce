//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose       : This file defines the classes for the time integration
//                 methods --- the "interface base class" along with the
//                 accompanying integration methods classes which can be
//                 used in the time integration algorithm.
//
// Special Notes :
//
// Creator       : Buddy Watts
//
// Creation Date : 6/1/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_TimeIntegrationMethods_h
#define Xyce_N_TIA_TimeIntegrationMethods_h

// ---------- Standard Includes ----------
#include <list>

// ----------   Xyce Includes   ----------
#include <N_ANP_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TIA_fwd.h>
#include <N_UTL_fwd.h>

namespace Xyce {
namespace TimeIntg {

//-----------------------------------------------------------------------------
// Class         : TimeIntegrationMethod
//
// Purpose       : This is the integration methods interface class, from which
//                 specific integration methods (such as Gear, trap, etc) are
//                 derrived.
//
// Special Notes :
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
class TimeIntegrationMethod
{
protected:
  TimeIntegrationMethod()
  {}

public:
  virtual ~TimeIntegrationMethod()
  {}

  virtual const char *getName() const = 0;

  // Predict solution at next time point.
  virtual void obtainPredictor() = 0;

  // Evaluate the predictor derivative formula.
  virtual void obtainPredictorDeriv() = 0;

  // Evaluate the corrector derivative formula.
  virtual void obtainCorrectorDeriv() = 0;

  // Compute an estimate of the error in the integration step.
  virtual double computeErrorEstimate() const = 0;

  // Interpolate solution, state or store approximation at prescribed time point.
  virtual bool interpolateSolution(double timepoint, Linear::Vector * tmpSolVectorPtr, std::vector<Linear::Vector*> & historyVec) = 0;

  // Print output using interpolation when order is high
  virtual bool printOutputSolution(
    Analysis::OutputMgrAdapter &        outputManagerAdapter,
    const TIAParams &                   tia_params,
    const double                        time,
    Linear::Vector *                    solnVecPtr,
    const bool                          doNotInterpolate,
    const std::vector<double> &         outputInterpolationTimes,
    bool                                skipPrintLineOutput) = 0;

  // Print MPDE output using local interpolation methods
  virtual bool printMPDEOutputSolution(
    Analysis::OutputMgrAdapter &        outputManagerAdapter,
    const double                        time,
    Linear::Vector *                    solnVecPtr,
    const std::vector<double> &         fastTimes)
  {
    return false;
  }

  // Print WaMPDE output using local interpolation methods
  virtual bool printWaMPDEOutputSolution(
    Analysis::OutputMgrAdapter &        outputManagerAdapter,
    const double                        time,
    Linear::Vector *                    solnVecPtr,
    const std::vector<double> &         fastTimes,
    const int                           phiGID)
  {
    return false;
  }

  // .SAVE output using interpolation when order is high
  virtual bool saveOutputSolution(
    Parallel::Machine                   comm,
    IO::InitialConditionsManager &      initial_conditions_manager,
    const NodeNameMap &                 node_name_map,
    const TIAParams &                   tia_params,
    Linear::Vector *                    solnVecPtr,
    const double                        saveTime,
    const bool                          doNotInterpolate) = 0;

  // Gets the time-integration order.
  virtual int getOrder() const = 0;
  virtual int getNumberOfSteps() const = 0;
  virtual int getUsedOrder() const = 0;
  virtual int getNscsco() const = 0;
  virtual int getMaxOrder() const = 0;                  ///< Return max order of method (this should obey user option maxorder)

  virtual void getInitialQnorm(TwoLevelError & tle) const = 0;
  virtual void getTwoLevelError(TwoLevelError & tle) const = 0;

  // This is for new-DAE.
  virtual void updateStateDeriv ()
  {}

  // calculates dQ/dt component of lead current vector and adds it to the lead current vector
  virtual void updateLeadCurrentVec ()
  {}

  // Returns the current partial time derivative for either the solution or
  // state vector.
  virtual double partialTimeDeriv() const = 0;

  // Gets the leading coefficient for the specified time-integration method.
  virtual double getLeadingCoeff() const = 0;

  // sets the leading coefficient for the specified time-integration method.
  virtual void setLeadingCoeff(double & LC) = 0;

  // Evaluate residual for nonlinear solver
  virtual void obtainResidual() = 0;

  // obtain residuals for transient sensitivity analysis
  virtual void obtainSensitivityResiduals() = 0;

  // obtain function derivatives for transient adjoint sensitivity analysis
  virtual void obtainFunctionDerivativesForTranAdjoint () = 0;

  // obtain function derivatives for transient adjoint sensitivity analysis, sparse version
  virtual void obtainSparseFunctionDerivativesForTranAdjoint () = 0;

  // obtain residuals for transient adjoint sensitivity analysis
  virtual void obtainAdjointSensitivityResidual () = 0;

  // Evaluate Jacobian for nonlinear solver
  virtual void obtainJacobian()
  {}

  // Apply Jacobian for nonlinear solver
  virtual void applyJacobian(const Linear::Vector& input, Linear::Vector& result)
  {}

  // Update history array after a successful step
  virtual void updateHistory()
  {}

  // Update history array for adjoints
  virtual void updateSensitivityHistoryAdjoint()
  {}
  virtual void updateSensitivityHistoryAdjoint2()
  {}

  // Update adjoint sensitivity history array after a successful step
  virtual void updateAdjointSensitivityHistory()
  {}


  // Restore history array after a failed step
  virtual void restoreHistory()
  {}

  // Update method coefficients
  virtual void updateCoeffs()
  {}

  // Update method coefficients, adjoint
  virtual void updateAdjointCoeffs()
  {}

  // Initialize method with initial solution & step-size
  virtual void initialize(const TIAParams &tia_params) = 0;

  // Initialize method for adjoint sensitivities
  virtual void initializeAdjoint (int index) = 0;

  // setup 2-level data.
  virtual void setTwoLevelTimeInfo() = 0;

  // Reject a step (this turns back history and chooses new order & step-size)
  virtual void rejectStep(const TIAParams &tia_params) = 0;

  // Reject a step, but only restore the history, as the new step will be
  // imposed by an outside program.
  virtual void rejectStepForHabanero()
  {}

  // Complete a step (this updates history and chooses new order & step-size)
  virtual void completeStep(const TIAParams &tia_params) = 0;

  // Complete an adjoint step (this updates history and chooses new order & step-size)
  virtual void completeAdjointStep(const TIAParams &tia_params) = 0;

  // Update "delta" vectors in the data store and local to this integration method.
  virtual void stepLinearCombo() = 0; 

  // Restart methods.
  // Save data local to this integration method.
  virtual bool getSolnVarData( const int & gid, std::vector<double> & varData ) = 0;
  virtual bool setSolnVarData( const int & gid, const std::vector<double> & varData ) = 0;
  virtual bool getStateVarData( const int & gid, std::vector<double> & varData ) = 0;
  virtual bool setStateVarData( const int & gid, const std::vector<double> & varData ) = 0;
  virtual bool getStoreVarData( const int & gid, std::vector<double> & varData ) = 0;
  virtual bool setStoreVarData( const int & gid, const std::vector<double> & varData ) = 0;

};

} // namespace TimeIntg
} // namespace Xyce

#endif // Xyce_N_TIA_TimeIntegrationMethods_h

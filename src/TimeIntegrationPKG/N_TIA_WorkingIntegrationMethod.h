//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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

#ifndef Xyce_N_TIA_WorkingIntegrationMethods_h
#define Xyce_N_TIA_WorkingIntegrationMethods_h

// ---------- Standard Includes ----------
#include <iostream>
#include <N_UTL_Math.h>
#include <list>

// ----------   Xyce Includes   ----------
#include <N_ANP_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_TIA_fwd.h>
#include <N_UTL_fwd.h>

#include <N_UTL_Stats.h>

namespace Xyce {
namespace TimeIntg {

typedef TimeIntegrationMethod *(*Factory)(const TIAParams &tia_params, StepErrorControl &step_error_control, DataStore &data_store);

void registerTimeIntegrationMethods();

void registerFactory(int type, const char *name, Factory factory);

template <class T>
inline
void
registerTimeIntegrationMethod()
{
  registerFactory(T::type, T::name, T::factory);
}


const char *getTimeIntegrationName(int type);

TimeIntegrationMethod *createTimeIntegrationMethod(int type, const TIAParams & tia_params, StepErrorControl & step_error_control, DataStore & data_store);

//-----------------------------------------------------------------------------
// Class         : WorkingIntegrationMethod
// Purpose       : This class provides a way for obtaining a specific
//                 working integration method and its associated data items.
// Special Notes :
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
class WorkingIntegrationMethod
{
public:
//  WorkingIntegrationMethod( );
  WorkingIntegrationMethod(Stats::Stat parent_stat);

  virtual ~WorkingIntegrationMethod();

  // Method which creates a time-integration method.
  void createTimeIntegMethod(
    int                 type,
    const TIAParams &   tia_params,
    StepErrorControl &  step_error_control,
    DataStore &         data_store);

  // accessors to the integration method object functions:
  bool isTimeIntegrationMethodCreated();

  double partialTimeDeriv() const;
  void obtainPredictor();
  void obtainPredictorDeriv();
  void obtainCorrectorDeriv();

  int getOrder() const;
  int getUsedOrder() const;
  int getMethod() const;
  int getNumberOfSteps() const;
  int getNscsco() const;
  void getInitialQnorm (TwoLevelError & tle) const;
  void getTwoLevelError(TwoLevelError & tle) const;
  void setTwoLevelTimeInfo();
  void updateCoeffs();
  void updateAdjointCoeffs();
  void rejectStepForHabanero ();
  void initialize(const TIAParams &tia_params);
  void initializeAdjoint(int index);
  void completeStep(const TIAParams &tia_params);
  void completeAdjointStep(const TIAParams &tia_params);
  void rejectStep(const TIAParams &tia_params);
  double computeErrorEstimate() const;
  void updateStateDeriv ();
  void updateLeadCurrent ();
  void updateLeadCurrentVec ();
  void obtainResidual();
  void obtainSensitivityResiduals();
  void loadFinalSensitivityDerivatives();
  void updateSensitivityHistoryAdjoint();
  void updateSensitivityHistoryAdjoint2();
  void obtainFunctionDerivativesForTranAdjoint ();
  void obtainSparseFunctionDerivativesForTranAdjoint ();
  void obtainAdjointSensitivityResidual ();
  void obtainJacobian();
  void applyJacobian(const Linear::Vector& input, Linear::Vector& result);

  bool printMPDEOutputSolution(
    Analysis::OutputMgrAdapter &        outputManagerAdapter,
    const double                        time,
    Linear::Vector *                      solnVecPtr,
    const std::vector<double> &         fastTimes );

  bool printWaMPDEOutputSolution(
    Analysis::OutputMgrAdapter &        outputManagerAdapter,
    const double                        time,
    Linear::Vector *                      solnVecPtr,
    const std::vector<double> &         fastTimes,
    const int                           phiGID );

  bool printOutputSolution(
    Analysis::OutputMgrAdapter &        outputManagerAdapter,
    const TIAParams &                   tia_params,
    const double                        time,
    Linear::Vector *                      solnVecPtr,
    const bool                          doNotInterpolate,
    const std::vector<double> &         outputInterpolationTimes,
    bool                                skipPrintLineOutput );

  bool saveOutputSolution(
    Parallel::Machine                   comm,
    IO::InitialConditionsManager &      initial_conditions_manager,
    const NodeNameMap &                 node_name_map,
    const TIAParams &                   tia_params,
    Linear::Vector *                    solnVecPtr,
    const double                        saveTime,
    const bool                          doNotInterpolate);

  // Update "delta" vectors in the data store an in each integration method.
  void stepLinearCombo ();

  bool getSolnVarData( const int & gid, std::vector<double> & varData ); 
  bool setSolnVarData( const int & gid, const std::vector<double> & varData );
  bool getStateVarData( const int & gid, std::vector<double> & varData );
  bool setStateVarData( const int & gid, const std::vector<double> & varData );
  bool getStoreVarData( const int & gid, std::vector<double> & varData );
  bool setStoreVarData( const int & gid, const std::vector<double> & varData );

private:
  TimeIntegrationMethod *       timeIntegrationMethod_;         ///< Pointer to the integration method.
  bool                          jacLimitFlag;
  double                        jacLimit;

  Stats::Stat                   timeIntegratorStat_;
  Stats::Stat                   predictorStat_; 
  Stats::Stat                   completeStepStat_;
  Stats::Stat                   rejectStepStat_;
  Stats::Stat                   updateCoefStat_;
  Stats::Stat                   residualStat_;
  Stats::Stat                   jacobianStat_;
  Stats::Stat                   initializeStat_;

  Stats::Stat                   updateLeadStat_; 
};

} // namespace TimeIntg
} // namespace Xyce

#endif // Xyce_N_TIA_WorkingIntegrationMethods_h

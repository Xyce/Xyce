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
// Purpose       : This file contains the functions which define the
//		             time integration methods classes.
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

#include <Xyce_config.h>

#include <N_ERH_ErrorMgr.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TIAParams.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_TIA_WorkingIntegrationMethod.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace TimeIntg {

namespace {

typedef std::map<int, std::pair<const char *, Factory> > Registry;

Registry &
getRegistry() 
{
  static Registry s_registry;

  return s_registry;
}

} // namespace <unnamed>

void
registerFactory(int type, const char *name, Factory factory)
{
  std::pair<Registry::iterator, bool> result = getRegistry().insert(Registry::value_type(type, std::pair<const char *, Factory>(name, factory)));
  if (!result.second && name != (*result.first).second.first)
    Report::DevelFatal0() << "Time integration factory " << type << " named " << name << " already registered with name " << (*result.first).second.first;
}

TimeIntegrationMethod *
createTimeIntegrationMethod(
  int                   type,
  const TIAParams &     tia_params,
  StepErrorControl &    step_error_control,
  DataStore &           data_store)
{
  Registry::const_iterator it = getRegistry().find(type);
  if (it == getRegistry().end())
    return 0;

  return (*(*it).second.second)(tia_params, step_error_control, data_store);
}

const char *
getTimeIntegrationName(int type) 
{
  Registry::const_iterator it = getRegistry().find(type);
  if (it == getRegistry().end())
    return "<none>";

  return (*it).second.first;
}

//-----------------------------------------------------------------------------
// Function      : WorkingIntegrationMethod::WorkingIntegrationMethod
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
WorkingIntegrationMethod::WorkingIntegrationMethod(Stats::Stat parent_stat )
  : timeIntegrationMethod_(0),
    jacLimitFlag(false),
    jacLimit(1.0e+17),
    timeIntegratorStat_("Time integrator", parent_stat),
    predictorStat_("Predictor", timeIntegratorStat_),
    completeStepStat_("Successful Step", timeIntegratorStat_),
    rejectStepStat_("Failed Step", timeIntegratorStat_),
    updateCoefStat_("Update Coef", timeIntegratorStat_),
    residualStat_("Load Residual", timeIntegratorStat_),
    jacobianStat_("Load Jacobian", timeIntegratorStat_),
    initializeStat_("Initialize",  timeIntegratorStat_),
    updateLeadStat_("Lead Currents",  timeIntegratorStat_)
{}

//-----------------------------------------------------------------------------
// Function      : WorkingIntegrationMethod::~WorkingIntegrationMethod()
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
WorkingIntegrationMethod::~WorkingIntegrationMethod()
{
  delete timeIntegrationMethod_;
}

//-----------------------------------------------------------------------------
// Function      : WorkingIntegrationMethod::createTimeIntegMethod
// Purpose       : Creates the time integration method class --- assigning a
//                 pointer and the Leading Coefficient value of the method.
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
void WorkingIntegrationMethod::createTimeIntegMethod(
  int                   type,
  const TIAParams &     tia_params,
  StepErrorControl &    step_error_control,
  DataStore &           data_store)
{
  jacLimitFlag = tia_params.jacLimitFlag;
  jacLimit = tia_params.jacLimit;

  delete timeIntegrationMethod_;
  timeIntegrationMethod_ = createTimeIntegrationMethod(type, tia_params, step_error_control, data_store);

  if (!timeIntegrationMethod_)
    Report::DevelFatal0().in("WorkingIntegrationMethod::createTimeIntegMethod") << "Invalid integration method " << type << " specified";

  if (VERBOSE_TIME)
    Xyce::lout() << "  Integration method = " << timeIntegrationMethod_->getName() << std::endl;
}

bool WorkingIntegrationMethod::isTimeIntegrationMethodCreated()
{
  return timeIntegrationMethod_ != 0;
}

double WorkingIntegrationMethod::partialTimeDeriv() const
{
  double pdt = timeIntegrationMethod_->partialTimeDeriv();
  if (jacLimitFlag && pdt > jacLimit)
    pdt = jacLimit;

  return pdt;
}

void WorkingIntegrationMethod::obtainPredictor()
{
//  Stats::TimeBlock x(predictorStat_);

  timeIntegrationMethod_->obtainPredictor();
}

void WorkingIntegrationMethod::obtainPredictorDeriv()
{
//  Stats::TimeBlock x(predictorStat_);
  timeIntegrationMethod_->obtainPredictorDeriv();
}

void WorkingIntegrationMethod::obtainCorrectorDeriv()
{
  timeIntegrationMethod_->obtainCorrectorDeriv();
}

int WorkingIntegrationMethod::getOrder() const
{
  return timeIntegrationMethod_->getOrder();
}

int WorkingIntegrationMethod::getUsedOrder() const
{
  return timeIntegrationMethod_->getUsedOrder();
}

int WorkingIntegrationMethod::getNumberOfSteps() const
{
  return timeIntegrationMethod_->getNumberOfSteps();
}

int WorkingIntegrationMethod::getNscsco() const
{
  return timeIntegrationMethod_->getNscsco();
}

void WorkingIntegrationMethod::getInitialQnorm(TwoLevelError & tle) const
{
  return timeIntegrationMethod_->getInitialQnorm (tle);
}

void WorkingIntegrationMethod::getTwoLevelError(TwoLevelError & tle) const
{
  return timeIntegrationMethod_->getTwoLevelError(tle);
}

void WorkingIntegrationMethod::setTwoLevelTimeInfo()
{
  return timeIntegrationMethod_->setTwoLevelTimeInfo();
}

void WorkingIntegrationMethod::updateCoeffs()
{
 
//  Stats::TimeBlock x(updateCoefStat_);
  return timeIntegrationMethod_->updateCoeffs();
}

void WorkingIntegrationMethod::updateAdjointCoeffs()
{
  return timeIntegrationMethod_->updateAdjointCoeffs();
}

void WorkingIntegrationMethod::rejectStepForHabanero ()
{
  return timeIntegrationMethod_->rejectStepForHabanero();
}

void WorkingIntegrationMethod::initialize(const TIAParams &tia_params)
{
//  Stats::TimeBlock x( initializeStat_ );
  return timeIntegrationMethod_->initialize(tia_params);
}


void WorkingIntegrationMethod::initializeAdjoint (int index)
{
  return timeIntegrationMethod_->initializeAdjoint(index);
}

void WorkingIntegrationMethod::completeStep(const TIAParams &tia_params)
{

//  Stats::TimeBlock x(completeStepStat_);

  return timeIntegrationMethod_->completeStep(tia_params);
}


void WorkingIntegrationMethod::completeAdjointStep(const TIAParams &tia_params)
{

//  Stats::TimeBlock x(completeStepStat_);

  return timeIntegrationMethod_->completeAdjointStep(tia_params);
}


void WorkingIntegrationMethod::rejectStep(const TIAParams &tia_params)
{
//  Stats::TimeBlock x(rejectStepStat_);
  return timeIntegrationMethod_->rejectStep(tia_params);
}

double WorkingIntegrationMethod::computeErrorEstimate() const
{

//  Stats::Stat ErrorStat_("error estimation", timeIntegratorStat_);
//  Stats::TimeBlock x(ErrorStat_);
  return timeIntegrationMethod_->computeErrorEstimate();
}

void WorkingIntegrationMethod::updateStateDeriv ()
{
  return timeIntegrationMethod_->updateStateDeriv ();
}

void WorkingIntegrationMethod::updateLeadCurrent ()
{

//  Stats::TimeBlock x( updateLeadStat_);
  return timeIntegrationMethod_->updateLeadCurrentVec ();
}

void WorkingIntegrationMethod::obtainResidual()
{
//  Stats::TimeBlock x(residualStat_);
  return timeIntegrationMethod_->obtainResidual();
}

void WorkingIntegrationMethod::obtainSensitivityResiduals()
{
  return timeIntegrationMethod_->obtainSensitivityResiduals();
}

void WorkingIntegrationMethod::updateSensitivityHistoryAdjoint()
{
  return timeIntegrationMethod_->updateSensitivityHistoryAdjoint();
}

void WorkingIntegrationMethod::updateSensitivityHistoryAdjoint2()
{
  return timeIntegrationMethod_->updateSensitivityHistoryAdjoint2();
}

void WorkingIntegrationMethod::obtainFunctionDerivativesForTranAdjoint()
{
  return timeIntegrationMethod_->obtainFunctionDerivativesForTranAdjoint();
}

void WorkingIntegrationMethod::obtainSparseFunctionDerivativesForTranAdjoint()
{
  return timeIntegrationMethod_->obtainSparseFunctionDerivativesForTranAdjoint();
}

void WorkingIntegrationMethod::obtainAdjointSensitivityResidual()
{
  return timeIntegrationMethod_->obtainAdjointSensitivityResidual();
}


void WorkingIntegrationMethod::obtainJacobian()
{

//  Stats::TimeBlock x(jacobianStat_);
  return timeIntegrationMethod_->obtainJacobian();
}

void WorkingIntegrationMethod::applyJacobian(const Linear::Vector& input, Linear::Vector& result)
{
//  Stats::TimeBlock x(jacobianStat_);
  return timeIntegrationMethod_->applyJacobian(input, result);
}

bool WorkingIntegrationMethod::printMPDEOutputSolution(
  Analysis::OutputMgrAdapter &  outputManagerAdapter,
  const double                  time,
  Linear::Vector *              solnVecPtr,
  const std::vector<double> &   fastTimes )
{
  return timeIntegrationMethod_->printMPDEOutputSolution(
    outputManagerAdapter, time, solnVecPtr, fastTimes );
}

bool WorkingIntegrationMethod::printWaMPDEOutputSolution(
  Analysis::OutputMgrAdapter &  outputManagerAdapter,
  const double                  time,
  Linear::Vector *              solnVecPtr,
  const std::vector<double> &   fastTimes,
  const int                     phiGID )
{
  return timeIntegrationMethod_->printWaMPDEOutputSolution(
    outputManagerAdapter, time, solnVecPtr, fastTimes, phiGID );
}

bool WorkingIntegrationMethod::printOutputSolution(
  Analysis::OutputMgrAdapter &  outputManagerAdapter,
  const TIAParams &             tia_params,
  const double                  time,
  Linear::Vector *              solnVecPtr,
  const bool                    doNotInterpolate,
  const std::vector<double> &   outputInterpolationTimes,
  bool                          skipPrintLineOutput )
{
//  Stats::TimeBlock x( othersStat_);
  return timeIntegrationMethod_->printOutputSolution(
    outputManagerAdapter, tia_params, time, solnVecPtr, doNotInterpolate, outputInterpolationTimes, skipPrintLineOutput) ;
}

bool WorkingIntegrationMethod::saveOutputSolution(
  Parallel::Machine                     comm,
  IO::InitialConditionsManager &        initial_conditions_manager,
  const NodeNameMap &                   node_name_map,
  const TIAParams &                     tia_params,
  Linear::Vector *                      solnVecPtr,
  const double                          saveTime,
  const bool                            doNotInterpolate)
{
  return timeIntegrationMethod_->saveOutputSolution(comm, initial_conditions_manager, node_name_map, tia_params, solnVecPtr, saveTime, doNotInterpolate);
}

void WorkingIntegrationMethod::stepLinearCombo()
{
  if (timeIntegrationMethod_)
  {
    timeIntegrationMethod_->stepLinearCombo();
  }
}

bool WorkingIntegrationMethod::getSolnVarData( const int & gid, std::vector<double> & varData ) 
{ 
  if (timeIntegrationMethod_)
  {
    return timeIntegrationMethod_->getSolnVarData( gid, varData );
  }
  return false;
}

bool WorkingIntegrationMethod::setSolnVarData( const int & gid, const std::vector<double> & varData ) 
{ 
  if (timeIntegrationMethod_)
  {
    return timeIntegrationMethod_->setSolnVarData( gid, varData );
  }
  return false;
}

bool WorkingIntegrationMethod::getStateVarData( const int & gid, std::vector<double> & varData ) 
{ 
  if (timeIntegrationMethod_)
  {
    return timeIntegrationMethod_->getStateVarData( gid, varData );
  }
  return false;
}

bool WorkingIntegrationMethod::setStateVarData( const int & gid, const std::vector<double> & varData ) 
{ 
  if (timeIntegrationMethod_)
  {
    return timeIntegrationMethod_->setStateVarData( gid, varData );
  }
  return false;
}

bool WorkingIntegrationMethod::getStoreVarData( const int & gid, std::vector<double> & varData ) 
{ 
  if (timeIntegrationMethod_)
  {
    return timeIntegrationMethod_->getStoreVarData( gid, varData );
  }
  return false;
}

bool WorkingIntegrationMethod::setStoreVarData( const int & gid, const std::vector<double> & varData ) 
{ 
  if (timeIntegrationMethod_)
  {
    return timeIntegrationMethod_->setStoreVarData( gid, varData );
  }
  return false;
}

} // namespace TimeIntg
} // namespace Xyce

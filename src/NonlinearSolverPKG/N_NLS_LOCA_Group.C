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

//-------------------------------------------------------------------------
//
// Purpose        : Interface to Xyce for LOCA continuation routines.
//
// Special Notes  :
//
// Creator        : Roger Pawlowski, SNL 9233
//
// Creation Date  : 02/17/03
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_OutputMgr.h>
#include <N_LAS_Builder.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_Vector.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_NLS_LOCA_Group.h>
#include <N_NLS_NOX_AugmentLinSys.h>
#include <N_NLS_NOX_SharedSystem.h>
#include <N_PDS_Comm.h>
#include <N_PDS_ParMap.h>
#include <N_UTL_FeatureTest.h>

// ----------   NOX Includes   ----------
#include <Teuchos_ParameterList.hpp>
#include <NOX_Abstract_Vector.H>
#include <LOCA_Parameter_Vector.H>

// using namespace Xyce::Nonlinear::N_NLS_NOX;

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_LOCA {

//-----------------------------------------------------------------------------
// Function      : Group::Group
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Group::Group(Teuchos::RCP<LOCA::GlobalData> gd,
             N_NLS_NOX::SharedSystem& s, Loader::NonlinearEquationLoader& l,
             IO::OutputMgr& o, Analysis::AnalysisManager & t) :
  N_NLS_NOX::Group(s),
  LOCA::Abstract::Group(gd),
  loader(l),
  outputMgr(o),
  anaInt(t),
  derivUtils(gd),
  tmpVectorPtr(0),
  scalingVecPtr(0),
  useAugmentLinSys_(false),
  nonContinuationSolve_(true)
{
  tmpVectorPtr = sharedSystemPtr_->getLasSystem()->builder().createVector();
}

//-----------------------------------------------------------------------------
// Function      : Group::Group
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Group::Group(const Group& source, NOX::CopyType type) :
  N_NLS_NOX::Group(source, type),
  LOCA::Abstract::Group(source, type),
  loader(source.loader),
  outputMgr(source.outputMgr),
  anaInt(source.anaInt),
  params(source.params),
  derivUtils(source.derivUtils),
  tmpVectorPtr(0),
  scalingVecPtr(source.scalingVecPtr),
  useAugmentLinSys_(source.useAugmentLinSys_),
  augmentLSStrategy_(source.augmentLSStrategy_),
  nonContinuationSolve_(source.nonContinuationSolve_)
{
  tmpVectorPtr = sharedSystemPtr_->getLasSystem()->builder().createVector();
}

//-----------------------------------------------------------------------------
// Function      : Group::~Group
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Group::~Group()
{
  delete tmpVectorPtr;
}

//-----------------------------------------------------------------------------
// Function      : Group::operator=
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
NOX::Abstract::Group&
Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&>(source));
}

//-----------------------------------------------------------------------------
// Function      : Group::operator=
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
N_NLS_NOX::Group&
Group::operator=(const N_NLS_NOX::Group& source)
{
  return operator=(dynamic_cast<const Group&>(source));
}

//-----------------------------------------------------------------------------
// Function      : Group::operator=
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
LOCA::Abstract::Group&
Group::operator=(const LOCA::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&>(source));
}

//-----------------------------------------------------------------------------
// Function      : Group::operator=
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Group&
Group::operator=(const Group& source)
{
  N_NLS_NOX::Group::operator=(source);
  params = source.params;
  derivUtils = source.derivUtils;
  if (source.scalingVecPtr != 0)
    scalingVecPtr = source.scalingVecPtr;
  useAugmentLinSys_ = source.useAugmentLinSys_;
  augmentLSStrategy_ = source.augmentLSStrategy_;
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : Group::copy
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void Group::copy(const NOX::Abstract::Group &source)
{
  *this = source;
}

//-----------------------------------------------------------------------------
// Function      : Group::clone
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Teuchos::RCP<NOX::Abstract::Group> Group::
clone(NOX::CopyType type) const
{
  Teuchos::RCP<Group> ptr =
    Teuchos::rcp(new Group(*this, type));
  return ptr;
}

//-----------------------------------------------------------------------------
// Function      : Group::computeF()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType Group::computeF()
{
  // Set continuation parameters, if neccessary
  // Note:  If the tranOP solve was a continuation, but the
  // transient is traditional Newton, make SURE that the
  // setParam calls are not happening for transient!
  if (!nonContinuationSolve_)
  {
    for (int i = 0; i < params.length(); ++i) 
    {
      std::string label = params.getLabel(i);
      loader.setParam(label, params.getValue(i));

      if (label == "GSTEPPING" && useAugmentLinSys_)
        augmentLSStrategy_->setProgressVariable(params.getValue(i));

    }
  }

  NOX::Abstract::Group::ReturnType status = N_NLS_NOX::Group::computeF();

  if (useAugmentLinSys_)
    augmentLSStrategy_->augmentResidual(xVec_.getNativeVectorPtr(),
					fVec_.getNativeVectorPtr());

  return status;
}

//-----------------------------------------------------------------------------
// Function      : Group::computeJacobian()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType Group::computeJacobian()
{
  // Set continuation parameters, if neccessary
  // Note:  If the tranOP solve was a continuation, but the
  // transient is traditional Newton, make SURE that the
  // setParam calls are not happening for transient!
  if (!nonContinuationSolve_)
  {
    for (int i = 0; i < params.length(); ++i) {
      std::string label = params.getLabel(i);
      loader.setParam(label, params.getValue(i));

      if (label == "GSTEPPING" && useAugmentLinSys_)
        augmentLSStrategy_->setProgressVariable(params.getValue(i));

    }
  }

  NOX::Abstract::Group::ReturnType status =
    N_NLS_NOX::Group::computeJacobian();

  // Augment jacobian for pseudo transient if enabled
  if (useAugmentLinSys_) {
    Linear::Matrix& jacobian =
      const_cast<Linear::Matrix&>(sharedSystemPtr_->getJacobian());
    augmentLSStrategy_->augmentJacobian(&jacobian);
  }

  if (DEBUG_NONLINEAR)
  {
    Linear::Matrix& jacobian = const_cast<Linear::Matrix&>(sharedSystemPtr_->getJacobian());
    sharedSystemPtr_->debugOutput1( jacobian, (*(fVec_.getNativeVectorPtr())));
  }

  return status;
}

//-----------------------------------------------------------------------------
// Function      : Group::computeDfDpMulti	
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType Group::computeDfDpMulti	
  (const std::vector< int > & paramIDs, 
   NOX::Abstract::MultiVector & dfdp, 
   bool isValidF)
{
  // Do it the old way (LOCA performs the finite difference DfDp):
  //if (true)
  if (false)
  {
    // obtain voltage limiting status:
    bool voltageLimterStatus =  loader.getVoltageLimiterStatus();

    // turn off voltage limiting:
    loader.setVoltageLimiterStatus(false);

    LOCA::Abstract::Group::computeDfDpMulti (paramIDs, dfdp, isValidF);

    // turn on voltage limiting (assuming it was on to begin with):
    loader.setVoltageLimiterStatus(voltageLimterStatus);
  }
  else 
  {
    // Do it the new way (Xyce provides DfDp)
    bool tmp = sharedSystemPtr_->computeDfDpMulti (paramIDs, dfdp, isValidF);
  }

  return NOX::Abstract::Group::Ok;
}

//-----------------------------------------------------------------------------
// Function      : Group::setParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void Group::setParams(const LOCA::ParameterVector& p)
{
  N_NLS_NOX::Group::resetIsValid_();
  params = p;
  return;
}

//-----------------------------------------------------------------------------
// Function      : Group::getParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
const LOCA::ParameterVector& Group::getParams() const
{
  return params;
}

//-----------------------------------------------------------------------------
// Function      : Group::setParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void Group::setParam(int paramID, double value)
{
  N_NLS_NOX::Group::resetIsValid_();
  params.setValue(paramID, value);
  return;
}

//-----------------------------------------------------------------------------
// Function      : Group::getParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
double Group::getParam(int paramID) const
{
  return params.getValue(paramID);
}

//-----------------------------------------------------------------------------
// Function      : Group::setParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void Group::setParam(std::string paramID, double value)
{
  N_NLS_NOX::Group::resetIsValid_();
  params.setValue(paramID, value);
  return;
}

//-----------------------------------------------------------------------------
// Function      : Group::getParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
double Group::getParam(std::string paramID) const
{
  return params.getValue(paramID);
}

//-----------------------------------------------------------------------------
// Function      : Group::setScaleVec
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void Group::setScaleVec(const NOX::Abstract::Vector& s)
{
  scalingVecPtr = &s;
}

//-----------------------------------------------------------------------------
// Function      : Group::getScaleVec()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
const NOX::Abstract::Vector& Group::getScaleVec() const
{
  if (scalingVecPtr == 0)
  {
    Xyce::Report::DevelFatal().in("Group::getScaleVec")
      << "scaling vector not set!";
  }

  return (*scalingVecPtr);
}

//-----------------------------------------------------------------------------
// Function      : Group::augmentJacobianForHomotopy
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType
Group::augmentJacobianForHomotopy(double conParamValue)
{

  Linear::Matrix& jacobian =
    const_cast<Linear::Matrix&>(sharedSystemPtr_->getJacobian());

  //jacobian.printPetraObject();
  jacobian.scale(conParamValue);
  jacobian.getDiagonal(*tmpVectorPtr);
  (*tmpVectorPtr).addScalar(1.0 - conParamValue);
  jacobian.replaceDiagonal(*tmpVectorPtr);
  //jacobian.printPetraObject();

  return NOX::Abstract::Group::Ok;
}

//-----------------------------------------------------------------------------
// Function      : Group::printSolution
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void Group::printSolution (const double conParam) const
{
  // ERK:  This is a KLUDGE!  I put it here because the "printSolution"
  // functions were the only place in this class that were called
  // only after a continuation step was successful.
  //
  // Note: homotopy output is now called from the time integrator,
  // so it is this function still results in the the solution being
  // printed.
  //
  // Note:  This will be replaced by "stepSuccess" (see below), once
  // stepSuccess completely works.
  anaInt.completeHomotopyStep(loader, params.getNamesVector(), params.getValuesVector(), xVec_.getNativeVectorPtr());

  return;
}

//-----------------------------------------------------------------------------
// Function      : Group::printSolution
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void Group::printSolution (const NOX::Abstract::Vector &x,
				       const double conParam) const
{
  // const N_NLS_NOX::Vector& noxVector =
  //   dynamic_cast<const N_NLS_NOX::Vector&>(x);

  // ERK: The function call below (completeHomotopyStep) handles both PRINT 
  // output , and also other data management that is necessary at the end of a 
  // successful homotopy step. I put it here because (at the time) the 
  // "printSolution" functions were the only place in this class that were called
  // only after a continuation step was successful.  That is no longer
  // true (see postPorcessContinuationStep, below), but there isn't
  // a compelling reason to change it at this point.
  anaInt.completeHomotopyStep(loader, params.getNamesVector(), params.getValuesVector(), xVec_.getNativeVectorPtr());

  return;
}

//-----------------------------------------------------------------------------
// Function      : Group::preProcessContinuationStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/20/2015
//-----------------------------------------------------------------------------
void Group::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
}

//-----------------------------------------------------------------------------
// Function      : Group::postProcessContinuationStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/20/2015
//-----------------------------------------------------------------------------
void Group::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  if (stepStatus == 0)
  {
    anaInt.failHomotopyStep(loader);
  }
}


//-----------------------------------------------------------------------------
// Function      : Group::setAugmentLinearSystem
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void Group::setAugmentLinearSystem(bool enable_value,
		const Teuchos::RCP<N_NLS_NOX::AugmentLinSys>& ls)
{
  useAugmentLinSys_ = enable_value;
  augmentLSStrategy_ = ls;
}

//-----------------------------------------------------------------------------
// Function      : Group::setNonContinuationFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void Group::setNonContinuationFlag (bool value)
{
  nonContinuationSolve_ = value;
}

//-----------------------------------------------------------------------------
// Function      : Group::getNonContinuationFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Group::getNonContinuationFlag ()
{
  return nonContinuationSolve_;
}

}}} // namespace N_NLS_LOCA

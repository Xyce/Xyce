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
// Purpose       : This file contains the functions which define the
//		   time integration methods classes.
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

#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_ERH_ErrorMgr.h>
#include <N_IO_InitialConditions.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_Vector.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_FilteredMultiVector.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_NoTimeIntegration.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TIAParams.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_TIA_TwoLevelError.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace TimeIntg {

const char *
NoTimeIntegration::name = "None";

//-----------------------------------------------------------------------------
// Function      : NoTimeIntegration::factory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/08/01
//-----------------------------------------------------------------------------
TimeIntegrationMethod *
NoTimeIntegration::factory(
    const TIAParams &   tia_params,
    StepErrorControl &  step_error_control,
    DataStore &         data_store)
{
  return new NoTimeIntegration(tia_params, step_error_control, data_store);
}

//-----------------------------------------------------------------------------
// Function      : NoTimeIntegration::NoTimeIntegration
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
NoTimeIntegration::NoTimeIntegration(
  const TIAParams & tia_params,
  StepErrorControl & secTmp,
  DataStore & dsTmp)
  : TimeIntegrationMethod(),
    ds(dsTmp),
    sec(secTmp),
    leadingCoeff(1.0)
{
  leadingCoeff = 1.0;
  alphas = -1.0;
  return;
}

//-----------------------------------------------------------------------------
// Function      : NoTimeIntegration::~NoTimeIntegration
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
NoTimeIntegration::~NoTimeIntegration() {}


void
NoTimeIntegration::completeStep(
  const TIAParams &     tia_params)
{
  sec.lastTime    = sec.currentTime;
  sec.currentTime = sec.nextTime;
}

void
NoTimeIntegration::rejectStep(
  const TIAParams &     tia_params)
{}

//-----------------------------------------------------------------------------
// Function      : NoTimeIntegration::obtainCorrectorDeriv
// Purpose       : Evaluate the Corrector Derivative Formula
// Special Notes : For "no integration" the derivatives should always be zero.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/10/01
//-----------------------------------------------------------------------------
void NoTimeIntegration::obtainCorrectorDeriv()
{
  ds.nextStateDerivPtr->putScalar(0.0);
}


//-----------------------------------------------------------------------------
// Function      : NoTimeIntegration::obtainResidual
//
// Purpose       : This function returns the residual for the steady state
//                 case.
//
// Special Notes : For "no integration" the derivatives should always be zero,
//                 so don't add in the dqdt term.
//
//                 This function is only called in the new-DAE case.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/09/04
//-----------------------------------------------------------------------------
void NoTimeIntegration::obtainResidual ()
{
  ds.RHSVectorPtr->linearCombo(+1.0,*ds.daeFVectorPtr,-1.0,*ds.daeBVectorPtr);

  // since the nonlinear solver is expecting a -f, scale by -1.0:
  ds.RHSVectorPtr->scale(-1.0);

  // if voltage limiting is on, add it in:
  if (ds.limiterFlag)
  {
    (ds.RHSVectorPtr)->daxpy(
      *(ds.RHSVectorPtr), +1.0, *(ds.dFdxdVpVectorPtr));
  }
}

//-----------------------------------------------------------------------------
// Function      : NoTimeIntegration::obtainSensitivityResiduals
//
// Purpose       : This function returns the residual for the steady state
//                 case.
//
// Special Notes : For "no integration" the derivatives should always be zero,
//                 so don't add in the dqdt term.
//
//                 This function is only called in the new-DAE case.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void NoTimeIntegration::obtainSensitivityResiduals ()
{
  ds.sensRHSPtrVector->linearCombo(+1.0, *(ds.nextDfdpPtrVector),
                                   -1.0, *(ds.nextDbdpPtrVector));
 
  // since the nonlinear solver is expecting a -f, scale by -1.0:
  ds.sensRHSPtrVector->scale(-1.0);

#ifdef DEBUG_SENS
    Xyce::dout() << "NoTimeIntegration: obtainSensitivityResiduals: " << std::endl; 
    ds.sensRHSPtrVector->printPetraObject(Xyce::dout());
#endif

}

//-----------------------------------------------------------------------------
// Function      : NoTimeIntegration::obtainFunctionDerivativesForTranAdjoint
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/13/2016
//-----------------------------------------------------------------------------
void NoTimeIntegration::obtainFunctionDerivativesForTranAdjoint ()
{
  ds.sensRHSPtrVector->linearCombo(+1.0, *(ds.nextDfdpPtrVector),
                                   -1.0, *(ds.nextDbdpPtrVector));

  // since the nonlinear solver is expecting a -f, scale by -1.0:
  ds.sensRHSPtrVector->scale(-1.0);
}

//-----------------------------------------------------------------------------
// Function      : NoTimeIntegration::obtainSparseFunctionDerivativesForTranAdjoint
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/9/2017
//-----------------------------------------------------------------------------
void NoTimeIntegration::obtainSparseFunctionDerivativesForTranAdjoint ()
{
  ds.sensRHSPtrVector->linearCombo(+1.0, *(ds.nextDfdpPtrVector),
                                   -1.0, *(ds.nextDbdpPtrVector));

  // since the nonlinear solver is expecting a -f, scale by -1.0:
  ds.sensRHSPtrVector->scale(-1.0);

  // Collect sparse information from sensRHSPtrVector.
  ds.sparseSensRHSMV->replaceValues( ds.masterIndexVector, *ds.sensRHSPtrVector );  
    
}

//-----------------------------------------------------------------------------
// Function      : NoTimeIntegration::obtainAdjointSensitivityResidual
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/23/2015
//-----------------------------------------------------------------------------
void NoTimeIntegration::obtainAdjointSensitivityResidual ()
{

}

//-----------------------------------------------------------------------------
// Function      : NoTimeIntegration::obtainJacobian
//
// Purpose       : Returns the full Jacobian matrix for the steady state
//                 case.  
//
// Special Notes : For "no integration" the derivatives should always be zero.
//                 However, to prevent singular matrices, the dQdt term is
//                 added in anyway, with an assumed very large time step.
//
//                 There may be a better way to handle this - the assumed
//                 large time step leads to a very small alpha/dt, which
//                 could (maybe?) have an adverse effect on the matrix
//                 conditioning.
//
//                 The singular matrix problem will happen for the case of
//                 capacitors in serial.  Any voltage node which is *only*
//                 connected to capacitors will not really be part of the
//                 system of equations.  Ideally nodes like this would just
//                 be removed from the system altogether.  No current is
//                 passing through them, and they are really just floating
//                 in space.
//
//                 For the time being, however, what we do is put some
//                 bogus C*alpha/dt terms into the Jacobian.  If dt is
//                 large, then the C*alpha/dt is very small, and doesn't
//                 significantly change any Jacobian entries, except for
//                 zero entries.
//
//                 This function is only called in the new-DAE case.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/09/04
//-----------------------------------------------------------------------------
void NoTimeIntegration::obtainJacobian ()
{
  Linear::Matrix & dQdx = *(ds.dQdxMatrixPtr);
  Linear::Matrix & dFdx = *(ds.dFdxMatrixPtr);
  Linear::Matrix & Jac = *(ds.JMatrixPtr);

  Jac.linearCombo( 1.0e-20, dQdx, 1.0, dFdx );
}

//-----------------------------------------------------------------------------
// Function      : NoTimeIntegration::applyJacobian
//
// Purpose       : Applies the Jacobian operator for the steady state
//                 case.  
//
// Special Notes : 
//                 This function is only called in the new-DAE HB (matrix-free) case.
//
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
void NoTimeIntegration::applyJacobian (const Linear::Vector& input, Linear::Vector& result)
{
  Linear::Vector & dQdxV = *(ds.dQdxVecVectorPtr);
  Linear::Vector & dFdxV = *(ds.dFdxVecVectorPtr);
  result.linearCombo( 1.0e-20, dQdxV, 1.0, dFdxV );
}

//-----------------------------------------------------------------------------
// Function      : NoTimeIntegration::getInitialQnorm
// Purpose       : Needed by 2-level solves.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/18/07
//-----------------------------------------------------------------------------
void NoTimeIntegration::getInitialQnorm(TwoLevelError & tle) const
{
  tle.q1HistorySum = 0.0;
}

//-----------------------------------------------------------------------------
// Function      : NoTimeIntegration::setupTwoLevelError
// Purpose       : Needed by 2-level solves.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
void NoTimeIntegration::getTwoLevelError(TwoLevelError & tle) const
{
  tle.xErrorSum    = 0.0;
  tle.qErrorSum    = 0.0;
  tle.xErrorSum_m1 = 0.0;
  tle.innerSize    = ds.newtonCorrectionPtr->globalLength();
}

//-----------------------------------------------------------------------------
bool NoTimeIntegration::printOutputSolution(
  Analysis::OutputMgrAdapter &  outputManagerAdapter, 
  const TIAParams &             tia_params,
  const double                  time,
  Linear::Vector *                solnVecPtr,
  const bool                    doNotInterpolate,
  const std::vector<double> &   outputInterpolationTimes,
  bool                          skipPrintLineOutput )
{
  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << "NoTimeIntegration::printOutputSolution" << std::endl;

    Xyce::dout() << "tmp lead current vector:" <<std::endl;
    ds.tmpLeadCurrentVectorPtr->printPetraObject(Xyce::dout());

    Xyce::dout() << "tmp lead deltaV vector:" <<std::endl;
    ds.tmpLeadDeltaVPtr->printPetraObject(Xyce::dout());

  }

  double dt = 0.0;
  outputManagerAdapter.tranOutput( 
      time, dt, sec.finalTime,
      *solnVecPtr, 
      *ds.currStatePtr, *ds.currStorePtr, 
      *ds.tmpLeadCurrentVectorPtr, *ds.tmpLeadDeltaVPtr, *ds.tmpLeadCurrentQDerivVectorPtr,
                                   ds.objectiveVec_, 
                                   ds.dOdpVec_, ds.dOdpAdjVec_, 
                                   ds.scaled_dOdpVec_, ds.scaled_dOdpAdjVec_, 
                                   skipPrintLineOutput);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : TimeIntegrationMethod::saveOutputSolution
// Purpose       : 
// Special Notes : For the old method functions (old-DAE/ODE) no interpolation
//                 is possible, so this function calls directly through to the
//                 output manager.
// Scope         : public
// Creator       : Eric Keiter, SNL, 1437
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
bool NoTimeIntegration::saveOutputSolution(
  Parallel::Machine                     comm,
  IO::InitialConditionsManager &        initial_conditions_manager,
  const NodeNameMap &                   node_name_map,
  const TIAParams &                     tia_params,
  Linear::Vector *                      solnVecPtr,
  const double                          saveTime,
  const bool                            doNotInterpolate) 
{
  // outputManagerAdapter.outputDCOP( *(solnVecPtr) );
  initial_conditions_manager.outputDCOP(comm, node_name_map, *solnVecPtr);

  return true;
}

} // namespace TimeIntg
} // namespace Xyce

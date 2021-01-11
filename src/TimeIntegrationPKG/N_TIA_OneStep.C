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
//		   Trap method, order 1-2, class.
//
// Special Notes :
//
// Creator       : Ting Mei
//
// Creation Date : 2/16/04
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_TIA_OneStep.h>

#include <N_ANP_OutputMgrAdapter.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_InitialConditions.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_Vector.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_FilteredMultiVector.h>
#include <N_PDS_Comm.h>
#include <N_PDS_Manager.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TIAParams.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_MachDepParams.h>

using std::min;
using std::max;
using std::abs;

namespace Xyce {
namespace TimeIntg {

const char *
OneStep::name = "Onestep: Trapezoidal";

//-----------------------------------------------------------------------------
// Function      : OneStep::factory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :  10/31/07
//-----------------------------------------------------------------------------
TimeIntegrationMethod *
OneStep::factory(
    const TIAParams &   tia_params,
    StepErrorControl &  step_error_control,
    DataStore &         data_store)
{
  return new OneStep(tia_params, step_error_control, data_store);
}


//-----------------------------------------------------------------------------
// Function      : OneStep::OneStep
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------
OneStep::OneStep(
  const TIAParams & tia_params,
  StepErrorControl & secTmp,
  DataStore & dsTmp)
  : TimeIntegrationMethod(),
    ds(dsTmp),
    sec(secTmp),
    leadingCoeff(1.0)
{
  leadingCoeff = 1;
  sec.maxOrder_=(std::min(2,tia_params.maxOrder));
  sec.minOrder_=(std::max(1,tia_params.minOrder));

  if (sec.minOrder_ > sec.maxOrder_)
  {
    sec.minOrder_ = sec.maxOrder_;
  }
  //  sec.maxOrder_ = 2;
  timept_ = -1.0;
  timeStepForHistory2_ = 0.0;

  sec.currentOrder_ = (std::min(sec.currentOrder_, sec.maxOrder_) );
}

//-----------------------------------------------------------------------------
// Function      : OneStep::obtainPredictor
// Purpose       : Calculate predictor
// Special Notes : stored in ds.xn0Ptr,qn0Ptr
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------
void OneStep::obtainPredictor()
{
  // evaluate predictor
  *ds.xn0Ptr = *(ds.xHistory[0]);
  *ds.qn0Ptr = *(ds.qHistory[0]);

  for (int i=1;i<=sec.currentOrder_;++i)
  {
    ds.xn0Ptr->update(sec.beta_[i],*(ds.xHistory[i]));
  }

  if (DEBUG_TIME && isActive(Diag::TIME_PREDICTOR))
  {
    Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() <<
      "  OneStep::obtainPredictor" << std::endl;
    Xyce::dout() << "\n currentOrder = " << sec.currentOrder_ << std::endl;
    Xyce::dout() << "\n sec.nscsco_: " << sec.nscsco_ << std::endl;
    for (int i=0; i<=sec.currentOrder_ ; ++i)
      Xyce::dout() << "\n sec.beta_[" << i << "] = " << sec.beta_[i] << "\n" << std::endl;
    for (int i=0; i<=sec.currentOrder_ ; ++i)
    {
      Xyce::dout() << "\n xHistory["<< i << "]: \n" << std::endl;
      (ds.xHistory[i])->print(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    for (int i=0; i<=sec.currentOrder_ ; ++i)
    {
      Xyce::dout() << "\n qHistory["<< i << "]: \n" << std::endl;
      (ds.qHistory[i])->print(Xyce::dout());
      Xyce::dout() << std::endl;
    }

    Xyce::dout() << "\n sHistory["<< 0 << "]: \n" << std::endl;
    (ds.sHistory[0])->print(Xyce::dout());
    Xyce::dout() << std::endl;

    for (int i=0; i<sec.currentOrder_ ; ++i)
    {
      Xyce::dout() << "\n stoHistory["<< i << "]: \n" << std::endl;
      (ds.stoHistory[i])->print(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    Xyce::dout() << "\n xn0: \n" << std::endl;
    ds.xn0Ptr->print(Xyce::dout());
    Xyce::dout() << std::endl;
    Xyce::dout() << "\n qn0: \n" << std::endl;
    ds.qn0Ptr->print(Xyce::dout());
    Xyce::dout() << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  // copy the prediction into the next solution:
  *(ds.nextSolutionPtr) = *(ds.xn0Ptr);

  obtainSensitivityPredictors();

  return;
}

//-----------------------------------------------------------------------------
// Function      : OneStep::obtainSensitivityPredictors
// Purpose       : Calculate predictor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void OneStep::obtainSensitivityPredictors()
{

  return;
}

//-----------------------------------------------------------------------------
// Function      : OneStep::obtainResidual
// Purpose       : Calculate Residual
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------
void OneStep::obtainResidual()
{
  // output: ds.RHSVectorPtr
  // Note:  ds.nextSolutionPtr is used to get Q,F,B in Analysis::AnalysisManager::loadRHS.
  ds.RHSVectorPtr->linearCombo(1.0,*ds.daeQVectorPtr,-1.0,*(ds.qHistory[0]));

  if (DEBUG_TIME && isActive(Diag::TIME_RESIDUAL))
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() <<
      "  OneStep::obtainResidual" << std::endl;
    Xyce::dout() << "\n t = " << sec.nextTime << "\n" << std::endl;
    Xyce::dout() << "\n solution: \n" << std::endl;
    ds.nextSolutionPtr->print(Xyce::dout());
    Xyce::dout() << "\n daeQVector: \n" << std::endl;
    ds.daeQVectorPtr->print(Xyce::dout());
    Xyce::dout() << "\n qn0: \n" << std::endl;
    ds.qn0Ptr->print(Xyce::dout());
    Xyce::dout() << "\n sec.alphas_/hn: " << sec.alphas_/sec.currentTimeStep << "\n" << std::endl;
    Xyce::dout() << "\n daeFVector: \n" << std::endl;
    ds.daeFVectorPtr->print(Xyce::dout());

    Xyce::dout() << "\n dQdt-vector: \n" << std::endl;
    ds.RHSVectorPtr->print(Xyce::dout());
    Xyce::dout() << std::endl;
  }

  if (sec.currentOrder_  == 2)
  {
    ds.RHSVectorPtr->update(1.0/2.0,*ds.daeFVectorPtr,-1.0/2.0,*ds.daeBVectorPtr,1.0/sec.currentTimeStep);
    ds.RHSVectorPtr->update(+1.0/2.0,*ds.qHistory[2]);
  }
  else
  {
    ds.RHSVectorPtr->update(+1.0,*ds.daeFVectorPtr,-1.0,*ds.daeBVectorPtr,1.0/sec.currentTimeStep);
  }
  // since the nonlinear solver is expecting a -f, scale by -1.0:
  ds.RHSVectorPtr->scale(-1.0);

  // if voltage limiting is on, add it in:
  if (ds.limiterFlag)
  {
    (ds.dQdxdVpVectorPtr)->scale( -sec.alphas_/sec.currentTimeStep );

    (ds.RHSVectorPtr)->update(+1.0, *(ds.dQdxdVpVectorPtr));

    double fscalar(1.0);

    if (sec.currentOrder_  == 2)
      fscalar =1.0/2.0;

    (ds.RHSVectorPtr)->update(fscalar, *(ds.dFdxdVpVectorPtr));
  }

  if (DEBUG_TIME && isActive(Diag::TIME_RESIDUAL))
  {
    Xyce::dout() << "\n Residual-vector: \n" << std::endl;
    ds.RHSVectorPtr->print(Xyce::dout());
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : OneStep::obtainSensitivityResiduals
// Purpose       : Calculate Sensitivity Residual
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/02/2014
//-----------------------------------------------------------------------------
void OneStep::obtainSensitivityResiduals()
{
  if (sec.currentOrder_ == 2)
  {
    ds.nextDqdpDerivPtrVector->linearCombo(1.0,*ds.nextDqdpPtrVector,-1.0, *(ds.dqdpHistory[0]));
    ds.nextDqdpDerivPtrVector->scale(1.0/sec.currentTimeStep);

    ds.sensRHSPtrVector->linearCombo(1.0,*ds.nextDqdpDerivPtrVector,+0.5,*ds.nextDfdpPtrVector);
    ds.sensRHSPtrVector->update(+0.5,*ds.dfdpHistory[0]);

    ds.sensRHSPtrVector->update(-0.5,*ds.nextDbdpPtrVector);
    ds.sensRHSPtrVector->update(-0.5,*ds.dbdpHistory[0]);
  }
  else
  {
    ds.nextDqdpDerivPtrVector->linearCombo(1.0,*ds.nextDqdpPtrVector,-1.0, *(ds.dqdpHistory[0]));
    ds.nextDqdpDerivPtrVector->scale(1.0/sec.currentTimeStep);

    ds.sensRHSPtrVector->linearCombo(1.0,*ds.nextDqdpDerivPtrVector,+1.0,*ds.nextDfdpPtrVector);
    ds.sensRHSPtrVector->update(-1.0,*ds.nextDbdpPtrVector);
  }

  // since the nonlinear solver is expecting a -dFdp, scale by -1.0:
  ds.sensRHSPtrVector->scale(-1.0);

  // correction terms
  double qscalar(1.0/sec.currentTimeStep);
  ds.sensRHSPtrVector->update(qscalar, *ds.currDQdxDXdpPtrVector);

  // second order "correction" term
  if (sec.currentOrder_ == 2)
  {
    ds.sensRHSPtrVector->update(-0.5, *ds.currDFdxDXdpPtrVector);
  }

#ifdef DEBUG_SENS
  Xyce::dout() << "OneStep: obtainSensitivityResiduals: RHS Vector : " << std::endl;;
  ds.sensRHSPtrVector->print(Xyce::dout());
#endif
}

//-----------------------------------------------------------------------------
// Function      : OneStep::obtainFunctionDerivativesForTranAdjoint
// Purpose       : 
//
// Special Notes : This is similar to the direct-sensitivity method, 
//                 OneStep::obtainSensitivityResiduals .  The main difference is
//                 that this function doesn't include the so-called cross terms
//                 (or chain rule terms).  It only loads the function (device) 
//                 parameter derivatives.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/13/2016
//-----------------------------------------------------------------------------
void OneStep::obtainFunctionDerivativesForTranAdjoint ()
{
  if (sec.currentOrder_ == 2)
  {
    ds.nextDqdpDerivPtrVector->linearCombo(1.0,*ds.nextDqdpPtrVector,-1.0, *(ds.dqdpHistory[0]));
    ds.nextDqdpDerivPtrVector->scale(1.0/sec.currentTimeStep);

    ds.sensRHSPtrVector->linearCombo(1.0,*ds.nextDqdpDerivPtrVector,+0.5,*ds.nextDfdpPtrVector);
    ds.sensRHSPtrVector->update(+0.5,*ds.dfdpHistory[0]);

    ds.sensRHSPtrVector->update(-0.5,*ds.nextDbdpPtrVector);
    ds.sensRHSPtrVector->update(-0.5,*ds.dbdpHistory[0]);
  }
  else
  {
    ds.nextDqdpDerivPtrVector->linearCombo(1.0,*ds.nextDqdpPtrVector,-1.0, *(ds.dqdpHistory[0]));
    ds.nextDqdpDerivPtrVector->scale(1.0/sec.currentTimeStep);

    ds.sensRHSPtrVector->linearCombo(1.0,*ds.nextDqdpDerivPtrVector,+1.0,*ds.nextDfdpPtrVector);
    ds.sensRHSPtrVector->update(-1.0,*ds.nextDbdpPtrVector);
  }

  // since the nonlinear solver is expecting a -dFdp, scale by -1.0:
  ds.sensRHSPtrVector->scale(-1.0);
}

//-----------------------------------------------------------------------------
// Function      : OneStep::obtainSparseFunctionDerivativesForTranAdjoint
// Purpose       : 
//
// Special Notes : This is similar to the direct-sensitivity method, 
//                 OneStep::obtainSensitivityResiduals .  The main difference is
//                 that this function doesn't include the so-called cross terms
//                 (or chain rule terms).  It only loads the function (device) 
//                 parameter derivatives.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/9/2017
//-----------------------------------------------------------------------------
void OneStep::obtainSparseFunctionDerivativesForTranAdjoint ()
{
  if (sec.currentOrder_ == 2)
  {
    ds.nextDqdpDerivPtrVector->linearCombo(1.0,*ds.nextDqdpPtrVector,-1.0, *(ds.dqdpHistory[0]));
    ds.nextDqdpDerivPtrVector->scale(1.0/sec.currentTimeStep);

    ds.sensRHSPtrVector->linearCombo(1.0,*ds.nextDqdpDerivPtrVector,+0.5,*ds.nextDfdpPtrVector);
    ds.sensRHSPtrVector->update(+0.5,*ds.dfdpHistory[0]);

    ds.sensRHSPtrVector->update(-0.5,*ds.nextDbdpPtrVector);
    ds.sensRHSPtrVector->update(-0.5,*ds.dbdpHistory[0]);
  }
  else
  {
    ds.nextDqdpDerivPtrVector->linearCombo(1.0,*ds.nextDqdpPtrVector,-1.0, *(ds.dqdpHistory[0]));
    ds.nextDqdpDerivPtrVector->scale(1.0/sec.currentTimeStep);

    ds.sensRHSPtrVector->linearCombo(1.0,*ds.nextDqdpDerivPtrVector,+1.0,*ds.nextDfdpPtrVector);
    ds.sensRHSPtrVector->update(-1.0,*ds.nextDbdpPtrVector);
  }

  // since the nonlinear solver is expecting a -dFdp, scale by -1.0:
  ds.sensRHSPtrVector->scale(-1.0);

  // Collect sparse information from sensRHSPtrVector.
  ds.sparseSensRHSMV->replaceValues( ds.masterIndexVector, *ds.sensRHSPtrVector );
}

//-----------------------------------------------------------------------------
// Function      : OneStep::obtainAdjointSensitivityResidual
// Purpose       : Calculate Residual for adjoint sensitivities
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/23/2015
//-----------------------------------------------------------------------------
void OneStep::obtainAdjointSensitivityResidual()
{
  Linear::Vector & RHSVec  = *(ds.RHSVectorPtr);

  Linear::Vector & currDQdxLambda = *(ds.currDQdxLambdaPtr);
  Linear::Vector & currDFdxLambda = *(ds.currDFdxLambdaPtr);

  Linear::Vector & currLambda = *(ds.currLambdaPtr);

  Linear::Matrix & dQdx = *(ds.dQdxMatrixPtr);
  Linear::Matrix & dFdx = *(ds.dFdxMatrixPtr);

  // This assumes that the matvecs have been properly stored.
  // The trap coefficients are constant.
  int it=ds.itAdjointIndex;
  int itmax=ds.timeHistory.size();

  if (it<itmax-1) 
  {
    bool Transpose = true;

    {
    double qscalar(1.0/sec.lastTimeStep);
    currDQdxLambda.putScalar(0.0);
    dQdx.matvec( Transpose , currLambda, currDQdxLambda);
    RHSVec.update(+qscalar, currDQdxLambda);
    }

    if (ds.orderHistory[it+1] != 1)
    {
      double fscalar(-0.5);
      currDFdxLambda.putScalar(0.0);
      dFdx.matvec( Transpose , currLambda, currDFdxLambda);
      RHSVec.update(+fscalar, currDFdxLambda);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OneStep::obtainJacobian
// Purpose       : Calculate Jacobian
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------
void OneStep::obtainJacobian()
{
  if (DEBUG_TIME && isActive(Diag::TIME_JACOBIAN))
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() <<
      "  OneStep::obtainJacobian" << std::endl;
  }

  // output: ds.JMatrixPtr

  // This function returns the following matrix:
  // $-(sec.alphas_/hn)dQdx(x)+dFdx$

  // Note:  ds.nextSolutionPtr is used to get dQdx,dFdx in Analysis::AnalysisManager::loadJacobian.

  Linear::Matrix & dQdx = *(ds.dQdxMatrixPtr);
  Linear::Matrix & dFdx = *(ds.dFdxMatrixPtr);
  Linear::Matrix & Jac = *(ds.JMatrixPtr);

  double qscalar(-sec.alphas_/sec.currentTimeStep);
  double fscalar(1.0);
  if (sec.currentOrder_  == 2)
    fscalar =1.0/2.0;

  Jac.linearCombo( qscalar, dQdx, fscalar, dFdx );

  if (DEBUG_TIME && isActive(Diag::TIME_JACOBIAN))
  {
    Xyce::dout() << "fscalar = " << fscalar << "  qscalar = " << qscalar << std::endl;

    Xyce::dout() << "\n dFdx:" <<std::endl;
    dFdx.print(Xyce::dout());
    Xyce::dout() << "\n dQdx:" <<std::endl;
    dQdx.print(Xyce::dout());
    Xyce::dout() << "\n Total Jacobian:" <<std::endl;
    Jac.print(Xyce::dout());
    //    for (int i=0;i<3;++i)
    //    {
    //      printf("[ %25.20g\t%25.20g\t%25.20g ]\n",Jac[i][0],Jac[i][1],Jac[i][2]);
    //    }

    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : OneStep::interpolateSolution
// Purpose       : Interpolate solution approximation at prescribed time point.
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------
bool OneStep::interpolateSolution(double timepoint,
    Linear::Vector * tmpSolVectorPtr, std::vector<Linear::Vector*> & historyVec)
{
  // this is a very course approximation to determine if we are too
  // close to the actual time step to do an interpolation.
  // it could use more work.
  double dtr = timepoint - sec.currentTime;  // the delta time requested.

  if( -dtr < 2 * Util::MachineDependentParams::MachinePrecision() )
  {
    *tmpSolVectorPtr = *(historyVec[0]);
    return false;
  }

  if( sec.usedOrder_ <= 2)
  {
    // do first order interpolation
    // X_interp = X + delta_t_requested * delta_X/delta_t[last step]
    dtr = dtr / sec.lastTimeStep;
    tmpSolVectorPtr->linearCombo(1.0,*(historyVec[0]),dtr,*(historyVec[1]));
  }
  else
  {
    *tmpSolVectorPtr = *(historyVec[0]);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OneStep::interpolateMPDESolution
// Purpose       : Interpolate solution approximation at prescribed time points.
// Special Notes : This routine computes the solution at the output
//               : timepoints by intepolation of the history using the order
//               : used for the most recent completed step, orderUsed.
//               : The output is put into provided Linear::Vector pointer.
//               : The interpolation is as follows:
//               : tmpSolVectorPtr->block(i) is interpolated at timepoint(i)
//               : Therefore, if you want them all interpolated at the same time,
//               : then use timepoint(i) = timepoint(0) forall i
//               : or use interpolateSolution.
// Scope         : public
// Creator       : Ting Mei, Eric Keiter, SNL
// Creation Date : 11/28/06
//-----------------------------------------------------------------------------
bool OneStep::interpolateMPDESolution(std::vector<double>& timepoint,
    Linear::Vector * tmpSolVectorPtr)
{
  Linear::BlockVector * blockTempSolVectorPtr =
    dynamic_cast<Linear::BlockVector*>(tmpSolVectorPtr);
  if (blockTempSolVectorPtr == NULL)
  {
    Xyce::Report::DevelFatal0().in("OneStep::interpolateMPDESolution")
      << "Linear::Vector tmpSolVectorPtr is not of type Linear::BlockVector";
    return(false);
  }

  double tfuzz;   // fuzz factor to check for valid output time
  double tp;      // approximately t_{n-1}
  int numblocks = timepoint.size();
  int blockCount = blockTempSolVectorPtr->blockCount();
  if (numblocks > blockCount)
  {
    Xyce::Report::DevelFatal0().in("OneStep::interpolateMPDESolution")
      << "Number of time points requested is greater than number of fast time points in MPDE block vector";
    return(false);
  }
  double delt;
  double c = 1.0;
  double gam;
  int kord;       // order of interpolation
  double tn = sec.currentTime;
  double hh = sec.currentTimeStep;
  double hused = sec.usedStep_;
  int kused = sec.usedOrder_;
  double uround = 0.0;  // unit round-off (set to zero for now)

  tfuzz = 100 * uround * (tn + hh);
  tp = tn - hused - tfuzz;
  for (int i=0; i<numblocks ; ++i)
  {
    if ( (timepoint[i] - tp)*hh < 0.0 )
      return false;
  }

  *tmpSolVectorPtr = *(ds.xHistory[0]);

  Linear::Vector * solVectorPtr;
  Linear::Vector * xHistoryVectorPtr;
  // Loop over blocks first so that maximal order can be maintained
  for (int i=0; i < numblocks ; ++i)
  {
    if ((kused == 0) || (timepoint[i] == tn)) { kord = 1; }
    else { kord = kused; }
    solVectorPtr = &(blockTempSolVectorPtr->block(i));
    c = 1.0;
    delt = timepoint[i] - tn;
    gam = delt/sec.psi_[0];
    for (int j=1 ; j <= kord ; ++j)
    {
      c = c*gam;
      gam = (delt + sec.psi_[j-1])/sec.psi_[j];
      Linear::BlockVector * blockXHistoryVectorPtr =
        dynamic_cast<Linear::BlockVector*>(ds.xHistory[j]);
      if (blockXHistoryVectorPtr == NULL)
      {
        Xyce::Report::DevelFatal0().in("OneStep::interpolateMPDESolution") << "Linear::Vector ds.xHistory[j] is not of type Linear::BlockVector\n j = " << j;
        return(false);
      }
      xHistoryVectorPtr = &(blockXHistoryVectorPtr->block(i));
      solVectorPtr->update(c,*xHistoryVectorPtr);
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : OneStep::printMPDEOutputSolution()
// Purpose       : Print transient output from MPDE simulation
// Special Notes : This routine uses interpolateMPDESolution.
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/28/06
//-----------------------------------------------------------------------------
bool OneStep::printMPDEOutputSolution(
  Analysis::OutputMgrAdapter & outputManagerAdapter,
  const double time,
  Linear::Vector * solnVecPtr,
  const std::vector<double> & fastTimes )
{
  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() <<
      "  OneStep::printMPDEOutputSolution" << std::endl;
  }

  double timestep = sec.lastAttemptedTimeStep;
  double lasttime = sec.currentTime - timestep;
  double tn = sec.currentTime;
  // Set these values up to read output time intervals.  FIXME
  double beg_of_output_time_interval = lasttime;
  double end_of_output_time_interval = tn;
  double start_time = max(lasttime,beg_of_output_time_interval);
  double stop_time = min(tn,end_of_output_time_interval);
  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << "timestep = " << timestep << std::endl;
    Xyce::dout() << "lasttime = " << lasttime << std::endl;
    Xyce::dout() << "tn = " << tn << std::endl;
    Xyce::dout() << "beg_of_output_time_interval = " << beg_of_output_time_interval << std::endl;
    Xyce::dout() << "end_of_output_time_interval = " << end_of_output_time_interval << std::endl;
    Xyce::dout() << "start_time = " << start_time << std::endl;
    Xyce::dout() << "stop_time = " << stop_time << std::endl;
  }

  // Allocate tmp vectors in data store an initialize them.
  ds.allocateWaMPDEVectors();

  Linear::BlockVector * blockTmpSolVectorPtr =
    dynamic_cast<Linear::BlockVector*>(ds.tmpSolVectorPtr);
  if (blockTmpSolVectorPtr == NULL)
  {
    Xyce::Report::DevelFatal0().in("OneStep::printMPDEOutputSolution")
      << "Linear::Vector ds.tmpSolVectorPtr is not of type Linear::BlockVector";
    return(false);
  }
  int blockCount = blockTmpSolVectorPtr->blockCount();

  // Create list of timepoints to interpolate (along characteristic curve)
  double T2 = fastTimes.back();
  //double charcross = start_time - floor(start_time/T2)*T2; // (start_time mod T2)
  double charcross = fmod(start_time,T2); // (start_time mod T2)
  int s_ind_0 = -1;
  // find s_ind_0 = first fast time point >= charcross.
  // This could probably be made faster: FIXME
  for (int i=0 ; i<=blockCount ; ++i)
  {
    if (fastTimes[i] >= charcross)
    {
      s_ind_0 = i;
      break;
    }
  }
  if (s_ind_0 == -1)
  {
    Xyce::Report::DevelFatal0().in("OneStep::printMPDEOutputSolution")
      << "Cannot find where characteristic curve crosses fast time slice at start_time";
    return(false);
  }
  std::vector<double> h2(blockCount,0);
  for (int j=0 ; j < blockCount ; ++j)
  {
    h2[j] = fastTimes[j+1] - fastTimes[j];
  }
  std::vector<double> ti;
  //double first_interp = floor(start_time/T2)*T2 + fastTimes[s_ind_0];
  double first_interp = start_time - charcross + fastTimes[s_ind_0];
  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << "first_interp = " << first_interp << std::endl;
  }

  if (s_ind_0 == blockCount) { s_ind_0 = 0; };
  // Don't interpolate the first point
  double eps = fabs(start_time)*1.0e-6;
  if ( fabs(first_interp-timept_) <= eps )
  {
    first_interp += h2[s_ind_0];
    s_ind_0++;
    if (s_ind_0 == blockCount) { s_ind_0 = 0; };
    if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
    {
      Xyce::dout() << "Moving first_interp forward to avoid duplicate outputs:  " << first_interp << std::endl;
    }
  }
  int sn = s_ind_0;
  double t = first_interp;
  while (t <= stop_time)
  {
    ti.push_back(t);
    t += h2[sn];
    sn++;
    if (sn >= blockCount) { sn = 0; }
  }
  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << "T2 = " << T2 << std::endl;
    Xyce::dout() << "charcross = " << charcross << std::endl;
    Xyce::dout() << "s_ind_0 = " << s_ind_0 << std::endl;
    Xyce::dout() << "Expecting to interpolate the following points:" << std::endl;
    unsigned int numinterp = ti.size();
    for (unsigned int i=0 ; i < numinterp ; ++i)
    {
      Xyce::dout() << ti[i] << std::endl;
    }
    Xyce::dout() << "Total of " << numinterp << " points" << std::endl;
  }

  timept_ = start_time;  // used later for interpolating stop_time
  unsigned int tinum = ti.size();
  int total_interp = 0;
  std::vector<double> timepoint_vec(blockCount,stop_time);
  int num_interp_this_cycle = 0;
  int s_ind = s_ind_0;
  for (unsigned int i=0; i < tinum ; ++i)
  {
    timepoint_vec[s_ind] = ti[i];
    num_interp_this_cycle++;
    s_ind++;
    if (s_ind >= blockCount) { s_ind = 0; };
    // If we're back to s_ind_0 or out of ti points, then interpolate:
    if ((s_ind == s_ind_0) || (i == tinum-1))
    {
      interpolateMPDESolution(timepoint_vec, ds.tmpSolVectorPtr);
      // Now we print these points out
      int s = s_ind_0;
      for (int j=0 ; j < num_interp_this_cycle ; ++j)
      {
        double dt = 0.0; 
        // ERK. dt=0 is a kludge.  I don't have time to figure this function 
        // out.  This dt is wrong, but it only matters for expression output 
        // that depends on SDT or DDT, which probably don't work anyway for 
        // this case of interpolated MPDE
        timept_ = timepoint_vec[s];
        outputManagerAdapter.tranOutput(
	  timept_, dt, sec.finalTime,
          blockTmpSolVectorPtr->block(s),
          *ds.tmpStaVectorPtr, *ds.tmpStoVectorPtr, *ds.tmpLeadCurrentVectorPtr, *ds.tmpLeadDeltaVPtr, *ds.tmpLeadCurrentQDerivVectorPtr,
          ds.objectiveVec_, ds.dOdpVec_, ds.dOdpAdjVec_,
          ds.scaled_dOdpVec_, ds.scaled_dOdpAdjVec_);
        total_interp++;
        s++;
        if (s >= blockCount) { s = 0; }
        if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
        {
          Xyce::dout() << "Interpolated to t = " << timept_ << std::endl;
        }
      }
      num_interp_this_cycle = 0;
    }
  }
  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << "Total of " << total_interp << " points" << std::endl;
  }

  // Now we interpolate stop_time unless its too close to the last timept interpolated.
  eps = fabs(stop_time)*1.0e-8;
  // fudge factor for printing, this should come from elsewhere FIXME
  if (fabs(timept_ - stop_time) >= eps)
  {
    if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
    {
      Xyce::dout() << "Previous timept = " << timept_ << std::endl;
      Xyce::dout() << "Expecting to interpolate the following point: " << stop_time << std::endl;
    }
    Linear::Vector* tmpSolnVecPtr = solnVecPtr;
    Linear::Vector* tmpVecPtr = ds.tmpXn0APtr;
    if (stop_time < tn)
    {
      interpolateSolution(stop_time,ds.tmpXn0BPtr, ds.xHistory);
      tmpSolnVecPtr = ds.tmpXn0BPtr;
    }
    Linear::BlockVector * blockTmpSolnVecPtr =
      dynamic_cast<Linear::BlockVector*>(tmpSolnVecPtr);
    Linear::BlockVector * blockTmpVecPtr =
      dynamic_cast<Linear::BlockVector*>(tmpVecPtr);
    if (blockTmpSolnVecPtr == NULL)
    {
      Xyce::Report::DevelFatal0().in("OneStep::printMPDEOutputSolution")
        << "Linear::Vector tmpSolnVecPtr is not of type Linear::BlockVector";
      return(false);
    }
    if (blockTmpVecPtr == NULL)
    {
      Xyce::Report::DevelFatal0().in("OneStep::printMPDEOutputSolution")
        << "Linear::Vector tmpVecPtr is not of type Linear::BlockVector";
      return(false);
    }
    // Interpolate where the caracteristic crosses stop_time (instead of start_time).
    //charcross = stop_time - floor(stop_time/T2)*T2;
    charcross = fmod(stop_time,T2);
    int s_ind_1 = -1;
    // Find index just before charcross:
    if( charcross < fastTimes[0] )
    {
      // by periodicity, the one before in this case is the one at the end
      s_ind_1 = blockCount-1;
    }
    else
    {
      for (int i=blockCount-1 ; i>=0 ; --i)
      {
        if (fastTimes[i] <= charcross)
        {
          s_ind_1 = i;
          break;
        }
      }
    }
    if (s_ind_1 == -1)
    {
      Xyce::Report::DevelFatal0().in("OneStep::printMPDEOutputSolution")
        << "Cannot find where characteristic curve crosses fast time slice at stop_time";
      return(false);
    }
    int sm = s_ind_1;
    int sp = s_ind_1+1;
    double coeff_sm = fastTimes[sp]-charcross;
    double coeff_sp = charcross-fastTimes[sm];
    if (sp == blockCount) { sp = 0; }
    double dt = h2[s_ind_1];
    timept_ = stop_time;
    if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
    {
      Xyce::dout() << "charcross = " << charcross << std::endl;
      Xyce::dout() << "s_ind_1 = " << s_ind_1 << std::endl;
      Xyce::dout() << "sp = " << sp << std::endl;
      Xyce::dout() << "sm = " << sm << std::endl;
      Xyce::dout() << "dt = " << dt << std::endl;
      Xyce::dout() << "timept = " << timept_ << std::endl;
      Xyce::dout() << "coeff_sm = " << coeff_sm << std::endl;
      Xyce::dout() << "coeff_sp = " << coeff_sp << std::endl;
    }
    blockTmpVecPtr->block(0).linearCombo(
      coeff_sm/dt, blockTmpSolnVecPtr->block(sm),
      coeff_sp/dt, blockTmpSolnVecPtr->block(sp)
                                         );
    outputManagerAdapter.tranOutput(
        timept_, dt, sec.finalTime,
        blockTmpVecPtr->block(0), *ds.tmpStaVectorPtr, *ds.tmpStoVectorPtr, 
        *ds.tmpLeadCurrentVectorPtr, *ds.tmpLeadDeltaVPtr, *ds.tmpLeadCurrentQDerivVectorPtr,
                                    ds.objectiveVec_, ds.dOdpVec_, ds.dOdpAdjVec_,
                                    ds.scaled_dOdpVec_, ds.scaled_dOdpAdjVec_);
    if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
    {
      Xyce::dout() << "Interpolated to t = " << timept_ << std::endl;
    }
  }
  else if(DEBUG_TIME && isActive(Diag::TIME_OUTPUT))  // no extra interpolation
  {
    Xyce::dout() << "No further interpolation required." << std::endl;
  }

  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OneStep::printWaMPDEOutputSolution()
// Purpose       : Print transient output from WaMPDE simulation
// Special Notes : This routine uses interpolateSolution.
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 12/15/06
//-----------------------------------------------------------------------------
bool OneStep::printWaMPDEOutputSolution(
    Analysis::OutputMgrAdapter & outputManagerAdapter,
    const double time,
    Linear::Vector * solnVecPtr,
    const std::vector<double> & fastTimes,
    const int phiGID )
{
  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() <<
      "  OneStep::printWaMPDEOutputSolution" << std::endl;
  }

  double timestep = sec.lastAttemptedTimeStep;
  double lasttime = sec.currentTime - timestep;
  double tn = sec.currentTime;
  // Set these values up to read output time intervals.  FIXME
  double beg_of_output_time_interval = lasttime;
  double end_of_output_time_interval = tn;
  double start_time = max(lasttime,beg_of_output_time_interval);
  double stop_time = min(tn,end_of_output_time_interval);
  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << "start_time = " << start_time << std::endl;
    Xyce::dout() << "stop_time = " << stop_time << std::endl;
  }

  // 12/15/06 tscoffe:
  // First break the interval up as printOutputSolution does:
  // hh = timestep/(sec.usedOrder_) and interpolate in the intervals:
  // [tn+(i-1)*hh,tn+i*hh], i=1..usedOrder
  // Assume phi(t_1) is linear in these intervals and approximate with:
  // phi(t) = (1/hh)*(phi(tn)(tn+hh-t)+phi(tn+hh)(t-tn))
  // Then the t_1 values we want to interpolate are:
  // n2 = number of fast time points.
  // T2 = fast time period
  // h2 = T2/n2 = average spacing on fast time scale
  // t1_vals = [tn:h2:tn+hh]
  // And the t_2 values we want to interpolate are:
  // t2_vals = phi(t1_vals) mod T2
  // Then take the N_LAS blocks and do 2D linear interpolation on the intervals:
  // (t1,s1), (t1,s2), (t2,s1), (t2,s2)
  // x(t) = x(t,s) approx =
  // (1/(t2-t1))(1/(s2-s1))[  x(t1,s1)(t2-t)(s2-s)
  //                         +x(t1,s2)(t2-t)(s-s1)
  //                         +x(t2,s1)(t-t1)(s2-s)
  //                         +x(t2,s2)(t-t1)(s-s1) ]
  // where t = t1_vals and s = t2_vals

  // Allocate tmp vectors in data store an initialize them.
  ds.allocateWaMPDEVectors();

  Linear::BlockVector * blockTmpSolVectorPtr =
    dynamic_cast<Linear::BlockVector*>(ds.tmpSolVectorPtr);
  Linear::BlockVector * blockTmpXn0APtr =
    dynamic_cast<Linear::BlockVector*>(ds.tmpXn0APtr);
  Linear::BlockVector * blockTmpXn0BPtr =
    dynamic_cast<Linear::BlockVector*>(ds.tmpXn0BPtr);
  if (blockTmpSolVectorPtr == NULL)
  {
    Xyce::Report::DevelFatal0().in("OneStep::printWaMPDEOutputSolution")
      << "Linear::Vector ds.tmpSolVectorPtr is not of type Linear::BlockVector";
    return(false);
  }
  if (blockTmpXn0APtr == NULL)
  {
    Xyce::Report::DevelFatal0().in("OneStep::printWaMPDEOutputSolution")
      << "Linear::Vector ds.tmpXn0APtr is not of type Linear::BlockVector";
    return(false);
  }
  if (blockTmpXn0BPtr == NULL)
  {
    Xyce::Report::DevelFatal0().in("OneStep::printWaMPDEOutputSolution")
      << "Linear::Vector ds.tmpXn0BPtr is not of type Linear::BlockVector";
    return(false);
  }
  int phiLID = blockTmpSolVectorPtr->pmap()->globalToLocalIndex(phiGID);
  double hh = timestep/(sec.usedOrder_);
  double timeA = -1.0;
  double timeB = -1.0;
  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << " sec.usedOrder_ = " << sec.usedOrder_ << std::endl;
    Xyce::dout() << " sec.currentTime_ = " << sec.currentTime << std::endl;
    Xyce::dout() << " lasttime = " << lasttime << std::endl;
  }

  for (int i=0 ; i < sec.usedOrder_ ; ++i)
  {
    if (i == 0)
    {
      bool junk;
      timeA = lasttime + hh*i;
      junk = interpolateSolution(timeA,ds.tmpXn0APtr, ds.xHistory);
      if (!junk)
      {
        Xyce::Report::DevelFatal0().in("OneStep::printWaMPDEOutputSolution")
          <<  "interpolateSolution returned false!";
      }
    }
    else
    {
      // We don't need to interpolate this again.
      *ds.tmpXn0APtr = *ds.tmpXn0BPtr;
      timeA = timeB;
    }
    timeB = lasttime + hh*(i+1);
    interpolateSolution(timeB,ds.tmpXn0BPtr, ds.xHistory);
    if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
    {
      Xyce::dout() << "Interpolating in [ " << timeA << ", " << timeB << " ]" << std::endl;
      Xyce::dout() << "timeA = " << timeA << std::endl;
      Xyce::dout() << "timeB = " << timeB << std::endl;
    }

    // Now we can interpolate [tmpXn0APtr,tmpXn0BPtr] in [timeA,timeB].
    std::vector<double> t1vals;
    double T2 = fastTimes.back();
    int blockCount = blockTmpSolVectorPtr->blockCount();
    double h2 = T2/blockCount; // Average mesh width in fast time-scale
    double tval = timeA+h2;
    while (tval <= timeB)
    {
      t1vals.push_back(tval);
      tval += h2;
    }
    // fudge factor for printing, this should come from elsewhere FIXME
    double eps = fabs(timeB)*1.0e-8;
    if ( (t1vals.size() == 0) || (fabs(t1vals.back() - timeB) >= eps) )
    {
      t1vals.push_back(timeB);
    }
    std::vector<double> t2vals, phiAB(2);
    std::vector<double> tmpPhiAB(2, 0.0);
    if (phiLID >= 0)
    {
      tmpPhiAB[0] = (*ds.tmpXn0APtr)[phiLID]; // Get from MPDE Manager
      tmpPhiAB[1] = (*ds.tmpXn0BPtr)[phiLID]; // Get from MPDE Manager
    }
    blockTmpSolVectorPtr->pmap()->pdsComm().sumAll( &tmpPhiAB[0], &phiAB[0], 2 );

    double phiA = phiAB[0], phiB = phiAB[1];
    for (unsigned int j=0 ; j<t1vals.size() ; ++j)
    {
      double phi = (1/(timeB-timeA))*(phiA*(timeB-t1vals[j])+phiB*(t1vals[j]-timeA));
      t2vals.push_back(fmod(phi,T2)); // phi(t1vals[j]) mod T2
    }
    if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
    {
      Xyce::dout() << "t1vals = " << std::endl;
      for (unsigned int j=0 ; j < t1vals.size() ; ++j)
      {
        Xyce::dout() << t1vals[j] << std::endl;
      }
      Xyce::dout() << "phi(" << timeA << ") = " << phiA << std::endl;
      Xyce::dout() << "phi(" << timeB << ") = " << phiB << std::endl;
      Xyce::dout() << "t2vals = " << std::endl;
      for (unsigned int j=0 ; j< t2vals.size() ; ++j)
      {
        Xyce::dout() << t2vals[j] << std::endl;
      }
    }

    // Now we can do our block 2D interpolations
    // loop through t1vals and move the fast time blocks as we go
    double t1 = timeA; // slow time indices
    double t2 = timeB;
    double t = t1vals[0]; // current time indices
    double s = t2vals[0];
    int b1,b2; // fast time block indices corresponding to s1,s2
    // Find the block surrounding s:
    b1 = -2;
    for (int j=0 ; j < blockCount ; ++j)
    {
      if ((fastTimes[j] <= s) && (s < fastTimes[j+1]))
      {
        b1 = j;
      }
    }
    b2 = b1+1;
    if (b2 == blockCount)
    {
      b2 = 0;
    }
    double s1 = fastTimes[b1];
    double s2 = fastTimes[b1+1]; // Note:  fastTimes[blockCount] = T2
    if ((s < s1) || (s > s2))
    {
      Xyce::Report::DevelFatal0().in("OneStep::printWaMPDEOutputSolution")
        << "  Interpolator cannot find a fast time block containing the first point  ";
    }
    if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
    {
      Xyce::dout() << "Found s = " << s << " in block " << b1 << " with boundary = [" << s1 << "," << s2 << "]" << std::endl;
    }

    for (unsigned int j=0 ; j < t1vals.size() ; ++j)
    {
      t = t1vals[j];
      s = t2vals[j];
      if (t > t2) break; // This should never occur
      // If s has moved outside our block, then increment block.
      if ( (s < s1) || (s > s2) )
      {
        if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
        {
          Xyce::dout() << "Incrementing fast time block for next interpolation." << std::endl;
        }

        b1++;
        if (b1 == blockCount)
        {
          b1 = 0;
        }
        b2 = b1+1;
        if (b2 == blockCount)
        {
          b2 = 0;
        }
        s1 = fastTimes[b1];
        s2 = fastTimes[b1+1];
      }
      // If s isn't in the next block, then search for it.
      if ((s < s1) || (s > s2))
      {
        if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
        {
          Xyce::dout() << "Searching for fast time block for next interpolation." << std::endl;
        }

        b1 = -2;
        for (int j2=0 ; j2 < blockCount ; ++j2)
        {
          if ((fastTimes[j2] <= s) && (s < fastTimes[j2+1]))
          {
            b1 = j2;
          }
        }
        b2 = b1+1;
        if (b2 == blockCount)
        {
          b2 = 0;
        }
        s1 = fastTimes[b1];
        s2 = fastTimes[b1+1];
      }
      // If a block surrounding s can't be found, then quit.
      if ((s < s1) || (s > s2))
      {
        Xyce::Report::DevelFatal0().in("OneStep::printWaMPDEOutputSolution")
          << "  Interpolator moved fast time block but new point is not in this block  ";
      }
      if (s > T2) break; // Just double checking...
      //blockTmpXn0APtr->block(b1) // x(t1,s1)
      //blockTmpXn0APtr->block(b2) // x(t1,s2)
      //blockTmpXn0BPtr->block(b1) // x(t2,s1)
      //blockTmpXn0BPtr->block(b2) // x(t2,s2)
      // (1/(t2-t1))(1/(s2-s1))[  x(t1,s1)(t2-t)(s2-s)
      //                         +x(t1,s2)(t2-t)(s-s1)
      //                         +x(t2,s1)(t-t1)(s2-s)
      //                         +x(t2,s2)(t-t1)(s-s1) ]
      if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
      {
        Xyce::dout() << "Interpolating in the block:" << std::endl;
        Xyce::dout() << "(t1,t2) = (" << t1 << "," << t2 << ")" << std::endl;
        Xyce::dout() << "(s1,s2) = (" << s1 << "," << s2 << ")" << std::endl;
      }

      double denom = (t2-t1)*(s2-s1);
      double coeff0 = (t2-t)*(s2-s)/denom;
      double coeff1 = (t2-t)*(s-s1)/denom;
      double coeff2 = (t-t1)*(s2-s)/denom;
      double coeff3 = (t-t1)*(s-s1)/denom;
      (blockTmpSolVectorPtr->block(b1)).linearCombo(
          coeff0, blockTmpXn0APtr->block(b1),
          coeff1, blockTmpXn0APtr->block(b2));
      (blockTmpSolVectorPtr->block(b1)).update(
          coeff2, blockTmpXn0BPtr->block(b1),
          coeff3, blockTmpXn0BPtr->block(b2), 1.0 );

      // erkeite 2/24/07. This is needed because currently the interpolation goes back to t=-1.0.
      if (t >= 0.0)
      {
        double dt=0.0;
        // ERK. dt=0 is a kludge.  I don't have time to figure this function 
        // out.  This dt is wrong, but it only matters for expression output 
        // that depends on SDT or DDT, which probably don't work anyway for 
        // this case of interpolated WaMPDE
        outputManagerAdapter.tranOutput(t, dt, sec.finalTime, blockTmpSolVectorPtr->block(b1),
            *ds.tmpStaVectorPtr, *ds.tmpStoVectorPtr, 
            *ds.tmpLeadCurrentVectorPtr, *ds.tmpLeadDeltaVPtr, *ds.tmpLeadCurrentQDerivVectorPtr,
            ds.objectiveVec_, ds.dOdpVec_, ds.dOdpAdjVec_,
            ds.scaled_dOdpVec_, ds.scaled_dOdpAdjVec_);
        if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
        {
          Xyce::dout() << "Interpolated to (t,phi(t)) = (" << t << "," << s << ")" << std::endl;
        }
      }
    }
  }

  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OneStep::printOutputSolution()
// Purpose       : Print output that is dumbed down in order.
// Special Notes : This routine picks smaller time steps to approximate first
//               : order integration from the perspective of the output.
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------
bool OneStep::printOutputSolution(
  Analysis::OutputMgrAdapter &   outputManagerAdapter,
  const TIAParams &              tia_params,
  const double                   time,
  Linear::Vector *                 solnVecPtr,
  const bool                     doNotInterpolate,
  const std::vector<double>     &outputInterpolationTimes,
  bool                           skipPrintLineOutput)
{
  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() <<
      "  OneStep::printOutputSolution" << std::endl;
    Xyce::dout() << "usedOrder_ = " << sec.usedOrder_ << std::endl;
  }

  double timestep = sec.lastAttemptedTimeStep;
  bool dointerp = true;
  double hh = timestep/(sec.usedOrder_);

  if (hh <= 10*sec.minTimeStep)
  {
    dointerp = false;
  }

  if (!tia_params.interpOutputFlag)
  {
    dointerp = false;
  }

  if (doNotInterpolate)
  {
    dointerp = false;
  }

  if (dointerp && !outputInterpolationTimes.empty())
  {
    double dt=0.0;
    for (unsigned int i=0;i<outputInterpolationTimes.size();++i)
    {
      interpolateSolution(outputInterpolationTimes[i], ds.tmpSolVectorPtr, ds.xHistory);  // interpolate solution
      interpolateSolution(outputInterpolationTimes[i], ds.tmpStaVectorPtr, ds.sHistory);  // interpolate state
      interpolateSolution(outputInterpolationTimes[i], ds.tmpStoVectorPtr, ds.stoHistory);// interpolate store

      if (ds.leadCurrentSize)
      {
        interpolateSolution(outputInterpolationTimes[i], ds.tmpLeadCurrentVectorPtr, ds.leadCurrentHistory);// interpolate lead current
        interpolateSolution(outputInterpolationTimes[i], ds.tmpLeadCurrentQDerivVectorPtr, ds.leadCurrentQDerivHistory);// interpolate lead current q comp.
        interpolateSolution(outputInterpolationTimes[i], ds.tmpLeadDeltaVPtr, ds.leadDeltaVHistory);// interpolate junction voltage
      }

      if (i>0)
      {
        dt = outputInterpolationTimes[i]-outputInterpolationTimes[i-1]; 
      }
      else
      {
        dt = 0.0;
      }
      outputManagerAdapter.tranOutput(
         outputInterpolationTimes[i], dt, sec.finalTime,
         *ds.tmpSolVectorPtr,
         *ds.tmpStaVectorPtr, *ds.tmpStoVectorPtr, 
         *ds.tmpLeadCurrentVectorPtr, *ds.tmpLeadDeltaVPtr, *ds.tmpLeadCurrentQDerivVectorPtr,
         ds.objectiveVec_, ds.dOdpVec_, ds.dOdpAdjVec_,
         ds.scaled_dOdpVec_, ds.scaled_dOdpAdjVec_,
         skipPrintLineOutput);
    }
  }

  // Either way, do an output on the actual computed time step, but only
  // if we weren't given a list of specific times *or* we were told not to
  // interpoloate.
  if (outputInterpolationTimes.empty() || doNotInterpolate)
  {
    outputManagerAdapter.tranOutput(time, timestep, sec.finalTime,
        *ds.currSolutionPtr,
        *ds.currStatePtr, *ds.currStorePtr, *ds.currLeadCurrentPtr, *ds.currLeadDeltaVPtr, *ds.currLeadCurrentQPtr,
        ds.objectiveVec_, ds.dOdpVec_, ds.dOdpAdjVec_,
        ds.scaled_dOdpVec_, ds.scaled_dOdpAdjVec_,
        skipPrintLineOutput);
  }

  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OneStep::saveOutputSolution
// Purpose       : This is similar to printOutputSolution, but is in support of
//                 the .SAVE capability, rather than .PRINT.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
bool OneStep::saveOutputSolution(
  Parallel::Machine                     comm,
  IO::InitialConditionsManager &        initial_conditions_manager,
  const NodeNameMap &                   node_name_map,
  const TIAParams &                     tia_params,
  Linear::Vector *                      solnVecPtr,
  const double                          saveTime,
  const bool                            doNotInterpolate)
{
  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() << "  OneStep::saveOutputSolution" << std::endl;
  }

  if (!doNotInterpolate)
  {
    interpolateSolution(saveTime, ds.tmpSolVectorPtr, ds.xHistory);  // interpolate solution
    initial_conditions_manager.outputDCOP(comm, node_name_map, *ds.tmpSolVectorPtr);
  }
  else 
  {
    initial_conditions_manager.outputDCOP(comm, node_name_map, *ds.currSolutionPtr);
  }

  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OneStep::updateHistory
// Purpose       : Update history array after a successful step
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------
void OneStep::updateHistory()
{
  if (DEBUG_TIME && isActive(Diag::TIME_HISTORY))
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() <<
      "  OneStep::updateHistory" << std::endl;
    Xyce::dout() << "\n Before updates \n" << std::endl;
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      Xyce::dout() << "\n xHistory["<< i << "]: \n" << std::endl;
      (ds.xHistory[i])->print(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      Xyce::dout() << "\n qHistory["<< i << "]: \n" << std::endl;
      (ds.qHistory[i])->print(Xyce::dout());
      Xyce::dout() << std::endl;
    }

    Xyce::dout() << "\n sHistory["<< 0 << "]: \n" << std::endl;
    (ds.sHistory[0])->print(Xyce::dout());
    Xyce::dout() << std::endl;
    
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  if (sec.currentOrder_ == 2)
  {
    *(ds.xHistory[2]) = *(ds.xHistory[1]);
    (ds.qHistory[2])->linearCombo(1.0, *(ds.daeFVectorPtr), -1.0, *(ds.daeBVectorPtr));
  }

  (ds.xHistory[1])->linearCombo(1.0, *ds.nextSolutionPtr, -1.0,*(ds.xHistory[0]));
  (ds.qHistory[1])->linearCombo(1.0, *ds.daeQVectorPtr, -1.0,*(ds.qHistory[0]));

  (ds.stoHistory[1])->linearCombo(1.0, *ds.nextStorePtr, -1.0,*(ds.stoHistory[0]));

  *(ds.xHistory[0]) = *ds.nextSolutionPtr;
  *(ds.qHistory[0]) =  *ds.daeQVectorPtr;
  *(ds.sHistory[0]) =  *ds.nextStatePtr;
  *(ds.stoHistory[0]) =  *ds.nextStorePtr;
  
  if (ds.leadCurrentSize)
  {
    (ds.leadCurrentHistory[1])->linearCombo(1.0, *ds.nextLeadCurrentPtr, -1.0,*(ds.leadCurrentHistory[0]));
    (ds.leadCurrentQHistory[1])->linearCombo(1.0, *ds.nextLeadCurrentQPtr, -1.0,*(ds.leadCurrentQHistory[0]));
    (ds.leadDeltaVHistory[1])->linearCombo(1.0, *ds.nextLeadDeltaVPtr, -1.0,*(ds.leadDeltaVHistory[0]));

    *(ds.leadCurrentHistory[0]) = *ds.nextLeadCurrentPtr;
    *(ds.leadCurrentQHistory[0]) = *ds.nextLeadCurrentQPtr;
    *(ds.leadDeltaVHistory[0]) = *ds.nextLeadDeltaVPtr;
  }

  if (DEBUG_TIME && isActive(Diag::TIME_HISTORY))
  {
    Xyce::dout() << "\n After updates \n" << std::endl;
    Xyce::dout() << "\n newtonCorrectionPtr: " << std::endl;
    ds.newtonCorrectionPtr->print(Xyce::dout());
    Xyce::dout() << "\n qNewtonCorrectionPtr: " << std::endl;
    ds.qNewtonCorrectionPtr->print(Xyce::dout());
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      Xyce::dout() << "\n xHistory["<< i << "]: \n" << std::endl;
      (ds.xHistory[i])->print(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      Xyce::dout() << "\n qHistory["<< i << "]: \n" << std::endl;
      (ds.qHistory[i])->print(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    Xyce::dout() << std::endl;
    Xyce::dout() << "\n nextStatePtr: " << std::endl;
    ds.nextStatePtr->print(Xyce::dout());
    Xyce::dout() << std::endl;
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      Xyce::dout() << "\n sHistory["<< i << "]: \n" << std::endl;
      (ds.sHistory[i])->print(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  updateSensitivityHistory();
}

//-----------------------------------------------------------------------------
// Function      : OneStep::updateSensitivityHistory
// Purpose       : Update sensitivity history array after a successful step
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void OneStep::updateSensitivityHistory()
{
  if (ds.numParams)
  {
    if (sec.currentOrder_ == 2)
    {
      (ds.dqdpHistory[2])->linearCombo(1.0, *ds.nextDfdpPtrVector, -1.0, *ds.nextDbdpPtrVector);
    }

    (ds.dqdpHistory[1])->linearCombo(1.0, *ds.nextDqdpPtrVector, -1.0, *(ds.dqdpHistory[0]));
    (ds.dfdpHistory[1])->linearCombo(1.0, *ds.nextDfdpPtrVector, -1.0, *(ds.dfdpHistory[0]));
    (ds.dbdpHistory[1])->linearCombo(1.0, *ds.nextDbdpPtrVector, -1.0, *(ds.dbdpHistory[0]));

    *(ds.dqdpHistory[0]) = *(ds.nextDqdpPtrVector);
    *(ds.dfdpHistory[0]) = *(ds.nextDfdpPtrVector);
    *(ds.dbdpHistory[0]) = *(ds.nextDbdpPtrVector);
  }
}


//-----------------------------------------------------------------------------
// Function      : OneStep::updateAdjointSensitivityHistory
// Purpose       : Update sensitivity history array after a successful step
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/23/2015
//-----------------------------------------------------------------------------
void OneStep::updateAdjointSensitivityHistory()
{
  if (ds.numParams)
  {
    *(ds.lastLambdaPtr) = *(ds.currLambdaPtr);
    *(ds.currLambdaPtr) = *(ds.nextLambdaPtr);

    *(ds.lastDQdxLambdaPtr) = *(ds.currDQdxLambdaPtr);
    *(ds.currDQdxLambdaPtr) = *(ds.nextDQdxLambdaPtr);

    *(ds.lastDFdxLambdaPtr) = *(ds.currDFdxLambdaPtr);
    *(ds.currDFdxLambdaPtr) = *(ds.nextDFdxLambdaPtr);
  }
}

//-----------------------------------------------------------------------------
// Function      : OneStep::restoreHistory
// Purpose       : Restore history array after a failed step
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
void OneStep::restoreHistory()
{
  for (int i=1;i<=sec.currentOrder_;++i)
  {
    sec.psi_[i-1] = sec.psi_[i];
  }
  if (DEBUG_TIME && isActive(Diag::TIME_HISTORY))
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() <<
      "  OneStep::restoreHistory" << std::endl;
    for (int i=1;i<=sec.currentOrder_;++i)
      Xyce::dout() << "\n sec.psi_[i] = " << sec.psi_[i] << std::endl;
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      Xyce::dout() << "\n xHistory["<< i << "]: \n" << std::endl;
      (ds.xHistory[i])->print(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      Xyce::dout() << "\n qHistory["<< i << "]: \n" << std::endl;
      (ds.qHistory[i])->print(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    for (int i=0; i<=sec.maxOrder_ ; ++i)
    {
      Xyce::dout() << "\n sHistory["<< i << "]: \n" << std::endl;
      (ds.sHistory[i])->print(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    Xyce::dout() << Xyce::section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : OneStep::updateCoeffs
// Purpose       : Update method coefficients
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
void OneStep::updateCoeffs()
{
  // synchronize with Step Error Control
  //  sec.psi_[0] = sec.currentTimeStep;
  if (DEBUG_TIME && isActive(Diag::TIME_COEFFICIENTS))
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() << "  OneStep::updateCoeffs" << std::endl;
    Xyce::dout() << "  currentTimeStep = " <<  sec.currentTimeStep << std::endl;
    Xyce::dout() << "  numberOfSteps_ = " <<  sec.numberOfSteps_ << std::endl;
    Xyce::dout() << "  currentOrder_ = " <<  sec.currentOrder_ << std::endl;
    Xyce::dout() << "  nscsco_ = " <<  sec.nscsco_ << std::endl;
    Xyce::dout() << "  psi_[0] = " <<  sec.psi_[0] << std::endl;
  }

  double temp1 = sec.currentTimeStep;

  if (sec.currentOrder_ == 2)
  {
    sec.psi_[2] = sec.psi_[1];
  }
  sec.psi_[1] = sec.psi_[0];
  sec.psi_[0] = temp1;

  sec.beta_[0] = 1.0;
  sec.alpha_[0] = 1.0;
  //    sec.ck_ = 1.0;
  sec.alphas_ = -1.0;

  if (sec.currentOrder_ == 2)
  {
    double temp2 = sec.psi_[1];
    sec.beta_[1] = temp1/temp2 + (temp1/temp2)*(temp1/temp2)/2;
    sec.beta_[2] = -1.0*temp1*temp1/temp2/sec.psi_[2]/2;
    sec.ck_ = (sec.currentTimeStep/sec.currentTimeStepSum)/3.0;
  }
  else
  {
    sec.beta_[1] = temp1/sec.psi_[1];
    sec.ck_ = sec.currentTimeStep/sec.currentTimeStepSum;
  }

  if (DEBUG_TIME && isActive(Diag::TIME_COEFFICIENTS))
  {
    Xyce::dout() << "  nscsco_ = " <<  sec.nscsco_ << std::endl;
    Xyce::dout() << "  beta_[0] = " <<  sec.beta_[0] << std::endl;
    Xyce::dout() << "  beta_[1] = " <<  sec.beta_[1] << std::endl;
    Xyce::dout() << "  beta_[2] = " <<  sec.beta_[2] << std::endl;
    Xyce::dout() << "  beta_[3] = " <<  sec.beta_[3] << std::endl;
    Xyce::dout() << "  beta_[4] = " <<  sec.beta_[4] << std::endl;
    Xyce::dout() << "  alpha_[0] = " <<  sec.alpha_[0] << std::endl;
    Xyce::dout() << "  alpha_[1] = " <<  sec.alpha_[1] << std::endl;
    Xyce::dout() << "  alpha_[2] = " <<  sec.alpha_[2] << std::endl;
    Xyce::dout() << "  alpha_[3] = " <<  sec.alpha_[3] << std::endl;
    Xyce::dout() << "  alpha_[4] = " <<  sec.alpha_[4] << std::endl;
    Xyce::dout() << "  alphas_ = " <<  sec.alphas_ << std::endl;
    Xyce::dout() << "  alpha0_ = " <<  sec.alpha0_ << std::endl;
    Xyce::dout() << "  gamma_[0] = " <<  sec.gamma_[0] << std::endl;
    Xyce::dout() << "  gamma_[1] = " <<  sec.gamma_[1] << std::endl;
    Xyce::dout() << "  gamma_[2] = " <<  sec.gamma_[2] << std::endl;
    Xyce::dout() << "  gamma_[3] = " <<  sec.gamma_[3] << std::endl;
    Xyce::dout() << "  gamma_[4] = " <<  sec.gamma_[4] << std::endl;
    Xyce::dout() << "  psi_[0] = " <<  sec.psi_[0] << std::endl;
    Xyce::dout() << "  psi_[1] = " <<  sec.psi_[1] << std::endl;
    Xyce::dout() << "  psi_[2] = " <<  sec.psi_[2] << std::endl;
    Xyce::dout() << "  psi_[3] = " <<  sec.psi_[3] << std::endl;
    Xyce::dout() << "  psi_[4] = " <<  sec.psi_[4] << std::endl;
    Xyce::dout() << "  sigma_[0] = " <<  sec.sigma_[0] << std::endl;
    Xyce::dout() << "  sigma_[1] = " <<  sec.sigma_[1] << std::endl;
    Xyce::dout() << "  sigma_[2] = " <<  sec.sigma_[2] << std::endl;
    Xyce::dout() << "  sigma_[3] = " <<  sec.sigma_[3] << std::endl;
    Xyce::dout() << "  sigma_[4] = " <<  sec.sigma_[4] << std::endl;
    Xyce::dout() << "  ck_ = " <<  sec.ck_ << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : OneStep::updateAdjointCoeffs
// Purpose       : Update method coefficients for adjoint.    (just order)
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void OneStep::updateAdjointCoeffs()
{
  int it=ds.itAdjointIndex;
  sec.currentOrder_  = ds.orderHistory[it] ;
}

//-----------------------------------------------------------------------------
// Function      : OneStep::initialize
// Purpose       : Initialize method with initial solution & step-size
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
void OneStep::initialize(const TIAParams &tia_params)
{
  // we assume the solution vector is available here
  // Note that I'm using currSolutionPtr instead of
  // nextSolutionPtr because this is the first step.

  // Update next stop time from StepErrorControl:
  // ERK.  Commenting this out, as it is already called from Analysis::AnalysisManager,
  // right before this initialize call.  It should not be called 2x, as
  // it is history dependent (unfortunately), so calling it 2x in a row changes
  // the stop time to a different number.
  // sec.updateStopTime();

  // Choose initial step-size
  double time_to_stop = sec.stopTime - sec.currentTime;
  double currentTimeStep;

  sec.TimeStepLimitedbyBP =  false;

  if (tia_params.constantTimeStepFlag)
  {
    currentTimeStep = 0.1 * time_to_stop;
    currentTimeStep = std::min(sec.startingTimeStep, currentTimeStep);
    sec.currentTimeStep = currentTimeStep;
  }
  else
  {
    // compute an initial step-size based on rate of change in the
    // solution initially
    double dnorm_q = ds.delta_x_errorNorm_q1();
    if (dnorm_q > 0.0)  // time-dependent DAE
    {
      if (sec.currentTime == sec.initialTime)
        currentTimeStep = std::min(sec.h0_max_factor_*abs(time_to_stop),sqrt(2.0)/(sec.h0_safety_*dnorm_q));
      else
        currentTimeStep = 0.1* std::min(sec.savedTimeStep, abs(time_to_stop));
    }
    else  // non-time-dependent DAE
    {
      if (sec.currentTime == sec.initialTime)
        currentTimeStep = sec.h0_max_factor_*abs(time_to_stop);
      else
        currentTimeStep = 0.1* std::min(sec.savedTimeStep, abs(time_to_stop));
    }
    // choose min of user specified value and our value:
    if (sec.startingTimeStep > 0.0 && (sec.currentTime == sec.initialTime))
      currentTimeStep = std::min(sec.startingTimeStep, currentTimeStep);
    // check for maximum step-size:
    double rh = abs(currentTimeStep)*sec.h_max_inv_;
    if (rh>1.0) currentTimeStep = currentTimeStep/rh;


    // Apply this new stepsize only if it is smaller than the one preceding
    // the breakpoint, but only do this if this is a non-DCOP breakpoint.
    if (sec.currentTime != sec.initialTime) // if not DCOP:
    {
      sec.currentTimeStep = std::min(sec.currentTimeStep, currentTimeStep);
    }
    else // else if DCOP:
    {
      sec.currentTimeStep = currentTimeStep;
    }
  }

  sec.currentTimeStepRatio = 1.0;
  sec.currentTimeStepSum   = 2.0*sec.currentTimeStep;

  sec.lastTimeStep      = sec.currentTimeStep;
  sec.lastTimeStepRatio = sec.currentTimeStepRatio;
  sec.lastTimeStepSum   = sec.currentTimeStepSum;

  sec.numberSuccessiveFailures = 0;
  sec.stepAttemptStatus        = true;

  if (VERBOSE_TIME && tia_params.errorAnalysisOption == TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES)
  {
    Xyce::dout() << "ERROROPTION=1:  DeltaT Grow = 2" << "\n" << std::endl
                 << "ERROROPTION=1:  DeltaT Cut = 0.125" << "\n" << std::endl
                 << "ERROROPTION=1:  NL MIN = " << tia_params.NLmin << "\n" << std::endl
                 << "ERROROPTION=1:  NL MAX = " << tia_params.NLmax << "\n" << std::endl
                 << "ERROROPTION=1:  DELMAX = " << sec.maxTimeStep << "\n" << std::endl;
  }

  //  sec.tolAimFac_ = 0.5;

  sec.nextTime = sec.currentTime + sec.currentTimeStep;

  // x history
  *(ds.xHistory[0]) = *(ds.currSolutionPtr);
  (ds.xHistory[1])->putScalar(0.0); // no need to multiply by dt here

  // q history
  *(ds.qHistory[0]) = *(ds.daeQVectorPtr);
  (ds.qHistory[1])->linearCombo(1.0, *(ds.daeFVectorPtr), -1.0, *(ds.daeBVectorPtr));
  (ds.qHistory[1])->scale(-sec.currentTimeStep);

  // state history
  *(ds.sHistory[0]) = *(ds.currStatePtr);
  (ds.sHistory[1])->putScalar(0.0);

  // store history
  *(ds.stoHistory[0]) = *(ds.currStorePtr);
  (ds.stoHistory[1])->putScalar(0.0);

  if (ds.leadCurrentSize)
  {
    // new lead current Q component history
    *(ds.leadCurrentHistory[0]) = *(ds.currLeadCurrentPtr);
    (ds.leadCurrentHistory[1])->putScalar(0.0);

    *(ds.leadCurrentQHistory[0]) = *(ds.currLeadCurrentQPtr);
    (ds.leadCurrentQHistory[1])->putScalar(0.0);

    *(ds.leadDeltaVHistory[0]) = *(ds.currLeadDeltaVPtr);
    (ds.leadDeltaVHistory[1])->putScalar(0.0);
  }

  // Coefficient initialization
  sec.numberOfSteps_ = 0;    // number of total time integration steps taken
  sec.currentOrder_ = 1;
  sec.usedOrder_ = 1;
  sec.psi_[0] = sec.currentTimeStep;
  sec.cj_ = 1/sec.psi_[0];
  sec.nscsco_ = 0;
  if (DEBUG_TIME && isActive(Diag::TIME_HISTORY))
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() <<
      "  OneStep::initialize" << std::endl;
    Xyce::dout() << "\n xHistory: \n" << std::endl;
    (ds.xHistory[0])->print(Xyce::dout());
    Xyce::dout() << std::endl;
    (ds.xHistory[1])->print(Xyce::dout());
    Xyce::dout() << std::endl;
    Xyce::dout() << "\n qHistory: \n" << std::endl;
    (ds.qHistory[0])->print(Xyce::dout());
    Xyce::dout() << std::endl;
    (ds.qHistory[1])->print(Xyce::dout());
    Xyce::dout() << std::endl;
    Xyce::dout() << "\n sHistory: \n" << std::endl;
    (ds.sHistory[0])->print(Xyce::dout());
    Xyce::dout() << std::endl;
    (ds.sHistory[1])->print(Xyce::dout());
    Xyce::dout() << std::endl;
    Xyce::dout() << "\n" << "currentTimeStep = " << currentTimeStep << "\n" << std::endl;
    Xyce::dout() << "\n" << "time_to_stop = " << time_to_stop << "\n" << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  initializeSensitivities();
}


//-----------------------------------------------------------------------------
// Function      : OneStep::initializeAdjoint 
// Purpose       : Initialize method for adjoint sensitivities
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/28/2016
//-----------------------------------------------------------------------------
void OneStep::initializeAdjoint (int index)
{
  ds.nextDQdxLambdaPtr->putScalar(0.0);
  ds.currDQdxLambdaPtr->putScalar(0.0);
  ds.lastDQdxLambdaPtr->putScalar(0.0);

  ds.nextDFdxLambdaPtr->putScalar(0.0);
  ds.currDFdxLambdaPtr->putScalar(0.0);
  ds.lastDFdxLambdaPtr->putScalar(0.0);

  ds.nextLambdaPtr->putScalar(0.0);
  ds.currLambdaPtr->putScalar(0.0);
  ds.lastLambdaPtr->putScalar(0.0);

  sec.currentTime = ds.timeHistory[index];
  sec.nextTime = sec.currentTime;

  sec.currentTimeStep = ds.dtHistory[index];

  // not right yet but testing:
  sec.lastTimeStep = ds.dtHistory[index];
  sec.oldeTimeStep = ds.dtHistory[index];

  sec.lastTime    = sec.currentTime;
}

//-----------------------------------------------------------------------------
// Function      : OneStep::initializeSensitivities
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void OneStep::initializeSensitivities()
{
  if (ds.numParams)
  {
    // dqdp history
    *(ds.dqdpHistory[0]) = *ds.currDqdpPtrVector;
    *(ds.dqdpHistory[1]) = *ds.currDqdpPtrVector;
  }
}

//-----------------------------------------------------------------------------
// Function      : OneStep::setTwoLevelTimeInfo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/01/07
//-----------------------------------------------------------------------------
void OneStep::setTwoLevelTimeInfo()
{
  // x history
  *ds.xHistory[0] = *ds.currSolutionPtr;
  ds.xHistory[1]->putScalar(0.0); // no need to multiply by dt here

  // q history
  *ds.qHistory[0] = *ds.daeQVectorPtr;
  ds.qHistory[1]->linearCombo(1.0, *ds.daeFVectorPtr, -1.0, *ds.daeBVectorPtr);
  ds.qHistory[1]->scale(-sec.currentTimeStep);

  // state history
  *ds.sHistory[0] = *ds.nextStatePtr;
  ds.sHistory[1]->putScalar(0.0);

  // Coefficient initialization
  sec.numberOfSteps_ = 0;    // number of total time integration steps taken
  sec.usedOrder_ = 1;
  sec.psi_[0] = sec.currentTimeStep;
  sec.cj_ = 1/sec.psi_[0];
  sec.nscsco_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : OneStep::checkReduceOrder()
// Purpose       : check whether to reduce order independent of local error test
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
void OneStep::checkReduceOrder()
{

}

//-----------------------------------------------------------------------------
// Function      : OneStep::rejectStep()
// Purpose       : code to restore history, choose new order/step-size
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
void OneStep::rejectStep(const TIAParams & tia_params)
{
  // This routine puts its output   newTimeStep_
  sec.TimeStepLimitedbyBP = false;

  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() << "  OneStep::rejectStep" << std::endl;
  }

  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = !tia_params.constantTimeStepFlag;

  sec.lastAttemptedTimeStep = sec.currentTimeStep;

  double newTimeStep_ = sec.currentTimeStep;
  double rr = 1.0; // step size ratio = new step / old step
  if ((sec.stepAttemptStatus == false) && (adjustStep))
  {
    if (tia_params.errorAnalysisOption == TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES)
    {
      newTimeStep_ = sec.currentTimeStep/8;
    }
    else
    {
      sec.initialPhase_ = false;
      sec.nef_++;
      restoreHistory();
      if (sec.nef_ >= sec.max_LET_fail_)
      {
        Xyce::Report::DevelFatal0().in("OneStep::rejectStep")
          << "  Maximum number of local error test failures.  ";
      }

      if ((sec.newtonConvergenceStatus <= 0))
      {
        /// 11/11/05 erkeite:  If the Newton solver fails, don't
        // rely on the error estimate - it may be full of Nan's.
        //        rr = sec.r_min_;

        newTimeStep_ = sec.currentTimeStep/8;
        sec.currentOrder_ = sec.minOrder_;
      }
      else
      {
        // 03/11/04 tscoffe:  Here is the block for choosing order &
        // step-size when the local error test FAILS (but Newton
        // succeeded).
        if (sec.nef_== 1)
        {

          //	sec.estOverTol_
          sec.Est_ = sec.estOverTol_;

          sec.currentOrder_ = sec.minOrder_;
          rr = sec.tolAimFac_/(sec.estOverTol_ + 0.0001);
          rr = pow(rr, 1.0/(sec.currentOrder_+1.0));
          rr = std::max(sec.r_min_,std::min(sec.r_max_,rr));

          newTimeStep_ = rr * sec.currentTimeStep;
        }
        else
        {
          rr = sec.r_min_;
          newTimeStep_ = rr * sec.currentTimeStep;

          sec.currentOrder_ = sec.minOrder_;
        }
      }

      if (DEBUG_TIME && isActive(Diag::TIME_STEP))
      {
        Xyce::dout() << "  currentTimeStep = " <<  sec.currentTimeStep << std::endl;
        Xyce::dout() << "  numberOfSteps_ = " <<  sec.numberOfSteps_ << std::endl;
        Xyce::dout() << "  currentOrder_ = " <<  sec.currentOrder_ << std::endl;
        Xyce::dout() << "  nscsco_ = " <<  sec.nscsco_ << std::endl;
        Xyce::dout() << "  alpha_[0] = " <<  sec.alpha_[0] << std::endl;
        Xyce::dout() << "  alpha_[1] = " <<  sec.alpha_[1] << std::endl;
        Xyce::dout() << "  alpha_[2] = " <<  sec.alpha_[2] << std::endl;
        Xyce::dout() << "  alpha_[3] = " <<  sec.alpha_[3] << std::endl;
        Xyce::dout() << "  alpha_[4] = " <<  sec.alpha_[4] << std::endl;
        Xyce::dout() << "  psi_[0] = " <<  sec.psi_[0] << std::endl;
        Xyce::dout() << "  psi_[1] = " <<  sec.psi_[1] << std::endl;
        Xyce::dout() << "  psi_[2] = " <<  sec.psi_[2] << std::endl;
        Xyce::dout() << "  psi_[3] = " <<  sec.psi_[3] << std::endl;
        Xyce::dout() << "  psi_[4] = " <<  sec.psi_[4] << std::endl;
        Xyce::dout() << "  sigma_[0] = " <<  sec.sigma_[0] << std::endl;
        Xyce::dout() << "  sigma_[1] = " <<  sec.sigma_[1] << std::endl;
        Xyce::dout() << "  sigma_[2] = " <<  sec.sigma_[2] << std::endl;
        Xyce::dout() << "  sigma_[3] = " <<  sec.sigma_[3] << std::endl;
        Xyce::dout() << "  rr = " <<  rr << std::endl;
        Xyce::dout() << "  Est_ = " <<  sec.Est_ << std::endl;
        Xyce::dout() << "  nef_ = " <<  sec.nef_ << std::endl;
        Xyce::dout() << "  newOrder_ = " <<  sec.newOrder_ << std::endl;
        Xyce::dout() << "  currentTimeStep = " <<  sec.currentTimeStep << std::endl;
        Xyce::dout() << "  newTimeStep_ = " <<  newTimeStep_ << std::endl;
      }
    }
  }
  else if ((sec.stepAttemptStatus == false) & (!adjustStep))
  {
    std::string tmp = "  OneStep:rejectStep: Warning: Local error test failed with constant step-size.\n";
    Xyce::dout() << tmp << std::endl;
  }

  // If the step needs to be adjusted:
  if (adjustStep)
  {
    newTimeStep_ = std::max(newTimeStep_, sec.minTimeStep);
    newTimeStep_ = std::min(newTimeStep_, sec.maxTimeStep);

    double nextTimePt = sec.currentTime + newTimeStep_;

    if (nextTimePt > sec.stopTime)
    {
      nextTimePt  = sec.stopTime;
      newTimeStep_ = sec.stopTime - sec.currentTime;
      sec.TimeStepLimitedbyBP = true;
    }

    sec.nextTime = nextTimePt;

    sec.currentTimeStepRatio = newTimeStep_/sec.lastTimeStep;
    sec.currentTimeStepSum   = newTimeStep_ + sec.lastTimeStep;

    if (DEBUG_TIME && isActive(Diag::TIME_STEP))
    {
      Xyce::dout() << "  newTimeStep_ = " <<  newTimeStep_ << std::endl;
      Xyce::dout() << "  nextTime = " <<  sec.nextTime << std::endl;
    }

    sec.currentTimeStep = newTimeStep_;
  }
  else // if time step is constant for this step:
  {
    double nextTimePt = sec.currentTime + sec.currentTimeStep;

    if (nextTimePt > sec.stopTime)
    {
      nextTimePt      = sec.stopTime;
      sec.currentTimeStep = sec.stopTime - sec.currentTime;
    }

    sec.currentTimeStepRatio = sec.currentTimeStep / sec.lastTimeStep;
    sec.currentTimeStepSum   = sec.currentTimeStep + sec.lastTimeStep;

    sec.nextTime = nextTimePt;
  }
  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : OneStep::rejectStepForHabanero
// Purpose       : step rejection, but from an outside program (Habanero API)
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 08/11/09
//-----------------------------------------------------------------------------
void OneStep::rejectStepForHabanero()
{
  restoreHistory();
  sec.setTimeStep(sec.currentTimeStep);
}

//-----------------------------------------------------------------------------
// Function      : OneStep::completeStep()
// Purpose       : code to update history, choose new order/step-size
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
void OneStep::completeStep(const TIAParams &tia_params)
{
  sec.TimeStepLimitedbyBP = false;

  sec.numberOfSteps_ ++;
  sec.nef_ = 0;
  sec.lastTime    = sec.currentTime;
  sec.currentTime = sec.nextTime;

  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;

    Xyce::dout() << "  OneStep::completeStep" << std::endl;
  }

  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = !tia_params.constantTimeStepFlag;

  sec.lastAttemptedTimeStep = sec.currentTimeStep;

  double newTimeStep_ = sec.currentTimeStep;
  double rr = 1.0; // step size ratio = new step / old step
  // 03/11/04 tscoffe:  Here is the block for choosing order & step-size when
  // the local error test PASSES (and Newton succeeded).

  // need the lastTimeStep for output interpolation.  Save it before
  // overwriting it with currentTimeStep;
  // This history for this step is saved in updateHistory(), called below
  timeStepForHistory2_=sec.lastTimeStep;

  sec.lastTimeStep = sec.currentTimeStep;
  sec.lastTimeStepRatio = sec.currentTimeStepRatio; // copied from calcTStep1
  sec.lastTimeStepSum   = sec.currentTimeStepSum; // copied from calcTStep1
  sec.usedOrder_ = sec.currentOrder_;
  sec.usedStep_ = sec.currentTimeStep;

  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << "  initialPhase_ = " <<  sec.initialPhase_ << std::endl;
    Xyce::dout() << "  rr = " <<  rr << std::endl;
    Xyce::dout() << "  currentTimeStep = " <<  sec.currentTimeStep << std::endl;
    Xyce::dout() << "  currentTime = " <<  sec.currentTime << std::endl;
    Xyce::dout() << "  nextTime = " <<  sec.nextTime << std::endl;
    Xyce::dout() << "  newTimeStep_ = " <<  newTimeStep_ << std::endl;
    Xyce::dout() << "  minTimeStep = " <<  sec.minTimeStep << std::endl;
    Xyce::dout() << "  maxTimeStep = " <<  sec.maxTimeStep << std::endl;
  }

  if (tia_params.errorAnalysisOption == TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES)
  {
    if (sec.numberOfSteps_ >= 2 &&  sec.maxOrder_ == 2)
    {
      sec.currentOrder_ = 2;
    }

    rr = 1;

    if (sec.nIterations <= tia_params.NLmin)
      rr = 2;

    if (sec.nIterations > tia_params.NLmax)
      rr = 1.0/8;

    newTimeStep_ = rr*sec.currentTimeStep;
  }
  else
  {
    rr = sec.tolAimFac_/(sec.estOverTol_ + 0.0001);
    rr = pow(rr, 1.0/(sec.currentOrder_+1.0));

    if (sec.numberOfSteps_ >= 2 && sec.maxOrder_ == 2)
    {
      if (sec.currentOrder_ == 1)
      {
        sec.currentOrder_ = 2;
        rr = sec.tolAimFac_/(sec.estOverTol_ + 0.0001);
        rr = pow(rr, 1.0/(sec.currentOrder_+1.0));

        if (rr <= 1.05)
        {

          sec.currentOrder_ = sec.minOrder_;
        }
      }
    }
    if (DEBUG_TIME && isActive(Diag::TIME_STEP))
    {
      Xyce::dout() << "  currentOrder_ = " <<  sec.currentOrder_ << std::endl;
      Xyce::dout() << "  r_hincr_ = " <<  sec.r_hincr_ << std::endl;
      Xyce::dout() << "  r_hincr_test_ = " <<  sec.r_hincr_test_ << std::endl;
      Xyce::dout() << "  Est = " <<  sec.Est_ << std::endl;
      Xyce::dout() << "  raw rr = " <<  rr << std::endl;
    }

    if (rr >= sec.r_hincr_test_)
    {
      rr = sec.r_hincr_;
      newTimeStep_ = rr*sec.currentTimeStep;
    }
    else if (rr <= 1)
    {
      rr = std::max(sec.r_min_,std::min(sec.r_max_,rr));
      newTimeStep_ = rr*sec.currentTimeStep;
    }
  }

  updateHistory();

  newTimeStep_ = std::max(newTimeStep_, sec.minTimeStep);
  newTimeStep_ = std::min(newTimeStep_, sec.maxTimeStep);

  // 12/01/05 tscoffe:  This is necessary to avoid currentTimeStep == 0 right
  // before a breakpoint.  So I'm checking to see if currentTime is identically
  // equal to stopTime, in which case we are right before a breakpoint and we
  // should not adjust currentStepSize because that would result in
  // currentStepSize == 0.
  if ((sec.stopTime - sec.currentTime) >= sec.minTimeStep)
  {
    // If the step needs to be adjusted:
    if (adjustStep)
    {
      double nextTimePt = sec.currentTime + newTimeStep_;
      if (nextTimePt > sec.stopTime)
      {

        sec.savedTimeStep = newTimeStep_;

        nextTimePt  = sec.stopTime;

        newTimeStep_ = sec.stopTime - sec.currentTime;
        sec.TimeStepLimitedbyBP = true;
      }

      sec.nextTime = nextTimePt;

      sec.currentTimeStepRatio = newTimeStep_/sec.lastTimeStep;
      sec.currentTimeStepSum   = newTimeStep_ + sec.lastTimeStep;

      if (DEBUG_TIME && isActive(Diag::TIME_STEP))
      {
        Xyce::dout() << "  nextTime = " <<  sec.nextTime << std::endl;
        Xyce::dout() << "  newTimeStep_ = " <<  newTimeStep_ << std::endl;
      }

    sec.currentTimeStep = newTimeStep_;
    }
    else // if time step is constant for this step:
    {
      double nextTimePt = sec.currentTime + sec.currentTimeStep;

      if (nextTimePt > sec.stopTime)
      {
        sec.savedTimeStep = sec.currentTimeStep;

        nextTimePt      = sec.stopTime;
        sec.currentTimeStep = sec.stopTime - sec.currentTime;
      }

      sec.currentTimeStepRatio = sec.currentTimeStep / sec.lastTimeStep;
      sec.currentTimeStepSum   = sec.currentTimeStep + sec.lastTimeStep;

      sec.nextTime = nextTimePt;
    }
  }

  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

}


//-----------------------------------------------------------------------------
// Function      : OneStep::completeAdjointStep()
// Purpose       : code to update history, choose new order/step-size for adjoint
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void OneStep::completeAdjointStep(const TIAParams &tia_params)
{
  sec.lastTime    = sec.currentTime;

  int timeSize = ds.timeHistory.size();

  if (ds.itAdjointIndex < timeSize)
  {
    sec.currentTime = ds.timeHistory[ds.itAdjointIndex];
  }
  else
  {
    sec.currentTime = ds.timeHistory[timeSize-1];
  }

  if (ds.itAdjointIndex > 0)
  {
    sec.nextTime = ds.timeHistory[ds.itAdjointIndex-1];
  }
  else
  {
    sec.nextTime = ds.timeHistory[0];
  }
    
  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::completeAdjointStep" << std::endl;
  }

  sec.lastAttemptedTimeStep = sec.currentTimeStep;

  sec.oldeTimeStep = sec.lastTimeStep;
  sec.lastTimeStep = sec.currentTimeStep;

  if (ds.itAdjointIndex>0)
  {
    sec.currentTimeStep = ds.dtHistory[ds.itAdjointIndex-1];
  }
  else
  {
    sec.currentTimeStep = ds.dtHistory[ds.itAdjointIndex];
  }

  sec.usedOrder_ = sec.currentOrder_;
  sec.usedStep_ = sec.currentTimeStep;

  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << "  initialPhase_ = " <<  sec.initialPhase_ << std::endl
      << "  currentTimeStep = " <<  sec.currentTimeStep << std::endl
      << "  currentTime = " <<  sec.currentTime << std::endl
      << "  nextTime = " <<  sec.nextTime << std::endl
      << "  minTimeStep = " <<  sec.minTimeStep << std::endl
      << "  maxTimeStep = " <<  sec.maxTimeStep << std::endl;
  }

  updateAdjointSensitivityHistory();

  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : OneStep::updateStateDeriv
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/17/05
//-----------------------------------------------------------------------------
void OneStep::updateStateDeriv ()
{
  // dS/dt = spn0 - (sec.alpha_/hn)(S(x)-sn0)
  ds.nextStateDerivPtr->linearCombo(1.0,*ds.nextStatePtr,-1.0,*(ds.sHistory[0]));

  if (sec.currentOrder_ == 1)
  {
    ds.nextStateDerivPtr->scale(1.0/sec.currentTimeStep);
  }
  else
  {
    ds.nextStateDerivPtr->update(-1.0,*(ds.currStateDerivPtr), 2.0/sec.currentTimeStep);
  }

  if (DEBUG_TIME && isActive(Diag::TIME_DUMP_SOLUTION_ARRAYS))
  {
    Xyce::dout() << "\n Next State Ptr: \n" << std::endl;
    ds.nextStatePtr->print(Xyce::dout());
    Xyce::dout() << std::endl;

    Xyce::dout() << "\n Next State Deriv: \n" << std::endl;
    ds.nextStateDerivPtr->print(Xyce::dout());
    Xyce::dout() << std::endl;
  }
}


//-----------------------------------------------------------------------------
// Function      : OneStep::updateLeadCurrentVec
// Purpose       : calculates lead currents in lead current vector with
//                 the leadCurrentQVec.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling, SNL
// Creation Date : 01/22/2015
//-----------------------------------------------------------------------------
void OneStep::updateLeadCurrentVec ()
{
  if (ds.leadCurrentSize)
  {
    ds.nextLeadCurrentQDerivPtr->linearCombo(1.0,*ds.nextLeadCurrentQPtr,-1.0,*ds.leadCurrentQHistory[0]);
    if (sec.currentOrder_ == 1)
    {
      ds.nextLeadCurrentQDerivPtr->scale(1.0/sec.currentTimeStep);
    }
    else
    {
      ds.nextLeadCurrentQDerivPtr->
        update(-1.0,*(ds.currLeadCurrentQDerivPtr),2.0/sec.currentTimeStep);
    }
    ds.nextLeadCurrentPtr->update(1.0,*ds.nextLeadCurrentQDerivPtr);

    if (DEBUG_TIME && isActive(Diag::TIME_STEP))
    {
      Xyce::dout() << "\n next lead current Ptr: \n" << std::endl;
      ds.nextLeadCurrentPtr->print(Xyce::dout());
      Xyce::dout() << std::endl;

      Xyce::dout() << "\n next lead current Q Ptr: \n" << std::endl;
      ds.nextLeadCurrentQPtr->print(Xyce::dout());
      Xyce::dout() << std::endl;

      Xyce::dout() << "\n curr lead current Q Ptr: \n" << std::endl;
      ds.leadCurrentQHistory[0]->print(Xyce::dout());
      Xyce::dout() << std::endl;
    }
  }
}
//-----------------------------------------------------------------------------
// Function      : OneStep::getInitialQnorm
// Purpose       : Needed by 2-level solves.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/18/07
//-----------------------------------------------------------------------------
void OneStep::getInitialQnorm(TwoLevelError & tle) const
{
  tle.q1HistorySum = ds.partialSum_q1();
}

//-----------------------------------------------------------------------------
// Function      : OneStep::setupTwoLevelError
// Purpose       : Needed by 2-level solves.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
void OneStep::getTwoLevelError(TwoLevelError & tle) const
{
  tle.xErrorSum    = ds.partialErrorNormSum ();
  tle.qErrorSum    = ds.partialQErrorNormSum ();
  tle.xErrorSum_m1 = ds.partialSum_m1 (sec.currentOrder_);
  tle.xErrorSum_p1 = ds.partialSum_p1 (sec.currentOrder_, sec.maxOrder_);
  tle.innerSize    = ds.newtonCorrectionPtr->globalLength();

  if (DEBUG_TIME && isActive(Diag::TIME_ERROR))
  {
    Xyce::dout() << tle;
  }
}

//-----------------------------------------------------------------------------
// Function      : OneStep::partialTimeDeriv
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 1/20/07
//-----------------------------------------------------------------------------
double OneStep::partialTimeDeriv() const
{
  if (sec.currentTimeStep < 1e-30)
  {
    Xyce::Report::UserWarning() << "Excessively small current time step, incorrectly returning with large value";

    return leadingCoeff * 1.e+30;
  }

  return leadingCoeff / sec.currentTimeStep;
}

//-----------------------------------------------------------------------------
// Function      : OneStep::getSolnVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/24/2016
//-----------------------------------------------------------------------------
bool OneStep::getSolnVarData( const int & gid, std::vector<double> & varData )
{
  // Dump DataStore information.
  int num = ds.getNumSolnVarData();
  bool ret = ds.getSolnVarData( gid, varData );

  // Determine order for OneStep dumping.
  int order = sec.currentOrder_;

  if (ret)
  {
    varData.resize( num + 2*(order+1) );
    for (int i=0; i <= order; ++i )
    {
      varData[num++] = ds.xHistory[i]->getElementByGlobalIndex ( gid );
      varData[num++] = ds.qHistory[i]->getElementByGlobalIndex ( gid );
    }
  }
  return ret;
}

//-----------------------------------------------------------------------------
// Function      : OneStep::setSolnVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/24/2016
//-----------------------------------------------------------------------------
bool OneStep::setSolnVarData( const int & gid, const std::vector<double> & varData )
{
  // Restore data store information
  int num = ds.getNumSolnVarData();
  bool ret = ds.setSolnVarData( gid, varData );

  // Determine order for OneStep restoring.
  int order = sec.currentOrder_;


//  if (sec.maxOrder_ < sec.currentOrder_)
//    order = sec.maxOrder_;
  if (ret)
  {
    for (int i=0; i <= order; ++i )
    {
      ds.xHistory[i]->setElementByGlobalIndex( gid, varData[num++] );
      ds.qHistory[i]->setElementByGlobalIndex( gid, varData[num++] );
    }
  }
  return ret;
}

//-----------------------------------------------------------------------------
// Function      : OneStep::getStateVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/24/2016
//-----------------------------------------------------------------------------
bool OneStep::getStateVarData( const int & gid, std::vector<double> & varData )
{
  // Dump DataStore information.
  int num = ds.getNumStateVarData();
  bool ret = ds.getStateVarData( gid, varData );

  if (ret)
  {
    varData.resize( num + 2 );
    varData[num++] = ds.sHistory[0]->getElementByGlobalIndex ( gid );
    varData[num++] = ds.sHistory[1]->getElementByGlobalIndex ( gid );
  }
  return ret;
}

//-----------------------------------------------------------------------------
// Function      : OneStep::setStateVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/24/2016
//-----------------------------------------------------------------------------
bool OneStep::setStateVarData( const int & gid, const std::vector<double> & varData )
{
  // Restore data store information
  int num = ds.getNumStateVarData();
  bool ret = ds.setStateVarData( gid, varData );

  if (ret)
  {
    ds.sHistory[0]->setElementByGlobalIndex( gid, varData[num++] );
    ds.sHistory[1]->setElementByGlobalIndex( gid, varData[num++] );
  }
  return ret;
}

//-----------------------------------------------------------------------------
// Function      : OneStep::getStoreVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/24/2016
//-----------------------------------------------------------------------------
bool OneStep::getStoreVarData( const int & gid, std::vector<double> & varData )
{
  // Dump DataStore information.
  int num = ds.getNumStoreVarData();
  bool ret = ds.getStoreVarData( gid, varData );

  if (ret)
  {
    varData.resize( num + 2 );
    varData[num++] = ds.stoHistory[0]->getElementByGlobalIndex ( gid );
    varData[num++] = ds.stoHistory[1]->getElementByGlobalIndex ( gid );
  }
  return ret;
}

//-----------------------------------------------------------------------------
// Function      : OneStep::setStoreVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/24/2016
//-----------------------------------------------------------------------------
bool OneStep::setStoreVarData( const int & gid, const std::vector<double> & varData )
{
  // Restore data store information
  int num = ds.getNumStoreVarData();
  bool ret = ds.setStoreVarData( gid, varData );

  if (ret)
  {
    ds.stoHistory[0]->setElementByGlobalIndex( gid, varData[num++] );
    ds.stoHistory[1]->setElementByGlobalIndex( gid, varData[num++] );
  }
  return ret;
}

} // namespace TimeIntg
} // namespace Xyce

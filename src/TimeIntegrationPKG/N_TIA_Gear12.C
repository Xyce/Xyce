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
// Purpose       : This file contains the functions which define the
//		             backward differentiation, order 1-2, class.
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
#include <N_TIA_Gear12.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TIAParams.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_MachDepParams.h>

using std::abs;

namespace Xyce {
namespace TimeIntg {

const char *
Gear12::name = "Gear 12";

//-----------------------------------------------------------------------------
// Function      : Gear12::factory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :  10/31/07 
//-----------------------------------------------------------------------------
TimeIntegrationMethod * Gear12::factory(
    const TIAParams &   tia_params,
    StepErrorControl &  step_error_control,
    DataStore &         data_store)
{
  return new Gear12(tia_params, step_error_control, data_store);
}

//-----------------------------------------------------------------------------
// Function      : Gear::Gear
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/10/12
//-----------------------------------------------------------------------------
Gear12::Gear12(
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


  sec.currentOrder_ = (std::min(sec.currentOrder_, sec.maxOrder_) );

}

//-----------------------------------------------------------------------------
// Function      : Gear12::obtainPredictor
// Purpose       : Calculate predictor 
// Special Notes : stored in ds.xn0Ptr,qn0Ptr
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/10/12
//-----------------------------------------------------------------------------
void Gear12::obtainPredictor()
{ 
  // evaluate predictor
  ds.xn0Ptr->putScalar(0.0);
  ds.qn0Ptr->putScalar(0.0);

  for (int i=0;i<=sec.currentOrder_;++i)
  {
    ds.xn0Ptr->update(sec.beta_[i],*(ds.xHistory[i]));
    ds.qn0Ptr->update(sec.beta_[i],*(ds.qHistory[i]));
  }

  if (DEBUG_TIME && isActive(Diag::TIME_PREDICTOR))
  {
    Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::obtainPredictor" << std::endl
      << "\n currentOrder = " << sec.currentOrder_ << std::endl
      << "\n sec.nscsco_: " << sec.nscsco_ << std::endl;
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
    for (int i=0; i<=sec.currentOrder_ ; ++i)
    {
      Xyce::dout() << "\n sHistory["<< i << "]: \n" << std::endl;
      (ds.sHistory[i])->print(Xyce::dout());
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
// Function      : Gear12::obtainSensitivityPredictors
// Purpose       : Calculate predictor 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void Gear12::obtainSensitivityPredictors()
{}


//-----------------------------------------------------------------------------
// Function      : Gear12::obtainResidual
// Purpose       : Calculate Residual
// Special Notes : 
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------
void Gear12::obtainResidual()
{
  // output: ds.RHSVectorPtr
  // Note:  ds.nextSolutionPtr is used to get Q,F,B in Analysis::AnalysisManager::loadRHS.
  ds.RHSVectorPtr->update(sec.alpha_[0],*ds.daeQVectorPtr, sec.alpha_[1],*(ds.qHistory[0]),0.0);

  if (DEBUG_TIME && isActive(Diag::TIME_RESIDUAL))
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::obtainResidual" << std::endl
      << "\n t = " << sec.nextTime << "\n" << std::endl
      << "\n solution: \n" << std::endl;
    ds.nextSolutionPtr->print(Xyce::dout());
    Xyce::dout() << "\n daeQVector: \n" << std::endl;
    ds.daeQVectorPtr->print(Xyce::dout());
    Xyce::dout() << "\n qn0: \n" << std::endl;
    ds.qn0Ptr->print(Xyce::dout());
    Xyce::dout() << "\n sec.alphas_/hn: " << sec.alphas_/sec.currentTimeStep << "\n" << std::endl
      << "\n daeFVector: \n" << std::endl;
    ds.daeFVectorPtr->print(Xyce::dout());

    Xyce::dout() << "\n dQdt-vector: \n" << std::endl;
    ds.RHSVectorPtr->print(Xyce::dout());
    Xyce::dout() << std::endl;
  }

  if (sec.currentOrder_  == 2)
  {
    ds.RHSVectorPtr->update(sec.alpha_[2],*(ds.qHistory[1]));
  }

  ds.RHSVectorPtr->update(+1.0,*ds.daeFVectorPtr,-1.0,*ds.daeBVectorPtr,1.0/sec.currentTimeStep);

  // since the nonlinear solver is expecting a -f, scale by -1.0:
  ds.RHSVectorPtr->scale(-1.0);

  // if voltage limiting is on, add it in:
  if (ds.limiterFlag)
  {
    (ds.dQdxdVpVectorPtr)->scale( sec.alpha_[0]/sec.currentTimeStep );
    //        double qscalar(sec.alpha_[0]/sec.currentTimeStep);

    (ds.RHSVectorPtr)->update(+1.0, *(ds.dQdxdVpVectorPtr));

    (ds.RHSVectorPtr)->update(+1.0, *(ds.dFdxdVpVectorPtr));
  }

  if (DEBUG_TIME && isActive(Diag::TIME_RESIDUAL))
  {
    Xyce::dout() << "\n Residual-vector: \n" << std::endl;
    ds.RHSVectorPtr->print(Xyce::dout());
    Xyce::dout() << Xyce::section_divider << std::endl
      << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Gear12::obtainSensitivityResiduals 
// Purpose       : Calculate sensitivity residual
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Gear12::obtainSensitivityResiduals()
{
  ds.sensRHSPtrVector->update(sec.alpha_[0], *(ds.nextDqdpPtrVector),
                              sec.alpha_[1], *(ds.dqdpHistory[0]), 0.0);

  if (sec.currentOrder_  == 2)
  {
    ds.sensRHSPtrVector->update(sec.alpha_[2],*(ds.dqdpHistory[1]));
  }

  ds.sensRHSPtrVector->update(+1.0, *(ds.nextDfdpPtrVector),
                              -1.0, *(ds.nextDbdpPtrVector),
                              +1.0/sec.currentTimeStep);

  // since the nonlinear solver is expecting a -f, scale by -1.0:
  ds.sensRHSPtrVector->scale(-1.0);

  double qscalar1(sec.alpha_[1]/sec.currentTimeStep);
  ds.sensRHSPtrVector->update(-qscalar1, *ds.currDQdxDXdpPtrVector);

  if (sec.currentOrder_  == 2)
  {
    double qscalar2(sec.alpha_[2]/sec.currentTimeStep);
    ds.sensRHSPtrVector->update(-qscalar2, *ds.lastDQdxDXdpPtrVector);
  }

#ifdef DEBUG_SENS
  Xyce::dout() << "Gear12: obtainSensitivityResiduals: RHS Vector : " << std::endl;
  ds.sensRHSPtrVector->print(Xyce::dout());
#endif
}


//-----------------------------------------------------------------------------
// Function      : Gear12::obtainFunctionDerivativesForTranAdjoint 
// Purpose       : 
//
// Special Notes : This is similar to the direct-sensitivity method, 
//                 Gear12::obtainSensitivityResiduals .  The main difference is
//                 that this function doesn't include the so-called cross terms
//                 (or chain rule terms).  It only loads the function (device) 
//                 parameter derivatives.
//
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 3/13/2016
//-----------------------------------------------------------------------------
void Gear12::obtainFunctionDerivativesForTranAdjoint()
{
  ds.sensRHSPtrVector->update(sec.alpha_[0], *(ds.nextDqdpPtrVector),
                              sec.alpha_[1], *(ds.dqdpHistory[0]), 0.0);

  if (sec.currentOrder_  == 2)
  {
    ds.sensRHSPtrVector->update(sec.alpha_[2],*(ds.dqdpHistory[1]));
  }

  ds.sensRHSPtrVector->update(+1.0, *(ds.nextDfdpPtrVector),
                              -1.0, *(ds.nextDbdpPtrVector),
                              1.0/sec.currentTimeStep);

  // since the nonlinear solver is expecting a -f, scale by -1.0:
  ds.sensRHSPtrVector->scale(-1.0);
}

//-----------------------------------------------------------------------------
// Function      : Gear12::obtainSparseFunctionDerivativesForTranAdjoint 
// Purpose       : 
//
// Special Notes : This is similar to the direct-sensitivity method, 
//                 Gear12::obtainSensitivityResiduals .  The main difference is
//                 that this function doesn't include the so-called cross terms
//                 (or chain rule terms).  It only loads the function (device) 
//                 parameter derivatives.
//
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 1/9/2017
//-----------------------------------------------------------------------------
void Gear12::obtainSparseFunctionDerivativesForTranAdjoint()
{
  ds.sensRHSPtrVector->update(sec.alpha_[0], *(ds.nextDqdpPtrVector),
                              sec.alpha_[1], *(ds.dqdpHistory[0]), 0.0);

  if (sec.currentOrder_  == 2)
  {
    ds.sensRHSPtrVector->update(sec.alpha_[2],*(ds.dqdpHistory[1]));
  }

  ds.sensRHSPtrVector->update(+1.0, *(ds.nextDfdpPtrVector),
                              -1.0, *(ds.nextDbdpPtrVector),
                              1.0/sec.currentTimeStep);

  // since the nonlinear solver is expecting a -f, scale by -1.0:
  ds.sensRHSPtrVector->scale(-1.0);

  // Collect sparse information from sensRHSPtrVector.
  ds.sparseSensRHSMV->replaceValues( ds.masterIndexVector, *ds.sensRHSPtrVector );
}

//-----------------------------------------------------------------------------
// Function      : Gear12::obtainAdjointSensitivityResidual 
// Purpose       : Calculate adjoint sensitivity residual.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 2/23/2015
//-----------------------------------------------------------------------------
void Gear12::obtainAdjointSensitivityResidual()
{
  Linear::Vector & RHSVec  = *(ds.RHSVectorPtr);

  Linear::Vector & currDQdxLambda = *(ds.currDQdxLambdaPtr);
  Linear::Vector & lastDQdxLambda = *(ds.lastDQdxLambdaPtr);

  Linear::Vector & currLambda = *(ds.currLambdaPtr);
  Linear::Vector & lastLambda = *(ds.lastLambdaPtr);

  Linear::Matrix & dQdx = *(ds.dQdxMatrixPtr);

  // This assumes that the matvecs have been properly stored.
  // The BDF2 coefficients (alpha0, alpha1 and alpha2) are not constant, 
  // and can change from step to step in the forward calculation.  For 
  // this version of the adjoint solve, it is necessary to compute the 
  // coefficients.

  int it=ds.itAdjointIndex;
  int itmax=ds.timeHistory.size();

  if (it<itmax-1) 
  {
    double alpha2=0.0;
    double alpha1=0.0;
    double alpha0=0.0;

      //orderHistory;
    if (ds.orderHistory[it+1] == 1)
    {
      alpha0 = 1.00000000e+00;
      alpha1 = -1.00000000e+00;
    }
    else
    {
      double dt1 = ds.dtHistory[it+1];
      double dtOld1 = ds.dtHistory[it];

      alpha2 = -dt1/dtOld1 * dt1/(2 * dt1 + dtOld1); 
      alpha1 = 1 - alpha2;
      alpha0 = -alpha1 - alpha2 * (1 + dtOld1/dt1);

      alpha2 = alpha2/alpha0;
      alpha1 = alpha1/alpha0;
      alpha0 = -1/alpha0;
    }

    double qscalar1(alpha1/sec.lastTimeStep);
    bool Transpose = true;
    currDQdxLambda.putScalar(0.0);
    dQdx.matvec( Transpose , currLambda, currDQdxLambda);
    currDQdxLambda.scale(-1);
    RHSVec.update(+qscalar1, currDQdxLambda);
  }

  if (it<itmax-2)
  {
    if (ds.orderHistory[it+2] == 1)
    {
      // don't do anything.
    }
    else
    {
      double alpha2=0.0;
      double alpha1=0.0;
      double alpha0=0.0;

      double dt1 = ds.dtHistory[it+2];
      double dtOld1 = ds.dtHistory[it+1];

      alpha2 = -dt1/dtOld1 * dt1/(2 * dt1 + dtOld1); 
      alpha1 = 1 - alpha2;
      alpha0 = -alpha1 - alpha2 * (1 + dtOld1/dt1);

      alpha2 = alpha2/alpha0;
      alpha1 = alpha1/alpha0;
      alpha0 = -1/alpha0;

      double qscalar2(alpha2/sec.oldeTimeStep);
      bool Transpose = true;
      lastDQdxLambda.putScalar(0.0);
      dQdx.matvec( Transpose , lastLambda, lastDQdxLambda);
      lastDQdxLambda.scale(-1);
      RHSVec.update(+qscalar2, lastDQdxLambda);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Gear12::obtainJacobian
// Purpose       : Calculate Jacobian
// Special Notes : 
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------
void Gear12::obtainJacobian()
{

  if (DEBUG_TIME && isActive(Diag::TIME_JACOBIAN))
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::obtainJacobian" << std::endl;
  }

  // output: ds.JMatrixPtr

  // This function returns the following matrix:
  // $-(sec.alphas_/hn)dQdx(x)+dFdx$

  // Note:  ds.nextSolutionPtr is used to get dQdx,dFdx in Analysis::AnalysisManager::loadJacobian.

  Linear::Matrix & dQdx = *(ds.dQdxMatrixPtr);
  Linear::Matrix & dFdx = *(ds.dFdxMatrixPtr);
  Linear::Matrix & Jac = *(ds.JMatrixPtr);

  double qscalar(sec.alpha_[0]/sec.currentTimeStep);
  double fscalar(1.0);

  Jac.linearCombo( qscalar, dQdx, fscalar, dFdx );

  if (DEBUG_TIME && isActive(Diag::TIME_JACOBIAN))
  {
    Xyce::dout() << "\n dFdx:" <<std::endl;
    dFdx.print(Xyce::dout());
    Xyce::dout() << "\n Total Jacobian:" <<std::endl;
    Jac.print(Xyce::dout());
    //    for (int i=0;i<3;++i)
    //    {
    //      printf("[ %25.20g\t%25.20g\t%25.20g ]\n",Jac[i][0],Jac[i][1],Jac[i][2]);
    //    }

    Xyce::dout() << Xyce::section_divider << std::endl << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Gear12::interpolateSolution
// Purpose       : Interpolate solution approximation at prescribed time point.
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------
bool Gear12::interpolateSolution(double timepoint, 
    Linear::Vector * tmpSolVectorPtr, std::vector<Linear::Vector*> & historyVec)

{
  // this is a very course approximation to determine if we are too 
  // close to the actual time step to do an interpolation.
  // it could use more work. 
  double dtr = timepoint - sec.currentTime;  // the delta time requested.
  if( -dtr < 100 * Util::MachineDependentParams::MachinePrecision() )
  {
    *tmpSolVectorPtr = *(historyVec[0]);
    return false;
  }

  tmpSolVectorPtr->update(1.0, *(historyVec[0]), -1.0, *(historyVec[1]), 0.0);

  if( sec.usedOrder_ <= 2)
  {
    // do first order interpolation
    // X_interp = X + delta_t_requested * delta_X/delta_t[last step]
    dtr = dtr / sec.lastTimeStep;
    tmpSolVectorPtr->update(1.0, *(historyVec[0]), dtr);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::interpolateMPDESolution
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
bool Gear12::interpolateMPDESolution(std::vector<double>& timepoint, 
    Linear::Vector * tmpSolVectorPtr)
{
  Linear::BlockVector * blockTempSolVectorPtr = 
    dynamic_cast<Linear::BlockVector*>(tmpSolVectorPtr);
  if (blockTempSolVectorPtr == NULL)
  {
    Xyce::Report::DevelFatal0().in("Gear12::interpolateMPDESolution")
      << "Linear::Vector tmpSolVectorPtr is not of type Linear::BlockVector";
    return(false);
  }

  double tfuzz;   // fuzz factor to check for valid output time
  double tp;      // approximately t_{n-1}
  int numblocks = timepoint.size();
  int blockCount = blockTempSolVectorPtr->blockCount();
  if (numblocks > blockCount)
  {
    Xyce::Report::DevelFatal0().in("Gear12::interpolateMPDESolution")
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
        Xyce::Report::DevelFatal0().in("Gear12::interpolateMPDESolution") << "Linear::Vector ds.xHistory[j] is not of type Linear::BlockVector\n j = " << j;
        return(false);
      }
      xHistoryVectorPtr = &(blockXHistoryVectorPtr->block(i));
      solVectorPtr->update(c,*xHistoryVectorPtr);
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::printMPDEOutputSolution()
// Purpose       : Print transient output from MPDE simulation
// Special Notes : This routine uses interpolateMPDESolution.
// Scope         : public
// Creator       : Ting Mei, SNL, 1414
// Creation Date : 11/28/06
//-----------------------------------------------------------------------------
bool Gear12::printMPDEOutputSolution(
    Analysis::OutputMgrAdapter & outputManagerAdapter,
    const double time,
    Linear::Vector * solnVecPtr,
    const std::vector<double> & fastTimes )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::printWaMPDEOutputSolution()
// Purpose       : Print transient output from WaMPDE simulation
// Special Notes : This routine uses interpolateSolution.
// Scope         : public
// Creator       : Ting Mei, SNL, 1414
// Creation Date : 12/15/06
//-----------------------------------------------------------------------------
bool Gear12::printWaMPDEOutputSolution(
    Analysis::OutputMgrAdapter & outputManagerAdapter,
    const double time,
    Linear::Vector * solnVecPtr,
    const std::vector<double> & fastTimes,
    const int phiGID )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::printOutputSolution()
// Purpose       : Print output that is dumbed down in order.
// Special Notes : This routine picks smaller time steps to approximate first
//               : order integration from the perspective of the output.
// Scope         : public
// Creator       : Ting Mei, SNL, 1414
// Creation Date : 11/16/07 
//-----------------------------------------------------------------------------
bool Gear12::printOutputSolution(
  Analysis::OutputMgrAdapter &  outputManagerAdapter,
  const TIAParams &             tia_params, 
  const double                  time,
  Linear::Vector *                solnVecPtr,
  const bool                    doNotInterpolate,
  const std::vector<double> &   outputInterpolationTimes,
  bool                          skipPrintLineOutput)
{
  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  Gear12::printOutputSolution" << std::endl
                 << "usedOrder_ = " << sec.usedOrder_ << std::endl;
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
    for (unsigned int i=0;i<outputInterpolationTimes.size();++i)
    {
      double dt=0.0;
      if (i>0) 
      {
        dt = outputInterpolationTimes[i]-outputInterpolationTimes[i-1];
      }
      else
      {
        dt = 0.0;
      }

      interpolateSolution(outputInterpolationTimes[i], ds.tmpSolVectorPtr, ds.xHistory);    // interpolate solution vector
      interpolateSolution(outputInterpolationTimes[i], ds.tmpStaVectorPtr, ds.sHistory);    // interpolate state vector
      interpolateSolution(outputInterpolationTimes[i], ds.tmpStoVectorPtr, ds.stoHistory);  // interpolate store vector
      if (ds.leadCurrentSize)
      {
        interpolateSolution(outputInterpolationTimes[i], ds.tmpLeadCurrentVectorPtr, ds.leadCurrentHistory);  // interpolate store vector
        interpolateSolution(outputInterpolationTimes[i], ds.tmpLeadDeltaVPtr, ds.leadDeltaVHistory);  // interpolate store vector
        interpolateSolution(outputInterpolationTimes[i], ds.tmpLeadCurrentQDerivVectorPtr, ds.leadCurrentQDerivHistory);  // interpolate store vector
      }
      outputManagerAdapter.tranOutput(outputInterpolationTimes[i], dt, sec.finalTime,
          *ds.tmpSolVectorPtr, *ds.tmpStaVectorPtr, *ds.tmpStoVectorPtr, 
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
    outputManagerAdapter.tranOutput( time, timestep, sec.finalTime,
        *ds.currSolutionPtr, *ds.currStatePtr, *ds.currStorePtr, 
        *ds.currLeadCurrentPtr, *ds.currLeadDeltaVPtr, *ds.tmpLeadCurrentQDerivVectorPtr,
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
// Function      : Gear12::saveOutputSolution
// Purpose       : This is similar to printOutputSolution, but is in support of
//                 the .SAVE capability, rather than .PRINT.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
bool Gear12::saveOutputSolution( 
  Parallel::Machine                 comm,
  IO::InitialConditionsManager &    initial_conditions_manager,
  const NodeNameMap &               node_name_map,
  const TIAParams &                 tia_params,
  Linear::Vector *                  solnVecPtr,
  const double                      saveTime,
  const bool                        doNotInterpolate)
{
  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::saveOutputSolution" << std::endl;
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
    Xyce::dout() << Xyce::section_divider << std::endl;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::updateHistory
// Purpose       : Update history array after a successful step 
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/16/07
//-----------------------------------------------------------------------------
void Gear12::updateHistory()
{
  if (DEBUG_TIME && isActive(Diag::TIME_OUTPUT))
  {
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  Gear12::updateHistory" << std::endl
                 << "\n Before updates \n" << std::endl;
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

  if (sec.currentOrder_ == 2)
  {
    *(ds.xHistory[2]) = *(ds.xHistory[1]);
  }
  *(ds.qHistory[1]) = *(ds.qHistory[0]);
  *(ds.xHistory[1]) = *(ds.xHistory[0]);

  *(ds.sHistory[1]) = *(ds.sHistory[0]);
  *(ds.stoHistory[1]) = *(ds.stoHistory[0]);

  if (ds.leadCurrentSize)
  {
    *(ds.leadCurrentHistory[1]) = *(ds.leadCurrentHistory[0]);
    *(ds.leadCurrentQHistory[1]) = *(ds.leadCurrentQHistory[0]);
    *(ds.leadDeltaVHistory[1]) = *(ds.leadDeltaVHistory[0]);

    *(ds.leadCurrentHistory[0]) = *ds.nextLeadCurrentPtr;
    *(ds.leadCurrentQHistory[0]) = *ds.nextLeadCurrentQPtr;
    *(ds.leadDeltaVHistory[0]) = *ds.nextLeadDeltaVPtr;
  }

  *(ds.xHistory[0]) = *ds.nextSolutionPtr;
  *(ds.qHistory[0]) =  *ds.daeQVectorPtr;
  *(ds.sHistory[0]) =  *ds.nextStatePtr;    
  *(ds.stoHistory[0]) = *ds.nextStorePtr;

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
// Function      : Gear12::updateSensitivityHistory
// Purpose       : Update sensitivity history array after a successful step 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void Gear12::updateSensitivityHistory()
{
  if (ds.numParams)
  {
    if (sec.currentOrder_ == 2)
    {
      *(ds.dqdpHistory[1]) = *(ds.dqdpHistory[0]);
    }
    *(ds.dqdpHistory[0]) = *ds.nextDqdpPtrVector;
  }
}

//-----------------------------------------------------------------------------
// Function      : Gear12::updateSensitivityHistoryAdjoint
// Purpose       : Update sensitivity history array after a successful step 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void Gear12::updateSensitivityHistoryAdjoint()
{
  if (ds.numParams)
  {
    // to keep things simple always store 2 steps back.
    *(ds.dqdpHistory[1]) = *(ds.dqdpHistory[0]);
    *(ds.dfdpHistory[1]) = *(ds.dfdpHistory[0]);

    *(ds.dqdpHistory[0]) = *(ds.nextDqdpPtrVector);
    *(ds.dfdpHistory[0]) = *(ds.nextDfdpPtrVector);
  }
}

//-----------------------------------------------------------------------------
// Function      : Gear12::updateSensitivityHistoryAdjoint2
// Purpose       : Update sensitivity history array after a successful step 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void Gear12::updateSensitivityHistoryAdjoint2()
{
  if (ds.numParams)
  {
    Linear::MultiVector * tmpMV = 0;

    tmpMV = ds.nextDqdpPtrVector;
    ds.nextDqdpPtrVector = ds.dqdpHistory[0];
    ds.dqdpHistory[0] = ds.dqdpHistory[1];
    ds.dqdpHistory[1] = tmpMV;

    tmpMV = ds.nextDfdpPtrVector;
    ds.nextDfdpPtrVector = ds.dfdpHistory[0];
    ds.dfdpHistory[0] = ds.dfdpHistory[1];
    ds.dqdpHistory[1] = tmpMV;
  }
}

//-----------------------------------------------------------------------------
// Function      : Gear12::updateAdjointSensitivityHistory
// Purpose       : Update sensitivity history array after a successful step 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/23/2015
//-----------------------------------------------------------------------------
void Gear12::updateAdjointSensitivityHistory()
{
  if (ds.numParams)
  {
    *(ds.lastLambdaPtr) = *(ds.currLambdaPtr);
    *(ds.currLambdaPtr) = *(ds.nextLambdaPtr);

    *(ds.lastDQdxLambdaPtr) = *(ds.currDQdxLambdaPtr);
    *(ds.currDQdxLambdaPtr) = *(ds.nextDQdxLambdaPtr);
  }
}

//-----------------------------------------------------------------------------
// Function      : Gear12::restoreHistory
// Purpose       : Restore history array after a failed step
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07 
//-----------------------------------------------------------------------------
void Gear12::restoreHistory()
{
  for (int i=1;i<=sec.currentOrder_;++i)
  {
    sec.psi_[i-1] = sec.psi_[i];
  }
  if (DEBUG_TIME && isActive(Diag::TIME_HISTORY))
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::restoreHistory" << std::endl;
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
// Function      : Gear12::updateCoeffs
// Purpose       : Update method coefficients
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
void Gear12::updateCoeffs()
{
  // synchronize with Step Error Control
  //  sec.psi_[0] = sec.currentTimeStep;
  if (DEBUG_TIME && isActive(Diag::TIME_COEFFICIENTS))
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::updateCoeffs" << std::endl
      << "  currentTimeStep = " <<  sec.currentTimeStep << std::endl
      << "  numberOfSteps_ = " <<  sec.numberOfSteps_ << std::endl
      << "  currentOrder_ = " <<  sec.currentOrder_ << std::endl
      << "  nscsco_ = " <<  sec.nscsco_ << std::endl
      << "  psi_[0] = " <<  sec.psi_[0] << std::endl;
  }

  double temp1 = sec.currentTimeStep;

  if (sec.currentOrder_ == 2)
  {
    sec.psi_[2] = sec.psi_[1];
  }
  sec.psi_[1] = sec.psi_[0];
  sec.psi_[0] = temp1;

  //    sec.beta_[0] = 1.0;
  //    sec.alpha_[0] = 1.0; 
  sec.ck_ = 1.0;
  sec.alphas_ = -1.0;

  if (sec.currentOrder_ == 2)
  {
    // the coeffs of predictor
    sec.beta_[2] = temp1/sec.psi_[2] * (temp1 + sec.psi_[1])/(sec.psi_[1] + sec.psi_[2]);
    sec.beta_[1] = -temp1/sec.psi_[1] - sec.beta_[2] * (sec.psi_[1] + sec.psi_[2])/sec.psi_[1];
    sec.beta_[0] = 1.0 - sec.beta_[2] - sec.beta_[1];

    sec.alpha_[2] = -temp1/sec.psi_[1] * temp1/(2 * temp1 + sec.psi_[1]); 
    sec.alpha_[1] = 1 - sec.alpha_[2];
    sec.alpha_[0] = -sec.alpha_[1] - sec.alpha_[2] * (1 + sec.psi_[1]/temp1);

    sec.alpha_[2] = sec.alpha_[2]/sec.alpha_[0];
    sec.alpha_[1] = sec.alpha_[1]/sec.alpha_[0];
    sec.alpha_[0] = -1/sec.alpha_[0];

    sec.ck_ = sec.currentTimeStep/(temp1 + sec.psi_[1] + sec.psi_[2]);
  }
  else
  {
    sec.beta_[0] = 1.0 + temp1/sec.psi_[1];
    sec.beta_[1] = -temp1/sec.psi_[1];
    sec.alpha_[0] = 1.0;
    sec.alpha_[1] = -1.0;

    sec.ck_ = sec.currentTimeStep/(temp1 + sec.psi_[1]);
  }

  if (DEBUG_TIME && isActive(Diag::TIME_COEFFICIENTS))
  {
    Xyce::dout() << "  nscsco_ = " <<  sec.nscsco_ << std::endl
      << "  beta_[0] = " <<  sec.beta_[0] << std::endl
      << "  beta_[1] = " <<  sec.beta_[1] << std::endl
      << "  beta_[2] = " <<  sec.beta_[2] << std::endl
      << "  beta_[3] = " <<  sec.beta_[3] << std::endl
      << "  beta_[4] = " <<  sec.beta_[4] << std::endl
      << "  alpha_[0] = " <<  sec.alpha_[0] << std::endl
      << "  alpha_[1] = " <<  sec.alpha_[1] << std::endl
      << "  alpha_[2] = " <<  sec.alpha_[2] << std::endl
      << "  alpha_[3] = " <<  sec.alpha_[3] << std::endl
      << "  alpha_[4] = " <<  sec.alpha_[4] << std::endl
      << "  alphas_ = " <<  sec.alphas_ << std::endl
      << "  alpha0_ = " <<  sec.alpha0_ << std::endl
      << "  gamma_[0] = " <<  sec.gamma_[0] << std::endl
      << "  gamma_[1] = " <<  sec.gamma_[1] << std::endl
      << "  gamma_[2] = " <<  sec.gamma_[2] << std::endl
      << "  gamma_[3] = " <<  sec.gamma_[3] << std::endl
      << "  gamma_[4] = " <<  sec.gamma_[4] << std::endl
      << "  psi_[0] = " <<  sec.psi_[0] << std::endl
      << "  psi_[1] = " <<  sec.psi_[1] << std::endl
      << "  psi_[2] = " <<  sec.psi_[2] << std::endl
      << "  psi_[3] = " <<  sec.psi_[3] << std::endl
      << "  psi_[4] = " <<  sec.psi_[4] << std::endl
      << "  sigma_[0] = " <<  sec.sigma_[0] << std::endl
      << "  sigma_[1] = " <<  sec.sigma_[1] << std::endl
      << "  sigma_[2] = " <<  sec.sigma_[2] << std::endl
      << "  sigma_[3] = " <<  sec.sigma_[3] << std::endl
      << "  sigma_[4] = " <<  sec.sigma_[4] << std::endl
      << "  ck_ = " <<  sec.ck_ << std::endl
      << Xyce::section_divider << std::endl;
  }
}


//-----------------------------------------------------------------------------
// Function      : Gear12::updateAdjointCoeffs
// Purpose       : Update method coefficients for adjoint.  
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/24/2015
//-----------------------------------------------------------------------------
void Gear12::updateAdjointCoeffs()
{
  int it=ds.itAdjointIndex;

  if (ds.orderHistory[it] == 1)
  {
    sec.alpha_[0] = 1.00000000e+00;
    sec.alpha_[1] = -1.00000000e+00;
  }
  else
  {
    double dt1 = ds.dtHistory[it];
    double dtOld1 = ds.dtHistory[it-1];

    sec.alpha_[2] = -dt1/dtOld1 * dt1/(2 * dt1 + dtOld1); 
    sec.alpha_[1] = 1 - sec.alpha_[2];
    sec.alpha_[0] = -sec.alpha_[1] - sec.alpha_[2] * (1 + dtOld1/dt1);

    sec.alpha_[2] = sec.alpha_[2]/sec.alpha_[0];
    sec.alpha_[1] = sec.alpha_[1]/sec.alpha_[0];
    sec.alpha_[0] = -1/sec.alpha_[0];
  }
}

//-----------------------------------------------------------------------------
// Function      : Gear12::initialize
// Purpose       : Initialize method with initial solution & step-size
// Special Notes : 
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
void Gear12::initialize(const TIAParams &tia_params)
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
    if (rh>1.0)
      currentTimeStep = currentTimeStep/rh;

    sec.currentTimeStep = currentTimeStep;
  }


  sec.currentTimeStep = std::max(sec.currentTimeStep, sec.minTimeStep);
  sec.currentTimeStep = std::min(sec.currentTimeStep, sec.maxTimeStep);

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

  //  if (sec.currentTime == sec.initialTime)
  {
    // x history
    *(ds.xHistory[0]) = *(ds.currSolutionPtr);
    *(ds.xHistory[1]) = *(ds.currSolutionPtr);

    // q history
    *(ds.qHistory[0]) = *(ds.daeQVectorPtr);
    //  *(ds.qHistory[1]) = *(ds.daeFVectorPtr);
    //  (ds.qHistory[1])->scale(-sec.currentTimeStep);
    //  (ds.qHistory[1])->putScalar(0.0);
    *(ds.qHistory[1]) = *(ds.daeQVectorPtr);

    // state history
    *(ds.sHistory[0]) = *(ds.currStatePtr);
    (ds.sHistory[1])->putScalar(0.0);

    // store history
    *(ds.stoHistory[0]) = *(ds.currStorePtr);
    (ds.stoHistory[1])->putScalar(0.0);

    // lead current history
    if (ds.leadCurrentSize)
    {
      *(ds.leadCurrentHistory[0]) = *(ds.currLeadCurrentPtr);
      (ds.leadCurrentHistory[1])->putScalar(0.0);
      *(ds.leadCurrentQHistory[0]) = *(ds.currLeadCurrentQPtr);
      *(ds.leadCurrentQHistory[1]) = *(ds.currLeadCurrentQPtr);
      *(ds.leadDeltaVHistory[0]) = *(ds.currLeadDeltaVPtr);
      (ds.leadDeltaVHistory[1])->putScalar(0.0);
    }
  }

  // Coefficient initialization 
  sec.numberOfSteps_ = 0;    // number of total time integration steps taken
  sec.currentOrder_ = 1;
  //  sec.usedOrder_ = 1;
  //  if (sec.currentTime == sec.initialTime)
  //  {
  sec.usedOrder_ = 1;
  sec.psi_[0] = sec.currentTimeStep;
  sec.cj_ = 1/sec.psi_[0];
  //  }
  sec.nscsco_ = 0;
  if (DEBUG_TIME && isActive(Diag::TIME_HISTORY))
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::initialize" << std::endl
      << "\n xHistory: \n" << std::endl;
    (ds.xHistory[0])->print(Xyce::dout());
    Xyce::dout() << std::endl;
    (ds.xHistory[1])->print(Xyce::dout());
    Xyce::dout() << std::endl
      << "\n qHistory: \n" << std::endl;
    (ds.qHistory[0])->print(Xyce::dout());
    Xyce::dout() << std::endl;
    (ds.qHistory[1])->print(Xyce::dout());
    Xyce::dout() << std::endl
      << "\n sHistory: \n" << std::endl;
    (ds.sHistory[0])->print(Xyce::dout());
    Xyce::dout() << std::endl;
    (ds.sHistory[1])->print(Xyce::dout());
    Xyce::dout() << std::endl
      << "\n" << "currentTimeStep = " << currentTimeStep << "\n" << std::endl
      << "\n" << "time_to_stop = " << time_to_stop << "\n" << std::endl
      << Xyce::section_divider << std::endl;
  }

  initializeSensitivities();
}


//-----------------------------------------------------------------------------
// Function      : Gear12::initializeAdjoint
// Purpose       : Initialize method with initial solution & step-size
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/28/2016
//-----------------------------------------------------------------------------
void Gear12::initializeAdjoint (int index)
{
  ds.nextDQdxLambdaPtr->putScalar(0.0);
  ds.currDQdxLambdaPtr->putScalar(0.0);
  ds.lastDQdxLambdaPtr->putScalar(0.0);

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
// Function      : Gear12::initializeSensitivities
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void Gear12::initializeSensitivities()
{
  if (ds.numParams)
  {
    // dqdp history
    *(ds.dqdpHistory[0]) = *ds.currDqdpPtrVector;
    *(ds.dqdpHistory[1]) = *ds.currDqdpPtrVector;
  }
}

//-----------------------------------------------------------------------------
// Function      : Gear12::setTwoLevelTimeInfo
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/01/07
//-----------------------------------------------------------------------------
void Gear12::setTwoLevelTimeInfo()
{
  // x history
  *(ds.xHistory[0]) = *(ds.currSolutionPtr);
  *(ds.xHistory[1]) = *(ds.currSolutionPtr); 

  // q history
  *(ds.qHistory[0]) = *(ds.daeQVectorPtr);
  *(ds.qHistory[1]) = *(ds.daeQVectorPtr);

  // state history
  *(ds.sHistory[0]) = *(ds.currStatePtr);
  (ds.sHistory[1])->putScalar(0.0); 

  // Coefficient initialization 
  sec.numberOfSteps_ = 0;    // number of total time integration steps taken
  sec.usedOrder_ = 1;
  sec.psi_[0] = sec.currentTimeStep;
  sec.cj_ = 1/sec.psi_[0];
  sec.nscsco_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::checkReduceOrder()
// Purpose       : check whether to reduce order independent of local error test
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07 
//-----------------------------------------------------------------------------
void Gear12::checkReduceOrder()
{

}

//-----------------------------------------------------------------------------
// Function      : Gear12::rejectStep()
// Purpose       : code to restore history, choose new order/step-size
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07
//-----------------------------------------------------------------------------
void Gear12::rejectStep(const TIAParams & tia_params)
{
  // This routine puts its output in newTimeStep_ and sec.newOrder_

  // This routine changes the following variables:
  //    lastAttemptedTimeStep, sec.initialPhase_, sec.nef_, sec.psi_, newTimeStep_,
  //    sec.newOrder_, sec.currentOrder_, currentTimeStep_, ds.xHistory,
  //    ds.qHistory, nextTimePt, nextTime, currentTimeStepRatio,
  //    currentTimeStepSum, nextTimePt

  // This routine reades but does not change the following variables:
  //    stepAttemptStatus, sec.r_factor_, sec.r_safety_, sec.Est_, sec.r_fudge_, sec.r_min_, sec.r_max_,
  //    minTimeStep, maxTimeStep, currentTime, stopTime, lastTimeStep

  sec.TimeStepLimitedbyBP = false;

  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::rejectStep" << std::endl;
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
      // restore sec.psi_
      //    for (int i=1;i<=sec.currentOrder_;++i)
      //      sec.psi_[i-1] = sec.psi_[i] - sec.currentTimeStep;

      if (sec.nef_ >= sec.maxNumfail_)
      {
        Xyce::Report::DevelFatal0().in("Gear12::rejectStep")
          << "  Maximum number of failures at time "  << sec.currentTime;
      }

      if ((sec.newtonConvergenceStatus <= 0))
      {
        /// 11/11/05 erkeite:  If the Newton solver fails, don't 
        // rely on the error estimate - it may be full of Nan's.
        //        rr = sec.r_min_;

        newTimeStep_ = sec.currentTimeStep/8;
        //        sec.currentOrder_ = 1;
        sec.currentOrder_ = sec.minOrder_;
        //      if (sec.nef_ > 2) sec.newOrder_ = 1;//consistent with block below.
      }
      else
      {
        // 03/11/04 tscoffe:  Here is the block for choosing order & 
        // step-size when the local error test FAILS (but Newton 
        // succeeded). 
        if (sec.nef_ == 1) // first local error test failure
        {
          //	sec.estOverTol_ 
          sec.Est_ = sec.estOverTol_;

          rr = sec.tolAimFac_/(sec.estOverTol_ + 0.0001);
          rr = pow(rr, 1.0/(sec.currentOrder_+1.0));
          rr = std::max(sec.r_min_,std::min(sec.r_max_,rr));

          newTimeStep_ = rr * sec.currentTimeStep;

          //          sec.currentOrder_ = 1; 
          //          sec.currentOrder_ = sec.minOrder_;
        }
        else // if (sec.nef_ == 2) // second Dae failure
        {
          rr = sec.r_min_;
          newTimeStep_ = rr * sec.currentTimeStep;

          //          sec.currentOrder_ = 1;
          sec.currentOrder_ = sec.minOrder_;
        }
      }
      if (DEBUG_TIME && isActive(Diag::TIME_STEP) && isActive(Diag::TIME_COEFFICIENTS))
      {
        Xyce::dout() << "  currentTimeStep = " <<  sec.currentTimeStep << std::endl
          << "  numberOfSteps_ = " <<  sec.numberOfSteps_ << std::endl
          << "  currentOrder_ = " <<  sec.currentOrder_ << std::endl
          << "  nscsco_ = " <<  sec.nscsco_ << std::endl
          << "  alpha_[0] = " <<  sec.alpha_[0] << std::endl
          << "  alpha_[1] = " <<  sec.alpha_[1] << std::endl
          << "  alpha_[2] = " <<  sec.alpha_[2] << std::endl
          << "  alpha_[3] = " <<  sec.alpha_[3] << std::endl
          << "  alpha_[4] = " <<  sec.alpha_[4] << std::endl
          << "  psi_[0] = " <<  sec.psi_[0] << std::endl
          << "  psi_[1] = " <<  sec.psi_[1] << std::endl
          << "  psi_[2] = " <<  sec.psi_[2] << std::endl
          << "  psi_[3] = " <<  sec.psi_[3] << std::endl
          << "  psi_[4] = " <<  sec.psi_[4] << std::endl
          << "  sigma_[0] = " <<  sec.sigma_[0] << std::endl
          << "  sigma_[1] = " <<  sec.sigma_[1] << std::endl
          << "  sigma_[2] = " <<  sec.sigma_[2] << std::endl
          << "  sigma_[3] = " <<  sec.sigma_[3] << std::endl
          << "  sigma_[4] = " <<  sec.sigma_[4] << std::endl
          << "  rr = " <<  rr << std::endl
          << "  r_factor_ = " <<  sec.r_factor_ << std::endl
          << "  r_safety_ = " <<  sec.r_safety_ << std::endl
          << "  Est_ = " <<  sec.Est_ << std::endl
          << "  r_fudge_ = " <<  sec.r_fudge_ << std::endl
          << "  newOrder_ = " <<  sec.newOrder_ << std::endl
          << "  currentTimeStep = " <<  sec.currentTimeStep << std::endl
          << "  newTimeStep_ = " <<  newTimeStep_ << std::endl;
      }
    }
  }
  else if ((sec.stepAttemptStatus == false) & (!adjustStep))
  {
    std::string tmp = "  Gear12:rejectStep: Warning: Local error test failed with constant step-size.\n";
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
      Xyce::dout() << "  newTimeStep_ = " <<  newTimeStep_ << std::endl
        << "  nextTime = " <<  sec.nextTime << std::endl;
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
// Function      : Gear12::rejectStepForHabanero
// Purpose       : step rejection, but from an outside program (Habanero API)
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 08/11/09
//-----------------------------------------------------------------------------
void Gear12::rejectStepForHabanero()
{
  restoreHistory();
  sec.setTimeStep(sec.currentTimeStep);
}

//-----------------------------------------------------------------------------
// Function      : Gear12::completeStep()
// Purpose       : code to update history, choose new order/step-size
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 10/21/07 
//-----------------------------------------------------------------------------
void Gear12::completeStep(const TIAParams &tia_params)
{
  sec.TimeStepLimitedbyBP = false;

  sec.numberOfSteps_ ++;
  sec.nef_ = 0;
  sec.lastTime    = sec.currentTime;
  sec.currentTime = sec.nextTime;

  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "  Gear12::completeStep" << std::endl;
  }

  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = !tia_params.constantTimeStepFlag;

  sec.lastAttemptedTimeStep = sec.currentTimeStep;

  double newTimeStep_ = sec.currentTimeStep;
  double rr = 1.0; // step size ratio = new step / old step
  // 03/11/04 tscoffe:  Here is the block for choosing order & step-size when
  // the local error test PASSES (and Newton succeeded). 
  sec.oldeTimeStep = sec.lastTimeStep;
  sec.lastTimeStep = sec.currentTimeStep;
  sec.lastTimeStepRatio = sec.currentTimeStepRatio; // copied from calcTStep1
  sec.lastTimeStepSum   = sec.currentTimeStepSum; // copied from calcTStep1
  sec.usedOrder_ = sec.currentOrder_;
  sec.usedStep_ = sec.currentTimeStep;

  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << "  initialPhase_ = " <<  sec.initialPhase_ << std::endl
      << "  rr = " <<  rr << std::endl
      << "  currentTimeStep = " <<  sec.currentTimeStep << std::endl
      << "  currentTime = " <<  sec.currentTime << std::endl
      << "  nextTime = " <<  sec.nextTime << std::endl
      << "  newTimeStep_ = " <<  newTimeStep_ << std::endl
      << "  minTimeStep = " <<  sec.minTimeStep << std::endl
      << "  maxTimeStep = " <<  sec.maxTimeStep << std::endl;
  }

  // 03/22/04 tscoffe:  Note that updating the history has nothing to do with
  // the step-size and everything to do with the newton correction vectors.


  /*   if (sec.numberOfSteps_ >= 2)
       {
       sec.currentOrder_ = 2;
  //     (ds.relErrTolPtr)->putScalar(1e-2);	
  } 
  */ 
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
          //            sec.currentOrder_ = 1; 
          sec.currentOrder_ = sec.minOrder_;
        }
      }
    } 
    if (DEBUG_TIME && isActive(Diag::TIME_STEP))
    {
      Xyce::dout() << "  currentOrder_ = " <<  sec.currentOrder_ << std::endl;
      Xyce::dout() << "  r_safety = " <<  sec.r_safety_ << std::endl;
      Xyce::dout() << "  r_fudge_ = " <<  sec.r_fudge_ << std::endl;
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
  //  if (sec.currentTime < sec.stopTime)
  if ((sec.stopTime - sec.currentTime) >= sec.minTimeStep) 
  {
    // If the step needs to be adjusted:
    if (adjustStep)
    {
      //      newTimeStep_ = Xycemax(newTimeStep_, sec.minTimeStep);
      //      newTimeStep_ = Xycemin(newTimeStep_, sec.maxTimeStep);

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

//  sec.currentTimeStep = newTimeStep_;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::completeAdjointStep()
// Purpose       : code to update history, choose new order/step-size for 
//                 adjoint calculations.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void Gear12::completeAdjointStep(const TIAParams &tia_params)
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
// Function      : Gear12::updateStateDeriv
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 11/17/05
//-----------------------------------------------------------------------------
void Gear12::updateStateDeriv ()
{
  ds.nextStateDerivPtr->
    update(sec.alpha_[0],*ds.nextStatePtr, sec.alpha_[1],*(ds.sHistory[0]),0.0);

  if (sec.currentOrder_ == 2)
  {
    ds.nextStateDerivPtr->
      update(sec.alpha_[2], *(ds.sHistory[1]));
  }

  ds.nextStateDerivPtr->scale(1.0/sec.currentTimeStep);

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
// Function      : Gear12::updateLeadCurrentVec
// Purpose       : calculates lead currents in lead current vector with 
//                 the leadCurrQVec. 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling, SNL
// Creation Date : 03/22/2013
//-----------------------------------------------------------------------------
void Gear12::updateLeadCurrentVec ()
{
  if (ds.leadCurrentSize)
  {
    ds.nextLeadCurrentQDerivPtr->update(
        sec.alpha_[0], *ds.nextLeadCurrentQPtr, 
        sec.alpha_[1], *(ds.leadCurrentQHistory[0]),0.0);

    if (sec.currentOrder_ == 2)
    {
      ds.nextLeadCurrentQDerivPtr->
        update(sec.alpha_[2], *(ds.leadCurrentQHistory[1]));
    }

    ds.nextLeadCurrentQDerivPtr->scale(1.0/sec.currentTimeStep);

    ds.nextLeadCurrentPtr->update(1.0,*ds.nextLeadCurrentQDerivPtr);

    if (DEBUG_TIME && isActive(Diag::TIME_DUMP_SOLUTION_ARRAYS))
    {
      Xyce::dout() << "\n next leadcurrent q Ptr: \n" << std::endl;
      ds.nextLeadCurrentPtr->print(Xyce::dout());
      Xyce::dout() << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Gear12::getInitialQnorm
// Purpose       : Needed by 2-level solves.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/18/07
//-----------------------------------------------------------------------------
void Gear12::getInitialQnorm(TwoLevelError & tle) const
{
  tle.q1HistorySum = ds.partialSum_q1();
}

//-----------------------------------------------------------------------------
// Function      : Gear12::setupTwoLevelError
// Purpose       : Needed by 2-level solves.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
void Gear12::getTwoLevelError(TwoLevelError & tle) const
{
  tle.xErrorSum    = ds.partialErrorNormSum ();
  tle.qErrorSum    = ds.partialQErrorNormSum ();
  tle.xErrorSum_m1 = ds.partialSum_m1 (sec.currentOrder_);
  tle.xErrorSum_p1 = ds.partialSum_p1 (sec.currentOrder_, sec.maxOrder_);
  tle.innerSize    = ds.newtonCorrectionPtr->globalLength();

  if (DEBUG_TIME && isActive(Diag::TIME_ERROR))
    Xyce::dout() << tle;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::partialTimeDeriv
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/24/2016
//-----------------------------------------------------------------------------
double Gear12::partialTimeDeriv() const
{
  if (sec.currentTimeStep < 1e-30) 
  {
    Xyce::Report::UserWarning() 
      << "Excessively small current time step, incorrectly returning with large value";

    return leadingCoeff * 1.e+30;
  }
  
  return leadingCoeff / sec.currentTimeStep;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::getSolnVarData
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/24/2016
//-----------------------------------------------------------------------------
bool Gear12::getSolnVarData( const int & gid, std::vector<double> & varData )
{
  // Dump DataStore information.
  int num = ds.getNumSolnVarData();
  bool ret = ds.getSolnVarData( gid, varData );

  // Determine order for Gear12 dumping.
  int order = sec.currentOrder_;
  if (ret)
  {
    varData.resize( num + order + 3 );
    for (int i=0; i <= order; ++i )
    {
      varData[num++] = ds.xHistory[i]->getElementByGlobalIndex ( gid );
    }
    varData[num++] = ds.qHistory[0]->getElementByGlobalIndex ( gid );
    varData[num++] = ds.qHistory[1]->getElementByGlobalIndex ( gid );
  }
  return ret;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::setSolnVarData
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/24/2016
//-----------------------------------------------------------------------------
bool Gear12::setSolnVarData( const int & gid, const std::vector<double> & varData )
{
  // Restore data store information
  int num = ds.getNumSolnVarData();
  bool ret = ds.setSolnVarData( gid, varData );

  // Determine order for Gear12 restoring.
  int order = sec.currentOrder_;
  if (ret)
  {
    for (int i=0; i <= order; ++i )
    {
      ds.xHistory[i]->setElementByGlobalIndex( gid, varData[num++] );
    }
    ds.qHistory[0]->setElementByGlobalIndex( gid, varData[num++] );
    ds.qHistory[1]->setElementByGlobalIndex( gid, varData[num++] );
  }
  return ret;
}

//-----------------------------------------------------------------------------
// Function      : Gear12::getStateVarData
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/24/2016
//-----------------------------------------------------------------------------
bool Gear12::getStateVarData( const int & gid, std::vector<double> & varData )
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
// Function      : Gear12::setStateVarData
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/24/2016
//-----------------------------------------------------------------------------
bool Gear12::setStateVarData( const int & gid, const std::vector<double> & varData )
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
// Function      : Gear12::getStoreVarData
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/24/2016
//-----------------------------------------------------------------------------
bool Gear12::getStoreVarData( const int & gid, std::vector<double> & varData )
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
// Function      : Gear12::setStoreVarData
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/24/2016
//-----------------------------------------------------------------------------
bool Gear12::setStoreVarData( const int & gid, const std::vector<double> & varData )
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

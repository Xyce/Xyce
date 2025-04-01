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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 04/30/02
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ----------   Standard Includes   ----------
#include <string>
#include <iostream>
#include <cstdio>
#include <sstream>
#include <ios>

// ----------   Xyce Includes   ----------
#include <N_DEV_NumericalJacobian.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceMgr.h>

#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_ExtendedString.h>

// ---------- Static Initializations ----------


namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::NumericalJacobian
// Purpose       : block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/30/02
//-----------------------------------------------------------------------------
NumericalJacobian::NumericalJacobian(
  MatrixLoadData & mlData1,
  const SolverState &ss1,
  const ExternData  &ed1,
  const DeviceOptions & do1)
  : mlData(mlData1),
    cols(mlData1.cols),
    vals(mlData1.vals),
    Qvals(mlData1.Qvals),
    val_local(mlData1.val_local),
    Qval_local(mlData1.Qval_local),
    col_local(mlData1.col_local),
    row_local(mlData1.row_local),
    internalFlag(mlData1.internalFlag),
    devOptions(do1),
    solState (ss1),
    extData  (ed1),
    maxCols(10) // guess
{}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::NumericalJacobian
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/30/02
//-----------------------------------------------------------------------------
NumericalJacobian::NumericalJacobian(const NumericalJacobian & right)
  : mlData(right.mlData),
    cols(right.cols),
    vals(right.vals),
    Qvals(right.Qvals),
    val_local(right.val_local),
    Qval_local(right.Qval_local),
    col_local(right.col_local),
    row_local(right.row_local),
    internalFlag(right.internalFlag),
    devOptions(right.devOptions),
    solState (right.solState),
    extData  (right.extData),
    maxCols(right.maxCols)
{

}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::~NumericalJacobian
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/30/00
//-----------------------------------------------------------------------------
NumericalJacobian::~NumericalJacobian()
{

}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::testDAEMatrices
//
// Purpose       : Performs a numerical Jacobian test on the dFdx and dQdx
//                 matrices.
//
// Special Notes : This is the main numerical Jacoabian test.  
//
//                 The test compares the matrix stamps (numerical vs. analytic) 
//                 for a single device.  The numerical Jacobian values computed 
//                 in this class are not used by the solver.  Once the comparison
//                 is done the original analytical values are put back 
//                 in their original place.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/15/06
//-----------------------------------------------------------------------------
bool NumericalJacobian::testDAEMatrices(DeviceInstance & instance, const std::vector<const std::string *> & nameVec)
{

  // Set up various references to indexing arrays
  const std::vector<int> & devLIDs               = instance.getDevLIDs();
  const std::vector<int> & devStateLIDs          = instance.getStaLIDVec();
  const std::vector< std::vector<int> > & devJacLIDs = instance.getDevJacLIDs();
  const std::vector< std::vector<int> > & jacStamp   = instance.jacobianStamp();


  // Set up references to temporary data structures.
  std::vector< std::vector<double> > & numJacF = mlData.numJac;
  std::vector< std::vector<double> > & saveJacF = mlData.saveJac;
  std::vector< std::vector<double> > & devJacF = mlData.devJac;
  std::vector< std::vector<double> > & diffJacF = mlData.diffJac;
  std::vector< std::vector<double> > & relJacF = mlData.relJac;

  std::vector< std::vector<double> > & numJacQ = mlData.numJacQ;
  std::vector< std::vector<double> > & saveJacQ = mlData.saveJacQ;
  std::vector< std::vector<double> > & devJacQ = mlData.devJacQ;
  std::vector< std::vector<double> > & diffJacQ = mlData.diffJacQ;
  std::vector< std::vector<double> > & relJacQ = mlData.relJacQ;

  std::vector< std::vector<int> > & statusF  = mlData.status;
  std::vector< std::vector<int> > & statusQ  = mlData.statusQ;
  std::vector< std::vector<int> > & stencil = mlData.stencil;

  std::vector<double> & saveF = mlData.saveF;
  std::vector<double> & pertF = mlData.pertF;
  std::vector<double> & origF = mlData.origF;
  std::vector<double> & saveQ = mlData.saveQ;
  std::vector<double> & pertQ = mlData.pertQ;
  std::vector<double> & origQ = mlData.origQ;

  std::vector<double> & saveSoln = mlData.saveSoln;
  std::vector<double> & pertSoln = mlData.pertSoln;
  std::vector<double> & saveCurrSoln = mlData.saveCurrSoln;

  std::vector<double> & saveLastState = mlData.saveLastState;
  std::vector<double> & saveCurrState = mlData.saveCurrState;
  std::vector<double> & saveNextState = mlData.saveNextState;
  std::vector<double> & saveStateDerivs = mlData.saveStateDerivs;

  // set up references to linear algebra objects.
  Linear::Vector & Fvec          = (*extData.daeFVectorPtr);
  Linear::Vector & Qvec          = (*extData.daeQVectorPtr);

  Linear::Vector & currSol      = (*extData.currSolVectorPtr);
  Linear::Vector & nextSol      = (*extData.nextSolVectorPtr);

  Linear::Vector & lastSta      = (*extData.lastStaVectorPtr);
  Linear::Vector & currSta      = (*extData.currStaVectorPtr);
  Linear::Vector & nextSta      = (*extData.nextStaVectorPtr);
  Linear::Vector & nextStaDeriv = (*extData.nextStaDerivVectorPtr);

  Linear::Matrix & dQdxMat      = (*extData.dQdxMatrixPtr);
  Linear::Matrix & dFdxMat      = (*extData.dFdxMatrixPtr);

  int numRows, numCols;
  numCols = devLIDs.size();
  numRows = jacStamp.size();

  int testSize= (numRows>numCols)?numRows:numCols;
  if (testSize > saveF.size())
  {
    mlData.resizeTestJacSolData(testSize);
    mlData.resizeTestJacQData(testSize);
  }

  int numState = devStateLIDs.size();
  if (numState > saveCurrState.size())
  {
    mlData.resizeTestJacStateData(numState);
  }

  int i, j, jCol;

  if(devJacLIDs.empty())
  {
    if (DEBUG_DEVICE)
    {
      Report::UserWarning() << instance.getName() << " does not have jacLIDs available";
    }

    return true;
  }

  if (instance.getOrigFlag() && numRows > 0 && numRows == jacStamp.size())
  {

    // Zero out all the mlData structures
    saveF.assign(saveF.size(),0.0);
    pertF.assign(pertF.size(),0.0);
    origF.assign(origF.size(),0.0);

    saveQ.assign(saveQ.size(),0.0);
    pertQ.assign(pertQ.size(),0.0);
    origQ.assign(origQ.size(),0.0);

    saveSoln.assign(saveSoln.size(),0.0);
    pertSoln.assign(pertSoln.size(),0.0);
    saveCurrSoln.assign(saveCurrSoln.size(),0.0);

    saveLastState.assign(saveLastState.size(),0.0);
    saveCurrState.assign(saveCurrState.size(),0.0);
    saveNextState.assign(saveNextState.size(),0.0);
    saveStateDerivs.assign(saveStateDerivs.size(),0.0);
    for (i=0;i<numJacF.size();++i)
    {
      numJacF[i].assign(numJacF[i].size(),0.0);
      saveJacF[i].assign(saveJacF[i].size(),0.0);
      devJacF[i].assign(devJacF[i].size(),0.0);
      diffJacF[i].assign(diffJacF[i].size(),0.0);
      statusF[i].assign(statusF[i].size(),-1);
      stencil[i].assign(stencil[i].size(),0);

      numJacQ[i].assign(numJacQ[i].size(),0.0);
      saveJacQ[i].assign(saveJacQ[i].size(),0.0);
      devJacQ[i].assign(devJacQ[i].size(),0.0);
      diffJacQ[i].assign(diffJacQ[i].size(),0.0);
      statusQ[i].assign(statusQ[i].size(),-1);
    }

    // Save Soln, RHS, and State for this device
    bool origFlag = instance.getOrigFlag();
    //for (i=0 ; i<numRows ; ++i)
    int tmpSize= (numRows>numCols)?numRows:numCols;
    double sqrtEta=devOptions.testJac_SqrtEta;
    for (i=0 ; i<tmpSize; ++i)
    {
      saveF[i]      = Fvec[devLIDs[i]];
      saveQ[i]      = Qvec[devLIDs[i]];

      saveSoln[i]     = nextSol[devLIDs[i]];
      saveCurrSoln[i] = currSol[devLIDs[i]];
      pertSoln[i]     = sqrtEta * (1.0 + fabs(saveSoln[i]));
    }

    for (i=0 ; i<numState ; ++i)
    {
      saveLastState[i] = lastSta[devStateLIDs[i]];
      saveCurrState[i] = currSta[devStateLIDs[i]];
      saveNextState[i] = nextSta[devStateLIDs[i]];
      saveStateDerivs[i] = nextStaDeriv[devStateLIDs[i]];
    }

    // Save the original matrix for later:
    for (i=0 ; i<numRows ; ++i)
    {
      jCol = devJacLIDs[i].size();
      for (j=0 ; j<jCol ; ++j)
      {
        double valF = dFdxMat[devLIDs[i]][devJacLIDs[i][j]];
        saveJacF[i][j] = valF;
        double valQ = dQdxMat[devLIDs[i]][devJacLIDs[i][j]];
        saveJacQ[i][j] = valQ;
      }
    }

    // Zeroing needs to be done after all saved values are
    // recorded because there can be multiple references
    // to the same element
    for (i=0 ; i<numRows ; ++i)
    {
      jCol = devJacLIDs[i].size();
      for (j=0 ; j<jCol ; ++j)
      {
        dFdxMat[devLIDs[i]][devJacLIDs[i][j]] = 0;
        dQdxMat[devLIDs[i]][devJacLIDs[i][j]] = 0;
      }
    }

    // Now that the original load has been zeroed out, re-load the
    // analytic contributions, to get the contributions from *just* this
    // device.
    instance.loadDAEdQdx ();
    instance.loadDAEdFdx ();

    for (i=0 ; i<numRows ; ++i)
    {
      devJacF[i].assign(devJacF[i].size(),0.0);
      devJacQ[i].assign(devJacQ[i].size(),0.0);
      stencil[i].assign(stencil[i].size(),0);

      jCol = devJacLIDs[i].size();
      for (j=0 ; j<jCol ; ++j)
      {
        double valF = dFdxMat[devLIDs[i]][devJacLIDs[i][j]];
        devJacF[i][jacStamp[i][j]] = valF;
        double valQ = dQdxMat[devLIDs[i]][devJacLIDs[i][j]];
        devJacQ[i][jacStamp[i][j]] = valQ;
        stencil[i][jacStamp[i][j]] = 1;
      }
    }

    // zero out the RHS, and re-load, so that we have only the
    // elements from one device present.
    for (i=0 ; i<numRows ; ++i)
    {
      Fvec[devLIDs[i]] = 0.0;
      Qvec[devLIDs[i]] = 0.0;
    }
    // re-load for just this instance:
    loadLocalDAEVectors (instance);

    // Save RHS for just this instance:
    for (i=0 ; i<numRows ; ++i)
    {
      origF[i] = Fvec[devLIDs[i]];
      origQ[i] = Qvec[devLIDs[i]];
    }

    for (i=0 ; i<numRows ; ++i)
    {
      jCol = devJacLIDs[i].size();
      for (j=0 ; j<jCol ; ++j)
      {
        statusF[i][jacStamp[i][j]] = -1;
        statusQ[i][jacStamp[i][j]] = -1;
      }
    }

    // These are the tolerances for determining jacobian agreement
    //      relT: similar to reltol, fractional error that is allowed
    //      absT: similar to abstol, error that is allowed for small derivatives
    double relTol=devOptions.testJac_relTol;
    double absTol=devOptions.testJac_absTol;

    double ndFdx, adFdx, ddFdx, relError_dFdx;
    double ndQdx, adQdx, ddQdx, relError_dQdx;

    bool failedTest = false;
    for (i=0 ; i<numCols ; ++i)
    {
      // Don't bother perturbing gnd.
      if (*nameVec[devLIDs[i]] == "gnd") continue;

      for (j=0 ; j<numRows ; ++j)
      {
        nextSol[devLIDs[j]] = saveSoln[j];
        currSol[devLIDs[j]] = saveCurrSoln[j];
        Fvec[devLIDs[j]] = 0.0;
        Qvec[devLIDs[j]] = 0.0;
      }

      for (j=0 ; j<numState ; ++j)
      {
        lastSta[devStateLIDs[j]] = saveLastState[j];
        currSta[devStateLIDs[j]] = saveCurrState[j];
        nextSta[devStateLIDs[j]] = saveNextState[j];
        nextStaDeriv[devStateLIDs[j]] = saveStateDerivs[j];
      }

      // Perturb the solution.
      double dX = pertSoln[i];
      nextSol[devLIDs[i]] += dX;

      // Re-load the F,Q vectors:
      loadLocalDAEVectors (instance);

      for (j=0 ; j<numRows ; ++j)
      {
        pertF[j] = Fvec[devLIDs[j]];
        pertQ[j] = Qvec[devLIDs[j]];
      }

#ifdef Xyce_DEBUG_TESTJAC
      testDebugHead (instance, nameVec, i, dX);
#endif

      for (j=0 ; j<numRows ; ++j)
      {
        // Don't bother with the gnd row.
        if (*nameVec[devLIDs[j]] == "gnd") continue;

        // if this derivative is not loaded analytically, don't bother
        // with it.
        if (stencil[j][i]!=1) continue;

        double dF = (pertF[j]-origF[j]);
        double dQ = (pertQ[j]-origQ[j]);
        numJacF[j][i] = dF/dX;
        numJacQ[j][i] = dQ/dX;

        ndFdx = numJacF[j][i];
        adFdx = devJacF[j][i];
        ddFdx = fabs(adFdx-ndFdx);
        relError_dFdx = ddFdx/(relTol*fabs(ndFdx)+absTol);

        diffJacF[j][i] = ddFdx;
        relJacF[j][i] = relError_dFdx;

        ndQdx = numJacQ[j][i];
        adQdx = devJacQ[j][i];
        ddQdx = fabs(adQdx-ndQdx);
        relError_dQdx = ddQdx/(relTol*fabs(ndQdx)+absTol);

        diffJacQ[j][i] = ddQdx;
        relJacQ[j][i] = relError_dQdx;

#ifdef Xyce_DEBUG_TESTJAC
        if (std::isnan(numJacF[j][i]))
        {
          testDebugOut (instance, nameVec, i, j);
        }
#endif
        // if the device is a Inductor, and IC= has been specified,
        // then skip this term as it is a special case.
        ExtendedString varNameI(*nameVec[devLIDs[i]]); varNameI.toUpper();
        ExtendedString varNameJ(*nameVec[devLIDs[j]]); varNameJ.toUpper();
        if ( ((solState.dcopFlag) && varNameI[0]=='L' && varNameJ[0]=='L') )
        {
          // For the inductor branch current, the matrix element will be
          // there whether IC= was specified or not.
          if (adFdx == 1 && ndFdx == 0)
          {
            statusF[j][i] = 3;
            statusQ[j][i] = 3;
          }
        }
        // if the device is a capacitor, and it has a branch current,
        // that means that IC= has been used, thus it is a "special case"
        // that should be skipped. The branch current will only be there
        // for IC=.
        else if((!(solState.dcopFlag) && varNameI[0]=='C') && varNameJ[0]=='C')
        {
          statusF[j][i] = 3;
          statusQ[j][i] = 3;
        }

        if ( statusF[j][i] != 3 )
        {
          if (relError_dFdx > 1.0) // failure
          {
            statusF[j][i] = -2;
            failedTest = true;
          }
          else // success
          {
            statusF[j][i] = 1;
          }

          if (relError_dQdx > 1.0) // failure
          {
            statusQ[j][i] = -2;
            failedTest = true;
          }
          else // success
          {
            statusQ[j][i] = 1;
          }
        }
      }

#ifdef Xyce_DEBUG_TESTJAC
      testDebugTail (instance, nameVec);
#endif
    }

    // Output Jacobians Differences. If debug enabled, always output. otherwise only output for failures.
    if (DEBUG_DEVICE)
      printJacobian_(dout(), instance, nameVec, failedTest);
    else if (failedTest)
      printJacobian_(lout(), instance, nameVec, failedTest);

    // Restore jacobian, RHS for this device
    instance.setOrigFlag(origFlag);
    tmpSize= (numRows>numCols)?numRows:numCols;
    for (i=0 ; i<tmpSize; ++i)
    {
      Fvec[devLIDs[i]] = saveF[i];
      Qvec[devLIDs[i]] = saveQ[i];
      nextSol[devLIDs[i]] = saveSoln[i];
      currSol[devLIDs[i]] = saveCurrSoln[i];
    }
    for (i=0 ; i<numState ; ++i)
    {
      lastSta[devStateLIDs[i]] = saveLastState[i];
      currSta[devStateLIDs[i]] = saveCurrState[i];
      nextSta[devStateLIDs[i]] = saveNextState[i];
      nextStaDeriv[devStateLIDs[i]] = saveStateDerivs[i];
    }


    for (i=0 ; i<numRows ; ++i)
    {
      jCol = devJacLIDs[i].size();
      for (j=0 ; j<jCol ; ++j)
      {
        dFdxMat[devLIDs[i]][devJacLIDs[i][j]] = saveJacF[i][j];
        dQdxMat[devLIDs[i]][devJacLIDs[i][j]] = saveJacQ[i][j];
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::loadLocalDAEVectors
// Purpose       :
// Special Notes : Note this function ignores the B-vector.
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/15/06
//-----------------------------------------------------------------------------
void NumericalJacobian::loadLocalDAEVectors (DeviceInstance & instance)
{
  Linear::Vector & currSta      = (*extData.currStaVectorPtr);
  Linear::Vector & nextSta      = (*extData.nextStaVectorPtr);
  Linear::Vector & nextStaDeriv = (*extData.nextStaDerivVectorPtr);
  Linear::Vector & nextSol      = (*extData.nextSolVectorPtr);

  const std::vector<int> & devStateLIDs = instance.getStaLIDVec();
  int numState = devStateLIDs.size();

  instance.updateGlobalAndDependentParameters(false,false,false); // this line necessary for expressions
  instance.updatePrimaryState ();

  // Assume backward euler integration, so that the time integrator
  // accessors are not needed.
  //anaIntPtr_->updateDerivs(devLIDs, devStateLIDs);
  for (int j=0 ; j<numState ; ++j)
  {
    nextStaDeriv[devStateLIDs[j]] =
      solState.pdt_ * (nextSta[devStateLIDs[j]]-currSta[devStateLIDs[j]]);
  }

  instance.updateSecondaryState ();
  instance.loadDAEQVector ();
  instance.loadDAEFVector ();

  return;
}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::loadLocalDAEVectorsIncludingB
// Purpose       :
// Special Notes : Note this function does NOT ignore the B-vector.
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/4/2018
//-----------------------------------------------------------------------------
void NumericalJacobian::loadLocalDAEVectorsIncludingB (DeviceInstance & instance)
{
  Linear::Vector & currSta      = (*extData.currStaVectorPtr);
  Linear::Vector & nextSta      = (*extData.nextStaVectorPtr);
  Linear::Vector & nextStaDeriv = (*extData.nextStaDerivVectorPtr);
  Linear::Vector & nextSol      = (*extData.nextSolVectorPtr);

  const std::vector<int> & devStateLIDs = instance.getStaLIDVec();
  int numState = devStateLIDs.size();

  instance.updateGlobalAndDependentParameters(false,false,false); // this line necessary for expressions
  instance.updatePrimaryState ();

  // Assume backward euler integration, so that the time integrator
  // accessors are not needed.
  //anaIntPtr_->updateDerivs(devLIDs, devStateLIDs);
  for (int j=0 ; j<numState ; ++j)
  {
    nextStaDeriv[devStateLIDs[j]] =
      solState.pdt_ * (nextSta[devStateLIDs[j]]-currSta[devStateLIDs[j]]);
  }

  instance.updateSecondaryState ();
  instance.loadDAEQVector ();
  instance.loadDAEFVector ();
  instance.loadDAEBVector ();

  return;
}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::printJacobian_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/12/06
//-----------------------------------------------------------------------------
void NumericalJacobian::printJacobian_(
  std::ostream &                                os,
  const DeviceInstance &                        instance,
  const std::vector<const std::string *> &      nameVec,
  bool                                          failed)
{
  bool NAflag = false;
  const std::vector<int> & devLIDs               = instance.getDevLIDs();

  // These curly brackets are for scoping only.
  {
    const std::vector< std::vector<double> > & numJac = mlData.numJac;
    const std::vector< std::vector<double> > & anaJac = mlData.devJac;
    const std::vector< std::vector<double> > & diffJac = mlData.diffJac;
    const std::vector< std::vector<double> > & relJac = mlData.relJac;
    const std::vector< std::vector<int> > & stencil = mlData.stencil;
    const std::vector< std::vector<int> > & status = mlData.status;


    os << Xyce::section_divider << std::endl;
    os << "dFdx matrix for " << instance.getName();

    if (failed)
    {
      os << ":  JACOBIAN TEST FAILURE";
    }
    else
    {
      os << ":  JACOBIAN TEST SUCCESS";
    }
    os << " at time = " << solState.currTime_;
    os << " at Niter = " << solState.newtonIter;
    os << std::endl;
    os << "       Numerical     Analytic      absDiff       relative Error  Status   Names  (row, col)"<<std::endl;

    int i,j;
    int numCols = devLIDs.size();
    int numRows = (status.size()<numCols)?status.size():numCols;

    for (i = 0; i < numRows; ++i)
    {
      if (*nameVec[devLIDs[i]] == "gnd") continue;

      for (j = 0; j < numCols; ++j)
      {
        if (*nameVec[devLIDs[j]] == "gnd") continue;

        // if this variable has not been tested for any reason, skip.
        if (status[i][j]==-1) continue;

        // if this derivative is not loaded analytically, skip.
        if (stencil[i][j]!=1) continue;

        // Note: JT=jacobian test is there to make this easy to grep.
        std::stringstream prefix("FT:");
        std::stringstream tempChar("");
        tempChar.width(12);
        tempChar.precision(4);
        tempChar.flags(tempChar.flags() | std::ios_base::scientific);

        tempChar << prefix.str() << "  " << numJac[i][j] << "  " << anaJac[i][j] << "  " << diffJac[i][j] << "  " << relJac[i][j];
        if (status[i][j]==-2)
        {
          tempChar << "      fail";
        }
        else if (status[i][j]==3)
        {
          NAflag = true;
          tempChar << "      NA  ";
        }
        else
        {
          tempChar << "          ";
        }

        os << tempChar.str();
        os << "     ("<< *nameVec[devLIDs[i]]
            << ", " << *nameVec[devLIDs[j]]
            << ") "
            << " row,col=[ " << i << ", " << j << "]" << std::endl;
      }
    }
  }

  const std::vector< std::vector<double> > & numJac = mlData.numJacQ;
  const std::vector< std::vector<double> > & anaJac = mlData.devJacQ;
  const std::vector< std::vector<double> > & diffJac = mlData.diffJacQ;
  const std::vector< std::vector<double> > & relJac = mlData.relJacQ;
  const std::vector< std::vector<int> > & stencil = mlData.stencil;
  const std::vector< std::vector<int> > & status = mlData.statusQ;

  os << "dQdx matrix for " << instance.getName();

  if (failed)
  {
    os << ":  JACOBIAN TEST FAILURE";
  }
  else
  {
    os << ":  JACOBIAN TEST SUCCESS";
  }
  os << " at time = " << solState.currTime_;
  os << std::endl;
  os << "       Numerical     Analytic      absDiff       relative Error  Status   Names  (row, col)"<<std::endl;

  int i,j;
  int numCols = devLIDs.size();
  int numRows = (status.size()<numCols)?status.size():numCols;

  for (i = 0; i < numRows; ++i)
  {
    if (*nameVec[devLIDs[i]] == "gnd") continue;

    for (j = 0; j < numCols; ++j)
    {
      if (*nameVec[devLIDs[j]] == "gnd") continue;

      // if this variable has not been tested for any reason, skip.
      if (status[i][j]==-1) continue;

      // if this derivative is not loaded analytically, skip.
      if (stencil[i][j]!=1) continue;

      // Note: QT=jacobian test is there to make this easy to grep.
      std::stringstream tempChar("");
      tempChar.width(12);
      tempChar.precision(4);
      tempChar.flags(tempChar.flags() | std::ios_base::scientific);
      tempChar << "QT:" << "  " << numJac[i][j] << "  " << anaJac[i][j] << "  " << diffJac[i][j] << "  " << relJac[i][j];
      if (status[i][j]==-2)
      {
        tempChar << "      fail";
      }
      else if (status[i][j]==3)
      {
        NAflag = true;
        tempChar << "      NA  ";
      }
      else
      {
        tempChar << "          ";
      }

      os << tempChar.str();
      os << "     ("<< *nameVec[devLIDs[i]]
          << ", " << *nameVec[devLIDs[j]]
          << ") "
          << " row,col=[ " << i << ", " << j << "]" << std::endl;
    }
  }

  if(NAflag)
    os << " Note:  NA = untestable special case, such as IC=, etc." << std::endl;
  os << Xyce::section_divider << std::endl;

  if (failed)
  {
    if (!(devOptions.testJacWarn))
    {
      Report::UserError() << "Numerical Jacobian test failure" << std::endl
                          << "If you want this failure to be a warning, rather than an error, "
                          << "then add .options device testjacwarn=1 to the netlist.";
    }
    else
    {
      Report::UserWarning() << "Numerical Jacobian test failure";
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::testDebugHead
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/14/06
//-----------------------------------------------------------------------------
void NumericalJacobian::testDebugHead( DeviceInstance & instance, const std::vector<const std::string *> & nameVec, int i, double dX)
{
  const std::vector<int> & devLIDs = instance.getDevLIDs();

  Xyce::dout() << Xyce::section_divider<<std::endl;
  Xyce::dout() << "Perturbing (LID="<<devLIDs[i]<<") " << *nameVec[devLIDs[i]] << " by " << dX << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::testDebugOut
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/14/06
//-----------------------------------------------------------------------------
void NumericalJacobian::testDebugOut( DeviceInstance & instance, const std::vector<const std::string *> & nameVec, int i, int j)
{
  const std::vector<int> & devLIDs = instance.getDevLIDs();
  const std::vector<double> & pertF = mlData.pertF;
  const std::vector<double> & origF = mlData.origF;

  const std::vector< std::vector<double> > & numJac = mlData.numJac;
  const std::vector< std::vector<double> > & relJac = mlData.relJac;

  Xyce::dout().width(15); Xyce::dout().precision(7); Xyce::dout().setf(std::ios::scientific);
  Xyce::dout() << "dFdX: ";
  Xyce::dout() << " (" << devLIDs[j] << ", " << devLIDs[i] << ") ";
  Xyce::dout() << numJac[j][i];
  Xyce::dout() << " Forig = " << origF[j];
  Xyce::dout() << " Fperturb = " << pertF[j];
  double dF = -(pertF[j]-origF[j]);
  Xyce::dout() << " dF = " << dF;
  Xyce::dout() << " (" << *nameVec[devLIDs[j]] << ", " << *nameVec[devLIDs[i]] << ") ";
  Xyce::dout() << std::endl;
  Xyce::dout() << "  relative error = " << relJac[j][i] << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : NumericalJacobian::testDebugTail
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/14/06
//-----------------------------------------------------------------------------
void NumericalJacobian::testDebugTail( DeviceInstance & instance, const std::vector<const std::string *> & nameVec)
{
  Xyce::dout() << Xyce::section_divider<<std::endl;
}

} // namespace Device
} // namespace Xyce

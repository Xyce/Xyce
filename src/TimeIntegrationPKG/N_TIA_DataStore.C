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
// Purpose       : This file creates & initializes the data arrays needed for
//                 the time integration algorithms.
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

// ---------- Standard Includes ----------
#include <iostream>

#include <N_TIA_DataStore.h>

#include <N_LAS_Builder.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_FilteredMultiVector.h>
#include <N_LAS_Vector.h>
#include <N_TIA_TIAParams.h>
#include <N_UTL_DeleteList.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_MachDepParams.h>
#include <N_MPDE_Builder.h>

namespace Xyce {
namespace TimeIntg {

//-----------------------------------------------------------------------------
// Function      : DataStore::DataStore
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
DataStore::DataStore(
  int                           max_order,
  const Linear::Builder &       linear_system_builder)
  : builder_(linear_system_builder),
    limiterFlag(false),
    maxOrder(max_order),
    solutionSize(0),
    stateSize(0),
    storeSize(0),
    leadCurrentSize(0),
    tmpSolVectorPtr(0),
    tmpStaVectorPtr(0),
    tmpStaDerivPtr(0),
    tmpStoVectorPtr(0),
    tmpLeadCurrentVectorPtr(0),
    tmpLeadDeltaVPtr(0),
    tmpLeadCurrentQDerivVectorPtr(0),
    xn0Ptr (0),
    currSolutionPtr(0),
    lastSolutionPtr(0),
    nextSolutionPtr(0),
    savedNextSolutionPtr(0),
    currStatePtr(0),
    lastStatePtr(0),
    nextStatePtr(0),
    currStorePtr(0),
    nextStorePtr(0),
    currLeadCurrentPtr(0),
    nextLeadCurrentPtr(0),
    currLeadDeltaVPtr(0),
    nextLeadDeltaVPtr(0),
    currLeadCurrentQPtr(0),
    nextLeadCurrentQPtr(0),
    numParams(0),
    sensRHSPtrVector(0),
    sparseSensRHSMV(0),
    nextDfdpPtrVector(0),
    currDqdpPtrVector(0),
    nextDqdpPtrVector(0),
    nextDbdpPtrVector(0),
    currDXdpPtrVector(0),
    nextDXdpPtrVector(0),
    currDQdxDXdpPtrVector(0),
    lastDQdxDXdpPtrVector(0),
    nextDQdxDXdpPtrVector(0),
    currDFdxDXdpPtrVector(0),
    lastDFdxDXdpPtrVector(0),
    nextDFdxDXdpPtrVector(0),
    currStateDerivPtr(0),
    nextStateDerivPtr(0),
    currLeadCurrentQDerivPtr(0),
    nextLeadCurrentQDerivPtr(0),
    errWtVecPtr(0),
    JMatrixPtr(0),
    RHSVectorPtr(0),
    newtonCorrectionPtr(0),
    qNewtonCorrectionPtr(0),
    deviceErrorWeightMask_(0),
    qErrWtVecPtr(0),
    daeQVectorPtr(0),
    daeFVectorPtr(0),
    daeBVectorPtr(0),
    dQdxMatrixPtr(0),
    dFdxMatrixPtr(0),
    dQdxVecVectorPtr(0),
    dFdxVecVectorPtr(0),
    dFdxdVpVectorPtr(0),
    dQdxdVpVectorPtr(0),
    itAdjointIndex(0),
    adjointDcop(true),
    qn0Ptr(0),
    delta_x(0),
    tmpXn0APtr(0),
    tmpXn0BPtr(0),
    nextSolPtrSwitched_(false),
    absErrTol_(0.0),
    relErrTol_(0.0),
    solsMaxValue(0.0),
    maxSolutionPtr(0),
    relSolutionPtr(0),
    index(0),
    allocateSensitivityArraysComplete_(false),
    includeTransientAdjoint_(false)
{
  // temporary vectors:
  tmpSolVectorPtr = linear_system_builder.createVector();
  tmpStaVectorPtr = linear_system_builder.createStateVector();
  tmpStaDerivPtr = linear_system_builder.createStateVector();
  tmpStoVectorPtr = linear_system_builder.createStoreVector();
  tmpLeadCurrentVectorPtr = linear_system_builder.createLeadCurrentVector();
  tmpLeadDeltaVPtr = linear_system_builder.createLeadCurrentVector();
  tmpLeadCurrentQDerivVectorPtr = linear_system_builder.createLeadCurrentVector();
  
  xn0Ptr = linear_system_builder.createVector();

  // solution vectors:
  currSolutionPtr = linear_system_builder.createVector();
  lastSolutionPtr = linear_system_builder.createVector();
  nextSolutionPtr = linear_system_builder.createVector();
  solutionSize = nextSolutionPtr->localLength();  // get local length

  // state vectors:
  currStatePtr    = linear_system_builder.createStateVector();
  lastStatePtr    = linear_system_builder.createStateVector();
  nextStatePtr    = linear_system_builder.createStateVector();
  stateSize = nextStatePtr->globalLength();

  // store vectors:
  currStorePtr    = linear_system_builder.createStoreVector();
  nextStorePtr    = linear_system_builder.createStoreVector();
  storeSize = nextStorePtr->globalLength();

  // lead current and delta V vectors (for power calc).
  currLeadCurrentPtr= linear_system_builder.createLeadCurrentVector();
  nextLeadCurrentPtr= linear_system_builder.createLeadCurrentVector();
  currLeadDeltaVPtr= linear_system_builder.createLeadCurrentVector();
  nextLeadDeltaVPtr= linear_system_builder.createLeadCurrentVector();
  currLeadCurrentQPtr= linear_system_builder.createLeadCurrentVector();
  nextLeadCurrentQPtr= linear_system_builder.createLeadCurrentVector();
  currLeadCurrentQDerivPtr= linear_system_builder.createLeadCurrentVector();
  nextLeadCurrentQDerivPtr= linear_system_builder.createLeadCurrentVector();
  leadCurrentSize = nextLeadCurrentPtr->globalLength();

  // state derivative vectors:
  currStateDerivPtr = linear_system_builder.createStateVector();
  nextStateDerivPtr = linear_system_builder.createStateVector();

  // error vectors:
  errWtVecPtr      = linear_system_builder.createVector();

  errWtVecPtr->putScalar(1.0);

  relSolutionPtr = linear_system_builder.createVector();

  // device mask pointer
  deviceErrorWeightMask_ = linear_system_builder.createVector();
  deviceErrorWeightMask_->putScalar(1.0);

  // nonlinear solution vectors:
  newtonCorrectionPtr = linear_system_builder.createVector();
  qNewtonCorrectionPtr = linear_system_builder.createVector();

  // new-DAE stuff:
  // Error Vectors
  qErrWtVecPtr      = linear_system_builder.createVector();

  // DAE formulation vectors
  daeQVectorPtr      = linear_system_builder.createVector();
  daeFVectorPtr      = linear_system_builder.createVector();
  daeBVectorPtr      = linear_system_builder.createVector();

  // DAE formulation matrices
  dQdxMatrixPtr = linear_system_builder.createMatrix();
  dFdxMatrixPtr = linear_system_builder.createMatrix();

  // History arrays
  for (int i = 0; i < max_order + 1; ++i)
  {
    xHistory.push_back(linear_system_builder.createVector());
    qHistory.push_back(linear_system_builder.createVector());
    sHistory.push_back(linear_system_builder.createStateVector());
    stoHistory.push_back(linear_system_builder.createStoreVector());
    leadCurrentHistory.push_back(linear_system_builder.createLeadCurrentVector());
    leadCurrentQHistory.push_back(linear_system_builder.createLeadCurrentVector());
    leadDeltaVHistory.push_back(linear_system_builder.createLeadCurrentVector());
    leadCurrentQDerivHistory.push_back(linear_system_builder.createLeadCurrentVector());
  }

  // Predictors
  qn0Ptr = linear_system_builder.createVector();
}

//-----------------------------------------------------------------------------
// Function      : DataStore::DataStore
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
DataStore::~DataStore()
{
  delete tmpSolVectorPtr;
  delete tmpStaVectorPtr;
  delete tmpStaDerivPtr;
  delete tmpStoVectorPtr;
  delete tmpLeadCurrentVectorPtr;
  delete tmpLeadDeltaVPtr;
  delete tmpLeadCurrentQDerivVectorPtr;
  delete xn0Ptr;

  delete currSolutionPtr;
  delete lastSolutionPtr;

  delete nextSolutionPtr;

  delete currStatePtr;
  delete lastStatePtr;

  delete nextStatePtr;

  delete currStorePtr;
  delete nextStorePtr;

  // Lead current and power vectors
  delete currLeadCurrentPtr;
  delete nextLeadCurrentPtr;
  delete currLeadDeltaVPtr;
  delete nextLeadDeltaVPtr;

  // for lead current calculations.  F component is
  // held in the store vector, Q component is here
  delete currLeadCurrentQPtr;
  delete nextLeadCurrentQPtr;
  delete currLeadCurrentQDerivPtr;
  delete nextLeadCurrentQDerivPtr;

  delete currStateDerivPtr;

  delete nextStateDerivPtr;

  delete relSolutionPtr;
  delete maxSolutionPtr;

  delete errWtVecPtr;

  delete deviceErrorWeightMask_;

  delete newtonCorrectionPtr;
  delete qNewtonCorrectionPtr;

  //new-DAE:
  // Error Vectors
  delete qErrWtVecPtr;

  // DAE formulation vectors
  delete daeQVectorPtr;
  delete daeFVectorPtr;
  delete daeBVectorPtr;

  // DAE formulation matrices
  delete dQdxMatrixPtr;
  delete dFdxMatrixPtr;

  // HB temporary vectors
  delete dQdxVecVectorPtr;
  delete dFdxVecVectorPtr;

  Xyce::deleteList(xHistory.begin(), xHistory.end());
  Xyce::deleteList(qHistory.begin(), qHistory.end());
  Xyce::deleteList(sHistory.begin(), sHistory.end());
  Xyce::deleteList(stoHistory.begin(), stoHistory.end());
  Xyce::deleteList(leadCurrentHistory.begin(), leadCurrentHistory.end());
  Xyce::deleteList(leadCurrentQHistory.begin(), leadCurrentQHistory.end());
  Xyce::deleteList(leadDeltaVHistory.begin(), leadDeltaVHistory.end());
  Xyce::deleteList(leadCurrentQDerivHistory.begin(), leadCurrentQDerivHistory.end());

  delete qn0Ptr;

  // Step-size selection temporary vectors
  delete delta_x;

  // Temporary vector for WaMPDE interpolation
  delete tmpXn0APtr;
  delete tmpXn0BPtr;

  // Delete data in the fast time storage for HB and MPDE
  Xyce::deleteList(fastTimeSolutionVec.begin(), fastTimeSolutionVec.end());
  Xyce::deleteList(fastTimeStateVec.begin(), fastTimeStateVec.end());
  Xyce::deleteList(fastTimeQVec.begin(), fastTimeQVec.end());
  Xyce::deleteList(fastTimeStoreVec.begin(), fastTimeStoreVec.end());

  // delete the sensitivity related stuff, if necessary
  deleteSensitivityArrays();
}

//-----------------------------------------------------------------------------
// Function      : DataStore::deleteSensitivityArrays
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/6/2014
//-----------------------------------------------------------------------------
void DataStore::deleteSensitivityArrays()
{
  // adjoint stuff  (this should eventually have its own if-statement)

  Xyce::deleteList( solutionHistory.begin(), solutionHistory.end());
  Xyce::deleteList( stateHistory.begin(), stateHistory.end());
  Xyce::deleteList( storeHistory.begin(), storeHistory.end());

  if ( allocateSensitivityArraysComplete_ && includeTransientAdjoint_)
  {
    delete nextLambdaPtr;
    delete currLambdaPtr;
    delete lastLambdaPtr;

    delete nextDQdxLambdaPtr;
    delete currDQdxLambdaPtr;
    delete lastDQdxLambdaPtr;

    delete nextDFdxLambdaPtr;
    delete currDFdxLambdaPtr;
    delete lastDFdxLambdaPtr;
  
    delete tmpMatrixPtr;
  }

  if (numParams)
  {
    delete sensRHSPtrVector;
    delete sparseSensRHSMV;

    delete nextDfdpPtrVector;
    delete currDqdpPtrVector;
    delete nextDqdpPtrVector;
    delete nextDbdpPtrVector;

    delete currDXdpPtrVector;
    delete nextDXdpPtrVector;

    if (includeTransientDirect_)
    {
      delete currDQdxDXdpPtrVector;
      delete lastDQdxDXdpPtrVector;
      delete nextDQdxDXdpPtrVector;

      delete currDFdxDXdpPtrVector;
      delete lastDFdxDXdpPtrVector;
      delete nextDFdxDXdpPtrVector;
    }

    delete nextDqdpDerivPtrVector;

    for (unsigned int i=0; i<functionSensitivityHistory.size(); ++i)
    {
      delete functionSensitivityHistory[i];
    }
    functionSensitivityHistory.clear();

    for (unsigned int i=0; i<sparseFunctionSensitivityHistory.size(); ++i)
    {
      delete sparseFunctionSensitivityHistory[i];
    }
    sparseFunctionSensitivityHistory.clear();

    for (int i = 0; i < dqdpHistory.size(); ++i)
    {
      delete dqdpHistory[i];
      delete dfdpHistory[i];
      delete dbdpHistory[i];
    }

    // history
    dqdpHistory.clear();
    dfdpHistory.clear();
    dbdpHistory.clear();
  }
}

//-----------------------------------------------------------------------------
// Function      : DataStore::allocateHBVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/13/2019
//-----------------------------------------------------------------------------
void DataStore::allocateHBVectors()
{
  // Temporary vectors HB matrix vector products
  if (!dQdxVecVectorPtr && !dFdxVecVectorPtr)
  {
    dQdxVecVectorPtr = builder_.createVector();
    dFdxVecVectorPtr = builder_.createVector();
  }
}
//-----------------------------------------------------------------------------
// Function      : DataStore::allocateWaMPDEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/13/2019
//-----------------------------------------------------------------------------
void DataStore::allocateWaMPDEVectors()
{
  // Temporary vector for MPDE & WaMPDE interpolation
  if (!tmpXn0APtr && !tmpXn0BPtr)
  {
    tmpXn0APtr = builder_.createVector();
    tmpXn0BPtr = builder_.createVector();
  }
}

//-----------------------------------------------------------------------------
// Function      : DataStore::allocateSensitivityArrays
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/6/2014
//-----------------------------------------------------------------------------
void DataStore::allocateSensitivityArrays(int nP, 
        bool includeTransientDirect, bool includeTransientAdjoint)
{
  numParams = nP;

  includeTransientAdjoint_ = includeTransientAdjoint;
  includeTransientDirect_ = includeTransientDirect;

  if (!allocateSensitivityArraysComplete_) // only allocate if not already allocated
  {
    sensRHSPtrVector = builder_.createMultiVector( numParams );

    nextDfdpPtrVector = builder_.createMultiVector( numParams );

    currDqdpPtrVector = builder_.createMultiVector( numParams );
    nextDqdpPtrVector = builder_.createMultiVector( numParams );

    nextDbdpPtrVector = builder_.createMultiVector( numParams );

    currDXdpPtrVector = builder_.createMultiVector( numParams );
    nextDXdpPtrVector = builder_.createMultiVector( numParams );

    if (includeTransientDirect)
    {
      currDQdxDXdpPtrVector = builder_.createMultiVector( numParams );
      lastDQdxDXdpPtrVector = builder_.createMultiVector( numParams );
      nextDQdxDXdpPtrVector = builder_.createMultiVector( numParams );

      currDFdxDXdpPtrVector = builder_.createMultiVector( numParams );
      lastDFdxDXdpPtrVector = builder_.createMultiVector( numParams );
      nextDFdxDXdpPtrVector = builder_.createMultiVector( numParams );
    }

    nextDqdpDerivPtrVector = builder_.createMultiVector( numParams );

    // history
    dbdpHistory.resize(maxOrder+1);
    dfdpHistory.resize(maxOrder+1);
    dqdpHistory.resize(maxOrder+1);

    for (int i = 0; i < maxOrder + 1; ++i)
    {
      dqdpHistory[i] = builder_.createMultiVector( numParams );
      dbdpHistory[i] = builder_.createMultiVector( numParams );
      dfdpHistory[i] = builder_.createMultiVector( numParams );
    }

    // adjoint stuff  
    if (includeTransientAdjoint)
    {
      // adjoint sparse storage experiment
      sparseSensRHSMV = new Linear::FilteredMultiVector( numParams );
      masterIndexVector.resize(numParams);
      masterIndexVectorSize.resize(numParams);


      nextLambdaPtr = builder_.createVector();
      currLambdaPtr = builder_.createVector();
      lastLambdaPtr = builder_.createVector();

      nextDQdxLambdaPtr = builder_.createVector();
      currDQdxLambdaPtr = builder_.createVector();
      lastDQdxLambdaPtr = builder_.createVector();

      nextDFdxLambdaPtr = builder_.createVector();
      currDFdxLambdaPtr = builder_.createVector();
      lastDFdxLambdaPtr = builder_.createVector();
  
      tmpMatrixPtr = builder_.createMatrix();
    }

    allocateSensitivityArraysComplete_ = true;
  }
}

//-----------------------------------------------------------------------------
// Function      : DataStore::setConstantHistory
// Purpose       : This function is called  after the operating point
//                 calculation has been called.  Once the operating point
//                 solution has been obtained, the code should regard that
//                 solution as having been the existing constant solution since
//                 the dawn of time.
// Special Notes : The most recent solution, etc., are in the "next" vectors.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/23/00
//-----------------------------------------------------------------------------
void DataStore::setConstantHistory()
{
  if (DEBUG_TIME && isActive(Diag::TIME_HISTORY))
    Xyce::dout() << "\nDataStore::setConstantHistory" << std::endl;

  // Solutions:
  *(lastSolutionPtr) = *(nextSolutionPtr);
  *(currSolutionPtr) = *(nextSolutionPtr);

  if (stateSize)
  {
    // States:
    *(lastStatePtr) = *(nextStatePtr);
    *(currStatePtr) = *(nextStatePtr);

    // Derivative of States:
    *(currStateDerivPtr) = *(nextStateDerivPtr);
  }

  if (storeSize)
  { 
    // Stores:
    *(currStorePtr) = *(nextStorePtr);
  }

  if (leadCurrentSize)
  {
    // lead current and junction voltage vectors 
    *(currLeadCurrentPtr) = *(nextLeadCurrentPtr);
    *(currLeadCurrentQPtr) = *(nextLeadCurrentQPtr);
    *(currLeadDeltaVPtr) = *(nextLeadDeltaVPtr);
    *(currLeadCurrentQDerivPtr) = *(nextLeadCurrentQDerivPtr);
  }

  setConstantSensitivityHistory();
}

//-----------------------------------------------------------------------------
// Function      : DataStore::setConstantSensitivityHistory
// Purpose       :
// Special Notes : The most recent solution, etc., are in the "next" vectors.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void DataStore::setConstantSensitivityHistory()
{
  if (numParams)
  {
    *currDqdpPtrVector = *nextDqdpPtrVector;
    *currDXdpPtrVector = *nextDXdpPtrVector;

    if (includeTransientDirect_)
    {
      *lastDQdxDXdpPtrVector = *nextDQdxDXdpPtrVector;
      *currDQdxDXdpPtrVector = *nextDQdxDXdpPtrVector;

      *lastDFdxDXdpPtrVector = *nextDFdxDXdpPtrVector;
      *currDFdxDXdpPtrVector = *nextDFdxDXdpPtrVector;
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : DataStore::setConstantHistoryAdjoint 
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/9/2016
//-----------------------------------------------------------------------------
void DataStore::setConstantHistoryAdjoint ()
{
  if (DEBUG_TIME && isActive(Diag::TIME_HISTORY))
    Xyce::dout() << "\nDataStore::setConstantHistoryAdjoint" << std::endl;

  int finalPoint = timeHistory.size() - 1;

  // Solutions:
  *(nextSolutionPtr) = *(solutionHistory[finalPoint]);
  *(lastSolutionPtr) = *(nextSolutionPtr);
  *(currSolutionPtr) = *(nextSolutionPtr);

  if (stateSize)
  {
    // States:
    *(nextStatePtr) = *(stateHistory[finalPoint]);
    *(lastStatePtr) = *(nextStatePtr);
    *(currStatePtr) = *(nextStatePtr);

    // Derivative of States:  Fix later
    //*(currStateDerivPtr) = *(nextStateDerivPtr);
  }

  if (storeSize)
  { 
    // Stores:
    *(nextStorePtr) = *(storeHistory[finalPoint] );
    *(currStorePtr) = *(nextStorePtr);
  }

  // don't bother with lead currents (for now? maybe ever)
#if 0
  if (leadCurrentSize)
  {
    // lead current and junction voltage vectors 
    *(currLeadCurrentPtr) = *(nextLeadCurrentPtr);
    *(currLeadCurrentQPtr) = *(nextLeadCurrentQPtr);
    *(currLeadDeltaVPtr) = *(nextLeadDeltaVPtr);
    *(currLeadCurrentQDerivPtr) = *(nextLeadCurrentQDerivPtr);
  }
#endif

}


//-----------------------------------------------------------------------------
// Function      : DataStore::resetAll
//
// Purpose       : This function resets everything so that a transient loop
//                 can be started from the beginning.
//
// Special Notes : This function was needed for HB.
//
// Scope         : public
// Creator       : T. Mei, SNL
// Creation Date : 02/26/09
//-----------------------------------------------------------------------------
bool DataStore::resetAll(
  double absolute_error_tolerance,
  double relative_error_tolerance)
{
  absErrTol_ = absolute_error_tolerance; // tiaParams_.absErrorTol);
  relErrTol_ = relative_error_tolerance; // tiaParams_.relErrorTol);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::resetFastTimeData
//
// Purpose       : This function deletes the information from all the vectors
//                 that store fast time data for HB and MPDE
//
// Special Notes : This function was needed for HB.
//
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 06/05/13
//-----------------------------------------------------------------------------
bool DataStore::resetFastTimeData()
{
  // Clear the time step vectors
  timeSteps.clear();
  timeStepsBreakpointFlag.clear();

  // Delete any stored up any solution or state info
  Xyce::deleteList(fastTimeSolutionVec.begin(), fastTimeSolutionVec.end());
  Xyce::deleteList(fastTimeStateVec.begin(), fastTimeStateVec.end());
  Xyce::deleteList(fastTimeQVec.begin(), fastTimeQVec.end());
  Xyce::deleteList(fastTimeStoreVec.begin(), fastTimeStoreVec.end());

  fastTimeSolutionVec.clear();
  fastTimeStateVec.clear();
  fastTimeQVec.clear();
  fastTimeStoreVec.clear();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::updateSolDataArrays
// Purpose       : Update the necessary integration data arrays for
//                 preparation of the next integration step. This is done after
//                 a successful step has been taken.
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
void DataStore::updateSolDataArrays()
{
  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << "\nDataStore::updateSolDataArrays " << std::endl;
  }

  Linear::Vector * tmpPtr = 0;
  Linear::MultiVector * tmpMVPtr = 0;

  // if the next solution has been switched out (probably because of NOX),
  // then reset it
  if (nextSolPtrSwitched_)
    unsetNextSolVectorPtr();

  // Solutions:
  tmpPtr = lastSolutionPtr;
  lastSolutionPtr = currSolutionPtr;
  currSolutionPtr = nextSolutionPtr;
  nextSolutionPtr = tmpPtr;

  if (stateSize)
  {
    // States:
    tmpPtr = lastStatePtr;
    lastStatePtr = currStatePtr;
    currStatePtr = nextStatePtr;
    nextStatePtr = tmpPtr;

    // Derivative of States:
    tmpPtr = currStateDerivPtr;
    currStateDerivPtr = nextStateDerivPtr;
    nextStateDerivPtr = tmpPtr;
  }

  if (storeSize)
  {
    // Stores:
    tmpPtr = currStorePtr;
    currStorePtr = nextStorePtr;
    nextStorePtr = tmpPtr;
  }

  if (leadCurrentSize)
  {
    // lead current and junction voltage vectors
    tmpPtr = currLeadCurrentPtr;
    currLeadCurrentPtr = nextLeadCurrentPtr;
    nextLeadCurrentPtr = tmpPtr;

    tmpPtr = currLeadCurrentQPtr;
    currLeadCurrentQPtr = nextLeadCurrentQPtr;
    nextLeadCurrentQPtr = tmpPtr;
 
    tmpPtr = currLeadCurrentQDerivPtr;
    currLeadCurrentQDerivPtr = nextLeadCurrentQDerivPtr;
    nextLeadCurrentQDerivPtr = tmpPtr;
 
    tmpPtr = currLeadDeltaVPtr;
    currLeadDeltaVPtr = nextLeadDeltaVPtr;
    nextLeadDeltaVPtr = tmpPtr;
  }

  if (numParams)
  {
    tmpMVPtr = currDqdpPtrVector;
    currDqdpPtrVector = nextDqdpPtrVector;
    nextDqdpPtrVector = tmpMVPtr;

    tmpMVPtr = currDXdpPtrVector;
    currDXdpPtrVector = nextDXdpPtrVector;
    nextDXdpPtrVector = tmpMVPtr;

    if (includeTransientDirect_)
    {
      tmpMVPtr = lastDQdxDXdpPtrVector;
      lastDQdxDXdpPtrVector = currDQdxDXdpPtrVector;
      currDQdxDXdpPtrVector = nextDQdxDXdpPtrVector;
      nextDQdxDXdpPtrVector = tmpMVPtr; 

      tmpMVPtr = lastDFdxDXdpPtrVector;
      lastDFdxDXdpPtrVector = currDFdxDXdpPtrVector;
      currDFdxDXdpPtrVector = nextDFdxDXdpPtrVector;
      nextDFdxDXdpPtrVector = tmpMVPtr; 
    }
  }

  // copy contents of "curr" into "next".  This is to insure
  // that at a minimum, the initial guess for the Newton solve
  // will at least be the results of the previous Newton solve.
  *(nextSolutionPtr) = *(currSolutionPtr);
  if (stateSize)
    *(nextStatePtr)    = *(currStatePtr);
  if (storeSize)
    *(nextStorePtr)    = *(currStorePtr);
  if (leadCurrentSize)
  {
    *(nextLeadCurrentPtr) = *(currLeadCurrentPtr);
    *(nextLeadDeltaVPtr) = *(currLeadDeltaVPtr);
  }
}

//-----------------------------------------------------------------------------
// Function      : DataStore::updateStateDataArrays
//
//
// Purpose       : Same as updateSolDataArrays, but this function only
//                 advances the state vector, and leaves the
//                 solution alone.
//
// Special Notes : The main usage of this function is LOCA.
//                 After each continuation step, LOCA needs to call
//                 this function.  LOCA keeps track of solution vectors
//                 on its own, which is why updateSolDataArrays
//                 is inappropriate for LOCA.
//
//                 This is necessary for voltage limiting to be
//                 consistent with LOCA.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, 9233.
// Creation Date : 3/06/05
//-----------------------------------------------------------------------------
bool DataStore::updateStateDataArrays()
{
  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
    Xyce::dout() << "\nDataStore::updateStateDataArrays " << std::endl;

  Linear::Vector * tmpPtr = 0;

  if (stateSize)
  {
    // States:
    tmpPtr = lastStatePtr;
    lastStatePtr = currStatePtr;
    currStatePtr = nextStatePtr;
    nextStatePtr = tmpPtr;

    // Derivative of States:
    tmpPtr = currStateDerivPtr;
    currStateDerivPtr = nextStateDerivPtr;
    nextStateDerivPtr = tmpPtr;
  }

  if (storeSize)
  {
    // Stores:
    tmpPtr = currStorePtr;
    currStorePtr = nextStorePtr;
    nextStorePtr = tmpPtr;
  }

  if (leadCurrentSize)
  {
    // lead current and junction voltage vectors
    tmpPtr = currLeadCurrentPtr;
    currLeadCurrentPtr = nextLeadCurrentPtr;
    nextLeadCurrentPtr = tmpPtr;
 
    tmpPtr = currLeadCurrentQPtr;
    currLeadCurrentQPtr = nextLeadCurrentQPtr;
    nextLeadCurrentQPtr = tmpPtr;
  
    tmpPtr = currLeadDeltaVPtr;
    currLeadDeltaVPtr = nextLeadDeltaVPtr;
    nextLeadDeltaVPtr = tmpPtr;
  }
 
  // Now, make the "next" stuff the same as the "curr" stuff.
  // This is done because at the end of the tranop, but before
  // the transient phase starts, the function setConstantHistory
  // will be called.  When it is called, the most recent values
  // for state variables need to be in the "next" vectors.
  //
  // As long as LOCA solves are never used for transient, and only
  // for DC and tranop solves, this is OK.  This is, of course,
  // a bit of a kludge.
  if (stateSize)
    *(nextStatePtr) = *(currStatePtr);
  if (storeSize)
    *(nextStorePtr) = *(currStorePtr);
  if (leadCurrentSize)
  {
    *(nextLeadCurrentPtr) = *(currLeadCurrentPtr);
    *(nextLeadDeltaVPtr) = *(currLeadDeltaVPtr);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::updateSolDataArraysAdjoint
// Purpose       : This function unrolls the previously computed solution/state
//                 history.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 3/9/16
//-----------------------------------------------------------------------------
void DataStore::updateSolDataArraysAdjoint(int timeIndex)
{
  if (DEBUG_TIME && isActive(Diag::TIME_STEP))
  {
    Xyce::dout() << "\nDataStore::updateSolDataArraysAdjoint " << std::endl;
  }

  if(timeIndex==0)
  {
    adjointDcop = true;
  }
  else
  {
    adjointDcop = false;
  }

  int size = solutionHistory.size();

  // Solutions:
  if (timeIndex < size)
  {
    *(nextSolutionPtr) = *(solutionHistory[timeIndex]);
  }

  if (timeIndex < size-1)
  {
    *(currSolutionPtr) = *(solutionHistory[timeIndex+1]);
  }
  else
  {
    *(currSolutionPtr) = *(nextSolutionPtr);
  }

  if (timeIndex < size-2)
  {
    *(lastSolutionPtr) = *(solutionHistory[timeIndex+2]);
  }
  else
  {
    *(lastSolutionPtr) = *(currSolutionPtr);
  }

  if (stateSize)
  {
    // States:
    if (timeIndex < size)
    {
      *(nextStatePtr) = *(stateHistory[timeIndex]);
    }

    if (timeIndex < size-1)
    {
      *(currStatePtr) = *(stateHistory[timeIndex+1]);
    }
    else
    {
      *(currStatePtr) = *(nextStatePtr);
    }

    if (timeIndex < size-2)
    {
      *(lastStatePtr) = *(stateHistory[timeIndex+2]);
    }
    else
    {
      *(lastStatePtr) = *(currStatePtr);
    }

#if 0
    // Derivative of States:  Skip for now
    tmpPtr = currStateDerivPtr;
    currStateDerivPtr = nextStateDerivPtr;
    nextStateDerivPtr = tmpPtr;
#endif
  }

  if (storeSize)
  {
    // Stores:
    if (timeIndex < size)
    {
      *(nextStorePtr) = *(storeHistory[timeIndex]);
    }

    if (timeIndex < size-1)
    {
      *(currStorePtr) = *(storeHistory[timeIndex+1]);
    }
    else
    {
      *(currStorePtr) = *(nextStorePtr);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DataStore::partialErrorNormSum
// Purpose       : Needed by 2-level solves.
//
// Special Notes : A weighted RMS norm is this:
//
//                  norm = sqrt ( 1/n * sum (x/w)^2 )
//
//                  What this function returns is: sum (x/w)^2
//
//                  It will later be summed with other partial sums
//                  to get the complete WRMS value.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
double DataStore::partialErrorNormSum ()
{
  double errorNorm = 0.0;
  newtonCorrectionPtr->wRMSNorm(*errWtVecPtr, &errorNorm);
  double sum = errorNorm*errorNorm;
  double length = newtonCorrectionPtr->globalLength();

  sum *= length;

  return sum;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::equateTmpVectors
// Purpose       : This function equates the 6 temporary vectors  with
//                 their "next" vector equivalents.  This function
//                 is neccessary for the nonlinear solver damping loop.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/24/02
//-----------------------------------------------------------------------------
bool DataStore::equateTmpVectors()
{
  // next solution vector:
  *tmpSolVectorPtr = *nextSolutionPtr;

  if (stateSize)
  {
    // next state vector:
    *tmpStaVectorPtr = *nextStatePtr;

    // next state derivative vector:
    *tmpStaDerivPtr  = *nextStateDerivPtr;
  }

  // next store vector:
  if (storeSize)
    *tmpStoVectorPtr = *nextStorePtr;

  // next lead currrent and junction voltage vectors
  if (leadCurrentSize)
  {
    *tmpLeadCurrentVectorPtr = *nextLeadCurrentPtr;
    *tmpLeadDeltaVPtr = *nextLeadDeltaVPtr;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::usePreviousSolAsPredictor
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/28/01
//-----------------------------------------------------------------------------
bool DataStore::usePreviousSolAsPredictor ()
{
  bool bsuccess = true;

  *nextSolutionPtr = *currSolutionPtr;
  if (stateSize)
    *nextStatePtr    = *currStatePtr;
  if (storeSize)
    *nextStorePtr    = *currStorePtr;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::setNextSolVectorPtr
// Purpose       :
// Special Notes : Only needed for NOX, and other solvers that prefer to
//                 own the solution vector.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/22/03
//-----------------------------------------------------------------------------
bool DataStore::setNextSolVectorPtr (Linear::Vector * solVecPtr)
{
  // only save the old pointer if it hasn't been switched yet.
  if (!nextSolPtrSwitched_)
  {
    savedNextSolutionPtr = nextSolutionPtr;
    nextSolPtrSwitched_ = true;
  }
  nextSolutionPtr = solVecPtr;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::setNextSolVectorPtr
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Timur Takhtaganov, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool DataStore::setNextSolVectorPtr (Linear::Vector & solVecPtr)
{
  // only save the old pointer if it hasn't been switched yet.
  if (!nextSolPtrSwitched_)
  {
    savedNextSolutionPtr = nextSolutionPtr;
    nextSolPtrSwitched_ = true;
  }
  *(nextSolutionPtr) = solVecPtr;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::unsetNextSolVectorPtr
// Purpose       :
// Special Notes : This also copies over the solution, in addition to the
//                 pointer.
//
//                 This is only called when it is time to rotate the
//                 pointers for the next time step.  If we have been
//                 running with NOX, or with some other solver that prefers
//                 to own the solution vector, this is neccessary, and the
//                 switch flag should be set.  Otherwise, not.
//
//                 Basically, if we've been running with NOX, then the next
//                 solution vector ptr has probably been switched out at least
//                 once.  We need to maintain the history, so we make a
//                 copy of this switched solution vector, and the
//                 restore the old pointer.
//
//                 This is kludgy, but will have to do for now.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/22/03
//-----------------------------------------------------------------------------
bool DataStore::unsetNextSolVectorPtr ()
{
  if (nextSolPtrSwitched_)
  {
    *savedNextSolutionPtr = *nextSolutionPtr;
    nextSolutionPtr = savedNextSolutionPtr;
    nextSolPtrSwitched_ = false;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::setZeroHistory
// Purpose       : Sets everything to zero.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/24/07
//-----------------------------------------------------------------------------
void DataStore::setZeroHistory()
{
  if (DEBUG_TIME && isActive(Diag::TIME_HISTORY))
    Xyce::dout() << "\nDataStore::setZeroHistory" << std::endl;

  // Solutions:
  nextSolutionPtr->putScalar(0.0);

  if (stateSize)
  {
    // States:
    nextStatePtr->putScalar(0.0);
    nextStateDerivPtr->putScalar(0.0);
  }

  if (storeSize)
  {
    nextStorePtr->putScalar(0.0);
  }
  
  if (leadCurrentSize)
  {
    nextLeadCurrentPtr->putScalar(0.0);
    nextLeadCurrentQPtr->putScalar(0.0);
    nextLeadDeltaVPtr->putScalar(0.0);
  }

  if (numParams)
  {
    nextDfdpPtrVector->putScalar(0.0);
    nextDqdpPtrVector->putScalar(0.0);
    nextDbdpPtrVector->putScalar(0.0);
    nextDXdpPtrVector->putScalar(0.0);

    if (includeTransientDirect_)
    {
      nextDQdxDXdpPtrVector->putScalar(0.0);
      nextDFdxDXdpPtrVector->putScalar(0.0);
    }

    nextDqdpDerivPtrVector->putScalar(0.0);
  }

  qErrWtVecPtr->putScalar(0.0);

  // DAE formulation vectors
  daeQVectorPtr->putScalar(0.0);
  daeFVectorPtr->putScalar(0.0);
  daeBVectorPtr->putScalar(0.0);

  // Predictors
  xn0Ptr->putScalar(0.0);
  qn0Ptr->putScalar(0.0);

  // Nonlinear solution vector:
  qNewtonCorrectionPtr->putScalar(0.0);

  // This just sets the "oldDAE" history vectors to zero.
  setConstantHistory ();

  // new-DAE history:
  for (int i = 0; i < maxOrder + 1; ++i)
  {
    xHistory[i]->putScalar(0.0);
    qHistory[i]->putScalar(0.0);
    sHistory[i]->putScalar(0.0);
    stoHistory[i]->putScalar(0.0);
    leadCurrentHistory[i]->putScalar(0.0);
    leadCurrentQHistory[i]->putScalar(0.0);
    leadDeltaVHistory[i]->putScalar(0.0);

    if (numParams)
    {
      dqdpHistory[i]->putScalar(0.0);
      dbdpHistory[i]->putScalar(0.0);
      dfdpHistory[i]->putScalar(0.0);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DataStore::setErrorWtVector
// Purpose       : Set the Error Weight Vector (defined in terms of the
//                 solution approximation and error tolerances).
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL,Parallel Computational Sciences
// Creation Date : 5/25/05
//-----------------------------------------------------------------------------
void DataStore::setErrorWtVector(
  const TIAParams &             tia_params,
  const std::vector<char> &     variable_type)
{
//tia_params.maskIVars = true;
  // presort the variables into types
  if (indexVVars.empty() && indexMaskedVars.empty() && indexIVars.empty() )
  {
    // now fill the index arrays
    for (int k = 0; k < solutionSize; ++k)
    {
      if ((*deviceErrorWeightMask_)[k] == 0.0)
        indexMaskedVars.push_back(k);

      else if (tia_params.maskIVars && (variable_type[k] == 'I' ))
 //      else if (variable_type[k] == 'I')
        indexIVars.push_back(k);
      else
        indexVVars.push_back(k);
    }
  }

  if (DEBUG_TIME && isActive(Diag::TIME_ERROR))
  {
    Xyce::dout() << Xyce::section_divider << std::endl
                 << "DataStore::setErrorWtVector, relErrTol = " << relErrTol_ 
                 << ", absErrTol = " << absErrTol_ << std::endl << std::endl
                 << "   errorWtVector    currSolution " << std::endl
                 << "   --------------  --------------"  << std::endl;
  }


  switch(tia_params.newLte )
  {
    case 0:    // point local
    {
      relSolutionPtr->absValue(*currSolutionPtr); 
      if (DEBUG_TIME && isActive(Diag::TIME_ERROR)) 
      {
        double currMaxValue = 0.0;
        currSolutionPtr->infNorm(&currMaxValue);
        std::vector<int> index(1, -1);
        currSolutionPtr->infNormIndex( &index[0] );
        Xyce::dout() << "currMaxValueoldLte = " << currMaxValue << ", currMaxValueIndex = " << index[0]  << std::endl;
        
        std::cout << " old LTE rel reference:"  << std::endl;
        relSolutionPtr->printPetraObject(Xyce::dout());

      }
    }
      break;


    case 1:   //point global
    {  
      double currMaxValue = 0.0;
      currSolutionPtr->infNorm(&currMaxValue);

      relSolutionPtr->putScalar(currMaxValue);
      if (DEBUG_TIME && isActive(Diag::TIME_ERROR))
      {
        std::vector<int> index(1, -1);
        currSolutionPtr->infNormIndex( &index[0] );
        Xyce::dout() << "currMaxValue = " << currMaxValue << ", currMaxValueIndex = " << index[0] << std::endl;

        std::cout << " lte =1 rel reference:"  << std::endl;
        relSolutionPtr->printPetraObject(Xyce::dout());

        std::cout << " lte = 1 currSolution :"  << std::endl;
        currSolutionPtr->printPetraObject(Xyce::dout());
      }
     } 
      break;
    

    case 2:
    {
      double currMaxValue = 0.0;
      currSolutionPtr->infNorm(&currMaxValue);

      std::vector<int> currIndex(1, -1);
      currSolutionPtr->infNormIndex( &currIndex[0] );

//      int index;
      if (currMaxValue > solsMaxValue)
      {
        solsMaxValue = currMaxValue;
        index = currIndex[0]; 
      }

      relSolutionPtr->putScalar( solsMaxValue);
//      errWtVecPtr->putScalar( solsMaxValue); 
 
      if (DEBUG_TIME && isActive(Diag::TIME_ERROR))
      {
        Xyce::dout() << "signal global MaxValue = " << solsMaxValue << ", signal global MaxValueIndex = " << index << std::endl;
        Xyce::dout() << " currMaxValue = " << currMaxValue << ", currMaxValueIndex = " <<  currIndex[0] << std::endl;
      }
    }
      break;

    case 3:  //sig local 
    {
      // Create maxSolutionPtr only for this LTE calculation
      if ( !maxSolutionPtr )
      {
        maxSolutionPtr = builder_.createVector();
        maxSolutionPtr->putScalar( 0.0 );
      } 

      for ( int i=0; i< solutionSize; i++)
      {
        if (fabs( (*currSolutionPtr)[i]) >  (*maxSolutionPtr)[i] )
          (*maxSolutionPtr)[i] = fabs( (*currSolutionPtr)[i]);

      }

//      errWtVecPtr =  maxSolutionPtr;

      for ( int i=0; i< solutionSize; i++)
      {
//        if (fabs( (*nextSolutionPtr)[i]) >  (*maxSolutionPtr)[i] )
//          (*relSolutionPtr)[i] = fabs( (*nextSolutionPtr)[i]);
//        else
          (*relSolutionPtr)[i] = (*maxSolutionPtr)[i];
      }

      if (DEBUG_TIME && isActive(Diag::TIME_ERROR))
      {
        std::cout << " LTE 3 rel reference:"  << std::endl; std::cout << " LTE 3 rel reference:"  << std::endl;
        relSolutionPtr->printPetraObject(Xyce::dout()); 

        std::cout << " LTE = 3 next solutions  :"  << std::endl;
        nextSolutionPtr->printPetraObject(Xyce::dout()); 
      }

    }
    break; 

    default:
       std::cout << "Unsupported new LTE options" << std::endl;
  }

  qErrWtVecPtr->absValue(*daeQVectorPtr);


  // Voltage variables
  for (std::vector<int>::const_iterator it = indexVVars.begin(), end = indexVVars.end(); it != end; ++it)
  {
    int i = *it;
    (*errWtVecPtr)[i] = relErrTol_ * (*relSolutionPtr)[i] + absErrTol_;
    (*qErrWtVecPtr)[i] = relErrTol_ * (*qErrWtVecPtr)[i] + absErrTol_;
//    (*errWtVecPtr)[i] = relErrTol_ * (*errWtVecPtr)[i] + absErrTol_;
  }

  // Current variables
  // if !fastTests, then I vars are treated with V vars above.

  for (std::vector<int>::const_iterator it = indexIVars.begin(), end = indexIVars.end(); it != end; ++it)
  {
    int i = *it;
    (*errWtVecPtr)[i] = (*qErrWtVecPtr)[i] = Util::MachineDependentParams::MachineBig();
  }

  // Masked variables
  for (std::vector<int>::const_iterator it = indexMaskedVars.begin(), end = indexMaskedVars.end(); it != end; ++it)
  {
    int i = *it;
    (*errWtVecPtr)[i] = (*qErrWtVecPtr)[i] = Util::MachineDependentParams::MachineBig();
  }


  if (DEBUG_TIME && isActive(Diag::TIME_ERROR))
  {
    for (int k = 0; k < solutionSize; ++k)
      Xyce::dout() << (*errWtVecPtr)[k] << " "
                   << (*currSolutionPtr)[k] << std::endl; 

    Xyce::dout() << ""  << std::endl
                 << Xyce::section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : DataStore::WRMS_errorNorm
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/25/05
//-----------------------------------------------------------------------------
double DataStore::WRMS_errorNorm()
{
  double errorNorm = 0.0, qErrorNorm = 0.0;
  newtonCorrectionPtr->wRMSNorm(*errWtVecPtr, &errorNorm);
  qNewtonCorrectionPtr->wRMSNorm(*qErrWtVecPtr, &qErrorNorm);

  if (DEBUG_TIME && isActive(Diag::TIME_ERROR))
  {
    Xyce::dout() << "DataStore::errorNorm = " << errorNorm << std::endl;
    Xyce::dout() << "DataStore::qErrorNorm = " << qErrorNorm << std::endl;
  }

  // This is for handling the total errorNorm for 2-level solves.
  //
  // Note:  This version of the function assumes that the error norm
  // and q-error norm are the same size.  (they have to be...)
  if ( !(innerErrorInfoVec.empty()) )
  {
    double upperSize = newtonCorrectionPtr->globalLength();
    int sumSize = innerErrorInfoVec.size();

    double totalSize = upperSize;
    double totalSum = errorNorm*errorNorm*upperSize;
    double totalQSum = qErrorNorm*qErrorNorm*upperSize;

    for (int i=0;i<sumSize;++i)
    {
      double innerSum = innerErrorInfoVec[i].xErrorSum;
      double innerQSum = innerErrorInfoVec[i].qErrorSum;
      double innerSize = innerErrorInfoVec[i].innerSize;

      if (DEBUG_TIME && isActive(Diag::TIME_ERROR))
      {
        Xyce::dout() << "DSdae:innerSum["<<i<<"] = " << innerSum <<std::endl;
        Xyce::dout() << "DSdae:innerQSum["<<i<<"] = " << innerQSum <<std::endl;
        Xyce::dout() << "DSdae:innerSize["<<i<<"] = " << innerSize <<std::endl;
      }

      totalSize += innerSize;
      totalSum += innerSum;
      totalQSum += innerQSum;
    }

    double recip = 1.0/totalSize;
    errorNorm = sqrt(recip*totalSum);
    qErrorNorm = sqrt(recip*totalQSum);

    if (DEBUG_TIME && isActive(Diag::TIME_ERROR))
    {
      Xyce::dout() << "DSdae:upperSize = " << upperSize << std::endl;
      Xyce::dout() << "DSdae:totalSum = " << totalSum << std::endl;
      Xyce::dout() << "DSdae:totalQSum = " << totalQSum << std::endl;
      Xyce::dout() << "DSdae:totalSize = " << totalSize << std::endl;
      Xyce::dout() << "DSdae:2-level errorNorm = " << errorNorm << std::endl;
      Xyce::dout() << "DSdae:2-level qErrorNorm = " << qErrorNorm << std::endl;
    }
  }

  return errorNorm;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::partialQErrorNormSum
// Purpose       : Needed by 2-level solves.  This is the Q-vector version
//                 of this function.
//
// Special Notes : A weighted RMS norm is this:
//
//                  norm = sqrt ( 1/n * sum (x/w)^2 )
//
//                  What this function returns is: sum (x/w)^2
//
//                  It will later be summed with other partial sums
//                  to get the complete WRMS value.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
double DataStore::partialQErrorNormSum ()
{
  double qErrorNorm = 0.0;
  qNewtonCorrectionPtr->wRMSNorm(*qErrWtVecPtr, &qErrorNorm);
  double sum = qErrorNorm*qErrorNorm;
  double length = qNewtonCorrectionPtr->globalLength();
  sum *= length;

  return sum;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::partialSum_m1
// Purpose       : Needed by 2-level solves.
//
// Special Notes : A weighted RMS norm is this:
//
//                  norm = sqrt ( 1/n * sum (x/w)^2 )
//
//                  What this function returns is: sum (x/w)^2
//
//                  It will later be summed with other partial sums
//                  to get the complete WRMS value.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double DataStore::partialSum_m1 (int currentOrder)
{
  double sum = 0.0;

  if (currentOrder>1)
  {
    if (!delta_x)
      delta_x = builder_.createVector();

    delta_x->linearCombo(1.0,*(xHistory[currentOrder]),1.0,*newtonCorrectionPtr);
    double norm = 0.0;
    delta_x->wRMSNorm(*errWtVecPtr, &norm);
    sum = norm*norm;
    double length = newtonCorrectionPtr->globalLength();
    sum *= length;
  }

  return sum;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::partialSum_p1
// Purpose       : Needed by 2-level solves.
//
// Special Notes : A weighted RMS norm is this:
//
//                  norm = sqrt ( 1/n * sum (x/w)^2 )
//
//                  What this function returns is: sum (x/w)^2
//
//                  It will later be summed with other partial sums
//                  to get the complete WRMS value.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double DataStore::partialSum_p1 (int currentOrder, int maxOrder)
{
  double sum = 0.0;

  if (currentOrder<maxOrder)
  {
    if (!delta_x)
      delta_x = builder_.createVector();

    delta_x->linearCombo(1.0,*newtonCorrectionPtr,-1.0,*(xHistory[currentOrder+1]));
    double norm = 0.0;
    delta_x->wRMSNorm(*errWtVecPtr, &norm);
    sum = norm*norm;
    double length = newtonCorrectionPtr->globalLength();
    sum *= length;
  }

  return sum;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::partialSum_q1
// Purpose       : Needed by 2-level solves.
//
// Special Notes : A weighted RMS norm is this:
//
//                  norm = sqrt ( 1/n * sum (x/w)^2 )
//
//                  What this function returns is: sum (x/w)^2
//
//                  It will later be summed with other partial sums
//                  to get the complete WRMS value.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double DataStore::partialSum_q1 ()
{
  double sum = 0.0;

  double norm = 0.0;
  (qHistory[1])->wRMSNorm(*qErrWtVecPtr, &norm);

  sum = norm*norm;
  double length = qNewtonCorrectionPtr->globalLength();
  sum *= length;

  return sum;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::delta_x_errorNorm_q1
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double DataStore::delta_x_errorNorm_q1()
{
  double errorNorm = 0.0;
  (qHistory[1])->wRMSNorm(*qErrWtVecPtr, &errorNorm);

  // This is for handling the total errorNorm for 2-level solves.
  // It really should get a separate overloaded function.
  if ( !(innerErrorInfoVec.empty()) )
  {
    double upperSize = newtonCorrectionPtr->globalLength();
    int sumSize = innerErrorInfoVec.size();

    double totalSize = upperSize;
    double totalSum = errorNorm*errorNorm*upperSize;

    for (int i=0;i<sumSize;++i)
    {
      double innerSum = innerErrorInfoVec[i].q1HistorySum;
      double innerSize = innerErrorInfoVec[i].innerSize;

      totalSize += innerSize;
      totalSum += innerSum;
    }

    double recip = 1.0/totalSize;
    errorNorm = sqrt(recip*totalSum);
  }

  return errorNorm;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::stepLinearCombo
// Purpose       : setup the newtonCorrection vectors.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 02/20/07
//-----------------------------------------------------------------------------
void DataStore::stepLinearCombo()
{
  // 03/16/04 tscoffe:  update the newton correction.  Note:  this should be
  // available from NOX, but for now I'm going to do the difference anyway.
  newtonCorrectionPtr->linearCombo (1.0,*nextSolutionPtr,-1.0,*xn0Ptr);

  // We need to compute the correction in Q here
  // I'm assuming dsDaePtr_->daeQVectorPtr will be fresh from the end of the
  // nonlinear solve.
  qNewtonCorrectionPtr->linearCombo (1.0,*daeQVectorPtr,-1.0,*qn0Ptr);

  if (DEBUG_TIME && isActive(Diag::TIME_ERROR))
  {
    Xyce::dout() << "\n newtonCorrection: \n" << std::endl;
    newtonCorrectionPtr->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
    Xyce::dout() << "\n qNewtonCorrection: \n" << std::endl;
    qNewtonCorrectionPtr->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : DataStore::getSolnVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/27/06
//-----------------------------------------------------------------------------
bool DataStore::getSolnVarData(
  const int &           gid,
  std::vector<double> & varData )
{
  varData.resize(getNumSolnVarData());
  int i=0;
  varData[i++] = tmpSolVectorPtr->getElementByGlobalIndex( gid );
  varData[i++] = currSolutionPtr->getElementByGlobalIndex( gid );
  varData[i++] = lastSolutionPtr->getElementByGlobalIndex( gid );
  varData[i++] = nextSolutionPtr->getElementByGlobalIndex( gid );
  varData[i++] = errWtVecPtr->getElementByGlobalIndex( gid );
  varData[i++] = qErrWtVecPtr->getElementByGlobalIndex ( gid );
  varData[i++] = daeQVectorPtr->getElementByGlobalIndex ( gid );
  varData[i++] = daeFVectorPtr->getElementByGlobalIndex ( gid );
  varData[i++] = daeBVectorPtr->getElementByGlobalIndex ( gid );
  varData[i++] = dFdxdVpVectorPtr->getElementByGlobalIndex ( gid );
  varData[i++] = dQdxdVpVectorPtr->getElementByGlobalIndex ( gid );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::getStateVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/27/06
//-----------------------------------------------------------------------------
bool DataStore::getStateVarData( const int & gid,
			 	       std::vector<double> & varData )
{
  int i=0;
  varData.resize(getNumStateVarData());
  varData[i++] = tmpStaVectorPtr->getElementByGlobalIndex( gid );
  varData[i++] = tmpStaDerivPtr->getElementByGlobalIndex( gid );
  varData[i++] = currStatePtr->getElementByGlobalIndex( gid );
  varData[i++] = lastStatePtr->getElementByGlobalIndex( gid );
  varData[i++] = nextStatePtr->getElementByGlobalIndex( gid );
  varData[i++] = currStateDerivPtr->getElementByGlobalIndex( gid );
  varData[i++] = nextStateDerivPtr->getElementByGlobalIndex( gid );
  return true;
}


//-----------------------------------------------------------------------------
// Function      : DataStore::getStoreVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date :
//-----------------------------------------------------------------------------
bool DataStore::getStoreVarData( const int & gid,
			 	       std::vector<double> & varData )
{
  int i=0;
  varData.resize(getNumStoreVarData());
  varData[i++] = tmpStoVectorPtr->getElementByGlobalIndex( gid );
  varData[i++] = currStorePtr->getElementByGlobalIndex( gid );
  varData[i++] = nextStorePtr->getElementByGlobalIndex( gid );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::setSolnVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/27/06
//-----------------------------------------------------------------------------
bool DataStore::setSolnVarData( const int & gid,
				      const std::vector<double> & varData )
{
  int i=0;
  tmpSolVectorPtr->setElementByGlobalIndex       ( gid, varData[i++] );
  currSolutionPtr->setElementByGlobalIndex       ( gid, varData[i++] );
  lastSolutionPtr->setElementByGlobalIndex       ( gid, varData[i++] );
  nextSolutionPtr->setElementByGlobalIndex       ( gid, varData[i++] );
  errWtVecPtr->setElementByGlobalIndex           ( gid, varData[i++] );
  qErrWtVecPtr->setElementByGlobalIndex          ( gid, varData[i++] );
  daeQVectorPtr->setElementByGlobalIndex         ( gid, varData[i++] );
  daeFVectorPtr->setElementByGlobalIndex         ( gid, varData[i++] );
  daeBVectorPtr->setElementByGlobalIndex         ( gid, varData[i++] );
  dFdxdVpVectorPtr->setElementByGlobalIndex      ( gid, varData[i++] );
  dQdxdVpVectorPtr->setElementByGlobalIndex      ( gid, varData[i++] );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::setStateVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/27/06
//-----------------------------------------------------------------------------
bool DataStore::setStateVarData( const int & gid,
			 	       const std::vector<double> & varData )
{
  int i=0;
  tmpStaVectorPtr->setElementByGlobalIndex    ( gid, varData[i++] );
  tmpStaDerivPtr->setElementByGlobalIndex     ( gid, varData[i++] );
  currStatePtr->setElementByGlobalIndex       ( gid, varData[i++] );
  lastStatePtr->setElementByGlobalIndex       ( gid, varData[i++] );
  nextStatePtr->setElementByGlobalIndex       ( gid, varData[i++] );
  currStateDerivPtr->setElementByGlobalIndex  ( gid, varData[i++] );
  nextStateDerivPtr->setElementByGlobalIndex  ( gid, varData[i++] );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DataStore::setStoreVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/27/06
//-----------------------------------------------------------------------------
bool DataStore::setStoreVarData( const int & gid,
			 	       const std::vector<double> & varData )
{
  int i=0;
  tmpStoVectorPtr->setElementByGlobalIndex    ( gid, varData[i++] );
  currStorePtr->setElementByGlobalIndex       ( gid, varData[i++] );
  nextStorePtr->setElementByGlobalIndex       ( gid, varData[i++] );
  return true;
}

} // namespace TimeIntg
} // namespace Xyce

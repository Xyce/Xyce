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


//-------------------------------------------------------------------------
//
// Purpose        : This is the test program for the time integration 
//                  pacakge.
//
// Special Notes  :
//
// Creator        : Buddy Watts
//
// Creation Date  : 5/26/00
//
//-------------------------------------------------------------------------

// ----------   Xyce Includes   ----------
#include <TIATest.h>

#ifdef TESTLA
  #include <N_LAS_Dummy.h>
#else
  #include <N_LAS_LAFactory.h>
  #include <N_LAS_System.h>
  #include <N_LAS_MultiVector.h>
  #include <N_LAS_Matrix.h>
#endif

#include <N_NLS_Manager.h>
#include <N_NLS_NonLinearSolver.h>

#include <N_DEV_DeviceMgr.h>

#include <N_LOA_LoaderMgr.h>
#include <N_LOA_Loader.h>

//-----------------------------------------------------------------------------
// Function      : TIATestor::setTiaParams
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/09/00
//-----------------------------------------------------------------------------
bool TIATestor::setTiaParams ()
{
  tiaParams_.RestartingIntegrationFromSSS = 0;
  tiaParams_.InitialTime = 1.0;
  tiaParams_.FinalTime = 5.0;
  tiaParams_.IntegrationMethod = 2;
  tiaParams_.StartingTimeStep = 0.5;
  tiaParams_.ConstantStepSize = 1;
  tiaParams_.solutionSize = 3;
  tiaParams_.stateSize = 3;

  //tiaParams_.xxx = 1.2;      //      Initial_Solution_(1)
  //tiaParams_.xxx = 3.4;      //      Initial_Solution_(2)
  //tiaParams_.xxx = 5.6;      //      Initial_Solution_(3)

  tiaParams_.ScalarTolerances = 1;
  tiaParams_.IntegRelErrorTol = 1.0e-4;
  tiaParams_.IntegAbsErrorThreshold = 1.0e-8;
  tiaParams_.NumDiscontinuityPts = 1;

  //tiaParams_.xxx = 3.5;      //      PtOfDiscontinuity

  return STATUS_SUCCESS;
}

//-----------------------------------------------------------------------------
// Function      : TIATestor::setNLParams
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/09/00
//-----------------------------------------------------------------------------
bool TIATestor::setNLParams ()
{
  nlParamsPtr_->deltaxTol     = 1.0e-8; 
  nlParamsPtr_->normRHS       = 0.0; 
  nlParamsPtr_->absTol        = 1.0e-8; 
  nlParamsPtr_->relTol        = 1.0e-8; 
  nlParamsPtr_->dampFactor    = 1.0; 
  nlParamsPtr_->maxChange     = 0.0; 
  nlParamsPtr_->newtonStep    = 0; 
  nlParamsPtr_->maxNewtonStep = 15; 
  nlParamsPtr_->dampStep      = 0; 
  nlParamsPtr_->maxDampStep   = 7; 
  nlParamsPtr_->normLevel     = 2; 

  return STATUS_SUCCESS;
}



//-----------------------------------------------------------------------------
// Function      : TIATestor::doAllocations
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/06/00
//-----------------------------------------------------------------------------
bool TIATestor::doAllocations ()
{
  bool bsuccess = STATUS_SUCCESS;

  lasSysPtr_ = new N_LAS_System ();
  nlsMgrPtr_ = new N_NLS_Manager ();
  nlsMgrPtr_->createNLS();
  nlsPtr_ = nlsMgrPtr_->getNonLinearSolver ();

  nlParamsPtr_ = new N_NLS_NLParams ();

  bsuccess = bsuccess && (lasSysPtr_ != NULL);
  bsuccess = bsuccess && (nlsMgrPtr_ != NULL);
  bsuccess = bsuccess && (nlsPtr_ != NULL);

  // do the loader:
  loaderMgrPtr_ = new N_LOA_LoaderMgr (); 
  loaderPtr_ = loaderMgrPtr_->createLoader (_CKTLOADER);

  bsuccess = bsuccess &&  (loaderMgrPtr_ != NULL);
  bsuccess = bsuccess &&  (loaderPtr_ != NULL);

  // device  manager:
  devPtr_ = N_DEV_DeviceMgr::factory ();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : TIATestor::doRegistrations
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/06/00
//-----------------------------------------------------------------------------
bool TIATestor::doRegistrations ()
{
  bool bsuccess = STATUS_SUCCESS;

  // time integration registrations:
  bsuccess = bsuccess && tia_.registerTIAParams    (tiaParams_);
  bsuccess = bsuccess && tia_.registerLinearSystem (lasSysPtr_);
  bsuccess = bsuccess && tia_.registerNLSolver     (nlsPtr_);
  bsuccess = bsuccess && tia_.registerLoader       (loaderPtr_);

  // nonlinear solver registrations:
  bsuccess = bsuccess && nlsMgrPtr_->registerNLParams    (*nlParamsPtr_);
  bsuccess = bsuccess && nlsMgrPtr_->registerLoader       (loaderPtr_);
  bsuccess = bsuccess && nlsMgrPtr_->registerLinearSystem (lasSysPtr_);
  bsuccess = bsuccess && nlsMgrPtr_->registerTimeIntegrator (&tia_);

  // loader registrations:
  bsuccess = bsuccess && loaderMgrPtr_->registerDeviceManager (devPtr_);

  // device registrations:
  bsuccess = bsuccess && devPtr_->registerLinearSystem (lasSysPtr_);
  bsuccess = bsuccess && devPtr_->registerTimeIntegrator (&tia_);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : TIATestor::doMatrixCreation
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/12/00
//-----------------------------------------------------------------------------
bool TIATestor::doMatrixCreation ()
{
  bool bsuccess = STATUS_SUCCESS;

  // jacobian matrix was already created in the dummy version of the
  // N_LAS_System constructor.

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : TIATestor::doInitializations
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/12/00
//-----------------------------------------------------------------------------
bool TIATestor::doInitializations ()
{
  bool bsuccess = STATUS_SUCCESS;

  bsuccess = bsuccess &&     tia_.initializeAll ();
  bsuccess = bsuccess && devPtr_->initializeAll ();
  bsuccess = bsuccess && nlsMgrPtr_->initializeAll ();
 
  return bsuccess; 
}

//-----------------------------------------------------------------------------
// Function      : TIATestor::doDeAllocations
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/12/00
//-----------------------------------------------------------------------------
bool TIATestor::doDeAllocations ()
{
  bool bsuccess = STATUS_SUCCESS;

  if (lasSysPtr_    != NULL) delete lasSysPtr_;
  if (nlsMgrPtr_    != NULL) delete nlsMgrPtr_;
  if (loaderMgrPtr_ != NULL) delete loaderMgrPtr_;
  if (loaderPtr_    != NULL) delete loaderPtr_;
  if (devPtr_       != NULL) delete devPtr_;
  if (nlParamsPtr_  != NULL) delete nlParamsPtr_;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : TIATestor::runTests
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/06/00
//-----------------------------------------------------------------------------
bool TIATestor::runTests(int iargs, char  *cargs[]) 
{
  bool bsuccess;

  bsuccess = doAllocations ();  

  setTiaParams ();
  setNLParams  ();

  bsuccess = bsuccess && doRegistrations ();

  bsuccess = bsuccess && doMatrixCreation ();

  bsuccess = bsuccess && doInitializations ();

  bsuccess = bsuccess && tia_.runIntegrator ();

  bsuccess = bsuccess && doDeAllocations ();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : main
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/06/00
//-----------------------------------------------------------------------------
int main(int iargs, char  *cargs[]) 
{
  TIATestor * TIAT_ptr = new TIATestor ();
  bool bsuccess = TIAT_ptr->runTests(iargs, cargs);
  delete TIAT_ptr;
  return 0;
}

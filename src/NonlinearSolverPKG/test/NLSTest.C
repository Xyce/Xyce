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
// File          : NLSTest.C
//
// Purpose       : This function is the test driver for the nonlinear solver
//                 package.
//
// Special Notes :
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------

// ---------- Standard Includes ----------
#include <iostream>
#include <vector>
#include <list>
#include <string>

// ----------   Xyce Includes   ----------
#include <NLSTest.h>
#include <N_LOA_LoaderMgr.h>
#include <N_LOA_Loader.h>

#ifdef TESTLA
#include <N_LAS_Dummy.h>
#else
#include <N_LAS_LAFactory.h>
#endif

#ifdef TESTTI
#include <N_TIA_Dummy.h>
#else
#include <N_TIA_TimeIntegrationAlgorithm.h>
#endif

//-----------------------------------------------------------------------------
// Function      : NLSTestor::NLSTestor
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
NLSTestor::NLSTestor()
{

}

//-----------------------------------------------------------------------------
// Function      : NLSTestor::~NLSTestor
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
NLSTestor::~NLSTestor()
{

}

//-----------------------------------------------------------------------------
// Function      : NLSTestor::doAllocations
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
bool NLSTestor::doAllocations()
{
  int bsuccess = STATUS_SUCCESS;

  cout << endl;
  cout << "Performing allocations" << endl;

  LAS_SolverPtr_ = new N_LAS_IterativeSolver(); // Allocate LAS class
  LAS_MatrixPtr_ = new N_LAS_Matrix();          // Allocate LAS matrix class
  LAS_RHSVecPtr_ = new N_LAS_MultiVector();     // Allocate LAS vector class
  LAS_SolVecPtr_ = new N_LAS_MultiVector();     // Allocate LAS vector class
  LOA_LoaderMgrPtr_ = new N_LOA_LoaderMgr();    // Allocate loader manager.

  // Allocate Loader
  LOA_LoaderPtr_ = LOA_LoaderMgrPtr_->createLoader(_CKTLOADER);

  // Allocate Time Integration:
  TIA_Ptr_       = new N_TIA_TimeIntegrationAlgorithm();

  ERH_Ptr_       = new N_ERH_ErrorMgr();       // Allocate Error handler
  NLS_Ptr_       = new N_NLS_Manager();        // Allocate Nonlinear Solver
  DEV_Ptr_       = N_DEV_DeviceMgr::factory(); // Allocate device manager

  bsuccess = bsuccess && (LAS_SolverPtr_    != NULL);
  bsuccess = bsuccess && (LAS_MatrixPtr_    != NULL);
  bsuccess = bsuccess && (LAS_RHSVecPtr_    != NULL);
  bsuccess = bsuccess && (LAS_SolVecPtr_    != NULL);
  bsuccess = bsuccess && (LOA_LoaderMgrPtr_ != NULL);
  bsuccess = bsuccess && (LOA_LoaderPtr_    != NULL);
  bsuccess = bsuccess && (TIA_Ptr_          != NULL);
  bsuccess = bsuccess && (ERH_Ptr_          != NULL);
  bsuccess = bsuccess && (DEV_Ptr_          != NULL);

  cout << "Done with allocations";
  cout << endl;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : NLSTestor::doRegistrations
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
bool NLSTestor::doRegistrations()

{
  bool bsuccess = STATUS_SUCCESS;

  cout << endl;
  cout << "Performing registrations";
  cout << endl;

  bsuccess = bsuccess && NLS_Ptr_->registerLinearSolver(LAS_SolverPtr_);
  bsuccess = bsuccess && NLS_Ptr_->registerLinearSystem(LAS_MatrixPtr_,
                                                        LAS_SolVecPtr_,
                                                        LAS_RHSVecPtr_);

  bsuccess = bsuccess && NLS_Ptr_->registerLoader(LOA_LoaderPtr_);

  bsuccess = bsuccess && NLS_Ptr_->registerTimeIntegrator(TIA_Ptr_);

  bsuccess = bsuccess && LOA_LoaderMgrPtr_->registerDeviceManager(DEV_Ptr_);

  cout << "Done with registrations";
  cout << endl;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : NLSTestor::doDeAllocations
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
bool NLSTestor::doDeAllocations()
{
  bool bsuccess = STATUS_SUCCESS;

  cout << endl;
  cout << "Performing de-allocations";
  cout << endl;

  delete LAS_SolverPtr_;
  delete LAS_MatrixPtr_;
  delete LAS_RHSVecPtr_;
  delete LAS_SolVecPtr_;
  delete LOA_LoaderMgrPtr_;
  delete LOA_LoaderPtr_;
  delete TIA_Ptr_;
  delete ERH_Ptr_;
  delete NLS_Ptr_;
  delete DEV_Ptr_;

  cout << "Done with de-allocations";
  cout << endl;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : NLTestor::doInitialization
// Purpose       : Tells the nonlinear solver to allocate one of the various
//                 nonlinear solvers.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
bool NLSTestor::doInitialization()
{
  int isolver;
  N_NLS_NLParams nlParams;

  cout << "Performing initializations" << endl;

  // Which solver is allocated depends upon the command line args.

  bool bsuccess = NLS_Ptr_->createNLS();
  bool btmp     = NLS_Ptr_->initNLS(nlParams);

  cout << "Done Initializations" << endl;
  cout << endl;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : NLSTestor::doSolve
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
bool NLSTestor::doSolve()
{
  cout << "Performing solve" << endl;

  bool bsuccess = NLS_Ptr_->solve();

  cout << "Done with doSolve" << endl;
  cout << endl;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : NLSTestor::runTests
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
bool NLSTestor::runTests(int iargs_tmp, char *cargs_tmp[])
{
  bool bsuccess;
  bool btmp;

  iargs = iargs_tmp;   cargs = cargs_tmp;

  cout << endl;
  cout << "Welcome to the Nonlinear solver testing program";
  cout << endl;

  bsuccess = doAllocations();
  btmp     = doRegistrations();  bsuccess = bsuccess && btmp;
  btmp     = doInitialization(); bsuccess = bsuccess && btmp;
  btmp     = doSolve();          bsuccess = bsuccess && btmp;
  btmp     = doDeAllocations();  bsuccess = bsuccess && btmp;

  cout << endl;
  cout << "Done testing the Nonlinear solver";
  cout << endl;
  cout << endl;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : main
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
int main(int iargs, char *cargs[])
{
  NLSTestor *NLSPtr = new NLSTestor();

  bool bsuccess = NLSPtr->runTests(iargs, cargs);

  delete NLSPtr;

  return 0;
}


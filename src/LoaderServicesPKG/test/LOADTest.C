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
// File          : LOADTest.C
//
// Purpose       : This function is the test driver for the Loader Services
//                 package.
//
// Special Notes :
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------

// ---------- Standard Includes ----------
#include <iostream>
#include <vector>
#include <list>
#include <string>

// ----------   Xyce Includes   ----------
#include <LOADTest.h>
#include <N_LOA_LoaderMgr.h>

//-----------------------------------------------------------------------------
// Function      : LOADTestor::LOADTestor
// Purpose       : constructor
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
LOADTestor::LOADTestor ()
{

}

//-----------------------------------------------------------------------------
// Function      : LOADTestor::~LOADTestor
// Purpose       : destructor
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
LOADTestor::~LOADTestor ()
{

}

//-----------------------------------------------------------------------------
// Function      : LOADTestor::doAllocations
// Purpose       : 
// Special Notes : 
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool LOADTestor::doAllocations ()
{
  bool bsuccess = STATUS_SUCCESS;

  cout << endl;
  cout << "Performing allocations" << endl;

  //ERH_Ptr_          = new N_ERH_ErrorMgr  ();     // Allocate Error handler
  LOA_LoaderMgrPtr_ = new N_LOA_LoaderMgr ();     // Allocate Loader manager
  DEV_DeviceMgrPtr_ = N_DEV_DeviceMgr::factory(); // Allocate Device manager.

  //bsuccess = bsuccess && (ERH_Ptr_          != NULL);
  bsuccess = bsuccess && (LOA_LoaderMgrPtr_ != NULL);
  bsuccess = bsuccess && (DEV_DeviceMgrPtr_ != NULL);

  cout << "Done with allocations";
  cout << endl;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : LOADTestor::doRegistrations
// Purpose       :
// Special Notes : 
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool LOADTestor::doRegistrations ()
{
  bool bsuccess = STATUS_SUCCESS;

  cout << endl;
  cout << "Performing registrations";
  cout << endl;

  LOA_LoaderMgrPtr_->registerDeviceManager (DEV_DeviceMgrPtr_);

  cout << "Done with registrations";
  cout << endl;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : LOADTestor::doDeAllocations
// Purpose       : 
// Special Notes : 
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool LOADTestor::doDeAllocations ()
{
  bool bsuccess = STATUS_SUCCESS;

  cout << endl;
  cout << "Performing de-allocations";
  cout << endl;

  delete LOA_LoaderMgrPtr_;
  delete DEV_DeviceMgrPtr_;
  //delete ERH_Ptr_;

  cout << "Done with de-allocations";
  cout << endl;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : LOADTestor::doInitialization
// Purpose       : Tells the loader manager to allocate one of the
//                 various loaders.
// Special Notes : 
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool LOADTestor::doInitialization ()
{
  bool bsuccess = STATUS_SUCCESS;

  cout << "Performing initializations" << endl;

  LOA_LoaderMgrPtr_->createLoader(0);

  cout << "Done Initializations" << endl;
  cout << endl;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : LOADTestor::doLoad
// Purpose       : Tells the loader manager to perform a RHS and a Jacobian
//                 load.
// Special Notes : 
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool LOADTestor::doLoad ()
{
  bool bsuccess = STATUS_SUCCESS;
  bool btmp;

  cout << "Performing load tests:\n";

  btmp = LOA_LoaderMgrPtr_->loadJacobian (); bsuccess = bsuccess && btmp;
  cout << "Done with Jacobian\n";
  btmp = LOA_LoaderMgrPtr_->loadRHS      (); bsuccess = bsuccess && btmp;
  cout << "Done with RHS\n";

  cout << "done performing load tests:\n";

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : LOADTestor::runTests
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool LOADTestor::runTests (int iargs_tmp, char *cargs_tmp[])
{
  bool bsuccess;
  bool btmp;

  iargs = iargs_tmp;   cargs = cargs_tmp;

  cout << endl;
  cout << "Welcome to the Loader Services testing program";
  cout << endl;

  bsuccess = doAllocations     ();
  btmp     = doRegistrations   (); bsuccess = bsuccess && btmp;
  btmp     = doInitialization  (); bsuccess = bsuccess && btmp;
  btmp     = doLoad            (); bsuccess = bsuccess && btmp;
  btmp     = doDeAllocations   (); bsuccess = bsuccess && btmp;

  cout << endl;
  cout << "Done testing the Loader Services";
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
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
int main (int iargs, char *cargs[])
{
  LOADTestor *LOADPtr = new LOADTestor ();

  bool bsuccess = LOADPtr->runTests (iargs, cargs);

  delete LOADPtr;

  return 0;
}


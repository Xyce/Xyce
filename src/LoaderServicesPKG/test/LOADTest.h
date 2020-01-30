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
// File          : LOADTest.h
//
// Purpose       : This function is the header file which contains class
//                 definitions for the nonlinear solver package test 
//                 program.
//
// Special Notes :
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------


#ifndef  _LOADTEST_H
#define  _LOADTEST_H

// ---------- Standard Includes ----------
#include <iostream>
#include <vector>
#include <list>
#include <string>

// ----------   Xyce Includes   ----------
#include <N_LOA_LoaderMgr.h>
#include <N_DEV_DeviceMgr.h>
#include <N_ERH_ErrorMgr.h>

using namespace std;

//-----------------------------------------------------------------------------
// Class         : LOADTestor
// Purpose       : This is the top level class for the Loader Services
//                 testing program.  The member function, RunTests, 
//                 is the "main" function, essentially.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
class LOADTestor
{
  // functions:
  public:
    LOADTestor                ();
    ~LOADTestor               ();

    bool  runTests           (int iargs, char *cargs[]);

  protected:

  private:
    bool  doAllocations      ();
    bool  doRegistrations    ();
    bool  doDeAllocations    ();

    bool  doInitialization   ();
    bool  doLoad             ();

  // attributes
  public:

  protected:

  private:
    N_DEV_DeviceMgr   * DEV_DeviceMgrPtr_;
    N_LOA_LoaderMgr   * LOA_LoaderMgrPtr_;
    N_ERH_ErrorMgr    * ERH_Ptr_;

    int iargs;
    char **cargs;
};

#endif



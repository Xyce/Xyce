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
// Purpose        : This is the test program for the time integration 
//                  pacakge.
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 6/06/00
//
//-------------------------------------------------------------------------

// ---------- Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_TIA_TimeIntegrationAlgorithm.h>
#include <N_ERH_ErrorMgr.h>

// ---------- Forward Declarations -----
class N_DEV_DeviceMgr;
class N_LOA_LoaderMgr;
class N_LOA_Loader;
class N_LAS_Solver;
class N_LAS_System;
class N_NLS_Manager;
class N_NLS_NonLinearSolver;
class N_NLS_NLParams;

//-----------------------------------------------------------------------------
// Class         : TIATestor
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/06/00
//-----------------------------------------------------------------------------
class TIATestor
{
  public:
    TIATestor   () {};
    ~TIATestor  () {};

    bool setTiaParams      ();
    bool setNLParams       ();
    bool doAllocations     ();
    bool doRegistrations   ();
    bool doInitializations ();
    bool doMatrixCreation  ();
    bool doDeAllocations   ();

    bool runTests (int iargs, char *cargs[]);

  protected:
  private:

  public:
  protected:
  private:
    N_LAS_System           * lasSysPtr_; 
    N_NLS_Manager          * nlsMgrPtr_; 
    N_NLS_NonLinearSolver  * nlsPtr_;
    N_LOA_LoaderMgr        * loaderMgrPtr_;
    N_LOA_Loader           * loaderPtr_;
    N_DEV_DeviceMgr        * devPtr_;

    N_TIA_TimeIntegrationAlgorithm tia_;
    N_TIA_TIAParams tiaParams_;
    N_NLS_NLParams * nlParamsPtr_;
};




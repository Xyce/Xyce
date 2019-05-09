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


////////////////////////////////////////////////////////////////////////////
// File:        DeviceTestor.h
// Author:      Eric Keiter
// Description:
//
//

#ifndef  _DEVTEST_H
#define  _DEVTEST_H

// ---------- Standard Includes ----------
#include <iostream>
#include <vector>
#include <list>
#include <string>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_ERH_ErrorMgr.h>

// ---------- Forward Declarations ----------
class N_LAS_Solver;
class N_LAS_MultiVector;
class N_LAS_Matrix;
class N_TIA_TimeIntegrationAlgorithm;

class N_LAS_System;
class N_LAS_LAFactory;

class N_PDS_Manager;
class N_PDS_LoadBalance;
class N_PDS_Comm;
class N_PDS_ParMap;

enum ElementParamTest {
    _RESISTOR_ELEMENT_TEST,   //
    _CAPACITOR_ELEMENT_TEST,  //
    _INDUCTOR_ELEMENT_TEST,   //
    _VSRC_ELEMENT_TEST,       //
    _ISRC_ELEMENT_TEST,       //
    _NUM_ELEMENT_TESTS
};

enum ModelParamTest {
    _DIODE_MODEL_TEST,   //
    _NUM_MODEL_TESTS
};

using namespace std;

//-----------------------------------------------------------------------------
// Class         : DeviceTestor
// Purpose       : This is the top level class for the device 
//                 testing program.  The member function, runTests, 
//                 is the "main" function, essentially.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class DeviceTestor
{
  // functions:
  public:
    DeviceTestor                ();
    ~DeviceTestor               ();

    bool  setupParMgr           ();
    bool  doAllocations         ();
    bool  doRegistrations       ();
    bool  doInitializations     ();
    bool  doDeAllocations       ();
    bool  createDevices         ();
    bool  deleteDevices         ();
    bool  addModels             ();
    bool  addInstances          ();
    bool  getPointers           ();
    bool  getTopologies         ();
    bool  performRHSLoads       ();
    bool  performJacLoads       ();

    bool  getDefaultElementInfo ();
    bool  getDefaultModelInfo   ();

    REAL  getValue              (const string & str) const;
    void  loadElementFields     (vector<string> & str,
                                 const int iElement);

    void  loadModelFields       (vector<string> & str,
                                 const int iModel);

    bool  outputMI              ();

    int   runTests              (int iargs, char *cargs[]);

  protected:

  private:


  // attributes
  public:

  protected:

  private:
    Xyce::Device::DeviceMgr *           DevMgrPtr_;

    //N_LAS_Matrix                   * LAS_MatrixPtr_;
    //N_LAS_MultiVector              * LAS_RHSVecPtr_;
    //N_LAS_MultiVector              * LAS_SolVecPtr_;
    //N_LAS_MultiVector              * LAS_TmpSolVecPtr_;

    N_LAS_MultiVector * currStatePtr;
    N_LAS_MultiVector * currSolutionPtr;
    N_LAS_MultiVector * nextStatePtr;
    N_LAS_MultiVector * nextSolutionPtr;
    N_LAS_MultiVector * tmpSolVectorPtr;
    N_LAS_MultiVector * tmpStaVectorPtr;
    N_LAS_MultiVector * errorEstimatePtr;
    N_LAS_MultiVector * nextSolutionDerivPtr;

    N_LAS_System                   * lasSysPtr_;
    N_LAS_LAFactory                * lasLAFactoryPtr_;

    N_PDS_Manager                  * parMgrPtr_;
    N_PDS_ParMap                   * parMapPtr_;
    N_PDS_LoadBalance              * lasLBPtr_;
    N_PDS_Comm                     * PDSCommPtr_;


    N_TIA_TimeIntegrationAlgorithm * TIAPtr_;

};

#endif

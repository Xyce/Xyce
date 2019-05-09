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
// File          : DMTest.C
//
// Purpose       : This function is the test driver for the DeviceModel
//                 package.
//
// Special Notes :
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date : 3/20/00
//-----------------------------------------------------------------------------


// ---------- Standard Includes ----------
#include <iostream>
#include <vector>
#include <list>
#include <string>

// ----------   Xyce Includes   ----------
#include <DMTest.h>

#include <N_DEV_DeviceMgr.h>
#include <N_DEV_Device.h>
#include <N_DEV_DeviceBlock.h>

#ifdef TESTTI
  #include <N_TIA_Dummy.h>
#else
  #include <N_TIA_TimeIntegrationAlgorithm.h>
#endif

#ifdef TESTLA
  #include <N_LAS_Dummy.h>
#else
  #include <N_LAS_LAFactory.h>
  #include <N_LAS_Matrix.h>
  #include <N_LAS_MultiVector.h>
  #include <N_LAS_System.h>
  #include <N_LAS_LAFactory.h>
#endif

#include <N_PDS_Manager.h>
#include <N_PDS_LoadBalance.h>
#include <N_PDS_ParMap.h>

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::DeviceTestor
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/20/00
//-----------------------------------------------------------------------------
DeviceTestor::DeviceTestor()
{

}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::~DeviceTestor
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/20/00
//-----------------------------------------------------------------------------
DeviceTestor::~DeviceTestor()
{

}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::setupParMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------

bool DeviceTestor::setupParMgr()

{
  bool bSuccess = true;

#ifdef Xyce_PARALLEL_MPI
  static string             lbMethod = ("PARMETIS");
  static list<string_param> lbParams;
#endif

  // Setup the Parallel Mgr. with a default load-balance based on the numProc
  // value.
  parMgrPtr_ = new N_PDS_Manager(1
#ifdef Xyce_PARALLEL_MPI
                                 , lbMethod, lbParams
#endif
                                 );

  lasLBPtr_   = parMgrPtr_->getLAS_LB();
  parMapPtr_  = lasLBPtr_->getParallelMap();
  PDSCommPtr_ = parMapPtr_->getPDSComm();

  cout << "lasLBPtr_ = " << lasLBPtr_ << endl;
  cout << "parMapPtr_ = " << parMapPtr_ << endl;
  cout << "PDSCommPtr_ = " << PDSCommPtr_ << endl;

  // Register this information with the error handler.
  N_ERH_ErrorMgr::registerComm(PDSCommPtr_);

  return bSuccess;

}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::doAllocations
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/20/00
//-----------------------------------------------------------------------------
bool DeviceTestor::doAllocations()
{
  bool bsuccess = STATUS_SUCCESS;

  string topotype = "Basic";

  DevMgrPtr_ = Device::DeviceMgr::factory();  // Allocate device manager:

  // Linear Algebra allocations:
  lasLAFactoryPtr_ = new Linear::LAFactory();
  lasSysPtr_       = new Linear::System();


  currStatePtr         = lasSysPtr_->createStateVector ();
  nextStatePtr         = lasSysPtr_->createStateVector ();
  tmpStaVectorPtr      = lasSysPtr_->createStateVector ();

  currSolutionPtr      = lasSysPtr_->createVector ();
  nextSolutionPtr      = lasSysPtr_->createVector ();
  tmpSolVectorPtr      = lasSysPtr_->createVector ();
  errorEstimatePtr     = lasSysPtr_->createVector ();
  nextSolutionDerivPtr = lasSysPtr_->createVector ();

  // Allocate Time Integration
  TIAPtr_    = new N_TIA_TimeIntegrationAlgorithm();

  bsuccess = bsuccess && (DevMgrPtr_ != NULL);

  bsuccess = bsuccess && (lasLAFactoryPtr_ != NULL);
  bsuccess = bsuccess && (lasSysPtr_ != NULL);


  bsuccess = bsuccess && (currStatePtr != NULL);
  bsuccess = bsuccess && (nextStatePtr != NULL);
  bsuccess = bsuccess && (tmpStaVectorPtr != NULL);

  bsuccess = bsuccess && (currSolutionPtr != NULL);
  bsuccess = bsuccess && (nextSolutionPtr != NULL);
  bsuccess = bsuccess && (tmpSolVectorPtr != NULL);
  bsuccess = bsuccess && (errorEstimatePtr != NULL);
  bsuccess = bsuccess && (nextSolutionDerivPtr != NULL);

  bsuccess = bsuccess && (TIAPtr_ != NULL);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::doDeAllocations
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/20/00
//-----------------------------------------------------------------------------
bool DeviceTestor::doDeAllocations()
{
  // de-allocate the device manager:

  delete DevMgrPtr_;
#if 0
  delete LAS_MatrixPtr_;
  delete LAS_RHSVecPtr_;
  delete LAS_SolVecPtr_;
  delete LAS_TmpSolVecPtr_;
#else
  delete lasSysPtr_;
  delete lasLAFactoryPtr_;

  delete currStatePtr;
  delete currSolutionPtr;
  delete nextStatePtr;
  delete nextSolutionPtr;
  delete tmpSolVectorPtr;
  delete tmpStaVectorPtr;
  delete errorEstimatePtr;
  delete nextSolutionDerivPtr;

#endif
  delete TIAPtr_;

  return STATUS_SUCCESS;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::doRegistrations
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/20/00
//-----------------------------------------------------------------------------
bool  DeviceTestor::doRegistrations()
{
  bool bsuccess = STATUS_SUCCESS;

#if 0
  bsuccess = bsuccess && DevMgrPtr_->registerLinearSystem (LAS_MatrixPtr_,
                                                           &LAS_SolVecPtr_,
                                                           LAS_TmpSolVecPtr_,
                                                           LAS_RHSVecPtr_);
#else
  bsuccess = bsuccess && DevMgrPtr_->registerLinearSystem (lasSysPtr_);
#endif

  bsuccess = bsuccess && DevMgrPtr_->registerTimeIntegrator (TIAPtr_);



  // linear algebra registrations:
  // register current:
  lasSysPtr_->registerCurrStaVector(&(currStatePtr));
  lasSysPtr_->registerCurrSolVector(&(currSolutionPtr));

  // register next:
  lasSysPtr_->registerNextStaVector(&(nextStatePtr));

  cout << "About to register nextSolutionPtr = ";
  cout  << nextSolutionPtr << endl;

  lasSysPtr_->registerNextSolVector(&(nextSolutionPtr));

  // register temporaries:
  lasSysPtr_->registerTmpSolVector(tmpSolVectorPtr);
  lasSysPtr_->registerTmpStaVector(tmpStaVectorPtr);

  // register error estimate:
  lasSysPtr_->registerErrorEstVector(errorEstimatePtr);

  // register next solution derivative:
  lasSysPtr_->registerNextSolDerivVector(&(nextSolutionDerivPtr));

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::doInitializations
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
bool  DeviceTestor::doInitializations ()
{
  bool bsuccess = DevMgrPtr_->initializeAll ();
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::createDevices
// Purpose       : This function tests the creation of each type
//                  of device.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/20/00
//-----------------------------------------------------------------------------
bool DeviceTestor::createDevices()
{

  // loop over all the device types:
  int i;
  Dev_index imax = _NUMDEV;

  bool test;
  bool test_total = STATUS_SUCCESS;

  for (i=1;i<imax;i++)
  {
     test = DevMgrPtr_->createDevice(i);
     cout << "test for Device index = " << i << " is  " << test;
     test_total = test_total && test;
     cout << ".  test_total = " << test_total<< endl;
  }

  return test_total;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::deleteDevices
// Purpose       : This function tests the creation of each type
//                  of device.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/20/00
//-----------------------------------------------------------------------------
bool DeviceTestor::deleteDevices()
{
  DevMgrPtr_->flushDevices();
  return STATUS_SUCCESS;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::addModels
// Purpose       : This function tests out adding parameters defined via
//                 .model statements (or their equivalent).
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/20/00
//-----------------------------------------------------------------------------
bool DeviceTestor::addModels()
{

  Device::ModelBlock MB;
  list<tagged_param>::iterator iter;

  //--------------------------------------------------------------
  MB.name   = "Resistor";
  MB.type   = "R";
  MB.level = 0;
  tagged_param TP("R",3.0);
  MB.params.push_back(TP);

  bool is1 = DevMgrPtr_->addDeviceModel(MB);

  //--------------------------------------------------------------
  MB.name   = "BigCapacitor";
  MB.type   = "C";
  MB.level = 0;
  iter = MB.params.begin();
  iter->tag   = "C";
  iter->param = 10.0;

  bool is2= DevMgrPtr_->addDeviceModel(MB);

  //--------------------------------------------------------------
  MB.name   = "Inductomatic";
  MB.type   = "L";
  MB.level = 0;
  iter = MB.params.begin();
  iter->tag = "L";
  iter->param = 0.001;

  bool is3= DevMgrPtr_->addDeviceModel(MB);

  //--------------------------------------------------------------
  return (is1 && is2 && is3);
}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::addInstances
// Purpose       : This function tests out the addition of single instances
//                 to the device manager.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/20/00
//-----------------------------------------------------------------------------
bool DeviceTestor::addInstances()
{
  bool bsuccess = STATUS_SUCCESS;

  //--------------------------------------------------------------
  Device::InstanceBlock IB;

  IB.name = "R1";
  IB.getModelName() = "";
  tagged_param TP("R",4.0);
  IB.params.push_back(TP);

  Device::DeviceInstance * DIPtr;

  DIPtr = DevMgrPtr_->addDeviceInstance(IB);

  if (DIPtr != NULL) bsuccess = STATUS_SUCCESS;
  else               bsuccess = STATUS_FAILURE;

  //--------------------------------------------------------------
  Device::InstanceBlock IB2;

  IB2.name = "R2";
  IB2.getModelName() = "Resistor";
  tagged_param TP2("R",2.0);
  IB2.params.push_back(TP2);

  DIPtr = 0x0;
  DIPtr = DevMgrPtr_->addDeviceInstance(IB2);

  if (DIPtr != NULL) bsuccess = bsuccess && STATUS_SUCCESS;
  else               bsuccess = bsuccess && STATUS_FAILURE;

  //--------------------------------------------------------------
  Device::InstanceBlock IB3;

  IB3.name = "C2";
  IB3.getModelName() = "BigCapacitor";
  tagged_param TP3("C",4.0);
  IB3.params.push_back(TP3);

  DIPtr = 0x0;
  DIPtr = DevMgrPtr_->addDeviceInstance(IB3);

  if (DIPtr != NULL) bsuccess = bsuccess && STATUS_SUCCESS;
  else               bsuccess = bsuccess && STATUS_FAILURE;

  //--------------------------------------------------------------
  Device::InstanceBlock IB4;

  IB4.name = "C1";
  IB4.getModelName() = "";
  tagged_param TP4("C",4.0);
  IB4.params.push_back(TP4);

  DIPtr = 0x0;
  DIPtr = DevMgrPtr_->addDeviceInstance(IB4);

  if (DIPtr != NULL) bsuccess = bsuccess && STATUS_SUCCESS;
  else               bsuccess = bsuccess && STATUS_FAILURE;

  //--------------------------------------------------------------
  Device::InstanceBlock IB5;

  IB5.name = "L4";
  IB5.getModelName()   = "Inductomatic";
  tagged_param TP5("L",4.0);
  IB5.params.push_back(TP5);

  DIPtr = 0x0;
  DIPtr = DevMgrPtr_->addDeviceInstance(IB5);

  if (DIPtr != NULL) bsuccess = bsuccess && STATUS_SUCCESS;
  else               bsuccess = bsuccess && STATUS_FAILURE;

  //--------------------------------------------------------------
  Device::InstanceBlock IB6;

  IB6.name = "R6";
  IB6.getModelName() = "";
  tagged_param TP6("R",15.0);
  IB6.params.push_back(TP6);

  DIPtr = 0x0;
  DIPtr = DevMgrPtr_->addDeviceInstance(IB6);

  if (DIPtr != NULL) bsuccess = bsuccess && STATUS_SUCCESS;
  else               bsuccess = bsuccess && STATUS_FAILURE;

  //--------------------------------------------------------------
  Device::InstanceBlock IB7;

  IB7.name = "L1";
  IB7.getModelName()   = "";
  tagged_param TP7("L",4.0);
  IB7.params.push_back(TP7);

  DIPtr = 0x0;
  DIPtr = DevMgrPtr_->addDeviceInstance(IB7);

  if (DIPtr != NULL) bsuccess = bsuccess && STATUS_SUCCESS;
  else               bsuccess = bsuccess && STATUS_FAILURE;
  //--------------------------------------------------------------

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::outputMI
// Purpose       : This function outputs all the model and instance lists
//                 that currently are set up in the device manager.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/20/00
//-----------------------------------------------------------------------------
bool DeviceTestor::outputMI()
{
  cout << endl;
  DevMgrPtr_->printOutLists();

  return 1;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::getPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/20/00
//-----------------------------------------------------------------------------
bool DeviceTestor::getPointers()
{
  int  i,imax;
  imax = _NUMDEV;

  string devname;

  Device::Device * DevicePtr;

  cout << "Testing the pointer return functionality" << endl;
  cout << "Number of Device types = " << imax << endl;

  for (i=0;i<imax;i++)
  {
      DevicePtr = DevMgrPtr_->returnDevicePtr(i);
      devname = DevicePtr->deviceName();

      cout << "Device # " << i << "  ";
      cout << devname << endl;
  }

  return STATUS_SUCCESS;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::getTopologies
// Purpose       : tests out the "internal" topology extraction functions.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/12/00
//-----------------------------------------------------------------------------
bool DeviceTestor::getTopologies()
{
  int  i,imax;
  imax = _NUMDEV;
  int inum;

  cout << endl << "Testing local topology functionality";
  cout << endl;
  cout << endl;

  cout << "First testing the internal vars functions:\n";
  cout << endl;

  string devname;

  Device::Device * DevicePtr;

  for (i=0;i<imax;i++)
  {
      DevicePtr = DevMgrPtr_->returnDevicePtr(i);
      devname = DevicePtr->deviceName();
      inum = DevicePtr->getNumIntVars ();

      cout << "Device # " << i << "  ";
      cout << devname << "Internal Vars: " << inum << endl;
  }


  cout << endl;
  cout << "Now testing the PrintOutTopologies functions";
  cout << endl;

  for (i=0;i<imax;i++)
  {
      DevicePtr = DevMgrPtr_->returnDevicePtr(i);
      DevicePtr->printOutTopologies ();
  }

  cout << endl;
  cout << "Now testing the ReturnLocalTopologies functions";
  cout << endl;


  return STATUS_SUCCESS;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::performRHSLoads
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/20/00
//-----------------------------------------------------------------------------
bool DeviceTestor::performRHSLoads()
{

  cout << "About to test out the RHS loader functions:" << endl;

  bool isuccess = STATUS_SUCCESS;

  isuccess = DevMgrPtr_->loadRHSVector ();

  return isuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::performJacLoads
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/20/00
//-----------------------------------------------------------------------------
bool DeviceTestor::performJacLoads()
{
  cout << "About to test out the Jacoabian loader functions:" << endl;

  bool isuccess = STATUS_SUCCESS;
  isuccess = DevMgrPtr_->loadJacobianMatrix ();

  return isuccess;
}


//-----------------------------------------------------------------------------
// Function      : DeviceTestor::getValue
// Purpose       : This function is the functional equivalent to the function
//                 PSN_Value in the class  ParseSpiceNetlist  which can be
//                 found in the DMFInterface package.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/27/00
//-----------------------------------------------------------------------------
REAL DeviceTestor::getValue (const string & str) const
{
  REAL value = atof(str.c_str());

  int j = str.find_first_not_of("0123456789.-+eE", 0);

  if (j == str.npos) return value;

  switch (str[j]) {
  case 'T' : case 't' :
      return value*1.0e12;
  case 'G' : case 'g' :
      return value*1.0e9;
  case 'M' :
      return value*1.0e6;
  case 'K' : case 'k' :
      return value*1.0e3;
  case 'm' :
      if (str[j+1] == 'i' && str[j+2] == 'l')
        return value*25.4e-6;
      else
        return value*1.0e-3;
  case 'u' : case 'U' :
      return value*1.0e-6;
  case 'n' : case 'N' :
      return value*1.0e-9;
  case 'p' : case 'P' :
      return value*1.0e-12;
  case 'f' : case 'F' :
      return value*1.0e-15;
  }
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::loadElementFields
// Purpose       : sets up the fields array for a bunch of different netlist
//                 specific devices.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/28/00
//-----------------------------------------------------------------------------
void DeviceTestor::loadElementFields(vector<string> & str, const int iElement)
{

  str.erase(str.begin(), str.end());

  switch (iElement)
  {
    case _RESISTOR_ELEMENT_TEST:
      str.reserve(4);
      str.push_back("R1");
      str.push_back("N1");
      str.push_back("N2");
      str.push_back("R");
      str.push_back("=");
      str.push_back("15k");
      break;
    case _CAPACITOR_ELEMENT_TEST:
      str.reserve(4);
      str.push_back("C1");
      str.push_back("N1");
      str.push_back("N2");
      str.push_back("C");
      str.push_back("=");
      str.push_back("15u");
      break;
    case _INDUCTOR_ELEMENT_TEST:
      str.reserve(4);
      str.push_back("L1");
      str.push_back("N1");
      str.push_back("N2");
      str.push_back("L");
      str.push_back("=");
      str.push_back("55g");
      break;
    case _VSRC_ELEMENT_TEST:
      str.reserve(13);
      str.push_back("V1");
      str.push_back("N1");
      str.push_back("N2");
      str.push_back("pulse");
      str.push_back("(");
      str.push_back("V1");
      str.push_back("=");
      str.push_back("0.0");    // V1
      str.push_back("V2");
      str.push_back("=");
      str.push_back("10.0");   // V2
      str.push_back("TD");
      str.push_back("=");
      str.push_back("5.0u");   // TD
      str.push_back("TR");
      str.push_back("=");
      str.push_back("5.0u");   // TR
      //str.push_back("5.0u");   // TF
      //str.push_back("20.0u");  // PW
      //str.push_back("20.0u");  // PER
      str.push_back(")");
      break;

    case _ISRC_ELEMENT_TEST:  // for this one, parameters are out-of-order.
      str.reserve(13);
      str.push_back("I1");
      str.push_back("N1");
      str.push_back("N2");
      str.push_back("sin");
      str.push_back("(");

      str.push_back("VO");
      str.push_back("=");
      str.push_back("0.0");

      str.push_back("THETA");
      str.push_back("=");
      str.push_back("5.0");

      str.push_back("FREQ");
      str.push_back("=");
      str.push_back("5.0u");

      str.push_back("TD");
      str.push_back("=");
      str.push_back("5.0u");

      str.push_back("VA");
      str.push_back("=");
      str.push_back("10.0");

      str.push_back(")");
      break;

    default:
      break;
  }

}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::getDefaultElementInfo
// Purpose       : The purpose  of this function is  to test out the
//                 getDefaultInstanceParams, etc. functions of the
//                 device manager.  These functions are used by the
//                 parser to obtain default tagged paramter lists to
//                 use in parsing the netlist.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/28/00
//-----------------------------------------------------------------------------
bool DeviceTestor::getDefaultElementInfo ()
{

  cout << "\nAbout to test getDefaultElementInfo\n\n";

  Device::InstanceBlock IB;
  list<string> nodelist;
  int isuccess = STATUS_SUCCESS;
  string left_paren("(");
  string right_paren(")");


  for  (int iElement = 0; iElement < _NUM_ELEMENT_TESTS; iElement++)
  {
    // first set up the vector of elements that would come in from the
    // netlist.  This assumes that an element line has already been
    // read in from the netlist and separated into fields.
    vector<string> field;

    loadElementFields(field, iElement);

    IB.clear();

    // Now process the various fields.
    // first get the name:
    IB.name  = field[0];

    // Now get the default paramter values from the device manager:
    DevMgrPtr_->getDefaultInstanceParams(IB);

    // now get the names of the nodes:
    int i,j;
    for (i=1; i<IB.iNumNodes+1; i++)
    {
      nodelist.push_back(field[i]);
    }

    // Determine if this is a "special" element - one requiring a model
    // statement or a source data block.  The type of device will already
    // have determined if these fields are required, and this information
    // is represented in the instance block by the modelFlag and sourceFlag
    // boolean variables.
    //
    // Note:  Some devices can have optional model statements.  Functionality
    // for optional model statements will be added later.

    if (IB.modelFlag) { IB.getModelName() = field[i]; i++; }

    // if this is a device which can use an "off" parameter", serach
    // through the various fields to find it.

    int iNumParams;
    iNumParams = field.size();

    if (IB.offFlag)
    {
       for (int itmp=0;itmp<iNumParams;itmp++)
       {
         ExtendedString tmpES = field[itmp];
         tmpES.toUpper();
         if  (tmpES == "OFF") IB.off = 1;
       }
    }

    // Now get the tagged paramters.  If this is a source, need some extra
    // info first.  In either case, need to check if the paramters have
    // been tagged in the netlist first.

    list<tagged_param>::iterator iter;
    list<tagged_param>::iterator begin;
    list<tagged_param>::iterator end;

    int istart, ifinish;

    if (IB.sourceFlag)
    {
      DevMgrPtr_->getDefaultSourceParams(field[i], IB);
      i++;
      iNumParams = field.size();

      // find the open and close parens:
      int ileft  = 0;
      int iright = 0;

      for (j=0;j<iNumParams;j++)
      {
         if (field[j] == "(") ileft = j;
         if (field[j] == ")") iright = j;
      }

      istart  = ileft+1;
      ifinish = iright;
    }
    else
    {
      istart = i;
      ifinish = field.size();
    }

    // Determine if these are tagged parameters or not.
    // Note that the parser will have already split out "=" signs into
    // separate fields.  Therefore the expression VO=5.0 will be 3 fields:
    //  field[i] = "VO"  field[i+1] = "="  and field[i+2] = "5.0"
    int itag = 0;
    for (i=istart;i<ifinish;i++) if(field[i] == "=") itag = 1;

    // Loop over the parameter fields, and load them into the instance
    // block.  For non-tagged fields assume that the paramters are in
    // the same order as the default tagged param list.  Any missing
    // non-tagged parameters are merely missing from the end of the list.

    begin = IB.params.begin();
    end   = IB.params.end();

    if (itag == 0)
    {
      if (istart < ifinish)
      {
        iter = begin;
        j    = istart;
        while (j < ifinish && iter != end)
        {
          iter->param = getValue(field[j]);
          j++; iter++;
        }
      }
    }
    else
    {
       for (iter=begin;iter != end; iter++)
       {
          for (i=istart;i<ifinish;i++)
          {
             ExtendedString tmpstring1 = field[i];
             ExtendedString tmpstring2 = iter->tag;
             tmpstring1.toUpper();
             tmpstring2.toUpper();
             if (tmpstring1 == tmpstring2) iter->param = getValue(field[i+2]);
          }
       }
    }

    cout << IB << endl;

  } // end of iElement loop.

  return isuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::loadModelFields
// Purpose       : sets up the fields array for a bunch of different netlist
//                 specific devices.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/28/00
//-----------------------------------------------------------------------------
void DeviceTestor::loadModelFields(vector<string> & str, const int iModel)
{

  str.erase(str.begin(), str.end());

  switch (iModel)
  {
    case _DIODE_MODEL_TEST:
      str.reserve(4);
      str.push_back("DIODE5");
      str.push_back("D");
      str.push_back("(");

      // junction area:
      str.push_back("AREA");
      str.push_back("=");
      str.push_back("10");

      // breakdown voltage
      str.push_back("BV");
      str.push_back("=");
      str.push_back("5");

      // explosion model parameter:
      str.push_back("EXPLI");
      str.push_back("=");
      str.push_back("10");

      // current at breakdown voltage: (?)
      str.push_back("IBV");
      str.push_back("=");
      str.push_back("5");

      // Forward knee current
      str.push_back("IK");
      str.push_back("=");
      str.push_back("5");

      // Reverse knee current
      str.push_back("IKR");
      str.push_back("=");
      str.push_back("5");

      // saturation current ?
      str.push_back("IS");
      str.push_back("=");
      str.push_back("5");

      // sidewall saturation current
      str.push_back("JSW");
      str.push_back("=");
      str.push_back("5");

      // length of the diode
      str.push_back("L");
      str.push_back("=");
      str.push_back("5");

      // level (1,2 or 3)
      str.push_back("LEVEL");
      str.push_back("=");
      str.push_back("3");

      // emmision coefficient
      str.push_back("N");
      str.push_back("=");
      str.push_back("5");

      // junction periphery
      str.push_back("PJ");
      str.push_back("=");
      str.push_back("5");

      // ohmic series resistance
      str.push_back("RS");
      str.push_back("=");
      str.push_back("5");

      // shrink factor
      str.push_back("SHRINK");
      str.push_back("=");
      str.push_back("5");

      // masking and/or etching effects:
      str.push_back("XW");
      str.push_back("=");
      str.push_back("5");


      // Junction capacitance parameters:




      str.push_back(")");

      break;
    default:
      break;
  }

}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::getDefaultModelInfo
// Purpose       : The purpose  of this function is
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/05/00
//-----------------------------------------------------------------------------
bool DeviceTestor::getDefaultModelInfo ()
{

  cout << "\nAbout to test getDefaultModelInfo\n\n";

  int i;
  Device::ModelBlock MB;

  list<string> nodelist;
  int isuccess = STATUS_SUCCESS;

  for  (int iModel = 0; iModel < _NUM_MODEL_TESTS; iModel++)
  {
    // first set up the vector of elements that would come in from the
    // netlist.  This assumes that an element line has already been
    // read in from the netlist and separated into fields.
    vector<string> field;

    loadModelFields(field, iModel);

    MB.clear();

    // Now process the various fields.
    // first get the name:
    i=0;
    MB.name  = field[i]; i++;
    MB.type  = field[i]; i++;

    // Note: parameter lists for model statements are *always* tagged!
    //
    // Check if the level has been specified.  If so, set the model block
    // level to that value.  If not, set it to 1.  It is neccessary to
    // do this now because the default parameters will dependend on the
    // level.
    list<tagged_param>::iterator iter;
    list<tagged_param>::iterator begin;
    list<tagged_param>::iterator end;

    int istart    = i;
    int ifinish   = field.size();

    MB.level = 1;
    for (i=istart;i<ifinish;i++)
    {
      ExtendedString fieldES = field[i];
      fieldES.toUpper();
      if(fieldES == "LEVEL") { MB.level = (int)(getValue(field[i+2])); }
    }

    // Now load the model block with the default paramter values from
    // the device manager:
    DevMgrPtr_->getDefaultModelParams(MB);

    // Loop over the parameter fields, and load them into the model block.
    begin = MB.params.begin();
    end   = MB.params.end();

    for (iter=begin;iter != end; iter++)
    {
       for (i=istart;i<ifinish;i++)
       {
          ExtendedString tmpstring1 = field[i];
          ExtendedString tmpstring2 = iter->tag;
          tmpstring1.toUpper();
          tmpstring2.toUpper();
          if (tmpstring1 == tmpstring2) iter->param = getValue(field[i+2]);
       }
    }

    cout << MB << endl;

  } // end of iElement loop.

  return isuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTestor::runTests
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/20/00
//-----------------------------------------------------------------------------
int DeviceTestor::runTests(int iargs, char *cargs[])
{
  int isuccess;

  cout << endl;
  cout << "Welcome to the DeviceModel Testing program." << endl;
  cout << endl;

  isuccess = setupParMgr ();

  //////////////////////////////////////////////////////////////
  // Allocate all the various packages:

  isuccess = doAllocations();

  if (isuccess == STATUS_SUCCESS)
  {
    Xyce::dout() << "Allocation was successful" << std::endl;
  }
  else
  {
    Xyce::Report::DevelFatal0() << "Allocation was NOT successful." << std::endl;
  }


  //////////////////////////////////////////////////////////////
  // Register the external package pointers:

  isuccess = doRegistrations();

  if (isuccess == STATUS_SUCCESS)
  {
    Xyce::dout() << "Registration was successful." << std::endl;
  }
  else
  {
    Xyce::Report::DevelFatal0() << "Registration was NOT successful." << std::endl;
  }

  isuccess = doInitializations ();

  //////////////////////////////////////////////////////////////
  // Test the  "dummy" pointers:
  isuccess = getPointers();

  //////////////////////////////////////////////////////////////
  // create one instance of each device type:
  isuccess = createDevices();

  if (isuccess == STATUS_SUCCESS)
  {
    Xyce::dout() << "Device creation was successful." << std::endl;
  }
  else
  {
    Xyce::Report::DevelFatal0() << "Device creation was NOT successful." << std::endl;
  }

  //////////////////////////////////////////////////////////////
  // Test the default parameter utilities, which would be used
  // by the netlist parser.
  isuccess = getDefaultElementInfo();

  if (isuccess == STATUS_SUCCESS)
  {
    Xyce::dout() << "getting default element info successful." << std::endl;
  }
  else
  {
    Xyce::Report::DevelFatal0() << "getting default element info NOT successful." << std::endl;
  }

  isuccess = getDefaultModelInfo();

  if (isuccess == STATUS_SUCCESS)
  {
    Xyce::dout() << "getting default model info successful." << std::endl;
  }
  else
  {
    Xyce::Report::DevelFatal0() << "getting default model info NOT successful." << std::endl;
  }

  //////////////////////////////////////////////////////////////
  // Test out if the device classes are singletons:
  isuccess = createDevices();

  if (isuccess == STATUS_SUCCESS)
  {
    Xyce::Report::DevelFatal0() << "Device creator does not enforce singletons." << std::endl;
  }
  else
  {
    Xyce::dout() << "Device creator does enforce singletons." << std::endl;
  }

  //////////////////////////////////////////////////////////////
  // Get device pointers of each type:
  isuccess = getPointers();

  //////////////////////////////////////////////////////////////
  // Flush the device array, test:
  //cout << "Right  before  deleting devices:"<< endl;
  //isuccess = deleteDevices();
  //isuccess = getPointers();
  //
  //cout << "Right  before  creating devices, again:"<< endl;
  //isuccess = createDevices();
  //isuccess = getPointers();
  //

  //////////////////////////////////////////////////////////////
  // Add models of each device type:
  isuccess = addModels();
  if (isuccess == STATUS_SUCCESS)
  {
    Xyce::dout() << "addModels was successful." << std::endl;
  }
  else
  {
    Xyce::Report::DevelFatal0() << "addModels was NOT successful." << std::endl;
  }


  //////////////////////////////////////////////////////////////
  // Add instances of each device type:
  isuccess = addInstances();

  if (isuccess == STATUS_SUCCESS)
  {
    Xyce::dout() << "addInstances was successful." << std::endl;
  }
  else
  {
    Xyce::Report::DevelFatal0() << "addInstances was NOT successful." << std::endl;
  }

  //////////////////////////////////////////////////////////////
  // Print out all the instances and models
  isuccess = outputMI();

  //////////////////////////////////////////////////////////////
  // Test Topological extraction:
  isuccess = getTopologies();

  //////////////////////////////////////////////////////////////
  // Test loader functions:
  isuccess = performRHSLoads();

  isuccess = performJacLoads();

  //////////////////////////////////////////////////////////////
  // Test State extraction:

  //////////////////////////////////////////////////////////////
  // Test ML creation:

  //////////////////////////////////////////////////////////////
  // de-allocate everything:

  isuccess = doDeAllocations();

  if (isuccess == STATUS_SUCCESS)
  {
    Xyce::dout() << "De-Allocation was successful." << std::endl;
  }
  else
  {
    Xyce::Report::DevelFatal0() << "De-Allocation was NOT successful." << std::endl;
  }

  //////////////////////////////////////////////////////////////
  return STATUS_SUCCESS;

}

//-----------------------------------------------------------------------------
// Function      : main
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/20/00
//-----------------------------------------------------------------------------
int main (int iargs, char *cargs[])
{
  DeviceTestor * DevTestPtr = new DeviceTestor();

  DevTestPtr->runTests(iargs, cargs);

  delete DevTestPtr;

  exit(0);
}


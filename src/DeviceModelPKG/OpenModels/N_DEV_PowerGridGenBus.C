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

//----------------------------------------------------------------------------
//
// Purpose        : PowerGrid classes: provides a device that calculates the
//                  steady state power flow from an Ideal Generator Bus.
//                  This device has a fixed output power (P) and a fixed 
//                  voltage magnitude (VM).
//
// Special Notes  : Experimental new device for an LDRD.
//
// Creator        : Pete Sholander, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 12/9/14
//
//
//
//
//----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_PowerGridGenBus.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Message.h>
#include <N_DEV_SolverState.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_IndentStreamBuf.h>

namespace Xyce {
namespace Device {
namespace PowerGridGenBus {

// enables debug output for just this device
//static const int DEBUG_DEVICE = 1;

void Traits::loadInstanceParameters(ParametricData<PowerGridGenBus::Instance> &p)
{
  p.addPar("AT", std::string("PQP"), &PowerGridGenBus::Instance::analysisTypeStr_)
    .setUnit(U_NONE)
    .setDescription("Analysis Type");

  p.addPar("P", 1.0, &PowerGridGenBus::Instance::power_)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_PERUNIT)
    .setDescription("Generator Output Power");

  p.addPar("VM", 1.0, &PowerGridGenBus::Instance::VMag_)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_PERUNIT)
    .setDescription("Voltage Magnitude");

  p.addPar("QMAX", 1.0, &PowerGridGenBus::Instance::QMax_)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_PERUNIT)
    .setDescription("Reactive Power Max Limit");

  p.addPar("QMIN", 0.0, &PowerGridGenBus::Instance::QMin_)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_PERUNIT)
    .setDescription("Reactive Power Min Limit");

  p.addPar("QLED", 0 , &PowerGridGenBus::Instance::QLimitEnforceDelay_)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_NONE)
    .setDescription("Q-Limit Enforcement Delay");
}

void Traits::loadModelParameters(ParametricData<PowerGridGenBus::Model> &p)
{}



std::vector< std::vector<int> > Instance::jacStamp;

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/9/14
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &               IB,
  Model & Miter,
  const FactoryBlock &  factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    analysisTypeStr_("PQP"),
    power_(1.0),
    VMag_(1.0),
    QMax_(1.0),
    QMin_(0.0),
    analysisType_(PQP),
    holdQatLower_(false),
    holdQatUpper_(false),
    QMaxGiven_(false),
    QMinGiven_(false),
    QLimitEnforceDelay_(0)
{
  // used by the Topology Checker.  See SON Bug 974 for more details.
  devConMap.resize(4);
  devConMap[0]=1;
  devConMap[1]=1;
  devConMap[2]=1;
  devConMap[3]=1;

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  processParams();
 
  // Finish setting up the device.  Process the atStr enum here, because
  // the calculations in processParams do not depend on the analysis type.
  ExtendedString atStr(analysisTypeStr_);
  atStr.toUpper();
  if ( atStr == "IV" )
  {
    // temporary error message since only PQP format works for the Generator Bus
    // model in Xyce 6.6
    UserError(*this) << "Only PQP Analysis Type is supported for PowerGridGenBus device.";
    analysisType_ = IV;
    analysisTypeStr_ = "IV";
  }
  else if ( atStr == "PQR" )
  {
    // temporary error message since only PQP format works for the Generator Bus
    // model in Xyce 6.6
    UserError(*this) << "Only PQP Analysis Type is supported for PowerGridGenBus device.";
    analysisType_ = PQR;
    analysisTypeStr_ = "PQR";
  }
  else if ( atStr == "PQP" )
  {
    analysisType_ = PQP;
    analysisTypeStr_ = "PQP";
  }
  // these are for test purposes for development, and were made "not reachable" for the
  // Xyce 6.6 release code.  
  //else if ( iter->stringValue() == "DVS" )
  //{
  //  analysisType_ = DVS;
  //  analysisTypeStr_ = "DVS";
  //}
  //else if ( iter->stringValue() == "VCS" )
  //{
  //  analysisType_ = VCS;
  //   analysisTypeStr_ = "VCS";
  //}
  //else if ( iter->stringValue() == "CCS" )
  //{
  //  analysisType_ = CCS;
  //  analysisTypeStr_ = "CCS";
  //}
  //else if ( iter->stringValue() == "RCB" )
  //{
  //  analysisType_ = RCB;
  //  analysisTypeStr_ = "RCB";
  //}
  else
  {
    UserError(*this) << "Analysis Type must be IV, PQR or PQP in power grid device.";
  }

  if (given("QMAX"))
  {
    QMaxGiven_ = true;
  }
  if (given("QMIN"))
  {
    QMinGiven_ = true;
  }

  if ( analysisType_ == DVS || analysisType_ == IV ||  analysisType_ == PQR ||
       analysisType_ == VCS || analysisType_ == RCB)
  {
    numIntVars   = 2;
    numExtVars   = 4;   
    numStateVars = 0;
  }
  else if (analysisType_ == PQP || analysisType_ == CCS) 
  {
    numIntVars   = 1;
    numExtVars   = 4;   
    numStateVars = 0;
  }
  else
  {
    UserError(*this) << "Analysis Type must be IV, PQR or PQP in power grid device.";
    return;
  }
 
  if( jacStamp.empty() )
  {
    if ( analysisType_ == DVS || analysisType_ == IV )
    {
      jacStamp.resize(6);
      jacStamp[0].resize(1);
      jacStamp[1].resize(1);
      jacStamp[2].resize(1);
      jacStamp[3].resize(1);
      jacStamp[4].resize(2);
      jacStamp[5].resize(2);

      jacStamp[0][0] = 4;  
      jacStamp[1][0] = 4;
      jacStamp[2][0] = 5;
      jacStamp[3][0] = 5;
      jacStamp[4][0] = 0;
      jacStamp[4][1] = 1;
      jacStamp[5][0] = 2;
      jacStamp[5][1] = 3;
    }
    else if (analysisType_ == VCS )
    { // this is for test purposes.  It is a voltage-current source, where the branch current
      // is a solution variable.  The VCS code is not reachable in the Xyce 6.4 release code.
      jacStamp.resize(6);
      jacStamp[0].resize(1);
      jacStamp[1].resize(1);
      jacStamp[2].resize(1);
      jacStamp[3].resize(1);
      jacStamp[4].resize(2);
      jacStamp[5].resize(1);

      jacStamp[0][0] = 4;  
      jacStamp[1][0] = 4;
      jacStamp[2][0] = 5;
      jacStamp[3][0] = 5;
      jacStamp[4][0] = 0;
      jacStamp[4][1] = 1;
      jacStamp[5][0] = 5;
    }
    else if (analysisType_ == PQR) 
    {
      jacStamp.resize(6);
      jacStamp[0].resize(1);
      jacStamp[1].resize(1);
      jacStamp[2].resize(3);
      jacStamp[3].resize(3);
      jacStamp[4].resize(2);
      jacStamp[5].resize(2);

      jacStamp[0][0] = 4;  
      jacStamp[1][0] = 4;
      jacStamp[2][0] = 0;
      jacStamp[2][1] = 1;
      jacStamp[2][2] = 5;
      jacStamp[3][0] = 0;
      jacStamp[3][1] = 1;
      jacStamp[3][2] = 5;
      jacStamp[4][0] = 0;
      jacStamp[4][1] = 1;
      jacStamp[5][0] = 2;
      jacStamp[5][1] = 3;
    }
    else if (analysisType_ == PQP ) 
    {
      jacStamp.resize(5);
      jacStamp[2].resize(1);
      jacStamp[3].resize(1);
      jacStamp[4].resize(3);

      jacStamp[2][0] = 4;
      jacStamp[3][0] = 4;
      jacStamp[4][0] = 2;
      jacStamp[4][1] = 3;
      jacStamp[4][2] = 4;
    }
    else if (analysisType_ == CCS ) 
    { // This is for test purposes.  It is a current-current source, where the branch current is a
      // solution variable.  The CCS code is not reachable in the Xyce 6.4 release code.
      jacStamp.resize(5);
      jacStamp[2].resize(1);
      jacStamp[3].resize(1);
      jacStamp[4].resize(1);

      jacStamp[2][0] = 4;
      jacStamp[3][0] = 4;
      jacStamp[4][0] = 4;
    }
    else if (analysisType_ == RCB )
    {
      //This is for test purposes.  It is a prototype of a remote controlled bus.  This device
      // is not functional yet.  The RCB code is not reachable in the Xyce 6.4 release code.
      jacStamp.resize(6);
      jacStamp[0].resize(1);
      jacStamp[1].resize(1);
      jacStamp[4].resize(2);
      jacStamp[5].resize(2);

      jacStamp[0][0] = 4;
      jacStamp[1][0] = 4;
      jacStamp[4][0] = 0;
      jacStamp[4][1] = 1;
      jacStamp[5][0] = 2;
      jacStamp[5][1] = 3;
    }
  }  
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/9/14
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 11/17/16
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  // Error checking for a bad instance line that will cause problems later
  if ( VMag_ <= 0.0 )
  {
    UserError(*this) << "VM must be positive for PowerGridGenBus bus device.";
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef )
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // These if-else blocks may be made into separate classes,
  // in the final version of this model. 
  if ( analysisType_ == DVS || analysisType_ == VCS || analysisType_ == IV)
  {
    li_VR1 = extLIDVec[0];
    li_VR2 = extLIDVec[1];
    li_VI1 = extLIDVec[2];
    li_VI2 = extLIDVec[3];
    li_BranCurrR = intLIDVec[0];
    li_BranCurrI = intLIDVec[1];
    
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << getName() << " LIDs"
        << Util::push << std::endl
           << "li_VR1 = " << li_VR1 << std::endl
           << "li_VI1 = " << li_VI1 << std::endl
           << "li_VR2 = " << li_VR2 << std::endl
           << "li_VI2 = " << li_VI2 << std::endl
           << "li_BranCurrR = " << li_BranCurrR << std::endl
           << "li_BranCurrI = " << li_BranCurrI << std::endl
           << std::endl;
    }
  }
  else if ( analysisType_ == PQR)
  {
    li_VR1 = extLIDVec[0];
    li_VR2 = extLIDVec[1];
    li_VI1 = extLIDVec[2];
    li_VI2 = extLIDVec[3];
    li_BranCurrP = intLIDVec[0];
    li_BranCurrQ = intLIDVec[1];
    
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << getName() << " LIDs"
        << Util::push << std::endl
           << "li_VR1 = " << li_VR1 << std::endl
           << "li_VI1 = " << li_VI1 << std::endl
           << "li_VR2 = " << li_VR2 << std::endl
           << "li_VI2 = " << li_VI2 << std::endl
           << "li_BranCurrP = " << li_BranCurrP << std::endl
           << "li_BranCurrQ = " << li_BranCurrQ << std::endl
           << std::endl;
    }
  }
  else if (analysisType_ == PQP || analysisType_ == CCS) 
  {
    li_Theta1 = extLIDVec[0];
    li_Theta2 = extLIDVec[1];
    li_VM1 = extLIDVec[2];
    li_VM2 = extLIDVec[3];
    li_BranCurrQ = intLIDVec[0];
  
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << getName() << " LIDs"
        << Util::push << std::endl
           << "li_Theta1 = " << li_Theta1 << std::endl
           << "li_VM1 = " << li_VM1 << std::endl
           << "li_Theta2 = " << li_Theta2 << std::endl
           << "li_VM2 = " << li_VM2 << std::endl
           << "li_BranCurrQ = " << li_BranCurrQ << std::endl
           << std::endl;
    }
  }
  else if (analysisType_ == RCB)
  {
    li_Theta1 = extLIDVec[0];
    li_Theta2 = extLIDVec[1];
    li_VM1 = extLIDVec[2];
    li_VM2 = extLIDVec[3];
    li_BranCurrP = intLIDVec[0];
    li_BranCurrQ = intLIDVec[1];
  }
  else
  {
    UserError(*this) << "Analysis Type must be IV, PQR or PQP in power grid device.";
    return;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadNodeSymbols
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
void Instance::loadNodeSymbols(Util::SymbolTable &symbol_table) const
{
  // These if-else blocks may be made into separate classes,
  // in the final version of this model
  if ( analysisType_ == DVS || analysisType_ == VCS || analysisType_ == IV )
  {
    addInternalNode(symbol_table, li_BranCurrR, getName(), "BranchCurrR");
    addInternalNode(symbol_table, li_BranCurrI, getName(), "BranchCurrI");
  }
  else if ( analysisType_ == PQR || analysisType_ == RCB )
  {
    addInternalNode(symbol_table, li_BranCurrP, getName(), "BranchCurrP");
    addInternalNode(symbol_table, li_BranCurrQ, getName(), "BranchCurrQ");
  }
  else if ( analysisType_ == PQP || analysisType_ == CCS)
  {
    addInternalNode(symbol_table, li_BranCurrQ, getName(), "BranchCurrQ");
  }
  else
  {
    UserError(*this) << "Analysis Type must be IV, PQR or PQP in power grid device: " << getName();
    return;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       : Set up offsets so we can store quantities that need to be
//                 differentiated.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs(const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/20/01
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/12/14
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  // // These if-else blocks may be made into separate classes,
  // in the final version of this model
  if (analysisType_ == DVS || analysisType_ == IV )
  {
    VR1_IR_Offset = jacLIDVec[0][0];
    VR2_IR_Offset = jacLIDVec[1][0];
    VI1_II_Offset = jacLIDVec[2][0];
    VI2_II_Offset = jacLIDVec[3][0];

    IR_VR1_Offset = jacLIDVec[4][0];
    IR_VR2_Offset = jacLIDVec[4][1];
    II_VI1_Offset = jacLIDVec[5][0];
    II_VI2_Offset = jacLIDVec[5][1];

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << getName() << ": In registerJacLIDs, the offsets are:" << std::endl
                   << "VR offsets are (" << VR1_IR_Offset << " , " << VR2_IR_Offset << ")" << std::endl
	           << "VI offsets are (" << VI1_II_Offset << " , " << VI1_II_Offset << ")" << std::endl
		   << "IR offsets are (" << IR_VR1_Offset << " , " << IR_VR2_Offset << ")" << std::endl
	           << "II offsets are (" << II_VI1_Offset << " , " << II_VI2_Offset << ")" << std::endl 
                   << std::endl;
    }
  }
  else if (analysisType_ == VCS)
  { // This is for test purposes.  It is a voltage-current source, where the branch current is
    // a solution variable for the current source.
    VR1_IR_Offset = jacLIDVec[0][0];
    VR2_IR_Offset = jacLIDVec[1][0];
    VI1_II_Offset = jacLIDVec[2][0];
    VI2_II_Offset = jacLIDVec[3][0];

    IR_VR1_Offset = jacLIDVec[4][0];
    IR_VR2_Offset = jacLIDVec[4][1];
    II_II_Offset = jacLIDVec[5][0];

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << getName() << ": In registerJacLIDs, the offsets are:" << std::endl
                   << "VR offsets are (" << VR1_IR_Offset << " , " << VR2_IR_Offset << ")" << std::endl
	           << "VI offsets are (" << VI1_II_Offset << " , " << VI1_II_Offset << ")" << std::endl
		   << "IR offsets are (" << IR_VR1_Offset << " , " << IR_VR2_Offset << ")" << std::endl
	           << "II offsets are (" << II_VI1_Offset << " , " << II_VI2_Offset << ")" << std::endl 
                   << std::endl;
    }
  }
  else if (analysisType_ == PQR)
  {
    VR1_P_Offset = jacLIDVec[0][0];  
    VR2_P_Offset = jacLIDVec[1][0];
    VI1_VR1_Offset = jacLIDVec[2][0];
    VI1_VR2_Offset = jacLIDVec[2][1];
    VI1_Q_Offset = jacLIDVec[2][2];
    VI2_VR1_Offset = jacLIDVec[3][0];
    VI2_VR2_Offset = jacLIDVec[3][1];
    VI2_Q_Offset = jacLIDVec[3][2];

    P_VR1_Offset = jacLIDVec[4][0];
    P_VR2_Offset = jacLIDVec[4][1];
    Q_VI1_Offset = jacLIDVec[5][0];
    Q_VI2_Offset = jacLIDVec[5][1];

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << getName() << ": In registerJacLIDs, the offsets are:" << std::endl
                   << "VR offsets are (" << VR1_P_Offset << " , " << VR2_P_Offset << ")" << std::endl
	           << "VI offsets are (" << VI1_Q_Offset << " , " << VI1_Q_Offset << ")" << std::endl
		   << "IR offsets are (" << P_VR1_Offset << " , " << P_VR2_Offset << ")" << std::endl
	           << "II offsets are (" << Q_VI1_Offset << " , " << Q_VI2_Offset << ")" << std::endl 
                   << std::endl; 
    }
  }
  else if (analysisType_ == PQP)
  {
    VM1_Q_Offset = jacLIDVec[2][0];
    VM2_Q_Offset = jacLIDVec[3][0];

    Q_VM1_Offset = jacLIDVec[4][0];
    Q_VM2_Offset = jacLIDVec[4][1];
    Q_Q_Offset = jacLIDVec[4][2];

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << getName() << ": In registerJacLIDs, the offsets are:" << std::endl
	           << "VM offsets are (" << VM1_Q_Offset << " , " << VM2_Q_Offset << ")" << std::endl
	           << "Q offsets are (" << Q_VM1_Offset << " , " << Q_VM2_Offset 
                   << " , " << Q_Q_Offset << ")" << std::endl << std::endl;     
    }
  }
  else if (analysisType_ == CCS)
  { // This is for test purposes.  It is a current-current source, where the branch current is 
    // a solution variable.
    VM1_Q_Offset = jacLIDVec[2][0];
    VM2_Q_Offset = jacLIDVec[3][0];

    Q_Q_Offset = jacLIDVec[4][0];

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << getName() << ": In registerJacLIDs, the offsets are:" << std::endl
	           << "VM offsets are (" << VM1_Q_Offset << " , " << VM2_Q_Offset << ")" << std::endl
	           << "Q offsets are (" << Q_Q_Offset << ")" << std::endl 
                   << std::endl;     
    }
  }
  else if (analysisType_ == RCB)
  {
     Th1_P_Offset = jacLIDVec[0][0];
     Th2_P_Offset = jacLIDVec[1][0];
     P_Th1_Offset = jacLIDVec[4][0];
     P_Th2_Offset = jacLIDVec[4][1];
     Q_VM1_Offset = jacLIDVec[5][0];
     Q_VM2_Offset = jacLIDVec[5][1];
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one power grid
//                 generator bus instance
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  double * solVec = extData.nextSolVectorRawPtr;

  // These if-else blocks may be made into separate classes,
  // in the final version of this model
  if (analysisType_ == DVS || analysisType_ == VCS || analysisType_ == IV)
  {
    VR1 = solVec[li_VR1];
    VR2 = solVec[li_VR2];
    VI1 = solVec[li_VI1];
    VI2 = solVec[li_VI2];

    Bran_CurrR = solVec[li_BranCurrR];
    Bran_CurrI = solVec[li_BranCurrI]; 

    srcDropR = VR1 - VR2;
    srcBCR = VMag_;
    srcVoltageR = srcDropR - srcBCR;

    srcDropI = VI1 - VI2;
    srcBCI = VMag_;
    srcVoltageI = srcDropI - srcBCI;
  }
  else if (analysisType_ == PQR)
  {
    // implement P as a current source, which means it has no state

    // VI will be implmented as voltage source that depends on VR
    VR1 = solVec[li_VR1];
    VR2 = solVec[li_VR2];
    VI1 = solVec[li_VI1];
    VI2 = solVec[li_VI2];

    Bran_CurrP = solVec[li_BranCurrR];
    Bran_CurrQ = solVec[li_BranCurrI];

    srcDropR = VR1 - VR2;
    srcDropI = VI1 - VI2;

    // enforce Voltage Magnitude constraint
    if (srcDropR > VMag_)
    {
      srcBCR = VMag_;
      srcBCI = 0;
    }
    else if (-srcDropR > VMag_)
    {
      srcBCR = -VMag_;
      srcBCI = 0;
    }
    else 
    {
      srcBCR = srcDropR;
      double signVal = (srcBCI < 0) ? -1.0 : 1.0;
      srcBCI = signVal * sqrt( VMag_*VMag_ - srcDropR * srcDropR);
    }
    
    srcVoltageR = srcDropR - srcBCR;
    srcVoltageI = srcDropI - srcBCI;

    /*if (Bran_CurrQ > QMax_)
    {
      reactPower_ = QMax_;
    }
    else if (Bran_CurrQ < QMin_)
    {
      reactPower_ = QMin_;
    }
    else
    {
      reactPower_ = Bran_CurrQ;
    }*/

    /*Q1 = Bran_CurrQ;
    if (Q1 < -QMax_)
    {
      Q1 = -QMax_;
    }
    else if (Q1 > -QMin_)
    {
      Q1 = -QMin_;
    }
    else
    { // Q is within limits.  So, enforce the voltage magnitude constraint

    }

    Q2 = -Q1;*/
  }
  else if (analysisType_ == PQP)
  {    
    // implement P as a current source, which means it has no state

    // implement VM as a voltage source
    VM1 = solVec[li_VM1];
    VM2 = solVec[li_VM2];
    Bran_CurrQ = solVec[li_BranCurrQ];
    if ( getSolverState().newtonIter > QLimitEnforceDelay_ )
    {
      if ( (Bran_CurrQ > QMax_) && QMaxGiven_ )
      { 
        holdQatUpper_ = true;
      }
      else if ( (Bran_CurrQ < QMin_) && QMinGiven_ )
      {
        holdQatLower_ = true; 
      }
    }

    //Xyce::dout() << getName() << ": Bran_CurrQ and NewtIter in UpdateIntermediate State is: " << Bran_CurrQ 
    //		 <<  " and " << getSolverState().newtonIter << std::endl;

    srcDropVM = VM1 - VM2;
    srcBCVM = VMag_;
    srcVoltageVM = srcDropVM - srcBCVM;
  }
  else if (analysisType_ == CCS)
  { // this is for test purposes.  It is a current-current source.    
    // implement P as a current source, which means it has no state

    // implement a current source with the branch current as a solution variable.
    Bran_CurrQ = solVec[li_BranCurrQ];
  }
  else if (analysisType_ == RCB)
  {
    // implement VM as a voltage source
    Th1 = solVec[li_Theta1];
    Th2 = solVec[li_Theta2];
    VM1 = solVec[li_VM1];
    VM2 = solVec[li_VM2];
    Bran_CurrP = solVec[li_BranCurrP];

    srcDropTh = Th1 - Th1;
    srcDropVM = VM1 - VM2;
    srcBCVM = VMag_;
  }
  else
  {
    UserError(*this) << "Analysis Type must be IV, PQR or PQP in power grid device: " << getName();
    return false;
  }
 
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  return updateIntermediateVars ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions 
//
// Special Notes :
//
//
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;
  
  // These if-else blocks may be made into separate classes,
  // in the final version of this model
  if (analysisType_ == DVS || analysisType_ == IV)
  {
    fVec[li_VR1] += Bran_CurrR;
    fVec[li_VR2] -= Bran_CurrR; 
    fVec[li_VI1] += Bran_CurrI; 
    fVec[li_VI2] -= Bran_CurrI;
    
    fVec[li_BranCurrR] += srcDropR;
    fVec[li_BranCurrI] += srcDropI; 
  }
    else if (analysisType_ == VCS )
  {
    fVec[li_VR1] += Bran_CurrR;
    fVec[li_VR2] -= Bran_CurrR; 
    fVec[li_VI1] -= Bran_CurrI; 
    fVec[li_VI2] += Bran_CurrI;
    
    fVec[li_BranCurrR] += srcDropR;
    fVec[li_BranCurrI] += Bran_CurrI; 
  }
  else if (analysisType_ == PQR)
  {
    fVec[li_VR1] -= power_;
    fVec[li_VR2] += power_; 
    fVec[li_VI1] -= Bran_CurrQ; 
    fVec[li_VI2] += Bran_CurrQ;  

    fVec[li_BranCurrP] -= srcDropR;
    fVec[li_BranCurrP] += srcBCR; 

    fVec[li_BranCurrQ] -= srcDropI;
    fVec[li_BranCurrQ] += srcBCI;

    Xyce::dout() << getName() << ": F Vector Load info is (P,Q): (" << power_ << ","
                 << Bran_CurrQ << ")" << std::endl
                                << "                       (srcDropR,srcBCR,srcDropI,srcBCI: (" 
				<< srcDropR << "," << srcBCR << "," <<srcDropI << "," 
                                << srcBCI << ")"<< std::endl << std::endl
		 << "                  Solution Vector  Currents are (P,Q): (" << Bran_CurrP << ","
		 << Bran_CurrQ << ")" << std::endl
                                << "                       (VR1,VR2,VI1,VI2): (" 
				<< VR1 << "," << VR2 << "," << VI1 << "," 
                                << VI2 << ")"<< std::endl << std::endl;    
    //fVec[li_VR1] += P1;
    //fVec[li_VR2] += P2; 
    //fVec[li_VI1] += Q1; 
    //fVec[li_VI2] += Q2;  
  }
  else if (analysisType_ == PQP)
  {
    // initial version does not have reactive power limiting.  
    
    // implement P as a current source, which means it has no F-Vector contribution.  
    // Note though that the P is a generator. So, the current polarity is reversed from 
    // a normal Xyce current source, as explained below. The polarity difference is
    // accounted for in the loadDAEBVector() function.
    

    // implement VM as a voltage source, when we're inside of the Q limits.  Implement a current source
    // for the VM terminals, otherwise.  This device is a generator.  So, the current polarities on the
    // branch current contributions, fVec[li_VM1] and  fVec[li_VM2], are reversed in both cases from a 
    // normal Xyce voltage source. Also note that the dFdx[li_VM1][VM1_Q_Offset] and 
    // dFdx[li_VM2][VM2_Q_Offset] terms in the dFdX equations in loadDAEdFdx () are also reversed
    // from a normal Xyce voltage source.  For a Xyce voltage source, "positive current" flows 
    // from the positive node through the source to the negative node.  For a Power Grid generator
    // (like this device), the convention is that positive reactive power (Q) flows into the power 
    // grid from the positive VM terminal.  A similar convention holds for the real power (P) flow
    // from the positive Theta terminal.
     if (holdQatLower_ || holdQatUpper_)
    {
      // Outside of the Q limits.  So, implement a current source with the branch current as a
      // solution variable, since the jacStamp has to remain the same size.   Reverse the signs though
      // on the voltage terms since it's a generator, and positive power flows out of the positive (VM1)
      // terminal.
      fVec[li_VM1] -= Bran_CurrQ;
      fVec[li_VM2] += Bran_CurrQ;
      fVec[li_BranCurrQ] += Bran_CurrQ; 

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
        Xyce::dout() << getName() << ": F Vector Load info is (Q,QMax,Qmin): (" << Bran_CurrQ << 
		     "," << QMax_ << "," << QMin_ << ")" << std::endl << std::endl;     
      }
    }
    else
    {
      // implement as a voltage source, since we're were within the Q limits.  Again, the signs on
      // VM terms are reversed since positive power flows out of the positive (VM1) terminal.
      fVec[li_VM1] -= Bran_CurrQ;
      fVec[li_VM2] += Bran_CurrQ; 
      if (DEBUG_DEVICE && isActive(Diag::DEVICE_LOAD_VECTOR))
      {
        Xyce::dout() << "Bran_CurrQ = " << Bran_CurrQ << std::endl;
      }

      fVec[li_BranCurrQ] += srcDropVM;

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
        Xyce::dout() << getName() << ": F Vector Load info is (P,Q,srcDropVM,srcBCVM): (" << power_ << " , "
                   << Bran_CurrQ << " , " << srcDropVM << " , " << srcBCVM << ")" << std::endl << std::endl;     
      }
    }
  }
  else if (analysisType_ == CCS)
  { // for test purposes
    fVec[li_VM1] -= Bran_CurrQ;
    fVec[li_VM2] += Bran_CurrQ;
    fVec[li_BranCurrQ] += Bran_CurrQ;
  }
  else if (analysisType_ == RCB)
  {
     fVec[li_Theta1] -= Bran_CurrP;
     fVec[li_Theta2] += Bran_CurrP; 
     fVec[li_BranCurrP] += srcDropTh;
     fVec[li_BranCurrQ] += srcDropVM;
  }
  else
  {
    UserError(*this) << "Analysis Type must be IV, PQR or PQP in power grid device: " << getName();
    return false;
  } 
 
  if (getSolverState().dcopFlag)
  {
    // no op placeholder
  }
  else
  {
    // no op placeholder
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEBVector
//
// Purpose       : Loads the B-vector contributions for a single
//                 power grid generator bus instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool Instance::loadDAEBVector ()
{
  double * bVec = extData.daeBVectorRawPtr;

  // These if-else blocks may be made into separate classes,
  // in the final version of this model
  if (analysisType_ == DVS || analysisType_ == IV )
  {
    bVec[li_BranCurrR] += srcBCR;
    bVec[li_BranCurrI] += srcBCI;
  }
  else if (analysisType_ == VCS)
  { // for test purposes.  A voltage-current source
    bVec[li_BranCurrR] += srcBCR;
    bVec[li_BranCurrI] += power_;
  }
  else if (analysisType_ == PQR)
  {
    bVec[li_BranCurrR] -= srcBCR;
    bVec[li_BranCurrI] -= srcBCI;
  }
  else if (analysisType_ == PQP)
  {
    // current source for P.   Note the P is a generator. So, the current polarity for 
    // the bVec[li_Theta1] and bVec[li_Theta2] contributions are reversed from a normal Xyce
    // current source.  The comments for the loadDAEFVector() function provide more explanation.
    bVec[li_Theta1] += power_;
    bVec[li_Theta2] -= power_;  

     // Implement a current source for Q, if we're outside of one of the Q limits.
    if (holdQatUpper_ )
    {
      bVec[li_BranCurrQ] += QMax_;
    }
    else if (holdQatLower_)
    {
      bVec[li_BranCurrQ] += QMin_;
    }        
    else   
    {
    // implement as a voltage source for VM since we're were within the Q limits
      bVec[li_BranCurrQ] += srcBCVM;
    }
  }
  else if (analysisType_ == CCS)
  { // for test purposes.  A current-current source
    bVec[li_Theta1] += power_;
    bVec[li_Theta2] -= power_;  
    bVec[li_BranCurrQ] += power_;
  }
  else if (analysisType_ == RCB)
  {
    bVec[li_BranCurrQ] += srcBCVM;
  }
  else
  {
    UserError(*this) << "Analysis Type must be IV, PQR or PQP in power grid device: " << getName();
    return false;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the dFdx-matrix contributions 
//
// Special Notes :
//
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  // These if-else blocks may be made into separate classes,
  // in the final version of this model
  if (analysisType_ == DVS || analysisType_ == IV)
  {
    dFdx[li_VR1][VR1_IR_Offset] += 1.0;
    dFdx[li_VR2][VR2_IR_Offset] -= 1.0;
    dFdx[li_VI1][VI1_II_Offset] += 1.0;
    dFdx[li_VI2][VI2_II_Offset] -= 1.0;

    dFdx[li_BranCurrR][IR_VR1_Offset] += 1.0;
    dFdx[li_BranCurrR][IR_VR2_Offset] -= 1.0;
    dFdx[li_BranCurrI][II_VI1_Offset] += 1.0;
    dFdx[li_BranCurrI][II_VI2_Offset] -= 1.0;
  }
    else if (analysisType_ == VCS)
  { // for test purposes.  A voltage-current source
    dFdx[li_VR1][VR1_IR_Offset] += 1.0;
    dFdx[li_VR2][VR2_IR_Offset] -= 1.0;
    dFdx[li_VI1][VI1_II_Offset] -= 1.0;
    dFdx[li_VI2][VI2_II_Offset] += 1.0;

    dFdx[li_BranCurrR][IR_VR1_Offset] += 1.0;
    dFdx[li_BranCurrR][IR_VR2_Offset] -= 1.0;
    dFdx[li_BranCurrI][II_II_Offset] += 1.0;
  }
  else if (analysisType_ == PQR)
  {
    dFdx[li_VR1][VR1_P_Offset] += 1.0;
    dFdx[li_VR2][VR2_P_Offset] -= 1.0;
    dFdx[li_VI1][VI1_Q_Offset] += 1.0;
    dFdx[li_VI2][VI2_Q_Offset] -= 1.0;

    // jacobiam term that come from the VMag constraint
    double jacTerm = (VMag_ - std::abs(VR1-VR2)) > 0 
                         ?  (VR2 - VR1) / sqrt(VMag_*VMag_ - (VR1-VR2)*(VR1-VR2)) : 1e10*(VR2 - VR1); 
    Xyce::dout() << "Jacobian term = " << jacTerm << std::endl; 
    dFdx[li_VI1][VI1_VR1_Offset] += jacTerm;
    dFdx[li_VI1][VI1_VR2_Offset] -= jacTerm;
    dFdx[li_VI2][VI2_VR1_Offset] -= jacTerm;
    dFdx[li_VI2][VI2_VR2_Offset] += jacTerm;

    //Jacobian terms from the voltage drop constraints
    dFdx[li_BranCurrP][P_VR1_Offset] -= 1.0;
    dFdx[li_BranCurrP][P_VR2_Offset] += 1.0;
    dFdx[li_BranCurrQ][Q_VI1_Offset] -= 1.0;
    dFdx[li_BranCurrQ][Q_VI2_Offset] += 1.0;           
  }
  else if (analysisType_ == PQP)
  {
    // P is implemented as a current source.  So, there are no dFdX terms for it

    // VM is implemented as a voltage source. This device is a generator.  So, the signs of the 
    // dFdx[li_VM1][VM1_Q_Offset] and dFdx[li_VM2][VM2_Q_Offset] terms are reversed from that of
    // a normal Xyce voltage source.  Also note that the signs of the branch current contributions,
    // fVec[li_VM1] and  fVec[li_VM2], in loadDAEFVector () are also reversed from a normal Xyce
    //  voltage source. 
     if (holdQatLower_ || holdQatUpper_)
    {
      // Outside of the Q limits.  So, this is a current-source with the branch current as a solution
      // variable. 
      dFdx[li_VM1][VM1_Q_Offset] -= 1.0;
      dFdx[li_VM2][VM2_Q_Offset] += 1.0;

      dFdx[li_BranCurrQ][Q_Q_Offset] += 1.0;
    } 
    else 
    {
      // implement as a voltage source, since we're were within the Q limits
      dFdx[li_VM1][VM1_Q_Offset] -= 1.0;
      dFdx[li_VM2][VM2_Q_Offset] += 1.0;

      dFdx[li_BranCurrQ][Q_VM1_Offset] += 1.0;
      dFdx[li_BranCurrQ][Q_VM2_Offset] -= 1.0;
    }
  }
  else if (analysisType_ == CCS)
  { // for test purposes.  A current-current souce
    dFdx[li_VM1][VM1_Q_Offset] -= 1.0;
    dFdx[li_VM2][VM2_Q_Offset] += 1.0;

    dFdx[li_BranCurrQ][Q_Q_Offset] += 1.0;
  }
  else if (analysisType_ == RCB)
  {
    dFdx[li_Theta1][Theta1_P_Offset] -= 1.0;
    dFdx[li_Theta2][Theta2_P_Offset] += 1.0; 
    
    dFdx[li_BranCurrP][P_Theta1_Offset] += 1.0;
    dFdx[li_BranCurrP][P_Theta2_Offset] -= 1.0;

    dFdx[li_BranCurrQ][Q_VM1_Offset] += 1.0;
    dFdx[li_BranCurrQ][Q_VM2_Offset] -= 1.0;
  }
  else
  {
    UserError(*this) << "Analysis Type must be IV, PQR or PQP in power grid device: " << getName();
    return false;
  }
   
  if (getSolverState().dcopFlag)
  {
    // no op placeholder code;
  }
  else
  {
    // no op placeholder code;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : "Model block" constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------


// These are all just placeholder functions, as there is nothing to the model
// class.
// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : "Model block" constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block)
{
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
Model::~Model ()
{
  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    delete (*iter);
  }
}

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << std::endl;
  os << "    name     model name  Parameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "      ";
    os << getName();
    os << std::endl;
  }

  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Model::forEachInstance
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date : 2/4/2014
//-----------------------------------------------------------------------------
/// Apply a device instance "op" to all instances associated with this
/// model
/// 
/// @param[in] op Operator to apply to all instances.
/// 
/// 
void Model::forEachInstance(DeviceInstanceOp &op) const /* override */ 
{
  for (std::vector<Instance *>::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
    op(*it);
}


Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new DeviceMaster<Traits>(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() || (deviceMap.find("POWERGRIDGENBUS")!=deviceMap.end())
                        || (deviceMap.find("PGGB")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("PowerGridGenBus", 1)
      .registerDevice("PGGB", 1)
      .registerModelType("PowerGridGenBus", 1);
  }
}

} // namespace PowerGridGenBus
} // namespace Device
} // namespace Xyce

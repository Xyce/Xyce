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

//----------------------------------------------------------------------------
//
// Purpose        : PowerGrid classes: provides a device that calculates the
//                  steady state power flow in a transmission grid branch
//
// Special Notes  : Experimental new device for an LDRD.
//
// Creator        : Pete Sholander, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 9/12/14
//
//
//
//
//----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_PowerGridBranch.h>
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
namespace PowerGridBranch {

// enables debug output for just this device
//static const int DEBUG_DEVICE = 1;

void Traits::loadInstanceParameters(ParametricData<PowerGridBranch::Instance> &p)
{
  p.addPar("AT", std::string("PQP"), &PowerGridBranch::Instance::analysisTypeStr_)
    .setUnit(U_NONE)
    .setDescription("Analysis Type");

  p.addPar("R", 0.0, &PowerGridBranch::Instance::branchResistance_)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_PERUNIT)
    .setDescription("Branch Resistance");

  p.addPar("X", 0.0, &PowerGridBranch::Instance::branchReactance_)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_PERUNIT)
    .setDescription("Branch Reactance");

  p.addPar("B", 0.0, &PowerGridBranch::Instance::branchSusceptance_)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_PERUNIT)
    .setDescription("Branch Shunt Susceptance");
}

void Traits::loadModelParameters(ParametricData<PowerGridBranch::Model> &p)
{}



std::vector< std::vector<int> > Instance::jacStamp;

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/16/14
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &               IB,
  Model & Miter,
  const FactoryBlock &  factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    analysisTypeStr_("PQP"),
    branchResistance_(0.0),
    branchReactance_(0.0),
    branchSusceptance_(0.0),
    analysisType_(PQP)
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
    analysisType_ = IV;
    analysisTypeStr_ = "IV";
  }
  else if ( atStr == "PQR" )
  {
    analysisType_ = PQR;
    analysisTypeStr_ = "PQR";
  }
  else if ( atStr == "PQP" )
  {
    analysisType_ = PQP;
    analysisTypeStr_ = "PQP";
  }
  else
  {
    UserError(*this) << "Analysis Type must be IV, PQR or PQP in power grid device.";
  }
  
  numIntVars   = 0;
  numExtVars   = 4;  
  numStateVars = 0; 
  
  if( jacStamp.empty() )
  {
    jacStamp.resize(4);
    jacStamp[0].resize(4);
    jacStamp[1].resize(4);
    jacStamp[2].resize(4);
    jacStamp[3].resize(4);

    jacStamp[0][0] = 0;  
    jacStamp[0][1] = 1;
    jacStamp[0][2] = 2;
    jacStamp[0][3] = 3;
    jacStamp[1][0] = 0;
    jacStamp[1][1] = 1;
    jacStamp[1][2] = 2;
    jacStamp[1][3] = 3;
    jacStamp[2][0] = 0;
    jacStamp[2][1] = 1;
    jacStamp[2][2] = 2;
    jacStamp[2][3] = 3;
    jacStamp[3][0] = 0;
    jacStamp[3][1] = 1;
    jacStamp[3][2] = 2;
    jacStamp[3][3] = 3;  
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
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
  // Error checking for a bad instance line that will cause some of the y parameters
  // to be infinite.
  if ( !given("R") && !given("X") )
  {
    UserError(*this) << "Either R or X must be specified for PowerGridBranch device.";
  }
  else if ( (branchResistance_ == 0.0) && (branchReactance_ == 0.0) )
  {
    UserError(*this) << "Either R or X must be non-zero for PowerGridBranch device.";
  }

  // Calculate the Y parameters from the branch data
  std::complex<double> zVal, selfVal; 

  // makes equations easier to read
  zVal=std::complex<double>(branchResistance_,branchReactance_);
  
  // off-diagonal elements are identical but specify both of them because it makes
  // the rest of the code (in the F and dF/dX functions) more readable
  y12 = (-1./zVal);
  y21 = (-1./zVal);

  // diagonal elements
  selfVal = (1./zVal);
  selfVal = std::complex<double>(selfVal.real(),selfVal.imag() + 0.5*branchSusceptance_);
  y11 = selfVal;
  y22 = selfVal;
  
  // get the real and imaginary parts.  It makes the rest of the code (in the F and dF/dX
  // functions) more readable
  g11 = y11.real();
  b11 = y11.imag();
  g12 = y12.real();
  b12 = y12.imag();
  g21 = y21.real();
  b21 = y21.imag();
  g22 = y22.real();
  b22 = y22.imag();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << getName() << " G's and B's are:"
      << Util::push << std::endl
	   << "(g11,b11) = (" << g11 << " , " << b11 << ")" << std::endl
           << "(g12,b12) = (" << g12 << " , " << b12 << ")" << std::endl
           << "(g21,b21) = (" << g21 << " , " << b21 << ")" << std::endl
           << "(g22,b22) = (" << g22 << " , " << b22 << ")" << std::endl
           << std::endl;
  }
  
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
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

  // in general, there will be 4 I/0 nodes for each branch.
  // These if-else blocks may be made into separate classes,
  // in the final version of this model.
  if ( analysisType_ == IV || analysisType_ == PQR)
  {
    li_VR1 = extLIDVec[0];
    li_VR2 = extLIDVec[1];
    li_VI1 = extLIDVec[2];
    li_VI2 = extLIDVec[3];
  
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_LIDS))
    {
      Xyce::dout() << getName() << " LIDs"
        << Util::push << std::endl
           << "li_VR1 = " << li_VR1 << std::endl
           << "li_VI1 = " << li_VI1 << std::endl
           << "li_VR2 = " << li_VR2 << std::endl
           << "li_VI2 = " << li_VI2 << std::endl
           << std::endl;
    }
  }
  else if (analysisType_ == PQP ) 
  {
    li_Theta1 = extLIDVec[0];
    li_Theta2 = extLIDVec[1];
    li_VM1 = extLIDVec[2];
    li_VM2 = extLIDVec[3];
  
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_LIDS))
    {
      Xyce::dout() << getName() << " LIDs"
        << Util::push << std::endl
           << "li_Theta1 = " << li_Theta1 << std::endl
           << "li_VM1 = " << li_VM1 << std::endl
           << "li_Theta2 = " << li_Theta2 << std::endl
           << "li_VM2 = " << li_VM2 << std::endl
           << std::endl;
    }
  }
  else
  {
    UserError(*this) << "Analysis Type must be IV, PQR or PQP in power grid device: " << getName();
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
{}

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

  // These if-else blocks may be made into separate classes,
  // in the final version of this model
  if (analysisType_ == IV || analysisType_ == PQR)
  {
    VR1_VR1_Offset = jacLIDVec[0][0];
    VR1_VR2_Offset = jacLIDVec[0][1];
    VR1_VI1_Offset = jacLIDVec[0][2];
    VR1_VI2_Offset = jacLIDVec[0][3];

    VR2_VR1_Offset = jacLIDVec[1][0];
    VR2_VR2_Offset = jacLIDVec[1][1];
    VR2_VI1_Offset = jacLIDVec[1][2];
    VR2_VI2_Offset = jacLIDVec[1][3];

    VI1_VR1_Offset = jacLIDVec[2][0];
    VI1_VR2_Offset = jacLIDVec[2][1];
    VI1_VI1_Offset = jacLIDVec[2][2];
    VI1_VI2_Offset = jacLIDVec[2][3];

    VI2_VR1_Offset = jacLIDVec[3][0];
    VI2_VR2_Offset = jacLIDVec[3][1];
    VI2_VI1_Offset = jacLIDVec[3][2];
    VI2_VI2_Offset = jacLIDVec[3][3];

    if (DEBUG_DEVICE && (isActive(Diag::DEVICE_JACSTAMP) || isActive(Diag::DEVICE_LIDS)))
    {
      Xyce::dout() << getName() << ": In registerJacLIDs, the offsets are:" << std::endl
                   << "VR1 offsets are (" << VR1_VR1_Offset << " , " << VR1_VR2_Offset << " , "
		   << VR1_VI1_Offset << " , " << VR1_VI2_Offset << ")" << std::endl
	           << "VR2 offsets are (" << VR2_VR1_Offset << " , " << VR2_VR2_Offset << " , "
		   << VR2_VI1_Offset << " , " << VR2_VI2_Offset << ")" << std::endl
		   << "VI1 offsets are (" << VI1_VR1_Offset << " , " << VI1_VR2_Offset << " , "
		   << VI1_VI1_Offset << " , " << VI1_VI2_Offset << ")" << std::endl
	           << "VI2 offsets are (" << VI2_VR1_Offset << " , " << VI2_VR2_Offset << " , "
		   << VI2_VI1_Offset << " , " << VI2_VI2_Offset << ")" << std::endl << std::endl;
    }
  }
  else if (analysisType_ == PQP)
  {
    Theta1_Theta1_Offset = jacLIDVec[0][0];
    Theta1_Theta2_Offset = jacLIDVec[0][1];
    Theta1_VM1_Offset = jacLIDVec[0][2];
    Theta1_VM2_Offset = jacLIDVec[0][3];
 
    Theta2_Theta1_Offset = jacLIDVec[1][0];
    Theta2_Theta2_Offset = jacLIDVec[1][1];
    Theta2_VM1_Offset = jacLIDVec[1][2];
    Theta2_VM2_Offset = jacLIDVec[1][3];
    
    VM1_Theta1_Offset = jacLIDVec[2][0];
    VM1_Theta2_Offset = jacLIDVec[2][1];
    VM1_VM1_Offset = jacLIDVec[2][2];
    VM1_VM2_Offset = jacLIDVec[2][3];

    VM2_Theta1_Offset = jacLIDVec[3][0];
    VM2_Theta2_Offset = jacLIDVec[3][1];
    VM2_VM1_Offset = jacLIDVec[3][2];
    VM2_VM2_Offset = jacLIDVec[3][3];

    if (DEBUG_DEVICE && (isActive(Diag::DEVICE_JACSTAMP) || isActive(Diag::DEVICE_LIDS)))
    {
      Xyce::dout() << getName() << ": In registerJacLIDs, the offsets are:" << std::endl
                   << "Theta1 offsets are (" << Theta1_Theta1_Offset << " , " << Theta1_Theta2_Offset << " , "
		   << Theta1_VM1_Offset << " , " << Theta1_VM2_Offset << ")" << std::endl
	           << "Theta2 offsets are (" << Theta2_Theta1_Offset << " , " << Theta2_Theta2_Offset << " , "
		   << Theta2_VM1_Offset << " , " << Theta2_VM2_Offset << ")" << std::endl
		   << "VM1 offsets are (" << VM1_Theta1_Offset << " , " << VM1_Theta2_Offset << " , "
		   << VM1_VM1_Offset << " , " << VM1_VM2_Offset << ")" << std::endl
	           << "VM2 offsets are ("  << VM2_Theta1_Offset << " , " << VM2_Theta2_Offset << " , "
		   << VM2_VM1_Offset << " , " << VM2_VM2_Offset << ")" << std::endl<< std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one power grid branch instance
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  double * solVec = extData.nextSolVectorRawPtr;

  // These if-else blocks may be split out into separate classes in the final version
  // of this model
  if (analysisType_ == IV)
  {
    VR1 = solVec[li_VR1];
    VR2 = solVec[li_VR2];
    VI1 = solVec[li_VI1];
    VI2 = solVec[li_VI2];

    IR1 = g11*VR1 +  g12*VR2 - b11*VI1 - b12*VI2;
    IR2 = g21*VR1 +  g22*VR2 - b21*VI1 - b22*VI2;
    II1 = b11*VR1 +  b12*VR2 + g11*VI1 + g12*VI2;
    II2 = b21*VR1 +  b22*VR2 + g21*VI1 + g22*VI2;
  }
  else if (analysisType_ == PQR)
  {
    VR1 = solVec[li_VR1];
    VR2 = solVec[li_VR2];
    VI1 = solVec[li_VI1];
    VI2 = solVec[li_VI2];

    P1 = g11*(VR1*VR1 + VI1*VI1) + VR1*(g12*VR2 - b12*VI2) + VI1*(b12*VR2 + g12*VI2);
    P2 = g22*(VR2*VR2 + VI2*VI2) + VR2*(g21*VR1 - b21*VI1) + VI2*(b21*VR1 + g21*VI1);

    Q1 = -1.0*b11*(VR1*VR1 + VI1*VI1) + VI1*(g12*VR2 - b12*VI2) - VR1*(b12*VR2 + g12*VI2);
    Q2 = -1.0*b22*(VR2*VR2 + VI2*VI2) + VI2*(g21*VR1 - b21*VI1) - VR2*(b21*VR1 + g21*VI1);
  }
  else if (analysisType_ == PQP)
  { 
    VM1 = solVec[li_VM1];
    VM2 = solVec[li_VM2];
    Theta1 = solVec[li_Theta1];
    Theta2 = solVec[li_Theta2];
   
    // these variables make the equations easier to read
    dSin12 = sin(Theta1 - Theta2);
    dSin21 = sin(Theta2 - Theta1);
    dCos12 = cos(Theta1 - Theta2);
    dCos21 = cos(Theta2 - Theta1);
  
    P1 = g11*VM1*VM1 + VM1*VM2*(g12*dCos12 + b12*dSin12);
    P2 = g22*VM2*VM2 + VM2*VM1*(g21*dCos21 + b21*dSin21);

    Q1 = -1*b11*VM1*VM1 + VM1*VM2*(g12*dSin12 - b12*dCos12);
    Q2 = -1*b22*VM2*VM2 + VM2*VM1*(g21*dSin21 - b21*dCos21);
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
  
  // These if-else blocks may be split out into separate classes in the final version
  // of this mode.
  if (analysisType_ == IV)
  {
    fVec[li_VR1] += IR1;
    fVec[li_VR2] += IR2; 
    fVec[li_VI1] += II1; 
    fVec[li_VI2] += II2;  
  }
  else if (analysisType_ == PQR)
  {
    fVec[li_VR1] += P1;
    fVec[li_VR2] += P2; 
    fVec[li_VI1] += Q1; 
    fVec[li_VI2] += Q2;  
  }
  else if (analysisType_ == PQP)
  {
    fVec[li_Theta1] += P1;
    fVec[li_Theta2] += P2; 
    fVec[li_VM1] += Q1; 
    fVec[li_VM2] += Q2;  

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_LOAD_VECTOR))
    {
      Xyce::dout() << getName() << ": F Vector Load info is: (P1,P2,Q1,Q2) load is: (" << P1 << " , "
                   << P2 << " , " << Q1 << " , " << Q2 << ")" << std::endl << std::endl;     
    }
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
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions.  This is a no-op for
//                 this model
//
// Special Notes :
//
//
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
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

  // These if-else blocks may be split out into separate classes in the final version
  // of this model
  if (analysisType_ == IV)
  {
    dFdx[li_VR1][VR1_VR1_Offset] += g11;
    dFdx[li_VR1][VR1_VR2_Offset] += g12;
    dFdx[li_VR1][VR1_VI1_Offset] -= b11;
    dFdx[li_VR1][VR1_VI2_Offset] -= b12;

    dFdx[li_VR2][VR2_VR1_Offset] += g21;
    dFdx[li_VR2][VR2_VR2_Offset] += g22;
    dFdx[li_VR2][VR2_VI1_Offset] -= b21;
    dFdx[li_VR2][VR2_VI2_Offset] -= b22;

    dFdx[li_VI1][VI1_VR1_Offset] += b11;
    dFdx[li_VI1][VI1_VR2_Offset] += b12;
    dFdx[li_VI1][VI1_VI1_Offset] += g11;
    dFdx[li_VI1][VI1_VI2_Offset] += g12;

    dFdx[li_VI2][VI2_VR1_Offset] += b21;
    dFdx[li_VI2][VI2_VR2_Offset] += b22;
    dFdx[li_VI2][VI2_VI1_Offset] += g21;
    dFdx[li_VI2][VI2_VI2_Offset] += g22;
  }
  else if (analysisType_ == PQR)
  {
    // VR1, VR2, VI1 and VI2 were calculated in updateIntermediateVars()
    dFdx[li_VR1][VR1_VR1_Offset] += 2*g11*VR1 + g12*VR2 - b12*VI2;
    dFdx[li_VR1][VR1_VR2_Offset] += g12*VR1 + b12*VI1;
    dFdx[li_VR1][VR1_VI1_Offset] += 2*g11*VI1 + b12*VR2 + g12*VI2;
    dFdx[li_VR1][VR1_VI2_Offset] += g12*VI1 - b12*VR1;

    dFdx[li_VR2][VR2_VR1_Offset] += g21*VR2 + b21*VI2;
    dFdx[li_VR2][VR2_VR2_Offset] += 2*g22*VR2 + g21*VR1 - b21*VI1;
    dFdx[li_VR2][VR2_VI1_Offset] += g21*VI2 - b21*VR2;
    dFdx[li_VR2][VR2_VI2_Offset] += 2*g22*VI2 + b21*VR1 + g21*VI1;

    dFdx[li_VI1][VI1_VR1_Offset] += -2*b11*VR1 - b12*VR2 - g12*VI2;
    dFdx[li_VI1][VI1_VR2_Offset] += g12*VI1 - b12*VR1;
    dFdx[li_VI1][VI1_VI1_Offset] += -2*b11*VI1 + g12*VR2 - b12*VI2;
    dFdx[li_VI1][VI1_VI2_Offset] -= g12*VR1 + b12*VI1;

    dFdx[li_VI2][VI2_VR1_Offset] += g21*VI2 - b21*VR2;
    dFdx[li_VI2][VI2_VR2_Offset] += -2*b22*VR2 - b21*VR1 - g12*VI1;
    dFdx[li_VI2][VI2_VI1_Offset] -= g21*VR2 + b21*VI2; 
    dFdx[li_VI2][VI2_VI2_Offset] += -2*b22*VI2 + g21*VR1 - b21*VI1; 
  }
  else if (analysisType_ == PQP)
  {
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_LOAD_VECTOR))
    {
      Xyce::dout() << getName() << ": dFdx info is: " << std::endl
	   << "(dSin12,dSin21,dCos12,dCos21 = (" << dSin12 << " , " << dSin21  
	   << " , " << dCos12 << " , " << dCos21 <<  ")" << std::endl
           << "(Theta1,Theta2) = (" << Theta1 << " , " << Theta2 << ")" << std::endl
           << "(VM1,VM2) = (" << VM1 << " , " << VM2 << ")" << std::endl;
    }

    // VM1, VM2, Theta1, Theta2, dSin12, dSin21, dCos12 and dCos21 were
    // calculated in updateIntermediateVars()  

    // dP / dTheta terms
    dFdx[li_Theta1][Theta1_Theta1_Offset] -=  VM1*VM2*(g12*dSin12 - b12*dCos12);
    dFdx[li_Theta1][Theta1_Theta2_Offset] +=  VM1*VM2*(g12*dSin12 - b12*dCos12);	 
    dFdx[li_Theta2][Theta2_Theta1_Offset] +=  VM1*VM2*(g21*dSin21 - b21*dCos21);
    dFdx[li_Theta2][Theta2_Theta2_Offset] -=  VM1*VM2*(g21*dSin21 - b21*dCos21);

    // dP / dV terms
    dFdx[li_Theta1][Theta1_VM1_Offset] += 2*VM1*g11 + VM2*(g12*dCos12 + b12*dSin12);  
    dFdx[li_Theta1][Theta1_VM2_Offset] += VM1*(g12*dCos12 + b12*dSin12);
    dFdx[li_Theta2][Theta2_VM1_Offset] += VM2*(g21*dCos21 + b21*dSin21); 
    dFdx[li_Theta2][Theta2_VM2_Offset] += 2*VM2*g22 + VM1*(g21*dCos21 + b21*dSin21);  

    // dQ / dTheta terms
    dFdx[li_VM1][VM1_Theta1_Offset] +=  VM1*VM2*(g12*dCos12 + b12*dSin12);
    dFdx[li_VM1][VM1_Theta2_Offset] -=  VM1*VM2*(g12*dCos12 + b12*dSin12);
    dFdx[li_VM2][VM2_Theta1_Offset] -=  VM1*VM2*(g21*dCos21 + b21*dSin21);
    dFdx[li_VM2][VM2_Theta2_Offset] +=  VM1*VM2*(g21*dCos21 + b21*dSin21);

    // dQ / dV terms
    dFdx[li_VM1][VM1_VM1_Offset] += -2*VM1*b11 + VM2*(g12*dSin12 - b12*dCos12);  
    dFdx[li_VM1][VM1_VM2_Offset] += VM1*(g12*dSin12 - b12*dCos12);
    dFdx[li_VM2][VM2_VM1_Offset] += VM2*(g21*dSin21 - b21*dCos21);
    dFdx[li_VM2][VM2_VM2_Offset] += -2*VM2*b22 + VM1*(g21*dSin21 - b21*dCos21);
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
// Function      : Instance::loadDAEdQdx ()
//
// Purpose       : Loads the dQdx-matrix contributions.  This is no-op for this
//                 model.
//
// Special Notes :
//
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
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
  if (deviceMap.empty() || (deviceMap.find("POWERGRIDBRANCH")!=deviceMap.end())
                        || (deviceMap.find("PGBR")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("PowerGridBranch", 1)
      .registerDevice("PGBR", 1)
      .registerModelType("PowerGridBranch", 1);
  }
}

} // namespace PowerGridBranch
} // namespace Device
} // namespace Xyce

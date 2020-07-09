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

//----------------------------------------------------------------------------
//
// Purpose        : PowerGrid classes: provides a device that calculates the
//                  steady state power flow in a transmission grid branch
//
// Special Notes  : Experimental new device for an LDRD.
//
// Creator        : Pete Sholander, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 2/22/15
//
//
//
//
//----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_PowerGridTransformer.h>
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
namespace PowerGridTransformer {

// enables debug output for just this device
//static const int DEBUG_DEVICE = 1;

void Traits::loadInstanceParameters(ParametricData<PowerGridTransformer::Instance> &p)
{
  p.addPar("AT", std::string("PQP"), &PowerGridTransformer::Instance::analysisTypeStr_)
    .setUnit(U_NONE)
    .setDescription("Analysis Type");

  p.addPar("R", 0.0, &PowerGridTransformer::Instance::resistance_)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_PERUNIT)
    .setDescription("Resistance");

  p.addPar("X", 0.0, &PowerGridTransformer::Instance::reactance_)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_PERUNIT)
    .setDescription("Reactance");

  p.addPar("TR", 1.0, &PowerGridTransformer::Instance::turnsRatio_)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_PERUNIT)
    .setDescription("Transformer Turns Ratio");

  p.addPar("PS", 0.0, &PowerGridTransformer::Instance::phaseShift_)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_RAD)
    .setDescription("Phase Shift");
  
  p.addPar("TT", std::string("FT"), &PowerGridTransformer::Instance::transTypeStr_)
    .setUnit(U_NONE)
    .setDescription("Transformer Type");
}

void Traits::loadModelParameters(ParametricData<PowerGridTransformer::Model> &p)
{}


// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 2/22/15
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &               IB,
  Model & Miter,
  const FactoryBlock &  factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    analysisTypeStr_("PQP"),
    transTypeStr_("FT"),
    resistance_(0.0),
    reactance_(0.0),
    turnsRatio_(1.0),
    phaseShift_(0.0),
    analysisType_(PQP),
    transType_(FT)
{
  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  processParams();
  
  // Finish setting up the device.  Process the atStr enum in processParams,
  // rather than here like for the branch, shunt and genbus, since the 
  // g's and b's for the transformer depend on the analysis type in 
  // some cases.
 
  // The Fixed Tap (FT) transformer uses fixed values for the turns ratio and phase shift.  Those
  // values are set on the instance line.  The Variable Tap (VT) transformer has the real part
  // of the turns ratio (N) as a solution variable, using an internal variable.  
  // The phase shift is fixed, and specified on the instance line for both FT
  // and VT transformers.  The Phase Shift (PS) transformer has the real part of the turns
  // ratio (N) fixed but the phase shift (Phi) is now a solution variable.
  if ( (transType_ == VT) || (transType_ == PS) )
  { 
    numExtVars = 5;

    // used by the Topology Checker.  See SON Bug 974 for more details.
    // The optional control-node terminal does not have a DC path to the
    // four input/output nodes.
    devConMap.resize(5);
    devConMap[0]=1;
    devConMap[1]=1;
    devConMap[2]=1;
    devConMap[3]=1;
    devConMap[4]=2;
  }
  else
  {
    numExtVars   = 4;

    devConMap.resize(4);
    devConMap[0]=1;
    devConMap[1]=1;
    devConMap[2]=1;
    devConMap[3]=1;
  }  
  numStateVars = 0; 

  // check for errors in the number of nodes on the instance line
  if (numExtVars != IB.numExtVars)
  {
    UserError(*this) << "Incorrect number of inputs in power grid device."
                      << " Found " << IB.numExtVars << ", should be " << numExtVars  
                      << " for transformer type " << transTypeStr_ << " power grid transformer.";
  }

  // Since there are multiple transformer types, each instance may have a 
  // different sized jacStamp matrix.  So, we must make this matrix for
  // each instance.  It has four rows for all three transformer types.
  jacStamp.resize(4); 
  
  // number of columns depends on the transformer type though.
  if ( (transType_ == VT) || (transType_ == PS) )
  {
    // the Jacobian has a fifth column, for the N variable for a
    // variable tap (VT) transformer or for the Phi variable for a phase
    // shifting (PS) transformer.  It still has four rows though.
    jacStamp[0].resize(5);
    jacStamp[1].resize(5);
    jacStamp[2].resize(5);
    jacStamp[3].resize(5);

    jacStamp[0][4] = 4;
    jacStamp[1][4] = 4;
    jacStamp[2][4] = 4;
    jacStamp[3][4] = 4;
  }     
  else
  {   
    // the default Jacobian, for a Fixed Tap (FT) transformer with a 
    // fixed turns ratio, is a 4x4 matrix.
    jacStamp[0].resize(4);
    jacStamp[1].resize(4);
    jacStamp[2].resize(4);
    jacStamp[3].resize(4);
  }

  // fill out the remaining entries, in the leftmost 4x4 portion of
  // the Jacobian matrix, that are common to all three transformer types.
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
    UserError(*this) << "R or X must be specified for PowerGridTransformer device.";
  }
  else if ( (resistance_ == 0.0) && (reactance_ == 0.0) )
  {
    UserError(*this) << "Either R or X must be non-zero for PowerGridTransformer device.";
  }
  else if ( turnsRatio_ == 0.0 )
  {
    UserError(*this) << "TR must be non-zero for PowerGridTransformer device.";
  }

  // process the analysisType_ and transformerType_ enums here because the analysisType_ 
  // and transType_ are used later to set the b's and g's for the transformer. 
  // This is done in the constructor for the Branch, GenBus and Shunt device models.
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

  ExtendedString ttStr(transTypeStr_);
  ttStr.toUpper();
  if ( ttStr == "FT" )
  {
    transType_ = FT;
    transTypeStr_ = "FT";
  }
  else if ( ttStr == "VT" )
  {
    transType_ = VT;
    transTypeStr_ = "VT";
  }
  else if ( ttStr == "PS" )
  {
    transType_ = PS;
    transTypeStr_ = "PS";
  }
  else
  {
    UserError(*this) << "Transformer Type must be FT, VT or PS in PowerGridTransformer device.";
  }

  // Calculate the Y parameters from the branch data. These initial calculations
  // omit the phase shift.  That is added below.
  std::complex<double> zVal; 

  // makes equations easier to read.  The turns ratio will be accounted for in
  // updateIntermediateVars() for all three solution formats.  This allows the
  // turns ratio to be an input for the Variable Tap (VT) transformer
  invTurnsRatio_ = (1./turnsRatio_);
  zVal=std::complex<double>(resistance_,reactance_);
  
  // off-diagonal elements are identical but specify both of them because it makes
  // the rest of the code (in the F and dF/dX functions) more readable
  y12 = (-1./zVal);
  y21 = (-1./zVal);

  // diagonal elements
  y11 = (1./zVal);
  y22 = (1./zVal);
  
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

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && transType_ == 1)
  {
    Xyce::dout() << getName() << " G's and B's without the phase shift correction are:"
      << Util::push << std::endl
	   << "(g11,b11) = (" << g11 << " , " << b11 << ")" << std::endl
           << "(g12,b12) = (" << g12 << " , " << b12 << ")" << std::endl
           << "(g21,b21) = (" << g21 << " , " << b21 << ")" << std::endl
           << "(g22,b22) = (" << g22 << " , " << b22 << ")" << std::endl
	   << "(Turns Ratio, Phase Shift) = (" << turnsRatio_ << "," << phaseShift_ << ")" << std::endl
           << std::endl;
  }

  // Apply the phase shift directly to the g's and b's for the I=YV and PQ Rectangular 
  // formulations for the FT and VT transformers.  The phase shift is added 
  // in updateIntermediateVars() for the PQ Polar format for all transformer types, 
  // and for the PS (phase shift) transformer for I=YV and PQ Rectangular formats.
  if ( (analysisType_ == IV || analysisType_ == PQR) && (phaseShift_ != 0.0) &&
       ( (transType_ == FT) || (transType_ == VT)) )
  {
    // see pp. 245-247 of Kundar's text book "Power System Stability and Control" for
    // a derivation of these equations for g12, b12, g21 and b21
    double prevG12 = g12;
    double prevG21 = g21;
    
    g12 = g12*cos(phaseShift_) - b12*sin(phaseShift_);
    b12 = b12*cos(phaseShift_) + prevG12*sin(phaseShift_);

    g21 = g21*cos(phaseShift_) + b21*sin(phaseShift_);
    b21 = b21*cos(phaseShift_) - prevG21*sin(phaseShift_);

    if ( DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
    {
      Xyce::dout() << getName() << ": for " << analysisTypeStr_ << " analysis, G's and B's after phase shift correction are:"
        << Util::push << std::endl
	   << "(g11,b11) = (" << g11 << " , " << b11 << ")" << std::endl
           << "(g12,b12) = (" << g12 << " , " << b12 << ")" << std::endl
           << "(g21,b21) = (" << g21 << " , " << b21 << ")" << std::endl
           << "(g22,b22) = (" << g22 << " , " << b22 << ")" << std::endl
	   << "(Turns Ratio, Phase Shift) = (" << turnsRatio_ << "," << phaseShift_ << ")" << std::endl
           << std::endl;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/12/14
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef )
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // There will be 4 I/0 nodes for a Fixed Tap (FT) Transformer, and 
  // one additional internal node for the N control input for a 
  // Variable Tap (VT) transformer or for the Phi control input for a
  // Phase Shift (PS) transformer
  if ( analysisType_ == IV || analysisType_ == PQR)
  {
    li_VR1 = extLIDVec[0];
    li_VR2 = extLIDVec[1];
    li_VI1 = extLIDVec[2];
    li_VI2 = extLIDVec[3];
    if ( transType_ == VT )  
    { 
      li_N = extLIDVec[4]; 
    }
    else if ( transType_ == PS )  
    { 
      li_Phi = extLIDVec[4]; 
    }
  
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_LIDS))
    {
      Xyce::dout() << getName() << " LIDs are"
         << Util::push << std::endl
         << "li_VR1 = " << li_VR1 << std::endl
         << "li_VI1 = " << li_VI1 << std::endl
         << "li_VR2 = " << li_VR2 << std::endl
	 << "li_VI2 = " << li_VI2 << std::endl;
      if (transType_ == VT)
      {
	Xyce::dout() << "li_N = " << li_N << std::endl;
      }
      else if (transType_ == PS)
      {
        Xyce::dout() << "li_Phi = " << li_Phi << std::endl;
      }

      Xyce::dout() << std::endl;
    }
  }
  else if (analysisType_ == PQP ) 
  {
    li_Theta1 = extLIDVec[0];
    li_Theta2 = extLIDVec[1];
    li_VM1 = extLIDVec[2];
    li_VM2 = extLIDVec[3];
    if ( transType_ == VT )  
    { 
      li_N = extLIDVec[4]; 
    }
    else if ( transType_ == PS )  
    { 
      li_Phi = extLIDVec[4]; 
    }
  
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_LIDS))
    {
      Xyce::dout() << getName() << " LIDs"
         << Util::push << std::endl
         << "li_Theta1 = " << li_Theta1 << std::endl
         << "li_VM1 = " << li_VM1 << std::endl
         << "li_Theta2 = " << li_Theta2 << std::endl
	 << "li_VM2 = " << li_VM2 << std::endl;

      if (transType_ == VT)
      {
        Xyce::dout() << "li_N = " << li_N << std::endl;
      }
      else if (transType_ == PS)
      {
        Xyce::dout() << "li_Phi = " << li_Phi << std::endl;
      }

      Xyce::dout() << std::endl;
    }
  }
  else
  {
    UserError(*this) << "Analysis Type must be IV, PQR or PQP in power grid device";
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
  // the Variable Tap (VT) and Phase Shift (PS) transformer have a fifth column
  // in their Jacobian for the control input (either N or Phi)
  // the extra node symbol for N or Phi does not need to be added since it is external
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

    if ( transType_ == VT )
    {
      VR1_N_Offset = jacLIDVec[0][4];
      VR2_N_Offset = jacLIDVec[1][4];
      VI1_N_Offset = jacLIDVec[2][4];
      VI2_N_Offset = jacLIDVec[3][4];
    }
    else if ( transType_ == PS )
    {
      VR1_Phi_Offset = jacLIDVec[0][4];
      VR2_Phi_Offset = jacLIDVec[1][4];
      VI1_Phi_Offset = jacLIDVec[2][4];
      VI2_Phi_Offset = jacLIDVec[3][4];
    }

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

    if ( transType_ == VT )
    {
      Theta1_N_Offset = jacLIDVec[0][4];
      Theta2_N_Offset = jacLIDVec[1][4];
      VM1_N_Offset = jacLIDVec[2][4];
      VM2_N_Offset = jacLIDVec[3][4];
    }
    else if ( transType_ == PS )
    {
      Theta1_Phi_Offset = jacLIDVec[0][4];
      Theta2_Phi_Offset = jacLIDVec[1][4];
      VM1_Phi_Offset = jacLIDVec[2][4];
      VM2_Phi_Offset = jacLIDVec[3][4];
    }

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
// Purpose       : update intermediate variables for one power grid
//                 transformer instance
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  double * solVec = extData.nextSolVectorRawPtr;

  // The g's and b's are "adjusted" because of the turns ratio and/or phase shift.  For 
  // a fixed tap (FT) transformer, the turns ratio (N) is fixed and based on the instance parameter (TR).
  // For a variable tap (VT) transformer, N is an input that may be time-varying.  The equations
  // in loadDAEdFdx() use the adjusted parameters.  The adjustment for N is the same
  // for all three solution formats.
  if (transType_ == FT)
  {
    // turns ratio is fixed
    ga11 = g11*(invTurnsRatio_*invTurnsRatio_);
    ba11 = b11*(invTurnsRatio_*invTurnsRatio_);
    ga12 = g12*(invTurnsRatio_);
    ba12 = b12*(invTurnsRatio_);
    ga21 = g21*(invTurnsRatio_);
    ba21 = b21*(invTurnsRatio_);
    ga22 = g22;
    ba22 = b22;
  }
  else if ( transType_ == VT)
  {
    // The value of the solution variable (N) may be zero at the first newton step.
    // So, invN_ is set to the inverse of the TR instance parameter for 
    // the first newton step. 
    //N_ = (getSolverState().newtonIter == 0) ? turnsRatio_ : solVec[li_N];
    N_ = solVec[li_N];
    invN_ = 1/N_;
    ga11 = g11*(invN_*invN_);
    ba11 = b11*(invN_*invN_);
    ga12 = g12*(invN_);
    ba12 = b12*(invN_);
    ga21 = g21*(invN_);
    ba21 = b21*(invN_);
    ga22 = g22;
    ba22 = b22;
  }
  else if ( transType_ == PS )
  {
    Phi_ = solVec[li_Phi];
    if ( (analysisType_ == IV) || (analysisType_ == PQR) )
    {
      ga11 = g11*(invTurnsRatio_*invTurnsRatio_);
      ba11 = b11*(invTurnsRatio_*invTurnsRatio_);
      ga12 = ( g12*cos(Phi_) - b12*sin(Phi_) )*invTurnsRatio_;
      ba12 = ( b12*cos(Phi_) + g12*sin(Phi_) )*invTurnsRatio_;
      ga21 = ( g21*cos(Phi_) + b21*sin(Phi_) )*invTurnsRatio_;
      ba21 = ( b21*cos(Phi_) - g21*sin(Phi_) )*invTurnsRatio_;
      ga22 = g22;
      ba22 = b22;
    }
    else if ( analysisType_ == PQP )
    {
      // phase shift (Phi) will be added later, then the power injections are
      // calculated
      ga11 = g11*(invTurnsRatio_*invTurnsRatio_);
      ba11 = b11*(invTurnsRatio_*invTurnsRatio_);
      ga12 = g12*(invTurnsRatio_);
      ba12 = b12*(invTurnsRatio_);
      ga21 = g21*(invTurnsRatio_);
      ba21 = b21*(invTurnsRatio_);
      ga22 = g22;
      ba22 = b22;
    }
  }
  else
  {
    UserError(*this) << "Transformer Type must be 1, 2 or 3: " << getName();
    return false;
  }
  
  // power flow equations use the "adjusted" values for g's and b's that account 
  // for the value of the turns ratio (N) that may be time-varying. 
  if (analysisType_ == IV) 
  {
    // For IV and PQ Rectangular formats, g12a, b12a, g21a and b21a were updated 
    // to include both the turnsRatio_ and the phaseShift_ 
    VR1 = solVec[li_VR1];
    VR2 = solVec[li_VR2];
    VI1 = solVec[li_VI1];
    VI2 = solVec[li_VI2];

    IR1 = ga11*VR1 +  ga12*VR2 - ba11*VI1 - ba12*VI2;
    IR2 = ga21*VR1 +  ga22*VR2 - ba21*VI1 - ba22*VI2;
    II1 = ba11*VR1 +  ba12*VR2 + ga11*VI1 + ga12*VI2;
    II2 = ba21*VR1 +  ba22*VR2 + ga21*VI1 + ga22*VI2;
  }
  else if (analysisType_ == PQR)
  {
    // For IV and PQ Rectangular formats, g12a, b12a, g21a and b21a were updated  
    // to include both the turnsRatio_ and the phaseShift_ 
    VR1 = solVec[li_VR1];
    VR2 = solVec[li_VR2];
    VI1 = solVec[li_VI1];
    VI2 = solVec[li_VI2];

    P1 = ga11*(VR1*VR1 + VI1*VI1) + VR1*(ga12*VR2 - ba12*VI2) + VI1*(ba12*VR2 + ga12*VI2);
    P2 = ga22*(VR2*VR2 + VI2*VI2) + VR2*(ga21*VR1 - ba21*VI1) + VI2*(ba21*VR1 + ga21*VI1);

    Q1 = -1.0*ba11*(VR1*VR1 + VI1*VI1) + VI1*(ga12*VR2 - ba12*VI2) - VR1*(ba12*VR2 + ga12*VI2);
    Q2 = -1.0*ba22*(VR2*VR2 + VI2*VI2) + VI2*(ga21*VR1 - ba21*VI1) - VR2*(ba21*VR1 + ga21*VI1);
  }
  else if (analysisType_ == PQP)
  {    
    // phaseShift_ is on the primary side of the transformer.  So, it is subtracted from dSin12
    // and dCos12.  It is added to the dSin21 and dCos21 terms.  The rest of the equations are
    // then unchanged.  For the derivation, see eqn.11.27 on pg. 272-274 of Milano's textbook
    // "Power System Modeling and Scripting"
    VM1 = solVec[li_VM1];
    VM2 = solVec[li_VM2];
    Theta1 = solVec[li_Theta1];
    Theta2 = solVec[li_Theta2];
   
    // these variables make the equations easier to read
    if ( (transType_ == FT) || (transType_ == VT) ) 
    {
      dSin12 = sin(Theta1 - Theta2 - phaseShift_);
      dSin21 = sin(Theta2 - Theta1 + phaseShift_);
      dCos12 = cos(Theta1 - Theta2 - phaseShift_);
      dCos21 = cos(Theta2 - Theta1 + phaseShift_);
    }
    else if ( transType_ == PS )  
    {
      dSin12 = sin(Theta1 - Theta2 - Phi_);
      dSin21 = sin(Theta2 - Theta1 + Phi_);
      dCos12 = cos(Theta1 - Theta2 - Phi_);
      dCos21 = cos(Theta2 - Theta1 + Phi_);
    }

    P1 = ga11*VM1*VM1 + VM1*VM2*(ga12*dCos12 + ba12*dSin12);
    P2 = ga22*VM2*VM2 + VM2*VM1*(ga21*dCos21 + ba21*dSin21);

    Q1 = -1*ba11*VM1*VM1 + VM1*VM2*(ga12*dSin12 - ba12*dCos12);
    Q2 = -1*ba22*VM2*VM2 + VM2*VM1*(ga21*dSin21 - ba21*dCos21);
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

  // These if-else blocks may be made into separate classes,
  // in the final version of this model
  if (analysisType_ == IV)
  {   
    // For IV and PQ Rectangular formats, g12, b12, g21 and b21 were updated in the 
    // constructor to include only the phaseShift_, but for the FT and VT transformers.
    dFdx[li_VR1][VR1_VR1_Offset] += ga11;
    dFdx[li_VR1][VR1_VR2_Offset] += ga12;
    dFdx[li_VR1][VR1_VI1_Offset] -= ba11;
    dFdx[li_VR1][VR1_VI2_Offset] -= ba12;

    dFdx[li_VR2][VR2_VR1_Offset] += ga21;
    dFdx[li_VR2][VR2_VR2_Offset] += ga22;
    dFdx[li_VR2][VR2_VI1_Offset] -= ba21;
    dFdx[li_VR2][VR2_VI2_Offset] -= ba22;

    dFdx[li_VI1][VI1_VR1_Offset] += ba11;
    dFdx[li_VI1][VI1_VR2_Offset] += ba12;
    dFdx[li_VI1][VI1_VI1_Offset] += ga11;
    dFdx[li_VI1][VI1_VI2_Offset] += ga12;

    dFdx[li_VI2][VI2_VR1_Offset] += ba21;
    dFdx[li_VI2][VI2_VR2_Offset] += ba22;
    dFdx[li_VI2][VI2_VI1_Offset] += ga21;
    dFdx[li_VI2][VI2_VI2_Offset] += ga22;

    // also need to load the last column with the derivatives wrt the turns ratio,
    // if this is a Variable Tap (VT) transformer.
    if (transType_ == VT)
    {
      dFdx[li_VR1][VR1_N_Offset] -= ( 2*ga11*VR1*invN_ + ga12*VR2*invN_ 
                                    - 2*ba11*VI1*invN_ - ba12*VI2*invN_ );
      dFdx[li_VR2][VR2_N_Offset] -= ( ga21*VR1*invN_ - ba21*VI1*invN_ );
      dFdx[li_VI1][VI1_N_Offset] -= ( 2*ba11*VR1*invN_ + ba12*VR2*invN_ 
                                    + 2*ga11*VI1*invN_ + ga12*VI2*invN_ );
      dFdx[li_VI2][VI2_N_Offset] -= ( ba21*VR1*invN_ + ga21*VI1*invN_ );
    }
    else if (transType_ == PS)
    {
      // note: this else-if block uses g12 rather than ga12, etc.  That's because ga12 
      // as defined in updateIntermediateVars() is a function of Phi_.  So, this may need to
      // be re-factored at some point.
      dFdx[li_VR1][VR1_Phi_Offset] -= invTurnsRatio_*(g12*VR2*sin(Phi_) + b12*VR2*cos(Phi_) 
						      - b12*VI2*sin(Phi_) + g12*VI2*cos(Phi_) );
      dFdx[li_VR2][VR2_Phi_Offset] -= invTurnsRatio_*( g21*VR2*sin(Phi_) - b21*VR2*cos(Phi_) 
                                                      - b21*VI1*sin(Phi_) - g21*VI1*cos(Phi_) );
      dFdx[li_VI1][VI1_Phi_Offset] -= invTurnsRatio_*( b12*VR2*sin(Phi_) - g12*VR2*cos(Phi_) 
						      + g12*VI2*cos(Phi_) + b12*VI2*sin(Phi_) );
      dFdx[li_VI2][VI2_Phi_Offset] -= invTurnsRatio_*( b21*VR1*sin(Phi_) + g21*VR1*cos(Phi_) 
                                                       + g12*VI1*cos(Phi_) - b12*VI1*sin(Phi_));
    }
  }
  else if (analysisType_ == PQR) 
  {
    // For IV and PQ Rectangular formats, g12, b12, g21 and b21 were updated in the 
    // constructor to include only the phaseShift_ for the FT and VT transformers.
    dFdx[li_VR1][VR1_VR1_Offset] += 2*ga11*VR1 + ga12*VR2 - ba12*VI2;
    dFdx[li_VR1][VR1_VR2_Offset] += ga12*VR1 + ba12*VI1;
    dFdx[li_VR1][VR1_VI1_Offset] += 2*ga11*VI1 + ba12*VR2 + ga12*VI2;
    dFdx[li_VR1][VR1_VI2_Offset] += ga12*VI1 - ba12*VR1;

    dFdx[li_VR2][VR2_VR1_Offset] += ga21*VR2 + ba21*VI2;
    dFdx[li_VR2][VR2_VR2_Offset] += 2*ga22*VR2 + ga21*VR1 - ba21*VI1;
    dFdx[li_VR2][VR2_VI1_Offset] += ga21*VI2 - ba21*VR2;
    dFdx[li_VR2][VR2_VI2_Offset] += 2*ga22*VI2 + ba21*VR1 + ga21*VI1;

    dFdx[li_VI1][VI1_VR1_Offset] += -2*ba11*VR1 - ba12*VR2 - ga12*VI2;
    dFdx[li_VI1][VI1_VR2_Offset] += ga12*VI1 - ba12*VR1;
    dFdx[li_VI1][VI1_VI1_Offset] += -2*ba11*VI1 + ga12*VR2 - ba12*VI2;
    dFdx[li_VI1][VI1_VI2_Offset] -= ga12*VR1 + ba12*VI1;

    dFdx[li_VI2][VI2_VR1_Offset] += ga21*VI2 - ba21*VR2;
    dFdx[li_VI2][VI2_VR2_Offset] += -2*ba22*VR2 - ba21*VR1 - ga12*VI1;
    dFdx[li_VI2][VI2_VI1_Offset] -= ga21*VR2 + ba21*VI2; 
    dFdx[li_VI2][VI2_VI2_Offset] += -2*ba22*VI2 + ga21*VR1 - ba21*VI1; 

    if (transType_ == VT)
    {
      dFdx[li_VR1][VR1_N_Offset] -= ( ga11*(2*invN_)*(VR1*VR1 + VI1*VI1)
                                      + ga12*invN_*VR1*VR2 - ba12*invN_*VR1*VI2
                                      + ba12*invN_*VI1*VR2 + ga12*invN_*VI1*VI2);

      dFdx[li_VR2][VR2_N_Offset] -= ( ga21*invN_*VR2*VR1 - ba21*invN_*VR2*VI1
                                      + ba21*invN_*VI2*VR1 + ga21*invN_*VI2*VI1);

      dFdx[li_VI1][VI1_N_Offset] -= ( -ba11*(2*invN_)*(VR1*VR1 + VI1*VI1)
                                      + ga12*invN_*VI1*VR2 - ba12*invN_*VI1*VI2
                                      + ba12*invN_*VR1*VR2 + ga12*invN_*VR1*VI2);

      dFdx[li_VI2][VI2_N_Offset] -= ( ga21*invN_*VI2*VR1 - ba21*invN_*VI2*VI1
                                      - ba21*invN_*VR2*VR1 - ga21*invN_*VR2*VI1);
    }
    else if (transType_ == PS)
    {
      dFdx[li_VR1][VR1_Phi_Offset] -= invTurnsRatio_*(VR1*VR2*(g12*sin(Phi_) + b12*cos(Phi_))
						     -VR1*VI2*(b12*sin(Phi_) - g12*cos(Phi_))
						     +VI1*VR2*(b12*sin(Phi_) - g12*cos(Phi_))
						     -VI1*VI2*(b12*sin(Phi_) + g12*cos(Phi_)));
      dFdx[li_VR2][VR2_Phi_Offset] -= invTurnsRatio_*(VR2*VR1*(g21*sin(Phi_) - b21*cos(Phi_))
						     -VR2*VI1*(b21*sin(Phi_) + g21*cos(Phi_))
						     +VI2*VR1*(b21*sin(Phi_) + g21*cos(Phi_))
						     -VI2*VI1*(b21*sin(Phi_) - g21*cos(Phi_)));
      dFdx[li_VI1][VI1_Phi_Offset] -= invTurnsRatio_*(VI1*VR2*(g12*sin(Phi_) + b12*cos(Phi_))
						     -VI1*VI2*(b12*sin(Phi_) - g12*cos(Phi_))
						     -VR1*VR2*(b12*sin(Phi_) - g12*cos(Phi_))
                                                     +VR1*VI2*(b12*sin(Phi_) + g12*cos(Phi_)));
      dFdx[li_VI2][VI2_Phi_Offset] -= invTurnsRatio_*(VI2*VR1*(g21*sin(Phi_) - b21*cos(Phi_))
						     -VI2*VI1*(b21*sin(Phi_) + g21*cos(Phi_))
						     -VR2*VR1*(b21*sin(Phi_) + g21*cos(Phi_))
						     +VR2*VI1*(b21*sin(Phi_) - g21*cos(Phi_)));
    }
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
    dFdx[li_Theta1][Theta1_Theta1_Offset] -=  VM1*VM2*(ga12*dSin12 - ba12*dCos12);
    dFdx[li_Theta1][Theta1_Theta2_Offset] +=  VM1*VM2*(ga12*dSin12 - ba12*dCos12);	 
    dFdx[li_Theta2][Theta2_Theta1_Offset] +=  VM1*VM2*(ga21*dSin21 - ba21*dCos21);
    dFdx[li_Theta2][Theta2_Theta2_Offset] -=  VM1*VM2*(ga21*dSin21 - ba21*dCos21);

    // dP / dV terms
    dFdx[li_Theta1][Theta1_VM1_Offset] += 2*VM1*ga11 + VM2*(ga12*dCos12 + ba12*dSin12);  
    dFdx[li_Theta1][Theta1_VM2_Offset] += VM1*(ga12*dCos12 + ba12*dSin12);
    dFdx[li_Theta2][Theta2_VM1_Offset] += VM2*(ga21*dCos21 + ba21*dSin21); 
    dFdx[li_Theta2][Theta2_VM2_Offset] += 2*VM2*ga22 + VM1*(ga21*dCos21 + ba21*dSin21);  

    // dQ / dTheta terms
    dFdx[li_VM1][VM1_Theta1_Offset] +=  VM1*VM2*(ga12*dCos12 + ba12*dSin12);
    dFdx[li_VM1][VM1_Theta2_Offset] -=  VM1*VM2*(ga12*dCos12 + ba12*dSin12);
    dFdx[li_VM2][VM2_Theta1_Offset] -=  VM1*VM2*(ga21*dCos21 + ba21*dSin21);
    dFdx[li_VM2][VM2_Theta2_Offset] +=  VM1*VM2*(ga21*dCos21 + ba21*dSin21);

    // dQ / dV terms
    dFdx[li_VM1][VM1_VM1_Offset] += -2*VM1*ba11 + VM2*(ga12*dSin12 - ba12*dCos12);  
    dFdx[li_VM1][VM1_VM2_Offset] += VM1*(ga12*dSin12 - ba12*dCos12);
    dFdx[li_VM2][VM2_VM1_Offset] += VM2*(ga21*dSin21 - ba21*dCos21);
    dFdx[li_VM2][VM2_VM2_Offset] += -2*VM2*ba22 + VM1*(ga21*dSin21 - ba21*dCos21);

    if (transType_ == VT)
    { 
      dFdx[li_Theta1][Theta1_N_Offset] -= ( 2*invN_*ga11*VM1*VM1 
				      + ga12*dCos12*invN_*VM1*VM2 
                                      + ba12*dSin12*invN_*VM1*VM2);

      dFdx[li_Theta2][Theta2_N_Offset] -= ( ga21*dCos21*invN_*VM2*VM1 
					    + ba21*dSin21*invN_*VM2*VM1);                             

      dFdx[li_VM1][VM1_N_Offset] -= ( -2*invN_*ba11*VM1*VM1
	                              + ga12*dSin12*invN_*VM1*VM2 
                                      - ba12*dCos12*invN_*VM1*VM2);
 
      dFdx[li_VM2][VM2_N_Offset] -= ( ga21*dSin21*invN_*VM2*VM1 
                                       - ba21*dCos21*invN_*VM2*VM1); 
    }
    else if (transType_ == PS)
    { 
      dFdx[li_Theta1][Theta1_Phi_Offset] += ( ga12*sin(Theta1-Theta2-Phi_)*VM1*VM2 
					      - ba12*cos(Theta1-Theta2-Phi_)*VM1*VM2);

      dFdx[li_Theta2][Theta2_Phi_Offset] -= ( ga21*dSin21*VM2*VM1 
					    - ba21*dCos21*VM2*VM1);                             

      dFdx[li_VM1][VM1_Phi_Offset] += ( ga12*cos(Theta1-Theta2-Phi_)*VM1*VM2 
				      + ba12*sin(Theta1-Theta2-Phi_)*VM1*VM2);
 
      dFdx[li_VM2][VM2_Phi_Offset] += ( ga21*dCos21*VM2*VM1 
                                       - ba21*dSin21*VM2*VM1); 
    }
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
// Purpose       : Loads the Q-vector contributions.  This is no-op for this
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
  if (deviceMap.empty() || (deviceMap.find("POWERGRIDTRANSFORMER")!=deviceMap.end())
                        || (deviceMap.find("PGTR")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("PowerGridTransformer", 1)
      .registerDevice("PGTR", 1)
      .registerModelType("PowerGridTransformer", 1);
  }
}

} // namespace PowerGridTransformer
} // namespace Device
} // namespace Xyce

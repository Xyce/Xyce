//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
//                  steady state power flow in a transmission grid.  This is
//                  a "monolithic model" for the entire power grid.  It is
//                  being replaced by individual device models (branch, bus
//                  shunt, load, sources) that parallelize better.
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
#include <fstream>

// ----------   Xyce Includes   ----------
#include <N_DEV_PowerGrid.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Message.h>
#include <N_DEV_SolverState.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

#include <N_UTL_FeatureTest.h>
#include <N_UTL_IndentStreamBuf.h>

namespace Xyce {
namespace Device {


namespace PowerGrid {

// enables debug output for just the digital devices
//static const int DEBUG_DEVICE = 1;

void Traits::loadInstanceParameters(ParametricData<PowerGrid::Instance> &p)
{
  p.addPar("NB", 1, &PowerGrid::Instance::numBus_)
    .setUnit(U_NONE)
    .setDescription("Number of Buses");

  p.addPar("AT", std::string("IV"), &PowerGrid::Instance::analysisType_)
    .setUnit(U_NONE)
    .setDescription("Analysis Type");

  p.addPar("BUF", std::string("busFile"), &PowerGrid::Instance::busFileName_)
    .setUnit(U_NONE)
    .setDescription("IC File Name");

  p.addPar("BRF", std::string("branchFile"), &PowerGrid::Instance::branchFileName_)
    .setUnit(U_NONE)
    .setDescription("BD File Name");
}

void Traits::loadModelParameters(ParametricData<PowerGrid::Model> &p)
{}



std::vector< std::vector<int> > Instance::jacStamp;

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &               IB,
  Model & Miter,
  const FactoryBlock &  factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    numBus_(1),
    analysisType_(""),
    busFileName_(""),
    branchFileName_("")
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl << "In device constructor for " << getName() << std::endl;
  }

  std::vector<Param>::const_iterator iter = IB.params.begin();
  std::vector<Param>::const_iterator last = IB.params.end();

  for ( ; iter != last; ++iter)
  {
    const std::string & tmpname = iter->tag();
    const bool & tmpgiven = iter->given();

    // copy instance parameter values into private variables
    if (tmpname == "NB" && tmpgiven == true)
    {
      numBus_ = iter->getImmutableValue<int>();
    }
    if (tmpname == "AT" && tmpgiven == true)
    {
       analysisType_ = iter->stringValue();
    }
    if (tmpname == "BUF" && tmpgiven == true)
    {
      busFileName_ = iter->stringValue();
    }
    if (tmpname == "BRF" && tmpgiven == true)
    {
      branchFileName_ =  iter->stringValue();
    }
  }

  // read input data
  loadBusdata();
  loadBranchData();

  // do sanity checking on input files and entered values
  if (numBus_ < 1)
  {
    UserError(*this) << "NB Instance Parameter must be an integer >=1";
  }

  if (magICmap_.size() != numBus_)
  {
     UserError(*this) << "Incorrect number of Initial Conditions for device "
                      << getName() << std::endl << "Device has NB=" << numBus_
		      << " and " << magICmap_.size() << " initial conditions in file "
                      << busFileName_;
    return;
  }

  // Build the admittance matrix from the branch data
  buildYMatrixMap();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_DUMP_VECTORS))
  {
    printYMatrixMap();
  }

  numIntVars   = 0;
  numExtVars   = 4;  // In general, this will be 4*numBus_
  numStateVars = 0;

  // temporarily make this into a 2-port device, between inputs 1 and 2
  // this is just to get the rest of the "Xyce plumbing" working

  // this code made the device into a 2-port resistor, as an intermediate step
  // set G to real part of G11, just for initial testing
  //twodComplexMap::iterator YMatrixMapIter=yMatrixMap_.begin();
  //std::complex<double> cVal = YMatrixMapIter->second;
  //G = std::abs(cVal);

  // this code needs to be generalized
  g11 = yMatrixMap_[std::make_pair(1,1)].real();
  b11 = yMatrixMap_[std::make_pair(1,1)].imag();
  g12 = yMatrixMap_[std::make_pair(1,2)].real();
  b12 = yMatrixMap_[std::make_pair(1,2)].imag();
  g21 = yMatrixMap_[std::make_pair(2,1)].real();
  b21 = yMatrixMap_[std::make_pair(2,1)].imag();
  g22 = yMatrixMap_[std::make_pair(2,2)].real();
  b22 = yMatrixMap_[std::make_pair(2,2)].imag();

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

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:

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

  // in general, there will be 4 I/0 nodes for each branch
  // these if-else blocks will be made into separate classes,
  // once they work
  if ( analysisType_ == "IV" )
  {
    li_VR1 = extLIDVec[0];
    li_VR2 = extLIDVec[1];
    li_VI1 = extLIDVec[2];
    li_VI2 = extLIDVec[3];

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
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
  else if (analysisType_ == "PQ" )
  {
    li_Theta1 = extLIDVec[0];
    li_Theta2 = extLIDVec[1];
    li_VM1 = extLIDVec[2];
    li_VM2 = extLIDVec[3];

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
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
    UserError(*this) << "Analysis Type must be IV or PQ in power grid device: " << getName();
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
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/27/01
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  // These will be split out into separate classes once they both work
  if (analysisType_ == "IV")
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
  }
  else if (analysisType_ == "PQ")
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
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one diode instance
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  double * solVec = extData.nextSolVectorRawPtr;

  // These will be split out into separate classes once they work
  if (analysisType_ == "IV")
  {
    IR1 = g11*solVec[li_VR1] +  g12*solVec[li_VR2] - b11*solVec[li_VI1] - b12*solVec[li_VI2];
    IR2 = g21*solVec[li_VR1] +  g22*solVec[li_VR2] - b21*solVec[li_VI1] - b22*solVec[li_VI2];
    II1 = b11*solVec[li_VR1] +  b12*solVec[li_VR2] + g11*solVec[li_VI1] + g12*solVec[li_VI2];
    II2 = b21*solVec[li_VR1] +  b22*solVec[li_VR2] + g21*solVec[li_VI1] + g22*solVec[li_VI2];
  }
  else if (analysisType_ == "PQ")
  {
    double dSin12 = sin(solVec[li_Theta1]-solVec[li_Theta2]);
    double dSin21 = sin(solVec[li_Theta2]-solVec[li_Theta1]);
    double dCos12 = cos(solVec[li_Theta1]-solVec[li_Theta2]);
    double dCos21 = cos(solVec[li_Theta2]-solVec[li_Theta1]);

    P1 = g11*solVec[li_VM1]*solVec[li_VM1] + solVec[li_VM1]*solVec[li_VM2]*(g12*dCos12 + b12*dSin12);
    P2 = g22*solVec[li_VM2]*solVec[li_VM2] + solVec[li_VM2]*solVec[li_VM1]*(g21*dCos21 + b21*dSin21);
    Q1 = -1*b11*solVec[li_VM1]*solVec[li_VM1] + solVec[li_VM1]*solVec[li_VM2]*(g12*dSin12 - b12*dCos12);
    Q2 = -1*b22*solVec[li_VM2]*solVec[li_VM2] + solVec[li_VM2]*solVec[li_VM1]*(g21*dSin21 - b21*dCos21);
  }
  else
  {
    UserError(*this) << "Analysis Type must be IV or PQ in power grid device: " << getName();
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
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::updatePrimaryState" <<std::endl;
  }

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

  // this needs to be vectorized in the general case.
  // These will be split out into separate classes once they work
  if (analysisType_ == "IV")
  {
    fVec[li_VR1] += IR1;
    fVec[li_VR2] += IR2;
    fVec[li_VI1] += II1;
    fVec[li_VI2] += II2;
  }
  else if (analysisType_ == "PQ")
  {
    fVec[li_Theta1] += P1;
    fVec[li_Theta2] += P2;
    fVec[li_VM1] += Q1;
    fVec[li_VM2] += Q2;
  }
  else
  {
    UserError(*this) << "Analysis Type must be IV or PQ in power grid device: " << getName();
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
// Purpose       : Loads the F-vector contributions
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

  if (analysisType_ == "IV")
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
  else if (analysisType_ == "PQ")
  {
    double * solVec = extData.nextSolVectorRawPtr;
    double dSin12 = sin(solVec[li_Theta1]-solVec[li_Theta2]);
    double dSin21 = sin(solVec[li_Theta2]-solVec[li_Theta1]);
    double dCos12 = cos(solVec[li_Theta1]-solVec[li_Theta2]);
    double dCos21 = cos(solVec[li_Theta2]-solVec[li_Theta1]);

    // dP / dTheta terms
    dFdx[li_Theta1][Theta1_Theta1_Offset] -=  solVec[li_VM1]*solVec[li_VM2]*(g12*dSin12 - b12*dCos12);
    dFdx[li_Theta1][Theta1_Theta2_Offset] +=  solVec[li_VM1]*solVec[li_VM2]*(g12*dSin12 - b12*dCos12);
    dFdx[li_Theta2][Theta2_Theta1_Offset] +=  solVec[li_VM1]*solVec[li_VM2]*(g21*dSin21 - b21*dCos21);
    dFdx[li_Theta2][Theta2_Theta2_Offset] -=  solVec[li_VM1]*solVec[li_VM2]*(g21*dSin21 - b21*dCos21);

    // dP / dV terms
    dFdx[li_Theta1][Theta1_VM1_Offset] += 2*solVec[li_VM1]*g11 + solVec[li_VM2]*(g12*dCos12 + b12*dSin12);
    dFdx[li_Theta1][Theta1_VM2_Offset] += solVec[li_VM2]*(g21*dCos21 + b21*dSin21);
    dFdx[li_Theta2][Theta2_VM1_Offset] += solVec[li_VM1]*(g12*dCos12 + b12*dSin12);
    dFdx[li_Theta2][Theta2_VM2_Offset] += 2*solVec[li_VM2]*g22 + solVec[li_VM1]*(g21*dCos21 + b21*dSin21);

    // dQ / dTheta terms
    dFdx[li_VM1][VM1_Theta1_Offset] -=  solVec[li_VM1]*solVec[li_VM2]*(g12*dCos12 + b12*dSin12);
    dFdx[li_VM1][VM1_Theta2_Offset] +=  solVec[li_VM1]*solVec[li_VM2]*(g12*dCos12 + b12*dSin12);
    dFdx[li_VM2][VM2_Theta1_Offset] +=  solVec[li_VM1]*solVec[li_VM2]*(g21*dCos21 + b21*dSin21);
    dFdx[li_VM2][VM2_Theta2_Offset] -=  solVec[li_VM1]*solVec[li_VM2]*(g21*dCos21 + b21*dSin21);

    // dQ / dV terms
    dFdx[li_VM1][VM1_VM1_Offset] += -2*solVec[li_VM1]*b11 + solVec[li_VM2]*(g12*dSin12 - b12*dCos12);
    dFdx[li_VM1][VM1_VM2_Offset] += solVec[li_VM2]*(g21*dSin21 - b21*dCos21);
    dFdx[li_VM2][VM2_VM1_Offset] += solVec[li_VM1]*(g21*dSin12 - b21*dCos12);
    dFdx[li_VM2][VM2_VM2_Offset] += -2*solVec[li_VM2]*b22 + solVec[li_VM1]*(g21*dSin21 - b21*dCos21);
  }
  else
  {
    UserError(*this) << "Analysis Type must be IV or PQ in power grid device: " << getName();
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
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadBusdata ()
//
// Purpose       : Loads IC and Bus data from a file.  This function will likely
//                 go away, once the monolithic PowerGrid model does.
//                 In that case, it may be more convenient to load
//                 the initial conditions for each bus from IC=
//                 statements on the instance lines
//
// Special Notes :
//
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
bool Instance::loadBusdata ()
{
  // If the file name is enclosed in double quotes, strip them.
  if (busFileName_[0] == '"' && busFileName_[busFileName_.length()-1] =='"')
  {
    busFileName_ = busFileName_.substr(1, busFileName_.length()-2);
  }

  // make the file exists
  std::ifstream busDataIn;
  busDataIn.open(busFileName_.c_str(), std::ios::in);
  if ( !busDataIn.is_open() )
  {
    Report::UserError() << "Could not find file " << busFileName_;
    return false;
  }

  // read in data for initial conditions for voltage magnitude and
  // phase at each bus
  int bus;
  while ( busDataIn >> bus )
  {
    if ( magICmap_.count(bus) == 0 )
    {
      if (busDataIn >> magICmap_[bus] >>  angleICmap_[bus]
	            >> busShuntConductance_[bus] >> busShuntSusceptance_[bus])
      {
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_DUMP_VECTORS))
	{
          Xyce::dout() << "Read line in file " << busFileName_ << " for bus "
	               << bus << ". (Mag,angle) ICs=(" << magICmap_[bus] << ","
                       << angleICmap_[bus] << "). Shunt (G,B)=("
                       << busShuntConductance_[bus] << ","
	               << busShuntSusceptance_[bus] << ")" << std::endl;
        }
      }
      else
      {
        Report::UserError() << "Problem reading file " << busFileName_
	     << std::endl << "File must use space separated values,"
             << "with 3 items per line.";
        busDataIn.close();
        return false;
      }
    }
    else
    {
      Report::UserError() << "Problem reading " << busFileName_ << std::endl
	      << "There is a duplicate entry for bus number " << bus;
      busDataIn.close();
      return false;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadBranchData ()
//
// Purpose       : Loads branch data from a file.  This function will likely
//                 go away, once the monolithic PowerGrid model does.
//                 In that case, it may be more convenient to load
//                 the branch data from instance parameters on each
//                 branch's instance line
//
// Special Notes :
//
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/14
//-----------------------------------------------------------------------------
bool Instance::loadBranchData ()
{
  // If the file name is enclosed in double quotes, strip them.
  if (branchFileName_[0] == '"' && branchFileName_[branchFileName_.length()-1] =='"')
  {
    branchFileName_ = branchFileName_.substr(1, branchFileName_.length()-2);
  }

  std::ifstream branchDataIn;
  branchDataIn.open(branchFileName_.c_str(), std::ios::in);
  if ( !branchDataIn.is_open() )
  {
    Report::UserError() << "Could not find file " << branchFileName_;
    return false;
  }

  // read in data for branch data
  int bus1, bus2;
  double resistance, reactance, susceptance,turnsRatio;
  twodKey fromToKey, toFromKey;
  while ( branchDataIn >> bus1 >> bus2)
  {
    if ( (bus1 < 1) || (bus1 > numBus_) || (bus2 < 1) || (bus2 > numBus_) )
    {
      Report::UserError() << "Invalid bus number in file " << branchFileName_ << std::endl
			  << "Bus numbers should be between 1 and " << numBus_;
      return false;
    }

    // key1 is to-from data from file.  key2 is from-to version.
    // input data should only have one of them for each possible branch
    fromToKey = std::make_pair(bus1,bus2);
    toFromKey = std::make_pair(bus2,bus1);

    // insert data if neither key is in the map
    if (( branchResistance_.count(fromToKey) == 0 ) &&
               ( branchResistance_.count(toFromKey) == 0 ))
    {
      if ( branchDataIn >> resistance >> reactance >> susceptance >> turnsRatio )
      {
        branchResistance_[fromToKey]=resistance;
	branchReactance_[fromToKey]=reactance;
	branchSusceptance_[fromToKey]=susceptance;
        // IEEE CDF format uses 0 to denote no transformer is present.
        // This is effectively a turns-ratio of 1.0
        turnsRatio_[fromToKey]= (turnsRatio > 0) ? turnsRatio : 1.0;
	//branchResistance_.insert(twodMap::value_type(key1,resistance));

        if (DEBUG_DEVICE && isActive(Diag::DEVICE_DUMP_VECTORS))
	{
          Xyce::dout() << "Read line in file " << branchFileName_ << " for buses "
	      << bus1 << " " << bus2 << ". (R,X,B)=(" << branchResistance_[fromToKey]
	      << "," << branchReactance_[fromToKey] << "," << branchSusceptance_[fromToKey]
	      << "). Turns Ratio=" << turnsRatio_[fromToKey] << std::endl;
        }
      }
      else
      {
        Report::UserError() << "Problem reading file " << branchFileName_
             << std::endl << "File must use space separated values,"
             << "with 3 items per line.";
        branchDataIn.close();
       return false;
      }
    }
    else
    {
       Report::UserError() << "Problem reading " << branchFileName_ << std::endl
	      << "There is a duplicate branch data entry for bus numbers "
			   << bus1 << " and " << bus2 << std::endl;
      branchDataIn.close();
      return false;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::buildYMatrixMap ()
//
// Purpose       : Builds the Y Matrix from the branch data.
//
// Special Notes :
//
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/15/14
//-----------------------------------------------------------------------------
bool Instance::buildYMatrixMap()
{
  // may combine these maps into one data structure
  twodMap::iterator iterRes;
  twodMap::iterator iterReact=branchReactance_.begin();
  twodMap::iterator iterSuscept=branchSusceptance_.begin();
  twodKey fromToKey, toFromKey, toToKey, fromFromKey;
  int toID,fromID;
  double invTurnsRatio; // makes equations easier to read
  std::complex<double> zVal,selfVal;

  if ( (branchReactance_.size() != branchResistance_.size()) ||
       (branchSusceptance_.size() != branchResistance_.size()) ||
       (turnsRatio_.size() != branchResistance_.size()) )
  {
    Report::UserError() << "Branch Resistance, Reactance and Susceptance Matrices not same size.";
    return false;
  }

  for (iterRes=branchResistance_.begin();iterRes!=branchResistance_.end();++iterRes)
  {
    fromToKey=iterRes->first;
    fromID=fromToKey.first;
    toID=fromToKey.second;
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_DUMP_VECTORS))
    {
      Xyce::dout() << "Processing Y Matrix for branch ID: " << fromID << "-" << toID << std::endl;
    }

    toFromKey = std::make_pair(toID,fromID);
    fromFromKey = std::make_pair(fromID,fromID);;
    toToKey = std::make_pair(toID,toID);
    // makes equations easier to read
    invTurnsRatio = (1./turnsRatio_[fromToKey]);

    // zVal is the impedance between the FROM and TO buses
    // (1/zVal)*invTurnsRatio is the Y[i,j] and Y[j,i] term in the Admittance Matrix.
    //zVal.real(branchResistance_[fromToKey]);
    //zVal.imag(branchReactance_[fromToKey]);
    zVal=std::complex<double>(branchResistance_[fromToKey],branchReactance_[fromToKey]);
    yMatrixMap_[fromToKey] = (-1./zVal) * invTurnsRatio;
    yMatrixMap_[toFromKey] = (-1./zVal) * invTurnsRatio;

    // addition to the Yii term is then Yij + half of the branch susceptance defined
    // in the IEEE CDF Format, since CDF gives the parameters for a PI model of
    // the branch.  Y[j,j] term is unaffected by turns ratio.  Y[i,i] term is.
    selfVal = (1./zVal);
    //selfVal.imag(selfVal.imag() + 0.5*branchSusceptance_[fromToKey]);
    selfVal=std::complex<double>(selfVal.real(),selfVal.imag() + 0.5*branchSusceptance_[fromToKey]);
    yMatrixMap_[fromFromKey] += selfVal*invTurnsRatio*invTurnsRatio;
    yMatrixMap_[toToKey] += selfVal;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_DUMP_VECTORS))
    {
      Xyce::dout() << "zVal, selfVal and inverse turns Ratio are: " << zVal << " " << selfVal
         << " " << invTurnsRatio << std::endl << "Updated Y Matrix values are: " << std::endl
         << "  Index " << fromID << " and " << toID << " = " << yMatrixMap_[fromToKey] << std::endl
         << "  Index " << toID << " and " << fromID << " = " << yMatrixMap_[toFromKey] << std::endl
         << "  Index " << fromID << " and " << fromID << " = " << yMatrixMap_[fromFromKey] << std::endl
	 << "  Index " << toID << " and " << toID << " = " << yMatrixMap_[toToKey] << std::endl;
    }
  }

  // now add in bus conductance and susceptance
  for (fromID=1; fromID <= numBus_; fromID++)
  {
    fromFromKey = std::make_pair(fromID,fromID);
    yMatrixMap_[fromFromKey]=std::complex<double>(yMatrixMap_[fromFromKey].real() +  busShuntConductance_[fromID],
                                      yMatrixMap_[fromFromKey].imag() + busShuntSusceptance_[fromID]);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::printYMatrixMap ()
//
// Purpose       : prints the Y Matrix Map
//
// Special Notes :
//
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/16/14
//-----------------------------------------------------------------------------
bool Instance::printYMatrixMap()
{
  twodComplexMap::iterator iter;
  twodKey key;
  int toID, fromID;
  std::complex<double> cVal;

  std::cout << "Y Matrix complex values are: " << std::endl;
  for (iter=yMatrixMap_.begin();iter!=yMatrixMap_.end();++iter)
  {
    key=iter->first;
    toID=key.first;
    fromID=key.second;
    cVal= iter->second;

    Xyce::dout() << "    (" << toID << "," << fromID << ") = " << cVal << std::endl;
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
void Model::forEachInstance(DeviceInstanceOp &op) const
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
  if (deviceMap.empty() || (deviceMap.find("POWERGRID")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("PowerGrid", 1)
      .registerModelType("PowerGrid", 1);
  }
}

} // namespace PowerGrid
} // namespace Device
} // namespace Xyce

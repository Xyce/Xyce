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

//-------------------------------------------------------------------------
//
// Purpose        : Implement lossless transmission line
//
// Special Notes  :
//
// Creator        : Tom Russo, SNL, Component Information and Models
//
// Creation Date  : 06/14/01
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceState.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_TRA.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

#include <N_UTL_BreakPoint.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Functors.h>
#include <N_UTL_MachDepParams.h>
#include <N_UTL_AssemblyTypes.h>

namespace Xyce {
namespace Device {
namespace TRA {

void Traits::loadInstanceParameters(ParametricData<TRA::Instance> &p)
{
  p.addPar("Z0", 0.0, &TRA::Instance::Z0)
    .setUnit(U_OHM)
    .setDescription("Characteristic Impedance");

  p.addPar("ZO", 0.0, &TRA::Instance::ZO)
    .setUnit(U_OHM)
    .setDescription("Characteristic Impedance");

  p.addPar("TD", 0.0, &TRA::Instance::td)
    .setUnit(U_SECOND)
    .setDescription("Time delay");

  p.addPar("F", 0.0, &TRA::Instance::freq)
    .setUnit(U_HZ)
    .setDescription("Frequency");

  p.addPar("NL", 0.25, &TRA::Instance::NL)
    .setDescription("Length in wavelengths");
}

void Traits::loadModelParameters(ParametricData<TRA::Model> &p)
{}

std::vector< std::vector<int> > Instance::jacStamp;

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 6/15/01
//-----------------------------------------------------------------------------

Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & instance_block,
  Model &               model,
  const FactoryBlock &  factory_block)
  : DeviceInstance(instance_block, configuration.getInstanceParameters(), factory_block),
    model_(model),
    Z0(0.0),
    G0(0.0),
    td(0.0),
    freq(0.0),
    NL(0.25),
    li_Pos1(-1),
    li_Neg1(-1),
    li_Int1(-1),
    li_Ibr1(-1),
    li_Pos2(-1),
    li_Neg2(-1),
    li_Int2(-1),
    li_Ibr2(-1),
    li_branch_data_1(-1),
    li_branch_data_2(-1),
    APos1EquPos1NodeOffset(-1),
    APos1EquInt1NodeOffset(-1),
    AInt1EquPos1NodeOffset(-1),
    AInt1EquInt1NodeOffset(-1),
    AInt1EquIbr1NodeOffset(-1),
    ANeg1EquIbr1NodeOffset(-1),
    AIbr1EquInt1NodeOffset(-1),
    AIbr1EquNeg1NodeOffset(-1),
    APos2EquPos2NodeOffset(-1),
    APos2EquInt2NodeOffset(-1),
    AInt2EquPos2NodeOffset(-1),
    AInt2EquInt2NodeOffset(-1),
    AInt2EquIbr2NodeOffset(-1),
    ANeg2EquIbr2NodeOffset(-1),
    AIbr2EquInt2NodeOffset(-1),
    AIbr2EquNeg2NodeOffset(-1),
    AIbr1EquPos2NodeOffset(-1),
    AIbr1EquNeg2NodeOffset(-1),
    AIbr1EquIbr2NodeOffset(-1),
    AIbr2EquPos1NodeOffset(-1),
    AIbr2EquNeg1NodeOffset(-1),
    AIbr2EquIbr1NodeOffset(-1),
    last_t(0.0),
    v1(0.0),
    v2(0.0),
    first_BP_call_done(false),
    timeOld(-1.0),
    newBreakPoint(false),
    newBreakPointTime(0.0)
{
  numIntVars   = 4;
  numExtVars   = 4;
  numStateVars = 0;
  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 2;    // this is the space to allocate if lead current or power is needed.

  devConMap.resize(4);
  devConMap[0] = 1;
  devConMap[1] = 1;
  devConMap[2] = 2;
  devConMap[3] = 2;

  if( jacStamp.empty() )
  {
    jacStamp.resize(8);
    jacStamp[0].resize(2);
    jacStamp[0][0]=0;
    jacStamp[0][1]=4;
    jacStamp[1].resize(1);
    jacStamp[1][0]=5;
    jacStamp[2].resize(2);
    jacStamp[2][0]=2;
    jacStamp[2][1]=6;
    jacStamp[3].resize(1);
    jacStamp[3][0]=7;
    jacStamp[4].resize(3);
    jacStamp[4][0]=0;
    jacStamp[4][1]=4;
    jacStamp[4][2]=5;
    jacStamp[5].resize(5);
    jacStamp[5][0]=1;
    jacStamp[5][1]=2;
    jacStamp[5][2]=3;
    jacStamp[5][3]=4;
    jacStamp[5][4]=7;
    jacStamp[6].resize(3);
    jacStamp[6][0]=2;
    jacStamp[6][1]=6;
    jacStamp[6][2]=7;
    jacStamp[7].resize(5);
    jacStamp[7][0]=0;
    jacStamp[7][1]=1;
    jacStamp[7][2]=3;
    jacStamp[7][3]=5;
    jacStamp[7][4]=6;
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (instance_block.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << " Z0 = " << Z0 << std::endl;
    Xyce::dout() << " td = " << td << std::endl;
    Xyce::dout() << " freq = " << freq << std::endl;
    Xyce::dout() << " NL = " << NL << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 7/02/04
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  bool bsuccess = true;

  if (!given("Z0"))
  {
    if (given("ZO"))
      Z0 = ZO;
    else
    {
      UserError(*this) << "Z0 not given.";
      bsuccess = false;
    }
  }
  if (Z0>0)
  {
    G0 = 1.0/Z0;
  }
  else
  {
    UserError(*this) << "Invalid (zero or negative) impedance (Z0) given.";
    bsuccess = false;
  }

  // Must give either TD or F.
  if (!given("TD") && !given("F"))
  {
    UserError(*this) << "Neither time delay (TD) nor frequency (F) given.";
    bsuccess = false;
  }
  if (given("TD") && given("F"))
  {
    UserError(*this) << "Both time delay (TD) and frequency (F) given.  Pick one.";
    bsuccess = false;
  }

  if ( !given("TD") )
  {
    if (freq <= 0)
    {
      UserError(*this) << "Invalid (zero or negative) frequency (F) given.";
      bsuccess = false;
    }
    else if (NL <= 0)
    {
      UserError(*this) << "Invalid (zero or negative) NL parameter given.";
      bsuccess = false;
    }
    else
    {
      td = NL/freq;
      if (td <= 0)
      {
        UserError(*this) << "Zero or negative time delay.";
        bsuccess = false;
      }
    }
  }
  else
  {
    if (td <= 0)
    {
      UserError(*this) << "Zero or negative time delay.";
      bsuccess = false;
    }
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       : function for registering, and setting up, local ID's.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                                      const std::vector<int> & extLIDVecRef )
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "In Instance::registerLIDs\n\n";
    Xyce::dout() << "name             = " << getName() << std::endl;
    Xyce::dout() << "number of internal variables: " << numIntVars << std::endl;
    Xyce::dout() << "number of external variables: " << numExtVars << std::endl;
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the linear algebra
  // entities.  This assumes an order.

  li_Pos1 = extLIDVec[0];
  li_Neg1 = extLIDVec[1];
  li_Pos2 = extLIDVec[2];
  li_Neg2 = extLIDVec[3];

  // Now do the internal variables

  li_Int1 = intLIDVec[0];
  li_Ibr1 = intLIDVec[1];
  li_Int2 = intLIDVec[2];
  li_Ibr2 = intLIDVec[3];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << " VARIABLE Indicies " << std::endl;
    Xyce::dout() << "li_Pos1 = " << li_Pos1 << std::endl;
    Xyce::dout() << "li_Neg1 = " << li_Neg1 << std::endl;
    Xyce::dout() << "li_Int1 = " << li_Int1 << std::endl;
    Xyce::dout() << "li_Ibr1 = " << li_Ibr1 << std::endl;
    Xyce::dout() << "li_Pos2 = " << li_Pos2 << std::endl;
    Xyce::dout() << "li_Neg2 = " << li_Neg2 << std::endl;
    Xyce::dout() << "li_Int2 = " << li_Int2 << std::endl;
    Xyce::dout() << "li_Ibr2 = " << li_Ibr2 << std::endl;

    Xyce::dout() << "\nEnd of Instance::register LIDs\n";
    Xyce::dout() << section_divider << std::endl;
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
  addInternalNode(symbol_table,  li_Int1, getName(), "int1");
  addInternalNode(symbol_table,  li_Int2, getName(), "int2");
  addInternalNode(symbol_table,  li_Ibr1, getName(), "i1");
  addInternalNode(symbol_table,  li_Ibr2, getName(), "i2");

  if (loadLeadCurrent)
  {
    addBranchDataNode(symbol_table, li_branch_data_1, getName(), "BRANCH_D1");
    addBranchDataNode(symbol_table, li_branch_data_2, getName(), "BRANCH_D2");
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef)
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       : One store var for device current.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/05/2013
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef )
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());
}


//-----------------------------------------------------------------------------
// Function      : Instance::registerBranchDataLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 12/21/15
//-----------------------------------------------------------------------------
void Instance::registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef)
{
  AssertLIDs(branchLIDVecRef.size() == getNumBranchDataVars());

  if (loadLeadCurrent)
  {
    li_branch_data_1 = branchLIDVecRef[0];
    li_branch_data_2 = branchLIDVecRef[1];
  }
}


//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 9/2/02
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
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 9/2/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  APos1EquPos1NodeOffset = jacLIDVec[0][0];
  APos1EquInt1NodeOffset = jacLIDVec[0][1];

  ANeg1EquIbr1NodeOffset = jacLIDVec[1][0];

  APos2EquPos2NodeOffset = jacLIDVec[2][0];
  APos2EquInt2NodeOffset = jacLIDVec[2][1];

  ANeg2EquIbr2NodeOffset = jacLIDVec[3][0];

  AInt1EquPos1NodeOffset = jacLIDVec[4][0];
  AInt1EquInt1NodeOffset = jacLIDVec[4][1];
  AInt1EquIbr1NodeOffset = jacLIDVec[4][2];

  AIbr1EquNeg1NodeOffset = jacLIDVec[5][0];
  AIbr1EquPos2NodeOffset = jacLIDVec[5][1];
  AIbr1EquNeg2NodeOffset = jacLIDVec[5][2];
  AIbr1EquInt1NodeOffset = jacLIDVec[5][3];
  AIbr1EquIbr2NodeOffset = jacLIDVec[5][4];

  AInt2EquPos2NodeOffset = jacLIDVec[6][0];
  AInt2EquInt2NodeOffset = jacLIDVec[6][1];
  AInt2EquIbr2NodeOffset = jacLIDVec[6][2];

  AIbr2EquPos1NodeOffset = jacLIDVec[7][0];
  AIbr2EquNeg1NodeOffset = jacLIDVec[7][1];
  AIbr2EquNeg2NodeOffset = jacLIDVec[7][2];
  AIbr2EquIbr1NodeOffset = jacLIDVec[7][3];
  AIbr2EquInt2NodeOffset = jacLIDVec[7][4];

}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 TRA instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/19/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess = true;

  double coef_pos1;
  double coef_neg1;
  double coef_int1;
  double coef_ibr1;
  double coef_pos2;
  double coef_neg2;
  double coef_int2;
  double coef_ibr2;

  double * solVec = extData.nextSolVectorRawPtr;
  double * fVec = extData.daeFVectorRawPtr;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "  Instance::loadDAEFVector" << std::endl;
    Xyce::dout() << "  name = " << getName() <<std::endl;
  }

  // Most of the work has already been done by uIVB.
  coef_pos1 = (Vpos1-Vint1)*G0;
  coef_neg1 = -Ibr1;
  coef_int1 = -(Vpos1-Vint1)*G0+Ibr1;
  coef_ibr1 = ((Vint1-Vneg1)-v1);
  coef_pos2 = (Vpos2-Vint2)*G0;
  coef_neg2 = -Ibr2;
  coef_int2 = -(Vpos2-Vint2)*G0+Ibr2;
  coef_ibr2 = ((Vint2-Vneg2)-v2);


  fVec[li_Pos1] += coef_pos1;
  fVec[li_Neg1] += coef_neg1;
  fVec[li_Int1] += coef_int1;
  fVec[li_Ibr1] += coef_ibr1;
  fVec[li_Pos2] += coef_pos2;
  fVec[li_Neg2] += coef_neg2;
  fVec[li_Int2] += coef_int2;
  fVec[li_Ibr2] += coef_ibr2;

  if( loadLeadCurrent )
  {
    // The convention is that current "into" the positive terminal of each port of
    // the device has a positive sign.  Current out of the positive terminal of
    // each port of the device has a negative sign.
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    double * junctionV = extData.nextJunctionVCompRawPtr;
    leadF[li_branch_data_1] = Ibr1;
    leadF[li_branch_data_2] = Ibr2;

    junctionV[li_branch_data_1] = solVec[li_Pos1] - solVec[li_Neg1];
    junctionV[li_branch_data_2] = solVec[li_Pos2] - solVec[li_Neg2];
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single device instance.
//
// Special Notes : The F-vector is an algebraic constaint.
//
//                 The special notes below are those that were taken from
//                 the old loadAnalyticJacobian header.  The matrix
//                 it describes is the full jacobian matrix:
//----------------------------------------------------------------------
// Special Notes : This is based on there being two two-node ports
//                 and a model of the following sort:
//
//  Pos1  o-----+          +------o Pos2
//              |          |
//              \          \
//              /          /
//              \Z0        \ Z0
//              /          /
//              |          |
//              oInt1      o Int2
//              |          |
//           ++++++      ++++++
//           | V1 |      | V2 |
//           ------      ------
//              |          |
//   Neg1 o-----+          +------o Neg2
//
// There are also two branch currents, Ibr1 and Ibr2 for left and right
// sides as well.
//
//  The matrix for this ends up being:
//            V_Pos1  V_Neg1  V_Int1  Ibr1  V_Pos2  V_Neg2  V_Int2  V_Ibr2
//            ------------------------------------------------------------
//  KCL Pos1    a               b
//  KCL Neg1                           c
//  KCL Int1    d               e      f
//  KCL Ibr1            g       h             i       j               k
//  KCL Pos2                                  l                m
//  KCL Neg2                                                          n
//  KCL Int2                                  o                p      q
//  KCL Ibr2    r       s              t              u        v
//
//  When doing time integration, i,j,k,r,s and t are zero, those dependences
//  are time-delayed, i.e. the equations for the output depend on time-delayed
//  of the input.  For DC calculations, i,j,k,r,s and t are non-zero.
//  Note:  For frequency-domain, i,j,k,r,s, and t are non-zero and complex.
//  The frequency-domain case is not loaded here.
//
// The right hand sides are:
// Pos1:  (V_int1-V_Pos1)*G0            Pos2: (V_Int2-V_Pos2)*G0
// Neg1:  -Ibr1                         Neg2: -Ibr2
// Int1:  (V_Pos1-V_Int1)*G0+Ibr1       Int2: (V_Pos2-V_Int2)*G0+Ibr2
// Ibr1:  (V_Int1-V_Neg1)-V1            Ibr2: (V_Int2-V_Neg2)-V2
//
// For transient operation, v1 and v2 depend on values of voltage and
//  current at delayed time at the opposite port:
//   V1 = DeltaV2(t-td)+Z0*Ibr2(t-td)
//   V2 = DeltaV1(t-td)+Z0*Ibr1(t-td)
//
// For DC operation V1=VPos2-Vneg2+Ibr2*Z0, V2 = Vpos1-Vneg1+Ibr1*Z0
//--------------------------------------------------------------------
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/19/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "  name             = " << getName() << std::endl;
  }


  dFdx[li_Pos1][APos1EquPos1NodeOffset] += G0;

  dFdx[li_Pos1][APos1EquInt1NodeOffset] -= G0;


  dFdx[li_Int1][AInt1EquPos1NodeOffset] -= G0;

  dFdx[li_Int1][AInt1EquInt1NodeOffset] += G0;

  dFdx[li_Int1][AInt1EquIbr1NodeOffset] += 1.0;


  dFdx[li_Neg1][ANeg1EquIbr1NodeOffset] -= 1.0;


  dFdx[li_Ibr1][AIbr1EquInt1NodeOffset] += 1.0;

  dFdx[li_Ibr1][AIbr1EquNeg1NodeOffset] -= 1.0;
  if( DCMODE )
  {

    dFdx[li_Ibr1][AIbr1EquPos2NodeOffset] -= 1.0;

    dFdx[li_Ibr1][AIbr1EquNeg2NodeOffset] += 1.0;

    dFdx[li_Ibr1][AIbr1EquIbr2NodeOffset] -= Z0;
  }


  dFdx[li_Pos2][APos2EquPos2NodeOffset] += G0;

  dFdx[li_Pos2][APos2EquInt2NodeOffset] -= G0;


  dFdx[li_Int2][AInt2EquPos2NodeOffset] -= G0;

  dFdx[li_Int2][AInt2EquInt2NodeOffset] += G0;

  dFdx[li_Int2][AInt2EquIbr2NodeOffset] += 1.0;


  dFdx[li_Neg2][ANeg2EquIbr2NodeOffset] -= 1.0;


  dFdx[li_Ibr2][AIbr2EquInt2NodeOffset] += 1.0;

  dFdx[li_Ibr2][AIbr2EquNeg2NodeOffset] -= 1.0;
  if( DCMODE )
  {

    dFdx[li_Ibr2][AIbr2EquPos1NodeOffset] -= 1.0;

    dFdx[li_Ibr2][AIbr2EquNeg1NodeOffset] += 1.0;

    dFdx[li_Ibr2][AIbr2EquIbr1NodeOffset] -= Z0;
  }


  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    Xyce::dout() << subsection_divider << std::endl;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       : update primary state for one TRA instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 1/10/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState()
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << std::endl << subsection_divider << std::endl;
    Xyce::dout() << "In TRA::updatePrimaryState\n";
    Xyce::dout() << " last_t is " << last_t << std::endl;
    Xyce::dout() << " v1 is " << v1 << std::endl;
    Xyce::dout() << " v2 is " << v2 << std::endl;
  }

  return updateIntermediateVars ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one TRA instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 1/10/01
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  double * solVec = extData.nextSolVectorRawPtr;    // the current guess

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << std::endl << subsection_divider << std::endl;
    Xyce::dout() << "  In ::updateIntermediateVars\n\n";
  }

  Vpos1 = Vpos2 = Vneg1 = Vneg2 = Vint1 = Vint2 = 0.0;

  Vpos1 = solVec[li_Pos1];
  Vneg1 = solVec[li_Neg1];
  Vint1 = solVec[li_Int1];
  Ibr1  = solVec[li_Ibr1];
  Vpos2 = solVec[li_Pos2];
  Vneg2 = solVec[li_Neg2];
  Vint2 = solVec[li_Int2];
  Ibr2  = solVec[li_Ibr2];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << " Vpos1 = " << Vpos1 << std::endl;
    Xyce::dout() << " Vneg1 = " << Vneg1 << std::endl;
    Xyce::dout() << " Vint1 = " << Vint1 << std::endl;
    Xyce::dout() << " Ibr1 = " << Ibr1 << std::endl;
    Xyce::dout() << " Vpos2 = " << Vpos2 << std::endl;
    Xyce::dout() << " Vneg2 = " << Vneg2 << std::endl;
    Xyce::dout() << " Vint2 = " << Vint2 << std::endl;
    Xyce::dout() << " Ibr2 = " << Ibr2 << std::endl;
  }

  // Test if we're doing DC or Transient
  if ((getSolverState().dcopFlag))
  {
    // DC operation
    DCMODE=true;
    v1 = (Vpos2-Vneg2)+Z0*Ibr2;
    v2 = (Vpos1-Vneg1)+Z0*Ibr1;
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
      Xyce::dout() << "DC Mode, V1 = " << v1 <<  ", V2 = " << v2  << std::endl;
  }
  else
  {
    double currentTime = getSolverState().currTime_;
    DCMODE=false;
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "Not DC, newtonIter = " << getSolverState().newtonIter;
      Xyce::dout() << " Time is " << currentTime << std::endl;
    }
    // Transient operation
    // Now determine if we're on the first newton step of an iteration
    if (getSolverState().newtonIter == 0 && (currentTime != timeOld))
    {
      timeOld = currentTime;
      // we are, so need to manipulate history and calculate v1,v2.
      // If we're the first time step, we need to initialize it
      if (getSolverState().initTranFlag_)
      {
        last_t = currentTime;
        v1 = (Vpos2-Vneg2)+Z0*Ibr2;
        v2 = (Vpos1-Vneg1)+Z0*Ibr1;

        history.clear();
        history.push_back(History(-2*td,v1,v2));
        history.push_back(History(-td,v1,v2));
        history.push_back(History(0,v1,v2));
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << "Transient, first time, T = ";
          Xyce::dout() << currentTime;
          Xyce::dout() << ", V1 = " << v1;
          Xyce::dout() << ", V2 = " << v2  << std::endl;
        }
      }
      else
      {
        double delayedTime = currentTime-td;

        // now get the values of v1 and v2 from the delayed-time
        // information
        InterpV1V2FromHistory(delayedTime, &v1, &v2);
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << "  Done with interpolation to delayedTime=";
          Xyce::dout() << delayedTime;
          Xyce::dout() << ", have v1="<<v1 << " and v2=" << v2 << std::endl;
          Xyce::dout() << " INTERP " << delayedTime << " " << v1 << " " << v2 << std::endl;
          Xyce::dout() << " Set last_t to " << currentTime << std::endl;
        }
        // now save the current time so we can have it next time
        // we get to this block (i.e. on the next time step)
        last_t = currentTime;
      }
    }
    else
    {
      // we're on the second iteration or later of the second time
      // step or later.  Re-use the values of v1 and v2 from the
      // first iteration  of this time step.  We don't care  what time
      // it is
      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "second or later iteration, t is " << currentTime;
        Xyce::dout() << " have last_t = " << last_t << " v1="<<v1 << " and v2=" << v2 << std::endl;
      }
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::pruneHistory
// Purpose       : sift through the transmission line state history and
//                 delete records that are so old they'll never be used again
// Special Notes :
// Scope         : private
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 6/15/2001
//-----------------------------------------------------------------------------

void Instance::pruneHistory(double t1)
{

  // The input t is the oldest time for which we'll ever interpolate again.
  // That means we only need two times in the history that are older than t1,
  // so this routine drops everything off the head but the most recent 2 that
  // are older than t.

  std::vector<History>::iterator first = history.begin();
  std::vector<History>::iterator it1;
  std::vector<History>::iterator last = history.end();
  int i;

  last--; // point to last stored item, not end of the list!
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() << "Pruning for time t1="<<t1 << std::endl;
    Xyce::dout() << " Oldest in list is t="<<first->t<<" v1 = "<<first->v1 <<
      " v2="<<first->v2 << std::endl;
    Xyce::dout() << " latest in list is t="<<last->t<<" v1 = "<<last->v1
         << " v2="<<last->v2 << std::endl;
  }
  // First find the first element for which the stored time is greater than t
  for (it1 = first, i = 0; it1->t < t1 && it1 != last; ++it1, ++i)
  {
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "i = " << i << " t = " << it1->t;
      Xyce::dout() << " v1 = " << it1->v1;
      Xyce::dout() << " v2 = " << it1->v2 << std::endl;
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "   i ="<<i << std::endl;
  }

  //   Now it1 points to the first element with t>t1
  if (i > 2)
    {
      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "Need to prune.  Keeping " << it1->t << std::endl;
      }
      // if i>2 we have  too many old ones.
      // back up 2
      it1--;
      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "                Keeping " << it1->t << std::endl;
      }

      it1--;
      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "                Keeping " << it1->t << std::endl;
      }
      // delete everything from the first to it1, not counting it1
      history.erase(first,it1);
    }
  // otherwise we don't need to do anything.
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
  }
}


//-----------------------------------------------------------------------------
// Function      : Instance::InterpV1V2FromHistory
// Purpose       : Use 3-point lagrange interpolation to determine
//                 v1(t) and v2(t) at a specified time in the past
// Special Notes :
// Scope         : private
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 6/15/2001
//-----------------------------------------------------------------------------
void Instance::InterpV1V2FromHistory(double t, double * v1p,
                                              double *v2p)
{
  std::vector<History>::iterator first = history.begin();
  std::vector<History>::iterator it1;
  std::vector<History>::iterator last = history.end();
  double t1,t2,t3;
  double dt1,dt2,dt3;
  double v11,v21,v12,v22,v13,v23;
  double dt12,dt13,dt23;
  double f1,f2,f3;    // interpolating functions

  if (history.size() <= 0)
  {
    DevelFatal(*this).in("Instance::InterpV1V2FromHistory")
      << " InterpV1V2FromHistory called but history list is"
      << " empty.";
  }

  last--; // point to the last stored item, not the tail of the list!
  // sanity clause (you canna foola me, I know they're ain'ta no sanity
  // clause!)
  //  if (t < first->t || t > last->t)
  if (t - first->t < -Util::MachineDependentParams::MachinePrecision())
  {
    UserError(*this) << "Cannot interpolate to a time (" << t << ") prior to oldest("
                     << first->t << ") in history";
    return;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << " interpolating for t = " << t  << std::endl;
  }

  // If we are within roundoff of the endpoints of the history, just use
  // the endpoints, otherwise interpolate to get it.
  if ( fabs(t-first->t)<Util::MachineDependentParams::MachinePrecision())
  {
    *v1p = first->v1;
    *v2p = first->v2;
  }
  else if ( fabs(t-last->t)<Util::MachineDependentParams::MachinePrecision())
  {
    *v1p = last->v1;
    *v2p = last->v2;
  }
  else
  {
    // If there are no elements of history with time later than t, there
    // is no point using lower_bound to search for one.

    it1=last;
    if (it1->t < t)
    {
      // just use this last point, which will cause us to extrapolate
    }
    else
    {
      LessThan<History,double> lessFunct;
      it1 = lower_bound(history.begin(),history.end(),t,lessFunct);
    }

    // Now it1 points to the first element with time > t (or last element of
    // list if extrapolating)
    t3 = it1->t;
    v13 = it1->v1;
    v23 = it1->v2;
    it1--;
    t2 = it1->t;
    v12 = it1->v1;
    v22 = it1->v2;
    it1--;
    t1 = it1->t;
    v11 = it1->v1;
    v21 = it1->v2;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "Using time t3="<<t3<<" v1(t3)="<<v13<<" v2(t3)="<<v23  << std::endl;
      Xyce::dout() << "Using time t2="<<t2<<" v1(t2)="<<v12<<" v2(t2)="<<v22  << std::endl;
      Xyce::dout() << "Using time t1="<<t1<<" v1(t1)="<<v11<<" v2(t1)="<<v21  << std::endl;
    }

    //  now we have three values of each function to be interpolated, and three
    // times.  t3 is after the desired time, t1 and t2 are before (t2 might be
    // equal to the desired time)
    // Set up the differences for lagrange interpolation:
    dt12 = t1-t2;
    dt13 = t1-t3;
    dt23 = t2-t3;
    dt1 = t-t1;
    dt2 = t-t2;
    dt3 = t-t3;
    // now we set up the lagrange interpolating functions
    // e.g. f1 = (t-t2)*(t-t3)/((t1-t2)*(t1-t3))
    // so that fi is 1 at ti and 0 at the other times.
    f1 = dt2*dt3;
    f2 = dt1*dt3;
    f3 = dt1*dt2;
    if (dt12 != 0)
    {
      f1 /= dt12;
      f2 /= -dt12;
    }
    else
    {
      f1 = f2 = 0.0;
    }
    if (dt13 != 0)
    {
      f1 /= dt13;
      f3 /= -dt13;
    }
    else
    {
      f1 = f2 = 0.0;
    }
    if (dt23 != 0)
    {
      f2 /= dt23;
      f3 /= -dt23;
    }
    else
    {
      f2 = f3 = 0.0;
    }
    // that's it, we have the interpolation functions evaluated at the time t,
    // and the values of v1 and v2 at the points, perform  the interpolation

    double d11=(v13-v12)/(t3-t2);
    double d21=(v12-v11)/(t2-t1);
    double d12=(v23-v22)/(t3-t2);
    double d22=(v22-v21)/(t2-t1);

    // If the derivatives are changing dramatically, don't do quadradic
    // interpolation, just do linear between t2 and t3
    // The conditions here are the same as the conditions that would
    // make us set a breakpoint
    if (fabs(d11-d21) >= .99*std::max(fabs(d11),fabs(d21))+1)
    {
      // linear
      if (fabs(v13-v12)<Util::MachineDependentParams::MachinePrecision())
      {
        // this is a really pathological case where the history
        // after a breakpoint is totally flat.  Either extrapolation or
        // interpolation should just be the average of the two
        *v1p = (v13+v12)/2.0;
      }
      else
      {
        *v1p = v12+d11*(t-t2);
      }
    }
    else
    {
      *v1p = f1*v11+f2*v12+f3*v13;
    }

    if (fabs(d12-d22) >= .99*std::max(fabs(d12),fabs(d22))+1)
    {
      // linear
      if (fabs(v23-v22)<Util::MachineDependentParams::MachinePrecision())
      {
        // this is a really pathological case where the history
        // after a breakpoint is totally flat.  Either extrapolation or
        // interpolation should just be the average of the two
        *v2p = (v23+v22)/2.0;
      }
      else
      {
        *v2p = v22+d12*(t-t2);
      }
    }
    else
    {
      *v2p = f1*v21+f2*v22+f3*v23;
    }
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::getInstanceBreakPoints
// Purpose       : This function adds break points to a vector of breakpoints.
//
//                 It does not bother to check them in any way, or put them
//                 in order.  It only adds them in.
//
// Special Notes : The guts of this has been moved to acceptStep, which
//                 actually computes the breakpoints if needed.  We only add
//                 them to the list here if necessary.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/08/01
//-----------------------------------------------------------------------------
bool Instance::getInstanceBreakPoints ( std::vector<Util::BreakPoint> & breakPointTimes )
{
  bool bsuccess = true;

  double currentTime = getSolverState().currTime_;
  int timeStep = getSolverState().timeStepNumber_;

  //  We're called once prior to any newton iterations, not even the
  // DC Op point.  Never do anything if first_BP_call_done is false.

  if (timeStep != 0 && first_BP_call_done)
  {
    if (newBreakPoint)
    {
      breakPointTimes.push_back(newBreakPointTime);
      newBreakPoint = false;
    }
  }
  else
  {
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << " In Instance::getBreakPoints "<<std::endl;
        Xyce::dout() << " First time step, I don't get to set breakpoints.  Time is ";
        Xyce::dout() << currentTime << std::endl;
      }
  }

  first_BP_call_done=true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance:isConverged ()
// Purpose       : Allow the transmission line to trump the nonlinear
//                 solver's convergence test to reject a step
// Special Notes : The nonlinear solver will not declare a solution converged
//                 unless all devices return true from this function.
//
//                 If the current time step exceeds the time delay *AND*
//                 a discontinuity occured at the last step (which we couldn't
//                 have known until now), then we must reject this time step
//                 and force the time step selection to roll it back.
//
//                 The logic here is similar to acceptStep's, only we don't
//                 need to search history.  The relevant points are always
//                 the last two out of history and the current step's.
//
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 02/24/2021
//-----------------------------------------------------------------------------

inline bool Instance::isConverged()
{
  bool converged = true;

  if ((!getSolverState().dcopFlag) &&
      !(getSolverState().initTranFlag_ &&  getSolverState().newtonIter == 0 ))
  {

    double currentTime = getSolverState().currTime_;

    Linear::Vector *theSolVectorPtr = extData.nextSolVectorPtr;

    if (theSolVectorPtr == 0)
    {
      // we are so early in a restarted run that the solution hasn't even been
      // set up yet.  So let's just say we're NOT converged?
      converged=false;
    }
    else
    {
      std::vector<History>::iterator last = history.end();

      double oVp1,oVp2,oVn1,oVn2,oI1,oI2;

      oVp1 = (*theSolVectorPtr)[li_Pos1];
      oVn1 = (*theSolVectorPtr)[li_Neg1];
      oI1  = (*theSolVectorPtr)[li_Ibr1];
      oVp2 = (*theSolVectorPtr)[li_Pos2];
      oVn2 = (*theSolVectorPtr)[li_Neg2];
      oI2  = (*theSolVectorPtr)[li_Ibr2];

      // we now need to see if the current time and the past two times
      // show a discontinuity.  We do that by computing derivatives and watching
      // for a dramatic change.

      double t3=currentTime;
      double v13=(oVp2-oVn2)+Z0*oI2;
      ;      double v23=(oVp1-oVn1)+Z0*oI1;;

      last--;
      double t2=last->t;
      double v12=last->v1;
      double v22=last->v2;
      last--;
      double t1=last->t;
      double v11=last->v1;
      double v21=last->v2;

      // slope from last time to this time for both ends of line
      double d11=(v13-v12)/(t3-t2);
      double d12=(v23-v22)/(t3-t2);

      // slope from second-to last time to last time for both ends of line
      double d21=(v12-v11)/(t2-t1);
      double d22=(v22-v21)/(t2-t1);

      // If either end of the line has shown a dramatic chaing in slope,
      // we've got a discontinuity.
      if ((fabs(d11-d21) >= .99*std::max(fabs(d11),fabs(d21))+1) ||
          (fabs(d12-d22) >= .99*std::max(fabs(d12),fabs(d22))+1))
      {
        // And if we have a discontinuity *AND* the current time step is
        // larger than the time delay of the line, we're not OK with this step.
        if ( (currentTime - (t2 + td) ) > getSolverState().bpTol_ )
        {
          converged = false;
        }
      }
    }
  }
  return converged;
}

//-----------------------------------------------------------------------------
// Function      : Instance::acceptStep
// Purpose       : This function saves the values of v1 and v2 along with
//                 the current time.  It is to be called ONLY at the point
//                 when the time integrator has determined we've got a
//                 converged, acceptable solution and is accepting it,
//                 but before it's updated its times and rotated vectors.
//
// Special Notes : In SPICE this same stuff was done in the "TRAaccept" function.
//
// Scope         : public
// Creator       : Tom Russo, SNL
// Creation Date : 01/23/07
//-----------------------------------------------------------------------------
void Instance::acceptStep()
{
  if (!getSolverState().dcopFlag)
  {
    double currentTime = getSolverState().currTime_;

    double d11, d21, d12, d22;
    Linear::Vector *theSolVectorPtr = extData.nextSolVectorPtr;// the accepted
    // values from this
    // step

    std::vector<History>::iterator last = history.end();

    last--;  // point to last item, not past last item.

    //  We're called once prior to any newton iterations, not even the
    // DC Op point.  Never do anything if first_BP_call_done is false.
    double oVp1,oVp2,oVn1,oVn2,oI1,oI2;
    double ov1,ov2;
    double tmp_v1,tmp_v2, tmp_t;

    double oVi1,oVi2;
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << " In Instance::acceptStep "<<std::endl;
      Xyce::dout() << "I want breakpoints.  Time is " << currentTime << std::endl;
      Xyce::dout() << "   timeOld is " << timeOld << std::endl;
    }

    // we're the end of a time step, the solution has been accepted.
    // clean up the history by deleting records of times so far
    // back that they'll never be used for interpolation again
    // never try to prune history for anything but times that have
    // been accepted already.
    // TVR: The goal of this was to prune the early history so we don't
    // get unbounded growth of the history vector, with the intent of making
    // the interpolation method faster.  Turns out that deleting these vector
    // elements is very expensive, much more expensive than using "lower_bound"
    // to find a value in the long list.  So I'm commenting this out.
    // double delayedTime;
    //  if (timeOld != -1)
    //  {
    //    delayedTime = timeOld-td;
    //    pruneHistory(delayedTime);
    //  }

    oVp1 = (*theSolVectorPtr)[li_Pos1];
    oVn1 = (*theSolVectorPtr)[li_Neg1];
    oI1  = (*theSolVectorPtr)[li_Ibr1];
    oVp2 = (*theSolVectorPtr)[li_Pos2];
    oVn2 = (*theSolVectorPtr)[li_Neg2];
    oI2  = (*theSolVectorPtr)[li_Ibr2];

    // Having the old values means we can calculate what v1 and v2
    // were for that time.
    ov1=(oVp2-oVn2)+Z0*oI2;
    ov2=(oVp1-oVn1)+Z0*oI1;

    if (DEBUG_DEVICE)
    {
      oVi1 = (*theSolVectorPtr)[li_Int1];
      oVi2 = (*theSolVectorPtr)[li_Int2];

      if (isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << " ----- New time step -----" << std::endl;
        Xyce::dout() << " Last solution : " << std::endl;
        Xyce::dout() << " vpos1 = " << oVp1 << std::endl;
        Xyce::dout() << " vneg1 = " << oVn1 << std::endl;
        Xyce::dout() << " vint1 = " << oVi1 << std::endl;
        Xyce::dout() << " ibr1 = " << oI1 << std::endl;
        Xyce::dout() << " vpos2 = " << oVp2 << std::endl;
        Xyce::dout() << " vneg2 = " << oVn2 << std::endl;
        Xyce::dout() << " vint2 = " << oVi2 << std::endl;
        Xyce::dout() << " ibr2 = " << oI2 << std::endl;
        Xyce::dout() << "in acceptStep, saving for time=" << currentTime << ", V1 = " << ov1 <<  ", V2 = " << ov2  << std::endl;
        Xyce::dout() << " V1V2DBG " << currentTime << " " << ov1 << " " << ov2 << std::endl;
      }
    }

    history.push_back(History(currentTime,ov1,ov2));

    last = history.end();
    last--; // point to last item, not past last item.
    // Now calculate derivatives based on history
    tmp_v1 = last->v1; tmp_v2 = last->v2; tmp_t = last->t;
    last--;
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "tmp_t="  << tmp_t  << " last->t =" << last->t << std::endl;
      Xyce::dout() << "tmp_v1=" << tmp_v1 << " last->v1=" << last->v1 << std::endl;
      Xyce::dout() << "tmp_v2=" << tmp_v2 << " last->v2=" << last->v2 << std::endl;
    }
    d11 = (tmp_v1-last->v1)/(tmp_t-last->t);
    d12 = (tmp_v2-last->v2)/(tmp_t-last->t);
    tmp_v1 = last->v1; tmp_v2 = last->v2; tmp_t = last->t;
    last--;
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "tmp_t="  << tmp_t  << " last->t =" << last->t << std::endl;
      Xyce::dout() << "tmp_v1=" << tmp_v1 << " last->v1=" << last->v1 << std::endl;
      Xyce::dout() << "tmp_v2=" << tmp_v2 << " last->v2=" << last->v2 << std::endl;
    }
    d21 = (tmp_v1-last->v1)/(tmp_t-last->t);
    d22 = (tmp_v2-last->v2)/(tmp_t-last->t);
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "Derivs are " << d11 << " " << d21 << std::endl;
      Xyce::dout() << "   and " << d12 << " " <<d22 << std::endl;
      Xyce::dout() << " fabs(d11-d21) = " << fabs(d11-d21) << std::endl;
      Xyce::dout() << " fabs(d12-d22) = " << fabs(d12-d22) << std::endl;
      Xyce::dout() << "D1D2DBG " << currentTime << " " << d11 << " " << d12 << std::endl;
    }

    if ((fabs(d11-d21) >= .99*std::max(fabs(d11),fabs(d21))+1) ||
        (fabs(d12-d22) >= .99*std::max(fabs(d12),fabs(d22))+1))
    {
      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "Derivative is changing enough, I want to set a break point ";
        Xyce::dout() << td << " ahead of discontinuity, which is ";
        Xyce::dout() << tmp_t+td<<std::endl;
      }
      newBreakPointTime = (tmp_t+td);

      // Only set a breakpoint if it's not effectively the current time
      if ( fabs(currentTime - newBreakPointTime) > getSolverState().bpTol_ )
      {
        newBreakPoint = true;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInternalState
// Purpose       : Generates an DeviceState object and populates
//                 it with the contents of the history vector for use by
//                 restarts
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 09/03/04
//-----------------------------------------------------------------------------

DeviceState * Instance::getInternalState()
{
  int hsize,i,j;
  // allocate object to return
  DeviceState * myState = new DeviceState;


  myState->ID=getName().getEncodedName();
  // We'll pack our history data into the single vector of doubles
  myState->data.resize(history.size()*3);
  hsize=history.size();
  for (i=0;i<hsize;++i)
  {
    j=i*3;
    myState->data[j]=history[i].t;
    myState->data[j+1]=history[i].v1;
    myState->data[j+2]=history[i].v2;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << Xyce::section_divider
         << std::endl;
    Xyce::dout() << " In Instance::getInternalState " << std::endl;
    Xyce::dout() << "   name=" << getName() << std::endl;
    Xyce::dout() << "   history size = " << hsize << std::endl;
    Xyce::dout() << "   history  data: " << std::endl;
    for (i = 0 ; i < hsize ; ++i)
    {
      Xyce::dout() << "   (" << history[i].t << ", " << history[i].v1 << ", "
           << history[i].v2 << ")"<< std::endl;
    }

    Xyce::dout() << "   DeviceState ID = " << myState->ID << std::endl;
    Xyce::dout() << "   DeviceState data size " << myState->data.size() << std::endl;
    Xyce::dout() << "   Device State data: " << std::endl;
    for (i = 0 ; i < myState->data.size() ; ++i)
    {
      Xyce::dout() << "    " << myState->data[i] << std::endl;
    }
    Xyce::dout() << Xyce::section_divider
         << std::endl;
  }

  return myState;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setInternalState
// Purpose       : Reload history data from restart
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 09/03/04
//-----------------------------------------------------------------------------
bool Instance::setInternalState(const DeviceState &state)
{
  int dsize=state.data.size();
  int hsize,i,j;
  if (getName().getEncodedName() != state.ID)
  {
    DevelFatal(*this).in("TRA::Instance::setInternal") << "ID(" << state.ID << ") from restart does not match my name (" << getName() << ")";
    return false;
  }

  if (dsize%3 != 0)
  {
    UserError(*this) << "Data size from restart (" << dsize << ") not a multiple of 3";
    return false;
  }

  hsize=dsize/3;
  history.clear();
  history.resize(hsize);
  for ( i=0; i<hsize; ++i)
  {
    j=i*3;
    history[i].t=state.data[j];
    history[i].v1=state.data[j+1];
    history[i].v2=state.data[j+2];
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << Xyce::section_divider
         << std::endl;
    Xyce::dout() << " In Instance::setInternalState " << std::endl;
    Xyce::dout() << "   name=" << getName() << std::endl;
    Xyce::dout() << "   history size = " << hsize << std::endl;
    Xyce::dout() << "   history  data: " << std::endl;
    for (i = 0 ; i < hsize ; ++i)
    {
      Xyce::dout() << "   (" << history[i].t << ", " << history[i].v1 << ", "
           << history[i].v2 << ")"<< std::endl;
    }

    Xyce::dout() << "   DeviceState ID = " << state.ID << std::endl;
    Xyce::dout() << "   DeviceState data size " << state.data.size() << std::endl;
    Xyce::dout() << "   Device State data: " << std::endl;
    for (i = 0 ; i < state.data.size() ; ++i)
    {
      Xyce::dout() << "    " << state.data[i] << std::endl;
    }
    Xyce::dout() << Xyce::section_divider
         << std::endl;
  }
  return true;
}

// Master class functions

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       : Perform top-level loop over all transmission lines,
//                 calling their updatePrimaryState() functions as needed
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 26 April 2018
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec,
                            int loadType)
{
  bool bsuccess = true;

  // the TRA device is flagged as a nonlinear device because of its
  // use of history --- can't do linear loads here
  // The primary purpose of this conditional is to prevent us
  // from loading if we're doing frequency domain.
  if (loadType == NONLINEAR || loadType == ALL)
  {
    bsuccess = DeviceMaster<Traits>::updateState(solVec,staVec,stoVec);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 26 April 2018
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double * qVec,
                               double * bVec, double * leadF, double * leadQ,
                               double * junctionV, int loadType)
{
  bool bsuccess = true;
  // The primary purpose of this conditional is to prevent us
  // from loading if we're doing frequency domain.
  if (loadType == NONLINEAR || loadType == ALL)
  {
    bsuccess = DeviceMaster<Traits>::loadDAEVectors(solVec,fVec,qVec,
                                                   bVec, leadF, leadQ,
                                                   junctionV);
  }

  return bsuccess;
}
//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       : Perform top-level loop over all transmission lines,
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 26 April 2018
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdX,
                              Linear::Matrix & dQdX,
                              int loadType)
{
  bool bsuccess = true;

  // The primary purpose of this conditional is to prevent us
  // from loading if we're doing frequency domain.
  if (loadType == NONLINEAR || loadType == ALL)
  {
    bsuccess = DeviceMaster<Traits>::loadDAEMatrices(dFdX, dQdX);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadFreqDAEVectors
// Purpose       : Load frequency domain DAE vectors from stored intermediates
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 26 April 2018
//-----------------------------------------------------------------------------
bool Master::loadFreqDAEVectors(double frequency, std::complex<double>* solVec,
                                  std::vector<Util::FreqVecEntry>& fVec,
                                  std::vector<Util::FreqVecEntry>& bVec)
{

  InstanceVector::const_iterator it, end;

  it = getInstanceBegin();
  end = getInstanceEnd();

  fVec.clear();
  bVec.clear();

  for ( ; it != end; ++it )
  {
    Instance & theInstance = *(*it);
    Util::FreqVecEntry tmpEntry;

    std::complex<double> g0 = std::complex<double>(theInstance.G0,0.0);
    std::complex<double> z0 = std::complex<double>(theInstance.Z0,0.0);
    double omega=2*M_PI*frequency;
    std::complex<double> coef(cos(-omega*theInstance.td),sin(-omega*theInstance.td));

    // Pos1
    tmpEntry.val = g0*(solVec[theInstance.li_Pos1] - solVec[theInstance.li_Int1]);
    tmpEntry.lid = theInstance.li_Pos1;
    fVec.push_back(tmpEntry);

    // Pos 2
    tmpEntry.val = g0*(solVec[theInstance.li_Pos2] - solVec[theInstance.li_Int2]);
    tmpEntry.lid = theInstance.li_Pos2;
    fVec.push_back(tmpEntry);

    // Int1
    tmpEntry.val = g0*(solVec[theInstance.li_Int1]-solVec[theInstance.li_Pos1])
      +solVec[theInstance.li_Ibr1];
    tmpEntry.lid = theInstance.li_Int1;
    fVec.push_back(tmpEntry);

    // Int2
    tmpEntry.val = g0*(solVec[theInstance.li_Int2]-solVec[theInstance.li_Pos2])
      +solVec[theInstance.li_Ibr2];
    tmpEntry.lid = theInstance.li_Int2;
    fVec.push_back(tmpEntry);

    // Ibr1
    tmpEntry.val = solVec[theInstance.li_Int1]-solVec[theInstance.li_Neg1]
      + coef*(solVec[theInstance.li_Neg2]-solVec[theInstance.li_Pos2])
      - coef*z0*solVec[theInstance.li_Ibr2];
    tmpEntry.lid = theInstance.li_Ibr1;
    fVec.push_back(tmpEntry);

    // Ibr2
    tmpEntry.val = solVec[theInstance.li_Int2]-solVec[theInstance.li_Neg2]
      + coef*(solVec[theInstance.li_Neg1]-solVec[theInstance.li_Pos1])
      - coef*z0*solVec[theInstance.li_Ibr1];
    tmpEntry.lid = theInstance.li_Ibr2;
    fVec.push_back(tmpEntry);

    // Neg1
    tmpEntry.val = -solVec[theInstance.li_Ibr1];
    tmpEntry.lid = theInstance.li_Neg1;
    fVec.push_back(tmpEntry);

    // Neg2
    tmpEntry.val = -solVec[theInstance.li_Ibr2];
    tmpEntry.lid = theInstance.li_Neg2;
    fVec.push_back(tmpEntry);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadFreqDAEMatrices
// Purpose       : Load frequency domain DAE matrices from stored intermediates
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 26 April 2018
//-----------------------------------------------------------------------------
bool Master::loadFreqDAEMatrices(double frequency, std::complex<double>* solVec,
                                   std::vector<Util::FreqMatEntry>& dFdx)
{
  InstanceVector::const_iterator it, end;

  it = getInstanceBegin();
  end = getInstanceEnd();

  dFdx.clear();

  for ( ; it != end; ++it )
  {
    Instance & theInstance = *(*it);
    Util::FreqMatEntry tmpEntry;

    double g0=theInstance.G0;
    double z0=theInstance.Z0;
    double omega=2*M_PI*frequency;
    std::complex<double> coef(cos(-omega*theInstance.td),sin(-omega*theInstance.td));

    // These are all the same as the time-domain dFdX:
    // First do all the ones that have positive G0
    tmpEntry.val=std::complex<double>(g0,0.0);
    tmpEntry.row_lid = theInstance.li_Pos1;
    tmpEntry.col_lid = theInstance.APos1EquPos1NodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.row_lid = theInstance.li_Int1;
    tmpEntry.col_lid = theInstance.AInt1EquInt1NodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.row_lid = theInstance.li_Pos2;
    tmpEntry.col_lid = theInstance.APos2EquPos2NodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.row_lid = theInstance.li_Int2;
    tmpEntry.col_lid = theInstance.AInt2EquInt2NodeOffset;
    dFdx.push_back(tmpEntry);

    // do all the ones that have negative G0
    tmpEntry.val=std::complex<double>(-g0,0.0);
    tmpEntry.row_lid = theInstance.li_Pos1;
    tmpEntry.col_lid = theInstance.APos1EquInt1NodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.row_lid = theInstance.li_Int1;
    tmpEntry.col_lid = theInstance.AInt1EquPos1NodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.row_lid = theInstance.li_Pos2;
    tmpEntry.col_lid = theInstance.APos2EquInt2NodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.row_lid = theInstance.li_Int2;
    tmpEntry.col_lid = theInstance.AInt2EquPos2NodeOffset;
    dFdx.push_back(tmpEntry);

    // Now all the +1s
    tmpEntry.val=std::complex<double>(1.0,0.0);
    tmpEntry.row_lid = theInstance.li_Int1;
    tmpEntry.col_lid = theInstance.AInt1EquIbr1NodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.row_lid = theInstance.li_Ibr1;
    tmpEntry.col_lid = theInstance.AIbr1EquInt1NodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.row_lid = theInstance.li_Int2;
    tmpEntry.col_lid = theInstance.AInt2EquIbr2NodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.row_lid = theInstance.li_Ibr2;
    tmpEntry.col_lid = theInstance.AIbr2EquInt2NodeOffset;
    dFdx.push_back(tmpEntry);

    // Now all the -1s
    tmpEntry.val=std::complex<double>(-1.0,0.0);
    tmpEntry.row_lid = theInstance.li_Neg1;
    tmpEntry.col_lid = theInstance.ANeg1EquIbr1NodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.row_lid = theInstance.li_Ibr1;
    tmpEntry.col_lid = theInstance.AIbr1EquNeg1NodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.row_lid = theInstance.li_Neg2;
    tmpEntry.col_lid = theInstance.ANeg2EquIbr2NodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.row_lid = theInstance.li_Ibr2;
    tmpEntry.col_lid = theInstance.AIbr2EquNeg2NodeOffset;
    dFdx.push_back(tmpEntry);

    // Now we have to do the ones that are not the same as TD or DC:
    // (Note, these are the same *locations* as the extra DC elements,
    //  but their values are multiplied by coef.
    tmpEntry.val = -coef;
    tmpEntry.row_lid = theInstance.li_Ibr1;
    tmpEntry.col_lid = theInstance.AIbr1EquPos2NodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.row_lid = theInstance.li_Ibr2;
    tmpEntry.col_lid = theInstance.AIbr2EquPos1NodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.val = coef;
    tmpEntry.row_lid = theInstance.li_Ibr1;
    tmpEntry.col_lid = theInstance.AIbr1EquNeg2NodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.row_lid = theInstance.li_Ibr2;
    tmpEntry.col_lid = theInstance.AIbr2EquNeg1NodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.val = -coef*z0;
    tmpEntry.row_lid = theInstance.li_Ibr1;
    tmpEntry.col_lid = theInstance.AIbr1EquIbr2NodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.row_lid = theInstance.li_Ibr2;
    tmpEntry.col_lid = theInstance.AIbr2EquIbr1NodeOffset;
    dFdx.push_back(tmpEntry);

  }
  return true;
}

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 9/25/02
//-----------------------------------------------------------------------------
bool Model::processParams ()
{
  // there are no model parameters to process.
  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirely, PSSI
// Creation Date : 03/23/06
//----------------------------------------------------------------------------
bool Model::processInstanceParams()
{

  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    (*iter)->processParams();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
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
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
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
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
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
    os << "Z0 = " << (*iter)->Z0 << std::endl;
    os << "G0 = " << (*iter)->G0 << std::endl;
    os << "TD = " << (*iter)->td << std::endl;
    os << "FREQ = " << (*iter)->freq << std::endl;
    os << "NL = " << (*iter)->NL << std::endl;

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


//-----------------------------------------------------------------------------
// Function      : Instance::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 8/01/01
//-----------------------------------------------------------------------------
double Instance::getMaxTimeStepSize ()
{
  return td;
}

// Additional Declarations

// History member (trivial) functions

//-----------------------------------------------------------------------------
// Function      : History::History
// Purpose       : default constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 06/14/01
//-----------------------------------------------------------------------------
History::History()
  : t(0),v1(0),v2(0)
{
}

//-----------------------------------------------------------------------------
// Function      : History::History
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 06/14/01
//-----------------------------------------------------------------------------
History::~History()
{
}
//-----------------------------------------------------------------------------
// Function      : History::History
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 06/14/01
//-----------------------------------------------------------------------------
History::History(const History &right)
  : t(right.t),v1(right.v1),v2(right.v2)
{
}

//-----------------------------------------------------------------------------
// Function      : History::History
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 06/14/01
//-----------------------------------------------------------------------------
History::History(double a, double b, double c)
  : t(a),v1(b),v2(c)
{
}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() || (deviceMap.find("T")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("t", 1);
  }
}

} // namespace TRA
} // namespace Device
} // namespace Xyce

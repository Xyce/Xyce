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
// Purpose        : lossy transmission line
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL
//
// Creation Date  : 06/16/10
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceState.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_LTRA.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Message.h>
#include <N_DEV_SolverState.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

#include <N_UTL_BreakPoint.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Functors.h>
#include <N_UTL_Math.h>
#include <N_UTL_MathSpecialFunctions.h>
#include <N_UTL_AssemblyTypes.h>

namespace Xyce {
namespace Device {
namespace LTRA {


void Traits::loadInstanceParameters(ParametricData<LTRA::Instance> &p)
{
  p.addPar("V1", 0.0, &LTRA::Instance::initVolt1)
    .setGivenMember(&LTRA::Instance::initVolt1Given)
    .setUnit(U_VOLT)
    .setDescription("Initial voltage at end 1");

  p.addPar("V2", 0.0, &LTRA::Instance::initVolt2)
    .setGivenMember(&LTRA::Instance::initVolt2Given)
    .setUnit(U_VOLT)
    .setDescription("Initial voltage at end 2");

  p.addPar("I1", 0.0, &LTRA::Instance::initCur1)
    .setGivenMember(&LTRA::Instance::initCur1Given)
    .setUnit(U_AMP)
    .setDescription("Initial current at end 1");

  p.addPar("I2", 0.0, &LTRA::Instance::initCur2)
    .setGivenMember(&LTRA::Instance::initCur2Given)
    .setUnit(U_AMP)
    .setDescription("Initial current at end 2");
}

void Traits::loadModelParameters(ParametricData<LTRA::Model> &p)
{
  p.addPar("R", 0.0, &LTRA::Model::resist)
    .setGivenMember(&LTRA::Model::resistGiven)
    .setUnit(U_OHMMM1)
    .setDescription("Resistance per unit length");

  p.addPar("L", 0.0, &LTRA::Model::induct)
    .setGivenMember(&LTRA::Model::inductGiven)
    .setUnit(U_HMM1)
    .setDescription("Inductance per unit length");

  p.addPar("G", 0.0, &LTRA::Model::conduct)
    .setGivenMember(&LTRA::Model::conductGiven)
    .setUnit(U_OHMM1MM1)
    .setDescription("Conductance per unit length");

  p.addPar("C", 0.0, &LTRA::Model::capac)
    .setGivenMember(&LTRA::Model::capacGiven)
    .setUnit(U_FARADMM1)
    .setDescription("Capacitance per unit length");

  p.addPar("LEN", 0.0, &LTRA::Model::length)
    .setGivenMember(&LTRA::Model::lengthGiven)
    .setUnit(U_METER)
    .setDescription("length of line");

  p.addPar("REL", 1.0, &LTRA::Model::reltol)
    .setGivenMember(&LTRA::Model::reltolGiven)
    .setDescription("Rel. rate of change of deriv. for bkpt");

  p.addPar("ABS", 1.0, &LTRA::Model::abstol)
    .setGivenMember(&LTRA::Model::abstolGiven)
    .setDescription("Abs. rate of change of deriv. for bkpt");

  p.addPar("STEPLIMIT", true, &LTRA::Model::stepLimit)
    .setGivenMember(&LTRA::Model::stepLimitGiven)
    .setUnit(U_LOGIC)
    .setDescription("limit timestep size based on the time constant of the line");

  p.addPar("NOSTEPLIMIT", false, &LTRA::Model::noStepLimit)
    .setGivenMember(&LTRA::Model::noStepLimitGiven)
    .setUnit(U_LOGIC)
    .setDescription("don't limit timestep size based on the time constant of the line");

  p.addPar("COMPLEXSTEPCONTROL", false, &LTRA::Model::lteTimeStepControl)
    .setGivenMember(&LTRA::Model::lteTimeStepControlGiven)
    .setUnit(U_LOGIC)
    .setDescription("do complex time step control using local truncation error estimation");

  p.addPar("LININTERP", false, &LTRA::Model::linInterp)
    .setGivenMember(&LTRA::Model::linInterpGiven)
    .setUnit(U_LOGIC)
    .setDescription("use linear interpolation");

  p.addPar("QUADINTERP", true, &LTRA::Model::quadInterp)
    .setGivenMember(&LTRA::Model::quadInterpGiven)
    .setUnit(U_LOGIC)
    .setDescription("use quadratic interpolation");

  p.addPar("MIXEDINTERP", false, &LTRA::Model::mixedInterp)
    .setGivenMember(&LTRA::Model::mixedInterpGiven)
    .setUnit(U_LOGIC)
    .setDescription("use linear interpolation if quadratic results look unacceptable");

  p.addPar("COMPACTREL", 1.0e-3, &LTRA::Model::stLineReltol)
    .setGivenMember(&LTRA::Model::stLineReltolGiven)
    .setDescription("special reltol for straight line checking");

  p.addPar("COMPACTABS", 1.0e-12, &LTRA::Model::stLineAbstol)
    .setGivenMember(&LTRA::Model::stLineAbstolGiven)
    .setDescription("special abstol for straight line checking");

  p.addPar("TRUNCNR", false, &LTRA::Model::truncNR)
    .setGivenMember(&LTRA::Model::truncNRGiven)
    .setUnit(U_LOGIC)
    .setDescription("use N-R iterations for step calculation in LTRAtrunc");

  p.addPar("TRUNCDONTCUT", false, &LTRA::Model::truncDontCut)
    .setGivenMember(&LTRA::Model::truncDontCutGiven)
    .setUnit(U_LOGIC)
    .setDescription("don't limit timestep to keep impulse response calculation errors low");
}



std::vector< std::vector<int> > Instance::jacStamp;


// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------

Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock&  instance_block,
  Model &               model,
  const FactoryBlock &  factory_block)
  : DeviceInstance(instance_block, configuration.getInstanceParameters(), factory_block),
    model_(model),
    input1(0.0),
    input2(0.0),
    initVolt1(0.0),
    initVolt2(0.0),
    initCur1(0.0),
    initCur2(0.0),
    listSize(0),

    initVolt1Given(false),
    initVolt2Given(false),
    initCur1Given(false),
    initCur2Given(false),

    li_Pos1(-1),
    li_Neg1(-1),

    li_Pos2(-1),
    li_Neg2(-1),

    li_Ibr1(-1),
    li_Ibr2(-1),

    APos1EquPos1NodeOffset(-1),
    APos1EquIbr1NodeOffset(-1),

    ANeg1EquNeg1NodeOffset(-1),
    ANeg1EquIbr1NodeOffset(-1),

    APos2EquPos2NodeOffset(-1),
    APos2EquIbr2NodeOffset(-1),

    ANeg2EquNeg2NodeOffset(-1),
    ANeg2EquIbr2NodeOffset(-1),

    AIbr1EquPos1NodeOffset(-1),
    AIbr1EquNeg1NodeOffset(-1),
    AIbr1EquPos2NodeOffset(-1),
    AIbr1EquNeg2NodeOffset(-1),
    AIbr1EquIbr1NodeOffset(-1),
    AIbr1EquIbr2NodeOffset(-1),

    AIbr2EquPos1NodeOffset(-1),
    AIbr2EquNeg1NodeOffset(-1),
    AIbr2EquPos2NodeOffset(-1),
    AIbr2EquNeg2NodeOffset(-1),
    AIbr2EquIbr1NodeOffset(-1),
    AIbr2EquIbr2NodeOffset(-1),

    pos1Pos1Ptr(0),
    pos1Ibr1Ptr(0),

    neg1Neg1Ptr(0),
    neg1Ibr1Ptr(0),

    pos2Pos2Ptr(0),
    pos2Ibr2Ptr(0),

    neg2Neg2Ptr(0),
    neg2Ibr2Ptr(0),

    ibr1Pos1Ptr(0),
    ibr1Neg1Ptr(0),
    ibr1Pos2Ptr(0),
    ibr1Neg2Ptr(0),
    ibr1Ibr1Ptr(0),
    ibr1Ibr2Ptr(0),

    ibr2Pos1Ptr(0),
    ibr2Neg1Ptr(0),
    ibr2Pos2Ptr(0),
    ibr2Neg2Ptr(0),
    ibr2Ibr1Ptr(0),
    ibr2Ibr2Ptr(0),

    first_BP_call_done(false),
    newBreakPoint(false),
    newBreakPointTime(0.0)
{
  numIntVars   = 2;
  numExtVars   = 4;
  numStateVars = 0;

  devConMap.resize(4);
  devConMap[0] = 1;
  devConMap[1] = 1;
  devConMap[2] = 2;
  devConMap[3] = 2;

  if( jacStamp.empty() )
  {
    jacStamp.resize(6);

    jacStamp[0].resize(2);    // Pos1 equ row
    jacStamp[0][0] = 0;       // pos1Pos1
    jacStamp[0][1] = 4;       // pos1Ibr1

    jacStamp[1].resize(2);    // Neg1 equ row
    jacStamp[1][0] = 1;       // neg1Neg1
    jacStamp[1][1] = 4;       // neg1Ibr1

    jacStamp[2].resize(2);    // Pos2 equ row
    jacStamp[2][0] = 2;       // pos2Pos2
    jacStamp[2][1] = 5;       // pos2Ibr2

    jacStamp[3].resize(2);    // Neg2 equ row
    jacStamp[3][0] = 3;       // neg2Neg2
    jacStamp[3][1] = 5;       // neg2Ibr2

    jacStamp[4].resize(6);    // Ibr1 var row
    jacStamp[4][0] = 0;       // ibr1Pos1
    jacStamp[4][1] = 1;       // ibr1Neg1
    jacStamp[4][2] = 2;       // ibr1Pos2
    jacStamp[4][3] = 3;       // ibr1Neg2
    jacStamp[4][4] = 4;       // ibr1Ibr1
    jacStamp[4][5] = 5;       // ibr1Ibr2

    jacStamp[5].resize(6);    // Ibr2 var row
    jacStamp[5][0] = 0;       // ibr2Pos1
    jacStamp[5][1] = 1;       // ibr2Neg1
    jacStamp[5][2] = 2;       // ibr2Pos2
    jacStamp[5][3] = 3;       // ibr2Neg2
    jacStamp[5][4] = 4;       // ibr2Ibr1
    jacStamp[5][5] = 5;       // ibr2Ibr2
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

}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  bool bsuccess = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/16/10
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
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int>& intLIDVecRef,
                                       const std::vector<int>& extLIDVecRef )
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "[LTRA-DBG-DEV] In Instance::registerLIDs\n\n";
    Xyce::dout() << "name             = " << getName() << std::endl;
    Xyce::dout() << "[LTRA-DBG-DEV] number of internal variables: "
                 << numIntVars << std::endl;
    Xyce::dout() << "[LTRA-DBG-DEV] number of external variables: "
                 << numExtVars << std::endl;
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

  li_Ibr1 = intLIDVec[0];
  li_Ibr2 = intLIDVec[1];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "[LTRA-DBG-DEV]  VARIABLE Indicies " << std::endl;
    Xyce::dout() << "li_Pos1 = " << li_Pos1 << std::endl;
    Xyce::dout() << "li_Neg1 = " << li_Neg1 << std::endl;
    Xyce::dout() << "li_Pos2 = " << li_Pos2 << std::endl;
    Xyce::dout() << "li_Neg2 = " << li_Neg2 << std::endl;
    Xyce::dout() << "li_Ibr1 = " << li_Ibr1 << std::endl;
    Xyce::dout() << "li_Ibr1 = " << li_Ibr1 << std::endl;
    Xyce::dout() << "li_Ibr2 = " << li_Ibr2 << std::endl;
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
  addInternalNode(symbol_table, li_Ibr1, getName(), "branch1");
  addInternalNode(symbol_table, li_Ibr2, getName(), "branch2");
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int>& staLIDVecRef)
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> >& Instance::jacobianStamp() const
{
  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> >& jacLIDVec )
{

  DeviceInstance::registerJacLIDs( jacLIDVec );

  APos1EquPos1NodeOffset = jacLIDVec[0][0];
  APos1EquIbr1NodeOffset = jacLIDVec[0][1];

  ANeg1EquNeg1NodeOffset = jacLIDVec[1][0];
  ANeg1EquIbr1NodeOffset = jacLIDVec[1][1];

  APos2EquPos2NodeOffset = jacLIDVec[2][0];
  APos2EquIbr2NodeOffset = jacLIDVec[2][1];

  ANeg2EquNeg2NodeOffset = jacLIDVec[3][0];
  ANeg2EquIbr2NodeOffset = jacLIDVec[3][1];

  AIbr1EquPos1NodeOffset = jacLIDVec[4][0];
  AIbr1EquNeg1NodeOffset = jacLIDVec[4][1];
  AIbr1EquPos2NodeOffset = jacLIDVec[4][2];
  AIbr1EquNeg2NodeOffset = jacLIDVec[4][3];
  AIbr1EquIbr1NodeOffset = jacLIDVec[4][4];
  AIbr1EquIbr2NodeOffset = jacLIDVec[4][5];

  AIbr2EquPos1NodeOffset = jacLIDVec[5][0];
  AIbr2EquNeg1NodeOffset = jacLIDVec[5][1];
  AIbr2EquPos2NodeOffset = jacLIDVec[5][2];
  AIbr2EquNeg2NodeOffset = jacLIDVec[5][3];
  AIbr2EquIbr1NodeOffset = jacLIDVec[5][4];
  AIbr2EquIbr2NodeOffset = jacLIDVec[5][5];

}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
//
// Purpose       : Sets up pointers!?
//
// Special Notes :
//
// Scope         : public
// Creator       : Gary Hennigan, SNL
// Creation Date : 10/11/2012
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
  Linear::Matrix& dFdx = *(extData.dFdxMatrixPtr);

  pos1Ibr1Ptr = &(dFdx[li_Pos1][APos1EquIbr1NodeOffset]);
  pos1Pos1Ptr = &(dFdx[li_Pos1][APos1EquPos1NodeOffset]);

  neg1Ibr1Ptr = &(dFdx[li_Neg1][ANeg1EquIbr1NodeOffset]);
  neg1Neg1Ptr = &(dFdx[li_Neg1][ANeg1EquNeg1NodeOffset]);

  pos2Ibr2Ptr = &(dFdx[li_Pos2][APos2EquIbr2NodeOffset]);
  pos2Pos2Ptr = &(dFdx[li_Pos2][APos2EquPos2NodeOffset]);

  neg2Ibr2Ptr = &(dFdx[li_Neg2][ANeg2EquIbr2NodeOffset]);
  neg2Neg2Ptr = &(dFdx[li_Neg2][ANeg2EquNeg2NodeOffset]);

  ibr1Ibr1Ptr = &(dFdx[li_Ibr1][AIbr1EquIbr1NodeOffset]);
  ibr1Ibr2Ptr = &(dFdx[li_Ibr1][AIbr1EquIbr2NodeOffset]);
  ibr1Pos1Ptr = &(dFdx[li_Ibr1][AIbr1EquPos1NodeOffset]);
  ibr1Neg1Ptr = &(dFdx[li_Ibr1][AIbr1EquNeg1NodeOffset]);
  ibr1Pos2Ptr = &(dFdx[li_Ibr1][AIbr1EquPos2NodeOffset]);
  ibr1Neg2Ptr = &(dFdx[li_Ibr1][AIbr1EquNeg2NodeOffset]);

  ibr2Ibr1Ptr = &(dFdx[li_Ibr2][AIbr2EquIbr1NodeOffset]);
  ibr2Ibr2Ptr = &(dFdx[li_Ibr2][AIbr2EquIbr2NodeOffset]);
  ibr2Pos1Ptr = &(dFdx[li_Ibr2][AIbr2EquPos1NodeOffset]);
  ibr2Neg1Ptr = &(dFdx[li_Ibr2][AIbr2EquNeg1NodeOffset]);
  ibr2Pos2Ptr = &(dFdx[li_Ibr2][AIbr2EquPos2NodeOffset]);
  ibr2Neg2Ptr = &(dFdx[li_Ibr2][AIbr2EquNeg2NodeOffset]);
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 LTRA instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       : update primary state for one LTRA instance
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------

bool Instance::updatePrimaryState ()
{
  return updateIntermediateVars ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       : update secondary state for one LTRA instance
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------

bool Instance::updateSecondaryState()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one LTRA instance
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  return true;
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
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Instance::getInstanceBreakPoints ( std::vector<Util::BreakPoint>& breakPointTimes )
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
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) &&
        getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "[LTRA-DBG-DEV]  In Instance::getBreakPoints "
                     <<std::endl;
        Xyce::dout() << " First time step, I don't get to set breakpoints.  Time is ";
        Xyce::dout() << currentTime << std::endl;
      }
  }

  first_BP_call_done=true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::acceptStep
// Purpose       :
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
void Instance::acceptStep()
{
  // This stores the voltage and current time history at the ports. Note
  // that both the dc-op and first time step have timeStepNumber 0 so we
  // have to distinguish between them. For the purposes of these
  // histories the dc-op is stored in index 0 and the first time step is
  // at index 1. This is consistent with NG-Spice.
  if (getSolverState().dcopFlag && listSize == 0)
  {
    listSize = 10;

    v1.resize(listSize);
    v2.resize(listSize);
    i1.resize(listSize);
    i2.resize(listSize);
  }
  else if (getSolverState().ltraTimeIndex_ >= listSize)
  {
    listSize += 10;

    v1.resize(listSize);
    v2.resize(listSize);
    i1.resize(listSize);
    i2.resize(listSize);
  }

  // because of DCOP at index 0
  v1[getSolverState().ltraTimeIndex_] = vpos1 - vneg1;
  v2[getSolverState().ltraTimeIndex_] = vpos2 - vneg2;

  i1[getSolverState().ltraTimeIndex_] = currp1;
  i2[getSolverState().ltraTimeIndex_] = currp2;

  // Allocate storage for time history entities
  if (getSolverState().initTranFlag_ && !getSolverState().dcopFlag)
  {
    model_.listSize = 10;

    model_.h1dashCoeffs.resize(model_.listSize);
    model_.h2Coeffs.resize(model_.listSize);
    model_.h3dashCoeffs.resize(model_.listSize);
  }
  else if (!getSolverState().dcopFlag && getSolverState().ltraTimeIndex_ >= model_.listSize)
  {
    model_.listSize += 10;

    model_.h1dashCoeffs.resize(model_.listSize);
    model_.h2Coeffs.resize(model_.listSize);
    model_.h3dashCoeffs.resize(model_.listSize);
  }

  bool compact = false;
  if (getDeviceOptions().tryToCompact && getSolverState().ltraTimeIndex_ >= 2) {

    // Figure out if the last 3 points line on a straight line for all
    // the termainal variables
    double t1 = getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_ - 2];
    double t2 = getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_ - 1];
    double t3 = getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_];

    compact = model_.straightLineCheck_(t1, v1[getSolverState().ltraTimeIndex_-2],
                                        t2, v1[getSolverState().ltraTimeIndex_-1],
                                        t3, v1[getSolverState().ltraTimeIndex_],
                                        model_.stLineReltol,
                                        model_.stLineAbstol);
    if (compact) {
      compact = model_.straightLineCheck_(t1, v2[getSolverState().ltraTimeIndex_-2],
                                          t2, v2[getSolverState().ltraTimeIndex_-1],
                                          t3, v2[getSolverState().ltraTimeIndex_],
                                          model_.stLineReltol,
                                          model_.stLineAbstol);
    }
    if (compact) {
      compact = model_.straightLineCheck_(t1, i1[getSolverState().ltraTimeIndex_-2],
                                          t2, i1[getSolverState().ltraTimeIndex_-1],
                                          t3, i1[getSolverState().ltraTimeIndex_],
                                          model_.stLineReltol,
                                          model_.stLineAbstol);
    }
    if (compact) {
      compact = model_.straightLineCheck_(t1, i2[getSolverState().ltraTimeIndex_-2],
                                          t2, i2[getSolverState().ltraTimeIndex_-1],
                                          t3, i2[getSolverState().ltraTimeIndex_],
                                          model_.stLineReltol,
                                          model_.stLineAbstol);
    }
  }

  if (getSolverState().ltraTimeIndex_ > 0)
  {
    double v1_ = (v1[getSolverState().ltraTimeIndex_] + i1[getSolverState().ltraTimeIndex_] *
                  model_.imped) * model_.attenuation;

    double v2_ = (v1[getSolverState().ltraTimeIndex_ - 1] +
                  i1[getSolverState().ltraTimeIndex_ - 1] *
                  model_.imped) * model_.attenuation;

    double v3_ = getSolverState().ltraTimeIndex_ < 2 ? v2_ : (v1[getSolverState().ltraTimeIndex_-2] +
                                                              i1[getSolverState().ltraTimeIndex_-2] *
                                                              model_.imped) * model_.attenuation;
    double v4_ = (v2[getSolverState().ltraTimeIndex_] +
                  i2[getSolverState().ltraTimeIndex_] *
                  model_.imped) * model_.attenuation;

    double v5_ = (v2[getSolverState().ltraTimeIndex_-1] +
                  i2[getSolverState().ltraTimeIndex_-1] *
                  model_.imped) * model_.attenuation;

    double v6_ = getSolverState().ltraTimeIndex_ < 2 ? v5_ : (v2[getSolverState().ltraTimeIndex_-2] +
                                                              i2[getSolverState().ltraTimeIndex_-2] *
                                                              model_.imped) * model_.attenuation;

    double d1_ = (v1_ - v2_) / (getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_]
                                - getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_-1]);

    double d2_ = getSolverState().ltraTimeIndex_ < 2 ? d1_ :
                 (v2_ - v3_) / (getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_-1]
                                - getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_-2]);

    double d3_ = (v4_ - v5_) / (getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_]
                                - getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_-1]);

    double d4_ = getSolverState().ltraTimeIndex_ < 2 ? d3_ :
                 (v5_ - v6_) / (getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_-1]
                                - getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_-2]);

#define CHECK(a,b,c) (std::max(std::max(a,b),c)-std::min(std::min(a,b),c) >= \
                      fabs(50.0*(getDeviceOptions().reltol/3.0*(a+b+c) + \
                                 getDeviceOptions().abstol)))

    bool tmp_test = (fabs(d1_ - d2_) > model_.reltol * std::max(fabs(d1_), fabs(d2_)) +
                     model_.abstol) && CHECK(v1_,v2_,v3_);

    if (tmp_test || ((fabs(d3_ - d4_)
                      >= model_.reltol * std::max(fabs(d3_), fabs(d4_)) +
                      model_.abstol) && CHECK(v4_,v5_,v6_)))
    {
      // Set breakpoint here
      newBreakPoint = true;
      newBreakPointTime = getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_-1] + model_.td;

      if (DEBUG_DEVICE)
      {
        Xyce::dout() << "[LTRA-DBG-DEV]: At simulation time " << getSolverState().currTime_
                     << " adding a breakpoint at time " << newBreakPointTime << std::endl;
      }
    }
  }

  if (getDeviceOptions().tryToCompact && compact && getSolverState().ltraTimeIndex_ >= 2)
  {
    v1[getSolverState().ltraTimeIndex_-1] = v1[getSolverState().ltraTimeIndex_];
    v2[getSolverState().ltraTimeIndex_-1] = v2[getSolverState().ltraTimeIndex_];
    i1[getSolverState().ltraTimeIndex_-1] = i1[getSolverState().ltraTimeIndex_];
    i2[getSolverState().ltraTimeIndex_-1] = i2[getSolverState().ltraTimeIndex_];

    getSolverState().ltraDoCompact_ = true;
  }

  calculateMaxTimeStep_();

  if (DEBUG_DEVICE)
  {
    Xyce::dout() << "[LTRA-DBG-DEV]: At time: " << getSolverState().currTime_
                 << " max time step set to: " << model_.maxTimeStep << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::calculateMaxTimeStep_
// Purpose       : Calculates a maximum safe time step to avoid excessive
//                 errors
//
// Special Notes :
//
// Scope         : private
// Creator       : Gary Hennigan, SNL
// Creation Date : 12/04/2012
//-----------------------------------------------------------------------------
void Instance::calculateMaxTimeStep_()
{

  Model& model = model_;
  model.maxTimeStep = 1.0e99;

  if (getSolverState().ltraTimeIndex_ < 2)
  {
    model.maxTimeStep = std::min(model.td, model.maxSafeStep);
    return;
  }

  switch (model.specialCase)
  {
    case LTRA_MOD_LC:
    case LTRA_MOD_RLC:

      if (model.stepLimitType == LTRA_MOD_STEPLIMIT)
      {
        model.maxTimeStep = model.td;
      }
      else
      {
        size_t ti = getSolverState().ltraTimeIndex_;

        // Approximate derivative to detect changing slope and adjust
        // time step accordingly
        double i1_ = (v2[ti] * model.admit + i2[ti]) * model.attenuation;
        double i2_ = (v2[ti-1] * model.admit + i2[ti-1]) * model.attenuation;
        double i3_ = (v2[ti-2] * model.admit + i2[ti-2]) * model.attenuation;

        double i4_ = (v1[ti] * model.admit + i1[ti]) * model.attenuation;
        double i5_ = (v1[ti-1] * model.admit + i1[ti-1]) * model.attenuation;
        double i6_ = (v1[ti-2] * model.admit + i1[ti-2]) * model.attenuation;

        double d1_ = (i1_ - i2_) /
          (getSolverState().ltraTimePoints_[ti] - getSolverState().ltraTimePoints_[ti-1]);

        double d2_ = (i2_ - i3_) /
          (getSolverState().ltraTimePoints_[ti-1] - getSolverState().ltraTimePoints_[ti-2]);

        double d3_ = (i4_ - i5_) /
          (getSolverState().ltraTimePoints_[ti] - getSolverState().ltraTimePoints_[ti-1]);

        double d4_ = (i5_ - i6_) /
          (getSolverState().ltraTimePoints_[ti-1] - getSolverState().ltraTimePoints_[ti-2]);

        if ((fabs(d1_-d2_) >= model.reltol * std::max(fabs(d1_), fabs(d2_)) + model.abstol) ||
            (fabs(d3_-d4_) >= model.reltol * std::max(fabs(d3_), fabs(d4_)) + model.abstol))
        {
          model.maxTimeStep = std::min(model.maxTimeStep, model.td);
        }

      }
      break;

    case LTRA_MOD_RC:
    case LTRA_MOD_RG:
      break;

    default:
      DevelFatal(*this).in("Instance::calculateMaxTimeStep_")
        << ": Error. Case not handled in calculateMaxTimeStep_() for LTRA model ";
        return;

  }

  //
  // the above was for the parts of the equations that resemble the
  // lossless equations. Now we need to estimate the local truncation
  // error in each of the three convolution equations, and if possible
  // adjust the timestep so that all of them remain within some bound.
  // Unfortunately, the expression for the LTE in a convolution
  // operation is complicated and costly to evaluate; in addition, no
  // explicit inverse exists.
  //
  // So what we do here (for the moment) is check to see the current
  // error is acceptable. If so, the timestep is not changed. If not,
  // then an estimate is made for the new timestep using a few
  // iterations of the newton-raphson method.
  //
  // modification: we change the timestep to half its previous value
  //
  if ((model.specialCase == LTRA_MOD_RLC) && !(model.truncDontCut))
  {
    model.maxTimeStep = std::min(model.maxTimeStep, model.maxSafeStep);
  }

  // NOTE-GLH: None of the following code has been tested. As far as I
  // can tell there is no user option in Spice3 to turn this bit of code
  // on and as such there is no Xyce option to turn it on. Just to
  // reiterate it has NOT been tested or even run.
  if (model.lteTimeStepControl) {

    double current_lte;
    double tolerance;
    switch (model.specialCase) {

      case LTRA_MOD_RLC:
      case LTRA_MOD_RC:
        tolerance = 7.0 *
                    (getDeviceOptions().reltol * (fabs(input1) + fabs(input2)) + getDeviceOptions().abstol);

        current_lte = model.lteCalculate_(*this, getSolverState().currTime_);

        if (current_lte >= tolerance) {
          if (model.truncNR) {

            double ti = getSolverState().ltraTimeIndex_;
            double x = getSolverState().ltraTimePoints_[ti];
            double y = current_lte;
            for (;;)
            {
              double deriv_delta = 0.01 * (x - getSolverState().ltraTimePoints_[ti-1]);

              if (DEBUG_DEVICE)
              {
                if (deriv_delta <= 0.0)
                Xyce::dout() << "LTRAtrunc: error: timestep is now less than zero" << std::endl;
              }
              double deriv = model.lteCalculate_(*this, x + deriv_delta) - y;

              deriv /= deriv_delta;
              double change = (tolerance - y) / deriv;
              x += change;

              int maxiter=2;
              int iterations=0;
              if (maxiter == 0) {
                if (fabs(change) <= fabs(deriv_delta))
                  break;
              } else {
                iterations++;
                if (iterations >= maxiter)
                  break;
              }
              y = model.lteCalculate_(*this, x);
            }

            double tmp = x - getSolverState().ltraTimePoints_[ti-1];
            model.maxTimeStep = std::min(model.maxTimeStep, tmp);
          }
          else
            model.maxTimeStep *= 0.5;
        }
        break;

      case LTRA_MOD_RG:
      case LTRA_MOD_LC:
        break;

      default:
        DevelFatal(*this).in("Instance::calculateMaxTimeStep_")
          << ": Error. Case not handled in calculateMaxTimeStep_() [2] for LTRA model ";
        return;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInternalState
// Purpose       : Generates an DeviceState object and populates
//                 it with the contents of the history vector for use by
//                 restarts.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
DeviceState * Instance::getInternalState()
{
  int i,j;

  // allocate obiect to return
  DeviceState * myState = new DeviceState;

  myState->ID=getName().getEncodedName();

  // stuff owned by the instance:
  myState->dataInt.resize(2);
  myState->dataInt[0] = listSize;

  int origSize = myState->data.size();
  myState->data.resize(origSize + 4*listSize + 6);

  myState->data[origSize  ]=input1;
  myState->data[origSize+1]=input2;
  myState->data[origSize+2]= initVolt1 ;
  myState->data[origSize+3]= initVolt2 ;
  myState->data[origSize+4]= initCur1 ;
  myState->data[origSize+5]= initCur2 ;

  if (DEBUG_RESTART)
    Xyce::dout() << "LTRA::getInternalState:  input1="<<input1<<" input2="<<input2<<std::endl
                 << "LTRA::getInternalState:  initVolt11="<<initVolt1<<" initVolt2="<<initVolt2<<std::endl
                 << "LTRA::getInternalState:  initCur11="<<initCur1<<" initCur2="<<initCur2<<std::endl;


  for (i=0;i<listSize;++i)
  {
    j=(origSize+6)+i*4;
    myState->data[j  ]=v1[i];
    myState->data[j+1]=v2[i];
    myState->data[j+2]=i1[i];
    myState->data[j+3]=i2[i];
    if (DEBUG_RESTART)
    Xyce::dout() << "LTRA::getInternalState:  v1["<<i<<"]="<<v1[i]
                 <<" v2["<<i<<"]="<<v2[i]
                 <<" i1["<<i<<"]="<<i1[i]
                 <<" i2["<<i<<"]="<<i2[i]
                 <<std::endl;
  }

  // stuff owned by the model:
  //if (!(model_.restartStoredFlag))
  //{
    myState->dataInt[1] = model_.listSize;

    origSize = myState->data.size();
    myState->data.resize(origSize+model_.listSize*3);
    for (i=0;i<model_.listSize;++i)
    {
      j=origSize+i*3;
      myState->data[j]=model_.h1dashCoeffs[i];
      myState->data[j+1]=model_.h2Coeffs[i];
      myState->data[j+2]=model_.h3dashCoeffs[i];
      if (DEBUG_RESTART)
        Xyce::dout() << "LTRA::getInternalState:  h1dashCoeffs["<<i<<"] =" << model_.h1dashCoeffs[i]
                     <<" h2Coeffs["<<i<<"] =" << model_.h2Coeffs[i]
                     <<" h3dashCoeffs["<<i<<"] =" << model_.h3dashCoeffs[i]<<std::endl;
    }

    //model_.restartStoredFlag=true;
  //}

  return myState;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setInternalState
// Purpose       : Reload history data from restart
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Instance::setInternalState(const DeviceState &state)
{
  int i,j;
  if (getName().getEncodedName() != state.ID)
  {
    DevelFatal(*this).in("LTRA::Instance::setInternalState") << "ID(" << state.ID << ") from restart does not match my name (" << getName() <<")";
    return false;
  }

  // stuff owned by the instance:
  listSize=state.dataInt[0];
  v1.clear(); v2.clear(); i1.clear(); i2.clear();
  v1.resize(listSize); v2.resize(listSize); i1.resize(listSize); i2.resize(listSize);

  input1=state.data[0];
  input2=state.data[1];
  initVolt1=state.data[2];
  initVolt2=state.data[3];
  initCur1=state.data[4];
  initCur2=state.data[5];

  if (DEBUG_RESTART)
    Xyce::dout() << "LTRA::setInternalState:  input1="<<input1
                 <<" input2="<<input2<<std::endl
                 << "LTRA::setInternalState:  initVolt11="<<initVolt1
                 <<" initVolt2="<<initVolt2<<std::endl
                 << "LTRA::setInternalState:  initCur11="<<initCur1
                 <<" initCur2="<<initCur2<<std::endl;

  for ( i=0; i<listSize; ++i)
  {
    j=6+i*4;
    v1[i]= state.data[j  ];
    v2[i]= state.data[j+1];
    i1[i]= state.data[j+2];
    i2[i]= state.data[j+3];
    if (DEBUG_RESTART)
      Xyce::dout() << "LTRA::setInternalState:  v1["<<i<<"]="<<v1[i]
                   <<" v2["<<i<<"]="<<v2[i]
                   <<" i1["<<i<<"]="<<i1[i]
                   <<" i2["<<i<<"]="<<i2[i]
                   <<std::endl;
  }

  // stuff owned by the model:
  model_.listSize=state.dataInt[1];

  model_.h1dashCoeffs.clear();
  model_.h2Coeffs.clear();
  model_.h3dashCoeffs.clear();

  model_.h1dashCoeffs.resize(model_.listSize);
  model_.h2Coeffs.resize(model_.listSize);
  model_.h3dashCoeffs.resize(model_.listSize);

  for ( i=0; i<model_.listSize; ++i)
  {
    j=(listSize*4+6)+i*3;
    model_.h1dashCoeffs[i]= state.data[j];
    model_.h2Coeffs[i]= state.data[j+1];
    model_.h3dashCoeffs[i]= state.data[j+2];
    if (DEBUG_RESTART)
      Xyce::dout() << "LTRA::setInternalState:  h1dashCoeffs["<<i<<"] ="
                   << model_.h1dashCoeffs[i]
                   <<" h2Coeffs["<<i<<"] =" << model_.h2Coeffs[i]
                   <<" h3dashCoeffs["<<i<<"] ="
                   << model_.h3dashCoeffs[i]<<std::endl;
  }

  return true;
}

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Model::processParams ()
{

  // If tolerances for the line interpolation aren't given set them to
  // the same as the tolerances for the device options
  if (stLineReltol == 0.0)
    stLineReltol = getDeviceOptions().reltol;
  if (stLineAbstol == 0.0)
    stLineAbstol = getDeviceOptions().abstol;

  // Initialize which case this is based on the nonzero user-specified
  // parameters.
  if ((resist == 0) && (conduct == 0) && (capac != 0) && (induct != 0))
  {
    specialCase = LTRA_MOD_LC;
  }
  else if ((resist != 0) && (conduct == 0) && (capac != 0) && (induct != 0))
  {
    specialCase = LTRA_MOD_RLC;
  }
  else if ((resist != 0) && (conduct == 0) && (capac != 0) && (induct == 0))
  {
    specialCase = LTRA_MOD_RC;
  }
  else if ((resist != 0) && (conduct != 0) && (capac == 0) && (induct == 0))
  {
    specialCase = LTRA_MOD_RG;
  }
  else if ((resist != 0) && (conduct != 0) && (capac != 0) && (induct != 0))
  {
    UserError(*this) << "RLCG line not supported. "
                 << "Modes supported: RC, RG, LC, RLC";
    specialCase = LTRA_MOD_LTRA;
  }
  else if ((conduct != 0) && ((capac != 0) || (induct != 0)))
  {
    UserError(*this) << "Nonzero G (except RG) transmission line not supported. "
                 << "Modes supported: RC, RG, LC, RLC";

    specialCase = LTRA_MOD_LTRA;
  }
  else if ((resist != 0) && (conduct == 0) && (capac == 0) && (induct != 0))
  {
    UserError(*this) << "RL transmission line not supported. "
                 << "Modes supported: RC, RG, LC, RLC";

    specialCase = LTRA_MOD_LTRA;
  }

  if ((resist == 0.0 ? 0 : 1) + (conduct == 0.0 ? 0 : 1) +
      (induct == 0.0 ? 0 : 1) + (capac == 0.0 ? 0 : 1) <= 1)
  {
    UserError(*this) << "Invalid specification. Specify at least "
                 << "two of R, L, G, or C with nonzero values. "
                 << "Modes supported: RC, RG, LC, RLC";

    specialCase = LTRA_MOD_LTRA;
  }

  // Override the interpolation, either default or user specified, if
  // the TRYTOCOMPACT option is specified.
  if (getDeviceOptions().tryToCompact) {
    howToInterp = LTRA_MOD_LININTERP;
  }

  if (stepLimit && noStepLimit)
  {
    UserWarning(*this) << "Conflicting options STEPLIMIT and NOSTEPLIMIT given. Using STEPLIMIT";
    stepLimitType = LTRA_MOD_STEPLIMIT;
  }
  else if (stepLimit || !noStepLimit)
  {
    stepLimitType = LTRA_MOD_STEPLIMIT;
  }
  else if (noStepLimit || !stepLimit)
  {
    stepLimitType = LTRA_MOD_NOSTEPLIMIT;
  }
  else
  { // default
    stepLimitType = LTRA_MOD_STEPLIMIT;
  }

  // Calculate some derived parameters
  switch (specialCase)
  {

  case LTRA_MOD_LC:
    imped = sqrt(induct / capac);
    admit = 1.0 / imped;
    td = sqrt(induct*capac) * length;
    attenuation = 1.0;
    break;

  case LTRA_MOD_RLC:
    imped = sqrt(induct / capac);
    admit = 1.0 / imped;
    td = sqrt(induct * capac) * length;
    alpha = 0.5 * (resist / induct);
    beta = alpha;
    attenuation = exp(-beta * td);

    if (alpha > 0.0)
    {
      intH1dash = -1.0;
      intH2 = 1.0 - attenuation;
      intH3dash = -attenuation;
    }
    else
    {
      intH1dash = intH2 = intH3dash = 0.0;
    }

    // Sanity check
    if (alpha < 0.0) {
      UserError(*this) << "Resistance and inductance must be greater than zero";
      return false;
    }

    // Calculate the time step size limit in order to keep
    // impulse-response errors low
    if (!truncDontCut) {
      double xbig, xsmall, xmid, y1big, y1small, y1mid;
      int done = 0, maxiter = 50, iters = 0;

      xbig = 10.0 * td;
      xsmall = td;
      xmid = 0.5 * (xbig + xsmall);
      y1small = rlcH2Func_(xsmall, td, alpha, beta);
      iters = 0;
      for (;;) {

        iters++;
        y1big = rlcH2Func_(xbig, td, alpha, beta);
        y1mid = rlcH2Func_(xmid, td, alpha, beta);
        done = straightLineCheck_(xbig, y1big, xmid, y1mid, xsmall,
                                    y1small, stLineReltol,
                                    stLineAbstol) +
          straightLineCheck_(xbig, y1big, xmid, y1mid, xsmall,
                               y1small, stLineReltol,
                               stLineAbstol);
        if ((done == 2) || (iters > maxiter))
          break; // out of "for (;;)"
        xbig = xmid;
        xmid = 0.5 * (xbig + xsmall);
      }
      maxSafeStep = xbig - td;
    }
    break;

  case LTRA_MOD_RC:
    cByR = capac / resist;
    rclsqr = resist * capac * length * length;
    intH1dash = 0.0;
    intH2 = 1.0;
    intH3dash = 0.0;
    break;

  case LTRA_MOD_RG:
    break;

  default:
    {
      UserError(*this) << "Unhandled LTRA special case encountered.";
    }
  }
  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
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
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock&     MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),

    h1dashFirstVal(0.0),
    h2FirstVal(0.0),
    h3dashFirstVal(0.0),
    h1dashFirstCoeff(0.0),
    h2FirstCoeff(0.0),
    h3dashFirstCoeff(0.0),
    listSize(0),
    resist(0.0),
    induct(0.0),
    conduct(0.0),
    capac(0.0),
    length(0.0),
    reltol(0.0),
    abstol(0.0),

    noStepLimit(false),
    stepLimit(true),
    stepLimitType(LTRA_MOD_STEPLIMIT),

    linInterp(false),
    quadInterp(true),
    mixedInterp(false),

    stLineReltol(0.0),
    stLineAbstol(0.0),

    lteTimeStepControl(false),

    truncNR(false),
    truncDontCut(false),

    resistGiven(false),
    inductGiven(false),
    conductGiven(false),
    capacGiven(false),
    lengthGiven(false),
    reltolGiven(false),
    abstolGiven(false),
    noStepLimitGiven(false),
    stepLimitGiven(false),
    linInterpGiven(false),
    quadInterpGiven(false),
    mixedInterpGiven(false),
    stLineReltolGiven(false),
    stLineAbstolGiven(false),
    truncNRGiven(false),
    truncDontCutGiven(false),

    td(0.0),
    imped(0.0),
    admit(0.0),
    alpha(0.0),
    beta(0.0),
    attenuation(0.0),
    cByR(0.0),
    rclsqr(0.0),
    intH1dash(0.0),
    intH2(0.0),
    intH3dash(0.0),

    coshlrootGR(0.0),
    rRsLrGRorG(0.0),
    rGsLrGRorR(0.0),
    auxIndex(0),
    chopReltol(0.0),
    chopAbstol(0.0),

    maxSafeStep(1.0e99),
    maxTimeStep(1.0e99),
    howToInterp(LTRA_MOD_QUADINTERP),
    printFlag(false),
    specialCase(0),
    tdover(false),
    restartStoredFlag(false)

{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // initial sanity check on parameter
  if ( (resist < 0) || (conduct < 0) || (capac < 0) ||
       (induct < 0) || (length < 0) )
  {
    UserError(*this) << "R, L, C, G or length (LEN) must not be negative";
  }

  processParams ();
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
bool Instance::setIC()
{

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
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
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
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
// Function      : Model::modelCalculations_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Model::modelCalculations_(int& isaved,
                                         double& qf1, double& qf2, double& qf3,
                                         double& lf2, double& lf3)
{
  double t1(0.0),t2(0.0),t3(0.0);
  double dummy1(0.0), dummy2(0.0);

  // Initialize index
  isaved = 0;

  if( getSolverState().dcopFlag)
  {
    switch (specialCase)
    {
      case LTRA_MOD_RG:
        dummy1 = length*sqrt(resist*conduct);
        dummy2 = exp(-dummy1);
        dummy1 = exp(dummy1); // warning: may overflow!
        coshlrootGR = 0.5*(dummy1 + dummy2);

        if (conduct <= 1.0e-10)
        { // hack!
          rRsLrGRorG = length*resist;
        }
        else
        {
          rRsLrGRorG = 0.5*(dummy1 - dummy2)*sqrt(resist/conduct);
        }

        if (resist <= 1.0e-10)
        { // hack!
          rGsLrGRorR = length*conduct;
        }
        else
        {
          rGsLrGRorR = 0.5*(dummy1 - dummy2)*sqrt(conduct/resist);
        }
        break;

      case LTRA_MOD_RC:
      case LTRA_MOD_LC:
      case LTRA_MOD_RLC:
        // simple resistor-like behaviour nothing to set up
        break;

      default:
        DevelFatal(*this).in("Instance::modelCalculations_")
          << ": Error. Case not handled in modelCalculations_() for LTRA model ";
        return false;
    }
  }
  else
  {
    switch (specialCase)
    {
      case LTRA_MOD_RLC:
      case LTRA_MOD_LC:

        if (getSolverState().currTime_ > td)
          tdover = true;
        else
          tdover = false;

      default:
        break;
    }

    switch (specialCase)
    {
      case LTRA_MOD_RLC:

        // set up lists of values of the functions at the
        // necessary timepoints.

        // set up coefficient lists LTRAh1dashCoeffs,
        // LTRAh2Coeffs, LTRAh3dashCoeffs for current
        // timepoint

        // NOTE: h1, h2 and h3 here actually refer to h1tilde,
        // h2tilde, h3tilde in the paper

        // Note: many function evaluations are saved by doing
        // the following all together in one procedure

        (void) rlcCoeffsSetup_( h1dashFirstCoeff, h2FirstCoeff, h3dashFirstCoeff,
                                h1dashCoeffs, h2Coeffs, h3dashCoeffs,
                                listSize,
                                td, alpha, beta,
                                getSolverState().currTime_,
                                getSolverState().ltraTimePoints_,
                                getSolverState().ltraTimeIndex_,
                                chopReltol,
                                &(auxIndex));

      case LTRA_MOD_LC:
        // setting up the coefficients for interpolation
        if (tdover)
        { // serious hack -fix!
          int i = 0;
          for (i = getSolverState().ltraTimeIndex_; i >= 0; i--)
          {
            if (getSolverState().ltraTimePoints_[i] < getSolverState().currTime_ - td)
              break;

          }
          if (DEBUG_DEVICE)
          {
            if (i == getSolverState().ltraTimeIndex_)
          {
            Xyce::dout() << "[LTRA-DBG-DEV] LTRAload: Warning: timestep larger than delay of line" << std::endl;
            Xyce::dout() << "\tTime now: " << getSolverState().currTime_ << std::endl << std::endl;
          }
          }

          if (i == getSolverState().ltraTimeIndex_)
            i--;

          if (i == -1)
          {
            if (DEBUG_DEVICE)
            {
              Xyce::dout() << "[LTRA-DBG-DEV] LTRAload: mistake: cannot find delayed timepoint" << std::endl;
            }
            DevelFatal(*this).in("Instance::modelCalculations_")
              << "Error. Delayed time point not found for LTRA model. "
              << "Zero length line is one possible cause.";
            return false;
          }

          isaved = i;

          t2 = getSolverState().ltraTimePoints_[i];
          t3 = getSolverState().ltraTimePoints_[i+1];

          // quadratic interpolation
          if ((i != 0) && ((howToInterp == LTRA_MOD_QUADINTERP)|| (howToInterp == LTRA_MOD_MIXEDINTERP)))
          {
            t1 = getSolverState().ltraTimePoints_[i-1];
            quadInterp_(getSolverState().currTime_-td, t1, t2, t3, qf1, qf2, qf3);
          }

          // linear interpolation
          if ( (i == 0) || (howToInterp == LTRA_MOD_MIXEDINTERP) || (howToInterp ==  LTRA_MOD_LININTERP))
          {
            linInterp_(getSolverState().currTime_-td, t2, t3, lf2, lf3);
          }
        }

        // interpolation coefficients set-up
        break;

      case LTRA_MOD_RC:

        //
        // set up lists of values of the coefficients at the
        // necessary timepoints.
        //

        //  set up coefficient lists LTRAh1dashCoeffs, LTRAh2Coeffs,
        //    LTRAh3dashCoeffs for current timepoint

        // Note: many function evaluations are saved by doing the
        // following all together in one procedure
        //

        (void)
          rcCoeffsSetup_(h1dashFirstCoeff, h2FirstCoeff, h3dashFirstCoeff,
                         h1dashCoeffs, h2Coeffs, h3dashCoeffs,
                         listSize,
                         cByR, rclsqr,
                         getSolverState().currTime_,
                         getSolverState().ltraTimePoints_,
                         getSolverState().ltraTimeIndex_,
                         chopReltol);

        break;

      case LTRA_MOD_RG:
        break;

      default:
        return false;
        //return(E_BADPARM);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::quadInterp_
// Purpose       :
//
// quadratic interpolation function
// t = timepoint where value wanted
// t1, t2, t3 are three timepoints where the value is known
// c1, c2, c3 are set to the proper coefficients by the function
// the interpolated value is c1*v1 + c2*v2 + c3*v3; this should be
// done in the calling program; (v1,v2,v3 are the known values at
// t1,t2,t3)
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
int Model::quadInterp_ (double t, double t1, double t2, double t3, double& c1, double& c2, double& c3)
{
  double f1, f2, f3;

  if (t == t1)
  {
    c1 = 1.0;
    c2 = 0.0;
    c3 = 0.0;
    return(0);
  }
  if (t == t2)
  {
    c1 = 0.0;
    c2 = 1.0;
    c3 = 0.0;
    return(0);
  }
  if (t == t3)
  {
    c1 = 0.0;
    c2 = 0.0;
    c3 = 1.0;
    return(0);
  }
  if( (t2-t1)==0  || (t3-t2) == 0 || (t1 - t3) ==0) return(1);

  f1 = (t - t2) * (t - t3) ;
  f2 = (t - t1) * (t - t3) ;
  f3 = (t - t1) * (t - t2) ;
  if((t2-t1)==0)
  { //  should never happen, but don't want
    //  to divide by zero, EVER...
    f1=0;
    f2=0;
  }
  else
  {
    f1 /= (t1-t2);
    f2 /= (t2-t1);
  }
  if((t3-t2)==0)
  { //  should never happen, but don't want
    //  to divide by zero, EVER...
    f2=0;
    f3=0;
  }
  else
  {
    f2 /= (t2-t3);
    f3 /= (t2-t3);
  }
  if((t3-t1)==0)
  {  // should never happen, but don't want
     // to divide by zero, EVER...
    f1=0;
    f2=0;
  }
  else
  {
    f1 /= (t1-t3);
    f3 /= (t1-t3);
  }
  c1 = f1;
  c2 = f2;
  c3 = f3;
  return(0);
}

//-----------------------------------------------------------------------------
// Function      : Model::linInterp_
// Purpose       : linear interpolation
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
int Model::linInterp_  (double t, double t1, double t2, double& c1, double& c2)
{
  double temp;

  if (t1 == t2) return(1);

  if (t==t1)
  {
    c1 = 1.0;
    c2 = 0.0;
    return(0);
  }

  if (t==t2)
  {
    c1 = 0.0;
    c2 = 1.0;
    return(0);
  }

  temp = (t-t1)/(t2-t1);
  c2 = temp;
  c1 = 1-temp;

  return(0);
}

//-----------------------------------------------------------------------------
// Function      : Model::intlinfunc_
// Purpose       :
//
// intlinfunc returns \int_lolimit^hilimit h(\tau) d \tau, where
// h(\tau) is assumed to be linear, with values lovalue and hivalue
// \tau = t1 and t2 respectively
// this is used only locally
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::intlinfunc_ (double lolimit, double hilimit,
                                     double lovalue, double hivalue,
                                     double t1, double t2)
{
  double width, m;

  width = t2 - t1;
  if (width == 0.0) return(0.0);
  m = (hivalue - lovalue)/width;

  return( (hilimit-lolimit)*lovalue + 0.5*m*((hilimit-t1)*(hilimit-t1)
    - (lolimit - t1)*(lolimit - t1)));
}

//-----------------------------------------------------------------------------
// Function      : Model::twiceintlinfunc_
// Purpose       :
//
// twiceintlinfunc returns \int_lolimit^hilimit \int_otherlolimit^\tau
// h(\tau') d \tau' d \tau , where
// h(\tau') is assumed to be linear, with values lovalue and hivalue
// \tau = t1 and t2 respectively
// this is used only locally
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::twiceintlinfunc_(double lolimit, double hilimit,
                                         double otherlolimit, double lovalue,
                                         double hivalue, double t1, double t2)
{
  double width, m, dummy;
  double temp1, temp2, temp3;

  width = t2 - t1;
  if (width == 0.0) return(0.0);
  m = (hivalue - lovalue)/width;

  temp1 = hilimit - t1;
  temp2 = lolimit - t1;
  temp3 = otherlolimit - t1;
  dummy = lovalue*((hilimit - otherlolimit)*(hilimit - otherlolimit) -
    (lolimit - otherlolimit)*(lolimit - otherlolimit));
  dummy += m*((temp1*temp1*temp1 - temp2*temp2*temp2)/3.0 -
    temp3*temp3*(hilimit - lolimit));
  return(dummy*0.5);
}


//-----------------------------------------------------------------------------
// Function      : Model::thriceintlinfunc_
// Purpose       :
//
// thriceintlinfunc returns \int_lolimit^hilimit \int_secondlolimit^\tau
// \int_thirdlolimit^\tau' h(\tau'') d \tau'' d \tau' d \tau , where
// h(\tau'') is assumed to be linear, with values lovalue and hivalue
// \tau = t1 and t2 respectively
// this is used only locally
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::thriceintlinfunc_(double lolimit, double hilimit,
                                          double secondlolimit, double thirdlolimit,
                                          double lovalue, double  hivalue, double t1, double t2)
{
  double width, m, dummy;
  double temp1, temp2, temp3, temp4;
  double temp5, temp6, temp7, temp8, temp9, temp10;


  width = t2 - t1;
  if (width == 0.0) return(0.0);
  m = (hivalue - lovalue)/width;

  temp1 = hilimit - t1;
  temp2 = lolimit - t1;
  temp3 = secondlolimit - t1;
  temp4 = thirdlolimit - t1;
  temp5 = hilimit - thirdlolimit;
  temp6 = lolimit - thirdlolimit;
  temp7 = secondlolimit - thirdlolimit;
  temp8 = hilimit - lolimit;
  temp9 = hilimit - secondlolimit;
  temp10 = lolimit - secondlolimit;
  dummy = lovalue*((temp5*temp5*temp5 - temp6*temp6*temp6)/3 -
    temp7*temp5*temp8);
  dummy += m*(((temp1*temp1*temp1*temp1 - temp2*temp2*temp2*temp2)*0.25 -
    temp3*temp3*temp3*temp8)/3 - temp4*temp4*0.5*(temp9*temp9 -
    temp10*temp10));
  return(dummy*0.5);
}

//-----------------------------------------------------------------------------
// Function      : Model::rlcH1dashFunc_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::rlcH1dashFunc_(double time, double T, double alpha, double beta)
{
  double besselarg, exparg, returnval;
  // T is not used in this function

  // result = alpha * e^{- beta*time} * {I_1(alpha*time) -
  // I_0(alpha*time)}
  //

  if (alpha == 0.0) return(0.0);

  exparg = - beta * time;
  besselarg = alpha*time;

  returnval = (Xyce::Util::besselI1(besselarg)-Xyce::Util::besselI0(besselarg))* alpha * exp(exparg);
  return(returnval);
}

//-----------------------------------------------------------------------------
// Function      : Model::rlcH2Func_
// Purpose       : first impulse response function
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::rlcH2Func_(double time, double T, double alpha, double beta)
{
  double besselarg, exparg, returnval;

  //
  // result = 0, time < T
  //      = (alpha*T*e^{-beta*time})/sqrt(t^2 - T^2) *
  //        I_1(alpha*sqrt(t^2 - T^2)), time >= T
  //

  if (alpha == 0.0) return(0.0);
  if (time < T) return(0.0);

  if (time != T) {
    besselarg = alpha*sqrt(time*time - T*T);
  } else {
    besselarg = 0.0;
  }
  exparg = -beta*time;

  returnval = alpha*alpha*T*exp(exparg)*Xyce::Util::besselI1xOverX(besselarg);
  return(returnval);
}

//-----------------------------------------------------------------------------
// Function      : Model::rlcH3dashFunc_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::rlcH3dashFunc_(double time, double T, double alpha, double beta)
{
  double exparg,besselarg,returnval;

  //
  // result = 0, time < T
  //      = alpha*e^{-beta*time}*(t/sqrt(t^2-T^2)*
  //      I_1(alpha*sqrt(t^2-T^2)) - I_0(alpha*sqrt(t^2-T^2)))
  //

  if (alpha == 0.0) return(0.0);
  if (time < T) return(0.0);

  exparg = - beta*time;
  if (time != T) {
    besselarg = alpha*sqrt(time*time - T*T);
  } else {
    besselarg = 0.0;
  }

  returnval = alpha*time*Xyce::Util::besselI1xOverX(besselarg) - Xyce::Util::besselI0(besselarg);
  returnval *= alpha*exp(exparg);
  return(returnval);
}

//-----------------------------------------------------------------------------
// Function      : Model::rlcH1dashTwiceIntFunc_
//
// Purpose       : Twice repeated integral of h1dash for the
//                 special case of G = 0
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::rlcH1dashTwiceIntFunc_(double time, double beta)
{
  double arg, returnval;

  //  result = time * e^{- beta*time} * {I_0(beta*time) +
  // I_1(beta*time)} - time
  //

  if (beta == 0.0) return(time);
  arg = beta*time;
  if (arg == 0.0) return(0.0);

  returnval = (Xyce::Util::besselI1(arg)+Xyce::Util::besselI0(arg))* time * exp(-arg) - time;
  return(returnval);
}

//-----------------------------------------------------------------------------
// Function      : Model::rlcH3dashIntFunc_
//
// Purpose       : twice repeated integral of h1dash for the
//                 special case of G = 0
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::rlcH3dashIntFunc_(double time, double T, double beta)
{
  double exparg, besselarg;
  double returnval;

  if (time <= T) return(0.0);
  if (beta == 0.0) return(0.0);
  exparg = -beta*time;
  besselarg = beta*sqrt(time*time-T*T);
  returnval = exp(exparg)* Xyce::Util::besselI0(besselarg) - exp(-beta*T);
  return(returnval);
}

//-----------------------------------------------------------------------------
// Function      : Model::rcH1dashTwiceIntFunc_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::rcH1dashTwiceIntFunc_(double time, double cbyr)
{
  return(sqrt(4*cbyr*time/M_PI));
}

//-----------------------------------------------------------------------------
// Function      : Model::rcH2TwiceIntFunc_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::rcH2TwiceIntFunc_(double time, double rclsqr)
{
  double temp(0.0);
  if (time != 0.0)
  {
    temp = rclsqr/(4*time);

    double erfc_res = Xyce::Util::erfc(sqrt(temp));

    return((time + rclsqr*0.5)*erfc_res - sqrt(time*rclsqr/M_PI)*exp(- temp));
  }
  else
  {
    return(0.0);
  }
}

//-----------------------------------------------------------------------------
// Function      : Model::rcH3dashTwiceIntFunc_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::rcH3dashTwiceIntFunc_(double time, double cbyr, double rclsqr)
{
  double temp;
  if (time != 0.0)
  {
    temp =  rclsqr/(4*time);

    // see note in rcH2TwiceIntFunc_ about intel compilers
    double erfc_res = Xyce::Util::erfc(sqrt(temp));

    temp = 2*sqrt(time/M_PI)*exp(-temp) - sqrt(rclsqr)*erfc_res;
    return(sqrt(cbyr)*temp);
  }
  else
  {
    return(0.0);
  }
}

// coefficient setups:
//-----------------------------------------------------------------------------
// Function      : Model::rcCoeffsSetup_
// Purpose       :
//
// Sets up the all coefficient lists for the special case where L=G=0
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
void Model::rcCoeffsSetup_(
    double& h1dashfirstcoeff,
    double& h2firstcoeff,
    double& h3dashfirstcoeff,
    std::vector<double>& h1dashcoeffs,
    std::vector<double>& h2coeffs,
    std::vector<double>& h3dashcoeffs,
    int listsize, double cbyr, double rclsqr, double curtime,
    const std::vector<double>& timelist, int timeindex, double reltol)
{
  double delta1;
  double h1dummy1, h1dummy2;
  double h2dummy1, h2dummy2;
  double h3dummy1=0.0, h3dummy2;
  double hilimit1;
  double h1lovalue1,h1hivalue1,h1hivalue2;
  double h2lovalue1,h2hivalue1,h2hivalue2;
  double h3lovalue1,h3hivalue1,h3hivalue2;
  double temp, temp2, temp3, temp4, temp5;
  double h1relval, h2relval, h3relval;
  int doh1=1, doh2=1, doh3=1;
  int i,auxindex;

  // coefflists should already have been allocated to the necessary size

  if (DEBUG_DEVICE)
  {
    if (listsize < timeindex) {
    Xyce::dout() << "[LTRA-DBG-DEV]: LTRAcoeffSetup: not enough space in coefflist" << std::endl;
  }
  }

  auxindex = timeindex;

  // the first coefficients

  delta1 = curtime - timelist[auxindex];
  hilimit1 = delta1;

  h1lovalue1 = 0.0;
  h1hivalue1 = // LTRArcH1dashTwiceIntFunc(hilimit1,cbyr);
        sqrt(4*cbyr*hilimit1/M_PI);
  h1dummy1 = h1hivalue1/delta1;
  h1dashfirstcoeff = h1dummy1;
  h1relval = fabs(h1dummy1*reltol);

  temp = rclsqr/(4*hilimit1);

  // see note in :rcH2TwiceIntFunc_ re intel compiler
  temp2 = (temp >= 100.0 ? 0.0 : Xyce::Util::erfc(sqrt(temp)));
  temp3 = exp(-temp);
  temp4 = sqrt(rclsqr);
  temp5 = sqrt(cbyr);

  h2lovalue1 = 0.0;
  h2hivalue1 = // LTRArcH2TwiceIntFunc(hilimit1,rclsqr);
  (hilimit1 != 0.0?  (hilimit1 + rclsqr*0.5)*temp2 - sqrt(hilimit1*rclsqr/M_PI)*temp3 : 0.0);


  h2dummy1 = h2hivalue1/delta1;
  h2firstcoeff = h2dummy1;
  h2relval = fabs(h2dummy1*reltol);

  h3lovalue1 = 0.0;
  h3hivalue1 = // LTRArcH3dashTwiceIntFunc(hilimit1,cbyr,rclsqr);
    (hilimit1 != 0.0? temp = 2*sqrt(hilimit1/M_PI)*temp3 - temp4*temp2, (temp5*temp): 0.0);

  h3dummy1 = h3hivalue1/delta1;
  h3dashfirstcoeff = h3dummy1;
  h3relval = fabs(h3dummy1*reltol);

  // the coefficients for the rest of the timepoints

  for (i=auxindex; i>0; i--)
  {
    delta1 = timelist[i] - timelist[i - 1];
    hilimit1 = curtime - timelist[i - 1];

    if (doh1)
    {
      h1hivalue2 = h1hivalue1; //previous hivalue1
      h1dummy2 = h1dummy1; // previous dummy1

      h1lovalue1 = h1hivalue2;
      h1hivalue1 = // LTRArcH1dashTwiceIntFunc(hilimit1,cbyr);
            sqrt(4*cbyr*hilimit1/M_PI);
      h1dummy1 = (h1hivalue1 - h1lovalue1)/delta1;
      h1dashcoeffs[i] = h1dummy1 - h1dummy2;
      if (fabs(h1dashcoeffs[i]) < h1relval) doh1=0;
    }
    else
    {
      h1dashcoeffs[i] = 0.0;
    }

    if (doh2 || doh3) {
    temp = rclsqr/(4*hilimit1);
    // see note in :rcH2TwiceIntFunc_ re intel compiler
    temp2 = (temp >= 100.0 ? 0.0 : Xyce::Util::erfc(sqrt(temp)));
    temp3 = exp(-temp);
    }

    if (doh2)
    {
      h2hivalue2 = h2hivalue1; // previous hivalue1
      h2dummy2 = h2dummy1; // previous dummy1

      h2lovalue1 = h2hivalue2;
      h2hivalue1 = // LTRArcH2TwiceIntFunc(hilimit1,rclsqr);
          (hilimit1 != 0.0?  (hilimit1 + rclsqr*0.5)*temp2 - sqrt(hilimit1*rclsqr/M_PI)*temp3 : 0.0);
      h2dummy1 = (h2hivalue1 - h2lovalue1)/delta1;
      h2coeffs[i] = h2dummy1 - h2dummy2;
      if (fabs(h2coeffs[i]) < h2relval) doh2=0;
    }
    else
    {
      h2coeffs[i] = 0.0;
    }

    if (doh3)
    {
      h3hivalue2 = h3hivalue1; // previous hivalue1
      h3dummy2 = h3dummy1; // previous dummy1

      h3lovalue1 = h3hivalue2;
      h3hivalue1 = // LTRArcH3dashTwiceIntFunc(hilimit1,cbyr,rclsqr);
          (hilimit1 != 0.0? temp = 2*sqrt(hilimit1/M_PI)*temp3 - temp4*temp2, (temp5*temp): 0.0);
      h3dummy1 = (h3hivalue1 - h3lovalue1)/delta1;
      h3dashcoeffs[i] = h3dummy1 - h3dummy2;
      if (fabs(h3dashcoeffs[i]) < h3relval) doh3=0;
    }
    else
    {
      h3dashcoeffs[i] = 0.0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Model::rlcCoeffsSetup_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
void Model::rlcCoeffsSetup_(
    double& h1dashfirstcoeff,
    double& h2firstcoeff,
    double& h3dashfirstcoeff,
    std::vector<double>& h1dashcoeffs,
    std::vector<double>& h2coeffs,
    std::vector<double>& h3dashcoeffs,
    int listsize,
    double T, double alpha, double beta, double curtime,
    const std::vector<double>& timelist, int timeindex, double reltol, int* auxindexptr)
{
  unsigned exact;
  double lolimit1,lolimit2,hilimit1,hilimit2;
  double delta1;

  double h1dummy1, h1dummy2;
  double h1lovalue1,h1hivalue1,h1hivalue2;

  double h2dummy1=0.0, h2dummy2;
  double h2lovalue1=0.0,h2lovalue2,h2hivalue1=0.0,h2hivalue2;

  double h3dummy1=0.0, h3dummy2;
  double h3lovalue1,h3hivalue1=0.0,h3hivalue2;

  double exparg, besselarg, expterm, bessi1overxterm, bessi0term;
  double expbetaTterm=0.0, alphasqTterm=0.0;
  double h1relval, h2relval=0.0, h3relval=0.0;
  int doh1=1, doh2=1, doh3=1;

  int i,auxindex;

  // coefflists should already have been allocated to the necessary size

  if (DEBUG_DEVICE)
  {
    if (listsize < timeindex) {
    Xyce::dout() << "[LTRA-DBG-DEV]: LTRArlcCoeffsSetup_: not enough space in coefflist" << std::endl;
  }
  }


  //
  // we assume a piecewise linear function, and we calculate the
  // coefficients using this assumption in the integration of the
  // function


  if (T == 0.0) {
    auxindex = timeindex;
  } else {

    if (curtime - T <= 0.0) {
      auxindex = 0;
    } else {
      exact = 0;
      for (i = timeindex; i>= 0; i--) {
        if (curtime - timelist[i] ==  T) {
          exact =1;
          break;
        }
        if (curtime - timelist[i] > T) break;
      }

      if (DEBUG_DEVICE)
      {
        if ((i < 0) || ((i==0) && (exact==1)))
        Xyce::dout() << "[LTRA-DBG-DEV]: LTRAcoeffSetup: i <= 0: some mistake!" << std::endl;
      }

      if (exact == 1) {
        auxindex = i-1;
      } else {
        auxindex = i;
      }
    }
  }
  // the first coefficient

  if (auxindex != 0)
  {
    lolimit1 = T;
    hilimit1 = curtime - timelist[auxindex];
    delta1 = hilimit1 - lolimit1;

    h2lovalue1 = rlcH2Func_(T,T,alpha,beta);
    besselarg = (hilimit1 > T) ? alpha*sqrt(hilimit1*hilimit1-T*T):0.0;
    exparg = -beta*hilimit1;
    expterm = exp(exparg);
    bessi1overxterm = Xyce::Util::besselI1xOverX(besselarg);
    alphasqTterm = alpha*alpha*T;
    h2hivalue1 = // LTRArlcH2Func(hilimit1,T,alpha,beta);
      ((alpha == 0.0) || (hilimit1 < T)) ? 0.0: alphasqTterm*expterm*bessi1overxterm;

    h2dummy1 = twiceintlinfunc_(lolimit1,hilimit1,lolimit1,h2lovalue1,
      h2hivalue1,lolimit1,hilimit1)/delta1;
    h2firstcoeff = h2dummy1;
    h2relval = fabs(reltol*h2dummy1);

    h3lovalue1 = 0.0; // E3dash should be consistent with this
    bessi0term = Xyce::Util::besselI0(besselarg);
    expbetaTterm = exp(-beta*T);
    h3hivalue1 = // LTRArlcH3dashIntFunc(hilimit1,T,beta);
    ((hilimit1 <= T) || (beta == 0.0)) ? 0.0: expterm* bessi0term-expbetaTterm;
    h3dummy1 = intlinfunc_(lolimit1,hilimit1,h3lovalue1,
      h3hivalue1,lolimit1,hilimit1)/delta1;
    h3dashfirstcoeff = h3dummy1;
    h3relval = fabs(h3dummy1*reltol);
  }
  else
  {
    h2firstcoeff = h3dashfirstcoeff = 0.0;
  }

  lolimit1 = 0.0;
  hilimit1 = curtime - timelist[timeindex];
  delta1 = hilimit1 - lolimit1;
  exparg = -beta*hilimit1;
  expterm = exp(exparg);

  h1lovalue1 = 0.0;
  h1hivalue1 = //LTRArlcH1dashTwiceIntFunc(hilimit1,beta);
    (beta == 0.0) ? hilimit1 : ((hilimit1 == 0.0) ? 0.0 :
                               (Xyce::Util::besselI1(-exparg)+Xyce::Util::besselI0(-exparg))* hilimit1 * expterm - hilimit1);
  h1dummy1 = h1hivalue1/delta1;
  h1dashfirstcoeff = h1dummy1;
  h1relval = fabs(h1dummy1*reltol);


  // the coefficients for the rest of the timepoints

  for (i=timeindex; i>0; i--)
  {
    if (doh1 || doh2 || doh3)
    {
      lolimit2 = lolimit1; // previous lolimit1
      hilimit2 = hilimit1; // previous hilimit1

      lolimit1 = hilimit2;
      hilimit1 = curtime - timelist[i - 1];
      delta1 = timelist[i] - timelist[i - 1];

      exparg = -beta*hilimit1;
      expterm = exp(exparg);
    }

    if (doh1)
    {
      h1hivalue2 = h1hivalue1; // previous hivalue1
      h1dummy2 = h1dummy1; // previous dummy1

      h1lovalue1 = h1hivalue2;
      h1hivalue1 = // LTRArlcH1dashTwiceIntFunc(hilimit1,beta);
        (beta == 0.0) ? hilimit1 : ((hilimit1 == 0.0) ? 0.0 :
                                    (Xyce::Util::besselI1(-exparg)+Xyce::Util::besselI0(-exparg))* hilimit1 * expterm - hilimit1);
      h1dummy1 = (h1hivalue1 - h1lovalue1)/delta1;

      h1dashcoeffs[i] = h1dummy1 - h1dummy2;
      if (fabs(h1dashcoeffs[i]) <= h1relval) doh1 = 0;
    }
    else
    {
      h1dashcoeffs[i] = 0.0;
    }

    if (i <= auxindex)
    {
      // if (i == auxindex) {
      // lolimit2 = T;
      // delta2 = hilimit2 - lolimit2;
      // }

      if (doh2 || doh3)
      {
        besselarg = (hilimit1 > T) ? alpha*sqrt(hilimit1*hilimit1-T*T):0.0;
      }

      if (doh2)
      {
        h2lovalue2 = h2lovalue1; // previous lovalue1
        h2hivalue2 = h2hivalue1; // previous hivalue1
        h2dummy2 = h2dummy1; // previous dummy1

        h2lovalue1 = h2hivalue2;
        bessi1overxterm = Xyce::Util::besselI1xOverX(besselarg);
        h2hivalue1 = // rlcH2Func(hilimit1,T,alpha,beta);
          ((alpha == 0.0) || (hilimit1 < T)) ? 0.0: alphasqTterm*expterm*bessi1overxterm;
        h2dummy1 = twiceintlinfunc_(lolimit1,hilimit1,lolimit1,
        h2lovalue1,h2hivalue1,lolimit1,hilimit1)/delta1;

        h2coeffs[i] = h2dummy1 - h2dummy2 + intlinfunc_(lolimit2,hilimit2,
        h2lovalue2,h2hivalue2,lolimit2,hilimit2);
        if (fabs(h2coeffs[i]) <= h2relval) doh2 = 0;
      }
      else
      {
        h2coeffs[i] = 0.0;
      }

      if (doh3)
      {
        h3hivalue2 = h3hivalue1; //previous hivalue1
        h3dummy2 = h3dummy1; // previous dummy1

        h3lovalue1 = h3hivalue2;
        bessi0term = Xyce::Util::besselI0(besselarg);
        h3hivalue1 = //LTRArlcH3dashIntFunc(hilimit1,T,beta);
        ((hilimit1 <= T) || (beta == 0.0)) ? 0.0: expterm* bessi0term-expbetaTterm;
        h3dummy1 = intlinfunc_(lolimit1,hilimit1,h3lovalue1,h3hivalue1,lolimit1,hilimit1)/delta1;

        h3dashcoeffs[i] = h3dummy1 - h3dummy2;
        if (fabs(h3dashcoeffs[i]) <= h3relval) doh3 = 0;
      }
      else
      {
        h3dashcoeffs[i] = 0.0;
      }
    }
  }
  *auxindexptr = auxindex;
}

//-----------------------------------------------------------------------------
// Function      : Model::straightLineCheck_
// Purpose       :
//
// takes the co-ordinates of three points,
// finds the area of the triangle enclosed by these points and
// compares this area with the area of the quadrilateral formed by
// the line between the first point and the third point, the
// perpendiculars from the first and third points to the x-axis, and
// the x-axis. If within reltol, then it returns 1, else 0. The
// purpose of this function is to determine if three points lie
// acceptably close to a straight line. This area criterion is used
// because it is related to integrals and convolution
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
bool Model::straightLineCheck_(double x1, double y1,
                                         double x2, double y2,
                                         double x3, double y3,
                                         double reltol, double abstol)
{
  // double asqr, bsqr, csqr, c, c1sqr;
  // double htsqr;
  double TRarea, QUADarea1,QUADarea2,QUADarea3, area;

  // asqr = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1);
  // bsqr = (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2);
  // csqr = (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1);
  // c = sqrt(csqr);
  // c1sqr = (asqr - bsqr + csqr)/(2*c);
  // c1sqr *= c1sqr;
  // htsqr = asqr - c1sqr;
  // TRarea = c*sqrt(htsqr)*0.5;

  // this should work if y1,y2,y3 all have the same sign and x1,x2,x3
  // are in increasing order

  QUADarea1 = (fabs(y2)+fabs(y1))*0.5*fabs(x2-x1);
  QUADarea2 = (fabs(y3)+fabs(y2))*0.5*fabs(x3-x2);
  QUADarea3 = (fabs(y3)+fabs(y1))*0.5*fabs(x3-x1);
  TRarea = fabs( QUADarea3 - QUADarea1 - QUADarea2);
  area = QUADarea1 + QUADarea2;
  if (area*reltol + abstol > TRarea)
    return(true);
  else
    return(false);
}

// i is the index of the latest value,
// a,b,c values correspond to values at t_{i-2}, t{i-1} and t_i
//
// ERK: Note: check curtime.
#define SECONDDERIV(i,a,b,c) \
  (oof = (i==getSolverState().ltraTimeIndex_?getSolverState().currTime_:                  \
          (getSolverState().ltraTimePoints_[i])),                                \
   (( c - b )/(oof-(getSolverState().ltraTimePoints_[i-1])) -                    \
    ( b - a )/((getSolverState().ltraTimePoints_[i-1])-                          \
               (getSolverState().ltraTimePoints_[i-2])))/(oof -                  \
                                                 (getSolverState().ltraTimePoints_[i-2])))


//-----------------------------------------------------------------------------
// Function      : Model::SECONDDERIV_
// Purpose       :
// Special Notes : see macro, above, modified from spice3
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::SECONDDERIV_(int i, double a, double b, double c)
{
  double oof=0.0;
  return SECONDDERIV(i,a,b,c);
}

//-----------------------------------------------------------------------------
// Function      : Model::lteCalculate_
// Purpose       :
//
// returns sum of the absolute values of the total
// local truncation error of the 2 equations for the LTRAline
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::lteCalculate_ (
  Instance& instance,
  double curtime
  )
{
  double h1dashTfirstCoeff;
  double h2TfirstCoeff=0.0;
  double h3dashTfirstCoeff=0.0;
  double dashdash;
  double hilimit1, lolimit1, hivalue1, lovalue1, f1i, g1i;
  double eq1LTE=0.0, eq2LTE=0.0;
  int auxindex, tdover, i, exact;

  switch(specialCase)
  {
    case LTRA_MOD_LC:
    case LTRA_MOD_RG:
      return(0.0);
      break;

    case LTRA_MOD_RLC:

      if (curtime > td)
      {
        tdover = 1;

        exact = 0;

        for (i=(getSolverState().ltraTimeIndex_-1); i >= 0; i--)
        {
          if (curtime - getSolverState().ltraTimePoints_[i] ==  td)
          {
            exact = 1;
            break;
          }
          if (curtime - getSolverState().ltraTimePoints_[i] > td)
          {
            break;
          }
        }

        if (DEBUG_DEVICE)
        {
          if ((i < 0) || ((i==0) && (exact==1)))
          Xyce::dout() << "[LTRA-DBG-DEV]: lteCalculate_: i <= 0: some mistake!" << std::endl;
        }

        if (exact == 1)
        {
          auxindex = i-1;
        }
        else
        {
          auxindex = i;
        }
      }
      else
      {
        tdover = 0;
      }

      hilimit1 = curtime - getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_-1];
      lolimit1 = 0.0;
      hivalue1 = rlcH1dashTwiceIntFunc_(hilimit1,beta);
      lovalue1 = 0.0;

      f1i = hivalue1;
      g1i = intlinfunc_(lolimit1,hilimit1,lovalue1,hivalue1,
                        lolimit1,hilimit1);
      h1dashTfirstCoeff = 0.5 * f1i *
        (curtime - getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_-1]) - g1i;

      if (tdover)
      {
        hilimit1 = curtime - getSolverState().ltraTimePoints_[auxindex];
        lolimit1 = getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_-1] - getSolverState().ltraTimePoints_[auxindex];
        lolimit1 = std::max(td,lolimit1);

        // are the following really doing the operations in the write-up?
        hivalue1 = rlcH2Func_(hilimit1,td,alpha,beta);
        lovalue1 = rlcH2Func_(lolimit1,td,alpha,beta);
        f1i = twiceintlinfunc_(lolimit1,hilimit1,lolimit1,lovalue1,hivalue1,lolimit1,
                               hilimit1);
        g1i = thriceintlinfunc_(lolimit1,hilimit1,lolimit1,lolimit1,lovalue1,
                                hivalue1,lolimit1,hilimit1);

        h2TfirstCoeff = 0.5*f1i*(curtime-td-getSolverState().ltraTimePoints_[auxindex]) - g1i;

        hivalue1 = rlcH3dashIntFunc_(hilimit1,td,beta);
        lovalue1 = rlcH3dashIntFunc_(lolimit1,td,beta);
        f1i = intlinfunc_(lolimit1,hilimit1,lovalue1,hivalue1,lolimit1,
                          hilimit1);
        g1i = twiceintlinfunc_(lolimit1,hilimit1,lolimit1,lovalue1,
                               hivalue1,lolimit1,hilimit1);
        h3dashTfirstCoeff = 0.5*f1i*(curtime-td-getSolverState().ltraTimePoints_[auxindex]) - g1i;
      }


      //  LTEs for convolution with v1
      //  get divided differences for v1 (2nd derivative estimates)

      //   no need to subtract operating point values because
      //   taking differences anyway
      //


      dashdash = SECONDDERIV_(getSolverState().ltraTimeIndex_,
                              instance.v1[getSolverState().ltraTimeIndex_-2],
                              instance.v1[getSolverState().ltraTimeIndex_-1],
                              instance.v1[getSolverState().ltraTimeIndex_]);
      eq1LTE += admit*fabs(dashdash * h1dashTfirstCoeff);

      // not bothering to interpolate since everything is approximate
      // anyway
      if (tdover)
      {
        dashdash = SECONDDERIV_(auxindex+1,
                                instance.v1[auxindex - 1],
                                instance.v1[auxindex],
                                instance.v1[auxindex + 1]) ;

        eq2LTE += admit*fabs(dashdash * h3dashTfirstCoeff);
      }
      // end LTEs for convolution with v1

      // LTEs for convolution with v2
      // get divided differences for v2 (2nd derivative estimates)

      dashdash = SECONDDERIV_(getSolverState().ltraTimeIndex_,
                              instance.v2[getSolverState().ltraTimeIndex_-2],
                              instance.v2[getSolverState().ltraTimeIndex_-1],
                              instance.v2[getSolverState().ltraTimeIndex_]);

      eq2LTE += admit*fabs(dashdash * h1dashTfirstCoeff);

      if (tdover)
      {
        dashdash = SECONDDERIV_(auxindex+1,
                                instance.v2[auxindex - 1],
                                instance.v2[auxindex],
                                instance.v2[auxindex + 1]);

        eq1LTE += admit*fabs(dashdash * h3dashTfirstCoeff);
      }

      // end LTEs for convolution with v2

      // LTE for convolution with i1
      // get divided differences for i1 (2nd derivative estimates)

      if (tdover)
      {
        dashdash = SECONDDERIV_(auxindex+1,
                                instance.i1[auxindex - 1],
                                instance.i1[auxindex],
                                instance.i1[auxindex + 1]) ;

        eq2LTE += fabs(dashdash * h2TfirstCoeff);
      }
      // end LTE for convolution with i1

      // LTE for convolution with i2
      // get divided differences for i2 (2nd derivative estimates)

      if (tdover)
      {
        dashdash = SECONDDERIV_(auxindex+1,
                                instance.i2[auxindex - 1],
                                instance.i2[auxindex],
                                instance.i2[auxindex + 1]) ;

        eq1LTE += fabs(dashdash * h2TfirstCoeff);
      }

      // end LTE for convolution with i1

      break;

    case LTRA_MOD_RC:

      hilimit1 = curtime - getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_-1];
      lolimit1 = 0.0;

      hivalue1 = rcH1dashTwiceIntFunc_(hilimit1,cByR);
      lovalue1 = 0.0;

      f1i = hivalue1;
      g1i = intlinfunc_(lolimit1,hilimit1,lovalue1,hivalue1,lolimit1,hilimit1);

      h1dashTfirstCoeff = 0.5*f1i*(curtime-getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_-1]) - g1i;

      hivalue1 = rcH2TwiceIntFunc_(hilimit1,rclsqr);
      lovalue1 = 0.0;

      f1i = hivalue1;
      g1i = intlinfunc_(lolimit1,hilimit1,lovalue1,hivalue1,lolimit1,hilimit1);
      h1dashTfirstCoeff = 0.5*f1i*(curtime-getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_-1]) - g1i;

      hivalue1 = rcH2TwiceIntFunc_(hilimit1,rclsqr);
      lovalue1 = 0.0;

      f1i = hivalue1;
      g1i = intlinfunc_(lolimit1,hilimit1,lovalue1,
                        hivalue1,lolimit1,hilimit1);
      h1dashTfirstCoeff = 0.5*f1i*(curtime-getSolverState().ltraTimePoints_[getSolverState().ltraTimeIndex_-1]) - g1i;

      // LTEs for convolution with v1
      // get divided differences for v1 (2nd derivative estimates)

      // no need to subtract operating point values because
      // taking differences anyway

      dashdash = SECONDDERIV_( getSolverState().ltraTimeIndex_,
                               instance.v1[getSolverState().ltraTimeIndex_-2],
                               instance.v1[getSolverState().ltraTimeIndex_-1],
                               instance.v1[getSolverState().ltraTimeIndex_] );

      eq1LTE += fabs(dashdash * h1dashTfirstCoeff);
      eq2LTE += fabs(dashdash * h3dashTfirstCoeff);

      // end LTEs for convolution with v1

      // LTEs for convolution with v2
      // get divided differences for v2 (2nd derivative estimates)

      dashdash = SECONDDERIV_( getSolverState().ltraTimeIndex_,
                               instance.v2[getSolverState().ltraTimeIndex_-2],
                               instance.v2[getSolverState().ltraTimeIndex_-1],
                               instance.v2[getSolverState().ltraTimeIndex_] );

      eq2LTE += fabs(dashdash * h1dashTfirstCoeff);
      eq1LTE += fabs(dashdash * h3dashTfirstCoeff);

      // end LTEs for convolution with v2

      // LTE for convolution with i1
      // get divided differences for i1 (2nd derivative estimates)

      dashdash = SECONDDERIV_( getSolverState().ltraTimeIndex_,
                               instance.i1[getSolverState().ltraTimeIndex_-2],
                               instance.i1[getSolverState().ltraTimeIndex_-1],
                               instance.i1[getSolverState().ltraTimeIndex_] );

      eq2LTE += fabs(dashdash * h2TfirstCoeff);

      // end LTE for convolution with i1

      // LTE for convolution with i2
      // get divided differences for i2 (2nd derivative estimates)

      dashdash = SECONDDERIV_( getSolverState().ltraTimeIndex_,
                               instance.i2[getSolverState().ltraTimeIndex_-2],
                               instance.i2[getSolverState().ltraTimeIndex_-1],
                               instance.i2[getSolverState().ltraTimeIndex_] );

      eq1LTE += fabs(dashdash * h2TfirstCoeff);

      // end LTE for convolution with i1

      break;

    default:
      return(1/*error*/);
  }

  if (DEBUG_DEVICE)
  {
    Xyce::dout() << "[LTRA-DBG-DEV] " << instance.getName() << ": LTE/input for Eq1 at time "
            << curtime << " is: " << eq1LTE/instance.input1 << std::endl;

  Xyce::dout() << "[LTRA-DBG-DEV] " << instance.getName() << ": LTE/input for Eq2 at time "
            << curtime << " is: " << eq2LTE/instance.input1 << std::endl;
  }

  return(fabs(eq1LTE) + fabs(eq2LTE));
}

//-----------------------------------------------------------------------------
// Function      : Instance::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
double Instance::getMaxTimeStepSize ()
{
  return model_.maxTimeStep;
}

// LTRA Master functions:

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes : load-type aware upper level version
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 05/03/2018
//-----------------------------------------------------------------------------
bool Master::updateState (double* solVec, double* staVec, double* stoVec,
                          int loadType)
{
  bool bsuccess = true;

  // This is a nonlinear device (as stated by our isLinear method), so
  // we should only load when the loader is asking for nonlinear device loads
  // in time domain.  Skip altogether in frequency domain.
  if (loadType == NONLINEAR || loadType == ALL)
  {
    bsuccess = updateState(solVec, staVec, stoVec);
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Master::updateState (double* solVec, double* staVec, double* stoVec)
{

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance& theInstance = *(*it);

    // Update current state
    theInstance.vpos1 = solVec[theInstance.li_Pos1];
    theInstance.vneg1 = solVec[theInstance.li_Neg1];

    theInstance.vpos2 = solVec[theInstance.li_Pos2];
    theInstance.vneg2 = solVec[theInstance.li_Neg2];

    theInstance.currp1 = solVec[theInstance.li_Ibr1];
    theInstance.currp2 = solVec[theInstance.li_Ibr2];

    // Initial state, generally the end result of the DC-OP calculation
    if (getSolverState().dcopFlag)
    {
      theInstance.initVolt1 = theInstance.vpos1 - theInstance.vneg1;
      theInstance.initVolt2 = theInstance.vpos2 - theInstance.vneg2;

      theInstance.initCur1 = theInstance.currp1;
      theInstance.initCur2 = theInstance.currp2;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes : load-type aware upper level version
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 05/03/2018
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,
                             double * bVec, double * leadF, double * leadQ,
                             double * junctionV, int loadType)
{
  bool bsuccess=true;

  // This is a nonlinear device (as stated by our isLinear method), so
  // we should only load when the loader is asking for nonlinear device loads
  // in time domain.  Skip altogether in frequency domain.
  if (loadType == NONLINEAR || loadType == ALL)
  {
    bsuccess = loadDAEVectors(solVec, fVec, qVec, bVec, leadF, leadQ,
                              junctionV);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, double * leadF, double * leadQ, double * junctionV)
{
  double max(0.0),min(0.0);
  double v1d(0.0), v2d(0.0), i1d(0.0), i2d(0.0);
  double dummy1(0.0), dummy2(0.0);
  std::ostringstream msg;

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance& theInstance = *(*it);

    if( getSolverState().dcopFlag || (theInstance.getModel().specialCase == LTRA_MOD_RG))
    {
      switch (theInstance.getModel().specialCase)
      {
        case LTRA_MOD_RG:
          dummy1 = theInstance.getModel().length *
            std::sqrt(theInstance.getModel().resist *
                      theInstance.getModel().conduct);
          dummy2 = exp(-dummy1);
          dummy1 = exp(dummy1);  // May overflow
          theInstance.getModel().coshlrootGR = 0.5 * (dummy1 + dummy2);

          if (theInstance.getModel().conduct <= 1.0e-10)
          {	// Spice3 hack!
            theInstance.getModel().rRsLrGRorG = theInstance.getModel().length *
              theInstance.getModel().resist;
          }
          else
          {
            theInstance.getModel().rRsLrGRorG =
              0.5 * (dummy1 - dummy2) * sqrt(theInstance.getModel().resist /
                                             theInstance.getModel().conduct);
          }

          if (theInstance.getModel().resist <= 1.0e-10)
          {	// Spice3 hack!
            theInstance.getModel().rGsLrGRorR = theInstance.getModel().length *
	      theInstance.getModel().conduct;
          }
          else
          {
            theInstance.getModel().rGsLrGRorR =
              0.5 * (dummy1 - dummy2) * sqrt(theInstance.getModel().conduct /
                                             theInstance.getModel().resist);
          }

          fVec[theInstance.li_Ibr1] += (theInstance.vpos1 -
                               theInstance.vneg1 -
                               theInstance.getModel().coshlrootGR *
                                        theInstance.vpos2 +
                               theInstance.getModel().coshlrootGR *
                                        theInstance.vneg2 +
                               (1.0 + getDeviceOptions().gmin) *
                                        theInstance.getModel().rRsLrGRorG *
                                        theInstance.currp2);

          fVec[theInstance.li_Ibr2] += (theInstance.getModel().coshlrootGR *
                                        theInstance.currp2 -
                               (1.0 + getDeviceOptions().gmin) *
                                        theInstance.getModel().rGsLrGRorR *
                                        theInstance.vpos2 +
                               (1.0 + getDeviceOptions().gmin) *
                                        theInstance.getModel().rGsLrGRorR *
                                        theInstance.vneg2 +
                               theInstance.currp1);

          break;

        // load a simple resistor (DC case, C=open, L=short). In the
        // lossless case (R=0.0) the port voltages are equal.
        case LTRA_MOD_LC:
        case LTRA_MOD_RC:
        case LTRA_MOD_RLC:

          // i_1 + i_2 = 0
          fVec[theInstance.li_Ibr1] += (theInstance.currp1 + theInstance.currp2);

          // v_(n1+) - v_(n2+) - R_ltra * i_1 = 0
          fVec[theInstance.li_Ibr2] += (theInstance.vpos1 - theInstance.vpos2 -
                               theInstance.currp1 *
                               theInstance.getModel().resist *
                               theInstance.getModel().length);

          break;

        default:
          UserError(theInstance) << "Unknown LTRA configuration.  Must be one of RG, LC, RC, or RLC.";
          return false;
      }

      // These are common for all DC cases. They are just the residuals
      // that enforce that the current out of the positive terminal of
      // the TL is equal to the current in to the negative terminal at
      // the same end of the TL.
      fVec[theInstance.li_Pos1] += theInstance.currp1;  // i_(n1+) = i_1
      fVec[theInstance.li_Neg1] += -theInstance.currp1; // i_(n1-) = -i_1

      fVec[theInstance.li_Pos2] += theInstance.currp2;  // i_(n2+) = i_2
      fVec[theInstance.li_Neg2] += -theInstance.currp2; // i_(n2-) = -i_2

    }
    else
    {
      // all cases other than DC or the RG case

      int isaved = 0;
      double qf1, qf2, qf3;
      double lf2, lf3;

      qf1 = qf2 = qf3 = 0.0;
      lf2 = lf3 = 0.0;

      theInstance.getModel().modelCalculations_(isaved, qf1, qf2, qf3, lf2, lf3);

      theInstance.input1 = theInstance.input2 = 0.0;

      switch (theInstance.getModel().specialCase)
      {
        case LTRA_MOD_LC:
        case LTRA_MOD_RLC:

          if (theInstance.getModel().tdover)
          {
            // have to interpolate values
            if ((isaved != 0) &&
                ((theInstance.getModel().howToInterp == LTRA_MOD_QUADINTERP) || (theInstance.getModel().howToInterp == LTRA_MOD_MIXEDINTERP)))
            {
              v1d = theInstance.v1[isaved-1] * qf1
                + theInstance.v1[isaved]   * qf2
                + theInstance.v1[isaved+1] * qf3;

              max = std::max(theInstance.v1[isaved-1], theInstance.v1[isaved]);
              max = std::max(max,theInstance.v1[isaved+1]);
              min = std::min(theInstance.v1[isaved-1], theInstance.v1[isaved]);
              min = std::min(min,theInstance.v1[isaved+1]);
            }

            if ((theInstance.getModel().howToInterp == LTRA_MOD_LININTERP) ||
                (isaved == 0) ||
                ((isaved != 0) &&
                 ((theInstance.getModel().howToInterp == LTRA_MOD_QUADINTERP) ||
                  (theInstance.getModel().howToInterp == LTRA_MOD_MIXEDINTERP)) &&
                 ((v1d > max) || (v1d < min))))
            {
              if ((isaved != 0) &&
                  (theInstance.getModel().howToInterp == LTRA_MOD_QUADINTERP))
              {
                if (DEBUG_DEVICE)
                {
                  Xyce::dout() << "[LTRA-DBG-DEV] load: warning: interpolated v1 is out of range after timepoint "
                          << getSolverState().ltraTimeIndex_ << std::endl;
                Xyce::dout() << "         values: "
                          << theInstance.v1[isaved-1] << "  "
                          << theInstance.v1[isaved] << "  "
                          << theInstance.v1[isaved+1] << "; interpolated: "
                          << v1d << std::endl;
                Xyce::dout() << "        timepoints are: "
                          << getSolverState().currTime_ - theInstance.getModel().td << std::endl;
                }
              }
              else
              {
                v1d = theInstance.v1[isaved] * lf2 + theInstance.v1[isaved+1] * lf3;
              }
            }

            if ((isaved != 0) &&
                ((theInstance.getModel().howToInterp == LTRA_MOD_QUADINTERP) ||
                 (theInstance.getModel().howToInterp == LTRA_MOD_MIXEDINTERP)))
            {
              i1d = theInstance.i1[isaved-1] * qf1
                + theInstance.i1[isaved] * qf2
                + theInstance.i1[isaved+1] * qf3;

              max = std::max(theInstance.i1[isaved-1], theInstance.i1[isaved]);
              max = std::max(max,theInstance.i1[isaved+1]);
              min = std::min(theInstance.i1[isaved-1], theInstance.i1[isaved]);
              min = std::min(min,theInstance.i1[isaved+1]);
            }

            if ((theInstance.getModel().howToInterp == LTRA_MOD_LININTERP) ||
                (isaved == 0) ||
                ((isaved != 0) &&
                 ((theInstance.getModel().howToInterp == LTRA_MOD_QUADINTERP) ||
                  (theInstance.getModel().howToInterp == LTRA_MOD_MIXEDINTERP)) &&
                 ((i1d > max) || (i1d < min))))
            {

              if ((isaved != 0) &&
                  (theInstance.getModel().howToInterp == LTRA_MOD_QUADINTERP))
              {
                if (DEBUG_DEVICE)
                {
                  Xyce::dout() << "[LTRA-DBG-DEV] load: warning: interpolated i1 is out of range after timepoint "
                          << getSolverState().ltraTimeIndex_ << std::endl;
                Xyce::dout() << "         values: "
                          << theInstance.i1[isaved-1] << "  "
                          << theInstance.i1[isaved] << "  "
                          << theInstance.i1[isaved+1] << "; interpolated: "
                          << i1d << std::endl;
                Xyce::dout() << "        timepoints are: "
                          << getSolverState().currTime_ - theInstance.getModel().td << std::endl;
                }
              }
              else
              {
                i1d = theInstance.i1[isaved] * lf2 + theInstance.i1[isaved+1] *
                  lf3;
              }
            }

            if ((isaved != 0) &&
                ((theInstance.getModel().howToInterp == LTRA_MOD_QUADINTERP) ||
                 (theInstance.getModel().howToInterp == LTRA_MOD_MIXEDINTERP)))
            {
              v2d = theInstance.v2[isaved-1] * qf1
                + theInstance.v2[isaved] * qf2
                + theInstance.v2[isaved+1] * qf3;

              max = std::max(theInstance.v2[isaved-1], theInstance.v2[isaved]);
              max = std::max(max,theInstance.v2[isaved+1]);
              min = std::min(theInstance.v2[isaved-1], theInstance.v2[isaved]);
              min = std::min(min,theInstance.v2[isaved+1]);
            }

            if ((theInstance.getModel().howToInterp ==
                 LTRA_MOD_LININTERP) || (isaved == 0) ||
                ((isaved != 0) &&
                 ((theInstance.getModel().howToInterp == LTRA_MOD_QUADINTERP) ||
                  (theInstance.getModel().howToInterp == LTRA_MOD_MIXEDINTERP)) &&
                 ((v2d > max) || (v2d < min))))
            {

              if ((isaved != 0) &&
                  (theInstance.getModel().howToInterp == LTRA_MOD_QUADINTERP))
              {
                if (DEBUG_DEVICE)
                {
                  Xyce::dout() << "[LTRA-DBG-DEV] load: warning: interpolated v2 is out of range after timepoint "
                          << getSolverState().ltraTimeIndex_ << std::endl;
                Xyce::dout() << "         values: "
                          << theInstance.v2[isaved-1] << "  "
                          << theInstance.v2[isaved] << "  "
                          << theInstance.v2[isaved+1] << "; interpolated: "
                          << v2d << std::endl;
                Xyce::dout() << "        timepoints are: "
                          << getSolverState().currTime_ - theInstance.getModel().td << std::endl;
                }
              }
              else
              {
                v2d = theInstance.v2[isaved] * lf2
                  + theInstance.v2[isaved+1] *
                  lf3;
              }
            }

            if ((isaved != 0) &&
                ((theInstance.getModel().howToInterp == LTRA_MOD_QUADINTERP) ||
                 (theInstance.getModel().howToInterp == LTRA_MOD_MIXEDINTERP)))
            {
              i2d = theInstance.i2[isaved-1] * qf1
                + theInstance.i2[isaved] * qf2
                + theInstance.i2[isaved+1] * qf3;

              max = std::max(theInstance.i2[isaved-1], theInstance.i2[isaved]);
              max = std::max(max,theInstance.i2[isaved+1]);
              min = std::min(theInstance.i2[isaved-1], theInstance.i2[isaved]);
              min = std::min(min,theInstance.i2[isaved+1]);
            }

            if ((theInstance.getModel().howToInterp == LTRA_MOD_LININTERP) ||
                (isaved == 0) ||
                ((isaved != 0) &&
                 ((theInstance.getModel().howToInterp == LTRA_MOD_QUADINTERP) ||
                  (theInstance.getModel().howToInterp == LTRA_MOD_MIXEDINTERP)) &&
                 ((i2d > max) || (i2d < min))))
            {
              if ((isaved != 0) &&
                  (theInstance.getModel().howToInterp == LTRA_MOD_QUADINTERP))
              {
                if (DEBUG_DEVICE)
                {
                  Xyce::dout() << "[LTRA-DBG-DEV] load: warning: interpolated i2 is out of range after timepoint "
                          << getSolverState().ltraTimeIndex_ << std::endl;
                Xyce::dout() << "         values: "
                          << theInstance.i2[isaved-1] << "  "
                          << theInstance.i2[isaved] << "  "
                          << theInstance.i2[isaved+1] << "; interpolated: "
                          << i2d << std::endl;
                Xyce::dout() << "        timepoints are: "
                          << getSolverState().currTime_ - theInstance.getModel().td << std::endl;
                }
              }
              else
              {
                i2d = theInstance.i2[isaved] * lf2 + theInstance.i2[isaved+1] * lf3;
              }
            }
          }

          // interpolation done
          break;

      case LTRA_MOD_RC:
        break;

      default:
        return false;
        // return(E_BADPARM);
      }

      switch (theInstance.getModel().specialCase)
      {
        case LTRA_MOD_RLC:

          // begin convolution parts

          // convolution of h1dash with v1 and v2
          // the matrix has already been loaded above

          dummy1 = dummy2 = 0.0;
          for (int j = getSolverState().ltraTimeIndex_; j > 0; j--)
          {
            if (theInstance.getModel().h1dashCoeffs[j] != 0.0)
            {
              dummy1 += theInstance.getModel().h1dashCoeffs[j] *
                (theInstance.v1[j] - theInstance.initVolt1);
              dummy2 += theInstance.getModel().h1dashCoeffs[j] *
                (theInstance.v2[j] - theInstance.initVolt2);
            }
          }

          dummy1 += theInstance.initVolt1 * theInstance.getModel().intH1dash;
          dummy2 += theInstance.initVolt2 * theInstance.getModel().intH1dash;

          dummy1 -= theInstance.initVolt1 *
            theInstance.getModel().h1dashFirstCoeff;
          dummy2 -= theInstance.initVolt2 *
            theInstance.getModel().h1dashFirstCoeff;

          theInstance.input1 -= dummy1 * theInstance.getModel().admit;
          theInstance.input2 -= dummy2 * theInstance.getModel().admit;

          // end convolution of h1dash with v1 and v2

          // convolution of h2 with i2 and i1

          dummy1 = dummy2 = 0.0;
          if (theInstance.getModel().tdover)
          {
            // the term for ckt->CKTtime - theInstance.getModel().td
            dummy1 = (i2d - theInstance.initCur2)*
              theInstance.getModel().h2FirstCoeff;
            dummy2 = (i1d - theInstance.initCur1)*
              theInstance.getModel().h2FirstCoeff;

            // the rest of the convolution

            for (int j= theInstance.getModel().auxIndex; j > 0; j--)
            {

              if (theInstance.getModel().h2Coeffs[j] != 0.0)
              {
                dummy1 += theInstance.getModel().h2Coeffs[j] *
                  (theInstance.i2[j] - theInstance.initCur2);
                dummy2 += theInstance.getModel().h2Coeffs[j] *
                  (theInstance.i1[j] - theInstance.initCur1);
              }
            }
          }

          // the initial-condition terms

          dummy1 += theInstance.initCur2 * theInstance.getModel().intH2;
          dummy2 += theInstance.initCur1 * theInstance.getModel().intH2;

          theInstance.input1 += dummy1;
          theInstance.input2 += dummy2;

          // end convolution of h2 with i2 and i1
          // convolution of h3dash with v2 and v1
          // the term for ckt->CKTtime - theInstance.getModel().td

          dummy1 = dummy2 = 0.0;
          if (theInstance.getModel().tdover)
          {
            dummy1 = (v2d - theInstance.initVolt2)*
              theInstance.getModel().h3dashFirstCoeff;
            dummy2 = (v1d - theInstance.initVolt1)*
              theInstance.getModel().h3dashFirstCoeff;

            // the rest of the convolution

            for (int j= theInstance.getModel().auxIndex; j > 0; j--)
            {
              if (theInstance.getModel().h3dashCoeffs[j] != 0.0)
              {
                dummy1 += theInstance.getModel().h3dashCoeffs[j] *
                  (theInstance.v2[j] - theInstance.initVolt2);
                dummy2 += theInstance.getModel().h3dashCoeffs[j] *
                  (theInstance.v1[j] - theInstance.initVolt1);
              }
            }
          }

          // the initial-condition terms

          dummy1 += theInstance.initVolt2 * theInstance.getModel().intH3dash;
          dummy2 += theInstance.initVolt1 * theInstance.getModel().intH3dash;

          theInstance.input1 += theInstance.getModel().admit*dummy1;
          theInstance.input2 += theInstance.getModel().admit*dummy2;

          // end convolution of h3dash with v2 and v1

          // NOTE: this switch passes through to following case

        case LTRA_MOD_LC:
          // begin lossless-like parts

          if (!theInstance.getModel().tdover)
          {
            theInstance.input1 += theInstance.getModel().attenuation *
              (theInstance.initVolt2*theInstance.getModel().admit +
               theInstance.initCur2);
            theInstance.input2 += theInstance.getModel().attenuation *
              (theInstance.initVolt1*theInstance.getModel().admit +
               theInstance.initCur1);
          }
          else
          {
            theInstance.input1 += theInstance.getModel().attenuation *
              (v2d*theInstance.getModel().admit + i2d);
            theInstance.input2 += theInstance.getModel().attenuation *
              (v1d*theInstance.getModel().admit + i1d);
          }

          // Residuals for the internal equations. These are for both
          // the RLC and LC case.
          fVec[theInstance.li_Ibr1] +=
            ((theInstance.getModel().admit *
              (theInstance.getModel().h1dashFirstCoeff + 1.0)) *
             (theInstance.vpos1-theInstance.vneg1) - theInstance.currp1) -
            theInstance.input1;

          fVec[theInstance.li_Ibr2] +=
            ((theInstance.getModel().admit *
              (theInstance.getModel().h1dashFirstCoeff + 1.0)) *
             (theInstance.vpos2-theInstance.vneg2) - theInstance.currp2) -
            theInstance.input2;

          // end lossless-like parts
          break;

        case LTRA_MOD_RC:

          // begin convolution parts

          // convolution of h1dash with v1 and v2
          // the matrix has already been loaded above

          dummy1 = 0.0;
          dummy2 = 0.0;
          for (int j = getSolverState().ltraTimeIndex_; j > 0; j--)
          {
            if (theInstance.getModel().h1dashCoeffs[j]!= 0.0)
            {
              dummy1 += theInstance.getModel().h1dashCoeffs[j] *
                (theInstance.v1[j] - theInstance.initVolt1);
              dummy2 += theInstance.getModel().h1dashCoeffs[j] *
                (theInstance.v2[j] - theInstance.initVolt2);
            }
          }

          // the initial condition terms

          dummy1 += theInstance.initVolt1 * theInstance.getModel().intH1dash;
          dummy2 += theInstance.initVolt2 * theInstance.getModel().intH1dash;

          // the constant contributed by the init
          // condition and the latest timepoint

          dummy1 -= theInstance.initVolt1*
            theInstance.getModel().h1dashFirstCoeff;
          dummy2 -= theInstance.initVolt2*
            theInstance.getModel().h1dashFirstCoeff;

          theInstance.input1 -= dummy1;
          theInstance.input2 -= dummy2;

          // end convolution of h1dash with v1 and v2
          // convolution of h2 with i2 and i1

          dummy1=dummy2=0.0;

          for (int j = getSolverState().ltraTimeIndex_; j > 0; j--)
          {
            if (theInstance.getModel().h2Coeffs[j] != 0.0)
            {
              dummy1 += theInstance.getModel().h2Coeffs[j] *
                (theInstance.i2[j] - theInstance.initCur2);
              dummy2 += theInstance.getModel().h2Coeffs[j] *
                (theInstance.i1[j] - theInstance.initCur1);
            }
          }

          // the initial-condition terms
          dummy1 += theInstance.initCur2 * theInstance.getModel().intH2;
          dummy2 += theInstance.initCur1 * theInstance.getModel().intH2;

          dummy1 -= theInstance.initCur2* theInstance.getModel().h2FirstCoeff;
          dummy2 -= theInstance.initCur1* theInstance.getModel().h2FirstCoeff;

          theInstance.input1 += dummy1;
          theInstance.input2 += dummy2;

          // end convolution of h2 with i2 and i1
          // convolution of h3dash with v2 and v1

          dummy1 = dummy2 = 0.0;

          for (int j=getSolverState().ltraTimeIndex_; j > 0; j--)
          {
            if (theInstance.getModel().h3dashCoeffs[j] != 0.0)
            {
              dummy1 += theInstance.getModel().h3dashCoeffs[j] *
                (theInstance.v2[j] - theInstance.initVolt2);
              dummy2 += theInstance.getModel().h3dashCoeffs[j] *
                (theInstance.v1[j] - theInstance.initVolt1);
            }
          }

          // the initial-condition terms

          dummy1 += theInstance.initVolt2 * theInstance.getModel().intH3dash;
          dummy2 += theInstance.initVolt1 * theInstance.getModel().intH3dash;

          dummy1 -= theInstance.initVolt2*
            theInstance.getModel().h3dashFirstCoeff;
          dummy2 -= theInstance.initVolt1*
            theInstance.getModel().h3dashFirstCoeff;

          theInstance.input1 += dummy1;
          theInstance.input2 += dummy2;

          // Residuales for the internal equations.
          fVec[theInstance.li_Ibr1] +=
            ((theInstance.getModel().h1dashFirstCoeff *
              (theInstance.vpos1-theInstance.vneg1) -
              theInstance.getModel().h3dashFirstCoeff *
              (theInstance.vpos2-theInstance.vneg2) -
              theInstance.getModel().h2FirstCoeff * theInstance.currp2 -
              theInstance.currp1) - theInstance.input1);

          fVec[theInstance.li_Ibr2] +=
            ((theInstance.getModel().h1dashFirstCoeff *
              (theInstance.vpos2-theInstance.vneg2) -
              theInstance.getModel().h3dashFirstCoeff *
              (theInstance.vpos1-theInstance.vneg1) -
              theInstance.getModel().h2FirstCoeff * theInstance.currp1 -
              theInstance.currp2) - theInstance.input2);

          // end convolution of h3dash with v2 and v1

          break;

        default:
          return false;
            //return(E_BADPARM);
      }

      // Residuals (KCL) for "normal" nodes and common between all cases
      fVec[theInstance.li_Pos1] += theInstance.currp1;
      fVec[theInstance.li_Neg1] -= theInstance.currp1;

      fVec[theInstance.li_Pos2] += theInstance.currp2;
      fVec[theInstance.li_Neg2] -= theInstance.currp2;

    }
  }

  return true;
}
//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes : load-type aware upper level version
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 05/03/2018
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix& dFdx, Linear::Matrix& dQdx,
                              int loadType)
{
  bool bsuccess=true;

  // This is a nonlinear device (as stated by our isLinear method), so
  // we should only load when the loader is asking for nonlinear device loads
  // in time domain.  Skip altogether in frequency domain.
  if (loadType == NONLINEAR || loadType == ALL)
  {
    bsuccess = loadDAEMatrices(dFdx, dQdx);
  }
  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix& dFdx, Linear::Matrix& dQdx)
{
  double dummy1(0.0);
  std::ostringstream msg;

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance& theInstance = *(*it);

    if( getSolverState().dcopFlag || (theInstance.getModel().specialCase == LTRA_MOD_RG))
    {
      switch (theInstance.getModel().specialCase)
      {
        case LTRA_MOD_RG:
          *(theInstance.ibr1Pos1Ptr) +=  1.0;
          *(theInstance.ibr1Neg1Ptr) += -1.0;
          *(theInstance.ibr1Pos2Ptr) += -theInstance.getModel().coshlrootGR;
          *(theInstance.ibr1Neg2Ptr) +=  theInstance.getModel().coshlrootGR;
          *(theInstance.ibr1Ibr2Ptr) +=  (1.0 + getDeviceOptions().gmin) *
            theInstance.getModel().rRsLrGRorG;

          *(theInstance.ibr2Ibr2Ptr) +=  theInstance.getModel().coshlrootGR;
          *(theInstance.ibr2Pos2Ptr) += -(1.0 + getDeviceOptions().gmin) *
            theInstance.getModel().rGsLrGRorR;
          *(theInstance.ibr2Neg2Ptr) +=  (1.0 + getDeviceOptions().gmin) *
            theInstance.getModel().rGsLrGRorR;
          *(theInstance.ibr2Ibr1Ptr) +=  1.0;

          *(theInstance.pos1Ibr1Ptr) +=  1.0;
          *(theInstance.neg1Ibr1Ptr) += -1.0;
          *(theInstance.pos2Ibr2Ptr) +=  1.0;
          *(theInstance.neg2Ibr2Ptr) += -1.0;

          break;

        case LTRA_MOD_LC:
        case LTRA_MOD_RLC:
        case LTRA_MOD_RC: // load a simple resistor

          *(theInstance.pos1Ibr1Ptr) +=  1.0;
          *(theInstance.neg1Ibr1Ptr) += -1.0;

          *(theInstance.pos2Ibr2Ptr) +=  1.0;
          *(theInstance.neg2Ibr2Ptr) += -1.0;

          *(theInstance.ibr1Ibr1Ptr) +=  1.0;
          *(theInstance.ibr1Ibr2Ptr) +=  1.0;

          *(theInstance.ibr2Pos1Ptr) +=  1.0;
          *(theInstance.ibr2Pos2Ptr) += -1.0;
          *(theInstance.ibr2Ibr1Ptr) += -theInstance.getModel().resist*
            theInstance.getModel().length;

          break;

        default:
          UserError(theInstance) << "Unknown LTRA configuration, " << theInstance.getModel().specialCase << ". Must be one of RG, LC, RC, or RLC.";

          return false;
      }

    }
    else
    {
      // all cases other than DC or the RG case

      // matrix loading - done every time load is called
      switch (theInstance.getModel().specialCase)
      {
        case LTRA_MOD_RLC:
          // loading for convolution parts' first terms

          dummy1 = theInstance.getModel().admit *
            theInstance.getModel().h1dashFirstCoeff;

          *(theInstance.ibr1Pos1Ptr) += dummy1;
          *(theInstance.ibr1Neg1Ptr) -= dummy1;

          *(theInstance.ibr2Pos2Ptr) += dummy1;
          *(theInstance.ibr2Neg2Ptr) -= dummy1;
          // end loading for convolution parts' first terms

          // NOTE: This case intentionally falls through to the next case

        case LTRA_MOD_LC:
          // this section loads for the parts of the equations that
          // resemble the lossless equations

          *(theInstance.ibr1Pos1Ptr) += theInstance.getModel().admit;
          *(theInstance.ibr1Neg1Ptr) -= theInstance.getModel().admit;

          *(theInstance.ibr1Ibr1Ptr) -= 1.0;

          *(theInstance.pos1Ibr1Ptr) += 1.0;
          *(theInstance.neg1Ibr1Ptr) -= 1.0;

          *(theInstance.ibr2Pos2Ptr) += theInstance.getModel().admit;
          *(theInstance.ibr2Neg2Ptr) -= theInstance.getModel().admit;

          *(theInstance.ibr2Ibr2Ptr) -= 1.0;

          *(theInstance.pos2Ibr2Ptr) += 1.0;
          *(theInstance.neg2Ibr2Ptr) -= 1.0;

          // loading for lossless-like parts over
          break;

        case LTRA_MOD_RC:

          // this section loads for the parts of the equations that
          // have no convolution 

          *(theInstance.ibr1Ibr1Ptr) -= 1.0;

          *(theInstance.pos1Ibr1Ptr) += 1.0;
          *(theInstance.neg1Ibr1Ptr) -= 1.0;

          *(theInstance.ibr2Ibr2Ptr) -= 1.0;

          *(theInstance.pos2Ibr2Ptr) += 1.0;
          *(theInstance.neg2Ibr2Ptr) -= 1.0;

          // loading for non-convolution parts over
          // loading for convolution parts' first terms

          dummy1 = theInstance.getModel().h1dashFirstCoeff;

          *(theInstance.ibr1Pos1Ptr) += dummy1;
          *(theInstance.ibr1Neg1Ptr) -= dummy1;

          *(theInstance.ibr2Pos2Ptr) += dummy1;
          *(theInstance.ibr2Neg2Ptr) -= dummy1;

          dummy1 = theInstance.getModel().h2FirstCoeff;

          *(theInstance.ibr1Ibr2Ptr) -= dummy1;
          *(theInstance.ibr2Ibr1Ptr) -= dummy1;

          dummy1 = theInstance.getModel().h3dashFirstCoeff;

          *(theInstance.ibr1Pos2Ptr) -= dummy1;
          *(theInstance.ibr1Neg2Ptr) += dummy1;

          *(theInstance.ibr2Pos1Ptr) -= dummy1;
          *(theInstance.ibr2Neg1Ptr) += dummy1;

          // end loading for convolution parts' first terms

          break;

        default:
          return false;
      }
    }
  }

  return true;
}
//-----------------------------------------------------------------------------
// Function      : Master::loadFreqDAEVectors
// Purpose       : Load frequency domain DAE vectors
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 3 May 2018
//-----------------------------------------------------------------------------
bool Master::loadFreqDAEVectors(double frequency, std::complex<double>* solVec,
                                  std::vector<Util::FreqVecEntry>& fVec,
                                  std::vector<Util::FreqVecEntry>& bVec)
{

  InstanceVector::const_iterator it, end;

  double omega=2*M_PI*frequency;

  it = getInstanceBegin();
  end = getInstanceEnd();

  fVec.clear();
  bVec.clear();

  for ( ; it != end; ++it )
  {
    Instance & theInstance = *(*it);
    Model & theModel = theInstance.getModel();
    Util::FreqVecEntry tmpEntry;
    std::complex<double> y0,lambda;

    // This should be *exactly* the same as the corresponding block
    // in loadFreqDAEMatrices
    // -- begin cut/paste
    // For RLC, LC, and RC, the Y0 and Lambda are just:
    // Y0 = sqrt( (G+i*omega*C)/(R+i*omega*L))
    // lambda = sqrt( (G+i*omega*C)*(R+i*omega*L))
    // G is always 0.
    // SPICE  has to do this all in real form, so it looks more
    // complicated than this.
    switch (theModel.specialCase)
    {
    case LTRA_MOD_LC:
      {
        std::complex<double> capReact(0,omega*theModel.capac);
        std::complex<double> indReact(0,omega*theModel.induct);
        // Don't do admit=sqrt(capReact/capInduct) here, because we
        // can be called with omega=0!
        y0 = theModel.admit;
        lambda = sqrt(capReact*indReact);
      }
      break;
    case LTRA_MOD_RLC:
      {
        std::complex<double> capReact(0,omega*theModel.capac);
        std::complex<double> indReact(0,omega*theModel.induct);

        y0 = sqrt(capReact/(theModel.resist+indReact));
        lambda = sqrt(capReact*(theModel.resist+indReact));
      }
      break;
    case LTRA_MOD_RC:
      {
        double temp;
        temp = sqrt(0.5*omega*theModel.cByR);
        y0 = std::complex<double>(temp,temp);
        temp = sqrt(0.5*omega*theModel.resist*theModel.capac);
        lambda = std::complex<double>(temp,temp);
      }
      break;
    default:
      // do nothing for RG
      break;
    }
    // -- end cut/paste

    // Now we must simply construct the vectors as dFdx*X:
    // All cases except RG are the same, but using different Y0 and Lambda
    // RG is special.
    switch(theModel.specialCase)
    {
    case LTRA_MOD_RC:
    case LTRA_MOD_LC:
    case LTRA_MOD_RLC:
      {
        std::complex<double> explambda = exp(-lambda*theModel.length);
        std::complex<double> y0exp = y0*explambda;
        // Pos1
        tmpEntry.val = solVec[theInstance.li_Ibr1];
        tmpEntry.lid = theInstance.li_Pos1;
        fVec.push_back(tmpEntry);
        // Neg1
        tmpEntry.val = -solVec[theInstance.li_Ibr1];
        tmpEntry.lid = theInstance.li_Neg1;
        fVec.push_back(tmpEntry);
        // Pos2
        tmpEntry.val = solVec[theInstance.li_Ibr2];
        tmpEntry.lid = theInstance.li_Pos2;
        fVec.push_back(tmpEntry);
        // Neg2
        tmpEntry.val = -solVec[theInstance.li_Ibr2];
        tmpEntry.lid = theInstance.li_Neg2;
        fVec.push_back(tmpEntry);

        // Ibr1
        tmpEntry.val =
          y0*(solVec[theInstance.li_Pos1]-solVec[theInstance.li_Neg1])
          - y0exp*(solVec[theInstance.li_Pos2]-solVec[theInstance.li_Neg2])
          - solVec[theInstance.li_Ibr1]
          - explambda*solVec[theInstance.li_Ibr2];
        tmpEntry.lid = theInstance.li_Ibr1;
        fVec.push_back(tmpEntry);
        // Ibr2
        tmpEntry.val =
          y0*(solVec[theInstance.li_Pos2]-solVec[theInstance.li_Neg2])
          - y0exp*(solVec[theInstance.li_Pos1]-solVec[theInstance.li_Neg1])
          - solVec[theInstance.li_Ibr2]
          - explambda*solVec[theInstance.li_Ibr1];
        tmpEntry.lid = theInstance.li_Ibr2;
        fVec.push_back(tmpEntry);
      }
      break;
    case LTRA_MOD_RG:
      // Pos1
      tmpEntry.val = solVec[theInstance.li_Ibr1];
      tmpEntry.lid = theInstance.li_Pos1;
      fVec.push_back(tmpEntry);
      // Neg1
      tmpEntry.val = -solVec[theInstance.li_Ibr1];
      tmpEntry.lid = theInstance.li_Neg1;
      fVec.push_back(tmpEntry);
      // Pos2
      tmpEntry.val = solVec[theInstance.li_Ibr2];
      tmpEntry.lid = theInstance.li_Pos2;
      fVec.push_back(tmpEntry);
      // Neg2
      tmpEntry.val = -solVec[theInstance.li_Ibr2];
      tmpEntry.lid = theInstance.li_Neg2;
      fVec.push_back(tmpEntry);

      // Ibr1
      tmpEntry.val =
        (solVec[theInstance.li_Pos1] - solVec[theInstance.li_Neg1]) -
        (theModel.coshlrootGR *
         (solVec[theInstance.li_Neg2]-solVec[theInstance.li_Pos2])) +
        ((1.0+getDeviceOptions().gmin)*theModel.rRsLrGRorG*
         solVec[theInstance.li_Ibr2]);
      tmpEntry.lid = theInstance.li_Ibr1;
      fVec.push_back(tmpEntry);

      // Ibr2
      tmpEntry.val =
        theModel.coshlrootGR*solVec[theInstance.li_Ibr2] -
        ((1.0+getDeviceOptions().gmin)*theModel.rGsLrGRorR*
         (solVec[theInstance.li_Neg2]-solVec[theInstance.li_Pos2])) +
        solVec[theInstance.li_Ibr1];
      tmpEntry.lid = theInstance.li_Ibr2;
      fVec.push_back(tmpEntry);

      break;
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadFreqDAEMatrices
// Purpose       : Load frequency domain DAE matrices
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 3 May 2018
//-----------------------------------------------------------------------------
/// load Frequency Domain matrices for the LTRA device
///
/// The FD matrix for the LTRA device in Xyce is exactly the same as the
/// small-signal AC matrix from the SPICE3F5 LTRA device, which has the
/// following comments:
///
/// loads for LTRA lines for the s.s. ac case
/// the equations are the following:
///
/// Y_0(s) * V_1(s) - I_1(s) = exp(-lambda(s)*length) * (Y_0(s) * V_2(s) + I_2(s))
/// Y_0(s) * V_2(s) - I_2(s) = exp(-lambda(s)*length) * (Y_0(s) * V_1(s) + I_1(s))
///
/// where Y_0(s) and lambda(s) are as follows:
///
/// Y_0(s) = sqrt( (sC+G)/(sL+R) )
/// lambda(s) = sqrt( (sC+G)*(sL+R) )
///
/// for the RC, RLC, and LC cases, G=0. The RG case is handled
/// exactly as the DC case, (and the above equations require
/// reformulation because they become identical for the DC case.)
bool Master::loadFreqDAEMatrices(double frequency, std::complex<double>* solVec,
                                   std::vector<Util::FreqMatEntry>& dFdx)
{
  InstanceVector::const_iterator it, end;

  double omega=2*M_PI*frequency;

  it = getInstanceBegin();
  end = getInstanceEnd();

  dFdx.clear();

  for ( ; it != end; ++it )
  {
    Instance & theInstance = *(*it);
    Model & theModel = theInstance.getModel();

    Util::FreqMatEntry tmpEntry;
    std::complex<double> y0,lambda;


    // For RLC, LC, and RC, the Y0 and Lambda are just:
    // Y0 = sqrt( (G+i*omega*C)/(R+i*omega*L))
    // lambda = sqrt( (G+i*omega*C)*(R+i*omega*L))
    // G is always 0.
    // SPICE  has to do this all in real form, so it looks more
    // complicated than this.
    switch (theModel.specialCase)
    {
    case LTRA_MOD_LC:
      {
        std::complex<double> capReact(0,omega*theModel.capac);
        std::complex<double> indReact(0,omega*theModel.induct);
        // Don't do admit=sqrt(capReact/capInduct) here, because we
        // can be called with omega=0!
        y0 = theModel.admit;
        lambda = sqrt(capReact*indReact);
      }
      break;
    case LTRA_MOD_RLC:
      {
        std::complex<double> capReact(0,omega*theModel.capac);
        std::complex<double> indReact(0,omega*theModel.induct);

        y0 = sqrt(capReact/(theModel.resist+indReact));
        lambda = sqrt(capReact*(theModel.resist+indReact));
      }
      break;
    case LTRA_MOD_RC:
      {
        double temp;
        temp = sqrt(0.5*omega*theModel.cByR);
        y0 = std::complex<double>(temp,temp);
        temp = sqrt(0.5*omega*theModel.resist*theModel.capac);
        lambda = std::complex<double>(temp,temp);
      }
      break;
    default:
      // do nothing for RG
      break;
    }

    // Now we can actually do the load
    // The LTRA matrix for the RLC, LC, and RC cases is exactly the
    // same as the SPICE3F5 AC load.
    switch(theModel.specialCase)
    {
    case LTRA_MOD_RC:
    case LTRA_MOD_LC:
    case LTRA_MOD_RLC:
      {
        std::complex<double> explambda = exp(-lambda*theModel.length);
        std::complex<double> y0exp = y0*explambda;

        // Ibr1 row
        tmpEntry.val = y0;
        tmpEntry.row_lid = theInstance.li_Ibr1;
        tmpEntry.col_lid = theInstance.AIbr1EquPos1NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = -y0;
        tmpEntry.row_lid = theInstance.li_Ibr1;
        tmpEntry.col_lid = theInstance.AIbr1EquNeg1NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = std::complex<double>(-1.0, 0.0);
        tmpEntry.row_lid = theInstance.li_Ibr1;
        tmpEntry.col_lid = theInstance.AIbr1EquIbr1NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = -y0exp;
        tmpEntry.row_lid = theInstance.li_Ibr1;
        tmpEntry.col_lid = theInstance.AIbr1EquPos2NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = y0exp;
        tmpEntry.row_lid = theInstance.li_Ibr1;
        tmpEntry.col_lid = theInstance.AIbr1EquNeg2NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = -explambda;
        tmpEntry.row_lid = theInstance.li_Ibr1;
        tmpEntry.col_lid = theInstance.AIbr1EquIbr2NodeOffset;
        dFdx.push_back(tmpEntry);

        // Ibr2 row
        tmpEntry.val = y0;
        tmpEntry.row_lid = theInstance.li_Ibr2;
        tmpEntry.col_lid = theInstance.AIbr2EquPos2NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = -y0;
        tmpEntry.row_lid = theInstance.li_Ibr2;
        tmpEntry.col_lid = theInstance.AIbr2EquNeg2NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = std::complex<double>(-1.0, 0.0);
        tmpEntry.row_lid = theInstance.li_Ibr2;
        tmpEntry.col_lid = theInstance.AIbr2EquIbr2NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = -y0exp;
        tmpEntry.row_lid = theInstance.li_Ibr2;
        tmpEntry.col_lid = theInstance.AIbr2EquPos1NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = y0exp;
        tmpEntry.row_lid = theInstance.li_Ibr2;
        tmpEntry.col_lid = theInstance.AIbr2EquNeg1NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = -explambda;
        tmpEntry.row_lid = theInstance.li_Ibr2;
        tmpEntry.col_lid = theInstance.AIbr2EquIbr1NodeOffset;
        dFdx.push_back(tmpEntry);

        // Pos1 row
        tmpEntry.val = std::complex<double>(1.0, 0.0);
        tmpEntry.row_lid = theInstance.li_Pos1;
        tmpEntry.col_lid = theInstance.APos1EquIbr1NodeOffset;
        dFdx.push_back(tmpEntry);
        // Neg1 row
        tmpEntry.val = std::complex<double>(-1.0, 0.0);
        tmpEntry.row_lid = theInstance.li_Neg1;
        tmpEntry.col_lid = theInstance.ANeg1EquIbr1NodeOffset;
        dFdx.push_back(tmpEntry);
        // Pos2 row
        tmpEntry.val = std::complex<double>(1.0, 0.0);
        tmpEntry.row_lid = theInstance.li_Pos2;
        tmpEntry.col_lid = theInstance.APos2EquIbr2NodeOffset;
        dFdx.push_back(tmpEntry);
        // Neg2 row
        tmpEntry.val = std::complex<double>(-1.0, 0.0);
        tmpEntry.row_lid = theInstance.li_Neg2;
        tmpEntry.col_lid = theInstance.ANeg2EquIbr2NodeOffset;
        dFdx.push_back(tmpEntry);
      }
      break;
    case LTRA_MOD_RG:
      {
        // Ibr1 row
        tmpEntry.val = std::complex<double>(1.0,0.0);
        tmpEntry.row_lid = theInstance.li_Ibr1;
        tmpEntry.col_lid = theInstance.AIbr1EquPos1NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = std::complex<double>(-1.0,0.0);
        tmpEntry.row_lid = theInstance.li_Ibr1;
        tmpEntry.col_lid = theInstance.AIbr1EquNeg1NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = std::complex<double>(-theModel.coshlrootGR,0.0);
        tmpEntry.row_lid = theInstance.li_Ibr1;
        tmpEntry.col_lid = theInstance.AIbr1EquPos2NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = std::complex<double>(theModel.coshlrootGR,0.0);
        tmpEntry.row_lid = theInstance.li_Ibr1;
        tmpEntry.col_lid = theInstance.AIbr1EquNeg2NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = std::complex<double>((1.0+getDeviceOptions().gmin)*
                                            theModel.rRsLrGRorG, 0.0);
        tmpEntry.row_lid = theInstance.li_Ibr1;
        tmpEntry.col_lid = theInstance.AIbr1EquIbr2NodeOffset;
        dFdx.push_back(tmpEntry);

        // Ibr2 row
        tmpEntry.val = std::complex<double>(theModel.coshlrootGR,0.0);
        tmpEntry.row_lid = theInstance.li_Ibr2;
        tmpEntry.col_lid = theInstance.AIbr2EquIbr2NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = std::complex<double>( -(1.0 + getDeviceOptions().gmin) *
                                             theModel.rGsLrGRorR,
                                             0.0);
        tmpEntry.row_lid = theInstance.li_Ibr2;
        tmpEntry.col_lid = theInstance.AIbr2EquPos2NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = std::complex<double>((1.0 + getDeviceOptions().gmin) *
                                            theModel.rGsLrGRorR,
                                            0.0);
        tmpEntry.row_lid = theInstance.li_Ibr2;
        tmpEntry.col_lid = theInstance.AIbr2EquNeg2NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = std::complex<double> (1.0,0.0);
        tmpEntry.row_lid = theInstance.li_Ibr2;
        tmpEntry.col_lid = theInstance.AIbr2EquIbr1NodeOffset;
        dFdx.push_back(tmpEntry);

        // One entry each on pos1, neg1, pos2, neg2 rows:
        tmpEntry.val = std::complex<double>(1.0,0.0);
        tmpEntry.row_lid = theInstance.li_Pos1;
        tmpEntry.col_lid = theInstance.APos1EquIbr1NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = std::complex<double>(-1.0,0.0);
        tmpEntry.row_lid = theInstance.li_Neg1;
        tmpEntry.col_lid = theInstance.ANeg1EquIbr1NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = std::complex<double>(1.0,0.0);
        tmpEntry.row_lid = theInstance.li_Pos2;
        tmpEntry.col_lid = theInstance.APos2EquIbr2NodeOffset;
        dFdx.push_back(tmpEntry);

        tmpEntry.val = std::complex<double>(-1.0,0.0);
        tmpEntry.row_lid = theInstance.li_Neg2;
        tmpEntry.col_lid = theInstance.ANeg2EquIbr2NodeOffset;
        dFdx.push_back(tmpEntry);
      }
      break;
    }
  }
  return true;
}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() || (deviceMap.find("O")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("o", 1)
      .registerModelType("ltra", 1);
  }
}

} // namespace LTRA
} // namespace Device
} // namespace Xyce

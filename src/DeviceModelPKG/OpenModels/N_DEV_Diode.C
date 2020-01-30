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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <N_UTL_Math.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_Diode.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_ExtendedString.h>
#include <N_ANP_NoiseData.h>

namespace Xyce {
namespace Device {
namespace Diode {

void Traits::loadInstanceParameters(ParametricData<Diode::Instance> &p)
{
  p.addPar ("AREA",  1.0, &Diode::Instance::Area)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Area scaling value (scales IS, ISR, IKF, RS, CJ0, and IBV)");

  p.addPar ("IC",    0.0, &Diode::Instance::InitCond)
    .setGivenMember(&Diode::Instance::InitCondGiven)
    .setCategory(CAT_NONE);

  p.addPar ("TEMP",  0.0, &Diode::Instance::Temp)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setDescription("Device temperature");

  p.addPar ("OFF", 0, &Diode::Instance::off)
    .setUnit(U_LOGIC)
    .setCategory(CAT_CONTROL)
    .setDescription("Initial voltage drop across device set to zero");
  p.addPar ("LAMBERTW", 0, &Diode::Instance::lambertWFlag)
    .setUnit(U_LOGIC)
    .setCategory(CAT_CONTROL)
    .setDescription("Option to solve diode equations with the Lambert-W function");
}

void Traits::loadModelParameters(ParametricData<Diode::Model> &p)
{
  p.addPar ("IS", 1.0e-14, &Diode::Model::IS)
    .setUnit(U_AMP)
    .setCategory(CAT_CURRENT)
    .setDescription("Saturation current")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  // synonym for IS
  p.addPar ("JS", 1.0e-14, &Diode::Model::IS)
    .setUnit(U_AMP)
    .setCategory(CAT_CURRENT)
    .setDescription("Saturation current")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("RS",  0.0, &Diode::Model::RS)
    .setExpressionAccess(ParameterType::MIN_RES)
    .setUnit(U_OHM)
    .setCategory(CAT_RES)
    .setDescription("Parasitic resistance")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("N",  1.0, &Diode::Model::N)
    .setCategory(CAT_PROCESS)
    .setOriginalValueStored(true)
    .setDescription("Emission coefficient")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("ISR", 0.0, &Diode::Model::ISR)
    .setUnit(U_AMP)
    .setCategory(CAT_CURRENT)
    .setDescription("Recombination current parameter (level 2)")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("NR", 2.0, &Diode::Model::NR)
    .setCategory(CAT_NONE)
    .setDescription("Emission coefficient for ISR (level 2)")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("IKF", 0.0, &Diode::Model::IKF)
    .setUnit(U_AMP)
    .setCategory(CAT_CURRENT)
    .setDescription("High-injection \"knee\" current (level 2)")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("TT", 0.0, &Diode::Model::TT)
    .setUnit(U_SECOND)
    .setCategory(CAT_PROCESS)
    .setDescription("Transit time")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("CJO", 0.0, &Diode::Model::CJO)
    .setExpressionAccess(ParameterType::MIN_CAP)
    .setUnit(U_FARAD)
    .setCategory(CAT_CAP)
    .setDescription("Zero-bias p-n depletion capacitance")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  // synonyms for CJO
  p.addPar ("CJ", 0.0, &Diode::Model::CJO)
    .setExpressionAccess(ParameterType::MIN_CAP)
    .setUnit(U_FARAD)
    .setCategory(CAT_CAP)
    .setDescription("Zero-bias p-n depletion capacitance")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("CJ0", 0.0, &Diode::Model::CJO)
    .setExpressionAccess(ParameterType::MIN_CAP)
    .setUnit(U_FARAD)
    .setCategory(CAT_CAP)
    .setDescription("Zero-bias p-n depletion capacitance")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("VJ", 1.0, &Diode::Model::VJ)
    .setUnit(U_VOLT)
    .setCategory(CAT_VOLT)
    .setDescription("Potential for p-n junction")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("M", 0.5, &Diode::Model::M)
    .setCategory(CAT_PROCESS)
    .setDescription("Grading parameter for p-n junction")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("EG", 1.11, &Diode::Model::EG)
    .setUnit(U_EV)
    .setCategory(CAT_PROCESS)
    .setDescription("Bandgap voltage (barrier height)")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("XTI", 3.0, &Diode::Model::XTI)
    .setCategory(CAT_TEMP)
    .setDescription("IS temperature exponent")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("TIKF", 0.0, &Diode::Model::TIKF)
    .setUnit(U_DEGCM1)
    .setCategory(CAT_TEMP)
    .setDescription("IKF temperature coefficient (linear) (level 2)")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);
////
  p.addPar ("TBV1", 0.0, &Diode::Model::TBV1)
    .setUnit(U_DEGCM1)
    .setCategory(CAT_TEMP)
    .setDescription("BV temperature coefficient (linear) (level 2)");

  p.addPar ("TBV2", 0.0, &Diode::Model::TBV2)
    .setUnit(U_DEGCM2)
    .setCategory(CAT_TEMP)
    .setDescription("BV temperature coefficient (quadratic) (level 2)");

  p.addPar ("TRS1", 0.0, &Diode::Model::TRS1)
    .setUnit(U_DEGCM1)
    .setCategory(CAT_TEMP)
    .setDescription("RS temperature coefficient (linear) (level 2)");

  p.addPar ("TRS2", 0.0, &Diode::Model::TRS2)
    .setUnit(U_DEGCM2)
    .setCategory(CAT_TEMP)
    .setDescription("RS temperature coefficient (quadratic) (level 2)");
////
  p.addPar ("FC", 0.5, &Diode::Model::FC)
    .setCategory(CAT_CAP)
    .setDescription("Forward-bias depletion capacitance coefficient")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("BV", 1E99, &Diode::Model::BV)
    .setGivenMember(&Diode::Model::BVGiven)
    .setUnit(U_VOLT)
    .setCategory(CAT_VOLT)
    .setDescription("Reverse breakdown \"knee\" voltage")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  // synonym for BV
  p.addPar ("VB", 1E99, &Diode::Model::BV)
    .setGivenMember(&Diode::Model::BVGiven)
    .setUnit(U_VOLT)
    .setCategory(CAT_VOLT)
    .setDescription("Reverse breakdown \"knee\" voltage")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("IBV", 1.0e-3, &Diode::Model::IBV)
    .setUnit(U_AMP)
    .setCategory(CAT_CURRENT)
    .setDescription("Reverse breakdown \"knee\" current")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("IRF", 1.0, &Diode::Model::IRF)
    .setGivenMember(&Diode::Model::IRFGiven)
    .setCategory(CAT_CURRENT)
    .setDescription("Reverse current fitting factor")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("NBV", 1.0, &Diode::Model::NBV)
    .setCategory(CAT_PROCESS)
    .setDescription("Reverse breakdown ideality factor (level 2)")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("IBVL", 0.0, &Diode::Model::IBVL)
    .setUnit(U_AMP)
    .setCategory(CAT_CURRENT)
    .setDescription("Low-level reverse breakdown \"knee\" current (level 2)")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("NBVL", 1.0, &Diode::Model::NBVL)
    .setCategory(CAT_PROCESS)
    .setDescription("Low-level reverse breakdown ideality factor (level 2)")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("TNOM", 0.0, &Diode::Model::TNOM)
    .setCategory(CAT_NONE)
    .setDescription("")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("KF", 0.0, &Diode::Model::KF)
    .setCategory(CAT_FLICKER)
    .setDescription("Flicker noise coefficient")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);

  p.addPar ("AF", 1.0, &Diode::Model::AF)
    .setCategory(CAT_FLICKER)
    .setDescription("Flicker noise exponent")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&diodeSens);
}



std::vector< std::vector<int> > Instance::jacStamp_RS;
std::vector< std::vector<int> > Instance::jacStamp;

std::vector<int> Instance::jacMap_RS;
std::vector<int> Instance::jacMap;

std::vector< std::vector<int> > Instance::jacMap2_RS;
std::vector< std::vector<int> > Instance::jacMap2;

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Instance::processParams()
{
  // Set any non-constant parameter defaults:
  if (!given("LAMBERTW"))
    lambertWFlag = getDeviceOptions().lambertWFlag;
  if (!given("TEMP"))
    Temp = getDeviceOptions().temp.getImmutableValue<double>();
  if ( (model_.RS == 0.0) && (lambertWFlag == 1) )
    model_.RS =1.0e-12;

  updateTemperature( Temp );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & instance_block,
  Model &               model,
  const FactoryBlock &  factory_block)
  : DeviceInstance(instance_block, configuration.getInstanceParameters(), factory_block),
    model_(model),
    off(0),
    Area(1.0),
    InitCond(0.0),
    Temp(getDeviceOptions().temp.getImmutableValue<double>()),
    lambertWFlag(0),
    InitCondGiven(false),
    tJctPot(0.0),
    tJctCap(0.0),
    tDepCap(0.0),
    tSatCur(0.0),
    tVcrit(0.0),
    tF1(0.0),
    tBrkdwnV(0.0),
    tSatCurR(0.0),
    tIKF(0.0),
    tRS(0.0),
    tIRF(1.0),
    Id(0.0),
    Gd(0.0),
    Cd(0.0),
    Gcd(0.0),
    Qd(0.0),
    Icd(0.0),
    Gspr(0.0),
    Vpp(0.0),
    Vp(0.0),
    Vn(0.0),
    Vc(0.0),
    Vd(0.0),
    Vd_old(0.0),
    Vd_orig(0.0),
    li_storevd(-1),
    li_branch_data(-1),
    li_Pos(-1),
    li_Neg(-1),
    li_Pri(-1),
    APosEquPosNodeOffset(-1),
    APosEquPriNodeOffset(-1),
    ANegEquNegNodeOffset(-1),
    ANegEquPriNodeOffset(-1),
    APriEquPosNodeOffset(-1),
    APriEquNegNodeOffset(-1),
    APriEquPriNodeOffset(-1)
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    ,
    fPosEquPosNodePtr(0),
    fPosEquPriNodePtr(0),
    fNegEquNegNodePtr(0),
    fNegEquPriNodePtr(0),
    fPriEquPosNodePtr(0),
    fPriEquNegNodePtr(0),
    fPriEquPriNodePtr(0),
    qPosEquPosNodePtr(0),
    qPosEquPriNodePtr(0),
    qNegEquNegNodePtr(0),
    qNegEquPriNodePtr(0),
    qPriEquPosNodePtr(0),
    qPriEquNegNodePtr(0),
    qPriEquPriNodePtr(0)
#endif
{
  numIntVars   = 1;
  numExtVars   = 2;
  numStateVars = 0;
  setNumStoreVars(1);
  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1;    // this is the space to allocate if lead current or power is needed.

  if( jacStamp.empty() )
  {
    jacStamp_RS.resize(3);
    jacStamp_RS[0].resize(2);
    jacStamp_RS[0][0]=0;
    jacStamp_RS[0][1]=2;
    jacStamp_RS[1].resize(2);
    jacStamp_RS[1][0]=1;
    jacStamp_RS[1][1]=2;
    jacStamp_RS[2].resize(3);
    jacStamp_RS[2][0]=0;
    jacStamp_RS[2][1]=1;
    jacStamp_RS[2][2]=2;

    jacMap_RS.clear();
    jacStampMap(jacStamp_RS, jacMap_RS, jacMap2_RS,
                jacStamp,    jacMap, jacMap2, 2, 0, 3);

  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (instance_block.params);

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

  if ( (model_.RS == 0.0) || (lambertWFlag == 1) )
    numIntVars = 0;
  if ( model_.CJO == 0.0 )
    numStateVars = 1;

  if ( (model_.IRFGiven) )
  {
    UserWarning(*this) << " IRF has been specified in the model card.  This usage is deprecated.  Please see the Reference Guide for details";
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::~Instance()
{
}

// Additional Declarations
//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       : function for registering, and setting up, local ID's.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "In Instance::register LIDs\n\n";
    Xyce::dout() << "name             = " << getName() << std::endl;
    Xyce::dout() << "number of internal variables: " << numIntVars << std::endl;
    Xyce::dout() << "number of external variables: " << numExtVars << std::endl;
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the linear algebra
  // entities.  This assumes an order.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];

  //Setup of Pri node indices
  //If RS=0, Pri=Pos Node
  //--------------------------
  if( model_.RS && (lambertWFlag != 1) )
    li_Pri = intLIDVec[0];
  else
    li_Pri = li_Pos;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "\nSolution and RHS variables:\n";
    Xyce::dout() << "\nli_Pos = ";
    Xyce::dout().width(4);
    Xyce::dout() << li_Pos << std::endl;

    Xyce::dout() << "\nli_Neg = ";
    Xyce::dout().width(4);
    Xyce::dout() << li_Neg << std::endl;

    Xyce::dout() << "\nli_Pri = ";
    Xyce::dout().width(4);
    Xyce::dout() << li_Pri << std::endl;

  }



  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
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
  if (li_Pos != li_Pri)
    addInternalNode(symbol_table, li_Pri, getName(), "internal");

  addStoreNode(symbol_table, li_storevd, getName(), "vd");

  if (loadLeadCurrent)
  {
    addBranchDataNode( symbol_table, li_branch_data, getName(), "BRANCH_D");
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // copy over the global ID lists:
  staLIDVec = staLIDVecRef;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs(
  const std::vector<int> & stoLIDVecRef )
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());

// copy over the global ID lists:
  stoLIDVec = stoLIDVecRef;
  li_storevd = stoLIDVec[0];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "li_storevd = " << li_storevd;
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerBranchDataLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 12/18/2012
//-----------------------------------------------------------------------------
/// In addition to state vector, Xyce maintains a separate datastructure
/// called a "branch data" vector.  As with other such vectors, the device
/// declares at construction time how many branch vector entries it needs,
/// and later Topology assigns locations for devices, returning LIDs.
///
/// These LIDs are stored in this method for later use.
///

void Instance::registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef)
{
  AssertLIDs(branchLIDVecRef.size() == getNumBranchDataVars());

  if (loadLeadCurrent)
  {
    li_branch_data= branchLIDVecRef[0];
  }
}


//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/2/02
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  if( model_.RS && (lambertWFlag != 1) )
    return jacStamp_RS;
  else
    return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/2/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  std::vector<int> map;
  std::vector< std::vector<int> > map2;

  if( model_.RS && (lambertWFlag != 1) )
  {
    map = jacMap_RS;
    map2 = jacMap2_RS;
  }
  else
  {
    map = jacMap;
    map2 = jacMap2;
  }

  APosEquPosNodeOffset = jacLIDVec[map[0]][map2[0][0]];
  APosEquPriNodeOffset = jacLIDVec[map[0]][map2[0][1]];

  ANegEquNegNodeOffset = jacLIDVec[map[1]][map2[1][0]];
  ANegEquPriNodeOffset = jacLIDVec[map[1]][map2[1][1]];

  APriEquPosNodeOffset = jacLIDVec[map[2]][map2[2][0]];
  APriEquNegNodeOffset = jacLIDVec[map[2]][map2[2][1]];
  APriEquPriNodeOffset = jacLIDVec[map[2]][map2[2][2]];

}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/30/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);

  fPosEquPosNodePtr = &(dFdx[li_Pos][APosEquPosNodeOffset]);
  fPosEquPriNodePtr = &(dFdx[li_Pos][APosEquPriNodeOffset]);
  fNegEquNegNodePtr = &(dFdx[li_Neg][ANegEquNegNodeOffset]);
  fNegEquPriNodePtr = &(dFdx[li_Neg][ANegEquPriNodeOffset]);
  fPriEquPosNodePtr = &(dFdx[li_Pri][APriEquPosNodeOffset]);
  fPriEquNegNodePtr = &(dFdx[li_Pri][APriEquNegNodeOffset]);
  fPriEquPriNodePtr = &(dFdx[li_Pri][APriEquPriNodeOffset]);

  qPosEquPosNodePtr = &(dQdx[li_Pos][APosEquPosNodeOffset]);
  qPosEquPriNodePtr = &(dQdx[li_Pos][APosEquPriNodeOffset]);
  qNegEquNegNodePtr = &(dQdx[li_Neg][ANegEquNegNodeOffset]);
  qNegEquPriNodePtr = &(dQdx[li_Neg][ANegEquPriNodeOffset]);
  qPriEquPosNodePtr = &(dQdx[li_Pri][APriEquPosNodeOffset]);
  qPriEquNegNodePtr = &(dQdx[li_Pri][APriEquNegNodeOffset]);
  qPriEquPriNodePtr = &(dQdx[li_Pri][APriEquPriNodeOffset]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 diode instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/26/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  // note: the only capacitor goes from the negative to the
  // positive-prime node, so there is not contribution in this
  // function to the positive node.

  // load in the KCL for the negative node:
  double coef = Qd;
  (extData.daeQVectorRawPtr)[li_Neg] -= coef;

  // load in the KCL for the positive prime node:
  coef *= -1.0; // -Qd
  (extData.daeQVectorRawPtr)[li_Pri] -= coef;

  // load the voltage limiter vector.
  if( getDeviceOptions().voltageLimiterFlag )
  {
    double Vd_diff = Vd - Vd_orig;
    double Cd_Jdxp = 0.0;
    Cd_Jdxp = -( Cd ) * Vd_diff;

    // Load the dQdxdVp vector
    (extData.dQdxdVpVectorRawPtr)[li_Neg] += Cd_Jdxp;
    (extData.dQdxdVpVectorRawPtr)[li_Pri] -= Cd_Jdxp;
  }

  if( loadLeadCurrent && (model_.CJO != 0.0))
  {
    (extData.nextLeadCurrQCompRawPtr)[li_branch_data] = Qd;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 diode instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/26/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  // 3f5 compatible currents
  // Including derivation of Vd_diff and Limiting Correction
  //---------------------------------------------------------
  double Ir = Gspr * (Vp - Vpp);

  // load in the KCL for the positive node:
  double coef = -Ir;
  (extData.daeFVectorRawPtr)[li_Pos] -= coef;

  // load in the KCL for the negative node:
  coef = Id;
  (extData.daeFVectorRawPtr)[li_Neg] -= coef;

  // load in the KCL for the positive prime node:
  coef *= -1;
  coef += Ir;
  (extData.daeFVectorRawPtr)[li_Pri] -= coef;

  // load the voltage limiter vector.
  if( getDeviceOptions().voltageLimiterFlag )
  {
    double Vd_diff = Vd - Vd_orig;
    double Gd_Jdxp = 0.0;
    Gd_Jdxp = -( Gd ) * Vd_diff;

    // Load the dFdxdVp vector
    (extData.dFdxdVpVectorRawPtr)[li_Neg] += Gd_Jdxp;
    (extData.dFdxdVpVectorRawPtr)[li_Pri] -= Gd_Jdxp;
  }

  if( loadLeadCurrent )
  {
    (extData.nextLeadCurrFCompRawPtr)[li_branch_data] = Id;
    (extData.nextJunctionVCompRawPtr)[li_branch_data] = (extData.nextSolVectorRawPtr)[li_Pos] - (extData.nextSolVectorRawPtr)[li_Neg];
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 diode instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);

// New Spice3f5 matched conductivities
// Only major difference in loads is support for Pos=Pri node when RS=0
// For this DAE dQdx load, the capacitance contribution (Gcd) is the
// only part used.
//---------------------------------------------------------------------

  dQdx[li_Neg][ANegEquNegNodeOffset] += Cd;
  dQdx[li_Neg][ANegEquPriNodeOffset] -= Cd;

  dQdx[li_Pri][APriEquNegNodeOffset] -= Cd;
  dQdx[li_Pri][APriEquPriNodeOffset] += Cd;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 diode instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

// New Spice3f5 matched conductivities
// Only major difference in loads is support for Pos=Pri node when RS=0
// For this DAE dFdx load, the capacitance contribution (Gcd) is removed.
//---------------------------------------------------------------------

  dFdx[li_Pos][APosEquPosNodeOffset] += Gspr;
  dFdx[li_Pos][APosEquPriNodeOffset] -= Gspr;

  dFdx[li_Neg][ANegEquNegNodeOffset] += Gd;
  dFdx[li_Neg][ANegEquPriNodeOffset] -= Gd;

  dFdx[li_Pri][APriEquPosNodeOffset] -= Gspr;
  dFdx[li_Pri][APriEquNegNodeOffset] -= Gd;
  dFdx[li_Pri][APriEquPriNodeOffset] += Gspr + Gd;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       : update primary state for one diode instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 1/10/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  double * stoVec = extData.nextStoVectorRawPtr;
  bool bsuccess = updateIntermediateVars ();
  stoVec[li_storevd] = Vd;
  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Diode::Instance::getNumNoiseSources
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/27/2014
//-----------------------------------------------------------------------------
int Instance::getNumNoiseSources () const
{
  return 3;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Diode::Instance::setupNoiseSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void Instance::setupNoiseSources (Xyce::Analysis::NoiseData & noiseData)
{
  int numSources=3;
  noiseData.numSources = numSources;
  noiseData.resize(numSources);

  noiseData.deviceName = getName().getEncodedName();

  // Note: the letter suffixes (e.g., rs) are used by the DNO() and DNI() 
  // operators for .PRINT NOISE
  noiseData.noiseNames[0] = "noise_" + getName().getEncodedName()+
    std::string("_rs");  // noise due to rs 
  noiseData.noiseNames[1] = "noise_" + getName().getEncodedName()+
    std::string("_id");  // noise due to id 
  noiseData.noiseNames[2] = "noise_" + getName().getEncodedName()+
    std::string("_fn"); // flicker (1/f) noise

  noiseData.li_Pos[0] = li_Pri;
  noiseData.li_Neg[0] = li_Pos;

  noiseData.li_Pos[1] = li_Pri;
  noiseData.li_Neg[1] = li_Neg;

  noiseData.li_Pos[2] = li_Pri;
  noiseData.li_Neg[2] = li_Neg;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Diode::Instance::getNoiseSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void Instance::getNoiseSources (Xyce::Analysis::NoiseData & noiseData)
{
  // thermal noise from RS
  devSupport.noiseSupport(
      noiseData.noiseDens[0], noiseData.lnNoiseDens[0], 
      THERMNOISE, Gspr, Temp);

  // shot noise from the diode depletion region
  devSupport.noiseSupport( 
      noiseData.noiseDens[1], noiseData.lnNoiseDens[1], 
      SHOTNOISE, Id, Temp);

  // flicker noise
  noiseData.noiseDens[2] = model_.KF * std::exp(model_.AF * std::log(std::max(fabs(Id),N_MINLOG))) / noiseData.freq;
  noiseData.lnNoiseDens[2] = std::log(std::max(noiseData.noiseDens[2],N_MINLOG));
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one diode instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 1/10/01
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;

  double * solVec = extData.nextSolVectorRawPtr;

  //diode parameters
  double M;       // grading parameter
  double BV;      // breakdown voltage
  double IBV;     // reverse breakdown current
  double NBV;     // reverse breakdown ideality factor
  double IBVL;    // low-level reverse breakdown current
  double NBVL;    // low-level reverse breakdown ideality factor
  double N;       // non-ideality factor.
  double NR;      // emission coeff. for ISR.
  double TT;      // transit time.
  double F2;      // capacitive polynomial factor
  double F3;      // capacitive polynomial factor
  double Inorm;   // normal diffusion current
  double Irec;    // recombination current
  double Kgen;    // generation factor
  double Khi;     // high-injection factor
  double Gd1, DKgen;
  double Gd2, DKhi;

  //Model Diode parameters
  M    = model_.M;
  BV   = model_.BV;
  IBV  = model_.IBV;
  NBV  = model_.NBV;
  IBVL = model_.IBVL;
  NBVL = model_.NBVL;
  N    = model_.N;
  NR   = model_.NR;
  TT   = model_.TT;
  F2   = model_.F2;
  F3   = model_.F3;

  // obtain voltage drop accross the capacitor:
  Vp  = Vn = Vpp = 0.0;
  Vpp = solVec[li_Pri];
  Vn  = solVec[li_Neg];
  Vp  = solVec[li_Pos];

  // Junction Voltage
  Vd = Vpp - Vn;

  double Isat  = tSatCur * Area;
  double IsatR = tSatCurR * Area;
  double Vt    = CONSTKoverQ * Temp;
  double Vte   = N * Vt;
  double VteR  = NR * Vt;

  Gspr = tCOND * Area;

  Vd_orig = Vd;
  origFlag = true;

  // Setup initial junction conditions if UIC enabled
  //------------------------------------------------
  if (getSolverState().newtonIter == 0)
  {
    //Vd_old = Vd;
    if (getSolverState().initJctFlag_ && getDeviceOptions().voltageLimiterFlag)
    {
      if (InitCondGiven)
      {
        Vd = InitCond;
        origFlag = false;
      }
      else if (off)
      {
        Vd = 0.0;
        origFlag = false;
      }
      else
      {
        if (getSolverState().inputOPFlag)
        {
          Linear::Vector * flagSolVectorPtr = extData.flagSolVectorPtr;
          if ((*flagSolVectorPtr)[li_Pos] == 0 ||
              (*flagSolVectorPtr)[li_Neg] == 0 ||
              (*flagSolVectorPtr)[li_Pri] == 0)
          {
            Vd=tVcrit;
            origFlag = false;
          }
        }
        else
        {
          Vd=tVcrit;
          origFlag = false;
        }
      }
    }

    // assume there is no history -- then check if the
    // state vector can overwrite this
    Vd_old = Vd;

    if (!(getSolverState().dcopFlag)||(getSolverState().locaEnabledFlag && getSolverState().dcopFlag))
    {
      Vd_old = (extData.currStoVectorRawPtr)[li_storevd];
    }
  }
  else  // just do this whenever it isn't the first iteration:
  {
    Vd_old = (extData.nextStoVectorRawPtr)[li_storevd];
  }

  // Voltage limiting based on mode of diode
  //---------------------------------------
  if (getDeviceOptions().voltageLimiterFlag)
  {
    int ichk = 0;

    if (getSolverState().newtonIter >= 0)
    {
      //Test if breakdown voltage given or not
      if (model_.BVGiven && (Vd < Xycemin(0.0, -BV + 10.0 * Vte)))
      {
        double Vdtmp = -( BV + Vd );
        Vdtmp = devSupport.pnjlim(Vdtmp, -(Vd_old+BV), Vte, tVcrit, &ichk);
        Vd    = -(Vdtmp + BV);
      }
      else
        Vd = devSupport.pnjlim(Vd, Vd_old, Vte, tVcrit, &ichk);

      if (ichk) origFlag = false;
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() << "Instance::updateIntermediateVars " << getName()<<std::endl;
  }

  // Current and Conductivity
  //----------------------------------------------

  int level = model_.getLevel();
  if(level == 1)
  {
    // Using LambertW function
    if (lambertWFlag == 1)
    {
      double RS = model_.RS;
      if (Vd >= -3.0 * Vte)
      {
        lambertWCurrent(Isat, Vte, RS);
      }
      // linear reverse bias
      else if ( !tBrkdwnV || (Vd >= -tBrkdwnV) )
      {
        lambertWLinearReverseBias(Isat, Vte, RS);
      }
      // reverse breakdown
      else
      {
        lambertWBreakdownCurrent(Isat, Vte, RS);
      }

      double Vrs;
      Vrs = (Id + Icd)*RS;
      Vc = Vd - Vrs;
    }
    else if (lambertWFlag == 2)
    {
      if (Vd >= -3.0 * Vte)
      {
        lambertWCurrent(Isat, Vte, 1.0e-15);
      }
      // linear reverse bias
      else if ( !tBrkdwnV || (Vd >= -tBrkdwnV) )
      {
        lambertWLinearReverseBias(Isat, Vte, 1.0e-15);
      }
      // reverse breakdown
      else
      {
        lambertWBreakdownCurrent(Isat, Vte, 1.0e-15);
      }
      Vc = Vd;
    }

    // Normal exponential
    else
    {
      // Adjustment for linear portion of reverse current
      double IRfactor;
      //if(Vd >= 0) IRfactor = 1.0;  //error?  shouldn't this be Vd>=-3*Vte????
      if(Vd >= -3.0 * Vte) IRfactor = 1.0;  //changing this to be consistent
      else IRfactor = tIRF;                 //with model in reference guide.
      Isat *= IRfactor;

      if (Vd >= -3.0 * Vte)
      {
        double arg1 = Vd / Vte;
        arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
        double evd = exp(arg1);

        Id = Isat * (evd - 1.0) + getDeviceOptions().gmin * Vd;
        Gd = Isat * evd / Vte + getDeviceOptions().gmin;
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
        {
          Xyce::dout()  << "Normal exponential regime." << std::endl;
          Xyce::dout()  << " Vd  = " << Vd << std::endl;
          Xyce::dout()  << " Vte = " << Vte << std::endl;
          Xyce::dout()  << " Id  = " << Id  << std::endl;
          Xyce::dout()  << " Gd  = " << Gd  << std::endl;
        }

      }

      // Linear reverse bias
      else if(!tBrkdwnV || (Vd >= -tBrkdwnV))
      {
        double arg = 3.0 * Vte / (Vd * CONSTe);
        arg = arg * arg * arg;
        Id = -Isat * (1.0 + arg) + getDeviceOptions().gmin * Vd;
        Gd = Isat * 3.0 * arg / Vd + getDeviceOptions().gmin;
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
        {
          Xyce::dout()  << "Linear reverse bias regime." << std::endl;
          Xyce::dout()  << " Vd       = " << Vd << std::endl;
          Xyce::dout()  << " tBrkdwnV = " << tBrkdwnV << std::endl;
          Xyce::dout()  << " Id       = " << Id  << std::endl;
          Xyce::dout()  << " Gd       = " << Gd  << std::endl;
        }
      }

      // Reverse breakdown
      else
      {
        double arg1 = -(tBrkdwnV + Vd) / Vte;
        arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
        double evrev = exp(arg1);

#ifdef Xyce_BREAKDOWN_ORIGINAL
        Id = -Isat * evrev + getDeviceOptions().gmin * Vd;
        Gd = Isat * evrev / Vte + getDeviceOptions().gmin;
#else
        //added by K. Santarelli 9/18/07 to account for change in tBrkdwnV
        //calculation.
        double arg2=3.0*Vte/(CONSTe*tBrkdwnV);
        arg2=arg2*arg2*arg2;
        double Isat_tBrkdwnV=Isat*(1-arg2);
        Id = -Isat_tBrkdwnV * evrev + getDeviceOptions().gmin * Vd;
        Gd = Isat_tBrkdwnV * evrev / Vte + getDeviceOptions().gmin;
#endif

        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
        {
          Xyce::dout()  << "Reverse breakdown regime." << std::endl;
          Xyce::dout()  << " Vd       = " << Vd << std::endl;
          Xyce::dout()  << " tBrkdwnV = " << tBrkdwnV << std::endl;
          Xyce::dout()  << " Id       = " << Id  << std::endl;
          Xyce::dout()  << " Gd       = " << Gd  << std::endl;
        }
      }
      Vc = Vd;
    }

  }
  else if(level == 2)
  {
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << " Level 2 diode code " << std::endl;
    }
    // Adjustment for linear portion of reverse current
    double IRfactor;
    //if(Vd >= 0) IRfactor = 1.0;    //error?  Shouldn't this be Vd >= -3*Vte?
    if(Vd >= -3.0*Vte) IRfactor=1.0; //changing it to be consistent with model
    else IRfactor = tIRF;            //in the reference manual.
    Isat *= IRfactor;

    if (Vd >= -3.0 * Vte)
    {
      double arg1 = Vd / Vte;
      arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
      double evd = exp(arg1);
      Inorm = Isat * (evd - 1.0) + getDeviceOptions().gmin * Vd;
      Gd1 = Isat*evd/Vte + getDeviceOptions().gmin;

      arg1 = Vd / VteR;
      arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
      evd  = exp(arg1);
      Irec = IsatR * (evd - 1.0);
      Gd2  = IsatR*evd/VteR;

      Khi = 1;
      DKhi = 0;
      if(tIKF > 0)
      {
        Khi = sqrt(tIKF/(tIKF+Inorm));
        DKhi = 0.5*Khi*Gd1/(tIKF+Inorm);
      }
      Kgen = 0;
      DKgen = 0;
      if(Irec != 0)
      {
        Kgen = sqrt( pow(((1-Vd/tJctPot)*(1-Vd/tJctPot) + 0.005),M) );
        DKgen = -M*(1-Vd/tJctPot)/(tJctPot*Kgen);
      }

      Id = Inorm*Khi + Irec*Kgen;
      Gd = Gd1*Khi + Inorm*DKhi + Gd2*Kgen + Irec*DKgen;
      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
      {
        Xyce::dout()  << "L2 Normal exponential regime." << std::endl;
        Xyce::dout()  << " Vd  = " << Vd << std::endl;
        Xyce::dout()  << " Vte = " << Vte << std::endl;
        Xyce::dout()  << " Id  = " << Id  << std::endl;
        Xyce::dout()  << " Irec= " << Irec  << std::endl;
        Xyce::dout()  << " Gd  = " << Gd  << std::endl;
        Xyce::dout()  << " Gd1 = " << Gd1 << std::endl;
        Xyce::dout()  << " Gd2 = " << Gd2 << std::endl;
        Xyce::dout()  << " Khi = " << Khi << std::endl;
        Xyce::dout()  << " Kgen=" << Kgen << std::endl;
        Xyce::dout()  << "DKgen=" <<DKgen << std::endl;
      }
    }

    // Linear reverse bias
    else if(!tBrkdwnV || (Vd >= -tBrkdwnV))
    {
      double arg = 3.0 * Vte / (Vd * CONSTe);
      arg = arg * arg * arg;
      Id = -Isat * (1.0 + arg) + getDeviceOptions().gmin * Vd;
      Gd = Isat * 3.0 * arg / Vd + getDeviceOptions().gmin;
      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
      {
        Xyce::dout()  << "L2 Linear reverse bias regime." << std::endl;
        Xyce::dout()  << " Vd       = " << Vd << std::endl;
        Xyce::dout()  << " tBrkdwnV = " << tBrkdwnV << std::endl;
        Xyce::dout()  << " Id       = " << Id  << std::endl;
        Xyce::dout()  << " Gd       = " << Gd  << std::endl;
      }
    }

    // Reverse breakdown
    else
    {
      double arg1 = -(tBrkdwnV + Vd) / Vte;
      arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
      double evrev = exp(arg1);

#ifdef Xyce_BREAKDOWN_ORIGINAL
      Id = -Isat * evrev + getDeviceOptions().gmin * Vd;
      Gd = Isat * evrev / Vte + getDeviceOptions().gmin;
#else
      //added 9/18/07 by K. Santarelli to account for change in tBrkdwnV
      //calculation.
      double arg2=3.0*Vte/(CONSTe*tBrkdwnV);
      arg2=arg2*arg2*arg2;
      double Isat_tBrkdwnV=Isat*(1-arg2);
      Id = -Isat_tBrkdwnV * evrev + getDeviceOptions().gmin * Vd;
      Gd = Isat_tBrkdwnV * evrev / Vte + getDeviceOptions().gmin;
#endif
      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
      {
        Xyce::dout()  << "L2 Reverse breakdown regime." << std::endl;
        Xyce::dout()  << " Vd       = " << Vd << std::endl;
        Xyce::dout()  << " tBrkdwnV = " << tBrkdwnV << std::endl;
        Xyce::dout()  << " Id       = " << Id  << std::endl;
        Xyce::dout()  << " Gd       = " << Gd  << std::endl;
      }
    }
    Vc = Vd;

  }  // level

  // Only compute if Capacitance is non-zero
  //---------------------------------------
  if (tJctCap != 0.0)
  {
    // Charge Storage
    double Czero = tJctCap * Area;
    if (Vc < tDepCap)
    {
      //double arg = 1.0 - Vd/VJ;
      double arg = 1.0 - Vc / tJctPot;
      double arg1 = -M * log(arg);
      arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
      double sarg = exp(arg1);

      //Qd = TT*Id + VJ*Czero*(1.0-arg*sarg)/(1.0-M);
      Qd = TT * Id + tJctPot * Czero * (1.0 - arg * sarg) / (1.0 - M);
      Cd = TT * Gd + Czero * sarg;
    }
    else
    {
      double Czof2 = Czero / F2;
      double MotJctPot = M / tJctPot;

      //Qd = TT*Id + Czero*tF1 + Czof2*(F3*(Vd-tDepCap)+(M/(VJ+VJ))
      Qd = TT * Id + Czero * tF1 +
           Czof2 * (F3 * (Vc - tDepCap) + (0.5 * MotJctPot) *
                    (Vc * Vc - tDepCap * tDepCap));
      //Cd = TT*Gd + Czof2*(F3+M*Vd/VJ);
      Cd = TT * Gd + Czof2 * (F3 + MotJctPot * Vc);
    }
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "Qd = " << Qd << std::endl;
      Xyce::dout() << "Cd = " << Qd << std::endl;
    }
  }
  else
  {
    Qd = 0.0;
    Cd = 0.0;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       : update intermediate variables for one diode instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 1/10/01
//-----------------------------------------------------------------------------
bool Instance::updateTemperature( const double & temp )
{

  double vtnom = CONSTKoverQ * model_.TNOM;

  double xfc = log( 1.0 - model_.FC );

  if( temp != -999.0 ) Temp = temp;
  double TNOM = model_.TNOM;

  double vt = CONSTKoverQ * Temp;
  double fact2 = Temp / CONSTREFTEMP;
  double egfet = CONSTEg0 - (CONSTalphaEg*Temp*Temp)/(Temp+CONSTbetaEg);
  double arg = -egfet/(2.0*CONSTboltz*Temp) +
               CONSTEg300/(CONSTboltz*(CONSTREFTEMP+CONSTREFTEMP));
  double pbfact = -2.0*vt*(1.5*log(fact2)+CONSTQ*arg);
  double egfet1 = CONSTEg0 - (CONSTalphaEg*model_.TNOM*model_.TNOM)/
                  (model_.TNOM+CONSTbetaEg);
  double arg1 = -egfet1/(2.0*CONSTboltz*model_.TNOM) +
		CONSTEg300/(2.0*CONSTboltz*CONSTREFTEMP);
  double fact1 = model_.TNOM/CONSTREFTEMP;
  double pbfact1 = -2.0*vtnom*(1.5*log(fact1)+CONSTQ*arg1);

  double pbo = (model_.VJ-pbfact1)/fact1;
  double gmaold = (model_.VJ-pbo)/pbo;

  tJctCap = model_.CJO/
            (1.0+model_.M*(4.0e-4*(model_.TNOM-CONSTREFTEMP) -gmaold));

  tJctPot = pbfact+fact2*pbo;

  double gmanew = (tJctPot-pbo)/pbo;

  tJctCap *= 1.0 + model_.M*(4.0e-4*(Temp-CONSTREFTEMP)-gmanew);

  tSatCur = model_.IS*exp(((Temp/model_.TNOM)-1.0)*
                          model_.EG/(model_.N*vt)+
                          model_.XTI/model_.N*log(Temp/model_.TNOM));

  tF1 = tJctPot*(1.0-exp((1.0-model_.M)*xfc))/(1.0-model_.M);

  tDepCap = model_.FC*tJctPot;

  double vte = model_.N*vt;
  double tempBV;

  tVcrit = vte*log(vte/(CONSTroot2*tSatCur));
  tRS   = model_.RS;
  tCOND = model_.COND;
  if (model_.IRFGiven)
  {
    tIRF  = model_.IRF*pow(fact2,1.6);
  }
  else
  {
    tIRF = model_.IRF;
  }

  int level = model_.getLevel();
  if(level == 2)   // this section is PSPICE compatible
  {
    tSatCurR = model_.ISR*exp((Temp/TNOM - 1.0)*
                              model_.EG/(model_.NR*vt)+
                              model_.XTI/model_.NR*log(Temp/TNOM));

    tIKF = model_.IKF*(1 + model_.TIKF*(Temp-TNOM));

    tempBV = model_.BV*(1 + (Temp-TNOM)*
                        ( model_.TBV1 + model_.TBV2*(Temp-TNOM) ));

    tRS = model_.RS*(1 + (Temp-TNOM)*
                     ( model_.TRS1 + model_.TRS2*(Temp-TNOM) ));

    tCOND = 0.0;
    if(tRS != 0.0) tCOND = 1.0/tRS;

    tJctPot = (model_.VJ - egfet1)*fact2 - 3*vt*log(fact2) + egfet;

    tJctCap = model_.CJO/(1.0 +
                          model_.M*(4.0e-4*(Temp-TNOM) + (1-tJctPot/model_.VJ)));
  }
  else
  {
    tempBV=model_.BV;
  }

#ifdef Xyce_BREAKDOWN_ORIGINAL
  //Changed 9/18/07, K. R. Santarelli.  This is the original (broken) version
  //of the breakdown voltage computation.  It has two known issues:
  //
  //1.  It assumes N (the emission coefficient) is always 1 (division by
  //    vt instead of vte).
  //2.  While the for loop terminates due to the fact that it's implemented via
  //    a counter, the value of xbv does not converge in many cases, so the
  //    value of xbv upon terminating is often garbage.
  //
  //For consistency with old, existing simulations, this version of the
  //breakdown voltage calculation (along with the appropriate code for
  //computing the breakdown current in the level 1 and level 2 models---see
  //the appropriate sections of
  //Instance::updateIntermediateVars) can be invoked by
  //specifiying the CPPFLAG "-DXyce_BREAKDOWN_ORIGINAL" when creating a Xyce
  //Build.  The default, however, is the code which follows #else.

  double reltol = 1.0e-3;
  if( model_.BVGiven )
  {
    double IRfactor = tIRF;
    double cbv = model_.IBV;
    double xbv, xcbv;
    if( cbv < IRfactor*tSatCur*tempBV/vt )
    {
      cbv = IRfactor*tSatCur*tempBV/vt;
      xbv = tempBV;
    }
    else
    {
      double tol = reltol*cbv;
      xbv = tempBV-vt*log(1.0+cbv/(IRfactor*tSatCur));
      for( int i = 0; i < 25; ++i )
      {
        xbv = tempBV-vt*log(cbv/(IRfactor*tSatCur)+1.0-xbv/vt);
        xcbv = IRfactor*tSatCur*(exp((tempBV-xbv)/vt)-1.0+xbv/vt);
        if(fabs(xcbv-cbv)<=tol) break;
      }
    }
    tBrkdwnV = xbv;
  }
#else
  if( model_.BVGiven)
  {
    double IRFactor=tIRF;
    double cbv = model_.IBV;
    double xbv;
    double cthreshlow; //lower threshold for IBV
    double cthreshhigh; //high threshold for IBV
    int iter_count;
    const int ITER_COUNT_MAX=8;
    double arg2;

    double arg1=3.0*vte/(CONSTe*tempBV);
    arg1=arg1*arg1*arg1;
    cthreshlow=tSatCur*IRFactor*(1-arg1);
    cthreshhigh=tSatCur*IRFactor*(1-1.0/(CONSTe*CONSTe*CONSTe)) *
                exp(-1.0*(3.0*vte-tempBV)/vte);

    if(cbv >= cthreshhigh)
    {                     //if IBV is too high, tBrkdwnV will go below 3NVt.
      tBrkdwnV=3.0*vte;   //Clip tBrkdwnV to 3*N*Vt in this case (and hence
    }                     //clip IBV to cthreshhigh).


    else if(cbv <= cthreshlow)
    {                          //if IBV is too low, tBrkdwnV will go above
      tBrkdwnV=tempBV;  //BV.  Clip tBrkdwnV to BV in this case (and
    }                          //hence clip IBV to cthreshlow).


    //If IBV is in an acceptable range, perform a Picard iteration to find
    //tBrkdwnV, starting with an initial guess of tBrkdwnV=tempBV, and
    //running through the iteration ITER_COUNT_MAX times.

    else
    {
      xbv=tempBV;
      for(iter_count=0; iter_count < ITER_COUNT_MAX; iter_count++)
      {
        arg2=3.0*vte/(CONSTe*xbv);
        arg2=arg2*arg2*arg2;
        xbv=tempBV-vte*log(cbv/(tSatCur*IRFactor))+vte *
            log(1-arg2);
      }
      tBrkdwnV=xbv;
    }

//Note that, not only is tBrkdwnV adjusted using this method, but the effective
//value of Is (tSatCur) is adjusted as well.  The code used to use Is before,
//but now it will use Is*IRFactor*(1-(3*N*Vt/(e*tBrkdwnV))^3) instead.  This is
//reflected in changes made to Instance::updateIntermediateVars
//(for "normal" level 1 and 2 diode models) for the reverse breakdown region.
//Changes have not been made to the associated LambertW functions as of yet
//since there seem to be other issues involved with those functions independent
// of this change.
  }
#endif


  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() << "Instance::UpdateTemperature" << getName() <<std::endl;
    Xyce::dout() << " IS = " << model_.IS << std::endl;
    Xyce::dout() << " vtnom   = " << vtnom << std::endl;
    Xyce::dout() << " xfc     = " << xfc << std::endl;
    Xyce::dout() << " TNOM    = " << TNOM << std::endl;
    Xyce::dout() << " vt      = " << vt << std::endl;
    Xyce::dout() << " fact2   = " << fact2 << std::endl;
    Xyce::dout() << " egfet   = " << egfet << std::endl;
    Xyce::dout() << " arg     = " << arg   << std::endl;
    Xyce::dout() << " pbfact  = " << pbfact   << std::endl;
    Xyce::dout() << " egfet1  = " << egfet1   << std::endl;
    Xyce::dout() << " arg1    = " << arg1     << std::endl;
    Xyce::dout() << " fact1   = " << fact1    << std::endl;
    Xyce::dout() << " pbfact1 = " << pbfact1  << std::endl;
    Xyce::dout() << " pbo     = " << pbo      << std::endl;
    Xyce::dout() << " gmaold  = " << gmaold   << std::endl;
    Xyce::dout() << " gmanew  = " << gmanew   << std::endl;
    Xyce::dout() << " tJctCap = " << tJctCap  << std::endl;
    Xyce::dout() << " tJctPot = " << tJctPot  << std::endl;
    Xyce::dout() << " tSatCur = " << tSatCur  << std::endl;
    Xyce::dout() << " tF1     = " << tF1      << std::endl;
    Xyce::dout() << " tDepCap = " << tDepCap  << std::endl;
    Xyce::dout() << " vte     = " << vte      << std::endl;
    Xyce::dout() << " tempBV  = " << tempBV   << std::endl;
    Xyce::dout() << " tVcrit  = " << tVcrit   << std::endl;
    Xyce::dout() << " tRS     = " << tRS      << std::endl;
    Xyce::dout() << " tCOND   = " << tCOND    << std::endl;
    Xyce::dout() << " tIRF    = " << tIRF     << std::endl;
    Xyce::dout() << " tBrkdwnV= " << tBrkdwnV << std::endl;
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Instance::lambertWCurrent
// Purpose       : Determine the diode current using Lambert-W function
// Special Notes :
// Scope         : public
// Creator       : Nick Johnson, Summer Intern
// Creation Date : 7/5/02
//-----------------------------------------------------------------------------
bool Instance::lambertWCurrent(double Isat, double Vte, double RS)
{
  // with capacitor current and using LambertW accros whole model
  /*  if (Cd != 0.0 && Icd != 0 && (lambertWFlag == 1))
      {
      double AA = Vte*(1.0 + RS*getSolverState().pdt*Cd);
      double XX = Icd - getSolverState().pdt*Qd;
      double arg1 = (Vd + RS*Isat - RS*XX)/AA;
      arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
      double evd = exp(arg1);
      double lambWArg = Isat*RS/AA * evd;
      double lambWReturn;
      int ierr;
      double lambWError;
      devSupport.lambertw(lambWArg, lambWReturn, ierr, lambWError);

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
      if (ierr == 0)
      Xyce::dout() << "Safe LambertW return" << std::endl;
      else if (ierr == 1)
      Xyce::dout() << "LambertW argument not in domain" << std::endl;
      else
      Xyce::dout() << "Arithmetic problems with LambertW" << std::endl;
      }

      double prefac = AA/RS;
      Id = -Isat + prefac*lambWReturn;
      Gd = lambWReturn / ((1 + lambWReturn)*RS);

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
      Xyce::dout() << "lambWArg   = " << lambWArg << std::endl;
      Xyce::dout() << "lambWError = " << lambWError << std::endl;
      Xyce::dout() << "Id         = " << Id << std::endl;
      Xyce::dout() << "Gd         = " << Gd << std::endl;
      Xyce::dout() << "Icd        = " << Icd << std::endl;
      Xyce::dout() << "Cd         = " << Cd  << std::endl;
      Xyce::dout() << "Qd         = " << Qd  << std::endl;
      Xyce::dout() << "AA         = " << AA << std::endl;
      Xyce::dout() << "XX         = " << XX << std::endl;
      Xyce::dout() << "Using lambertwCurrent w/ capacitance" << std::endl;
      }
      #endif
      }

      // when capacitor current=0, there is no capacitance, or using LambertW
      // only across the diode element
      else
      { */
  double arg1 = (Vd + Isat*RS)/Vte;
  arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
  double evd = exp(arg1);
  double lambWArg = Isat*RS*evd/Vte;
  double lambWReturn;
  int ierr;
  double lambWError;
  devSupport.lambertw(lambWArg, lambWReturn, ierr, lambWError);

  Id = -Isat+Vte*(lambWReturn)/RS;
  Gd = lambWReturn / ((1 + lambWReturn)*RS);

  // }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::lambertWLinearReverseBias
// Purpose       : Function to maintain continuity between forward bias and
//                 breakdown voltages of diode
// Special Notes :
// Scope         : public
// Creator       : Nick Johnson, Summer Intern
// Creation Date : 7/30/02
//-----------------------------------------------------------------------------
bool Instance::lambertWLinearReverseBias(double Isat, double Vte, double RS)
{
  double FF1 = (Vd + tBrkdwnV)/(-3.0 * Vte + tBrkdwnV);
  double FF2 = (1/2 - FF1)*(-2);
  double arg = Vte/RS;
  double arg1 = (Isat/arg - 3);
  double arg2 = FF1*arg1;
  arg1 = Xycemin(CONSTMAX_EXP_ARG, arg2);
  double evd = exp(arg2);
  double lambWArg = Isat*evd/arg;
  double lambWReturn;
  int ierr;
  double lambWError;
  devSupport.lambertw(lambWArg, lambWReturn, ierr, lambWError);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    if (ierr == 0)
      Xyce::dout() << "Safe LambertW return" << std::endl;
    else if (ierr == 1)
      Xyce::dout() << "LambertW argument not in domain" << std::endl;
    else
      Xyce::dout() << "Arithmetic problems with LambertW" << std::endl;
  }

  Id = -Isat*FF1 + FF2*lambWReturn*arg;

  double GdFF1 = 1/(-3.0*Vte + tBrkdwnV);
  double GdFF2 = 2 * GdFF1;
  double GdW = arg*lambWReturn*GdFF1*arg1/(1 + lambWReturn);
  Gd = -Isat*GdFF1 + GdFF2*arg*lambWReturn + FF2*GdW;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "lambWArg   = " << lambWArg << std::endl;
    Xyce::dout() << "lambWError = " << lambWError << std::endl;
    Xyce::dout() << "lambWReturn= " << lambWReturn << std::endl;
    Xyce::dout() << "Id         = " << Id << std::endl;
    Xyce::dout() << "Gd         = " << Gd << std::endl;
    Xyce::dout() << "Using lambertwReverseBias" << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::lambertWBreakdownCurrent
// Purpose       : Determine the diode breakdown current using Lambert-W function
// Special Notes :
// Scope         : public
// Creator       : Nick Johnson, Summer Intern
// Creation Date : 7/11/02
//-----------------------------------------------------------------------------
bool Instance::lambertWBreakdownCurrent(double Isat, double Vte, double RS)
{
  // with capacitor current and applying LambertW accros whole diode model
  /*  if (Cd != 0.0 && Icd != 0 && (lambertWFlag == 1))
      {
      double AA = Vte*(1 + RS*getSolverState().pdt*Cd);
      double XX = Icd - getSolverState().pdt*Qd;
      double arg1 = (RS*XX - Vd)/AA - tBrkdwnV/Vte;
      arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
      double evd = exp(arg1);
      double lambWArg = Isat*RS*evd/AA;
      double lambWReturn;
      int ierr;
      double lambWError;
      devSupport.lambertw(lambWArg, lambWReturn, ierr, lambWError);

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
      if (ierr == 0)
      Xyce::dout() << "Safe LambertW return" << std::endl;
      else if (ierr == 1)
      Xyce::dout() << "LambertW argument not in domain" << std::endl;
      else
      Xyce::dout() << "Arithmetic problems with LambertW" << std::endl;
      }
      #endif

      Id = -AA*lambWReturn/RS;
      Gd = (lambWReturn / (1 + lambWReturn)) * (1/RS);

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
      Xyce::dout() << "lambWArg   = " << lambWArg << std::endl;
      Xyce::dout() << "lambWError = " << lambWError << std::endl;
      Xyce::dout() << "Id         = " << Id << std::endl;
      Xyce::dout() << "Gd         = " << Gd << std::endl;
      Xyce::dout() << "Icd        = " << Icd << std::endl;
      Xyce::dout() << "Cd         = " << Cd  << std::endl;
      Xyce::dout() << "Qd         = " << Qd  << std::endl;
      Xyce::dout() << "AA         = " << AA << std::endl;
      Xyce::dout() << "XX         = " << XX << std::endl;
      Xyce::dout() << "Using lambertwBreakdown w/ capacitance" << std::endl;
      }
      #endif
      }

      // without capacitor current or applying lambertW just accros diode element
      else
      {*/
  double arg1 = (-Vd - tBrkdwnV)/ Vte;
  arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
  double evd = exp(arg1);
  double lambWArg = Isat*RS*evd/Vte;
  double lambWReturn;
  int ierr;
  double lambWError;
  devSupport.lambertw(lambWArg, lambWReturn, ierr, lambWError);

  Id = -Vte*lambWReturn/RS;
  Gd = lambWReturn / ((1 + lambWReturn)*RS);

  // }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Model::processParams ()
{
  //limit grading coeff
  if( M > 0.9 ) M = 0.9;

  //limit activation energy
  if( EG < 0.1 ) EG = 0.1;

  //limit depl cap coeff
  if( FC > 0.95 ) FC = 0.95;

  if( RS==0.0 )
    COND = 0.0;
  else
    COND = 1.0/RS;

  double xfc = log(1.0-FC);
  F2 = exp((1.0+M)*xfc);
  F3 = 1.0-FC*(1.0+M);

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
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    IS(1.0e-14),
    RS(0.0),
    COND(0.0),
    N(1.0),
    ISR(0.0),
    NR(2.0),
    IKF(0.0),
    TT(0.0),
    CJO(0.0),
    VJ(1.0),
    M(0.5),
    EG(1.11),
    XTI(3.0),
    TIKF(0.0),
    TBV1(0.0),
    TBV2(0.0),
    TRS1(0.0),
    TRS2(0.0),
    FC(0.5),
    BV(1e99),
    IBV(1.0e-10),
    IRF(1.0),
    NBV(1.0),
    IBVL(0.0),
    NBVL(1.0),
    F2(0.0),
    F3(0.0),
    TNOM(27),
    KF(0.0),
    AF(1.0),
    BVGiven(false),
    IRFGiven(false)
{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
  if (!given("TNOM"))
    TNOM = getDeviceOptions().tnom;

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

  // Note: Level 2 only params are: ISR, NR, IKF, NBV, IBVL, NBVL, TIKF, TBV1, TBV2, TRS1, TRS2
  if (getLevel() == 1)
  {
    std::string bad_parameters;

    if (given("ISR"))
      bad_parameters += " ISR";
    if (given("NR"))
      bad_parameters += " NR";
    if (given("IKF"))
      bad_parameters += " IKF";
    if (given("NBV"))
      bad_parameters += " NBV";
    if (given("IBVL"))
      bad_parameters += " IBVL";
    if (given("NBVL"))
      bad_parameters += " NBVL";
    if (given("TIKF"))
      bad_parameters += " TIKF";
    if (given("TBV1"))
      bad_parameters += " TBV1";
    if (given("TBV2"))
      bad_parameters += " TBV2";
    if (given("TRS1"))
      bad_parameters += " TRS1";
    if (given("TRS2"))
      bad_parameters += " TRS2";
    if (!bad_parameters.empty())
    {
      UserError(*this) << "Illegal parameter(s) given for level 1 diode:" << bad_parameters;
    }
  }

  processParams ();
}


//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Model::~Model()
{
  std::vector<Instance*>::iterator iterI;
  std::vector<Instance*>::iterator firstI = instanceContainer.begin ();
  std::vector<Instance*>::iterator lastI  = instanceContainer.end ();

  // loop over instances:
  for (iterI = firstI; iterI != lastI; ++iterI)
  {
    delete (*iterI);
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
    os << "AREA   = " << (*iter)->Area << std::endl;
    os << "  IC   = " << (*iter)->InitCond   << std::endl;
    os << "TEMP   = " << (*iter)->Temp << std::endl;
    os << " off   = " << (*iter)->off  << std::endl;

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


//-----------------------------------------------------------------------------
// Diode Master functions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  bool bsuccess = true;

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & di = *(*it);
    bool btmp = di.updateIntermediateVars ();
    stoVec[di.li_storevd] = di.Vd;
    bsuccess = bsuccess && btmp;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, double * leadF, double * leadQ, double * junctionV)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & di = *(*it);
    // load F:
    double Ir = di.Gspr * (di.Vp - di.Vpp);

    fVec[di.li_Pos] -= -Ir;
    fVec[di.li_Neg] -= di.Id;
    fVec[di.li_Pri] -= (-di.Id + Ir);

    // load Q:
    qVec[di.li_Neg] -= di.Qd;
    qVec[di.li_Pri] -= -di.Qd;

    // voltage limiter vectors.
    if( getDeviceOptions().voltageLimiterFlag )
    {
      double Vd_diff = di.Vd - di.Vd_orig;
      double Cd_Jdxp = -( di.Cd ) * Vd_diff;
      double Gd_Jdxp = -( di.Gd ) * Vd_diff;

      double * dFdxdVp = di.extData.dFdxdVpVectorRawPtr;
      // dFdxdVp vector
      dFdxdVp[di.li_Neg] += Gd_Jdxp;
      dFdxdVp[di.li_Pri] -= Gd_Jdxp;

      double * dQdxdVp = di.extData.dQdxdVpVectorRawPtr;
      // dQdxdVp vector
      dQdxdVp[di.li_Neg] += Cd_Jdxp;
      dQdxdVp[di.li_Pri] -= Cd_Jdxp;
    }

    if( di.loadLeadCurrent )
    {
      if (di.model_.CJO != 0.0)
      {
        leadQ[di.li_branch_data] = di.Qd;
      }
      leadF[di.li_branch_data] = di.Id;
      junctionV[di.li_branch_data] = solVec[di.li_Pos] - solVec[di.li_Neg];
 
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & di = *(*it);

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    // F-matrix:
    *di.fPosEquPosNodePtr += di.Gspr;
    *di.fPosEquPriNodePtr -= di.Gspr;
    *di.fNegEquNegNodePtr += di.Gd;
    *di.fNegEquPriNodePtr -= di.Gd;
    *di.fPriEquPosNodePtr -= di.Gspr;
    *di.fPriEquNegNodePtr -= di.Gd;
    *di.fPriEquPriNodePtr += di.Gspr + di.Gd;

    // Q-matrix:
    *di.qNegEquNegNodePtr += di.Cd;
    *di.qNegEquPriNodePtr -= di.Cd;
    *di.qPriEquNegNodePtr -= di.Cd;
    *di.qPriEquPriNodePtr += di.Cd;
#else
    // F-matrix:
    dFdx[di.li_Pos][di.APosEquPosNodeOffset] += di.Gspr;
    dFdx[di.li_Pos][di.APosEquPriNodeOffset] -= di.Gspr;
    dFdx[di.li_Neg][di.ANegEquNegNodeOffset] += di.Gd;
    dFdx[di.li_Neg][di.ANegEquPriNodeOffset] -= di.Gd;
    dFdx[di.li_Pri][di.APriEquPosNodeOffset] -= di.Gspr;
    dFdx[di.li_Pri][di.APriEquNegNodeOffset] -= di.Gd;
    dFdx[di.li_Pri][di.APriEquPriNodeOffset] += di.Gspr + di.Gd;

    // Q-matrix:
    dQdx[di.li_Neg][di.ANegEquNegNodeOffset] += di.Cd;
    dQdx[di.li_Neg][di.ANegEquPriNodeOffset] -= di.Cd;
    dQdx[di.li_Pri][di.APriEquNegNodeOffset] -= di.Cd;
    dQdx[di.li_Pri][di.APriEquPriNodeOffset] += di.Cd;
#endif
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
  static bool initialized = false;

  if (!initialized && (deviceMap.empty() ||
      ((deviceMap.find("D")!=deviceMap.end()) && 
      ((levelSet.find(1)!=levelSet.end()) || (levelSet.find(2)!=levelSet.end())))))
  {
    initialized = true;

    Config<Traits>::addConfiguration()
      .registerDevice("d", 1)
      .registerDevice("d", 2)
      .registerModelType("d", 1)
      .registerModelType("d", 2);
  }
}

//-----------------------------------------------------------------------------
// sensitivity related code.
//
// This is not the final implementation, as it is too klunky.  But, it works 
// for now as a preliminary experiment.


//-----------------------------------------------------------------------------
// Function      : processParams
// Purpose       : patterned after the model version of this function.
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename ScalarT> 
bool processParams (
    ScalarT & M,
    ScalarT & EG,
    ScalarT & FC,
    const ScalarT & RS,
    ScalarT & COND,
    ScalarT & F2,
    ScalarT & F3
    )
{
  //limit grading coeff
  if( M > 0.9 ) M = 0.9;

  //limit activation energy
  if( EG < 0.1 ) EG = 0.1;

  //limit depl cap coeff
  if( FC > 0.95 ) FC = 0.95;

  if( RS==0.0 )
    COND = 0.0;
  else
    COND = 1.0/RS;

  ScalarT xfc = log(1.0-FC);
  F2 = exp((1.0+M)*xfc);
  F3 = 1.0-FC*(1.0+M);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : updateTemperature
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename ScalarT> 
bool updateTemperature
  (
   const double & temp, 

   // instance variables/params:
   ScalarT & Temp,
   ScalarT & tJctCap,
   ScalarT & tJctPot,
   ScalarT & tDepCap,
   ScalarT & tF1,
   ScalarT & tSatCur,
   ScalarT & tSatCurR,
   ScalarT & tVcrit,
   ScalarT & tRS,
   ScalarT & tCOND,
   ScalarT & tIRF,
   ScalarT & tIKF,
   ScalarT & tBrkdwnV,

   // model variables/params:
   const ScalarT & TNOM,
   const ScalarT & VJ,
   const ScalarT & CJO,
   const ScalarT & M,
   const ScalarT & N,
   const ScalarT & IS,
   const ScalarT & EG,
   const ScalarT & XTI,
   const ScalarT & RS,
   const ScalarT & COND,
   const ScalarT & IRF,
   const ScalarT & NR,
   const ScalarT & IKF,
   const ScalarT & TIKF,
   const ScalarT & ISR,
   const ScalarT & IBV,
   const ScalarT & BV,

   const bool & BVGiven,
   const bool & IRFGiven,

   const ScalarT & TBV1,
   const ScalarT & TBV2,
   const ScalarT & TRS1,
   const ScalarT & TRS2,
   const ScalarT & FC,

   const int  level

   )
{
  ScalarT KoverQ = static_cast<ScalarT>(CONSTKoverQ);
  ScalarT REFTEMP = static_cast<ScalarT>(CONSTREFTEMP);
  ScalarT Eg0 = static_cast<ScalarT>(CONSTEg0);
  ScalarT alphaEg = static_cast<ScalarT>(CONSTalphaEg);
  ScalarT betaEg = static_cast<ScalarT>(CONSTbetaEg);
  ScalarT boltz = static_cast<ScalarT>(CONSTboltz);
  ScalarT Q = static_cast<ScalarT>(CONSTQ);
  ScalarT Eg300 = static_cast<ScalarT>(CONSTEg300);
  ScalarT root2 = static_cast<ScalarT>(CONSTroot2);
  ScalarT e = static_cast<ScalarT>(CONSTe);

  ScalarT vtnom = KoverQ * TNOM;

  ScalarT xfc = log( 1.0 - FC );

  if( temp != -999.0 ) Temp = temp;
  //double TNOM = TNOM;

  ScalarT vt = KoverQ * Temp;
  ScalarT fact2 = Temp / REFTEMP;
  ScalarT egfet = Eg0 - (alphaEg*Temp*Temp)/(Temp+betaEg);
  ScalarT arg = -egfet/(2.0*boltz*Temp) +
               Eg300/(boltz*(REFTEMP+REFTEMP));
  ScalarT pbfact = -2.0*vt*(1.5*log(fact2)+Q*arg);
  ScalarT egfet1 = Eg0 - (alphaEg*TNOM*TNOM)/
                  (TNOM+betaEg);
  ScalarT arg1 = -egfet1/(2.0*boltz*TNOM) +
		Eg300/(2.0*boltz*REFTEMP);
  ScalarT fact1 = TNOM/REFTEMP;
  ScalarT pbfact1 = -2.0*vtnom*(1.5*log(fact1)+Q*arg1);

  ScalarT pbo = (VJ-pbfact1)/fact1;
  ScalarT gmaold = (VJ-pbo)/pbo;

  tJctCap = CJO/
            (1.0+M*(4.0e-4*(TNOM-REFTEMP) -gmaold));

  tJctPot = pbfact+fact2*pbo;

  ScalarT gmanew = (tJctPot-pbo)/pbo;

  tJctCap *= 1.0 + M*(4.0e-4*(Temp-REFTEMP)-gmanew);

  tSatCur = IS*exp(((Temp/TNOM)-1.0)*
                          EG/(N*vt)+
                          XTI/N*log(Temp/TNOM));

  tF1 = tJctPot*(1.0-exp((1.0-M)*xfc))/(1.0-M);

  tDepCap = FC*tJctPot;

  ScalarT vte = N*vt;
  ScalarT tempBV;

  tVcrit = vte*log(vte/(root2*tSatCur));
  tRS   = RS;
  tCOND = COND;
  if (IRFGiven)
  {
    tIRF  = IRF*pow(fact2,1.6);
  }
  else
  {
    tIRF = IRF;
  }

  //int level = getLevel();
  if(level == 2)   // this section is PSPICE compatible
  {
    tSatCurR = ISR*exp((Temp/TNOM - 1.0)*
                         EG/(NR*vt)+
                         XTI/NR*log(Temp/TNOM));

    tIKF = IKF*(1 + TIKF*(Temp-TNOM));

    tempBV = BV*(1 + (Temp-TNOM)*
                        ( TBV1 + TBV2*(Temp-TNOM) ));

    tRS = RS*(1 + (Temp-TNOM)*
                     ( TRS1 + TRS2*(Temp-TNOM) ));

    tCOND = 0.0;
    if(tRS != 0.0) tCOND = 1.0/tRS;

    tJctPot = (VJ - egfet1)*fact2 - 3*vt*log(fact2) + egfet;

    tJctCap = CJO/(1.0 +
                          M*(4.0e-4*(Temp-TNOM) + (1-tJctPot/VJ)));
  }
  else
  {
    tempBV=BV;
  }

#ifdef Xyce_BREAKDOWN_ORIGINAL
  //Changed 9/18/07, K. R. Santarelli.  This is the original (broken) version
  //of the breakdown voltage computation.  It has two known issues:
  //
  //1.  It assumes N (the emission coefficient) is always 1 (division by
  //    vt instead of vte).
  //2.  While the for loop terminates due to the fact that it's implemented via
  //    a counter, the value of xbv does not converge in many cases, so the
  //    value of xbv upon terminating is often garbage.
  //
  //For consistency with old, existing simulations, this version of the
  //breakdown voltage calculation (along with the appropriate code for
  //computing the breakdown current in the level 1 and level 2 models---see
  //the appropriate sections of
  //Instance::updateIntermediateVars) can be invoked by
  //specifiying the CPPFLAG "-DXyce_BREAKDOWN_ORIGINAL" when creating a Xyce
  //Build.  The default, however, is the code which follows #else.

  double reltol = 1.0e-3;
  if( BVGiven )
  {
    ScalarT IRfactor = tIRF;
    ScalarT cbv = IBV;
    ScalarT xbv, xcbv;
    if( cbv < IRfactor*tSatCur*tempBV/vt )
    {
      cbv = IRfactor*tSatCur*tempBV/vt;
      xbv = tempBV;
    }
    else
    {
      ScalarT tol = reltol*cbv;
      xbv = tempBV-vt*log(1.0+cbv/(IRfactor*tSatCur));
      for( int i = 0; i < 25; ++i )
      {
        xbv = tempBV-vt*log(cbv/(IRfactor*tSatCur)+1.0-xbv/vt);
        xcbv = IRfactor*tSatCur*(exp((tempBV-xbv)/vt)-1.0+xbv/vt);
        if(fabs(xcbv-cbv)<=tol) break;
      }
    }
    tBrkdwnV = xbv;
  }
#else
  if( BVGiven)
  {
    ScalarT IRFactor=tIRF;
    ScalarT cbv = IBV;
    ScalarT xbv;
    ScalarT cthreshlow; //lower threshold for IBV
    ScalarT cthreshhigh; //high threshold for IBV
    int iter_count;
    const int ITER_COUNT_MAX=8;
    ScalarT arg2;

    ScalarT arg1=3.0*vte/(e*tempBV);
    arg1=arg1*arg1*arg1;
    cthreshlow=tSatCur*IRFactor*(1-arg1);
    cthreshhigh=tSatCur*IRFactor*(1-1.0/(e*e*e)) *
                exp(-1.0*(3.0*vte-tempBV)/vte);

    if(cbv >= cthreshhigh)
    {                     //if IBV is too high, tBrkdwnV will go below 3NVt.
      tBrkdwnV=3.0*vte;   //Clip tBrkdwnV to 3*N*Vt in this case (and hence
    }                     //clip IBV to cthreshhigh).


    else if(cbv <= cthreshlow)
    {                          //if IBV is too low, tBrkdwnV will go above
      tBrkdwnV=tempBV;  //BV.  Clip tBrkdwnV to BV in this case (and
    }                          //hence clip IBV to cthreshlow).


    //If IBV is in an acceptable range, perform a Picard iteration to find
    //tBrkdwnV, starting with an initial guess of tBrkdwnV=tempBV, and
    //running through the iteration ITER_COUNT_MAX times.

    else
    {
      xbv=tempBV;
      for(iter_count=0; iter_count < ITER_COUNT_MAX; iter_count++)
      {
        arg2=3.0*vte/(e*xbv);
        arg2=arg2*arg2*arg2;
        xbv=tempBV-vte*log(cbv/(tSatCur*IRFactor))+vte *
            log(1-arg2);
      }
      tBrkdwnV=xbv;
    }

//Note that, not only is tBrkdwnV adjusted using this method, but the effective
//value of Is (tSatCur) is adjusted as well.  The code used to use Is before,
//but now it will use Is*IRFactor*(1-(3*N*Vt/(e*tBrkdwnV))^3) instead.  This is
//reflected in changes made to Instance::updateIntermediateVars
//(for "normal" level 1 and 2 diode models) for the reverse breakdown region.
//Changes have not been made to the associated LambertW functions as of yet
//since there seem to be other issues involved with those functions independent
// of this change.
  }
#endif

  return true;
}


//-----------------------------------------------------------------------------
// Function      : updateIntermediateVars
// Purpose       : update intermediate variables for one diode instance
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename ScalarT> 
bool updateIntermediateVars ( 

  // inputs:
   const ScalarT & Vp,
   const ScalarT & Vpp,
   const ScalarT & Vn,
   const ScalarT & Vd,

  // instance params:
   const ScalarT & Temp,
   const ScalarT & tJctCap,
   const ScalarT & tJctPot,
   const ScalarT & tDepCap,
   const ScalarT & tF1,
   const ScalarT & tSatCur,
   const ScalarT & tSatCurR,
   const ScalarT & tVcrit,
   const ScalarT & tRS,
   const ScalarT & tCOND,
   const ScalarT & tIRF,
   const ScalarT & tIKF,
   const ScalarT & tBrkdwnV,

   // instance variables:
   const ScalarT & Area,
   const int & lambertWFlag,
   const double & gmin,

  // model params:
  const ScalarT M   , // grading parameter
  const ScalarT BV  , // breakdown voltage
  const ScalarT IBV , // reverse breakdown current
  const ScalarT NBV , // reverse breakdown ideality factor
  const ScalarT IBVL, // low-level reverse breakdown current
  const ScalarT NBVL, // low-level reverse breakdown ideality factor
  const ScalarT N   , // non-ideality factor.
  const ScalarT NR  , // emission coeff. for ISR.
  const ScalarT TT  , // transit time.
  const ScalarT F2  , // capacitive polynomial factor
  const ScalarT F3  , // capacitive polynomial factor

  const int  level,

  // outputs:
  ScalarT & Id,
  ScalarT & Gd,
  ScalarT & Qd,
  ScalarT & Cd,
  ScalarT & Gspr

    )
{
  bool bsuccess = true;



  ScalarT Inorm;   // normal diffusion current
  ScalarT Irec;    // recombination current
  ScalarT Kgen;    // generation factor
  ScalarT Khi;     // high-injection factor
  ScalarT Gd1, DKgen;
  ScalarT Gd2, DKhi;

  ScalarT Vc;

  ScalarT Icd = 0.0;// NOT USED, must be deprecated.

  ScalarT Isat  = tSatCur * Area;
  ScalarT IsatR = tSatCurR * Area;
  ScalarT KoverQ = static_cast<ScalarT>(CONSTKoverQ);
  ScalarT Vt    = KoverQ * Temp;
  ScalarT Vte   = N * Vt;
  ScalarT VteR  = NR * Vt;

  Gspr = tCOND * Area;

  // Current and Conductivity
  //----------------------------------------------

  if(level == 1)
  {
    // Adjustment for linear portion of reverse current
    ScalarT IRfactor;
    //if(Vd >= 0) IRfactor = 1.0;  //error?  shouldn't this be Vd>=-3*Vte????
    if(Vd >= -3.0 * Vte) IRfactor = 1.0;  //changing this to be consistent
    else IRfactor = tIRF;                 //with model in reference guide.
    Isat *= IRfactor;

    if (Vd >= -3.0 * Vte)
    {
      ScalarT arg1 = Vd / Vte;
      arg1 = Xycemin(static_cast<fadType>(CONSTMAX_EXP_ARG), arg1);
      ScalarT evd = exp(arg1);

      Id = Isat * (evd - 1.0) + gmin * Vd;
      Gd = Isat * evd / Vte + gmin;
    }
    // Linear reverse bias
    else if(!tBrkdwnV || (Vd >= -tBrkdwnV))
    {
      ScalarT arg = 3.0 * Vte / (Vd * CONSTe);
      arg = arg * arg * arg;
      Id = -Isat * (1.0 + arg) + gmin * Vd;
      Gd = Isat * 3.0 * arg / Vd + gmin;
    }
    // Reverse breakdown
    else
    {
      ScalarT arg1 = -(tBrkdwnV + Vd) / Vte;
      arg1 = Xycemin(static_cast<fadType>(CONSTMAX_EXP_ARG), arg1);
      ScalarT evrev = exp(arg1);

#ifdef Xyce_BREAKDOWN_ORIGINAL
      Id = -Isat * evrev + gmin * Vd;
      Gd = Isat * evrev / Vte + gmin;
#else
      //added by K. Santarelli 9/18/07 to account for change in tBrkdwnV
      //calculation.
      ScalarT arg2=3.0*Vte/(CONSTe*tBrkdwnV);
      arg2=arg2*arg2*arg2;
      ScalarT Isat_tBrkdwnV=Isat*(1-arg2);
      Id = -Isat_tBrkdwnV * evrev + gmin * Vd;
      Gd = Isat_tBrkdwnV * evrev / Vte + gmin;
#endif
    }
    Vc = Vd;
  }
  else if(level == 2)
  {
    // Adjustment for linear portion of reverse current
    ScalarT IRfactor;
    //if(Vd >= 0) IRfactor = 1.0;    //error?  Shouldn't this be Vd >= -3*Vte?
    if(Vd >= -3.0*Vte) IRfactor=1.0; //changing it to be consistent with model
    else IRfactor = tIRF;            //in the reference manual.
    Isat *= IRfactor;

    if (Vd >= -3.0 * Vte)
    {
      ScalarT arg1 = Vd / Vte;
      arg1 = Xycemin(static_cast<fadType>(CONSTMAX_EXP_ARG), arg1);
      ScalarT evd = exp(arg1);
      Inorm = Isat * (evd - 1.0) + gmin * Vd;
      Gd1 = Isat*evd/Vte + gmin;

      arg1 = Vd / VteR;
      arg1 = Xycemin(static_cast<fadType>(CONSTMAX_EXP_ARG), arg1);
      evd  = exp(arg1);
      Irec = IsatR * (evd - 1.0);
      Gd2  = IsatR*evd/VteR;

      Khi = 1;
      DKhi = 0;
      if(tIKF > 0)
      {
        Khi = sqrt(tIKF/(tIKF+Inorm));
        DKhi = 0.5*Khi*Gd1/(tIKF+Inorm);
      }
      Kgen = 0;
      DKgen = 0;
      if(Irec != 0)
      {
        Kgen = sqrt( pow(((1-Vd/tJctPot)*(1-Vd/tJctPot) + 0.005),M) );
        DKgen = -M*(1-Vd/tJctPot)/(tJctPot*Kgen);
      }

      Id = Inorm*Khi + Irec*Kgen;
      Gd = Gd1*Khi + Inorm*DKhi + Gd2*Kgen + Irec*DKgen;
    }
    // Linear reverse bias
    else if(!tBrkdwnV || (Vd >= -tBrkdwnV))
    {
      ScalarT arg = 3.0 * Vte / (Vd * CONSTe);
      arg = arg * arg * arg;
      Id = -Isat * (1.0 + arg) + gmin * Vd;
      Gd = Isat * 3.0 * arg / Vd + gmin;
    }
    // Reverse breakdown
    else
    {
      ScalarT arg1 = -(tBrkdwnV + Vd) / Vte;
      arg1 = Xycemin(static_cast<fadType>(CONSTMAX_EXP_ARG), arg1);
      ScalarT evrev = exp(arg1);

#ifdef Xyce_BREAKDOWN_ORIGINAL
      Id = -Isat * evrev + gmin * Vd;
      Gd = Isat * evrev / Vte + gmin;
#else
      //added 9/18/07 by K. Santarelli to account for change in tBrkdwnV
      //calculation.
      ScalarT arg2=3.0*Vte/(CONSTe*tBrkdwnV);
      arg2=arg2*arg2*arg2;
      ScalarT Isat_tBrkdwnV=Isat*(1-arg2);
      Id = -Isat_tBrkdwnV * evrev + gmin * Vd;
      Gd = Isat_tBrkdwnV * evrev / Vte + gmin;
#endif
    }
    Vc = Vd;
  }  // level

  // Only compute if Capacitance is non-zero
  //---------------------------------------
  if (tJctCap != 0.0)
  {
    // Charge Storage
    ScalarT Czero = tJctCap * Area;
    if (Vc < tDepCap)
    {
      //ScalarT arg = 1.0 - Vd/VJ;
      ScalarT arg = 1.0 - Vc / tJctPot;
      ScalarT arg1 = -M * log(arg);
      arg1 = Xycemin(static_cast<fadType>(CONSTMAX_EXP_ARG), arg1);
      ScalarT sarg = exp(arg1);

      //Qd = TT*Id + VJ*Czero*(1.0-arg*sarg)/(1.0-M);
      Qd = TT * Id + tJctPot * Czero * (1.0 - arg * sarg) / (1.0 - M);
      Cd = TT * Gd + Czero * sarg;
    }
    else
    {
      ScalarT Czof2 = Czero / F2;
      ScalarT MotJctPot = M / tJctPot;

      //Qd = TT*Id + Czero*tF1 + Czof2*(F3*(Vd-tDepCap)+(M/(VJ+VJ))
      Qd = TT * Id + Czero * tF1 +
           Czof2 * (F3 * (Vc - tDepCap) + (0.5 * MotJctPot) *
                    (Vc * Vc - tDepCap * tDepCap));
      //Cd = TT*Gd + Czof2*(F3+M*Vd/VJ);
      Cd = TT * Gd + Czof2 * (F3 + MotJctPot * Vc);
    }
  }
  else
  {
    Qd = 0.0;
    Cd = 0.0;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : applyLimiters
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename ScalarT> 
bool applyLimiters
( 
   DeviceSupport & devSupport,

   // inputs:
   const ScalarT & Vp,
   const ScalarT & Vpp,
   const ScalarT & Vn,

   // parameters:
   const ScalarT & tVcrit,
   const ScalarT & Vte,
   const ScalarT & BV,

   // output and book-keeping
   ScalarT & Vd,
   ScalarT & Vd_orig,
   ScalarT & Vd_old,

   const ScalarT & currVd_old, // (extData.currStoVectorRawPtr)[li_storevd];
   const ScalarT & nextVd_old, // (extData.nextStoVectorRawPtr)[li_storevd];

   const double InitCond,
   const bool InitCondGiven,
   const bool BVGiven,
   const int  off, // instance variable
   bool & origFlag, // instance variable
   const bool dotICapplies, // check all the "flagSol". if any == 1, then true

  // solver state variables:
   const int & newtonIter,
   const bool initJctFlag,
   const bool voltageLimiterFlag,
   const bool dcopFlag,
   const bool locaEnabledFlag
  )
{
  // Junction Voltage
  Vd = Vpp - Vn;
  Vd_orig = Vd;
  origFlag = true;

  // Setup initial junction conditions if UIC enabled
  //------------------------------------------------
  if (newtonIter == 0)
  {
    if (initJctFlag && voltageLimiterFlag && !dotICapplies)
    {
      if (InitCondGiven)
      {
        Vd = InitCond;
        origFlag = false;
      }
      else if (off)
      {
        Vd = 0.0;
        origFlag = false;
      }
      else
      {
        Vd=tVcrit;
        origFlag = false;
      }
    }

    // assume there is no history -- then check if the
    // state vector can overwrite this
    Vd_old = Vd;

    if (!(dcopFlag)|| (locaEnabledFlag && dcopFlag))
    {
      Vd_old = currVd_old;
    }
  }
  else  // just do this whenever it isn't the first iteration:
  {
    Vd_old = nextVd_old;
  }

  // Voltage limiting based on mode of diode
  //---------------------------------------
  if (voltageLimiterFlag)
  {
    int ichk = 0;

    if (newtonIter >= 0)
    {
      //Test if breakdown voltage given or not
      if (BVGiven && (Vd < Xycemin(0.0, -BV + 10.0 * Vte)))
      {
        double Vdtmp = -( BV + Vd );
        Vdtmp = devSupport.pnjlim(Vdtmp, -(Vd_old+BV), Vte, tVcrit, &ichk);
        Vd    = -(Vdtmp + BV);
      }
      else
        Vd = devSupport.pnjlim(Vd, Vd_old, Vte, tVcrit, &ichk);

      if (ichk) origFlag = false;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : diodeSensitivity::operator
// Purpose       : produces df/dp and dq/dp, where p= any diode parameter.  
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 7/18/2014
//-----------------------------------------------------------------------------
void diodeSensitivity::operator()(
    const ParameterBase &entity,
    const std::string &name,
    std::vector<double> & dfdp, 
    std::vector<double> & dqdp, 
    std::vector<double> & dbdp, 
    std::vector<int> & Findices,
    std::vector<int> & Qindices,
    std::vector<int> & Bindices
    ) const
{
  const ParameterBase * e1 = &entity;
  const Model & mod = *(dynamic_cast<const Model *> (e1));

  int sizeInstance = mod.instanceContainer.size();

  dfdp.resize(3*sizeInstance);
  dqdp.resize(3*sizeInstance);
  Findices.resize(3*sizeInstance);
  Qindices.resize(3*sizeInstance);

  // model params:
  fadType TNOM	 = mod.TNOM;
  fadType VJ	 = mod.VJ;
  fadType CJO	 = mod.CJO;
  fadType M	 = mod.M;
  fadType N	 = mod.N;
  fadType IS	 = mod.IS;
  fadType EG	 = mod.EG;
  fadType XTI	 = mod.XTI;
  fadType RS	 = mod.RS;
  fadType COND	 = mod.COND;
  fadType IRF	 = mod.IRF;
  fadType NR	 = mod.NR;
  fadType IKF	 = mod.IKF;
  fadType TIKF	 = mod.TIKF;
  fadType ISR	 = mod.ISR;
  fadType IBV	 = mod.IBV;
  fadType IBVL	 = mod.IBVL;
  fadType NBV	 = mod.NBV;
  fadType NBVL	 = mod.NBVL;
  fadType BV	 = mod.BV;
  fadType TT	 = mod.TT;
  fadType F2	 = mod.F2;
  fadType F3	 = mod.F3;

  fadType TBV1 = mod.TBV1;
  fadType TBV2 = mod.TBV2;
  fadType TRS1 = mod.TRS1;
  fadType TRS2 = mod.TRS2;
  fadType FC = mod.FC;
  fadType KF = mod.KF;
  fadType AF = mod.AF;

  std::string paramName = ExtendedString( name ).toUpper();

  if      (paramName=="VJ")   { VJ.diff(0,1);   }
  else if (paramName=="CJO")  { CJO.diff(0,1);  }
  else if (paramName=="CJ")   { CJO.diff(0,1);  } // CJ is a synonym for CJO
  else if (paramName=="CJ0")  { CJO.diff(0,1);  } // CJ0 is a synonym for CJO
  else if (paramName=="M")    { M.diff(0,1);    }
  else if (paramName=="N")    { N.diff(0,1);    }
  else if (paramName=="IS")   { IS.diff(0,1);   }
  else if (paramName=="JS")   { IS.diff(0,1);   } // JS is a synonym for IS
  else if (paramName=="EG")   { EG.diff(0,1);   }
  else if (paramName=="XTI")  { XTI.diff(0,1);  }
  else if (paramName=="RS")   { RS.diff(0,1);   }
  else if (paramName=="COND") { COND.diff(0,1); }
  else if (paramName=="IRF")  { IRF.diff(0,1);  }
  else if (paramName=="NR")   { NR.diff(0,1);   }
  else if (paramName=="IKF")  { IKF.diff(0,1);  }
  else if (paramName=="TIKF") { TIKF.diff(0,1); }
  else if (paramName=="ISR")  { ISR.diff(0,1);  }
  else if (paramName=="IBV")  { IBV.diff(0,1);  }
  else if (paramName=="IBVL") { IBVL.diff(0,1); }
  else if (paramName=="NBV")  { NBV.diff(0,1);  }
  else if (paramName=="NBVL") { NBVL.diff(0,1); }
  else if (paramName=="BV")   { BV.diff(0,1);   }
  else if (paramName=="VB")   { BV.diff(0,1);   } // VB is a synonym for BV
  else if (paramName=="TT")   { TT.diff(0,1);   }
  else if (paramName=="FC")   { FC.diff(0,1);   }
  else if (paramName=="KF")   { KF.diff(0,1);   }
  else if (paramName=="AF")   { AF.diff(0,1);   }
  else if (paramName=="TNOM") { TNOM.diff(0,1); }

  else if (paramName=="TBV1") { TBV1.diff(0,1); }
  else if (paramName=="TBV2") { TBV2.diff(0,1); }
  else if (paramName=="TRS1") { TRS1.diff(0,1); }
  else if (paramName=="TRS2") { TRS2.diff(0,1); }

  processParams( M, EG, FC, RS, COND, F2, F3);

  int inst=0;

  for (std::vector<Instance *>::const_iterator in = mod.instanceContainer.begin(); 
      in != mod.instanceContainer.end(); ++in, ++inst)
  {
    double * solVec = (*in)->extData.nextSolVectorRawPtr;

    fadType Vpp = solVec[(*in)->li_Pri];
    fadType Vn  = solVec[(*in)->li_Neg];
    fadType Vp  = solVec[(*in)->li_Pos];
    fadType Vd = Vpp - Vn;

    // instance variables:
    fadType tJctCap = 0.0;
    fadType tJctPot = 0.0;
    fadType tDepCap = 0.0;
    fadType tF1 = 0.0;
    fadType tSatCur = 0.0;
    fadType tSatCurR = 0.0;
    fadType tVcrit = 0.0;
    fadType tRS = 0.0;
    fadType tCOND = 0.0;
    fadType tIRF = 0.0;
    fadType tIKF = 0.0;
    fadType tBrkdwnV = 0.0;

    // instance parameters:
    fadType Temp = (*in)->Temp;
    fadType Area = (*in)->Area;

    updateTemperature(
       (*in)->Temp,
       Temp, tJctCap, tJctPot, tDepCap, tF1, tSatCur, tSatCurR, 
       tVcrit, tRS, tCOND, tIRF, tIKF, tBrkdwnV,
       TNOM, VJ, CJO, M, N, IS, EG, XTI, RS, COND, IRF, 
       NR, IKF, TIKF, ISR, IBV, BV,
       mod.BVGiven,
       mod.IRFGiven,
       TBV1, TBV2, TRS1, TRS2, FC,
       mod.getLevel()
     );

    fadType Id = 0.0;
    fadType Gd = 0.0;
    fadType Qd = 0.0;
    fadType Cd = 0.0;
    fadType Gspr = 0.0;

    updateIntermediateVars (
      // inputs:
      Vp, Vpp, Vn, Vd,
      // instance params:
      Temp, tJctCap, tJctPot, tDepCap, tF1,
      tSatCur, tSatCurR, tVcrit, tRS, tCOND,
      tIRF, tIKF, tBrkdwnV,
      // instance variables:
      Area, (*in)->lambertWFlag, (*in)->getDeviceOptions().gmin,
      // model params:
      M, BV, IBV, NBV, IBVL, NBVL, N, NR, TT, F2, F3,
      mod.getLevel(),
      // outputs:
      Id, Gd, Qd, Cd, Gspr
    ); 

     int iPos=0+inst*3;
     int iNeg=1+inst*3;
     int iPri=2+inst*3;

     fadType Ir = Gspr * (Vp - Vpp);
     dfdp[iPos] -= -Ir.dx(0);
     dfdp[iNeg] -= Id.dx(0);
     dfdp[iPri] -= (-Id.dx(0) + Ir.dx(0));

     dqdp[iPos] = 0.0;
     dqdp[iNeg] -= Qd.dx(0);
     dqdp[iPri] -= -Qd.dx(0);

     Findices[iPos] = (*in)->li_Pos;
     Findices[iNeg] = (*in)->li_Neg;
     Findices[iPri] = (*in)->li_Pri;

     Qindices[iPos] = (*in)->li_Pos;
     Qindices[iNeg] = (*in)->li_Neg;
     Qindices[iPri] = (*in)->li_Pri;
  }
}

} // namespace Diode
} // namespace Device
} // namespace Xyce

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

//-------------------------------------------------------------------------
//
// Purpose        : Implement the MOSFET Level 3 static model
//
// Special Notes  : BADMOS3 option not supported.
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

// ---------- Xyce Includes ----------
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MOSFET3.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_DEV_Message.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_DEV_MOSFET1.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Math.h>
#include <N_ANP_NoiseData.h>

namespace Xyce {
namespace Device {

namespace MOSFET3 {

void Traits::loadInstanceParameters(ParametricData<MOSFET3::Instance> &p)
{
  p.addPar ("TEMP",0.0,&MOSFET3::Instance::temp)
   .setExpressionAccess(ParameterType::TIME_DEP)
   .setUnit(STANDARD)
   .setCategory(CAT_NONE)
   .setDescription("Device temperature"),

  p.addPar ("L",0.0,&MOSFET3::Instance::l)
   .setOriginalValueStored(true)
   .setUnit(U_METER)
   .setCategory(CAT_GEOMETRY)
   .setDescription("Channel length")
   .setLengthScaling(true);

  p.addPar ("W",0.0,&MOSFET3::Instance::w)
   .setOriginalValueStored(true)
   .setUnit(U_METER)
   .setCategory(CAT_GEOMETRY)
   .setDescription("Channel width")
   .setLengthScaling(true);

  p.addPar ("AD",0.0,&MOSFET3::Instance::drainArea)
   .setUnit(U_METER2)
   .setCategory(CAT_GEOMETRY)
   .setDescription("Drain diffusion area")
   .setAreaScaling(true);

  p.addPar ("AS",0.0,&MOSFET3::Instance::sourceArea)
   .setUnit(U_METER2)
   .setCategory(CAT_GEOMETRY)
   .setDescription("Source diffusion area")
   .setAreaScaling(true);

  p.addPar ("NRD",1.0,&MOSFET3::Instance::drainSquares)
   .setUnit(U_SQUARES)
   .setCategory(CAT_GEOMETRY)
   .setDescription("Multiplier for RSH to yield parasitic resistance of drain");

  p.addPar ("NRS",1.0,&MOSFET3::Instance::sourceSquares)
   .setUnit(U_SQUARES)
   .setCategory(CAT_GEOMETRY)
   .setDescription("Multiplier for RSH to yield parasitic resistance of source");

  p.addPar ("PD",0.0,&MOSFET3::Instance::drainPerimeter)
   .setUnit(U_METER)
   .setCategory(CAT_GEOMETRY)
   .setDescription("Drain diffusion perimeter")
   .setLengthScaling(true);

  p.addPar ("PS",0.0,&MOSFET3::Instance::sourcePerimeter)
   .setUnit(U_METER)
   .setCategory(CAT_GEOMETRY)
   .setDescription("Source diffusion perimeter")
   .setLengthScaling(true);

  p.addPar ("M",1.0,&MOSFET3::Instance::numberParallel)
   .setUnit(U_NONE)
   .setCategory(CAT_CONTROL)
   .setDescription("Multiplier for M devices connected in parallel");

  // Initial conditions
  p.addPar ("IC1",0.0,&MOSFET3::Instance::icVDS)
   .setGivenMember(&MOSFET3::Instance::IC_GIVEN)
   .setUnit(U_VOLT)
   .setCategory(CAT_INITIAL)
   .setDescription("Initial condition on Drain-Source voltage");

  p.addPar ("IC2",0.0,&MOSFET3::Instance::icVGS)
   .setGivenMember(&MOSFET3::Instance::IC_GIVEN)
   .setUnit(U_VOLT)
   .setCategory(CAT_INITIAL)
   .setDescription("Initial condition on Gate-Source voltage");

  p.addPar ("IC3",0.0,&MOSFET3::Instance::icVBS)
   .setGivenMember(&MOSFET3::Instance::IC_GIVEN)
   .setUnit(U_VOLT)
   .setCategory(CAT_INITIAL)
   .setDescription("Initial condition on Bulk-Source voltage");

  p.makeVector ("IC",3);

  // Set up non-double precision variables:
  p.addPar ("OFF",false, &MOSFET3::Instance::OFF)
   .setUnit(U_LOGIC)
   .setCategory(CAT_VOLT)
   .setDescription("Initial condition of no voltage drops across device");
}

void Traits::loadModelParameters(ParametricData<MOSFET3::Model> &p)
{
  // Set up double precision variables:
  p.addPar ("L",1e-4,&MOSFET3::Model::model_l)
   .setUnit(U_METER)
   .setCategory(CAT_GEOMETRY)
   .setDescription("Default channel length");

  p.addPar ("W",1e-4,&MOSFET3::Model::model_w)
   .setUnit(U_METER)
   .setCategory(CAT_GEOMETRY)
   .setDescription("Default channel width");

  p.addPar ("VTO",0.0,&MOSFET3::Model::vt0)
   .setUnit(U_VOLT)
   .setCategory(CAT_VOLT)
   .setDescription("Zero-bias threshold voltage");

  p.addPar ("VT0",0.0,&MOSFET3::Model::vt0)
   .setUnit(U_VOLT)
   .setCategory(CAT_VOLT)
   .setDescription("Zero-bias threshold voltage (alias for VTO)");

  p.addPar ("KP",2e-5,&MOSFET3::Model::transconductance)
   .setUnit(U_AMPVM2)
   .setCategory(CAT_PROCESS)
   .setDescription("Transconductance coefficient");

  p.addPar ("GAMMA",0.0,&MOSFET3::Model::gamma)
   .setUnit(U_VOLTH)
   .setCategory(CAT_PROCESS)
   .setDescription("Bulk threshold parameter");

  p.addPar ("PHI",0.6,&MOSFET3::Model::phi)
   .setUnit(U_VOLT)
   .setCategory(CAT_PROCESS)
   .setDescription("Surface potential");

  p.addPar ("RD",0.0,&MOSFET3::Model::drainResistance)
   .setExpressionAccess(ParameterType::MIN_RES)
   .setUnit(U_OHM)
   .setCategory(CAT_RES)
   .setDescription("Drain ohmic resistance");

  p.addPar ("RS",0.0,&MOSFET3::Model::sourceResistance)
   .setExpressionAccess(ParameterType::MIN_RES)
   .setUnit(U_OHM)
   .setCategory(CAT_RES)
   .setDescription("Source ohmic resistance");

  p.addPar ("CBD",0.0,&MOSFET3::Model::capBD)
   .setExpressionAccess(ParameterType::MIN_CAP)
   .setGivenMember(&MOSFET3::Model::capBDGiven)
   .setUnit(U_FARAD)
   .setCategory(CAT_CAP)
   .setDescription("Zero-bias bulk-drain p-n capacitance");

  p.addPar ("CBS",0.0,&MOSFET3::Model::capBS)
   .setExpressionAccess(ParameterType::MIN_CAP)
   .setGivenMember(&MOSFET3::Model::capBSGiven)
   .setUnit(U_FARAD)
   .setCategory(CAT_CAP)
   .setDescription("Zero-bias bulk-source p-n capacitance");

  p.addPar ("IS",1e-14, &MOSFET3::Model::jctSatCur)
   .setUnit(U_AMP)
   .setCategory(CAT_CURRENT)
   .setDescription("Bulk p-n saturation current");

  p.addPar ("PB",0.8,&MOSFET3::Model::bulkJctPotential)
   .setUnit(U_VOLT)
   .setCategory(CAT_VOLT)
   .setDescription("Bulk p-n bottom potential");

  p.addPar ("CGSO",0.0,&MOSFET3::Model::gateSourceOverlapCapFactor)
   .setUnit(U_FARADMM1)
   .setCategory(CAT_CAP)
   .setDescription("Gate-source overlap capacitance/channel width");

  p.addPar ("CGDO",0.0,&MOSFET3::Model::gateDrainOverlapCapFactor)
   .setUnit(U_FARADMM1)
   .setCategory(CAT_CAP)
   .setDescription("Gate-drain overlap capacitance/channel width");

  p.addPar ("CGBO",0.0,&MOSFET3::Model::gateBulkOverlapCapFactor)
   .setUnit(U_FARADMM1)
   .setCategory(CAT_CAP)
   .setDescription("Gate-bulk overlap capacitance/channel length");

  p.addPar ("RSH",0.0,&MOSFET3::Model::sheetResistance)
   .setUnit(U_OHM)
   .setCategory(CAT_RES)
   .setDescription("Drain,source diffusion sheet resistance");

  p.addPar ("CJ",0.0,&MOSFET3::Model::bulkCapFactor)
   .setGivenMember(&MOSFET3::Model::bulkCapFactorGiven)
   .setUnit(U_FARADMM2)
   .setCategory(CAT_CAP)
   .setDescription("Bulk p-n zero-bias bottom capacitance/area");

  p.addPar ("MJ",0.5,&MOSFET3::Model::bulkJctBotGradingCoeff)
   .setUnit(U_NONE)
   .setCategory(CAT_DOPING)
   .setDescription("Bulk p-n bottom grading coefficient");

  p.addPar ("CJSW",0.0,&MOSFET3::Model::sideWallCapFactor)
   .setGivenMember(&MOSFET3::Model::sideWallCapFactorGiven)
   .setUnit(U_FARADMM2)
   .setCategory(CAT_CAP)
   .setDescription("Bulk p-n zero-bias sidewall capacitance/area");

  p.addPar ("MJSW",0.33,&MOSFET3::Model::bulkJctSideGradingCoeff)
   .setUnit(U_NONE)
   .setCategory(CAT_DOPING)
   .setDescription("Bulk p-n sidewall grading coefficient");

  p.addPar ("JS",0.0,&MOSFET3::Model::jctSatCurDensity)
   .setUnit(U_AMPMM2)
   .setCategory(CAT_PROCESS)
   .setDescription("Bulk p-n saturation current density");

  p.addPar ("TOX",1e-7,&MOSFET3::Model::oxideThickness)
   .setOriginalValueStored(true)
   .setUnit(U_METER)
   .setCategory(CAT_GEOMETRY)
   .setDescription("Gate oxide thickness");

  p.addPar ("LD",0.0,&MOSFET3::Model::latDiff)
   .setUnit(U_METER)
   .setCategory(CAT_DOPING)
   .setDescription("Lateral diffusion length");

  p.addPar ("UO",600.0,&MOSFET3::Model::surfaceMobility)
   .setUnit(U_CMM2VM1SM1)
   .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
   .setDescription("Surface mobility");

  p.addPar ("U0",600.0,&MOSFET3::Model::surfaceMobility)
   .setUnit(U_CMM2VM1SM1)
   .setCategory(CAT_PROCESS)
   .setDescription("Surface mobility (alias for UO)");

  p.addPar ("FC",0.5,&MOSFET3::Model::fwdCapDepCoeff)
   .setUnit(U_NONE)
   .setCategory(CAT_CAP)
   .setDescription("Bulk p-n forward-bias capacitance coefficient");

  p.addPar ("NSUB",0.0,&MOSFET3::Model::substrateDoping)
   .setUnit(U_CMM3)
   .setCategory(CAT_DOPING)
   .setDescription("Substrate doping density");

  p.addPar ("NSS",0.0,&MOSFET3::Model::surfaceStateDensity)
   .setUnit(U_CMM2)
   .setCategory(CAT_PROCESS)
   .setDescription("Surface state density");

  p.addPar ("ETA",0.0,&MOSFET3::Model::eta)
   .setUnit(U_NONE)
   .setCategory(CAT_PROCESS)
   .setDescription("Static feedback");

  p.addPar ("DELTA",0.0,&MOSFET3::Model::delta)
   .setUnit(U_NONE)
   .setCategory(CAT_PROCESS)
   .setDescription("Width effect on threshold");

  p.addPar ("NFS",0.0,&MOSFET3::Model::fastSurfaceStateDensity)
   .setUnit(U_CMM2)
   .setCategory(CAT_PROCESS)
   .setDescription("Fast surface state density");

  p.addPar ("THETA",0.0,&MOSFET3::Model::theta)
   .setUnit(U_VOLTM1)
   .setCategory(CAT_PROCESS)
   .setDescription("Mobility modulation");

  p.addPar ("VMAX",0.0,&MOSFET3::Model::maxDriftVel)
   .setUnit(U_MSM1)
   .setCategory(CAT_PROCESS)
   .setDescription("Maximum drift velocity");

  p.addPar ("KAPPA",0.2,&MOSFET3::Model::kappa)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Saturation field factor");

  p.addPar ("XJ",0.0,&MOSFET3::Model::junctionDepth)
   .setUnit(U_METER)
   .setCategory(CAT_GEOMETRY)
   .setDescription("Metallurgical junction depth");

  p.addPar ("TNOM",27.0,&MOSFET3::Model::tnom)
   .setUnit(STANDARD)
   .setCategory(CAT_NONE)
   .setDescription("Parameter measurement temperature");

  p.addPar ("KF",0.0,&MOSFET3::Model::fNcoef)
   .setUnit(U_NONE)
   .setCategory(CAT_FLICKER)
   .setDescription("Flicker noise coefficient");

  p.addPar ("AF",1.0,&MOSFET3::Model::fNexp)
   .setUnit(U_NONE)
   .setCategory(CAT_FLICKER)
   .setDescription("Flicker noise exponent");

  // Set up non-double precision variables:
  p.addPar ("TPG",1,&MOSFET3::Model::gateType)
   .setUnit(U_NONE)
   .setCategory(CAT_MATERIAL)
   .setDescription("Gate material type (-1 = same as substrate,0 = aluminum,1 = opposite of substrate)");

  DeviceModel::initThermalModel(p);
}

std::vector< std::vector<int> > Instance::jacStamp_DC_SC;
std::vector< std::vector<int> > Instance::jacStamp_DC;
std::vector< std::vector<int> > Instance::jacStamp_SC;
std::vector< std::vector<int> > Instance::jacStamp;

std::vector<int> Instance::jacMap_DC_SC;
std::vector<int> Instance::jacMap_DC;
std::vector<int> Instance::jacMap_SC;
std::vector<int> Instance::jacMap;

std::vector< std::vector<int> > Instance::jacMap2_DC_SC;
std::vector< std::vector<int> > Instance::jacMap2_DC;
std::vector< std::vector<int> > Instance::jacMap2_SC;
std::vector< std::vector<int> > Instance::jacMap2;

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::applyScale
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/2021
//-----------------------------------------------------------------------------
bool Instance::applyScale ()
{
  // apply scale
  if (getDeviceOptions().lengthScale != 1.0)
  {
    if (given("L")) { l *= getDeviceOptions().lengthScale; }
    if (given("W")) { w *= getDeviceOptions().lengthScale; } 
    if (given("AS")) { sourceArea *= getDeviceOptions().lengthScale * getDeviceOptions().lengthScale ; } 
    if (given("AD")) { drainArea *= getDeviceOptions().lengthScale * getDeviceOptions().lengthScale ; } 
    if (given("PD")) { drainPerimeter *= getDeviceOptions().lengthScale; } 
    if (given("PS")) { sourcePerimeter *= getDeviceOptions().lengthScale; }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  // Set any non-constant parameter defaults:
  if (!given("TEMP"))
    temp = getDeviceOptions().temp.getImmutableValue<double>();
  if (!given("L"))
    l =model_.model_l;
  if (!given("W"))
    w = model_.model_w;

  if (!given("AD"))
    drainArea = getDeviceOptions().defad;
  if (!given("AS"))
    sourceArea = getDeviceOptions().defas;


  // now set the temperature related stuff
  // This *must* be done prior to the computations below, due to the
  // possibility of temperature interpolation of models.  
  updateTemperature(temp);
  
  // process source/drain series resistance (from mos1temp)
  if(model_.drainResistance != 0)
  {
    drainConductance = 1/model_.drainResistance;
  }
  else if (model_.given("RSH"))
  {
    if(model_.sheetResistance != 0)
    {
      drainConductance =
        1/(model_.sheetResistance*drainSquares);
    }
    else
    {
      drainConductance = 0;
    }
  }
  else
  {
    drainConductance = 0;
  }
  if(model_.sourceResistance != 0)
  {
    sourceConductance = 1/model_.sourceResistance;
  }
  else if (model_.given("RSH"))
  {
    if(model_.sheetResistance != 0)
    {
      sourceConductance =
        1/(model_.sheetResistance*sourceSquares);
    }
    else
    {
      sourceConductance = 0;
    }
  }
  else
  {
    sourceConductance = 0;
  }

  // calculate dependent (ie computed) params and check for errors:

  if(l - 2 * model_.latDiff <=0)
  {
    UserError(*this) << "Effective channel length less than zero.";
  }

  EffectiveLength=l - 2*model_.latDiff;
  GateSourceOverlapCap = model_.gateSourceOverlapCapFactor * w;
  GateDrainOverlapCap = model_.gateDrainOverlapCapFactor * w;
  GateBulkOverlapCap = model_.gateBulkOverlapCapFactor * EffectiveLength;
  OxideCap = model_.oxideCapFactor * EffectiveLength * w;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 3/21/01
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & IB,
   Model & Miter,
  const FactoryBlock &  factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    dNode(0),
    gNode(0),
    sNode(0),
    bNode(0),
    dNodePrime(0),
    sNodePrime(0),
    OFF(false),
    l(getDeviceOptions().defl),
    w(getDeviceOptions().defw),
    drainArea(getDeviceOptions().defad),
    sourceArea(getDeviceOptions().defas),
    drainSquares(1.0),
    sourceSquares(1.0),
    drainPerimeter(0.0),
    sourcePerimeter(0.0),
    sourceConductance(0.0),
    drainConductance(0.0),
  temp(getDeviceOptions().temp.getImmutableValue<double>()),
    numberParallel(1.0),
    tTransconductance(0.0),
    tSurfMob(0.0),
    tPhi(0.0),
    tVto(0.0),
    tSatCur(0.0),
    tSatCurDens(0.0),
    tCbd(0.0),
    tCbs(0.0),
    tCj(0.0),
    tCjsw(0.0),
    tBulkPot(0.0),
    tDepCap(0.0),
    tVbi(0.0),
    icVBS(0.0),
    icVDS(0.0),
    icVGS(0.0),
    von(0.0),
    vdsat(0.0),
    sourceVcrit(0.0),
    drainVcrit(0.0),
    cd(0.0),
    cbs(0.0),
    cbd(0.0),
    gmbs(0.0),
    gm(0.0),
    gds(0.0),
    gbd(0.0),
    gbs(0.0),
    capbd(0.0),
    capbs(0.0),
    Cbd(0.0),
    Cbdsw(0.0),
    Cbs(0.0),
    Cbssw(0.0),
    f2d(0.0),
    f3d(0.0),
    f4d(0.0),
    f2s(0.0),
    f3s(0.0),
    f4s(0.0),
    mode(1),
    mode_low(0.0),
    mode_high(0.0),
    limitedFlag(false),
    IC_GIVEN(false),
    //calculated quantities
    EffectiveLength(0),
    DrainSatCur(0),
    SourceSatCur(0),
    GateSourceOverlapCap(0),
    GateDrainOverlapCap(0),
    GateBulkOverlapCap(0),
    OxideCap(0),
    Isource(0.0),
    Idrain(0.0),
    // local indices
    li_Drain(-1),
    li_DrainPrime(-1),
    li_Source(-1),
    li_SourcePrime(-1),
    li_Gate(-1),
    li_Bulk(-1),
    // matrix and vector pointers:
    // jacobian:
    //  drain row
    ADrainEquDrainNodeOffset(-1),
    ADrainEquDrainPrimeNodeOffset(-1),
    //  gate row
    AGateEquGateNodeOffset(-1),
    AGateEquBulkNodeOffset(-1),
    AGateEquDrainPrimeNodeOffset(-1),
    AGateEquSourcePrimeNodeOffset(-1),
    //  source row
    ASourceEquSourceNodeOffset(-1),
    ASourceEquSourcePrimeNodeOffset(-1),
    //  bulk row
    ABulkEquGateNodeOffset(-1),
    ABulkEquBulkNodeOffset(-1),
    ABulkEquDrainPrimeNodeOffset(-1),
    ABulkEquSourcePrimeNodeOffset(-1),
    // drain' row
    ADrainPrimeEquDrainNodeOffset(-1),
    ADrainPrimeEquGateNodeOffset(-1),
    ADrainPrimeEquBulkNodeOffset(-1),
    ADrainPrimeEquDrainPrimeNodeOffset(-1),
    ADrainPrimeEquSourcePrimeNodeOffset(-1),
    // source' row
    ASourcePrimeEquGateNodeOffset(-1),
    ASourcePrimeEquSourceNodeOffset(-1),
    ASourcePrimeEquBulkNodeOffset(-1),
    ASourcePrimeEquDrainPrimeNodeOffset(-1),
    ASourcePrimeEquSourcePrimeNodeOffset(-1),

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // F-vector pointers:
  // V_d Row:
    f_DrainEquDrainNodePtr(0),             // a
    f_DrainEquDrainPrimeNodePtr(0),        // b
  // V_g Row:
    f_GateEquGateNodePtr(0),               // c
    f_GateEquBulkNodePtr(0),               // d
    f_GateEquDrainPrimeNodePtr(0),         // e
    f_GateEquSourcePrimeNodePtr(0),        // f
  // V_s Row:
    f_SourceEquSourceNodePtr(0),           // g
    f_SourceEquSourcePrimeNodePtr(0),      // h
  // V_b Row:
    f_BulkEquGateNodePtr(0),               // i
    f_BulkEquBulkNodePtr(0),               // j
    f_BulkEquDrainPrimeNodePtr(0),         // k
    f_BulkEquSourcePrimeNodePtr(0),        // l
  // V_d' Row:
    f_DrainPrimeEquDrainNodePtr(0),        // m
    f_DrainPrimeEquGateNodePtr(0),         // n
    f_DrainPrimeEquBulkNodePtr(0),         // o
    f_DrainPrimeEquDrainPrimeNodePtr(0),   // p
    f_DrainPrimeEquSourcePrimeNodePtr(0),  // q
  // V_s' Row:
    f_SourcePrimeEquGateNodePtr(0),        // r
    f_SourcePrimeEquSourceNodePtr(0),      // s
    f_SourcePrimeEquBulkNodePtr(0),        // t
    f_SourcePrimeEquDrainPrimeNodePtr(0),  // u
    f_SourcePrimeEquSourcePrimeNodePtr(0), // v
  // Q-vector pointers:
  // V_d Row:
    q_DrainEquDrainNodePtr(0),             // a
    q_DrainEquDrainPrimeNodePtr(0),        // b
  // V_g Row:
    q_GateEquGateNodePtr(0),               // c
    q_GateEquBulkNodePtr(0),               // d
    q_GateEquDrainPrimeNodePtr(0),         // e
    q_GateEquSourcePrimeNodePtr(0),        // f
  // V_s Row:
    q_SourceEquSourceNodePtr(0),           // g
    q_SourceEquSourcePrimeNodePtr(0),      // h
  // V_b Row:
    q_BulkEquGateNodePtr(0),               // i
    q_BulkEquBulkNodePtr(0),               // j
    q_BulkEquDrainPrimeNodePtr(0),         // k
    q_BulkEquSourcePrimeNodePtr(0),        // l
  // V_d' Row:
    q_DrainPrimeEquDrainNodePtr(0),        // m
    q_DrainPrimeEquGateNodePtr(0),         // n
    q_DrainPrimeEquBulkNodePtr(0),         // o
    q_DrainPrimeEquDrainPrimeNodePtr(0),   // p
    q_DrainPrimeEquSourcePrimeNodePtr(0),  // q
  // V_s' Row:
    q_SourcePrimeEquGateNodePtr(0),        // r
    q_SourcePrimeEquSourceNodePtr(0),      // s
    q_SourcePrimeEquBulkNodePtr(0),        // t
    q_SourcePrimeEquDrainPrimeNodePtr(0),  // u
    q_SourcePrimeEquSourcePrimeNodePtr(0), // v
#endif

    vbd(0.0),
    vbs(0.0),
    vgs(0.0),
    vds(0.0),
    vgs_orig(0.0),
    vds_orig(0.0),
    vbs_orig(0.0),
    vbd_orig(0.0),
    vgs_old(0.0),
    vds_old(0.0),
    vbs_old(0.0),
    vbd_old(0.0),
    capgs(0.0),
    qgs(0.0),
    capgd(0.0),
    qgd(0.0),
    capgb(0.0),
    qgb(0.0),
    qbd(0.0),
    qbs(0.0),
    // local indices
    li_store_vbd(-1),
    li_store_vbs(-1),
    li_store_vgs(-1),
    li_store_vds(-1),
    li_store_von(-1),
    li_store_gm (-1),
    li_branch_dev_id(-1),
    li_branch_dev_ig(-1),
    li_branch_dev_is(-1),
    li_branch_dev_ib(-1),
    li_state_capgs(-1),
    li_state_capgd(-1),
    li_state_capgb(-1),
    li_state_qgs(-1),
    li_state_qgd(-1),
    li_state_qgb(-1),
    li_state_qbd(-1),
    li_state_qbs(-1)
{
  numIntVars   = 2;
  numExtVars   = 4;

  setNumStoreVars(6);
  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 4;    // this is the space to allocate if lead current or power is needed.
  numStateVars = 8;

  devConMap.resize(4);
  devConMap[0] = 1;
  devConMap[1] = 2;
  devConMap[2] = 1;
  devConMap[3] = 3;

  if( jacStamp.empty() )
  {
    // stamp for RS!=0, RD!=0
    jacStamp_DC_SC.resize(6);
    jacStamp_DC_SC[0].resize(2);  // Drain row
    jacStamp_DC_SC[0][0]=0;       // d-d
    jacStamp_DC_SC[0][1]=4;       // d-d'
    jacStamp_DC_SC[1].resize(4);  // Gate row
    jacStamp_DC_SC[1][0]=1;       // g-g
    jacStamp_DC_SC[1][1]=3;       // g-b
    jacStamp_DC_SC[1][2]=4;       // g-d'
    jacStamp_DC_SC[1][3]=5;       // g-s'
    jacStamp_DC_SC[2].resize(2);  // Source row
    jacStamp_DC_SC[2][0]=2;       // s-s
    jacStamp_DC_SC[2][1]=5;       // s-s'
    jacStamp_DC_SC[3].resize(4);  // Bulk row
    jacStamp_DC_SC[3][0]=1;       // b-g
    jacStamp_DC_SC[3][1]=3;       // b-b
    jacStamp_DC_SC[3][2]=4;       // b-d'
    jacStamp_DC_SC[3][3]=5;       // b-s'
    jacStamp_DC_SC[4].resize(5);  // Drain' row
    jacStamp_DC_SC[4][0]=0;       // d'-d
    jacStamp_DC_SC[4][1]=1;       // d'-g
    jacStamp_DC_SC[4][2]=3;       // d'-b
    jacStamp_DC_SC[4][3]=4;       // d'-d'
    jacStamp_DC_SC[4][4]=5;       // d'-s'
    jacStamp_DC_SC[5].resize(5);  // Source' row
    jacStamp_DC_SC[5][0]=1;       // s'-g
    jacStamp_DC_SC[5][1]=2;       // s'-s
    jacStamp_DC_SC[5][2]=3;       // s'-b
    jacStamp_DC_SC[5][3]=4;       // s'-d'
    jacStamp_DC_SC[5][4]=5;       // s'-s'

    jacMap_DC_SC.clear();
    jacStampMap(jacStamp_DC_SC, jacMap_DC_SC, jacMap2_DC_SC,
                jacStamp_DC,    jacMap_DC, jacMap2_DC, 5, 2, 6);

    jacStampMap(jacStamp_DC_SC, jacMap_DC_SC, jacMap2_DC_SC,
                jacStamp_SC,    jacMap_SC, jacMap2_SC, 4, 0, 6);

    jacStampMap(jacStamp_DC, jacMap_DC, jacMap2_DC,
                jacStamp,    jacMap, jacMap2, 4, 0, 6);
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);


  // Calculate any parameters specified as expressions:

  updateDependentParameters();

  applyScale ();
  processParams ();

  // NOTE: The five lines below are exact duplicates of the last 5
  // lines of processParams.  Logic dictates that they should be
  // completely unnecessary.  Yet failing to do them *HERE* in the
  // Level 1 causes one test case in the Xyce regression suite to fail
  // randomly with timestep-too-small in parallel.  See joseki bug
  // 647.  These lines are left here for superstition, in case the same
  // problem exists in the level 3.
  //   TVR 16 Sep 2015
  EffectiveLength=l - 2*model_.latDiff;
  GateSourceOverlapCap = model_.gateSourceOverlapCapFactor * w;
  GateDrainOverlapCap = model_.gateDrainOverlapCapFactor * w;
  GateBulkOverlapCap = model_.gateBulkOverlapCapFactor * EffectiveLength;
  OxideCap = model_.oxideCapFactor * EffectiveLength * w;

  numIntVars = (((sourceConductance == 0.0)?0:1)+((drainConductance == 0.0) ? 0:1));
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

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                                          const std::vector<int> & extLIDVecRef )
{
  numIntVars = (((sourceConductance == 0.0)?0:1)+((drainConductance == 0.0) ? 0:1));

  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "  In Instance::register LIDs\n\n";
    Xyce::dout() << "  name             = " << getName() << std::endl;
    Xyce::dout() << "  number of internal variables: " << numIntVars << std::endl;
    Xyce::dout() << "  number of external variables: " << numExtVars << std::endl;
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Drain = extLIDVec[0];
  li_Gate = extLIDVec[1];
  li_Source = extLIDVec[2];
  li_Bulk = extLIDVec[3];

  int intLoc = 0;

  if( drainConductance )
    li_DrainPrime = intLIDVec[intLoc++];
  else
    li_DrainPrime = li_Drain;

  if( sourceConductance )
    li_SourcePrime = intLIDVec[intLoc];
  else
    li_SourcePrime = li_Source;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "\n variable local indices:\n";
    Xyce::dout() << "  li_Drain       = " << li_Drain << std::endl;
    Xyce::dout() << "  li_DrainPrime  = " << li_DrainPrime << std::endl;
    Xyce::dout() << "  li_Source      = " << li_Source << std::endl;
    Xyce::dout() << "  li_SourcePrime = " << li_SourcePrime << std::endl;
    Xyce::dout() << "  li_Gate        = " << li_Gate << std::endl;
    Xyce::dout() << "  li_Bulk        = " << li_Bulk << std::endl;

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
  if ( li_DrainPrime != li_Drain )
    addInternalNode(symbol_table, li_DrainPrime, getName(), "drainprime");

  if ( li_SourcePrime != li_Source )
    addInternalNode(symbol_table, li_SourcePrime, getName(), "sourceprime");

  if (loadLeadCurrent)
  {
    addBranchDataNode(symbol_table, li_branch_dev_id, getName(), "BRANCH_DD");
    addBranchDataNode(symbol_table, li_branch_dev_is, getName(), "BRANCH_DS");
    addBranchDataNode(symbol_table, li_branch_dev_ig, getName(), "BRANCH_DG");
    addBranchDataNode(symbol_table, li_branch_dev_ib, getName(), "BRANCH_DB");
  }
  addStoreNode(symbol_table,  li_store_gm, getName().getEncodedName() + ":gm");
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "  In Instance::registerStateLIDs\n\n";
    Xyce::dout() << "  name             = " << getName() << std::endl;
    Xyce::dout() << "  Number of State LIDs: " << numStateVars << std::endl;
  }

  // Copy over the global ID lists:
  staLIDVec = staLIDVecRef;

  int lid=0;

  li_state_qgs = staLIDVec[lid++];
  li_state_qgd = staLIDVec[lid++];
  li_state_qgb = staLIDVec[lid++];

  li_state_capgs = staLIDVec[lid++];
  li_state_capgd = staLIDVec[lid++];
  li_state_capgb = staLIDVec[lid++];

  li_state_qbd = staLIDVec[lid++];
  li_state_qbs = staLIDVec[lid++];


  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "  State local indices:" << std::endl;
    Xyce::dout() << std::endl;
    Xyce::dout() << "  li_state_qgs           = " << li_state_qgs ;
    Xyce::dout() << "  li_state_capgs         = " << li_state_capgs;
    Xyce::dout() << "  li_state_capgd         = " << li_state_capgd;
    Xyce::dout() << "  li_state_capgb         = " << li_state_capgb;
    Xyce::dout() << "  li_state_qgd           = " << li_state_qgd;
    Xyce::dout() << "  li_state_qgb           = " << li_state_qgb;
    Xyce::dout() << "  li_state_qbs           = " << li_state_qbs;
    Xyce::dout() << "  li_state_qbd           = " << li_state_qbd;
    Xyce::dout() << std::endl;
    Xyce::dout() << section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/11/11
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs( const std::vector<int> & stoLIDVecRef )
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());

  // Copy over the global ID lists:
  stoLIDVec = stoLIDVecRef;

  int lid=0;
  li_store_vbd = stoLIDVec[lid++];
  li_store_vbs = stoLIDVec[lid++];
  li_store_vgs = stoLIDVec[lid++];
  li_store_vds = stoLIDVec[lid++];
  li_store_von = stoLIDVec[lid++];
  li_store_gm  = stoLIDVec[lid++];
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
    li_branch_dev_id = branchLIDVecRef[0];
    li_branch_dev_ig = branchLIDVecRef[1];
    li_branch_dev_is = branchLIDVecRef[2];
    li_branch_dev_ib = branchLIDVecRef[3];
  }
}


//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, Computational Sciences
// Creation Date : 9/3/02
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  if( drainConductance != 0.0 && sourceConductance != 0.0 )
    return jacStamp_DC_SC;
  else if( drainConductance != 0.0 && sourceConductance == 0.0 )
    return jacStamp_DC;
  else if( drainConductance == 0.0 && sourceConductance != 0.0 )
    return jacStamp_SC;

  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, Computational Sciences
// Creation Date : 9/3/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  std::vector<int> map;
  std::vector< std::vector<int> > map2;

  if (drainConductance != 0.0)
  {
    if (sourceConductance != 0.0)
    {
      map = jacMap_DC_SC;
      map2 = jacMap2_DC_SC;
    }
    else
    {
      map = jacMap_DC;
      map2 = jacMap2_DC;
    }
  }
  else
  {
    if (sourceConductance != 0.0)
    {
      map = jacMap_SC;
      map2 = jacMap2_SC;
    }
    else
    {
      map = jacMap;
      map2 = jacMap2;
    }
  }
  ADrainEquDrainNodeOffset             = jacLIDVec[map[0]][map2[0][0]];
  ADrainEquDrainPrimeNodeOffset        = jacLIDVec[map[0]][map2[0][1]];

  AGateEquGateNodeOffset               = jacLIDVec[map[1]][map2[1][0]];
  AGateEquBulkNodeOffset               = jacLIDVec[map[1]][map2[1][1]];
  AGateEquDrainPrimeNodeOffset         = jacLIDVec[map[1]][map2[1][2]];
  AGateEquSourcePrimeNodeOffset        = jacLIDVec[map[1]][map2[1][3]];

  ASourceEquSourceNodeOffset           = jacLIDVec[map[2]][map2[2][0]];
  ASourceEquSourcePrimeNodeOffset      = jacLIDVec[map[2]][map2[2][1]];

  ABulkEquGateNodeOffset               = jacLIDVec[map[3]][map2[3][0]];
  ABulkEquBulkNodeOffset               = jacLIDVec[map[3]][map2[3][1]];
  ABulkEquDrainPrimeNodeOffset         = jacLIDVec[map[3]][map2[3][2]];
  ABulkEquSourcePrimeNodeOffset        = jacLIDVec[map[3]][map2[3][3]];

  ADrainPrimeEquDrainNodeOffset        = jacLIDVec[map[4]][map2[4][0]];
  ADrainPrimeEquGateNodeOffset         = jacLIDVec[map[4]][map2[4][1]];
  ADrainPrimeEquBulkNodeOffset         = jacLIDVec[map[4]][map2[4][2]];
  ADrainPrimeEquDrainPrimeNodeOffset   = jacLIDVec[map[4]][map2[4][3]];
  ADrainPrimeEquSourcePrimeNodeOffset  = jacLIDVec[map[4]][map2[4][4]];

  ASourcePrimeEquGateNodeOffset        = jacLIDVec[map[5]][map2[5][0]];
  ASourcePrimeEquSourceNodeOffset      = jacLIDVec[map[5]][map2[5][1]];
  ASourcePrimeEquBulkNodeOffset        = jacLIDVec[map[5]][map2[5][2]];
  ASourcePrimeEquDrainPrimeNodeOffset  = jacLIDVec[map[5]][map2[5][3]];
  ASourcePrimeEquSourcePrimeNodeOffset = jacLIDVec[map[5]][map2[5][4]];
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/06/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);

  // F-pointers:
  f_DrainEquDrainNodePtr             = 	  &(dFdx[li_Drain][ADrainEquDrainNodeOffset]);
  f_DrainEquDrainPrimeNodePtr        = 	  &(dFdx[li_Drain][ADrainEquDrainPrimeNodeOffset]);

  f_GateEquGateNodePtr               = 	  &(dFdx[li_Gate][AGateEquGateNodeOffset]);
  f_GateEquBulkNodePtr               = 	  &(dFdx[li_Gate][AGateEquBulkNodeOffset]);
  f_GateEquDrainPrimeNodePtr         = 	  &(dFdx[li_Gate][AGateEquDrainPrimeNodeOffset]);
  f_GateEquSourcePrimeNodePtr        = 	  &(dFdx[li_Gate][AGateEquSourcePrimeNodeOffset]);

  f_SourceEquSourceNodePtr           = 	  &(dFdx[li_Source][ASourceEquSourceNodeOffset]);
  f_SourceEquSourcePrimeNodePtr      = 	  &(dFdx[li_Source][ASourceEquSourcePrimeNodeOffset]);

  f_BulkEquGateNodePtr               = 	  &(dFdx[li_Bulk][ABulkEquGateNodeOffset]);
  f_BulkEquBulkNodePtr               = 	  &(dFdx[li_Bulk][ABulkEquBulkNodeOffset]);
  f_BulkEquDrainPrimeNodePtr         = 	  &(dFdx[li_Bulk][ABulkEquDrainPrimeNodeOffset]);
  f_BulkEquSourcePrimeNodePtr        = 	  &(dFdx[li_Bulk][ABulkEquSourcePrimeNodeOffset]);

  f_DrainPrimeEquDrainNodePtr        = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset]);
  f_DrainPrimeEquGateNodePtr         = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset]);
  f_DrainPrimeEquBulkNodePtr         = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset]);
  f_DrainPrimeEquDrainPrimeNodePtr   = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]);
  f_DrainPrimeEquSourcePrimeNodePtr  = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]);

  f_SourcePrimeEquGateNodePtr        = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset]);
  f_SourcePrimeEquSourceNodePtr      = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset]);
  f_SourcePrimeEquBulkNodePtr        = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset]);
  f_SourcePrimeEquDrainPrimeNodePtr  = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]);
  f_SourcePrimeEquSourcePrimeNodePtr = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]);

  // Q-pointers:
  q_DrainEquDrainNodePtr             = 	  &(dQdx[li_Drain][ADrainEquDrainNodeOffset]);
  q_DrainEquDrainPrimeNodePtr        = 	  &(dQdx[li_Drain][ADrainEquDrainPrimeNodeOffset]);

  q_GateEquGateNodePtr               = 	  &(dQdx[li_Gate][AGateEquGateNodeOffset]);
  q_GateEquBulkNodePtr               = 	  &(dQdx[li_Gate][AGateEquBulkNodeOffset]);
  q_GateEquDrainPrimeNodePtr         = 	  &(dQdx[li_Gate][AGateEquDrainPrimeNodeOffset]);
  q_GateEquSourcePrimeNodePtr        = 	  &(dQdx[li_Gate][AGateEquSourcePrimeNodeOffset]);

  q_SourceEquSourceNodePtr           = 	  &(dQdx[li_Source][ASourceEquSourceNodeOffset]);
  q_SourceEquSourcePrimeNodePtr      = 	  &(dQdx[li_Source][ASourceEquSourcePrimeNodeOffset]);

  q_BulkEquGateNodePtr               = 	  &(dQdx[li_Bulk][ABulkEquGateNodeOffset]);
  q_BulkEquBulkNodePtr               = 	  &(dQdx[li_Bulk][ABulkEquBulkNodeOffset]);
  q_BulkEquDrainPrimeNodePtr         = 	  &(dQdx[li_Bulk][ABulkEquDrainPrimeNodeOffset]);
  q_BulkEquSourcePrimeNodePtr        = 	  &(dQdx[li_Bulk][ABulkEquSourcePrimeNodeOffset]);

  q_DrainPrimeEquDrainNodePtr        = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset]);
  q_DrainPrimeEquGateNodePtr         = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset]);
  q_DrainPrimeEquBulkNodePtr         = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset]);
  q_DrainPrimeEquDrainPrimeNodePtr   = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]);
  q_DrainPrimeEquSourcePrimeNodePtr  = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]);

  q_SourcePrimeEquGateNodePtr        = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset]);
  q_SourcePrimeEquSourceNodePtr      = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset]);
  q_SourcePrimeEquBulkNodePtr        = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset]);
  q_SourcePrimeEquDrainPrimeNodePtr  = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]);
  q_SourcePrimeEquSourcePrimeNodePtr = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]);
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
//  This is similar to the loadRHS function, only this function
//  only loads capacitor charges, and loads them into the daeQ vector.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/06/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  double * qVec = extData.daeQVectorRawPtr;
  double coef(0.0);

  double Qeqbs(0.0),Qeqbd(0.0),Qeqgb(0.0), Qeqgs(0.0), Qeqgd(0.0);
  //double ceqbs,ceqbd,ceqgb, ceqgs, ceqgd; // 3f5 vars
  int Dtype=model_.dtype;

  // do the same Dtype corrections on the charges that
  // are performed on the currents in the loadRHS function.

  // What is cbs and cbd?  They are diode currents (from exponentials),
  // so they are left out of this function.
  Qeqbs = Dtype*(qbs);
  Qeqbd = Dtype*(qbd);
  // These need "Dtype" here because we use them later *without*
  // Dtype, where SPICE uses it *with*
  Qeqgb = Dtype*(qgb);
  Qeqgs = Dtype*(qgs);
  Qeqgd = Dtype*(qgd);

  // 2 KCL for gate node
  coef = (Qeqgs+Qeqgd+Qeqgb);
  qVec[li_Gate] += coef*numberParallel;

  // 4 KCL for bulk node
  coef = Qeqbs + Qeqbd - Qeqgb;
  qVec[li_Bulk] += coef*numberParallel;

  // 5 KCL for drain' node
  coef = -(Qeqbd + Qeqgd);
  qVec[li_DrainPrime] += coef*numberParallel;

  // 6 KCL for source' node
  coef = -(Qeqbs + Qeqgs);
  qVec[li_SourcePrime] += coef*numberParallel;

  if( loadLeadCurrent )
  {
    double * leadQ = extData.nextLeadCurrQCompRawPtr;
    leadQ[li_branch_dev_id] = (-(Qeqbd + Qeqgd))*numberParallel;
    leadQ[li_branch_dev_is] = (-(Qeqbs + Qeqgs))*numberParallel;
    leadQ[li_branch_dev_ig] = (Qeqgs+Qeqgd+Qeqgb)*numberParallel;
    leadQ[li_branch_dev_ib] = (Qeqbs + Qeqbd - Qeqgb)*numberParallel;
  }

  // Same as for the loadRHS function, but with capacitive terms:
  //  gcgd = Capgd;
  //  gcgs = Capgs;
  //  gcgb = Capgb;
  //  gcbs = capbs;
  //  gcbd = capbd;
  if(!origFlag)
  {
    // The setup of gcgd, etc. is the same as in the loadDAEdQdxVector
    // function.
    double gcgd, gcgs, gcgb, gcbs, gcbd;
    if (getSolverState().tranopFlag || getSolverState().acopFlag || getSolverState().transientFlag)
    {
      gcgd = Capgd;
      gcgs = Capgs;
      gcgb = Capgb;
      // get at the two parasitic caps the same way
      gcbs = capbs;
      gcbd = capbd;
    }
    else
    {
      gcgd = 0.0; gcgs = 0.0; gcgb = 0.0; gcbs = 0.0; gcbd = 0.0;
    }

    // KCL 2
    double coef_Jdxp2 = Dtype*(gcgd*(vgd-vgd_orig)+gcgs*(vgs-vgs_orig)+
                        gcgb*(vgs-vgs_orig-vbs+vbs_orig));

    // 4 KCL for bulk node
    double coef_Jdxp4 = Dtype*(
                       - (gcgb)*(vgs-vgs_orig-vbs+vbs_orig)
                       + (gcgb)*(vbd-vbd_orig)
                       + (gcbs)*(vbs-vbs_orig));

    // 5 KCL for drain' node
    double coef_Jdxp5 = Dtype*(
                        -(gcgd)*(vgd-vgd_orig)
                        -(gcbd)*(vbd-vbd_orig));

    // 6 KCL for source' node
    double coef_Jdxp6 = Dtype*(-gcgs*(vgs-vgs_orig)-(gcbs)*(vbs-vbs_orig));

    double * dQdxdVp = extData.dQdxdVpVectorRawPtr;
    dQdxdVp[li_Gate       ] += coef_Jdxp2*numberParallel;
    dQdxdVp[li_Bulk       ] += coef_Jdxp4*numberParallel;
    dQdxdVp[li_DrainPrime ] += coef_Jdxp5*numberParallel;
    dQdxdVp[li_SourcePrime] += coef_Jdxp6*numberParallel;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 diode instance.
//
// Special Notes : See the special notes for loadDAEFVector.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/06/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * solVec = extData.nextSolVectorRawPtr;
  double * fVec = extData.daeFVectorRawPtr;
  double coef(0.0);
  double ceqbs(0.0),ceqbd(0.0),ceqgb(0.0), ceqgs(0.0), ceqgd(0.0); // 3f5 vars
  double gmin1 = getDeviceOptions().gmin;

  int Dtype=model_.dtype;

  // The next few lines are the same as for loadRHS, except
  // the capcitor current terms have been set to zero here.
  ceqbs = Dtype*(cbs);
  ceqbd = Dtype*(cbd);
  // These need "Dtype" here because we use them later *without*
  // Dtype, where SPICE uses it *with*
  ceqgb = 0.0;
  ceqgs = 0.0;
  ceqgd = 0.0;

  // 1 KCL for drain node
  if (drainConductance != 0.0)
  {
    coef = Idrain;
    fVec[li_Drain] += coef*numberParallel;
  }

  // 2 KCL for gate node
  coef = (ceqgs+ceqgd+ceqgb);
  fVec[li_Gate] += coef*numberParallel;

  // 3 KCL for source node
  if (sourceConductance != 0.0)
  {
    coef = Isource;
    fVec[li_Source] += coef*numberParallel;
  }

  // 4 KCL for bulk node
  coef = ceqbs + ceqbd - ceqgb;
  fVec[li_Bulk] += coef*numberParallel;

  // 5 KCL for drain' node
  coef = -Idrain-(ceqbd - cdreq + ceqgd);
  fVec[li_DrainPrime] += coef*numberParallel;

  // 6 KCL for source' node
  coef = -Isource-(ceqbs + cdreq + ceqgs);
  fVec[li_SourcePrime] += coef*numberParallel;

  // Same as for the loadRHS function, but without capacitive terms:
  //  gcgd = Capgd;
  //  gcgs = Capgs;
  //  gcgb = Capgb;
  //  gcbs = capbs;
  //  gcbd = capbd;
  if (!origFlag)
  {
    // 4 KCL for bulk node
    double coef_Jdxp4 = Dtype*(
                       + ((gbd-gmin1))*(vbd-vbd_orig)
                       + ((gbs-gmin1))*(vbs-vbs_orig));

    // 5 KCL for drain' node
    double coef_Jdxp5 = Dtype*(
                        -((gbd-gmin1))*(vbd-vbd_orig)
                        +gds*(vds-vds_orig)
                        +Gm*((mode>0)?(vgs-vgs_orig):(vgd-vgd_orig))
                        +Gmbs*((mode>0)?(vbs-vbs_orig):(vbd-vbd_orig)));

    // 6 KCL for source' node
    double coef_Jdxp6 = Dtype*(
                        -((gbs-gmin1))*(vbs-vbs_orig)
                        -gds*(vds-vds_orig)
                        -Gm*((mode>0)?(vgs-vgs_orig):(vgd-vgd_orig))
                        -Gmbs*((mode>0)?(vbs-vbs_orig):(vbd-vbd_orig)));

    double * dFdxdVp = extData.dFdxdVpVectorRawPtr;
    dFdxdVp[li_Bulk       ] += coef_Jdxp4*numberParallel;
    dFdxdVp[li_DrainPrime ] += coef_Jdxp5*numberParallel;
    dFdxdVp[li_SourcePrime] += coef_Jdxp6*numberParallel;
  }

  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    double * junctionV = extData.nextJunctionVCompRawPtr;
    leadF[li_branch_dev_id] = (-(ceqbd - cdreq + ceqgd))*numberParallel;
    leadF[li_branch_dev_is] = (-(ceqbs + cdreq + ceqgs))*numberParallel;
    leadF[li_branch_dev_ig] = (ceqgs+ceqgd+ceqgb)*numberParallel;
    leadF[li_branch_dev_ib] = (ceqbs + ceqbd - ceqgb)*numberParallel;

    junctionV[li_branch_dev_id] = solVec[li_Drain] - solVec[li_Source];
    junctionV[li_branch_dev_ig] = solVec[li_Gate] - solVec[li_Source];
    junctionV[li_branch_dev_is] = 0.0;
    junctionV[li_branch_dev_ib] = 0.0 ;
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
// Creation Date : 03/06/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);

  double gcgd(0.0);  // d(cqgd)/dVgd
  double gcgs(0.0);  // d(cqgs)/dVgs
  double gcgb(0.0);  // d(cqgb)/dVgb
  double gcbs(0.0);  // d(cqbs)/dVbs
  double gcbd(0.0);  // d(cqbd)/dVbd

  // get at the "conductances" for the gate capacitors with this trick
  //      gcgd = model_.dtype*Capgd;
  //      gcgs = model_.dtype*Capgs;
  //      gcgb = model_.dtype*Capgb;
  //
  //      In the loadRHS function, these would all be multiplied by
  //      getSolverState().pdt.  Here, for dQdx, the pdt term is left out.
  if (getSolverState().tranopFlag || getSolverState().acopFlag || getSolverState().transientFlag)
  {
    gcgd = Capgd;
    gcgs = Capgs;
    gcgb = Capgb;
    // get at the two parasitic caps the same way
    gcbs = capbs;
    gcbd = capbd;
  }

  dQdx[li_Gate][AGateEquGateNodeOffset] += (gcgd+gcgs+gcgb)
    *numberParallel;
  dQdx[li_Gate][AGateEquBulkNodeOffset] -= gcgb*numberParallel;
  dQdx[li_Gate][AGateEquDrainPrimeNodeOffset] -= gcgd*numberParallel;
  dQdx[li_Gate][AGateEquSourcePrimeNodeOffset] -= gcgs*numberParallel;

  dQdx[li_Bulk][ABulkEquGateNodeOffset] -= gcgb*numberParallel;
  dQdx[li_Bulk][ABulkEquBulkNodeOffset] += +(gcbs+gcbd+gcgb)
    *numberParallel;
  dQdx[li_Bulk][ABulkEquDrainPrimeNodeOffset] -= +gcbd*numberParallel;
  dQdx[li_Bulk][ABulkEquSourcePrimeNodeOffset] -=
    +gcbs*numberParallel;

  dQdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset] +=
    -gcgd*numberParallel;
  dQdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset] +=
    -gcbd*numberParallel;
  dQdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset] +=
    +(gcbd+gcgd)*numberParallel;

  dQdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset] -=
    gcgs*numberParallel;
  dQdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset] -=
    +gcbs*numberParallel;
  dQdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]+=
    +(gcbs+gcgs)*numberParallel;

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
// Creation Date : 03/06/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Drain][ADrainEquDrainNodeOffset] +=
    drainConductance*numberParallel;
  dFdx[li_Drain][ADrainEquDrainPrimeNodeOffset] -=
    drainConductance*numberParallel;

  dFdx[li_Source][ASourceEquSourceNodeOffset] +=
    sourceConductance*numberParallel;
  dFdx[li_Source][ASourceEquSourcePrimeNodeOffset] -=
    sourceConductance*numberParallel;

  dFdx[li_Bulk][ABulkEquBulkNodeOffset] +=
    (gbs+gbd)*numberParallel;
  dFdx[li_Bulk][ABulkEquDrainPrimeNodeOffset] -= gbd*numberParallel;
  dFdx[li_Bulk][ABulkEquSourcePrimeNodeOffset] -= gbs*numberParallel;

  dFdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset] -=
    drainConductance*numberParallel;
  dFdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset] +=
    Gm*numberParallel;
  dFdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset] +=
    (-gbd+Gmbs)*numberParallel;
  dFdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset] +=
    (drainConductance+gds+gbd+revsum)*numberParallel;
  dFdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset] +=
    (-gds-nrmsum)*numberParallel;

  dFdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset] -=
    Gm*numberParallel;
  dFdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset] -=
    sourceConductance*numberParallel;
  dFdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset] -=
    (gbs+Gmbs)*numberParallel;
  dFdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset] -=
    (gds+revsum)*numberParallel;
  dFdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset] +=
    (sourceConductance+gds+gbs+nrmsum)*numberParallel;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 03/01/01
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;
  // 3f5 likes to use the same variable names in local variables and in
  // structures.  Messes with us!  Define some local versions with capitals
  // instead
  double Von;
  double Vdsat;
  double Beta;
  //
  double evbs;
  double evbd;
  double sarg;
  double sargsw;
  //  double vgs1;
  //  double vgd1;
  //  double vgb1;
  double arg;
  int Check = 1;

  double capgs_old;
  double capgd_old;
  double capgb_old;

  // temporary storage for vgs/vgd/vds so that homotopy doesn't impact
  // voltage limiting "jdxp" terms

  double vgs_save;
  double vgd_save;
  double vds_save;

  // This is one of the vars that are set up at the top of mos3load that
  // should *not* be moved to the instance constructor!  tTransconductance
  // is set by updateTemperature.  Perhaps I should remove it from the
  // instance variables, too, since now it's pretty much a local thing.
  // same goes for the other things that depend on the t* variables!

  if( (tSatCurDens == 0) || (drainArea == 0) || (sourceArea == 0))
  {
    DrainSatCur = tSatCur;
    SourceSatCur = tSatCur;
  }
  else
  {
    DrainSatCur = tSatCurDens * drainArea;
    SourceSatCur = tSatCurDens * sourceArea;
  }
  Beta = tTransconductance * w/EffectiveLength;

  // don't do this anymore, use the one from the manager
  //  newtonIter = nlsMgrPtr->getNonLinearIter ();

  //  we need our solution variables for any of this stuff

  Vd = extData.nextSolVectorRawPtr[li_Drain];
  Vg = extData.nextSolVectorRawPtr[li_Gate];
  Vs = extData.nextSolVectorRawPtr[li_Source];
  Vb = extData.nextSolVectorRawPtr[li_Bulk];
  Vsp = extData.nextSolVectorRawPtr[li_SourcePrime];
  Vdp = extData.nextSolVectorRawPtr[li_DrainPrime];

  // now we need voltage drops
  Vddp = Vd - Vdp;
  Vssp = Vs - Vsp;
  Vbsp = Vb - Vsp;
  Vbdp = Vb - Vdp;
  Vgsp = Vg - Vsp;
  Vgdp = Vg - Vdp;
  Vgb  = Vg - Vb;
  Vdpsp = Vdp - Vsp;

  // Now the things that the 3f5 code really uses (from mos3load's
  // "general iteration" part  at lines 276-295
  vbs = model_.dtype * Vbsp;
  vgs = model_.dtype * Vgsp;
  vds = model_.dtype * Vdpsp;

  vbd = vbs-vds;
  vgd = vgs-vds;

  origFlag = 1;
  limitedFlag=false;
  vgs_orig = vgs;
  vds_orig = vds;
  vbs_orig = vbs;
  vbd_orig = vbd;
  vgd_orig = vgd;

  if (getSolverState().initJctFlag_ && !OFF && getDeviceOptions().voltageLimiterFlag)
  {
    if (IC_GIVEN)
    {
      vds = model_.dtype*icVDS;
      vgs = model_.dtype*icVGS;
      vbs = model_.dtype*icVBS;
      vbd = vbs - vds;
      vgd = vgs - vds;
      origFlag = false;
    }
    else
    {
      if (getSolverState().inputOPFlag)
      {
        Linear::Vector * flagSolVectorPtr = extData.flagSolVectorPtr;
        if ((*flagSolVectorPtr)[li_Drain] == 0 || (*flagSolVectorPtr)[li_Gate] == 0 ||
            (*flagSolVectorPtr)[li_Source] == 0 || (*flagSolVectorPtr)[li_SourcePrime] ||
            (*flagSolVectorPtr)[li_DrainPrime] || (*flagSolVectorPtr)[li_Bulk] )
        {
          vbs = -1;
          vgs = model_.dtype*tVto;
          vds = 0;
          vbd = vbs-vds;
          vgd = vgs-vds;
        }
      }
      else
      {
        vbs = -1;
        vgs = model_.dtype*tVto;
        vds = 0;
        vbd = vbs-vds;
        vgd = vgs-vds;
      }
    }
  }
  else if ((getSolverState().initFixFlag || getSolverState().initJctFlag_) && OFF)
  {
    vbs = vgs = vds = 0;
    vbd = vgd = 0;
  }

  if (getSolverState().newtonIter == 0)
  {

    if (!(getSolverState().dcopFlag)||(getSolverState().locaEnabledFlag && getSolverState().dcopFlag))
    {
      vbs_old = (*extData.currStoVectorPtr)[li_store_vbs];
      vbd_old = (*extData.currStoVectorPtr)[li_store_vbd];
      vgs_old = (*extData.currStoVectorPtr)[li_store_vgs];
      vds_old = (*extData.currStoVectorPtr)[li_store_vds];
      Von = model_.dtype *
                (*extData.currStoVectorPtr)[li_store_von];
    }
    else
    { // otherwise there is no history
      vbs_old = vbs;
      vbd_old = vbd;
      vgs_old = vgs;
      vds_old = vds;
      Von = 0.0;
    }
    vgd_old = vgs_old-vds_old;
  }
  else
  {
    vbs_old = (*extData.nextStoVectorPtr)[li_store_vbs];
    vbd_old = (*extData.nextStoVectorPtr)[li_store_vbd];
    vgs_old = (*extData.nextStoVectorPtr)[li_store_vgs];
    vds_old = (*extData.nextStoVectorPtr)[li_store_vds];
    Von = model_.dtype *
              (*extData.nextStoVectorPtr)[li_store_von];
    vgd_old = vgs_old-vds_old;
  }

  ////////////////////////////////////////////
  // SPICE-type Voltage Limiting
  ////////////////////////////////////////////
  if (getDeviceOptions().voltageLimiterFlag)
  {
    // Do not do limiting if mode initfix and OFF:
    if (! (getSolverState().initFixFlag && OFF))
    {
      if (vds_old >= 0)
      {
        vgs = devSupport.fetlim( vgs, vgs_old, Von);
        vds = vgs - vgd;
        vds = devSupport.limvds( vds,  vds_old);
        vgd = vgs - vds;
      }
      else
      {
        vgd = devSupport.fetlim( vgd, vgd_old, Von);
        vds = vgs - vgd;
        vds = -devSupport.limvds( -vds, -vds_old );
        vgs = vgd + vds;
      }

      if (vds >= 0.0)
      {
        vbs = devSupport.pnjlim( vbs, vbs_old, vt, sourceVcrit, &Check);
        vbd = vbs - vds;
      }
      else
      {
        vbd = devSupport.pnjlim( vbd, vbd_old, vt, drainVcrit, &Check);
        vbs = vbd + vds;
      }

      // for convergence:
      if (Check == 1) limitedFlag=true;

    }
  }

  ////
  // now all the preliminaries are over - we can start doing the
  //  real work
  ////
  vbd = vbs - vds;
  vgd = vgs - vds;
  Vgb = vgs - vbs;

  // Now set the origFlag
  if (vgs_orig != vgs || vds_orig != vds ||
      vbs_orig != vbs || vbd_orig != vbd || vgd_orig != vgd) origFlag = 0;


  ////
  //  bulk-source and bulk-drain diodes
  //  here we just evaluate the ideal diode current and the
  //   corresponding derivative (conductance).
  ////
  if(vbs <= 0)
  {
    gbs = SourceSatCur/vt;
    gbs += getDeviceOptions().gmin;
    cbs = gbs*vbs;
  }
  else
  {
    evbs = exp(std::min(CONSTMAX_EXP_ARG,vbs/vt));
    gbs = (SourceSatCur*evbs/vt + getDeviceOptions().gmin);
    cbs = (SourceSatCur * (evbs-1) + getDeviceOptions().gmin*vbs);
  }
  if(vbd <= 0)
  {
    gbd = DrainSatCur/vt;
    gbd += getDeviceOptions().gmin;
    cbd = gbd *vbd;
  }
  else
  {
    evbd = exp(std::min(CONSTMAX_EXP_ARG,vbd/vt));
    gbd = (DrainSatCur*evbd/vt + getDeviceOptions().gmin);
    cbd = (DrainSatCur *(evbd-1) + getDeviceOptions().gmin*vbd);
  }

  // 3f5 does this simple stuff
  if (vds >= 0)
    mode = 1;
  else
    mode = -1;

  {
    // here 3f5 seems to have inserted a c-translation of a fortran routine
    // Sadly, there are GOTO's and labeled statements all over this thing.
    // I'll do what I can to eliminate them all, but I'll start by
    // leaving it VERBATIM and ugly
    //
    // evaluate the drain current, its derivatives and the
    // charges associated with the gate, channel and bulk for mosfets based on
    // semi-empirical equations

    ////
    // * subroutine moseq3(vds,vbs,vgs,gm,gds,gmbs,
    // * qg,qc,qb,cggb,cgdb,cgsb,cbgb,cbdb,cbsb)
    ////

    ////
    //     *     this routine evaluates the drain current, its derivatives and
    //     *     the charges associated with the gate, channel and bulk
    //     *     for mosfets based on semi-empirical equations
    ////

    //
    //      common /mosarg/ vto,beta,gamma,phi,phib,cox,xnsub,xnfs,xd,xj,xld,
    //      1   xlamda,uo,uexp,vbp,utra,vmax,xneff,xl,xw,vbi,von,vdsat,qspof,
    //      2   beta0,beta1,cdrain,xqco,xqc,fnarrw,fshort,lev
    //      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
    //      1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
    //      2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
    //      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
    //      1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
    //      2   pivtol,pivrel
    //

    //    // equivalence (xlamda,alpha),(vbp,theta),(uexp,eta),(utra,xkappa)

    double coeff0 = 0.0631353e0;
    double coeff1 = 0.8013292e0;
    double coeff2 = -0.01110777e0;
    double oneoverxl;  //  1/effective length
    double eta; // eta from model after length factor */
    double phibs;   // phi - vbs */
    double sqphbs;  // square root of phibs */
    double dsqdvb;  //  */
    double sqphis;  // square root of phi */
    double sqphs3;  // square root of phi cubed */
    double wps;
    double oneoverxj;   // 1/junction depth */
    double xjonxl;  // junction depth/effective length */
    double djonxj;
    double wponxj;
    double arga;
    double argb;
    double argc;
    double dwpdvb;
    double dadvb;
    double dbdvb;
    double gammas;
    double fbodys;
    double fbody;
    double onfbdy;
    double qbonco;
    double vbix;
    double wconxj;
    double dfsdvb;
    double dfbdvb;
    double dqbdvb;
    double vth;
    double dvtdvb;
    double csonco;
    double cdonco;
    double dxndvb;
    double dvodvb;
    double dvodvd;
    double vgsx;
    double dvtdvd;
    double onfg;
    double fgate;
    double us;
    double dfgdvg;
    double dfgdvd;
    double dfgdvb;
    double dvsdvg;
    double dvsdvb;
    double dvsdvd;
    double xn;
    double vdsc;
    double onvdsc;
    double dvsdga;
    double vdsx;
    double dcodvb;
    double cdnorm;
    double cdo;
    double cd1;
    double fdrain;
    double fd2;
    double dfddvg;
    double dfddvb;
    double dfddvd;
    double gdsat;
    double cdsat;
    double gdoncd;
    double gdonfd;
    double gdonfg;
    double dgdvg;
    double dgdvd;
    double dgdvb;
    double emax;
    double emongd;
    double demdvg;
    double demdvd;
    double demdvb;
    double delxl;
    double dldvd;
    double dldem;
    double ddldvg;
    double ddldvd;
    double ddldvb;
    double dlonxl;
    double xlfact;
    double diddl;
    double gds0;
    double emoncd;
    double ondvt;
    double onxn;
    double wfact;
    double gms;
    double gmw;
    double fshort;


    // Begin block of mosfet continuation code.
    // This idea is based, loosely, on a paper by Jaijeet
    // Roychowdhury.
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "HOMOTOPY INFO: gainscale   = " << getSolverState().gainScale_ << std::endl
                   << "HOMOTOPY INFO: before vds  = " << vds << std::endl
                   << "HOMOTOPY INFO: before vgs  = " << vgs << std::endl;
    }

    // Save these before allowing homotopy to tweak them.  It
    // is important to restore them before moving on to
    // calculate RHS, because then the Jdxp terms will attempt to force
    // the external circuit to make these voltage drops the real thing!
    vds_save=vds;
    vgs_save=vgs;
    vgd_save=vgd;

    if (getSolverState().artParameterFlag_)
    {
      double alpha = getSolverState().gainScale_;
      double vgstConst = getDeviceOptions().vgstConst;

      vds = devSupport.contVds (vds, getSolverState().nltermScale_, getDeviceOptions().vdsScaleMin);

      if (mode==1)
      {
        vgs = devSupport.contVgst (vgs, alpha, vgstConst);
      }
      else
      {
        vgd = devSupport.contVgst (vgd, alpha, vgstConst);
      }
    }

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "HOMOTOPY INFO: after vds   = " << vds << std::endl;
      Xyce::dout() << "HOMOTOPY INFO: after vgs   = " << vgs << std::endl;
    }
    // End of block of mosfet continuation code.

    ////
    // bypasses the computation of charges
    ////

    ////
    // reference cdrain equations to source and
    // charge equations to bulk
    ////
    Vdsat = 0.0;
    oneoverxl = 1.0/EffectiveLength;

    // TVR: Massobrio and Antognetti call this "sigma"
    eta = model_.eta * 8.15e-22/(model_.oxideCapFactor*
                                       EffectiveLength*EffectiveLength*EffectiveLength);
    ////
    //.....square root term
    ////
    if ( (mode==1?vbs:vbd) <=  0.0 )
    {
      if (tPhi > 0)
        sqphis = sqrt(tPhi);
      else
        sqphis = 0;
      sqphs3 = tPhi*sqphis;
      phibs  =  tPhi-(mode==1?vbs:vbd);
      sqphbs  =  sqrt(phibs);
      dsqdvb  =  -0.5/sqphbs;
    }
    else
    {
      sqphis = sqrt(tPhi);
      sqphs3 = tPhi*sqphis;
      sqphbs = sqphis/(1.0+(mode==1?vbs:vbd)/
                       (tPhi+tPhi));
      phibs = sqphbs*sqphbs;
      dsqdvb = -phibs/(sqphs3+sqphs3);
    }
    ////
    //.....short channel effect factor
    ////
    if ( (model_.junctionDepth != 0.0) &&
         (model_.coeffDepLayWidth != 0.0) )
    {
      wps = model_.coeffDepLayWidth*sqphbs;
      oneoverxj = 1.0/model_.junctionDepth;
      xjonxl = model_.junctionDepth*oneoverxl;
      djonxj = model_.latDiff*oneoverxj;
      wponxj = wps*oneoverxj;
      wconxj = coeff0+coeff1*wponxj+coeff2*wponxj*wponxj;
      arga = wconxj+djonxj;
      argc = wponxj/(1.0+wponxj);
      argb = sqrt(1.0-argc*argc);
      fshort = 1.0-xjonxl*(arga*argb-djonxj);
      dwpdvb = model_.coeffDepLayWidth*dsqdvb;
      dadvb = (coeff1+coeff2*(wponxj+wponxj))*dwpdvb*oneoverxj;
      dbdvb = -argc*argc*(1.0-argc)*dwpdvb/(argb*wps);
      dfsdvb = -xjonxl*(dadvb*argb+arga*dbdvb);
    }
    else
    {
      fshort = 1.0;
      dfsdvb = 0.0;
    }
    ////
    //.....body effect
    //
    gammas = model_.gamma*fshort;
    fbodys = 0.5*gammas/(sqphbs+sqphbs);
    fbody = fbodys+model_.narrowFactor/w;

    onfbdy = 1.0/(1.0+fbody);
    dfbdvb = -fbodys*dsqdvb/sqphbs+fbodys*dfsdvb/fshort;
    qbonco =gammas*sqphbs+model_.narrowFactor*phibs/w;
    dqbdvb = gammas*dsqdvb+model_.gamma*dfsdvb*sqphbs-
      model_.narrowFactor/w;
    ////
    //.....static feedback effect
    //
    vbix = tVbi*model_.dtype-eta*(mode*vds);
    ////
    //.....threshold voltage
    //
    vth = vbix+qbonco;
    dvtdvd = -eta;
    dvtdvb = dqbdvb;
    ////
    //.....joint weak inversion and strong inversion
    //
    Von = vth;

    if ( model_.fastSurfaceStateDensity != 0.0 )
    {
      csonco = CONSTQ*model_.fastSurfaceStateDensity *
        1e4  *EffectiveLength*w/OxideCap;  // 1e4 to convert (cm**2/m**2)
      cdonco = qbonco/(phibs+phibs);
      xn = 1.0+csonco+cdonco;
      Von = vth+vt*xn;
      dxndvb = dqbdvb/(phibs+phibs)-qbonco*dsqdvb/(phibs*sqphbs);
      dvodvd = dvtdvd;
      dvodvb = dvtdvb+vt*dxndvb;
    }
    else
    {
      ////
      //.....cutoff region
      //

      if ( (mode==1?vgs:vgd) <= Von ) {
        cdrain = 0.0;
        gm = 0.0;
        gds = 0.0;
        gmbs = 0.0;
        goto innerline1000;
      }
    }
    ////
    //     *.....device is on
    ////

    vgsx = std::max((mode==1?vgs:vgd),Von);

    ////
    //.....mobility modulation by gate voltage
    ////
    onfg = 1.0+model_.theta*(vgsx-vth);
    fgate = 1.0/onfg;
    us = tSurfMob * 1e-4  *fgate;  // 1e4 to convert (m**2/cm**2)
    dfgdvg = -model_.theta*fgate*fgate;
    dfgdvd = -dfgdvg*dvtdvd;
    dfgdvb = -dfgdvg*dvtdvb;

    ////
    //     *.....saturation voltage
    ////
    Vdsat = (vgsx-vth)*onfbdy;

    if ( model_.maxDriftVel <= 0.0 )
    {
      dvsdvg = onfbdy;
      dvsdvd = -dvsdvg*dvtdvd;
      dvsdvb = -dvsdvg*dvtdvb-Vdsat*dfbdvb*onfbdy;
    }
    else
    {
      vdsc = EffectiveLength*model_.maxDriftVel/us;
      onvdsc = 1.0/vdsc;
      arga = (vgsx-vth)*onfbdy;
      argb = sqrt(arga*arga+vdsc*vdsc);
      Vdsat = arga+vdsc-argb;
      dvsdga = (1.0-arga/argb)*onfbdy;
      dvsdvg = dvsdga-(1.0-vdsc/argb)*vdsc*dfgdvg*onfg;
      dvsdvd = -dvsdvg*dvtdvd;
      dvsdvb = -dvsdvg*dvtdvb-arga*dvsdga*dfbdvb;
    }
    ////
    //     *.....current factors in linear region
    ////
    vdsx = std::min((mode*vds),Vdsat);
    if ( vdsx == 0.0 ) goto line900;

    cdo = vgsx-vth-0.5*(1.0+fbody)*vdsx;
    dcodvb = -dvtdvb-0.5*dfbdvb*vdsx;

    ////
    //     *.....normalized drain current
    ////
    cdnorm = cdo*vdsx;

    gm = vdsx;
    gds = vgsx-vth-(1.0+fbody+dvtdvd)*vdsx;
    gmbs = dcodvb*vdsx;
    ////
    //     *.....drain current without velocity saturation effect
    ////
    cd1 = Beta*cdnorm;
    Beta = Beta*fgate;
    cdrain = Beta*cdnorm;
    gm = Beta*gm+dfgdvg*cd1;
    gds = Beta*gds+dfgdvd*cd1;
    gmbs = Beta*gmbs;

    ////
    //     *.....velocity saturation factor
    ////
    if ( model_.maxDriftVel != 0.0 )
    {
      fdrain = 1.0/(1.0+vdsx*onvdsc);
      fd2 = fdrain*fdrain;
      arga = fd2*vdsx*onvdsc*onfg;
      dfddvg = -dfgdvg*arga;
      dfddvd = -dfgdvd*arga-fd2*onvdsc;
      dfddvb = -dfgdvb*arga;
      ////
      //       *.....drain current
      ////
      gm = fdrain*gm+dfddvg*cdrain;
      gds = fdrain*gds+dfddvd*cdrain;
      gmbs = fdrain*gmbs+dfddvb*cdrain;
      cdrain = fdrain*cdrain;
      Beta = Beta*fdrain;
    }

    ////
    //     *.....channel length modulation
    ////
    if ( (mode*vds) <= Vdsat ) goto line700;
    if ( model_.maxDriftVel <= 0.0 ) goto line510;
    if (model_.alpha == 0.0) goto line700;
    cdsat = cdrain;
    gdsat = cdsat*(1.0-fdrain)*onvdsc;
    gdsat = std::max(1.0e-12,gdsat);
    gdoncd = gdsat/cdsat;
    gdonfd = gdsat/(1.0-fdrain);
    gdonfg = gdsat*onfg;
    dgdvg = gdoncd*gm-gdonfd*dfddvg+gdonfg*dfgdvg;
    dgdvd = gdoncd*gds-gdonfd*dfddvd+gdonfg*dfgdvd;
    dgdvb = gdoncd*gmbs-gdonfd*dfddvb+gdonfg*dfgdvb;

    // to have this condition make sense, the option BADMOS3 must be
    // recognized on the options line.  A problem with the KAPPA parameter
    // was detected and fixed in 3f2.  BADMOS3 enables the unfixed deal,
    // just in case parameter fitting is affected.
    // It has never been implemented in Xyce
    //#ifdef BADMOS3_IMPLEMENTED
    //    if (ckt->CKTbadMos3)
    //      emax = cdsat*oneoverxl/gdsat;
    //    else
    //#endif
    emax = model_.kappa * cdsat*oneoverxl/gdsat;
    emoncd = emax/cdsat;
    emongd = emax/gdsat;
    demdvg = emoncd*gm-emongd*dgdvg;
    demdvd = emoncd*gds-emongd*dgdvd;
    demdvb = emoncd*gmbs-emongd*dgdvb;

    arga = 0.5*emax*model_.alpha;
    argc = model_.kappa*model_.alpha;
    argb = sqrt(arga*arga+argc*((mode*vds)-Vdsat));
    delxl = argb-arga;
    dldvd = argc/(argb+argb);
    dldem = 0.5*(arga/argb-1.0)*model_.alpha;
    ddldvg = dldem*demdvg;
    ddldvd = dldem*demdvd-dldvd;
    ddldvb = dldem*demdvb;
    goto line520;
  line510:
    delxl = sqrt(model_.kappa*((mode*vds)-Vdsat)*
                 model_.alpha);
    dldvd = 0.5*delxl/((mode*vds)-Vdsat);
    ddldvg = 0.0;
    ddldvd = -dldvd;
    ddldvb = 0.0;
    ////
    //     *.....punch through approximation
    ////
  line520:
    if ( delxl > (0.5*EffectiveLength) )
    {
      delxl = EffectiveLength-(EffectiveLength*EffectiveLength/
                               (4.0*delxl));
      arga = 4.0*(EffectiveLength-delxl)*(EffectiveLength-delxl)/
        (EffectiveLength*EffectiveLength);
      ddldvg = ddldvg*arga;
      ddldvd = ddldvd*arga;
      ddldvb = ddldvb*arga;
      dldvd =  dldvd*arga;
    }
    ////
    //     *.....saturation region
    ////
    dlonxl = delxl*oneoverxl;
    xlfact = 1.0/(1.0-dlonxl);
    cdrain = cdrain*xlfact;
    diddl = cdrain/(EffectiveLength-delxl);
    gm = gm*xlfact+diddl*ddldvg;
    gds0 = gds*xlfact+diddl*ddldvd;
    gmbs = gmbs*xlfact+diddl*ddldvb;
    gm = gm+gds0*dvsdvg;
    gmbs = gmbs+gds0*dvsdvb;
    gds = gds0*dvsdvd+diddl*dldvd;

    ////
    //     *.....finish strong inversion case
    ////
  line700:
    if ( (mode==1?vgs:vgd) < Von )
    {
      ////
      //       *.....weak inversion
      ////
      onxn = 1.0/xn;
      ondvt = onxn/vt;
      wfact = exp( ((mode==1?vgs:vgd)-Von)*ondvt );
      cdrain = cdrain*wfact;
      gms = gm*wfact;
      gmw = cdrain*ondvt;
      gm = gmw;

      if ((mode*vds) > Vdsat)
      {
        gm = gm+gds0*dvsdvg*wfact;
      }
      gds = gds*wfact+(gms-gmw)*dvodvd;
      gmbs = gmbs*wfact+(gms-gmw)*dvodvb-gmw*
        ((mode==1?vgs:vgd)-Von)*onxn*dxndvb;
    }
    ////
    //     *.....charge computation
    ////
    goto innerline1000;
    ////
    //     *.....special case of vds = 0.0d0
    ////
  line900:

    Beta = Beta*fgate;
    cdrain = 0.0;
    gm = 0.0;
    gds = Beta*(vgsx-vth);
    gmbs = 0.0;
    if ( (model_.fastSurfaceStateDensity != 0.0) &&
         ((mode==1?vgs:vgd) < Von) )
    {
      gds *=exp(((mode==1?vgs:vgd)-Von)/(vt*xn));
    }
  innerline1000:;
    ////
    //     *.....done
    ////
    // end of moseq3
    ////
  }

  // now deal with n vs p polarity

  von = model_.dtype * Von;
  vdsat = model_.dtype * Vdsat;

  ////
  //   *  COMPUTE EQUIVALENT DRAIN CURRENT SOURCE
  ////

  cd = mode * cdrain - cbd;
  // in 3f5 this is all in a block conditioned on CKTmode, but since
  // it's valid for MODETRAN and MODETRANOP we'll just always do it

  ////
  // * now we do the hard part of the bulk-drain and bulk-source
  // * diode - we evaluate the non-linear capacitance and
  // * charge
  // *
  // * the basic equations are not hard, but the implementation
  // * is somewhat long in an attempt to avoid log/exponential
  // * evaluations
  ////
  ////
  // *  charge storage elements
  // *
  // *.. bulk-drain and bulk-source depletion capacitances
  ////
  // I took out all the CAPBYPASS stuff, and the
  // unnecessary curly braces that wind up there if you do

  // can't bypass the diode capacitance calculations
  if(Cbs != 0 || Cbssw != 0 )
  {
    if (vbs < tDepCap)
    {
      arg=1-vbs/tBulkPot;
      ////
      // * the following block looks somewhat long and messy,
      // * but since most users use the default grading
      // * coefficients of .5, and sqrt is MUCH faster than an
      // * exp(log()) we use this special case code to buy time.
      // * (as much as 10% of total job time!)
      ////
      if(model_.bulkJctBotGradingCoeff==model_.bulkJctSideGradingCoeff)
      {
        if(model_.bulkJctBotGradingCoeff == .5)
        {
          sarg = sargsw = 1/sqrt(arg);
        }
        else
        {
          sarg = sargsw = exp(-model_.bulkJctBotGradingCoeff*log(arg));
        }
      }
      else
      {
        if(model_.bulkJctBotGradingCoeff == .5)
        {
          sarg = 1/sqrt(arg);
        }
        else
        {
          sarg = exp(-model_.bulkJctBotGradingCoeff*log(arg));
        }
        if(model_.bulkJctSideGradingCoeff == .5)
        {
          sargsw = 1/sqrt(arg);
        }
        else
        {
          sargsw =exp(-model_.bulkJctSideGradingCoeff* log(arg));
        }
      }
      qbs = tBulkPot*(Cbs*(1-arg*sarg)/(1-model_.bulkJctBotGradingCoeff)
                      +Cbssw*(1-arg*sargsw)/(1-model_.bulkJctSideGradingCoeff));
      capbs=Cbs*sarg+ Cbssw*sargsw;

    }
    else
    {
      qbs = f4s + vbs*(f2s+vbs*(f3s/2));
      capbs=f2s+f3s*vbs;
    }
  }
  else
  {
    qbs = 0;
    capbs=0;
  }

  //// can't bypass the diode capacitance calculations
  if(Cbd != 0 || Cbdsw != 0 )
  {

    if (vbd < tDepCap)
    {
      arg=1-vbd/tBulkPot;
      ////
      // * the following block looks somewhat long and messy,
      // * but since most users use the default grading
      // * coefficients of .5, and sqrt is MUCH faster than an
      // * exp(log()) we use this special case code to buy time.
      // * (as much as 10% of total job time!)
      ////
      if(model_.bulkJctBotGradingCoeff == .5 &&
         model_.bulkJctSideGradingCoeff == .5)
      {
        sarg = sargsw = 1/sqrt(arg);
      }
      else
      {
        if(model_.bulkJctBotGradingCoeff == .5)
        {
          sarg = 1/sqrt(arg);
        }
        else
        {
          sarg = exp(-model_.bulkJctBotGradingCoeff*log(arg));
        }
        if(model_.bulkJctSideGradingCoeff == .5)
        {
          sargsw = 1/sqrt(arg);
        }
        else
        {
          sargsw =exp(-model_.bulkJctSideGradingCoeff*log(arg));
        }
      }
      qbd =
        tBulkPot*(Cbd*
                  (1-arg*sarg)
                  /(1-model_.bulkJctBotGradingCoeff)
                  +Cbdsw*
                  (1-arg*sargsw)
                  /(1-model_.bulkJctSideGradingCoeff));
      capbd=Cbd*sarg+
        Cbdsw*sargsw;
    }
    else
    {
      qbd = f4d +
        vbd * (f2d + vbd * f3d/2);
      capbd=f2d + vbd * f3d;
    }
  }
  else
  {
    qbd = 0;
    capbd = 0;
  }

  // Now  after a mess of convergence stuff that seems not to apply to us
  // 3f5 saves the vbs, vbd, etc. in the state vector (we don't, it gets
  // saved in updatePrimaryState)

  // Then 3f5 calculates meyer capacitances, capgs, capgb and capgd
  // Careful!  They use the local von and vdsat, which haven't got dtype
  // multiplying them!

  ////
  // *     calculate meyer's capacitors
  ////
  ////
  // * new cmeyer - this just evaluates at the current time,
  // * expects you to remember values from previous time
  // * returns 1/2 of non-constant portion of capacitance
  // * you must add in the other half from previous time
  // * and the constant part
  ////

  if (mode > 0)
  {
    devSupport.qmeyer (vgs,vgd,Vgb,Von,Vdsat,
                       capgs, capgd, capgb, tPhi,OxideCap);
  }
  else
  {
    devSupport.qmeyer (vgd,vgs,Vgb,Von,Vdsat,
                       capgd, capgs, capgb, tPhi,OxideCap);
  }


  //  vgs1 = vgs_old;
  //  vgd1 = vgs1 - vds_old;
  //   vgb1 = vgs1 - vbs_old;

  ////
  // TVR:  Note!  Caution with Capgs vs. capgs.  If I used the instance var
  // capgs on the left hand side, it's OK as long as we never try to implement
  // meyer back averaging, since capgs is going to be recalculated each
  // time through.  But if we do the averaging, the old one will have the
  // constant part and old old part added in, which is not what we wanted.
  // So I'll continue 3f5's way of having a local Capgs that's actually used
  // for the charge computation, and an instance capgs that's saved as state.
  ////

  if((getSolverState().dcopFlag))
//  if(!(getSolverState().tranopFlag))
  {
    Capgs =  2 * capgs + GateSourceOverlapCap ;
    Capgd =  2 * capgd + GateDrainOverlapCap ;
    Capgb =  2 * capgb + GateBulkOverlapCap ;
  }
  else
  {
    capgs_old = (*extData.currStaVectorPtr)[li_state_capgs];
    capgd_old = (*extData.currStaVectorPtr)[li_state_capgd];
    capgb_old = (*extData.currStaVectorPtr)[li_state_capgb];

    Capgs = ( capgs+
              capgs_old +
              GateSourceOverlapCap );
    Capgd = ( capgd+
              capgd_old +
              GateDrainOverlapCap );
    Capgb = ( capgb+
              capgb_old +
              GateBulkOverlapCap );
  }

  // Eric does this kludge in level 1, so I'll do it to make sure capacitances
  // are positive
  Capgs *= (Capgs < 0.0)?-1.0:1.0;
  Capgd *= (Capgd < 0.0)?-1.0:1.0;
  Capgb *= (Capgb < 0.0)?-1.0:1.0;


  // in the sequel, I'll refer to all derivatives of currents w.r.t. voltages
  // as conductances, whether they really are or not.

  // when we get here, cdrain has the channel current
  // cbd and cbs have the currents through the diodes
  // cd has the drain current minus the diode current and seems only to be used
  // for 3f5 convergence stuff
  // qbs and qbd have the charges on the parasitic capacitors
  // capbs and capbd have the capacitances of the parasitics.

  // none of the charges for the gate capacitors have been calculated yet.
  // We've saved the capacitances, so we can get the charges in
  // updatePrimaryState later.

  // Conductances:
  // gbd: the bulk-drain' conductance without the capacitor components
  //      We'll need to get the capacitor contribution in the actual load
  //      using C*dt
  // gbs: bulk-source' without capacitor
  // gm = derivative of channel current w.r.t. gate-source voltage --- gotta
  //     account for mode=normal or mode=reverse when using this!
  // gmbs = derivative of channel current w.r.t bulk-source voltage
  // gds = derivative of channel current w.r.t. drain-source voltage

  // the variables gcgs, gcgb, gcgd should be the conductances for the gate
  // capacitors, but we won't do those here (vide supra), we'll do them
  // in the jacobian load given the capacitances.

  // Now 3f5 doesn't do the resistor currents in the RHS load, because of
  // how they do their numerics.  We do, so let's save those here.

  Idrain = drainConductance * Vddp;
  Isource = sourceConductance * Vssp;

  if (mode >= 0)   // Normal mode
  {
    Gm = gm;      // (xnrm-xrev)*gm  in 3f5
    Gmbs = gmbs;  // (xnrm-xrev)*gmbs in 3f5
    nrmsum = Gm+Gmbs; // xnrm*(gm+gmbs)
    revsum = 0;       // xrev*(gm+gmbs)
    cdreq = model_.dtype*cdrain;
  }
  else
  {
    Gm = -gm;
    Gmbs = -gmbs;
    nrmsum = 0;
    revsum = -(Gm+Gmbs);  // because Gm and Gmbs already have - in them!
    cdreq = -(model_.dtype)*cdrain;
  }

  // It is now essential to restore the vds/vgs/vgd variables that might
  //   have been tweaked by homotopy, lest they have an effect on RHS
  // Jdxp terms.

  vds=vds_save;
  vgs=vgs_save;
  vgd=vgd_save;

  /// CURRENTS to load into RHS:

  // so at this point:

  // current out of drain is
  // Idrain

  // current out of gate:
  // dtype*( (deriv of qgs) + (deriv of qgd) + (deriv of qgb))

  //  the current *out of* the source should be simply
  // Isource.

  // current out of bulk is
  // dtype*(deriv of qbd) + dtype*cbd + dtype*cbs + dtype*(deriv of qbs)
  //  - dtype*(deriv of qgb)

  // current out of drain' is
  // -Idrain - dtype*(deriv of qgd) - (deriv of qbd) - dtype*cbd +
  //  mode*dtype*cdrain

  // the current out of the source' is
  //  -Isource - dtype*(deriv of qgs) - dtype*cbs - (deriv of qbs) -
  //   mode*dtype*cdrain

  //////Conductances to load into Jacobian as they relate to things here:
  /// all of the places where I say Cap/dt I really mean dtype*Cap/dt for
  // the meyer capacitors.  No dtype for the parasitics, though

  // 3f5 handles the mode by doing:
  //        where xnrm=1, xrev=0 if mode>=0, xnrm=0,xrev=1 if mode<0

  // drain-drain = a = drainConductance
  // drain-drain' = b = -drainConductance

  // gate-gate = c = "gcgd+gcgs+gcgb" = Capgd/dt+Capgs/dt+Capgb/dt
  // gate-bulk = d = -gcgb = -Capgb/dt
  // gate-drain' = e = -gcgd = -Capgd/dt
  // gate-source' = f = -gcgs = -Capgs/dt

  // source-source = g = sourceConductance
  // source-source' = h = -sourceConductance

  // bulk-gate = i = -gcgb = -Capgb/dt
  // bulk-bulk = j = gbs+gbd+gcgb+parasitics=gbs+gbd+Capgb/dt+capbs/dt+capbd/dt
  // bulk-drain' = k= -gbd-capbd/dt
  // bulk-source' = l= -gbs-capbs/dt

  // drain'-drain = m = -drainConductance
  // drain'-gate = n = (xnrm-xrev)*gm-Capgd/dt
  // drain'-bulk = o = -gbd-capbd/dt+(xnrm-xrev)*gmbs
  // drain'-drain' = p = drainConductance+gds+gbd+capbd/dt+
  //                   xrev*(gm+gmbs)+ Capgd/dt
  // drain'-source' = q = -gds-xnrm*(gm+gmbs)

  // source'-gate = r = -(xnrm-xrev)*gm-Capgs/dt
  // source'-source = s = -sourceConductance
  // source'-bulk = t = -gbs-capbs/dt-(xnrm-xrev)*gmbs
  // source'-drain' = u= -gds-xrev*(gm+gmbs)
  // source'-source' = v = sourceConductance+gds+gbs+capbs/dt+xnrm*(gm+gmbs)+
  //                       Capgs/dt

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 02/27/01
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp_tmp)
{
  // mos3temp vars
  double czbd;    // zero voltage bulk-drain capacitance
  double czbdsw;  // zero voltage bulk-drain sidewall capacitance
  double czbs;    // zero voltage bulk-source capacitance
  double czbssw;  // zero voltage bulk-source sidewall capacitance
  double arg;     // 1 - fc
  double sarg;    // (1-fc) ^^ (-mj)
  double sargsw;  // (1-fc) ^^ (-mjsw)
  double ratio,ratio4;
  double fact2;
  double kt;
  double egfet;
  double pbfact;
  double capfact;
  double phio;
  double pbo;
  double gmanew,gmaold;
  // end of mos3temp stuff

  double tnom;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "  Instance::Begin of updateTemperature. \n";
    Xyce::dout() <<" name = " << getName() << std::endl;
    Xyce::dout() << std::endl;
  }

  // first set the instance temperature to the new temperature:
  if (temp_tmp != -999.0) temp = temp_tmp;

  if (model_.interpolateTNOM(temp))
  {
    model_.processParams();
  }

  tnom = model_.tnom;
  ratio = temp/tnom;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Temperature = "<< temp << std::endl;
    Xyce::dout() << "tnom = " << tnom << std::endl;
    Xyce::dout() << "ratio = " << ratio << std::endl;
  }

  vt = temp * CONSTKoverQ;
  ratio = temp/tnom;
  fact2 = temp/CONSTREFTEMP;
  kt = temp * CONSTboltz;
  egfet = 1.16-(7.02e-4*temp*temp)/(temp+1108);
  arg = -egfet/(kt+kt)+1.1150877/(CONSTboltz*(CONSTREFTEMP+CONSTREFTEMP));
  pbfact = -2*vt *(1.5*log(fact2)+CONSTQ*arg);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "vt = " << vt << std::endl;
    Xyce::dout() << "ratio = " << ratio << std::endl;
    Xyce::dout() << "fact2 = " << fact2 << std::endl;
    Xyce::dout() << "kt = " << kt << std::endl;
    Xyce::dout() << "egfet = " << egfet << std::endl;
    Xyce::dout() << "arg = " << arg << std::endl;
    Xyce::dout() << "pbfact = " << pbfact << std::endl;
  }

  // in mos3temp 3f5 does a bunch of parameter defaulting (lines 155-163)
  //  but we assume that our various parameters have been set already in
  //  the constructors

  // lines 164-203 of mos3temp moved to instance block constructor

  // Here's the entire guts of the mos3temp instance loop, with obvious
  // modifications (here->MOS3 goes away, model->MOS3 turns into model_.)

  ratio4 = ratio * sqrt(ratio);
  tTransconductance = model_.transconductance / ratio4;
  tSurfMob = model_.surfaceMobility/ratio4;
  phio= (model_.phi-model_.pbfact1)/model_.fact1;
  tPhi = fact2 * phio + pbfact;
  tVbi = model_.vt0 - model_.dtype *
    (model_.gamma* sqrt(model_.phi))+.5*(model_.egfet1-egfet)
    + model_.dtype*.5* (tPhi-model_.phi);
  tVto = tVbi + model_.dtype * model_.gamma * sqrt(tPhi);
  tSatCur = model_.jctSatCur* exp(-egfet/vt+model_.egfet1/model_.vtnom);
  tSatCurDens = model_.jctSatCurDensity * exp(-egfet/vt+model_.egfet1/model_.vtnom);
  pbo = (model_.bulkJctPotential - model_.pbfact1)/model_.fact1;
  gmaold = (model_.bulkJctPotential-pbo)/pbo;
  capfact = 1/(1+model_.bulkJctBotGradingCoeff*
               (4e-4*(model_.tnom-CONSTREFTEMP)-gmaold));
  tCbd = model_.capBD * capfact;
  tCbs = model_.capBS * capfact;
  tCj = model_.bulkCapFactor * capfact;
  capfact = 1/(1+model_.bulkJctSideGradingCoeff*
               (4e-4*(model_.tnom-CONSTREFTEMP)-gmaold));
  tCjsw = model_.sideWallCapFactor * capfact;
  tBulkPot = fact2 * pbo+pbfact;
  gmanew = (tBulkPot-pbo)/pbo;
  capfact = (1+model_.bulkJctBotGradingCoeff*(4e-4*(temp-CONSTREFTEMP)-gmanew));
  tCbd *= capfact;
  tCbs *= capfact;
  tCj *= capfact;
  capfact = (1+model_.bulkJctSideGradingCoeff*(4e-4*(temp-CONSTREFTEMP)-gmanew));
  tCjsw *= capfact;
  tDepCap = model_.fwdCapDepCoeff * tBulkPot;

  if( (model_.jctSatCurDensity == 0) || (drainArea == 0) ||
      (sourceArea == 0) )
  {
    sourceVcrit = drainVcrit =
      vt*log(vt/(CONSTroot2*model_.jctSatCur));
  }
  else
  {
    drainVcrit = vt * log( vt / (CONSTroot2 *
                                 model_.jctSatCurDensity * drainArea));
    sourceVcrit = vt * log( vt / (CONSTroot2 *
                                  model_.jctSatCurDensity * sourceArea));
  }
  if(model_.capBDGiven)
  {
    czbd = tCbd;
  }
  else
  {
    if(model_.bulkCapFactorGiven)
    {
      czbd=tCj*drainArea;
    }
    else
    {
      czbd=0;
    }
  }
  if(model_.sideWallCapFactorGiven)
  {
    czbdsw= tCjsw * drainPerimeter;
  }
  else
  {
    czbdsw=0;
  }
  arg = 1-model_.fwdCapDepCoeff;
  sarg = exp( (-model_.bulkJctBotGradingCoeff) * log(arg) );
  sargsw = exp( (-model_.bulkJctSideGradingCoeff) * log(arg) );
  Cbd = czbd;
  Cbdsw = czbdsw;

  // The following lines were taken verbatim from the SPICE level 3 mosfet.
  // There is a mistake in f3d and f4d, which incorrectly use bulkJctPotential
  // instead of tBulkPot.  The result of this mistake is to have an incorrect
  // discontinuity in the charge calculations --- the derivative of charge
  // is right as the voltage drop approaches the discontinuity from either side
  // but when numerical differentiation is done to get currents one obtains a
  // massive contribution to the right hand side.  Since spice never checks
  // the RHS residual norm this seems not to have any effect on SPICE, but
  // Xyce pukes.
  //  f2d = czbd*(1-model_.fwdCapDepCoeff*
  //              (1+model_.bulkJctBotGradingCoeff))* sarg/arg
  //    +  czbdsw*(1-model_.fwdCapDepCoeff*
  //               (1+model_.bulkJctSideGradingCoeff))*
  //    sargsw/arg;
  //  f3d = czbd * model_.bulkJctBotGradingCoeff * sarg/arg/
  //    model_.bulkJctPotential
  //    + czbdsw * model_.bulkJctSideGradingCoeff * sargsw/arg /
  //    model_.bulkJctPotential;
  //  f4d = czbd*model_.bulkJctPotential*(1-arg*sarg)/
  //    (1-model_.bulkJctBotGradingCoeff)
  //    + czbdsw*model_.bulkJctPotential*(1-arg*sargsw)/
  //    (1-model_.bulkJctSideGradingCoeff)
  //    -f3d/2*
  //    (tDepCap*tDepCap)
  //    -tDepCap * f2d;

  // These lines were taken from the equivalent section of the level 1 mosfet
  // and do not have the mistake in the lines above.
  f2d = czbd*(1-model_.fwdCapDepCoeff*
              (1+model_.bulkJctBotGradingCoeff))* sarg/arg
    +  czbdsw*(1-model_.fwdCapDepCoeff*
               (1+model_.bulkJctSideGradingCoeff))*
    sargsw/arg;
  f3d = czbd * model_.bulkJctBotGradingCoeff * sarg/arg/
    tBulkPot
    + czbdsw * model_.bulkJctSideGradingCoeff * sargsw/arg /
    tBulkPot;
  f4d = czbd*tBulkPot*(1-arg*sarg)/
    (1-model_.bulkJctBotGradingCoeff)
    + czbdsw*tBulkPot*(1-arg*sargsw)/
    (1-model_.bulkJctSideGradingCoeff)
    -f3d/2*
    (tDepCap*tDepCap)
    -tDepCap * f2d;
  if(model_.capBSGiven)
  {
    czbs=tCbs;
  }
  else
  {
    if(model_.bulkCapFactorGiven)
    {
      czbs=tCj*sourceArea;
    }
    else
    {
      czbs=0;
    }
  }
  if(model_.sideWallCapFactorGiven)
  {
    czbssw = tCjsw * sourcePerimeter;
  }
  else
  {
    czbssw=0;
  }
  arg = 1-model_.fwdCapDepCoeff;
  sarg = exp( (-model_.bulkJctBotGradingCoeff) * log(arg) );
  sargsw = exp( (-model_.bulkJctSideGradingCoeff) * log(arg) );
  Cbs = czbs;
  Cbssw = czbssw;

  // See the comments above regarding f3d and f4d --- the same error
  // exists in the SPICE mosfet3 code
  //  f2s = czbs*(1-model_.fwdCapDepCoeff*
  //              (1+model_.bulkJctBotGradingCoeff))* sarg/arg
  //    +  czbssw*(1-model_.fwdCapDepCoeff*
  //               (1+model_.bulkJctSideGradingCoeff))*
  //    sargsw/arg;
  //  f3s = czbs * model_.bulkJctBotGradingCoeff * sarg/arg/
  //    model_.bulkJctPotential
  //    + czbssw * model_.bulkJctSideGradingCoeff * sargsw/arg /
  //    model_.bulkJctPotential;
  //  f4s = czbs*model_.bulkJctPotential*(1-arg*sarg)/
  //    (1-model_.bulkJctBotGradingCoeff)
  //    + czbssw*model_.bulkJctPotential*(1-arg*sargsw)/
  //    (1-model_.bulkJctSideGradingCoeff)
  //    -f3s/2*
  //    (tBulkPot*tBulkPot)
  //    -tBulkPot * f2s;

  // so we use these lines from the MOSFET 1 instead.
  f2s = czbs*(1-model_.fwdCapDepCoeff*
              (1+model_.bulkJctBotGradingCoeff))* sarg/arg
    +  czbssw*(1-model_.fwdCapDepCoeff*
               (1+model_.bulkJctSideGradingCoeff))*
    sargsw/arg;
  f3s = czbs * model_.bulkJctBotGradingCoeff * sarg/arg/
    tBulkPot
    + czbssw * model_.bulkJctSideGradingCoeff * sargsw/arg /
    tBulkPot;
  f4s = czbs*tBulkPot*(1-arg*sarg)/
    (1-model_.bulkJctBotGradingCoeff)
    + czbssw*tBulkPot*(1-arg*sargsw)/
    (1-model_.bulkJctSideGradingCoeff)
    -f3s/2*
    (tDepCap*tDepCap)
    -tDepCap * f2s;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 02/28/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  double * staVector = extData.nextStaVectorRawPtr;
  double * oldstaVector = extData.currStaVectorRawPtr;
  double * stoVector = extData.nextStoVectorRawPtr;
  double * oldstoVector = extData.currStoVectorRawPtr;

  double vgs1, vgd1, vbs1,vgb1, vds1;

  bool bsuccess = updateIntermediateVars ();

  // voltage drops:
  stoVector[li_store_vbd] = vbd;
  stoVector[li_store_vbs] = vbs;
  stoVector[li_store_vgs] = vgs;
  stoVector[li_store_vds] = vds;
  stoVector[li_store_von] = von;

  // save the transconductance
  stoVector[li_store_gm] = Gm;

  // now the meyer capacitances
  // we didn't calculate these charges in update IntermediateVars
  // but we did calculate the voltage drops and capacitances.
  // first store the capacitances themselves:
  staVector[li_state_capgs] = capgs;
  staVector[li_state_capgd] = capgd;
  staVector[li_state_capgb] = capgb;

  // now the charges
  // BE CAREFUL!  We can only do Q=CV for DCOP!  Otherwise it's
  // supposed to be *INTEGRATED*:
  // Q = int(t0,t1)C(V)*dV --- and we approximate that by
  // Q(t1)-Q(t0) = CBar*(V(t1)-V(t0)) where CBar is the average.
  // Now with Meyer back averaging, Capxx is the average between the last
  // time step and this one.  So we gotta do the right thing for non-DCOP
  // when backaverage is on.

  if((getSolverState().dcopFlag))
//  if(!(getSolverState().tranopFlag))
  {
    qgs = Capgs*vgs;
    qgd = Capgd*vgd;
    qgb = Capgb*Vgb;
  }
  else
  {
    // get the ones from last time step
    qgs = oldstaVector[li_state_qgs];
    qgd = oldstaVector[li_state_qgd];
    qgb = oldstaVector[li_state_qgb];
    // get the voltage drops, too
    vgs1 = oldstoVector[li_store_vgs];
    vbs1 = oldstoVector[li_store_vbs];
    vds1 = oldstoVector[li_store_vds];

    vgb1 = vgs1-vbs1;
    vgd1 = vgs1-vds1;

    // NOW we can calculate the charge update
    qgs += Capgs*(vgs-vgs1);
    qgd += Capgd*(vgd-vgd1);
    qgb += Capgb*((vgs-vbs)-vgb1);
  }

  staVector[li_state_qgs] = qgs;
  staVector[li_state_qgd] = qgd;
  staVector[li_state_qgb] = qgb;

  // and the diode parasitic capacitors
  // these charges were set in updateIntermediateVars
  staVector[li_state_qbd] = qbd;
  staVector[li_state_qbs] = qbs;

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MOSFET3::Instance::getNumNoiseSources
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/27/2014
//-----------------------------------------------------------------------------
int Instance::getNumNoiseSources () const
{
  return 4;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MOSFET3::Instance::setupNoiseSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void Instance::setupNoiseSources (Xyce::Analysis::NoiseData & noiseData)
{
  int numSources=4;
  noiseData.numSources = numSources;
  noiseData.resize(numSources);

  noiseData.deviceName = getName().getEncodedName();

  // Note: the letter suffixes (e.g., rd) are used by the DNO() and DNI() 
  // operators for .PRINT NOISE
  noiseData.noiseNames[0] = "noise_" + getName().getEncodedName()+
    std::string("_rd");  // noise due to rd 
  noiseData.noiseNames[1] = "noise_" + getName().getEncodedName()+
    std::string("_rs");  // noise due to rs 
  noiseData.noiseNames[2] = "noise_" + getName().getEncodedName()+
    std::string("_id");  // noise due to id 
  noiseData.noiseNames[3] = "noise_" + getName().getEncodedName()+
    std::string("_fn"); // flicker (1/f) noise

  // RD thermal:
  noiseData.li_Pos[0] = li_DrainPrime;
  noiseData.li_Neg[0] = li_Drain;

  // RS thermal:
  noiseData.li_Pos[1] = li_SourcePrime;
  noiseData.li_Neg[1] = li_Source;

  // ID thermal:
  noiseData.li_Pos[2] = li_DrainPrime;
  noiseData.li_Neg[2] = li_SourcePrime;

  // flicker:
  noiseData.li_Pos[3] = li_DrainPrime;
  noiseData.li_Neg[3] = li_SourcePrime;

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MOSFET3::Instance::getNoiseSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void Instance::getNoiseSources (Xyce::Analysis::NoiseData & noiseData)
{
  // thermal noise, RD:
  devSupport.noiseSupport(
      noiseData.noiseDens[0], noiseData.lnNoiseDens[0], THERMNOISE, 
                  drainConductance*numberParallel,
                  temp);

  // thermal noise, RS:
  devSupport.noiseSupport(
      noiseData.noiseDens[1], noiseData.lnNoiseDens[1], THERMNOISE, 
                  sourceConductance*numberParallel,
                  temp);

  // thermal noise, ID:
  devSupport.noiseSupport(
      noiseData.noiseDens[2], noiseData.lnNoiseDens[2], THERMNOISE, 
                  (2.0/3.0 * fabs(gm))*numberParallel, temp);

  // flicker noise 
  noiseData.noiseDens[3] = numberParallel*model_.fNcoef * std::exp(model_.fNexp *
				 std::log(std::max(fabs(cd),N_MINLOG))) /
				 (noiseData.freq * w * (l - 2*model_.latDiff) *
				 model_.oxideCapFactor * model_.oxideCapFactor);

  noiseData.lnNoiseDens[3] = std::log(std::max(noiseData.noiseDens[3],N_MINLOG));
}

// Additional Declarations

// Class Model

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
  double wkfngs;
  double wkfng;
  double fermig;
  double fermis;
  double vfb;
  double kt1;
  double arg1;

  fact1 = tnom/CONSTREFTEMP;
  vtnom = tnom*CONSTKoverQ;
  kt1 = CONSTboltz * tnom;
  egfet1 = 1.16-(7.02e-4*tnom*tnom)/(tnom+1108);
  arg1 = -egfet1/(kt1+kt1)+1.1150877/(CONSTboltz*(CONSTREFTEMP+CONSTREFTEMP));
  pbfact1 = -2*vtnom *(1.5*log(fact1)+CONSTQ*arg1);

  if(oxideThickness == 0)
  {
    UserError(*this) << getName() << " has TOX=0";
  }
  else
  {
    oxideCapFactor = 3.9 * 8.854214871e-12/oxideThickness;
  }
  if(!given("U0") && !given("UO")) surfaceMobility=600;
  if(!given("KP"))
    transconductance = surfaceMobility * oxideCapFactor * 1e-4;
  if(given("NSUB"))
  {
    if(substrateDoping*1e6 >1.45e16)
    {
      if(!given("PHI"))
      {
        phi = 2*vtnom*
          log(substrateDoping*1e6/1.45e16);
        phi = std::max(0.1,phi);
      }
      fermis = dtype * .5 * phi;
      wkfng = 3.2;
      if(!given("TPG")) gateType=1;
      if(gateType != 0)
      {
        fermig = dtype *gateType*.5*egfet1;
        wkfng = 3.25 + .5 * egfet1 - fermig;
      }
      wkfngs = wkfng - (3.25 + .5 * egfet1 +fermis);
      if(!given("GAMMA"))
      {
        gamma = sqrt(2 * 11.70 * CONSTperm0 *
                     CONSTQ * substrateDoping*1e6)/
          oxideCapFactor;
      }
      if(!given("VTO") && !given("VT0"))
      {
        if(!given("NSS"))
          surfaceStateDensity=0;
        vfb = wkfngs - surfaceStateDensity*1e4*CONSTQ/oxideCapFactor;
        vt0 = vfb + dtype * (gamma * sqrt(phi)+ phi);
      }
      else
      {
        vfb = vt0 - dtype * (gamma*sqrt(phi)+phi);
      }
      alpha = ((11.7*CONSTperm0)+(11.7*CONSTperm0))/
        (CONSTQ*substrateDoping*1e6); //(cm**3/m**3)
      coeffDepLayWidth = sqrt(alpha);
    }
    else
    {
      UserError(*this) << "Nsub < Ni";
    }
  }

  // now model parameter preprocessing
  double mpi = M_PI;
  narrowFactor = delta * 0.5 * mpi * (11.7*CONSTperm0) /oxideCapFactor ;

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
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 2/26/01
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    dtype(CONSTNMOS),
    tnom(getDeviceOptions().tnom),
    latDiff(0.0),
    jctSatCurDensity(0.0),
    jctSatCur(0.0),
    drainResistance(0.0),
    sourceResistance(0.0),
    sheetResistance(0.0),
    transconductance(0.0),
    gateSourceOverlapCapFactor(0.0),
    gateDrainOverlapCapFactor(0.0),
    gateBulkOverlapCapFactor(0.0),
    oxideCapFactor(0.0),
    vt0(0.0),
    capBD(0.0),
    capBS(0.0),
    bulkCapFactor(0.0),
    sideWallCapFactor(0.0),
    bulkJctPotential(0.0),
    bulkJctBotGradingCoeff(0.0),
    bulkJctSideGradingCoeff(0.0),
    fwdCapDepCoeff(0.0),
    phi(0.0),
    gamma(0.0),
    substrateDoping(0.0),
    gateType(0),
    surfaceStateDensity(0.0),
    oxideThickness(0.0),
    surfaceMobility(0.0),
    eta(0.0),
    junctionDepth(0.0),
    coeffDepLayWidth(0.0),
    narrowFactor(0.0),
    delta(0.0),
    fastSurfaceStateDensity(0.0),
    theta(0.0),
    maxDriftVel(0.0),
    alpha(0.0),
    kappa(0.0),
    fNcoef(0.0),
    fNexp(0.0),
    capBDGiven(0),
    capBSGiven(0),
    bulkCapFactorGiven(0),
    sideWallCapFactorGiven(0)
{
  if (getType() != "")
  {
    if (getType() == "NMOS") {
      dtype = CONSTNMOS;
    }
    else if (getType() == "PMOS") {
      dtype = CONSTPMOS;
    }
    else
    {
      UserError(*this) << "Could not recognize the type for model " << getName();
    }
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
  if (!given("L"))
    model_l=getDeviceOptions().defl;
  if (!given("W"))
    model_w=getDeviceOptions().defw;
  if (!given("TNOM"))
    tnom = getDeviceOptions().tnom;
  if (capBD != 0)
    capBDGiven = true;
  if (capBS != 0)
    capBSGiven = true;

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
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
    os << "  " << i << ": " << (*iter)->getName() << "\t";
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


//-----------------------------------------------------------------------------
// MOSFET3 Master functions:
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
    Instance & mi = *(*it);

    double * oldstaVectorPtr = mi.extData.currStaVectorRawPtr;
    double * stoVec = mi.extData.nextStoVectorRawPtr;
    double * oldstoVec = mi.extData.currStoVectorRawPtr;

    double vgs1(0.0), vgd1(0.0), vbs1(0.0),vgb1(0.0), vds1(0.0);

    bool btmp = mi.updateIntermediateVars ();
    bsuccess = bsuccess && btmp;

    // voltage drops:
    stoVec[mi.li_store_vbd] = mi.vbd;
    stoVec[mi.li_store_vbs] = mi.vbs;
    stoVec[mi.li_store_vgs] = mi.vgs;
    stoVec[mi.li_store_vds] = mi.vds;
    stoVec[mi.li_store_von] = mi.von;

    // save the transconductance
    stoVec[mi.li_store_gm]  = mi.Gm;

    // now the meyer capacitances
    // we didn't calculate these charges in update IntermediateVars
    // but we did calculate the voltage drops and capacitances.
    // first store the capacitances themselves:
    staVec[mi.li_state_capgs] = mi.capgs;
    staVec[mi.li_state_capgd] = mi.capgd;
    staVec[mi.li_state_capgb] = mi.capgb;

    // now the charges
    // BE CAREFUL!  We can only do Q=CV for DCOP!  Otherwise it's
    // supposed to be *INTEGRATED*:
    // Q = int(t0,t1)C(V)*dV --- and we approximate that by
    // Q(t1)-Q(t0) = CBar*(V(t1)-V(t0)) where CBar is the average.
    // Now with Meyer back averaging, Capxx is the average between the last
    // time step and this one.  So we gotta do the right thing for non-DCOP
    // when backaverage is on.

    if((getSolverState().dcopFlag))
//    if(!(getSolverState().tranopFlag))
    {
      mi.qgs = mi.Capgs*mi.vgs;
      mi.qgd = mi.Capgd*mi.vgd;
      mi.qgb = mi.Capgb*mi.Vgb;
    }
    else
    {
      // get the ones from last time step
      mi.qgs = (oldstaVectorPtr)[mi.li_state_qgs];
      mi.qgd = (oldstaVectorPtr)[mi.li_state_qgd];
      mi.qgb = (oldstaVectorPtr)[mi.li_state_qgb];
      // get the voltage drops, too
      vgs1 = oldstoVec[mi.li_store_vgs];
      vbs1 = oldstoVec[mi.li_store_vbs];
      vds1 = oldstoVec[mi.li_store_vds];

      vgb1 = vgs1-vbs1;
      vgd1 = vgs1-vds1;

      // NOW we can calculate the charge update
      mi.qgs += mi.Capgs*(mi.vgs-vgs1);
      mi.qgd += mi.Capgd*(mi.vgd-vgd1);
      mi.qgb += mi.Capgb*((mi.vgs-mi.vbs)-vgb1);
    }

    staVec[mi.li_state_qgs] = mi.qgs;
    staVec[mi.li_state_qgd] = mi.qgd;
    staVec[mi.li_state_qgb] = mi.qgb;

    // and the diode parasitic capacitors
    // these charges were set in updateIntermediateVars
    staVec[mi.li_state_qbd] = mi.qbd;
    staVec[mi.li_state_qbs] = mi.qbs;
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
  double gmin1 = getDeviceOptions().gmin;

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & mi = *(*it);

    int Dtype=mi.getModel().dtype;
    double ceqbs(0.0),ceqbd(0.0),ceqgb(0.0), ceqgs(0.0), ceqgd(0.0);
    double Qeqbs(0.0),Qeqbd(0.0),Qeqgb(0.0), Qeqgs(0.0), Qeqgd(0.0);
    double coef(0.0);

    // F-Vector:
    ceqbs = Dtype*(mi.cbs);
    ceqbd = Dtype*(mi.cbd);
    // These need "Dtype" here because we use them later *without*
    // Dtype, where SPICE uses it *with*
    ceqgb = 0.0;
    ceqgs = 0.0;
    ceqgd = 0.0;

    if (mi.drainConductance != 0.0)
    {
      fVec[mi.li_Drain] += mi.Idrain*mi.numberParallel;
    }

    coef = (ceqgs+ceqgd+ceqgb);

    fVec[mi.li_Gate] += coef*mi.numberParallel;

    if (mi.sourceConductance != 0.0)
    {
      fVec[mi.li_Source] += mi.Isource*mi.numberParallel;
    }

    coef = ceqbs + ceqbd - ceqgb;

    fVec[mi.li_Bulk] += coef*mi.numberParallel;

    coef = -mi.Idrain-(ceqbd - mi.cdreq + ceqgd);

    fVec[mi.li_DrainPrime] += coef*mi.numberParallel;

    coef = -mi.Isource-(ceqbs + mi.cdreq + ceqgs);

    fVec[mi.li_SourcePrime] += coef*mi.numberParallel;

    // Q-Vector:
    Qeqbs = Dtype*(mi.qbs);
    Qeqbd = Dtype*(mi.qbd);
    // These need "Dtype" here because we use them later *without*
    // Dtype, where SPICE uses it *with*
    Qeqgb = Dtype*(mi.qgb);
    Qeqgs = Dtype*(mi.qgs);
    Qeqgd = Dtype*(mi.qgd);

    coef = (Qeqgs+Qeqgd+Qeqgb);

    qVec[mi.li_Gate] += coef*mi.numberParallel;

    coef = Qeqbs + Qeqbd - Qeqgb;

    qVec[mi.li_Bulk] += coef*mi.numberParallel;

    coef = -(Qeqbd + Qeqgd);

    qVec[mi.li_DrainPrime] += coef*mi.numberParallel;

    coef = -(Qeqbs + Qeqgs);

    qVec[mi.li_SourcePrime] += coef*mi.numberParallel;

    // voltage limiters:
    if (!mi.origFlag)
    {
      // F-limiters:
      double coef_Jdxp4 = Dtype*(
            + ((mi.gbd-gmin1))*(mi.vbd-mi.vbd_orig)
            + ((mi.gbs-gmin1))*(mi.vbs-mi.vbs_orig));

      double coef_Jdxp5 = Dtype*(
            -((mi.gbd-gmin1))*(mi.vbd-mi.vbd_orig)
            +mi.gds*(mi.vds-mi.vds_orig)
            +mi.Gm*((mi.mode>0)?(mi.vgs-mi.vgs_orig):(mi.vgd-mi.vgd_orig))
            +mi.Gmbs*((mi.mode>0)?(mi.vbs-mi.vbs_orig):(mi.vbd-mi.vbd_orig)));

      double coef_Jdxp6 = Dtype*(
            -((mi.gbs-gmin1))*(mi.vbs-mi.vbs_orig)
            -mi.gds*(mi.vds-mi.vds_orig)
            -mi.Gm*((mi.mode>0)?(mi.vgs-mi.vgs_orig):(mi.vgd-mi.vgd_orig))
            -mi.Gmbs*((mi.mode>0)?(mi.vbs-mi.vbs_orig):(mi.vbd-mi.vbd_orig)));

      double * dFdxdVp = mi.extData.dFdxdVpVectorRawPtr;
      dFdxdVp[mi.li_Bulk       ] += coef_Jdxp4*mi.numberParallel;
      dFdxdVp[mi.li_DrainPrime ] += coef_Jdxp5*mi.numberParallel;
      dFdxdVp[mi.li_SourcePrime] += coef_Jdxp6*mi.numberParallel;

      // Q-limiters:
      {
        double gcgd(0.0), gcgs(0.0), gcgb(0.0), gcbs(0.0), gcbd(0.0);
        if (getSolverState().tranopFlag || getSolverState().acopFlag || getSolverState().transientFlag)
        {
          gcgd = mi.Capgd;
          gcgs = mi.Capgs;
          gcgb = mi.Capgb;
          // get at the two parasitic caps the same way
          gcbs = mi.capbs;
          gcbd = mi.capbd;
        }
        else
        {
          gcgd = 0.0; gcgs = 0.0; gcgb = 0.0; gcbs = 0.0; gcbd = 0.0;
        }

        double coef_Jdxp2 =
          Dtype*(gcgd*(mi.vgd-mi.vgd_orig)+gcgs*(mi.vgs-mi.vgs_orig)+
          gcgb*(mi.vgs-mi.vgs_orig-mi.vbs+mi.vbs_orig));

        double coef_Jdxp4 = Dtype*(
            - (gcgb)*(mi.vgs-mi.vgs_orig-mi.vbs+mi.vbs_orig)
            + (gcgb)*(mi.vbd-mi.vbd_orig)
            + (gcbs)*(mi.vbs-mi.vbs_orig));

        double coef_Jdxp5 = Dtype*(
            -(gcgd)*(mi.vgd-mi.vgd_orig)
            -(gcbd)*(mi.vbd-mi.vbd_orig));

        // 6 KCL for source' node
        double coef_Jdxp6 = Dtype*
          (-gcgs*(mi.vgs-mi.vgs_orig)-(gcbs)*(mi.vbs-mi.vbs_orig));

        double * dQdxdVp = mi.extData.dQdxdVpVectorRawPtr;
        dQdxdVp[mi.li_Gate       ] += coef_Jdxp2*mi.numberParallel;
        dQdxdVp[mi.li_Bulk       ] += coef_Jdxp4*mi.numberParallel;
        dQdxdVp[mi.li_DrainPrime ] += coef_Jdxp5*mi.numberParallel;
        dQdxdVp[mi.li_SourcePrime] += coef_Jdxp6*mi.numberParallel;
      }
    }

    if( mi.loadLeadCurrent )
    {
      leadF[mi.li_branch_dev_id] = (-(ceqbd - mi.cdreq + ceqgd))*mi.numberParallel;
      leadQ[mi.li_branch_dev_id] = (-(Qeqbd + Qeqgd))*mi.numberParallel;
      leadF[mi.li_branch_dev_is] = (-(ceqbs + mi.cdreq + ceqgs))*mi.numberParallel;
      leadQ[mi.li_branch_dev_is] = (-(Qeqbs + Qeqgs))*mi.numberParallel;
      leadF[mi.li_branch_dev_ig] = (ceqgs+ceqgd+ceqgb)*mi.numberParallel;
      leadQ[mi.li_branch_dev_ig] = (Qeqgs+Qeqgd+Qeqgb)*mi.numberParallel;
      leadF[mi.li_branch_dev_ib] = (ceqbs + ceqbd - ceqgb)*mi.numberParallel;
      leadQ[mi.li_branch_dev_ib] = (Qeqbs + Qeqbd - Qeqgb)*mi.numberParallel;

      junctionV[mi.li_branch_dev_id] = solVec[mi.li_Drain] - solVec[mi.li_Source];
      junctionV[mi.li_branch_dev_ig] = solVec[mi.li_Gate] - solVec[mi.li_Source];
      junctionV[mi.li_branch_dev_is] = 0.0;
      junctionV[mi.li_branch_dev_ib] = 0.0 ; 
    }
  }

  return true;
}

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
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
    Instance & mi = *(*it);

    // F-matrix:

    *mi.f_DrainEquDrainNodePtr +=
    mi.drainConductance*mi.numberParallel;

    *mi.f_DrainEquDrainPrimeNodePtr -=
    mi.drainConductance*mi.numberParallel;


    *mi.f_SourceEquSourceNodePtr +=
    mi.sourceConductance*mi.numberParallel;

    *mi.f_SourceEquSourcePrimeNodePtr -=
    mi.sourceConductance*mi.numberParallel;


    *mi.f_BulkEquBulkNodePtr +=
    (mi.gbs+mi.gbd)*mi.numberParallel;

    *mi.f_BulkEquDrainPrimeNodePtr -= mi.gbd*mi.numberParallel;

    *mi.f_BulkEquSourcePrimeNodePtr -= mi.gbs*mi.numberParallel;


    *mi.f_DrainPrimeEquDrainNodePtr -=
    mi.drainConductance*mi.numberParallel;

    *mi.f_DrainPrimeEquGateNodePtr +=
    mi.Gm*mi.numberParallel;

    *mi.f_DrainPrimeEquBulkNodePtr +=
    (-mi.gbd+mi.Gmbs)*mi.numberParallel;

    *mi.f_DrainPrimeEquDrainPrimeNodePtr +=
    (mi.drainConductance+mi.gds+mi.gbd+mi.revsum)*mi.numberParallel;

    *mi.f_DrainPrimeEquSourcePrimeNodePtr +=
    (-mi.gds-mi.nrmsum)*mi.numberParallel;


    *mi.f_SourcePrimeEquGateNodePtr -=
    mi.Gm*mi.numberParallel;

    *mi.f_SourcePrimeEquSourceNodePtr -=
    mi.sourceConductance*mi.numberParallel;

    *mi.f_SourcePrimeEquBulkNodePtr -=
    (mi.gbs+mi.Gmbs)*mi.numberParallel;

    *mi.f_SourcePrimeEquDrainPrimeNodePtr -=
    (mi.gds+mi.revsum)*mi.numberParallel;

    *mi.f_SourcePrimeEquSourcePrimeNodePtr +=
    (mi.sourceConductance+mi.gds+mi.gbs+mi.nrmsum)*mi.numberParallel;

    // Q-matrix:
    double gcgd(0.0);  // d(cqgd)/dVgd
    double gcgs(0.0);  // d(cqgs)/dVgs
    double gcgb(0.0);  // d(cqgb)/dVgb
    double gcbs(0.0);  // d(cqbs)/dVbs
    double gcbd(0.0);  // d(cqbd)/dVbd

    // get at the "conductances" for the gate capacitors with this trick
    //      gcgd = model_.dtype*Capgd;
    //      gcgs = model_.dtype*Capgs;
    //      gcgb = model_.dtype*Capgb;
    //
    //      In the loadRHS function, these would all be multiplied by
    //      getSolverState().pdt.  Here, for *mi.q_, the pdt term is left out.
    if (getSolverState().tranopFlag || getSolverState().acopFlag || getSolverState().transientFlag)
    {
      gcgd = mi.Capgd;
      gcgs = mi.Capgs;
      gcgb = mi.Capgb;
      // get at the two parasitic caps the same way
      gcbs = mi.capbs;
      gcbd = mi.capbd;
    }


    *mi.q_GateEquGateNodePtr +=
    (gcgd+gcgs+gcgb)*mi.numberParallel;

    *mi.q_GateEquBulkNodePtr -= gcgb*mi.numberParallel;

    *mi.q_GateEquDrainPrimeNodePtr -= gcgd*mi.numberParallel;

    *mi.q_GateEquSourcePrimeNodePtr -= gcgs*mi.numberParallel;


    *mi.q_BulkEquGateNodePtr -= gcgb*mi.numberParallel;

    *mi.q_BulkEquBulkNodePtr +=
    (+gcbs+gcbd+gcgb)*mi.numberParallel;

    *mi.q_BulkEquDrainPrimeNodePtr -= +gcbd*mi.numberParallel;

    *mi.q_BulkEquSourcePrimeNodePtr -=
    +gcbs*mi.numberParallel;


    *mi.q_DrainPrimeEquGateNodePtr +=
    -gcgd*mi.numberParallel;

    *mi.q_DrainPrimeEquBulkNodePtr +=
    -gcbd*mi.numberParallel;

    *mi.q_DrainPrimeEquDrainPrimeNodePtr +=
    (+gcbd+gcgd)*mi.numberParallel;


    *mi.q_SourcePrimeEquGateNodePtr -=
    gcgs*mi.numberParallel;

    *mi.q_SourcePrimeEquBulkNodePtr -=
    +gcbs*mi.numberParallel;

    *mi.q_SourcePrimeEquSourcePrimeNodePtr+=
    (+gcbs+gcgs)*mi.numberParallel;
  }
  return true;
}
#else
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
    Instance & mi = *(*it);

    // F-matrix:
    dFdx[mi.li_Drain][mi.ADrainEquDrainNodeOffset] +=
    mi.drainConductance*mi.numberParallel;

    dFdx[mi.li_Drain][mi.ADrainEquDrainPrimeNodeOffset] -=
    mi.drainConductance*mi.numberParallel;

    dFdx[mi.li_Source][mi.ASourceEquSourceNodeOffset] +=
    mi.sourceConductance*mi.numberParallel;

    dFdx[mi.li_Source][mi.ASourceEquSourcePrimeNodeOffset] -=
    mi.sourceConductance*mi.numberParallel;

    dFdx[mi.li_Bulk][mi.ABulkEquBulkNodeOffset] +=
    (mi.gbs+mi.gbd)*mi.numberParallel;

    dFdx[mi.li_Bulk][mi.ABulkEquDrainPrimeNodeOffset] -= mi.gbd*mi.numberParallel;
    dFdx[mi.li_Bulk][mi.ABulkEquSourcePrimeNodeOffset] -= mi.gbs*mi.numberParallel;


    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquDrainNodeOffset] -=
    mi.drainConductance*mi.numberParallel;

    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquGateNodeOffset] +=
    mi.Gm*mi.numberParallel;

    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquBulkNodeOffset] +=
    (-mi.gbd+mi.Gmbs)*mi.numberParallel;

    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquDrainPrimeNodeOffset] +=
    (mi.drainConductance+mi.gds+mi.gbd+mi.revsum)*mi.numberParallel;

    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquSourcePrimeNodeOffset] +=
    (-mi.gds-mi.nrmsum)*mi.numberParallel;


    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquGateNodeOffset] -=
    mi.Gm*mi.numberParallel;

    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquSourceNodeOffset] -=
    mi.sourceConductance*mi.numberParallel;

    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquBulkNodeOffset] -=
    (mi.gbs+mi.Gmbs)*mi.numberParallel;

    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquDrainPrimeNodeOffset] -=
    (mi.gds+mi.revsum)*mi.numberParallel;

    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquSourcePrimeNodeOffset] +=
    (mi.sourceConductance+mi.gds+mi.gbs+mi.nrmsum)*mi.numberParallel;

    // Q-matrix:
    double gcgd(0.0);  // d(cqgd)/dVgd
    double gcgs(0.0);  // d(cqgs)/dVgs
    double gcgb(0.0);  // d(cqgb)/dVgb
    double gcbs(0.0);  // d(cqbs)/dVbs
    double gcbd(0.0);  // d(cqbd)/dVbd

    // get at the "conductances" for the gate capacitors with this trick
    //      gcgd = model_.dtype*Capgd;
    //      gcgs = model_.dtype*Capgs;
    //      gcgb = model_.dtype*Capgb;
    //
    //      In the loadRHS function, these would all be multiplied by
    //      getSolverState().pdt.  Here, for dQdx, the pdt term is left out.
    if (getSolverState().tranopFlag || getSolverState().acopFlag || getSolverState().transientFlag)
    {
      gcgd = mi.Capgd;
      gcgs = mi.Capgs;
      gcgb = mi.Capgb;
      // get at the two parasitic caps the same way
      gcbs = mi.capbs;
      gcbd = mi.capbd;
    }


    dQdx[mi.li_Gate][mi.AGateEquGateNodeOffset] +=
    (gcgd+gcgs+gcgb)*mi.numberParallel;

    dQdx[mi.li_Gate][mi.AGateEquBulkNodeOffset] -= gcgb*mi.numberParallel;
    dQdx[mi.li_Gate][mi.AGateEquDrainPrimeNodeOffset] -= gcgd*mi.numberParallel;
    dQdx[mi.li_Gate][mi.AGateEquSourcePrimeNodeOffset] -= gcgs*mi.numberParallel;
    dQdx[mi.li_Bulk][mi.ABulkEquGateNodeOffset] -= gcgb*mi.numberParallel;
    dQdx[mi.li_Bulk][mi.ABulkEquBulkNodeOffset] +=
    (+gcbs+gcbd+gcgb)*mi.numberParallel;
    dQdx[mi.li_Bulk][mi.ABulkEquDrainPrimeNodeOffset] -= +gcbd*mi.numberParallel;
    dQdx[mi.li_Bulk][mi.ABulkEquSourcePrimeNodeOffset] -=
    +gcbs*mi.numberParallel;
    dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquGateNodeOffset] +=
    -gcgd*mi.numberParallel;
    dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquBulkNodeOffset] +=
    -gcbd*mi.numberParallel;
    dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquDrainPrimeNodeOffset] +=
    (+gcbd+gcgd)*mi.numberParallel;
    dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquGateNodeOffset] -=
    gcgs*mi.numberParallel;
    dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquBulkNodeOffset] -=
    +gcbs*mi.numberParallel;
    dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquSourcePrimeNodeOffset]+=
    (+gcbs+gcgs)*mi.numberParallel;
  }
  return true;
}
#endif

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() ||
      ((deviceMap.find("M")!=deviceMap.end()) && (levelSet.find(3)!=levelSet.end())))
  {
    MOSFET1::registerDevice();

    Config<Traits>::addConfiguration()
      .registerDevice("m", 3)
      .registerModelType("pmos", 3)
      .registerModelType("nmos", 3);
  }
}

} // namespace MOSFET3
} // namespace Device
} // namespace Xyce

//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Purpose        : Implement the Vertical double-diffused power MOSFET
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
#include <climits>

// ---------- Xyce Includes ----------
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MOSFET1.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Message.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_VDMOS.h>
#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_MachDepParams.h>

namespace Xyce {
namespace Device {


namespace VDMOS {


void Traits::loadInstanceParameters(ParametricData<VDMOS::Instance> &p)
{
// Set up map for normal (double) param variables:
    p.addPar ("L",0.0,&VDMOS::Instance::l)
     .setUnit(U_METER)
     .setCategory(CAT_GEOMETRY)
     .setDescription("Channel length")
     .setLengthScaling(true);

    p.addPar ("W",0.0,&VDMOS::Instance::w)
     .setUnit(U_METER)
     .setCategory(CAT_GEOMETRY)
     .setDescription("Channel width")
     .setLengthScaling(true);

    p.addPar ("AD",0.0,&VDMOS::Instance::drainArea)
     .setUnit(U_METER2)
     .setCategory(CAT_GEOMETRY)
     .setDescription("Drain diffusion area")
     .setAreaScaling(true);

    p.addPar ("AS",0.0,&VDMOS::Instance::sourceArea)
     .setUnit(U_METER2)
     .setCategory(CAT_GEOMETRY)
     .setDescription("Source diffusion area")
     .setAreaScaling(true);

    p.addPar ("NRD",1.0,&VDMOS::Instance::drainSquares)
     .setUnit(U_SQUARES)
     .setCategory(CAT_GEOMETRY)
     .setDescription("Multiplier for RSH to yield parasitic resistance of drain");

    p.addPar ("NRS",1.0,&VDMOS::Instance::sourceSquares)
     .setUnit(U_SQUARES)
     .setCategory(CAT_GEOMETRY)
     .setDescription("Multiplier for RSH to yield parasitic resistance of source");

    p.addPar ("PD",0.0,&VDMOS::Instance::drainPerimeter)
     .setUnit(U_METER)
     .setCategory(CAT_GEOMETRY)
     .setDescription("Drain diffusion perimeter")
     .setLengthScaling(true);

    p.addPar ("PS",0.0,&VDMOS::Instance::sourcePerimeter)
     .setUnit(U_METER)
     .setCategory(CAT_GEOMETRY)
     .setDescription("Source diffusion perimeter")
     .setLengthScaling(true);

    p.addPar ("M",1.0,&VDMOS::Instance::numberParallel)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Multiplier for M devices connected in parallel");

    p.addPar ("TEMP",0.0,&VDMOS::Instance::temp)
     .setExpressionAccess(ParameterType::TIME_DEP)
     .setUnit(STANDARD)
     .setCategory(CAT_NONE)
     .setDescription("Device temperature");
}

void Traits::loadModelParameters(ParametricData<VDMOS::Model> &p)
{
    p.addPar ("L0",0.0,&VDMOS::Model::l0)
     .setUnit(U_METER)
     .setCategory(CAT_GEOMETRY)
     .setDescription("Gate length of nominal device");

    p.addPar ("W0",0.0,&VDMOS::Model::w0)
     .setUnit(U_METER)
     .setCategory(CAT_GEOMETRY)
     .setDescription("Gate width of nominal device");

    p.addPar ("VTO",0.0,&VDMOS::Model::vt0)
     .setUnit(U_VOLT)
     .setCategory(CAT_VOLT)
     .setDescription("Zero-bias threshold voltage");

    p.addPar ("PHI",0.6,&VDMOS::Model::phi)
     .setUnit(U_VOLT)
     .setCategory(CAT_PROCESS)
     .setDescription("Surface potential");

    p.addPar ("RD",0.0,&VDMOS::Model::drainResistance)
        .setExpressionAccess(ParameterType::MIN_RES)
      .setUnit(U_OHM)
     .setCategory(CAT_RES)
     .setDescription("Drain ohmic resistance");

    p.addPar ("RG",0.0,&VDMOS::Model::gateResistance)
     .setUnit(U_OHM)
     .setCategory(CAT_RES)
     .setDescription("Gate ohmic resistance");

    p.addPar ("RS",0.0,&VDMOS::Model::sourceResistance)
     .setExpressionAccess(ParameterType::MIN_RES)
     .setUnit(U_OHM)
     .setCategory(CAT_RES)
     .setDescription("Source ohmic resistance");

    p.addPar ("CBD",0.0,&VDMOS::Model::capBD)
     .setGivenMember(&VDMOS::Model::capBDGiven)
     .setExpressionAccess(ParameterType::MIN_CAP)
     .setUnit(U_FARAD)
     .setCategory(CAT_CAP)
     .setDescription("Zero-bias bulk-drain p-n capacitance");

    p.addPar ("CBS",0.0,&VDMOS::Model::capBS)
     .setGivenMember(&VDMOS::Model::capBSGiven)
     .setExpressionAccess(ParameterType::MIN_CAP)
     .setUnit(U_FARAD)
     .setCategory(CAT_CAP)
     .setDescription("Zero-bias bulk-source p-n capacitance");

    p.addPar ("TS",0.0,&VDMOS::Model::timeScale)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("IS",1e-14,&VDMOS::Model::jctSatCur)
     .setUnit(U_AMP)
     .setCategory(CAT_CURRENT)
     .setDescription("Bulk p-n saturation current");

    p.addPar ("PB",0.8,&VDMOS::Model::bulkJctPotential)
     .setUnit(U_VOLT)
     .setCategory(CAT_VOLT)
     .setDescription("Bulk p-n bottom potential");

    p.addPar ("CGSO",0.0,&VDMOS::Model::gateSourceOverlapCapFactor)
     .setUnit(U_FARADMM1)
     .setCategory(CAT_CAP)
     .setDescription("Gate-source overlap capacitance/channel width");

    p.addPar ("CGDO",0.0,&VDMOS::Model::gateDrainOverlapCapFactor)
     .setUnit(U_FARADMM1)
     .setCategory(CAT_CAP)
     .setDescription("Gate-drain overlap capacitance/channel width");

    p.addPar ("CGBO",0.0,&VDMOS::Model::gateBulkOverlapCapFactor)
     .setUnit(U_FARADMM1)
     .setCategory(CAT_CAP)
     .setDescription("Gate-bulk overlap capacitance/channel length");

    p.addPar ("RSH",0.0,&VDMOS::Model::sheetResistance)
     .setUnit(U_OHM)
     .setCategory(CAT_RES)
     .setDescription("Drain,source diffusion sheet resistance");

    p.addPar ("CJ",0.0,&VDMOS::Model::bulkCapFactor)
     .setGivenMember(&VDMOS::Model::bulkCapFactorGiven)
     .setUnit(U_FARADMM2)
     .setCategory(CAT_CAP)
     .setDescription("Bulk p-n zero-bias bottom capacitance/area");

    p.addPar ("MJ",0.5,&VDMOS::Model::bulkJctBotGradingCoeff)
     .setUnit(U_NONE)
     .setCategory(CAT_DOPING)
     .setDescription("Bulk p-n bottom grading coefficient");

    p.addPar ("CJSW",0.0,&VDMOS::Model::sideWallCapFactor)
     .setGivenMember(&VDMOS::Model::sideWallCapFactorGiven)
     .setUnit(U_FARADMM2)
     .setCategory(CAT_CAP)
     .setDescription("Bulk p-n zero-bias sidewall capacitance/area");

    p.addPar ("MJSW",0.5,&VDMOS::Model::bulkJctSideGradingCoeff)
     .setUnit(U_NONE)
     .setCategory(CAT_DOPING)
     .setDescription("Bulk p-n sidewall grading coefficient");

    p.addPar ("JS",0.0,&VDMOS::Model::jctSatCurDensity)
     .setUnit(U_AMPMM2)
     .setCategory(CAT_PROCESS)
     .setDescription("Bulk p-n saturation current density");

    p.addPar ("TOX",1e-7,&VDMOS::Model::oxideThickness)
     .setUnit(U_METER)
     .setCategory(CAT_GEOMETRY)
     .setDescription("Gate oxide thickness");

    p.addPar ("LD",0.0,&VDMOS::Model::latDiff)
      .setUnit(U_METER)
     .setCategory(CAT_DOPING)
     .setDescription("Lateral diffusion length");

    p.addPar ("UO",280., &VDMOS::Model::surfaceMobility)
     .setUnit(U_CMM2VM1SM1)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("Surface mobility");

    p.addPar ("U0",280., &VDMOS::Model::surfaceMobility0)
     .setUnit(U_CMM2VM1SM1)
     .setCategory(CAT_PROCESS)
     .setDescription("Surface mobility");

    p.addPar ("FC",0.5,&VDMOS::Model::fwdCapDepCoeff)
     .setUnit(U_NONE)
     .setCategory(CAT_CAP)
     .setDescription("Coefficient for forward-bias depletion capacitance formula");

    p.addPar ("NSUB",0.0,&VDMOS::Model::substrateDoping)
     .setUnit(U_CMM3)
     .setCategory(CAT_DOPING)
     .setDescription("Substrate doping density");

    p.addPar ("NSS",0.0,&VDMOS::Model::surfaceStateDensity)
      .setUnit(U_CMM2)
     .setCategory(CAT_PROCESS)
     .setDescription("Surface state density");

    p.addPar ("TNOM",0.0,&VDMOS::Model::tnom)
     .setUnit(STANDARD)
     .setCategory(CAT_NONE)
     .setDescription("Parameter measurement temperature");

    p.addPar ("VMAX",4e+4,&VDMOS::Model::maxDriftVel)
     .setUnit(U_MSM1)
     .setCategory(CAT_PROCESS)
     .setDescription("Maximum drift velocity for carriers");

    p.addPar ("XJ",0.0,&VDMOS::Model::junctionDepth)
     .setUnit(U_METER)
     .setCategory(CAT_GEOMETRY)
     .setDescription("Metallurgical junction depth");

    p.addPar ("LAMBDA",0.048,&VDMOS::Model::lambda)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_PROCESS)
     .setDescription("Output conductance parameter");

    p.addPar ("DELTA",5.0,&VDMOS::Model::delta)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Transition width parameter");

    p.addPar ("ETA",1.32,&VDMOS::Model::eta)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Subthreshold ideality factor");

    p.addPar ("M",4.0,&VDMOS::Model::m)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Knee shape parameter");

    p.addPar ("MC",3.0,&VDMOS::Model::mc)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("SIGMA0",0.048,&VDMOS::Model::sigma0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("DIBL parameter");

    p.addPar ("VSIGMAT",1.7,&VDMOS::Model::vsigmat)
     .setUnit(U_VOLT)
     .setCategory(CAT_VOLT)
     .setDescription("DIBL parameter");

    p.addPar ("VSIGMA",0.2,&VDMOS::Model::vsigma)
     .setUnit(U_VOLT)
     .setCategory(CAT_VOLT)
     .setDescription("DIBL parameter");

    p.addPar ("THETA",0.0,&VDMOS::Model::theta)
     .setUnit(U_MVM1)
     .setCategory(CAT_PROCESS)
     .setDescription("Mobility degradation parameter");

    p.addPar ("GAMMAS0",0.5,&VDMOS::Model::gammas0)
     .setUnit(U_VOLTMH)
     .setCategory(CAT_NONE)
     .setDescription("Body effect constant in front of square root term");

    p.addPar ("GAMMAL0",0.0,&VDMOS::Model::gammal0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body effect constant in front of linear term");

    p.addPar ("LGAMMAS",0.0,&VDMOS::Model::lgammas)
     .setUnit(U_VOLTMH)
     .setCategory(CAT_NONE)
     .setDescription("Sensitivity of gS on device length");

    p.addPar ("WGAMMAS",0.0,&VDMOS::Model::wgammas)
     .setUnit(U_VOLTMH)
     .setCategory(CAT_NONE)
     .setDescription("Sensitivity of gS on device width");

    p.addPar ("LGAMMAL",0.0,&VDMOS::Model::lgammal_)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Sensitivity of gL on device length");

    p.addPar ("WGAMMAL",0.0,&VDMOS::Model::wgammal)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Sensitivity of gL on device width");

    p.addPar ("K",0.0,&VDMOS::Model::k)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("KVT",0.0,&VDMOS::Model::kvt)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("MDTEMP",0.0,&VDMOS::Model::mdtemp)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("KVS",0.0,&VDMOS::Model::kvs)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("TVS",0.0,&VDMOS::Model::tvs)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("MTH",0.0,&VDMOS::Model::mth)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("ARTD",0.0,&VDMOS::Model::artd)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("BRTD",0.035,&VDMOS::Model::brtd)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("CRTD",0.1472,&VDMOS::Model::crtd)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("DRTD",0.0052,&VDMOS::Model::drtd)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("NRTD",0.115,&VDMOS::Model::nrtd)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("N2",1.0,&VDMOS::Model::n2)
      .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("XQC",0.6,&VDMOS::Model::xqc)
     .setUnit(U_NONE)
     .setCategory(CAT_PROCESS)
     .setDescription("Charge partitioning factor");

    p.addPar ("MCV",10.0,&VDMOS::Model::mcv)
     .setUnit(U_NONE)
     .setCategory(CAT_GEOMETRY)
     .setDescription("Transition width parameter used by the charge partitioning scheme");

    p.addPar ("VFB",0.0,&VDMOS::Model::vfb)
     .setUnit(U_VOLT)
     .setCategory(CAT_VOLT)
     .setDescription("Flat band voltage");

    p.addPar ("ALPHA",0.0,&VDMOS::Model::alpha)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Parameter accounting for the threshold dependence on the channel potential");

    p.addPar ("LS",35e-9,&VDMOS::Model::ls)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("RSUB",0.0,&VDMOS::Model::rsub)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("VP",0.0,&VDMOS::Model::vp)
     .setGivenMember(&VDMOS::Model::vpGiven)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("AI",2e+9,&VDMOS::Model::ai)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("BI",8e+8,&VDMOS::Model::bi)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("DELMAX",0.9,&VDMOS::Model::delmax)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("MD",2.0,&VDMOS::Model::md)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    p.addPar ("DRIFTPARAMA",0.08,&VDMOS::Model::driftParamA)
     .setUnit(U_OHM)
     .setCategory(CAT_RES)
     .setDescription("Drift region resistance intercept parameter");

    p.addPar ("DRIFTPARAMB",0.013,&VDMOS::Model::driftParamB)
     .setUnit(U_OHMPV)
     .setCategory(CAT_RES)
     .setDescription("Drift region resistance slope parameter");

    p.addPar ("RDSSHUNT",0.0,&VDMOS::Model::rdsshunt)
     .setUnit(U_OHM)
     .setCategory(CAT_RES)
     .setDescription("Drain-source shunt resistance");

    p.addPar ("D1IS",1e-14,&VDMOS::Model::D1DIOsatCur)
     .setUnit(U_AMP)
     .setCategory(CAT_CURRENT)
     .setDescription("Drain-Source diode saturation current");

    p.addPar ("D1RS",0.0,&VDMOS::Model::D1DIOresist)
     .setUnit(U_OHM)
     .setCategory(CAT_RES)
     .setDescription("Drain-source diode ohmic resistance");

    p.addPar ("D1N",1.0,&VDMOS::Model::D1DIOemissionCoeff)
     .setUnit(U_NONE)
     .setCategory(CAT_PROCESS)
     .setDescription("Drain-source diode emission coefficient");

    p.addPar ("D1TT",0.0,&VDMOS::Model::D1DIOtransitTime)
     .setUnit(U_SECOND)
     .setCategory(CAT_PROCESS)
     .setDescription("Drain-source diode transit time");

    p.addPar ("D1CJO",0.0,&VDMOS::Model::D1DIOjunctionCap)
     .setUnit(U_FARAD)
     .setCategory(CAT_CAP)
     .setDescription("Drain-source diode junction capacitance");

    p.addPar ("D1VJ",1.0,&VDMOS::Model::D1DIOjunctionPot)
     .setUnit(U_VOLT)
     .setCategory(CAT_VOLT)
     .setDescription("Drain-source diode junction potential");

    p.addPar ("D1M",0.5,&VDMOS::Model::D1DIOgradingCoeff)
     .setUnit(U_NONE)
     .setCategory(CAT_PROCESS)
     .setDescription("Drain-source diode grading coefficient");

    p.addPar ("D1EG",1.11,&VDMOS::Model::D1DIOactivationEnergy)
     .setUnit(U_EV)
     .setCategory(CAT_PROCESS)
     .setDescription("Drain-source diode activation energy");

    p.addPar ("D1XTI",3.0,&VDMOS::Model::D1DIOsaturationCurrentExp)
     .setUnit(U_NONE)
     .setCategory(CAT_PROCESS)
     .setDescription("Drain-source diode sat. current temperature exponent");

    p.addPar ("D1FC",0.5,&VDMOS::Model::D1DIOdepletionCapCoeff)
     .setUnit(U_NONE)
     .setCategory(CAT_CAP)
     .setDescription("Drain-source diode forward bias depletion capacitance");

    p.addPar ("D1BV",1E99,&VDMOS::Model::D1DIObreakdownVoltage)
     .setGivenMember(&VDMOS::Model::D1DIObreakdownVoltageGiven)
     .setUnit(U_VOLT)
     .setCategory(CAT_VOLT)
     .setDescription("Drain-source diode reverse breakdown voltage");

    p.addPar ("D1IBV",1e-3,&VDMOS::Model::D1DIObreakdownCurrent)
     .setUnit(U_AMP)
     .setCategory(CAT_CURRENT)
     .setDescription("Drain-source diode current at breakdown voltage");

    p.addPar ("D1IKF",0.0,&VDMOS::Model::D1DIOikf)
     .setUnit(U_AMP)
     .setCategory(CAT_CURRENT)
     .setDescription("Drain-source diode high injection knee currrent");

    p.addPar ("D1ISR",0.0,&VDMOS::Model::D1DIOisr)
     .setUnit(U_AMP)
     .setCategory(CAT_CURRENT)
     .setDescription("Drain-source diode recombination saturation current");

    p.addPar ("D1NR",2.0,&VDMOS::Model::D1DIOnr)
     .setUnit(U_NONE)
     .setCategory(CAT_PROCESS)
     .setDescription("Drain-source diode recombination emission coefficient");

    p.addPar ("D1KF",0.0,&VDMOS::Model::D1DIOfNcoef)
     .setUnit(U_NONE)
     .setCategory(CAT_FLICKER)
     .setDescription("Drain-source diode flicker noise coefficient");

    p.addPar ("D1AF",1.0,&VDMOS::Model::D1DIOfNexp)
     .setUnit(U_NONE)
     .setCategory(CAT_FLICKER)
     .setDescription("Drain-source diode flicker noise exponent");

    p.addPar ("D1TNOM",300.15,&VDMOS::Model::D1DIOnomTemp)
     .setUnit(U_DEGC)
     .setCategory(CAT_TEMP)
     .setDescription("Drain-source diode nominal temperature");

    // Set up non-double precision variables:
    p.addPar ("TPG",1,&VDMOS::Model::gateType)
     .setUnit(U_NONE)
     .setCategory(CAT_MATERIAL)
     .setDescription("Gate material type (-1 = same as substrate, 0 = aluminum,1 = opposite of substrate)");

    p.addPar ("CV",1,&VDMOS::Model::cv)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Charge model storage selector");

    p.addPar ("CVE",1,&VDMOS::Model::cve)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Meyer-like capacitor model selector");

    p.addPar ("FPE",1,&VDMOS::Model::fpe)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Charge partitioning scheme selector");

    p.addPar ("ISUBMOD",0,&VDMOS::Model::isubmod)
     .setUnit(U_NONE)
     .setCategory(ParameterCategory(CAT_PROCESS | UNDOCUMENTED))
     .setDescription("");

    DeviceModel::initThermalModel(p);
}

std::vector< std::vector<int> > Instance::jacStamp_DC_SC_GC;
std::vector< std::vector<int> > Instance::jacStamp_DC_GC;
std::vector< std::vector<int> > Instance::jacStamp_SC_GC;
std::vector< std::vector<int> > Instance::jacStamp_DC_SC;
std::vector< std::vector<int> > Instance::jacStamp_GC;
std::vector< std::vector<int> > Instance::jacStamp_SC;
std::vector< std::vector<int> > Instance::jacStamp_DC;
std::vector< std::vector<int> > Instance::jacStamp;

std::vector<int> Instance::jacMap_DC_SC_GC;
std::vector<int> Instance::jacMap_DC_GC;
std::vector<int> Instance::jacMap_SC_GC;
std::vector<int> Instance::jacMap_DC_SC;
std::vector<int> Instance::jacMap_GC;
std::vector<int> Instance::jacMap_SC;
std::vector<int> Instance::jacMap_DC;
std::vector<int> Instance::jacMap;

std::vector< std::vector<int> > Instance::jacMap2_DC_SC_GC;
std::vector< std::vector<int> > Instance::jacMap2_DC_GC;
std::vector< std::vector<int> > Instance::jacMap2_SC_GC;
std::vector< std::vector<int> > Instance::jacMap2_DC_SC;
std::vector< std::vector<int> > Instance::jacMap2_GC;
std::vector< std::vector<int> > Instance::jacMap2_SC;
std::vector< std::vector<int> > Instance::jacMap2_DC;
std::vector< std::vector<int> > Instance::jacMap2;

// duplicating, but with RD1RS nonzero
std::vector< std::vector<int> > Instance::jacStamp_D1C_DC_SC_GC;
std::vector< std::vector<int> > Instance::jacStamp_D1C_DC_GC;
std::vector< std::vector<int> > Instance::jacStamp_D1C_SC_GC;
std::vector< std::vector<int> > Instance::jacStamp_D1C_DC_SC;
std::vector< std::vector<int> > Instance::jacStamp_D1C_GC;
std::vector< std::vector<int> > Instance::jacStamp_D1C_SC;
std::vector< std::vector<int> > Instance::jacStamp_D1C_DC;
std::vector< std::vector<int> > Instance::jacStamp_D1C;

std::vector<int> Instance::jacMap_D1C_DC_SC_GC;
std::vector<int> Instance::jacMap_D1C_DC_GC;
std::vector<int> Instance::jacMap_D1C_SC_GC;
std::vector<int> Instance::jacMap_D1C_DC_SC;
std::vector<int> Instance::jacMap_D1C_GC;
std::vector<int> Instance::jacMap_D1C_SC;
std::vector<int> Instance::jacMap_D1C_DC;
std::vector<int> Instance::jacMap_D1C;

std::vector< std::vector<int> > Instance::jacMap2_D1C_DC_SC_GC;
std::vector< std::vector<int> > Instance::jacMap2_D1C_DC_GC;
std::vector< std::vector<int> > Instance::jacMap2_D1C_SC_GC;
std::vector< std::vector<int> > Instance::jacMap2_D1C_DC_SC;
std::vector< std::vector<int> > Instance::jacMap2_D1C_GC;
std::vector< std::vector<int> > Instance::jacMap2_D1C_SC;
std::vector< std::vector<int> > Instance::jacMap2_D1C_DC;
std::vector< std::vector<int> > Instance::jacMap2_D1C;

//--------------------- Class Model ----------------------------

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
    gateType(0),

    l0(2e-6),
    w0(2e-5),
    tnom(getDeviceOptions().tnom),
    latDiff(0.0),
    jctSatCurDensity(0.0),
    jctSatCur(0.0),
    drainResistance(0.0),
    gateResistance(0.0),
    sourceResistance(0.0),
    sheetResistance(0.0),
    gateSourceOverlapCapFactor(0.0),
    gateDrainOverlapCapFactor(0.0),
    gateBulkOverlapCapFactor(0.0),
    oxideCapFactor(0.0),
    vt0(0.0),
    capBD(0.0),
    capBS(0.0),
    timeScale(0.0),
    bulkCapFactor(0.0),
    sideWallCapFactor(0.0),
    bulkJctPotential(0.0),
    bulkJctBotGradingCoeff(0.0),
    bulkJctSideGradingCoeff(0.0),
    fwdCapDepCoeff(0.0),
    phi(0.0),
    gamma(0.0),
    lambda(0.0),
    substrateDoping(0.0),
    surfaceStateDensity(0.0),
    oxideThickness(0.0),
    surfaceMobility(0.0),
    surfaceMobility0(0.0),

    capBDGiven(0),
    capBSGiven(0),
    bulkCapFactorGiven(0),
    sideWallCapFactorGiven(0),
    vpGiven(0),

    maxDriftVel(0.0),
    junctionDepth(0.0),
    rdi(0.0),
    rsi(0.0),

    delta(5.0),
    eta(1.32),
    m(4.0),
    mc(3.0),
    sigma0(0.048),
    vsigmat(1.7),
    vsigma(0.2),
    theta(0.0),
    gammas0(0.5),
    gammal0(0.0),
    lgammas(0.0),
    wgammas(0.0),
    lgammal_(0.0),
    wgammal(0.0),
    kacc(0.0),
    gb(0.0),
    knit(0.0),
    nitd(0.0),
    nits(0.0),
    mm(0.5),
    k(0.0),
    deltaSqr(0.0),
    kvt(0.0),
    mdtemp(0.0),
    kvs(0.0),
    tvs(0.0),
    mth(0.0),
    artd(0.0),
    brtd(0.035),
    crtd(0.1472),
    drtd(0.0052),
    nrtd(0.115),
    n2(1.0),

    xqc(0.6),
    mcv(10.0),
    vfb(0.0),
    fpe(1),
    alpha(1.05),

    cv(1),
    cve(1),
    ls(35e-9),
    rsub(0.0),
    vp(0.0),
    ai(2e9),
    bi(8e8),
    delmax(0.9),
    md(2.0),
    isubmod(0),

    kaccd(0.0),
    kaccs(0.0),
    invm(0.0),
    invmc(0.0),
    invmd(0.0),

    fact1(0.0),
    vtnom(0.0),
    egfet1(0.0),
    pbfact1(0.0),

    driftParamA(0.08),
    driftParamB(0.013),
    rdsshunt(0.0),

    D1DIOsatCur(0.0),
    D1DIOresist(0.0),
    D1DIOconductance(0.0),
    D1DIOemissionCoeff(0.0),
    D1DIOtransitTime(0.0),
    D1DIOjunctionCap(0.0),
    D1DIOjunctionPot(0.0),
    D1DIOgradingCoeff(0.0),
    D1DIOactivationEnergy(0.0),
    D1DIOsaturationCurrentExp(0.0),
    D1DIOdepletionCapCoeff(0.0),
    D1DIObreakdownVoltage(0.0),
    D1DIObreakdownCurrent(0.0),
    D1DIOf2(0.0),
    D1DIOf3(0.0),
    D1DIOnomTemp(0.0),
    D1DIOfNcoef(0.0),
    D1DIOfNexp(0.0),
    D1DIOikf(0.0),
    D1DIOisr(0.0),
    D1DIOnr(0.0),
    D1DIObreakdownVoltageGiven(0)

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
  if (!given("TNOM"))
    tnom = getDeviceOptions().tnom;
  if (!given("L0"))
    l0 = getDeviceOptions().defl;
  if (!given("W0"))
    w0 = getDeviceOptions().defw;
  if (!given("ALPHA"))
    alpha = (dtype == CONSTNMOS) ? 1.05 : 1.3;
  if (capBD != 0)
    capBDGiven = true;
  if (capBS != 0)
    capBSGiven = true;

  // Calculate any parameters specified as expressions:

  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  if (given("U0"))
  {
    if (given("UO"))
      UserError(*this) << "Both uo and u0 have been specified and, which is not allowed";
    else
      UserWarning(*this) << "Surface mobility has been specified as u0 instead of uo, uo is the preferred syntax";

    surfaceMobility = surfaceMobility0;
  }

  if(eta == 0)
  {
    UserError(*this) << "ETA cannot be zero for level 18";
  }

  if(driftParamA == 0 && driftParamB == 0)
  {
    UserError(*this) << "Both driftParamA and driftParamB cannot be zero";
  }

  if(cv != 1 && cv != 2)
  {
    UserError(*this) << "Model error: use cv=1 (Meyer's model) or cv=2 (Meyer-like Model)";
  }

  if(fpe != 1 && fpe != 2 && fpe != 3)
  {
    UserError(*this) << "Model error: use fpe=1 (and xqc), fpe=2 (and xqc,mcv) or fpe=3";
  }

  if(mcv <= 1)
  {
    UserError(*this) << "Model error: charge conservation requires mcv > 1";
  }

  if(xqc > 1 || xqc < 0.5)
  {
    UserError(*this) << "Model error: charge conservation requires 0.5 <= xqc <= 1.0";
  }

  // New version of VDMOSModel. We have to set lambda to zero
  if(isubmod) lambda = 0.0;

  // limit diode #1 parameters
  // limit grading coeff to max of 0.9
  if(D1DIOgradingCoeff > 0.9)
  {
    UserWarning(*this) << "grading coefficient too large, limited to 0.9";
    D1DIOgradingCoeff = 0.9;
  }

  // limit activation energy to min of 0.1
  if(D1DIOactivationEnergy < 0.1)
  {
    UserWarning(*this) << "activation energy too small, limited to 0.1";
    D1DIOactivationEnergy = 0.1;
  }

  // limit depletion cap coeff to max of 0.95
  if(D1DIOdepletionCapCoeff > 0.95)
  {
    UserWarning(*this) << "coefficient Fc too large, limited to 0.95";
    D1DIOdepletionCapCoeff = 0.95;
  }

  processParams ();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "   " << std::endl;
    Xyce::dout() << "l0 " << l0 << std::endl;
    Xyce::dout() << "w0 " << w0 << std::endl;
    Xyce::dout() << "vt0 = " << vt0 << std::endl;
    Xyce::dout() << "gamma = " << gamma << std::endl;
    Xyce::dout() << "lambda = " << lambda << std::endl;
    Xyce::dout() << "phi = " << phi << std::endl;
    Xyce::dout() << "drainResistance = " << drainResistance << std::endl;
    Xyce::dout() << "gateResistance  = " << gateResistance << std::endl;
    Xyce::dout() << "sourceResistance = " << sourceResistance << std::endl;
    Xyce::dout() << "capBD = " << capBD << std::endl;
    Xyce::dout() << "capBS = " << capBS << std::endl;
    Xyce::dout() << "timeScale = " << timeScale << std::endl;
    Xyce::dout() << "jctSatCur = " << jctSatCur << std::endl;
    Xyce::dout() << "bulkJctPotential = " << bulkJctPotential << std::endl;
    Xyce::dout() << "gateSourceOverlapCapFactor = " << gateSourceOverlapCapFactor << std::endl;
    Xyce::dout() << "gateDrainOverlapCapFactor = " << gateDrainOverlapCapFactor << std::endl;
    Xyce::dout() << "gateBulkOverlapCapFactor = " << gateBulkOverlapCapFactor << std::endl;
    Xyce::dout() << "sheetResistance = " << sheetResistance << std::endl;
    Xyce::dout() << "bulkCapFactor   = " << bulkCapFactor   << std::endl;
    Xyce::dout() << "bulkJctBotGradingCoeff   = " << bulkJctBotGradingCoeff   << std::endl;
    Xyce::dout() << "sideWallCapFactor   = " << sideWallCapFactor   << std::endl;
    Xyce::dout() << "bulkJctSideGradingCoeff   = " << bulkJctSideGradingCoeff   << std::endl;
    Xyce::dout() << "jctSatCurDensity   = " << jctSatCurDensity   << std::endl;
    Xyce::dout() << "oxideThickness   = " << oxideThickness   << std::endl;
    Xyce::dout() << "oxideCapFactor = " << oxideCapFactor << std::endl;
    Xyce::dout() << "latDiff   = " << latDiff   << std::endl;
    Xyce::dout() << "surfaceMobility   = " << surfaceMobility   << std::endl;
    Xyce::dout() << "fwdCapDepCoeff   = " << fwdCapDepCoeff   << std::endl;
    Xyce::dout() << "substrateDoping   = " << substrateDoping   << std::endl;
    Xyce::dout() << "gateType   = " << gateType   << std::endl;
    Xyce::dout() << "surfaceStateDensity   = " << surfaceStateDensity   << std::endl;
    Xyce::dout() << "tnom   = " << tnom   << std::endl;
    Xyce::dout() << "driftParamA = " << driftParamA << std::endl;
    Xyce::dout() << "driftParamB = " << driftParamB << std::endl;
    Xyce::dout() << "rdsshunt = " << rdsshunt << std::endl;
  }

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
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : pmc
// Creation Date : 2/17/2004
//-----------------------------------------------------------------------------
bool Model::processParams ()
{
  double wkfngs(0.0);
  double wkfng(0.0);
  double fermig(0.0);
  double fermis(0.0);
  double kt1(0.0);
  double arg1(0.0);

  fact1  = tnom/CONSTREFTEMP;
  vtnom  = tnom*CONSTKoverQ;
  kt1    = CONSTboltz * tnom;
  egfet1 = 1.16-(7.02e-4*tnom*tnom)/(tnom+1108);
  arg1   = -egfet1/(kt1+kt1)+1.1150877/(CONSTboltz*(CONSTREFTEMP+CONSTREFTEMP));
  pbfact1= -2*vtnom *(1.5*log(fact1)+CONSTQ*arg1);

// assuming silicon - make definition for epsilon of silicon
//#define EPSSIL (11.7 * 8.854214871e-12) = CONSTEPSSI (the Xyce value)

  oxideCapFactor = 3.9 * 8.854214871e-12/oxideThickness;
  if(given("NSUB"))
  {
    if(substrateDoping*1e6 >1.45e16)
    {
      if(!given("PHI"))
      {
        phi = 2*vtnom* log(substrateDoping*1e6/1.45e16);
        phi = std::max(.1,phi);
      }
      fermis = dtype*0.5*phi;
      wkfng = 3.2;
      if(gateType != 0)
      {
        fermig = dtype *gateType*0.5*egfet1;
        wkfng = 3.25 + 0.5 * egfet1 - fermig;
      }
      wkfngs = wkfng - (3.25 + 0.5 * egfet1 + fermis);
      gamma = sqrt(2 * 11.70 * 8.854214871e-12 * CONSTQ
                                      * substrateDoping*1e6)/oxideCapFactor;
      if(!given("GAMMAS0")) gammas0 = gamma;

      if(!given("VTO"))
      {
        if(!given("NSS")) surfaceStateDensity=0;
        if(!given("VFB"))
        {
          vfb = wkfngs - surfaceStateDensity*1e4*CONSTQ/oxideCapFactor;
        }
        vt0 = vfb + dtype * (gamma * sqrt(phi)+ phi);
      }
      else
      {
        if(!given("VFB")) vfb = vt0 - dtype * (gamma * sqrt(phi)+ phi);
      }
    }
    else
    {
      UserError(*this) << "SubstrateDoping is 0";

      substrateDoping = 0;
    }
  }

  if(!given("CJ"))
  {
    bulkCapFactor = sqrt(CONSTEPSSI*CONSTQ*substrateDoping*1e6
           /(2*bulkJctPotential));
  }
  if(!given("VFB"))
  {
    if(fabs(vfb) < 1e-12) 
    {
      vfb = 1e-12;
    }
    else 
    {
      //vfb = vfb;
    }
  }

  kaccd = kacc/(1+pow(k*nitd,mm));
  kaccs = kacc/(1+pow(k*nits,mm));
  deltaSqr = delta*delta;

  if(D1DIOresist == 0)
    D1DIOconductance = 0;
  else
    D1DIOconductance = 1/D1DIOresist;

  double D1xfc = log(1-D1DIOdepletionCapCoeff);
  D1DIOf2 = exp((1 + D1DIOgradingCoeff)*D1xfc);
  D1DIOf3 = 1 - D1DIOdepletionCapCoeff*(1 + D1DIOgradingCoeff);

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

//------------------ Class Instance ----------------------------


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
  const InstanceBlock &         IB,
  Model &                       Miter,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    dNode(0),
    gNode(0),
    sNode(0),
    bNode(0),
    dNodePrime(0),
    sNodePrime(0),
    dDriftNode(0),
    l(getDeviceOptions().defl),
    w(getDeviceOptions().defw),
    drainArea(getDeviceOptions().defad),
    sourceArea(getDeviceOptions().defas),
    drainSquares(1.0),
    sourceSquares(1.0),
    drainPerimeter(0.0),
    sourcePerimeter(0.0),
    sourceCond(0.0),
    gateCond(0.0),
    drainCond(0.0),
    draindriftCond(0.0),
    numberParallel(1),
    vt(0.0),
    temp(getDeviceOptions().temp.getImmutableValue<double>()),
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

    von(0.0),
    vdsat(0.0),
    vddsat(0.0),
    sourceVcrit(0.0),
    drainVcrit(0.0),
    draindriftVcrit(0.0),
    cdd(0.0),
    cd(0.0),
    gmbs(0.0),
    gm(0.0),
    gddd(0.0),
    dIdd_dVd(0.0),
    gds(0.0),
    gdds(0.0),
    gdsshunt(0.0),
    gbs(0.0),
    gbd(0.0),
    cbd(0.0),
    Cbd(0.0),
    Cbdsw(0.0),
    cbs(0.0),
    Cbs(0.0),
    Cbssw(0.0),
    f2d(0.0),
    f3d(0.0),
    f4d(0.0),
    f2s(0.0),
    f3s(0.0),
    f4s(0.0),
    gammas(0.0),
    gammal(0.0),
    gchi0(0.0),
    vtoo(0.0),
    vthLimit(0.0),

    mode(1),
    mode_low(0.0),
    mode_high(0.0),
    off(0),
    dNodePrimeSet(0),
    sNodePrimeSet(0),
    limitedFlag(false),
    //calculated quantities
    EffectiveLength(0),
    DrainSatCur(0),
    SourceSatCur(0),
    GateSourceOverlapCap(0),
    GateDrainOverlapCap(0),
    GateBulkOverlapCap(0),
    OxideCap(0),

    Vd    (0.0),
    Vs    (0.0),
    Vg    (0.0),
    Vb    (0.0),
    Vdp   (0.0),
    Vgp   (0.0),
    Vsp   (0.0),
    Vdd   (0.0),
    Vddp  (0.0),
    Vddd  (0.0),
    Vdddp (0.0),
    Vssp  (0.0),
    Vbsp  (0.0),
    Vbdp  (0.0),
    Vggp  (0.0),
    Vgpsp (0.0),
    Vgpdp (0.0),
    Vgpb  (0.0),
    Vdpsp (0.0),
    Vbdd  (0.0),
    D1vd  (0.0),
    Capgs (0.0),
    Capgdd(0.0),
    Capgb (0.0),
    Isource(0.0),
    Igate(0.0),
    Idrain(0.0),
    Idraindrift(0.0),
    Irdsshunt(0.0),
    Ird1rs(0.0),
    mm1(0.0),
    dmm1vgs(0.0),
    dmm1vds(0.0),
    dmm1vbs(0.0),

    cdrain(0.0),
    cdraindrift(0.0),
    Gm    (0.0),
    Gmbs  (0.0),
    revsum(0.0),
    nrmsum(0.0),
    D1DIOoff1(0),
    D1DIOarea(1.0),
    D1DIOinitCond(0.0),
    D1DIOtemp(0.0),
    D1DIOtJctPot(0.0),
    D1DIOtJctCap(0.0),
    D1DIOtDepCap(0.0),
    D1DIOtSatCur(0.0),
    D1DIOtSatRCur(0.0),
    D1DIOtVcrit(0.0),
    D1DIOtF1(0.0),
    D1DIOtBrkdwnV(0.0),
    D1gspr(0.0),
    D1gd(0.0),
    D1cdeq(0.0),
    D1vt(0.0),
    D1vte(0.0),


    // solution,rhs vector:
    // local indices
    li_Drain                    (-1),
    li_DrainPrime               (-1),
    li_Source                   (-1),
    li_SourcePrime              (-1),
    li_Gate                     (-1),
    li_GatePrime                (-1),
    li_Bulk                     (-1),
    li_DrainDrift               (-1),
    li_D1Prime                  (-1),

    // jacobian matrix:
    // matrix and vector pointers:
    //  drain row
    ADrainEquDrainNodeOffset              (-1),
    ADrainEquSourceNodeOffset             (-1),
    ADrainEquDrainDriftNodeOffset         (-1),
    ADrainEquD1PrimeNodeOffset            (-1),
    //  gate row
    AGateEquGateNodeOffset                (-1),
    AGateEquGatePrimeNodeOffset           (-1),
    //  source row
    ASourceEquDrainNodeOffset             (-1),
    ASourceEquSourceNodeOffset            (-1),
    ASourceEquSourcePrimeNodeOffset       (-1),
    //  bulk row
    ABulkEquBulkNodeOffset                (-1),
    ABulkEquDrainPrimeNodeOffset          (-1),
    ABulkEquGatePrimeNodeOffset           (-1),
    ABulkEquSourcePrimeNodeOffset         (-1),
    // drain' row
    ADrainPrimeEquBulkNodeOffset          (-1),
    ADrainPrimeEquDrainPrimeNodeOffset    (-1),
    ADrainPrimeEquGatePrimeNodeOffset     (-1),
    ADrainPrimeEquSourcePrimeNodeOffset   (-1),
    ADrainPrimeEquDrainDriftNodeOffset    (-1),
    //  gate' row
    AGatePrimeEquGateNodeOffset           (-1),
    AGatePrimeEquBulkNodeOffset           (-1),
    AGatePrimeEquDrainPrimeNodeOffset     (-1),
    AGatePrimeEquGatePrimeNodeOffset      (-1),
    AGatePrimeEquSourcePrimeNodeOffset    (-1),
    // source' row
    ASourcePrimeEquSourceNodeOffset       (-1),
    ASourcePrimeEquBulkNodeOffset         (-1),
    ASourcePrimeEquDrainPrimeNodeOffset   (-1),
    ASourcePrimeEquGatePrimeNodeOffset    (-1),
    ASourcePrimeEquSourcePrimeNodeOffset  (-1),
    // drain drift row
    ADrainDriftEquDrainNodeOffset         (-1),
    ADrainDriftEquDrainPrimeNodeOffset    (-1),
    ADrainDriftEquDrainDriftNodeOffset    (-1),
    // D1Prime row
    AD1PrimeEquDrainNodeOffset            (-1),
    AD1PrimeEquSourceNodeOffset           (-1),
    AD1PrimeEquD1PrimeNodeOffset          (-1),

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    // F-matrix pointers:
    f_DrainEquDrainNodePtr(0),
    f_DrainEquSourceNodePtr(0),
    f_DrainEquDrainDriftNodePtr(0),
    f_DrainEquD1PrimeNodePtr(0),

    f_GateEquGateNodePtr(0),
    f_GateEquGatePrimeNodePtr(0),

    f_SourceEquDrainNodePtr(0),
    f_SourceEquSourceNodePtr(0),
    f_SourceEquSourcePrimeNodePtr(0),
    f_SourceEquD1PrimeNodePtr(0),

    f_BulkEquBulkNodePtr(0),
    f_BulkEquDrainPrimeNodePtr(0),
    f_BulkEquGatePrimeNodePtr(0),
    f_BulkEquSourcePrimeNodePtr(0),

    f_DrainPrimeEquBulkNodePtr(0),
    f_DrainPrimeEquDrainPrimeNodePtr(0),
    f_DrainPrimeEquGatePrimeNodePtr(0),
    f_DrainPrimeEquSourcePrimeNodePtr(0),
    f_DrainPrimeEquDrainDriftNodePtr(0),

    f_GatePrimeEquGateNodePtr(0),
    f_GatePrimeEquBulkNodePtr(0),
    f_GatePrimeEquDrainPrimeNodePtr(0),
    f_GatePrimeEquGatePrimeNodePtr(0),
    f_GatePrimeEquSourcePrimeNodePtr(0),

    f_SourcePrimeEquSourceNodePtr(0),
    f_SourcePrimeEquBulkNodePtr(0),
    f_SourcePrimeEquDrainPrimeNodePtr(0),
    f_SourcePrimeEquGatePrimeNodePtr(0),
    f_SourcePrimeEquSourcePrimeNodePtr(0),

    f_DrainDriftEquDrainNodePtr(0),
    f_DrainDriftEquDrainPrimeNodePtr(0),
    f_DrainDriftEquDrainDriftNodePtr(0),

    f_D1PrimeEquDrainNodePtr(0),
    f_D1PrimeEquSourceNodePtr(0),
    f_D1PrimeEquD1PrimeNodePtr(0),

    // Q-matrix pointers:
    q_DrainEquDrainNodePtr(0),
    q_DrainEquSourceNodePtr(0),
    q_DrainEquDrainDriftNodePtr(0),
    q_DrainEquD1PrimeNodePtr(0),

    q_GateEquGateNodePtr(0),
    q_GateEquGatePrimeNodePtr(0),

    q_SourceEquDrainNodePtr(0),
    q_SourceEquSourceNodePtr(0),
    q_SourceEquSourcePrimeNodePtr(0),
    q_SourceEquD1PrimeNodePtr(0),

    q_BulkEquBulkNodePtr(0),
    q_BulkEquDrainPrimeNodePtr(0),
    q_BulkEquGatePrimeNodePtr(0),
    q_BulkEquSourcePrimeNodePtr(0),

    q_DrainPrimeEquBulkNodePtr(0),
    q_DrainPrimeEquDrainPrimeNodePtr(0),
    q_DrainPrimeEquGatePrimeNodePtr(0),
    q_DrainPrimeEquSourcePrimeNodePtr(0),
    q_DrainPrimeEquDrainDriftNodePtr(0),

    q_GatePrimeEquGateNodePtr(0),
    q_GatePrimeEquBulkNodePtr(0),
    q_GatePrimeEquDrainPrimeNodePtr(0),
    q_GatePrimeEquGatePrimeNodePtr(0),
    q_GatePrimeEquSourcePrimeNodePtr(0),

    q_SourcePrimeEquSourceNodePtr(0),
    q_SourcePrimeEquBulkNodePtr(0),
    q_SourcePrimeEquDrainPrimeNodePtr(0),
    q_SourcePrimeEquGatePrimeNodePtr(0),
    q_SourcePrimeEquSourcePrimeNodePtr(0),

    q_DrainDriftEquDrainNodePtr(0),
    q_DrainDriftEquDrainPrimeNodePtr(0),
    q_DrainDriftEquDrainDriftNodePtr(0),

    q_D1PrimeEquDrainNodePtr(0),
    q_D1PrimeEquSourceNodePtr(0),
    q_D1PrimeEquD1PrimeNodePtr(0),
#endif

    vbdd(0.0),
    vbs(0.0),
    vgpdd(0.0),
    vgps(0.0),
    vdds(0.0),

    vbdd_orig(0.0),
    vbs_orig(0.0),
    vgpdd_orig(0.0),
    vgps_orig(0.0),
    vdds_orig(0.0),
    D1vd_orig(0.0),

    vbdd_old(0.0),
    vbs_old(0.0),
    vgpdd_old(0.0),
    vgps_old(0.0),
    vdds_old(0.0),
    D1vd_old(0.0),

    capgs(0.0),
    qgs(0.0),
    cqgs(0.0),
    capgdd(0.0),
    qgdd(0.0),
    cqgdd(0.0),
    capgb(0.0),
    qgb(0.0),
    cqgb(0.0),
    capbd(0.0),
    qbd(0.0),
    cqbd(0.0),
    capbs(0.0),
    qbs(0.0),
    cqbs(0.0),
    D1DIOcapCharge(0.0),
    D1DIOcapCurrent(0.0),
    D1capd(0.0),

    // local indices
    li_state_vbdd(-1),
    li_state_vbs(-1),
    li_state_vgps(-1),
    li_state_vdds(-1),
    li_state_D1vd(-1),
    li_state_capgs(-1),
    li_state_capgdd(-1),
    li_state_capgb(-1),
    li_state_qgs(-1),
    li_state_qgdd(-1),
    li_state_qgb(-1),
    li_state_qbd(-1),
    li_state_qbs(-1),

    li_state_D1DIOcapCharge(-1),
    li_state_von(-1),
    li_branch_data_d(-1),
    li_branch_data_g(-1),
    li_branch_data_s(-1),
    li_branch_data_b(-1)
{
  numStateVars = 15;
  numExtVars   = 4;
  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 4;    // this is the space to allocate if lead current or power is needed.
  devConMap.resize(4);
  devConMap[0] = 1;
  devConMap[1] = 2;
  devConMap[2] = 1;
  devConMap[3] = 3;

  if( jacStamp.empty() )
  {
    // stamp for RS!=0, RD!=0, RG!=0, RD1RS!=0
    jacStamp_D1C_DC_SC_GC.resize(9);
    jacStamp_D1C_DC_SC_GC[0].resize(4);  // Drain row
    jacStamp_D1C_DC_SC_GC[0][0]=0;       // d-d
    jacStamp_D1C_DC_SC_GC[0][1]=2;       // d-s
    jacStamp_D1C_DC_SC_GC[0][2]=7;       // d-dd
    jacStamp_D1C_DC_SC_GC[0][3]=8;       // d-d1'
    jacStamp_D1C_DC_SC_GC[1].resize(2);  // Gate row
    jacStamp_D1C_DC_SC_GC[1][0]=1;       // g-g
    jacStamp_D1C_DC_SC_GC[1][1]=5;       // g-g'
    jacStamp_D1C_DC_SC_GC[2].resize(4);  // Source row
    jacStamp_D1C_DC_SC_GC[2][0]=0;       // s-d
    jacStamp_D1C_DC_SC_GC[2][1]=2;       // s-s
    jacStamp_D1C_DC_SC_GC[2][2]=6;       // s-s'
    jacStamp_D1C_DC_SC_GC[2][3]=8;       // s-d1'
    jacStamp_D1C_DC_SC_GC[3].resize(4);  // Bulk row
    jacStamp_D1C_DC_SC_GC[3][0]=3;       // b-b
    jacStamp_D1C_DC_SC_GC[3][1]=4;       // b-d'
    jacStamp_D1C_DC_SC_GC[3][2]=5;       // b-g'
    jacStamp_D1C_DC_SC_GC[3][3]=6;       // b-s'
    jacStamp_D1C_DC_SC_GC[4].resize(5);  // Drain' row
    jacStamp_D1C_DC_SC_GC[4][0]=3;       // d'-b
    jacStamp_D1C_DC_SC_GC[4][1]=4;       // d'-d'
    jacStamp_D1C_DC_SC_GC[4][2]=5;       // d'-g'
    jacStamp_D1C_DC_SC_GC[4][3]=6;       // d'-s'
    jacStamp_D1C_DC_SC_GC[4][4]=7;       // d'-dd
    jacStamp_D1C_DC_SC_GC[5].resize(5);  // Gate' row
    jacStamp_D1C_DC_SC_GC[5][0]=1;       // g'-g
    jacStamp_D1C_DC_SC_GC[5][1]=3;       // g'-b
    jacStamp_D1C_DC_SC_GC[5][2]=4;       // g'-d'
    jacStamp_D1C_DC_SC_GC[5][3]=5;       // g'-g'
    jacStamp_D1C_DC_SC_GC[5][4]=6;       // g'-s'
    jacStamp_D1C_DC_SC_GC[6].resize(5);  // Source' row
    jacStamp_D1C_DC_SC_GC[6][0]=2;       // s'-s
    jacStamp_D1C_DC_SC_GC[6][1]=3;       // s'-b
    jacStamp_D1C_DC_SC_GC[6][2]=4;       // s'-d'
    jacStamp_D1C_DC_SC_GC[6][3]=5;       // s'-g'
    jacStamp_D1C_DC_SC_GC[6][4]=6;       // s'-s'
    jacStamp_D1C_DC_SC_GC[7].resize(3);  // DrainDrift row
    jacStamp_D1C_DC_SC_GC[7][0]=0;       // dd-d
    jacStamp_D1C_DC_SC_GC[7][1]=4;       // dd-d'
    jacStamp_D1C_DC_SC_GC[7][2]=7;       // dd-dd
    jacStamp_D1C_DC_SC_GC[8].resize(3);  // D1'pos row
    jacStamp_D1C_DC_SC_GC[8][0]=0;       // d1'-d
    jacStamp_D1C_DC_SC_GC[8][1]=2;       // d1'-s
    jacStamp_D1C_DC_SC_GC[8][2]=8;       // d1'-d1'
    jacMap_D1C_DC_SC_GC.clear();

    // First, generate the stamp for when RD1RS is zero --- it will be used
    // for all the other calculations for that case.
    // When RD1RS is zero, d1' is just S
    jacStampMap(jacStamp_D1C_DC_SC_GC, jacMap_D1C_DC_SC_GC,
                jacMap2_D1C_DC_SC_GC,
                jacStamp_DC_SC_GC,    jacMap_DC_SC_GC, jacMap2_DC_SC_GC,
                8, 2, 9);

    // Now do the cases where RD1RS is nonzero, but one of the others is
    jacStampMap(jacStamp_D1C_DC_SC_GC, jacMap_D1C_DC_SC_GC,
                jacMap2_D1C_DC_SC_GC,
                jacStamp_D1C_DC_GC, jacMap_D1C_DC_GC, jacMap2_D1C_DC_GC,
                6, 2, 9); // s' becomes same as s
    jacStampMap(jacStamp_D1C_DC_SC_GC, jacMap_D1C_DC_SC_GC,
                jacMap2_D1C_DC_SC_GC,
                jacStamp_D1C_SC_GC,    jacMap_D1C_SC_GC, jacMap2_D1C_SC_GC,
                7, 4, 9); // dd becomes same as d'
    jacStampMap(jacStamp_D1C_DC_SC_GC, jacMap_D1C_DC_SC_GC,
                jacMap2_D1C_DC_SC_GC,
                jacStamp_D1C_DC_SC,    jacMap_D1C_DC_SC, jacMap2_D1C_DC_SC,
                5, 1, 9); // g' becomes same as g

    // s' is already same as s here, so dd has become 6 instead of 7, and
    // needs to be mapped onto d' (which is still 4)
    jacStampMap(jacStamp_D1C_DC_GC, jacMap_D1C_DC_GC, jacMap2_D1C_DC_GC,
                jacStamp_D1C_GC,    jacMap_D1C_GC, jacMap2_D1C_GC, 6, 4, 9);

    // g' is same as g here, making s' 5 instead of 6, map it onto s now:
    jacStampMap(jacStamp_D1C_DC_SC, jacMap_D1C_DC_SC, jacMap2_D1C_DC_SC,
                jacStamp_D1C_DC,    jacMap_D1C_DC, jacMap2_D1C_DC, 5, 2, 9);
    // g' is same as g here, making dd 6 instead of 7, map it onto d' now:
    jacStampMap(jacStamp_D1C_DC_SC, jacMap_D1C_DC_SC, jacMap2_D1C_DC_SC,
                jacStamp_D1C_SC,    jacMap_D1C_SC, jacMap2_D1C_SC, 6, 4, 9);

    // g'=g and s'=s, so dd is 5 instead of 7, map it onto d'
    jacStampMap(jacStamp_D1C_DC, jacMap_D1C_DC, jacMap2_D1C_DC,
                jacStamp_D1C,    jacMap_D1C, jacMap2_D1C, 5, 4, 9);

    // Now do the ones where RD1RS is zero, starting from the case where
    // everything *but* RD1RS is nonzero.  These follow exactly the same
    // pattern as above, but starting from _DC_SC_GC instead of _D1C_DC_SC_GC
    // This only works out this way because d1' is the last row, and nothing
    // gets moved up when it goes away.
    // I won't repeat the comments for each block.
    // only one row goes away:
    jacStampMap(jacStamp_DC_SC_GC, jacMap_DC_SC_GC, jacMap2_DC_SC_GC,
                jacStamp_DC_GC,    jacMap_DC_GC, jacMap2_DC_GC, 6, 2, 9);
    jacStampMap(jacStamp_DC_SC_GC, jacMap_DC_SC_GC, jacMap2_DC_SC_GC,
                jacStamp_SC_GC,    jacMap_SC_GC, jacMap2_SC_GC, 7, 4, 9);
    jacStampMap(jacStamp_DC_SC_GC, jacMap_DC_SC_GC, jacMap2_DC_SC_GC,
                jacStamp_DC_SC,    jacMap_DC_SC, jacMap2_DC_SC, 5, 1, 9);

    // s'==s already
    jacStampMap(jacStamp_DC_GC, jacMap_DC_GC, jacMap2_DC_GC,
                jacStamp_GC,    jacMap_GC, jacMap2_GC, 6, 4, 9);

    // g'=g already
    jacStampMap(jacStamp_DC_SC, jacMap_DC_SC, jacMap2_DC_SC,
                jacStamp_DC,    jacMap_DC, jacMap2_DC, 5, 2, 9);
    jacStampMap(jacStamp_DC_SC, jacMap_DC_SC, jacMap2_DC_SC,
                jacStamp_SC,    jacMap_SC, jacMap2_SC, 6, 4, 9);

    // s'=s, g'=g already
    jacStampMap(jacStamp_DC, jacMap_DC, jacMap2_DC,
                jacStamp,    jacMap, jacMap2, 5, 4, 9);
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

  // if options scale has been set in the netlist, apply it.
  applyScale ();

  // call updateTemp which may change model parameters if
  // temperature effects are present
  processParams ();

  numIntVars = 1 + (((sourceCond == 0.0)?0:1)+((drainCond == 0.0) ? 0:1)
                 +  ((gateCond == 0.0) ? 0:1))
                 +  ((model_.D1DIOconductance == 0.0) ? 0:1);

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
  li_DrainPrime = intLIDVec[intLoc++];

  if( gateCond )
    li_GatePrime = intLIDVec[intLoc++];
  else
    li_GatePrime = li_Gate;

  if( sourceCond )
    li_SourcePrime = intLIDVec[intLoc++];
  else
    li_SourcePrime = li_Source;

  if( drainCond )
    li_DrainDrift = intLIDVec[intLoc++];
  else
    li_DrainDrift = li_DrainPrime;

  if( model_.D1DIOconductance )
    li_D1Prime = intLIDVec[intLoc++];
  else
    li_D1Prime = li_Source;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "\n variable local indices:\n";
    Xyce::dout() << "  li_Drain       = " << li_Drain << std::endl;
    Xyce::dout() << "  li_Source      = " << li_Source << std::endl;
    Xyce::dout() << "  li_Gate        = " << li_Gate << std::endl;
    Xyce::dout() << "  li_Bulk        = " << li_Bulk << std::endl;
    Xyce::dout() << "  li_DrainPrime  = " << li_DrainPrime << std::endl;
    Xyce::dout() << "  li_GatePrime   = " << li_GatePrime << std::endl;
    Xyce::dout() << "  li_SourcePrime = " << li_SourcePrime << std::endl;
    Xyce::dout() << "  li_DrainDrift  = " << li_DrainDrift << std::endl;
    Xyce::dout() << "  li_D1Prime       = " << li_D1Prime << std::endl;

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
  addInternalNode(symbol_table, li_DrainPrime, getName(), "drainprime");

  if (li_GatePrime != li_Gate )
    addInternalNode(symbol_table, li_GatePrime, getName(), "gateprime");

  if (li_SourcePrime != li_Source )
    addInternalNode(symbol_table, li_SourcePrime, getName(), "sourceprime");

  if (li_DrainDrift != li_DrainPrime )
    addInternalNode(symbol_table, li_DrainDrift, getName(), "draindrift");

  if (li_D1Prime != li_Source )
    addInternalNode(symbol_table, li_D1Prime, getName(), "d1pos");

  if (loadLeadCurrent)
  {
    addBranchDataNode(symbol_table, li_branch_data_d, getName(), "BRANCH_DD");
    addBranchDataNode(symbol_table, li_branch_data_s, getName(), "BRANCH_DS");
    addBranchDataNode(symbol_table, li_branch_data_g, getName(), "BRANCH_DG");
    addBranchDataNode(symbol_table, li_branch_data_b, getName(), "BRANCH_DB");
  }
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
  li_state_vbdd = staLIDVec[lid++];
  li_state_vbs  = staLIDVec[lid++];
  li_state_vgps = staLIDVec[lid++];
  li_state_vdds = staLIDVec[lid++];
  li_state_D1vd = staLIDVec[lid++];

  li_state_qgs   = staLIDVec[lid++];
  //li_state_cqgs  = staLIDVec[lid++];
  li_state_qgdd  = staLIDVec[lid++];
  //li_state_cqgdd = staLIDVec[lid++];
  li_state_qgb   = staLIDVec[lid++];
  //li_state_cqgb  = staLIDVec[lid++];

  li_state_capgs  = staLIDVec[lid++];
  li_state_capgdd = staLIDVec[lid++];
  li_state_capgb  = staLIDVec[lid++];

  li_state_qbd = staLIDVec[lid++];
  //li_state_cqbd = staLIDVec[lid++];
  li_state_qbs = staLIDVec[lid++];
  //li_state_cqbs = staLIDVec[lid++];

  li_state_D1DIOcapCharge  = staLIDVec[lid++];
  //li_state_D1DIOcapCurrent = staLIDVec[lid++];

  li_state_von = staLIDVec[lid++];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "  State local indices:" << std::endl;
    Xyce::dout() << std::endl;

    Xyce::dout() << "  li_state_vbdd          = " << li_state_vbdd << "\n";
    Xyce::dout() << "  li_state_vbs           = " << li_state_vbs << "\n";
    Xyce::dout() << "  li_state_vgps          = " << li_state_vgps << "\n";
    Xyce::dout() << "  li_state_vdds          = " << li_state_vdds << "\n";
    Xyce::dout() << "  li_state_qgs           = " << li_state_qgs  << "\n";
    //Xyce::dout() << "  li_state_cqgs          = " << li_state_cqgs << "\n";
    Xyce::dout() << "  li_state_capgs         = " << li_state_capgs << "\n";
    Xyce::dout() << "  li_state_capgdd        = " << li_state_capgdd << "\n";
    Xyce::dout() << "  li_state_capgb         = " << li_state_capgb << "\n";
    Xyce::dout() << "  li_state_qgdd          = " << li_state_qgdd << "\n";
    //Xyce::dout() << "  li_state_cqgdd         = " << li_state_cqgdd << "\n";
    Xyce::dout() << "  li_state_qgb           = " << li_state_qgb << "\n";
    //Xyce::dout() << "  li_state_cqgb          = " << li_state_cqgb << "\n";
    Xyce::dout() << "  li_state_qbs           = " << li_state_qbs << "\n";
    //Xyce::dout() << "  li_state_cqbs          = " << li_state_cqbs << "\n";
    Xyce::dout() << "  li_state_qbd           = " << li_state_qbd << "\n";
    //Xyce::dout() << "  li_state_cqbd          = " << li_state_cqbd << "\n";
    Xyce::dout() << "  li_state_D1DIOcapCharge  = " << li_state_D1DIOcapCharge << "\n";
    //Xyce::dout() << "  li_state_D1DIOcapCurrent = " << li_state_D1DIOcapCurrent << "\n";
    Xyce::dout() << "  li_state_von           = " << li_state_von << "\n";

    Xyce::dout() << std::endl;
    Xyce::dout() << section_divider << std::endl;
  }

}


//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 4/23/2013
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs( const std::vector<int> & stoLIDVecRef )
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
    li_branch_data_d = branchLIDVecRef[0];
    li_branch_data_g = branchLIDVecRef[1];
    li_branch_data_s = branchLIDVecRef[2];
    li_branch_data_b = branchLIDVecRef[3];
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
  if( drainCond != 0.0 && gateCond != 0.0 && sourceCond != 0.0)
    return (model_.D1DIOconductance != 0)?jacStamp_D1C_DC_SC_GC:jacStamp_DC_SC_GC;
  else if( drainCond != 0.0 && gateCond != 0.0 && sourceCond == 0.0 )
    return (model_.D1DIOconductance != 0)?jacStamp_D1C_DC_GC:jacStamp_DC_GC;
  else if( drainCond == 0.0 && gateCond != 0.0 && sourceCond != 0.0 )
    return (model_.D1DIOconductance != 0)?jacStamp_D1C_SC_GC:jacStamp_SC_GC;
  else if( drainCond != 0.0 && gateCond == 0.0 && sourceCond != 0.0 )
    return (model_.D1DIOconductance != 0)?jacStamp_D1C_DC_SC:jacStamp_DC_SC;
  else if( drainCond == 0.0 && gateCond != 0.0 && sourceCond == 0.0 )
    return (model_.D1DIOconductance != 0)?jacStamp_D1C_GC:jacStamp_GC;
  else if( drainCond != 0.0 && gateCond == 0.0 && sourceCond == 0.0 )
    return (model_.D1DIOconductance != 0)?jacStamp_D1C_DC:jacStamp_DC;
  else if( drainCond == 0.0 && gateCond == 0.0 && sourceCond != 0.0 )
    return (model_.D1DIOconductance != 0)?jacStamp_D1C_SC:jacStamp_SC;
  else if( drainCond == 0.0 && gateCond == 0.0 && sourceCond == 0.0 )
    return (model_.D1DIOconductance != 0)?jacStamp_D1C:jacStamp;
  else
    return (model_.D1DIOconductance != 0)?jacStamp_D1C:jacStamp;
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
  double d1dcond=(model_.D1DIOconductance);
  if (gateCond != 0.0)
  {
    if (drainCond != 0.0)
    {
      if (sourceCond != 0.0)
      {
        if (d1dcond != 0)
        {
          map =  jacMap_D1C_DC_SC_GC;
          map2 = jacMap2_D1C_DC_SC_GC;
        }
        else
        {
          map =  jacMap_DC_SC_GC;
          map2 = jacMap2_DC_SC_GC;
        }
      }
      else
      {
        if (d1dcond != 0)
        {
          map = jacMap_D1C_DC_GC;
          map2 = jacMap2_D1C_DC_GC;
        }
        else
        {
          map = jacMap_DC_GC;
          map2 = jacMap2_DC_GC;
        }
      }
    }
    else
    {
      if (sourceCond != 0.0)
      {
        if (d1dcond != 0)
        {
          map = jacMap_D1C_SC_GC;
          map2 = jacMap2_D1C_SC_GC;
        }
        else
        {
          map = jacMap_SC_GC;
          map2 = jacMap2_SC_GC;
        }
      }
      else
      {
        if (d1dcond != 0)
        {
          map = jacMap_D1C_GC;
          map2 = jacMap2_D1C_GC;
        }
        else
        {
          map = jacMap_GC;
          map2 = jacMap2_GC;
        }
      }
    }
  }
  else
  {
    if (drainCond != 0.0)
    {
      if (sourceCond != 0.0)
      {
        if (d1dcond != 0)
        {
          map = jacMap_D1C_DC_SC;
          map2 = jacMap2_D1C_DC_SC;
        }
        else
        {
          map = jacMap_DC_SC;
          map2 = jacMap2_DC_SC;
        }
      }
      else
      {
        if (d1dcond != 0)
        {
          map = jacMap_D1C_DC;
          map2 = jacMap2_D1C_DC;
        }
        else
        {
          map = jacMap_DC;
          map2 = jacMap2_DC;
        }
      }
    }
    else
    {
      if (sourceCond != 0.0)
      {
        if (d1dcond != 0)
        {
          map = jacMap_D1C_SC;
          map2 = jacMap2_D1C_SC;
        }
        else
        {
          map = jacMap_SC;
          map2 = jacMap2_SC;
        }
      }
      else
      {
        if (d1dcond != 0)
        {
          map = jacMap_D1C;
          map2 = jacMap2_D1C;
        }
        else
        {
          map = jacMap;
          map2 = jacMap2;
        }
      }
    }
  }

  ADrainEquDrainNodeOffset             = jacLIDVec[map[0]][map2[0][0]];
  ADrainEquSourceNodeOffset            = jacLIDVec[map[0]][map2[0][1]];
  ADrainEquDrainDriftNodeOffset        = jacLIDVec[map[0]][map2[0][2]];
  ADrainEquD1PrimeNodeOffset           = jacLIDVec[map[0]][map2[0][3]];

  AGateEquGateNodeOffset               = jacLIDVec[map[1]][map2[1][0]];
  AGateEquGatePrimeNodeOffset          = jacLIDVec[map[1]][map2[1][1]];

  ASourceEquDrainNodeOffset            = jacLIDVec[map[2]][map2[2][0]];
  ASourceEquSourceNodeOffset           = jacLIDVec[map[2]][map2[2][1]];
  ASourceEquSourcePrimeNodeOffset      = jacLIDVec[map[2]][map2[2][2]];
  ASourceEquD1PrimeNodeOffset          = jacLIDVec[map[2]][map2[2][3]];

  ABulkEquBulkNodeOffset               = jacLIDVec[map[3]][map2[3][0]];
  ABulkEquDrainPrimeNodeOffset         = jacLIDVec[map[3]][map2[3][1]];
  ABulkEquGatePrimeNodeOffset          = jacLIDVec[map[3]][map2[3][2]];
  ABulkEquSourcePrimeNodeOffset        = jacLIDVec[map[3]][map2[3][3]];

  ADrainPrimeEquBulkNodeOffset         = jacLIDVec[map[4]][map2[4][0]];
  ADrainPrimeEquDrainPrimeNodeOffset   = jacLIDVec[map[4]][map2[4][1]];
  ADrainPrimeEquGatePrimeNodeOffset    = jacLIDVec[map[4]][map2[4][2]];
  ADrainPrimeEquSourcePrimeNodeOffset  = jacLIDVec[map[4]][map2[4][3]];
  ADrainPrimeEquDrainDriftNodeOffset   = jacLIDVec[map[4]][map2[4][4]];

  AGatePrimeEquGateNodeOffset          = jacLIDVec[map[5]][map2[5][0]];
  AGatePrimeEquBulkNodeOffset          = jacLIDVec[map[5]][map2[5][1]];
  AGatePrimeEquDrainPrimeNodeOffset    = jacLIDVec[map[5]][map2[5][2]];
  AGatePrimeEquGatePrimeNodeOffset     = jacLIDVec[map[5]][map2[5][3]];
  AGatePrimeEquSourcePrimeNodeOffset   = jacLIDVec[map[5]][map2[5][4]];

  ASourcePrimeEquSourceNodeOffset      = jacLIDVec[map[6]][map2[6][0]];
  ASourcePrimeEquBulkNodeOffset        = jacLIDVec[map[6]][map2[6][1]];
  ASourcePrimeEquDrainPrimeNodeOffset  = jacLIDVec[map[6]][map2[6][2]];
  ASourcePrimeEquGatePrimeNodeOffset   = jacLIDVec[map[6]][map2[6][3]];
  ASourcePrimeEquSourcePrimeNodeOffset = jacLIDVec[map[6]][map2[6][4]];

  ADrainDriftEquDrainNodeOffset        = jacLIDVec[map[7]][map2[7][0]];
  ADrainDriftEquDrainPrimeNodeOffset   = jacLIDVec[map[7]][map2[7][1]];
  ADrainDriftEquDrainDriftNodeOffset   = jacLIDVec[map[7]][map2[7][2]];

  AD1PrimeEquDrainNodeOffset           = jacLIDVec[map[8]][map2[8][0]];
  AD1PrimeEquSourceNodeOffset          = jacLIDVec[map[8]][map2[8][1]];
  AD1PrimeEquD1PrimeNodeOffset         = jacLIDVec[map[8]][map2[8][2]];
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

  f_DrainEquDrainNodePtr             = 	  &(dFdx[li_Drain][ADrainEquDrainNodeOffset             ]);
  f_DrainEquSourceNodePtr            = 	  &(dFdx[li_Drain][ADrainEquSourceNodeOffset            ]);
  f_DrainEquDrainDriftNodePtr        = 	  &(dFdx[li_Drain][ADrainEquDrainDriftNodeOffset        ]);
  f_DrainEquD1PrimeNodePtr           = 	  &(dFdx[li_Drain][ADrainEquD1PrimeNodeOffset           ]);

  f_GateEquGateNodePtr               = 	  &(dFdx[li_Gate][AGateEquGateNodeOffset               ]);
  f_GateEquGatePrimeNodePtr          = 	  &(dFdx[li_Gate][AGateEquGatePrimeNodeOffset          ]);

  f_SourceEquDrainNodePtr            = 	  &(dFdx[li_Source][ASourceEquDrainNodeOffset            ]);
  f_SourceEquSourceNodePtr           = 	  &(dFdx[li_Source][ASourceEquSourceNodeOffset           ]);
  f_SourceEquSourcePrimeNodePtr      = 	  &(dFdx[li_Source][ASourceEquSourcePrimeNodeOffset      ]);
  f_SourceEquD1PrimeNodePtr          = 	  &(dFdx[li_Source][ASourceEquD1PrimeNodeOffset          ]);

  f_BulkEquBulkNodePtr               = 	  &(dFdx[li_Bulk][ABulkEquBulkNodeOffset               ]);
  f_BulkEquDrainPrimeNodePtr         = 	  &(dFdx[li_Bulk][ABulkEquDrainPrimeNodeOffset         ]);
  f_BulkEquGatePrimeNodePtr          = 	  &(dFdx[li_Bulk][ABulkEquGatePrimeNodeOffset          ]);
  f_BulkEquSourcePrimeNodePtr        = 	  &(dFdx[li_Bulk][ABulkEquSourcePrimeNodeOffset        ]);

  f_DrainPrimeEquBulkNodePtr         = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset         ]);
  f_DrainPrimeEquDrainPrimeNodePtr   = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset   ]);
  f_DrainPrimeEquGatePrimeNodePtr    = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquGatePrimeNodeOffset    ]);
  f_DrainPrimeEquSourcePrimeNodePtr  = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset  ]);
  f_DrainPrimeEquDrainDriftNodePtr   = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquDrainDriftNodeOffset   ]);

  f_GatePrimeEquGateNodePtr          = 	  &(dFdx[li_GatePrime][AGatePrimeEquGateNodeOffset          ]);
  f_GatePrimeEquBulkNodePtr          = 	  &(dFdx[li_GatePrime][AGatePrimeEquBulkNodeOffset          ]);
  f_GatePrimeEquDrainPrimeNodePtr    = 	  &(dFdx[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset    ]);
  f_GatePrimeEquGatePrimeNodePtr     = 	  &(dFdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset     ]);
  f_GatePrimeEquSourcePrimeNodePtr   = 	  &(dFdx[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset   ]);

  f_SourcePrimeEquSourceNodePtr      = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset      ]);
  f_SourcePrimeEquBulkNodePtr        = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset        ]);
  f_SourcePrimeEquDrainPrimeNodePtr  = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset  ]);
  f_SourcePrimeEquGatePrimeNodePtr   = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquGatePrimeNodeOffset   ]);
  f_SourcePrimeEquSourcePrimeNodePtr = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset ]);

  f_DrainDriftEquDrainNodePtr        = 	  &(dFdx[li_DrainDrift][ADrainDriftEquDrainNodeOffset        ]);
  f_DrainDriftEquDrainPrimeNodePtr   = 	  &(dFdx[li_DrainDrift][ADrainDriftEquDrainPrimeNodeOffset   ]);
  f_DrainDriftEquDrainDriftNodePtr   = 	  &(dFdx[li_DrainDrift][ADrainDriftEquDrainDriftNodeOffset   ]);

  f_D1PrimeEquDrainNodePtr           = 	  &(dFdx[li_D1Prime][AD1PrimeEquDrainNodeOffset           ]);
  f_D1PrimeEquSourceNodePtr          = 	  &(dFdx[li_D1Prime][AD1PrimeEquSourceNodeOffset          ]);
  f_D1PrimeEquD1PrimeNodePtr         = 	  &(dFdx[li_D1Prime][AD1PrimeEquD1PrimeNodeOffset         ]);



  q_DrainEquDrainNodePtr             = 	  &(dQdx[li_Drain][ADrainEquDrainNodeOffset             ]);
  q_DrainEquSourceNodePtr            = 	  &(dQdx[li_Drain][ADrainEquSourceNodeOffset            ]);
  q_DrainEquDrainDriftNodePtr        = 	  &(dQdx[li_Drain][ADrainEquDrainDriftNodeOffset        ]);
  q_DrainEquD1PrimeNodePtr           = 	  &(dQdx[li_Drain][ADrainEquD1PrimeNodeOffset           ]);

  q_GateEquGateNodePtr               = 	  &(dQdx[li_Gate][AGateEquGateNodeOffset               ]);
  q_GateEquGatePrimeNodePtr          = 	  &(dQdx[li_Gate][AGateEquGatePrimeNodeOffset          ]);

  q_SourceEquDrainNodePtr            = 	  &(dQdx[li_Source][ASourceEquDrainNodeOffset            ]);
  q_SourceEquSourceNodePtr           = 	  &(dQdx[li_Source][ASourceEquSourceNodeOffset           ]);
  q_SourceEquSourcePrimeNodePtr      = 	  &(dQdx[li_Source][ASourceEquSourcePrimeNodeOffset      ]);
  q_SourceEquD1PrimeNodePtr          = 	  &(dQdx[li_Source][ASourceEquD1PrimeNodeOffset          ]);

  q_BulkEquBulkNodePtr               = 	  &(dQdx[li_Bulk][ABulkEquBulkNodeOffset               ]);
  q_BulkEquDrainPrimeNodePtr         = 	  &(dQdx[li_Bulk][ABulkEquDrainPrimeNodeOffset         ]);
  q_BulkEquGatePrimeNodePtr          = 	  &(dQdx[li_Bulk][ABulkEquGatePrimeNodeOffset          ]);
  q_BulkEquSourcePrimeNodePtr        = 	  &(dQdx[li_Bulk][ABulkEquSourcePrimeNodeOffset        ]);

  q_DrainPrimeEquBulkNodePtr         = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset         ]);
  q_DrainPrimeEquDrainPrimeNodePtr   = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset   ]);
  q_DrainPrimeEquGatePrimeNodePtr    = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquGatePrimeNodeOffset    ]);
  q_DrainPrimeEquSourcePrimeNodePtr  = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset  ]);
  q_DrainPrimeEquDrainDriftNodePtr   = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquDrainDriftNodeOffset   ]);

  q_GatePrimeEquGateNodePtr          = 	  &(dQdx[li_GatePrime][AGatePrimeEquGateNodeOffset          ]);
  q_GatePrimeEquBulkNodePtr          = 	  &(dQdx[li_GatePrime][AGatePrimeEquBulkNodeOffset          ]);
  q_GatePrimeEquDrainPrimeNodePtr    = 	  &(dQdx[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset    ]);
  q_GatePrimeEquGatePrimeNodePtr     = 	  &(dQdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset     ]);
  q_GatePrimeEquSourcePrimeNodePtr   = 	  &(dQdx[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset   ]);

  q_SourcePrimeEquSourceNodePtr      = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset      ]);
  q_SourcePrimeEquBulkNodePtr        = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset        ]);
  q_SourcePrimeEquDrainPrimeNodePtr  = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset  ]);
  q_SourcePrimeEquGatePrimeNodePtr   = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquGatePrimeNodeOffset   ]);
  q_SourcePrimeEquSourcePrimeNodePtr = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset ]);

  q_DrainDriftEquDrainNodePtr        = 	  &(dQdx[li_DrainDrift][ADrainDriftEquDrainNodeOffset        ]);
  q_DrainDriftEquDrainPrimeNodePtr   = 	  &(dQdx[li_DrainDrift][ADrainDriftEquDrainPrimeNodeOffset   ]);
  q_DrainDriftEquDrainDriftNodePtr   = 	  &(dQdx[li_DrainDrift][ADrainDriftEquDrainDriftNodeOffset   ]);

  q_D1PrimeEquDrainNodePtr           = 	  &(dQdx[li_D1Prime][AD1PrimeEquDrainNodeOffset           ]);
  q_D1PrimeEquSourceNodePtr          = 	  &(dQdx[li_D1Prime][AD1PrimeEquSourceNodeOffset          ]);
  q_D1PrimeEquD1PrimeNodePtr         = 	  &(dQdx[li_D1Prime][AD1PrimeEquD1PrimeNodeOffset         ]);




#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 voltage source instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) + B(t) = 0
//
//                 The "Q" vector contains charges and fluxes, mostly.
//                 The voltage source will not make any contributions to Q,
//                 so this function does nothing.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/02/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  double * qVec = extData.daeQVectorRawPtr;
  double * dQdxdVp = extData.dQdxdVpVectorRawPtr;

  double coef_Jdxp(0.0);

  double ceqbs(0.0), ceqbd(0.0), ceqgb(0.0),  ceqgs(0.0),  ceqgdd(0.0);  // 3f5 vars
  double Dtype =model_.dtype;

  ceqbs = Dtype*qbs;
  ceqbd = Dtype*qbd;
  ceqgb = Dtype*qgb;
  ceqgs = Dtype*qgs;
  ceqgdd = Dtype*qgdd;

  qVec[li_Bulk]        += (ceqbs + ceqbd - ceqgb);
  qVec[li_DrainPrime]  += -(ceqbd + ceqgdd);
  qVec[li_GatePrime]   += (ceqgs+ceqgdd+ceqgb);
  qVec[li_SourcePrime] += (-(ceqbs + ceqgs));
  qVec[li_D1Prime] +=  D1DIOcapCharge;
  qVec[li_Drain]   += -D1DIOcapCharge;
  // voltlim section:
  if (!origFlag)
  {
    // bulk
    coef_Jdxp = Dtype*(-Capgb*(vgps-vgps_orig-vbs+vbs_orig)
                       + (+Capgb)*(vbdd-vbdd_orig)
                       + (+capbs)*(vbs-vbs_orig));
    dQdxdVp[li_Bulk] += coef_Jdxp;

    // drain-prime
    coef_Jdxp = Dtype*(-Capgdd*(vgpdd-vgpdd_orig)-capbd*(vbdd-vbdd_orig));
    dQdxdVp[li_DrainPrime] += coef_Jdxp;

    // gate-prime
    coef_Jdxp = Dtype*(Capgdd*(vgpdd-vgpdd_orig)+Capgs*(vgps-vgps_orig)+
                       Capgb*(vgps-vgps_orig-vbs+vbs_orig));
    dQdxdVp[li_GatePrime] += coef_Jdxp;

    // source-prime
    coef_Jdxp = Dtype*(-Capgs*(vgps-vgps_orig)-capbs*(vbs-vbs_orig));
    dQdxdVp[li_SourcePrime] += coef_Jdxp;

    coef_Jdxp = -D1capd*(D1vd-D1vd_orig);
    dQdxdVp[li_D1Prime] -= coef_Jdxp;
    dQdxdVp[li_Drain]   += coef_Jdxp;
  }

  if( loadLeadCurrent )
  {
    double * leadQ = extData.nextLeadCurrQCompRawPtr;
    leadQ[li_branch_data_d] = -D1DIOcapCharge;
    leadQ[li_branch_data_s] = 0;
    leadQ[li_branch_data_g] = 0;
    leadQ[li_branch_data_b] = (ceqbs + ceqbd - ceqgb);
    // case where optional nodes become external nodes:
    if( !gateCond )
    {
      // G' is G
      leadQ[li_branch_data_g] += (ceqgs+ceqgdd+ceqgb);
    }

    if( !sourceCond )
    {
      // S' is S
      leadQ[li_branch_data_s] += (-(ceqbs + ceqgs));
    }

    if( !drainCond )
    {
      // Ddrift is D'
    }

    if( !model_.D1DIOconductance )
    {
      // D1' is S
      leadQ[li_branch_data_s] += D1DIOcapCharge;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 VDMOS instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/02/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * solVec = extData.nextSolVectorRawPtr;
  double * fVec = extData.daeFVectorRawPtr;
  double * dFdxdVp = extData.dFdxdVpVectorRawPtr;

  double coef_Jdxp(0.0), gd_Jdxp(0.0);
  double gmin1 = getDeviceOptions().gmin;
  double Dtype =model_.dtype;

  double ceqbs = Dtype*(cbs+cqbs);
  double ceqbd = Dtype*(cbd+cqbd);

  double D1current = Dtype*D1cdeq; // don't add in the capacitor stuff
  fVec[li_Drain] += (Idraindrift + Irdsshunt - D1current);

  if (Igate != 0.0)
  {
    fVec[li_Gate]  += Igate;
    fVec[li_GatePrime]   += -Igate;
  }
  fVec[li_Source]      += (Isource - Irdsshunt + Ird1rs);
  fVec[li_Bulk]        += (ceqbs + ceqbd);
  fVec[li_DrainPrime]  += (-Idrain-(ceqbd - cdreq));
  fVec[li_SourcePrime] += (-Isource-(ceqbs + cdreq));
  fVec[li_DrainDrift]  += (Idrain - Idraindrift);
  fVec[li_D1Prime]     += (D1current - Ird1rs);

  // limiter terms:
  if (!origFlag)
  {
    // bulk
    coef_Jdxp = Dtype*( + ((gbd-gmin1))*(vbdd-vbdd_orig)
                        + ((gbs-gmin1))*(vbs-vbs_orig));
    dFdxdVp[li_Bulk] += coef_Jdxp;

    // drain-prime
    coef_Jdxp = Dtype*(-((gbd-gmin1))*
                       (vbdd-vbdd_orig)+gdds*(vdds-vdds_orig)
                       +Gm*((mode>0)?(vgps-vgps_orig):(vgpdd-vgpdd_orig))
                       +Gmbs*((mode>0)?(vbs-vbs_orig):(vbdd-vbdd_orig)));
    dFdxdVp[li_DrainPrime] += coef_Jdxp;

    // source prime
    coef_Jdxp = Dtype*(-((gbs-gmin1))*(vbs-vbs_orig)
                       -gdds*(vdds-vdds_orig)
                       -Gm*((mode>0)?(vgps-vgps_orig):(vgpdd-vgpdd_orig))
                       -Gmbs*((mode>0)?(vbs-vbs_orig):(vbdd-vbdd_orig)));
    dFdxdVp[li_SourcePrime] += coef_Jdxp;

    gd_Jdxp = -D1gd * (D1vd-D1vd_orig);
    dFdxdVp[li_Drain] += gd_Jdxp;
    dFdxdVp[li_D1Prime] -= gd_Jdxp;
  }

  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    double * junctionV = extData.nextJunctionVCompRawPtr;
    leadF[li_branch_data_d] = (Idraindrift + Irdsshunt - D1current);
    leadF[li_branch_data_s] = (Isource - Irdsshunt + Ird1rs);
    leadF[li_branch_data_g] = 0;
    leadF[li_branch_data_b] = (ceqbs + ceqbd);
    // case where optional nodes become external nodes:
    if (Igate != 0.0)
    {
      leadF[li_branch_data_g]  += Igate;
    }
    if( !gateCond )
    {
      // G' is G
      leadF[li_branch_data_g] += -Igate;
    }

    if( !sourceCond )
    {
      // S' is S
      leadF[li_branch_data_s] += (-Isource-(ceqbs + cdreq));
    }

    if( !drainCond )
    {
      // Ddrift is D'
    }

    if( !model_.D1DIOconductance )
    {
      // D1' is S
      leadF[li_branch_data_s] += (D1current - Ird1rs);
    }

    junctionV[li_branch_data_d] = solVec[li_Drain] - solVec[li_Source];
    junctionV[li_branch_data_g] = solVec[li_Gate] - solVec[li_Source];
    junctionV[li_branch_data_s] = 0.0;
    junctionV[li_branch_data_b] = 0.0 ;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 VDMOS instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/02/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);

  if (getSolverState().dcopFlag) return true;

  dQdx[li_Bulk][ABulkEquBulkNodeOffset] += +capbs+capbd+Capgb;
  dQdx[li_Bulk][ABulkEquDrainPrimeNodeOffset] -= +capbd;
  dQdx[li_Bulk][ABulkEquGatePrimeNodeOffset] -= Capgb;
  dQdx[li_Bulk][ABulkEquSourcePrimeNodeOffset] -= +capbs;

  dQdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset] += -capbd;
  dQdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset] += +capbd+Capgdd;
  dQdx[li_DrainPrime][ADrainPrimeEquGatePrimeNodeOffset] += -Capgdd;

  dQdx[li_GatePrime][AGatePrimeEquBulkNodeOffset] -= Capgb;
  dQdx[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset] -= Capgdd;
  dQdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset] += +Capgdd+Capgs+Capgb;
  dQdx[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset] -= Capgs;

  dQdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset] -= +capbs;
  dQdx[li_SourcePrime][ASourcePrimeEquGatePrimeNodeOffset] -= Capgs;
  dQdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset] += +capbs+Capgs;

  dQdx[li_Drain][ADrainEquDrainNodeOffset] += D1capd;
  dQdx[li_Drain][ADrainEquD1PrimeNodeOffset] -= D1capd;

  dQdx[li_D1Prime][AD1PrimeEquDrainNodeOffset] -= D1capd;
  dQdx[li_D1Prime][AD1PrimeEquD1PrimeNodeOffset] += D1capd;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 instance.
//
// Special Notes : The F-vector is an algebraic constaint.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/02/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Drain][ADrainEquDrainNodeOffset] += dIdd_dVd + gdsshunt + D1gd;
  dFdx[li_Drain][ADrainEquSourceNodeOffset] -= gdsshunt;
  dFdx[li_Drain][ADrainEquDrainDriftNodeOffset] -= dIdd_dVd;
  dFdx[li_Drain][ADrainEquD1PrimeNodeOffset] -= D1gd;

  if (gateCond != 0)
  {
    dFdx[li_Gate][AGateEquGateNodeOffset] += gateCond;
    dFdx[li_Gate][AGateEquGatePrimeNodeOffset] -= gateCond;
  }

  dFdx[li_Source][ASourceEquDrainNodeOffset] -= gdsshunt;
  dFdx[li_Source][ASourceEquSourceNodeOffset] += sourceCond + gdsshunt + D1gspr;
  dFdx[li_Source][ASourceEquD1PrimeNodeOffset] += -D1gspr;

  if (sourceCond != 0.0)
  {
    dFdx[li_Source][ASourceEquSourcePrimeNodeOffset] -= sourceCond;
  }

  dFdx[li_Bulk][ABulkEquBulkNodeOffset] += gbs+gbd;
  dFdx[li_Bulk][ABulkEquDrainPrimeNodeOffset] -= gbd;

  dFdx[li_Bulk][ABulkEquSourcePrimeNodeOffset] -= gbs;

  dFdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset] += -gbd+Gmbs;
  dFdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset] += drainCond+gdds+gbd+revsum;
  dFdx[li_DrainPrime][ADrainPrimeEquGatePrimeNodeOffset] += Gm;
  dFdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset] += -gdds-nrmsum;
  if (drainCond != 0.0)
  {
    dFdx[li_DrainPrime][ADrainPrimeEquDrainDriftNodeOffset] -= drainCond;
  }

  if (gateCond != 0)
  {
    dFdx[li_GatePrime][AGatePrimeEquGateNodeOffset] -= gateCond;
    dFdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset] += gateCond;
  }

  if (sourceCond != 0.0)
  {
    dFdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset] -= sourceCond;
  }
  dFdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset] -= gbs+Gmbs;
  dFdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset] -= gdds+revsum;
  dFdx[li_SourcePrime][ASourcePrimeEquGatePrimeNodeOffset] -= Gm;
  dFdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset] += sourceCond+gdds+gbs+nrmsum;

  dFdx[li_DrainDrift][ADrainDriftEquDrainNodeOffset] -= dIdd_dVd;

  if (drainCond != 0.0)
  {
    dFdx[li_DrainDrift][ADrainDriftEquDrainPrimeNodeOffset] -= drainCond;
  }
  dFdx[li_DrainDrift][ADrainDriftEquDrainDriftNodeOffset] += dIdd_dVd + drainCond;

  dFdx[li_D1Prime][AD1PrimeEquDrainNodeOffset] -= D1gd;
  dFdx[li_D1Prime][AD1PrimeEquSourceNodeOffset] += -D1gspr;
  dFdx[li_D1Prime][AD1PrimeEquD1PrimeNodeOffset] += D1gd + D1gspr;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : pmc
// Creation Date : 2/17/2004
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;

  double * solVec = extData.nextSolVectorRawPtr;

  // 3f5 likes to use the same variable names in local variables and in
  // structures.  Messes with us!  Define some local versions with capitals
  // instead
#define ISUBMOD  model_.isubmod

  double Von(0.0);
  double Vddsat(0.0);
  //
  double evbs(0.0);
  double evbdd(0.0);
  double sarg(0.0);
  double sargsw(0.0);
  double arg(0.0);
  int Check(1);

  double capgs_old(0.0);
  double capgdd_old(0.0);
  double capgb_old(0.0);

  // diode #1
  double D1temp(0.0);
  double D1power(0.0);
  double D1arg(0.0);
  double D1cd(0.0);
  double D1cdr(0.0);
  double D1csat(0.0);
  double D1csatr(0.0);
  double D1czero(0.0);
  double D1czof2(0.0);
  double D1evd(0.0);
  double D1evr(0.0);
  double D1isr(0.0);
  double D1evrev(0.0);
  double D1gdr(0.0);
  double D1sarg(0.0);
  double D1vdtemp(0.0);
  double D1vtr(0.0);
  int    D1Check(0);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() <<"  Instance::updateIntermediateVars.\n"<<std::endl;
    Xyce::dout() <<"  name = " << getName() << std::endl;
    Xyce::dout() <<"  Model name = " << model_.getName() << std::endl;
    Xyce::dout() <<"  dtype is " << model_.dtype << std::endl;
    Xyce::dout() << std::endl;
    Xyce::dout().width(25); Xyce::dout().precision(17); Xyce::dout().setf(std::ios::scientific);
  }

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

  //  we need our solution variables for any of this stuff

  Vd  = 0.0;
  Vs  = 0.0;
  Vg  = 0.0;
  Vb  = 0.0;
  Vdp = 0.0;
  Vgp = 0.0;
  Vsp = 0.0;
  Vdd = 0.0;

  Vd  = solVec[li_Drain];
  Vg  = solVec[li_Gate];
  Vs  = solVec[li_Source];
  Vb  = solVec[li_Bulk];
  Vsp = solVec[li_SourcePrime];
  Vgp = solVec[li_GatePrime];
  Vdp = solVec[li_DrainPrime];
  Vdd = solVec[li_DrainDrift];

  // node for diode #1
  Vd1p = solVec[li_D1Prime];

  D1vt  = CONSTKoverQ * D1DIOtemp;
  D1vte = model_.D1DIOemissionCoeff * D1vt;
  D1vtr = model_.D1DIOnr * D1vt;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "   "  << std::endl;
    Xyce::dout() << " Solution vector: " << std::endl;
    Xyce::dout() << " Vg  = " << Vg << std::endl;
    Xyce::dout() << " Vd  = " << Vd << std::endl;
    Xyce::dout() << " Vs  = " << Vs << std::endl;
    Xyce::dout() << " Vb  = " << Vb << std::endl;
    Xyce::dout() << " Vdp = " << Vdp << std::endl;
    Xyce::dout() << " Vgp = " << Vgp << std::endl;
    Xyce::dout() << " Vsp = " << Vsp << std::endl;
    Xyce::dout() << " Vdd = " << Vdd << std::endl;
    Xyce::dout() << " Vd1p= " << Vd1p << std::endl;
  }

  // voltage drops
  Vddp  = Vd  - Vdp;
  Vddd  = Vd  - Vdd;
  Vdddp = Vdd - Vdp;
  Vssp  = Vs  - Vsp;
  Vbsp  = Vb  - Vsp;
  Vbdp  = Vb  - Vdp;

  Vggp  = Vg  - Vgp;
  Vgpsp = Vgp - Vsp;
  Vgpdp = Vgp - Vdp;
  Vgpb  = Vgp - Vb;
  Vdpsp = Vdp - Vsp;

  // setup corrections for dtype.
  vbs  = model_.dtype * Vbsp;
  vgps = model_.dtype * Vgpsp;
  vdds = model_.dtype * Vdpsp;
  D1vd = model_.dtype * (Vd1p  - Vd);

  vbdd  = vbs-vdds;
  vgpdd = vgps-vdds;

  // set up the orig voltages.
  origFlag  = 1;
  limitedFlag = false;
  vgps_orig = vgps;
  vdds_orig = vdds;
  vbs_orig  = vbs;
  vbdd_orig = vbdd;
  vgpdd_orig = vgpdd;
  D1vd_orig = D1vd;

  if (getSolverState().newtonIter == 0)
  {

    if (getSolverState().initJctFlag_ && getDeviceOptions().voltageLimiterFlag)
    {
      if (getSolverState().inputOPFlag)
      {
        Linear::Vector * flagSolVectorPtr = extData.flagSolVectorPtr;
        if ((*flagSolVectorPtr)[li_Drain] == 0 || (*flagSolVectorPtr)[li_Gate] == 0 ||
            (*flagSolVectorPtr)[li_Source] == 0 || (*flagSolVectorPtr)[li_SourcePrime] ||
            (*flagSolVectorPtr)[li_DrainDrift] == 0 || (*flagSolVectorPtr)[li_GatePrime] ||
            (*flagSolVectorPtr)[li_DrainPrime] || (*flagSolVectorPtr)[li_Bulk] )
        {
          vbs  = -1;
          vgps = model_.dtype*tVto;
          vdds = 0;
          vbdd = vbs-vdds;
          vgpdd = vgps-vdds;
        }
      }
      else
      {
        vbs  = -1;
        vgps = model_.dtype*tVto;
        vdds = 0;
        vbdd = vbs-vdds;
        vgpdd = vgps-vdds;
      }
    }

    ////////////////////////////////////////////////////////////////////////
    // Note:  the "old" variables should be the values for the previous
    //          Newton iteration.   That previous newton iteration could
    //          have happened in a previous time step...
    if (!(getSolverState().dcopFlag)||(getSolverState().locaEnabledFlag && getSolverState().dcopFlag))
    {
      vbs_old  = (*extData.currStaVectorPtr)[li_state_vbs];
      vbdd_old = (*extData.currStaVectorPtr)[li_state_vbdd];
      vgps_old = (*extData.currStaVectorPtr)[li_state_vgps];
      vdds_old = (*extData.currStaVectorPtr)[li_state_vdds];
      D1vd_old = (*extData.currStaVectorPtr)[li_state_D1vd];
      Von = model_.dtype *
                (*extData.currStaVectorPtr)[li_state_von];
    }
    else
    { // there's no history
      vbs_old = vbs;
      vbdd_old = vbdd;
      vgps_old = vgps;
      vdds_old = vdds;
      D1vd_old = D1vd;
    }
  }
  else
  {
    vbs_old  = (*extData.nextStaVectorPtr)[li_state_vbs];
    vbdd_old = (*extData.nextStaVectorPtr)[li_state_vbdd];
    vgps_old = (*extData.nextStaVectorPtr)[li_state_vgps];
    vdds_old = (*extData.nextStaVectorPtr)[li_state_vdds];
    D1vd_old = (*extData.nextStaVectorPtr)[li_state_D1vd];
    Von = model_.dtype *
              (*extData.nextStaVectorPtr)[li_state_von];
    vgpdd_old = vgps_old-vdds_old;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "After Von first set, ";
    Xyce::dout() << " von = " << von << " Von = " << Von << std::endl;
  }

  ////////////////////////////////////////////
  // SPICE-type Voltage Limiting (PINNING)
  ////////////////////////////////////////////
  if (getDeviceOptions().voltageLimiterFlag)
  {
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << " checking whether to limit voltages  "<< std::endl;
      Xyce::dout() << "  Von = " << Von << std::endl;
      Xyce::dout() << " before limiting: " << std::endl;
      Xyce::dout() << " vgpdd = " << vgpdd <<  " vgpdd_old = " << vgpdd_old << std::endl;
      Xyce::dout() << " vgps  = " << vgps <<  " vgps_old = " << vgps_old << std::endl;
      Xyce::dout() << " vdds = " << vdds <<  " vdds_old = " << vdds_old << std::endl;
      Xyce::dout() << " vbs  = " << vbs <<  " vbs_old = " << vbs_old << std::endl;
      Xyce::dout() << " vbdd = " << vbdd <<  " vbdd_old = " << vbdd_old << std::endl;
    }

    if (vdds_old >= 0)
    {
      vgps = devSupport.fetlim( vgps, vgps_old, Von);
      vdds = vgps - vgpdd;
      vdds = devSupport.limvds( vdds,  vdds_old);
      vgpdd = vgps - vdds;
    }
    else
    {
      vgpdd = devSupport.fetlim( vgpdd, vgpdd_old, Von);
      vdds = vgps - vgpdd;
      vdds = -devSupport.limvds( -vdds, -vdds_old );
      vgps = vgpdd + vdds;
    }

    if (vdds >= 0.0)
    {
      vbs  = devSupport.pnjlim( vbs, vbs_old, vt, sourceVcrit, &Check);
      vbdd = vbs - vdds;
    }
    else
    {
      vbdd = devSupport.pnjlim( vbdd, vbdd_old, vt, drainVcrit, &Check);
      vbs  = vbdd + vdds;
    }

    if (Check == 1) limitedFlag=true;
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << " After limiting: " << std::endl;
      Xyce::dout() << " vgpdd = " << vgpdd << std::endl;
      Xyce::dout() << " vgps  = " << vgps << std::endl;
      Xyce::dout() << " vdds = " << vdds << std::endl;
      Xyce::dout() << " vbs  = " << vbs << std::endl;
      Xyce::dout() << " vbdd = " << vbdd << std::endl;
    }

    //   limit diode junction voltage
    bool tmpGiven = false;
    tmpGiven = (model_.D1DIObreakdownVoltageGiven);
    if (tmpGiven && (D1vd < std::min(0.0, 10.0*D1vte-D1DIOtBrkdwnV) ))
    {
        D1vdtemp = -(D1vd + D1DIOtBrkdwnV);
        D1vdtemp = devSupport.pnjlim(D1vdtemp,
                    -(D1vd_old + D1DIOtBrkdwnV), D1vte,
                                 D1DIOtVcrit, &D1Check);
        D1vd = -(D1vdtemp + D1DIOtBrkdwnV);
    }
    else
    {
        D1vd = devSupport.pnjlim(D1vd, D1vd_old, D1vte,
                                                   D1DIOtVcrit, &D1Check);
    }

    if (D1Check == 1) limitedFlag=true;
  } // voltage limiter flag

  ////
  // now all the preliminaries are over - we can start doing the
  //  real work
  ////
  vbdd  = vbs - vdds;
  vgpdd = vgps - vdds;
  Vgpb  = vgps - vbs;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << " vbs  = " << vbs << std::endl;
    Xyce::dout() << " vgps = " << vgps << std::endl;
    Xyce::dout() << " vdds = " << vdds << std::endl;
    Xyce::dout() << " vbdd = " << vbdd << std::endl;
    Xyce::dout() << " vgpdd= " << vgpdd << std::endl;

    Xyce::dout() << " Vddp = " << Vddp << std::endl;
    Xyce::dout() << " Vddd = " << Vddd << std::endl;
    Xyce::dout() << " Vdddp= " << Vdddp << std::endl;
    Xyce::dout() << " Vssp = " << Vssp << std::endl;
    Xyce::dout() << " Vbsp = " << Vbsp << std::endl;
    Xyce::dout() << " Vbdp = " << Vbdp << std::endl;
    Xyce::dout() << " Vggp = " << Vggp << std::endl;
    Xyce::dout() << " Vgpsp = " << Vgpsp << std::endl;
    Xyce::dout() << " Vgpdp = " << Vgpdp << std::endl;
    Xyce::dout() << " Vgpb  = " << Vgpb  << std::endl;
    Xyce::dout() << " Vdpsp= " << Vdpsp << std::endl;

  }

  // Now set the origFlag
  if (vgps_orig != vgps || vdds_orig != vdds ||
      vbs_orig != vbs || vbdd_orig != vbdd || vgpdd_orig != vgpdd) origFlag = 0;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    if (origFlag == 0)
    {
      Xyce::dout() << " Something modified the voltages. " << std::endl;
      Xyce::dout() << " Voltage       before      after        diff " << std::endl;
      Xyce::dout() << " vgps       " << vgps_orig << "  " << vgps << " " << vgps-vgps_orig << std::endl;
      Xyce::dout() << " vdds       " << vdds_orig << "  " << vdds << " " << vdds-vdds_orig << std::endl;
      Xyce::dout() << " vbs       " << vbs_orig << "  " << vbs << " " << vbs-vbs_orig << std::endl;
      Xyce::dout() << " vbdd       " << vbdd_orig << "  " << vbdd << " " << vbdd-vbdd_orig << std::endl;
      Xyce::dout() << " vgpdd       " << vgpdd_orig << "  " << vgpdd << " " << vgpdd-vgpdd_orig << std::endl;
    }
  }


  ////
  //   bulk-source and bulk-drain diodes
  //   here we just evaluate the ideal diode current and the
  //   corresponding derivative (conductance).
  ////

  double pi(0.0);
  double end10(0.0);
  double end1(0.0);
  double end20(0.0);
  double end2(0.0);
  double end3(0.0);
  double j1(0.0);
  double fr1(0.0);
  double fr2(0.0);
  double d1j1(0.0);
  double fr3(0.0);
  double d2j1(0.0);
  double dj1(0.0);

  if(model_.artd != 0.0 && sourceArea != 0.0)
  {
    end10=model_.brtd-model_.crtd+model_.nrtd*vbs;
    end1=1+exp(std::min(CONSTMAX_EXP_ARG,end10/vt));
    end20=model_.brtd-model_.crtd-model_.nrtd*vbs;
    end2=1+exp(std::min(CONSTMAX_EXP_ARG,end20/vt));
    end3=(model_.crtd-model_.nrtd*vbs)/model_.drtd;
    pi = 3.1415927;
    j1=model_.artd*log(end1/end2)*(pi/2+atan(end3));
    j1=j1*sourceArea;

    fr1=(end1-1)*model_.nrtd/(end1*vt);
    fr2=(end2-1)*model_.nrtd/(end2*vt);
    d1j1=(fr1-fr2)*(pi/2+atan(end3));
    fr3=(-model_.nrtd/model_.drtd)/(1+end3*end3);
    d2j1=log(end1/end2)*fr3;
    dj1=model_.artd*(d1j1+d2j1);
    dj1=dj1*sourceArea;
  }
  else
  {
    j1  = 0.0;
    dj1 = 0.0;
  }

  if(vbs <= 0)
  {
    gbs = numberParallel*(SourceSatCur*model_.n2/vt+dj1);
    gbs += numberParallel*getDeviceOptions().gmin;
    cbs = gbs*vbs+j1;
// PMC  for D1 diode testing
//    gbs = 0.0;
//    cbs = 0.0;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "*******Setting cbs for vbs<=0 ******" << std::endl;
      Xyce::dout() << "         vbs  = " << vbs << std::endl;
      Xyce::dout() << "          vt  = " << vt  << std::endl;
      Xyce::dout() << "          n2  = " << model_.n2  << std::endl;
      Xyce::dout() << "         SSC  = " << SourceSatCur << std::endl;
      Xyce::dout() << "         gbs  = " << gbs  << std::endl;
      Xyce::dout() << "         cbs  = " << cbs << std::endl;
      Xyce::dout() << "         j1   = " << j1 << std::endl;
    }
  }
  else
  {
    evbs = exp(std::min(CONSTMAX_EXP_ARG,model_.n2*vbs/vt));
    gbs  = numberParallel*(SourceSatCur*model_.n2*evbs/vt
                           + getDeviceOptions().gmin + dj1);
    cbs  = numberParallel*(SourceSatCur*(evbs-1) + j1);
// PMC  for D1 diode testing
//    gbs = 0.0;
//    cbs = 0.0;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "*******Setting cbs for vbs>0 ******" << std::endl;
      Xyce::dout() << "     vbs  = " << vbs << std::endl;
      Xyce::dout() << "      vt  = " << vt  << std::endl;
      Xyce::dout() << "      n2  = " << model_.n2 << std::endl;
      Xyce::dout() << "  vbs/vt  = " << vbs/vt  << std::endl;
      Xyce::dout() << "  Maxarg  = " << CONSTMAX_EXP_ARG << std::endl;
      Xyce::dout() << "     arg  = " << std::min(CONSTMAX_EXP_ARG,
                                         model_.n2*vbs/vt) << std::endl;
      Xyce::dout() << "    evbs  = " << evbs << std::endl;
      Xyce::dout() << "     gbs  = " << gbs << std::endl;
      Xyce::dout() << "     cbs  = " << cbs << std::endl;
    }
  }
  if(vbdd <= 0)
  {
    gbd = numberParallel*(DrainSatCur*model_.n2/vt);
    gbd += numberParallel*getDeviceOptions().gmin;
    cbd = gbd *vbdd;
// PMC  for D1 diode testing
//    gbd = 0.0;
//    cbd = 0.0;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "*******Setting cbd for vbdd<=0 ******" << std::endl;
      Xyce::dout() << "         vbdd = " << vbdd << std::endl;
      Xyce::dout() << "         SSC  = " << SourceSatCur << std::endl;
      Xyce::dout() << "          vt  = " << vt  << std::endl;
      Xyce::dout() << "          n2  = " << model_.n2 << std::endl;
      Xyce::dout() << "         gbd  = " << gbd << std::endl;
      Xyce::dout() << "         cbd  = " << cbd << std::endl;
    }

  }
  else
  {
    evbdd = exp(std::min(CONSTMAX_EXP_ARG,model_.n2*vbdd/vt));
    gbd   = numberParallel*(DrainSatCur*model_.n2*evbdd/vt
                            + getDeviceOptions().gmin);
    cbd   = numberParallel*DrainSatCur*(evbdd-1);
// PMC  for D1 diode testing
//    gbd = 0.0;
//    cbd = 0.0;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "*******Setting cbd for vbdd>0 ******" << std::endl;
      Xyce::dout() << "     vbdd  = " << vbdd << std::endl;
      Xyce::dout() << "      vt   = " << vt  << std::endl;
      Xyce::dout() << "  vbdd/vt  = " << vbdd/vt  << std::endl;
      Xyce::dout() << "  Maxarg   = " << CONSTMAX_EXP_ARG << std::endl;
      Xyce::dout() << "     arg   = " << std::min(CONSTMAX_EXP_ARG,
                                        model_.n2*vbdd/vt) << std::endl;
      Xyce::dout() << "    evbdd  = " << evbdd << std::endl;
      Xyce::dout() << "     gbd   = " << gbd << std::endl;
      Xyce::dout() << "     cbd   = " << cbd << std::endl;
    }
  }

  if (vdds >= 0)
    mode = 1;
  else
    mode = -1;

  double dvonvbs(0.0);
  double lvbs = mode == 1 ? vbs : vbdd;
  UCCMcvon(lvbs, von, dvonvbs);
  Von   = model_.dtype * von;   // von contains dtype, Von does not

  if(!ISUBMOD)
  {
    // Old version of model with corrected Jacobian terms
    UCCMmosa1(vdds>0?vgps:vgpdd,mode*vdds,dvonvbs, cdraindrift, vddsat);
  }
  else
  {
    // Version of model with substrate currents
    UCCMmosa2(vdds>0?vgps:vgpdd,mode*vdds,dvonvbs, cdraindrift, vddsat);
  }
  Vddsat = model_.dtype * vddsat;

  ////
  //     COMPUTE EQUIVALENT DRAIN CURRENT SOURCE
  ////

    // Substrate current calculations
  ISUB            = cdraindrift*mm1;
  GMSUB           = gm*mm1   + cdraindrift*dmm1vgs;
  GDDSSUB         = gdds*mm1 + cdraindrift*dmm1vds;
  GBSSUB          = gmbs*mm1 + cdraindrift*dmm1vbs;

// PMC  ISUBMOD not fully implemented
//                 - probably we don't have the right model parameters
//  cdrain = cdraindrift + ISUB;

  cdrain = cdraindrift;
  cd = mode*cdrain - cbd;

  // diode #1 values
  D1csat  = D1DIOtSatCur * D1DIOarea;
  D1csatr = D1DIOtSatRCur * D1DIOarea;
  D1gspr  = model_.D1DIOconductance * D1DIOarea;
  D1vtr   = model_.D1DIOnr * D1vt;
  D1evd   = D1arg = D1evrev = 0.0;

  //   compute diode dc current and derivatives
  if(D1csatr != 0)
  {
    D1evr   = exp(D1vd/D1vtr);
    D1temp  = 1 - D1vd/D1DIOtJctPot;
    D1arg   = D1temp*D1temp;
    D1power = pow(D1arg + 0.001,model_.D1DIOgradingCoeff/2);
    D1isr   = D1csatr*D1power;
    D1cdr   = D1isr*(D1evr - 1);
    D1gdr   = -D1cdr*model_.D1DIOgradingCoeff*D1temp/
                         (D1DIOtJctPot*(D1arg + 0.001)) + D1isr*D1evr/D1vtr;
  } else
  {
    D1cdr = 0;
    D1gdr = 0;
  }

  // add in diode bulk resistance
  if(D1vd >= -3*D1vte)
  {
    D1evd = exp(D1vd/D1vte);
    D1cd = D1csat*(D1evd-1)+D1cdr + getDeviceOptions().gmin*D1vd;
    D1gd = D1csat*D1evd/D1vte+D1gdr + getDeviceOptions().gmin;
    if(model_.D1DIOikf > 0) {
      D1arg = sqrt(model_.D1DIOikf/(model_.D1DIOikf + D1cd));
      D1gd = D1arg*D1gd*(1 - 0.5*D1cd/(model_.D1DIOikf + D1cd));
      D1cd = D1arg*D1cd;
    }
  }
  else if( (!(D1DIOtBrkdwnV)) || (D1vd >= -D1DIOtBrkdwnV) )
  {
    D1arg = 3*D1vte/(D1vd*CONSTe);
    D1arg = D1arg * D1arg * D1arg;
    D1cd  = -D1csat*(1 + D1arg) + getDeviceOptions().gmin*D1vd + D1cdr;
    D1gd  = D1csat*3*D1arg/D1vd + getDeviceOptions().gmin + D1gdr;
  }
  else
  {
    D1evrev= exp(-(D1DIOtBrkdwnV + D1vd)/D1vte);
    D1cd   = -D1csat*D1evrev + getDeviceOptions().gmin*D1vd + D1cdr;
    D1gd   = D1csat*D1evrev/D1vte + getDeviceOptions().gmin + D1gdr;
  }
  D1cdeq = D1cd;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "  "  << std::endl;
    Xyce::dout() << " dtype   = " << model_.dtype << std::endl;
    Xyce::dout() << " D1vd    = " << D1vd << std::endl;
    Xyce::dout() << " D1vtr   = " << D1vtr << std::endl;
    Xyce::dout() << " D1csat  = " << D1csat << std::endl;
    Xyce::dout() << " D1csatr = " << D1csatr << std::endl;
    Xyce::dout() << " D1gspr  = " << D1gspr << std::endl;
    Xyce::dout() << " D1cdr   = " << D1cdr << std::endl;
    Xyce::dout() << " D1gdr   = " << D1gdr << std::endl;
    Xyce::dout() << " D1vte   = " << D1vte << std::endl;
    Xyce::dout() << " D1evd   = " << D1evd << std::endl;
    Xyce::dout() << " D1arg   = " << D1arg << std::endl;
    Xyce::dout() << " D1evrev = " << D1evrev << std::endl;
    Xyce::dout() << " D1DIOtBrkdwnV = " << D1DIOtBrkdwnV << std::endl << std::endl;
    Xyce::dout() << " D1cd    = " << D1cd << std::endl;
    Xyce::dout() << " D1gd    = " << D1gd << std::endl;
    Xyce::dout() << " D1cdeq  = " << D1cdeq << std::endl;
  }

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
      capbs=Cbs*sarg + Cbssw*sargsw;
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

    if (vbdd < tDepCap)
    {
      arg=1-vbdd/tBulkPot;
      ////
      // * the following block looks somewhat long and messy,
      // * but since most users use the default grading
      // * coefficients of .5, and sqrt is MUCH faster than an
      // * exp(log()) we use this special case code to buy time.
      // * (as much as 10% of total job time!)
      ////
      if(model_.bulkJctBotGradingCoeff  == .5 &&
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
        tBulkPot*(Cbd*(1-arg*sarg)
                  /(1-model_.bulkJctBotGradingCoeff)
                  +Cbdsw*(1-arg*sargsw)
                  /(1-model_.bulkJctSideGradingCoeff));
      capbd=Cbd*sarg + Cbdsw*sargsw;
    }
    else
    {
      qbd  = f4d + vbdd * (f2d + vbdd * f3d/2);
      capbd= f2d + vbdd * f3d;
    }
  }
  else
  {
    qbd = 0;
    capbd = 0;
  }

  //   charge storage elements
  D1czero = D1DIOtJctCap*D1DIOarea;
  if (D1vd < D1DIOtDepCap)
  {
    D1arg = 1-D1vd/model_.D1DIOjunctionPot;
    D1sarg = exp(-model_.D1DIOgradingCoeff*log(D1arg));
    D1DIOcapCharge = model_.D1DIOtransitTime*D1cd +
                  model_.D1DIOjunctionPot*D1czero*
                  (1-D1arg*D1sarg)/(1-model_.D1DIOgradingCoeff);
    D1capd = model_.D1DIOtransitTime*D1gd + D1czero*D1sarg;
  }
  else
  {
    D1czof2 = D1czero/model_.D1DIOf2;
    D1DIOcapCharge = model_.D1DIOtransitTime*D1cd+D1czero*
                        D1DIOtF1+D1czof2*(model_.D1DIOf3*
                     (D1vd-D1DIOtDepCap) + (model_.D1DIOgradingCoeff/
                        (model_.D1DIOjunctionPot +
                        model_.D1DIOjunctionPot))*
                        (D1vd*D1vd-D1DIOtDepCap*D1DIOtDepCap));
    D1capd = model_.D1DIOtransitTime*D1gd+D1czof2*
             (model_.D1DIOf3 + model_.D1DIOgradingCoeff*D1vd/
                                 model_.D1DIOjunctionPot);
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "  "  << std::endl;
    Xyce::dout() << " Going into qmeyer..." << std::endl;
    Xyce::dout() << " Mode is " << mode << std::endl;
    Xyce::dout() << " Args are vgps = " << vgps << " vgpdd = " << vgpdd << std::endl;
    Xyce::dout() << " Vgpb = " << Vgpb << " Von = " << Von << " Vddsat = " << Vddsat << std::endl;
    Xyce::dout() << " tPhi = " << tPhi << " OxideCap = " << OxideCap << std::endl;
  }
  if (mode > 0)
  {
    if(model_.cve == 1)
    {
      UCCMqmeyer (vgps,vgpdd,Vgpb,Von,Vddsat,
                       capgs, capgdd, capgb, tPhi,OxideCap);
    }
    else
    {
      // Ward-like model
      UCCMMeyercap(vgps,vgpdd,Vgpb, capgs, capgdd, capgb);
    }
  }
  else
  {
    if(model_.cve == 1)
    {
      UCCMqmeyer (vgpdd,vgps,Vgpb,Von,Vddsat,
                       capgdd, capgs, capgb, tPhi,OxideCap);
    }
    else
    {
      // Ward-like model
      UCCMMeyercap(vgpdd,vgps,Vgpb, capgs, capgdd, capgb);
    }
  }

  if((getSolverState().dcopFlag))
  {
    Capgs  =  2 * capgs  + GateSourceOverlapCap ;
    Capgdd =  2 * capgdd + GateDrainOverlapCap ;
    Capgb  =  2 * capgb  + GateBulkOverlapCap ;
  }
  else
  {
    capgs_old  = (*extData.currStaVectorPtr)[li_state_capgs];
    capgdd_old = (*extData.currStaVectorPtr)[li_state_capgdd];
    capgb_old  = (*extData.currStaVectorPtr)[li_state_capgb];

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "Doing meyer back averaging..."<< std::endl;
      Xyce::dout() << " capgs  = " << capgs  << " capgs_old  = " << capgs_old << std::endl;
      Xyce::dout() << " capgdd = " << capgdd << " capgdd_old = " << capgdd_old << std::endl;
      Xyce::dout() << " capgb  = " << capgb  << " capgb_old  = " << capgb_old << std::endl;
    }
    Capgs  = ( capgs+capgs_old   + GateSourceOverlapCap );
    Capgdd = ( capgdd+capgdd_old + GateDrainOverlapCap );
    Capgb  = ( capgb+capgb_old   + GateBulkOverlapCap );
  }
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Capgs  = " << Capgs << std::endl;
    Xyce::dout() << "Capgdd = " << Capgdd << std::endl;
    Xyce::dout() << "Capgb  = " << Capgb << std::endl;
    Xyce::dout() << "capgs  = " << capgs << std::endl;
    Xyce::dout() << "capgdd = " << capgdd << std::endl;
    Xyce::dout() << "capgb  = " << capgb << std::endl;
  }
  Capgs  *= (Capgs  < 0.0)?-1.0:1.0;
  Capgdd *= (Capgdd < 0.0)?-1.0:1.0;
  Capgb  *= (Capgb  < 0.0)?-1.0:1.0;

  // Voltage-dependent drain-drift conductance
  double absV = fabs(Vddd);
  gddd = (model_.driftParamA + model_.driftParamB*absV);
  if(gddd != 0) gddd = 1.0/gddd;
  draindriftCond = gddd;

  // need a more precise derivative than gddd for the jacobian.
  // w.r.t. Vd and Vdd.
  double dVddd_dVd (1.0);
  double d_absVddd_dVd (0.0);
  if (Vddd > 0.0)
  {
    d_absVddd_dVd = 1.0;
  }
  else if (Vddd < 0.0)
  {
    d_absVddd_dVd = -1.0;
  }

  dIdd_dVd = dVddd_dVd * gddd -
              Vddd * (gddd*gddd)*(model_.driftParamB)*d_absVddd_dVd;

  // now do the shundt resistor
  if(model_.rdsshunt == 0.0)
  {
     gdsshunt = 0.0;
  }
  else
  {
    gdsshunt = 1.0/model_.rdsshunt + getDeviceOptions().gmin;
  }

  Idrain  = drainCond * Vdddp;
  Igate   = gateCond * Vggp;
  Isource = sourceCond * Vssp;
  Idraindrift = draindriftCond * Vddd;
  Irdsshunt = gdsshunt * (Vd - Vs);
  // Rd1rs current:
  Ird1rs = D1gspr*(Vs - Vd1p);

  if (mode >= 0)   // Normal mode
  {
    Gm     = gm;      // (xnrm-xrev)*gm  in 3f5
    Gmbs   = gmbs;    // (xnrm-xrev)*gmbs in 3f5
    nrmsum = Gm+Gmbs; // xnrm*(gm+gmbs)
    revsum = 0;       // xrev*(gm+gmbs)
    cdreq  = model_.dtype*cdrain;
  }
  else
  {
    Gm     = -gm;
    Gmbs   = -gmbs;
    nrmsum = 0;
    revsum = -(Gm+Gmbs);  // because Gm and Gmbs already have - in them!
    cdreq  = -(model_.dtype)*cdrain;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << " Done with Instance::updateIntermediateVars " << std::endl;
    Xyce::dout() << "  mode         = " << mode << std::endl;
    Xyce::dout() << "  Idrain       = " << Idrain << std::endl;
    Xyce::dout() << "  Igate        = " << Igate << std::endl;
    Xyce::dout() << "  Isource      = " << Isource << std::endl;
    Xyce::dout() << "  Idraindrift  = " << Idraindrift << std::endl;
    Xyce::dout() << "  Ird1rs  = " << Ird1rs << std::endl;
    Xyce::dout() << "  dIdd_dVd     = " << dIdd_dVd << std::endl;
    Xyce::dout() << "  gddd         = " << gddd << std::endl;
    Xyce::dout() << "  Irdsshunt    = " << Irdsshunt << std::endl;
    Xyce::dout() << "  D1DIOcapCurrent = " << D1DIOcapCurrent << std::endl;
    Xyce::dout() << "  cbd          = " << cbd << std::endl;
    Xyce::dout() << "  cbs          = " << cbs  << std::endl;
    Xyce::dout() << "  qbd          = " << qbd  << std::endl;
    Xyce::dout() << "  qbs          = " << qbs << std::endl;
    Xyce::dout() << "  cdrain       = " << cdrain  << std::endl;
    Xyce::dout() << "  cdraindrift  = " << cdraindrift  << std::endl;
    Xyce::dout() << "  cdreq        = " << cdreq  << std::endl;
    Xyce::dout() << "  gdds         = " << gdds << std::endl;
    Xyce::dout() << "  gdsshunt     = " << gdsshunt << std::endl;
    Xyce::dout() << "  gm           = " << gm << std::endl;
    Xyce::dout() << "  gmbs         = " << gmbs << std::endl;
    Xyce::dout() << "  Gm           = " << Gm << std::endl;
    Xyce::dout() << "  Gmbs         = " << Gmbs << std::endl;
  }

  /// CURRENTS to load into RHS:

  // current out of drain is
  // Idraindrift + Irdsshunt - Dtype*(D1cdeq + D1DIOcapCurrent)

  // current out of gate:
  // dtype*( (deriv of qgs) + (deriv of qgdd) + (deriv of qgb))

  //  the current *out of* the source is
  // Isource - Irdsshunt + Ird1rs

  // current out of bulk is
  // dtype*(deriv of qbd) + dtype*cbd + dtype*cbs + dtype*(deriv of qbs)
  //  - dtype*(deriv of qgb)

  // current out of drain' is
  // -Idrain - dtype*(deriv of qgd) - (deriv of qbd) - dtype*cbd +
  //  mode*dtype*cdrain

  // the current out of the source' is
  //  -Isource - dtype*(deriv of qgs) - dtype*cbs - (deriv of qbs) -
  //   mode*dtype*cdrain -Irdsshunt

  // the current out of the draindrift node is
  // Idrain - IdraindriftDtype*(D1cdeq + D1DIOcapCurrent)

 // the current out of the d1' pos node is
  // Dtype*(D1cdeq + D1DIOcapCurrent) - Ird1rs

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : pmc
// Creation Date : 2/17/2004
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp_tmp)
{
  bool bsuccess = true;
  // mos3temp vars
  double czbd(0.0);    // zero voltage bulk-drain capacitance
  double czbdsw(0.0);  // zero voltage bulk-drain sidewall capacitance
  double czbs(0.0);    // zero voltage bulk-source capacitance
  double czbssw(0.0);  // zero voltage bulk-source sidewall capacitance
  double arg(0.0);     // 1 - fc
  double sarg(0.0);    // (1-fc) ^^ (-mj)
  double sargsw(0.0);  // (1-fc) ^^ (-mjsw)
  double ratio,ratio4(0.0);
  double fact2(0.0);
  double kt(0.0);
  double egfet(0.0);
  double pbfact(0.0);
  double capfact(0.0);
  double phio(0.0);
  double pbo(0.0);
  double gmanew,gmaold(0.0);
  // end of mos3temp stuff

//variables related to the diode 1
  double D1xfc(0.0);
  double D1vte_loc(0.0);
  double D1cbv(0.0);
  double D1xbv(0.0);
  double D1xcbv(0.0);
  double D1tol(0.0);
  double D1vt_loc(0.0);
  int D1iter(0);
  double D1egfet1(0.0),D1arg1(0.0),D1fact1(0.0),D1pbfact1(0.0),D1pbo(0.0),D1gmaold(0.0);
  double D1fact2(0.0),D1pbfact(0.0),D1arg(0.0),D1egfet(0.0),D1gmanew(0.0);
  double reltol(0.0);

  double tnom(0.0);
  double VMAX(0.0);

// oxide dielectric permitivity
//#define EPSSIO2    3.453130e-11 == CONSTEPSOX (Xyce value)

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
    // make sure interpolation doesn't take any resistance negative
    if(model_.sheetResistance  < 0) model_.sheetResistance = 0;
    if(model_.drainResistance  < 0) model_.drainResistance = 0;
    if(model_.gateResistance   < 0) model_.gateResistance  = 0;
    if(model_.sourceResistance < 0) model_.sourceResistance= 0;
    if(model_.D1DIOresist      < 0) model_.D1DIOresist = 0;

    // some params may have changed during interpolation
    model_.processParams();
  }

  VMAX = model_.maxDriftVel;
  tnom = model_.tnom;
  ratio = temp/tnom;
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << " Temperature = "<< temp << std::endl;
    Xyce::dout() << " tnom = " << tnom << std::endl;
    Xyce::dout() << " ratio = " << ratio << std::endl;
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
    Xyce::dout() << " fact1  = " << model_.fact1 << std::endl;
    Xyce::dout() << " vtnom  = " << model_.vtnom << std::endl;
    Xyce::dout() << " egfet1 = " << model_.egfet1 << std::endl;
    Xyce::dout() << " pbfact1= " << model_.pbfact1 << std::endl;
    Xyce::dout() << " vt     = " << vt << std::endl;
    Xyce::dout() << " ratio  = " << ratio << std::endl;
    Xyce::dout() << " fact2  = " << fact2 << std::endl;
    Xyce::dout() << " kt     = " << kt << std::endl;
    Xyce::dout() << " egfet  = " << egfet << std::endl;
    Xyce::dout() << " arg    = " << arg << std::endl;
    Xyce::dout() << " pbfact = " << pbfact << std::endl;
  }

  ratio4 = ratio * sqrt(ratio);
  tSurfMob = 1.e-4*model_.surfaceMobility/ratio4;
  phio= (model_.phi-model_.pbfact1)/model_.fact1;
  tPhi = fact2 * phio + pbfact;
  tVbi = model_.vt0 - model_.dtype *
    (model_.gamma* sqrt(model_.phi))+.5*(model_.egfet1-egfet)
    + model_.dtype*.5* (tPhi-model_.phi);
  tVto = tVbi + model_.dtype * model_.gamma * sqrt(tPhi);
  if(model_.mdtemp != 0) tVto = model_.vt0 -
                                       model_.kvt*(temp-tnom);
  tSatCur = model_.jctSatCur *
                   exp(-egfet/vt+model_.egfet1/model_.vtnom);
  tSatCurDens = model_.jctSatCurDensity *
                   exp(-egfet/vt+model_.egfet1/model_.vtnom);
  pbo = (model_.bulkJctPotential-model_.pbfact1)/model_.fact1;
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
  capfact = (1+model_.bulkJctBotGradingCoeff *
                                           (4e-4*(temp-CONSTREFTEMP)-gmanew));
  tCbd *= capfact;
  tCbs *= capfact;
  tCj *= capfact;
  capfact = (1+model_.bulkJctSideGradingCoeff *
                                            (4e-4*(temp-CONSTREFTEMP)-gmanew));
  tCjsw *= capfact;
  tDepCap = model_.fwdCapDepCoeff * tBulkPot;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << " ratio4 = " << ratio4 << std::endl;
    Xyce::dout() << " tSurfMob = " << tSurfMob << std::endl;
    Xyce::dout() << " phio = " << phio << std::endl;
    Xyce::dout() << " tPhi = " << tPhi << std::endl;
    Xyce::dout() << " tVbi = " << tVbi << std::endl;
    Xyce::dout() << " tVto = " << tVto << std::endl;
    Xyce::dout() << " tSatCur = " << tSatCur << std::endl;
    Xyce::dout() << " tSatCurDens = " << tSatCurDens << std::endl;
    Xyce::dout() << " pbo = " << pbo << std::endl;
    Xyce::dout() << " gmaold = " << gmaold << std::endl;
    Xyce::dout() << " tBulkPot = " << tBulkPot << std::endl;
    Xyce::dout() << " gmanew = " << gmanew << std::endl;
    Xyce::dout() << " capfact = " << capfact << std::endl;
    Xyce::dout() << " tCbd = " << tCbd << std::endl;
    Xyce::dout() << " tCbs = " << tCbs << std::endl;
    Xyce::dout() << " tCj = " << tCj << std::endl;
    Xyce::dout() << " capfact = " << capfact << std::endl;
    Xyce::dout() << " tCjsw = " << tCjsw << std::endl;
    Xyce::dout() << " tDepCap = " << tDepCap << std::endl;
  }

  if( (model_.jctSatCurDensity == 0) || (drainArea == 0) ||
      (sourceArea == 0) )
  {
    sourceVcrit = drainVcrit =
      vt*log(vt/(CONSTroot2*model_.jctSatCur));
  }
  else
  {
    drainVcrit  = vt * log( vt / (CONSTroot2 *
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
  f2d = czbd*(1-model_.fwdCapDepCoeff*
              (1+model_.bulkJctBotGradingCoeff))* sarg/arg
               + czbdsw*(1-model_.fwdCapDepCoeff*
               (1+model_.bulkJctSideGradingCoeff))*sargsw/arg;
  f3d = czbd * model_.bulkJctBotGradingCoeff * sarg/arg/tBulkPot
               + czbdsw * model_.bulkJctSideGradingCoeff *
                sargsw/arg/tBulkPot;
  f4d = czbd*tBulkPot*(1-arg*sarg)/(1-model_.bulkJctBotGradingCoeff) +
                czbdsw*tBulkPot*(1-arg*sargsw)/
                (1-model_.bulkJctSideGradingCoeff) -
                 f3d/2*(tDepCap*tDepCap) - tDepCap * f2d;
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
  f2s = czbs*(1-model_.fwdCapDepCoeff*
              (1+model_.bulkJctBotGradingCoeff))*sarg/arg +
               czbssw*(1-model_.fwdCapDepCoeff*
               (1+model_.bulkJctSideGradingCoeff))*sargsw/arg;
  f3s = czbs * model_.bulkJctBotGradingCoeff * sarg/arg/tBulkPot
               + czbssw * model_.bulkJctSideGradingCoeff *
               sargsw/arg /tBulkPot;
  f4s = czbs*tBulkPot*(1-arg*sarg)/(1-model_.bulkJctBotGradingCoeff)
               + czbssw*tBulkPot*(1-arg*sargsw)/
               (1-model_.bulkJctSideGradingCoeff)-f3s/2*
               (tDepCap*tDepCap) - tDepCap * f2s;

  n0 = CONSTEPSOX*model_.eta*vt/2/CONSTQ/model_.oxideThickness;

  if(!model_.vpGiven)
  {
    vp = VMAX*(l-2*model_.latDiff)/tSurfMob;
    if(model_.dtype == CONSTNMOS) vp *= 2.0;
  } else
    vp = model_.vp;

  gchi0 = CONSTQ*w/(l-2*model_.latDiff);
  gammas =
    model_.gammas0+model_.lgammas*(1-model_.l0/l)+
    model_.wgammas*(1-model_.w0/w);
  gammal =
    model_.gammal0+model_.lgammal_*(1-model_.l0/l)+
    model_.wgammal*(1-model_.w0/w);
  if(gammal != 0.0) {
    vthLimit = gammas/(2*gammal);
    vthLimit = vthLimit*vthLimit;
  } else
    vthLimit = Util::MachineDependentParams::DoubleMax();
  vtoo = tVto*model_.dtype+ gammal*tPhi- gammas*sqrt(tPhi);

  // quantities related to diode #1

    D1DIOtemp = temp;
    model_.D1DIOnomTemp = tnom;
    D1vt_loc = CONSTKoverQ * D1DIOtemp;
    D1xfc = log(1-model_.D1DIOdepletionCapCoeff);

// this part gets really ugly - I don't know how to explain these equations

    D1fact2 = D1DIOtemp/CONSTREFTEMP;
    D1egfet = 1.16-(7.02e-4*D1DIOtemp*D1DIOtemp)/(D1DIOtemp+1108);
    D1arg = -D1egfet/(2*CONSTboltz*D1DIOtemp) + 1.1150877/
                                (CONSTboltz*(CONSTREFTEMP+CONSTREFTEMP));
    D1pbfact = -2*D1vt_loc*(1.5*log(D1fact2)+CONSTQ*D1arg);
    D1egfet1 = 1.16 - (7.02e-4*model_.D1DIOnomTemp*
                 model_.D1DIOnomTemp)/(model_.D1DIOnomTemp+1108);
    D1arg1 = -D1egfet1/(CONSTboltz*2*model_.D1DIOnomTemp) +
                                       1.1150877/(2*CONSTboltz*CONSTREFTEMP);
    D1fact1 = model_.D1DIOnomTemp/CONSTREFTEMP;
    D1pbfact1 = -2*CONSTKoverQ*model_.D1DIOnomTemp*
                                          (1.5*log(D1fact1)+CONSTQ*D1arg1);
    D1pbo = (model_.D1DIOjunctionPot-D1pbfact1)/D1fact1;
    D1gmaold = (model_.D1DIOjunctionPot -D1pbo)/D1pbo;
    D1DIOtJctCap = model_.D1DIOjunctionCap/
                   (1+model_.D1DIOgradingCoeff*
                   (400e-6*(model_.D1DIOnomTemp-CONSTREFTEMP)-D1gmaold) );
    D1DIOtJctPot = D1pbfact+D1fact2*D1pbo;
    D1gmanew = (D1DIOtJctPot-D1pbo)/D1pbo;
    D1DIOtJctCap *= 1+model_.D1DIOgradingCoeff*
                    (400e-6*(D1DIOtemp-CONSTREFTEMP)-D1gmanew);
    D1DIOtSatCur = model_.D1DIOsatCur*exp( ((D1DIOtemp/
                   model_.D1DIOnomTemp)-1)*
                   model_.D1DIOactivationEnergy/
                   (model_.D1DIOemissionCoeff*D1vt_loc)+
                   model_.D1DIOsaturationCurrentExp/
                   model_.D1DIOemissionCoeff*
                   log(D1DIOtemp/model_.D1DIOnomTemp) );
    D1DIOtSatRCur = model_.D1DIOisr;

    // the defintion of f1, just recompute after temperature adjusting
    // all the variables used in it
    D1DIOtF1=D1DIOtJctPot*(1-exp((1-model_.D1DIOgradingCoeff)*D1xfc))/
                    (1-model_.D1DIOgradingCoeff);

    // same for Depletion Capacitance
    D1DIOtDepCap=model_.D1DIOdepletionCapCoeff*D1DIOtJctPot;

    // and Vcrit
    D1vte_loc=model_.D1DIOemissionCoeff*D1vt_loc;
    D1DIOtVcrit=D1vte_loc*log(D1vte_loc/(CONSTroot2*D1DIOtSatCur));

    // and now to copute the breakdown voltage, again using
    // temperature adjusted basic parameters
    if (model_.D1DIObreakdownVoltageGiven)
    {
        D1cbv=model_.D1DIObreakdownCurrent;
        if (D1cbv < D1DIOtSatCur*model_.D1DIObreakdownVoltage/D1vt_loc)
        {
          D1cbv=D1DIOtSatCur*model_.D1DIObreakdownVoltage/D1vt_loc;
          Xyce::dout() << " breakdown current increased to " << D1cbv <<
                    "to resolve incompatability " <<
                    "with specified saturation current" << std::endl;
          D1xbv=model_.D1DIObreakdownVoltage;
        } else
        {
          reltol = 1e-3;
          D1tol=reltol*D1cbv;
          D1xbv=model_.D1DIObreakdownVoltage-
                                        D1vt_loc*log(1+D1cbv/D1DIOtSatCur);
          for(D1iter=0; D1iter<25; ++D1iter)
          {
              D1xbv=model_.D1DIObreakdownVoltage-D1vt_loc*log(D1cbv/
                              D1DIOtSatCur+1-D1xbv/D1vt_loc);
              D1xcbv=D1DIOtSatCur*(exp((model_.D1DIObreakdownVoltage
                              -D1xbv)/D1vt_loc)-1+D1xbv/D1vt_loc);
              if (fabs(D1xcbv-D1cbv) <= D1tol) goto matched;
          }
          Xyce::dout() << " unable to match forward and reverse diode regions: D1bv = "
              << D1xbv << " D1ibv = " << D1xcbv << std::endl;
        }
        matched:
        D1DIOtBrkdwnV = D1xbv;
    }

  return bsuccess;
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
  bool bsuccess = true;
  double vgs1(0.0), vgdd1(0.0), vbs1(0.0),vgb1(0.0), vdds1(0.0);
  double * staVector = extData.nextStaVectorRawPtr;
  double * currStaVector = extData.currStaVectorRawPtr;

  bool tmpBool = updateIntermediateVars ();
  bsuccess = bsuccess && tmpBool;

  // voltage drops:
  staVector[li_state_vbdd] = vbdd;
  staVector[li_state_vbs]  = vbs;
  staVector[li_state_vgps] = vgps;
  staVector[li_state_vdds] = vdds;
  staVector[li_state_D1vd] = D1vd;

  // now the meyer capacitances
  // we didn't calculate these charges in update IntermediateVars
  // but we did calculate the voltage drops and capacitances.
  // first store the capacitances themselves:
  staVector[li_state_capgs]  = capgs;
  staVector[li_state_capgdd] = capgdd;
  staVector[li_state_capgb]  = capgb;

  // now the charges
  // BE CAREFUL!  We can only do Q=CV for DCOP!  Otherwise it's
  // supposed to be *INTEGRATED*:
  // Q = int(t0,t1)C(V)*dV --- and we approximate that by
  // Q(t1)-Q(t0) = CBar*(V(t1)-V(t0)) where CBar is the average.
  // Now with Meyer back averaging, Capxx is the average between the last
  // time step and this one.  So we gotta do the right thing for non-DCOP
  // when backaverage is on.

  if((getSolverState().dcopFlag))
  {
    qgs  = Capgs *vgps;
    qgdd = Capgdd*vgpdd;
    qgb  = Capgb *Vgpb;
  }
  else
  {
    // get the ones from last time step
    qgs   = currStaVector[li_state_qgs];
    qgdd  = currStaVector[li_state_qgdd];
    qgb   = currStaVector[li_state_qgb];
    // get the voltage drops, too
    vgs1  = currStaVector[li_state_vgps];
    vbs1  = currStaVector[li_state_vbs];
    vdds1 = currStaVector[li_state_vdds];

    vgb1  = vgs1-vbs1;
    vgdd1 = vgs1-vdds1;

    // NOW we can calculate the charge update
    qgs  += Capgs*(vgps-vgs1);
    qgdd += Capgdd*(vgpdd-vgdd1);
    qgb  += Capgb*((vgps-vbs)-vgb1);
  }

  staVector[li_state_qgs]  = qgs;
  staVector[li_state_qgdd] = qgdd;
  staVector[li_state_qgb]  = qgb;

  // and the diode parasitic capacitors
  // these charges were set in updateIntermediateVars
  staVector[li_state_qbd] = qbd;
  staVector[li_state_qbs] = qbs;

// diode #1 quatntities

  staVector[li_state_D1DIOcapCharge]   = D1DIOcapCharge;

  // In the case of a charge, need to handle this:
  // Ensure dQ/dt = 0 for initial step after dcOP by setting the history equal
  // to the current value
  if( !getSolverState().dcopFlag && getSolverState().initTranFlag_ && !getSolverState().newtonIter )
  {
    currStaVector[li_state_D1DIOcapCharge]   = D1DIOcapCharge;
  }

  staVector[ li_state_von ] = von;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/09/01
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  double * staDerivVector = (extData.nextStaDerivVectorRawPtr);

  cqgs  = staDerivVector[li_state_qgs];
  cqgdd = staDerivVector[li_state_qgdd];
  cqgb  = staDerivVector[li_state_qgb];
  cqbd  = staDerivVector[li_state_qbd];
  cqbs  = staDerivVector[li_state_qbs];

  D1DIOcapCurrent  = staDerivVector[li_state_D1DIOcapCharge];

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::applyScale
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
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
    l = model_.l0;
  if (!given("W"))
    w = model_.w0;

  if(model_.drainResistance != 0)
  {
    drainCond = numberParallel/model_.drainResistance;
  }
  else if (model_.given("RSH"))
  {
    if(model_.sheetResistance != 0)
    {
      drainCond =
        numberParallel/(model_.sheetResistance*drainSquares);
    }
    else
    {
      drainCond = 0;
    }
  }
  else
  {
    drainCond = 0;
  }
  if(model_.sourceResistance != 0)
  {
    sourceCond = numberParallel/model_.sourceResistance;
  }
  else if (model_.given("RSH"))
  {
    if(model_.sheetResistance != 0)
    {
      sourceCond =
        numberParallel/(model_.sheetResistance*sourceSquares);
    }
    else
    {
      sourceCond = 0;
    }
  }
  else
  {
    sourceCond = 0;
  }

  if(model_.given("RG"))
  {
    if(model_.gateResistance != 0)
    {
      gateCond = numberParallel/model_.gateResistance;
    }
    else
    {
      gateCond = 0;
    }
  }
  else
  {
    gateCond = 0;
  }

  w *= numberParallel;

  gddd = model_.driftParamA;
  if(gddd != 0) gddd = 1.0/gddd;
  draindriftCond = gddd;

  if(l - 2 * model_.latDiff <=0)
  {
    UserError(*this) << "Effective channel length less than zero.";
  }

  // values for diode #1
  D1DIOarea = 1;

  EffectiveLength = l - 2*model_.latDiff;
  GateSourceOverlapCap = model_.gateSourceOverlapCapFactor * w;
  GateDrainOverlapCap = model_.gateDrainOverlapCapFactor * w;
  GateBulkOverlapCap = model_.gateBulkOverlapCapFactor * EffectiveLength;
  OxideCap        = model_.oxideCapFactor * EffectiveLength * w;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << " L = " << l << std::endl;
    Xyce::dout() << " W = " << w<< std::endl;
    Xyce::dout() << " drainArea = " << drainArea<< std::endl;
    Xyce::dout() << " sourceArea = " << sourceArea<< std::endl;
    Xyce::dout() << " drainSquares = " << drainSquares<< std::endl;
    Xyce::dout() << " sourceSquares = " << sourceSquares<< std::endl;
    Xyce::dout() << " drainPerimeter = " << drainPerimeter<< std::endl;
    Xyce::dout() << " sourcePerimeter = " << sourcePerimeter<< std::endl;
    Xyce::dout() << " drainCond = " << drainCond<< std::endl;
    Xyce::dout() << " sourceCond = " << sourceCond << std::endl;
    Xyce::dout() << " draindriftCond = " << draindriftCond<< std::endl;
    Xyce::dout() << " temp = " << temp<< std::endl;
  }

  // now set the temperature related stuff
  updateTemperature(temp);

  return true;
}

// Additional Declarations


//-----------------------------------------------------------------------------
// Function      : Instance::UCCMqmeyer
// Purpose       : Compute the MOS overlap capacitances as functions of the
//                 device terminal voltages
//
// Special Notes :
//
// Scope         : public
// Creator       : pmc
// Creation Date : 02/17/2004
//-----------------------------------------------------------------------------
bool Instance::UCCMqmeyer (
   double xvgs,    // initial voltage gate-source (mode > 0)
   double xvgdd,   // initial voltage gate-drain (mode > 0)
   double vgb,     // initial voltage gate-bulk
   double von_local,     // threshold voltage
   double vddsat_local,  // saturation drain voltage
   double & capgs_local, // non-constant portion of g-s overlap capacitance
   double & capgdd_local,// non-constant portion of g-d overlap capacitance
   double & capgb_local, // non-constant portion of g-b overlap capacitance
   double phi,
   double cox      // oxide capactiance
 )
{
  double vdds_local(0.0);
  double vddif(0.0);
  double vddif1(0.0);
  double vddif2(0.0);
  double vgst(0.0);
  double etavt(0.0);
  double cgc(0.0);
  double x(0.0);

#define W        w
#define L        EffectiveLength
#define ETA      model_.eta
#define TOX      model_.oxideThickness

  vgst = xvgs-von_local;
  if (vgst <= -phi)
  {
    capgb_local  = cox/2;
    capgs_local  = 0;
    capgdd_local = 0;
  }
  else if (vgst <= -phi/2)
  {
    capgb_local  = -vgst*cox/(2*phi);
    capgs_local  = 0;
    capgdd_local = 0;
  }
  else if (vgst <= 0)
  {
    capgb_local  = -vgst*cox/(2*phi);
    capgs_local  = vgst*cox/(1.5*phi)+cox/3;
    capgdd_local = 0;
  }
  else
  {
    vdds_local = xvgs-xvgdd;
    capgb_local = 0;
    if (vddsat_local <= vdds_local)
    {
      capgs_local  = cox/3;
      capgdd_local = 0;
    }
    else
    {
      vddif   = 2.0*vddsat_local-vdds_local;
      vddif1  = vddsat_local-vdds_local/*-1.0e-12*/;
      vddif2  = vddif*vddif;
      capgdd_local = cox*(1.0-vddsat_local*vddsat_local/vddif2)/3;
      capgs_local  = cox*(1.0-vddif1*vddif1/vddif2)/3;
    }
  }

  // at this point we have unmodified Xyce values
  // the following is the UCCM adjustment

  if(model_.cv == 2 && vddsat_local != 0)
  {
    vdds_local    = fabs(xvgs-xvgdd);
    vdds_local    = vdds_local/pow(1+pow(vdds_local/vddsat_local,model_.mc),
                                            1.0/model_.mc);
    vddif   = 2.0*vddsat_local-vdds_local;
    vddif1  = vddsat_local-vdds_local;
    vddif2  = vddif*vddif;
    etavt   = ETA*vt;
    x       = vgst/etavt;
    cgc     = L*W/(TOX/CONSTEPSSI+etavt/(CONSTQ*n0)*exp(-x));
    capgdd_local = cgc*(1.0-vddsat_local*vddsat_local/vddif2)/3;
    capgs_local  = cgc*(1.0-vddif1*vddif1/vddif2)/3;
  }

  return  true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::UCCMMeyercap
// Purpose       :
//
//
// Special Notes :
//
// Scope         : public
// Creator       : pmc
// Creation Date : 02/17/2004
//-----------------------------------------------------------------------------
bool Instance::UCCMMeyercap (
   double xvgs,    // initial voltage gate-source  (mode > 0)
   double xvgdd,   // initial voltage gate-drain  (mode > 0)
   double vgb,     // initial voltage gate-bulk
   double & cgs,
   double & cgd,
   double & cgb
  )
{
   double vbs_local(0.0);
   double vdds_local(0.0);
   double zetanb(0.0);
   double xi(0.0);
   double xisqrt(0.0);
   double nsd(0.0);
   double nss(0.0);
   double tnss(0.0);
   double tnsd(0.0);
   double DeltaVT(0.0);
   double DnsdVgs(0.0);
   double DnsdVgd(0.0);
   double DnsdVgb(0.0);
   double DnssVgs(0.0);
   double DnssVgd(0.0);
   double DnssVgb(0.0);
   double DndepVgs(0.0);
   double DndepVgd(0.0);
   double DndepVgb(0.0);
   double UI(0.0);
   double VI(0.0);
   double DxiVgs(0.0);
   double DxiVgd(0.0);
   double DxiVgb(0.0);
   double DUIVgs(0.0);
   double DUIVgd(0.0);
   double DUIVgb(0.0);
   double DVIVgs(0.0);
   double DVIVgd(0.0);
   double DVIVgb(0.0);     // With use of the chain rule, this routine may
   double DqiVgs(0.0);     // be used to calculate the Ward capacitors.
   double DqiVgd(0.0);     //
   double DqiVgb(0.0);     //  Cs, gs=dqs/dVgs= DqiVgs*fp+qiVgs*dfp
   double DqbVgs(0.0);     //  fp: partitioning factor
   double DqbVgd(0.0);
   double DqbVgb(0.0);
   double DqgVgs(0.0);
   double DqgVgd(0.0);
   double DqgVgb(0.0);
   double SigmaVgs(0.0);   // This Sigma-parameter acounts for the DIBL
   double SigmaVgd(0.0);   // effect on the capasitors. The physical
   double etavt(0.0);
   double VFB(0.0);
   double mqWL(0.0);
   double gamma(0.0);
   double A(0.0);
   double ALPHA(0.0);

   vbs_local   = xvgs-vgb;
   vdds_local  = xvgs-xvgdd;
   A     = CONSTQ/model_.oxideCapFactor;
   VFB   =model_.vfb;
   ALPHA =model_.alpha;

   if(xvgs > VFB+vbs_local)
   {
      gamma= model_.gammas0;
      mqWL=(-CONSTQ)*EffectiveLength*(w);
      etavt=model_.eta*vt;
      nss=2.0*n0*log(1+0.5*exp((xvgs-von)/etavt));
      nsd=2.0*n0*log(1+0.5*exp((xvgs-von-ALPHA*vdds_local)/etavt));

      if (nss<1e-36)  {nss=1.0e-36;}
      if (nsd<1e-36)  {nsd=1.0e-36;}
      zetanb = (1-ALPHA)/ALPHA;
      xi = gamma*gamma/(A*A)+4/(A*A)*(vgb-VFB)-4.0/A*nss;
      xisqrt=sqrt(xi);
      tnsd = etavt/nsd+A;
      tnss = etavt/nss+A;
      if (vbs_local <= 0) { DeltaVT= 0.5*gamma/sqrt(tPhi-vbs_local); }
      if ((vbs_local > 0)&&(vbs_local <= 2*tPhi)) { DeltaVT=0.5*gamma/sqrt(tPhi); }
      if ((vbs_local > 0)&&(vbs_local > 2*tPhi)) { DeltaVT=0; }

      DnsdVgs=(1+DeltaVT-ALPHA+SigmaVgs)/tnsd;
      DnsdVgd=(-SigmaVgd+ALPHA)/tnsd;
      DnsdVgb=-DeltaVT/tnsd;

      DnssVgs= (1+DeltaVT+SigmaVgs)/tnss;
      DnssVgd=-SigmaVgd/tnss;
      DnssVgb=-DeltaVT/tnss;

      UI= etavt/2.0*(nsd*nsd-nss*nss) + A/3.0*(nsd*nsd*nsd-nss*nss*nss);
      VI= etavt*(nsd-nss) + A/2.0*(nsd*nsd-nss*nss);
      DUIVgs = etavt*(nsd*DnsdVgs-nss*DnssVgs) +
                             A*(nsd*nsd*DnsdVgs-nss*nss*DnssVgs);
      DUIVgd= etavt*(nsd*DnsdVgd-nss*DnssVgd) +
                             A*(nsd*nsd*DnsdVgd-nss*nss*DnssVgd);
      DUIVgb= etavt*(nsd*DnsdVgb-nss*DnssVgb) +
                             A*(nsd*nsd*DnsdVgb-nss*nss*DnssVgb);

      DVIVgs= etavt*(DnsdVgs-DnssVgs) + A*(nsd*DnsdVgs-nss*DnssVgs);
      DVIVgd= etavt*(DnsdVgd-DnssVgd) + A*(nsd*DnsdVgd-nss*DnssVgd);
      DVIVgb= etavt*(DnsdVgb-DnssVgb) + A*(nsd*DnsdVgb-nss*DnssVgb);
      if (VI != 0)
      {
         DqiVgs= mqWL*(DUIVgs*VI-DVIVgs*UI)/(VI*VI);
         DqiVgd= mqWL*(DUIVgd*VI-DVIVgd*UI)/(VI*VI);
         DqiVgb= mqWL*(DUIVgb*VI-DVIVgb*UI)/(VI*VI);
      }
      else
      {
         DqiVgs=mqWL*DnssVgs;
         DqiVgd=mqWL*DnssVgd;
         DqiVgb=mqWL*DnssVgb;
      }

      DxiVgs= -4.0/A*DnssVgs;
      DxiVgd= -4.0/A*DnssVgd;
      DxiVgb= 4.0/(A*A)-4/A*DnssVgb;

      if (xisqrt !=0)
      {DndepVgs= gamma/4.0*1.0/xisqrt*DxiVgs;
       DndepVgd= gamma/4.0*1.0/xisqrt*DxiVgd;
       DndepVgb= gamma/4.0*1.0/xisqrt*DxiVgb;
      }
      else
      {
       DndepVgs=0;
       DndepVgd=0;
       DndepVgb=0;
      }

      DqbVgs = mqWL*(DndepVgs-zetanb*DnssVgs)+zetanb*DqiVgs;
      DqbVgd = mqWL*(DndepVgd-zetanb*DnssVgd)+zetanb*DqiVgd;
      DqbVgb = mqWL*(DndepVgb-zetanb*DnssVgb)+zetanb*DqiVgb;

      DqgVgs= -DqiVgs-DqbVgs;
      DqgVgd= -DqiVgd-DqbVgd;
      DqgVgb= -DqiVgb-DqbVgb;
   }
   else
   {
      // Below flat band
      DqgVgs=0.0;
      DqgVgd=0.0;
      DqgVgb=model_.oxideCapFactor * EffectiveLength * w;
   }

   if (DqgVgs<0)  {DqgVgs=0;}   // Because of truncation error cgb may turn
   if (DqgVgd<0)  {DqgVgd=0;}   // negative in extreme inversion, while cgd
   if (DqgVgb<0)  {DqgVgb=0;}   // and cgs may turn negative in accumulation
   cgs= 0.5*DqgVgs;
   cgd= 0.5*DqgVgd;
   cgb= 0.5*DqgVgb;

   return true;
}

//----------------------------------------------------------------------------
// Function      : Instance::UCCMCharges
// Purpose       :
//
//
// Special Notes :
//
// Scope         : public
// Creator       : pmc
// Creation Date : 02/17/2004
//----------------------------------------------------------------------------
bool Instance::UCCMCharges( double vgs, double vgdd,
                           double vgb, double & qD, double & qS, double & qB)
{
  double VT(0.0);     // Threshold voltage
  double nss(0.0);    // Inversion charge density at source
  double nsd(0.0);    // Inversion charge density at drain
  double etavth(0.0);
  double mqWL(0.0);
  double nsdsqr(0.0);
  double nsssqr(0.0);
  double arg1(0.0);
  double arg2(0.0);
  double qn(0.0);     // Inversion charge
  double VFB(0.0);    // Flat-band voltage
  double xisqrt(0.0);
  double ndeps(0.0);  // Depletion charge density
  double zetanb(0.0);
  double fp1(0.0);    // Partitioning factor
  double gamma(0.0);
  double ALPHA(0.0);
  double A(0.0);
  double cox(0.0);    // oxide capasity per unit area
  double vdds_local(0.0);

  cox   = model_.oxideCapFactor;
  gamma = model_.gammas0;
  A     = CONSTQ/cox;
  ETA   = model_.eta;
  ALPHA = model_.alpha;
  VFB   = model_.vfb;

  vdds_local=vgs-vgdd;

  if (vgb > VFB)
  {
     VT = von;
     etavth = ETA*temp*CONSTKoverQ;

     nsd=2*n0*log(1+0.5*exp((vgs-VT-ALPHA*vdds_local)/etavth));
     nss=2*n0*log(1+0.5*exp((vgs-VT)/etavth));

     nsdsqr=nsd*nsd;
     nsssqr=nss*nss;
     mqWL=-CONSTQ*EffectiveLength*(w);

     arg1= 0.5*etavth*(nsdsqr-nsssqr)+1/3.0*A*(nsdsqr*nsd-nsssqr*nss);
     arg2= etavth*(nsd-nss) +0.5*A*(nsdsqr-nsssqr);

     if (arg2==0)        //||(fabs(nss-nsd)<1e-12))
        { qn=mqWL*nss;}
     else {qn = mqWL*arg1/arg2;}
     xisqrt= sqrt(gamma*gamma/(A*A)+4.0/(A*A)*(vgb-VFB)-4.0/A*nss);
     ndeps = -gamma*gamma/(2.0*A)+ gamma/2.0*xisqrt;
     zetanb =(1-ALPHA)/ALPHA;
     mqWL=-CONSTQ*EffectiveLength*(w);
     qB =mqWL*(ndeps-zetanb*nss)+zetanb*qn;

    switch(model_.fpe)
     {
      case 1:
           fp1=model_.xqc;
         break;

      case 2:
        {
          double vsat(0.0);
          double beta(0.0);
          double a0(0.0);
          double vdsabs(0.0);
          double m(0.0);
          double fp1denum(0.0);

          vdsabs=fabs(vdds_local);
          vsat=vddsat;
          if (vgs-VT > 1e-36)
          {
              if((vsat > 1e-36)&&(vdds_local>1e-36))
                    { beta=vdsabs/vsat;}else{beta=1e-36;}
              m  = model_.mcv;
              a0 = (model_.xqc - 0.5)/vsat;
              fp1denum = exp(1/m*log(1+exp(m*log(beta))));
              if (fp1denum > 1e-36)
                  {fp1= 0.5+ a0*vdsabs/fp1denum;}
              else {fp1=0.5;}

          } else {fp1=0.5;}
        }
         break;
      case 3:
        {
          double nd0(0.0);
          double ns0(0.0);
          double arg1_loc(0.0);
          double arg2_loc(0.0);
          double arg3_loc(0.0);

          nd0=nss/(A*etavth);
          ns0=nss/(A*etavth);
          arg1_loc= nd0*nd0*nd0*(1/3.0+nd0*(3.0/8.0+1/10.0*nd0))
                - (1.0/2.0+1.0/3.0*nd0)*(1.0+1/2.0*ns0)*nd0*nd0*ns0
                + ns0*ns0*ns0*(1.0/6.0+ns0*(5.0/24.0+1.0/15.0*ns0));
          arg2_loc = 1.0/2.0*(nd0*nd0-ns0*ns0)+1.0/3.0*(nd0*nd0*nd0-ns0*ns0*ns0);
          arg3_loc = nd0-ns0+1.0/2.0*(nd0*nd0-ns0*ns0);
          if ((arg2_loc==0)||(arg3_loc==0)) {fp1=0.5;}
          else { fp1 = 1-arg1_loc/(arg2_loc*arg3_loc);}
        }
          break;
      default:
        {
          UserWarning(*this) << "Partitioning model does not exist";
          return true;
        }
     }
  }
  else
  {    //  Below flatband
        qB = -EffectiveLength*(w)*cox*(vgb-VFB);
        qn  = 0;
  }

  qS = -fp1*qn;
  qD = -(1-fp1)*qn;
  qB = -qB;

   return true;
}

//----------------------------------------------------------------------------
// Function      : Instance::UCCMcvon
// Purpose       : Compute the threshold voltage
//
// Special Notes : ERK.  The "local" named variables (von_local, dvonvbs_local)
//                 are named that way to avoid conflicts with the instance
//                 variables of the same name.
//
// Scope         : public
// Creator       : pmc
// Creation Date : 02/17/2004
//----------------------------------------------------------------------------
bool Instance::UCCMcvon(double vbs_local,
                                       double & von_local, double & dvonvbs_local)
{
   double PhiMinVbs = tPhi - vbs_local;
   double sarg(0.0);
   double dsrgdb(0.0);

// vtoo, vthLimit calculated in updateTemp

   if(vtoo-vbs_local > vthLimit) {
     von_local = vtoo + gammal*vthLimit;
     dvonvbs_local = 0.0;
     return true;
   }
   if(PhiMinVbs > 0.0)
   {
     sarg = sqrt(PhiMinVbs);
     dsrgdb = -0.5/sarg;
   } else
   {
     sarg = 0;
     dsrgdb = 0;
   }
   von_local = vtoo + gammas*sarg - gammal*PhiMinVbs;
   dvonvbs_local = gammas*dsrgdb + gammal;

   return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::UCCMmosa1
// Purpose       : Compute currents and conductances in drain region
//
//
// Special Notes :
//
// Scope         : public
// Creator       : pmc
// Creation Date : 02/17/2004
//-----------------------------------------------------------------------------

bool Instance::UCCMmosa1(double xvgs, double xvdds,
                      double dvonvbs, double & cdraindrift_loc, double & vsate)
{
  double etavt(0.0);
  double vgt(0.0);
  double vgt0(0.0);
  double sigma(0.0);
  double vgte(0.0);
  double isat(0.0);
  double mu(0.0);
  double ns(0.0);
  double a(0.0);
  double b(0.0);
  double d(0.0);
  double g(0.0);
  double h(0.0);
  double q(0.0);
  double r(0.0);
  double t(0.0);
  double u(0.0);
  double x(0.0),y(0.0),z(0.0);

  double gch(0.0);
  double gchi(0.0);
  double rt(0.0);
  double vl(0.0);
  double vl2(0.0);
  double ichoo(0.0);

  double dichoodvds(0.0);
  double dichoodvgs(0.0);
  double dichoodvbs(0.0);
  double dichooisat(0.0);
  double dichoodgch(0.0);
  double delgchgchi(0.0);
  double disatdvds(0.0);
  double disatdvgs(0.0);
  double disatdvbs(0.0);
  double dgchidvds(0.0);
  double dgchidvgs(0.0);
  double dgchidvbs(0.0);
  double disatvgte(0.0);
  double disatgchi(0.0);
  double dvgtedvgt(0.0);
  double dvgtdvds(0.0);
  double dvgtdvgs(0.0);
  double dvgtdvbs(0.0);
  double dnsdvgt(0.0);
  double dnsdvds(0.0);
  double dnsdvgs(0.0);
  double dnsdvbs(0.0);
  double dmudvgte(0.0);
  double dmudvds(0.0);
  double dmudvgs(0.0);
  double dmudvbs(0.0);

  static int output=0;

#define ETA      model_.eta
#define RS       model_.sourceResistance
#define RD       model_.drainResistance
#define VSIGMAT  model_.vsigmat
#define VSIGMA   model_.vsigma
#define SIGMA0   model_.sigma0
#define THETA    model_.theta
#define LAMBDA   model_.lambda
#define DELTA    model_.delta
#define VMAX     model_.maxDriftVel
#define TOX      model_.oxideThickness

#define EXP_MAX     150.0

  if(output == 1) {
    Xyce::dout() << "                 " << std::endl;
    Xyce::dout() << "ETA     = " << ETA << std::endl;
    Xyce::dout() << "RS      = " << RS << std::endl;
    Xyce::dout() << "RD      = " << RD << std::endl;
    Xyce::dout() << "VSIGMAT = " << VSIGMAT << std::endl;
    Xyce::dout() << "VSIGMA  = " << VSIGMA << std::endl;
    Xyce::dout() << "SIGMA0  = " << SIGMA0 << std::endl;
    Xyce::dout() << "THETA   = " << THETA << std::endl;
    Xyce::dout() << "LAMBDA  = " << LAMBDA << std::endl;
    Xyce::dout() << "DELTA   = " << DELTA << std::endl;
    Xyce::dout() << "VMAX    = " << VMAX << std::endl;
    Xyce::dout() << "TOX     = " << TOX << std::endl;
    output=0;
  }

  etavt    = ETA*vt;
  rt       = RS+RD;
  vgt0     = xvgs - von;
  a        = exp((vgt0-VSIGMAT)/VSIGMA);
  sigma    = SIGMA0/(1+a);
  vgt      = vgt0+sigma*xvdds;
  b        = 0.5*vgt/vt-1;
  q        = sqrt(model_.deltaSqr+b*b);
  vgte     = vt*(1+b+1+q);
  u        = 1+THETA*(vgte+2*von)/TOX;
  mu       = tSurfMob/u;

  x        = vgt/etavt;
  if (x > 50.0)
  {
     ns = n0*2.0*x;
  }
  else if (x < -30)
  {
     ns = n0*exp(x);
  }
  else
  {
     ns  = 2.0*n0*log(1+0.5*exp(x));
  }

  if(ns < 1.0e-38) {
    cdraindrift_loc = 0.0;
    vsate = 0.0;
    gm   = 0.0;
    gdds = 0.0;
    gmbs = 0.0;
    mm1  = 0.0;
    dmm1vgs = 0.0;
    dmm1vds = 0.0;
    dmm1vbs = 0.0;
    return true;
  }

  gchi       = gchi0*mu*ns;
  gch        = gchi/(1+gchi*rt);
  t          = VMAX*L;
//  vl         = t/mu;
  vl         = t/tSurfMob;
  vl2        = vl*vl;
  d          = sqrt(1+2*gchi*RS + vgte*vgte/vl2);
  r          = gchi*vgte;
  h          = 1+gchi*RS+d;
  isat       = r/h;
  vsate      = isat/gch;
  y          = xvdds/(vsate);

  if(fabs(y) > EXP_MAX) z = 0;
  else    z  = 1/(cosh(y)*cosh(y));

  // drain current
  ichoo      = isat*(1+LAMBDA*xvdds)*tanh(y);

  dvgtedvgt  = 0.5*(1+b/q);
  dnsdvgt    = n0/(etavt*(exp(-x)+0.5));
  dmudvgte   = -tSurfMob*THETA/(TOX*u*u);

  dvgtdvds   = sigma;
  dvgtdvgs   = 1 - SIGMA0*a*xvdds/(VSIGMA*(1+a)*(1+a));
  dvgtdvbs   = -dvonvbs*dvgtdvgs;

  dnsdvds    = dnsdvgt*dvgtdvds;
  dnsdvgs    = dnsdvgt*dvgtdvgs;
  dnsdvbs    = dnsdvgt*dvgtdvbs;
  dmudvds    = (dmudvgte*dvgtedvgt)*dvgtdvds;;
  dmudvgs    = (dmudvgte*dvgtedvgt)*dvgtdvgs;;
  dmudvbs    = (dmudvgte*dvgtedvgt)*dvgtdvbs + 2*dmudvgte*dvonvbs;

  dgchidvds  = gchi0*(mu*dnsdvds + ns*dmudvds);
  dgchidvgs  = gchi0*(mu*dnsdvgs + ns*dmudvgs);
  dgchidvbs  = gchi0*(mu*dnsdvbs + ns*dmudvbs);

  disatvgte  = gchi/h - gchi*vgte*vgte/vl2/(d*h*h);
  disatgchi  = vgte/h - gchi*vgte*RS*(1+1/d)/(h*h);

  disatdvds  = (disatvgte*dvgtedvgt)*dvgtdvds + disatgchi*dgchidvds;
  disatdvgs  = (disatvgte*dvgtedvgt)*dvgtdvgs + disatgchi*dgchidvgs;
  disatdvbs  = (disatvgte*dvgtedvgt)*dvgtdvbs + disatgchi*dgchidvbs;

  dichoodgch = (1+LAMBDA*xvdds)*xvdds*z;
  dichooisat = (1+LAMBDA*xvdds)*(tanh(y) - gch*xvdds*z/isat);
  g          = 1+gchi*rt;
  delgchgchi = 1/(g*g);

  dichoodvds = dichooisat*disatdvds + (dichoodgch*delgchgchi)*dgchidvds
                       + isat*LAMBDA*tanh(y) + (1+LAMBDA*xvdds)*gch*z;
  dichoodvgs = dichooisat*disatdvgs + (dichoodgch*delgchgchi)*dgchidvgs;
  dichoodvbs = dichooisat*disatdvbs + (dichoodgch*delgchgchi)*dgchidvbs;


  cdraindrift_loc = ichoo;
  gm      = dichoodvgs;
  gdds    = dichoodvds;
  gmbs    = dichoodvbs;

// the mm1 factors are for the ISUB calc which we don't use
  mm1     = 0.0;
  dmm1vgs = 0.0;
  dmm1vds = 0.0;
  dmm1vbs = 0.0;

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Instance::UCCMmosa2
// Purpose       : Compute currents and conductances in
//               :                            drain and substrate regions
//
//
// Special Notes :
//
// Scope         : public
// Creator       : pmc
// Creation Date : 02/17/2004
//-----------------------------------------------------------------------------

bool Instance::UCCMmosa2(double xvgs, double xvdds,
                      double dvonvbs, double & cdraindrift_loc, double & vsate)
{
  double etavt(0.0);
  double vgt(0.0);
  double vgt0(0.0);
  double sigma(0.0);
  double vgte(0.0);
  double isat(0.0);
  double mu(0.0);
  double ns(0.0);
  double a(0.0);
  double b(0.0);
  double c(0.0);
  double d(0.0);
  double e(0.0);
  double f(0.0);
  double g(0.0);
  double h(0.0);
  double p(0.0);
  double q(0.0);
  double r(0.0);
  double s(0.0);
  double t(0.0);
  double u(0.0);
  double x,y,z(0.0);
  double gch(0.0);
  double gchi(0.0);
  double rt(0.0);
  double vl(0.0);
  double vl2(0.0);
  double ichoo(0.0);
  double gmoo(0.0);
  double gdsoo(0.0);
  double gbsoo(0.0);
  double icho(0.0);
  double gmo(0.0);
  double gdso(0.0);
  double gbso(0.0);
  double vsat(0.0);
  double vdso(0.0);
  double vdse(0.0);
  double cox(0.0);
  double temp1(0.0);
  double qs(0.0);
  double dqsvgt(0.0);
  double dqsvbs(0.0);
  double delidgch(0.0);
  double delgchgchi(0.0);
  double delgchins(0.0);
  double dvgtevgt(0.0);
  double delidvsate(0.0);
  double delvsateisat(0.0);
  double delisatvgte(0.0);
  double delisatgchi(0.0);
  double delisatvl(0.0);
  double dgchivgt(0.0);
  double delvsategch(0.0);
  double delidvds(0.0);
  double dsigmavgs(0.0);
  double delnsvgt(0.0);
  double dvsatevgt(0.0);
  double dvgtvgs(0.0);
  double dmuvon(0.0);
  double dmuvgt(0.0);
  double dvgtvon(0.0);
  double delvsatqs(0.0);
  double delvsatmu(0.0);
  double dvsatvgt(0.0);
  double dvsatvbs(0.0);
  double delvdsevdso(0.0);
  double delvdsevsat(0.0);
  double dvdsevgs(0.0);
  double dvdsevds(0.0);
  double dvdsevbs(0.0);

  double del(0.0);
  double delmm1vdso(0.0);
  double ddelvgs(0.0);
  double ddelvds(0.0);
  double ddelvbs(0.0);

  // These "loc" variables are named such to avoid name conflicts with
  // instance variables.
  double mm1_loc(0.0);
  double dmm1vgs_loc(0.0);
  double dmm1vds_loc(0.0);
  double dmm1vbs_loc(0.0);

  static int output(0);

#define ETA      model_.eta
#define RS       model_.sourceResistance
#define RD       model_.drainResistance
#define VSIGMAT  model_.vsigmat
#define VSIGMA   model_.vsigma
#define SIGMA0   model_.sigma0
#define THETA    model_.theta
#define LAMBDA   model_.lambda
#define DELTA    model_.delta
#define ALPHA    model_.alpha
#define VMAX     model_.maxDriftVel
#define TOX      model_.oxideThickness
#define W        w
#define L        EffectiveLength
#define LS       model_.ls
#define M        model_.m
#define INVM     1.0/model_.m
#define MC       model_.mc
#define INVMC    1.0/model_.mc
#define RSUB     model_.rsub
#define AI       model_.ai
#define BI       model_.bi
#define MD       model_.md
#define INVMD    1.0/model_.md
#define DELMAX   model_.delmax

   if(output == 1) {
     Xyce::dout() << "                 " << std::endl;
     Xyce::dout() << "ETA     = " << ETA << std::endl;
     Xyce::dout() << "RS      = " << RS << std::endl;
     Xyce::dout() << "RD      = " << RD << std::endl;
     Xyce::dout() << "VSIGMAT = " << VSIGMAT << std::endl;
     Xyce::dout() << "VSIGMA  = " << VSIGMA << std::endl;
     Xyce::dout() << "SIGMA0  = " << SIGMA0 << std::endl;
     Xyce::dout() << "THETA   = " << THETA << std::endl;
     Xyce::dout() << "LAMBDA  = " << LAMBDA << std::endl;
     Xyce::dout() << "DELTA   = " << DELTA << std::endl;
     Xyce::dout() << "ALPHA   = " << ALPHA << std::endl;
     Xyce::dout() << "VMAX    = " << VMAX << std::endl;
     Xyce::dout() << "TOX     = " << TOX << std::endl;
     Xyce::dout() << "W       = " << W << std::endl;
     Xyce::dout() << "L       = " << L << std::endl;
     Xyce::dout() << "LS      = " << LS << std::endl;
     Xyce::dout() << "M       = " << M << std::endl;
     Xyce::dout() << "INVM    = " << INVM << std::endl;
     Xyce::dout() << "MC      = " << MC << std::endl;
     Xyce::dout() << "INVMC   = " << INVMC << std::endl;
     Xyce::dout() << "RSUB    = " << RSUB << std::endl;
     Xyce::dout() << "AI      = " << AI << std::endl;
     Xyce::dout() << "BI      = " << BI << std::endl;
     Xyce::dout() << "MD      = " << MD << std::endl;
     Xyce::dout() << "INVMD   = " << INVMD << std::endl;
     Xyce::dout() << "DELMAX  = " << DELMAX << std::endl;
     Xyce::dout() << "VP      = " << vp << std::endl;
     output=0;
   }


  etavt    = ETA*vt;
  rt       = RS+RD;
  vgt0     = xvgs - von;
  a        = exp((vgt0-VSIGMAT)/VSIGMA);
  sigma    = SIGMA0/(1+a);
  vgt      = vgt0+sigma*xvdds;
  b        = 0.5*vgt/vt-1;
  q        = sqrt(model_.deltaSqr+b*b);
  vgte     = vt*(1+b+1+q);
  u        = 1+THETA*(vgte+2*von)/TOX;
  mu       = tSurfMob/u;

  x        = vgt/etavt;
  if (x > 50.0)
  {
     ns = n0*2.0*x;
  }
  else if (x < -30)
  {
     ns = n0*exp(x);
  }
  else
  {
     ns  = 2.0*n0*log(1+0.5*exp(x));
  }
  //ns  = 2.0*n0*log(1+0.5*c);

  if(ns < 1.0e-38)
  {
    cdraindrift_loc = 0.0;
    vsate = 0.0;
    gm   = 0.0;
    gdds = 0.0;
    gmbs = 0.0;

    mm1 = mm1_loc  = 0.0;
    dmm1vgs = dmm1vgs_loc = 0.0;
    dmm1vds = dmm1vds_loc = 0.0;
    dmm1vbs = dmm1vbs_loc = 0.0;
    return true;
  }

  gchi     = gchi0*mu*ns;
  gch      = gchi/(1+gchi*rt);
  t        = VMAX*L;
  vl       = t/mu;
  vl2      = vl*vl;
  d        = sqrt(1+2*gchi*RS + vgte*vgte/vl2);
  r        = gchi*vgte;
  isat     = r/(1+gchi*RS+d);
  vsate   = isat/gch;
  y        = xvdds/(vsate);
  z        = 0.5*(cosh(2*y)+1);

//  Old model
//  e        = pow(xvdds/(*vsate),M);
//  f        = pow(1+e,INVM);
//  delidgch = xvdds*(1+LAMBDA*xvdds)/f;

  ichoo        = isat*(1+LAMBDA*xvdds)*tanh(y);
  delidgch     = ichoo/gch;
  g            = 1+gchi*rt;
  delgchgchi   = 1/(g*g);
  delgchins    = gchi0*mu;
  x            = vgt/etavt;
  delnsvgt     = n0/(etavt*(exp(-x) + 0.5));
  dvgtevgt     = 0.5*(1+b/q);

//  delidvds     = gch*(1+2*LAMBDA*xvdds)/f-
//                 ichoo*pow(vdds/(*vsate),M-1)/((*vsate)*(1+e));
//  delidvsate   = ichoo*e/((*vsate)*(1+e));
  delidvds     = ichoo*LAMBDA/(1+LAMBDA*xvdds) + gch*(1+LAMBDA*xvdds)/z;
  delidvsate   = gch*(1+LAMBDA*xvdds)*xvdds/(vsate)/z;

  delvsateisat = 1/gch;
  delvsategch  = -(vsate)/gch;
  h            = 1+gchi*RS+d;
  delisatgchi  = (vgte*h - gchi*vgte*RS*(1+1/d))/(h*h);
  delisatvgte  = (gchi*h - gchi*vgte*vgte/(vl2*d))/(h*h);
  delisatvl    = isat*isat*vgte/(gchi*d*vl2*vl);
  s            = THETA/(tSurfMob*TOX);
  dgchivgt     = delgchins*delnsvgt-gchi*mu*s*dvgtevgt;
  dsigmavgs    = -SIGMA0*a/(VSIGMA*((1+a)*(1+a)));
  p            = delgchgchi*dgchivgt;
  dvsatevgt    = delvsateisat*(delisatgchi*dgchivgt+
                 (delisatvgte+delisatvl*t*s)*dvgtevgt)+ delvsategch*p;
  g            = delidgch*p+delidvsate*dvsatevgt;
  dvgtvgs      = (1+xvdds*dsigmavgs);
  gmoo         = g*dvgtvgs;
  gdsoo        = delidvds+g*sigma;

  dmuvon       = -(2-dvgtevgt*dvgtvgs)*mu*mu*s;
  dmuvgt       = -mu*mu*THETA*dvgtevgt/(tSurfMob*TOX);
  dvgtvon      = -dvgtvgs;

  if(THETA == 0.0)
  {
    gbsoo = -(gmoo)*dvonvbs;
  }
  else
  {
    double dgchivon = -delgchins*delnsvgt*dvgtvgs+gchi*dmuvon/mu;
    p        = delgchgchi*dgchivon;
    double dvsatevon = delvsateisat*(delisatgchi*dgchivon+
      delisatvgte*dvgtevgt*dvgtvon+delisatvl*(-vl2/t)*dmuvon)+delvsategch*p;
    gbsoo = (delidgch*p+delidvsate*dvsatevon)*dvonvbs;
  }

//  mosa2: New Version

  cox    = CONSTEPSOX/TOX;
  temp1  = cox*RS/CONSTQ;
  qs     = CONSTQ*(ns-ichoo*temp1);
  dqsvgt = CONSTQ*(delnsvgt-g*temp1);
  dqsvbs = CONSTQ*(delnsvgt*dvgtvon*dvonvbs-gbsoo*temp1);

  // A more precise vsat model for deltal and isub calculations
  if(model_.dtype == CONSTNMOS)
  {
    a    = 2*qs*VMAX*L;
    b    = qs*mu+2*cox*ALPHA*VMAX*L;
    temp1 = b*b;
    vsat = a/b;
    delvsatqs = 2*VMAX*L*(b-qs*mu)/temp1;
    delvsatmu = -a*qs/temp1;
    dvsatvgt  = delvsatqs*dqsvgt+delvsatmu*dmuvgt;
    dvsatvbs  = delvsatqs*dqsvbs+delvsatmu*dmuvon*dvonvbs;
  }
  else
  {
    a = VMAX*L/mu;
    b = 2/(cox*ALPHA*VMAX*L);
    c = sqrt(1+qs*mu*b);
    vsat = a*(c-1);
    delvsatqs = a*(0.5*b*mu/c);
    delvsatmu = -vsat/mu+a*(0.5*b*qs/c);
    dvsatvgt  = delvsatqs*dqsvgt+delvsatmu*dmuvgt;
    dvsatvbs  = delvsatqs*dqsvbs+delvsatmu*dmuvon*dvonvbs;
  }

  // Effective intrinsic drain-source voltage calculation
  vdso = xvdds-ichoo*rt;
  a    = pow(vdso/vsat,MC);
  b    = pow(1+a,INVMC);
  vdse = vdso/b;
  delvdsevdso  = (1-a/(1+a))/b;
  delvdsevsat = vdse*a/(vsat*(1+a));

  // DeltaL calculation
  double deltal(0.0);
  double deldeltalvdse(0.0);
  double deldeltalvdso(0.0);
  double deldeltalmu(0.0);
  double ddeltalvgs(0.0);
  double ddeltalvds(0.0);
  double ddeltalvbs(0.0);

  a      = 1+(vdso-vdse)/vp;
  b      = W*mu*cox*vdse*RS/L;
  c      = 1+vdse/vp+b;
  d      = LS*log10(a);
  deltal = d/c;
  e      = 1/(1-deltal/L);
  icho   = ichoo*e;
  f      = log(10.0);
  deldeltalvdse = (-LS/(f*a*vp)*c-d*(1/vp+W*mu*cox*RS/L))/(c*c);
  deldeltalvdso = LS/(c*a*f*vp);
  deldeltalmu   = -d*W*cox*vdse*RS/(L*c*c);
  temp1         = -gmoo*rt;
  dvdsevgs      = delvdsevsat*dvsatvgt*dvgtvgs+delvdsevdso*temp1;
  ddeltalvgs    = deldeltalvdse*dvdsevgs + deldeltalmu*dmuvgt*dvgtvgs
                       + deldeltalvdso*temp1;
  gmo           = e*(gmoo+ichoo*ddeltalvgs*e/L);
  temp1         = 1-gdsoo*rt;
  dvdsevds      = delvdsevsat*dvsatvgt*sigma+delvdsevdso*temp1;
  ddeltalvds    = deldeltalvdse*dvdsevds + deldeltalmu*dmuvgt*sigma
                       + deldeltalvdso*temp1;
  gdso          = e*(gdsoo+ichoo*ddeltalvds*e/L);
  temp          = -gbsoo*rt;
  dvdsevbs      = delvdsevsat*dvsatvbs+delvdsevdso*temp1;
  ddeltalvbs    = deldeltalvdse*dvdsevbs + deldeltalmu*dmuvon*dvonvbs
                       + deldeltalvdso*temp1;;
  gbso          = e*(gbsoo+ichoo*ddeltalvbs*e/L);

  a   = AI/BI;
  b   = vdso-vdse+vt;
  c   = exp(-LS*BI/b);
  mm1_loc = a*b*c;
  if(RSUB != 0.0)
  {
    d   = 0.5*W*(ALPHA-1)*cox*RSUB/L;
    del = d*mm1_loc*mu*vdse;
    f   = pow(del/DELMAX,MD);
    g   = pow(1+f,INVMD);
    h   = (1-f/(1+f))/g;
    del = del/g;
    e   = 1/(1-del);//1+del;
  }
  else
  {
    d = 0;
    h = 1;
    e = 1;
  }

  cdraindrift_loc = icho*e;
  delmm1vdso = a*c+mm1_loc*LS*BI/(b*b);
  dmm1vgs_loc    = delmm1vdso*(-gmoo*rt-dvdsevgs);
  ddelvgs    = h*d*(mu*vdse*dmm1vgs_loc + mm1_loc*vdse*dmuvgt*dvgtvgs
                  + mu*mm1_loc*dvdsevgs);
  gm       = e*(gmo+icho*ddelvgs*e);      //gmo*e+icho*ddelvgs;
  dmm1vds_loc  = delmm1vdso*(1-gdsoo*rt-dvdsevds);
  ddelvds  = h*d*(mu*vdse*dmm1vds_loc + mm1_loc*vdse*dmuvgt*sigma
                  + mu*mm1_loc*dvdsevds);
  gdds     = e*(gdso+icho*ddelvds*e);     //gdso*e+icho*ddelvds;
  dmm1vbs_loc  = delmm1vdso*(-gbsoo*rt-dvdsevbs);
  ddelvbs  = h*d*(mu*vdse*dmm1vbs_loc + mm1_loc*vdse*dmuvon*dvonvbs
                  + mu*mm1_loc*dvdsevbs);
  gmbs     = e*(gbso+icho*ddelvbs*e);     //gbso*e+icho*ddelvbs;


  // saving some quantities for substrate current calculations
  mm1     = mm1_loc;
  dmm1vgs = dmm1vgs_loc;
  dmm1vds = dmm1vds_loc;
  dmm1vbs = dmm1vbs_loc;

   return true;
}


//-----------------------------------------------------------------------------
// VDMOS Master functions:
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

      double * currStaVector = mi.extData.currStaVectorRawPtr;

      bool btmp = mi.updateIntermediateVars ();
      bsuccess = bsuccess && btmp;

      double vgs1(0.0), vgdd1(0.0), vbs1(0.0),vgb1(0.0), vdds1(0.0);

      // voltage drops:
      staVec[mi.li_state_vbdd] = mi.vbdd;
      staVec[mi.li_state_vbs]  = mi.vbs;
      staVec[mi.li_state_vgps] = mi.vgps;
      staVec[mi.li_state_vdds] = mi.vdds;
      staVec[mi.li_state_D1vd] = mi.D1vd;

      // now the meyer capacitances
      // we didn't calculate these charges in update IntermediateVars
      // but we did calculate the voltage drops and capacitances.
      // first store the capacitances themselves:
      staVec[mi.li_state_capgs]  = mi.capgs;
      staVec[mi.li_state_capgdd] = mi.capgdd;
      staVec[mi.li_state_capgb]  = mi.capgb;

      // now the charges
      // BE CAREFUL!  We can only do Q=CV for DCOP!  Otherwise it's
      // supposed to be *INTEGRATED*:
      // Q = int(t0,t1)C(V)*dV --- and we approximate that by
      // Q(t1)-Q(t0) = CBar*(V(t1)-V(t0)) where CBar is the average.
      // Now with Meyer back averaging, Capxx is the average between the last
      // time step and this one.  So we gotta do the right thing for non-DCOP
      // when backaverage is on.

      if((getSolverState().dcopFlag))
      {
        mi.qgs  = mi.Capgs *mi.vgps;
        mi.qgdd = mi.Capgdd*mi.vgpdd;
        mi.qgb  = mi.Capgb *mi.Vgpb;
      }
      else
      {
        // get the ones from last time step
        mi.qgs   = currStaVector[mi.li_state_qgs];
        mi.qgdd  = currStaVector[mi.li_state_qgdd];
        mi.qgb   = currStaVector[mi.li_state_qgb];
        // get the voltage drops, too
        vgs1  = currStaVector[mi.li_state_vgps];
        vbs1  = currStaVector[mi.li_state_vbs];
        vdds1 = currStaVector[mi.li_state_vdds];

        vgb1  = vgs1-vbs1;
        vgdd1 = vgs1-vdds1;

        // NOW we can calculate the charge update
        mi.qgs  += mi.Capgs*(mi.vgps-vgs1);
        mi.qgdd += mi.Capgdd*(mi.vgpdd-vgdd1);
        mi.qgb  += mi.Capgb*((mi.vgps-mi.vbs)-vgb1);
      }

      staVec[mi.li_state_qgs]  = mi.qgs;
      staVec[mi.li_state_qgdd] = mi.qgdd;
      staVec[mi.li_state_qgb]  = mi.qgb;

      // and the diode parasitic capacitors
      // these charges were set in updateIntermediateVars
      staVec[mi.li_state_qbd] = mi.qbd;
      staVec[mi.li_state_qbs] = mi.qbs;

    // diode #1 quatntities

      staVec[mi.li_state_D1DIOcapCharge]   = mi.D1DIOcapCharge;

      // In the case of a charge, need to handle this:
      // Ensure dQ/dt = 0 for initial step after dcOP by setting the history equal
      // to the current value
      if( !getSolverState().dcopFlag && getSolverState().initTranFlag_ && !getSolverState().newtonIter )
      {
        currStaVector[mi.li_state_D1DIOcapCharge]   = mi.D1DIOcapCharge;
      }

      staVec[ mi.li_state_von ] = mi.von;

    }
//  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateSecondaryState ( double * staDerivVec, double * stoVec )
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
    {
      Instance & mi = *(*it);

      mi.cqgs  = staDerivVec[mi.li_state_qgs];
      mi.cqgdd = staDerivVec[mi.li_state_qgdd];
      mi.cqgb  = staDerivVec[mi.li_state_qgb];
      mi.cqbd  = staDerivVec[mi.li_state_qbd];
      mi.cqbs  = staDerivVec[mi.li_state_qbs];

      mi.D1DIOcapCurrent  = staDerivVec[mi.li_state_D1DIOcapCharge];

    }

  return true;
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
    Instance & mi = *(*it);

    // F-vector:
      double * dFdxdVp = mi.extData.dFdxdVpVectorRawPtr;

    double coef_Jdxp(0.0), gd_Jdxp(0.0);
    double gmin1 = getDeviceOptions().gmin;
    double Dtype =mi.getModel().dtype;

    double ceqbs = Dtype*(mi.cbs+mi.cqbs);
    double ceqbd = Dtype*(mi.cbd+mi.cqbd);

    double D1current = Dtype*mi.D1cdeq; // don't add in the capacitor stuff


    fVec[mi.li_Drain] += (mi.Idraindrift + mi.Irdsshunt - D1current);

    if (mi.Igate != 0.0)
    {
      fVec[mi.li_Gate]  += mi.Igate;
      fVec[mi.li_GatePrime]   += -mi.Igate;
    }


    fVec[mi.li_Source]      += (mi.Isource - mi.Irdsshunt + mi.Ird1rs);
    fVec[mi.li_Bulk]        += (ceqbs + ceqbd);
    fVec[mi.li_DrainPrime]  += (-mi.Idrain-(ceqbd - mi.cdreq));
    fVec[mi.li_SourcePrime] += (-mi.Isource-(ceqbs + mi.cdreq));
    fVec[mi.li_DrainDrift]  += (mi.Idrain - mi.Idraindrift);
    fVec[mi.li_D1Prime]     += (D1current - mi.Ird1rs);

    // limiter terms:
    if (!mi.origFlag)
    {
      // bulk
      coef_Jdxp = Dtype*( + ((mi.gbd-gmin1))*(mi.vbdd-mi.vbdd_orig)
                          + ((mi.gbs-gmin1))*(mi.vbs-mi.vbs_orig));
      dFdxdVp[mi.li_Bulk] += coef_Jdxp;

      // drain-prime
      coef_Jdxp = Dtype*(-((mi.gbd-gmin1))*
                        (mi.vbdd-mi.vbdd_orig)+mi.gdds*(mi.vdds-mi.vdds_orig)
                        +mi.Gm*((mi.mode>0)?(mi.vgps-mi.vgps_orig):(mi.vgpdd-mi.vgpdd_orig))
                        +mi.Gmbs*((mi.mode>0)?(mi.vbs-mi.vbs_orig):(mi.vbdd-mi.vbdd_orig)));
      dFdxdVp[mi.li_DrainPrime] += coef_Jdxp;

      // source prime
      coef_Jdxp = Dtype*(-((mi.gbs-gmin1))*(mi.vbs-mi.vbs_orig)
                        -mi.gdds*(mi.vdds-mi.vdds_orig)
                        -mi.Gm*((mi.mode>0)?(mi.vgps-mi.vgps_orig):(mi.vgpdd-mi.vgpdd_orig))
                        -mi.Gmbs*((mi.mode>0)?(mi.vbs-mi.vbs_orig):(mi.vbdd-mi.vbdd_orig)));


      dFdxdVp[mi.li_SourcePrime] += coef_Jdxp;
      gd_Jdxp = -mi.D1gd * (mi.D1vd-mi.D1vd_orig);
      dFdxdVp[mi.li_Drain] += gd_Jdxp;
      dFdxdVp[mi.li_D1Prime] -= gd_Jdxp;
    }

    // Q-vector:
    double * dQdxdVp = mi.extData.dQdxdVpVectorRawPtr;

    double qcoef_Jdxp(0.0);

    double Qeqbs = Dtype*mi.qbs;
    double Qeqbd = Dtype*mi.qbd;
    double Qeqgb = Dtype*mi.qgb;
    double Qeqgs = Dtype*mi.qgs;
    double Qeqgdd = Dtype*mi.qgdd;


    qVec[mi.li_Bulk]        += (Qeqbs + Qeqbd - Qeqgb);
    qVec[mi.li_DrainPrime]  += -(Qeqbd + Qeqgdd);
    qVec[mi.li_GatePrime]   += (Qeqgs+Qeqgdd+Qeqgb);
    qVec[mi.li_SourcePrime] += (-(Qeqbs + Qeqgs));
    qVec[mi.li_D1Prime] +=  mi.D1DIOcapCharge;
    qVec[mi.li_Drain]   += -mi.D1DIOcapCharge;

    // voltlim section:
    if (!mi.origFlag)
    {
      // bulk
      qcoef_Jdxp = Dtype*(-mi.Capgb*(mi.vgps-mi.vgps_orig-mi.vbs+mi.vbs_orig)
                        + (+mi.Capgb)*(mi.vbdd-mi.vbdd_orig)
                        + (+mi.capbs)*(mi.vbs-mi.vbs_orig));
      dQdxdVp[mi.li_Bulk] += qcoef_Jdxp;

      // drain-prime
      qcoef_Jdxp = Dtype*(-mi.Capgdd*(mi.vgpdd-mi.vgpdd_orig)-mi.capbd*(mi.vbdd-mi.vbdd_orig));
      dQdxdVp[mi.li_DrainPrime] += qcoef_Jdxp;

      // gate-prime
      qcoef_Jdxp = Dtype*(mi.Capgdd*(mi.vgpdd-mi.vgpdd_orig)+mi.Capgs*(mi.vgps-mi.vgps_orig)+
                        mi.Capgb*(mi.vgps-mi.vgps_orig-mi.vbs+mi.vbs_orig));
      dQdxdVp[mi.li_GatePrime] += qcoef_Jdxp;

      // source-prime
      qcoef_Jdxp = Dtype*(-mi.Capgs*(mi.vgps-mi.vgps_orig)-mi.capbs*(mi.vbs-mi.vbs_orig));
      dQdxdVp[mi.li_SourcePrime] += qcoef_Jdxp;
      qcoef_Jdxp = -mi.D1capd*(mi.D1vd-mi.D1vd_orig);
      dQdxdVp[mi.li_D1Prime] -= qcoef_Jdxp;
      dQdxdVp[mi.li_Drain]   += qcoef_Jdxp;
    }

    if( mi.loadLeadCurrent )
    {
      leadF[mi.li_branch_data_d] = (mi.Idraindrift + mi.Irdsshunt - D1current);
      leadF[mi.li_branch_data_s] = (mi.Isource - mi.Irdsshunt + mi.Ird1rs);
      leadF[mi.li_branch_data_g] = 0;
      leadF[mi.li_branch_data_b] = (ceqbs + ceqbd);
      // case where optional nodes become external nodes:
      if (mi.Igate != 0.0)
      {
        leadF[mi.li_branch_data_g]  += mi.Igate;
      }
      if( !mi.gateCond )
      {
        // G' is G
        leadF[mi.li_branch_data_g] += -mi.Igate;
      }

      if( !mi.sourceCond )
      {
        // S' is S
        leadF[mi.li_branch_data_s] += (-mi.Isource-(ceqbs + mi.cdreq));
      }

      if( !mi.drainCond )
      {
        // Ddrift is D'
      }

      if( !mi.model_.D1DIOconductance )
      {
        // D1' is S
        leadF[mi.li_branch_data_s] += (D1current - mi.Ird1rs);
      }

      leadQ[mi.li_branch_data_d] = -mi.D1DIOcapCharge;
      leadQ[mi.li_branch_data_s] = 0;
      leadQ[mi.li_branch_data_g] = 0;
      leadQ[mi.li_branch_data_b] = (Qeqbs + Qeqbd - Qeqgb);
      // case where optional nodes become external nodes:
      if( !mi.gateCond )
      {
        // G' is G
        leadQ[mi.li_branch_data_g] += (Qeqgs+Qeqgdd+Qeqgb);
      }

      if( !mi.sourceCond )
      {
        // S' is S
        leadQ[mi.li_branch_data_s] += (-(Qeqbs + Qeqgs));
      }

      if( !mi.drainCond )
      {
        // Ddrift is D'
      }

      if( !mi.model_.D1DIOconductance )
      {
        // D1' is S
        leadQ[mi.li_branch_data_s] += mi.D1DIOcapCharge;
      }

      junctionV[mi.li_branch_data_d] = solVec[mi.li_Drain] - solVec[mi.li_Source];
      junctionV[mi.li_branch_data_g] = solVec[mi.li_Gate] - solVec[mi.li_Source];
      junctionV[mi.li_branch_data_s] = 0.0;
      junctionV[mi.li_branch_data_b] = 0.0 ; 
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

    *mi.f_DrainEquDrainNodePtr += mi.dIdd_dVd + mi.gdsshunt + mi.D1gd;

    *mi.f_DrainEquSourceNodePtr -= mi.gdsshunt;

    *mi.f_DrainEquDrainDriftNodePtr -= mi.dIdd_dVd;

    *mi.f_DrainEquD1PrimeNodePtr -= mi.D1gd;

    if (mi.gateCond != 0)
    {

      *mi.f_GateEquGateNodePtr += mi.gateCond;

      *mi.f_GateEquGatePrimeNodePtr -= mi.gateCond;
    }


    *mi.f_SourceEquDrainNodePtr -= mi.gdsshunt;

    *mi.f_SourceEquSourceNodePtr += mi.sourceCond + mi.gdsshunt + mi.D1gspr;

    *mi.f_SourceEquD1PrimeNodePtr += -mi.D1gspr;

    if (mi.sourceCond != 0.0)
    {

      *mi.f_SourceEquSourcePrimeNodePtr -= mi.sourceCond;
    }


    *mi.f_BulkEquBulkNodePtr += mi.gbs+mi.gbd;

    *mi.f_BulkEquDrainPrimeNodePtr -= mi.gbd;


    *mi.f_BulkEquSourcePrimeNodePtr -= mi.gbs;


    *mi.f_DrainPrimeEquBulkNodePtr += -mi.gbd+mi.Gmbs;

    *mi.f_DrainPrimeEquDrainPrimeNodePtr += mi.drainCond+mi.gdds+mi.gbd+mi.revsum;

    *mi.f_DrainPrimeEquGatePrimeNodePtr += mi.Gm;

    *mi.f_DrainPrimeEquSourcePrimeNodePtr += -mi.gdds-mi.nrmsum;
    if (mi.drainCond != 0.0)
    {

      *mi.f_DrainPrimeEquDrainDriftNodePtr -= mi.drainCond;
    }

    if (mi.gateCond != 0)
    {

      *mi.f_GatePrimeEquGateNodePtr -= mi.gateCond;

      *mi.f_GatePrimeEquGatePrimeNodePtr += mi.gateCond;
    }

    if (mi.sourceCond != 0.0)
    {

      *mi.f_SourcePrimeEquSourceNodePtr -= mi.sourceCond;
    }

    *mi.f_SourcePrimeEquBulkNodePtr -= mi.gbs+mi.Gmbs;

    *mi.f_SourcePrimeEquDrainPrimeNodePtr -= mi.gdds+mi.revsum;

    *mi.f_SourcePrimeEquGatePrimeNodePtr -= mi.Gm;

    *mi.f_SourcePrimeEquSourcePrimeNodePtr += mi.sourceCond+mi.gdds+mi.gbs+mi.nrmsum;


    *mi.f_DrainDriftEquDrainNodePtr -= mi.dIdd_dVd;

    if (mi.drainCond != 0.0)
    {

      *mi.f_DrainDriftEquDrainPrimeNodePtr -= mi.drainCond;
    }

    *mi.f_DrainDriftEquDrainDriftNodePtr += mi.dIdd_dVd + mi.drainCond;


    *mi.f_D1PrimeEquDrainNodePtr -= mi.D1gd;

    *mi.f_D1PrimeEquSourceNodePtr += -mi.D1gspr;

    *mi.f_D1PrimeEquD1PrimeNodePtr += mi.D1gd + mi.D1gspr;

    // Q-matrix:
    if (!getSolverState().dcopFlag)
    {

      *mi.q_BulkEquBulkNodePtr += +mi.capbs+mi.capbd+mi.Capgb;

      *mi.q_BulkEquDrainPrimeNodePtr -= +mi.capbd;

      *mi.q_BulkEquGatePrimeNodePtr -= mi.Capgb;

      *mi.q_BulkEquSourcePrimeNodePtr -= +mi.capbs;


      *mi.q_DrainPrimeEquBulkNodePtr += -mi.capbd;

      *mi.q_DrainPrimeEquDrainPrimeNodePtr += +mi.capbd+mi.Capgdd;

      *mi.q_DrainPrimeEquGatePrimeNodePtr += -mi.Capgdd;


      *mi.q_GatePrimeEquBulkNodePtr -= mi.Capgb;

      *mi.q_GatePrimeEquDrainPrimeNodePtr -= mi.Capgdd;

      *mi.q_GatePrimeEquGatePrimeNodePtr += +mi.Capgdd+mi.Capgs+mi.Capgb;

      *mi.q_GatePrimeEquSourcePrimeNodePtr -= mi.Capgs;


      *mi.q_SourcePrimeEquBulkNodePtr -= +mi.capbs;

      *mi.q_SourcePrimeEquGatePrimeNodePtr -= mi.Capgs;

      *mi.q_SourcePrimeEquSourcePrimeNodePtr += +mi.capbs+mi.Capgs;


      *mi.q_DrainEquDrainNodePtr += mi.D1capd;

      *mi.q_DrainEquD1PrimeNodePtr -= mi.D1capd;


      *mi.q_D1PrimeEquDrainNodePtr -= mi.D1capd;

      *mi.q_D1PrimeEquD1PrimeNodePtr += mi.D1capd;
    }
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

    dFdx[mi.li_Drain][mi.ADrainEquDrainNodeOffset] += mi.dIdd_dVd + mi.gdsshunt + mi.D1gd;

    dFdx[mi.li_Drain][mi.ADrainEquSourceNodeOffset] -= mi.gdsshunt;

    dFdx[mi.li_Drain][mi.ADrainEquDrainDriftNodeOffset] -= mi.dIdd_dVd;

    dFdx[mi.li_Drain][mi.ADrainEquD1PrimeNodeOffset] -= mi.D1gd;

    if (mi.gateCond != 0)
    {

      dFdx[mi.li_Gate][mi.AGateEquGateNodeOffset] += mi.gateCond;

      dFdx[mi.li_Gate][mi.AGateEquGatePrimeNodeOffset] -= mi.gateCond;
    }


    dFdx[mi.li_Source][mi.ASourceEquDrainNodeOffset] -= mi.gdsshunt;

    dFdx[mi.li_Source][mi.ASourceEquSourceNodeOffset] += mi.sourceCond + mi.gdsshunt + mi.D1gspr;

    dFdx[mi.li_Source][mi.ASourceEquD1PrimeNodeOffset] += -mi.D1gspr;

    if (mi.sourceCond != 0.0)
    {

      dFdx[mi.li_Source][mi.ASourceEquSourcePrimeNodeOffset] -= mi.sourceCond;
    }


    dFdx[mi.li_Bulk][mi.ABulkEquBulkNodeOffset] += mi.gbs+mi.gbd;

    dFdx[mi.li_Bulk][mi.ABulkEquDrainPrimeNodeOffset] -= mi.gbd;

    dFdx[mi.li_Bulk][mi.ABulkEquSourcePrimeNodeOffset] -= mi.gbs;


    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquBulkNodeOffset] += -mi.gbd+mi.Gmbs;

    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquDrainPrimeNodeOffset] += mi.drainCond+mi.gdds+mi.gbd+mi.revsum;

    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquGatePrimeNodeOffset] += mi.Gm;

    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquSourcePrimeNodeOffset] += -mi.gdds-mi.nrmsum;
    if (mi.drainCond != 0.0)
    {

      dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquDrainDriftNodeOffset] -= mi.drainCond;
    }

    if (mi.gateCond != 0)
    {

      dFdx[mi.li_GatePrime][mi.AGatePrimeEquGateNodeOffset] -= mi.gateCond;

      dFdx[mi.li_GatePrime][mi.AGatePrimeEquGatePrimeNodeOffset] += mi.gateCond;
    }

    if (mi.sourceCond != 0.0)
    {

      dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquSourceNodeOffset] -= mi.sourceCond;
    }

    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquBulkNodeOffset] -= mi.gbs+mi.Gmbs;

    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquDrainPrimeNodeOffset] -= mi.gdds+mi.revsum;

    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquGatePrimeNodeOffset] -= mi.Gm;

    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquSourcePrimeNodeOffset] += mi.sourceCond+mi.gdds+mi.gbs+mi.nrmsum;


    dFdx[mi.li_DrainDrift][mi.ADrainDriftEquDrainNodeOffset] -= mi.dIdd_dVd;

    if (mi.drainCond != 0.0)
    {

      dFdx[mi.li_DrainDrift][mi.ADrainDriftEquDrainPrimeNodeOffset] -= mi.drainCond;
    }

    dFdx[mi.li_DrainDrift][mi.ADrainDriftEquDrainDriftNodeOffset] += mi.dIdd_dVd + mi.drainCond;


    dFdx[mi.li_D1Prime][mi.AD1PrimeEquDrainNodeOffset] -= mi.D1gd;

    dFdx[mi.li_D1Prime][mi.AD1PrimeEquSourceNodeOffset] += -mi.D1gspr;

    dFdx[mi.li_D1Prime][mi.AD1PrimeEquD1PrimeNodeOffset] += mi.D1gd + mi.D1gspr;

    // Q-matrix:
    if (!getSolverState().dcopFlag)
    {

      dQdx[mi.li_Bulk][mi.ABulkEquBulkNodeOffset] += +mi.capbs+mi.capbd+mi.Capgb;

      dQdx[mi.li_Bulk][mi.ABulkEquDrainPrimeNodeOffset] -= +mi.capbd;

      dQdx[mi.li_Bulk][mi.ABulkEquGatePrimeNodeOffset] -= mi.Capgb;

      dQdx[mi.li_Bulk][mi.ABulkEquSourcePrimeNodeOffset] -= +mi.capbs;


      dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquBulkNodeOffset] += -mi.capbd;

      dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquDrainPrimeNodeOffset] += +mi.capbd+mi.Capgdd;

      dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquGatePrimeNodeOffset] += -mi.Capgdd;


      dQdx[mi.li_GatePrime][mi.AGatePrimeEquBulkNodeOffset] -= mi.Capgb;

      dQdx[mi.li_GatePrime][mi.AGatePrimeEquDrainPrimeNodeOffset] -= mi.Capgdd;

      dQdx[mi.li_GatePrime][mi.AGatePrimeEquGatePrimeNodeOffset] += +mi.Capgdd+mi.Capgs+mi.Capgb;

      dQdx[mi.li_GatePrime][mi.AGatePrimeEquSourcePrimeNodeOffset] -= mi.Capgs;


      dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquBulkNodeOffset] -= +mi.capbs;

      dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquGatePrimeNodeOffset] -= mi.Capgs;

      dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquSourcePrimeNodeOffset] += +mi.capbs+mi.Capgs;


      dQdx[mi.li_Drain][mi.ADrainEquDrainNodeOffset] += mi.D1capd;

      dQdx[mi.li_Drain][mi.ADrainEquD1PrimeNodeOffset] -= mi.D1capd;


      dQdx[mi.li_D1Prime][mi.AD1PrimeEquDrainNodeOffset] -= mi.D1capd;

      dQdx[mi.li_D1Prime][mi.AD1PrimeEquD1PrimeNodeOffset] += mi.D1capd;
    }

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
      ((deviceMap.find("M")!=deviceMap.end()) && (levelSet.find(18)!=levelSet.end())))
  {
    MOSFET1::registerDevice();

    Config<Traits>::addConfiguration()
      .registerDevice("m", 18)
      .registerModelType("pmos", 18)
      .registerModelType("nmos", 18);
  }
}

} // namespace VDMOS
} // namespace Device
} // namespace Xyce

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
// Purpose        : This file implements the BSIM3 MOSFET model.  It
//                  is intended to be compatible with the Berkeley SPICE
//                  (3f5) version, BSIM3 version 3.2.2.
//
// Special Notes  :
//
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 01/12/01
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#ifndef Xyce_USE_BSIM3_CONST
#include <N_DEV_Const.h>
#endif

#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MOSFET_B3.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_DEV_MOSFET1.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Math.h>
#include <N_ANP_NoiseData.h>

// ---------- BSIM3 constants ---------------
// Normally, these should be obtained from N_DEV_Const.h
// I have them here, for purposes of comparison to spice.
#ifdef Xyce_USE_BSIM3_CONST

#define CONSTMAX_EXP 5.834617425e+14
#define CONSTMIN_EXP 1.713908431e-15

#define CONSTEXP_THRESHOLD 34.0
#define CONSTEPSOX 3.453133e-11
#define CONSTEPSSI 1.03594e-10

#define CONSTQ 1.60219e-19

#define CONSTDELTA_1 0.02
#define CONSTDELTA_2 0.02
#define CONSTDELTA_3 0.02
#define CONSTDELTA_4 0.02

#define CONSTroot2    sqrt(2.0)
#define CONSTCtoK     (273.15)
#define CONSTREFTEMP  (300.15)

#define CONSTboltz    1.3806226e-23
#define CONSTKoverQ   8.617087e-5  // Kb / q  where q = 1.60219e-19

#define CONSTvt0     (CONSTboltz * (27.0 +CONSTCtoK)/CONSTQ)

#define CONSTEg300   (1.1150877) // band gap for Si at T=300.15K (room temp)
#define CONSTEg0     (1.16)      // band gap for Si at T=0K. (eV)
#define CONSTalphaEg (7.02e-4)   // (eV/K)
#define CONSTbetaEg  (1108.0)    // (K)

#define CONSTNi0     (1.45e10)   // carrier concentration at room temp.

#define CONSTNMOS 1
#define CONSTPMOS -1

#endif

namespace Xyce {
namespace Device {

namespace MOSFET_B3 {

void Traits::loadInstanceParameters(ParametricData<MOSFET_B3::Instance> &p)
{
  p.addPar("TEMP",0.0,&MOSFET_B3::Instance::temp)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(STANDARD)
    .setCategory(CAT_NONE)
    .setDescription("Device temperature")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3InstanceSens);

  p.addPar("L",0.0,&MOSFET_B3::Instance::l)
    .setOriginalValueStored(true)
    .setUnit(U_METER)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Channel length")
    .setLengthScaling(true)
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3InstanceSens);

  p.addPar("W",0.0,&MOSFET_B3::Instance::w)
    .setOriginalValueStored(true)
    .setUnit(U_METER)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Channel width")
    .setLengthScaling(true)
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3InstanceSens);

  p.addPar("AD",0.0,&MOSFET_B3::Instance::drainArea)
    .setUnit(U_METER2)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Drain diffusion area")
    .setAreaScaling(true)
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3InstanceSens);

  p.addPar("AS",0.0,&MOSFET_B3::Instance::sourceArea)
    .setUnit(U_METER2)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Source diffusion area")
    .setAreaScaling(true)
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3InstanceSens);

  p.addPar("NRD",1.0,&MOSFET_B3::Instance::drainSquares)
    .setUnit(U_SQUARES)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Multiplier for RSH to yield parasitic resistance of drain")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3InstanceSens);

  p.addPar("NRS",1.0,&MOSFET_B3::Instance::sourceSquares)
    .setUnit(U_SQUARES)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Multiplier for RSH to yield parasitic resistance of source")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3InstanceSens);

  p.addPar("PD",0.0,&MOSFET_B3::Instance::drainPerimeter)
    .setUnit(U_METER)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Drain diffusion perimeter")
    .setLengthScaling(true)
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3InstanceSens);

  p.addPar("PS",0.0,&MOSFET_B3::Instance::sourcePerimeter)
    .setUnit(U_METER)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Source diffusion perimeter")
    .setLengthScaling(true)
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3InstanceSens);

  p.addPar("M",1.0,&MOSFET_B3::Instance::numberParallel)
    .setUnit(U_NONE)
    .setCategory(CAT_CONTROL)
    .setDescription("Multiplier for M devices connected in parallel")
    .setAnalyticSensitivityAvailable(false);

  p.addPar("IC1",0.0,&MOSFET_B3::Instance::icVDS)
    .setGivenMember(&MOSFET_B3::Instance::icVDSGiven)
    .setUnit(U_VOLT)
    .setCategory(CAT_VOLT)
    .setDescription("Initial condition on Vds")
    .setAnalyticSensitivityAvailable(false);

  p.addPar("IC2",0.0,&MOSFET_B3::Instance::icVGS)
    .setGivenMember(&MOSFET_B3::Instance::icVGSGiven)
    .setUnit(U_VOLT)
    .setCategory(CAT_VOLT)
    .setDescription("Initial condition on Vgs")
    .setAnalyticSensitivityAvailable(false);

  p.addPar("IC3",0.0,&MOSFET_B3::Instance::icVBS)
    .setGivenMember(&MOSFET_B3::Instance::icVBSGiven)
    .setUnit(U_VOLT)
    .setCategory(CAT_VOLT)
    .setDescription("Initial condition on Vbs")
    .setAnalyticSensitivityAvailable(false);

  // Set up non-double precision variables:
  p.addPar("NQSMOD",0,&MOSFET_B3::Instance::nqsMod)
    .setUnit(U_NONE)
    .setCategory(CAT_CONTROL)
    .setDescription("Flag for NQS model")
    .setAnalyticSensitivityAvailable(false);

  p.addPar("OFF",false,&MOSFET_B3::Instance::OFF)
    .setUnit(U_LOGIC)
    .setCategory(CAT_VOLT)
    .setDescription("Initial condition of no voltage drops accross device")
    .setAnalyticSensitivityAvailable(false);

  p.addPar("DTEMP", 0.0, &MOSFET_B3::Instance::dtemp)
    .setGivenMember(&MOSFET_B3::Instance::dtempGiven)
    .setUnit(U_DEGC)
    .setCategory(CAT_TEMP)
    .setDescription("Device delta temperature")
    .setAnalyticSensitivityAvailable(false);

  // This tells the parser that IC1,IC2,and IC3 are to be input as a vector of "IC"
  p.makeVector ("IC",3);
}

void Traits::loadModelParameters(ParametricData<MOSFET_B3::Model> &p)
{
  // Set up map for normal (double) param variables:
  p.addPar("TOX",150.e-10,&MOSFET_B3::Model::tox)
    .setOriginalValueStored(true)
    .setUnit(U_METER)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Gate oxide thickness")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("TOXM",0.0,&MOSFET_B3::Model::toxm)
    .setUnit(U_METER)
    .setCategory(CAT_PROCESS)
    .setDescription("Gate oxide thickness used in extraction")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("CDSC",2.4e-4,&MOSFET_B3::Model::cdsc)
    .setUnit(U_FARADMM2)
    .setCategory(CAT_DC)
    .setDescription("Drain/source to channel coupling capacitance")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("CDSCB",0.0,&MOSFET_B3::Model::cdscb)
    .setUnit(U_FVM1MM2)
    .setCategory(CAT_DC)
    .setDescription("Body-bias sensitivity of CDSC")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("CDSCD",0.0,&MOSFET_B3::Model::cdscd)
    .setUnit(U_FVM1MM2)
    .setCategory(CAT_DC)
    .setDescription("Drain-bias sensitivity of CDSC")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("CIT",0.0,&MOSFET_B3::Model::cit)
    .setUnit(U_FARADMM2)
    .setCategory(CAT_DC)
    .setDescription("Interface trap capacitance")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("NFACTOR",1.0,&MOSFET_B3::Model::nfactor)
    .setUnit(U_NONE)
    .setCategory(CAT_DC)
    .setDescription("Subthreshold swing factor")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("XJ",0.15e-6,&MOSFET_B3::Model::xj)
    .setUnit(U_METER)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Junction depth")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("VSAT",8.0e4,&MOSFET_B3::Model::vsat)
    .setUnit(U_MSM1)
    .setCategory(CAT_DC)
    .setDescription("Saturation velocity at temp = TNOM")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("AT",3.3e4,&MOSFET_B3::Model::at)
    .setUnit(U_MSM1)
    .setCategory(CAT_TEMP)
    .setDescription("Temperature coefficient for saturation velocity")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("A0",1.0,&MOSFET_B3::Model::a0)
    .setUnit(U_NONE)
    .setCategory(CAT_DC)
    .setDescription("Bulk charge effect coefficient for channel length")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("AGS",0.0,&MOSFET_B3::Model::ags)
    .setUnit(U_VOLTM1)
    .setCategory(CAT_DC)
    .setDescription("Gate-bias coefficient of abulk")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("A1",0.0,&MOSFET_B3::Model::a1)
    .setUnit(U_VOLTM1)
    .setCategory(CAT_DC)
    .setDescription("First non-saturation effect parameter")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("A2",1.0,&MOSFET_B3::Model::a2)
    .setUnit(U_NONE)
    .setCategory(CAT_DC)
    .setDescription("Second non-saturation factor")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("KETA",-0.047,&MOSFET_B3::Model::keta)
    .setUnit(U_VOLTM1)
    .setCategory(CAT_DC)
    .setDescription("Body-bias coefficient of bulk charge effect")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("NSUB",6.0e16,&MOSFET_B3::Model::nsub)
    .setGivenMember(&MOSFET_B3::Model::nsubGiven)
    .setUnit(U_CMM3)
    .setCategory(CAT_DOPING)
    .setDescription("Substrate doping density")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("NCH",1.7e17,&MOSFET_B3::Model::npeak)
    .setGivenMember(&MOSFET_B3::Model::npeakGiven)
    .setUnit(U_CMM3)
    .setCategory(CAT_PROCESS)
    .setDescription("Channel doping concentration")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("NGATE",0.0,&MOSFET_B3::Model::ngate)
    .setUnit(U_CMM3)
    .setCategory(CAT_DC)
    .setDescription("Poly gate doping concentration")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("GAMMA1",0.0,&MOSFET_B3::Model::gamma1)
    .setGivenMember(&MOSFET_B3::Model::gamma1Given)
    .setUnit(U_VOLTH)
    .setCategory(CAT_PROCESS)
    .setDescription("Body effect coefficient near the surface")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("GAMMA2",0.0,&MOSFET_B3::Model::gamma2)
    .setGivenMember(&MOSFET_B3::Model::gamma2Given)
    .setUnit(U_VOLTH)
    .setCategory(CAT_PROCESS)
    .setDescription("Body effect coefficient in the bulk")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("VBX",0.0,&MOSFET_B3::Model::vbx)
    .setGivenMember(&MOSFET_B3::Model::vbxGiven)
    .setUnit(U_VOLT)
    .setCategory(CAT_PROCESS)
    .setDescription("Vbs at which the depetion region = XT")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("VBM",-3.0,&MOSFET_B3::Model::vbm)
    .setGivenMember(&MOSFET_B3::Model::vbmGiven)
    .setUnit(U_VOLT)
    .setCategory(CAT_DC)
    .setDescription("Maximum applied body-bias in threshold voltage calculation")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("XT",1.55e-7,&MOSFET_B3::Model::xt)
    .setGivenMember(&MOSFET_B3::Model::xtGiven)
    .setUnit(U_METER)
    .setCategory(CAT_PROCESS)
    .setDescription("Doping depth")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("K1",0.0,&MOSFET_B3::Model::k1)
    .setGivenMember(&MOSFET_B3::Model::k1Given)
    .setUnit(U_VOLTH)
    .setCategory(CAT_DC)
    .setDescription("First-order body effect coefficient")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("KT1",-0.11,&MOSFET_B3::Model::kt1)
    .setUnit(U_VOLT)
    .setCategory(CAT_TEMP)
    .setDescription("Temperature coefficient for threshold voltage")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("KT1L",0.0,&MOSFET_B3::Model::kt1l)
    .setUnit(U_VM)
    .setCategory(CAT_TEMP)
    .setDescription("Channel length dependence of the temerature coefficient for the threshold voltage")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("KT2",0.022,&MOSFET_B3::Model::kt2)
    .setUnit(U_NONE)
    .setCategory(CAT_TEMP)
    .setDescription("Body-bias coefficient fo the threshold voltage temperature effect")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("K2",0.0,&MOSFET_B3::Model::k2)
    .setGivenMember(&MOSFET_B3::Model::k2Given)
    .setUnit(U_NONE)
    .setCategory(CAT_DC)
    .setDescription("second-order body effect coefficient")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("K3",80.0,&MOSFET_B3::Model::k3)
    .setUnit(U_NONE)
    .setCategory(CAT_DC)
    .setDescription("Narrow width coefficient")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("K3B",0.0,&MOSFET_B3::Model::k3b)
    .setUnit(U_VOLTM1)
    .setCategory(CAT_DC)
    .setDescription("Body effect coefficient of K3")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("W0",2.5e-6,&MOSFET_B3::Model::w0)
    .setUnit(U_METER)
    .setCategory(CAT_DC)
    .setDescription("Narrow-width paameter")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("NLX",1.74e-7,&MOSFET_B3::Model::nlx)
    .setUnit(U_METER)
    .setCategory(CAT_DC)
    .setDescription("Lateral non-uniform doping parameter")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("DVT0",2.2,&MOSFET_B3::Model::dvt0)
    .setUnit(U_NONE)
    .setCategory(CAT_DC)
    .setDescription("First coefficient of short-channel effect effect on threshold voltage")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("DVT1",0.53,&MOSFET_B3::Model::dvt1)
    .setUnit(U_NONE)
    .setCategory(CAT_DC)
    .setDescription("Second coefficient of short-channel effect effect on threshold voltage")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("DVT2",-0.032,&MOSFET_B3::Model::dvt2)
    .setUnit(U_VOLTM1)
    .setCategory(CAT_DC)
    .setDescription("Body-bias coefficient of short-channel effect effect on threshold voltage")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("DVT0W",0.0,&MOSFET_B3::Model::dvt0w)
    .setUnit(U_METERM1)
    .setCategory(CAT_DC)
    .setDescription("First coefficient of narrow-width effect effect on threshold voltage for small channel length")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("DVT1W",5.3e6,&MOSFET_B3::Model::dvt1w)
    .setUnit(U_METERM1)
    .setCategory(CAT_DC)
    .setDescription("Second coefficient of narrow-width effect effect on threshold voltage for small channel length")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("DVT2W",-0.032,&MOSFET_B3::Model::dvt2w)
    .setUnit(U_VOLTM1)
    .setCategory(CAT_DC)
    .setDescription("Body-bias coefficient of narrow-width effect effect on threshold voltage for small channel length")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("DROUT",0.56,&MOSFET_B3::Model::drout)
    .setUnit(U_NONE)
    .setCategory(CAT_DC)
    .setDescription("L-depedance Coefficient of the DIBL correction parameter in Rout")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("DSUB",0.0,&MOSFET_B3::Model::dsub)
    .setUnit(U_NONE)
    .setCategory(CAT_DC)
    .setDescription("DIBL coefficient exponent in subthreshhold region")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("VTH0",0.0,&MOSFET_B3::Model::vth0)
    .setGivenMember(&MOSFET_B3::Model::vth0Given)
    .setUnit(U_VOLT)
    .setCategory(CAT_DC)
    .setDescription("Threshold voltage at Vbs = 0 for large L")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("UA",2.25e-9,&MOSFET_B3::Model::ua)
    .setUnit(U_MVM1)
    .setCategory(CAT_DC)
    .setDescription("First-order mobility degradation coefficient")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("UA1",4.31e-9,&MOSFET_B3::Model::ua1)
    .setUnit(U_MVM1)
    .setCategory(CAT_TEMP)
    .setDescription("Temperature coefficient for UA")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("UB",5.87e-19,&MOSFET_B3::Model::ub)
    .setUnit(U_M2VM2)
    .setCategory(CAT_DC)
    .setDescription("First-order mobility degradation coefficient")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("UB1",-7.61e-18,&MOSFET_B3::Model::ub1)
    .setUnit(U_M2VM2)
    .setCategory(CAT_TEMP)
    .setDescription("Temperature coefficient for UB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("UC",0.0,&MOSFET_B3::Model::uc)
    .setUnit(U_MVM2)
    .setCategory(CAT_DC)
    .setDescription("Body effect of mobility degridation coefficient")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("UC1",0.0,&MOSFET_B3::Model::uc1)
    .setUnit(U_MVM2DEGCM1)
    .setCategory(CAT_TEMP)
    .setDescription("Temperature coefficient for UC")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("U0",0.0,&MOSFET_B3::Model::u0)
    .setUnit(U_CMM2VM1SM1)
    .setCategory(CAT_PROCESS)
    .setDescription("Surface mobility")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("UTE",-1.5,&MOSFET_B3::Model::ute)
    .setUnit(U_NONE)
    .setCategory(CAT_TEMP)
    .setDescription("Mobility temerature exponent")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("VOFF",-0.08,&MOSFET_B3::Model::voff)
    .setUnit(U_VOLT)
    .setCategory(CAT_DC)
    .setDescription("Offset voltage in the subthreshold region at large W and L")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("RDSW",0.0,&MOSFET_B3::Model::rdsw)
    .setUnit(U_OHMMICRON)
    .setCategory(CAT_DC)
    .setDescription("Parasitic resistance per unit width")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PRWG",0.0,&MOSFET_B3::Model::prwg)
    .setUnit(U_VOLTM1)
    .setCategory(CAT_DC)
    .setDescription("Gate-bias effect coefficient of RDSW")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PRWB",0.0,&MOSFET_B3::Model::prwb)
    .setUnit(U_VOLTMH)
    .setCategory(CAT_DC)
    .setDescription("Body effect coefficient of RDSW")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PRT",0.0,&MOSFET_B3::Model::prt)
    .setUnit(U_OHMMICRON)
    .setCategory(CAT_TEMP)
    .setDescription("Temerature coefficient for RDSW")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("ETA0",0.08,&MOSFET_B3::Model::eta0)
    .setUnit(U_NONE)
    .setCategory(CAT_DC)
    .setDescription("DIBL coefficient in subthreshold region")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("ETAB",-0.07,&MOSFET_B3::Model::etab)
    .setUnit(U_VOLTM1)
    .setCategory(CAT_DC)
    .setDescription("Body-bias coefficient for the subthreshold DIBL effect")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PCLM",1.3,&MOSFET_B3::Model::pclm)
    .setUnit(U_NONE)
    .setCategory(CAT_DC)
    .setDescription("Channel length modulation parameter")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PDIBLC1",0.39,&MOSFET_B3::Model::pdibl1)
    .setUnit(U_NONE)
    .setCategory(CAT_DC)
    .setDescription("First output resistance DIBL effect correction parameter")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PDIBLC2",0.0086,&MOSFET_B3::Model::pdibl2)
    .setUnit(U_NONE)
    .setCategory(CAT_DC)
    .setDescription("Second output resistance DIBL effect correction parameter")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PDIBLCB",0.0,&MOSFET_B3::Model::pdiblb)
    .setUnit(U_VOLTM1)
    .setCategory(CAT_DC)
    .setDescription("Body effect coefficient of DIBL correction parameter")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PSCBE1",4.24e8,&MOSFET_B3::Model::pscbe1)
    .setUnit(U_VMM1)
    .setCategory(CAT_DC)
    .setDescription("First substrate current body effect parameter")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PSCBE2",1.0e-5,&MOSFET_B3::Model::pscbe2)
    .setUnit(U_VMM1)
    .setCategory(CAT_DC)
    .setDescription("second substrate current body effect parameter")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PVAG",0.0,&MOSFET_B3::Model::pvag)
    .setUnit(U_NONE)
    .setCategory(CAT_DC)
    .setDescription("Gate dependence of early voltage")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("DELTA",0.01,&MOSFET_B3::Model::delta)
    .setUnit(U_VOLT)
    .setCategory(CAT_DC)
    .setDescription("Effective Vds parameter")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WR",1.0,&MOSFET_B3::Model::wr)
    .setUnit(U_NONE)
    .setCategory(CAT_DC)
    .setDescription("Width offset from Weff for Rds Calculation")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("DWG",0.0,&MOSFET_B3::Model::dwg)
    .setUnit(U_MVMH)
    .setCategory(CAT_DC)
    .setDescription("Coefficient of gate depedence of Weff")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("DWB",0.0,&MOSFET_B3::Model::dwb)
    .setUnit(U_MVMH)
    .setCategory(CAT_DC)
    .setDescription("Coefficient of substrate body bias dependence of Weff")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("B0",0.0,&MOSFET_B3::Model::b0)
    .setUnit(U_METER)
    .setCategory(CAT_DC)
    .setDescription("Bulk charge effect coefficient for channel width")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("B1",0.0,&MOSFET_B3::Model::b1)
    .setUnit(U_METER)
    .setCategory(CAT_DC)
    .setDescription("Bulk charge effect offset")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("ALPHA0",0.0,&MOSFET_B3::Model::alpha0)
    .setUnit(U_MVM1)
    .setCategory(CAT_DC)
    .setDescription("First parameter of impact-ionization current")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("ALPHA1",0.0,&MOSFET_B3::Model::alpha1)
    .setUnit(U_VOLTM1)
    .setCategory(CAT_DC)
    .setDescription("Isub parameter for length scaling")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("BETA0",30.0,&MOSFET_B3::Model::beta0)
    .setUnit(U_VOLT)
    .setCategory(CAT_DC)
    .setDescription("Second parameter of impact-ionization current")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("IJTH",0.1,&MOSFET_B3::Model::ijth)
    .setUnit(U_AMP)
    .setCategory(CAT_DC)
    .setDescription("Diode limiting current")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("VFB",0.0,&MOSFET_B3::Model::vfb)
    .setGivenMember(&MOSFET_B3::Model::vfbGiven)
    .setUnit(U_VOLT)
    .setCategory(CAT_DC)
    .setDescription("Flat-band voltage")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("ELM",5.0,&MOSFET_B3::Model::elm)
    .setUnit(U_NONE)
    .setCategory(CAT_NQS)
    .setDescription("Elmore constant of the channel")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("CGSL",0.0,&MOSFET_B3::Model::cgsl)
    .setUnit(U_FARADMM1)
    .setCategory(CAT_CAP)
    .setDescription("Light-doped source-gate region overlap capacitance")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("CGDL",0.0,&MOSFET_B3::Model::cgdl)
    .setUnit(U_FARADMM1)
    .setCategory(CAT_CAP)
    .setDescription("Light-doped drain-gate region overlap capacitance")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("CKAPPA",0.6,&MOSFET_B3::Model::ckappa)
    .setUnit(U_FARADMM1)
    .setCategory(CAT_CAP)
    .setDescription("Coefficient for lightly doped region overlap capacitance fireing field capacitance")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("CF",0.0,&MOSFET_B3::Model::cf)
    .setUnit(U_FARADMM1)
    .setCategory(CAT_CAP)
    .setDescription("Firing field capacitance")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("VFBCV",-1.0,&MOSFET_B3::Model::vfbcv)
    .setUnit(U_VOLT)
    .setCategory(CAT_CAP)
    .setDescription("Flat-band voltage parameter (for CAPMOD = 0 only)")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("CLC",0.1e-6,&MOSFET_B3::Model::clc)
    .setUnit(U_METER)
    .setCategory(CAT_CAP)
    .setDescription("Constant term for short-channel model")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("CLE",0.6,&MOSFET_B3::Model::cle)
    .setUnit(U_NONE)
    .setCategory(CAT_CAP)
    .setDescription("Exponetial term for the short-channel model")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("DWC",0.0,&MOSFET_B3::Model::dwc)
    .setUnit(U_METER)
    .setCategory(CAT_CAP)
    .setDescription("Width offset fitting parameter from C-V")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("DLC",0.0,&MOSFET_B3::Model::dlc)
    .setUnit(U_METER)
    .setCategory(CAT_CAP)
    .setDescription("Length offset fitting parameter from C-V")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("NOFF",1.0,&MOSFET_B3::Model::noff)
    .setUnit(U_NONE)
    .setCategory(CAT_CAP)
    .setDescription("CV parameter in Vgsteff,CV for weak to strong inversion")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("VOFFCV",0.0,&MOSFET_B3::Model::voffcv)
    .setUnit(U_VOLT)
    .setCategory(CAT_CAP)
    .setDescription("CV parameter in Vgsteff,CV for weak to strong inversion")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("ACDE",1.0,&MOSFET_B3::Model::acde)
    .setUnit(U_MVM1)
    .setCategory(CAT_CAP)
    .setDescription("Exponetial coefficient for charge thickness in capmod = 3 for accumulation and depletion regions")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("MOIN",15.0,&MOSFET_B3::Model::moin)
    .setUnit(U_NONE)
    .setCategory(CAT_CAP)
    .setDescription("Coefficient for the gate-bias dependent surface potential")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("TCJ",0.0,&MOSFET_B3::Model::tcj)
    .setUnit(U_KM1)
    .setCategory(CAT_TEMP)
    .setDescription("Temperature coefficient of Cj")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("TCJSW",0.0,&MOSFET_B3::Model::tcjsw)
    .setUnit(U_KM1)
    .setCategory(CAT_TEMP)
    .setDescription("Temperature coefficient of Cswj")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("TCJSWG",0.0,&MOSFET_B3::Model::tcjswg)
    .setUnit(U_KM1)
    .setCategory(CAT_TEMP)
    .setDescription("Temperature coefficient of Cjswg")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("TPB",0.0,&MOSFET_B3::Model::tpb)
    .setUnit(U_VKM1)
    .setCategory(CAT_TEMP)
    .setDescription("Temperature coefficient of Pb")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("TPBSW",0.0,&MOSFET_B3::Model::tpbsw)
    .setUnit(U_VKM1)
    .setCategory(CAT_TEMP)
    .setDescription("Temperature coefficient of Pbsw")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("TPBSWG",0.0,&MOSFET_B3::Model::tpbswg)
    .setUnit(U_VKM1)
    .setCategory(CAT_TEMP)
    .setDescription("Temperature coefficient of Pbswg")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LCDSC",0.0,&MOSFET_B3::Model::lcdsc)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of cdsc")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LCDSCB",0.0,&MOSFET_B3::Model::lcdscb)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of cdscb")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LCDSCD",0.0,&MOSFET_B3::Model::lcdscd)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of cdscd")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LCIT",0.0,&MOSFET_B3::Model::lcit)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of cit")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LNFACTOR",0.0,&MOSFET_B3::Model::lnfactor)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of nfactor")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LXJ",0.0,&MOSFET_B3::Model::lxj)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of xj")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LVSAT",0.0,&MOSFET_B3::Model::lvsat)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of vsat")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LAT",0.0,&MOSFET_B3::Model::lat)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of at")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LA0",0.0,&MOSFET_B3::Model::la0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of a0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LAGS",0.0,&MOSFET_B3::Model::lags)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of ags")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LA1",0.0,&MOSFET_B3::Model::la1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of a1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LA2",0.0,&MOSFET_B3::Model::la2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of a2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LKETA",0.0,&MOSFET_B3::Model::lketa)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of keta")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LNSUB",0.0,&MOSFET_B3::Model::lnsub)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of nsub")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LNCH",0.0,&MOSFET_B3::Model::lnpeak)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of nch")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LNGATE",0.0,&MOSFET_B3::Model::lngate)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of ngate")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LGAMMA1",0.0,&MOSFET_B3::Model::lgamma1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of gamma1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LGAMMA2",0.0,&MOSFET_B3::Model::lgamma2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of gamma2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LVBX",0.0,&MOSFET_B3::Model::lvbx)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of VBX")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LVBM",0.0,&MOSFET_B3::Model::lvbm)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of VBM")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LXT",0.0,&MOSFET_B3::Model::lxt)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of XT")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LK1",0.0,&MOSFET_B3::Model::lk1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of K1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LKT1",0.0,&MOSFET_B3::Model::lkt1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of KT1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LKT1L",0.0,&MOSFET_B3::Model::lkt1l)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of KT1L")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LKT2",0.0,&MOSFET_B3::Model::lkt2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of KT2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LK2",0.0,&MOSFET_B3::Model::lk2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of K2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LK3",0.0,&MOSFET_B3::Model::lk3)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of K3")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LK3B",0.0,&MOSFET_B3::Model::lk3b)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of K3B")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LW0",0.0,&MOSFET_B3::Model::lw0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of W0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LNLX",0.0,&MOSFET_B3::Model::lnlx)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of NLX")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LDVT0",0.0,&MOSFET_B3::Model::ldvt0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of DVT0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LDVT1",0.0,&MOSFET_B3::Model::ldvt1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of DVT1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LDVT2",0.0,&MOSFET_B3::Model::ldvt2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of DVT2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LDVT0W",0.0,&MOSFET_B3::Model::ldvt0w)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of DVT0W")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LDVT1W",0.0,&MOSFET_B3::Model::ldvt1w)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of DVT1W")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LDVT2W",0.0,&MOSFET_B3::Model::ldvt2w)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of DVT2W")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LDROUT",0.0,&MOSFET_B3::Model::ldrout)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of DROUT")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LDSUB",0.0,&MOSFET_B3::Model::ldsub)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of LDSUB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LVTH0",0.0,&MOSFET_B3::Model::lvth0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of VT0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LUA",0.0,&MOSFET_B3::Model::lua)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of UA")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LUA1",0.0,&MOSFET_B3::Model::lua1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of UA1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LUB",0.0,&MOSFET_B3::Model::lub)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of UB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LUB1",0.0,&MOSFET_B3::Model::lub1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of UB1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LUC",0.0,&MOSFET_B3::Model::luc)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of UC")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LUC1",0.0,&MOSFET_B3::Model::luc1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of UC1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LU0",0.0,&MOSFET_B3::Model::lu0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of U0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LUTE",0.0,&MOSFET_B3::Model::lute)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of UTE")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LVOFF",0.0,&MOSFET_B3::Model::lvoff)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of VOFF")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LRDSW",0.0,&MOSFET_B3::Model::lrdsw)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of RDSW")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LPRWG",0.0,&MOSFET_B3::Model::lprwg)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of PRWG")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LPRWB",0.0,&MOSFET_B3::Model::lprwb)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of PRWB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LPRT",0.0,&MOSFET_B3::Model::lprt)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of PRT")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LETA0",0.0,&MOSFET_B3::Model::leta0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of ETA0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LETAB",0.0,&MOSFET_B3::Model::letab)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of ETAB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LPCLM",0.0,&MOSFET_B3::Model::lpclm)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of PCLM")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LPDIBLC1",0.0,&MOSFET_B3::Model::lpdibl1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of PDIBLC1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LPDIBLC2",0.0,&MOSFET_B3::Model::lpdibl2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of PDIBLC2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LPDIBLCB",0.0,&MOSFET_B3::Model::lpdiblb)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of PDIBLCB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LPSCBE1",0.0,&MOSFET_B3::Model::lpscbe1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of PSCBE1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LPSCBE2",0.0,&MOSFET_B3::Model::lpscbe2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of PSCBE2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LPVAG",0.0,&MOSFET_B3::Model::lpvag)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of PVAG")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LDELTA",0.0,&MOSFET_B3::Model::ldelta)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of DELTA")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LWR",0.0,&MOSFET_B3::Model::lwr)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of WR")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LDWG",0.0,&MOSFET_B3::Model::ldwg)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of DWG")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LDWB",0.0,&MOSFET_B3::Model::ldwb)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of DWB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LB0",0.0,&MOSFET_B3::Model::lb0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of B0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LB1",0.0,&MOSFET_B3::Model::lb1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of B1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LALPHA0",0.0,&MOSFET_B3::Model::lalpha0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of ALPHA0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LALPHA1",0.0,&MOSFET_B3::Model::lalpha1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of ALPHA1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LBETA0",0.0,&MOSFET_B3::Model::lbeta0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of BETA0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LVFB",0.0,&MOSFET_B3::Model::lvfb)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of VFB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LELM",0.0,&MOSFET_B3::Model::lelm)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of ELM")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LCGSL",0.0,&MOSFET_B3::Model::lcgsl)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of CGSL")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LCGDL",0.0,&MOSFET_B3::Model::lcgdl)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of CGDL")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LCKAPPA",0.0,&MOSFET_B3::Model::lckappa)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of CKAPPA")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LCF",0.0,&MOSFET_B3::Model::lcf)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of CF")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LCLC",0.0,&MOSFET_B3::Model::lclc)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of CLC")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LCLE",0.0,&MOSFET_B3::Model::lcle)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of CLE")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LVFBCV",0.0,&MOSFET_B3::Model::lvfbcv)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of VFBCV")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LNOFF",0.0,&MOSFET_B3::Model::lnoff)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of NOFF")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LVOFFCV",0.0,&MOSFET_B3::Model::lvoffcv)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of VOFFCV")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LACDE",0.0,&MOSFET_B3::Model::lacde)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of ACDE")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LMOIN",0.0,&MOSFET_B3::Model::lmoin)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Length dependence of MOIN")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WCDSC",0.0,&MOSFET_B3::Model::wcdsc)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of CDSC")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WCDSCB",0.0,&MOSFET_B3::Model::wcdscb)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of CDSCB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WCDSCD",0.0,&MOSFET_B3::Model::wcdscd)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of CDSCD")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WCIT",0.0,&MOSFET_B3::Model::wcit)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of CIT")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WNFACTOR",0.0,&MOSFET_B3::Model::wnfactor)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of NFACTOR")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WXJ",0.0,&MOSFET_B3::Model::wxj)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of XJ")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WVSAT",0.0,&MOSFET_B3::Model::wvsat)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of VSAT")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WAT",0.0,&MOSFET_B3::Model::wat)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of AT")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WA0",0.0,&MOSFET_B3::Model::wa0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of A0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WAGS",0.0,&MOSFET_B3::Model::wags)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of AGS")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WA1",0.0,&MOSFET_B3::Model::wa1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of A1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WA2",0.0,&MOSFET_B3::Model::wa2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of A2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WKETA",0.0,&MOSFET_B3::Model::wketa)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of KETA")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WNSUB",0.0,&MOSFET_B3::Model::wnsub)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of NSUB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WNCH",0.0,&MOSFET_B3::Model::wnpeak)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of NCH")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WNGATE",0.0,&MOSFET_B3::Model::wngate)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of NGATE")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WGAMMA1",0.0,&MOSFET_B3::Model::wgamma1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of GAMMA1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WGAMMA2",0.0,&MOSFET_B3::Model::wgamma2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of GAMMA2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WVBX",0.0,&MOSFET_B3::Model::wvbx)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of VBX")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WVBM",0.0,&MOSFET_B3::Model::wvbm)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of VBM")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WXT",0.0,&MOSFET_B3::Model::wxt)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of XT")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WK1",0.0,&MOSFET_B3::Model::wk1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of K1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WKT1",0.0,&MOSFET_B3::Model::wkt1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of KT1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WKT1L",0.0,&MOSFET_B3::Model::wkt1l)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of KT1L")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WKT2",0.0,&MOSFET_B3::Model::wkt2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of KT2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WK2",0.0,&MOSFET_B3::Model::wk2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of K2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WK3",0.0,&MOSFET_B3::Model::wk3)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of K3")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WK3B",0.0,&MOSFET_B3::Model::wk3b)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of K3B")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WW0",0.0,&MOSFET_B3::Model::ww0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of W0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WNLX",0.0,&MOSFET_B3::Model::wnlx)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of NLX")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WDVT0",0.0,&MOSFET_B3::Model::wdvt0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of DVT0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WDVT1",0.0,&MOSFET_B3::Model::wdvt1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of DVT1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WDVT2",0.0,&MOSFET_B3::Model::wdvt2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of DVT2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WDVT0W",0.0,&MOSFET_B3::Model::wdvt0w)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of DVT0W")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WDVT1W",0.0,&MOSFET_B3::Model::wdvt1w)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of DVT1W")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WDVT2W",0.0,&MOSFET_B3::Model::wdvt2w)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of DVT2W")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WDROUT",0.0,&MOSFET_B3::Model::wdrout)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of DROUT")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WDSUB",0.0,&MOSFET_B3::Model::wdsub)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of DSUB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WVTH0",0.0,&MOSFET_B3::Model::wvth0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of VTO")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WUA",0.0,&MOSFET_B3::Model::wua)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of UA")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WUA1",0.0,&MOSFET_B3::Model::wua1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of UA1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WUB",0.0,&MOSFET_B3::Model::wub)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of UB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WUB1",0.0,&MOSFET_B3::Model::wub1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of UB1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WUC",0.0,&MOSFET_B3::Model::wuc)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of UC")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WUC1",0.0,&MOSFET_B3::Model::wuc1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of UC1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WU0",0.0,&MOSFET_B3::Model::wu0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of U0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WUTE",0.0,&MOSFET_B3::Model::wute)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of UTE")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WVOFF",0.0,&MOSFET_B3::Model::wvoff)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of VOFF")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WRDSW",0.0,&MOSFET_B3::Model::wrdsw)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of RDSW")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WPRWG",0.0,&MOSFET_B3::Model::wprwg)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of PRWG")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WPRWB",0.0,&MOSFET_B3::Model::wprwb)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of PRWB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WPRT",0.0,&MOSFET_B3::Model::wprt)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of PRT")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WETA0",0.0,&MOSFET_B3::Model::weta0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of ETA0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WETAB",0.0,&MOSFET_B3::Model::wetab)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of ETAB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WPCLM",0.0,&MOSFET_B3::Model::wpclm)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of PCLM")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WPDIBLC1",0.0,&MOSFET_B3::Model::wpdibl1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of PDIBLC1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WPDIBLC2",0.0,&MOSFET_B3::Model::wpdibl2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of PDIBLC2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WPDIBLCB",0.0,&MOSFET_B3::Model::wpdiblb)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of PDIBLCB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WPSCBE1",0.0,&MOSFET_B3::Model::wpscbe1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of PSCBE1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WPSCBE2",0.0,&MOSFET_B3::Model::wpscbe2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of PSCBE2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WPVAG",0.0,&MOSFET_B3::Model::wpvag)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of PVAG")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WDELTA",0.0,&MOSFET_B3::Model::wdelta)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of DELTA")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WWR",0.0,&MOSFET_B3::Model::wwr)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of WR")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WDWG",0.0,&MOSFET_B3::Model::wdwg)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of WG")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WDWB",0.0,&MOSFET_B3::Model::wdwb)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of DWB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WB0",0.0,&MOSFET_B3::Model::wb0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of B0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WB1",0.0,&MOSFET_B3::Model::wb1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of B1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WALPHA0",0.0,&MOSFET_B3::Model::walpha0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of ALPHA0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WALPHA1",0.0,&MOSFET_B3::Model::walpha1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of ALPHA1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WBETA0",0.0,&MOSFET_B3::Model::wbeta0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of BETA0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WVFB",0.0,&MOSFET_B3::Model::wvfb)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of VFB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WELM",0.0,&MOSFET_B3::Model::welm)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of ELM")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WCGSL",0.0,&MOSFET_B3::Model::wcgsl)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of CGSL")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WCGDL",0.0,&MOSFET_B3::Model::wcgdl)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of CGDL")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WCKAPPA",0.0,&MOSFET_B3::Model::wckappa)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of CKAPPA")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WCF",0.0,&MOSFET_B3::Model::wcf)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of CF")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WCLC",0.0,&MOSFET_B3::Model::wclc)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of CLC")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WCLE",0.0,&MOSFET_B3::Model::wcle)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of CLE")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WVFBCV",0.0,&MOSFET_B3::Model::wvfbcv)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of VFBCV")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WNOFF",0.0,&MOSFET_B3::Model::wnoff)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of NOFF")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WVOFFCV",0.0,&MOSFET_B3::Model::wvoffcv)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of VOFFCV")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WACDE",0.0,&MOSFET_B3::Model::wacde)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of ACDE")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WMOIN",0.0,&MOSFET_B3::Model::wmoin)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Width dependence of MOIN")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PCDSC",0.0,&MOSFET_B3::Model::pcdsc)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of CDSC")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PCDSCB",0.0,&MOSFET_B3::Model::pcdscb)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of CDSCB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PCDSCD",0.0,&MOSFET_B3::Model::pcdscd)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of CDSCD")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PCIT",0.0,&MOSFET_B3::Model::pcit)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of CIT")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PNFACTOR",0.0,&MOSFET_B3::Model::pnfactor)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of NFACTOR")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PXJ",0.0,&MOSFET_B3::Model::pxj)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of XJ")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PVSAT",0.0,&MOSFET_B3::Model::pvsat)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of VSAT")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PAT",0.0,&MOSFET_B3::Model::pat)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of AT")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PA0",0.0,&MOSFET_B3::Model::pa0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of A0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PAGS",0.0,&MOSFET_B3::Model::pags)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of AGS")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PA1",0.0,&MOSFET_B3::Model::pa1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of A1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PA2",0.0,&MOSFET_B3::Model::pa2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of A2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PKETA",0.0,&MOSFET_B3::Model::pketa)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of KETA")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PNSUB",0.0,&MOSFET_B3::Model::pnsub)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of NSUB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PNCH",0.0,&MOSFET_B3::Model::pnpeak)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of NCH")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PNGATE",0.0,&MOSFET_B3::Model::pngate)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of NGATE")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PGAMMA1",0.0,&MOSFET_B3::Model::pgamma1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of GAMMA1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PGAMMA2",0.0,&MOSFET_B3::Model::pgamma2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of GAMMA2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PVBX",0.0,&MOSFET_B3::Model::pvbx)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of VBX")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PVBM",0.0,&MOSFET_B3::Model::pvbm)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of VBM")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PXT",0.0,&MOSFET_B3::Model::pxt)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of XT")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PK1",0.0,&MOSFET_B3::Model::pk1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of K1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PKT1",0.0,&MOSFET_B3::Model::pkt1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of KT1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PKT1L",0.0,&MOSFET_B3::Model::pkt1l)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of KT1L")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PKT2",0.0,&MOSFET_B3::Model::pkt2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of KT2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PK2",0.0,&MOSFET_B3::Model::pk2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of K2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PK3",0.0,&MOSFET_B3::Model::pk3)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of K3")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PK3B",0.0,&MOSFET_B3::Model::pk3b)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of K3B")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PW0",0.0,&MOSFET_B3::Model::pw0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of W0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PNLX",0.0,&MOSFET_B3::Model::pnlx)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of NLX")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PDVT0",0.0,&MOSFET_B3::Model::pdvt0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of DVT0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PDVT1",0.0,&MOSFET_B3::Model::pdvt1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of DVT1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PDVT2",0.0,&MOSFET_B3::Model::pdvt2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of DVT2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PDVT0W",0.0,&MOSFET_B3::Model::pdvt0w)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of DVT0W")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PDVT1W",0.0,&MOSFET_B3::Model::pdvt1w)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of DVT1W")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PDVT2W",0.0,&MOSFET_B3::Model::pdvt2w)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of DVT2W")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PDROUT",0.0,&MOSFET_B3::Model::pdrout)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of DROUT")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PDSUB",0.0,&MOSFET_B3::Model::pdsub)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of DSUB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PVTH0",0.0,&MOSFET_B3::Model::pvth0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of VT0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PUA",0.0,&MOSFET_B3::Model::pua)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of UA")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PUA1",0.0,&MOSFET_B3::Model::pua1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of UA1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PUB",0.0,&MOSFET_B3::Model::pub)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of UB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PUB1",0.0,&MOSFET_B3::Model::pub1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of UB1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PUC",0.0,&MOSFET_B3::Model::puc)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of UC")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PUC1",0.0,&MOSFET_B3::Model::puc1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of UC1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PU0",0.0,&MOSFET_B3::Model::pu0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of U0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PUTE",0.0,&MOSFET_B3::Model::pute)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of UTE")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PVOFF",0.0,&MOSFET_B3::Model::pvoff)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of VOFF")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PRDSW",0.0,&MOSFET_B3::Model::prdsw)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of RDSW")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PPRWG",0.0,&MOSFET_B3::Model::pprwg)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of PRWG")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PPRWB",0.0,&MOSFET_B3::Model::pprwb)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of PRWB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PPRT",0.0,&MOSFET_B3::Model::pprt)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of PRT")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PETA0",0.0,&MOSFET_B3::Model::peta0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of ETA0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PETAB",0.0,&MOSFET_B3::Model::petab)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of ETAB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PPCLM",0.0,&MOSFET_B3::Model::ppclm)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of PCLM")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PPDIBLC1",0.0,&MOSFET_B3::Model::ppdibl1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of PDIBLC1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PPDIBLC2",0.0,&MOSFET_B3::Model::ppdibl2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of PDIBLC2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PPDIBLCB",0.0,&MOSFET_B3::Model::ppdiblb)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of PDIBLCB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PPSCBE1",0.0,&MOSFET_B3::Model::ppscbe1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of PSCBE1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PPSCBE2",0.0,&MOSFET_B3::Model::ppscbe2)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of PSCBE2")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PPVAG",0.0,&MOSFET_B3::Model::ppvag)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of PVAG")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PDELTA",0.0,&MOSFET_B3::Model::pdelta)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of DELTA")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PWR",0.0,&MOSFET_B3::Model::pwr)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of WR")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PDWG",0.0,&MOSFET_B3::Model::pdwg)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of DWG")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PDWB",0.0,&MOSFET_B3::Model::pdwb)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of DWB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PB0",0.0,&MOSFET_B3::Model::pb0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of B0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PB1",0.0,&MOSFET_B3::Model::pb1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of B1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PALPHA0",0.0,&MOSFET_B3::Model::palpha0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of ALPHA0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PALPHA1",0.0,&MOSFET_B3::Model::palpha1)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of ALPHA1")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PBETA0",0.0,&MOSFET_B3::Model::pbeta0)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of BETA0")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PVFB",0.0,&MOSFET_B3::Model::pvfb)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of VFB")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PELM",0.0,&MOSFET_B3::Model::pelm)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of ELM")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PCGSL",0.0,&MOSFET_B3::Model::pcgsl)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of CGSL")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PCGDL",0.0,&MOSFET_B3::Model::pcgdl)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of CGDL")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PCKAPPA",0.0,&MOSFET_B3::Model::pckappa)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of CKAPPA")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PCF",0.0,&MOSFET_B3::Model::pcf)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of CF")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PCLC",0.0,&MOSFET_B3::Model::pclc)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of CLC")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PCLE",0.0,&MOSFET_B3::Model::pcle)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of CLE")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PVFBCV",0.0,&MOSFET_B3::Model::pvfbcv)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of VFBCV")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PNOFF",0.0,&MOSFET_B3::Model::pnoff)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of NOFF")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PVOFFCV",0.0,&MOSFET_B3::Model::pvoffcv)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of VOFFCV")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PACDE",0.0,&MOSFET_B3::Model::pacde)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of ACDE")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PMOIN",0.0,&MOSFET_B3::Model::pmoin)
    .setUnit(U_INVALID)
    .setCategory(CAT_DEPENDENCY)
    .setDescription("Cross-term dependence of MOIN")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("TNOM",0.0,&MOSFET_B3::Model::tnom)
    .setUnit(STANDARD)
    .setCategory(CAT_NONE)
    .setDescription("Parameter measurement temperature")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("CGSO",0.0,&MOSFET_B3::Model::cgso)
    .setUnit(U_FARADMM1)
    .setCategory(CAT_CAP)
    .setDescription("Non-LLD region source-gate overlap capacitance per unit channel length")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("CGDO",0.0,&MOSFET_B3::Model::cgdo)
    .setUnit(U_FARADMM1)
    .setCategory(CAT_CAP)
    .setDescription("Non-LLD region drain-gate overlap capacitance per unit channel length")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("CGBO",0.0,&MOSFET_B3::Model::cgbo)
    .setUnit(U_FARADMM1)
    .setCategory(CAT_CAP)
    .setDescription("Gate-bulk overlap capacitance per unit channel length")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("XPART",0.0,&MOSFET_B3::Model::xpart)
    .setUnit(U_NONE)
    .setCategory(CAT_CAP)
    .setDescription("Charge partitioning rate flag")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("RSH",0.0,&MOSFET_B3::Model::sheetResistance)
    .setExpressionAccess(ParameterType::MIN_RES)
    .setUnit(U_OHM)
    .setCategory(CAT_RES)
    .setDescription("Drain,source diffusion sheet resistance")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("JS",1.e-4,&MOSFET_B3::Model::jctSatCurDensity)
    .setUnit(U_AMPMM2)
    .setCategory(CAT_PROCESS)
    .setDescription("Bulk p-n saturation current density")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("JSW",0.0,&MOSFET_B3::Model::jctSidewallSatCurDensity)
    .setUnit(U_AMPMM1)
    .setCategory(CAT_DC)
    .setDescription("Sidewall saturation current per unit length")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PB",1.0,&MOSFET_B3::Model::bulkJctPotential)
    .setUnit(U_VOLT)
    .setCategory(CAT_VOLT)
    .setDescription("Bulk p-n bottom potential")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("MJ",0.5,&MOSFET_B3::Model::bulkJctBotGradingCoeff)
    .setUnit(U_NONE)
    .setCategory(CAT_DOPING)
    .setDescription("Bulk p-n bottom grading coefficient")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PBSW",1.0,&MOSFET_B3::Model::sidewallJctPotential)
    .setUnit(U_VOLT)
    .setCategory(CAT_CAP)
    .setDescription("Source/drain side junction built-in potential")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("PBSWG",0.0,&MOSFET_B3::Model::GatesidewallJctPotential)
    .setUnit(U_VOLT)
    .setCategory(CAT_CAP)
    .setDescription("Source/drain gate sidewall junction built-in potential")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("MJSW",0.33,&MOSFET_B3::Model::bulkJctSideGradingCoeff)
    .setUnit(U_NONE)
    .setCategory(CAT_DOPING)
    .setDescription("Bulk p-n sidewall grading coefficient")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("CJ",5.e-4,&MOSFET_B3::Model::unitAreaJctCap)
    .setUnit(U_FARADMM2)
    .setCategory(CAT_CAP)
    .setDescription("Bulk p-n zero-bias bottom capacitance/area")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("CJSW",5.e-10,&MOSFET_B3::Model::unitLengthSidewallJctCap)
    .setUnit(U_FARADMM2)
    .setCategory(CAT_CAP)
    .setDescription("Bulk p-n zero-bias sidewall capacitance/area")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("MJSWG",0.0,&MOSFET_B3::Model::bulkJctGateSideGradingCoeff)
    .setUnit(U_NONE)
    .setCategory(CAT_CAP)
    .setDescription("Source/grain gate sidewall junction capacitance grading coeficient")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("CJSWG",0.0,&MOSFET_B3::Model::unitLengthGateSidewallJctCap)
    .setUnit(U_FARADMM1)
    .setCategory(CAT_CAP)
    .setDescription("Source/grain gate sidewall junction capacitance per unit width")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("NJ",1.0,&MOSFET_B3::Model::jctEmissionCoeff)
    .setUnit(U_NONE)
    .setCategory(CAT_TEMP)
    .setDescription("Emission coefficient of junction")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("XTI",3.0,&MOSFET_B3::Model::jctTempExponent)
    .setUnit(U_NONE)
    .setCategory(CAT_TEMP)
    .setDescription("Junction current temperature exponent coefficient")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("NOIA",0.0,&MOSFET_B3::Model::oxideTrapDensityA)
    .setUnit(U_NONE)
    .setCategory(CAT_FLICKER)
    .setDescription("Noise parameter a")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("NOIB",0.0,&MOSFET_B3::Model::oxideTrapDensityB)
    .setUnit(U_NONE)
    .setCategory(CAT_FLICKER)
    .setDescription("Noise parameter b")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("NOIC",0.0,&MOSFET_B3::Model::oxideTrapDensityC)
    .setUnit(U_NONE)
    .setCategory(CAT_FLICKER)
    .setDescription("Noise parameter c")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("EM",4.1e7,&MOSFET_B3::Model::em)
    .setUnit(U_VMM1)
    .setCategory(CAT_FLICKER)
    .setDescription("Saturation field")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("EF",1.0,&MOSFET_B3::Model::ef)
    .setUnit(U_NONE)
    .setCategory(CAT_FLICKER)
    .setDescription("Flicker exponent")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("AF",1.0,&MOSFET_B3::Model::af)
    .setUnit(U_NONE)
    .setCategory(CAT_FLICKER)
    .setDescription("Flicker noise exponent")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("KF",0.0,&MOSFET_B3::Model::kf)
    .setUnit(U_NONE)
    .setCategory(CAT_FLICKER)
    .setDescription("Flicker noise coefficient")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LINTNOI",0.0,&MOSFET_B3::Model::lintnoi)
    .setUnit(U_METER)
    .setCategory(CAT_DC)
    .setDescription("lint offset for noise calculation")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LINT",0.0,&MOSFET_B3::Model::Lint)
    .setUnit(U_METER)
    .setCategory(CAT_DC)
    .setDescription("Length of offset fiting parameter from I-V without bias")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LL",0.0,&MOSFET_B3::Model::Ll)
    .setUnit(U_MEXPLL)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Coefficient of length dependence for length offset")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LLC",0.0,&MOSFET_B3::Model::Llc)
    .setUnit(U_MEXPLL)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Coefficient of length dependence for CV channel length offset")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LLN",0.0,&MOSFET_B3::Model::Lln)
    .setUnit(U_NONE)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Power of length dependence for length offset")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LW",0.0,&MOSFET_B3::Model::Lw)
    .setUnit(U_MEXPLW)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Coefficient of width dependence for length offset")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LWC",0.0,&MOSFET_B3::Model::Lwc)
    .setUnit(U_MEXPLW)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Coefficient of width dependence for channel length offset")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LWN",0.0,&MOSFET_B3::Model::Lwn)
    .setUnit(U_NONE)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Power of width dependence for length offset")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LWL",0.0,&MOSFET_B3::Model::Lwl)
    .setUnit(U_MEXPLLLW)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Coefficient of length and width cross term for length offset")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LWLC",0.0,&MOSFET_B3::Model::Lwlc)
    .setUnit(U_MEXPLLLW)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Coefficient of length and width dependence for CV channel length offset")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WINT",0.0,&MOSFET_B3::Model::Wint)
    .setUnit(U_METER)
    .setCategory(CAT_DC)
    .setDescription("Width-offset fitting parameter from I-V without bias")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WL",0.0,&MOSFET_B3::Model::Wl)
    .setUnit(U_MEXPWL)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Coefficient of length dependence for width offset")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WLC",0.0,&MOSFET_B3::Model::Wlc)
    .setUnit(U_MEXPWL)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Coefficient of length dependence for CV channel width offset")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WLN",0.0,&MOSFET_B3::Model::Wln)
    .setUnit(U_NONE)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Power of length dependece of width offset")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WW",0.0,&MOSFET_B3::Model::Ww)
    .setUnit(U_MEXPWW)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Coefficient of width dependence for width offset")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WWC",0.0,&MOSFET_B3::Model::Wwc)
    .setUnit(U_MEXPWW)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Coefficient of width dependence for CV channel width offset")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WWN",0.0,&MOSFET_B3::Model::Wwn)
    .setUnit(U_NONE)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Power of width dependence of width offset")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WWL",0.0,&MOSFET_B3::Model::Wwl)
    .setUnit(U_MEXPWLWW)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Coefficient of length and width cross term for width offset")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WWLC",0.0,&MOSFET_B3::Model::Wwlc)
    .setUnit(U_MEXPWLWW)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Coefficient of length and width dependence for CV channel width offset")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("L",5.e-6,&MOSFET_B3::Model::model_l)
    .setUnit(U_METER)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Channel length")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("W",5.e-6,&MOSFET_B3::Model::model_w)
    .setUnit(U_METER)
    .setCategory(CAT_GEOMETRY)
    .setDescription("Channel width")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LMAX",1.0,&MOSFET_B3::Model::Lmax)
    .setUnit(U_METER)
    .setCategory(CAT_BIN)
    .setDescription("Maximum channel length")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("LMIN",0.0,&MOSFET_B3::Model::Lmin)
    .setUnit(U_METER)
    .setCategory(CAT_BIN)
    .setDescription("Minimum channel length")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WMAX",1.0,&MOSFET_B3::Model::Wmax)
    .setUnit(U_METER)
    .setCategory(CAT_BIN)
    .setDescription("Maximum channel width")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("WMIN",0.0,&MOSFET_B3::Model::Wmin)
    .setUnit(U_METER)
    .setCategory(CAT_BIN)
    .setDescription("Minimum channel width")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&bsim3ModelSens);

  p.addPar("MOBMOD",1,&MOSFET_B3::Model::mobMod)
    .setUnit(U_NONE)
    .setCategory(CAT_CONTROL)
    .setDescription("Mobility model selector")
    .setAnalyticSensitivityAvailable(false);

  p.addPar("BINUNIT",1,&MOSFET_B3::Model::binUnit)
    .setUnit(U_NONE)
    .setCategory(CAT_CONTROL)
    .setDescription("Binning unit selector")
    .setAnalyticSensitivityAvailable(false);

  p.addPar("CAPMOD",3,&MOSFET_B3::Model::capMod)
    .setUnit(U_NONE)
    .setCategory(CAT_CONTROL)
    .setDescription("Flag for capacitance models")
    .setAnalyticSensitivityAvailable(false);

  p.addPar("PARAMCHK",0,&MOSFET_B3::Model::paramChk)
    .setUnit(U_NONE)
    .setCategory(CAT_CONTROL)
    .setDescription("Parameter value check")
    .setAnalyticSensitivityAvailable(false);

  p.addPar("NOIMOD",1,&MOSFET_B3::Model::noiMod)
    .setUnit(U_NONE)
    .setCategory(CAT_CONTROL)
    .setDescription("Flag for noise models")
    .setAnalyticSensitivityAvailable(false);

  p.addPar("VERSION",std::string("3.2.2"),&MOSFET_B3::Model::version)
    .setUnit(U_NONE)
    .setCategory(CAT_CONTROL)
    .setDescription("Version number")
    .setAnalyticSensitivityAvailable(false);
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
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
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
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  // Set any non-constant parameter defaults:
  if (!given("TEMP"))
  {
    temp = getDeviceOptions().temp.getImmutableValue<double>();
    if  (!dtempGiven)
      dtemp = 0.0;
  }
  else
  {
    dtemp = 0.0;
    if  (dtempGiven)
    {
      UserWarning(*this) << "Instance temperature specified, dtemp ignored";
    }
  }

  if (!given("L"))
    l =model_.model_l;
  if (!given("W"))
    w = model_.model_w;
  if (!given("AD"))
    drainArea = getDeviceOptions().defad;
  if (!given("AS"))
    sourceArea = getDeviceOptions().defas;


  // process source/drain series resistance
  drainConductance = model_.sheetResistance * drainSquares;

  if (drainConductance > 0.0)
    drainConductance = 1.0 / drainConductance;
  else
    drainConductance = 0.0;

  sourceConductance = model_.sheetResistance * sourceSquares;

  if (sourceConductance > 0.0)
    sourceConductance = 1.0 / sourceConductance;
  else
    sourceConductance = 0.0;

  if (given("NQSMOD"))
  {
    UserWarning(*this) << 
      "  nsqMod = 1.  Not allowed yet.  Setting to 0.";
    nqsMod = 0;
  }

  if (getDeviceOptions().verboseLevel > 0 && (l > model_.Lmax || l < model_.Lmin))
  {
    UserWarning(*this) << "Channel length out of range";
  }

  if (getDeviceOptions().verboseLevel > 0 && (w > model_.Wmax || w < model_.Wmin))
  {
    UserWarning(*this) << "Channel width out of range";
  }
  
  // If there are any time dependent parameters, set their values at for
  // the current time.

  // now set the temperature related stuff.
  updateTemperature(temp);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/00
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & IB,
  Model &model,
  const FactoryBlock &  factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_                          (model),
    dNode                             (0),
    gNode                             (0),
    sNode                             (0),
    bNode                             (0),
    dNodePrime                        (0),
    sNodePrime                        (0),
    qNode                             (0),
    ueff                              (0.0),
    thetavth                          (0.0),
    von                               (0.0),
    vdsat                             (0.0),
    cgdo                              (0.0),
    cgso                              (0.0),
    vjsm                              (0.0),
    IsEvjsm                           (0.0),
    vjdm                              (0.0),
    IsEvjdm                           (0.0),
    l                                 (getDeviceOptions().defl),
    w                                 (getDeviceOptions().defw),
    numberParallel                    (1.0),
    drainArea                         (getDeviceOptions().defad),
    sourceArea                        (getDeviceOptions().defas),
    drainSquares                      (1.0),
    sourceSquares                     (1.0),
    drainPerimeter                    (0.0),
    sourcePerimeter                   (0.0),
    sourceConductance                 (0.0),
    drainConductance                  (0.0),
    icVBS                             (0.0),
    icVDS                             (0.0),
    icVGS                             (0.0),
    OFF                               (false),
    mode                              (0),
    nqsMod                            (0),
    qinv                              (0.0),
    cd                                (0.0),
    cbs                               (0.0),
    cbd                               (0.0),
    csub                              (0.0),
    cdrain                            (0.0),
    gm                                (0.0),
    gds                               (0.0),
    gmbs                              (0.0),
    gbd                               (0.0),
    gbs                               (0.0),
    gbbs                              (0.0),
    gbgs                              (0.0),
    gbds                              (0.0),
    cggb                              (0.0),
    cgdb                              (0.0),
    cgsb                              (0.0),
    cbgb                              (0.0),
    cbdb                              (0.0),
    cbsb                              (0.0),
    cdgb                              (0.0),
    cddb                              (0.0),
    cdsb                              (0.0),
    capbd                             (0.0),
    capbs                             (0.0),
    cqgb                              (0.0),
    cqdb                              (0.0),
    cqsb                              (0.0),
    cqbb                              (0.0),
    qgate                             (0.0),
    qbulk                             (0.0),
    qdrn                              (0.0),
    gtau                              (0.0),
    gtg                               (0.0),
    gtd                               (0.0),
    gts                               (0.0),
    gtb                               (0.0),
    rds                               (0.0),
    Vgsteff                           (0.0),
    Vdseff                            (0.0),
    Abulk                             (0.0),
    AbovVgst2Vtm                      (0.0),
    limitedFlag                       (false),
    paramPtr                          (NULL),
    icVBSGiven                        (0),
    icVDSGiven                        (0),
    icVGSGiven                        (0),
    ChargeComputationNeeded           (true),
    gcbdb                             (0.0),
    gcbgb                             (0.0),
    gcbsb                             (0.0),
    gcddb                             (0.0),
    gcdgb                             (0.0),
    gcdsb                             (0.0),
    gcgdb                             (0.0),
    gcggb                             (0.0),
    gcgsb                             (0.0),
    gcsdb                             (0.0),
    gcsgb                             (0.0),
    gcssb                             (0.0),
    qgd                               (0.0),
    qgs                               (0.0),
    qgb                               (0.0),
    qgdo                              (0.0),
    qgso                              (0.0),
    qsrc                              (0.0),
    CoxWL                             (0.0),
    Cgg                               (0.0),
    Cgd                               (0.0),
    Cgb                               (0.0),
    Cdg                               (0.0),
    Cdd                               (0.0),
    Cds                               (0.0),
    Csg                               (0.0),
    Csd                               (0.0),
    Css                               (0.0),
    Csb                               (0.0),
    Cbg                               (0.0),
    Cbd                               (0.0),
    Cbb                               (0.0),
    CAPcggb                           (0.0),
    CAPcgdb                           (0.0),
    CAPcgsb                           (0.0),
    CAPcbgb                           (0.0),
    CAPcbdb                           (0.0),
    CAPcbsb                           (0.0),
    CAPcdgb                           (0.0),
    CAPcddb                           (0.0),
    CAPcdsb                           (0.0),
    CAPcsgb                           (0.0),
    CAPcsdb                           (0.0),
    CAPcssb                           (0.0),
    Qeqqd_Jdxp                        (0.0),
    Qeqqb_Jdxp                        (0.0),
    Qeqqg_Jdxp                        (0.0),
    dxpart                            (0.0),
    sxpart                            (0.0),
    ggtg                              (0.0),
    ggtd                              (0.0),
    ggts                              (0.0),
    ggtb                              (0.0),
    ddxpart_dVd                       (0.0),
    ddxpart_dVg                       (0.0),
    ddxpart_dVb                       (0.0),
    ddxpart_dVs                       (0.0),
    dsxpart_dVd                       (0.0),
    dsxpart_dVg                       (0.0),
    dsxpart_dVb                       (0.0),
    dsxpart_dVs                       (0.0),
    gbspsp                            (0.0),
    gbbdp                             (0.0),
    gbbsp                             (0.0),
    gbspg                             (0.0),
    gbspb                             (0.0),
    gbspdp                            (0.0),
    gbdpdp                            (0.0),
    gbdpg                             (0.0),
    gbdpb                             (0.0),
    gbdpsp                            (0.0),
    cdreq                             (0.0),
    ceqbd                             (0.0),
    ceqbs                             (0.0),
    cdreq_Jdxp        (0.0),
    ceqbd_Jdxp        (0.0),
    ceqbs_Jdxp        (0.0),
    Gm                                (0.0),
    Gmbs                              (0.0),
    FwdSum                            (0.0),
    RevSum                            (0.0),
    T1global                          (0.0),
    dVgst_dVg                         (0.0),
    dVgst_dVb                         (0.0),
    dVgs_eff_dVg                      (0.0),
    dDeltaPhi_dVg                     (0.0),
    dDeltaPhi_dVd                     (0.0),
    dDeltaPhi_dVb                     (0.0),
    gqdef             (0.0),
    gcqdb             (0.0),
    gcqsb             (0.0),
    gcqgb             (0.0),
    gcqbb             (0.0),
    ScalingFactor     (0.0),
    cqgate            (0.0),
    cqbulk            (0.0),
    cqdrn             (0.0),
    vtm                               (0.0),
    jctTempSatCurDensity              (0.0),
    jctSidewallTempSatCurDensity      (0.0),
    unitAreaJctCapTemp                (0.0),
    unitLengthSidewallJctCapTemp      (0.0),
    unitLengthGateSidewallJctCapTemp  (0.0),
    PhiBTemp                          (0.0),
    PhiBSWTemp                        (0.0),
    PhiBSWGTemp                       (0.0),
    temp(getDeviceOptions().temp.getImmutableValue<double>()),
    dtemp(0.0),
    dtempGiven(false),
    Vd                                (0.0),
    Vs                                (0.0),
    Vg                                (0.0),
    Vb                                (0.0),
    Vsp                               (0.0),
    Vdp                               (0.0),
    Qtotal                            (0.0),
    Vddp                              (0.0),
    Vssp                              (0.0),
    Vbsp                              (0.0),
    Vbdp                              (0.0),
    Vgsp                              (0.0),
    Vgdp                              (0.0),
    Vgb                               (0.0),
    Vdpsp                             (0.0),
    Idrain                            (0.0),
    Isource                           (0.0),
    df1dVdp                           (0.0),
    df2dVdp                           (0.0),
    df1dVsp                           (0.0),
    df2dVsp                           (0.0),
    df1dVg                            (0.0),
    df2dVg                            (0.0),
    df1dVb                            (0.0),
    df2dVb                            (0.0),
    vgb               (0.0),
    vgd               (0.0),
    cqdef             (0.0),
    cqdef_Jdxp        (0.0),
    ceqqd_Jdxp        (0.0),
    ceqqb_Jdxp        (0.0),
    ceqqg_Jdxp        (0.0),
    vbd               (0.0),
    vbs               (0.0),
    vgs               (0.0),
    vds               (0.0),
    vbd_old                           (0.0),
    vbs_old                           (0.0),
    vgs_old                           (0.0),
    vds_old                           (0.0),
    vgs_orig                          (0.0),
    vds_orig                          (0.0),
    vbs_orig                          (0.0),
    vbd_orig                          (0.0),
    vgd_orig                          (0.0),
    qb                (0.0),
    qg                (0.0),
    qd                (0.0),
    qbs               (0.0),
    qbd               (0.0),
    qcheq             (0.0),
    qcdump            (0.0),
    qdef              (0.0),
// matrix and vectors indices:
// state vector: (local indices)
    li_store_vbd             (-1),
    li_store_vbs             (-1),
    li_store_vgs             (-1),
    li_store_vds             (-1),
    li_store_von             (-1),
    li_store_gm              (-1),
    li_store_Vds             (-1),
    li_store_Vgs             (-1),
    li_store_Vbs             (-1),
    li_store_Vdsat           (-1),
    li_store_Vth             (-1),
    li_store_Gds     (-1),
    li_store_Cgs     (-1),
    li_store_Cgd     (-1),
    Vds(0.0),
    Vgs(0.0), 
    Vbs(0.0),
    Vdsat(0.0),
    Vth(0.0),
    li_branch_dev_id         (-1),
    li_branch_dev_ig         (-1),
    li_branch_dev_is         (-1),
    li_branch_dev_ib         (-1),
    li_state_qb              (-1),
    li_state_qg              (-1),
    li_state_qd              (-1),
    li_state_qbs             (-1),
    li_state_qbd             (-1),
    li_state_qcheq           (-1),
    li_state_qcdump          (-1),
    li_state_qdef            (-1),
// solution vector: (local indices)
    li_Drain                 (-1),
    li_Gate                  (-1),
    li_Source                (-1),
    li_Bulk                  (-1),
    li_DrainPrime            (-1),
    li_SourcePrime           (-1),
    li_Charge                (-1),
    li_Ibs                   (-1),
    li_Ids                   (-1),
    li_Igs                   (-1),
// matrix offsets:
// jacobian:
//  drain row
    ADrainEquDrainNodeOffset             (-1),
    ADrainEquDrainPrimeNodeOffset        (-1),
    ADrainEquIdsOffset                   (-1),
//  gate row
    AGateEquGateNodeOffset               (-1),
    AGateEquBulkNodeOffset               (-1),
    AGateEquDrainPrimeNodeOffset         (-1),
    AGateEquSourcePrimeNodeOffset        (-1),
    AGateEquChargeVarOffset              (-1),
    AGateEquIgsOffset                    (-1),
//  source row
    ASourceEquSourceNodeOffset           (-1),
    ASourceEquSourcePrimeNodeOffset      (-1),
    ASourceEquIbsOffset                  (-1),
    ASourceEquIdsOffset                  (-1),
    ASourceEquIgsOffset                  (-1),
//  bulk row
    ABulkEquGateNodeOffset               (-1),
    ABulkEquBulkNodeOffset               (-1),
    ABulkEquDrainPrimeNodeOffset         (-1),
    ABulkEquSourcePrimeNodeOffset        (-1),
    ABulkEquChargeVarOffset              (-1),
    ABulkEquIbsOffset                    (-1),
// drain' row
    ADrainPrimeEquDrainNodeOffset        (-1),
    ADrainPrimeEquGateNodeOffset         (-1),
    ADrainPrimeEquBulkNodeOffset         (-1),
    ADrainPrimeEquDrainPrimeNodeOffset   (-1),
    ADrainPrimeEquSourcePrimeNodeOffset  (-1),
    ADrainPrimeEquChargeVarOffset        (-1),
// source' row
    ASourcePrimeEquGateNodeOffset        (-1),
    ASourcePrimeEquSourceNodeOffset      (-1),
    ASourcePrimeEquBulkNodeOffset        (-1),
    ASourcePrimeEquDrainPrimeNodeOffset  (-1),
    ASourcePrimeEquSourcePrimeNodeOffset (-1),
    ASourcePrimeEquChargeVarOffset       (-1),
// Charge row
    AChargeEquChargeVarOffset            (-1),
    AChargeEquDrainPrimeNodeOffset       (-1),
    AChargeEquGateNodeOffset             (-1),
    AChargeEquSourcePrimeNodeOffset      (-1),
    AChargeEquBulkNodeOffset             (-1),
    icVBSEquVsOffset                     (-1),
    icVBSEquVbOffset                     (-1),
    icVBSEquIbsOffset                    (-1),
    icVDSEquVdOffset                     (-1),
    icVDSEquVsOffset                     (-1),
    icVDSEquIdsOffset                    (-1),
    icVGSEquVgOffset                     (-1),
    icVGSEquVsOffset                     (-1),
    icVGSEquIgsOffset                    (-1),
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
//
  // V_d Row:
    f_DrainEquDrainNodePtr(0),             // a
    f_DrainEquDrainPrimeNodePtr(0),        // b
    f_DrainEquIdsPtr(0),                   // i1

  // V_g Row:
    f_GateEquGateNodePtr(0),               // c
    f_GateEquBulkNodePtr(0),               // d
    f_GateEquDrainPrimeNodePtr(0),         // e
    f_GateEquSourcePrimeNodePtr(0),        // f
    f_GateEquChargeVarPtr(0),              // 1
    f_GateEquIgsPtr(0),                    // i2

  // V_s Row:
    f_SourceEquSourceNodePtr(0),           // g
    f_SourceEquSourcePrimeNodePtr(0),      // h
    f_SourceEquIbsPtr(0),                  // i3
    f_SourceEquIdsPtr(0),                  // i4
    f_SourceEquIgsPtr(0),                  // i5

  // V_b Row:
    f_BulkEquGateNodePtr(0),               // i
    f_BulkEquBulkNodePtr(0),               // j
    f_BulkEquDrainPrimeNodePtr(0),         // k
    f_BulkEquSourcePrimeNodePtr(0),        // l
    f_BulkEquChargeVarPtr(0),              // 2
    f_BulkEquIbsPtr(0),                    // i6

  // V_d' Row:
    f_DrainPrimeEquDrainNodePtr(0),        // m
    f_DrainPrimeEquGateNodePtr(0),         // n
    f_DrainPrimeEquBulkNodePtr(0),         // o
    f_DrainPrimeEquDrainPrimeNodePtr(0),   // p
    f_DrainPrimeEquSourcePrimeNodePtr(0),  // q
    f_DrainPrimeEquChargeVarPtr(0),        // 3

  // V_s' Row:
    f_SourcePrimeEquGateNodePtr(0),        // r
    f_SourcePrimeEquSourceNodePtr(0),      // s
    f_SourcePrimeEquBulkNodePtr(0),        // t
    f_SourcePrimeEquDrainPrimeNodePtr(0),  // u
    f_SourcePrimeEquSourcePrimeNodePtr(0), // v
    f_SourcePrimeEquChargeVarPtr(0),       // 4

  // MOSFET charge (Q) Row:
    f_ChargeEquChargeVarPtr(0),            // 5
    f_ChargeEquDrainPrimeNodePtr(0),       // 6
    f_ChargeEquGateNodePtr(0),             // 7
    f_ChargeEquSourcePrimeNodePtr(0),      // 8
    f_ChargeEquBulkNodePtr(0),             // 9

    // icVBS
    f_icVBSEquVsPtr(0),                   // i7
    f_icVBSEquVbPtr(0),                   // i8
    f_icVBSEquIbsPtr(0),                  // i9

    // icVDS
    f_icVDSEquVdPtr(0),                   // i10
    f_icVDSEquVsPtr(0),                   // i11
    f_icVDSEquIdsPtr(0),                  // i12

    // icVGS
    f_icVGSEquVgPtr(0),                   // i13
    f_icVGSEquVsPtr(0),                   // i14
    f_icVGSEquIgsPtr(0),                  // i15

  // V_d Row:
    q_DrainEquDrainNodePtr(0),             // a
    q_DrainEquDrainPrimeNodePtr(0),        // b
    q_DrainEquIdsPtr(0),                   // i1

  // V_g Row:
    q_GateEquGateNodePtr(0),               // c
    q_GateEquBulkNodePtr(0),               // d
    q_GateEquDrainPrimeNodePtr(0),         // e
    q_GateEquSourcePrimeNodePtr(0),        // f
    q_GateEquChargeVarPtr(0),              // 1
    q_GateEquIgsPtr(0),                    // i2

  // V_s Row:
    q_SourceEquSourceNodePtr(0),           // g
    q_SourceEquSourcePrimeNodePtr(0),      // h
    q_SourceEquIbsPtr(0),                  // i3
    q_SourceEquIdsPtr(0),                  // i4
    q_SourceEquIgsPtr(0),                  // i5

  // V_b Row:
    q_BulkEquGateNodePtr(0),               // i
    q_BulkEquBulkNodePtr(0),               // j
    q_BulkEquDrainPrimeNodePtr(0),         // k
    q_BulkEquSourcePrimeNodePtr(0),        // l
    q_BulkEquChargeVarPtr(0),              // 2
    q_BulkEquIbsPtr(0),                    // i6

  // V_d' Row:
    q_DrainPrimeEquDrainNodePtr(0),        // m
    q_DrainPrimeEquGateNodePtr(0),         // n
    q_DrainPrimeEquBulkNodePtr(0),         // o
    q_DrainPrimeEquDrainPrimeNodePtr(0),   // p
    q_DrainPrimeEquSourcePrimeNodePtr(0),  // q
    q_DrainPrimeEquChargeVarPtr(0),        // 3

  // V_s' Row:
    q_SourcePrimeEquGateNodePtr(0),        // r
    q_SourcePrimeEquSourceNodePtr(0),      // s
    q_SourcePrimeEquBulkNodePtr(0),        // t
    q_SourcePrimeEquDrainPrimeNodePtr(0),  // u
    q_SourcePrimeEquSourcePrimeNodePtr(0), // v
    q_SourcePrimeEquChargeVarPtr(0),       // 4

  // MOSFET charge (Q) Row:
    q_ChargeEquChargeVarPtr(0),            // 5
    q_ChargeEquDrainPrimeNodePtr(0),       // 6
    q_ChargeEquGateNodePtr(0),             // 7
    q_ChargeEquSourcePrimeNodePtr(0),      // 8
    q_ChargeEquBulkNodePtr(0),             // 9

    // icVBS
    q_icVBSEquVsPtr(0),                   // i7
    q_icVBSEquVbPtr(0),                   // i8
    q_icVBSEquIbsPtr(0),                  // i9

    // icVDS
    q_icVDSEquVdPtr(0),                   // i10
    q_icVDSEquVsPtr(0),                   // i11
    q_icVDSEquIdsPtr(0),                  // i12

    // icVGS
    q_icVGSEquVgPtr(0),                   // i13
    q_icVGSEquVsPtr(0),                   // i14
    q_icVGSEquIgsPtr(0),                  // i15
//
#endif
    updateTemperatureCalled_ (false)
{
  numIntVars   = 3;
  numExtVars   = 4;
  numStateVars = 12;
  setNumStoreVars(14);
  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 4;    // this is the space to allocate if lead current or power is needed.

  devConMap.resize(4);
  devConMap[0] = 1;
  devConMap[1] = 2;
  devConMap[2] = 1;
  devConMap[3] = 3;

  if( jacStamp.empty() )
  {
    jacStamp_DC_SC.resize(6);
    jacStamp_DC_SC[0].resize(2);
    jacStamp_DC_SC[0][0]=0;
    jacStamp_DC_SC[0][1]=4;
    jacStamp_DC_SC[1].resize(4);
    jacStamp_DC_SC[1][0]=1;
    jacStamp_DC_SC[1][1]=3;
    jacStamp_DC_SC[1][2]=4;
    jacStamp_DC_SC[1][3]=5;
    jacStamp_DC_SC[2].resize(2);
    jacStamp_DC_SC[2][0]=2;
    jacStamp_DC_SC[2][1]=5;
    jacStamp_DC_SC[3].resize(4);
    jacStamp_DC_SC[3][0]=1;
    jacStamp_DC_SC[3][1]=3;
    jacStamp_DC_SC[3][2]=4;
    jacStamp_DC_SC[3][3]=5;
    jacStamp_DC_SC[4].resize(5);
    jacStamp_DC_SC[4][0]=0;
    jacStamp_DC_SC[4][1]=1;
    jacStamp_DC_SC[4][2]=3;
    jacStamp_DC_SC[4][3]=4;
    jacStamp_DC_SC[4][4]=5;
    jacStamp_DC_SC[5].resize(5);
    jacStamp_DC_SC[5][0]=1;
    jacStamp_DC_SC[5][1]=2;
    jacStamp_DC_SC[5][2]=3;
    jacStamp_DC_SC[5][3]=4;
    jacStamp_DC_SC[5][4]=5;

    jacMap_DC_SC.clear();
    jacStampMap(jacStamp_DC_SC, jacMap_DC_SC, jacMap2_DC_SC,
                jacStamp_DC,    jacMap_DC, jacMap2_DC, 5, 2, 6);

    jacStampMap(jacStamp_DC_SC, jacMap_DC_SC, jacMap2_DC_SC,
                jacStamp_SC,    jacMap_SC, jacMap2_SC, 4, 0, 6);

    jacStampMap(jacStamp_DC, jacMap_DC, jacMap2_DC,
                jacStamp,    jacMap, jacMap2, 4, 0, 6);

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "Instance::Instance jacStampMap_DS_SC" << std::endl;
      for (int k = 0; k<jacMap_DC_SC.size(); ++k )
      {
        Xyce::dout() << "jacStamp_DC_SC[ " << jacMap_DC_SC[k] << " ] = { ";
        for (int q = 0; q < jacMap2_DC_SC[k].size(); ++q )
        {
          Xyce::dout() << jacStamp_DC_SC[jacMap_DC_SC[k]][jacMap2_DC_SC[k][q]] << "  ";
        }
        Xyce::dout() << "}" << std::endl;
      }

      Xyce::dout() << "Instance::Instance jacStampMap_DS" << std::endl;
      for (int k = 0; k<jacMap_DC.size(); ++k )
      {
        Xyce::dout() << "jacStamp_DC[ " << jacMap_DC[k] << " ] = { ";
        for (int q = 0; q < jacMap2_DC[k].size(); ++q )
        {
          Xyce::dout() << jacStamp_DC[jacMap_DC[k]][jacMap2_DC[k][q]] << "  ";
        }
        Xyce::dout() << "}" << std::endl;
      }

      Xyce::dout() << "Instance::Instance jacStampMap_SC" << std::endl;
      for (int k = 0; k<jacMap_SC.size(); ++k )
      {
        Xyce::dout() << "jacStamp_SC[ " << jacMap_SC[k] << " ] = { ";
        for (int q = 0; q < jacMap2_SC[k].size(); ++q )
        {
          Xyce::dout() << jacStamp_SC[jacMap_SC[k]][jacMap2_SC[k][q]] << "  ";
        }
        Xyce::dout() << "}" << std::endl;
      }

      Xyce::dout() << "Instance::Instance jacStampMap" << std::endl;
      for (int k = 0; k<jacMap.size(); ++k )
      {
        Xyce::dout() << "jacStamp[ " << jacMap[k] << " ] = { ";
        for (int q = 0; q < jacMap2[k].size(); ++q )
        {
          Xyce::dout() << jacStamp[jacMap[k]][jacMap2[k][q]] << "  ";
        }
        Xyce::dout() << "}" << std::endl;
      }
    }
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // if options scale has been set in the netlist, apply it.
  applyScale ();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

  numIntVars = 0;
  if ( sourceConductance!=0.0 ) ++numIntVars;
  if ( drainConductance!=0.0 ) ++numIntVars;
  if ( nqsMod ) ++numIntVars;

  if (icVDSGiven) ++numIntVars;
  if (icVGSGiven) ++numIntVars;
  if (icVBSGiven) ++numIntVars;


  // if we need to, build the jacobian stamp and operating point
  // jacobian stamp specificly needed for this device's initial
  // conditions

  // There are 4 optional internal variables that may be included in
  // this device: charge Q, Vbs, Vds and Vgs.  We need to assign these
  // variables to columns in our jacobian depending on what's there.
  // For example if they're all required than Q is the 7th variable
  // and goes in column 6 (numbering from zero), Vbs is 8th in column
  // 7, Vds is 9th in column, Vgs is 10th in column 9.  If any
  // variable is not needed, those in higher columns shift
  // down. I.e. if Q isn't needed, then Vbs is in column 6, Vds in
  // column 7 and Vgs in column 8.

  int currentCol = 6;
  int numExtraCol = 0;
  int qCol = -1, icVBSCol = -1, icVDSCol = -1, icVGSCol = -1;
  if( nqsMod )
  {
    qCol = currentCol;
    ++currentCol;
    ++numExtraCol;
  }
  if( icVBSGiven )
  {
    icVBSCol = currentCol;
    ++currentCol;
    ++numExtraCol;
  }
  if( icVDSGiven )
  {
    icVDSCol = currentCol;
    ++currentCol;
    ++numExtraCol;
  }
  if( icVGSGiven )
  {
    icVGSCol = currentCol;
    ++currentCol;
    ++numExtraCol;
  }

  // ok now build this instance's jacobian stamp
  if( nqsMod || icVBSGiven || icVDSGiven || icVGSGiven )
  {
    // we need a special version of the jacStamp for the
    jacStampSpecial.resize(6 + numExtraCol);
    // row KCL d
    int numNonZeros = 2 + (icVDSGiven ? 1:0);
    jacStampSpecial[0].resize( numNonZeros );
    jacStampSpecial[0][0]=0;
    jacStampSpecial[0][1]=4;
    currentCol = 2;
    if( icVDSGiven )
    {
      jacStampSpecial[0][currentCol]=icVDSCol;
      ++currentCol;
    }

    // row KCL g
    numNonZeros = 4 + (nqsMod ? 1:0 ) + (icVGSGiven ? 1:0);
    jacStampSpecial[1].resize( numNonZeros );
    jacStampSpecial[1][0]=1;
    jacStampSpecial[1][1]=3;
    jacStampSpecial[1][2]=4;
    jacStampSpecial[1][3]=5;
    currentCol = 4;
    if( nqsMod )
    {
      jacStampSpecial[1][currentCol]=qCol;
      ++currentCol;
    }
    if( icVGSGiven )
    {
      jacStampSpecial[1][currentCol]=icVGSCol;
    }

    // row KCL s
    numNonZeros = 2 + (icVBSGiven ? 1:0) + (icVDSGiven ? 1:0) + (icVGSGiven ? 1:0);
    jacStampSpecial[2].resize( numNonZeros );
    jacStampSpecial[2][0]=2;
    jacStampSpecial[2][1]=5;
    currentCol = 2;
    if ( icVBSGiven )
    {
      jacStampSpecial[2][currentCol] = icVBSCol;
      ++currentCol;
    }
    if( icVDSGiven )
    {
      jacStampSpecial[2][currentCol] = icVDSCol;
      ++currentCol;
    }
    if( icVGSGiven )
    {
      jacStampSpecial[2][currentCol] = icVGSCol;
      ++currentCol;
    }

    // row KCL b
    numNonZeros = 4 + (nqsMod ? 1:0 ) + (icVBSGiven ? 1:0);
    jacStampSpecial[3].resize( numNonZeros );
    jacStampSpecial[3][0]=1;
    jacStampSpecial[3][1]=3;
    jacStampSpecial[3][2]=4;
    jacStampSpecial[3][3]=5;
    currentCol = 4;
    if( nqsMod )
    {
      jacStampSpecial[3][currentCol] = qCol;
      ++currentCol;
    }
    if( icVBSGiven )
    {
      jacStampSpecial[3][currentCol] = icVBSCol;
      ++currentCol;
    }

    // row KCL d'
    numNonZeros = 5 + (nqsMod ? 1:0 );
    jacStampSpecial[4].resize( numNonZeros );
    jacStampSpecial[4][0]=0;
    jacStampSpecial[4][1]=1;
    jacStampSpecial[4][2]=3;
    jacStampSpecial[4][3]=4;
    jacStampSpecial[4][4]=5;
    currentCol = 5;
    if( nqsMod )
    {
      jacStampSpecial[4][currentCol] = qCol;
      ++currentCol;
    }

    // row KCL s'
    numNonZeros = 5 + (nqsMod ? 1:0 );
    jacStampSpecial[5].resize( numNonZeros );
    jacStampSpecial[5][0]=1;
    jacStampSpecial[5][1]=2;
    jacStampSpecial[5][2]=3;
    jacStampSpecial[5][3]=4;
    jacStampSpecial[5][4]=5;
    currentCol = 5;
    if( nqsMod )
    {
      jacStampSpecial[5][currentCol] = qCol;
      ++currentCol;
    }

    int currentRow = 6;

    // Q row if we need it
    if( nqsMod )
    {
      // add in charge row
      jacStampSpecial[currentRow].resize( 5 );
      jacStampSpecial[currentRow][0] = 1;
      jacStampSpecial[currentRow][1] = 3;
      jacStampSpecial[currentRow][2] = 4;
      jacStampSpecial[currentRow][3] = 5;
      jacStampSpecial[currentRow][4] = qCol;
      ++currentRow;
    }

    // icVBS row if we need it
    if( icVBSGiven )
    {
      jacStampSpecial[currentRow].resize( 3 );
      jacStampSpecial[currentRow][0] = 2;
      jacStampSpecial[currentRow][1] = 3;
      jacStampSpecial[currentRow][2] = icVBSCol;
      ++currentRow;
    }

    // icVDS row if we need it
    if( icVDSGiven )
    {
      jacStampSpecial[currentRow].resize( 3 );
      jacStampSpecial[currentRow][0] = 0;
      jacStampSpecial[currentRow][1] = 2;
      jacStampSpecial[currentRow][2] = icVDSCol;
      ++currentRow;
    }

    // icVGS row if we need it
    if( icVGSGiven )
    {
      jacStampSpecial[currentRow].resize( 3 );
      jacStampSpecial[currentRow][0] = 1;
      jacStampSpecial[currentRow][1] = 2;
      jacStampSpecial[currentRow][2] = icVGSCol;
      ++currentRow;
    }

    // now we just need to merge Vd' and or Vs' nodes if that's called
    // for in this device
    if ( (drainConductance == 0.0) && (sourceConductance == 0.0) )
    {
      // temporary to hold intermediate results
      std::vector< std::vector<int> > jacStampSpecialMergedTemp;
      std::vector<int>           jacSpecialMapTemp;
      std::vector< std::vector<int> > jacSpecialMapTemp2;

      jacStampMap( jacStampSpecial,           jacSpecialMap,     jacSpecialMap2,
               jacStampSpecialMergedTemp, jacSpecialMapTemp, jacSpecialMapTemp2,
               5, 2, jacStampSpecial.size() );

      jacStampMap( jacStampSpecialMergedTemp, jacSpecialMapTemp,   jacSpecialMapTemp2,
               jacStampSpecialMerged,     jacSpecialMergedMap, jacSpecialMergedMap2,
               4, 0, jacStampSpecial.size() );

    }
    else if (drainConductance == 0.0)
    {
      jacStampMap( jacStampSpecial,          jacSpecialMap,          jacSpecialMap2,
               jacStampSpecialMerged,    jacSpecialMergedMap,    jacSpecialMergedMap2,
               4, 0, jacStampSpecial.size() );

    }
    else if (sourceConductance == 0.0)
    {
      jacStampMap( jacStampSpecial,       jacSpecialMap,       jacSpecialMap2,
               jacStampSpecialMerged, jacSpecialMergedMap, jacSpecialMergedMap2,
               5, 2, jacStampSpecial.size() );
    }
    else
    {
      // no rows or columns were merged, but we need to initialize
      // jacSpecialMap and jacSpecialMap2 as these will be used to
      // index into the jacobian in registerJacLIDs()
      // copied from DeviceInstance::jacStampMap initialization

      if (jacSpecialMap.size() == 0)
      {
        jacSpecialMap.resize(jacStampSpecial.size());
        jacSpecialMap2.resize(jacStampSpecial.size());
        for (int i=0 ; i<jacStampSpecial.size() ; ++i)
        {
          jacSpecialMap[i] = i;
          jacSpecialMap2[i].resize(jacStampSpecial[i].size());
          for (int j=0 ; j<jacStampSpecial[i].size() ; ++j)
          {
            jacSpecialMap2[i][j] = j;
          }
        }
      }
    }
  }
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Instance::Instance jacStampSpecial" << std::endl;
    for (int k = 0; k<jacSpecialMap.size(); ++k )
    {
      Xyce::dout() << "jacSpecialMap[" << jacSpecialMap[k] << " ] = { ";
      for (int q = 0; q < jacSpecialMap2[k].size(); ++q )
      {
        Xyce::dout() << jacStampSpecial[jacSpecialMap[k]][jacSpecialMap2[k][q]] <<"  ";
      }
      Xyce::dout() << "}" << std::endl;
    }

    Xyce::dout() << "Instance::Instance jacStampSpecialMerged" << std::endl;
    for (int k = 0; k<jacSpecialMergedMap.size(); ++k )
    {
      Xyce::dout() << "jacSpecialMap[" << jacSpecialMergedMap[k] << " ] = { ";
      for (int q = 0; q < jacSpecialMergedMap2[k].size(); ++q )
      {
        Xyce::dout() << jacStampSpecialMerged[jacSpecialMergedMap[k]][jacSpecialMergedMap2[k][q]] <<"  ";
      }
      Xyce::dout() << "}" << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/00
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                                            const std::vector<int> & extLIDVecRef )
{
  numIntVars = 0;
  if ( sourceConductance!=0.0 ) ++numIntVars;
  if ( drainConductance!=0.0 ) ++numIntVars;
  if ( nqsMod ) ++numIntVars;
  if (icVBSGiven) ++numIntVars;
  if (icVDSGiven) ++numIntVars;
  if (icVGSGiven) ++numIntVars;

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

  if( drainConductance != 0.0 )
    li_DrainPrime = intLIDVec[intLoc++];
  else
    li_DrainPrime = li_Drain;

  if( sourceConductance != 0.0 )
    li_SourcePrime = intLIDVec[intLoc++];
  else
    li_SourcePrime = li_Source;

  if( nqsMod )
    li_Charge = intLIDVec[intLoc++];

  if( icVBSGiven )
  {
    if( li_Bulk == li_Source )
    {
      DevelFatal(*this).in("Instance::registerLIDs")
        << "Tried to specify an initial condition on V_Bulk_Source when Bulk and Source nodes are the same node";
    }
    li_Ibs = intLIDVec[intLoc++];
  }

  if( icVDSGiven )
  {
    if( li_Drain == li_Source )
    {
      DevelFatal(*this).in("Instance::registerLIDs")
        << "Tried to specify an initial condition on V_Drain_Source when Drain and Source nodes are the same node";
    }
    li_Ids = intLIDVec[intLoc++];
  }

  if( icVGSGiven )
  {
    if( li_Gate == li_Source )
    {
      DevelFatal(*this).in("Instance::registerLIDs")
        << "Tried to specify an initial condition on V_Gate_Source when Gate and Source nodes are the same node";
    }
   li_Igs = intLIDVec[intLoc++];
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "\n  local variable indices:\n";
    Xyce::dout() << "  li_Drain         = " << li_Drain << std::endl;
    Xyce::dout() << "  li_Gate          = " << li_Gate << std::endl;
    Xyce::dout() << "  li_Source        = " << li_Source << std::endl;
    Xyce::dout() << "  li_Bulk          = " << li_Bulk << std::endl;

    if (drainConductance)
      Xyce::dout() << "  li_DrainPrime  = " << li_DrainPrime << std::endl;
    if (sourceConductance)
      Xyce::dout() << "  li_SourcePrime = " << li_SourcePrime << std::endl;

    if (nqsMod)
      Xyce::dout() << "  li_Charge      = " << li_Charge << std::endl;

    if (icVBSGiven)
      Xyce::dout() << "  li_Ibs         = " << li_Ibs << std::endl;

    if (icVDSGiven)
      Xyce::dout() << "  li_Ids         = " << li_Ids << std::endl;

    if (icVGSGiven)
      Xyce::dout() << "  li_Igs         = " << li_Igs << std::endl;
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
  if (drainConductance != 0.0)
    addInternalNode(symbol_table, li_DrainPrime, getName(), "drainprime");

  if (sourceConductance != 0.0)
    addInternalNode(symbol_table, li_SourcePrime, getName(), "sourceprime");

  if (icVDSGiven)
    addInternalNode(symbol_table, li_Ids, getName(), "branch_DS");

  if (icVGSGiven)
    addInternalNode(symbol_table, li_Igs, getName(), "branch_GS");

  if (icVBSGiven)
    addInternalNode(symbol_table, li_Ibs, getName(), "branch_BS");

  if (loadLeadCurrent)
  {
    addBranchDataNode(symbol_table, li_branch_dev_id, getName(), "BRANCH_DD");
    addBranchDataNode(symbol_table, li_branch_dev_ig, getName(), "BRANCH_DG");
    addBranchDataNode(symbol_table, li_branch_dev_is, getName(), "BRANCH_DS");
    addBranchDataNode(symbol_table, li_branch_dev_ib, getName(), "BRANCH_DB");
  }

  addStoreNode(symbol_table,  li_store_gm, getName().getEncodedName() + ":gm");
  addStoreNode(symbol_table,  li_store_Vds, getName().getEncodedName() + ":Vds");
  addStoreNode(symbol_table,  li_store_Vgs, getName().getEncodedName() + ":Vgs");
  addStoreNode(symbol_table,  li_store_Vbs, getName().getEncodedName() + ":Vbs");
  addStoreNode(symbol_table,  li_store_Vdsat, getName().getEncodedName() + ":Vdsat");
  addStoreNode(symbol_table,  li_store_Vth, getName().getEncodedName() + ":Vth");
  addStoreNode(symbol_table,  li_store_Gds, getName().getEncodedName() + ":Gds");
  addStoreNode(symbol_table,  li_store_Cgs, getName().getEncodedName() + ":Cgs");
  addStoreNode(symbol_table,  li_store_Cgd, getName().getEncodedName() + ":Cgd");
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
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
  // Intrinsic capacitors:
  li_state_qb        = staLIDVec[lid++];
  li_state_qg        = staLIDVec[lid++];
  li_state_qd        = staLIDVec[lid++];

  // Parasitic capacitors:
  li_state_qbs       = staLIDVec[lid++];
  li_state_qbd       = staLIDVec[lid++];

  // state variables, cheq
  li_state_qcheq     = staLIDVec[lid++];

  // state variables, cdump
  li_state_qcdump    = staLIDVec[lid++];

  // state variable, qdef
  li_state_qdef      = staLIDVec[lid++];


  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "  Local State indices:" << std::endl;
    Xyce::dout() << std::endl;
    Xyce::dout() << "  li_state_qb            = " << li_state_qb  << std::endl;
    Xyce::dout() << "  li_state_qg            = " << li_state_qg << std::endl;
    Xyce::dout() << "  li_state_qd            = " << li_state_qd << std::endl;
    Xyce::dout() << "  li_state_qbs           = " << li_state_qbs << std::endl;
    Xyce::dout() << "  li_state_qbd           = " << li_state_qbd << std::endl;
    Xyce::dout() << "  li_state_qcheq         = " << li_state_qcheq << std::endl;
    Xyce::dout() << "  li_state_qcdump        = " << li_state_qcdump << std::endl;
    Xyce::dout() << "  li_state_qdef          = " << li_state_qdef << std::endl;
    Xyce::dout() << section_divider << std::endl;
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/9/11
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs( const std::vector<int> & stoLIDVecRef )
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());

  // Copy over the global ID lists:
  stoLIDVec = stoLIDVecRef;

  int lid=0;
  // Voltage drops:
  li_store_vbd       = stoLIDVec[lid++];
  li_store_vbs       = stoLIDVec[lid++];
  li_store_vgs       = stoLIDVec[lid++];
  li_store_vds       = stoLIDVec[lid++];
  li_store_von       = stoLIDVec[lid++];

  // transconductance, and other outputs:
  li_store_gm        = stoLIDVec[lid++];
  li_store_Vds       = stoLIDVec[lid++];
  li_store_Vgs       = stoLIDVec[lid++];
  li_store_Vbs       = stoLIDVec[lid++];
  li_store_Vdsat     = stoLIDVec[lid++];
  li_store_Vth       = stoLIDVec[lid++];

  li_store_Gds   = stoLIDVec[lid++];
  li_store_Cgs   = stoLIDVec[lid++];
  li_store_Cgd   = stoLIDVec[lid++];
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
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/4/02
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  if (icVBSGiven || icVDSGiven || icVGSGiven || nqsMod )
  {
    if( !drainConductance  || !sourceConductance )
    {
      return jacStampSpecialMerged;
    }
    else
    {
      return jacStampSpecial;
    }
  }
  else
  {
    if( drainConductance && sourceConductance && !nqsMod )
      return jacStamp_DC_SC;
    else if( drainConductance && !sourceConductance && !nqsMod )
      return jacStamp_DC;
    else if( !drainConductance && sourceConductance && !nqsMod )
      return jacStamp_SC;
    else if( !drainConductance && !sourceConductance && !nqsMod )
      return jacStamp;
    else
      DevelFatal(*this).in("Instance::jacobianStamp")
        << "NQSMOD not supported for DIRECT MATRIX ACCESS\n";
  }

DevelFatal(*this).in("Instance::jacobianStamp") << 
          "should not get here!";

  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/4/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  std::vector<int> map;
  std::vector< std::vector<int> > map2;


  if (icVBSGiven || icVDSGiven || icVGSGiven || nqsMod )
  {
    if( drainConductance && sourceConductance )
    {
      map = jacSpecialMap;
      map2 = jacSpecialMap2;
    }
    else
    {
      map = jacSpecialMergedMap;
      map2 = jacSpecialMergedMap2;
    }
  }
  else
  {
    if (drainConductance)
    {
      if (sourceConductance)
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
      if (sourceConductance)
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
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "Instance::registerJacLIDs map selected" << std::endl;
      for (int k = 0; k<map.size(); ++k )
        {
          Xyce::dout() << "map[ " << k << "] = " << map[k] << "  map2[] = { ";
          for (int q = 0; q < map2[k].size(); ++q )
            {
              Xyce::dout() << map2[k][q] <<"  ";
            }
          Xyce::dout() << "}" << std::endl;
        }

      for(int k = 0; k<jacLIDVec.size(); ++k )
        {
          Xyce::dout() << "jacLIDVec[ " << k << "] = { ";
          for (int q = 0; q < jacLIDVec[k].size(); ++q )
            {
              Xyce::dout() << jacLIDVec[k][q] <<"  ";
            }
          Xyce::dout() << "}" << std::endl;
        }
    }

  int nextColumn = 0;

  // V_d row
  ADrainEquDrainNodeOffset             = jacLIDVec[map[0]][map2[0][0]];
  ADrainEquDrainPrimeNodeOffset        = jacLIDVec[map[0]][map2[0][1]];
  if( icVDSGiven )
  {
    ADrainEquIdsOffset           = jacLIDVec[map[0]][map2[0][2]];
  }

  // V_g row
  AGateEquGateNodeOffset               = jacLIDVec[map[1]][map2[1][0]];
  AGateEquBulkNodeOffset               = jacLIDVec[map[1]][map2[1][1]];
  AGateEquDrainPrimeNodeOffset         = jacLIDVec[map[1]][map2[1][2]];
  AGateEquSourcePrimeNodeOffset        = jacLIDVec[map[1]][map2[1][3]];
  nextColumn = 4;
  if( nqsMod )
  {
    AGateEquChargeVarOffset          = jacLIDVec[map[1]][map2[1][nextColumn]];
    ++nextColumn;
  }
  if( icVGSGiven )
  {
    AGateEquIgsOffset            = jacLIDVec[map[1]][map2[1][nextColumn]];
    ++nextColumn;
  }

  // V_s row
  ASourceEquSourceNodeOffset           = jacLIDVec[map[2]][map2[2][0]];
  ASourceEquSourcePrimeNodeOffset      = jacLIDVec[map[2]][map2[2][1]];
  nextColumn = 2;
  if( icVBSGiven )
  {
    ASourceEquIbsOffset          = jacLIDVec[map[2]][map2[2][nextColumn]];
    ++nextColumn;
  }
  if( icVDSGiven )
  {
    ASourceEquIdsOffset          = jacLIDVec[map[2]][map2[2][nextColumn]];
    ++nextColumn;
  }
  if( icVGSGiven )
  {
    ASourceEquIgsOffset          = jacLIDVec[map[2]][map2[2][nextColumn]];
    ++nextColumn;
  }

  // V_b row
  ABulkEquGateNodeOffset               = jacLIDVec[map[3]][map2[3][0]];
  ABulkEquBulkNodeOffset               = jacLIDVec[map[3]][map2[3][1]];
  ABulkEquDrainPrimeNodeOffset         = jacLIDVec[map[3]][map2[3][2]];
  ABulkEquSourcePrimeNodeOffset        = jacLIDVec[map[3]][map2[3][3]];
  nextColumn = 4;
  if( nqsMod )
  {
    ABulkEquChargeVarOffset          = jacLIDVec[map[3]][map2[3][nextColumn]];
    ++nextColumn;
  }
  if( icVBSGiven )
  {
    ABulkEquIbsOffset            = jacLIDVec[map[3]][map2[3][nextColumn]];
    ++nextColumn;
  }

  // V_d'
  ADrainPrimeEquDrainNodeOffset        = jacLIDVec[map[4]][map2[4][0]];
  ADrainPrimeEquGateNodeOffset         = jacLIDVec[map[4]][map2[4][1]];
  ADrainPrimeEquBulkNodeOffset         = jacLIDVec[map[4]][map2[4][2]];
  ADrainPrimeEquDrainPrimeNodeOffset   = jacLIDVec[map[4]][map2[4][3]];
  ADrainPrimeEquSourcePrimeNodeOffset  = jacLIDVec[map[4]][map2[4][4]];
  if( nqsMod )
  {
    ADrainPrimeEquChargeVarOffset    = jacLIDVec[map[4]][map2[4][5]];
  }


  // V_s'
  ASourcePrimeEquGateNodeOffset        = jacLIDVec[map[5]][map2[5][0]];
  ASourcePrimeEquSourceNodeOffset      = jacLIDVec[map[5]][map2[5][1]];
  ASourcePrimeEquBulkNodeOffset        = jacLIDVec[map[5]][map2[5][2]];
  ASourcePrimeEquDrainPrimeNodeOffset  = jacLIDVec[map[5]][map2[5][3]];
  ASourcePrimeEquSourcePrimeNodeOffset = jacLIDVec[map[5]][map2[5][4]];
  if( nqsMod )
  {
    ASourcePrimeEquChargeVarOffset   = jacLIDVec[map[5]][map2[5][5]];
  }

  int nextRow = 6;
  if( nqsMod )
  {
    AChargeEquChargeVarOffset         = jacLIDVec[map[nextRow]][map2[nextRow][0]];
    AChargeEquDrainPrimeNodeOffset    = jacLIDVec[map[nextRow]][map2[nextRow][1]];
    AChargeEquGateNodeOffset          = jacLIDVec[map[nextRow]][map2[nextRow][2]];
    AChargeEquSourcePrimeNodeOffset   = jacLIDVec[map[nextRow]][map2[nextRow][3]];
    AChargeEquBulkNodeOffset          = jacLIDVec[map[nextRow]][map2[nextRow][4]];
    ++nextRow;
  }


  if( icVBSGiven )
  {
    icVBSEquVbOffset                  = jacLIDVec[map[nextRow]][map2[nextRow][0]];
    icVBSEquVsOffset                  = jacLIDVec[map[nextRow]][map2[nextRow][1]];
    icVBSEquIbsOffset                 = jacLIDVec[map[nextRow]][map2[nextRow][2]];
    ++nextRow;
  }

  if( icVDSGiven )
  {
    icVDSEquVdOffset                  = jacLIDVec[map[nextRow]][map2[nextRow][0]];
    icVDSEquVsOffset                  = jacLIDVec[map[nextRow]][map2[nextRow][1]];
    icVDSEquIdsOffset                 = jacLIDVec[map[nextRow]][map2[nextRow][2]];
    ++nextRow;
  }

  if( icVGSGiven )
  {
    icVGSEquVgOffset                  = jacLIDVec[map[nextRow]][map2[nextRow][0]];
    icVGSEquVsOffset                  = jacLIDVec[map[nextRow]][map2[nextRow][1]];
    icVGSEquIgsOffset                 = jacLIDVec[map[nextRow]][map2[nextRow][2]];
    ++nextRow;
  }

  if (nqsMod)
  {
    DevelFatal(*this).in("Instance::registerJacLIDs")
      <<" NQSMOD not supported.";
  }
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

  // V_d row
  f_DrainEquDrainNodePtr             = 	  &(dFdx[li_Drain][ADrainEquDrainNodeOffset]);
  f_DrainEquDrainPrimeNodePtr        = 	  &(dFdx[li_Drain][ADrainEquDrainPrimeNodeOffset]);
  if( icVDSGiven )
  {
    f_DrainEquIdsPtr           = 	    &(dFdx[li_Drain][ADrainEquIdsOffset]);
  }

  // V_g row
  f_GateEquGateNodePtr               = 	  &(dFdx[li_Gate][AGateEquGateNodeOffset]);
  f_GateEquBulkNodePtr               = 	  &(dFdx[li_Gate][AGateEquBulkNodeOffset]);
  f_GateEquDrainPrimeNodePtr         = 	  &(dFdx[li_Gate][AGateEquDrainPrimeNodeOffset]);
  f_GateEquSourcePrimeNodePtr        = 	  &(dFdx[li_Gate][AGateEquSourcePrimeNodeOffset]);
  if( nqsMod )
  {
    f_GateEquChargeVarPtr          = 	    &(dFdx[li_Gate][AGateEquChargeVarOffset]);

  }
  if( icVGSGiven )
  {
    f_GateEquIgsPtr            = 	    &(dFdx[li_Gate][AGateEquIgsOffset]);

  }

  // V_s row	  // V_s row
  f_SourceEquSourceNodePtr           = 	  &(dFdx[li_Source][ASourceEquSourceNodeOffset]);
  f_SourceEquSourcePrimeNodePtr      = 	  &(dFdx[li_Source][ASourceEquSourcePrimeNodeOffset]);
  if( icVBSGiven )
  {
    f_SourceEquIbsPtr          = 	    &(dFdx[li_Source][ASourceEquIbsOffset]);

  }
  if( icVDSGiven )
  {
    f_SourceEquIdsPtr          = 	    &(dFdx[li_Source][ASourceEquIdsOffset]);

  }
  if( icVGSGiven )
  {
    f_SourceEquIgsPtr          = 	    &(dFdx[li_Source][ASourceEquIgsOffset]);

  }

  // V_b row
  f_BulkEquGateNodePtr               = 	  &(dFdx[li_Bulk][ABulkEquGateNodeOffset]);
  f_BulkEquBulkNodePtr               = 	  &(dFdx[li_Bulk][ABulkEquBulkNodeOffset]);
  f_BulkEquDrainPrimeNodePtr         = 	  &(dFdx[li_Bulk][ABulkEquDrainPrimeNodeOffset]);
  f_BulkEquSourcePrimeNodePtr        = 	  &(dFdx[li_Bulk][ABulkEquSourcePrimeNodeOffset]);
  if( nqsMod )
  {
    f_BulkEquChargeVarPtr          = 	    &(dFdx[li_Bulk][ABulkEquChargeVarOffset]);

  }
  if( icVBSGiven )
  {
    f_BulkEquIbsPtr            = 	    &(dFdx[li_Bulk][ABulkEquIbsOffset]);

  }

  // V_d'
  f_DrainPrimeEquDrainNodePtr        = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset]);
  f_DrainPrimeEquGateNodePtr         = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset]);
  f_DrainPrimeEquBulkNodePtr         = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset]);
  f_DrainPrimeEquDrainPrimeNodePtr   = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]);
  f_DrainPrimeEquSourcePrimeNodePtr  = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]);
  if( nqsMod )
  {
    f_DrainPrimeEquChargeVarPtr    = 	    &(dFdx[li_DrainPrime][ADrainPrimeEquChargeVarOffset]);
  }


  // V_s'
  f_SourcePrimeEquGateNodePtr        = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset]);
  f_SourcePrimeEquSourceNodePtr      = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset]);
  f_SourcePrimeEquBulkNodePtr        = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset]);
  f_SourcePrimeEquDrainPrimeNodePtr  = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]);
  f_SourcePrimeEquSourcePrimeNodePtr = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]);
  if( nqsMod )
  {
    f_SourcePrimeEquChargeVarPtr   = 	    &(dFdx[li_SourcePrime][ASourcePrimeEquChargeVarOffset]);
  }

  if( nqsMod )
  {
    f_ChargeEquChargeVarPtr         = 	    &(dFdx[li_SourcePrime][AChargeEquChargeVarOffset]);
    f_ChargeEquDrainPrimeNodePtr    = 	    &(dFdx[li_SourcePrime][AChargeEquDrainPrimeNodeOffset]);
    f_ChargeEquGateNodePtr          = 	    &(dFdx[li_SourcePrime][AChargeEquGateNodeOffset]);
    f_ChargeEquSourcePrimeNodePtr   = 	    &(dFdx[li_SourcePrime][AChargeEquSourcePrimeNodeOffset]);
    f_ChargeEquBulkNodePtr          = 	    &(dFdx[li_SourcePrime][AChargeEquBulkNodeOffset]);
  }


  if( icVBSGiven )
  {
    f_icVBSEquVbPtr                  = 	    &(dFdx[li_Ibs][icVBSEquVbOffset]);
    f_icVBSEquVsPtr                  = 	    &(dFdx[li_Ibs][icVBSEquVsOffset]);
    f_icVBSEquIbsPtr                 = 	    &(dFdx[li_Ibs][icVBSEquIbsOffset]);
  }

  if( icVDSGiven )
  {
    f_icVDSEquVdPtr                  = 	    &(dFdx[li_Ids][icVDSEquVdOffset]);
    f_icVDSEquVsPtr                  = 	    &(dFdx[li_Ids][icVDSEquVsOffset]);
    f_icVDSEquIdsPtr                 = 	    &(dFdx[li_Ids][icVDSEquIdsOffset]);
  }

  if( icVGSGiven )
  {
    f_icVGSEquVgPtr                  = 	    &(dFdx[li_Igs][icVGSEquVgOffset]);
    f_icVGSEquVsPtr                  = 	    &(dFdx[li_Igs][icVGSEquVsOffset]);
    f_icVGSEquIgsPtr                 = 	    &(dFdx[li_Igs][icVGSEquIgsOffset]);
  }



  // V_d row
  q_DrainEquDrainNodePtr             = 	  &(dQdx[li_Drain][ADrainEquDrainNodeOffset]);
  q_DrainEquDrainPrimeNodePtr        = 	  &(dQdx[li_Drain][ADrainEquDrainPrimeNodeOffset]);
  if( icVDSGiven )
  {
    q_DrainEquIdsPtr           = 	    &(dQdx[li_Drain][ADrainEquIdsOffset]);
  }

  // V_g row
  q_GateEquGateNodePtr               = 	  &(dQdx[li_Gate][AGateEquGateNodeOffset]);
  q_GateEquBulkNodePtr               = 	  &(dQdx[li_Gate][AGateEquBulkNodeOffset]);
  q_GateEquDrainPrimeNodePtr         = 	  &(dQdx[li_Gate][AGateEquDrainPrimeNodeOffset]);
  q_GateEquSourcePrimeNodePtr        = 	  &(dQdx[li_Gate][AGateEquSourcePrimeNodeOffset]);
  if( nqsMod )
  {
    q_GateEquChargeVarPtr          = 	    &(dQdx[li_Gate][AGateEquChargeVarOffset]);

  }
  if( icVGSGiven )
  {
    q_GateEquIgsPtr            = 	    &(dQdx[li_Gate][AGateEquIgsOffset]);

  }

  // V_s row	  // V_s row
  q_SourceEquSourceNodePtr           = 	  &(dQdx[li_Source][ASourceEquSourceNodeOffset]);
  q_SourceEquSourcePrimeNodePtr      = 	  &(dQdx[li_Source][ASourceEquSourcePrimeNodeOffset]);
  if( icVBSGiven )
  {
    q_SourceEquIbsPtr          = 	    &(dQdx[li_Source][ASourceEquIbsOffset]);

  }
  if( icVDSGiven )
  {
    q_SourceEquIdsPtr          = 	    &(dQdx[li_Source][ASourceEquIdsOffset]);

  }
  if( icVGSGiven )
  {
    q_SourceEquIgsPtr          = 	    &(dQdx[li_Source][ASourceEquIgsOffset]);

  }

  // V_b row
  q_BulkEquGateNodePtr               = 	  &(dQdx[li_Bulk][ABulkEquGateNodeOffset]);
  q_BulkEquBulkNodePtr               = 	  &(dQdx[li_Bulk][ABulkEquBulkNodeOffset]);
  q_BulkEquDrainPrimeNodePtr         = 	  &(dQdx[li_Bulk][ABulkEquDrainPrimeNodeOffset]);
  q_BulkEquSourcePrimeNodePtr        = 	  &(dQdx[li_Bulk][ABulkEquSourcePrimeNodeOffset]);
  if( nqsMod )
  {
    q_BulkEquChargeVarPtr          = 	    &(dQdx[li_Bulk][ABulkEquChargeVarOffset]);

  }
  if( icVBSGiven )
  {
    q_BulkEquIbsPtr            = 	    &(dQdx[li_Bulk][ABulkEquIbsOffset]);

  }

  // V_d'
  q_DrainPrimeEquDrainNodePtr        = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset]);
  q_DrainPrimeEquGateNodePtr         = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset]);
  q_DrainPrimeEquBulkNodePtr         = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset]);
  q_DrainPrimeEquDrainPrimeNodePtr   = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]);
  q_DrainPrimeEquSourcePrimeNodePtr  = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]);
  if( nqsMod )
  {
    q_DrainPrimeEquChargeVarPtr    = 	    &(dQdx[li_DrainPrime][ADrainPrimeEquChargeVarOffset]);
  }


  // V_s'
  q_SourcePrimeEquGateNodePtr        = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset]);
  q_SourcePrimeEquSourceNodePtr      = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset]);
  q_SourcePrimeEquBulkNodePtr        = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset]);
  q_SourcePrimeEquDrainPrimeNodePtr  = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]);
  q_SourcePrimeEquSourcePrimeNodePtr = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]);
  if( nqsMod )
  {
    q_SourcePrimeEquChargeVarPtr   = 	    &(dQdx[li_SourcePrime][ASourcePrimeEquChargeVarOffset]);
  }

  if( nqsMod )
  {
    q_ChargeEquChargeVarPtr         = 	    &(dQdx[li_SourcePrime][AChargeEquChargeVarOffset]);
    q_ChargeEquDrainPrimeNodePtr    = 	    &(dQdx[li_SourcePrime][AChargeEquDrainPrimeNodeOffset]);
    q_ChargeEquGateNodePtr          = 	    &(dQdx[li_SourcePrime][AChargeEquGateNodeOffset]);
    q_ChargeEquSourcePrimeNodePtr   = 	    &(dQdx[li_SourcePrime][AChargeEquSourcePrimeNodeOffset]);
    q_ChargeEquBulkNodePtr          = 	    &(dQdx[li_SourcePrime][AChargeEquBulkNodeOffset]);
  }


  if( icVBSGiven )
  {
    q_icVBSEquVbPtr                  = 	    &(dQdx[li_Ibs][icVBSEquVbOffset]);
    q_icVBSEquVsPtr                  = 	    &(dQdx[li_Ibs][icVBSEquVsOffset]);
    q_icVBSEquIbsPtr                 = 	    &(dQdx[li_Ibs][icVBSEquIbsOffset]);
  }

  if( icVDSGiven )
  {
    q_icVDSEquVdPtr                  = 	    &(dQdx[li_Ids][icVDSEquVdOffset]);
    q_icVDSEquVsPtr                  = 	    &(dQdx[li_Ids][icVDSEquVsOffset]);
    q_icVDSEquIdsPtr                 = 	    &(dQdx[li_Ids][icVDSEquIdsOffset]);
  }

  if( icVGSGiven )
  {
    q_icVGSEquVgPtr                  = 	    &(dQdx[li_Igs][icVGSEquVgOffset]);
    q_icVGSEquVsPtr                  = 	    &(dQdx[li_Igs][icVGSEquVsOffset]);
    q_icVGSEquIgsPtr                 = 	    &(dQdx[li_Igs][icVGSEquIgsOffset]);
  }

#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       : This updates all the instance-owned paramters which
//                 are temperature dependent.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/22/00
//-----------------------------------------------------------------------------
bool Instance::updateTemperature (const double & temp_tmp)
{
  char msg[128];

  double tmp, tmp1, tmp2, tmp3, Eg;
  double T0, T1, T2, T3, T4, T5, Ldrn, Wdrn;
  double delTemp, TRatio, Inv_L, Inv_W, Inv_LW;
  //double Dw, Dl;
  double Tnom;
  double Nvtm, SourceSatCurrent, DrainSatCurrent;

  // stuff from model paramters:
  double Eg0 = model_.Eg0;
  double ni = model_.ni;
  double Vtm0 = model_.Vtm0;

  bool bsuccess = true;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << std::endl << subsection_divider << std::endl;
    Xyce::dout() << "Instance::updateTemperature\n";
    Xyce::dout() << "name = " << getName() << std::endl;
  }

  // first set the instance temperature to the new temperature:
  if (temp_tmp != -999.0) 
  {
    temp = temp_tmp;
    temp += dtemp;
  }

  Tnom = model_.tnom;
  TRatio = temp/Tnom;

  vtm = CONSTKoverQ * temp;
  Eg = CONSTEg0 - CONSTalphaEg * temp * temp / (temp + CONSTbetaEg);

  if (temp != Tnom)
  {
    T0 = Eg0 / Vtm0 - Eg/vtm + model_.jctTempExponent*log(temp/Tnom);
    T1 = exp(T0 / model_.jctEmissionCoeff);
    jctTempSatCurDensity         = model_.jctSatCurDensity* T1;
    jctSidewallTempSatCurDensity = model_.jctSidewallSatCurDensity * T1;
  }
  else
  {
    jctTempSatCurDensity         = model_.jctSatCurDensity;
    jctSidewallTempSatCurDensity = model_.jctSidewallSatCurDensity;
  }

  if (jctTempSatCurDensity < 0.0)         jctTempSatCurDensity = 0.0;
  if (jctSidewallTempSatCurDensity < 0.0) jctSidewallTempSatCurDensity = 0.0;


  // Temperature dependence of D/B and S/B diode capacitance begins
  delTemp = temp - Tnom;
  T0 = model_.tcj * delTemp;

  if (T0 >= -1.0)
  {
    unitAreaJctCapTemp = model_.unitAreaJctCap *(1.0 + T0);
  }
  else if (unitAreaJctCapTemp > 0.0)
  {
    unitAreaJctCapTemp = 0.0;

    lout() << "Temperature effect has caused cj to be negative. Cj clamped to zero.\n" << std::endl;
  }

  T0 = model_.tcjsw * delTemp;

  if (T0 >= -1.0)
  {
    unitLengthSidewallJctCapTemp =
      model_.unitLengthSidewallJctCap *(1.0 + T0);
  }
  else if (unitLengthSidewallJctCapTemp > 0.0)
  {
    unitLengthSidewallJctCapTemp = 0.0;
    lout() << "Temperature effect has caused cjsw to be negative. Cjsw clamped to zero.\n" << std::endl;
  }

  T0 = model_.tcjswg * delTemp;

  if (T0 >= -1.0)
  {
    unitLengthGateSidewallJctCapTemp =
      model_.unitLengthGateSidewallJctCap *(1.0 + T0);
  }
  else if (unitLengthGateSidewallJctCapTemp > 0.0)
  {
    unitLengthGateSidewallJctCapTemp = 0.0;
    lout() << "Temperature effect has caused cjswg to be negative. Cjswg clamped to zero.\n" << std::endl;
  }

  PhiBTemp = model_.bulkJctPotential - model_.tpb * delTemp;

  if (PhiBTemp < 0.01)
  {
    PhiBTemp = 0.01;
    lout() << "Temperature effect has caused pb to be < 0.01. Pb clamped to 0.01.\n" << std::endl;
  }

  PhiBSWTemp = model_.sidewallJctPotential - model_.tpbsw * delTemp;

  if (PhiBSWTemp <= 0.01)
  {
    PhiBSWTemp = 0.01;
    lout() << "Temperature effect has caused pbsw to be < 0.01. Pbsw clamped to 0.01.\n" << std::endl;
  }

  PhiBSWGTemp = model_.GatesidewallJctPotential - model_.tpbswg * delTemp;

  if (PhiBSWGTemp <= 0.01)
  {
    PhiBSWGTemp = 0.01;
    lout() << "Temperature effect has caused pbswg to be < 0.01. Pbswg clamped to 0.01.\n" << std::endl;
  }
  // End of junction capacitance


  // This next block determines whether  or not to use a previously allocated
  // set of size dependent parameters. These are stored in a list that is
  // owned by the model.  If the values for length and width match those of
  // a previously allocated set, then use the old set.  If not, allocate a new set.

  std::list<SizeDependParam<double> *>::iterator it_dpL =
    model_.sizeDependParamList.begin();
  std::list<SizeDependParam<double> *>::iterator end_dpL =
    model_.sizeDependParamList.end();

  paramPtr = NULL;

  for( ; it_dpL != end_dpL; ++it_dpL )
    if( (*it_dpL)->Length == l && (*it_dpL)->Width == w  && (*it_dpL)->referenceTemperature == temp_tmp)
      paramPtr = (*it_dpL);

  if ( paramPtr != NULL )
  {
  }
  else
  {
    paramPtr = new SizeDependParam<double> ();

    model_.sizeDependParamList.push_back( paramPtr );
    paramPtr->referenceTemperature = temp_tmp;

    Ldrn = l;
    Wdrn = w;
    paramPtr->Length = Ldrn;
    paramPtr->Width = Wdrn;

    T0 = pow(Ldrn, model_.Lln);
    T1 = pow(Wdrn, model_.Lwn);

    tmp1 = model_.Ll / T0 + model_.Lw / T1
           + model_.Lwl / (T0 * T1);

    paramPtr->dl = model_.Lint + tmp1;

    tmp2 = model_.Llc / T0 + model_.Lwc / T1
           + model_.Lwlc / (T0 * T1);

    paramPtr->dlc = model_.dlc + tmp2;

    T2 = pow(Ldrn, model_.Wln);
    T3 = pow(Wdrn, model_.Wwn);

    tmp1 = model_.Wl / T2 + model_.Ww / T3
           + model_.Wwl / (T2 * T3);

    paramPtr->dw = model_.Wint + tmp1;
    tmp2 = model_.Wlc / T2 + model_.Wwc / T3
           + model_.Wwlc / (T2 * T3);

    paramPtr->dwc = model_.dwc + tmp2;

    paramPtr->leff = l - 2.0 * paramPtr->dl;
    if (paramPtr->leff <= 0.0)
    {
      UserWarning(*this) << "Effective channel length <= 0";
    }

    paramPtr->weff = w - 2.0 * paramPtr->dw;
    if (paramPtr->weff <= 0.0)
    {
      UserWarning(*this) << "Effective channel width <= 0";
    }

    paramPtr->leffCV = l - 2.0 * paramPtr->dlc;
    if (paramPtr->leffCV <= 0.0)
    {
      UserWarning(*this) << "Effective channel length for C-V <= 0";
    }

    paramPtr->weffCV = w - 2.0 * paramPtr->dwc;
    if (paramPtr->weffCV <= 0.0)
    {
      UserWarning(*this) << "Effective channel width for C-V <= 0";
    }


    if (model_.binUnit == 1)
    {
      Inv_L = 1.0e-6 / paramPtr->leff;
      Inv_W = 1.0e-6 / paramPtr->weff;
      Inv_LW = 1.0e-12 / (paramPtr->leff * paramPtr->weff);
    }
    else
    {
      Inv_L = 1.0 / paramPtr->leff;
      Inv_W = 1.0 / paramPtr->weff;
      Inv_LW = 1.0 / (paramPtr->leff * paramPtr->weff);
    }

    paramPtr->cdsc = model_.cdsc
                     + model_.lcdsc * Inv_L
                     + model_.wcdsc * Inv_W
                     + model_.pcdsc * Inv_LW;

    paramPtr->cdscb = model_.cdscb
                      + model_.lcdscb * Inv_L
                      + model_.wcdscb * Inv_W
                      + model_.pcdscb * Inv_LW;

    paramPtr->cdscd = model_.cdscd
                      + model_.lcdscd * Inv_L
                      + model_.wcdscd * Inv_W
                      + model_.pcdscd * Inv_LW;

    paramPtr->cit = model_.cit
                    + model_.lcit * Inv_L
                    + model_.wcit * Inv_W
                    + model_.pcit * Inv_LW;

    paramPtr->nfactor = model_.nfactor
                        + model_.lnfactor * Inv_L
                        + model_.wnfactor * Inv_W
                        + model_.pnfactor * Inv_LW;

    paramPtr->xj = model_.xj
                   + model_.lxj * Inv_L
                   + model_.wxj * Inv_W
                   + model_.pxj * Inv_LW;

    paramPtr->vsat = model_.vsat
                     + model_.lvsat * Inv_L
                     + model_.wvsat * Inv_W
                     + model_.pvsat * Inv_LW;

    paramPtr->at = model_.at
                   + model_.lat * Inv_L
                   + model_.wat * Inv_W
                   + model_.pat * Inv_LW;

    paramPtr->a0 = model_.a0
                   + model_.la0 * Inv_L
                   + model_.wa0 * Inv_W
                   + model_.pa0 * Inv_LW;

    paramPtr->ags = model_.ags
                    + model_.lags * Inv_L
                    + model_.wags * Inv_W
                    + model_.pags * Inv_LW;

    paramPtr->a1 = model_.a1
                   + model_.la1 * Inv_L
                   + model_.wa1 * Inv_W
                   + model_.pa1 * Inv_LW;

    paramPtr->a2 = model_.a2
                   + model_.la2 * Inv_L
                   + model_.wa2 * Inv_W
                   + model_.pa2 * Inv_LW;

    paramPtr->keta = model_.keta
                     + model_.lketa * Inv_L
                     + model_.wketa * Inv_W
                     + model_.pketa * Inv_LW;

    paramPtr->nsub = model_.nsub
                     + model_.lnsub * Inv_L
                     + model_.wnsub * Inv_W
                     + model_.pnsub * Inv_LW;

    paramPtr->npeak = model_.npeak
                      + model_.lnpeak * Inv_L
                      + model_.wnpeak * Inv_W
                      + model_.pnpeak * Inv_LW;

    paramPtr->ngate = model_.ngate
                      + model_.lngate * Inv_L
                      + model_.wngate * Inv_W
                      + model_.pngate * Inv_LW;

    paramPtr->gamma1 = model_.gamma1
                       + model_.lgamma1 * Inv_L
                       + model_.wgamma1 * Inv_W
                       + model_.pgamma1 * Inv_LW;

    paramPtr->gamma2 = model_.gamma2
                       + model_.lgamma2 * Inv_L
                       + model_.wgamma2 * Inv_W
                       + model_.pgamma2 * Inv_LW;

    paramPtr->vbx = model_.vbx
                    + model_.lvbx * Inv_L
                    + model_.wvbx * Inv_W
                    + model_.pvbx * Inv_LW;

    paramPtr->vbm = model_.vbm
                    + model_.lvbm * Inv_L
                    + model_.wvbm * Inv_W
                    + model_.pvbm * Inv_LW;

    paramPtr->xt = model_.xt
                   + model_.lxt * Inv_L
                   + model_.wxt * Inv_W
                   + model_.pxt * Inv_LW;

    paramPtr->vfb = model_.vfb
                    + model_.lvfb * Inv_L
                    + model_.wvfb * Inv_W
                    + model_.pvfb * Inv_LW;

    paramPtr->k1 = model_.k1
                   + model_.lk1 * Inv_L
                   + model_.wk1 * Inv_W
                   + model_.pk1 * Inv_LW;

    paramPtr->kt1 = model_.kt1
                    + model_.lkt1 * Inv_L
                    + model_.wkt1 * Inv_W
                    + model_.pkt1 * Inv_LW;

    paramPtr->kt1l = model_.kt1l
                     + model_.lkt1l * Inv_L
                     + model_.wkt1l * Inv_W
                     + model_.pkt1l * Inv_LW;

    paramPtr->k2 = model_.k2
                   + model_.lk2 * Inv_L
                   + model_.wk2 * Inv_W
                   + model_.pk2 * Inv_LW;

    paramPtr->kt2 = model_.kt2
                    + model_.lkt2 * Inv_L
                    + model_.wkt2 * Inv_W
                    + model_.pkt2 * Inv_LW;

    paramPtr->k3 = model_.k3
                   + model_.lk3 * Inv_L
                   + model_.wk3 * Inv_W
                   + model_.pk3 * Inv_LW;

    paramPtr->k3b = model_.k3b
                    + model_.lk3b * Inv_L
                    + model_.wk3b * Inv_W
                    + model_.pk3b * Inv_LW;

    paramPtr->w0 = model_.w0
                   + model_.lw0 * Inv_L
                   + model_.ww0 * Inv_W
                   + model_.pw0 * Inv_LW;

    paramPtr->nlx = model_.nlx
                    + model_.lnlx * Inv_L
                    + model_.wnlx * Inv_W
                    + model_.pnlx * Inv_LW;

    paramPtr->dvt0 = model_.dvt0
                     + model_.ldvt0 * Inv_L
                     + model_.wdvt0 * Inv_W
                     + model_.pdvt0 * Inv_LW;

    paramPtr->dvt1 = model_.dvt1
                     + model_.ldvt1 * Inv_L
                     + model_.wdvt1 * Inv_W
                     + model_.pdvt1 * Inv_LW;

    paramPtr->dvt2 = model_.dvt2
                     + model_.ldvt2 * Inv_L
                     + model_.wdvt2 * Inv_W
                     + model_.pdvt2 * Inv_LW;

    paramPtr->dvt0w = model_.dvt0w
                      + model_.ldvt0w * Inv_L
                      + model_.wdvt0w * Inv_W
                      + model_.pdvt0w * Inv_LW;

    paramPtr->dvt1w = model_.dvt1w
                      + model_.ldvt1w * Inv_L
                      + model_.wdvt1w * Inv_W
                      + model_.pdvt1w * Inv_LW;

    paramPtr->dvt2w = model_.dvt2w
                      + model_.ldvt2w * Inv_L
                      + model_.wdvt2w * Inv_W
                      + model_.pdvt2w * Inv_LW;

    paramPtr->drout = model_.drout
                      + model_.ldrout * Inv_L
                      + model_.wdrout * Inv_W
                      + model_.pdrout * Inv_LW;

    paramPtr->dsub = model_.dsub
                     + model_.ldsub * Inv_L
                     + model_.wdsub * Inv_W
                     + model_.pdsub * Inv_LW;

    paramPtr->vth0 = model_.vth0
                     + model_.lvth0 * Inv_L
                     + model_.wvth0 * Inv_W
                     + model_.pvth0 * Inv_LW;

    paramPtr->ua = model_.ua
                   + model_.lua * Inv_L
                   + model_.wua * Inv_W
                   + model_.pua * Inv_LW;

    paramPtr->ua1 = model_.ua1
                    + model_.lua1 * Inv_L
                    + model_.wua1 * Inv_W
                    + model_.pua1 * Inv_LW;

    paramPtr->ub = model_.ub
                   + model_.lub * Inv_L
                   + model_.wub * Inv_W
                   + model_.pub * Inv_LW;

    paramPtr->ub1 = model_.ub1
                    + model_.lub1 * Inv_L
                    + model_.wub1 * Inv_W
                    + model_.pub1 * Inv_LW;

    paramPtr->uc = model_.uc
                   + model_.luc * Inv_L
                   + model_.wuc * Inv_W
                   + model_.puc * Inv_LW;

    paramPtr->uc1 = model_.uc1
                    + model_.luc1 * Inv_L
                    + model_.wuc1 * Inv_W
                    + model_.puc1 * Inv_LW;

    paramPtr->u0 = model_.u0
                   + model_.lu0 * Inv_L
                   + model_.wu0 * Inv_W
                   + model_.pu0 * Inv_LW;

    paramPtr->ute = model_.ute
                    + model_.lute * Inv_L
                    + model_.wute * Inv_W
                    + model_.pute * Inv_LW;

    paramPtr->voff = model_.voff
                     + model_.lvoff * Inv_L
                     + model_.wvoff * Inv_W
                     + model_.pvoff * Inv_LW;

    paramPtr->delta = model_.delta
                      + model_.ldelta * Inv_L
                      + model_.wdelta * Inv_W
                      + model_.pdelta * Inv_LW;

    paramPtr->rdsw = model_.rdsw
                     + model_.lrdsw * Inv_L
                     + model_.wrdsw * Inv_W
                     + model_.prdsw * Inv_LW;

    paramPtr->prwg = model_.prwg
                     + model_.lprwg * Inv_L
                     + model_.wprwg * Inv_W
                     + model_.pprwg * Inv_LW;

    paramPtr->prwb = model_.prwb
                     + model_.lprwb * Inv_L
                     + model_.wprwb * Inv_W
                     + model_.pprwb * Inv_LW;

    paramPtr->prt = model_.prt
                    + model_.lprt * Inv_L
                    + model_.wprt * Inv_W
                    + model_.pprt * Inv_LW;

    paramPtr->eta0 = model_.eta0
                     + model_.leta0 * Inv_L
                     + model_.weta0 * Inv_W
                     + model_.peta0 * Inv_LW;

    paramPtr->etab = model_.etab
                     + model_.letab * Inv_L
                     + model_.wetab * Inv_W
                     + model_.petab * Inv_LW;

    paramPtr->pclm = model_.pclm
                     + model_.lpclm * Inv_L
                     + model_.wpclm * Inv_W
                     + model_.ppclm * Inv_LW;

    paramPtr->pdibl1 = model_.pdibl1
                       + model_.lpdibl1 * Inv_L
                       + model_.wpdibl1 * Inv_W
                       + model_.ppdibl1 * Inv_LW;

    paramPtr->pdibl2 = model_.pdibl2
                       + model_.lpdibl2 * Inv_L
                       + model_.wpdibl2 * Inv_W
                       + model_.ppdibl2 * Inv_LW;

    paramPtr->pdiblb = model_.pdiblb
                       + model_.lpdiblb * Inv_L
                       + model_.wpdiblb * Inv_W
                       + model_.ppdiblb * Inv_LW;

    paramPtr->pscbe1 = model_.pscbe1
                       + model_.lpscbe1 * Inv_L
                       + model_.wpscbe1 * Inv_W
                       + model_.ppscbe1 * Inv_LW;

    paramPtr->pscbe2 = model_.pscbe2
                       + model_.lpscbe2 * Inv_L
                       + model_.wpscbe2 * Inv_W
                       + model_.ppscbe2 * Inv_LW;

    paramPtr->pvag = model_.pvag
                     + model_.lpvag * Inv_L
                     + model_.wpvag * Inv_W
                     + model_.ppvag * Inv_LW;

    paramPtr->wr = model_.wr
                   + model_.lwr * Inv_L
                   + model_.wwr * Inv_W
                   + model_.pwr * Inv_LW;

    paramPtr->dwg = model_.dwg
                    + model_.ldwg * Inv_L
                    + model_.wdwg * Inv_W
                    + model_.pdwg * Inv_LW;

    paramPtr->dwb = model_.dwb
                    + model_.ldwb * Inv_L
                    + model_.wdwb * Inv_W
                    + model_.pdwb * Inv_LW;

    paramPtr->b0 = model_.b0
                   + model_.lb0 * Inv_L
                   + model_.wb0 * Inv_W
                   + model_.pb0 * Inv_LW;

    paramPtr->b1 = model_.b1
                   + model_.lb1 * Inv_L
                   + model_.wb1 * Inv_W
                   + model_.pb1 * Inv_LW;

    paramPtr->alpha0 = model_.alpha0
                       + model_.lalpha0 * Inv_L
                       + model_.walpha0 * Inv_W
                       + model_.palpha0 * Inv_LW;

    paramPtr->alpha1 = model_.alpha1
                       + model_.lalpha1 * Inv_L
                       + model_.walpha1 * Inv_W
                       + model_.palpha1 * Inv_LW;

    paramPtr->beta0 = model_.beta0
                      + model_.lbeta0 * Inv_L
                      + model_.wbeta0 * Inv_W
                      + model_.pbeta0 * Inv_LW;

    // CV model
    paramPtr->elm = model_.elm
                    + model_.lelm * Inv_L
                    + model_.welm * Inv_W
                    + model_.pelm * Inv_LW;

    paramPtr->cgsl = model_.cgsl
                     + model_.lcgsl * Inv_L
                     + model_.wcgsl * Inv_W
                     + model_.pcgsl * Inv_LW;

    paramPtr->cgdl = model_.cgdl
                     + model_.lcgdl * Inv_L
                     + model_.wcgdl * Inv_W
                     + model_.pcgdl * Inv_LW;

    paramPtr->ckappa = model_.ckappa
                       + model_.lckappa * Inv_L
                       + model_.wckappa * Inv_W
                       + model_.pckappa * Inv_LW;

    paramPtr->cf = model_.cf
                   + model_.lcf * Inv_L
                   + model_.wcf * Inv_W
                   + model_.pcf * Inv_LW;

    paramPtr->clc = model_.clc
                    + model_.lclc * Inv_L
                    + model_.wclc * Inv_W
                    + model_.pclc * Inv_LW;

    paramPtr->cle = model_.cle
                    + model_.lcle * Inv_L
                    + model_.wcle * Inv_W
                    + model_.pcle * Inv_LW;

    paramPtr->vfbcv = model_.vfbcv
                      + model_.lvfbcv * Inv_L
                      + model_.wvfbcv * Inv_W
                      + model_.pvfbcv * Inv_LW;

    paramPtr->acde = model_.acde
                     + model_.lacde * Inv_L
                     + model_.wacde * Inv_W
                     + model_.pacde * Inv_LW;

    paramPtr->moin = model_.moin
                     + model_.lmoin * Inv_L
                     + model_.wmoin * Inv_W
                     + model_.pmoin * Inv_LW;

    paramPtr->noff = model_.noff
                     + model_.lnoff * Inv_L
                     + model_.wnoff * Inv_W
                     + model_.pnoff * Inv_LW;

    paramPtr->voffcv = model_.voffcv
                       + model_.lvoffcv * Inv_L
                       + model_.wvoffcv * Inv_W
                       + model_.pvoffcv * Inv_LW;

    paramPtr->abulkCVfactor = 1.0
                              + pow((paramPtr->clc / paramPtr->leffCV), paramPtr->cle);

    T0 = (TRatio - 1.0);

    paramPtr->ua = paramPtr->ua + paramPtr->ua1 * T0;
    paramPtr->ub = paramPtr->ub + paramPtr->ub1 * T0;
    paramPtr->uc = paramPtr->uc + paramPtr->uc1 * T0;

    if (paramPtr->u0 > 1.0) paramPtr->u0 = paramPtr->u0 / 1.0e4;

    paramPtr->u0temp = paramPtr->u0 * pow(TRatio, paramPtr->ute);

    paramPtr->vsattemp = paramPtr->vsat - paramPtr->at * T0;

    paramPtr->rds0 = (paramPtr->rdsw + paramPtr->prt * T0)
                     / pow(paramPtr->weff * 1E6, paramPtr->wr);

      paramPtr->cgdo = (model_.cgdo + paramPtr->cf) * paramPtr->weffCV;
      paramPtr->cgso = (model_.cgso + paramPtr->cf) * paramPtr->weffCV;
      paramPtr->cgbo = model_.cgbo * paramPtr->leffCV;

      T0 = paramPtr->leffCV * paramPtr->leffCV;

      paramPtr->tconst = paramPtr->u0temp * paramPtr->elm / (model_.cox
                                                             * paramPtr->weffCV * paramPtr->leffCV * T0);

      if (!model_.npeakGiven && model_.gamma1Given)
      {
        T0 = paramPtr->gamma1 * model_.cox;
        paramPtr->npeak = 3.021E22 * T0 * T0;
      }

      paramPtr->phi     = 2.0 * Vtm0 * log(paramPtr->npeak / ni);
      paramPtr->sqrtPhi = sqrt(paramPtr->phi);
      paramPtr->phis3   = paramPtr->sqrtPhi * paramPtr->phi;

      paramPtr->Xdep0 = sqrt(2.0 * CONSTEPSSI / (CONSTQ * paramPtr->npeak * 1.0e6))
                        * paramPtr->sqrtPhi;

      paramPtr->sqrtXdep0 = sqrt(paramPtr->Xdep0);
      paramPtr->litl = sqrt(3.0 * paramPtr->xj * model_.tox);

      paramPtr->vbi = Vtm0 * log(1.0e20 * paramPtr->npeak / (ni * ni));

      paramPtr->cdep0 = sqrt(CONSTQ * CONSTEPSSI * paramPtr->npeak * 1.0e6 / 2.0
                             / paramPtr->phi);

      paramPtr->ldeb = sqrt(CONSTEPSSI * Vtm0 / (CONSTQ
                                                 * paramPtr->npeak * 1.0e6)) / 3.0;

      paramPtr->acde *= pow((paramPtr->npeak / 2.0e16), -0.25);


      if (model_.k1Given || model_.k2Given)
      {
        if (!model_.k1Given)
        {
          UserWarning(*this) << "k1 should be specified with k2.";
          paramPtr->k1 = 0.53;
        }

        if (!model_.k2Given)
        {
          UserWarning(*this) << "k2 should be specified with k1.";
          paramPtr->k2 = -0.0186;
        }

        if (model_.nsubGiven)
        {
          UserWarning(*this) << "nsub is ignored because k1 or k2 is given.";
        }

        if (model_.xtGiven)
        {
          UserWarning(*this) << "xt is ignored because k1 or k2 is given.";
        }

        if (model_.vbxGiven)
        {
          UserWarning(*this) << "vbx is ignored because k1 or k2 is given.";
        }

        if (model_.gamma1Given)
        {
          UserWarning(*this) << "gamma1 is ignored because k1 or k2 is given.";
        }

        if (model_.gamma2Given)
        {
          UserWarning(*this) << "gamma2 is ignored because k1 or k2 is given.";
        }
      }
      else
      {
        if (!model_.vbxGiven)
          paramPtr->vbx = paramPtr->phi - 7.7348e-4 * paramPtr->npeak
                          * paramPtr->xt * paramPtr->xt;

        if (paramPtr->vbx > 0.0)
          paramPtr->vbx = -paramPtr->vbx;

        if (paramPtr->vbm > 0.0)
          paramPtr->vbm = -paramPtr->vbm;

        if (!model_.gamma1Given)
          paramPtr->gamma1 = 5.753e-12 * sqrt(paramPtr->npeak) / model_.cox;

        if (!model_.gamma2Given)
          paramPtr->gamma2 = 5.753e-12 * sqrt(paramPtr->nsub) / model_.cox;

        T0 = paramPtr->gamma1 - paramPtr->gamma2;
        T1 = sqrt(paramPtr->phi - paramPtr->vbx) - paramPtr->sqrtPhi;
        T2 = sqrt(paramPtr->phi * (paramPtr->phi - paramPtr->vbm)) - paramPtr->phi;

        paramPtr->k2 = T0 * T1 / (2.0 * T2 + paramPtr->vbm);
        paramPtr->k1 = paramPtr->gamma2 - 2.0 * paramPtr->k2 * sqrt(paramPtr->phi
                                                                    - paramPtr->vbm);
      }

      if (paramPtr->k2 < 0.0)
      {
        T0 = 0.5 * paramPtr->k1 / paramPtr->k2;
        paramPtr->vbsc = 0.9 * (paramPtr->phi - T0 * T0);

        if (paramPtr->vbsc > -3.0) paramPtr->vbsc = -3.0;
        else if (paramPtr->vbsc < -30.0) paramPtr->vbsc = -30.0;
      }
      else
      {
        paramPtr->vbsc = -30.0;
      }

      if (paramPtr->vbsc > paramPtr->vbm) paramPtr->vbsc = paramPtr->vbm;

      if (!model_.vfbGiven)
      {
        if (model_.vth0Given)
        {
          paramPtr->vfb = model_.dtype * paramPtr->vth0
                          - paramPtr->phi - paramPtr->k1 * paramPtr->sqrtPhi;
        }
        else
        {   paramPtr->vfb = -1.0;
        }
      }

      if (!model_.vth0Given)
      {
        paramPtr->vth0 = model_.dtype
                         * (paramPtr->vfb + paramPtr->phi + paramPtr->k1
                            * paramPtr->sqrtPhi);
      }

      paramPtr->k1ox = paramPtr->k1 * model_.tox / model_.toxm;
      paramPtr->k2ox = paramPtr->k2 * model_.tox / model_.toxm;

      T1 = sqrt(CONSTEPSSI / CONSTEPSOX * model_.tox * paramPtr->Xdep0);
      T0 = exp(-0.5 * paramPtr->dsub * paramPtr->leff / T1);

      paramPtr->theta0vb0 = (T0 + 2.0 * T0 * T0);

      T0 = exp(-0.5 * paramPtr->drout * paramPtr->leff / T1);
      T2 = (T0 + 2.0 * T0 * T0);

      paramPtr->thetaRout = paramPtr->pdibl1 * T2 + paramPtr->pdibl2;

      tmp = sqrt(paramPtr->Xdep0);
      tmp1 = paramPtr->vbi - paramPtr->phi;
      tmp2 = model_.factor1 * tmp;

      T0 = -0.5 * paramPtr->dvt1w * paramPtr->weff * paramPtr->leff / tmp2;

      if (T0 > -CONSTEXP_THRESHOLD)
      {
        T1 = exp(T0);
        T2 = T1 * (1.0 + 2.0 * T1);
      }
      else
      {
        T1 = CONSTMIN_EXP;
        T2 = T1 * (1.0 + 2.0 * T1);
      }
      T0 = paramPtr->dvt0w * T2;
      T2 = T0 * tmp1;

      T0 = -0.5 * paramPtr->dvt1 * paramPtr->leff / tmp2;

      if (T0 > -CONSTEXP_THRESHOLD)
      {
        T1 = exp(T0);
        T3 = T1 * (1.0 + 2.0 * T1);
      }
      else
      {
        T1 = CONSTMIN_EXP;
        T3 = T1 * (1.0 + 2.0 * T1);
      }

      T3 = paramPtr->dvt0 * T3 * tmp1;

      T4 = model_.tox * paramPtr->phi / (paramPtr->weff + paramPtr->w0);

      T0 = sqrt(1.0 + paramPtr->nlx / paramPtr->leff);
      T5 = paramPtr->k1ox * (T0 - 1.0) * paramPtr->sqrtPhi
           + (paramPtr->kt1 + paramPtr->kt1l / paramPtr->leff) * (TRatio - 1.0);

      tmp3 = model_.dtype * paramPtr->vth0 - T2 - T3 + paramPtr->k3 * T4 + T5;

      paramPtr->vfbzb = tmp3 - paramPtr->phi - paramPtr->k1 * paramPtr->sqrtPhi;

    } // End of vfbzb

    cgso = paramPtr->cgso;
    cgdo = paramPtr->cgdo;

    Nvtm = vtm * model_.jctEmissionCoeff;

    if ((sourceArea <= 0.0) &&
        (sourcePerimeter <= 0.0))
    {
      SourceSatCurrent = 1.0e-14;
    }
    else
    {
      SourceSatCurrent = sourceArea * jctTempSatCurDensity
                         + sourcePerimeter
                         * jctSidewallTempSatCurDensity;
    }

    if ((SourceSatCurrent > 0.0) && (model_.ijth > 0.0))
    {
      vjsm = Nvtm * log(model_.ijth / SourceSatCurrent + 1.0);
      IsEvjsm = SourceSatCurrent * exp(vjsm / Nvtm);
    }

    if ((drainArea <= 0.0) &&
        (drainPerimeter <= 0.0))
    {
      DrainSatCurrent = 1.0e-14;
    }
    else
    {
      DrainSatCurrent = drainArea * jctTempSatCurDensity
                        + drainPerimeter
                        * jctSidewallTempSatCurDensity;
    }

    if ((DrainSatCurrent > 0.0) && (model_.ijth > 0.0))
    {
      vjdm = Nvtm * log(model_.ijth / DrainSatCurrent + 1.0);
      IsEvjdm = DrainSatCurrent * exp(vjdm / Nvtm);
    }

    updateTemperatureCalled_ = true;

    return bsuccess;
  }

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/09/01
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;

// begin the b3ld.c parameters:
  double SourceSatCurrent(0.0), DrainSatCurrent(0.0);
  double vgdo(0.0);

  double VgstNVt(0.0), ExpVgst(0.0);

  double czbd(0.0), czbdsw(0.0), czbdswg(0.0), czbs(0.0), czbssw(0.0), czbsswg(0.0);
  double evbd(0.0), evbs(0.0), arg(0.0), sarg(0.0);

  double Vfbeff(0.0), dVfbeff_dVg(0.0), dVfbeff_dVb(0.0), V3(0.0), V4(0.0);

  double MJ(0.0), MJSW(0.0), MJSWG(0.0);

  double qinoi(0.0);
  //double Vds(0.0);  // made into instance variables
  //double Vgs(0.0), Vbs(0.0);
  double Vgs_eff(0.0), Vfb(0.0);
  double Phis(0.0), dPhis_dVb(0.0), sqrtPhis(0.0), dsqrtPhis_dVb(0.0);
  //double Vth(0.0);  // made into instance variable
  double dVth_dVb(0.0), dVth_dVd(0.0);
  double Vgst(0.0);

  double Nvtm(0.0);
  double Vtm(0.0);
  double n(0.0), dn_dVb(0.0), dn_dVd(0.0), voffcv(0.0), noff(0.0), dnoff_dVd(0.0), dnoff_dVb(0.0);
  double ExpArg(0.0), V0(0.0), CoxWLcen(0.0), QovCox(0.0), LINK(0.0);
  double DeltaPhi(0.0);

  double Cox(0.0), Tox(0.0), Tcen(0.0), dTcen_dVg(0.0), dTcen_dVd(0.0), dTcen_dVb(0.0);
  double Ccen(0.0), Coxeff(0.0), dCoxeff_dVg(0.0), dCoxeff_dVd(0.0), dCoxeff_dVb(0.0);
  double Denomi(0.0), dDenomi_dVg(0.0), dDenomi_dVd(0.0), dDenomi_dVb(0.0);

  double dueff_dVg(0.0), dueff_dVd(0.0), dueff_dVb(0.0);
  double Esat(0.0);

  //double Vdsat(0.0); // made into instance variable

  double EsatL(0.0), dEsatL_dVg(0.0), dEsatL_dVd(0.0), dEsatL_dVb(0.0);

  double dVdsat_dVg(0.0), dVdsat_dVb(0.0), dVdsat_dVd(0.0), Vasat(0.0), dAlphaz_dVg(0.0), dAlphaz_dVb(0.0);
  double dVasat_dVg(0.0), dVasat_dVb(0.0), dVasat_dVd(0.0), Va(0.0);

  double dVa_dVd(0.0), dVa_dVg(0.0), dVa_dVb(0.0);
  double Vbseff(0.0), dVbseff_dVb(0.0), VbseffCV(0.0), dVbseffCV_dVb(0.0);
  double Arg1(0.0);

  double One_Third_CoxWL(0.0), Two_Third_CoxWL(0.0), Alphaz(0.0);

  double T0(0.0), dT0_dVg(0.0), dT0_dVd(0.0), dT0_dVb(0.0);
  double T1(0.0), dT1_dVg(0.0), dT1_dVd(0.0), dT1_dVb(0.0);
  double T2(0.0), dT2_dVg(0.0), dT2_dVd(0.0), dT2_dVb(0.0);
  double T3(0.0), dT3_dVg(0.0), dT3_dVd(0.0), dT3_dVb(0.0);
  double T4(0.0);

  double T5(0.0);
  double T6(0.0);
  double T7(0.0);
  double T8(0.0);
  double T9(0.0);
  double T10(0.0);
  double T11(0.0), T12(0.0);

  double tmp(0.0); 

  //double Abulk(0.0);  // needs to be instance var for noise

  double dAbulk_dVb(0.0), Abulk0(0.0), dAbulk0_dVb(0.0);

  double VACLM(0.0), dVACLM_dVg(0.0), dVACLM_dVd(0.0), dVACLM_dVb(0.0);
  double VADIBL(0.0), dVADIBL_dVg(0.0), dVADIBL_dVd(0.0), dVADIBL_dVb(0.0);

  double Xdep(0.0), dXdep_dVb(0.0), lt1(0.0), dlt1_dVb(0.0), ltw(0.0), dltw_dVb(0.0);
  double Delt_vth(0.0), dDelt_vth_dVb(0.0);

  double Theta0(0.0), dTheta0_dVb(0.0);

  double TempRatio(0.0), tmp1(0.0), tmp2(0.0), tmp3(0.0), tmp4(0.0);

  double DIBL_Sft(0.0), dDIBL_Sft_dVd(0.0);

  double Lambda(0.0), dLambda_dVg(0.0);

  double a1(0.0);

  //double Vgsteff(0.0);   // needs to be an instance var for noise.
  
  double dVgsteff_dVg(0.0), dVgsteff_dVd(0.0), dVgsteff_dVb(0.0);
  //double Vdseff(0.0);  // needs to be an instance var for noise.
  double dVdseff_dVg(0.0), dVdseff_dVd(0.0), dVdseff_dVb(0.0);
  double VdseffCV(0.0), dVdseffCV_dVg(0.0), dVdseffCV_dVd(0.0), dVdseffCV_dVb(0.0);
  double diffVds(0.0);

  double dAbulk_dVg(0.0);
  double beta(0.0), dbeta_dVg(0.0), dbeta_dVd(0.0), dbeta_dVb(0.0);
  double gche(0.0), dgche_dVg(0.0), dgche_dVd(0.0), dgche_dVb(0.0);
  double fgche1(0.0), dfgche1_dVg(0.0), dfgche1_dVd(0.0), dfgche1_dVb(0.0);
  double fgche2(0.0), dfgche2_dVg(0.0), dfgche2_dVd(0.0), dfgche2_dVb(0.0);
  double Idl(0.0), dIdl_dVg(0.0), dIdl_dVd(0.0), dIdl_dVb(0.0);
  double Idsa(0.0), dIdsa_dVg(0.0), dIdsa_dVd(0.0), dIdsa_dVb(0.0);
  double Ids(0.0);

  double Gds(0.0), Gmb(0.0);
  double Isub(0.0);

  double Gbd(0.0), Gbg(0.0), Gbb(0.0);
  double VASCBE(0.0), dVASCBE_dVg(0.0), dVASCBE_dVd(0.0), dVASCBE_dVb(0.0);
  double CoxWovL(0.0);
  double Rds(0.0), dRds_dVg(0.0), dRds_dVb(0.0), WVCox(0.0), WVCoxRds(0.0);
  double Vgst2Vtm(0.0), VdsatCV(0.0);

  double dVdsatCV_dVg(0.0), dVdsatCV_dVb(0.0);
  double Leff(0.0), Weff(0.0), dWeff_dVg(0.0), dWeff_dVb(0.0);
  double AbulkCV(0.0), dAbulkCV_dVb(0.0);

  double gtau_diff(0.0), gtau_drift(0.0);
  // these shadow member varables and then are unitialized
  // when used in later calculations
  // double qcheq(0.0), cqcheq(0.0), qdef(0.0);

  double Cgg1(0.0), Cgb1(0.0), Cgd1(0.0), Cbg1(0.0), Cbb1(0.0), Cbd1(0.0);

  double Qac0(0.0), Qsub0(0.0);
  double dQac0_dVg(0.0), dQac0_dVb(0.0), dQsub0_dVg(0.0), dQsub0_dVd(0.0), dQsub0_dVb(0.0);
  double von_local(0.0);

  ScalingFactor = 1.0e-9;

  // Don't do charge computations in DC sweeps.
  if (getSolverState().tranopFlag || getSolverState().acopFlag || getSolverState().transientFlag)
  {
    ChargeComputationNeeded = true;
  }
  else
  {
    ChargeComputationNeeded = false;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "  Instance::updateIntermediateVars\n";
    Xyce::dout() << "  name = " << getName();
    Xyce::dout() << "  model name = " << model_.getName();
    Xyce::dout() <<"   dtype is " << model_.dtype << std::endl;
    Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << "  " << std::endl;
  }

  int Check = 1;

  // The first block of code in b3ld.c basically sets up, locally,
  // what the load function should use as values for the various solution
  // variables.  There is a series of IF statements which are dependent
  // upon the mode.  (transient, initializing transient, operating point,
  // small signal, etc.).  Xyce treats the operating point and transient
  // calculation in the same way, from the device's point of view, and
  // we don't support any of the other modes.  Therefore most of these
  // mode options are not here - only the transient mode stuff.

  // First get some of the needed solution variables:
  Vd     = 0.0;
  Vs     = 0.0;
  Vg     = 0.0;
  Vb     = 0.0;
  Vsp    = 0.0;
  Vdp    = 0.0;
  Qtotal = 0.0;

  Vd = (extData.nextSolVectorRawPtr)[li_Drain];
  Vg = (extData.nextSolVectorRawPtr)[li_Gate];
  Vs = (extData.nextSolVectorRawPtr)[li_Source];
  Vb = (extData.nextSolVectorRawPtr)[li_Bulk];
  Vsp = (extData.nextSolVectorRawPtr)[li_SourcePrime];
  Vdp = (extData.nextSolVectorRawPtr)[li_DrainPrime];
  if( nqsMod )
  {
    Qtotal = (extData.nextSolVectorRawPtr)[li_Charge];
  }
  else
  {
    Qtotal = 0.0;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "     Vg = " << Vg << std::endl;
    Xyce::dout() << "     Vb = " << Vb << std::endl;
    Xyce::dout() << "     Vs = " << Vs << std::endl;
    Xyce::dout() << "     Vd = " << Vd << std::endl;
    Xyce::dout() << "    Vsp = " << Vsp << std::endl;
    Xyce::dout() << "    Vdp = " << Vdp << std::endl;
    Xyce::dout() << " Qtotal = " << Qtotal << std::endl;
  }

  Vddp  = Vd   - Vdp;
  Vssp  = Vs   - Vsp;
  Vbsp  = Vb   - Vsp;
  Vbdp  = Vb   - Vdp;
  Vgsp  = Vg   - Vsp;
  Vgdp  = Vg   - Vdp;
  Vgb   = Vg   - Vb;

  Vdpsp = Vdp  - Vsp;

  // modified from b3ld:  (see lines 221-230)
  vbs  = model_.dtype * Vbsp;
  vgs  = model_.dtype * Vgsp;
  vds  = model_.dtype * Vdpsp;

  qdef = model_.dtype * Qtotal;

  vbd = vbs - vds;
  vgd = vgs - vds;

  origFlag = 1;
  limitedFlag = false;
  vgs_orig = vgs;
  vds_orig = vds;
  vbs_orig = vbs;
  vbd_orig = vbd;
  vgd_orig = vgd;

  // What follows is a block of code designed to impose some  limits,
  //  or initial conditions on the junction voltages.  Initial conditions
  //  should only be imposed on the first Newton step of an operating point.
  //
  // The first possible limit on the  junction voltages has to do with
  // limiting the percent change of junction voltages between  Newton
  // iterations.  The second has to do with avoiding extra floating point
  // operations in the event that the device has in some sense converged
  // (aka BYPASS).  Although the primary point of BYPASS is to reduce
  // neccessary work, it also seems to reduce the number of Newton iterations.
  //
  // NOTE:  We do not support BYPASS.
  //
  // The "old" variables should be the values for the previous
  // Newton iteration, if indeed there was a previous Newton
  // iteration.  If not, just set the  old values equal to
  // the current ones.
  //

  // set an initial condition if appropriate:
  if (getSolverState().initJctFlag_ && !OFF && getDeviceOptions().voltageLimiterFlag)
  {
    if (getSolverState().inputOPFlag)
    {
      Linear::Vector * flagSolVectorPtr = extData.flagSolVectorPtr;
      if ((*flagSolVectorPtr)[li_Drain] == 0 || (*flagSolVectorPtr)[li_Gate] == 0 ||
          (*flagSolVectorPtr)[li_Source] == 0 || (*flagSolVectorPtr)[li_SourcePrime] == 0 ||
          (*flagSolVectorPtr)[li_DrainPrime] == 0 || (*flagSolVectorPtr)[li_Bulk] == 0 )
      {
        vbs = 0.0;
        vgs = model_.dtype * paramPtr->vth0 + 0.1;
        vds = 0.1;
        origFlag = 0;
      }
    }
    else
    {
      vbs = 0.0;
      vgs = model_.dtype * paramPtr->vth0 + 0.1;
      vds = 0.1;
      origFlag = 0;
    }
    vbd = vbs - vds;
    vgd = vgs - vds;
    //origFlag = 0;
  }
  else if ((getSolverState().initFixFlag || getSolverState().initJctFlag_) && OFF)
  {
    qdef = vbs = vgs = vds = 0;
  }


  if (getSolverState().newtonIter == 0)
  {

    if (!getSolverState().dcopFlag || (getSolverState().locaEnabledFlag && getSolverState().dcopFlag))
      // ie, first newton step of a transient time step or DCOP continuation step.
    {
      vbs_old = (extData.currStoVectorRawPtr)[li_store_vbs];
      vbd_old = (extData.currStoVectorRawPtr)[li_store_vbd];
      vgs_old = (extData.currStoVectorRawPtr)[li_store_vgs];
      vds_old = (extData.currStoVectorRawPtr)[li_store_vds];
      von_local = (extData.currStoVectorRawPtr)[li_store_von];
    }
    else
    {  // no history
      vbs_old = vbs;
      vbd_old = vbd;
      vgs_old = vgs;
      vds_old = vds;
      von_local = 0.0;
    }
  }
  else
  {
    vbs_old = (extData.nextStoVectorRawPtr)[li_store_vbs];
    vbd_old = (extData.nextStoVectorRawPtr)[li_store_vbd];
    vgs_old = (extData.nextStoVectorRawPtr)[li_store_vgs];
    vds_old = (extData.nextStoVectorRawPtr)[li_store_vds];
    von_local = (extData.nextStoVectorRawPtr)[li_store_von];
  }

  vgdo = vgs_old - vds_old;

  // This next block performs checks on the junction voltages and
  // imposes limits on them if they are too big.
  // Note:  In the level=1 von is multiplied by dtype.  Here it is not.  They
  // are both right.

  if (getDeviceOptions().voltageLimiterFlag && !(getSolverState().initFixFlag && OFF))
  {

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << "  von_local = " << von_local << std::endl;
      Xyce::dout() << "  CONSTvt0  = " << CONSTvt0 << std::endl;
      Xyce::dout() << "  vcrit     = " << model_.vcrit << std::endl;
      Xyce::dout().width(3);
      Xyce::dout() << getSolverState().newtonIter;
      Xyce::dout().width(5);Xyce::dout() << getName();
      Xyce::dout() << " old :";
      Xyce::dout()<<" vgs:";Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vgs_old;
      Xyce::dout()<<" vds:";Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vds_old;
      Xyce::dout()<<" vbs:";Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vbs_old;
      Xyce::dout()<<" vbd:";Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vbd_old << std::endl;
      Xyce::dout().width(3);
      Xyce::dout() << getSolverState().newtonIter;
      Xyce::dout().width(5);Xyce::dout() << getName();
      Xyce::dout() << " Blim:";
      Xyce::dout()<<" vgs:";Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vgs;
      Xyce::dout()<<" vds:";Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vds;
      Xyce::dout()<<" vbs:";Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vbs;
      Xyce::dout()<<" vbd:";Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vbd << std::endl;
      Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
    }

    // only do this if we are beyond the first Newton iteration.  On the
    // first newton iteration, the "old" values are from a previous time
    // step.

    if (getSolverState().newtonIter >= 0)
    {
      if (vds_old >= 0.0)
      {
        vgs = devSupport.fetlim( vgs, vgs_old, von_local);
        vds = vgs - vgd;
        vds = devSupport.limvds( vds,  vds_old);
        vgd = vgs - vds;
      }
      else
      {
        vgd = devSupport.fetlim( vgd, vgdo, von_local);
        vds = vgs - vgd;
        vds = -devSupport.limvds( -vds, -vds_old );
        vgs = vgd + vds;
      }

      if (vds >= 0.0)
      {
        vbs = devSupport.pnjlim( vbs, vbs_old, CONSTvt0,
                                 model_.vcrit, &Check);
        vbd = vbs - vds;
      }
      else
      {
        vbd = devSupport.pnjlim( vbd, vbd_old, CONSTvt0,
                                 model_.vcrit, &Check);
        vbs = vbd + vds;
      }
    }

    // set the origFlag:
#ifdef Xyce_NEW_ORIG_TEST
    double vgs_diff = fabs(vgs - vgs_orig);
    double vds_diff = fabs(vds - vds_orig);
    double vbs_diff = fabs(vbs - vbs_orig);
    double vgd_diff = fabs(vgd - vgd_orig);

    bool noOrigFlag_vgs = 0;
    bool noOrigFlag_vds = 0;
    bool noOrigFlag_vbs = 0;
    bool noOrigFlag_vgd = 0;

    if (vgs_diff != 0.0)
    {
      if (vgs_orig != 0.0)
      {
        if ( fabs(vgs_diff/vgs_orig) > reltol) noOrigFlag_vgs = 1;
      }
      else
      {
        if ( fabs(vgs_diff) > voltTol ) noOrigFlag_vgs = 1;
      }
    }

    if (vds_diff != 0.0)
    {
      if (vds_orig != 0.0)
      {
        if ( fabs(vds_diff/vds_orig) > reltol) noOrigFlag_vds = 1;
      }
      else
      {
        if ( fabs(vds_diff) > voltTol ) noOrigFlag_vds = 1;
      }
    }

    if (vbs_diff != 0.0)
    {
      if (vbs_orig != 0.0)
      {
        if ( fabs(vbs_diff/vbs_orig) > reltol) noOrigFlag_vbs = 1;
      }
      else
      {
        if ( fabs(vbs_diff) > voltTol ) noOrigFlag_vbs = 1;
      }
    }

    if (vgd_diff != 0.0)
    {
      if (vgd_orig != 0.0)
      {
        if ( fabs(vgd_diff/vgd_orig) > reltol) noOrigFlag_vgd = 1;
      }
      else
      {
        if ( fabs(vgd_diff) > voltTol ) noOrigFlag_vgd = 1;
      }
    }

    origFlag = !( noOrigFlag_vgs || noOrigFlag_vds ||
                  noOrigFlag_vbs || noOrigFlag_vgd);

#else
    if (vgs_orig != vgs || vds_orig != vds ||
        vbs_orig != vbs || vbd_orig != vbd || vgd_orig != vgd) origFlag = 0;
#endif

    // for convergence testing:
    if (Check == 1) limitedFlag=true;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout().width(3);
      Xyce::dout() << getSolverState().newtonIter;
      Xyce::dout().width(5);Xyce::dout() << getName();
      Xyce::dout() << " Alim:";
      Xyce::dout()<<" vgs:";Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vgs;
      Xyce::dout()<<" vds:";Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vds;
      Xyce::dout()<<" vbs:";Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vbs;
      Xyce::dout()<<" vbd:";Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vbd;
      if (origFlag) Xyce::dout() << " SAME";
      else          Xyce::dout() << " DIFF";
      Xyce::dout() << std::endl;
      Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
    }

  } // getDeviceOptions().voltageLimiterFlag

  // Finished with what would have been the series of CKTmode
  // IF statements...

  // file: b3ld.c   line: 345
  // determine DC current and derivatives
  vbd = vbs - vds;
  vgd = vgs - vds;
  vgb = vgs - vbs;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << " mod type = " << model_.modType << std::endl;
    Xyce::dout() << "dtype = " << model_.dtype << std::endl;
    Xyce::dout() << "  vbs = " << vbs << std::endl;
    Xyce::dout() << "  vds = " << vds << std::endl;
    Xyce::dout() << "  vgs = " << vgs << std::endl;

    Xyce::dout() << "  vbd = " << vbd << std::endl;
    Xyce::dout() << "  vgd = " << vgd << std::endl;
    Xyce::dout() << "  vgb = " << vgs << std::endl;
    Xyce::dout() << " qdef = " << qdef << std::endl;
  }

  // Source/drain junction diode DC model begins
  Nvtm = vtm * model_.jctEmissionCoeff;

  if ((sourceArea <= 0.0) && (sourcePerimeter <= 0.0))
  {
    SourceSatCurrent = 1.0e-14;
  }
  else
  {
    SourceSatCurrent = sourceArea
                       * jctTempSatCurDensity
                       + sourcePerimeter
                       * jctSidewallTempSatCurDensity;
  }

  if (SourceSatCurrent <= 0.0)
  {
    gbs = getDeviceOptions().gmin;
    cbs = gbs * vbs;
  }
  else
  {
    if (model_.ijth == 0.0)
    {
      evbs = exp(vbs / Nvtm);
      gbs = SourceSatCurrent * evbs / Nvtm + getDeviceOptions().gmin;
      cbs = SourceSatCurrent * (evbs - 1.0) + getDeviceOptions().gmin * vbs;
    }
    else
    {
      if (vbs < vjsm)
      {
        evbs = exp(vbs / Nvtm);
        gbs = SourceSatCurrent * evbs / Nvtm + getDeviceOptions().gmin;
        cbs = SourceSatCurrent * (evbs - 1.0) + getDeviceOptions().gmin * vbs;
      }
      else
      {
        T0 = IsEvjsm / Nvtm;
        gbs = T0 + getDeviceOptions().gmin;
        cbs = IsEvjsm - SourceSatCurrent
              + T0 * (vbs - vjsm)
              + getDeviceOptions().gmin * vbs;
      }
    }
  }

  if ((drainArea <= 0.0) && (drainPerimeter <= 0.0))
  {
    DrainSatCurrent = 1.0e-14;
  }
  else
  {
    DrainSatCurrent = drainArea
                      * jctTempSatCurDensity
                      + drainPerimeter
                      * jctSidewallTempSatCurDensity;
  }

  if (DrainSatCurrent <= 0.0)
  {
    gbd = getDeviceOptions().gmin;
    cbd = gbd * vbd;
  }
  else
  {
    if (model_.ijth == 0.0)
    {
      evbd = exp(vbd / Nvtm);
      gbd = DrainSatCurrent * evbd / Nvtm + getDeviceOptions().gmin;
      cbd = DrainSatCurrent * (evbd - 1.0) + getDeviceOptions().gmin * vbd;
    }
    else
    {
      if (vbd < vjdm)
      {
        evbd = exp(vbd / Nvtm);
        gbd = DrainSatCurrent * evbd / Nvtm + getDeviceOptions().gmin;
        cbd = DrainSatCurrent * (evbd - 1.0) + getDeviceOptions().gmin * vbd;
      }
      else
      {
        T0 = IsEvjdm / Nvtm;
        gbd = T0 + getDeviceOptions().gmin;
        cbd = IsEvjdm - DrainSatCurrent
              + T0 * (vbd - vjdm)
              + getDeviceOptions().gmin * vbd;
      }
    }
  }
  // End of diode DC model

  if (vds >= 0.0)
  {   // normal mode
    mode = 1;
    Vds = vds;
    Vgs = vgs;
    Vbs = vbs;
  }
  else
  {   // inverse mode
    mode = -1;
    Vds = -vds;
    Vgs = vgd;
    Vbs = vbd;
  }

  // mosfet continuation.
  // This idea is based, loosely, on a paper by Jaijeet
  // Rosychowdhury.  If the artificial parameter flag has been enabled,
  // modify Vds and Vgs.
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "HOMOTOPY INFO: gainscale   = " << getSolverState().gainScale_ << std::endl
                 << "HOMOTOPY INFO: before vds  = " << Vds << std::endl
                 << "HOMOTOPY INFO: before vgst = " << Vgs << std::endl;
  }
  if (getSolverState().artParameterFlag_)
  {

    double alpha = getSolverState().gainScale_;
    double vgstConst = getDeviceOptions().vgstConst;

    Vds = devSupport.contVds (Vds,getSolverState().nltermScale_, getDeviceOptions().vdsScaleMin);
    Vgs = devSupport.contVgst(Vgs, alpha, vgstConst);
  }
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "HOMOTOPY INFO: after vds   = " << Vds << std::endl;
    Xyce::dout() << "HOMOTOPY INFO: after vgst  = " << Vgs << std::endl;
  }
  // end of mosfet continuation block.

  T0 = Vbs - paramPtr->vbsc - 0.001;
  T1global = sqrt(T0 * T0 - 0.004 * paramPtr->vbsc);
  Vbseff = paramPtr->vbsc + 0.5 * (T0 + T1global);
  dVbseff_dVb = 0.5 * (1.0 + T0 / T1global);

  if (Vbseff < Vbs) Vbseff = Vbs;

  if (Vbseff > 0.0)
  {
    T0 = paramPtr->phi / (paramPtr->phi + Vbseff);
    Phis = paramPtr->phi * T0;
    dPhis_dVb = -T0 * T0;
    sqrtPhis = paramPtr->phis3 / (paramPtr->phi + 0.5 * Vbseff);
    dsqrtPhis_dVb = -0.5 * sqrtPhis * sqrtPhis / paramPtr->phis3;
  }
  else
  {
    Phis = paramPtr->phi - Vbseff;
    dPhis_dVb = -1.0;
    sqrtPhis = sqrt(Phis);
    dsqrtPhis_dVb = -0.5 / sqrtPhis;
  }

  Xdep = paramPtr->Xdep0 * sqrtPhis / paramPtr->sqrtPhi;
  dXdep_dVb = (paramPtr->Xdep0 / paramPtr->sqrtPhi) * dsqrtPhis_dVb;

  Leff = paramPtr->leff;
  Vtm = vtm;
  // Vth Calculation
  T3 = sqrt(Xdep);
  V0 = paramPtr->vbi - paramPtr->phi;

  T0 = paramPtr->dvt2 * Vbseff;
  if (T0 >= - 0.5)
  {
    T1 = 1.0 + T0;
    T2 = paramPtr->dvt2;
  }
  else // Added to avoid any discontinuity problems caused by dvt2
  {
    T4 = 1.0 / (3.0 + 8.0 * T0);
    T1 = (1.0 + 3.0 * T0) * T4;
    T2 = paramPtr->dvt2 * T4 * T4;
  }

  lt1 = model_.factor1 * T3 * T1;
  dlt1_dVb = model_.factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

  T0 = paramPtr->dvt2w * Vbseff;
  if (T0 >= - 0.5)
  {
    T1 = 1.0 + T0;
    T2 = paramPtr->dvt2w;
  }
  else // Added to avoid any discontinuity problems caused by dvt2w
  {
    T4 = 1.0 / (3.0 + 8.0 * T0);
    T1 = (1.0 + 3.0 * T0) * T4;
    T2 = paramPtr->dvt2w * T4 * T4;
  }

  ltw = model_.factor1 * T3 * T1;
  dltw_dVb = model_.factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

  T0 = -0.5 * paramPtr->dvt1 * Leff / lt1;
  if (T0 > -CONSTEXP_THRESHOLD)
  {
    T1 = exp(T0);
    Theta0 = T1 * (1.0 + 2.0 * T1);
    dT1_dVb = -T0 / lt1 * T1 * dlt1_dVb;
    dTheta0_dVb = (1.0 + 4.0 * T1) * dT1_dVb;
  }
  else
  {
    T1 = CONSTMIN_EXP;
    Theta0 = T1 * (1.0 + 2.0 * T1);
    dTheta0_dVb = 0.0;
  }

  thetavth = paramPtr->dvt0 * Theta0;
  Delt_vth = thetavth * V0;
  dDelt_vth_dVb = paramPtr->dvt0 * dTheta0_dVb * V0;

  T0 = -0.5 * paramPtr->dvt1w * paramPtr->weff * Leff / ltw;
  if (T0 > -CONSTEXP_THRESHOLD)
  {
    T1 = exp(T0);
    T2 = T1 * (1.0 + 2.0 * T1);
    dT1_dVb = -T0 / ltw * T1 * dltw_dVb;
    dT2_dVb = (1.0 + 4.0 * T1) * dT1_dVb;
  }
  else
  {
    T1 = CONSTMIN_EXP;
    T2 = T1 * (1.0 + 2.0 * T1);
    dT2_dVb = 0.0;
  }

  T0 = paramPtr->dvt0w * T2;
  T2 = T0 * V0;
  dT2_dVb = paramPtr->dvt0w * dT2_dVb * V0;

  TempRatio =  temp / model_.tnom - 1.0;
  T0 = sqrt(1.0 + paramPtr->nlx / Leff);
  T1 = paramPtr->k1ox * (T0 - 1.0) * paramPtr->sqrtPhi
       + (paramPtr->kt1 + paramPtr->kt1l / Leff
          +  paramPtr->kt2 * Vbseff) * TempRatio;

  tmp2 = model_.tox * paramPtr->phi / (paramPtr->weff + paramPtr->w0);

  T3 = paramPtr->eta0 + paramPtr->etab * Vbseff;
  if (T3 < 1.0e-4) // avoid  discontinuity problems caused by etab
  {
    T9 = 1.0 / (3.0 - 2.0e4 * T3);
    T3 = (2.0e-4 - T3) * T9;
    T4 = T9 * T9;
  }
  else
  {
    T4 = 1.0;
  }

  dDIBL_Sft_dVd = T3 * paramPtr->theta0vb0;
  DIBL_Sft = dDIBL_Sft_dVd * Vds;

  Vth = model_.dtype * paramPtr->vth0 - paramPtr->k1
        * paramPtr->sqrtPhi + paramPtr->k1ox * sqrtPhis
        - paramPtr->k2ox * Vbseff - Delt_vth - T2 + (paramPtr->k3
                                                     + paramPtr->k3b * Vbseff) * tmp2 + T1 - DIBL_Sft;

  von = Vth;

  dVth_dVb = paramPtr->k1ox * dsqrtPhis_dVb - paramPtr->k2ox
             - dDelt_vth_dVb - dT2_dVb + paramPtr->k3b * tmp2
             - paramPtr->etab * Vds * paramPtr->theta0vb0 * T4
             + paramPtr->kt2 * TempRatio;

  dVth_dVd = -dDIBL_Sft_dVd;

  // Calculate n
  tmp2 = paramPtr->nfactor * CONSTEPSSI / Xdep;
  tmp3 = paramPtr->cdsc + paramPtr->cdscb * Vbseff
         + paramPtr->cdscd * Vds;
  tmp4 = (tmp2 + tmp3 * Theta0 + paramPtr->cit) / model_.cox;

  if (tmp4 >= -0.5)
  {
    n = 1.0 + tmp4;
    dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
              + paramPtr->cdscb * Theta0) / model_.cox;
    dn_dVd = paramPtr->cdscd * Theta0 / model_.cox;
  }
  else // avoid  discontinuity problems caused by tmp4
  {
    T0 = 1.0 / (3.0 + 8.0 * tmp4);
    n = (1.0 + 3.0 * tmp4) * T0;
    T0 *= T0;
    dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
              + paramPtr->cdscb * Theta0) / model_.cox * T0;
    dn_dVd = paramPtr->cdscd * Theta0 / model_.cox * T0;
  }

  // Poly Gate Si Depletion Effect
  T0 = paramPtr->vfb + paramPtr->phi;

  // added to avoid the problem caused by ngate
  if ((paramPtr->ngate > 1.e18) && (paramPtr->ngate < 1.e25) && (Vgs > T0))
  {
    T1 = 1.0e6 * CONSTQ * CONSTEPSSI * paramPtr->ngate
         / (model_.cox * model_.cox);
    T4 = sqrt(1.0 + 2.0 * (Vgs - T0) / T1);

    T2 = T1 * (T4 - 1.0);
    T3 = 0.5 * T2 * T2 / T1; // T3 = Vpoly
    T7 = 1.12 - T3 - 0.05;
    T6 = sqrt(T7 * T7 + 0.224);
    T5 = 1.12 - 0.5 * (T7 + T6);
    Vgs_eff = Vgs - T5;
    dVgs_eff_dVg = 1.0 - (0.5 - 0.5 / T4) * (1.0 + T7 / T6);
  }
  else
  {
    Vgs_eff = Vgs;
    dVgs_eff_dVg = 1.0;
  }
  Vgst = Vgs_eff - Vth;

  // Effective Vgst (Vgsteff) Calculation
  T10 = 2.0 * n * Vtm;
  VgstNVt = Vgst / T10;
  ExpArg = (2.0 * paramPtr->voff - Vgst) / T10;

  // MCJ: Very small Vgst
  if (VgstNVt > CONSTEXP_THRESHOLD)
  {
    Vgsteff = Vgst;
    dVgsteff_dVg = dVgs_eff_dVg;
    dVgsteff_dVd = -dVth_dVd;
    dVgsteff_dVb = -dVth_dVb;
  }
  else if (ExpArg > CONSTEXP_THRESHOLD)
  {
    T0 = (Vgst - paramPtr->voff) / (n * Vtm);
    ExpVgst = exp(T0);
    Vgsteff = Vtm * paramPtr->cdep0 / model_.cox * ExpVgst;
    dVgsteff_dVg = Vgsteff / (n * Vtm);
    dVgsteff_dVd = -dVgsteff_dVg * (dVth_dVd + T0 * Vtm * dn_dVd);
    dVgsteff_dVb = -dVgsteff_dVg * (dVth_dVb + T0 * Vtm * dn_dVb);
    dVgsteff_dVg *= dVgs_eff_dVg;
  }
  else
  {
    ExpVgst = exp(VgstNVt);
    T1 = T10 * log(1.0 + ExpVgst);
    dT1_dVg = ExpVgst / (1.0 + ExpVgst);
    dT1_dVb = -dT1_dVg * (dVth_dVb + Vgst / n * dn_dVb) + T1 / n * dn_dVb;
    dT1_dVd = -dT1_dVg * (dVth_dVd + Vgst / n * dn_dVd) + T1 / n * dn_dVd;

    dT2_dVg = -model_.cox / (Vtm * paramPtr->cdep0) * exp(ExpArg);
    T2 = 1.0 - T10 * dT2_dVg;

    dT2_dVd = -dT2_dVg * (dVth_dVd - 2.0 * Vtm * ExpArg * dn_dVd)
              + (T2 - 1.0) / n * dn_dVd;

    dT2_dVb = -dT2_dVg * (dVth_dVb - 2.0 * Vtm * ExpArg * dn_dVb)
              + (T2 - 1.0) / n * dn_dVb;

    Vgsteff = T1 / T2;
    T3 = T2 * T2;

    dVgsteff_dVg = (T2 * dT1_dVg - T1 * dT2_dVg) / T3 * dVgs_eff_dVg;
    dVgsteff_dVd = (T2 * dT1_dVd - T1 * dT2_dVd) / T3;
    dVgsteff_dVb = (T2 * dT1_dVb - T1 * dT2_dVb) / T3;
  }

  // Calculate Effective Channel Geometry
  T9 = sqrtPhis - paramPtr->sqrtPhi;
  Weff = paramPtr->weff -2.0 *(paramPtr->dwg * Vgsteff + paramPtr->dwb * T9);
  dWeff_dVg = -2.0 * paramPtr->dwg;
  dWeff_dVb = -2.0 * paramPtr->dwb * dsqrtPhis_dVb;

  if (Weff < 2.0e-8) // to avoid the discontinuity problem due to Weff
  {
    T0 = 1.0 / (6.0e-8 - 2.0 * Weff);
    Weff = 2.0e-8 * (4.0e-8 - Weff) * T0;
    T0 *= T0 * 4.0e-16;
    dWeff_dVg *= T0;
    dWeff_dVb *= T0;
  }

  T0 = paramPtr->prwg * Vgsteff + paramPtr->prwb * T9;
  if (T0 >= -0.9)
  {
    Rds = paramPtr->rds0 * (1.0 + T0);
    dRds_dVg = paramPtr->rds0 * paramPtr->prwg;
    dRds_dVb = paramPtr->rds0 * paramPtr->prwb * dsqrtPhis_dVb;
  }
  else // to avoid the discontinuity problem due to prwg and prwb
  {
    T1 = 1.0 / (17.0 + 20.0 * T0);
    Rds = paramPtr->rds0 * (0.8 + T0) * T1;
    T1 *= T1;
    dRds_dVg = paramPtr->rds0 * paramPtr->prwg * T1;
    dRds_dVb = paramPtr->rds0 * paramPtr->prwb * dsqrtPhis_dVb * T1;
  }

  // Calculate Abulk
  T1 = 0.5 * paramPtr->k1ox / sqrtPhis;
  dT1_dVb = -T1 / sqrtPhis * dsqrtPhis_dVb;

  T9 = sqrt(paramPtr->xj * Xdep);
  tmp1 = Leff + 2.0 * T9;
  T5 = Leff / tmp1;
  tmp2 = paramPtr->a0 * T5;
  tmp3 = paramPtr->weff + paramPtr->b1;
  tmp4 = paramPtr->b0 / tmp3;
  T2 = tmp2 + tmp4;
  dT2_dVb = -T9 / tmp1 / Xdep * dXdep_dVb;
  T6 = T5 * T5;
  T7 = T5 * T6;

  Abulk0 = 1.0 + T1 * T2;
  dAbulk0_dVb = T1 * tmp2 * dT2_dVb + T2 * dT1_dVb;

  T8 = paramPtr->ags * paramPtr->a0 * T7;
  dAbulk_dVg = -T1 * T8;
  Abulk = Abulk0 + dAbulk_dVg * Vgsteff;
  dAbulk_dVb = dAbulk0_dVb - T8 * Vgsteff * (dT1_dVb + 3.0 * T1 * dT2_dVb);

  if (Abulk0 < 0.1) // added to avoid the problems caused by Abulk0
  {
    T9 = 1.0 / (3.0 - 20.0 * Abulk0);
    Abulk0 = (0.2 - Abulk0) * T9;
    dAbulk0_dVb *= T9 * T9;
  }

  if (Abulk < 0.1) // added to avoid the problems caused by Abulk
  {
    T9 = 1.0 / (3.0 - 20.0 * Abulk);
    Abulk = (0.2 - Abulk) * T9;
    T10 = T9 * T9;
    dAbulk_dVb *= T10;
    dAbulk_dVg *= T10;
  }

  T2 = paramPtr->keta * Vbseff;
  if (T2 >= -0.9)
  {
    T0 = 1.0 / (1.0 + T2);
    dT0_dVb = -paramPtr->keta * T0 * T0;
  }
  else // added to avoid the problems caused by Keta
  {
    T1 = 1.0 / (0.8 + T2);
    T0 = (17.0 + 20.0 * T2) * T1;
    dT0_dVb = -paramPtr->keta * T1 * T1;
  }

  dAbulk_dVg *= T0;
  dAbulk_dVb = dAbulk_dVb * T0 + Abulk * dT0_dVb;
  dAbulk0_dVb = dAbulk0_dVb * T0 + Abulk0 * dT0_dVb;
  Abulk *= T0;
  Abulk0 *= T0;

  // Mobility calculation
  if (model_.mobMod == 1)
  {
    T0 = Vgsteff + Vth + Vth;
    T2 = paramPtr->ua + paramPtr->uc * Vbseff;
    T3 = T0 / model_.tox;
    T5 = T3 * (T2 + paramPtr->ub * T3);
    dDenomi_dVg = (T2 + 2.0 * paramPtr->ub * T3) / model_.tox;
    dDenomi_dVd = dDenomi_dVg * 2.0 * dVth_dVd;
    dDenomi_dVb = dDenomi_dVg * 2.0 * dVth_dVb + paramPtr->uc * T3;
  }
  else if (model_.mobMod == 2)
  {
    T5 = Vgsteff / model_.tox * (paramPtr->ua
                                 + paramPtr->uc * Vbseff + paramPtr->ub * Vgsteff / model_.tox);

    dDenomi_dVg = (paramPtr->ua + paramPtr->uc * Vbseff
                   + 2.0 * paramPtr->ub * Vgsteff / model_.tox) / model_.tox;

    dDenomi_dVd = 0.0;
    dDenomi_dVb = Vgsteff * paramPtr->uc / model_.tox;
  }
  else
  {
    T0 = Vgsteff + Vth + Vth;
    T2 = 1.0 + paramPtr->uc * Vbseff;
    T3 = T0 / model_.tox;
    T4 = T3 * (paramPtr->ua + paramPtr->ub * T3);
    T5 = T4 * T2;

    dDenomi_dVg = (paramPtr->ua + 2.0 * paramPtr->ub * T3) * T2 /model_.tox;
    dDenomi_dVd = dDenomi_dVg * 2.0 * dVth_dVd;
    dDenomi_dVb = dDenomi_dVg * 2.0 * dVth_dVb + paramPtr->uc * T4;
  }

  if (T5 >= -0.8)
  {
    Denomi = 1.0 + T5;
  }
  else // Added to avoid the discontinuity problem caused by ua and ub
  {
    T9 = 1.0 / (7.0 + 10.0 * T5);
    Denomi = (0.6 + T5) * T9;
    T9 *= T9;
    dDenomi_dVg *= T9;
    dDenomi_dVd *= T9;
    dDenomi_dVb *= T9;
  }

  ueff = paramPtr->u0temp / Denomi;
  T9 = -ueff / Denomi;
  dueff_dVg = T9 * dDenomi_dVg;
  dueff_dVd = T9 * dDenomi_dVd;
  dueff_dVb = T9 * dDenomi_dVb;

  // Saturation Drain Voltage  Vdsat
  WVCox = Weff * paramPtr->vsattemp * model_.cox;
  WVCoxRds = WVCox * Rds;

  Esat = 2.0 * paramPtr->vsattemp / ueff;
  EsatL = Esat * Leff;
  T0 = -EsatL /ueff;
  dEsatL_dVg = T0 * dueff_dVg;
  dEsatL_dVd = T0 * dueff_dVd;
  dEsatL_dVb = T0 * dueff_dVb;

  // Sqrt()
  a1 = paramPtr->a1;
  if (a1 == 0.0)
  {
    Lambda = paramPtr->a2;
    dLambda_dVg = 0.0;
  }
  else if (a1 > 0.0) // Added to avoid the discontinuity problem
    // caused by a1 and a2 (Lambda)
  {
    T0 = 1.0 - paramPtr->a2;
    T1 = T0 - paramPtr->a1 * Vgsteff - 0.0001;
    T2 = sqrt(T1 * T1 + 0.0004 * T0);
    Lambda = paramPtr->a2 + T0 - 0.5 * (T1 + T2);
    dLambda_dVg = 0.5 * paramPtr->a1 * (1.0 + T1 / T2);
  }
  else
  {
    T1 = paramPtr->a2 + paramPtr->a1 * Vgsteff - 0.0001;
    T2 = sqrt(T1 * T1 + 0.0004 * paramPtr->a2);
    Lambda = 0.5 * (T1 + T2);
    dLambda_dVg = 0.5 * paramPtr->a1 * (1.0 + T1 / T2);
  }

  Vgst2Vtm = Vgsteff + 2.0 * Vtm;
  if (Rds > 0)
  {
    tmp2 = dRds_dVg / Rds + dWeff_dVg / Weff;
    tmp3 = dRds_dVb / Rds + dWeff_dVb / Weff;
  }
  else
  {
    tmp2 = dWeff_dVg / Weff;
    tmp3 = dWeff_dVb / Weff;
  }

  if ((Rds == 0.0) && (Lambda == 1.0))
  {
    T0 = 1.0 / (Abulk * EsatL + Vgst2Vtm);
    tmp1 = 0.0;
    T1 = T0 * T0;
    T2 = Vgst2Vtm * T0;
    T3 = EsatL * Vgst2Vtm;
    Vdsat = T3 * T0;

    dT0_dVg = -(Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 1.0) * T1;
    dT0_dVd = -(Abulk * dEsatL_dVd) * T1;
    dT0_dVb = -(Abulk * dEsatL_dVb + dAbulk_dVb * EsatL) * T1;

    dVdsat_dVg = T3 * dT0_dVg + T2 * dEsatL_dVg + EsatL * T0;
    dVdsat_dVd = T3 * dT0_dVd + T2 * dEsatL_dVd;
    dVdsat_dVb = T3 * dT0_dVb + T2 * dEsatL_dVb;
  }
  else
  {
    tmp1 = dLambda_dVg / (Lambda * Lambda);
    T9 = Abulk * WVCoxRds;
    T8 = Abulk * T9;
    T7 = Vgst2Vtm * T9;
    T6 = Vgst2Vtm * WVCoxRds;
    T0 = 2.0 * Abulk * (T9 - 1.0 + 1.0 / Lambda);
    dT0_dVg = 2.0 * (T8 * tmp2 - Abulk * tmp1
                     + (2.0 * T9 + 1.0 / Lambda - 1.0) * dAbulk_dVg);

    dT0_dVb = 2.0 * (T8 * (2.0 / Abulk * dAbulk_dVb + tmp3)
                     + (1.0 / Lambda - 1.0) * dAbulk_dVb);
    dT0_dVd = 0.0;
    T1 = Vgst2Vtm * (2.0 / Lambda - 1.0) + Abulk * EsatL + 3.0 * T7;

    dT1_dVg = (2.0 / Lambda - 1.0) - 2.0 * Vgst2Vtm * tmp1
              + Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 3.0 * (T9
                                                                 + T7 * tmp2 + T6 * dAbulk_dVg);

    dT1_dVb = Abulk * dEsatL_dVb + EsatL * dAbulk_dVb
              + 3.0 * (T6 * dAbulk_dVb + T7 * tmp3);

    dT1_dVd = Abulk * dEsatL_dVd;

    T2 = Vgst2Vtm * (EsatL + 2.0 * T6);
    dT2_dVg = EsatL + Vgst2Vtm * dEsatL_dVg
              + T6 * (4.0 + 2.0 * Vgst2Vtm * tmp2);

    dT2_dVb = Vgst2Vtm * (dEsatL_dVb + 2.0 * T6 * tmp3);
    dT2_dVd = Vgst2Vtm * dEsatL_dVd;

    T3 = sqrt(T1 * T1 - 2.0 * T0 * T2);
    Vdsat = (T1 - T3) / T0;

    dT3_dVg = (T1 * dT1_dVg - 2.0 * (T0 * dT2_dVg + T2 * dT0_dVg)) / T3;
    dT3_dVd = (T1 * dT1_dVd - 2.0 * (T0 * dT2_dVd + T2 * dT0_dVd)) / T3;
    dT3_dVb = (T1 * dT1_dVb - 2.0 * (T0 * dT2_dVb + T2 * dT0_dVb)) / T3;

    dVdsat_dVg = (dT1_dVg - (T1 * dT1_dVg - dT0_dVg * T2
                             - T0 * dT2_dVg) / T3 - Vdsat * dT0_dVg) / T0;

    dVdsat_dVb = (dT1_dVb - (T1 * dT1_dVb - dT0_dVb * T2
                             - T0 * dT2_dVb) / T3 - Vdsat * dT0_dVb) / T0;

    dVdsat_dVd = (dT1_dVd - (T1 * dT1_dVd - T0 * dT2_dVd) / T3) / T0;
  }
  vdsat = Vdsat;

  // Effective Vds (Vdseff) Calculation
  T1 = Vdsat - Vds - paramPtr->delta;
  dT1_dVg = dVdsat_dVg;
  dT1_dVd = dVdsat_dVd - 1.0;
  dT1_dVb = dVdsat_dVb;

  T2 = sqrt(T1 * T1 + 4.0 * paramPtr->delta * Vdsat);
  T0 = T1 / T2;
  T3 = 2.0 * paramPtr->delta / T2;
  dT2_dVg = T0 * dT1_dVg + T3 * dVdsat_dVg;
  dT2_dVd = T0 * dT1_dVd + T3 * dVdsat_dVd;
  dT2_dVb = T0 * dT1_dVb + T3 * dVdsat_dVb;

  Vdseff = Vdsat - 0.5 * (T1 + T2);
  dVdseff_dVg = dVdsat_dVg - 0.5 * (dT1_dVg + dT2_dVg);
  dVdseff_dVd = dVdsat_dVd - 0.5 * (dT1_dVd + dT2_dVd);
  dVdseff_dVb = dVdsat_dVb - 0.5 * (dT1_dVb + dT2_dVb);

  // Added to eliminate non-zero Vdseff at Vds=0.0
  if (Vds == 0.0)
  {
    Vdseff = 0.0;
    dVdseff_dVg = 0.0;
    dVdseff_dVb = 0.0;
  }

  // Calculate VAsat
  tmp4 = 1.0 - 0.5 * Abulk * Vdsat / Vgst2Vtm;
  T9 = WVCoxRds * Vgsteff;
  T8 = T9 / Vgst2Vtm;
  T0 = EsatL + Vdsat + 2.0 * T9 * tmp4;

  T7 = 2.0 * WVCoxRds * tmp4;
  dT0_dVg = dEsatL_dVg + dVdsat_dVg + T7 * (1.0 + tmp2 * Vgsteff)
            - T8 * (Abulk * dVdsat_dVg - Abulk * Vdsat / Vgst2Vtm
                    + Vdsat * dAbulk_dVg);

  dT0_dVb = dEsatL_dVb + dVdsat_dVb + T7 * tmp3 * Vgsteff
            - T8 * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
  dT0_dVd = dEsatL_dVd + dVdsat_dVd - T8 * Abulk * dVdsat_dVd;

  T9 = WVCoxRds * Abulk;
  T1 = 2.0 / Lambda - 1.0 + T9;
  dT1_dVg = -2.0 * tmp1 +  WVCoxRds * (Abulk * tmp2 + dAbulk_dVg);
  dT1_dVb = dAbulk_dVb * WVCoxRds + T9 * tmp3;

  Vasat = T0 / T1;
  dVasat_dVg = (dT0_dVg - Vasat * dT1_dVg) / T1;
  dVasat_dVb = (dT0_dVb - Vasat * dT1_dVb) / T1;
  dVasat_dVd = dT0_dVd / T1;

  if (Vdseff > Vds) Vdseff = Vds;

  diffVds = Vds - Vdseff;

  // Calculate VACLM
  if ((paramPtr->pclm > 0.0) && (diffVds > 1.0e-10))
  {
    T0 = 1.0 / (paramPtr->pclm * Abulk * paramPtr->litl);
    dT0_dVb = -T0 / Abulk * dAbulk_dVb;
    dT0_dVg = -T0 / Abulk * dAbulk_dVg;

    T2 = Vgsteff / EsatL;
    T1 = Leff * (Abulk + T2);
    dT1_dVg = Leff * ((1.0 - T2 * dEsatL_dVg) / EsatL + dAbulk_dVg);
    dT1_dVb = Leff * (dAbulk_dVb - T2 * dEsatL_dVb / EsatL);
    dT1_dVd = -T2 * dEsatL_dVd / Esat;

    T9 = T0 * T1;
    VACLM = T9 * diffVds;
    dVACLM_dVg = T0 * dT1_dVg * diffVds - T9 * dVdseff_dVg
                 + T1 * diffVds * dT0_dVg;

    dVACLM_dVb = (dT0_dVb * T1 + T0 * dT1_dVb) * diffVds
                 - T9 * dVdseff_dVb;

    dVACLM_dVd = T0 * dT1_dVd * diffVds + T9 * (1.0 - dVdseff_dVd);
  }
  else
  {
    VACLM = CONSTMAX_EXP;
    dVACLM_dVd = dVACLM_dVg = dVACLM_dVb = 0.0;
  }

  // Calculate VADIBL
  if (paramPtr->thetaRout > 0.0)
  {
    T8 = Abulk * Vdsat;
    T0 = Vgst2Vtm * T8;
    dT0_dVg = Vgst2Vtm * Abulk * dVdsat_dVg + T8
              + Vgst2Vtm * Vdsat * dAbulk_dVg;

    dT0_dVb = Vgst2Vtm * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
    dT0_dVd = Vgst2Vtm * Abulk * dVdsat_dVd;

    T1 = Vgst2Vtm + T8;
    dT1_dVg = 1.0 + Abulk * dVdsat_dVg + Vdsat * dAbulk_dVg;
    dT1_dVb = Abulk * dVdsat_dVb + dAbulk_dVb * Vdsat;
    dT1_dVd = Abulk * dVdsat_dVd;

    T9 = T1 * T1;
    T2 = paramPtr->thetaRout;

    VADIBL = (Vgst2Vtm - T0 / T1) / T2;
    dVADIBL_dVg = (1.0 - dT0_dVg / T1 + T0 * dT1_dVg / T9) / T2;
    dVADIBL_dVb = (-dT0_dVb / T1 + T0 * dT1_dVb / T9) / T2;
    dVADIBL_dVd = (-dT0_dVd / T1 + T0 * dT1_dVd / T9) / T2;

    T7 = paramPtr->pdiblb * Vbseff;
    if (T7 >= -0.9)
    {
      T3 = 1.0 / (1.0 + T7);
      VADIBL *= T3;
      dVADIBL_dVg *= T3;
      dVADIBL_dVb = (dVADIBL_dVb - VADIBL * paramPtr->pdiblb) * T3;
      dVADIBL_dVd *= T3;
    }
    else // Added to avoid the discontinuity problem caused by pdiblcb
    {
      T4 = 1.0 / (0.8 + T7);
      T3 = (17.0 + 20.0 * T7) * T4;
      dVADIBL_dVg *= T3;
      dVADIBL_dVb = dVADIBL_dVb * T3 - VADIBL * paramPtr->pdiblb * T4 * T4;

      dVADIBL_dVd *= T3;
      VADIBL *= T3;
    }
  }
  else
  {
    VADIBL = CONSTMAX_EXP;
    dVADIBL_dVd = dVADIBL_dVg = dVADIBL_dVb = 0.0;
  }

  // Calculate VA
  T8 = paramPtr->pvag / EsatL;
  T9 = T8 * Vgsteff;
  if (T9 > -0.9)
  {
    T0 = 1.0 + T9;
    dT0_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL);
    dT0_dVb = -T9 * dEsatL_dVb / EsatL;
    dT0_dVd = -T9 * dEsatL_dVd / EsatL;
  }
  else /* Added to avoid the discontinuity problems caused by pvag */
  {
    T1 = 1.0 / (17.0 + 20.0 * T9);
    T0 = (0.8 + T9) * T1;
    T1 *= T1;
    dT0_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL) * T1;

    T9 *= T1 / EsatL;
    dT0_dVb = -T9 * dEsatL_dVb;
    dT0_dVd = -T9 * dEsatL_dVd;
  }

  tmp1 = VACLM * VACLM;
  tmp2 = VADIBL * VADIBL;
  tmp3 = VACLM + VADIBL;

  T1 = VACLM * VADIBL / tmp3;
  tmp3 *= tmp3;
  dT1_dVg = (tmp1 * dVADIBL_dVg + tmp2 * dVACLM_dVg) / tmp3;
  dT1_dVd = (tmp1 * dVADIBL_dVd + tmp2 * dVACLM_dVd) / tmp3;
  dT1_dVb = (tmp1 * dVADIBL_dVb + tmp2 * dVACLM_dVb) / tmp3;

  Va = Vasat + T0 * T1;
  dVa_dVg = dVasat_dVg + T1 * dT0_dVg + T0 * dT1_dVg;
  dVa_dVd = dVasat_dVd + T1 * dT0_dVd + T0 * dT1_dVd;
  dVa_dVb = dVasat_dVb + T1 * dT0_dVb + T0 * dT1_dVb;

  // Calculate VASCBE
  if (paramPtr->pscbe2 > 0.0)
  {
    if (diffVds > paramPtr->pscbe1 * paramPtr->litl / CONSTEXP_THRESHOLD)
    {
      T0 =  paramPtr->pscbe1 * paramPtr->litl / diffVds;
      VASCBE = Leff * exp(T0) / paramPtr->pscbe2;
      T1 = T0 * VASCBE / diffVds;
      dVASCBE_dVg = T1 * dVdseff_dVg;
      dVASCBE_dVd = -T1 * (1.0 - dVdseff_dVd);
      dVASCBE_dVb = T1 * dVdseff_dVb;
    }
    else
    {
      VASCBE = CONSTMAX_EXP * Leff/paramPtr->pscbe2;
      dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
    }
  }
  else
  {
    VASCBE = CONSTMAX_EXP;
    dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
  }

  // Calculate Ids
  CoxWovL = model_.cox * Weff / Leff;
  beta = ueff * CoxWovL;
  dbeta_dVg = CoxWovL * dueff_dVg + beta * dWeff_dVg / Weff;
  dbeta_dVd = CoxWovL * dueff_dVd;
  dbeta_dVb = CoxWovL * dueff_dVb + beta * dWeff_dVb / Weff;

  T0 = 1.0 - 0.5 * Abulk * Vdseff / Vgst2Vtm;
  dT0_dVg = -0.5 * (Abulk * dVdseff_dVg
                    - Abulk * Vdseff / Vgst2Vtm + Vdseff * dAbulk_dVg) / Vgst2Vtm;
  dT0_dVd = -0.5 * Abulk * dVdseff_dVd / Vgst2Vtm;
  dT0_dVb = -0.5 * (Abulk * dVdseff_dVb + dAbulk_dVb * Vdseff) / Vgst2Vtm;

  fgche1 = Vgsteff * T0;
  dfgche1_dVg = Vgsteff * dT0_dVg + T0;
  dfgche1_dVd = Vgsteff * dT0_dVd;
  dfgche1_dVb = Vgsteff * dT0_dVb;

  T9 = Vdseff / EsatL;
  fgche2 = 1.0 + T9;
  dfgche2_dVg = (dVdseff_dVg - T9 * dEsatL_dVg) / EsatL;
  dfgche2_dVd = (dVdseff_dVd - T9 * dEsatL_dVd) / EsatL;
  dfgche2_dVb = (dVdseff_dVb - T9 * dEsatL_dVb) / EsatL;

  gche = beta * fgche1 / fgche2;
  dgche_dVg = (beta * dfgche1_dVg + fgche1 * dbeta_dVg
               - gche * dfgche2_dVg) / fgche2;

  dgche_dVd = (beta * dfgche1_dVd + fgche1 * dbeta_dVd
               - gche * dfgche2_dVd) / fgche2;

  dgche_dVb = (beta * dfgche1_dVb + fgche1 * dbeta_dVb
               - gche * dfgche2_dVb) / fgche2;

  T0 = 1.0 + gche * Rds;
  T9 = Vdseff / T0;
  Idl = gche * T9;

  dIdl_dVg = (gche * dVdseff_dVg + T9 * dgche_dVg) / T0
             - Idl * gche / T0 * dRds_dVg ;

  dIdl_dVd = (gche * dVdseff_dVd + T9 * dgche_dVd) / T0;
  dIdl_dVb = (gche * dVdseff_dVb + T9 * dgche_dVb
              - Idl * dRds_dVb * gche) / T0;

  T9 =  diffVds / Va;
  T0 =  1.0 + T9;
  Idsa = Idl * T0;
  dIdsa_dVg = T0 * dIdl_dVg - Idl * (dVdseff_dVg + T9 * dVa_dVg) / Va;
  dIdsa_dVd = T0 * dIdl_dVd + Idl * (1.0 - dVdseff_dVd
                                     - T9 * dVa_dVd) / Va;

  dIdsa_dVb = T0 * dIdl_dVb - Idl * (dVdseff_dVb + T9 * dVa_dVb) / Va;

  T9 = diffVds / VASCBE;
  T0 = 1.0 + T9;
  Ids = Idsa * T0;

  Gm  = T0 * dIdsa_dVg - Idsa * (dVdseff_dVg + T9 * dVASCBE_dVg) / VASCBE;
  Gds = T0 * dIdsa_dVd + Idsa * (1.0 - dVdseff_dVd
                                 - T9 * dVASCBE_dVd) / VASCBE;
  Gmb = T0 * dIdsa_dVb - Idsa * (dVdseff_dVb
                                 + T9 * dVASCBE_dVb) / VASCBE;

  Gds += Gm * dVgsteff_dVd;
  Gmb += Gm * dVgsteff_dVb;
  Gm *= dVgsteff_dVg;
  Gmb *= dVbseff_dVb;

  // Substrate current begins
  tmp = paramPtr->alpha0 + paramPtr->alpha1 * Leff;
  if ((tmp <= 0.0) || (paramPtr->beta0 <= 0.0))
  {
    Isub = Gbd = Gbb = Gbg = 0.0;
  }
  else
  {
    T2 = tmp / Leff;
    if (diffVds > paramPtr->beta0 / CONSTEXP_THRESHOLD)
    {
      T0 = -paramPtr->beta0 / diffVds;
      T1 = T2 * diffVds * exp(T0);
      T3 = T1 / diffVds * (T0 - 1.0);
      dT1_dVg = T3 * dVdseff_dVg;
      dT1_dVd = T3 * (dVdseff_dVd - 1.0);
      dT1_dVb = T3 * dVdseff_dVb;
    }
    else
    {
      T3 = T2 * CONSTMIN_EXP;
      T1 = T3 * diffVds;
      dT1_dVg = -T3 * dVdseff_dVg;
      dT1_dVd = T3 * (1.0 - dVdseff_dVd);
      dT1_dVb = -T3 * dVdseff_dVb;
    }
    Isub = T1 * Idsa;
    Gbg = T1 * dIdsa_dVg + Idsa * dT1_dVg;
    Gbd = T1 * dIdsa_dVd + Idsa * dT1_dVd;
    Gbb = T1 * dIdsa_dVb + Idsa * dT1_dVb;

    Gbd += Gbg * dVgsteff_dVd;
    Gbb += Gbg * dVgsteff_dVb;
    Gbg *= dVgsteff_dVg;
    Gbb *= dVbseff_dVb; // bug fixing
  }

  // copy over local drain (channel) current vars to instance vars:
  cdrain = Ids;
  gds = Gds;
  gm = Gm;
  gmbs = Gmb;

  // copy over local substrate current vars to instance vars:
  gbbs = Gbb;
  gbgs = Gbg;
  gbds = Gbd;

  csub = Isub;

  //  thermal noise Qinv calculated from all capMod
  //  * 0, 1, 2 & 3 stored in iterI->qinv 1/1998
  if ((model_.xpart < 0) || (!ChargeComputationNeeded))
  {
    qgate  = qdrn = qsrc = qbulk = 0.0;
    cggb = cgsb = cgdb = 0.0;
    cdgb = cdsb = cddb = 0.0;
    cbgb = cbsb = cbdb = 0.0;
    cqdb = cqsb = cqgb = cqbb = 0.0;

    gtau = 0.0;
    goto finished;
  }
  else if (model_.capMod == 0)
  {
    if (Vbseff < 0.0)
    {
      Vbseff = Vbs;
      dVbseff_dVb = 1.0;
    }
    else
    {
      Vbseff = paramPtr->phi - Phis;
      dVbseff_dVb = -dPhis_dVb;
    }

    Vfb = paramPtr->vfbcv;
    Vth = Vfb + paramPtr->phi + paramPtr->k1ox * sqrtPhis;
    Vgst = Vgs_eff - Vth;
    dVth_dVb = paramPtr->k1ox * dsqrtPhis_dVb;
    dVgst_dVb = -dVth_dVb;
    dVgst_dVg = dVgs_eff_dVg;

    CoxWL = model_.cox * paramPtr->weffCV * paramPtr->leffCV;
    Arg1 = Vgs_eff - Vbseff - Vfb;

    if (Arg1 <= 0.0)
    {
      qgate = CoxWL * Arg1;
      qbulk = -qgate;
      qdrn = 0.0;

      cggb = CoxWL * dVgs_eff_dVg;
      cgdb = 0.0;
      cgsb = CoxWL * (dVbseff_dVb - dVgs_eff_dVg);

      cdgb = 0.0;
      cddb = 0.0;
      cdsb = 0.0;

      cbgb = -CoxWL * dVgs_eff_dVg;
      cbdb = 0.0;
      cbsb = -cgsb;
      qinv = 0.0;
    }
    else if (Vgst <= 0.0)
    {
      T1 = 0.5 * paramPtr->k1ox;
      T2 = sqrt(T1 * T1 + Arg1);
      qgate = CoxWL * paramPtr->k1ox * (T2 - T1);
      qbulk = -qgate;
      qdrn = 0.0;

      T0 = CoxWL * T1 / T2;
      cggb = T0 * dVgs_eff_dVg;
      cgdb = 0.0;
      cgsb = T0 * (dVbseff_dVb - dVgs_eff_dVg);

      cdgb = 0.0;
      cddb = 0.0;
      cdsb = 0.0;

      cbgb = -cggb;
      cbdb = 0.0;
      cbsb = -cgsb;
      qinv = 0.0;
    }
    else
    {
      One_Third_CoxWL = CoxWL / 3.0;
      Two_Third_CoxWL = 2.0 * One_Third_CoxWL;

      AbulkCV = Abulk0 * paramPtr->abulkCVfactor;
      dAbulkCV_dVb = paramPtr->abulkCVfactor * dAbulk0_dVb;

      Vdsat = Vgst / AbulkCV;
      dVdsat_dVg = dVgs_eff_dVg / AbulkCV;
      dVdsat_dVb = - (Vdsat * dAbulkCV_dVb + dVth_dVb)/ AbulkCV;

      if (model_.xpart > 0.5)
      { // 0/100 Charge partition model
        if (Vdsat <= Vds)
        {   // saturation region
          T1 = Vdsat / 3.0;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi - T1);

          T2 = -Two_Third_CoxWL * Vgst;
          qbulk = -(qgate + T2);
          qdrn = 0.0;

          cggb = One_Third_CoxWL * (3.0 - dVdsat_dVg)* dVgs_eff_dVg;
          T2 = -One_Third_CoxWL * dVdsat_dVb;

          cgsb = -(cggb + T2);
          cgdb = 0.0;

          cdgb = 0.0;
          cddb = 0.0;
          cdsb = 0.0;

          cbgb = -(cggb - Two_Third_CoxWL * dVgs_eff_dVg);
          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);

          cbsb = -(cbgb + T3);
          cbdb = 0.0;
          qinv = -(qgate + qbulk);
        }
        else
        {   // linear region
          Alphaz = Vgst / Vdsat;
          T1 = 2.0 * Vdsat - Vds;
          T2 = Vds / (3.0 * T1);
          T3 = T2 * Vds;
          T9 = 0.25 * CoxWL;
          T4 = T9 * Alphaz;
          T7 = 2.0 * Vds - T1 - 3.0 * T3;
          T8 = T3 - T1 - 2.0 * Vds;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi - 0.5 * (Vds - T3));

          T10 = T4 * T8;
          qdrn = T4 * T7;
          qbulk = -(qgate + qdrn + T10);

          T5 = T3 / T1;
          cggb = CoxWL * (1.0 - T5 * dVdsat_dVg) * dVgs_eff_dVg;

          T11 = -CoxWL * T5 * dVdsat_dVb;
          cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
          cgsb = -(cggb + T11 + cgdb);

          T6 = 1.0 / Vdsat;
          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

          T7 = T9 * T7;
          T8 = T9 * T8;
          T9 = 2.0 * T4 * (1.0 - 3.0 * T5);
          cdgb = (T7 * dAlphaz_dVg - T9* dVdsat_dVg) * dVgs_eff_dVg;

          T12 = T7 * dAlphaz_dVb - T9 * dVdsat_dVb;
          cddb = T4 * (3.0 - 6.0 * T2 - 3.0 * T5);
          cdsb = -(cdgb + T12 + cddb);

          T9  = 2.0 * T4 * (1.0 + T5);
          T10 = (T8 * dAlphaz_dVg - T9 * dVdsat_dVg) * dVgs_eff_dVg;
          T11 = T8 * dAlphaz_dVb - T9 * dVdsat_dVb;
          T12 = T4 * (2.0 * T2 + T5 - 1.0);
          T0  = -(T10 + T11 + T12);


          cbgb = -(cggb + cdgb + T10);
          cbdb = -(cgdb + cddb + T12);
          cbsb = -(cgsb + cdsb + T0);
          qinv = -(qgate + qbulk);
        }
      }
      else if (model_.xpart < 0.5)
      {   // 40/60 Charge partition model
        if (Vds >= Vdsat)
        { // saturation region
          T1 = Vdsat / 3.0;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi - T1);

          T2 = -Two_Third_CoxWL * Vgst;
          qbulk = -(qgate + T2);
          qdrn = 0.4 * T2;

          cggb = One_Third_CoxWL* (3.0 - dVdsat_dVg) * dVgs_eff_dVg;

          T2 = -One_Third_CoxWL * dVdsat_dVb;
          cgsb = -(cggb + T2);
          cgdb = 0.0;

          T3 = 0.4 * Two_Third_CoxWL;
          cdgb = -T3 * dVgs_eff_dVg;
          cddb = 0.0;

          T4 = T3 * dVth_dVb;
          cdsb = -(T4 + cdgb);

          cbgb = -(cggb - Two_Third_CoxWL * dVgs_eff_dVg);
          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
          cbsb = -(cbgb + T3);
          cbdb = 0.0;
          qinv = -(qgate + qbulk);
        }
        else
        {   // linear region
          Alphaz = Vgst / Vdsat;
          T1 = 2.0 * Vdsat - Vds;
          T2 = Vds / (3.0 * T1);
          T3 = T2 * Vds;
          T9 = 0.25 * CoxWL;
          T4 = T9 * Alphaz;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi - 0.5 * (Vds - T3));

          T5 = T3 / T1;
          cggb = CoxWL * (1.0 - T5 * dVdsat_dVg) * dVgs_eff_dVg;
          tmp = -CoxWL * T5 * dVdsat_dVb;
          cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
          cgsb = -(cggb + cgdb + tmp);

          T6 = 1.0 / Vdsat;
          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

          T6 = 8.0 * Vdsat * Vdsat - 6.0 * Vdsat * Vds + 1.2 * Vds * Vds;
          T8 = T2 / T1;
          T7 = Vds - T1 - T8 * T6;
          qdrn = T4 * T7;
          T7 *= T9;
          tmp = T8 / T1;
          tmp1 = T4*(2.0 - 4.0 * tmp * T6 + T8 *(16.0 * Vdsat - 6.0 *Vds));

          cdgb = (T7 *dAlphaz_dVg - tmp1 *dVdsat_dVg) *dVgs_eff_dVg;

          T10 = T7 * dAlphaz_dVb - tmp1 * dVdsat_dVb;
          cddb = T4 * (2.0 - (1.0 / (3.0 * T1
                                     * T1) + 2.0 * tmp) * T6 + T8
                       * (6.0 * Vdsat - 2.4 * Vds));

          cdsb = -(cdgb + T10 + cddb);

          T7 = 2.0 * (T1 + T3);
          qbulk = -(qgate - T4 * T7);
          T7 *= T9;
          T0 = 4.0 * T4 * (1.0 - T5);
          T12 = (-T7 * dAlphaz_dVg - cdgb
                 - T0 * dVdsat_dVg) * dVgs_eff_dVg;
          T11 = -T7 * dAlphaz_dVb - T10 - T0 * dVdsat_dVb;
          T10 = -4.0 * T4 * (T2 - 0.5 + 0.5 * T5)
                - cddb;

          tmp = -(T10 + T11 + T12);

          cbgb = -(cggb + cdgb + T12);
          cbdb = -(cgdb + cddb + T11);
          cbsb = -(cgsb + cdsb + tmp);
          qinv = -(qgate + qbulk);
        }
      }
      else
      {   // 50/50 partitioning
        if (Vds >= Vdsat)
        {   // saturation region
          T1 = Vdsat / 3.0;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi - T1);
          T2 = -Two_Third_CoxWL * Vgst;
          qbulk = -(qgate + T2);
          qdrn = 0.5 * T2;

          cggb = One_Third_CoxWL * (3.0 -dVdsat_dVg) * dVgs_eff_dVg;
          T2 = -One_Third_CoxWL * dVdsat_dVb;
          cgsb = -(cggb + T2);
          cgdb = 0.0;

          cdgb = -One_Third_CoxWL * dVgs_eff_dVg;
          cddb = 0.0;
          T4 = One_Third_CoxWL * dVth_dVb;
          cdsb = -(T4 + cdgb);

          cbgb = -(cggb - Two_Third_CoxWL * dVgs_eff_dVg);
          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
          cbsb = -(cbgb + T3);
          cbdb = 0.0;
          qinv = -(qgate + qbulk);
        }
        else
        {   // linear region
          Alphaz = Vgst / Vdsat;
          T1 = 2.0 * Vdsat - Vds;
          T2 = Vds / (3.0 * T1);
          T3 = T2 * Vds;
          T9 = 0.25 * CoxWL;
          T4 = T9 * Alphaz;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi - 0.5 * (Vds - T3))
                  ;

          T5 = T3 / T1;
          cggb = CoxWL * (1.0 - T5 * dVdsat_dVg) * dVgs_eff_dVg;

          tmp = -CoxWL * T5 * dVdsat_dVb;
          cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
          cgsb = -(cggb + cgdb + tmp);

          T6 = 1.0 / Vdsat;
          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

          T7 = T1 + T3;
          qdrn = -T4 * T7;
          qbulk = - (qgate + qdrn + qdrn);
          T7 *= T9;
          T0 = T4 * (2.0 * T5 - 2.0);

          cdgb = (T0 * dVdsat_dVg - T7 *dAlphaz_dVg) *dVgs_eff_dVg;
          T12 = T0 * dVdsat_dVb - T7 * dAlphaz_dVb;
          cddb = T4 * (1.0 - 2.0 * T2 - T5);
          cdsb = -(cdgb + T12 + cddb);

          cbgb = -(cggb + 2.0 * cdgb);
          cbdb = -(cgdb + 2.0 * cddb);
          cbsb = -(cgsb + 2.0 * cdsb);
          qinv = -(qgate + qbulk);
        }
      }
    }
  }
  else
  {
    if (Vbseff < 0.0)
    {
      VbseffCV = Vbseff;
      dVbseffCV_dVb = 1.0;
    }
    else
    {
      VbseffCV = paramPtr->phi - Phis;
      dVbseffCV_dVb = -dPhis_dVb;
    }

    CoxWL = model_.cox * paramPtr->weffCV * paramPtr->leffCV;

    // Seperate VgsteffCV with noff and voffcv
    noff = n * paramPtr->noff;
    dnoff_dVd = paramPtr->noff * dn_dVd;
    dnoff_dVb = paramPtr->noff * dn_dVb;
    T0 = Vtm * noff;
    voffcv = paramPtr->voffcv;
    VgstNVt = (Vgst - voffcv) / T0;

    if (VgstNVt > CONSTEXP_THRESHOLD)
    {
      Vgsteff = Vgst - voffcv;
      dVgsteff_dVg = dVgs_eff_dVg;
      dVgsteff_dVd = -dVth_dVd;
      dVgsteff_dVb = -dVth_dVb;
    }
    else if (VgstNVt < -CONSTEXP_THRESHOLD)
    {
      Vgsteff = T0 * log(1.0 + CONSTMIN_EXP);
      dVgsteff_dVg = 0.0;
      dVgsteff_dVd = Vgsteff / noff;
      dVgsteff_dVb = dVgsteff_dVd * dnoff_dVb;
      dVgsteff_dVd *= dnoff_dVd;
    }
    else
    {
      ExpVgst = exp(VgstNVt);
      Vgsteff = T0 * log(1.0 + ExpVgst);
      dVgsteff_dVg = ExpVgst / (1.0 + ExpVgst);
      dVgsteff_dVd = -dVgsteff_dVg * (dVth_dVd + (Vgst - voffcv)
                                      / noff * dnoff_dVd) + Vgsteff / noff * dnoff_dVd;
      dVgsteff_dVb = -dVgsteff_dVg * (dVth_dVb + (Vgst - voffcv)
                                      / noff * dnoff_dVb) + Vgsteff / noff * dnoff_dVb;
      dVgsteff_dVg *= dVgs_eff_dVg;
    } // End of VgsteffCV

    if (model_.capMod == 1)
    {
      Vfb = paramPtr->vfbzb;
      Arg1 = Vgs_eff - VbseffCV - Vfb - Vgsteff;

      if (Arg1 <= 0.0)
      {
        qgate = CoxWL * Arg1;
        Cgg = CoxWL * (dVgs_eff_dVg - dVgsteff_dVg);
        Cgd = -CoxWL * dVgsteff_dVd;
        Cgb = -CoxWL * (dVbseffCV_dVb + dVgsteff_dVb);
      }
      else
      {
        T0 = 0.5 * paramPtr->k1ox;
        T1 = sqrt(T0 * T0 + Arg1);
        T2 = CoxWL * T0 / T1;

        qgate = CoxWL * paramPtr->k1ox * (T1 - T0);

        Cgg = T2 * (dVgs_eff_dVg - dVgsteff_dVg);
        Cgd = -T2 * dVgsteff_dVd;
        Cgb = -T2 * (dVbseffCV_dVb + dVgsteff_dVb);
      }
      qbulk = -qgate;
      Cbg = -Cgg;
      Cbd = -Cgd;
      Cbb = -Cgb;

      One_Third_CoxWL = CoxWL / 3.0;
      Two_Third_CoxWL = 2.0 * One_Third_CoxWL;
      AbulkCV = Abulk0 * paramPtr->abulkCVfactor;
      dAbulkCV_dVb = paramPtr->abulkCVfactor * dAbulk0_dVb;
      VdsatCV = Vgsteff / AbulkCV;

      if (VdsatCV < Vds)
      {
        dVdsatCV_dVg = 1.0 / AbulkCV;
        dVdsatCV_dVb = -VdsatCV * dAbulkCV_dVb / AbulkCV;
        T0 = Vgsteff - VdsatCV / 3.0;
        dT0_dVg = 1.0 - dVdsatCV_dVg / 3.0;
        dT0_dVb = -dVdsatCV_dVb / 3.0;
        qgate += CoxWL * T0;
        Cgg1 = CoxWL * dT0_dVg;
        Cgb1 = CoxWL * dT0_dVb + Cgg1 * dVgsteff_dVb;
        Cgd1 = Cgg1 * dVgsteff_dVd;
        Cgg1 *= dVgsteff_dVg;
        Cgg += Cgg1;
        Cgb += Cgb1;
        Cgd += Cgd1;

        T0 = VdsatCV - Vgsteff;
        dT0_dVg = dVdsatCV_dVg - 1.0;
        dT0_dVb = dVdsatCV_dVb;
        qbulk += One_Third_CoxWL * T0;
        Cbg1 = One_Third_CoxWL * dT0_dVg;
        Cbb1 = One_Third_CoxWL * dT0_dVb + Cbg1 * dVgsteff_dVb;
        Cbd1 = Cbg1 * dVgsteff_dVd;
        Cbg1 *= dVgsteff_dVg;
        Cbg += Cbg1;
        Cbb += Cbb1;
        Cbd += Cbd1;

        if (model_.xpart > 0.5)      T0 = -Two_Third_CoxWL;
        else if (model_.xpart < 0.5) T0 = -0.4 * CoxWL;
        else                         T0 = -One_Third_CoxWL;

        qsrc = T0 * Vgsteff;
        Csg = T0 * dVgsteff_dVg;
        Csb = T0 * dVgsteff_dVb;
        Csd = T0 * dVgsteff_dVd;
        Cgb *= dVbseff_dVb;
        Cbb *= dVbseff_dVb;
        Csb *= dVbseff_dVb;
      }
      else
      {
        T0 = AbulkCV * Vds;
        T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.e-20);
        T2 = Vds / T1;
        T3 = T0 * T2;
        dT3_dVg = -12.0 * T2 * T2 * AbulkCV;
        dT3_dVd = 6.0 * T0 * (4.0 * Vgsteff - T0) / T1 / T1 - 0.5;
        dT3_dVb = 12.0 * T2 * T2 * dAbulkCV_dVb * Vgsteff;

        qgate += CoxWL * (Vgsteff - 0.5 * Vds + T3);
        Cgg1 = CoxWL * (1.0 + dT3_dVg);
        Cgb1 = CoxWL * dT3_dVb + Cgg1 * dVgsteff_dVb;
        Cgd1 = CoxWL * dT3_dVd + Cgg1 * dVgsteff_dVd;
        Cgg1 *= dVgsteff_dVg;
        Cgg += Cgg1;
        Cgb += Cgb1;
        Cgd += Cgd1;

        qbulk += CoxWL * (1.0 - AbulkCV) * (0.5 * Vds - T3);
        Cbg1 = -CoxWL * ((1.0 - AbulkCV) * dT3_dVg);
        Cbb1 = -CoxWL * ((1.0 - AbulkCV) * dT3_dVb
                         + (0.5 * Vds - T3) * dAbulkCV_dVb)
               + Cbg1 * dVgsteff_dVb;
        Cbd1 = -CoxWL * (1.0 - AbulkCV) * dT3_dVd
               + Cbg1 * dVgsteff_dVd;
        Cbg1 *= dVgsteff_dVg;
        Cbg += Cbg1;
        Cbb += Cbb1;
        Cbd += Cbd1;

        if (model_.xpart > 0.5)
        {   // 0/100 Charge petition model
          T1 = T1 + T1;
          qsrc = -CoxWL * (0.5 * Vgsteff + 0.25 * T0 - T0 * T0 / T1);
          Csg = -CoxWL * (0.5 + 24.0 * T0 * Vds / T1 / T1 * AbulkCV);
          Csb = -CoxWL * (0.25 * Vds * dAbulkCV_dVb
                          - 12.0 * T0 * Vds / T1 / T1 * (4.0 * Vgsteff - T0)
                          * dAbulkCV_dVb) + Csg * dVgsteff_dVb;
          Csd = -CoxWL * (0.25 * AbulkCV - 12.0 * AbulkCV * T0
                          / T1 / T1 * (4.0 * Vgsteff - T0))
                + Csg * dVgsteff_dVd;
          Csg *= dVgsteff_dVg;
        }
        else if (model_.xpart < 0.5)
        {   // 40/60 Charge petition model
          T1 = T1 / 12.0;
          T2 = 0.5 * CoxWL / (T1 * T1);
          T3 = Vgsteff * (2.0 * T0 * T0 / 3.0 + Vgsteff
                          * (Vgsteff - 4.0 * T0 / 3.0))
               - 2.0 * T0 * T0 * T0 / 15.0;
          qsrc = -T2 * T3;
          T4 = 4.0 / 3.0 * Vgsteff * (Vgsteff - T0)
               + 0.4 * T0 * T0;
          Csg = -2.0 * qsrc / T1 - T2 * (Vgsteff * (3.0
                                                    * Vgsteff - 8.0 * T0 / 3.0)
                                         + 2.0 * T0 * T0 / 3.0);
          Csb = (qsrc / T1 * Vds + T2 * T4 * Vds) * dAbulkCV_dVb
                + Csg * dVgsteff_dVb;
          Csd = (qsrc / T1 + T2 * T4) * AbulkCV
                + Csg * dVgsteff_dVd;
          Csg *= dVgsteff_dVg;
        }
        else
        {   // 50/50 Charge petition model
          qsrc = -0.5 * (qgate + qbulk);
          Csg = -0.5 * (Cgg1 + Cbg1);
          Csb = -0.5 * (Cgb1 + Cbb1);
          Csd = -0.5 * (Cgd1 + Cbd1);
        }
        Cgb *= dVbseff_dVb;
        Cbb *= dVbseff_dVb;
        Csb *= dVbseff_dVb;
      }
      qdrn = -(qgate + qbulk + qsrc);
      cggb = Cgg;
      cgsb = -(Cgg + Cgd + Cgb);
      cgdb = Cgd;
      cdgb = -(Cgg + Cbg + Csg);
      cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb + Csg + Csd + Csb);
      cddb = -(Cgd + Cbd + Csd);
      cbgb = Cbg;
      cbsb = -(Cbg + Cbd + Cbb);
      cbdb = Cbd;
      qinv = -(qgate + qbulk);
    }
    else if (model_.capMod == 2)
    {
      Vfb = paramPtr->vfbzb;
      V3 = Vfb - Vgs_eff + VbseffCV - CONSTDELTA_3;
      if (Vfb <= 0.0)
      {
        T0 = sqrt(V3 * V3 - 4.0 * CONSTDELTA_3 * Vfb);
        T2 = -CONSTDELTA_3 / T0;
      }
      else
      {
        T0 = sqrt(V3 * V3 + 4.0 * CONSTDELTA_3 * Vfb);
        T2 = CONSTDELTA_3 / T0;
      }

      T1 = 0.5 * (1.0 + V3 / T0);
      Vfbeff = Vfb - 0.5 * (V3 + T0);
      dVfbeff_dVg = T1 * dVgs_eff_dVg;
      dVfbeff_dVb = -T1 * dVbseffCV_dVb;
      Qac0 = CoxWL * (Vfbeff - Vfb);
      dQac0_dVg = CoxWL * dVfbeff_dVg;
      dQac0_dVb = CoxWL * dVfbeff_dVb;

      T0 = 0.5 * paramPtr->k1ox;
      T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
      if (paramPtr->k1ox == 0.0)
      {
        T1 = 0.0;
        T2 = 0.0;
      }
      else if (T3 < 0.0)
      {
        T1 = T0 + T3 / paramPtr->k1ox;
        T2 = CoxWL;
      }
      else
      {
        T1 = sqrt(T0 * T0 + T3);
        T2 = CoxWL * T0 / T1;
      }

      Qsub0 = CoxWL * paramPtr->k1ox * (T1 - T0);

      dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg);
      dQsub0_dVd = -T2 * dVgsteff_dVd;
      dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb + dVgsteff_dVb);

      AbulkCV = Abulk0 * paramPtr->abulkCVfactor;
      dAbulkCV_dVb = paramPtr->abulkCVfactor * dAbulk0_dVb;
      VdsatCV = Vgsteff / AbulkCV;

      V4 = VdsatCV - Vds - CONSTDELTA_4;
      T0 = sqrt(V4 * V4 + 4.0 * CONSTDELTA_4 * VdsatCV);
      VdseffCV = VdsatCV - 0.5 * (V4 + T0);
      T1 = 0.5 * (1.0 + V4 / T0);
      T2 = CONSTDELTA_4 / T0;
      T3 = (1.0 - T1 - T2) / AbulkCV;
      dVdseffCV_dVg = T3;
      dVdseffCV_dVd = T1;
      dVdseffCV_dVb = -T3 * VdsatCV * dAbulkCV_dVb;
      // Added to eliminate non-zero VdseffCV at Vds=0.0
      if (Vds == 0.0)
      {
        VdseffCV = 0.0;
        dVdseffCV_dVg = 0.0;
        dVdseffCV_dVb = 0.0;
      }

      T0 = AbulkCV * VdseffCV;
      T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1e-20);
      T2 = VdseffCV / T1;
      T3 = T0 * T2;

      T4 = (1.0 - 12.0 * T2 * T2 * AbulkCV);
      T5 = (6.0 * T0 * (4.0 * Vgsteff - T0) / (T1 * T1) - 0.5);
      T6 = 12.0 * T2 * T2 * Vgsteff;

      qinoi = -CoxWL * (Vgsteff - 0.5 * T0 + AbulkCV * T3);
      qgate = CoxWL * (Vgsteff - 0.5 * VdseffCV + T3);
      Cgg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
      Cgd1 = CoxWL * T5 * dVdseffCV_dVd + Cgg1 * dVgsteff_dVd;
      Cgb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
             + Cgg1 * dVgsteff_dVb;
      Cgg1 *= dVgsteff_dVg;

      T7 = 1.0 - AbulkCV;
      qbulk = CoxWL * T7 * (0.5 * VdseffCV - T3);
      T4 = -T7 * (T4 - 1.0);
      T5 = -T7 * T5;
      T6 = -(T7 * T6 + (0.5 * VdseffCV - T3));
      Cbg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
      Cbd1 = CoxWL * T5 * dVdseffCV_dVd + Cbg1 * dVgsteff_dVd;
      Cbb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
             + Cbg1 * dVgsteff_dVb;
      Cbg1 *= dVgsteff_dVg;

      if (model_.xpart > 0.5)
      {   // 0/100 Charge petition model
        T1 = T1 + T1;
        qsrc = -CoxWL * (0.5 * Vgsteff + 0.25 * T0 - T0 * T0 / T1);
        T7 = (4.0 * Vgsteff - T0) / (T1 * T1);
        T4 = -(0.5 + 24.0 * T0 * T0 / (T1 * T1));
        T5 = -(0.25 * AbulkCV - 12.0 * AbulkCV * T0 * T7);
        T6 = -(0.25 * VdseffCV - 12.0 * T0 * VdseffCV * T7);
        Csg = CoxWL * (T4 + T5 * dVdseffCV_dVg);
        Csd = CoxWL * T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
        Csb = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
              + Csg * dVgsteff_dVb;
        Csg *= dVgsteff_dVg;
      }
      else if (model_.xpart < 0.5)
      {   // 40/60 Charge petition model
        T1 = T1 / 12.0;
        T2 = 0.5 * CoxWL / (T1 * T1);
        T3 = Vgsteff *(2.0 * T0 *T0/3.0 +Vgsteff *(Vgsteff - 4.0 *T0/ 3.0))
             - 2.0 * T0 * T0 * T0 / 15.0;
        qsrc = -T2 * T3;
        T7 = 4.0 / 3.0 * Vgsteff * (Vgsteff - T0)
             + 0.4 * T0 * T0;
        T4 = -2.0 * qsrc / T1 - T2 * (Vgsteff * (3.0
                                                 * Vgsteff - 8.0 * T0 / 3.0)
                                      + 2.0 * T0 * T0 / 3.0);
        T5 = (qsrc / T1 + T2 * T7) * AbulkCV;
        T6 = (qsrc / T1 * VdseffCV + T2 * T7 * VdseffCV);
        Csg = (T4 + T5 * dVdseffCV_dVg);
        Csd = T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
        Csb = (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
              + Csg * dVgsteff_dVb;
        Csg *= dVgsteff_dVg;
      }
      else
      {   // 50/50 Charge petition model
        qsrc = -0.5 * (qgate + qbulk);
        Csg = -0.5 * (Cgg1 + Cbg1);
        Csb = -0.5 * (Cgb1 + Cbb1);
        Csd = -0.5 * (Cgd1 + Cbd1);
      }

      qgate += Qac0 + Qsub0;
      qbulk -= (Qac0 + Qsub0);
      qdrn = -(qgate + qbulk + qsrc);

      Cgg = dQac0_dVg + dQsub0_dVg + Cgg1;
      Cgd = dQsub0_dVd + Cgd1;
      Cgb = dQac0_dVb + dQsub0_dVb + Cgb1;

      Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
      Cbd = Cbd1 - dQsub0_dVd;
      Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;

      Cgb *= dVbseff_dVb;
      Cbb *= dVbseff_dVb;
      Csb *= dVbseff_dVb;

      cggb = Cgg;
      cgsb = -(Cgg + Cgd + Cgb);
      cgdb = Cgd;
      cdgb = -(Cgg + Cbg + Csg);
      cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb + Csg + Csd + Csb);
      cddb = -(Cgd + Cbd + Csd);
      cbgb = Cbg;
      cbsb = -(Cbg + Cbd + Cbb);
      cbdb = Cbd;
      qinv = qinoi;
    }

    // New Charge-Thickness capMod (CTM) begins
    else if (model_.capMod == 3)
    {
      V3 = paramPtr->vfbzb - Vgs_eff + VbseffCV - CONSTDELTA_3;
      if (paramPtr->vfbzb <= 0.0)
      {
        T0 = sqrt(V3 * V3 - 4.0 * CONSTDELTA_3 * paramPtr->vfbzb);
        T2 = -CONSTDELTA_3 / T0;
      }
      else
      {
        T0 = sqrt(V3 * V3 + 4.0 * CONSTDELTA_3 * paramPtr->vfbzb);
        T2 = CONSTDELTA_3 / T0;
      }

      T1 = 0.5 * (1.0 + V3 / T0);
      Vfbeff = paramPtr->vfbzb - 0.5 * (V3 + T0);
      dVfbeff_dVg = T1 * dVgs_eff_dVg;
      dVfbeff_dVb = -T1 * dVbseffCV_dVb;

      Cox = model_.cox;
      Tox = 1.0e8 * model_.tox;
      T0 = (Vgs_eff - VbseffCV - paramPtr->vfbzb) / Tox;
      dT0_dVg = dVgs_eff_dVg / Tox;
      dT0_dVb = -dVbseffCV_dVb / Tox;

      tmp = T0 * paramPtr->acde;
      if ((-CONSTEXP_THRESHOLD < tmp) && (tmp < CONSTEXP_THRESHOLD))
      {
        Tcen = paramPtr->ldeb * exp(tmp);
        dTcen_dVg = paramPtr->acde * Tcen;
        dTcen_dVb = dTcen_dVg * dT0_dVb;
        dTcen_dVg *= dT0_dVg;
      }
      else if (tmp <= -CONSTEXP_THRESHOLD)
      {
        Tcen = paramPtr->ldeb * CONSTMIN_EXP;
        dTcen_dVg = dTcen_dVb = 0.0;
      }
      else
      {
        Tcen = paramPtr->ldeb * CONSTMAX_EXP;
        dTcen_dVg = dTcen_dVb = 0.0;
      }

      LINK = 1.0e-3 * model_.tox;
      V3 = paramPtr->ldeb - Tcen - LINK;
      V4 = sqrt(V3 * V3 + 4.0 * LINK * paramPtr->ldeb);
      Tcen = paramPtr->ldeb - 0.5 * (V3 + V4);
      T1 = 0.5 * (1.0 + V3 / V4);
      dTcen_dVg *= T1;
      dTcen_dVb *= T1;

      Ccen = CONSTEPSSI / Tcen;
      T2 = Cox / (Cox + Ccen);
      Coxeff = T2 * Ccen;
      T3 = -Ccen / Tcen;
      dCoxeff_dVg = T2 * T2 * T3;
      dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
      dCoxeff_dVg *= dTcen_dVg;
      CoxWLcen = CoxWL * Coxeff / Cox;

      Qac0 = CoxWLcen * (Vfbeff - paramPtr->vfbzb);
      QovCox = Qac0 / Coxeff;
      dQac0_dVg = CoxWLcen * dVfbeff_dVg + QovCox * dCoxeff_dVg;
      dQac0_dVb = CoxWLcen * dVfbeff_dVb + QovCox * dCoxeff_dVb;

      T0 = 0.5 * paramPtr->k1ox;
      T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
      if (paramPtr->k1ox == 0.0)
      {
        T1 = 0.0;
        T2 = 0.0;
      }
      else if (T3 < 0.0)
      {
        T1 = T0 + T3 / paramPtr->k1ox;
        T2 = CoxWLcen;
      }
      else
      {
        T1 = sqrt(T0 * T0 + T3);
        T2 = CoxWLcen * T0 / T1;
      }

      Qsub0 = CoxWLcen * paramPtr->k1ox * (T1 - T0);
      QovCox = Qsub0 / Coxeff;
      dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg)
                   + QovCox * dCoxeff_dVg;
      dQsub0_dVd = -T2 * dVgsteff_dVd;
      dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb + dVgsteff_dVb)
                   + QovCox * dCoxeff_dVb;

      // Gate-bias dependent delta Phis begins */
      if (paramPtr->k1ox <= 0.0)
      {
        Denomi = 0.25 * paramPtr->moin * Vtm;
        T0 = 0.5 * paramPtr->sqrtPhi;
      }
      else
      {
        Denomi = paramPtr->moin * Vtm * paramPtr->k1ox * paramPtr->k1ox;
        T0 = paramPtr->k1ox * paramPtr->sqrtPhi;
      }
      T1 = 2.0 * T0 + Vgsteff;

      DeltaPhi = Vtm * log(1.0 + T1 * Vgsteff / Denomi);
      dDeltaPhi_dVg = 2.0 * Vtm * (T1 -T0) / (Denomi + T1 * Vgsteff);
      dDeltaPhi_dVd = dDeltaPhi_dVg * dVgsteff_dVd;
      dDeltaPhi_dVb = dDeltaPhi_dVg * dVgsteff_dVb;
      // End of delta Phis

      T3 = 4.0 * (Vth - paramPtr->vfbzb - paramPtr->phi);
      Tox += Tox;
      if (T3 >= 0.0)
      {  T0 = (Vgsteff + T3) / Tox;
        dT0_dVd = (dVgsteff_dVd + 4.0 * dVth_dVd) / Tox;
        dT0_dVb = (dVgsteff_dVb + 4.0 * dVth_dVb) / Tox;
      }
      else
      {  T0 = (Vgsteff + 1.0e-20) / Tox;
        dT0_dVd = dVgsteff_dVd / Tox;
        dT0_dVb = dVgsteff_dVb / Tox;
      }
      tmp = exp(0.7 * log(T0));
      T1 = 1.0 + tmp;
      T2 = 0.7 * tmp / (T0 * Tox);
      Tcen = 1.9e-9 / T1;
      dTcen_dVg = -1.9e-9 * T2 / T1 /T1;
      dTcen_dVd = Tox * dTcen_dVg;
      dTcen_dVb = dTcen_dVd * dT0_dVb;
      dTcen_dVd *= dT0_dVd;
      dTcen_dVg *= dVgsteff_dVg;

      Ccen = CONSTEPSSI / Tcen;
      T0 = Cox / (Cox + Ccen);
      Coxeff = T0 * Ccen;
      T1 = -Ccen / Tcen;
      dCoxeff_dVg = T0 * T0 * T1;
      dCoxeff_dVd = dCoxeff_dVg * dTcen_dVd;
      dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
      dCoxeff_dVg *= dTcen_dVg;
      CoxWLcen = CoxWL * Coxeff / Cox;

      AbulkCV = Abulk0 * paramPtr->abulkCVfactor;
      dAbulkCV_dVb = paramPtr->abulkCVfactor * dAbulk0_dVb;
      VdsatCV = (Vgsteff - DeltaPhi) / AbulkCV;
      V4 = VdsatCV - Vds - CONSTDELTA_4;
      T0 = sqrt(V4 * V4 + 4.0 * CONSTDELTA_4 * VdsatCV);
      VdseffCV = VdsatCV - 0.5 * (V4 + T0);
      T1 = 0.5 * (1.0 + V4 / T0);
      T2 = CONSTDELTA_4 / T0;
      T3 = (1.0 - T1 - T2) / AbulkCV;
      T4 = T3 * ( 1.0 - dDeltaPhi_dVg);
      dVdseffCV_dVg = T4;
      dVdseffCV_dVd = T1;
      dVdseffCV_dVb = -T3 * VdsatCV * dAbulkCV_dVb;

      // Added to eliminate non-zero VdseffCV at Vds=0.0
      if (Vds == 0.0)
      {  VdseffCV = 0.0;
        dVdseffCV_dVg = 0.0;
        dVdseffCV_dVb = 0.0;
      }

      T0 = AbulkCV * VdseffCV;
      T1 = Vgsteff - DeltaPhi;
      T2 = 12.0 * (T1 - 0.5 * T0 + 1.0e-20);
      T3 = T0 / T2;
      T4 = 1.0 - 12.0 * T3 * T3;
      T5 = AbulkCV * (6.0 * T0 * (4.0 * T1 - T0) / (T2 * T2) - 0.5);
      T6 = T5 * VdseffCV / AbulkCV;

      qgate = qinoi = CoxWLcen * (T1 - T0 * (0.5 - T3));
      QovCox = qgate / Coxeff;
      Cgg1 = CoxWLcen * (T4 * (1.0 - dDeltaPhi_dVg) + T5 * dVdseffCV_dVg);
      Cgd1 = CoxWLcen * T5 * dVdseffCV_dVd + Cgg1
             * dVgsteff_dVd + QovCox * dCoxeff_dVd;
      Cgb1 = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
             + Cgg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
      Cgg1 = Cgg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;


      T7 = 1.0 - AbulkCV;
      T8 = T2 * T2;
      T9 = 12.0 * T7 * T0 * T0 / (T8 * AbulkCV);
      T10 = T9 * (1.0 - dDeltaPhi_dVg);
      T11 = -T7 * T5 / AbulkCV;
      T12 = -(T9 * T1 / AbulkCV + VdseffCV * (0.5 - T0 / T2));

      qbulk = CoxWLcen * T7 * (0.5 * VdseffCV - T0 * VdseffCV / T2);
      QovCox = qbulk / Coxeff;
      Cbg1 = CoxWLcen * (T10 + T11 * dVdseffCV_dVg);
      Cbd1 = CoxWLcen * T11 * dVdseffCV_dVd + Cbg1
             * dVgsteff_dVd + QovCox * dCoxeff_dVd;
      Cbb1 = CoxWLcen * (T11 * dVdseffCV_dVb + T12 * dAbulkCV_dVb)
             + Cbg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
      Cbg1 = Cbg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;

      if (model_.xpart > 0.5)
      {   // 0/100 partition
        qsrc = -CoxWLcen * (T1 / 2.0 + T0 / 4.0 - 0.5 * T0 * T0 / T2);
        QovCox = qsrc / Coxeff;
        T2 += T2;
        T3 = T2 * T2;
        T7 = -(0.25 - 12.0 * T0 * (4.0 * T1 - T0) / T3);
        T4 = -(0.5 + 24.0 * T0 * T0 / T3) * (1.0 - dDeltaPhi_dVg);
        T5 = T7 * AbulkCV;
        T6 = T7 * VdseffCV;

        Csg = CoxWLcen * (T4 + T5 * dVdseffCV_dVg);
        Csd = CoxWLcen * T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd
              + QovCox * dCoxeff_dVd;
        Csb = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
              + Csg * dVgsteff_dVb + QovCox * dCoxeff_dVb;
        Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
      }
      else if (model_.xpart < 0.5)
      {   // 40/60 partition
        T2 = T2 / 12.0;
        T3 = 0.5 * CoxWLcen / (T2 * T2);
        T4 = T1 * (2.0 * T0 * T0 / 3.0 + T1 * (T1 - 4.0
                                               * T0 / 3.0)) - 2.0 * T0 * T0 * T0 / 15.0;
        qsrc = -T3 * T4;
        QovCox = qsrc / Coxeff;
        T8 = 4.0 / 3.0 * T1 * (T1 - T0) + 0.4 * T0 * T0;
        T5 = -2.0 * qsrc / T2 - T3 * (T1 * (3.0 * T1 - 8.0
                                            * T0 / 3.0) + 2.0 * T0 * T0 / 3.0);
        T6 = AbulkCV * (qsrc / T2 + T3 * T8);
        T7 = T6 * VdseffCV / AbulkCV;

        Csg = T5 * (1.0 - dDeltaPhi_dVg) + T6 * dVdseffCV_dVg;
        Csd = Csg * dVgsteff_dVd + T6 * dVdseffCV_dVd
              + QovCox * dCoxeff_dVd;
        Csb = Csg * dVgsteff_dVb + T6 * dVdseffCV_dVb
              + T7 * dAbulkCV_dVb + QovCox * dCoxeff_dVb;
        Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
      }
      else
      {   // 50/50 partition
        qsrc = -0.5 * qgate;
        Csg = -0.5 * Cgg1;
        Csd = -0.5 * Cgd1;
        Csb = -0.5 * Cgb1;
      }

      qgate += Qac0 + Qsub0 - qbulk;
      qbulk -= (Qac0 + Qsub0);
      qdrn = -(qgate + qbulk + qsrc);

      Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
      Cbd = Cbd1 - dQsub0_dVd;
      Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;

      Cgg = Cgg1 - Cbg;
      Cgd = Cgd1 - Cbd;
      Cgb = Cgb1 - Cbb;

      Cgb *= dVbseff_dVb;
      Cbb *= dVbseff_dVb;
      Csb *= dVbseff_dVb;

      cggb = Cgg;
      cgsb = -(Cgg + Cgd + Cgb);
      cgdb = Cgd;
      cdgb = -(Cgg + Cbg + Csg);
      cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb + Csg + Csd + Csb);
      cddb = -(Cgd + Cbd + Csd);
      cbgb = Cbg;
      cbsb = -(Cbg + Cbd + Cbb);
      cbdb = Cbd;
      qinv = -qinoi;
    }  // End of CTM
  }

  finished:
  // Returning Values to Calling Routine
  // COMPUTE EQUIVALENT DRAIN CURRENT SOURCE

  // copy local "cdrain" variable over to instance variable cd.
  cd = cdrain;

  // charge storage elements:
  // bulk-drain and bulk-source depletion capacitances
  // czbd : zero bias drain junction capacitance
  // czbs : zero bias source junction capacitance
  // czbdsw: zero bias drain junction sidewall capacitance
  //         along field oxide
  // czbssw: zero bias source junction sidewall capacitance
  //         along field oxide
  // czbdswg: zero bias drain junction sidewall capacitance
  //          along gate side
  // czbsswg: zero bias source junction sidewall capacitance
  //          along gate side
  if (ChargeComputationNeeded)
  {
    czbd = unitAreaJctCapTemp * drainArea;
    czbs = unitAreaJctCapTemp * sourceArea;
    if (drainPerimeter < paramPtr->weff)
    {
      czbdswg = unitLengthGateSidewallJctCapTemp * drainPerimeter;
      czbdsw = 0.0;
    }
    else
    {
      czbdsw = unitLengthSidewallJctCapTemp
               * (drainPerimeter - paramPtr->weff);

      czbdswg = unitLengthGateSidewallJctCapTemp *  paramPtr->weff;
    }
    if (sourcePerimeter < paramPtr->weff)
    {
      czbssw = 0.0;
      czbsswg = unitLengthGateSidewallJctCapTemp * sourcePerimeter;
    }
    else
    {
      czbssw = unitLengthSidewallJctCapTemp
               * (sourcePerimeter - paramPtr->weff);
      czbsswg = unitLengthGateSidewallJctCapTemp *  paramPtr->weff;
    }

    MJ    = model_.bulkJctBotGradingCoeff;
    MJSW  = model_.bulkJctSideGradingCoeff;
    MJSWG = model_.bulkJctGateSideGradingCoeff;

    // Source Bulk Junction
    if (vbs == 0.0)
    {
      //*(ckt->CKTstate0 + iterI->qbs) = 0.0;
      qbs = 0.0;
      capbs = czbs + czbssw + czbsswg;
    }
    else if (vbs < 0.0)
    {
      if (czbs > 0.0)
      {
        arg = 1.0 - vbs / PhiBTemp;

        if (MJ == 0.5) sarg = 1.0 / sqrt(arg);
        else           sarg = exp(-MJ * log(arg));

        //*(ckt->CKTstate0 + iterI->qbs) =
        qbs = PhiBTemp * czbs * (1.0 - arg * sarg) / (1.0 - MJ);

        capbs = czbs * sarg;
      }
      else
      {
        //*(ckt->CKTstate0 + iterI->qbs) = 0.0;
        qbs = 0.0;
        capbs = 0.0;
      }

      if (czbssw > 0.0)
      {
        arg = 1.0 - vbs / PhiBSWTemp;
        if (MJSW == 0.5) sarg = 1.0 / sqrt(arg);
        else             sarg = exp(-MJSW * log(arg));

        //*(ckt->CKTstate0 + iterI->qbs) +=
        qbs += PhiBSWTemp * czbssw * (1.0 - arg * sarg) / (1.0 - MJSW);

        capbs += czbssw * sarg;
      }

      if (czbsswg > 0.0)
      {
        arg = 1.0 - vbs / PhiBSWGTemp;
        if (MJSWG == 0.5) sarg = 1.0 / sqrt(arg);
        else              sarg = exp(-MJSWG * log(arg));

        //*(ckt->CKTstate0 + iterI->qbs) +=
        qbs += PhiBSWGTemp * czbsswg * (1.0 - arg * sarg) / (1.0 - MJSWG);

        capbs += czbsswg * sarg;
      }

    }
    else
    {
      T0 = czbs + czbssw + czbsswg;
      T1 = vbs * (czbs * MJ / PhiBTemp + czbssw * MJSW
                  / PhiBSWTemp + czbsswg * MJSWG / PhiBSWGTemp);

      //*(ckt->CKTstate0 + iterI->
      qbs = vbs * (T0 + 0.5 * T1);
      capbs = T0 + T1;
    }

    // Drain Bulk Junction
    if (vbd == 0.0)
    {
      //*(ckt->CKTstate0 + iterI->qbd) = 0.0;
      qbd = 0.0;
      capbd = czbd + czbdsw + czbdswg;
    }
    else if (vbd < 0.0)
    {
      if (czbd > 0.0)
      {
        arg = 1.0 - vbd / PhiBTemp;
        if (MJ == 0.5) sarg = 1.0 / sqrt(arg);
        else           sarg = exp(-MJ * log(arg));

        //*(ckt->CKTstate0 + iterI->qbd) =
        qbd = PhiBTemp * czbd * (1.0 - arg * sarg) / (1.0 - MJ);
        capbd = czbd * sarg;
      }
      else
      {
        //*(ckt->CKTstate0 + iterI->qbd) = 0.0;
        qbd = 0.0;
        capbd = 0.0;
      }

      if (czbdsw > 0.0)
      {
        arg = 1.0 - vbd / PhiBSWTemp;
        if (MJSW == 0.5) sarg = 1.0 / sqrt(arg);
        else             sarg = exp(-MJSW * log(arg));

        //*(ckt->CKTstate0 + iterI->qbd) +=
        qbd += PhiBSWTemp * czbdsw * (1.0 - arg * sarg) / (1.0 - MJSW);
        capbd += czbdsw * sarg;
      }

      if (czbdswg > 0.0)
      {
        arg = 1.0 - vbd / PhiBSWGTemp;
        if (MJSWG == 0.5) sarg = 1.0 / sqrt(arg);
        else              sarg = exp(-MJSWG * log(arg));

        //*(ckt->CKTstate0 + iterI->qbd) +=
        qbd += PhiBSWGTemp * czbdswg * (1.0 - arg * sarg) / (1.0 - MJSWG);
        capbd += czbdswg * sarg;
      }
    }
    else
    {
      T0 = czbd + czbdsw + czbdswg;
      T1 = vbd * (czbd * MJ / PhiBTemp + czbdsw * MJSW
                  / PhiBSWTemp + czbdswg * MJSWG / PhiBSWGTemp);

      //*(ckt->CKTstate0 + iterI->qbd) = vbd * (T0 + 0.5 * T1);
      qbd = vbd * (T0 + 0.5 * T1);
      capbd = T0 + T1;
    }
  }

  // There is a spice3f5 convergence check that would happen here.
  // (line 2404) skipping...

  // In 3f5, loading a bunch of things into the state vector at this point.
  // (line 2433) skipping...

  // bulk and channel charge plus overlaps

  if (ChargeComputationNeeded)
  {
    // NQS begins
    if (nqsMod)
    {
      qcheq = -(qbulk + qgate);

      cqgb = -(cggb + cbgb);
      cqdb = -(cgdb + cbdb);
      cqsb = -(cgsb + cbsb);
      cqbb = -(cqgb + cqdb + cqsb);

      gtau_drift = fabs(paramPtr->tconst * qcheq) * ScalingFactor;
      T0 = paramPtr->leffCV * paramPtr->leffCV;
      gtau_diff = 16.0 * paramPtr->u0temp * model_.vtm / T0 * ScalingFactor;

      gtau =  gtau_drift + gtau_diff;
    }

    if (model_.capMod == 0)
    {
      if (vgd < 0.0)
      {
        cgdo = paramPtr->cgdo;
        qgdo = paramPtr->cgdo * vgd;
      }
      else
      {
        cgdo = paramPtr->cgdo;
        qgdo =  paramPtr->cgdo * vgd;
      }

      if (vgs < 0.0)
      {
        cgso = paramPtr->cgso;
        qgso = paramPtr->cgso * vgs;
      }
      else
      {
        cgso = paramPtr->cgso;
        qgso =  paramPtr->cgso * vgs;
      }
    }
    else if (model_.capMod == 1)
    {
      if (vgd < 0.0)
      {
        T1 = sqrt(1.0 - 4.0 * vgd / paramPtr->ckappa);
        cgdo = paramPtr->cgdo + paramPtr->weffCV * paramPtr->cgdl / T1;

        qgdo = paramPtr->cgdo * vgd - paramPtr->weffCV * 0.5
               * paramPtr->cgdl * paramPtr->ckappa * (T1 - 1.0);
      }
      else
      {
        cgdo = paramPtr->cgdo + paramPtr->weffCV * paramPtr->cgdl;
        qgdo = (paramPtr->weffCV * paramPtr->cgdl + paramPtr->cgdo) * vgd;
      }

      if (vgs < 0.0)
      {
        T1 = sqrt(1.0 - 4.0 * vgs / paramPtr->ckappa);
        cgso = paramPtr->cgso + paramPtr->weffCV * paramPtr->cgsl / T1;
        qgso = paramPtr->cgso * vgs - paramPtr->weffCV * 0.5
               * paramPtr->cgsl * paramPtr->ckappa * (T1 - 1.0);
      }
      else
      {
        cgso = paramPtr->cgso + paramPtr->weffCV * paramPtr->cgsl;
        qgso = (paramPtr->weffCV * paramPtr->cgsl + paramPtr->cgso) * vgs;
      }
    }
    else
    {
      T0 = vgd + CONSTDELTA_1;
      T1 = sqrt(T0 * T0 + 4.0 * CONSTDELTA_1);
      T2 = 0.5 * (T0 - T1);

      T3 = paramPtr->weffCV * paramPtr->cgdl;
      T4 = sqrt(1.0 - 4.0 * T2 / paramPtr->ckappa);
      cgdo = paramPtr->cgdo + T3 - T3 * (1.0 - 1.0 / T4) * (0.5 - 0.5 * T0 / T1);

      qgdo = (paramPtr->cgdo + T3) * vgd - T3 * (T2
                                                 + 0.5 * paramPtr->ckappa * (T4 - 1.0));

      T0 = vgs + CONSTDELTA_1;
      T1 = sqrt(T0 * T0 + 4.0 * CONSTDELTA_1);
      T2 = 0.5 * (T0 - T1);
      T3 = paramPtr->weffCV * paramPtr->cgsl;
      T4 = sqrt(1.0 - 4.0 * T2 / paramPtr->ckappa);
      cgso = paramPtr->cgso + T3 - T3 * (1.0 - 1.0 / T4) * (0.5 - 0.5 * T0 / T1);
      qgso = (paramPtr->cgso + T3) * vgs - T3 * (T2
                                                 + 0.5 * paramPtr->ckappa * (T4 - 1.0));
    }

    //cgdo = cgdo;
    //cgso = cgso;

    setupCapacitors_oldDAE ();
    setupCapacitors_newDAE ();

    //cqdef = cqcheq = 0.0;
    cqdef = 0.0;

    // set some state variables:
    qg = qgate;
    qd = qdrn - qbd;
    qb = qbulk + qbd + qbs;
    if (nqsMod) qcdump = qdef * ScalingFactor;

  } // end of ChargeComputationNeeded if statement.

  // store small signal parameters
  //if (ckt->CKTmode & MODEINITSMSIG) goto line1000;
  // Note: in 3f5, line1000 is at the end of the load, after
  //        the loads to the rhs and the matrix.  So it looks
  //        like this goto essentially means return.


  // Setting up a few currents for the RHS load:
  Idrain = drainConductance * Vddp;
  Isource = sourceConductance * Vssp;

  // Put this kludge in because the matrix load needs T1 but it is used
  // all over the place:
  T1global = T1;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
//
// Purpose       : This function sets up the primaray state variables into
//                 the primary state vector.
//
//                 These variables include qb, qg, qd and, in the
//                 event that nqsMod=1, qcdump and qcheq.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/09/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  bool bsuccess = true;
  double * staVec = extData.nextStaVectorRawPtr;
  double * stoVec = extData.nextStoVectorRawPtr;

  bsuccess = updateIntermediateVars ();

  // voltage drops:
  stoVec[li_store_vbs] = vbs;
  stoVec[li_store_vgs] = vgs;
  stoVec[li_store_vds] = vds;
  stoVec[li_store_vbd] = vbd;
  stoVec[li_store_von] = von;

  // transconductance:
  if (mode >= 0)
  {
    stoVec[li_store_gm] = gm;
  }
  else
  {
    stoVec[li_store_gm] = -gm;
  }
  stoVec[li_store_Vds] = Vds;
  stoVec[li_store_Vgs] = Vgs;
  stoVec[li_store_Vbs] = Vbs;
  stoVec[li_store_Vdsat] = Vdsat;
  stoVec[li_store_Vth] = Vth;
  stoVec[li_store_Gds] = gds;
  stoVec[li_store_Cgs] = CAPcgsb;
  stoVec[li_store_Cgd] = CAPcgdb;

  // intrinsic capacitors:
  staVec[li_state_qb] = qb;
  staVec[li_state_qg] = qg;
  staVec[li_state_qd] = qd;

  // parasitic capacitors:
  staVec[li_state_qbs] = qbs;
  staVec[li_state_qbd] = qbd;

  if( nqsMod )
  {
    staVec[li_state_qcheq] = qcheq;
    staVec[li_state_qcdump] = qcdump;
  }

  // if this is the first newton step of the first time step
  // of the transient simulation, we need to enforce that the
  // time derivatives w.r.t. charge are zero.  This is to maintain 3f5
  // compatibility.  ERK.

  // Note:  I think this kind of thing is enforced (or should be enforced,
  // anyway) at the time integration level.  So I'm not sure this step is
  // really needed, at least for new-DAE.  Derivatives out of the DCOP
  // are supposed to be zero at the first newton step.

  if (!(getSolverState().dcopFlag) && getSolverState().initTranFlag_ && getSolverState().newtonIter==0)
  {
    // re-set the state vector pointer that we are using to the "current"
    // pointer, rather than the "next" pointer.
    double * currStaVec = extData.currStaVectorRawPtr;

    // intrinsic capacitors:
    currStaVec[li_state_qb] = qb;
    currStaVec[li_state_qg] = qg;
    currStaVec[li_state_qd] = qd;

    // parasitic capacitors:
    currStaVec[li_state_qbs] = qbs;
    currStaVec[li_state_qbd] = qbd;

    if( nqsMod )
    {
      currStaVec[li_state_qcheq] = qcheq;
      currStaVec[li_state_qcdump] = qcdump;
    }
  }

  return bsuccess;
}

// Noise functions

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MOSFET_B3::Instance::StrongInversionNoiseEval
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/13/2015
//-----------------------------------------------------------------------------
double Instance::StrongInversionNoiseEval(double Vds, double freq, double temp)
{
  double cdAbs, esat, DelClm, EffFreq, N0, Nl, Vgst, Leff, Leffsq;
  double T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, Ssi;

  cdAbs = fabs(cd);
  Leff = paramPtr->leff - 2.0 * model_.lintnoi;
  Leffsq = Leff * Leff;
  esat = 2.0 * paramPtr->vsattemp / ueff;
  if(model_.em<=0.0) 
  {
    DelClm = 0.0; 
  }
  else 
  { 
    T0 = ((((Vds - Vdseff) / paramPtr->litl) 
          + model_.em) / esat);
    DelClm = paramPtr->litl * std::log (std::max(T0, N_MINLOG));
    if (DelClm < 0.0)       DelClm = 0.0;  /* bugfix */
  }

  EffFreq = pow(freq, model_.ef);
  T1 = CONSTQ * CONSTQ * 8.62e-5 * cdAbs * temp * ueff;
  T2 = 1.0e8 * EffFreq * Abulk * model_.cox * Leffsq;
  N0 = model_.cox * Vgsteff / CONSTQ;
  Nl = model_.cox * Vgsteff * (1.0 - AbovVgst2Vtm * Vdseff) / CONSTQ;

  T3 = model_.oxideTrapDensityA * std::log(std::max(((N0 + 2.0e14) / (Nl + 2.0e14)), N_MINLOG));
  T4 = model_.oxideTrapDensityB * (N0 - Nl);
  T5 = model_.oxideTrapDensityC * 0.5 * (N0 * N0 - Nl * Nl);

  T6 = 8.62e-5 * temp * cdAbs * cdAbs;
  T7 = 1.0e8 * EffFreq * Leffsq * paramPtr->weff;
  T8 = model_.oxideTrapDensityA + model_.oxideTrapDensityB * Nl 
     + model_.oxideTrapDensityC * Nl * Nl;
  
  T9 = (Nl + 2.0e14) * (Nl + 2.0e14);

  Ssi = T1 / T2 * (T3 + T4 + T5) + T6 / T7 * DelClm * T8 / T9;
  return Ssi;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MOSFET_B3::Instance::getNumNoiseSources
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/13/2015
//-----------------------------------------------------------------------------
int Instance::getNumNoiseSources () const
{
  return 4;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MOSFET_B3::Instance::setupNoiseSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 1/13/2015
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
// Function      : Xyce::Device::MOSFET_B3::Instance::getNoiseSources
// Purpose       :
// Special Notes :
//
//   Channel thermal and flicker noises are calculated based on the value
//   of model->BSIM3noiMod.
//   If model->BSIM3noiMod = 1,
//     Channel thermal noise = SPICE2 model
//     Flicker noise         = SPICE2 model
//   If model->BSIM3noiMod = 2,
//     Channel thermal noise = BSIM3 model
//     Flicker noise         = BSIM3 model
//   If model->BSIM3noiMod = 3,
//     Channel thermal noise = SPICE2 model
//     Flicker noise         = BSIM3 model
//   If model->BSIM3noiMod = 4,
//     Channel thermal noise = BSIM3 model
//     Flicker noise         = SPICE2 model
//   If model->BSIM3noiMod = 5,
//     Channel thermal noise = SPICE2 model with linear/sat fix
//     Flicker noise         = SPICE2 model
//   If model->BSIM3noiMod = 6,
//     Channel thermal noise = SPICE2 model with linear/sat fix
//     Flicker noise         = BSIM3 model
//
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 1/13/2015
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
  switch(model_.noiMod)
  {
    case 1:
    case 3:
      devSupport.noiseSupport(
         noiseData.noiseDens[2], noiseData.lnNoiseDens[2], THERMNOISE,
         numberParallel * 2.0 * fabs(gm + gds + gmbs) / 3.0, temp);

      break;
    case 5:
    case 6:
      {
      double vdsLocal = std::min( vds, vdsat);
      devSupport.noiseSupport(
         noiseData.noiseDens[2], noiseData.lnNoiseDens[2], THERMNOISE,
         numberParallel*(
            (3.0 - vdsLocal / vdsat) * fabs(gm + gds + gmbs) / 3.0
                         )
         , temp);
      }
      break;
    case 2:
    case 4:
      devSupport.noiseSupport(
         noiseData.noiseDens[2], noiseData.lnNoiseDens[2], THERMNOISE,
         numberParallel*(
            (ueff * fabs(qinv) / (paramPtr->leff * paramPtr->leff + ueff *fabs(qinv) * rds))
                         )
            , temp);
      break;
  }

  // flicker noise
  switch(model_.noiMod )
  {
    case 1:
    case 4:
    case 5:
      noiseData.noiseDens[3] = numberParallel*model_.kf *
        std::exp(model_.af * std::log(std::max(fabs(cd), N_MINLOG)))
        / (std::pow(noiseData.freq,model_.ef) * paramPtr->leff * paramPtr->leff *model_.cox);
      break;
    case 2:
    case 3:
    case 6:
      {
        double vdsLocal = vds;
        if (vdsLocal < 0.0)
        {
          vdsLocal = -vdsLocal;
        }
        double Ssi,Swi,T10,T11,T1;
        Ssi = StrongInversionNoiseEval(vdsLocal, noiseData.freq, temp);
        T10 =model_.oxideTrapDensityA * 8.62e-5 * temp;
        T11 = paramPtr->weff * paramPtr->leff * std::pow(noiseData.freq,model_.ef) * 4.0e36;
        Swi = T10 / T11 * cd * cd;
        T1 = Swi + Ssi;
        if (T1 > 0.0)
        {
          noiseData.noiseDens[3] = numberParallel*(Ssi * Swi) / T1;
        }
        else
        {
          noiseData.noiseDens[3] = 0.0;
        }
      }
      break;
  }
  noiseData.lnNoiseDens[3] = std::log(std::max(noiseData.noiseDens[3],N_MINLOG));
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 bsim3 instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/07/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  double * qVec = extData.daeQVectorRawPtr;
  double * dQdxdVp = extData.dQdxdVpVectorRawPtr;

  auxChargeCalculations ();

  double Qeqqg = 0.0;   // gate charge
  double Qeqqb = 0.0;   // bulk charge
  double Qeqqd = 0.0;   // drain charge
  double Qqdef = 0.0;   // nqs-related charge.
  double Qqcheq = 0.0;  // nqs-related charge.

  // These 3 vars are class variables, and are set up elsewhere.
  //double Qeqqg_Jdxp = 0.0; // limiter, related to gate  cap.
  //double Qeqqb_Jdxp = 0.0; // limiter, related to bulk  cap.
  //double Qeqqd_Jdxp = 0.0; // limiter, related to drain cap.

  if (model_.dtype > 0)
  {
    Qeqqg  = qg;
    Qeqqb  = qb;
    Qeqqd  = qd;
    Qqdef  = qcdump; // this needs to be fixed...
    Qqcheq = qcheq;
  }
  else  // need to convert these to charges.
  {
    Qeqqg  = -qg;
    Qeqqb  = -qb;
    Qeqqd  = -qd;
    Qqdef  = -qcdump;
    Qqcheq = -qcheq;
  }

  qVec[li_Gate] += Qeqqg*numberParallel;
  qVec[li_Bulk] += (Qeqqb)*numberParallel;
  qVec[li_DrainPrime] += (-(-Qeqqd))*numberParallel;
  qVec[li_SourcePrime] += (-(+ Qeqqg + Qeqqb + Qeqqd))*numberParallel;

  if( loadLeadCurrent )
  {
    double * leadQ = extData.nextLeadCurrQCompRawPtr;
    if (drainConductance == 0.0)
    {
      leadQ[li_branch_dev_id] = (-(-Qeqqd))*numberParallel;
    }
    if (sourceConductance == 0.0)
    {
      leadQ[li_branch_dev_is] = (-(Qeqqg + Qeqqb + Qeqqd))*numberParallel;
    }
    leadQ[li_branch_dev_ig] = Qeqqg*numberParallel;
    leadQ[li_branch_dev_ib] = (Qeqqb)*numberParallel;
  }

  if (nqsMod)
  {
    // 7 equ. for nqs modification. charge equation.
    qVec[li_Charge] += -(Qqcheq - Qqdef)*numberParallel;
  }

  //////////////////////////////////////////////////
  // limiting section:
  if (getDeviceOptions().voltageLimiterFlag)
  {
    // Need the following:
    //  Qeqqg_Jdxp
    //  Qeqqb_Jdxp
    //  Qeqqd_Jdxp
    if (model_.dtype < 0)
    {
      Qeqqg_Jdxp = -Qeqqg_Jdxp;
      Qeqqb_Jdxp = -Qeqqb_Jdxp;
      Qeqqd_Jdxp = -Qeqqd_Jdxp;
    }

    if (!origFlag)
    {
      dQdxdVp[li_Gate] += -Qeqqg_Jdxp*numberParallel;
      dQdxdVp[li_Bulk] += -(+Qeqqb_Jdxp)*numberParallel;
      dQdxdVp[li_DrainPrime] += (-Qeqqd_Jdxp)*numberParallel;
      dQdxdVp[li_SourcePrime] += (+Qeqqg_Jdxp+Qeqqb_Jdxp+Qeqqd_Jdxp) *numberParallel;
    } // orig flag.
  } // limiter flag

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::auxChargeCalculations
//
// Purpose       : This function does some final "cleanup" calculations
//                 having to do with the capacitors.
//
// Special Notes : About all this function really does is set up some
//                 voltlim terms, and some unused nqs stuff.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/18/04
//-----------------------------------------------------------------------------
bool Instance::auxChargeCalculations ()
{
  double T0, T1;

  if (!ChargeComputationNeeded)
  {
    sxpart = (1.0 - (dxpart = (mode > 0) ? 0.4 : 0.6));
    ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
    dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;

    if (nqsMod)
    {
      gtau = 16.0 * paramPtr->u0temp * model_.vtm
             / paramPtr->leffCV / paramPtr->leffCV * ScalingFactor;
    }
    else
    {
      gtau = 0.0;
    }
  }
  else  // ChargeComputation is needed
  {
    double vgb_orig = vgs_orig - vbs_orig;

    Qeqqg_Jdxp = 0.0;
    if (!origFlag)
    {
      Qeqqg_Jdxp = - CAPcggb * (vgb-vgb_orig)
                   + CAPcgdb * (vbd-vbd_orig)
                   + CAPcgsb * (vbs-vbs_orig);
    }

    Qeqqb_Jdxp = 0.0;
    if (!origFlag)
    {
      Qeqqb_Jdxp = - CAPcbgb * (vgb-vgb_orig)
                   + CAPcbdb * (vbd-vbd_orig)
                   + CAPcbsb * (vbs-vbs_orig);
    }

    Qeqqd_Jdxp = 0.0;
    if (!origFlag)
    {
      Qeqqd_Jdxp = - CAPcdgb * (vgb-vgb_orig)
                   + CAPcddb * (vbd-vbd_orig)
                   + CAPcdsb * (vbs-vbs_orig);
    }

    // Note:  nqs stuff is not yet finished, especially the voltage
    // limiting aspect.  For limiting, need to re-do T0 and the term
    // added to ceqqd.
    if (nqsMod)
    {
      DevelFatal(*this).in("Instance::auxChargeCalculations")
        << "Instance::auxChargeCalculations ()" << std::endl
        << " nqsMod=1 is not ready yet.  Re-run with nqsMod=0";

      T0 = ggtg * vgb - ggtd * vbd - ggts * vbs;
      ceqqg += T0;
      T1 = qdef * gtau;
      ceqqd -= dxpart * T0 + T1 * (ddxpart_dVg * vgb - ddxpart_dVd
                                   * vbd - ddxpart_dVs * vbs);

      //cqdef = cqcdump - gqdef * qdef;

      //if (!origFlag)
      //{
      //double tmp =   - (gcqgb * (vgb-vgb_orig)
      //- gcqdb * (vbd-vbd_orig)
      //- gcqsb * (vbs-vbs_orig)) + T0;
      //cqcheq += tmp;
      //cqcheq_Jdxp = tmp;
      //}
    }
  } // !ChargeComputationNeeded

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupCapacitors_newDAE ()
//
// Purpose       : This takes a lot of the individual capacitive terms and
//                 sums them together for loading into the Jacobian.
//
// Special Notes : This was extracted from updateIntermediateVars.
//                 Different variables are used, and nothing is multiplied
//                 by ag0 = solState.pdt.  The new-dae formulation handles
//                 all the 1/dt - related stuff up in the time integrator.
//
//                 NOTE: nqs not even close to being supported.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/23/04
//-----------------------------------------------------------------------------
bool Instance::setupCapacitors_newDAE ()
{
  if (mode > 0)
  {
    if (nqsMod == 0)
    {
      CAPcggb = (cggb + cgdo + cgso + paramPtr->cgbo );
      CAPcgdb = (cgdb - cgdo);
      CAPcgsb = (cgsb - cgso);

      CAPcdgb = (cdgb - cgdo);
      CAPcddb = (cddb + capbd + cgdo);
      CAPcdsb = cdsb;

      CAPcsgb = -(cggb + cbgb + cdgb + cgso);
      CAPcsdb = -(cgdb + cbdb + cddb);
      CAPcssb = (capbs + cgso - (cgsb + cbsb + cdsb));

      CAPcbgb = (cbgb - paramPtr->cgbo);
      CAPcbdb = (cbdb - capbd);
      CAPcbsb = (cbsb - capbs);
    }
    else  // nqsMode != 0
    {

    } // nqsMod
  }
  else
  {
    if (nqsMod == 0)
    {
      CAPcggb = (cggb + cgdo + cgso + paramPtr->cgbo );
      CAPcgdb = (cgsb - cgdo);
      CAPcgsb = (cgdb - cgso);

      CAPcdgb = -(cggb + cbgb + cdgb + cgdo);
      CAPcddb = (capbd + cgdo - (cgsb + cbsb + cdsb));
      CAPcdsb = -(cgdb + cbdb + cddb);

      CAPcsgb = (cdgb - cgso);
      CAPcsdb = cdsb;
      CAPcssb = (cddb + capbs + cgso);

      CAPcbgb = (cbgb - paramPtr->cgbo);
      CAPcbdb = (cbsb - capbd);
      CAPcbsb = (cbdb - capbs);
    }
    else  // nqsMode != 0
    {

    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupCapacitors_oldDAE ()
// Purpose       : Same as new-DAE version, but including pdt, essentially.
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/17/06
//-----------------------------------------------------------------------------
bool Instance::setupCapacitors_oldDAE ()
{
  double ag0 = getSolverState().pdt_;
  double T0 = 0.0;

  // ERK. 12/17/2006.
  // It is necessary to set ag0=0.0, because for the first time step out of
  // the DCOP, all the time derivatives are forced to be zero.  Thus, all
  // their derivatives should also be zero.  If it wasn't for that, then ag0
  // could always be pdt.  (it used to be, before the -jacobian_test capability).
  if (!(getSolverState().dcopFlag) && getSolverState().initTranFlag_ && getSolverState().newtonIter==0)
  {
    ag0 = 0.0;
  }

  if (mode > 0)
  {
    if (nqsMod == 0)
    {
      gcggb = (cggb + cgdo + cgso + paramPtr->cgbo ) * ag0;
      gcgdb = (cgdb - cgdo) * ag0;
      gcgsb = (cgsb - cgso) * ag0;

      gcdgb = (cdgb - cgdo) * ag0;
      gcddb = (cddb + capbd + cgdo) * ag0;
      gcdsb = cdsb * ag0;

      gcsgb = -(cggb + cbgb + cdgb + cgso) * ag0;
      gcsdb = -(cgdb + cbdb + cddb) * ag0;
      gcssb = (capbs + cgso - (cgsb + cbsb + cdsb)) * ag0;

      gcbgb = (cbgb - paramPtr->cgbo) * ag0;
      gcbdb = (cbdb - capbd) * ag0;
      gcbsb = (cbsb - capbs) * ag0;

      qgd = qgdo;
      qgs = qgso;
      qgb = paramPtr->cgbo * vgb;
      qgate += qgd + qgs + qgb;
      qbulk -= qgb;
      qdrn -= qgd;
      qsrc = -(qgate + qbulk + qdrn);

      ggtg = ggtd = ggtb = ggts = 0.0;
      sxpart = 0.6;
      dxpart = 0.4;
      ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
      dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
    }
    else  // nqsMode != 0
    {
      if (qcheq > 0.0)
      {
        T0 = paramPtr->tconst * qdef * ScalingFactor;
      }
      else
      {
        T0 = -paramPtr->tconst * qdef * ScalingFactor;
      }

      ggtg = gtg = T0 * cqgb;
      ggtd = gtd = T0 * cqdb;
      ggts = gts = T0 * cqsb;
      ggtb = gtb = T0 * cqbb;
      gqdef = ScalingFactor * ag0;

      gcqgb = cqgb * ag0;
      gcqdb = cqdb * ag0;
      gcqsb = cqsb * ag0;
      gcqbb = cqbb * ag0;

      gcggb = (cgdo + cgso + paramPtr->cgbo ) * ag0;
      gcgdb = -cgdo * ag0;
      gcgsb = -cgso * ag0;

      gcdgb = -cgdo * ag0;
      gcddb = (capbd + cgdo) * ag0;
      gcdsb = 0.0;

      gcsgb = -cgso * ag0;
      gcsdb = 0.0;
      gcssb = (capbs + cgso) * ag0;

      gcbgb = -paramPtr->cgbo * ag0;
      gcbdb = -capbd * ag0;
      gcbsb = -capbs * ag0;

      CoxWL = model_.cox * paramPtr->weffCV * paramPtr->leffCV;

      if (fabs(qcheq) <= 1.0e-5 * CoxWL)
      {
        if      (model_.xpart < 0.5) dxpart = 0.4;
        else if (model_.xpart > 0.5) dxpart = 0.0;
        else                               dxpart = 0.5;

        ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
      }
      else
      {
        dxpart = qdrn / qcheq;
        Cdd = cddb;
        Csd = -(cgdb + cddb + cbdb);
        ddxpart_dVd = (Cdd - dxpart * (Cdd + Csd)) / qcheq;
        Cdg = cdgb;
        Csg = -(cggb + cdgb + cbgb);
        ddxpart_dVg = (Cdg - dxpart * (Cdg + Csg)) / qcheq;

        Cds = cdsb;
        Css = -(cgsb + cdsb + cbsb);
        ddxpart_dVs = (Cds - dxpart * (Cds + Css)) / qcheq;

        ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);
      }
      sxpart = 1.0 - dxpart;
      dsxpart_dVd = -ddxpart_dVd;
      dsxpart_dVg = -ddxpart_dVg;
      dsxpart_dVs = -ddxpart_dVs;
      dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);

      qgd = qgdo;
      qgs = qgso;
      qgb = paramPtr->cgbo * vgb;
      qgate = qgd + qgs + qgb;
      qbulk = -qgb;
      qdrn = -qgd;
      qsrc = -(qgate + qbulk + qdrn);
    } // nqsMod
  }
  else
  {
    if (nqsMod == 0)
    {
      gcggb = (cggb + cgdo + cgso + paramPtr->cgbo ) * ag0;
      gcgdb = (cgsb - cgdo) * ag0;
      gcgsb = (cgdb - cgso) * ag0;

      gcdgb = -(cggb + cbgb + cdgb + cgdo) * ag0;
      gcddb = (capbd + cgdo - (cgsb + cbsb + cdsb)) * ag0;
      gcdsb = -(cgdb + cbdb + cddb) * ag0;

      gcsgb = (cdgb - cgso) * ag0;
      gcsdb = cdsb * ag0;
      gcssb = (cddb + capbs + cgso) * ag0;

      gcbgb = (cbgb - paramPtr->cgbo) * ag0;
      gcbdb = (cbsb - capbd) * ag0;
      gcbsb = (cbdb - capbs) * ag0;

      qgd = qgdo;
      qgs = qgso;
      qgb = paramPtr->cgbo * vgb;
      qgate += qgd + qgs + qgb;
      qbulk -= qgb;
      qsrc = qdrn - qgs;
      qdrn = -(qgate + qbulk + qsrc);

      ggtg = ggtd = ggtb = ggts = 0.0;
      sxpart = 0.4;
      dxpart = 0.6;
      ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
      dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
    }
    else  // nqsMode != 0
    {
      if (qcheq > 0.0)
      {
        T0 = paramPtr->tconst * qdef * ScalingFactor;
      }
      else
      {
        T0 = -paramPtr->tconst * qdef * ScalingFactor;
      }

      ggtg = gtg = T0 * cqgb;
      ggts = gtd = T0 * cqdb;
      ggtd = gts = T0 * cqsb;
      ggtb = gtb = T0 * cqbb;
      gqdef = ScalingFactor * ag0;

      gcqgb = cqgb * ag0;
      gcqdb = cqsb * ag0;
      gcqsb = cqdb * ag0;
      gcqbb = cqbb * ag0;

      gcggb = (cgdo + cgso + paramPtr->cgbo) * ag0;
      gcgdb = -cgdo * ag0;
      gcgsb = -cgso * ag0;

      gcdgb = -cgdo * ag0;
      gcddb = (capbd + cgdo) * ag0;
      gcdsb = 0.0;

      gcsgb = -cgso * ag0;
      gcsdb = 0.0;
      gcssb = (capbs + cgso) * ag0;

      gcbgb = -paramPtr->cgbo * ag0;
      gcbdb = -capbd * ag0;
      gcbsb = -capbs * ag0;

      CoxWL = model_.cox * paramPtr->weffCV * paramPtr->leffCV;

      if (fabs(qcheq) <= 1.0e-5 * CoxWL)
      {
        if      (model_.xpart < 0.5) sxpart = 0.4;
        else if (model_.xpart > 0.5) sxpart = 0.0;
        else                         sxpart = 0.5;

        dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
      }
      else
      {
        sxpart = qdrn / qcheq;
        Css = cddb;
        Cds = -(cgdb + cddb + cbdb);
        dsxpart_dVs = (Css - sxpart * (Css + Cds)) / qcheq;
        Csg = cdgb;
        Cdg = -(cggb + cdgb + cbgb);
        dsxpart_dVg = (Csg - sxpart * (Csg + Cdg)) / qcheq;

        Csd = cdsb;
        Cdd = -(cgsb + cdsb + cbsb);
        dsxpart_dVd = (Csd - sxpart * (Csd + Cdd)) / qcheq;

        dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);
      }

      dxpart = 1.0 - sxpart;
      ddxpart_dVd = -dsxpart_dVd;
      ddxpart_dVg = -dsxpart_dVg;
      ddxpart_dVs = -dsxpart_dVs;
      ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);

      qgd = qgdo;
      qgs = qgso;
      qgb = paramPtr->cgbo * vgb;
      qgate = qgd + qgs + qgb;
      qbulk = -qgb;
      qsrc = -qgs;
      qdrn = -(qgate + qbulk + qsrc);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 bsim3 instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/07/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * solVec = extData.nextSolVectorRawPtr;
  double * fVec = extData.daeFVectorRawPtr;
  double * dFdxdVp = extData.dFdxdVpVectorRawPtr;

  double coef(0.0);

  cdreq_Jdxp = 0.0;
  ceqbd_Jdxp = 0.0;
  ceqbs_Jdxp = 0.0;

  // Do a few auxilliary calculations, derived from  3f5.
  // load current vector
  if (mode >= 0)
  {
    Gm = gm;
    Gmbs = gmbs;
    FwdSum = Gm + Gmbs;
    RevSum = 0.0;

    cdreq =  model_.dtype * (cd);
    ceqbd = -model_.dtype * (csub);

    ceqbs = 0.0;

    gbbdp = -gbds;
    gbbsp = (gbds + gbgs + gbbs);

    gbdpg = gbgs;
    gbdpdp = gbds;
    gbdpb = gbbs;
    gbdpsp = -(gbdpg + gbdpdp + gbdpb);

    gbspg = 0.0;
    gbspdp = 0.0;
    gbspb = 0.0;
    gbspsp = 0.0;
  }
  else
  {
    Gm = -gm;
    Gmbs = -gmbs;
    FwdSum = 0.0;
    RevSum = -(Gm + Gmbs);

    cdreq = -model_.dtype * (cd);
    ceqbs = -model_.dtype * (csub);

    ceqbd = 0.0;

    gbbsp = -gbds;
    gbbdp = (gbds + gbgs + gbbs);

    gbdpg = 0.0;
    gbdpsp = 0.0;
    gbdpb = 0.0;
    gbdpdp = 0.0;

    gbspg = gbgs;
    gbspsp = gbds;
    gbspb = gbbs;
    gbspdp = -(gbspg + gbspsp + gbspb);
  }

  if (model_.dtype > 0)
  {
    ceqbs += (cbs);
    ceqbd += (cbd);
  }
  else
  {
    ceqbs -= (cbs);
    ceqbd -= (cbd);
  }

  if (drainConductance != 0.0)
  {
    fVec[li_Drain] += Idrain*numberParallel;
  }
  if (sourceConductance != 0.0)
  {
    fVec[li_Source] += Isource*numberParallel;
  }
  fVec[li_Bulk] += (ceqbs + ceqbd)*numberParallel;
  fVec[li_DrainPrime] += (-(ceqbd - cdreq)-Idrain)*numberParallel;
  fVec[li_SourcePrime] += (-(cdreq + ceqbs)-Isource)*numberParallel;

  // lead current support
  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    double * junctionV = extData.nextJunctionVCompRawPtr;
    if (drainConductance != 0.0)
    {
      leadF[li_branch_dev_id] = Idrain*numberParallel;
    }
    else
    {
      leadF[li_branch_dev_id] = (-(ceqbd - cdreq)-Idrain)*numberParallel;
    }
    if (sourceConductance != 0.0)
    {
      leadF[li_branch_dev_is] = Isource*numberParallel;
    }
    else
    {
      leadF[li_branch_dev_is] = (-(cdreq + ceqbs)-Isource)*numberParallel;
    }
    leadF[li_branch_dev_ig] = 0.0;
    leadF[li_branch_dev_ib] = (ceqbs + ceqbd)*numberParallel;

    junctionV[li_branch_dev_id] = solVec[li_Drain] - solVec[li_Source];
    junctionV[li_branch_dev_ig] = solVec[li_Gate] - solVec[li_Source];
    junctionV[li_branch_dev_is] = 0.0;
    junctionV[li_branch_dev_ib] = 0.0 ;
  }

  // Initial condition support

  if( getSolverState().dcopFlag && icVDSGiven )
  {
    coef = extData.nextSolVectorRawPtr[li_Ids];
    fVec[li_Drain] += coef;
    fVec[li_Source] += -coef;
    if( loadLeadCurrent )
    {
      double * leadF = extData.nextLeadCurrFCompRawPtr;
      leadF[li_branch_dev_id]= coef;
      leadF[li_branch_dev_is]= -coef;
    }
  }

  if( getSolverState().dcopFlag && icVGSGiven )
  {
    coef = extData.nextSolVectorRawPtr[li_Igs];
    fVec[li_Gate] += coef;
    fVec[li_Source] += -coef;
    if( loadLeadCurrent )
    {
      double * leadF = extData.nextLeadCurrFCompRawPtr;
      leadF[li_branch_dev_ig]= coef;
      leadF[li_branch_dev_is]= -coef;
    }
  }


  if( getSolverState().dcopFlag && icVBSGiven )
  {
    coef = extData.nextSolVectorRawPtr[li_Ibs];
    fVec[li_Bulk] += coef;
    fVec[li_Source] += -coef;
    if( loadLeadCurrent )
    {
      double * leadF = extData.nextLeadCurrFCompRawPtr;
      leadF[li_branch_dev_ib]= coef;
      leadF[li_branch_dev_is]= -coef;
    }
  }



  //////////////////////////////////////////////////
  // limiting section:
  if (getDeviceOptions().voltageLimiterFlag)
  {
    if (!origFlag)
    {
      if (mode >= 0)
      {
        // option 1
        double tmp = model_.dtype * (-gds * (vds-vds_orig) -
                                     Gm * (vgs-vgs_orig) -
                                     Gmbs * (vbs-vbs_orig));

        cdreq_Jdxp += tmp;
        cdreq      += tmp;

        tmp = -model_.dtype * (-gbds * (vds-vds_orig) -
                               gbgs * (vgs-vgs_orig) -
                               gbbs * (vbs-vbs_orig));
        ceqbd_Jdxp += tmp;
        ceqbd      += tmp;
      }
      else
      {
        // option 2
        double tmp = -model_.dtype * (gds * (vds-vds_orig) +
                                      Gm * (vgd-vgd_orig) +
                                      Gmbs * (vbd-vbd_orig));
        cdreq_Jdxp += tmp;
        cdreq      += tmp;

        tmp = -model_.dtype * (gbds * (vds-vds_orig) -
                               gbgs * (vgd-vgd_orig) -
                               gbbs * (vbd-vbd_orig));
        ceqbd_Jdxp += tmp;
        ceqbd      += tmp;
      }


      if (model_.dtype > 0)
      {
        ceqbs_Jdxp += (-gbs*(vbs-vbs_orig));
        ceqbs      += (-gbs*(vbs-vbs_orig));

        ceqbd_Jdxp += (-gbd*(vbd-vbd_orig));
        ceqbd      += (-gbd*(vbd-vbd_orig));
      }
      else
      {
        ceqbs_Jdxp -= (-gbs*(vbs-vbs_orig) );
        ceqbs      -= (-gbs*(vbs-vbs_orig) );

        ceqbd_Jdxp -= (-gbd*(vbd-vbd_orig) );
        ceqbd      -= (-gbd*(vbd-vbd_orig) );
      }
      dFdxdVp[li_Bulk] += -(ceqbs_Jdxp+ceqbd_Jdxp)*numberParallel;
      dFdxdVp[li_DrainPrime] += (ceqbd_Jdxp-cdreq_Jdxp)*numberParallel;
      dFdxdVp[li_SourcePrime] += (cdreq_Jdxp+ceqbs_Jdxp) *numberParallel;
    } // orig flag.
  } // voltage limiter flag

  // Row associated with icVBS
  if( getSolverState().dcopFlag && icVBSGiven )
  {
    // get the voltage drop from the previous solution
    double cVs = extData.nextSolVectorRawPtr[li_Source];
    double cVb = extData.nextSolVectorRawPtr[li_Bulk];
    fVec[li_Ibs] += (cVb - cVs - icVBS);
  }

  // Row associated with icVDS
  if( getSolverState().dcopFlag && icVDSGiven )
  {
    // get the voltage drop from the previous solution
    double cVd = extData.nextSolVectorRawPtr[li_Drain];
    double cVs = extData.nextSolVectorRawPtr[li_Source];
    fVec[li_Ids] += (cVd - cVs - icVDS);
  }

  // Row associated with icVGS
  if( getSolverState().dcopFlag && icVGSGiven )
  {
    // get the voltage drop from the previous solution
    double cVg = extData.nextSolVectorRawPtr[li_Gate];
    double cVs = extData.nextSolVectorRawPtr[li_Source];

    fVec[li_Igs] += (cVg - cVs - icVGS);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 bsim3 instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/07/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  {
    Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);

    // Row corresponding to the KCL for the drain node: NOTHING

    // Row corresponding to the KCL for the source node: NOTHING

    // Row corresponding to the KCL for the gate node:
    // Check this later.   ERK.  See the comments in the function
    // loadDAEdQdx, regarding ggtg, ggtb, etc.
    //
    // For now I am leaving out the gg terms, as they are zero when
    // nqsMod=0, which is always true.
    //
    dQdx[li_Gate][AGateEquGateNodeOffset]
      += (CAPcggb )*numberParallel;
    dQdx[li_Gate][AGateEquBulkNodeOffset]
      -= (CAPcggb + CAPcgdb + CAPcgsb )*numberParallel;
    dQdx[li_Gate][AGateEquDrainPrimeNodeOffset]
      += (CAPcgdb )*numberParallel;
    dQdx[li_Gate][AGateEquSourcePrimeNodeOffset]
      += (CAPcgsb )*numberParallel;

    // Row corresponding to the KCL for the bulk node:
    dQdx[li_Bulk][ABulkEquGateNodeOffset]
      += (CAPcbgb)*numberParallel;

    dQdx[li_Bulk][ABulkEquBulkNodeOffset]
      += (- CAPcbgb - CAPcbdb - CAPcbsb)*numberParallel;

    dQdx[li_Bulk][ABulkEquDrainPrimeNodeOffset]
      += (CAPcbdb)*numberParallel;

    dQdx[li_Bulk][ABulkEquSourcePrimeNodeOffset]
      += (CAPcbsb)*numberParallel;


    // Row corresponding to the KCL for the drain prime node:
    dQdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset]
      -= (+ CAPcdgb + CAPcddb + CAPcdsb )*numberParallel;

    dQdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset]
      += (CAPcdgb) *numberParallel;

    dQdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]
      += (+ CAPcddb )*numberParallel;

    dQdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]
      -= (- CAPcdsb) *numberParallel;

    // Row corresponding to the KCL for the source prime node:
    dQdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset]
      += (CAPcsgb) *numberParallel;

    dQdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset]
      -= (+ CAPcsgb + CAPcsdb + CAPcssb) *numberParallel;

    dQdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]
      -= (- CAPcsdb) *numberParallel;

    dQdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]
      += (+ CAPcssb) *numberParallel;

    // Row associated with the charge equation
    // This is currently not supported.
    if (nqsMod)
    {
      DevelFatal(*this).in("Instance::loadDAEdQdX")
        << " nqsMod=1 is not ready yet.  Re-run with nqsMod=0.";
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 bsim3 instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/07/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  // Row corresponding to the KCL for the drain node:
  dFdx [li_Drain][ADrainEquDrainNodeOffset]
    += drainConductance*numberParallel;
  dFdx [li_Drain][ADrainEquDrainPrimeNodeOffset]
    -= drainConductance*numberParallel;

  // Extra term for initial conditions on Vds in operating point
  if( getSolverState().dcopFlag && icVDSGiven )
  {
    dFdx[li_Drain][ADrainEquIdsOffset] += 1.0;
  }

  // Row corresponding to the KCL for the source node:
  dFdx [li_Source][ASourceEquSourceNodeOffset]
    += sourceConductance*numberParallel;
  dFdx [li_Source][ASourceEquSourcePrimeNodeOffset]
    -= sourceConductance*numberParallel;

  // Extra term for initial conditions on Vbs in operating point
  if( getSolverState().dcopFlag && icVBSGiven )
  {
    dFdx[li_Source][ASourceEquIbsOffset] -= 1.0;
  }
  // Extra term for initial conditions on Vds in operating point
  if( getSolverState().dcopFlag && icVDSGiven )
  {
    dFdx[li_Source][ASourceEquIdsOffset] -= 1.0;
  }
  // Extra term for initial conditions on Vgs in operating point
  if( getSolverState().dcopFlag && icVGSGiven )
  {
    dFdx[li_Source][ASourceEquIgsOffset] -= 1.0;
  }

  // Row corresponding to the KCL for the gate node: NOTHING
  // Check this later.   ERK.
  //
  //  The terms beginning with "gc" (gcggb, etc.) definately do NOT
  //  belong here.  I'm not sure aboug the gg terms.  On one hand, the
  //  rhs vector component for the gate node ONLY seems to take
  //  capacitive currents, which implies that all of these are capacitive
  //  conductances.  On the other hand, the gg terms do not appear to
  //  have been created by multiplying by ag0 = pdt = 1/dt.  Generally
  //  capacitive conductances are of the form g = C/dt, and the gg terms
  //  do not have this form.
  //
  //  For now, the gg issue is moot b/c those terms are only nonzero
  //  if nqsMod = 1, which is not a supported option.

  //  However, the gg
  //  terms (ggtg, ggtb, ggtd and ggts)
  //
  //(*JMatPtr)[li_Gate][AGateEquGateNodeOffset]
  //  += (gcggb - ggtg)*numberParallel;
  //(*JMatPtr)[li_Gate][AGateEquBulkNodeOffset]
  //  -= (gcggb + gcgdb + gcgsb + ggtb)*numberParallel;
  //(*JMatPtr)[li_Gate][AGateEquDrainPrimeNodeOffset]
  //  += (gcgdb - ggtd)*numberParallel;
  //(*JMatPtr)[li_Gate][AGateEquSourcePrimeNodeOffset]
  //  += (gcgsb - ggts)*numberParallel;

  // initial conditions on gate node
  // Extra term for initial conditions on Vgs in operating point
  if( getSolverState().dcopFlag && icVGSGiven )
  {
    dFdx[li_Gate][AGateEquIgsOffset] += 1.0;
  }

  // Row corresponding to the KCL for the bulk node:
  dFdx [li_Bulk][ABulkEquGateNodeOffset]
    += (- gbgs)*numberParallel;

  dFdx [li_Bulk][ABulkEquBulkNodeOffset]
    += (gbd + gbs - gbbs)*numberParallel;

  dFdx [li_Bulk][ABulkEquDrainPrimeNodeOffset]
    += (- gbd + gbbdp)*numberParallel;

  dFdx [li_Bulk][ABulkEquSourcePrimeNodeOffset]
    += (- gbs + gbbsp)*numberParallel;

  // Extra term for initial conditions on Vbs in operating point
  if( getSolverState().dcopFlag && icVBSGiven )
  {
    dFdx[li_Bulk][ABulkEquIbsOffset] += 1.0;
  }

  // Row corresponding to the KCL for the drain prime node:
  dFdx [li_DrainPrime][ADrainPrimeEquDrainNodeOffset]
    -= drainConductance*numberParallel;

  dFdx [li_DrainPrime][ADrainPrimeEquBulkNodeOffset]
    -= (gbd - Gmbs - dxpart*ggtb
        - T1global*ddxpart_dVb - gbdpb)*numberParallel;

  dFdx [li_DrainPrime][ADrainPrimeEquGateNodeOffset]
    += (Gm + dxpart*ggtg + T1global*ddxpart_dVg + gbdpg)
    *numberParallel;

  dFdx [li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]
    += (drainConductance + gds + gbd + RevSum + dxpart*ggtd
        + T1global*ddxpart_dVd + gbdpdp)*numberParallel;

  dFdx [li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]
    -= (gds + FwdSum - dxpart*ggts - T1global*ddxpart_dVs - gbdpsp)
    *numberParallel;

  // Row corresponding to the KCL for the source prime node:
  dFdx [li_SourcePrime][ASourcePrimeEquGateNodeOffset]
    += (- Gm + sxpart*ggtg + T1global*dsxpart_dVg + gbspg)
    *numberParallel;

  dFdx [li_SourcePrime][ASourcePrimeEquBulkNodeOffset]
    -= (gbs + Gmbs - sxpart*ggtb
        - T1global*dsxpart_dVb - gbspb)*numberParallel;

  dFdx [li_SourcePrime][ASourcePrimeEquSourceNodeOffset]
    -= sourceConductance*numberParallel;

  dFdx [li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]
    -= (gds + RevSum - sxpart*ggtd - T1global*dsxpart_dVd - gbspdp)
    *numberParallel;

  dFdx [li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]
    += (sourceConductance + gds + gbs + FwdSum + sxpart*ggts
        + T1global*dsxpart_dVs + gbspsp)*numberParallel;

  // Row associated with the charge equation
  if (nqsMod)
  {
    DevelFatal(*this).in("Instance::loadDAEdFdx")
      << " nqsMod=1 is not ready yet.  Re-run with nqsMod=0.";
  }

  // Initial condition rows
  // Row associated with icVBS
  if( icVBSGiven )
  {
    if( getSolverState().dcopFlag  )
    {
      dFdx[li_Ibs][icVBSEquVbOffset] += 1.0;
      dFdx[li_Ibs][icVBSEquVsOffset] -= 1.0;
    }
    else
    {
      dFdx[li_Ibs][icVBSEquIbsOffset] += 1.0;
    }
  }

  // Row associated with icVDS
  if( icVDSGiven )
  {
    if( getSolverState().dcopFlag  )
    {
      dFdx[li_Ids][icVDSEquVdOffset] += 1.0;
      dFdx[li_Ids][icVDSEquVsOffset] -= 1.0;
    }
    else
    {
      dFdx[li_Ids][icVDSEquIdsOffset] += 1.0;
    }
  }

  // Row associated with icVGS
  if( icVGSGiven )
  {
    if( getSolverState().dcopFlag  )
    {
      dFdx[li_Igs][icVGSEquVgOffset] += 1.0;
      dFdx[li_Igs][icVGSEquVsOffset] -= 1.0;
    }
    else
    {
      dFdx[li_Igs][icVGSEquIgsOffset] += 1.0;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/10/02
//-----------------------------------------------------------------------------
bool Instance::setIC ()
{
  bool bsuccess = true;

  if( icVBSGiven )
  {
    extData.currStoVectorRawPtr[li_store_vbs] = icVBS;
    extData.nextStoVectorRawPtr[li_store_vbs] = icVBS;
  }

  if( icVDSGiven )
  {
    extData.currStoVectorRawPtr[li_store_vds] = icVDS;
    extData.nextStoVectorRawPtr[li_store_vds] = icVDS;
  }

  if( icVGSGiven )
  {
    extData.currStoVectorRawPtr[li_store_vgs] = icVGS;
    extData.nextStoVectorRawPtr[li_store_vgs] = icVGS;
  }

  return bsuccess;
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
  std::string msg;
  cox = 3.453133e-11 / tox;
  if (!given("TOXM"))    toxm = tox;
  if (!given("DSUB"))  dsub = drout;
  if (!given("LLC"))  Llc = Ll;
  if (!given("LWC"))  Lwc = Lw;
  if (!given("LWLC")) Lwlc = Lwl;
  if (!given("WLC"))  Wlc = Wl;
  if (!given("WWL")) Wwlc = Wwl;
  if (!given("WWLC")) Wwlc = Wwl;
  if (!given("DWC"))  dwc = Wint;
  if (!given("DLC"))  dlc = Lint;

  if (!given("CF"))
  {
    double C1 = 2.0 * CONSTEPSOX;
    double C5 = M_PI;
    double C2 = 1.0 + (0.4e-6 / tox);
    double C3 = log(C2);
    cf = C1*C3/C5;
  }

  if (!given("CGDO"))
  {
    if (given("DLC") && (dlc > 0.0)) cgdo = dlc * cox - cgdl ;
    else                         cgdo = 0.6 * xj * cox;
  }

  if (!given("CGSO"))
  {
    if (given("DLC") && (dlc > 0.0)) cgso = dlc * cox - cgsl ;
    else                         cgso = 0.6 * xj * cox;
  }

  if (!given("CGBO")) cgbo = 2.0 * dwc * cox;

  if (!given("CJSWG"))
    unitLengthGateSidewallJctCap = unitLengthSidewallJctCap ;

  if (!given("PBSWG"))
    GatesidewallJctPotential = sidewallJctPotential;

  if (!given("MJSWG"))
    bulkJctGateSideGradingCoeff = bulkJctSideGradingCoeff;


  // More initializations:  taken from b3temp.c:
  if (bulkJctPotential < 0.1)
  {
    bulkJctPotential = 0.1;
    UserWarning(*this) << "Given pb is less than 0.1. Pb is set to 0.1.";
  }

  if (sidewallJctPotential < 0.1)
  {
    sidewallJctPotential = 0.1;
    UserWarning(*this) << "Given pbsw is less than 0.1. Pbsw is set to 0.1.";
  }

  if (GatesidewallJctPotential < 0.1)
  {
    GatesidewallJctPotential = 0.1;
    UserWarning(*this) << "Given pbswg is less than 0.1. Pbswg is set to 0.1.";
  }

  vcrit   = CONSTvt0 * log(CONSTvt0 / (CONSTroot2 * 1.0e-14));
  factor1 = sqrt(CONSTEPSSI / CONSTEPSOX * tox);

  Vtm0 = CONSTKoverQ * tnom;
  Eg0  = CONSTEg0 - CONSTalphaEg * tnom * tnom / (tnom + CONSTbetaEg);
  ni   = CONSTNi0 * (tnom / CONSTREFTEMP) * sqrt(tnom / CONSTREFTEMP)
         * exp(21.5565981 - Eg0 / (2.0 * Vtm0));

  // If there are any time dependent parameters, set their values at for
  // the current time.

  // We have changed model parameters (maybe) and so all size dependent params
  // we may have stored may be invalid.  Clear them
  clearTemperatureData();
  
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
// Creation Date : 11/14/00
//-----------------------------------------------------------------------------
  Model::Model(
  const Configuration & configuration,
    const ModelBlock &    MB,
    const FactoryBlock &  factory_block)
    : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    modType       (0),
    dtype         (CONSTNMOS),
    mobMod        (0),
    capMod        (0),
    noiMod        (1),
    binUnit       (0),
    paramChk      (0),
    version       ("3.2.2"),
    model_l  (0.0),
    model_w  (0.0),
    tox           (0.0),
    toxm          (0.0),
    cdsc          (0.0),
    cdscb         (0.0),
    cdscd         (0.0),
    cit           (0.0),
    nfactor       (0.0),
    xj            (0.0),
    vsat          (0.0),
    at            (0.0),
    a0            (0.0),
    ags           (0.0),
    a1            (0.0),
    a2            (0.0),
    keta          (0.0),
    nsub          (0.0),
    npeak         (0.0),
    ngate         (0.0),
    gamma1        (0.0),
    gamma2        (0.0),
    vbx           (0.0),
    vbm           (0.0),
    xt            (0.0),
    k1            (0.0),
    kt1           (0.0),
    kt1l          (0.0),
    kt2           (0.0),
    k2            (0.0),
    k3            (0.0),
    k3b           (0.0),
    w0            (0.0),
    nlx           (0.0),
    dvt0          (0.0),
    dvt1          (0.0),
    dvt2          (0.0),
    dvt0w         (0.0),
    dvt1w         (0.0),
    dvt2w         (0.0),
    drout         (0.0),
    dsub          (0.0),
    vth0          (0.0),
    ua            (0.0),
    ua1           (0.0),
    ub            (0.0),
    ub1           (0.0),
    uc            (0.0),
    uc1           (0.0),
    u0            (0.0),
    ute           (0.0),
    voff          (0.0),
    delta         (0.0),
    rdsw          (0.0),
    prwg          (0.0),
    prwb          (0.0),
    prt           (0.0),
    eta0          (0.0),
    etab          (0.0),
    pclm          (0.0),
    pdibl1        (0.0),
    pdibl2        (0.0),
    pdiblb        (0.0),
    pscbe1        (0.0),
    pscbe2        (0.0),
    pvag          (0.0),
    wr            (0.0),
    dwg           (0.0),
    dwb           (0.0),
    b0            (0.0),
    b1            (0.0),
    alpha0        (0.0),
    alpha1        (0.0),
    beta0         (0.0),
    ijth          (0.0),
    vfb           (0.0),
    elm           (0.0),
    cgsl          (0.0),
    cgdl          (0.0),
    ckappa        (0.0),
    cf            (0.0),
    vfbcv         (0.0),
    clc           (0.0),
    cle           (0.0),
    dwc           (0.0),
    dlc           (0.0),
    noff          (0.0),
    voffcv        (0.0),
    acde          (0.0),
    moin          (0.0),
    tcj           (0.0),
    tcjsw         (0.0),
    tcjswg        (0.0),
    tpb           (0.0),
    tpbsw         (0.0),
    tpbswg        (0.0),
    lcdsc         (0.0),
    lcdscb        (0.0),
    lcdscd        (0.0),
    lcit          (0.0),
    lnfactor      (0.0),
    lxj           (0.0),
    lvsat         (0.0),
    lat           (0.0),
    la0           (0.0),
    lags          (0.0),
    la1           (0.0),
    la2           (0.0),
    lketa         (0.0),
    lnsub         (0.0),
    lnpeak        (0.0),
    lngate        (0.0),
    lgamma1       (0.0),
    lgamma2       (0.0),
    lvbx          (0.0),
    lvbm          (0.0),
    lxt           (0.0),
    lk1           (0.0),
    lkt1          (0.0),
    lkt1l         (0.0),
    lkt2          (0.0),
    lk2           (0.0),
    lk3           (0.0),
    lk3b          (0.0),
    lw0           (0.0),
    lnlx          (0.0),
    ldvt0         (0.0),
    ldvt1         (0.0),
    ldvt2         (0.0),
    ldvt0w        (0.0),
    ldvt1w        (0.0),
    ldvt2w        (0.0),
    ldrout        (0.0),
    ldsub         (0.0),
    lvth0         (0.0),
    lua           (0.0),
    lua1          (0.0),
    lub           (0.0),
    lub1          (0.0),
    luc           (0.0),
    luc1          (0.0),
    lu0           (0.0),
    lute          (0.0),
    lvoff         (0.0),
    ldelta        (0.0),
    lrdsw         (0.0),
    lprwg         (0.0),
    lprwb         (0.0),
    lprt          (0.0),
    leta0         (0.0),
    letab         (0.0),
    lpclm         (0.0),
    lpdibl1       (0.0),
    lpdibl2       (0.0),
    lpdiblb       (0.0),
    lpscbe1       (0.0),
    lpscbe2       (0.0),
    lpvag         (0.0),
    lwr           (0.0),
    ldwg          (0.0),
    ldwb          (0.0),
    lb0           (0.0),
    lb1           (0.0),
    lalpha0       (0.0),
    lalpha1       (0.0),
    lbeta0        (0.0),
    lvfb          (0.0),
    lelm          (0.0),
    lcgsl         (0.0),
    lcgdl         (0.0),
    lckappa       (0.0),
    lcf           (0.0),
    lclc          (0.0),
    lcle          (0.0),
    lvfbcv        (0.0),
    lnoff         (0.0),
    lvoffcv       (0.0),
    lacde         (0.0),
    lmoin         (0.0),
    wcdsc         (0.0),
    wcdscb        (0.0),
    wcdscd        (0.0),
    wcit          (0.0),
    wnfactor      (0.0),
    wxj           (0.0),
    wvsat         (0.0),
    wat           (0.0),
    wa0           (0.0),
    wags          (0.0),
    wa1           (0.0),
    wa2           (0.0),
    wketa         (0.0),
    wnsub         (0.0),
    wnpeak        (0.0),
    wngate        (0.0),
    wgamma1       (0.0),
    wgamma2       (0.0),
    wvbx          (0.0),
    wvbm          (0.0),
    wxt           (0.0),
    wk1           (0.0),
    wkt1          (0.0),
    wkt1l         (0.0),
    wkt2          (0.0),
    wk2           (0.0),
    wk3           (0.0),
    wk3b          (0.0),
    ww0           (0.0),
    wnlx          (0.0),
    wdvt0         (0.0),
    wdvt1         (0.0),
    wdvt2         (0.0),
    wdvt0w        (0.0),
    wdvt1w        (0.0),
    wdvt2w        (0.0),
    wdrout        (0.0),
    wdsub         (0.0),
    wvth0         (0.0),
    wua           (0.0),
    wua1          (0.0),
    wub           (0.0),
    wub1          (0.0),
    wuc           (0.0),
    wuc1          (0.0),
    wu0           (0.0),
    wute          (0.0),
    wvoff         (0.0),
    wdelta        (0.0),
    wrdsw         (0.0),
    wprwg         (0.0),
    wprwb         (0.0),
    wprt          (0.0),
    weta0         (0.0),
    wetab         (0.0),
    wpclm         (0.0),
    wpdibl1       (0.0),
    wpdibl2       (0.0),
    wpdiblb       (0.0),
    wpscbe1       (0.0),
    wpscbe2       (0.0),
    wpvag         (0.0),
    wwr           (0.0),
    wdwg          (0.0),
    wdwb          (0.0),
    wb0           (0.0),
    wb1           (0.0),
    walpha0       (0.0),
    walpha1       (0.0),
    wbeta0        (0.0),
    wvfb          (0.0),
    welm          (0.0),
    wcgsl         (0.0),
    wcgdl         (0.0),
    wckappa       (0.0),
    wcf           (0.0),
    wclc          (0.0),
    wcle          (0.0),
    wvfbcv        (0.0),
    wnoff         (0.0),
    wvoffcv       (0.0),
    wacde         (0.0),
    wmoin         (0.0),
    pcdsc         (0.0),
    pcdscb        (0.0),
    pcdscd        (0.0),
    pcit          (0.0),
    pnfactor      (0.0),
    pxj           (0.0),
    pvsat         (0.0),
    pat           (0.0),
    pa0           (0.0),
    pags          (0.0),
    pa1           (0.0),
    pa2           (0.0),
    pketa         (0.0),
    pnsub         (0.0),
    pnpeak        (0.0),
    pngate        (0.0),
    pgamma1       (0.0),
    pgamma2       (0.0),
    pvbx          (0.0),
    pvbm          (0.0),
    pxt           (0.0),
    pk1           (0.0),
    pkt1          (0.0),
    pkt1l         (0.0),
    pkt2          (0.0),
    pk2           (0.0),
    pk3           (0.0),
    pk3b          (0.0),
    pw0           (0.0),
    pnlx          (0.0),
    pdvt0         (0.0),
    pdvt1         (0.0),
    pdvt2         (0.0),
    pdvt0w        (0.0),
    pdvt1w        (0.0),
    pdvt2w        (0.0),
    pdrout        (0.0),
    pdsub         (0.0),
    pvth0         (0.0),
    pua           (0.0),
    pua1          (0.0),
    pub           (0.0),
    pub1          (0.0),
    puc           (0.0),
    puc1          (0.0),
    pu0           (0.0),
    pute          (0.0),
    pvoff         (0.0),
    pdelta        (0.0),
    prdsw         (0.0),
    pprwg         (0.0),
    pprwb         (0.0),
    pprt          (0.0),
    peta0         (0.0),
    petab         (0.0),
    ppclm         (0.0),
    ppdibl1       (0.0),
    ppdibl2       (0.0),
    ppdiblb       (0.0),
    ppscbe1       (0.0),
    ppscbe2       (0.0),
    ppvag         (0.0),
    pwr           (0.0),
    pdwg          (0.0),
    pdwb          (0.0),
    pb0           (0.0),
    pb1           (0.0),
    palpha0       (0.0),
    palpha1       (0.0),
    pbeta0        (0.0),
    pvfb          (0.0),
    pelm          (0.0),
    pcgsl         (0.0),
    pcgdl         (0.0),
    pckappa       (0.0),
    pcf           (0.0),
    pclc          (0.0),
    pcle          (0.0),
    pvfbcv        (0.0),
    pnoff         (0.0),
    pvoffcv       (0.0),
    pacde         (0.0),
    pmoin         (0.0),
    tnom          (getDeviceOptions().tnom),
    cgso          (0.0),
    cgdo          (0.0),
    cgbo          (0.0),
    xpart         (0.0),
    cFringOut     (0.0),
    cFringMax     (0.0),
    sheetResistance              (0.0),
    jctSatCurDensity             (0.0),
    jctSidewallSatCurDensity     (0.0),
    bulkJctPotential             (0.0),
    bulkJctBotGradingCoeff       (0.0),
    bulkJctSideGradingCoeff      (0.0),
    bulkJctGateSideGradingCoeff  (0.0),
    sidewallJctPotential         (0.0),
    GatesidewallJctPotential     (0.0),
    unitAreaJctCap               (0.0),
    unitLengthSidewallJctCap     (0.0),
    unitLengthGateSidewallJctCap (0.0),
    jctEmissionCoeff             (0.0),
    jctTempExponent              (0.0),
    Lint     (0.0),
    Ll       (0.0),
    Llc      (0.0),
    Lln      (0.0),
    Lw       (0.0),
    Lwc      (0.0),
    Lwn      (0.0),
    Lwl      (0.0),
    Lwlc     (0.0),
    Lmin     (0.0),
    Lmax     (0.0),
    Wint     (0.0),
    Wl       (0.0),
    Wlc      (0.0),
    Wln      (0.0),
    Ww       (0.0),
    Wwc      (0.0),
    Wwn      (0.0),
    Wwl      (0.0),
    Wwlc     (0.0),
    Wmin     (0.0),
    Wmax     (0.0),
    vtm      (0.0),
    cox      (0.0),
    cof1     (0.0),
    cof2     (0.0),
    cof3     (0.0),
    cof4     (0.0),
    vcrit    (0.0),
    factor1  (0.0),
    PhiB     (0.0),
    PhiBSW   (0.0),
    PhiBSWG  (0.0),
    oxideTrapDensityA (0.0),
    oxideTrapDensityB (0.0),
    oxideTrapDensityC (0.0),
    em       (0.0),
    ef       (0.0),
    af       (0.0),
    kf       (0.0),
    lintnoi  (0.0),
    npeakGiven     (0.0),
    gamma1Given    (0.0),
    gamma2Given    (0.0),
    k1Given        (0.0),
    k2Given        (0.0),
    nsubGiven      (0.0),
    xtGiven        (0.0),
    vbxGiven       (0.0),
    vbmGiven       (0.0),
    vfbGiven       (0.0),
    vth0Given      (0.0),
    Vtm0            (0.0),
    Eg0             (0.0),
    ni              (0.0)
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
  setDefaultParams();

  // Set params according to .model line and constant defaults from metadata:
  setModParams(MB.params);

  // Set any non-constant parameter defaults:
  if (!given("TNOM"))
    tnom = getDeviceOptions().tnom;

  if (!given("TNOM") && tnom == 0.0)
    UserFatal(*this) << "TNOM is zero";

// Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  if (!given("VTH0"))
    vth0 = (dtype == CONSTNMOS) ? 0.7 : -0.7;
  if (!given("UC"))
    uc  = (mobMod == 3) ? -0.0465 : -0.0465e-9;
  if (!given("UC1"))
    uc1 = (mobMod == 3) ? -0.056 : -0.056e-9;
  if (!given("U0"))
    u0  = (dtype == CONSTNMOS) ? 0.067 : 0.025;
  if (!given("NOIA"))
  {
    if (dtype == CONSTNMOS)
      oxideTrapDensityA = 1e20;
    else
      oxideTrapDensityA = 9.9e18;
  }
  if (!given("NOIB"))
  {
    if (dtype == CONSTNMOS)
      oxideTrapDensityB = 5e4;
    else
      oxideTrapDensityB = 2.4e3;
  }
  if (!given("NOIC"))
  {
    if (dtype == CONSTNMOS)
      oxideTrapDensityC = -1.4e-12;
    else
      oxideTrapDensityC = 1.4e-12;
  }

  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/00
//-----------------------------------------------------------------------------
Model::~Model ()
{
  std::list<SizeDependParam<double> *>::iterator it_dpL =
    sizeDependParamList.begin();
  std::list<SizeDependParam<double> *>::iterator end_dpL =
    sizeDependParamList.end();
  for( ; it_dpL != end_dpL; ++it_dpL )
    delete (*it_dpL);

  sizeDependParamList.clear ();

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
std::ostream &Model::printOutInstances (std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << std::endl;
  os << "    name     modelName  Parameters" << std::endl;

  for (i=0, iter=first; iter!=last; ++iter,++i)
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



//----------------------------------------------------------------------------
// Function       : Model::clearTemperatureData
//
// Purpose        : This is mainly here to delete rid of the size
//                  dependent parameters, which are also temperature dependent.
//
// Special Notes  : This is called right before the circuit temperature is
//                  changed.
//
// Scope          : public
// Creator        : Eric R. Keiter, 9233, computation sciences
// Creation Date  : 10/26/2004
//----------------------------------------------------------------------------
bool Model::clearTemperatureData ()
{
  std::list<SizeDependParam<double> *>::iterator it_dpL =
    sizeDependParamList.begin();
  std::list<SizeDependParam<double> *>::iterator end_dpL =
    sizeDependParamList.end();
  for( ; it_dpL != end_dpL; ++it_dpL )
    delete (*it_dpL);

  sizeDependParamList.clear ();

  return true;
}

//-----------------------------------------------------------------------------
// MOSFET_B3 Master functions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/02/08
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  bool bsuccess = true;
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & mi = *(*it);

    bool btmp = mi.updateIntermediateVars ();
    bsuccess = bsuccess && btmp;

    // voltage drops:
    double * stoVec = mi.extData.nextStoVectorRawPtr;
    stoVec[mi.li_store_vbs] = mi.vbs;
    stoVec[mi.li_store_vgs] = mi.vgs;
    stoVec[mi.li_store_vds] = mi.vds;
    stoVec[mi.li_store_vbd] = mi.vbd;
    stoVec[mi.li_store_von] = mi.von;

    // transconductance:
    if (mi.mode >= 0)
    {
      stoVec[mi.li_store_gm] = mi.gm;
    }
    else
    {
      stoVec[mi.li_store_gm] = -mi.gm;
    }
    stoVec[mi.li_store_Vds] = mi.Vds;
    stoVec[mi.li_store_Vgs] = mi.Vgs;
    stoVec[mi.li_store_Vbs] = mi.Vbs;
    stoVec[mi.li_store_Vdsat] = mi.Vdsat;
    stoVec[mi.li_store_Vth] = mi.Vth;

    stoVec[mi.li_store_Gds] = mi.gds;
    stoVec[mi.li_store_Cgs] = mi.CAPcgsb;
    stoVec[mi.li_store_Cgd] = mi.CAPcgdb;

    // intrinsic capacitors:
    staVec[mi.li_state_qb] = mi.qb;
    staVec[mi.li_state_qg] = mi.qg;
    staVec[mi.li_state_qd] = mi.qd;

    // parasitic capacitors:
    staVec[mi.li_state_qbs] = mi.qbs;
    staVec[mi.li_state_qbd] = mi.qbd;

    if( mi.nqsMod )
    {
      staVec[mi.li_state_qcheq] = mi.qcheq;
      staVec[mi.li_state_qcdump] = mi.qcdump;
    }

    // if this is the first newton step of the first time step
    // of the transient simulation, we need to enforce that the
    // time derivatives w.r.t. charge are zero.  This is to maintain 3f5
    // compatibility.  ERK.

    // Note:  I think this kind of thing is enforced (or should be enforced,
    // anyway) at the time integration level.  So I'm not sure this step is
    // really needed, at least for new-DAE.  Derivatives out of the DCOP
    // are supposed to be zero at the first newton step.

    if (!(getSolverState().dcopFlag) && getSolverState().initTranFlag_ && getSolverState().newtonIter==0)
    {
      // re-set the state vector pointer that we are using to the "current"
      // pointer, rather than the "next" pointer.
      double * currStaVec = mi.extData.currStaVectorRawPtr;

      // intrinsic capacitors:
      currStaVec[mi.li_state_qb] = mi.qb;
      currStaVec[mi.li_state_qg] = mi.qg;
      currStaVec[mi.li_state_qd] = mi.qd;

      // parasitic capacitors:
      currStaVec[mi.li_state_qbs] = mi.qbs;
      currStaVec[mi.li_state_qbd] = mi.qbd;

      if( mi.nqsMod )
      {
        currStaVec[mi.li_state_qcheq] = mi.qcheq;
        currStaVec[mi.li_state_qcdump] = mi.qcdump;
      }
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/02/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, double * leadF, double * leadQ, double * junctionV)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & mi = *(*it);

    double * dFdxdVp = mi.extData.dFdxdVpVectorRawPtr;
    double * dQdxdVp = mi.extData.dQdxdVpVectorRawPtr;

    double coef(0.0);
    // F-vector:
    mi.cdreq_Jdxp = 0.0;
    mi.ceqbd_Jdxp = 0.0;
    mi.ceqbs_Jdxp = 0.0;

    // Do a few auxilliary calculations, derived from  3f5.
    // load current vector
    if (mi.mode >= 0)
    {
      mi.Gm = mi.gm;
      mi.Gmbs = mi.gmbs;
      mi.FwdSum = mi.Gm + mi.Gmbs;
      mi.RevSum = 0.0;

      mi.cdreq =  mi.model_.dtype * (mi.cd);
      mi.ceqbd = -mi.model_.dtype * (mi.csub);

      mi.ceqbs = 0.0;

      mi.gbbdp = -mi.gbds;
      mi.gbbsp = (mi.gbds + mi.gbgs + mi.gbbs);

      mi.gbdpg = mi.gbgs;
      mi.gbdpdp = mi.gbds;
      mi.gbdpb = mi.gbbs;
      mi.gbdpsp = -(mi.gbdpg + mi.gbdpdp + mi.gbdpb);

      mi.gbspg = 0.0;
      mi.gbspdp = 0.0;
      mi.gbspb = 0.0;
      mi.gbspsp = 0.0;
    }
    else
    {
      mi.Gm = -mi.gm;
      mi.Gmbs = -mi.gmbs;
      mi.FwdSum = 0.0;
      mi.RevSum = -(mi.Gm + mi.Gmbs);

      mi.cdreq = -mi.model_.dtype * (mi.cd);
      mi.ceqbs = -mi.model_.dtype * (mi.csub);

      mi.ceqbd = 0.0;

      mi.gbbsp = -mi.gbds;
      mi.gbbdp = (mi.gbds + mi.gbgs + mi.gbbs);

      mi.gbdpg = 0.0;
      mi.gbdpsp = 0.0;
      mi.gbdpb = 0.0;
      mi.gbdpdp = 0.0;

      mi.gbspg = mi.gbgs;
      mi.gbspsp = mi.gbds;
      mi.gbspb = mi.gbbs;
      mi.gbspdp = -(mi.gbspg + mi.gbspsp + mi.gbspb);
    }

    if (mi.model_.dtype > 0)
    {
      mi.ceqbs += (mi.cbs);
      mi.ceqbd += (mi.cbd);
    }
    else
    {
      mi.ceqbs -= (mi.cbs);
      mi.ceqbd -= (mi.cbd);
    }

    if (mi.drainConductance != 0.0)
    {
      fVec[mi.li_Drain] += mi.Idrain*mi.numberParallel;
    }
    if (mi.sourceConductance != 0.0)
    {
      fVec[mi.li_Source] += mi.Isource*mi.numberParallel;
    }

    fVec[mi.li_Bulk] += (mi.ceqbs + mi.ceqbd)*mi.numberParallel;
    fVec[mi.li_DrainPrime] += (-(mi.ceqbd - mi.cdreq)-mi.Idrain)*mi.numberParallel;
    fVec[mi.li_SourcePrime] += (-(mi.cdreq + mi.ceqbs)-mi.Isource)*mi.numberParallel;

    if( mi.loadLeadCurrent )
    {
      if (mi.drainConductance != 0.0)
      {
        leadF[mi.li_branch_dev_id] = mi.Idrain*mi.numberParallel;
      }
      else
      {
        leadF[mi.li_branch_dev_id] = (-(mi.ceqbd - mi.cdreq)-mi.Idrain)*mi.numberParallel;
      }
      if (mi.sourceConductance != 0.0)
      {
        leadF[mi.li_branch_dev_is] = mi.Isource*mi.numberParallel;
      }
      else
      {
        leadF[mi.li_branch_dev_is] = (-(mi.cdreq + mi.ceqbs)-mi.Isource)*mi.numberParallel;
      }
      leadF[mi.li_branch_dev_ig] = 0.0;
      leadF[mi.li_branch_dev_ib] = (mi.ceqbs + mi.ceqbd)*mi.numberParallel;

      junctionV[mi.li_branch_dev_id] = solVec[mi.li_Drain] - solVec[mi.li_Source];
      junctionV[mi.li_branch_dev_ig] = solVec[mi.li_Gate] - solVec[mi.li_Source];
      junctionV[mi.li_branch_dev_is] = 0.0;
      junctionV[mi.li_branch_dev_ib] = 0.0 ; 
    }

    // Initial condition support

    if( getSolverState().dcopFlag && mi.icVDSGiven )
    {
      coef = mi.extData.nextSolVectorRawPtr[mi.li_Ids];
      fVec[mi.li_Drain] += coef;
      fVec[mi.li_Source] += -coef;
      if( mi.loadLeadCurrent )
      {
        leadF[mi.li_branch_dev_id]= coef;
        leadF[mi.li_branch_dev_is]= -coef;
      }
    }

    if( getSolverState().dcopFlag && mi.icVGSGiven )
    {
      coef = mi.extData.nextSolVectorRawPtr[mi.li_Igs];
      fVec[mi.li_Gate] += coef;
      fVec[mi.li_Source] += -coef;
      if( mi.loadLeadCurrent )
      {
        leadF[mi.li_branch_dev_ig]= coef;
        leadF[mi.li_branch_dev_is]= -coef;
      }
    }


    if( getSolverState().dcopFlag && mi.icVBSGiven )
    {
      coef = mi.extData.nextSolVectorRawPtr[mi.li_Ibs];
      fVec[mi.li_Bulk] += coef;
      fVec[mi.li_Source] += -coef;
      if( mi.loadLeadCurrent )
      {
        leadF[mi.li_branch_dev_ib]= coef;
        leadF[mi.li_branch_dev_is]= -coef;
      }
    }

    //////////////////////////////////////////////////
    // limiting section:
    if (getDeviceOptions().voltageLimiterFlag)
    {
      if (!mi.origFlag)
      {
        if (mi.mode >= 0)
        {
          // option 1
          double tmp = mi.model_.dtype * (-mi.gds * (mi.vds-mi.vds_orig) -
                                          mi.Gm * (mi.vgs-mi.vgs_orig) -
                                          mi.Gmbs * (mi.vbs-mi.vbs_orig));

          mi.cdreq_Jdxp += tmp;
          mi.cdreq      += tmp;

          tmp = -mi.model_.dtype * (-mi.gbds * (mi.vds-mi.vds_orig) -
                                    mi.gbgs * (mi.vgs-mi.vgs_orig) -
                                    mi.gbbs * (mi.vbs-mi.vbs_orig));
          mi.ceqbd_Jdxp += tmp;
          mi.ceqbd      += tmp;
        }
        else
        {
          // option 2
          double tmp = -mi.model_.dtype * (mi.gds * (mi.vds-mi.vds_orig) +
                                           mi.Gm * (mi.vgd-mi.vgd_orig) +
                                           mi.Gmbs * (mi.vbd-mi.vbd_orig));
          mi.cdreq_Jdxp += tmp;
          mi.cdreq      += tmp;

          tmp = -mi.model_.dtype * (mi.gbds * (mi.vds-mi.vds_orig) -
                                    mi.gbgs * (mi.vgd-mi.vgd_orig) -
                                    mi.gbbs * (mi.vbd-mi.vbd_orig));
          mi.ceqbd_Jdxp += tmp;
          mi.ceqbd      += tmp;
        }


        if (mi.model_.dtype > 0)
        {
          mi.ceqbs_Jdxp += (-mi.gbs*(mi.vbs-mi.vbs_orig));
          mi.ceqbs      += (-mi.gbs*(mi.vbs-mi.vbs_orig));

          mi.ceqbd_Jdxp += (-mi.gbd*(mi.vbd-mi.vbd_orig));
          mi.ceqbd      += (-mi.gbd*(mi.vbd-mi.vbd_orig));
        }
        else
        {
          mi.ceqbs_Jdxp -= (-mi.gbs*(mi.vbs-mi.vbs_orig) );
          mi.ceqbs      -= (-mi.gbs*(mi.vbs-mi.vbs_orig) );

          mi.ceqbd_Jdxp -= (-mi.gbd*(mi.vbd-mi.vbd_orig) );
          mi.ceqbd      -= (-mi.gbd*(mi.vbd-mi.vbd_orig) );
        }

        dFdxdVp[mi.li_Bulk] += -(mi.ceqbs_Jdxp+mi.ceqbd_Jdxp)*mi.numberParallel;
        dFdxdVp[mi.li_DrainPrime] += (mi.ceqbd_Jdxp-mi.cdreq_Jdxp)*mi.numberParallel;
        dFdxdVp[mi.li_SourcePrime] += (mi.cdreq_Jdxp+mi.ceqbs_Jdxp)*mi.numberParallel;

      } // orig flag.
    } // voltage limiter flag

    // Row associated with icVBS
    if( getSolverState().dcopFlag && mi.icVBSGiven )
    {
      // get the voltage drop from the previous solution
      double cVs = mi.extData.nextSolVectorRawPtr[mi.li_Source];
      double cVb = mi.extData.nextSolVectorRawPtr[mi.li_Bulk];

      fVec[mi.li_Ibs] += (cVb - cVs - mi.icVBS);
    }

    // Row associated with icVDS
    if( getSolverState().dcopFlag && mi.icVDSGiven )
    {
      // get the voltage drop from the previous solution
      double cVd = mi.extData.nextSolVectorRawPtr[mi.li_Drain];
      double cVs = mi.extData.nextSolVectorRawPtr[mi.li_Source];

      fVec[mi.li_Ids] += (cVd - cVs - mi.icVDS);
    }

    // Row associated with icVGS
    if( getSolverState().dcopFlag && mi.icVGSGiven )
    {
      // get the voltage drop from the previous solution
      double cVg = mi.extData.nextSolVectorRawPtr[mi.li_Gate];
      double cVs = mi.extData.nextSolVectorRawPtr[mi.li_Source];

      fVec[mi.li_Igs] += (cVg - cVs - mi.icVGS);
    }

    // Q-vector:

    mi.auxChargeCalculations ();

    double Qeqqg = 0.0;   // gate charge
    double Qeqqb = 0.0;   // bulk charge
    double Qeqqd = 0.0;   // drain charge
    double Qqdef = 0.0;   // nqs-related charge.
    double Qqcheq = 0.0;  // nqs-related charge.

    // These 3 vars are class variables, and are set up elsewhere.
    //double Qeqqg_Jdxp = 0.0; // limiter, related to gate  cap.
    //double Qeqqb_Jdxp = 0.0; // limiter, related to bulk  cap.
    //double Qeqqd_Jdxp = 0.0; // limiter, related to drain cap.

    if (mi.model_.dtype > 0)
    {
      Qeqqg  = mi.qg;
      Qeqqb  = mi.qb;
      Qeqqd  = mi.qd;
      Qqdef  = mi.qcdump; // this needs to be fixed...
      Qqcheq = mi.qcheq;
    }
    else  // need to convert these to charges.
    {
      Qeqqg  = -mi.qg;
      Qeqqb  = -mi.qb;
      Qeqqd  = -mi.qd;
      Qqdef  = -mi.qcdump;
      Qqcheq = -mi.qcheq;
    }

    qVec[mi.li_Gate] += Qeqqg*mi.numberParallel;
    qVec[mi.li_Bulk] += (Qeqqb)*mi.numberParallel;
    qVec[mi.li_DrainPrime] += (-(-Qeqqd))*mi.numberParallel;
    qVec[mi.li_SourcePrime] += (-(+ Qeqqg + Qeqqb + Qeqqd))*mi.numberParallel;

    if (mi.nqsMod)
    {
      // 7 equ. for nqs modification. charge equation.
      qVec[mi.li_Charge] += -(Qqcheq - Qqdef)*mi.numberParallel;
    }

    if( mi.loadLeadCurrent )
    {
      if (mi.drainConductance == 0.0)
      {
        leadQ[mi.li_branch_dev_id] = (-(-Qeqqd))*mi.numberParallel;
      }
      if (mi.sourceConductance == 0.0)
      {
        leadQ[mi.li_branch_dev_is] = (-(Qeqqg + Qeqqb + Qeqqd))*mi.numberParallel;
      }
      leadQ[mi.li_branch_dev_ig] = Qeqqg*mi.numberParallel;
      leadQ[mi.li_branch_dev_ib] = (Qeqqb)*mi.numberParallel;
    }

    //////////////////////////////////////////////////
    // limiting section:
    if (getDeviceOptions().voltageLimiterFlag)
    {
      // Need the following:
      //  Qeqqg_Jdxp
      //  Qeqqb_Jdxp
      //  Qeqqd_Jdxp
      if (mi.model_.dtype < 0)
      {
        mi.Qeqqg_Jdxp = -mi.Qeqqg_Jdxp;
        mi.Qeqqb_Jdxp = -mi.Qeqqb_Jdxp;
        mi.Qeqqd_Jdxp = -mi.Qeqqd_Jdxp;
      }

      if (!mi.origFlag)
      {
        dQdxdVp[mi.li_Gate] += -mi.Qeqqg_Jdxp*mi.numberParallel;
        dQdxdVp[mi.li_Bulk] += -(+mi.Qeqqb_Jdxp)*mi.numberParallel;
        dQdxdVp[mi.li_DrainPrime] += (-mi.Qeqqd_Jdxp)*mi.numberParallel;
        dQdxdVp[mi.li_SourcePrime] += (+mi.Qeqqg_Jdxp+mi.Qeqqb_Jdxp+mi.Qeqqd_Jdxp) *mi.numberParallel;
      } // orig flag.
    } // limiter flag
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
// Creation Date : 12/02/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & mi = *(*it);

    // F-matrix:
    // Row corresponding to the KCL for the drain node:

    *mi.f_DrainEquDrainNodePtr
      += mi.drainConductance*mi.numberParallel;

    *mi.f_DrainEquDrainPrimeNodePtr
      -= mi.drainConductance*mi.numberParallel;

    // Extra term for initial conditions on Vds in operating point
    if( getSolverState().dcopFlag && mi.icVDSGiven )
    {
      *mi.f_DrainEquIdsPtr += 1.0;
    }

    // Row corresponding to the KCL for the source node:

    *mi.f_SourceEquSourceNodePtr
      += mi.sourceConductance*mi.numberParallel;

    *mi.f_SourceEquSourcePrimeNodePtr
      -= mi.sourceConductance*mi.numberParallel;

    // Extra term for initial conditions on Vbs in operating point
    if( getSolverState().dcopFlag && mi.icVBSGiven )
    {
      *mi.f_SourceEquIbsPtr -= 1.0;
    }
    // Extra term for initial conditions on Vds in operating point
    if( getSolverState().dcopFlag && mi.icVDSGiven )
    {
      *mi.f_SourceEquIdsPtr -= 1.0;
    }
    // Extra term for initial conditions on Vgs in operating point
    if( getSolverState().dcopFlag && mi.icVGSGiven )
    {
      *mi.f_SourceEquIgsPtr -= 1.0;
    }

    // Row corresponding to the KCL for the gate node: NOTHING
    // Check this later.   ERK.
    //
    //  The terms beginning with "gc" (gcggb, etc.) definately do NOT
    //  belong here.  I'm not sure aboug the gg terms.  On one hand, the
    //  rhs vector component for the gate node ONLY seems to take
    //  capacitive currents, which implies that all of these are capacitive
    //  conductances.  On the other hand, the gg terms do not appear to
    //  have been created by multiplying by ag0 = pdt = 1/dt.  Generally
    //  capacitive conductances are of the form g = C/dt, and the gg terms
    //  do not have this form.
    //
    //  For now, the gg issue is moot b/c those terms are only nonzero
    //  if mi.nqsMod = 1, which is not a supported option.

    //  However, the gg
    //  terms (mi.ggtg, mi.ggtb, ggtd and ggts)
    //
    //(*JMatPtr)[GateEquGateNodePtr
    //  += (gcggb - mi.ggtg)*mi.numberParallel;
    //(*JMatPtr)[GateEquBulkNodePtr
    //  -= (gcggb + gcgdb + gcgsb + mi.ggtb)*mi.numberParallel;
    //(*JMatPtr)[GateEquDrainPrimeNodePtr
    //  += (gcgdb - ggtd)*mi.numberParallel;
    //(*JMatPtr)[GateEquSourcePrimeNodePtr
    //  += (gcgsb - ggts)*mi.numberParallel;

    // initial conditions on gate node
    // Extra term for initial conditions on Vgs in operating point
    if( getSolverState().dcopFlag && mi.icVGSGiven )
    {
      *mi.f_GateEquIgsPtr += 1.0;
    }

    // Row corresponding to the KCL for the bulk node:

    *mi.f_BulkEquGateNodePtr
      += (- mi.gbgs)*mi.numberParallel;


    *mi.f_BulkEquBulkNodePtr
      += (mi.gbd + mi.gbs - mi.gbbs)*mi.numberParallel;


    *mi.f_BulkEquDrainPrimeNodePtr
      += (- mi.gbd + mi.gbbdp)*mi.numberParallel;


    *mi.f_BulkEquSourcePrimeNodePtr
      += (- mi.gbs + mi.gbbsp)*mi.numberParallel;

    // Extra term for initial conditions on Vbs in operating point
    if( getSolverState().dcopFlag && mi.icVBSGiven )
    {
      *mi.f_BulkEquIbsPtr += 1.0;
    }

    // Row corresponding to the KCL for the drain prime node:

    *mi.f_DrainPrimeEquDrainNodePtr
      -= mi.drainConductance*mi.numberParallel;

    *mi.f_DrainPrimeEquBulkNodePtr
      -= (mi.gbd - mi.Gmbs - mi.dxpart*mi.ggtb
          - mi.T1global*mi.ddxpart_dVb - mi.gbdpb)*mi.numberParallel;

    *mi.f_DrainPrimeEquGateNodePtr
      += (mi.Gm + mi.dxpart*mi.ggtg + mi.T1global*mi.ddxpart_dVg + mi.gbdpg)
      *mi.numberParallel;

    *mi.f_DrainPrimeEquDrainPrimeNodePtr
      += (mi.drainConductance + mi.gds + mi.gbd + mi.RevSum + mi.dxpart*mi.ggtd
          + mi.T1global*mi.ddxpart_dVd + mi.gbdpdp)*mi.numberParallel;

    *mi.f_DrainPrimeEquSourcePrimeNodePtr
      -= (mi.gds + mi.FwdSum - mi.dxpart*mi.ggts - mi.T1global*mi.ddxpart_dVs - mi.gbdpsp)
      *mi.numberParallel;

    // Row corresponding to the KCL for the source prime node:
    *mi.f_SourcePrimeEquGateNodePtr
      += (- mi.Gm + mi.sxpart*mi.ggtg + mi.T1global*mi.dsxpart_dVg + mi.gbspg)
      *mi.numberParallel;

    *mi.f_SourcePrimeEquBulkNodePtr
      -= (mi.gbs + mi.Gmbs - mi.sxpart*mi.ggtb
          - mi.T1global*mi.dsxpart_dVb - mi.gbspb)*mi.numberParallel;

    *mi.f_SourcePrimeEquSourceNodePtr
      -= mi.sourceConductance*mi.numberParallel;

    *mi.f_SourcePrimeEquDrainPrimeNodePtr
      -= (mi.gds + mi.RevSum - mi.sxpart*mi.ggtd - mi.T1global*mi.dsxpart_dVd - mi.gbspdp)
      *mi.numberParallel;

    *mi.f_SourcePrimeEquSourcePrimeNodePtr
      += (mi.sourceConductance + mi.gds + mi.gbs + mi.FwdSum + mi.sxpart*mi.ggts
          + mi.T1global*mi.dsxpart_dVs + mi.gbspsp)*mi.numberParallel;

    // Row associated with the charge equation
    if (mi.nqsMod)
    {
      DevelFatal(*this) << " nqsMod=1 is not ready yet.  Re-run with nqsMod=0.";
    }

    // Initial condition rows
    // Row associated with mi.icVBS
    if( mi.icVBSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        *mi.f_icVBSEquVbPtr += 1.0;
        *mi.f_icVBSEquVsPtr -= 1.0;
      }
      else
      {
        *mi.f_icVBSEquIbsPtr += 1.0;
      }
    }

    // Row associated with mi.icVDS
    if( mi.icVDSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        *mi.f_icVDSEquVdPtr += 1.0;
        *mi.f_icVDSEquVsPtr -= 1.0;
      }
      else
      {
        *mi.f_icVDSEquIdsPtr += 1.0;
      }
    }

    // Row associated with mi.icVGS
    if( mi.icVGSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        *mi.f_icVGSEquVgPtr += 1.0;
        *mi.f_icVGSEquVsPtr -= 1.0;
      }
      else
      {
        *mi.f_icVGSEquIgsPtr += 1.0;
      }
    }


    {
      // Row corresponding to the KCL for the drain node: NOTHING

      // Row corresponding to the KCL for the source node: NOTHING

      // Row corresponding to the KCL for the gate node:
      // Check this later.   ERK.  See the comments in the function
      // loadDAE*mi.q_, regarding ggtg, ggtb, etc.
      //
      // For now I am leaving out the gg terms, as they are zero when
      // nqsMod=0, which is always true.
      //
      *mi.q_GateEquGateNodePtr
        += (mi.CAPcggb )*mi.numberParallel;
      *mi.q_GateEquBulkNodePtr
        -= (mi.CAPcggb + mi.CAPcgdb + mi.CAPcgsb )*mi.numberParallel;
      *mi.q_GateEquDrainPrimeNodePtr
        += (mi.CAPcgdb )*mi.numberParallel;
      *mi.q_GateEquSourcePrimeNodePtr
        += (mi.CAPcgsb )*mi.numberParallel;

      // Row corresponding to the KCL for the bulk node:
      *mi.q_BulkEquGateNodePtr
        += (mi.CAPcbgb)*mi.numberParallel;

      *mi.q_BulkEquBulkNodePtr
        += (- mi.CAPcbgb - mi.CAPcbdb - mi.CAPcbsb)*mi.numberParallel;

      *mi.q_BulkEquDrainPrimeNodePtr
        += (mi.CAPcbdb)*mi.numberParallel;

      *mi.q_BulkEquSourcePrimeNodePtr
        += (mi.CAPcbsb)*mi.numberParallel;


      // Row corresponding to the KCL for the drain prime node:
      *mi.q_DrainPrimeEquBulkNodePtr
        -= (+ mi.CAPcdgb + mi.CAPcddb + mi.CAPcdsb )*mi.numberParallel;

      *mi.q_DrainPrimeEquGateNodePtr
        += (mi.CAPcdgb) *mi.numberParallel;

      *mi.q_DrainPrimeEquDrainPrimeNodePtr
        += (+ mi.CAPcddb )*mi.numberParallel;

      *mi.q_DrainPrimeEquSourcePrimeNodePtr
        -= (- mi.CAPcdsb) *mi.numberParallel;

      // Row corresponding to the KCL for the source prime node:
      *mi.q_SourcePrimeEquGateNodePtr
        += (mi.CAPcsgb) *mi.numberParallel;

      *mi.q_SourcePrimeEquBulkNodePtr
        -= (+ mi.CAPcsgb + mi.CAPcsdb + mi.CAPcssb) *mi.numberParallel;

      *mi.q_SourcePrimeEquDrainPrimeNodePtr
        -= (- mi.CAPcsdb) *mi.numberParallel;

      *mi.q_SourcePrimeEquSourcePrimeNodePtr
        += (+ mi.CAPcssb) *mi.numberParallel;

      // Row associated with the charge equation
      // This is currently not supported.
      if (mi.nqsMod)
      {
        DevelFatal(*this).in("Master::loadDAEMatrices")
          << " nqsMod=1 is not ready yet.  Re-run with nqsMod=0\n";
      }
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
// Creation Date : 12/02/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & mi = *(*it);

    int count = 0;
    // F-matrix:
    // Row corresponding to the KCL for the drain node:

    dFdx [mi.li_Drain][mi.ADrainEquDrainNodeOffset]
      += mi.drainConductance*mi.numberParallel;

    dFdx [mi.li_Drain][mi.ADrainEquDrainPrimeNodeOffset]
      -= mi.drainConductance*mi.numberParallel;

    // Extra term for initial conditions on Vds in operating point
    if( getSolverState().dcopFlag && mi.icVDSGiven )
    {
      dFdx[mi.li_Drain][mi.ADrainEquIdsOffset] += 1.0;
    }

    // Row corresponding to the KCL for the source node:

    dFdx [mi.li_Source][mi.ASourceEquSourceNodeOffset]
      += mi.sourceConductance*mi.numberParallel;

    dFdx [mi.li_Source][mi.ASourceEquSourcePrimeNodeOffset]
      -= mi.sourceConductance*mi.numberParallel;

    // Extra term for initial conditions on Vbs in operating point
    if( getSolverState().dcopFlag && mi.icVBSGiven )
    {
      dFdx[mi.li_Source][mi.ASourceEquIbsOffset] -= 1.0;
    }
    // Extra term for initial conditions on Vds in operating point
    if( getSolverState().dcopFlag && mi.icVDSGiven )
    {
      dFdx[mi.li_Source][mi.ASourceEquIdsOffset] -= 1.0;
    }
    // Extra term for initial conditions on Vgs in operating point
    if( getSolverState().dcopFlag && mi.icVGSGiven )
    {
      dFdx[mi.li_Source][mi.ASourceEquIgsOffset] -= 1.0;
    }

    // Row corresponding to the KCL for the gate node: NOTHING
    // Check this later.   ERK.
    //
    //  The terms beginning with "gc" (gcggb, etc.) definately do NOT
    //  belong here.  I'm not sure aboug the gg terms.  On one hand, the
    //  rhs vector component for the gate node ONLY seems to take
    //  capacitive currents, which implies that all of these are capacitive
    //  conductances.  On the other hand, the gg terms do not appear to
    //  have been created by multiplying by ag0 = pdt = 1/dt.  Generally
    //  capacitive conductances are of the form g = C/dt, and the gg terms
    //  do not have this form.
    //
    //  For now, the gg issue is moot b/c those terms are only nonzero
    //  if mi.nqsMod = 1, which is not a supported option.

    //  However, the gg
    //  terms (mi.ggtg, mi.ggtb, ggtd and ggts)
    //
    //(*JMatPtr)[mi.li_Gate][mi.AGateEquGateNodeOffset]
    //  += (gcggb - mi.ggtg)*mi.numberParallel;
    //(*JMatPtr)[mi.li_Gate][mi.AGateEquBulkNodeOffset]
    //  -= (gcggb + gcgdb + gcgsb + mi.ggtb)*mi.numberParallel;
    //(*JMatPtr)[mi.li_Gate][mi.AGateEquDrainPrimeNodeOffset]
    //  += (gcgdb - ggtd)*mi.numberParallel;
    //(*JMatPtr)[mi.li_Gate][mi.AGateEquSourcePrimeNodeOffset]
    //  += (gcgsb - ggts)*mi.numberParallel;

    // initial conditions on gate node
    // Extra term for initial conditions on Vgs in operating point
    if( getSolverState().dcopFlag && mi.icVGSGiven )
    {
      dFdx[mi.li_Gate][mi.AGateEquIgsOffset] += 1.0;
    }

    // Row corresponding to the KCL for the bulk node:

    dFdx [mi.li_Bulk][mi.ABulkEquGateNodeOffset]
      += (- mi.gbgs)*mi.numberParallel;


    dFdx [mi.li_Bulk][mi.ABulkEquBulkNodeOffset]
      += (mi.gbd + mi.gbs - mi.gbbs)*mi.numberParallel;


    dFdx [mi.li_Bulk][mi.ABulkEquDrainPrimeNodeOffset]
      += (- mi.gbd + mi.gbbdp)*mi.numberParallel;


    dFdx [mi.li_Bulk][mi.ABulkEquSourcePrimeNodeOffset]
      += (- mi.gbs + mi.gbbsp)*mi.numberParallel;

    // Extra term for initial conditions on Vbs in operating point
    if( getSolverState().dcopFlag && mi.icVBSGiven )
    {
      dFdx[mi.li_Bulk][mi.ABulkEquIbsOffset] += 1.0;
    }

    // Row corresponding to the KCL for the drain prime node:

    dFdx [mi.li_DrainPrime][mi.ADrainPrimeEquDrainNodeOffset]
      -= mi.drainConductance*mi.numberParallel;


    dFdx [mi.li_DrainPrime][mi.ADrainPrimeEquBulkNodeOffset]
      -= (mi.gbd - mi.Gmbs - mi.dxpart*mi.ggtb
          - mi.T1global*mi.ddxpart_dVb - mi.gbdpb)*mi.numberParallel;


    dFdx [mi.li_DrainPrime][mi.ADrainPrimeEquGateNodeOffset]
      += (mi.Gm + mi.dxpart*mi.ggtg + mi.T1global*mi.ddxpart_dVg + mi.gbdpg)
      *mi.numberParallel;


    dFdx [mi.li_DrainPrime][mi.ADrainPrimeEquDrainPrimeNodeOffset]
      += (mi.drainConductance + mi.gds + mi.gbd + mi.RevSum + mi.dxpart*mi.ggtd
          + mi.T1global*mi.ddxpart_dVd + mi.gbdpdp)*mi.numberParallel;


    dFdx [mi.li_DrainPrime][mi.ADrainPrimeEquSourcePrimeNodeOffset]
      -= (mi.gds + mi.FwdSum - mi.dxpart*mi.ggts - mi.T1global*mi.ddxpart_dVs - mi.gbdpsp)
      *mi.numberParallel;

    // Row corresponding to the KCL for the source prime node:

    dFdx [mi.li_SourcePrime][mi.ASourcePrimeEquGateNodeOffset]
      += (- mi.Gm + mi.sxpart*mi.ggtg + mi.T1global*mi.dsxpart_dVg + mi.gbspg)
      *mi.numberParallel;


    dFdx [mi.li_SourcePrime][mi.ASourcePrimeEquBulkNodeOffset]
      -= (mi.gbs + mi.Gmbs - mi.sxpart*mi.ggtb
          - mi.T1global*mi.dsxpart_dVb - mi.gbspb)*mi.numberParallel;


    dFdx [mi.li_SourcePrime][mi.ASourcePrimeEquSourceNodeOffset]
      -= mi.sourceConductance*mi.numberParallel;


    dFdx [mi.li_SourcePrime][mi.ASourcePrimeEquDrainPrimeNodeOffset]
      -= (mi.gds + mi.RevSum - mi.sxpart*mi.ggtd - mi.T1global*mi.dsxpart_dVd - mi.gbspdp)
      *mi.numberParallel;


    dFdx [mi.li_SourcePrime][mi.ASourcePrimeEquSourcePrimeNodeOffset]
      += (mi.sourceConductance + mi.gds + mi.gbs + mi.FwdSum + mi.sxpart*mi.ggts
          + mi.T1global*mi.dsxpart_dVs + mi.gbspsp)*mi.numberParallel;

    // Row associated with the charge equation
    if (mi.nqsMod)
    {
      DevelFatal(*this).in("Master::loadDAEMatrices")
        << " nqsMod=1 is not ready yet.  Re-run with nqsMod=0";
    }

    // Initial condition rows
    // Row associated with mi.icVBS
    if( mi.icVBSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        dFdx[mi.li_Ibs][mi.icVBSEquVbOffset] += 1.0;
        dFdx[mi.li_Ibs][mi.icVBSEquVsOffset] -= 1.0;
      }
      else
      {
        dFdx[mi.li_Ibs][mi.icVBSEquIbsOffset] += 1.0;
      }
    }

    // Row associated with mi.icVDS
    if( mi.icVDSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        dFdx[mi.li_Ids][mi.icVDSEquVdOffset] += 1.0;
        dFdx[mi.li_Ids][mi.icVDSEquVsOffset] -= 1.0;
      }
      else
      {
        dFdx[mi.li_Ids][mi.icVDSEquIdsOffset] += 1.0;
      }
    }

    // Row associated with mi.icVGS
    if( mi.icVGSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        dFdx[mi.li_Igs][mi.icVGSEquVgOffset] += 1.0;
        dFdx[mi.li_Igs][mi.icVGSEquVsOffset] -= 1.0;
      }
      else
      {
        dFdx[mi.li_Igs][mi.icVGSEquIgsOffset] += 1.0;
      }
    }

    {
      // Row corresponding to the KCL for the drain node: NOTHING

      // Row corresponding to the KCL for the source node: NOTHING

      // Row corresponding to the KCL for the gate node:
      // Check this later.   ERK.  See the comments in the function
      // loadDAEdQdx, regarding ggtg, ggtb, etc.
      //
      // For now I am leaving out the gg terms, as they are zero when
      // nqsMod=0, which is always true.
      //
      dQdx[mi.li_Gate][mi.AGateEquGateNodeOffset]
        += (mi.CAPcggb )*mi.numberParallel;
      dQdx[mi.li_Gate][mi.AGateEquBulkNodeOffset]
        -= (mi.CAPcggb + mi.CAPcgdb + mi.CAPcgsb )*mi.numberParallel;
      dQdx[mi.li_Gate][mi.AGateEquDrainPrimeNodeOffset]
        += (mi.CAPcgdb )*mi.numberParallel;
      dQdx[mi.li_Gate][mi.AGateEquSourcePrimeNodeOffset]
        += (mi.CAPcgsb )*mi.numberParallel;

      // Row corresponding to the KCL for the bulk node:
      dQdx[mi.li_Bulk][mi.ABulkEquGateNodeOffset]
        += (mi.CAPcbgb)*mi.numberParallel;

      dQdx[mi.li_Bulk][mi.ABulkEquBulkNodeOffset]
        += (- mi.CAPcbgb - mi.CAPcbdb - mi.CAPcbsb)*mi.numberParallel;

      dQdx[mi.li_Bulk][mi.ABulkEquDrainPrimeNodeOffset]
        += (mi.CAPcbdb)*mi.numberParallel;

      dQdx[mi.li_Bulk][mi.ABulkEquSourcePrimeNodeOffset]
        += (mi.CAPcbsb)*mi.numberParallel;


      // Row corresponding to the KCL for the drain prime node:
      dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquBulkNodeOffset]
        -= (+ mi.CAPcdgb + mi.CAPcddb + mi.CAPcdsb )*mi.numberParallel;

      dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquGateNodeOffset]
        += (mi.CAPcdgb) *mi.numberParallel;

      dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquDrainPrimeNodeOffset]
        += (+ mi.CAPcddb )*mi.numberParallel;

      dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquSourcePrimeNodeOffset]
        -= (- mi.CAPcdsb) *mi.numberParallel;

      // Row corresponding to the KCL for the source prime node:
      dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquGateNodeOffset]
        += (mi.CAPcsgb) *mi.numberParallel;

      dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquBulkNodeOffset]
        -= (+ mi.CAPcsgb + mi.CAPcsdb + mi.CAPcssb) *mi.numberParallel;

      dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquDrainPrimeNodeOffset]
        -= (- mi.CAPcsdb) *mi.numberParallel;

      dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquSourcePrimeNodeOffset]
        += (+ mi.CAPcssb) *mi.numberParallel;

      // Row associated with the charge equation
      // This is currently not supported.
      if (mi.nqsMod)
      {
        DevelFatal(*this).in("Master::loadDAEMatrices")
          << " nqsMod=1 is not ready yet.  Re-run with nqsMod=0";
      }
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
      ((deviceMap.find("M")!=deviceMap.end()) 
      && ((levelSet.find(9)!=levelSet.end()) || (levelSet.find(49)!=levelSet.end()))))
  {
    MOSFET1::registerDevice();

    Config<Traits>::addConfiguration()
      .registerDevice("m", 9)
      .registerDevice("m", 49)
      .registerModelType("pmos", 9)
      .registerModelType("nmos", 9)
      .registerModelType("pmos", 49)
      .registerModelType("nmos", 49);
  }
}


//-----------------------------------------------------------------------------
// sensitivity related code.


//-----------------------------------------------------------------------------
template <typename ScalarT> 
bool updateTemperature (
  const ScalarT & temp_tmp,
  const ScalarT & dtemp,
  const ScalarT Eg0,
  const ScalarT ni,
  const ScalarT Vtm0,
  const ScalarT model_tnom,
  const ScalarT model_jctTempExponent,
  const ScalarT model_jctEmissionCoeff,
  const ScalarT model_jctSatCurDensity,
  const ScalarT model_jctSidewallSatCurDensity,
  const ScalarT model_tcj,
  const ScalarT model_unitAreaJctCap,
  const ScalarT model_tcjsw,
  const ScalarT model_unitLengthSidewallJctCap, 
  const ScalarT model_tcjswg,
  const ScalarT model_unitLengthGateSidewallJctCap,
  const ScalarT model_bulkJctPotential,
  const ScalarT model_tpb,
  ScalarT PhiBSWTemp,
  const ScalarT model_sidewallJctPotential,
  const ScalarT model_tpbsw,
  ScalarT PhiBSWGTemp,
  const ScalarT model_GatesidewallJctPotential,
  const ScalarT model_tpbswg,
  const ScalarT l,
  const ScalarT w,
  const ScalarT model_Lln,
  const ScalarT model_Lwn,
  const ScalarT model_Ll,
  const ScalarT model_Lw,
  const ScalarT model_Lwl,
  const ScalarT model_Lint,
  const ScalarT model_Llc,
  const ScalarT model_Lwc,
  const ScalarT model_Lwlc,
  const ScalarT model_dlc,
  const ScalarT model_Wln,
  const ScalarT model_Wwn,
  const ScalarT model_Wl,
  const ScalarT model_Ww,
  const ScalarT model_Wwl,
  const ScalarT model_Wint,
  const ScalarT model_Wlc,
  const ScalarT model_Wwc,
  const ScalarT model_Wwlc,
  const ScalarT model_dwc,
  const ScalarT model_cgdo,
  const ScalarT model_cgso,
  const ScalarT model_cgbo,
  const ScalarT model_cox,
  const bool model_npeakGiven,
  const bool model_gamma1Given,
  const ScalarT model_ijth,
  const ScalarT sourceArea, 
  const ScalarT sourcePerimeter,
  const ScalarT drainArea,
  const ScalarT drainPerimeter,

  const ScalarT model_cdsc,
  const ScalarT model_lcdsc,
  const ScalarT model_wcdsc,
  const ScalarT model_pcdsc,
  const ScalarT model_cdscb,
  const ScalarT model_lcdscb,
  const ScalarT model_wcdscb,
  const ScalarT model_pcdscb,
  const ScalarT model_cdscd,
  const ScalarT model_lcdscd,
  const ScalarT model_wcdscd,
  const ScalarT model_pcdscd,
  const ScalarT model_cit,
  const ScalarT model_lcit,
  const ScalarT model_wcit,
  const ScalarT model_pcit,
  const ScalarT model_nfactor,
  const ScalarT model_lnfactor,
  const ScalarT model_wnfactor,
  const ScalarT model_pnfactor,
  const ScalarT model_xj,
  const ScalarT model_lxj,
  const ScalarT model_wxj,
  const ScalarT model_pxj,
  const ScalarT model_vsat,
  const ScalarT model_lvsat,
  const ScalarT model_wvsat,
  const ScalarT model_pvsat,
  const ScalarT model_at,
  const ScalarT model_lat,
  const ScalarT model_wat,
  const ScalarT model_pat,
  const ScalarT model_a0,
  const ScalarT model_la0,
  const ScalarT model_wa0,
  const ScalarT model_pa0,
  const ScalarT model_ags,
  const ScalarT model_lags,
  const ScalarT model_wags,
  const ScalarT model_pags,
  const ScalarT model_a1,
  const ScalarT model_la1,
  const ScalarT model_wa1,
  const ScalarT model_pa1,
  const ScalarT model_a2,
  const ScalarT model_la2,
  const ScalarT model_wa2,
  const ScalarT model_pa2,
  const ScalarT model_keta,
  const ScalarT model_lketa,
  const ScalarT model_wketa,
  const ScalarT model_pketa,
  const ScalarT model_nsub,
  const ScalarT model_lnsub,
  const ScalarT model_wnsub,
  const ScalarT model_pnsub,
  const ScalarT model_npeak,
  const ScalarT model_lnpeak,
  const ScalarT model_wnpeak,
  const ScalarT model_pnpeak,
  const ScalarT model_ngate,
  const ScalarT model_lngate,
  const ScalarT model_wngate,
  const ScalarT model_pngate,
  const ScalarT model_gamma1,
  const ScalarT model_lgamma1,
  const ScalarT model_wgamma1,
  const ScalarT model_pgamma1,
  const ScalarT model_gamma2,
  const ScalarT model_lgamma2,
  const ScalarT model_wgamma2,
  const ScalarT model_pgamma2,
  const ScalarT model_vbx,
  const ScalarT model_lvbx,
  const ScalarT model_wvbx,
  const ScalarT model_pvbx,
  const ScalarT model_vbm,
  const ScalarT model_lvbm,
  const ScalarT model_wvbm,
  const ScalarT model_pvbm,
  const ScalarT model_xt,
  const ScalarT model_lxt,
  const ScalarT model_wxt,
  const ScalarT model_pxt,
  const ScalarT model_vfb,
  const ScalarT model_lvfb,
  const ScalarT model_wvfb,
  const ScalarT model_pvfb,
  const ScalarT model_k1,
  const ScalarT model_lk1,
  const ScalarT model_wk1,
  const ScalarT model_pk1,
  const ScalarT model_kt1,
  const ScalarT model_lkt1,
  const ScalarT model_wkt1,
  const ScalarT model_pkt1,
  const ScalarT model_kt1l,
  const ScalarT model_lkt1l,
  const ScalarT model_wkt1l,
  const ScalarT model_pkt1l,
  const ScalarT model_k2,
  const ScalarT model_lk2,
  const ScalarT model_wk2,
  const ScalarT model_pk2,
  const ScalarT model_kt2,
  const ScalarT model_lkt2,
  const ScalarT model_wkt2,
  const ScalarT model_pkt2,
  const ScalarT model_k3,
  const ScalarT model_lk3,
  const ScalarT model_wk3,
  const ScalarT model_pk3,
  const ScalarT model_k3b,
  const ScalarT model_lk3b,
  const ScalarT model_wk3b,
  const ScalarT model_pk3b,
  const ScalarT model_w0,
  const ScalarT model_lw0,
  const ScalarT model_ww0,
  const ScalarT model_pw0,
  const ScalarT model_nlx,
  const ScalarT model_lnlx,
  const ScalarT model_wnlx,
  const ScalarT model_pnlx,
  const ScalarT model_dvt0,
  const ScalarT model_ldvt0,
  const ScalarT model_wdvt0,
  const ScalarT model_pdvt0,
  const ScalarT model_dvt1,
  const ScalarT model_ldvt1,
  const ScalarT model_wdvt1,
  const ScalarT model_pdvt1,
  const ScalarT model_dvt2,
  const ScalarT model_ldvt2,
  const ScalarT model_wdvt2,
  const ScalarT model_pdvt2,
  const ScalarT model_dvt0w,
  const ScalarT model_ldvt0w,
  const ScalarT model_wdvt0w,
  const ScalarT model_pdvt0w,
  const ScalarT model_dvt1w,
  const ScalarT model_ldvt1w,
  const ScalarT model_wdvt1w,
  const ScalarT model_pdvt1w,
  const ScalarT model_dvt2w,
  const ScalarT model_ldvt2w,
  const ScalarT model_wdvt2w,
  const ScalarT model_pdvt2w,
  const ScalarT model_drout,
  const ScalarT model_ldrout,
  const ScalarT model_wdrout,
  const ScalarT model_pdrout,
  const ScalarT model_dsub,
  const ScalarT model_ldsub,
  const ScalarT model_wdsub,
  const ScalarT model_pdsub,
  const ScalarT model_vth0,
  const ScalarT model_lvth0,
  const ScalarT model_wvth0,
  const ScalarT model_pvth0,
  const ScalarT model_ua,
  const ScalarT model_lua,
  const ScalarT model_wua,
  const ScalarT model_pua,
  const ScalarT model_ua1,
  const ScalarT model_lua1,
  const ScalarT model_wua1,
  const ScalarT model_pua1,
  const ScalarT model_ub,
  const ScalarT model_lub,
  const ScalarT model_wub,
  const ScalarT model_pub,
  const ScalarT model_ub1,
  const ScalarT model_lub1,
  const ScalarT model_wub1,
  const ScalarT model_pub1,
  const ScalarT model_uc,
  const ScalarT model_luc,
  const ScalarT model_wuc,
  const ScalarT model_puc,
  const ScalarT model_uc1,
  const ScalarT model_luc1,
  const ScalarT model_wuc1,
  const ScalarT model_puc1,
  const ScalarT model_u0,
  const ScalarT model_lu0,
  const ScalarT model_wu0,
  const ScalarT model_pu0,
  const ScalarT model_ute,
  const ScalarT model_lute,
  const ScalarT model_wute,
  const ScalarT model_pute,
  const ScalarT model_voff,
  const ScalarT model_lvoff,
  const ScalarT model_wvoff,
  const ScalarT model_pvoff,
  const ScalarT model_delta,
  const ScalarT model_ldelta,
  const ScalarT model_wdelta,
  const ScalarT model_pdelta,
  const ScalarT model_rdsw,
  const ScalarT model_lrdsw,
  const ScalarT model_wrdsw,
  const ScalarT model_prdsw,
  const ScalarT model_prwg,
  const ScalarT model_lprwg,
  const ScalarT model_wprwg,
  const ScalarT model_pprwg,
  const ScalarT model_prwb,
  const ScalarT model_lprwb,
  const ScalarT model_wprwb,
  const ScalarT model_pprwb,
  const ScalarT model_prt,
  const ScalarT model_lprt,
  const ScalarT model_wprt,
  const ScalarT model_pprt,
  const ScalarT model_eta0,
  const ScalarT model_leta0,
  const ScalarT model_weta0,
  const ScalarT model_peta0,
  const ScalarT model_etab,
  const ScalarT model_letab,
  const ScalarT model_wetab,
  const ScalarT model_petab,
  const ScalarT model_pclm,
  const ScalarT model_lpclm,
  const ScalarT model_wpclm,
  const ScalarT model_ppclm,
  const ScalarT model_pdibl1,
  const ScalarT model_lpdibl1,
  const ScalarT model_wpdibl1,
  const ScalarT model_ppdibl1,
  const ScalarT model_pdibl2,
  const ScalarT model_lpdibl2,
  const ScalarT model_wpdibl2,
  const ScalarT model_ppdibl2,
  const ScalarT model_pdiblb,
  const ScalarT model_lpdiblb,
  const ScalarT model_wpdiblb,
  const ScalarT model_ppdiblb,
  const ScalarT model_pscbe1,
  const ScalarT model_lpscbe1,
  const ScalarT model_wpscbe1,
  const ScalarT model_ppscbe1,
  const ScalarT model_pscbe2,
  const ScalarT model_lpscbe2,
  const ScalarT model_wpscbe2,
  const ScalarT model_ppscbe2,
  const ScalarT model_pvag,
  const ScalarT model_lpvag,
  const ScalarT model_wpvag,
  const ScalarT model_ppvag,
  const ScalarT model_wr,
  const ScalarT model_lwr,
  const ScalarT model_wwr,
  const ScalarT model_pwr,
  const ScalarT model_dwg,
  const ScalarT model_ldwg,
  const ScalarT model_wdwg,
  const ScalarT model_pdwg,
  const ScalarT model_dwb,
  const ScalarT model_ldwb,
  const ScalarT model_wdwb,
  const ScalarT model_pdwb,
  const ScalarT model_b0,
  const ScalarT model_lb0,
  const ScalarT model_wb0,
  const ScalarT model_pb0,
  const ScalarT model_b1,
  const ScalarT model_lb1,
  const ScalarT model_wb1,
  const ScalarT model_pb1,
  const ScalarT model_alpha0,
  const ScalarT model_lalpha0,
  const ScalarT model_walpha0,
  const ScalarT model_palpha0,
  const ScalarT model_alpha1,
  const ScalarT model_lalpha1,
  const ScalarT model_walpha1,
  const ScalarT model_palpha1,
  const ScalarT model_beta0,
  const ScalarT model_lbeta0,
  const ScalarT model_wbeta0,
  const ScalarT model_pbeta0,
  const ScalarT model_elm,
  const ScalarT model_lelm,
  const ScalarT model_welm,
  const ScalarT model_pelm,
  const ScalarT model_cgsl,
  const ScalarT model_lcgsl,
  const ScalarT model_wcgsl,
  const ScalarT model_pcgsl,
  const ScalarT model_cgdl,
  const ScalarT model_lcgdl,
  const ScalarT model_wcgdl,
  const ScalarT model_pcgdl,
  const ScalarT model_ckappa,
  const ScalarT model_lckappa,
  const ScalarT model_wckappa,
  const ScalarT model_pckappa,
  const ScalarT model_cf,
  const ScalarT model_lcf,
  const ScalarT model_wcf,
  const ScalarT model_pcf,
  const ScalarT model_clc,
  const ScalarT model_lclc,
  const ScalarT model_wclc,
  const ScalarT model_pclc,
  const ScalarT model_cle,
  const ScalarT model_lcle,
  const ScalarT model_wcle,
  const ScalarT model_pcle,
  const ScalarT model_vfbcv,
  const ScalarT model_lvfbcv,
  const ScalarT model_wvfbcv,
  const ScalarT model_pvfbcv,
  const ScalarT model_acde,
  const ScalarT model_lacde,
  const ScalarT model_wacde,
  const ScalarT model_pacde,
  const ScalarT model_moin,
  const ScalarT model_lmoin,
  const ScalarT model_wmoin,
  const ScalarT model_pmoin,
  const ScalarT model_noff,
  const ScalarT model_lnoff,
  const ScalarT model_wnoff,
  const ScalarT model_pnoff,
  const ScalarT model_voffcv,
  const ScalarT model_lvoffcv,
  const ScalarT model_wvoffcv,
  const ScalarT model_pvoffcv,
  const int model_binUnit,
  const ScalarT model_tox,
  const bool model_k1Given,
  const bool model_k2Given,
  const bool model_nsubGiven,
  const bool model_xtGiven,
  const bool model_vbxGiven,
  const bool model_gamma2Given,
  const bool model_vfbGiven,
  const bool model_vth0Given,
  const int model_dtype,
  const ScalarT model_toxm,
  const ScalarT model_factor1,
  // outputs
  ScalarT & vtm, 
  ScalarT & instance_temp,
  ScalarT & jctTempSatCurDensity,
  ScalarT & jctSidewallTempSatCurDensity,
  ScalarT & unitAreaJctCapTemp,
  ScalarT & unitLengthSidewallJctCapTemp,
  ScalarT & unitLengthGateSidewallJctCapTemp,
  ScalarT & PhiBTemp,
  ScalarT & cgso,
  ScalarT & cgdo,
  ScalarT & vjsm,
  ScalarT & IsEvjsm,
  ScalarT & vjdm,
  ScalarT & IsEvjdm,
  bool & updateTemperatureCalled_,
  SizeDependParam<ScalarT> & sizeDepParams
  )
{
  ScalarT tmp, tmp1, tmp2, tmp3, Eg;
  ScalarT T0, T1, T2, T3, T4, T5, Ldrn, Wdrn;
  ScalarT delTemp, TRatio, Inv_L, Inv_W, Inv_LW;
  ScalarT Tnom;
  ScalarT Nvtm, SourceSatCurrent, DrainSatCurrent;

  bool bsuccess = true;

  // first set the instance temperature to the new temperature:
  if (temp_tmp != -999.0) 
  {
    instance_temp = temp_tmp;
    instance_temp += dtemp;
  }

  Tnom = model_tnom;
  TRatio = instance_temp/Tnom;

  vtm = CONSTKoverQ * instance_temp;
  Eg = CONSTEg0 - CONSTalphaEg * instance_temp * instance_temp / (instance_temp + CONSTbetaEg);

  if (instance_temp != Tnom)
  {
    T0 = Eg0 / Vtm0 - Eg/vtm + model_jctTempExponent*log(instance_temp/Tnom);
    T1 = exp(T0 / model_jctEmissionCoeff);
    jctTempSatCurDensity         = model_jctSatCurDensity* T1;
    jctSidewallTempSatCurDensity = model_jctSidewallSatCurDensity * T1;
  }
  else
  {
    jctTempSatCurDensity         = model_jctSatCurDensity;
    jctSidewallTempSatCurDensity = model_jctSidewallSatCurDensity;
  }

  if (jctTempSatCurDensity < 0.0)         jctTempSatCurDensity = 0.0;
  if (jctSidewallTempSatCurDensity < 0.0) jctSidewallTempSatCurDensity = 0.0;


  // Temperature dependence of D/B and S/B diode capacitance begins
  delTemp = instance_temp - Tnom;
  T0 = model_tcj * delTemp;

  if (T0 >= -1.0)
  {
    unitAreaJctCapTemp = model_unitAreaJctCap *(1.0 + T0);
  }
  else if (unitAreaJctCapTemp > 0.0)
  {
    unitAreaJctCapTemp = 0.0;

    lout() << "Temperature effect has caused cj to be negative. Cj clamped to zero.\n" << std::endl;
  }

  T0 = model_tcjsw * delTemp;

  if (T0 >= -1.0)
  {
    unitLengthSidewallJctCapTemp =
      model_unitLengthSidewallJctCap *(1.0 + T0);
  }
  else if (unitLengthSidewallJctCapTemp > 0.0)
  {
    unitLengthSidewallJctCapTemp = 0.0;
    lout() << "Temperature effect has caused cjsw to be negative. Cjsw clamped to zero.\n" << std::endl;
  }

  T0 = model_tcjswg * delTemp;

  if (T0 >= -1.0)
  {
    unitLengthGateSidewallJctCapTemp =
      model_unitLengthGateSidewallJctCap *(1.0 + T0);
  }
  else if (unitLengthGateSidewallJctCapTemp > 0.0)
  {
    unitLengthGateSidewallJctCapTemp = 0.0;
    lout() << "Temperature effect has caused cjswg to be negative. Cjswg clamped to zero.\n" << std::endl;
  }

  PhiBTemp = model_bulkJctPotential - model_tpb * delTemp;

  if (PhiBTemp < 0.01)
  {
    PhiBTemp = 0.01;
    lout() << "Temperature effect has caused pb to be < 0.01. Pb clamped to 0.01.\n" << std::endl;
  }

  PhiBSWTemp = model_sidewallJctPotential - model_tpbsw * delTemp;

  if (PhiBSWTemp <= 0.01)
  {
    PhiBSWTemp = 0.01;
    lout() << "Temperature effect has caused pbsw to be < 0.01. Pbsw clamped to 0.01.\n" << std::endl;
  }

  PhiBSWGTemp = model_GatesidewallJctPotential - model_tpbswg * delTemp;

  if (PhiBSWGTemp <= 0.01)
  {
    PhiBSWGTemp = 0.01;
    lout() << "Temperature effect has caused pbswg to be < 0.01. Pbswg clamped to 0.01.\n" << std::endl;
  }
  // End of junction capacitance

    // for sensitivity calculation, simply assume that we have to create a new size dep param, no matter what.
    // This is not efficient, but there you go.
  {
    sizeDepParams.referenceTemperature = temp_tmp;

    Ldrn = l;
    Wdrn = w;
    sizeDepParams.Length = Ldrn;
    sizeDepParams.Width = Wdrn;

    T0 = pow(Ldrn, model_Lln);
    T1 = pow(Wdrn, model_Lwn);

    tmp1 = model_Ll / T0 + model_Lw / T1
      + model_Lwl / (T0 * T1);

    sizeDepParams.dl = model_Lint + tmp1;

    tmp2 = model_Llc / T0 + model_Lwc / T1
      + model_Lwlc / (T0 * T1);

    sizeDepParams.dlc = model_dlc + tmp2;

    T2 = pow(Ldrn, model_Wln);
    T3 = pow(Wdrn, model_Wwn);

    tmp1 = model_Wl / T2 + model_Ww / T3
      + model_Wwl / (T2 * T3);

    sizeDepParams.dw = model_Wint + tmp1;
    tmp2 = model_Wlc / T2 + model_Wwc / T3
      + model_Wwlc / (T2 * T3);

    sizeDepParams.dwc = model_dwc + tmp2;

    sizeDepParams.leff = l - 2.0 * sizeDepParams.dl;
    sizeDepParams.weff = w - 2.0 * sizeDepParams.dw;
    sizeDepParams.leffCV = l - 2.0 * sizeDepParams.dlc;
    sizeDepParams.weffCV = w - 2.0 * sizeDepParams.dwc;

    if (model_binUnit == 1)
    {
      Inv_L = 1.0e-6 / sizeDepParams.leff;
      Inv_W = 1.0e-6 / sizeDepParams.weff;
      Inv_LW = 1.0e-12 / (sizeDepParams.leff * sizeDepParams.weff);
    }
    else
    {
      Inv_L = 1.0 / sizeDepParams.leff;
      Inv_W = 1.0 / sizeDepParams.weff;
      Inv_LW = 1.0 / (sizeDepParams.leff * sizeDepParams.weff);
    }

    sizeDepParams.cdsc = model_cdsc
      + model_lcdsc * Inv_L
      + model_wcdsc * Inv_W
      + model_pcdsc * Inv_LW;

    sizeDepParams.cdscb = model_cdscb
      + model_lcdscb * Inv_L
      + model_wcdscb * Inv_W
      + model_pcdscb * Inv_LW;

    sizeDepParams.cdscd = model_cdscd
      + model_lcdscd * Inv_L
      + model_wcdscd * Inv_W
      + model_pcdscd * Inv_LW;

    sizeDepParams.cit = model_cit
      + model_lcit * Inv_L
      + model_wcit * Inv_W
      + model_pcit * Inv_LW;

    sizeDepParams.nfactor = model_nfactor
      + model_lnfactor * Inv_L
      + model_wnfactor * Inv_W
      + model_pnfactor * Inv_LW;

    sizeDepParams.xj = model_xj
      + model_lxj * Inv_L
      + model_wxj * Inv_W
      + model_pxj * Inv_LW;

    sizeDepParams.vsat = model_vsat
      + model_lvsat * Inv_L
      + model_wvsat * Inv_W
      + model_pvsat * Inv_LW;

    sizeDepParams.at = model_at
      + model_lat * Inv_L
      + model_wat * Inv_W
      + model_pat * Inv_LW;

    sizeDepParams.a0 = model_a0
      + model_la0 * Inv_L
      + model_wa0 * Inv_W
      + model_pa0 * Inv_LW;

    sizeDepParams.ags = model_ags
      + model_lags * Inv_L
      + model_wags * Inv_W
      + model_pags * Inv_LW;

    sizeDepParams.a1 = model_a1
      + model_la1 * Inv_L
      + model_wa1 * Inv_W
      + model_pa1 * Inv_LW;

    sizeDepParams.a2 = model_a2
      + model_la2 * Inv_L
      + model_wa2 * Inv_W
      + model_pa2 * Inv_LW;

    sizeDepParams.keta = model_keta
      + model_lketa * Inv_L
      + model_wketa * Inv_W
      + model_pketa * Inv_LW;

    sizeDepParams.nsub = model_nsub
      + model_lnsub * Inv_L
      + model_wnsub * Inv_W
      + model_pnsub * Inv_LW;

    sizeDepParams.npeak = model_npeak
      + model_lnpeak * Inv_L
      + model_wnpeak * Inv_W
      + model_pnpeak * Inv_LW;

    sizeDepParams.ngate = model_ngate
      + model_lngate * Inv_L
      + model_wngate * Inv_W
      + model_pngate * Inv_LW;

    sizeDepParams.gamma1 = model_gamma1
      + model_lgamma1 * Inv_L
      + model_wgamma1 * Inv_W
      + model_pgamma1 * Inv_LW;

    sizeDepParams.gamma2 = model_gamma2
      + model_lgamma2 * Inv_L
      + model_wgamma2 * Inv_W
      + model_pgamma2 * Inv_LW;

    sizeDepParams.vbx = model_vbx
      + model_lvbx * Inv_L
      + model_wvbx * Inv_W
      + model_pvbx * Inv_LW;

    sizeDepParams.vbm = model_vbm
      + model_lvbm * Inv_L
      + model_wvbm * Inv_W
      + model_pvbm * Inv_LW;

    sizeDepParams.xt = model_xt
      + model_lxt * Inv_L
      + model_wxt * Inv_W
      + model_pxt * Inv_LW;

    sizeDepParams.vfb = model_vfb
      + model_lvfb * Inv_L
      + model_wvfb * Inv_W
      + model_pvfb * Inv_LW;

    sizeDepParams.k1 = model_k1
      + model_lk1 * Inv_L
      + model_wk1 * Inv_W
      + model_pk1 * Inv_LW;

    sizeDepParams.kt1 = model_kt1
      + model_lkt1 * Inv_L
      + model_wkt1 * Inv_W
      + model_pkt1 * Inv_LW;

    sizeDepParams.kt1l = model_kt1l
      + model_lkt1l * Inv_L
      + model_wkt1l * Inv_W
      + model_pkt1l * Inv_LW;

    sizeDepParams.k2 = model_k2
      + model_lk2 * Inv_L
      + model_wk2 * Inv_W
      + model_pk2 * Inv_LW;

    sizeDepParams.kt2 = model_kt2
      + model_lkt2 * Inv_L
      + model_wkt2 * Inv_W
      + model_pkt2 * Inv_LW;

    sizeDepParams.k3 = model_k3
      + model_lk3 * Inv_L
      + model_wk3 * Inv_W
      + model_pk3 * Inv_LW;

    sizeDepParams.k3b = model_k3b
      + model_lk3b * Inv_L
      + model_wk3b * Inv_W
      + model_pk3b * Inv_LW;

    sizeDepParams.w0 = model_w0
      + model_lw0 * Inv_L
      + model_ww0 * Inv_W
      + model_pw0 * Inv_LW;

    sizeDepParams.nlx = model_nlx
      + model_lnlx * Inv_L
      + model_wnlx * Inv_W
      + model_pnlx * Inv_LW;

    sizeDepParams.dvt0 = model_dvt0
      + model_ldvt0 * Inv_L
      + model_wdvt0 * Inv_W
      + model_pdvt0 * Inv_LW;

    sizeDepParams.dvt1 = model_dvt1
      + model_ldvt1 * Inv_L
      + model_wdvt1 * Inv_W
      + model_pdvt1 * Inv_LW;

    sizeDepParams.dvt2 = model_dvt2
      + model_ldvt2 * Inv_L
      + model_wdvt2 * Inv_W
      + model_pdvt2 * Inv_LW;

    sizeDepParams.dvt0w = model_dvt0w
      + model_ldvt0w * Inv_L
      + model_wdvt0w * Inv_W
      + model_pdvt0w * Inv_LW;

    sizeDepParams.dvt1w = model_dvt1w
      + model_ldvt1w * Inv_L
      + model_wdvt1w * Inv_W
      + model_pdvt1w * Inv_LW;

    sizeDepParams.dvt2w = model_dvt2w
      + model_ldvt2w * Inv_L
      + model_wdvt2w * Inv_W
      + model_pdvt2w * Inv_LW;

    sizeDepParams.drout = model_drout
      + model_ldrout * Inv_L
      + model_wdrout * Inv_W
      + model_pdrout * Inv_LW;

    sizeDepParams.dsub = model_dsub
      + model_ldsub * Inv_L
      + model_wdsub * Inv_W
      + model_pdsub * Inv_LW;

    sizeDepParams.vth0 = model_vth0
      + model_lvth0 * Inv_L
      + model_wvth0 * Inv_W
      + model_pvth0 * Inv_LW;

    sizeDepParams.ua = model_ua
      + model_lua * Inv_L
      + model_wua * Inv_W
      + model_pua * Inv_LW;

    sizeDepParams.ua1 = model_ua1
      + model_lua1 * Inv_L
      + model_wua1 * Inv_W
      + model_pua1 * Inv_LW;

    sizeDepParams.ub = model_ub
      + model_lub * Inv_L
      + model_wub * Inv_W
      + model_pub * Inv_LW;

    sizeDepParams.ub1 = model_ub1
      + model_lub1 * Inv_L
      + model_wub1 * Inv_W
      + model_pub1 * Inv_LW;

    sizeDepParams.uc = model_uc
      + model_luc * Inv_L
      + model_wuc * Inv_W
      + model_puc * Inv_LW;

    sizeDepParams.uc1 = model_uc1
      + model_luc1 * Inv_L
      + model_wuc1 * Inv_W
      + model_puc1 * Inv_LW;

    sizeDepParams.u0 = model_u0
      + model_lu0 * Inv_L
      + model_wu0 * Inv_W
      + model_pu0 * Inv_LW;

    sizeDepParams.ute = model_ute
      + model_lute * Inv_L
      + model_wute * Inv_W
      + model_pute * Inv_LW;

    sizeDepParams.voff = model_voff
      + model_lvoff * Inv_L
      + model_wvoff * Inv_W
      + model_pvoff * Inv_LW;

    sizeDepParams.delta = model_delta
      + model_ldelta * Inv_L
      + model_wdelta * Inv_W
      + model_pdelta * Inv_LW;

    sizeDepParams.rdsw = model_rdsw
      + model_lrdsw * Inv_L
      + model_wrdsw * Inv_W
      + model_prdsw * Inv_LW;

    sizeDepParams.prwg = model_prwg
      + model_lprwg * Inv_L
      + model_wprwg * Inv_W
      + model_pprwg * Inv_LW;

    sizeDepParams.prwb = model_prwb
      + model_lprwb * Inv_L
      + model_wprwb * Inv_W
      + model_pprwb * Inv_LW;

    sizeDepParams.prt = model_prt
      + model_lprt * Inv_L
      + model_wprt * Inv_W
      + model_pprt * Inv_LW;

    sizeDepParams.eta0 = model_eta0
      + model_leta0 * Inv_L
      + model_weta0 * Inv_W
      + model_peta0 * Inv_LW;

    sizeDepParams.etab = model_etab
      + model_letab * Inv_L
      + model_wetab * Inv_W
      + model_petab * Inv_LW;

    sizeDepParams.pclm = model_pclm
      + model_lpclm * Inv_L
      + model_wpclm * Inv_W
      + model_ppclm * Inv_LW;

    sizeDepParams.pdibl1 = model_pdibl1
      + model_lpdibl1 * Inv_L
      + model_wpdibl1 * Inv_W
      + model_ppdibl1 * Inv_LW;

    sizeDepParams.pdibl2 = model_pdibl2
      + model_lpdibl2 * Inv_L
      + model_wpdibl2 * Inv_W
      + model_ppdibl2 * Inv_LW;

    sizeDepParams.pdiblb = model_pdiblb
      + model_lpdiblb * Inv_L
      + model_wpdiblb * Inv_W
      + model_ppdiblb * Inv_LW;

    sizeDepParams.pscbe1 = model_pscbe1
      + model_lpscbe1 * Inv_L
      + model_wpscbe1 * Inv_W
      + model_ppscbe1 * Inv_LW;

    sizeDepParams.pscbe2 = model_pscbe2
      + model_lpscbe2 * Inv_L
      + model_wpscbe2 * Inv_W
      + model_ppscbe2 * Inv_LW;

    sizeDepParams.pvag = model_pvag
      + model_lpvag * Inv_L
      + model_wpvag * Inv_W
      + model_ppvag * Inv_LW;

    sizeDepParams.wr = model_wr
      + model_lwr * Inv_L
      + model_wwr * Inv_W
      + model_pwr * Inv_LW;

    sizeDepParams.dwg = model_dwg
      + model_ldwg * Inv_L
      + model_wdwg * Inv_W
      + model_pdwg * Inv_LW;

    sizeDepParams.dwb = model_dwb
      + model_ldwb * Inv_L
      + model_wdwb * Inv_W
      + model_pdwb * Inv_LW;

    sizeDepParams.b0 = model_b0
      + model_lb0 * Inv_L
      + model_wb0 * Inv_W
      + model_pb0 * Inv_LW;

    sizeDepParams.b1 = model_b1
      + model_lb1 * Inv_L
      + model_wb1 * Inv_W
      + model_pb1 * Inv_LW;

    sizeDepParams.alpha0 = model_alpha0
      + model_lalpha0 * Inv_L
      + model_walpha0 * Inv_W
      + model_palpha0 * Inv_LW;

    sizeDepParams.alpha1 = model_alpha1
      + model_lalpha1 * Inv_L
      + model_walpha1 * Inv_W
      + model_palpha1 * Inv_LW;

    sizeDepParams.beta0 = model_beta0
      + model_lbeta0 * Inv_L
      + model_wbeta0 * Inv_W
      + model_pbeta0 * Inv_LW;

    // CV model
    sizeDepParams.elm = model_elm
      + model_lelm * Inv_L
      + model_welm * Inv_W
      + model_pelm * Inv_LW;

    sizeDepParams.cgsl = model_cgsl
      + model_lcgsl * Inv_L
      + model_wcgsl * Inv_W
      + model_pcgsl * Inv_LW;

    sizeDepParams.cgdl = model_cgdl
      + model_lcgdl * Inv_L
      + model_wcgdl * Inv_W
      + model_pcgdl * Inv_LW;

    sizeDepParams.ckappa = model_ckappa
      + model_lckappa * Inv_L
      + model_wckappa * Inv_W
      + model_pckappa * Inv_LW;

    sizeDepParams.cf = model_cf
      + model_lcf * Inv_L
      + model_wcf * Inv_W
      + model_pcf * Inv_LW;

    sizeDepParams.clc = model_clc
      + model_lclc * Inv_L
      + model_wclc * Inv_W
      + model_pclc * Inv_LW;

    sizeDepParams.cle = model_cle
      + model_lcle * Inv_L
      + model_wcle * Inv_W
      + model_pcle * Inv_LW;

    sizeDepParams.vfbcv = model_vfbcv
      + model_lvfbcv * Inv_L
      + model_wvfbcv * Inv_W
      + model_pvfbcv * Inv_LW;

    sizeDepParams.acde = model_acde
      + model_lacde * Inv_L
      + model_wacde * Inv_W
      + model_pacde * Inv_LW;

    sizeDepParams.moin = model_moin
      + model_lmoin * Inv_L
      + model_wmoin * Inv_W
      + model_pmoin * Inv_LW;

    sizeDepParams.noff = model_noff
      + model_lnoff * Inv_L
      + model_wnoff * Inv_W
      + model_pnoff * Inv_LW;

    sizeDepParams.voffcv = model_voffcv
      + model_lvoffcv * Inv_L
      + model_wvoffcv * Inv_W
      + model_pvoffcv * Inv_LW;

    sizeDepParams.abulkCVfactor = 1.0
      + pow((sizeDepParams.clc / sizeDepParams.leffCV), sizeDepParams.cle);

    T0 = (TRatio - 1.0);

    sizeDepParams.ua = sizeDepParams.ua + sizeDepParams.ua1 * T0;
    sizeDepParams.ub = sizeDepParams.ub + sizeDepParams.ub1 * T0;
    sizeDepParams.uc = sizeDepParams.uc + sizeDepParams.uc1 * T0;

    if (sizeDepParams.u0 > 1.0) sizeDepParams.u0 = sizeDepParams.u0 / 1.0e4;

    sizeDepParams.u0temp = sizeDepParams.u0 * pow(TRatio, sizeDepParams.ute);

    sizeDepParams.vsattemp = sizeDepParams.vsat - sizeDepParams.at * T0;

    sizeDepParams.rds0 = (sizeDepParams.rdsw + sizeDepParams.prt * T0)
      / pow(sizeDepParams.weff * 1E6, sizeDepParams.wr);

    sizeDepParams.cgdo = (model_cgdo + sizeDepParams.cf) * sizeDepParams.weffCV;
    sizeDepParams.cgso = (model_cgso + sizeDepParams.cf) * sizeDepParams.weffCV;
    sizeDepParams.cgbo = model_cgbo * sizeDepParams.leffCV;

    T0 = sizeDepParams.leffCV * sizeDepParams.leffCV;

    sizeDepParams.tconst = sizeDepParams.u0temp * sizeDepParams.elm / (model_cox
        * sizeDepParams.weffCV * sizeDepParams.leffCV * T0);

    if (!model_npeakGiven && model_gamma1Given)
    {
      T0 = sizeDepParams.gamma1 * model_cox;
      sizeDepParams.npeak = 3.021E22 * T0 * T0;
    }

    sizeDepParams.phi     = 2.0 * Vtm0 * log(sizeDepParams.npeak / ni);
    sizeDepParams.sqrtPhi = sqrt(sizeDepParams.phi);
    sizeDepParams.phis3   = sizeDepParams.sqrtPhi * sizeDepParams.phi;

    sizeDepParams.Xdep0 = sqrt(2.0 * CONSTEPSSI / (CONSTQ * sizeDepParams.npeak * 1.0e6))
      * sizeDepParams.sqrtPhi;

    sizeDepParams.sqrtXdep0 = sqrt(sizeDepParams.Xdep0);
    sizeDepParams.litl = sqrt(3.0 * sizeDepParams.xj * model_tox);

    sizeDepParams.vbi = Vtm0 * log(1.0e20 * sizeDepParams.npeak / (ni * ni));

    sizeDepParams.cdep0 = sqrt(CONSTQ * CONSTEPSSI * sizeDepParams.npeak * 1.0e6 / 2.0
        / sizeDepParams.phi);

    sizeDepParams.ldeb = sqrt(CONSTEPSSI * Vtm0 / (CONSTQ
          * sizeDepParams.npeak * 1.0e6)) / 3.0;

    sizeDepParams.acde *= pow((sizeDepParams.npeak / 2.0e16), -0.25);


    if (model_k1Given || model_k2Given)
    {
      if (!model_k1Given)
      {
        sizeDepParams.k1 = 0.53;
      }

      if (!model_k2Given)
      {
        sizeDepParams.k2 = -0.0186;
      }
    }
    else
    {
      if (!model_vbxGiven)
        sizeDepParams.vbx = sizeDepParams.phi - 7.7348e-4 * sizeDepParams.npeak
          * sizeDepParams.xt * sizeDepParams.xt;

      if (sizeDepParams.vbx > 0.0)
        sizeDepParams.vbx = -sizeDepParams.vbx;

      if (sizeDepParams.vbm > 0.0)
        sizeDepParams.vbm = -sizeDepParams.vbm;

      if (!model_gamma1Given)
        sizeDepParams.gamma1 = 5.753e-12 * sqrt(sizeDepParams.npeak) / model_cox;

      if (!model_gamma2Given)
        sizeDepParams.gamma2 = 5.753e-12 * sqrt(sizeDepParams.nsub) / model_cox;

      T0 = sizeDepParams.gamma1 - sizeDepParams.gamma2;
      T1 = sqrt(sizeDepParams.phi - sizeDepParams.vbx) - sizeDepParams.sqrtPhi;
      T2 = sqrt(sizeDepParams.phi * (sizeDepParams.phi - sizeDepParams.vbm)) - sizeDepParams.phi;

      sizeDepParams.k2 = T0 * T1 / (2.0 * T2 + sizeDepParams.vbm);
      sizeDepParams.k1 = sizeDepParams.gamma2 - 2.0 * sizeDepParams.k2 * sqrt(sizeDepParams.phi
          - sizeDepParams.vbm);
    }

    if (sizeDepParams.k2 < 0.0)
    {
      T0 = 0.5 * sizeDepParams.k1 / sizeDepParams.k2;
      sizeDepParams.vbsc = 0.9 * (sizeDepParams.phi - T0 * T0);

      if (sizeDepParams.vbsc > -3.0) sizeDepParams.vbsc = -3.0;
      else if (sizeDepParams.vbsc < -30.0) sizeDepParams.vbsc = -30.0;
    }
    else
    {
      sizeDepParams.vbsc = -30.0;
    }

    if (sizeDepParams.vbsc > sizeDepParams.vbm) sizeDepParams.vbsc = sizeDepParams.vbm;

    if (!model_vfbGiven)
    {
      if (model_vth0Given)
      {
        sizeDepParams.vfb = model_dtype * sizeDepParams.vth0
          - sizeDepParams.phi - sizeDepParams.k1 * sizeDepParams.sqrtPhi;
      }
      else
      {   sizeDepParams.vfb = -1.0;
      }
    }

    if (!model_vth0Given)
    {
      sizeDepParams.vth0 = model_dtype
        * (sizeDepParams.vfb + sizeDepParams.phi + sizeDepParams.k1
            * sizeDepParams.sqrtPhi);
    }

    sizeDepParams.k1ox = sizeDepParams.k1 * model_tox / model_toxm;
    sizeDepParams.k2ox = sizeDepParams.k2 * model_tox / model_toxm;

    T1 = sqrt(CONSTEPSSI / CONSTEPSOX * model_tox * sizeDepParams.Xdep0);
    T0 = exp(-0.5 * sizeDepParams.dsub * sizeDepParams.leff / T1);

    sizeDepParams.theta0vb0 = (T0 + 2.0 * T0 * T0);

    T0 = exp(-0.5 * sizeDepParams.drout * sizeDepParams.leff / T1);
    T2 = (T0 + 2.0 * T0 * T0);

    sizeDepParams.thetaRout = sizeDepParams.pdibl1 * T2 + sizeDepParams.pdibl2;

    tmp = sqrt(sizeDepParams.Xdep0);
    tmp1 = sizeDepParams.vbi - sizeDepParams.phi;
    tmp2 = model_factor1 * tmp;

    T0 = -0.5 * sizeDepParams.dvt1w * sizeDepParams.weff * sizeDepParams.leff / tmp2;

    if (T0 > -CONSTEXP_THRESHOLD)
    {
      T1 = exp(T0);
      T2 = T1 * (1.0 + 2.0 * T1);
    }
    else
    {
      T1 = CONSTMIN_EXP;
      T2 = T1 * (1.0 + 2.0 * T1);
    }
    T0 = sizeDepParams.dvt0w * T2;
    T2 = T0 * tmp1;

    T0 = -0.5 * sizeDepParams.dvt1 * sizeDepParams.leff / tmp2;

    if (T0 > -CONSTEXP_THRESHOLD)
    {
      T1 = exp(T0);
      T3 = T1 * (1.0 + 2.0 * T1);
    }
    else
    {
      T1 = CONSTMIN_EXP;
      T3 = T1 * (1.0 + 2.0 * T1);
    }

    T3 = sizeDepParams.dvt0 * T3 * tmp1;

    T4 = model_tox * sizeDepParams.phi / (sizeDepParams.weff + sizeDepParams.w0);

    T0 = sqrt(1.0 + sizeDepParams.nlx / sizeDepParams.leff);
    T5 = sizeDepParams.k1ox * (T0 - 1.0) * sizeDepParams.sqrtPhi
      + (sizeDepParams.kt1 + sizeDepParams.kt1l / sizeDepParams.leff) * (TRatio - 1.0);

    tmp3 = model_dtype * sizeDepParams.vth0 - T2 - T3 + sizeDepParams.k3 * T4 + T5;

    sizeDepParams.vfbzb = tmp3 - sizeDepParams.phi - sizeDepParams.k1 * sizeDepParams.sqrtPhi;

  } // End of vfbzb

  cgso = sizeDepParams.cgso;
  cgdo = sizeDepParams.cgdo;

  Nvtm = vtm * model_jctEmissionCoeff;

  if ((sourceArea <= 0.0) &&
      (sourcePerimeter <= 0.0))
  {
    SourceSatCurrent = 1.0e-14;
  }
  else
  {
    SourceSatCurrent = sourceArea * jctTempSatCurDensity
      + sourcePerimeter
      * jctSidewallTempSatCurDensity;
  }

  if ((SourceSatCurrent > 0.0) && (model_ijth > 0.0))
  {
    vjsm = Nvtm * log(model_ijth / SourceSatCurrent + 1.0);
    IsEvjsm = SourceSatCurrent * exp(vjsm / Nvtm);
  }

  if ((drainArea <= 0.0) &&
      (drainPerimeter <= 0.0))
  {
    DrainSatCurrent = 1.0e-14;
  }
  else
  {
    DrainSatCurrent = drainArea * jctTempSatCurDensity
      + drainPerimeter
      * jctSidewallTempSatCurDensity;
  }

  if ((DrainSatCurrent > 0.0) && (model_ijth > 0.0))
  {
    vjdm = Nvtm * log(model_ijth / DrainSatCurrent + 1.0);
    IsEvjdm = DrainSatCurrent * exp(vjdm / Nvtm);
  }

  updateTemperatureCalled_ = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
template <typename ScalarT> 
bool processParams (
  const DeviceOptions & devOptions,

   bool given_TEMP,
   bool given_DTEMP,
   bool given_L,
   bool given_W,
   bool given_AD,
   bool given_AS,
   bool given_NQSMOD,
   // inputs
   const ScalarT & model_l,
   const ScalarT & model_w,
   const ScalarT & model_sheetResistance,
   const ScalarT & drainSquares,
   const ScalarT & sourceSquares,
   // outputs
   ScalarT & temp,
   ScalarT & dtemp,
   ScalarT & l,
   ScalarT & w,
   ScalarT & drainArea,
   ScalarT & sourceArea,
   ScalarT & drainConductance,
   ScalarT & sourceConductance,
   int & nqsMod
    )
{
  // Set any non-constant parameter defaults:
  if (!given_TEMP)
  {
    temp = devOptions.temp.getImmutableValue<double>();
    if (!given_DTEMP)
      dtemp = 0.0;
  }
  else
  {
    dtemp = 0.0;
    if (given_DTEMP)
    {
      // don't issue a warning in this templated version of the function, as it is only called for .sens
    }
  }

  if (!given_L)
    l =model_l;
  if (!given_W)
    w = model_w;
  if (!given_AD)
    drainArea = devOptions.defad;
  if (!given_AS)
    sourceArea = devOptions.defas;

  // process source/drain series resistance
  drainConductance = model_sheetResistance * drainSquares;

  if (drainConductance > 0.0)
    drainConductance = 1.0 / drainConductance;
  else
    drainConductance = 0.0;

  sourceConductance = model_sheetResistance * sourceSquares;

  if (sourceConductance > 0.0)
    sourceConductance = 1.0 / sourceConductance;
  else
    sourceConductance = 0.0;

  if (given_NQSMOD)
  {
    nqsMod = 0;
  }

  // If there are any time dependent parameters, set their values at for
  // the current time.

  // now set the temperature related stuff.
  //updateTemperature(temp);
  //
  // this will be called separately, so I don't have to re-do the argument list

  return true;
}

//-----------------------------------------------------------------------------
// Function      : processModelParams
// Purpose       : patterned after the model version of this function.
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename ScalarT> 
bool processModelParams (
      // givens:
      bool given_TOXM, bool given_DSUB, bool given_LLC, bool given_LWC, bool given_LWLC, bool given_WLC,
      bool given_WWL, bool given_WWLC, bool given_DWC, bool given_DLC, bool given_CF, bool given_CGDO,
      bool given_CGSO, bool given_CGBO, bool given_CJSWG, bool given_PBSWG, bool given_MJSWG, 
      // inputs:
      const ScalarT tox, const ScalarT drout, const ScalarT Ll, const ScalarT Lw, const ScalarT Lwl,
      const ScalarT Wl, const ScalarT Wwl, const ScalarT Wint, const ScalarT Lint, const ScalarT tnom,
      const ScalarT cgdl, const ScalarT cgsl, const ScalarT unitLengthSidewallJctCap,
      const ScalarT bulkJctSideGradingCoeff, const ScalarT xj,

      //outputs:
      ScalarT & toxm, ScalarT & cox, ScalarT & dsub, ScalarT & Llc, ScalarT & Lwc, ScalarT & Lwlc,
      ScalarT & Wlc, ScalarT & Wwlc, ScalarT & dwc, ScalarT & dlc, 
      ScalarT & cf, ScalarT & cgdo, ScalarT & cgso, ScalarT & cgbo,
      ScalarT & unitLengthGateSidewallJctCap, ScalarT & GatesidewallJctPotential, ScalarT & bulkJctGateSideGradingCoeff,
      ScalarT & bulkJctPotential, ScalarT & sidewallJctPotential,
      ScalarT & vcrit, ScalarT & factor1, ScalarT & Vtm0, ScalarT & Eg0, ScalarT & ni  
      )
{
  cox = 3.453133e-11 / tox;
  if (!given_TOXM)    toxm = tox;
  if (!given_DSUB)  dsub = drout;
  if (!given_LLC)  Llc = Ll;
  if (!given_LWC)  Lwc = Lw;
  if (!given_LWLC) Lwlc = Lwl;
  if (!given_WLC)  Wlc = Wl;
  if (!given_WWL) Wwlc = Wwl;
  if (!given_WWLC) Wwlc = Wwl;
  if (!given_DWC)  dwc = Wint;
  if (!given_DLC)  dlc = Lint;

  if (!given_CF)
  {
    ScalarT C1 = 2.0 * CONSTEPSOX;
    ScalarT C5 = M_PI;
    ScalarT C2 = 1.0 + (0.4e-6 / tox);
    ScalarT C3 = log(C2);
    cf = C1*C3/C5;
  }

  if (!given_CGDO)
  {
    if (given_DLC && (dlc > 0.0)) cgdo = dlc * cox - cgdl ;
    else                         cgdo = 0.6 * xj * cox;
  }

  if (!given_CGSO)
  {
    if (given_DLC && (dlc > 0.0)) cgso = dlc * cox - cgsl ;
    else                         cgso = 0.6 * xj * cox;
  }

  if (!given_CGBO) cgbo = 2.0 * dwc * cox;

  if (!given_CJSWG)
    unitLengthGateSidewallJctCap = unitLengthSidewallJctCap ;

  if (!given_PBSWG)
    GatesidewallJctPotential = sidewallJctPotential;

  if (!given_MJSWG)
    bulkJctGateSideGradingCoeff = bulkJctSideGradingCoeff;


  // More initializations:  taken from b3temp.c:
  if (bulkJctPotential < 0.1)
  {
    bulkJctPotential = 0.1;
  }

  if (sidewallJctPotential < 0.1)
  {
    sidewallJctPotential = 0.1;
  }

  if (GatesidewallJctPotential < 0.1)
  {
    GatesidewallJctPotential = 0.1;
  }

  vcrit   = CONSTvt0 * log(CONSTvt0 / (CONSTroot2 * 1.0e-14));
  factor1 = sqrt(CONSTEPSSI / CONSTEPSOX * tox);

  Vtm0 = CONSTKoverQ * tnom;
  Eg0  = CONSTEg0 - CONSTalphaEg * tnom * tnom / (tnom + CONSTbetaEg);
  ni   = CONSTNi0 * (tnom / CONSTREFTEMP) * sqrt(tnom / CONSTREFTEMP)
    * exp(21.5565981 - Eg0 / (2.0 * Vtm0));

  // If there are any time dependent parameters, set their values at for
  // the current time.

  // We have changed model parameters (maybe) and so all size dependent params
  // we may have stored may be invalid.  Clear them
  //
  // For the sensitivity calculations, nothing to clear.
  //
  //clearTemperatureData();

  return true;
}


//-----------------------------------------------------------------------------
// Function      : setupCapacitors_oldDAE ()
//
// Purpose       : Same as new-DAE version, but including pdt, essentially.
//
// Special Notes : The original version of this function is still present in the
//                 model because it computes several terms used in the Jacobian
//                 loads.
//
//                 Most of these terms are completely useless for sensitivity 
//                 calculations, which is what this version of the function
//                 is for.  So, many variables are local to the function that 
//                 would otherwise have been function arguments to be used later
//                 in subsequent calculations.  As such, some parts of this 
//                 function could probably go away for efficiency.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/29/2018
//-----------------------------------------------------------------------------
template <typename ScalarT> 
bool setupCapacitors_oldDAE 
 (
  const SolverState & solverState,
  int mode,
  int nqsMod,
  const SizeDependParam<ScalarT> & sizeDepParams,
  const ScalarT &  cggb,
  const ScalarT &  cgdo,
  const ScalarT &  cgso,
  const ScalarT &  cgdb,
  const ScalarT &  cgsb,
  const ScalarT &  cdgb,
  const ScalarT &  cddb,
  const ScalarT &  capbd,
  const ScalarT &  cdsb,
  const ScalarT &  cbgb,
  const ScalarT &  cbdb,
  const ScalarT &  capbs,
  const ScalarT &  cbsb,
  const ScalarT &  qgdo,
  const ScalarT &  qgso,
  const ScalarT &  vgb,
  const ScalarT &  qcheq,
  const ScalarT &  qdef,
  const ScalarT &  ScalingFactor,
  const ScalarT & cqgb,
  const ScalarT & cqdb,
  const ScalarT & cqsb,
  const ScalarT & cqbb,
  const ScalarT & model_cox,
  const ScalarT & model_xpart,
  // outputs
  ScalarT & qgate,
  ScalarT & qbulk,
  ScalarT & qdrn,
  ScalarT & qsrc,
  ScalarT & CoxWL
  )
{
  ScalarT ag0 = solverState.pdt_;
  ScalarT T0 = 0.0;

  // ---- variables that originally were class variables, but are now local 
  // (b/c they aren't as important for sensitivities)
  //
  ScalarT qgb; // temp variable that is used
  ScalarT gtg; // pointless??  seems redundant with ggtg
  ScalarT gtd; // pointless??  seems redundant with ggtd
  ScalarT gts; // pointless??  seems redundant with ggts
  ScalarT gtb; // pointless??  seems redundant with ggtb
  ScalarT gqdef;  // doens't seem to be used anywhere


  // a lot of these "gc" terms used to be used in the Jacobian, but are not any more.  
  ScalarT   gcggb;
  ScalarT   gcgdb;
  ScalarT   gcgsb;
  ScalarT   gcdgb;
  ScalarT   gcddb;
  ScalarT   gcdsb;
  ScalarT   gcsgb;
  ScalarT   gcsdb;
  ScalarT   gcssb;
  ScalarT   gcbgb;
  ScalarT   gcbdb;
  ScalarT   gcbsb;
  ScalarT   qgd;
  ScalarT   qgs;


  // these are used in Jacobian but not elsewhere, AFAIK
  ScalarT ggtg;
  ScalarT ggtd;
  ScalarT ggtb;
  ScalarT ggts;
  ScalarT sxpart;
  ScalarT dxpart;
  ScalarT ddxpart_dVd;
  ScalarT ddxpart_dVg;
  ScalarT ddxpart_dVb;
  ScalarT ddxpart_dVs;
  ScalarT dsxpart_dVd;
  ScalarT dsxpart_dVg;
  ScalarT dsxpart_dVb;
  ScalarT dsxpart_dVs;
  ScalarT gcqgb;
  ScalarT gcqdb;
  ScalarT gcqsb;
  ScalarT gcqbb;

  ScalarT Cdd;
  ScalarT Csd;
  ScalarT Cdg;
  ScalarT Csg;
  ScalarT Cds;
  ScalarT Css;

  // ---- end of variables that originally were class variables, but are now local 

  // It is necessary to set ag0=0.0, because for the first time step out of
  // the DCOP, all the time derivatives are forced to be zero.  Thus, all
  // their derivatives should also be zero.  If it wasn't for that, then ag0
  // could always be pdt.  (it used to be, before the -jacobian_test capability).
  if (!(solverState.dcopFlag) && solverState.initTranFlag_ && solverState.newtonIter==0)
  {
    ag0 = 0.0;
  }

  if (mode > 0)
  {
    if (nqsMod == 0)
    {
      gcggb = (cggb + cgdo + cgso + sizeDepParams.cgbo ) * ag0;
      gcgdb = (cgdb - cgdo) * ag0;
      gcgsb = (cgsb - cgso) * ag0;

      gcdgb = (cdgb - cgdo) * ag0;
      gcddb = (cddb + capbd + cgdo) * ag0;
      gcdsb = cdsb * ag0;

      gcsgb = -(cggb + cbgb + cdgb + cgso) * ag0;
      gcsdb = -(cgdb + cbdb + cddb) * ag0;
      gcssb = (capbs + cgso - (cgsb + cbsb + cdsb)) * ag0;

      gcbgb = (cbgb - sizeDepParams.cgbo) * ag0;
      gcbdb = (cbdb - capbd) * ag0;
      gcbsb = (cbsb - capbs) * ag0;

      qgd = qgdo;
      qgs = qgso;
      qgb = sizeDepParams.cgbo * vgb;
      qgate += qgd + qgs + qgb;
      qbulk -= qgb;
      qdrn -= qgd;
      qsrc = -(qgate + qbulk + qdrn);

      ggtg = ggtd = ggtb = ggts = 0.0;
      sxpart = 0.6;
      dxpart = 0.4;
      ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
      dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
    }
    else  // nqsMode != 0
    {
      if (qcheq > 0.0)
      {
        T0 = sizeDepParams.tconst * qdef * ScalingFactor;
      }
      else
      {
        T0 = -sizeDepParams.tconst * qdef * ScalingFactor;
      }

      ggtg = gtg = T0 * cqgb;
      ggtd = gtd = T0 * cqdb;
      ggts = gts = T0 * cqsb;
      ggtb = gtb = T0 * cqbb;
      gqdef = ScalingFactor * ag0;

      gcqgb = cqgb * ag0;
      gcqdb = cqdb * ag0;
      gcqsb = cqsb * ag0;
      gcqbb = cqbb * ag0;

      gcggb = (cgdo + cgso + sizeDepParams.cgbo ) * ag0;
      gcgdb = -cgdo * ag0;
      gcgsb = -cgso * ag0;

      gcdgb = -cgdo * ag0;
      gcddb = (capbd + cgdo) * ag0;
      gcdsb = 0.0;

      gcsgb = -cgso * ag0;
      gcsdb = 0.0;
      gcssb = (capbs + cgso) * ag0;

      gcbgb = -sizeDepParams.cgbo * ag0;
      gcbdb = -capbd * ag0;
      gcbsb = -capbs * ag0;

      CoxWL = model_cox * sizeDepParams.weffCV * sizeDepParams.leffCV;

      if (fabs(qcheq) <= 1.0e-5 * CoxWL)
      {
        if      (model_xpart < 0.5) dxpart = 0.4;
        else if (model_xpart > 0.5) dxpart = 0.0;
        else                               dxpart = 0.5;

        ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
      }
      else
      {
        dxpart = qdrn / qcheq;
        Cdd = cddb;
        Csd = -(cgdb + cddb + cbdb);
        ddxpart_dVd = (Cdd - dxpart * (Cdd + Csd)) / qcheq;
        Cdg = cdgb;
        Csg = -(cggb + cdgb + cbgb);
        ddxpart_dVg = (Cdg - dxpart * (Cdg + Csg)) / qcheq;

        Cds = cdsb;
        Css = -(cgsb + cdsb + cbsb);
        ddxpart_dVs = (Cds - dxpart * (Cds + Css)) / qcheq;

        ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);
      }
      sxpart = 1.0 - dxpart;
      dsxpart_dVd = -ddxpart_dVd;
      dsxpart_dVg = -ddxpart_dVg;
      dsxpart_dVs = -ddxpart_dVs;
      dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);

      qgd = qgdo;
      qgs = qgso;
      qgb = sizeDepParams.cgbo * vgb;
      qgate = qgd + qgs + qgb;
      qbulk = -qgb;
      qdrn = -qgd;
      qsrc = -(qgate + qbulk + qdrn);
    } // nqsMod
  }
  else
  {
    if (nqsMod == 0)
    {
      gcggb = (cggb + cgdo + cgso + sizeDepParams.cgbo ) * ag0;
      gcgdb = (cgsb - cgdo) * ag0;
      gcgsb = (cgdb - cgso) * ag0;

      gcdgb = -(cggb + cbgb + cdgb + cgdo) * ag0;
      gcddb = (capbd + cgdo - (cgsb + cbsb + cdsb)) * ag0;
      gcdsb = -(cgdb + cbdb + cddb) * ag0;

      gcsgb = (cdgb - cgso) * ag0;
      gcsdb = cdsb * ag0;
      gcssb = (cddb + capbs + cgso) * ag0;

      gcbgb = (cbgb - sizeDepParams.cgbo) * ag0;
      gcbdb = (cbsb - capbd) * ag0;
      gcbsb = (cbdb - capbs) * ag0;

      qgd = qgdo;
      qgs = qgso;
      qgb = sizeDepParams.cgbo * vgb;
      qgate += qgd + qgs + qgb;
      qbulk -= qgb;
      qsrc = qdrn - qgs;
      qdrn = -(qgate + qbulk + qsrc);

      ggtg = ggtd = ggtb = ggts = 0.0;
      sxpart = 0.4;
      dxpart = 0.6;
      ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
      dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
    }
    else  // nqsMode != 0
    {
      if (qcheq > 0.0)
      {
        T0 = sizeDepParams.tconst * qdef * ScalingFactor;
      }
      else
      {
        T0 = -sizeDepParams.tconst * qdef * ScalingFactor;
      }

      ggtg = gtg = T0 * cqgb;
      ggts = gtd = T0 * cqdb;
      ggtd = gts = T0 * cqsb;
      ggtb = gtb = T0 * cqbb;
      gqdef = ScalingFactor * ag0;

      gcqgb = cqgb * ag0;
      gcqdb = cqsb * ag0;
      gcqsb = cqdb * ag0;
      gcqbb = cqbb * ag0;

      gcggb = (cgdo + cgso + sizeDepParams.cgbo) * ag0;
      gcgdb = -cgdo * ag0;
      gcgsb = -cgso * ag0;

      gcdgb = -cgdo * ag0;
      gcddb = (capbd + cgdo) * ag0;
      gcdsb = 0.0;

      gcsgb = -cgso * ag0;
      gcsdb = 0.0;
      gcssb = (capbs + cgso) * ag0;

      gcbgb = -sizeDepParams.cgbo * ag0;
      gcbdb = -capbd * ag0;
      gcbsb = -capbs * ag0;

      CoxWL = model_cox * sizeDepParams.weffCV * sizeDepParams.leffCV;

      if (fabs(qcheq) <= 1.0e-5 * CoxWL)
      {
        if      (model_xpart < 0.5) sxpart = 0.4;
        else if (model_xpart > 0.5) sxpart = 0.0;
        else                         sxpart = 0.5;

        dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
      }
      else
      {
        sxpart = qdrn / qcheq;
        Css = cddb;
        Cds = -(cgdb + cddb + cbdb);
        dsxpart_dVs = (Css - sxpart * (Css + Cds)) / qcheq;
        Csg = cdgb;
        Cdg = -(cggb + cdgb + cbgb);
        dsxpart_dVg = (Csg - sxpart * (Csg + Cdg)) / qcheq;

        Csd = cdsb;
        Cdd = -(cgsb + cdsb + cbsb);
        dsxpart_dVd = (Csd - sxpart * (Csd + Cdd)) / qcheq;

        dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);
      }

      dxpart = 1.0 - sxpart;
      ddxpart_dVd = -dsxpart_dVd;
      ddxpart_dVg = -dsxpart_dVg;
      ddxpart_dVs = -dsxpart_dVs;
      ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);

      qgd = qgdo;
      qgs = qgso;
      qgb = sizeDepParams.cgbo * vgb;
      qgate = qgd + qgs + qgb;
      qbulk = -qgb;
      qsrc = -qgs;
      qdrn = -(qgate + qbulk + qsrc);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Purpose       : setupCapacitors_newDAE
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/29/2018
//-----------------------------------------------------------------------------
template <typename ScalarT> 
bool setupCapacitors_newDAE (
    int mode,
    int nqsMod,
    const SizeDependParam<ScalarT> & sizeDepParams,
    const ScalarT & cggb,
    const ScalarT & cgdo,
    const ScalarT & cgso,
    const ScalarT & cgdb,
    const ScalarT & cgsb,
    const ScalarT & cdgb,
    const ScalarT & cddb,
    const ScalarT & capbd,
    const ScalarT & cdsb,
    const ScalarT & cbgb,
    const ScalarT & cbdb,
    const ScalarT & capbs,
    const ScalarT & cbsb,
    //outputs
    ScalarT & CAPcggb,
    ScalarT & CAPcgdb,
    ScalarT & CAPcgsb,
    ScalarT & CAPcdgb,
    ScalarT & CAPcddb,
    ScalarT & CAPcdsb,
    ScalarT & CAPcsgb,
    ScalarT & CAPcsdb,
    ScalarT & CAPcssb,
    ScalarT & CAPcbgb,
    ScalarT & CAPcbdb,
    ScalarT & CAPcbsb
    )
{
  if (mode > 0)
  {
    if (nqsMod == 0)
    {
      CAPcggb = (cggb + cgdo + cgso + sizeDepParams.cgbo );
      CAPcgdb = (cgdb - cgdo);
      CAPcgsb = (cgsb - cgso);

      CAPcdgb = (cdgb - cgdo);
      CAPcddb = (cddb + capbd + cgdo);
      CAPcdsb = cdsb;

      CAPcsgb = -(cggb + cbgb + cdgb + cgso);
      CAPcsdb = -(cgdb + cbdb + cddb);
      CAPcssb = (capbs + cgso - (cgsb + cbsb + cdsb));

      CAPcbgb = (cbgb - sizeDepParams.cgbo);
      CAPcbdb = (cbdb - capbd);
      CAPcbsb = (cbsb - capbs);
    }
    else  // nqsMode != 0
    {

    } // nqsMod
  }
  else
  {
    if (nqsMod == 0)
    {
      CAPcggb = (cggb + cgdo + cgso + sizeDepParams.cgbo );
      CAPcgdb = (cgsb - cgdo);
      CAPcgsb = (cgdb - cgso);

      CAPcdgb = -(cggb + cbgb + cdgb + cgdo);
      CAPcddb = (capbd + cgdo - (cgsb + cbsb + cdsb));
      CAPcdsb = -(cgdb + cbdb + cddb);

      CAPcsgb = (cdgb - cgso);
      CAPcsdb = cdsb;
      CAPcssb = (cddb + capbs + cgso);

      CAPcbgb = (cbgb - sizeDepParams.cgbo);
      CAPcbdb = (cbsb - capbd);
      CAPcbsb = (cbdb - capbs);
    }
    else  // nqsMode != 0
    {

    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/29/2018
//-----------------------------------------------------------------------------
template <typename ScalarT> 
bool updateIntermediateVars (
  const SolverState & solverState,
  const DeviceOptions & deviceOptions,
  const SizeDependParam<ScalarT> & sizeDepParams,

  int model_dtype,
  int model_mobMod,
  int model_capMod,

  const ScalarT & vbs,
  const ScalarT & vds,
  const ScalarT & vgs,
  const ScalarT & vtm,
  const ScalarT & model_jctEmissionCoeff,
  const ScalarT & model_ijth,
  const ScalarT & model_factor1,
  const ScalarT & model_tnom,
  const ScalarT & model_tox,
  const ScalarT & model_cox,
  const ScalarT & model_xpart,
  const ScalarT & model_bulkJctBotGradingCoeff,
  const ScalarT & model_bulkJctSideGradingCoeff,
  const ScalarT & model_bulkJctGateSideGradingCoeff,

  const ScalarT & sourceArea,
  const ScalarT & sourcePerimeter,
  const ScalarT & drainArea,
  const ScalarT & drainPerimeter,
  const ScalarT & jctTempSatCurDensity,
  const ScalarT & jctSidewallTempSatCurDensity,
  const ScalarT & vjsm,
  const ScalarT & IsEvjsm,
  const ScalarT & vjdm,
  const ScalarT & IsEvjdm,

  int mode,

  ScalarT &  T1global,
  ScalarT &  thetavth,
  const ScalarT &  temp,
  ScalarT &  Vth,
  ScalarT &  von,
  ScalarT &  dVgs_eff_dVg,
  ScalarT &  Vgsteff,
  ScalarT &  Abulk,
  ScalarT &  ueff,
  ScalarT &  Vdsat,
  ScalarT &  vdsat,
  ScalarT &  Vdseff,
  ScalarT &  Gm,
  ScalarT &  cdrain,
  ScalarT &  gds,
  ScalarT &  gm,

  ScalarT &  gmbs,
  ScalarT &  gbbs,
  ScalarT &  gbgs,
  ScalarT &  gbds,
  ScalarT &  csub,
  ScalarT &  qgate,
  ScalarT &  qdrn,
  ScalarT &  qsrc,
  ScalarT &  qbulk,
  ScalarT &  cggb,
  ScalarT &  cgsb,
  ScalarT &  cgdb,
  ScalarT &  cdgb,
  ScalarT &  cdsb,
  ScalarT &  cddb,
  ScalarT &  cbgb,
  ScalarT &  cbsb,
  ScalarT &  cbdb,
  ScalarT &  cqdb,
  ScalarT &  cqsb,
  ScalarT &  cqgb,
  ScalarT &  cqbb,
  ScalarT &  gtau,

  ScalarT &  dVgst_dVb,
  ScalarT &  dVgst_dVg,
  ScalarT &  CoxWL,
  ScalarT &  qinv,
  ScalarT &  Cgg,
  ScalarT &  Cgd,
  ScalarT &  Cgb,
  ScalarT &  Cbg,
  ScalarT &  Cbd,
  ScalarT &  Cbb,
  ScalarT &  Csg,
  ScalarT &  Csb,
  ScalarT &  Csd,
  ScalarT &  dDeltaPhi_dVg,
  ScalarT &  dDeltaPhi_dVd,
  ScalarT &  dDeltaPhi_dVb,
  ScalarT &  cd,

  const ScalarT &  unitAreaJctCapTemp,
  const ScalarT &  unitLengthGateSidewallJctCapTemp,
  const ScalarT &  unitLengthSidewallJctCapTemp,

  ScalarT &  qbs,
  ScalarT &  capbs,

  const ScalarT &  PhiBTemp,
  const ScalarT &  PhiBSWTemp,
  const ScalarT &  PhiBSWGTemp,

  ScalarT &  qbd,
  ScalarT &  capbd,
  int nqsMod,
  ScalarT &  qcheq,

  const ScalarT &  model_vtm,

  ScalarT &  cgdo,
  ScalarT &  qgdo,
  ScalarT &  cgso,
  ScalarT &  qgso,
  ScalarT &  cqdef,
  ScalarT &  qg,
  ScalarT &  qd,
  ScalarT &  qb,
  ScalarT &  qcdump,

  const ScalarT &  qdef,
  ScalarT &  Idrain,
  const ScalarT &  drainConductance,
  const ScalarT &  Vddp,
  ScalarT &  Isource,
  const ScalarT &  sourceConductance,
  const ScalarT &  Vssp,

  //outputs
  ScalarT & gbs,
  ScalarT & cbs,
  ScalarT & gbd,
  ScalarT & cbd,

  ScalarT & CAPcggb,
  ScalarT & CAPcgdb,
  ScalarT & CAPcgsb,
  ScalarT & CAPcdgb,
  ScalarT & CAPcddb,
  ScalarT & CAPcdsb,
  ScalarT & CAPcsgb,
  ScalarT & CAPcsdb,
  ScalarT & CAPcssb,
  ScalarT & CAPcbgb,
  ScalarT & CAPcbdb,
  ScalarT & CAPcbsb,
  bool & ChargeComputationNeeded
    )
{
  bool bsuccess = true;

  // begin the b3ld.c parameters:
  ScalarT SourceSatCurrent(0.0), DrainSatCurrent(0.0);
  ScalarT vgdo(0.0);

  ScalarT VgstNVt(0.0), ExpVgst(0.0);

  ScalarT czbd(0.0), czbdsw(0.0), czbdswg(0.0), czbs(0.0), czbssw(0.0), czbsswg(0.0);
  ScalarT evbd(0.0), evbs(0.0), arg(0.0), sarg(0.0);

  ScalarT Vfbeff(0.0), dVfbeff_dVg(0.0), dVfbeff_dVb(0.0), V3(0.0), V4(0.0);

  ScalarT MJ(0.0), MJSW(0.0), MJSWG(0.0);

  ScalarT qinoi(0.0);

  ScalarT Vds(0.0);
  ScalarT Vgs(0.0), Vbs(0.0);

  ScalarT Vgs_eff(0.0), Vfb(0.0);
  ScalarT Phis(0.0), dPhis_dVb(0.0), sqrtPhis(0.0), dsqrtPhis_dVb(0.0);
  //ScalarT Vth(0.0);  // made into instance variable
  ScalarT dVth_dVb(0.0), dVth_dVd(0.0);
  ScalarT Vgst(0.0);

  ScalarT Nvtm(0.0);
  ScalarT Vtm(0.0);
  ScalarT n(0.0), dn_dVb(0.0), dn_dVd(0.0), voffcv(0.0), noff(0.0), dnoff_dVd(0.0), dnoff_dVb(0.0);
  ScalarT ExpArg(0.0), V0(0.0), CoxWLcen(0.0), QovCox(0.0), LINK(0.0);
  ScalarT DeltaPhi(0.0);

  ScalarT Cox(0.0), Tox(0.0), Tcen(0.0), dTcen_dVg(0.0), dTcen_dVd(0.0), dTcen_dVb(0.0);
  ScalarT Ccen(0.0), Coxeff(0.0), dCoxeff_dVg(0.0), dCoxeff_dVd(0.0), dCoxeff_dVb(0.0);
  ScalarT Denomi(0.0), dDenomi_dVg(0.0), dDenomi_dVd(0.0), dDenomi_dVb(0.0);

  ScalarT dueff_dVg(0.0), dueff_dVd(0.0), dueff_dVb(0.0);
  ScalarT Esat(0.0);

  //ScalarT Vdsat(0.0); // made into instance variable

  ScalarT EsatL(0.0), dEsatL_dVg(0.0), dEsatL_dVd(0.0), dEsatL_dVb(0.0);

  ScalarT dVdsat_dVg(0.0), dVdsat_dVb(0.0), dVdsat_dVd(0.0), Vasat(0.0), dAlphaz_dVg(0.0), dAlphaz_dVb(0.0);
  ScalarT dVasat_dVg(0.0), dVasat_dVb(0.0), dVasat_dVd(0.0), Va(0.0);

  ScalarT dVa_dVd(0.0), dVa_dVg(0.0), dVa_dVb(0.0);
  ScalarT Vbseff(0.0), dVbseff_dVb(0.0), VbseffCV(0.0), dVbseffCV_dVb(0.0);
  ScalarT Arg1(0.0);

  ScalarT One_Third_CoxWL(0.0), Two_Third_CoxWL(0.0), Alphaz(0.0);

  ScalarT T0(0.0), dT0_dVg(0.0), dT0_dVd(0.0), dT0_dVb(0.0);
  ScalarT T1(0.0), dT1_dVg(0.0), dT1_dVd(0.0), dT1_dVb(0.0);
  ScalarT T2(0.0), dT2_dVg(0.0), dT2_dVd(0.0), dT2_dVb(0.0);
  ScalarT T3(0.0), dT3_dVg(0.0), dT3_dVd(0.0), dT3_dVb(0.0);
  ScalarT T4(0.0);

  ScalarT T5(0.0);
  ScalarT T6(0.0);
  ScalarT T7(0.0);
  ScalarT T8(0.0);
  ScalarT T9(0.0);
  ScalarT T10(0.0);
  ScalarT T11(0.0), T12(0.0);

  ScalarT tmp(0.0); 

  //ScalarT Abulk(0.0);  // needs to be instance var for noise

  ScalarT dAbulk_dVb(0.0), Abulk0(0.0), dAbulk0_dVb(0.0);

  ScalarT VACLM(0.0), dVACLM_dVg(0.0), dVACLM_dVd(0.0), dVACLM_dVb(0.0);
  ScalarT VADIBL(0.0), dVADIBL_dVg(0.0), dVADIBL_dVd(0.0), dVADIBL_dVb(0.0);

  ScalarT Xdep(0.0), dXdep_dVb(0.0), lt1(0.0), dlt1_dVb(0.0), ltw(0.0), dltw_dVb(0.0);
  ScalarT Delt_vth(0.0), dDelt_vth_dVb(0.0);

  ScalarT Theta0(0.0), dTheta0_dVb(0.0);

  ScalarT TempRatio(0.0), tmp1(0.0), tmp2(0.0), tmp3(0.0), tmp4(0.0);

  ScalarT DIBL_Sft(0.0), dDIBL_Sft_dVd(0.0);

  ScalarT Lambda(0.0), dLambda_dVg(0.0);

  ScalarT a1(0.0);

  //ScalarT Vgsteff(0.0);   // needs to be an instance var for noise.

  ScalarT dVgsteff_dVg(0.0), dVgsteff_dVd(0.0), dVgsteff_dVb(0.0);
  //ScalarT Vdseff(0.0);  // needs to be an instance var for noise.
  ScalarT dVdseff_dVg(0.0), dVdseff_dVd(0.0), dVdseff_dVb(0.0);
  ScalarT VdseffCV(0.0), dVdseffCV_dVg(0.0), dVdseffCV_dVd(0.0), dVdseffCV_dVb(0.0);
  ScalarT diffVds(0.0);

  ScalarT dAbulk_dVg(0.0);
  ScalarT beta(0.0), dbeta_dVg(0.0), dbeta_dVd(0.0), dbeta_dVb(0.0);
  ScalarT gche(0.0), dgche_dVg(0.0), dgche_dVd(0.0), dgche_dVb(0.0);
  ScalarT fgche1(0.0), dfgche1_dVg(0.0), dfgche1_dVd(0.0), dfgche1_dVb(0.0);
  ScalarT fgche2(0.0), dfgche2_dVg(0.0), dfgche2_dVd(0.0), dfgche2_dVb(0.0);
  ScalarT Idl(0.0), dIdl_dVg(0.0), dIdl_dVd(0.0), dIdl_dVb(0.0);
  ScalarT Idsa(0.0), dIdsa_dVg(0.0), dIdsa_dVd(0.0), dIdsa_dVb(0.0);
  ScalarT Ids(0.0);

  ScalarT Gds(0.0), Gmb(0.0);
  ScalarT Isub(0.0);

  ScalarT Gbd(0.0), Gbg(0.0), Gbb(0.0);
  ScalarT VASCBE(0.0), dVASCBE_dVg(0.0), dVASCBE_dVd(0.0), dVASCBE_dVb(0.0);
  ScalarT CoxWovL(0.0);
  ScalarT Rds(0.0), dRds_dVg(0.0), dRds_dVb(0.0), WVCox(0.0), WVCoxRds(0.0);
  ScalarT Vgst2Vtm(0.0), VdsatCV(0.0);

  ScalarT dVdsatCV_dVg(0.0), dVdsatCV_dVb(0.0);
  ScalarT Leff(0.0), Weff(0.0), dWeff_dVg(0.0), dWeff_dVb(0.0);
  ScalarT AbulkCV(0.0), dAbulkCV_dVb(0.0);

  ScalarT gtau_diff(0.0), gtau_drift(0.0);
  // these shadow member varables and then are unitialized
  // when used in later calculations
  // ScalarT qcheq(0.0), cqcheq(0.0), qdef(0.0);

  ScalarT Cgg1(0.0), Cgb1(0.0), Cgd1(0.0), Cbg1(0.0), Cbb1(0.0), Cbd1(0.0);

  ScalarT Qac0(0.0), Qsub0(0.0);
  ScalarT dQac0_dVg(0.0), dQac0_dVb(0.0), dQsub0_dVg(0.0), dQsub0_dVd(0.0), dQsub0_dVb(0.0);
  ScalarT von_local(0.0);

  ScalarT ScalingFactor = 1.0e-9;

  //bool ChargeComputationNeeded = true;

  // Don't do charge computations in DC sweeps.
  if (solverState.tranopFlag || solverState.acopFlag || solverState.transientFlag)
  {
    ChargeComputationNeeded = true;
  }
  else
  {
    ChargeComputationNeeded = false;
  }

  int Check = 1;

// ------------------
// get the voltage drops ...

// ------------------

  // determine DC current and derivatives
  ScalarT vbd = vbs - vds;
  ScalarT vgd = vgs - vds;
  ScalarT vgb = vgs - vbs;

  // Source/drain junction diode DC model begins
  Nvtm = vtm * model_jctEmissionCoeff;

  if ((sourceArea <= 0.0) && (sourcePerimeter <= 0.0))
  {
    SourceSatCurrent = 1.0e-14;
  }
  else
  {
    SourceSatCurrent = sourceArea
      * jctTempSatCurDensity
      + sourcePerimeter
      * jctSidewallTempSatCurDensity;
  }

  if (SourceSatCurrent <= 0.0)
  {
    gbs = deviceOptions.gmin;
    cbs = gbs * vbs;
  }
  else
  {
    if (model_ijth == 0.0)
    {
      evbs = exp(vbs / Nvtm);
      gbs = SourceSatCurrent * evbs / Nvtm + deviceOptions.gmin;
      cbs = SourceSatCurrent * (evbs - 1.0) + deviceOptions.gmin * vbs;
    }
    else
    {
      if (vbs < vjsm)
      {
        evbs = exp(vbs / Nvtm);
        gbs = SourceSatCurrent * evbs / Nvtm + deviceOptions.gmin;
        cbs = SourceSatCurrent * (evbs - 1.0) + deviceOptions.gmin * vbs;
      }
      else
      {
        T0 = IsEvjsm / Nvtm;
        gbs = T0 + deviceOptions.gmin;
        cbs = IsEvjsm - SourceSatCurrent
          + T0 * (vbs - vjsm)
          + deviceOptions.gmin * vbs;
      }
    }
  }

  if ((drainArea <= 0.0) && (drainPerimeter <= 0.0))
  {
    DrainSatCurrent = 1.0e-14;
  }
  else
  {
    DrainSatCurrent = drainArea
      * jctTempSatCurDensity
      + drainPerimeter
      * jctSidewallTempSatCurDensity;
  }

  if (DrainSatCurrent <= 0.0)
  {
    gbd = deviceOptions.gmin;
    cbd = gbd * vbd;
  }
  else
  {
    if (model_ijth == 0.0)
    {
      evbd = exp(vbd / Nvtm);
      gbd = DrainSatCurrent * evbd / Nvtm + deviceOptions.gmin;
      cbd = DrainSatCurrent * (evbd - 1.0) + deviceOptions.gmin * vbd;
    }
    else
    {
      if (vbd < vjdm)
      {
        evbd = exp(vbd / Nvtm);
        gbd = DrainSatCurrent * evbd / Nvtm + deviceOptions.gmin;
        cbd = DrainSatCurrent * (evbd - 1.0) + deviceOptions.gmin * vbd;
      }
      else
      {
        T0 = IsEvjdm / Nvtm;
        gbd = T0 + deviceOptions.gmin;
        cbd = IsEvjdm - DrainSatCurrent
          + T0 * (vbd - vjdm)
          + deviceOptions.gmin * vbd;
      }
    }
  }
  // End of diode DC model

  if (vds >= 0.0)
  {   // normal mode
    mode = 1;
    Vds = vds;
    Vgs = vgs;
    Vbs = vbs;
  }
  else
  {   // inverse mode
    mode = -1;
    Vds = -vds;
    Vgs = vgd;
    Vbs = vbd;
  }

  // exclude the continuation stuff.  hmmm.    This seems correct for traditional sensitivities, but 
  // not sure if it will work with arclength continuation ...  check later.

  T0 = Vbs - sizeDepParams.vbsc - 0.001;
  T1global = sqrt(T0 * T0 - 0.004 * sizeDepParams.vbsc);
  Vbseff = sizeDepParams.vbsc + 0.5 * (T0 + T1global);
  dVbseff_dVb = 0.5 * (1.0 + T0 / T1global);

  if (Vbseff < Vbs) Vbseff = Vbs;

  if (Vbseff > 0.0)
  {
    T0 = sizeDepParams.phi / (sizeDepParams.phi + Vbseff);
    Phis = sizeDepParams.phi * T0;
    dPhis_dVb = -T0 * T0;
    sqrtPhis = sizeDepParams.phis3 / (sizeDepParams.phi + 0.5 * Vbseff);
    dsqrtPhis_dVb = -0.5 * sqrtPhis * sqrtPhis / sizeDepParams.phis3;
  }
  else
  {
    Phis = sizeDepParams.phi - Vbseff;
    dPhis_dVb = -1.0;
    sqrtPhis = sqrt(Phis);
    dsqrtPhis_dVb = -0.5 / sqrtPhis;
  }

  Xdep = sizeDepParams.Xdep0 * sqrtPhis / sizeDepParams.sqrtPhi;
  dXdep_dVb = (sizeDepParams.Xdep0 / sizeDepParams.sqrtPhi) * dsqrtPhis_dVb;

  Leff = sizeDepParams.leff;
  Vtm = vtm;
  // Vth Calculation
  T3 = sqrt(Xdep);
  V0 = sizeDepParams.vbi - sizeDepParams.phi;

  T0 = sizeDepParams.dvt2 * Vbseff;
  if (T0 >= - 0.5)
  {
    T1 = 1.0 + T0;
    T2 = sizeDepParams.dvt2;
  }
  else // Added to avoid any discontinuity problems caused by dvt2
  {
    T4 = 1.0 / (3.0 + 8.0 * T0);
    T1 = (1.0 + 3.0 * T0) * T4;
    T2 = sizeDepParams.dvt2 * T4 * T4;
  }

  lt1 = model_factor1 * T3 * T1;
  dlt1_dVb = model_factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

  T0 = sizeDepParams.dvt2w * Vbseff;
  if (T0 >= - 0.5)
  {
    T1 = 1.0 + T0;
    T2 = sizeDepParams.dvt2w;
  }
  else // Added to avoid any discontinuity problems caused by dvt2w
  {
    T4 = 1.0 / (3.0 + 8.0 * T0);
    T1 = (1.0 + 3.0 * T0) * T4;
    T2 = sizeDepParams.dvt2w * T4 * T4;
  }

  ltw = model_factor1 * T3 * T1;
  dltw_dVb = model_factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

  T0 = -0.5 * sizeDepParams.dvt1 * Leff / lt1;
  if (T0 > -CONSTEXP_THRESHOLD)
  {
    T1 = exp(T0);
    Theta0 = T1 * (1.0 + 2.0 * T1);
    dT1_dVb = -T0 / lt1 * T1 * dlt1_dVb;
    dTheta0_dVb = (1.0 + 4.0 * T1) * dT1_dVb;
  }
  else
  {
    T1 = CONSTMIN_EXP;
    Theta0 = T1 * (1.0 + 2.0 * T1);
    dTheta0_dVb = 0.0;
  }

  thetavth = sizeDepParams.dvt0 * Theta0;
  Delt_vth = thetavth * V0;
  dDelt_vth_dVb = sizeDepParams.dvt0 * dTheta0_dVb * V0;

  T0 = -0.5 * sizeDepParams.dvt1w * sizeDepParams.weff * Leff / ltw;
  if (T0 > -CONSTEXP_THRESHOLD)
  {
    T1 = exp(T0);
    T2 = T1 * (1.0 + 2.0 * T1);
    dT1_dVb = -T0 / ltw * T1 * dltw_dVb;
    dT2_dVb = (1.0 + 4.0 * T1) * dT1_dVb;
  }
  else
  {
    T1 = CONSTMIN_EXP;
    T2 = T1 * (1.0 + 2.0 * T1);
    dT2_dVb = 0.0;
  }

  T0 = sizeDepParams.dvt0w * T2;
  T2 = T0 * V0;
  dT2_dVb = sizeDepParams.dvt0w * dT2_dVb * V0;

  TempRatio =  temp / model_tnom - 1.0;
  T0 = sqrt(1.0 + sizeDepParams.nlx / Leff);
  T1 = sizeDepParams.k1ox * (T0 - 1.0) * sizeDepParams.sqrtPhi
    + (sizeDepParams.kt1 + sizeDepParams.kt1l / Leff
        +  sizeDepParams.kt2 * Vbseff) * TempRatio;

  tmp2 = model_tox * sizeDepParams.phi / (sizeDepParams.weff + sizeDepParams.w0);

  T3 = sizeDepParams.eta0 + sizeDepParams.etab * Vbseff;
  if (T3 < 1.0e-4) // avoid  discontinuity problems caused by etab
  {
    T9 = 1.0 / (3.0 - 2.0e4 * T3);
    T3 = (2.0e-4 - T3) * T9;
    T4 = T9 * T9;
  }
  else
  {
    T4 = 1.0;
  }

  dDIBL_Sft_dVd = T3 * sizeDepParams.theta0vb0;
  DIBL_Sft = dDIBL_Sft_dVd * Vds;

  Vth = model_dtype * sizeDepParams.vth0 - sizeDepParams.k1
    * sizeDepParams.sqrtPhi + sizeDepParams.k1ox * sqrtPhis
    - sizeDepParams.k2ox * Vbseff - Delt_vth - T2 + (sizeDepParams.k3
        + sizeDepParams.k3b * Vbseff) * tmp2 + T1 - DIBL_Sft;

  von = Vth;

  dVth_dVb = sizeDepParams.k1ox * dsqrtPhis_dVb - sizeDepParams.k2ox
    - dDelt_vth_dVb - dT2_dVb + sizeDepParams.k3b * tmp2
    - sizeDepParams.etab * Vds * sizeDepParams.theta0vb0 * T4
    + sizeDepParams.kt2 * TempRatio;

  dVth_dVd = -dDIBL_Sft_dVd;

  // Calculate n
  tmp2 = sizeDepParams.nfactor * CONSTEPSSI / Xdep;
  tmp3 = sizeDepParams.cdsc + sizeDepParams.cdscb * Vbseff
    + sizeDepParams.cdscd * Vds;
  tmp4 = (tmp2 + tmp3 * Theta0 + sizeDepParams.cit) / model_cox;

  if (tmp4 >= -0.5)
  {
    n = 1.0 + tmp4;
    dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
        + sizeDepParams.cdscb * Theta0) / model_cox;
    dn_dVd = sizeDepParams.cdscd * Theta0 / model_cox;
  }
  else // avoid  discontinuity problems caused by tmp4
  {
    T0 = 1.0 / (3.0 + 8.0 * tmp4);
    n = (1.0 + 3.0 * tmp4) * T0;
    T0 *= T0;
    dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
        + sizeDepParams.cdscb * Theta0) / model_cox * T0;
    dn_dVd = sizeDepParams.cdscd * Theta0 / model_cox * T0;
  }

  // Poly Gate Si Depletion Effect
  T0 = sizeDepParams.vfb + sizeDepParams.phi;

  // added to avoid the problem caused by ngate
  if ((sizeDepParams.ngate > 1.e18) && (sizeDepParams.ngate < 1.e25) && (Vgs > T0))
  {
    T1 = 1.0e6 * CONSTQ * CONSTEPSSI * sizeDepParams.ngate
      / (model_cox * model_cox);
    T4 = sqrt(1.0 + 2.0 * (Vgs - T0) / T1);

    T2 = T1 * (T4 - 1.0);
    T3 = 0.5 * T2 * T2 / T1; // T3 = Vpoly
    T7 = 1.12 - T3 - 0.05;
    T6 = sqrt(T7 * T7 + 0.224);
    T5 = 1.12 - 0.5 * (T7 + T6);
    Vgs_eff = Vgs - T5;
    dVgs_eff_dVg = 1.0 - (0.5 - 0.5 / T4) * (1.0 + T7 / T6);
  }
  else
  {
    Vgs_eff = Vgs;
    dVgs_eff_dVg = 1.0;
  }
  Vgst = Vgs_eff - Vth;

  // Effective Vgst (Vgsteff) Calculation
  T10 = 2.0 * n * Vtm;
  VgstNVt = Vgst / T10;
  ExpArg = (2.0 * sizeDepParams.voff - Vgst) / T10;

  // MCJ: Very small Vgst
  if (VgstNVt > CONSTEXP_THRESHOLD)
  {
    Vgsteff = Vgst;
    dVgsteff_dVg = dVgs_eff_dVg;
    dVgsteff_dVd = -dVth_dVd;
    dVgsteff_dVb = -dVth_dVb;
  }
  else if (ExpArg > CONSTEXP_THRESHOLD)
  {
    T0 = (Vgst - sizeDepParams.voff) / (n * Vtm);
    ExpVgst = exp(T0);
    Vgsteff = Vtm * sizeDepParams.cdep0 / model_cox * ExpVgst;
    dVgsteff_dVg = Vgsteff / (n * Vtm);
    dVgsteff_dVd = -dVgsteff_dVg * (dVth_dVd + T0 * Vtm * dn_dVd);
    dVgsteff_dVb = -dVgsteff_dVg * (dVth_dVb + T0 * Vtm * dn_dVb);
    dVgsteff_dVg *= dVgs_eff_dVg;
  }
  else
  {
    ExpVgst = exp(VgstNVt);
    T1 = T10 * log(1.0 + ExpVgst);
    dT1_dVg = ExpVgst / (1.0 + ExpVgst);
    dT1_dVb = -dT1_dVg * (dVth_dVb + Vgst / n * dn_dVb) + T1 / n * dn_dVb;
    dT1_dVd = -dT1_dVg * (dVth_dVd + Vgst / n * dn_dVd) + T1 / n * dn_dVd;

    dT2_dVg = -model_cox / (Vtm * sizeDepParams.cdep0) * exp(ExpArg);
    T2 = 1.0 - T10 * dT2_dVg;

    dT2_dVd = -dT2_dVg * (dVth_dVd - 2.0 * Vtm * ExpArg * dn_dVd)
      + (T2 - 1.0) / n * dn_dVd;

    dT2_dVb = -dT2_dVg * (dVth_dVb - 2.0 * Vtm * ExpArg * dn_dVb)
      + (T2 - 1.0) / n * dn_dVb;

    Vgsteff = T1 / T2;
    T3 = T2 * T2;

    dVgsteff_dVg = (T2 * dT1_dVg - T1 * dT2_dVg) / T3 * dVgs_eff_dVg;
    dVgsteff_dVd = (T2 * dT1_dVd - T1 * dT2_dVd) / T3;
    dVgsteff_dVb = (T2 * dT1_dVb - T1 * dT2_dVb) / T3;
  }

  // Calculate Effective Channel Geometry
  T9 = sqrtPhis - sizeDepParams.sqrtPhi;
  Weff = sizeDepParams.weff -2.0 *(sizeDepParams.dwg * Vgsteff + sizeDepParams.dwb * T9);
  dWeff_dVg = -2.0 * sizeDepParams.dwg;
  dWeff_dVb = -2.0 * sizeDepParams.dwb * dsqrtPhis_dVb;

  if (Weff < 2.0e-8) // to avoid the discontinuity problem due to Weff
  {
    T0 = 1.0 / (6.0e-8 - 2.0 * Weff);
    Weff = 2.0e-8 * (4.0e-8 - Weff) * T0;
    T0 *= T0 * 4.0e-16;
    dWeff_dVg *= T0;
    dWeff_dVb *= T0;
  }

  T0 = sizeDepParams.prwg * Vgsteff + sizeDepParams.prwb * T9;
  if (T0 >= -0.9)
  {
    Rds = sizeDepParams.rds0 * (1.0 + T0);
    dRds_dVg = sizeDepParams.rds0 * sizeDepParams.prwg;
    dRds_dVb = sizeDepParams.rds0 * sizeDepParams.prwb * dsqrtPhis_dVb;
  }
  else // to avoid the discontinuity problem due to prwg and prwb
  {
    T1 = 1.0 / (17.0 + 20.0 * T0);
    Rds = sizeDepParams.rds0 * (0.8 + T0) * T1;
    T1 *= T1;
    dRds_dVg = sizeDepParams.rds0 * sizeDepParams.prwg * T1;
    dRds_dVb = sizeDepParams.rds0 * sizeDepParams.prwb * dsqrtPhis_dVb * T1;
  }

  // Calculate Abulk
  T1 = 0.5 * sizeDepParams.k1ox / sqrtPhis;
  dT1_dVb = -T1 / sqrtPhis * dsqrtPhis_dVb;

  T9 = sqrt(sizeDepParams.xj * Xdep);
  tmp1 = Leff + 2.0 * T9;
  T5 = Leff / tmp1;
  tmp2 = sizeDepParams.a0 * T5;
  tmp3 = sizeDepParams.weff + sizeDepParams.b1;
  tmp4 = sizeDepParams.b0 / tmp3;
  T2 = tmp2 + tmp4;
  dT2_dVb = -T9 / tmp1 / Xdep * dXdep_dVb;
  T6 = T5 * T5;
  T7 = T5 * T6;

  Abulk0 = 1.0 + T1 * T2;
  dAbulk0_dVb = T1 * tmp2 * dT2_dVb + T2 * dT1_dVb;

  T8 = sizeDepParams.ags * sizeDepParams.a0 * T7;
  dAbulk_dVg = -T1 * T8;
  Abulk = Abulk0 + dAbulk_dVg * Vgsteff;
  dAbulk_dVb = dAbulk0_dVb - T8 * Vgsteff * (dT1_dVb + 3.0 * T1 * dT2_dVb);

  if (Abulk0 < 0.1) // added to avoid the problems caused by Abulk0
  {
    T9 = 1.0 / (3.0 - 20.0 * Abulk0);
    Abulk0 = (0.2 - Abulk0) * T9;
    dAbulk0_dVb *= T9 * T9;
  }

  if (Abulk < 0.1) // added to avoid the problems caused by Abulk
  {
    T9 = 1.0 / (3.0 - 20.0 * Abulk);
    Abulk = (0.2 - Abulk) * T9;
    T10 = T9 * T9;
    dAbulk_dVb *= T10;
    dAbulk_dVg *= T10;
  }

  T2 = sizeDepParams.keta * Vbseff;
  if (T2 >= -0.9)
  {
    T0 = 1.0 / (1.0 + T2);
    dT0_dVb = -sizeDepParams.keta * T0 * T0;
  }
  else // added to avoid the problems caused by Keta
  {
    T1 = 1.0 / (0.8 + T2);
    T0 = (17.0 + 20.0 * T2) * T1;
    dT0_dVb = -sizeDepParams.keta * T1 * T1;
  }

  dAbulk_dVg *= T0;
  dAbulk_dVb = dAbulk_dVb * T0 + Abulk * dT0_dVb;
  dAbulk0_dVb = dAbulk0_dVb * T0 + Abulk0 * dT0_dVb;
  Abulk *= T0;
  Abulk0 *= T0;

  // Mobility calculation
  if (model_mobMod == 1)
  {
    T0 = Vgsteff + Vth + Vth;
    T2 = sizeDepParams.ua + sizeDepParams.uc * Vbseff;
    T3 = T0 / model_tox;
    T5 = T3 * (T2 + sizeDepParams.ub * T3);
    dDenomi_dVg = (T2 + 2.0 * sizeDepParams.ub * T3) / model_tox;
    dDenomi_dVd = dDenomi_dVg * 2.0 * dVth_dVd;
    dDenomi_dVb = dDenomi_dVg * 2.0 * dVth_dVb + sizeDepParams.uc * T3;
  }
  else if (model_mobMod == 2)
  {
    T5 = Vgsteff / model_tox * (sizeDepParams.ua
        + sizeDepParams.uc * Vbseff + sizeDepParams.ub * Vgsteff / model_tox);

    dDenomi_dVg = (sizeDepParams.ua + sizeDepParams.uc * Vbseff
        + 2.0 * sizeDepParams.ub * Vgsteff / model_tox) / model_tox;

    dDenomi_dVd = 0.0;
    dDenomi_dVb = Vgsteff * sizeDepParams.uc / model_tox;
  }
  else
  {
    T0 = Vgsteff + Vth + Vth;
    T2 = 1.0 + sizeDepParams.uc * Vbseff;
    T3 = T0 / model_tox;
    T4 = T3 * (sizeDepParams.ua + sizeDepParams.ub * T3);
    T5 = T4 * T2;

    dDenomi_dVg = (sizeDepParams.ua + 2.0 * sizeDepParams.ub * T3) * T2 /model_tox;
    dDenomi_dVd = dDenomi_dVg * 2.0 * dVth_dVd;
    dDenomi_dVb = dDenomi_dVg * 2.0 * dVth_dVb + sizeDepParams.uc * T4;
  }

  if (T5 >= -0.8)
  {
    Denomi = 1.0 + T5;
  }
  else // Added to avoid the discontinuity problem caused by ua and ub
  {
    T9 = 1.0 / (7.0 + 10.0 * T5);
    Denomi = (0.6 + T5) * T9;
    T9 *= T9;
    dDenomi_dVg *= T9;
    dDenomi_dVd *= T9;
    dDenomi_dVb *= T9;
  }

  ueff = sizeDepParams.u0temp / Denomi;
  T9 = -ueff / Denomi;
  dueff_dVg = T9 * dDenomi_dVg;
  dueff_dVd = T9 * dDenomi_dVd;
  dueff_dVb = T9 * dDenomi_dVb;

  // Saturation Drain Voltage  Vdsat
  WVCox = Weff * sizeDepParams.vsattemp * model_cox;
  WVCoxRds = WVCox * Rds;

  Esat = 2.0 * sizeDepParams.vsattemp / ueff;
  EsatL = Esat * Leff;
  T0 = -EsatL /ueff;
  dEsatL_dVg = T0 * dueff_dVg;
  dEsatL_dVd = T0 * dueff_dVd;
  dEsatL_dVb = T0 * dueff_dVb;

  // Sqrt()
  a1 = sizeDepParams.a1;
  if (a1 == 0.0)
  {
    Lambda = sizeDepParams.a2;
    dLambda_dVg = 0.0;
  }
  else if (a1 > 0.0) // Added to avoid the discontinuity problem
    // caused by a1 and a2 (Lambda)
  {
    T0 = 1.0 - sizeDepParams.a2;
    T1 = T0 - sizeDepParams.a1 * Vgsteff - 0.0001;
    T2 = sqrt(T1 * T1 + 0.0004 * T0);
    Lambda = sizeDepParams.a2 + T0 - 0.5 * (T1 + T2);
    dLambda_dVg = 0.5 * sizeDepParams.a1 * (1.0 + T1 / T2);
  }
  else
  {
    T1 = sizeDepParams.a2 + sizeDepParams.a1 * Vgsteff - 0.0001;
    T2 = sqrt(T1 * T1 + 0.0004 * sizeDepParams.a2);
    Lambda = 0.5 * (T1 + T2);
    dLambda_dVg = 0.5 * sizeDepParams.a1 * (1.0 + T1 / T2);
  }

  Vgst2Vtm = Vgsteff + 2.0 * Vtm;
  if (Rds > 0)
  {
    tmp2 = dRds_dVg / Rds + dWeff_dVg / Weff;
    tmp3 = dRds_dVb / Rds + dWeff_dVb / Weff;
  }
  else
  {
    tmp2 = dWeff_dVg / Weff;
    tmp3 = dWeff_dVb / Weff;
  }

  if ((Rds == 0.0) && (Lambda == 1.0))
  {
    T0 = 1.0 / (Abulk * EsatL + Vgst2Vtm);
    tmp1 = 0.0;
    T1 = T0 * T0;
    T2 = Vgst2Vtm * T0;
    T3 = EsatL * Vgst2Vtm;
    Vdsat = T3 * T0;

    dT0_dVg = -(Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 1.0) * T1;
    dT0_dVd = -(Abulk * dEsatL_dVd) * T1;
    dT0_dVb = -(Abulk * dEsatL_dVb + dAbulk_dVb * EsatL) * T1;

    dVdsat_dVg = T3 * dT0_dVg + T2 * dEsatL_dVg + EsatL * T0;
    dVdsat_dVd = T3 * dT0_dVd + T2 * dEsatL_dVd;
    dVdsat_dVb = T3 * dT0_dVb + T2 * dEsatL_dVb;
  }
  else
  {
    tmp1 = dLambda_dVg / (Lambda * Lambda);
    T9 = Abulk * WVCoxRds;
    T8 = Abulk * T9;
    T7 = Vgst2Vtm * T9;
    T6 = Vgst2Vtm * WVCoxRds;
    T0 = 2.0 * Abulk * (T9 - 1.0 + 1.0 / Lambda);
    dT0_dVg = 2.0 * (T8 * tmp2 - Abulk * tmp1
        + (2.0 * T9 + 1.0 / Lambda - 1.0) * dAbulk_dVg);

    dT0_dVb = 2.0 * (T8 * (2.0 / Abulk * dAbulk_dVb + tmp3)
        + (1.0 / Lambda - 1.0) * dAbulk_dVb);
    dT0_dVd = 0.0;
    T1 = Vgst2Vtm * (2.0 / Lambda - 1.0) + Abulk * EsatL + 3.0 * T7;

    dT1_dVg = (2.0 / Lambda - 1.0) - 2.0 * Vgst2Vtm * tmp1
      + Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 3.0 * (T9
          + T7 * tmp2 + T6 * dAbulk_dVg);

    dT1_dVb = Abulk * dEsatL_dVb + EsatL * dAbulk_dVb
      + 3.0 * (T6 * dAbulk_dVb + T7 * tmp3);

    dT1_dVd = Abulk * dEsatL_dVd;

    T2 = Vgst2Vtm * (EsatL + 2.0 * T6);
    dT2_dVg = EsatL + Vgst2Vtm * dEsatL_dVg
      + T6 * (4.0 + 2.0 * Vgst2Vtm * tmp2);

    dT2_dVb = Vgst2Vtm * (dEsatL_dVb + 2.0 * T6 * tmp3);
    dT2_dVd = Vgst2Vtm * dEsatL_dVd;

    T3 = sqrt(T1 * T1 - 2.0 * T0 * T2);
    Vdsat = (T1 - T3) / T0;

    dT3_dVg = (T1 * dT1_dVg - 2.0 * (T0 * dT2_dVg + T2 * dT0_dVg)) / T3;
    dT3_dVd = (T1 * dT1_dVd - 2.0 * (T0 * dT2_dVd + T2 * dT0_dVd)) / T3;
    dT3_dVb = (T1 * dT1_dVb - 2.0 * (T0 * dT2_dVb + T2 * dT0_dVb)) / T3;

    dVdsat_dVg = (dT1_dVg - (T1 * dT1_dVg - dT0_dVg * T2
          - T0 * dT2_dVg) / T3 - Vdsat * dT0_dVg) / T0;

    dVdsat_dVb = (dT1_dVb - (T1 * dT1_dVb - dT0_dVb * T2
          - T0 * dT2_dVb) / T3 - Vdsat * dT0_dVb) / T0;

    dVdsat_dVd = (dT1_dVd - (T1 * dT1_dVd - T0 * dT2_dVd) / T3) / T0;
  }
  vdsat = Vdsat;

  // Effective Vds (Vdseff) Calculation
  T1 = Vdsat - Vds - sizeDepParams.delta;
  dT1_dVg = dVdsat_dVg;
  dT1_dVd = dVdsat_dVd - 1.0;
  dT1_dVb = dVdsat_dVb;

  T2 = sqrt(T1 * T1 + 4.0 * sizeDepParams.delta * Vdsat);
  T0 = T1 / T2;
  T3 = 2.0 * sizeDepParams.delta / T2;
  dT2_dVg = T0 * dT1_dVg + T3 * dVdsat_dVg;
  dT2_dVd = T0 * dT1_dVd + T3 * dVdsat_dVd;
  dT2_dVb = T0 * dT1_dVb + T3 * dVdsat_dVb;

  Vdseff = Vdsat - 0.5 * (T1 + T2);
  dVdseff_dVg = dVdsat_dVg - 0.5 * (dT1_dVg + dT2_dVg);
  dVdseff_dVd = dVdsat_dVd - 0.5 * (dT1_dVd + dT2_dVd);
  dVdseff_dVb = dVdsat_dVb - 0.5 * (dT1_dVb + dT2_dVb);

  // Added to eliminate non-zero Vdseff at Vds=0.0
  if (Vds == 0.0)
  {
    Vdseff = 0.0;
    dVdseff_dVg = 0.0;
    dVdseff_dVb = 0.0;
  }

  // Calculate VAsat
  tmp4 = 1.0 - 0.5 * Abulk * Vdsat / Vgst2Vtm;
  T9 = WVCoxRds * Vgsteff;
  T8 = T9 / Vgst2Vtm;
  T0 = EsatL + Vdsat + 2.0 * T9 * tmp4;

  T7 = 2.0 * WVCoxRds * tmp4;
  dT0_dVg = dEsatL_dVg + dVdsat_dVg + T7 * (1.0 + tmp2 * Vgsteff)
    - T8 * (Abulk * dVdsat_dVg - Abulk * Vdsat / Vgst2Vtm
        + Vdsat * dAbulk_dVg);

  dT0_dVb = dEsatL_dVb + dVdsat_dVb + T7 * tmp3 * Vgsteff
    - T8 * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
  dT0_dVd = dEsatL_dVd + dVdsat_dVd - T8 * Abulk * dVdsat_dVd;

  T9 = WVCoxRds * Abulk;
  T1 = 2.0 / Lambda - 1.0 + T9;
  dT1_dVg = -2.0 * tmp1 +  WVCoxRds * (Abulk * tmp2 + dAbulk_dVg);
  dT1_dVb = dAbulk_dVb * WVCoxRds + T9 * tmp3;

  Vasat = T0 / T1;
  dVasat_dVg = (dT0_dVg - Vasat * dT1_dVg) / T1;
  dVasat_dVb = (dT0_dVb - Vasat * dT1_dVb) / T1;
  dVasat_dVd = dT0_dVd / T1;

  if (Vdseff > Vds) Vdseff = Vds;

  diffVds = Vds - Vdseff;

  // Calculate VACLM
  if ((sizeDepParams.pclm > 0.0) && (diffVds > 1.0e-10))
  {
    T0 = 1.0 / (sizeDepParams.pclm * Abulk * sizeDepParams.litl);
    dT0_dVb = -T0 / Abulk * dAbulk_dVb;
    dT0_dVg = -T0 / Abulk * dAbulk_dVg;

    T2 = Vgsteff / EsatL;
    T1 = Leff * (Abulk + T2);
    dT1_dVg = Leff * ((1.0 - T2 * dEsatL_dVg) / EsatL + dAbulk_dVg);
    dT1_dVb = Leff * (dAbulk_dVb - T2 * dEsatL_dVb / EsatL);
    dT1_dVd = -T2 * dEsatL_dVd / Esat;

    T9 = T0 * T1;
    VACLM = T9 * diffVds;
    dVACLM_dVg = T0 * dT1_dVg * diffVds - T9 * dVdseff_dVg
      + T1 * diffVds * dT0_dVg;

    dVACLM_dVb = (dT0_dVb * T1 + T0 * dT1_dVb) * diffVds
      - T9 * dVdseff_dVb;

    dVACLM_dVd = T0 * dT1_dVd * diffVds + T9 * (1.0 - dVdseff_dVd);
  }
  else
  {
    VACLM = CONSTMAX_EXP;
    dVACLM_dVd = dVACLM_dVg = dVACLM_dVb = 0.0;
  }

  // Calculate VADIBL
  if (sizeDepParams.thetaRout > 0.0)
  {
    T8 = Abulk * Vdsat;
    T0 = Vgst2Vtm * T8;
    dT0_dVg = Vgst2Vtm * Abulk * dVdsat_dVg + T8
      + Vgst2Vtm * Vdsat * dAbulk_dVg;

    dT0_dVb = Vgst2Vtm * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
    dT0_dVd = Vgst2Vtm * Abulk * dVdsat_dVd;

    T1 = Vgst2Vtm + T8;
    dT1_dVg = 1.0 + Abulk * dVdsat_dVg + Vdsat * dAbulk_dVg;
    dT1_dVb = Abulk * dVdsat_dVb + dAbulk_dVb * Vdsat;
    dT1_dVd = Abulk * dVdsat_dVd;

    T9 = T1 * T1;
    T2 = sizeDepParams.thetaRout;

    VADIBL = (Vgst2Vtm - T0 / T1) / T2;
    dVADIBL_dVg = (1.0 - dT0_dVg / T1 + T0 * dT1_dVg / T9) / T2;
    dVADIBL_dVb = (-dT0_dVb / T1 + T0 * dT1_dVb / T9) / T2;
    dVADIBL_dVd = (-dT0_dVd / T1 + T0 * dT1_dVd / T9) / T2;

    T7 = sizeDepParams.pdiblb * Vbseff;
    if (T7 >= -0.9)
    {
      T3 = 1.0 / (1.0 + T7);
      VADIBL *= T3;
      dVADIBL_dVg *= T3;
      dVADIBL_dVb = (dVADIBL_dVb - VADIBL * sizeDepParams.pdiblb) * T3;
      dVADIBL_dVd *= T3;
    }
    else // Added to avoid the discontinuity problem caused by pdiblcb
    {
      T4 = 1.0 / (0.8 + T7);
      T3 = (17.0 + 20.0 * T7) * T4;
      dVADIBL_dVg *= T3;
      dVADIBL_dVb = dVADIBL_dVb * T3 - VADIBL * sizeDepParams.pdiblb * T4 * T4;

      dVADIBL_dVd *= T3;
      VADIBL *= T3;
    }
  }
  else
  {
    VADIBL = CONSTMAX_EXP;
    dVADIBL_dVd = dVADIBL_dVg = dVADIBL_dVb = 0.0;
  }

  // Calculate VA
  T8 = sizeDepParams.pvag / EsatL;
  T9 = T8 * Vgsteff;
  if (T9 > -0.9)
  {
    T0 = 1.0 + T9;
    dT0_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL);
    dT0_dVb = -T9 * dEsatL_dVb / EsatL;
    dT0_dVd = -T9 * dEsatL_dVd / EsatL;
  }
  else /* Added to avoid the discontinuity problems caused by pvag */
  {
    T1 = 1.0 / (17.0 + 20.0 * T9);
    T0 = (0.8 + T9) * T1;
    T1 *= T1;
    dT0_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL) * T1;

    T9 *= T1 / EsatL;
    dT0_dVb = -T9 * dEsatL_dVb;
    dT0_dVd = -T9 * dEsatL_dVd;
  }

  tmp1 = VACLM * VACLM;
  tmp2 = VADIBL * VADIBL;
  tmp3 = VACLM + VADIBL;

  T1 = VACLM * VADIBL / tmp3;
  tmp3 *= tmp3;
  dT1_dVg = (tmp1 * dVADIBL_dVg + tmp2 * dVACLM_dVg) / tmp3;
  dT1_dVd = (tmp1 * dVADIBL_dVd + tmp2 * dVACLM_dVd) / tmp3;
  dT1_dVb = (tmp1 * dVADIBL_dVb + tmp2 * dVACLM_dVb) / tmp3;

  Va = Vasat + T0 * T1;
  dVa_dVg = dVasat_dVg + T1 * dT0_dVg + T0 * dT1_dVg;
  dVa_dVd = dVasat_dVd + T1 * dT0_dVd + T0 * dT1_dVd;
  dVa_dVb = dVasat_dVb + T1 * dT0_dVb + T0 * dT1_dVb;

  // Calculate VASCBE
  if (sizeDepParams.pscbe2 > 0.0)
  {
    if (diffVds > sizeDepParams.pscbe1 * sizeDepParams.litl / CONSTEXP_THRESHOLD)
    {
      T0 =  sizeDepParams.pscbe1 * sizeDepParams.litl / diffVds;
      VASCBE = Leff * exp(T0) / sizeDepParams.pscbe2;
      T1 = T0 * VASCBE / diffVds;
      dVASCBE_dVg = T1 * dVdseff_dVg;
      dVASCBE_dVd = -T1 * (1.0 - dVdseff_dVd);
      dVASCBE_dVb = T1 * dVdseff_dVb;
    }
    else
    {
      VASCBE = CONSTMAX_EXP * Leff/sizeDepParams.pscbe2;
      dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
    }
  }
  else
  {
    VASCBE = CONSTMAX_EXP;
    dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
  }

  // Calculate Ids
  CoxWovL = model_cox * Weff / Leff;
  beta = ueff * CoxWovL;
  dbeta_dVg = CoxWovL * dueff_dVg + beta * dWeff_dVg / Weff;
  dbeta_dVd = CoxWovL * dueff_dVd;
  dbeta_dVb = CoxWovL * dueff_dVb + beta * dWeff_dVb / Weff;

  T0 = 1.0 - 0.5 * Abulk * Vdseff / Vgst2Vtm;
  dT0_dVg = -0.5 * (Abulk * dVdseff_dVg
      - Abulk * Vdseff / Vgst2Vtm + Vdseff * dAbulk_dVg) / Vgst2Vtm;
  dT0_dVd = -0.5 * Abulk * dVdseff_dVd / Vgst2Vtm;
  dT0_dVb = -0.5 * (Abulk * dVdseff_dVb + dAbulk_dVb * Vdseff) / Vgst2Vtm;

  fgche1 = Vgsteff * T0;
  dfgche1_dVg = Vgsteff * dT0_dVg + T0;
  dfgche1_dVd = Vgsteff * dT0_dVd;
  dfgche1_dVb = Vgsteff * dT0_dVb;

  T9 = Vdseff / EsatL;
  fgche2 = 1.0 + T9;
  dfgche2_dVg = (dVdseff_dVg - T9 * dEsatL_dVg) / EsatL;
  dfgche2_dVd = (dVdseff_dVd - T9 * dEsatL_dVd) / EsatL;
  dfgche2_dVb = (dVdseff_dVb - T9 * dEsatL_dVb) / EsatL;

  gche = beta * fgche1 / fgche2;
  dgche_dVg = (beta * dfgche1_dVg + fgche1 * dbeta_dVg
      - gche * dfgche2_dVg) / fgche2;

  dgche_dVd = (beta * dfgche1_dVd + fgche1 * dbeta_dVd
      - gche * dfgche2_dVd) / fgche2;

  dgche_dVb = (beta * dfgche1_dVb + fgche1 * dbeta_dVb
      - gche * dfgche2_dVb) / fgche2;

  T0 = 1.0 + gche * Rds;
  T9 = Vdseff / T0;
  Idl = gche * T9;

  dIdl_dVg = (gche * dVdseff_dVg + T9 * dgche_dVg) / T0
    - Idl * gche / T0 * dRds_dVg ;

  dIdl_dVd = (gche * dVdseff_dVd + T9 * dgche_dVd) / T0;
  dIdl_dVb = (gche * dVdseff_dVb + T9 * dgche_dVb
      - Idl * dRds_dVb * gche) / T0;

  T9 =  diffVds / Va;
  T0 =  1.0 + T9;
  Idsa = Idl * T0;
  dIdsa_dVg = T0 * dIdl_dVg - Idl * (dVdseff_dVg + T9 * dVa_dVg) / Va;
  dIdsa_dVd = T0 * dIdl_dVd + Idl * (1.0 - dVdseff_dVd
      - T9 * dVa_dVd) / Va;

  dIdsa_dVb = T0 * dIdl_dVb - Idl * (dVdseff_dVb + T9 * dVa_dVb) / Va;

  T9 = diffVds / VASCBE;
  T0 = 1.0 + T9;
  Ids = Idsa * T0;

  Gm  = T0 * dIdsa_dVg - Idsa * (dVdseff_dVg + T9 * dVASCBE_dVg) / VASCBE;
  Gds = T0 * dIdsa_dVd + Idsa * (1.0 - dVdseff_dVd
      - T9 * dVASCBE_dVd) / VASCBE;
  Gmb = T0 * dIdsa_dVb - Idsa * (dVdseff_dVb
      + T9 * dVASCBE_dVb) / VASCBE;

  Gds += Gm * dVgsteff_dVd;
  Gmb += Gm * dVgsteff_dVb;
  Gm *= dVgsteff_dVg;
  Gmb *= dVbseff_dVb;

  // Substrate current begins
  tmp = sizeDepParams.alpha0 + sizeDepParams.alpha1 * Leff;
  if ((tmp <= 0.0) || (sizeDepParams.beta0 <= 0.0))
  {
    Isub = Gbd = Gbb = Gbg = 0.0;
  }
  else
  {
    T2 = tmp / Leff;
    if (diffVds > sizeDepParams.beta0 / CONSTEXP_THRESHOLD)
    {
      T0 = -sizeDepParams.beta0 / diffVds;
      T1 = T2 * diffVds * exp(T0);
      T3 = T1 / diffVds * (T0 - 1.0);
      dT1_dVg = T3 * dVdseff_dVg;
      dT1_dVd = T3 * (dVdseff_dVd - 1.0);
      dT1_dVb = T3 * dVdseff_dVb;
    }
    else
    {
      T3 = T2 * CONSTMIN_EXP;
      T1 = T3 * diffVds;
      dT1_dVg = -T3 * dVdseff_dVg;
      dT1_dVd = T3 * (1.0 - dVdseff_dVd);
      dT1_dVb = -T3 * dVdseff_dVb;
    }
    Isub = T1 * Idsa;
    Gbg = T1 * dIdsa_dVg + Idsa * dT1_dVg;
    Gbd = T1 * dIdsa_dVd + Idsa * dT1_dVd;
    Gbb = T1 * dIdsa_dVb + Idsa * dT1_dVb;

    Gbd += Gbg * dVgsteff_dVd;
    Gbb += Gbg * dVgsteff_dVb;
    Gbg *= dVgsteff_dVg;
    Gbb *= dVbseff_dVb; // bug fixing
  }

  // copy over local drain (channel) current vars to instance vars:
  cdrain = Ids;
  gds = Gds;
  gm = Gm;
  gmbs = Gmb;

  // copy over local substrate current vars to instance vars:
  gbbs = Gbb;
  gbgs = Gbg;
  gbds = Gbd;

  csub = Isub;

  //  thermal noise Qinv calculated from all capMod
  //  * 0, 1, 2 & 3 stored in iterI->qinv 1/1998
  if ((model_xpart < 0) || (!ChargeComputationNeeded))
  {
    qgate  = qdrn = qsrc = qbulk = 0.0;
    cggb = cgsb = cgdb = 0.0;
    cdgb = cdsb = cddb = 0.0;
    cbgb = cbsb = cbdb = 0.0;
    cqdb = cqsb = cqgb = cqbb = 0.0;

    gtau = 0.0;
    goto finished;
  }
  else if (model_capMod == 0)
  {
    if (Vbseff < 0.0)
    {
      Vbseff = Vbs;
      dVbseff_dVb = 1.0;
    }
    else
    {
      Vbseff = sizeDepParams.phi - Phis;
      dVbseff_dVb = -dPhis_dVb;
    }

    Vfb = sizeDepParams.vfbcv;
    Vth = Vfb + sizeDepParams.phi + sizeDepParams.k1ox * sqrtPhis;
    Vgst = Vgs_eff - Vth;
    dVth_dVb = sizeDepParams.k1ox * dsqrtPhis_dVb;
    dVgst_dVb = -dVth_dVb;
    dVgst_dVg = dVgs_eff_dVg;

    CoxWL = model_cox * sizeDepParams.weffCV * sizeDepParams.leffCV;
    Arg1 = Vgs_eff - Vbseff - Vfb;

    if (Arg1 <= 0.0)
    {
      qgate = CoxWL * Arg1;
      qbulk = -qgate;
      qdrn = 0.0;

      cggb = CoxWL * dVgs_eff_dVg;
      cgdb = 0.0;
      cgsb = CoxWL * (dVbseff_dVb - dVgs_eff_dVg);

      cdgb = 0.0;
      cddb = 0.0;
      cdsb = 0.0;

      cbgb = -CoxWL * dVgs_eff_dVg;
      cbdb = 0.0;
      cbsb = -cgsb;
      qinv = 0.0;
    }
    else if (Vgst <= 0.0)
    {
      T1 = 0.5 * sizeDepParams.k1ox;
      T2 = sqrt(T1 * T1 + Arg1);
      qgate = CoxWL * sizeDepParams.k1ox * (T2 - T1);
      qbulk = -qgate;
      qdrn = 0.0;

      T0 = CoxWL * T1 / T2;
      cggb = T0 * dVgs_eff_dVg;
      cgdb = 0.0;
      cgsb = T0 * (dVbseff_dVb - dVgs_eff_dVg);

      cdgb = 0.0;
      cddb = 0.0;
      cdsb = 0.0;

      cbgb = -cggb;
      cbdb = 0.0;
      cbsb = -cgsb;
      qinv = 0.0;
    }
    else
    {
      One_Third_CoxWL = CoxWL / 3.0;
      Two_Third_CoxWL = 2.0 * One_Third_CoxWL;

      AbulkCV = Abulk0 * sizeDepParams.abulkCVfactor;
      dAbulkCV_dVb = sizeDepParams.abulkCVfactor * dAbulk0_dVb;

      Vdsat = Vgst / AbulkCV;
      dVdsat_dVg = dVgs_eff_dVg / AbulkCV;
      dVdsat_dVb = - (Vdsat * dAbulkCV_dVb + dVth_dVb)/ AbulkCV;

      if (model_xpart > 0.5)
      { // 0/100 Charge partition model
        if (Vdsat <= Vds)
        {   // saturation region
          T1 = Vdsat / 3.0;
          qgate = CoxWL * (Vgs_eff - Vfb - sizeDepParams.phi - T1);

          T2 = -Two_Third_CoxWL * Vgst;
          qbulk = -(qgate + T2);
          qdrn = 0.0;

          cggb = One_Third_CoxWL * (3.0 - dVdsat_dVg)* dVgs_eff_dVg;
          T2 = -One_Third_CoxWL * dVdsat_dVb;

          cgsb = -(cggb + T2);
          cgdb = 0.0;

          cdgb = 0.0;
          cddb = 0.0;
          cdsb = 0.0;

          cbgb = -(cggb - Two_Third_CoxWL * dVgs_eff_dVg);
          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);

          cbsb = -(cbgb + T3);
          cbdb = 0.0;
          qinv = -(qgate + qbulk);
        }
        else
        {   // linear region
          Alphaz = Vgst / Vdsat;
          T1 = 2.0 * Vdsat - Vds;
          T2 = Vds / (3.0 * T1);
          T3 = T2 * Vds;
          T9 = 0.25 * CoxWL;
          T4 = T9 * Alphaz;
          T7 = 2.0 * Vds - T1 - 3.0 * T3;
          T8 = T3 - T1 - 2.0 * Vds;
          qgate = CoxWL * (Vgs_eff - Vfb - sizeDepParams.phi - 0.5 * (Vds - T3));

          T10 = T4 * T8;
          qdrn = T4 * T7;
          qbulk = -(qgate + qdrn + T10);

          T5 = T3 / T1;
          cggb = CoxWL * (1.0 - T5 * dVdsat_dVg) * dVgs_eff_dVg;

          T11 = -CoxWL * T5 * dVdsat_dVb;
          cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
          cgsb = -(cggb + T11 + cgdb);

          T6 = 1.0 / Vdsat;
          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

          T7 = T9 * T7;
          T8 = T9 * T8;
          T9 = 2.0 * T4 * (1.0 - 3.0 * T5);
          cdgb = (T7 * dAlphaz_dVg - T9* dVdsat_dVg) * dVgs_eff_dVg;

          T12 = T7 * dAlphaz_dVb - T9 * dVdsat_dVb;
          cddb = T4 * (3.0 - 6.0 * T2 - 3.0 * T5);
          cdsb = -(cdgb + T12 + cddb);

          T9  = 2.0 * T4 * (1.0 + T5);
          T10 = (T8 * dAlphaz_dVg - T9 * dVdsat_dVg) * dVgs_eff_dVg;
          T11 = T8 * dAlphaz_dVb - T9 * dVdsat_dVb;
          T12 = T4 * (2.0 * T2 + T5 - 1.0);
          T0  = -(T10 + T11 + T12);


          cbgb = -(cggb + cdgb + T10);
          cbdb = -(cgdb + cddb + T12);
          cbsb = -(cgsb + cdsb + T0);
          qinv = -(qgate + qbulk);
        }
      }
      else if (model_xpart < 0.5)
      {   // 40/60 Charge partition model
        if (Vds >= Vdsat)
        { // saturation region
          T1 = Vdsat / 3.0;
          qgate = CoxWL * (Vgs_eff - Vfb - sizeDepParams.phi - T1);

          T2 = -Two_Third_CoxWL * Vgst;
          qbulk = -(qgate + T2);
          qdrn = 0.4 * T2;

          cggb = One_Third_CoxWL* (3.0 - dVdsat_dVg) * dVgs_eff_dVg;

          T2 = -One_Third_CoxWL * dVdsat_dVb;
          cgsb = -(cggb + T2);
          cgdb = 0.0;

          T3 = 0.4 * Two_Third_CoxWL;
          cdgb = -T3 * dVgs_eff_dVg;
          cddb = 0.0;

          T4 = T3 * dVth_dVb;
          cdsb = -(T4 + cdgb);

          cbgb = -(cggb - Two_Third_CoxWL * dVgs_eff_dVg);
          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
          cbsb = -(cbgb + T3);
          cbdb = 0.0;
          qinv = -(qgate + qbulk);
        }
        else
        {   // linear region
          Alphaz = Vgst / Vdsat;
          T1 = 2.0 * Vdsat - Vds;
          T2 = Vds / (3.0 * T1);
          T3 = T2 * Vds;
          T9 = 0.25 * CoxWL;
          T4 = T9 * Alphaz;
          qgate = CoxWL * (Vgs_eff - Vfb - sizeDepParams.phi - 0.5 * (Vds - T3));

          T5 = T3 / T1;
          cggb = CoxWL * (1.0 - T5 * dVdsat_dVg) * dVgs_eff_dVg;
          tmp = -CoxWL * T5 * dVdsat_dVb;
          cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
          cgsb = -(cggb + cgdb + tmp);

          T6 = 1.0 / Vdsat;
          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

          T6 = 8.0 * Vdsat * Vdsat - 6.0 * Vdsat * Vds + 1.2 * Vds * Vds;
          T8 = T2 / T1;
          T7 = Vds - T1 - T8 * T6;
          qdrn = T4 * T7;
          T7 *= T9;
          tmp = T8 / T1;
          tmp1 = T4*(2.0 - 4.0 * tmp * T6 + T8 *(16.0 * Vdsat - 6.0 *Vds));

          cdgb = (T7 *dAlphaz_dVg - tmp1 *dVdsat_dVg) *dVgs_eff_dVg;

          T10 = T7 * dAlphaz_dVb - tmp1 * dVdsat_dVb;
          cddb = T4 * (2.0 - (1.0 / (3.0 * T1
                  * T1) + 2.0 * tmp) * T6 + T8
              * (6.0 * Vdsat - 2.4 * Vds));

          cdsb = -(cdgb + T10 + cddb);

          T7 = 2.0 * (T1 + T3);
          qbulk = -(qgate - T4 * T7);
          T7 *= T9;
          T0 = 4.0 * T4 * (1.0 - T5);
          T12 = (-T7 * dAlphaz_dVg - cdgb
              - T0 * dVdsat_dVg) * dVgs_eff_dVg;
          T11 = -T7 * dAlphaz_dVb - T10 - T0 * dVdsat_dVb;
          T10 = -4.0 * T4 * (T2 - 0.5 + 0.5 * T5)
            - cddb;

          tmp = -(T10 + T11 + T12);

          cbgb = -(cggb + cdgb + T12);
          cbdb = -(cgdb + cddb + T11);
          cbsb = -(cgsb + cdsb + tmp);
          qinv = -(qgate + qbulk);
        }
      }
      else
      {   // 50/50 partitioning
        if (Vds >= Vdsat)
        {   // saturation region
          T1 = Vdsat / 3.0;
          qgate = CoxWL * (Vgs_eff - Vfb - sizeDepParams.phi - T1);
          T2 = -Two_Third_CoxWL * Vgst;
          qbulk = -(qgate + T2);
          qdrn = 0.5 * T2;

          cggb = One_Third_CoxWL * (3.0 -dVdsat_dVg) * dVgs_eff_dVg;
          T2 = -One_Third_CoxWL * dVdsat_dVb;
          cgsb = -(cggb + T2);
          cgdb = 0.0;

          cdgb = -One_Third_CoxWL * dVgs_eff_dVg;
          cddb = 0.0;
          T4 = One_Third_CoxWL * dVth_dVb;
          cdsb = -(T4 + cdgb);

          cbgb = -(cggb - Two_Third_CoxWL * dVgs_eff_dVg);
          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
          cbsb = -(cbgb + T3);
          cbdb = 0.0;
          qinv = -(qgate + qbulk);
        }
        else
        {   // linear region
          Alphaz = Vgst / Vdsat;
          T1 = 2.0 * Vdsat - Vds;
          T2 = Vds / (3.0 * T1);
          T3 = T2 * Vds;
          T9 = 0.25 * CoxWL;
          T4 = T9 * Alphaz;
          qgate = CoxWL * (Vgs_eff - Vfb - sizeDepParams.phi - 0.5 * (Vds - T3))
            ;

          T5 = T3 / T1;
          cggb = CoxWL * (1.0 - T5 * dVdsat_dVg) * dVgs_eff_dVg;

          tmp = -CoxWL * T5 * dVdsat_dVb;
          cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
          cgsb = -(cggb + cgdb + tmp);

          T6 = 1.0 / Vdsat;
          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

          T7 = T1 + T3;
          qdrn = -T4 * T7;
          qbulk = - (qgate + qdrn + qdrn);
          T7 *= T9;
          T0 = T4 * (2.0 * T5 - 2.0);

          cdgb = (T0 * dVdsat_dVg - T7 *dAlphaz_dVg) *dVgs_eff_dVg;
          T12 = T0 * dVdsat_dVb - T7 * dAlphaz_dVb;
          cddb = T4 * (1.0 - 2.0 * T2 - T5);
          cdsb = -(cdgb + T12 + cddb);

          cbgb = -(cggb + 2.0 * cdgb);
          cbdb = -(cgdb + 2.0 * cddb);
          cbsb = -(cgsb + 2.0 * cdsb);
          qinv = -(qgate + qbulk);
        }
      }
    }
  }
  else
  {
    if (Vbseff < 0.0)
    {
      VbseffCV = Vbseff;
      dVbseffCV_dVb = 1.0;
    }
    else
    {
      VbseffCV = sizeDepParams.phi - Phis;
      dVbseffCV_dVb = -dPhis_dVb;
    }

    CoxWL = model_cox * sizeDepParams.weffCV * sizeDepParams.leffCV;

    // Seperate VgsteffCV with noff and voffcv
    noff = n * sizeDepParams.noff;
    dnoff_dVd = sizeDepParams.noff * dn_dVd;
    dnoff_dVb = sizeDepParams.noff * dn_dVb;
    T0 = Vtm * noff;
    voffcv = sizeDepParams.voffcv;
    VgstNVt = (Vgst - voffcv) / T0;

    if (VgstNVt > CONSTEXP_THRESHOLD)
    {
      Vgsteff = Vgst - voffcv;
      dVgsteff_dVg = dVgs_eff_dVg;
      dVgsteff_dVd = -dVth_dVd;
      dVgsteff_dVb = -dVth_dVb;
    }
    else if (VgstNVt < -CONSTEXP_THRESHOLD)
    {
      Vgsteff = T0 * log(1.0 + CONSTMIN_EXP);
      dVgsteff_dVg = 0.0;
      dVgsteff_dVd = Vgsteff / noff;
      dVgsteff_dVb = dVgsteff_dVd * dnoff_dVb;
      dVgsteff_dVd *= dnoff_dVd;
    }
    else
    {
      ExpVgst = exp(VgstNVt);
      Vgsteff = T0 * log(1.0 + ExpVgst);
      dVgsteff_dVg = ExpVgst / (1.0 + ExpVgst);
      dVgsteff_dVd = -dVgsteff_dVg * (dVth_dVd + (Vgst - voffcv)
          / noff * dnoff_dVd) + Vgsteff / noff * dnoff_dVd;
      dVgsteff_dVb = -dVgsteff_dVg * (dVth_dVb + (Vgst - voffcv)
          / noff * dnoff_dVb) + Vgsteff / noff * dnoff_dVb;
      dVgsteff_dVg *= dVgs_eff_dVg;
    } // End of VgsteffCV

    if (model_capMod == 1)
    {
      Vfb = sizeDepParams.vfbzb;
      Arg1 = Vgs_eff - VbseffCV - Vfb - Vgsteff;

      if (Arg1 <= 0.0)
      {
        qgate = CoxWL * Arg1;
        Cgg = CoxWL * (dVgs_eff_dVg - dVgsteff_dVg);
        Cgd = -CoxWL * dVgsteff_dVd;
        Cgb = -CoxWL * (dVbseffCV_dVb + dVgsteff_dVb);
      }
      else
      {
        T0 = 0.5 * sizeDepParams.k1ox;
        T1 = sqrt(T0 * T0 + Arg1);
        T2 = CoxWL * T0 / T1;

        qgate = CoxWL * sizeDepParams.k1ox * (T1 - T0);

        Cgg = T2 * (dVgs_eff_dVg - dVgsteff_dVg);
        Cgd = -T2 * dVgsteff_dVd;
        Cgb = -T2 * (dVbseffCV_dVb + dVgsteff_dVb);
      }
      qbulk = -qgate;
      Cbg = -Cgg;
      Cbd = -Cgd;
      Cbb = -Cgb;

      One_Third_CoxWL = CoxWL / 3.0;
      Two_Third_CoxWL = 2.0 * One_Third_CoxWL;
      AbulkCV = Abulk0 * sizeDepParams.abulkCVfactor;
      dAbulkCV_dVb = sizeDepParams.abulkCVfactor * dAbulk0_dVb;
      VdsatCV = Vgsteff / AbulkCV;

      if (VdsatCV < Vds)
      {
        dVdsatCV_dVg = 1.0 / AbulkCV;
        dVdsatCV_dVb = -VdsatCV * dAbulkCV_dVb / AbulkCV;
        T0 = Vgsteff - VdsatCV / 3.0;
        dT0_dVg = 1.0 - dVdsatCV_dVg / 3.0;
        dT0_dVb = -dVdsatCV_dVb / 3.0;
        qgate += CoxWL * T0;
        Cgg1 = CoxWL * dT0_dVg;
        Cgb1 = CoxWL * dT0_dVb + Cgg1 * dVgsteff_dVb;
        Cgd1 = Cgg1 * dVgsteff_dVd;
        Cgg1 *= dVgsteff_dVg;
        Cgg += Cgg1;
        Cgb += Cgb1;
        Cgd += Cgd1;

        T0 = VdsatCV - Vgsteff;
        dT0_dVg = dVdsatCV_dVg - 1.0;
        dT0_dVb = dVdsatCV_dVb;
        qbulk += One_Third_CoxWL * T0;
        Cbg1 = One_Third_CoxWL * dT0_dVg;
        Cbb1 = One_Third_CoxWL * dT0_dVb + Cbg1 * dVgsteff_dVb;
        Cbd1 = Cbg1 * dVgsteff_dVd;
        Cbg1 *= dVgsteff_dVg;
        Cbg += Cbg1;
        Cbb += Cbb1;
        Cbd += Cbd1;

        if (model_xpart > 0.5)      T0 = -Two_Third_CoxWL;
        else if (model_xpart < 0.5) T0 = -0.4 * CoxWL;
        else                         T0 = -One_Third_CoxWL;

        qsrc = T0 * Vgsteff;
        Csg = T0 * dVgsteff_dVg;
        Csb = T0 * dVgsteff_dVb;
        Csd = T0 * dVgsteff_dVd;
        Cgb *= dVbseff_dVb;
        Cbb *= dVbseff_dVb;
        Csb *= dVbseff_dVb;
      }
      else
      {
        T0 = AbulkCV * Vds;
        T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.e-20);
        T2 = Vds / T1;
        T3 = T0 * T2;
        dT3_dVg = -12.0 * T2 * T2 * AbulkCV;
        dT3_dVd = 6.0 * T0 * (4.0 * Vgsteff - T0) / T1 / T1 - 0.5;
        dT3_dVb = 12.0 * T2 * T2 * dAbulkCV_dVb * Vgsteff;

        qgate += CoxWL * (Vgsteff - 0.5 * Vds + T3);
        Cgg1 = CoxWL * (1.0 + dT3_dVg);
        Cgb1 = CoxWL * dT3_dVb + Cgg1 * dVgsteff_dVb;
        Cgd1 = CoxWL * dT3_dVd + Cgg1 * dVgsteff_dVd;
        Cgg1 *= dVgsteff_dVg;
        Cgg += Cgg1;
        Cgb += Cgb1;
        Cgd += Cgd1;

        qbulk += CoxWL * (1.0 - AbulkCV) * (0.5 * Vds - T3);
        Cbg1 = -CoxWL * ((1.0 - AbulkCV) * dT3_dVg);
        Cbb1 = -CoxWL * ((1.0 - AbulkCV) * dT3_dVb
            + (0.5 * Vds - T3) * dAbulkCV_dVb)
          + Cbg1 * dVgsteff_dVb;
        Cbd1 = -CoxWL * (1.0 - AbulkCV) * dT3_dVd
          + Cbg1 * dVgsteff_dVd;
        Cbg1 *= dVgsteff_dVg;
        Cbg += Cbg1;
        Cbb += Cbb1;
        Cbd += Cbd1;

        if (model_xpart > 0.5)
        {   // 0/100 Charge petition model
          T1 = T1 + T1;
          qsrc = -CoxWL * (0.5 * Vgsteff + 0.25 * T0 - T0 * T0 / T1);
          Csg = -CoxWL * (0.5 + 24.0 * T0 * Vds / T1 / T1 * AbulkCV);
          Csb = -CoxWL * (0.25 * Vds * dAbulkCV_dVb
              - 12.0 * T0 * Vds / T1 / T1 * (4.0 * Vgsteff - T0)
              * dAbulkCV_dVb) + Csg * dVgsteff_dVb;
          Csd = -CoxWL * (0.25 * AbulkCV - 12.0 * AbulkCV * T0
              / T1 / T1 * (4.0 * Vgsteff - T0))
            + Csg * dVgsteff_dVd;
          Csg *= dVgsteff_dVg;
        }
        else if (model_xpart < 0.5)
        {   // 40/60 Charge petition model
          T1 = T1 / 12.0;
          T2 = 0.5 * CoxWL / (T1 * T1);
          T3 = Vgsteff * (2.0 * T0 * T0 / 3.0 + Vgsteff
              * (Vgsteff - 4.0 * T0 / 3.0))
            - 2.0 * T0 * T0 * T0 / 15.0;
          qsrc = -T2 * T3;
          T4 = 4.0 / 3.0 * Vgsteff * (Vgsteff - T0)
            + 0.4 * T0 * T0;
          Csg = -2.0 * qsrc / T1 - T2 * (Vgsteff * (3.0
                * Vgsteff - 8.0 * T0 / 3.0)
              + 2.0 * T0 * T0 / 3.0);
          Csb = (qsrc / T1 * Vds + T2 * T4 * Vds) * dAbulkCV_dVb
            + Csg * dVgsteff_dVb;
          Csd = (qsrc / T1 + T2 * T4) * AbulkCV
            + Csg * dVgsteff_dVd;
          Csg *= dVgsteff_dVg;
        }
        else
        {   // 50/50 Charge petition model
          qsrc = -0.5 * (qgate + qbulk);
          Csg = -0.5 * (Cgg1 + Cbg1);
          Csb = -0.5 * (Cgb1 + Cbb1);
          Csd = -0.5 * (Cgd1 + Cbd1);
        }
        Cgb *= dVbseff_dVb;
        Cbb *= dVbseff_dVb;
        Csb *= dVbseff_dVb;
      }
      qdrn = -(qgate + qbulk + qsrc);
      cggb = Cgg;
      cgsb = -(Cgg + Cgd + Cgb);
      cgdb = Cgd;
      cdgb = -(Cgg + Cbg + Csg);
      cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb + Csg + Csd + Csb);
      cddb = -(Cgd + Cbd + Csd);
      cbgb = Cbg;
      cbsb = -(Cbg + Cbd + Cbb);
      cbdb = Cbd;
      qinv = -(qgate + qbulk);
    }
    else if (model_capMod == 2)
    {
      Vfb = sizeDepParams.vfbzb;
      V3 = Vfb - Vgs_eff + VbseffCV - CONSTDELTA_3;
      if (Vfb <= 0.0)
      {
        T0 = sqrt(V3 * V3 - 4.0 * CONSTDELTA_3 * Vfb);
        T2 = -CONSTDELTA_3 / T0;
      }
      else
      {
        T0 = sqrt(V3 * V3 + 4.0 * CONSTDELTA_3 * Vfb);
        T2 = CONSTDELTA_3 / T0;
      }

      T1 = 0.5 * (1.0 + V3 / T0);
      Vfbeff = Vfb - 0.5 * (V3 + T0);
      dVfbeff_dVg = T1 * dVgs_eff_dVg;
      dVfbeff_dVb = -T1 * dVbseffCV_dVb;
      Qac0 = CoxWL * (Vfbeff - Vfb);
      dQac0_dVg = CoxWL * dVfbeff_dVg;
      dQac0_dVb = CoxWL * dVfbeff_dVb;

      T0 = 0.5 * sizeDepParams.k1ox;
      T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
      if (sizeDepParams.k1ox == 0.0)
      {
        T1 = 0.0;
        T2 = 0.0;
      }
      else if (T3 < 0.0)
      {
        T1 = T0 + T3 / sizeDepParams.k1ox;
        T2 = CoxWL;
      }
      else
      {
        T1 = sqrt(T0 * T0 + T3);
        T2 = CoxWL * T0 / T1;
      }

      Qsub0 = CoxWL * sizeDepParams.k1ox * (T1 - T0);

      dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg);
      dQsub0_dVd = -T2 * dVgsteff_dVd;
      dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb + dVgsteff_dVb);

      AbulkCV = Abulk0 * sizeDepParams.abulkCVfactor;
      dAbulkCV_dVb = sizeDepParams.abulkCVfactor * dAbulk0_dVb;
      VdsatCV = Vgsteff / AbulkCV;

      V4 = VdsatCV - Vds - CONSTDELTA_4;
      T0 = sqrt(V4 * V4 + 4.0 * CONSTDELTA_4 * VdsatCV);
      VdseffCV = VdsatCV - 0.5 * (V4 + T0);
      T1 = 0.5 * (1.0 + V4 / T0);
      T2 = CONSTDELTA_4 / T0;
      T3 = (1.0 - T1 - T2) / AbulkCV;
      dVdseffCV_dVg = T3;
      dVdseffCV_dVd = T1;
      dVdseffCV_dVb = -T3 * VdsatCV * dAbulkCV_dVb;
      // Added to eliminate non-zero VdseffCV at Vds=0.0
      if (Vds == 0.0)
      {
        VdseffCV = 0.0;
        dVdseffCV_dVg = 0.0;
        dVdseffCV_dVb = 0.0;
      }

      T0 = AbulkCV * VdseffCV;
      T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1e-20);
      T2 = VdseffCV / T1;
      T3 = T0 * T2;

      T4 = (1.0 - 12.0 * T2 * T2 * AbulkCV);
      T5 = (6.0 * T0 * (4.0 * Vgsteff - T0) / (T1 * T1) - 0.5);
      T6 = 12.0 * T2 * T2 * Vgsteff;

      qinoi = -CoxWL * (Vgsteff - 0.5 * T0 + AbulkCV * T3);
      qgate = CoxWL * (Vgsteff - 0.5 * VdseffCV + T3);
      Cgg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
      Cgd1 = CoxWL * T5 * dVdseffCV_dVd + Cgg1 * dVgsteff_dVd;
      Cgb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
        + Cgg1 * dVgsteff_dVb;
      Cgg1 *= dVgsteff_dVg;

      T7 = 1.0 - AbulkCV;
      qbulk = CoxWL * T7 * (0.5 * VdseffCV - T3);
      T4 = -T7 * (T4 - 1.0);
      T5 = -T7 * T5;
      T6 = -(T7 * T6 + (0.5 * VdseffCV - T3));
      Cbg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
      Cbd1 = CoxWL * T5 * dVdseffCV_dVd + Cbg1 * dVgsteff_dVd;
      Cbb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
        + Cbg1 * dVgsteff_dVb;
      Cbg1 *= dVgsteff_dVg;

      if (model_xpart > 0.5)
      {   // 0/100 Charge petition model
        T1 = T1 + T1;
        qsrc = -CoxWL * (0.5 * Vgsteff + 0.25 * T0 - T0 * T0 / T1);
        T7 = (4.0 * Vgsteff - T0) / (T1 * T1);
        T4 = -(0.5 + 24.0 * T0 * T0 / (T1 * T1));
        T5 = -(0.25 * AbulkCV - 12.0 * AbulkCV * T0 * T7);
        T6 = -(0.25 * VdseffCV - 12.0 * T0 * VdseffCV * T7);
        Csg = CoxWL * (T4 + T5 * dVdseffCV_dVg);
        Csd = CoxWL * T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
        Csb = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
          + Csg * dVgsteff_dVb;
        Csg *= dVgsteff_dVg;
      }
      else if (model_xpart < 0.5)
      {   // 40/60 Charge petition model
        T1 = T1 / 12.0;
        T2 = 0.5 * CoxWL / (T1 * T1);
        T3 = Vgsteff *(2.0 * T0 *T0/3.0 +Vgsteff *(Vgsteff - 4.0 *T0/ 3.0))
          - 2.0 * T0 * T0 * T0 / 15.0;
        qsrc = -T2 * T3;
        T7 = 4.0 / 3.0 * Vgsteff * (Vgsteff - T0)
          + 0.4 * T0 * T0;
        T4 = -2.0 * qsrc / T1 - T2 * (Vgsteff * (3.0
              * Vgsteff - 8.0 * T0 / 3.0)
            + 2.0 * T0 * T0 / 3.0);
        T5 = (qsrc / T1 + T2 * T7) * AbulkCV;
        T6 = (qsrc / T1 * VdseffCV + T2 * T7 * VdseffCV);
        Csg = (T4 + T5 * dVdseffCV_dVg);
        Csd = T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
        Csb = (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
          + Csg * dVgsteff_dVb;
        Csg *= dVgsteff_dVg;
      }
      else
      {   // 50/50 Charge petition model
        qsrc = -0.5 * (qgate + qbulk);
        Csg = -0.5 * (Cgg1 + Cbg1);
        Csb = -0.5 * (Cgb1 + Cbb1);
        Csd = -0.5 * (Cgd1 + Cbd1);
      }

      qgate += Qac0 + Qsub0;
      qbulk -= (Qac0 + Qsub0);
      qdrn = -(qgate + qbulk + qsrc);

      Cgg = dQac0_dVg + dQsub0_dVg + Cgg1;
      Cgd = dQsub0_dVd + Cgd1;
      Cgb = dQac0_dVb + dQsub0_dVb + Cgb1;

      Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
      Cbd = Cbd1 - dQsub0_dVd;
      Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;

      Cgb *= dVbseff_dVb;
      Cbb *= dVbseff_dVb;
      Csb *= dVbseff_dVb;

      cggb = Cgg;
      cgsb = -(Cgg + Cgd + Cgb);
      cgdb = Cgd;
      cdgb = -(Cgg + Cbg + Csg);
      cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb + Csg + Csd + Csb);
      cddb = -(Cgd + Cbd + Csd);
      cbgb = Cbg;
      cbsb = -(Cbg + Cbd + Cbb);
      cbdb = Cbd;
      qinv = qinoi;
    }

    // New Charge-Thickness capMod (CTM) begins
    else if (model_capMod == 3)
    {
      V3 = sizeDepParams.vfbzb - Vgs_eff + VbseffCV - CONSTDELTA_3;
      if (sizeDepParams.vfbzb <= 0.0)
      {
        T0 = sqrt(V3 * V3 - 4.0 * CONSTDELTA_3 * sizeDepParams.vfbzb);
        T2 = -CONSTDELTA_3 / T0;
      }
      else
      {
        T0 = sqrt(V3 * V3 + 4.0 * CONSTDELTA_3 * sizeDepParams.vfbzb);
        T2 = CONSTDELTA_3 / T0;
      }

      T1 = 0.5 * (1.0 + V3 / T0);
      Vfbeff = sizeDepParams.vfbzb - 0.5 * (V3 + T0);
      dVfbeff_dVg = T1 * dVgs_eff_dVg;
      dVfbeff_dVb = -T1 * dVbseffCV_dVb;

      Cox = model_cox;
      Tox = 1.0e8 * model_tox;
      T0 = (Vgs_eff - VbseffCV - sizeDepParams.vfbzb) / Tox;
      dT0_dVg = dVgs_eff_dVg / Tox;
      dT0_dVb = -dVbseffCV_dVb / Tox;

      tmp = T0 * sizeDepParams.acde;
      if ((-CONSTEXP_THRESHOLD < tmp) && (tmp < CONSTEXP_THRESHOLD))
      {
        Tcen = sizeDepParams.ldeb * exp(tmp);
        dTcen_dVg = sizeDepParams.acde * Tcen;
        dTcen_dVb = dTcen_dVg * dT0_dVb;
        dTcen_dVg *= dT0_dVg;
      }
      else if (tmp <= -CONSTEXP_THRESHOLD)
      {
        Tcen = sizeDepParams.ldeb * CONSTMIN_EXP;
        dTcen_dVg = dTcen_dVb = 0.0;
      }
      else
      {
        Tcen = sizeDepParams.ldeb * CONSTMAX_EXP;
        dTcen_dVg = dTcen_dVb = 0.0;
      }

      LINK = 1.0e-3 * model_tox;
      V3 = sizeDepParams.ldeb - Tcen - LINK;
      V4 = sqrt(V3 * V3 + 4.0 * LINK * sizeDepParams.ldeb);
      Tcen = sizeDepParams.ldeb - 0.5 * (V3 + V4);
      T1 = 0.5 * (1.0 + V3 / V4);
      dTcen_dVg *= T1;
      dTcen_dVb *= T1;

      Ccen = CONSTEPSSI / Tcen;
      T2 = Cox / (Cox + Ccen);
      Coxeff = T2 * Ccen;
      T3 = -Ccen / Tcen;
      dCoxeff_dVg = T2 * T2 * T3;
      dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
      dCoxeff_dVg *= dTcen_dVg;
      CoxWLcen = CoxWL * Coxeff / Cox;

      Qac0 = CoxWLcen * (Vfbeff - sizeDepParams.vfbzb);
      QovCox = Qac0 / Coxeff;
      dQac0_dVg = CoxWLcen * dVfbeff_dVg + QovCox * dCoxeff_dVg;
      dQac0_dVb = CoxWLcen * dVfbeff_dVb + QovCox * dCoxeff_dVb;

      T0 = 0.5 * sizeDepParams.k1ox;
      T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
      if (sizeDepParams.k1ox == 0.0)
      {
        T1 = 0.0;
        T2 = 0.0;
      }
      else if (T3 < 0.0)
      {
        T1 = T0 + T3 / sizeDepParams.k1ox;
        T2 = CoxWLcen;
      }
      else
      {
        T1 = sqrt(T0 * T0 + T3);
        T2 = CoxWLcen * T0 / T1;
      }

      Qsub0 = CoxWLcen * sizeDepParams.k1ox * (T1 - T0);
      QovCox = Qsub0 / Coxeff;
      dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg)
        + QovCox * dCoxeff_dVg;
      dQsub0_dVd = -T2 * dVgsteff_dVd;
      dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb + dVgsteff_dVb)
        + QovCox * dCoxeff_dVb;

      // Gate-bias dependent delta Phis begins */
      if (sizeDepParams.k1ox <= 0.0)
      {
        Denomi = 0.25 * sizeDepParams.moin * Vtm;
        T0 = 0.5 * sizeDepParams.sqrtPhi;
      }
      else
      {
        Denomi = sizeDepParams.moin * Vtm * sizeDepParams.k1ox * sizeDepParams.k1ox;
        T0 = sizeDepParams.k1ox * sizeDepParams.sqrtPhi;
      }
      T1 = 2.0 * T0 + Vgsteff;

      DeltaPhi = Vtm * log(1.0 + T1 * Vgsteff / Denomi);
      dDeltaPhi_dVg = 2.0 * Vtm * (T1 -T0) / (Denomi + T1 * Vgsteff);
      dDeltaPhi_dVd = dDeltaPhi_dVg * dVgsteff_dVd;
      dDeltaPhi_dVb = dDeltaPhi_dVg * dVgsteff_dVb;
      // End of delta Phis

      T3 = 4.0 * (Vth - sizeDepParams.vfbzb - sizeDepParams.phi);
      Tox += Tox;
      if (T3 >= 0.0)
      {  T0 = (Vgsteff + T3) / Tox;
        dT0_dVd = (dVgsteff_dVd + 4.0 * dVth_dVd) / Tox;
        dT0_dVb = (dVgsteff_dVb + 4.0 * dVth_dVb) / Tox;
      }
      else
      {  T0 = (Vgsteff + 1.0e-20) / Tox;
        dT0_dVd = dVgsteff_dVd / Tox;
        dT0_dVb = dVgsteff_dVb / Tox;
      }
      tmp = exp(0.7 * log(T0));
      T1 = 1.0 + tmp;
      T2 = 0.7 * tmp / (T0 * Tox);
      Tcen = 1.9e-9 / T1;
      dTcen_dVg = -1.9e-9 * T2 / T1 /T1;
      dTcen_dVd = Tox * dTcen_dVg;
      dTcen_dVb = dTcen_dVd * dT0_dVb;
      dTcen_dVd *= dT0_dVd;
      dTcen_dVg *= dVgsteff_dVg;

      Ccen = CONSTEPSSI / Tcen;
      T0 = Cox / (Cox + Ccen);
      Coxeff = T0 * Ccen;
      T1 = -Ccen / Tcen;
      dCoxeff_dVg = T0 * T0 * T1;
      dCoxeff_dVd = dCoxeff_dVg * dTcen_dVd;
      dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
      dCoxeff_dVg *= dTcen_dVg;
      CoxWLcen = CoxWL * Coxeff / Cox;

      AbulkCV = Abulk0 * sizeDepParams.abulkCVfactor;
      dAbulkCV_dVb = sizeDepParams.abulkCVfactor * dAbulk0_dVb;
      VdsatCV = (Vgsteff - DeltaPhi) / AbulkCV;
      V4 = VdsatCV - Vds - CONSTDELTA_4;
      T0 = sqrt(V4 * V4 + 4.0 * CONSTDELTA_4 * VdsatCV);
      VdseffCV = VdsatCV - 0.5 * (V4 + T0);
      T1 = 0.5 * (1.0 + V4 / T0);
      T2 = CONSTDELTA_4 / T0;
      T3 = (1.0 - T1 - T2) / AbulkCV;
      T4 = T3 * ( 1.0 - dDeltaPhi_dVg);
      dVdseffCV_dVg = T4;
      dVdseffCV_dVd = T1;
      dVdseffCV_dVb = -T3 * VdsatCV * dAbulkCV_dVb;

      // Added to eliminate non-zero VdseffCV at Vds=0.0
      if (Vds == 0.0)
      {  VdseffCV = 0.0;
        dVdseffCV_dVg = 0.0;
        dVdseffCV_dVb = 0.0;
      }

      T0 = AbulkCV * VdseffCV;
      T1 = Vgsteff - DeltaPhi;
      T2 = 12.0 * (T1 - 0.5 * T0 + 1.0e-20);
      T3 = T0 / T2;
      T4 = 1.0 - 12.0 * T3 * T3;
      T5 = AbulkCV * (6.0 * T0 * (4.0 * T1 - T0) / (T2 * T2) - 0.5);
      T6 = T5 * VdseffCV / AbulkCV;

      qgate = qinoi = CoxWLcen * (T1 - T0 * (0.5 - T3));
      QovCox = qgate / Coxeff;
      Cgg1 = CoxWLcen * (T4 * (1.0 - dDeltaPhi_dVg) + T5 * dVdseffCV_dVg);
      Cgd1 = CoxWLcen * T5 * dVdseffCV_dVd + Cgg1
        * dVgsteff_dVd + QovCox * dCoxeff_dVd;
      Cgb1 = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
        + Cgg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
      Cgg1 = Cgg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;


      T7 = 1.0 - AbulkCV;
      T8 = T2 * T2;
      T9 = 12.0 * T7 * T0 * T0 / (T8 * AbulkCV);
      T10 = T9 * (1.0 - dDeltaPhi_dVg);
      T11 = -T7 * T5 / AbulkCV;
      T12 = -(T9 * T1 / AbulkCV + VdseffCV * (0.5 - T0 / T2));

      qbulk = CoxWLcen * T7 * (0.5 * VdseffCV - T0 * VdseffCV / T2);
      QovCox = qbulk / Coxeff;
      Cbg1 = CoxWLcen * (T10 + T11 * dVdseffCV_dVg);
      Cbd1 = CoxWLcen * T11 * dVdseffCV_dVd + Cbg1
        * dVgsteff_dVd + QovCox * dCoxeff_dVd;
      Cbb1 = CoxWLcen * (T11 * dVdseffCV_dVb + T12 * dAbulkCV_dVb)
        + Cbg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
      Cbg1 = Cbg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;

      if (model_xpart > 0.5)
      {   // 0/100 partition
        qsrc = -CoxWLcen * (T1 / 2.0 + T0 / 4.0 - 0.5 * T0 * T0 / T2);
        QovCox = qsrc / Coxeff;
        T2 += T2;
        T3 = T2 * T2;
        T7 = -(0.25 - 12.0 * T0 * (4.0 * T1 - T0) / T3);
        T4 = -(0.5 + 24.0 * T0 * T0 / T3) * (1.0 - dDeltaPhi_dVg);
        T5 = T7 * AbulkCV;
        T6 = T7 * VdseffCV;

        Csg = CoxWLcen * (T4 + T5 * dVdseffCV_dVg);
        Csd = CoxWLcen * T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd
          + QovCox * dCoxeff_dVd;
        Csb = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
          + Csg * dVgsteff_dVb + QovCox * dCoxeff_dVb;
        Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
      }
      else if (model_xpart < 0.5)
      {   // 40/60 partition
        T2 = T2 / 12.0;
        T3 = 0.5 * CoxWLcen / (T2 * T2);
        T4 = T1 * (2.0 * T0 * T0 / 3.0 + T1 * (T1 - 4.0
              * T0 / 3.0)) - 2.0 * T0 * T0 * T0 / 15.0;
        qsrc = -T3 * T4;
        QovCox = qsrc / Coxeff;
        T8 = 4.0 / 3.0 * T1 * (T1 - T0) + 0.4 * T0 * T0;
        T5 = -2.0 * qsrc / T2 - T3 * (T1 * (3.0 * T1 - 8.0
              * T0 / 3.0) + 2.0 * T0 * T0 / 3.0);
        T6 = AbulkCV * (qsrc / T2 + T3 * T8);
        T7 = T6 * VdseffCV / AbulkCV;

        Csg = T5 * (1.0 - dDeltaPhi_dVg) + T6 * dVdseffCV_dVg;
        Csd = Csg * dVgsteff_dVd + T6 * dVdseffCV_dVd
          + QovCox * dCoxeff_dVd;
        Csb = Csg * dVgsteff_dVb + T6 * dVdseffCV_dVb
          + T7 * dAbulkCV_dVb + QovCox * dCoxeff_dVb;
        Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
      }
      else
      {   // 50/50 partition
        qsrc = -0.5 * qgate;
        Csg = -0.5 * Cgg1;
        Csd = -0.5 * Cgd1;
        Csb = -0.5 * Cgb1;
      }

      qgate += Qac0 + Qsub0 - qbulk;
      qbulk -= (Qac0 + Qsub0);
      qdrn = -(qgate + qbulk + qsrc);

      Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
      Cbd = Cbd1 - dQsub0_dVd;
      Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;

      Cgg = Cgg1 - Cbg;
      Cgd = Cgd1 - Cbd;
      Cgb = Cgb1 - Cbb;

      Cgb *= dVbseff_dVb;
      Cbb *= dVbseff_dVb;
      Csb *= dVbseff_dVb;

      cggb = Cgg;
      cgsb = -(Cgg + Cgd + Cgb);
      cgdb = Cgd;
      cdgb = -(Cgg + Cbg + Csg);
      cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb + Csg + Csd + Csb);
      cddb = -(Cgd + Cbd + Csd);
      cbgb = Cbg;
      cbsb = -(Cbg + Cbd + Cbb);
      cbdb = Cbd;
      qinv = -qinoi;
    }  // End of CTM
  }

finished:
  // Returning Values to Calling Routine
  // COMPUTE EQUIVALENT DRAIN CURRENT SOURCE

  // copy local "cdrain" variable over to instance variable cd.
  cd = cdrain;

  // charge storage elements:
  // bulk-drain and bulk-source depletion capacitances
  // czbd : zero bias drain junction capacitance
  // czbs : zero bias source junction capacitance
  // czbdsw: zero bias drain junction sidewall capacitance
  //         along field oxide
  // czbssw: zero bias source junction sidewall capacitance
  //         along field oxide
  // czbdswg: zero bias drain junction sidewall capacitance
  //          along gate side
  // czbsswg: zero bias source junction sidewall capacitance
  //          along gate side
  if (ChargeComputationNeeded)
  {
    czbd = unitAreaJctCapTemp * drainArea;
    czbs = unitAreaJctCapTemp * sourceArea;
    if (drainPerimeter < sizeDepParams.weff)
    {
      czbdswg = unitLengthGateSidewallJctCapTemp * drainPerimeter;
      czbdsw = 0.0;
    }
    else
    {
      czbdsw = unitLengthSidewallJctCapTemp
        * (drainPerimeter - sizeDepParams.weff);

      czbdswg = unitLengthGateSidewallJctCapTemp *  sizeDepParams.weff;
    }
    if (sourcePerimeter < sizeDepParams.weff)
    {
      czbssw = 0.0;
      czbsswg = unitLengthGateSidewallJctCapTemp * sourcePerimeter;
    }
    else
    {
      czbssw = unitLengthSidewallJctCapTemp
        * (sourcePerimeter - sizeDepParams.weff);
      czbsswg = unitLengthGateSidewallJctCapTemp *  sizeDepParams.weff;
    }

    MJ    = model_bulkJctBotGradingCoeff;
    MJSW  = model_bulkJctSideGradingCoeff;
    MJSWG = model_bulkJctGateSideGradingCoeff;

    // Source Bulk Junction
    if (vbs == 0.0)
    {
      //*(ckt->CKTstate0 + iterI->qbs) = 0.0;
      qbs = 0.0;
      capbs = czbs + czbssw + czbsswg;
    }
    else if (vbs < 0.0)
    {
      if (czbs > 0.0)
      {
        arg = 1.0 - vbs / PhiBTemp;

        if (MJ == 0.5) sarg = 1.0 / sqrt(arg);
        else           sarg = exp(-MJ * log(arg));

        //*(ckt->CKTstate0 + iterI->qbs) =
        qbs = PhiBTemp * czbs * (1.0 - arg * sarg) / (1.0 - MJ);

        capbs = czbs * sarg;
      }
      else
      {
        //*(ckt->CKTstate0 + iterI->qbs) = 0.0;
        qbs = 0.0;
        capbs = 0.0;
      }

      if (czbssw > 0.0)
      {
        arg = 1.0 - vbs / PhiBSWTemp;
        if (MJSW == 0.5) sarg = 1.0 / sqrt(arg);
        else             sarg = exp(-MJSW * log(arg));

        //*(ckt->CKTstate0 + iterI->qbs) +=
        qbs += PhiBSWTemp * czbssw * (1.0 - arg * sarg) / (1.0 - MJSW);

        capbs += czbssw * sarg;
      }

      if (czbsswg > 0.0)
      {
        arg = 1.0 - vbs / PhiBSWGTemp;
        if (MJSWG == 0.5) sarg = 1.0 / sqrt(arg);
        else              sarg = exp(-MJSWG * log(arg));

        //*(ckt->CKTstate0 + iterI->qbs) +=
        qbs += PhiBSWGTemp * czbsswg * (1.0 - arg * sarg) / (1.0 - MJSWG);

        capbs += czbsswg * sarg;
      }

    }
    else
    {
      T0 = czbs + czbssw + czbsswg;
      T1 = vbs * (czbs * MJ / PhiBTemp + czbssw * MJSW
          / PhiBSWTemp + czbsswg * MJSWG / PhiBSWGTemp);

      //*(ckt->CKTstate0 + iterI->
      qbs = vbs * (T0 + 0.5 * T1);
      capbs = T0 + T1;
    }

    // Drain Bulk Junction
    if (vbd == 0.0)
    {
      //*(ckt->CKTstate0 + iterI->qbd) = 0.0;
      qbd = 0.0;
      capbd = czbd + czbdsw + czbdswg;
    }
    else if (vbd < 0.0)
    {
      if (czbd > 0.0)
      {
        arg = 1.0 - vbd / PhiBTemp;
        if (MJ == 0.5) sarg = 1.0 / sqrt(arg);
        else           sarg = exp(-MJ * log(arg));

        //*(ckt->CKTstate0 + iterI->qbd) =
        qbd = PhiBTemp * czbd * (1.0 - arg * sarg) / (1.0 - MJ);
        capbd = czbd * sarg;
      }
      else
      {
        //*(ckt->CKTstate0 + iterI->qbd) = 0.0;
        qbd = 0.0;
        capbd = 0.0;
      }

      if (czbdsw > 0.0)
      {
        arg = 1.0 - vbd / PhiBSWTemp;
        if (MJSW == 0.5) sarg = 1.0 / sqrt(arg);
        else             sarg = exp(-MJSW * log(arg));

        //*(ckt->CKTstate0 + iterI->qbd) +=
        qbd += PhiBSWTemp * czbdsw * (1.0 - arg * sarg) / (1.0 - MJSW);
        capbd += czbdsw * sarg;
      }

      if (czbdswg > 0.0)
      {
        arg = 1.0 - vbd / PhiBSWGTemp;
        if (MJSWG == 0.5) sarg = 1.0 / sqrt(arg);
        else              sarg = exp(-MJSWG * log(arg));

        //*(ckt->CKTstate0 + iterI->qbd) +=
        qbd += PhiBSWGTemp * czbdswg * (1.0 - arg * sarg) / (1.0 - MJSWG);
        capbd += czbdswg * sarg;
      }
    }
    else
    {
      T0 = czbd + czbdsw + czbdswg;
      T1 = vbd * (czbd * MJ / PhiBTemp + czbdsw * MJSW
          / PhiBSWTemp + czbdswg * MJSWG / PhiBSWGTemp);

      //*(ckt->CKTstate0 + iterI->qbd) = vbd * (T0 + 0.5 * T1);
      qbd = vbd * (T0 + 0.5 * T1);
      capbd = T0 + T1;
    }
  }

  // There is a spice3f5 convergence check that would happen here.
  // (line 2404) skipping...

  // In 3f5, loading a bunch of things into the state vector at this point.
  // (line 2433) skipping...

  // bulk and channel charge plus overlaps

  if (ChargeComputationNeeded)
  {
    // NQS begins
    if (nqsMod)
    {
      qcheq = -(qbulk + qgate);

      cqgb = -(cggb + cbgb);
      cqdb = -(cgdb + cbdb);
      cqsb = -(cgsb + cbsb);
      cqbb = -(cqgb + cqdb + cqsb);

      gtau_drift = fabs(sizeDepParams.tconst * qcheq) * ScalingFactor;
      T0 = sizeDepParams.leffCV * sizeDepParams.leffCV;
      gtau_diff = 16.0 * sizeDepParams.u0temp * model_vtm / T0 * ScalingFactor;

      gtau =  gtau_drift + gtau_diff;
    }

    if (model_capMod == 0)
    {
      if (vgd < 0.0)
      {
        cgdo = sizeDepParams.cgdo;
        qgdo = sizeDepParams.cgdo * vgd;
      }
      else
      {
        cgdo = sizeDepParams.cgdo;
        qgdo =  sizeDepParams.cgdo * vgd;
      }

      if (vgs < 0.0)
      {
        cgso = sizeDepParams.cgso;
        qgso = sizeDepParams.cgso * vgs;
      }
      else
      {
        cgso = sizeDepParams.cgso;
        qgso =  sizeDepParams.cgso * vgs;
      }
    }
    else if (model_capMod == 1)
    {
      if (vgd < 0.0)
      {
        T1 = sqrt(1.0 - 4.0 * vgd / sizeDepParams.ckappa);
        cgdo = sizeDepParams.cgdo + sizeDepParams.weffCV * sizeDepParams.cgdl / T1;

        qgdo = sizeDepParams.cgdo * vgd - sizeDepParams.weffCV * 0.5
          * sizeDepParams.cgdl * sizeDepParams.ckappa * (T1 - 1.0);
      }
      else
      {
        cgdo = sizeDepParams.cgdo + sizeDepParams.weffCV * sizeDepParams.cgdl;
        qgdo = (sizeDepParams.weffCV * sizeDepParams.cgdl + sizeDepParams.cgdo) * vgd;
      }

      if (vgs < 0.0)
      {
        T1 = sqrt(1.0 - 4.0 * vgs / sizeDepParams.ckappa);
        cgso = sizeDepParams.cgso + sizeDepParams.weffCV * sizeDepParams.cgsl / T1;
        qgso = sizeDepParams.cgso * vgs - sizeDepParams.weffCV * 0.5
          * sizeDepParams.cgsl * sizeDepParams.ckappa * (T1 - 1.0);
      }
      else
      {
        cgso = sizeDepParams.cgso + sizeDepParams.weffCV * sizeDepParams.cgsl;
        qgso = (sizeDepParams.weffCV * sizeDepParams.cgsl + sizeDepParams.cgso) * vgs;
      }
    }
    else
    {
      T0 = vgd + CONSTDELTA_1;
      T1 = sqrt(T0 * T0 + 4.0 * CONSTDELTA_1);
      T2 = 0.5 * (T0 - T1);

      T3 = sizeDepParams.weffCV * sizeDepParams.cgdl;
      T4 = sqrt(1.0 - 4.0 * T2 / sizeDepParams.ckappa);
      cgdo = sizeDepParams.cgdo + T3 - T3 * (1.0 - 1.0 / T4) * (0.5 - 0.5 * T0 / T1);

      qgdo = (sizeDepParams.cgdo + T3) * vgd - T3 * (T2
          + 0.5 * sizeDepParams.ckappa * (T4 - 1.0));

      T0 = vgs + CONSTDELTA_1;
      T1 = sqrt(T0 * T0 + 4.0 * CONSTDELTA_1);
      T2 = 0.5 * (T0 - T1);
      T3 = sizeDepParams.weffCV * sizeDepParams.cgsl;
      T4 = sqrt(1.0 - 4.0 * T2 / sizeDepParams.ckappa);
      cgso = sizeDepParams.cgso + T3 - T3 * (1.0 - 1.0 / T4) * (0.5 - 0.5 * T0 / T1);
      qgso = (sizeDepParams.cgso + T3) * vgs - T3 * (T2
          + 0.5 * sizeDepParams.ckappa * (T4 - 1.0));
    }

    setupCapacitors_oldDAE (
        // inputs
  solverState,
  mode,
  nqsMod,
  sizeDepParams,
    cggb,
    cgdo,
    cgso,
    cgdb,
    cgsb,
    cdgb,
    cddb,
    capbd,
    cdsb,
    cbgb,
    cbdb,
    capbs,
    cbsb,
    qgdo,
    qgso,
    vgb,
    qcheq,
    qdef,
    ScalingFactor,
   cqgb,
   cqdb,
   cqsb,
   cqbb,
   model_cox,
   model_xpart,
      // outputs

       qgate,
       qbulk,
       qdrn,
       qsrc,
       CoxWL
  );

    setupCapacitors_newDAE (
      mode, nqsMod, sizeDepParams,
      cggb, cgdo, cgso, cgdb, cgsb, cdgb, cddb,
      capbd, cdsb, cbgb, cbdb, capbs, cbsb,
      //outputs
      CAPcggb, CAPcgdb, CAPcgsb, CAPcdgb, CAPcddb,
      CAPcdsb, CAPcsgb, CAPcsdb, CAPcssb, CAPcbgb,
      CAPcbdb, CAPcbsb);

    cqdef = 0.0;

    // set some state variables:
    qg = qgate;
    qd = qdrn - qbd;
    qb = qbulk + qbd + qbs;
    if (nqsMod) qcdump = qdef * ScalingFactor;

  } // end of ChargeComputationNeeded if statement.

  // store small signal parameters
  //if (ckt->CKTmode & MODEINITSMSIG) goto line1000;
  // Note: in 3f5, line1000 is at the end of the load, after
  //        the loads to the rhs and the matrix.  So it looks
  //        like this goto essentially means return.


  // Setting up a few currents for the RHS load:
  Idrain = drainConductance * Vddp;
  Isource = sourceConductance * Vssp;

  // Put this kludge in because the matrix load needs T1 but it is used
  // all over the place:
  T1global = T1;

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : auxChargeCalculations
//
// Purpose       : This function does some final "cleanup" calculations
//                 having to do with the capacitors.
//
// Special Notes : About all this function really does is set up some
//                 voltlim terms, and some unused nqs stuff.  
//                 For sensitivities, this is useless ; but I implemented b4 
//                 I realized that.  If the origFlag = true, then limiting is not
//                 applied.  So most of the below code is then a no-op anyway.
//
//                 All this does at this point is zero out gtau.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/30/2018
//-----------------------------------------------------------------------------
template <typename ScalarT> 
bool auxChargeCalculations (
       bool ChargeComputationNeeded,
       int mode,
       int nqsMod,
       const ScalarT & model_vtm,
       const SizeDependParam<ScalarT> & sizeDepParams,
    const ScalarT & ScalingFactor,
       ScalarT & gtau
    )
{
  ScalarT T0, T1;

  if (!ChargeComputationNeeded)
  {
    if (nqsMod) // useless as nqsMod isn't implemented
    {
      gtau = 16.0 * sizeDepParams.u0temp * model_vtm
             / sizeDepParams.leffCV / sizeDepParams.leffCV * ScalingFactor;
    }
    else
    {
      gtau = 0.0;
    }
  }
  else  // ChargeComputation is needed
  {
    // everything that was in here in the original function was voltlim related.
    // voltlim isn't used during senstivity calculation, so not relevant, and deleted.

    // Note:  nqs stuff is not finished, so excluded here.

  } // !ChargeComputationNeeded

  return true;
}

//-----------------------------------------------------------------------------
// Function      : bsim3InstanceSensitivity::operator
// Purpose       : produces df/dp and dq/dp, where p= any bsim3 instance parameter.  
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/20/2018
//-----------------------------------------------------------------------------
void bsim3InstanceSensitivity::operator()(
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
  const Instance & inst = *(dynamic_cast<const Instance *> (e1));
  const Xyce::Device::MOSFET_B3::Model & mod = inst.model_;

  dfdp.resize(10);
  dqdp.resize(10);
  Findices.resize(10);
  Qindices.resize(10);

  // model params
  // sensitivity supported, (doubles):
  fadType modelPar_tox=mod.tox;	  bool modelPar_given_TOX=mod.given("TOX");
  fadType modelPar_toxm=mod.toxm;	  bool modelPar_given_TOXM=mod.given("TOXM");
  fadType modelPar_cdsc=mod.cdsc;	  bool modelPar_given_CDSC=mod.given("CDSC");
  fadType modelPar_cdscb=mod.cdscb;	  bool modelPar_given_CDSCB=mod.given("CDSCB");
  fadType modelPar_cdscd=mod.cdscd;	  bool modelPar_given_CDSCD=mod.given("CDSCD");
  fadType modelPar_cit=mod.cit;	  bool modelPar_given_CIT=mod.given("CIT");
  fadType modelPar_nfactor=mod.nfactor;	  bool modelPar_given_NFACTOR=mod.given("NFACTOR");
  fadType modelPar_xj=mod.xj;	  bool modelPar_given_XJ=mod.given("XJ");
  fadType modelPar_vsat=mod.vsat;	  bool modelPar_given_VSAT=mod.given("VSAT");
  fadType modelPar_at=mod.at;	  bool modelPar_given_AT=mod.given("AT");
  fadType modelPar_a0=mod.a0;	  bool modelPar_given_A0=mod.given("A0");
  fadType modelPar_ags=mod.ags;	  bool modelPar_given_AGS=mod.given("AGS");
  fadType modelPar_a1=mod.a1;	  bool modelPar_given_A1=mod.given("A1");
  fadType modelPar_a2=mod.a2;	  bool modelPar_given_A2=mod.given("A2");
  fadType modelPar_keta=mod.keta;	  bool modelPar_given_KETA=mod.given("KETA");
  fadType modelPar_nsub=mod.nsub;	  bool modelPar_given_NSUB=mod.given("NSUB");
  fadType modelPar_npeak=mod.npeak;	  bool modelPar_given_NCH=mod.given("NCH");
  fadType modelPar_ngate=mod.ngate;	  bool modelPar_given_NGATE=mod.given("NGATE");
  fadType modelPar_gamma1=mod.gamma1;	  bool modelPar_given_GAMMA1=mod.given("GAMMA1");
  fadType modelPar_gamma2=mod.gamma2;	  bool modelPar_given_GAMMA2=mod.given("GAMMA2");
  fadType modelPar_vbx=mod.vbx;	  bool modelPar_given_VBX=mod.given("VBX");
  fadType modelPar_vbm=mod.vbm;	  bool modelPar_given_VBM=mod.given("VBM");
  fadType modelPar_xt=mod.xt;	  bool modelPar_given_XT=mod.given("XT");
  fadType modelPar_k1=mod.k1;	  bool modelPar_given_K1=mod.given("K1");
  fadType modelPar_kt1=mod.kt1;	  bool modelPar_given_KT1=mod.given("KT1");
  fadType modelPar_kt1l=mod.kt1l;	  bool modelPar_given_KT1L=mod.given("KT1L");
  fadType modelPar_kt2=mod.kt2;	  bool modelPar_given_KT2=mod.given("KT2");
  fadType modelPar_k2=mod.k2;	  bool modelPar_given_K2=mod.given("K2");
  fadType modelPar_k3=mod.k3;	  bool modelPar_given_K3=mod.given("K3");
  fadType modelPar_k3b=mod.k3b;	  bool modelPar_given_K3B=mod.given("K3B");
  fadType modelPar_w0=mod.w0;	  bool modelPar_given_W0=mod.given("W0");
  fadType modelPar_nlx=mod.nlx;	  bool modelPar_given_NLX=mod.given("NLX");
  fadType modelPar_dvt0=mod.dvt0;	  bool modelPar_given_DVT0=mod.given("DVT0");
  fadType modelPar_dvt1=mod.dvt1;	  bool modelPar_given_DVT1=mod.given("DVT1");
  fadType modelPar_dvt2=mod.dvt2;	  bool modelPar_given_DVT2=mod.given("DVT2");
  fadType modelPar_dvt0w=mod.dvt0w;	  bool modelPar_given_DVT0W=mod.given("DVT0W");
  fadType modelPar_dvt1w=mod.dvt1w;	  bool modelPar_given_DVT1W=mod.given("DVT1W");
  fadType modelPar_dvt2w=mod.dvt2w;	  bool modelPar_given_DVT2W=mod.given("DVT2W");
  fadType modelPar_drout=mod.drout;	  bool modelPar_given_DROUT=mod.given("DROUT");
  fadType modelPar_dsub=mod.dsub;	  bool modelPar_given_DSUB=mod.given("DSUB");
  fadType modelPar_vth0=mod.vth0;	  bool modelPar_given_VTH0=mod.given("VTH0");
  fadType modelPar_ua=mod.ua;	  bool modelPar_given_UA=mod.given("UA");
  fadType modelPar_ua1=mod.ua1;	  bool modelPar_given_UA1=mod.given("UA1");
  fadType modelPar_ub=mod.ub;	  bool modelPar_given_UB=mod.given("UB");
  fadType modelPar_ub1=mod.ub1;	  bool modelPar_given_UB1=mod.given("UB1");
  fadType modelPar_uc=mod.uc;	  bool modelPar_given_UC=mod.given("UC");
  fadType modelPar_uc1=mod.uc1;	  bool modelPar_given_UC1=mod.given("UC1");
  fadType modelPar_u0=mod.u0;	  bool modelPar_given_U0=mod.given("U0");
  fadType modelPar_ute=mod.ute;	  bool modelPar_given_UTE=mod.given("UTE");
  fadType modelPar_voff=mod.voff;	  bool modelPar_given_VOFF=mod.given("VOFF");
  fadType modelPar_rdsw=mod.rdsw;	  bool modelPar_given_RDSW=mod.given("RDSW");
  fadType modelPar_prwg=mod.prwg;	  bool modelPar_given_PRWG=mod.given("PRWG");
  fadType modelPar_prwb=mod.prwb;	  bool modelPar_given_PRWB=mod.given("PRWB");
  fadType modelPar_prt=mod.prt;	  bool modelPar_given_PRT=mod.given("PRT");
  fadType modelPar_eta0=mod.eta0;	  bool modelPar_given_ETA0=mod.given("ETA0");
  fadType modelPar_etab=mod.etab;	  bool modelPar_given_ETAB=mod.given("ETAB");
  fadType modelPar_pclm=mod.pclm;	  bool modelPar_given_PCLM=mod.given("PCLM");
  fadType modelPar_pdibl1=mod.pdibl1;	  bool modelPar_given_PDIBLC1=mod.given("PDIBLC1");
  fadType modelPar_pdibl2=mod.pdibl2;	  bool modelPar_given_PDIBLC2=mod.given("PDIBLC2");
  fadType modelPar_pdiblb=mod.pdiblb;	  bool modelPar_given_PDIBLCB=mod.given("PDIBLCB");
  fadType modelPar_pscbe1=mod.pscbe1;	  bool modelPar_given_PSCBE1=mod.given("PSCBE1");
  fadType modelPar_pscbe2=mod.pscbe2;	  bool modelPar_given_PSCBE2=mod.given("PSCBE2");
  fadType modelPar_pvag=mod.pvag;	  bool modelPar_given_PVAG=mod.given("PVAG");
  fadType modelPar_delta=mod.delta;	  bool modelPar_given_DELTA=mod.given("DELTA");
  fadType modelPar_wr=mod.wr;	  bool modelPar_given_WR=mod.given("WR");
  fadType modelPar_dwg=mod.dwg;	  bool modelPar_given_DWG=mod.given("DWG");
  fadType modelPar_dwb=mod.dwb;	  bool modelPar_given_DWB=mod.given("DWB");
  fadType modelPar_b0=mod.b0;	  bool modelPar_given_B0=mod.given("B0");
  fadType modelPar_b1=mod.b1;	  bool modelPar_given_B1=mod.given("B1");
  fadType modelPar_alpha0=mod.alpha0;	  bool modelPar_given_ALPHA0=mod.given("ALPHA0");
  fadType modelPar_alpha1=mod.alpha1;	  bool modelPar_given_ALPHA1=mod.given("ALPHA1");
  fadType modelPar_beta0=mod.beta0;	  bool modelPar_given_BETA0=mod.given("BETA0");
  fadType modelPar_ijth=mod.ijth;	  bool modelPar_given_IJTH=mod.given("IJTH");
  fadType modelPar_vfb=mod.vfb;	  bool modelPar_given_VFB=mod.given("VFB");
  fadType modelPar_elm=mod.elm;	  bool modelPar_given_ELM=mod.given("ELM");
  fadType modelPar_cgsl=mod.cgsl;	  bool modelPar_given_CGSL=mod.given("CGSL");
  fadType modelPar_cgdl=mod.cgdl;	  bool modelPar_given_CGDL=mod.given("CGDL");
  fadType modelPar_ckappa=mod.ckappa;	  bool modelPar_given_CKAPPA=mod.given("CKAPPA");
  fadType modelPar_cf=mod.cf;	  bool modelPar_given_CF=mod.given("CF");
  fadType modelPar_vfbcv=mod.vfbcv;	  bool modelPar_given_VFBCV=mod.given("VFBCV");
  fadType modelPar_clc=mod.clc;	  bool modelPar_given_CLC=mod.given("CLC");
  fadType modelPar_cle=mod.cle;	  bool modelPar_given_CLE=mod.given("CLE");
  fadType modelPar_dwc=mod.dwc;	  bool modelPar_given_DWC=mod.given("DWC");
  fadType modelPar_dlc=mod.dlc;	  bool modelPar_given_DLC=mod.given("DLC");
  fadType modelPar_noff=mod.noff;	  bool modelPar_given_NOFF=mod.given("NOFF");
  fadType modelPar_voffcv=mod.voffcv;	  bool modelPar_given_VOFFCV=mod.given("VOFFCV");
  fadType modelPar_acde=mod.acde;	  bool modelPar_given_ACDE=mod.given("ACDE");
  fadType modelPar_moin=mod.moin;	  bool modelPar_given_MOIN=mod.given("MOIN");
  fadType modelPar_tcj=mod.tcj;	  bool modelPar_given_TCJ=mod.given("TCJ");
  fadType modelPar_tcjsw=mod.tcjsw;	  bool modelPar_given_TCJSW=mod.given("TCJSW");
  fadType modelPar_tcjswg=mod.tcjswg;	  bool modelPar_given_TCJSWG=mod.given("TCJSWG");
  fadType modelPar_tpb=mod.tpb;	  bool modelPar_given_TPB=mod.given("TPB");
  fadType modelPar_tpbsw=mod.tpbsw;	  bool modelPar_given_TPBSW=mod.given("TPBSW");
  fadType modelPar_tpbswg=mod.tpbswg;	  bool modelPar_given_TPBSWG=mod.given("TPBSWG");
  fadType modelPar_lcdsc=mod.lcdsc;	  bool modelPar_given_LCDSC=mod.given("LCDSC");
  fadType modelPar_lcdscb=mod.lcdscb;	  bool modelPar_given_LCDSCB=mod.given("LCDSCB");
  fadType modelPar_lcdscd=mod.lcdscd;	  bool modelPar_given_LCDSCD=mod.given("LCDSCD");
  fadType modelPar_lcit=mod.lcit;	  bool modelPar_given_LCIT=mod.given("LCIT");
  fadType modelPar_lnfactor=mod.lnfactor;	  bool modelPar_given_LNFACTOR=mod.given("LNFACTOR");
  fadType modelPar_lxj=mod.lxj;	  bool modelPar_given_LXJ=mod.given("LXJ");
  fadType modelPar_lvsat=mod.lvsat;	  bool modelPar_given_LVSAT=mod.given("LVSAT");
  fadType modelPar_lat=mod.lat;	  bool modelPar_given_LAT=mod.given("LAT");
  fadType modelPar_la0=mod.la0;	  bool modelPar_given_LA0=mod.given("LA0");
  fadType modelPar_lags=mod.lags;	  bool modelPar_given_LAGS=mod.given("LAGS");
  fadType modelPar_la1=mod.la1;	  bool modelPar_given_LA1=mod.given("LA1");
  fadType modelPar_la2=mod.la2;	  bool modelPar_given_LA2=mod.given("LA2");
  fadType modelPar_lketa=mod.lketa;	  bool modelPar_given_LKETA=mod.given("LKETA");
  fadType modelPar_lnsub=mod.lnsub;	  bool modelPar_given_LNSUB=mod.given("LNSUB");
  fadType modelPar_lnpeak=mod.lnpeak;	  bool modelPar_given_LNCH=mod.given("LNCH");
  fadType modelPar_lngate=mod.lngate;	  bool modelPar_given_LNGATE=mod.given("LNGATE");
  fadType modelPar_lgamma1=mod.lgamma1;	  bool modelPar_given_LGAMMA1=mod.given("LGAMMA1");
  fadType modelPar_lgamma2=mod.lgamma2;	  bool modelPar_given_LGAMMA2=mod.given("LGAMMA2");
  fadType modelPar_lvbx=mod.lvbx;	  bool modelPar_given_LVBX=mod.given("LVBX");
  fadType modelPar_lvbm=mod.lvbm;	  bool modelPar_given_LVBM=mod.given("LVBM");
  fadType modelPar_lxt=mod.lxt;	  bool modelPar_given_LXT=mod.given("LXT");
  fadType modelPar_lk1=mod.lk1;	  bool modelPar_given_LK1=mod.given("LK1");
  fadType modelPar_lkt1=mod.lkt1;	  bool modelPar_given_LKT1=mod.given("LKT1");
  fadType modelPar_lkt1l=mod.lkt1l;	  bool modelPar_given_LKT1L=mod.given("LKT1L");
  fadType modelPar_lkt2=mod.lkt2;	  bool modelPar_given_LKT2=mod.given("LKT2");
  fadType modelPar_lk2=mod.lk2;	  bool modelPar_given_LK2=mod.given("LK2");
  fadType modelPar_lk3=mod.lk3;	  bool modelPar_given_LK3=mod.given("LK3");
  fadType modelPar_lk3b=mod.lk3b;	  bool modelPar_given_LK3B=mod.given("LK3B");
  fadType modelPar_lw0=mod.lw0;	  bool modelPar_given_LW0=mod.given("LW0");
  fadType modelPar_lnlx=mod.lnlx;	  bool modelPar_given_LNLX=mod.given("LNLX");
  fadType modelPar_ldvt0=mod.ldvt0;	  bool modelPar_given_LDVT0=mod.given("LDVT0");
  fadType modelPar_ldvt1=mod.ldvt1;	  bool modelPar_given_LDVT1=mod.given("LDVT1");
  fadType modelPar_ldvt2=mod.ldvt2;	  bool modelPar_given_LDVT2=mod.given("LDVT2");
  fadType modelPar_ldvt0w=mod.ldvt0w;	  bool modelPar_given_LDVT0W=mod.given("LDVT0W");
  fadType modelPar_ldvt1w=mod.ldvt1w;	  bool modelPar_given_LDVT1W=mod.given("LDVT1W");
  fadType modelPar_ldvt2w=mod.ldvt2w;	  bool modelPar_given_LDVT2W=mod.given("LDVT2W");
  fadType modelPar_ldrout=mod.ldrout;	  bool modelPar_given_LDROUT=mod.given("LDROUT");
  fadType modelPar_ldsub=mod.ldsub;	  bool modelPar_given_LDSUB=mod.given("LDSUB");
  fadType modelPar_lvth0=mod.lvth0;	  bool modelPar_given_LVTH0=mod.given("LVTH0");
  fadType modelPar_lua=mod.lua;	  bool modelPar_given_LUA=mod.given("LUA");
  fadType modelPar_lua1=mod.lua1;	  bool modelPar_given_LUA1=mod.given("LUA1");
  fadType modelPar_lub=mod.lub;	  bool modelPar_given_LUB=mod.given("LUB");
  fadType modelPar_lub1=mod.lub1;	  bool modelPar_given_LUB1=mod.given("LUB1");
  fadType modelPar_luc=mod.luc;	  bool modelPar_given_LUC=mod.given("LUC");
  fadType modelPar_luc1=mod.luc1;	  bool modelPar_given_LUC1=mod.given("LUC1");
  fadType modelPar_lu0=mod.lu0;	  bool modelPar_given_LU0=mod.given("LU0");
  fadType modelPar_lute=mod.lute;	  bool modelPar_given_LUTE=mod.given("LUTE");
  fadType modelPar_lvoff=mod.lvoff;	  bool modelPar_given_LVOFF=mod.given("LVOFF");
  fadType modelPar_lrdsw=mod.lrdsw;	  bool modelPar_given_LRDSW=mod.given("LRDSW");
  fadType modelPar_lprwg=mod.lprwg;	  bool modelPar_given_LPRWG=mod.given("LPRWG");
  fadType modelPar_lprwb=mod.lprwb;	  bool modelPar_given_LPRWB=mod.given("LPRWB");
  fadType modelPar_lprt=mod.lprt;	  bool modelPar_given_LPRT=mod.given("LPRT");
  fadType modelPar_leta0=mod.leta0;	  bool modelPar_given_LETA0=mod.given("LETA0");
  fadType modelPar_letab=mod.letab;	  bool modelPar_given_LETAB=mod.given("LETAB");
  fadType modelPar_lpclm=mod.lpclm;	  bool modelPar_given_LPCLM=mod.given("LPCLM");
  fadType modelPar_lpdibl1=mod.lpdibl1;	  bool modelPar_given_LPDIBLC1=mod.given("LPDIBLC1");
  fadType modelPar_lpdibl2=mod.lpdibl2;	  bool modelPar_given_LPDIBLC2=mod.given("LPDIBLC2");
  fadType modelPar_lpdiblb=mod.lpdiblb;	  bool modelPar_given_LPDIBLCB=mod.given("LPDIBLCB");
  fadType modelPar_lpscbe1=mod.lpscbe1;	  bool modelPar_given_LPSCBE1=mod.given("LPSCBE1");
  fadType modelPar_lpscbe2=mod.lpscbe2;	  bool modelPar_given_LPSCBE2=mod.given("LPSCBE2");
  fadType modelPar_lpvag=mod.lpvag;	  bool modelPar_given_LPVAG=mod.given("LPVAG");
  fadType modelPar_ldelta=mod.ldelta;	  bool modelPar_given_LDELTA=mod.given("LDELTA");
  fadType modelPar_lwr=mod.lwr;	  bool modelPar_given_LWR=mod.given("LWR");
  fadType modelPar_ldwg=mod.ldwg;	  bool modelPar_given_LDWG=mod.given("LDWG");
  fadType modelPar_ldwb=mod.ldwb;	  bool modelPar_given_LDWB=mod.given("LDWB");
  fadType modelPar_lb0=mod.lb0;	  bool modelPar_given_LB0=mod.given("LB0");
  fadType modelPar_lb1=mod.lb1;	  bool modelPar_given_LB1=mod.given("LB1");
  fadType modelPar_lalpha0=mod.lalpha0;	  bool modelPar_given_LALPHA0=mod.given("LALPHA0");
  fadType modelPar_lalpha1=mod.lalpha1;	  bool modelPar_given_LALPHA1=mod.given("LALPHA1");
  fadType modelPar_lbeta0=mod.lbeta0;	  bool modelPar_given_LBETA0=mod.given("LBETA0");
  fadType modelPar_lvfb=mod.lvfb;	  bool modelPar_given_LVFB=mod.given("LVFB");
  fadType modelPar_lelm=mod.lelm;	  bool modelPar_given_LELM=mod.given("LELM");
  fadType modelPar_lcgsl=mod.lcgsl;	  bool modelPar_given_LCGSL=mod.given("LCGSL");
  fadType modelPar_lcgdl=mod.lcgdl;	  bool modelPar_given_LCGDL=mod.given("LCGDL");
  fadType modelPar_lckappa=mod.lckappa;	  bool modelPar_given_LCKAPPA=mod.given("LCKAPPA");
  fadType modelPar_lcf=mod.lcf;	  bool modelPar_given_LCF=mod.given("LCF");
  fadType modelPar_lclc=mod.lclc;	  bool modelPar_given_LCLC=mod.given("LCLC");
  fadType modelPar_lcle=mod.lcle;	  bool modelPar_given_LCLE=mod.given("LCLE");
  fadType modelPar_lvfbcv=mod.lvfbcv;	  bool modelPar_given_LVFBCV=mod.given("LVFBCV");
  fadType modelPar_lnoff=mod.lnoff;	  bool modelPar_given_LNOFF=mod.given("LNOFF");
  fadType modelPar_lvoffcv=mod.lvoffcv;	  bool modelPar_given_LVOFFCV=mod.given("LVOFFCV");
  fadType modelPar_lacde=mod.lacde;	  bool modelPar_given_LACDE=mod.given("LACDE");
  fadType modelPar_lmoin=mod.lmoin;	  bool modelPar_given_LMOIN=mod.given("LMOIN");
  fadType modelPar_wcdsc=mod.wcdsc;	  bool modelPar_given_WCDSC=mod.given("WCDSC");
  fadType modelPar_wcdscb=mod.wcdscb;	  bool modelPar_given_WCDSCB=mod.given("WCDSCB");
  fadType modelPar_wcdscd=mod.wcdscd;	  bool modelPar_given_WCDSCD=mod.given("WCDSCD");
  fadType modelPar_wcit=mod.wcit;	  bool modelPar_given_WCIT=mod.given("WCIT");
  fadType modelPar_wnfactor=mod.wnfactor;	  bool modelPar_given_WNFACTOR=mod.given("WNFACTOR");
  fadType modelPar_wxj=mod.wxj;	  bool modelPar_given_WXJ=mod.given("WXJ");
  fadType modelPar_wvsat=mod.wvsat;	  bool modelPar_given_WVSAT=mod.given("WVSAT");
  fadType modelPar_wat=mod.wat;	  bool modelPar_given_WAT=mod.given("WAT");
  fadType modelPar_wa0=mod.wa0;	  bool modelPar_given_WA0=mod.given("WA0");
  fadType modelPar_wags=mod.wags;	  bool modelPar_given_WAGS=mod.given("WAGS");
  fadType modelPar_wa1=mod.wa1;	  bool modelPar_given_WA1=mod.given("WA1");
  fadType modelPar_wa2=mod.wa2;	  bool modelPar_given_WA2=mod.given("WA2");
  fadType modelPar_wketa=mod.wketa;	  bool modelPar_given_WKETA=mod.given("WKETA");
  fadType modelPar_wnsub=mod.wnsub;	  bool modelPar_given_WNSUB=mod.given("WNSUB");
  fadType modelPar_wnpeak=mod.wnpeak;	  bool modelPar_given_WNCH=mod.given("WNCH");
  fadType modelPar_wngate=mod.wngate;	  bool modelPar_given_WNGATE=mod.given("WNGATE");
  fadType modelPar_wgamma1=mod.wgamma1;	  bool modelPar_given_WGAMMA1=mod.given("WGAMMA1");
  fadType modelPar_wgamma2=mod.wgamma2;	  bool modelPar_given_WGAMMA2=mod.given("WGAMMA2");
  fadType modelPar_wvbx=mod.wvbx;	  bool modelPar_given_WVBX=mod.given("WVBX");
  fadType modelPar_wvbm=mod.wvbm;	  bool modelPar_given_WVBM=mod.given("WVBM");
  fadType modelPar_wxt=mod.wxt;	  bool modelPar_given_WXT=mod.given("WXT");
  fadType modelPar_wk1=mod.wk1;	  bool modelPar_given_WK1=mod.given("WK1");
  fadType modelPar_wkt1=mod.wkt1;	  bool modelPar_given_WKT1=mod.given("WKT1");
  fadType modelPar_wkt1l=mod.wkt1l;	  bool modelPar_given_WKT1L=mod.given("WKT1L");
  fadType modelPar_wkt2=mod.wkt2;	  bool modelPar_given_WKT2=mod.given("WKT2");
  fadType modelPar_wk2=mod.wk2;	  bool modelPar_given_WK2=mod.given("WK2");
  fadType modelPar_wk3=mod.wk3;	  bool modelPar_given_WK3=mod.given("WK3");
  fadType modelPar_wk3b=mod.wk3b;	  bool modelPar_given_WK3B=mod.given("WK3B");
  fadType modelPar_ww0=mod.ww0;	  bool modelPar_given_WW0=mod.given("WW0");
  fadType modelPar_wnlx=mod.wnlx;	  bool modelPar_given_WNLX=mod.given("WNLX");
  fadType modelPar_wdvt0=mod.wdvt0;	  bool modelPar_given_WDVT0=mod.given("WDVT0");
  fadType modelPar_wdvt1=mod.wdvt1;	  bool modelPar_given_WDVT1=mod.given("WDVT1");
  fadType modelPar_wdvt2=mod.wdvt2;	  bool modelPar_given_WDVT2=mod.given("WDVT2");
  fadType modelPar_wdvt0w=mod.wdvt0w;	  bool modelPar_given_WDVT0W=mod.given("WDVT0W");
  fadType modelPar_wdvt1w=mod.wdvt1w;	  bool modelPar_given_WDVT1W=mod.given("WDVT1W");
  fadType modelPar_wdvt2w=mod.wdvt2w;	  bool modelPar_given_WDVT2W=mod.given("WDVT2W");
  fadType modelPar_wdrout=mod.wdrout;	  bool modelPar_given_WDROUT=mod.given("WDROUT");
  fadType modelPar_wdsub=mod.wdsub;	  bool modelPar_given_WDSUB=mod.given("WDSUB");
  fadType modelPar_wvth0=mod.wvth0;	  bool modelPar_given_WVTH0=mod.given("WVTH0");
  fadType modelPar_wua=mod.wua;	  bool modelPar_given_WUA=mod.given("WUA");
  fadType modelPar_wua1=mod.wua1;	  bool modelPar_given_WUA1=mod.given("WUA1");
  fadType modelPar_wub=mod.wub;	  bool modelPar_given_WUB=mod.given("WUB");
  fadType modelPar_wub1=mod.wub1;	  bool modelPar_given_WUB1=mod.given("WUB1");
  fadType modelPar_wuc=mod.wuc;	  bool modelPar_given_WUC=mod.given("WUC");
  fadType modelPar_wuc1=mod.wuc1;	  bool modelPar_given_WUC1=mod.given("WUC1");
  fadType modelPar_wu0=mod.wu0;	  bool modelPar_given_WU0=mod.given("WU0");
  fadType modelPar_wute=mod.wute;	  bool modelPar_given_WUTE=mod.given("WUTE");
  fadType modelPar_wvoff=mod.wvoff;	  bool modelPar_given_WVOFF=mod.given("WVOFF");
  fadType modelPar_wrdsw=mod.wrdsw;	  bool modelPar_given_WRDSW=mod.given("WRDSW");
  fadType modelPar_wprwg=mod.wprwg;	  bool modelPar_given_WPRWG=mod.given("WPRWG");
  fadType modelPar_wprwb=mod.wprwb;	  bool modelPar_given_WPRWB=mod.given("WPRWB");
  fadType modelPar_wprt=mod.wprt;	  bool modelPar_given_WPRT=mod.given("WPRT");
  fadType modelPar_weta0=mod.weta0;	  bool modelPar_given_WETA0=mod.given("WETA0");
  fadType modelPar_wetab=mod.wetab;	  bool modelPar_given_WETAB=mod.given("WETAB");
  fadType modelPar_wpclm=mod.wpclm;	  bool modelPar_given_WPCLM=mod.given("WPCLM");
  fadType modelPar_wpdibl1=mod.wpdibl1;	  bool modelPar_given_WPDIBLC1=mod.given("WPDIBLC1");
  fadType modelPar_wpdibl2=mod.wpdibl2;	  bool modelPar_given_WPDIBLC2=mod.given("WPDIBLC2");
  fadType modelPar_wpdiblb=mod.wpdiblb;	  bool modelPar_given_WPDIBLCB=mod.given("WPDIBLCB");
  fadType modelPar_wpscbe1=mod.wpscbe1;	  bool modelPar_given_WPSCBE1=mod.given("WPSCBE1");
  fadType modelPar_wpscbe2=mod.wpscbe2;	  bool modelPar_given_WPSCBE2=mod.given("WPSCBE2");
  fadType modelPar_wpvag=mod.wpvag;	  bool modelPar_given_WPVAG=mod.given("WPVAG");
  fadType modelPar_wdelta=mod.wdelta;	  bool modelPar_given_WDELTA=mod.given("WDELTA");
  fadType modelPar_wwr=mod.wwr;	  bool modelPar_given_WWR=mod.given("WWR");
  fadType modelPar_wdwg=mod.wdwg;	  bool modelPar_given_WDWG=mod.given("WDWG");
  fadType modelPar_wdwb=mod.wdwb;	  bool modelPar_given_WDWB=mod.given("WDWB");
  fadType modelPar_wb0=mod.wb0;	  bool modelPar_given_WB0=mod.given("WB0");
  fadType modelPar_wb1=mod.wb1;	  bool modelPar_given_WB1=mod.given("WB1");
  fadType modelPar_walpha0=mod.walpha0;	  bool modelPar_given_WALPHA0=mod.given("WALPHA0");
  fadType modelPar_walpha1=mod.walpha1;	  bool modelPar_given_WALPHA1=mod.given("WALPHA1");
  fadType modelPar_wbeta0=mod.wbeta0;	  bool modelPar_given_WBETA0=mod.given("WBETA0");
  fadType modelPar_wvfb=mod.wvfb;	  bool modelPar_given_WVFB=mod.given("WVFB");
  fadType modelPar_welm=mod.welm;	  bool modelPar_given_WELM=mod.given("WELM");
  fadType modelPar_wcgsl=mod.wcgsl;	  bool modelPar_given_WCGSL=mod.given("WCGSL");
  fadType modelPar_wcgdl=mod.wcgdl;	  bool modelPar_given_WCGDL=mod.given("WCGDL");
  fadType modelPar_wckappa=mod.wckappa;	  bool modelPar_given_WCKAPPA=mod.given("WCKAPPA");
  fadType modelPar_wcf=mod.wcf;	  bool modelPar_given_WCF=mod.given("WCF");
  fadType modelPar_wclc=mod.wclc;	  bool modelPar_given_WCLC=mod.given("WCLC");
  fadType modelPar_wcle=mod.wcle;	  bool modelPar_given_WCLE=mod.given("WCLE");
  fadType modelPar_wvfbcv=mod.wvfbcv;	  bool modelPar_given_WVFBCV=mod.given("WVFBCV");
  fadType modelPar_wnoff=mod.wnoff;	  bool modelPar_given_WNOFF=mod.given("WNOFF");
  fadType modelPar_wvoffcv=mod.wvoffcv;	  bool modelPar_given_WVOFFCV=mod.given("WVOFFCV");
  fadType modelPar_wacde=mod.wacde;	  bool modelPar_given_WACDE=mod.given("WACDE");
  fadType modelPar_wmoin=mod.wmoin;	  bool modelPar_given_WMOIN=mod.given("WMOIN");
  fadType modelPar_pcdsc=mod.pcdsc;	  bool modelPar_given_PCDSC=mod.given("PCDSC");
  fadType modelPar_pcdscb=mod.pcdscb;	  bool modelPar_given_PCDSCB=mod.given("PCDSCB");
  fadType modelPar_pcdscd=mod.pcdscd;	  bool modelPar_given_PCDSCD=mod.given("PCDSCD");
  fadType modelPar_pcit=mod.pcit;	  bool modelPar_given_PCIT=mod.given("PCIT");
  fadType modelPar_pnfactor=mod.pnfactor;	  bool modelPar_given_PNFACTOR=mod.given("PNFACTOR");
  fadType modelPar_pxj=mod.pxj;	  bool modelPar_given_PXJ=mod.given("PXJ");
  fadType modelPar_pvsat=mod.pvsat;	  bool modelPar_given_PVSAT=mod.given("PVSAT");
  fadType modelPar_pat=mod.pat;	  bool modelPar_given_PAT=mod.given("PAT");
  fadType modelPar_pa0=mod.pa0;	  bool modelPar_given_PA0=mod.given("PA0");
  fadType modelPar_pags=mod.pags;	  bool modelPar_given_PAGS=mod.given("PAGS");
  fadType modelPar_pa1=mod.pa1;	  bool modelPar_given_PA1=mod.given("PA1");
  fadType modelPar_pa2=mod.pa2;	  bool modelPar_given_PA2=mod.given("PA2");
  fadType modelPar_pketa=mod.pketa;	  bool modelPar_given_PKETA=mod.given("PKETA");
  fadType modelPar_pnsub=mod.pnsub;	  bool modelPar_given_PNSUB=mod.given("PNSUB");
  fadType modelPar_pnpeak=mod.pnpeak;	  bool modelPar_given_PNCH=mod.given("PNCH");
  fadType modelPar_pngate=mod.pngate;	  bool modelPar_given_PNGATE=mod.given("PNGATE");
  fadType modelPar_pgamma1=mod.pgamma1;	  bool modelPar_given_PGAMMA1=mod.given("PGAMMA1");
  fadType modelPar_pgamma2=mod.pgamma2;	  bool modelPar_given_PGAMMA2=mod.given("PGAMMA2");
  fadType modelPar_pvbx=mod.pvbx;	  bool modelPar_given_PVBX=mod.given("PVBX");
  fadType modelPar_pvbm=mod.pvbm;	  bool modelPar_given_PVBM=mod.given("PVBM");
  fadType modelPar_pxt=mod.pxt;	  bool modelPar_given_PXT=mod.given("PXT");
  fadType modelPar_pk1=mod.pk1;	  bool modelPar_given_PK1=mod.given("PK1");
  fadType modelPar_pkt1=mod.pkt1;	  bool modelPar_given_PKT1=mod.given("PKT1");
  fadType modelPar_pkt1l=mod.pkt1l;	  bool modelPar_given_PKT1L=mod.given("PKT1L");
  fadType modelPar_pkt2=mod.pkt2;	  bool modelPar_given_PKT2=mod.given("PKT2");
  fadType modelPar_pk2=mod.pk2;	  bool modelPar_given_PK2=mod.given("PK2");
  fadType modelPar_pk3=mod.pk3;	  bool modelPar_given_PK3=mod.given("PK3");
  fadType modelPar_pk3b=mod.pk3b;	  bool modelPar_given_PK3B=mod.given("PK3B");
  fadType modelPar_pw0=mod.pw0;	  bool modelPar_given_PW0=mod.given("PW0");
  fadType modelPar_pnlx=mod.pnlx;	  bool modelPar_given_PNLX=mod.given("PNLX");
  fadType modelPar_pdvt0=mod.pdvt0;	  bool modelPar_given_PDVT0=mod.given("PDVT0");
  fadType modelPar_pdvt1=mod.pdvt1;	  bool modelPar_given_PDVT1=mod.given("PDVT1");
  fadType modelPar_pdvt2=mod.pdvt2;	  bool modelPar_given_PDVT2=mod.given("PDVT2");
  fadType modelPar_pdvt0w=mod.pdvt0w;	  bool modelPar_given_PDVT0W=mod.given("PDVT0W");
  fadType modelPar_pdvt1w=mod.pdvt1w;	  bool modelPar_given_PDVT1W=mod.given("PDVT1W");
  fadType modelPar_pdvt2w=mod.pdvt2w;	  bool modelPar_given_PDVT2W=mod.given("PDVT2W");
  fadType modelPar_pdrout=mod.pdrout;	  bool modelPar_given_PDROUT=mod.given("PDROUT");
  fadType modelPar_pdsub=mod.pdsub;	  bool modelPar_given_PDSUB=mod.given("PDSUB");
  fadType modelPar_pvth0=mod.pvth0;	  bool modelPar_given_PVTH0=mod.given("PVTH0");
  fadType modelPar_pua=mod.pua;	  bool modelPar_given_PUA=mod.given("PUA");
  fadType modelPar_pua1=mod.pua1;	  bool modelPar_given_PUA1=mod.given("PUA1");
  fadType modelPar_pub=mod.pub;	  bool modelPar_given_PUB=mod.given("PUB");
  fadType modelPar_pub1=mod.pub1;	  bool modelPar_given_PUB1=mod.given("PUB1");
  fadType modelPar_puc=mod.puc;	  bool modelPar_given_PUC=mod.given("PUC");
  fadType modelPar_puc1=mod.puc1;	  bool modelPar_given_PUC1=mod.given("PUC1");
  fadType modelPar_pu0=mod.pu0;	  bool modelPar_given_PU0=mod.given("PU0");
  fadType modelPar_pute=mod.pute;	  bool modelPar_given_PUTE=mod.given("PUTE");
  fadType modelPar_pvoff=mod.pvoff;	  bool modelPar_given_PVOFF=mod.given("PVOFF");
  fadType modelPar_prdsw=mod.prdsw;	  bool modelPar_given_PRDSW=mod.given("PRDSW");
  fadType modelPar_pprwg=mod.pprwg;	  bool modelPar_given_PPRWG=mod.given("PPRWG");
  fadType modelPar_pprwb=mod.pprwb;	  bool modelPar_given_PPRWB=mod.given("PPRWB");
  fadType modelPar_pprt=mod.pprt;	  bool modelPar_given_PPRT=mod.given("PPRT");
  fadType modelPar_peta0=mod.peta0;	  bool modelPar_given_PETA0=mod.given("PETA0");
  fadType modelPar_petab=mod.petab;	  bool modelPar_given_PETAB=mod.given("PETAB");
  fadType modelPar_ppclm=mod.ppclm;	  bool modelPar_given_PPCLM=mod.given("PPCLM");
  fadType modelPar_ppdibl1=mod.ppdibl1;	  bool modelPar_given_PPDIBLC1=mod.given("PPDIBLC1");
  fadType modelPar_ppdibl2=mod.ppdibl2;	  bool modelPar_given_PPDIBLC2=mod.given("PPDIBLC2");
  fadType modelPar_ppdiblb=mod.ppdiblb;	  bool modelPar_given_PPDIBLCB=mod.given("PPDIBLCB");
  fadType modelPar_ppscbe1=mod.ppscbe1;	  bool modelPar_given_PPSCBE1=mod.given("PPSCBE1");
  fadType modelPar_ppscbe2=mod.ppscbe2;	  bool modelPar_given_PPSCBE2=mod.given("PPSCBE2");
  fadType modelPar_ppvag=mod.ppvag;	  bool modelPar_given_PPVAG=mod.given("PPVAG");
  fadType modelPar_pdelta=mod.pdelta;	  bool modelPar_given_PDELTA=mod.given("PDELTA");
  fadType modelPar_pwr=mod.pwr;	  bool modelPar_given_PWR=mod.given("PWR");
  fadType modelPar_pdwg=mod.pdwg;	  bool modelPar_given_PDWG=mod.given("PDWG");
  fadType modelPar_pdwb=mod.pdwb;	  bool modelPar_given_PDWB=mod.given("PDWB");
  fadType modelPar_pb0=mod.pb0;	  bool modelPar_given_PB0=mod.given("PB0");
  fadType modelPar_pb1=mod.pb1;	  bool modelPar_given_PB1=mod.given("PB1");
  fadType modelPar_palpha0=mod.palpha0;	  bool modelPar_given_PALPHA0=mod.given("PALPHA0");
  fadType modelPar_palpha1=mod.palpha1;	  bool modelPar_given_PALPHA1=mod.given("PALPHA1");
  fadType modelPar_pbeta0=mod.pbeta0;	  bool modelPar_given_PBETA0=mod.given("PBETA0");
  fadType modelPar_pvfb=mod.pvfb;	  bool modelPar_given_PVFB=mod.given("PVFB");
  fadType modelPar_pelm=mod.pelm;	  bool modelPar_given_PELM=mod.given("PELM");
  fadType modelPar_pcgsl=mod.pcgsl;	  bool modelPar_given_PCGSL=mod.given("PCGSL");
  fadType modelPar_pcgdl=mod.pcgdl;	  bool modelPar_given_PCGDL=mod.given("PCGDL");
  fadType modelPar_pckappa=mod.pckappa;	  bool modelPar_given_PCKAPPA=mod.given("PCKAPPA");
  fadType modelPar_pcf=mod.pcf;	  bool modelPar_given_PCF=mod.given("PCF");
  fadType modelPar_pclc=mod.pclc;	  bool modelPar_given_PCLC=mod.given("PCLC");
  fadType modelPar_pcle=mod.pcle;	  bool modelPar_given_PCLE=mod.given("PCLE");
  fadType modelPar_pvfbcv=mod.pvfbcv;	  bool modelPar_given_PVFBCV=mod.given("PVFBCV");
  fadType modelPar_pnoff=mod.pnoff;	  bool modelPar_given_PNOFF=mod.given("PNOFF");
  fadType modelPar_pvoffcv=mod.pvoffcv;	  bool modelPar_given_PVOFFCV=mod.given("PVOFFCV");
  fadType modelPar_pacde=mod.pacde;	  bool modelPar_given_PACDE=mod.given("PACDE");
  fadType modelPar_pmoin=mod.pmoin;	  bool modelPar_given_PMOIN=mod.given("PMOIN");
  fadType modelPar_tnom=mod.tnom;	  bool modelPar_given_TNOM=mod.given("TNOM");
  fadType modelPar_cgso=mod.cgso;	  bool modelPar_given_CGSO=mod.given("CGSO");
  fadType modelPar_cgdo=mod.cgdo;	  bool modelPar_given_CGDO=mod.given("CGDO");
  fadType modelPar_cgbo=mod.cgbo;	  bool modelPar_given_CGBO=mod.given("CGBO");
  fadType modelPar_xpart=mod.xpart;	  bool modelPar_given_XPART=mod.given("XPART");
  fadType modelPar_sheetResistance=mod.sheetResistance;	  bool modelPar_given_RSH=mod.given("RSH");
  fadType modelPar_jctSatCurDensity=mod.jctSatCurDensity;	  bool modelPar_given_JS=mod.given("JS");
  fadType modelPar_jctSidewallSatCurDensity=mod.jctSidewallSatCurDensity;	  bool modelPar_given_JSW=mod.given("JSW");
  fadType modelPar_bulkJctPotential=mod.bulkJctPotential;	  bool modelPar_given_PB=mod.given("PB");
  fadType modelPar_bulkJctBotGradingCoeff=mod.bulkJctBotGradingCoeff;	  bool modelPar_given_MJ=mod.given("MJ");
  fadType modelPar_sidewallJctPotential=mod.sidewallJctPotential;	  bool modelPar_given_PBSW=mod.given("PBSW");
  fadType modelPar_GatesidewallJctPotential=mod.GatesidewallJctPotential;	  bool modelPar_given_PBSWG=mod.given("PBSWG");
  fadType modelPar_bulkJctSideGradingCoeff=mod.bulkJctSideGradingCoeff;	  bool modelPar_given_MJSW=mod.given("MJSW");
  fadType modelPar_unitAreaJctCap=mod.unitAreaJctCap;	  bool modelPar_given_CJ=mod.given("CJ");
  fadType modelPar_unitLengthSidewallJctCap=mod.unitLengthSidewallJctCap;	  bool modelPar_given_CJSW=mod.given("CJSW");
  fadType modelPar_bulkJctGateSideGradingCoeff=mod.bulkJctGateSideGradingCoeff;	  bool modelPar_given_MJSWG=mod.given("MJSWG");
  fadType modelPar_unitLengthGateSidewallJctCap=mod.unitLengthGateSidewallJctCap;	  bool modelPar_given_CJSWG=mod.given("CJSWG");
  fadType modelPar_jctEmissionCoeff=mod.jctEmissionCoeff;	  bool modelPar_given_NJ=mod.given("NJ");
  fadType modelPar_jctTempExponent=mod.jctTempExponent;	  bool modelPar_given_XTI=mod.given("XTI");
  fadType modelPar_oxideTrapDensityA=mod.oxideTrapDensityA;	  bool modelPar_given_NOIA=mod.given("NOIA");
  fadType modelPar_oxideTrapDensityB=mod.oxideTrapDensityB;	  bool modelPar_given_NOIB=mod.given("NOIB");
  fadType modelPar_oxideTrapDensityC=mod.oxideTrapDensityC;	  bool modelPar_given_NOIC=mod.given("NOIC");
  fadType modelPar_em=mod.em;	  bool modelPar_given_EM=mod.given("EM");
  fadType modelPar_ef=mod.ef;	  bool modelPar_given_EF=mod.given("EF");
  fadType modelPar_af=mod.af;	  bool modelPar_given_AF=mod.given("AF");
  fadType modelPar_kf=mod.kf;	  bool modelPar_given_KF=mod.given("KF");
  fadType modelPar_lintnoi=mod.lintnoi;	  bool modelPar_given_LINTNOI=mod.given("LINTNOI");
  fadType modelPar_Lint=mod.Lint;	  bool modelPar_given_LINT=mod.given("LINT");
  fadType modelPar_Ll=mod.Ll;	  bool modelPar_given_LL=mod.given("LL");
  fadType modelPar_Llc=mod.Llc;	  bool modelPar_given_LLC=mod.given("LLC");
  fadType modelPar_Lln=mod.Lln;	  bool modelPar_given_LLN=mod.given("LLN");
  fadType modelPar_Lw=mod.Lw;	  bool modelPar_given_LW=mod.given("LW");
  fadType modelPar_Lwc=mod.Lwc;	  bool modelPar_given_LWC=mod.given("LWC");
  fadType modelPar_Lwn=mod.Lwn;	  bool modelPar_given_LWN=mod.given("LWN");
  fadType modelPar_Lwl=mod.Lwl;	  bool modelPar_given_LWL=mod.given("LWL");
  fadType modelPar_Lwlc=mod.Lwlc;	  bool modelPar_given_LWLC=mod.given("LWLC");
  fadType modelPar_Wint=mod.Wint;	  bool modelPar_given_WINT=mod.given("WINT");
  fadType modelPar_Wl=mod.Wl;	  bool modelPar_given_WL=mod.given("WL");
  fadType modelPar_Wlc=mod.Wlc;	  bool modelPar_given_WLC=mod.given("WLC");
  fadType modelPar_Wln=mod.Wln;	  bool modelPar_given_WLN=mod.given("WLN");
  fadType modelPar_Ww=mod.Ww;	  bool modelPar_given_WW=mod.given("WW");
  fadType modelPar_Wwc=mod.Wwc;	  bool modelPar_given_WWC=mod.given("WWC");
  fadType modelPar_Wwn=mod.Wwn;	  bool modelPar_given_WWN=mod.given("WWN");
  fadType modelPar_Wwl=mod.Wwl;	  bool modelPar_given_WWL=mod.given("WWL");
  fadType modelPar_Wwlc=mod.Wwlc;	  bool modelPar_given_WWLC=mod.given("WWLC");
  fadType modelPar_model_l=mod.model_l;	  bool modelPar_given_L=mod.given("L");
  fadType modelPar_model_w=mod.model_w;	  bool modelPar_given_W=mod.given("W");
  fadType modelPar_Lmax=mod.Lmax;	  bool modelPar_given_LMAX=mod.given("LMAX");
  fadType modelPar_Lmin=mod.Lmin;	  bool modelPar_given_LMIN=mod.given("LMIN");
  fadType modelPar_Wmax=mod.Wmax;	  bool modelPar_given_WMAX=mod.given("WMAX");
  fadType modelPar_Wmin=mod.Wmin;	  bool modelPar_given_WMIN=mod.given("WMIN");

  // model params NOT sensitivity supported:
  int  dtype = mod.dtype;

  int mobMod = mod.mobMod; // int
  int binUnit = mod.binUnit; // int
  int capMod = mod.capMod; // int
  int paramChk = mod.paramChk; // int
  int noiMod = mod.noiMod; // int
  std::string version = mod.version; // string

  unordered_map <std::string,fadType*,HashNoCase,EqualNoCase> inParamMap;

  // instance params
  
  fadType instancePar_temp=inst.temp;	  
  bool instancePar_given_temp=inst.given("TEMP");
  inParamMap["TEMP"] = &instancePar_temp;

  fadType instancePar_dtemp=inst.dtemp;	  
  bool instancePar_given_dtemp=inst.given("DTEMP");
  inParamMap["DTEMP"] = &instancePar_dtemp;
  
  fadType instancePar_l=inst.l;	  
  bool instancePar_given_l=inst.given("L");
  inParamMap["L"] = &instancePar_l;
  
  fadType instancePar_w=inst.w;	  
  bool instancePar_given_w=inst.given("W");
  inParamMap["W"] = &instancePar_w;
  
  fadType instancePar_drainArea=inst.drainArea;	  
  bool instancePar_given_drainArea=inst.given("AD");
  inParamMap["AD"] = &instancePar_drainArea;
  
  fadType instancePar_sourceArea=inst.sourceArea;	  
  bool instancePar_given_sourceArea=inst.given("AS");
  inParamMap["AS"] = &instancePar_sourceArea;
  
  fadType instancePar_drainSquares=inst.drainSquares;	  
  bool instancePar_given_drainSquares=inst.given("NRD");
  inParamMap["NRD"] = &instancePar_drainSquares;
  
  fadType instancePar_sourceSquares=inst.sourceSquares;	  
  bool instancePar_given_sourceSquares=inst.given("NRS");
  inParamMap["NRS"] = &instancePar_sourceSquares;
  
  fadType instancePar_drainPerimeter=inst.drainPerimeter;	  
  bool instancePar_given_drainPerimeter=inst.given("PD");
  inParamMap["PD"] = &instancePar_drainPerimeter;
  
  fadType instancePar_sourcePerimeter=inst.sourcePerimeter;	  
  bool instancePar_given_sourcePerimeter=inst.given("PS");
  inParamMap["PS"] = &instancePar_sourcePerimeter;

  bool instancePar_given_nqsMod = inst.given("NQSMOD");

  inParamMap[name]->diff(0,1);

  // instance params NOT sensitivity supported:  
  double numberParallel = inst.numberParallel ;
  double icVDS = inst.icVDS ;
  double icVGS = inst.icVGS ;
  double icVBS = inst.icVBS ;
  int nqsMod = inst.nqsMod ;
  bool OFF = inst.OFF ;

  // other variables
  fadType drainConductance, sourceConductance;

  fadType modelPar_cox = mod.cox;
  fadType modelPar_vcrit = mod.vcrit;
  fadType modelPar_factor1 = mod.factor1;
  fadType modelPar_Vtm0 = mod.Vtm0;
  fadType modelPar_Eg0 = mod.Eg0;
  fadType modelPar_ni  = mod.ni;

  // model process params:
  processModelParams (
      // givens:
      modelPar_given_TOXM, modelPar_given_DSUB, modelPar_given_LLC, modelPar_given_LWC,
      modelPar_given_LWLC, modelPar_given_WLC, modelPar_given_WWL, modelPar_given_WWLC,
      modelPar_given_DWC, modelPar_given_DLC, modelPar_given_CF, modelPar_given_CGDO,
      modelPar_given_CGSO, modelPar_given_CGBO, modelPar_given_CJSWG, modelPar_given_PBSWG,
      modelPar_given_MJSWG, 
      // inputs:
      modelPar_tox, modelPar_drout, modelPar_Ll, modelPar_Lw,
      modelPar_Lwl, modelPar_Wl, modelPar_Wwl, modelPar_Wint,
      modelPar_Lint, modelPar_tnom, 
      modelPar_cgdl, modelPar_cgsl,
      modelPar_unitLengthSidewallJctCap, modelPar_bulkJctSideGradingCoeff, modelPar_xj, 
      //outputs:
      modelPar_toxm, modelPar_cox, modelPar_dsub, modelPar_Llc, modelPar_Lwc,
      modelPar_Lwlc, modelPar_Wlc, modelPar_Wwlc, modelPar_dwc, modelPar_dlc,
      modelPar_cf, modelPar_cgdo, modelPar_cgso, modelPar_cgbo,
      modelPar_unitLengthGateSidewallJctCap, modelPar_GatesidewallJctPotential, modelPar_bulkJctGateSideGradingCoeff,
      modelPar_bulkJctPotential, modelPar_sidewallJctPotential,
      modelPar_vcrit, modelPar_factor1, modelPar_Vtm0, modelPar_Eg0, modelPar_ni 
      );

  // instance variables:
  processParams (
    inst.getDeviceOptions(),

    instancePar_given_temp,
    instancePar_given_dtemp,
    instancePar_given_l,
    instancePar_given_w,
  instancePar_given_drainArea,
  instancePar_given_sourceArea,
  instancePar_given_nqsMod,

   // inputs
    modelPar_model_l,
    modelPar_model_w,

    modelPar_sheetResistance,

    instancePar_drainSquares,
    instancePar_sourceSquares,

   // outputs
    instancePar_temp,
    instancePar_dtemp,
    instancePar_l,
    instancePar_w,
    instancePar_drainArea,
    instancePar_sourceArea,
      drainConductance,
      sourceConductance,
     nqsMod
    );

  fadType instance_PhiBSWTemp = inst.PhiBSWTemp;
  fadType instance_PhiBSWGTemp = inst.PhiBSWGTemp;

  fadType instance_vtm = inst.vtm;
  fadType model_vtm = mod.vtm;

  fadType instance_jctTempSatCurDensity = inst.jctTempSatCurDensity;
  fadType instance_jctSidewallTempSatCurDensity = inst.jctSidewallTempSatCurDensity;
  fadType instance_unitAreaJctCapTemp = inst.unitAreaJctCapTemp;
  fadType instance_unitLengthSidewallJctCapTemp = inst.unitLengthSidewallJctCapTemp;
  fadType instance_unitLengthGateSidewallJctCapTemp = inst.unitLengthGateSidewallJctCapTemp;

  fadType instance_PhiBTemp = inst.PhiBTemp;

  fadType instance_cgso = inst.cgso;
  fadType instance_cgdo = inst.cgdo;
  fadType instance_vjsm = inst.vjsm;
  fadType instance_IsEvjsm = inst.IsEvjsm;
  fadType instance_vjdm = inst.vjdm;
  fadType instance_IsEvjdm = inst.IsEvjdm;

  bool instance_updateTemperatureCalled_ = inst.updateTemperatureCalled_;  // seems useless  ... 

  // ----------------------------------------------------------------------------
  SizeDependParam<fadType> sizeDepParams;

  updateTemperature (
    instancePar_temp,
    instancePar_dtemp,
  modelPar_Eg0,
  modelPar_ni,
  modelPar_Vtm0,
    modelPar_tnom,
    modelPar_jctTempExponent,
    modelPar_jctEmissionCoeff,
    modelPar_jctSatCurDensity,
    modelPar_jctSidewallSatCurDensity,
    modelPar_tcj,
    modelPar_unitAreaJctCap,
    modelPar_tcjsw,
    modelPar_unitLengthSidewallJctCap, 
    modelPar_tcjswg,
    modelPar_unitLengthGateSidewallJctCap,
    modelPar_bulkJctPotential,
    modelPar_tpb,
  instance_PhiBSWTemp, // non-CONST !!
    modelPar_sidewallJctPotential,
    modelPar_tpbsw,
  instance_PhiBSWGTemp, // non-CONST !!
    modelPar_GatesidewallJctPotential,
    modelPar_tpbswg,
    instancePar_l,
    instancePar_w,
    modelPar_Lln,
    modelPar_Lwn,
    modelPar_Ll,
    modelPar_Lw,
    modelPar_Lwl,
    modelPar_Lint,
    modelPar_Llc,
    modelPar_Lwc,
    modelPar_Lwlc,
    modelPar_dlc,
    modelPar_Wln,
    modelPar_Wwn,
    modelPar_Wl,
    modelPar_Ww,
    modelPar_Wwl,
    modelPar_Wint,
    modelPar_Wlc,
    modelPar_Wwc,
    modelPar_Wwlc,
    modelPar_dwc,
    modelPar_cgdo,
    modelPar_cgso,
    modelPar_cgbo,
    modelPar_cox,
    modelPar_given_NCH, //modelPar_npeakGiven,
    modelPar_given_GAMMA1, //modelPar_gamma1Given,
    modelPar_ijth,
    instancePar_sourceArea, 
    instancePar_sourcePerimeter,
    instancePar_drainArea,
    instancePar_drainPerimeter,
    modelPar_cdsc,
    modelPar_lcdsc,
    modelPar_wcdsc,
    modelPar_pcdsc,
    modelPar_cdscb,
    modelPar_lcdscb,
    modelPar_wcdscb,
    modelPar_pcdscb,
    modelPar_cdscd,
    modelPar_lcdscd,
    modelPar_wcdscd,
    modelPar_pcdscd,
    modelPar_cit,
    modelPar_lcit,
    modelPar_wcit,
    modelPar_pcit,
    modelPar_nfactor,
    modelPar_lnfactor,
    modelPar_wnfactor,
    modelPar_pnfactor,
    modelPar_xj,
    modelPar_lxj,
    modelPar_wxj,
    modelPar_pxj,
    modelPar_vsat,
    modelPar_lvsat,
    modelPar_wvsat,
    modelPar_pvsat,
    modelPar_at,
    modelPar_lat,
    modelPar_wat,
    modelPar_pat,
    modelPar_a0,
    modelPar_la0,
    modelPar_wa0,
    modelPar_pa0,
    modelPar_ags,
    modelPar_lags,
    modelPar_wags,
    modelPar_pags,
    modelPar_a1,
    modelPar_la1,
    modelPar_wa1,
    modelPar_pa1,
    modelPar_a2,
    modelPar_la2,
    modelPar_wa2,
    modelPar_pa2,
    modelPar_keta,
    modelPar_lketa,
    modelPar_wketa,
    modelPar_pketa,
    modelPar_nsub,
    modelPar_lnsub,
    modelPar_wnsub,
    modelPar_pnsub,
    modelPar_npeak,
    modelPar_lnpeak,
    modelPar_wnpeak,
    modelPar_pnpeak,
    modelPar_ngate,
    modelPar_lngate,
    modelPar_wngate,
    modelPar_pngate,
    modelPar_gamma1,
    modelPar_lgamma1,
    modelPar_wgamma1,
    modelPar_pgamma1,
    modelPar_gamma2,
    modelPar_lgamma2,
    modelPar_wgamma2,
    modelPar_pgamma2,
    modelPar_vbx,
    modelPar_lvbx,
    modelPar_wvbx,
    modelPar_pvbx,
    modelPar_vbm,
    modelPar_lvbm,
    modelPar_wvbm,
    modelPar_pvbm,
    modelPar_xt,
    modelPar_lxt,
    modelPar_wxt,
    modelPar_pxt,
    modelPar_vfb,
    modelPar_lvfb,
    modelPar_wvfb,
    modelPar_pvfb,
    modelPar_k1,
    modelPar_lk1,
    modelPar_wk1,
    modelPar_pk1,
    modelPar_kt1,
    modelPar_lkt1,
    modelPar_wkt1,
    modelPar_pkt1,
    modelPar_kt1l,
    modelPar_lkt1l,
    modelPar_wkt1l,
    modelPar_pkt1l,
    modelPar_k2,
    modelPar_lk2,
    modelPar_wk2,
    modelPar_pk2,
    modelPar_kt2,
    modelPar_lkt2,
    modelPar_wkt2,
    modelPar_pkt2,
    modelPar_k3,
    modelPar_lk3,
    modelPar_wk3,
    modelPar_pk3,
    modelPar_k3b,
    modelPar_lk3b,
    modelPar_wk3b,
    modelPar_pk3b,
    modelPar_w0,
    modelPar_lw0,
    modelPar_ww0,
    modelPar_pw0,
    modelPar_nlx,
    modelPar_lnlx,
    modelPar_wnlx,
    modelPar_pnlx,
    modelPar_dvt0,
    modelPar_ldvt0,
    modelPar_wdvt0,
    modelPar_pdvt0,
    modelPar_dvt1,
    modelPar_ldvt1,
    modelPar_wdvt1,
    modelPar_pdvt1,
    modelPar_dvt2,
    modelPar_ldvt2,
    modelPar_wdvt2,
    modelPar_pdvt2,
    modelPar_dvt0w,
    modelPar_ldvt0w,
    modelPar_wdvt0w,
    modelPar_pdvt0w,
    modelPar_dvt1w,
    modelPar_ldvt1w,
    modelPar_wdvt1w,
    modelPar_pdvt1w,
    modelPar_dvt2w,
    modelPar_ldvt2w,
    modelPar_wdvt2w,
    modelPar_pdvt2w,
    modelPar_drout,
    modelPar_ldrout,
    modelPar_wdrout,
    modelPar_pdrout,
    modelPar_dsub,
    modelPar_ldsub,
    modelPar_wdsub,
    modelPar_pdsub,
    modelPar_vth0,
    modelPar_lvth0,
    modelPar_wvth0,
    modelPar_pvth0,
    modelPar_ua,
    modelPar_lua,
    modelPar_wua,
    modelPar_pua,
    modelPar_ua1,
    modelPar_lua1,
    modelPar_wua1,
    modelPar_pua1,
    modelPar_ub,
    modelPar_lub,
    modelPar_wub,
    modelPar_pub,
    modelPar_ub1,
    modelPar_lub1,
    modelPar_wub1,
    modelPar_pub1,
    modelPar_uc,
    modelPar_luc,
    modelPar_wuc,
    modelPar_puc,
    modelPar_uc1,
    modelPar_luc1,
    modelPar_wuc1,
    modelPar_puc1,
    modelPar_u0,
    modelPar_lu0,
    modelPar_wu0,
    modelPar_pu0,
    modelPar_ute,
    modelPar_lute,
    modelPar_wute,
    modelPar_pute,
    modelPar_voff,
    modelPar_lvoff,
    modelPar_wvoff,
    modelPar_pvoff,
    modelPar_delta,
    modelPar_ldelta,
    modelPar_wdelta,
    modelPar_pdelta,
    modelPar_rdsw,
    modelPar_lrdsw,
    modelPar_wrdsw,
    modelPar_prdsw,
    modelPar_prwg,
    modelPar_lprwg,
    modelPar_wprwg,
    modelPar_pprwg,
    modelPar_prwb,
    modelPar_lprwb,
    modelPar_wprwb,
    modelPar_pprwb,
    modelPar_prt,
    modelPar_lprt,
    modelPar_wprt,
    modelPar_pprt,
    modelPar_eta0,
    modelPar_leta0,
    modelPar_weta0,
    modelPar_peta0,
    modelPar_etab,
    modelPar_letab,
    modelPar_wetab,
    modelPar_petab,
    modelPar_pclm,
    modelPar_lpclm,
    modelPar_wpclm,
    modelPar_ppclm,
    modelPar_pdibl1,
    modelPar_lpdibl1,
    modelPar_wpdibl1,
    modelPar_ppdibl1,
    modelPar_pdibl2,
    modelPar_lpdibl2,
    modelPar_wpdibl2,
    modelPar_ppdibl2,
    modelPar_pdiblb,
    modelPar_lpdiblb,
    modelPar_wpdiblb,
    modelPar_ppdiblb,
    modelPar_pscbe1,
    modelPar_lpscbe1,
    modelPar_wpscbe1,
    modelPar_ppscbe1,
    modelPar_pscbe2,
    modelPar_lpscbe2,
    modelPar_wpscbe2,
    modelPar_ppscbe2,
    modelPar_pvag,
    modelPar_lpvag,
    modelPar_wpvag,
    modelPar_ppvag,
    modelPar_wr,
    modelPar_lwr,
    modelPar_wwr,
    modelPar_pwr,
    modelPar_dwg,
    modelPar_ldwg,
    modelPar_wdwg,
    modelPar_pdwg,
    modelPar_dwb,
    modelPar_ldwb,
    modelPar_wdwb,
    modelPar_pdwb,
    modelPar_b0,
    modelPar_lb0,
    modelPar_wb0,
    modelPar_pb0,
    modelPar_b1,
    modelPar_lb1,
    modelPar_wb1,
    modelPar_pb1,
    modelPar_alpha0,
    modelPar_lalpha0,
    modelPar_walpha0,
    modelPar_palpha0,
    modelPar_alpha1,
    modelPar_lalpha1,
    modelPar_walpha1,
    modelPar_palpha1,
    modelPar_beta0,
    modelPar_lbeta0,
    modelPar_wbeta0,
    modelPar_pbeta0,
    modelPar_elm,
    modelPar_lelm,
    modelPar_welm,
    modelPar_pelm,
    modelPar_cgsl,
    modelPar_lcgsl,
    modelPar_wcgsl,
    modelPar_pcgsl,
    modelPar_cgdl,
    modelPar_lcgdl,
    modelPar_wcgdl,
    modelPar_pcgdl,
    modelPar_ckappa,
    modelPar_lckappa,
    modelPar_wckappa,
    modelPar_pckappa,
    modelPar_cf,
    modelPar_lcf,
    modelPar_wcf,
    modelPar_pcf,
    modelPar_clc,
    modelPar_lclc,
    modelPar_wclc,
    modelPar_pclc,
    modelPar_cle,
    modelPar_lcle,
    modelPar_wcle,
    modelPar_pcle,
    modelPar_vfbcv,
    modelPar_lvfbcv,
    modelPar_wvfbcv,
    modelPar_pvfbcv,
    modelPar_acde,
    modelPar_lacde,
    modelPar_wacde,
    modelPar_pacde,
    modelPar_moin,
    modelPar_lmoin,
    modelPar_wmoin,
    modelPar_pmoin,
    modelPar_noff,
    modelPar_lnoff,
    modelPar_wnoff,
    modelPar_pnoff,
    modelPar_voffcv,
    modelPar_lvoffcv,
    modelPar_wvoffcv,
    modelPar_pvoffcv,
    binUnit, // this is an int //modelPar_binUnit,
    modelPar_tox,
    modelPar_given_K1, //modelPar_k1Given,
    modelPar_given_K2, //modelPar_k2Given,
    modelPar_given_NSUB, //modelPar_nsubGiven,
    modelPar_given_XT, //modelPar_xtGiven,
    modelPar_given_VBX, //modelPar_vbxGiven,
    modelPar_given_GAMMA2, //modelPar_gamma2Given,
    modelPar_given_VFB, //modelPar_vfbGiven,
    modelPar_given_VTH0, //modelPar_vth0Given,
    dtype, // this is an int //modelPar_dtype,
    modelPar_toxm,
    modelPar_factor1,
    // outputs
    instance_vtm, 
    instancePar_temp, // this is goofy ...  //temp,
    instance_jctTempSatCurDensity,
    instance_jctSidewallTempSatCurDensity,
    instance_unitAreaJctCapTemp,
    instance_unitLengthSidewallJctCapTemp,
    instance_unitLengthGateSidewallJctCapTemp,
    instance_PhiBTemp,
    instance_cgso,
    instance_cgdo,
    instance_vjsm,
    instance_IsEvjsm,
    instance_vjdm,
    instance_IsEvjdm,
    instance_updateTemperatureCalled_,
    sizeDepParams
    );

  // updateIntermediateVars  stuff
  ///////////////////////////////////////////
  // obtain voltages:
  double * solVec = inst.extData.nextSolVectorRawPtr;

  fadType Vd = (solVec)[inst.li_Drain];
  fadType Vg = (solVec)[inst.li_Gate];
  fadType Vs = (solVec)[inst.li_Source];
  fadType Vb = (solVec)[inst.li_Bulk];
  fadType Vsp = (solVec)[inst.li_SourcePrime];
  fadType Vdp = (solVec)[inst.li_DrainPrime];
  fadType Qtotal = 0.0;
  if( inst.nqsMod ) { Qtotal = (solVec)[inst.li_Charge]; }

  fadType Vddp  = Vd   - Vdp;
  fadType Vssp  = Vs   - Vsp;
  fadType Vbsp  = Vb   - Vsp;
  fadType Vbdp  = Vb   - Vdp;
  fadType Vgsp  = Vg   - Vsp;
  fadType Vgdp  = Vg   - Vdp;
  fadType Vgb   = Vg   - Vb;
  fadType Vdpsp = Vdp  - Vsp;

  fadType vbs  = mod.dtype * Vbsp;
  fadType vgs  = mod.dtype * Vgsp;
  fadType vds  = mod.dtype * Vdpsp;
  fadType qdef = mod.dtype * Qtotal;
  fadType vbd = vbs - vds;
  fadType vgd = vgs - vds;

  fadType instance_T1global = inst.T1global;
  fadType instance_thetavth = inst.thetavth;

  fadType instance_Vth = inst.Vth;
  fadType instance_von = inst.von;
  fadType instance_dVgs_eff_dVg = inst.dVgs_eff_dVg;
  fadType instance_Vgsteff = inst.Vgsteff;
  fadType instance_Abulk = inst.Abulk;
  fadType instance_ueff = inst.ueff;

  fadType instance_Vdsat = inst.Vdsat;
  fadType instance_vdsat = inst.vdsat;
  fadType instance_Vdseff = inst.Vdseff;
  fadType instance_Gm = inst.Gm;
  fadType instance_cdrain = inst.cdrain;
  fadType instance_gds = inst.gds;
  fadType instance_gm = inst.gm;
  fadType instance_gmbs = inst.gmbs;
  fadType instance_gbbs = inst.gbbs;
  fadType instance_gbgs = inst.gbgs;
  fadType instance_gbds = inst.gbds;
  fadType instance_csub = inst.csub;
  fadType instance_qgate = inst.qgate;
  fadType instance_qdrn = inst.qdrn;
  fadType instance_qsrc = inst.qsrc;
  fadType instance_qbulk = inst.qbulk;
  fadType instance_cggb = inst.cggb;
  fadType instance_cgsb = inst.cgsb;
  fadType instance_cgdb = inst.cgdb;
  fadType instance_cdgb = inst.cdgb;
  fadType instance_cdsb = inst.cdsb;
  fadType instance_cddb = inst.cddb;
  fadType instance_cbgb = inst.cbgb;
  fadType instance_cbsb = inst.cbsb;
  fadType instance_cbdb = inst.cbdb;
  fadType instance_cqdb = inst.cqdb;
  fadType instance_cqsb = inst.cqsb;
  fadType instance_cqgb = inst.cqgb;
  fadType instance_cqbb = inst.cqbb;
  fadType instance_gtau = inst.gtau;

  fadType instance_dVgst_dVb = inst.dVgst_dVb;
  fadType instance_dVgst_dVg = inst.dVgst_dVg;
  fadType instance_CoxWL = inst.CoxWL;
  fadType instance_qinv = inst.qinv;
  fadType instance_Cgg = inst.Cgg;
  fadType instance_Cgd = inst.Cgd;
  fadType instance_Cgb = inst.Cgb;
  fadType instance_Cbg = inst.Cbg;
  fadType instance_Cbd = inst.Cbd;
  fadType instance_Cbb = inst.Cbb;
  fadType instance_Csg = inst.Csg;
  fadType instance_Csb = inst.Csb;
  fadType instance_Csd = inst.Csd;
  fadType instance_dDeltaPhi_dVg = inst.dDeltaPhi_dVg;
  fadType instance_dDeltaPhi_dVd = inst.dDeltaPhi_dVd;
  fadType instance_dDeltaPhi_dVb = inst.dDeltaPhi_dVb;
  fadType instance_cd = inst.cd;

  fadType instance_qbs = inst.qbs;
  fadType instance_capbs = inst.capbs;

  fadType instance_qbd = inst.qbd;
  fadType instance_capbd = inst.capbd;
  //fadType instance_nqsMod = inst.nqsMod;
  fadType instance_qcheq = inst.qcheq;

  fadType instance_qgso = inst.qgso;
  fadType instance_qgdo = inst.qgdo;
  fadType instance_cqdef = inst.cqdef;
  fadType instance_qg = inst.qg;
  fadType instance_qd = inst.qd;
  fadType instance_qb = inst.qb;
  fadType instance_qcdump = inst.qcdump;
  fadType instance_Idrain = inst.Idrain;
  fadType instance_drainConductance = inst.drainConductance;
  fadType instance_Isource = inst.Isource;
  fadType instance_sourceConductance = inst.sourceConductance;

  fadType instance_gbs = inst.gbs;
  fadType instance_cbs = inst.cbs;
  fadType instance_gbd = inst.gbd;
  fadType instance_cbd = inst.cbd;

  fadType instance_CAPcggb = inst.CAPcggb;
  fadType instance_CAPcgdb = inst.CAPcgdb;
  fadType instance_CAPcgsb = inst.CAPcgsb;
  fadType instance_CAPcdgb = inst.CAPcdgb;
  fadType instance_CAPcddb = inst.CAPcddb;
  fadType instance_CAPcdsb = inst.CAPcdsb;
  fadType instance_CAPcsgb = inst.CAPcsgb;
  fadType instance_CAPcsdb = inst.CAPcsdb;
  fadType instance_CAPcssb = inst.CAPcssb;
  fadType instance_CAPcbgb = inst.CAPcbgb;
  fadType instance_CAPcbdb = inst.CAPcbdb;
  fadType instance_CAPcbsb = inst.CAPcbsb;

  bool ChargeComputationNeeded;

  bool UIVsuccess = updateIntermediateVars<fadType> 
    (
      inst.getSolverState(),
      inst.getDeviceOptions(),
      sizeDepParams,

      mod.dtype,
      mobMod,
      capMod,

       vbs,
       vds,
       vgs,
       instance_vtm,
       modelPar_jctEmissionCoeff,
       modelPar_ijth,
       modelPar_factor1,
       modelPar_tnom,
       modelPar_tox,
       modelPar_cox,
       modelPar_xpart,
       modelPar_bulkJctBotGradingCoeff,
       modelPar_bulkJctSideGradingCoeff,
       modelPar_bulkJctGateSideGradingCoeff,


        instancePar_sourceArea, 
        instancePar_sourcePerimeter,
        instancePar_drainArea,
        instancePar_drainPerimeter,

       instance_jctTempSatCurDensity,
       instance_jctSidewallTempSatCurDensity,
        instance_vjsm,
        instance_IsEvjsm,
        instance_vjdm,
        instance_IsEvjdm,

        inst.mode,

        instance_T1global,
        instance_thetavth,
        instancePar_temp,
        instance_Vth,
        instance_von,
        instance_dVgs_eff_dVg,
        instance_Vgsteff,
        instance_Abulk,
        instance_ueff,
        instance_Vdsat,
        instance_vdsat,
        instance_Vdseff,
        instance_Gm,
        instance_cdrain,
        instance_gds,
        instance_gm,

        instance_gmbs,
        instance_gbbs,
        instance_gbgs,
        instance_gbds,
        instance_csub,
        instance_qgate,
        instance_qdrn,
        instance_qsrc,
        instance_qbulk,
        instance_cggb,
        instance_cgsb,
        instance_cgdb,
        instance_cdgb,
        instance_cdsb,
        instance_cddb,
        instance_cbgb,
        instance_cbsb,
        instance_cbdb,
        instance_cqdb,
        instance_cqsb,
        instance_cqgb,
        instance_cqbb,
        instance_gtau,

        instance_dVgst_dVb,
        instance_dVgst_dVg,
        instance_CoxWL,
        instance_qinv,
        instance_Cgg,
        instance_Cgd,
        instance_Cgb,
        instance_Cbg,
        instance_Cbd,
        instance_Cbb,
        instance_Csg,
        instance_Csb,
        instance_Csd,
        instance_dDeltaPhi_dVg,
        instance_dDeltaPhi_dVd,
        instance_dDeltaPhi_dVb,
        instance_cd,

        instance_unitAreaJctCapTemp,
        instance_unitLengthGateSidewallJctCapTemp,
        instance_unitLengthSidewallJctCapTemp,

        instance_qbs,
        instance_capbs,

        instance_PhiBTemp,
        instance_PhiBSWTemp,
        instance_PhiBSWGTemp,

        instance_qbd,
        instance_capbd,
        inst.nqsMod, //instance_nqsMod,
        instance_qcheq,

        model_vtm,

        instance_cgdo,
        instance_qgdo,
        instance_cgso,
        instance_qgso,
        instance_cqdef,
        instance_qg,
        instance_qd,
        instance_qb,
        instance_qcdump,

        qdef,
        instance_Idrain,
        instance_drainConductance,
        Vddp,
        instance_Isource,
        instance_sourceConductance,
        Vssp,

      //outputs
       instance_gbs,
       instance_cbs,
       instance_gbd,
       instance_cbd,

       instance_CAPcggb,
       instance_CAPcgdb,
       instance_CAPcgsb,
       instance_CAPcdgb,
       instance_CAPcddb,
       instance_CAPcdsb,
       instance_CAPcsgb,
       instance_CAPcsdb,
       instance_CAPcssb,
       instance_CAPcbgb,
       instance_CAPcbdb,
       instance_CAPcbsb,
       ChargeComputationNeeded
    );

  // loads, and some related calcs
  fadType Gm, Gmbs, FwdSum, RevSum, cdreq, ceqbd, ceqbs, gbbdp, gbbsp;
  fadType gbdpg, gbdpdp, gbdpb, gbdpsp, gbspg, gbspdp, gbspb, gbspsp;
  if (inst.mode >= 0)
  {
    Gm = instance_gm;
    Gmbs = instance_gmbs;
    FwdSum = Gm + Gmbs;
    RevSum = 0.0;

    cdreq =  mod.dtype * (instance_cd);
    ceqbd = -mod.dtype * (instance_csub);

    ceqbs = 0.0;

    gbbdp = -instance_gbds;
    gbbsp = (instance_gbds + instance_gbgs + instance_gbbs);

    gbdpg = instance_gbgs;
    gbdpdp = instance_gbds;
    gbdpb = instance_gbbs;
    gbdpsp = -(gbdpg + gbdpdp + gbdpb);

    gbspg = 0.0;
    gbspdp = 0.0;
    gbspb = 0.0;
    gbspsp = 0.0;
  }
  else
  {
    Gm = -instance_gm;
    Gmbs = -instance_gmbs;
    FwdSum = 0.0;
    RevSum = -(Gm + Gmbs);

    cdreq = -mod.dtype * (instance_cd);
    ceqbs = -mod.dtype * (instance_csub);

    ceqbd = 0.0;

    gbbsp = -instance_gbds;
    gbbdp = (instance_gbds + instance_gbgs + instance_gbbs);

    gbdpg = 0.0;
    gbdpsp = 0.0;
    gbdpb = 0.0;
    gbdpdp = 0.0;

    gbspg = instance_gbgs;
    gbspsp = instance_gbds;
    gbspb = instance_gbbs;
    gbspdp = -(gbspg + gbspsp + gbspb);
  }

  if (mod.dtype > 0)
  {
    ceqbs += (instance_cbs);
    ceqbd += (instance_cbd);
  }
  else
  {
    ceqbs -= (instance_cbs);
    ceqbd -= (instance_cbd);
  }

  int local_Drain = 0;
  int local_Gate = 1;
  int local_Source = 2;
  int local_Bulk = 3;
  int local_DrainPrime = 4;
  int local_SourcePrime = 5;
  int local_Charge = 6;
  int local_Ibs = 7;
  int local_Ids = 8;
  int local_Igs = 9;

  if (inst.drainConductance != 0.0) { dfdp[local_Drain] += instance_Idrain.dx(0)*inst.numberParallel; }
  if (inst.sourceConductance != 0.0) { dfdp[local_Source] += instance_Isource.dx(0)*inst.numberParallel; }

  dfdp[local_Bulk] += (ceqbs.dx(0) + ceqbd.dx(0))*inst.numberParallel;
  dfdp[local_DrainPrime] += (-(ceqbd.dx(0) - cdreq.dx(0))-instance_Idrain.dx(0))*inst.numberParallel;
  dfdp[local_SourcePrime] += (-(cdreq.dx(0) + ceqbs.dx(0))-instance_Isource.dx(0))*inst.numberParallel;

// don't bother with the IC stuff

  fadType ScalingFactor = 1.0e-9; // originally declared in updateIntermediateVars.  This is a copy

  auxChargeCalculations (
       ChargeComputationNeeded,
       inst.mode,
       inst.nqsMod,
       model_vtm,
       sizeDepParams,
    ScalingFactor,
       instance_gtau
      );

  fadType Qeqqg = 0.0;   // gate charge
  fadType Qeqqb = 0.0;   // bulk charge
  fadType Qeqqd = 0.0;   // drain charge
  fadType Qqdef = 0.0;   // nqs-related charge.
  fadType Qqcheq = 0.0;  // nqs-related charge.
  if (inst.model_.dtype > 0)
  {
    Qeqqg  = instance_qg;
    Qeqqb  = instance_qb;
    Qeqqd  = instance_qd;
    Qqdef  = instance_qcdump; // this needs to be fixed...
    Qqcheq = instance_qcheq;
  }
  else  // need to convert these to charges.
  {
    Qeqqg  = -instance_qg;
    Qeqqb  = -instance_qb;
    Qeqqd  = -instance_qd;
    Qqdef  = -instance_qcdump;
    Qqcheq = -instance_qcheq;
  }

  dqdp[local_Gate] += Qeqqg.dx(0)*inst.numberParallel;
  dqdp[local_Bulk] += (Qeqqb.dx(0))*inst.numberParallel;
  dqdp[local_DrainPrime] += (-(-Qeqqd.dx(0)))*inst.numberParallel;
  dqdp[local_SourcePrime] += (-(+ Qeqqg.dx(0) + Qeqqb.dx(0) + Qeqqd.dx(0)))*inst.numberParallel;

  if (inst.nqsMod)
  {
    // 7 equ. for nqs modification. charge equation.
    dqdp[local_Charge] += -(Qqcheq.dx(0) - Qqdef.dx(0))*inst.numberParallel;
  }

  Findices[local_Gate] = inst.li_Gate;
  Findices[local_Source] = inst.li_Source;
  Findices[local_Bulk] = inst.li_Bulk;
  Findices[local_DrainPrime] = inst.li_DrainPrime;
  Findices[local_SourcePrime] = inst.li_SourcePrime;
  Findices[local_Charge] = inst.li_Charge;
  Findices[local_Ibs] = inst.li_Ibs;
  Findices[local_Ids] = inst.li_Ids;
  Findices[local_Igs] = inst.li_Igs;

  Qindices[local_Drain] = inst.li_Drain;
  Qindices[local_Gate] = inst.li_Gate;
  Qindices[local_Source] = inst.li_Source;
  Qindices[local_Bulk] = inst.li_Bulk;
  Qindices[local_DrainPrime] = inst.li_DrainPrime;
  Qindices[local_SourcePrime] = inst.li_SourcePrime;
  Qindices[local_Charge] = inst.li_Charge;
  Qindices[local_Ibs] = inst.li_Ibs;
  Qindices[local_Ids] = inst.li_Ids;
  Qindices[local_Igs] = inst.li_Igs;
}

//-----------------------------------------------------------------------------
// Function      : bsim3ModelSensitivity::operator
// Purpose       : produces df/dp and dq/dp, where p= any bsim3 instance parameter.  
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/30/2018
//-----------------------------------------------------------------------------
void bsim3ModelSensitivity::operator()(
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

  std::vector<Instance *>::const_iterator inTmp = mod.instanceContainer.begin(); 
  const Instance & instTmp = *((*inTmp));

  dfdp.resize(10*sizeInstance);
  dqdp.resize(10*sizeInstance);
  Findices.resize(10*sizeInstance);
  Qindices.resize(10*sizeInstance);

  unordered_map <std::string,fadType*,HashNoCase,EqualNoCase> modParamMap;

  // model params
  // sensitivity supported, (doubles):
  fadType modelPar_tox=mod.tox; 
  bool modelPar_given_TOX=mod.given("TOX"); 
  modParamMap["TOX"] = &modelPar_tox;
  
  fadType modelPar_toxm=mod.toxm; 
  bool modelPar_given_TOXM=mod.given("TOXM"); 
  modParamMap["TOXM"] = &modelPar_toxm;
  
  fadType modelPar_cdsc=mod.cdsc; 
  bool modelPar_given_CDSC=mod.given("CDSC"); 
  modParamMap["CDSC"] = &modelPar_cdsc;
  
  fadType modelPar_cdscb=mod.cdscb; 
  bool modelPar_given_CDSCB=mod.given("CDSCB"); 
  modParamMap["CDSCB"] = &modelPar_cdscb;
  
  fadType modelPar_cdscd=mod.cdscd; 
  bool modelPar_given_CDSCD=mod.given("CDSCD"); 
  modParamMap["CDSCD"] = &modelPar_cdscd;
  
  fadType modelPar_cit=mod.cit; 
  bool modelPar_given_CIT=mod.given("CIT"); 
  modParamMap["CIT"] = &modelPar_cit;
  
  fadType modelPar_nfactor=mod.nfactor; 
  bool modelPar_given_NFACTOR=mod.given("NFACTOR"); 
  modParamMap["NFACTOR"] = &modelPar_nfactor;
  
  fadType modelPar_xj=mod.xj; 
  bool modelPar_given_XJ=mod.given("XJ"); 
  modParamMap["XJ"] = &modelPar_xj;
  
  fadType modelPar_vsat=mod.vsat; 
  bool modelPar_given_VSAT=mod.given("VSAT"); 
  modParamMap["VSAT"] = &modelPar_vsat;
  
  fadType modelPar_at=mod.at; 
  bool modelPar_given_AT=mod.given("AT"); 
  modParamMap["AT"] = &modelPar_at;
  
  fadType modelPar_a0=mod.a0; 
  bool modelPar_given_A0=mod.given("A0"); 
  modParamMap["A0"] = &modelPar_a0;
  
  fadType modelPar_ags=mod.ags; 
  bool modelPar_given_AGS=mod.given("AGS"); 
  modParamMap["AGS"] = &modelPar_ags;
  
  fadType modelPar_a1=mod.a1; 
  bool modelPar_given_A1=mod.given("A1"); 
  modParamMap["A1"] = &modelPar_a1;
  
  fadType modelPar_a2=mod.a2; 
  bool modelPar_given_A2=mod.given("A2"); 
  modParamMap["A2"] = &modelPar_a2;
  
  fadType modelPar_keta=mod.keta; 
  bool modelPar_given_KETA=mod.given("KETA"); 
  modParamMap["KETA"] = &modelPar_keta;
  
  fadType modelPar_nsub=mod.nsub; 
  bool modelPar_given_NSUB=mod.given("NSUB"); 
  modParamMap["NSUB"] = &modelPar_nsub;
  
  fadType modelPar_npeak=mod.npeak; 
  bool modelPar_given_NCH=mod.given("NCH"); 
  modParamMap["NCH"] = &modelPar_npeak;
  
  fadType modelPar_ngate=mod.ngate; 
  bool modelPar_given_NGATE=mod.given("NGATE"); 
  modParamMap["NGATE"] = &modelPar_ngate;
  
  fadType modelPar_gamma1=mod.gamma1; 
  bool modelPar_given_GAMMA1=mod.given("GAMMA1"); 
  modParamMap["GAMMA1"] = &modelPar_gamma1;
  
  fadType modelPar_gamma2=mod.gamma2; 
  bool modelPar_given_GAMMA2=mod.given("GAMMA2"); 
  modParamMap["GAMMA2"] = &modelPar_gamma2;
  
  fadType modelPar_vbx=mod.vbx; 
  bool modelPar_given_VBX=mod.given("VBX"); 
  modParamMap["VBX"] = &modelPar_vbx;
  
  fadType modelPar_vbm=mod.vbm; 
  bool modelPar_given_VBM=mod.given("VBM"); 
  modParamMap["VBM"] = &modelPar_vbm;
  
  fadType modelPar_xt=mod.xt; 
  bool modelPar_given_XT=mod.given("XT"); 
  modParamMap["XT"] = &modelPar_xt;
  
  fadType modelPar_k1=mod.k1; 
  bool modelPar_given_K1=mod.given("K1"); 
  modParamMap["K1"] = &modelPar_k1;
  
  fadType modelPar_kt1=mod.kt1; 
  bool modelPar_given_KT1=mod.given("KT1"); 
  modParamMap["KT1"] = &modelPar_kt1;
  
  fadType modelPar_kt1l=mod.kt1l; 
  bool modelPar_given_KT1L=mod.given("KT1L"); 
  modParamMap["KT1L"] = &modelPar_kt1l;
  
  fadType modelPar_kt2=mod.kt2; 
  bool modelPar_given_KT2=mod.given("KT2"); 
  modParamMap["KT2"] = &modelPar_kt2;
  
  fadType modelPar_k2=mod.k2; 
  bool modelPar_given_K2=mod.given("K2"); 
  modParamMap["K2"] = &modelPar_k2;
  
  fadType modelPar_k3=mod.k3; 
  bool modelPar_given_K3=mod.given("K3"); 
  modParamMap["K3"] = &modelPar_k3;
  
  fadType modelPar_k3b=mod.k3b; 
  bool modelPar_given_K3B=mod.given("K3B"); 
  modParamMap["K3B"] = &modelPar_k3b;
  
  fadType modelPar_w0=mod.w0; 
  bool modelPar_given_W0=mod.given("W0"); 
  modParamMap["W0"] = &modelPar_w0;
  
  fadType modelPar_nlx=mod.nlx; 
  bool modelPar_given_NLX=mod.given("NLX"); 
  modParamMap["NLX"] = &modelPar_nlx;
  
  fadType modelPar_dvt0=mod.dvt0; 
  bool modelPar_given_DVT0=mod.given("DVT0"); 
  modParamMap["DVT0"] = &modelPar_dvt0;
  
  fadType modelPar_dvt1=mod.dvt1; 
  bool modelPar_given_DVT1=mod.given("DVT1"); 
  modParamMap["DVT1"] = &modelPar_dvt1;
  
  fadType modelPar_dvt2=mod.dvt2; 
  bool modelPar_given_DVT2=mod.given("DVT2"); 
  modParamMap["DVT2"] = &modelPar_dvt2;
  
  fadType modelPar_dvt0w=mod.dvt0w; 
  bool modelPar_given_DVT0W=mod.given("DVT0W"); 
  modParamMap["DVT0W"] = &modelPar_dvt0w;
  
  fadType modelPar_dvt1w=mod.dvt1w; 
  bool modelPar_given_DVT1W=mod.given("DVT1W"); 
  modParamMap["DVT1W"] = &modelPar_dvt1w;
  
  fadType modelPar_dvt2w=mod.dvt2w; 
  bool modelPar_given_DVT2W=mod.given("DVT2W"); 
  modParamMap["DVT2W"] = &modelPar_dvt2w;
  
  fadType modelPar_drout=mod.drout; 
  bool modelPar_given_DROUT=mod.given("DROUT"); 
  modParamMap["DROUT"] = &modelPar_drout;
  
  fadType modelPar_dsub=mod.dsub; 
  bool modelPar_given_DSUB=mod.given("DSUB"); 
  modParamMap["DSUB"] = &modelPar_dsub;
  
  fadType modelPar_vth0=mod.vth0; 
  bool modelPar_given_VTH0=mod.given("VTH0"); 
  modParamMap["VTH0"] = &modelPar_vth0;
  
  fadType modelPar_ua=mod.ua; 
  bool modelPar_given_UA=mod.given("UA"); 
  modParamMap["UA"] = &modelPar_ua;
  
  fadType modelPar_ua1=mod.ua1; 
  bool modelPar_given_UA1=mod.given("UA1"); 
  modParamMap["UA1"] = &modelPar_ua1;
  
  fadType modelPar_ub=mod.ub; 
  bool modelPar_given_UB=mod.given("UB"); 
  modParamMap["UB"] = &modelPar_ub;
  
  fadType modelPar_ub1=mod.ub1; 
  bool modelPar_given_UB1=mod.given("UB1"); 
  modParamMap["UB1"] = &modelPar_ub1;
  
  fadType modelPar_uc=mod.uc; 
  bool modelPar_given_UC=mod.given("UC"); 
  modParamMap["UC"] = &modelPar_uc;
  
  fadType modelPar_uc1=mod.uc1; 
  bool modelPar_given_UC1=mod.given("UC1"); 
  modParamMap["UC1"] = &modelPar_uc1;
  
  fadType modelPar_u0=mod.u0; 
  bool modelPar_given_U0=mod.given("U0"); 
  modParamMap["U0"] = &modelPar_u0;
  
  fadType modelPar_ute=mod.ute; 
  bool modelPar_given_UTE=mod.given("UTE"); 
  modParamMap["UTE"] = &modelPar_ute;
  
  fadType modelPar_voff=mod.voff; 
  bool modelPar_given_VOFF=mod.given("VOFF"); 
  modParamMap["VOFF"] = &modelPar_voff;
  
  fadType modelPar_rdsw=mod.rdsw; 
  bool modelPar_given_RDSW=mod.given("RDSW"); 
  modParamMap["RDSW"] = &modelPar_rdsw;
  
  fadType modelPar_prwg=mod.prwg; 
  bool modelPar_given_PRWG=mod.given("PRWG"); 
  modParamMap["PRWG"] = &modelPar_prwg;
  
  fadType modelPar_prwb=mod.prwb; 
  bool modelPar_given_PRWB=mod.given("PRWB"); 
  modParamMap["PRWB"] = &modelPar_prwb;
  
  fadType modelPar_prt=mod.prt; 
  bool modelPar_given_PRT=mod.given("PRT"); 
  modParamMap["PRT"] = &modelPar_prt;
  
  fadType modelPar_eta0=mod.eta0; 
  bool modelPar_given_ETA0=mod.given("ETA0"); 
  modParamMap["ETA0"] = &modelPar_eta0;
  
  fadType modelPar_etab=mod.etab; 
  bool modelPar_given_ETAB=mod.given("ETAB"); 
  modParamMap["ETAB"] = &modelPar_etab;
  
  fadType modelPar_pclm=mod.pclm; 
  bool modelPar_given_PCLM=mod.given("PCLM"); 
  modParamMap["PCLM"] = &modelPar_pclm;
  
  fadType modelPar_pdibl1=mod.pdibl1; 
  bool modelPar_given_PDIBLC1=mod.given("PDIBLC1"); 
  modParamMap["PDIBLC1"] = &modelPar_pdibl1;
  
  fadType modelPar_pdibl2=mod.pdibl2; 
  bool modelPar_given_PDIBLC2=mod.given("PDIBLC2"); 
  modParamMap["PDIBLC2"] = &modelPar_pdibl2;
  
  fadType modelPar_pdiblb=mod.pdiblb; 
  bool modelPar_given_PDIBLCB=mod.given("PDIBLCB"); 
  modParamMap["PDIBLCB"] = &modelPar_pdiblb;
  
  fadType modelPar_pscbe1=mod.pscbe1; 
  bool modelPar_given_PSCBE1=mod.given("PSCBE1"); 
  modParamMap["PSCBE1"] = &modelPar_pscbe1;
  
  fadType modelPar_pscbe2=mod.pscbe2; 
  bool modelPar_given_PSCBE2=mod.given("PSCBE2"); 
  modParamMap["PSCBE2"] = &modelPar_pscbe2;
  
  fadType modelPar_pvag=mod.pvag; 
  bool modelPar_given_PVAG=mod.given("PVAG"); 
  modParamMap["PVAG"] = &modelPar_pvag;
  
  fadType modelPar_delta=mod.delta; 
  bool modelPar_given_DELTA=mod.given("DELTA"); 
  modParamMap["DELTA"] = &modelPar_delta;
  
  fadType modelPar_wr=mod.wr; 
  bool modelPar_given_WR=mod.given("WR"); 
  modParamMap["WR"] = &modelPar_wr;
  
  fadType modelPar_dwg=mod.dwg; 
  bool modelPar_given_DWG=mod.given("DWG"); 
  modParamMap["DWG"] = &modelPar_dwg;
  
  fadType modelPar_dwb=mod.dwb; 
  bool modelPar_given_DWB=mod.given("DWB"); 
  modParamMap["DWB"] = &modelPar_dwb;
  
  fadType modelPar_b0=mod.b0; 
  bool modelPar_given_B0=mod.given("B0"); 
  modParamMap["B0"] = &modelPar_b0;
  
  fadType modelPar_b1=mod.b1; 
  bool modelPar_given_B1=mod.given("B1"); 
  modParamMap["B1"] = &modelPar_b1;
  
  fadType modelPar_alpha0=mod.alpha0; 
  bool modelPar_given_ALPHA0=mod.given("ALPHA0"); 
  modParamMap["ALPHA0"] = &modelPar_alpha0;
  
  fadType modelPar_alpha1=mod.alpha1; 
  bool modelPar_given_ALPHA1=mod.given("ALPHA1"); 
  modParamMap["ALPHA1"] = &modelPar_alpha1;
  
  fadType modelPar_beta0=mod.beta0; 
  bool modelPar_given_BETA0=mod.given("BETA0"); 
  modParamMap["BETA0"] = &modelPar_beta0;
  
  fadType modelPar_ijth=mod.ijth; 
  bool modelPar_given_IJTH=mod.given("IJTH"); 
  modParamMap["IJTH"] = &modelPar_ijth;
  
  fadType modelPar_vfb=mod.vfb; 
  bool modelPar_given_VFB=mod.given("VFB"); 
  modParamMap["VFB"] = &modelPar_vfb;
  
  fadType modelPar_elm=mod.elm; 
  bool modelPar_given_ELM=mod.given("ELM"); 
  modParamMap["ELM"] = &modelPar_elm;
  
  fadType modelPar_cgsl=mod.cgsl; 
  bool modelPar_given_CGSL=mod.given("CGSL"); 
  modParamMap["CGSL"] = &modelPar_cgsl;
  
  fadType modelPar_cgdl=mod.cgdl; 
  bool modelPar_given_CGDL=mod.given("CGDL"); 
  modParamMap["CGDL"] = &modelPar_cgdl;
  
  fadType modelPar_ckappa=mod.ckappa; 
  bool modelPar_given_CKAPPA=mod.given("CKAPPA"); 
  modParamMap["CKAPPA"] = &modelPar_ckappa;
  
  fadType modelPar_cf=mod.cf; 
  bool modelPar_given_CF=mod.given("CF"); 
  modParamMap["CF"] = &modelPar_cf;
  
  fadType modelPar_vfbcv=mod.vfbcv; 
  bool modelPar_given_VFBCV=mod.given("VFBCV"); 
  modParamMap["VFBCV"] = &modelPar_vfbcv;
  
  fadType modelPar_clc=mod.clc; 
  bool modelPar_given_CLC=mod.given("CLC"); 
  modParamMap["CLC"] = &modelPar_clc;
  
  fadType modelPar_cle=mod.cle; 
  bool modelPar_given_CLE=mod.given("CLE"); 
  modParamMap["CLE"] = &modelPar_cle;
  
  fadType modelPar_dwc=mod.dwc; 
  bool modelPar_given_DWC=mod.given("DWC"); 
  modParamMap["DWC"] = &modelPar_dwc;
  
  fadType modelPar_dlc=mod.dlc; 
  bool modelPar_given_DLC=mod.given("DLC"); 
  modParamMap["DLC"] = &modelPar_dlc;
  
  fadType modelPar_noff=mod.noff; 
  bool modelPar_given_NOFF=mod.given("NOFF"); 
  modParamMap["NOFF"] = &modelPar_noff;
  
  fadType modelPar_voffcv=mod.voffcv; 
  bool modelPar_given_VOFFCV=mod.given("VOFFCV"); 
  modParamMap["VOFFCV"] = &modelPar_voffcv;
  
  fadType modelPar_acde=mod.acde; 
  bool modelPar_given_ACDE=mod.given("ACDE"); 
  modParamMap["ACDE"] = &modelPar_acde;
  
  fadType modelPar_moin=mod.moin; 
  bool modelPar_given_MOIN=mod.given("MOIN"); 
  modParamMap["MOIN"] = &modelPar_moin;
  
  fadType modelPar_tcj=mod.tcj; 
  bool modelPar_given_TCJ=mod.given("TCJ"); 
  modParamMap["TCJ"] = &modelPar_tcj;
  
  fadType modelPar_tcjsw=mod.tcjsw; 
  bool modelPar_given_TCJSW=mod.given("TCJSW"); 
  modParamMap["TCJSW"] = &modelPar_tcjsw;
  
  fadType modelPar_tcjswg=mod.tcjswg; 
  bool modelPar_given_TCJSWG=mod.given("TCJSWG"); 
  modParamMap["TCJSWG"] = &modelPar_tcjswg;
  
  fadType modelPar_tpb=mod.tpb; 
  bool modelPar_given_TPB=mod.given("TPB"); 
  modParamMap["TPB"] = &modelPar_tpb;
  
  fadType modelPar_tpbsw=mod.tpbsw; 
  bool modelPar_given_TPBSW=mod.given("TPBSW"); 
  modParamMap["TPBSW"] = &modelPar_tpbsw;
  
  fadType modelPar_tpbswg=mod.tpbswg; 
  bool modelPar_given_TPBSWG=mod.given("TPBSWG"); 
  modParamMap["TPBSWG"] = &modelPar_tpbswg;
  
  fadType modelPar_lcdsc=mod.lcdsc; 
  bool modelPar_given_LCDSC=mod.given("LCDSC"); 
  modParamMap["LCDSC"] = &modelPar_lcdsc;
  
  fadType modelPar_lcdscb=mod.lcdscb; 
  bool modelPar_given_LCDSCB=mod.given("LCDSCB"); 
  modParamMap["LCDSCB"] = &modelPar_lcdscb;
  
  fadType modelPar_lcdscd=mod.lcdscd; 
  bool modelPar_given_LCDSCD=mod.given("LCDSCD"); 
  modParamMap["LCDSCD"] = &modelPar_lcdscd;
  
  fadType modelPar_lcit=mod.lcit; 
  bool modelPar_given_LCIT=mod.given("LCIT"); 
  modParamMap["LCIT"] = &modelPar_lcit;
  
  fadType modelPar_lnfactor=mod.lnfactor; 
  bool modelPar_given_LNFACTOR=mod.given("LNFACTOR"); 
  modParamMap["LNFACTOR"] = &modelPar_lnfactor;
  
  fadType modelPar_lxj=mod.lxj; 
  bool modelPar_given_LXJ=mod.given("LXJ"); 
  modParamMap["LXJ"] = &modelPar_lxj;
  
  fadType modelPar_lvsat=mod.lvsat; 
  bool modelPar_given_LVSAT=mod.given("LVSAT"); 
  modParamMap["LVSAT"] = &modelPar_lvsat;
  
  fadType modelPar_lat=mod.lat; 
  bool modelPar_given_LAT=mod.given("LAT"); 
  modParamMap["LAT"] = &modelPar_lat;
  
  fadType modelPar_la0=mod.la0; 
  bool modelPar_given_LA0=mod.given("LA0"); 
  modParamMap["LA0"] = &modelPar_la0;
  
  fadType modelPar_lags=mod.lags; 
  bool modelPar_given_LAGS=mod.given("LAGS"); 
  modParamMap["LAGS"] = &modelPar_lags;
  
  fadType modelPar_la1=mod.la1; 
  bool modelPar_given_LA1=mod.given("LA1"); 
  modParamMap["LA1"] = &modelPar_la1;
  
  fadType modelPar_la2=mod.la2; 
  bool modelPar_given_LA2=mod.given("LA2"); 
  modParamMap["LA2"] = &modelPar_la2;
  
  fadType modelPar_lketa=mod.lketa; 
  bool modelPar_given_LKETA=mod.given("LKETA"); 
  modParamMap["LKETA"] = &modelPar_lketa;
  
  fadType modelPar_lnsub=mod.lnsub; 
  bool modelPar_given_LNSUB=mod.given("LNSUB"); 
  modParamMap["LNSUB"] = &modelPar_lnsub;
  
  fadType modelPar_lnpeak=mod.lnpeak; 
  bool modelPar_given_LNCH=mod.given("LNCH"); 
  modParamMap["LNCH"] = &modelPar_lnpeak;
  
  fadType modelPar_lngate=mod.lngate; 
  bool modelPar_given_LNGATE=mod.given("LNGATE"); 
  modParamMap["LNGATE"] = &modelPar_lngate;
  
  fadType modelPar_lgamma1=mod.lgamma1; 
  bool modelPar_given_LGAMMA1=mod.given("LGAMMA1"); 
  modParamMap["LGAMMA1"] = &modelPar_lgamma1;
  
  fadType modelPar_lgamma2=mod.lgamma2; 
  bool modelPar_given_LGAMMA2=mod.given("LGAMMA2"); 
  modParamMap["LGAMMA2"] = &modelPar_lgamma2;
  
  fadType modelPar_lvbx=mod.lvbx; 
  bool modelPar_given_LVBX=mod.given("LVBX"); 
  modParamMap["LVBX"] = &modelPar_lvbx;
  
  fadType modelPar_lvbm=mod.lvbm; 
  bool modelPar_given_LVBM=mod.given("LVBM"); 
  modParamMap["LVBM"] = &modelPar_lvbm;
  
  fadType modelPar_lxt=mod.lxt; 
  bool modelPar_given_LXT=mod.given("LXT"); 
  modParamMap["LXT"] = &modelPar_lxt;
  
  fadType modelPar_lk1=mod.lk1; 
  bool modelPar_given_LK1=mod.given("LK1"); 
  modParamMap["LK1"] = &modelPar_lk1;
  
  fadType modelPar_lkt1=mod.lkt1; 
  bool modelPar_given_LKT1=mod.given("LKT1"); 
  modParamMap["LKT1"] = &modelPar_lkt1;
  
  fadType modelPar_lkt1l=mod.lkt1l; 
  bool modelPar_given_LKT1L=mod.given("LKT1L"); 
  modParamMap["LKT1L"] = &modelPar_lkt1l;
  
  fadType modelPar_lkt2=mod.lkt2; 
  bool modelPar_given_LKT2=mod.given("LKT2"); 
  modParamMap["LKT2"] = &modelPar_lkt2;
  
  fadType modelPar_lk2=mod.lk2; 
  bool modelPar_given_LK2=mod.given("LK2"); 
  modParamMap["LK2"] = &modelPar_lk2;
  
  fadType modelPar_lk3=mod.lk3; 
  bool modelPar_given_LK3=mod.given("LK3"); 
  modParamMap["LK3"] = &modelPar_lk3;
  
  fadType modelPar_lk3b=mod.lk3b; 
  bool modelPar_given_LK3B=mod.given("LK3B"); 
  modParamMap["LK3B"] = &modelPar_lk3b;
  
  fadType modelPar_lw0=mod.lw0; 
  bool modelPar_given_LW0=mod.given("LW0"); 
  modParamMap["LW0"] = &modelPar_lw0;
  
  fadType modelPar_lnlx=mod.lnlx; 
  bool modelPar_given_LNLX=mod.given("LNLX"); 
  modParamMap["LNLX"] = &modelPar_lnlx;
  
  fadType modelPar_ldvt0=mod.ldvt0; 
  bool modelPar_given_LDVT0=mod.given("LDVT0"); 
  modParamMap["LDVT0"] = &modelPar_ldvt0;
  
  fadType modelPar_ldvt1=mod.ldvt1; 
  bool modelPar_given_LDVT1=mod.given("LDVT1"); 
  modParamMap["LDVT1"] = &modelPar_ldvt1;
  
  fadType modelPar_ldvt2=mod.ldvt2; 
  bool modelPar_given_LDVT2=mod.given("LDVT2"); 
  modParamMap["LDVT2"] = &modelPar_ldvt2;
  
  fadType modelPar_ldvt0w=mod.ldvt0w; 
  bool modelPar_given_LDVT0W=mod.given("LDVT0W"); 
  modParamMap["LDVT0W"] = &modelPar_ldvt0w;
  
  fadType modelPar_ldvt1w=mod.ldvt1w; 
  bool modelPar_given_LDVT1W=mod.given("LDVT1W"); 
  modParamMap["LDVT1W"] = &modelPar_ldvt1w;
  
  fadType modelPar_ldvt2w=mod.ldvt2w; 
  bool modelPar_given_LDVT2W=mod.given("LDVT2W"); 
  modParamMap["LDVT2W"] = &modelPar_ldvt2w;
  
  fadType modelPar_ldrout=mod.ldrout; 
  bool modelPar_given_LDROUT=mod.given("LDROUT"); 
  modParamMap["LDROUT"] = &modelPar_ldrout;
  
  fadType modelPar_ldsub=mod.ldsub; 
  bool modelPar_given_LDSUB=mod.given("LDSUB"); 
  modParamMap["LDSUB"] = &modelPar_ldsub;
  
  fadType modelPar_lvth0=mod.lvth0; 
  bool modelPar_given_LVTH0=mod.given("LVTH0"); 
  modParamMap["LVTH0"] = &modelPar_lvth0;
  
  fadType modelPar_lua=mod.lua; 
  bool modelPar_given_LUA=mod.given("LUA"); 
  modParamMap["LUA"] = &modelPar_lua;
  
  fadType modelPar_lua1=mod.lua1; 
  bool modelPar_given_LUA1=mod.given("LUA1"); 
  modParamMap["LUA1"] = &modelPar_lua1;
  
  fadType modelPar_lub=mod.lub; 
  bool modelPar_given_LUB=mod.given("LUB"); 
  modParamMap["LUB"] = &modelPar_lub;
  
  fadType modelPar_lub1=mod.lub1; 
  bool modelPar_given_LUB1=mod.given("LUB1"); 
  modParamMap["LUB1"] = &modelPar_lub1;
  
  fadType modelPar_luc=mod.luc; 
  bool modelPar_given_LUC=mod.given("LUC"); 
  modParamMap["LUC"] = &modelPar_luc;
  
  fadType modelPar_luc1=mod.luc1; 
  bool modelPar_given_LUC1=mod.given("LUC1"); 
  modParamMap["LUC1"] = &modelPar_luc1;
  
  fadType modelPar_lu0=mod.lu0; 
  bool modelPar_given_LU0=mod.given("LU0"); 
  modParamMap["LU0"] = &modelPar_lu0;
  
  fadType modelPar_lute=mod.lute; 
  bool modelPar_given_LUTE=mod.given("LUTE"); 
  modParamMap["LUTE"] = &modelPar_lute;
  
  fadType modelPar_lvoff=mod.lvoff; 
  bool modelPar_given_LVOFF=mod.given("LVOFF"); 
  modParamMap["LVOFF"] = &modelPar_lvoff;
  
  fadType modelPar_lrdsw=mod.lrdsw; 
  bool modelPar_given_LRDSW=mod.given("LRDSW"); 
  modParamMap["LRDSW"] = &modelPar_lrdsw;
  
  fadType modelPar_lprwg=mod.lprwg; 
  bool modelPar_given_LPRWG=mod.given("LPRWG"); 
  modParamMap["LPRWG"] = &modelPar_lprwg;
  
  fadType modelPar_lprwb=mod.lprwb; 
  bool modelPar_given_LPRWB=mod.given("LPRWB"); 
  modParamMap["LPRWB"] = &modelPar_lprwb;
  
  fadType modelPar_lprt=mod.lprt; 
  bool modelPar_given_LPRT=mod.given("LPRT"); 
  modParamMap["LPRT"] = &modelPar_lprt;
  
  fadType modelPar_leta0=mod.leta0; 
  bool modelPar_given_LETA0=mod.given("LETA0"); 
  modParamMap["LETA0"] = &modelPar_leta0;
  
  fadType modelPar_letab=mod.letab; 
  bool modelPar_given_LETAB=mod.given("LETAB"); 
  modParamMap["LETAB"] = &modelPar_letab;
  
  fadType modelPar_lpclm=mod.lpclm; 
  bool modelPar_given_LPCLM=mod.given("LPCLM"); 
  modParamMap["LPCLM"] = &modelPar_lpclm;
  
  fadType modelPar_lpdibl1=mod.lpdibl1; 
  bool modelPar_given_LPDIBLC1=mod.given("LPDIBLC1"); 
  modParamMap["LPDIBLC1"] = &modelPar_lpdibl1;
  
  fadType modelPar_lpdibl2=mod.lpdibl2; 
  bool modelPar_given_LPDIBLC2=mod.given("LPDIBLC2"); 
  modParamMap["LPDIBLC2"] = &modelPar_lpdibl2;
  
  fadType modelPar_lpdiblb=mod.lpdiblb; 
  bool modelPar_given_LPDIBLCB=mod.given("LPDIBLCB"); 
  modParamMap["LPDIBLCB"] = &modelPar_lpdiblb;
  
  fadType modelPar_lpscbe1=mod.lpscbe1; 
  bool modelPar_given_LPSCBE1=mod.given("LPSCBE1"); 
  modParamMap["LPSCBE1"] = &modelPar_lpscbe1;
  
  fadType modelPar_lpscbe2=mod.lpscbe2; 
  bool modelPar_given_LPSCBE2=mod.given("LPSCBE2"); 
  modParamMap["LPSCBE2"] = &modelPar_lpscbe2;
  
  fadType modelPar_lpvag=mod.lpvag; 
  bool modelPar_given_LPVAG=mod.given("LPVAG"); 
  modParamMap["LPVAG"] = &modelPar_lpvag;
  
  fadType modelPar_ldelta=mod.ldelta; 
  bool modelPar_given_LDELTA=mod.given("LDELTA"); 
  modParamMap["LDELTA"] = &modelPar_ldelta;
  
  fadType modelPar_lwr=mod.lwr; 
  bool modelPar_given_LWR=mod.given("LWR"); 
  modParamMap["LWR"] = &modelPar_lwr;
  
  fadType modelPar_ldwg=mod.ldwg; 
  bool modelPar_given_LDWG=mod.given("LDWG"); 
  modParamMap["LDWG"] = &modelPar_ldwg;
  
  fadType modelPar_ldwb=mod.ldwb; 
  bool modelPar_given_LDWB=mod.given("LDWB"); 
  modParamMap["LDWB"] = &modelPar_ldwb;
  
  fadType modelPar_lb0=mod.lb0; 
  bool modelPar_given_LB0=mod.given("LB0"); 
  modParamMap["LB0"] = &modelPar_lb0;
  
  fadType modelPar_lb1=mod.lb1; 
  bool modelPar_given_LB1=mod.given("LB1"); 
  modParamMap["LB1"] = &modelPar_lb1;
  
  fadType modelPar_lalpha0=mod.lalpha0; 
  bool modelPar_given_LALPHA0=mod.given("LALPHA0"); 
  modParamMap["LALPHA0"] = &modelPar_lalpha0;
  
  fadType modelPar_lalpha1=mod.lalpha1; 
  bool modelPar_given_LALPHA1=mod.given("LALPHA1"); 
  modParamMap["LALPHA1"] = &modelPar_lalpha1;
  
  fadType modelPar_lbeta0=mod.lbeta0; 
  bool modelPar_given_LBETA0=mod.given("LBETA0"); 
  modParamMap["LBETA0"] = &modelPar_lbeta0;
  
  fadType modelPar_lvfb=mod.lvfb; 
  bool modelPar_given_LVFB=mod.given("LVFB"); 
  modParamMap["LVFB"] = &modelPar_lvfb;
  
  fadType modelPar_lelm=mod.lelm; 
  bool modelPar_given_LELM=mod.given("LELM"); 
  modParamMap["LELM"] = &modelPar_lelm;
  
  fadType modelPar_lcgsl=mod.lcgsl; 
  bool modelPar_given_LCGSL=mod.given("LCGSL"); 
  modParamMap["LCGSL"] = &modelPar_lcgsl;
  
  fadType modelPar_lcgdl=mod.lcgdl; 
  bool modelPar_given_LCGDL=mod.given("LCGDL"); 
  modParamMap["LCGDL"] = &modelPar_lcgdl;
  
  fadType modelPar_lckappa=mod.lckappa; 
  bool modelPar_given_LCKAPPA=mod.given("LCKAPPA"); 
  modParamMap["LCKAPPA"] = &modelPar_lckappa;
  
  fadType modelPar_lcf=mod.lcf; 
  bool modelPar_given_LCF=mod.given("LCF"); 
  modParamMap["LCF"] = &modelPar_lcf;
  
  fadType modelPar_lclc=mod.lclc; 
  bool modelPar_given_LCLC=mod.given("LCLC"); 
  modParamMap["LCLC"] = &modelPar_lclc;
  
  fadType modelPar_lcle=mod.lcle; 
  bool modelPar_given_LCLE=mod.given("LCLE"); 
  modParamMap["LCLE"] = &modelPar_lcle;
  
  fadType modelPar_lvfbcv=mod.lvfbcv; 
  bool modelPar_given_LVFBCV=mod.given("LVFBCV"); 
  modParamMap["LVFBCV"] = &modelPar_lvfbcv;
  
  fadType modelPar_lnoff=mod.lnoff; 
  bool modelPar_given_LNOFF=mod.given("LNOFF"); 
  modParamMap["LNOFF"] = &modelPar_lnoff;
  
  fadType modelPar_lvoffcv=mod.lvoffcv; 
  bool modelPar_given_LVOFFCV=mod.given("LVOFFCV"); 
  modParamMap["LVOFFCV"] = &modelPar_lvoffcv;
  
  fadType modelPar_lacde=mod.lacde; 
  bool modelPar_given_LACDE=mod.given("LACDE"); 
  modParamMap["LACDE"] = &modelPar_lacde;
  
  fadType modelPar_lmoin=mod.lmoin; 
  bool modelPar_given_LMOIN=mod.given("LMOIN"); 
  modParamMap["LMOIN"] = &modelPar_lmoin;
  
  fadType modelPar_wcdsc=mod.wcdsc; 
  bool modelPar_given_WCDSC=mod.given("WCDSC"); 
  modParamMap["WCDSC"] = &modelPar_wcdsc;
  
  fadType modelPar_wcdscb=mod.wcdscb; 
  bool modelPar_given_WCDSCB=mod.given("WCDSCB"); 
  modParamMap["WCDSCB"] = &modelPar_wcdscb;
  
  fadType modelPar_wcdscd=mod.wcdscd; 
  bool modelPar_given_WCDSCD=mod.given("WCDSCD"); 
  modParamMap["WCDSCD"] = &modelPar_wcdscd;
  
  fadType modelPar_wcit=mod.wcit; 
  bool modelPar_given_WCIT=mod.given("WCIT"); 
  modParamMap["WCIT"] = &modelPar_wcit;
  
  fadType modelPar_wnfactor=mod.wnfactor; 
  bool modelPar_given_WNFACTOR=mod.given("WNFACTOR"); 
  modParamMap["WNFACTOR"] = &modelPar_wnfactor;
  
  fadType modelPar_wxj=mod.wxj; 
  bool modelPar_given_WXJ=mod.given("WXJ"); 
  modParamMap["WXJ"] = &modelPar_wxj;
  
  fadType modelPar_wvsat=mod.wvsat; 
  bool modelPar_given_WVSAT=mod.given("WVSAT"); 
  modParamMap["WVSAT"] = &modelPar_wvsat;
  
  fadType modelPar_wat=mod.wat; 
  bool modelPar_given_WAT=mod.given("WAT"); 
  modParamMap["WAT"] = &modelPar_wat;
  
  fadType modelPar_wa0=mod.wa0; 
  bool modelPar_given_WA0=mod.given("WA0"); 
  modParamMap["WA0"] = &modelPar_wa0;
  
  fadType modelPar_wags=mod.wags; 
  bool modelPar_given_WAGS=mod.given("WAGS"); 
  modParamMap["WAGS"] = &modelPar_wags;
  
  fadType modelPar_wa1=mod.wa1; 
  bool modelPar_given_WA1=mod.given("WA1"); 
  modParamMap["WA1"] = &modelPar_wa1;
  
  fadType modelPar_wa2=mod.wa2; 
  bool modelPar_given_WA2=mod.given("WA2"); 
  modParamMap["WA2"] = &modelPar_wa2;
  
  fadType modelPar_wketa=mod.wketa; 
  bool modelPar_given_WKETA=mod.given("WKETA"); 
  modParamMap["WKETA"] = &modelPar_wketa;
  
  fadType modelPar_wnsub=mod.wnsub; 
  bool modelPar_given_WNSUB=mod.given("WNSUB"); 
  modParamMap["WNSUB"] = &modelPar_wnsub;
  
  fadType modelPar_wnpeak=mod.wnpeak; 
  bool modelPar_given_WNCH=mod.given("WNCH"); 
  modParamMap["WNCH"] = &modelPar_wnpeak;
  
  fadType modelPar_wngate=mod.wngate; 
  bool modelPar_given_WNGATE=mod.given("WNGATE"); 
  modParamMap["WNGATE"] = &modelPar_wngate;
  
  fadType modelPar_wgamma1=mod.wgamma1; 
  bool modelPar_given_WGAMMA1=mod.given("WGAMMA1"); 
  modParamMap["WGAMMA1"] = &modelPar_wgamma1;
  
  fadType modelPar_wgamma2=mod.wgamma2; 
  bool modelPar_given_WGAMMA2=mod.given("WGAMMA2"); 
  modParamMap["WGAMMA2"] = &modelPar_wgamma2;
  
  fadType modelPar_wvbx=mod.wvbx; 
  bool modelPar_given_WVBX=mod.given("WVBX"); 
  modParamMap["WVBX"] = &modelPar_wvbx;
  
  fadType modelPar_wvbm=mod.wvbm; 
  bool modelPar_given_WVBM=mod.given("WVBM"); 
  modParamMap["WVBM"] = &modelPar_wvbm;
  
  fadType modelPar_wxt=mod.wxt; 
  bool modelPar_given_WXT=mod.given("WXT"); 
  modParamMap["WXT"] = &modelPar_wxt;
  
  fadType modelPar_wk1=mod.wk1; 
  bool modelPar_given_WK1=mod.given("WK1"); 
  modParamMap["WK1"] = &modelPar_wk1;
  
  fadType modelPar_wkt1=mod.wkt1; 
  bool modelPar_given_WKT1=mod.given("WKT1"); 
  modParamMap["WKT1"] = &modelPar_wkt1;
  
  fadType modelPar_wkt1l=mod.wkt1l; 
  bool modelPar_given_WKT1L=mod.given("WKT1L"); 
  modParamMap["WKT1L"] = &modelPar_wkt1l;
  
  fadType modelPar_wkt2=mod.wkt2; 
  bool modelPar_given_WKT2=mod.given("WKT2"); 
  modParamMap["WKT2"] = &modelPar_wkt2;
  
  fadType modelPar_wk2=mod.wk2; 
  bool modelPar_given_WK2=mod.given("WK2"); 
  modParamMap["WK2"] = &modelPar_wk2;
  
  fadType modelPar_wk3=mod.wk3; 
  bool modelPar_given_WK3=mod.given("WK3"); 
  modParamMap["WK3"] = &modelPar_wk3;
  
  fadType modelPar_wk3b=mod.wk3b; 
  bool modelPar_given_WK3B=mod.given("WK3B"); 
  modParamMap["WK3B"] = &modelPar_wk3b;
  
  fadType modelPar_ww0=mod.ww0; 
  bool modelPar_given_WW0=mod.given("WW0"); 
  modParamMap["WW0"] = &modelPar_ww0;
  
  fadType modelPar_wnlx=mod.wnlx; 
  bool modelPar_given_WNLX=mod.given("WNLX"); 
  modParamMap["WNLX"] = &modelPar_wnlx;
  
  fadType modelPar_wdvt0=mod.wdvt0; 
  bool modelPar_given_WDVT0=mod.given("WDVT0"); 
  modParamMap["WDVT0"] = &modelPar_wdvt0;
  
  fadType modelPar_wdvt1=mod.wdvt1; 
  bool modelPar_given_WDVT1=mod.given("WDVT1"); 
  modParamMap["WDVT1"] = &modelPar_wdvt1;
  
  fadType modelPar_wdvt2=mod.wdvt2; 
  bool modelPar_given_WDVT2=mod.given("WDVT2"); 
  modParamMap["WDVT2"] = &modelPar_wdvt2;
  
  fadType modelPar_wdvt0w=mod.wdvt0w; 
  bool modelPar_given_WDVT0W=mod.given("WDVT0W"); 
  modParamMap["WDVT0W"] = &modelPar_wdvt0w;
  
  fadType modelPar_wdvt1w=mod.wdvt1w; 
  bool modelPar_given_WDVT1W=mod.given("WDVT1W"); 
  modParamMap["WDVT1W"] = &modelPar_wdvt1w;
  
  fadType modelPar_wdvt2w=mod.wdvt2w; 
  bool modelPar_given_WDVT2W=mod.given("WDVT2W"); 
  modParamMap["WDVT2W"] = &modelPar_wdvt2w;
  
  fadType modelPar_wdrout=mod.wdrout; 
  bool modelPar_given_WDROUT=mod.given("WDROUT"); 
  modParamMap["WDROUT"] = &modelPar_wdrout;
  
  fadType modelPar_wdsub=mod.wdsub; 
  bool modelPar_given_WDSUB=mod.given("WDSUB"); 
  modParamMap["WDSUB"] = &modelPar_wdsub;
  
  fadType modelPar_wvth0=mod.wvth0; 
  bool modelPar_given_WVTH0=mod.given("WVTH0"); 
  modParamMap["WVTH0"] = &modelPar_wvth0;
  
  fadType modelPar_wua=mod.wua; 
  bool modelPar_given_WUA=mod.given("WUA"); 
  modParamMap["WUA"] = &modelPar_wua;
  
  fadType modelPar_wua1=mod.wua1; 
  bool modelPar_given_WUA1=mod.given("WUA1"); 
  modParamMap["WUA1"] = &modelPar_wua1;
  
  fadType modelPar_wub=mod.wub; 
  bool modelPar_given_WUB=mod.given("WUB"); 
  modParamMap["WUB"] = &modelPar_wub;
  
  fadType modelPar_wub1=mod.wub1; 
  bool modelPar_given_WUB1=mod.given("WUB1"); 
  modParamMap["WUB1"] = &modelPar_wub1;
  
  fadType modelPar_wuc=mod.wuc; 
  bool modelPar_given_WUC=mod.given("WUC"); 
  modParamMap["WUC"] = &modelPar_wuc;
  
  fadType modelPar_wuc1=mod.wuc1; 
  bool modelPar_given_WUC1=mod.given("WUC1"); 
  modParamMap["WUC1"] = &modelPar_wuc1;
  
  fadType modelPar_wu0=mod.wu0; 
  bool modelPar_given_WU0=mod.given("WU0"); 
  modParamMap["WU0"] = &modelPar_wu0;
  
  fadType modelPar_wute=mod.wute; 
  bool modelPar_given_WUTE=mod.given("WUTE"); 
  modParamMap["WUTE"] = &modelPar_wute;
  
  fadType modelPar_wvoff=mod.wvoff; 
  bool modelPar_given_WVOFF=mod.given("WVOFF"); 
  modParamMap["WVOFF"] = &modelPar_wvoff;
  
  fadType modelPar_wrdsw=mod.wrdsw; 
  bool modelPar_given_WRDSW=mod.given("WRDSW"); 
  modParamMap["WRDSW"] = &modelPar_wrdsw;
  
  fadType modelPar_wprwg=mod.wprwg; 
  bool modelPar_given_WPRWG=mod.given("WPRWG"); 
  modParamMap["WPRWG"] = &modelPar_wprwg;
  
  fadType modelPar_wprwb=mod.wprwb; 
  bool modelPar_given_WPRWB=mod.given("WPRWB"); 
  modParamMap["WPRWB"] = &modelPar_wprwb;
  
  fadType modelPar_wprt=mod.wprt; 
  bool modelPar_given_WPRT=mod.given("WPRT"); 
  modParamMap["WPRT"] = &modelPar_wprt;
  
  fadType modelPar_weta0=mod.weta0; 
  bool modelPar_given_WETA0=mod.given("WETA0"); 
  modParamMap["WETA0"] = &modelPar_weta0;
  
  fadType modelPar_wetab=mod.wetab; 
  bool modelPar_given_WETAB=mod.given("WETAB"); 
  modParamMap["WETAB"] = &modelPar_wetab;
  
  fadType modelPar_wpclm=mod.wpclm; 
  bool modelPar_given_WPCLM=mod.given("WPCLM"); 
  modParamMap["WPCLM"] = &modelPar_wpclm;
  
  fadType modelPar_wpdibl1=mod.wpdibl1; 
  bool modelPar_given_WPDIBLC1=mod.given("WPDIBLC1"); 
  modParamMap["WPDIBLC1"] = &modelPar_wpdibl1;
  
  fadType modelPar_wpdibl2=mod.wpdibl2; 
  bool modelPar_given_WPDIBLC2=mod.given("WPDIBLC2"); 
  modParamMap["WPDIBLC2"] = &modelPar_wpdibl2;
  
  fadType modelPar_wpdiblb=mod.wpdiblb; 
  bool modelPar_given_WPDIBLCB=mod.given("WPDIBLCB"); 
  modParamMap["WPDIBLCB"] = &modelPar_wpdiblb;
  
  fadType modelPar_wpscbe1=mod.wpscbe1; 
  bool modelPar_given_WPSCBE1=mod.given("WPSCBE1"); 
  modParamMap["WPSCBE1"] = &modelPar_wpscbe1;
  
  fadType modelPar_wpscbe2=mod.wpscbe2; 
  bool modelPar_given_WPSCBE2=mod.given("WPSCBE2"); 
  modParamMap["WPSCBE2"] = &modelPar_wpscbe2;
  
  fadType modelPar_wpvag=mod.wpvag; 
  bool modelPar_given_WPVAG=mod.given("WPVAG"); 
  modParamMap["WPVAG"] = &modelPar_wpvag;
  
  fadType modelPar_wdelta=mod.wdelta; 
  bool modelPar_given_WDELTA=mod.given("WDELTA"); 
  modParamMap["WDELTA"] = &modelPar_wdelta;
  
  fadType modelPar_wwr=mod.wwr; 
  bool modelPar_given_WWR=mod.given("WWR"); 
  modParamMap["WWR"] = &modelPar_wwr;
  
  fadType modelPar_wdwg=mod.wdwg; 
  bool modelPar_given_WDWG=mod.given("WDWG"); 
  modParamMap["WDWG"] = &modelPar_wdwg;
  
  fadType modelPar_wdwb=mod.wdwb; 
  bool modelPar_given_WDWB=mod.given("WDWB"); 
  modParamMap["WDWB"] = &modelPar_wdwb;
  
  fadType modelPar_wb0=mod.wb0; 
  bool modelPar_given_WB0=mod.given("WB0"); 
  modParamMap["WB0"] = &modelPar_wb0;
  
  fadType modelPar_wb1=mod.wb1; 
  bool modelPar_given_WB1=mod.given("WB1"); 
  modParamMap["WB1"] = &modelPar_wb1;
  
  fadType modelPar_walpha0=mod.walpha0; 
  bool modelPar_given_WALPHA0=mod.given("WALPHA0"); 
  modParamMap["WALPHA0"] = &modelPar_walpha0;
  
  fadType modelPar_walpha1=mod.walpha1; 
  bool modelPar_given_WALPHA1=mod.given("WALPHA1"); 
  modParamMap["WALPHA1"] = &modelPar_walpha1;
  
  fadType modelPar_wbeta0=mod.wbeta0; 
  bool modelPar_given_WBETA0=mod.given("WBETA0"); 
  modParamMap["WBETA0"] = &modelPar_wbeta0;
  
  fadType modelPar_wvfb=mod.wvfb; 
  bool modelPar_given_WVFB=mod.given("WVFB"); 
  modParamMap["WVFB"] = &modelPar_wvfb;
  
  fadType modelPar_welm=mod.welm; 
  bool modelPar_given_WELM=mod.given("WELM"); 
  modParamMap["WELM"] = &modelPar_welm;
  
  fadType modelPar_wcgsl=mod.wcgsl; 
  bool modelPar_given_WCGSL=mod.given("WCGSL"); 
  modParamMap["WCGSL"] = &modelPar_wcgsl;
  
  fadType modelPar_wcgdl=mod.wcgdl; 
  bool modelPar_given_WCGDL=mod.given("WCGDL"); 
  modParamMap["WCGDL"] = &modelPar_wcgdl;
  
  fadType modelPar_wckappa=mod.wckappa; 
  bool modelPar_given_WCKAPPA=mod.given("WCKAPPA"); 
  modParamMap["WCKAPPA"] = &modelPar_wckappa;
  
  fadType modelPar_wcf=mod.wcf; 
  bool modelPar_given_WCF=mod.given("WCF"); 
  modParamMap["WCF"] = &modelPar_wcf;
  
  fadType modelPar_wclc=mod.wclc; 
  bool modelPar_given_WCLC=mod.given("WCLC"); 
  modParamMap["WCLC"] = &modelPar_wclc;
  
  fadType modelPar_wcle=mod.wcle; 
  bool modelPar_given_WCLE=mod.given("WCLE"); 
  modParamMap["WCLE"] = &modelPar_wcle;
  
  fadType modelPar_wvfbcv=mod.wvfbcv; 
  bool modelPar_given_WVFBCV=mod.given("WVFBCV"); 
  modParamMap["WVFBCV"] = &modelPar_wvfbcv;
  
  fadType modelPar_wnoff=mod.wnoff; 
  bool modelPar_given_WNOFF=mod.given("WNOFF"); 
  modParamMap["WNOFF"] = &modelPar_wnoff;
  
  fadType modelPar_wvoffcv=mod.wvoffcv; 
  bool modelPar_given_WVOFFCV=mod.given("WVOFFCV"); 
  modParamMap["WVOFFCV"] = &modelPar_wvoffcv;
  
  fadType modelPar_wacde=mod.wacde; 
  bool modelPar_given_WACDE=mod.given("WACDE"); 
  modParamMap["WACDE"] = &modelPar_wacde;
  
  fadType modelPar_wmoin=mod.wmoin; 
  bool modelPar_given_WMOIN=mod.given("WMOIN"); 
  modParamMap["WMOIN"] = &modelPar_wmoin;
  
  fadType modelPar_pcdsc=mod.pcdsc; 
  bool modelPar_given_PCDSC=mod.given("PCDSC"); 
  modParamMap["PCDSC"] = &modelPar_pcdsc;
  
  fadType modelPar_pcdscb=mod.pcdscb; 
  bool modelPar_given_PCDSCB=mod.given("PCDSCB"); 
  modParamMap["PCDSCB"] = &modelPar_pcdscb;
  
  fadType modelPar_pcdscd=mod.pcdscd; 
  bool modelPar_given_PCDSCD=mod.given("PCDSCD"); 
  modParamMap["PCDSCD"] = &modelPar_pcdscd;
  
  fadType modelPar_pcit=mod.pcit; 
  bool modelPar_given_PCIT=mod.given("PCIT"); 
  modParamMap["PCIT"] = &modelPar_pcit;
  
  fadType modelPar_pnfactor=mod.pnfactor; 
  bool modelPar_given_PNFACTOR=mod.given("PNFACTOR"); 
  modParamMap["PNFACTOR"] = &modelPar_pnfactor;
  
  fadType modelPar_pxj=mod.pxj; 
  bool modelPar_given_PXJ=mod.given("PXJ"); 
  modParamMap["PXJ"] = &modelPar_pxj;
  
  fadType modelPar_pvsat=mod.pvsat; 
  bool modelPar_given_PVSAT=mod.given("PVSAT"); 
  modParamMap["PVSAT"] = &modelPar_pvsat;
  
  fadType modelPar_pat=mod.pat; 
  bool modelPar_given_PAT=mod.given("PAT"); 
  modParamMap["PAT"] = &modelPar_pat;
  
  fadType modelPar_pa0=mod.pa0; 
  bool modelPar_given_PA0=mod.given("PA0"); 
  modParamMap["PA0"] = &modelPar_pa0;
  
  fadType modelPar_pags=mod.pags; 
  bool modelPar_given_PAGS=mod.given("PAGS"); 
  modParamMap["PAGS"] = &modelPar_pags;
  
  fadType modelPar_pa1=mod.pa1; 
  bool modelPar_given_PA1=mod.given("PA1"); 
  modParamMap["PA1"] = &modelPar_pa1;
  
  fadType modelPar_pa2=mod.pa2; 
  bool modelPar_given_PA2=mod.given("PA2"); 
  modParamMap["PA2"] = &modelPar_pa2;
  
  fadType modelPar_pketa=mod.pketa; 
  bool modelPar_given_PKETA=mod.given("PKETA"); 
  modParamMap["PKETA"] = &modelPar_pketa;
  
  fadType modelPar_pnsub=mod.pnsub; 
  bool modelPar_given_PNSUB=mod.given("PNSUB"); 
  modParamMap["PNSUB"] = &modelPar_pnsub;
  
  fadType modelPar_pnpeak=mod.pnpeak; 
  bool modelPar_given_PNCH=mod.given("PNCH"); 
  modParamMap["PNCH"] = &modelPar_pnpeak;
  
  fadType modelPar_pngate=mod.pngate; 
  bool modelPar_given_PNGATE=mod.given("PNGATE"); 
  modParamMap["PNGATE"] = &modelPar_pngate;
  
  fadType modelPar_pgamma1=mod.pgamma1; 
  bool modelPar_given_PGAMMA1=mod.given("PGAMMA1"); 
  modParamMap["PGAMMA1"] = &modelPar_pgamma1;
  
  fadType modelPar_pgamma2=mod.pgamma2; 
  bool modelPar_given_PGAMMA2=mod.given("PGAMMA2"); 
  modParamMap["PGAMMA2"] = &modelPar_pgamma2;
  
  fadType modelPar_pvbx=mod.pvbx; 
  bool modelPar_given_PVBX=mod.given("PVBX"); 
  modParamMap["PVBX"] = &modelPar_pvbx;
  
  fadType modelPar_pvbm=mod.pvbm; 
  bool modelPar_given_PVBM=mod.given("PVBM"); 
  modParamMap["PVBM"] = &modelPar_pvbm;
  
  fadType modelPar_pxt=mod.pxt; 
  bool modelPar_given_PXT=mod.given("PXT"); 
  modParamMap["PXT"] = &modelPar_pxt;
  
  fadType modelPar_pk1=mod.pk1; 
  bool modelPar_given_PK1=mod.given("PK1"); 
  modParamMap["PK1"] = &modelPar_pk1;
  
  fadType modelPar_pkt1=mod.pkt1; 
  bool modelPar_given_PKT1=mod.given("PKT1"); 
  modParamMap["PKT1"] = &modelPar_pkt1;
  
  fadType modelPar_pkt1l=mod.pkt1l; 
  bool modelPar_given_PKT1L=mod.given("PKT1L"); 
  modParamMap["PKT1L"] = &modelPar_pkt1l;
  
  fadType modelPar_pkt2=mod.pkt2; 
  bool modelPar_given_PKT2=mod.given("PKT2"); 
  modParamMap["PKT2"] = &modelPar_pkt2;
  
  fadType modelPar_pk2=mod.pk2; 
  bool modelPar_given_PK2=mod.given("PK2"); 
  modParamMap["PK2"] = &modelPar_pk2;
  
  fadType modelPar_pk3=mod.pk3; 
  bool modelPar_given_PK3=mod.given("PK3"); 
  modParamMap["PK3"] = &modelPar_pk3;
  
  fadType modelPar_pk3b=mod.pk3b; 
  bool modelPar_given_PK3B=mod.given("PK3B"); 
  modParamMap["PK3B"] = &modelPar_pk3b;
  
  fadType modelPar_pw0=mod.pw0; 
  bool modelPar_given_PW0=mod.given("PW0"); 
  modParamMap["PW0"] = &modelPar_pw0;
  
  fadType modelPar_pnlx=mod.pnlx; 
  bool modelPar_given_PNLX=mod.given("PNLX"); 
  modParamMap["PNLX"] = &modelPar_pnlx;
  
  fadType modelPar_pdvt0=mod.pdvt0; 
  bool modelPar_given_PDVT0=mod.given("PDVT0"); 
  modParamMap["PDVT0"] = &modelPar_pdvt0;
  
  fadType modelPar_pdvt1=mod.pdvt1; 
  bool modelPar_given_PDVT1=mod.given("PDVT1"); 
  modParamMap["PDVT1"] = &modelPar_pdvt1;
  
  fadType modelPar_pdvt2=mod.pdvt2; 
  bool modelPar_given_PDVT2=mod.given("PDVT2"); 
  modParamMap["PDVT2"] = &modelPar_pdvt2;
  
  fadType modelPar_pdvt0w=mod.pdvt0w; 
  bool modelPar_given_PDVT0W=mod.given("PDVT0W"); 
  modParamMap["PDVT0W"] = &modelPar_pdvt0w;
  
  fadType modelPar_pdvt1w=mod.pdvt1w; 
  bool modelPar_given_PDVT1W=mod.given("PDVT1W"); 
  modParamMap["PDVT1W"] = &modelPar_pdvt1w;
  
  fadType modelPar_pdvt2w=mod.pdvt2w; 
  bool modelPar_given_PDVT2W=mod.given("PDVT2W"); 
  modParamMap["PDVT2W"] = &modelPar_pdvt2w;
  
  fadType modelPar_pdrout=mod.pdrout; 
  bool modelPar_given_PDROUT=mod.given("PDROUT"); 
  modParamMap["PDROUT"] = &modelPar_pdrout;
  
  fadType modelPar_pdsub=mod.pdsub; 
  bool modelPar_given_PDSUB=mod.given("PDSUB"); 
  modParamMap["PDSUB"] = &modelPar_pdsub;
  
  fadType modelPar_pvth0=mod.pvth0; 
  bool modelPar_given_PVTH0=mod.given("PVTH0"); 
  modParamMap["PVTH0"] = &modelPar_pvth0;
  
  fadType modelPar_pua=mod.pua; 
  bool modelPar_given_PUA=mod.given("PUA"); 
  modParamMap["PUA"] = &modelPar_pua;
  
  fadType modelPar_pua1=mod.pua1; 
  bool modelPar_given_PUA1=mod.given("PUA1"); 
  modParamMap["PUA1"] = &modelPar_pua1;
  
  fadType modelPar_pub=mod.pub; 
  bool modelPar_given_PUB=mod.given("PUB"); 
  modParamMap["PUB"] = &modelPar_pub;
  
  fadType modelPar_pub1=mod.pub1; 
  bool modelPar_given_PUB1=mod.given("PUB1"); 
  modParamMap["PUB1"] = &modelPar_pub1;
  
  fadType modelPar_puc=mod.puc; 
  bool modelPar_given_PUC=mod.given("PUC"); 
  modParamMap["PUC"] = &modelPar_puc;
  
  fadType modelPar_puc1=mod.puc1; 
  bool modelPar_given_PUC1=mod.given("PUC1"); 
  modParamMap["PUC1"] = &modelPar_puc1;
  
  fadType modelPar_pu0=mod.pu0; 
  bool modelPar_given_PU0=mod.given("PU0"); 
  modParamMap["PU0"] = &modelPar_pu0;
  
  fadType modelPar_pute=mod.pute; 
  bool modelPar_given_PUTE=mod.given("PUTE"); 
  modParamMap["PUTE"] = &modelPar_pute;
  
  fadType modelPar_pvoff=mod.pvoff; 
  bool modelPar_given_PVOFF=mod.given("PVOFF"); 
  modParamMap["PVOFF"] = &modelPar_pvoff;
  
  fadType modelPar_prdsw=mod.prdsw; 
  bool modelPar_given_PRDSW=mod.given("PRDSW"); 
  modParamMap["PRDSW"] = &modelPar_prdsw;
  
  fadType modelPar_pprwg=mod.pprwg; 
  bool modelPar_given_PPRWG=mod.given("PPRWG"); 
  modParamMap["PPRWG"] = &modelPar_pprwg;
  
  fadType modelPar_pprwb=mod.pprwb; 
  bool modelPar_given_PPRWB=mod.given("PPRWB"); 
  modParamMap["PPRWB"] = &modelPar_pprwb;
  
  fadType modelPar_pprt=mod.pprt; 
  bool modelPar_given_PPRT=mod.given("PPRT"); 
  modParamMap["PPRT"] = &modelPar_pprt;
  
  fadType modelPar_peta0=mod.peta0; 
  bool modelPar_given_PETA0=mod.given("PETA0"); 
  modParamMap["PETA0"] = &modelPar_peta0;
  
  fadType modelPar_petab=mod.petab; 
  bool modelPar_given_PETAB=mod.given("PETAB"); 
  modParamMap["PETAB"] = &modelPar_petab;
  
  fadType modelPar_ppclm=mod.ppclm; 
  bool modelPar_given_PPCLM=mod.given("PPCLM"); 
  modParamMap["PPCLM"] = &modelPar_ppclm;
  
  fadType modelPar_ppdibl1=mod.ppdibl1; 
  bool modelPar_given_PPDIBLC1=mod.given("PPDIBLC1"); 
  modParamMap["PPDIBLC1"] = &modelPar_ppdibl1;
  
  fadType modelPar_ppdibl2=mod.ppdibl2; 
  bool modelPar_given_PPDIBLC2=mod.given("PPDIBLC2"); 
  modParamMap["PPDIBLC2"] = &modelPar_ppdibl2;
  
  fadType modelPar_ppdiblb=mod.ppdiblb; 
  bool modelPar_given_PPDIBLCB=mod.given("PPDIBLCB"); 
  modParamMap["PPDIBLCB"] = &modelPar_ppdiblb;
  
  fadType modelPar_ppscbe1=mod.ppscbe1; 
  bool modelPar_given_PPSCBE1=mod.given("PPSCBE1"); 
  modParamMap["PPSCBE1"] = &modelPar_ppscbe1;
  
  fadType modelPar_ppscbe2=mod.ppscbe2; 
  bool modelPar_given_PPSCBE2=mod.given("PPSCBE2"); 
  modParamMap["PPSCBE2"] = &modelPar_ppscbe2;
  
  fadType modelPar_ppvag=mod.ppvag; 
  bool modelPar_given_PPVAG=mod.given("PPVAG"); 
  modParamMap["PPVAG"] = &modelPar_ppvag;
  
  fadType modelPar_pdelta=mod.pdelta; 
  bool modelPar_given_PDELTA=mod.given("PDELTA"); 
  modParamMap["PDELTA"] = &modelPar_pdelta;
  
  fadType modelPar_pwr=mod.pwr; 
  bool modelPar_given_PWR=mod.given("PWR"); 
  modParamMap["PWR"] = &modelPar_pwr;
  
  fadType modelPar_pdwg=mod.pdwg; 
  bool modelPar_given_PDWG=mod.given("PDWG"); 
  modParamMap["PDWG"] = &modelPar_pdwg;
  
  fadType modelPar_pdwb=mod.pdwb; 
  bool modelPar_given_PDWB=mod.given("PDWB"); 
  modParamMap["PDWB"] = &modelPar_pdwb;
  
  fadType modelPar_pb0=mod.pb0; 
  bool modelPar_given_PB0=mod.given("PB0"); 
  modParamMap["PB0"] = &modelPar_pb0;
  
  fadType modelPar_pb1=mod.pb1; 
  bool modelPar_given_PB1=mod.given("PB1"); 
  modParamMap["PB1"] = &modelPar_pb1;
  
  fadType modelPar_palpha0=mod.palpha0; 
  bool modelPar_given_PALPHA0=mod.given("PALPHA0"); 
  modParamMap["PALPHA0"] = &modelPar_palpha0;
  
  fadType modelPar_palpha1=mod.palpha1; 
  bool modelPar_given_PALPHA1=mod.given("PALPHA1"); 
  modParamMap["PALPHA1"] = &modelPar_palpha1;
  
  fadType modelPar_pbeta0=mod.pbeta0; 
  bool modelPar_given_PBETA0=mod.given("PBETA0"); 
  modParamMap["PBETA0"] = &modelPar_pbeta0;
  
  fadType modelPar_pvfb=mod.pvfb; 
  bool modelPar_given_PVFB=mod.given("PVFB"); 
  modParamMap["PVFB"] = &modelPar_pvfb;
  
  fadType modelPar_pelm=mod.pelm; 
  bool modelPar_given_PELM=mod.given("PELM"); 
  modParamMap["PELM"] = &modelPar_pelm;
  
  fadType modelPar_pcgsl=mod.pcgsl; 
  bool modelPar_given_PCGSL=mod.given("PCGSL"); 
  modParamMap["PCGSL"] = &modelPar_pcgsl;
  
  fadType modelPar_pcgdl=mod.pcgdl; 
  bool modelPar_given_PCGDL=mod.given("PCGDL"); 
  modParamMap["PCGDL"] = &modelPar_pcgdl;
  
  fadType modelPar_pckappa=mod.pckappa; 
  bool modelPar_given_PCKAPPA=mod.given("PCKAPPA"); 
  modParamMap["PCKAPPA"] = &modelPar_pckappa;
  
  fadType modelPar_pcf=mod.pcf; 
  bool modelPar_given_PCF=mod.given("PCF"); 
  modParamMap["PCF"] = &modelPar_pcf;
  
  fadType modelPar_pclc=mod.pclc; 
  bool modelPar_given_PCLC=mod.given("PCLC"); 
  modParamMap["PCLC"] = &modelPar_pclc;
  
  fadType modelPar_pcle=mod.pcle; 
  bool modelPar_given_PCLE=mod.given("PCLE"); 
  modParamMap["PCLE"] = &modelPar_pcle;
  
  fadType modelPar_pvfbcv=mod.pvfbcv; 
  bool modelPar_given_PVFBCV=mod.given("PVFBCV"); 
  modParamMap["PVFBCV"] = &modelPar_pvfbcv;
  
  fadType modelPar_pnoff=mod.pnoff; 
  bool modelPar_given_PNOFF=mod.given("PNOFF"); 
  modParamMap["PNOFF"] = &modelPar_pnoff;
  
  fadType modelPar_pvoffcv=mod.pvoffcv; 
  bool modelPar_given_PVOFFCV=mod.given("PVOFFCV"); 
  modParamMap["PVOFFCV"] = &modelPar_pvoffcv;
  
  fadType modelPar_pacde=mod.pacde; 
  bool modelPar_given_PACDE=mod.given("PACDE"); 
  modParamMap["PACDE"] = &modelPar_pacde;
  
  fadType modelPar_pmoin=mod.pmoin; 
  bool modelPar_given_PMOIN=mod.given("PMOIN"); 
  modParamMap["PMOIN"] = &modelPar_pmoin;
  
  fadType modelPar_tnom=mod.tnom; 
  bool modelPar_given_TNOM=mod.given("TNOM"); 
  modParamMap["TNOM"] = &modelPar_tnom;
  
  fadType modelPar_cgso=mod.cgso; 
  bool modelPar_given_CGSO=mod.given("CGSO"); 
  modParamMap["CGSO"] = &modelPar_cgso;
  
  fadType modelPar_cgdo=mod.cgdo; 
  bool modelPar_given_CGDO=mod.given("CGDO"); 
  modParamMap["CGDO"] = &modelPar_cgdo;
  
  fadType modelPar_cgbo=mod.cgbo; 
  bool modelPar_given_CGBO=mod.given("CGBO"); 
  modParamMap["CGBO"] = &modelPar_cgbo;
  
  fadType modelPar_xpart=mod.xpart; 
  bool modelPar_given_XPART=mod.given("XPART"); 
  modParamMap["XPART"] = &modelPar_xpart;
  
  fadType modelPar_sheetResistance=mod.sheetResistance; 
  bool modelPar_given_RSH=mod.given("RSH"); 
  modParamMap["RSH"] = &modelPar_sheetResistance;
  
  fadType modelPar_jctSatCurDensity=mod.jctSatCurDensity; 
  bool modelPar_given_JS=mod.given("JS"); 
  modParamMap["JS"] = &modelPar_jctSatCurDensity;
  
  fadType modelPar_jctSidewallSatCurDensity=mod.jctSidewallSatCurDensity; 
  bool modelPar_given_JSW=mod.given("JSW"); 
  modParamMap["JSW"] = &modelPar_jctSidewallSatCurDensity;
  
  fadType modelPar_bulkJctPotential=mod.bulkJctPotential; 
  bool modelPar_given_PB=mod.given("PB"); 
  modParamMap["PB"] = &modelPar_bulkJctPotential;
  
  fadType modelPar_bulkJctBotGradingCoeff=mod.bulkJctBotGradingCoeff; 
  bool modelPar_given_MJ=mod.given("MJ"); 
  modParamMap["MJ"] = &modelPar_bulkJctBotGradingCoeff;
  
  fadType modelPar_sidewallJctPotential=mod.sidewallJctPotential; 
  bool modelPar_given_PBSW=mod.given("PBSW"); 
  modParamMap["PBSW"] = &modelPar_sidewallJctPotential;
  
  fadType modelPar_GatesidewallJctPotential=mod.GatesidewallJctPotential; 
  bool modelPar_given_PBSWG=mod.given("PBSWG"); 
  modParamMap["PBSWG"] = &modelPar_GatesidewallJctPotential;
  
  fadType modelPar_bulkJctSideGradingCoeff=mod.bulkJctSideGradingCoeff; 
  bool modelPar_given_MJSW=mod.given("MJSW"); 
  modParamMap["MJSW"] = &modelPar_bulkJctSideGradingCoeff;
  
  fadType modelPar_unitAreaJctCap=mod.unitAreaJctCap; 
  bool modelPar_given_CJ=mod.given("CJ"); 
  modParamMap["CJ"] = &modelPar_unitAreaJctCap;
  
  fadType modelPar_unitLengthSidewallJctCap=mod.unitLengthSidewallJctCap; 
  bool modelPar_given_CJSW=mod.given("CJSW"); 
  modParamMap["CJSW"] = &modelPar_unitLengthSidewallJctCap;
  
  fadType modelPar_bulkJctGateSideGradingCoeff=mod.bulkJctGateSideGradingCoeff; 
  bool modelPar_given_MJSWG=mod.given("MJSWG"); 
  modParamMap["MJSWG"] = &modelPar_bulkJctGateSideGradingCoeff;
  
  fadType modelPar_unitLengthGateSidewallJctCap=mod.unitLengthGateSidewallJctCap; 
  bool modelPar_given_CJSWG=mod.given("CJSWG"); 
  modParamMap["CJSWG"] = &modelPar_unitLengthGateSidewallJctCap;
  
  fadType modelPar_jctEmissionCoeff=mod.jctEmissionCoeff; 
  bool modelPar_given_NJ=mod.given("NJ"); 
  modParamMap["NJ"] = &modelPar_jctEmissionCoeff;
  
  fadType modelPar_jctTempExponent=mod.jctTempExponent; 
  bool modelPar_given_XTI=mod.given("XTI"); 
  modParamMap["XTI"] = &modelPar_jctTempExponent;
  
  fadType modelPar_oxideTrapDensityA=mod.oxideTrapDensityA; 
  bool modelPar_given_NOIA=mod.given("NOIA"); 
  modParamMap["NOIA"] = &modelPar_oxideTrapDensityA;
  
  fadType modelPar_oxideTrapDensityB=mod.oxideTrapDensityB; 
  bool modelPar_given_NOIB=mod.given("NOIB"); 
  modParamMap["NOIB"] = &modelPar_oxideTrapDensityB;
  
  fadType modelPar_oxideTrapDensityC=mod.oxideTrapDensityC; 
  bool modelPar_given_NOIC=mod.given("NOIC"); 
  modParamMap["NOIC"] = &modelPar_oxideTrapDensityC;
  
  fadType modelPar_em=mod.em; 
  bool modelPar_given_EM=mod.given("EM"); 
  modParamMap["EM"] = &modelPar_em;
  
  fadType modelPar_ef=mod.ef; 
  bool modelPar_given_EF=mod.given("EF"); 
  modParamMap["EF"] = &modelPar_ef;
  
  fadType modelPar_af=mod.af; 
  bool modelPar_given_AF=mod.given("AF"); 
  modParamMap["AF"] = &modelPar_af;
  
  fadType modelPar_kf=mod.kf; 
  bool modelPar_given_KF=mod.given("KF"); 
  modParamMap["KF"] = &modelPar_kf;
  
  fadType modelPar_lintnoi=mod.lintnoi; 
  bool modelPar_given_LINTNOI=mod.given("LINTNOI"); 
  modParamMap["LINTNOI"] = &modelPar_lintnoi;
  
  fadType modelPar_Lint=mod.Lint; 
  bool modelPar_given_LINT=mod.given("LINT"); 
  modParamMap["LINT"] = &modelPar_Lint;
  
  fadType modelPar_Ll=mod.Ll; 
  bool modelPar_given_LL=mod.given("LL"); 
  modParamMap["LL"] = &modelPar_Ll;
  
  fadType modelPar_Llc=mod.Llc; 
  bool modelPar_given_LLC=mod.given("LLC"); 
  modParamMap["LLC"] = &modelPar_Llc;
  
  fadType modelPar_Lln=mod.Lln; 
  bool modelPar_given_LLN=mod.given("LLN"); 
  modParamMap["LLN"] = &modelPar_Lln;
  
  fadType modelPar_Lw=mod.Lw; 
  bool modelPar_given_LW=mod.given("LW"); 
  modParamMap["LW"] = &modelPar_Lw;
  
  fadType modelPar_Lwc=mod.Lwc; 
  bool modelPar_given_LWC=mod.given("LWC"); 
  modParamMap["LWC"] = &modelPar_Lwc;
  
  fadType modelPar_Lwn=mod.Lwn; 
  bool modelPar_given_LWN=mod.given("LWN"); 
  modParamMap["LWN"] = &modelPar_Lwn;
  
  fadType modelPar_Lwl=mod.Lwl; 
  bool modelPar_given_LWL=mod.given("LWL"); 
  modParamMap["LWL"] = &modelPar_Lwl;
  
  fadType modelPar_Lwlc=mod.Lwlc; 
  bool modelPar_given_LWLC=mod.given("LWLC"); 
  modParamMap["LWLC"] = &modelPar_Lwlc;
  
  fadType modelPar_Wint=mod.Wint; 
  bool modelPar_given_WINT=mod.given("WINT"); 
  modParamMap["WINT"] = &modelPar_Wint;
  
  fadType modelPar_Wl=mod.Wl; 
  bool modelPar_given_WL=mod.given("WL"); 
  modParamMap["WL"] = &modelPar_Wl;
  
  fadType modelPar_Wlc=mod.Wlc; 
  bool modelPar_given_WLC=mod.given("WLC"); 
  modParamMap["WLC"] = &modelPar_Wlc;
  
  fadType modelPar_Wln=mod.Wln; 
  bool modelPar_given_WLN=mod.given("WLN"); 
  modParamMap["WLN"] = &modelPar_Wln;
  
  fadType modelPar_Ww=mod.Ww; 
  bool modelPar_given_WW=mod.given("WW"); 
  modParamMap["WW"] = &modelPar_Ww;
  
  fadType modelPar_Wwc=mod.Wwc; 
  bool modelPar_given_WWC=mod.given("WWC"); 
  modParamMap["WWC"] = &modelPar_Wwc;
  
  fadType modelPar_Wwn=mod.Wwn; 
  bool modelPar_given_WWN=mod.given("WWN"); 
  modParamMap["WWN"] = &modelPar_Wwn;
  
  fadType modelPar_Wwl=mod.Wwl; 
  bool modelPar_given_WWL=mod.given("WWL"); 
  modParamMap["WWL"] = &modelPar_Wwl;
  
  fadType modelPar_Wwlc=mod.Wwlc; 
  bool modelPar_given_WWLC=mod.given("WWLC"); 
  modParamMap["WWLC"] = &modelPar_Wwlc;
  
  fadType modelPar_model_l=mod.model_l; 
  bool modelPar_given_L=mod.given("L"); 
  modParamMap["L"] = &modelPar_model_l;
  
  fadType modelPar_model_w=mod.model_w; 
  bool modelPar_given_W=mod.given("W"); 
  modParamMap["W"] = &modelPar_model_w;
  
  fadType modelPar_Lmax=mod.Lmax; 
  bool modelPar_given_LMAX=mod.given("LMAX"); 
  modParamMap["LMAX"] = &modelPar_Lmax;
  
  fadType modelPar_Lmin=mod.Lmin; 
  bool modelPar_given_LMIN=mod.given("LMIN"); 
  modParamMap["LMIN"] = &modelPar_Lmin;
  
  fadType modelPar_Wmax=mod.Wmax; 
  bool modelPar_given_WMAX=mod.given("WMAX"); 
  modParamMap["WMAX"] = &modelPar_Wmax;
  
  fadType modelPar_Wmin=mod.Wmin; 
  bool modelPar_given_WMIN=mod.given("WMIN"); 
  modParamMap["WMIN"] = &modelPar_Wmin;

  // model params NOT sensitivity supported:
  int  dtype = mod.dtype;

  int mobMod = mod.mobMod; // int
  int binUnit = mod.binUnit; // int
  int capMod = mod.capMod; // int
  int paramChk = mod.paramChk; // int
  int noiMod = mod.noiMod; // int
  std::string version = mod.version; // string

  modParamMap[name]->diff(0,1);

  // Now loop over all instances and do the deed
  int inst=0;
  for (std::vector<Instance*>::const_iterator in_it=mod.instanceContainer.begin(); in_it != mod.instanceContainer.end(); ++in_it,++inst)
  {
    Instance & in=*(*in_it);

  // instance params
  fadType instancePar_temp=in.temp;	  
  bool instancePar_given_temp=in.given("TEMP");

  fadType instancePar_dtemp=in.dtemp;	  
  bool instancePar_given_dtemp=in.given("DTEMP");
  
  fadType instancePar_l=in.l;	  
  bool instancePar_given_l=in.given("L");
  
  fadType instancePar_w=in.w;	  
  bool instancePar_given_w=in.given("W");
  
  fadType instancePar_drainArea=in.drainArea;	  
  bool instancePar_given_drainArea=in.given("AD");
  
  fadType instancePar_sourceArea=in.sourceArea;	  
  bool instancePar_given_sourceArea=in.given("AS");
  
  fadType instancePar_drainSquares=in.drainSquares;	  
  bool instancePar_given_drainSquares=in.given("NRD");
  
  fadType instancePar_sourceSquares=in.sourceSquares;	  
  bool instancePar_given_sourceSquares=in.given("NRS");
  
  fadType instancePar_drainPerimeter=in.drainPerimeter;	  
  bool instancePar_given_drainPerimeter=in.given("PD");
  
  fadType instancePar_sourcePerimeter=in.sourcePerimeter;	  
  bool instancePar_given_sourcePerimeter=in.given("PS");

  bool instancePar_given_nqsMod = in.given("NQSMOD");

  // instance params NOT sensitivity supported:  
  double numberParallel = in.numberParallel ;
  double icVDS = in.icVDS ;
  double icVGS = in.icVGS ;
  double icVBS = in.icVBS ;
  int nqsMod = in.nqsMod ;
  bool OFF = in.OFF ;

  // other variables
  fadType drainConductance, sourceConductance;

  fadType modelPar_cox = mod.cox;
  fadType modelPar_vcrit = mod.vcrit;
  fadType modelPar_factor1 = mod.factor1;
  fadType modelPar_Vtm0 = mod.Vtm0;
  fadType modelPar_Eg0 = mod.Eg0;
  fadType modelPar_ni  = mod.ni;

  // model process params:
  processModelParams (
      // givens:
      modelPar_given_TOXM, modelPar_given_DSUB, modelPar_given_LLC, modelPar_given_LWC,
      modelPar_given_LWLC, modelPar_given_WLC, modelPar_given_WWL, modelPar_given_WWLC,
      modelPar_given_DWC, modelPar_given_DLC, modelPar_given_CF, modelPar_given_CGDO,
      modelPar_given_CGSO, modelPar_given_CGBO, modelPar_given_CJSWG, modelPar_given_PBSWG,
      modelPar_given_MJSWG, 
      // inputs:
      modelPar_tox, modelPar_drout, modelPar_Ll, modelPar_Lw,
      modelPar_Lwl, modelPar_Wl, modelPar_Wwl, modelPar_Wint,
      modelPar_Lint, modelPar_tnom, 
      modelPar_cgdl, modelPar_cgsl,
      modelPar_unitLengthSidewallJctCap, modelPar_bulkJctSideGradingCoeff, modelPar_xj, 
      //outputs:
      modelPar_toxm, modelPar_cox, modelPar_dsub, modelPar_Llc, modelPar_Lwc,
      modelPar_Lwlc, modelPar_Wlc, modelPar_Wwlc, modelPar_dwc, modelPar_dlc,
      modelPar_cf, modelPar_cgdo, modelPar_cgso, modelPar_cgbo,
      modelPar_unitLengthGateSidewallJctCap, modelPar_GatesidewallJctPotential, modelPar_bulkJctGateSideGradingCoeff,
      modelPar_bulkJctPotential, modelPar_sidewallJctPotential,
      modelPar_vcrit, modelPar_factor1, modelPar_Vtm0, modelPar_Eg0, modelPar_ni 
      );

  // instance variables:
  processParams (
    in.getDeviceOptions(),

    instancePar_given_temp,
    instancePar_given_dtemp,
    instancePar_given_l,
    instancePar_given_w,
  instancePar_given_drainArea,
  instancePar_given_sourceArea,
  instancePar_given_nqsMod,

   // inputs
    modelPar_model_l,
    modelPar_model_w,

    modelPar_sheetResistance,

    instancePar_drainSquares,
    instancePar_sourceSquares,

   // outputs
    instancePar_temp,
    instancePar_dtemp,
    instancePar_l,
    instancePar_w,
    instancePar_drainArea,
    instancePar_sourceArea,
      drainConductance,
      sourceConductance,
     nqsMod
    );

  fadType instance_PhiBSWTemp = in.PhiBSWTemp;
  fadType instance_PhiBSWGTemp = in.PhiBSWGTemp;

  fadType instance_vtm = in.vtm;
  fadType model_vtm = mod.vtm;

  fadType instance_jctTempSatCurDensity = in.jctTempSatCurDensity;
  fadType instance_jctSidewallTempSatCurDensity = in.jctSidewallTempSatCurDensity;
  fadType instance_unitAreaJctCapTemp = in.unitAreaJctCapTemp;
  fadType instance_unitLengthSidewallJctCapTemp = in.unitLengthSidewallJctCapTemp;
  fadType instance_unitLengthGateSidewallJctCapTemp = in.unitLengthGateSidewallJctCapTemp;

  fadType instance_PhiBTemp = in.PhiBTemp;

  fadType instance_cgso = in.cgso;
  fadType instance_cgdo = in.cgdo;
  fadType instance_vjsm = in.vjsm;
  fadType instance_IsEvjsm = in.IsEvjsm;
  fadType instance_vjdm = in.vjdm;
  fadType instance_IsEvjdm = in.IsEvjdm;

  bool instance_updateTemperatureCalled_ = in.updateTemperatureCalled_;  // seems useless  ... 

  // ----------------------------------------------------------------------------
  SizeDependParam<fadType> sizeDepParams;

  updateTemperature (
    instancePar_temp,
    instancePar_dtemp,
  modelPar_Eg0,
  modelPar_ni,
  modelPar_Vtm0,
    modelPar_tnom,
    modelPar_jctTempExponent,
    modelPar_jctEmissionCoeff,
    modelPar_jctSatCurDensity,
    modelPar_jctSidewallSatCurDensity,
    modelPar_tcj,
    modelPar_unitAreaJctCap,
    modelPar_tcjsw,
    modelPar_unitLengthSidewallJctCap, 
    modelPar_tcjswg,
    modelPar_unitLengthGateSidewallJctCap,
    modelPar_bulkJctPotential,
    modelPar_tpb,
  instance_PhiBSWTemp, // non-CONST !!
    modelPar_sidewallJctPotential,
    modelPar_tpbsw,
  instance_PhiBSWGTemp, // non-CONST !!
    modelPar_GatesidewallJctPotential,
    modelPar_tpbswg,
    instancePar_l,
    instancePar_w,
    modelPar_Lln,
    modelPar_Lwn,
    modelPar_Ll,
    modelPar_Lw,
    modelPar_Lwl,
    modelPar_Lint,
    modelPar_Llc,
    modelPar_Lwc,
    modelPar_Lwlc,
    modelPar_dlc,
    modelPar_Wln,
    modelPar_Wwn,
    modelPar_Wl,
    modelPar_Ww,
    modelPar_Wwl,
    modelPar_Wint,
    modelPar_Wlc,
    modelPar_Wwc,
    modelPar_Wwlc,
    modelPar_dwc,
    modelPar_cgdo,
    modelPar_cgso,
    modelPar_cgbo,
    modelPar_cox,
    modelPar_given_NCH, //modelPar_npeakGiven,
    modelPar_given_GAMMA1, //modelPar_gamma1Given,
    modelPar_ijth,
    instancePar_sourceArea, 
    instancePar_sourcePerimeter,
    instancePar_drainArea,
    instancePar_drainPerimeter,
    modelPar_cdsc,
    modelPar_lcdsc,
    modelPar_wcdsc,
    modelPar_pcdsc,
    modelPar_cdscb,
    modelPar_lcdscb,
    modelPar_wcdscb,
    modelPar_pcdscb,
    modelPar_cdscd,
    modelPar_lcdscd,
    modelPar_wcdscd,
    modelPar_pcdscd,
    modelPar_cit,
    modelPar_lcit,
    modelPar_wcit,
    modelPar_pcit,
    modelPar_nfactor,
    modelPar_lnfactor,
    modelPar_wnfactor,
    modelPar_pnfactor,
    modelPar_xj,
    modelPar_lxj,
    modelPar_wxj,
    modelPar_pxj,
    modelPar_vsat,
    modelPar_lvsat,
    modelPar_wvsat,
    modelPar_pvsat,
    modelPar_at,
    modelPar_lat,
    modelPar_wat,
    modelPar_pat,
    modelPar_a0,
    modelPar_la0,
    modelPar_wa0,
    modelPar_pa0,
    modelPar_ags,
    modelPar_lags,
    modelPar_wags,
    modelPar_pags,
    modelPar_a1,
    modelPar_la1,
    modelPar_wa1,
    modelPar_pa1,
    modelPar_a2,
    modelPar_la2,
    modelPar_wa2,
    modelPar_pa2,
    modelPar_keta,
    modelPar_lketa,
    modelPar_wketa,
    modelPar_pketa,
    modelPar_nsub,
    modelPar_lnsub,
    modelPar_wnsub,
    modelPar_pnsub,
    modelPar_npeak,
    modelPar_lnpeak,
    modelPar_wnpeak,
    modelPar_pnpeak,
    modelPar_ngate,
    modelPar_lngate,
    modelPar_wngate,
    modelPar_pngate,
    modelPar_gamma1,
    modelPar_lgamma1,
    modelPar_wgamma1,
    modelPar_pgamma1,
    modelPar_gamma2,
    modelPar_lgamma2,
    modelPar_wgamma2,
    modelPar_pgamma2,
    modelPar_vbx,
    modelPar_lvbx,
    modelPar_wvbx,
    modelPar_pvbx,
    modelPar_vbm,
    modelPar_lvbm,
    modelPar_wvbm,
    modelPar_pvbm,
    modelPar_xt,
    modelPar_lxt,
    modelPar_wxt,
    modelPar_pxt,
    modelPar_vfb,
    modelPar_lvfb,
    modelPar_wvfb,
    modelPar_pvfb,
    modelPar_k1,
    modelPar_lk1,
    modelPar_wk1,
    modelPar_pk1,
    modelPar_kt1,
    modelPar_lkt1,
    modelPar_wkt1,
    modelPar_pkt1,
    modelPar_kt1l,
    modelPar_lkt1l,
    modelPar_wkt1l,
    modelPar_pkt1l,
    modelPar_k2,
    modelPar_lk2,
    modelPar_wk2,
    modelPar_pk2,
    modelPar_kt2,
    modelPar_lkt2,
    modelPar_wkt2,
    modelPar_pkt2,
    modelPar_k3,
    modelPar_lk3,
    modelPar_wk3,
    modelPar_pk3,
    modelPar_k3b,
    modelPar_lk3b,
    modelPar_wk3b,
    modelPar_pk3b,
    modelPar_w0,
    modelPar_lw0,
    modelPar_ww0,
    modelPar_pw0,
    modelPar_nlx,
    modelPar_lnlx,
    modelPar_wnlx,
    modelPar_pnlx,
    modelPar_dvt0,
    modelPar_ldvt0,
    modelPar_wdvt0,
    modelPar_pdvt0,
    modelPar_dvt1,
    modelPar_ldvt1,
    modelPar_wdvt1,
    modelPar_pdvt1,
    modelPar_dvt2,
    modelPar_ldvt2,
    modelPar_wdvt2,
    modelPar_pdvt2,
    modelPar_dvt0w,
    modelPar_ldvt0w,
    modelPar_wdvt0w,
    modelPar_pdvt0w,
    modelPar_dvt1w,
    modelPar_ldvt1w,
    modelPar_wdvt1w,
    modelPar_pdvt1w,
    modelPar_dvt2w,
    modelPar_ldvt2w,
    modelPar_wdvt2w,
    modelPar_pdvt2w,
    modelPar_drout,
    modelPar_ldrout,
    modelPar_wdrout,
    modelPar_pdrout,
    modelPar_dsub,
    modelPar_ldsub,
    modelPar_wdsub,
    modelPar_pdsub,
    modelPar_vth0,
    modelPar_lvth0,
    modelPar_wvth0,
    modelPar_pvth0,
    modelPar_ua,
    modelPar_lua,
    modelPar_wua,
    modelPar_pua,
    modelPar_ua1,
    modelPar_lua1,
    modelPar_wua1,
    modelPar_pua1,
    modelPar_ub,
    modelPar_lub,
    modelPar_wub,
    modelPar_pub,
    modelPar_ub1,
    modelPar_lub1,
    modelPar_wub1,
    modelPar_pub1,
    modelPar_uc,
    modelPar_luc,
    modelPar_wuc,
    modelPar_puc,
    modelPar_uc1,
    modelPar_luc1,
    modelPar_wuc1,
    modelPar_puc1,
    modelPar_u0,
    modelPar_lu0,
    modelPar_wu0,
    modelPar_pu0,
    modelPar_ute,
    modelPar_lute,
    modelPar_wute,
    modelPar_pute,
    modelPar_voff,
    modelPar_lvoff,
    modelPar_wvoff,
    modelPar_pvoff,
    modelPar_delta,
    modelPar_ldelta,
    modelPar_wdelta,
    modelPar_pdelta,
    modelPar_rdsw,
    modelPar_lrdsw,
    modelPar_wrdsw,
    modelPar_prdsw,
    modelPar_prwg,
    modelPar_lprwg,
    modelPar_wprwg,
    modelPar_pprwg,
    modelPar_prwb,
    modelPar_lprwb,
    modelPar_wprwb,
    modelPar_pprwb,
    modelPar_prt,
    modelPar_lprt,
    modelPar_wprt,
    modelPar_pprt,
    modelPar_eta0,
    modelPar_leta0,
    modelPar_weta0,
    modelPar_peta0,
    modelPar_etab,
    modelPar_letab,
    modelPar_wetab,
    modelPar_petab,
    modelPar_pclm,
    modelPar_lpclm,
    modelPar_wpclm,
    modelPar_ppclm,
    modelPar_pdibl1,
    modelPar_lpdibl1,
    modelPar_wpdibl1,
    modelPar_ppdibl1,
    modelPar_pdibl2,
    modelPar_lpdibl2,
    modelPar_wpdibl2,
    modelPar_ppdibl2,
    modelPar_pdiblb,
    modelPar_lpdiblb,
    modelPar_wpdiblb,
    modelPar_ppdiblb,
    modelPar_pscbe1,
    modelPar_lpscbe1,
    modelPar_wpscbe1,
    modelPar_ppscbe1,
    modelPar_pscbe2,
    modelPar_lpscbe2,
    modelPar_wpscbe2,
    modelPar_ppscbe2,
    modelPar_pvag,
    modelPar_lpvag,
    modelPar_wpvag,
    modelPar_ppvag,
    modelPar_wr,
    modelPar_lwr,
    modelPar_wwr,
    modelPar_pwr,
    modelPar_dwg,
    modelPar_ldwg,
    modelPar_wdwg,
    modelPar_pdwg,
    modelPar_dwb,
    modelPar_ldwb,
    modelPar_wdwb,
    modelPar_pdwb,
    modelPar_b0,
    modelPar_lb0,
    modelPar_wb0,
    modelPar_pb0,
    modelPar_b1,
    modelPar_lb1,
    modelPar_wb1,
    modelPar_pb1,
    modelPar_alpha0,
    modelPar_lalpha0,
    modelPar_walpha0,
    modelPar_palpha0,
    modelPar_alpha1,
    modelPar_lalpha1,
    modelPar_walpha1,
    modelPar_palpha1,
    modelPar_beta0,
    modelPar_lbeta0,
    modelPar_wbeta0,
    modelPar_pbeta0,
    modelPar_elm,
    modelPar_lelm,
    modelPar_welm,
    modelPar_pelm,
    modelPar_cgsl,
    modelPar_lcgsl,
    modelPar_wcgsl,
    modelPar_pcgsl,
    modelPar_cgdl,
    modelPar_lcgdl,
    modelPar_wcgdl,
    modelPar_pcgdl,
    modelPar_ckappa,
    modelPar_lckappa,
    modelPar_wckappa,
    modelPar_pckappa,
    modelPar_cf,
    modelPar_lcf,
    modelPar_wcf,
    modelPar_pcf,
    modelPar_clc,
    modelPar_lclc,
    modelPar_wclc,
    modelPar_pclc,
    modelPar_cle,
    modelPar_lcle,
    modelPar_wcle,
    modelPar_pcle,
    modelPar_vfbcv,
    modelPar_lvfbcv,
    modelPar_wvfbcv,
    modelPar_pvfbcv,
    modelPar_acde,
    modelPar_lacde,
    modelPar_wacde,
    modelPar_pacde,
    modelPar_moin,
    modelPar_lmoin,
    modelPar_wmoin,
    modelPar_pmoin,
    modelPar_noff,
    modelPar_lnoff,
    modelPar_wnoff,
    modelPar_pnoff,
    modelPar_voffcv,
    modelPar_lvoffcv,
    modelPar_wvoffcv,
    modelPar_pvoffcv,
    binUnit, // this is an int //modelPar_binUnit,
    modelPar_tox,
    modelPar_given_K1, //modelPar_k1Given,
    modelPar_given_K2, //modelPar_k2Given,
    modelPar_given_NSUB, //modelPar_nsubGiven,
    modelPar_given_XT, //modelPar_xtGiven,
    modelPar_given_VBX, //modelPar_vbxGiven,
    modelPar_given_GAMMA2, //modelPar_gamma2Given,
    modelPar_given_VFB, //modelPar_vfbGiven,
    modelPar_given_VTH0, //modelPar_vth0Given,
    dtype, // this is an int //modelPar_dtype,
    modelPar_toxm,
    modelPar_factor1,
    // outputs
    instance_vtm, 
    instancePar_temp, // this is goofy ...  //temp,
    instance_jctTempSatCurDensity,
    instance_jctSidewallTempSatCurDensity,
    instance_unitAreaJctCapTemp,
    instance_unitLengthSidewallJctCapTemp,
    instance_unitLengthGateSidewallJctCapTemp,
    instance_PhiBTemp,
    instance_cgso,
    instance_cgdo,
    instance_vjsm,
    instance_IsEvjsm,
    instance_vjdm,
    instance_IsEvjdm,
    instance_updateTemperatureCalled_,
    sizeDepParams
    );
  //
    double * solVec = in.extData.nextSolVectorRawPtr;

    fadType Vd = (solVec)[in.li_Drain];
    fadType Vg = (solVec)[in.li_Gate];
    fadType Vs = (solVec)[in.li_Source];
    fadType Vb = (solVec)[in.li_Bulk];
    fadType Vsp = (solVec)[in.li_SourcePrime];
    fadType Vdp = (solVec)[in.li_DrainPrime];
    fadType Qtotal = 0.0;
    if( in.nqsMod ) { Qtotal = (solVec)[in.li_Charge]; }

    fadType Vddp  = Vd   - Vdp;
    fadType Vssp  = Vs   - Vsp;
    fadType Vbsp  = Vb   - Vsp;
    fadType Vbdp  = Vb   - Vdp;
    fadType Vgsp  = Vg   - Vsp;
    fadType Vgdp  = Vg   - Vdp;
    fadType Vgb   = Vg   - Vb;
    fadType Vdpsp = Vdp  - Vsp;

    fadType vbs  = mod.dtype * Vbsp;
    fadType vgs  = mod.dtype * Vgsp;
    fadType vds  = mod.dtype * Vdpsp;
    fadType qdef = mod.dtype * Qtotal;
    fadType vbd = vbs - vds;
    fadType vgd = vgs - vds;

    fadType instance_T1global = in.T1global;
    fadType instance_thetavth = in.thetavth;

    fadType instance_Vth = in.Vth;
    fadType instance_von = in.von;
    fadType instance_dVgs_eff_dVg = in.dVgs_eff_dVg;
    fadType instance_Vgsteff = in.Vgsteff;
    fadType instance_Abulk = in.Abulk;
    fadType instance_ueff = in.ueff;

    fadType instance_Vdsat = in.Vdsat;
    fadType instance_vdsat = in.vdsat;
    fadType instance_Vdseff = in.Vdseff;
    fadType instance_Gm = in.Gm;
    fadType instance_cdrain = in.cdrain;
    fadType instance_gds = in.gds;
    fadType instance_gm = in.gm;
    fadType instance_gmbs = in.gmbs;
    fadType instance_gbbs = in.gbbs;
    fadType instance_gbgs = in.gbgs;
    fadType instance_gbds = in.gbds;
    fadType instance_csub = in.csub;
    fadType instance_qgate = in.qgate;
    fadType instance_qdrn = in.qdrn;
    fadType instance_qsrc = in.qsrc;
    fadType instance_qbulk = in.qbulk;
    fadType instance_cggb = in.cggb;
    fadType instance_cgsb = in.cgsb;
    fadType instance_cgdb = in.cgdb;
    fadType instance_cdgb = in.cdgb;
    fadType instance_cdsb = in.cdsb;
    fadType instance_cddb = in.cddb;
    fadType instance_cbgb = in.cbgb;
    fadType instance_cbsb = in.cbsb;
    fadType instance_cbdb = in.cbdb;
    fadType instance_cqdb = in.cqdb;
    fadType instance_cqsb = in.cqsb;
    fadType instance_cqgb = in.cqgb;
    fadType instance_cqbb = in.cqbb;
    fadType instance_gtau = in.gtau;

    fadType instance_dVgst_dVb = in.dVgst_dVb;
    fadType instance_dVgst_dVg = in.dVgst_dVg;
    fadType instance_CoxWL = in.CoxWL;
    fadType instance_qinv = in.qinv;
    fadType instance_Cgg = in.Cgg;
    fadType instance_Cgd = in.Cgd;
    fadType instance_Cgb = in.Cgb;
    fadType instance_Cbg = in.Cbg;
    fadType instance_Cbd = in.Cbd;
    fadType instance_Cbb = in.Cbb;
    fadType instance_Csg = in.Csg;
    fadType instance_Csb = in.Csb;
    fadType instance_Csd = in.Csd;
    fadType instance_dDeltaPhi_dVg = in.dDeltaPhi_dVg;
    fadType instance_dDeltaPhi_dVd = in.dDeltaPhi_dVd;
    fadType instance_dDeltaPhi_dVb = in.dDeltaPhi_dVb;
    fadType instance_cd = in.cd;

    fadType instance_qbs = in.qbs;
    fadType instance_capbs = in.capbs;

    fadType instance_qbd = in.qbd;
    fadType instance_capbd = in.capbd;
    //fadType instance_nqsMod = in.nqsMod;
    fadType instance_qcheq = in.qcheq;

    fadType instance_qgso = in.qgso;
    fadType instance_qgdo = in.qgdo;
    fadType instance_cqdef = in.cqdef;
    fadType instance_qg = in.qg;
    fadType instance_qd = in.qd;
    fadType instance_qb = in.qb;
    fadType instance_qcdump = in.qcdump;
    fadType instance_Idrain = in.Idrain;
    fadType instance_drainConductance = in.drainConductance;
    fadType instance_Isource = in.Isource;
    fadType instance_sourceConductance = in.sourceConductance;

    fadType instance_gbs = in.gbs;
    fadType instance_cbs = in.cbs;
    fadType instance_gbd = in.gbd;
    fadType instance_cbd = in.cbd;

    fadType instance_CAPcggb = in.CAPcggb;
    fadType instance_CAPcgdb = in.CAPcgdb;
    fadType instance_CAPcgsb = in.CAPcgsb;
    fadType instance_CAPcdgb = in.CAPcdgb;
    fadType instance_CAPcddb = in.CAPcddb;
    fadType instance_CAPcdsb = in.CAPcdsb;
    fadType instance_CAPcsgb = in.CAPcsgb;
    fadType instance_CAPcsdb = in.CAPcsdb;
    fadType instance_CAPcssb = in.CAPcssb;
    fadType instance_CAPcbgb = in.CAPcbgb;
    fadType instance_CAPcbdb = in.CAPcbdb;
    fadType instance_CAPcbsb = in.CAPcbsb;

    bool ChargeComputationNeeded;

    bool UIVsuccess = updateIntermediateVars<fadType> 
      (
        in.getSolverState(),
        in.getDeviceOptions(),
        sizeDepParams,

        mod.dtype,
        mobMod,
        capMod,

         vbs,
         vds,
         vgs,
         instance_vtm,
         modelPar_jctEmissionCoeff,
         modelPar_ijth,
         modelPar_factor1,
         modelPar_tnom,
         modelPar_tox,
         modelPar_cox,
         modelPar_xpart,
         modelPar_bulkJctBotGradingCoeff,
         modelPar_bulkJctSideGradingCoeff,
         modelPar_bulkJctGateSideGradingCoeff,


          instancePar_sourceArea, 
          instancePar_sourcePerimeter,
          instancePar_drainArea,
          instancePar_drainPerimeter,

         instance_jctTempSatCurDensity,
         instance_jctSidewallTempSatCurDensity,
          instance_vjsm,
          instance_IsEvjsm,
          instance_vjdm,
          instance_IsEvjdm,

          in.mode,

          instance_T1global,
          instance_thetavth,
          instancePar_temp,
          instance_Vth,
          instance_von,
          instance_dVgs_eff_dVg,
          instance_Vgsteff,
          instance_Abulk,
          instance_ueff,
          instance_Vdsat,
          instance_vdsat,
          instance_Vdseff,
          instance_Gm,
          instance_cdrain,
          instance_gds,
          instance_gm,

          instance_gmbs,
          instance_gbbs,
          instance_gbgs,
          instance_gbds,
          instance_csub,
          instance_qgate,
          instance_qdrn,
          instance_qsrc,
          instance_qbulk,
          instance_cggb,
          instance_cgsb,
          instance_cgdb,
          instance_cdgb,
          instance_cdsb,
          instance_cddb,
          instance_cbgb,
          instance_cbsb,
          instance_cbdb,
          instance_cqdb,
          instance_cqsb,
          instance_cqgb,
          instance_cqbb,
          instance_gtau,

          instance_dVgst_dVb,
          instance_dVgst_dVg,
          instance_CoxWL,
          instance_qinv,
          instance_Cgg,
          instance_Cgd,
          instance_Cgb,
          instance_Cbg,
          instance_Cbd,
          instance_Cbb,
          instance_Csg,
          instance_Csb,
          instance_Csd,
          instance_dDeltaPhi_dVg,
          instance_dDeltaPhi_dVd,
          instance_dDeltaPhi_dVb,
          instance_cd,

          instance_unitAreaJctCapTemp,
          instance_unitLengthGateSidewallJctCapTemp,
          instance_unitLengthSidewallJctCapTemp,

          instance_qbs,
          instance_capbs,

          instance_PhiBTemp,
          instance_PhiBSWTemp,
          instance_PhiBSWGTemp,

          instance_qbd,
          instance_capbd,
          in.nqsMod, //instance_nqsMod,
          instance_qcheq,

          model_vtm,

          instance_cgdo,
          instance_qgdo,
          instance_cgso,
          instance_qgso,
          instance_cqdef,
          instance_qg,
          instance_qd,
          instance_qb,
          instance_qcdump,

          qdef,
          instance_Idrain,
          instance_drainConductance,
          Vddp,
          instance_Isource,
          instance_sourceConductance,
          Vssp,

        //outputs
         instance_gbs,
         instance_cbs,
         instance_gbd,
         instance_cbd,

         instance_CAPcggb,
         instance_CAPcgdb,
         instance_CAPcgsb,
         instance_CAPcdgb,
         instance_CAPcddb,
         instance_CAPcdsb,
         instance_CAPcsgb,
         instance_CAPcsdb,
         instance_CAPcssb,
         instance_CAPcbgb,
         instance_CAPcbdb,
         instance_CAPcbsb,
         ChargeComputationNeeded
      );

    // loads, and some related calcs
    fadType Gm, Gmbs, FwdSum, RevSum, cdreq, ceqbd, ceqbs, gbbdp, gbbsp;
    fadType gbdpg, gbdpdp, gbdpb, gbdpsp, gbspg, gbspdp, gbspb, gbspsp;
    if (in.mode >= 0)
    {
      Gm = instance_gm;
      Gmbs = instance_gmbs;
      FwdSum = Gm + Gmbs;
      RevSum = 0.0;

      cdreq =  mod.dtype * (instance_cd);
      ceqbd = -mod.dtype * (instance_csub);

      ceqbs = 0.0;

      gbbdp = -instance_gbds;
      gbbsp = (instance_gbds + instance_gbgs + instance_gbbs);

      gbdpg = instance_gbgs;
      gbdpdp = instance_gbds;
      gbdpb = instance_gbbs;
      gbdpsp = -(gbdpg + gbdpdp + gbdpb);

      gbspg = 0.0;
      gbspdp = 0.0;
      gbspb = 0.0;
      gbspsp = 0.0;
    }
    else
    {
      Gm = -instance_gm;
      Gmbs = -instance_gmbs;
      FwdSum = 0.0;
      RevSum = -(Gm + Gmbs);

      cdreq = -mod.dtype * (instance_cd);
      ceqbs = -mod.dtype * (instance_csub);

      ceqbd = 0.0;

      gbbsp = -instance_gbds;
      gbbdp = (instance_gbds + instance_gbgs + instance_gbbs);

      gbdpg = 0.0;
      gbdpsp = 0.0;
      gbdpb = 0.0;
      gbdpdp = 0.0;

      gbspg = instance_gbgs;
      gbspsp = instance_gbds;
      gbspb = instance_gbbs;
      gbspdp = -(gbspg + gbspsp + gbspb);
    }

    if (mod.dtype > 0)
    {
      ceqbs += (instance_cbs);
      ceqbd += (instance_cbd);
    }
    else
    {
      ceqbs -= (instance_cbs);
      ceqbd -= (instance_cbd);
    }

    int local_Drain = 0+inst*10;
    int local_Gate = 1+inst*10;
    int local_Source = 2+inst*10;
    int local_Bulk = 3+inst*10;
    int local_DrainPrime = 4+inst*10;
    int local_SourcePrime = 5+inst*10;
    int local_Charge = 6+inst*10;
    int local_Ibs = 7+inst*10;
    int local_Ids = 8+inst*10;
    int local_Igs = 9+inst*10;

    if (in.drainConductance != 0.0) { dfdp[local_Drain] += instance_Idrain.dx(0)*in.numberParallel; }
    if (in.sourceConductance != 0.0) { dfdp[local_Source] += instance_Isource.dx(0)*in.numberParallel; }

    dfdp[local_Bulk] += (ceqbs.dx(0) + ceqbd.dx(0))*in.numberParallel;
    dfdp[local_DrainPrime] += (-(ceqbd.dx(0) - cdreq.dx(0))-instance_Idrain.dx(0))*in.numberParallel;
    dfdp[local_SourcePrime] += (-(cdreq.dx(0) + ceqbs.dx(0))-instance_Isource.dx(0))*in.numberParallel;

  // don't bother with the IC stuff

    fadType ScalingFactor = 1.0e-9; // originally declared in updateIntermediateVars.  This is a copy

    auxChargeCalculations (
         ChargeComputationNeeded,
         in.mode,
         in.nqsMod,
         model_vtm,
         sizeDepParams,
      ScalingFactor,
         instance_gtau
        );

    fadType Qeqqg = 0.0;   // gate charge
    fadType Qeqqb = 0.0;   // bulk charge
    fadType Qeqqd = 0.0;   // drain charge
    fadType Qqdef = 0.0;   // nqs-related charge.
    fadType Qqcheq = 0.0;  // nqs-related charge.
    if (in.model_.dtype > 0)
    {
      Qeqqg  = instance_qg;
      Qeqqb  = instance_qb;
      Qeqqd  = instance_qd;
      Qqdef  = instance_qcdump; // this needs to be fixed...
      Qqcheq = instance_qcheq;
    }
    else  // need to convert these to charges.
    {
      Qeqqg  = -instance_qg;
      Qeqqb  = -instance_qb;
      Qeqqd  = -instance_qd;
      Qqdef  = -instance_qcdump;
      Qqcheq = -instance_qcheq;
    }

    dqdp[local_Gate] += Qeqqg.dx(0)*in.numberParallel;
    dqdp[local_Bulk] += (Qeqqb.dx(0))*in.numberParallel;
    dqdp[local_DrainPrime] += (-(-Qeqqd.dx(0)))*in.numberParallel;
    dqdp[local_SourcePrime] += (-(+ Qeqqg.dx(0) + Qeqqb.dx(0) + Qeqqd.dx(0)))*in.numberParallel;

    if (in.nqsMod)
    {
      // 7 equ. for nqs modification. charge equation.
      dqdp[local_Charge] += -(Qqcheq.dx(0) - Qqdef.dx(0))*in.numberParallel;
    }

    Findices[local_Gate] = in.li_Gate;
    Findices[local_Source] = in.li_Source;
    Findices[local_Bulk] = in.li_Bulk;
    Findices[local_DrainPrime] = in.li_DrainPrime;
    Findices[local_SourcePrime] = in.li_SourcePrime;
    Findices[local_Charge] = in.li_Charge;
    Findices[local_Ibs] = in.li_Ibs;
    Findices[local_Ids] = in.li_Ids;
    Findices[local_Igs] = in.li_Igs;

    Qindices[local_Drain] = in.li_Drain;
    Qindices[local_Gate] = in.li_Gate;
    Qindices[local_Source] = in.li_Source;
    Qindices[local_Bulk] = in.li_Bulk;
    Qindices[local_DrainPrime] = in.li_DrainPrime;
    Qindices[local_SourcePrime] = in.li_SourcePrime;
    Qindices[local_Charge] = in.li_Charge;
    Qindices[local_Ibs] = in.li_Ibs;
    Qindices[local_Ids] = in.li_Ids;
    Qindices[local_Igs] = in.li_Igs;

  //
  }
}


} // namespace MOSFET_B3
} // namespace Device
} // namespace Xyce


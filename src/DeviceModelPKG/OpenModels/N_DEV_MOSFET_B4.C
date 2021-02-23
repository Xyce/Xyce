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
// Purpose        : This file implements the BSIM4 MOSFET model.  It
//                  is intended to be compatible with the Berkeley SPICE
//                  (3f5) version, BSIM4 version 4.6.1
//
// Special Notes  : Updated from BSIM 4.5.0 to version 4.6.1 by TVR, August 07
//
//
// Creator        : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 11/25/06
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MOSFET1.h>
#include <N_DEV_MOSFET_B4.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Message.h>
#include <N_DEV_SolverState.h>
#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_MachDepParams.h>
#include <N_UTL_Math.h>
#include <N_ANP_NoiseData.h>

// ---------- BSIM4 constants ---------------
// Many constants are obtained from N_DEV_Const.h
// A few overrides (2) are below, to make this device as close
// to the original spice3f5 code as possible.

#define CONSTEPS0 8.85418e-12
#if 0
#define CONSTMAX_EXP 5.834617425e+14
#define CONSTMIN_EXP 1.713908431e-15

#define CONSTEXP_THRESHOLD 34.0

#define CONSTDELTA_1 0.02
#define CONSTDELTA_2 0.02
#define CONSTDELTA_3 0.02
#define CONSTDELTA_4 0.02
#endif

#define Charge_q     (1.60219e-19) // electron charge, used in the
                                   // updateTemperature function instead of
                                   // of CONSTQ.  Having 2 constants for
                                   // the same quantity (with different
                                   // precision) doesn't make sense, but
                                   // it is what is in the spice3f5 bsim4.

#define CONSTKboQ 8.617087e-5  // another updateTemperature constant, which is
                               // not really necessary but I am keeping for
                               // compatibility with the spice3f5 bsim4.

#define CONSTvt0     (CONSTboltz * (27.0 +CONSTCtoK)/CONSTQ)

#define CONSTMM  3  // smooth coeff

#define DEXP(A,B,C) {                                                         \
        if (A > CONSTEXP_THRESHOLD) {                                         \
            B = CONSTMAX_EXP*(1.0+(A)-CONSTEXP_THRESHOLD);                    \
            C = CONSTMAX_EXP;                                                 \
        } else if (A < -CONSTEXP_THRESHOLD)  {                                \
            B = CONSTMIN_EXP;                                                 \
            C = 0;                                                            \
        } else   {                                                            \
            B = exp(A);                                                       \
            C = B;                                                            \
        }                                                                     \
    }


#define DELTA  1.0E-9
#define DEXP2(A,B) {                                                       \
        if (A > CONSTEXP_THRESHOLD) {                                      \
            B = CONSTMAX_EXP*(1.0+(A)-CONSTEXP_THRESHOLD);                 \
        } else if (A < -CONSTEXP_THRESHOLD)  {                             \
            B = CONSTMIN_EXP;                                              \
        } else   {                                                         \
            B = exp(A);                                                    \
        }                                                                  \
    }


namespace Xyce {
namespace Device {

namespace MOSFET_B4 {

void Traits::loadInstanceParameters(ParametricData<MOSFET_B4::Instance> &p)
{
    p.addPar ("TEMP",0.0,&MOSFET_B4::Instance::temp)
     .setGivenMember(&MOSFET_B4::Instance::TEMPgiven)
     .setExpressionAccess(ParameterType::TIME_DEP)
     .setUnit(STANDARD)
     .setCategory(CAT_NONE)
     .setDescription("Device temperature");

    p.addPar ("L",5.0e-6,&MOSFET_B4::Instance::l)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length");

    p.addPar ("W",5.0e-6,&MOSFET_B4::Instance::w)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width");

    p.addPar ("NF",1.0,&MOSFET_B4::Instance::nf)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Number of fingers");

    p.addPar ("SA",0.0,&MOSFET_B4::Instance::sa)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("distance between  OD edge to poly of one side ");

    p.addPar ("SB",0.0,&MOSFET_B4::Instance::sb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("distance between  OD edge to poly of the other side");

    p.addPar ("SD",0.0,&MOSFET_B4::Instance::sd)
     .setGivenMember(&MOSFET_B4::Instance::SDgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("distance between neighbour fingers");

    p.addPar ("SCA",0.0,&MOSFET_B4::Instance::sca)
     .setGivenMember(&MOSFET_B4::Instance::scaGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Integral of the first distribution function for scattered well dopant");

    p.addPar ("SCB",0.0,&MOSFET_B4::Instance::scb)
     .setGivenMember(&MOSFET_B4::Instance::scbGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Integral of the second distribution function for scattered well dopant");

    p.addPar ("SCC",0.0,&MOSFET_B4::Instance::scc)
     .setGivenMember(&MOSFET_B4::Instance::sccGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Integral of the third distribution function for scattered well dopant");

    p.addPar ("SC",0.0,&MOSFET_B4::Instance::sc)
     .setGivenMember(&MOSFET_B4::Instance::scGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Distance to a single well edge ");


    p.addPar ("AD",0.0,&MOSFET_B4::Instance::drainArea)
     .setGivenMember(&MOSFET_B4::Instance::drainAreaGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain area");

    p.addPar ("AS",0.0,&MOSFET_B4::Instance::sourceArea)
     .setGivenMember(&MOSFET_B4::Instance::sourceAreaGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source area");

    p.addPar ("PD",0.0,&MOSFET_B4::Instance::drainPerimeter)
     .setGivenMember(&MOSFET_B4::Instance::drainPerimeterGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain perimeter");


    p.addPar ("PS",0.0,&MOSFET_B4::Instance::sourcePerimeter)
     .setGivenMember(&MOSFET_B4::Instance::sourcePerimeterGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source perimeter");

    p.addPar ("NRD",1.0,&MOSFET_B4::Instance::drainSquares)
     .setGivenMember(&MOSFET_B4::Instance::drainSquaresGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_ASYMRDS)
     .setDescription("Number of squares in drain");

    p.addPar ("NRS",1.0,&MOSFET_B4::Instance::sourceSquares)
     .setGivenMember(&MOSFET_B4::Instance::sourceSquaresGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_ASYMRDS)
     .setDescription("Number of squares in source");


    p.addPar ("RBDB",0.0,&MOSFET_B4::Instance::rbdb)
     .setGivenMember(&MOSFET_B4::Instance::RBDBgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance");

    p.addPar ("RBSB",0.0,&MOSFET_B4::Instance::rbsb)
     .setGivenMember(&MOSFET_B4::Instance::RBSBgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance");

    p.addPar ("RBPB",0.0,&MOSFET_B4::Instance::rbpb)
     .setGivenMember(&MOSFET_B4::Instance::RBPBgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance");

    p.addPar ("RBPS",0.0,&MOSFET_B4::Instance::rbps)
     .setGivenMember(&MOSFET_B4::Instance::RBPSgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance");

    p.addPar ("RBPD",0.0,&MOSFET_B4::Instance::rbpd)
     .setGivenMember(&MOSFET_B4::Instance::RBPDgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance");

    p.addPar ("DELVTO",0.0,&MOSFET_B4::Instance::delvto)
     .setUnit(U_VOLT)
     .setCategory(CAT_BASIC)
     .setDescription("Zero bias threshold voltage variation");

    p.addPar ("XGW",0.0,&MOSFET_B4::Instance::xgw)
     .setGivenMember(&MOSFET_B4::Instance::XGWgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Distance from gate contact center to device edge");

    p.addPar ("NGCON",0.0,&MOSFET_B4::Instance::ngcon)
     .setGivenMember(&MOSFET_B4::Instance::NGCONgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Number of gate contacts");

    p.addPar ("M",1.0,&MOSFET_B4::Instance::numberParallel)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Number of parallel copies");

    p.addPar ("IC1",0.0,&MOSFET_B4::Instance::icVDS)
     .setGivenMember(&MOSFET_B4::Instance::icVDSGiven)
     .setUnit(U_VOLT)
     .setCategory(CAT_VOLT)
     .setDescription("Vector of initial values: Vds,Vgs,Vbs");

    p.addPar ("IC2",0.0,&MOSFET_B4::Instance::icVGS)
     .setGivenMember(&MOSFET_B4::Instance::icVGSGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("");

    p.addPar ("IC3",0.0,&MOSFET_B4::Instance::icVBS)
     .setGivenMember(&MOSFET_B4::Instance::icVBSGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("");

    // Set up non-double precision variables:
    p.addPar ("TRNQSMOD",0,&MOSFET_B4::Instance::trnqsMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Transient NQS model selector");

    p.addPar ("ACNQSMOD",0,&MOSFET_B4::Instance::acnqsMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("AC NQS model selector");

    p.addPar ("RBODYMOD",0,&MOSFET_B4::Instance::rbodyMod)
     .setGivenMember(&MOSFET_B4::Instance::RBODYMODgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Distributed body R model selector");

    p.addPar ("RGATEMOD",0,&MOSFET_B4::Instance::rgateMod)
     .setGivenMember(&MOSFET_B4::Instance::RGATEMODgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Gate resistance model selector");

    p.addPar ("GEOMOD",0,&MOSFET_B4::Instance::geoMod)
     .setGivenMember(&MOSFET_B4::Instance::GEOMODgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Geometry dependent parasitics model selector");

    p.addPar ("RGEOMOD",0,&MOSFET_B4::Instance::rgeoMod)
     .setGivenMember(&MOSFET_B4::Instance::RGEOMODgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("S/D resistance and contact model selector");

    p.addPar ("MIN",0,&MOSFET_B4::Instance::min)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Minimize either D or S");

    p.addPar ("OFF",false,&MOSFET_B4::Instance::OFF)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Device is initially off");

    // This tells the parser that IC1,IC2,and IC3 are to be input as a vector of "IC"
    p.makeVector ("IC",3);
}

void Traits::loadModelParameters(ParametricData<MOSFET_B4::Model> &p)
{
    p.addPar ("EOT",15.0e-10,&MOSFET_B4::Model::eot)
     .setUnit(U_METER)
     .setCategory(CAT_PROCESS)
     .setDescription("Equivalent gate oxide thickness in meters");

    p.addPar ("VDDEOT",1.5,&MOSFET_B4::Model::vddeot)
     .setUnit(U_VOLT)
     .setCategory(CAT_BASIC)
     .setDescription("Voltage for extraction of equivalent gate oxide thickness");

    p.addPar ("ADOS",1.0,&MOSFET_B4::Model::ados)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("Charge centroid parameter");

    p.addPar ("BDOS",1.0,&MOSFET_B4::Model::bdos)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("Charge centroid parameter");

    p.addPar ("TOXE",30.0e-10,&MOSFET_B4::Model::toxe)
     .setGivenMember(&MOSFET_B4::Model::toxeGiven)
     .setUnit(U_METER)
     .setCategory(CAT_PROCESS)
     .setDescription("Electrical gate oxide thickness in meters");

    p.addPar ("TOXP",30.0e-10,&MOSFET_B4::Model::toxp)
     .setGivenMember(&MOSFET_B4::Model::toxpGiven)
     .setUnit(U_METER)
     .setCategory(CAT_PROCESS)
     .setDescription("Physical gate oxide thickness in meters");

    p.addPar ("TOXM",30.0e-10,&MOSFET_B4::Model::toxm)
     .setUnit(U_METER)
     .setCategory(CAT_PROCESS)
     .setDescription("Gate oxide thickness at which parameters are extracted");

    p.addPar ("TOXREF",30.0e-10,&MOSFET_B4::Model::toxref)
     .setUnit(U_METER)
     .setCategory(CAT_TUNNEL)
     .setDescription("Target tox value");

    p.addPar ("DTOX",0.0,&MOSFET_B4::Model::dtox)
     .setGivenMember(&MOSFET_B4::Model::dtoxGiven)
     .setUnit(U_METER)
     .setCategory(CAT_PROCESS)
     .setDescription("Defined as (toxe - toxp) ");

    p.addPar ("EPSROX",3.9,&MOSFET_B4::Model::epsrox)
     .setUnit(U_NONE)
     .setCategory(CAT_PROCESS)
     .setDescription("Dielectric constant of the gate oxide relative to vacuum");


    p.addPar ("CDSC",2.4e-4,&MOSFET_B4::Model::cdsc)
     .setUnit(U_FARADMM2)
     .setCategory(CAT_BASIC)
     .setDescription("Drain/Source and channel coupling capacitance");

    p.addPar ("CDSCB",0.0,&MOSFET_B4::Model::cdscb)
     .setUnit(U_FVM1MM2)
     .setCategory(CAT_BASIC)
     .setDescription("Body-bias dependence of cdsc");

    p.addPar ("CDSCD",0.0,&MOSFET_B4::Model::cdscd)
     .setUnit(U_FVM1MM2)
     .setCategory(CAT_BASIC)
     .setDescription("Drain-bias dependence of cdsc");

    p.addPar ("CIT",0.0,&MOSFET_B4::Model::cit)
     .setUnit(U_FARADMM2)
     .setCategory(CAT_BASIC)
     .setDescription("Interface state capacitance");

    p.addPar ("NFACTOR",1.0,&MOSFET_B4::Model::nfactor)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("Subthreshold swing Coefficient");

    p.addPar ("XJ",0.15e-6,&MOSFET_B4::Model::xj)
     .setUnit(U_METER)
     .setCategory(CAT_PROCESS)
     .setDescription("Junction depth in meters");

    p.addPar ("VSAT",8.0e4,&MOSFET_B4::Model::vsat)
     .setUnit(U_MSM1)
     .setCategory(CAT_BASIC)
     .setDescription("Saturation velocity at tnom");

    p.addPar ("AT",3.3e4,&MOSFET_B4::Model::at)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient of vsat");

    p.addPar ("A0",1.0,&MOSFET_B4::Model::a0)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("Non-uniform depletion width effect coefficient.");

    p.addPar ("AGS",0.0,&MOSFET_B4::Model::ags)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_BASIC)
     .setDescription("Gate bias  coefficient of Abulk.");

    p.addPar ("A1",0.0,&MOSFET_B4::Model::a1)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_BASIC)
     .setDescription("Non-saturation effect coefficient");

    p.addPar ("A2",1.0,&MOSFET_B4::Model::a2)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("Non-saturation effect coefficient");

    p.addPar ("KETA",-0.047,&MOSFET_B4::Model::keta)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_BASIC)
     .setDescription("Body-bias coefficient of non-uniform depletion width effect.");

    p.addPar ("PHIG",4.05,&MOSFET_B4::Model::phig)
     .setGivenMember(&MOSFET_B4::Model::phigGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Work Function of gate");

    p.addPar ("EPSRGATE",11.7,&MOSFET_B4::Model::epsrgate)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Dielectric constant of gate relative to vacuum");

    p.addPar ("EASUB",4.05,&MOSFET_B4::Model::easub)
     .setUnit(U_VOLT)
     .setCategory(CAT_BASIC)
     .setDescription("Electron affinity of substrate");

    p.addPar ("EPSRSUB",11.7,&MOSFET_B4::Model::epsrsub)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("Dielectric constant of substrate relative to vacuum");

    p.addPar ("NI0SUB",1.45e10,&MOSFET_B4::Model::ni0sub)
     .setUnit(U_CMM3)
     .setCategory(CAT_BASIC)
     .setDescription("Intrinsic carrier concentration of substrate at 300.15K");

    p.addPar ("BG0SUB",1.16,&MOSFET_B4::Model::bg0sub)
     .setUnit(U_EV)
     .setCategory(CAT_BASIC)
     .setDescription("Band-gap of substrate at T=0K");

    p.addPar ("TBGASUB",7.02e-4,&MOSFET_B4::Model::tbgasub)
     .setUnit(U_EVDEGKM1)
     .setCategory(CAT_BASIC)
     .setDescription("First parameter of band-gap change due to temperature");


    p.addPar ("TBGBSUB",1108.0,&MOSFET_B4::Model::tbgbsub)
     .setUnit(U_DEGK)
     .setCategory(CAT_BASIC)
     .setDescription("Second parameter of band-gap change due to temperature");

    p.addPar ("NSUB",6.0e16,&MOSFET_B4::Model::nsub)
     .setGivenMember(&MOSFET_B4::Model::nsubGiven)
     .setUnit(U_CMM3)
     .setCategory(CAT_PROCESS)
     .setDescription("Substrate doping concentration");

    p.addPar ("NDEP",1.7e17,&MOSFET_B4::Model::ndep)
     .setGivenMember(&MOSFET_B4::Model::ndepGiven)
     .setUnit(U_CMM3)
     .setCategory(CAT_PROCESS)
     .setDescription("Channel doping concentration at the depletion edge");

    p.addPar ("NSD",1.0e20,&MOSFET_B4::Model::nsd)
     .setUnit(U_CMM3)
     .setCategory(CAT_PROCESS)
     .setDescription("S/D doping concentration");

    p.addPar ("PHIN",0.0,&MOSFET_B4::Model::phin)
     .setUnit(U_VOLT)
     .setCategory(CAT_BASIC)
     .setDescription("Adjusting parameter for surface potential due to non-uniform vertical doping");

    p.addPar ("NGATE",0.0,&MOSFET_B4::Model::ngate)
     .setUnit(U_CMM3)
     .setCategory(CAT_PROCESS)
     .setDescription("Poly-gate doping concentration");

    p.addPar ("GAMMA1",0.0,&MOSFET_B4::Model::gamma1)
      .setGivenMember(&MOSFET_B4::Model::gamma1Given)
     .setUnit(U_VOLTH)
     .setCategory(CAT_PROCESS)
     .setDescription("Vth body coefficient");

    p.addPar ("GAMMA2",0.0,&MOSFET_B4::Model::gamma2)
      .setGivenMember(&MOSFET_B4::Model::gamma2Given)
     .setUnit(U_VOLTH)
     .setCategory(CAT_PROCESS)
     .setDescription("Vth body coefficient");

    p.addPar ("VBX",0.0,&MOSFET_B4::Model::vbx)
      .setGivenMember(&MOSFET_B4::Model::vbxGiven)
     .setUnit(U_VOLT)
     .setCategory(CAT_PROCESS)
     .setDescription("Vth transition body Voltage");

    p.addPar ("VBM",-3.0,&MOSFET_B4::Model::vbm)
     .setUnit(U_VOLT)
     .setCategory(CAT_BASIC)
     .setDescription("Maximum body voltage");

    p.addPar ("XT",1.55e-7,&MOSFET_B4::Model::xt)
     .setGivenMember(&MOSFET_B4::Model::xtGiven)
     .setUnit(U_METER)
     .setCategory(CAT_PROCESS)
     .setDescription("Doping depth");

    p.addPar ("K1",0.0,&MOSFET_B4::Model::k1)
     .setGivenMember(&MOSFET_B4::Model::k1Given)
     .setUnit(U_VOLTMH)
     .setCategory(CAT_BASIC)
     .setDescription("Bulk effect coefficient 1");

    p.addPar ("KT1",-0.11,&MOSFET_B4::Model::kt1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient of Vth");

    p.addPar ("KT1L",0.0,&MOSFET_B4::Model::kt1l)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient of Vth");

    p.addPar ("KT2",0.022,&MOSFET_B4::Model::kt2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body-coefficient of kt1");

    p.addPar ("K2",0.0,&MOSFET_B4::Model::k2)
     .setGivenMember(&MOSFET_B4::Model::k2Given)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("Bulk effect coefficient 2");

    p.addPar ("K3",80.0,&MOSFET_B4::Model::k3)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("Narrow width effect coefficient");

    p.addPar ("K3B",0.0,&MOSFET_B4::Model::k3b)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body effect coefficient of k3");

    p.addPar ("W0",2.5e-6,&MOSFET_B4::Model::w0)
     .setUnit(U_METER)
     .setCategory(CAT_BASIC)
     .setDescription("Narrow width effect parameter");

    p.addPar ("DVTP0",0.0,&MOSFET_B4::Model::dvtp0)
     .setUnit(U_METER)
     .setCategory(CAT_BASIC)
     .setDescription("First parameter for Vth shift due to pocket");

    p.addPar ("DVTP1",0.0,&MOSFET_B4::Model::dvtp1)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_BASIC)
     .setDescription("Second parameter for Vth shift due to pocket");

    p.addPar ("LPE0",1.74e-7,&MOSFET_B4::Model::lpe0)
     .setUnit(U_METER)
     .setCategory(CAT_BASIC)
     .setDescription("Equivalent length of pocket region at zero bias");

    p.addPar ("LPEB",0.0,&MOSFET_B4::Model::lpeb)
     .setUnit(U_METER)
     .setCategory(CAT_BASIC)
     .setDescription("Equivalent length of pocket region accounting for body bias");


    p.addPar ("DVT0",2.2,&MOSFET_B4::Model::dvt0)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("Short channel effect coeff. 0");

    p.addPar ("DVT1",0.53,&MOSFET_B4::Model::dvt1)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("Short channel effect coeff. 1");

    p.addPar ("DVT2",-0.032,&MOSFET_B4::Model::dvt2)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_BASIC)
     .setDescription("Short channel effect coeff. 2");

    p.addPar ("DVT0W",0.0,&MOSFET_B4::Model::dvt0w)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("Narrow Width coeff. 0");

    p.addPar ("DVT1W",5.3e6,&MOSFET_B4::Model::dvt1w)
     .setUnit(U_METERM1)
     .setCategory(CAT_BASIC)
     .setDescription("Narrow Width effect coeff. 1");

    p.addPar ("DVT2W",-0.032,&MOSFET_B4::Model::dvt2w)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_BASIC)
     .setDescription("Narrow Width effect coeff. 2");

    p.addPar ("DROUT",0.56,&MOSFET_B4::Model::drout)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("DIBL coefficient of output resistance");

    p.addPar ("DSUB",0.0,&MOSFET_B4::Model::dsub)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("DIBL coefficient in the subthreshold region");

    p.addPar ("VTH0",0.0,&MOSFET_B4::Model::vth0)
     .setGivenMember(&MOSFET_B4::Model::vth0Given)
     .setUnit(U_VOLT)
     .setCategory(CAT_BASIC)
     .setDescription("");

    p.addPar ("UA",0.0,&MOSFET_B4::Model::ua)
     .setUnit(U_MVM1)
     .setCategory(CAT_BASIC)
     .setDescription("Linear gate dependence of mobility");

    p.addPar ("UA1",1.0e-9,&MOSFET_B4::Model::ua1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient of ua");

    p.addPar ("UB",1.0e-19,&MOSFET_B4::Model::ub)
     .setUnit(U_M2VM2)
     .setCategory(CAT_BASIC)
     .setDescription("Quadratic gate dependence of mobility");

    p.addPar ("UB1",-1.0e-18,&MOSFET_B4::Model::ub1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient of ub");

    p.addPar ("UC",0.0,&MOSFET_B4::Model::uc)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_BASIC)
     .setDescription("Body-bias dependence of mobility");

    p.addPar ("UC1",0.0,&MOSFET_B4::Model::uc1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient of uc");

    p.addPar ("UD",0.0,&MOSFET_B4::Model::ud)
     .setUnit(U_METERM2)
     .setCategory(CAT_BASIC)
     .setDescription("Coulomb scattering factor of mobility");

    p.addPar ("UD1",0.0,&MOSFET_B4::Model::ud1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient of ud");

    p.addPar ("UP",0.0,&MOSFET_B4::Model::up)
     .setUnit(U_METERM2)
     .setCategory(CAT_BASIC)
     .setDescription("Channel length linear factor of mobility");

    p.addPar ("LP",1.0e-8,&MOSFET_B4::Model::lp)
     .setUnit(U_METER)
     .setCategory(CAT_BASIC)
     .setDescription("Channel length exponential factor of mobility");

    p.addPar ("U0",0.0,&MOSFET_B4::Model::u0)
     .setUnit(U_M2VM1SM1)
     .setCategory(CAT_BASIC)
     .setDescription("Low-field mobility at Tnom");

    p.addPar ("EU",0.0,&MOSFET_B4::Model::eu)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("Mobility exponent");

    p.addPar ("UTE",-1.5,&MOSFET_B4::Model::ute)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient of mobility");

    p.addPar ("VOFF",-0.08,&MOSFET_B4::Model::voff)
     .setUnit(U_VOLT)
     .setCategory(CAT_BASIC)
     .setDescription("Threshold voltage offset");

    p.addPar ("MINV",0.0,&MOSFET_B4::Model::minv)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("Fitting parameter for moderate inversion in Vgsteff");

    p.addPar ("MINVCV",0.0,&MOSFET_B4::Model::minvcv)
     .setUnit(U_NONE)
     .setCategory(CAT_CAP)
     .setDescription("Fitting parameter for moderate inversion in Vgsteffcv");

    p.addPar ("VOFFL",0.0,&MOSFET_B4::Model::voffl)
     .setUnit(U_VOLT)
     .setCategory(CAT_BASIC)
     .setDescription("Length dependence parameter for Vth offset");

    p.addPar ("VOFFCVL",0.0,&MOSFET_B4::Model::voffcvl)
     .setUnit(U_NONE)
     .setCategory(CAT_CAP)
     .setDescription("Length dependence parameter for Vth offset in CV");

    p.addPar ("TNOM",0.0,&MOSFET_B4::Model::tnom)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Parameter measurement temperature");

    p.addPar ("CGSO",0.0,&MOSFET_B4::Model::cgso)
     .setGivenMember(&MOSFET_B4::Model::cgsoGiven)
     .setUnit(U_FARADMM1)
     .setCategory(CAT_CAP)
     .setDescription("Gate-source overlap capacitance per width");

    p.addPar ("CGDO",0.0,&MOSFET_B4::Model::cgdo)
     .setGivenMember(&MOSFET_B4::Model::cgdoGiven)
     .setUnit(U_FARADMM1)
     .setCategory(CAT_CAP)
     .setDescription("Gate-drain overlap capacitance per width");

    p.addPar ("CGBO",0.0,&MOSFET_B4::Model::cgbo)
     .setGivenMember(&MOSFET_B4::Model::cgboGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_CAP)
     .setDescription("Gate-bulk overlap capacitance per length");

    p.addPar ("XPART",0.0,&MOSFET_B4::Model::xpart)
     .setUnit(U_FARADMM1)
     .setCategory(CAT_CAP)
     .setDescription("Channel charge partitioning");

    p.addPar ("DELTA",0.01,&MOSFET_B4::Model::delta)
     .setUnit(U_VOLT)
     .setCategory(CAT_BASIC)
     .setDescription("Effective Vds parameter");

    p.addPar ("RSH",0.0,&MOSFET_B4::Model::sheetResistance)
     .setUnit(U_OSQM1)
     .setCategory(CAT_PROCESS)
     .setDescription("Source-drain sheet resistance");

    p.addPar ("RDSW",200.0,&MOSFET_B4::Model::rdsw)
     .setUnit(U_OHMMICRON)
     .setCategory(CAT_ASYMRDS)
     .setDescription("Source-drain resistance per width");

    p.addPar ("RDSWMIN",0.0,&MOSFET_B4::Model::rdswmin)
     .setUnit(U_OHMMICRON)
     .setCategory(CAT_ASYMRDS)
     .setDescription("Source-drain resistance per width at high Vg");

    p.addPar ("RSW",100.0,&MOSFET_B4::Model::rsw)
     .setUnit(U_OHMMICRON)
     .setCategory(CAT_ASYMRDS)
     .setDescription("Source resistance per width");

    p.addPar ("RDW",100.0,&MOSFET_B4::Model::rdw)
     .setUnit(U_OHMMICRON)
     .setCategory(CAT_ASYMRDS)
     .setDescription("Drain resistance per width");

    p.addPar ("RDWMIN",0.0,&MOSFET_B4::Model::rdwmin)
     .setUnit(U_OHMMICRON)
     .setCategory(CAT_ASYMRDS)
     .setDescription("Drain resistance per width at high Vg");

    p.addPar ("RSWMIN",0.0,&MOSFET_B4::Model::rswmin)
     .setUnit(U_OHMMICRON)
     .setCategory(CAT_ASYMRDS)
     .setDescription("Source resistance per width at high Vg");

    p.addPar ("PRWG",1.0,&MOSFET_B4::Model::prwg)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_ASYMRDS)
     .setDescription("Gate-bias effect on parasitic resistance ");

    p.addPar ("PRWB",0.0,&MOSFET_B4::Model::prwb)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_ASYMRDS)
     .setDescription("Body-effect on parasitic resistance ");

    p.addPar ("PRT",0.0,&MOSFET_B4::Model::prt)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient of parasitic resistance ");

    p.addPar ("ETA0",0.08,&MOSFET_B4::Model::eta0)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("Subthreshold region DIBL coefficient");

    p.addPar ("ETAB",-0.07,&MOSFET_B4::Model::etab)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_BASIC)
     .setDescription("Subthreshold region DIBL coefficient");

    p.addPar ("PCLM",1.3,&MOSFET_B4::Model::pclm)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("Channel length modulation Coefficient");

    // check these:  are they PDIBL1,etc?
    p.addPar ("PDIBLC1",0.39,&MOSFET_B4::Model::pdibl1)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("Drain-induced barrier lowering coefficient");

    p.addPar ("PDIBLC2",0.0086,&MOSFET_B4::Model::pdibl2)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription("Drain-induced barrier lowering coefficient");

    p.addPar ("PDIBLCB",0.0,&MOSFET_B4::Model::pdiblb)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_BASIC)
     .setDescription("Body-effect on drain-induced barrier lowering");

    p.addPar ("FPROUT",0.0,&MOSFET_B4::Model::fprout)
     .setUnit(U_VMMH)
     .setCategory(CAT_BASIC)
     .setDescription("Rout degradation coefficient for pocket devices");

    p.addPar ("PDITS",0.0,&MOSFET_B4::Model::pdits)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_BASIC)
     .setDescription("Coefficient for drain-induced Vth shifts");

    p.addPar ("PDITSL",0.0,&MOSFET_B4::Model::pditsl)
     .setUnit(U_METERM1)
     .setCategory(CAT_BASIC)
     .setDescription("Length dependence of drain-induced Vth shifts");

    p.addPar ("PDITSD",0.0,&MOSFET_B4::Model::pditsd)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_BASIC)
     .setDescription("Vds dependence of drain-induced Vth shifts");

    p.addPar ("PSCBE1",4.24e8,&MOSFET_B4::Model::pscbe1)
     .setUnit(U_VMM1)
     .setCategory(CAT_BASIC)
     .setDescription("Substrate current body-effect coefficient");

    p.addPar ("PSCBE2",1.0e-5,&MOSFET_B4::Model::pscbe2)
     .setUnit(U_MVM1)
     .setCategory(CAT_BASIC)
     .setDescription("Substrate current body-effect coefficient");

    p.addPar ("PVAG",0.0,&MOSFET_B4::Model::pvag)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Gate dependence of output resistance parameter");

    p.addPar ("JSS",1e-4,&MOSFET_B4::Model::SjctSatCurDensity)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Bottom source junction reverse saturation current density");

    p.addPar ("JSWS",0.0,&MOSFET_B4::Model::SjctSidewallSatCurDensity)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Isolation edge sidewall source junction reverse saturation current density");

    p.addPar ("JSWGS",0.0,&MOSFET_B4::Model::SjctGateSidewallSatCurDensity)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Gate edge source junction reverse saturation current density");

    p.addPar ("PBS",1.0,&MOSFET_B4::Model::SbulkJctPotential)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source junction built-in potential");

    p.addPar ("NJS",1.0,&MOSFET_B4::Model::SjctEmissionCoeff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source junction emission coefficient");

    p.addPar ("XTIS",3.0,&MOSFET_B4::Model::SjctTempExponent)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source junction current temperature exponent");

    p.addPar ("MJS",0.5,&MOSFET_B4::Model::SbulkJctBotGradingCoeff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source bottom junction capacitance grading coefficient");

    p.addPar ("PBSWS",1.0,&MOSFET_B4::Model::SsidewallJctPotential)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source sidewall junction capacitance built in potential");

    p.addPar ("MJSWS",0.33,&MOSFET_B4::Model::SbulkJctSideGradingCoeff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source sidewall junction capacitance grading coefficient");

    p.addPar ("PBSWGS",0.0,&MOSFET_B4::Model::SGatesidewallJctPotential)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source (gate side) sidewall junction capacitance built in potential");

    p.addPar ("MJSWGS",0.33,&MOSFET_B4::Model::SbulkJctGateSideGradingCoeff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source (gate side) sidewall junction capacitance grading coefficient");

    p.addPar ("CJS",5e-4,&MOSFET_B4::Model::SunitAreaJctCap)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source bottom junction capacitance per unit area");

    p.addPar ("CJSWS",5e-10,&MOSFET_B4::Model::SunitLengthSidewallJctCap)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source sidewall junction capacitance per unit periphery");

    p.addPar ("CJSWGS",0.0,&MOSFET_B4::Model::SunitLengthGateSidewallJctCap)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source (gate side) sidewall junction capacitance per unit width");

    p.addPar ("JSD",1e-4,&MOSFET_B4::Model::DjctSatCurDensity)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Bottom drain junction reverse saturation current density");

    p.addPar ("JSWD",0.0,&MOSFET_B4::Model::DjctSidewallSatCurDensity)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Isolation edge sidewall drain junction reverse saturation current density");

    p.addPar ("JSWGD",0.0,&MOSFET_B4::Model::DjctGateSidewallSatCurDensity)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Gate edge drain junction reverse saturation current density");

    p.addPar ("PBD",1.0,&MOSFET_B4::Model::DbulkJctPotential)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain junction built-in potential");

    p.addPar ("NJD",1.0,&MOSFET_B4::Model::DjctEmissionCoeff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain junction emission coefficient");

    p.addPar ("XTID",3.0,&MOSFET_B4::Model::DjctTempExponent)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drainjunction current temperature exponent");

    p.addPar ("MJD",0.5,&MOSFET_B4::Model::DbulkJctBotGradingCoeff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain bottom junction capacitance grading coefficient");

    p.addPar ("PBSWD",1.0,&MOSFET_B4::Model::DsidewallJctPotential )
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain sidewall junction capacitance built in potential");

    p.addPar ("MJSWD",0.33,&MOSFET_B4::Model::DbulkJctSideGradingCoeff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain sidewall junction capacitance grading coefficient");

    p.addPar ("PBSWGD",0.0,&MOSFET_B4::Model::DGatesidewallJctPotential)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain (gate side) sidewall junction capacitance built in potential");

    p.addPar ("MJSWGD",0.33,&MOSFET_B4::Model::DbulkJctGateSideGradingCoeff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain (gate side) sidewall junction capacitance grading coefficient");

    p.addPar ("CJD",5e-4,&MOSFET_B4::Model::DunitAreaJctCap)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain bottom junction capacitance per unit area");

    p.addPar ("CJSWD",5e-10,&MOSFET_B4::Model::DunitLengthSidewallJctCap)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain sidewall junction capacitance per unit periphery");

    p.addPar ("CJSWGD",0.0,&MOSFET_B4::Model::DunitLengthGateSidewallJctCap)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain (gate side) sidewall junction capacitance per unit width");

    p.addPar ("VFBCV",-1.0,&MOSFET_B4::Model::vfbcv)
     .setUnit(U_VOLT)
     .setCategory(CAT_CAP)
     .setDescription("Flat Band Voltage parameter for capmod=0 only");

    p.addPar ("VFB",-1.0,&MOSFET_B4::Model::vfb)
     .setGivenMember(&MOSFET_B4::Model::vfbGiven)
     .setUnit(U_VOLT)
     .setCategory(CAT_BASIC)
     .setDescription("Flat Band Voltage");

    p.addPar ("TPB",0.0,&MOSFET_B4::Model::tpb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient of pb");

    p.addPar ("TCJ",0.0,&MOSFET_B4::Model::tcj)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient of cj");

    p.addPar ("TPBSW",0.0,&MOSFET_B4::Model::tpbsw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient of pbsw");

    p.addPar ("TCJSW",0.0,&MOSFET_B4::Model::tcjsw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient of cjsw");

    p.addPar ("TPBSWG",0.0,&MOSFET_B4::Model::tpbswg)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient of pbswg");

    p.addPar ("TCJSWG",0.0,&MOSFET_B4::Model::tcjswg)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient of cjswg");

    p.addPar ("ACDE",1.0,&MOSFET_B4::Model::acde)
     .setUnit(U_MVM1)
     .setCategory(CAT_CAP)
     .setDescription("Exponential coefficient for finite charge thickness");

    p.addPar ("MOIN",15.0,&MOSFET_B4::Model::moin)
     .setUnit(U_NONE)
     .setCategory(CAT_CAP)
     .setDescription("Coefficient for gate-bias dependent surface potential");

    p.addPar ("NOFF",1.0,&MOSFET_B4::Model::noff)
     .setUnit(U_NONE)
     .setCategory(CAT_CAP)
     .setDescription("C-V turn-on/off parameter");

    p.addPar ("VOFFCV",0.0,&MOSFET_B4::Model::voffcv)
     .setUnit(U_VOLT)
     .setCategory(CAT_CAP)
     .setDescription("C-V lateral-shift parameter");

    p.addPar ("DMCG",0.0,&MOSFET_B4::Model::dmcg)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Distance of Mid-Contact to Gate edge");

    p.addPar ("DMCI",0.0,&MOSFET_B4::Model::dmci)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Distance of Mid-Contact to Isolation");

    p.addPar ("DMDG",0.0,&MOSFET_B4::Model::dmdg)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Distance of Mid-Diffusion to Gate edge");

    p.addPar ("DMCGT",0.0,&MOSFET_B4::Model::dmcgt)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Distance of Mid-Contact to Gate edge in Test structures");

    p.addPar ("XGW",0.0,&MOSFET_B4::Model::xgw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Distance from gate contact center to device edge");

    p.addPar ("XGL",0.0,&MOSFET_B4::Model::xgl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Variation in Ldrawn");

    p.addPar ("RSHG",0.1,&MOSFET_B4::Model::rshg)
     .setUnit(U_OSQM1)
     .setCategory(CAT_PROCESS)
     .setDescription("Gate sheet resistance");

    p.addPar ("NGCON",1.0,&MOSFET_B4::Model::ngcon)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Number of gate contacts");

    p.addPar ("XRCRG1",12.0,&MOSFET_B4::Model::xrcrg1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("First fitting parameter the bias-dependent Rg");

    p.addPar ("XRCRG2",1.0,&MOSFET_B4::Model::xrcrg2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Second fitting parameter the bias-dependent Rg");

    p.addPar ("LAMBDA",0.0,&MOSFET_B4::Model::lambda)
     .setGivenMember(&MOSFET_B4::Model::lambdaGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription(" Velocity overshoot parameter");

    p.addPar ("VTL",2.0e5,&MOSFET_B4::Model::vtl)
     .setGivenMember(&MOSFET_B4::Model::vtlGiven)
     .setUnit(U_MSM1)
     .setCategory(CAT_BASIC)
     .setDescription(" thermal velocity");

    p.addPar ("LC",5.0e-9,&MOSFET_B4::Model::lc)
     .setUnit(U_METER)
     .setCategory(CAT_BASIC)
     .setDescription(" back scattering parameter");

    p.addPar ("XN",3.0,&MOSFET_B4::Model::xn)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setDescription(" back scattering parameter");

    p.addPar ("VFBSDOFF",0.0,&MOSFET_B4::Model::vfbsdoff)
     .setUnit(U_VOLT)
     .setCategory(CAT_TUNNEL)
     .setDescription("S/D flatband voltage offset");

    p.addPar ("TVFBSDOFF",0.0,&MOSFET_B4::Model::tvfbsdoff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature parameter for vfbsdoff");

    p.addPar ("TVOFF",0.0,&MOSFET_B4::Model::tvoff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature parameter for voff");

    p.addPar ("LINTNOI",0.0,&MOSFET_B4::Model::lintnoi)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("lint offset for noise calculation");

    p.addPar ("LINT",0.0,&MOSFET_B4::Model::Lint)
     .setUnit(U_METER)
     .setCategory(CAT_BASIC)
     .setDescription("Length reduction parameter");

    p.addPar ("LL",0.0,&MOSFET_B4::Model::Ll)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length reduction parameter");

    p.addPar ("LLC",0.0,&MOSFET_B4::Model::Llc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length reduction parameter for CV");

    p.addPar ("LLN",1.0,&MOSFET_B4::Model::Lln)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length reduction parameter");

    p.addPar ("LW",0.0,&MOSFET_B4::Model::Lw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length reduction parameter");

    p.addPar ("LWC",0.0,&MOSFET_B4::Model::Lwc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length reduction parameter for CV");

    p.addPar ("LWN",1.0,&MOSFET_B4::Model::Lwn)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length reduction parameter");

    p.addPar ("LWL",0.0,&MOSFET_B4::Model::Lwl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length reduction parameter");

    p.addPar ("LWLC",0.0,&MOSFET_B4::Model::Lwlc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length reduction parameter for CV");

    p.addPar ("LMIN",0.0,&MOSFET_B4::Model::Lmin)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Minimum length for the model");

    p.addPar ("LMAX",1.0,&MOSFET_B4::Model::Lmax)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Maximum length for the model");

    p.addPar ("WR",1.0,&MOSFET_B4::Model::wr)
     .setUnit(U_NONE)
     .setCategory(CAT_ASYMRDS)
     .setDescription("Width dependence of rds");

    p.addPar ("WINT",0.0,&MOSFET_B4::Model::Wint)
     .setUnit(U_METER)
     .setCategory(CAT_BASIC)
     .setDescription("Width reduction parameter");

    p.addPar ("DWG",0.0,&MOSFET_B4::Model::dwg)
     .setUnit(U_MVM1)
     .setCategory(CAT_BASIC)
     .setDescription("Width reduction parameter");

    p.addPar ("DWB",0.0,&MOSFET_B4::Model::dwb)
     .setUnit(U_MVMH)
     .setCategory(CAT_BASIC)
     .setDescription("Width reduction parameter");

    p.addPar ("WL",0.0,&MOSFET_B4::Model::Wl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width reduction parameter");

    p.addPar ("WLC",0.0,&MOSFET_B4::Model::Wlc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width reduction parameter for CV");

    p.addPar ("WLN",1.0,&MOSFET_B4::Model::Wln)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width reduction parameter");

    p.addPar ("WW",0.0,&MOSFET_B4::Model::Ww)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width reduction parameter");

    p.addPar ("WWC",0.0,&MOSFET_B4::Model::Wwc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width reduction parameter for CV");

    p.addPar ("WWN",1.0,&MOSFET_B4::Model::Wwn)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width reduction parameter");

    p.addPar ("WWL",0.0,&MOSFET_B4::Model::Wwl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width reduction parameter");

    p.addPar ("WWLC",0.0,&MOSFET_B4::Model::Wwlc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width reduction parameter for CV");

    p.addPar ("WMIN",0.0,&MOSFET_B4::Model::Wmin)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Minimum width for the model");

    p.addPar ("WMAX",1.0,&MOSFET_B4::Model::Wmax)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Maximum width for the model");

    p.addPar ("B0",0.0,&MOSFET_B4::Model::b0)
     .setUnit(U_METER)
     .setCategory(CAT_BASIC)
     .setDescription("Abulk narrow width parameter");

    p.addPar ("B1",0.0,&MOSFET_B4::Model::b1)
     .setUnit(U_METER)
     .setCategory(CAT_BASIC)
     .setDescription("Abulk narrow width parameter");

    p.addPar ("CGSL",0.0,&MOSFET_B4::Model::cgsl)
     .setUnit(U_FARADMM1)
     .setCategory(CAT_CAP)
     .setDescription("New C-V model parameter");

    p.addPar ("CGDL",0.0,&MOSFET_B4::Model::cgdl)
     .setUnit(U_FARADMM1)
     .setCategory(CAT_CAP)
     .setDescription("New C-V model parameter");

    p.addPar ("CKAPPAS",0.6,&MOSFET_B4::Model::ckappas)
     .setUnit(U_VOLT)
     .setCategory(CAT_CAP)
     .setDescription("S/G overlap C-V parameter ");

    p.addPar ("CKAPPAD",0.6,&MOSFET_B4::Model::ckappad)
     .setUnit(U_VOLT)
     .setCategory(CAT_CAP)
     .setDescription("D/G overlap C-V parameter");

    p.addPar ("CF",0.0,&MOSFET_B4::Model::cf)
     .setUnit(U_FARADMM1)
     .setCategory(CAT_CAP)
     .setDescription("Fringe capacitance parameter");

    p.addPar ("CLC",0.1e-6,&MOSFET_B4::Model::clc)
     .setUnit(U_METER)
     .setCategory(CAT_CAP)
     .setDescription("Vdsat parameter for C-V model");

    p.addPar ("CLE",0.6,&MOSFET_B4::Model::cle)
     .setUnit(U_NONE)
     .setCategory(CAT_CAP)
     .setDescription("Vdsat parameter for C-V model");

    p.addPar ("DWC",0.0,&MOSFET_B4::Model::dwc)
     .setUnit(U_METER)
     .setCategory(CAT_CAP)
     .setDescription("Delta W for C-V model");

    p.addPar ("DLC",0.0,&MOSFET_B4::Model::dlc)
     .setGivenMember(&MOSFET_B4::Model::dlcGiven)
     .setUnit(U_METER)
     .setCategory(CAT_CAP)
     .setDescription("Delta L for C-V model");

    p.addPar ("XW",0.0,&MOSFET_B4::Model::xw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("W offset for channel width due to mask/etch effect");

    p.addPar ("XL",0.0,&MOSFET_B4::Model::xl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("L offset for channel length due to mask/etch effect");

    p.addPar ("DLCIG",0.0,&MOSFET_B4::Model::dlcig)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Delta L for Ig model");

    p.addPar ("DLCIGD",0.0,&MOSFET_B4::Model::dlcigd)
     .setUnit(U_METER)
     .setCategory(CAT_TUNNEL)
     .setDescription("Delta L for Ig model drain side");

    p.addPar ("DWJ",0.0,&MOSFET_B4::Model::dwj)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Delta W for S/D junctions");

    p.addPar ("ALPHA0",0.0,&MOSFET_B4::Model::alpha0)
     .setUnit(U_MVM1)
     .setCategory(CAT_IMPACT)
     .setDescription("substrate current model parameter");

    p.addPar ("ALPHA1",0.0,&MOSFET_B4::Model::alpha1)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_IMPACT)
     .setDescription("substrate current model parameter");

    p.addPar ("BETA0",0.0,&MOSFET_B4::Model::beta0)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_IMPACT)
     .setDescription("substrate current model parameter");

    p.addPar ("AGIDL",0.0,&MOSFET_B4::Model::agidl)
     .setUnit(U_OHMM1)
     .setCategory(CAT_GDLEAKAGE)
     .setDescription("Pre-exponential constant for GIDL");

    p.addPar ("BGIDL",2.3e9,&MOSFET_B4::Model::bgidl)
     .setUnit(U_VMM1)
     .setCategory(CAT_GDLEAKAGE)
     .setDescription("Exponential constant for GIDL");

    p.addPar ("CGIDL",0.5,&MOSFET_B4::Model::cgidl)
     .setUnit(U_VOLT3)
     .setCategory(CAT_GDLEAKAGE)
     .setDescription("Parameter for body-bias dependence of GIDL");

    p.addPar ("EGIDL",0.8,&MOSFET_B4::Model::egidl)
     .setUnit(U_VOLT)
     .setCategory(CAT_GDLEAKAGE)
     .setDescription("Fitting parameter for Bandbending");

    p.addPar ("AGISL",0.0,&MOSFET_B4::Model::agisl)
     .setUnit(U_OHMM1)
     .setCategory(CAT_GDLEAKAGE)
     .setDescription("Pre-exponential constant for GISL");

    p.addPar ("BGISL",2.3e-9,&MOSFET_B4::Model::bgisl)
     .setUnit(U_VMM1)
     .setCategory(CAT_GDLEAKAGE)
     .setDescription("Exponential constant for GISL");

    p.addPar ("CGISL",0.5,&MOSFET_B4::Model::cgisl)
     .setUnit(U_VOLT3)
     .setCategory(CAT_GDLEAKAGE)
     .setDescription("Parameter for body-bias dependence of GISL");

    p.addPar ("EGISL",0.8,&MOSFET_B4::Model::egisl)
     .setUnit(U_VOLT)
     .setCategory(CAT_GDLEAKAGE)
     .setDescription("Fitting parameter for Bandbending");

    // These are type-dependent. For the table,assuming NMOS.
    p.addPar ("AIGC",1.36e-2,&MOSFET_B4::Model::aigc)
     .setUnit(U_FS2HGMHMM1)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igc");

    p.addPar ("BIGC",1.71e-3,&MOSFET_B4::Model::bigc)
     .setUnit(U_FS2HGMHMM1VM1)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igc");

    p.addPar ("CIGC",0.075,&MOSFET_B4::Model::cigc)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igc");

    p.addPar ("AIGSD",1.36e-2,&MOSFET_B4::Model::aigsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Parameter for Igs,d");

    p.addPar ("BIGSD",1.71e-3,&MOSFET_B4::Model::bigsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Parameter for Igs,d");

    p.addPar ("CIGSD",0.075,&MOSFET_B4::Model::cigsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Parameter for Igs,d");

    p.addPar ("AIGS",1.36e-2,&MOSFET_B4::Model::aigs)
     .setUnit(U_FS2HGMHMM1)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igs");

    p.addPar ("BIGS",1.71e-3,&MOSFET_B4::Model::bigs)
     .setUnit(U_FS2HGMHMM1VM1)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igs");

    p.addPar ("CIGS",0.075,&MOSFET_B4::Model::cigs)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igs");

    p.addPar ("AIGD",1.36e-2,&MOSFET_B4::Model::aigd)
     .setUnit(U_FS2HGMHMM1)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igd");

    p.addPar ("BIGD",1.71e-3,&MOSFET_B4::Model::bigd)
     .setUnit(U_FS2HGMHMM1VM1)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igd");

    p.addPar ("CIGD",0.075,&MOSFET_B4::Model::cigd)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igd");

    p.addPar ("AIGBACC",1.36e-2,&MOSFET_B4::Model::aigbacc)
     .setUnit(U_FS2HGMHMM1)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igb");

    p.addPar ("BIGBACC",1.71e-3,&MOSFET_B4::Model::bigbacc)
     .setUnit(U_FS2HGMHMM1VM1)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igb");

    p.addPar ("CIGBACC",0.075,&MOSFET_B4::Model::cigbacc)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igb");

    p.addPar ("AIGBINV",1.11e-2,&MOSFET_B4::Model::aigbinv)
     .setUnit(U_FS2HGMHMM1)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igb");

    p.addPar ("BIGBINV",9.49e-4,&MOSFET_B4::Model::bigbinv)
     .setUnit(U_FS2HGMHMM1VM1)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igb");

    p.addPar ("CIGBINV",0.006,&MOSFET_B4::Model::cigbinv)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igb");

    p.addPar ("NIGC",1.0,&MOSFET_B4::Model::nigc)
     .setUnit(U_NONE)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igc slope");

    p.addPar ("NIGBINV",3.0,&MOSFET_B4::Model::nigbinv)
     .setUnit(U_NONE)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igbinv slope");

    p.addPar ("NIGBACC",1.0,&MOSFET_B4::Model::nigbacc)
     .setUnit(U_NONE)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igbacc slope");

    p.addPar ("NTOX",1.0,&MOSFET_B4::Model::ntox)
     .setUnit(U_NONE)
     .setCategory(CAT_TUNNEL)
     .setDescription("Exponent for Tox ratio");

    p.addPar ("EIGBINV",1.1,&MOSFET_B4::Model::eigbinv)
     .setUnit(U_VOLT)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for the Si bandgap for Igbinv");

    p.addPar ("PIGCD",1.0,&MOSFET_B4::Model::pigcd)
     .setGivenMember(&MOSFET_B4::Model::pigcdGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_TUNNEL)
     .setDescription("Parameter for Igc partition");

    p.addPar ("POXEDGE",1.0,&MOSFET_B4::Model::poxedge)
     .setUnit(U_NONE)
     .setCategory(CAT_TUNNEL)
     .setDescription("Factor for the gate edge Tox");

    // By default,the drain values are set to the source values here:
    p.addPar ("IJTHDFWD",0.1,&MOSFET_B4::Model::ijthdfwd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Forward drain diode forward limiting current");

    p.addPar ("IJTHSFWD",0.1,&MOSFET_B4::Model::ijthsfwd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Forward source diode forward limiting current");

    p.addPar ("IJTHDREV",0.1,&MOSFET_B4::Model::ijthdrev)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Reverse drain diode forward limiting current");

    p.addPar ("IJTHSREV",0.1,&MOSFET_B4::Model::ijthsrev)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Reverse source diode forward limiting current");

    p.addPar ("XJBVD",1.0,&MOSFET_B4::Model::xjbvd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Fitting parameter for drain diode breakdown current");

    p.addPar ("XJBVS",1.0,&MOSFET_B4::Model::xjbvs)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Fitting parameter for source diode breakdown current");

    p.addPar ("BVD",10.0,&MOSFET_B4::Model::bvd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain diode breakdown voltage");

    p.addPar ("BVS",10.0,&MOSFET_B4::Model::bvs)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source diode breakdown voltage");

    p.addPar ("JTSS",0.0,&MOSFET_B4::Model::jtss)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source bottom trap-assisted saturation current density");

    p.addPar ("JTSD",0.0,&MOSFET_B4::Model::jtsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain bottom trap-assisted saturation current density");

    p.addPar ("JTSSWS",0.0,&MOSFET_B4::Model::jtssws)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source STI sidewall trap-assisted saturation current density");

    p.addPar ("JTSSWD",0.0,&MOSFET_B4::Model::jtsswd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain STI sidewall trap-assisted saturation current density");

    p.addPar ("JTSSWGS",0.0,&MOSFET_B4::Model::jtsswgs)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source gate-edge sidewall trap-assisted saturation current density");

    p.addPar ("JTSSWGD",0.0,&MOSFET_B4::Model::jtsswgd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain gate-edge sidewall trap-assisted saturation current density");

    p.addPar ("NJTS",20.0,&MOSFET_B4::Model::njts)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Non-ideality factor for bottom junction");

    p.addPar ("NJTSSW",20.0,&MOSFET_B4::Model::njtssw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Non-ideality factor for STI sidewall junction");

    p.addPar ("NJTSSWG",20.0,&MOSFET_B4::Model::njtsswg)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Non-ideality factor for gate-edge sidewall junction");

    p.addPar ("NJTSD",20.0,&MOSFET_B4::Model::njtsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Non-ideality factor for bottom junction drain side");

    p.addPar ("NJTSSWD",20.0,&MOSFET_B4::Model::njtsswd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Non-ideality factor for STI sidewall junction drain side");

    p.addPar ("NJTSSWGD",20.0,&MOSFET_B4::Model::njtsswgd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Non-ideality factor for gate-edge sidewall junction drain side");

    p.addPar ("XTSS",0.02,&MOSFET_B4::Model::xtss)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Power dependence of JTSS on temperature");

    p.addPar ("XTSD",0.02,&MOSFET_B4::Model::xtsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Power dependence of JTSD on temperature");

    p.addPar ("XTSSWS",0.02,&MOSFET_B4::Model::xtssws)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Power dependence of JTSSWS on temperature");

    p.addPar ("XTSSWD",0.02,&MOSFET_B4::Model::xtsswd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Power dependence of JTSSWD on temperature");

    p.addPar ("XTSSWGS",0.02,&MOSFET_B4::Model::xtsswgs)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Power dependence of JTSSWGS on temperature");

    p.addPar ("XTSSWGD",0.02,&MOSFET_B4::Model::xtsswgd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Power dependence of JTSSWGD on temperature");

    p.addPar ("TNJTS",0.0,&MOSFET_B4::Model::tnjts)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient for NJTS");

    p.addPar ("TNJTSSW",0.0,&MOSFET_B4::Model::tnjtssw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient for NJTSSW");

    p.addPar ("TNJTSSWG",0.0,&MOSFET_B4::Model::tnjtsswg)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient for NJTSSWG");

    p.addPar ("TNJTSD",0.0,&MOSFET_B4::Model::tnjtsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient for NJTSD");

    p.addPar ("TNJTSSWD",0.0,&MOSFET_B4::Model::tnjtsswd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient for NJTSSWD");

    p.addPar ("TNJTSSWGD",0.0,&MOSFET_B4::Model::tnjtsswgd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient for NJTSSWGD");

    p.addPar ("VTSS",10.0,&MOSFET_B4::Model::vtss)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source bottom trap-assisted voltage dependent parameter");

    p.addPar ("VTSD",10.0,&MOSFET_B4::Model::vtsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain bottom trap-assisted voltage dependent parameter");

    p.addPar ("VTSSWS",10.0,&MOSFET_B4::Model::vtssws)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source STI sidewall trap-assisted voltage dependent parameter");

    p.addPar ("VTSSWD",10.0,&MOSFET_B4::Model::vtsswd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain STI sidewall trap-assisted voltage dependent parameter");

    p.addPar ("VTSSWGS",10.0,&MOSFET_B4::Model::vtsswgs)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source gate-edge sidewall trap-assisted voltage dependent parameter");

    p.addPar ("VTSSWGD",10.0,&MOSFET_B4::Model::vtsswgd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain gate-edge sidewall trap-assisted voltage dependent parameter");

    p.addPar ("GBMIN",1.0e-12,&MOSFET_B4::Model::gbmin)
     .setUnit(U_OHMM1)
     .setCategory(CAT_NONE)
     .setDescription("Minimum body conductance");

    p.addPar ("RBDB",50.0,&MOSFET_B4::Model::rbdb)
     .setUnit(U_OHM)
     .setCategory(CAT_NONE)
     .setDescription("Resistance between bNode and dbNode");

    p.addPar ("RBPB",50.0,&MOSFET_B4::Model::rbpb)
     .setUnit(U_OHM)
     .setCategory(CAT_NONE)
     .setDescription("Resistance between bNodePrime and bNode");

    p.addPar ("RBSB",50.0,&MOSFET_B4::Model::rbsb)
     .setUnit(U_OHM)
     .setCategory(CAT_NONE)
     .setDescription("Resistance between bNode and sbNode");

    p.addPar ("RBPS",50.0,&MOSFET_B4::Model::rbps)
     .setUnit(U_OHM)
     .setCategory(CAT_NONE)
     .setDescription("Resistance between bNodePrime and sbNode");

    p.addPar ("RBPD",50.0,&MOSFET_B4::Model::rbpd)
     .setUnit(U_OHM)
     .setCategory(CAT_NONE)
     .setDescription("Resistance between bNodePrime and bNode");

    p.addPar ("RBPS0",50.0,&MOSFET_B4::Model::rbps0)
     .setGivenMember(&MOSFET_B4::Model::rbps0Given)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBPS scaling");

    p.addPar ("RBPSL",0.0,&MOSFET_B4::Model::rbpsl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBPS L scaling");

    p.addPar ("RBPSW",0.0,&MOSFET_B4::Model::rbpsw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBPS W scaling");

    p.addPar ("RBPSNF",0.0,&MOSFET_B4::Model::rbpsnf)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBPS NF scaling");

    p.addPar ("RBPD0",50.0,&MOSFET_B4::Model::rbpd0)
     .setGivenMember(&MOSFET_B4::Model::rbpd0Given)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBPD scaling");

    p.addPar ("RBPDL",0.0,&MOSFET_B4::Model::rbpdl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBPD L scaling");

    p.addPar ("RBPDW",0.0,&MOSFET_B4::Model::rbpdw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBPD W scaling");

    p.addPar ("RBPDNF",0.0,&MOSFET_B4::Model::rbpdnf)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBPD NF scaling");

    p.addPar ("RBPBX0",100.0,&MOSFET_B4::Model::rbpbx0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBPBX  scaling");

    p.addPar ("RBPBXL",0.0,&MOSFET_B4::Model::rbpbxl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBPBX L scaling");

    p.addPar ("RBPBXW",0.0,&MOSFET_B4::Model::rbpbxw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBPBX W scaling");

    p.addPar ("RBPBXNF",0.0,&MOSFET_B4::Model::rbpbxnf)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBPBX NF scaling");

    p.addPar ("RBPBY0",100.0,&MOSFET_B4::Model::rbpby0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBPBY  scaling");

    p.addPar ("RBPBYL",0.0,&MOSFET_B4::Model::rbpbyl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBPBY L scaling");

    p.addPar ("RBPBYW",0.0,&MOSFET_B4::Model::rbpbyw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBPBY W scaling");

    p.addPar ("RBPBYNF",0.0,&MOSFET_B4::Model::rbpbynf)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBPBY NF scaling");

    p.addPar ("RBSBX0",100.0,&MOSFET_B4::Model::rbsbx0)
     .setGivenMember(&MOSFET_B4::Model::rbsbx0Given)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBSBX  scaling");

    p.addPar ("RBSBY0",100.0,&MOSFET_B4::Model::rbsby0)
     .setGivenMember(&MOSFET_B4::Model::rbsby0Given)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBSBY  scaling");

    p.addPar ("RBDBX0",100.0,&MOSFET_B4::Model::rbdbx0)
     .setGivenMember(&MOSFET_B4::Model::rbdbx0Given)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBDBX  scaling");

    p.addPar ("RBDBY0",100.0,&MOSFET_B4::Model::rbdby0)
     .setGivenMember(&MOSFET_B4::Model::rbdby0Given)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBDBY  scaling");

    p.addPar ("RBSDBXL",0.0,&MOSFET_B4::Model::rbsdbxl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBSDBX L scaling");

    p.addPar ("RBSDBXW",0.0,&MOSFET_B4::Model::rbsdbxw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBSDBX W scaling");

    p.addPar ("RBSDBXNF",0.0,&MOSFET_B4::Model::rbsdbxnf)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBSDBX NF scaling");

    p.addPar ("RBSDBYL",0.0,&MOSFET_B4::Model::rbsdbyl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBSDBY L scaling");

    p.addPar ("RBSDBYW",0.0,&MOSFET_B4::Model::rbsdbyw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBSDBY W scaling");

    p.addPar ("RBSDBYNF",0.0,&MOSFET_B4::Model::rbsdbynf)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Body resistance RBSDBY NF scaling");

    p.addPar ("LCDSC",0.0,&MOSFET_B4::Model::lcdsc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of cdsc");

    p.addPar ("LCDSCB",0.0,&MOSFET_B4::Model::lcdscb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of cdscb");

    p.addPar ("LCDSCD",0.0,&MOSFET_B4::Model::lcdscd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of cdscd");

    p.addPar ("LCIT",0.0,&MOSFET_B4::Model::lcit)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of cit");

    p.addPar ("LNFACTOR",0.0,&MOSFET_B4::Model::lnfactor)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of nfactor");

    p.addPar ("LXJ",0.0,&MOSFET_B4::Model::lxj)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of xj");

    p.addPar ("LVSAT",0.0,&MOSFET_B4::Model::lvsat)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of vsat");

    p.addPar ("LAT",0.0,&MOSFET_B4::Model::lat)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of at");

    p.addPar ("LA0",0.0,&MOSFET_B4::Model::la0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of a0");

    p.addPar ("LAGS",0.0,&MOSFET_B4::Model::lags)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of ags");

    p.addPar ("LA1",0.0,&MOSFET_B4::Model::la1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of a1");

    p.addPar ("LA2",0.0,&MOSFET_B4::Model::la2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of a2");

    p.addPar ("LKETA",0.0,&MOSFET_B4::Model::lketa)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of keta");

    p.addPar ("LNSUB",0.0,&MOSFET_B4::Model::lnsub)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of nsub");

    p.addPar ("LNDEP",0.0,&MOSFET_B4::Model::lndep)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of ndep");

    p.addPar ("LNSD",0.0,&MOSFET_B4::Model::lnsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of nsd");

    p.addPar ("LPHIN",0.0,&MOSFET_B4::Model::lphin)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of phin");

    p.addPar ("LNGATE",0.0,&MOSFET_B4::Model::lngate)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of ngate");

    p.addPar ("LGAMMA1",0.0,&MOSFET_B4::Model::lgamma1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of gamma1");

    p.addPar ("LGAMMA2",0.0,&MOSFET_B4::Model::lgamma2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of gamma2");

    p.addPar ("LVBX",0.0,&MOSFET_B4::Model::lvbx)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of vbx");

    p.addPar ("LVBM",0.0,&MOSFET_B4::Model::lvbm)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of vbm");

    p.addPar ("LXT",0.0,&MOSFET_B4::Model::lxt)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of xt");

    p.addPar ("LK1",0.0,&MOSFET_B4::Model::lk1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of k1");

    p.addPar ("LKT1",0.0,&MOSFET_B4::Model::lkt1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of kt1");

    p.addPar ("LKT1L",0.0,&MOSFET_B4::Model::lkt1l)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of kt1l");

    p.addPar ("LKT2",0.0,&MOSFET_B4::Model::lkt2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of kt2");

    p.addPar ("LK2",0.0,&MOSFET_B4::Model::lk2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of k2");

    p.addPar ("LK3",0.0,&MOSFET_B4::Model::lk3)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of k3");

    p.addPar ("LK3B",0.0,&MOSFET_B4::Model::lk3b)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of k3b");

    p.addPar ("LW0",0.0,&MOSFET_B4::Model::lw0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of w0");

    p.addPar ("LDVTP0",0.0,&MOSFET_B4::Model::ldvtp0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of dvtp0");

    p.addPar ("LDVTP1",0.0,&MOSFET_B4::Model::ldvtp1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of dvtp1");

    p.addPar ("LLPE0",0.0,&MOSFET_B4::Model::llpe0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of lpe0");

    p.addPar ("LLPEB",0.0,&MOSFET_B4::Model::llpeb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of lpeb");

    p.addPar ("LDVT0",0.0,&MOSFET_B4::Model::ldvt0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of dvt0");

    p.addPar ("LDVT1",0.0,&MOSFET_B4::Model::ldvt1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of dvt1");

    p.addPar ("LDVT2",0.0,&MOSFET_B4::Model::ldvt2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of dvt2");

    p.addPar ("LDVT0W",0.0,&MOSFET_B4::Model::ldvt0w)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of dvt0w");

    p.addPar ("LDVT1W",0.0,&MOSFET_B4::Model::ldvt1w)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of dvt1w");

    p.addPar ("LDVT2W",0.0,&MOSFET_B4::Model::ldvt2w)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of dvt2w");

    p.addPar ("LDROUT",0.0,&MOSFET_B4::Model::ldrout)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of drout");

    p.addPar ("LDSUB",0.0,&MOSFET_B4::Model::ldsub)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of dsub");

    p.addPar ("LVTH0",0.0,&MOSFET_B4::Model::lvth0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("");

    p.addPar ("LUA",0.0,&MOSFET_B4::Model::lua)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of ua");

    p.addPar ("LUA1",0.0,&MOSFET_B4::Model::lua1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of ua1");

    p.addPar ("LUB",0.0,&MOSFET_B4::Model::lub)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of ub");

    p.addPar ("LUB1",0.0,&MOSFET_B4::Model::lub1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of ub1");

    p.addPar ("LUC",0.0,&MOSFET_B4::Model::luc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of uc");

    p.addPar ("LUC1",0.0,&MOSFET_B4::Model::luc1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of uc1");

    p.addPar ("LUD",0.0,&MOSFET_B4::Model::lud)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of ud");

    p.addPar ("LUD1",0.0,&MOSFET_B4::Model::lud1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of ud1");

    p.addPar ("LUP",0.0,&MOSFET_B4::Model::lup)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of up");

    p.addPar ("LLP",0.0,&MOSFET_B4::Model::llp)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of lp");

    p.addPar ("LU0",0.0,&MOSFET_B4::Model::lu0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of u0");

    p.addPar ("LUTE",0.0,&MOSFET_B4::Model::lute)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of ute");

    p.addPar ("LVOFF",0.0,&MOSFET_B4::Model::lvoff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of voff");

    p.addPar ("LMINV",0.0,&MOSFET_B4::Model::lminv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of minv");

    p.addPar ("LMINVCV",0.0,&MOSFET_B4::Model::lminvcv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of minvcv");

    p.addPar ("LDELTA",0.0,&MOSFET_B4::Model::ldelta)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of delta");

    p.addPar ("LRDSW",0.0,&MOSFET_B4::Model::lrdsw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of rdsw ");

    p.addPar ("LRSW",0.0,&MOSFET_B4::Model::lrsw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of rsw");

    p.addPar ("LRDW",0.0,&MOSFET_B4::Model::lrdw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of rdw");

    p.addPar ("LPRWG",0.0,&MOSFET_B4::Model::lprwg)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of prwg ");

    p.addPar ("LPRWB",0.0,&MOSFET_B4::Model::lprwb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of prwb ");

    p.addPar ("LPRT",0.0,&MOSFET_B4::Model::lprt)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of prt ");

    p.addPar ("LETA0",0.0,&MOSFET_B4::Model::leta0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of eta0");

    p.addPar ("LETAB",0.0,&MOSFET_B4::Model::letab)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of etab");

    p.addPar ("LPCLM",0.0,&MOSFET_B4::Model::lpclm)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of pclm");

    p.addPar ("LPDIBLC1",0.0,&MOSFET_B4::Model::lpdibl1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of pdiblc1");

    p.addPar ("LPDIBLC2",0.0,&MOSFET_B4::Model::lpdibl2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of pdiblc2");

    p.addPar ("LPDIBLCB",0.0,&MOSFET_B4::Model::lpdiblb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of pdiblcb");

    p.addPar ("LFPROUT",0.0,&MOSFET_B4::Model::lfprout)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of pdiblcb");

    p.addPar ("LPDITS",0.0,&MOSFET_B4::Model::lpdits)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of pdits");

    p.addPar ("LPDITSD",0.0,&MOSFET_B4::Model::lpditsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of pditsd");

    p.addPar ("LPSCBE1",0.0,&MOSFET_B4::Model::lpscbe1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of pscbe1");

    p.addPar ("LPSCBE2",0.0,&MOSFET_B4::Model::lpscbe2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of pscbe2");

    p.addPar ("LPVAG",0.0,&MOSFET_B4::Model::lpvag)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of pvag");

    p.addPar ("LWR",0.0,&MOSFET_B4::Model::lwr)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of wr");

    p.addPar ("LDWG",0.0,&MOSFET_B4::Model::ldwg)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of dwg");

    p.addPar ("LDWB",0.0,&MOSFET_B4::Model::ldwb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of dwb");

    p.addPar ("LB0",0.0,&MOSFET_B4::Model::lb0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of b0");

    p.addPar ("LB1",0.0,&MOSFET_B4::Model::lb1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of b1");

    p.addPar ("LCGSL",0.0,&MOSFET_B4::Model::lcgsl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of cgsl");

    p.addPar ("LCGDL",0.0,&MOSFET_B4::Model::lcgdl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of cgdl");

    p.addPar ("LCKAPPAS",0.0,&MOSFET_B4::Model::lckappas)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of ckappas");

    p.addPar ("LCKAPPAD",0.0,&MOSFET_B4::Model::lckappad)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of ckappad");

    p.addPar ("LCF",0.0,&MOSFET_B4::Model::lcf)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of cf");

    p.addPar ("LCLC",0.0,&MOSFET_B4::Model::lclc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of clc");

    p.addPar ("LCLE",0.0,&MOSFET_B4::Model::lcle)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of cle");

    p.addPar ("LALPHA0",0.0,&MOSFET_B4::Model::lalpha0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of alpha0");

    p.addPar ("LALPHA1",0.0,&MOSFET_B4::Model::lalpha1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of alpha1");

    p.addPar ("LBETA0",0.0,&MOSFET_B4::Model::lbeta0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of beta0");

    p.addPar ("LAGIDL",0.0,&MOSFET_B4::Model::lagidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of agidl");

    p.addPar ("LBGIDL",0.0,&MOSFET_B4::Model::lbgidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of bgidl");

    p.addPar ("LCGIDL",0.0,&MOSFET_B4::Model::lcgidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of cgidl");

    p.addPar ("LEGIDL",0.0,&MOSFET_B4::Model::legidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of egidl");

    p.addPar ("LAGISL",0.0,&MOSFET_B4::Model::lagisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of agisl");

    p.addPar ("LBGISL",0.0,&MOSFET_B4::Model::lbgisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of bgisl");

    p.addPar ("LCGISL",0.0,&MOSFET_B4::Model::lcgisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of cgisl");

    p.addPar ("LEGISL",0.0,&MOSFET_B4::Model::legisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of egisl");

    p.addPar ("LAIGC",0.0,&MOSFET_B4::Model::laigc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of aigc");

    p.addPar ("LBIGC",0.0,&MOSFET_B4::Model::lbigc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of bigc");

    p.addPar ("LCIGC",0.0,&MOSFET_B4::Model::lcigc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of cigc");

    p.addPar ("LAIGSD",0.0,&MOSFET_B4::Model::laigsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of aigsd");

    p.addPar ("LBIGSD",0.0,&MOSFET_B4::Model::lbigsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of bigsd");

    p.addPar ("LCIGSD",0.0,&MOSFET_B4::Model::lcigsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of cigsd");


    p.addPar ("LAIGS",0.0,&MOSFET_B4::Model::laigs)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of aigs");

    p.addPar ("LBIGS",0.0,&MOSFET_B4::Model::lbigs)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of bigs");

    p.addPar ("LCIGS",0.0,&MOSFET_B4::Model::lcigs)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of cigs");

    p.addPar ("LAIGD",0.0,&MOSFET_B4::Model::laigd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of aigd");

    p.addPar ("LBIGD",0.0,&MOSFET_B4::Model::lbigd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of bigd");

    p.addPar ("LCIGD",0.0,&MOSFET_B4::Model::lcigd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of cigd");

    p.addPar ("LAIGBACC",0.0,&MOSFET_B4::Model::laigbacc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of aigbacc");

    p.addPar ("LBIGBACC",0.0,&MOSFET_B4::Model::lbigbacc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of bigbacc");

    p.addPar ("LCIGBACC",0.0,&MOSFET_B4::Model::lcigbacc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of cigbacc");

    p.addPar ("LAIGBINV",0.0,&MOSFET_B4::Model::laigbinv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of aigbinv");

    p.addPar ("LBIGBINV",0.0,&MOSFET_B4::Model::lbigbinv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of bigbinv");

    p.addPar ("LCIGBINV",0.0,&MOSFET_B4::Model::lcigbinv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of cigbinv");

    p.addPar ("LNIGC",0.0,&MOSFET_B4::Model::lnigc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of nigc");

    p.addPar ("LNIGBINV",0.0,&MOSFET_B4::Model::lnigbinv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of nigbinv");

    p.addPar ("LNIGBACC",0.0,&MOSFET_B4::Model::lnigbacc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of nigbacc");

    p.addPar ("LNTOX",0.0,&MOSFET_B4::Model::lntox)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of ntox");

    p.addPar ("LEIGBINV",0.0,&MOSFET_B4::Model::leigbinv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence for eigbinv");

    p.addPar ("LPIGCD",0.0,&MOSFET_B4::Model::lpigcd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence for pigcd");

    p.addPar ("LPOXEDGE",0.0,&MOSFET_B4::Model::lpoxedge)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence for poxedge");

    p.addPar ("LVFBCV",0.0,&MOSFET_B4::Model::lvfbcv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of vfbcv");

    p.addPar ("LVFB",0.0,&MOSFET_B4::Model::lvfb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of vfb");

    p.addPar ("LACDE",0.0,&MOSFET_B4::Model::lacde)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of acde");

    p.addPar ("LMOIN",0.0,&MOSFET_B4::Model::lmoin)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of moin");

    p.addPar ("LNOFF",0.0,&MOSFET_B4::Model::lnoff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of noff");

    p.addPar ("LVOFFCV",0.0,&MOSFET_B4::Model::lvoffcv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of voffcv");

    p.addPar ("LXRCRG1",0.0,&MOSFET_B4::Model::lxrcrg1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of xrcrg1");

    p.addPar ("LXRCRG2",0.0,&MOSFET_B4::Model::lxrcrg2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of xrcrg2");

    p.addPar ("LLAMBDA",0.0,&MOSFET_B4::Model::llambda)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of lambda");

    p.addPar ("LVTL",0.0,&MOSFET_B4::Model::lvtl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription(" Length dependence of vtl");

    p.addPar ("LXN",0.0,&MOSFET_B4::Model::lxn)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription(" Length dependence of xn");

    p.addPar ("LEU",0.0,&MOSFET_B4::Model::leu)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription(" Length dependence of eu");

    p.addPar ("LVFBSDOFF",0.0,&MOSFET_B4::Model::lvfbsdoff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of vfbsdoff");

    p.addPar ("LTVFBSDOFF",0.0,&MOSFET_B4::Model::ltvfbsdoff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of tvfbsdoff");

    p.addPar ("LTVOFF",0.0,&MOSFET_B4::Model::ltvoff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of tvoff");

    p.addPar ("WCDSC",0.0,&MOSFET_B4::Model::wcdsc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of cdsc");

    p.addPar ("WCDSCB",0.0,&MOSFET_B4::Model::wcdscb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of cdscb");

    p.addPar ("WCDSCD",0.0,&MOSFET_B4::Model::wcdscd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of cdscd");

    p.addPar ("WCIT",0.0,&MOSFET_B4::Model::wcit)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of cit");

    p.addPar ("WNFACTOR",0.0,&MOSFET_B4::Model::wnfactor)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of nfactor");

    p.addPar ("WXJ",0.0,&MOSFET_B4::Model::wxj)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of xj");

    p.addPar ("WVSAT",0.0,&MOSFET_B4::Model::wvsat)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of vsat");

    p.addPar ("WAT",0.0,&MOSFET_B4::Model::wat)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of at");

    p.addPar ("WA0",0.0,&MOSFET_B4::Model::wa0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of a0");

    p.addPar ("WAGS",0.0,&MOSFET_B4::Model::wags)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of ags");

    p.addPar ("WA1",0.0,&MOSFET_B4::Model::wa1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of a1");

    p.addPar ("WA2",0.0,&MOSFET_B4::Model::wa2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of a2");

    p.addPar ("WKETA",0.0,&MOSFET_B4::Model::wketa)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of keta");

    p.addPar ("WNSUB",0.0,&MOSFET_B4::Model::wnsub)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of nsub");

    p.addPar ("WNDEP",0.0,&MOSFET_B4::Model::wndep)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of ndep");

    p.addPar ("WNSD",0.0,&MOSFET_B4::Model::wnsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of nsd");

    p.addPar ("WPHIN",0.0,&MOSFET_B4::Model::wphin)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of phin");

    p.addPar ("WNGATE",0.0,&MOSFET_B4::Model::wngate)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of ngate");

    p.addPar ("WGAMMA1",0.0,&MOSFET_B4::Model::wgamma1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of gamma1");

    p.addPar ("WGAMMA2",0.0,&MOSFET_B4::Model::wgamma2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of gamma2");

    p.addPar ("WVBX",0.0,&MOSFET_B4::Model::wvbx)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of vbx");

    p.addPar ("WVBM",0.0,&MOSFET_B4::Model::wvbm)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of vbm");

    p.addPar ("WXT",0.0,&MOSFET_B4::Model::wxt)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of xt");

    p.addPar ("WK1",0.0,&MOSFET_B4::Model::wk1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of k1");

    p.addPar ("WKT1",0.0,&MOSFET_B4::Model::wkt1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of kt1");

    p.addPar ("WKT1L",0.0,&MOSFET_B4::Model::wkt1l)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of kt1l");

    p.addPar ("WKT2",0.0,&MOSFET_B4::Model::wkt2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of kt2");

    p.addPar ("WK2",0.0,&MOSFET_B4::Model::wk2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of k2");

    p.addPar ("WK3",0.0,&MOSFET_B4::Model::wk3)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of k3");

    p.addPar ("WK3B",0.0,&MOSFET_B4::Model::wk3b)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of k3b");

    p.addPar ("WW0",0.0,&MOSFET_B4::Model::ww0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of w0");

    p.addPar ("WDVTP0",0.0,&MOSFET_B4::Model::wdvtp0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of dvtp0");

    p.addPar ("WDVTP1",0.0,&MOSFET_B4::Model::wdvtp1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of dvtp1");

    p.addPar ("WLPE0",0.0,&MOSFET_B4::Model::wlpe0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of lpe0");

    p.addPar ("WLPEB",0.0,&MOSFET_B4::Model::wlpeb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of lpeb");

    p.addPar ("WDVT0",0.0,&MOSFET_B4::Model::wdvt0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of dvt0");

    p.addPar ("WDVT1",0.0,&MOSFET_B4::Model::wdvt1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of dvt1");

    p.addPar ("WDVT2",0.0,&MOSFET_B4::Model::wdvt2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of dvt2");

    p.addPar ("WDVT0W",0.0,&MOSFET_B4::Model::wdvt0w)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of dvt0w");

    p.addPar ("WDVT1W",0.0,&MOSFET_B4::Model::wdvt1w)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of dvt1w");

    p.addPar ("WDVT2W",0.0,&MOSFET_B4::Model::wdvt2w)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of dvt2w");

    p.addPar ("WDROUT",0.0,&MOSFET_B4::Model::wdrout)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of drout");

    p.addPar ("WDSUB",0.0,&MOSFET_B4::Model::wdsub)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of dsub");

    p.addPar ("WVTH0",0.0,&MOSFET_B4::Model::wvth0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("");

    p.addPar ("WUA",0.0,&MOSFET_B4::Model::wua)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of ua");

    p.addPar ("WUA1",0.0,&MOSFET_B4::Model::wua1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of ua1");

    p.addPar ("WUB",0.0,&MOSFET_B4::Model::wub)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of ub");

    p.addPar ("WUB1",0.0,&MOSFET_B4::Model::wub1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of ub1");

    p.addPar ("WUC",0.0,&MOSFET_B4::Model::wuc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of uc");

    p.addPar ("WUC1",0.0,&MOSFET_B4::Model::wuc1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of uc1");

    p.addPar ("WUD",0.0,&MOSFET_B4::Model::wud)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of ud");

    p.addPar ("WUD1",0.0,&MOSFET_B4::Model::wud1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of ud1");

    p.addPar ("WUP",0.0,&MOSFET_B4::Model::wup)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of up");

    p.addPar ("WLP",0.0,&MOSFET_B4::Model::wlp)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of lp");

    p.addPar ("WU0",0.0,&MOSFET_B4::Model::wu0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of u0");

    p.addPar ("WUTE",0.0,&MOSFET_B4::Model::wute)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of ute");

    p.addPar ("WVOFF",0.0,&MOSFET_B4::Model::wvoff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of voff");

    p.addPar ("WMINV",0.0,&MOSFET_B4::Model::wminv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of minv");

    p.addPar ("WMINVCV",0.0,&MOSFET_B4::Model::wminvcv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of minvcv");

    p.addPar ("WDELTA",0.0,&MOSFET_B4::Model::wdelta)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of delta");

    p.addPar ("WRDSW",0.0,&MOSFET_B4::Model::wrdsw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of rdsw ");

    p.addPar ("WRSW",0.0,&MOSFET_B4::Model::wrsw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of rsw");

    p.addPar ("WRDW",0.0,&MOSFET_B4::Model::wrdw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of rdw");

    p.addPar ("WPRWG",0.0,&MOSFET_B4::Model::wprwg)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of prwg ");

    p.addPar ("WPRWB",0.0,&MOSFET_B4::Model::wprwb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of prwb ");

    p.addPar ("WPRT",0.0,&MOSFET_B4::Model::wprt)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of prt");

    p.addPar ("WETA0",0.0,&MOSFET_B4::Model::weta0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of eta0");

    p.addPar ("WETAB",0.0,&MOSFET_B4::Model::wetab)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of etab");

    p.addPar ("WPCLM",0.0,&MOSFET_B4::Model::wpclm)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of pclm");

    p.addPar ("WPDIBLC1",0.0,&MOSFET_B4::Model::wpdibl1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of pdiblc1");

    p.addPar ("WPDIBLC2",0.0,&MOSFET_B4::Model::wpdibl2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of pdiblc2");

    p.addPar ("WPDIBLCB",0.0,&MOSFET_B4::Model::wpdiblb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of pdiblcb");

    p.addPar ("WFPROUT",0.0,&MOSFET_B4::Model::wfprout)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of pdiblcb");

    p.addPar ("WPDITS",0.0,&MOSFET_B4::Model::wpdits)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of pdits");

    p.addPar ("WPDITSD",0.0,&MOSFET_B4::Model::wpditsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of pditsd");

    p.addPar ("WPSCBE1",0.0,&MOSFET_B4::Model::wpscbe1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of pscbe1");

    p.addPar ("WPSCBE2",0.0,&MOSFET_B4::Model::wpscbe2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of pscbe2");

    p.addPar ("WPVAG",0.0,&MOSFET_B4::Model::wpvag)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of pvag");

    p.addPar ("WWR",0.0,&MOSFET_B4::Model::wwr)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of wr");

    p.addPar ("WDWG",0.0,&MOSFET_B4::Model::wdwg)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of dwg");

    p.addPar ("WDWB",0.0,&MOSFET_B4::Model::wdwb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of dwb");

    p.addPar ("WB0",0.0,&MOSFET_B4::Model::wb0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of b0");

    p.addPar ("WB1",0.0,&MOSFET_B4::Model::wb1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of b1");

    p.addPar ("WCGSL",0.0,&MOSFET_B4::Model::wcgsl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of cgsl");

    p.addPar ("WCGDL",0.0,&MOSFET_B4::Model::wcgdl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of cgdl");

    p.addPar ("WCKAPPAS",0.0,&MOSFET_B4::Model::wckappas)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of ckappas");

    p.addPar ("WCKAPPAD",0.0,&MOSFET_B4::Model::wckappad)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of ckappad");

    p.addPar ("WCF",0.0,&MOSFET_B4::Model::wcf)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of cf");

    p.addPar ("WCLC",0.0,&MOSFET_B4::Model::wclc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of clc");

    p.addPar ("WCLE",0.0,&MOSFET_B4::Model::wcle)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of cle");

    p.addPar ("WALPHA0",0.0,&MOSFET_B4::Model::walpha0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of alpha0");

    p.addPar ("WALPHA1",0.0,&MOSFET_B4::Model::walpha1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of alpha1");

    p.addPar ("WBETA0",0.0,&MOSFET_B4::Model::wbeta0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of beta0");

    p.addPar ("WAGIDL",0.0,&MOSFET_B4::Model::wagidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of agidl");

    p.addPar ("WBGIDL",0.0,&MOSFET_B4::Model::wbgidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of bgidl");

    p.addPar ("WCGIDL",0.0,&MOSFET_B4::Model::wcgidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of cgidl");

    p.addPar ("WEGIDL",0.0,&MOSFET_B4::Model::wegidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of egidl");

    p.addPar ("WAGISL",0.0,&MOSFET_B4::Model::wagisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of agisl");

    p.addPar ("WBGISL",0.0,&MOSFET_B4::Model::wbgisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of bgisl");

    p.addPar ("WCGISL",0.0,&MOSFET_B4::Model::wcgisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of cgisl");

    p.addPar ("WEGISL",0.0,&MOSFET_B4::Model::wegisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of egisl");

    p.addPar ("WAIGC",0.0,&MOSFET_B4::Model::waigc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of aigc");

    p.addPar ("WBIGC",0.0,&MOSFET_B4::Model::wbigc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of bigc");

    p.addPar ("WCIGC",0.0,&MOSFET_B4::Model::wcigc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of cigc");

    p.addPar ("WAIGSD",0.0,&MOSFET_B4::Model::waigsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of aigsd");

    p.addPar ("WBIGSD",0.0,&MOSFET_B4::Model::wbigsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of bigsd");

    p.addPar ("WCIGSD",0.0,&MOSFET_B4::Model::wcigsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of cigsd");

    p.addPar ("WAIGS",0.0,&MOSFET_B4::Model::waigs )
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of aigs");

    p.addPar ("WBIGS",0.0,&MOSFET_B4::Model::wbigs )
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of bigs");

    p.addPar ("WCIGS",0.0,&MOSFET_B4::Model::wcigs )
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of cigs");

    p.addPar ("WAIGD",0.0,&MOSFET_B4::Model::waigd )
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of aigd");

    p.addPar ("WBIGD",0.0,&MOSFET_B4::Model::wbigd )
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of bigd");

    p.addPar ("WCIGD",0.0,&MOSFET_B4::Model::wcigd )
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of cigd");

    p.addPar ("WAIGBACC",0.0,&MOSFET_B4::Model::waigbacc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of aigbacc");

    p.addPar ("WBIGBACC",0.0,&MOSFET_B4::Model::wbigbacc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of bigbacc");

    p.addPar ("WCIGBACC",0.0,&MOSFET_B4::Model::wcigbacc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of cigbacc");

    p.addPar ("WAIGBINV",0.0,&MOSFET_B4::Model::waigbinv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of aigbinv");

    p.addPar ("WBIGBINV",0.0,&MOSFET_B4::Model::wbigbinv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of bigbinv");

    p.addPar ("WCIGBINV",0.0,&MOSFET_B4::Model::wcigbinv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of cigbinv");

    p.addPar ("WNIGC",0.0,&MOSFET_B4::Model::wnigc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of nigc");

    p.addPar ("WNIGBINV",0.0,&MOSFET_B4::Model::wnigbinv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of nigbinv");

    p.addPar ("WNIGBACC",0.0,&MOSFET_B4::Model::wnigbacc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of nigbacc");

    p.addPar ("WNTOX",0.0,&MOSFET_B4::Model::wntox)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of ntox");

    p.addPar ("WEIGBINV",0.0,&MOSFET_B4::Model::weigbinv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence for eigbinv");

    p.addPar ("WPIGCD",0.0,&MOSFET_B4::Model::wpigcd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence for pigcd");

    p.addPar ("WPOXEDGE",0.0,&MOSFET_B4::Model::wpoxedge)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence for poxedge");

    p.addPar ("WVFBCV",0.0,&MOSFET_B4::Model::wvfbcv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of vfbcv");

    p.addPar ("WVFB",0.0,&MOSFET_B4::Model::wvfb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of vfb");

    p.addPar ("WACDE",0.0,&MOSFET_B4::Model::wacde)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of acde");

    p.addPar ("WMOIN",0.0,&MOSFET_B4::Model::wmoin)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of moin");

    p.addPar ("WNOFF",0.0,&MOSFET_B4::Model::wnoff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of noff");

    p.addPar ("WVOFFCV",0.0,&MOSFET_B4::Model::wvoffcv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of voffcv");

    p.addPar ("WXRCRG1",0.0,&MOSFET_B4::Model::wxrcrg1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of xrcrg1");

    p.addPar ("WXRCRG2",0.0,&MOSFET_B4::Model::wxrcrg2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of xrcrg2");

    p.addPar ("WLAMBDA",0.0,&MOSFET_B4::Model::wlambda)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of lambda");

    p.addPar ("WVTL",0.0,&MOSFET_B4::Model::wvtl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of vtl");

    p.addPar ("WXN",0.0,&MOSFET_B4::Model::wxn)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of xn");

    p.addPar ("WEU",0.0,&MOSFET_B4::Model::weu)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of eu");

    p.addPar ("WVFBSDOFF",0.0,&MOSFET_B4::Model::wvfbsdoff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of vfbsdoff");

    p.addPar ("WTVFBSDOFF",0.0,&MOSFET_B4::Model::wtvfbsdoff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of tvfbsdoff");

    p.addPar ("WTVOFF",0.0,&MOSFET_B4::Model::wtvoff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of tvoff");

    p.addPar ("PCDSC",0.0,&MOSFET_B4::Model::pcdsc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of cdsc");

    p.addPar ("PCDSCB",0.0,&MOSFET_B4::Model::pcdscb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of cdscb");

    p.addPar ("PCDSCD",0.0,&MOSFET_B4::Model::pcdscd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of cdscd");

    p.addPar ("PCIT",0.0,&MOSFET_B4::Model::pcit)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of cit");

    p.addPar ("PNFACTOR",0.0,&MOSFET_B4::Model::pnfactor)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of nfactor");

    p.addPar ("PXJ",0.0,&MOSFET_B4::Model::pxj)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of xj");

    p.addPar ("PVSAT",0.0,&MOSFET_B4::Model::pvsat)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of vsat");

    p.addPar ("PAT",0.0,&MOSFET_B4::Model::pat)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of at");

    p.addPar ("PA0",0.0,&MOSFET_B4::Model::pa0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of a0");

    p.addPar ("PAGS",0.0,&MOSFET_B4::Model::pags)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of ags");

    p.addPar ("PA1",0.0,&MOSFET_B4::Model::pa1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of a1");

    p.addPar ("PA2",0.0,&MOSFET_B4::Model::pa2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of a2");

    p.addPar ("PKETA",0.0,&MOSFET_B4::Model::pketa)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of keta");

    p.addPar ("PNSUB",0.0,&MOSFET_B4::Model::pnsub)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of nsub");

    p.addPar ("PNDEP",0.0,&MOSFET_B4::Model::pndep)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of ndep");

    p.addPar ("PNSD",0.0,&MOSFET_B4::Model::pnsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of nsd");

    p.addPar ("PPHIN",0.0,&MOSFET_B4::Model::pphin)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of phin");

    p.addPar ("PNGATE",0.0,&MOSFET_B4::Model::pngate)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of ngate");

    p.addPar ("PGAMMA1",0.0,&MOSFET_B4::Model::pgamma1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of gamma1");

    p.addPar ("PGAMMA2",0.0,&MOSFET_B4::Model::pgamma2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of gamma2");

    p.addPar ("PVBX",0.0,&MOSFET_B4::Model::pvbx)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of vbx");

    p.addPar ("PVBM",0.0,&MOSFET_B4::Model::pvbm)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of vbm");

    p.addPar ("PXT",0.0,&MOSFET_B4::Model::pxt)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of xt");

    p.addPar ("PK1",0.0,&MOSFET_B4::Model::pk1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of k1");

    p.addPar ("PKT1",0.0,&MOSFET_B4::Model::pkt1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of kt1");

    p.addPar ("PKT1L",0.0,&MOSFET_B4::Model::pkt1l)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of kt1l");

    p.addPar ("PKT2",0.0,&MOSFET_B4::Model::pkt2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of kt2");

    p.addPar ("PK2",0.0,&MOSFET_B4::Model::pk2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of k2");

    p.addPar ("PK3",0.0,&MOSFET_B4::Model::pk3)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of k3");

    p.addPar ("PK3B",0.0,&MOSFET_B4::Model::pk3b)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of k3b");

    p.addPar ("PW0",0.0,&MOSFET_B4::Model::pw0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of w0");

    p.addPar ("PDVTP0",0.0,&MOSFET_B4::Model::pdvtp0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of dvtp0");

    p.addPar ("PDVTP1",0.0,&MOSFET_B4::Model::pdvtp1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of dvtp1");

    p.addPar ("PLPE0",0.0,&MOSFET_B4::Model::plpe0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of lpe0");

    p.addPar ("PLPEB",0.0,&MOSFET_B4::Model::plpeb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of lpeb");

    p.addPar ("PDVT0",0.0,&MOSFET_B4::Model::pdvt0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of dvt0");

    p.addPar ("PDVT1",0.0,&MOSFET_B4::Model::pdvt1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of dvt1");

    p.addPar ("PDVT2",0.0,&MOSFET_B4::Model::pdvt2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of dvt2");

    p.addPar ("PDVT0W",0.0,&MOSFET_B4::Model::pdvt0w)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of dvt0w");

    p.addPar ("PDVT1W",0.0,&MOSFET_B4::Model::pdvt1w)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of dvt1w");

    p.addPar ("PDVT2W",0.0,&MOSFET_B4::Model::pdvt2w)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of dvt2w");

    p.addPar ("PDROUT",0.0,&MOSFET_B4::Model::pdrout)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of drout");

    p.addPar ("PDSUB",0.0,&MOSFET_B4::Model::pdsub)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of dsub");

    p.addPar ("PVTH0",0.0,&MOSFET_B4::Model::pvth0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("");

    p.addPar ("PUA",0.0,&MOSFET_B4::Model::pua)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of ua");

    p.addPar ("PUA1",0.0,&MOSFET_B4::Model::pua1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of ua1");

    p.addPar ("PUB",0.0,&MOSFET_B4::Model::pub)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of ub");

    p.addPar ("PUB1",0.0,&MOSFET_B4::Model::pub1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of ub1");

    p.addPar ("PUC",0.0,&MOSFET_B4::Model::puc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of uc");

    p.addPar ("PUC1",0.0,&MOSFET_B4::Model::puc1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of uc1");

    p.addPar ("PUD",0.0,&MOSFET_B4::Model::pud)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of ud");

    p.addPar ("PUD1",0.0,&MOSFET_B4::Model::pud1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of ud1");

    p.addPar ("PUP",0.0,&MOSFET_B4::Model::pup)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of up");

    p.addPar ("PLP",0.0,&MOSFET_B4::Model::plp)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of lp");

    p.addPar ("PU0",0.0,&MOSFET_B4::Model::pu0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of u0");

    p.addPar ("PUTE",0.0,&MOSFET_B4::Model::pute)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of ute");

    p.addPar ("PVOFF",0.0,&MOSFET_B4::Model::pvoff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of voff");

    p.addPar ("PMINV",0.0,&MOSFET_B4::Model::pminv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of minv");

    p.addPar ("PMINVCV",0.0,&MOSFET_B4::Model::pminvcv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of minvcv");

    p.addPar ("PDELTA",0.0,&MOSFET_B4::Model::pdelta)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of delta");

    p.addPar ("PRDSW",0.0,&MOSFET_B4::Model::prdsw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of rdsw ");

    p.addPar ("PRSW",0.0,&MOSFET_B4::Model::prsw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of rsw");

    p.addPar ("PRDW",0.0,&MOSFET_B4::Model::prdw)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of rdw");

    p.addPar ("PPRWG",0.0,&MOSFET_B4::Model::pprwg)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of prwg ");

    p.addPar ("PPRWB",0.0,&MOSFET_B4::Model::pprwb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of prwb ");

    p.addPar ("PPRT",0.0,&MOSFET_B4::Model::pprt)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of prt ");

    p.addPar ("PETA0",0.0,&MOSFET_B4::Model::peta0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of eta0");

    p.addPar ("PETAB",0.0,&MOSFET_B4::Model::petab)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of etab");

    p.addPar ("PPCLM",0.0,&MOSFET_B4::Model::ppclm)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of pclm");

    p.addPar ("PPDIBLC1",0.0,&MOSFET_B4::Model::ppdibl1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of pdiblc1");

    p.addPar ("PPDIBLC2",0.0,&MOSFET_B4::Model::ppdibl2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of pdiblc2");

    p.addPar ("PPDIBLCB",0.0,&MOSFET_B4::Model::ppdiblb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of pdiblcb");

    p.addPar ("PFPROUT",0.0,&MOSFET_B4::Model::pfprout)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of pdiblcb");

    p.addPar ("PPDITS",0.0,&MOSFET_B4::Model::ppdits)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of pdits");

    p.addPar ("PPDITSD",0.0,&MOSFET_B4::Model::ppditsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of pditsd");

    p.addPar ("PPSCBE1",0.0,&MOSFET_B4::Model::ppscbe1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of pscbe1");

    p.addPar ("PPSCBE2",0.0,&MOSFET_B4::Model::ppscbe2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of pscbe2");

    p.addPar ("PPVAG",0.0,&MOSFET_B4::Model::ppvag)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of pvag");

    p.addPar ("PWR",0.0,&MOSFET_B4::Model::pwr)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of wr");

    p.addPar ("PDWG",0.0,&MOSFET_B4::Model::pdwg)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of dwg");

    p.addPar ("PDWB",0.0,&MOSFET_B4::Model::pdwb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of dwb");

    p.addPar ("PB0",0.0,&MOSFET_B4::Model::pb0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of b0");

    p.addPar ("PB1",0.0,&MOSFET_B4::Model::pb1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of b1");

    p.addPar ("PCGSL",0.0,&MOSFET_B4::Model::pcgsl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of cgsl");

    p.addPar ("PCGDL",0.0,&MOSFET_B4::Model::pcgdl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of cgdl");

    p.addPar ("PCKAPPAS",0.0,&MOSFET_B4::Model::pckappas)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of ckappas");

    p.addPar ("PCKAPPAD",0.0,&MOSFET_B4::Model::pckappad)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of ckappad");

    p.addPar ("PCF",0.0,&MOSFET_B4::Model::pcf)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of cf");

    p.addPar ("PCLC",0.0,&MOSFET_B4::Model::pclc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of clc");

    p.addPar ("PCLE",0.0,&MOSFET_B4::Model::pcle)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of cle");

    p.addPar ("PALPHA0",0.0,&MOSFET_B4::Model::palpha0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of alpha0");

    p.addPar ("PALPHA1",0.0,&MOSFET_B4::Model::palpha1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of alpha1");

    p.addPar ("PBETA0",0.0,&MOSFET_B4::Model::pbeta0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of beta0");

    p.addPar ("PAGIDL",0.0,&MOSFET_B4::Model::pagidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of agidl");

    p.addPar ("PBGIDL",0.0,&MOSFET_B4::Model::pbgidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of bgidl");

    p.addPar ("PCGIDL",0.0,&MOSFET_B4::Model::pcgidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of cgidl");

    p.addPar ("PEGIDL",0.0,&MOSFET_B4::Model::pegidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of egidl");

    p.addPar ("PAGISL",0.0,&MOSFET_B4::Model::pagisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of agisl");

    p.addPar ("PBGISL",0.0,&MOSFET_B4::Model::pbgisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of bgisl");

    p.addPar ("PCGISL",0.0,&MOSFET_B4::Model::pcgisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of cgisl");

    p.addPar ("PEGISL",0.0,&MOSFET_B4::Model::pegisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of egisl");

    p.addPar ("PAIGC",0.0,&MOSFET_B4::Model::paigc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of aigc");

    p.addPar ("PBIGC",0.0,&MOSFET_B4::Model::pbigc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of bigc");

    p.addPar ("PCIGC",0.0,&MOSFET_B4::Model::pcigc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of cigc");

    p.addPar ("PAIGSD",0.0,&MOSFET_B4::Model::paigsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of aigsd");

    p.addPar ("PBIGSD",0.0,&MOSFET_B4::Model::pbigsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of bigsd");

    p.addPar ("PCIGSD",0.0,&MOSFET_B4::Model::pcigsd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of cigsd");

    p.addPar ("PAIGS",0.0,&MOSFET_B4::Model::paigs)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of aigs");

    p.addPar ("PBIGS",0.0,&MOSFET_B4::Model::pbigs)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of bigs");

    p.addPar ("PCIGS",0.0,&MOSFET_B4::Model::pcigs)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of cigs");

    p.addPar ("PAIGD",0.0,&MOSFET_B4::Model::paigd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of aigd");

    p.addPar ("PBIGD",0.0,&MOSFET_B4::Model::pbigd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of bigd");

    p.addPar ("PCIGD",0.0,&MOSFET_B4::Model::pcigd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of cigd");

    p.addPar ("PAIGBACC",0.0,&MOSFET_B4::Model::paigbacc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of aigbacc");

    p.addPar ("PBIGBACC",0.0,&MOSFET_B4::Model::pbigbacc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of bigbacc");

    p.addPar ("PCIGBACC",0.0,&MOSFET_B4::Model::pcigbacc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of cigbacc");

    p.addPar ("PAIGBINV",0.0,&MOSFET_B4::Model::paigbinv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of aigbinv");

    p.addPar ("PBIGBINV",0.0,&MOSFET_B4::Model::pbigbinv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of bigbinv");

    p.addPar ("PCIGBINV",0.0,&MOSFET_B4::Model::pcigbinv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of cigbinv");

    p.addPar ("PNIGC",0.0,&MOSFET_B4::Model::pnigc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of nigc");

    p.addPar ("PNIGBINV",0.0,&MOSFET_B4::Model::pnigbinv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of nigbinv");

    p.addPar ("PNIGBACC",0.0,&MOSFET_B4::Model::pnigbacc)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of nigbacc");

    p.addPar ("PNTOX",0.0,&MOSFET_B4::Model::pntox)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of ntox");

    p.addPar ("PEIGBINV",0.0,&MOSFET_B4::Model::peigbinv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence for eigbinv");

    p.addPar ("PPIGCD",0.0,&MOSFET_B4::Model::ppigcd)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence for pigcd");

    p.addPar ("PPOXEDGE",0.0,&MOSFET_B4::Model::ppoxedge)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence for poxedge");

    p.addPar ("PVFBCV",0.0,&MOSFET_B4::Model::pvfbcv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of vfbcv");

    p.addPar ("PVFB",0.0,&MOSFET_B4::Model::pvfb)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of vfb");

    p.addPar ("PACDE",0.0,&MOSFET_B4::Model::pacde)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of acde");

    p.addPar ("PMOIN",0.0,&MOSFET_B4::Model::pmoin)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of moin");

    p.addPar ("PNOFF",0.0,&MOSFET_B4::Model::pnoff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of noff");

    p.addPar ("PVOFFCV",0.0,&MOSFET_B4::Model::pvoffcv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of voffcv");

    p.addPar ("PXRCRG1",0.0,&MOSFET_B4::Model::pxrcrg1)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of xrcrg1");

    p.addPar ("PXRCRG2",0.0,&MOSFET_B4::Model::pxrcrg2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of xrcrg2");

    p.addPar ("PLAMBDA",0.0,&MOSFET_B4::Model::plambda)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of lambda");

    p.addPar ("PVTL",0.0,&MOSFET_B4::Model::pvtl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of vtl");

    p.addPar ("PXN",0.0,&MOSFET_B4::Model::pxn)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of xn");

    p.addPar ("PEU",0.0,&MOSFET_B4::Model::peu)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of eu");

    p.addPar ("PVFBSDOFF",0.0,&MOSFET_B4::Model::pvfbsdoff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of vfbsdoff");

    p.addPar ("PTVFBSDOFF",0.0,&MOSFET_B4::Model::ptvfbsdoff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of tvfbsdoff");

    p.addPar ("PTVOFF",0.0,&MOSFET_B4::Model::ptvoff)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of tvoff");

    // stress effect
    p.addPar ("SAREF",1.0e-6,&MOSFET_B4::Model::saref)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Reference distance between OD edge to poly of one side");

    p.addPar ("SBREF",1.0e-6,&MOSFET_B4::Model::sbref)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Reference distance between OD edge to poly of the other side");

    p.addPar ("WLOD",0.0,&MOSFET_B4::Model::wlod)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width parameter for stress effect");

    p.addPar ("KU0",0.0,&MOSFET_B4::Model::ku0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Mobility degradation/enhancement coefficient for LOD");

    p.addPar ("KVSAT",0.0,&MOSFET_B4::Model::kvsat)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Saturation velocity degradation/enhancement parameter for LOD");

    p.addPar ("KVTH0",0.0,&MOSFET_B4::Model::kvth0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Threshold degradation/enhancement parameter for LOD");

    p.addPar ("TKU0",0.0,&MOSFET_B4::Model::tku0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient of KU0");

    p.addPar ("LLODKU0",0.0,&MOSFET_B4::Model::llodku0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length parameter for u0 LOD effect");

    p.addPar ("WLODKU0",0.0,&MOSFET_B4::Model::wlodku0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width parameter for u0 LOD effect");

    p.addPar ("LLODVTH",0.0,&MOSFET_B4::Model::llodvth)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length parameter for vth LOD effect");

    p.addPar ("WLODVTH",0.0,&MOSFET_B4::Model::wlodvth)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width parameter for vth LOD effect");

    p.addPar ("LKU0",0.0,&MOSFET_B4::Model::lku0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of ku0");

    p.addPar ("WKU0",0.0,&MOSFET_B4::Model::wku0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of ku0");

    p.addPar ("PKU0",0.0,&MOSFET_B4::Model::pku0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of ku0");

    p.addPar ("LKVTH0",0.0,&MOSFET_B4::Model::lkvth0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of kvth0");

    p.addPar ("WKVTH0",0.0,&MOSFET_B4::Model::wkvth0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of kvth0");

    p.addPar ("PKVTH0",0.0,&MOSFET_B4::Model::pkvth0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of kvth0");

    p.addPar ("STK2",0.0,&MOSFET_B4::Model::stk2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("K2 shift factor related to stress effect on vth");

    p.addPar ("LODK2",1.0,&MOSFET_B4::Model::lodk2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("K2 shift modification factor for stress effect");

    p.addPar ("STETA0",0.0,&MOSFET_B4::Model::steta0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("eta0 shift factor related to stress effect on vth");

    p.addPar ("LODETA0",1.0,&MOSFET_B4::Model::lodeta0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("eta0 shift modification factor for stress effect");

    // Well Proximity Effect
    p.addPar ("WEB",0.0,&MOSFET_B4::Model::web)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Coefficient for SCB");

    p.addPar ("WEC",0.0,&MOSFET_B4::Model::wec)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Coefficient for SCC");

    p.addPar ("KVTH0WE",0.0,&MOSFET_B4::Model::kvth0we)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Threshold shift factor for well proximity effect");

    p.addPar ("K2WE",0.0,&MOSFET_B4::Model::k2we)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription(" K2 shift factor for well proximity effect ");

    p.addPar ("KU0WE",0.0,&MOSFET_B4::Model::ku0we)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription(" Mobility degradation factor for well proximity effect ");

    p.addPar ("SCREF",1.0e-6,&MOSFET_B4::Model::scref)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription(" Reference distance to calculate SCA,SCB and SCC");

    p.addPar ("LKVTH0WE",0.0,&MOSFET_B4::Model::lkvth0we)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length dependence of kvth0we");

    p.addPar ("LK2WE",0.0,&MOSFET_B4::Model::lk2we)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription(" Length dependence of k2we ");

    p.addPar ("LKU0WE",0.0,&MOSFET_B4::Model::lku0we)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription(" Length dependence of ku0we ");

    p.addPar ("WKVTH0WE",0.0,&MOSFET_B4::Model::wkvth0we)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width dependence of kvth0we");

    p.addPar ("WK2WE",0.0,&MOSFET_B4::Model::wk2we)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription(" Width dependence of k2we ");

    p.addPar ("WKU0WE",0.0,&MOSFET_B4::Model::wku0we)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription(" Width dependence of ku0we ");

    p.addPar ("PKVTH0WE",0.0,&MOSFET_B4::Model::pkvth0we)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Cross-term dependence of kvth0we");

    p.addPar ("PK2WE",0.0,&MOSFET_B4::Model::pk2we)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription(" Cross-term dependence of k2we ");

    p.addPar ("PKU0WE",0.0,&MOSFET_B4::Model::pku0we)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription(" Cross-term dependence of ku0we ");

    p.addPar ("NOIA",0.0,&MOSFET_B4::Model::oxideTrapDensityA)
     .setUnit(U_NONE)
     .setCategory(CAT_FLICKER)
     .setDescription("Flicker Noise parameter a");

    p.addPar ("NOIB",0.0,&MOSFET_B4::Model::oxideTrapDensityB)
     .setUnit(U_NONE)
     .setCategory(CAT_FLICKER)
     .setDescription("Flicker Noise parameter b");

    p.addPar ("NOIC",0.0,&MOSFET_B4::Model::oxideTrapDensityC)
     .setUnit(U_NONE)
     .setCategory(CAT_FLICKER)
     .setDescription("Flicker Noise parameter c");

    p.addPar ("TNOIA",1.5,&MOSFET_B4::Model::tnoia)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Thermal noise parameter");

    p.addPar ("TNOIB",3.5,&MOSFET_B4::Model::tnoib)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Thermal noise parameter");

    p.addPar ("RNOIA",0.577,&MOSFET_B4::Model::rnoia)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Thermal noise coefficient");

    p.addPar ("RNOIB",0.5164,&MOSFET_B4::Model::rnoib)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Thermal noise coefficient");

    p.addPar ("NTNOI",1.0,&MOSFET_B4::Model::ntnoi)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Thermal noise parameter");

    p.addPar ("EM",4.1e7,&MOSFET_B4::Model::em)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Flicker noise parameter");

    p.addPar ("EF",1.0,&MOSFET_B4::Model::ef)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Flicker noise frequency exponent");

    p.addPar ("AF",1.0,&MOSFET_B4::Model::af)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Flicker noise exponent");

    p.addPar ("KF",0.0,&MOSFET_B4::Model::kf)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Flicker noise coefficient");

    p.addPar ("WPEMOD",0.0,&MOSFET_B4::Model::wpemod)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription(" Flag for WPE model (WPEMOD=1 to activate this model) ");

    // Set up non-double precision variables:
    p.addPar ("CVCHARGEMOD",0,&MOSFET_B4::Model::cvchargeMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Capacitance charge model selector");

    p.addPar ("CAPMOD",2,&MOSFET_B4::Model::capMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Capacitance model selector");

    p.addPar ("DIOMOD",1,&MOSFET_B4::Model::dioMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Diode IV model selector");

    p.addPar ("RDSMOD",0,&MOSFET_B4::Model::rdsMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Bias-dependent S/D resistance model selector");

    p.addPar ("TRNQSMOD",0,&MOSFET_B4::Model::trnqsMod)
     .setGivenMember(&MOSFET_B4::Instance::TRNQSMODgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Transient NQS model selector");

    p.addPar ("ACNQSMOD",0,&MOSFET_B4::Model::acnqsMod)
     .setGivenMember(&MOSFET_B4::Instance::ACNQSMODgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("AC NQS model selector");

    p.addPar ("MOBMOD",0,&MOSFET_B4::Model::mobMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Mobility model selector");

    p.addPar ("RBODYMOD",0,&MOSFET_B4::Model::rbodyMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Distributed body R model selector");

    p.addPar ("RGATEMOD",0,&MOSFET_B4::Model::rgateMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Gate R model selector");

    p.addPar ("PERMOD",1,&MOSFET_B4::Model::perMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Pd and Ps model selector");

    p.addPar ("GEOMOD",0,&MOSFET_B4::Model::geoMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Geometry dependent parasitics model selector");

    p.addPar ("FNOIMOD",1,&MOSFET_B4::Model::fnoiMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Flicker noise model selector");

    p.addPar ("TNOIMOD",0,&MOSFET_B4::Model::tnoiMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Thermal noise model selector");

    p.addPar ("MTRLMOD",0,&MOSFET_B4::Model::mtrlMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("parameter for nonm-silicon substrate or metal gate selector");

    p.addPar ("IGCMOD",0,&MOSFET_B4::Model::igcMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Gate-to-channel Ig model selector");

    p.addPar ("IGBMOD",0,&MOSFET_B4::Model::igbMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Gate-to-body Ig model selector");

    p.addPar ("TEMPMOD",0,&MOSFET_B4::Model::tempMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Temperature model selector");

    p.addPar ("PARAMCHK",1,&MOSFET_B4::Model::paramChk)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Model parameter checking selector");

    p.addPar ("BINUNIT",1,&MOSFET_B4::Model::binUnit)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Bin  unit  selector");

    p.addPar ("VERSION",std::string("4.6.1"),&MOSFET_B4::Model::version)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("parameter for model version");
}

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  double Rtot;

  // apply scale
  if (getDeviceOptions().lengthScale != 1.0)
  {
    if (given("L")) { l *= getDeviceOptions().lengthScale; }
    if (given("W")) { w *= getDeviceOptions().lengthScale; }
    if (sourceAreaGiven) { sourceArea *= getDeviceOptions().lengthScale * getDeviceOptions().lengthScale ; }
    if (drainAreaGiven) { drainArea *= getDeviceOptions().lengthScale * getDeviceOptions().lengthScale ; }
    if (drainPerimeterGiven) { drainPerimeter *= getDeviceOptions().lengthScale; }
    if (sourcePerimeterGiven) { sourcePerimeter *= getDeviceOptions().lengthScale; }
  }

  // Set any non-constant parameter defaults:
  if (!RBDBgiven)
    rbdb = model_.rbdb;
  if (!RBSBgiven)
    rbsb = model_.rbsb;
  if (!RBPBgiven)
    rbpb = model_.rbpb;
  if (!RBPSgiven)
    rbps = model_.rbps;
  if (!RBPDgiven)
    rbpd = model_.rbpd;
  if (!XGWgiven)
    xgw = model_.xgw;
  if (!NGCONgiven)
    ngcon = model_.ngcon;
  if (!SDgiven)
    sd = 2.0 * model_.dmcg;
  if (!TEMPgiven)
    temp = getDeviceOptions().temp.getImmutableValue<double>();
  if (!drainAreaGiven)
    drainArea = getDeviceOptions().defad;
  if (!sourceAreaGiven)
    sourceArea = getDeviceOptions().defas;

  // Process instance (*M_iter) selectors, some
  // may override their global counterparts
  //
  if (!RBODYMODgiven)
  {
    rbodyMod = model_.rbodyMod;
  }
  else if ((rbodyMod != 0) && (rbodyMod != 1) && (rbodyMod != 2))
  {
    rbodyMod = model_.rbodyMod;
    UserWarning(*this) << "rbodyMod has been set to its global value: ";
  }

  if (!RGATEMODgiven)
  {
    rgateMod = model_.rgateMod;
  }
  else if ((rgateMod != 0) && (rgateMod != 1) && (rgateMod != 2) && (rgateMod != 3))
  {
    rgateMod = model_.rgateMod;
    UserWarning(*this) << "rgateMod has been set to its global value: ";
  }

  if (!GEOMODgiven)
  {
    geoMod = model_.geoMod;
  }
  if (!RGEOMODgiven)
  {
    rgeoMod = 0;
  }
  if (!TRNQSMODgiven)
  {
    trnqsMod = model_.trnqsMod;
  }
  else if ((trnqsMod != 0) && (trnqsMod != 1))
  {
    trnqsMod = model_.trnqsMod;
    UserWarning(*this) << "trnqsMod has been set to its global value: ";
  }

  if (!ACNQSMODgiven)
  {
    acnqsMod = model_.acnqsMod;
  }
  else if ((acnqsMod != 0) && (acnqsMod != 1))
  {
    acnqsMod = model_.acnqsMod;
    UserWarning(*this) << "acnqsMod has been set to its global value: ";
  }

  // now set the temperature related stuff.
  updateTemperature(temp);

  // Note:  Xyce does not support noise analysis, so the noiseAnalGive
  // flag is always false.  However, I've left it in here to show consistency
  // with the original spice3f5 code.

  bool noiseAnalGiven=false;

  // process drain series resistance
  int createNode = 0;
  if ( (model_.rdsMod != 0) || (model_.tnoiMod != 0 && noiseAnalGiven))
  {
    createNode = 1;
  }
  else if (model_.sheetResistance > 0)
  {
    if (drainSquaresGiven && drainSquares > 0)
    {
      createNode = 1;
    }
    else if (!drainSquaresGiven && (rgeoMod != 0))
    {
      RdseffGeo(nf, geoMod, rgeoMod, min,
              w, model_.sheetResistance,
              DMCGeff, DMCIeff, DMDGeff, 0, Rtot);

      if(Rtot > 0)
      {
        createNode = 1;
      }
    }
  }

  if ( createNode != 0 )
  {
    drainMOSFET_B4Exists = true;
  }
  else
  {
    drainMOSFET_B4Exists = false;
  }

  // process source series resistance
  createNode = 0;
  if ( (model_.rdsMod != 0) || (model_.tnoiMod != 0 && noiseAnalGiven))
  {
    createNode = 1;
  }
  else if (model_.sheetResistance > 0)
  {
    if (sourceSquaresGiven && sourceSquares > 0)
    {
      createNode = 1;
    }
    else if (!sourceSquaresGiven && (rgeoMod != 0))
    {
      RdseffGeo(nf, geoMod, rgeoMod, min,
              w, model_.sheetResistance,
              DMCGeff, DMCIeff, DMDGeff, 1, Rtot);

      if(Rtot > 0)
      {
        createNode = 1;
      }
    }
  }

  if ( createNode != 0 )
  {
    sourceMOSFET_B4Exists = true;
  }
  else
  {
    sourceMOSFET_B4Exists = false;
  }


  // set up numIntVars:
  numIntVars = 0;

  if (drainMOSFET_B4Exists) ++ numIntVars;
  if (sourceMOSFET_B4Exists) ++ numIntVars;

  if (rgateMod == 1 || rgateMod == 2) ++numIntVars;
  else if (rgateMod == 3) numIntVars+=2;

  if ( trnqsMod ) ++numIntVars;
  if ( rbodyMod ) numIntVars+=3;

  if (icVBSGiven) ++numIntVars;
  if (icVDSGiven) ++numIntVars;
  if (icVGSGiven) ++numIntVars;

  // set up numStateVars
  numStateVars = 3;
  setNumStoreVars(22);

  if (rgateMod == 3)
  {
    numStateVars += 1;
  }

  // parasitic capacitors:
  if (rbodyMod)
  {
    numStateVars += 2;
  }

  if (trnqsMod)
  {
    numStateVars += 2;
  }

  // If there are any time dependent parameters, set their values at for
  // the current time.

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::debugJacStampOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/30/06
//-----------------------------------------------------------------------------
void Instance::debugJacStampOutput ()
{
  Xyce::dout() << "Jacobian stamp:" << std::endl;
  for( int rw=0; rw < jacStamp.size() ; ++rw  )
  {
    Xyce::dout() << "jacStamp[ " << rw << "] = { " ;
    for( int cl=0; cl < jacStamp[rw].size(); ++cl )
    {
      Xyce::dout() << jacStamp[rw][cl];
      if( cl != (jacStamp[rw].size()-1) )
      {
        Xyce::dout() << ", ";
      }
    }
    Xyce::dout() << "}" << std::endl;
  }
  Xyce::dout() << std::endl;

  Xyce::dout() << "And as viewed through the maps"  << std::endl;
  for( int rw=0; rw < jacMap.size() ; ++rw  )
  {
    Xyce::dout() << "jacStamp[ " << rw << "] mapped to jacStamp[ " << jacMap[rw] << "] = { " ;
    for( int cl=0; cl < jacMap2[rw].size(); ++cl )
    {
      Xyce::dout() << jacStamp[jacMap[rw]][jacMap2[rw][cl]];
      if( cl != (jacMap2[rw].size()-1) )
      {
        Xyce::dout() << ", ";
      }
    }
    Xyce::dout() << "}" << std::endl;
  }
  Xyce::dout() << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------

Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & IB,
  Model & Miter,
  const FactoryBlock &  factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    ueff  (0.0),
    thetavth  (0.0),
    von  (0.0),
    vdsat  (0.0),
    cgdo  (0.0),
    qgdo  (0.0),
    cgso  (0.0),
    qgso  (0.0),
    grbsb  (0.0),
    grbdb  (0.0),
    grbpb  (0.0),
    grbps  (0.0),
    grbpd  (0.0),
    vjsmFwd  (0.0),
    vjsmRev  (0.0),
    vjdmFwd  (0.0),
    vjdmRev  (0.0),
    XExpBVS  (0.0),
    XExpBVD  (0.0),
    SslpFwd  (0.0),
    SslpRev  (0.0),
    DslpFwd  (0.0),
    DslpRev  (0.0),
    IVjsmFwd  (0.0),
    IVjsmRev  (0.0),
    IVjdmFwd  (0.0),
    IVjdmRev  (0.0),
    grgeltd  (0.0),
    Pseff  (0.0),
    Pdeff  (0.0),
    Aseff  (0.0),
    Adeff  (0.0),
    l                (getDeviceOptions().defl),
    w                (getDeviceOptions().defw),
    numberParallel   (1.0),
    drainArea        (getDeviceOptions().defad),
    sourceArea       (getDeviceOptions().defas),
    drainSquares(0.0),
    sourceSquares(0.0),
    drainPerimeter(0.0),
    sourcePerimeter(0.0),
    sourceConductance(0.0),
    drainConductance(0.0),
    sa(0.0),
    sb(0.0),
    sd(0.0),
    SDgiven(false),
    sca(0.0),
    scb(0.0),
    scc(0.0),
    sc(0.0),
    rbdb(0.0),
    rbsb(0.0),
    rbpb(0.0),
    rbps(0.0),
    rbpd(0.0),

    RBDBgiven(false),
    RBSBgiven(false),
    RBPBgiven(false),
    RBPSgiven(false),
    RBPDgiven(false),

    delvto(0.0),
    xgw(0.0),
    ngcon(0.0),

    XGWgiven(false),
    NGCONgiven(false),

    u0temp(0.0),
    vsattemp(0.0),
    vth0(0.0),
    vfb(0.0),
    vfbzb(0.0),
    vtfbphi1(0.0),
    vtfbphi2(0.0),
    k2(0.0),
    vbsc(0.0),
    k2ox(0.0),
    eta0(0.0),
    icVDS(0.0),
    icVGS(0.0),
    icVBS(0.0),
    nf(0.0),
    OFF(0),
    mode(0),
    trnqsMod(0),
    acnqsMod(0),
    rbodyMod(0),
    rgateMod(0),
    geoMod(0),
    rgeoMod(0),
    min(0),

    RBODYMODgiven(false),
    RGATEMODgiven(false),
    GEOMODgiven(false),
    RGEOMODgiven(false),
    TRNQSMODgiven(false),
    ACNQSMODgiven(false),

    Vgsteff(0.0),
    Vgsteff_forNoise(0.0),
    vgs_eff(0.0),
    vgd_eff(0.0),
    dvgs_eff_dvg(0.0),
    dvgd_eff_dvg(0.0),
    Vdseff(0.0),
    Vdseff_forNoise(0.0),
    nstar(0.0),
    Abulk(0.0),
    Abulk_forNoise (0.0),
    EsatL(0.0),
    AbovVgst2Vtm(0.0),
    qinv(0.0),
    cd(0.0),
    cbs(0.0),
    cbd(0.0),
    csub(0.0),
    Igidl(0.0),
    Igisl(0.0),
    gm(0.0),
    gds(0.0),
    gmbs(0.0),
    gbd(0.0),
    gbs(0.0),
    gbbs(0.0),
    gbgs(0.0),
    gbds(0.0),
    ggidld(0.0),
    ggidlg(0.0),
    ggidls(0.0),
    ggidlb(0.0),
    ggisld(0.0),
    ggislg(0.0),
    ggisls(0.0),
    ggislb(0.0),
    Igcs(0.0),
    gIgcsg(0.0),
    gIgcsd(0.0),
    gIgcss(0.0),
    gIgcsb(0.0),
    Igcd(0.0),
    gIgcdg(0.0),
    gIgcdd(0.0),
    gIgcds(0.0),
    gIgcdb(0.0),
    Igs(0.0),
    gIgsg(0.0),
    gIgss(0.0),
    Igd(0.0),
    gIgdg(0.0),
    gIgdd(0.0),
    Igb(0.0),
    gIgbg(0.0),
    gIgbd(0.0),
    gIgbs(0.0),
    gIgbb(0.0),
    grdsw(0.0),
    IdovVds(0.0),
    gcrg(0.0),
    gcrgd(0.0),
    gcrgg(0.0),
    gcrgs(0.0),
    gcrgb(0.0),
    gstot(0.0),
    gstotd(0.0),
    gstotg(0.0),
    gstots(0.0),
    gstotb(0.0),
    gdtot(0.0),
    gdtotd(0.0),
    gdtotg(0.0),
    gdtots(0.0),
    gdtotb(0.0),
    cggb(0.0),
    cgdb(0.0),
    cgsb(0.0),
    cbgb(0.0),
    cbdb(0.0),
    cbsb(0.0),
    cdgb(0.0),
    cddb(0.0),
    cdsb(0.0),
    csgb(0.0),
    csdb(0.0),
    cssb(0.0),
    cgbb(0.0),
    cdbb(0.0),
    csbb(0.0),
    cbbb(0.0),
    capbd(0.0),
    capbs(0.0),
    cqgb(0.0),
    cqdb(0.0),
    cqsb(0.0),
    cqbb(0.0),
    qgate(0.0),
    qbulk(0.0),
    qdrn(0.0),
    qsrc(0.0),
    qchqs(0.0),
    taunet(0.0),
    gtau(0.0),
    gtg(0.0),
    gtd(0.0),
    gts(0.0),
    gtb(0.0),
    SjctTempRevSatCur(0.0),
    DjctTempRevSatCur(0.0),
    SswTempRevSatCur(0.0),
    DswTempRevSatCur(0.0),
    SswgTempRevSatCur(0.0),
    DswgTempRevSatCur(0.0),
    limitedFlag(false),
    paramPtr                          (NULL),
    icVBSGiven                        (false),
    icVDSGiven                        (false),
    icVGSGiven                        (false),
    scaGiven(false),
    scbGiven(false),
    sccGiven(false),
    scGiven(false),
    sourcePerimeterGiven(false),
    drainPerimeterGiven(false),
    sourceAreaGiven(false),
    drainAreaGiven(false),
    sourceSquaresGiven(false),
    drainSquaresGiven(false),
    drainMOSFET_B4Exists(false),
    sourceMOSFET_B4Exists(false),
    ChargeComputationNeeded           (true),
    temp                              (getDeviceOptions().temp.getImmutableValue<double>()),
    TEMPgiven(false),
    // missing stuff...
    Vd (0.0),
    Vs (0.0),
    Vb (0.0),
    Vdp(0.0),
    Vsp(0.0),
    Vgp(0.0),
    Vbp(0.0),
    Vge(0.0),
    Vgm(0.0),
    Vdb(0.0),
    Vsb(0.0),
    Vds(0.0),
    Vgs(0.0),
    Vbs(0.0),

    Qtotal                            (0.0),
    Vddp                              (0.0),
    Vssp                              (0.0),
    Vbsp                              (0.0),
    Vbdp                              (0.0),
    Vgsp                              (0.0),
    Vgdp                              (0.0),
    Vgb                               (0.0),
    Vdpsp                             (0.0),
    Vdbb                              (0.0),
    Vdbbp                             (0.0),
    Vbpb                              (0.0),
    Vsbb                              (0.0),
    Vsbbp                             (0.0),
    Idrain                            (0.0),
    Isource                           (0.0),
    Idbb                              (0.0),
    Idbbp                             (0.0),
    Ibpb                              (0.0),
    Isbb                              (0.0),
    Isbbp                             (0.0),
    df1dVdp                           (0.0),
    df2dVdp                           (0.0),
    df1dVsp                           (0.0),
    df2dVsp                           (0.0),
    df1dVg                            (0.0),
    df2dVg                            (0.0),
    df1dVb                            (0.0),
    df2dVb                            (0.0),

    vgb (0.0),
    vgd (0.0),
    cqdef (0.0),
    ceqqd (0.0),
    ceqqb (0.0),
    ceqqg (0.0),

    cqdef_Jdxp (0.0),
    ceqqd_Jdxp (0.0),
    ceqqb_Jdxp (0.0),
    ceqqg_Jdxp (0.0),

    vbd(0.0),
    vbs(0.0),
    vgs(0.0),
    vds(0.0),
    vges(0.0),
    vgms(0.0),
    vdes(0.0),
    vses(0.0),
    vdbs(0.0),
    vsbs(0.0),
    vdbd(0.0),
    vged(0.0),
    vgmd(0.0),

    vbd_old(0.0),
    vbs_old(0.0),
    vgs_old(0.0),
    vds_old(0.0),
    vges_old(0.0),
    vgms_old(0.0),
    vdes_old(0.0),
    vses_old(0.0),
    vdbs_old(0.0),
    vsbs_old(0.0),
    vdbd_old(0.0),
    vged_old(0.0),
    vgmd_old(0.0),

    vbd_orig(0.0),
    vbs_orig(0.0),
    vgs_orig(0.0),
    vds_orig(0.0),
    vgd_orig(0.0),
    vges_orig(0.0),
    vgms_orig(0.0),
    vdes_orig(0.0),
    vses_orig(0.0),
    vdbs_orig(0.0),
    vsbs_orig(0.0),
    vdbd_orig(0.0),
    vbs_jct_orig(0.0),
    vbd_jct_orig(0.0),
    vgmb_orig(0.0),
    vgb_orig(0.0),
    vged_orig(0.0),
    vgmd_orig(0.0),

    Gm(0.0),
    Gmbs(0.0),
    FwdSum(0.0),
    RevSum(0.0),
    ceqdrn(0.0),
    cdrain(0.0),
    ceqbd(0.0),
    ceqbs(0.0),
    cqgate(0.0),
    cqbody(0.0),
    cqdrn(0.0),
    gbbdp(0.0),
    gbbsp(0.0),
    gbdpg(0.0),
    gbdpdp(0.0),
    gbdpb(0.0),
    gbdpsp(0.0),
    gbspg(0.0),
    gbspdp(0.0),
    gbspb(0.0),
    gbspsp(0.0),
    Istoteq(0.0),
    gIstotg(0.0),
    gIstotd(0.0),
    gIstots(0.0),
    gIstotb(0.0),
    Idtoteq(0.0),
    gIdtotg(0.0),
    gIdtotd(0.0),
    gIdtots(0.0),
    gIdtotb(0.0),
    Ibtoteq(0.0),
    gIbtotg(0.0),
    gIbtotd(0.0),
    gIbtots(0.0),
    gIbtotb(0.0),
    Igtoteq(0.0),
    gIgtotg(0.0),
    gIgtotd(0.0),
    gIgtots(0.0),
    gIgtotb(0.0),
    ceqgcrg(0.0),
    ceqgstot(0.0),
    ceqgdtot(0.0),
    ceqjs(0.0),
    ceqjd(0.0),
    vbs_jct(0.0),
    vbd_jct(0.0),

    ceqqjs(0.0),
    ceqqjd(0.0),
    ceqqgmid(0.0),
    gjbd(0.0),
    gjbs(0.0),
    gdpr(0.0),
    gspr(0.0),
    geltd(0.0),
    gcggb(0.0),
    ggtg(0.0),
    gcgdb(0.0),
    ggtd(0.0),
    gcgsb(0.0),
    ggts(0.0),
    gcgbb(0.0),
    ggtb(0.0),
    dxpart(0.0),
    sxpart(0.0),
    gqdef(0.0),
    ddxpart_dVd(0.0),
    ddxpart_dVg(0.0),
    ddxpart_dVb(0.0),
    ddxpart_dVs(0.0),
    dsxpart_dVd(0.0),
    dsxpart_dVg(0.0),
    dsxpart_dVb(0.0),
    dsxpart_dVs(0.0),
    vgmb(0.0),
    CoxWL(0.0),
    ScalingFactor(0.0),
    DMCGeff(0.0),
    DMCIeff(0.0),
    DMDGeff(0.0),

    ceqqgmid_Jdxp(0.0),
    ceqqjs_Jdxp(0.0),
    ceqqjd_Jdxp(0.0),

    qgmb(0.0),
    qgb(0.0),
    Cgg(0.0),
    Cgd(0.0),
    Cgb(0.0),
    Cdg(0.0),
    Cdd(0.0),
    Cds(0.0),
    Csg(0.0),
    Csd(0.0),
    Css(0.0),
    Csb(0.0),
    Cbg(0.0),
    Cbd(0.0),
    Cbb(0.0),


    CAPcggb(0.0),
    CAPcgdb(0.0),
    CAPcgsb(0.0),
    CAPcbgb(0.0),
    CAPcbdb(0.0),
    CAPcbsb(0.0),
    CAPcdgb(0.0),
    CAPcdgmb(0.0),
    CAPcddb(0.0),
    CAPcdbdb(0.0),
    CAPcdsb(0.0),
    CAPcsgb(0.0),
    CAPcsdb(0.0),
    CAPcssb(0.0),
    CAPcgmdb(0.0),
    CAPcgmsb(0.0),
    CAPcgmgmb(0.0),
    CAPcbgmb(0.0),
    CAPcsbsb(0.0),
    CAPcqgb(0.0),
    CAPcqdb(0.0),
    CAPcqsb(0.0),

    CAPcgmbb(0.0),
    CAPcsgmb(0.0),
    CAPcgbb(0.0),
    CAPcdbb(0.0),
    CAPcsbb(0.0),
    CAPcbbb(0.0),
    CAPcqbb(0.0),

    Qeqqd_Jdxp(0.0),
    Qeqqb_Jdxp(0.0),
    Qeqqg_Jdxp(0.0),

    Qeqqgmid_Jdxp(0.0),
    Qeqqjs_Jdxp(0.0),
    Qeqqjd_Jdxp(0.0),
    Qqcheq_Jdxp(0.0),

    Igate(0.0),
    IgateMid(0.0),
    Vgegp(0.0),
    Vgegm(0.0),
    Vgmgp(0.0),

    qb                (0.0),
    qg                (0.0),
    qd                (0.0),
    qgmid             (0.0),
    qbs               (0.0),
    qbd               (0.0),
    qcheq             (0.0),
    cqcheq            (0.0),
    cqcheq_Jdxp       (0.0),
    cqcdump           (0.0),
    qdef(0.0),

// matrix and vectors indices:
// state vector: (local indices)
    li_store_vbd(-1),
    li_store_vbs(-1),
    li_store_vgs(-1),
    li_store_vds(-1),
    li_store_vges(-1),
    li_store_vgms(-1),
    li_store_vdes(-1),
    li_store_vses(-1),
    li_store_vdbs(-1),
    li_store_vsbs(-1),
    li_store_vdbd(-1),
    li_store_vged(-1),
    li_store_vgmd(-1),
    li_store_gm  (-1),
    li_store_Vds     (-1),
    li_store_Vgs     (-1),
    li_store_Vbs     (-1),
    li_store_Vdsat   (-1),
    li_store_Vth     (-1),
    li_store_Gds     (-1),
    li_store_Cgs     (-1),
    li_store_Cgd     (-1),
    Vdsat(0.0),
    Vth(0.0),

    li_state_qb              (-1),
    li_state_qg              (-1),
    li_state_qd              (-1),
    li_state_qgmid           (-1),
    li_state_qbs             (-1),
    li_state_qbd             (-1),
    li_state_qcheq           (-1),
    li_state_qcdump          (-1),
    li_state_qdef            (-1),
    li_branch_dev_id         (-1),
    li_branch_dev_ig         (-1),
    li_branch_dev_is         (-1),
    li_branch_dev_ib         (-1),
// solution vector: (local indices)
    li_Drain                 (-1),
    li_GateExt               (-1),
    li_Source                (-1),
    li_Body                  (-1),
    li_DrainPrime            (-1),
    li_GatePrime             (-1),
    li_GateMid               (-1),
    li_SourcePrime           (-1),
    li_BodyPrime             (-1),
    li_DrainBody             (-1),
    li_SourceBody            (-1),
    li_Charge                (-1),
    li_Ibs                   (-1),
    li_Ids                   (-1),
    li_Igs                   (-1),
    // jacobian matrix offsets:
    GEge(-1),
    GEgp(-1),
    GEdp(-1),
    GEsp(-1),
    GEbp(-1),
    GEgm(-1),
    GEigs(-1),
    GPge(-1),
    GPgp(-1),
    GPdp(-1),
    GPsp(-1),
    GPbp(-1),
    GPq(-1),
    GPgm(-1),
    GMge(-1),
    GMgp(-1),
    GMdp(-1),
    GMsp(-1),
    GMbp(-1),
    GMgm(-1),
    DPgm(-1),
    DPdp(-1),
    DPd(-1),
    DPgp(-1),
    DPsp(-1),
    DPbp(-1),
    DPdb(-1),
    DPq(-1),
    Dd(-1),
    Dgp(-1),
    Ddp(-1),
    Dsp(-1),
    Dbp(-1),
    Dids(-1),
    SPgm(-1),
    SPdp(-1),
    SPgp(-1),
    SPsp(-1),
    SPs(-1),
    SPbp(-1),
    SPsb(-1),
    SPq(-1),
    Ss(-1),
    Sdp(-1),
    Sgp(-1),
    Ssp(-1),
    Sbp(-1),
    Sibs(-1),
    Sids(-1),
    Sigs(-1),
    BPgm(-1),
    BPdp(-1),
    BPgp(-1),
    BPsp(-1),
    BPb(-1),
    BPbp(-1),
    BPdb(-1),
    BPsb(-1),
    DBdp(-1),
    DBdb(-1),
    DBbp(-1),
    DBb(-1),
    SBsp(-1),
    SBbp(-1),
    SBb(-1),
    SBsb(-1),
    Bdb(-1),
    Bbp(-1),
    Bsb(-1),
    Bb(-1),
    Bibs(-1),
    Qq(-1),
    Qgp(-1),
    Qdp(-1),
    Qsp(-1),
    Qbp(-1),

    IBSb(-1),
    IBSs(-1),
    IBSibs(-1),

    IDSd(-1),
    IDSs(-1),
    IDSids(-1),

    IGSg(-1),
    IGSs(-1),
    IGSigs(-1),

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    f_DdPtr  (0),	  q_DdPtr  (0),
    f_DdpPtr (0),	  q_DdpPtr (0),
    f_DspPtr (0),	  q_DspPtr (0),
    f_DgpPtr (0),	  q_DgpPtr (0),
    f_DbpPtr (0),	  q_DbpPtr (0),
    f_DidsPtr(0),

    f_GEgePtr (0),	  q_GEgePtr (0),
    f_GEdpPtr (0),	  q_GEdpPtr (0),
    f_GEspPtr (0),	  q_GEspPtr (0),
    f_GEgpPtr (0),	  q_GEgpPtr (0),
    f_GEgmPtr (0),	  q_GEgmPtr (0),
    f_GEbpPtr (0),	  q_GEbpPtr (0),
    f_GEigsPtr (0),

    f_SsPtr  (0),	  q_SsPtr  (0),
    f_SdpPtr (0),	  q_SdpPtr (0),
    f_SspPtr (0),	  q_SspPtr (0),
    f_SgpPtr (0),	  q_SgpPtr (0),
    f_SbpPtr (0),	  q_SbpPtr (0),
    f_SibsPtr(0),
    f_SidsPtr(0),
    f_SigsPtr(0),

    f_BbPtr  (0),	  q_BbPtr  (0),
    f_BbpPtr (0),	  q_BbpPtr (0),
    f_BsbPtr (0),	  q_BsbPtr (0),
    f_BdbPtr (0),	  q_BdbPtr (0),
    f_BibsPtr(0),

    f_DPdPtr  (0),	  q_DPdPtr  (0),
    f_DPdpPtr (0),	  q_DPdpPtr (0),
    f_DPspPtr (0),	  q_DPspPtr (0),
    f_DPgpPtr (0),	  q_DPgpPtr (0),
    f_DPgmPtr (0),	  q_DPgmPtr (0),
    f_DPbpPtr (0),	  q_DPbpPtr (0),
    f_DPdbPtr (0),	  q_DPdbPtr (0),


    f_DPqPtr (0),	    q_DPqPtr (0),


    f_SPsPtr  (0),	  q_SPsPtr  (0),
    f_SPdpPtr (0),	  q_SPdpPtr (0),
    f_SPspPtr (0),	  q_SPspPtr (0),
    f_SPgpPtr (0),	  q_SPgpPtr (0),
    f_SPgmPtr (0),	  q_SPgmPtr (0),
    f_SPbpPtr (0),	  q_SPbpPtr (0),
    f_SPsbPtr (0),	  q_SPsbPtr (0),

    f_SPqPtr (0),	    q_SPqPtr (0),

    f_GPgePtr(0),	  q_GPgePtr(0),
    f_GPdpPtr(0),	  q_GPdpPtr(0),
    f_GPspPtr(0),	  q_GPspPtr(0),
    f_GPgpPtr(0),	  q_GPgpPtr(0),
    f_GPgmPtr(0),	  q_GPgmPtr(0),
    f_GPbpPtr(0),	  q_GPbpPtr(0),

    f_GPqPtr(0),	    q_GPqPtr(0),

    f_GMgePtr(0),	  q_GMgePtr(0),
    f_GMdpPtr(0),	  q_GMdpPtr(0),
    f_GMspPtr(0),	  q_GMspPtr(0),
    f_GMgpPtr(0),	  q_GMgpPtr(0),
    f_GMgmPtr(0),	  q_GMgmPtr(0),
    f_GMbpPtr(0),	  q_GMbpPtr(0),

    f_BPbPtr (0),	  q_BPbPtr (0),
    f_BPdpPtr(0),	  q_BPdpPtr(0),
    f_BPspPtr(0),	  q_BPspPtr(0),
    f_BPgpPtr(0),	  q_BPgpPtr(0),
    f_BPgmPtr(0),	  q_BPgmPtr(0),
    f_BPbpPtr(0),	  q_BPbpPtr(0),
    f_BPsbPtr(0),	  q_BPsbPtr(0),
    f_BPdbPtr(0),	  q_BPdbPtr(0),

    f_SBbPtr (0),	  q_SBbPtr (0),
    f_SBspPtr(0),	  q_SBspPtr(0),
    f_SBbpPtr(0),	  q_SBbpPtr(0),
    f_SBsbPtr(0),	  q_SBsbPtr(0),

    f_DBbPtr (0),	  q_DBbPtr (0),
    f_DBdpPtr(0),	  q_DBdpPtr(0),
    f_DBbpPtr(0),	  q_DBbpPtr(0),
    f_DBdbPtr(0),	  q_DBdbPtr(0),

    f_QdpPtr(0),	    q_QdpPtr(0),
    f_QspPtr(0),	    q_QspPtr(0),
    f_QgpPtr(0),	    q_QgpPtr(0),
    f_QbpPtr(0),	    q_QbpPtr(0),
    f_QqPtr (0),	    q_QqPtr (0),
    f_IBSbPtr(0),
    f_IBSsPtr(0),
    f_IBSibsPtr(0),

    f_IDSdPtr(0),
    f_IDSsPtr(0),
    f_IDSidsPtr(0),
    f_IGSgPtr(0),
    f_IGSsPtr(0),
    f_IGSigsPtr(0),
#endif

    updateTemperatureCalled_ (false),

    jacStamp(0, std::vector<int>(0)),
    jacMap(0),
    jacMap2(0, std::vector<int>(0)),

    blockHomotopyID                      (0),
    randomPerturb                        (0.0),

    ceqdrn_Jdxp(0.0),
    ceqbd_Jdxp(0.0),
    ceqbs_Jdxp(0.0),
    Istoteq_Jdxp(0.0),
    Idtoteq_Jdxp(0.0),
    Ibtoteq_Jdxp(0.0),
    Igtoteq_Jdxp(0.0) ,
    ceqgcrg_Jdxp(0.0),
    ceqgstot_Jdxp(0.0),
    ceqgdtot_Jdxp(0.0),
    ceqjs_Jdxp(0.0),
    ceqjd_Jdxp(0.0),
    T0(0.0)
{
  numIntVars   = 3;
  numExtVars   = 4;
  numStateVars = 17;
  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 4;    // this is the space to allocate if lead current or power is needed.

  devConMap.resize(4);
  devConMap[0] = 1;
  devConMap[1] = 2;
  devConMap[2] = 1;
  devConMap[3] = 3;

  blockHomotopyID =
    devSupport.getGainScaleBlockID(getDeviceOptions().numGainScaleBlocks);
  randomPerturb =
    devSupport.getRandomPerturbation();


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

  if (given("TRNQSMOD"))
  {
    UserWarning(*this) << "nsqMod = 1.  Not allowed yet.  Setting to 0";
  }

  if (getDeviceOptions().verboseLevel > 0 && (l > model_.Lmax || l < model_.Lmin))
  {
    UserWarning(*this) << "Channel length out of range";
  }

  if (getDeviceOptions().verboseLevel > 0 && (w > model_.Wmax || w < model_.Wmin))
  {
    UserWarning(*this) << "Channel width out of range";
  }

  // Note important difference from the BSIM3 --- we do NOT keep any
  // static jacstamps for special cases, we build the entire real jacstamp
  // for this particular device, starting from the full general case, then
  // mapping away what's unneeded.  There are simply too many special cases
  // to do it the old way.

  if( jacStamp.empty() )
  {
    jacStamp.resize(11);
    jacStamp[0].resize(5);   // Drain Row
    jacStamp[0][0]=0;        // Dd    (drain-drain)
    jacStamp[0][1]=4;        // Ddp   (drain-drain')
    jacStamp[0][2]=5;        // Dsp   (drain-source')
    jacStamp[0][3]=6;        // Dgp   (drain-gate')
    jacStamp[0][4]=8;        // Dbp   (drain-body')

    jacStamp[1].resize(6);   // GateExt Row
    jacStamp[1][0]= 1;       // GEge
    jacStamp[1][1]= 4;       // GEdp
    jacStamp[1][2]= 5;       // GEsp
    jacStamp[1][3]= 6;       // GEgp
    jacStamp[1][4]= 7;       // GEgm
    jacStamp[1][5]= 8;       // GEbp

    jacStamp[2].resize(5);   // Source Row
    jacStamp[2][0]=2;        // Ss    (source-source)
    jacStamp[2][1]=4;        // Sdp   (source-drain')
    jacStamp[2][2]=5;        // Ssp   (source-source')
    jacStamp[2][3]=6;        // Sgp   (source-gate')
    jacStamp[2][4]=8;        // Sbp   (source-body')

    jacStamp[3].resize(4);   // Body Row
    jacStamp[3][0]= 3;       // Bb
    jacStamp[3][1]= 8;       // Bbp
    jacStamp[3][2]= 9;       // Bsb
    jacStamp[3][3]= 10;      // Bdb

    // Optional terminal resistors (rdsMod):
    jacStamp[4].resize(7);   // Drain' Row
    jacStamp[4][0]=0;        // DPd   (drain'-drain)
    jacStamp[4][1]=4;        // DPdp  (drain'-drain')
    jacStamp[4][2]=5;        // DPsp  (drain'-source')
    jacStamp[4][3]=6;        // DPgp  (drain'-gate')
    jacStamp[4][4]=7;        // DPgm  (drain'-gateMid)
    jacStamp[4][5]=8;        // DPbp  (drain'-body')
    jacStamp[4][6]=10;       // DPdb  (drain'-drainBody)

    jacStamp[5].resize(7);   // Source' Row
    jacStamp[5][0]=2;        // SPs   (source'-source)
    jacStamp[5][1]=4;        // SPdp  (source'-drain')
    jacStamp[5][2]=5;        // SPsp  (source'-source')
    jacStamp[5][3]=6;        // SPgp  (source'-gate')
    jacStamp[5][4]=7;        // SPgm  (source'-gateMid)
    jacStamp[5][5]=8;        // SPbp  (source'-body')
    jacStamp[5][6]=9;        // SPsb  (source'-sourceBody)

    // Optional gate resistors (rgateMod)
    jacStamp[6].resize(6);   // Gate' Row (needed of rgateMod > 0)
    jacStamp[6][0]=1;        // GPge  (gate'-gateExt)
    jacStamp[6][1]=4;        // GPdp  (gate'-drain')
    jacStamp[6][2]=5;        // GPsp  (gate'-source')
    jacStamp[6][3]=6;        // GPgp  (gate'-gate')
    jacStamp[6][4]=7;        // GPgm  (gate'-gateMid)
    jacStamp[6][5 ]=8;        // GPbp  (gate'-body')

    jacStamp[7].resize(6);   // GateMid Row (needed of rgateMod == 3)
    jacStamp[7][0]= 1;        // GMge
    jacStamp[7][1]= 4;        // GMdp
    jacStamp[7][2]= 5;        // GMsp
    jacStamp[7][3]= 6;        // GMgp
    jacStamp[7][4]= 7;        // GMgm
    jacStamp[7][5]= 8;        // GMbp

    // Optional rbodyMod variables:
    jacStamp[8].resize(8);   // Body' Row
    jacStamp[8][0]=3;        // BPb   (body'-body)
    jacStamp[8][1]=4;        // BPdp  (body'-drain')
    jacStamp[8][2]=5;        // BPsp  (body'-source')
    jacStamp[8][3]=6;        // BPgp  (body'-gate')
    jacStamp[8][4]=7;        // BPgm  (body'-gateMid)
    jacStamp[8][5]=8;        // BPbp  (body'-body')
    jacStamp[8][6]=9;        // BPsb  (body'-sourceBody)
    jacStamp[8][7]=10;       // BPdb  (body'-drainBody)

    jacStamp[9].resize(4);   // SourceBody Row
    jacStamp[9][0] = 3;      // SBb
    jacStamp[9][1] = 5;      // SBsp
    jacStamp[9][2] = 8;      // SBbp
    jacStamp[9][3] = 9;      // SBsb

    jacStamp[10].resize(4);   // DrainBody Row
    jacStamp[10][0] = 3;      // DBb
    jacStamp[10][1] = 4;      // DBdp
    jacStamp[10][2] = 8;      // DBbp
    jacStamp[10][3] = 10;     // DBdb

    // Optional nqs row:
    if( trnqsMod )
    {
      int oldSize = jacStamp.size();
      int newMaxCol=oldSize;
      int newSize = oldSize+1;
      jacStamp.resize(newSize);
      jacStamp[newMaxCol].resize(5);   // Charge Row
      jacStamp[newMaxCol][0] = 4;      // Qdp
      jacStamp[newMaxCol][1] = 5;      // Qsp
      jacStamp[newMaxCol][2] = 6;      // Qgp
      jacStamp[newMaxCol][3] = 8;      // Qbp
      jacStamp[newMaxCol][4] = newMaxCol;     // Qq

      // extra columns in other rows:
      jacStamp[6].resize(7);           // Gate' Row
      jacStamp[6][6]=newMaxCol;        // GPq   (gate'-charge)

      jacStamp[4].resize(8);           // Drain' Row
      jacStamp[4][7]=newMaxCol;        // DPq   (drain'-charge)

      jacStamp[5].resize(8);           // Source' Row
      jacStamp[5][7]=newMaxCol;        // SPq   (source'-charge)
    }

    if (icVBSGiven)
    {
      int oldSize=jacStamp.size();
      int icVBSCol=oldSize;
      int newSize=oldSize+1;

      jacStamp.resize(newSize);
      // New Ibs row:
      jacStamp[icVBSCol].resize(3);
      jacStamp[icVBSCol][0]=2;      // Ibs-Source
      jacStamp[icVBSCol][1]=3;      // Ibs-Body
      jacStamp[icVBSCol][2]=icVBSCol; // Ibs-Ibs


      // Fix up Source row:
      int sourceSize=jacStamp[2].size();
      int newSourceCol=sourceSize;
      sourceSize++;
      jacStamp[2].resize(sourceSize);
      jacStamp[2][newSourceCol]=icVBSCol;  // Source-IBS

      // Fix up Body row:
      int bodySize=jacStamp[3].size();
      int newBodyCol=bodySize;
      bodySize++;
      jacStamp[3].resize(bodySize);
      jacStamp[3][newBodyCol]=icVBSCol;    // Body-IBS
    }

    if (icVDSGiven)
    {
      int oldSize=jacStamp.size();
      int icVDSCol=oldSize;
      int newSize=oldSize+1;

      jacStamp.resize(newSize);
      // New Ids row:
      jacStamp[icVDSCol].resize(3);
      jacStamp[icVDSCol][0]=0;      // Ids-Drain
      jacStamp[icVDSCol][1]=2;      // Ids-Source
      jacStamp[icVDSCol][2]=icVDSCol; // Ids-Ids


      // Fix up Source row:
      int sourceSize=jacStamp[2].size();
      int newSourceCol=sourceSize;
      sourceSize++;
      jacStamp[2].resize(sourceSize);
      jacStamp[2][newSourceCol]=icVDSCol; // Source-IDS

      // Fix up Drain row:
      int drainSize=jacStamp[0].size();
      int newDrainCol=drainSize;
      drainSize++;
      jacStamp[0].resize(drainSize);
      jacStamp[0][newDrainCol]=icVDSCol;  // Drain-IDS
    }

    if (icVGSGiven)
    {
      int oldSize=jacStamp.size();
      int icVGSCol=oldSize;
      int newSize=oldSize+1;

      jacStamp.resize(newSize);
      // New Igs row:
      jacStamp[icVGSCol].resize(3);
      jacStamp[icVGSCol][0]=1;      // Ids-GateExt
      jacStamp[icVGSCol][1]=2;      // Ids-Source
      jacStamp[icVGSCol][2]=icVGSCol; // Igs-Igs


      // Fix up Source row:
      int sourceSize=jacStamp[2].size();
      int newSourceCol=sourceSize;
      sourceSize++;
      jacStamp[2].resize(sourceSize);
      jacStamp[2][newSourceCol]=icVGSCol;  // Source-IGS

      // Fix up gate row:
      int gateSize=jacStamp[1].size();
      int newGateCol=gateSize;
      gateSize++;
      jacStamp[1].resize(gateSize);
      jacStamp[1][newGateCol]=icVGSCol;    // Gate-IGS
    }

    jacMap.clear();
    jacMap2.clear();
    jacMap.resize(jacStamp.size());
    jacMap2.resize(jacStamp.size());

    int mapSize = jacMap.size();
    for (int i=0;i<mapSize;++i)
    {
      jacMap[i]=i;
      jacMap2[i].resize(jacStamp[i].size());
      for (int j=0;j<jacStamp[i].size();++j)
      {
        jacMap2[i][j] = j;
      }
    }

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "About to remap away optional nodes from the jacStamp!" << std::endl;
    debugJacStampOutput ();
    }

    // Selectively map away bits and pieces of the
    // jacobian stamp based on absent resistances

    // temporary stamps and maps
    std::vector< std::vector<int> > tempStamp;
    std::vector<int> tempMap;
    std::vector< std::vector<int> > tempMap2;

    int OriginalSize = jacMap.size();

    // Body mod:
    if (!rbodyMod) // remove body, sourcebody, drainbody
    {
      // Map away the drainBody (originally 10)
      tempStamp.clear(); tempMap.clear(); tempMap2.clear();
      jacStampMap(jacStamp, jacMap, jacMap2,
                  tempStamp, tempMap, tempMap2,
                  jacMap[10] , jacMap[3],
                  OriginalSize);

      // now move the new stuff into the old place
      jacStamp = tempStamp; jacMap = tempMap; jacMap2 = tempMap2;

      // Map away the sourceBody (originally 9)
      tempStamp.clear(); tempMap.clear(); tempMap2.clear();
      jacStampMap(jacStamp, jacMap, jacMap2,
                  tempStamp, tempMap, tempMap2,
                  jacMap[9] , jacMap[3],
                  OriginalSize);

      // now move the new stuff into the old place
      jacStamp = tempStamp; jacMap = tempMap; jacMap2 = tempMap2;

      // Map  Body' (originally 8) onto Body
      tempStamp.clear(); tempMap.clear(); tempMap2.clear();
      jacStampMap(jacStamp, jacMap, jacMap2,
                  tempStamp, tempMap, tempMap2,
                  jacMap[8] , jacMap[3],
                  OriginalSize);

      // now move the new stuff into the old place
      jacStamp = tempStamp; jacMap = tempMap; jacMap2 = tempMap2;
    }

    // rgateMod.  The gate resistor has 4 options
    // rgateMod==0   no gate resistor.
    // rgateMod==1   linear gate resistor
    // rgateMod==2   nonlinear gate resistor
    // rgateMod==3   2 gate resistors, in series.:

    if (rgateMod != 3) // remove gateMid.
    {
      // Map away the gateMid (originally 7)
      tempStamp.clear(); tempMap.clear(); tempMap2.clear();
      jacStampMap(jacStamp, jacMap, jacMap2,
                  tempStamp, tempMap, tempMap2,
                  jacMap[7] , jacMap[1],
                  OriginalSize);

      // now move the new stuff into the old place
      jacStamp = tempStamp; jacMap = tempMap; jacMap2 = tempMap2;
    }

    if (rgateMod < 1) // remove gatePrime
    {
      // Map away the gatePrime (originally 6)
      tempStamp.clear(); tempMap.clear(); tempMap2.clear();
      jacStampMap(jacStamp, jacMap, jacMap2,
                  tempStamp, tempMap, tempMap2,
                  jacMap[6] , jacMap[1],
                  OriginalSize);

      // now move the new stuff into the old place
      jacStamp = tempStamp; jacMap = tempMap; jacMap2 = tempMap2;
    }

    // drain and source terminal resistors.
    if(!sourceMOSFET_B4Exists)
    {
      // Map away the source' (originally 5), into the source (2):
      tempStamp.clear(); tempMap.clear(); tempMap2.clear();
      jacStampMap(jacStamp, jacMap, jacMap2,
                  tempStamp, tempMap, tempMap2,
                  jacMap[5], jacMap[2],
                  OriginalSize);
      // now move the new stuff into the old place
      jacStamp = tempStamp; jacMap = tempMap; jacMap2 = tempMap2;
    }


    if (!drainMOSFET_B4Exists)
    {
      // Map away the drain' (always 4), into the drain (0):
      tempStamp.clear(); tempMap.clear(); tempMap2.clear();
      jacStampMap(jacStamp, jacMap, jacMap2,
                  tempStamp, tempMap, tempMap2,
                  jacMap[4] , jacMap[0],
                  OriginalSize);

      // now move the new stuff into the old place
      jacStamp = tempStamp; jacMap = tempMap; jacMap2 = tempMap2;
    }


    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "Done remap away optional nodes from the jacStamp!" << std::endl;
    debugJacStampOutput ();
    }
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                                            const std::vector<int> & extLIDVecRef )
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << section_divider << std::endl
                 << "  In Instance::register LIDs" << std::endl
                 << "  name             = " << getName() << std::endl
                 << "  number of internal variables: " << intLIDVecRef.size() << std::endl
                 << "  number of external variables: " << extLIDVecRef.size() << std::endl;
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Drain    = extLIDVec[0];
  li_GateExt  = extLIDVec[1];
  li_Source   = extLIDVec[2];
  li_Body     = extLIDVec[3];

  int intLoc = 0;

  if (drainMOSFET_B4Exists)
  {
    li_DrainPrime = intLIDVec[intLoc++];
  }
  else
  {
    li_DrainPrime = li_Drain;
  }

  if (sourceMOSFET_B4Exists)
  {
    li_SourcePrime = intLIDVec[intLoc++];
  }
  else
  {
    li_SourcePrime = li_Source;
  }

  if (rgateMod>0)
  {
    li_GatePrime = intLIDVec[intLoc++];
  }
  else
  {
    li_GatePrime = li_GateExt;
  }

  if (rgateMod==3)
  {
    li_GateMid   = intLIDVec[intLoc++];
  }
  else
  {
    li_GateMid = li_GateExt;
  }

  if ( rbodyMod )
  {
    li_BodyPrime  = intLIDVec[intLoc++];
    li_SourceBody = intLIDVec[intLoc++];
    li_DrainBody  = intLIDVec[intLoc++];
  }
  else
  {
    li_BodyPrime  = li_Body;
    li_SourceBody = li_Body;
    li_DrainBody  = li_Body;
  }

  if( trnqsMod )
  {
    li_Charge = intLIDVec[intLoc++];
  }

  if( icVBSGiven )
  {
    if( li_Body == li_Source )
    {
      UserError(*this) << "Tried to specify an initial condition on V_Bulk_Source when Bulk and Source nodes are the same node";
    }
    li_Ibs = intLIDVec[intLoc++];
  }

  if( icVDSGiven )
  {
    if( li_Drain == li_Source )
    {
      UserError(*this) << "Tried to specify an initial condition on V_Drain_Source when Drain and Source nodes are the same node";
    }
    li_Ids = intLIDVec[intLoc++];
  }

  if( icVGSGiven )
  {
    if( li_GateExt == li_Source )
    {
      UserError(*this) << "Tried to specify an initial condition on V_Gate_Source when Gate and Source nodes are the same node";
    }
   li_Igs = intLIDVec[intLoc++];
  }


  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "\n  local variable indices:\n";
    Xyce::dout() << " li_Drain       = " << li_Drain << std::endl;
    Xyce::dout() << " li_GateExt     = " << li_GateExt << std::endl;
    Xyce::dout() << " li_Source      = " << li_Source << std::endl;
    Xyce::dout() << " li_Body        = " << li_Body << std::endl;
    Xyce::dout() << " li_DrainPrime  = " << li_DrainPrime << std::endl;
    Xyce::dout() << " li_GatePrime   = " << li_GatePrime << std::endl;
    Xyce::dout() << " li_GateMid     = " << li_GateMid << std::endl;
    Xyce::dout() << " li_SourcePrime = " << li_SourcePrime << std::endl;
    Xyce::dout() << " li_BodyPrime   = " << li_BodyPrime << std::endl;
    Xyce::dout() << " li_DrainBody   = " << li_DrainBody << std::endl;
    Xyce::dout() << " li_SourceBody  = " << li_SourceBody << std::endl;
    Xyce::dout() << " li_Charge      = " << li_Charge << std::endl;
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
  if (drainMOSFET_B4Exists)
    addInternalNode(symbol_table, li_DrainPrime, getName(), "drainprime");

  if (sourceMOSFET_B4Exists)
    addInternalNode(symbol_table, li_SourcePrime, getName(), "sourceprime");

  if (rgateMod>0)
    addInternalNode(symbol_table, li_GatePrime, getName(), "GatePrime");

  if (rgateMod==3)
    addInternalNode(symbol_table, li_GateMid, getName(), "MidGate");

  if (rbodyMod)
  {
    addInternalNode(symbol_table, li_BodyPrime, getName(), "BodyPrime");
    addInternalNode(symbol_table, li_SourceBody, getName(), "SourceBody");
    addInternalNode(symbol_table, li_DrainBody, getName(), "DrainBody");
  }

  if (trnqsMod)
    addInternalNode(symbol_table, li_Charge, getName(), "charge");

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
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl
                 << section_divider << std::endl
                 << "  In Instance::registerStateLIDs" << std::endl
                 << "  name             = " << getName() << std::endl
                 << "  Number of State LIDs: " << staLIDVecRef.size() << std::endl;
  }

  // Copy over the global ID lists:
  staLIDVec = staLIDVecRef;

  int lid=0;
  // Intrinsic capacitors:
  li_state_qb        = staLIDVec[lid++];
  li_state_qg        = staLIDVec[lid++];
  li_state_qd        = staLIDVec[lid++];

  if (rgateMod == 3)
  {
    li_state_qgmid = staLIDVec[lid++];
  }

  // Parasitic capacitors:
  if (rbodyMod)
  {
    li_state_qbs        = staLIDVec[lid++];
    li_state_qbd        = staLIDVec[lid++];
  }

  if( trnqsMod )
  {
    li_state_qcheq     = staLIDVec[lid++];
    li_state_qcdump       = staLIDVec[lid++];
  }

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
// Creation Date : 12/9/11
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs( const std::vector<int> & stoLIDVecRef )
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());

  // Copy over the global ID lists:
  stoLIDVec = stoLIDVecRef;

  int lid=0;
  // Voltage drops:
  li_store_vbd  = stoLIDVec[lid++];
  li_store_vbs  = stoLIDVec[lid++];
  li_store_vgs  = stoLIDVec[lid++];
  li_store_vds  = stoLIDVec[lid++];
  li_store_vges = stoLIDVec[lid++];
  li_store_vgms = stoLIDVec[lid++];
  li_store_vdes = stoLIDVec[lid++];
  li_store_vses = stoLIDVec[lid++];
  li_store_vdbs = stoLIDVec[lid++];
  li_store_vsbs = stoLIDVec[lid++];
  li_store_vdbd = stoLIDVec[lid++];
  li_store_vged = stoLIDVec[lid++];
  li_store_vgmd = stoLIDVec[lid++];

  // transconductance, and other outputs:
  li_store_gm    = stoLIDVec[lid++];
  li_store_Vds   = stoLIDVec[lid++];
  li_store_Vgs   = stoLIDVec[lid++];
  li_store_Vbs   = stoLIDVec[lid++];
  li_store_Vdsat = stoLIDVec[lid++];
  li_store_Vth   = stoLIDVec[lid++];

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
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
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
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  std::vector<int> & map = jacMap;
  std::vector< std::vector<int> > & map2 = jacMap2;;

  Dd   =  jacLIDVec[map[0]][map2[0][0]];
  Ddp  =  jacLIDVec[map[0]][map2[0][1]];
  Dsp  =  jacLIDVec[map[0]][map2[0][2]];
  Dgp  =  jacLIDVec[map[0]][map2[0][3]];
  Dbp  =  jacLIDVec[map[0]][map2[0][4]];
  if (icVDSGiven)
  {
    Dids = jacLIDVec[map[0]][map2[0][5]];
  }

  GEge =  jacLIDVec[map[1]][map2[1][0]];
  GEdp =  jacLIDVec[map[1]][map2[1][1]];
  GEsp =  jacLIDVec[map[1]][map2[1][2]];
  GEgp =  jacLIDVec[map[1]][map2[1][3]];
  GEgm =  jacLIDVec[map[1]][map2[1][4]];
  GEbp =  jacLIDVec[map[1]][map2[1][5]];
  if (icVGSGiven)
  {
    GEigs = jacLIDVec[map[1]][map2[1][6]];
  }

  Ss   =  jacLIDVec[map[2]][map2[2][0]];
  Sdp  =  jacLIDVec[map[2]][map2[2][1]];
  Ssp  =  jacLIDVec[map[2]][map2[2][2]];
  Sgp  =  jacLIDVec[map[2]][map2[2][3]];
  Sbp  =  jacLIDVec[map[2]][map2[2][4]];
  int currentCol=5;
  // these must be assigned in the same order that they were set up
  // in the jacstamp set-up in the constructor
  if (icVBSGiven)
  {
    Sibs = jacLIDVec[map[2]][map2[2][currentCol++]];
  }
  if (icVDSGiven)
  {
    Sids = jacLIDVec[map[2]][map2[2][currentCol++]];
  }
  if (icVGSGiven)
  {
    Sigs = jacLIDVec[map[2]][map2[2][currentCol++]];
  }

  Bb   =  jacLIDVec[map[3]][map2[3][0]];
  Bbp  =  jacLIDVec[map[3]][map2[3][1]];
  Bsb  =  jacLIDVec[map[3]][map2[3][2]];
  Bdb  =  jacLIDVec[map[3]][map2[3][3]];
  if (icVBSGiven)
  {
    Bibs = jacLIDVec[map[3]][map2[3][4]];
  }

  DPd  =  jacLIDVec[map[4]][map2[4][0]];
  DPdp =  jacLIDVec[map[4]][map2[4][1]];
  DPsp =  jacLIDVec[map[4]][map2[4][2]];
  DPgp =  jacLIDVec[map[4]][map2[4][3]];
  DPgm =  jacLIDVec[map[4]][map2[4][4]];
  DPbp =  jacLIDVec[map[4]][map2[4][5]];
  DPdb =  jacLIDVec[map[4]][map2[4][6]];
  if( trnqsMod )
  {
    DPq  =  jacLIDVec[map[4]][map2[4][7]];
  }

  SPs  =  jacLIDVec[map[5]][map2[5][0]];
  SPdp =  jacLIDVec[map[5]][map2[5][1]];
  SPsp =  jacLIDVec[map[5]][map2[5][2]];
  SPgp =  jacLIDVec[map[5]][map2[5][3]];
  SPgm =  jacLIDVec[map[5]][map2[5][4]];
  SPbp =  jacLIDVec[map[5]][map2[5][5]];
  SPsb =  jacLIDVec[map[5]][map2[5][6]];

  if( trnqsMod )
  {
    SPq  =  jacLIDVec[map[5]][map2[5][7]];
  }

  GPge =  jacLIDVec[map[6]][map2[6][0]];
  GPdp =  jacLIDVec[map[6]][map2[6][1]];
  GPsp =  jacLIDVec[map[6]][map2[6][2]];
  GPgp =  jacLIDVec[map[6]][map2[6][3]];
  GPgm =  jacLIDVec[map[6]][map2[6][4]];
  GPbp =  jacLIDVec[map[6]][map2[6][5]];

  if( trnqsMod )
  {
    GPq  =  jacLIDVec[map[6]][map2[6][6]];
  }

  GMge =  jacLIDVec[map[7]][map2[7][0]];
  GMdp =  jacLIDVec[map[7]][map2[7][1]];
  GMsp =  jacLIDVec[map[7]][map2[7][2]];
  GMgp =  jacLIDVec[map[7]][map2[7][3]];
  GMgm =  jacLIDVec[map[7]][map2[7][4]];
  GMbp =  jacLIDVec[map[7]][map2[7][5]];

  BPb  =  jacLIDVec[map[8]][map2[8][0]];
  BPdp =  jacLIDVec[map[8]][map2[8][1]];
  BPsp =  jacLIDVec[map[8]][map2[8][2]];
  BPgp =  jacLIDVec[map[8]][map2[8][3]];
  BPgm =  jacLIDVec[map[8]][map2[8][4]];
  BPbp =  jacLIDVec[map[8]][map2[8][5]];
  BPsb =  jacLIDVec[map[8]][map2[8][6]];
  BPdb =  jacLIDVec[map[8]][map2[8][7]];

  SBb  =  jacLIDVec[map[9]][map2[9][0]];
  SBsp =  jacLIDVec[map[9]][map2[9][1]];
  SBbp =  jacLIDVec[map[9]][map2[9][2]];
  SBsb =  jacLIDVec[map[9]][map2[9][3]];

  DBb  =  jacLIDVec[map[10]][map2[10][0]];
  DBdp =  jacLIDVec[map[10]][map2[10][1]];
  DBbp =  jacLIDVec[map[10]][map2[10][2]];
  DBdb =  jacLIDVec[map[10]][map2[10][3]];

  int currentRow=11;
  // Optional nqs row:
  if( trnqsMod )
  {
    Qdp  =  jacLIDVec[map[currentRow]][map2[currentRow][0]];
    Qsp  =  jacLIDVec[map[currentRow]][map2[currentRow][1]];
    Qgp  =  jacLIDVec[map[currentRow]][map2[currentRow][2]];
    Qbp  =  jacLIDVec[map[currentRow]][map2[currentRow][3]];
    Qq   =  jacLIDVec[map[currentRow]][map2[currentRow][4]];
    currentRow++;
  }

  if (icVBSGiven)
  {
    IBSs = jacLIDVec[map[currentRow]][map2[currentRow][0]];
    IBSb = jacLIDVec[map[currentRow]][map2[currentRow][1]];
    IBSibs = jacLIDVec[map[currentRow]][map2[currentRow][2]];
    currentRow++;
  }

  if (icVDSGiven)
  {
    IDSd = jacLIDVec[map[currentRow]][map2[currentRow][0]];
    IDSs = jacLIDVec[map[currentRow]][map2[currentRow][1]];
    IDSids = jacLIDVec[map[currentRow]][map2[currentRow][2]];
    currentRow++;
  }

  if (icVGSGiven)
  {
    IGSg = jacLIDVec[map[currentRow]][map2[currentRow][0]];
    IGSs = jacLIDVec[map[currentRow]][map2[currentRow][1]];
    IGSigs = jacLIDVec[map[currentRow]][map2[currentRow][2]];
    currentRow++;
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/01/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);

  f_DdPtr  =  &(dFdx[li_Drain][Dd  ]);	  q_DdPtr  =  &(dQdx[li_Drain][Dd  ]);
  f_DdpPtr =  &(dFdx[li_Drain][Ddp ]);	  q_DdpPtr =  &(dQdx[li_Drain][Ddp ]);
  f_DspPtr =  &(dFdx[li_Drain][Dsp ]);	  q_DspPtr =  &(dQdx[li_Drain][Dsp ]);
  f_DgpPtr =  &(dFdx[li_Drain][Dgp ]);	  q_DgpPtr =  &(dQdx[li_Drain][Dgp ]);
  f_DbpPtr =  &(dFdx[li_Drain][Dbp ]);	  q_DbpPtr =  &(dQdx[li_Drain][Dbp ]);
  if (icVDSGiven)
  {
    f_DidsPtr = &(dFdx[li_Drain][Dids]);
  }

  f_GEgePtr =  &(dFdx[li_GateExt][GEge ]);	  q_GEgePtr =  &(dQdx[li_GateExt][GEge ]);
  f_GEdpPtr =  &(dFdx[li_GateExt][GEdp ]);	  q_GEdpPtr =  &(dQdx[li_GateExt][GEdp ]);
  f_GEspPtr =  &(dFdx[li_GateExt][GEsp ]);	  q_GEspPtr =  &(dQdx[li_GateExt][GEsp ]);
  f_GEgpPtr =  &(dFdx[li_GateExt][GEgp ]);	  q_GEgpPtr =  &(dQdx[li_GateExt][GEgp ]);
  f_GEgmPtr =  &(dFdx[li_GateExt][GEgm ]);	  q_GEgmPtr =  &(dQdx[li_GateExt][GEgm ]);
  f_GEbpPtr =  &(dFdx[li_GateExt][GEbp ]);	  q_GEbpPtr =  &(dQdx[li_GateExt][GEbp ]);
  if (icVGSGiven)
  {
    f_GEigsPtr = &(dFdx[li_GateExt][GEigs]);
  }

  f_SsPtr  =  &(dFdx[li_Source][Ss  ]);	  q_SsPtr  =  &(dQdx[li_Source][Ss  ]);
  f_SdpPtr =  &(dFdx[li_Source][Sdp ]);	  q_SdpPtr =  &(dQdx[li_Source][Sdp ]);
  f_SspPtr =  &(dFdx[li_Source][Ssp ]);	  q_SspPtr =  &(dQdx[li_Source][Ssp ]);
  f_SgpPtr =  &(dFdx[li_Source][Sgp ]);	  q_SgpPtr =  &(dQdx[li_Source][Sgp ]);
  f_SbpPtr =  &(dFdx[li_Source][Sbp ]);	  q_SbpPtr =  &(dQdx[li_Source][Sbp ]);
  if (icVBSGiven)
  {
    f_SibsPtr = &(dFdx[li_Source][Sibs]);
  }
  if (icVDSGiven)
  {
    f_SidsPtr = &(dFdx[li_Source][Sids]);
  }
  if (icVGSGiven)
  {
    f_SigsPtr = &(dFdx[li_Source][Sigs]);
  }

  f_BbPtr  =  &(dFdx[li_Body][Bb  ]);	  q_BbPtr  =  &(dQdx[li_Body][Bb  ]);
  f_BbpPtr =  &(dFdx[li_Body][Bbp ]);	  q_BbpPtr =  &(dQdx[li_Body][Bbp ]);
  f_BsbPtr =  &(dFdx[li_Body][Bsb ]);	  q_BsbPtr =  &(dQdx[li_Body][Bsb ]);
  f_BdbPtr =  &(dFdx[li_Body][Bdb ]);	  q_BdbPtr =  &(dQdx[li_Body][Bdb ]);
  if (icVBSGiven)
  {
    f_BibsPtr = &(dFdx[li_Body][Bibs]);
  }

  f_DPdPtr  =  &(dFdx[li_DrainPrime][DPd  ]);	  q_DPdPtr  =  &(dQdx[li_DrainPrime][DPd  ]);
  f_DPdpPtr =  &(dFdx[li_DrainPrime][DPdp ]);	  q_DPdpPtr =  &(dQdx[li_DrainPrime][DPdp ]);
  f_DPspPtr =  &(dFdx[li_DrainPrime][DPsp ]);	  q_DPspPtr =  &(dQdx[li_DrainPrime][DPsp ]);
  f_DPgpPtr =  &(dFdx[li_DrainPrime][DPgp ]);	  q_DPgpPtr =  &(dQdx[li_DrainPrime][DPgp ]);
  f_DPgmPtr =  &(dFdx[li_DrainPrime][DPgm ]);	  q_DPgmPtr =  &(dQdx[li_DrainPrime][DPgm ]);
  f_DPbpPtr =  &(dFdx[li_DrainPrime][DPbp ]);	  q_DPbpPtr =  &(dQdx[li_DrainPrime][DPbp ]);
  f_DPdbPtr =  &(dFdx[li_DrainPrime][DPdb ]);	  q_DPdbPtr =  &(dQdx[li_DrainPrime][DPdb ]);
  if( trnqsMod )
  {
    f_DPqPtr =  &(dFdx[li_DrainPrime][DPq ]);	    q_DPqPtr =  &(dQdx[li_DrainPrime][DPq ]);
  }

  f_SPsPtr  =  &(dFdx[li_SourcePrime][SPs  ]);	  q_SPsPtr  =  &(dQdx[li_SourcePrime][SPs  ]);
  f_SPdpPtr =  &(dFdx[li_SourcePrime][SPdp ]);	  q_SPdpPtr =  &(dQdx[li_SourcePrime][SPdp ]);
  f_SPspPtr =  &(dFdx[li_SourcePrime][SPsp ]);	  q_SPspPtr =  &(dQdx[li_SourcePrime][SPsp ]);
  f_SPgpPtr =  &(dFdx[li_SourcePrime][SPgp ]);	  q_SPgpPtr =  &(dQdx[li_SourcePrime][SPgp ]);
  f_SPgmPtr =  &(dFdx[li_SourcePrime][SPgm ]);	  q_SPgmPtr =  &(dQdx[li_SourcePrime][SPgm ]);
  f_SPbpPtr =  &(dFdx[li_SourcePrime][SPbp ]);	  q_SPbpPtr =  &(dQdx[li_SourcePrime][SPbp ]);
  f_SPsbPtr =  &(dFdx[li_SourcePrime][SPsb ]);	  q_SPsbPtr =  &(dQdx[li_SourcePrime][SPsb ]);

  if( trnqsMod )
  {
    f_SPqPtr = &(dFdx[li_SourcePrime][SPq ]);	    q_SPqPtr = &(dQdx[li_SourcePrime][SPq ]);
  }

  f_GPgePtr = &(dFdx[li_GatePrime][GPge ]);	  q_GPgePtr = &(dQdx[li_GatePrime][GPge ]);
  f_GPdpPtr = &(dFdx[li_GatePrime][GPdp ]);	  q_GPdpPtr = &(dQdx[li_GatePrime][GPdp ]);
  f_GPspPtr = &(dFdx[li_GatePrime][GPsp ]);	  q_GPspPtr = &(dQdx[li_GatePrime][GPsp ]);
  f_GPgpPtr = &(dFdx[li_GatePrime][GPgp ]);	  q_GPgpPtr = &(dQdx[li_GatePrime][GPgp ]);
  f_GPgmPtr = &(dFdx[li_GatePrime][GPgm ]);	  q_GPgmPtr = &(dQdx[li_GatePrime][GPgm ]);
  f_GPbpPtr = &(dFdx[li_GatePrime][GPbp ]);	  q_GPbpPtr = &(dQdx[li_GatePrime][GPbp ]);

  if( trnqsMod )
  {
    f_GPqPtr = &(dFdx[li_GatePrime][GPq ]);	    q_GPqPtr = &(dQdx[li_GatePrime][GPq ]);
  }

  f_GMgePtr = &(dFdx[li_GateMid][GMge ]);	  q_GMgePtr = &(dQdx[li_GateMid][GMge ]);
  f_GMdpPtr = &(dFdx[li_GateMid][GMdp ]);	  q_GMdpPtr = &(dQdx[li_GateMid][GMdp ]);
  f_GMspPtr = &(dFdx[li_GateMid][GMsp ]);	  q_GMspPtr = &(dQdx[li_GateMid][GMsp ]);
  f_GMgpPtr = &(dFdx[li_GateMid][GMgp ]);	  q_GMgpPtr = &(dQdx[li_GateMid][GMgp ]);
  f_GMgmPtr = &(dFdx[li_GateMid][GMgm ]);	  q_GMgmPtr = &(dQdx[li_GateMid][GMgm ]);
  f_GMbpPtr = &(dFdx[li_GateMid][GMbp ]);	  q_GMbpPtr = &(dQdx[li_GateMid][GMbp ]);

  f_BPbPtr  = &(dFdx[li_BodyPrime][BPb  ]);	  q_BPbPtr  = &(dQdx[li_BodyPrime][BPb  ]);
  f_BPdpPtr = &(dFdx[li_BodyPrime][BPdp ]);	  q_BPdpPtr = &(dQdx[li_BodyPrime][BPdp ]);
  f_BPspPtr = &(dFdx[li_BodyPrime][BPsp ]);	  q_BPspPtr = &(dQdx[li_BodyPrime][BPsp ]);
  f_BPgpPtr = &(dFdx[li_BodyPrime][BPgp ]);	  q_BPgpPtr = &(dQdx[li_BodyPrime][BPgp ]);
  f_BPgmPtr = &(dFdx[li_BodyPrime][BPgm ]);	  q_BPgmPtr = &(dQdx[li_BodyPrime][BPgm ]);
  f_BPbpPtr = &(dFdx[li_BodyPrime][BPbp ]);	  q_BPbpPtr = &(dQdx[li_BodyPrime][BPbp ]);
  f_BPsbPtr = &(dFdx[li_BodyPrime][BPsb ]);	  q_BPsbPtr = &(dQdx[li_BodyPrime][BPsb ]);
  f_BPdbPtr = &(dFdx[li_BodyPrime][BPdb ]);	  q_BPdbPtr = &(dQdx[li_BodyPrime][BPdb ]);

  f_SBbPtr  = &(dFdx[li_SourceBody][SBb  ]);	  q_SBbPtr  = &(dQdx[li_SourceBody][SBb  ]);
  f_SBspPtr = &(dFdx[li_SourceBody][SBsp ]);	  q_SBspPtr = &(dQdx[li_SourceBody][SBsp ]);
  f_SBbpPtr = &(dFdx[li_SourceBody][SBbp ]);	  q_SBbpPtr = &(dQdx[li_SourceBody][SBbp ]);
  f_SBsbPtr = &(dFdx[li_SourceBody][SBsb ]);	  q_SBsbPtr = &(dQdx[li_SourceBody][SBsb ]);

  f_DBbPtr  = &(dFdx[li_DrainBody][DBb  ]);	  q_DBbPtr  = &(dQdx[li_DrainBody][DBb  ]);
  f_DBdpPtr = &(dFdx[li_DrainBody][DBdp ]);	  q_DBdpPtr = &(dQdx[li_DrainBody][DBdp ]);
  f_DBbpPtr = &(dFdx[li_DrainBody][DBbp ]);	  q_DBbpPtr = &(dQdx[li_DrainBody][DBbp ]);
  f_DBdbPtr = &(dFdx[li_DrainBody][DBdb ]);	  q_DBdbPtr = &(dQdx[li_DrainBody][DBdb ]);

  // Optional nqs row:
  if( trnqsMod )
  {
    f_QdpPtr = &(dFdx[li_Charge][Qdp ]);	    q_QdpPtr = &(dQdx[li_Charge][Qdp ]);
    f_QspPtr = &(dFdx[li_Charge][Qsp ]);	    q_QspPtr = &(dQdx[li_Charge][Qsp ]);
    f_QgpPtr = &(dFdx[li_Charge][Qgp ]);	    q_QgpPtr = &(dQdx[li_Charge][Qgp ]);
    f_QbpPtr = &(dFdx[li_Charge][Qbp ]);	    q_QbpPtr = &(dQdx[li_Charge][Qbp ]);
    f_QqPtr  = &(dFdx[li_Charge][Qq  ]);	    q_QqPtr  = &(dQdx[li_Charge][Qq  ]);
  }
  if (icVBSGiven)
  {
    f_IBSbPtr = &(dFdx[li_Ibs][IBSb]);
    f_IBSsPtr = &(dFdx[li_Ibs][IBSs]);
    f_IBSibsPtr = &(dFdx[li_Ibs][IBSibs]);
  }
  if (icVDSGiven)
  {
    f_IDSdPtr = &(dFdx[li_Ids][IDSd]);
    f_IDSsPtr = &(dFdx[li_Ids][IDSs]);
    f_IDSidsPtr = &(dFdx[li_Ids][IDSids]);
  }
  if (icVGSGiven)
  {
    f_IGSgPtr = &(dFdx[li_Igs][IGSg]);
    f_IGSsPtr = &(dFdx[li_Igs][IGSs]);
    f_IGSigsPtr = &(dFdx[li_Igs][IGSigs]);
  }

#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       : This updates all the instance-owned paramters which
//                 are temperature dependent.
//
// Special Notes : Annoyingly, some model-owned parameters need to be
//                 tweaked here because of how the SPICE code is set up.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
bool Instance::updateTemperature (const double & temp_tmp)
{
  std::string msg="";

  double tmp(0.0), tmp1(0.0), tmp2(0.0), tmp3(0.0), Eg(0.0), Eg0(0.0), ni,epssub;
  double T0(0.0), T1(0.0);
  double T2(0.0), T3(0.0), T4(0.0), T5(0.0), T6(0.0), T7(0.0), T8(0.0), T9(0.0), Lnew(0.0), Wnew(0.0);
  double delTemp(0.0), TRatio(0.0), Inv_L(0.0), Inv_W(0.0), Inv_LW(0.0), Vtm0, Tnom(0.0);
  double dumPs(0.0), dumPd(0.0), dumAs(0.0), dumAd(0.0), PowWeffWr(0.0);
  double Nvtms(0.0), Nvtmd(0.0), SourceSatCurrent(0.0), DrainSatCurrent(0.0);
  double T10(0.0);
  double Inv_saref(0.0), Inv_sbref(0.0), Inv_sa(0.0), Inv_sb(0.0), rho(0.0), Ldrn(0.0), dvth0_lod(0.0);
  double W_tmp(0.0), Inv_ODeff(0.0), OD_offset(0.0), dk2_lod(0.0), deta0_lod(0.0);
  double lnl(0.0), lnw(0.0), lnnf(0.0), rbpbx(0.0), rbpby(0.0), rbsbx(0.0), rbsby(0.0), rbdbx(0.0), rbdby(0.0),bodymode(0.0);
  double kvsat(0.0), wlod(0.0), sceff(0.0), Wdrn(0.0);
  double V0, lt1, ltw, Theta0, Delt_vth, TempRatio, Vth_NarrowW, Lpe_Vb; 
  //double Vth; // converted to instance variable
  double n, Vgsteff, Vgs_eff, toxpf, toxpi, Tcen, toxe, epsrox, vddeot;

  int niter;

  bool bsuccess = true;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << std::endl << subsection_divider << std::endl;
    Xyce::dout() << "Instance::updateTemperature\n";
    Xyce::dout() << "name = " << getName() << std::endl;
  }

  // first set the instance temperature to the new temperature:
  if (temp_tmp != -999.0) temp = temp_tmp;

  Tnom = model_.tnom;
  TRatio = temp/Tnom;

  ///////////////////////////////////////////////////////////////////////////////
  // Model-specific stuff:
  // This is kludgey - the model-specific stuff should be handled in a model function.
  // Some of this used to be in the model class's processParams, but was
  // moved back here because it makes updating from new spice BSIM4 code
  // less painful, even though it makes no sense here.
  //

  if (model_.mtrlMod == 0)
  {
    if ((model_.toxeGiven) && (model_.toxpGiven) &&
        (model_.dtoxGiven) &&
        (model_.toxe != (model_.toxp +model_.dtox)))
    {
      UserWarning(*this) << "toxe, toxp and dtox all given and toxe != toxp + dtox; dtox ignored";
    }
    else if ((model_.toxeGiven) && (!model_.toxpGiven))
    {
      model_.toxp = model_.toxe - model_.dtox;
    }
    else if ((!model_.toxeGiven) && (model_.toxpGiven))
    {
      model_.toxe = model_.toxp + model_.dtox;
    }
  }

  if (model_.mtrlMod)
  {
    epsrox = 3.9;
    toxe = model_.eot;
    epssub = CONSTEPS0 * model_.epsrsub;
  }
  else
  {
    epsrox = model_.epsrox;
    toxe = model_.toxe;
    epssub = CONSTEPSSI;
  }


  model_.coxe = epsrox * CONSTEPS0 / toxe;
  if (model_.mtrlMod == 0)
    model_.coxp = model_.epsrox * CONSTEPS0 / model_.toxp;

  if (!model_.cgdoGiven)
  {
    if (model_.dlcGiven && (model_.dlc > 0.0))
    {
      model_.cgdo = model_.dlc * model_.coxe - model_.cgdl ;
    }
    else
    {
      model_.cgdo = 0.6 * model_.xj * model_.coxe;
    }
  }

  if (!model_.cgsoGiven)
  {
    if (model_.dlcGiven && (model_.dlc > 0.0))
    {
      model_.cgso = model_.dlc * model_.coxe - model_.cgsl ;
    }
    else
    {
      model_.cgso = 0.6 * model_.xj * model_.coxe;
    }
  }
  if (!model_.cgboGiven)
  {
    model_.cgbo = 2.0 * model_.dwc * model_.coxe;
  }

  model_.vcrit = CONSTvt0 * log(CONSTvt0 / (CONSTroot2 * 1.0e-14));
  model_.factor1 = sqrt(epssub / (epsrox * CONSTEPS0) * toxe);

  Vtm0 = model_.vtm0 = CONSTKboQ * Tnom;

  if (model_.mtrlMod == 0)
  {
    Eg0 = 1.16 - 7.02e-4 * Tnom * Tnom / (Tnom + 1108.0);
    ni = 1.45e10 * (Tnom / 300.15) * sqrt(Tnom / 300.15)
      * exp(21.5565981 - Eg0 / (2.0 * Vtm0));
  }
  else
  {
    Eg0 = model_.bg0sub - model_.tbgasub * Tnom * Tnom
      / (Tnom + model_.tbgbsub);
    T0 =  model_.bg0sub - model_.tbgasub * 90090.0225
      / (300.15 + model_.tbgbsub);
    ni = model_.ni0sub * (Tnom / 300.15) * sqrt(Tnom / 300.15)
      * exp((T0 - Eg0) / (2.0 * Vtm0));
  }

  model_.Eg0 = Eg0;
  model_.vtm = CONSTKboQ * temp;
  if (model_.mtrlMod == 0)
  {
    Eg = 1.16 - 7.02e-4 * temp * temp / (temp + 1108.0);
  }
  else
  {
    Eg = model_.bg0sub - model_.tbgasub * temp * temp
      / (temp + model_.tbgbsub);
  }

  if (temp != Tnom)
  {
    T0 = Eg0 / Vtm0 - Eg / model_.vtm;
    T1 = log(temp / Tnom);
    T2 = T0 + model_.SjctTempExponent * T1;
    T3 = exp(T2 / model_.SjctEmissionCoeff);
    model_.SjctTempSatCurDensity = model_.SjctSatCurDensity * T3;

    model_.SjctSidewallTempSatCurDensity
         = model_.SjctSidewallSatCurDensity * T3;

    model_.SjctGateSidewallTempSatCurDensity
                     = model_.SjctGateSidewallSatCurDensity * T3;

    T2 = T0 + model_.DjctTempExponent * T1;
    T3 = exp(T2 / model_.DjctEmissionCoeff);

    model_.DjctTempSatCurDensity = model_.DjctSatCurDensity * T3;

    model_.DjctSidewallTempSatCurDensity
           = model_.DjctSidewallSatCurDensity * T3;

    model_.DjctGateSidewallTempSatCurDensity
           = model_.DjctGateSidewallSatCurDensity * T3;

  }
  else
  {
   model_.SjctTempSatCurDensity = model_.SjctSatCurDensity;
   model_.SjctSidewallTempSatCurDensity = model_.SjctSidewallSatCurDensity;

   model_.SjctGateSidewallTempSatCurDensity
                      = model_.SjctGateSidewallSatCurDensity;

   model_.DjctTempSatCurDensity = model_.DjctSatCurDensity;

   model_.DjctSidewallTempSatCurDensity
                      = model_.DjctSidewallSatCurDensity;

   model_.DjctGateSidewallTempSatCurDensity
                      = model_.DjctGateSidewallSatCurDensity;
  }

  if (model_.SjctTempSatCurDensity < 0.0)
     model_.SjctTempSatCurDensity = 0.0;

  if (model_.SjctSidewallTempSatCurDensity < 0.0)
     model_.SjctSidewallTempSatCurDensity = 0.0;

  if (model_.SjctGateSidewallTempSatCurDensity < 0.0)
     model_.SjctGateSidewallTempSatCurDensity = 0.0;

  if (model_.DjctTempSatCurDensity < 0.0)
     model_.DjctTempSatCurDensity = 0.0;

  if (model_.DjctSidewallTempSatCurDensity < 0.0)
     model_.DjctSidewallTempSatCurDensity = 0.0;

  if (model_.DjctGateSidewallTempSatCurDensity < 0.0)
     model_.DjctGateSidewallTempSatCurDensity = 0.0;

  // Temperature dependence of D/B and S/B diode capacitance begins
  delTemp = temp - Tnom;
  T0 = model_.tcj * delTemp;
  if (T0 >= -1.0)
  {   model_.SunitAreaTempJctCap = model_.SunitAreaJctCap *(1.0 + T0); //bug_fix -JX
     model_.DunitAreaTempJctCap = model_.DunitAreaJctCap *(1.0 + T0);
  }
  else
  {
   if (model_.SunitAreaJctCap > 0.0)
   {
     model_.SunitAreaTempJctCap = 0.0;
     UserWarning(*this) << "Temperature effect has caused cjs to be negative. Cjs is clamped to zero";

   }
   if (model_.DunitAreaJctCap > 0.0)
   {
     model_.DunitAreaTempJctCap = 0.0;
     UserWarning(*this) << "Temperature effect has caused cjd to be negative. Cjd is clamped to zero";
   }
  }

  T0 = model_.tcjsw * delTemp;
  if (T0 >= -1.0)
  {
   model_.SunitLengthSidewallTempJctCap = model_.SunitLengthSidewallJctCap *(1.0 + T0);
   model_.DunitLengthSidewallTempJctCap = model_.DunitLengthSidewallJctCap *(1.0 + T0);
  }
  else
  {
   if (model_.SunitLengthSidewallJctCap > 0.0)
   {
     model_.SunitLengthSidewallTempJctCap = 0.0;
     UserWarning(*this) << "Temperature effect has caused cjsws to be negative. Cjsws is clamped to zero";
   }
   if (model_.DunitLengthSidewallJctCap > 0.0)
   {
     model_.DunitLengthSidewallTempJctCap = 0.0;
     UserWarning(*this) << "Temperature effect has caused cjswd to be negative. Cjswd is clamped to zero";
   }
  }

  T0 = model_.tcjswg * delTemp;
  if (T0 >= -1.0)
  {
   model_.SunitLengthGateSidewallTempJctCap = model_.SunitLengthGateSidewallJctCap *(1.0 + T0);
   model_.DunitLengthGateSidewallTempJctCap = model_.DunitLengthGateSidewallJctCap *(1.0 + T0);
  }
  else
  {
   if (model_.SunitLengthGateSidewallJctCap > 0.0)
   {
     model_.SunitLengthGateSidewallTempJctCap = 0.0;
     UserWarning(*this) << "Temperature effect has caused cjswgs to be negative. Cjswgs is clamped to zero";
   }
   if (model_.DunitLengthGateSidewallJctCap > 0.0)
   {
     model_.DunitLengthGateSidewallTempJctCap = 0.0;
     UserWarning(*this) << "Temperature effect has caused cjswgd to be negative. Cjswgd is clamped to zero";
   }
  }

  model_.PhiBS = model_.SbulkJctPotential - model_.tpb * delTemp;

  if (model_.PhiBS < 0.01)
  {
   model_.PhiBS = 0.01;
   UserWarning(*this) << "Temperature effect has caused pbs to be less than 0.01. Pbs is clamped to 0.01";
  }

  model_.PhiBD = model_.DbulkJctPotential - model_.tpb * delTemp;
  if (model_.PhiBD < 0.01)
  {
   model_.PhiBD = 0.01;
   UserWarning(*this) << "Temperature effect has caused pbd to be less than 0.01. Pbd is clamped to 0.01";
  }

  model_.PhiBSWS = model_.SsidewallJctPotential - model_.tpbsw * delTemp;
  if (model_.PhiBSWS <= 0.01)
  {
   model_.PhiBSWS = 0.01;
   UserWarning(*this) << "Temperature effect has caused pbsws to be less than 0.01. Pbsws is clamped to 0.01";
  }

  model_.PhiBSWD = model_.DsidewallJctPotential - model_.tpbsw * delTemp;
  if (model_.PhiBSWD <= 0.01)
  {
   model_.PhiBSWD = 0.01;
   UserWarning(*this) << "Temperature effect has caused pbswd to be less than 0.01. Pbswd is clamped to 0.01";
  }

  model_.PhiBSWGS = model_.SGatesidewallJctPotential - model_.tpbswg * delTemp;
  if (model_.PhiBSWGS <= 0.01)
  {
   model_.PhiBSWGS = 0.01;
   UserWarning(*this) << "Temperature effect has caused pbswgs to be less than 0.01. Pbswgs is clamped to 0.01";
  }

  model_.PhiBSWGD = model_.DGatesidewallJctPotential - model_.tpbswg * delTemp;
  if (model_.PhiBSWGD <= 0.01)
  {
   model_.PhiBSWGD = 0.01;
   UserWarning(*this) << "Temperature effect has caused pbswgd to be less than 0.01. Pbswgd is clamped to 0.01";
  } // End of junction capacitance


  if (model_.ijthdfwd <= 0.0)
  {
   model_.ijthdfwd = 0.1;
   UserWarning(*this) << "Ijthdfwd reset to " << model_.ijthdfwd;
  }
  if (model_.ijthsfwd <= 0.0)
  {
   model_.ijthsfwd = 0.1;
   UserWarning(*this) << "Ijthsfwd reset to " << model_.ijthsfwd;
  }
  if (model_.ijthdrev <= 0.0)
  {
   model_.ijthdrev = 0.1;
   UserWarning(*this) << "Ijthdrev reset to " << model_.ijthdrev;
  }
  if (model_.ijthsrev <= 0.0)
  {
   model_.ijthsrev = 0.1;
   UserWarning(*this) << "Ijthsrev reset to " << model_.ijthsrev;
  }

  if ((model_.xjbvd <= 0.0) && (model_.dioMod == 2))
  {
   model_.xjbvd = 1.0;
   UserWarning(*this) << "Xjbvd reset to " << model_.xjbvd;
  }
  else if ((model_.xjbvd < 0.0) && (model_.dioMod == 0))
  {
   model_.xjbvd = 1.0;
   UserWarning(*this) << "Xjbvd reset to " << model_.xjbvd;
  }

  if (model_.bvd <= 0.0)
  {
   model_.bvd = 10.0;
   UserWarning(*this) << "BVD reset to " << model_.bvd;
  }

  if ((model_.xjbvs <= 0.0) && (model_.dioMod == 2))
  {
   model_.xjbvs = 1.0;
   UserWarning(*this) << "Xjbvs reset to " << model_.xjbvs;
  }
  else if ((model_.xjbvs < 0.0) && (model_.dioMod == 0))
  {
   model_.xjbvs = 1.0;
   UserWarning(*this) << "Xjbvs reset to " << model_.xjbvs;
  }

  if (model_.bvs <= 0.0)
  {
   model_.bvs = 10.0;
   UserWarning(*this) << "BVS reset to " << model_.bvs;
  }


  ///////////////////////////////////////////////////////////////////////////////
  // Instance stuff:
  // (loop through all the instances of the model)

  // stress effect
  Ldrn = l;
  Wdrn = w / nf;

  // This next block determines whether  or not to use a previously allocated
  // set of size dependent parameters. These are stored in a list that is
  // owned by the model.  If the values for length and width match those of
  // a previously allocated set, then use the old set.  If not, allocate a new set.

  std::list<SizeDependParam*>::iterator it_dpL =
    model_.sizeDependParamList.begin();
  std::list<SizeDependParam*>::iterator end_dpL =
    model_.sizeDependParamList.end();

  paramPtr = NULL;

  for( ; it_dpL != end_dpL; ++it_dpL )
  {
    if( ((*it_dpL)->Length  == l)
     && ((*it_dpL)->Width   == w)
     && ((*it_dpL)->NFinger == nf) 
     && ((*it_dpL)->referenceTemperature == temp_tmp))
    {
      paramPtr = (*it_dpL);
    }
  }

  // This was inside of the "Size_Not_Found" if-statement, but that
  // won't work here - it winds up being uninitialized whenever the
  // size pointer is found
  Lnew = l  + model_.xl ;
  Wnew = w / nf + model_.xw;

  if ( paramPtr != NULL )
  {
  }
  else
  {
    paramPtr = new SizeDependParam ();

    model_.sizeDependParamList.push_back( paramPtr );
    paramPtr->referenceTemperature = temp_tmp;

    //paramPtr->pNext = NULL;

    paramPtr->Length = l;
    paramPtr->Width = w;
    paramPtr->NFinger = nf;
    //Lnew = l  + model_.xl ;
    //Wnew = w / nf + model_.xw;

    T0 = pow(Lnew, model_.Lln);
    T1 = pow(Wnew, model_.Lwn);
    tmp1 = model_.Ll / T0 + model_.Lw / T1
         + model_.Lwl / (T0 * T1);
    paramPtr->dl = model_.Lint + tmp1;
    tmp2 = model_.Llc / T0 + model_.Lwc / T1
         + model_.Lwlc / (T0 * T1);
    paramPtr->dlc = model_.dlc + tmp2;

    T2 = pow(Lnew, model_.Wln);
    T3 = pow(Wnew, model_.Wwn);
    tmp1 = model_.Wl / T2 + model_.Ww / T3
         + model_.Wwl / (T2 * T3);
    paramPtr->dw = model_.Wint + tmp1;
    tmp2 = model_.Wlc / T2 + model_.Wwc / T3
         + model_.Wwlc / (T2 * T3);
    paramPtr->dwc = model_.dwc + tmp2;
    paramPtr->dwj = model_.dwj + tmp2;

    paramPtr->leff = Lnew - 2.0 * paramPtr->dl;

    if (paramPtr->leff <= 0.0)
    {
      UserError(*this) << "mosfet " << getName() << " model " << model_.getName()
                        << "  Effective channel length <= 0";
    }

    paramPtr->weff = Wnew - 2.0 * paramPtr->dw;
    if (paramPtr->weff <= 0.0)
    {
      UserError(*this) << "mosfet " << getName() << " model " << model_.getName()
                        << "  Effective channel width <= 0";
    }

    paramPtr->leffCV = Lnew - 2.0 * paramPtr->dlc;
    if (paramPtr->leffCV <= 0.0)
    {
      UserError(*this) << "mosfet " << getName() << " model " << model_.getName()
                        << "  Effective channel length for C-V <= 0";
    }

    paramPtr->weffCV = Wnew - 2.0 * paramPtr->dwc;
    if (paramPtr->weffCV <= 0.0)
    {
      UserError(*this) << "mosfet " << getName() << " model " << model_.getName()
                        << "  Effective channel width for C-V <= 0";
    }

    paramPtr->weffCJ = Wnew - 2.0 * paramPtr->dwj;
    if (paramPtr->weffCJ <= 0.0)
    {
      UserError(*this) << "mosfet " << getName() << " model " << model_.getName()
                        << "  Effective channel width for S/D junctions <= 0";
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
    paramPtr->ndep = model_.ndep
                  + model_.lndep * Inv_L
                  + model_.wndep * Inv_W
                  + model_.pndep * Inv_LW;
    paramPtr->nsd = model_.nsd
                     + model_.lnsd * Inv_L
                     + model_.wnsd * Inv_W
                     + model_.pnsd * Inv_LW;
    paramPtr->phin = model_.phin
                      + model_.lphin * Inv_L
                      + model_.wphin * Inv_W
                      + model_.pphin * Inv_LW;
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
    paramPtr->lpe0 = model_.lpe0
                      + model_.llpe0 * Inv_L
                      + model_.wlpe0 * Inv_W
                      + model_.plpe0 * Inv_LW;
    paramPtr->lpeb = model_.lpeb
                      + model_.llpeb * Inv_L
                      + model_.wlpeb * Inv_W
                      + model_.plpeb * Inv_LW;
    paramPtr->dvtp0 = model_.dvtp0
                       + model_.ldvtp0 * Inv_L
                       + model_.wdvtp0 * Inv_W
                       + model_.pdvtp0 * Inv_LW;
    paramPtr->dvtp1 = model_.dvtp1
                       + model_.ldvtp1 * Inv_L
                       + model_.wdvtp1 * Inv_W
                       + model_.pdvtp1 * Inv_LW;
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
    paramPtr->ud = model_.ud
                      + model_.lud * Inv_L
                      + model_.wud * Inv_W
                      + model_.pud * Inv_LW;
    paramPtr->ud1 = model_.ud1
                      + model_.lud1 * Inv_L
                      + model_.wud1 * Inv_W
                      + model_.pud1 * Inv_LW;
    paramPtr->up = model_.up
                      + model_.lup * Inv_L
                      + model_.wup * Inv_W
                      + model_.pup * Inv_LW;
    paramPtr->lp = model_.lp
                      + model_.llp * Inv_L
                      + model_.wlp * Inv_W
                      + model_.plp * Inv_LW;
    paramPtr->eu = model_.eu
                    + model_.leu * Inv_L
                    + model_.weu * Inv_W
                    + model_.peu * Inv_LW;
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
    paramPtr->tvoff = model_.tvoff
                    + model_.ltvoff * Inv_L
                    + model_.wtvoff * Inv_W
                    + model_.ptvoff * Inv_LW;
    paramPtr->minv = model_.minv
                      + model_.lminv * Inv_L
                      + model_.wminv * Inv_W
                      + model_.pminv * Inv_LW;
    paramPtr->minvcv = model_.minvcv
                      + model_.lminvcv * Inv_L
                      + model_.wminvcv * Inv_W
                      + model_.pminvcv * Inv_LW;
    paramPtr->fprout = model_.fprout
                       + model_.lfprout * Inv_L
                       + model_.wfprout * Inv_W
                       + model_.pfprout * Inv_LW;
    paramPtr->pdits = model_.pdits
                       + model_.lpdits * Inv_L
                       + model_.wpdits * Inv_W
                       + model_.ppdits * Inv_LW;
    paramPtr->pditsd = model_.pditsd
                        + model_.lpditsd * Inv_L
                        + model_.wpditsd * Inv_W
                        + model_.ppditsd * Inv_LW;
    paramPtr->delta = model_.delta
                        + model_.ldelta * Inv_L
                        + model_.wdelta * Inv_W
                        + model_.pdelta * Inv_LW;
    paramPtr->rdsw = model_.rdsw
                        + model_.lrdsw * Inv_L
                        + model_.wrdsw * Inv_W
                        + model_.prdsw * Inv_LW;
    paramPtr->rdw = model_.rdw
                      + model_.lrdw * Inv_L
                      + model_.wrdw * Inv_W
                      + model_.prdw * Inv_LW;
    paramPtr->rsw = model_.rsw
                      + model_.lrsw * Inv_L
                      + model_.wrsw * Inv_W
                      + model_.prsw * Inv_LW;
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
    paramPtr->agidl = model_.agidl
                       + model_.lagidl * Inv_L
                       + model_.wagidl * Inv_W
                       + model_.pagidl * Inv_LW;
    paramPtr->bgidl = model_.bgidl
                       + model_.lbgidl * Inv_L
                       + model_.wbgidl * Inv_W
                       + model_.pbgidl * Inv_LW;
    paramPtr->cgidl = model_.cgidl
                       + model_.lcgidl * Inv_L
                       + model_.wcgidl * Inv_W
                       + model_.pcgidl * Inv_LW;
    paramPtr->egidl = model_.egidl
                       + model_.legidl * Inv_L
                       + model_.wegidl * Inv_W
                       + model_.pegidl * Inv_LW;
    paramPtr->agisl = model_.agisl
                         + model_.lagisl * Inv_L
                         + model_.wagisl * Inv_W
                         + model_.pagisl * Inv_LW;
    paramPtr->bgisl = model_.bgisl
                         + model_.lbgisl * Inv_L
                         + model_.wbgisl * Inv_W
                         + model_.pbgisl * Inv_LW;
    paramPtr->cgisl = model_.cgisl
                         + model_.lcgisl * Inv_L
                         + model_.wcgisl * Inv_W
                         + model_.pcgisl * Inv_LW;
    paramPtr->egisl = model_.egisl
                         + model_.legisl * Inv_L
                         + model_.wegisl * Inv_W
                         + model_.pegisl * Inv_LW;
    paramPtr->aigc = model_.aigc
                       + model_.laigc * Inv_L
                       + model_.waigc * Inv_W
                       + model_.paigc * Inv_LW;
    paramPtr->bigc = model_.bigc
                       + model_.lbigc * Inv_L
                       + model_.wbigc * Inv_W
                       + model_.pbigc * Inv_LW;
    paramPtr->cigc = model_.cigc
                       + model_.lcigc * Inv_L
                       + model_.wcigc * Inv_W
                       + model_.pcigc * Inv_LW;
    paramPtr->aigs = model_.aigs
                       + model_.laigs * Inv_L
                       + model_.waigs * Inv_W
                       + model_.paigs * Inv_LW;
    paramPtr->bigs = model_.bigs
                       + model_.lbigs * Inv_L
                       + model_.wbigs * Inv_W
                       + model_.pbigs * Inv_LW;
    paramPtr->cigs = model_.cigs
                       + model_.lcigs * Inv_L
                       + model_.wcigs * Inv_W
                       + model_.pcigs * Inv_LW;
    paramPtr->aigd = model_.aigd
                       + model_.laigd * Inv_L
                       + model_.waigd * Inv_W
                       + model_.paigd * Inv_LW;
    paramPtr->bigd = model_.bigd
                       + model_.lbigd * Inv_L
                       + model_.wbigd * Inv_W
                       + model_.pbigd * Inv_LW;
    paramPtr->cigd = model_.cigd
                       + model_.lcigd * Inv_L
                       + model_.wcigd * Inv_W
                       + model_.pcigd * Inv_LW;
    paramPtr->aigbacc = model_.aigbacc
                         + model_.laigbacc * Inv_L
                         + model_.waigbacc * Inv_W
                         + model_.paigbacc * Inv_LW;
    paramPtr->bigbacc = model_.bigbacc
                         + model_.lbigbacc * Inv_L
                         + model_.wbigbacc * Inv_W
                         + model_.pbigbacc * Inv_LW;
    paramPtr->cigbacc = model_.cigbacc
                         + model_.lcigbacc * Inv_L
                         + model_.wcigbacc * Inv_W
                         + model_.pcigbacc * Inv_LW;
    paramPtr->aigbinv = model_.aigbinv
                         + model_.laigbinv * Inv_L
                         + model_.waigbinv * Inv_W
                         + model_.paigbinv * Inv_LW;
    paramPtr->bigbinv = model_.bigbinv
                         + model_.lbigbinv * Inv_L
                         + model_.wbigbinv * Inv_W
                         + model_.pbigbinv * Inv_LW;
    paramPtr->cigbinv = model_.cigbinv
                         + model_.lcigbinv * Inv_L
                         + model_.wcigbinv * Inv_W
                         + model_.pcigbinv * Inv_LW;
    paramPtr->nigc = model_.nigc
                         + model_.lnigc * Inv_L
                         + model_.wnigc * Inv_W
                         + model_.pnigc * Inv_LW;
    paramPtr->nigbacc = model_.nigbacc
                         + model_.lnigbacc * Inv_L
                         + model_.wnigbacc * Inv_W
                         + model_.pnigbacc * Inv_LW;
    paramPtr->nigbinv = model_.nigbinv
                         + model_.lnigbinv * Inv_L
                         + model_.wnigbinv * Inv_W
                         + model_.pnigbinv * Inv_LW;
    paramPtr->ntox = model_.ntox
                      + model_.lntox * Inv_L
                      + model_.wntox * Inv_W
                      + model_.pntox * Inv_LW;
    paramPtr->eigbinv = model_.eigbinv
                         + model_.leigbinv * Inv_L
                         + model_.weigbinv * Inv_W
                         + model_.peigbinv * Inv_LW;
    paramPtr->pigcd = model_.pigcd
                       + model_.lpigcd * Inv_L
                       + model_.wpigcd * Inv_W
                       + model_.ppigcd * Inv_LW;
    paramPtr->poxedge = model_.poxedge
                         + model_.lpoxedge * Inv_L
                         + model_.wpoxedge * Inv_W
                         + model_.ppoxedge * Inv_LW;
    paramPtr->xrcrg1 = model_.xrcrg1
                        + model_.lxrcrg1 * Inv_L
                        + model_.wxrcrg1 * Inv_W
                        + model_.pxrcrg1 * Inv_LW;
    paramPtr->xrcrg2 = model_.xrcrg2
                        + model_.lxrcrg2 * Inv_L
                        + model_.wxrcrg2 * Inv_W
                        + model_.pxrcrg2 * Inv_LW;
    paramPtr->lambda = model_.lambda
                        + model_.llambda * Inv_L
                        + model_.wlambda * Inv_W
                        + model_.plambda * Inv_LW;
    paramPtr->vtl = model_.vtl
                        + model_.lvtl * Inv_L
                        + model_.wvtl * Inv_W
                        + model_.pvtl * Inv_LW;
    paramPtr->xn = model_.xn
                        + model_.lxn * Inv_L
                        + model_.wxn * Inv_W
                        + model_.pxn * Inv_LW;
    paramPtr->vfbsdoff = model_.vfbsdoff
                        + model_.lvfbsdoff * Inv_L
                        + model_.wvfbsdoff * Inv_W
                        + model_.pvfbsdoff * Inv_LW;
    paramPtr->tvfbsdoff = model_.tvfbsdoff
                        + model_.ltvfbsdoff * Inv_L
                        + model_.wtvfbsdoff * Inv_W
                        + model_.ptvfbsdoff * Inv_LW;

    paramPtr->cgsl = model_.cgsl
                        + model_.lcgsl * Inv_L
                        + model_.wcgsl * Inv_W
                        + model_.pcgsl * Inv_LW;
    paramPtr->cgdl = model_.cgdl
                        + model_.lcgdl * Inv_L
                        + model_.wcgdl * Inv_W
                        + model_.pcgdl * Inv_LW;
    paramPtr->ckappas = model_.ckappas
                        + model_.lckappas * Inv_L
                        + model_.wckappas * Inv_W
                        + model_.pckappas * Inv_LW;
    paramPtr->ckappad = model_.ckappad
                         + model_.lckappad * Inv_L
                         + model_.wckappad * Inv_W
                         + model_.pckappad * Inv_LW;
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
    paramPtr->kvth0we = model_.kvth0we
                        + model_.lkvth0we * Inv_L
                        + model_.wkvth0we * Inv_W
                        + model_.pkvth0we * Inv_LW;
    paramPtr->k2we = model_.k2we
                        + model_.lk2we * Inv_L
                        + model_.wk2we * Inv_W
                        + model_.pk2we * Inv_LW;
    paramPtr->ku0we = model_.ku0we
                        + model_.lku0we * Inv_L
                        + model_.wku0we * Inv_W
                        + model_.pku0we * Inv_LW;

    paramPtr->abulkCVfactor = 1.0 + pow((paramPtr->clc / paramPtr->leffCV), paramPtr->cle);

    T0 = (TRatio - 1.0);

    PowWeffWr = pow(paramPtr->weffCJ * 1.0e6, paramPtr->wr) * nf;

    T1 = T2 = T3 = T4 = 0.0;
    if(model_.tempMod == 0)
    {
      paramPtr->ua = paramPtr->ua + paramPtr->ua1 * T0;
      paramPtr->ub = paramPtr->ub + paramPtr->ub1 * T0;
      paramPtr->uc = paramPtr->uc + paramPtr->uc1 * T0;
      paramPtr->ud = paramPtr->ud + paramPtr->ud1 * T0;
        paramPtr->vsattemp = paramPtr->vsat - paramPtr->at * T0;
      T10 = paramPtr->prt * T0;
      if(model_.rdsMod)
      {
          // External Rd(V)
          T1 = paramPtr->rdw + T10;
            T2 = model_.rdwmin + T10;
          // External Rs(V)
          T3 = paramPtr->rsw + T10;
            T4 = model_.rswmin + T10;
      }
      // Internal Rds(V) in IV
      paramPtr->rds0 = (paramPtr->rdsw + T10) * nf / PowWeffWr;
      paramPtr->rdswmin = (model_.rdswmin + T10) * nf / PowWeffWr;
    }
    else
    { // tempMod = 1, 2
      paramPtr->ua = paramPtr->ua * (1.0 + paramPtr->ua1 * delTemp) ;
      paramPtr->ub = paramPtr->ub * (1.0 + paramPtr->ub1 * delTemp);
      paramPtr->uc = paramPtr->uc * (1.0 + paramPtr->uc1 * delTemp);
      paramPtr->ud = paramPtr->ud * (1.0 + paramPtr->ud1 * delTemp);
      paramPtr->vsattemp = paramPtr->vsat * (1.0 - paramPtr->at * delTemp);
      T10 = 1.0 + paramPtr->prt * delTemp;
      if(model_.rdsMod)
      {
        // External Rd(V)
        T1 = paramPtr->rdw * T10;
        T2 = model_.rdwmin * T10;
        // External Rs(V)
        T3 = paramPtr->rsw * T10;
        T4 = model_.rswmin * T10;
      }
      // Internal Rds(V) in IV
      paramPtr->rds0 = paramPtr->rdsw * T10 * nf / PowWeffWr;
      paramPtr->rdswmin = model_.rdswmin * T10 * nf / PowWeffWr;
    }
    if (T1 < 0.0)
    {
      T1 = 0.0;
      UserWarning(*this) << "Rdw at current temperature is negative; set to 0";
    }
    if (T2 < 0.0)
    {
      T2 = 0.0;
      UserWarning(*this) << "Rdwmin at current temperature is negative; set to 0";
    }
    paramPtr->rd0 = T1 / PowWeffWr;
    paramPtr->rdwmin = T2 / PowWeffWr;
    if (T3 < 0.0)
    {
      T3 = 0.0;
      UserWarning(*this) << "Rsw at current temperature is negative; set to 0";
    }
    if (T4 < 0.0)
    {
      T4 = 0.0;
      UserWarning(*this) << "Rswmin at current temperature is negative; set to 0";
    }
    paramPtr->rs0 = T3 / PowWeffWr;
    paramPtr->rswmin = T4 / PowWeffWr;

    if (paramPtr->u0 > 1.0)
        paramPtr->u0 = paramPtr->u0 / 1.0e4;

    // mobility channel length dependence
    T5 = 1.0 - paramPtr->up * exp( - paramPtr->leff / paramPtr->lp);
    paramPtr->u0temp = paramPtr->u0 * T5
    * pow(TRatio, paramPtr->ute);
    if (paramPtr->eu < 0.0)
    {
      paramPtr->eu = 0.0;
      UserWarning(*this) << "eu has been negative; reset to 0.0";
    }

    paramPtr->vfbsdoff = paramPtr->vfbsdoff * (1.0 + paramPtr->tvfbsdoff * delTemp);
    paramPtr->voff = paramPtr->voff * (1.0 + paramPtr->tvoff * delTemp);

    // Source End Velocity Limit
    if((model_.vtlGiven) && (model_.vtl > 0.0) )
    {
       if(model_.lc < 0.0) paramPtr->lc = 0.0;
       else   paramPtr->lc = model_.lc ;
       T0 = paramPtr->leff / (paramPtr->xn * paramPtr->leff + paramPtr->lc);
       paramPtr->tfactor = (1.0 - T0) / (1.0 + T0 );
    }

    paramPtr->cgdo = (model_.cgdo + paramPtr->cf) * paramPtr->weffCV;
    paramPtr->cgso = (model_.cgso + paramPtr->cf) * paramPtr->weffCV;
    paramPtr->cgbo = model_.cgbo * paramPtr->leffCV * nf;

    if (!model_.ndepGiven && model_.gamma1Given)
    {
      T0 = paramPtr->gamma1 * model_.coxe;
      paramPtr->ndep = 3.01248e22 * T0 * T0;
    }

    paramPtr->phi = Vtm0 * log(paramPtr->ndep / ni) + paramPtr->phin + 0.4;

    paramPtr->sqrtPhi = sqrt(paramPtr->phi);
    paramPtr->phis3 = paramPtr->sqrtPhi * paramPtr->phi;

    paramPtr->Xdep0 = sqrt(2.0 * epssub / (Charge_q
                        * paramPtr->ndep * 1.0e6))
                        * paramPtr->sqrtPhi;

    paramPtr->sqrtXdep0 = sqrt(paramPtr->Xdep0);

    if (model_.mtrlMod == 0)
    {
      paramPtr->litl = sqrt(3.0 * paramPtr->xj * toxe);
    }
    else
    {
      paramPtr->litl = sqrt(model_.epsrsub/epsrox * paramPtr->xj * toxe);
    }

    paramPtr->vbi = Vtm0 * log(paramPtr->nsd * paramPtr->ndep / (ni * ni));

    if (model_.mtrlMod == 0)
    {
      if (paramPtr->ngate > 0.0)
      {   paramPtr->vfbsd = Vtm0 * log(paramPtr->ngate
                                          / paramPtr->nsd);
      }
      else
        paramPtr->vfbsd = 0.0;
    }
    else
    {
      T0 = Vtm0 * log(paramPtr->nsd/ni);
      T1 = 0.5 * Eg0;
      if(T0 > T1)
        T0 = T1;
      T2 = model_.easub + T1 - model_.dtype * T0;
      paramPtr->vfbsd = model_.phig - T2;
    }

    paramPtr->cdep0 = sqrt(Charge_q * epssub
                        * paramPtr->ndep * 1.0e6 / 2.0 / paramPtr->phi);

    paramPtr->ToxRatio = exp(paramPtr->ntox
                            * log(model_.toxref / toxe))
                            / toxe / toxe;

    paramPtr->ToxRatioEdge = exp(paramPtr->ntox
                              * log(model_.toxref
                              / (toxe * paramPtr->poxedge)))
                              / toxe / toxe
                              / paramPtr->poxedge / paramPtr->poxedge;

    paramPtr->Aechvb = (model_.dtype == CONSTNMOS) ? 4.97232e-7 : 3.42537e-7;
    paramPtr->Bechvb = (model_.dtype == CONSTNMOS) ? 7.45669e11 : 1.16645e12;
    paramPtr->AechvbEdgeS = paramPtr->Aechvb * paramPtr->weff
                            * model_.dlcig * paramPtr->ToxRatioEdge;
    paramPtr->AechvbEdgeD = paramPtr->Aechvb * paramPtr->weff
                            * model_.dlcigd * paramPtr->ToxRatioEdge;
    paramPtr->BechvbEdge = -paramPtr->Bechvb
                            * model_.toxe * paramPtr->poxedge;
    paramPtr->Aechvb *= paramPtr->weff * paramPtr->leff
                          * paramPtr->ToxRatio;
    paramPtr->Bechvb *= -toxe;


    paramPtr->mstar = 0.5 + atan(paramPtr->minv) / M_PI;
    paramPtr->mstarcv = 0.5 + atan(paramPtr->minvcv) / M_PI;
    paramPtr->voffcbn =  paramPtr->voff + model_.voffl / paramPtr->leff;
    paramPtr->voffcbncv =  paramPtr->voffcv + model_.voffcvl / paramPtr->leff;

    paramPtr->ldeb = sqrt(epssub * Vtm0 / (Charge_q
                                               * paramPtr->ndep * 1.0e6)) / 3.0;
    paramPtr->acde *= pow((paramPtr->ndep / 2.0e16), -0.25);


    if (model_.k1Given || model_.k2Given)
    {
      if (!model_.k1Given)
      {
        UserWarning(*this) << "k1 should be specified with k2";
        paramPtr->k1 = 0.53;
      }
      if (!model_.k2Given)
      {
        UserWarning(*this) << "k2 should be specified with k1";
        paramPtr->k2 = -0.0186;
      }
      if (model_.nsubGiven)
      {
        UserWarning(*this) << "nsub is ignored because k1 or k2 is given";
      }
      if (model_.xtGiven)
      {
        UserWarning(*this) << "xt is ignored because k1 or k2 is given";
      }
      if (model_.vbxGiven)
      {
        UserWarning(*this) << "vbx is ignored because k1 or k2 is given";
      }
      if (model_.gamma1Given)
      {
        UserWarning(*this) << "gamma1 is ignored because k1 or k2 is given";
      }
      if (model_.gamma2Given)
      {
        UserWarning(*this) << "gamma2 is ignored because k1 or k2 is given";
      }
    }
    else
    {
      if (!model_.vbxGiven)
            paramPtr->vbx = paramPtr->phi - 7.7348e-4
                             * paramPtr->ndep * paramPtr->xt * paramPtr->xt;

      if (paramPtr->vbx > 0.0)
            paramPtr->vbx = -paramPtr->vbx;
      if (paramPtr->vbm > 0.0)
            paramPtr->vbm = -paramPtr->vbm;

      if (!model_.gamma1Given)
            paramPtr->gamma1 = 5.753e-12
                                * sqrt(paramPtr->ndep) / model_.coxe;

      if (!model_.gamma2Given)
            paramPtr->gamma2 = 5.753e-12
                              * sqrt(paramPtr->nsub) / model_.coxe;

      T0 = paramPtr->gamma1 - paramPtr->gamma2;

      T1 = sqrt(paramPtr->phi - paramPtr->vbx) - paramPtr->sqrtPhi;

      T2 = sqrt(paramPtr->phi * (paramPtr->phi - paramPtr->vbm)) - paramPtr->phi;

      paramPtr->k2 = T0 * T1 / (2.0 * T2 + paramPtr->vbm);

      paramPtr->k1 = paramPtr->gamma2 - 2.0
                    * paramPtr->k2 * sqrt(paramPtr->phi - paramPtr->vbm);
    }

    if (!model_.vfbGiven)
    {
      if (model_.vth0Given)
      {   paramPtr->vfb = model_.dtype * paramPtr->vth0
                             - paramPtr->phi - paramPtr->k1
                             * paramPtr->sqrtPhi;
      }
      else
      {
        if ((model_.mtrlMod) && (model_.phigGiven) &&
            (model_.nsubGiven))
        {
          T0 = Vtm0 * log(paramPtr->nsub/ni);
          T1 = 0.5 * Eg0;
          if(T0 > T1)
            T0 = T1;
          T2 = model_.easub + T1 + model_.dtype * T0;
          paramPtr->vfb = model_.phig - T2;
        }
        else
        {
          paramPtr->vfb = -1.0;
        }
      }
    }
    if (!model_.vth0Given)
    {
      paramPtr->vth0 = model_.dtype * (paramPtr->vfb
                          + paramPtr->phi + paramPtr->k1
                          * paramPtr->sqrtPhi);
    }

    paramPtr->k1ox = paramPtr->k1 * toxe / model_.toxm;

    tmp = sqrt(epssub / (epsrox * CONSTEPS0)
        * toxe * paramPtr->Xdep0);

    T0 = paramPtr->dsub * paramPtr->leff / tmp;
    if (T0 < CONSTEXP_THRESHOLD)
    {
      T1 = exp(T0);
      T2 = T1 - 1.0;
      T3 = T2 * T2;
      T4 = T3 + 2.0 * T1 * CONSTMIN_EXP;
      paramPtr->theta0vb0 = T1 / T4;
    }
    else
    {
      paramPtr->theta0vb0 = 1.0 / (CONSTMAX_EXP - 2.0);
    }

    T0 = paramPtr->drout * paramPtr->leff / tmp;
    if (T0 < CONSTEXP_THRESHOLD)
    {
      T1 = exp(T0);
      T2 = T1 - 1.0;
      T3 = T2 * T2;
      T4 = T3 + 2.0 * T1 * CONSTMIN_EXP;
      T5 = T1 / T4;
    }
    else
        T5 = 1.0 / (CONSTMAX_EXP - 2.0); // 3.0 * CONSTMIN_EXP omitted

    paramPtr->thetaRout = paramPtr->pdibl1 * T5 + paramPtr->pdibl2;

    tmp = sqrt(paramPtr->Xdep0);
    tmp1 = paramPtr->vbi - paramPtr->phi;
    tmp2 = model_.factor1 * tmp;

    T0 = paramPtr->dvt1w * paramPtr->weff * paramPtr->leff / tmp2;

    if (T0 < CONSTEXP_THRESHOLD)
    {
      T1 = exp(T0);
      T2 = T1 - 1.0;
      T3 = T2 * T2;
      T4 = T3 + 2.0 * T1 * CONSTMIN_EXP;
      T8 = T1 / T4;
    }
    else
    {
      T8 = 1.0 / (CONSTMAX_EXP - 2.0);
    }

    T0 = paramPtr->dvt0w * T8;
    T8 = T0 * tmp1;

    T0 = paramPtr->dvt1 * paramPtr->leff / tmp2;
    if (T0 < CONSTEXP_THRESHOLD)
    {
      T1 = exp(T0);
      T2 = T1 - 1.0;
      T3 = T2 * T2;
      T4 = T3 + 2.0 * T1 * CONSTMIN_EXP;
      T9 = T1 / T4;
    }
    else
      T9 = 1.0 / (CONSTMAX_EXP - 2.0);
    T9 = paramPtr->dvt0 * T9 * tmp1;

    T4 = toxe * paramPtr->phi
       / (paramPtr->weff + paramPtr->w0);

    T0 = sqrt(1.0 + paramPtr->lpe0 / paramPtr->leff);
    if((model_.tempMod == 1) || (model_.tempMod == 0))
      T3 = (paramPtr->kt1 + paramPtr->kt1l / paramPtr->leff)
          * (TRatio - 1.0);
    if(model_.tempMod == 2)
          T3 = - paramPtr->kt1 * (TRatio - 1.0);

    T5 = paramPtr->k1ox * (T0 - 1.0) * paramPtr->sqrtPhi + T3;

    paramPtr->vfbzbfactor
      = - T8 - T9 + paramPtr->k3 * T4 + T5 - paramPtr->phi - paramPtr->k1 * paramPtr->sqrtPhi;

    // stress effect

    wlod = model_.wlod;
    if (model_.wlod < 0.0)
    {
      UserWarning(*this) << "WLOD =is less than 0. 0.0 is used";
      wlod = 0.0;
    }
    T0 = pow(Lnew, model_.llodku0);
    W_tmp = Wnew + wlod;
    T1 = pow(W_tmp, model_.wlodku0);
    tmp1 = model_.lku0 / T0 + model_.wku0 / T1
           + model_.pku0 / (T0 * T1);
    paramPtr->ku0 = 1.0 + tmp1;

    T0 = pow(Lnew, model_.llodvth);
    T1 = pow(W_tmp, model_.wlodvth);
    tmp1 = model_.lkvth0 / T0 + model_.wkvth0 / T1
         + model_.pkvth0 / (T0 * T1);
    paramPtr->kvth0 = 1.0 + tmp1;
    paramPtr->kvth0 = sqrt(paramPtr->kvth0*paramPtr->kvth0 + DELTA);

    T0 = (TRatio - 1.0);
    paramPtr->ku0temp = paramPtr->ku0 * (1.0 + model_.tku0 *T0) + DELTA;

    Inv_saref = 1.0/(model_.saref + 0.5*Ldrn);
    Inv_sbref = 1.0/(model_.sbref + 0.5*Ldrn);
    paramPtr->inv_od_ref = Inv_saref + Inv_sbref;
    paramPtr->rho_ref = model_.ku0 / paramPtr->ku0temp * paramPtr->inv_od_ref;

  } // End of size if-statement

  //  stress effect
  if( (sa > 0.0) && (sb > 0.0) &&
      ((nf == 1.0) || ((nf > 1.0) && (sd > 0.0))) )
  {
    Inv_sa = 0;
    Inv_sb = 0;

    kvsat = model_.kvsat;
    if (model_.kvsat < -1.0 )
    {
      UserWarning(*this) << "KVSAT is too small; -1.0 is used";
      kvsat = -1.0;
    }
    if (model_.kvsat > 1.0)
    {
      UserWarning(*this) << "KVSAT is too big; 1.0 is used";
      kvsat = 1.0;
    }

    int i=0;
    for(i = 0; i < nf; i++)
    {
      T0 = 1.0 / nf / (sa + 0.5*Ldrn + i * (sd +Ldrn));
      T1 = 1.0 / nf / (sb + 0.5*Ldrn + i * (sd +Ldrn));
      Inv_sa += T0;
      Inv_sb += T1;
    }
    Inv_ODeff = Inv_sa + Inv_sb;
    rho = model_.ku0 / paramPtr->ku0temp * Inv_ODeff;
    T0 = (1.0 + rho)/(1.0 + paramPtr->rho_ref);
    u0temp = paramPtr->u0temp * T0;

    T1 = (1.0 + kvsat * rho)/(1.0 + kvsat * paramPtr->rho_ref);
    vsattemp = paramPtr->vsattemp * T1;

    OD_offset = Inv_ODeff - paramPtr->inv_od_ref;
    dvth0_lod = model_.kvth0 / paramPtr->kvth0 * OD_offset;
    dk2_lod = model_.stk2 / pow(paramPtr->kvth0, model_.lodk2) *
                           OD_offset;
    deta0_lod = model_.steta0 / pow(paramPtr->kvth0, model_.lodeta0) *
                             OD_offset;
    vth0 = paramPtr->vth0 + dvth0_lod;

    eta0 = paramPtr->eta0 + deta0_lod;
    k2 = paramPtr->k2 + dk2_lod;
  }
  else
  {
    u0temp = paramPtr->u0temp;
    vth0 = paramPtr->vth0;
    vsattemp = paramPtr->vsattemp;
    eta0 = paramPtr->eta0;
    k2 = paramPtr->k2;
  }

  //  Well Proximity Effect
  if (model_.wpemod)
  {
    if( (!scaGiven) && (!scbGiven) && (!sccGiven) )
    {
      if((scGiven) && (sc > 0.0) )
      {
        T1 = sc + Wdrn;
        T2 = 1.0 / model_.scref;
        sca = model_.scref * model_.scref / (sc * T1);
        scb = ( (0.1 * sc + 0.01 * model_.scref)
            * exp(-10.0 * sc * T2)
            - (0.1 * T1 + 0.01 * model_.scref)
            * exp(-10.0 * T1 * T2) ) / Wdrn;
        scc = ( (0.05 * sc + 0.0025 * model_.scref)
                                  * exp(-20.0 * sc * T2)
                                  - (0.05 * T1 + 0.0025 * model_.scref)
                                  * exp(-20.0 * T1 * T2) ) / Wdrn;
      }
      else
      {
        UserWarning(*this) << "No WPE as none of SCA, SCB, SCC, SC is given and/or SC not positive";
      }
    }
    sceff = sca + model_.web * scb
                + model_.wec * scc;
    vth0 += paramPtr->kvth0we * sceff;
    k2 +=  paramPtr->k2we * sceff;
    T3 =  1.0 + paramPtr->ku0we * sceff;
    if (T3 <= 0.0)
    {
      T3 = 0.0;
      UserWarning(*this) << "ku0we = %g is negatively too high. Negative mobility!";
    }
    u0temp *= T3;
  }

  // adding delvto
  vth0 += delvto;
  vfb = paramPtr->vfb + model_.dtype * delvto;

  // Instance variables calculation
  T3 = model_.dtype * vth0
         - vfb - paramPtr->phi;
  T4 = T3 + T3;
  T5 = 2.5 * T3;
  vtfbphi1 = (model_.dtype == CONSTNMOS) ? T4 : T5;
  if (vtfbphi1 < 0.0)
      vtfbphi1 = 0.0;

  vtfbphi2 = 4.0 * T3;
  if (vtfbphi2 < 0.0)
      vtfbphi2 = 0.0;

  if (k2 < 0.0)
  {
    T0 = 0.5 * paramPtr->k1 / k2;
    vbsc = 0.9 * (paramPtr->phi - T0 * T0);
    if (vbsc > -3.0)
              vbsc = -3.0;
    else if (vbsc < -30.0)
              vbsc = -30.0;
  }
  else
    vbsc = -30.0;
  if (vbsc > paramPtr->vbm)
      vbsc = paramPtr->vbm;
  k2ox = k2 * toxe
                        / model_.toxm;

  vfbzb = paramPtr->vfbzbfactor
                    +  model_.dtype * vth0 ;

  cgso = paramPtr->cgso;
  cgdo = paramPtr->cgdo;

  lnl = log(paramPtr->leff * 1.0e6);
  lnw = log(paramPtr->weff * 1.0e6);
  lnnf = log(nf);

  bodymode = 5;
  if( ( !model_.rbps0Given) || ( !model_.rbpd0Given) )
    bodymode = 1;
  else
    if( (!model_.rbsbx0Given && !model_.rbsby0Given) ||
        (!model_.rbdbx0Given && !model_.rbdby0Given) )
    bodymode = 3;

  if(rbodyMod == 2)
  {
    if (bodymode == 5)
    {
      rbsbx =  exp( log(model_.rbsbx0) + model_.rbsdbxl * lnl +
        model_.rbsdbxw * lnw + model_.rbsdbxnf * lnnf );
      rbsby =  exp( log(model_.rbsby0) + model_.rbsdbyl * lnl +
        model_.rbsdbyw * lnw + model_.rbsdbynf * lnnf );
      rbsb = rbsbx * rbsby / (rbsbx + rbsby);


      rbdbx =  exp( log(model_.rbdbx0) + model_.rbsdbxl * lnl +
        model_.rbsdbxw * lnw + model_.rbsdbxnf * lnnf );
      rbdby =  exp( log(model_.rbdby0) + model_.rbsdbyl * lnl +
        model_.rbsdbyw * lnw + model_.rbsdbynf * lnnf );
      rbdb = rbdbx * rbdby / (rbdbx + rbdby);
    }

    if ((bodymode == 3)|| (bodymode == 5))
    {
      rbps = exp( log(model_.rbps0) + model_.rbpsl * lnl +
           model_.rbpsw * lnw + model_.rbpsnf * lnnf );
      rbpd = exp( log(model_.rbpd0) + model_.rbpdl * lnl +
           model_.rbpdw * lnw + model_.rbpdnf * lnnf );
    }

    rbpbx =  exp( log(model_.rbpbx0) + model_.rbpbxl * lnl +
      model_.rbpbxw * lnw + model_.rbpbxnf * lnnf );
    rbpby =  exp( log(model_.rbpby0) + model_.rbpbyl * lnl +
      model_.rbpbyw * lnw + model_.rbpbynf * lnnf );
    rbpb = rbpbx*rbpby/(rbpbx + rbpby);
  }


  if ((rbodyMod == 1 ) || ((rbodyMod == 2 ) && (bodymode == 5)) )
  {
    if (rbdb < 1.0e-3) grbdb = 1.0e3; // in mho
    else               grbdb = model_.gbmin + 1.0 / rbdb;

    if (rbpb < 1.0e-3) grbpb = 1.0e3;
    else               grbpb = model_.gbmin + 1.0 / rbpb;

    if (rbps < 1.0e-3) grbps = 1.0e3;
    else               grbps = model_.gbmin + 1.0 / rbps;

    if (rbsb < 1.0e-3) grbsb = 1.0e3;
    else               grbsb = model_.gbmin + 1.0 / rbsb;

    if (rbpd < 1.0e-3) grbpd = 1.0e3;
    else               grbpd = model_.gbmin + 1.0 / rbpd;

  }

  if((rbodyMod == 2) && (bodymode == 3))
  {
    grbdb = grbsb = model_.gbmin;

    if (rbpb < 1.0e-3) grbpb = 1.0e3;
    else               grbpb = model_.gbmin + 1.0 / rbpb;

    if (rbps < 1.0e-3) grbps = 1.0e3;
    else               grbps = model_.gbmin + 1.0 / rbps;

    if (rbpd < 1.0e-3) grbpd = 1.0e3;
    else               grbpd = model_.gbmin + 1.0 / rbpd;
  }

  if((rbodyMod == 2) && (bodymode == 1))
  {
    grbdb = grbsb = model_.gbmin;
    grbps = grbpd = 1.0e3;
    if (rbpb < 1.0e-3)
    {
      grbpb = 1.0e3;
    }
    else
    {
      grbpb = model_.gbmin + 1.0 / rbpb;
    }
  }

  // Process geomertry dependent parasitics
  grgeltd = model_.rshg * (xgw
        + paramPtr->weffCJ / 3.0 / ngcon) /
        (ngcon * nf * (Lnew - model_.xgl));

  if (grgeltd > 0.0)
  {
     grgeltd = 1.0 / grgeltd;
  }
  else
  {
   grgeltd = 1.0e3; // mho
   if (rgateMod != 0)
   {
     UserWarning(*this) << "The gate conductance reset to 1.0e3 mho";
   }
  }

  DMCGeff = model_.dmcg - model_.dmcgt;
  DMCIeff = model_.dmci;
  DMDGeff = model_.dmdg - model_.dmcgt;

  if (sourcePerimeterGiven)
  {
   if (model_.perMod == 0)
      Pseff = sourcePerimeter;
   else
      Pseff = sourcePerimeter - paramPtr->weffCJ * nf;
  }
  else
  {
    PAeffGeo(nf, geoMod, min,
             paramPtr->weffCJ, DMCGeff, DMCIeff, DMDGeff,
            (Pseff), dumPd, dumAs, dumAd);
  }

  if (drainPerimeterGiven)
  {
   if (model_.perMod == 0)
     Pdeff = drainPerimeter;
   else
     Pdeff = drainPerimeter - paramPtr->weffCJ * nf;
  }
  else
  {
    PAeffGeo(nf, geoMod, min,
              paramPtr->weffCJ, DMCGeff, DMCIeff, DMDGeff,
              dumPs, (Pdeff), dumAs, dumAd);
  }

  if (sourceAreaGiven)
   Aseff = sourceArea;
  else
   PAeffGeo(nf, geoMod, min,
              paramPtr->weffCJ, DMCGeff, DMCIeff, DMDGeff,
              dumPs, dumPd, (Aseff), dumAd);

  if (drainAreaGiven)
    Adeff = drainArea;
  else
    PAeffGeo(nf, geoMod, min,
              paramPtr->weffCJ, DMCGeff, DMCIeff, DMDGeff,
              dumPs, dumPd, dumAs, (Adeff));

  // Processing S/D resistance and conductance below
  if (model_.rdsMod || sourceSquaresGiven)
  {
    sourceConductance = 0.0;
    if(sourceSquaresGiven)
    {
      sourceConductance = model_.sheetResistance * sourceSquares;
    }
    else if (rgeoMod > 0)
    {
      RdseffGeo(nf, geoMod,
        rgeoMod, min,
        paramPtr->weffCJ, model_.sheetResistance,
        DMCGeff, DMCIeff, DMDGeff, 1, (sourceConductance));
    }
    else
    {
      sourceConductance = 0.0;
    }

    if (sourceConductance > 0.0)
    {
      sourceConductance = 1.0 / sourceConductance;
    }
    else
    {
      sourceConductance = 1.0e3; // mho
      UserWarning(*this) << "Source conductance reset to 1.0e3 mho";
    }
  }
  else
  {
    sourceConductance = 0.0;
  }

  if (model_.rdsMod || drainSquaresGiven)
  {
    drainConductance = 0.0;
    if(drainSquaresGiven)
    {
      drainConductance = model_.sheetResistance * drainSquares;
    }
    else if (rgeoMod > 0)
    {
      RdseffGeo(nf, geoMod, rgeoMod, min,
        paramPtr->weffCJ, model_.sheetResistance,
        DMCGeff, DMCIeff, DMDGeff, 0, (drainConductance));
    }
    else
    {
      drainConductance = 0.0;
    }

    if (drainConductance > 0.0)
    {
      drainConductance = 1.0 / drainConductance;
    }
    else
    {
      drainConductance = 1.0e3; // mho
      UserWarning(*this) << "Drain conductance reset to 1.0e3 mho";
    }
  }
  else
  {
    drainConductance = 0.0;
  }

  // End of Rsd processing

  Nvtms = model_.vtm * model_.SjctEmissionCoeff;
  if ((Aseff <= 0.0) && (Pseff <= 0.0))
  {
    SourceSatCurrent = 1.0e-14;
  }
  else
  {
    SourceSatCurrent = Aseff * model_.SjctTempSatCurDensity
                       + Pseff * model_.SjctSidewallTempSatCurDensity
                       + paramPtr->weffCJ * nf
                       * model_.SjctGateSidewallTempSatCurDensity;
  }

  if (SourceSatCurrent > 0.0)
  {
    switch(model_.dioMod)
    {
      case 0:
        if ((model_.bvs / Nvtms) > CONSTEXP_THRESHOLD)
          XExpBVS = model_.xjbvs * CONSTMIN_EXP;
        else
          XExpBVS = model_.xjbvs * exp(-model_.bvs / Nvtms);
        break;
      case 1:
        DioIjthVjmEval(Nvtms, model_.ijthsfwd, SourceSatCurrent,
                    0.0, (vjsmFwd));
        IVjsmFwd = SourceSatCurrent * exp(vjsmFwd / Nvtms);
        break;
      case 2:
        if ((model_.bvs / Nvtms) > CONSTEXP_THRESHOLD)
        {
          XExpBVS = model_.xjbvs * CONSTMIN_EXP;
          tmp = CONSTMIN_EXP;
        }
        else
        {
          XExpBVS = exp(-model_.bvs / Nvtms);
          tmp = XExpBVS;
          XExpBVS *= model_.xjbvs;
        }

        DioIjthVjmEval(Nvtms, model_.ijthsfwd, SourceSatCurrent,
                                XExpBVS, (vjsmFwd));
        T0 = exp(vjsmFwd / Nvtms);
        IVjsmFwd = SourceSatCurrent * (T0 - XExpBVS / T0
                            + XExpBVS - 1.0);
        SslpFwd = SourceSatCurrent
                             * (T0 + XExpBVS / T0) / Nvtms;

        T2 = model_.ijthsrev / SourceSatCurrent;
        if (T2 < 1.0)
        {
          T2 = 10.0;
          UserWarning(*this) << "ijthsrev too small and set to 10 times IsbSat";
        }
        vjsmRev = -model_.bvs
                             - Nvtms * log((T2 - 1.0) / model_.xjbvs);
        T1 = model_.xjbvs * exp(-(model_.bvs
            + vjsmRev) / Nvtms);
        IVjsmRev = SourceSatCurrent * (1.0 + T1);
        SslpRev = -SourceSatCurrent * T1 / Nvtms;
        break;
      default:
        UserError(*this) << "Specified dioMod not matched.  dioMod = " << model_.dioMod;
    }
  }

  Nvtmd = model_.vtm * model_.DjctEmissionCoeff;
  if ((Adeff <= 0.0) && (Pdeff <= 0.0))
  {
    DrainSatCurrent = 1.0e-14;
  }
  else
  {
    DrainSatCurrent = Adeff * model_.DjctTempSatCurDensity
                    + Pdeff * model_.DjctSidewallTempSatCurDensity
                    + paramPtr->weffCJ * nf
                    * model_.DjctGateSidewallTempSatCurDensity;
  }

  if (DrainSatCurrent > 0.0)
  {
    switch(model_.dioMod)
    {
      case 0:
        if ((model_.bvd / Nvtmd) > CONSTEXP_THRESHOLD)
          XExpBVD = model_.xjbvd * CONSTMIN_EXP;
        else
          XExpBVD = model_.xjbvd * exp(-model_.bvd / Nvtmd);
        break;
      case 1:
        DioIjthVjmEval(Nvtmd, model_.ijthdfwd, DrainSatCurrent,
                            0.0, (vjdmFwd));
        IVjdmFwd = DrainSatCurrent * exp(vjdmFwd / Nvtmd);
        break;
      case 2:
        if ((model_.bvd / Nvtmd) > CONSTEXP_THRESHOLD)
        {   XExpBVD = model_.xjbvd * CONSTMIN_EXP;
            tmp = CONSTMIN_EXP;
        }
        else
        {   XExpBVD = exp(-model_.bvd / Nvtmd);
            tmp = XExpBVD;
            XExpBVD *= model_.xjbvd;
        }

        DioIjthVjmEval(Nvtmd, model_.ijthdfwd, DrainSatCurrent,
                            XExpBVD, (vjdmFwd));
        T0 = exp(vjdmFwd / Nvtmd);
        IVjdmFwd = DrainSatCurrent * (T0 - XExpBVD / T0 + XExpBVD - 1.0);
        DslpFwd = DrainSatCurrent * (T0 + XExpBVD / T0) / Nvtmd;

        T2 = model_.ijthdrev / DrainSatCurrent;
        if (T2 < 1.0)
        {
          T2 = 10.0;
          UserWarning(*this) << "ijthdrev too small and set to 10 times IdbSat";
        }
        vjdmRev = -model_.bvd
                           - Nvtmd * log((T2 - 1.0) / model_.xjbvd); // bugfix
        T1 = model_.xjbvd * exp(-(model_.bvd
           + vjdmRev) / Nvtmd);
        IVjdmRev = DrainSatCurrent * (1.0 + T1);
        DslpRev = -DrainSatCurrent * T1 / Nvtmd;
        break;
      default:
        UserError(*this) << "Specified dioMod not matched.  dioMod = " << model_.dioMod;
    }
  }

  // GEDL current reverse bias
  T0 = (TRatio - 1.0);
  model_.njtsstemp = model_.njts * (1.0 + model_.tnjts * T0);
  model_.njtsswstemp = model_.njtssw * (1.0 + model_.tnjtssw * T0);
  model_.njtsswgstemp = model_.njtsswg * (1.0 + model_.tnjtsswg * T0);
  model_.njtsdtemp = model_.njtsd * (1.0 + model_.tnjtsd * T0);
  model_.njtsswdtemp = model_.njtsswd * (1.0 + model_.tnjtsswd * T0);
  model_.njtsswgdtemp = model_.njtsswgd * (1.0 + model_.tnjtsswgd * T0);

  T7 = Eg0 / model_.vtm * T0;

  T9 = model_.xtss * T7;
  DEXP2(T9, T1);
  T9 = model_.xtsd * T7;
  DEXP2(T9, T2);
  T9 = model_.xtssws * T7;
  DEXP2(T9, T3);
  T9 = model_.xtsswd * T7;
  DEXP2(T9, T4);
  T9 = model_.xtsswgs * T7;
  DEXP2(T9, T5);
  T9 = model_.xtsswgd * T7;
  DEXP2(T9, T6);

  T10 = paramPtr->weffCJ * nf;
  SjctTempRevSatCur = T1 * Aseff * model_.jtss;
  DjctTempRevSatCur = T2 * Adeff * model_.jtsd;
  SswTempRevSatCur = T3 * Pseff * model_.jtssws;
  DswTempRevSatCur = T4 * Pdeff * model_.jtsswd;
  SswgTempRevSatCur = T5 * T10 * model_.jtsswgs;
  DswgTempRevSatCur = T6 * T10 * model_.jtsswgd;

  if(model_.mtrlMod)
  {
    /* Calculate TOXP from EOT */

    /* Calculate Vgs_eff @ Vgs = VDD with Poly Depletion Effect */
    tmp2 = vfb + paramPtr->phi;
    vddeot = model_.dtype * model_.vddeot;
    T0 = model_.epsrgate * CONSTEPS0;
    if ((paramPtr->ngate > 1.0e18) && (paramPtr->ngate < 1.0e25)
        && (vddeot > tmp2) && (T0!=0))
    {
      T1 = 1.0e6 * CONSTQ * T0 * paramPtr->ngate /
        (model_.coxe * model_.coxe);
      T8 = vddeot - tmp2;
      T4 = sqrt(1.0 + 2.0 * T8 / T1);
      T2 = 2.0 * T8 / (T4 + 1.0);
      T3 = 0.5 * T2 * T2 / T1;
      T7 = 1.12 - T3 - 0.05;
      T6 = sqrt(T7 * T7 + 0.224);
      T5 = 1.12 - 0.5 * (T7 + T6);
      Vgs_eff = vddeot - T5;
    }
    else
      Vgs_eff = vddeot;

    /* Calculate Vth @ Vds=Vbs=0 */
    V0 = paramPtr->vbi - paramPtr->phi;
    lt1 = model_.factor1* paramPtr->sqrtXdep0;
    ltw = lt1;
    T0 = paramPtr->dvt1 * paramPtr->leff / lt1;
    if (T0 < CONSTEXP_THRESHOLD)
    {
      T1 = exp(T0);
      T2 = T1 - 1.0;
      T3 = T2 * T2;
      T4 = T3 + 2.0 * T1 * CONSTMIN_EXP;
      Theta0 = T1 / T4;
    }
    else
      Theta0 = 1.0 / (CONSTMAX_EXP - 2.0);
    Delt_vth = paramPtr->dvt0 * Theta0 * V0;
    T0 = paramPtr->dvt1w * paramPtr->weff * paramPtr->leff / ltw;
    if (T0 < CONSTEXP_THRESHOLD)
    {   T1 = exp(T0);
      T2 = T1 - 1.0;
      T3 = T2 * T2;
      T4 = T3 + 2.0 * T1 * CONSTMIN_EXP;
      T5 = T1 / T4;
    }
    else
      T5 = 1.0 / (CONSTMAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
    T2 = paramPtr->dvt0w * T5 * V0;
    TempRatio =  temp / model_.tnom - 1.0;
    T0 = sqrt(1.0 + paramPtr->lpe0 / paramPtr->leff);
    T1 = paramPtr->k1ox * (T0 - 1.0) * paramPtr->sqrtPhi
      + (paramPtr->kt1 + paramPtr->kt1l / paramPtr->leff) * TempRatio;
    Vth_NarrowW = toxe * paramPtr->phi
      / (paramPtr->weff + paramPtr->w0);
    Lpe_Vb = sqrt(1.0 + paramPtr->lpeb / paramPtr->leff);
    Vth = model_.dtype * vth0 +
      (paramPtr->k1ox - paramPtr->k1)*paramPtr->sqrtPhi*Lpe_Vb
      - Delt_vth - T2 + paramPtr->k3 * Vth_NarrowW + T1;

    /* Calculate n */
    tmp1 = epssub / paramPtr->Xdep0;
    nstar = model_.vtm / Charge_q *
      (model_.coxe	+ tmp1 + paramPtr->cit);
    tmp2 = paramPtr->nfactor * tmp1;
    tmp3 = (tmp2 + paramPtr->cdsc * Theta0 + paramPtr->cit) / model_.coxe;
    if (tmp3 >= -0.5)
      n = 1.0 + tmp3;
    else
    {
      T0 = 1.0 / (3.0 + 8.0 * tmp3);
			n = (1.0 + 3.0 * tmp3) * T0;
    }

    /* Vth correction for Pocket implant */
    if (paramPtr->dvtp0 > 0.0)
    {
      T3 = paramPtr->leff + paramPtr->dvtp0 * 2.0;
      if (model_.tempMod < 2)
        T4 = model_.vtm * log(paramPtr->leff / T3);
      else
        T4 = model_.vtm0 * log(paramPtr->leff / T3);
      Vth -= n * T4;
    }
    Vgsteff = Vgs_eff-Vth;
    /* calculating Toxp */
    niter = 0;
    toxpf = toxe;
    do
    {
      toxpi = toxpf;
      tmp2 = 2.0e8 * toxpf;
      T0 = (Vgsteff + vtfbphi2) / tmp2;
      T1 = 1.0 + exp(model_.bdos * 0.7 * log(T0));
      Tcen = model_.ados * 1.9e-9 / T1;
      toxpf = toxe - epsrox/model_.epsrsub * Tcen;
      niter++;
    } while ((niter<=4)&&(fabs(toxpf-toxpi)>1e-12));
    model_.toxp = toxpf;
    model_.coxp = epsrox * CONSTEPS0 / model_.toxp;
  }

#if 0
    if (checkModel(model, here, ckt))
    {   IFuid namarray[2];
        namarray[0] = model_.name;
        namarray[1] = name;
        (*(SPfrontEnd->IFerror)) (ERR_FATAL, "Fatal error(s) detected during .5.0 parameter checking for %s in model %s", namarray);
        return(E_BADPARM);
    }
#endif

  ///////////////////////////////////////////////////////////////////////////////

  updateTemperatureCalled_ = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;

  // begin the b4ld.c parameters:
  double dgstot_dvd(0.0), dgstot_dvg(0.0), dgstot_dvs(0.0), dgstot_dvb(0.0);
  double dgdtot_dvd(0.0), dgdtot_dvg(0.0), dgdtot_dvs(0.0), dgdtot_dvb(0.0);
  double Rs(0.0), Rd(0.0);
  double dRs_dvg(0.0), dRd_dvg(0.0),  dRs_dvb(0.0), dRd_dvb(0.0);
  double dT0_dvg(0.0), dT1_dvb(0.0), dT3_dvg(0.0), dT3_dvb(0.0);

  double vgd_old(0.0), vsbd_old(0.0), vsbd(0.0);
  double SourceSatCurrent(0.0), DrainSatCurrent(0.0);
  double VgstNVt(0.0), ExpVgst(0.0);
  double czbd(0.0), czbdsw(0.0), czbdswg(0.0),
         czbs(0.0), czbssw(0.0), czbsswg(0.0);
  double evbd(0.0), evbs(0.0),  arg(0.0), sarg(0.0);
  double Vfbeff(0.0), dVfbeff_dVg(0.0), dVfbeff_dVb(0.0), V3(0.0), V4(0.0);
  double MJD(0.0), MJSWD(0.0), MJSWGD(0.0);
  double MJS(0.0), MJSWS(0.0), MJSWGS(0.0);
  double Ggidld(0.0), Ggidlg(0.0), Ggidlb(0.0);
  double Voxacc(0.0), dVoxacc_dVg(0.0), dVoxacc_dVb(0.0);
  double Voxdepinv(0.0);
  double dVoxdepinv_dVg(0.0), dVoxdepinv_dVd(0.0), dVoxdepinv_dVb(0.0);
  double VxNVt(0.0), ExpVxNVt(0.0);
  double Vaux(0.0), dVaux_dVg(0.0), dVaux_dVd(0.0), dVaux_dVb(0.0);
  double Igc(0.0), dIgc_dVg(0.0), dIgc_dVd(0.0), dIgc_dVb(0.0);
  double dIgcs_dVg(0.0), dIgcs_dVd(0.0), dIgcs_dVb(0.0);
  double  dIgcd_dVg(0.0), dIgcd_dVd(0.0), dIgcd_dVb(0.0);
  double dIgs_dVg(0.0), dIgs_dVs(0.0), dIgd_dVg(0.0), dIgd_dVd(0.0);
  double Igbacc(0.0), dIgbacc_dVg(0.0);
  double dIgbacc_dVb(0.0);
  double Igbinv(0.0), dIgbinv_dVg(0.0), dIgbinv_dVd(0.0), dIgbinv_dVb(0.0);
  double Pigcd(0.0), dPigcd_dVg(0.0), dPigcd_dVd(0.0), dPigcd_dVb(0.0);

  double Vgs_eff(0.0), Vfb(0.0);
  double Vth_NarrowW(0.0);
  double Phis(0.0), dPhis_dVb(0.0), sqrtPhis(0.0), dsqrtPhis_dVb(0.0);
  //double Vth(0.0);  // converted to instance variable
  double dVth_dVb(0.0), dVth_dVd(0.0);
  double Vgst(0.0), dVgst_dVg(0.0), dVgst_dVb(0.0), dVgs_eff_dVg(0.0);
  double Nvtms(0.0), Nvtmd(0.0);
  double Vtm(0.0), Vtm0(0.0);
  double n(0.0), dn_dVb(0.0), dn_dVd(0.0), voffcv (0.0);
  double noff(0.0), dnoff_dVd(0.0), dnoff_dVb(0.0);
  double CoxWLcen(0.0), QovCox(0.0), LINK(0.0), V0(0.0);
  double DeltaPhi(0.0), dDeltaPhi_dVg(0.0), VgDP(0.0), dVgDP_dVg(0.0);
  double Cox(0.0), Tox(0.0);
  double Tcen(0.0), dTcen_dVg(0.0), dTcen_dVd(0.0), dTcen_dVb(0.0);
  double Ccen(0.0);
  double Coxeff(0.0), dCoxeff_dVd(0.0), dCoxeff_dVg(0.0), dCoxeff_dVb(0.0);
  double Denomi(0.0), dDenomi_dVg(0.0), dDenomi_dVd(0.0), dDenomi_dVb(0.0);
  double dueff_dVg(0.0), dueff_dVd(0.0), dueff_dVb(0.0);
  double Esat(0.0);
  //double Vdsat(0.0); // made into an instance variable
  double dEsatL_dVg(0.0), dEsatL_dVd(0.0), dEsatL_dVb(0.0);
  double dVdsat_dVg(0.0), dVdsat_dVb(0.0);
  double dVdsat_dVd(0.0), Vasat(0.0), dAlphaz_dVg(0.0), dAlphaz_dVb(0.0);
  double dVasat_dVg(0.0), dVasat_dVb(0.0);
  double dVasat_dVd(0.0), Va(0.0), dVa_dVd(0.0), dVa_dVg(0.0), dVa_dVb(0.0);
  double Vbseff(0.0), dVbseff_dVb(0.0), VbseffCV(0.0), dVbseffCV_dVb(0.0);
  double Arg1(0.0), One_Third_CoxWL(0.0), Two_Third_CoxWL(0.0), Alphaz(0.0);

  double T0,dT0_dVg(0.0), dT0_dVd(0.0), dT0_dVb(0.0);
  double T1,dT1_dVg(0.0), dT1_dVd(0.0), dT1_dVb(0.0);

  double T2(0.0), dT2_dVg(0.0), dT2_dVd(0.0), dT2_dVb(0.0);
  double T3(0.0), dT3_dVg(0.0), dT3_dVd(0.0), dT3_dVb(0.0);
  double T4(0.0),               dT4_dVd(0.0), dT4_dVb(0.0);
  double T5(0.0), dT5_dVg(0.0), dT5_dVd(0.0), dT5_dVb(0.0);
  double T6(0.0), dT6_dVg(0.0), dT6_dVd(0.0), dT6_dVb(0.0);
  double T7(0.0), dT7_dVg(0.0), dT7_dVd(0.0), dT7_dVb(0.0);
  double T8(0.0), dT8_dVg(0.0), dT8_dVd(0.0), dT8_dVb(0.0);
  double T9(0.0), dT9_dVg(0.0), dT9_dVd(0.0), dT9_dVb(0.0);
  double T10(0.0), dT10_dVg(0.0), dT10_dVb(0.0), dT10_dVd(0.0);
  double T11(0.0), T12(0.0), T13(0.0), T14(0.0);
  double tmp(0.0);
  double dAbulk_dVb(0.0), Abulk0(0.0), dAbulk0_dVb(0.0);
  double Cclm(0.0), dCclm_dVg(0.0), dCclm_dVd(0.0), dCclm_dVb(0.0);
  double FP(0.0), dFP_dVg(0.0);
  double PvagTerm(0.0), dPvagTerm_dVg(0.0);
  double dPvagTerm_dVd(0.0), dPvagTerm_dVb(0.0);
  double VADITS(0.0), dVADITS_dVg(0.0), dVADITS_dVd(0.0);
  double Lpe_Vb(0.0);
  double dDITS_Sft_dVb(0.0), dDITS_Sft_dVd(0.0);
  double VACLM(0.0), dVACLM_dVg(0.0), dVACLM_dVd(0.0), dVACLM_dVb(0.0);
  double VADIBL(0.0), dVADIBL_dVg(0.0), dVADIBL_dVd(0.0), dVADIBL_dVb(0.0);
  double Xdep(0.0), dXdep_dVb(0.0);
  double lt1(0.0), dlt1_dVb(0.0), ltw(0.0), dltw_dVb(0.0);
  double Delt_vth(0.0), dDelt_vth_dVb(0.0);
  double Theta0(0.0), dTheta0_dVb(0.0);

  double TempRatio(0.0), tmp1(0.0), tmp2(0.0), tmp3(0.0), tmp4(0.0);
  double DIBL_Sft(0.0), dDIBL_Sft_dVd(0.0);
  double Lambda(0.0), dLambda_dVg(0.0);
  double a1(0.0);

  double dVgsteff_dVg(0.0), dVgsteff_dVd(0.0), dVgsteff_dVb(0.0);
  double dVdseff_dVg(0.0), dVdseff_dVd(0.0), dVdseff_dVb(0.0);
  double VdseffCV(0.0),
         dVdseffCV_dVg(0.0), dVdseffCV_dVd(0.0), dVdseffCV_dVb(0.0);
  double diffVds(0.0);
  double dAbulk_dVg(0.0);
  double beta(0.0), dbeta_dVg(0.0), dbeta_dVd(0.0), dbeta_dVb(0.0);
  double gche(0.0), dgche_dVg(0.0), dgche_dVd(0.0), dgche_dVb(0.0);
  double fgche1(0.0), dfgche1_dVg(0.0), dfgche1_dVd(0.0), dfgche1_dVb(0.0);
  double fgche2(0.0), dfgche2_dVg(0.0), dfgche2_dVd(0.0), dfgche2_dVb(0.0);
  double Idl(0.0), dIdl_dVg(0.0), dIdl_dVd(0.0), dIdl_dVb(0.0);
  double Idsa(0.0), dIdsa_dVg(0.0), dIdsa_dVd(0.0), dIdsa_dVb(0.0);
  double Ids(0.0), Gmb(0.0);
  double devbs_dvb(0.0), devbd_dvb(0.0);
  double Isub(0.0), Gbd(0.0), Gbg(0.0), Gbb(0.0), Gds(0.0);
  double VASCBE(0.0), dVASCBE_dVg(0.0), dVASCBE_dVd(0.0), dVASCBE_dVb(0.0);
  double CoxeffWovL(0.0);
  double Rds(0.0), dRds_dVg(0.0), dRds_dVb(0.0), WVCox(0.0), WVCoxRds(0.0);
  double Vgst2Vtm(0.0), VdsatCV(0.0);
  double Leff(0.0), Weff(0.0), dWeff_dVg(0.0), dWeff_dVb(0.0);
  double AbulkCV(0.0), dAbulkCV_dVb(0.0);

  double Cgg1(0.0), Cgb1(0.0), Cgd1(0.0), Cbg1(0.0), Cbb1(0.0), Cbd1(0.0);
  double Qac0(0.0), Qsub0(0.0);
  double dQac0_dVg(0.0), dQac0_dVb(0.0);
  double dQsub0_dVg(0.0), dQsub0_dVd(0.0), dQsub0_dVb(0.0);
  double Ggislg(0.0), Ggislb(0.0), Ggisls(0.0);
  double Nvtmrss(0.0), Nvtmrssws(0.0), Nvtmrsswgs(0.0);
  double Nvtmrsd(0.0), Nvtmrsswd(0.0), Nvtmrsswgd(0.0);

  double vs(0.0), Fsevl(0.0);
  double dvs_dVg(0.0), dvs_dVd(0.0), dvs_dVb(0.0), dFsevl_dVg(0.0);
  double dFsevl_dVd(0.0), dFsevl_dVb(0.0);
  double vgdx(0.0), vgsx(0.0),epssub(0.0),toxe(0.0),epsrox(0.0);
  double von_local(0.0);

  // end b4ld.c parameters

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

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag && isActive(Diag::TIME_STEP) )
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "  Instance::updateIntermediateVars\n";
    Xyce::dout() << "  name = " << getName();
    Xyce::dout() << "  model name = " << model_.getName();
    Xyce::dout() <<"   dtype is " << model_.dtype << std::endl;
    Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << "  " << std::endl;
  }

  int Check = 0;
  int Check1 = 0;
  int Check2 = 0;

  limitedFlag=false;

  // The first block of code in b4ld.c basically sets up, locally,
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
  Vb     = 0.0;
  Vsp    = 0.0;
  Vdp    = 0.0;
  Vgp    = 0.0;
  Vbp    = 0.0;
  Vge    = 0.0;
  Vgm    = 0.0;
  Vdb    = 0.0;
  Vsb    = 0.0;
  Qtotal = 0.0;

  Vd  = (extData.nextSolVectorRawPtr)[li_Drain] ;
  Vs  = (extData.nextSolVectorRawPtr)[li_Source] ;
  Vb  = (extData.nextSolVectorRawPtr)[li_Body] ;
  Vsp = (extData.nextSolVectorRawPtr)[li_SourcePrime] ;
  Vdp = (extData.nextSolVectorRawPtr)[li_DrainPrime] ;

  Vgp = (extData.nextSolVectorRawPtr)[li_GatePrime] ;
  Vbp = (extData.nextSolVectorRawPtr)[li_BodyPrime] ;
  Vge = (extData.nextSolVectorRawPtr)[li_GateExt] ;

  if (li_GateMid >= 0) // only true for rgateMod==3
  {
    Vgm = (extData.nextSolVectorRawPtr)[li_GateMid] ;
  }

  Vdb = (extData.nextSolVectorRawPtr)[li_DrainBody] ;
  Vsb = (extData.nextSolVectorRawPtr)[li_SourceBody] ;

  if (trnqsMod)
  {
    Qtotal = (extData.nextSolVectorRawPtr)[li_Charge];
  }
  else
  {
    Qtotal = 0.0;
  }

  Vddp  = Vd   - Vdp;
  Vssp  = Vs   - Vsp;
  //Vbsp  = Vb   - Vsp;
  //Vbdp  = Vb   - Vdp;
  //Vgsp  = Vg   - Vsp;
  //Vgdp  = Vg   - Vdp;
  //Vgb   = Vg   - Vb;

  //Vdpsp = Vdp  - Vsp;

  // substrate network:
  Vdbb  = Vdb - Vb;
  Vdbbp = Vdb - Vbp;
  Vsbb  = Vsb - Vb;
  Vsbbp = Vsb - Vbp;
  Vbpb  = Vbp - Vb;

  // modified from b4ld:
  vds  = model_.dtype * (Vdp - Vsp);
  vgs  = model_.dtype * (Vgp - Vsp);
  vbs  = model_.dtype * (Vbp - Vsp);
  vges = model_.dtype * (Vge - Vsp);
  vgms = model_.dtype * (Vgm - Vsp);
  vdbs = model_.dtype * (Vdb - Vsp);
  vsbs = model_.dtype * (Vsb - Vsp);
  vses = model_.dtype * (Vs - Vsp);
  vdes = model_.dtype * (Vd - Vsp);
  qdef = model_.dtype * (Qtotal);

  vbd = vbs - vds;
  vgd = vgs - vds;
  vgb = vgs - vbs;
  vged = vges - vds;
  vgmd = vgms - vds;
  vgmb = vgms - vbs;
  vdbd = vdbs - vds;

  vbs_jct = (!rbodyMod) ? vbs : vsbs;
  vbd_jct = (!rbodyMod) ? vbd : vdbd;

  // Set up the linear resistors.  We need the exact drops, not munged by
  // type.
  Vgegp = Vge - Vgp;
  Vgegm = Vge - Vgm;
  Vgmgp = Vgm - Vgp;

  origFlag = 1;

  vbd_orig = vbd;
  vbs_orig = vbs;
  vgs_orig = vgs;
  vds_orig = vds;
  vgd_orig = vgd;
  vges_orig = vges;
  vgms_orig = vgms;
  vdes_orig = vdes;
  vses_orig = vses;
  vdbs_orig = vdbs;
  vsbs_orig = vsbs;
  vdbd_orig = vdbd;
  vged_orig = vged;
  vgmd_orig = vgmd;
  vbs_jct_orig = vbs_jct;
  vbd_jct_orig = vbd_jct;
  vgmb_orig = vgmb;
  vgb_orig = vgb;

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
      if ((*flagSolVectorPtr)[li_Drain]       == 0 ||
          (*flagSolVectorPtr)[li_GateExt]     == 0 ||
          (*flagSolVectorPtr)[li_Source]      == 0 ||
          (*flagSolVectorPtr)[li_Body]        == 0 ||
          (*flagSolVectorPtr)[li_DrainPrime]  == 0 ||
          (*flagSolVectorPtr)[li_GatePrime]   == 0 ||
          (*flagSolVectorPtr)[li_GateMid]     == 0 ||
          (*flagSolVectorPtr)[li_SourcePrime] == 0 ||
          (*flagSolVectorPtr)[li_BodyPrime]   == 0 ||
          (*flagSolVectorPtr)[li_DrainBody]   == 0 ||
          (*flagSolVectorPtr)[li_SourceBody]  == 0 ||
          (*flagSolVectorPtr)[li_Charge]      == 0 )
      {
        vds = 0.1;
        vdes = 0.11;
        vses = -0.01;
        vgs = vges = vgms = model_.dtype * vth0 + 0.1;
        origFlag = 0;
      }
    }
    else
    {
      vds = 0.1;
      vdes = 0.11;
      vses = -0.01;
      vgs = vges = vgms = model_.dtype * vth0 + 0.1;
      origFlag = 0;
    }
    vbs = vdbs = vsbs = 0.0;
    vbd = vbs - vds;
    vdbd = vdbs - vds;
    vgd = vgs - vds;
    vged = vges-vds;
    vgmd = vgms-vds;
    //origFlag = 0;
  }
  else if ((getSolverState().initFixFlag || getSolverState().initJctFlag_) && OFF)
  {
    vds = vgs = vbs = vges = vgms = 0.0;
    vds = vsbs = vdes = vses = qdef = 0.0;
  }


  if (getSolverState().newtonIter == 0)
  {

    if (!getSolverState().dcopFlag || (getSolverState().locaEnabledFlag && getSolverState().dcopFlag))
    // ie, first newton step of a transient time step or DCOP continuation step.
    {
      vbd_old = (extData.currStoVectorRawPtr)[li_store_vbd];
      vbs_old = (extData.currStoVectorRawPtr)[li_store_vbs];
      vgs_old = (extData.currStoVectorRawPtr)[li_store_vgs];
      vds_old = (extData.currStoVectorRawPtr)[li_store_vds];
      vges_old = (extData.currStoVectorRawPtr)[li_store_vges];
      vgms_old = (extData.currStoVectorRawPtr)[li_store_vgms];
      vdes_old = (extData.currStoVectorRawPtr)[li_store_vdes];
      vses_old = (extData.currStoVectorRawPtr)[li_store_vses];
      vdbs_old = (extData.currStoVectorRawPtr)[li_store_vdbs];
      vsbs_old = (extData.currStoVectorRawPtr)[li_store_vsbs];
      vdbd_old = (extData.currStoVectorRawPtr)[li_store_vdbd];
      vged_old = (extData.currStoVectorRawPtr)[li_store_vged];
      vgmd_old = (extData.currStoVectorRawPtr)[li_store_vgmd];
    }
    else
    {  // no history
      vbd_old = vbd;
      vbs_old = vbs;
      vgs_old = vgs;
      vds_old = vds;
      vges_old = vges;
      vgms_old = vgms;
      vdes_old = vdes;
      vses_old = vses;
      vdbs_old = vdbs;
      vsbs_old = vsbs;
      vdbd_old = vdbd;
      vged_old = vged;
      vgmd_old = vgmd;
    }
  }
  else
  {
    vbd_old = (extData.nextStoVectorRawPtr)[li_store_vbd];
    vbs_old = (extData.nextStoVectorRawPtr)[li_store_vbs];
    vgs_old = (extData.nextStoVectorRawPtr)[li_store_vgs];
    vds_old = (extData.nextStoVectorRawPtr)[li_store_vds];
    vges_old = (extData.nextStoVectorRawPtr)[li_store_vges];
    vgms_old = (extData.nextStoVectorRawPtr)[li_store_vgms];
    vdes_old = (extData.nextStoVectorRawPtr)[li_store_vdes];
    vses_old = (extData.nextStoVectorRawPtr)[li_store_vses];
    vdbs_old = (extData.nextStoVectorRawPtr)[li_store_vdbs];
    vsbs_old = (extData.nextStoVectorRawPtr)[li_store_vsbs];
    vdbd_old = (extData.nextStoVectorRawPtr)[li_store_vdbd];
    vged_old = (extData.nextStoVectorRawPtr)[li_store_vged];
    vgmd_old = (extData.nextStoVectorRawPtr)[li_store_vgmd];
  }

  vgd_old = vgs_old - vds_old;

  // This next block performs checks on the junction voltages and
  // imposes limits on them if they are too big.
  // Note:  In the level=1 von is multiplied by dtype.  Here it is not.  They
  // are both right.

  if (getDeviceOptions().voltageLimiterFlag && !(getSolverState().initFixFlag && OFF))
  {

#if 0
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << "  CONSTvt0  = " << CONSTvt0 << std::endl;
      Xyce::dout() << "  vcrit     = " << model_.vcrit << std::endl;
#if 0
      Xyce::dout().width(3);
      Xyce::dout() << getSolverState().newtonIter;
      Xyce::dout().width(5);Xyce::dout() << getName();
      Xyce::dout() << " old :";
      Xyce::dout()<<" vgs:";//Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vgs_old;
      Xyce::dout()<<" vds:";//Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vds_old;
      Xyce::dout()<<" vbs:";//Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vbs_old;
      Xyce::dout()<<" vbd:";//Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vbd_old;
      Xyce::dout()<<" vges:";//Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vges_old;
      Xyce::dout()<<" vgms:";//Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vgms_old;
      Xyce::dout()<<" vged:";//Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vged_old;
      Xyce::dout()<<" vgmd:";//Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vgmd_old << std::endl;
#endif
      Xyce::dout().width(3);
      Xyce::dout() << getSolverState().newtonIter;
      Xyce::dout().width(5);Xyce::dout() << getName();
      Xyce::dout() << " Blim:";
      Xyce::dout()<<" vgs:";//Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vgs;
      Xyce::dout()<<" vds:";//Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vds;
      Xyce::dout()<<" vbs:";//Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vbs;
      Xyce::dout()<<" vbd:";//Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vbd;
      Xyce::dout()<<" vges:";//Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vges;
      Xyce::dout()<<" vgms:";//Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vgms;
      Xyce::dout()<<" vged:";//Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vged;
      Xyce::dout()<<" vgmd:";//Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vgmd << std::endl;
      Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
    }
#endif

    // only do this if we are beyond the first Newton iteration.  On the
    // first newton iteration, the "old" values are from a previous time
    // step.
    von_local=von;
    if (getSolverState().newtonIter >= 0 && !(getSolverState().initJctFlag_))
    {
      if (vds_old >= 0.0)
      {
        vgs = devSupport.fetlim(vgs, vgs_old, von_local);
        vds = vgs - vgd;
        vds = devSupport.limvds(vds, vds_old);
        vgd = vgs - vds;
        if (rgateMod == 3)
        {
         vges = devSupport.fetlim(vges, vges_old, von_local);
         vgms = devSupport.fetlim(vgms, vgms_old, von_local);
         vged = vges - vds;
         vgmd = vgms - vds;
        }
        else if ((rgateMod == 1) || (rgateMod == 2))
        {
          vges = devSupport.fetlim(vges, vges_old, von_local);
          vged = vges - vds;
        }

        if (model_.rdsMod)
        {
          vdes = devSupport.limvds(vdes, vdes_old);
          vses = -devSupport.limvds(-vses, -vses_old);
        }
      }
      else
      {
        vgd = devSupport.fetlim(vgd, vgd_old, von_local);
        vds = vgs - vgd;
        vds = -devSupport.limvds(-vds, -vds_old);
        vgs = vgd + vds;

        if (rgateMod == 3)
        {
          vged = devSupport.fetlim(vged, vged_old, von_local);
          vges = vged + vds;
          vgmd = devSupport.fetlim(vgmd, vgmd_old, von_local);
          vgms = vgmd + vds;
        }
        if ((rgateMod == 1) || (rgateMod == 2))
        {
          vged = devSupport.fetlim(vged, vged_old, von_local);
          vges = vged + vds;
        }

        if (model_.rdsMod)
        {
          vdes = -devSupport.limvds(-vdes, -vdes_old);
          vses = devSupport.limvds(vses, vses_old);
        }
      }

      if (vds >= 0.0)
      {
        vbs = devSupport.pnjlim(vbs, vbs_old,
                 CONSTvt0, model_.vcrit, &Check);
        vbd = vbs - vds;
        if (rbodyMod)
        {
          vdbs = devSupport.pnjlim(vdbs, vdbs_old,
                      CONSTvt0, model_.vcrit, &Check1);
          vdbd = vdbs - vds;
          vsbs = devSupport.pnjlim(vsbs, vsbs_old,
                      CONSTvt0, model_.vcrit, &Check2);
          if ((Check1 != 0) || (Check2 != 0))
          {
            Check = 1;
          }
        }
      }
      else
      {
        vbd = devSupport.pnjlim(vbd, vbd_old,
                   CONSTvt0, model_.vcrit, &Check);
        vbs = vbd + vds;
        if (rbodyMod)
        {
          vdbd = devSupport.pnjlim(vdbd, vdbd_old,
                        CONSTvt0, model_.vcrit, &Check1);
          vdbs = vdbd + vds;
          vsbd_old = vsbs_old - vds_old;
          vsbd = vsbs - vds;
          vsbd = devSupport.pnjlim(vsbd, vsbd_old, CONSTvt0, model_.vcrit, &Check2);
          vsbs = vsbd + vds;
          if ((Check1 != 0) || (Check2 != 0))
          {
            Check = 1;
          }
        }
      }
    }

    // for convergence testing:
    if (Check == 1) limitedFlag=true;

  } // getDeviceOptions().voltageLimiterFlag

  // Calculate DC currents and their derivatives
  vbd = vbs - vds;
  vgd = vgs - vds;
  vgb = vgs - vbs;
  vged = vges - vds;
  vgmd = vgms - vds;
  vgmb = vgms - vbs;
  vdbd = vdbs - vds;

  vbs_jct = (!rbodyMod) ? vbs : vsbs;
  vbd_jct = (!rbodyMod) ? vbd : vdbd;

  if (getDeviceOptions().voltageLimiterFlag)
  {
    double machprec= Util::MachineDependentParams::MachinePrecision();

    if (
       fabs( vbs_orig - vbs) > machprec ||
       fabs( vgs_orig - vgs) > machprec ||
       fabs( vds_orig - vds) > machprec ||
       fabs( vges_orig - vges) > machprec ||
       fabs( vgms_orig - vgms) > machprec ||
       fabs( vdes_orig - vdes) > machprec ||
       fabs( vses_orig - vses) > machprec ||
       fabs( vdbs_orig - vdbs) > machprec ||
       fabs( vsbs_orig - vsbs) > machprec 
       )
    {
      origFlag = 0;
    }
  }


  if (DEBUG_DEVICE && getDeviceOptions().voltageLimiterFlag)
  {
    if (isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      if (!origFlag)
      {
        double machprec = Util::MachineDependentParams::MachinePrecision();
#if 0
        double vbd_update = vbd_orig-vbd;
        double vbs_update = vbs_orig-vbs;
        double vgs_update = vgs_orig-vgs;
        double vds_update = vds_orig-vds;
        double vgd_update = vgd_orig-vgd;
        double vges_update = vges_orig-vges;
        double vgms_update = vgms_orig-vgms;
        double vdes_update = vdes_orig-vdes;
        double vses_update = vses_orig-vses;
        double vdbs_update = vdbs_orig-vdbs;
        double vsbs_update = vsbs_orig-vsbs;
        double vdbd_update = vdbd_orig-vdbd;
        double vbs_jct_update = vbs_jct_orig-vbs_jct;
        double vbd_jct_update = vbd_jct_orig-vbd_jct;
        double vgmb_update = vgmb_orig-vgmb;
        double vgb_update = vgb_orig-vgb;
        double vged_update = vged_orig-vged;
        double vgmd_update = vgmd_orig-vgmd;
#endif
#define SET_UPDATE(VAR) double VAR ## _update = VAR ## _orig -  VAR;
#define ORIG_OUTPUT(VAR) if ( fabs(VAR ## _update) > machprec ) Xyce::dout()<<" "#VAR ":"	  << ((VAR ## _orig>=0)?"+":"")  << VAR ## _orig;
#define LIMIT_OUTPUT(VAR) if ( fabs(VAR ## _update) > machprec ) Xyce::dout()<<" "#VAR ":"	  << ((VAR >=0)?"+":"")  << VAR;
#define DIFF_OUTPUT(VAR) if ( fabs(VAR ## _update) > machprec ) Xyce::dout()<<" "#VAR ":"	  << ((VAR ## _update >=0)?"+":"")  << VAR ## _update;

        SET_UPDATE(vbd) SET_UPDATE(vbs) SET_UPDATE(vgs) SET_UPDATE(vds) SET_UPDATE(vgd) SET_UPDATE(vges)
        SET_UPDATE(vgms) SET_UPDATE(vdes) SET_UPDATE(vses) SET_UPDATE(vdbs) SET_UPDATE(vsbs) SET_UPDATE(vdbd)
        SET_UPDATE(vbs_jct) SET_UPDATE(vbd_jct) SET_UPDATE(vgmb) SET_UPDATE(vged) SET_UPDATE(vgmd)
        SET_UPDATE(vgb)

        Xyce::dout().width(3);
        Xyce::dout() << getSolverState().newtonIter;
        Xyce::dout().width(5);Xyce::dout() << getName();
        Xyce::dout() << " Blim:";
        Xyce::dout().width(12); Xyce::dout().precision(4); Xyce::dout().setf(std::ios::scientific);
        LIMIT_OUTPUT(vbd) LIMIT_OUTPUT(vbs) LIMIT_OUTPUT(vgs) LIMIT_OUTPUT(vds) LIMIT_OUTPUT(vgd) LIMIT_OUTPUT(vges)
        LIMIT_OUTPUT(vgms) LIMIT_OUTPUT(vdes) LIMIT_OUTPUT(vses) LIMIT_OUTPUT(vdbs) LIMIT_OUTPUT(vsbs) LIMIT_OUTPUT(vdbd)
        LIMIT_OUTPUT(vbs_jct) LIMIT_OUTPUT(vbd_jct) LIMIT_OUTPUT(vgmb) LIMIT_OUTPUT(vged) LIMIT_OUTPUT(vgmd)
        LIMIT_OUTPUT(vgb)
        Xyce::dout() << std::endl;

        Xyce::dout().width(3);
        Xyce::dout() << getSolverState().newtonIter;
        Xyce::dout().width(5);Xyce::dout() << getName();
        Xyce::dout() << " Alim:";
        Xyce::dout().width(12); Xyce::dout().precision(4); Xyce::dout().setf(std::ios::scientific);
        LIMIT_OUTPUT(vbd) LIMIT_OUTPUT(vbs) LIMIT_OUTPUT(vgs) LIMIT_OUTPUT(vds) LIMIT_OUTPUT(vgd) LIMIT_OUTPUT(vges)
        LIMIT_OUTPUT(vgms) LIMIT_OUTPUT(vdes) LIMIT_OUTPUT(vses) LIMIT_OUTPUT(vdbs) LIMIT_OUTPUT(vsbs) LIMIT_OUTPUT(vdbd)
        LIMIT_OUTPUT(vbs_jct) LIMIT_OUTPUT(vbd_jct) LIMIT_OUTPUT(vgmb) LIMIT_OUTPUT(vged) LIMIT_OUTPUT(vgmd)
        LIMIT_OUTPUT(vgb)
        Xyce::dout()  <<std::endl;

        Xyce::dout().width(3);
        Xyce::dout() << getSolverState().newtonIter;
        Xyce::dout().width(5);Xyce::dout() << getName();
        Xyce::dout() << " Dlim:";
        Xyce::dout().width(12); Xyce::dout().precision(4); Xyce::dout().setf(std::ios::scientific);
        DIFF_OUTPUT(vbd) DIFF_OUTPUT(vbs) DIFF_OUTPUT(vgs) DIFF_OUTPUT(vds) DIFF_OUTPUT(vgd) DIFF_OUTPUT(vges)
        DIFF_OUTPUT(vgms) DIFF_OUTPUT(vdes) DIFF_OUTPUT(vses) DIFF_OUTPUT(vdbs) DIFF_OUTPUT(vsbs) DIFF_OUTPUT(vdbd)
        DIFF_OUTPUT(vbs_jct) DIFF_OUTPUT(vbd_jct) DIFF_OUTPUT(vgmb) DIFF_OUTPUT(vged) DIFF_OUTPUT(vgmd)
        DIFF_OUTPUT(vgb)
        Xyce::dout()  <<std::endl;
      }
    }
  }

  // Source/drain junction diode DC model begins
  Nvtms = model_.vtm * model_.SjctEmissionCoeff;
  if ((Aseff <= 0.0) && (Pseff <= 0.0))
  {
    SourceSatCurrent = 1.0e-14;
  }
  else
  {
    SourceSatCurrent = Aseff * model_.SjctTempSatCurDensity
                     + Pseff * model_.SjctSidewallTempSatCurDensity
                     + paramPtr->weffCJ * nf
                     * model_.SjctGateSidewallTempSatCurDensity;
  }

  if (SourceSatCurrent <= 0.0)
  {
    gbs = getDeviceOptions().gmin;
    cbs = gbs * vbs_jct;
  }
  else
  {
    switch(model_.dioMod)
    {
      case 0:
          evbs = exp(vbs_jct / Nvtms);
          T1 = model_.xjbvs * exp(-(model_.bvs + vbs_jct) / Nvtms);
          // WDLiu: Magic T1 in this form; different from  beta.
          gbs = SourceSatCurrent * (evbs + T1) / Nvtms + getDeviceOptions().gmin;
          cbs = SourceSatCurrent * (evbs + XExpBVS
                                    - T1 - 1.0) + getDeviceOptions().gmin * vbs_jct;
          break;
      case 1:
        T2 = vbs_jct / Nvtms;
        if (T2 < -CONSTEXP_THRESHOLD)
        {
          gbs = getDeviceOptions().gmin;
          cbs = SourceSatCurrent * (CONSTMIN_EXP - 1.0)
                + getDeviceOptions().gmin * vbs_jct;
        }
        else if (vbs_jct <= vjsmFwd)
        {
          evbs = exp(T2);
          gbs = SourceSatCurrent * evbs / Nvtms + getDeviceOptions().gmin;
          cbs = SourceSatCurrent * (evbs - 1.0)
                + getDeviceOptions().gmin * vbs_jct;
        }
        else
        {
          T0 = IVjsmFwd / Nvtms;
          gbs = T0 + getDeviceOptions().gmin;
          cbs = IVjsmFwd - SourceSatCurrent + T0
                * (vbs_jct - vjsmFwd) + getDeviceOptions().gmin * vbs_jct;
        }
        break;
      case 2:
        if (vbs_jct < vjsmRev)
        {
          T0 = vbs_jct / Nvtms;
          if (T0 < -CONSTEXP_THRESHOLD)
          {
            evbs = CONSTMIN_EXP;
            devbs_dvb = 0.0;
          }
          else
          {
            evbs = exp(T0);
            devbs_dvb = evbs / Nvtms;
          }

          T1 = evbs - 1.0;
          T2 = IVjsmRev + SslpRev * (vbs_jct - vjsmRev);
          gbs = devbs_dvb * T2 + T1 * SslpRev + getDeviceOptions().gmin;
          cbs = T1 * T2 + getDeviceOptions().gmin * vbs_jct;
        }
        else if (vbs_jct <= vjsmFwd)
        {
          T0 = vbs_jct / Nvtms;
          if (T0 < -CONSTEXP_THRESHOLD)
          {
            evbs = CONSTMIN_EXP;
            devbs_dvb = 0.0;
          }
          else
          {
            evbs = exp(T0);
            devbs_dvb = evbs / Nvtms;
          }

          T1 = (model_.bvs + vbs_jct) / Nvtms;
          if (T1 > CONSTEXP_THRESHOLD)
          {
            T2 = CONSTMIN_EXP;
            T3 = 0.0;
          }
          else
          {
            T2 = exp(-T1);
            T3 = -T2 /Nvtms;
          }
          gbs = SourceSatCurrent * (devbs_dvb - model_.xjbvs * T3)
                + getDeviceOptions().gmin;
          cbs = SourceSatCurrent * (evbs + XExpBVS - 1.0
                                    - model_.xjbvs * T2)
                + getDeviceOptions().gmin * vbs_jct;
        }
        else
        {
          gbs = SslpFwd + getDeviceOptions().gmin;
          cbs = IVjsmFwd + SslpFwd * (vbs_jct
                                      - vjsmFwd) + getDeviceOptions().gmin * vbs_jct;
        }
        break;
    default: break;
    }
  }

  Nvtmd = model_.vtm * model_.DjctEmissionCoeff;
  if ((Adeff <= 0.0) && (Pdeff <= 0.0))
  {
    DrainSatCurrent = 1.0e-14;
  }
  else
  {
    DrainSatCurrent = Adeff * model_.DjctTempSatCurDensity
                    + Pdeff * model_.DjctSidewallTempSatCurDensity
                    + paramPtr->weffCJ * nf
                    * model_.DjctGateSidewallTempSatCurDensity;
  }

  if (DrainSatCurrent <= 0.0)
  {
    gbd = getDeviceOptions().gmin;
    cbd = gbd * vbd_jct;
  }
  else
  {
    switch(model_.dioMod)
    {
      case 0:
          evbd = exp(vbd_jct / Nvtmd);
          T1 = model_.xjbvd * exp(-(model_.bvd + vbd_jct) / Nvtmd);
          // WDLiu: Magic T1 in this form; different from  beta.
          gbd = DrainSatCurrent * (evbd + T1) / Nvtmd + getDeviceOptions().gmin;
          cbd = DrainSatCurrent * (evbd + XExpBVD
                                   - T1 - 1.0) + getDeviceOptions().gmin * vbd_jct;
          break;
      case 1:
          T2 = vbd_jct / Nvtmd;
          if (T2 < -CONSTEXP_THRESHOLD)
          {
            gbd = getDeviceOptions().gmin;
            cbd = DrainSatCurrent * (CONSTMIN_EXP - 1.0)
                  + getDeviceOptions().gmin * vbd_jct;
          }
          else if (vbd_jct <= vjdmFwd)
          {
            evbd = exp(T2);
            gbd = DrainSatCurrent * evbd / Nvtmd + getDeviceOptions().gmin;
            cbd = DrainSatCurrent * (evbd - 1.0)
                  + getDeviceOptions().gmin * vbd_jct;
          }
          else
          {
            T0 = IVjdmFwd / Nvtmd;
            gbd = T0 + getDeviceOptions().gmin;
            cbd = IVjdmFwd - DrainSatCurrent + T0
                  * (vbd_jct - vjdmFwd) + getDeviceOptions().gmin * vbd_jct;
          }
          break;
      case 2:
          if (vbd_jct < vjdmRev)
          {
            T0 = vbd_jct / Nvtmd;
            if (T0 < -CONSTEXP_THRESHOLD)
            {
              evbd = CONSTMIN_EXP;
              devbd_dvb = 0.0;
            }
            else
            {
              evbd = exp(T0);
              devbd_dvb = evbd / Nvtmd;
            }

            T1 = evbd - 1.0;
            T2 = IVjdmRev + DslpRev * (vbd_jct - vjdmRev);
            gbd = devbd_dvb * T2 + T1 * DslpRev + getDeviceOptions().gmin;
            cbd = T1 * T2 + getDeviceOptions().gmin * vbd_jct;
          }
          else if (vbd_jct <= vjdmFwd)
          {
            T0 = vbd_jct / Nvtmd;
            if (T0 < -CONSTEXP_THRESHOLD)
            {
              evbd = CONSTMIN_EXP;
              devbd_dvb = 0.0;
            }
            else
            {
              evbd = exp(T0);
              devbd_dvb = evbd / Nvtmd;
            }

            T1 = (model_.bvd + vbd_jct) / Nvtmd;
            if (T1 > CONSTEXP_THRESHOLD)
            {
              T2 = CONSTMIN_EXP;
              T3 = 0.0;
            }
            else
            {
              T2 = exp(-T1);
              T3 = -T2 /Nvtmd;
            }
            gbd = DrainSatCurrent * (devbd_dvb - model_.xjbvd * T3)
                  + getDeviceOptions().gmin;
            cbd = DrainSatCurrent * (evbd + XExpBVD - 1.0
                                     - model_.xjbvd * T2) + getDeviceOptions().gmin * vbd_jct;
          }
          else
          {
            gbd = DslpFwd + getDeviceOptions().gmin;
            cbd = IVjdmFwd + DslpFwd * (vbd_jct
                                        - vjdmFwd) + getDeviceOptions().gmin * vbd_jct;
          }
          break;
      default: break;
    }
  }

  /* trap-assisted tunneling and recombination current for reverse bias  */
  Nvtmrssws = model_.vtm0 * model_.njtsswstemp;
  Nvtmrsswgs = model_.vtm0 * model_.njtsswgstemp;
  Nvtmrss = model_.vtm0 * model_.njtsstemp;
  Nvtmrsswd = model_.vtm0 * model_.njtsswdtemp;
  Nvtmrsswgd = model_.vtm0 * model_.njtsswgdtemp;
  Nvtmrsd = model_.vtm0 * model_.njtsdtemp;

  if ((model_.vtss - vbs_jct) < (model_.vtss * 1e-3))
  {
    T9 = 1.0e3;
    T0 = - vbs_jct / Nvtmrss * T9;
    DEXP(T0, T1, T10);
    dT1_dVb = T10 / Nvtmrss * T9;
  }
  else
  {
    T9 = 1.0 / (model_.vtss - vbs_jct);
    T0 = -vbs_jct / Nvtmrss * model_.vtss * T9;
    dT0_dVb = model_.vtss / Nvtmrss * (T9 + vbs_jct * T9 * T9) ;
    DEXP(T0, T1, T10);
    dT1_dVb = T10 * dT0_dVb;
  }

  if ((model_.vtsd - vbd_jct) < (model_.vtsd * 1e-3) )
  {
    T9 = 1.0e3;
    T0 = -vbd_jct / Nvtmrsd * T9;
    DEXP(T0, T2, T10);
    dT2_dVb = T10 / Nvtmrsd * T9;
  }
  else
  {
    T9 = 1.0 / (model_.vtsd - vbd_jct);
    T0 = -vbd_jct / Nvtmrsd * model_.vtsd * T9;
    dT0_dVb = model_.vtsd / Nvtmrsd * (T9 + vbd_jct * T9 * T9) ;
    DEXP(T0, T2, T10);
    dT2_dVb = T10 * dT0_dVb;
  }

  if ((model_.vtssws - vbs_jct) < (model_.vtssws * 1e-3) )
  {
    T9 = 1.0e3;
    T0 = -vbs_jct / Nvtmrssws * T9;
    DEXP(T0, T3, T10);
    dT3_dVb = T10 / Nvtmrssws * T9;
  }
  else
  {
    T9 = 1.0 / (model_.vtssws - vbs_jct);
    T0 = -vbs_jct / Nvtmrssws * model_.vtssws * T9;
    dT0_dVb = model_.vtssws / Nvtmrssws * (T9 + vbs_jct * T9 * T9) ;
    DEXP(T0, T3, T10);
    dT3_dVb = T10 * dT0_dVb;
  }

  if ((model_.vtsswd - vbd_jct) < (model_.vtsswd * 1e-3) )
  {
    T9 = 1.0e3;
    T0 = -vbd_jct / Nvtmrsswd * T9;
    DEXP(T0, T4, T10);
    dT4_dVb = T10 / Nvtmrsswd * T9;
  }
  else
  {
    T9 = 1.0 / (model_.vtsswd - vbd_jct);
    T0 = -vbd_jct / Nvtmrsswd * model_.vtsswd * T9;
    dT0_dVb = model_.vtsswd / Nvtmrsswd * (T9 + vbd_jct * T9 * T9) ;
    DEXP(T0, T4, T10);
    dT4_dVb = T10 * dT0_dVb;
  }

  if ((model_.vtsswgs - vbs_jct) < (model_.vtsswgs * 1e-3) )
  {
    T9 = 1.0e3;
    T0 = -vbs_jct / Nvtmrsswgs * T9;
    DEXP(T0, T5, T10);
    dT5_dVb = T10 / Nvtmrsswgs * T9;
  }
  else
  {
    T9 = 1.0 / (model_.vtsswgs - vbs_jct);
    T0 = -vbs_jct / Nvtmrsswgs * model_.vtsswgs * T9;
    dT0_dVb = model_.vtsswgs / Nvtmrsswgs * (T9 + vbs_jct * T9 * T9) ;
    DEXP(T0, T5, T10);
    dT5_dVb = T10 * dT0_dVb;
  }

  if ((model_.vtsswgd - vbd_jct) < (model_.vtsswgd * 1e-3) )
  {
    T9 = 1.0e3;
    T0 = -vbd_jct / Nvtmrsswgd * T9;
    DEXP(T0, T6, T10);
    dT6_dVb = T10 / Nvtmrsswgd * T9;
  }
  else
  {
    T9 = 1.0 / (model_.vtsswgd - vbd_jct);
    T0 = -vbd_jct / Nvtmrsswgd * model_.vtsswgd * T9;
    dT0_dVb = model_.vtsswgd / Nvtmrsswgd * (T9 + vbd_jct * T9 * T9) ;
    DEXP(T0, T6, T10);
    dT6_dVb = T10 * dT0_dVb;
  }

  gbs += SjctTempRevSatCur * dT1_dVb
        + SswTempRevSatCur * dT3_dVb
        + SswgTempRevSatCur * dT5_dVb;
  cbs -= SjctTempRevSatCur * (T1 - 1.0)
        + SswTempRevSatCur * (T3 - 1.0)
        + SswgTempRevSatCur * (T5 - 1.0);
  gbd += DjctTempRevSatCur * dT2_dVb
        + DswTempRevSatCur * dT4_dVb
        + DswgTempRevSatCur * dT6_dVb;
  cbd -= DjctTempRevSatCur * (T2 - 1.0)
        + DswTempRevSatCur * (T4 - 1.0)
        + DswgTempRevSatCur * (T6 - 1.0);

  // End of diode DC model

  if (vds >= 0.0)
  {
    mode = 1;
    Vds = vds;
    Vgs = vgs;
    Vbs = vbs;
    Vdb = vds - vbs;  // WDLiu: for GIDL
  }
  else
  {
    mode = -1;
    Vds = -vds;
    Vgs = vgd;
    Vbs = vbd;
    Vdb = -vbs;
  }

  // dunga
  if(model_.mtrlMod)
  {
    epsrox = 3.9;
    toxe = model_.eot;
    epssub = CONSTEPS0 * model_.epsrsub;
  }
  else
  {
    epsrox = model_.epsrox;
    toxe = model_.toxe;
    epssub = CONSTEPSSI;
  }

  /////////////////////////////////////////////////////////////////////////////
  // mosfet continuation.
  // This idea is based, loosely, on a paper by Jaijeet
  // Rosychowdhury.  If the artificial parameter flag has been enabled,
  // modify Vds and Vgs.
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    if ( getSolverState().dcopFlag || getSolverState().tranopFlag )
    {
      Xyce::dout() << "HOMOTOPY INFO: gainscale   = " << getSolverState().gainScale_[blockHomotopyID] << std::endl
                   << "HOMOTOPY INFO: before vds  = " << Vds << std::endl
                   << "HOMOTOPY INFO: before vgst = " << Vgs << std::endl;
    }
  }
  if (getSolverState().artParameterFlag_)
  {

    double alpha = getSolverState().gainScale_[blockHomotopyID];
    if (getDeviceOptions().staggerGainScale)
    {
      alpha *= (0.3 * randomPerturb + 1.0);
      if (alpha > 1.0)
      {
        alpha = 1.0;
      }
    }
    double vgstConst = getDeviceOptions().vgstConst;
    if (getDeviceOptions().randomizeVgstConst)
    {
      vgstConst *= randomPerturb;
    }

    Vds = devSupport.contVds (Vds,getSolverState().nltermScale_, getDeviceOptions().vdsScaleMin);
    Vgs = devSupport.contVgst(Vgs, alpha, vgstConst);
  }
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag && isActive(Diag::TIME_STEP))
  {
    if ( getSolverState().dcopFlag || getSolverState().tranopFlag )
    {
      Xyce::dout() << "HOMOTOPY INFO: after vds   = " << Vds << std::endl;
      Xyce::dout() << "HOMOTOPY INFO: after vgst  = " << Vgs << std::endl;
    }
  }
  // end of mosfet continuation block.
  /////////////////////////////////////////////////////////////////////////////

  T0 = Vbs - vbsc - 0.001;
  T1 = sqrt(T0 * T0 - 0.004 * vbsc);
  if (T0 >= 0.0)
  {
    Vbseff = vbsc + 0.5 * (T0 + T1);
    dVbseff_dVb = 0.5 * (1.0 + T0 / T1);
  }
  else
  {
    T2 = -0.002 / (T1 - T0);
    Vbseff = vbsc * (1.0 + T2);
    dVbseff_dVb = T2 * vbsc / T1;
  }

  // JX: Correction to forward body bias
  T9 = 0.95 * paramPtr->phi;
  T0 = T9 - Vbseff - 0.001;
  T1 = sqrt(T0 * T0 + 0.004 * T9);
  Vbseff = T9 - 0.5 * (T0 + T1);
  dVbseff_dVb *= 0.5 * (1.0 + T0 / T1);

  Phis = paramPtr->phi - Vbseff;
  dPhis_dVb = -1.0;
  sqrtPhis = sqrt(Phis);
  dsqrtPhis_dVb = -0.5 / sqrtPhis;

  Xdep = paramPtr->Xdep0 * sqrtPhis / paramPtr->sqrtPhi;
  dXdep_dVb = (paramPtr->Xdep0 / paramPtr->sqrtPhi) * dsqrtPhis_dVb;

  Leff = paramPtr->leff;
  Vtm = model_.vtm;
  Vtm0 = model_.vtm0;

  // Vth Calculation
  T3 = sqrt(Xdep);
  V0 = paramPtr->vbi - paramPtr->phi;

  T0 = paramPtr->dvt2 * Vbseff;
  if (T0 >= - 0.5)
  {
    T1 = 1.0 + T0;
    T2 = paramPtr->dvt2;
  }
  else
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
  else
  {
    T4 = 1.0 / (3.0 + 8.0 * T0);
    T1 = (1.0 + 3.0 * T0) * T4;
    T2 = paramPtr->dvt2w * T4 * T4;
  }
  ltw = model_.factor1 * T3 * T1;
  dltw_dVb = model_.factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

  T0 = paramPtr->dvt1 * Leff / lt1;
  if (T0 < CONSTEXP_THRESHOLD)
  {
    T1 = exp(T0);
    T2 = T1 - 1.0;
    T3 = T2 * T2;
    T4 = T3 + 2.0 * T1 * CONSTMIN_EXP;
    Theta0 = T1 / T4;
    dT1_dVb = -T0 * T1 * dlt1_dVb / lt1;
    dTheta0_dVb = dT1_dVb * (T4 - 2.0 * T1 * (T2 + CONSTMIN_EXP)) / T4 / T4;
  }
  else
  {
    Theta0 = 1.0 / (CONSTMAX_EXP - 2.0); // 3.0 * CONSTMIN_EXP omitted
    dTheta0_dVb = 0.0;
  }
  thetavth = paramPtr->dvt0 * Theta0;
  Delt_vth = thetavth * V0;
  dDelt_vth_dVb = paramPtr->dvt0 * dTheta0_dVb * V0;

  T0 = paramPtr->dvt1w * paramPtr->weff * Leff / ltw;
  if (T0 < CONSTEXP_THRESHOLD)
  {
    T1 = exp(T0);
    T2 = T1 - 1.0;
    T3 = T2 * T2;
    T4 = T3 + 2.0 * T1 * CONSTMIN_EXP;
    T5 = T1 / T4;
    dT1_dVb = -T0 * T1 * dltw_dVb / ltw;
    dT5_dVb = dT1_dVb * (T4 - 2.0 * T1 * (T2 + CONSTMIN_EXP)) / T4 / T4;
  }
  else
  {
    T5 = 1.0 / (CONSTMAX_EXP - 2.0); // 3.0 * CONSTMIN_EXP omitted
    dT5_dVb = 0.0;
  }

  T0 = paramPtr->dvt0w * T5;
  T2 = T0 * V0;
  dT2_dVb = paramPtr->dvt0w * dT5_dVb * V0;

  TempRatio =  temp / model_.tnom - 1.0;
  T0 = sqrt(1.0 + paramPtr->lpe0 / Leff);
  T1 = paramPtr->k1ox * (T0 - 1.0) * paramPtr->sqrtPhi
     + (paramPtr->kt1 + paramPtr->kt1l / Leff
     + paramPtr->kt2 * Vbseff) * TempRatio;
  Vth_NarrowW = toxe * paramPtr->phi
              / (paramPtr->weff + paramPtr->w0);

  T3 = eta0 + paramPtr->etab * Vbseff;
  if (T3 < 1.0e-4)
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

  Lpe_Vb = sqrt(1.0 + paramPtr->lpeb / Leff);

  Vth = model_.dtype * vth0 + (paramPtr->k1ox * sqrtPhis
      - paramPtr->k1 * paramPtr->sqrtPhi) * Lpe_Vb
      - k2ox * Vbseff - Delt_vth - T2 + (paramPtr->k3
      + paramPtr->k3b * Vbseff) * Vth_NarrowW + T1 - DIBL_Sft;

  dVth_dVb = Lpe_Vb * paramPtr->k1ox * dsqrtPhis_dVb - k2ox
           - dDelt_vth_dVb - dT2_dVb + paramPtr->k3b * Vth_NarrowW
           - paramPtr->etab * Vds * paramPtr->theta0vb0 * T4
           + paramPtr->kt2 * TempRatio;
  dVth_dVd = -dDIBL_Sft_dVd;


  // Calculate n
  tmp1 = epssub / Xdep;
  nstar = model_.vtm / CONSTQ * (model_.coxe + tmp1 + paramPtr->cit);
  tmp2 = paramPtr->nfactor * tmp1;
  tmp3 = paramPtr->cdsc + paramPtr->cdscb * Vbseff
       + paramPtr->cdscd * Vds;
  tmp4 = (tmp2 + tmp3 * Theta0 + paramPtr->cit) / model_.coxe;
  if (tmp4 >= -0.5)
  {
    n = 1.0 + tmp4;
    dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
           + paramPtr->cdscb * Theta0) / model_.coxe;
    dn_dVd = paramPtr->cdscd * Theta0 / model_.coxe;
  }
  else
  {
    T0 = 1.0 / (3.0 + 8.0 * tmp4);
    n = (1.0 + 3.0 * tmp4) * T0;
    T0 *= T0;
    dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
           + paramPtr->cdscb * Theta0) / model_.coxe * T0;
    dn_dVd = paramPtr->cdscd * Theta0 / model_.coxe * T0;
  }


  // Vth correction for Pocket implant
  if (paramPtr->dvtp0 > 0.0)
  {
    T0 = -paramPtr->dvtp1 * Vds;
    if (T0 < -CONSTEXP_THRESHOLD)
    {
      T2 = CONSTMIN_EXP;
      dT2_dVd = 0.0;
    }
    else
    {
      T2 = exp(T0);
      dT2_dVd = -paramPtr->dvtp1 * T2;
    }

    T3 = Leff + paramPtr->dvtp0 * (1.0 + T2);
    dT3_dVd = paramPtr->dvtp0 * dT2_dVd;
    if (model_.tempMod < 2)
    {
      T4 = Vtm * log(Leff / T3);
      dT4_dVd = -Vtm * dT3_dVd / T3;
    }
    else
    {
      T4 = model_.vtm0 * log(Leff / T3);
      dT4_dVd = -model_.vtm0 * dT3_dVd / T3;
    }
    dDITS_Sft_dVd = dn_dVd * T4 + n * dT4_dVd;
    dDITS_Sft_dVb = T4 * dn_dVb;

    Vth -= n * T4;
    dVth_dVd -= dDITS_Sft_dVd;
    dVth_dVb -= dDITS_Sft_dVb;
  }
  von = Vth;

  // Poly Gate Si Depletion Effect
  T0 = vfb + paramPtr->phi;
  if(model_.mtrlMod == 0)
    T1 = CONSTEPSSI;
  else
    T1 = model_.epsrgate * CONSTEPS0;

  polyDepletion(T0, paramPtr->ngate, T1, model_.coxe,
                vgs, vgs_eff, dvgs_eff_dvg);

  polyDepletion(T0, paramPtr->ngate, T1, model_.coxe,
                vgd, vgd_eff, dvgd_eff_dvg);

  if(mode>0)
  {
    Vgs_eff = vgs_eff;
    dVgs_eff_dVg = dvgs_eff_dvg;
  }
  else
  {
    Vgs_eff = vgd_eff;
    dVgs_eff_dVg = dvgd_eff_dvg;
  }

  Vgst = Vgs_eff - Vth;

  // Calculate Vgsteff
  T0 = n * Vtm;
  T1 = paramPtr->mstar * Vgst;
  T2 = T1 / T0;
  if (T2 > CONSTEXP_THRESHOLD)
  {
    T10 = T1;
    dT10_dVg = paramPtr->mstar * dVgs_eff_dVg;
    dT10_dVd = -dVth_dVd * paramPtr->mstar;
    dT10_dVb = -dVth_dVb * paramPtr->mstar;
  }
  else if (T2 < -CONSTEXP_THRESHOLD)
  {
    T10 = Vtm * log(1.0 + CONSTMIN_EXP);
    dT10_dVg = 0.0;
    dT10_dVd = T10 * dn_dVd;
    dT10_dVb = T10 * dn_dVb;
    T10 *= n;
  }
  else
  {
    ExpVgst = exp(T2);
    T3 = Vtm * log(1.0 + ExpVgst);
    T10 = n * T3;
    dT10_dVg = paramPtr->mstar * ExpVgst / (1.0 + ExpVgst);
    dT10_dVb = T3 * dn_dVb - dT10_dVg * (dVth_dVb + Vgst * dn_dVb / n);
    dT10_dVd = T3 * dn_dVd - dT10_dVg * (dVth_dVd + Vgst * dn_dVd / n);
    dT10_dVg *= dVgs_eff_dVg;
  }

  T1 = paramPtr->voffcbn - (1.0 - paramPtr->mstar) * Vgst;
  T2 = T1 / T0;

  if (T2 < -CONSTEXP_THRESHOLD)
  {
    T3 = model_.coxe * CONSTMIN_EXP / paramPtr->cdep0;
    T9 = paramPtr->mstar + T3 * n;
    dT9_dVg = 0.0;
    dT9_dVd = dn_dVd * T3;
    dT9_dVb = dn_dVb * T3;
  }
  else if (T2 > CONSTEXP_THRESHOLD)
  {
    T3 = model_.coxe * CONSTMAX_EXP / paramPtr->cdep0;
    T9 = paramPtr->mstar + T3 * n;
    dT9_dVg = 0.0;
    dT9_dVd = dn_dVd * T3;
    dT9_dVb = dn_dVb * T3;
  }
  else
  {
    ExpVgst = exp(T2);
    T3 = model_.coxe / paramPtr->cdep0;
    T4 = T3 * ExpVgst;
    T5 = T1 * T4 / T0;
    T9 = paramPtr->mstar + n * T4;
    dT9_dVg = T3 * (paramPtr->mstar - 1.0) * ExpVgst / Vtm;
    dT9_dVb = T4 * dn_dVb - dT9_dVg * dVth_dVb - T5 * dn_dVb;
    dT9_dVd = T4 * dn_dVd - dT9_dVg * dVth_dVd - T5 * dn_dVd;
    dT9_dVg *= dVgs_eff_dVg;
  }

  Vgsteff = T10 / T9;
  Vgsteff_forNoise = Vgsteff;
  T11 = T9 * T9;
  dVgsteff_dVg = (T9 * dT10_dVg - T10 * dT9_dVg) / T11;
  dVgsteff_dVd = (T9 * dT10_dVd - T10 * dT9_dVd) / T11;
  dVgsteff_dVb = (T9 * dT10_dVb - T10 * dT9_dVb) / T11;

  // Calculate Effective Channel Geometry
  T9 = sqrtPhis - paramPtr->sqrtPhi;
  Weff = paramPtr->weff - 2.0 * (paramPtr->dwg * Vgsteff
      + paramPtr->dwb * T9);
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

  if (model_.rdsMod == 1)
  {
    Rds = dRds_dVg = dRds_dVb = 0.0;
  }
  else
  {
    T0 = 1.0 + paramPtr->prwg * Vgsteff;
    dT0_dVg = -paramPtr->prwg / T0 / T0;
    T1 = paramPtr->prwb * T9;
    dT1_dVb = paramPtr->prwb * dsqrtPhis_dVb;

    T2 = 1.0 / T0 + T1;
    T3 = T2 + sqrt(T2 * T2 + 0.01); // 0.01 = 4.0 * 0.05 * 0.05
    dT3_dVg = 1.0 + T2 / (T3 - T2);
    dT3_dVb = dT3_dVg * dT1_dVb;
    dT3_dVg *= dT0_dVg;

    T4 = paramPtr->rds0 * 0.5;
    Rds = paramPtr->rdswmin + T3 * T4;
    dRds_dVg = T4 * dT3_dVg;
    dRds_dVb = T4 * dT3_dVb;

    if (Rds > 0.0)
    {
     grdsw = 1.0 / Rds;
    }
    else
    {
     grdsw = 0.0;
    }
  }

  // Calculate Abulk
  T9 = 0.5 * paramPtr->k1ox * Lpe_Vb / sqrtPhis;
  T1 = T9 + k2ox - paramPtr->k3b * Vth_NarrowW;
  dT1_dVb = -T9 / sqrtPhis * dsqrtPhis_dVb;

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

  if (Abulk < 0.1)
  {
    T9 = 1.0 / (3.0 - 20.0 * Abulk);
    Abulk = (0.2 - Abulk) * T9;
    T10 = T9 * T9;
    dAbulk_dVb *= T10;
    dAbulk_dVg *= T10;
  }

  Abulk_forNoise = Abulk; // ERK. this is a bit screwy.  But in spice3/ngspice, Abulk is saved as an instance variable and used later in noise calculations.  But, within the load function the "local" copy of Abulk gets modified further afterwards.  So it is then "wrong" for noise at that point.  As we've made Abulk a class variable, there needs to be an extra copy for noise, saved at the right time.

  T2 = paramPtr->keta * Vbseff;
  if (T2 >= -0.9)
  {
    T0 = 1.0 / (1.0 + T2);
    dT0_dVb = -paramPtr->keta * T0 * T0;
  }
  else
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
  if (model_.mtrlMod)
    T14 = 2.0 * model_.dtype *(model_.phig - model_.easub - 0.5*model_.Eg0 + 0.45);
  else
    T14 = 0.0;

  if (model_.mobMod == 0)
  { T0 = Vgsteff + Vth + Vth - T14;
    T2 = paramPtr->ua + paramPtr->uc * Vbseff;
    T3 = T0 / toxe;
    T12 = sqrt(Vth * Vth + 0.0001);
    T9 = 1.0/(Vgsteff + 2*T12);
    T10 = T9*toxe;
    T8 = paramPtr->ud * T10 * T10 * Vth;
    T6 = T8 * Vth;
    T5 = T3 * (T2 + paramPtr->ub * T3) + T6;
    T7 = - 2.0 * T6 * T9;
    T11 = T7 * Vth/T12;
    dDenomi_dVg = (T2 + 2.0 * paramPtr->ub * T3) / toxe;
    T13 = 2.0 * (dDenomi_dVg + T11 + T8);
    dDenomi_dVd = T13 * dVth_dVd;
    dDenomi_dVb = T13 * dVth_dVb + paramPtr->uc * T3;
    dDenomi_dVg+= T7;
  }
  else if (model_.mobMod == 1)
  {   T0 = Vgsteff + Vth + Vth - T14;
    T2 = 1.0 + paramPtr->uc * Vbseff;
    T3 = T0 / toxe;
    T4 = T3 * (paramPtr->ua + paramPtr->ub * T3);
    T12 = sqrt(Vth * Vth + 0.0001);
    T9 = 1.0/(Vgsteff + 2*T12);
    T10 = T9*toxe;
    T8 = paramPtr->ud * T10 * T10 * Vth;
    T6 = T8 * Vth;
    T5 = T4 * T2 + T6;
    T7 = - 2.0 * T6 * T9;
    T11 = T7 * Vth/T12;
    dDenomi_dVg = (paramPtr->ua + 2.0 * paramPtr->ub * T3) * T2 / toxe;
    T13 = 2.0 * (dDenomi_dVg + T11 + T8);
    dDenomi_dVd = T13 * dVth_dVd;
    dDenomi_dVb = T13 * dVth_dVb + paramPtr->uc * T4;
    dDenomi_dVg+= T7;
  }
  else
  {   T0 = (Vgsteff + vtfbphi1) / toxe;
    T1 = exp(paramPtr->eu * log(T0));
    dT1_dVg = T1 * paramPtr->eu / T0 / toxe;
    T2 = paramPtr->ua + paramPtr->uc * Vbseff;
    T3 = T0 / toxe;
    T12 = sqrt(Vth * Vth + 0.0001);
    T9 = 1.0/(Vgsteff + 2*T12);
    T10 = T9*toxe;
    T8 = paramPtr->ud * T10 * T10 * Vth;
    T6 = T8 * Vth;
    T5 = T1 * T2 + T6;
    T7 = - 2.0 * T6 * T9;
    T11 = T7 * Vth/T12;
    dDenomi_dVg = T2 * dT1_dVg + T7;
    T13 = 2.0 * (T11 + T8);
    dDenomi_dVd = T13 * dVth_dVd;
    dDenomi_dVb = T13 * dVth_dVb + T1 * paramPtr->uc;
  }

  if (T5 >= -0.8)
  {
    Denomi = 1.0 + T5;
  }
  else
  {
    T9 = 1.0 / (7.0 + 10.0 * T5);
    Denomi = (0.6 + T5) * T9;
    T9 *= T9;
    dDenomi_dVg *= T9;
    dDenomi_dVd *= T9;
    dDenomi_dVb *= T9;
  }

  ueff = ueff = u0temp / Denomi;
  T9 = -ueff / Denomi;
  dueff_dVg = T9 * dDenomi_dVg;
  dueff_dVd = T9 * dDenomi_dVd;
  dueff_dVb = T9 * dDenomi_dVb;

  // Saturation Drain Voltage  Vdsat
  WVCox = Weff * vsattemp * model_.coxe;
  WVCoxRds = WVCox * Rds;

  Esat = 2.0 * vsattemp / ueff;
  EsatL = EsatL = Esat * Leff;
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
  else if (a1 > 0.0)
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

    dT3_dVg = (T1 * dT1_dVg - 2.0 * (T0 * dT2_dVg + T2 * dT0_dVg))
     / T3;
    dT3_dVd = (T1 * dT1_dVd - 2.0 * (T0 * dT2_dVd + T2 * dT0_dVd))
     / T3;
    dT3_dVb = (T1 * dT1_dVb - 2.0 * (T0 * dT2_dVb + T2 * dT0_dVb))
     / T3;

    dVdsat_dVg = (dT1_dVg - (T1 * dT1_dVg - dT0_dVg * T2
      - T0 * dT2_dVg) / T3 - Vdsat * dT0_dVg) / T0;
    dVdsat_dVb = (dT1_dVb - (T1 * dT1_dVb - dT0_dVb * T2
      - T0 * dT2_dVb) / T3 - Vdsat * dT0_dVb) / T0;
    dVdsat_dVd = (dT1_dVd - (T1 * dT1_dVd - T0 * dT2_dVd) / T3) / T0;
  }
  vdsat = Vdsat;

  // Calculate Vdseff
  T1 = Vdsat - Vds - paramPtr->delta;
  dT1_dVg = dVdsat_dVg;
  dT1_dVd = dVdsat_dVd - 1.0;
  dT1_dVb = dVdsat_dVb;

  T2 = sqrt(T1 * T1 + 4.0 * paramPtr->delta * Vdsat);
  T0 = T1 / T2;
  T9 = 2.0 * paramPtr->delta;
  T3 = T9 / T2;
  dT2_dVg = T0 * dT1_dVg + T3 * dVdsat_dVg;
  dT2_dVd = T0 * dT1_dVd + T3 * dVdsat_dVd;
  dT2_dVb = T0 * dT1_dVb + T3 * dVdsat_dVb;

  if (T1 >= 0.0)
  {
    Vdseff = Vdsat - 0.5 * (T1 + T2);
    dVdseff_dVg = dVdsat_dVg - 0.5 * (dT1_dVg + dT2_dVg);
    dVdseff_dVd = dVdsat_dVd - 0.5 * (dT1_dVd + dT2_dVd);
    dVdseff_dVb = dVdsat_dVb - 0.5 * (dT1_dVb + dT2_dVb);
  }
  else
  {
    T4 = T9 / (T2 - T1);
    T5 = 1.0 - T4;
    T6 = Vdsat * T4 / (T2 - T1);
    Vdseff = Vdsat * T5;
    dVdseff_dVg = dVdsat_dVg * T5 + T6 * (dT2_dVg - dT1_dVg);
    dVdseff_dVd = dVdsat_dVd * T5 + T6 * (dT2_dVd - dT1_dVd);
    dVdseff_dVb = dVdsat_dVb * T5 + T6 * (dT2_dVb - dT1_dVb);
  }

  if (Vds == 0.0)
  {
    Vdseff = 0.0;
    dVdseff_dVg = 0.0;
    dVdseff_dVb = 0.0;
  }

  if (Vdseff > Vds)
  {
    Vdseff = Vds;
  }

  diffVds = Vds - Vdseff;
  Vdseff_forNoise = Vdseff;

  // Velocity Overshoot
  if((model_.lambdaGiven) && (model_.lambda > 0.0) )
  {
    T1 =  Leff * ueff;
    T2 = paramPtr->lambda / T1;
    T3 = -T2 / T1 * Leff;
    dT2_dVd = T3 * dueff_dVd;
    dT2_dVg = T3 * dueff_dVg;
    dT2_dVb = T3 * dueff_dVb;
    T5 = 1.0 / (Esat * paramPtr->litl);
    T4 = -T5 / EsatL;
    dT5_dVg = dEsatL_dVg * T4;
    dT5_dVd = dEsatL_dVd * T4;
    dT5_dVb = dEsatL_dVb * T4;
    T6 = 1.0 + diffVds  * T5;
    dT6_dVg = dT5_dVg * diffVds - dVdseff_dVg * T5;
    dT6_dVd = dT5_dVd * diffVds + (1.0 - dVdseff_dVd) * T5;
    dT6_dVb = dT5_dVb * diffVds - dVdseff_dVb * T5;
    T7 = 2.0 / (T6 * T6 + 1.0);
    T8 = 1.0 - T7;
    T9 = T6 * T7 * T7;
    dT8_dVg = T9 * dT6_dVg;
    dT8_dVd = T9 * dT6_dVd;
    dT8_dVb = T9 * dT6_dVb;
    T10 = 1.0 + T2 * T8;
    dT10_dVg = dT2_dVg * T8 + T2 * dT8_dVg;
    dT10_dVd = dT2_dVd * T8 + T2 * dT8_dVd;
    dT10_dVb = dT2_dVb * T8 + T2 * dT8_dVb;
    if(T10 == 1.0)
    {
      dT10_dVg = dT10_dVd = dT10_dVb = 0.0;
    }

    dEsatL_dVg *= T10;
    dEsatL_dVg += EsatL * dT10_dVg;
    dEsatL_dVd *= T10;
    dEsatL_dVd += EsatL * dT10_dVd;
    dEsatL_dVb *= T10;
    dEsatL_dVb += EsatL * dT10_dVb;
    EsatL *= T10;
  }

  // Calculate Vasat
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

  // Calculate Idl first
  tmp1 = vtfbphi2;
  tmp2 = 2.0e8 * model_.toxp;
  dT0_dVg = 1.0 / tmp2;
  T0 = (Vgsteff + tmp1) * dT0_dVg;

  tmp3 = exp(model_.bdos * 0.7 * log(T0));
  T1 = 1.0 + tmp3;
  T2 = model_.bdos * 0.7 * tmp3 / T0;
  Tcen = model_.ados * 1.9e-9 / T1;
  dTcen_dVg = -Tcen * T2 * dT0_dVg / T1;

  Coxeff = epssub * model_.coxp / (epssub + model_.coxp * Tcen);
  dCoxeff_dVg = -Coxeff * Coxeff * dTcen_dVg / epssub;

  CoxeffWovL = Coxeff * Weff / Leff;
  beta = ueff * CoxeffWovL;
  T3 = ueff / Leff;
  dbeta_dVg = CoxeffWovL * dueff_dVg + T3
          * (Weff * dCoxeff_dVg + Coxeff * dWeff_dVg);
  dbeta_dVd = CoxeffWovL * dueff_dVd;
  dbeta_dVb = CoxeffWovL * dueff_dVb + T3 * Coxeff * dWeff_dVb;

  AbovVgst2Vtm = Abulk / Vgst2Vtm;
  T0 = 1.0 - 0.5 * Vdseff * AbovVgst2Vtm;
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
  Idl = gche / T0;
  T1 = (1.0 - Idl * Rds) / T0;
  T2 = Idl * Idl;
  dIdl_dVg = T1 * dgche_dVg - T2 * dRds_dVg;
  dIdl_dVd = T1 * dgche_dVd;
  dIdl_dVb = T1 * dgche_dVb - T2 * dRds_dVb;

  // Calculate degradation factor due to pocket implant
  if (paramPtr->fprout <= 0.0)
  {
    FP = 1.0;
    dFP_dVg = 0.0;
  }
  else
  {
    T9 = paramPtr->fprout * sqrt(Leff) / Vgst2Vtm;
    FP = 1.0 / (1.0 + T9);
    dFP_dVg = FP * FP * T9 / Vgst2Vtm;
  }

  // Calculate VACLM
  T8 = paramPtr->pvag / EsatL;
  T9 = T8 * Vgsteff;
  if (T9 > -0.9)
  {
    PvagTerm = 1.0 + T9;
    dPvagTerm_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL);
    dPvagTerm_dVb = -T9 * dEsatL_dVb / EsatL;
    dPvagTerm_dVd = -T9 * dEsatL_dVd / EsatL;
  }
  else
  {
    T4 = 1.0 / (17.0 + 20.0 * T9);
    PvagTerm = (0.8 + T9) * T4;
    T4 *= T4;
    dPvagTerm_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL) * T4;
    T9 *= T4 / EsatL;
    dPvagTerm_dVb = -T9 * dEsatL_dVb;
    dPvagTerm_dVd = -T9 * dEsatL_dVd;
  }

  if ((paramPtr->pclm > CONSTMIN_EXP) && (diffVds > 1.0e-10))
  {
    T0 = 1.0 + Rds * Idl;
    dT0_dVg = dRds_dVg * Idl + Rds * dIdl_dVg;
    dT0_dVd = Rds * dIdl_dVd;
    dT0_dVb = dRds_dVb * Idl + Rds * dIdl_dVb;

    T2 = Vdsat / Esat;
    T1 = Leff + T2;
    dT1_dVg = (dVdsat_dVg - T2 * dEsatL_dVg / Leff) / Esat;
    dT1_dVd = (dVdsat_dVd - T2 * dEsatL_dVd / Leff) / Esat;
    dT1_dVb = (dVdsat_dVb - T2 * dEsatL_dVb / Leff) / Esat;

    Cclm = FP * PvagTerm * T0 * T1 / (paramPtr->pclm * paramPtr->litl);
    dCclm_dVg = Cclm * (dFP_dVg / FP + dPvagTerm_dVg / PvagTerm
              + dT0_dVg / T0 + dT1_dVg / T1);
    dCclm_dVb = Cclm * (dPvagTerm_dVb / PvagTerm + dT0_dVb / T0
              + dT1_dVb / T1);
    dCclm_dVd = Cclm * (dPvagTerm_dVd / PvagTerm + dT0_dVd / T0
              + dT1_dVd / T1);
    VACLM = Cclm * diffVds;

    dVACLM_dVg = dCclm_dVg * diffVds - dVdseff_dVg * Cclm;
    dVACLM_dVb = dCclm_dVb * diffVds - dVdseff_dVb * Cclm;
    dVACLM_dVd = dCclm_dVd * diffVds + (1.0 - dVdseff_dVd) * Cclm;
  }
  else
  {
    VACLM = Cclm = CONSTMAX_EXP;
    dVACLM_dVd = dVACLM_dVg = dVACLM_dVb = 0.0;
    dCclm_dVd = dCclm_dVg = dCclm_dVb = 0.0;
  }

  // Calculate VADIBL
  if (paramPtr->thetaRout > CONSTMIN_EXP)
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
    else
    {
      T4 = 1.0 / (0.8 + T7);
      T3 = (17.0 + 20.0 * T7) * T4;
      dVADIBL_dVg *= T3;
      dVADIBL_dVb = dVADIBL_dVb * T3
        - VADIBL * paramPtr->pdiblb * T4 * T4;
      dVADIBL_dVd *= T3;
      VADIBL *= T3;
    }

    dVADIBL_dVg = dVADIBL_dVg * PvagTerm + VADIBL * dPvagTerm_dVg;
    dVADIBL_dVb = dVADIBL_dVb * PvagTerm + VADIBL * dPvagTerm_dVb;
    dVADIBL_dVd = dVADIBL_dVd * PvagTerm + VADIBL * dPvagTerm_dVd;
    VADIBL *= PvagTerm;
  }
  else
  {
    VADIBL = CONSTMAX_EXP;
    dVADIBL_dVd = dVADIBL_dVg = dVADIBL_dVb = 0.0;
  }

  // Calculate Va
  Va = Vasat + VACLM;
  dVa_dVg = dVasat_dVg + dVACLM_dVg;
  dVa_dVb = dVasat_dVb + dVACLM_dVb;
  dVa_dVd = dVasat_dVd + dVACLM_dVd;

  // Calculate VADITS
  T0 = paramPtr->pditsd * Vds;
  if (T0 > CONSTEXP_THRESHOLD)
  {
    T1 = CONSTMAX_EXP;
    dT1_dVd = 0;
  }
  else
  {
    T1 = exp(T0);
    dT1_dVd = T1 * paramPtr->pditsd;
  }

  if (paramPtr->pdits > CONSTMIN_EXP)
  {
    T2 = 1.0 + model_.pditsl * Leff;
    VADITS = (1.0 + T2 * T1) / paramPtr->pdits;
    dVADITS_dVg = VADITS * dFP_dVg;
    dVADITS_dVd = FP * T2 * dT1_dVd / paramPtr->pdits;
    VADITS *= FP;
  }
  else
  {
    VADITS = CONSTMAX_EXP;
    dVADITS_dVg = dVADITS_dVd = 0;
  }

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

  // Add DIBL to Ids
  T9 = diffVds / VADIBL;
  T0 = 1.0 + T9;
  Idsa = Idl * T0;
  dIdsa_dVg = T0 * dIdl_dVg - Idl * (dVdseff_dVg + T9 * dVADIBL_dVg) / VADIBL;
  dIdsa_dVd = T0 * dIdl_dVd + Idl
          * (1.0 - dVdseff_dVd - T9 * dVADIBL_dVd) / VADIBL;
  dIdsa_dVb = T0 * dIdl_dVb - Idl * (dVdseff_dVb + T9 * dVADIBL_dVb) / VADIBL;

  // Add DITS to Ids
  T9 = diffVds / VADITS;
  T0 = 1.0 + T9;
  dIdsa_dVg = T0 * dIdsa_dVg - Idsa * (dVdseff_dVg + T9 * dVADITS_dVg) / VADITS;
  dIdsa_dVd = T0 * dIdsa_dVd + Idsa * (1.0 - dVdseff_dVd - T9 * dVADITS_dVd) / VADITS;
  dIdsa_dVb = T0 * dIdsa_dVb - Idsa * dVdseff_dVb / VADITS;
  Idsa *= T0;

  // Add CLM to Ids
  T0 = log(Va / Vasat);
  dT0_dVg = dVa_dVg / Va - dVasat_dVg / Vasat;
  dT0_dVb = dVa_dVb / Va - dVasat_dVb / Vasat;
  dT0_dVd = dVa_dVd / Va - dVasat_dVd / Vasat;
  T1 = T0 / Cclm;
  T9 = 1.0 + T1;
  dT9_dVg = (dT0_dVg - T1 * dCclm_dVg) / Cclm;
  dT9_dVb = (dT0_dVb - T1 * dCclm_dVb) / Cclm;
  dT9_dVd = (dT0_dVd - T1 * dCclm_dVd) / Cclm;

  dIdsa_dVg = dIdsa_dVg * T9 + Idsa * dT9_dVg;
  dIdsa_dVb = dIdsa_dVb * T9 + Idsa * dT9_dVb;
  dIdsa_dVd = dIdsa_dVd * T9 + Idsa * dT9_dVd;
  Idsa *= T9;

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
    T4 = Idsa * Vdseff;
    Isub = T1 * T4;
    Gbg = T1 * (dIdsa_dVg * Vdseff + Idsa * dVdseff_dVg)
        + T4 * dT1_dVg;
    Gbd = T1 * (dIdsa_dVd * Vdseff + Idsa * dVdseff_dVd)
        + T4 * dT1_dVd;
    Gbb = T1 * (dIdsa_dVb * Vdseff + Idsa * dVdseff_dVb)
        + T4 * dT1_dVb;

    Gbd += Gbg * dVgsteff_dVd;
    Gbb += Gbg * dVgsteff_dVb;
    Gbg *= dVgsteff_dVg;
    Gbb *= dVbseff_dVb;
  }
  csub = Isub;
  gbbs = Gbb;
  gbgs = Gbg;
  gbds = Gbd;

  // Add SCBE to Ids
  T9 = diffVds / VASCBE;
  T0 = 1.0 + T9;
  Ids = Idsa * T0;

  Gm = T0 * dIdsa_dVg - Idsa
   * (dVdseff_dVg + T9 * dVASCBE_dVg) / VASCBE;
  Gds = T0 * dIdsa_dVd + Idsa
    * (1.0 - dVdseff_dVd - T9 * dVASCBE_dVd) / VASCBE;
  Gmb = T0 * dIdsa_dVb - Idsa
    * (dVdseff_dVb + T9 * dVASCBE_dVb) / VASCBE;


  tmp1 = Gds + Gm * dVgsteff_dVd;
  tmp2 = Gmb + Gm * dVgsteff_dVb;
  tmp3 = Gm;

  Gm = (Ids * dVdseff_dVg + Vdseff * tmp3) * dVgsteff_dVg;
  Gds = Ids * (dVdseff_dVd + dVdseff_dVg * dVgsteff_dVd)
    + Vdseff * tmp1;
  Gmb = (Ids * (dVdseff_dVb + dVdseff_dVg * dVgsteff_dVb)
    + Vdseff * tmp2) * dVbseff_dVb;

  cdrain = Ids * Vdseff;

  // Source End Velocity Limit
  if((model_.vtlGiven) && (model_.vtl > 0.0) )
  {
    T12 = 1.0 / Leff / CoxeffWovL;
    T11 = T12 / Vgsteff;
    T10 = -T11 / Vgsteff;
    vs = cdrain * T11; // vs
    dvs_dVg = Gm * T11 + cdrain * T10 * dVgsteff_dVg;
    dvs_dVd = Gds * T11 + cdrain * T10 * dVgsteff_dVd;
    dvs_dVb = Gmb * T11 + cdrain * T10 * dVgsteff_dVb;
    T0 = 2 * CONSTMM;
    T1 = vs / (paramPtr->vtl * paramPtr->tfactor);

    if(T1 > 0.0)
    {
      T2 = 1.0 + exp(T0 * log(T1));
      T3 = (T2 - 1.0) * T0 / vs;
      Fsevl = 1.0 / exp(log(T2)/ T0);
      dT2_dVg = T3 * dvs_dVg;
      dT2_dVd = T3 * dvs_dVd;
      dT2_dVb = T3 * dvs_dVb;
      T4 = -1.0 / T0 * Fsevl / T2;
      dFsevl_dVg = T4 * dT2_dVg;
      dFsevl_dVd = T4 * dT2_dVd;
      dFsevl_dVb = T4 * dT2_dVb;
    }
    else
    {
      Fsevl = 1.0;
      dFsevl_dVg = 0.0;
      dFsevl_dVd = 0.0;
      dFsevl_dVb = 0.0;
    }
    Gm *=Fsevl;
    Gm += cdrain * dFsevl_dVg;
    Gmb *=Fsevl;
    Gmb += cdrain * dFsevl_dVb;
    Gds *=Fsevl;
    Gds += cdrain * dFsevl_dVd;

    cdrain *= Fsevl;
  }

  gds = Gds;
  gm = Gm;
  gmbs = Gmb;
  IdovVds = Ids;
  if( IdovVds <= 1.0e-9) IdovVds = 1.0e-9;

  // Calculate Rg
  if ((rgateMod > 1) || (trnqsMod != 0) || (acnqsMod != 0))
  {
    T9 = paramPtr->xrcrg2 * model_.vtm;
    T0 = T9 * beta;
    dT0_dVd = (dbeta_dVd + dbeta_dVg * dVgsteff_dVd) * T9;
    dT0_dVb = (dbeta_dVb + dbeta_dVg * dVgsteff_dVb) * T9;
    dT0_dVg = dbeta_dVg * T9;

    gcrg = paramPtr->xrcrg1 * ( T0 + Ids);
    gcrgd = paramPtr->xrcrg1 * (dT0_dVd + tmp1);
    gcrgb = paramPtr->xrcrg1 * (dT0_dVb + tmp2)
                     * dVbseff_dVb;
    gcrgg = paramPtr->xrcrg1 * (dT0_dVg + tmp3)
                     * dVgsteff_dVg;

    if (nf != 1.0)
    {
      gcrg *= nf;
      gcrgg *= nf;
      gcrgd *= nf;
      gcrgb *= nf;
    }

    if (rgateMod == 2)
    {
      T10 = grgeltd * grgeltd;
      T11 = grgeltd + gcrg;
      gcrg = grgeltd * gcrg / T11;
      T12 = T10 / T11 / T11;
      gcrgg *= T12;
      gcrgd *= T12;
      gcrgb *= T12;
    }
    gcrgs = -(gcrgg + gcrgd + gcrgb);
  }

  // Calculate bias-dependent external S/D resistance
  if (model_.rdsMod)
  {   // Rs(V)
    T0 = vgs - paramPtr->vfbsd;
    T1 = sqrt(T0 * T0 + 1.0e-4);
    vgs_eff = 0.5 * (T0 + T1);
    dvgs_eff_dvg = vgs_eff / T1;

    T0 = 1.0 + paramPtr->prwg * vgs_eff;
    dT0_dvg = -paramPtr->prwg / T0 / T0 * dvgs_eff_dvg;
    T1 = -paramPtr->prwb * vbs;
    dT1_dvb = -paramPtr->prwb;

    T2 = 1.0 / T0 + T1;
    T3 = T2 + sqrt(T2 * T2 + 0.01);
    dT3_dvg = T3 / (T3 - T2);
    dT3_dvb = dT3_dvg * dT1_dvb;
    dT3_dvg *= dT0_dvg;

    T4 = paramPtr->rs0 * 0.5;
    Rs = paramPtr->rswmin + T3 * T4;
    dRs_dvg = T4 * dT3_dvg;
    dRs_dvb = T4 * dT3_dvb;

    T0 = 1.0 + sourceConductance * Rs;
    gstot = sourceConductance / T0;
    T0 = -gstot * gstot;
    dgstot_dvd = 0.0; // place holder
    dgstot_dvg = T0 * dRs_dvg;
    dgstot_dvb = T0 * dRs_dvb;
    dgstot_dvs = -(dgstot_dvg + dgstot_dvb + dgstot_dvd);

    // Rd(V)
    T0 = vgd - paramPtr->vfbsd;
    T1 = sqrt(T0 * T0 + 1.0e-4);
    vgd_eff = 0.5 * (T0 + T1);
    dvgd_eff_dvg = vgd_eff / T1;

    T0 = 1.0 + paramPtr->prwg * vgd_eff;
    dT0_dvg = -paramPtr->prwg / T0 / T0 * dvgd_eff_dvg;
    T1 = -paramPtr->prwb * vbd;
    dT1_dvb = -paramPtr->prwb;

    T2 = 1.0 / T0 + T1;
    T3 = T2 + sqrt(T2 * T2 + 0.01);
    dT3_dvg = T3 / (T3 - T2);
    dT3_dvb = dT3_dvg * dT1_dvb;
    dT3_dvg *= dT0_dvg;

    T4 = paramPtr->rd0 * 0.5;
    Rd = paramPtr->rdwmin + T3 * T4;
    dRd_dvg = T4 * dT3_dvg;
    dRd_dvb = T4 * dT3_dvb;

    T0 = 1.0 + drainConductance * Rd;
    gdtot = drainConductance / T0;
    T0 = -gdtot * gdtot;
    dgdtot_dvs = 0.0;
    dgdtot_dvg = T0 * dRd_dvg;
    dgdtot_dvb = T0 * dRd_dvb;
    dgdtot_dvd = -(dgdtot_dvg + dgdtot_dvb + dgdtot_dvs);

    gstotd = vses * dgstot_dvd;
    gstotg = vses * dgstot_dvg;
    gstots = vses * dgstot_dvs;
    gstotb = vses * dgstot_dvb;

    T2 = vdes - vds;
    gdtotd = T2 * dgdtot_dvd;
    gdtotg = T2 * dgdtot_dvg;
    gdtots = T2 * dgdtot_dvs;
    gdtotb = T2 * dgdtot_dvb;
  }
  else // WDLiu: for bypass
  {
    gstot = gstotd = gstotg = 0.0;
    gstots = gstotb = 0.0;
    gdtot = gdtotd = gdtotg = 0.0;
    gdtots = gdtotb = 0.0;
  }

  // Calculate GIDL current
  if (model_.mtrlMod == 0)
  {
    T0 = 3.0 * toxe;
  }
  else
  {
    T0 = model_.epsrsub * toxe / epsrox;
  }

  if(model_.mtrlMod ==0)
    T1 = (vds - vgs_eff - paramPtr->egidl ) / T0;
  else
    T1 = (vds - vgs_eff - paramPtr->egidl + paramPtr->vfbsd) / T0;

  if ((paramPtr->agidl <= 0.0) || (paramPtr->bgidl <= 0.0)
                || (T1 <= 0.0) || (paramPtr->cgidl <= 0.0) || (vbd > 0.0))
  {
    Igidl = Ggidld = Ggidlg = Ggidlb = 0.0;
  }
  else
  {
    dT1_dVd = 1.0 / T0;
    dT1_dVg = -dvgs_eff_dvg * dT1_dVd;
    T2 = paramPtr->bgidl / T1;
    if (T2 < 100.0)
    {
      Igidl = paramPtr->agidl * paramPtr->weffCJ * T1 * exp(-T2);
      T3 = Igidl * (1.0 + T2) / T1;
      Ggidld = T3 * dT1_dVd;
      Ggidlg = T3 * dT1_dVg;
    }
    else
    {
      Igidl = paramPtr->agidl * paramPtr->weffCJ * 3.720075976e-44;
      Ggidld = Igidl * dT1_dVd;
      Ggidlg = Igidl * dT1_dVg;
      Igidl *= T1;
    }

    T4 = vbd * vbd;
    T5 = -vbd * T4;
    T6 = paramPtr->cgidl + T5;
    T7 = T5 / T6;
    T8 = 3.0 * paramPtr->cgidl * T4 / T6 / T6;
    Ggidld = Ggidld * T7 + Igidl * T8;
    Ggidlg = Ggidlg * T7;
    Ggidlb = -Igidl * T8;
    Igidl *= T7;
  }
  ggidld = Ggidld;
  ggidlg = Ggidlg;
  ggidlb = Ggidlb;

  // Calculate GISL current

  if (model_.mtrlMod == 0)
  {
    T1 = (-vds - vgd_eff - paramPtr->egisl ) / T0;
  }
  else
  {
    T1 = (-vds - vgd_eff - paramPtr->egisl + paramPtr->vfbsd ) / T0;
  }

  if ((paramPtr->agisl <= 0.0) || (paramPtr->bgisl <= 0.0)
    || (T1 <= 0.0) || (paramPtr->cgisl <= 0.0) || (vbs > 0.0))
  {
    Igisl = Ggisls = Ggislg = Ggislb = 0.0;
  }
  else
  {
    dT1_dVd = 1.0 / T0;
    dT1_dVg = -dvgd_eff_dvg * dT1_dVd;
    T2 = paramPtr->bgisl / T1;
    if (T2 < 100.0)
    {
      Igisl = paramPtr->agisl * paramPtr->weffCJ * T1 * exp(-T2);
      T3 = Igisl * (1.0 + T2) / T1;
      Ggisls = T3 * dT1_dVd;
      Ggislg = T3 * dT1_dVg;
    }
    else
    {
      Igisl = paramPtr->agisl * paramPtr->weffCJ * 3.720075976e-44;
      Ggisls = Igisl * dT1_dVd;
      Ggislg = Igisl * dT1_dVg;
      Igisl *= T1;
    }

    T4 = vbs * vbs;
    T5 = -vbs * T4;
    T6 = paramPtr->cgisl + T5;
    T7 = T5 / T6;
    T8 = 3.0 * paramPtr->cgisl * T4 / T6 / T6;
    Ggisls = Ggisls * T7 + Igisl * T8;
    Ggislg = Ggislg * T7;
    Ggislb = -Igisl * T8;
    Igisl *= T7;
  }
  ggisls = Ggisls;
  ggislg = Ggislg;
  ggislb = Ggislb;


  // Calculate gate tunneling current
  if ((model_.igcMod != 0) || (model_.igbMod != 0))
  {
    Vfb = vfbzb;
    V3 = Vfb - Vgs_eff + Vbseff - CONSTDELTA_3;
    if (Vfb <= 0.0)
        T0 = sqrt(V3 * V3 - 4.0 * CONSTDELTA_3 * Vfb);
    else
        T0 = sqrt(V3 * V3 + 4.0 * CONSTDELTA_3 * Vfb);
    T1 = 0.5 * (1.0 + V3 / T0);
    Vfbeff = Vfb - 0.5 * (V3 + T0);
    dVfbeff_dVg = T1 * dVgs_eff_dVg;
    dVfbeff_dVb = -T1; // WDLiu: -No surprise? No. -Good!

    Voxacc = Vfb - Vfbeff;
    dVoxacc_dVg = -dVfbeff_dVg;
    dVoxacc_dVb = -dVfbeff_dVb;
    if (Voxacc < 0.0) // WDLiu: Avoiding numerical instability.
    {
      Voxacc = dVoxacc_dVg = dVoxacc_dVb = 0.0;
    }

    T0 = 0.5 * paramPtr->k1ox;
    T3 = Vgs_eff - Vfbeff - Vbseff - Vgsteff;
    if (paramPtr->k1ox == 0.0)
    {
      Voxdepinv = dVoxdepinv_dVg = dVoxdepinv_dVd = dVoxdepinv_dVb = 0.0;
    }
    else if (T3 < 0.0)
    {
      Voxdepinv = -T3;
      dVoxdepinv_dVg = -dVgs_eff_dVg + dVfbeff_dVg + dVgsteff_dVg;
      dVoxdepinv_dVd = dVgsteff_dVd;
      dVoxdepinv_dVb = dVfbeff_dVb + 1.0 + dVgsteff_dVb;
    }
    else
    {
      T1 = sqrt(T0 * T0 + T3);
      T2 = T0 / T1;
      Voxdepinv = paramPtr->k1ox * (T1 - T0);
      dVoxdepinv_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg);
      dVoxdepinv_dVd = -T2 * dVgsteff_dVd;
      dVoxdepinv_dVb = -T2 * (dVfbeff_dVb + 1.0 + dVgsteff_dVb);
    }

    Voxdepinv += Vgsteff;
    dVoxdepinv_dVg += dVgsteff_dVg;
    dVoxdepinv_dVd += dVgsteff_dVd;
    dVoxdepinv_dVb += dVgsteff_dVb;
  }

  if(model_.tempMod < 2)
  {
    tmp = Vtm;
  }
  else // model_.tempMod = 2
  {
    tmp = Vtm0;
  }

  if (model_.igcMod)
  {
    T0 = tmp * paramPtr->nigc;
    if(model_.igcMod == 1)
    {
      VxNVt = (Vgs_eff - model_.dtype * vth0) / T0;
      if (VxNVt > CONSTEXP_THRESHOLD)
      {
        Vaux = Vgs_eff - model_.dtype * vth0;
        dVaux_dVg = dVgs_eff_dVg;
        dVaux_dVd = 0.0;
        dVaux_dVb = 0.0;
      }
    }
    else if (model_.igcMod == 2)
    {
      VxNVt = (Vgs_eff - von) / T0;
      if (VxNVt > CONSTEXP_THRESHOLD)
      {
        Vaux = Vgs_eff - von;
        dVaux_dVg = dVgs_eff_dVg;
        dVaux_dVd = -dVth_dVd;
        dVaux_dVb = -dVth_dVb;
      }
    }

    if (VxNVt < -CONSTEXP_THRESHOLD)
    {
      Vaux = T0 * log(1.0 + CONSTMIN_EXP);
      dVaux_dVg = dVaux_dVd = dVaux_dVb = 0.0;
    }
    else if ((VxNVt >= -CONSTEXP_THRESHOLD) && (VxNVt <= CONSTEXP_THRESHOLD))
    {
      ExpVxNVt = exp(VxNVt);
      Vaux = T0 * log(1.0 + ExpVxNVt);
      dVaux_dVg = ExpVxNVt / (1.0 + ExpVxNVt);
      if(model_.igcMod == 1)
      {
        dVaux_dVd = 0.0;
        dVaux_dVb = 0.0;
      }
      else if (model_.igcMod == 2)
      {
        dVaux_dVd = -dVgs_eff_dVg * dVth_dVd;
        dVaux_dVb = -dVgs_eff_dVg * dVth_dVb;
      }
      dVaux_dVg *= dVgs_eff_dVg;
    }

    T2 = Vgs_eff * Vaux;
    dT2_dVg = dVgs_eff_dVg * Vaux + Vgs_eff * dVaux_dVg;
    dT2_dVd = Vgs_eff * dVaux_dVd;
    dT2_dVb = Vgs_eff * dVaux_dVb;

    T11 = paramPtr->Aechvb;
    T12 = paramPtr->Bechvb;
    T3 = paramPtr->aigc * paramPtr->cigc
       - paramPtr->bigc;
    T4 = paramPtr->bigc * paramPtr->cigc;
    T5 = T12 * (paramPtr->aigc + T3 * Voxdepinv
       - T4 * Voxdepinv * Voxdepinv);

    if (T5 > CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMAX_EXP;
      dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
    }
    else if (T5 < -CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMIN_EXP;
      dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
    }
    else
    {
      T6 = exp(T5);
      dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxdepinv);
      dT6_dVd = dT6_dVg * dVoxdepinv_dVd;
      dT6_dVb = dT6_dVg * dVoxdepinv_dVb;
      dT6_dVg *= dVoxdepinv_dVg;
    }

    Igc = T11 * T2 * T6;
    dIgc_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
    dIgc_dVd = T11 * (T2 * dT6_dVd + T6 * dT2_dVd);
    dIgc_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);

    if (model_.pigcdGiven)
    {
      Pigcd = paramPtr->pigcd;
      dPigcd_dVg = dPigcd_dVd = dPigcd_dVb = 0.0;
    }
    else
    {
      T11 = paramPtr->Bechvb * toxe;
      T12 = Vgsteff + 1.0e-20;
      T13 = T11 / T12 / T12;
      T14 = -T13 / T12;
      Pigcd = T13 * (1.0 - 0.5 * Vdseff / T12);
      dPigcd_dVg = T14 * (2.0 + 0.5 * (dVdseff_dVg
                    - 3.0 * Vdseff / T12));
      dPigcd_dVd = 0.5 * T14 * dVdseff_dVd;
      dPigcd_dVb = 0.5 * T14 * dVdseff_dVb;
    }

    T7 = -Pigcd * Vdseff; // bugfix
    dT7_dVg = -Vdseff * dPigcd_dVg - Pigcd * dVdseff_dVg;
    dT7_dVd = -Vdseff * dPigcd_dVd - Pigcd * dVdseff_dVd + dT7_dVg * dVgsteff_dVd;
    dT7_dVb = -Vdseff * dPigcd_dVb - Pigcd * dVdseff_dVb + dT7_dVg * dVgsteff_dVb;
    dT7_dVg *= dVgsteff_dVg;
    dT7_dVb *= dVbseff_dVb;
    T8 = T7 * T7 + 2.0e-4;
    dT8_dVg = 2.0 * T7;
    dT8_dVd = dT8_dVg * dT7_dVd;
    dT8_dVb = dT8_dVg * dT7_dVb;
    dT8_dVg *= dT7_dVg;

    if (T7 > CONSTEXP_THRESHOLD)
    {
      T9 = CONSTMAX_EXP;
      dT9_dVg = dT9_dVd = dT9_dVb = 0.0;
    }
    else if (T7 < -CONSTEXP_THRESHOLD)
    {
      T9 = CONSTMIN_EXP;
      dT9_dVg = dT9_dVd = dT9_dVb = 0.0;
    }
    else
    {
      T9 = exp(T7);
      dT9_dVg = T9 * dT7_dVg;
      dT9_dVd = T9 * dT7_dVd;
      dT9_dVb = T9 * dT7_dVb;
    }

    T0 = T8 * T8;
    T1 = T9 - 1.0 + 1.0e-4;
    T10 = (T1 - T7) / T8;
    dT10_dVg = (dT9_dVg - dT7_dVg - T10 * dT8_dVg) / T8;
    dT10_dVd = (dT9_dVd - dT7_dVd - T10 * dT8_dVd) / T8;
    dT10_dVb = (dT9_dVb - dT7_dVb - T10 * dT8_dVb) / T8;

    Igcs = Igc * T10;
    dIgcs_dVg = dIgc_dVg * T10 + Igc * dT10_dVg;
    dIgcs_dVd = dIgc_dVd * T10 + Igc * dT10_dVd;
    dIgcs_dVb = dIgc_dVb * T10 + Igc * dT10_dVb;

    T1 = T9 - 1.0 - 1.0e-4;
    T10 = (T7 * T9 - T1) / T8;
    dT10_dVg = (dT7_dVg * T9 + (T7 - 1.0) * dT9_dVg
             - T10 * dT8_dVg) / T8;
    dT10_dVd = (dT7_dVd * T9 + (T7 - 1.0) * dT9_dVd
             - T10 * dT8_dVd) / T8;
    dT10_dVb = (dT7_dVb * T9 + (T7 - 1.0) * dT9_dVb
             - T10 * dT8_dVb) / T8;
    Igcd = Igc * T10;
    dIgcd_dVg = dIgc_dVg * T10 + Igc * dT10_dVg;
    dIgcd_dVd = dIgc_dVd * T10 + Igc * dT10_dVd;
    dIgcd_dVb = dIgc_dVb * T10 + Igc * dT10_dVb;

    //    Igcs = Igcs;
    gIgcsg = dIgcs_dVg;
    gIgcsd = dIgcs_dVd;
    gIgcsb =  dIgcs_dVb * dVbseff_dVb;
    //    Igcd = Igcd;
    gIgcdg = dIgcd_dVg;
    gIgcdd = dIgcd_dVd;
    gIgcdb = dIgcd_dVb * dVbseff_dVb;

    T0 = vgs - (paramPtr->vfbsd + paramPtr->vfbsdoff);
    vgs_eff = sqrt(T0 * T0 + 1.0e-4);
    dvgs_eff_dvg = T0 / vgs_eff;

    T2 = vgs * vgs_eff;
    dT2_dVg = vgs * dvgs_eff_dvg + vgs_eff;
    T11 = paramPtr->AechvbEdgeS;
    T12 = paramPtr->BechvbEdge;
    T3 = paramPtr->aigs * paramPtr->cigs
       - paramPtr->bigs;
    T4 = paramPtr->bigs * paramPtr->cigs;
    T5 = T12 * (paramPtr->aigs + T3 * vgs_eff
       - T4 * vgs_eff * vgs_eff);
    if (T5 > CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMAX_EXP;
      dT6_dVg = 0.0;
    }
    else if (T5 < -CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMIN_EXP;
      dT6_dVg = 0.0;
    }
    else
    {
      T6 = exp(T5);
      dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * vgs_eff)
                * dvgs_eff_dvg;
    }
    Igs = T11 * T2 * T6;
    dIgs_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
    dIgs_dVs = -dIgs_dVg;

    T0 = vgd - (paramPtr->vfbsd + paramPtr->vfbsdoff);
    vgd_eff = sqrt(T0 * T0 + 1.0e-4);
    dvgd_eff_dvg = T0 / vgd_eff;

    T2 = vgd * vgd_eff;
    dT2_dVg = vgd * dvgd_eff_dvg + vgd_eff;
    T11 = paramPtr->AechvbEdgeD;
    T3 = paramPtr->aigd * paramPtr->cigd
      - paramPtr->bigd;
    T4 = paramPtr->bigd * paramPtr->cigd;
    T5 = T12 * (paramPtr->aigd + T3 * vgd_eff
       - T4 * vgd_eff * vgd_eff);
    if (T5 > CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMAX_EXP;
      dT6_dVg = 0.0;
    }
    else if (T5 < -CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMIN_EXP;
      dT6_dVg = 0.0;
    }
    else
    {
      T6 = exp(T5);
      dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * vgd_eff) * dvgd_eff_dvg;
    }
    Igd = T11 * T2 * T6;
    dIgd_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
    dIgd_dVd = -dIgd_dVg;

    //Igs = Igs;
    gIgsg = dIgs_dVg;
    gIgss = dIgs_dVs;
    //Igd = Igd;
    gIgdg = dIgd_dVg;
    gIgdd = dIgd_dVd;

  }
  else
  {
    Igcs = gIgcsg = gIgcsd = gIgcsb = 0.0;
    Igcd = gIgcdg = gIgcdd = gIgcdb = 0.0;
    Igs = gIgsg = gIgss = 0.0;
    Igd = gIgdg = gIgdd = 0.0;
  }

  if (model_.igbMod)
  {
    T0 = tmp * paramPtr->nigbacc;
    T1 = -Vgs_eff + Vbseff + Vfb;
    VxNVt = T1 / T0;
    if (VxNVt > CONSTEXP_THRESHOLD)
    {
      Vaux = T1;
      dVaux_dVg = -dVgs_eff_dVg;
      dVaux_dVb = 1.0;
    }
    else if (VxNVt < -CONSTEXP_THRESHOLD)
    {
      Vaux = T0 * log(1.0 + CONSTMIN_EXP);
      dVaux_dVg = dVaux_dVb = 0.0;
    }
    else
    {
      ExpVxNVt = exp(VxNVt);
      Vaux = T0 * log(1.0 + ExpVxNVt);
      dVaux_dVb = ExpVxNVt / (1.0 + ExpVxNVt);
      dVaux_dVg = -dVaux_dVb * dVgs_eff_dVg;
    }

    T2 = (Vgs_eff - Vbseff) * Vaux;
    dT2_dVg = dVgs_eff_dVg * Vaux + (Vgs_eff - Vbseff) * dVaux_dVg;
    dT2_dVb = -Vaux + (Vgs_eff - Vbseff) * dVaux_dVb;

    T11 = 4.97232e-7 * paramPtr->weff
            * paramPtr->leff * paramPtr->ToxRatio;
    T12 = -7.45669e11 * toxe;
    T3 = paramPtr->aigbacc * paramPtr->cigbacc
           - paramPtr->bigbacc;
    T4 = paramPtr->bigbacc * paramPtr->cigbacc;
    T5 = T12 * (paramPtr->aigbacc + T3 * Voxacc
         - T4 * Voxacc * Voxacc);

    if (T5 > CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMAX_EXP;
      dT6_dVg = dT6_dVb = 0.0;
    }
    else if (T5 < -CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMIN_EXP;
      dT6_dVg = dT6_dVb = 0.0;
    }
    else
    {
      T6 = exp(T5);
      dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxacc);
      dT6_dVb = dT6_dVg * dVoxacc_dVb;
      dT6_dVg *= dVoxacc_dVg;
    }

    Igbacc = T11 * T2 * T6;
    dIgbacc_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
    dIgbacc_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);


    T0 = tmp * paramPtr->nigbinv;
    T1 = Voxdepinv - paramPtr->eigbinv;
    VxNVt = T1 / T0;
    if (VxNVt > CONSTEXP_THRESHOLD)
    {
      Vaux = T1;
      dVaux_dVg = dVoxdepinv_dVg;
      dVaux_dVd = dVoxdepinv_dVd;
      dVaux_dVb = dVoxdepinv_dVb;
    }
    else if (VxNVt < -CONSTEXP_THRESHOLD)
    {
      Vaux = T0 * log(1.0 + CONSTMIN_EXP);
      dVaux_dVg = dVaux_dVd = dVaux_dVb = 0.0;
    }
    else
    {
      ExpVxNVt = exp(VxNVt);
      Vaux = T0 * log(1.0 + ExpVxNVt);
      dVaux_dVg = ExpVxNVt / (1.0 + ExpVxNVt);
      dVaux_dVd = dVaux_dVg * dVoxdepinv_dVd;
      dVaux_dVb = dVaux_dVg * dVoxdepinv_dVb;
      dVaux_dVg *= dVoxdepinv_dVg;
    }

    T2 = (Vgs_eff - Vbseff) * Vaux;
    dT2_dVg = dVgs_eff_dVg * Vaux + (Vgs_eff - Vbseff) * dVaux_dVg;
    dT2_dVd = (Vgs_eff - Vbseff) * dVaux_dVd;
    dT2_dVb = -Vaux + (Vgs_eff - Vbseff) * dVaux_dVb;

    T11 *= 0.75610;
    T12 *= 1.31724;
    T3 = paramPtr->aigbinv * paramPtr->cigbinv
       - paramPtr->bigbinv;
    T4 = paramPtr->bigbinv * paramPtr->cigbinv;
    T5 = T12 * (paramPtr->aigbinv + T3 * Voxdepinv
       - T4 * Voxdepinv * Voxdepinv);

    if (T5 > CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMAX_EXP;
      dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
    }
    else if (T5 < -CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMIN_EXP;
      dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
    }
    else
    {
      T6 = exp(T5);
      dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxdepinv);
      dT6_dVd = dT6_dVg * dVoxdepinv_dVd;
      dT6_dVb = dT6_dVg * dVoxdepinv_dVb;
      dT6_dVg *= dVoxdepinv_dVg;
    }

    Igbinv = T11 * T2 * T6;
    dIgbinv_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
    dIgbinv_dVd = T11 * (T2 * dT6_dVd + T6 * dT2_dVd);
    dIgbinv_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);

    Igb = Igbinv + Igbacc;
    gIgbg = dIgbinv_dVg + dIgbacc_dVg;
    gIgbd = dIgbinv_dVd;
    gIgbb = (dIgbinv_dVb + dIgbacc_dVb) * dVbseff_dVb;
  }
  else
  {
    Igb = gIgbg = gIgbd = gIgbs = gIgbb = 0.0;
  } // End of Gate current

  if (nf != 1.0)
  {
    cdrain *= nf;
    gds *= nf;
    gm *= nf;
    gmbs *= nf;
    IdovVds *= nf;

    gbbs *= nf;
    gbgs *= nf;
    gbds *= nf;
    csub *= nf;

    Igidl *= nf;
    ggidld *= nf;
    ggidlg *= nf;
    ggidlb *= nf;

    Igisl *= nf;
    ggisls *= nf;
    ggislg *= nf;
    ggislb *= nf;

    Igcs *= nf;
    gIgcsg *= nf;
    gIgcsd *= nf;
    gIgcsb *= nf;
    Igcd *= nf;
    gIgcdg *= nf;
    gIgcdd *= nf;
    gIgcdb *= nf;

    Igs *= nf;
    gIgsg *= nf;
    gIgss *= nf;
    Igd *= nf;
    gIgdg *= nf;
    gIgdd *= nf;

    Igb *= nf;
    gIgbg *= nf;
    gIgbd *= nf;
    gIgbb *= nf;
  }

  ggidls = -(ggidld + ggidlg + ggidlb);
  ggisld = -(ggisls + ggislg + ggislb);
  gIgbs = -(gIgbg + gIgbd + gIgbb);
  gIgcss = -(gIgcsg + gIgcsd + gIgcsb);
  gIgcds = -(gIgcdg + gIgcdd + gIgcdb);
  cd = cdrain;

  if (model_.tnoiMod == 0)
  {
    Abulk = Abulk0 * paramPtr->abulkCVfactor;
    Vdsat = Vgsteff / Abulk;
    T0 = Vdsat - Vds - CONSTDELTA_4;
    T1 = sqrt(T0 * T0 + 4.0 * CONSTDELTA_4 * Vdsat);
    if (T0 >= 0.0)
    {
      Vdseff = Vdsat - 0.5 * (T0 + T1);
    }
    else
    {
      T3 = (CONSTDELTA_4 + CONSTDELTA_4) / (T1 - T0);
      T4 = 1.0 - T3;
      T5 = Vdsat * T3 / (T1 - T0);
      Vdseff = Vdsat * T4;
    }
    if (Vds == 0.0)
    {
      Vdseff = 0.0;
    }

    T0 = Abulk * Vdseff;
    T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.0e-20);
    T2 = Vdseff / T1;
    T3 = T0 * T2;
    qinv = Coxeff * paramPtr->weffCV * nf
                    * paramPtr->leffCV
                    * (Vgsteff - 0.5 * T0 + Abulk * T3);
  }

  // C-V begins

  if ((model_.xpart < 0) || (!ChargeComputationNeeded))
  {
    qgate  = qdrn = qsrc = qbulk = 0.0;
    cggb = cgsb = cgdb = 0.0;
    cdgb = cdsb = cddb = 0.0;
    cbgb = cbsb = cbdb = 0.0;
    csgb = cssb = csdb = 0.0;
    cgbb = csbb = cdbb = cbbb = 0.0;
    cqdb = cqsb = cqgb = cqbb = 0.0;
    gtau = 0.0;
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

    CoxWL = model_.coxe * paramPtr->weffCV * paramPtr->leffCV * nf;
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
    } // Arg1 <= 0.0, end of accumulation
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
    } // Vgst <= 0.0, end of depletion
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
      {
        // 0/100 Charge partition model
        if (Vdsat <= Vds)
        {   // saturation region
          T1 = Vdsat / 3.0;
          qgate = CoxWL * (Vgs_eff - Vfb
              - paramPtr->phi - T1);
          T2 = -Two_Third_CoxWL * Vgst;
          qbulk = -(qgate + T2);
          qdrn = 0.0;

          cggb = One_Third_CoxWL * (3.0
            - dVdsat_dVg) * dVgs_eff_dVg;
          T2 = -One_Third_CoxWL * dVdsat_dVb;
          cgsb = -(cggb + T2);
          cgdb = 0.0;

          cdgb = 0.0;
          cddb = 0.0;
          cdsb = 0.0;

          cbgb = -(cggb
            - Two_Third_CoxWL * dVgs_eff_dVg);
          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
          cbsb = -(cbgb + T3);
          cbdb = 0.0;
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
          qgate = CoxWL * (Vgs_eff - Vfb
                - paramPtr->phi - 0.5 * (Vds - T3));
          T10 = T4 * T8;
          qdrn = T4 * T7;
          qbulk = -(qgate + qdrn + T10);

          T5 = T3 / T1;
          cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
                          * dVgs_eff_dVg;
          T11 = -CoxWL * T5 * dVdsat_dVb;
          cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
          cgsb = -(cggb + T11 + cgdb);
          T6 = 1.0 / Vdsat;
          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);
          T7 = T9 * T7;
          T8 = T9 * T8;
          T9 = 2.0 * T4 * (1.0 - 3.0 * T5);
          cdgb = (T7 * dAlphaz_dVg - T9 * dVdsat_dVg) * dVgs_eff_dVg;
          T12 = T7 * dAlphaz_dVb - T9 * dVdsat_dVb;
          cddb = T4 * (3.0 - 6.0 * T2 - 3.0 * T5);
          cdsb = -(cdgb + T12 + cddb);

          T9 = 2.0 * T4 * (1.0 + T5);
          T10 = (T8 * dAlphaz_dVg - T9 * dVdsat_dVg)
              * dVgs_eff_dVg;
          T11 = T8 * dAlphaz_dVb - T9 * dVdsat_dVb;
          T12 = T4 * (2.0 * T2 + T5 - 1.0);
          T0 = -(T10 + T11 + T12);

          cbgb = -(cggb + cdgb + T10);
          cbdb = -(cgdb + cddb + T12);
          cbsb = -(cgsb + cdsb + T0);
        }
      }
      else if (model_.xpart < 0.5)
      {   // 40/60 Charge partition model
        if (Vds >= Vdsat)
        {   // saturation region
          T1 = Vdsat / 3.0;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi - T1);
          T2 = -Two_Third_CoxWL * Vgst;
          qbulk = -(qgate + T2);
          qdrn = 0.4 * T2;

          cggb = One_Third_CoxWL * (3.0 - dVdsat_dVg) * dVgs_eff_dVg;
          T2 = -One_Third_CoxWL * dVdsat_dVb;
          cgsb = -(cggb + T2);
          cgdb = 0.0;

          T3 = 0.4 * Two_Third_CoxWL;
          cdgb = -T3 * dVgs_eff_dVg;
          cddb = 0.0;
          T4 = T3 * dVth_dVb; cdsb = -(T4 + cdgb);

          cbgb = -(cggb - Two_Third_CoxWL * dVgs_eff_dVg);
          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
          cbsb = -(cbgb + T3);
          cbdb = 0.0;
        }
        else
        {   // linear region
          Alphaz = Vgst / Vdsat;
          T1 = 2.0 * Vdsat - Vds;
          T2 = Vds / (3.0 * T1);
          T3 = T2 * Vds;
          T9 = 0.25 * CoxWL;
          T4 = T9 * Alphaz;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi
          - 0.5 * (Vds - T3));

          T5 = T3 / T1;
          cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
              * dVgs_eff_dVg;
          tmp = -CoxWL * T5 * dVdsat_dVb;
          cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
          cgsb = -(cggb
              + cgdb + tmp);

          T6 = 1.0 / Vdsat;
          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

          T6 = 8.0 * Vdsat * Vdsat - 6.0 * Vdsat * Vds
             + 1.2 * Vds * Vds;
          T8 = T2 / T1;
          T7 = Vds - T1 - T8 * T6;
          qdrn = T4 * T7;
          T7 *= T9;
          tmp = T8 / T1;
          tmp1 = T4 * (2.0 - 4.0 * tmp * T6
               + T8 * (16.0 * Vdsat - 6.0 * Vds));

          cdgb = (T7 * dAlphaz_dVg - tmp1
            * dVdsat_dVg) * dVgs_eff_dVg;
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
          T10 = -4.0 * T4 * (T2 - 0.5 + 0.5 * T5) - cddb;
          tmp = -(T10 + T11 + T12);

          cbgb = -(cggb + cdgb + T12);
          cbdb = -(cgdb + cddb + T10);
          cbsb = -(cgsb + cdsb + tmp);
        }
      }
      else
      {   // 50/50 partitioning
        if (Vds >= Vdsat)
        {   // saturation region
          T1 = Vdsat / 3.0;
          qgate = CoxWL * (Vgs_eff - Vfb
                - paramPtr->phi - T1);
          T2 = -Two_Third_CoxWL * Vgst;
          qbulk = -(qgate + T2);
          qdrn = 0.5 * T2;

          cggb = One_Third_CoxWL * (3.0
                          - dVdsat_dVg) * dVgs_eff_dVg;
          T2 = -One_Third_CoxWL * dVdsat_dVb;
          cgsb = -(cggb + T2);
          cgdb = 0.0;

          cdgb = -One_Third_CoxWL * dVgs_eff_dVg;
          cddb = 0.0;
          T4 = One_Third_CoxWL * dVth_dVb;
          cdsb = -(T4 + cdgb);

          cbgb = -(cggb
                            - Two_Third_CoxWL * dVgs_eff_dVg);
          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
          cbsb = -(cbgb + T3);
          cbdb = 0.0;
        }
        else
        {   // linear region
          Alphaz = Vgst / Vdsat;
          T1 = 2.0 * Vdsat - Vds;
          T2 = Vds / (3.0 * T1);
          T3 = T2 * Vds;
          T9 = 0.25 * CoxWL;
          T4 = T9 * Alphaz;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi
          - 0.5 * (Vds - T3));

          T5 = T3 / T1;
          cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
              * dVgs_eff_dVg;
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

          cdgb = (T0 * dVdsat_dVg - T7 * dAlphaz_dVg) * dVgs_eff_dVg;
          T12 = T0 * dVdsat_dVb - T7 * dAlphaz_dVb;
          cddb = T4 * (1.0 - 2.0 * T2 - T5);
          cdsb = -(cdgb + T12 + cddb);

          cbgb = -(cggb + 2.0 * cdgb);
          cbdb = -(cgdb + 2.0 * cddb);
          cbsb = -(cgsb + 2.0 * cdsb);
        } // end of linear region
      } // end of 50/50 partition
    } // end of inversion
  } // end of capMod=0
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

    CoxWL = model_.coxe * paramPtr->weffCV
          * paramPtr->leffCV * nf;

    if (model_.cvchargeMod == 0)
    {
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
      }
      // End of VgsteffCV for cvchargeMod = 0
    }
    else
    {
      T0 = n * Vtm;
      T1 = paramPtr->mstarcv * Vgst;
      T2 = T1 / T0;
      if (T2 > CONSTEXP_THRESHOLD)
      {
		      T10 = T1;
		      dT10_dVg = paramPtr->mstarcv * dVgs_eff_dVg;
		      dT10_dVd = -dVth_dVd * paramPtr->mstarcv;
		      dT10_dVb = -dVth_dVb * paramPtr->mstarcv;
      }
      else if (T2 < -CONSTEXP_THRESHOLD)
      {
        T10 = Vtm * log(1.0 + CONSTMIN_EXP);
        dT10_dVg = 0.0;
        dT10_dVd = T10 * dn_dVd;
        dT10_dVb = T10 * dn_dVb;
        T10 *= n;
      }
      else
      {
        ExpVgst = exp(T2);
        T3 = Vtm * log(1.0 + ExpVgst);
        T10 = n * T3;
        dT10_dVg = paramPtr->mstarcv * ExpVgst / (1.0 + ExpVgst);
        dT10_dVb = T3 * dn_dVb - dT10_dVg * (dVth_dVb + Vgst * dn_dVb / n);
        dT10_dVd = T3 * dn_dVd - dT10_dVg * (dVth_dVd + Vgst * dn_dVd / n);
        dT10_dVg *= dVgs_eff_dVg;
      }

      T1 = paramPtr->voffcbncv - (1.0 - paramPtr->mstarcv) * Vgst;
      T2 = T1 / T0;
      if (T2 < -CONSTEXP_THRESHOLD)
      {
        T3 = model_.coxe * CONSTMIN_EXP / paramPtr->cdep0;
        T9 = paramPtr->mstarcv + T3 * n;
        dT9_dVg = 0.0;
        dT9_dVd = dn_dVd * T3;
        dT9_dVb = dn_dVb * T3;
      }
      else if (T2 > CONSTEXP_THRESHOLD)
      {
        T3 = model_.coxe * CONSTMAX_EXP / paramPtr->cdep0;
        T9 = paramPtr->mstarcv + T3 * n;
        dT9_dVg = 0.0;
        dT9_dVd = dn_dVd * T3;
        dT9_dVb = dn_dVb * T3;
      }
      else
      {
        ExpVgst = exp(T2);
        T3 = model_.coxe / paramPtr->cdep0;
        T4 = T3 * ExpVgst;
        T5 = T1 * T4 / T0;
        T9 = paramPtr->mstarcv + n * T4;
        dT9_dVg = T3 * (paramPtr->mstarcv - 1.0) * ExpVgst / Vtm;
        dT9_dVb = T4 * dn_dVb - dT9_dVg * dVth_dVb - T5 * dn_dVb;
        dT9_dVd = T4 * dn_dVd - dT9_dVg * dVth_dVd - T5 * dn_dVd;
        dT9_dVg *= dVgs_eff_dVg;
      }

      Vgsteff = T10 / T9;
      T11 = T9 * T9;
      dVgsteff_dVg = (T9 * dT10_dVg - T10 * dT9_dVg) / T11;
      dVgsteff_dVd = (T9 * dT10_dVd - T10 * dT9_dVd) / T11;
      dVgsteff_dVb = (T9 * dT10_dVb - T10 * dT9_dVb) / T11;
      // End of VgsteffCV for cvchargeMod = 1
    }



    if (model_.capMod == 1)
    {
      Vfb = vfbzb;
      V3 = Vfb - Vgs_eff + VbseffCV - CONSTDELTA_3;
      if (Vfb <= 0.0)
      {
        T0 = sqrt(V3 * V3 - 4.0 * CONSTDELTA_3 * Vfb);
      }
      else
      {
        T0 = sqrt(V3 * V3 + 4.0 * CONSTDELTA_3 * Vfb);
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
      dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb
               + dVgsteff_dVb);

      AbulkCV = Abulk0 * paramPtr->abulkCVfactor;
      dAbulkCV_dVb = paramPtr->abulkCVfactor * dAbulk0_dVb;
      VdsatCV = Vgsteff / AbulkCV;

      T0 = VdsatCV - Vds - CONSTDELTA_4;
      dT0_dVg = 1.0 / AbulkCV;
      dT0_dVb = -VdsatCV * dAbulkCV_dVb / AbulkCV;
      T1 = sqrt(T0 * T0 + 4.0 * CONSTDELTA_4 * VdsatCV);
      dT1_dVg = (T0 + CONSTDELTA_4 + CONSTDELTA_4) / T1;
      dT1_dVd = -T0 / T1;
      dT1_dVb = dT1_dVg * dT0_dVb;
      dT1_dVg *= dT0_dVg;
      if (T0 >= 0.0)
      {
        VdseffCV = VdsatCV - 0.5 * (T0 + T1);
        dVdseffCV_dVg = 0.5 * (dT0_dVg - dT1_dVg);
        dVdseffCV_dVd = 0.5 * (1.0 - dT1_dVd);
        dVdseffCV_dVb = 0.5 * (dT0_dVb - dT1_dVb);
      }
      else
      {
        T3 = (CONSTDELTA_4 + CONSTDELTA_4) / (T1 - T0);
        T4 = 1.0 - T3;
        T5 = VdsatCV * T3 / (T1 - T0);
        VdseffCV = VdsatCV * T4;
        dVdseffCV_dVg = dT0_dVg * T4 + T5 * (dT1_dVg - dT0_dVg);
        dVdseffCV_dVd = T5 * (dT1_dVd + 1.0);
        dVdseffCV_dVb = dT0_dVb * (T4 - T5) + T5 * dT1_dVb;
      }

      if (Vds == 0.0)
      {
        VdseffCV = 0.0;
        dVdseffCV_dVg = 0.0;
        dVdseffCV_dVb = 0.0;
      }

      T0 = AbulkCV * VdseffCV;
      T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.0e-20);
      T2 = VdseffCV / T1;
      T3 = T0 * T2;

      T4 = (1.0 - 12.0 * T2 * T2 * AbulkCV);
      T5 = (6.0 * T0 * (4.0 * Vgsteff - T0) / (T1 * T1) - 0.5);
      T6 = 12.0 * T2 * T2 * Vgsteff;

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
        qsrc = -CoxWL * (0.5 * Vgsteff + 0.25 * T0
           - T0 * T0 / T1);
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
      { // 40/60 Charge petition model
        T1 = T1 / 12.0;
        T2 = 0.5 * CoxWL / (T1 * T1);
        T3 = Vgsteff * (2.0 * T0 * T0 / 3.0 + Vgsteff
           * (Vgsteff - 4.0 * T0 / 3.0))
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
      { // 50/50 Charge petition model
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
      cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb
        + Csg + Csd + Csb);
      cddb = -(Cgd + Cbd + Csd);
      cbgb = Cbg;
      cbsb = -(Cbg + Cbd + Cbb);
      cbdb = Cbd;
    }
    // Charge-Thickness capMod (CTM) begins
    else if (model_.capMod == 2)
    {
      V3 = vfbzb - Vgs_eff + VbseffCV - CONSTDELTA_3;
      if (vfbzb <= 0.0)
          T0 = sqrt(V3 * V3 - 4.0 * CONSTDELTA_3 * vfbzb);
      else
          T0 = sqrt(V3 * V3 + 4.0 * CONSTDELTA_3 * vfbzb);

      T1 = 0.5 * (1.0 + V3 / T0);
      Vfbeff = vfbzb - 0.5 * (V3 + T0);
      dVfbeff_dVg = T1 * dVgs_eff_dVg;
      dVfbeff_dVb = -T1 * dVbseffCV_dVb;

      Cox = model_.coxp;
      Tox = 1.0e8 * model_.toxp;
      T0 = (Vgs_eff - VbseffCV - vfbzb) / Tox;
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

      LINK = 1.0e-3 * model_.toxp;
      V3 = paramPtr->ldeb - Tcen - LINK;
      V4 = sqrt(V3 * V3 + 4.0 * LINK * paramPtr->ldeb);
      Tcen = paramPtr->ldeb - 0.5 * (V3 + V4);
      T1 = 0.5 * (1.0 + V3 / V4);
      dTcen_dVg *= T1;
      dTcen_dVb *= T1;

      Ccen = epssub / Tcen;
      T2 = Cox / (Cox + Ccen);
      Coxeff = T2 * Ccen;
      T3 = -Ccen / Tcen;
      dCoxeff_dVg = T2 * T2 * T3;
      dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
      dCoxeff_dVg *= dTcen_dVg;
      CoxWLcen = CoxWL * Coxeff / model_.coxe;

      Qac0 = CoxWLcen * (Vfbeff - vfbzb);
      QovCox = Qac0 / Coxeff;
      dQac0_dVg = CoxWLcen * dVfbeff_dVg
                + QovCox * dCoxeff_dVg;
      dQac0_dVb = CoxWLcen * dVfbeff_dVb
                + QovCox * dCoxeff_dVb;

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

      // Gate-bias dependent delta Phis begins
      if (paramPtr->k1ox <= 0.0)
      {   Denomi = 0.25 * paramPtr->moin * Vtm;
                      T0 = 0.5 * paramPtr->sqrtPhi;
      }
      else
      {   Denomi = paramPtr->moin * Vtm
           * paramPtr->k1ox * paramPtr->k1ox;
          T0 = paramPtr->k1ox * paramPtr->sqrtPhi;
      }
      T1 = 2.0 * T0 + Vgsteff;

      DeltaPhi = Vtm * log(1.0 + T1 * Vgsteff / Denomi);
      dDeltaPhi_dVg = 2.0 * Vtm * (T1 -T0) / (Denomi + T1 * Vgsteff);
      // End of delta Phis

      // VgDP = Vgsteff - DeltaPhi
      T0 = Vgsteff - DeltaPhi - 0.001;
      dT0_dVg = 1.0 - dDeltaPhi_dVg;
      T1 = sqrt(T0 * T0 + Vgsteff * 0.004);
      VgDP = 0.5 * (T0 + T1);
      dVgDP_dVg = 0.5 * (dT0_dVg + (T0 * dT0_dVg + 0.002) / T1);

      Tox += Tox; // WDLiu: Tcen reevaluated below due to different Vgsteff
      T0 = (Vgsteff + vtfbphi2) / Tox;
      tmp = exp(model_.bdos * 0.7 * log(T0));
      T1 = 1.0 + tmp;
      T2 = model_.bdos * 0.7 * tmp / (T0 * Tox);
      Tcen = model_.ados * 1.9e-9 / T1;
      dTcen_dVg = -Tcen * T2 / T1;
      dTcen_dVd = dTcen_dVg * dVgsteff_dVd;
      dTcen_dVb = dTcen_dVg * dVgsteff_dVb;
      dTcen_dVg *= dVgsteff_dVg;

      Ccen = epssub / Tcen;
      T0 = Cox / (Cox + Ccen);
      Coxeff = T0 * Ccen;
      T1 = -Ccen / Tcen;
      dCoxeff_dVg = T0 * T0 * T1;
      dCoxeff_dVd = dCoxeff_dVg * dTcen_dVd;
      dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
      dCoxeff_dVg *= dTcen_dVg;
      CoxWLcen = CoxWL * Coxeff / model_.coxe;

      AbulkCV = Abulk0 * paramPtr->abulkCVfactor;
      dAbulkCV_dVb = paramPtr->abulkCVfactor * dAbulk0_dVb;
      VdsatCV = VgDP / AbulkCV;

      T0 = VdsatCV - Vds - CONSTDELTA_4;
      dT0_dVg = dVgDP_dVg / AbulkCV;
      dT0_dVb = -VdsatCV * dAbulkCV_dVb / AbulkCV;
      T1 = sqrt(T0 * T0 + 4.0 * CONSTDELTA_4 * VdsatCV);
      dT1_dVg = (T0 + CONSTDELTA_4 + CONSTDELTA_4) / T1;
      dT1_dVd = -T0 / T1;
      dT1_dVb = dT1_dVg * dT0_dVb;
      dT1_dVg *= dT0_dVg;
      if (T0 >= 0.0)
      {
        VdseffCV = VdsatCV - 0.5 * (T0 + T1);
        dVdseffCV_dVg = 0.5 * (dT0_dVg - dT1_dVg);
        dVdseffCV_dVd = 0.5 * (1.0 - dT1_dVd);
        dVdseffCV_dVb = 0.5 * (dT0_dVb - dT1_dVb);
      }
      else
      {
        T3 = (CONSTDELTA_4 + CONSTDELTA_4) / (T1 - T0);
        T4 = 1.0 - T3;
        T5 = VdsatCV * T3 / (T1 - T0);
        VdseffCV = VdsatCV * T4;
        dVdseffCV_dVg = dT0_dVg * T4 + T5 * (dT1_dVg - dT0_dVg);
        dVdseffCV_dVd = T5 * (dT1_dVd + 1.0);
        dVdseffCV_dVb = dT0_dVb * (T4 - T5) + T5 * dT1_dVb;
      }

      if (Vds == 0.0)
      {
        VdseffCV = 0.0;
        dVdseffCV_dVg = 0.0;
        dVdseffCV_dVb = 0.0;
      }

      T0 = AbulkCV * VdseffCV;
      T1 = VgDP;
      T2 = 12.0 * (T1 - 0.5 * T0 + 1.0e-20);
      T3 = T0 / T2;
      T4 = 1.0 - 12.0 * T3 * T3;
      T5 = AbulkCV * (6.0 * T0 * (4.0 * T1 - T0) / (T2 * T2) - 0.5);
      T6 = T5 * VdseffCV / AbulkCV;

      qgate = CoxWLcen * (T1 - T0 * (0.5 - T3));
      QovCox = qgate / Coxeff;
      Cgg1 = CoxWLcen * (T4 * dVgDP_dVg
           + T5 * dVdseffCV_dVg);
      Cgd1 = CoxWLcen * T5 * dVdseffCV_dVd + Cgg1
           * dVgsteff_dVd + QovCox * dCoxeff_dVd;
      Cgb1 = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
           + Cgg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
      Cgg1 = Cgg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;


      T7 = 1.0 - AbulkCV;
      T8 = T2 * T2;
      T9 = 12.0 * T7 * T0 * T0 / (T8 * AbulkCV);
      T10 = T9 * dVgDP_dVg;
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
      {
        // 0/100 partition
        qsrc = -CoxWLcen * (T1 / 2.0 + T0 / 4.0
               - 0.5 * T0 * T0 / T2);
        QovCox = qsrc / Coxeff;
        T2 += T2;
        T3 = T2 * T2;
        T7 = -(0.25 - 12.0 * T0 * (4.0 * T1 - T0) / T3);
        T4 = -(0.5 + 24.0 * T0 * T0 / T3) * dVgDP_dVg;
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
      {
        // 40/60 partition
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

        Csg = T5 * dVgDP_dVg + T6 * dVdseffCV_dVg;
        Csd = Csg * dVgsteff_dVd + T6 * dVdseffCV_dVd
            + QovCox * dCoxeff_dVd;
        Csb = Csg * dVgsteff_dVb + T6 * dVdseffCV_dVb
            + T7 * dAbulkCV_dVb + QovCox * dCoxeff_dVb;
        Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
      }
      else
      {
        // 50/50 partition
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
      cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb
                        + Csg + Csd + Csb);
      cddb = -(Cgd + Cbd + Csd);
      cbgb = Cbg;
      cbsb = -(Cbg + Cbd + Cbb);
      cbdb = Cbd;
    }  // End of CTM
  }

  csgb = - cggb - cdgb - cbgb;
  csdb = - cgdb - cddb - cbdb;
  cssb = - cgsb - cdsb - cbsb;
  cgbb = - cgdb - cggb - cgsb;
  cdbb = - cddb - cdgb - cdsb;
  cbbb = - cbgb - cbdb - cbsb;
  csbb = - cgbb - cdbb - cbbb;
  // These three lines are commented out because in SPICE they're
  // here->BSIM4qgate = qgate
  // and used only for bypass.  It's pointless to do this here:
  //  qgate = qgate;
  //  qbulk = qbulk;
  //  qdrn = qdrn;
  // This line in spice is actually:
  // here->BSIM4qsrc = -(qgate + qbulk + qdrn);
  // and that saved version is *never used* --- let's not overwrite
  // our regular qsrc this way.  Don't think it matters, but it's not
  // what spice is doing here.
  //  qsrc = -(qgate + qbulk + qdrn);

  // NQS begins
  if ((trnqsMod) || (acnqsMod))
  {
    qchqs = qcheq = -(qbulk + qgate);
    cqgb = -(cggb + cbgb);
    cqdb = -(cgdb + cbdb);
    cqsb = -(cgsb + cbsb);
    cqbb = -(cqgb + cqdb + cqsb);

    CoxWL = model_.coxe * paramPtr->weffCV * nf
          * paramPtr->leffCV;
    T1 = gcrg / CoxWL; // 1 / tau
    gtau = T1 * ScalingFactor;

    if (acnqsMod)
    {
      taunet = 1.0 / T1;
    }

#if 0
    *(ckt->CKTstate0 + qcheq) = qcheq;
    if (ckt->CKTmode & MODEINITTRAN)
    {
      *(ckt->CKTstate1 + qcheq) = *(ckt->CKTstate0 + qcheq);
    }
    if (trnqsMod)
    {
      error = NIintegrate(ckt, &geq, &ceq, 0.0, qcheq);
      if (error)
      {
        return(error);
      }
    }
#endif
  }

  // Calculate junction C-V
  if (ChargeComputationNeeded)
  {
    czbd = model_.DunitAreaTempJctCap * Adeff; // bug fix
    czbs = model_.SunitAreaTempJctCap * Aseff;
    czbdsw = model_.DunitLengthSidewallTempJctCap * Pdeff;
    czbdswg = model_.DunitLengthGateSidewallTempJctCap
            * paramPtr->weffCJ * nf;
    czbssw = model_.SunitLengthSidewallTempJctCap * Pseff;
    czbsswg = model_.SunitLengthGateSidewallTempJctCap
            * paramPtr->weffCJ * nf;

    MJS = model_.SbulkJctBotGradingCoeff;
    MJSWS = model_.SbulkJctSideGradingCoeff;
    MJSWGS = model_.SbulkJctGateSideGradingCoeff;

    MJD = model_.DbulkJctBotGradingCoeff;
    MJSWD = model_.DbulkJctSideGradingCoeff;
    MJSWGD = model_.DbulkJctGateSideGradingCoeff;

    // Source Bulk Junction
    if (vbs_jct == 0.0)
    {
      qbs = 0.0;
      capbs = czbs + czbssw + czbsswg;
    }
    else if (vbs_jct < 0.0)
    {
      if (czbs > 0.0)
      {
        arg = 1.0 - vbs_jct / model_.PhiBS;
        if (MJS == 0.5)
        {
          sarg = 1.0 / sqrt(arg);
        }
        else
        {
          sarg = exp(-MJS * log(arg));
        }
        qbs = model_.PhiBS * czbs * (1.0 - arg * sarg) / (1.0 - MJS);
        capbs = czbs * sarg;
      }
      else
      {
        qbs = 0.0;
        capbs = 0.0;
      }
      if (czbssw > 0.0)
      {
        arg = 1.0 - vbs_jct / model_.PhiBSWS;
        if (MJSWS == 0.5)
        {
          sarg = 1.0 / sqrt(arg);
        }
        else
        {
          sarg = exp(-MJSWS * log(arg));
        }
        qbs += model_.PhiBSWS * czbssw
           * (1.0 - arg * sarg) / (1.0 - MJSWS);
        capbs += czbssw * sarg;
      }
      if (czbsswg > 0.0)
      {
        arg = 1.0 - vbs_jct / model_.PhiBSWGS;
        if (MJSWGS == 0.5)
        {
          sarg = 1.0 / sqrt(arg);
        }
        else
        {
          sarg = exp(-MJSWGS * log(arg));
        }
        qbs += model_.PhiBSWGS * czbsswg * (1.0 - arg * sarg) / (1.0 - MJSWGS);
        capbs += czbsswg * sarg;
      }
    }
    else
    {
      T0 = czbs + czbssw + czbsswg;
      T1 = vbs_jct * (czbs * MJS / model_.PhiBS + czbssw * MJSWS
           / model_.PhiBSWS + czbsswg * MJSWGS / model_.PhiBSWGS);

      qbs = vbs_jct * (T0 + 0.5 * T1);
      capbs = T0 + T1;
    }

    // Drain Bulk Junction
    if (vbd_jct == 0.0)
    {
      qbd = 0.0;
      capbd = czbd + czbdsw + czbdswg;
    }
    else if (vbd_jct < 0.0)
    {
      if (czbd > 0.0)
      {
        arg = 1.0 - vbd_jct / model_.PhiBD;
        if (MJD == 0.5)
        {
          sarg = 1.0 / sqrt(arg);
        }
        else
        {
          sarg = exp(-MJD * log(arg));
        }
        qbd = model_.PhiBD* czbd * (1.0 - arg * sarg) / (1.0 - MJD);
        capbd = czbd * sarg;
      }
      else
      {
        qbd = 0.0;
        capbd = 0.0;
      }
      if (czbdsw > 0.0)
      {
        arg = 1.0 - vbd_jct / model_.PhiBSWD;
        if (MJSWD == 0.5)
        {
          sarg = 1.0 / sqrt(arg);
        }
        else
        {
          sarg = exp(-MJSWD * log(arg));
        }
        qbd += model_.PhiBSWD * czbdsw
         * (1.0 - arg * sarg) / (1.0 - MJSWD);
        capbd += czbdsw * sarg;
      }
      if (czbdswg > 0.0)
      {
        arg = 1.0 - vbd_jct / model_.PhiBSWGD;
        if (MJSWGD == 0.5)
            sarg = 1.0 / sqrt(arg);
        else
            sarg = exp(-MJSWGD * log(arg));
        qbd += model_.PhiBSWGD * czbdswg
           * (1.0 - arg * sarg) / (1.0 - MJSWGD);
        capbd += czbdswg * sarg;
      }
    }
    else
    {
      T0 = czbd + czbdsw + czbdswg;
      T1 = vbd_jct * (czbd * MJD / model_.PhiBD + czbdsw * MJSWD
         / model_.PhiBSWD + czbdswg * MJSWGD / model_.PhiBSWGD);
      qbd = vbd_jct * (T0 + 0.5 * T1);
      capbd = T0 + T1;
    }
  } // ChargeComputation

  if (rgateMod == 3)
  {
    vgdx = vgmd;
    vgsx = vgms;
  }
  else  // For rgateMod == 0, 1 and 2
  {
    vgdx = vgd;
    vgsx = vgs;
  }

  // gate resistor model currents.  It is not necessary to calculate these
  // directly in spice3f5, but it is necessary in Xyce.
  Igate = IgateMid = 0.0;
  if(rgateMod == 1)
  {
    Igate = grgeltd * (Vgegp);
  }
  else if(rgateMod == 2)
  {
    Igate = (gcrg) * (Vgegp);
  }
  else if(rgateMod == 3)
  {
    Igate = grgeltd * (Vgegm);
    IgateMid = gcrg * (Vgmgp);
  }

  if (model_.capMod == 0)
  {
    cgdo = paramPtr->cgdo;
    qgdo = paramPtr->cgdo * vgdx;
    cgso = paramPtr->cgso;
    qgso = paramPtr->cgso * vgsx;
  }
  else // For both capMod == 1 and 2
  {
    T0 = vgdx + CONSTDELTA_1;
    T1 = sqrt(T0 * T0 + 4.0 * CONSTDELTA_1);
    T2 = 0.5 * (T0 - T1);

    T3 = paramPtr->weffCV * paramPtr->cgdl;
    T4 = sqrt(1.0 - 4.0 * T2 / paramPtr->ckappad);
    cgdo = paramPtr->cgdo + T3 - T3 * (1.0 - 1.0 / T4)
           * (0.5 - 0.5 * T0 / T1);
    qgdo = (paramPtr->cgdo + T3) * vgdx - T3 * (T2
           + 0.5 * paramPtr->ckappad * (T4 - 1.0));

    T0 = vgsx + CONSTDELTA_1;
    T1 = sqrt(T0 * T0 + 4.0 * CONSTDELTA_1);
    T2 = 0.5 * (T0 - T1);
    T3 = paramPtr->weffCV * paramPtr->cgsl;
    T4 = sqrt(1.0 - 4.0 * T2 / paramPtr->ckappas);
    cgso = paramPtr->cgso + T3 - T3 * (1.0 - 1.0 / T4)
         * (0.5 - 0.5 * T0 / T1);
    qgso = (paramPtr->cgso + T3) * vgsx - T3 * (T2
         + 0.5 * paramPtr->ckappas * (T4 - 1.0));
  }

  if (nf != 1.0)
  {
    cgdo *= nf;
    cgso *= nf;
    qgdo *= nf;
    qgso *= nf;
  }
  // This silliness unnecessary, conversion from spice's
  // here->BSIM4cgdo = cgdo;
  // that stuff only needed for bypass
  //  cgdo = cgdo;
  //  qgdo = qgdo;
  //  cgso = cgso;
  //  qgso = qgso;

  setupCapacitors_oldDAE();
  setupCapacitors_newDAE ();

  // Setting up a few currents for the RHS load:
  if (model_.rdsMod == 1)
  {
    Idrain = gdtot * Vddp;
    Isource = gstot * Vssp;
  }
  else
  {
    Idrain = drainConductance * Vddp;
    Isource = sourceConductance * Vssp;
  }

  // More terms that Spice leaves out because of its formulation, but which
  // Xyce absolutely needs in the RHS.
  if (model_.rbodyMod != 0)
  {
    Idbb = grbdb * Vdbb;
    Idbbp = grbpd * Vdbbp;
    Isbb = grbsb * Vsbb;
    Isbbp = grbps * Vsbbp;
    Ibpb = grbpb * Vbpb;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
//
// Purpose       : This function sets up the primaray state variables into
//                 the primary state vector.
//
//                 These variables include qbulk, qgate, qdrn and, in the
//                 event that nqsMod=1, qcdump and qcheq.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  bool bsuccess = true;

  double * staVec = extData.nextStaVectorRawPtr;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "  Begin of updatePrimaryState. \n";
    Xyce::dout() << std::endl;
  }

  bsuccess = updateIntermediateVars ();

  // voltage drops:
  double * stoVec = extData.nextStoVectorRawPtr;
  stoVec[li_store_vbd]  = vbd;
  stoVec[li_store_vbs]  = vbs;
  stoVec[li_store_vgs]  = vgs;
  stoVec[li_store_vds]  = vds;
  stoVec[li_store_vges] = vges;
  stoVec[li_store_vgms] = vgms;
  stoVec[li_store_vdes] = vdes;
  stoVec[li_store_vses] = vses;
  stoVec[li_store_vdbs] = vdbs;
  stoVec[li_store_vsbs] = vsbs;
  stoVec[li_store_vdbd] = vdbd;

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
  // Note the weirdness --- we have a "qg", "qb" and "qd" state variable,
  // but no corresponding instance variable --- they are all calculated from
  // other quantities that are NOT stored in the instance.  We use them
  // only in their derivative forms, cqg, cqb, and cqd.
  qg = staVec[li_state_qg   ] = qgate;
  qd = staVec[li_state_qd   ] = qdrn-qbd;

  if (!rbodyMod)
  {
    qb =staVec[li_state_qb   ] = qbulk+qbd+qbs;
  }
  else
  {
    qb = staVec[li_state_qb   ] = qbulk;
  }

  if (rgateMod == 3)
  {
    staVec[li_state_qgmid] = qgmid;
  }

  // parasitic capacitors:
  if (rbodyMod)
  {
    staVec[li_state_qbs] = qbs;
    staVec[li_state_qbd] = qbd;
  }

  if( trnqsMod )
  {
    staVec[li_state_qcheq] = qcheq;
    staVec[li_state_qcdump] = qdef * ScalingFactor;
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
    currStaVec[li_state_qg   ] = qgate;
    currStaVec[li_state_qd   ] = qdrn-qbd;
    if (!rbodyMod)
    {
      currStaVec[li_state_qb   ] = qbulk+qbd+qbs;
    }
    else
    {
      currStaVec[li_state_qb   ] = qbulk;
    }

    if (rgateMod == 3)
    {
      currStaVec[li_state_qgmid] = qgmid;
    }

    // parasitic capacitors:
    if (rbodyMod)
    {
      currStaVec[li_state_qbs] = qbs;
      currStaVec[li_state_qbd] = qbd;
    }

    if( trnqsMod )
    {
      currStaVec[li_state_qcheq] = qcheq;
      currStaVec[li_state_qcdump] = qdef * ScalingFactor;
    }
  }

  return bsuccess;
}


// Noise functions

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MOSFET_B4::Instance::Eval1ovFNoise
// Purpose       : 
// Special Notes :
//
// WDL: 1/f noise model has been smoothed out and enhanced with
// bulk charge effect as well as physical N* equ. and necessary
// conversion into the SI unit system.S
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/13/2018
//-----------------------------------------------------------------------------
double Instance::Eval1ovFNoise(double Vds, double freq, double temp)
{
  double cdLocal, esat, DelClm, EffFreq, N0, Nl, Leff, Leffsq;
  double T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, Ssi;

  cdLocal = fabs(cd);
  Leff = paramPtr->leff - 2.0 * model_.lintnoi;
  Leffsq = Leff * Leff;
  esat = 2.0 * vsattemp / ueff;
  if(model_.em<=0.0) 
  {
    DelClm = 0.0; /* flicker noise modified -JX  */
  }
  else 
  {
    T0 = ((((Vds - Vdseff_forNoise) / paramPtr->litl) + model_.em) / esat);
    DelClm = paramPtr->litl * log (std::max(T0, N_MINLOG));
    if (DelClm < 0.0)	
    {
      DelClm = 0.0;  /* bugfix */
    }
  }
  EffFreq = pow(freq, model_.ef);
  T1 = CONSTQ * CONSTQ * CONSTboltz * cdLocal * temp * ueff;
  T2 = 1.0e10 * EffFreq * Abulk_forNoise * model_.coxe * Leffsq;
  N0 = model_.coxe * Vgsteff_forNoise / CONSTQ;
  Nl = model_.coxe * Vgsteff_forNoise * (1.0 - AbovVgst2Vtm * Vdseff_forNoise) / CONSTQ;

  T3 = model_.oxideTrapDensityA * log(std::max(((N0 + nstar) / (Nl + nstar)), N_MINLOG));
  T4 = model_.oxideTrapDensityB * (N0 - Nl);
  T5 = model_.oxideTrapDensityC * 0.5 * (N0 * N0 - Nl * Nl);

  T6 = CONSTboltz * temp * cdLocal * cdLocal;
  T7 = 1.0e10 * EffFreq * Leffsq * paramPtr->weff * nf;
  T8 = model_.oxideTrapDensityA + model_.oxideTrapDensityB * Nl + model_.oxideTrapDensityC * Nl * Nl;
  T9 = (Nl + nstar) * (Nl + nstar);
  Ssi = T1 / T2 * (T3 + T4 + T5) + T6 / T7 * DelClm * T8 / T9;

  return Ssi;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MOSFET_B4::Instance::getNumNoiseSources
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/13/2018
//-----------------------------------------------------------------------------
int Instance::getNumNoiseSources () const
{
  return NUMNOIZ;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MOSFET_B4::Instance::setupNoiseSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 11/13/2018
//-----------------------------------------------------------------------------
void Instance::setupNoiseSources (Xyce::Analysis::NoiseData & noiseData)
{
  int numSources=NUMNOIZ;
  noiseData.numSources = numSources;
  noiseData.resize(numSources);

  noiseData.deviceName = getName().getEncodedName();

  // Note: the letter suffixes (e.g., rd) are used by the DNO() and DNI()
  // operators for .PRINT NOISE
  noiseData.noiseNames[RDNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_rd");  // noise due to rd 
  noiseData.noiseNames[RSNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_rs");  // noise due to rs 
  noiseData.noiseNames[RGNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_rg");  // noise due to rg 

  noiseData.noiseNames[RBPSNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_rbps");
  noiseData.noiseNames[RBPDNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_rbpd");
  noiseData.noiseNames[RBPBNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_rbpb");
  noiseData.noiseNames[RBSBNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_rbsb");
  noiseData.noiseNames[RBDBNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_rbdb");

  noiseData.noiseNames[IDNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_id");
  noiseData.noiseNames[FLNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_fl");
  noiseData.noiseNames[IGSNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_igs");
  noiseData.noiseNames[IGDNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_igd");
  noiseData.noiseNames[IGBNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_igb");

  // RD thermal:
  noiseData.li_Pos[RDNOIZ] = li_DrainPrime;
  noiseData.li_Neg[RDNOIZ] = li_Drain;

  // RS thermal:
  noiseData.li_Pos[RSNOIZ] = li_SourcePrime;
  noiseData.li_Neg[RSNOIZ] = li_Source;

  // RG thermal:
  //       gNodePrime, gNodeExt,
  // Xyce: li_GatePrime, 
  noiseData.li_Pos[RGNOIZ] = li_GatePrime;
  noiseData.li_Neg[RGNOIZ] = li_GateExt;

  int bodymode = 5;
  if (rbodyMod == 2)
  {	
    if( ( !model_.rbps0Given) || ( !model_.rbpd0Given) )
    {
      bodymode = 1;
    }
    else if( (!model_.rbsbx0Given && !model_.rbsby0Given) || (!model_.rbdbx0Given && !model_.rbdby0Given) )
    {
      bodymode = 3;
    }
  }

  if (rbodyMod)
  {
    if(bodymode == 5) // all 5 noise sources (RBPS, RBPD, RBPB, RBSB, RBDB)
    {
      // RBPS nodes and conductance: bNodePrime, sbNode, grbps
      //                       Xyce:  li_BodyPrime, li_SourceBody
      noiseData.li_Pos[RBPSNOIZ] = li_BodyPrime; 
      noiseData.li_Neg[RBPSNOIZ] = li_SourceBody;

      // RBPD nodes and conductance: bNodePrime, dbNode, grbpd
      //                       Xyce:  li_BodyPrime, li_DrainBody
      noiseData.li_Pos[RBPDNOIZ] = li_BodyPrime;
      noiseData.li_Neg[RBPDNOIZ] = li_DrainBody;

      // RBPB nodes and conductance: bNodePrime, bNode, grbpb
      //                       Xyce:  li_BodyPrime, li_Body
      noiseData.li_Pos[RBPBNOIZ] = li_BodyPrime;
      noiseData.li_Neg[RBPBNOIZ] = li_Body;

      // RBSB nodes and conductance: bNode, sbNode, grbsb
      //                       Xyce:  li_Body, li_SourceBody
      noiseData.li_Pos[RBSBNOIZ] = li_Body;
      noiseData.li_Neg[RBSBNOIZ] = li_SourceBody;

      // RBDB nodes and conductance: bNode, dbNode, grbdb
      //                       Xyce:  li_Body, li_DrainBody
      noiseData.li_Pos[RBDBNOIZ] = li_Body;
      noiseData.li_Neg[RBDBNOIZ] = li_DrainBody;

    }
    if(bodymode == 3) // 3 out of 5 noise sources (RBPS, RBPD, RBPB) (RBSB and RBDB excluded)
    {
      // RBPS nodes and conductance: bNodePrime, sbNode, grbps
      //                       Xyce:  li_BodyPrime, li_SourceBody
      noiseData.li_Pos[RBPSNOIZ] = li_BodyPrime; 
      noiseData.li_Neg[RBPSNOIZ] = li_SourceBody;

      // RBPD nodes and conductance: bNodePrime, dbNode, grbpd
      //                       Xyce:  li_BodyPrime, li_DrainBody
      noiseData.li_Pos[RBPDNOIZ] = li_BodyPrime;
      noiseData.li_Neg[RBPDNOIZ] = li_DrainBody;

      // RBPB nodes and conductance: bNodePrime, bNode, grbpb
      //                       Xyce:  li_BodyPrime, li_Body
      noiseData.li_Pos[RBPBNOIZ] = li_BodyPrime;
      noiseData.li_Neg[RBPBNOIZ] = li_Body;

      noiseData.li_Pos[RBSBNOIZ] = -1;
      noiseData.li_Neg[RBSBNOIZ] = -1;
      noiseData.li_Pos[RBDBNOIZ] = -1;
      noiseData.li_Neg[RBDBNOIZ] = -1;
    }
    if(bodymode == 1) // 1 out of 5 noise sources (only RBPB)
    {
      // RBPB nodes and conductance: bNodePrime, bNode, grbpb
      //                       Xyce:  li_BodyPrime, li_Body
      noiseData.li_Pos[RBPBNOIZ] = li_BodyPrime;
      noiseData.li_Neg[RBPBNOIZ] = li_Body;

      noiseData.li_Pos[RBPSNOIZ] = -1;
      noiseData.li_Neg[RBPSNOIZ] = -1;
      noiseData.li_Pos[RBPDNOIZ] = -1;
      noiseData.li_Neg[RBPDNOIZ] = -1;
      noiseData.li_Pos[RBSBNOIZ] = -1;
      noiseData.li_Neg[RBSBNOIZ] = -1;
      noiseData.li_Pos[RBDBNOIZ] = -1;
      noiseData.li_Neg[RBDBNOIZ] = -1;
    }
  }
  else
  {
    noiseData.li_Pos[RBPSNOIZ] = -1;
    noiseData.li_Neg[RBPSNOIZ] = -1;
    noiseData.li_Pos[RBPDNOIZ] = -1;
    noiseData.li_Neg[RBPDNOIZ] = -1;
    noiseData.li_Pos[RBPBNOIZ] = -1;
    noiseData.li_Neg[RBPBNOIZ] = -1;
    noiseData.li_Pos[RBSBNOIZ] = -1;
    noiseData.li_Neg[RBSBNOIZ] = -1;
    noiseData.li_Pos[RBDBNOIZ] = -1;
    noiseData.li_Neg[RBDBNOIZ] = -1;
  }

  //  dNodePrime, sNodePrime,
  //  Xyce:  li_DrainPrime, li_SourcePrime
  noiseData.li_Pos[IDNOIZ] = li_DrainPrime;
  noiseData.li_Neg[IDNOIZ] = li_SourcePrime;

  //dNodePrime, sNodePrime, 
  //  Xyce:  li_DrainPrime, li_SourcePrime
  noiseData.li_Pos[FLNOIZ] = li_DrainPrime;
  noiseData.li_Neg[FLNOIZ] = li_SourcePrime;

  //gNodePrime, sNodePrime,
  //  Xyce:  li_GatePrime, li_SourcePrime
  noiseData.li_Pos[IGSNOIZ] = li_GatePrime;
  noiseData.li_Neg[IGSNOIZ] = li_SourcePrime;

  //gNodePrime, dNodePrime,
  //  Xyce:  li_GatePrime, li_DrainPrime
  noiseData.li_Pos[IGDNOIZ] = li_GatePrime;
  noiseData.li_Neg[IGDNOIZ] = li_DrainPrime;

  //gNodePrime, bNodePrime,
  //  Xyce:  li_GatePrime, li_BodyPrime;
  noiseData.li_Pos[IGBNOIZ] = li_GatePrime;
  noiseData.li_Neg[IGBNOIZ] = li_BodyPrime;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MOSFET_B4::Instance::getNoiseSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 11/13/2018
//-----------------------------------------------------------------------------
void Instance::getNoiseSources (Xyce::Analysis::NoiseData & noiseData)
{
  double tmp = 0.0;
  double T1 = 0.0;
  double T2 = 0.0;
  double T5 = 0.0;
  double T10 = 0.0;
  double T11 = 0.0;
  double Ssi = 0.0;
  double Swi = 0.0;
  double npart_beta = 0.0;
  double npart_theta = 0.0;
  double igsquare = 0.0;

  if (model_.tnoiMod == 0)
  {
    if (model_.rdsMod == 0)
    {   
      gspr = sourceConductance;
      gdpr = drainConductance;

      if (grdsw > 0.0)
        tmp = 1.0 / grdsw; /* tmp used below */ 
      else
        tmp = 0.0;
    }
    else
    {   
      gspr = gstot;
      gdpr = gdtot;
      tmp = 0.0;
    }
  }
  else
  {
    T5 = Vgsteff_forNoise / EsatL;
    T5 *= T5;
    npart_beta = model_.rnoia * (1.0 + T5 * model_.tnoia * paramPtr->leff);
    npart_theta = model_.rnoib * (1.0 + T5 * model_.tnoib * paramPtr->leff);

    if (model_.rdsMod == 0)
    {   
      gspr = sourceConductance;
      gdpr = drainConductance;
    }
    else
    {   
      gspr = gstot;
      gdpr = gdtot;
    }

    if (vds >= 0.0)
    {
      gspr = gspr / (1.0 + npart_theta * npart_theta * gspr / IdovVds);  /* bugfix */
    }
    else
    {
      gdpr = gdpr / (1.0 + npart_theta * npart_theta * gdpr / IdovVds);
    }
  } 

  devSupport.noiseSupport(noiseData.noiseDens[RDNOIZ],
   noiseData.lnNoiseDens[RDNOIZ], THERMNOISE,
   //dNodePrime, dNode,
   gdpr,temp);

  devSupport.noiseSupport(noiseData.noiseDens[RSNOIZ],
      noiseData.lnNoiseDens[RSNOIZ],  THERMNOISE,
      //sNodePrime, sNode,
      gspr,temp);

  if ((rgateMod == 1) || (rgateMod == 2))
  {   
    devSupport.noiseSupport(noiseData.noiseDens[RGNOIZ],
      noiseData.lnNoiseDens[RGNOIZ],  THERMNOISE,
      //gNodePrime, gNodeExt,
      grgeltd,temp);
  }
  else if (rgateMod == 3)
  {   
    devSupport.noiseSupport(noiseData.noiseDens[RGNOIZ],
      noiseData.lnNoiseDens[RGNOIZ],  THERMNOISE,
      //gNodeMid, gNodeExt,
      grgeltd,temp);
  }
  else
  {    
    noiseData.noiseDens[RGNOIZ] = 0.0;
    noiseData.lnNoiseDens[RGNOIZ] = std::log(std::max(noiseData.noiseDens[RGNOIZ], N_MINLOG));
  }

  int bodymode = 5;
  if (rbodyMod == 2)
  {	
    if( ( !model_.rbps0Given) || ( !model_.rbpd0Given) )
    {
      bodymode = 1;
    }
    else if( (!model_.rbsbx0Given && !model_.rbsby0Given) || (!model_.rbdbx0Given && !model_.rbdby0Given) )
    {
      bodymode = 3;
    }
  }

  if (rbodyMod)
  { 
    if(bodymode == 5)
    {
      devSupport.noiseSupport(noiseData.noiseDens[RBPSNOIZ],
          noiseData.lnNoiseDens[RBPSNOIZ],  THERMNOISE,
          //bNodePrime, sbNode,
          grbps,temp);

      devSupport.noiseSupport(noiseData.noiseDens[RBPDNOIZ],
          noiseData.lnNoiseDens[RBPDNOIZ],  THERMNOISE,
          //bNodePrime, dbNode,
          grbpd,temp);

      devSupport.noiseSupport(noiseData.noiseDens[RBPBNOIZ],
          noiseData.lnNoiseDens[RBPBNOIZ],  THERMNOISE,
          //bNodePrime, bNode,
          grbpb,temp);

      devSupport.noiseSupport(noiseData.noiseDens[RBSBNOIZ],
          noiseData.lnNoiseDens[RBSBNOIZ],  THERMNOISE,
          //bNode, sbNode,
          grbsb,temp);

      devSupport.noiseSupport(noiseData.noiseDens[RBDBNOIZ],
          noiseData.lnNoiseDens[RBDBNOIZ],  THERMNOISE,
          //bNode, dbNode,
          grbdb,temp);
    }
    if(bodymode == 3)
    {
      devSupport.noiseSupport(noiseData.noiseDens[RBPSNOIZ],
          noiseData.lnNoiseDens[RBPSNOIZ],  THERMNOISE,
          //bNodePrime, sbNode,
          grbps,temp);

      devSupport.noiseSupport(noiseData.noiseDens[RBPDNOIZ],
          noiseData.lnNoiseDens[RBPDNOIZ],  THERMNOISE,
          //bNodePrime, dbNode,
          grbpd,temp);

      devSupport.noiseSupport(noiseData.noiseDens[RBPBNOIZ],
          noiseData.lnNoiseDens[RBPBNOIZ],  THERMNOISE,
          //bNodePrime, bNode,
          grbpb,temp);

      noiseData.noiseDens[RBSBNOIZ] = 0.0;
      noiseData.noiseDens[RBDBNOIZ] = 0.0;
      noiseData.lnNoiseDens[RBSBNOIZ] = std::log(std::max(noiseData.noiseDens[RBSBNOIZ], N_MINLOG));
      noiseData.lnNoiseDens[RBDBNOIZ] = std::log(std::max(noiseData.noiseDens[RBDBNOIZ], N_MINLOG));				     
    }
    if(bodymode == 1)
    {
      devSupport.noiseSupport(noiseData.noiseDens[RBPBNOIZ],
          noiseData.lnNoiseDens[RBPBNOIZ],  THERMNOISE,
          //bNodePrime, bNode,
          grbpb,temp);				    

      noiseData.noiseDens[RBPSNOIZ] = 0.0;
      noiseData.noiseDens[RBPDNOIZ] = 0.0;                                    

      noiseData.noiseDens[RBSBNOIZ] = 0.0;
      noiseData.noiseDens[RBDBNOIZ] = 0.0;
      noiseData.lnNoiseDens[RBPSNOIZ] = std::log(std::max(noiseData.noiseDens[RBPSNOIZ], N_MINLOG));
      noiseData.lnNoiseDens[RBPDNOIZ] = std::log(std::max(noiseData.noiseDens[RBPDNOIZ], N_MINLOG));
      noiseData.lnNoiseDens[RBSBNOIZ] = std::log(std::max(noiseData.noiseDens[RBSBNOIZ], N_MINLOG));
      noiseData.lnNoiseDens[RBDBNOIZ] = std::log(std::max(noiseData.noiseDens[RBDBNOIZ], N_MINLOG));
    }
  }
  else
  {
    noiseData.noiseDens[RBPSNOIZ] = 0.0;
    noiseData.noiseDens[RBPDNOIZ] = 0.0;   
    noiseData.noiseDens[RBPBNOIZ] = 0.0;
    noiseData.noiseDens[RBSBNOIZ] = 0.0;
    noiseData.noiseDens[RBDBNOIZ] = 0.0;
    noiseData.lnNoiseDens[RBPSNOIZ] = std::log(std::max(noiseData.noiseDens[RBPSNOIZ], N_MINLOG));
    noiseData.lnNoiseDens[RBPDNOIZ] = std::log(std::max(noiseData.noiseDens[RBPDNOIZ], N_MINLOG));
    noiseData.lnNoiseDens[RBPBNOIZ] = std::log(std::max(noiseData.noiseDens[RBPBNOIZ], N_MINLOG));
    noiseData.lnNoiseDens[RBSBNOIZ] = std::log(std::max(noiseData.noiseDens[RBSBNOIZ], N_MINLOG));
    noiseData.lnNoiseDens[RBDBNOIZ] = std::log(std::max(noiseData.noiseDens[RBDBNOIZ], N_MINLOG));
  }

  switch(model_.tnoiMod)
  {
    case 0:
      T0 = ueff * fabs(qinv);
      T1 = T0 * tmp + paramPtr->leff * paramPtr->leff;

      devSupport.noiseSupport(noiseData.noiseDens[IDNOIZ], noiseData.lnNoiseDens[IDNOIZ], THERMNOISE, 
          //dNodePrime, sNodePrime,
          (T0 / T1) * model_.ntnoi,temp);

      break;
    case 1:
      T0 = gm + gmbs + gds;
      T0 *= T0;
      igsquare = npart_theta * npart_theta * T0 / IdovVds;
      T1 = npart_beta * (gm + gmbs) + gds;
      T2 = T1 * T1 / IdovVds;

      devSupport.noiseSupport(noiseData.noiseDens[IDNOIZ], noiseData.lnNoiseDens[IDNOIZ], THERMNOISE, 
          //dNodePrime, sNodePrime, 
          (T2 - igsquare),temp);

      break;
  }

  switch(model_.fnoiMod)
  {
    case 0:
      noiseData.noiseDens[FLNOIZ] *= model_.kf * exp(model_.af * std::log(std::max(fabs(cd), N_MINLOG)))
        / (pow(noiseData.freq, model_.ef) * paramPtr->leff * paramPtr->leff * model_.coxe);

      break;
    case 1:
      double VdsLocal = vds;
      if (VdsLocal < 0.0)
      {
        VdsLocal = -VdsLocal;
      }

      Ssi = Eval1ovFNoise(VdsLocal, noiseData.freq, temp);
      T10 = model_.oxideTrapDensityA * CONSTboltz * temp;
      T11 = paramPtr->weff * nf * paramPtr->leff * pow(noiseData.freq, model_.ef) * 1.0e10 * nstar * nstar;
      Swi = T10 / T11 * cd * cd;
      T1 = Swi + Ssi;

      if (T1 > 0.0)
      {
        noiseData.noiseDens[FLNOIZ] *= (Ssi * Swi) / T1;
      }
      else
      {
        noiseData.noiseDens[FLNOIZ] *= 0.0;
      }
      break;
  }

  noiseData.lnNoiseDens[FLNOIZ] = std::log(std::max(noiseData.noiseDens[FLNOIZ], N_MINLOG));

  if(mode >= 0) 
  {  /* bugfix  */
    devSupport.noiseSupport(noiseData.noiseDens[IGSNOIZ],
        noiseData.lnNoiseDens[IGSNOIZ],  SHOTNOISE,
        //gNodePrime, sNodePrime,
        (Igs + Igcs),temp);

    devSupport.noiseSupport(noiseData.noiseDens[IGDNOIZ],
        noiseData.lnNoiseDens[IGDNOIZ],  SHOTNOISE,
        //gNodePrime, dNodePrime,
        (Igd + Igcd),temp);
  } 
  else 
  {
    devSupport.noiseSupport(noiseData.noiseDens[IGSNOIZ],
        noiseData.lnNoiseDens[IGSNOIZ],  SHOTNOISE,
        //gNodePrime, sNodePrime,
        (Igs + Igcd),temp);

    devSupport.noiseSupport(noiseData.noiseDens[IGDNOIZ],
        noiseData.lnNoiseDens[IGDNOIZ],  SHOTNOISE,
        //gNodePrime, dNodePrime,
        (Igd + Igcs),temp);
  }
  devSupport.noiseSupport(noiseData.noiseDens[IGBNOIZ],
      noiseData.lnNoiseDens[IGBNOIZ],  SHOTNOISE,
      //gNodePrime, bNodePrime,
      Igb,temp);
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 bsim4 instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  bool bsuccess=true;

  double * qVec = extData.daeQVectorRawPtr;
  double * dQdxdVp = extData.dQdxdVpVectorRawPtr;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << " MOSFET BSIM4 loadDAEQVector " << std::endl;
    Xyce::dout() << "  name             = " << getName() << std::endl;
    Xyce::dout().width(28); Xyce::dout().precision(20); Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << "  " << std::endl;
  }

  auxChargeCalculations ();

  double Qeqqg    = 0.0;   // gate charge
  double Qeqqb    = 0.0;   // bulk charge
  double Qeqqd    = 0.0;   // drain charge
  double Qeqqgmid = 0.0;   //
  double Qeqqjs   = 0.0;   // source-junction charge
  double Qeqqjd   = 0.0;   // drain-junction charge
  double Qqdef    = 0.0;   // nqs-related charge.
  double Qqcheq   = 0.0;   // nqs-related charge.

  if (model_.dtype > 0)
  {
    Qeqqg = qg;
    Qeqqd = qd;
    Qeqqb = qb;

    if (trnqsMod)
    {
     Qqdef = qdef;
     Qqcheq = qcheq;
    }

    if (rbodyMod)
    {
     Qeqqjs = qbs;
     Qeqqjd = qbd;
    }

    if (rgateMod == 3)
    {
     Qeqqgmid = qgmid;
    }
  }
  else
  {
    Qeqqg = -qg;
    Qeqqd = -qd;
    Qeqqb = -qb;

    if (trnqsMod)
    {
     Qqdef = -qdef;
     Qqcheq = -qcheq;
    }

    if (rbodyMod)
    {
     Qeqqjs = -qbs;
     Qeqqjd = -qbd;
    }

    if (rgateMod == 3)
    {
     Qeqqgmid = -qgmid;
    }
  }

  // Loading q-vector:
  qVec[li_DrainPrime] += -(-Qeqqd)*numberParallel;
  qVec[li_GatePrime] -= -(Qeqqg)*numberParallel;

  if (rgateMod == 3)
  {
    qVec[li_GateMid] -= -(+Qeqqgmid)*numberParallel;
  }

  if (!rbodyMod)
  {
    qVec[li_BodyPrime] += -(-Qeqqb)*numberParallel;
    qVec[li_SourcePrime] += -(+Qeqqg + Qeqqb + Qeqqd + Qeqqgmid)*numberParallel;
  }
  else
  {
    qVec[li_DrainBody] -= -(Qeqqjd)*numberParallel;
    qVec[li_BodyPrime] += -(-Qeqqb)*numberParallel;
    qVec[li_SourceBody] -= -(Qeqqjs)*numberParallel;
    qVec[li_SourcePrime] += -(Qeqqd + Qeqqg + Qeqqb + Qeqqjd + Qeqqjs + Qeqqgmid)*numberParallel;
  }

  if (trnqsMod)
  {
    qVec[li_Charge] += -(Qqcheq - Qqdef)*numberParallel;
  }

  // limiter section
  if (getDeviceOptions().voltageLimiterFlag && !origFlag)
  {
    dQdxdVp[li_DrainPrime] += (-Qeqqd_Jdxp)*numberParallel;
    dQdxdVp[li_GatePrime] -= (Qeqqg_Jdxp)*numberParallel;

    if (rgateMod == 3)
    {
      dQdxdVp[li_GateMid] -= (+Qeqqgmid_Jdxp)*numberParallel;
    }

    if (!rbodyMod)
    {
      dQdxdVp[li_BodyPrime] += (-Qeqqb_Jdxp)*numberParallel;
      dQdxdVp[li_SourcePrime] += (+Qeqqg_Jdxp + Qeqqb_Jdxp + Qeqqd_Jdxp + Qeqqgmid_Jdxp)*numberParallel;
    }
    else
    {
      dQdxdVp[li_DrainBody] -= (Qeqqjd_Jdxp)*numberParallel;
      dQdxdVp[li_BodyPrime] += (-Qeqqb_Jdxp)*numberParallel;
      dQdxdVp[li_SourceBody] -= (+Qeqqjs_Jdxp)*numberParallel;
      dQdxdVp[li_SourcePrime] += (+Qeqqd_Jdxp + Qeqqg_Jdxp + Qeqqb_Jdxp + Qeqqjd_Jdxp + Qeqqjs_Jdxp + Qeqqgmid_Jdxp)*numberParallel;
    }

    if (trnqsMod)
    {
      dQdxdVp[li_Charge] += (Qqcheq_Jdxp)*numberParallel;
    }
  }
    
  if( loadLeadCurrent )
  {
    double * leadQ = extData.nextLeadCurrQCompRawPtr;
    
    leadQ[li_branch_dev_id] = (Qeqqd)*numberParallel; 
    leadQ[li_branch_dev_ig] = (Qeqqg)*numberParallel; 
    leadQ[li_branch_dev_ib] = (Qeqqb)*numberParallel;
    
    if (!rbodyMod)
    {
      leadQ[li_branch_dev_is] = -(+Qeqqg + Qeqqb + Qeqqd + Qeqqgmid)*numberParallel;
    }
    else
    {
      leadQ[li_branch_dev_is] = -(Qeqqd + Qeqqg + Qeqqb + Qeqqjd + Qeqqjs + Qeqqgmid)*numberParallel;
    }
  }

  return bsuccess;
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
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
bool Instance::auxChargeCalculations ()
{
  double T0, T1;

  if (!ChargeComputationNeeded)
  {
    sxpart = (1.0 - (dxpart = (mode > 0) ? 0.4 : 0.6));
    ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
    dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
    Qeqqg_Jdxp=0.0;
    Qeqqd_Jdxp=0.0;
    Qeqqb_Jdxp=0.0;
    Qeqqgmid_Jdxp=0.0;
    Qeqqjs_Jdxp = 0.0;
    Qeqqjd_Jdxp = 0.0;

    if (trnqsMod)
    {
      CoxWL = model_.coxe * paramPtr->weffCV * nf
            * paramPtr->leffCV;
      T1 = gcrg / CoxWL;
      gtau = T1 * ScalingFactor;
    }
    else
    {
      gtau = 0.0;
    }
  }
  else  // ChargeComputation is needed
  {

    Qeqqg_Jdxp = 0.0;
    if (!origFlag)
    {
      Qeqqg_Jdxp = - CAPcggb * (vgb-vgb_orig)
                  + CAPcgdb * (vbd-vbd_orig)
                  + CAPcgsb * (vbs-vbs_orig);
    }

    Qeqqd_Jdxp = 0.0;
    if (!origFlag)
    {
      Qeqqd_Jdxp = - CAPcdgb  * (vgb-vgb_orig)
                   - CAPcdgmb * (vgmb-vgmb_orig)
                   + (CAPcddb + CAPcdbdb) * (vbd-vbd_orig)
                   - CAPcdbdb * (vbd_jct-vbd_jct_orig)
                   + CAPcdsb * (vbs-vbs_orig);
    }

    Qeqqb_Jdxp = 0.0;
    if (!origFlag)
    {
      Qeqqb_Jdxp = - CAPcbgb  * (vgb-vgb_orig)
                   - CAPcbgmb * (vgmb-vgmb_orig)
                   + CAPcbdb  * (vbd-vbd_orig)
                   + CAPcbsb  * (vbs-vbs_orig);
    }

    if (rgateMod == 3)
    {
      Qeqqgmid_Jdxp = 0.0;
      if (!origFlag)
      {
        Qeqqgmid_Jdxp = + CAPcgmdb  * (vbd-vbd_orig)
                        + CAPcgmsb  * (vbs-vbs_orig)
                        - CAPcgmgmb * (vgmb-vgmb_orig);
      }
    }
    else
    {
      Qeqqgmid_Jdxp = 0.0;
    }

    if (rbodyMod)
    {
      Qeqqjs_Jdxp = 0.0;
      Qeqqjd_Jdxp = 0.0;
      if (!origFlag)
      {
        Qeqqjs_Jdxp = CAPcsbsb * (vbs_jct-vbs_jct_orig);
        Qeqqjd_Jdxp = CAPcdbdb * (vbd_jct-vbd_jct_orig);
      }
    }

    if (trnqsMod)
    {
      // Not sure if these nqs-related terms are correct.
      T0 = ggtg * (vgb-vgb_orig) - ggtd * (vbd-vbd_orig) - ggts * (vbs-vbs_orig);

      //ceqqg += 0.0;
      Qeqqg_Jdxp += T0;
      T1 = qdef * gtau;

      //ceqqd -= 0.0;
      Qeqqd_Jdxp -= dxpart * T0
              + T1 * (ddxpart_dVg * (vgb-vgb_orig)
                    - ddxpart_dVd * (vbd-vbd_orig)
                    - ddxpart_dVs * (vbs-vbs_orig));

      cqdef = cqcdump - gqdef * qdef;

      //cqcheq = cqcheq; // redundant..
      Qqcheq_Jdxp = -(
            CAPcqgb * (vgb-vgb_orig)
          - CAPcqdb * (vbd-vbd_orig)
          - CAPcqsb * (vbs-vbs_orig)) + T0;
    }

#if 0
    if (ckt->CKTmode & MODEINITTRAN)
    {
      *(ckt->CKTstate1 + cqb) = *(ckt->CKTstate0 + cqb);
      *(ckt->CKTstate1 + cqg) = *(ckt->CKTstate0 + cqg);
      *(ckt->CKTstate1 + cqd) = *(ckt->CKTstate0 + cqd);

      if (rgateMod == 3)
      {
          *(ckt->CKTstate1 + cqgmid) = *(ckt->CKTstate0 + cqgmid);
      }

      if (rbodyMod)
      {
        *(ckt->CKTstate1 + cqbs) = *(ckt->CKTstate0 + cqbs);
        *(ckt->CKTstate1 + cqbd) = *(ckt->CKTstate0 + cqbd);
      }
    }
#endif

  } // !ChargeComputationNeeded

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupCapacitors_newDAE ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
bool Instance::setupCapacitors_newDAE ()
{
  if (mode > 0)
  {
    if (trnqsMod == 0)
    {
      if (rgateMod == 3)
      {
        CAPcgmgmb = (cgdo + cgso + paramPtr->cgbo) ;
        CAPcgmdb = -cgdo ;
        CAPcgmsb = -cgso ;
        CAPcgmbb = -paramPtr->cgbo ;

        CAPcdgmb = CAPcgmdb;
        CAPcsgmb = CAPcgmsb;
        CAPcbgmb = CAPcgmbb;

        CAPcggb = cggb ;
        CAPcgdb = cgdb ;
        CAPcgsb = cgsb ;
        CAPcgbb = -(CAPcggb + CAPcgdb + CAPcgsb);

        CAPcdgb = cdgb ;
        CAPcsgb = -(cggb + cbgb + cdgb) ;
        CAPcbgb = cbgb ;
      }
      else
      {
        CAPcggb = (cggb + cgdo + cgso + paramPtr->cgbo ) ;
        CAPcgdb = (cgdb - cgdo) ;
        CAPcgsb = (cgsb - cgso) ;
        CAPcgbb = -(CAPcggb + CAPcgdb + CAPcgsb);

        CAPcdgb = (cdgb - cgdo) ;
        CAPcsgb = -(cggb + cbgb + cdgb + cgso) ;
        CAPcbgb = (cbgb - paramPtr->cgbo) ;

        CAPcdgmb = CAPcsgmb = CAPcbgmb = 0.0;

      }
      CAPcddb = (cddb + capbd + cgdo) ;
      CAPcdsb = cdsb ;

      CAPcsdb = -(cgdb + cbdb + cddb) ;
      CAPcssb = (capbs + cgso - (cgsb + cbsb + cdsb)) ;

      if (!rbodyMod)
      {
        CAPcdbb = -(CAPcdgb + CAPcddb + CAPcdsb + CAPcdgmb);
        CAPcsbb = -(CAPcsgb + CAPcsdb + CAPcssb + CAPcsgmb);
        CAPcbdb = (cbdb - capbd) ;
        CAPcbsb = (cbsb - capbs) ;
        CAPcdbdb = 0.0; CAPcsbsb = 0.0;
      }
      else
      {
        CAPcdbb  = -(cddb + cdgb + cdsb) ;
        CAPcsbb = -(CAPcsgb + CAPcsdb + CAPcssb + CAPcsgmb) + capbs ;
        CAPcbdb = cbdb ;
        CAPcbsb = cbsb ;

        CAPcdbdb = -capbd ;
        CAPcsbsb = -capbs ;
      }
      CAPcbbb = -(CAPcbdb + CAPcbgb + CAPcbsb + CAPcbgmb);

    }
    else
    {
      CAPcqgb = cqgb ;
      CAPcqdb = cqdb ;
      CAPcqsb = cqsb ;
      CAPcqbb = cqbb ;


      if (rgateMod == 3)
      {
        CAPcgmgmb = (cgdo + cgso + paramPtr->cgbo) ;
        CAPcgmdb = -cgdo ;
        CAPcgmsb = -cgso ;
        CAPcgmbb = -paramPtr->cgbo ;

        CAPcdgmb = CAPcgmdb;
        CAPcsgmb = CAPcgmsb;
        CAPcbgmb = CAPcgmbb;

        CAPcdgb = CAPcsgb = CAPcbgb = 0.0;
        CAPcggb = CAPcgdb = CAPcgsb = CAPcgbb = 0.0;

      }
      else
      {
        CAPcggb = (cgdo + cgso + paramPtr->cgbo ) ;
        CAPcgdb = -cgdo ;
        CAPcgsb = -cgso ;
        CAPcgbb = -paramPtr->cgbo ;

        CAPcdgb = CAPcgdb;
        CAPcsgb = CAPcgsb;
        CAPcbgb = CAPcgbb;
        CAPcdgmb = CAPcsgmb = CAPcbgmb = 0.0;

      }

      CAPcddb = (capbd + cgdo) ;
      CAPcdsb = CAPcsdb = 0.0;
      CAPcssb = (capbs + cgso) ;

      if (!rbodyMod)
      {
        CAPcdbb = -(CAPcdgb + CAPcddb + CAPcdgmb);
        CAPcsbb = -(CAPcsgb + CAPcssb + CAPcsgmb);
        CAPcbdb = -capbd ;
        CAPcbsb = -capbs ;
        CAPcdbdb = 0.0; CAPcsbsb = 0.0;
      }
      else
      {
        CAPcdbb = CAPcsbb = CAPcbdb = CAPcbsb = 0.0;
        CAPcdbdb = -capbd ;
        CAPcsbsb = -capbs ;
      }
      CAPcbbb = -(CAPcbdb + CAPcbgb + CAPcbsb + CAPcbgmb);
    }
  }
  else
  {
    if (trnqsMod == 0)
    {
      if (rgateMod == 3)
      {
        CAPcgmgmb = (cgdo + cgso + paramPtr->cgbo) ;
        CAPcgmdb = -cgdo ;
        CAPcgmsb = -cgso ;
        CAPcgmbb = -paramPtr->cgbo ;

        CAPcdgmb = CAPcgmdb;
        CAPcsgmb = CAPcgmsb;
        CAPcbgmb = CAPcgmbb;

        CAPcggb = cggb ;
        CAPcgdb = cgsb ;
        CAPcgsb = cgdb ;
        CAPcgbb = -(CAPcggb + CAPcgdb + CAPcgsb);

        CAPcdgb = -(cggb + cbgb + cdgb) ;
        CAPcsgb = cdgb ;
        CAPcbgb = cbgb ;

      }
      else
      {
        CAPcggb = (cggb + cgdo + cgso + paramPtr->cgbo ) ;
        CAPcgdb = (cgsb - cgdo) ;
        CAPcgsb = (cgdb - cgso) ;
        CAPcgbb = -(CAPcggb + CAPcgdb + CAPcgsb);

        CAPcdgb = -(cggb + cbgb + cdgb + cgdo) ;
        CAPcsgb = (cdgb - cgso) ;
        CAPcbgb = (cbgb - paramPtr->cgbo) ;

        CAPcdgmb = CAPcsgmb = CAPcbgmb = 0.0;

      }
      CAPcddb = (capbd + cgdo - (cgsb + cbsb + cdsb)) ;
      CAPcdsb = -(cgdb + cbdb + cddb) ;

      CAPcsdb = cdsb ;
      CAPcssb = (cddb + capbs + cgso) ;

      if (!rbodyMod)
      {
        CAPcdbb = -(CAPcdgb + CAPcddb + CAPcdsb + CAPcdgmb);
        CAPcsbb = -(CAPcsgb + CAPcsdb + CAPcssb + CAPcsgmb);
        CAPcbdb = (cbsb - capbd) ;
        CAPcbsb = (cbdb - capbs) ;
        CAPcdbdb = 0.0; CAPcsbsb = 0.0;
      }
      else
      {
        CAPcdbb = -(CAPcdgb + CAPcddb + CAPcdsb + CAPcdgmb) + capbd ;
        CAPcsbb = -(cddb + cdgb + cdsb) ;
        CAPcbdb = cbsb ;
        CAPcbsb = cbdb ;
        CAPcdbdb = -capbd ;
        CAPcsbsb = -capbs ;
      }
      CAPcbbb = -(CAPcbgb + CAPcbdb + CAPcbsb + CAPcbgmb);
    }
    else
    {

      CAPcqgb = cqgb ;
      CAPcqdb = cqsb ;
      CAPcqsb = cqdb ;
      CAPcqbb = cqbb ;


      if (rgateMod == 3)
      {
        CAPcgmgmb = (cgdo + cgso + paramPtr->cgbo) ;
        CAPcgmdb = -cgdo ;
        CAPcgmsb = -cgso ;
        CAPcgmbb = -paramPtr->cgbo ;

        CAPcdgmb = CAPcgmdb;
        CAPcsgmb = CAPcgmsb;
        CAPcbgmb = CAPcgmbb;

        CAPcdgb = CAPcsgb = CAPcbgb = 0.0;
        CAPcggb = CAPcgdb = CAPcgsb = CAPcgbb = 0.0;

      }
      else
      {
        CAPcggb = (cgdo + cgso + paramPtr->cgbo ) ;
        CAPcgdb = -cgdo ;
        CAPcgsb = -cgso ;
        CAPcgbb = -paramPtr->cgbo ;

        CAPcdgb = CAPcgdb;
        CAPcsgb = CAPcgsb;
        CAPcbgb = CAPcgbb;
        CAPcdgmb = CAPcsgmb = CAPcbgmb = 0.0;

      }

      CAPcddb = (capbd + cgdo) ;
      CAPcdsb = CAPcsdb = 0.0;
      CAPcssb = (capbs + cgso) ;
      if (!rbodyMod)
      {
        CAPcdbb = -(CAPcdgb + CAPcddb + CAPcdgmb);
        CAPcsbb = -(CAPcsgb + CAPcssb + CAPcsgmb);
        CAPcbdb = -capbd ;
        CAPcbsb = -capbs ;
        CAPcdbdb = 0.0; CAPcsbsb = 0.0;
      }
      else
      {
        CAPcdbb = CAPcsbb = CAPcbdb = CAPcbsb = 0.0;
        CAPcdbdb = -capbd ;
        CAPcsbsb = -capbs ;
      }
      CAPcbbb = -(CAPcbdb + CAPcbgb + CAPcbsb + CAPcbgmb);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupCapacitors_oldDAE ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/17/06
//-----------------------------------------------------------------------------
bool Instance::setupCapacitors_oldDAE ()
{
  double ag0 = getSolverState().pdt_;
  double T0;

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
    if (trnqsMod == 0)
    {
      qdrn -= qgdo;
      if (rgateMod == 3)
      {
        qgmb = paramPtr->cgbo * vgmb;
        qgmid = qgdo + qgso + qgmb;
        qbulk -= qgmb;
        qsrc = -(qgate + qgmid + qbulk + qdrn);
      }
      else
      {
        qgb = paramPtr->cgbo * vgb;
        qgate += qgdo + qgso + qgb;
        qbulk -= qgb;
        qsrc = -(qgate + qbulk + qdrn);
      }

      ggtg = ggtd = ggtb = ggts = 0.0;
      sxpart = 0.6;
      dxpart = 0.4;
      ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
      dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
    }
    else
    {
      qcheq = qchqs;
      CoxWL = model_.coxe * paramPtr->weffCV * nf
        * paramPtr->leffCV;
      T0 = qdef * ScalingFactor / CoxWL;

      ggtg = gtg = T0 * gcrgg;
      ggtd = gtd = T0 * gcrgd;
      ggts = gts = T0 * gcrgs;
      ggtb = gtb = T0 * gcrgb;
      gqdef = ScalingFactor * ag0;

      if (fabs(qcheq) <= 1.0e-5 * CoxWL)
      {
        if (model_.xpart < 0.5)
        {
          dxpart = 0.4;
        }
        else if (model_.xpart > 0.5)
        {
          dxpart = 0.0;
        }
        else
        {
          dxpart = 0.5;
        }
        ddxpart_dVd = ddxpart_dVg = ddxpart_dVb
                    = ddxpart_dVs = 0.0;
      }
      else
      {
        dxpart = qdrn / qcheq;
        Cdd = cddb;
        Csd = -(cgdb + cddb
            + cbdb);
        ddxpart_dVd = (Cdd - dxpart * (Cdd + Csd)) / qcheq;
        Cdg = cdgb;
        Csg = -(cggb + cdgb
            + cbgb);
        ddxpart_dVg = (Cdg - dxpart * (Cdg + Csg)) / qcheq;

        Cds = cdsb;
        Css = -(cgsb + cdsb
            + cbsb);
        ddxpart_dVs = (Cds - dxpart * (Cds + Css)) / qcheq;

        ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);
      }
      sxpart = 1.0 - dxpart;
      dsxpart_dVd = -ddxpart_dVd;
      dsxpart_dVg = -ddxpart_dVg;
      dsxpart_dVs = -ddxpart_dVs;
      dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);

      if (rgateMod == 3)
      {
        qgmb = paramPtr->cgbo * vgmb;
        qgmid = qgdo + qgso + qgmb;
        qgate = 0.0;
        qbulk = -qgmb;
        qdrn = -qgdo;
        qsrc = -(qgmid + qbulk + qdrn);
      }
      else
      {
        qgb = paramPtr->cgbo * vgb;
        qgate = qgdo + qgso + qgb;
        qbulk = -qgb;
        qdrn = -qgdo;
        qsrc = -(qgate + qbulk + qdrn);
      }
    }
  }
  else
  {
    if (trnqsMod == 0)
    {
      qsrc = qdrn - qgso;
      if (rgateMod == 3)
      {
        qgmb = paramPtr->cgbo * vgmb;
        qgmid = qgdo + qgso + qgmb;
        qbulk -= qgmb;
        qdrn = -(qgate + qgmid + qbulk + qsrc);
      }
      else
      {
        qgb = paramPtr->cgbo * vgb;
        qgate += qgdo + qgso + qgb;
        qbulk -= qgb;
        qdrn = -(qgate + qbulk + qsrc);
      }

      ggtg = ggtd = ggtb = ggts = 0.0;
      sxpart = 0.4;
      dxpart = 0.6;
      ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
      dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
    }
    else
    {
      qcheq = qchqs;
      CoxWL = model_.coxe * paramPtr->weffCV * nf
            * paramPtr->leffCV;
      T0 = qdef * ScalingFactor / CoxWL;
      ggtg = gtg = T0 * gcrgg;
      ggts = gts = T0 * gcrgd;
      ggtd = gtd = T0 * gcrgs;
      ggtb = gtb = T0 * gcrgb;

      if (fabs(qcheq) <= 1.0e-5 * CoxWL)
      {
        if (model_.xpart < 0.5)
        {
          sxpart = 0.4;
        }
        else if (model_.xpart > 0.5)
        {
          sxpart = 0.0;
        }
        else
        {
          sxpart = 0.5;
        }
        dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
      }
      else
      {
        sxpart = qdrn / qcheq;
        Css = cddb;
        Cds = -(cgdb + cddb
            + cbdb);
        dsxpart_dVs = (Css - sxpart * (Css + Cds)) / qcheq;
        Csg = cdgb;
        Cdg = -(cggb + cdgb
            + cbgb);
        dsxpart_dVg = (Csg - sxpart * (Csg + Cdg)) / qcheq;

        Csd = cdsb;
        Cdd = -(cgsb + cdsb
            + cbsb);
        dsxpart_dVd = (Csd - sxpart * (Csd + Cdd)) / qcheq;

        dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);
      }
      dxpart = 1.0 - sxpart;
      ddxpart_dVd = -dsxpart_dVd;
      ddxpart_dVg = -dsxpart_dVg;
      ddxpart_dVs = -dsxpart_dVs;
      ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);

      if (rgateMod == 3)
      {
        qgmb = paramPtr->cgbo * vgmb;
        qgmid = qgdo + qgso + qgmb;
        qgate = 0.0;
        qbulk = -qgmb;
        qdrn = -qgdo;
        qsrc = -qgso;
      }
      else
      {
        qgb = paramPtr->cgbo * vgb;
        qgate = qgdo + qgso + qgb;
        qbulk = -qgb;
        qdrn = -qgdo;
        qsrc = -qgso;
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 bsim4 instance.
//
// Special Notes : See the special notes for loadDAEFVector.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess = true;

  double * solVec = extData.nextSolVectorRawPtr;
  double * fVec = extData.daeFVectorRawPtr;
  double * dFdxdVp = extData.dFdxdVpVectorRawPtr;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << std::endl << subsection_divider << std::endl;
    Xyce::dout() << "  Begin of Instance::loadDAEFVector. ";
    Xyce::dout() << "  origFlag = " << origFlag << "  name = " << getName() << std::endl;

    Xyce::dout().width(28); Xyce::dout().precision(20); Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << "  " << std::endl;
  }

  setupFVectorVars ();

  // Loading F-vector
  fVec[li_DrainPrime] += -(ceqjd - ceqbd - ceqdrn + Idtoteq)*numberParallel;
  fVec[li_GatePrime] -= -(-ceqgcrg + Igtoteq)*numberParallel;

  if (rgateMod == 1)
  {
    fVec[li_GateExt  ] += (Igate)*numberParallel;
    fVec[li_GatePrime] -= (Igate)*numberParallel;
  }
  else if (rgateMod == 2)
  {
    fVec[li_GateExt] += (Igate + ceqgcrg)*numberParallel;
    fVec[li_GatePrime] -= (Igate)*numberParallel;
  }
  else if (rgateMod == 3)
  {
    fVec[li_GateExt] += (Igate)*numberParallel;
    fVec[li_GateMid] += (IgateMid - Igate + ceqgcrg)*numberParallel;
    fVec[li_GatePrime] -= (IgateMid)*numberParallel;
  }

  if (!rbodyMod)
  {
    fVec[li_BodyPrime] += -(ceqbd + ceqbs - ceqjd - ceqjs + Ibtoteq)*numberParallel;
    fVec[li_SourcePrime] += -(ceqdrn - ceqbs + ceqjs + Istoteq)*numberParallel;
  }
  else
  {
    fVec[li_DrainBody] -= -(ceqjd + Idbb + Idbbp)*numberParallel;
    fVec[li_BodyPrime] += -(ceqbd + ceqbs + Ibtoteq +
                                     Idbbp + Isbbp - Ibpb)*numberParallel;
    fVec[li_Body] += - (Isbb + Idbb + Ibpb)*numberParallel;
    fVec[li_SourceBody] -= -(ceqjs + Isbb + Isbbp)*numberParallel;
    fVec[li_SourcePrime] += -(ceqdrn - ceqbs + ceqjs + Istoteq)*numberParallel;
  }

  if (model_.rdsMod)
  {
    fVec[li_Drain]  += -(-ceqgdtot)*numberParallel;
    fVec[li_Source] +=  -(ceqgstot)*numberParallel;
    fVec[li_DrainPrime]  +=  -(ceqgdtot)*numberParallel;
    fVec[li_SourcePrime] += -(-ceqgstot)*numberParallel;
  }

  // Idrain, Isource are linear terminal resistor currents
  if (drainMOSFET_B4Exists)
  {
    fVec[li_Drain]  += -(-Idrain)*numberParallel;
    fVec[li_DrainPrime]  += -(Idrain)*numberParallel;
  }

  if (sourceMOSFET_B4Exists)
  {
    fVec[li_Source] += -(-Isource)*numberParallel;
    fVec[li_SourcePrime] += -(+Isource)*numberParallel;
  }

  // Initial condition support
  if (getSolverState().dcopFlag && icVBSGiven)
  {
    double coef = extData.nextSolVectorRawPtr[li_Ibs];
    fVec[li_Body] += coef;
    fVec[li_Source] += -coef;
    double cVs = extData.nextSolVectorRawPtr[li_Source];
    double cVb = extData.nextSolVectorRawPtr[li_Body];
    fVec[li_Ibs] += (cVb-cVs-icVBS);
  }

  if (getSolverState().dcopFlag && icVDSGiven)
  {
    double coef = extData.nextSolVectorRawPtr[li_Ids];
    fVec[li_Drain] += coef;
    fVec[li_Source] += -coef;
    double cVs = extData.nextSolVectorRawPtr[li_Source];
    double cVd = extData.nextSolVectorRawPtr[li_Drain];
    fVec[li_Ids] += (cVd-cVs-icVDS);
  }

  if (getSolverState().dcopFlag && icVGSGiven)
  {
    double coef = extData.nextSolVectorRawPtr[li_Igs];
    fVec[li_GateExt] += coef;
    fVec[li_Source] += -coef;
    double cVs = extData.nextSolVectorRawPtr[li_Source];
    double cVg = extData.nextSolVectorRawPtr[li_GateExt];
    fVec[li_Igs] += (cVg-cVs-icVGS);
  }

  // limiter section
  if (getDeviceOptions().voltageLimiterFlag && !origFlag)
  {
    dFdxdVp[li_DrainPrime] += (ceqjd_Jdxp - ceqbd_Jdxp - ceqdrn_Jdxp + Idtoteq_Jdxp)*numberParallel;
    dFdxdVp[li_GatePrime] -= (- ceqgcrg_Jdxp + Igtoteq_Jdxp)*numberParallel;

    if (rgateMod == 2)
    {
      dFdxdVp[li_GateExt] += (-ceqgcrg_Jdxp)*numberParallel;
    }
    else if (rgateMod == 3)
    {
      dFdxdVp[li_GateMid] += (-ceqgcrg_Jdxp)*numberParallel;
    }

    if (!rbodyMod)
    {
      dFdxdVp[li_BodyPrime] += (ceqbd_Jdxp + ceqbs_Jdxp - ceqjd_Jdxp - ceqjs_Jdxp + Ibtoteq_Jdxp)*numberParallel;
      dFdxdVp[li_SourcePrime] += (ceqdrn_Jdxp - ceqbs_Jdxp + ceqjs_Jdxp + Istoteq_Jdxp)*numberParallel;
    }
    else
    {
      dFdxdVp[li_DrainBody] -= (ceqjd_Jdxp)*numberParallel;
      dFdxdVp[li_BodyPrime] += (ceqbd_Jdxp + ceqbs_Jdxp + Ibtoteq_Jdxp)*numberParallel;
      dFdxdVp[li_SourceBody] -= (ceqjs_Jdxp )*numberParallel;
      dFdxdVp[li_SourcePrime] += (ceqdrn_Jdxp - ceqbs_Jdxp + ceqjs_Jdxp + Istoteq_Jdxp)*numberParallel;
    }

    if (model_.rdsMod)
    {
      dFdxdVp[li_Drain] -= (ceqgdtot_Jdxp)*numberParallel;
      dFdxdVp[li_Source] += (ceqgstot_Jdxp)*numberParallel;
      dFdxdVp[li_DrainPrime] += (ceqgdtot_Jdxp)*numberParallel;
      dFdxdVp[li_SourcePrime] -= (ceqgstot_Jdxp)*numberParallel;
    }
  }
  
  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    double * junctionV = extData.nextJunctionVCompRawPtr;

    leadF[li_branch_dev_id] = -(ceqjd - ceqbd - ceqdrn + Idtoteq)*numberParallel;
    leadF[li_branch_dev_is] = Isource*numberParallel;
    leadF[li_branch_dev_ig] = (-ceqgcrg + Igtoteq)*numberParallel;
    leadF[li_branch_dev_ib] = 0;

    if (rgateMod == 1)
    {
      leadF[li_branch_dev_ig] += (Igate)*numberParallel;
    }
    else if (rgateMod == 2)
    {
      leadF[li_branch_dev_ig] += (Igate)*numberParallel;
    }
    else if (rgateMod == 3)
    {
      leadF[li_branch_dev_ig] += (IgateMid)*numberParallel;
    }

    if (!rbodyMod)
    {
      leadF[li_branch_dev_ib] += -(ceqbd + ceqbs - ceqjd - ceqjs + Ibtoteq)*numberParallel;
      leadF[li_branch_dev_is] += -(ceqdrn - ceqbs + ceqjs + Istoteq)*numberParallel;
    }
    else
    {
      leadF[li_branch_dev_ib] = - (Isbb + Idbb + Ibpb)*numberParallel;
      leadF[li_branch_dev_is] += -(ceqdrn - ceqbs + ceqjs + Istoteq)*numberParallel;
    }

    if (model_.rdsMod)
    {
      leadF[li_branch_dev_id]  += -(-ceqgdtot)*numberParallel;
      leadF[li_branch_dev_is]  +=  -(ceqgstot)*numberParallel;
    }

    junctionV[li_branch_dev_id] = solVec[li_Drain] - solVec[li_Source];
    junctionV[li_branch_dev_ig] = solVec[li_GateExt] - solVec[li_Source];
    junctionV[li_branch_dev_is] = 0.0;
    junctionV[li_branch_dev_ib] = 0.0 ;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupFVectorVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/01/08
//-----------------------------------------------------------------------------
void Instance::setupFVectorVars ()
{
  ceqdrn_Jdxp=0.0; ceqbd_Jdxp=0.0; ceqbs_Jdxp=0.0;
  Istoteq_Jdxp=0.0; Idtoteq_Jdxp=0.0;
  Ibtoteq_Jdxp=0.0; Igtoteq_Jdxp=0.0;
  ceqgcrg_Jdxp=0.0; ceqgstot_Jdxp=0.0; ceqgdtot_Jdxp=0.0;
  ceqjs_Jdxp=0.0; ceqjd_Jdxp=0.0;
  T0=0.0;

  if (mode >= 0)
  {
    Gm = gm;
    Gmbs = gmbs;
    FwdSum = Gm + Gmbs;
    RevSum = 0.0;

    ceqdrn = model_.dtype * cdrain;
    ceqdrn_Jdxp = model_.dtype *
             (-gds  * (vds-vds_orig)
              -Gm   * (vgs-vgs_orig)
              -Gmbs * (vbs-vbs_orig));

    ceqbd = model_.dtype * (csub + Igidl);
    ceqbd_Jdxp = model_.dtype * (
               - (gbds + ggidld) * (vds-vds_orig)
               - (gbgs + ggidlg) * (vgs-vgs_orig)
               - (gbbs + ggidlb) * (vbs-vbs_orig));

    ceqbs = model_.dtype * Igisl;
    ceqbs_Jdxp = model_.dtype * (
               + ggisls * (vds-vds_orig)
               - ggislg * (vgd-vgd_orig)
               - ggislb * (vbd-vbd_orig));

    gbbdp = -(gbds);
    gbbsp = gbds + gbgs + gbbs;

    gbdpg = gbgs;
    gbdpdp = gbds;
    gbdpb = gbbs;
    gbdpsp = -(gbdpg + gbdpdp + gbdpb);

    gbspg = 0.0;
    gbspdp = 0.0;
    gbspb = 0.0;
    gbspsp = 0.0;

    if (model_.igcMod)
    {
      gIstotg = gIgsg + gIgcsg;
      gIstotd = gIgcsd;
      gIstots = gIgss + gIgcss;
      gIstotb = gIgcsb;
      Istoteq = model_.dtype * (Igs + Igcs);
      Istoteq_Jdxp = model_.dtype * (
          - gIstotg * (vgs-vgs_orig)
          - gIgcsd  * (vds-vds_orig)
          - gIgcsb  * (vbs-vbs_orig));

      gIdtotg = gIgdg + gIgcdg;
      gIdtotd = gIgdd + gIgcdd;
      gIdtots = gIgcds;
      gIdtotb = gIgcdb;
      Idtoteq = model_.dtype * (Igd + Igcd);
      Idtoteq_Jdxp = model_.dtype * (
          - gIgdg  * (vgd-vgd_orig)
          - gIgcdg * (vgs-vgs_orig)
          - gIgcdd * (vds-vds_orig)
          - gIgcdb * (vbs-vbs_orig));
    }
    else
    {
      gIstotg = gIstotd = gIstots = gIstotb = Istoteq = 0.0;
      gIdtotg = gIdtotd = gIdtots = gIdtotb = Idtoteq = 0.0;
      Istoteq_Jdxp = 0.0;
      Idtoteq_Jdxp = 0.0;
    }

    if (model_.igbMod)
    {
      gIbtotg = gIgbg;
      gIbtotd = gIgbd;
      gIbtots = gIgbs;
      gIbtotb = gIgbb;
      Ibtoteq = model_.dtype * Igb;
      Ibtoteq_Jdxp = model_.dtype * (
          - gIgbg * (vgs-vgs_orig)
          - gIgbd * (vds-vds_orig)
          - gIgbb * (vbs-vbs_orig));
    }
    else
    {
      gIbtotg = gIbtotd = gIbtots = gIbtotb = Ibtoteq = 0.0;
      Ibtoteq_Jdxp = 0.0;
    }

    if ((model_.igcMod != 0) || (model_.igbMod != 0))
    {
      gIgtotg = gIstotg + gIdtotg + gIbtotg;
      gIgtotd = gIstotd + gIdtotd + gIbtotd ;
      gIgtots = gIstots + gIdtots + gIbtots;
      gIgtotb = gIstotb + gIdtotb + gIbtotb;
      Igtoteq = Istoteq + Idtoteq + Ibtoteq;
      Igtoteq_Jdxp = Istoteq_Jdxp + Idtoteq_Jdxp + Ibtoteq_Jdxp;
    }
    else
    {
      gIgtotg = gIgtotd = gIgtots = gIgtotb = Igtoteq = 0.0;
      Igtoteq_Jdxp = 0.0;
    }


    if (rgateMod == 2)
    {
      T0 = vges - vgs;
    }
    else if (rgateMod == 3)
    {
      T0 = vgms - vgs;
    }

    if (rgateMod > 1)
    {
      gcrgd = gcrgd * T0;
      gcrgg = gcrgg * T0;
      gcrgs = gcrgs * T0;
      gcrgb = gcrgb * T0;
      ceqgcrg = 0.0;
      ceqgcrg_Jdxp = -(
                        gcrgd * (vds-vds_orig)
                      + gcrgg * (vgs-vgs_orig)
                      + gcrgb * (vbs-vbs_orig));
      gcrgg -= gcrg;
      //gcrg = gcrg;
    }
    else
    {
      ceqgcrg = gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;
      ceqgcrg_Jdxp = 0.0;
    }
  }
  else
  {
    Gm = -gm;
    Gmbs = -gmbs;
    FwdSum = 0.0;
    RevSum = -(Gm + Gmbs);

    ceqdrn = -model_.dtype * cdrain;
    ceqdrn_Jdxp = -model_.dtype * (
        + gds  * (vds-vds_orig)
        + Gm   * (vgd-vgd_orig)
        + Gmbs * (vbd-vbd_orig));

    ceqbs = model_.dtype * (csub + Igisl);
    ceqbs_Jdxp = model_.dtype * (
          + (gbds + ggisls) * (vds-vds_orig)
          - (gbgs + ggislg) * (vgd-vgd_orig)
          - (gbbs + ggislb) * (vbd-vbd_orig));
    ceqbd = model_.dtype * Igidl;
    ceqbd_Jdxp = model_.dtype * (
          - ggidld * (vds-vds_orig)
          - ggidlg * (vgs-vgs_orig)
          - ggidlb * (vbs-vbs_orig));

    gbbsp = -(gbds);
    gbbdp = gbds + gbgs + gbbs;

    gbdpg = 0.0;
    gbdpsp = 0.0;
    gbdpb = 0.0;
    gbdpdp = 0.0;

    gbspg = gbgs;
    gbspsp = gbds;
    gbspb = gbbs;
    gbspdp = -(gbspg + gbspsp + gbspb);

    if (model_.igcMod)
    {
      gIstotg = gIgsg + gIgcdg;
      gIstotd = gIgcds;
      gIstots = gIgss + gIgcdd;
      gIstotb = gIgcdb;
      Istoteq = model_.dtype * (Igs + Igcd);
      Istoteq_Jdxp = model_.dtype * (
          - gIgsg  * (vgs-vgs_orig)
          - gIgcdg * (vgd-vgd_orig)
          + gIgcdd * (vds-vds_orig)
          - gIgcdb * (vbd-vbd_orig));

      gIdtotg = gIgdg + gIgcsg;
      gIdtotd = gIgdd + gIgcss;
      gIdtots = gIgcsd;
      gIdtotb = gIgcsb;
      Idtoteq = model_.dtype * (Igd + Igcs);
      Idtoteq_Jdxp = model_.dtype * (
            - (gIgdg + gIgcsg) * (vgd-vgd_orig)
            + gIgcsd * (vds-vds_orig)
            - gIgcsb * (vbd-vbd_orig));
    }
    else
    {
      gIstotg = gIstotd = gIstots = gIstotb = Istoteq = 0.0;
      gIdtotg = gIdtotd = gIdtots = gIdtotb = Idtoteq = 0.0;
      Istoteq_Jdxp = 0.0;
      Idtoteq_Jdxp = 0.0;
    }

    if (model_.igbMod)
    {
      gIbtotg = gIgbg;
      gIbtotd = gIgbs;
      gIbtots = gIgbd;
      gIbtotb = gIgbb;
      Ibtoteq = model_.dtype * Igb;
      Ibtoteq_Jdxp = model_.dtype * (
          - gIgbg * (vgd-vgd_orig)
          + gIgbd * (vds-vds_orig)
          - gIgbb * (vbd-vbd_orig));
    }
    else
    {
      gIbtotg = gIbtotd = gIbtots = gIbtotb = Ibtoteq = 0.0;
      Ibtoteq_Jdxp = 0.0;
    }

    if ((model_.igcMod != 0) || (model_.igbMod != 0))
    {
      gIgtotg = gIstotg + gIdtotg + gIbtotg;
      gIgtotd = gIstotd + gIdtotd + gIbtotd ;
      gIgtots = gIstots + gIdtots + gIbtots;
      gIgtotb = gIstotb + gIdtotb + gIbtotb;
      Igtoteq = Istoteq + Idtoteq + Ibtoteq;
      Igtoteq_Jdxp = Istoteq_Jdxp + Idtoteq_Jdxp + Ibtoteq_Jdxp;
    }
    else
    {
      gIgtotg = gIgtotd = gIgtots = gIgtotb = Igtoteq = 0.0;
      Igtoteq_Jdxp = 0.0;
    }

    if (rgateMod == 2)
    {
      T0 = vges - vgs;
    }
    else if (rgateMod == 3)
    {
      T0 = vgms - vgs;
    }

    if (rgateMod > 1)
    {
    double tmp_gcrgd = gcrgd;
      gcrgd = gcrgs * T0;
      gcrgg = gcrgg * T0;
      gcrgs = tmp_gcrgd * T0;
      gcrgb = gcrgb * T0;
      ceqgcrg = 0.0;
      ceqgcrg_Jdxp = -(
            gcrgg * (vgd-vgd_orig)
          - gcrgs * (vds-vds_orig)
          + gcrgb * (vbd-vbd_orig));
      gcrgg -= gcrg;
      //gcrg = gcrg;
    }
    else
    {
      ceqgcrg = gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;
    }
  }

  if (model_.rdsMod == 1)
  {
    ceqgstot = 0.0;
    ceqgstot_Jdxp = model_.dtype * (
          gstotd * (vds-vds_orig)
        + gstotg * (vgs-vgs_orig)
        + gstotb * (vbs-vbs_orig));
    gstots = gstots - gstot;

    ceqgdtot = 0.0;
    ceqgdtot_Jdxp = -model_.dtype * (
          gdtotd * (vds-vds_orig)
        + gdtotg * (vgs-vgs_orig)
        + gdtotb * (vbs-vbs_orig));
    gdtotd = gdtotd - gdtot;
  }
  else
  {
    gstot = gstotd = gstotg = gstots = gstotb = ceqgstot = 0.0;
    gdtot = gdtotd = gdtotg = gdtots = gdtotb = ceqgdtot = 0.0;
    ceqgstot_Jdxp = 0.0;
    ceqgdtot_Jdxp = 0.0;
  }

  if (model_.dtype > 0)
  {
    ceqjs = (cbs);
    ceqjs_Jdxp = (- gbs * (vbs_jct-vbs_jct_orig));
    ceqjd = (cbd);
    ceqjd_Jdxp = (- gbd * (vbd_jct-vbd_jct_orig));
  }
  else
  {
    ceqjs = -(cbs);
    ceqjs_Jdxp = (gbs * (vbs_jct-vbs_jct_orig));
    ceqjd = -(cbd);
    ceqjd_Jdxp = (gbd * (vbd_jct-vbd_jct_orig));
    ceqgcrg = -ceqgcrg;

    ceqgcrg_Jdxp = -ceqgcrg_Jdxp;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 bsim4 instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;

  {
    Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << std::endl << subsection_divider << std::endl;
      Xyce::dout() << "  name = " << getName() << std::endl;
    }

    if (rgateMod == 1)
    {
      dQdx[li_GatePrime][GPgp] += (CAPcggb)*numberParallel;
      dQdx[li_GatePrime][GPdp] += (CAPcgdb)*numberParallel;
      dQdx[li_GatePrime][GPsp] += (CAPcgsb)*numberParallel;
      dQdx[li_GatePrime][GPbp] += (CAPcgbb)*numberParallel;
    } // WDLiu: CAPcrg already subtracted from all CAPcrgg below
    else if (rgateMod == 2)
    {
      dQdx[li_GatePrime][GPgp] += (CAPcggb)*numberParallel;
      dQdx[li_GatePrime][GPdp] += (CAPcgdb)*numberParallel;
      dQdx[li_GatePrime][GPsp] += (CAPcgsb)*numberParallel;
      dQdx[li_GatePrime][GPbp] += (CAPcgbb)*numberParallel;
    }
    else if (rgateMod == 3)
    {
      dQdx[li_GateMid][GMgm] += (+ CAPcgmgmb)*numberParallel;

      dQdx[li_GateMid][GMdp] += (CAPcgmdb)*numberParallel;
      dQdx[li_GateMid][GMsp] += (CAPcgmsb)*numberParallel;
      dQdx[li_GateMid][GMbp] += (CAPcgmbb)*numberParallel;

      dQdx[li_DrainPrime][DPgm] += (CAPcdgmb)*numberParallel;
      dQdx[li_SourcePrime][SPgm] += (CAPcsgmb)*numberParallel;
      dQdx[li_BodyPrime][BPgm] += (CAPcbgmb)*numberParallel;

      dQdx[li_GatePrime][GPgp] += (CAPcggb)*numberParallel;
      dQdx[li_GatePrime][GPdp] += (CAPcgdb)*numberParallel;
      dQdx[li_GatePrime][GPsp] += (CAPcgsb)*numberParallel;
      dQdx[li_GatePrime][GPbp] += (CAPcgbb)*numberParallel;
    }
    else
    {
      dQdx[li_GatePrime][GPgp] += (CAPcggb)*numberParallel;
      dQdx[li_GatePrime][GPdp] += (CAPcgdb)*numberParallel;
      dQdx[li_GatePrime][GPsp] += (CAPcgsb)*numberParallel;
      dQdx[li_GatePrime][GPbp] += (CAPcgbb)*numberParallel;
    }

    dQdx[li_DrainPrime][DPdp] += (CAPcddb)*numberParallel;

    dQdx[li_DrainPrime][DPgp] += (+ CAPcdgb)*numberParallel;

    dQdx[li_DrainPrime][DPsp] -= (- CAPcdsb)*numberParallel;

    dQdx[li_DrainPrime][DPbp] -= (- CAPcdbb)*numberParallel;


    dQdx[li_SourcePrime][SPdp] -= (- CAPcsdb)*numberParallel;

    dQdx[li_SourcePrime][SPgp] += (CAPcsgb)*numberParallel;

    dQdx[li_SourcePrime][SPsp] += (CAPcssb)*numberParallel;

    dQdx[li_SourcePrime][SPbp] -= (- CAPcsbb)*numberParallel;

    dQdx[li_BodyPrime][BPdp] += (CAPcbdb)*numberParallel;
    dQdx[li_BodyPrime][BPgp] += (CAPcbgb)*numberParallel;
    dQdx[li_BodyPrime][BPsp] += (CAPcbsb)*numberParallel;
    dQdx[li_BodyPrime][BPbp] += (CAPcbbb)*numberParallel;

    if (rbodyMod)
    {
      dQdx[li_DrainPrime][DPdb] += (CAPcdbdb)*numberParallel;
      dQdx[li_SourcePrime][SPsb] -= (- CAPcsbsb)*numberParallel;

      dQdx[li_DrainBody][DBdp] += (CAPcdbdb)*numberParallel;
      dQdx[li_DrainBody][DBdb] += (- CAPcdbdb)*numberParallel;

      dQdx[li_SourceBody][SBsp] += (CAPcsbsb)*numberParallel;
      dQdx[li_SourceBody][SBsb] += (- CAPcsbsb)*numberParallel;
    }

    if (trnqsMod)
    {
      dQdx[li_Charge][Qgp] += (- CAPcqgb)*numberParallel;
      dQdx[li_Charge][Qdp] += (- CAPcqdb)*numberParallel;
      dQdx[li_Charge][Qsp] += (- CAPcqsb)*numberParallel;
      dQdx[li_Charge][Qbp] += (- CAPcqbb)*numberParallel;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 bsim4 instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << std::endl << subsection_divider << std::endl;
    Xyce::dout() << "Instance::loadDAEdFdx";
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  if (!rbodyMod)
  {
    gjbd = gbd;
    gjbs = gbs;
  }
  else
  {
    gjbd = gjbs = 0.0;
  }

  if (!model_.rdsMod)
  {
    gdpr = drainConductance;
    gspr = sourceConductance;
  }
  else
  {
    gdpr = gspr = 0.0;
  }

  geltd = grgeltd;

  double T1 = qdef * gtau;

  if (rgateMod == 1)
  {
    dFdx[li_GateExt][GEge] += (geltd)*numberParallel;
    dFdx[li_GateExt][GEgp] -= (geltd)*numberParallel;
    dFdx[li_GatePrime][GPge] -= (geltd)*numberParallel;
    dFdx[li_GatePrime][GPgp] += (+ geltd - ggtg + gIgtotg)*numberParallel;
    dFdx[li_GatePrime][GPdp] += (- ggtd + gIgtotd)*numberParallel;
    dFdx[li_GatePrime][GPsp] += (- ggts + gIgtots)*numberParallel;
    dFdx[li_GatePrime][GPbp] += (- ggtb + gIgtotb)*numberParallel;
  } // WDLiu: gcrg already subtracted from all gcrgg below
  else if (rgateMod == 2)
  {
    dFdx[li_GateExt][GEge] += (gcrg)*numberParallel;
    dFdx[li_GateExt][GEgp] += (gcrgg)*numberParallel;
    dFdx[li_GateExt][GEdp] += (gcrgd)*numberParallel;
    dFdx[li_GateExt][GEsp] += (gcrgs)*numberParallel;
    dFdx[li_GateExt][GEbp] += (gcrgb)*numberParallel;

    dFdx[li_GatePrime][GPge] -= (gcrg)*numberParallel;
    dFdx[li_GatePrime][GPgp] += (- gcrgg - ggtg + gIgtotg)*numberParallel;
    dFdx[li_GatePrime][GPdp] += (- gcrgd - ggtd + gIgtotd)*numberParallel;
    dFdx[li_GatePrime][GPsp] += (- gcrgs - ggts + gIgtots)*numberParallel;
    dFdx[li_GatePrime][GPbp] += (- gcrgb - ggtb + gIgtotb)*numberParallel;
  }
  else if (rgateMod == 3)
  {
    dFdx[li_GateExt][GEge] += (geltd)*numberParallel;
    dFdx[li_GateExt][GEgm] -= (geltd)*numberParallel;
    dFdx[li_GateMid][GMge] -= (geltd)*numberParallel;
    dFdx[li_GateMid][GMgm] += (geltd + gcrg)*numberParallel;

    dFdx[li_GateMid][GMdp] += (gcrgd)*numberParallel;
    dFdx[li_GateMid][GMgp] += (gcrgg)*numberParallel;
    dFdx[li_GateMid][GMsp] += (gcrgs)*numberParallel;
    dFdx[li_GateMid][GMbp] += (gcrgb)*numberParallel;

    dFdx[li_GatePrime][GPgm] -= (gcrg)*numberParallel;

    dFdx[li_GatePrime][GPgp] += (- gcrgg - ggtg + gIgtotg)*numberParallel;
    dFdx[li_GatePrime][GPdp] += (- gcrgd - ggtd + gIgtotd)*numberParallel;
    dFdx[li_GatePrime][GPsp] += (- gcrgs - ggts + gIgtots)*numberParallel;
    dFdx[li_GatePrime][GPbp] += (- gcrgb - ggtb + gIgtotb)*numberParallel;
  }
  else
  {
    dFdx[li_GatePrime][GPgp] += (- ggtg + gIgtotg)*numberParallel;
    dFdx[li_GatePrime][GPdp] += (- ggtd + gIgtotd)*numberParallel;
    dFdx[li_GatePrime][GPsp] += (- ggts + gIgtots)*numberParallel;
    dFdx[li_GatePrime][GPbp] += (- ggtb + gIgtotb)*numberParallel;
  }

  if (model_.rdsMod)
  {
    dFdx[li_Drain][Dgp] += (gdtotg)*numberParallel;
    dFdx[li_Drain][Dsp] += (gdtots)*numberParallel;
    dFdx[li_Drain][Dbp] += (gdtotb)*numberParallel;
    dFdx[li_Source][Sdp] += (gstotd)*numberParallel;
    dFdx[li_Source][Sgp] += (gstotg)*numberParallel;
    dFdx[li_Source][Sbp] += (gstotb)*numberParallel;
  }

  dFdx[li_DrainPrime][DPdp] += (gdpr + gds + gbd + T1 * ddxpart_dVd
             - gdtotd + RevSum + gbdpdp + dxpart * ggtd - gIdtotd)*numberParallel;

  dFdx[li_DrainPrime][DPd] -= (gdpr + gdtot)*numberParallel;
  dFdx[li_DrainPrime][DPgp] += (Gm - gdtotg + gbdpg - gIdtotg
             + dxpart * ggtg + T1 * ddxpart_dVg)*numberParallel;

  dFdx[li_DrainPrime][DPsp] -= (gds + gdtots - dxpart * ggts + gIdtots
             - T1 * ddxpart_dVs + FwdSum - gbdpsp)*numberParallel;

  dFdx[li_DrainPrime][DPbp] -= (gjbd + gdtotb - Gmbs - gbdpb + gIdtotb
             - T1 * ddxpart_dVb - dxpart * ggtb)*numberParallel;

  dFdx[li_Drain][Ddp] -= (gdpr - gdtotd)*numberParallel;
  dFdx[li_Drain][Dd] += (gdpr + gdtot)*numberParallel;

  dFdx[li_SourcePrime][SPdp] -= (gds + gstotd + RevSum - gbspdp
                         - T1 * dsxpart_dVd - sxpart * ggtd + gIstotd)*numberParallel;

  dFdx[li_SourcePrime][SPgp] += (- Gm - gstotg + gbspg + sxpart * ggtg
                       + T1 * dsxpart_dVg - gIstotg)*numberParallel;

  dFdx[li_SourcePrime][SPsp] += (gspr + gds + gbs + T1 * dsxpart_dVs
                         - gstots + FwdSum + gbspsp + sxpart * ggts - gIstots)*numberParallel;

  dFdx[li_SourcePrime][SPs] -= (gspr + gstot)*numberParallel;

  dFdx[li_SourcePrime][SPbp] -= (gjbs + gstotb + Gmbs - gbspb - sxpart * ggtb
                         - T1 * dsxpart_dVb + gIstotb)*numberParallel;

  dFdx[li_Source][Ssp] -= (gspr - gstots)*numberParallel;
  dFdx[li_Source][Ss] += (gspr + gstot)*numberParallel;

  dFdx[li_BodyPrime][BPdp] += (- gjbd + gbbdp - gIbtotd)*numberParallel;
  dFdx[li_BodyPrime][BPgp] += (- gbgs - gIbtotg)*numberParallel;
  dFdx[li_BodyPrime][BPsp] += (- gjbs + gbbsp - gIbtots)*numberParallel;
  dFdx[li_BodyPrime][BPbp] += (gjbd + gjbs - gbbs - gIbtotb)*numberParallel;

  //ggidld = (ggidld)*numberParallel;
  //ggidlg = (ggidlg)*numberParallel;
  //ggidlb = (ggidlb)*numberParallel;
  //ggislg = (ggislg)*numberParallel;
  //ggisls = (ggisls)*numberParallel;
  //ggislb = (ggislb)*numberParallel;

  // stamp gidl
  dFdx[li_DrainPrime][DPdp] += (ggidld)*numberParallel;
  dFdx[li_DrainPrime][DPgp] += (ggidlg)*numberParallel;
  dFdx[li_DrainPrime][DPsp] -= ((ggidlg + ggidld + ggidlb))*numberParallel;
  dFdx[li_DrainPrime][DPbp] += (ggidlb)*numberParallel;
  dFdx[li_BodyPrime][BPdp] -= (ggidld)*numberParallel;
  dFdx[li_BodyPrime][BPgp] -= (ggidlg)*numberParallel;
  dFdx[li_BodyPrime][BPsp] += ((ggidlg + ggidld + ggidlb))*numberParallel;
  dFdx[li_BodyPrime][BPbp] -= (ggidlb)*numberParallel;
  // stamp gisl
  dFdx[li_SourcePrime][SPdp] -= ((ggisls + ggislg + ggislb))*numberParallel;
  dFdx[li_SourcePrime][SPgp] += (ggislg)*numberParallel;
  dFdx[li_SourcePrime][SPsp] += (ggisls)*numberParallel;
  dFdx[li_SourcePrime][SPbp] += (ggislb)*numberParallel;
  dFdx[li_BodyPrime][BPdp] += ((ggislg + ggisls + ggislb))*numberParallel;
  dFdx[li_BodyPrime][BPgp] -= (ggislg)*numberParallel;
  dFdx[li_BodyPrime][BPsp] -= (ggisls)*numberParallel;
  dFdx[li_BodyPrime][BPbp] -= (ggislb)*numberParallel;


  if (rbodyMod)
  {
    dFdx[li_DrainPrime][DPdb] += (- gbd)*numberParallel;
    dFdx[li_SourcePrime][SPsb] -= (gbs)*numberParallel;

    dFdx[li_DrainBody][DBdp] += (- gbd)*numberParallel;
    dFdx[li_DrainBody][DBdb] += (gbd + grbpd + grbdb)*numberParallel;
    dFdx[li_DrainBody][DBbp] -= (grbpd)*numberParallel;
    dFdx[li_DrainBody][DBb] -= (grbdb)*numberParallel;

    dFdx[li_BodyPrime][BPdb] -= (grbpd)*numberParallel;
    dFdx[li_BodyPrime][BPb] -= (grbpb)*numberParallel;
    dFdx[li_BodyPrime][BPsb] -= (grbps)*numberParallel;
    dFdx[li_BodyPrime][BPbp] += (grbpd + grbps + grbpb)*numberParallel;
    // WDLiu: (gcbbb - gbbs) already added to BPbpPtr

    dFdx[li_SourceBody][SBsp] += (- gbs)*numberParallel;
    dFdx[li_SourceBody][SBbp] -= (grbps)*numberParallel;
    dFdx[li_SourceBody][SBb] -= (grbsb)*numberParallel;
    dFdx[li_SourceBody][SBsb] += (gbs + grbps + grbsb)*numberParallel;

    dFdx[li_Body][Bdb] -= (grbdb)*numberParallel;
    dFdx[li_Body][Bbp] -= (grbpb)*numberParallel;
    dFdx[li_Body][Bsb] -= (grbsb)*numberParallel;
    dFdx[li_Body][Bb] += (grbsb + grbdb + grbpb)*numberParallel;
  }

  if (trnqsMod)
  {
    dFdx[li_Charge][Qq] += (gqdef + gtau)*numberParallel;
    dFdx[li_Charge][Qgp] += (ggtg)*numberParallel;
    dFdx[li_Charge][Qdp] += (ggtd)*numberParallel;
    dFdx[li_Charge][Qsp] += (ggts)*numberParallel;
    dFdx[li_Charge][Qbp] += (ggtb)*numberParallel;

    dFdx[li_DrainPrime][DPq] += (dxpart * gtau)*numberParallel;
    dFdx[li_SourcePrime][SPq] += (sxpart * gtau)*numberParallel;
    dFdx[li_GatePrime][GPq] -= (gtau)*numberParallel;
  }


  // Initial Conditions:
  if (icVBSGiven)
  {
    if (getSolverState().dcopFlag)
    {
      dFdx[li_Body][Bibs] += 1.0;
      dFdx[li_Source][Sibs] += -1.0;
      dFdx[li_Ibs][IBSb] += 1.0;
      dFdx[li_Ibs][IBSs] += -1.0;
    }
    else
    {
      dFdx[li_Ibs][IBSibs] = 1.0;
    }
  }

  if (icVDSGiven)
  {
    if (getSolverState().dcopFlag)
    {
      dFdx[li_Drain][Dids] += 1.0;
      dFdx[li_Source][Sids] += -1.0;
      dFdx[li_Ids][IDSd] += 1.0;
      dFdx[li_Ids][IDSs] += -1.0;
    }
    else
    {
      dFdx[li_Ids][IDSids] = 1.0;
    }
  }

  if (icVGSGiven)
  {
    if (getSolverState().dcopFlag)
    {dFdx[li_GateExt][GEigs] += 1.0;
      dFdx[li_Source][Sigs] += -1.0;
      dFdx[li_Igs][IGSg] += 1.0;
      dFdx[li_Igs][IGSs] += -1.0;
    }
    else
    {
      dFdx[li_Igs][IGSigs] = 1.0;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
bool Instance::setIC ()
{
  bool bsuccess = true;


  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance:polyDepletion
// Purpose       : Function to compute poly depletion effect .
//
// Special Notes : Vgs_arg is named to distinguish it from the instance
//                 variable, Vgs.
//
// Scope         : public
// Creator       : Eric Keiter,SNL
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
int Instance::polyDepletion( double  phi, double  ngate,
                                            double epsgate,
                                            double  coxe, double  Vgs_arg,
                                            double & Vgs_eff,
                                            double & dVgs_eff_dVg)
{
  double T1(0.0), T2(0.0), T3(0.0), T4(0.0), T5(0.0), T6(0.0), T7(0.0), T8(0.0);

  // CONSTQ = 1.6021918e-19
  // CONSTEPSSI = 1.03594e-10
  // Poly Gate Si Depletion Effect
  if ((ngate > 1.0e18) && (ngate < 1.0e25) && (Vgs_arg > phi) && (epsgate!=0))
  {
    T1 = 1.0e6 * CONSTQ * epsgate * ngate / (coxe * coxe);
    T8 = Vgs_arg - phi;
    T4 = sqrt(1.0 + 2.0 * T8 / T1);
    T2 = 2.0 * T8 / (T4 + 1.0);
    T3 = 0.5 * T2 * T2 / T1; // T3 = Vpoly
    T7 = 1.12 - T3 - 0.05;
    T6 = sqrt(T7 * T7 + 0.224);
    T5 = 1.12 - 0.5 * (T7 + T6);
    Vgs_eff = Vgs_arg - T5;
    dVgs_eff_dVg = 1.0 - (0.5 - 0.5 / T4) * (1.0 + T7 / T6);
  }
  else
  {
    Vgs_eff = Vgs_arg;
    dVgs_eff_dVg = 1.0;
  }
  return(0);
}

//-----------------------------------------------------------------------------
// Function      : Instance:DioIjthVjmEval
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter,SNL
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
int Instance::DioIjthVjmEval
  (double Nvtm, double Ijth, double Isb, double XExpBV, double & Vjm)
{
  double Tb(0.0), Tc(0.0), EVjmovNv(0.0);

  Tc = XExpBV;
  Tb = 1.0 + Ijth / Isb - Tc;
  EVjmovNv = 0.5 * (Tb + sqrt(Tb * Tb + 4.0 * Tc));
  Vjm = Nvtm * log(EVjmovNv);

  return 0;
}

//
// WDLiu:
// This subroutine is a special module to process the geometry dependent
// parasitics for BSIM4, which calculates Ps, Pd, As, Ad, and Rs and  Rd
// for multi-fingers and varous GEO and RGEO options.
//
//-----------------------------------------------------------------------------
// Function      : Instance::NumFingerDiff
// Purpose       :
//
// Special Notes : nf_arg is named to distinguish it from the
//                 instance variable nf.
//
// Scope         : public
// Creator       : Eric Keiter,SNL
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
int Instance::NumFingerDiff
   (double nf_arg, int minSD,
    double & nuIntD, double & nuEndD, double & nuIntS, double & nuEndS)
{
  int NF = static_cast<int>(nf_arg);

	if ((NF%2) != 0)
	{
    nuEndD = nuEndS = 1.0;
    nuIntD = nuIntS = 2.0 * std::max((nf_arg - 1.0) / 2.0, 0.0);
	}
	else
  {
    if (minSD == 1) // minimize # of source
    {
      nuEndD = 2.0;
      nuIntD = 2.0 * std::max((nf_arg / 2.0 - 1.0), 0.0);
      nuEndS = 0.0;
      nuIntS = nf_arg;
    }
    else
    {
      nuEndD = 0.0;
      nuIntD = nf_arg;
      nuEndS = 2.0;
      nuIntS = 2.0 * std::max((nf_arg / 2.0 - 1.0), 0.0);
    }
	}
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : Instance::PAeffGeo
// Purpose       :
//
// Special Notes : nf_arg is named to distinguish it from the
//                 instance variable nf.
//
// Scope         : public
// Creator       : Eric Keiter,SNL
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
int Instance::PAeffGeo
  (double nf_arg, int geo, int minSD,
   double Weffcj, double DMCG, double DMCI, double DMDG,
   double & Ps, double & Pd, double & As, double & Ad)
{
  double T0(0.0), T1(0.0), T2(0.0);
  double ADiso(0.0), ADsha(0.0), ADmer(0.0), ASiso(0.0), ASsha(0.0), ASmer(0.0);
  double PDiso(0.0), PDsha(0.0), PDmer(0.0), PSiso(0.0), PSsha(0.0), PSmer(0.0);
  double nuIntD (0.0), nuEndD (0.0), nuIntS (0.0), nuEndS (0.0);
  
  if (geo < 9) // For geo = 9 and 10, the numbers of S/D diffusions already known
    NumFingerDiff(nf_arg, minSD, nuIntD, nuEndD, nuIntS, nuEndS);
  
  T0 = DMCG + DMCI;
  T1 = DMCG + DMCG;
  T2 = DMDG + DMDG;
  
  PSiso = PDiso = T0 + T0 + Weffcj;
  PSsha = PDsha = T1;
  PSmer = PDmer = T2;
  
  ASiso = ADiso = T0 * Weffcj;
  ASsha = ADsha = DMCG * Weffcj;
  ASmer = ADmer = DMDG * Weffcj;
  
  switch(geo)
  {
    case 0:
        Ps = nuEndS * PSiso + nuIntS * PSsha;
        Pd = nuEndD * PDiso + nuIntD * PDsha;
        As = nuEndS * ASiso + nuIntS * ASsha;
        Ad = nuEndD * ADiso + nuIntD * ADsha;
        break;
    case 1:
        Ps = nuEndS * PSiso + nuIntS * PSsha;
        Pd = (nuEndD + nuIntD) * PDsha;
        As = nuEndS * ASiso + nuIntS * ASsha;
        Ad = (nuEndD + nuIntD) * ADsha;
        break;
    case 2:
        Ps = (nuEndS + nuIntS) * PSsha;
        Pd = nuEndD * PDiso + nuIntD * PDsha;
        As = (nuEndS + nuIntS) * ASsha;
        Ad = nuEndD * ADiso + nuIntD * ADsha;
        break;
    case 3:
        Ps = (nuEndS + nuIntS) * PSsha;
        Pd = (nuEndD + nuIntD) * PDsha;
        As = (nuEndS + nuIntS) * ASsha;
        Ad = (nuEndD + nuIntD) * ADsha;
        break;
    case 4:
        Ps = nuEndS * PSiso + nuIntS * PSsha;
        Pd = nuEndD * PDmer + nuIntD * PDsha;
        As = nuEndS * ASiso + nuIntS * ASsha;
        Ad = nuEndD * ADmer + nuIntD * ADsha;
        break;
    case 5:
        Ps = (nuEndS + nuIntS) * PSsha;
        Pd = nuEndD * PDmer + nuIntD * PDsha;
        As = (nuEndS + nuIntS) * ASsha;
        Ad = nuEndD * ADmer + nuIntD * ADsha;
        break;
    case 6:
        Ps = nuEndS * PSmer + nuIntS * PSsha;
        Pd = nuEndD * PDiso + nuIntD * PDsha;
        As = nuEndS * ASmer + nuIntS * ASsha;
        Ad = nuEndD * ADiso + nuIntD * ADsha;
        break;
    case 7:
        Ps = nuEndS * PSmer + nuIntS * PSsha;
        Pd = (nuEndD + nuIntD) * PDsha;
        As = nuEndS * ASmer + nuIntS * ASsha;
        Ad = (nuEndD + nuIntD) * ADsha;
        break;
    case 8:
        Ps = nuEndS * PSmer + nuIntS * PSsha;
        Pd = nuEndD * PDmer + nuIntD * PDsha;
        As = nuEndS * ASmer + nuIntS * ASsha;
        Ad = nuEndD * ADmer + nuIntD * ADsha;
        break;
    case 9: // geo = 9 and 10 happen only when nf_arg = even
        Ps = PSiso + (nf_arg - 1.0) * PSsha;
        Pd = nf_arg * PDsha;
        As = ASiso + (nf_arg - 1.0) * ASsha;
        Ad = nf_arg * ADsha;
        break;
    case 10:
        Ps = nf_arg * PSsha;
        Pd = PDiso + (nf_arg - 1.0) * PDsha;
        As = nf_arg * ASsha;
        Ad = ADiso + (nf_arg - 1.0) * ADsha;
        break;
    default:
      UserWarning(*this) << "Specified GEO not matched\n";
  }
  
  return 0;
}


//-----------------------------------------------------------------------------
// Function      : Instance::RdseffGeo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter,SNL
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
int Instance::RdseffGeo
  (double nf_arg,
   int geo, int rgeo, int minSD,
   double Weffcj, double Rsh, double DMCG, double DMCI, double DMDG,
   int Type, double & Rtot)
{
  std::string msg="";
  double Rint(0.0), Rend (0.0);
  double nuIntD (0.0), nuEndD (0.0), nuIntS (0.0), nuEndS (0.0);

  if (geo < 9) // since geo = 9 and 10 only happen when nf_arg = even
  {
    NumFingerDiff(nf_arg, minSD, nuIntD, nuEndD, nuIntS, nuEndS);

    // Internal S/D resistance -- assume shared S or D and all wide contacts
    if (Type == 1)
    {
      if (nuIntS == 0.0)
        Rint = 0.0;
      else
        Rint = Rsh * DMCG / ( Weffcj * nuIntS);
    }
    else
    {
      if (nuIntD == 0.0)
        Rint = 0.0;
      else
        Rint = Rsh * DMCG / ( Weffcj * nuIntD);
    }
  }

  // End S/D resistance  -- geo dependent
  switch(geo)
  {   case 0:
      if (Type == 1) RdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndS, rgeo, 1, Rend);
      else           RdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndD, rgeo, 0, Rend);
        break;
    case 1:
        if (Type == 1) RdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndS, rgeo, 1, Rend);
        else           RdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG, static_cast<int>(nuEndD), rgeo, 0, Rend);
        break;
    case 2:
        if (Type == 1) RdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG, static_cast<int>(nuEndS), rgeo, 1, Rend);
        else           RdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndD, rgeo, 0, Rend);
        break;
    case 3:
        if (Type == 1) RdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG, static_cast<int>(nuEndS), rgeo, 1, Rend);
        else           RdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG, static_cast<int>(nuEndD), rgeo, 0, Rend);
        break;
    case 4:
        if (Type == 1) RdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndS, rgeo, 1, Rend);
        else           Rend = Rsh * DMDG / Weffcj;
        break;
    case 5:
        if (Type == 1) RdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG, static_cast<int>(nuEndS), rgeo, 1, Rend);
        else           Rend = Rsh * DMDG / (Weffcj * nuEndD);
        break;
    case 6:
        if (Type == 1) Rend = Rsh * DMDG / Weffcj;
        else           RdsEndIso(Weffcj, Rsh, DMCG, DMCI, DMDG, nuEndD, rgeo, 0, Rend);
        break;
    case 7:
        if (Type == 1) Rend = Rsh * DMDG / (Weffcj * nuEndS);
        else           RdsEndSha(Weffcj, Rsh, DMCG, DMCI, DMDG, static_cast<int>(nuEndD), rgeo, 0, Rend);
        break;
    case 8:
        Rend = Rsh * DMDG / Weffcj;
        break;
    case 9: // all wide contacts assumed for geo = 9 and 10
        if (Type == 1)
        {
          Rend = 0.5 * Rsh * DMCG / Weffcj;
          if (nf_arg == 2.0)
            Rint = 0.0;
          else
            Rint = Rsh * DMCG / (Weffcj * (nf_arg - 2.0));
        }
        else
        {
          Rend = 0.0;
          Rint = Rsh * DMCG / (Weffcj * nf_arg);
        }
        break;
    case 10:
        if (Type == 1)
        {
          Rend = 0.0;
          Rint = Rsh * DMCG / (Weffcj * nf_arg);
        }
        else
        {
          Rend = 0.5 * Rsh * DMCG / Weffcj;;
          if (nf_arg == 2.0)
            Rint = 0.0;
          else
            Rint = Rsh * DMCG / (Weffcj * (nf_arg - 2.0));
        }
        break;
    default:
        UserWarning(*this) << "Specified GEO not matched\n";
  }

  if (Rint <= 0.0)
    Rtot = Rend;
  else if (Rend <= 0.0)
    Rtot = Rint;
  else
    Rtot = Rint * Rend / (Rint + Rend);

  if(Rtot==0.0)
  {
    UserWarning(*this) << "Zero resistance returned from RdseffGeo\n";
  }

  return 0;
}


//-----------------------------------------------------------------------------
// Function      : Instance::RdsEndIso
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter,SNL
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
int Instance::RdsEndIso
  (double Weffcj, double Rsh, double DMCG, double DMCI, double DMDG,
   double nuEnd, int rgeo, int Type, double & Rend)
{
  std::string msg="";
	if (Type == 1)
	{
    switch(rgeo)
    {
      case 1:
      case 2:
      case 5:
          if (nuEnd == 0.0)
              Rend = 0.0;
          else
              Rend = Rsh * DMCG / (Weffcj * nuEnd);
          break;
      case 3:
      case 4:
      case 6:
          if ((DMCG + DMCI) == 0.0)
            msg = "(DMCG + DMCI) can not be equal to zero\n";
          if (nuEnd == 0.0)
            Rend = 0.0;
          else
            Rend = Rsh * Weffcj / (3.0 * nuEnd * (DMCG + DMCI));
          break;
      default:
          UserWarning(*this) << "Specified RGEO not matched\n";
    }
	}
	else
	{
    switch(rgeo)
    {
      case 1:
      case 3:
      case 7:
        if (nuEnd == 0.0)
          Rend = 0.0;
        else
          Rend = Rsh * DMCG / (Weffcj * nuEnd);
        break;
      case 2:
      case 4:
      case 8:
        if ((DMCG + DMCI) == 0.0)
          msg = "(DMCG + DMCI) can not be equal to zero\n";
        if (nuEnd == 0.0)
          Rend = 0.0;
        else
          Rend = Rsh * Weffcj / (3.0 * nuEnd * (DMCG + DMCI));
        break;
      default:
        UserWarning(*this) << "Specified RGEO not matched\n";
            }
	}
  return 0;
}


//-----------------------------------------------------------------------------
// Function      : Instance::RdsEndSha
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter,SNL
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
int Instance::RdsEndSha
   (double Weffcj, double Rsh, double DMCG, double DMCI, double DMDG,
    int rgeo, int Type, double nuEnd, double & Rend)
{
  std::string msg = "";
  if (Type == 1)
  {
    switch(rgeo)
    {
      case 1:
      case 2:
      case 5:
          if (nuEnd == 0.0)
              Rend = 0.0;
          else
              Rend = Rsh * DMCG / (Weffcj * nuEnd);
          break;
      case 3:
      case 4:
      case 6:
          if (DMCG == 0.0)
              msg = "DMCG can not be equal to zero\n";
          if (nuEnd == 0.0)
              Rend = 0.0;
          else
              Rend = Rsh * Weffcj / (6.0 * nuEnd * DMCG);
          break;
      default:
          UserWarning(*this) << "Specified RGEO not matched\n";
    }
  }
  else
  {
    switch(rgeo)
    {
      case 1:
      case 3:
      case 7:
          if (nuEnd == 0.0)
              Rend = 0.0;
          else
              Rend = Rsh * DMCG / (Weffcj * nuEnd);
          break;
      case 2:
      case 4:
      case 8:
          if (DMCG == 0.0)
              msg = "DMCG can not be equal to zero\n";
          if (nuEnd == 0.0)
              Rend = 0.0;
          else
              Rend = Rsh * Weffcj / (6.0 * nuEnd * DMCG);
          break;
      default:
          UserWarning(*this) << "Specified RGEO = %d not matched\n";
    }
  }
  return 0;
}

// Additional Declarations

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
bool Model::processParams ()
{
  std::string msg;

  if (SbulkJctPotential < 0.1)
  {
    SbulkJctPotential = 0.1;
    UserWarning(*this) << "Given pbs is less than 0.1. Pbs is set to 0.1";
  }
  if (SsidewallJctPotential < 0.1)
  {
    SsidewallJctPotential = 0.1;
     UserWarning(*this) << "Given pbsws is less than 0.1. Pbsws is set to 0.1";
  }
  if (SGatesidewallJctPotential < 0.1)
  {
    SGatesidewallJctPotential = 0.1;
     UserWarning(*this) << "Given pbswgs is less than 0.1. Pbswgs is set to 0.1";
  }

  if (DbulkJctPotential < 0.1)
  {
    DbulkJctPotential = 0.1;
     UserWarning(*this) << "Given pbd is less than 0.1. Pbd is set to 0.1";
  }
  if (DsidewallJctPotential < 0.1)
  {
    DsidewallJctPotential = 0.1;
     UserWarning(*this) << "Given pbswd is less than 0.1. Pbswd is set to 0.1";
  }
  if (DGatesidewallJctPotential < 0.1)
  {
    DGatesidewallJctPotential = 0.1;
    UserWarning(*this) << "Given pbswgd is less than 0.1. Pbswgd is set to 0.1";
  }

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
// Creator       : Eric Keiter, SNL
// Creation Date : 11/25/06
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
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    instanceContainer(0),
    modType       (0),
    dtype         (CONSTNMOS),
    mobMod(0),
    cvchargeMod(0),
    capMod(0),
    dioMod(0),
    trnqsMod(0),
    acnqsMod(0),
    fnoiMod(0),
    tnoiMod(0),
    rdsMod(0),
    rbodyMod(0),
    rgateMod(0),
    perMod(0),
    geoMod(0),
    mtrlMod(0),
    igcMod(0),
    igbMod(0),
    tempMod(0),
    binUnit(0),
    paramChk(0),
    version("4.6.1"),
    eot(0.0),
    vddeot(0.0),
    ados(0.0),
    bdos(0.0),
    toxe(0.0),
    toxp(0.0),
    toxm(0.0),
    dtox(0.0),
    epsrox(0.0),
    cdsc(0.0),
    cdscb(0.0),
    cdscd(0.0),
    cit(0.0),
    nfactor(0.0),
    xj(0.0),
    vsat(0.0),
    at(0.0),
    a0(0.0),
    ags(0.0),
    a1(0.0),
    a2(0.0),
    keta(0.0),
    nsub(0.0),
    phig(0.0),
    epsrgate(0.0),
    easub(0.0),
    epsrsub(0.0),
    ni0sub(0.0),
    bg0sub(0.0),
    tbgasub(0.0),
    tbgbsub(0.0),
    ndep(0.0),
    nsd(0.0),
    phin(0.0),
    ngate(0.0),
    gamma1(0.0),
    gamma2(0.0),
    vbx(0.0),
    vbm(0.0),
    xt(0.0),
    k1(0.0),
    kt1(0.0),
    kt1l(0.0),
    kt2(0.0),
    k2(0.0),
    k3(0.0),
    k3b(0.0),
    w0(0.0),
    dvtp0(0.0),
    dvtp1(0.0),
    lpe0(0.0),
    lpeb(0.0),
    dvt0(0.0),
    dvt1(0.0),
    dvt2(0.0),
    dvt0w(0.0),
    dvt1w(0.0),
    dvt2w(0.0),
    drout(0.0),
    dsub(0.0),
    vth0(0.0),
    eu(0.0),
    ua(0.0),
    ua1(0.0),
    ub(0.0),
    ub1(0.0),
    uc(0.0),
    uc1(0.0),
    ud(0.0),
    ud1(0.0),
    up(0.0),
    lp(0.0),
    u0(0.0),
    ute(0.0),
    voff(0.0),
    tvoff(0.0),
    minv(0.0),
    minvcv(0.0),
    voffl(0.0),
    voffcvl(0.0),
    delta(0.0),
    rdsw(0.0),
    rdswmin(0.0),
    rdwmin(0.0),
    rswmin(0.0),
    rsw(0.0),
    rdw(0.0),
    prwg(0.0),
    prwb(0.0),
    prt(0.0),
    eta0(0.0),
    etab(0.0),
    pclm(0.0),
    pdibl1(0.0),
    pdibl2(0.0),
    pdiblb(0.0),
    fprout(0.0),
    pdits(0.0),
    pditsd(0.0),
    pditsl(0.0),
    pscbe1(0.0),
    pscbe2(0.0),
    pvag(0.0),
    wr(0.0),
    dwg(0.0),
    dwb(0.0),
    b0(0.0),
    b1(0.0),
    alpha0(0.0),
    alpha1(0.0),
    beta0(0.0),
    agidl(0.0),
    bgidl(0.0),
    cgidl(0.0),
    egidl(0.0),
    agisl(0.0),
    bgisl(0.0),
    cgisl(0.0),
    egisl(0.0),
    aigc(0.0),
    bigc(0.0),
    cigc(0.0),
    aigsd(0.0),
    bigsd(0.0),
    cigsd(0.0),
    aigs(0.0),
    bigs(0.0),
    cigs(0.0),
    aigd(0.0),
    bigd(0.0),
    cigd(0.0),
    aigbacc(0.0),
    bigbacc(0.0),
    cigbacc(0.0),
    aigbinv(0.0),
    bigbinv(0.0),
    cigbinv(0.0),
    nigc(0.0),
    nigbacc(0.0),
    nigbinv(0.0),
    ntox(0.0),
    eigbinv(0.0),
    pigcd(0.0),
    poxedge(0.0),
    toxref(0.0),
    ijthdfwd(0.0),
    ijthsfwd(0.0),
    ijthdrev(0.0),
    ijthsrev(0.0),
    xjbvd(0.0),
    xjbvs(0.0),
    bvd(0.0),
    bvs(0.0),
    jtss(0.0),
    jtsd(0.0),
    jtssws(0.0),
    jtsswd(0.0),
    jtsswgs(0.0),
    jtsswgd(0.0),
    njts(0.0),
    njtssw(0.0),
    njtsswg(0.0),
    njtsd(0.0),
    njtsswd(0.0),
    njtsswgd(0.0),
    xtss(0.0),
    xtsd(0.0),
    xtssws(0.0),
    xtsswd(0.0),
    xtsswgs(0.0),
    xtsswgd(0.0),
    tnjts(0.0),
    tnjtssw(0.0),
    tnjtsswg(0.0),
    tnjtsd(0.0),
    tnjtsswd(0.0),
    tnjtsswgd(0.0),
    vtss(0.0),
    vtsd(0.0),
    vtssws(0.0),
    vtsswd(0.0),
    vtsswgs(0.0),
    vtsswgd(0.0),

    xrcrg1(0.0),
    xrcrg2(0.0),
    lambda(0.0),
    vtl(0.0),
    lc(0.0),
    xn(0.0),
    vfbsdoff(0.0),  // S/D flatband offset voltage
    lintnoi(0.0),  // lint offset for noise calculation
    tvfbsdoff(0.0),

    vfb(0.0),
    gbmin(0.0),
    rbdb(0.0),
    rbsb(0.0),
    rbpb(0.0),
    rbps(0.0),
    rbpd(0.0),

    rbps0(0.0),
    rbpsl(0.0),
    rbpsw(0.0),
    rbpsnf(0.0),

    rbpd0(0.0),
    rbpdl(0.0),
    rbpdw(0.0),
    rbpdnf(0.0),

    rbpbx0(0.0),
    rbpbxl(0.0),
    rbpbxw(0.0),
    rbpbxnf(0.0),
    rbpby0(0.0),
    rbpbyl(0.0),
    rbpbyw(0.0),
    rbpbynf(0.0),

    rbsbx0(0.0),
    rbsby0(0.0),
    rbdbx0(0.0),
    rbdby0(0.0),

    rbsdbxl(0.0),
    rbsdbxw(0.0),
    rbsdbxnf(0.0),
    rbsdbyl(0.0),
    rbsdbyw(0.0),
    rbsdbynf(0.0),

    tnoia(0.0),
    tnoib(0.0),
    rnoia(0.0),
    rnoib(0.0),
    ntnoi(0.0),

    // CV model and Parasitics
    cgsl(0.0),
    cgdl(0.0),
    ckappas(0.0),
    ckappad(0.0),
    cf(0.0),
    vfbcv(0.0),
    clc(0.0),
    cle(0.0),
    dwc(0.0),
    dlc(0.0),
    xw(0.0),
    xl(0.0),
    dlcig(0.0),
    dlcigd(0.0),
    dwj(0.0),
    noff(0.0),
    voffcv(0.0),
    acde(0.0),
    moin(0.0),
    tcj(0.0),
    tcjsw(0.0),
    tcjswg(0.0),
    tpb(0.0),
    tpbsw(0.0),
    tpbswg(0.0),
    dmcg(0.0),
    dmci(0.0),
    dmdg(0.0),
    dmcgt(0.0),
    xgw(0.0),
    xgl(0.0),
    rshg(0.0),
    ngcon(0.0),

    // Length Dependence
    lcdsc(0.0),
    lcdscb(0.0),
    lcdscd(0.0),
    lcit(0.0),
    lnfactor(0.0),
    lxj(0.0),
    lvsat(0.0),
    lat(0.0),
    la0(0.0),
    lags(0.0),
    la1(0.0),
    la2(0.0),
    lketa(0.0),
    lnsub(0.0),
    lndep(0.0),
    lnsd(0.0),
    lphin(0.0),
    lngate(0.0),
    lgamma1(0.0),
    lgamma2(0.0),
    lvbx(0.0),
    lvbm(0.0),
    lxt(0.0),
    lk1(0.0),
    lkt1(0.0),
    lkt1l(0.0),
    lkt2(0.0),
    lk2(0.0),
    lk3(0.0),
    lk3b(0.0),
    lw0(0.0),
    ldvtp0(0.0),
    ldvtp1(0.0),
    llpe0(0.0),
    llpeb(0.0),
    ldvt0(0.0),
    ldvt1(0.0),
    ldvt2(0.0),
    ldvt0w(0.0),
    ldvt1w(0.0),
    ldvt2w(0.0),
    ldrout(0.0),
    ldsub(0.0),
    lvth0(0.0),
    lua(0.0),
    lua1(0.0),
    lub(0.0),
    lub1(0.0),
    luc(0.0),
    luc1(0.0),
    lud(0.0),
    lud1(0.0),
    lup(0.0),
    llp(0.0),
    lu0(0.0),
    leu(0.0),
    lute(0.0),
    lvoff(0.0),
    ltvoff(0.0),
    lminv(0.0),
    lminvcv(0.0),
    ldelta(0.0),
    lrdsw(0.0),
    lrsw(0.0),
    lrdw(0.0),
    lprwg(0.0),
    lprwb(0.0),
    lprt(0.0),
    leta0(0.0),
    letab(0.0),
    lpclm(0.0),
    lpdibl1(0.0),
    lpdibl2(0.0),
    lpdiblb(0.0),
    lfprout(0.0),
    lpdits(0.0),
    lpditsd(0.0),
    lpscbe1(0.0),
    lpscbe2(0.0),
    lpvag(0.0),
    lwr(0.0),
    ldwg(0.0),
    ldwb(0.0),
    lb0(0.0),
    lb1(0.0),
    lalpha0(0.0),
    lalpha1(0.0),
    lbeta0(0.0),
    lvfb(0.0),
    lagidl(0.0),
    lbgidl(0.0),
    lcgidl(0.0),
    legidl(0.0),
    lagisl(0.0),
    lbgisl(0.0),
    lcgisl(0.0),
    legisl(0.0),
    laigc(0.0),
    lbigc(0.0),
    lcigc(0.0),
    laigsd(0.0),
    lbigsd(0.0),
    lcigsd(0.0),
    laigs(0.0),
    lbigs(0.0),
    lcigs(0.0),
    laigd(0.0),
    lbigd(0.0),
    lcigd(0.0),
    laigbacc(0.0),
    lbigbacc(0.0),
    lcigbacc(0.0),
    laigbinv(0.0),
    lbigbinv(0.0),
    lcigbinv(0.0),
    lnigc(0.0),
    lnigbacc(0.0),
    lnigbinv(0.0),
    lntox(0.0),
    leigbinv(0.0),
    lpigcd(0.0),
    lpoxedge(0.0),
    lxrcrg1(0.0),
    lxrcrg2(0.0),
    llambda(0.0),
    lvtl(0.0),
    lxn(0.0),
    lvfbsdoff(0.0),
    ltvfbsdoff(0.0),

    // CV model
    lcgsl(0.0),
    lcgdl(0.0),
    lckappas(0.0),
    lckappad(0.0),
    lcf(0.0),
    lclc(0.0),
    lcle(0.0),
    lvfbcv(0.0),
    lnoff(0.0),
    lvoffcv(0.0),
    lacde(0.0),
    lmoin(0.0),

    // Width Dependence
    wcdsc(0.0),
    wcdscb(0.0),
    wcdscd(0.0),
    wcit(0.0),
    wnfactor(0.0),
    wxj(0.0),
    wvsat(0.0),
    wat(0.0),
    wa0(0.0),
    wags(0.0),
    wa1(0.0),
    wa2(0.0),
    wketa(0.0),
    wnsub(0.0),
    wndep(0.0),
    wnsd(0.0),
    wphin(0.0),
    wngate(0.0),
    wgamma1(0.0),
    wgamma2(0.0),
    wvbx(0.0),
    wvbm(0.0),
    wxt(0.0),
    wk1(0.0),
    wkt1(0.0),
    wkt1l(0.0),
    wkt2(0.0),
    wk2(0.0),
    wk3(0.0),
    wk3b(0.0),
    ww0(0.0),
    wdvtp0(0.0),
    wdvtp1(0.0),
    wlpe0(0.0),
    wlpeb(0.0),
    wdvt0(0.0),
    wdvt1(0.0),
    wdvt2(0.0),
    wdvt0w(0.0),
    wdvt1w(0.0),
    wdvt2w(0.0),
    wdrout(0.0),
    wdsub(0.0),
    wvth0(0.0),
    wua(0.0),
    wua1(0.0),
    wub(0.0),
    wub1(0.0),
    wuc(0.0),
    wuc1(0.0),
    wud(0.0),
    wud1(0.0),
    wup(0.0),
    wlp(0.0),
    wu0(0.0),
    weu(0.0),
    wute(0.0),
    wvoff(0.0),
    wtvoff(0.0),
    wminv(0.0),
    wminvcv(0.0),
    wdelta(0.0),
    wrdsw(0.0),
    wrsw(0.0),
    wrdw(0.0),
    wprwg(0.0),
    wprwb(0.0),
    wprt(0.0),
    weta0(0.0),
    wetab(0.0),
    wpclm(0.0),
    wpdibl1(0.0),
    wpdibl2(0.0),
    wpdiblb(0.0),
    wfprout(0.0),
    wpdits(0.0),
    wpditsd(0.0),
    wpscbe1(0.0),
    wpscbe2(0.0),
    wpvag(0.0),
    wwr(0.0),
    wdwg(0.0),
    wdwb(0.0),
    wb0(0.0),
    wb1(0.0),
    walpha0(0.0),
    walpha1(0.0),
    wbeta0(0.0),
    wvfb(0.0),
    wagidl(0.0),
    wbgidl(0.0),
    wcgidl(0.0),
    wegidl(0.0),
    wagisl(0.0),
    wbgisl(0.0),
    wcgisl(0.0),
    wegisl(0.0),
    waigc(0.0),
    wbigc(0.0),
    wcigc(0.0),
    waigsd(0.0),
    wbigsd(0.0),
    wcigsd(0.0),
    waigs(0.0),
    wbigs(0.0),
    wcigs(0.0),
    waigd(0.0),
    wbigd(0.0),
    wcigd(0.0),
    waigbacc(0.0),
    wbigbacc(0.0),
    wcigbacc(0.0),
    waigbinv(0.0),
    wbigbinv(0.0),
    wcigbinv(0.0),
    wnigc(0.0),
    wnigbacc(0.0),
    wnigbinv(0.0),
    wntox(0.0),
    weigbinv(0.0),
    wpigcd(0.0),
    wpoxedge(0.0),
    wxrcrg1(0.0),
    wxrcrg2(0.0),
    wlambda(0.0),
    wvtl(0.0),
    wxn(0.0),
    wvfbsdoff(0.0),
    wtvfbsdoff(0.0),

    // CV model
    wcgsl(0.0),
    wcgdl(0.0),
    wckappas(0.0),
    wckappad(0.0),
    wcf(0.0),
    wclc(0.0),
    wcle(0.0),
    wvfbcv(0.0),
    wnoff(0.0),
    wvoffcv(0.0),
    wacde(0.0),
    wmoin(0.0),

    // Cross-term Dependence
    pcdsc(0.0),
    pcdscb(0.0),
    pcdscd(0.0),
    pcit(0.0),
    pnfactor(0.0),
    pxj(0.0),
    pvsat(0.0),
    pat(0.0),
    pa0(0.0),
    pags(0.0),
    pa1(0.0),
    pa2(0.0),
    pketa(0.0),
    pnsub(0.0),
    pndep(0.0),
    pnsd(0.0),
    pphin(0.0),
    pngate(0.0),
    pgamma1(0.0),
    pgamma2(0.0),
    pvbx(0.0),
    pvbm(0.0),
    pxt(0.0),
    pk1(0.0),
    pkt1(0.0),
    pkt1l(0.0),
    pkt2(0.0),
    pk2(0.0),
    pk3(0.0),
    pk3b(0.0),
    pw0(0.0),
    pdvtp0(0.0),
    pdvtp1(0.0),
    plpe0(0.0),
    plpeb(0.0),
    pdvt0(0.0),
    pdvt1(0.0),
    pdvt2(0.0),
    pdvt0w(0.0),
    pdvt1w(0.0),
    pdvt2w(0.0),
    pdrout(0.0),
    pdsub(0.0),
    pvth0(0.0),
    pua(0.0),
    pua1(0.0),
    pub(0.0),
    pub1(0.0),
    puc(0.0),
    puc1(0.0),
    pud(0.0),
    pud1(0.0),
    pup(0.0),
    plp(0.0),
    pu0(0.0),
    peu(0.0),
    pute(0.0),
    pvoff(0.0),
    ptvoff(0.0),
    pminv(0.0),
    pminvcv(0.0),
    pdelta(0.0),
    prdsw(0.0),
    prsw(0.0),
    prdw(0.0),
    pprwg(0.0),
    pprwb(0.0),
    pprt(0.0),
    peta0(0.0),
    petab(0.0),
    ppclm(0.0),
    ppdibl1(0.0),
    ppdibl2(0.0),
    ppdiblb(0.0),
    pfprout(0.0),
    ppdits(0.0),
    ppditsd(0.0),
    ppscbe1(0.0),
    ppscbe2(0.0),
    ppvag(0.0),
    pwr(0.0),
    pdwg(0.0),
    pdwb(0.0),
    pb0(0.0),
    pb1(0.0),
    palpha0(0.0),
    palpha1(0.0),
    pbeta0(0.0),
    pvfb(0.0),
    pagidl(0.0),
    pbgidl(0.0),
    pcgidl(0.0),
    pegidl(0.0),
    pagisl(0.0),
    pbgisl(0.0),
    pcgisl(0.0),
    pegisl(0.0),
    paigc(0.0),
    pbigc(0.0),
    pcigc(0.0),
    paigsd(0.0),
    pbigsd(0.0),
    pcigsd(0.0),
    paigs(0.0),
    pbigs(0.0),
    pcigs(0.0),
    paigd(0.0),
    pbigd(0.0),
    pcigd(0.0),
    paigbacc(0.0),
    pbigbacc(0.0),
    pcigbacc(0.0),
    paigbinv(0.0),
    pbigbinv(0.0),
    pcigbinv(0.0),
    pnigc(0.0),
    pnigbacc(0.0),
    pnigbinv(0.0),
    pntox(0.0),
    peigbinv(0.0),
    ppigcd(0.0),
    ppoxedge(0.0),
    pxrcrg1(0.0),
    pxrcrg2(0.0),
    plambda(0.0),
    pvtl(0.0),
    pxn(0.0),
    pvfbsdoff(0.0),
    ptvfbsdoff(0.0),

    // CV model
    pcgsl(0.0),
    pcgdl(0.0),
    pckappas(0.0),
    pckappad(0.0),
    pcf(0.0),
    pclc(0.0),
    pcle(0.0),
    pvfbcv(0.0),
    pnoff(0.0),
    pvoffcv(0.0),
    pacde(0.0),
    pmoin(0.0),

    tnom(0.0),
    cgso(0.0),
    cgdo(0.0),
    cgbo(0.0),
    xpart(0.0),
    cFringOut(0.0),
    cFringMax(0.0),

    sheetResistance(0.0),
    SjctSatCurDensity(0.0),
    DjctSatCurDensity(0.0),
    SjctSidewallSatCurDensity(0.0),
    DjctSidewallSatCurDensity(0.0),
    SjctGateSidewallSatCurDensity(0.0),
    DjctGateSidewallSatCurDensity(0.0),
    SbulkJctPotential(0.0),
    DbulkJctPotential(0.0),
    SbulkJctBotGradingCoeff(0.0),
    DbulkJctBotGradingCoeff(0.0),
    SbulkJctSideGradingCoeff(0.0),
    DbulkJctSideGradingCoeff(0.0),
    SbulkJctGateSideGradingCoeff(0.0),
    DbulkJctGateSideGradingCoeff(0.0),
    SsidewallJctPotential(0.0),
    DsidewallJctPotential(0.0),
    SGatesidewallJctPotential(0.0),
    DGatesidewallJctPotential(0.0),
    SunitAreaJctCap(0.0),
    DunitAreaJctCap(0.0),
    SunitLengthSidewallJctCap(0.0),
    DunitLengthSidewallJctCap(0.0),
    SunitLengthGateSidewallJctCap(0.0),
    DunitLengthGateSidewallJctCap(0.0),
    SjctEmissionCoeff(0.0),
    DjctEmissionCoeff(0.0),
    SjctTempExponent(0.0),
    DjctTempExponent(0.0),
    njtsstemp(0.0),
    njtsswstemp(0.0),
    njtsswgstemp(0.0),
    njtsdtemp(0.0),
    njtsswdtemp(0.0),
    njtsswgdtemp(0.0),


    Lint(0.0),
    Ll(0.0),
    Llc(0.0),
    Lln(0.0),
    Lw(0.0),
    Lwc(0.0),
    Lwn(0.0),
    Lwl(0.0),
    Lwlc(0.0),
    Lmin(0.0),
    Lmax(0.0),

    Wint(0.0),
    Wl(0.0),
    Wlc(0.0),
    Wln(0.0),
    Ww(0.0),
    Wwc(0.0),
    Wwn(0.0),
    Wwl(0.0),
    Wwlc(0.0),
    Wmin(0.0),
    Wmax(0.0),

    // added for stress effect
    saref(0.0),
    sbref(0.0),
    wlod(0.0),
    ku0(0.0),
    kvsat(0.0),
    kvth0(0.0),
    tku0(0.0),
    llodku0(0.0),
    wlodku0(0.0),
    llodvth(0.0),
    wlodvth(0.0),
    lku0(0.0),
    wku0(0.0),
    pku0(0.0),
    lkvth0(0.0),
    wkvth0(0.0),
    pkvth0(0.0),
    stk2(0.0),
    lodk2(0.0),
    steta0(0.0),
    lodeta0(0.0),

    web(0.0),
    wec(0.0),
    kvth0we(0.0),
    k2we(0.0),
    ku0we(0.0),
    scref(0.0),
    wpemod(0.0),
    lkvth0we(0.0),
    lk2we(0.0),
    lku0we(0.0),
    wkvth0we(0.0),
    wk2we(0.0),
    wku0we(0.0),
    pkvth0we(0.0),
    pk2we(0.0),
    pku0we(0.0),

// Pre-calculated constants
// move to size-dependent param
    Eg0(0.0),
    vtm(0.0),
    vtm0(0.0),
    coxe(0.0),
    coxp(0.0),
    cof1(0.0),
    cof2(0.0),
    cof3(0.0),
    cof4(0.0),
    vcrit(0.0),
    factor1(0.0),
    PhiBS(0.0),
    PhiBSWS(0.0),
    PhiBSWGS(0.0),
    SjctTempSatCurDensity(0.0),
    SjctSidewallTempSatCurDensity(0.0),
    SjctGateSidewallTempSatCurDensity(0.0),
    PhiBD(0.0),
    PhiBSWD(0.0),
    PhiBSWGD(0.0),
    DjctTempSatCurDensity(0.0),
    DjctSidewallTempSatCurDensity(0.0),
    DjctGateSidewallTempSatCurDensity(0.0),
    SunitAreaTempJctCap(0.0),
    DunitAreaTempJctCap(0.0),
    SunitLengthSidewallTempJctCap(0.0),
    DunitLengthSidewallTempJctCap(0.0),
    SunitLengthGateSidewallTempJctCap(0.0),
    DunitLengthGateSidewallTempJctCap(0.0),

    oxideTrapDensityA(0.0),
    oxideTrapDensityB(0.0),
    oxideTrapDensityC(0.0),
    em(0.0),
    ef(0.0),
    af(0.0),
    kf(0.0),
    ni(0.0),
    Vtm0(0.0),

    vtlGiven(false),
    ndepGiven(false),
    gamma1Given(false),
    k1Given(false),
    k2Given(false),
    nsubGiven(false),
    phigGiven(false),
    xtGiven(false),
    vbxGiven(false),
    gamma2Given(false),
    vfbGiven(false),
    vth0Given(false),
    rbps0Given(false),
    rbpd0Given(false),
    rbsbx0Given(false),
    rbsby0Given(false),
    rbdbx0Given(false),
    rbdby0Given(false),
    lambdaGiven(false),
    pigcdGiven(false),
    toxeGiven(false),
    toxpGiven(false),
    dtoxGiven(false),
    cgdoGiven(false),
    dlcGiven(false),
    cgsoGiven(false),
    cgboGiven(false),
    sizeDependParamList()

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

  // Calculate any parameters specified as expressions:

  updateDependentParameters();

  // calculate dependent (ie computed) params
  if (!given("TOXP") )
    toxp = toxe;
  if (!given("TOXM") )
    toxm = toxe;
  if (!given("DSUB") )
    dsub = drout;

  if (!given("VTH0") )
    vth0 = (dtype == CONSTNMOS) ? 0.7 : -0.7;

  if (!given("VDDEOT"))
    vddeot= (dtype == CONSTNMOS)?1.5:-1.5;
  if (!given("EU"))
    eu =(dtype == CONSTNMOS) ? 1.67 : 1.0;;
  if (!given("UA"))
    ua =(mobMod == 2) ? 1.0E-15 : 1.0E-9; // UNIT M/V
  if (!given("UC"))
    uc = (mobMod == 1) ? -0.0465 : -0.0465E-9;
  if (!given("UC1"))
    uc1 =(mobMod == 1) ? -0.056 : -0.056E-9;
  if (!given("U0"))
    u0 = (dtype == CONSTNMOS) ? 0.067 : 0.025;
  if (!given("AIGC"))
    aigc =(dtype == CONSTNMOS) ? 1.36E-2 : 9.80E-3;
  if (!given("BIGC"))
    bigc =(dtype == CONSTNMOS) ? 1.71E-3 : 7.59E-4;
  if (!given("CIGC"))
    cigc = (dtype == CONSTNMOS) ? 0.075 : 0.03;
  if (given("AIGSD"))
  {
    aigs = aigd = aigsd;
  }
  else
  {
    aigsd = (dtype == CONSTNMOS) ? 1.36E-2 : 9.80E-3;
    if (!given("AIGS"))
      aigs = aigsd;
    if (!given("AIGD"))
      aigd = aigsd;
  }

  if (given("BIGSD"))
  {
    bigs = bigd = bigsd;
  }
  else
  {
    bigsd = (dtype == CONSTNMOS) ? 1.71E-3 : 7.59E-4;
    if (!given("BIGS"))
      bigs = bigsd;
    if (!given("BIGD"))
      bigd = bigsd;
  }
  if (given("CIGSD"))
  {
    cigs=cigd=cigsd;
  }
  else
  {
    cigsd = (dtype == CONSTNMOS) ? 0.075 : 0.03;
    if (!given("CIGS"))
      cigs=cigsd;
    if (!given("CIGD"))
      cigd=cigsd;
  }
  if (!given("IJTHDFWD"))
    ijthdfwd = ijthsfwd;
  if (!given("IJTHDREV"))
    ijthdrev = ijthsrev;
  if (!given("XJBVD"))
    xjbvd = xjbvs;
  if (!given("BVD"))
    bvd = bvs;
  if (!given("CKAPPAD"))
    ckappad = ckappas;
  if (!given("DMCI"))
    dmci = dmcg;
  if (!given("LLC"))
    Llc = Ll;
  if (!given("LWC"))
    Lwc = Lw;
  if (!given("LWLC"))
    Lwlc = Lwl;
  if (!given("WLC"))
    Wlc = Wl;
  if (!given("WWC"))
    Wwc = Ww;
  if (!given("WWLC"))
    Wwlc = Wwl;
  if (!given("DWC"))
    dwc = Wint;
  if (!given("DLC"))
    dlc = Lint;

  if (!given("AGISL"))
  {
    if (given("AGIDL"))
      agisl=agidl;
    // otherwise our initializer took care of it already
  }
  if (!given("BGISL"))
  {
    if (given("BGIDL"))
      bgisl=bgidl;
  }
  if (!given("CGISL"))
  {
    if (given("CGIDL"))
      cgisl=cgidl;
  }
  if (!given("EGISL"))
  {
    if (given("EGIDL"))
      egisl=egidl;
  }

  if (!given("DLCIG"))
    dlcig = Lint;
  if (!given("DLCIGD"))
  {
    if (!given("DLCIG"))
      dlcigd = Lint;
    else
      dlcigd = dlcig;
  }
  if (!given("DWJ"))
    dwj = dwc;




  // The seemingly unnecessary and absurd parentheses around CONSTEPS0/M_PI
  // turned out to be necessary for getting this statement to evaluate properly
  // on Windows with Intel compilers.  Bleah.
  if (!given("CF"))
    cf = 2.0 * epsrox * ((CONSTEPS0) / (M_PI)) * log(1.0 + 0.4E-6 / toxe);

  if (!given("JSD"))
    DjctSatCurDensity=SjctSatCurDensity;
  if (!given("JSWD"))
    DjctSidewallSatCurDensity=SjctSidewallSatCurDensity;
  if (!given("JSWGD"))
    DjctGateSidewallSatCurDensity=SjctGateSidewallSatCurDensity;
  if (!given("PBD"))
    DbulkJctPotential=SbulkJctPotential;
  if (!given("NJD"))
    DjctEmissionCoeff = SjctEmissionCoeff;
  if (!given("XTID"))
    DjctTempExponent = SjctTempExponent;
  if (!given("MJD"))
    DbulkJctBotGradingCoeff = SbulkJctBotGradingCoeff;
  if (!given("MJSWD"))
    DbulkJctSideGradingCoeff = SbulkJctSideGradingCoeff;
  if (!given("MJSWGS"))
    SbulkJctGateSideGradingCoeff = SbulkJctSideGradingCoeff;
  if (!given("MJSWGD"))
    DbulkJctGateSideGradingCoeff = SbulkJctGateSideGradingCoeff;
  if (!given("PBSWD"))
    DsidewallJctPotential = SsidewallJctPotential;
  if (!given("PBSWGS"))
    SGatesidewallJctPotential = SsidewallJctPotential;
  if (!given("PBSWGD"))
    DGatesidewallJctPotential = DsidewallJctPotential;
  if (!given("CJD"))
    DunitAreaJctCap=SunitAreaJctCap;
  if (!given("CJSWD"))
    DunitLengthSidewallJctCap = SunitLengthSidewallJctCap;
  if (!given("CJSWGS"))
    SunitLengthGateSidewallJctCap = SunitLengthSidewallJctCap ;
  if (!given("CJSWGD"))
    DunitLengthGateSidewallJctCap = SunitLengthGateSidewallJctCap;

  if (!given("JTSD"))
    jtsd = jtss;
  if (!given("JTSSWD"))
    jtsswd = jtssws;
  if (!given("JTSSWGD"))
    jtsswgd = jtsswgs;

  if (!given("NJTSD"))
  {
    if (given("NJTS"))
      njtsd =  njts;
  }
  if (!given("NJTSSWD"))
  {
    if (given("NJTSSW"))
      njtsswd =  njtssw;
  }
  if (!given("NJTSSWGD"))
  {
    if (given("NJTSSWG"))
      njtsswgd =  njtsswg;
  }

  if (!given("XTSD"))
    xtsd = xtss;
  if (!given("XTSSWD"))
    xtsswd = xtssws;
  if (!given("XTSSWGD"))
    xtsswgd = xtsswgs;

  if (!given("TNJTSD"))
  {
    if (given("TNJTS"))
      tnjtsd =  tnjts;
  }
  if (!given("TNJTSSWD"))
  {
    if (given("TNJTSSW"))
      tnjtsswd =  tnjtssw;
  }
  if (!given("TNJTSSWGD"))
  {
    if (given("TNJTSSWG"))
      tnjtsswgd =  tnjtsswg;
  }

  if (!given("VTSD"))
    vtsd = vtss;
  if (!given("VTSSWD"))
    vtsswd = vtssws;
  if (!given("VTSSWGD"))
    vtsswgd = vtsswgs;

  if (!given("LAGISL"))
  {
    if (given("LAGIDL"))
      lagisl = lagidl;
  }
  if (!given("LBGISL"))
  {
    if (given("LBGIDL"))
      lbgisl = lbgidl;
  }
  if (!given("LCGISL"))
  {
    if (given("LCGIDL"))
      lcgisl = lcgidl;
  }
  if (!given("LEGISL"))
  {
    if (given("LEGIDL"))
      legisl = legidl;
  }

  // This is ugly, ugly, ugly.
  // This stuff is all the "else" clauses from the spice code, which are
  // all like this:
  //  if (given("xxxsd") && (given(xxxs) || given(xxxd)))
  //  {
  //    set not-given things to constants
  //  }
  //  else
  //  {
  //     set the s and d versions to the sd version
  //  }
  //  But we set the constants in the ParTable stuff, so all that is left is
  // the else.  So I negated the conditional and left only the stuff in the
  // else.
  // TVR 1 Aug 07
  //
  if (!(!given("AIGSD") && (given("AIGS") || given("AIGD"))))
  {
    laigs = laigd = laigsd;
  }
  if (!(!given("BIGSD") && (given("BIGS") || given("BIGD"))))
  {
    lbigs = lbigd = lbigsd;
  }
  if (!(!given("CIGSD") && (given("CIGS") || given("CIGD"))))
  {
    lcigs = lcigd = lcigsd;
  }

  if (!given("WAGISL"))
  {
    if (given("WAGIDL"))
      wagisl = wagidl;
  }
  if (!given("WBGISL"))
  {
    if (given("WBGIDL"))
      wbgisl = wbgidl;
  }
  if (!given("WCGISL"))
  {
    if (given("WCGIDL"))
      wcgisl = wcgidl;
  }
  if (!given("WEGISL"))
  {
    if (given("WEGIDL"))
      wegisl = wegidl;
  }

  // See above, under "ugly, ugly, ugly"
  if (!(!given("AIGSD") && (given("AIGS") || given("AIGD"))))
  {
    waigs = waigd = waigsd;
  }
  if (!(!given("BIGSD") && (given("BIGS") || given("BIGD"))))
  {
    wbigs = wbigd = wbigsd;
  }
  if (!(!given("CIGSD") && (given("CIGS") || given("CIGD"))))
  {
    wcigs = wcigd = wcigsd;
  }

  if (!given("PAGISL"))
  {
    if (given("PAGIDL"))
      pagisl = pagidl;
  }
  if (!given("PBGISL"))
  {
    if (given("PBGIDL"))
      pbgisl = pbgidl;
  }
  if (!given("PCGISL"))
  {
    if (given("PCGIDL"))
      pcgisl = pcgidl;
  }
  if (!given("PEGISL"))
  {
    if (given("PEGIDL"))
      pegisl = pegidl;
  }

  // Vide supra, re "ugly"
  if (!(!given("AIGSD") && (given("AIGS") || given("AIGD"))))
  {
    paigs = paigd = paigsd;
  }
  if (!(!given("BIGSD") && (given("BIGS") || given("BIGD"))))
  {
    pbigs = pbigd = pbigsd;
  }
  if (!(!given("CIGSD") && (given("CIGS") || given("CIGD"))))
  {
    pcigs = pcigd = pcigsd;
  }

  if (!given("NOIA"))
  {
    if (dtype == CONSTNMOS)
      oxideTrapDensityA = 6.25e41;
    else
      oxideTrapDensityA= 6.188e40;
  }
  if (!given("NOIB"))
  {
    if (dtype == CONSTNMOS)
      oxideTrapDensityB = 3.125e26;
    else
      oxideTrapDensityB = 1.5e25;
  }
  if (!given("NOIC"))
  {
     oxideTrapDensityC = 8.75e9;
  }

  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
Model::~Model ()
{
  std::list<SizeDependParam*>::iterator it_dpL =
   sizeDependParamList.begin();
  std::list<SizeDependParam*>::iterator end_dpL =
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
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << std::endl;
  os << "    name     model name  Parameters" << std::endl;

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
// Creator        : Eric R. Keiter, SNL
// Creation Date  : 11/25/06
//----------------------------------------------------------------------------
bool Model::clearTemperatureData ()
{
  std::list<SizeDependParam*>::iterator it_dpL =
   sizeDependParamList.begin();
  std::list<SizeDependParam*>::iterator end_dpL =
   sizeDependParamList.end();
  for( ; it_dpL != end_dpL; ++it_dpL )
    delete (*it_dpL);

  sizeDependParamList.clear ();

  return true;
}

//-----------------------------------------------------------------------------
// MOSFET_B4 Master functions:
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

    bool btmp = mi.updateIntermediateVars ();
    bsuccess = bsuccess && btmp;

    // voltage drops:
    double * stoVec = mi.extData.nextStoVectorRawPtr;
    stoVec[mi.li_store_vbd]  = mi.vbd;
    stoVec[mi.li_store_vbs]  = mi.vbs;
    stoVec[mi.li_store_vgs]  = mi.vgs;
    stoVec[mi.li_store_vds]  = mi.vds;
    stoVec[mi.li_store_vges] = mi.vges;
    stoVec[mi.li_store_vgms] = mi.vgms;
    stoVec[mi.li_store_vdes] = mi.vdes;
    stoVec[mi.li_store_vses] = mi.vses;
    stoVec[mi.li_store_vdbs] = mi.vdbs;
    stoVec[mi.li_store_vsbs] = mi.vsbs;
    stoVec[mi.li_store_vdbd] = mi.vdbd;

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
    // Note the wierdness --- we have a "qg", "qb" and "qd" state variable,
    // but no corresponding instance variable --- they are all calculated from
    // other quantities that are NOT stored in the instance.  We use them
    // only in their derivative forms, cqg, cqb, and cqd.
    mi.qg = staVec[mi.li_state_qg   ] = mi.qgate;
    mi.qd = staVec[mi.li_state_qd   ] = mi.qdrn-mi.qbd;

    if (!mi.rbodyMod)
    {
      mi.qb =staVec[mi.li_state_qb   ] = mi.qbulk+mi.qbd+mi.qbs;
    }
    else
    {
      mi.qb = staVec[mi.li_state_qb   ] = mi.qbulk;
    }

    if (mi.rgateMod == 3)
    {
      staVec[mi.li_state_qgmid] = mi.qgmid;
    }

    // parasitic capacitors:
    if (mi.rbodyMod)
    {
      staVec[mi.li_state_qbs] = mi.qbs;
      staVec[mi.li_state_qbd] = mi.qbd;
    }

    if( mi.trnqsMod )
    {
      staVec[mi.li_state_qcheq] = mi.qcheq;
      staVec[mi.li_state_qcdump] = mi.qdef * mi.ScalingFactor;
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
      currStaVec[mi.li_state_qg   ] = mi.qgate;
      currStaVec[mi.li_state_qd   ] = mi.qdrn-mi.qbd;
      if (!mi.rbodyMod)
      {
        currStaVec[mi.li_state_qb   ] = mi.qbulk+mi.qbd+mi.qbs;
      }
      else
      {
        currStaVec[mi.li_state_qb   ] = mi.qbulk;
      }

      if (mi.rgateMod == 3)
      {
        currStaVec[mi.li_state_qgmid] = mi.qgmid;
      }

      // parasitic capacitors:
      if (mi.rbodyMod)
      {
        currStaVec[mi.li_state_qbs] = mi.qbs;
        currStaVec[mi.li_state_qbd] = mi.qbd;
      }

      if( mi.trnqsMod )
      {
        currStaVec[mi.li_state_qcheq] = mi.qcheq;
        currStaVec[mi.li_state_qcdump] = mi.qdef * mi.ScalingFactor;
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
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, double * leadF, double * leadQ, double * junctionV)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & mi = *(*it);

    double * dFdxdVp = mi.extData.dFdxdVpVectorRawPtr;
    double * dQdxdVp = mi.extData.dQdxdVpVectorRawPtr;
    double coef(0.0);

    mi.setupFVectorVars ();

    // Loading F-vector
    fVec[mi.li_DrainPrime] += -(mi.ceqjd - mi.ceqbd - mi.ceqdrn + mi.Idtoteq)*mi.numberParallel;
    fVec[mi.li_GatePrime] -= -(-mi.ceqgcrg + mi.Igtoteq)*mi.numberParallel;

    if (mi.rgateMod == 1)
    {
      fVec[mi.li_GateExt  ] += (mi.Igate)*mi.numberParallel;
      fVec[mi.li_GatePrime] -= (mi.Igate)*mi.numberParallel;
    }
    else if (mi.rgateMod == 2)
    {
      fVec[mi.li_GateExt] += (mi.Igate + mi.ceqgcrg)*mi.numberParallel;
      fVec[mi.li_GatePrime] -= +(mi.Igate)*mi.numberParallel;
    }
    else if (mi.rgateMod == 3)
    {
      fVec[mi.li_GateExt] += (mi.Igate)*mi.numberParallel;
      fVec[mi.li_GateMid] += (mi.IgateMid - mi.Igate + mi.ceqgcrg)*mi.numberParallel;
      fVec[mi.li_GatePrime] -= -(mi.IgateMid)*mi.numberParallel;
    }

    if (!mi.rbodyMod)
    {
      fVec[mi.li_BodyPrime] += -(mi.ceqbd + mi.ceqbs - mi.ceqjd - mi.ceqjs + mi.Ibtoteq)*mi.numberParallel;
      fVec[mi.li_SourcePrime] += -(mi.ceqdrn - mi.ceqbs + mi.ceqjs + mi.Istoteq)*mi.numberParallel;
    }
    else
    {
      fVec[mi.li_DrainBody] -= -(mi.ceqjd + mi.Idbb + mi.Idbbp)*mi.numberParallel;
      fVec[mi.li_BodyPrime] += -(mi.ceqbd + mi.ceqbs + mi.Ibtoteq +
                                      mi.Idbbp + mi.Isbbp - mi.Ibpb)*mi.numberParallel;
      fVec[mi.li_Body] += - (mi.Isbb + mi.Idbb + mi.Ibpb)*mi.numberParallel;
      fVec[mi.li_SourceBody] -= -(mi.ceqjs + mi.Isbb + mi.Isbbp)*mi.numberParallel;
      fVec[mi.li_SourcePrime] += -(mi.ceqdrn - mi.ceqbs + mi.ceqjs + mi.Istoteq)*mi.numberParallel;
    }

    if (mi.getModel().rdsMod)
    {
      fVec[mi.li_Drain]  += -(-mi.ceqgdtot)*mi.numberParallel;
      fVec[mi.li_Source] +=  -(mi.ceqgstot)*mi.numberParallel;
      fVec[mi.li_DrainPrime]  +=  -(mi.ceqgdtot)*mi.numberParallel;
      fVec[mi.li_SourcePrime] += -(-mi.ceqgstot)*mi.numberParallel;
    }

    // Idrain, Isource are linear terminal resistor currents
    if (mi.drainMOSFET_B4Exists)
    {
      fVec[mi.li_Drain]  += -(-mi.Idrain)*mi.numberParallel;
      fVec[mi.li_DrainPrime]  += -(mi.Idrain)*mi.numberParallel;
    }

    if (mi.sourceMOSFET_B4Exists)
    {
      fVec[mi.li_Source] += -(-mi.Isource)*mi.numberParallel;
      fVec[mi.li_SourcePrime] += -(+mi.Isource)*mi.numberParallel;
    }

    // Initial condition support:
    if (getSolverState().dcopFlag && mi.icVBSGiven)
    {
      coef = mi.extData.nextSolVectorRawPtr[mi.li_Ibs];
      fVec[mi.li_Body] += coef;
      fVec[mi.li_Source] += -coef;
      double cVs = mi.extData.nextSolVectorRawPtr[mi.li_Source];
      double cVb = mi.extData.nextSolVectorRawPtr[mi.li_Body];
      fVec[mi.li_Ibs] += (cVb-cVs-mi.icVBS);
    }

    if (getSolverState().dcopFlag && mi.icVDSGiven)
    {
      coef = mi.extData.nextSolVectorRawPtr[mi.li_Ids];
      fVec[mi.li_Drain] += coef;
      fVec[mi.li_Source] += -coef;
      double cVs = mi.extData.nextSolVectorRawPtr[mi.li_Source];
      double cVd = mi.extData.nextSolVectorRawPtr[mi.li_Drain];
      fVec[mi.li_Ids] += (cVd-cVs-mi.icVDS);
    }

    if (getSolverState().dcopFlag && mi.icVGSGiven)
    {
      coef = mi.extData.nextSolVectorRawPtr[mi.li_Igs];
      fVec[mi.li_GateExt] += coef;
      fVec[mi.li_Source] += -coef;
      double cVs = mi.extData.nextSolVectorRawPtr[mi.li_Source];
      double cVg = mi.extData.nextSolVectorRawPtr[mi.li_GateExt];
      fVec[mi.li_Igs] += (cVg-cVs-mi.icVGS);
    }

    // limiter section
    if (getDeviceOptions().voltageLimiterFlag && !mi.origFlag)
    {
      dFdxdVp[mi.li_DrainPrime] += (mi.ceqjd_Jdxp - mi.ceqbd_Jdxp - mi.ceqdrn_Jdxp + mi.Idtoteq_Jdxp)*mi.numberParallel;
      dFdxdVp[mi.li_GatePrime] -= (- mi.ceqgcrg_Jdxp + mi.Igtoteq_Jdxp)*mi.numberParallel;

      if (mi.rgateMod == 2)
      {
        dFdxdVp[mi.li_GateExt] += (- mi.ceqgcrg_Jdxp)*mi.numberParallel;
      }
      else if (mi.rgateMod == 3)
      {
        dFdxdVp[mi.li_GateMid] += (-mi.ceqgcrg_Jdxp)*mi.numberParallel;
      }

      if (!mi.rbodyMod)
      {
        dFdxdVp[mi.li_BodyPrime] += (mi.ceqbd_Jdxp + mi.ceqbs_Jdxp - mi.ceqjd_Jdxp - mi.ceqjs_Jdxp + mi.Ibtoteq_Jdxp)*mi.numberParallel;
        dFdxdVp[mi.li_SourcePrime] += (mi.ceqdrn_Jdxp - mi.ceqbs_Jdxp + mi.ceqjs_Jdxp + mi.Istoteq_Jdxp)*mi.numberParallel;
      }
      else
      {
        dFdxdVp[mi.li_DrainBody] -= (mi.ceqjd_Jdxp)*mi.numberParallel;
        dFdxdVp[mi.li_BodyPrime] += (mi.ceqbd_Jdxp + mi.ceqbs_Jdxp + mi.Ibtoteq_Jdxp)*mi.numberParallel;
        dFdxdVp[mi.li_SourceBody] -= (mi.ceqjs_Jdxp )*mi.numberParallel;
        dFdxdVp[mi.li_SourcePrime] += (mi.ceqdrn_Jdxp - mi.ceqbs_Jdxp + mi.ceqjs_Jdxp + mi.Istoteq_Jdxp)*mi.numberParallel;
      }

      if (mi.getModel().rdsMod)
      {
        dFdxdVp[mi.li_Drain] -= (mi.ceqgdtot_Jdxp)*mi.numberParallel;
        dFdxdVp[mi.li_Source] += (mi.ceqgstot_Jdxp)*mi.numberParallel;
        dFdxdVp[mi.li_DrainPrime] += (mi.ceqgdtot_Jdxp)*mi.numberParallel;
        dFdxdVp[mi.li_SourcePrime] -= (mi.ceqgstot_Jdxp)*mi.numberParallel;
      }
    }

    // Loading Q-vector
    mi.auxChargeCalculations ();

    double Qeqqg    = 0.0;   // gate charge
    double Qeqqb    = 0.0;   // bulk charge
    double Qeqqd    = 0.0;   // drain charge
    double Qeqqgmid = 0.0;   //
    double Qeqqjs   = 0.0;   // source-junction charge
    double Qeqqjd   = 0.0;   // drain-junction charge
    double Qqdef    = 0.0;   // nqs-related charge.
    double Qqcheq   = 0.0;   // nqs-related charge.

    if (mi.getModel().dtype > 0)
    {
      Qeqqg = mi.qg;
      Qeqqd = mi.qd;
      Qeqqb = mi.qb;

      if (mi.trnqsMod)
      {
        Qqdef = mi.qdef;
        Qqcheq = mi.qcheq;
      }

      if (mi.rbodyMod)
      {
        Qeqqjs = mi.qbs;
        Qeqqjd = mi.qbd;
      }

      if (mi.rgateMod == 3)
      {
        Qeqqgmid = mi.qgmid;
      }
    }
    else
    {
      Qeqqg = -mi.qg;
      Qeqqd = -mi.qd;
      Qeqqb = -mi.qb;

      if (mi.trnqsMod)
      {
        Qqdef = -mi.qdef;
        Qqcheq = -mi.qcheq;
      }

      if (mi.rbodyMod)
      {
        Qeqqjs = -mi.qbs;
        Qeqqjd = -mi.qbd;
      }

      if (mi.rgateMod == 3)
      {
        Qeqqgmid = -mi.qgmid;
      }
    }

    // Loading q-vector:

    qVec[mi.li_DrainPrime] += -(-Qeqqd)*mi.numberParallel;

    qVec[mi.li_GatePrime] -= -(Qeqqg)*mi.numberParallel;

    if (mi.rgateMod == 3)
    {
      qVec[mi.li_GateMid] -= -(+Qeqqgmid)*mi.numberParallel;
    }

    if (!mi.rbodyMod)
    {
      qVec[mi.li_BodyPrime] += -(-Qeqqb)*mi.numberParallel;
      qVec[mi.li_SourcePrime] += -(+Qeqqg + Qeqqb + Qeqqd + Qeqqgmid)*mi.numberParallel;
    }
    else
    {
      qVec[mi.li_DrainBody] -= -(Qeqqjd)*mi.numberParallel;
      qVec[mi.li_BodyPrime] += -(-Qeqqb)*mi.numberParallel;
      qVec[mi.li_SourceBody] -= -(Qeqqjs)*mi.numberParallel;
      qVec[mi.li_SourcePrime] += -(Qeqqd + Qeqqg + Qeqqb + Qeqqjd + Qeqqjs + Qeqqgmid)*mi.numberParallel;
    }

    if (mi.trnqsMod)
    {
      qVec[mi.li_Charge] += -(Qqcheq - Qqdef)*mi.numberParallel;
    }

    // limiter section
    if (getDeviceOptions().voltageLimiterFlag && !mi.origFlag)
    {
      dQdxdVp[mi.li_DrainPrime] += (-mi.Qeqqd_Jdxp)*mi.numberParallel;
      dQdxdVp[mi.li_GatePrime] -= (mi.Qeqqg_Jdxp)*mi.numberParallel;

      if (mi.rgateMod == 3)
      {
        dQdxdVp[mi.li_GateMid] -= (+mi.Qeqqgmid_Jdxp)*mi.numberParallel;
      }

      if (!mi.rbodyMod)
      {
        dQdxdVp[mi.li_BodyPrime] += (-mi.Qeqqb_Jdxp)*mi.numberParallel;
        dQdxdVp[mi.li_SourcePrime] += (+mi.Qeqqg_Jdxp + mi.Qeqqb_Jdxp + mi.Qeqqd_Jdxp + mi.Qeqqgmid_Jdxp)*mi.numberParallel;
      }
      else
      {
        dQdxdVp[mi.li_DrainBody] -= (mi.Qeqqjd_Jdxp)*mi.numberParallel;
        dQdxdVp[mi.li_BodyPrime] += (-mi.Qeqqb_Jdxp)*mi.numberParallel;
        dQdxdVp[mi.li_SourceBody] -= (+mi.Qeqqjs_Jdxp)*mi.numberParallel;
        dQdxdVp[mi.li_SourcePrime] += (+mi.Qeqqd_Jdxp + mi.Qeqqg_Jdxp + mi.Qeqqb_Jdxp + mi.Qeqqjd_Jdxp + mi.Qeqqjs_Jdxp + mi.Qeqqgmid_Jdxp)*mi.numberParallel;
      }

      if (mi.trnqsMod)
      {
        dQdxdVp[mi.li_Charge] += (mi.Qqcheq_Jdxp)*mi.numberParallel;
      }
    }
    
    if( mi.loadLeadCurrent )
    {	
      leadQ[mi.li_branch_dev_id] = (Qeqqd)*mi.numberParallel; 
      leadQ[mi.li_branch_dev_ig] = (Qeqqg)*mi.numberParallel; 
      leadQ[mi.li_branch_dev_ib] = (Qeqqb)*mi.numberParallel;
    
      if (!mi.rbodyMod)
      {
        leadQ[mi.li_branch_dev_is] = -(+Qeqqg + Qeqqb + Qeqqd + Qeqqgmid)*mi.numberParallel;
      }
      else
      {
        leadQ[mi.li_branch_dev_is] = -(Qeqqd + Qeqqg + Qeqqb + Qeqqjd + Qeqqjs + Qeqqgmid)*mi.numberParallel;
      }
      
      leadF[mi.li_branch_dev_id] = -(mi.ceqjd - mi.ceqbd - mi.ceqdrn + mi.Idtoteq)*mi.numberParallel;
      leadF[mi.li_branch_dev_is] = mi.Isource*mi.numberParallel;
      leadF[mi.li_branch_dev_ig] = (-mi.ceqgcrg + mi.Igtoteq)*mi.numberParallel;
      leadF[mi.li_branch_dev_ib] = 0;
    
      if (mi.rgateMod == 1)
      {
        leadF[mi.li_branch_dev_ig] += (mi.Igate)*mi.numberParallel;
      }
      else if (mi.rgateMod == 2)
      {
        leadF[mi.li_branch_dev_ig] += (mi.Igate)*mi.numberParallel;
      }
      else if (mi.rgateMod == 3)
      {
        leadF[mi.li_branch_dev_ig] += (mi.IgateMid)*mi.numberParallel;
      }
    
      if (!mi.rbodyMod)
      {
        leadF[mi.li_branch_dev_ib] += -(mi.ceqbd + mi.ceqbs - mi.ceqjd - mi.ceqjs + mi.Ibtoteq)*mi.numberParallel;
        leadF[mi.li_branch_dev_is] += -(mi.ceqdrn - mi.ceqbs + mi.ceqjs + mi.Istoteq)*mi.numberParallel;
      }
      else
      {
        leadF[mi.li_branch_dev_ib] = - (mi.Isbb + mi.Idbb + mi.Ibpb)*mi.numberParallel;
        leadF[mi.li_branch_dev_is] += -(mi.ceqdrn - mi.ceqbs + mi.ceqjs + mi.Istoteq)*mi.numberParallel;
      }
    
      if (mi.model_.rdsMod)
      {
        leadF[mi.li_branch_dev_id]  += (mi.ceqgdtot)*mi.numberParallel;
        leadF[mi.li_branch_dev_is]  +=  -(mi.ceqgstot)*mi.numberParallel;
      }
      
      junctionV[mi.li_branch_dev_id] = solVec[mi.li_Drain] - solVec[mi.li_Source];
      junctionV[mi.li_branch_dev_ig] = solVec[mi.li_GateExt] - solVec[mi.li_Source];
      junctionV[mi.li_branch_dev_is] = 0.0;
      junctionV[mi.li_branch_dev_ib] = 0.0 ; 
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
    Instance & mi = *(*it);

    // F-matrix:
    if (!mi.rbodyMod)
    {
      mi.gjbd = mi.gbd;
      mi.gjbs = mi.gbs;
    }
    else
    {
      mi.gjbd = mi.gjbs = 0.0;
    }

    if (!mi.getModel().rdsMod)
    {
      mi.gdpr = mi.drainConductance;
      mi.gspr = mi.sourceConductance;
    }
    else
    {
      mi.gdpr = mi.gspr = 0.0;
    }

    mi.geltd = mi.grgeltd;

    double T1 = mi.qdef * mi.gtau;

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    if (mi.rgateMod == 1)
    {
      *mi.f_GEgePtr += (mi.geltd)*mi.numberParallel;
      *mi.f_GEgpPtr -= (mi.geltd)*mi.numberParallel;
      *mi.f_GPgePtr -= (mi.geltd)*mi.numberParallel;
      *mi.f_GPgpPtr += (+ mi.geltd - mi.ggtg + mi.gIgtotg)*mi.numberParallel;
      *mi.f_GPdpPtr += (- mi.ggtd + mi.gIgtotd)*mi.numberParallel;
      *mi.f_GPspPtr += (- mi.ggts + mi.gIgtots)*mi.numberParallel;
      *mi.f_GPbpPtr += (- mi.ggtb + mi.gIgtotb)*mi.numberParallel;
    } // WDLiu: gcrg already subtracted from all gcrgg below
    else if (mi.rgateMod == 2)
    {
      *mi.f_GEgePtr += (mi.gcrg)*mi.numberParallel;
      *mi.f_GEgpPtr += (mi.gcrgg)*mi.numberParallel;
      *mi.f_GEdpPtr += (mi.gcrgd)*mi.numberParallel;
      *mi.f_GEspPtr += (mi.gcrgs)*mi.numberParallel;
      *mi.f_GEbpPtr += (mi.gcrgb)*mi.numberParallel;

      *mi.f_GPgePtr -= (mi.gcrg)*mi.numberParallel;
      *mi.f_GPgpPtr += (- mi.gcrgg - mi.ggtg + mi.gIgtotg)*mi.numberParallel;
      *mi.f_GPdpPtr += (- mi.gcrgd - mi.ggtd + mi.gIgtotd)*mi.numberParallel;
      *mi.f_GPspPtr += (- mi.gcrgs - mi.ggts + mi.gIgtots)*mi.numberParallel;
      *mi.f_GPbpPtr += (- mi.gcrgb - mi.ggtb + mi.gIgtotb)*mi.numberParallel;
    }
    else if (mi.rgateMod == 3)
    {
      *mi.f_GEgePtr += (mi.geltd)*mi.numberParallel;
      *mi.f_GEgmPtr -= (mi.geltd)*mi.numberParallel;
      *mi.f_GMgePtr -= (mi.geltd)*mi.numberParallel;
      *mi.f_GMgmPtr += (mi.geltd + mi.gcrg)*mi.numberParallel;

      *mi.f_GMdpPtr += (mi.gcrgd)*mi.numberParallel;
      *mi.f_GMgpPtr += (mi.gcrgg)*mi.numberParallel;
      *mi.f_GMspPtr += (mi.gcrgs)*mi.numberParallel;
      *mi.f_GMbpPtr += (mi.gcrgb)*mi.numberParallel;

      *mi.f_GPgmPtr -= (mi.gcrg)*mi.numberParallel;

      *mi.f_GPgpPtr += (- mi.gcrgg - mi.ggtg + mi.gIgtotg)*mi.numberParallel;
      *mi.f_GPdpPtr += (- mi.gcrgd - mi.ggtd + mi.gIgtotd)*mi.numberParallel;
      *mi.f_GPspPtr += (- mi.gcrgs - mi.ggts + mi.gIgtots)*mi.numberParallel;
      *mi.f_GPbpPtr += (- mi.gcrgb - mi.ggtb + mi.gIgtotb)*mi.numberParallel;
    }
    else
    {
      *mi.f_GPgpPtr += (- mi.ggtg + mi.gIgtotg)*mi.numberParallel;
      *mi.f_GPdpPtr += (- mi.ggtd + mi.gIgtotd)*mi.numberParallel;
      *mi.f_GPspPtr += (- mi.ggts + mi.gIgtots)*mi.numberParallel;
      *mi.f_GPbpPtr += (- mi.ggtb + mi.gIgtotb)*mi.numberParallel;
    }

    if (mi.getModel().rdsMod)
    {
      *mi.f_DgpPtr += (mi.gdtotg)*mi.numberParallel;
      *mi.f_DspPtr += (mi.gdtots)*mi.numberParallel;
      *mi.f_DbpPtr += (mi.gdtotb)*mi.numberParallel;
      *mi.f_SdpPtr += (mi.gstotd)*mi.numberParallel;
      *mi.f_SgpPtr += (mi.gstotg)*mi.numberParallel;
      *mi.f_SbpPtr += (mi.gstotb)*mi.numberParallel;
    }


    *mi.f_DPdpPtr += (mi.gdpr + mi.gds + mi.gbd + T1 * mi.ddxpart_dVd
              - mi.gdtotd + mi.RevSum + mi.gbdpdp + mi.dxpart * mi.ggtd - mi.gIdtotd)*mi.numberParallel;


    *mi.f_DPdPtr -= (mi.gdpr + mi.gdtot)*mi.numberParallel;

    *mi.f_DPgpPtr += (mi.Gm - mi.gdtotg + mi.gbdpg - mi.gIdtotg
              + mi.dxpart * mi.ggtg + T1 * mi.ddxpart_dVg)*mi.numberParallel;


    *mi.f_DPspPtr -= (mi.gds + mi.gdtots - mi.dxpart * mi.ggts + mi.gIdtots
              - T1 * mi.ddxpart_dVs + mi.FwdSum - mi.gbdpsp)*mi.numberParallel;


    *mi.f_DPbpPtr -= (mi.gjbd + mi.gdtotb - mi.Gmbs - mi.gbdpb + mi.gIdtotb
              - T1 * mi.ddxpart_dVb - mi.dxpart * mi.ggtb)*mi.numberParallel;


    *mi.f_DdpPtr -= (mi.gdpr - mi.gdtotd)*mi.numberParallel;

    *mi.f_DdPtr += (mi.gdpr + mi.gdtot)*mi.numberParallel;


    *mi.f_SPdpPtr -= (mi.gds + mi.gstotd + mi.RevSum - mi.gbspdp
                          - T1 * mi.dsxpart_dVd - mi.sxpart * mi.ggtd + mi.gIstotd)*mi.numberParallel;


    *mi.f_SPgpPtr += (- mi.Gm - mi.gstotg + mi.gbspg + mi.sxpart * mi.ggtg
                        + T1 * mi.dsxpart_dVg - mi.gIstotg)*mi.numberParallel;


    *mi.f_SPspPtr += (mi.gspr + mi.gds + mi.gbs + T1 * mi.dsxpart_dVs
                          - mi.gstots + mi.FwdSum + mi.gbspsp + mi.sxpart * mi.ggts - mi.gIstots)*mi.numberParallel;


    *mi.f_SPsPtr -= (mi.gspr + mi.gstot)*mi.numberParallel;


    *mi.f_SPbpPtr -= (mi.gjbs + mi.gstotb + mi.Gmbs - mi.gbspb - mi.sxpart * mi.ggtb
                          - T1 * mi.dsxpart_dVb + mi.gIstotb)*mi.numberParallel;


    *mi.f_SspPtr -= (mi.gspr - mi.gstots)*mi.numberParallel;

    *mi.f_SsPtr += (mi.gspr + mi.gstot)*mi.numberParallel;


    *mi.f_BPdpPtr += (- mi.gjbd + mi.gbbdp - mi.gIbtotd)*mi.numberParallel;

    *mi.f_BPgpPtr += (- mi.gbgs - mi.gIbtotg)*mi.numberParallel;

    *mi.f_BPspPtr += (- mi.gjbs + mi.gbbsp - mi.gIbtots)*mi.numberParallel;

    *mi.f_BPbpPtr += (mi.gjbd + mi.gjbs - mi.gbbs - mi.gIbtotb)*mi.numberParallel;

    //ggidld = (ggidld)*mi.numberParallel;
    //ggidlg = (ggidlg)*mi.numberParallel;
    //ggidlb = (ggidlb)*mi.numberParallel;
    //ggislg = (ggislg)*mi.numberParallel;
    //ggisls = (ggisls)*mi.numberParallel;
    //ggislb = (ggislb)*mi.numberParallel;

    // stamp gidl

    *mi.f_DPdpPtr += (mi.ggidld)*mi.numberParallel;

    *mi.f_DPgpPtr += (mi.ggidlg)*mi.numberParallel;

    *mi.f_DPspPtr -= ((mi.ggidlg + mi.ggidld + mi.ggidlb))*mi.numberParallel;

    *mi.f_DPbpPtr += (mi.ggidlb)*mi.numberParallel;

    *mi.f_BPdpPtr -= (mi.ggidld)*mi.numberParallel;

    *mi.f_BPgpPtr -= (mi.ggidlg)*mi.numberParallel;

    *mi.f_BPspPtr += ((mi.ggidlg + mi.ggidld + mi.ggidlb))*mi.numberParallel;

    *mi.f_BPbpPtr -= (mi.ggidlb)*mi.numberParallel;
    // stamp gisl

    *mi.f_SPdpPtr -= ((mi.ggisls + mi.ggislg + mi.ggislb))*mi.numberParallel;

    *mi.f_SPgpPtr += (mi.ggislg)*mi.numberParallel;

    *mi.f_SPspPtr += (mi.ggisls)*mi.numberParallel;

    *mi.f_SPbpPtr += (mi.ggislb)*mi.numberParallel;

    *mi.f_BPdpPtr += ((mi.ggislg + mi.ggisls + mi.ggislb))*mi.numberParallel;

    *mi.f_BPgpPtr -= (mi.ggislg)*mi.numberParallel;

    *mi.f_BPspPtr -= (mi.ggisls)*mi.numberParallel;

    *mi.f_BPbpPtr -= (mi.ggislb)*mi.numberParallel;


    if (mi.rbodyMod)
    {
      *mi.f_DPdbPtr += (- mi.gbd)*mi.numberParallel;
      *mi.f_SPsbPtr -= (mi.gbs)*mi.numberParallel;

      *mi.f_DBdpPtr += (- mi.gbd)*mi.numberParallel;
      *mi.f_DBdbPtr += (mi.gbd + mi.grbpd + mi.grbdb)*mi.numberParallel;
      *mi.f_DBbpPtr -= (mi.grbpd)*mi.numberParallel;
      *mi.f_DBbPtr -= (mi.grbdb)*mi.numberParallel;

      *mi.f_BPdbPtr -= (mi.grbpd)*mi.numberParallel;
      *mi.f_BPbPtr -= (mi.grbpb)*mi.numberParallel;
      *mi.f_BPsbPtr -= (mi.grbps)*mi.numberParallel;
      *mi.f_BPbpPtr += (mi.grbpd + mi.grbps + mi.grbpb)*mi.numberParallel;
      // WDLiu: (gcbbb - gbbs) already added to mi.BPbpPtr

      *mi.f_SBspPtr += (- mi.gbs)*mi.numberParallel;
      *mi.f_SBbpPtr -= (mi.grbps)*mi.numberParallel;
      *mi.f_SBbPtr -= (mi.grbsb)*mi.numberParallel;
      *mi.f_SBsbPtr += (mi.gbs + mi.grbps + mi.grbsb)*mi.numberParallel;

      *mi.f_BdbPtr -= (mi.grbdb)*mi.numberParallel;
      *mi.f_BbpPtr -= (mi.grbpb)*mi.numberParallel;
      *mi.f_BsbPtr -= (mi.grbsb)*mi.numberParallel;
      *mi.f_BbPtr += (mi.grbsb + mi.grbdb + mi.grbpb)*mi.numberParallel;
    }

    if (mi.trnqsMod)
    {
      *mi.f_QqPtr += (mi.gqdef + mi.gtau)*mi.numberParallel;
      *mi.f_QgpPtr += (mi.ggtg)*mi.numberParallel;
      *mi.f_QdpPtr += (mi.ggtd)*mi.numberParallel;
      *mi.f_QspPtr += (mi.ggts)*mi.numberParallel;
      *mi.f_QbpPtr += (mi.ggtb)*mi.numberParallel;

      *mi.f_DPqPtr += (mi.dxpart * mi.gtau)*mi.numberParallel;
      *mi.f_SPqPtr += (mi.sxpart * mi.gtau)*mi.numberParallel;
      *mi.f_GPqPtr -= (mi.gtau)*mi.numberParallel;
    }

    // Initial Conditions:
    if (mi.icVBSGiven)
    {
      if (getSolverState().dcopFlag)
      {
        *mi.f_BibsPtr += 1.0;
        *mi.f_SibsPtr += -1.0;
        *mi.f_IBSbPtr += 1.0;
        *mi.f_IBSsPtr += -1.0;
      }
      else
      {
        *mi.f_IBSibsPtr = 1.0;
      }
    }

    if (mi.icVDSGiven)
    {
      if (getSolverState().dcopFlag)
      {
        *mi.f_DidsPtr += 1.0;
        *mi.f_SidsPtr += -1.0;
        *mi.f_IDSdPtr += 1.0;
        *mi.f_IDSsPtr += -1.0;
      }
      else
      {
        *mi.f_IDSidsPtr = 1.0;
      }
    }

    if (mi.icVGSGiven)
    {
      if (getSolverState().dcopFlag)
      {
        *mi.f_GEigsPtr += 1.0;
        *mi.f_SigsPtr += -1.0;
        *mi.f_IGSgPtr += 1.0;
        *mi.f_IGSsPtr += -1.0;
      }
      else
      {
        *mi.f_IGSigsPtr = 1.0;
      }
    }

#else
    if (mi.rgateMod == 1)
    {
      dFdx[mi.li_GateExt][mi.GEge] += (mi.geltd)*mi.numberParallel;
      dFdx[mi.li_GateExt][mi.GEgp] -= (mi.geltd)*mi.numberParallel;
      dFdx[mi.li_GatePrime][mi.GPge] -= (mi.geltd)*mi.numberParallel;
      dFdx[mi.li_GatePrime][mi.GPgp] += (+ mi.geltd - mi.ggtg + mi.gIgtotg)*mi.numberParallel;
      dFdx[mi.li_GatePrime][mi.GPdp] += (- mi.ggtd + mi.gIgtotd)*mi.numberParallel;
      dFdx[mi.li_GatePrime][mi.GPsp] += (- mi.ggts + mi.gIgtots)*mi.numberParallel;
      dFdx[mi.li_GatePrime][mi.GPbp] += (- mi.ggtb + mi.gIgtotb)*mi.numberParallel;
    } // WDLiu: gcrg already subtracted from all gcrgg below
    else if (mi.rgateMod == 2)
    {
      dFdx[mi.li_GateExt][mi.GEge] += (mi.gcrg)*mi.numberParallel;
      dFdx[mi.li_GateExt][mi.GEgp] += (mi.gcrgg)*mi.numberParallel;
      dFdx[mi.li_GateExt][mi.GEdp] += (mi.gcrgd)*mi.numberParallel;
      dFdx[mi.li_GateExt][mi.GEsp] += (mi.gcrgs)*mi.numberParallel;
      dFdx[mi.li_GateExt][mi.GEbp] += (mi.gcrgb)*mi.numberParallel;

      dFdx[mi.li_GatePrime][mi.GPge] -= (mi.gcrg)*mi.numberParallel;
      dFdx[mi.li_GatePrime][mi.GPgp] += (- mi.gcrgg - mi.ggtg + mi.gIgtotg)*mi.numberParallel;
      dFdx[mi.li_GatePrime][mi.GPdp] += (- mi.gcrgd - mi.ggtd + mi.gIgtotd)*mi.numberParallel;
      dFdx[mi.li_GatePrime][mi.GPsp] += (- mi.gcrgs - mi.ggts + mi.gIgtots)*mi.numberParallel;
      dFdx[mi.li_GatePrime][mi.GPbp] += (- mi.gcrgb - mi.ggtb + mi.gIgtotb)*mi.numberParallel;
    }
    else if (mi.rgateMod == 3)
    {
      dFdx[mi.li_GateExt][mi.GEge] += (mi.geltd)*mi.numberParallel;
      dFdx[mi.li_GateExt][mi.GEgm] -= (mi.geltd)*mi.numberParallel;
      dFdx[mi.li_GateMid][mi.GMge] -= (mi.geltd)*mi.numberParallel;
      dFdx[mi.li_GateMid][mi.GMgm] += (mi.geltd + mi.gcrg)*mi.numberParallel;

      dFdx[mi.li_GateMid][mi.GMdp] += (mi.gcrgd)*mi.numberParallel;
      dFdx[mi.li_GateMid][mi.GMgp] += (mi.gcrgg)*mi.numberParallel;
      dFdx[mi.li_GateMid][mi.GMsp] += (mi.gcrgs)*mi.numberParallel;
      dFdx[mi.li_GateMid][mi.GMbp] += (mi.gcrgb)*mi.numberParallel;

      dFdx[mi.li_GatePrime][mi.GPgm] -= (mi.gcrg)*mi.numberParallel;

      dFdx[mi.li_GatePrime][mi.GPgp] += (- mi.gcrgg - mi.ggtg + mi.gIgtotg)*mi.numberParallel;
      dFdx[mi.li_GatePrime][mi.GPdp] += (- mi.gcrgd - mi.ggtd + mi.gIgtotd)*mi.numberParallel;
      dFdx[mi.li_GatePrime][mi.GPsp] += (- mi.gcrgs - mi.ggts + mi.gIgtots)*mi.numberParallel;
      dFdx[mi.li_GatePrime][mi.GPbp] += (- mi.gcrgb - mi.ggtb + mi.gIgtotb)*mi.numberParallel;
    }
    else
    {
      dFdx[mi.li_GatePrime][mi.GPgp] += (- mi.ggtg + mi.gIgtotg)*mi.numberParallel;
      dFdx[mi.li_GatePrime][mi.GPdp] += (- mi.ggtd + mi.gIgtotd)*mi.numberParallel;
      dFdx[mi.li_GatePrime][mi.GPsp] += (- mi.ggts + mi.gIgtots)*mi.numberParallel;
      dFdx[mi.li_GatePrime][mi.GPbp] += (- mi.ggtb + mi.gIgtotb)*mi.numberParallel;
    }

    if (mi.getModel().rdsMod)
    {
      dFdx[mi.li_Drain][mi.Dgp] += (mi.gdtotg)*mi.numberParallel;
      dFdx[mi.li_Drain][mi.Dsp] += (mi.gdtots)*mi.numberParallel;
      dFdx[mi.li_Drain][mi.Dbp] += (mi.gdtotb)*mi.numberParallel;
      dFdx[mi.li_Source][mi.Sdp] += (mi.gstotd)*mi.numberParallel;
      dFdx[mi.li_Source][mi.Sgp] += (mi.gstotg)*mi.numberParallel;
      dFdx[mi.li_Source][mi.Sbp] += (mi.gstotb)*mi.numberParallel;
    }


    dFdx[mi.li_DrainPrime][mi.DPdp] += (mi.gdpr + mi.gds + mi.gbd + T1 * mi.ddxpart_dVd
              - mi.gdtotd + mi.RevSum + mi.gbdpdp + mi.dxpart * mi.ggtd - mi.gIdtotd)*mi.numberParallel;


    dFdx[mi.li_DrainPrime][mi.DPd] -= (mi.gdpr + mi.gdtot)*mi.numberParallel;

    dFdx[mi.li_DrainPrime][mi.DPgp] += (mi.Gm - mi.gdtotg + mi.gbdpg - mi.gIdtotg
              + mi.dxpart * mi.ggtg + T1 * mi.ddxpart_dVg)*mi.numberParallel;


    dFdx[mi.li_DrainPrime][mi.DPsp] -= (mi.gds + mi.gdtots - mi.dxpart * mi.ggts + mi.gIdtots
              - T1 * mi.ddxpart_dVs + mi.FwdSum - mi.gbdpsp)*mi.numberParallel;


    dFdx[mi.li_DrainPrime][mi.DPbp] -= (mi.gjbd + mi.gdtotb - mi.Gmbs - mi.gbdpb + mi.gIdtotb
              - T1 * mi.ddxpart_dVb - mi.dxpart * mi.ggtb)*mi.numberParallel;


    dFdx[mi.li_Drain][mi.Ddp] -= (mi.gdpr - mi.gdtotd)*mi.numberParallel;

    dFdx[mi.li_Drain][mi.Dd] += (mi.gdpr + mi.gdtot)*mi.numberParallel;


    dFdx[mi.li_SourcePrime][mi.SPdp] -= (mi.gds + mi.gstotd + mi.RevSum - mi.gbspdp
                          - T1 * mi.dsxpart_dVd - mi.sxpart * mi.ggtd + mi.gIstotd)*mi.numberParallel;


    dFdx[mi.li_SourcePrime][mi.SPgp] += (- mi.Gm - mi.gstotg + mi.gbspg + mi.sxpart * mi.ggtg
                        + T1 * mi.dsxpart_dVg - mi.gIstotg)*mi.numberParallel;


    dFdx[mi.li_SourcePrime][mi.SPsp] += (mi.gspr + mi.gds + mi.gbs + T1 * mi.dsxpart_dVs
                          - mi.gstots + mi.FwdSum + mi.gbspsp + mi.sxpart * mi.ggts - mi.gIstots)*mi.numberParallel;


    dFdx[mi.li_SourcePrime][mi.SPs] -= (mi.gspr + mi.gstot)*mi.numberParallel;


    dFdx[mi.li_SourcePrime][mi.SPbp] -= (mi.gjbs + mi.gstotb + mi.Gmbs - mi.gbspb - mi.sxpart * mi.ggtb
                          - T1 * mi.dsxpart_dVb + mi.gIstotb)*mi.numberParallel;


    dFdx[mi.li_Source][mi.Ssp] -= (mi.gspr - mi.gstots)*mi.numberParallel;

    dFdx[mi.li_Source][mi.Ss] += (mi.gspr + mi.gstot)*mi.numberParallel;


    dFdx[mi.li_BodyPrime][mi.BPdp] += (- mi.gjbd + mi.gbbdp - mi.gIbtotd)*mi.numberParallel;

    dFdx[mi.li_BodyPrime][mi.BPgp] += (- mi.gbgs - mi.gIbtotg)*mi.numberParallel;

    dFdx[mi.li_BodyPrime][mi.BPsp] += (- mi.gjbs + mi.gbbsp - mi.gIbtots)*mi.numberParallel;

    dFdx[mi.li_BodyPrime][mi.BPbp] += (mi.gjbd + mi.gjbs - mi.gbbs - mi.gIbtotb)*mi.numberParallel;

    //ggidld = (ggidld)*mi.numberParallel;
    //ggidlg = (ggidlg)*mi.numberParallel;
    //ggidlb = (ggidlb)*mi.numberParallel;
    //ggislg = (ggislg)*mi.numberParallel;
    //ggisls = (ggisls)*mi.numberParallel;
    //ggislb = (ggislb)*mi.numberParallel;

    // stamp gidl

    dFdx[mi.li_DrainPrime][mi.DPdp] += (mi.ggidld)*mi.numberParallel;

    dFdx[mi.li_DrainPrime][mi.DPgp] += (mi.ggidlg)*mi.numberParallel;

    dFdx[mi.li_DrainPrime][mi.DPsp] -= ((mi.ggidlg + mi.ggidld + mi.ggidlb))*mi.numberParallel;

    dFdx[mi.li_DrainPrime][mi.DPbp] += (mi.ggidlb)*mi.numberParallel;

    dFdx[mi.li_BodyPrime][mi.BPdp] -= (mi.ggidld)*mi.numberParallel;

    dFdx[mi.li_BodyPrime][mi.BPgp] -= (mi.ggidlg)*mi.numberParallel;

    dFdx[mi.li_BodyPrime][mi.BPsp] += ((mi.ggidlg + mi.ggidld + mi.ggidlb))*mi.numberParallel;

    dFdx[mi.li_BodyPrime][mi.BPbp] -= (mi.ggidlb)*mi.numberParallel;
    // stamp gisl

    dFdx[mi.li_SourcePrime][mi.SPdp] -= ((mi.ggisls + mi.ggislg + mi.ggislb))*mi.numberParallel;

    dFdx[mi.li_SourcePrime][mi.SPgp] += (mi.ggislg)*mi.numberParallel;

    dFdx[mi.li_SourcePrime][mi.SPsp] += (mi.ggisls)*mi.numberParallel;

    dFdx[mi.li_SourcePrime][mi.SPbp] += (mi.ggislb)*mi.numberParallel;

    dFdx[mi.li_BodyPrime][mi.BPdp] += ((mi.ggislg + mi.ggisls + mi.ggislb))*mi.numberParallel;

    dFdx[mi.li_BodyPrime][mi.BPgp] -= (mi.ggislg)*mi.numberParallel;

    dFdx[mi.li_BodyPrime][mi.BPsp] -= (mi.ggisls)*mi.numberParallel;

    dFdx[mi.li_BodyPrime][mi.BPbp] -= (mi.ggislb)*mi.numberParallel;


    if (mi.rbodyMod)
    {
      dFdx[mi.li_DrainPrime][mi.DPdb] += (- mi.gbd)*mi.numberParallel;
      dFdx[mi.li_SourcePrime][mi.SPsb] -= (mi.gbs)*mi.numberParallel;

      dFdx[mi.li_DrainBody][mi.DBdp] += (- mi.gbd)*mi.numberParallel;
      dFdx[mi.li_DrainBody][mi.DBdb] += (mi.gbd + mi.grbpd + mi.grbdb)*mi.numberParallel;
      dFdx[mi.li_DrainBody][mi.DBbp] -= (mi.grbpd)*mi.numberParallel;
      dFdx[mi.li_DrainBody][mi.DBb] -= (mi.grbdb)*mi.numberParallel;

      dFdx[mi.li_BodyPrime][mi.BPdb] -= (mi.grbpd)*mi.numberParallel;
      dFdx[mi.li_BodyPrime][mi.BPb] -= (mi.grbpb)*mi.numberParallel;
      dFdx[mi.li_BodyPrime][mi.BPsb] -= (mi.grbps)*mi.numberParallel;
      dFdx[mi.li_BodyPrime][mi.BPbp] += (mi.grbpd + mi.grbps + mi.grbpb)*mi.numberParallel;
      // WDLiu: (gcbbb - gbbs) already added to mi.BPbpPtr

      dFdx[mi.li_SourceBody][mi.SBsp] += (- mi.gbs)*mi.numberParallel;
      dFdx[mi.li_SourceBody][mi.SBbp] -= (mi.grbps)*mi.numberParallel;
      dFdx[mi.li_SourceBody][mi.SBb] -= (mi.grbsb)*mi.numberParallel;
      dFdx[mi.li_SourceBody][mi.SBsb] += (mi.gbs + mi.grbps + mi.grbsb)*mi.numberParallel;

      dFdx[mi.li_Body][mi.Bdb] -= (mi.grbdb)*mi.numberParallel;
      dFdx[mi.li_Body][mi.Bbp] -= (mi.grbpb)*mi.numberParallel;
      dFdx[mi.li_Body][mi.Bsb] -= (mi.grbsb)*mi.numberParallel;
      dFdx[mi.li_Body][mi.Bb] += (mi.grbsb + mi.grbdb + mi.grbpb)*mi.numberParallel;
    }

    if (mi.trnqsMod)
    {
      dFdx[mi.li_Charge][mi.Qq] += (mi.gqdef + mi.gtau)*mi.numberParallel;
      dFdx[mi.li_Charge][mi.Qgp] += (mi.ggtg)*mi.numberParallel;
      dFdx[mi.li_Charge][mi.Qdp] += (mi.ggtd)*mi.numberParallel;
      dFdx[mi.li_Charge][mi.Qsp] += (mi.ggts)*mi.numberParallel;
      dFdx[mi.li_Charge][mi.Qbp] += (mi.ggtb)*mi.numberParallel;

      dFdx[mi.li_DrainPrime][mi.DPq] += (mi.dxpart * mi.gtau)*mi.numberParallel;
      dFdx[mi.li_SourcePrime][mi.SPq] += (mi.sxpart * mi.gtau)*mi.numberParallel;
      dFdx[mi.li_GatePrime][mi.GPq] -= (mi.gtau)*mi.numberParallel;
    }

    // Initial Conditions:
    if (mi.icVBSGiven)
    {
      if (getSolverState().dcopFlag)
      {
        dFdx[mi.li_Body][mi.Bibs] += 1.0;
        dFdx[mi.li_Source][mi.Sibs] += -1.0;
        dFdx[mi.li_Ibs][mi.IBSb] += 1.0;
        dFdx[mi.li_Ibs][mi.IBSs] += -1.0;
      }
      else
      {
        dFdx[mi.li_Ibs][mi.IBSibs] = 1.0;
      }
    }

    if (mi.icVDSGiven)
    {
      if (getSolverState().dcopFlag)
      {
        dFdx[mi.li_Drain][mi.Dids] += 1.0;
        dFdx[mi.li_Source][mi.Sids] += -1.0;
        dFdx[mi.li_Ids][mi.IDSd] += 1.0;
        dFdx[mi.li_Ids][mi.IDSs] += -1.0;
      }
      else
      {
        dFdx[mi.li_Ids][mi.IDSids] = 1.0;
      }
    }

    if (mi.icVGSGiven)
    {
      if (getSolverState().dcopFlag)
      {
        dFdx[mi.li_GateExt][mi.GEigs] += 1.0;
        dFdx[mi.li_Source][mi.Sigs] += -1.0;
        dFdx[mi.li_Igs][mi.IGSg] += 1.0;
        dFdx[mi.li_Igs][mi.IGSs] += -1.0;
      }
      else
      {
        dFdx[mi.li_Igs][mi.IGSigs] = 1.0;
      }
    }

#endif
    {
      // These are only used for the nqsMod variation
      // which isn't currently implemented

#ifndef Xyce_NONPOINTER_MATRIX_LOAD

      if (mi.rgateMod == 1)
      {
        *mi.q_GPgpPtr += (mi.CAPcggb)*mi.numberParallel;
        *mi.q_GPdpPtr += (mi.CAPcgdb)*mi.numberParallel;
        *mi.q_GPspPtr += (mi.CAPcgsb)*mi.numberParallel;
        *mi.q_GPbpPtr += (mi.CAPcgbb)*mi.numberParallel;
      } // WDLiu: CAPcrg already subtracted from all CAPcrgg below
      else if (mi.rgateMod == 2)
      {
        *mi.q_GPgpPtr += (mi.CAPcggb)*mi.numberParallel;
        *mi.q_GPdpPtr += (mi.CAPcgdb)*mi.numberParallel;
        *mi.q_GPspPtr += (mi.CAPcgsb)*mi.numberParallel;
        *mi.q_GPbpPtr += (mi.CAPcgbb)*mi.numberParallel;
      }
      else if (mi.rgateMod == 3)
      {
        *mi.q_GMgmPtr += (+ mi.CAPcgmgmb)*mi.numberParallel;

        *mi.q_GMdpPtr += (mi.CAPcgmdb)*mi.numberParallel;
        *mi.q_GMspPtr += (mi.CAPcgmsb)*mi.numberParallel;
        *mi.q_GMbpPtr += (mi.CAPcgmbb)*mi.numberParallel;

        *mi.q_DPgmPtr += (mi.CAPcdgmb)*mi.numberParallel;
        *mi.q_SPgmPtr += (mi.CAPcsgmb)*mi.numberParallel;
        *mi.q_BPgmPtr += (mi.CAPcbgmb)*mi.numberParallel;

        *mi.q_GPgpPtr += (mi.CAPcggb)*mi.numberParallel;
        *mi.q_GPdpPtr += (mi.CAPcgdb)*mi.numberParallel;
        *mi.q_GPspPtr += (mi.CAPcgsb)*mi.numberParallel;
        *mi.q_GPbpPtr += (mi.CAPcgbb)*mi.numberParallel;
      }
      else
      {
        *mi.q_GPgpPtr += (mi.CAPcggb)*mi.numberParallel;
        *mi.q_GPdpPtr += (mi.CAPcgdb)*mi.numberParallel;
        *mi.q_GPspPtr += (mi.CAPcgsb)*mi.numberParallel;
        *mi.q_GPbpPtr += (mi.CAPcgbb)*mi.numberParallel;
      }

      *mi.q_DPdpPtr += (mi.CAPcddb)*mi.numberParallel;
      *mi.q_DPgpPtr += (+ mi.CAPcdgb)*mi.numberParallel;
      *mi.q_DPspPtr -= (- mi.CAPcdsb)*mi.numberParallel;
      *mi.q_DPbpPtr -= (- mi.CAPcdbb)*mi.numberParallel;

      *mi.q_SPdpPtr -= (- mi.CAPcsdb)*mi.numberParallel;
      *mi.q_SPgpPtr += (mi.CAPcsgb)*mi.numberParallel;
      *mi.q_SPspPtr += (mi.CAPcssb)*mi.numberParallel;
      *mi.q_SPbpPtr -= (- mi.CAPcsbb)*mi.numberParallel;

      *mi.q_BPdpPtr += (mi.CAPcbdb)*mi.numberParallel;
      *mi.q_BPgpPtr += (mi.CAPcbgb)*mi.numberParallel;
      *mi.q_BPspPtr += (mi.CAPcbsb)*mi.numberParallel;
      *mi.q_BPbpPtr += (mi.CAPcbbb)*mi.numberParallel;

      if (mi.rbodyMod)
      {
        *mi.q_DPdbPtr += (mi.CAPcdbdb)*mi.numberParallel;
        *mi.q_SPsbPtr -= (- mi.CAPcsbsb)*mi.numberParallel;

        *mi.q_DBdpPtr += (mi.CAPcdbdb)*mi.numberParallel;
        *mi.q_DBdbPtr += (- mi.CAPcdbdb)*mi.numberParallel;

        *mi.q_SBspPtr += (mi.CAPcsbsb)*mi.numberParallel;
        *mi.q_SBsbPtr += (- mi.CAPcsbsb)*mi.numberParallel;
      }

      if (mi.trnqsMod)
      {
        *mi.q_QgpPtr += (- mi.CAPcqgb)*mi.numberParallel;
        *mi.q_QdpPtr += (- mi.CAPcqdb)*mi.numberParallel;
        *mi.q_QspPtr += (- mi.CAPcqsb)*mi.numberParallel;
        *mi.q_QbpPtr += (- mi.CAPcqbb)*mi.numberParallel;
      }

#else
      if (mi.rgateMod == 1)
      {
        dQdx[mi.li_GatePrime][mi.GPgp] += (mi.CAPcggb)*mi.numberParallel;
        dQdx[mi.li_GatePrime][mi.GPdp] += (mi.CAPcgdb)*mi.numberParallel;
        dQdx[mi.li_GatePrime][mi.GPsp] += (mi.CAPcgsb)*mi.numberParallel;
        dQdx[mi.li_GatePrime][mi.GPbp] += (mi.CAPcgbb)*mi.numberParallel;
      } // WDLiu: CAPcrg already subtracted from all CAPcrgg below
      else if (mi.rgateMod == 2)
      {
        dQdx[mi.li_GatePrime][mi.GPgp] += (mi.CAPcggb)*mi.numberParallel;
        dQdx[mi.li_GatePrime][mi.GPdp] += (mi.CAPcgdb)*mi.numberParallel;
        dQdx[mi.li_GatePrime][mi.GPsp] += (mi.CAPcgsb)*mi.numberParallel;
        dQdx[mi.li_GatePrime][mi.GPbp] += (mi.CAPcgbb)*mi.numberParallel;
      }
      else if (mi.rgateMod == 3)
      {
        dQdx[mi.li_GateMid][mi.GMgm] += (+ mi.CAPcgmgmb)*mi.numberParallel;

        dQdx[mi.li_GateMid][mi.GMdp] += (mi.CAPcgmdb)*mi.numberParallel;
        dQdx[mi.li_GateMid][mi.GMsp] += (mi.CAPcgmsb)*mi.numberParallel;
        dQdx[mi.li_GateMid][mi.GMbp] += (mi.CAPcgmbb)*mi.numberParallel;

        dQdx[mi.li_DrainPrime][mi.DPgm] += (mi.CAPcdgmb)*mi.numberParallel;
        dQdx[mi.li_SourcePrime][mi.SPgm] += (mi.CAPcsgmb)*mi.numberParallel;
        dQdx[mi.li_BodyPrime][mi.BPgm] += (mi.CAPcbgmb)*mi.numberParallel;

        dQdx[mi.li_GatePrime][mi.GPgp] += (mi.CAPcggb)*mi.numberParallel;
        dQdx[mi.li_GatePrime][mi.GPdp] += (mi.CAPcgdb)*mi.numberParallel;
        dQdx[mi.li_GatePrime][mi.GPsp] += (mi.CAPcgsb)*mi.numberParallel;
        dQdx[mi.li_GatePrime][mi.GPbp] += (mi.CAPcgbb)*mi.numberParallel;
      }
      else
      {
        dQdx[mi.li_GatePrime][mi.GPgp] += (mi.CAPcggb)*mi.numberParallel;
        dQdx[mi.li_GatePrime][mi.GPdp] += (mi.CAPcgdb)*mi.numberParallel;
        dQdx[mi.li_GatePrime][mi.GPsp] += (mi.CAPcgsb)*mi.numberParallel;
        dQdx[mi.li_GatePrime][mi.GPbp] += (mi.CAPcgbb)*mi.numberParallel;
      }

      dQdx[mi.li_DrainPrime][mi.DPdp] += (mi.CAPcddb)*mi.numberParallel;
      dQdx[mi.li_DrainPrime][mi.DPgp] += (+ mi.CAPcdgb)*mi.numberParallel;
      dQdx[mi.li_DrainPrime][mi.DPsp] -= (- mi.CAPcdsb)*mi.numberParallel;
      dQdx[mi.li_DrainPrime][mi.DPbp] -= (- mi.CAPcdbb)*mi.numberParallel;

      dQdx[mi.li_SourcePrime][mi.SPdp] -= (- mi.CAPcsdb)*mi.numberParallel;
      dQdx[mi.li_SourcePrime][mi.SPgp] += (mi.CAPcsgb)*mi.numberParallel;
      dQdx[mi.li_SourcePrime][mi.SPsp] += (mi.CAPcssb)*mi.numberParallel;
      dQdx[mi.li_SourcePrime][mi.SPbp] -= (- mi.CAPcsbb)*mi.numberParallel;

      dQdx[mi.li_BodyPrime][mi.BPdp] += (mi.CAPcbdb)*mi.numberParallel;
      dQdx[mi.li_BodyPrime][mi.BPgp] += (mi.CAPcbgb)*mi.numberParallel;
      dQdx[mi.li_BodyPrime][mi.BPsp] += (mi.CAPcbsb)*mi.numberParallel;
      dQdx[mi.li_BodyPrime][mi.BPbp] += (mi.CAPcbbb)*mi.numberParallel;

      if (mi.rbodyMod)
      {
        dQdx[mi.li_DrainPrime][mi.DPdb] += (mi.CAPcdbdb)*mi.numberParallel;
        dQdx[mi.li_SourcePrime][mi.SPsb] -= (- mi.CAPcsbsb)*mi.numberParallel;

        dQdx[mi.li_DrainBody][mi.DBdp] += (mi.CAPcdbdb)*mi.numberParallel;
        dQdx[mi.li_DrainBody][mi.DBdb] += (- mi.CAPcdbdb)*mi.numberParallel;

        dQdx[mi.li_SourceBody][mi.SBsp] += (mi.CAPcsbsb)*mi.numberParallel;
        dQdx[mi.li_SourceBody][mi.SBsb] += (- mi.CAPcsbsb)*mi.numberParallel;
      }

      if (mi.trnqsMod)
      {
        dQdx[mi.li_Charge][mi.Qgp] += (- mi.CAPcqgb)*mi.numberParallel;
        dQdx[mi.li_Charge][mi.Qdp] += (- mi.CAPcqdb)*mi.numberParallel;
        dQdx[mi.li_Charge][mi.Qsp] += (- mi.CAPcqsb)*mi.numberParallel;
        dQdx[mi.li_Charge][mi.Qbp] += (- mi.CAPcqbb)*mi.numberParallel;
      }

#endif
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
  if (deviceMap.empty() ||
      ((deviceMap.find("M")!=deviceMap.end())
      && ((levelSet.find(14)!=levelSet.end()) || (levelSet.find(54)!=levelSet.end()))))
  {
    MOSFET1::registerDevice();

    Config<Traits>::addConfiguration()
      .registerDevice("m", 14)
      .registerDevice("m", 54)
      .registerModelType("pmos", 14)
      .registerModelType("nmos", 14)
      .registerModelType("pmos", 54)
      .registerModelType("nmos", 54);
  }
}

} // namespace MOSFET_B4
} // namespace Device
} // namespace Xyce

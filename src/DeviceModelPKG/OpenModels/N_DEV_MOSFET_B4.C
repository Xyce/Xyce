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


namespace Xyce {
namespace Device {

namespace MOSFET_B4 {

// Utility functions
namespace {

//-----------------------------------------------------------------------------
// Function      : convertVersToDouble
// Purpose       : converts a version string like "4.6.1" to a CMC-style
//                 double (4.61).  Having internal representation of version
//                 as a double eases work on support of multiple versions
// Special Notes :
// Scope         : (unnamed namespace under MOSFET_B4)
// Creator       : Tom Russo, SNL 1355
//-----------------------------------------------------------------------------
double convertVersToDouble(const std::string &versionString)
{
  double retval=0.0;

  std::string::size_type first_dot = versionString.find_first_of(".");
  std::string::size_type last_dot = versionString.find_last_of(".");

  if (first_dot == std::string::npos || first_dot == last_dot)
  {
    // then we are trivially convertible to double
    retval=std::stod(versionString);
  }
  else
  {
    retval=std::stod(versionString.substr(0,first_dot))
      + 0.1*convertVersToDouble(versionString.substr(first_dot+1));
  }

  return retval;
}
} // namespace (unnamed)

void Traits::loadInstanceParameters(ParametricData<MOSFET_B4::Instance> &p)
{
    p.addPar ("TEMP",0.0,&MOSFET_B4::Instance::temp)
     .setGivenMember(&MOSFET_B4::Instance::TEMPgiven)
     .setExpressionAccess(ParameterType::TIME_DEP)
     .setUnit(U_DEGC)
     .setCategory(CAT_NONE)
     .setDescription("Device temperature");

    p.addPar ("L",5.0e-6,&MOSFET_B4::Instance::l)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Length")
     .setLengthScaling(true);

    p.addPar ("W",5.0e-6,&MOSFET_B4::Instance::w)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Width")
     .setLengthScaling(true);

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
     .setDescription("Drain area")
     .setAreaScaling(true);

    p.addPar ("AS",0.0,&MOSFET_B4::Instance::sourceArea)
     .setGivenMember(&MOSFET_B4::Instance::sourceAreaGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source area")
     .setAreaScaling(true);

    p.addPar ("PD",0.0,&MOSFET_B4::Instance::drainPerimeter)
     .setGivenMember(&MOSFET_B4::Instance::drainPerimeterGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Drain perimeter")
     .setLengthScaling(true);

    p.addPar ("PS",0.0,&MOSFET_B4::Instance::sourcePerimeter)
     .setGivenMember(&MOSFET_B4::Instance::sourcePerimeterGiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Source perimeter")
     .setLengthScaling(true);

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

    p.addPar ("DELVT0",0.0,&MOSFET_B4::Instance::delvto)
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

    p.addPar ("DTEMP",0.0,&MOSFET_B4::Instance::dtemp)
     .setGivenMember(&MOSFET_B4::Instance::dtempGiven)
     .setUnit(U_DEGC)
     .setCategory(CAT_NONE)
     .setDescription("Device delta temperature");

    // This tells the parser that IC1,IC2,and IC3 are to be input as a vector of "IC"
    p.makeVector ("IC",3);
}

void Traits::loadModelParameters(ParametricData<MOSFET_B4::Model> &p)
{
    p.addPar ("EOT",15.0e-10,&MOSFET_B4::Model::eot)
     .setGivenMember(&MOSFET_B4::Model::eotGiven)
     .setUnit(U_METER)
     .setCategory(CAT_PROCESS)
     .setDescription("Equivalent gate oxide thickness in meters");

    p.addPar ("VDDEOT",1.5,&MOSFET_B4::Model::vddeot)
     .setUnit(U_VOLT)
     .setCategory(CAT_BASIC)
     .setDescription("Voltage for extraction of equivalent gate oxide thickness");

    p.addPar ("TEMPEOT",300.15,&MOSFET_B4::Model::tempeot)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setMinimumVersion(4.70)
     .setDescription("Temperature for extraction of EOT");

    p.addPar ("LEFFEOT",1e-6,&MOSFET_B4::Model::leffeot)
     .setUnit(U_METER)
     .setCategory(CAT_BASIC)
     .setMinimumVersion(4.70)
     .setDescription("Effective length for extraction of EOT");

    p.addPar ("WEFFEOT",10e-6,&MOSFET_B4::Model::weffeot)
     .setUnit(U_METER)
     .setCategory(CAT_BASIC)
     .setMinimumVersion(4.70)
     .setDescription("Effective width for extraction of EOT");

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
     .setGivenMember(&MOSFET_B4::Model::toxmGiven)
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

    p.addPar ("DVTP2",0.0,&MOSFET_B4::Model::dvtp2)
     .setUnit(U_VMX)
     .setCategory(CAT_BASIC)
     .setMinimumVersion(4.70)
     .setDescription("3rd parameter for Vth shift due to pocket");

    p.addPar ("DVTP3",0.0,&MOSFET_B4::Model::dvtp3)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setMinimumVersion(4.70)
     .setDescription("4th parameter for Vth shift due to pocket");

    p.addPar ("DVTP4",0.0,&MOSFET_B4::Model::dvtp4)
     .setUnit(U_VOLTM1)
     .setCategory(CAT_BASIC)
     .setMinimumVersion(4.70)
     .setDescription("5th parameter for Vth shift due to pocket");

    p.addPar ("DVTP5",0.0,&MOSFET_B4::Model::dvtp4)
     .setUnit(U_VOLT)
     .setCategory(CAT_BASIC)
     .setMinimumVersion(4.70)
     .setDescription("6th parameter for Vth shift due to pocket");

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

    p.addPar ("UCS",1.67,&MOSFET_B4::Model::ucs)
     .setUnit(U_NONE)
     .setCategory(CAT_BASIC)
     .setMinimumVersion(4.70)
     .setDescription("Colombic scattering exponent");

    p.addPar ("UTE",-1.5,&MOSFET_B4::Model::ute)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Temperature coefficient of mobility");

    p.addPar ("UCSTE",-4.775e-3,&MOSFET_B4::Model::ucste)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Temperature coefficient of colombic mobility");

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

    p.addPar ("TNFACTOR",0.0,&MOSFET_B4::Model::tnfactor)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Temperature parameter for nfactor");

    p.addPar ("TETA0",0.0,&MOSFET_B4::Model::teta0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Temperature parameter for eta0");

    p.addPar ("TVOFFCV",0.0,&MOSFET_B4::Model::tvoffcv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Temperature parameter for tvoffcv");

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
     .setGivenMember(&MOSFET_B4::Model::cfGiven)
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

    p.addPar ("RGIDL",1.0,&MOSFET_B4::Model::rgidl)
     .setUnit(U_NONE)
     .setCategory(CAT_GDLEAKAGE)
     .setMinimumVersion(4.70)
     .setDescription("GIDL vg parameter");

    p.addPar ("KGIDL",0.0,&MOSFET_B4::Model::kgidl)
     .setUnit(U_VOLT)
     .setCategory(CAT_GDLEAKAGE)
     .setMinimumVersion(4.70)
     .setDescription("GIDL vb parameter");

    p.addPar ("FGIDL",0.0,&MOSFET_B4::Model::fgidl)
     .setUnit(U_VOLT)
     .setCategory(CAT_GDLEAKAGE)
     .setMinimumVersion(4.70)
     .setDescription("GIDL vb parameter");

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

    p.addPar ("RGISL",1.0,&MOSFET_B4::Model::rgisl)
     .setUnit(U_NONE)
     .setCategory(CAT_GDLEAKAGE)
     .setMinimumVersion(4.70)
     .setDescription("Parameter for GISL gate bias dependence");

    p.addPar ("KGISL",0.0,&MOSFET_B4::Model::kgisl)
     .setUnit(U_VOLT)
     .setCategory(CAT_GDLEAKAGE)
     .setMinimumVersion(4.70)
     .setDescription("Parameter for GISL body bias dependence");

    p.addPar ("FGISL",0.0,&MOSFET_B4::Model::fgisl)
     .setUnit(U_VOLT)
     .setCategory(CAT_GDLEAKAGE)
     .setMinimumVersion(4.70)
     .setDescription("Parameter for GISL body bias dependence");

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

    p.addPar ("JTWEFF",0.0,&MOSFET_B4::Model::jtweff)
     .setUnit(U_METER)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("TAT current width dependence");

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

    p.addPar ("LDVTP2",0.0,&MOSFET_B4::Model::ldvtp2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Length dependence of dvtp2");

    p.addPar ("LDVTP3",0.0,&MOSFET_B4::Model::ldvtp3)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Length dependence of dvtp3");

    p.addPar ("LDVTP4",0.0,&MOSFET_B4::Model::ldvtp4)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Length dependence of dvtp4");

    p.addPar ("LDVTP5",0.0,&MOSFET_B4::Model::ldvtp5)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Length dependence of dvtp5");

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

    p.addPar ("LUCSTE",0.0,&MOSFET_B4::Model::lucste)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Length dependence of ucste");

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

    p.addPar ("LRGIDL",0.0,&MOSFET_B4::Model::lrgidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Length dependence of rgidl");

    p.addPar ("LKGIDL",0.0,&MOSFET_B4::Model::lkgidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Length dependence of kgidl");

    p.addPar ("LFGIDL",0.0,&MOSFET_B4::Model::lfgidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Length dependence of fgidl");

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

    p.addPar ("LRGISL",0.0,&MOSFET_B4::Model::lrgisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Length dependence of rgisl");

    p.addPar ("LKGISL",0.0,&MOSFET_B4::Model::lkgisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Length dependence of kgisl");

    p.addPar ("LFGISL",0.0,&MOSFET_B4::Model::lfgisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Length dependence of fgisl");

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

    p.addPar ("LUCS",0.0,&MOSFET_B4::Model::lucs)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription(" Length dependence of ucs");

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

    p.addPar ("LTNFACTOR",0.0,&MOSFET_B4::Model::ltnfactor)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Length dependence of tnfactor");

    p.addPar ("LTETA0",0.0,&MOSFET_B4::Model::lteta0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Length dependence of teta0");

    p.addPar ("LTVOFFCV",0.0,&MOSFET_B4::Model::ltvoffcv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Length dependence of tvoffcv");

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

    p.addPar ("WDVTP2",0.0,&MOSFET_B4::Model::wdvtp2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Width dependence of dvtp2");

    p.addPar ("WDVTP3",0.0,&MOSFET_B4::Model::wdvtp3)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Width dependence of dvtp3");

    p.addPar ("WDVTP4",0.0,&MOSFET_B4::Model::wdvtp4)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Width dependence of dvtp4");

    p.addPar ("WDVTP5",0.0,&MOSFET_B4::Model::wdvtp5)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Width dependence of dvtp5");

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

    p.addPar ("WUCSTE",0.0,&MOSFET_B4::Model::wucste)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Width dependence of ucste");

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

    p.addPar ("WRGIDL",0.0,&MOSFET_B4::Model::wrgidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Width dependence of rgidl");

    p.addPar ("WKGIDL",0.0,&MOSFET_B4::Model::wkgidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Width dependence of kgidl");

    p.addPar ("WFGIDL",0.0,&MOSFET_B4::Model::wfgidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Width dependence of fgidl");

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

    p.addPar ("WRGISL",0.0,&MOSFET_B4::Model::wrgisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Width dependence of rgisl");

    p.addPar ("WKGISL",0.0,&MOSFET_B4::Model::wkgisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Width dependence of kgisl");

    p.addPar ("WFGISL",0.0,&MOSFET_B4::Model::wfgisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Width dependence of fgisl");

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

    p.addPar ("WUCS",0.0,&MOSFET_B4::Model::wucs)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Width dependence of ucs");

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

    p.addPar ("WTNFACTOR",0.0,&MOSFET_B4::Model::wtnfactor)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Width dependence of tnfactor");

    p.addPar ("WTETA0",0.0,&MOSFET_B4::Model::wteta0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Width dependence of teta0");

    p.addPar ("WTVOFFCV",0.0,&MOSFET_B4::Model::wtvoffcv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Width dependence of tvoffcv");

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

    p.addPar ("PDVTP2",0.0,&MOSFET_B4::Model::pdvtp2)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Cross-term dependence of dvtp2");

    p.addPar ("PDVTP3",0.0,&MOSFET_B4::Model::pdvtp3)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Cross-term dependence of dvtp3");

    p.addPar ("PDVTP4",0.0,&MOSFET_B4::Model::pdvtp4)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Cross-term dependence of dvtp4");

    p.addPar ("PDVTP5",0.0,&MOSFET_B4::Model::pdvtp5)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Cross-term dependence of dvtp5");

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

    p.addPar ("PUCSTE",0.0,&MOSFET_B4::Model::pucste)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Cross-term dependence of ucste");

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

    p.addPar ("PRGIDL",0.0,&MOSFET_B4::Model::prgidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Cross-term dependence of rgidl");

    p.addPar ("PKGIDL",0.0,&MOSFET_B4::Model::pkgidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Cross-term dependence of kgidl");

    p.addPar ("PFGIDL",0.0,&MOSFET_B4::Model::pfgidl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Cross-term dependence of fgidl");

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

    p.addPar ("PRGISL",0.0,&MOSFET_B4::Model::prgisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Cross-term dependence of rgisl");

    p.addPar ("PKGISL",0.0,&MOSFET_B4::Model::pkgisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Cross-term dependence of kgisl");

    p.addPar ("PFGISL",0.0,&MOSFET_B4::Model::pfgisl)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Cross-term dependence of fgisl");

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

    p.addPar ("PUCS",0.0,&MOSFET_B4::Model::pucs)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Cross-term dependence of ucs");

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

    p.addPar ("PTNFACTOR",0.0,&MOSFET_B4::Model::ptnfactor)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Cross-term dependence of tnfactor");

    p.addPar ("PTETA0",0.0,&MOSFET_B4::Model::pteta0)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Cross-term dependence of teta0");

    p.addPar ("PTVOFFCV",0.0,&MOSFET_B4::Model::ptvoffcv)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Cross-term dependence of tvoffcv");

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

    p.addPar ("TNOIC",0.0,&MOSFET_B4::Model::tnoic)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Thermal noise parameter");

    p.addPar ("RNOIA",0.577,&MOSFET_B4::Model::rnoia)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Thermal noise coefficient");

    p.addPar ("RNOIB",0.5164,&MOSFET_B4::Model::rnoib)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Thermal noise coefficient");

    p.addPar ("RNOIC",0.395,&MOSFET_B4::Model::rnoic)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setMinimumVersion(4.70)
     .setDescription("Thermal noise coefficient");

    p.addPar ("GIDLCLAMP",-1.0e-5,&MOSFET_B4::Model::gidlclamp)
      .setUnit(U_NONE)
      .setCategory(CAT_NONE)
      .setMinimumVersion(4.82)
      .setDescription("gidl clamp value");

    // NOTE:  This is NOT a typo.  The BSIM4.8.2 model defines a model
    // parameter "IDOVVDS" but associates it with an internal variable
    // idovvdsc.
    p.addPar ("IDOVVDS",1e-9,&MOSFET_B4::Model::idovvdsc)
      .setUnit(U_NONE)
      .setCategory(CAT_NONE)
      .setMinimumVersion(4.82)
      .setDescription("noise clamping limit parameter");

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

    p.addPar ("RGEOMOD",0,&MOSFET_B4::Model::rgeoMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("S/D resistance and contact model selector");

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

    p.addPar ("MTRLCOMPATMOD",0,&MOSFET_B4::Model::mtrlCompatMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setMinimumVersion(4.70)
     .setDescription("New material Mod backward compatibility selector");

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

    p.addPar ("GIDLMOD",0,&MOSFET_B4::Model::gidlMod)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setMinimumVersion(4.70)
     .setDescription("parameter for GIDL selector");

    p.addPar ("PARAMCHK",1,&MOSFET_B4::Model::paramChk)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Model parameter checking selector");

    p.addPar ("BINUNIT",1,&MOSFET_B4::Model::binUnit)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Bin  unit  selector");

    p.addPar ("VERSION",std::string("4.8.2"),&MOSFET_B4::Model::version)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("parameter for model version");
}

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
    if (sourceAreaGiven) { sourceArea *= getDeviceOptions().lengthScale * getDeviceOptions().lengthScale ; }
    if (drainAreaGiven) { drainArea *= getDeviceOptions().lengthScale * getDeviceOptions().lengthScale ; }
    if (drainPerimeterGiven) { drainPerimeter *= getDeviceOptions().lengthScale; }
    if (sourcePerimeterGiven) { sourcePerimeter *= getDeviceOptions().lengthScale; }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       : wrapper for version-specific processParams calls
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, 1355
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  return (this->*processParamsPtr_)();
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
// Function      : setupVersionPointers_
// Purpose       : Set up the various member-function pointers needed
//                 to make us perform version-specific operations
// Special Notes :
// Scope         : private
// Creator       : Tom Russo, SNL, 1355
// Creation Date : 24 Aug 2022
//-----------------------------------------------------------------------------
void Instance::setupVersionPointers_()
{
  // set up member function pointers based on version:
  // nonsensicial if structure is a placeholder for when we actually
  // have more than one version.
  if (model_.versionDouble == 4.61)
  {
    processParamsPtr_ = &Instance::processParams4p61_;
    updateTemperaturePtr_ = &Instance::updateTemperature4p61_;
    updateIntermediateVarsPtr_ = &Instance::updateIntermediateVars4p61_;
    setupNoiseSourcesPtr_ = &Instance::setupNoiseSources4p61_;
    getNoiseSourcesPtr_ = &Instance::getNoiseSources4p61_;
    RdsEndIsoPtr_ = &Instance::RdsEndIso4p61_;
  }
  else   if (model_.versionDouble == 4.70)
  {
    processParamsPtr_ = &Instance::processParams4p70_;
    updateTemperaturePtr_ = &Instance::updateTemperature4p70_;
    updateIntermediateVarsPtr_ = &Instance::updateIntermediateVars4p70_;
    setupNoiseSourcesPtr_ = &Instance::setupNoiseSources4p70_;
    getNoiseSourcesPtr_ = &Instance::getNoiseSources4p70_;
    RdsEndIsoPtr_ = &Instance::RdsEndIso4p70_;
  }
  else if (model_.versionDouble >= 4.80)
  {
    processParamsPtr_ = &Instance::processParams4p82_;
    updateTemperaturePtr_ = &Instance::updateTemperature4p82_;
    updateIntermediateVarsPtr_ = &Instance::updateIntermediateVars4p82_;
    setupNoiseSourcesPtr_ = &Instance::setupNoiseSources4p82_;
    getNoiseSourcesPtr_ = &Instance::getNoiseSources4p82_;
    RdsEndIsoPtr_ = &Instance::RdsEndIso4p82_;
  }
  else
  {
    processParamsPtr_ = &Instance::processParams4p82_;
    updateTemperaturePtr_ = &Instance::updateTemperature4p82_;
    updateIntermediateVarsPtr_ = &Instance::updateIntermediateVars4p82_;
    setupNoiseSourcesPtr_ = &Instance::setupNoiseSources4p82_;
    getNoiseSourcesPtr_ = &Instance::getNoiseSources4p82_;
    RdsEndIsoPtr_ = &Instance::RdsEndIso4p82_;
  }
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
    processParamsPtr_(&Instance::processParams4p70_),
    updateTemperaturePtr_(&Instance::updateTemperature4p70_),
    updateIntermediateVarsPtr_(&Instance::updateIntermediateVars4p70_),
    setupNoiseSourcesPtr_(&Instance::setupNoiseSources4p70_),
    getNoiseSourcesPtr_(&Instance::getNoiseSources4p70_),
    RdsEndIsoPtr_(&Instance::RdsEndIso4p70_),
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
    toxp(0.0),
    coxp(0.0),
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
    noiGd0(0.0),
    Coxeff(0.0),
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
    dtemp(0.0),
    dtempGiven(false),
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

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  checkParamVersions(model_.versionDouble);
  setupVersionPointers_();

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // if options scale has been set in the netlist, apply it.
  applyScale ();

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
// Special Notes : Wrapper for calls through a member function pointer
// Creator       : Tom Russo, SNL, 1355
// Creation Date : 24 Aug 2022
//-----------------------------------------------------------------------------
bool Instance::updateTemperature (const double & temp_tmp)
{
  return (this->*updateTemperaturePtr_)(temp_tmp);
}


//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : Wrapper for real version-specific updateIntermediateVars call
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, 1355
// Creation Date : 24 Aug 2022
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  return (this->*updateIntermediateVarsPtr_)();
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
// Function      : setupNoiseSources
// Purpose       : Version-generic setup of noise sources
// Special Notes :  wrapper for version-specific code
// Scope         : public
// Creator       : Tom Russo, SNL, 1355
// Creation Date : 24 Aug 2022
//-----------------------------------------------------------------------------
void Instance::setupNoiseSources (Xyce::Analysis::NoiseData & noiseData)
{
  (this->*setupNoiseSourcesPtr_)(noiseData);
}

//-----------------------------------------------------------------------------
// Function      : getNoiseSources
// Purpose       : Version-generic getNoiseSources
// Special Notes :
// Scope         : private
// Creator       : Tom Russo, SNL, 1355
// Creation Date : 24 Aug 2022
//-----------------------------------------------------------------------------
void Instance::getNoiseSources (Xyce::Analysis::NoiseData & noiseData)
{
  (this->*getNoiseSourcesPtr_)(noiseData);
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
  return (this->*RdsEndIsoPtr_)(Weffcj, Rsh, DMCG, DMCI, DMDG,
                                nuEnd, rgeo, Type, Rend);
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
// Function      : checkAndFixVersion_
// Purpose       : check the version string given on the model line (if any)
//                 set the double precision version, and reset as needed to
//                 match a version that actually exists
// Special Notes :
// Scope         : private
// Creator       : Tom Russo, SNL 1355
//-----------------------------------------------------------------------------
void Model::checkAndFixVersion_()
{
  if (given("version"))
    versionDouble = convertVersToDouble(version);
  else
    versionDouble=4.82;     // not strictly necessary, because we set this
                            // in the constructor initializers already, but
                            // let's play it safe

  if (versionDouble < 4.61)
  {
    UserWarning(*this) << "Model card specifies BSIM4 version " << version
         << " which is older than the oldest version supported in Xyce (4.6.1). "
                       << " Using oldest version available."
                       << std::endl;
    versionDouble=4.61;
  }
  else if (versionDouble < 4.70)
  {
    if (versionDouble != 4.61)
    {
      UserWarning(*this) << "Model card specifies BSIM4 version " << version
         << " not supported by Xyce. "
         << " Using version 4.6.1, the supported version prior to the requested version "
                       << std::endl;
    }
    versionDouble=4.61;
  }
  else if (versionDouble < 4.80)
  {
    if (versionDouble != 4.70)
    {
      UserWarning(*this) << "Model card specifies BSIM4 version " << version
                         << " not supported by Xyce. "
                         << " Using 4.7.0 instead."
                         << std::endl;
    }
    versionDouble=4.70;
  }
  else if (versionDouble >= 4.80)
  {
    if (versionDouble < 4.82)
    {
      UserWarning(*this) << "Model card specifies BSIM4 version " << version
                         << " not supported by Xyce. "
                         << " Using 4.8.2 instead."
                         << std::endl;
    }
    if (versionDouble > 4.82)
    {
      UserWarning(*this) << "Model card specifies BSIM4 version " << version
                         << " which is newer than the latest version supported in Xyce (4.8.2)"
                         << " Using 4.8.2 instead."
                         << std::endl;
    }
    // DON'T reset the version value, though, because 4.8.2 has some
    // conditionals that make it work differently if version <= 4.80
  }
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes : Wrapper for the "guts" that are hidden in a function pointer
// Scope         : public
// Creator       : Tom Russo, SNL, 01355
// Creation Date : 24 Aug 2022
//-----------------------------------------------------------------------------
bool Model::processParams()
{
  return (this->*processParamsPtr_)();
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
// Function      : setupVersionPointers_
// Purpose       : Set up the various member-function pointers needed
//                 to make us perform version-specific operations
// Special Notes :
// Scope         : private
// Creator       : Tom Russo, SNL, 1355
// Creation Date : 24 Aug 2022
//-----------------------------------------------------------------------------
void Model::setupVersionPointers_()
{
  // set up member function pointers based on version:
  // nonsensicial if structure is a placeholder for when we actually
  // have more than one version.
  if (versionDouble == 4.61)
  {
    processParamsPtr_ = &Model::processParams4p61_;
  }
  else if (versionDouble == 4.70)
  {
    processParamsPtr_ = &Model::processParams4p70_;
  }
  else if (versionDouble >= 4.80)
  {
    processParamsPtr_ = &Model::processParams4p82_;
  }
  else
  {
    processParamsPtr_ = &Model::processParams4p82_;
  }
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
    processParamsPtr_(&Model::processParams4p70_),
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
    rgeoMod(0),
    mtrlMod(0),
    mtrlCompatMod(0),
    igcMod(0),
    igbMod(0),
    tempMod(0),
    gidlMod(0),
    binUnit(0),
    paramChk(0),
    version("4.8.2"),
    versionDouble(4.82),
    eot(0.0),
    vddeot(0.0),
    tempeot(0.0),
    leffeot(0.0),
    weffeot(0.0),
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
    dvtp2(0.0),
    dvtp3(0.0),
    dvtp4(0.0),
    dvtp5(0.0),
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
    ucs(0.0),
    ute(0.0),
    ucste(0.0),
    voff(0.0),
    tvoff(0.0),
    tnfactor(0.0),
    teta0(0.0),
    tvoffcv(0.0),
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
    rgidl(0.0),
    kgidl(0.0),
    fgidl(0.0),
    egidl(0.0),
    agisl(0.0),
    bgisl(0.0),
    cgisl(0.0),
    rgisl(0.0),
    kgisl(0.0),
    fgisl(0.0),
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
    jtweff(0.0),
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
    tnoic(0.0),
    rnoia(0.0),
    rnoib(0.0),
    rnoic(0.0),
    gidlclamp(-1.0e-5),
    idovvdsc(1e-9),
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
    ldvtp2(0.0),
    ldvtp3(0.0),
    ldvtp4(0.0),
    ldvtp5(0.0),
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
    lucste(0.0),
    lvoff(0.0),
    ltvoff(0.0),
    ltnfactor(0.0),
    lteta0(0.0),
    ltvoffcv(0.0),
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
    lrgidl(0.0),
    lkgidl(0.0),
    lfgidl(0.0),
    legidl(0.0),
    lagisl(0.0),
    lbgisl(0.0),
    lcgisl(0.0),
    lrgisl(0.0),
    lkgisl(0.0),
    lfgisl(0.0),
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
    wdvtp2(0.0),
    wdvtp3(0.0),
    wdvtp4(0.0),
    wdvtp5(0.0),
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
    wucs(0.0),
    wute(0.0),
    wucste(0.0),
    wvoff(0.0),
    wtvoff(0.0),
    wtnfactor(0.0),
    wteta0(0.0),
    wtvoffcv(0.0),
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
    wrgidl(0.0),
    wkgidl(0.0),
    wfgidl(0.0),
    wegidl(0.0),
    wagisl(0.0),
    wbgisl(0.0),
    wcgisl(0.0),
    wrgisl(0.0),
    wkgisl(0.0),
    wfgisl(0.0),
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
    pdvtp2(0.0),
    pdvtp3(0.0),
    pdvtp4(0.0),
    pdvtp5(0.0),
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
    pucs(0.0),
    pute(0.0),
    pucste(0.0),
    pvoff(0.0),
    ptvoff(0.0),
    ptnfactor(0.0),
    pteta0(0.0),
    ptvoffcv(0.0),
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
    prgidl(0.0),
    pkgidl(0.0),
    pfgidl(0.0),
    pegidl(0.0),
    pagisl(0.0),
    pbgisl(0.0),
    pcgisl(0.0),
    prgisl(0.0),
    pkgisl(0.0),
    pfgisl(0.0),
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
    eotGiven(false),
    toxeGiven(false),
    toxpGiven(false),
    toxmGiven(false),
    dtoxGiven(false),
    cgdoGiven(false),
    dlcGiven(false),
    cgsoGiven(false),
    cgboGiven(false),
    cfGiven(false),
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

  checkAndFixVersion_();
  checkParamVersions(versionDouble);
  setupVersionPointers_();

  // Set any non-constant parameter defaults:
  if (!given("TNOM"))
    tnom = getDeviceOptions().tnom;

  // Calculate any parameters specified as expressions:

  updateDependentParameters();

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

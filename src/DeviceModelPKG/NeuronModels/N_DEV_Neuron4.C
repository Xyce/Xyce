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
// Creator        : Richard Schiek, Electrical and Microsytem Modeling
//
// Creation Date  : 06/10/09
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------


// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Neuron4.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_DEV_Neuron.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Math.h>

namespace Xyce {
namespace Device {


namespace Neuron4 {


void Traits::loadInstanceParameters(ParametricData<Neuron4::Instance> &p)
{
  p.addPar ("R",0.0,&Neuron4::Instance::rInt)
   .setGivenMember(&Neuron4::Instance::rIntGiven)
   .setUnit(U_OHMMM1)
   .setCategory(CAT_NONE)
   .setDescription("Intracellular resistivity");

  p.addPar ("A",0.0,&Neuron4::Instance::radius)
   .setGivenMember(&Neuron4::Instance::radiusGiven)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Segment radius");

  p.addPar ("L",0.0,&Neuron4::Instance::length)
   .setGivenMember(&Neuron4::Instance::lengthGiven)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Cable length");

  p.addPar ("RPS",1.0e-6,&Neuron4::Instance::rIntPrevious)
   .setGivenMember(&Neuron4::Instance::rIntPreviousGiven)
   .setUnit(U_OHMMM1)
   .setCategory(CAT_NONE)
   .setDescription("Previous segment,intracellular resistivity");

  p.addPar ("APS",0.0,&Neuron4::Instance::radiusPrevious)
   .setGivenMember(&Neuron4::Instance::radiusPreviousGiven)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Previous segment,segment radius");

  p.addPar ("LPS",0.0,&Neuron4::Instance::lengthPrevious)
   .setGivenMember(&Neuron4::Instance::lengthPreviousGiven)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Previous segment length");

  p.addPar ("RNS",1.0e-6,&Neuron4::Instance::rIntNext)
   .setGivenMember(&Neuron4::Instance::rIntNextGiven)
   .setUnit(U_OHMMM1)
   .setCategory(CAT_NONE)
   .setDescription("Next segment,intracellular resistivity");

  p.addPar ("ANS",0.0,&Neuron4::Instance::radiusNext)
   .setGivenMember(&Neuron4::Instance::radiusNextGiven)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Next segment,segment radius");

  p.addPar ("LNS",0.0,&Neuron4::Instance::lengthNext)
   .setGivenMember(&Neuron4::Instance::lengthNextGiven)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Next segment length");

  p.addPar("N",0,&Neuron4::Instance::nSeg)
   .setGivenMember(&Neuron4::Instance::nSegGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Number of segments");
}

void Traits::loadModelParameters(ParametricData<Neuron4::Model> &p)
{
  p.addPar ("CMEM",0.0,&Neuron4::Model::cMem)
   .setGivenMember(&Neuron4::Model::cMemGiven)
   .setUnit(U_FARAD)
   .setCategory(CAT_NONE)
   .setDescription("Membrane capacitance");

  p.addPar ("GMEM",0.0,&Neuron4::Model::gMem)
   .setGivenMember(&Neuron4::Model::gMemGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("Membrane conductance");

  p.addPar ("VREST",0.0,&Neuron4::Model::vRest)
   .setGivenMember(&Neuron4::Model::vRestGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Resting potential");

  p.addPar ("EK",0.0,&Neuron4::Model::eK)
   .setGivenMember(&Neuron4::Model::eKGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Potassium resting potential");

  p.addPar ("GK",0.0,&Neuron4::Model::gK)
   .setGivenMember(&Neuron4::Model::gKGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("Potassium base conductance");

  p.addPar ("ENA",0.0,&Neuron4::Model::eNa)
   .setGivenMember(&Neuron4::Model::eNaGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Sodium resting potential");

  p.addPar ("GNA",0.0,&Neuron4::Model::gNa)
   .setGivenMember(&Neuron4::Model::gNaGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("Sodium base conductance");

  p.addPar ("EA",0.0,&Neuron4::Model::eA)
   .setGivenMember(&Neuron4::Model::eAGiven)
   .setUnit(U_CM)
   .setCategory(CAT_NONE)
   .setDescription("a-current rest potential");

  p.addPar ("GA",0.0,&Neuron4::Model::gA)
   .setGivenMember(&Neuron4::Model::gAGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("a-current base conductance");

  p.addPar ("ECA",0.0,&Neuron4::Model::eCa)
   .setGivenMember(&Neuron4::Model::eCaGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Calcium rest potential");

  p.addPar ("GCA",0.0,&Neuron4::Model::gCa)
   .setGivenMember(&Neuron4::Model::gCaGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("Calcium base conductance");

  p.addPar ("EKCA",0.0,&Neuron4::Model::eKCa)
   .setGivenMember(&Neuron4::Model::eKCaGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("Potassium-calcium rest potential");

  p.addPar ("GKCA",0.0,&Neuron4::Model::gKCa)
   .setGivenMember(&Neuron4::Model::gKCaGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("Potassium-calcium base conductance");

  p.addPar ("CAINIT",0.0,&Neuron4::Model::CaInit)
   .setGivenMember(&Neuron4::Model::CaInitGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("initial intra-cellular calcium concentration");

  p.addPar ("CAGAMMA",0.0,&Neuron4::Model::CaGamma)
   .setGivenMember(&Neuron4::Model::CaGammaGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("calcium current to concentration multiplier");

  p.addPar ("CATAU",0.0,&Neuron4::Model::CaTau)
   .setGivenMember(&Neuron4::Model::CaTauGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("calcium removal time constant");

  p.addPar ("R",0.0,&Neuron4::Model::rInt)
   .setGivenMember(&Neuron4::Model::rIntGiven)
   .setUnit(U_OHMMM1)
   .setCategory(CAT_NONE)
   .setDescription("Intracellular resistivity");

  p.addPar ("A",0.0,&Neuron4::Model::radius)
   .setGivenMember(&Neuron4::Model::radiusGiven)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Segment radius");

  p.addPar ("L",0.0,&Neuron4::Model::length)
   .setGivenMember(&Neuron4::Model::lengthGiven)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Cable length");

  p.addPar ("RPS",1.0e-6,&Neuron4::Model::rIntPrevious)
   .setGivenMember(&Neuron4::Model::rIntPreviousGiven)
   .setUnit(U_OHMMM1)
   .setCategory(CAT_NONE)
   .setDescription("Previous segment,intracellular resistivity");

  p.addPar ("APS",0.0,&Neuron4::Model::radiusPrevious)
   .setGivenMember(&Neuron4::Model::radiusPreviousGiven)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Previous segment,segment radius");

  p.addPar ("LPS",0.0,&Neuron4::Model::lengthPrevious)
   .setGivenMember(&Neuron4::Model::lengthPreviousGiven)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Previous segment length");

  p.addPar ("RNS",1.0e-6,&Neuron4::Model::rIntNext)
   .setGivenMember(&Neuron4::Model::rIntNextGiven)
   .setUnit(U_OHMMM1)
   .setCategory(CAT_NONE)
   .setDescription("Next segment,intracellular resistivity");

  p.addPar ("ANS",0.0,&Neuron4::Model::radiusNext)
   .setGivenMember(&Neuron4::Model::radiusNextGiven)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Next segment,segment radius");

  p.addPar ("LNS",0.0,&Neuron4::Model::lengthNext)
   .setGivenMember(&Neuron4::Model::lengthNextGiven)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Next segment length");

  p.addPar("N",0,&Neuron4::Model::nSeg)
   .setGivenMember(&Neuron4::Model::nSegGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Number of segments");
}


// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::processParams()
{
  // If there are any time dependent parameters, set their values at for
  // the current time.

  // pull unspecified params out of the model if they weren't specified here
  if( !nSegGiven )
  {
    if ( model_.nSegGiven )
    {
      nSeg = model_.nSeg;
    }
    else
    {
      nSeg = 10;
    }
  }

  if( !rIntGiven && model_.rIntGiven )
  {
    rInt = model_.rInt;
  }
  if( !radiusGiven && model_.radiusGiven )
  {
    radius = model_.radius;
  }
  if( !lengthGiven && model_.lengthGiven )
  {
    length = model_.length;
  }

  // if nSeg is still unknown then estimate it via lambda-d rule (TO DO)

  // now that nSeg, length and radius are set calculate segment area
  segArea = 2.0 * M_PI * radius * length / nSeg;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp)
{
  bool bsuccess = true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       Miter,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    rInt(0.0),
    radius(0.0),
    length(0.0),
    segArea(0.0),
    nSeg(0),
    rIntGiven(false),
    radiusGiven(false),
    lengthGiven(false),
    nSegGiven(false),
    rIntPrevious(0.0),
    radiusPrevious(0.0),
    lengthPrevious(0.0),
    rIntNext(0.0),
    radiusNext(0.0),
    lengthNext(0.0),
    rIntPreviousGiven(false),
    radiusPreviousGiven(false),
    lengthPreviousGiven(false),
    rIntNextGiven(false),
    radiusNextGiven(false),
    lengthNextGiven(false),
    kcl1Fvalue(0.0),
    kcl2Fvalue(0.0),
    dkcl1F_dVin(0.0),
    dkcl1F_dVs0(0.0),
    dkcl2F_dVout(0.0),
    dkcl2F_dVsn(0.0),
    li_Pos(0),
    li_Neg(0),
    APosEquPosNodeOffset(0),
    APosEquNextNodeOffset(0),
    ANegEquNegNodeOffset(0),
    ANegEquLastNodeOffset(0)
{
  // we have pased the model and instance parameters so now we can calculate
  // the number of internal vars
  numExtVars   = 2;  // input and output voltage

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line
  setParams (IB.params);

  // Set any non-constant parameter defaults:

  //if (!given("TEMP"))
  //  temp = getDeviceOptions().temp.dVal();

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

  // now we can calculate the number of internal vars
  numIntVars   = nSeg * 10;  // ion channel vars + one internal voltage node var for each segment
  numStateVars = nSeg * 2;   // two currents per segment
  int numVars = numExtVars + numIntVars;

  //
  // i(n) - I(n)/A - g(n,n+1) * (V(n+1) - V(n)) - g(n,n-1) * (V(n-1) - V(n)) + Cm dV(n)/d(t) = 0
  //
  // Vin   : g(0,1) * (V(1) - Vin ) = 0
  // Vout  : g(n,n-1) * (V(n-1) - Vout) = 0
  // Vnode : i(n) - I(n)/A - g(n,n+1) * (V(n+1) - V(n)) - g(n,n-1) * (V(n-1) - V(n)) + Cm dV(n)/d(t) = 0
  //        plus node supporting equations (a, b, m)
  //
  // jacobian format
  //             Vin      Vout      V1   n   m   h  a  b  M  H  c  Ca  V2   n   m   h    V(nSeg)   n   m   h
  // kcl Vin    -g(0,1)            g(0,1)
  // kcl Vout           -g(n,n-1)                                      g(n,n-1)
  // kcl V1     yes                yes yes yes yes  y  y  y  y  y     yes
  // n                             yes yes
  // m                             yes     yes
  // h                             yes         yes
  // a                             yes             yes
  // b                             yes                yes
  // M                             yes                   yes
  // H                             yes                     yes
  // c                             yes                        yes yes
  // Ca                            yes                   yes yes  yes
  //
  //
  // jacobian element count by row:
  // vin:  2
  // vout: 2
  // v1:   11
  // n1:   2
  // m1:   2
  // h1:   2
  // a1:   2
  // b1:   2
  // M1:   2
  // H1:   2
  // c1:   3
  // Ca1:  4
  //
  // set up jacStamp
  if( jacStamp.empty() )       // redundant as jacStamp is not static for this device
  {                            // it can't be as each cable may have a different number of nodes
    jacStamp.resize(numVars);
    jacStamp[0].resize(2);
    jacStamp[0][0] = 0;
    jacStamp[0][1] = 1;
    jacStamp[1].resize(2);
    jacStamp[1][0] = 1;
    jacStamp[1][1] = numVars - 10;
    for( int i=2; i<numVars; i+=10)
    {
      jacStamp[i].resize(11);
      if( i == 2 )
      {
        jacStamp[i][0] = 0;     // v_in
      }
      else
      {
        jacStamp[i][0] = i-10;  // v_prev
      }
      jacStamp[i][1] = i;    // v
      jacStamp[i][2] = i+1;  // n
      jacStamp[i][3] = i+2;  // m
      jacStamp[i][4] = i+3;  // h
      jacStamp[i][5] = i+4;  // a
      jacStamp[i][6] = i+5;  // b
      jacStamp[i][7] = i+6;  // M
      jacStamp[i][8] = i+7;  // H
      jacStamp[i][9] = i+8;  // c
      if( i==(numVars-10) )
      {
        jacStamp[i][10] = 1;  // v_out
      }
      else
      {
        jacStamp[i][10] = i+10; // v_next
      }
      jacStamp[i+1].resize(2);   // n
      jacStamp[i+1][0] = i;
      jacStamp[i+1][1] = i+1;
      jacStamp[i+2].resize(2);   // m
      jacStamp[i+2][0] = i;
      jacStamp[i+2][1] = i+2;
      jacStamp[i+3].resize(2);   // h
      jacStamp[i+3][0] = i;
      jacStamp[i+3][1] = i+3;
      jacStamp[i+4].resize(2);   // a
      jacStamp[i+4][0] = i;
      jacStamp[i+4][1] = i+4;
      jacStamp[i+5].resize(2);   // b
      jacStamp[i+5][0] = i;
      jacStamp[i+5][1] = i+5;
      jacStamp[i+6].resize(2);   // M
      jacStamp[i+6][0] = i;
      jacStamp[i+6][1] = i+6;
      jacStamp[i+7].resize(2);   // H
      jacStamp[i+7][0] = i;
      jacStamp[i+7][1] = i+7;
      jacStamp[i+8].resize(3);   // c
      jacStamp[i+8][0] = i;
      jacStamp[i+8][1] = i+8;
      jacStamp[i+8][2] = i+9;
      jacStamp[i+9].resize(4);   // ca
      jacStamp[i+9][0] = i;
      jacStamp[i+9][1] = i+6;
      jacStamp[i+9][2] = i+7;
      jacStamp[i+9][3] = i+9;

    }

  }

  /*
  // print out jacStamp
  int numRows = jacStamp.size();
  for( int i=0; i< numRows; i++ )
  {
  int numCol = jacStamp[i].size();
  Xyce::dout() << "jacStamp[ " << i << " ] = { ";
  for(int j=0; j<numCol; j++)
  {
  Xyce::dout() << jacStamp[i][j] << "  ";
  }
  Xyce::dout() << " } " <<  std::endl;
  }
  Xyce::dout() << std::endl;
  */

  // calculate segment conductivities used in load calls:
  // Note: conductivity to the previous and next segment is not symmetric if the segment length and/or radius is not are not equal
  // the full formula is:
  //   g(n,n') = radius * (radius')^2 / ( rInt segLength * ( segLength * (radius')^2 + segLength' * (radius)^2 ) )
  //
  // equation 6.30 Theoretical neuroscience: computational and mathematical modeling of neural systems, Peter Dayan and L. F. Abbot 2001
  // If we allow variable segment lengths and radii then we'll need to update this calculation
  //
  // Note: the above equation has the wrong units unless rInt (which is Ohm Length) is rLong (which is Ohm/Length).
  //
  gForward.resize(nSeg);
  gBackward.resize(nSeg);
  double segLength = length / nSeg;
  double rLong = rInt / (M_PI * radius * radius);  // longitudinal resistivity (ohm/length)
  gBackward[0] = radius * (radiusPrevious * radiusPrevious) / (rLong * segLength * ( segLength * radiusPrevious * radiusPrevious + lengthPrevious * radius * radius ));
  gForward[0] = radius * (radius * radius) / (rLong * segLength * ( segLength * radius * radius + segLength * radius * radius ));
  gBackward[nSeg-1] = radius * (radius * radius) / (rLong * segLength * ( segLength * radius * radius + segLength * radius * radius ));
  gForward[nSeg-1] = radius * (radiusNext * radiusNext) / (rLong * segLength * ( segLength * radiusNext * radiusNext + lengthNext * radius * radius ));
  for(int i=1; i<(nSeg-1); i++)
  {
    gBackward[i] = radius * (radius * radius) / (rLong * segLength * ( segLength * radius * radius + segLength * radius * radius ));
    gForward[i] = gBackward[i];
  }
  // gForward.resize(nSeg);
  // gBackward.resize(nSeg);
  // double segLength = length / nSeg;
  // gBackward[0] = radius * (radiusPrevious * radiusPrevious) / (rInt * segLength * ( segLength * radiusPrevious * radiusPrevious + lengthPrevious * radius * radius ));
  // gForward[0] = radius * (radius * radius) / (rInt * segLength * ( segLength * radius * radius + segLength * radius * radius ));
  // gBackward[nSeg-1] = radius * (radius * radius) / (rInt * segLength * ( segLength * radius * radius + segLength * radius * radius ));
  // gForward[nSeg-1] = radius * (radiusNext * radiusNext) / (rInt * segLength * ( segLength * radiusNext * radiusNext + lengthNext * radius * radius ));
  // for(int i=1; i<(nSeg-1); i++)
  // {
  //   gBackward[i] = radius * (radius * radius) / (rInt * segLength * ( segLength * radius * radius + segLength * radius * radius ));
  //   gForward[i] = gBackward[i];
  // }

  // allocate space for load and jacobian terms per segment
  // variable indecies loads
  li_Vol.resize(nSeg);
  li_nPro.resize(nSeg);
  li_mPro.resize(nSeg);
  li_hPro.resize(nSeg);
  li_aPro.resize(nSeg);
  li_bPro.resize(nSeg);
  li_MPro.resize(nSeg);
  li_HPro.resize(nSeg);
  li_cPro.resize(nSeg);
  li_CaPro.resize(nSeg);
  li_KCurrentState.resize(nSeg);
  li_NaCurrentState.resize(nSeg);
  segFvalue.resize(nSeg);
  segQvalue.resize(nSeg);
  segNEquFvalue.resize(nSeg);
  segNEquQvalue.resize(nSeg);
  segMEquFvalue.resize(nSeg);
  segMEquQvalue.resize(nSeg);
  segHEquFvalue.resize(nSeg);
  segHEquQvalue.resize(nSeg);
  segAEquFvalue.resize(nSeg);
  segAEquQvalue.resize(nSeg);
  segBEquFvalue.resize(nSeg);
  segBEquQvalue.resize(nSeg);
  segM_EquFvalue.resize(nSeg);
  segM_EquQvalue.resize(nSeg);
  segH_EquFvalue.resize(nSeg);
  segH_EquQvalue.resize(nSeg);
  segCEquFvalue.resize(nSeg);
  segCEquQvalue.resize(nSeg);
  segCaEquFvalue.resize(nSeg);
  segCaEquQvalue.resize(nSeg);

  // jacobian elements
  segF_dVp.resize(nSeg);
  segF_dV.resize(nSeg);
  segF_dVn.resize(nSeg);
  segF_dn.resize(nSeg);
  segF_dm.resize(nSeg);
  segF_dh.resize(nSeg);
  segF_da.resize(nSeg);
  segF_db.resize(nSeg);
  segF_dM.resize(nSeg);
  segF_dH.resize(nSeg);
  segF_dc.resize(nSeg);
  segQ_dV.resize(nSeg);
  dnF_dV.resize(nSeg);
  dnF_dn.resize(nSeg);
  dnQ_dn.resize(nSeg);
  dmF_dV.resize(nSeg);
  dmF_dm.resize(nSeg);
  dmQ_dm.resize(nSeg);
  dhF_dV.resize(nSeg);
  dhF_dh.resize(nSeg);
  dhQ_dh.resize(nSeg);
  daF_dV.resize(nSeg);
  daF_da.resize(nSeg);
  daQ_da.resize(nSeg);
  dbF_dV.resize(nSeg);
  dbF_db.resize(nSeg);
  dbQ_db.resize(nSeg);
  dMF_dV.resize(nSeg);
  dMF_dM.resize(nSeg);
  dMQ_dM.resize(nSeg);
  dHF_dV.resize(nSeg);
  dHF_dH.resize(nSeg);
  dHQ_dH.resize(nSeg);
  dcF_dV.resize(nSeg);
  dcF_dc.resize(nSeg);
  dcF_dCa.resize(nSeg);
  dcQ_dc.resize(nSeg);
  dCaF_dV.resize(nSeg);
  dCaF_dM.resize(nSeg);
  dCaF_dH.resize(nSeg);
  dCaF_dCa.resize(nSeg);
  dCaQ_dCa.resize(nSeg);

  // state vars
  potassiumCurrent.resize(nSeg);
  sodiumCurrent.resize(nSeg);

}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
Instance::~Instance()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::registerLIDs(const std::vector<int> & intLIDVecRef,
                            const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  Instance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  li_Pos  = extLIDVec[0];
  li_Neg  = extLIDVec[1];
  for( int i=0, j=0; i<nSeg; i++, j+=10)
  {
    li_Vol[i]  = intLIDVec[j];
    li_nPro[i] = intLIDVec[j+1];
    li_mPro[i] = intLIDVec[j+2];
    li_hPro[i] = intLIDVec[j+3];
    li_aPro[i] = intLIDVec[j+4];
    li_bPro[i] = intLIDVec[j+5];
    li_MPro[i] = intLIDVec[j+6];
    li_HPro[i] = intLIDVec[j+7];
    li_cPro[i] = intLIDVec[j+8];
    li_CaPro[i] = intLIDVec[j+9];
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl
         << "  li_Neg = " << li_Neg << std::endl;
    for( int i=0; i<nSeg; i++ )
    {
      Xyce::dout() << "  li_Vol[ " << i << " ] = " << li_Vol[i] << std::endl
           << "  li_nPro[ " << i << " ] = " << li_nPro[i] << std::endl
           << "  li_mPro[ " << i << " ] = " << li_mPro[i] << std::endl
           << "  li_hPro[ " << i << " ] = " << li_hPro[i] << std::endl
           << "  li_aPro[ " << i << " ] = " << li_aPro[i] << std::endl
           << "  li_bPro[ " << i << " ] = " << li_bPro[i] << std::endl
           << "  li_MPro[ " << i << " ] = " << li_MPro[i] << std::endl
           << "  li_HPro[ " << i << " ] = " << li_HPro[i] << std::endl
           << "  li_cPro[ " << i << " ] = " << li_cPro[i] << std::endl
           << "  li_CaPro[ " << i << " ] = " << li_CaPro[i] << std::endl;
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
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
  for (int i = 0; i < nSeg; i++)
  {
    std::ostringstream segNumber;
    segNumber << i;
    addInternalNode(symbol_table, li_Vol[i], getName(), "V" + segNumber.str());
    addInternalNode(symbol_table, li_nPro[i], getName(), "N" + segNumber.str());
    addInternalNode(symbol_table, li_mPro[i], getName(), "M" + segNumber.str());
    addInternalNode(symbol_table, li_hPro[i], getName(), "H" + segNumber.str());
    addInternalNode(symbol_table, li_aPro[i], getName(), "A" + segNumber.str());
    addInternalNode(symbol_table, li_bPro[i], getName(), "B" + segNumber.str());
    addInternalNode(symbol_table, li_MPro[i], getName(), "M_" + segNumber.str());
    addInternalNode(symbol_table, li_HPro[i], getName(), "H_" + segNumber.str());
    addInternalNode(symbol_table, li_cPro[i], getName(), "C" + segNumber.str());
    addInternalNode(symbol_table, li_CaPro[i], getName(), "Ca" + segNumber.str());
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // copy over the global ID lists.
  staLIDVec = staLIDVecRef;

  for( int i=0, j=0; i<nSeg; i++, j+=2)
  {
    li_KCurrentState[i] = staLIDVec[j];
    li_NaCurrentState[i] = staLIDVec[j+1];
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
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
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  // external terminals
  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNextNodeOffset = jacLIDVec[0][1];
  ANegEquNegNodeOffset = jacLIDVec[1][0];
  ANegEquLastNodeOffset = jacLIDVec[1][1];
  /*
    Xyce::dout() << "APosEquPosNodeOffset = " << APosEquPosNodeOffset << std::endl;
    Xyce::dout() << "APosEquNextNodeOffset = " << APosEquNextNodeOffset << std::endl;
    Xyce::dout() << "ANegEquNegNodeOffset = " << ANegEquNegNodeOffset << std::endl;
    Xyce::dout() << "ANegEquLastNodeOffset = " << ANegEquLastNodeOffset << std::endl;
  */

  // internal variables
  SegVEqnVpreOffset.resize(nSeg);
  SegVEqnVsegOffset.resize(nSeg);
  SegVEqnVnexOffset.resize(nSeg);
  SegVEqnNOffset.resize(nSeg);
  SegVEqnMOffset.resize(nSeg);
  SegVEqnHOffset.resize(nSeg);
  NEquVNodeOffset.resize(nSeg);
  NEquNNodeOffset.resize(nSeg);
  MEquVNodeOffset.resize(nSeg);
  MEquMNodeOffset.resize(nSeg);
  HEquVNodeOffset.resize(nSeg);
  HEquHNodeOffset.resize(nSeg);
  AEquVNodeOffset.resize(nSeg);
  AEquANodeOffset.resize(nSeg);
  BEquVNodeOffset.resize(nSeg);
  BEquBNodeOffset.resize(nSeg);
  M_EquVNodeOffset.resize(nSeg);
  M_EquM_NodeOffset.resize(nSeg);
  H_EquVNodeOffset.resize(nSeg);
  H_EquH_NodeOffset.resize(nSeg);
  CEquVNodeOffset.resize(nSeg);
  CEquCNodeOffset.resize(nSeg);
  CEquCaNodeOffset.resize(nSeg);
  CaEquVNodeOffset.resize(nSeg);
  CaEquM_NodeOffset.resize(nSeg);
  CaEquH_NodeOffset.resize(nSeg);
  CaEquCaNodeOffset.resize(nSeg);

  for(int i=0, j=2; i<nSeg; i++, j+=10 )
  {
    SegVEqnVpreOffset[i] = jacLIDVec[j][0];
    SegVEqnVsegOffset[i] = jacLIDVec[j][1];
    SegVEqnNOffset[i] = jacLIDVec[j][2];
    SegVEqnMOffset[i] = jacLIDVec[j][3];
    SegVEqnHOffset[i] = jacLIDVec[j][4];
    SegVEqnVnexOffset[i] = jacLIDVec[j][5];

    NEquVNodeOffset[i] = jacLIDVec[j+1][0];
    NEquNNodeOffset[i] = jacLIDVec[j+1][1];
    MEquVNodeOffset[i] = jacLIDVec[j+2][0];
    MEquMNodeOffset[i] = jacLIDVec[j+2][1];
    HEquVNodeOffset[i] = jacLIDVec[j+3][0];
    HEquHNodeOffset[i] = jacLIDVec[j+3][1];
    AEquVNodeOffset[i] = jacLIDVec[j+4][0];
    AEquANodeOffset[i] = jacLIDVec[j+4][1];
    BEquVNodeOffset[i] = jacLIDVec[j+5][0];
    BEquBNodeOffset[i] = jacLIDVec[j+5][1];
    M_EquVNodeOffset[i] = jacLIDVec[j+6][0];
    M_EquM_NodeOffset[i] = jacLIDVec[j+6][1];
    H_EquVNodeOffset[i] = jacLIDVec[j+7][0];
    H_EquH_NodeOffset[i] = jacLIDVec[j+7][1];
    CEquVNodeOffset[i] = jacLIDVec[j+8][0];
    CEquCNodeOffset[i] = jacLIDVec[j+8][1];
    CEquCaNodeOffset[i] = jacLIDVec[j+8][2];
    CaEquVNodeOffset[i] = jacLIDVec[j+9][0];
    CaEquM_NodeOffset[i] = jacLIDVec[j+9][1];
    CaEquH_NodeOffset[i] = jacLIDVec[j+9][2];
    CaEquCaNodeOffset[i] = jacLIDVec[j+9][3];

    /*
      Xyce::dout() <<  SegVEqnVpreOffset[i] << ", "
      << SegVEqnVsegOffset[i] << ", "
      << SegVEqnNOffset[i] << ", "
      << SegVEqnMOffset[i] << ", "
      << SegVEqnHOffset[i] << ", "
      << SegVEqnVnexOffset[i] << ", "
      << NEquVNodeOffset[i] << ", "
      << NEquNNodeOffset[i] << ", "
      << MEquVNodeOffset[i] << ", "
      << MEquMNodeOffset[i] << ", "
      << HEquVNodeOffset[i] << ", "
      << HEquHNodeOffset[i] << std::endl;
    */
  }
}


//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;

  // here we take the current solutions for V1, V2, n, m and h
  // and use those to calculate all the terms needed for the next
  // load cycle (F, Q, dFdX, dQdX)
  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;

  double vIn = (*solVectorPtr)[li_Pos];
  double vOut = (*solVectorPtr)[li_Neg];

  // take care of the input and output nodes as they are different
  kcl1Fvalue = gForward[0] * ((*solVectorPtr)[li_Vol[0]] - vIn );
  dkcl1F_dVin = -gForward[0];
  dkcl1F_dVs0 =  gForward[0];
  kcl2Fvalue = gBackward[nSeg-1] * ((*solVectorPtr)[li_Vol[nSeg-1]] - vOut );
  dkcl2F_dVout = -gBackward[nSeg-1];
  dkcl2F_dVsn = gBackward[nSeg-1];

  // loop over segments getting all the load and jacobian terms for each segment
  for( int i=0; i<nSeg; i++ )
  {
    // for this segment get the values of the local vars
    double vSeg  = (*solVectorPtr)[li_Vol[i]];
    double vNext = 0.0;
    if (i == (nSeg - 1))
    {
      vNext = vOut;
    }
    else
    {
      vNext = (*solVectorPtr)[li_Vol[i+1]];
    }
    double vPrev = 0.0;
    if (i == 0 )
    {
      vPrev = vIn;
    }
    else
    {
      vPrev = (*solVectorPtr)[li_Vol[i-1]];
    }
    double nVarSeg = (*solVectorPtr)[li_nPro[i]];
    double mVarSeg = (*solVectorPtr)[li_mPro[i]];
    double hVarSeg = (*solVectorPtr)[li_hPro[i]];
    double aVarSeg = (*solVectorPtr)[li_aPro[i]];
    double bVarSeg = (*solVectorPtr)[li_bPro[i]];
    double M_VarSeg = (*solVectorPtr)[li_MPro[i]];
    double H_VarSeg = (*solVectorPtr)[li_HPro[i]];
    double cVarSeg = (*solVectorPtr)[li_cPro[i]];
    double CaVarSeg = (*solVectorPtr)[li_CaPro[i]];

    // do the voltage equation for this node
    // get function and derivative values as we go.
    {
      // F part
      // use scoping to avoid lots of similar variable names
      const int numDeriv = 11;
      Sacado::Fad::SFad<double,11> vVar ( numDeriv,  0, vSeg );
      Sacado::Fad::SFad<double,11> vVpr ( numDeriv,  1, vPrev );
      Sacado::Fad::SFad<double,11> vVne ( numDeriv,  2, vNext );
      Sacado::Fad::SFad<double,11> nVar ( numDeriv,  3, nVarSeg );
      Sacado::Fad::SFad<double,11> mVar ( numDeriv,  4, mVarSeg );
      Sacado::Fad::SFad<double,11> hVar ( numDeriv,  5, hVarSeg );
      Sacado::Fad::SFad<double,11> aVar ( numDeriv,  6, aVarSeg );
      Sacado::Fad::SFad<double,11> bVar ( numDeriv,  7, bVarSeg );
      Sacado::Fad::SFad<double,11> M_Var( numDeriv,  8, M_VarSeg );
      Sacado::Fad::SFad<double,11> H_Var( numDeriv,  9, H_VarSeg );
      Sacado::Fad::SFad<double,11> cVar ( numDeriv, 10, cVarSeg );

      // parameters
      Sacado::Fad::SFad<double,11> gPrev( gBackward[i] );
      Sacado::Fad::SFad<double,11> gNext( gForward[i] );
      Sacado::Fad::SFad<double,11> gMemVar( model_.gMem * segArea );
      Sacado::Fad::SFad<double,11> vRestVar( model_.vRest );
      Sacado::Fad::SFad<double,11> gKVar( model_.gK * segArea );
      Sacado::Fad::SFad<double,11> eKVar( model_.eK );
      Sacado::Fad::SFad<double,11> gNaVar( model_.gNa * segArea );
      Sacado::Fad::SFad<double,11> eNaVar( model_.eNa );
      Sacado::Fad::SFad<double,11> gAVar( model_.gA * segArea );
      Sacado::Fad::SFad<double,11> eAVar( model_.eA );
      Sacado::Fad::SFad<double,11> gCaVar( model_.gCa * segArea );
      Sacado::Fad::SFad<double,11> eCaVar( model_.eCa );
      Sacado::Fad::SFad<double,11> gKCaVar( model_.gKCa * segArea );
      Sacado::Fad::SFad<double,11> CaInitVar( model_.CaInit );
      Sacado::Fad::SFad<double,11> CaGammaVar( model_.CaGamma );
      Sacado::Fad::SFad<double,11> CaTauVar( model_.CaTau );

      // compute the vaud and derivative terms for KCL 1 F
      Sacado::Fad::SFad<double,11> resultFad;
      resultFad = kcl1EquF( vVar, vVpr, vVne, nVar, mVar, hVar, aVar, bVar, M_Var, H_Var, cVar, gPrev, gNext, gMemVar, vRestVar, gKVar, eKVar, gNaVar, eNaVar, gAVar, eAVar, gCaVar, eCaVar, gKCaVar);

      segFvalue[i] = resultFad.val();
      segF_dV[i]   = resultFad.dx(0);
      segF_dVp[i]  = resultFad.dx(1);
      segF_dVn[i]  = resultFad.dx(2);
      segF_dn[i]   = resultFad.dx(3);
      segF_dm[i]   = resultFad.dx(4);
      segF_dh[i]   = resultFad.dx(5);
      segF_da[i]   = resultFad.dx(6);
      segF_db[i]   = resultFad.dx(7);
      segF_dM[i]   = resultFad.dx(8);
      segF_dH[i]   = resultFad.dx(9);
      segF_dc[i]   = resultFad.dx(10);
    }
    {
      // Q part
      const int numDeriv = 1;
      Sacado::Fad::SFad<double,1> vVar( numDeriv, 0, vSeg );

      // parameters
      Sacado::Fad::SFad<double,1> cMemVar( model_.cMem * segArea );

      Sacado::Fad::SFad<double,1> resultFad;
      resultFad    = kcl1EquQ( vVar, cMemVar );
      segQvalue[i] = resultFad.val();
      segQ_dV[i]   = resultFad.dx(0);

    }

    // n - equation
    {
      const int numDeriv = 2;
      Sacado::Fad::SFad<double,2> vVar( numDeriv, 0, vSeg );
      Sacado::Fad::SFad<double,2> nVar( numDeriv, 1, nVarSeg );
      // parameter
      Sacado::Fad::SFad<double,2> vRestVar( model_.vRest );

      Sacado::Fad::SFad<double,2> resultFad = nEquF( vVar, nVar, vRestVar);
      segNEquFvalue[i] = resultFad.val();
      dnF_dV[i]        = resultFad.dx(0);
      dnF_dn[i]        = resultFad.dx(1);
    }
    {
      const int numDeriv = 1;
      Sacado::Fad::SFad<double,1> nVar( numDeriv, 0, nVarSeg );

      Sacado::Fad::SFad<double,1> resultFad = nEquQ( nVar );
      segNEquQvalue[i] = resultFad.val();
      dnQ_dn[i]        = resultFad.dx(0);
    }

    // m - equation
    {
      const int numDeriv = 2;
      Sacado::Fad::SFad<double,2> vVar( numDeriv, 0, vSeg );
      Sacado::Fad::SFad<double,2> mVar( numDeriv, 1, mVarSeg );
      // parameter
      Sacado::Fad::SFad<double,2> vRestVar( model_.vRest );

      Sacado::Fad::SFad<double,2> resultFad = mEquF( vVar, mVar, vRestVar );
      segMEquFvalue[i] = resultFad.val();
      dmF_dV[i]        = resultFad.dx(0);
      dmF_dm[i]        = resultFad.dx(1);
    }
    {
      const int numDeriv = 1;
      Sacado::Fad::SFad<double,1> mVar( numDeriv, 0, mVarSeg );

      Sacado::Fad::SFad<double,1> resultFad = mEquQ( mVar );
      segMEquQvalue[i] = resultFad.val();
      dmQ_dm[i]        = resultFad.dx(0);
    }

    // h - equation
    {
      const int numDeriv = 2;
      Sacado::Fad::SFad<double,2> vVar( numDeriv, 0, vSeg );
      Sacado::Fad::SFad<double,2> hVar( numDeriv, 1, hVarSeg );
      // parameter
      Sacado::Fad::SFad<double,2> vRestVar( model_.vRest );

      Sacado::Fad::SFad<double,2> resultFad = hEquF( vVar, hVar, vRestVar );
      segHEquFvalue[i] = resultFad.val();
      dhF_dV[i]        = resultFad.dx(0);
      dhF_dh[i]        = resultFad.dx(1);
    }
    {
      const int numDeriv = 1;
      Sacado::Fad::SFad<double,1> hVar( numDeriv, 0, hVarSeg );

      Sacado::Fad::SFad<double,1> resultFad = hEquQ( hVar );
      segHEquQvalue[i] = resultFad.val();
      dhQ_dh[i]        = resultFad.dx(0);
    }

    // a - equation
    {
      const int numDeriv = 2;
      Sacado::Fad::SFad<double,2> vVar( numDeriv, 0, vSeg );
      Sacado::Fad::SFad<double,2> aVar( numDeriv, 1, aVarSeg );
      // parameter
      Sacado::Fad::SFad<double,2> vRestVar( model_.vRest );

      Sacado::Fad::SFad<double,2> resultFad = aEquF( vVar, aVar, vRestVar );
      segAEquFvalue[i] = resultFad.val();
      daF_dV[i]        = resultFad.dx(0);
      daF_da[i]        = resultFad.dx(1);
    }
    {
      const int numDeriv = 1;
      Sacado::Fad::SFad<double,1> aVar( numDeriv, 0, aVarSeg );

      Sacado::Fad::SFad<double,1> resultFad = aEquQ( aVar );
      segAEquQvalue[i] = resultFad.val();
      daQ_da[i]        = resultFad.dx(0);
    }

    // b - equation
    {
      const int numDeriv = 2;
      Sacado::Fad::SFad<double,2> vVar( numDeriv, 0, vSeg );
      Sacado::Fad::SFad<double,2> bVar( numDeriv, 1, bVarSeg );
      // parameter
      Sacado::Fad::SFad<double,2> vRestVar( model_.vRest );

      Sacado::Fad::SFad<double,2> resultFad = bEquF( vVar, bVar, vRestVar );
      segBEquFvalue[i] = resultFad.val();
      dbF_dV[i]        = resultFad.dx(0);
      dbF_db[i]        = resultFad.dx(1);
    }
    {
      const int numDeriv = 1;
      Sacado::Fad::SFad<double,1> bVar( numDeriv, 0, bVarSeg );

      Sacado::Fad::SFad<double,1> resultFad = bEquQ( bVar );
      segBEquQvalue[i] = resultFad.val();
      dbQ_db[i]        = resultFad.dx(0);
    }

    // M - equation
    {
      const int numDeriv = 2;
      Sacado::Fad::SFad<double,2> vVar( numDeriv, 0, vSeg );
      Sacado::Fad::SFad<double,2> M_Var( numDeriv, 1, M_VarSeg );
      // parameter
      Sacado::Fad::SFad<double,2> vRestVar( model_.vRest );

      Sacado::Fad::SFad<double,2> resultFad = M_EquF( vVar, M_Var, vRestVar );
      segM_EquFvalue[i] = resultFad.val();
      dMF_dV[i]        = resultFad.dx(0);
      dMF_dM[i]        = resultFad.dx(1);
    }
    {
      const int numDeriv = 1;
      Sacado::Fad::SFad<double,1> M_Var( numDeriv, 0, M_VarSeg );

      Sacado::Fad::SFad<double,1> resultFad = M_EquQ( M_Var );
      segM_EquQvalue[i] = resultFad.val();
      dMQ_dM[i]        = resultFad.dx(0);
    }

    // H - equation
    {
      const int numDeriv = 2;
      Sacado::Fad::SFad<double,2> vVar( numDeriv, 0, vSeg );
      Sacado::Fad::SFad<double,2> H_Var( numDeriv, 1, H_VarSeg );
      // parameter
      Sacado::Fad::SFad<double,2> vRestVar( model_.vRest );

      Sacado::Fad::SFad<double,2> resultFad = H_EquF( vVar, H_Var, vRestVar );
      segH_EquFvalue[i] = resultFad.val();
      dHF_dV[i]        = resultFad.dx(0);
      dHF_dH[i]        = resultFad.dx(1);
    }
    {
      const int numDeriv = 1;
      Sacado::Fad::SFad<double,1> H_Var( numDeriv, 0, H_VarSeg );

      Sacado::Fad::SFad<double,1> resultFad = H_EquQ( H_Var );
      segH_EquQvalue[i] = resultFad.val();
      dHQ_dH[i]        = resultFad.dx(0);
    }

    // c - equation
    {
      const int numDeriv = 3;
      Sacado::Fad::SFad<double,3> vVar( numDeriv, 0, vSeg );
      Sacado::Fad::SFad<double,3> cVar( numDeriv, 1, cVarSeg );
      Sacado::Fad::SFad<double,3> CaVar( numDeriv, 2, CaVarSeg );
      // parameter
      Sacado::Fad::SFad<double,3> vRestVar( model_.vRest );

      Sacado::Fad::SFad<double,3> resultFad = C_EquF( vVar, cVar, CaVar, vRestVar );
      segCEquFvalue[i] = resultFad.val();
      dcF_dV[i]        = resultFad.dx(0);
      dcF_dc[i]        = resultFad.dx(1);
      dcF_dCa[i]       = resultFad.dx(2);
    }
    {
      const int numDeriv = 1;
      Sacado::Fad::SFad<double,1> cVar( numDeriv, 0, cVarSeg );

      Sacado::Fad::SFad<double,1> resultFad = C_EquQ( cVar );
      segCEquQvalue[i] = resultFad.val();
      dcQ_dc[i]        = resultFad.dx(0);
    }

    // Ca - equation
    {
      const int numDeriv = 4;
      Sacado::Fad::SFad<double,4> vVar( numDeriv, 0, vSeg );
      Sacado::Fad::SFad<double,4> M_Var( numDeriv, 1, M_VarSeg );
      Sacado::Fad::SFad<double,4> H_Var( numDeriv, 2, H_VarSeg );
      Sacado::Fad::SFad<double,4> CaVar( numDeriv, 3, CaVarSeg );
      // parameter
      Sacado::Fad::SFad<double,4> gCaVar( model_.gCa );
      Sacado::Fad::SFad<double,4> eCaVar( model_.gCa );
      Sacado::Fad::SFad<double,4> CaGammaVar( model_.CaGamma );
      Sacado::Fad::SFad<double,4> CaTauVar( model_.CaTau );

      Sacado::Fad::SFad<double,4> resultFad = Ca_EquF( vVar, M_Var, H_Var, CaVar, gCaVar, eCaVar, CaGammaVar, CaTauVar );
      segCaEquFvalue[i] = resultFad.val();
      dCaF_dV[i]        = resultFad.dx(0);
      dCaF_dM[i]        = resultFad.dx(1);
      dCaF_dH[i]        = resultFad.dx(2);
      dCaF_dCa[i]       = resultFad.dx(3);
    }
    {
      const int numDeriv = 1;
      Sacado::Fad::SFad<double,1> CaVar( numDeriv, 0, CaVarSeg );

      Sacado::Fad::SFad<double,1> resultFad = Ca_EquQ( CaVar );
      segCaEquQvalue[i] = resultFad.val();
      dCaQ_dCa[i]       = resultFad.dx(0);
    }


  }

#if 0
  //Xyce::dout() << "Instance::updateIntermediateVars()" << std::endl;

  bool bsuccess = true;

  // here we take the current solutions for V1, V2, n, m and h
  // and use those to calculate all the terms needed for the next
  // load cycle (F, Q, dFdX, dQdX)

  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;

  // use suffix "now" to to clarify that this for the latest solution
  double v1Now  = (*solVectorPtr)[li_Pos];
  double v2Now  = (*solVectorPtr)[li_Neg];
  double nNow   = (*solVectorPtr)[li_nPro];
  double mNow   = (*solVectorPtr)[li_mPro];
  double hNow   = (*solVectorPtr)[li_hPro];
  double aNow   = (*solVectorPtr)[li_aPro];
  double bNow   = (*solVectorPtr)[li_bPro];
  double M_Now  = (*solVectorPtr)[li_M_Pro];
  double H_Now  = (*solVectorPtr)[li_H_Pro];
  double cNow   = (*solVectorPtr)[li_cPro];
  double CaNow  = (*solVectorPtr)[li_CaPro];

  // get function and derivative values
  // independent variables
  // use scoping to avoid lots of similar variable names
  {
    Sacado::Fad::SFad<double,10> v1Var( 10, 0, v1Now );
    Sacado::Fad::SFad<double,10> v2Var( 10, 1, v2Now );
    Sacado::Fad::SFad<double,10> nVar( 10, 2, nNow );
    Sacado::Fad::SFad<double,10> mVar( 10, 3, mNow );
    Sacado::Fad::SFad<double,10> hVar( 10, 4, hNow );
    Sacado::Fad::SFad<double,10> aVar( 10, 5, aNow );
    Sacado::Fad::SFad<double,10> bVar( 10, 6, bNow );
    Sacado::Fad::SFad<double,10> M_Var( 10, 7, M_Now );
    Sacado::Fad::SFad<double,10> H_Var( 10, 8, H_Now );
    Sacado::Fad::SFad<double,10> cVar( 10, 9, cNow );

    // parameters from the model that we'll need.
    Sacado::Fad::SFad<double,10> gMemVar( model_.gMem );
    Sacado::Fad::SFad<double,10> vRestVar( model_.vRest );
    Sacado::Fad::SFad<double,10> gKVar( model_.gK );
    Sacado::Fad::SFad<double,10> eKVar( model_.eK );
    Sacado::Fad::SFad<double,10> gNaVar( model_.gNa );
    Sacado::Fad::SFad<double,10> eNaVar( model_.eNa );
    Sacado::Fad::SFad<double,10> gAVar( model_.gA );
    Sacado::Fad::SFad<double,10> eAVar( model_.eA );
    Sacado::Fad::SFad<double,10> gCaVar( model_.gCa );
    Sacado::Fad::SFad<double,10> eCaVar( model_.gCa );
    Sacado::Fad::SFad<double,10> gKCaVar( model_.gKCa );
    Sacado::Fad::SFad<double,10> CaInitVar( model_.CaInit );
    Sacado::Fad::SFad<double,10> CaGammaVar( model_.CaGamma );
    Sacado::Fad::SFad<double,10> CaTauVar( model_.CaTau );

    // compute the vaud and derivative terms for KCL 1 F
    Sacado::Fad::SFad<double,10> resultFad;
    resultFad = kcl1EquF( v1Var, v2Var, nVar, mVar, hVar, aVar, bVar, M_Var, H_Var, cVar, gMemVar, vRestVar, gKVar, eKVar, gNaVar, eNaVar, gAVar, eAVar, gCaVar, eCaVar, gKCaVar);
    kcl1Fvalue = resultFad.val();
    dkcl1F_dV1 = resultFad.dx(0);
    dkcl1F_dV2 = resultFad.dx(1);
    dkcl1F_dn  = resultFad.dx(2);
    dkcl1F_dm  = resultFad.dx(3);
    dkcl1F_dh  = resultFad.dx(4);
    dkcl1F_da  = resultFad.dx(5);
    dkcl1F_db  = resultFad.dx(6);
    dkcl1F_dM  = resultFad.dx(7);
    dkcl1F_dH  = resultFad.dx(8);
    dkcl1F_dc  = resultFad.dx(9);

    // compute the vaud and derivative terms for KCL 2 F
    resultFad = kcl2EquF( v1Var, v2Var, nVar, mVar, hVar, aVar, bVar, M_Var, H_Var, cVar, gMemVar, vRestVar, gKVar, eKVar, gNaVar, eNaVar, gAVar, eAVar, gCaVar, eCaVar, gKCaVar);
    kcl2Fvalue = resultFad.val();
    dkcl2F_dV1 = resultFad.dx(0);
    dkcl2F_dV2 = resultFad.dx(1);
    dkcl2F_dn  = resultFad.dx(2);
    dkcl2F_dm  = resultFad.dx(3);
    dkcl2F_dh  = resultFad.dx(4);
    dkcl2F_da  = resultFad.dx(5);
    dkcl2F_db  = resultFad.dx(6);
    dkcl2F_dM  = resultFad.dx(7);
    dkcl2F_dH  = resultFad.dx(8);
    dkcl2F_dc  = resultFad.dx(9);
  }

  {
    Sacado::Fad::SFad<double,2> v1Var( 2, 0, v1Now );
    Sacado::Fad::SFad<double,2> v2Var( 2, 1, v2Now );

    // parameters
    Sacado::Fad::SFad<double,2> cMemVar( model_.cMem );

    Sacado::Fad::SFad<double,2> resultFad;
    resultFad = kcl1EquQ( v1Var, v2Var, cMemVar );
    kcl1Qvalue = resultFad.val();
    dkcl1Q_dV1 = resultFad.dx(0);
    dkcl1Q_dV2 = resultFad.dx(1);

    resultFad = kcl2EquQ( v1Var, v2Var, cMemVar );
    kcl2Qvalue = resultFad.val();
    dkcl2Q_dV1 = resultFad.dx(0);
    dkcl2Q_dV2 = resultFad.dx(1);
  }

  // n - equation
  {
    Sacado::Fad::SFad<double,2> v1Var( 2, 0, v1Now );
    Sacado::Fad::SFad<double,2> nVar( 2, 1, nNow );
    // parameter
    Sacado::Fad::SFad<double,2> vRestVar( model_.vRest );

    Sacado::Fad::SFad<double,2> resultFad = nEquF( v1Var, nVar, vRestVar);
    nEquFvalue = resultFad.val();
    dnF_dV1 = resultFad.dx(0);
    dnF_dn  = resultFad.dx(1);
  }

  {
    Sacado::Fad::SFad<double,1> nVar( 1, 0, nNow );

    Sacado::Fad::SFad<double,1> resultFad = nEquQ( nVar );
    nEquQvalue = resultFad.val();
    dnQ_dn     = resultFad.dx(0);
  }

  // m - equation
  {
    Sacado::Fad::SFad<double,2> v1Var( 2, 0, v1Now );
    Sacado::Fad::SFad<double,2> mVar( 2, 1, mNow );
    // parameter
    Sacado::Fad::SFad<double,2> vRestVar( model_.vRest );

    Sacado::Fad::SFad<double,2> resultFad = mEquF( v1Var, mVar, vRestVar );
    mEquFvalue = resultFad.val();
    dmF_dV1 = resultFad.dx(0);
    dmF_dm  = resultFad.dx(1);
  }
  {
    Sacado::Fad::SFad<double,1> mVar( 1, 0, mNow );

    Sacado::Fad::SFad<double,1> resultFad = mEquQ( mVar );
    mEquQvalue = resultFad.val();
    dmQ_dm     = resultFad.dx(0);
  }

  // h - equation
  {
    Sacado::Fad::SFad<double,2> v1Var( 2, 0, v1Now );
    Sacado::Fad::SFad<double,2> hVar( 2, 1, hNow );
    // parameter
    Sacado::Fad::SFad<double,2> vRestVar( model_.vRest );

    Sacado::Fad::SFad<double,2> resultFad = hEquF( v1Var, hVar, vRestVar );
    hEquFvalue = resultFad.val();
    dhF_dV1 = resultFad.dx(0);
    dhF_dh  = resultFad.dx(1);
  }
  {
    Sacado::Fad::SFad<double,1> hVar( 1, 0, hNow );

    Sacado::Fad::SFad<double,1> resultFad = hEquQ( hVar );
    hEquQvalue = resultFad.val();
    dhQ_dh     = resultFad.dx(0);
  }

  // a - equation
  {
    Sacado::Fad::SFad<double,2> v1Var( 2, 0, v1Now );
    Sacado::Fad::SFad<double,2> aVar( 2, 1, aNow );
    // parameter
    Sacado::Fad::SFad<double,2> vRestVar( model_.vRest );

    Sacado::Fad::SFad<double,2> resultFad = aEquF( v1Var, aVar, vRestVar );
    aEquFvalue = resultFad.val();
    daF_dV1 = resultFad.dx(0);
    daF_da  = resultFad.dx(1);
  }
  {
    Sacado::Fad::SFad<double,1> aVar( 1, 0, aNow );

    Sacado::Fad::SFad<double,1> resultFad = aEquQ( aVar );
    aEquQvalue = resultFad.val();
    daQ_da     = resultFad.dx(0);
  }

  // b - equation
  {
    Sacado::Fad::SFad<double,2> v1Var( 2, 0, v1Now );
    Sacado::Fad::SFad<double,2> bVar( 2, 1, bNow );
    // parameter
    Sacado::Fad::SFad<double,2> vRestVar( model_.vRest );

    Sacado::Fad::SFad<double,2> resultFad = bEquF( v1Var, bVar, vRestVar );
    bEquFvalue = resultFad.val();
    dbF_dV1 = resultFad.dx(0);
    dbF_db  = resultFad.dx(1);
  }
  {
    Sacado::Fad::SFad<double,1> bVar( 1, 0, bNow );

    Sacado::Fad::SFad<double,1> resultFad = aEquQ( bVar );
    bEquQvalue = resultFad.val();
    dbQ_db     = resultFad.dx(0);
  }

  // M - equation
  {
    Sacado::Fad::SFad<double,2> v1Var( 2, 0, v1Now );
    Sacado::Fad::SFad<double,2> M_Var( 2, 1, M_Now );
    // parameter
    Sacado::Fad::SFad<double,2> vRestVar( model_.vRest );

    Sacado::Fad::SFad<double,2> resultFad = M_EquF( v1Var, M_Var, vRestVar );
    M_EquFvalue = resultFad.val();
    dMF_dV1 = resultFad.dx(0);
    dMF_dM  = resultFad.dx(1);
  }
  {
    Sacado::Fad::SFad<double,1> M_Var( 1, 0, M_Now );

    Sacado::Fad::SFad<double,1> resultFad = aEquQ( M_Var );
    M_EquQvalue = resultFad.val();
    dMQ_dM     = resultFad.dx(0);
  }

  // H - equation
  {
    Sacado::Fad::SFad<double,2> v1Var( 2, 0, v1Now );
    Sacado::Fad::SFad<double,2> H_Var( 2, 1, H_Now );
    // parameter
    Sacado::Fad::SFad<double,2> vRestVar( model_.vRest );

    Sacado::Fad::SFad<double,2> resultFad = H_EquF( v1Var, H_Var, vRestVar );
    H_EquFvalue = resultFad.val();
    dHF_dV1 = resultFad.dx(0);
    dHF_dH  = resultFad.dx(1);
  }
  {
    Sacado::Fad::SFad<double,1> H_Var( 1, 0, H_Now );

    Sacado::Fad::SFad<double,1> resultFad = H_EquQ( H_Var );
    H_EquQvalue = resultFad.val();
    dHQ_dH     = resultFad.dx(0);
  }

  // c - equation
  {
    Sacado::Fad::SFad<double,3> v1Var( 3, 0, v1Now );
    Sacado::Fad::SFad<double,3> cVar( 3, 1, cNow );
    Sacado::Fad::SFad<double,3> CaVar( 3, 2, CaNow );
    // parameter
    Sacado::Fad::SFad<double,3> vRestVar( model_.vRest );

    Sacado::Fad::SFad<double,3> resultFad = C_EquF( v1Var, cVar, CaVar, vRestVar );
    cEquFvalue = resultFad.val();
    dcF_dV1 = resultFad.dx(0);
    dcF_dc  = resultFad.dx(1);
    dcF_dCa = resultFad.dx(2);
  }
  {
    Sacado::Fad::SFad<double,1> cVar( 1, 0, cNow );

    Sacado::Fad::SFad<double,1> resultFad = C_EquQ( cVar );
    cEquQvalue = resultFad.val();
    dcQ_dc     = resultFad.dx(0);
  }

  // Ca - equation
  {
    Sacado::Fad::SFad<double,5> v1Var( 5, 0, v1Now );
    Sacado::Fad::SFad<double,5> v2Var( 5, 1, v2Now );
    Sacado::Fad::SFad<double,5> M_Var( 5, 2, M_Now );
    Sacado::Fad::SFad<double,5> H_Var( 5, 3, H_Now );
    Sacado::Fad::SFad<double,5> CaVar( 5, 4, CaNow );

    // parameter
    Sacado::Fad::SFad<double,5> gCaVar( model_.gCa );
    Sacado::Fad::SFad<double,5> eCaVar( model_.gCa );
    Sacado::Fad::SFad<double,5> CaGammaVar( model_.CaGamma );
    Sacado::Fad::SFad<double,5> CaTauVar( model_.CaTau );

    Sacado::Fad::SFad<double,5> resultFad = Ca_EquF( v1Var, v2Var, M_Var, H_Var, CaVar, gCaVar, eCaVar, CaGammaVar, CaTauVar );
    CaEquFvalue = resultFad.val();
    dCaF_dV1 = resultFad.dx(0);
    dCaF_dV2 = resultFad.dx(1);
    dCaF_dM  = resultFad.dx(2);
    dCaF_dH  = resultFad.dx(3);
    dCaF_dCa = resultFad.dx(4);
  }
  {
    Sacado::Fad::SFad<double,1> CaVar( 1, 0, CaNow );

    Sacado::Fad::SFad<double,1> resultFad = Ca_EquQ( CaVar );
    CaEquQvalue = resultFad.val();
    dCaQ_dCa    = resultFad.dx(0);
  }
  // calculate potassium current
  potassiumCurrent = model_.gK * pow(nNow, 4.0) * (v1Now - v2Now - model_.eK);

  // calculate sodium current
  sodiumCurrent = model_.gNa * pow(mNow, 3.0) * hNow * (v1Now - v2Now - model_.eNa);


  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::updateIntermediateVars()" << std::endl
              << "v1 = " << v1Now << std::endl
              << "v2 = " << v2Now << std::endl
              << "nNow = " << nNow << std::endl
              << "mNow = " << mNow << std::endl
              << "hNow = " << hNow << std::endl
              << "aNow = " << aNow << std::endl
              << "bNow = " <<  bNow << std::endl
              << "M_Now = " << M_Now << std::endl
              << "H_Now = " << H_Now << std::endl
              << "cNow = " << cNow << std::endl
              << "CaNow = " << CaNow << std::endl
              << "kcl1Fvalue = " << kcl1Fvalue << std::endl
              << "dkcl1F_dV1 = " << dkcl1F_dV1 << std::endl
              << "dkcl1F_dV2 = " << dkcl1F_dV2 << std::endl
              << "dkcl1F_dn = " << dkcl1F_dn << std::endl
              << "dkcl1F_dm = " << dkcl1F_dm << std::endl
              << "dkcl1F_dh = " << dkcl1F_dh << std::endl
              << "kcl2Fvalue = " << kcl2Fvalue << std::endl
              << "dkcl2F_dV1 = " << dkcl2F_dV1 << std::endl
              << "dkcl2F_dV2 = " << dkcl2F_dV2 << std::endl
              << "dkcl2F_dn = " << dkcl2F_dn << std::endl
              << "dkcl2F_dm = " << dkcl2F_dm << std::endl
              << "dkcl2F_dh = " << dkcl2F_dh << std::endl
              << "alphaN = " << alphaN<double>( v1Now ) << std::endl
              << "betaN = " << betaN<double>( v1Now ) << std::endl
              << "nEquFvalue = " << nEquFvalue << std::endl
              << "dnF_dV1 = " << dnF_dV1 << std::endl
              << "dnF_dn = " << dnF_dn << std::endl
              << "nEquQvalue = " << nEquQvalue << std::endl
              << "dnQ_dn = " << dnQ_dn << std::endl
              << "alphaM = " << alphaM<double>( v1Now )  << std::endl
              << "betaM = " << betaM<double>( v1Now ) << std::endl
              << "mEquFvalue = " << mEquFvalue << std::endl
              << "dmF_dV1 = " << dmF_dV1 << std::endl
              << "dmF_dm = " << dmF_dm << std::endl
              << "mEquQvalue = " << mEquQvalue << std::endl
              << "dmQ_dm = " << dmQ_dm << std::endl
              << "alphaH = " << alphaH<double>( v1Now ) << std::endl
              << "betaH = " << betaH<double>( v1Now ) << std::endl
              << "hEquFvalue = " << hEquFvalue << std::endl
              << "dhF_dV1 = " << dhF_dV1 << std::endl
              << "dhF_dh = " << dhF_dh << std::endl
              << "hEquQvalue = " << hEquQvalue << std::endl
              << "dhQ_dh = " << dhQ_dh << std::endl

              << "aInf = " << aInf<double>( v1Now ) << std::endl
              << "aTau = " << aTau<double>( v1Now ) << std::endl
              << "aEquFvalue = " << aEquFvalue << std::endl
              << "daF_dV1 = " << daF_dV1 << std::endl
              << "daF_da = " << daF_da << std::endl
              << "aEquQvalue = " << aEquQvalue << std::endl
              << "daQ_da = " << daQ_da << std::endl

              << "bInf = " << bInf<double>( v1Now ) << std::endl
              << "bTau = " << bTau<double>( v1Now ) << std::endl
              << "bEquFvalue = " << bEquFvalue << std::endl
              << "dbF_dV1 = " << dbF_dV1 << std::endl
              << "dbF_db = " << dbF_db << std::endl
              << "bEquQvalue = " << bEquQvalue << std::endl
              << "dbQ_db = " << dbQ_db << std::endl

              << "M_Inf = " << M_Inf<double>( v1Now ) << std::endl
              << "M_Tau = " << M_Tau<double>( v1Now ) << std::endl
              << "M_EquFvalue = " << M_EquFvalue << std::endl
              << "dMF_dV1 = " << dMF_dV1 << std::endl
              << "dMF_dM = " << dMF_dM << std::endl
              << "M_EquQvalue = " << M_EquQvalue << std::endl
              << "dMQ_dM = " << dMQ_dM << std::endl

              << "H_Inf = " << H_Inf<double>( v1Now ) << std::endl
              << "H_Tau = " << H_Tau<double>( v1Now ) << std::endl
              << "H_EquFvalue = " << H_EquFvalue << std::endl
              << "dHF_dV1 = " << dHF_dV1 << std::endl
              << "dHF_dH = " << dHF_dH << std::endl
              << "H_EquQvalue = " << H_EquQvalue << std::endl
              << "dHQ_dH = " << dHQ_dH << std::endl

              << "cEquFvalue = " << cEquFvalue << std::endl
              << "dcF_dV1 = " << dcF_dV1 << std::endl
              << "dcF_dc = " << dcF_dc << std::endl
              << "cEquQvalue = " << cEquQvalue << std::endl
              << "dcQ_dc = " << dcQ_dc << std::endl

              << "CaEquFvalue = " << CaEquFvalue << std::endl
              << "dCaF_dV1 = " << dCaF_dV1 << std::endl
              << "dCaF_dV2 = " << dCaF_dV2 << std::endl
              << "dCaF_dM = " << dCaF_dM << std::endl
              << "dCaF_dH = " << dCaF_dH << std::endl
              << "dCaF_dCa = " << dCaF_dCa << std::endl
              << "CaEquQvalue = " << CaEquQvalue << std::endl
              << "dCaQ_dCa = " << dCaQ_dCa << std::endl

              << std::endl;
  }

#endif
  return bsuccess;
}
//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState()
{
  bool bsuccess = true;

  updateIntermediateVars();

  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;
  Linear::Vector * staVectorPtr = extData.nextStaVectorPtr;

  for( int i=0; i<nSeg; i++)
  {
    (*staVectorPtr)[li_KCurrentState[i]]  = potassiumCurrent[i];
    (*staVectorPtr)[li_NaCurrentState[i]] = sodiumCurrent[i];
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  bool bsuccess = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 Neuron 4 instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  bool bsuccess = true;

  Linear::Vector * daeQVecPtr = extData.daeQVectorPtr;

  for( int i=0; i<nSeg ; i++)
  {
    (*daeQVecPtr)[li_Vol[i]]  += segQvalue[i];
    (*daeQVecPtr)[li_nPro[i]] += segNEquQvalue[i];
    (*daeQVecPtr)[li_mPro[i]] += segMEquQvalue[i];
    (*daeQVecPtr)[li_hPro[i]] += segHEquQvalue[i];
    (*daeQVecPtr)[li_aPro[i]] += segAEquQvalue[i];
    (*daeQVecPtr)[li_bPro[i]] += segBEquQvalue[i];
    (*daeQVecPtr)[li_MPro[i]] += segM_EquQvalue[i];
    (*daeQVecPtr)[li_HPro[i]] += segH_EquQvalue[i];
    (*daeQVecPtr)[li_cPro[i]] += segCEquQvalue[i];
    (*daeQVecPtr)[li_CaPro[i]] += segCaEquQvalue[i];
  }

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess=true;

  Linear::Vector * daeFVecPtr = extData.daeFVectorPtr;

  (*daeFVecPtr)[li_Pos]  += kcl1Fvalue;
  (*daeFVecPtr)[li_Neg]  += kcl2Fvalue;

  for( int i=0; i<nSeg ; i++)
  {
    (*daeFVecPtr)[li_Vol[i]]  += segFvalue[i];
    (*daeFVecPtr)[li_nPro[i]] += segNEquFvalue[i];
    (*daeFVecPtr)[li_mPro[i]] += segMEquFvalue[i];
    (*daeFVecPtr)[li_hPro[i]] += segHEquFvalue[i];
    (*daeFVecPtr)[li_aPro[i]] += segAEquFvalue[i];
    (*daeFVecPtr)[li_bPro[i]] += segBEquFvalue[i];
    (*daeFVecPtr)[li_MPro[i]] += segM_EquFvalue[i];
    (*daeFVecPtr)[li_HPro[i]] += segH_EquFvalue[i];
    (*daeFVecPtr)[li_cPro[i]] += segCEquFvalue[i];
    (*daeFVecPtr)[li_CaPro[i]] += segCaEquFvalue[i];
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 instance.
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;

  Linear::Matrix * dQdxMatPtr = extData.dQdxMatrixPtr;

  for( int i=0; i<nSeg ; i++)
  {
    (*dQdxMatPtr)[li_Vol[i]][SegVEqnVsegOffset[i]] += segQ_dV[i];
    (*dQdxMatPtr)[li_nPro[i]][NEquNNodeOffset[i]]  += dnQ_dn[i];
    (*dQdxMatPtr)[li_mPro[i]][MEquMNodeOffset[i]]  += dmQ_dm[i];
    (*dQdxMatPtr)[li_hPro[i]][HEquHNodeOffset[i]]  += dhQ_dh[i];
    (*dQdxMatPtr)[li_aPro[i]][AEquANodeOffset[i]] += daQ_da[i];
    (*dQdxMatPtr)[li_bPro[i]][BEquBNodeOffset[i]] += dbQ_db[i];
    (*dQdxMatPtr)[li_MPro[i]][M_EquM_NodeOffset[i]] += dMQ_dM[i];
    (*dQdxMatPtr)[li_HPro[i]][H_EquH_NodeOffset[i]] += dHQ_dH[i];
    (*dQdxMatPtr)[li_cPro[i]][CEquCNodeOffset[i]] += dcQ_dc[i];
    (*dQdxMatPtr)[li_CaPro[i]][CaEquCaNodeOffset[i]] += dCaQ_dCa[i];
  }

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 Neuron 4 instance.
//
// Special Notes : This is an algebraic constaint.
//
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;

  Linear::Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  (*dFdxMatPtr)[li_Pos][APosEquPosNodeOffset] += dkcl1F_dVin;
  (*dFdxMatPtr)[li_Pos][APosEquNextNodeOffset] += dkcl1F_dVs0;

  (*dFdxMatPtr)[li_Neg][ANegEquNegNodeOffset] += dkcl2F_dVout;
  (*dFdxMatPtr)[li_Neg][ANegEquLastNodeOffset] += dkcl2F_dVsn;

  for( int i=0; i<nSeg ; i++)
  {
    (*dFdxMatPtr)[li_Vol[i]][SegVEqnVpreOffset[i]] += segF_dVp[i];
    (*dFdxMatPtr)[li_Vol[i]][SegVEqnVsegOffset[i]] += segF_dV[i];
    (*dFdxMatPtr)[li_Vol[i]][SegVEqnVnexOffset[i]] += segF_dVn[i];
    (*dFdxMatPtr)[li_Vol[i]][SegVEqnNOffset[i]]    += segF_dn[i];
    (*dFdxMatPtr)[li_Vol[i]][SegVEqnMOffset[i]]    += segF_dm[i];
    (*dFdxMatPtr)[li_Vol[i]][SegVEqnHOffset[i]]    += segF_dh[i];

    (*dFdxMatPtr)[li_nPro[i]][NEquVNodeOffset[i]]  += dnF_dV[i];
    (*dFdxMatPtr)[li_nPro[i]][NEquNNodeOffset[i]]  += dnF_dn[i];
    (*dFdxMatPtr)[li_mPro[i]][MEquVNodeOffset[i]]  += dmF_dV[i];
    (*dFdxMatPtr)[li_mPro[i]][MEquMNodeOffset[i]]  += dmF_dm[i];
    (*dFdxMatPtr)[li_hPro[i]][HEquVNodeOffset[i]]  += dhF_dV[i];
    (*dFdxMatPtr)[li_hPro[i]][HEquHNodeOffset[i]]  += dhF_dh[i];

    (*dFdxMatPtr)[li_aPro[i]][AEquVNodeOffset[i]] += daF_dV[i];
    (*dFdxMatPtr)[li_aPro[i]][AEquANodeOffset[i]]   += daF_da[i];

    (*dFdxMatPtr)[li_bPro[i]][BEquVNodeOffset[i]] += dbF_dV[i];
    (*dFdxMatPtr)[li_bPro[i]][BEquBNodeOffset[i]]   += dbF_db[i];

    (*dFdxMatPtr)[li_MPro[i]][M_EquVNodeOffset[i]] += dMF_dV[i];
    (*dFdxMatPtr)[li_MPro[i]][M_EquM_NodeOffset[i]]  += dMF_dM[i];

    (*dFdxMatPtr)[li_HPro[i]][H_EquVNodeOffset[i]] += dHF_dV[i];
    (*dFdxMatPtr)[li_HPro[i]][H_EquH_NodeOffset[i]]  += dHF_dH[i];

    (*dFdxMatPtr)[li_cPro[i]][CEquVNodeOffset[i]]   += dcF_dV[i];
    (*dFdxMatPtr)[li_cPro[i]][CEquCNodeOffset[i]]     += dcF_dc[i];
    (*dFdxMatPtr)[li_cPro[i]][CEquCaNodeOffset[i]]    += dcF_dCa[i];

    (*dFdxMatPtr)[li_CaPro[i]][CaEquVNodeOffset[i]] += dCaF_dV[i];
    (*dFdxMatPtr)[li_CaPro[i]][CaEquM_NodeOffset[i]]  += dCaF_dM[i];
    (*dFdxMatPtr)[li_CaPro[i]][CaEquH_NodeOffset[i]]  += dCaF_dH[i];
    (*dFdxMatPtr)[li_CaPro[i]][CaEquCaNodeOffset[i]]  += dCaF_dCa[i];

  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::setIC ()
{
  bool bsuccess = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  //varTypeVec.resize(1);
  //varTypeVec[0] = 'I';
}


//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Model::processParams ()
{
  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
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
// Purpose       : block constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block)
{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
  //if (!given("TNOM"))
  //  tnom = getDeviceOptions().tnom;

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
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
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

// additional Declarations

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i, isize;
  isize = instanceContainer.size();

  os << std::endl;
  os << "Number of Neuron instances: " << isize << std::endl;
  os << "    name=\t\tmodelName\tParameters" << std::endl;
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


Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new DeviceMaster<Traits>( configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() ||
      ((deviceMap.find("NEURON")!=deviceMap.end()) && (levelSet.find(4)!=levelSet.end())))
  {
    Neuron::registerDevice();

    Config<Traits>::addConfiguration()
      .registerDevice("neuron", 4)
      .registerModelType("neuron", 4);
  }
}

} // namespace Neuron4
} // namespace Device
} // namespace Xyce

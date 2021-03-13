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
#include <N_DEV_MembraneCS.h>
#include <N_DEV_MembraneHH.h>
#include <N_DEV_MembranePassive.h>
#include <N_DEV_MembraneUserDefined.h>
#include <N_DEV_Neuron6.h>
#include <N_DEV_Neuron_CommonEquations.h>
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

namespace Neuron6 {

void Traits::loadInstanceParameters(ParametricData<Neuron6::Instance> &p)
{
  p.addPar ("R",1.0,&Neuron6::Instance::rInt)
   .setGivenMember(&Neuron6::Instance::rIntGiven)
   .setUnit(U_OHMM)
   .setCategory(CAT_NONE)
   .setDescription("Intracellular resistivity");
  // typical values 1-3 kOhm-mm; use 1 Kohm-mm which is the same as 1 Ohm-m
  p.addPar ("A",0.00025,&Neuron6::Instance::radius)
   .setGivenMember(&Neuron6::Instance::radiusGiven)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Segment radius");
  // 250 microns,based on NEURON default diameter of 500 microns
  p.addPar ("L",0.0001,&Neuron6::Instance::length)
   .setGivenMember(&Neuron6::Instance::lengthGiven)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Cable length");
  // 100 microns,based on NEURON default length

  p.addPar("N",1,&Neuron6::Instance::nSeg)
   .setGivenMember(&Neuron6::Instance::nSegGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Number of segments");
}

void Traits::loadModelParameters(ParametricData<Neuron6::Model> &p)
{
  p.addPar ("CMEM",0.0,&Neuron6::Model::cMem)
   .setGivenMember(&Neuron6::Model::cMemGiven)
   .setUnit(U_FARADMM2)
   .setCategory(CAT_NONE)
   .setDescription("Membrane capacitance");

  p.addPar ("GMEM",0.0,&Neuron6::Model::gMem)
   .setGivenMember(&Neuron6::Model::gMemGiven)
   .setUnit(U_OHMM1MM2)
   .setCategory(CAT_NONE)
   .setDescription("Membrane conductance");

  p.addPar ("VREST",0.0,&Neuron6::Model::vRest)
   .setGivenMember(&Neuron6::Model::vRestGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Resting potential");

  p.addPar ("EK",0.0,&Neuron6::Model::eK)
   .setGivenMember(&Neuron6::Model::eKGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Potassium resting potential");

  p.addPar ("GK",0.0,&Neuron6::Model::gK)
   .setGivenMember(&Neuron6::Model::gKGiven)
   .setUnit(U_OHMM1MM2)
   .setCategory(CAT_NONE)
   .setDescription("Potassium base conductance");

  p.addPar ("ENA",0.0,&Neuron6::Model::eNa)
   .setGivenMember(&Neuron6::Model::eNaGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Sodium resting potential");

  p.addPar ("GNA",0.0,&Neuron6::Model::gNa)
   .setGivenMember(&Neuron6::Model::gNaGiven)
   .setUnit(U_OHMM1MM2)
   .setCategory(CAT_NONE)
   .setDescription("Sodium base conductance");

  p.addPar ("EA",0.0,&Neuron6::Model::eA)
   .setGivenMember(&Neuron6::Model::eAGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("a-current rest potential");

  p.addPar ("GA",0.0,&Neuron6::Model::gA)
   .setGivenMember(&Neuron6::Model::gAGiven)
   .setUnit(U_OHMM1MM2)
   .setCategory(CAT_NONE)
   .setDescription("a-current base conductance");

  p.addPar ("ECA",0.0,&Neuron6::Model::eCa)
   .setGivenMember(&Neuron6::Model::eCaGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Calcium rest potential");

  p.addPar ("GCA",0.0,&Neuron6::Model::gCa)
   .setGivenMember(&Neuron6::Model::gCaGiven)
   .setUnit(U_OHMM1MM2)
   .setCategory(CAT_NONE)
   .setDescription("Calcium base conductance");

  p.addPar ("EKCA",0.0,&Neuron6::Model::eKCa)
   .setGivenMember(&Neuron6::Model::eKCaGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Potassium-calcium rest potential");

  p.addPar ("GKCA",0.0,&Neuron6::Model::gKCa)
   .setGivenMember(&Neuron6::Model::gKCaGiven)
   .setUnit(U_OHMM1MM2)
   .setCategory(CAT_NONE)
   .setDescription("Potassium-calcium base conductance");

  p.addPar ("CAINIT",0.0,&Neuron6::Model::CaInit)
   .setGivenMember(&Neuron6::Model::CaInitGiven)
   .setUnit(U_MOLAR)
   .setCategory(CAT_NONE)
   .setDescription("initial intra-cellular calcium concentration");

  p.addPar ("CAGAMMA",0.0,&Neuron6::Model::CaGamma)
   .setGivenMember(&Neuron6::Model::CaGammaGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("calcium current to concentration multiplier");

  p.addPar ("CATAU",0.0,&Neuron6::Model::CaTau)
   .setGivenMember(&Neuron6::Model::CaTauGiven)
   .setUnit(U_SECOND)
   .setCategory(CAT_NONE)
   .setDescription("calcium removal time constant");

  p.addPar ("R",1.0,&Neuron6::Model::rInt)
   .setGivenMember(&Neuron6::Model::rIntGiven)
   .setUnit(U_OHMM)
   .setCategory(CAT_NONE)
   .setDescription("Intracellular resistivity");
  // typical values 1-3 kOhm-mm; use 1 Kohm-mm which is the same as 1 Ohm-m

  p.addPar ("A",0.00025,&Neuron6::Model::radius)
   .setGivenMember(&Neuron6::Model::radiusGiven)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Segment radius");
  // 250 microns,based on NEURON default diameter of 500 microns

  p.addPar ("L",0.0001,&Neuron6::Model::length)
   .setGivenMember(&Neuron6::Model::lengthGiven)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Cable length");
  // 100 microns,based on NEURON default length

  p.addPar ("I",0.0,&Neuron6::Model::I)
   .setExpressionAccess(ParameterType::SOLN_DEP)
   .setUnit(U_AMP)
   .setCategory(CAT_NONE)
   .setDescription("Current for user-defined current equation");

  p.addPar("IONCHANNELMODEL","",&Neuron6::Model::ionChannelModel)
   .setGivenMember(&Neuron6::Model::ionChannelModelGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Neuron6::Model to use for ion channels");

  p.addPar("N",1,&Neuron6::Model::nSeg)
   .setGivenMember(&Neuron6::Model::nSegGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Number of segments");

  p.addPar("MM_CURRENT",std::vector<std::string>(),&Neuron6::Model::membraneCurrentEqus)
   .setGivenMember(&Neuron6::Model::membraneCurrentEqusGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Contribution to membrane current");

  p.addPar("MM_INDVARS",std::vector<std::string>(),&Neuron6::Model::membraneIndpVars)
   .setGivenMember(&Neuron6::Model::membraneIndpVarsGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Independant variables for ion channel equations");

  p.addPar("MM_INDFEQUS",std::vector<std::string>(),&Neuron6::Model::membraneIndpFEqus)
   .setGivenMember(&Neuron6::Model::membraneIndpFEqusGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Independant variables: F equations");

  p.addPar("MM_INDQEQUS",std::vector<std::string>(),&Neuron6::Model::membraneIndpQEqus)
   .setGivenMember(&Neuron6::Model::membraneIndpQEqusGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Independant variables: Q equations");

  p.addPar("MM_FUNCTIONS",std::vector<std::string>(),&Neuron6::Model::membraneFunctions)
   .setGivenMember(&Neuron6::Model::membraneFunctionsGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Functions for membrane Neuron6::Model");

  p.addPar("MM_PARAMETERS",std::vector<std::string>(),&Neuron6::Model::membraneParameters)
   .setGivenMember(&Neuron6::Model::membraneParametersGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Parameters for membrane Neuron6::Model");
}


// Class Instance

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
    gSeg(0.0),
    numIntVarsPerSegment(0),
    numStateVarsPerSegment(0),
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
  // by default we'll have a voltage at each node (nSeg)
  // and no state vars.  If the user has an ion channel on, then we'll add in vars for that
  // ask the membrane model for the number of vars it has.
  numIntVarsPerSegment = model_.membraneModel_->numVars();
  numStateVarsPerSegment = 0;

  /*
    if( model_.ConnorStevensOn_ ) {
    numIntVarsPerSegment += 9;
    numStateVarsPerSegment += 2;   // two currents per segment
    }
  */

  numIntVars = numIntVarsPerSegment*nSeg;
  numStateVars = numStateVarsPerSegment*nSeg;

  // total up number of vars.
  int numVars = numExtVars + numIntVars;

  //
  // i(n) - I(n)/A - g(n,n+1) * (V(n+1) - V(n)) - g(n,n-1) * (V(n-1) - V(n)) + Cm dV(n)/d(t) = 0 LHS needs leak term memG ( vSeg - vRest )
  //
  // Vin   : i(n) - I(n)/A - g(n,n+1) * (V(n+1) - V(n)) + Cm dV(n)/d(t) = 0
  // Vout  : i(n) - I(n)/A - g(n,n-1) * (V(n-1) - V(n)) + Cm dV(n)/d(t) = 0
  // Vnode : i(n) - I(n)/A - g(n,n+1) * (V(n+1) - V(n)) - g(n,n-1) * (V(n-1) - V(n)) + Cm dV(n)/d(t) = 0
  //        plus node supporting equations (a, b, m)
  //
  // jacobian format for just the passive cable
  //             Vin      Vout      V1        V(nSeg)
  // kcl Vin    yes               g(0,1)
  // kcl Vout             yes         g(n,n-1)
  // kcl V1     yes                yes        yes
  //

  // if the membrane model includes additional internal variables, the above changes to something
  // along these lines:
  //             Vin      Vout      V1     x1     y1    V(nSeg)     x(nSeg)     y(nSeg)
  // kcl Vin    yes               g(0,1)
  // kcl Vout             yes         g(n,n-1)
  // kcl V1     yes                yes     yes   yes     yes
  // note that V1's dependence on internal variables for segment 1 comes before dependence on
  // Vnext, except for the last segment, in which case Vnext is Vout

  // set up jacStamp.  This is dependant on the membrane model.  The only part this
  // constructor really knows about is the external variables Vin and vOut and internal node
  // voltages

  if( jacStamp.empty() )       // redundant as jacStamp is not static for this device
  {                            // it can't be as each cable may have a different number of nodes
    jacStamp.resize(numVars);
    jacStamp[0].resize(2);
    jacStamp[0][0] = 0;                               // Vin
    jacStamp[0][1] = 2;                               // Vseg[0]
    jacStamp[1].resize(2);
    jacStamp[1][0] = 1;                               // Vout
    jacStamp[1][1] = numVars - numIntVarsPerSegment;  // VnSeg[nSeg-1]

    // now loop over each segment and have the membrane model fill that instanceContainer
    for( int i=0; i<nSeg; i++ )
    {
      // the cable model should take care of the Vpre, V, Vnext dependence as that
      // is part of the cable equation.  Let the membraneModel_ handle what happens
      // at the membrane level
      int offset = numExtVars + i * numIntVarsPerSegment;	// row for vSeg
      jacStamp[offset].resize( numIntVarsPerSegment + 2 );   // + 2 for Vin and Vout

      // have to handle a few special cases which can alter the ordering of Vprev, Vseg, Vnext
      // for each segment number, save the offsets for the previous, current, and next segment so
      // we don't have to rethink these special cases every time.
      if( nSeg == 1 )
      {
        // here the ordering is Vin  Vout  Vseg
        prevMap[i] = 0;
        nextMap[i] = 1;
        segMap[i] = 2;
        jacStamp[offset][prevMap[i]] = 0;                                     // Vin
        jacStamp[offset][nextMap[i]] = 1;                                     // Vout
        jacStamp[offset][segMap[i]] = 2;                                      // Vseg

      }
      else if( i==0 )
      {
        // ordering is Vin Vseg (seg int vars)  Vnext
        prevMap[i] = 0;
        segMap[i] = 1;
        nextMap[i] = numIntVarsPerSegment+1;
        jacStamp[offset][prevMap[i]] = 0;                                     // Vin
        jacStamp[offset][segMap[i]] = offset;                                 // Vseg
        jacStamp[offset][nextMap[i]] = offset + numIntVarsPerSegment;         // Vnext
      }
      else if( i==(nSeg-1) )
      {
        // ordering is Vout Vprev Vseg
        nextMap[i] = 0;
        prevMap[i] = 1;
        segMap[i] = 2;
        jacStamp[offset][nextMap[i]] = 1;                                     // Vout
        jacStamp[offset][prevMap[i]] = offset - numIntVarsPerSegment;         // Vprev
        jacStamp[offset][segMap[i]] = offset;                                 // Vseg
      }
      else
      {
        // ordering is Vprev Vseg (seg int vars) Vnext
        prevMap[i] = 0;
        segMap[i] = 1;
        nextMap[i] = numIntVarsPerSegment+1;
        jacStamp[offset][prevMap[i]] = offset - numIntVarsPerSegment;         // Vprev
        jacStamp[offset][segMap[i]] = offset;                                 // Vseg
        jacStamp[offset][nextMap[i]] = offset + numIntVarsPerSegment;         // Vnext
      }

      // pass the membraneModel_ enough information for it to construct its part of the jacobian
      model_.membraneModel_->setJacStamp( numExtVars, i, segMap[i], jacStamp );
    }

  }

  /*
  // print out jacStamp
  Xyce::dout() << "jacStamp for Neuron6" << std::endl;
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
  // equation 6.30 Theoretical neuroscience: computational and mathematical modeling of neural systems, Peter Dayan and L. F. Abbot 2001
  //   g(n,n') = radius * (radius')^2 / ( rInt segLength * ( segLength * (radius')^2 + segLength' * (radius)^2 ) )
  // this equation is in terms of current/area, in which case conductivity to the previous and next
  // segment is not symmetric if the segment length and/or radius is not are not equal
  // But since we work in current rather than current density, this is not a problem
  double segLength = length / nSeg;
  double rLong = rInt * segLength / (M_PI * radius * radius);  // longitudinal resistance (ohm)

  // watch out for divide-by-0 cases; if resistivity is close to 0, just set conductance to some large value
  // TODO:  Really should be testing for anything close enough to zero to cause problems.  Is there a better 'large' value to use?
  if (rLong == 0.0)
  {
    gSeg = 1000.0;
  }
  else
  {
    gSeg = 1.0 / rLong;
  }

  // variable indecies loads
  if( model_.ConnorStevensOn_ )
  {
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
  else
  {
    /*
      li_Vol.resize(nSeg);
      segFvalue.resize(nSeg);
      segQvalue.resize(nSeg);

      // jacobian elements
      segF_dVp.resize(nSeg);
      segF_dV.resize(nSeg);
      segF_dVn.resize(nSeg);
      segQ_dV.resize(nSeg);
    */
  }

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
    else    //  estimate it via lambda-d rule
    {
    // equations from The Neuron Book, pgs 122-123
    //     says d_lambda value of 0.1 is adequate for most purposes;
    //     use a smaller value if membrane time constant (tau_m) is shorter than 8 ms
    //     but generally uses directly specify nSeg if the default rule doesn't apply
    double d_lambda = 0.1;
    double f = 100.0;	// frequency
    // In equations from Neuron book, C was given in uF; we use F
    double cMem = model_.cMem * 1.0e6;
    // NEURON version of lambda_f equation took d in microns, other
    // distances in cm, and returned lambda_f in microns:
    //   lambda_f = 1.0e5 * sqrt(2*radius/(4*M_PI*f*rInt*cMem));
    //   nSeg = int((length/(d_lambda*lambda_f)+0.9)/2)*2 + 1;
    // I modified the coefficient in the lambda_f equation
    // to keep distance units in cm
    // (should also work in other units as long as the units are consistent)
    double lambda_f = 1.0e3 * sqrt(2*radius/(4*M_PI*f*rInt*cMem));
    nSeg = int((length/(d_lambda*lambda_f)+0.9)/2)*2 + 1;
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

  // resize our storage location for the internal vars.
  li_internalVars.resize( numIntVars );

  // now copy in the local ids
  for( int i=0; i<numIntVars; i++ )
  {
    li_internalVars[i] = intLIDVec[i];
  }
  /*
    if( model_.ConnorStevensOn_ )
    {
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
    }
    else
    {
    for( int i=0, j=0; i<nSeg; i++, j+=1)
    {
    li_Vol[i]  = intLIDVec[j];
    }
    }

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
    {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl
    << "  li_Neg = " << li_Neg << std::endl;
    for( int i=0; i<nSeg; i++ )
    {
    Xyce::dout() << "  li_Vol[ " << i << " ] = " << li_Vol[i] << std::endl;

    if( model_.ConnorStevensOn_ )
    {
    Xyce::dout() << "  li_nPro[ " << i << " ] = " << li_nPro[i] << std::endl
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
    }
    #endif

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
    {
    Xyce::dout() << section_divider << std::endl;
    }
    #endif
  */
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
  for (int i = 0; i < nSeg ; i++)
  {
    std::ostringstream segNumber;
    segNumber << i;
    addInternalNode(symbol_table, li_internalVars[i* numIntVarsPerSegment ], getName(), "V" + segNumber.str());
    if( numIntVarsPerSegment > 1 ) 
    {
      // need to handle this better as the HH and CS (and userdefined) optinos will have differnet internal variables
      addInternalNode(symbol_table, li_internalVars[i*numIntVarsPerSegment +1], getName(), "N" + segNumber.str());
      addInternalNode(symbol_table, li_internalVars[i*numIntVarsPerSegment +2], getName(), "M" + segNumber.str());
      addInternalNode(symbol_table, li_internalVars[i*numIntVarsPerSegment +3], getName(), "H" + segNumber.str());
    }

    if( model_.ConnorStevensOn_ )
    {
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

  if( model_.ConnorStevensOn_ )
  {
    for( int i=0, j=0; i<nSeg; i++, j+=2)
    {
      li_KCurrentState[i] = staLIDVec[j];
      li_NaCurrentState[i] = staLIDVec[j+1];
    }
  }
  else
  {
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

  // resize our storage location and store the results

  int numRows = jacLIDVec.size();
  jacobianOffsets.resize( numRows );
  for( int i=0; i< numRows; i++ )
  {
    int numCol = jacLIDVec[i].size();
    jacobianOffsets[i].resize( numCol );
    for( int j=0; j< numCol; j++ )
    {
      jacobianOffsets[i][j] = jacLIDVec[i][j];
    }
  }

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

  if( model_.ConnorStevensOn_ )
  {
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
  else
  {
    // internal variables
    SegVEqnVpreOffset.resize(nSeg);
    SegVEqnVsegOffset.resize(nSeg);
    SegVEqnVnexOffset.resize(nSeg);

    for(int i=0, j=2; i<nSeg; i++, j+=numIntVarsPerSegment )
    {
      // Xyce::dout() << " i = " << i << " j = " << j << " jacLIDVec[ " << j << " ].size() = " << jacLIDVec[j].size() << std::endl;
      SegVEqnVpreOffset[i] = jacLIDVec[j][0];
      SegVEqnVsegOffset[i] = jacLIDVec[j][1];
      SegVEqnVnexOffset[i] = jacLIDVec[j][2];

      /*
        Xyce::dout() <<  SegVEqnVpreOffset[i] << ", "
        << SegVEqnVsegOffset[i] << ", "
        << SegVEqnVnexOffset[i] << std::endl;
      */
    }
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
bool Instance::updatePrimaryState ()
{
  bool bsuccess = true;

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
//                 Neuron 6 instance.
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

  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;
  Linear::Vector * daeQVecPtr = extData.daeQVectorPtr;

  // no Q component for the cable component of this devcie

  // now let the membrane model load it's Q component (memCap dV/dt, etc.)
  for( int i=0; i<nSeg ; i++)
  {
    model_.membraneModel_->loadDAEQVector ( i, li_internalVars, solVectorPtr, daeQVecPtr, segArea );
    /*
      (*daeQVecPtr)[li_Vol[i]]  += segQvalue[i];
      if( model_.ConnorStevensOn_ )
      {
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
    */
  }

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 Neuron 6 instance.
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

  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;
  Linear::Vector * daeFVecPtr = extData.daeFVectorPtr;

  double vIn = (*solVectorPtr)[li_Pos];
  double vOut = (*solVectorPtr)[li_Neg];

  // take care of the input and output nodes as they are different
  (*daeFVecPtr)[li_Pos]  += -2.0 * gSeg * ((*solVectorPtr)[li_internalVars[0]] - vIn );
  (*daeFVecPtr)[li_Neg]  += -2.0 * gSeg * ((*solVectorPtr)[li_internalVars[(nSeg-1)*numIntVarsPerSegment]] - vOut );

  for( int i=0; i<nSeg ; i++)
  {
    // for this segment get the values of the local vars
    double vSeg  = (*solVectorPtr)[li_internalVars[i*numIntVarsPerSegment]];
    double vNext = 0.0;
    double gNext = 0.0;
    if (i == (nSeg - 1))
    {
      vNext = vOut;
      gNext = gSeg * 2.0;
    }
    else
    {
      vNext = (*solVectorPtr)[li_internalVars[(i+1)*numIntVarsPerSegment]];
      gNext = gSeg;
    }
    double vPrev = 0.0;
    double gPrev = 0.0;
    if (i == 0 )
    {
      vPrev = vIn;
      gPrev = gSeg * 2.0;
    }
    else
    {
      vPrev = (*solVectorPtr)[li_internalVars[(i-1)*numIntVarsPerSegment]];
      gPrev = gSeg;
    }

    // li_internalVars lists numIntVarsPerSegment for each segment.  V for each segment will
    // be first.  So, V's offset in li_internalVars is i * numIntVarsPerSegment

    (*daeFVecPtr)[li_internalVars[i*numIntVarsPerSegment]] += - gPrev * (vPrev - vSeg) - gNext * (vNext - vSeg);

    model_.membraneModel_->loadDAEFVector ( i, li_internalVars, solVectorPtr, daeFVecPtr, segArea );

  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 Neuron 6 instance.
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;

  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;
  Linear::Matrix * dQdxMatPtr = extData.dQdxMatrixPtr;

  for( int i=0; i<nSeg ; i++)
  {
    // let the membrane model load it's part
    model_.membraneModel_->loadDAEdQdx ( i, segMap[i], li_internalVars, jacobianOffsets, solVectorPtr, dQdxMatPtr, segArea );

    /*
      if( model_.ConnorStevensOn_ )
      {
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
    */
  }

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 Neuron 6 instance.
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

  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;
  Linear::Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  (*dFdxMatPtr)[li_Pos][APosEquPosNodeOffset]  +=  2.0 * gSeg;
  (*dFdxMatPtr)[li_Pos][APosEquNextNodeOffset] += -2.0 * gSeg;

  (*dFdxMatPtr)[li_Neg][ANegEquNegNodeOffset]  +=  2.0 * gSeg;
  (*dFdxMatPtr)[li_Neg][ANegEquLastNodeOffset] += -2.0 * gSeg;

  for( int i=0; i<nSeg ; i++)
  {
    int offset = i * numIntVarsPerSegment;

    int row = numExtVars + i * numIntVarsPerSegment;

    double gPrev = gSeg;
    double gNext = gSeg;
    if (i == 0)
    {
      gPrev = gSeg * 2.0;
    }
    if (i == (nSeg-1))
    {
      gNext = gSeg * 2.0;
    }

    (*dFdxMatPtr)[li_internalVars[offset]][jacobianOffsets[row][prevMap[i]]]           +=  -gPrev;         // Vpre
    (*dFdxMatPtr)[li_internalVars[offset]][jacobianOffsets[row][segMap[i]]]     +=  gPrev + gNext;  // Vseg
    (*dFdxMatPtr)[li_internalVars[offset]][jacobianOffsets[row][nextMap[i]]] +=  -gNext;         // Vnext

    // now let the membrane model load it's part
    model_.membraneModel_->loadDAEdFdx ( i, segMap[i], li_internalVars, jacobianOffsets, solVectorPtr, dFdxMatPtr, segArea );

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
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    cMem(0.0),
    gMem(0.0),
    vRest(0.0),
    eNa(0.0),
    gNa(0.0),
    eK(0.0),
    gK(0.0),
    eA(0.0),
    gA(0.0),
    eCa(0.0),
    gCa(0.0),
    eKCa(0.0),
    gKCa(0.0),
    CaInit(0.0),
    CaGamma(0.0),
    CaTau(0.0),
    rInt(0.0),
    radius(0.0),
    length(0.0),
    nSeg(0.0),
    rIntGiven(false),
    radiusGiven(false),
    lengthGiven(false),
    ionChannelModelGiven(false),
    nSegGiven(false),
    cMemGiven(false),
    gMemGiven(false),
    vRestGiven(false),
    eNaGiven(false),
    gNaGiven(false),
    eKGiven(false),
    gKGiven(false),
    eAGiven(false),
    gAGiven(false),
    eCaGiven(false),
    gCaGiven(false),
    eKCaGiven(false),
    gKCaGiven(false),
    CaInitGiven(false),
    CaGammaGiven(false),
    CaTauGiven(false),
    membraneIndpVarsGiven(false),
    membraneIndpFEqusGiven(false),
    membraneIndpQEqusGiven(false),
    hodgenHuxleyOn_(false),
    ConnorStevensOn_(false),
    sodiumOn_(false),
    potassiumOn_(false),
    aCurrentOn_(false),
    calciumOn_(false)
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

  // check the specified model type and allocate the needed membrane model
  if( ionChannelModelGiven )
  {
    // should change the case on ionChannelModel and simplify this set of clauses
    if( ionChannelModel == "passive" || ionChannelModel == "PASSIVE" )
    {
      membraneModel_ = rcp( new MembranePassive( getSolverState(), cMem, gMem, vRest) );

    }
    else if( ionChannelModel == "hh" || ionChannelModel == "HH" )
    {
      membraneModel_ = rcp( new MembraneHH( getSolverState(), cMem, gMem, vRest, eK, gK, eNa, gNa) );
    }
    else if( ionChannelModel == "cs" || ionChannelModel == "CS")
    {
      membraneModel_ = rcp( new MembraneCS( getSolverState() ) );
    }
    else if( ionChannelModel == "ud" || ionChannelModel == "UD")
    {
      membraneModel_ = rcp( new MembraneUserDefined( getSolverState(), cMem, gMem, vRest,
                                                           membraneCurrentEqus, membraneIndpVars, membraneIndpFEqus, membraneIndpQEqus, membraneFunctions, membraneParameters) );
    }
    else
    {
      UserFatal(*this) << "Model: unknown ion channel model given \"" << ionChannelModel << "\"";
    }

  }
  else
  {
    // no model given so assume passive
    membraneModel_ = rcp( new MembranePassive( getSolverState(), cMem, gMem, vRest) );
  }
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

  return new DeviceMaster<Traits>(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() ||
      ((deviceMap.find("NEURON")!=deviceMap.end()) && (levelSet.find(6)!=levelSet.end())))
  {
    Neuron::registerDevice();

    Config<Traits>::addConfiguration()
      .registerDevice("neuron", 6)
      .registerModelType("neuron", 6);
  }
}

} // namespace Neuron6
} // namespace Device
} // namespace Xyce

//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Creator        : Richard Schiek, Electrical Systems Modeling
//
// Creation Date  : 03/08/2012
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Neuron8.h>
#include <N_DEV_Neuron_CommonEquations.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_DEV_Neuron.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {

namespace Neuron8 {

void Traits::loadInstanceParameters(ParametricData<Neuron8::Instance> &p)
{
  p.addPar ("MEMC",0.0,&Neuron8::Instance::memCap)
   .setGivenMember(&Neuron8::Instance::memCapGiven)
   .setUnit(U_FARAD)
   .setCategory(CAT_NONE)
   .setDescription("Membrane capacitance");

  p.addPar ("VT",0.0,&Neuron8::Instance::Vt)
   .setGivenMember(&Neuron8::Instance::VtGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Instantaneous threshold voltage");

  p.addPar ("VR",0.0,&Neuron8::Instance::Vr)
   .setGivenMember(&Neuron8::Instance::VrGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Resting membrane potential");

  p.addPar ("VP",0.0,&Neuron8::Instance::Vpeak)
   .setGivenMember(&Neuron8::Instance::VpeakGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Peak voltage");

  p.addPar ("K",0.0,&Neuron8::Instance::k)
   .setGivenMember(&Neuron8::Instance::kGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("modeling parameter");

  p.addPar ("A",0.0,&Neuron8::Instance::a)
   .setGivenMember(&Neuron8::Instance::aGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("modeling parameter");

  p.addPar ("B",0.0,&Neuron8::Instance::b)
   .setGivenMember(&Neuron8::Instance::bGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("modeling parameter");

  p.addPar ("C",0.0,&Neuron8::Instance::c)
   .setGivenMember(&Neuron8::Instance::cGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("modeling parameter");

  p.addPar ("D",0.0,&Neuron8::Instance::d)
   .setGivenMember(&Neuron8::Instance::dGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("modeling parameter");

  p.addPar ("USCALE",1.0e-9,&Neuron8::Instance::uscale)
   .setGivenMember(&Neuron8::Instance::uscaleGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("scaling for u variable");

  p.addPar ("FALLRATE",1.0e3,&Neuron8::Instance::fallRate)
   .setGivenMember(&Neuron8::Instance::fallRateGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("recovery rate");
}

void Traits::loadModelParameters(ParametricData<Neuron8::Model> &p)
{
  p.addPar ("MEMC",0.0,&Neuron8::Model::memCap)
   .setGivenMember(&Neuron8::Model::memCapGiven)
   .setUnit(U_FARAD)
   .setCategory(CAT_NONE)
   .setDescription("Membrane capacitance");

  p.addPar ("VT",0.0,&Neuron8::Model::Vt)
   .setGivenMember(&Neuron8::Model::VtGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Instantaneous threshold voltage");

  p.addPar ("VR",0.0,&Neuron8::Model::Vr)
   .setGivenMember(&Neuron8::Model::VrGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Resting membrane potential");

  p.addPar ("VP",0.0,&Neuron8::Model::Vpeak)
   .setGivenMember(&Neuron8::Model::VpeakGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Peak voltage");

  p.addPar ("K",0.0,&Neuron8::Model::k)
   .setGivenMember(&Neuron8::Model::kGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Neuron8::Modeling parameter");

  p.addPar ("A",0.0,&Neuron8::Model::a)
   .setGivenMember(&Neuron8::Model::aGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Neuron8::Modeling parameter");

  p.addPar ("B",0.0,&Neuron8::Model::b)
   .setGivenMember(&Neuron8::Model::bGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Neuron8::Modeling parameter");

  p.addPar ("C",0.0,&Neuron8::Model::c)
   .setGivenMember(&Neuron8::Model::cGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Neuron8::Modeling parameter");

  p.addPar ("D",0.0,&Neuron8::Model::d)
   .setGivenMember(&Neuron8::Model::dGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Neuron8::Modeling parameter");

  p.addPar ("USCALE",1.0e-9,&Neuron8::Model::uscale)
   .setGivenMember(&Neuron8::Model::uscaleGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("scaling for u variable");

  p.addPar ("FALLRATE",1.0e3,&Neuron8::Model::fallRate)
   .setGivenMember(&Neuron8::Model::fallRateGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("recovery rate");
}

//
// static class member inits
//
std::vector< std::vector<int> > Instance::jacStamp;

// Class Instance

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       Miter,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    memCap(0.0),
    Vt(0.0),
    Vr(0.0),
    Vpeak(0.0),
    k(0.0),
    a(0.0),
    b(0.0),
    c(0.0),
    d(0.0),
    uscale(1.0e-9),
    fallRate(1.0e3),
    memCapGiven(false),
    VtGiven(false),
    VrGiven(false),
    VpeakGiven(false),
    kGiven(false),
    aGiven(false),
    bGiven(false),
    cGiven(false),
    dGiven(false),
    uscaleGiven(false),
    fallRateGiven(false),
    vEquFvalue(0.0),
    vEquQvalue(0.0),
    vEqudFdv(0.0),
    vEqudFdu(0.0),
    vEqudQdv(0.0),
    uEquFvalue(0.0),
    uEquQvalue(0.0),
    uEqudFdv(0.0),
    uEqudFdu(0.0),
    uEqudQdu(0.0),
    resetting(false),
    uPeak(0.0),
    li_V(-1),
    li_U(-1),
    vEquVOffset(-1),
    vEquUOffset(-1),
    uEquVOffset(-1),
    uEquUOffset(-1)
{
  numExtVars = 1;  // membrane voltage

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

  /*
    Xyce::dout() << "Instance::Instance" << std::endl
    << "memC  = " << memCap << std::endl
    << "Vt    = " << Vt << std::endl
    << "Vr    = " << Vr << std::endl
    << "Vp    = " << Vpeak << std::endl
    << "k     = " << k << std::endl
    << "a     = " << a << std::endl
    << "b     = " << b << std::endl
    << "c     = " << c << std::endl
    << "d     = " << d << std::endl
    << "uscale = " << uscale << std::endl
    << "fallRate = " << fallRate << std::endl;
  */

  numIntVars = 1;  // u
  numStateVars = 0;

  // total up number of vars.
  int numVars = numExtVars + numIntVars;


  //
  // equations are:
  // v equation
  //  k * (v-Vr) * (v-Vt) - u + I - C dv/dt = 0
  // u euqtion
  //  a * ( b * (v-Vr) - u ) - du/dt = 0

  //  df/dx
  //            V                u
  // Vequ  2k V - Vt - Vr       -1
  // Uequ      a b              -a
  //
  // dq/dx
  //            V                u
  // Vequ    -memCap             0
  // Uequ       0               -1
  //
  // so jacStamp is dense.


  // set up jacStamp
  if( jacStamp.empty() )
  {
    jacStamp.resize(2);
    jacStamp[0].resize(2);
    jacStamp[0][0] = 0;
    jacStamp[0][1] = 1;
    jacStamp[1].resize(2);
    jacStamp[1][0] = 0;
    jacStamp[1][1] = 1;
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

}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
//-----------------------------------------------------------------------------
Instance::~Instance()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
//-----------------------------------------------------------------------------
bool Instance::processParams()
{
  // If there are any time dependent parameters, set their values at for
  // the current time.

  // pull unspecified params out of the model if they weren't specified here
  if( !memCapGiven && model_.memCapGiven )
  {
    memCap = model_.memCap;
  }
  if( !VtGiven && model_.VtGiven )
  {
    Vt = model_.Vt;
  }
  if( !VrGiven && model_.VrGiven )
  {
    Vr = model_.Vr;
  }
  if( !VpeakGiven && model_.VpeakGiven )
  {
    Vpeak = model_.Vpeak;
  }
  if( !kGiven && model_.kGiven )
  {
    k = model_.k;
  }
  if( !aGiven && model_.aGiven )
  {
    a = model_.a;
  }
  if( !bGiven && model_.bGiven )
  {
    b = model_.b;
  }
  if( !cGiven && model_.cGiven )
  {
    c = model_.c;
  }
  if( !dGiven && model_.dGiven )
  {
    d = model_.d;
  }
  if( !uscaleGiven && model_.uscaleGiven )
  {
    uscale = model_.uscale;
  }
  if( !fallRateGiven && model_.fallRateGiven )
  {
    fallRate = model_.fallRate;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
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
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
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

  li_V  = extLIDVec[0];
  li_U  = intLIDVec[0];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << "  li_V = " << li_V << std::endl
         << "  li_U = " << li_U << std::endl;

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
  addInternalNode(symbol_table, li_U, getName(), "U");
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
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
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  vEquVOffset = jacLIDVec[0][0];
  vEquUOffset = jacLIDVec[0][1];
  uEquVOffset = jacLIDVec[1][0];
  uEquUOffset = jacLIDVec[1][1];
}


//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars()
{
  bool bsuccess = true;

  // here we take the current solutions for V1, V2, n, m and h
  // and use those to calculate all the terms needed for the next
  // load cycle (F, Q, dFdX, dQdX)
  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;

  double vVal = (*solVectorPtr)[li_V];
  double uVal = (*solVectorPtr)[li_U];

  // not sure if I need to do this as it should just
  // fall out of the solution when d/dt = 0.  Possibly
  // it doesn't because there are two solutions when d/dt = 0,
  // either vVal = Vr or vVal = Vt we want vVal=Vr.
  if( getSolverState().dcopFlag  )
  {
    // vVal - Vr = 0 and uVal=0
    vEquFvalue = vVal - Vr;
    vEquQvalue = - memCap * vVal;
    vEqudFdv   = 1.0;
    vEqudFdu   = 0.0;
    vEqudQdv   = - memCap;
    uEquFvalue = uVal*uscale;
    uEquQvalue = -uVal*uscale;
    uEqudFdv   = 0.0;
    uEqudFdu   = 1.0*uscale;
    uEqudQdu   = uscale;
  }
  else
  {
    // need to cacth case where v> Vpeak
    if( !resetting && (vVal >= Vpeak) )
    {
      resetting = true;
      uPeak = uVal;
    }
    else
    {
      // only break out of the resetting phase at the start of a newton solve
      if( (getSolverState().newtonIter == 0 ) && (vVal <= (c+0.0001) ) && (uVal >= (uPeak + (d-1.0e-13)/uscale) ))
      {
        // reset is done, go back to old form
        resetting = false;
      }
    }

    if ( resetting )
    {
      //Xyce::dout() << "In resetting section. (uPeak*uscale + d) = " << (uPeak + d/uscale) << std::endl;
      vEquFvalue = -fallRate*(vVal - c) - uVal*uscale;
      vEquQvalue = -memCap* vVal;
      vEqudFdv   = -fallRate;
      vEqudFdu   = -uscale;
      vEqudQdv   = -memCap;

      uEquFvalue = -fallRate*(uVal - (uPeak + d/uscale) );
      uEquQvalue = -uVal;
      uEqudFdv   =  0.0;
      uEqudFdu   = -fallRate;
      uEqudQdu   = -1.0;
    }
    else
    {
      //Xyce::dout() << "In normal section." << std::endl;
      vEquFvalue = k * (vVal - Vr) * (vVal - Vt) - uVal*uscale;
      vEquQvalue = - memCap * vVal;
      vEqudFdv   = k * ( 2*vVal - Vt - Vr);
      vEqudFdu   = -uscale;
      vEqudQdv   = -memCap;

      uEquFvalue = a * ( b * (vVal - Vr)/uscale - uVal );
      uEquQvalue = -uVal;
      uEqudFdv   =  a * b / uscale;
      uEqudFdu   = -a;
      uEqudQdu   = -1.0;
    }


  }

  /*
    Xyce::dout() << "Instance::updateIntermediateVars()" << std::endl
    << "vEquFvalue = " <<  vEquFvalue << std::endl
    << "vEquQvalue = " << vEquQvalue << std::endl
    << "vEqudFdv   = " << vEqudFdv << std::endl
    << "vEqudFdu   = " << vEqudFdu << std::endl
    << "vEqudQdv   = " << vEqudQdv << std::endl
    << "uEquFvalue = " << uEquFvalue << std::endl
    << "uEquQvalue = " << uEquQvalue << std::endl
    << "uEqudFdv   = " << uEqudFdv << std::endl
    << "uEqudFdu   = " << uEqudFdu << std::endl
    << "uEqudQdu   = " << uEqudQdu << std::endl
    << "vVal       = " << vVal << std::endl
    << "uVal       = " << uVal << std::endl;
  */
  return bsuccess;
}

/*
  this approach can be problematic if many of these devices are out of
  phase and issueing breakpoints all the time.  Thus, not very workable -- RLS
//-----------------------------------------------------------------------------
// Function      : Instance::getInstanceBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
//-----------------------------------------------------------------------------
bool Instance::getInstanceBreakPoints(
std::vector<Util::BreakPoint> &breakPointTimes)
{
breakPointTimes.push_back(breakPoint);
return true;
}
*/

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState()
{
  bool bsuccess = true;

  updateIntermediateVars();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
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
//                 Neuron 7 instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  bool bsuccess = true;

  Linear::Vector * daeQVecPtr = extData.daeQVectorPtr;
  (*daeQVecPtr)[li_V] += vEquQvalue;
  (*daeQVecPtr)[li_U] += uEquQvalue;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 Neuron 7 instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess=true;

  Linear::Vector * daeFVecPtr = extData.daeFVectorPtr;

  (*daeFVecPtr)[li_V]  += vEquFvalue;
  (*daeFVecPtr)[li_U]  += uEquFvalue;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 Neuron 7 instance.
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;

  Linear::Matrix * dQdxMatPtr = extData.dQdxMatrixPtr;

  (*dQdxMatPtr)[li_V][vEquVOffset] += vEqudQdv;
  (*dQdxMatPtr)[li_U][uEquUOffset] += uEqudQdu;

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 Neuron 7 instance.
//
// Special Notes : This is an algebraic constaint.
//
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;

  Linear::Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  (*dFdxMatPtr)[li_V][vEquVOffset] += vEqudFdv;
  (*dFdxMatPtr)[li_V][vEquUOffset] += vEqudFdu;

  (*dFdxMatPtr)[li_U][uEquVOffset] += uEqudFdv;
  (*dFdxMatPtr)[li_U][uEquUOffset] += uEqudFdu;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
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
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
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
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    memCap(0.0),
    Vt(0.0),
    Vr(0.0),
    Vpeak(0.0),
    k(0.0),
    a(0.0),
    b(0.0),
    c(0.0),
    d(0.0),
    uscale(1.0e-9),
    fallRate(1.0e3),
    memCapGiven(false),
    VtGiven(false),
    VrGiven(false),
    VpeakGiven(false),
    kGiven(false),
    aGiven(false),
    bGiven(false),
    cGiven(false),
    dGiven(false),
    uscaleGiven(false),
    fallRateGiven(false)
{

  /*
    Xyce::dout() << "Model::Model" << std::endl
    << "memC  = " << memCap << std::endl
    << "Vt    = " << Vt << std::endl
    << "Vr    = " << Vr << std::endl
    << "Vp    = " << Vpeak << std::endl
    << "k     = " << k << std::endl
    << "a     = " << a << std::endl
    << "b     = " << b << std::endl
    << "c     = " << c << std::endl
    << "d     = " << d << std::endl
    << "uscale = " << uscale << std::endl
    << "fallRate = " << fallRate << std::endl;
  */

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
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
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
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
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
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
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
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/08/2012
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
      ((deviceMap.find("NEURON")!=deviceMap.end()) && (levelSet.find(8)!=levelSet.end())))
  {
    Neuron::registerDevice();

    Config<Traits>::addConfiguration()
      .registerDevice("neuron", 8)
      .registerModelType("neuron", 8);
  }
}

} // namespace Neuron8
} // namespace Device
} // namespace Xyce

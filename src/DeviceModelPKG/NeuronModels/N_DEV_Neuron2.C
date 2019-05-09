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
// Creator        : Richard Schiek, Electrical and Microsytem Modeling
//
// Creation Date  : 02/28/00
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
#include <N_DEV_Neuron2.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {

namespace Neuron2 {

void Traits::loadInstanceParameters(ParametricData<Neuron2::Instance> &p)
{
}

void Traits::loadModelParameters(ParametricData<Neuron2::Model> &p)
{
  p.addPar ("CMEM",0.0,&Neuron2::Model::cMem)
   .setGivenMember(&Neuron2::Model::cMemGiven)
   .setUnit(U_FARAD)
   .setCategory(CAT_NONE)
   .setDescription("Membrane capacitance");

  p.addPar ("GMEM",0.0,&Neuron2::Model::gMem)
   .setGivenMember(&Neuron2::Model::gMemGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("Membrane conductance");

  p.addPar ("VREST",0.0,&Neuron2::Model::vRest)
   .setGivenMember(&Neuron2::Model::vRestGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Resting potential");

  p.addPar ("EK",0.0,&Neuron2::Model::eK)
   .setGivenMember(&Neuron2::Model::eKGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Potassium resting potential");

  p.addPar ("GK",0.0,&Neuron2::Model::gK)
   .setGivenMember(&Neuron2::Model::gKGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("Potassium base conductance");

  p.addPar ("ENA",0.0,&Neuron2::Model::eNa)
   .setGivenMember(&Neuron2::Model::eNaGiven)
   .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Sodium resting potential");

  p.addPar ("GNA",0.0,&Neuron2::Model::gNa)
   .setGivenMember(&Neuron2::Model::gNaGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("Sodium base conductance");

  p.addPar ("EA",0.0,&Neuron2::Model::eA)
   .setGivenMember(&Neuron2::Model::eAGiven)
   .setUnit(U_CM)
   .setCategory(CAT_NONE)
   .setDescription("a-current rest potential");

  p.addPar ("GA",0.0,&Neuron2::Model::gA)
   .setGivenMember(&Neuron2::Model::gAGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("a-current base conductance");

  p.addPar ("ECA",0.0,&Neuron2::Model::eCa)
   .setGivenMember(&Neuron2::Model::eCaGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Calcium rest potential");

  p.addPar ("GCA",0.0,&Neuron2::Model::gCa)
   .setGivenMember(&Neuron2::Model::gCaGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("Calcium base conductance");

  p.addPar ("EKCA",0.0,&Neuron2::Model::eKCa)
   .setGivenMember(&Neuron2::Model::eKCaGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("Potassium-calcium rest potential");

  p.addPar ("GKCA",0.0,&Neuron2::Model::gKCa)
   .setGivenMember(&Neuron2::Model::gKCaGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("Potassium-calcium base conductance");

  p.addPar ("CAINIT",0.0,&Neuron2::Model::CaInit)
   .setGivenMember(&Neuron2::Model::CaInitGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("initial intra-cellular calcium concentration");

  p.addPar ("CAGAMMA",0.0,&Neuron2::Model::CaGamma)
   .setGivenMember(&Neuron2::Model::CaGammaGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("calcium current to concentration multiplier");

  p.addPar ("CATAU",0.0,&Neuron2::Model::CaTau)
   .setGivenMember(&Neuron2::Model::CaTauGiven)
   .setUnit(U_OHMM1)
   .setCategory(CAT_NONE)
   .setDescription("calcium removal time constant");
}

//
// static class member inits
//
std::vector< std::vector<int> > Instance::jacStamp;

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 01/02/08
//-----------------------------------------------------------------------------
bool Instance::processParams()
{
  // If there are any time dependent parameters, set their values at for
  // the current time.

  // now set the temperature related stuff.
  //updateTemperature(temp);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 01/02/08
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
// Creation Date : 01/02/08
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       Miter,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    kcl1Fvalue(0.0),
    kcl1Qvalue(0.0),
    kcl2Fvalue(0.0),
    kcl2Qvalue(0.0),
    nEquFvalue(0.0),
    nEquQvalue(0.0),
    mEquFvalue(0.0),
    mEquQvalue(0.0),
    hEquFvalue(0.0),
    hEquQvalue(0.0),
    dkcl1F_dV1(0.0),
    dkcl1F_dV2(0.0),
    dkcl1F_dn(0.0),
    dkcl1F_dm(0.0),
    dkcl1F_dh(0.0),
    dkcl1Q_dV1(0.0),
    dkcl1Q_dV2(0.0),
    dkcl2F_dV1(0.0),
    dkcl2F_dV2(0.0),
    dkcl2F_dn(0.0),
    dkcl2F_dm(0.0),
    dkcl2F_dh(0.0),
    dkcl2Q_dV1(0.0),
    dkcl2Q_dV2(0.0),
    dnF_dV1(0.0),
    dnF_dn(0.0),
    dnQ_dn(0.0),
    dmF_dV1(0.0),
    dmF_dm(0.0),
    dmQ_dm(0.0),
    dhF_dV1(0.0),
    dhF_dh(0.0),
    dhQ_dh(0.0)
{
  numExtVars   = 2;
  numIntVars   = 9;
  numStateVars = 2;

  // set up jacStamp
  if( jacStamp.empty() )
  {
    jacStamp.resize(11);
    jacStamp[0].resize(10);  // kcl 1
    jacStamp[0][0] = 0;
    jacStamp[0][1] = 1;
    jacStamp[0][2] = 2;
    jacStamp[0][3] = 3;
    jacStamp[0][4] = 4;
    jacStamp[0][5] = 5;
    jacStamp[0][6] = 6;
    jacStamp[0][7] = 7;
    jacStamp[0][8] = 8;
    jacStamp[0][9] = 9;
    jacStamp[1].resize(10);  // kcl 2
    jacStamp[1][0] = 0;
    jacStamp[1][1] = 1;
    jacStamp[1][2] = 2;
    jacStamp[1][3] = 3;
    jacStamp[1][4] = 4;
    jacStamp[1][5] = 5;
    jacStamp[1][6] = 6;
    jacStamp[1][7] = 7;
    jacStamp[1][8] = 8;
    jacStamp[1][9] = 9;
    jacStamp[2].resize(2);   // n equation
    jacStamp[2][0] = 0;
    jacStamp[2][1] = 2;
    jacStamp[3].resize(2);   // m equation
    jacStamp[3][0] = 0;
    jacStamp[3][1] = 3;
    jacStamp[4].resize(2);   // h equation
    jacStamp[4][0] = 0;
    jacStamp[4][1] = 4;
    jacStamp[5].resize(2);   // a equation
    jacStamp[5][0] = 0;
    jacStamp[5][1] = 5;
    jacStamp[6].resize(2);   // b equation
    jacStamp[6][0] = 0;
    jacStamp[6][1] = 6;
    jacStamp[7].resize(2);   // M equation
    jacStamp[7][0] = 0;
    jacStamp[7][1] = 7;
    jacStamp[8].resize(2);   // H equation
    jacStamp[8][0] = 0;
    jacStamp[8][1] = 8;
    jacStamp[9].resize(3);   // c equation
    jacStamp[9][0] = 0;
    jacStamp[9][1] = 9;
    jacStamp[9][2] = 10;
    jacStamp[10].resize(5);   // Ca equation
    jacStamp[10][0] = 0;
    jacStamp[10][1] = 1;
    jacStamp[10][2] = 7;
    jacStamp[10][3] = 8;
    jacStamp[10][4] = 10;
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  // there are no params for this device as yet, so calling setParams causes
  // the code to exit.
  // setParams (IB.params);

  // Set any non-constant parameter defaults:

  //if (!given("TEMP"))
  //  temp = getDeviceOptions().temp.dVal();

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 01/02/08
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
// Creation Date : 01/02/08
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
  li_nPro = intLIDVec[0];
  li_mPro = intLIDVec[1];
  li_hPro = intLIDVec[2];
  li_aPro = intLIDVec[3];
  li_bPro = intLIDVec[4];
  li_M_Pro = intLIDVec[5];
  li_H_Pro = intLIDVec[6];
  li_cPro = intLIDVec[7];
  li_CaPro = intLIDVec[8];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl
      << "  li_Neg = " << li_Neg << std::endl
      << "  li_nPro = " << li_nPro << std::endl
      << "  li_mPro = " << li_mPro << std::endl
      << "  li_hPro = " << li_hPro << std::endl
      << "  li_aPro = " << li_aPro << std::endl
      << "  li_bPro = " << li_bPro << std::endl
      << "  li_M_Pro = " << li_M_Pro << std::endl
      << "  li_H_Pro = " << li_H_Pro << std::endl
      << "  li_cPro = " << li_cPro << std::endl
      << "  li_CaPro = " << li_CaPro << std::endl;
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
    addInternalNode(symbol_table, li_nPro, getName(), "N");
    addInternalNode(symbol_table, li_mPro, getName(), "M");
    addInternalNode(symbol_table, li_hPro, getName(), "H");
    addInternalNode(symbol_table, li_aPro, getName(), "A");
    addInternalNode(symbol_table, li_bPro, getName(), "B");
    addInternalNode(symbol_table, li_M_Pro, getName(), "M_");
    addInternalNode(symbol_table, li_H_Pro, getName(), "H_");
    addInternalNode(symbol_table, li_cPro, getName(), "C");
    addInternalNode(symbol_table, li_CaPro, getName(), "Ca");

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 01/02/08
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // copy over the global ID lists.
  staLIDVec = staLIDVecRef;

  li_KCurrentState = staLIDVec[0];
  li_NaCurrentState = staLIDVec[1];

}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 01/02/08
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
// Creation Date : 01/02/08
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNegNodeOffset = jacLIDVec[0][1];
  APosEquNNodeOffset   = jacLIDVec[0][2];
  APosEquMNodeOffset   = jacLIDVec[0][3];
  APosEquHNodeOffset   = jacLIDVec[0][4];
  APosEquANodeOffset   = jacLIDVec[0][5];
  APosEquBNodeOffset   = jacLIDVec[0][6];
  APosEquM_NodeOffset  = jacLIDVec[0][7];
  APosEquH_NodeOffset  = jacLIDVec[0][8];
  APosEquCNodeOffset   = jacLIDVec[0][9];

  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];
  ANegEquNNodeOffset   = jacLIDVec[1][2];
  ANegEquMNodeOffset   = jacLIDVec[1][3];
  ANegEquHNodeOffset   = jacLIDVec[1][4];
  ANegEquANodeOffset   = jacLIDVec[1][5];
  ANegEquBNodeOffset   = jacLIDVec[1][6];
  ANegEquM_NodeOffset  = jacLIDVec[1][7];
  ANegEquH_NodeOffset  = jacLIDVec[1][8];
  ANegEquCNodeOffset   = jacLIDVec[1][9];

  ANEquPosNodeOffset   = jacLIDVec[2][0];
  ANEquNNodeOffset     = jacLIDVec[2][1];

  AMEquPosNodeOffset   = jacLIDVec[3][0];
  AMEquMNodeOffset     = jacLIDVec[3][1];

  AHEquPosNodeOffset   = jacLIDVec[4][0];
  AHEquHNodeOffset     = jacLIDVec[4][1];

  AAEquPosNodeOffset   = jacLIDVec[5][0];
  AAEquANodeOffset     = jacLIDVec[5][1];

  ABEquPosNodeOffset   = jacLIDVec[6][0];
  ABEquBNodeOffset     = jacLIDVec[6][1];

  AM_EquPosNodeOffset  = jacLIDVec[7][0];
  AM_EquM_NodeOffset   = jacLIDVec[7][1];

  AH_EquPosNodeOffset  = jacLIDVec[8][0];
  AH_EquH_NodeOffset   = jacLIDVec[8][1];

  ACEquPosNodeOffset   = jacLIDVec[9][0];
  ACEquCNodeOffset     = jacLIDVec[9][1];
  ACEquCaNodeOffset    = jacLIDVec[9][2];

  ACaEquPosNodeOffset  = jacLIDVec[10][0];
  ACaEquNegNodeOffset  = jacLIDVec[10][1];
  ACaEquM_NodeOffset   = jacLIDVec[10][2];
  ACaEquH_NodeOffset   = jacLIDVec[10][3];
  ACaEquCaNodeOffset   = jacLIDVec[10][4];

}


//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 01/02/08
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars()
{
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
    Sacado::Fad::SFad<double,10> eCaVar( model_.eCa );
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


  return bsuccess;
}
//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 01/02/08
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState()
{
  bool bsuccess = true;

  updateIntermediateVars();

  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;
  Linear::Vector * staVectorPtr = extData.nextStaVectorPtr;

  (*staVectorPtr)[li_KCurrentState]  = potassiumCurrent;
  (*staVectorPtr)[li_NaCurrentState] = sodiumCurrent;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 01/02/08
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
//                 Neuron2 instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 01/02/08
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  bool bsuccess = true;

  Linear::Vector * daeQVecPtr = extData.daeQVectorPtr;

  (*daeQVecPtr)[li_Pos]  += kcl1Qvalue;
  (*daeQVecPtr)[li_Neg]  += kcl2Qvalue;
  (*daeQVecPtr)[li_nPro] += nEquQvalue;
  (*daeQVecPtr)[li_mPro] += mEquQvalue;
  (*daeQVecPtr)[li_hPro] += hEquQvalue;
  (*daeQVecPtr)[li_aPro] += aEquQvalue;
  (*daeQVecPtr)[li_bPro] += bEquQvalue;
  (*daeQVecPtr)[li_M_Pro] += M_EquQvalue;
  (*daeQVecPtr)[li_H_Pro] += H_EquQvalue;
  (*daeQVecPtr)[li_cPro] += cEquQvalue;
  (*daeQVecPtr)[li_CaPro] += CaEquQvalue;

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 Neuron2 instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 01/02/08
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess=true;

  Linear::Vector * daeFVecPtr = extData.daeFVectorPtr;

  (*daeFVecPtr)[li_Pos]  += kcl1Fvalue;
  (*daeFVecPtr)[li_Neg]  += kcl2Fvalue;
  (*daeFVecPtr)[li_nPro] += nEquFvalue;
  (*daeFVecPtr)[li_mPro] += mEquFvalue;
  (*daeFVecPtr)[li_hPro] += hEquFvalue;
  (*daeFVecPtr)[li_aPro] += aEquFvalue;
  (*daeFVecPtr)[li_bPro] += bEquFvalue;
  (*daeFVecPtr)[li_M_Pro] += M_EquFvalue;
  (*daeFVecPtr)[li_H_Pro] += H_EquFvalue;
  (*daeFVecPtr)[li_cPro] += cEquFvalue;
  (*daeFVecPtr)[li_CaPro] += CaEquFvalue;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 Neuron 2 instance.
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 01/02/08
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;

  Linear::Matrix * dQdxMatPtr = extData.dQdxMatrixPtr;

  (*dQdxMatPtr)[li_Pos][APosEquPosNodeOffset] += dkcl1Q_dV1;
  (*dQdxMatPtr)[li_Pos][APosEquNegNodeOffset] += dkcl1Q_dV2;

  (*dQdxMatPtr)[li_Neg][ANegEquPosNodeOffset] += dkcl2Q_dV1;
  (*dQdxMatPtr)[li_Neg][ANegEquNegNodeOffset] += dkcl2Q_dV2;

  (*dQdxMatPtr)[li_nPro][ANEquNNodeOffset] += dnQ_dn;
  (*dQdxMatPtr)[li_mPro][AMEquMNodeOffset] += dmQ_dm;
  (*dQdxMatPtr)[li_hPro][AHEquHNodeOffset] += dhQ_dh;
  (*dQdxMatPtr)[li_aPro][AAEquANodeOffset] += daQ_da;
  (*dQdxMatPtr)[li_bPro][ABEquBNodeOffset] += dbQ_db;
  (*dQdxMatPtr)[li_M_Pro][AM_EquM_NodeOffset] += dMQ_dM;
  (*dQdxMatPtr)[li_H_Pro][AH_EquH_NodeOffset] += dHQ_dH;
  (*dQdxMatPtr)[li_cPro][ACEquCNodeOffset] += dcQ_dc;
  (*dQdxMatPtr)[li_CaPro][ACaEquCaNodeOffset] += dCaQ_dCa;

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 Neuron 2 instance.
//
// Special Notes : This is an algebraic constaint.
//
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 01/02/08
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;

  Linear::Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  (*dFdxMatPtr)[li_Pos][APosEquPosNodeOffset] += dkcl1F_dV1;
  (*dFdxMatPtr)[li_Pos][APosEquNegNodeOffset] += dkcl1F_dV2;
  (*dFdxMatPtr)[li_Pos][APosEquNNodeOffset]   += dkcl1F_dn;
  (*dFdxMatPtr)[li_Pos][APosEquMNodeOffset]   += dkcl1F_dm;
  (*dFdxMatPtr)[li_Pos][APosEquHNodeOffset]   += dkcl1F_dh;
  (*dFdxMatPtr)[li_Pos][APosEquANodeOffset]   += dkcl1F_da;
  (*dFdxMatPtr)[li_Pos][APosEquBNodeOffset]   += dkcl1F_db;
  (*dFdxMatPtr)[li_Pos][APosEquM_NodeOffset]  += dkcl1F_dM;
  (*dFdxMatPtr)[li_Pos][APosEquH_NodeOffset]  += dkcl1F_dH;
  (*dFdxMatPtr)[li_Pos][APosEquCNodeOffset]   += dkcl1F_dc;

  (*dFdxMatPtr)[li_Neg][ANegEquPosNodeOffset] += dkcl2F_dV1;
  (*dFdxMatPtr)[li_Neg][ANegEquNegNodeOffset] += dkcl2F_dV2;
  (*dFdxMatPtr)[li_Neg][ANegEquNNodeOffset]   += dkcl2F_dn;
  (*dFdxMatPtr)[li_Neg][ANegEquMNodeOffset]   += dkcl2F_dm;
  (*dFdxMatPtr)[li_Neg][ANegEquHNodeOffset]   += dkcl2F_dh;
  (*dFdxMatPtr)[li_Neg][ANegEquANodeOffset]   += dkcl2F_da;
  (*dFdxMatPtr)[li_Neg][ANegEquBNodeOffset]   += dkcl2F_db;
  (*dFdxMatPtr)[li_Neg][ANegEquM_NodeOffset]  += dkcl2F_dM;
  (*dFdxMatPtr)[li_Neg][ANegEquH_NodeOffset]  += dkcl2F_dH;
  (*dFdxMatPtr)[li_Neg][ANegEquCNodeOffset]   += dkcl2F_dc;

  (*dFdxMatPtr)[li_nPro][ANEquPosNodeOffset] += dnF_dV1;
  (*dFdxMatPtr)[li_nPro][ANEquNNodeOffset]   += dnF_dn;

  (*dFdxMatPtr)[li_mPro][AMEquPosNodeOffset] += dmF_dV1;
  (*dFdxMatPtr)[li_mPro][AMEquMNodeOffset]   += dmF_dm;

  (*dFdxMatPtr)[li_hPro][AHEquPosNodeOffset] += dhF_dV1;
  (*dFdxMatPtr)[li_hPro][AHEquHNodeOffset]   += dhF_dh;

  (*dFdxMatPtr)[li_aPro][AAEquPosNodeOffset] += daF_dV1;
  (*dFdxMatPtr)[li_aPro][AAEquANodeOffset]   += daF_da;

  (*dFdxMatPtr)[li_bPro][ABEquPosNodeOffset] += dbF_dV1;
  (*dFdxMatPtr)[li_bPro][ABEquBNodeOffset]   += dbF_db;

  (*dFdxMatPtr)[li_M_Pro][AM_EquPosNodeOffset] += dMF_dV1;
  (*dFdxMatPtr)[li_M_Pro][AM_EquM_NodeOffset]  += dMF_dM;

  (*dFdxMatPtr)[li_H_Pro][AH_EquPosNodeOffset] += dHF_dV1;
  (*dFdxMatPtr)[li_H_Pro][AH_EquH_NodeOffset]  += dHF_dH;

  (*dFdxMatPtr)[li_cPro][ACEquPosNodeOffset]   += dcF_dV1;
  (*dFdxMatPtr)[li_cPro][ACEquCNodeOffset]     += dcF_dc;
  (*dFdxMatPtr)[li_cPro][ACEquCaNodeOffset]    += dcF_dCa;

  (*dFdxMatPtr)[li_CaPro][ACaEquPosNodeOffset] += dCaF_dV1;
  (*dFdxMatPtr)[li_CaPro][ACaEquNegNodeOffset] += dCaF_dV2;
  (*dFdxMatPtr)[li_CaPro][ACaEquM_NodeOffset]  += dCaF_dM;
  (*dFdxMatPtr)[li_CaPro][ACaEquH_NodeOffset]  += dCaF_dH;
  (*dFdxMatPtr)[li_CaPro][ACaEquCaNodeOffset]  += dCaF_dCa;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 01/02/08
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
// Creation Date : 01/02/08
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
// Creation Date : 01/02/08
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
// Creation Date : 01/02/08
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
// Creation Date : 01/02/08
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
// Creation Date : 01/02/08
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
// Creation Date : 01/02/08
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
      ((deviceMap.find("NEURON")!=deviceMap.end()) && (levelSet.find(2)!=levelSet.end())))
  {
    Neuron::registerDevice();

    Config<Traits>::addConfiguration()
      .registerDevice("neuron", 2)
      .registerModelType("neuron", 2);
  }
}

} // namespace Neuron2
} // namespace Device
} // namespace Xyce

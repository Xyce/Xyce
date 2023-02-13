//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
#include <N_DEV_ExternData.h>
#include <N_DEV_MESFET.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {
namespace MESFET {

void Traits::loadInstanceParameters(ParametricData<MESFET::Instance> &p)
{
  p.addPar("TEMP", 0.0, &MESFET::Instance::temp)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setDescription("Device temperature");

  p.addPar("AREA", 1.0, &MESFET::Instance::area)
    .setUnit(U_METER2)
    .setCategory(CAT_GEOMETRY)
    .setDescription("device area");
}

void Traits::loadModelParameters(ParametricData<MESFET::Model> &p)
{
  p.addPar("AF", 1.0, &MESFET::Model::AF)
    .setUnit(U_NONE)
    .setCategory(CAT_FLICKER)
    .setDescription("Flicker noise exponent");

  p.addPar("B", 0.3, &MESFET::Model::B)
    .setUnit(U_VOLTM1)
    .setCategory(CAT_PROCESS)
    .setDescription("Doping tail parameter");

  p.addPar("BETA", 2.5e-3, &MESFET::Model::BETA)
    .setUnit(U_AMPVM2)
    .setCategory(CAT_PROCESS)
    .setDescription("Transconductance parameter");

  p.addPar("ALPHA", 2.0, &MESFET::Model::ALPHA)
    .setUnit(U_VOLTM1)
    .setCategory(CAT_PROCESS)
    .setDescription("Saturation voltage parameter");

  p.addPar("CGS", 0.0, &MESFET::Model::CGS)
    .setExpressionAccess(ParameterType::MIN_CAP)
    .setUnit(U_FARAD)
    .setCategory(CAT_CAP)
    .setDescription("Zero-bias gate-source junction capacitance");

  p.addPar("CGD", 0.0, &MESFET::Model::CGD)
    .setExpressionAccess(ParameterType::MIN_CAP)
    .setUnit(U_FARAD)
    .setCategory(CAT_CAP)
    .setDescription("Zero-bias gate-drain junction capacitance");

  p.addPar("FC", 0.5, &MESFET::Model::FC)
    .setUnit(U_FARAD)
    .setCategory(CAT_CAP)
    .setDescription("Coefficient for forward-bias depletion capacitance");

  p.addPar("IS", 1e-14, &MESFET::Model::IS)
    .setUnit(U_AMP)
    .setCategory(CAT_CURRENT)
    .setDescription("Gate junction saturation current");

  p.addPar("KF", 0.05, &MESFET::Model::KF)
    .setUnit(U_NONE)
    .setCategory(CAT_FLICKER)
    .setDescription("Flicker noise coefficient");

  p.addPar("LAMBDA",0.0, &MESFET::Model::LAMBDA)
    .setUnit(U_VOLTM1)
    .setCategory(CAT_VOLT)
    .setDescription("Channel length modulation");

  p.addPar("PB", 1.0, &MESFET::Model::PB)
    .setUnit(U_VOLT)
    .setCategory(CAT_VOLT)
    .setDescription("Gate junction potential");

  p.addPar("RD", 0.0, &MESFET::Model::RD)
    .setExpressionAccess(ParameterType::MIN_RES)
    .setUnit(U_OHM)
    .setCategory(CAT_RES)
    .setDescription("Drain ohmic resistance");

  p.addPar("RS", 0.0, &MESFET::Model::RS)
    .setExpressionAccess(ParameterType::MIN_RES)
    .setUnit(U_OHM)
    .setCategory(CAT_RES)
    .setDescription("Source ohmic resistance");

  p.addPar("TNOM", 0.0, &MESFET::Model::TNOM)
   .setUnit(STANDARD)
   .setCategory(CAT_NONE)
   .setDescription("Parameter measurement temperature");

  p.addPar("VTO", 0.0, &MESFET::Model::VTO)
    .setUnit(U_VOLT)
    .setCategory(CAT_VOLT)
    .setDescription("Threshold voltage");

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

//------------------- Class Model ---------------------------------
//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : pmc
// Creation Date : 11/16/2003
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
        AF(1.0),
        B(0.3),
        ALPHA(2.0),
        BETA(2.5e-3),
        CGS(0.0),
        CGD(0.0),
        FC(0.5),
        IS(1.0e-14),
        KF(0.0),
        LAMBDA(0.0),
        PB(1.0),
        RD(0.0),
        RS(0.0),
        TNOM(CONSTREFTEMP),
        VTO(-2.0),
        fNcoef(0.0),
        fNexp(1.0),
        dtype(CONSTNMOS)
{
  if (getType() != "")
  {
    if (getType() == "NMF") {
      dtype = CONSTNMOS;
    }
    else if (getType() == "PMF") {
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
    TNOM = getDeviceOptions().tnom;

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
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i,isize;
  isize = instanceContainer.size();
  os << std::endl;
  os << "Number of MESFET Instances: " << isize << std::endl;
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
// Creation Date : 11/16/2003
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

//------------------------ Class Instance -------------------------
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : pmc
// Creation Date : 11/16/2003
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & instance_block,
  Model &               model,
  const FactoryBlock &  factory_block)
  : DeviceInstance(instance_block, configuration.getInstanceParameters(), factory_block),
    model_(model),
    limitedFlag(false),
    off(0),
    ic(0),
    area(1.0),
    ic_vds(0.0),
    ic_vgs(0.0),
    temp(getDeviceOptions().temp.getImmutableValue<double>()),
    drainCond(0.0),
    sourceCond(0.0),
    tCGS(0.0),
    tCGD(0.0),
    tIS(0.0),
    tPB(0.0),
    tMESb(0.0),
    tBeta(0.0),
    tvt0(0.0),
    tLambda(0.0),
    tAlpha(0.0),
    tRD(0.0),
    tRS(0.0),
    dNode(0),
    gNode(0),
    sNode(0),
    dpNode(0),
    spNode(0),
    Vgs(0.0),
    Vgd(0.0),
    gm(0.0),
    gds(0.0),
    ggs(0.0),
    ggd(0.0),
    Bfac(0.0),
    p(0.0),
    // Solution variables and intermediate quantities
    // drain,source,gate, drainprime and sourceprime voltages
    Vd(0.0),
    Vs(0.0),
    Vg(0.0),
    Vdp(0.0),
    Vsp(0.0),
    // vector local indices
    li_Drain(-1),
    li_DrainPrime(-1),
    li_Source(-1),
    li_SourcePrime(-1),
    li_Gate(-1),
  // Jacobian Matrix
   // Jacobian Matrix Offset:
  // V_d Row:
    ADrainEquDrainNodeOffset(-1),
    ADrainEquDrainPrimeNodeOffset(-1),
  // V_g Row:
    AGateEquGateNodeOffset(-1),
    AGateEquDrainPrimeNodeOffset(-1),
    AGateEquSourcePrimeNodeOffset(-1),
  // V_s Row:
    ASourceEquSourceNodeOffset(-1),
    ASourceEquSourcePrimeNodeOffset(-1),
  // V_d' Row:
    ADrainPrimeEquDrainNodeOffset(-1),
    ADrainPrimeEquGateNodeOffset(-1),
    ADrainPrimeEquDrainPrimeNodeOffset(-1),
    ADrainPrimeEquSourcePrimeNodeOffset(-1),
 // V_s' Row:
    ASourcePrimeEquGateNodeOffset(-1),
    ASourcePrimeEquSourceNodeOffset(-1),
    ASourcePrimeEquDrainPrimeNodeOffset(-1),
    ASourcePrimeEquSourcePrimeNodeOffset(-1),

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
   // dFdx Matrix Ptr:
  // V_d Row:
    f_DrainEquDrainNodePtr(0),
    f_DrainEquDrainPrimeNodePtr(0),
  // V_g Row:
    f_GateEquGateNodePtr(0),
    f_GateEquDrainPrimeNodePtr(0),
    f_GateEquSourcePrimeNodePtr(0),
  // V_s Row:
    f_SourceEquSourceNodePtr(0),
    f_SourceEquSourcePrimeNodePtr(0),
  // V_d' Row:
    f_DrainPrimeEquDrainNodePtr(0),
    f_DrainPrimeEquGateNodePtr(0),
    f_DrainPrimeEquDrainPrimeNodePtr(0),
    f_DrainPrimeEquSourcePrimeNodePtr(0),
 // V_s' Row:
    f_SourcePrimeEquGateNodePtr(0),
    f_SourcePrimeEquSourceNodePtr(0),
    f_SourcePrimeEquDrainPrimeNodePtr(0),
    f_SourcePrimeEquSourcePrimeNodePtr(0),

   // dQdx Matrix Ptr:
  // V_d Row:
    q_DrainEquDrainNodePtr(0),
    q_DrainEquDrainPrimeNodePtr(0),
  // V_g Row:
    q_GateEquGateNodePtr(0),
    q_GateEquDrainPrimeNodePtr(0),
    q_GateEquSourcePrimeNodePtr(0),
  // V_s Row:
    q_SourceEquSourceNodePtr(0),
    q_SourceEquSourcePrimeNodePtr(0),
  // V_d' Row:
    q_DrainPrimeEquDrainNodePtr(0),
    q_DrainPrimeEquGateNodePtr(0),
    q_DrainPrimeEquDrainPrimeNodePtr(0),
    q_DrainPrimeEquSourcePrimeNodePtr(0),
 // V_s' Row:
    q_SourcePrimeEquGateNodePtr(0),
    q_SourcePrimeEquSourceNodePtr(0),
    q_SourcePrimeEquDrainPrimeNodePtr(0),
    q_SourcePrimeEquSourcePrimeNodePtr(0),
#endif
    vgs(0.0),
    vgd(0.0),
    vgs_orig(0.0),
    vgd_orig(0.0),
    vgs_old(0.0),
    vgd_old(0.0),
    vds_old(0.0),

    mode(1),
    
    capgs(0.0),
    qgs(0.0),
    cqgs(0.0),
    capgd(0.0),
    qgd(0.0),
    cqgd(0.0),
    // local indices
    li_store_vgs(-1),
    li_store_vgd(-1),
    li_branch_dev_id(-1),
    li_branch_dev_ig(-1),
    li_branch_dev_is(-1),
    li_state_qgs(-1),
    li_state_gcgs(-1),
    li_state_qgd(-1),
    li_state_gcgd(-1)
{
  numIntVars   = 2;
  numExtVars   = 3;
  numStateVars = 4;
  setNumStoreVars(2);
  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 3;    // this is the space to allocate if lead current or power is needed.

  devConMap.resize(3);
  devConMap[0] = 1;
  devConMap[1] = 2;
  devConMap[2] = 1;

  if( jacStamp.empty() )
  {
    // stamp for RS!=0, RD!=0
    jacStamp_DC_SC.resize(5);
    jacStamp_DC_SC[0].resize(2);  // Drain row
    jacStamp_DC_SC[0][0]=0;       // d-d
    jacStamp_DC_SC[0][1]=3;       // d-d'
    jacStamp_DC_SC[1].resize(3);  // Gate row
    jacStamp_DC_SC[1][0]=1;       // g-g
    jacStamp_DC_SC[1][1]=3;       // g-d'
    jacStamp_DC_SC[1][2]=4;       // g-s'
    jacStamp_DC_SC[2].resize(2);  // Source row
    jacStamp_DC_SC[2][0]=2;       // s-s
    jacStamp_DC_SC[2][1]=4;       // s-s'
    jacStamp_DC_SC[3].resize(4);  // Drain' row
    jacStamp_DC_SC[3][0]=0;       // d'-d
    jacStamp_DC_SC[3][1]=1;       // d'-g
    jacStamp_DC_SC[3][2]=3;       // d'-d'
    jacStamp_DC_SC[3][3]=4;       // d'-s'
    jacStamp_DC_SC[4].resize(4);  // Source' row
    jacStamp_DC_SC[4][0]=1;       // s'-g
    jacStamp_DC_SC[4][1]=2;       // s'-s
    jacStamp_DC_SC[4][2]=3;       // s'-d'
    jacStamp_DC_SC[4][3]=4;       // s'-s'

    jacMap_DC_SC.clear();
    jacStampMap(jacStamp_DC_SC, jacMap_DC_SC, jacMap2_DC_SC,
                jacStamp_DC,    jacMap_DC, jacMap2_DC, 4, 2, 5);

    jacStampMap(jacStamp_DC_SC, jacMap_DC_SC, jacMap2_DC_SC,
                jacStamp_SC,    jacMap_SC, jacMap2_SC, 3, 0, 5);

    jacStampMap(jacStamp_DC, jacMap_DC, jacMap2_DC,
                jacStamp,    jacMap, jacMap2, 3, 0, 5);

  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams(instance_block.params);

  // Set any non-constant parameter defaults:
  if (!given("TEMP"))
    temp = getDeviceOptions().temp.getImmutableValue<double>();

  updateDependentParameters();

  // Calculate any parameters specified as expressions:
  processParams ();

  numIntVars = (((sourceCond == 0.0)?0:1)+((drainCond == 0.0)?0:1));

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

//----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : pmc
// Creation Date : 11/16/2003
//----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                                       const std::vector<int> & extLIDVecRef )
{
  numIntVars = (((sourceCond == 0.0)?0:1)+((drainCond == 0.0)?0:1));

  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  Instance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
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
  li_Gate  = extLIDVec[1];
  li_Source = extLIDVec[2];

  int intLoc = 0;

  if( drainCond )
    li_DrainPrime = intLIDVec[intLoc++];
  else
    li_DrainPrime = li_Drain;

  if( sourceCond )
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
  if (drainCond != 0.0)
    addInternalNode(symbol_table, li_DrainPrime, getName(), "drainprime");

  if (sourceCond != 0.0)
    addInternalNode(symbol_table, li_SourcePrime, getName(), "sourceprime");

  if (loadLeadCurrent)
  {
    addBranchDataNode(symbol_table, li_branch_dev_id, getName(), "BRANCH_DD");
    addBranchDataNode(symbol_table, li_branch_dev_is, getName(), "BRANCH_DS");
    addBranchDataNode(symbol_table, li_branch_dev_ig, getName(), "BRANCH_DG");
  }
}

//----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//----------------------------------------------------------------------------
void Instance::registerStateLIDs(const std::vector<int> & staLIDVecRef)
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
  li_state_qgs  = staLIDVec[lid++];
  li_state_gcgs = staLIDVec[lid++];

  li_state_qgd  = staLIDVec[lid++];
  li_state_gcgd = staLIDVec[lid++];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "  State local indices:" << std::endl;
    Xyce::dout() << std::endl;

    Xyce::dout() << "  li_state_qgs       = " << li_state_qgs << std::endl;
    Xyce::dout() << "  li_state_gcgs      = " << li_state_gcgs;
    Xyce::dout() << "  li_state_qgd       = " << li_state_qgd;
    Xyce::dout() << "  li_state_gcgd      = " << li_state_gcgd << std::endl;;

    Xyce::dout() << section_divider << std::endl;
  }

}

//----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/9/11
//----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef)
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());

  // Copy over the global ID lists:
  stoLIDVec = stoLIDVecRef;

  int lid=0;
  li_store_vgs = stoLIDVec[lid++];
  li_store_vgd = stoLIDVec[lid++];
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
    li_branch_dev_id =  branchLIDVecRef[0];
    li_branch_dev_ig =  branchLIDVecRef[1];
    li_branch_dev_is =  branchLIDVecRef[2];
  }
}


//----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  if( drainCond != 0.0 && sourceCond != 0.0 )
    return jacStamp_DC_SC;
  else if( drainCond != 0.0 && sourceCond == 0.0 )
    return jacStamp_DC;
  else if( drainCond == 0.0 && sourceCond != 0.0 )
    return jacStamp_SC;
  else if( drainCond == 0.0 && sourceCond == 0.0 )
    return jacStamp;
  else
    return jacStamp;
}

//----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  std::vector<int> map;
  std::vector< std::vector<int> > map2;

  if (drainCond != 0.0)
  {
    if (sourceCond != 0.0)
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
    if (sourceCond != 0.0)
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
  AGateEquDrainPrimeNodeOffset         = jacLIDVec[map[1]][map2[1][1]];
  AGateEquSourcePrimeNodeOffset        = jacLIDVec[map[1]][map2[1][2]];

  ASourceEquSourceNodeOffset           = jacLIDVec[map[2]][map2[2][0]];
  ASourceEquSourcePrimeNodeOffset      = jacLIDVec[map[2]][map2[2][1]];

  ADrainPrimeEquDrainNodeOffset        = jacLIDVec[map[3]][map2[3][0]];
  ADrainPrimeEquGateNodeOffset         = jacLIDVec[map[3]][map2[3][1]];
  ADrainPrimeEquDrainPrimeNodeOffset   = jacLIDVec[map[3]][map2[3][2]];
  ADrainPrimeEquSourcePrimeNodeOffset  = jacLIDVec[map[3]][map2[3][3]];

  ASourcePrimeEquGateNodeOffset        = jacLIDVec[map[4]][map2[4][0]];
  ASourcePrimeEquSourceNodeOffset      = jacLIDVec[map[4]][map2[4][1]];
  ASourcePrimeEquDrainPrimeNodeOffset  = jacLIDVec[map[4]][map2[4][2]];
  ASourcePrimeEquSourcePrimeNodeOffset = jacLIDVec[map[4]][map2[4][3]];
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

  // F-matrix:
  f_DrainEquDrainNodePtr             = 	&(dFdx[li_Drain][ADrainEquDrainNodeOffset]);
  f_DrainEquDrainPrimeNodePtr        = 	&(dFdx[li_Drain][ADrainEquDrainPrimeNodeOffset]);

  f_GateEquGateNodePtr               = 	&(dFdx[li_Gate][AGateEquGateNodeOffset]);
  f_GateEquDrainPrimeNodePtr         = 	&(dFdx[li_Gate][AGateEquDrainPrimeNodeOffset]);
  f_GateEquSourcePrimeNodePtr        = 	&(dFdx[li_Gate][AGateEquSourcePrimeNodeOffset]);

  f_SourceEquSourceNodePtr           = 	&(dFdx[li_Source][ASourceEquSourceNodeOffset]);
  f_SourceEquSourcePrimeNodePtr      = 	&(dFdx[li_Source][ASourceEquSourcePrimeNodeOffset]);

  f_DrainPrimeEquDrainNodePtr        = 	&(dFdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset]);
  f_DrainPrimeEquGateNodePtr         = 	&(dFdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset]);
  f_DrainPrimeEquDrainPrimeNodePtr   = 	&(dFdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]);
  f_DrainPrimeEquSourcePrimeNodePtr  = 	&(dFdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]);

  f_SourcePrimeEquGateNodePtr        = 	&(dFdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset]);
  f_SourcePrimeEquSourceNodePtr      = 	&(dFdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset]);
  f_SourcePrimeEquDrainPrimeNodePtr  = 	&(dFdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]);
  f_SourcePrimeEquSourcePrimeNodePtr = 	&(dFdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]);

  // Q-matrix:
  q_DrainEquDrainNodePtr             = 	&(dQdx[li_Drain][ADrainEquDrainNodeOffset]);
  q_DrainEquDrainPrimeNodePtr        = 	&(dQdx[li_Drain][ADrainEquDrainPrimeNodeOffset]);

  q_GateEquGateNodePtr               = 	&(dQdx[li_Gate][AGateEquGateNodeOffset]);
  q_GateEquDrainPrimeNodePtr         = 	&(dQdx[li_Gate][AGateEquDrainPrimeNodeOffset]);
  q_GateEquSourcePrimeNodePtr        = 	&(dQdx[li_Gate][AGateEquSourcePrimeNodeOffset]);

  q_SourceEquSourceNodePtr           = 	&(dQdx[li_Source][ASourceEquSourceNodeOffset]);
  q_SourceEquSourcePrimeNodePtr      = 	&(dQdx[li_Source][ASourceEquSourcePrimeNodeOffset]);

  q_DrainPrimeEquDrainNodePtr        = 	&(dQdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset]);
  q_DrainPrimeEquGateNodePtr         = 	&(dQdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset]);
  q_DrainPrimeEquDrainPrimeNodePtr   = 	&(dQdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]);
  q_DrainPrimeEquSourcePrimeNodePtr  = 	&(dQdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]);

  q_SourcePrimeEquGateNodePtr        = 	&(dQdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset]);
  q_SourcePrimeEquSourceNodePtr      = 	&(dQdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset]);
  q_SourcePrimeEquDrainPrimeNodePtr  = 	&(dQdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]);
  q_SourcePrimeEquSourcePrimeNodePtr = 	&(dQdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]);

#endif
}

//----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : pmc
// Creation Date : 11/16/2003
//----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  double * staVec = extData.nextStaVectorRawPtr;
  bool bsuccess = updateIntermediateVars ();

  double * stoVec = extData.nextStoVectorRawPtr;
  stoVec[li_store_vgs] = vgs;
  stoVec[li_store_vgd] = vgd;
  staVec[li_state_qgs] = qgs;
  staVec[li_state_qgd] = qgd;

  return  bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : pmc
// Creation Date : 11/16/2003
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  double * solVec = extData.nextSolVectorRawPtr;
  double * currStaVec = extData.currStaVectorRawPtr;

  int    dtype;
  double csat, betap;
  double vgst, vgdt;
  double evgs, evgd;
  double sarg, vtf;;

// from the spice jfet
  double czgd, czgs;
  double czgdf2, czgsf2;
  double fcpb2;
  double twop;
  int    icheck, ichk1;

// for the Shockley version
//  double A, B, C, B12, C12, D, Vdsat;
//  double delta;
  double prod, denom, invdenom, afact, lfact;

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

  icheck = 1;
  dtype  = model_.dtype;

  //  we need our solution variables for any of this stuff
  Vd  = 0.0;
  Vs  = 0.0;
  Vg  = 0.0;
  Vdp = 0.0;
  Vsp = 0.0;

  Vd  = solVec[li_Drain];
  Vg  = solVec[li_Gate];
  Vs  = solVec[li_Source];
  Vsp = solVec[li_SourcePrime];
  Vdp = solVec[li_DrainPrime];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << " " << std::endl;
    Xyce::dout() << " Vg  = " << Vg << std::endl;
    Xyce::dout() << " Vd  = " << Vd << std::endl;
    Xyce::dout() << " Vs  = " << Vs << std::endl;
    Xyce::dout() << " Vdp = " << Vdp << std::endl;
    Xyce::dout() << " Vsp = " << Vsp << std::endl;
  }

  // now we need voltage drops
  Vddp  = Vd - Vdp;
  Vssp  = Vs - Vsp;
  Vgsp  = Vg - Vsp;
  Vgdp  = Vg - Vdp;
  Vdpsp = Vdp - Vsp;

  // Now the things that the 3f5 code really uses
  vgs = dtype * Vgsp;
  vgd = dtype * Vgdp;
  vds = vgs-vgd;

  origFlag = 1;
  limitedFlag = false;
  vgs_orig = vgs;
  vgd_orig = vgd;
  vds_orig = vds;

  if (getSolverState().newtonIter == 0)
  {
    if (getSolverState().initJctFlag_ && getDeviceOptions().voltageLimiterFlag)
    {
      if (getSolverState().inputOPFlag)
      {
        Linear::Vector * flagSolVectorPtr = extData.flagSolVectorPtr;
        if ((*flagSolVectorPtr)[li_Drain] == 0 || (*flagSolVectorPtr)[li_Gate] == 0 ||
            (*flagSolVectorPtr)[li_Source] == 0 || (*flagSolVectorPtr)[li_SourcePrime] ||
            (*flagSolVectorPtr)[li_DrainPrime] )
        {
          vgs = 0;
          vgd = 0;
          vds = vgs-vgd;
        }
      }
      else
      {
        vgs = 0;
        vgd = 0;
        vds = vgs-vgd;
      }
    }
    if (!(getSolverState().dcopFlag)||(getSolverState().locaEnabledFlag && getSolverState().dcopFlag))
    {
      double * currStoVec = extData.currStoVectorRawPtr;
      vgs_old = currStoVec[li_store_vgs];
      vgd_old = currStoVec[li_store_vgd];
    }
    else
    { // there is no history
      vgs_old = vgs;
      vgd_old = vgd;
    }
  }
  else
  {
    double *stoVec = extData.nextStoVectorRawPtr;
    vgs_old = stoVec[li_store_vgs];
    vgd_old = stoVec[li_store_vgd];
  }

  // SPICE-type Voltage Limiting
  ///////////////////////////////

  if (getDeviceOptions().voltageLimiterFlag)
  {
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << " before limiting: " << std::endl;
      Xyce::dout() << " vgs = " << vgs <<  "   vgs_old = " << vgs_old << std::endl;
      Xyce::dout() << " vgd = " << vgd <<  "   vgd_old = " << vgd_old << std::endl;
    }

    ichk1=1;
    vgs = devSupport.pnjlim(vgs, vgs_old, vt, vcrit, &icheck);
    vgd = devSupport.pnjlim(vgd, vgd_old, vt, vcrit, &ichk1);

    if (ichk1 == 1) {icheck=1;}
    if (icheck == 1) limitedFlag=true;

    vgs = devSupport.fetlim(vgs, vgs_old, tvt0);
    vgd = devSupport.fetlim(vgd, vgd_old, tvt0);
    vds = vgs-vgd;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << " After limiting: " << std::endl;
      Xyce::dout() << " vgs = " << vgs << std::endl;
      Xyce::dout() << " vgd = " << vgd << std::endl;
      Xyce::dout() << " " << std::endl;
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "vgs   = " << vgs << std::endl;
    Xyce::dout() << "vgd   = " << vgd << std::endl;
    Xyce::dout() << "vds   = " << vds << std::endl;
    Xyce::dout() << "Vddp  = " << Vddp << std::endl;
    Xyce::dout() << "Vssp  = " << Vssp << std::endl;
    Xyce::dout() << "Vgsp  = " << Vgsp << std::endl;
    Xyce::dout() << "Vgdp  = " << Vgdp << std::endl;
    Xyce::dout() << "Vdpsp = " << Vdpsp << std::endl;
    Xyce::dout() << " " << std::endl;
  }
  // Now set the origFlag
  if (vgs_orig != vgs || vds_orig != vds || vgd_orig != vgd) origFlag = 0;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    if (origFlag == 0)
    {
      Xyce::dout() << " Something modified the voltages. " << std::endl;
      Xyce::dout() << " Voltage       before                    after                   diff " << std::endl;
      Xyce::dout() << " vgs  " << vgs_orig << "  " << vgs << " " << vgs-vgs_orig << std::endl;
      Xyce::dout() << " vgd  " << vgd_orig << "  " << vgd << " " << vgd-vgd_orig << std::endl;
      Xyce::dout() << " vds  " << vds_orig << "  " << vds << " " << vds-vds_orig << std::endl;
      Xyce::dout() << " " << std::endl;
    }
  }

  //
  //  the following block of code evaluates the dc current and its
  //  derivatives and the charges associated with the gate and
  //  channel
  //

  // vt set in updateTemperature
  vtf = 5.0*vt;
  csat = tIS;
  if (vgs <= -vtf)
  {
    ggs = -csat/vgs + getDeviceOptions().gmin;
    cg  = ggs*vgs;
  }
  else
  {
    evgs = exp(vgs/vt);
    ggs = csat*evgs/vt + getDeviceOptions().gmin;
    cg = csat*(evgs-1) + getDeviceOptions().gmin*vgs;
  }
  if (vgd <= -vtf)
  {
    ggd = -csat/vgd + getDeviceOptions().gmin;
    cgd = ggd*vgd;
  }
  else
  {
    evgd = exp(vgd/vt);
    ggd = csat*evgd/vt  + getDeviceOptions().gmin;
    cgd = csat*(evgd-1) + getDeviceOptions().gmin*vgd;
  }
  cg = cg + cgd;

  // 3f5 does this simple stuff
  if (vds >= 0)
    mode = 1;
  else
    mode = -1;

  if (vds >= 0)  // normal mode
  {
    vgst = vgs-tvt0;
    if (vgst <= 0)
    {
      //
      //   normal mode, cutoff region
      //
      cdrain = 0;
      gm = 0;
      gds = 0;
    }
    else
    {
      prod = 1 + tLambda*vds;
      betap = tBeta*prod;
      denom = 1 + tMESb*vgst;
      invdenom = 1/denom;
      if (vds >= ( 3/tAlpha ) )
      {
        //
        //   normal mode, saturation region
        //
        cdrain = betap*vgst*vgst*invdenom;
        gm = betap*vgst*(1 + denom)*invdenom*invdenom;
        gds = tLambda*tBeta*vgst*vgst*invdenom;
      }
      else
      {
        //
        //   normal mode, linear region
        //
        afact = 1 - tAlpha*vds/3;
        lfact = 1 - afact*afact*afact;
        cdrain = betap*vgst*vgst*invdenom*lfact;
        gm = betap*vgst*(1 + denom)*invdenom*invdenom*lfact;
        gds = tBeta*vgst*vgst*invdenom*(tAlpha*afact*afact*prod + lfact*tLambda);
      }
    }
  }
  else   // inverse mode
  {
    vgdt = vgd - tvt0;
    if (vgdt <= 0)
    {
      //
      //   inverse mode, cutoff region
      //
      cdrain = 0;
      gm = 0;
      gds = 0;
    }
    else
    {
      //
      //   inverse mode, saturation region
      //
      prod = 1 - tLambda*vds;
      betap = tBeta*prod;
      denom = 1 + tMESb*vgdt;
      invdenom = 1/denom;
      if ( -vds >= 3/tAlpha )
      {
        cdrain = -betap*vgdt*vgdt*invdenom;
        gm = -betap*vgdt*(1 + denom)*invdenom*invdenom;
        gds = tLambda*tBeta*vgdt*vgdt*invdenom - gm;
      }
      else
      {
        //
        //  inverse mode, linear region
        //
        afact = 1 + tAlpha*vds/3;
        lfact = 1 - afact*afact*afact;
        cdrain = -betap*vgdt*vgdt*invdenom*lfact;
        gm = -betap*vgdt*(1 + denom)*invdenom*invdenom*lfact;
        gds = tBeta*vgdt*vgdt*invdenom*(tAlpha*afact*afact*prod
                                                    + lfact*tLambda) - gm;
      }
    }
  }
  cd = cdrain-cgd;

  //
  //     charge storage elements
  //
  twop  = 2.0*tPB;
  fcpb2 = corDepCap*corDepCap;
  czgs  = tCGS;
  czgd  = tCGD;
  if(czgs != 0)
  {
    czgsf2=czgs/f2;
    if (vgs < corDepCap)
    {
      sarg=sqrt(1-vgs/tPB);
      qgs = twop*czgs*(1-sarg);
      capgs=czgs/sarg;
    }
    else
    {
      qgs = czgs*f1 + czgsf2*(f3 *(vgs - corDepCap)
           +(vgs*vgs - fcpb2)/(2*twop));
      capgs=czgsf2*(f3 + vgs/twop);
    }
  }
  else
  {
    qgs=0.0;
    capgs=0.0;
  }

  if(czgd != 0)
  {
    czgdf2=czgd/f2;
    if (vgd < corDepCap)
    {
      sarg=sqrt(1-vgd/tPB);
      qgd = twop*czgd*(1-sarg);
      capgd=czgd/sarg;
    }
    else
    {
      qgd = czgd*f1 + czgdf2*( f3*(vgd - corDepCap)
          +(vgd*vgd - fcpb2)/(2*twop) );
      capgd=czgdf2*(f3 + vgd/twop);
    }
  }
  else
  {
    qgd=0.0;
    capgd=0.0;
  }

  Idrain  = drainCond  * Vddp;
  Isource = sourceCond * Vssp;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << " Done with Instance::updateIntermediateVars." << std::endl;
    Xyce::dout() << "  mode    = " << mode << std::endl;
    Xyce::dout() << "  tBeta   = " << tBeta << std::endl;
    Xyce::dout() << "  Idrain  = " << Idrain << std::endl;
    Xyce::dout() << "  Isource = " << Isource << std::endl;
    Xyce::dout() << "  gds     = " << gds << std::endl;
    Xyce::dout() << "  gm      = " << gm << std::endl;
  }

  /// CURRENTS to load into RHS:

  // so at this point:

  // current out of drain is
  // Idrain

  // current out of gate:
  // dtype*( d/dt(qgs) + d/dt(qgd) )

  //  the current *out of* the source should be simply
  // Isource

  // current out of drain' is
  // -Idrain - dtype*( d/dt(qgd) -  cdrain )

  // the current out of the source' is
  //  -Isource - dtype*( d/dt(qgs) +  cdrain )

  return true;
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
//    from updateSecondaryState:
//
//    ggd = ggd + capgd*(getSolverState().pdt);
//    ggs = ggs + capgs*(getSolverState().pdt);
//
//    // Sum the capacitor currents into the DC currents.
//    cg = cg + cqgs + cqgd;
//    cd  = cd  - cqgd;
//    cgd = cgd + cqgd;
//
//    So:
//
//    replace ggd with capgd.
//    replace ggs with capgs
//
//    replace cg with qgs+qgd
//    replace cd with -qgd
//    replace cgd with qgd.
//
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/01/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  double * qVec = extData.daeQVectorRawPtr;
  double * dQdxdVp = extData.dQdxdVpVectorRawPtr;

  // set up the final load variables:
  int Dtype = model_.dtype;
  double ceqgd = Dtype*(qgd);
  double ceqgs = Dtype*(((qgs+qgd)-qgd));
  double cdreq = Dtype*(((-qgd)+qgd));

  double ceqgd_Jdxp = -Dtype*(capgd*(vgd-vgd_orig));
  double ceqgs_Jdxp = -Dtype*(capgs*(vgs-vgs_orig));
  double cdreq_Jdxp = 0.0;

  qVec[li_Gate       ] += ( ceqgs+ceqgd);
  qVec[li_DrainPrime ] -= (-cdreq+ceqgd);
  qVec[li_SourcePrime] -= ( cdreq+ceqgs);

  if (!origFlag)
  {
    dQdxdVp[li_Gate       ] -= ( ceqgs_Jdxp+ceqgd_Jdxp);
    dQdxdVp[li_DrainPrime ] += (-cdreq_Jdxp+ceqgd_Jdxp);
    dQdxdVp[li_SourcePrime] += ( cdreq_Jdxp+ceqgs_Jdxp);
  }
  
  if( loadLeadCurrent )
  {
    double * leadQ = extData.nextLeadCurrQCompRawPtr;
    if (drainCond == 0.0)
    {
      leadQ[li_branch_dev_id] = -(-cdreq+ceqgd);
    }
    if (sourceCond == 0.0)
    {
      leadQ[li_branch_dev_is] = -( cdreq+ceqgs);
    }
    leadQ[li_branch_dev_ig] = (ceqgs+ceqgd);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 MESFET instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/01/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * solVec = extData.nextSolVectorRawPtr;
  double * fVec = extData.daeFVectorRawPtr;
  double * dFdxdVp = extData.dFdxdVpVectorRawPtr;

  // set up the final load variables:
  int Dtype = model_.dtype;
  double ceqgd = Dtype*(cgd);
  double ceqgs = Dtype*((cg-cgd));
  double cdreq = Dtype*((cd+cgd));

  double ceqgd_Jdxp = -Dtype*(ggd*(vgd-vgd_orig));
  double ceqgs_Jdxp = -Dtype*(ggs*(vgs-vgs_orig));
  double cdreq_Jdxp = -Dtype*(gds*(vds-vds_orig)+gm*(vgs-vgs_orig));

  // optional load resistors:
  if (drainCond  != 0.0)
  {
    fVec[li_Drain ] += Idrain;
  }
  if (sourceCond != 0.0)
  {
    fVec[li_Source] += Isource;
  }

  fVec[li_Gate       ] += (ceqgs+ceqgd);
  fVec[li_DrainPrime ] -= (Idrain +(-cdreq+ceqgd));
  fVec[li_SourcePrime] -= (Isource+(cdreq+ceqgs));

  if (!origFlag)
  {
    dFdxdVp[li_Gate       ] -= ( ceqgs_Jdxp+ceqgd_Jdxp);
    dFdxdVp[li_DrainPrime ] += (-cdreq_Jdxp+ceqgd_Jdxp);
    dFdxdVp[li_SourcePrime] += ( cdreq_Jdxp+ceqgs_Jdxp);
  }
  
  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    double * junctionV = extData.nextJunctionVCompRawPtr;
    if (drainCond != 0.0)
    {
      leadF[li_branch_dev_id] = Idrain;
    }
    else
    {
      leadF[li_branch_dev_id] = -(Idrain +(-cdreq+ceqgd));
    }
    if (sourceCond != 0.0)
    {
      leadF[li_branch_dev_is] = Isource;
    }
    else
    {
      leadF[li_branch_dev_is] = -(Isource+(cdreq+ceqgs));
    }
    leadF[li_branch_dev_ig] = (ceqgs+ceqgd);

    junctionV[li_branch_dev_id] = solVec[li_Drain] - solVec[li_Source];
    junctionV[li_branch_dev_ig] = solVec[li_Gate] - solVec[li_Source];
    junctionV[li_branch_dev_is] = 0.0;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 MESFET instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/01/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);

  dQdx[li_Gate       ][AGateEquGateNodeOffset        ] += capgd+capgs;
  dQdx[li_Gate       ][AGateEquDrainPrimeNodeOffset        ] -= capgd;
  dQdx[li_Gate       ][AGateEquSourcePrimeNodeOffset       ] -= capgs;
  dQdx[li_DrainPrime ][ADrainPrimeEquGateNodeOffset        ] -= capgd;
  dQdx[li_DrainPrime ][ADrainPrimeEquDrainPrimeNodeOffset  ] += capgd;
  dQdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset       ] -= capgs;
  dQdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset] += capgs;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 resistor  instance.
//
// Special Notes : The F-vector is an algebraic constaint.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/01/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Drain][ADrainEquDrainNodeOffset] += drainCond;
  dFdx[li_Drain][ADrainEquDrainPrimeNodeOffset] -= drainCond;

  dFdx[li_Gate][AGateEquGateNodeOffset] += ggd+ggs;
  dFdx[li_Gate][AGateEquDrainPrimeNodeOffset] -= ggd;
  dFdx[li_Gate][AGateEquSourcePrimeNodeOffset] -= ggs;

  dFdx[li_Source][ASourceEquSourceNodeOffset] += sourceCond;
  dFdx[li_Source][ASourceEquSourcePrimeNodeOffset] -= sourceCond;

  dFdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset] -= drainCond;
  dFdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset] += gm-ggd;
  dFdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset] +=
     drainCond+gds+ggd;
  dFdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset] += -gds-gm;

  dFdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset] -= gm+ggs;
  dFdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset] -= sourceCond;
  dFdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset] -= gds;
  dFdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]
     += sourceCond+gds+gm+ggs;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp_tmp)
{
  bool bsuccess = true;
  double tnom, ratio;
  double arg, arg1;
  double ratio1;
  double fact1, fact2;
  double kt, kt1;
  double vtnom;
  double egfet, egfet1;
  double pbfact;
  double cjfact, cjfact1;
  double gmanew, gmaold;
  double pbo;
  double xfc;
  double Pb;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
//    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "  Instance::Begin of updateTemperature. \n";
    Xyce::dout() << "  name = " << getName() << std::endl;
    Xyce::dout() << std::endl;
  }

  // first set the instance temperature to the new temperature:
  if (temp_tmp != -999.0) temp = temp_tmp;
  if (model_.interpolateTNOM(temp))
  {
    // make sure interpolation doesn't take any resistance negative
    if(model_.RD < 0) model_.RD = 0;
    if(model_.RS < 0) model_.RS = 0;

    // some params may have changed during interpolation
    // model_.processParams();
  }

  Pb   = model_.PB;
  tnom = model_.TNOM;
  ratio = temp/tnom;

  //  first do the model stuff

  vtnom  = tnom*CONSTKoverQ;
  fact1  = tnom/CONSTREFTEMP;
  kt1    = CONSTboltz*tnom;
  egfet1 = 1.16 - (7.02e-4*tnom*tnom)/(tnom + 1108);
  arg1   = -egfet1/(2.0*kt1) + 1.1150877/(CONSTboltz*2.0*CONSTREFTEMP);
  pbfact = -2.0*vtnom*(1.5*log(fact1) + CONSTQ*arg1);
  pbo    = (Pb - pbfact)/fact1;
  gmaold = (Pb - pbo)/pbo;
  cjfact = 1.0/(1.0 + 0.5*(4e-4*(tnom - CONSTREFTEMP) - gmaold));

  if(model_.FC >.95) {
      Xyce::dout() << "Depletion cap. coeff. FC too large, limited to .95";
      Xyce::dout() << std::endl;
      model_.FC = .95;
  }
  xfc = log(1.0 - model_.FC);
  f2  = exp(1.5*xfc);
  f3  = 1.0 - 1.5*model_.FC;
  // skip  bFac

  //  now do the instance stuff

  vt = temp*CONSTKoverQ;
  kt = temp*CONSTboltz;
  fact2 = temp/CONSTREFTEMP;
  ratio1 = ratio - 1.0;
  tIS = model_.IS*exp(ratio1*1.11/vt)*area;

  tCGS  = model_.CGS*cjfact*area;
  tCGD  = model_.CGD*cjfact*area;
  egfet = 1.16 - (7.02e-4*temp*temp)/(temp + 1108);
  arg   = -egfet/(2.0*kt) + 1.1150877/(CONSTboltz*2.0*CONSTREFTEMP);
  pbfact = -2.0*vt*(1.5*log(fact2) + CONSTQ*arg);
  tPB    = fact2*pbo + pbfact;
  gmanew = (tPB - pbo)/pbo;
  cjfact1 = 1.0 + 0.5*(4e-4*(temp - CONSTREFTEMP) - gmanew);
  tCGS *= cjfact1;
  tCGD *= cjfact1;

  corDepCap = model_.FC*tPB;
  f1    = tPB*(1.0 - exp((0.5)*xfc))/(0.5);
  vcrit = vt * log(vt/(CONSTroot2 * tIS));

  // the following parameters have no temperature dependence in Spice 3f5
  //
  tBeta   =  model_.BETA*area;   // transconductance parameter
  tvt0    =  model_.VTO;         // threshold voltage
  tLambda =  model_.LAMBDA;      // channel-length modulation
  tAlpha  =  model_.ALPHA;       // saturation voltage parameter
  tRD     =  model_.RD/area;     // drain ohmic resistance
  tRS     =  model_.RS/area;     // source ohmic resistance
  tMESb   =  model_.B;           // dopinng tail parameter

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "temp   = "<< temp << std::endl;
    Xyce::dout() << "tnom   = " << tnom << std::endl;
    Xyce::dout() << "ratio  = " << ratio << std::endl;
    Xyce::dout() << "vt     = " << vt << std::endl;
    Xyce::dout() << "kt     = " << kt << std::endl;
    Xyce::dout() << "fact2  = " << fact2 << std::endl;
    Xyce::dout() << "egfet  = " << egfet << std::endl;
    Xyce::dout() << "arg    = " << arg << std::endl;
    Xyce::dout() << "pbfact = " << pbfact << std::endl;
    Xyce::dout() << "PB     = " << Pb << std::endl;
    Xyce::dout() << "pbo    = " << pbo << std::endl;
    Xyce::dout() << "f2     = " << f2 << std::endl;
    Xyce::dout() << "f3     = " << f3 << std::endl;
    Xyce::dout() << "corDepCap= " << corDepCap << std::endl;
    Xyce::dout() << "tBeta   = " << tBeta << std::endl;
    Xyce::dout() << "tvt0    = " << tvt0 << std::endl;
    Xyce::dout() << "tPB     = " << tPB << std::endl;
    Xyce::dout() << "tMESb   = " << tMESb << std::endl;
    Xyce::dout() << "tLambda = " << tLambda << std::endl;
    //Xyce::dout() << "tTheta  = " << tTheta << std::endl;
    Xyce::dout() << "tRD     = " << tRD << std::endl;
    Xyce::dout() << "tRS     = " << tRS << std::endl;
    Xyce::dout() << " " << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  // process source/drain series resistance
  drainCond = 0;
  if (model_.RD != 0)
    drainCond = area/model_.RD;
  sourceCond = 0;
  if (model_.RS != 0)
    sourceCond = area/model_.RS;

  updateTemperature(temp);

  return true;
}

// MESFET Master functions:

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
    Instance & ji = *(*it);

    bool btmp = ji.updateIntermediateVars ();
    bsuccess = bsuccess && btmp;

    double * stoVec = ji.extData.nextStoVectorRawPtr;
    stoVec[ji.li_store_vgs] = ji.vgs;
    stoVec[ji.li_store_vgd] = ji.vgd;
    staVec[ji.li_state_qgs] = ji.qgs;
    staVec[ji.li_state_qgd] = ji.qgd;
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
    Instance & ji = *(*it);

    // F-vector:
    double * dFdxdVp = ji.extData.dFdxdVpVectorRawPtr;

    // set up the final load variables:
    int Dtype = ji.getModel().dtype;
    double f_ceqgd = Dtype*(ji.cgd);
    double f_ceqgs = Dtype*((ji.cg-ji.cgd));
    double f_cdreq = Dtype*((ji.cd+ji.cgd));

    double f_ceqgd_Jdxp = -Dtype*(ji.ggd*(ji.vgd-ji.vgd_orig));
    double f_ceqgs_Jdxp = -Dtype*(ji.ggs*(ji.vgs-ji.vgs_orig));
    double f_cdreq_Jdxp = -Dtype*(ji.gds*(ji.vds-ji.vds_orig)+ji.gm*(ji.vgs-ji.vgs_orig));

    // optional load resistors:
    if (ji.drainCond  != 0.0)
    {
      fVec[ji.li_Drain ] += ji.Idrain;
    }
    if (ji.sourceCond != 0.0)
    {
      fVec[ji.li_Source] += ji.Isource;
    }
    fVec[ji.li_Gate       ] += (f_ceqgs+f_ceqgd);
    fVec[ji.li_DrainPrime ] -= (ji.Idrain +(-f_cdreq+f_ceqgd));
    fVec[ji.li_SourcePrime] -= (ji.Isource+(f_cdreq+f_ceqgs));

    if (!ji.origFlag)
    {
      dFdxdVp[ji.li_Gate       ] -= ( f_ceqgs_Jdxp+f_ceqgd_Jdxp);
      dFdxdVp[ji.li_DrainPrime ] += (-f_cdreq_Jdxp+f_ceqgd_Jdxp);
      dFdxdVp[ji.li_SourcePrime] += ( f_cdreq_Jdxp+f_ceqgs_Jdxp);
    }

    // Q-vector:
    double * dQdxdVp = ji.extData.dQdxdVpVectorRawPtr;

    // set up the final load variables:
    double q_ceqgd = Dtype*(ji.qgd);
    double q_ceqgs = Dtype*(((ji.qgs+ji.qgd)-ji.qgd));
    double q_cdreq = Dtype*(((-ji.qgd)+ji.qgd));

    double q_ceqgd_Jdxp = -Dtype*(ji.capgd*(ji.vgd-ji.vgd_orig));
    double q_ceqgs_Jdxp = -Dtype*(ji.capgs*(ji.vgs-ji.vgs_orig));
    double q_cdreq_Jdxp = 0.0;

    qVec[ji.li_Gate       ] += ( q_ceqgs+q_ceqgd);
    qVec[ji.li_DrainPrime ] -= (-q_cdreq+q_ceqgd);
    qVec[ji.li_SourcePrime] -= ( q_cdreq+q_ceqgs);

    if (!ji.origFlag)
    {
      dQdxdVp[ji.li_Gate       ] -= ( q_ceqgs_Jdxp+q_ceqgd_Jdxp);
      dQdxdVp[ji.li_DrainPrime ] += (-q_cdreq_Jdxp+q_ceqgd_Jdxp);
      dQdxdVp[ji.li_SourcePrime] += ( q_cdreq_Jdxp+q_ceqgs_Jdxp);
    }

    if( ji.loadLeadCurrent )
    {
      if (ji.drainCond != 0.0)
      {
        leadF[ji.li_branch_dev_id] = ji.Idrain;
      }
      else
      {
        leadF[ji.li_branch_dev_id] = -(ji.Idrain +(-f_cdreq+f_ceqgd));
        leadQ[ji.li_branch_dev_id] = -(-q_cdreq+q_ceqgd);
      }
      if (ji.sourceCond != 0.0)
      {
        leadF[ji.li_branch_dev_is] = ji.Isource;
      }
      else
      {
        leadF[ji.li_branch_dev_is] = -(ji.Isource+(f_cdreq+f_ceqgs));
        leadQ[ji.li_branch_dev_is] = -( q_cdreq+q_ceqgs);
      }
      leadF[ji.li_branch_dev_ig] = (f_ceqgs+f_ceqgd);
      leadQ[ji.li_branch_dev_ig] = (q_ceqgs+q_ceqgd);

      junctionV[ji.li_branch_dev_id] = solVec[ji.li_Drain] - solVec[ji.li_Source];
      junctionV[ji.li_branch_dev_ig] = solVec[ji.li_Gate] - solVec[ji.li_Source];
      junctionV[ji.li_branch_dev_is] = 0.0;
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
// Creation Date : 12/12/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & ji = *(*it);

    // F-matrix:

    *ji.f_DrainEquDrainNodePtr += ji.drainCond;

    *ji.f_DrainEquDrainPrimeNodePtr -= ji.drainCond;


    *ji.f_GateEquGateNodePtr += ji.ggd+ji.ggs;

    *ji.f_GateEquDrainPrimeNodePtr -= ji.ggd;

    *ji.f_GateEquSourcePrimeNodePtr -= ji.ggs;


    *ji.f_SourceEquSourceNodePtr += ji.sourceCond;

    *ji.f_SourceEquSourcePrimeNodePtr -= ji.sourceCond;


    *ji.f_DrainPrimeEquDrainNodePtr -= ji.drainCond;

    *ji.f_DrainPrimeEquGateNodePtr += ji.gm-ji.ggd;

    *ji.f_DrainPrimeEquDrainPrimeNodePtr += ji.drainCond+ji.gds+ji.ggd;

    *ji.f_DrainPrimeEquSourcePrimeNodePtr += -ji.gds-ji.gm;


    *ji.f_SourcePrimeEquGateNodePtr -= ji.gm+ji.ggs;

    *ji.f_SourcePrimeEquSourceNodePtr -= ji.sourceCond;

    *ji.f_SourcePrimeEquDrainPrimeNodePtr -= ji.gds;

    *ji.f_SourcePrimeEquSourcePrimeNodePtr += ji.sourceCond+ji.gds+ji.gm+ji.ggs;

    // Q-matrix:

    *ji.q_GateEquGateNodePtr         += ji.capgd+ji.capgs;

    *ji.q_GateEquDrainPrimeNodePtr         -= ji.capgd;

    *ji.q_GateEquSourcePrimeNodePtr        -= ji.capgs;

    *ji.q_DrainPrimeEquGateNodePtr         -= ji.capgd;

    *ji.q_DrainPrimeEquDrainPrimeNodePtr   += ji.capgd;

    *ji.q_SourcePrimeEquGateNodePtr        -= ji.capgs;

    *ji.q_SourcePrimeEquSourcePrimeNodePtr += ji.capgs;
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
// Creation Date : 12/12/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & ji = *(*it);

    // dFdx matrix:

    dFdx[ji.li_Drain][ji.ADrainEquDrainNodeOffset] += ji.drainCond;

    dFdx[ji.li_Drain][ji.ADrainEquDrainPrimeNodeOffset] -= ji.drainCond;


    dFdx[ji.li_Gate][ji.AGateEquGateNodeOffset] += ji.ggd+ji.ggs;

    dFdx[ji.li_Gate][ji.AGateEquDrainPrimeNodeOffset] -= ji.ggd;

    dFdx[ji.li_Gate][ji.AGateEquSourcePrimeNodeOffset] -= ji.ggs;


    dFdx[ji.li_Source][ji.ASourceEquSourceNodeOffset] += ji.sourceCond;

    dFdx[ji.li_Source][ji.ASourceEquSourcePrimeNodeOffset] -= ji.sourceCond;


    dFdx[ji.li_DrainPrime][ji.ADrainPrimeEquDrainNodeOffset] -= ji.drainCond;

    dFdx[ji.li_DrainPrime][ji.ADrainPrimeEquGateNodeOffset] += ji.gm-ji.ggd;

    dFdx[ji.li_DrainPrime][ji.ADrainPrimeEquDrainPrimeNodeOffset] +=
      ji.drainCond+ji.gds+ji.ggd;

    dFdx[ji.li_DrainPrime][ji.ADrainPrimeEquSourcePrimeNodeOffset] += -ji.gds-ji.gm;


    dFdx[ji.li_SourcePrime][ji.ASourcePrimeEquGateNodeOffset] -= ji.gm+ji.ggs;

    dFdx[ji.li_SourcePrime][ji.ASourcePrimeEquSourceNodeOffset] -= ji.sourceCond;

    dFdx[ji.li_SourcePrime][ji.ASourcePrimeEquDrainPrimeNodeOffset] -= ji.gds;

    dFdx[ji.li_SourcePrime][ji.ASourcePrimeEquSourcePrimeNodeOffset]
      += ji.sourceCond+ji.gds+ji.gm+ji.ggs;

    // dQdx matrix:

    dQdx[ji.li_Gate       ][ji.AGateEquGateNodeOffset        ] += ji.capgd+ji.capgs;

    dQdx[ji.li_Gate       ][ji.AGateEquDrainPrimeNodeOffset        ] -= ji.capgd;

    dQdx[ji.li_Gate       ][ji.AGateEquSourcePrimeNodeOffset       ] -= ji.capgs;

    dQdx[ji.li_DrainPrime ][ji.ADrainPrimeEquGateNodeOffset        ] -= ji.capgd;

    dQdx[ji.li_DrainPrime ][ji.ADrainPrimeEquDrainPrimeNodeOffset  ] += ji.capgd;

    dQdx[ji.li_SourcePrime][ji.ASourcePrimeEquGateNodeOffset       ] -= ji.capgs;

    dQdx[ji.li_SourcePrime][ji.ASourcePrimeEquSourcePrimeNodeOffset] += ji.capgs;
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
  if (deviceMap.empty() || (deviceMap.find("Z")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("z", 1)
      .registerModelType("nmf", 1)
      .registerModelType("pmf", 1);
  }
}

} // namespace MESFET
} // namespace Device
} // namespace Xyce

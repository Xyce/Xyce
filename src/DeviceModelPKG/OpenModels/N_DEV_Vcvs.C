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

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Vcvs.h>
#include <N_DEV_Bsrc.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {

namespace Vcvs {

void Traits::loadInstanceParameters(ParametricData<Vcvs::Instance> &p)
{
// Set up double precision variables:
   p.addPar ("G",0.0,&Vcvs::Instance::Gain)
    .setUnit(U_NONE)
    .setCategory(CAT_NONE)
    .setDescription("Gain");
}

void Traits::loadModelParameters(ParametricData<Vcvs::Model> &p)
{
}

std::vector< std::vector<int> > Instance::jacStamp;
std::vector< std::vector<int> > Instance::jacStampPDE;

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 12/20/00
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       Viter,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Viter),
    Gain(1.0),
    li_Pos(-1),
    li_Neg(-1),
    li_Bra(-1),
    li_ContPos(-1),
    li_ContNeg(-1),
    li_branch_data(0),
    ABraEquPosNodeOffset(-1),
    ABraEquNegNodeOffset(-1),
    ABraEquContPosNodeOffset(-1),
    ABraEquContNegNodeOffset(-1),
    APosEquBraVarOffset(-1),
    ANegEquBraVarOffset(-1),

    f_BraEquPosNodePtr(0),
    f_BraEquNegNodePtr(0),
    f_BraEquContPosNodePtr(0),
    f_BraEquContNegNodePtr(0),
    f_PosEquBraVarPtr(0),
    f_NegEquBraVarPtr(0)
{
  numIntVars   = 1;
  numExtVars   = 4;
  numStateVars = 0;
  setNumBranchDataVars(0);          // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1; // this is the space to allocate if lead current or power is needed.

  devConMap.resize(4);
  devConMap[0] = 1;
  devConMap[1] = 1;
  devConMap[2] = 2;
  devConMap[3] = 2;

  if( jacStamp.empty() )
  {
    jacStamp.resize(5);
    jacStamp[0].resize(1);
    jacStamp[0][0]=4;
    jacStamp[1].resize(1);
    jacStamp[1][0]=4;
    jacStamp[4].resize(4);
    jacStamp[4][0]=0;
    jacStamp[4][1]=1;
    jacStamp[4][2]=2;
    jacStamp[4][3]=3;

    // PDE supporting stampPDE.  This includes diagonal elements, needed by the
    // 2-level Newton.
    jacStampPDE.resize(5);
    jacStampPDE[0].resize(2);
    jacStampPDE[0][0]=0;
    jacStampPDE[0][1]=4;
    jacStampPDE[1].resize(2);
    jacStampPDE[1][0]=1;
    jacStampPDE[1][1]=4;

    jacStampPDE[2].resize(1);
    jacStampPDE[1][0]=2;

    jacStampPDE[3].resize(1);
    jacStampPDE[1][0]=3;

    jacStampPDE[4].resize(5);
    jacStampPDE[4][0]=0;
    jacStampPDE[4][1]=1;
    jacStampPDE[4][2]=2;
    jacStampPDE[4][3]=3;
    jacStampPDE[4][4]=4;
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (!given("G"))
  {
    UserError(*this) << "Could not find Gain parameter in instance.";
  }
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
// Creator       : Robert Hoekstra, SNL, Computational Science
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs ( const std::vector<int> & intLIDVecRef,
                                        const std::vector<int> & extLIDVecRef )
{
  std::string msg;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  VcvsInstance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

  if (numInt != numIntVars)
  {
    DevelFatal(*this).in("Instance::registerLIDs")
      <<  "numInt != numIntVars";
  }

  if (numExt != numExtVars)
  {
    DevelFatal(*this).in("Instance::registerLIDs")
      <<  "numExt != numExtVars";
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];
  li_ContPos = extLIDVec[2];
  li_ContNeg = extLIDVec[3];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl;
    Xyce::dout() << "  li_Neg = " << li_Neg << std::endl;
    Xyce::dout() << "  li_ContPos = " << li_ContPos << std::endl;
    Xyce::dout() << "  li_ContNeg = " << li_ContNeg << std::endl;
  }

  li_Bra = intLIDVec[0];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    Xyce::dout() << "  li_Bra = " << li_Bra << std::endl;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    Xyce::dout() << section_divider << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerBranchDataLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 12/18/2012
//-----------------------------------------------------------------------------
/// Register the local store IDs
///
/// In addition to state vector, Xyce maintains a separate datastructure
/// called a "branch data" vector.  As with other such vectors, the device
/// declares at construction time how many branch vector entries it needs,
/// and later Topology assigns locations for devices, returning LIDs.
///
/// These LIDs are stored in this method for later use.
///
/// The Resistor device uses exactly one "branch data vector" element, where
///
/// @param stoLIDVecRef Store variable local IDs
///
/// @author Richard Schiek, Electrical Systems Modeling
/// @date   12/18/2012

void Instance::registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef)
{
  AssertLIDs(branchLIDVecRef.size() == getNumBranchDataVars());

  if (loadLeadCurrent)
  {
    li_branch_data= branchLIDVecRef[0];
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
  addInternalNode(symbol_table, li_Bra, getName(), "branch");

  if (loadLeadCurrent)
  {
    addBranchDataNode( symbol_table, li_branch_data, getName(), "BRANCH_D");
  }
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
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/2/02
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  if (getSolverState().isPDESystem_)
    return jacStampPDE;

  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/2/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  if (getSolverState().isPDESystem_)
  {
    APosEquBraVarOffset = jacLIDVec[0][1];
    ANegEquBraVarOffset = jacLIDVec[1][1];
    ABraEquPosNodeOffset = jacLIDVec[4][0];
    ABraEquNegNodeOffset = jacLIDVec[4][1];
    ABraEquContPosNodeOffset = jacLIDVec[4][2];
    ABraEquContNegNodeOffset = jacLIDVec[4][3];
  }
  else
  {
    APosEquBraVarOffset = jacLIDVec[0][0];
    ANegEquBraVarOffset = jacLIDVec[1][0];
    ABraEquPosNodeOffset = jacLIDVec[4][0];
    ABraEquNegNodeOffset = jacLIDVec[4][1];
    ABraEquContPosNodeOffset = jacLIDVec[4][2];
    ABraEquContNegNodeOffset = jacLIDVec[4][3];
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/12/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  f_PosEquBraVarPtr = &(dFdx[li_Pos][ APosEquBraVarOffset ]);
  f_NegEquBraVarPtr = &(dFdx[li_Neg][ ANegEquBraVarOffset ]);
  f_BraEquPosNodePtr = &(dFdx[li_Bra][ ABraEquPosNodeOffset ]);
  f_BraEquNegNodePtr = &(dFdx[li_Bra][ ABraEquNegNodeOffset ]);
  f_BraEquContPosNodePtr = &(dFdx[li_Bra][ ABraEquContPosNodeOffset ]);
  f_BraEquContNegNodePtr = &(dFdx[li_Bra][ ABraEquContNegNodeOffset ]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/29/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 Vcvs instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/19/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;
  double * solVec = extData.nextSolVectorRawPtr;

  double v_pos = solVec[li_Pos];
  double v_neg = solVec[li_Neg];
  double v_cont_pos = solVec[li_ContPos];
  double v_cont_neg = solVec[li_ContNeg];

  double i_bra = solVec[li_Bra];

  double c_tmp = i_bra;
  double v_tmp = (v_pos-v_neg);

  fVec[li_Pos] += i_bra;
  fVec[li_Neg] += -i_bra;

  double src = Gain * ( v_cont_pos - v_cont_neg ) - v_pos + v_neg;
  fVec[li_Bra] += -src;

  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    leadF[li_branch_data] = c_tmp;
    double * junctionV = extData.nextJunctionVCompRawPtr;
    junctionV[li_branch_data] = v_tmp; 
  }

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
// Creation Date : 02/19/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Pos][APosEquBraVarOffset] += 1.0;
  dFdx[li_Neg][ANegEquBraVarOffset] -= 1.0;

  dFdx[li_Bra][ABraEquPosNodeOffset] += 1.0;
  dFdx[li_Bra][ABraEquNegNodeOffset] -= 1.0;
  dFdx[li_Bra][ABraEquContPosNodeOffset] -= Gain;
  dFdx[li_Bra][ABraEquContNegNodeOffset] += Gain;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/17/04
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  varTypeVec.resize(1);
  varTypeVec[0] = 'I';
}

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : "Model block" constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 12/21/00
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block)
{
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
    os << "  " << i << ": " << (*iter)->getName() << "      ";
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


// Vcvs Master functions:

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
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, 
                             double * leadF, double * leadQ, double * junctionV, int loadType)
{
  InstanceVector::const_iterator it, end;

  if (loadType == LINEAR_FREQ)
    loadType = LINEAR;

  if (!separateInstances_ && ( loadType == LINEAR || loadType == NONLINEAR ))
  {
    separateInstanceTypes(linearInstances_, nonlinearInstances_);
    separateInstances_ = true;
  }

  if (loadType == ALL)
  {
    it = getInstanceBegin();
    end = getInstanceEnd();
  }
  else if (loadType == LINEAR)
  {
    it = linearInstances_.begin();
    end = linearInstances_.end();
  }
  else
  {
    it = nonlinearInstances_.begin();
    end = nonlinearInstances_.end();
  }

  for ( ; it != end; ++it )
  {
    Instance & vi = *(*it);

    double v_pos = solVec[vi.li_Pos];
    double v_neg = solVec[vi.li_Neg];
    double v_cont_pos = solVec[vi.li_ContPos];
    double v_cont_neg = solVec[vi.li_ContNeg];

    double i_bra = solVec[vi.li_Bra];

    fVec[vi.li_Pos] += i_bra;
    fVec[vi.li_Neg] += -i_bra;

    double src = vi.Gain * ( v_cont_pos - v_cont_neg ) - v_pos + v_neg;
    fVec[vi.li_Bra] += -src;

    if( vi.loadLeadCurrent )
      {
        leadF[vi.li_branch_data] = i_bra;
        junctionV[vi.li_branch_data] = v_pos-v_neg; 
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
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx, int loadType)
{
  InstanceVector::const_iterator it, end;

  if (loadType == LINEAR_FREQ)
    loadType = LINEAR;

  if (!separateInstances_ && ( loadType == LINEAR || loadType == NONLINEAR ))
  {
    separateInstanceTypes(linearInstances_, nonlinearInstances_);
    separateInstances_ = true;
  }

  if (loadType == ALL)
  {
    it = getInstanceBegin();
    end = getInstanceEnd();
  }
  else if (loadType == LINEAR)
  {
    it = linearInstances_.begin();
    end = linearInstances_.end();
  }
  else
  {
    it = nonlinearInstances_.begin();
    end = nonlinearInstances_.end();
  }

  for ( ; it != end; ++it )
  {
    Instance & vi = *(*it);

#ifndef Xyce_NONPOINTER_MATRIX_LOAD

    *vi.f_PosEquBraVarPtr += 1.0;
    *vi.f_NegEquBraVarPtr -= 1.0;

    *vi.f_BraEquPosNodePtr += 1.0;
    *vi.f_BraEquNegNodePtr -= 1.0;
    *vi.f_BraEquContPosNodePtr -= vi.Gain;
    *vi.f_BraEquContNegNodePtr += vi.Gain;
#else

    dFdx[vi.li_Pos][vi.APosEquBraVarOffset] += 1.0;
    dFdx[vi.li_Neg][vi.ANegEquBraVarOffset] -= 1.0;

    dFdx[vi.li_Bra][vi.ABraEquPosNodeOffset] += 1.0;
    dFdx[vi.li_Bra][vi.ABraEquNegNodeOffset] -= 1.0;
    dFdx[vi.li_Bra][vi.ABraEquContPosNodeOffset] -= vi.Gain;
    dFdx[vi.li_Bra][vi.ABraEquContNegNodeOffset] += vi.Gain;
#endif
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
  if (deviceMap.empty() || (deviceMap.find("E")!=deviceMap.end()))
  {
    Bsrc::registerDevice();

    Config<Traits>::addConfiguration()
      .registerDevice("e", 1)
      .registerModelType("e", 1);
  }
}

} // namespace Vcvs
} // namespace Device
} // namespace Xyce

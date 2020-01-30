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
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
//
// Creation Date : 2/1/10
//
//
//
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

#include <memory>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Resistor3.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_DEV_Resistor.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {
namespace Resistor3 {

void Traits::loadInstanceParameters(ParametricData<Resistor3::Instance> &p)
{
  p.addPar ("R",1000.0,&Resistor3::Instance::R)
   .setExpressionAccess(ParameterType::TIME_DEP)
   .setUnit(U_OHM)
   .setCategory(CAT_NONE)
   .setDescription("Resistance");

   p.addPar("M", 1.0, &Resistor3::Instance::multiplicityFactor)
    .setUnit(U_NONE)
    .setDescription("Multiplicity Factor");

  p.addPar ("L",0.0,&Resistor3::Instance::length)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Length");

  p.addPar ("W",0.0,&Resistor3::Instance::width)
   .setUnit(U_METER)
   .setCategory(CAT_NONE)
   .setDescription("Width");

  p.addPar ("TEMP",0.0,&Resistor3::Instance::temp)
   .setExpressionAccess(ParameterType::TIME_DEP)
   .setUnit(U_DEGC)
   .setCategory(CAT_NONE)
   .setDescription("Device temperature");

  p.addPar ("TC1",0.0,&Resistor3::Instance::tempCoeff1)
   .setGivenMember(&Resistor3::Instance::tempCoeff1Given)
   .setUnit(U_DEGCM1)
   .setCategory(CAT_NONE)
   .setDescription("Linear Temperature Coefficient");

  p.addPar ("TC2",0.0,&Resistor3::Instance::tempCoeff2)
   .setGivenMember(&Resistor3::Instance::tempCoeff2Given)
   .setUnit(U_DEGCM2)
   .setCategory(CAT_NONE)
   .setDescription("Quadratic Temperature Coefficient");

  // This call tells the parameter handling code that TC can be specified
  // as a vector with up to two elements as in TC=a,b.  It then translates
  // TC=a,b into TC1=a TC2=b.  Likewise,TC=a will translate into TC1=a
  p.makeVector ("TC",2);

  p.addPar("TCE", 0.0, &Resistor3::Instance::tempCoeffExp)
    .setGivenMember(&Resistor3::Instance::tempCoeffExpGiven)
    .setUnit(U_PERCENTDEGCM1)
    .setDescription("Exponential Temperature Coefficient");

  p.addPar ("DTEMP",0.0,&Resistor3::Instance::dtemp)
   .setGivenMember(&Resistor3::Instance::dtempGiven)
   .setUnit(U_DEGC)
   .setCategory(CAT_NONE)
   .setDescription("Device Temperature -- For compatibility only. Parameter is NOT used");
}

void Traits::loadModelParameters(ParametricData<Resistor3::Model> &p)
{}

std::vector< std::vector<int> > Instance::jacStamp;
std::vector< std::vector<int> > Instance::jacStampPDE;


// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       Viter,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Viter),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    tempCoeffExp(0.0),
    dtemp(0.0),
    tempCoeff1Given(false),
    tempCoeff2Given(false),
    tempCoeffExpGiven(false),
    dtempGiven(false),
    li_Pos(-1),
    li_Neg(-1),
    li_Bra(-1),

    ABraEquPosNodeOffset(-1),
    ABraEquNegNodeOffset(-1),
    APosEquBraVarOffset(-1),
    ANegEquBraVarOffset(-1),
    APosEquPosNodeOffset(-1),
    ANegEquNegNodeOffset(-1),
    ABraEquBraVarOffset(-1)

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  ,fBraEquPosNodePtr(0),
  fBraEquNegNodePtr(0),
  fPosEquBraVarPtr(0),
  fNegEquBraVarPtr(0),
  fPosEquPosNodePtr(0),
  fNegEquNegNodePtr(0),
  fBraEquBraVarPtr(0)
#endif
{
  numIntVars   = 1;
  numExtVars   = 2;
  numStateVars = 0;

  if( jacStamp.empty() )
  {
    jacStamp.resize(3);
    jacStamp[0].resize(1);
    jacStamp[0][0] = 2;
    jacStamp[1].resize(1);
    jacStamp[1][0] = 2;
    jacStamp[2].resize(2);
    jacStamp[2][0] = 0;
    jacStamp[2][1] = 1;

    // PDE supporting stamp.  This includes diagonal elements, needed by the
    // 2-level Newton.
    jacStampPDE.resize(3);
    jacStampPDE[0].resize(2);
    jacStampPDE[0][0] = 0;
    jacStampPDE[0][1] = 2;
    jacStampPDE[1].resize(2);
    jacStampPDE[1][0] = 1;
    jacStampPDE[1][1] = 2;
    jacStampPDE[2].resize(3);
    jacStampPDE[2][0] = 0;
    jacStampPDE[2][1] = 1;
    jacStampPDE[2][2] = 2;
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  processParams();

  // Calculate any parameters specified as expressions:
  updateDependentParameters();
  processParams();
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
void Instance::registerLIDs ( const std::vector<int> & intLIDVecRef,
	                                const std::vector<int> & extLIDVecRef)
{
  std::string msg;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  Instance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

  if (numInt != numIntVars)
  {
    DevelFatal(*this).in("Instance::registerLIDs")
      << "numInt != numIntVars";
  }

  if (numExt != numExtVars)
  {
    DevelFatal(*this).in("Instance::registerLIDs")
      << "numExt != numExtVars";
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];
  li_Bra = intLIDVec[0];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl;
    Xyce::dout() << "  li_Neg = " << li_Neg << std::endl;
    Xyce::dout() << "  li_Bra = " << li_Bra << std::endl;
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
  addInternalNode(symbol_table, li_Bra, getName(), "branch");
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
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
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
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
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  if (getSolverState().isPDESystem_)
  {
    APosEquBraVarOffset  = jacLIDVec[0][1];
    ANegEquBraVarOffset  = jacLIDVec[1][1];
    ABraEquPosNodeOffset = jacLIDVec[2][0];
    ABraEquNegNodeOffset = jacLIDVec[2][1];
  }
  else
  {
    APosEquBraVarOffset  = jacLIDVec[0][0];
    ANegEquBraVarOffset  = jacLIDVec[1][0];
    ABraEquPosNodeOffset = jacLIDVec[2][0];
    ABraEquNegNodeOffset = jacLIDVec[2][1];
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  fPosEquBraVarPtr = &(dFdx[li_Pos][APosEquBraVarOffset]);
  fNegEquBraVarPtr = &(dFdx[li_Neg][ANegEquBraVarOffset]);
  fBraEquPosNodePtr = &(dFdx[li_Bra][ABraEquPosNodeOffset]);
  fBraEquNegNodePtr = &(dFdx[li_Bra][ABraEquNegNodeOffset]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 vsrc instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  // Using values determined during the loadRHS function call.
  // (srcCurrent, srcVoltage).

  double * fVec = extData.daeFVectorRawPtr;
  double * solVec = extData.nextSolVectorRawPtr;

  fVec[li_Pos] += solVec[li_Bra];
  fVec[li_Neg] -= solVec[li_Bra];
  fVec[li_Bra] += (solVec[li_Pos]-solVec[li_Neg]);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 resistor  instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Pos][APosEquBraVarOffset] += 1.0;
  dFdx[li_Neg][ANegEquBraVarOffset] -= 1.0;
  dFdx[li_Bra][ABraEquPosNodeOffset] += 1.0;
  dFdx[li_Bra][ABraEquNegNodeOffset] -= 1.0;

  return true;
}

// end of new-DAE functions

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  varTypeVec.resize(1);
  varTypeVec[0] = 'I';
}

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    DC_TRAN (0)
{
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
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

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
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



//-----------------------------------------------------------------------------
// Resistor3 Master functions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, double * leadF, double * leadQ, double * junctionV)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & vi = *(*it);

    fVec[vi.li_Pos] += solVec[vi.li_Bra];
    fVec[vi.li_Neg] -= solVec[vi.li_Bra];
    fVec[vi.li_Bra] += (solVec[vi.li_Pos]-solVec[vi.li_Neg]);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & vi = *(*it);

#ifndef Xyce_NONPOINTER_MATRIX_LOAD

    *(vi.fPosEquBraVarPtr) += 1.0;

    *(vi.fNegEquBraVarPtr) -= 1.0;

    *(vi.fBraEquPosNodePtr) += 1.0;

    *(vi.fBraEquNegNodePtr) -= 1.0;
#else

    dFdx[vi.li_Pos][vi.APosEquBraVarOffset] += 1.0;

    dFdx[vi.li_Neg][vi.ANegEquBraVarOffset] -= 1.0;

    dFdx[vi.li_Bra][vi.ABraEquPosNodeOffset] += 1.0;

    dFdx[vi.li_Bra][vi.ABraEquNegNodeOffset] -= 1.0;
#endif

  }

  return true;
}

Device *
Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{
  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() || (deviceMap.find("R")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("r", 3)
      .registerModelType("r", 3);
  }
}

} // namespace Resistor3
} // namespace Device
} // namespace Xyce

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
#include <algorithm>
#include <fstream>
#include <set>
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_Inductor.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Expression.h>

namespace Xyce {
namespace Device {
namespace Inductor {

void Traits::loadInstanceParameters(ParametricData<Inductor::Instance> &p)
{
  p.addPar("L",    0.0, &Inductor::Instance::baseL)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_HENRY)
    .setDescription("Inductance")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&indSens);

  p.addPar("M", 1.0, &Inductor::Instance::multiplicityFactor)
    .setUnit(U_NONE)
    .setDescription("Multiplicity Factor");

  p.addPar("IC",   0.0, &Inductor::Instance::IC)
    .setGivenMember(&Inductor::Instance::ICGiven)
    .setUnit(U_AMP)
    .setDescription("Initial current through device");

  p.addPar("TEMP", 0.0, &Inductor::Instance::temp)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setGivenMember(&Inductor::Instance::tempGiven)
    .setUnit(U_DEGC)
    .setCategory(CAT_MATERIAL)
    .setDescription("Device temperature");

  p.addPar("TC1", 0.0, &Inductor::Instance::tempCoeff1)
    .setGivenMember(&Inductor::Instance::tempCoeff1Given)
    .setUnit(U_DEGCM1)
    .setDescription("Linear Temperature Coefficient");

  p.addPar("TC2", 0.0, &Inductor::Instance::tempCoeff2)
    .setGivenMember(&Inductor::Instance::tempCoeff2Given)
    .setUnit(U_DEGCM2)
    .setDescription("Quadratic Temperature Coefficient");

  // This call tells the parameter handling code that TC can be specified
  // as a vector with up to two elements as in TC=a,b.  It then translates
  // TC=a,b into TC1=a TC2=b.  Likewise, TC=a will translate into TC1=a
  p.makeVector ("TC", 2);
}

void Traits::loadModelParameters(ParametricData<Inductor::Model> &p)
{
  // Set up double precision variables:
  p.addPar("L", 1.0, &Inductor::Model::inductanceMultiplier)
    .setUnit(U_NONE)
    .setDescription("Inductance Multiplier");

  p.addPar("IC", 0.0, &Inductor::Model::IC)
    .setUnit(U_AMP)
    .setDescription("Initial current through device");

  p.addPar("TNOM", 27.0, &Inductor::Model::tnom)
    .setUnit(U_DEGC)
    .setCategory(CAT_MATERIAL)
    .setDescription("Reference temperature");

  p.addPar("TC1",0.0, &Inductor::Model::tempCoeff1)
    .setUnit(U_DEGCM1)
    .setCategory(CAT_MATERIAL)
    .setDescription("First order temperature coeff.");

  p.addPar("TC2", 0.0, &Inductor::Model::tempCoeff2)
    .setUnit(U_DEGCM2)
    .setCategory(CAT_MATERIAL)
    .setDescription("Second order temperature coeff.");
}

//
// static class member inits
//
std::vector< std::vector<int> > Instance::jacStamp_BASE;

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Instance::processParams()
{
  if (!given("TEMP"))
    temp = getDeviceOptions().temp.getImmutableValue<double>();

  if (!tempCoeff1Given)
    tempCoeff1=model_.tempCoeff1;
  if (!tempCoeff2Given)
    tempCoeff2=model_.tempCoeff2;

  // M must be non-negative
  if (multiplicityFactor <= 0)
  {
    UserError(*this) << "Multiplicity Factor (M) must be non-negative" << std::endl;
  }

  // If there are any time dependent parameters, set their values at for
  // the current time.

  // now set the temperature related stuff.
  updateTemperature(temp);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 02/27/01
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp)
{
  double difference = temp - model_.tnom;
  // Apply both the inductanceMultipler (from the model card) and the multiplicityFactor M 
  // (from the instance line).
  double factor = model_.inductanceMultiplier*(1.0 + tempCoeff1*difference +
                         tempCoeff2*difference*difference)/multiplicityFactor;
  L = baseL*factor;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & instance_block,
  Model &               model,
  const FactoryBlock &  factory_block)
  : DeviceInstance(instance_block, configuration.getInstanceParameters(), factory_block),
    model_(model),
    L(0),
    multiplicityFactor(1.0),
    IC(0),
    ICGiven(false),
    baseL(0.0),
    temp(getDeviceOptions().temp.getImmutableValue<double>()),
    tempGiven(0),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    tempCoeff1Given(false),
    tempCoeff2Given(false),
    li_fstate(-1),
    li_Pos(-1),
    li_Neg(-1),
    li_Bra(-1),
    li_branch_data(0),
    ABraEquPosNodeOffset(-1),
    ABraEquNegNodeOffset(-1),
    ABraEquBraVarOffset(-1),
    APosEquBraVarOffset(-1),
    ANegEquBraVarOffset(-1)
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  ,
    fPosEquBraVarPtr(0),
    fNegEquBraVarPtr(0),
    fBraEquPosNodePtr(0),
    fBraEquNegNodePtr(0),
    fBraEquBraVarPtr(0),
    qBraEquBraVarPtr(0)
#endif
{
  numExtVars   = 2;
  numIntVars   = 1;
  numStateVars = 1;
  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1;    // this is the space to allocate if lead current or power is needed.

  if( jacStamp_BASE.empty() )
  {
    jacStamp_BASE.resize(3);

    jacStamp_BASE[0].resize(1);
    jacStamp_BASE[0][0] = 2;

    jacStamp_BASE[1].resize(1);
    jacStamp_BASE[1][0] = 2;

    jacStamp_BASE[2].resize(3);
    jacStamp_BASE[2][0] = 0;
    jacStamp_BASE[2][1] = 1;
    jacStamp_BASE[2][2] = 2;
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (instance_block.params);

  // Set any non-constant parameter defaults:
  if (!given("L"))
  {
    UserError(*this) << "Could not find L parameter in instance.";
  }

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

  // set up numIntVars:
  numIntVars = 1;

  // set up numStateVars:
  numStateVars = 2;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::~Instance()
{
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

  fPosEquBraVarPtr  = &(dFdx[li_Pos][APosEquBraVarOffset]);
  fNegEquBraVarPtr  = &(dFdx[li_Neg][ANegEquBraVarOffset]);
  fBraEquPosNodePtr = &(dFdx[li_Bra][ABraEquPosNodeOffset]);
  fBraEquNegNodePtr = &(dFdx[li_Bra][ABraEquNegNodeOffset]);
  fBraEquBraVarPtr  = &(dFdx[li_Bra][ABraEquBraVarOffset]);

  qBraEquBraVarPtr = &(dQdx[li_Bra][ABraEquBraVarOffset]);

#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::isLinearDevice
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool Instance::isLinearDevice() const
{
  if( loadLeadCurrent )
  {
    return false;
  }

  const std::vector<Depend> & depVec = const_cast<Xyce::Device::Inductor::Instance*>(this)->getDependentParams();
  if ( depVec.size() )
  {
    std::vector<Depend>::const_iterator d;
    std::vector<Depend>::const_iterator begin=depVec.begin();
    std::vector<Depend>::const_iterator end=depVec.end();

    for (d=begin; d!=end; ++d)
    {
      int expNumVars = d->numVars;
      int expNumGlobal = d->numGlobals;
      Util::Expression* expPtr = d->expr;

      if (expNumVars > 0 || expPtr->isTimeDependent() || expNumGlobal > 0 )
      { 
        return false;
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs(const std::vector<int> & intLIDVecRef,
                                          const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];
  li_Bra = intLIDVec[0];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << section_divider << std::endl;

    Xyce::dout() << "::registerLIDs:\n";
    Xyce::dout() << "  name = " << getName() << std::endl;

    Xyce::dout() << "\nlocal solution indices:\n";
    Xyce::dout() << "  li_Pos = "<< li_Pos << std::endl;
    Xyce::dout() << "  li_Neg = "<< li_Neg << std::endl;
    Xyce::dout() << "  li_Bra = "<< li_Bra << std::endl;

    Xyce::dout() << section_divider << std::endl;
  }

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
// Creation Date : 6/22/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // copy over the global ID lists.
  staLIDVec = staLIDVecRef;

  li_fstate = staLIDVec[0];
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL
// Creation Date : 08/21/02
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp_BASE;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL
// Creation Date : 08/27/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  APosEquBraVarOffset = jacLIDVec[0][0];
  ANegEquBraVarOffset = jacLIDVec[1][0];
  ABraEquPosNodeOffset = jacLIDVec[2][0];
  ABraEquNegNodeOffset = jacLIDVec[2][1];
  ABraEquBraVarOffset = jacLIDVec[2][2];
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/26/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  double * qVec = extData.daeQVectorRawPtr;
  double * solVec = extData.nextSolVectorRawPtr;
  double * staVec = extData.nextStaVectorRawPtr;

  double current = solVec[li_Bra];
  if( (getSolverState().dcopFlag) && ICGiven )
    current = IC;

  f0 = L*current;
  staVec[li_fstate] = f0;

  qVec[li_Bra] += f0;

  return true;
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
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/26/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * solVec = extData.nextSolVectorRawPtr;
  double * fVec = extData.daeFVectorRawPtr;

  double vind = solVec[li_Pos]-solVec[li_Neg];

  // In the case that an initial condition is specified, the inductor is set up
  // like a current source for the DC operating point. We don't deal with the
  // node voltages in that case, so set coef to 0.
  double current = solVec[li_Bra];
  double coef = -vind;

  if (getSolverState().dcopFlag && ICGiven)
  {
    current = IC;
    coef = 0.0;
  }

  // load the current into the two KCL rhs vector elements
  fVec[li_Pos] += current;
  fVec[li_Neg] += -current;
  fVec[li_Bra] += coef;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 instance.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);
  dQdx[li_Bra][ABraEquBraVarOffset] += L;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  if ( getSolverState().dcopFlag && ICGiven )
  {
    // In the case that an initial condition is specified for an
    // inductor, the DC op should be set up like a current source just
    // for the operating point calculation.
    dFdx[li_Pos][APosEquBraVarOffset]  += 0.0;
    dFdx[li_Neg][ANegEquBraVarOffset]  += 0.0;
    dFdx[li_Bra][ABraEquPosNodeOffset] += 0.0;
    dFdx[li_Bra][ABraEquNegNodeOffset] += 0.0;
    dFdx[li_Bra][ABraEquBraVarOffset]  += 1.0;
  }
  else
  {
    dFdx[li_Pos][APosEquBraVarOffset]  += 1.0;
    dFdx[li_Neg][ANegEquBraVarOffset]  -= 1.0;
    dFdx[li_Bra][ABraEquPosNodeOffset] -= 1.0;
    dFdx[li_Bra][ABraEquNegNodeOffset] += 1.0;
    dFdx[li_Bra][ABraEquBraVarOffset]  += 0.0;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 1/11/02
//-----------------------------------------------------------------------------
bool Instance::setIC ()
{
  int i_bra_sol;
  int i_f_state;
  double * nextStaVector = extData.nextStaVectorRawPtr;
  double * currStaVector = extData.currStaVectorRawPtr;

  double * nextSolVector = extData.nextSolVectorRawPtr;
  double * currSolVector = extData.currSolVectorRawPtr;

  if (ICGiven)
  {
    f0 = L*IC;
    currStaVector[li_fstate] = f0;
    nextStaVector[li_fstate] = f0;

    currSolVector[li_Bra] = IC;
    nextSolVector[li_Bra] = IC;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/17/04
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  varTypeVec.resize(1);
  varTypeVec[0] = 'I';
}

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
// Purpose       : block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    inductanceMultiplier(1.0),
    IC(0.0),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    tnom(getDeviceOptions().tnom),
    tnomGiven(0)
{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
  if (!given("TNOM"))
    tnom = getDeviceOptions().tnom;

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
// Creation Date : 4/04/00
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i, isize;
  isize = instanceContainer.size();

  os << std::endl;
  os << "Number of Inductor instances: " << isize << std::endl;
  os << "    name=\t\tmodelName\tParameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << getName();
    os << "\t\tL = " << (*iter)->L;
    os << "\tIC = " << (*iter)->IC;
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


// Inductor Master functions:

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

  if (!separateInstances_ && ( loadType == LINEAR || loadType == NONLINEAR))
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
    Instance & inst = *(*it);

    double coef = 0.0;
    double current = solVec[inst.li_Bra];

    if( (getSolverState().dcopFlag) && inst.ICGiven )
      current = inst.IC;

    inst.f0 = inst.L*current;
    double vind = solVec[inst.li_Pos]-solVec[inst.li_Neg];

    // In the case that an initial condition is specified, the inductor is set up
    // like a current source for the DC operating point. We don't deal with the
    // node voltages in that case, so set coef to 0.
    if (getSolverState().dcopFlag && inst.ICGiven)
    {
      solVec[inst.li_Bra] = current;
    }
    else
    {
      coef = -vind;
    }

    // load the current into the two KCL rhs vector elements
    fVec[inst.li_Pos] += current;
    fVec[inst.li_Neg] += -current;
    fVec[inst.li_Bra] += coef;
    qVec[inst.li_Bra] += inst.f0;
    if( inst.loadLeadCurrent )
    {
      leadF[inst.li_branch_data] = current;
      junctionV[inst.li_branch_data] = vind;
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

  if (!separateInstances_ && ( loadType == LINEAR || loadType == NONLINEAR))
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
    Instance & inst = *(*it);

    if ( getSolverState().dcopFlag && inst.ICGiven )
    {
      // In the case that an initial condition is specified for an
      // inductor, the DC op should be set up like a current source just
      // for the operating point calculation.
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
      *inst.fBraEquBraVarPtr  += 1.0;
#else
      dFdx[inst.li_Bra][inst.ABraEquBraVarOffset]  += 1.0;
#endif
    }
    else
    {
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
      *inst.fPosEquBraVarPtr  += 1.0;
      *inst.fNegEquBraVarPtr  -= 1.0;
      *inst.fBraEquPosNodePtr -= 1.0;
      *inst.fBraEquNegNodePtr += 1.0;
#else
      dFdx[inst.li_Pos][inst.APosEquBraVarOffset]  += 1.0;
      dFdx[inst.li_Neg][inst.ANegEquBraVarOffset]  -= 1.0;
      dFdx[inst.li_Bra][inst.ABraEquPosNodeOffset] -= 1.0;
      dFdx[inst.li_Bra][inst.ABraEquNegNodeOffset] += 1.0;
#endif
    }

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    *inst.qBraEquBraVarPtr += inst.L;
#else
    dQdx[inst.li_Bra][inst.ABraEquBraVarOffset] += inst.L;
#endif
  }
  return true;
}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void registerDevice() 
{
  Config<Traits>::addConfiguration()
    .registerDevice("l", 1)
    .registerModelType("l", 1)
    .registerModelType("ind", 1);
}

//-----------------------------------------------------------------------------
// Function      : indSensitivity::operator
// Purpose       : produces df/dp and dq/dp, where p=L.  
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 7/31/2014
//-----------------------------------------------------------------------------
void indSensitivity::operator()(
    const ParameterBase &entity,
    const std::string & name,
    std::vector<double> & dfdp, 
    std::vector<double> & dqdp, 
    std::vector<double> & dbdp, 
    std::vector<int> & Findices,
    std::vector<int> & Qindices,
    std::vector<int> & Bindices
    ) const
{
  const ParameterBase * e1 = &entity;
  const Instance * in = dynamic_cast<const Instance *> (e1);

  double * solVec = in->extData.nextSolVectorRawPtr;
  double current = solVec[in->li_Bra];
  if( (in->getSolverState().dcopFlag) && in->ICGiven )
  {
    current = in->IC;
  }

  double dqdpLoc = current;

  dqdp.resize(1);
  dqdp[0] = dqdpLoc;

  Qindices.resize(1);
  Qindices[0] = in->li_Bra;
}

} // namespace Inductor
} // namespace Device
} // namespace Xyce

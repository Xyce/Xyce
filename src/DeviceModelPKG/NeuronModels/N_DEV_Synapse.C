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

//----------------------------------------------------------------------------
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
//----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
//
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Synapse.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {

namespace Synapse {


void Traits::loadInstanceParameters(ParametricData<Synapse::Instance> &p)
{
  // Set up configuration constants:
}

void Traits::loadModelParameters(ParametricData<Synapse::Model> &p)
{
  // Default values are taken from Destexhe 98, values for AMPA receptor
  p.addPar ("GMAX", 1.0E-9, &Synapse::Model::gMax)
    .setUnit(U_OHMM1)
    .setDescription("Maximal Synaptic Conductance");

  p.addPar ("EREV", 0.0, &Synapse::Model::eRev)
    .setUnit(U_VOLT)
    .setDescription("Reversal Potential");

  p.addPar ("ALPHA", 1.1E-6, &Synapse::Model::alpha)
    .setUnit(U_SECM1)
    .setDescription("Forward rate constant for receptor opening");

  p.addPar ("BETA", 190.0, &Synapse::Model::beta)
    .setUnit(U_SECM1)
    .setDescription("Backward rate constant for receptor opening");

  p.addPar ("VP", 0.002, &Synapse::Model::vP)
    .setUnit(U_VOLT)
    .setDescription("Presynaptic voltage at which neurotransmitter concentration is half-maximal");

  // KP should NOT be 0
  p.addPar ("KP", 0.005, &Synapse::Model::kP)
    .setUnit(U_VOLT)
    .setDescription("Steepness parameter for neurotransmitter concentration");

  p.addPar ("TMAX", 0.001, &Synapse::Model::tMax)
    .setDescription("Maximal neurotransmitter concentration");

  // NOTE:  not including concentration units - TMAX should be in moles/liter, and
  // alpha has concentration in the denominator; they cancel out
}



std::vector< std::vector<int> > Instance::jacStamp;

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  // now set the temperature related stuff.
  //updateTemperature(temp);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & IB,
  Model & Riter,
  const FactoryBlock &  factory_block)

  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Riter),
    li_Prev(-1),
    li_Post(-1),
    APostEquPostNodeOffset(-1),
    APostEquRNodeOffset(-1),
    AREquPostNodeOffset(-1),
    AREquRNodeOffset(-1),
    f_PostEquPostNodePtr(0),
    f_PostEquRNodePtr(0),
    f_REquPostNodePtr(0),
    f_REquRNodePtr(0)
{
  numIntVars   = 1;
  numExtVars   = 2;
  numStateVars = 0;

  if( jacStamp.empty() )
  {
    // Vpre, Vpost, r
    jacStamp.resize(3);
    jacStamp[0].resize(0);
    jacStamp[1].resize(2);
    jacStamp[2].resize(2);
    jacStamp[1][0] = 1;  // Vpost  PostPostVar
    jacStamp[1][1] = 2;  // r      PostRVar
    jacStamp[2][0] = 0;  // Vpre   RVpre
    jacStamp[2][1] = 2;  // r      RR
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:

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
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/12/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef )
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  SynapseInstance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  li_Prev = extLIDVec[0];
  li_Post = extLIDVec[1];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << "  li_Prev = " << li_Prev << std::endl;
    Xyce::dout() << "  li_Post = " << li_Post << std::endl;
  }

  li_rVar = intLIDVec[0];

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
  addInternalNode(symbol_table, li_rVar, getName(), "r");
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       : Note that the Synapse does not have any state vars.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/12/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs(const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/20/01
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
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/27/01
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  APostEquPostNodeOffset = jacLIDVec[1][0];
  APostEquRNodeOffset = jacLIDVec[1][1];
  AREquPostNodeOffset = jacLIDVec[2][0];
  AREquRNodeOffset = jacLIDVec[2][1];
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
  f_PostEquPostNodePtr = &(dFdx[li_Post][APostEquPostNodeOffset]);
  f_PostEquRNodePtr = &(dFdx[li_Post][APostEquRNodeOffset]);
  f_REquPostNodePtr = &(dFdx[li_rVar][AREquPostNodeOffset]);
  f_REquRNodePtr = &(dFdx[li_rVar][AREquRNodeOffset]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one diode instance
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, Dept. 9233.
// Creation Date : 3/05/04
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars()
{
  bool bsuccess = true;

  Linear::Vector *  solVecPtr = extData.nextSolVectorPtr;

  double vPre  = (*solVecPtr)[li_Prev];
  double vPost  = (*solVecPtr)[li_Post];
  double rVar  = (*solVecPtr)[li_rVar];

  // Load RHS vector element for the positive circuit node KCL equ.
  {
    Sacado::Fad::SFad<double,2> vPostVar( 2, 0, vPost );
    Sacado::Fad::SFad<double,2> rVarS( 2, 1, rVar);

    // parameters
    Sacado::Fad::SFad<double,2> gMaxVar( model_.gMax );
    Sacado::Fad::SFad<double,2> eRevVar( model_.eRev );

    Sacado::Fad::SFad<double,2> resultFad;
    resultFad = PostCurrentEqu( vPostVar, rVarS, gMaxVar, eRevVar );
    ipost = resultFad.val();
    didVpost = resultFad.dx(0);
    didr = resultFad.dx(1);
  }
  {
    Sacado::Fad::SFad<double,2> vPreVar( 2, 0, vPre );
    Sacado::Fad::SFad<double,2> rVarS( 2, 1, rVar);

    // parameters
    Sacado::Fad::SFad<double,2> alphaVar( model_.alpha );
    Sacado::Fad::SFad<double,2> betaVar( model_.beta );
    Sacado::Fad::SFad<double,2> tMaxVar( model_.tMax );
    Sacado::Fad::SFad<double,2> vPVar( model_.vP );
    Sacado::Fad::SFad<double,2> kPVar( model_.kP );

    Sacado::Fad::SFad<double,2> resultFad;
    resultFad = rEquF( vPreVar, rVarS, alphaVar, betaVar, tMaxVar, vPVar, kPVar);
    rFval = resultFad.val();
    drFdVpre = resultFad.dx(0);
    drFdr = resultFad.dx(1);
  }
  return bsuccess;
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
  bool bsuccess = updateIntermediateVars();
  return  bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/29/01
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState()
{
  return  true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 synapse
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 8/31/2010
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  Linear::Vector *  qVecPtr = extData.daeQVectorPtr;
  Linear::Vector *  solVecPtr = extData.nextSolVectorPtr;
  double rVar  = (*solVecPtr)[li_rVar];
  (*qVecPtr)[li_rVar] -= rVar;

  return true;
}
//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 Synapse  instance.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/24/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  Linear::Vector *  fVecPtr = extData.daeFVectorPtr;
  (*fVecPtr)[li_Prev] += 0.0;
  (*fVecPtr)[li_Post] += ipost;
  (*fVecPtr)[li_rVar] += rFval;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the dQdx-matrix contributions for a single
//                 synapse instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);
  dQdx[li_rVar][AREquRNodeOffset] += -1.0;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
// Purpose       : Loads the F-vector contributions for a single
//                 Synapse  instance.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  dFdx[li_Post][APostEquPostNodeOffset] += didVpost;
  dFdx[li_Post][APostEquRNodeOffset]    += didr;
  dFdx[li_rVar][AREquPostNodeOffset]    += drFdVpre;
  dFdx[li_rVar][AREquRNodeOffset]       += drFdr;

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
bool Instance::updateTemperature ( const double & temp_tmp)
{
  bool bsuccess = true;
  return bsuccess;
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
  if (kP == 0.0)
    kP = 0.001;

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
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    gMax(0.0),
    eRev(0.0),
    alpha(0.0),
    beta(0.0),
    vP(0.0),
    kP(0.0),
    tMax(0.0)
{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams();
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

  int i,isize;
  isize = instanceContainer.size();
  os << std::endl;
  os << "Number of Synapse Instances: " << isize << std::endl;
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


// Synapse Master functions:

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
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    (*it)->updateIntermediateVars();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateSecondaryState ( double * staDerivVec, double * stoVec )
{
  return true;
}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  static bool initialized = false;

  if (!initialized && (deviceMap.empty() ||
      ((deviceMap.find("SYNAPSE")!=deviceMap.end()) && (levelSet.find(1)!=levelSet.end()))))
  {
    initialized = true;

    Config<Traits>::addConfiguration()
      .registerDevice("synapse", 1)
      .registerModelType("synapse", 1);
  }
}

} // namespace Synapse
} // namespace Device
} // namespace Xyce

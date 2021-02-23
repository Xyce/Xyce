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
// Purpose        :  Provide a linking device for coupling to Alegra
//                   (and other) simulations.
//
// Special Notes  :  The external simulator is responsible for passing
//                   in the number of internal variables in an initialization
//                   step, and for providing F, Q, B (and their derivatives)
//                   each time Xyce is asked to simulate.
//
// Creator        : Tom Russo, SNL, Electrical Models and Simulation
//
// Creation Date  : 02/28/2017
//
//
//
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------

#include <sstream>

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_GeneralExternal.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_AssemblyTypes.h>

#include <N_DEV_VectorComputeInterface.h>

namespace Xyce {
namespace Device {

// These manage the data for the DPARAMS vector-composite parameter
template<>
ParametricData<GenExtDoubleData>::ParametricData()
{
  addPar ( "NAME", "PARAM0", &GenExtDoubleData::name_);
  addPar ( "VALUE", 0.0, &GenExtDoubleData::value_);
}

ParametricData<GenExtDoubleData> &GenExtDoubleData::getParametricData()
{
  static ParametricData<GenExtDoubleData> parMap;
  return parMap;
}

//-----------------------------------------------------------------------------
// Function      : GenExtDoubleData::GenExtDoubleData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
GenExtDoubleData::GenExtDoubleData()
  : CompositeParam(getParametricData()),
    name_("param"),
    value_(0.0)
{}

//-----------------------------------------------------------------------------
// Function      : GenExtDoubleData::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
void GenExtDoubleData::processParams ()
{
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : "<<" operator for GenExtDoubleData
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const GenExtDoubleData & gEdd)
{
  os << " GenExtDoubleData for: name = " << gEdd.getName() <<
    " Value=" << gEdd.getValue() <<
    std::endl;

  return os;
}


// These manage the data for the IPARAMS vector-composite parameter
template<>
ParametricData<GenExtIntData>::ParametricData()
{
  addPar ( "NAME", "PARAM0", &GenExtIntData::name_);
  addPar ( "VALUE", 0, &GenExtIntData::value_);
}

ParametricData<GenExtIntData> &GenExtIntData::getParametricData()
{
  static ParametricData<GenExtIntData> parMap;
  return parMap;
}

//-----------------------------------------------------------------------------
// Function      : GenExtIntData::GenExtIntData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
GenExtIntData::GenExtIntData()
  : CompositeParam(getParametricData()),
    name_("param"),
    value_(0.0)
{}

//-----------------------------------------------------------------------------
// Function      : GenExtIntData::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
void GenExtIntData::processParams ()
{
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : "<<" operator for GenExtIntData
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const GenExtIntData & gEid)
{
  os << " GenExtIntData for: name = " << gEid.getName() <<
    " Value=" << gEid.getValue() <<
    std::endl;

  return os;
}

// These manage the data for the BPARAMS vector-composite parameter
template<>
ParametricData<GenExtBoolData>::ParametricData()
{
  addPar ( "NAME", "PARAM0", &GenExtBoolData::name_);
  addPar ( "VALUE", false, &GenExtBoolData::value_);
}

ParametricData<GenExtBoolData> &GenExtBoolData::getParametricData()
{
  static ParametricData<GenExtBoolData> parMap;
  return parMap;
}

//-----------------------------------------------------------------------------
// Function      : GenExtBoolData::GenExtBoolData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
GenExtBoolData::GenExtBoolData()
  : CompositeParam(getParametricData()),
    name_("param"),
    value_(false)
{}

//-----------------------------------------------------------------------------
// Function      : GenExtBoolData::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
void GenExtBoolData::processParams ()
{
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : "<<" operator for GenExtBoolData
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const GenExtBoolData & gEbd)
{
  os << " GenExtBoolData for: name = " << gEbd.getName() <<
    " Value=" << gEbd.getValue() <<
    std::endl;

  return os;
}

// These manage the data for the SPARAMS vector-composite parameter
template<>
ParametricData<GenExtStringData>::ParametricData()
{
  addPar ( "NAME", "PARAM0", &GenExtStringData::name_);
  addPar ( "VALUE", "", &GenExtStringData::value_);
}

ParametricData<GenExtStringData> &GenExtStringData::getParametricData()
{
  static ParametricData<GenExtStringData> parMap;
  return parMap;
}

//-----------------------------------------------------------------------------
// Function      : GenExtStringData::GenExtStringData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
GenExtStringData::GenExtStringData()
  : CompositeParam(getParametricData()),
    name_("param"),
    value_("")
{}

//-----------------------------------------------------------------------------
// Function      : GenExtStringData::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
void GenExtStringData::processParams ()
{
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : "<<" operator for GenExtStringData
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const GenExtStringData & gEsd)
{
  os << " GenExtStringData for: name = " << gEsd.getName() <<
    " Value=" << gEsd.getValue() <<
    std::endl;

  return os;
}

namespace GeneralExternal {

void Traits::loadInstanceParameters(ParametricData<GeneralExternal::Instance> &p)
{
  // Vector-composite for double params
  p.addComposite("DPARAMS", GenExtDoubleData::getParametricData(), &GeneralExternal::Instance::doubleDataMap_);
  // Vector-composite for int params
  p.addComposite("IPARAMS", GenExtIntData::getParametricData(), &GeneralExternal::Instance::intDataMap_);
  p.addComposite("BPARAMS", GenExtBoolData::getParametricData(), &GeneralExternal::Instance::boolDataMap_);
  p.addComposite("SPARAMS", GenExtStringData::getParametricData(), &GeneralExternal::Instance::stringDataMap_);
}

void Traits::loadModelParameters(ParametricData<GeneralExternal::Model> &p)
{}


// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
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
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
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
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       Miter,
  const FactoryBlock &          factory_block)

  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    vciPtr_(0)
{
  numExtVars   = IB.numExtVars;  // we have as many as were specified on the
                                 // instance line
  numIntVars=0;
  numStateVars=0;

  // Do not allocate lead current vectors by default:
  setNumBranchDataVars(0);
  // but if requested, allocate one for each external node
  numBranchDataVarsIfAllocated=numExtVars;

  leadCurrentF.resize(numBranchDataVarsIfAllocated);
  leadCurrentQ.resize(numBranchDataVarsIfAllocated);

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

  // Now do some stuff with vector-composite params
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    if (!doubleDataVec_.empty())
    {
      Xyce::dout() << getName() << " was given " << doubleDataVec_.size()
                   << " double params." << std::endl;
      for (int i=0; i< doubleDataVec_.size(); i++)
      {
        Xyce::dout() << "Parameter " << doubleDataVec_[i]->getName()
                     << " has value " << doubleDataVec_[i]->getValue()
                     << std::endl;
      }
    }
    else
    {
      Xyce::dout() << getName() << " was not given any double params."
                   << std::endl;
    }
    if (!intDataVec_.empty())
    {
      Xyce::dout() << getName() << " was given " << intDataVec_.size()
                   << " int params." << std::endl;
      for (int i=0; i< intDataVec_.size(); i++)
      {
        Xyce::dout() << "Parameter " << intDataVec_[i]->getName()
                     << " has value " << intDataVec_[i]->getValue()
                     << std::endl;
      }
    }
    else
    {
      Xyce::dout() << getName() << " was not given any int params."
                   << std::endl;
    }

    if (!boolDataVec_.empty())
    {
      Xyce::dout() << getName() << " was given " << boolDataVec_.size()
                   << " bool params." << std::endl;
      for (int i=0; i< boolDataVec_.size(); i++)
      {
        Xyce::dout() << "Parameter " << boolDataVec_[i]->getName()
                     << " has value " << boolDataVec_[i]->getValue()
                     << std::endl;
      }
    }
    else
    {
      Xyce::dout() << getName() << " was not given any bool params."
                   << std::endl;
    }

    if (!stringDataVec_.empty())
    {
      Xyce::dout() << getName() << " was given " << stringDataVec_.size()
                   << " string params." << std::endl;
      for (int i=0; i< stringDataVec_.size(); i++)
      {
        Xyce::dout() << "Parameter " << stringDataVec_[i]->getName()
                     << " has value " << stringDataVec_[i]->getValue()
                     << std::endl;
      }
    }
    else
    {
      Xyce::dout() << getName() << " was not given any string params."
                   << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
Instance::~Instance()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupDenseJacStamp_
// Purpose       : Utility function to set up the jacobian stamp
// Special Notes : Just makes a dense jacobian stamp for now
// Scope         : private
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
void Instance::setupDenseJacStamp_()
{
  int numNodes = numExtVars+numIntVars;

  jacStamp_.resize(numNodes);
  for (int i=0; i<numNodes; ++i)
  {
    jacStamp_[i].resize(numNodes);
    for (int j=0; j<numNodes; ++j)
    {
      jacStamp_[i][j]=j;
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
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

  int numNodes = numExtVars + numIntVars;

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  li_Nodes_.resize(numNodes);

  int extVar=0;
  int intVar=0;
  int theVar=0;

  // First comes the external nodes
  for (theVar=0; theVar<numExtVars; ++theVar)
  {
    li_Nodes_[theVar] = extLIDVec[theVar];
  }
  // now the internals

  for (theVar=0; theVar<+numIntVars; ++theVar)
  {
      li_Nodes_[theVar+numExtVars] = intLIDVec[theVar];
  }


  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    for (int i=0; i<numIntVars+numExtVars;++i)
    {
      Xyce::dout() << "LID var " << i << " ";
      if (i>=numExtVars)
        Xyce::dout() << " (internal) ";
      Xyce::dout() << ": "<< li_Nodes_[i] << std::endl;
    }
    Xyce::dout() << section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerBranchDataLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL
// Creation Date : 11/15/2017
//-----------------------------------------------------------------------------
void Instance::registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef)
{
  AssertLIDs(branchLIDVecRef.size() == getNumBranchDataVars());

  if (loadLeadCurrent)
  {
    li_Branches_.resize(numBranchDataVarsIfAllocated);
    for (int i=0;i<numBranchDataVarsIfAllocated;i++)
    {
      li_Branches_[i] = branchLIDVecRef[i];
    }
  }
}
//-----------------------------------------------------------------------------
// Function      : Instance::loadNodeSymbols
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
void Instance::loadNodeSymbols(Util::SymbolTable &symbol_table) const
{
  if (numIntVars != 0)
  {
    for (int i = 0; i < numIntVars; ++i)
    {
      std::ostringstream ost;
      ost << "InternalNode_" << i;
      addInternalNode(symbol_table, li_Nodes_[numExtVars+i],
                      getName(), ost.str());
    }
  }
  if (loadLeadCurrent)
  {
    for (int i=0; i< numBranchDataVarsIfAllocated; i++)
    {
      std::ostringstream ost;
      ost << "BRANCH_D" << i+1;
      addBranchDataNode(symbol_table, li_Branches_[i], getName(), ost.str());
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
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
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp_;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  if (jacStamp_.empty())
  {
    DevelFatal(*this) << "Improper setup.  No Jacobian stamp has been set."
                      << std::endl
                      << "The General External device can only be used through the "
                      << "GenCouplingSimulator interface, and cannot be used in normal Xyce runs. "
                      << std::endl
                      << "If you are calling this device through the GenCouplingInterface, be sure "
                      << "that you are setting the number of interal variables with setNumInternalVariables"
                      << "and associating  a vector loader with setVectorLoader." << std::endl;
  }

  DeviceInstance::registerJacLIDs( jacLIDVec );

  int numVars=numExtVars+numIntVars;
   A_Equ_NodeOffsets_.resize(numVars);
  for (int equ=0; equ < numVars; ++equ)
  {
    A_Equ_NodeOffsets_[equ].resize(jacStamp_[equ].size());
    for (int node=0; node < jacStamp_[equ].size(); ++node)
    {
      A_Equ_NodeOffsets_[equ][node] = jacLIDVec[equ][node];
    }
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState()
{
  bool bsuccess = true;

  bsuccess = updateIntermediateVars();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  bool bsuccess = true;

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{

  bool bsuccess=true;
  int numVars=numIntVars+numExtVars;
  double * solVector = extData.nextSolVectorRawPtr;
  if (vciPtr_)
  {
    if (solutionVars_.empty())
      solutionVars_.resize(numVars);

    // Copy out the solution variables we need.  The ordering
    // is already correct, because of the way we constructed
    // li_Nodes
    for (int i=0;i<numVars;i++)
      solutionVars_[i] = solVector[li_Nodes_[i]];

    // Call back to the computation object and get this devices
    // contributions.  IT IS THAT FUNCTION'S RESPONSIBILITY TO SIZE THESE
    // VECTORS APPROPRIATELY!
    bsuccess=vciPtr_->computeXyceVectors(solutionVars_,
                                         getSolverState().currTime_,
                                         FVec_, QVec_, BVec_,
                                         dFdXMat_, dQdXMat_);

    if (!bsuccess)
    {
      UserError(*this) << "Time Domain vector compute function returned false!";
    }

    if (loadLeadCurrent)
    {
      for (int i=0; i< numBranchDataVarsIfAllocated; i++)
      {
        leadCurrentF[i] = leadCurrentQ[i] = 0.0;
      }

      if (!FVec_.empty())
      {
        for (int i=0;i<numExtVars;i++)
        {
          leadCurrentF[i] += FVec_[i];
        }
      }
      if (!QVec_.empty())
      {
        for (int i=0;i<numExtVars;i++)
        {
          leadCurrentQ[i] += QVec_[i];
        }
      }
      if (!BVec_.empty())
      {
        for (int i=0;i<numExtVars;i++)
        {
          leadCurrentF[i] -= BVec_[i];
        }
      }
    }

    // NOTE:  This sanity checking could be costly, and should
    // only be enabled with DEBUG_DEVICE turned on.
    if (DEBUG_DEVICE)
    {
      if (bsuccess)
      {
        // Sanity check sizes.  Each vector MUST be either 0 length (indicating
        // no contribution to that vector) *OR* numVars long.  And at least ONE
        // of the vectors MUST be non-zero length.
        if ((FVec_.empty() && QVec_.empty() && BVec_.empty())
            || (FVec_.size()!=numVars && QVec_.size()!=numVars
                && BVec_.size()!= numVars)
            )
        {
          UserError(*this) << "At least one of the F, Q, or B matrices must have length " << numVars;
          bsuccess=false;
        }
        else
        {
          if ( !((FVec_.size() == numVars || FVec_.size() == 0) &&
                 (QVec_.size() == numVars || QVec_.size() == 0) &&
                 (BVec_.size() == numVars || BVec_.size() == 0)))
          {
            UserError(*this) << "F, Q, or B matrices must have length "
                             << numVars << " or zero. "
                             << " F length=" << FVec_.size()
                             << " Q length=" << QVec_.size()
                             << " B length=" << BVec_.size();
            bsuccess=false;
          }
          else
          {
            if ( dFdXMat_.size() != FVec_.size() || dQdXMat_.size() != QVec_.size())
            {
              UserError(*this) << "dFdX matrix has " << dFdXMat_.size()
                               << " rows "
                               << " and should match F vector length "
                               << FVec_.size()
                               << ", dQdX matrix has " << dQdXMat_.size()
                               << " rows "
                               << " and should match Q vector length "
                               << QVec_.size();
              bsuccess=false;
            }
            else
            {
              // Here's the potentially costly one
              // All the vectors are the right sizes, all the matrices
              // have the right rows --- do all the rows have the right
              // lengths?
              for (int i=0; i<dFdXMat_.size();i++)
              {
                if (dFdXMat_[i].size() != FVec_.size())
                {
                  UserError(*this) << " F matrix row " << i << " has "
                                   << dFdXMat_[i].size()
                                   << " must have " << FVec_.size();
                  bsuccess=false;
                }
              }
              for (int i=0; i<dQdXMat_.size();i++)
              {
                if (dQdXMat_[i].size() != QVec_.size())
                {
                  UserError(*this) << " Q matrix row " << i << " has "
                                   << dQdXMat_[i].size()
                                   << " must have " << QVec_.size();
                  bsuccess=false;
                }
              }
            }
          }
        }
      }
    }
  }
  else
  {
    bsuccess=false;
    DevelFatal(*this) << "Improper setup.  No vector loader for this has been assigned."
                      << std::endl
                      << "The General External device can only be used through the"
                      << "GenCouplingSimulator interface, and cannot be used in normal Xyce runs."
                      << std::endl
                      << "If you are calling this device through the GenCouplingInterface, be sure"
                      << "that you are setting the number of interal variables with setNumInternalVariables "
                      << "and associating  a vector loader with setVectorLoader." << std::endl;

  }

  return bsuccess;
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
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  bool bsuccess = true;
  int numVars = numExtVars+numIntVars;

  Linear::Vector * daeQVecPtr = extData.daeQVectorPtr;

  if (QVec_.size() != 0)
  {
    for (int i=0; i<numVars; i++)
      (*daeQVecPtr)[li_Nodes_[i]] += QVec_[i];
  }

  if (loadLeadCurrent)
  {
    double * leadQ = extData.nextLeadCurrQCompRawPtr;
    for (int i=0; i<numBranchDataVarsIfAllocated; i++)
    {
      leadQ[li_Branches_[i]] = leadCurrentQ[i];
    }
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
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess=true;

  Linear::Vector * daeFVecPtr = extData.daeFVectorPtr;
  int numVars = numExtVars+numIntVars;

  if (FVec_.size() != 0)
  {
    for (int i=0; i<numVars; i++)
      (*daeFVecPtr)[li_Nodes_[i]] += FVec_[i];
  }

  if (loadLeadCurrent)
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    for (int i=0; i<numBranchDataVarsIfAllocated; i++)
    {
      leadF[li_Branches_[i]] = leadCurrentF[i];
    }
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEBVector
//
// Purpose       : Loads the B-vector contributions for a single
//                 instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
bool Instance::loadDAEBVector ()
{
  bool bsuccess=true;

  Linear::Vector * daeBVecPtr = extData.daeBVectorPtr;
  int numVars = numExtVars+numIntVars;

  if (BVec_.size() != 0)
  {
    for (int i=0; i<numVars; i++)
      (*daeBVecPtr)[li_Nodes_[i]] += BVec_[i];
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 instance.
//
// Special Notes:  We are assuming a dense jacstamp here!
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;
  Linear::Matrix * dQdxMatPtr = extData.dQdxMatrixPtr;
  int numVars = numExtVars+numIntVars;

  if (QVec_.size() != 0)
  {
    for (int row=0; row<numVars; row++)
    {
      for (int col=0; col<A_Equ_NodeOffsets_[row].size(); col++)
      {
        (*dQdxMatPtr)[li_Nodes_[row]][A_Equ_NodeOffsets_[row][col]]
          += dQdXMat_[row][jacStamp_[row][col]];
      }
    }
  }

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 instance.
//
// Special Notes:  We are assuming a dense jacstamp here!
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;

  Linear::Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;
  int numVars = numExtVars+numIntVars;

  if (FVec_.size() != 0)
  {
    for (int row=0; row<numVars; row++)
    {
      for (int col=0; col<A_Equ_NodeOffsets_[row].size(); col++)
      {
        (*dFdxMatPtr)[li_Nodes_[row]][A_Equ_NodeOffsets_[row][col]]
          += dFdXMat_[row][jacStamp_[row][col]];
      }
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getSolution
//
// Purpose       : Copies this device's solution variables into provide
//                 vector
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
void Instance::getSolution(std::vector<double> &sV)
{
  double * solVector = extData.nextSolVectorRawPtr;
  int numVars = numExtVars+numIntVars;

  if (solutionVars_.empty())
    solutionVars_.resize(numVars);

  // Copy out the solution variables we need.  The ordering
  // is already correct, because of the way we constructed
  // li_Nodes
  for (int i=0;i<numVars;i++)
    solutionVars_[i] = solVector[li_Nodes_[i]];

  sV = solutionVars_;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getDParams
//
// Purpose       : Copies name/value pairs of params from DPARAMS into provided
//                 vectors
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
void Instance::getDParams(std::vector<std::string> &names,
                          std::vector<double> &values)
{
  names.clear();
  values.clear();
  for (int i=0; i< doubleDataVec_.size(); i++)
  {
    names.push_back(doubleDataVec_[i]->getName());
    values.push_back(doubleDataVec_[i]->getValue());
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::getIParams
//
// Purpose       : Copies name/value pairs of params from IPARAMS into provided
//                 vectors
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
void Instance::getIParams(std::vector<std::string> &names,
                          std::vector<int> &values)
{
  names.clear();
  values.clear();
  for (int i=0; i< intDataVec_.size(); i++)
  {
    names.push_back(intDataVec_[i]->getName());
    values.push_back(intDataVec_[i]->getValue());
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::getBParams
//
// Purpose       : Copies name/value pairs of params from BPARAMS into provided
//                 vectors
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
void Instance::getBParams(std::vector<std::string> &names,
                          std::vector<bool> &values)
{
  names.clear();
  values.clear();
  for (int i=0; i< boolDataVec_.size(); i++)
  {
    names.push_back(boolDataVec_[i]->getName());
    values.push_back(boolDataVec_[i]->getValue());
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::getSParams
//
// Purpose       : Copies name/value pairs of params from SPARAMS into provided
//                 vectors
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
void Instance::getSParams(std::vector<std::string> &names,
                          std::vector<std::string> &values)
{
  names.clear();
  values.clear();
  for (int i=0; i< stringDataVec_.size(); i++)
  {
    names.push_back(stringDataVec_[i]->getName());
    values.push_back(stringDataVec_[i]->getValue());
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::constructComposite
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
CompositeParam *Instance::constructComposite(const std::string & cName, const std::string & pName)
{
  if (cName == "DPARAMS")
  {
    GenExtDoubleData *gEdd = new GenExtDoubleData ();
    doubleDataVec_.push_back(gEdd);
    return (static_cast<CompositeParam *> (gEdd));
  }
  else if (cName == "IPARAMS")
  {
    GenExtIntData *gEid = new GenExtIntData ();
    intDataVec_.push_back(gEid);
    return (static_cast<CompositeParam *> (gEid));
  }
  else if (cName == "BPARAMS")
  {
    GenExtBoolData *gEbd = new GenExtBoolData ();
    boolDataVec_.push_back(gEbd);
    return (static_cast<CompositeParam *> (gEbd));
  }
  else if (cName == "SPARAMS")
  {
    GenExtStringData *gEsd = new GenExtStringData ();
    stringDataVec_.push_back(gEsd);
    return (static_cast<CompositeParam *> (gEsd));
  }
  else
  {
    DevelFatal(*this).in("Instance::constructComposite")
      << "unrecognized composite name: "
      << cName;
  }
  // never reached
  return NULL;
}

// Master class functions

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       : Perform top-level loop over all genext devices,
//                 calling their updatePrimaryState() functions as needed
// Special Notes : A general external device must implement the time-domain
//                 vector loader, but may or may not implement a frequency
//                 domain loader.  If we're doing FD loads, we MUST use
//                 any FD loaders provided, and refrain from loading TD.
//                 If only a TD loader is provided, we use it in FD.
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 22 March 2018
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec,
                            int loadType)
{
  bool bsuccess = true;

  for (InstanceVector::const_iterator it = getInstanceBegin();
       it != getInstanceEnd();
       ++it)
  {
    // Here we are explicitly treating all general devices as
    // nonlinear even if the details of the vector loader are implemented
    // so that they could have been treated as linear.
    if ((( loadType == NONLINEAR )
         || (( loadType == NONLINEAR_FREQ ) &&
             !((*it)->vciPtr_->haveFDLoads()) )
         || ( loadType == ALL ) ))
    {
      bsuccess = (*it)->updatePrimaryState() && bsuccess;
    }
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       : Perform top-level loop over all genext devices,
//                 calling their loadDAE*Vector() functions as needed
// Special Notes : A general external device must implement the time-domain
//                 vector loader, but may or may not implement a frequency
//                 domain loader.  If we're doing FD loads, we MUST use
//                 any FD loaders provided, and refrain from loading TD.
//                 If only a TD loader is provided, we use it in FD.
//
//                 Note that the various vectors passed in are never actually
//                 used elsewhere --- that's the way it is in the
//                 default DeviceMaster::loadDAEVectors function, too.
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 22 March 2018
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double * qVec,
                               double * bVec, double * leadF, double * leadQ,
                               double * junctionV, int loadType)
{
  bool bsuccess = true;
  for (InstanceVector::const_iterator it = getInstanceBegin();
       it != getInstanceEnd();
       ++it)
  {
    // Here we are explicitly treating all general devices as
    // nonlinear even if the details of the vector loader are implemented
    // so that they could have been treated as linear.
    if ((( loadType == NONLINEAR )
         || (( loadType == NONLINEAR_FREQ ) &&
             !((*it)->vciPtr_->haveFDLoads()) )
         || ( loadType == ALL ) ))
    {
      bool tmpBool = (*it)->loadDAEFVector();
      bsuccess = bsuccess && tmpBool;

      tmpBool = (*it)->loadDAEQVector();
      bsuccess = bsuccess && tmpBool;

      tmpBool = (*it)->loadDAEBVector();
      bsuccess = bsuccess && tmpBool;
    }
  }
  return bsuccess;
}
//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       : Perform top-level loop over all genext devices,
//                 calling their loadDAE*Vector() functions as needed
// Special Notes : A general external device must implement the time-domain
//                 vector loader, but may or may not implement a frequency
//                 domain loader.  If we're doing FD loads, we MUST use
//                 any FD loaders provided, and refrain from loading TD.
//                 If only a TD loader is provided, we use it in FD.
//
//                 Note that the various vectors passed in are never actually
//                 used elsewhere --- that's the way it is in the
//                 default DeviceMaster::loadDAEMatrices function, too.
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 22 March 2018
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdX,
                              Linear::Matrix & dQdX,
                              int loadType)
{
  bool bsuccess = true;
  for (InstanceVector::const_iterator it = getInstanceBegin();
       it != getInstanceEnd();
       ++it)
  {
    // Here we are explicitly treating all general devices as
    // nonlinear even if the details of the vector loader are implemented
    // so that they could have been treated as linear.
    if ((( loadType == NONLINEAR )
         || (( loadType == NONLINEAR_FREQ ) &&
             !((*it)->vciPtr_->haveFDLoads()) )
         || ( loadType == ALL ) ))
    {
      bool tmpBool = (*it)->loadDAEdFdx();
      bsuccess = bsuccess && tmpBool;

      tmpBool = (*it)->loadDAEdQdx();
      bsuccess = bsuccess && tmpBool;
    }
  }
  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Master::updateFDIntermediateVars
// Purpose       : Calls vector loader and sets the frequency domain
//                 intermediate variables.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo and Ting Mei
// Creation Date : 20 March 2018
//-----------------------------------------------------------------------------

bool Master::updateFDIntermediateVars(double frequency,
                                      std::complex<double> * solVec)
{
  bool bsuccess=true;
  InstanceVector::const_iterator it, end;

  end = getInstanceEnd();

  for (  it = getInstanceBegin(); it != end; ++it)
  {
    Instance & theInstance = *(*it);
    int numVars=theInstance.numIntVars+theInstance.numExtVars;

    theInstance.solutionFDVars_.resize(numVars);
    for (size_t i = 0; i< numVars; i++)
    {
      theInstance.solutionFDVars_[i] = ( solVec[theInstance.li_Nodes_[i]] );
    }

    if (theInstance.vciPtr_)
    {
      // This conditional not really necessary as the default implemenation
      // of computeXyceFDVectors is actually an empty function, but
      // we'll check anyway
      if (theInstance.vciPtr_->haveFDLoads())
      {
        bsuccess =
          theInstance.vciPtr_->computeXyceFDVectors(theInstance.solutionFDVars_,
                                                    frequency,
                                                    theInstance.fDFVec_,
                                                    theInstance.fDBVec_,
                                                    theInstance.dFdXFDMat_)
          && bsuccess;
      }
    }
  }
  return bsuccess;
}
//-----------------------------------------------------------------------------
// Function      : Master::loadFreqDAEVectors
// Purpose       : Load frequency domain DAE vectors from stored intermediates
// Special Notes :
// Scope         : public
// Creator       : Tom Russo and Ting Mei
// Creation Date : 22 March 2018
//-----------------------------------------------------------------------------
bool Master::loadFreqDAEVectors(double frequency, std::complex<double>* solVec,
                                  std::vector<Util::FreqVecEntry>& fVec,
                                  std::vector<Util::FreqVecEntry>& bVec)
{
  InstanceVector::const_iterator it, end;

  it = getInstanceBegin();
  end = getInstanceEnd();

  fVec.clear();
  bVec.clear();

  for ( ; it != end; ++it )
  {
    Instance & theInstance = *(*it);
    Util::FreqVecEntry tmpEntry;

    for (int i =0; i<theInstance.fDFVec_.size(); i++)
    {
      tmpEntry.val = theInstance.fDFVec_[i];
      tmpEntry.lid = theInstance.li_Nodes_[i];
      fVec.push_back(tmpEntry);
    }
    for (int i =0; i<theInstance.fDBVec_.size(); i++)
    {
      tmpEntry.val = theInstance.fDBVec_[i];
      tmpEntry.lid = theInstance.li_Nodes_[i];
      bVec.push_back(tmpEntry);
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadFreqDAEMatrices
// Purpose       : Load frequency domain DAE matrices from stored intermediates
// Special Notes :
// Scope         : public
// Creator       : Tom Russo and Ting Mei
// Creation Date : 22 March 2018
//-----------------------------------------------------------------------------
bool Master::loadFreqDAEMatrices(double frequency, std::complex<double>* solVec,
                                   std::vector<Util::FreqMatEntry>& dFdx)
{
  InstanceVector::const_iterator it, end;

  it = getInstanceBegin();
  end = getInstanceEnd();

  dFdx.clear();

  for ( ; it != end; ++it )
  {
    Instance & theInstance = *(*it);
    Util::FreqMatEntry tmpEntry;

    for (int row =0; row<theInstance.dFdXFDMat_.size(); row++)
    {
      for (int col =0; col<theInstance.A_Equ_NodeOffsets_[row].size(); col++)
      {
        tmpEntry.val = theInstance.dFdXFDMat_[row][theInstance.jacStamp_[row][col]];
        tmpEntry.row_lid = theInstance.li_Nodes_[row];
        tmpEntry.col_lid = theInstance.A_Equ_NodeOffsets_[row][col];
        dFdx.push_back(tmpEntry);
      }
    }
  }
  return true;
}

// Model class functions

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
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
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
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
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
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
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
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
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i, isize;
  isize = instanceContainer.size();

  os << std::endl;
  os << "Number of GenExt instances: " << isize << std::endl;
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
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
/// Apply a device instance "op" to all instances associated with this
/// model
///
/// @param[in] op Operator to apply to all instances.
///
///
void Model::forEachInstance(DeviceInstanceOp &op) const
{
  for (std::vector<Instance *>::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
    op(*it);
}


Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() || (deviceMap.find("GENEXT")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("genext", 1)
      .registerModelType("genext", 1);
  }
}

} // namespace GeneralExternal
} // namespace Device
} // namespace Xyce

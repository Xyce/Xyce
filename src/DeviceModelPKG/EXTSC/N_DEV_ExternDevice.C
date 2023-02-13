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
// Creation Date  : 04/09/05
//
//
//
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

#include <sstream>
#include <N_UTL_Math.h>

#include <N_DEV_1D_Interface.h>
#include <N_DEV_CharonInterface.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_ExternDevice.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Message.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_XyceInterface.h>
#include <N_DEV_fwd.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_PDS_Comm.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<VoltageNode>::ParametricData()
{
  // Set up map for double precision variables:
  addPar("INITVAL", 0.0, &VoltageNode::initVal);
  addPar("LIMITHIGH", 0.0, &VoltageNode::limValHigh);
  addPar("LIMITLOW", 0.0, &VoltageNode::limValLow);

  // Set up exceptions (ie variables that are not doubles):
  addPar("NAME", "VCONNECT0000", &VoltageNode::vsrcName);
}

ParametricData<VoltageNode> &VoltageNode::getParametricData() {
  static ParametricData<VoltageNode> parMap;

  return parMap;
}


//-----------------------------------------------------------------------------
// Function      : VoltageNode::VoltageNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/03/06
//-----------------------------------------------------------------------------
VoltageNode::VoltageNode()
  : CompositeParam(getParametricData()),
    vsrcName(""),
    initVal(0.0),
    limValHigh(0.15),
    limValLow(0.10)
{
}

//-----------------------------------------------------------------------------
// Function      : VoltageNode::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/03/06
//-----------------------------------------------------------------------------
void VoltageNode::processParams ()
{}


namespace ExternDevice {

void Traits::loadInstanceParameters(ParametricData<ExternDevice::Instance> &p)
{
  // Set up configuration constants:
// Set up map for double precision variables:

  // Set up map for non-double precision variables:
  p.addPar("EXTERNCODE", std::string("xyce"), &ExternDevice::Instance::externCode_);
  p.addPar("NETLIST", std::string("input.cir"), &ExternDevice::Instance::netlistFilename_);
  p.addPar("VOLTLIM", false, &ExternDevice::Instance:: voltageLimiterFlag);

  p.addComposite("NODE", VoltageNode::getParametricData(), &ExternDevice::Instance::voltageNodeVec);
}

void Traits::loadModelParameters(ParametricData<ExternDevice::Model> &p)
{}


//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Instance::processParams()
{

  std::string msg;

  if ( externCode_ == "xyce" )
  {
    if (extCodePtr_==0)
    {
      extCodePtr_ = new XyceInterface(getDeviceOptions(), commandLine_, netlistFilename_);
    }
  }
  else if ( externCode_ == "1d")
  {
    UserError(*this) << "1D interface not suppported";
  }
  else if ( externCode_ == "charon")
  {
    if (extCodePtr_==0)
    {
      extCodePtr_ = new CharonInterface(getDeviceOptions(), netlistFilename_, getSolverState());
    }
  }
  else
  {
    UserError(*this) << externCode_ << " is not a recognized external device";
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &       instance_block,
  Model &                     model,
  const FactoryBlock &        factory_block)
  : DeviceInstance(instance_block, configuration.getInstanceParameters(), factory_block),
    externCode_("xyce"),
    netlistFilename_("input.cir"),
    extCodePtr_(0),
    model_(model),
    initializeFlag_(false),
    innerSolveStatus_(false),
    voltageLimiterFlag(false),
    initJctGiven_ (false),
    nodeGiven_(false),
    locallyConnected_(true),
    owningProc_(-1),
    comm_(0),
    commandLine_(factory_block.commandLine_),
    instanceBlock_(instance_block)
{

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "name             = " << getName() << std::endl;
  }

  // allow for a variable numbers of nodes
  numExtVars = instance_block.numExtVars;
  numStateVars = instance_block.numExtVars;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "numExtVars = " << numExtVars << std::endl;

    std::vector<Param>::const_iterator iter_t;
    std::vector<Param>::const_iterator begin_t = instance_block.params.begin();
    std::vector<Param>::const_iterator end_t   = instance_block.params.end();

    for (iter_t = begin_t; iter_t !=  end_t;  ++iter_t)
    {
      Xyce::dout() << "Tag: " << iter_t->tag();
      Xyce::dout() << "  Value = " << iter_t->stringValue() << "  Given: " << iter_t->given() << std::endl;
    }
  }

  int size = numExtVars;
  if( jacStamp.empty() )
  {
    jacStamp.resize(size);
    for( int i=0; i<size ; i++ )
    {
      jacStamp[i].resize(size);
      for( int j=0; j<size; j++ )
      {
        jacStamp[i][j] = j;
      }
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    for( int i=0; i<size; i++)
    {
      for( int j=0; j<size; j++)
      {
        Xyce::dout() << "jacStamp[ " << i << ", " << j << "] = "
             << jacStamp[i][j] << std::endl;
      }
    }
    Xyce::dout() << std::endl;
  }

  currentOutputVector_.resize(size,0.0);
  conductanceJacobian_.resize(size);
  for (int i=0;i<size;++i)
  {
    conductanceJacobian_[i].resize(size,0.0);

    // This name MUST be all caps!
    std::ostringstream oss;
    oss << "VCONNECT" << std::setw(4) << std::setfill('0') << i;
    voltageInputMap_[oss.str()] = 0.0;
  }

  // voltlim vectors:
  voltageOld_.resize      (size,0.0);
  voltageOrig_.resize     (size,0.0);
  voltageLastCall_.resize (size,0.0);
  voltageDiff_.resize     (size,0.0);
  voltageFactor_.resize   (size,0.0);
  voltageStateID_.resize  (size,0);


  // Set params to constant default values:

  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:

  setParams (instance_block.params);

  // Set any non-constant parameter defaults:

  if (!given("EXTERNCODE"))
    externCode_ = "xyce";
  if (!given("NETLIST"))
    netlistFilename_ = "input.cir";

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
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
Instance::~Instance()
{
  if (extCodePtr_!=0)
  {
    delete extCodePtr_ ;
    extCodePtr_ = 0;
  }

  // delete the pointers that are in the voltage node map.
  std::vector<VoltageNode *>::iterator iterNode = voltageNodeVec.begin();
  std::vector<VoltageNode *>::iterator endNode  = voltageNodeVec.end ();
  for (;iterNode!=endNode;++iterNode)
  {
    delete (*iterNode);
    (*iterNode) = 0;
  }
}

// Additional Declarations
//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       : function for registering, and setting up, local ID's.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter , SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "In Instance::register LIDs\n\n";
    Xyce::dout() << "name             = " << getName() << std::endl;
    Xyce::dout() << "number of internal variables: " << numIntVars << std::endl;
    Xyce::dout() << "number of external variables: " << numExtVars << std::endl;
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "\nEnd of Instance::register LIDs\n";
    Xyce::dout() << section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter , SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  int i=0;
  for ( ; i<numStateVars;++i)
  {
    voltageStateID_[i] = staLIDVecRef[i];
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter , SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
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
// Creator       : Eric R. Keiter , SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  jacLIDs = jacLIDVec;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 external device instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 external device instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  Linear::Vector * daeFVecPtr = extData.daeFVectorPtr;
  Linear::Vector * dFdxdVpPtr = extData.dFdxdVpVectorPtr;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "In ::loadDAEFVector " << getName() << std::endl;
  }

  if( locallyConnected_ )
  {
    // sum the current vector into the RHS:
    for (int i=0;i<currentOutputVector_.size();++i)
    {
      int iRow = extLIDVec[i];

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
      {
        Xyce::dout() <<"Instance: iRow = " << iRow << "  current = ";
        Xyce::dout() << currentOutputVector_[i] << std::endl;
      }
      (*daeFVecPtr)[iRow] -= currentOutputVector_[i];

      if (voltageLimiterFlag && getDeviceOptions().voltageLimiterFlag && !origFlag)
      {
        (*dFdxdVpPtr)[iRow] += voltageFactor_[i];
      }
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Done with Instance::loadDAEFVector "<< std::endl;
  }


  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 external device instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 external device instance.
//
// Special Notes : This is an algebraic constaint, and as such the external device
//                 does make a contribution to it.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;
  Linear::Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  if( locallyConnected_ )
  {
    // Load in the conductance Jacobian:
    int iRowLoc,iColLoc;
    int size = conductanceJacobian_.size();
    for (iRowLoc=0;iRowLoc < size; ++iRowLoc)
    {
      for (iColLoc=0;iColLoc < size; ++iColLoc)
      {
        int iRow = extLIDVec[iRowLoc];
        int iCol = jacLIDs[iRowLoc][iColLoc];
        double val = conductanceJacobian_[iRowLoc][iColLoc];
        (*dFdxMatPtr)[iRow][iCol] += val;
      }
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       : update primary state for one external device
//                 instance
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState()
{
  bool bsuccess=true;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "  Instance::updatePrimaryState" << std::endl;
    Xyce::dout() << "  name   = " << getName() << std::endl;
  }

//  bsuccess = bsuccess && updateIntermediateVars();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       : update secondary state for one external device
//                 instance
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/06/04
//-----------------------------------------------------------------------------
bool Instance::getBreakPoints(
    std::vector<Util::BreakPoint> & breakPointTimes,
    std::vector<Util::BreakPoint> & pauseBreakPointTimes)
{
  // if necessary, initialize
  initialize();
  bool bsuccess = extCodePtr_->getBreakPoints(breakPointTimes, pauseBreakPointTimes);
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::runExternalDevice
// Purpose       : actually run the solve of the external device
//
// Special Notes : updateIntermediateVars is called here, instead of
//                 in updatePrimaryState, b/c these devices need to do
//                 synchronized loads in parallel.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Elec. & Micro Systems
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool Instance::runExternalDevice()
{
  bool bsuccess = true;
  bsuccess = bsuccess && updateIntermediateVars();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupVoltageLimiting_
// Purpose       :
// Special Notes : The limits are hardwired for now, but should be made to
//                 be more flexible later.  Specifically, the limits should
//                 take advantage of gradient information.
//
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 03/24/06
//-----------------------------------------------------------------------------
bool Instance::setupVoltageLimiting_ ()
{
  int i(0);

  if (DEBUG_DEVICE)
    Xyce::dout() << "In setupVoltageLimiting" << std::endl;

  origFlag=true;

  std::map<std::string,double>::iterator iterM = voltageInputMap_.begin();
  std::map<std::string,double>::iterator endM  = voltageInputMap_.end ();

  int voltSize = voltageOrig_.size();
  // save the orig, unlimited voltages:
  for(i=0;i<voltSize;++i,++iterM)
  {
    voltageOrig_[i] = iterM->second;
    if (DEBUG_DEVICE)
      Xyce::dout() << "voltageOrig_["<<i<<"] = " << voltageOrig_[i] << std::endl;
  }

  if (DEBUG_DEVICE)
  {

    Xyce::dout() << getName() << "  Before Limiting:"<< std::endl;
    iterM = voltageInputMap_.begin();
    for( ; iterM != endM; ++iterM)
    {
      Xyce::dout() << "Blim: ";
      Xyce::dout() << iterM->first;
      Xyce::dout() << "\t";
      Xyce::dout() << iterM->second;
      Xyce::dout() << std::endl;
    }
  }

  // grab the "old" voltages from either the "last call" vector, or the
  // state vector.
  if (getSolverState().newtonIter == 0)
  {
    if (DEBUG_DEVICE)
    {
      Xyce::dout() << "Inside the newtonIter==0... " << std::endl;
    }

    if (!(getSolverState().dcopFlag)||(getSolverState().locaEnabledFlag && getSolverState().dcopFlag))
    {
      if (DEBUG_DEVICE)
      {
        Xyce::dout() << "Inside !dcopFlag || (locaEnabledFlago && dcopFlag) " << std::endl;
      }
      for(i=0;i<voltSize;++i)
      {
        int li_state = voltageStateID_[i];
        //voltageOld_[i] = (*(*extData.currStaVectorPtrPtr))[li_state];
        voltageOld_[i] = (*extData.currStaVectorPtr)[li_state];
      }
    }
    else  // no history, just set old voltages to the current ones.
    {
      if (DEBUG_DEVICE)
      {
        Xyce::dout() << "Inside no-history" << std::endl;
      }
      iterM = voltageInputMap_.begin();
      for(i=0;i<voltSize;++i,++iterM)
      {
        voltageOld_[i] = iterM->second;
      }
    }
  }
  else 
  {
    if (DEBUG_DEVICE)
    {
      Xyce::dout() << "Inside lastCall. " << std::endl;
    }
    for(i=0;i<voltSize;++i)
    {
      voltageOld_[i] = voltageLastCall_[i];
    }
  }

  bool setInitJct = false;
  if (getSolverState().initJctFlag_ && getDeviceOptions().voltageLimiterFlag)
  {
    if (DEBUG_DEVICE)
    {
      Xyce::dout() << "Inside the initJctFlag==true." << std::endl;
    }
    if (initJctGiven_)
    {
      if (DEBUG_DEVICE)
      {
        Xyce::dout() << "Inside the local initJctGiven_==true." << std::endl;
      }
      origFlag=false;
      setInitJct = true;

      iterM = voltageInputMap_.begin();
      for(i=0;i<voltSize;++i,++iterM)
      {
        iterM->second = voltageNodeVec[i]->initVal;
        if (DEBUG_DEVICE)
        {
          Xyce::dout() << "initVal : " << iterM->first << " = ";
          Xyce::dout() << voltageNodeVec[i]->initVal << std::endl;
        }
      }
    }
  }

  if (!setInitJct)
  {
    iterM = voltageInputMap_.begin();
    endM = voltageInputMap_.end ();
    for(i=0; iterM != endM; ++iterM, ++i)
    {
      double v1     = iterM->second;
      double v1_old = voltageOld_[i];
      double delV1 = v1 - v1_old;

      // Note: For now we are using hardwired values.
      // In the future, this should
      // be a more sophisticated limiter(s).
      //double limValHigh = voltageNodeVec[i].limValHigh;
      double limValHigh = 0.15;
      if ( delV1 > +limValHigh)
      {
        v1 = v1_old + limValHigh;
        origFlag=false;
      }
      //double limValLow = voltageNodeVec[i].limValLow;
      double limValLow = 0.10;
      if ( delV1 < -limValLow)
      {
        v1 = v1_old - limValLow;
        origFlag=false;
      }
      iterM->second = v1;
    }
  }

  // save the junction voltages in case the newton iteration is about
  // to be declared as successful.  If so, these "lastcall" variables will
  // get copied into their "old" versions.
  iterM = voltageInputMap_.begin();
  for(i=0;i<voltSize;++i,++iterM)
  {
    voltageLastCall_[i] = iterM->second;
  }

  if (DEBUG_DEVICE)
  {
    Xyce::dout() << getName() << "  After Limiting:"<< std::endl;
    iterM = voltageInputMap_.begin();
    for( ; iterM != endM; ++iterM)
    {
      Xyce::dout() << "Alim: ";
      Xyce::dout() << iterM->first;
      Xyce::dout() << "\t";
      Xyce::dout() << iterM->second;
      Xyce::dout() << std::endl;
    }
    if (origFlag)
      Xyce::dout() << "origFlag is TRUE" << std::endl;
    else
      Xyce::dout() << "origFlag is FALSE" << std::endl;
  }

  iterM = voltageInputMap_.begin();
  for(i=0;i<voltSize;++i,++iterM)
  {
    voltageDiff_[i] = iterM->second - voltageOrig_[i];
  }

  // Finally, save voltages in the state vector:
  Linear::Vector * staVectorPtr;
  //staVectorPtr = *(extData.nextStaVectorPtrPtr);
  staVectorPtr = extData.nextStaVectorPtr;
  for(i=0;i<voltSize;++i)
  {
    int li_state = voltageStateID_[i];
    (*staVectorPtr)[li_state] = voltageOld_[i];
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcVoltLimFactors_
//
// Purpose       : This function calculates the factors that will need to
//                 be loaded into the RHS vector to support voltage limiting.
//
// Special Notes : This function does a very small matvec, to produce
//                 jdxp terms.
//
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 03/24/06
//-----------------------------------------------------------------------------
bool Instance::calcVoltLimFactors_ ()
{
  bool bsuccess = true;

  int i=0;
  int j=0;
  int voltSize = voltageDiff_.size();
  for (;i<voltSize;++i)
  {
    double fac = 0.0;
    for (j=0;j<voltSize;++j)
    {
      fac += conductanceJacobian_[i][j] * voltageDiff_[j];
    }
    voltageFactor_[i] = fac;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::initialize
// Purpose       : If this device hasn't been initialized yet, do it here.
//                 If it has already been initialized, then this function is
//                 a no-op.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Elec. & Micro Systems
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool
Instance::initialize()  
{
  // If necessary, initialize the external device/code.
  if (!initializeFlag_)
  {
    if( locallyConnected_ )
    {
      // Do some initializations related to the voltageNodeVec.
      // I'm not sure where else to put this, so put it here.
      // bool nodeGiven = given("NODE");
      // vector composites don't seem to use the given function correctly...
      if (nodeGiven_)
      {
        // Check that the voltageNodeVec is the right size.
        int sizeVoltMap = voltageNodeVec.size();
        if (sizeVoltMap != numExtVars)
        {
          UserError(*this) << "Number of specified nodes != numExtVars";
        }

        // The name parameter was required, so force the voltageInputMap
        // to have the same keys as the voltageNodeVec.
        voltageInputMap_.clear();
        int i=0; int vnSize=voltageNodeVec.size();
        for (i=0;i<vnSize;++i)
        {
          voltageInputMap_[voltageNodeVec[i]->vsrcName] = 0.0;
        }

        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
        {
          std::map<std::string,double>::const_iterator iterVI = voltageInputMap_.begin();
          std::map<std::string,double>::const_iterator endVI = voltageInputMap_.end ();
          for (;iterVI!=endVI;++iterVI)
          {
            Xyce::dout() << "voltageInputMap["<<iterVI->first<<"] = ";
            Xyce::dout() << iterVI->second << std::endl;
          }
        }

        // Now check if initial values were specified.
        bool allInitSet = true;
        bool oneOrMoreInitSet = false;
        for (i=0;i<vnSize;++i)
        {
          VoltageNode & vn = *(voltageNodeVec[i]);
          bool tmpBool = vn.given("INITVAL");
          allInitSet = allInitSet && tmpBool;
          oneOrMoreInitSet = oneOrMoreInitSet || tmpBool;
        }

        if (oneOrMoreInitSet && !allInitSet)
        {
          UserError(*this) << "Initial value set on at least one node, but not all of them.  You must either set all or none.";
        }

        initJctGiven_ = allInitSet;
      }
      else // if the "node" vector-composite was *not* specified, we need to use the
           // defaults, and make sure that they are set up.
      {
        int sizeVoltMap = voltageNodeVec.size();
        if (sizeVoltMap != numExtVars)
        {
          int i=0;
          voltageNodeVec.resize(numExtVars);
          for (i=0;i<numExtVars;++i)
          {
            std::ostringstream oss;
            oss << "VCONNECT" << std::setw(4) << std::setfill('0') << i;
            voltageNodeVec[i] = new VoltageNode();
            voltageNodeVec[i]->vsrcName = oss.str();
          }
        }
      }
    }

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "Initializing internal Xyce" << std::endl;
    }
    extCodePtr_->initialize(comm_);

    // startTimeStep(
    //   beginIntegrationFlag,
    //   nextTimeStep,
    //   nextTime,
    //   currentOrder);

    initializeFlag_ = true;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::constructComposite
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/03/06
//-----------------------------------------------------------------------------
CompositeParam *Instance::constructComposite(const std::string & cName, const std::string & pName)
{
  CompositeParam * cpPtr = NULL;

  if (cName == "NODE")
  {
    nodeGiven_ = true;
    VoltageNode *vNodePtr = new VoltageNode ();
    voltageNodeVec.push_back(vNodePtr);
    cpPtr = (static_cast<CompositeParam *> (vNodePtr));
  }
  else
  {
    UserError(*this) << "Unrecognized composite name " << cName;
  }

  return cpPtr;
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
}


//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one external device
//                 instance
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars()
{
  Linear::Vector * solPtr = extData.nextSolVectorPtr;
  int i=0;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::updateIntermediateVars:" << std::endl;

    int vnSize = voltageNodeVec.size();

    for (i=0;i<vnSize;++i)
    {
      Xyce::dout() << voltageNodeVec[i];
    }
  }

  // If necessary, initialize the external device/code.
  initialize();

  innerSolveStatus_ = false;

  if (Parallel::is_parallel_run(comm_->comm()))
    if (comm_->procID() != owningProc_)
      locallyConnected_ = false;

  if( locallyConnected_ )
  {
    // Determine the input variables, etc., and place them in
    // the voltageInputMap
    for (i=0;i<extLIDVec.size();++i)
    {
      voltageInputMap_[voltageNodeVec[i]->vsrcName] = (*solPtr)[extLIDVec[i]];
    }

    if (voltageLimiterFlag && getDeviceOptions().voltageLimiterFlag)
    {
      setupVoltageLimiting_ ();
    }
  }

  if (Parallel::is_parallel_run(comm_->comm())) {

    if( comm_->procID() != owningProc_ ) locallyConnected_ = false;

    std::map<std::string,double>::iterator iterM = voltageInputMap_.begin();
    std::map<std::string,double>::iterator endM  = voltageInputMap_.end ();

    for( ; iterM != endM; ++iterM )
    {
      double val = iterM->second;
      double result = 0.0;
      if( !locallyConnected_ ) val = 0.0;
      comm_->sumAll( &val, &result, 1 );
      iterM->second = result;
    }
  }

  // Now, have the code perform a short simulation, and return currents
  // and conductances.
  //
  // had to change the first two arguments from
  //   currtime-currTimeStep, currTime
  // to
  //   currTime, currTime + currTimeStep
  // to avoid negative intial times and cases where this is
  // no time difference between the two arguments.
  //

  innerSolveStatus_ = extCodePtr_->simulateStep (
    getSolverState(),
    voltageInputMap_,
    currentOutputVector_,
    conductanceJacobian_,
    tlError_
                                                 );

  if( locallyConnected_ )
  {
    // Do some voltlim stuff here, if necessary.
    if (voltageLimiterFlag && getDeviceOptions().voltageLimiterFlag)
    {
      calcVoltLimFactors_ ();
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       : update intermediate variables for one external
//                 device instance
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Instance::updateTemperature( const double & temp )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::homotopyStepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void Instance::homotopyStepSuccess
(const std::vector<std::string> & paramNames,
 const std::vector<double> & paramVals)
{
  extCodePtr_->homotopyStepSuccess ( paramNames, paramVals);
}

//-----------------------------------------------------------------------------
// Function      : Instance::homotopyStepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void Instance::homotopyStepFailure ()
{
  extCodePtr_->homotopyStepFailure ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::stepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void Instance::stepSuccess(Analysis::TwoLevelMode analysis)
{
  extCodePtr_->stepSuccess(analysis);
}

//-----------------------------------------------------------------------------
// Function      : Instance::stepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void Instance::stepFailure(Analysis::TwoLevelMode analysis)
{
  extCodePtr_->stepFailure(analysis);
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInitialQnorm
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool Instance::getInitialQnorm (TimeIntg::TwoLevelError & tle)
{
  extCodePtr_->getInitialQnorm (tlError_);
  tle = tlError_;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInnerLoopErrorSum
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool Instance::getInnerLoopErrorSum (TimeIntg::TwoLevelError & tle)
{
  tle = tlError_;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateStateArrays
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool Instance::updateStateArrays ()
{
  bool bsuccess = true;
  if (extCodePtr_)
  {
    bool bs1 = extCodePtr_->updateStateArrays ();
    bsuccess = bsuccess && bs1;
  }
  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::startTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool Instance::startTimeStep(
  bool          beginIntegrationFlag,
  double        nextTimeStep,
  double        nextTime,
  int           currentOrder)
  
{
  bool bsuccess = true;
  if (extCodePtr_)
  {
    bool bs1 = extCodePtr_->startTimeStep(
      beginIntegrationFlag,
      nextTimeStep,
      nextTime,
      currentOrder);

    bsuccess = bsuccess && bs1;
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setInternalParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool Instance::setInternalParam(const std::string & name1, double val)
{
  bool bsuccess = true;
  initialize();
  if (extCodePtr_)
  {
    bool bs1 = extCodePtr_->setInternalParam(name1,val);
    bsuccess = bsuccess && bs1;
  }
  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Model::processParams ()
{


  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    model_block,
  const FactoryBlock &  factory_block)
  : DeviceModel(model_block, configuration.getModelParameters(), factory_block)
{
  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
Model::~Model()
{
  std::vector<Instance*>::iterator iterI;
  std::vector<Instance*>::iterator firstI = instanceContainer.begin ();
  std::vector<Instance*>::iterator lastI  = instanceContainer.end ();

  // loop over instances:
  for (iterI = firstI; iterI != lastI; ++iterI)
  {
    delete (*iterI);
  }


}

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
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


Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new DeviceMaster<Traits>(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("ext", 1)
    .registerModelType("ext", 1);
}

} // namespace ExternDevice
} // namespace Device
} // namespace Xyce


#include <iostream>
#include <unordered_map>
#include <string>

#include <mainXyceExpressionGroup.h>
#include <ast.h>
#include <newExpression.h>

#include <N_TIA_DataStore.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TOP_Topology.h>
#include <N_LAS_Vector.h>
#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Manager.h>
#include <N_PDS_Serial.h>

#include <N_ANP_AnalysisManager.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_Op.h>
#include <N_DEV_Const.h>

#include <N_ERH_ErrorMgr.h>

#include <N_UTL_DeviceNameConverters.h>

namespace Xyce {
namespace Util {

//-------------------------------------------------------------------------- 
//  THIS IS COPIED FROM OpBuilders.C.  Fix later ...
//-------------------------------------------------------------------------- 
// Function      : findNodeIndex 
// Purpose       : Determine whether a requested node name is "valid".
//                 That includes node names in the solution vector, the
//                 aliasNodeMap and also the Ground node (0).
//                 The return values are:
//                   a) -2 if the node is not found on this processor.
//                   b) -1 if it is a Ground node.
//                   c) the node index (0 ...N), otherwise.
// Special Notes :  
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
int findNodeIndex(
  const std::string &       name,
  const NodeNameMap &       node_map,
  const IO::AliasNodeMap &      alias_map)
{
  // The return value will be -2 if the specified node name is not
  // found AND the specified node name is not Ground (0). A value 
  // of -2 may be returned even for a valid node name, in parallel,
  // because of how node_map is used in parallel.
  int nodeIndex = -2;

  // Handle Ground (0).  If .PREPROCESS REPLACEGROUND was used,
  // then this function assumes that all of the GND nodes have
  // been replaced with 0 by this point in Xyce startup.
  if ( name == "0")
  { 
    // A node index of -1 will be used in the various get() functions
    // to denote the Ground node.
    nodeIndex = -1;
  } 
  else
  {
    NodeNameMap::const_iterator node_it = node_map.find(name);
    if (node_it == node_map.end()) 
    {
      // If the specified node name is not found in the node_map then
      // look for it in the AliasNodeMap. An example where this can 
      // happen is for a subcircuit interface node.
      IO::AliasNodeMap::const_iterator alias_node_it = alias_map.find(name);
      if (alias_node_it != alias_map.end())
      {
        // (*alias_node_it).second will be the "real name" of the alias node.
        node_it = node_map.find((*alias_node_it).second);
      }
    }

    // get the node index if it exists on this processor.
    if (node_it != node_map.end()) 
    {
      nodeIndex = (*node_it).second;
    } 
  }

  return nodeIndex;
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::mainXyceExpressionGroup 
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-------------------------------------------------------------------------------
mainXyceExpressionGroup::mainXyceExpressionGroup ( 
 N_PDS_Comm & comm, Topo::Topology & top,
 Analysis::AnalysisManager &analysis_manager,
 Device::DeviceMgr & device_manager,
 IO::OutputMgr &output_manager
 ) :
 comm_(comm),
 top_(top),
 analysisManager_(analysis_manager),
 deviceManager_(device_manager),
 outputManager_(output_manager),
 time_(0.0), temp_(0.0), VT_(0.0), freq_(0.0), dt_(0.0), alpha_(0.0)
{

}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::mainXyceExpressionGroup 
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-------------------------------------------------------------------------------
mainXyceExpressionGroup::~mainXyceExpressionGroup ()
{
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getSolutionGID_
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020
//-------------------------------------------------------------------------------
int mainXyceExpressionGroup::getSolutionGID_(const std::string & nodeName)
{
  int tmpGID=-1;
  std::vector<int> svGIDList1, dummyList;
  char type1;

  std::string nodeNameUpper = nodeName;
  Xyce::Util::toUpper(nodeNameUpper);

  bool foundLocal = top_.getNodeSVarGIDs(NodeID(nodeNameUpper, Xyce::_VNODE), svGIDList1, dummyList, type1);
  bool found = static_cast<int>(foundLocal);
  Xyce::Parallel::AllReduce(comm_.comm(), MPI_LOR, &found, 1);

  // if looking for this as a voltage node failed, try a "device" (i.e. current) node.  I(Vsrc)
  bool foundLocal2 = false;
  if (!found)
  {
    foundLocal2 = top_.getNodeSVarGIDs(NodeID(nodeNameUpper, Xyce::_DNODE), svGIDList1, dummyList, type1);
  }
  bool found2 = static_cast<int>(foundLocal2);
  Xyce::Parallel::AllReduce(comm_.comm(), MPI_LOR, &found2, 1);

  // Check if this is a subcircuit interface node name, which would be found in the aliasNodeMap.
  // If it is then get the GID for the corresponding "real node name". See SRN Bug 1962 for 
  // more details but, as an example, these netlist lines:
  //
  //   X1 1 2 MYSUB
  //   .SUBCKT MYSUB a c 
  //   R1   a b 0.5
  //   R2   b c 0.5
  //  .ENDS
  //
  // would produce key-value pairs of <"X1:C","2"> and <"X1:A","1"> in the aliasNodeMap 
  bool foundAliasNodeLocal = false;
  if (!found && !found2)
  {
    IO::AliasNodeMap::const_iterator alias_it = aliasNodeMap_.find( nodeNameUpper );
    if (alias_it != aliasNodeMap_.end())
    {      
      foundAliasNodeLocal = top_.getNodeSVarGIDs(NodeID((*alias_it).second, Xyce::_VNODE), svGIDList1, dummyList, type1);
    }
  }
  bool foundAliasNode = static_cast<int>(foundAliasNodeLocal);
  Xyce::Parallel::AllReduce(comm_.comm(), MPI_LOR, &foundAliasNode, 1);

  if(svGIDList1.size()==1)
  {
    tmpGID = svGIDList1.front();
  }

  return tmpGID;
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getSolutionVal
// Purpose       : 
// Special Notes : double precision version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 3/20/2020
//-------------------------------------------------------------------------------
bool mainXyceExpressionGroup::getSolutionVal(const std::string & nodeName, double & retval )
{
  retval = 0.0;

  int tmpGID = getSolutionGID_(nodeName);
  if (tmpGID >= 0)
  {
    const TimeIntg::DataStore & dataStore_ = *(analysisManager_.getDataStore());
    retval = dataStore_.nextSolutionPtr->getElementByGlobalIndex(tmpGID, 0);
  }
  Xyce::Parallel::AllReduce(comm_.comm(), MPI_SUM, &retval, 1);
  return (tmpGID>=0);
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getSolutionVal
// Purpose       : 
// Special Notes : std::complex<double> precision version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 3/20/2020
//-------------------------------------------------------------------------------
bool mainXyceExpressionGroup::getSolutionVal(const std::string & nodeName, std::complex<double> & retval)
{
  double real_val=0.0;
  double imag_val=0.0;

  int tmpGID = getSolutionGID_(nodeName);
  if (tmpGID >= 0)
  {
    const TimeIntg::DataStore & dataStore_ = *(analysisManager_.getDataStore());
    real_val = dataStore_.nextSolutionPtr->getElementByGlobalIndex(tmpGID, 0);
  }

  Xyce::Parallel::AllReduce(comm_.comm(), MPI_SUM, &real_val, 1);

  // ERK.  To Do:
  // need to add logic to see if this is frequency domain situation or not.
  // If not, don't bother getting the imaginary part.
  //
  // Xyce currently uses a real equivalent form for everything, rather than std::complex.
  //
  // For now, however, just set imag to zero.

  retval = std::complex<double>(real_val,imag_val);
  return (tmpGID>=0);
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getGlobalParameterVal
//
// Purpose       : retrieve the value of a parameter that has been 
//                 declared to be a "var" via the make_var function.
//
// Special Notes : double precision version
//
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020
//-------------------------------------------------------------------------------
bool mainXyceExpressionGroup::getGlobalParameterVal(const std::string &paramName, double & retval)
{
  bool success=true;
  Device::getParamAndReduce(comm_.comm(), deviceManager_, paramName, retval);
  return success;
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getGlobalParameterVal
//
// Purpose       : retrieve the value of a parameter that has been 
//                 declared to be a "var" via the make_var function.
//
// Special Notes : std::complex<double> version
//
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020
//-------------------------------------------------------------------------------
bool mainXyceExpressionGroup::getGlobalParameterVal (const std::string & paramName, std::complex<double> & retval)
{
  bool success=true;
  double tmpval;
  Device::getParamAndReduce(comm_.comm(), deviceManager_, paramName, tmpval);
  retval = std::complex<double>(tmpval,0.0);
  return success;
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getTimeStep
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getTimeStep ()
{
  dt_ = deviceManager_.getSolverState().currTimeStep_;
  return dt_;
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getTime
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getTime() 
{ 
  // I would have preferred to use this but as the code is currently written it
  // is not safe.  The earliest call I would need to make to getTime happens before 
  // the stepErrorControl class has been allocated.  Unfortunately, the analysis
  // manager accessor returns an invalid reference in that case, which I can't 
  // really test for.
  //
  //const TimeIntg::StepErrorControl & secControl_ = (analysisManager_.getStepErrorControl());
  //time_ = secControl_.nextTime;
  
  time_ = deviceManager_.getSolverState().currTime_;
  return time_;
} 

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getTemp
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getTemp() 
{ 
  temp_ = deviceManager_.getDeviceOptions().temp.getImmutableValue<double>() - CONSTCtoK;
  return temp_;
} 

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getVT  () 
{ 
  VT_ = (deviceManager_.getDeviceOptions().temp.getImmutableValue<double>())*CONSTKoverQ;
  return VT_;
} 

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getFreq() 
{ 
  freq_ = deviceManager_.getSolverState().currFreq_;
  return freq_;
} 

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getGmin() 
{ 
  gmin_ = deviceManager_.getGmin();
  return gmin_;
} 

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getBpTol()
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getBpTol()
{
  return deviceManager_.getSolverState().bpTol_;
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getStartingTimeStep
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/27/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getStartingTimeStep()
{
  return deviceManager_.getSolverState().startingTimeStep_;
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getFinalTime()
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/27/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getFinalTime()
{
  return deviceManager_.getSolverState().finalTime_;
}

} // Util
} // Xyce

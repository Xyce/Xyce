
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

namespace Xyce {
namespace Util {

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
 Device::DeviceMgr & device_manager
 ) :
 comm_(comm),
 top_(top),
 analysisManager_(analysis_manager),
 deviceManager_(device_manager),
 //tempOp_(new Device::ArtificialParameterOp("TEMP", deviceManager_, *(*deviceManager_.getArtificialParameterMap().find("TEMP")).second, "TEMP")),
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
  //delete tempOp_;
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
    //IO::AliasNodeMap::const_iterator alias_it = output_manager.getAliasNodeMap().find( nodeNameUpper );
    //IO::AliasNodeMap::const_iterator alias_it = aliasNodeMap_.find( nodeNameUpper );
    IO::AliasNodeMap::const_iterator alias_it = aliasNodeMapPtr_->find( nodeNameUpper );
    //if (alias_it != output_manager.getAliasNodeMap().end())
    //if (alias_it != aliasNodeMap_.end())
    if (alias_it != aliasNodeMapPtr_->end())
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
#if 0
  std::string nodeNameUpper = nodeName;
  Xyce::Util::toUpper(nodeNameUpper);

  // ERK.  Fix this so that these don't have to be looked up over and over.
  // get the GID for this variable
  {
    std::vector<int> svGIDList1, dummyList;
    char type1;

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
      //IO::AliasNodeMap::const_iterator alias_it = output_manager.getAliasNodeMap().find( nodeNameUpper );
      //IO::AliasNodeMap::const_iterator alias_it = aliasNodeMap_.find( nodeNameUpper );
      IO::AliasNodeMap::const_iterator alias_it = aliasNodeMapPtr_->find( nodeNameUpper );
      //if (alias_it != output_manager.getAliasNodeMap().end())
      //if (alias_it != aliasNodeMap_.end())
      if (alias_it != aliasNodeMapPtr_->end())
      {
        foundAliasNodeLocal = top_.getNodeSVarGIDs(NodeID((*alias_it).second, Xyce::_VNODE), svGIDList1, dummyList, type1);
      }
    }
    bool foundAliasNode = static_cast<int>(foundAliasNodeLocal);
    Xyce::Parallel::AllReduce(comm_.comm(), MPI_LOR, &foundAliasNode, 1);

    // once we have the GID, get value
    if(svGIDList1.size()==1)
    {
      tmpGID = svGIDList1.front();
    }
  }
#else
  int tmpGID = getSolutionGID_(nodeName);
#endif

  retval = 0.0;
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
  int tmpGID = getSolutionGID_(nodeName);

  double real_val=0.0;
  double imag_val=0.0;
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
#if 0
  // may not need to bother with the tempOp.  I copied it from the OutputMgrAdapter
  Util::Op::OpData op_data;
  temp_ = (*tempOp_)(comm_.comm(), op_data).real();
#else
  temp_ = deviceManager_.getDeviceOptions().temp.getImmutableValue<double>();
#endif
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
double mainXyceExpressionGroup::getBpTol()
{
  return deviceManager_.getSolverState().bpTol_;
}

} // Util
} // Xyce

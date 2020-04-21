
#include <iostream>
#include <unordered_map>
#include <string>

#include <mainXyceExpressionGroup.h>
#include <ast.h>
#include <newExpression.h>

#include <N_TIA_DataStore.h>
#include <N_TOP_Topology.h>
#include <N_LAS_Vector.h>
#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Manager.h>
#include <N_PDS_Serial.h>

#include <N_ANP_AnalysisManager.h>

namespace Xyce {
namespace Util {

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getSolutionVal
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 3/20/2020
//-------------------------------------------------------------------------------
bool mainXyceExpressionGroup::getSolutionVal(const std::string & nodeName, double & retval )
{
  bool success=true;
  int tmpGID=-1;
  retval = 0.0;
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

  if (tmpGID >= 0)
  {
    const TimeIntg::DataStore & dataStore_ = *(analysisManager_.getDataStore());
    retval = dataStore_.nextSolutionPtr->getElementByGlobalIndex(tmpGID, 0);
  }
  else 
  {
    retval = 0.0;
  }
  Xyce::Parallel::AllReduce(comm_.comm(), MPI_SUM, &retval, 1);

  return success; // FIX THIS
}

}
}

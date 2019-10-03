//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/20/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef N_TOP_Topology_h
#define N_TOP_Topology_h

// ---------- Standard Includes ----------
#include <string>
#include <list>
#include <map>
#include <set>
#include <iosfwd>
#include <vector>

#include <N_ANP_fwd.h>
#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TOP_fwd.h>
#include <N_LAS_fwd.h>
#include <N_UTL_fwd.h>

#include <N_DEV_DeviceBlock.h>
#include <N_UTL_Misc.h>
#include <N_UTL_NoCase.h>
#include <N_UTL_NodeSymbols.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : Topology
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class Topology
{
  friend std::ostream & operator << (std::ostream & os, const Topology & topo);

public:
  Topology(const IO::CmdParse &cp, IO::HangingResistor &hanging_resistor, N_PDS_Manager &pds_manager);
  ~Topology();

private:
  Topology(const Topology & right);
  Topology &operator=(const Topology & right);

public:
  // Setup linear system information.
  bool setupGlobalIndices();

  // Current add methods for devices tailored to work with device ghosting.
  void addDevice(Device::DeviceMgr &device_manager, const NodeDevBlock & nb);

  // Main call to register ext global id's with device nodes.
  void registerGIDs();

  // Main call to register local id's with device instances.
  void registerLIDswithDevs();

  // Resolution of Secondary Dependencies for devices.
  void resolveDependentVars();

  // Output to screen, debug, traversals of graph.
  void OutputBFSGraphLists() const;

  // Generate Ordered node list using BFS.
  void generateOrderedNodeList() const;

  const CktNodeList &getOrderedNodeList() const
  {
    return *orderedNodeListPtr_;
  }

  CktNodeList &getOrderedNodeList() 
  {
    return *orderedNodeListPtr_;
  }

  // merge the off processor superNodeList_ and communicate the same list
  // to all procs so topology reduction is the same on all procs
  void mergeOffProcTaggedNodesAndDevices();

  // this functions builds up the supernode list.  Called at the end of
  // parsing from n_cir_xyce
  void verifyNodesAndDevices(Device::DeviceMgr &device_manager);

 // remove any nodes and devices that were tagged for removal during parsing
  void removeTaggedNodesAndDevices();

  // Used for delayed instantiation of devices
  void instantiateDevices();

  // Return list of solution var indices for named node.  Returns false if node
  // not owned or not local.
  bool getNodeSVarGIDs( const NodeID &id,
                        std::vector<int> & sVarGIDList,
                        std::vector<int> & extSVarGIDList,
                        char & type ) const;

  // Accessor to get utility class for linear solver.
  Linear::QueryUtil *getLinearSolverUtility()
  {
    return linearSolverUtility_;
  }

  void regenerateGIDNodeMap();

  bool generateICLoader(Device::DeviceMgr &device_manager);

  // Restart capability.
  bool getRestartNodes(Analysis::AnalysisManager &analysis_manager, std::vector< IO::RestartNode * > & nodeVec);
  bool restoreRestartNodes(Analysis::AnalysisManager &analysis_manager, const std::vector< IO::RestartNode * > & nodeVec);

  // Solution variable name output.
  bool outputNameFile(Parallel::Machine comm, const std::string &path, bool overRideOutput = false);

  // These calls generate and return node maps.
  // NOTE:  The node type is known, this map will have unique keys.
  void loadNodeSymbols() const;

  const std::vector<char> &getVarTypes() const;

  const NodeNameMap &getSolutionNodeNameMap() const {
    loadNodeSymbols();
    return nodeSymbols_[Util::SOLUTION_SYMBOL];
  }

  const NodeNameMap &getStateNodeNameMap() const {
    loadNodeSymbols();
    return nodeSymbols_[Util::STATE_SYMBOL];
  }

  const NodeNameMap &getStoreNodeNameMap() const {
    loadNodeSymbols();
    return nodeSymbols_[Util::STORE_SYMBOL];
  }

  const NodeNameMap &getExternalNodeNameMap() const {
    loadNodeSymbols();
    return nodeSymbols_[Util::EXTERN_SYMBOL];
  }

  const NodeNameMap &getBranchVarsNodeNameMap() const {
    loadNodeSymbols();
    return nodeSymbols_[Util::BRANCH_SYMBOL];
  }

  // used to get the devices names (that have noise sources) 
  // from the Symbol Table
  const NodeNameMap &getNoiseDeviceNameMap() const {
    loadNodeSymbols();
    return nodeSymbols_[Util::NOISE_DEVICE_SYMBOL];
  }

  // used to get the noise type names (where each device may have
  // multiple noise sources/types) from the Symbol Table
  const NodeNameMap &getNoiseTypeNameMap() const {
    loadNodeSymbols();
    return nodeSymbols_[Util::NOISE_TYPE_SYMBOL];
  }
  
  // These functions are added to augment a copy of the netlist file to
  // include resistors which connect ground to "dangling" nodes.
  // NOTE: These IDs are assumed to be _VNODE, so no need to use NodeID.
  void writeNetlistAddResistors( unordered_set<std::string> & connToOneTermIDs,
                                 unordered_set<std::string> & noDCPathIDs );

  // Print out warnings if topological issues are found with the circuit.
  void outputTopoWarnings(unordered_set<std::string> & oneTermName,
                          unordered_set<std::string> & noDCName );

  // Returns adj IDs to the given ID
  void returnAdjIDs( const NodeID& id, std::vector<NodeID>& adj_ids ) const;

  // Returns adj GIDs to the given GID, including ground (-1)
  void returnAdjGIDsWithGround( int gid, std::vector<int>& adj_gids ) const;

  // Returns the number of adjacent nodes, without their GIDs.
  int numAdjNodes( int gid ) const;

  const std::map< int, int > &getDepSolnGIDMap() const {
    return depSolnGIDMap_;
  }

  N_PDS_Manager &getPDSManager() {
    return pdsManager_;
  }

  Util::SymbolTable &getNodeSymbols() {
    return nodeSymbols_;
  }

  const std::vector<const std::string *> &getSolutionNodeNames() const {
    return solutionNodeNames_;
  }
  
  CktGraph &getMainGraph() {
    return *mainGraphPtr_;
  }

  const CktGraph &getMainGraph() const {
    return *mainGraphPtr_;
  }

  // Returns pointer to specified ckt node.
  const CktNode *findCktNode(const NodeID& cnID) const;

private:
  const IO::CmdParse &                  commandLine_;           ///< Command line
  IO::HangingResistor &                 hangingResistor_;       ///< Hanging resistor repair
  Linear::QueryUtil *                   linearSolverUtility_;   ///< Utility class to extract alloc info for linear solver.

  N_PDS_Manager &                       pdsManager_;            ///< Parallel services manager object

  CktGraph *                            mainGraphPtr_;

  // Pointer to the time-integration manager object.
  mutable CktNodeList *                 orderedNodeListPtr_;

  std::map< int, int >                  depSolnGIDMap_;

  // this is a list of nodes to be supernoded
  // format is: nodeToBeReplaced, nodeToBeKept in each pair.
  // NOTE: this memory is cleared after the nodes are removed from the ordered node list. 
  std::vector< std::pair<NodeID, NodeID> >              superNodeList_;
  std::vector< NodeID >                                 badDeviceList_;

  mutable Util::SymbolTable                             nodeSymbols_;
  mutable std::vector<char>                             variableTypes_;
  mutable std::vector<const std::string *>              solutionNodeNames_;
  const std::string                                     gnd_;

  // Methods for writing out netlist file.
  void addResistors(const unordered_set<std::string> & inputVec, const std::string & resValue,
		    bool oneTermNotNoDCPath);

  void appendEndStatement();

  // Collect the set of names on processor 0.
  void generateGlobalNameSet( unordered_set<std::string> & local_names );
};

bool registerPkgOptionsMgr(Topology &topology, IO::PkgOptionsMgr &options_manager);

} // namespace Topo
} // namespace Xyce

#endif

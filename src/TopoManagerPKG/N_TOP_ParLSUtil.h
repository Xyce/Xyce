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
// Creation Date  : 05/26/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef N_TOP_ParLSUtil_h
#define N_TOP_ParLSUtil_h 1

// ---------- Standard Includes ----------

#include <list>
#include <iosfwd>
#include <vector>
#include <set>
#include <map>

#include <N_IO_fwd.h>
#include <N_TOP_fwd.h>
#include <N_PDS_fwd.h>

#include <N_LAS_QueryUtil.h>
#include <N_TOP_Topology.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : ParLSUtil
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/26/00
//-----------------------------------------------------------------------------
class ParLSUtil : public Linear::QueryUtil
{
  friend std::ostream & operator << (std::ostream & os, const ParLSUtil & tlsu);

public:
  ParLSUtil(Topology &topology, N_PDS_Manager &pds_manager);

  // Destructor
  ~ParLSUtil() {}

  // Generation of Linear system data:

  // Setup row/col data for linear solver including reorder.
  bool setupRowCol();

  // Accessor methods.

  const std::vector<int> & vnodeGIDVec() const { return vnodeGIDVector_; }
  const std::vector<int> & vsrcGIDVec() const { return vsrcGIDVector_; }
  const std::vector<int> & nonlinGIDVec() const { return nonlinGIDVector_; }

  int numGlobalNodes() const { return numGlobalNodes_; }
  int numLocalNodes() const { return numLocalNodes_; }

  int numGlobalRows() const { return numGlobalRows_; }
  int numLocalRows() const { return numLocalRows_; }
  int numExternRows() const { return numExternRows_; }
  int numGlobalExternRows() const { return numGlobalExternRows_; }
  int baseRowGID() const { return baseRowGID_; }
  const std::vector<int> & rowList_GID() const { return rowList_GID_; }
  const std::vector< std::pair<int,int> > & rowList_ExternGID() const { return rowList_ExternGID_; }

  int numGlobalStateVars() const { return numGlobalStateVars_; }
  int numLocalStateVars() const { return numLocalStateVars_; }
  int baseStateVarGID() const { return baseStateVarGID_; }
  const std::vector<int> & rowList_StateGID() const { return rowList_StateGID_; }

  int numGlobalStoreVars() const { return numGlobalStoreVars_; }
  int numLocalStoreVars() const { return numLocalStoreVars_; }
  int baseStoreVarGID() const { return baseStoreVarGID_; }
  const std::vector<int> & rowList_StoreGID() const { return rowList_StoreGID_; }
  
  int numGlobalLeadCurrentVars() const { return numGlobalLeadCurrentVars_; }
  int numLocalLeadCurrentVars() const { return numLocalLeadCurrentVars_; }
  int baseLeadCurrentVarGID() const { return baseLeadCurrentVarGID_; }
  const std::vector<int> & rowList_LeadCurrentGID() const { return rowList_LeadCurrentGID_; }

  int numGlobalNZs() const { return numGlobalNZs_; }
  int numLocalNZs() const { return numLocalNZs_; }
  const std::vector<int> & rowList_NumNZs() const { return rowList_NumNZs_; }
  const std::vector< std::vector<int> > & rowList_ColList() const { return rowList_ColList_; }

  const std::vector<char> & rowList_VarType() const { return topology_.getVarTypes(); }

  void cleanRowLists();

private:

  // Generate Ordering and Var GIDs.
  bool setupNodeGIDs();

  // Generate Ordering and Var GIDs
  bool setupSolnAndStateGIDs();

  // Generate row/col data for linear solver.
  void generateRowColData();

  // Iterate over the ordered list and extract all local and external nodes.
  void extractNodeGIDs();

  // Iterate over the ordered list from topology and fill all the local GID vectors.
  void extractAllGIDsFromTopology();

private:
  Topology &            topology_;                      ///< Topology
  N_PDS_Manager &       pdsManager_;                    ///< Parallel services manager object.

  // Number of global (across all processors) nodes in the topology.
  int numGlobalNodes_;
  // Number of local (on processor) nodes in the topology.
  int numLocalNodes_;
  int baseNodeGID_;
  std::vector<int> nodeList_GID_;
  std::vector< std::pair<int,int> > nodeList_ExternGID_;

  // Number of global (across all processors) rows in the linear system.
  int numGlobalRows_;
  // Number of local (on processor) rows in the linear system.
  int numLocalRows_;
  int numExternRows_;
  int numGlobalExternRows_;
  int baseRowGID_;
  std::vector<int> rowList_GID_;
  std::vector< std::pair<int,int> > rowList_ExternGID_;

  // Number of global (across all processors) state-variables associated with
  // the linear system.
  int numGlobalStateVars_;
  // Number of local (on processor) state-variables associated with the linear
  // system.
  int numLocalStateVars_;
  int baseStateVarGID_;
  std::vector<int> rowList_StateGID_;

  // Number of global (across all processors) store-variables associated with
  // the linear system.
  int numGlobalStoreVars_;
  // Number of local (on processor) store-variables associated with the linear
  // system.
  int numLocalStoreVars_;
  int baseStoreVarGID_;
  std::vector<int> rowList_StoreGID_;

  // Number of global (across all processors) leadcurrent-variables associated with
  // the linear system.
  int numGlobalLeadCurrentVars_;
  // Number of local (on processor) leadcurrent-variables associated with the linear
  // system.
  int numLocalLeadCurrentVars_;
  int baseLeadCurrentVarGID_;
  std::vector<int> rowList_LeadCurrentGID_;

  //Graph info for jacobian matrix
  int numGlobalNZs_;
  int numLocalNZs_;
  std::vector<int> rowList_NumNZs_;
  std::vector< std::vector<int> > rowList_ColList_;

  std::vector<int> vnodeGIDVector_;
  std::vector<int> vsrcGIDVector_;
  std::vector<int> nonlinGIDVector_;

  //Adding these lists to detect the IDs of nodes for which there is
  //no DC path to ground or for which the node is only connected to one device
  //terminal.
  unordered_set<std::string> noDCPathIDs_;
  unordered_set<std::string> connToOneTermIDs_;

private:

  //testing routine for problems with voltage node connectivity
  bool testVoltageNodeConnectivity_();
  void comm_boundaries (std::map<int, std::vector<int> > & gid_map,
                      std::vector<int> & actual_buf_in, std::vector<int> & actual_buf_out,
                      std::vector<int> & buf_len, std::vector<int> & buf_dest,
                      std::vector<int *> & buf_in, std::vector<int *> & buf_out, int mode);

  // Don't allow copy construction or assignment.
  // Copy constructor (private)
  ParLSUtil(const ParLSUtil & right);
  // Assignment operator (private).
  ParLSUtil & operator = (const ParLSUtil & right);
};

} // namespace Topo
} // namespace Xyce

#endif

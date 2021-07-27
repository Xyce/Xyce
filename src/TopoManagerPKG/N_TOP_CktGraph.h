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

#ifndef N_TOP_CktGraph_h
#define N_TOP_CktGraph_h 1

#include <string>
#include <ostream>
#include <list>
#include <vector>
#include <N_TOP_fwd.h>

#include <N_UTL_Misc.h>
#include <N_UTL_Graph.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : CktGraph
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class CktGraph
{
public:
  typedef Util::Graph<NodeID, CktNode *, int> Graph;

  // Default constructor.
  CktGraph()
    : id_()
  {}

  // Constructor
  CktGraph(const std::string& cgID)
  : id_(cgID)
  {}


private:
  CktGraph(const CktGraph &);
  CktGraph &operator=(const CktGraph &);

public:
  virtual ~CktGraph()
  {}

  // Inserts graph node for ckt node if does not exist (abstract).
  virtual void InsertNode(CktNode* node,
                          const std::vector<NodeID>& neighborList) = 0;

  // Returns pointer to specified ckt node (abstract).
  virtual CktNode* FindCktNode(const NodeID& cnID) = 0;

  // Produces list of ckt nodes in breadth-first traversal order (abstract).
  virtual CktNodeList* getBFSNodeList() = 0;

  // Returns the node list from the circuit graph without any specific ordering
  virtual const Graph::Data1Map &getNodeList() = 0;

  // Loop over nodes and register int and ext global ids with each device
  // (abstract).
  virtual void registerGIDs() = 0;

  // Loop over nodes and register int and ext global ids with each device
  virtual void registerLIDswithDevs( Indexor & indexor ) = 0;
  // Loop over nodes and register state local ids with each device
  virtual void registerStateLIDswithDevs( Indexor & indexor ) = 0;
  virtual void registerStoreLIDswithDevs( Indexor & indexor ) = 0;

  // Loop over nodes and register all dep. local ids with each device
  virtual void registerDepLIDswithDevs( Indexor & indexor ) = 0;
  // Loop over nodes and register branch local ids with each device
  virtual void registerBranchDataLIDswithDevs( Indexor & indexor ) = 0;

  // Registration of jacobian local offsets with devices
  virtual void registerJacLIDswithDevs( Indexor & indexor ) = 0;

  // Returns adj IDs to the given ID
  virtual void returnAdjIDs( const NodeID& id, std::vector<NodeID>& adj_ids,
                             bool withGnd = false ) = 0;

  // Returns adj GIDs to the given GID, including ground
  virtual void returnAdjGIDsWithGround( int gid, std::vector<int>& adj_gids ) = 0;

  virtual int numAdjNodes( int gid ) = 0;
  virtual int numAdjNodesWithGround( int gid ) = 0;

  // Redo global index node map (abstract).
  virtual void regenerateGIDNodeMap() = 0;
  
  // supernode given nodes
  virtual CktNode * replaceNode( const NodeID nodeToBeReplaced, const NodeID nodeToKeep ) = 0;
  virtual void removeRedundantDevices( std::vector< NodeID > & devicesToBeRemoved, std::vector< CktNode * > & removedDevices) = 0;
  virtual void removeNodes( const std::vector< NodeID > & nodesToBeRemoved, std::vector< CktNode * > & removedNodes ) = 0;

  virtual std::vector< Xyce::NodeID > analyzeDeviceNodeGraph(std::ostream & os) = 0;

  const std::string& get_id() const { return id_; }

protected:

  std::string id_;

  virtual std::ostream & put(std::ostream & os) const = 0;

  friend std::ostream& operator<<(std::ostream& os, const CktGraph& cg);
};

} // namespace Topo
} // namespace Xyce

#endif

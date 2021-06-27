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
// Creator        : Robert J. Hoekstra, SNL, Electrical & Microsystems
//
// Creation Date  : 08/10/06
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef N_TOP_CktGraphBasic_h
#define N_TOP_CktGraphBasic_h 1

#include <iosfwd>
#include <string>

#include <N_TOP_fwd.h>
#include <N_TOP_CktGraph.h>
#include <N_UTL_Misc.h>
#include <N_UTL_Graph.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : CktGraphBasic
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Electrical & Microsystems
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
class CktGraphBasic : public CktGraph
{
public:
  typedef Util::Graph<NodeID, CktNode *, int> Graph;

  // Default constructor.
  CktGraphBasic();

  // Constructor
  CktGraphBasic(const std::string & cgID);

private:
  CktGraphBasic(const CktGraphBasic &);
  CktGraphBasic &operator=(const CktGraphBasic &);

public:
  // Destructor
  ~CktGraphBasic();

  // Inserts graph node for ckt node if does not exist
  void InsertNode(CktNode * cktnode,
                  const std::vector<NodeID> & neighborList);

  // Returns pointer to specified ckt node
  CktNode * FindCktNode(const NodeID & cnID);

  // Breadth first, depth first, and basic traversals of graph
  //----------------------------------------------------------

  // Produces list of ckt nodes in breadth-first traversal order
  CktNodeList* getBFSNodeList();

  // Returns the node list from the circuit graph without any specific ordering
  const Graph::Data1Map &getNodeList() { return cktgph_.getData1Map(); }

  // assigns global ids to device int. and ext. vars
  //----------------------------------------------------------

  // Loop over nodes and register ext global ids with each device node
  void registerGIDs();

  // Loop over nodes and register int and ext global ids with each device
  void registerStateGIDswithDevs();

  // Loop over nodes and register int and ext local ids with each device
  void registerLIDswithDevs( Indexor & indexor );
  // Loop over nodes and register state local ids with each device
  void registerStateLIDswithDevs( Indexor & indexor );
  void registerStoreLIDswithDevs( Indexor & indexor );

  // Loop over nodes and register dep. local ids with each device
  void registerDepLIDswithDevs( Indexor & indexor );
  // Loop over nodes and register branch local ids with each device
  void registerBranchDataLIDswithDevs( Indexor & indexor );

  //Loop over nodes and register jacobian offsets with each device
  void registerJacLIDswithDevs( Indexor & indexor );

  // Returns vector of adj ids
  void returnAdjIDs( const NodeID & id, std::vector<NodeID> & adj_ids );

  // Returns adj GIDs to the given GID, including ground
  void returnAdjGIDsWithGround( int gid, std::vector<int>& adj_gids );
  
  // Returns number of adj gids
  int numAdjNodes( int gid );

  // Returns number of adj gids, including ground node
  int numAdjNodesWithGround( int gid );

  // Redo global index node map
  void regenerateGIDNodeMap();

  // supernode given nodes
  CktNode * replaceNode( const NodeID nodeToBeReplaced, const NodeID nodeToKeep );
  void removeRedundantDevices( std::vector< NodeID > & devicesToBeRemoved, std::vector< CktNode * > & removedDevices);
  void removeNodes( const std::vector< NodeID > & nodesToBeRemoved, std::vector< CktNode * > & removedNodes );

private:
  Graph                 cktgph_;                ///< Circuit graph pair = <id, node type>, CktNode = data
  CktNodeList           BFSNodeList_;           ///< List of ckt nodes in breadth-first traversal order.
  bool                  isModified_;            ///< Don't update traversals if not modified - flag.

  unordered_map<int,int>     indexToGID_;
  unordered_map<int,int>     gIDtoIndex_;
 
public:
  std::ostream & put(std::ostream & os) const;
};

} // namespace Topo
} // namespace Xyce

#endif

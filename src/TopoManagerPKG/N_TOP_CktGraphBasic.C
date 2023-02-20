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
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 08/10/06
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_DEV_DeviceBlock.h>
#include <N_ERH_ErrorMgr.h>
#include <N_PDS_Manager.h>
#include <N_TOP_CktGraphBasic.h>
#include <N_TOP_CktNode.h>
#include <N_TOP_CktNode_Dev.h>
#include <N_TOP_Indexor.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Functors.h>

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::CktGraphBasic
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
CktGraphBasic::CktGraphBasic()
 : CktGraph(""),
   isModified_(true)
{
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::CktGraphBasic
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
CktGraphBasic::CktGraphBasic(const std::string &cgID) 
 : CktGraph(cgID),
   isModified_(true)
{
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::~CktGraphBasic()
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
CktGraphBasic::~CktGraphBasic()
{
  CktNodeList::iterator it_nL = BFSNodeList_.begin();
  CktNodeList::iterator it_nL_end = BFSNodeList_.end();
  for( ; it_nL != it_nL_end; ++it_nL )
    delete *it_nL;
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::InsertNode
// Purpose       : Inserts graph node for ckt node if does not exist
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::InsertNode(CktNode* cktnode,
                               const std::vector<NodeID> &neighborList)
{
  if (DEBUG_TOPOLOGY)
    Xyce::dout() << "Inserting Node: " << cktnode->get_nodeID() << std::endl
                 << *cktnode;

  bool inserted = cktgph_.insertNode( cktnode->get_nodeID(),
                                      neighborList, cktnode );
  
  if( !inserted )
  {
    delete cktnode;
  }
  else
  {
    //------ Graph has been changed so traversals will change
    isModified_ = true;
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::FindCktNode
// Purpose       : Returns ptr to specified ckt node.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
CktNode* CktGraphBasic::FindCktNode (const NodeID &cnID)
{
  // If the node type cannot be determined before hand, try to see if there is
  // a unique node associated with the string in cnID.
  if (cnID.second != -1)
  {
    if( cktgph_.checkKey(cnID) )
      return cktgph_.getData(cnID);
  }
  else
  {
    // Loop through all the possible node types to see if this node exists.
    // Only return it if there is one unique node.  Otherwise, this operation
    // is not well defined.
    int numFound       = 0;
    CktNode* foundNode = 0;
    for (int i         = 0; i<_NUM_NODE_TYPES; ++i)
    {
      NodeID tmpID( cnID.first, i );
      if (cktgph_.checkKey(cnID) )
      {
        foundNode      = cktgph_.getData(cnID);
        numFound++;
      }
    }
    if (numFound == 1)
      return foundNode;
  }

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::getBFSNodeList
// Purpose       : Produces list of ckt nodes in breadth first traversal order
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
CktNodeList *
CktGraphBasic::getBFSNodeList()
{
  //------- Check if list not created or graph modified

  if( BFSNodeList_.empty() || isModified_ )
  {
    //------- Include disconnected parts for graph in traversal
    int numNodes = cktgph_.numNodes();

    if( numNodes )
    {
      int numLevels = cktgph_.generateBFT();

      if (DEBUG_TOPOLOGY)
        Xyce::dout() << "Number of levels in BFT : " << numLevels << std::endl;
    
      const std::vector<int> & bfsVec = cktgph_.getBFT();

      BFSNodeList_.clear();
      BFSNodeList_.reserve( bfsVec.size() );
      std::vector<int>::const_reverse_iterator bfs_rit = bfsVec.rbegin();
      std::vector<int>::const_reverse_iterator bfs_rit_end = bfsVec.rend();
      for( ; bfs_rit != bfs_rit_end; ++bfs_rit )
        BFSNodeList_.push_back( cktgph_.getData(*bfs_rit) );
    }

    isModified_ = false;
  }

  return &BFSNodeList_;
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::returnAdjIDs
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::returnAdjIDs( const NodeID & id, std::vector<NodeID> & adj_ids,
                                  bool withGnd )
{
  adj_ids.clear();

  const std::vector<NodeID> &adjIDs = cktgph_.getAdjacent(id);
  std::vector<NodeID>::const_iterator adj_it = adjIDs.begin();
  std::vector<NodeID>::const_iterator adj_end = adjIDs.end();

  for( ; adj_it != adj_end; adj_it++ )
  {
    if( withGnd )
      adj_ids.push_back( *adj_it );
    else if( adj_it->first != "0" ) 
      adj_ids.push_back( *adj_it );
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::returnAdjGIDsWithGround
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 8/14/18
//-----------------------------------------------------------------------------
void CktGraphBasic::returnAdjGIDsWithGround( int gid,
                                             std::vector<int> & adj_gids )
{
  adj_gids.clear();
  
  const std::vector<int>& adjIndices = cktgph_.getAdjacentRow( gIDtoIndex_[ gid ] );
  std::vector<int>::const_iterator adj_it = adjIndices.begin();
  std::vector<int>::const_iterator adj_end = adjIndices.end();

  for( ; adj_it != adj_end; adj_it++ )
  {
    adj_gids.push_back( indexToGID_[ *adj_it ] );
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::numAdjNodes
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
int CktGraphBasic::numAdjNodes( int gid )
{
  int count = 0; 

  unordered_map<int,int>::iterator it = gIDtoIndex_.find( gid ); 
  if ( it != gIDtoIndex_.end() )
  {    
    const std::vector<int>& adjIndices = cktgph_.getAdjacentRow( it->second );
    std::vector<int>::const_iterator adj_it = adjIndices.begin();
    std::vector<int>::const_iterator adj_end = adjIndices.end();

    for( ; adj_it != adj_end; adj_it++ )
    {
      if( indexToGID_[ *adj_it ] != -1 ) 
        count++;
    }
  }
  return count;
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::numAdjNodesWithGround
// Purpose       : This also counts the ground node as an adjacent node
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/7/2020
//-----------------------------------------------------------------------------
int CktGraphBasic::numAdjNodesWithGround( int gid )
{
  int count = 0;

  // Return number of adjacent nodes before or after GIDs are determined.
  int id = gid;
  if (gIDtoIndex_.size())
  {
    unordered_map<int,int>::iterator it = gIDtoIndex_.find( gid );
    if ( it!= gIDtoIndex_.end() )
      id = it->second;
    else
      id = -1;
  } 

  if ( id != -1 )
  {
    const std::vector<int>& adjIndices = cktgph_.getAdjacentRow( id );
    count = adjIndices.size();
  }

  return count;
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::registerGIDs
// Purpose       : Loop over nodes and register ext global ids
//                 with each device node
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::registerGIDs()
{
  const CktNodeList * tmpNodeList = getBFSNodeList();
  CktNodeList::const_iterator it_nL = tmpNodeList->begin();
  CktNodeList::const_iterator it_nL_end = tmpNodeList->end();

  const std::vector<int> & tmpBfsVec = cktgph_.getBFT();
  std::vector<int>::const_reverse_iterator it_bfs = tmpBfsVec.rbegin();

  // loop over nodes:
  for ( ; it_nL != it_nL_end; ++it_nL, ++it_bfs)
  {
    int type = (*it_nL)->type();
    if (type == _DNODE)
    {
      CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_nL);
      std::vector<int> & svGIDList = cktNodeDevPtr->get_ExtSolnVarGIDList();

      //------  Generate global id list for external variables

      // initialize the local svGIDList.
      const std::vector<int> &adjIDs = cktgph_.getAdjacentRow( *it_bfs );
      std::vector<int>::const_iterator adj_it = adjIDs.begin();
      std::vector<int>::const_iterator adj_end = adjIDs.end();

      for ( ; adj_it != adj_end; adj_it++ )
      {
        CktNode * cktnode = cktgph_.getData(*adj_it);
        svGIDList.insert( svGIDList.end(), cktnode->get_SolnVarGIDList().begin(), cktnode->get_SolnVarGIDList().end() );
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::registerLIDswithDevs
// Purpose       : Loop over nodes and register int and ext local ids
//                 with each device
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::registerLIDswithDevs( Indexor & indexor )
{
  std::vector<int> svGIDList;
  std::vector<int>::const_iterator it_iL, end_iL;
  std::vector<int> intVec, extVec;

  CktNodeList * tmpNodeList = getBFSNodeList();

  CktNodeList::const_iterator it_nL = tmpNodeList->begin();
  CktNodeList::const_iterator it_nL_end = tmpNodeList->end();

  // loop over nodes:
  for( ; it_nL != it_nL_end; ++it_nL )
  {
    int type = (*it_nL)->type();
    if (type == _DNODE)
    {
      CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_nL);

      //------- Clear lists for each node
      svGIDList.clear();

      //------- Push list of int. and ext. local id's to ckt node

      // first arg. is internal variable list
      // 2nd   arg. is external variable list

      it_iL = cktNodeDevPtr->get_SolnVarGIDList().begin();
      end_iL = cktNodeDevPtr->get_SolnVarGIDList().end();
      intVec.assign( it_iL, end_iL );

      bool success = indexor.globalToLocal(Parallel::SOLUTION_OVERLAP_GND, intVec);

      it_iL = cktNodeDevPtr->get_ExtSolnVarGIDList().begin();
      end_iL = cktNodeDevPtr->get_ExtSolnVarGIDList().end();
      extVec.assign( it_iL, end_iL );

      success = success && indexor.globalToLocal(Parallel::SOLUTION_OVERLAP_GND, extVec);

      cktNodeDevPtr->registerLIDswithDev( intVec, extVec );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::registerStateLIDswithDevs
// Purpose       : Loop over nodes and register int and ext global ids
//                 with each device
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::registerStateLIDswithDevs( Indexor & indexor )
{
  std::vector<int>::const_iterator it_iL, end_iL;
  std::vector<int> stateVec;

  CktNodeList * tmpNodeList = getBFSNodeList();

  CktNodeList::const_iterator it_nL = tmpNodeList->begin();
  CktNodeList::const_iterator it_nL_end = tmpNodeList->end();

  // loop over nodes:
  for( ; it_nL != it_nL_end; ++it_nL )
  {
    int type = (*it_nL)->type();
    if (type == _DNODE)
    {
      CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_nL);

      it_iL = cktNodeDevPtr->get_StateVarGIDList().begin();
      end_iL = cktNodeDevPtr->get_StateVarGIDList().end();

      stateVec.assign( it_iL, end_iL );

      indexor.globalToLocal( Parallel::STATE, stateVec );

      //------- Push list of state global id's to ckt node
      cktNodeDevPtr->registerStateLIDswithDev( stateVec );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::registerStoreLIDswithDevs
// Purpose       : Loop over nodes and register int and ext global ids
//                 with each device
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void CktGraphBasic::registerStoreLIDswithDevs( Indexor & indexor )
{
  std::vector<int>::const_iterator it_iL, end_iL;
  std::vector<int> storeVec;

  CktNodeList * tmpNodeList = getBFSNodeList();

  CktNodeList::const_iterator it_nL = tmpNodeList->begin();
  CktNodeList::const_iterator it_nL_end = tmpNodeList->end();

  // loop over nodes:
  for( ; it_nL != it_nL_end; ++it_nL )
  {
    int type = (*it_nL)->type();
    if (type == _DNODE)
    {
      CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_nL);

      it_iL = cktNodeDevPtr->get_StoreVarGIDList().begin();
      end_iL = cktNodeDevPtr->get_StoreVarGIDList().end();

      storeVec.assign( it_iL, end_iL );

      indexor.globalToLocal(Parallel::STORE, storeVec);

      //------- Push list of store global id's to ckt node
      cktNodeDevPtr->registerStoreLIDswithDev( storeVec );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::registerBranchDataLIDswithDevs
// Purpose       : Loop over nodes and register lids for branch data vec.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void CktGraphBasic::registerBranchDataLIDswithDevs( Indexor & indexor )
{
  std::vector<int>::const_iterator it_iL, end_iL;
  std::vector<int> branchDataVec;

  CktNodeList * tmpNodeList = getBFSNodeList();

  CktNodeList::const_iterator it_nL = tmpNodeList->begin();
  CktNodeList::const_iterator it_nL_end = tmpNodeList->end();

  // loop over nodes:
  for( ; it_nL != it_nL_end; ++it_nL )
  {
    int type = (*it_nL)->type();
    if (type == _DNODE)
    {
      CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_nL);

      it_iL = cktNodeDevPtr->get_LeadCurrentVarGIDList().begin();
      end_iL = cktNodeDevPtr->get_LeadCurrentVarGIDList().end();

      branchDataVec.assign( it_iL, end_iL );

      indexor.globalToLocal(Parallel::LEADCURRENT, branchDataVec);

      //------- Push list of store global id's to ckt node
      cktNodeDevPtr->registerLeadCurrentLIDswithDev( branchDataVec );
    }
  }
}



//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::registerDepLIDswithDevs
// Purpose       : Loop over nodes and register dep local ids
//                 with each device
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::registerDepLIDswithDevs( Indexor & indexor )
{
  std::vector< std::vector<int> > indexVec;

  CktNodeList * tmpNodeList = getBFSNodeList();

  CktNodeList::iterator it_nL = tmpNodeList->begin();
  CktNodeList::iterator it_nL_end = tmpNodeList->end();

  // loop over nodes:
  for( ; it_nL != it_nL_end; ++it_nL )
  {
    int type = (*it_nL)->type();
    if (type == _DNODE)
    {
      CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_nL);

      // Register solution LIDs.
      if ( cktNodeDevPtr->depSolnVarCount() )
      {
        cktNodeDevPtr->get_DepSolnGIDVec( indexVec );

        for( unsigned int i = 0; i < indexVec.size(); ++i )
          indexor.globalToLocal(Parallel::SOLUTION_OVERLAP_GND, indexVec[i] );

        cktNodeDevPtr->registerDepLIDswithDev( indexVec );
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::registerJacLIDswithDevs
// Purpose       : Loop over nodes and register jacobian local ids
//                 with each device
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::registerJacLIDswithDevs( Indexor & indexor )
{
  std::vector< std::vector<int> > stampVec;

  indexor.setupAcceleratedMatrixIndexing(Parallel::JACOBIAN_OVERLAP);

  CktNodeList * tmpNodeList = getBFSNodeList();

  CktNodeList::const_iterator it_nL = tmpNodeList->begin();
  CktNodeList::const_iterator it_nL_end = tmpNodeList->end();

  // loop over nodes:
  for( ; it_nL != it_nL_end; ++it_nL )
  {
    int type = (*it_nL)->type();
    if (type == _DNODE)
    {
      CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_nL);

      const std::vector<int> & intGIDs = cktNodeDevPtr->get_SolnVarGIDList();
      const std::vector<int> & extGIDs = cktNodeDevPtr->get_ExtSolnVarGIDList();
      const std::vector<int> & depGIDs = cktNodeDevPtr->get_DepSolnGIDJacVec();
      std::vector<int> gids( intGIDs.size() + extGIDs.size() + depGIDs.size() );
      std::copy( extGIDs.begin(), extGIDs.end(), gids.begin() );
      std::copy( intGIDs.begin(), intGIDs.end(), gids.begin() + extGIDs.size() );
      std::copy( depGIDs.begin(), depGIDs.end(), gids.begin() + extGIDs.size() + intGIDs.size() );

      stampVec = cktNodeDevPtr->jacobianStamp();

      int numRows = stampVec.size();
      for( int i = 0; i < numRows; ++i )
      {
        int numCols = stampVec[i].size();
        for( int j = 0; j < numCols; ++j )
          stampVec[i][j] = gids[ stampVec[i][j] ];
      }

      std::vector<int> counts(3);
      counts[0] = extGIDs.size();
      counts[1] = intGIDs.size();
      counts[2] = depGIDs.size();

      indexor.matrixGlobalToLocal(Parallel::JACOBIAN_OVERLAP, gids, stampVec );

      cktNodeDevPtr->registerJacLIDswithDev( stampVec );
    }
  }

  indexor.deleteAcceleratedMatrixIndexing();
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::regenerateGIDNodeMap
// Purpose       : Redo global index node map
// Special Notes : Should find a way to automate this
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::regenerateGIDNodeMap()
{
  // Clear the map, just in case this was done before.
  indexToGID_.clear();
  gIDtoIndex_.clear();

  CktNodeList::const_iterator it_cnL = BFSNodeList_.begin();
  CktNodeList::const_iterator it_cnL_end = BFSNodeList_.end();
  for( ; it_cnL != it_cnL_end ; ++it_cnL )
  {
    int index = cktgph_.getIndex( (*it_cnL)->get_nodeID() );
    int gid = (*it_cnL)->get_gID();
    indexToGID_[ index ] = gid;
    gIDtoIndex_[ gid ] = index;
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::replaceNode
// Purpose       : replace supernode two nodes
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/10/10
//-----------------------------------------------------------------------------
CktNode * CktGraphBasic::replaceNode( const NodeID nodeToBeReplaced, 
                                      const NodeID nodeToKeep )
{
  // procedure is as follows:
  // 1. get the adjacency of both nodes
  // 2. combine to the two adjacencies
  // 3. set nodeToKeep's adjacency to what we found in step 2.
  // 4. tell the graph to replace any other adjacencies that refer to nodeToBeReplaced with nodeToKeep
  // 5. Remove nodeToBeReplaced and return it's associated CktNode object
  
  // look up key for these nodes in graph.
  CktNode * nodeToBeReplacedCktNodePtr = FindCktNode( nodeToBeReplaced );
 
  if (nodeToBeReplacedCktNodePtr)
  { 
    // 1. get the adjacency of both nodes
    std::vector<NodeID> adjNodeToBeReplaced;
    returnAdjIDs( nodeToBeReplaced, adjNodeToBeReplaced );
  
    // 2. set nodeToKeep's adjacency to include nodeToBeReplaced adjacencies.
    cktgph_.addToAdjacent( nodeToBeReplaced, nodeToKeep, adjNodeToBeReplaced );
 
    // 3. tell the graph to replace any other adjacencies that refer to nodeToBeReplaced with nodeToKeep
    cktgph_.replaceAdjacent( nodeToBeReplaced, nodeToKeep );

    // 4. Remove nodeToBeReplaced and return it's associated CktNode object
    cktgph_.removeKey( nodeToBeReplaced );

    // changed ordering so set isModified_ flag
    isModified_ = true;
  } 

  return nodeToBeReplacedCktNodePtr;
}


//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::removeRedundantDevices
// Purpose       : remove any device nodes that are only connected to one ckt node
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/10/10
//-----------------------------------------------------------------------------
void CktGraphBasic::removeRedundantDevices( std::vector< NodeID > & devicesToBeRemoved, 
                                            std::vector< CktNode * > & removedDevices)
{
  std::ostringstream outputStringStream;
  // Collect the deviceIDs for all the devices to be removed from the graph
  std::vector<NodeID> deviceIDs;
  const Graph::Data1Map &dataMap = cktgph_.getData1Map();

  Graph::Data1Map::const_iterator currentCktNodeItr = dataMap.begin();
  Graph::Data1Map::const_iterator endCktNodeItr = dataMap.end();
  while( currentCktNodeItr != endCktNodeItr )
  {
    if( ((*currentCktNodeItr).second)->type() == _DNODE )
    {
      CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>((*currentCktNodeItr).second);
      const Device::InstanceBlock *deviceInstanceBlockPtr = cktNodeDevPtr->deviceInstanceBlock();
      if (deviceInstanceBlockPtr)
      {
        // have a valid device.
        int deviceIndex = cktgph_.getIndex( (*currentCktNodeItr).first );
        const std::vector<int> &adjacentIDs = cktgph_.getAdjacentRow( deviceIndex );
        int numAdjIDs = adjacentIDs.size(); 
 
        bool removeDevice=false;
        if( numAdjIDs == 2 )
        {
          if( adjacentIDs[0] == adjacentIDs[1] )
            removeDevice=true;
        }
        else if( numAdjIDs == 3 )
        {
          if( ( adjacentIDs[0] == adjacentIDs[1] ) && (adjacentIDs[1] == adjacentIDs[2]) )
            removeDevice=true;
        }
        else if( numAdjIDs == 4 )
        {
          if( ( adjacentIDs[0] == adjacentIDs[1] ) && (adjacentIDs[1] == adjacentIDs[2])
            && (adjacentIDs[2] == adjacentIDs[3]) )
            removeDevice=true;
        }
        if( removeDevice )
        {
          deviceIDs.push_back((*currentCktNodeItr).first);
          devicesToBeRemoved.push_back((*currentCktNodeItr).first);
          removedDevices.push_back((*currentCktNodeItr).second);
        }
      }
      else
      {
        // issue fatal error as this case shouldn't occur
        Report::DevelFatal() << "Topology::removeRedundantDevices() null Device Instance Block pointer.";
      }
    }
    currentCktNodeItr++;
  }

  // If there are any devices to remove, remove them from the graph and mark the object as modified.
  if( removedDevices.size() > 0 )
  {
    cktgph_.removeKeys( deviceIDs );

    // After the devices are removed, check the graph to see if any singletons were created.
    // These are most likely ghost nodes that connected the removed devices to the rest of the circuit.
    // NOTE:  This can be a real issue in practice, so don't remove this search.
    std::vector< NodeID> singletonIDs = cktgph_.getSingletons();
    if (singletonIDs.size() > 0)
    {
      std::vector< NodeID> singletonDevs;

      cktgph_.removeKeys( singletonIDs );
   
      // Check if the singleton was a device
      std::vector< NodeID>::iterator sIDs_it = singletonIDs.begin();
      for ( ; sIDs_it != singletonIDs.end(); ++sIDs_it )
      {
        if( Xyce::get_node_type( *sIDs_it ) == _DNODE )
          singletonDevs.push_back( *sIDs_it );
      }

      // If there are singleton devices, find the device object to add to removedDevices
      std::vector< NodeID>::iterator sdIDs_it = singletonDevs.begin();
      for ( ; sdIDs_it != singletonDevs.end(); ++sdIDs_it )
      {
        Graph::Data1Map::const_iterator currentCktNodeItr = dataMap.begin();
        Graph::Data1Map::const_iterator endCktNodeItr = dataMap.end();
        while( currentCktNodeItr != endCktNodeItr )
        {
          if( (((*currentCktNodeItr).second)->type() == _DNODE) && ((*currentCktNodeItr).first == *sdIDs_it) )
          {
            devicesToBeRemoved.push_back((*currentCktNodeItr).first);
            removedDevices.push_back((*currentCktNodeItr).second);
          }
          currentCktNodeItr++;
        }
      }
    }

    isModified_=true;
  }

}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::put
// Purpose       : Allows virtual override of operator<<
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::removeNodes( const std::vector< NodeID > & nodesToBeRemoved, 
                                 std::vector< CktNode * > & removedNodes )
{ 
  std::vector< NodeID >::const_iterator it = nodesToBeRemoved.begin();
  for ( ; it != nodesToBeRemoved.end(); ++it )
  {
    CktNode * nodeToBeReplacedCktNodePtr = FindCktNode(*it);
    if (nodeToBeReplacedCktNodePtr)
      removedNodes.push_back( nodeToBeReplacedCktNodePtr );
  }  
 
  cktgph_.removeKeys( nodesToBeRemoved ); 
  isModified_=true; 
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::put
// Purpose       : Allows virtual override of operator<<
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
std::ostream& CktGraphBasic::put(std::ostream& os) const
{
  int i = 0;

  CktNodeList::const_iterator it_cnL = BFSNodeList_.begin();
  CktNodeList::const_iterator it_cnL_end = BFSNodeList_.end();
  for( ; it_cnL != it_cnL_end ; ++it_cnL, ++i )
  {
    os << "[" << i << "]:" << **it_cnL << std::endl;
  }

  // Print out information about the circuit graph.
  cktgph_.print( os );

  return os;
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::analyzeDeviceNodeGraph
// Purpose       : Find floating nodes in device node graph for removal or rebalance
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 7/21/2021
//-----------------------------------------------------------------------------
std::vector< Xyce::NodeID > CktGraphBasic::analyzeDeviceNodeGraph(std::ostream & os) 
{
  // The device node graph is the distance-2 graph
  const CktGraph::Graph::Index1Map& indexMap = cktgph_.getIndex1Map();
  
  // Get the local ID for the ground node, which we will ignore when generating the device graph.
  Xyce::NodeID gnd( "0", _VNODE );
  CktGraph::Graph::Index1Map::const_iterator currentIndexItr = indexMap.find( gnd );
  int gndIdx = -1;
  if (currentIndexItr != indexMap.end())
    gndIdx = currentIndexItr->second; 
 
  unordered_map<int, std::vector<int> > deviceAdj;

  int maxIndex = 0; 
  currentIndexItr = indexMap.begin();
  CktGraph::Graph::Index1Map::const_iterator endIndexItr = indexMap.end();
  while( currentIndexItr != endIndexItr )
  { 
    if( Xyce::get_node_type((*currentIndexItr).first) == _DNODE )
    { 
      // Get the voltage nodes adjacent to this device
      const std::vector<int>& vAdj = cktgph_.getAdjacentRow( currentIndexItr->second );

      // Loop over the voltage nodes and get their adjacencies, which are devices
      std::vector<int>::const_iterator vIdx = vAdj.begin();
      for ( ; vIdx != vAdj.end(); ++vIdx )
      {
        if (*vIdx != gndIdx)
        {
          const std::vector<int>& dAdj = cktgph_.getAdjacentRow( *vIdx ); 
          deviceAdj[ currentIndexItr->second ].insert( deviceAdj[ currentIndexItr->second ].end(),
                                                       dAdj.begin(), dAdj.end() );
        }
      }

      // Get the maximum device index
      if ( currentIndexItr->second > maxIndex )
        maxIndex = currentIndexItr->second;
    }
    currentIndexItr++;
  } 

  int totalEntries = 0; 
  std::vector< Xyce::NodeID > floatingDevs;

  os << "-------------------- Device Graph ----------------------------\n";
  os << deviceAdj.size() << " " << maxIndex << std::endl;
  unordered_map<int, std::vector<int> >::iterator dIdx = deviceAdj.begin();
  for( ; dIdx != deviceAdj.end(); ++dIdx )
  {
    // Make sure the entries for this device are unique
    std::sort( (*dIdx).second.begin(), (*dIdx).second.end() );
    (*dIdx).second.erase(std::unique((*dIdx).second.begin(), (*dIdx).second.end()), (*dIdx).second.end());

    // Write out the device id and its connected devices to the output stream
    os << "[ " << dIdx->first << ", " << (cktgph_.getKey( dIdx->first )).first << " ] : ";
    for( size_t j = 0; j < (*dIdx).second.size(); ++j )
    {
      if (((*dIdx).second)[j] != dIdx->first)
        os << " " << ((*dIdx).second)[j]; 
    }
    os << std::endl;
    if ((*dIdx).second.size() == 1)
      floatingDevs.push_back( cktgph_.getKey( dIdx->first ) );
  }

  return floatingDevs;
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::streamCircuitGraph
// Purpose       : Send the circuit graph to an output stream
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 7/21/2021
//-----------------------------------------------------------------------------
void CktGraphBasic::streamCircuitGraph(std::ostream & os)
{
  // The device node graph is the distance-2 graph
  const CktGraph::Graph::Index1Map& indexMap = cktgph_.getIndex1Map();

  CktGraph::Graph::Index1Map::const_iterator currentIndexItr = indexMap.begin();
  CktGraph::Graph::Index1Map::const_iterator endIndexItr = indexMap.end();

  os << "-------------------- Circuit Graph ----------------------------\n";
  os << indexMap.size() << std::endl;

  while( currentIndexItr != endIndexItr )
  { 
    // Get the voltage nodes adjacent to this device
    const std::vector<int>& vAdj = cktgph_.getAdjacentRow( currentIndexItr->second );

    // Loop over the voltage nodes and get their adjacencies, which are devices
    std::vector<int>::const_iterator vIdx = vAdj.begin();

    // Write out the device id and its connected devices to the output stream
    os << "[ " << currentIndexItr->second <<  ", " << currentIndexItr->first << " ] : ";
    for( ; vIdx != vAdj.end(); ++vIdx )
    {
        os << " " << *vIdx;
    }
    os << std::endl;

    currentIndexItr++;
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::removeFloatingNodes
// Purpose       : Find unattached device or voltage nodes in graph for removal.
// Special Notes : This cleanup is necessary after all graph analyses have been performed.
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 8/3/2021
//-----------------------------------------------------------------------------
void CktGraphBasic::removeUnattachedNodes()
{
  const CktGraph::Graph::Index1Map& indexMap = cktgph_.getIndex1Map();

  CktGraph::Graph::Index1Map::const_iterator currentIndexItr = indexMap.begin();
  CktGraph::Graph::Index1Map::const_iterator endIndexItr = indexMap.end();
  std::vector< Xyce::NodeID > nodesToBeRemoved;
  std::vector< CktNode * > removedNodes;
 
  while( currentIndexItr != endIndexItr )
  {
    // Get the nodes adjacent to this node
    const std::vector<int>& vAdj = cktgph_.getAdjacentRow( currentIndexItr->second );
    if (!vAdj.size())
    {
      nodesToBeRemoved.push_back( currentIndexItr->first );
  
      CktNode * nodeToBeReplacedCktNodePtr = FindCktNode( currentIndexItr->first );
      if (nodeToBeReplacedCktNodePtr)
        removedNodes.push_back( nodeToBeReplacedCktNodePtr );
    }       
    currentIndexItr++;
  }    

  Xyce::dout() << "removeUnattachedNodes found " << nodesToBeRemoved.size() << " nodes to remove!" << std::endl;
  cktgph_.removeKeys( nodesToBeRemoved );

  // now it's safe to delete the device nodes
  for (std::vector< CktNode * >::iterator it = removedNodes.begin(); it != removedNodes.end(); ++it)
    delete *it;  

  isModified_=true;
}


} // namespace Topo
} // namespace Xyce

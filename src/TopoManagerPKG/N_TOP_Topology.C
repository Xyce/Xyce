//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
// Creation Date  : 03/20/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <sstream>
#include <fstream>

#include <N_UTL_fwd.h>

#include <N_ANP_AnalysisManager.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_DeviceState.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CmdParse.h>
#include <N_IO_HangingResistor.h>
#include <N_IO_RestartNode.h>
#include <N_PDS_Comm.h>
#include <N_PDS_Manager.h>
#include <N_TOP_CktGraph.h>
#include <N_TOP_CktGraphBasic.h>
#include <N_TOP_CktNode.h>
#include <N_TOP_CktNode_Dev.h>
#include <N_TOP_CktNode_V.h>
#include <N_TOP_Directory.h>
#include <N_TOP_Indexor.h>
#include <N_TOP_NodeDevBlock.h>
#include <N_TOP_LSUtilFactory.h>
#include <N_TOP_Topology.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Functors.h>
#include <N_UTL_Marshal.h>
#include <N_UTL_OptionBlock.h>
#include <N_LAS_QueryUtil.h>

#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>


namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : Topology::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool registerPkgOptionsMgr(Topology &topology, IO::PkgOptionsMgr &options_manager)
{
  registerPkgOptionsMgr(*topology.getLinearSolverUtility(), options_manager);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Topology::Topology
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
Topology::Topology(
  const IO::CmdParse &          cp,
  IO::HangingResistor &         hanging_resistor,
  N_PDS_Manager &               pds_manager)
  : commandLine_(cp),
    hangingResistor_(hanging_resistor),
    linearSolverUtility_(Topo::LSUtilFactory::newLSUtil(*this, pds_manager)),
    pdsManager_(pds_manager),
    mainGraphPtr_( new CktGraphBasic() ),
    nodeSymbols_(Util::SYMBOL_TYPE_END),
    gnd_("gnd")
{}

//-----------------------------------------------------------------------------
// Function      : Topology::~Topology
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
Topology::~Topology()
{
  delete linearSolverUtility_;
  delete mainGraphPtr_;
}

//-----------------------------------------------------------------------------
// Returns adj IDs to the given ID
//-----------------------------------------------------------------------------
void Topology::returnAdjIDs( const NodeID& id, std::vector<NodeID>& adj_ids ) const
{
  return mainGraphPtr_->returnAdjIDs(id, adj_ids);
}

//-----------------------------------------------------------------------------
// Returns adj GIDs to the given GID, including the ground node (GID of -1)
//-----------------------------------------------------------------------------
void Topology::returnAdjGIDsWithGround( int gid, std::vector<int>& adj_gids ) const
{
  return mainGraphPtr_->returnAdjGIDsWithGround(gid, adj_gids);
}

//-----------------------------------------------------------------------------
// Returns number of adj nodes
//-----------------------------------------------------------------------------
int Topology::numAdjNodes( int gid ) const
{
  return mainGraphPtr_->numAdjNodes( gid );
}

//-----------------------------------------------------------------------------
// Returns number of adj nodes, including the ground node
//-----------------------------------------------------------------------------
int Topology::numAdjNodesWithGround( int gid ) const
{
  return mainGraphPtr_->numAdjNodesWithGround( gid );
}

//-----------------------------------------------------------------------------
// Returns pointer to specified ckt device node.
//-----------------------------------------------------------------------------
const CktNode * Topology::findCktNode(const NodeID& cnID) const
{
  return mainGraphPtr_->FindCktNode(cnID);
}

//-----------------------------------------------------------------------------
// Function      : Topology::setupGlobalIndices
// Purpose       : Lin system data setup
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/1/01
//-----------------------------------------------------------------------------
bool Topology::setupGlobalIndices()
{
  return linearSolverUtility_->setupRowCol();
}

//-----------------------------------------------------------------------------
// Function      : Topology::addDevice
// Purpose       : New dev-node instantiator (planarized ckts only)
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/20/00
//-----------------------------------------------------------------------------
void
Topology::addDevice(
  Device::DeviceMgr &           device_manager,
  const NodeDevBlock &             node_block)
{
  std::vector<NodeID> emptyNLList, nlList;
  //------ Add connected voltage nodes

  std::vector<std::string>::const_iterator it_nlL = node_block.get_NodeList().begin();
  std::vector<std::string>::const_iterator end_nlL = node_block.get_NodeList().end();

  for( ; it_nlL != end_nlL ; ++it_nlL )
  {
    //----- insert each v-node: Unless the v-node already exists,
    //----- it will be instantiated but not owned and ProcNum_
    //----- identifies the owning processor
    mainGraphPtr_->InsertNode(new CktNode_V(*it_nlL), emptyNLList);

    nlList.push_back( NodeID( *it_nlL, _VNODE ) );
  }

  //------- Now instantiate the device node,  DistribMgr has
  //------- already set ownership

  mainGraphPtr_->InsertNode(new CktNode_Dev(node_block.getDevBlock(), device_manager), nlList);
}

//-----------------------------------------------------------------------------
// Function      : Topology::registerGIDs
// Purpose       : register the ext. global ids stored by
//                 the cktnodes with their respective devices
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/22/00
//-----------------------------------------------------------------------------
void Topology::registerGIDs()
{
  mainGraphPtr_->registerGIDs();
}

//-----------------------------------------------------------------------------
// Function      : Topology::registerLIDswithDevs
// Purpose       : register the int. and ext. local ids stored by
//                 the cktnodes with their respective devices
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/12/02
//-----------------------------------------------------------------------------
void Topology::registerLIDswithDevs()
{
  Indexor indexor( pdsManager_ );

  mainGraphPtr_->registerLIDswithDevs( indexor );
  mainGraphPtr_->registerStateLIDswithDevs( indexor );
  mainGraphPtr_->registerStoreLIDswithDevs( indexor );
  mainGraphPtr_->registerBranchDataLIDswithDevs( indexor );

  mainGraphPtr_->registerDepLIDswithDevs( indexor );
  mainGraphPtr_->registerJacLIDswithDevs( indexor );
}

//-----------------------------------------------------------------------------
// Function      : Topology::resolveDependentVars
// Purpose       : loop through devices and resolve their secondary
//                 dependencies
// Special Notes : Used to resolve Expressions and Current Dependencies
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/05/01
//-----------------------------------------------------------------------------
void Topology::resolveDependentVars()
{
  int solnCount = 0;
  std::vector<int> solnLocVec;
  std::vector<NodeID> solnNidVec;
  std::vector<NodeID> idVec;

  // Find out if any device nodes have dependent variables that need to be resolved.
  for (CktNodeList::iterator it_cnL = mainGraphPtr_->getBFSNodeList()->begin(),
      it_cnL_end = mainGraphPtr_->getBFSNodeList()->end(); it_cnL != it_cnL_end; ++it_cnL )
  {
    if ( (*it_cnL)->type() == _DNODE )
    {
      CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_cnL);

      // Get dependent solution variables
      cktNodeDevPtr->getDepSolnVars( idVec );
      if( !idVec.empty() )
      {
        solnLocVec.push_back( solnCount );
        solnNidVec.insert( solnNidVec.end(), idVec.begin(), idVec.end() );
        solnCount += idVec.size();
      }
      else
        solnLocVec.push_back( solnCount );
    }
  }
  solnLocVec.push_back( solnCount );

  // Determine if we need to generate the parallel directory by summing all counts.
  int tSolnCount;
  pdsManager_.getPDSComm()->sumAll( &solnCount, &tSolnCount, 1 );

  // If there is work to do, then create the directory and resolve the GIDs.
  if (tSolnCount)
  {
    // Generate the directory.
    Directory directory(*this, *pdsManager_.getPDSComm());
    directory.generateDirectory();

    std::vector< std::vector<int> > gidVec, indexVec;
    std::vector<int> procVec;

    directory.getSolnGIDs( solnNidVec, gidVec, procVec );

    int count = 0;
    for (CktNodeList::iterator it_cnL = mainGraphPtr_->getBFSNodeList()->begin(),
         it_cnL_end = mainGraphPtr_->getBFSNodeList()->end(); it_cnL != it_cnL_end; ++it_cnL )
    {
      if ( (*it_cnL)->type() == _DNODE )
      {
        CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_cnL);
        indexVec.assign( gidVec.begin()+solnLocVec[count],
                         gidVec.begin()+solnLocVec[count+1] );
        if (!indexVec.empty()) 
        {
          cktNodeDevPtr->set_DepSolnGIDVec( indexVec );
          cktNodeDevPtr->registerDepSolnGIDs( indexVec );
        }
        ++count;
      }
    }

    if (!pdsManager_.getPDSComm()->isSerial())
    {
      depSolnGIDMap_.clear();
      for( unsigned int i = 0; i < gidVec.size(); ++i )
        for( unsigned int ii = 0; ii < gidVec[i].size(); ++ii )
          depSolnGIDMap_[ gidVec[i][ii] ] = procVec[i];
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::OutputBFSGraphLists
// Purpose       : Output to Xyce::dout() BFS node list for debugging
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
void Topology::OutputBFSGraphLists() const
{
  Xyce::dout() << "BFS Node Listing for Graphs" << std::endl;

  for (CktNodeList::const_iterator it_cnL = mainGraphPtr_->getBFSNodeList()->begin(),
      end_cnL = mainGraphPtr_->getBFSNodeList()->end(); it_cnL != end_cnL; ++it_cnL )
  {
    Xyce::dout() << *(*it_cnL) << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::generateOrderedNodeList
// Purpose       : Currently sets orderedNodeListPtr_ attribute to
//                 BFS traversal of main ckt.
// Special Notes :
// Scope         : public
// Creator       : Rob  Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
void Topology::generateOrderedNodeList() const
{
  orderedNodeListPtr_ = mainGraphPtr_->getBFSNodeList();
}

//-----------------------------------------------------------------------------
// Function      : verifyNodesAndDevices -- scans existing topology
//                 and asks device manager to verify each device.
//                 if the device manager returns false on any verification,
//                 then the nodes on that device are added to a list of nodes
//                 to be supernoded and the device will later be removed as redundant
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems modeling
// Creation Date : 2/18/2010
//-----------------------------------------------------------------------------
void Topology::verifyNodesAndDevices(
  Device::DeviceMgr &   device_manager)
{
  if( linearSolverUtility_->supernodeFlag() )
  {
    int badDeviceCount=0;

    const CktGraph::Graph::Data1Map & dataMap = mainGraphPtr_->getNodeList();
    CktGraph::Graph::Data1Map::const_iterator currentCktNodeItr = dataMap.begin();
    CktGraph::Graph::Data1Map::const_iterator endCktNodeItr = dataMap.end();
    while( currentCktNodeItr != endCktNodeItr )
    {
      if( ((*currentCktNodeItr).second)->type() == _DNODE )
      {
        CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>((*currentCktNodeItr).second);
        const Device::InstanceBlock *deviceInstanceBlockPtr = cktNodeDevPtr->deviceInstanceBlock();
        if (deviceInstanceBlockPtr)
        {
          bool deviceInstanceOk = device_manager.verifyDeviceInstance( *deviceInstanceBlockPtr );
          if( !deviceInstanceOk )
          {
            badDeviceCount++;
            // place the nodes for this device in the superNodeList_
            // collected nodes so that they can be combined (supernoded)
            // It's ok to let the device node get created and inserted into the
            // topology.  we can remove it later
            const NodeID& deviceID = (*currentCktNodeItr).second->get_nodeID();
            badDeviceList_.push_back( deviceID );

            std::vector< NodeID > adjacentIDs;
            mainGraphPtr_->returnAdjIDs( deviceID, adjacentIDs );

            std::vector< NodeID >::iterator currentID = adjacentIDs.begin();
            std::vector< NodeID >::iterator endID = adjacentIDs.end();
            std::vector< NodeID >::iterator nextID = currentID;
            if (currentID != endID)
              nextID++;

            while( nextID != endID )
            {
              
              if ( (*currentID) != (*nextID) )
              {
                // these nodes are to be supernoded.  Take the lexically smaller one
                // as the node to keep (i.e. A < B )
                if( (*currentID) < (*nextID) )
                {
                  // format of pair is nodeToReplace, nodeToKeep
                  superNodeList_.push_back( make_pair( *nextID, *currentID ) );
                }
                else
                {
                  // format of pair is nodeToReplace, nodeToKeep
                  superNodeList_.push_back( make_pair( *currentID, *nextID ) );
                }
             }
              currentID++;
              nextID++;
            }
          }
        }
        else
        {
          // issue fatal error as this case shouldn't occur
          Report::DevelFatal().in("Topology::verifyNodesAndDevices") 
            << "null Device Instance Block pointer";
        }
      }
      currentCktNodeItr++;
    }

    N_PDS_Comm * commPtr = pdsManager_.getPDSComm();
    int totalBadDevices = 0;
    commPtr->sumAll( &badDeviceCount, &totalBadDevices, 1 );

    Report::UserWarning0() << "Device verification found " << totalBadDevices << " device(s) to remove";
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::removeTaggedNodesAndDevices
// Purpose       : Remove devices and nodes that were tagged for removal
//                 during parsing.  Node removal is done through supernoding,
//                 .i.e. replacing one node globally with another.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems modeling
// Creation Date : 2/2/2010
//-----------------------------------------------------------------------------
void Topology::removeTaggedNodesAndDevices()
{
  if( linearSolverUtility_->supernodeFlag() )
  {
    if (DEBUG_TOPOLOGY)
      Xyce::dout() << "Topology::removeTaggedNodesAndDevices" << std::endl
                   << *this << std::endl;

    // First remove bad devices for this processor
    std::vector< CktNode * > removedDevices;
    mainGraphPtr_->removeNodes( badDeviceList_, removedDevices );

    // storage for oldNode that we'll delete when done with this routine
    // use a set because we can get the same oldNode CktNode pointer
    // multiple times
    CktNodeList oldNodeList;

    // print out current state of supernode list
    std::vector< std::pair<NodeID, NodeID> >::iterator currentNodePair = superNodeList_.begin();
    std::vector< std::pair<NodeID, NodeID> >::iterator endNodePair = superNodeList_.end();
    unordered_set<NodeID> nodesReplaced;
    while ( currentNodePair != endNodePair )
    {
      NodeID nodeToBeReplaced( currentNodePair->first );
      NodeID replacementNode( currentNodePair->second );
      if( nodeToBeReplaced != replacementNode)
      {
        unordered_set<NodeID>::iterator nodesReplacedEnd = nodesReplaced.end();
        unordered_set<NodeID>::iterator nodeReplacedLoc = nodesReplaced.find( nodeToBeReplaced );
        if( nodeReplacedLoc == nodesReplacedEnd )
        {
          if (DEBUG_TOPOLOGY && nodeToBeReplaced < replacementNode)
            Xyce::dout() << "Ordering is wrong on nodes!" << std::endl;

          // the replacement node might be from another processor, so it is not in this processors' graph; add it now
          CktNode * newNode = mainGraphPtr_->FindCktNode( replacementNode );
          if ( !newNode )         
          {
            std::vector<NodeID> emptyNLList;
            newNode = new CktNode_V( Xyce::get_node_id( replacementNode ) );
            newNode->set_IsOwned( false );
            mainGraphPtr_->InsertNode( newNode, emptyNLList );
          }

          CktNode * oldNode = mainGraphPtr_->replaceNode( nodeToBeReplaced, replacementNode );

          // there is a possibility that this node may not be found, so ignore it and move on.
          if ( oldNode )
          {
            nodesReplaced.insert( nodeToBeReplaced );

            // If we delete this old node now, we'll make the orderedNodeListPtr_ untraversable.
            // this would be ok as we can regenerate it, but that takes time and we would only
            // invalidate it when we do our next delete.  So, store up the oldNode so we can
            // delete them when were done
            oldNodeList.push_back( oldNode );

            // now that we've replaced "nodeToBeReplaced" with "replacementNode" we need to
            // search the superNodeList_ from this point on also doing this same substitution
            //
            // for example if our super node list was:
            //
            //    B   A
            //    C   B
            //
            // If after the B->A substitution we didn't update our list, then we would next bring
            // back the B's with C->B.
            //
            std::vector< std::pair<NodeID, NodeID> >::iterator nextNodePair = currentNodePair;
            nextNodePair++;
            while ( nextNodePair != endNodePair )
            {
              if( nextNodePair->first == nodeToBeReplaced )
              {
                if( replacementNode < nextNodePair->second )
                {
                  // need to swap on insert
                  *nextNodePair = make_pair( nextNodePair->second, replacementNode );
                }
                else
                {
                  // just insert in order, new pair is ( replacementNode, nextNodePair->second )
                  nextNodePair->first = replacementNode;
                }
              }
              else if( nextNodePair->second == nodeToBeReplaced )
              {
                if( replacementNode < nextNodePair->first )
                {
                  // no swap needed, new pair is ( nextNodePair->first, replacementNode )
                  nextNodePair->second = replacementNode;
                }
                else
                {
                  // need to swap
                  *nextNodePair = make_pair( replacementNode, nextNodePair->first );
                }
              }
              nextNodePair++;
            }
          }
        }
      }
      currentNodePair++;
    }

    {
      //Stats::StatTop _topoStat("Topology Remove Redundant Devices");
      //Stats::TimeBlock _topoTimer(_topoStat);
      mainGraphPtr_->removeRedundantDevices(removedDevices);
    }

    // now it's safe to delete the old nodes and device nodes
    for (CktNodeList::iterator it = oldNodeList.begin(), end = oldNodeList.end(); it != end; ++it)
    {
      delete *it;
    }

    if (!removedDevices.empty())
    {
      // delete old devices that were removed
      for (std::vector<CktNode *>::iterator it = removedDevices.begin(), end = removedDevices.end(); it != end; ++it)
      {
        CktNode_Dev *cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it);

        if (cktNodeDevPtr && cktNodeDevPtr->deviceInstanceBlock())
        {
          delete cktNodeDevPtr;
        }
      }
    }

    // Clear supernode storage.
    superNodeList_.clear();
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::mergeOffProcTaggedNodesAndDevices
// Purpose       : Merge the off processor superNodeList_ and communicate
//                 the same list to all procs so topology reduction is the
//                 same on all procs
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems modeling
// Creation Date : 2/2/2010
//-----------------------------------------------------------------------------
void Topology::mergeOffProcTaggedNodesAndDevices()
{
  if( linearSolverUtility_->supernodeFlag() && !(pdsManager_.getPDSComm())->isSerial() )
  {
    N_PDS_Comm * commPtr = pdsManager_.getPDSComm();
    int numProcs = commPtr->numProc();
    int thisProc = commPtr->procID();

    int localNodes = superNodeList_.size();
    std::vector< std::pair<NodeID, NodeID> > externalSuperNodeList;

    // Count the bytes for packing these strings
    int byteCount = 0;

    // First count the size of the local superNodeList_
    byteCount += sizeof(int);

    // Now count all the NodeIDs in the list
    unordered_set<NodeID> nodesReplaced;
    for (std::vector< std::pair<NodeID, NodeID> >::iterator nodePair = superNodeList_.begin(); 
        nodePair != superNodeList_.end(); nodePair++) 
    {
      byteCount += Xyce::packedByteCount(nodePair->first);
      byteCount += Xyce::packedByteCount(nodePair->second);
      nodesReplaced.insert( nodePair->first );
    }

    for( int p = 0; p < numProcs; ++p )
    {
      commPtr->barrier();

      // Broadcast the buffer size for this processor.
      int bsize=0;
      if (p==thisProc) { bsize = byteCount; }
      commPtr->bcast( &bsize, 1, p );

      // Create buffer.
      int pos = 0;
      char * superNodeBuffer = new char[bsize];

      if (p==thisProc) 
      {
        // Pack number of supernodes on this processor and place in externalSuperNodeList
        commPtr->pack( &localNodes, 1, superNodeBuffer, bsize, pos );
        for (std::vector< std::pair<NodeID, NodeID> >::iterator nodePair = superNodeList_.begin(); 
            nodePair != superNodeList_.end(); nodePair++) 
        {
          // Pack first entry of pair.
          Xyce::pack(nodePair->first, superNodeBuffer, bsize, pos, commPtr);

          // Pack second entry of pair.
          Xyce::pack(nodePair->second, superNodeBuffer, bsize, pos, commPtr);
        }

        // Broadcast packed buffer.
        commPtr->bcast( superNodeBuffer, bsize, p );
      }
      else 
      {
        // Unpack buffer and place in externalSuperNodeList
        commPtr->bcast( superNodeBuffer, bsize, p );

        // Get number of supernodes from that processor
        int numSuperNodes = 0;
        commPtr->unpack( superNodeBuffer, bsize, pos, &numSuperNodes, 1 );

        // Extract supernode pairs and push to the back of the externalSuperNodeList
        std::pair<NodeID,NodeID> tmpPair;

        for (int i=0; i<numSuperNodes; ++i) {

          // Unpack first pair.
          Xyce::unpack(tmpPair.first, superNodeBuffer, bsize, pos, commPtr); 

          // Unpack second pair.
          Xyce::unpack(tmpPair.second, superNodeBuffer, bsize, pos, commPtr); 

          // Push to back of externalSuperNodeList
          externalSuperNodeList.push_back( tmpPair );
        }
      }
      // Clean up.
      delete [] superNodeBuffer;
    }


    // Now go through the local superNodeList_ and tack on any boundary cases, where adjacencies may effect device removal.
    int i = 0, nodeListSize = superNodeList_.size();
    while( i < nodeListSize )
    {
      NodeID nodeToBeReplaced = superNodeList_[i].first;
      NodeID replacementNode = superNodeList_[i].second;
      std::vector< std::pair<NodeID, NodeID> >::iterator currentGlobalSN = externalSuperNodeList.begin();
      std::vector< std::pair<NodeID, NodeID> >::iterator endGlobalSN = externalSuperNodeList.end();
      while ( currentGlobalSN != endGlobalSN )
      {
        // Add pair if replacement node is node to be removed by another processor.
        if (currentGlobalSN->first == replacementNode)
        {
          superNodeList_.push_back( *currentGlobalSN );
          nodesReplaced.insert( currentGlobalSN->first );
          nodeListSize++;
        }
        // Add pair if node to be removed is a replacement node on another processor.
        else if (currentGlobalSN->second == nodeToBeReplaced)
        {  
          superNodeList_.push_back( *currentGlobalSN );
          nodesReplaced.insert( currentGlobalSN->first );
          nodeListSize++;
        }
        currentGlobalSN++;
      }
      i++;
    }

    // Condensing external super node list so that any supernodes found in the next step do not need iteration.
    // Ex.  Processor 2 may have a supernode ( A, 0 ) to ( B, 0 ) that hasn't been accounted for, but Processor 3 has ( B, 0 ) to ( C, 0 ).
    //      This condensing should mean that ( A, 0 ) to ( C, 0 ) is the supernode that is picked out below.
    //      Need to consider all global supernodes because they could have been collected in any order, keep iterating over them
    //      until nothing is modified on this processor.
    int numModified=1;
    while (numModified > 0)
    {
      std::vector< std::pair<NodeID, NodeID> >::iterator currentGlobalSN = externalSuperNodeList.begin();
      std::vector< std::pair<NodeID, NodeID> >::iterator endGlobalSN = externalSuperNodeList.end();
      numModified = 0;
      for ( ; currentGlobalSN != endGlobalSN; ++currentGlobalSN )
      {
        NodeID nodeToBeReplaced = currentGlobalSN->first;
        NodeID replacementNode = currentGlobalSN->second;
        std::vector< std::pair<NodeID, NodeID> >::iterator nodePair = externalSuperNodeList.begin();
        for ( ; nodePair != endGlobalSN; ++nodePair )
        {
          if( (nodePair != currentGlobalSN) )
          {
            if( nodePair->first == replacementNode )
            {
              currentGlobalSN->second = nodePair->second;
              replacementNode = nodePair->second;
              numModified++;
            }
            else if( nodePair->second == nodeToBeReplaced )
            {
              nodePair->second = replacementNode;
              numModified++;
            }
          }
        }
      }
    } 

    // There might be a chance that the node to be removed on another processor is a node on this processor, but not a supernode.
    std::vector< std::pair<NodeID, NodeID> >::iterator currentGlobalSN = externalSuperNodeList.begin();
    std::vector< std::pair<NodeID, NodeID> >::iterator endGlobalSN = externalSuperNodeList.end();
    for ( ; currentGlobalSN != endGlobalSN; ++currentGlobalSN )
      if ( (nodesReplaced.find( currentGlobalSN->first ) == nodesReplaced.end()) && mainGraphPtr_->FindCktNode( currentGlobalSN->first ) )
        superNodeList_.push_back( *currentGlobalSN );

  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::instantiateDevices
// Purpose       : Delayed instantiation of devices
// Special Notes :
// Scope         : public
// Creator       : Rob  Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/03
//-----------------------------------------------------------------------------
void Topology::instantiateDevices()
{
  generateOrderedNodeList();

  for (CktNodeList::iterator it = orderedNodeListPtr_->begin(), 
      end = orderedNodeListPtr_->end(); it != end; ++it )
  {
    if( (*it)->type() == _DNODE )
      (dynamic_cast<CktNode_Dev*>(*it))->instantiate();
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::getNodeSVarGIDs
// Purpose       : Return list of solution var indices for named node.
// Special Notes : returns false if node not owned or not local.
// Scope         : public
// Creator       : Rob  Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/11/00
//-----------------------------------------------------------------------------
bool Topology::getNodeSVarGIDs( const NodeID& id,
		                std::vector<int> & sVarGIDList,
		                std::vector<int> & extSVarGIDList,
                                char & type ) const
{
  CktNode * cnPtr = mainGraphPtr_->FindCktNode( id );

  if( cnPtr != NULL )
  {
    sVarGIDList.assign( cnPtr->get_SolnVarGIDList().begin(), cnPtr->get_SolnVarGIDList().end() );

    if( cnPtr->type() == _DNODE )
    {
      type = 'D';
      CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(cnPtr);
      extSVarGIDList.assign( cktNodeDevPtr->get_ExtSolnVarGIDList().begin(), cktNodeDevPtr->get_ExtSolnVarGIDList().end() );
    }
    else
    { 
      type = 'V';
    }

    if( cnPtr->get_IsOwned() )
    {
      return true;
    }
    else
    {
      sVarGIDList.clear();
      return false;
    }
  }
  else
    return false;
}

//-----------------------------------------------------------------------------
// Function      : Topology::regenerateGIDNodeMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/6/01
//-----------------------------------------------------------------------------
void Topology::regenerateGIDNodeMap()
{
  mainGraphPtr_->regenerateGIDNodeMap();
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : implementation for debugging
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const Topology &topology)
{
  os << topology.getMainGraph() << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Topology::getRestartNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/29/01
//-----------------------------------------------------------------------------
bool Topology::getRestartNodes(Analysis::AnalysisManager &analysis_manager, std::vector<IO::RestartNode*> & nodeVec )
{
  if (!orderedNodeListPtr_)
    return false;

  int count = 0;
  for (CktNodeList::const_iterator iterCN = orderedNodeListPtr_->begin(), 
      endCN  = orderedNodeListPtr_->end(); iterCN != endCN; ++iterCN )
  {
    if( (*iterCN)->get_IsOwned() && (*iterCN)->get_gID() != -1 )
    {
      ++count;
    }
  }

  nodeVec.resize(count);

  count = 0;
  for (CktNodeList::iterator iterCN = orderedNodeListPtr_->begin(), 
      endCN  = orderedNodeListPtr_->end(); iterCN != endCN; ++iterCN )
  {
    if( (*iterCN)->get_IsOwned() && (*iterCN)->get_gID() != -1 )
    {
      CktNode *cnP = *iterCN;

      nodeVec[count] = new IO::RestartNode( cnP->get_id(), cnP->type() );

      int i = 0;
      const std::vector<int> & gidList = cnP->get_SolnVarGIDList();
      nodeVec[count]->solnVarData.resize( gidList.size() );
      for( std::vector<int>::const_iterator iterIC = gidList.begin(); iterIC != gidList.end(); ++iterIC, ++i )
        analysis_manager.getSolnVarData( *iterIC, nodeVec[count]->solnVarData[i] );

      if( cnP->type() == _DNODE )
      {
        CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(cnP);

        i = 0;
        int numStateVars = cktNodeDevPtr->stateVarCount();
        nodeVec[count]->stateVarData.resize( numStateVars );
        if (numStateVars)
        {
          const std::vector<int> & gidList2 = cktNodeDevPtr->get_StateVarGIDList();
          for( std::vector<int>::const_iterator iterIC = gidList2.begin(); iterIC != gidList2.end(); ++iterIC, ++i )
            analysis_manager.getStateVarData( *iterIC, nodeVec[count]->stateVarData[i] );
        }

        i = 0;
        int numStoreVars = cktNodeDevPtr->storeVarCount();
        nodeVec[count]->storeVarData.resize( numStoreVars );
        if (numStoreVars)
        {
          const std::vector<int> & gidList3 = cktNodeDevPtr->get_StoreVarGIDList();
          for( std::vector<int>::const_iterator iterIC = gidList3.begin(); iterIC != gidList3.end(); ++iterIC, ++i )
            analysis_manager.getStoreVarData( *iterIC, nodeVec[count]->storeVarData[i] );
        }

        nodeVec[count]->devState = cktNodeDevPtr->getDevState();
      }
      ++count;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Topology::restoreRestartNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/29/01
//-----------------------------------------------------------------------------
bool Topology::restoreRestartNodes(
    Analysis::AnalysisManager &analysis_manager, 
    const std::vector<IO::RestartNode*> & nodeVec)
{
  CktNode * cnP;

  for( unsigned int i = 0; i < nodeVec.size(); ++i )
  {
    cnP = mainGraphPtr_->FindCktNode( NodeID(nodeVec[i]->ID,nodeVec[i]->type) );

    if( cnP != NULL )
    {
      if (DEBUG_RESTART)
      {
        Xyce::dout() << "Restoring Node: " << nodeVec[i]->ID << std::endl;
      }

      const std::vector<int> & gidList = cnP->get_SolnVarGIDList();
      int pos = 0;
      // Check if the new solution variable data is the same size as the old.
      // Do not set old solution variable data if they do not match, model has changed.
      if (gidList.size() == nodeVec[i]->solnVarData.size())
      {
        for( std::vector<int>::const_iterator iterIC = gidList.begin();
             iterIC != gidList.end(); ++iterIC, ++pos )
        {
          analysis_manager.setSolnVarData( *iterIC, nodeVec[i]->solnVarData[pos] );
        }
      }
      else
      {
        Report::UserWarning() << "Cannot restore solution variables for node : " << nodeVec[i]->ID << std::endl;
      }

      if ( cnP->type() == _DNODE )
      {
        CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(cnP);
        if (cktNodeDevPtr->stateVarCount())
        {
          const std::vector<int> & gidList2 = cktNodeDevPtr->get_StateVarGIDList();
          pos = 0;
          // Check if the new state variable data is the same size as the old.
          // Do not set old state variable data if they do not match, model has changed.
          if (gidList2.size() == nodeVec[i]->stateVarData.size())
          {
            for( std::vector<int>::const_iterator iterIC = gidList2.begin();
                 iterIC != gidList2.end(); ++iterIC, ++pos )
            {
              analysis_manager.setStateVarData( *iterIC, nodeVec[i]->stateVarData[pos] );
            }
          }
          else
          {
            Report::UserWarning() << "Cannot restore state variables for node : " << nodeVec[i]->ID << std::endl;
          }
        }

        if (cktNodeDevPtr->storeVarCount())
        {
          const std::vector<int> & gidList3 = cktNodeDevPtr->get_StoreVarGIDList();
          pos = 0;
          // Check if the new store variable data is the same size as the old.
          // Do not set old store variable data if they do not match, model has changed.
          if (gidList3.size() == nodeVec[i]->storeVarData.size())
          {
            for( std::vector<int>::const_iterator iterIC = gidList3.begin();
                 iterIC != gidList3.end(); ++iterIC, ++pos )
            {
              analysis_manager.setStoreVarData( *iterIC, nodeVec[i]->storeVarData[pos] );
            }
          }
          else
          {
            Report::UserWarning() << "Cannot restore store variables for node : " << nodeVec[i]->ID << std::endl;
          }
        }

        if( nodeVec[i]->devState != NULL && nodeVec[i]->devState != 0 )
        {
          cktNodeDevPtr->setDevState( *nodeVec[i]->devState );
        }
      }
    }
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Topology::outputNameFile
// Purpose       : This is a kludgy function designed to output all the
//                 solution variable indices and their respective names
//                 to a file.
//
//                 The original point was to create a file to compare
//                 with SPICE, so the names needed to be as similar
//                 as possible to SPICE's naming convention.
//
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/07/01
//-----------------------------------------------------------------------------
bool Topology::outputNameFile(
  Parallel::Machine     comm,
  const std::string &   path,
  bool                  override_output)
{
  loadNodeSymbols();


  if (linearSolverUtility_->namesFileFlag() || override_output)
  {
    Indexor indexor(pdsManager_);

    std::ostringstream oss;
    std::vector<int> id(1);
    for (int i = 0, size = solutionNodeNames_.size(); i < size; ++i)
    {
      if (*solutionNodeNames_[i] != gnd_) {
        id[0] = i;
        std::string s = *solutionNodeNames_[i];
        Util::toLower(s);
        indexor.localToGlobal(Parallel::SOLUTION, id);

        oss << "\t" << id[0] << "\t" << std::setw(12) << s << std::endl;
      }
    }

    Util::Marshal mout;
    mout << oss.str();

    std::vector<std::string> dest;
    Parallel::GatherV(comm, 0, mout.str(), dest);

    if (Parallel::rank(comm) == 0) {
      std::ofstream output_stream(path.c_str(), std::ios_base::out);

      if (output_stream.fail())
      {
        Report::UserWarning() << "Unable to open names file" <<std::endl;
      }
      else
      {
        output_stream << "HEADER" << std::endl;
        for (int p = 0; p < Parallel::size(comm); ++p) {
          Util::Marshal min(dest[p]);

          std::string s;
          min >> s;
          output_stream << s;
        }
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void Topology::loadNodeSymbols() const
{
#if 1
  std::cout  << "In Topology::loadNodeSymbols. "
    << "  nodeSymbols_[Util::SOLUTION_SYMBOL].size = " << nodeSymbols_[Util::SOLUTION_SYMBOL].size() 
    << "  nodeSymbols_[Util::STATE_SYMBOL].size = " << nodeSymbols_[Util::STATE_SYMBOL].size() 
    << "  nodeSymbols_[Util::STORE_SYMBOL].size = " << nodeSymbols_[Util::STORE_SYMBOL].size() 
    << "  nodeSymbols_[Util::EXTERN_SYMBOL].size = " << nodeSymbols_[Util::EXTERN_SYMBOL].size() 
    << "  nodeSymbols_[Util::BRANCH_SYMBOL].size = " << nodeSymbols_[Util::BRANCH_SYMBOL].size() 
    << std::endl;
#endif

  generateOrderedNodeList();

  if (nodeSymbols_[Util::SOLUTION_SYMBOL].empty()) {
    for (CktNodeList::const_iterator it = orderedNodeListPtr_->begin(), end  = orderedNodeListPtr_->end(); it != end; ++it)
    {
      if ((*it)->get_IsOwned())
      {
        (*it)->loadNodeSymbols(const_cast<Topology &>(*this));
      }
    }
    // Resize the solution node names list and add one more entry for the ground node.
    solutionNodeNames_.resize(nodeSymbols_[Util::SOLUTION_SYMBOL].size() + 1, &gnd_);
    for (NodeNameMap::const_iterator it = nodeSymbols_[Util::SOLUTION_SYMBOL].begin(), end = nodeSymbols_[Util::SOLUTION_SYMBOL].end(); it != end; ++it) {
      ThrowRequire((*it).second < solutionNodeNames_.size());

      solutionNodeNames_[(*it).second] = &(*it).first;
    }
  }
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
const std::vector<char> & Topology::getVarTypes() const
{
  generateOrderedNodeList();

  if (variableTypes_.empty())
  {
    for (CktNodeList::iterator it = orderedNodeListPtr_->begin(), 
        end = orderedNodeListPtr_->end(); it != end; ++it )
    {
      if( (*it)->get_IsOwned() && ( (*it)->get_gID() != -1 ) )
      {
          (*it)->varTypeList(variableTypes_); // variableTypes_.push_back( 'V' );
      }
    }
  }

  return variableTypes_;
}

//-----------------------------------------------------------------------------
// Function      : Topology::addResistors
// Purpose       : Adds resistors (between ground and nodes which are connected
//                 to only one device terminal) to a copy of the netlist file.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 1/1/18
//-----------------------------------------------------------------------------
void Topology::writeNetlistAddResistors
(unordered_set<std::string> &        connToOneTermIDs,
 unordered_set<std::string> &        noDCPathIDs)
{
  N_PDS_Comm * commPtr = pdsManager_.getPDSComm();
  int procID = commPtr->procID();

  if (procID == 0 && hangingResistor_.getNetlistCopy())
  {
    // We let 'connected to one terminal' take precedence
    // over 'no DC path.'  (Don't want to label a node as
    // both* being connected to only one terminal and having no
    // DC path to ground since this will cause the addition of two resistors.
    unordered_set<std::string>::iterator it = connToOneTermIDs.begin();
    unordered_set<std::string>::iterator it_end = connToOneTermIDs.end();
    for ( ; it != it_end; ++it )
    {
      unordered_set<std::string>::iterator it_noDC = noDCPathIDs.find( *it );
      if (it_noDC != noDCPathIDs.end())
        noDCPathIDs.erase( it_noDC );
    }

    bool oneTermNotNoDCPath = true;

    //We use this boolean to print a different banner depending upon whether
    //resistors are being added because they are connected to only one device
    //terminal, or if they are being added because they have no DC path to
    //ground.

    //append resistors to nodes with only one terminal connection.

    if (!connToOneTermIDs.empty())
    {
      std::string onetermres(hangingResistor_.getOneTermRes());
      addResistors(connToOneTermIDs,onetermres,oneTermNotNoDCPath);
    }


    //append resistors to nodes with no dc path to ground.
    if (!noDCPathIDs.empty())
    {
      std::string nodcpathres(hangingResistor_.getNoDCPathRes());
      addResistors(noDCPathIDs,nodcpathres,!oneTermNotNoDCPath);

    }

    //if we've requested to produce a netlist file with resistors between
    //dangling nodes and ground, but it turns out that there aren't any
    //dangling nodes, we just need to add a ".END" to the end of the netlist
    //file copy that we produced in the I/O Interface Package (not critical).
    appendEndStatement();
    
  }
} 

//-----------------------------------------------------------------------------
// Function      : Topology::addResistors
// Purpose       : Adds resistors (between ground and nodes which are connected
//                 to only one device terminal) to a copy of the netlist file.
// Special Notes :
// Scope         : private
// Creator       : Keith Santarelli, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/11/07
//-----------------------------------------------------------------------------
void Topology::addResistors
(const unordered_set<std::string> &     inputVec,
 const std::string &                    resValue,
 bool                                   oneTermNotNoDCPath)
{
  std::string netlistCopy(commandLine_.getArgumentValue("netlist"));
  netlistCopy += "_xyce.cir";
  std::ofstream copyFile;
  copyFile.open(netlistCopy.c_str(),std::ios::app);  //put in append mode

  std::string msg("");

  if (DEBUG_IO) {
    if (oneTermNotNoDCPath)
    {
      msg = "Adding resistors of value ";
      //msg += resvalue;
      msg += " between ground and nodes connected to only one device";
      msg += "terminal in file ";
      msg += netlistCopy;
    }
    else
    {
      msg = "Adding resistors of value ";
      //msg += resvalue;
      msg += " between ground and nodes with no DC path to ground in file ";
      msg += netlistCopy;
    }

    Xyce::dout() << msg << std::endl;
  }

  //Some error checking in case we can't open the file.
  if(copyFile.fail())
  {
    if (oneTermNotNoDCPath)
    {
      Report::UserError() 
        << "Error adding resistors between ground and nodes connected to only one device terminal: cannot open file " 
        << netlistCopy;
    }
    else
    {
      Report::UserError() 
        << "Error adding resistors between ground and nodes with no DC path to ground: cannot open file " 
        << netlistCopy;
    }

    return;
  }

  std::string banner("");

  if (oneTermNotNoDCPath)
  {
    banner = "*XYCE-GENERATED OUTPUT:  Adding resistors between ground and ";
    banner += "nodes connected to only 1 device terminal:";
  }
  else
  {
    banner = "*XYCE-GENERATED OUTPUT:  Adding resistors between ground and ";
    banner += "nodes with no DC path to ground:";
  }

  copyFile << std::endl << std::endl << banner << std::endl << std::endl;

  //Now, loop through the ids in inputVec and add the resistors.
  unordered_set<std::string>::const_iterator inputit = inputVec.begin();
  unordered_set<std::string>::const_iterator inputend = inputVec.end();
  int count = 0;

  for ( ; inputit != inputend; ++inputit, ++count )
  {
    std::string resname("R");

    if (oneTermNotNoDCPath)
    {
      resname += "ONETERM";
    }
    else
    {
      resname += "NODCPATH";
    }

    copyFile << resname;
    copyFile << count+1 << " " << *inputit << " 0 " << resValue;
    copyFile << std::endl;
  }
  copyFile.close();
}

//-----------------------------------------------------------------------------
// Function      : Topology::appendEndStatement
// Purpose       : Adds ".END" to a copy of the netlist file.
// Special Notes :
// Scope         : private
// Creator       : Keith Santarelli, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/11/07
//-----------------------------------------------------------------------------
void Topology::appendEndStatement()
{
  std::string netlistCopy(commandLine_.getArgumentValue("netlist"));
  netlistCopy += "_xyce.cir";
  std::ofstream copyFile;
  copyFile.open(netlistCopy.c_str(),std::ios::app);  //put in append mode

  //Some error checking in case we can't open the file.
  if(copyFile.fail())
  {
    Report::UserError() << "Attempt to append .END statement as part of netlist copy procedure:  Cannot open file " << netlistCopy;
  }

  copyFile << std::endl << ".END" << std::endl;
  copyFile.close();
}

//-----------------------------------------------------------------------------
// Function      : Topology::outputTopoWarnings
// Purpose       : Output warnings from node connectivity checks
// Special Notes : Warnings are output based on preprocess settings.
// Scope         : public 
// Creator       : Dave Shirley, Heidi Thornquist, SNL
// Creation Date : 06/21/05
//-----------------------------------------------------------------------------
void Topology::outputTopoWarnings(unordered_set<std::string> & oneTermName,
                                  unordered_set<std::string> & noDCName)
{
  N_PDS_Comm & comm = *(pdsManager_.getPDSComm());
  int procID = comm.procID();

  std::string warn1("connected to only 1 device Terminal");
  std::string warn2("does not have a DC path to ground");

  if (!comm.isSerial())
  {
    generateGlobalNameSet( oneTermName );  
    generateGlobalNameSet( noDCName );
  }

  if (procID == 0 && !oneTermName.empty())
  {
    unordered_set<std::string>::const_iterator it = oneTermName.begin();
    unordered_set<std::string>::const_iterator it_end = oneTermName.end();
    for ( ; it != it_end; ++it )
    {
      std::string msg("Voltage Node (" + *it + ") " + warn1);
      Report::UserWarning0() << msg;
    }
  }

  if (procID == 0 && !noDCName.empty())
  {
    unordered_set<std::string>::const_iterator it = noDCName.begin();
    unordered_set<std::string>::const_iterator it_end = noDCName.end();
    for ( ; it != it_end; ++it )
    {
      std::string msg("Voltage Node (" + *it + ") " + warn2);
      Report::UserWarning0() << msg;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::generateGlobalNameSet
// Purpose       : Collect an unordered set of strings on one processor.
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley, Heidi Thornquist, SNL
// Creation Date : 06/21/05
//-----------------------------------------------------------------------------
void Topology::generateGlobalNameSet( unordered_set<std::string> & nodeName )
{
  N_PDS_Comm & comm = *(pdsManager_.getPDSComm());
  int procCnt = comm.numProc();
  int procID = comm.procID();

  std::vector<int> candidates_all (procCnt,0);
  std::vector<int> candidates_sum (procCnt,0);

  int bs;
  candidates_all[procID] = nodeName.size();
  comm.sumAll(&candidates_all[0], &candidates_sum[0], procCnt);

  unordered_set<std::string>::const_iterator it = nodeName.begin();
  unordered_set<std::string>::const_iterator it_end = nodeName.end();

  if (procID != 0)
  {
    for ( ; it != it_end; ++it )
    {
      bs = (*it).size();
      comm.send (&bs, 1, 0);
      comm.send (&(*it)[0], bs, 0);
    }
  }
  else
  {
    std::string buf;

    for (int i=1 ; i<procCnt ; ++i)
    {
      for (int j=0 ; j<candidates_sum[i] ; ++j)
      {
        comm.recv (&bs, 1, i);
        buf.resize(bs);
        comm.recv (&buf[0], bs, i);
        nodeName.insert( buf );
      }
    }
  }
}

} // namespace Topo
} // namespace Xyce

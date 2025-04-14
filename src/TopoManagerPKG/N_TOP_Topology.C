//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
#include <numeric>

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
  Parallel::Manager &           pds_manager)
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
  std::vector<bool> solnFndGndVec;
  std::vector<NodeID> solnNidVec;
  std::vector<NodeID> idVec;
  NodeID gndNode( "0", _VNODE );

  // Find out if any device nodes have dependent variables that need to be resolved.
  for (CktNodeList::iterator it_cnL = mainGraphPtr_->getBFSNodeList()->begin(),
      it_cnL_end = mainGraphPtr_->getBFSNodeList()->end(); it_cnL != it_cnL_end; ++it_cnL )
  {
    if ( (*it_cnL)->type() == _DNODE )
    {
      CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_cnL);

      // Get dependent solution variables
      cktNodeDevPtr->getDepSolnVars( idVec );

      solnLocVec.push_back( solnCount );
      solnFndGndVec.push_back( false );

      if( !idVec.empty() )
      {
        for (int i=0; i<idVec.size(); ++i)
        {
          unordered_map<std::string, std::string>::iterator alias_it = aliasMap_.find( Xyce::get_node_id( idVec[i] ) );
          if ( idVec[i] == gndNode )
          {
            (solnFndGndVec.back()) = true;
          }
          else if ( alias_it != aliasMap_.end() )
          {
            solnNidVec.push_back( Xyce::NodeID( alias_it->second, _VNODE ) );
            solnCount++;  
          }
          else
          {
            solnNidVec.push_back( idVec[i] );
            solnCount++;
          }
        }
      }
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
    std::vector<int> gndGID( 1, -1 );

    // Get the solution GIDs and check for bad solution nodes.
    std::vector<NodeID> badSolnNodes = directory.getSolnGIDs( solnNidVec, gidVec, procVec );

    // Compute global number of bad solution nodes
    int local_bSN_size = badSolnNodes.size(), global_bSN_size = 0;
    pdsManager_.getPDSComm()->sumAll(&local_bSN_size, &global_bSN_size, 1);

    // Output error messages for bad solution nodes
    checkForBadSolnNodes( badSolnNodes );

    if (global_bSN_size == 0)
    {
      int count = 0;
      for (CktNodeList::iterator it_cnL = mainGraphPtr_->getBFSNodeList()->begin(),
         it_cnL_end = mainGraphPtr_->getBFSNodeList()->end(); it_cnL != it_cnL_end; ++it_cnL )
      {
        if ( (*it_cnL)->type() == _DNODE )
        {
          CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_cnL);

          if (solnLocVec[count] != solnLocVec[count+1])
          {
            
            indexVec.assign( gidVec.begin()+solnLocVec[count],
                             gidVec.begin()+solnLocVec[count+1] );

            if (solnFndGndVec[count] == true)
            {
              // Get dependent solution variables
              cktNodeDevPtr->getDepSolnVars( idVec );

              for (int i=0; i<idVec.size(); ++i)
              {
                if ( idVec[i] == gndNode )
                {
                  indexVec.insert( indexVec.begin()+i, gndGID );
                }
              } 
            }
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
        depSolnGIDMap_[-1] = pdsManager_.getPDSComm()->procID();
      }
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
// Function      : Topology::removeFloatingNodes
// Purpose       : Takes a look at the current devices in the graph to detect
//               : imbalance before completing topology setup
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 7/1/21
//-----------------------------------------------------------------------------
void Topology::removeFloatingNodes(Device::DeviceMgr &device_manager)
{
  if( linearSolverUtility_->floatingnodeFlag() ) 
  {
    int numProcs = pdsManager_.getPDSComm()->numProc();
    int myProc = pdsManager_.getPDSComm()->procID();

    // First identify any devices (on this processor) that are 
    // disconnected from the rest of the circuit.
    std::ostringstream oss;
    std::vector< Xyce::NodeID > floatingDevs = mainGraphPtr_->analyzeDeviceNodeGraph( oss );

    std::vector< int > floatingDevNodes;
    floatingDevNodes.push_back(0); // For the start pointer

    // Identify the voltage nodes attached to these disconnected devices, 
    // which are also floating and might need to be removed.
    std::vector< Xyce::NodeID > floatingNodes;
    for ( std::vector< Xyce::NodeID >::iterator it = floatingDevs.begin(); it != floatingDevs.end(); ++it )
    {
      std::vector< Xyce::NodeID > v_ids;
      mainGraphPtr_->returnAdjIDs(*it, v_ids);
      floatingNodes.insert( floatingNodes.end(), v_ids.begin(), v_ids.end() );
      floatingDevNodes.push_back( floatingDevNodes.back() + v_ids.size() );
    }

    int tNumFloatingDevs = floatingDevs.size();
    int tNumFloatingNodes = floatingNodes.size();

    unordered_map< Xyce::NodeID, int > devicesToSend;

    // In parallel, a floating device/node on one processor might be connected to devices 
    // on other processors.  Communicate the voltage nodes connected to these devices to
    // find out.
    if (!pdsManager_.getPDSComm()->isSerial())
    {
      std::vector<int> tFNodes( numProcs, 0 ), fNodes( numProcs, 0 );

      //Xyce::dout() << "Processor " << myProc << " has " << floatingDevs.size() << " floating devices to consider " << std::endl;
      //Xyce::dout() << "Processor " << myProc << " has " << floatingNodes.size() << " floating nodes to consider " << std::endl;

      fNodes[ myProc ] = floatingNodes.size();
      pdsManager_.getPDSComm()->sumAll( &fNodes[0], &tFNodes[0], numProcs );

      // Count the bytes for packing these strings
      int byteCount = 0;
      for ( std::vector< Xyce::NodeID >::iterator it = floatingNodes.begin(); it != floatingNodes.end(); ++it )
        byteCount += Xyce::packedByteCount( *it );

      // Collect finalFloating nodes and devices for local removal. 
      std::vector< Xyce::NodeID > finalFloatingDevs, finalFloatingNodes;

      for (int proc = 0; proc < numProcs; ++proc)
      {
        // Broadcast the buffer size for this processor.
        int bsize=0;
        if (proc==myProc) { bsize = byteCount; }
        pdsManager_.getPDSComm()->bcast( &bsize, 1, proc );

        // Create buffer.
        int pos = 0;
        char * floatingNodeBuffer = new char[bsize];

        // IMPORTANT DETAIL:  Set the flag to -(numProcs+1) initially. 
        //  - If a processor also has that non-floating node, it will set the flag to (myProc+1).  
        //  - If that node is also a floating node on this processor, set it to -(myProc+1).
        // That way, a processor with a non-floating node will be chosen to send the device to.
        // If multiple processors have the same floating node, they could be consolidated on one
        // processor to reduce the communication.  
        std::vector< int > nodeCnt( tFNodes[proc], -(numProcs+1) ), tmpNodeCnt( tFNodes[proc], -(numProcs+1) );

        if ( proc == myProc )
        {
          // Pack the floating node buffer
          for (std::vector< Xyce::NodeID >::iterator node = floatingNodes.begin();
            node != floatingNodes.end(); node++)
          {
            // Pack the node
            Xyce::pack(*node, floatingNodeBuffer, bsize, pos, pdsManager_.getPDSComm());
          }

          // Broadcast packed buffer.
          pdsManager_.getPDSComm()->bcast( floatingNodeBuffer, bsize, proc );

          // Wait for a response from the other processors on whether they have 
          // any of these floating nodes.
          pdsManager_.getPDSComm()->maxAll( tmpNodeCnt.data(), nodeCnt.data(), tFNodes[proc] );

          std::vector< Xyce::NodeID >::iterator it = floatingDevs.begin();
          for (int i=0; it != floatingDevs.end(); ++it, ++i)
          {
            int sendToProc = -1;
            for ( int j=floatingDevNodes[i]; j<floatingDevNodes[i+1]; ++j )
            {
              // This node has been found on another processor, and is not a floating node.
              if ( nodeCnt[j] > 0 )
                sendToProc = nodeCnt[j]-1;  // Remove the 1-based processor offset
            }
            // If this device is truly floating or is going to be sent to another processor,
            // prepare for its removal.
            if (sendToProc > -1)
            {
              devicesToSend[ *it ] = sendToProc;
            }
            else
            {
              finalFloatingDevs.push_back( *it );
              for ( int j=floatingDevNodes[i]; j<floatingDevNodes[i+1]; ++j )
                finalFloatingNodes.push_back( *(floatingNodes.begin()+j) );
            }
          }

          //Xyce::dout() << "Processor " << myProc << " needs to send " << devicesToSend.size() << " device(s)!" << std::endl;
        }
        else
        {
          // Unpack buffer.
          pdsManager_.getPDSComm()->bcast( floatingNodeBuffer, bsize, proc );

          // Extract floating node pairs.
          std::vector< NodeID > tmpNode( tFNodes[proc] );

          // Unpack each node.
          for (int i=0; i<tFNodes[proc]; ++i)
            Xyce::unpack(tmpNode[i], floatingNodeBuffer, bsize, pos, pdsManager_.getPDSComm());

          // Check if this processor has any of these floating nodes in the graph.
          for (int i=0; i<tFNodes[proc]; ++i)
          {
            // Check if this node is on this processor
            if ( mainGraphPtr_->FindCktNode( tmpNode[i] ) != 0 )
              tmpNodeCnt[i] = myProc+1;
          }

          // Now sum the tmpNodeCnt to see which voltage nodes are on other processors. 
          pdsManager_.getPDSComm()->maxAll( tmpNodeCnt.data(), nodeCnt.data(), tFNodes[proc] );
        }

        // Clean up.
        delete [] floatingNodeBuffer;
      }

      floatingDevs = finalFloatingDevs;
      floatingNodes = finalFloatingNodes;

      int numFloatingDevs = floatingDevs.size();
      int numFloatingNodes = floatingNodes.size();

      pdsManager_.getPDSComm()->sumAll( &numFloatingNodes, &tNumFloatingNodes, 1 );
      pdsManager_.getPDSComm()->sumAll( &numFloatingDevs, &tNumFloatingDevs, 1 );
    }

    // Now local floating devices and nodes can be removed.
    if (tNumFloatingDevs)
      Xyce::dout() << "Device verification found " << tNumFloatingDevs << " disconnected device(s) to remove" << std::endl;

    if (tNumFloatingNodes)
      Xyce::dout() << "Device verification found " << tNumFloatingNodes << " disconnected nodes(s) to remove" << std::endl;

    // First remove floating devices for this processor
    std::vector< CktNode * > removedDevices;
    mainGraphPtr_->removeNodes( floatingDevs, removedDevices );

    // now it's safe to delete the device nodes
    for (std::vector< CktNode * >::iterator it = removedDevices.begin(); it != removedDevices.end(); ++it)
      delete *it;

    // Now remove floating voltage nodes for this processor
    std::vector< CktNode * > removedNodes;
    mainGraphPtr_->removeNodes( floatingNodes, removedNodes );

    // now it's safe to delete the voltage nodes
    for (std::vector< CktNode * >::iterator it = removedNodes.begin(); it != removedNodes.end(); ++it)
      delete *it;
  }
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
void Topology::verifyNodesAndDevices(Device::DeviceMgr & device_manager)
{
  if( linearSolverUtility_->supernodeFlag() )
  {
    unordered_set<NodeID> nodesReplaced, redundantNodes;
    std::vector< std::pair<NodeID, NodeID> > redundantSupernodes;

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
                  std::pair<unordered_set<NodeID>::iterator,bool> ret = nodesReplaced.insert( *nextID );
                  if (!ret.second)
                  {
                    //Xyce::dout() << "Voltage node " << *nextID << " is being replaced twice " << std::endl;
                    redundantSupernodes.push_back( std::make_pair( *nextID, *currentID ) );
                    redundantNodes.insert( *nextID );
                  }
                  else
                  {
                    // format of pair is nodeToReplace, nodeToKeep
                    superNodeList_.push_back( std::make_pair( *nextID, *currentID ) );
                  }
                }
                else
                {
                  std::pair<unordered_set<NodeID>::iterator,bool> ret = nodesReplaced.insert( *currentID );
                  if (!ret.second)
                  {
                    //Xyce::dout() << "Voltage node " << *currentID << " is being replaced again" << std::endl;
                    redundantSupernodes.push_back( std::make_pair( *currentID, *nextID ) );
                    redundantNodes.insert( *currentID );
                  }
                  else
                  {
                    // format of pair is nodeToReplace, nodeToKeep
                    superNodeList_.push_back( std::make_pair( *currentID, *nextID ) );
                  }
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
            << "null device instance block pointer";
        }
      }
      currentCktNodeItr++;
    }

    //Xyce::dout() << "Redundant supernodes: " << redundantSupernodes.size() << std::endl;
    // Find the redundant supernodes on this processor and collapse them all.
    unordered_set<Xyce::NodeID>::iterator currNodeID = redundantNodes.begin();
    unordered_set<Xyce::NodeID>::iterator endNodeID = redundantNodes.end();
    while( currNodeID != endNodeID )
    {
      std::set<Xyce::NodeID> setKeep; // This is an ordered set, so the first is the lexically smallest 
      Xyce::NodeID replaceNode = *currNodeID;

      // Collect all the kept nodes in setKeep and then replace them with the smallest node.
      // Additionally, new supernodes need to be added to collapse the redundant nodes.
      for (int i=0; i < redundantSupernodes.size(); ++i)
      {
        if( replaceNode == redundantSupernodes[i].first )
          setKeep.insert( redundantSupernodes[i].second );
      }
      for (int i=0; i < superNodeList_.size(); ++i)
      {
        if( replaceNode == superNodeList_[i].first )
        {
          setKeep.insert( superNodeList_[i].second );
          if( superNodeList_[i].second != *(setKeep.begin()) )
          {
            superNodeList_.push_back( std::make_pair( superNodeList_[i].second, *(setKeep.begin()) ) );
            //Xyce::dout() << "Adding supernode: " << superNodeList_[i].second << ", " << *(setKeep.begin()) << std::endl;
            superNodeList_[i].second = *(setKeep.begin());
          }
          break;
        }
      }
      // Collapse redundant nodes
      for (int i=0; i < redundantSupernodes.size(); ++i)
      {
        if( replaceNode == redundantSupernodes[i].first )
        {
          if( redundantSupernodes[i].second != *(setKeep.begin()) )
          {
            superNodeList_.push_back( std::make_pair( redundantSupernodes[i].second, *(setKeep.begin()) ) );
            //Xyce::dout() << "Adding supernode: " << redundantSupernodes[i].second << ", " << *(setKeep.begin()) << std::endl;
          }
        }
      }

      currNodeID++;
    }
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
unordered_map<std::string, std::string> 
Topology::removeTaggedNodesAndDevices(Device::DeviceMgr &device_manager)
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

    std::vector< std::pair<NodeID, NodeID> >::iterator currentNodePair = superNodeList_.begin();
    std::vector< std::pair<NodeID, NodeID> >::iterator endNodePair = superNodeList_.end();
    unordered_set<NodeID> nodesReplaced;
    while ( currentNodePair != endNodePair )
    {
      NodeID nodeToBeReplaced( currentNodePair->first );
      NodeID replacementNode( currentNodePair->second );
      //Xyce::dout() << "Replacing node " << nodeToBeReplaced << " with " << replacementNode << std::endl;
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
                  *nextNodePair = std::make_pair( nextNodePair->second, replacementNode );
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
                  *nextNodePair = std::make_pair( replacementNode, nextNodePair->first );
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
      mainGraphPtr_->removeRedundantDevices( badDeviceList_, removedDevices );
    }

    int tNumBadDevices = 0;
    int numBadDevices = badDeviceList_.size();
    pdsManager_.getPDSComm()->sumAll( &numBadDevices, &tNumBadDevices, 1 );
    if (tNumBadDevices)
    {
      Xyce::dout() << "Device verification found " << tNumBadDevices << " device(s) to remove" << std::endl;
    }

    // Communicate the node aliases for the removed nodes, which is held by the superNodeList_
    std::vector< std::pair<NodeID, NodeID> >::iterator it = superNodeList_.begin();
    std::vector< std::pair<NodeID, NodeID> >::iterator it_end = superNodeList_.end();
    for ( ; it != it_end; ++it )
      aliasMap_[ Xyce::get_node_id( (*it).first ) ] = Xyce::get_node_id( (*it).second );

    // now it's safe to delete the old nodes and device nodes
    for (CktNodeList::iterator it = oldNodeList.begin(), end = oldNodeList.end(); it != end; ++it)
    {
      delete *it;
    }

    if (!removedDevices.empty())
    {
      std::vector< std::string > badDeviceNames;

      // delete old devices that were removed
      for (std::vector<CktNode *>::iterator it = removedDevices.begin(), end = removedDevices.end(); it != end; ++it)
      {
        CktNode_Dev *cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it);

        if (cktNodeDevPtr && cktNodeDevPtr->deviceInstanceBlock())
        {
          // Collect information from the device instance block
          badDeviceNames.push_back( cktNodeDevPtr->deviceInstanceBlock()->getInstanceName().getEncodedName() );  
          delete cktNodeDevPtr;
        }
      }
     
      // Communicate removed devices to the device manager
      device_manager.registerRemovedDevices( badDeviceNames );
    }
  }

  return aliasMap_;
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
    Parallel::Communicator * commPtr = pdsManager_.getPDSComm();
    int numProcs = commPtr->numProc();
    int thisProc = commPtr->procID();

    int localNodes = superNodeList_.size();
    std::vector< std::pair<NodeID, NodeID> > externalSuperNodeList;

    // Count the bytes for packing these strings
    int byteCount = 0;

    // First count the size of the local superNodeList_
    byteCount += sizeof(int);

    // Now count all the NodeIDs in the list
    unordered_set<NodeID> nodesReplaced, replacementNodes;
    for (std::vector< std::pair<NodeID, NodeID> >::iterator nodePair = superNodeList_.begin();
        nodePair != superNodeList_.end(); nodePair++)
    { 
      byteCount += Xyce::packedByteCount(nodePair->first);
      byteCount += Xyce::packedByteCount(nodePair->second);
      nodesReplaced.insert( nodePair->first );
      replacementNodes.insert( nodePair->second );
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

    std::set< NodeID > checkReplacedNodes;
    unordered_map< NodeID, std::set<NodeID> > redundantSN;
    std::vector< std::pair<NodeID, NodeID> >::iterator currentGlobalSN = externalSuperNodeList.begin();
    std::vector< std::pair<NodeID, NodeID> >::iterator endGlobalSN = externalSuperNodeList.end();
    while ( currentGlobalSN != endGlobalSN )
    {
      //Xyce::dout() << "Analyzing global supernode " << currentGlobalSN->first << " with " << currentGlobalSN->second << std::endl; 
      // Add pair if replacement node is node to be removed by another processor.
      if (replacementNodes.find( currentGlobalSN->first ) != replacementNodes.end())
      {
        //Xyce::dout() << "Found matching replacement node and external supernode: " << currentGlobalSN->first << std::endl;
        superNodeList_.push_back( *currentGlobalSN );
      }
      // Add pair if node to be removed is a replacement node on another processor.
      else if (nodesReplaced.find( currentGlobalSN->second ) != nodesReplaced.end())
      {
        //Xyce::dout() << "Found matching node replaced and external supernode: " << currentGlobalSN->second << std::endl;
        superNodeList_.push_back( *currentGlobalSN );
      }
      else if (nodesReplaced.find( currentGlobalSN->first ) != nodesReplaced.end())
      {
        //Xyce::dout() << "Another processor is replacing node " << currentGlobalSN->first << " with " << currentGlobalSN->second << std::endl;
        redundantSN[ currentGlobalSN->first ].insert( currentGlobalSN->second );
      }
      // There might be a chance that the node to be removed on another processor
      // is a node on this processor, but not a supernode.
      else if ( (nodesReplaced.find( currentGlobalSN->first ) == nodesReplaced.end())
               && mainGraphPtr_->FindCktNode( currentGlobalSN->first ) )
      {
        //Xyce::dout() << "A non-supernode is being replaced " << currentGlobalSN->first << " with " << currentGlobalSN->second << std::endl;
        superNodeList_.push_back( *currentGlobalSN );
      }
      else if ( (replacementNodes.find( currentGlobalSN->second ) == replacementNodes.end())
                && mainGraphPtr_->FindCktNode( currentGlobalSN->second ) )
      {
        checkReplacedNodes.insert( currentGlobalSN->second );
        //Xyce::dout() << "A local node is a replacement node, might want to check if the node being replaced is collapsed " << currentGlobalSN->first << " with " << currentGlobalSN->second << std::endl;
      }
      currentGlobalSN++;
    }

    // Check if there are redundant supernodes from different processors that changes the replacement node.
    std::vector< std::pair< NodeID, NodeID > > redSuperNodes;
    unordered_map< NodeID, std::set<NodeID> >::iterator redSN_iter = redundantSN.begin();
    unordered_map< NodeID, std::set<NodeID> >::iterator redSN_iter_end = redundantSN.end();
    while ( redSN_iter != redSN_iter_end )
    { 
      currentGlobalSN = superNodeList_.begin();
      endGlobalSN = superNodeList_.end();
      while ( currentGlobalSN != endGlobalSN )
      {
        if ( currentGlobalSN->first == redSN_iter->first )
        {
          std::set<NodeID>::iterator redSN_set_iter = (redSN_iter->second).begin();
          std::set<NodeID>::iterator redSN_set_iter_end = (redSN_iter->second).end();
          NodeID replacementNode = *redSN_set_iter;
          if ( replacementNode < currentGlobalSN->second )
          {
            //Xyce::dout() << "Creating new supernode " << currentGlobalSN->second << ", " << *((redSN_iter->second).begin()) << std::endl;
            redSuperNodes.push_back( std::make_pair( currentGlobalSN->second, *((redSN_iter->second).begin()) ) ); 
          }
          redSN_set_iter++;
          while ( redSN_set_iter != redSN_set_iter_end )
          {
            if ( (nodesReplaced.find( *redSN_set_iter ) == nodesReplaced.end())
                 && mainGraphPtr_->FindCktNode( currentGlobalSN->first )
                 && (replacementNode < *redSN_set_iter) )
            {
              //Xyce::dout() << "Creating new supernode " << *redSN_set_iter << ", " << replacementNode << std::endl;
              redSuperNodes.push_back( std::make_pair( *redSN_set_iter, replacementNode ) );
            }
            redSN_set_iter++;
          } 
        }
        currentGlobalSN++;
      }

      redSN_iter++;
    }

    if (redSuperNodes.size())
    {
      currentGlobalSN = superNodeList_.begin();
      endGlobalSN = superNodeList_.end();
      while ( currentGlobalSN != endGlobalSN )
      {
        std::vector< std::pair<NodeID, NodeID> >::iterator iterSN = redSuperNodes.begin();
        std::vector< std::pair<NodeID, NodeID> >::iterator iterSN_end = redSuperNodes.end();
        while ( iterSN != iterSN_end )
        {
          if ( currentGlobalSN->second == iterSN->first )
          {
            currentGlobalSN->second = iterSN->second;
            break;
          }
          iterSN++;
        }
        //Xyce::dout() << "Supernode " << currentGlobalSN->first << ", " << currentGlobalSN->second << std::endl;
        currentGlobalSN++;
      }

      // Add new supernodes to old supernode list.
      superNodeList_.insert( superNodeList_.end(), redSuperNodes.begin(), redSuperNodes.end() );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::finalOutput
// Purpose       : Print output topology details by request
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 4/16/22
//-----------------------------------------------------------------------------
void Topology::finalOutput()
{
  if (linearSolverUtility_->outputGraphFlag() == 1)
  {
    // Stream the circuit graph 
    std::ostringstream oss;
    mainGraphPtr_->streamCircuitGraph( oss );

    // Write out device graph
    std::string filename(commandLine_.getArgumentValue("netlist"));
    filename += "_circuitgraph";
    std::ofstream output_stream(filename.c_str(), std::ios_base::out);

    if (output_stream.fail())
      Report::UserWarning() << "Unable to open circuit graph file" <<std::endl;

    output_stream << oss.str() << std::endl;
  }
  if (linearSolverUtility_->outputGraphFlag() == 2)
  {
    // Compute the device graph (distance-2 graph of circuit graph)
    std::ostringstream oss;
    mainGraphPtr_->analyzeDeviceNodeGraph( oss );

    // Write out device graph
    std::string filename(commandLine_.getArgumentValue("netlist"));
    filename += "_devicegraph";
    std::ofstream output_stream(filename.c_str(), std::ios_base::out);

    if (output_stream.fail())
      Report::UserWarning() << "Unable to open device graph file" <<std::endl;

    output_stream << oss.str() << std::endl;
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
// Purpose       : This is a function designed to output all the
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
  Parallel::Communicator * commPtr = pdsManager_.getPDSComm();
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
  Parallel::Communicator & comm = *(pdsManager_.getPDSComm());
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
// Function      : Topology::checkForBadSolnNodes
// Purpose       : Check if any solution nodes were not found in resolveDependentVars 
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist, SNL
// Creation Date : 08/16/05
//-----------------------------------------------------------------------------
void Topology::checkForBadSolnNodes( const std::vector<NodeID>& badSolnNodes )
{
  if ( badSolnNodes.size() )
  {
    CktNodeList::iterator it_cnL = mainGraphPtr_->getBFSNodeList()->begin();
    CktNodeList::iterator it_cnL_end = mainGraphPtr_->getBFSNodeList()->end(); 
    std::vector<NodeID>::const_iterator it_bSN = badSolnNodes.begin();

    std::vector<NodeID> idVec;
    NodeID gndNode( "0", _VNODE );

    for ( ; it_bSN != badSolnNodes.end(); ++it_bSN )
    {
      bool foundNode = false;

      // Find out if any device nodes have dependent variables that need to be resolved.
      for ( ; it_cnL != it_cnL_end && !foundNode; ++it_cnL )
      { 
        if ( (*it_cnL)->type() == _DNODE )
        { 
          CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_cnL);
      
          // Get dependent solution variables
          cktNodeDevPtr->getDepSolnVars( idVec );
      
          if( !idVec.empty() )
          { 
            for (int i=0; i<idVec.size() && !foundNode; ++i)
            { 
              if ( idVec[i] != gndNode )
              { 
                unordered_map<std::string, std::string>::iterator alias_it = aliasMap_.find( Xyce::get_node_id( idVec[i] ) );
                if ( alias_it != aliasMap_.end() )
                { 
                  if ( Xyce::get_node_id(*it_bSN) == alias_it->second )
                  {
                    Report::UserError() << "Device " << (*it_cnL)->get_id() << " refers to unknown solution node " << Xyce::get_node_id(*it_bSN) << std::endl;
                    foundNode = true;
                  }
                }
                else if ( *it_bSN == idVec[i] )
                { 
                  Report::UserError() << "Device " << (*it_cnL)->get_id() << " refers to unknown solution node " << Xyce::get_node_id(*it_bSN) << std::endl;
                  foundNode = true;
                } 
              }
            }
          }
        }
      }
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
  Parallel::Communicator & comm = *(pdsManager_.getPDSComm());
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

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
// Creation Date  : 07/02/01
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <algorithm>

#include <N_TOP_Directory.h>

#include <N_TOP_Topology.h>
#include <N_TOP_CktGraph.h>
#include <N_TOP_CktNode.h>
#include <N_TOP_CktNode_Dev.h>

#include <N_UTL_Functors.h>

#include <N_PDS_Migrate.h>
#include <N_PDS_Comm.h>

#include <N_PDS_Directory.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : DirectoryData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/26/04
//-----------------------------------------------------------------------------
struct DirectoryData
{
  typedef RCP<ParNode> NodePtr;
  typedef std::map<NodeID,NodePtr> NodeContainer;
                                                                                           
  typedef Xyce::Parallel::Hash<NodeID> NodeIDHash;
  typedef Xyce::Parallel::Migrate<NodeID,ParNode> NodeMigrate;

  typedef Xyce::Parallel::Directory< NodeID,
                                     ParNode,
                                     NodeIDHash,
                                     NodeContainer,
                                     NodeMigrate >  NodeDir;

  DirectoryData( Parallel::Communicator & comm )
  : hash(comm.numProc()),
    migrate(comm),
    directory(migrate,hash)
  {}

  ~DirectoryData() {}

  NodeIDHash hash;
  NodeMigrate migrate;

  NodeDir directory;
};
  

//-----------------------------------------------------------------------------
// Function      : Directory::~Directory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/02/01
//-----------------------------------------------------------------------------
Directory::~Directory()
{
  delete data_;
}

//-----------------------------------------------------------------------------
// Function      : Directory::generateDirectory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/02/01
//-----------------------------------------------------------------------------
bool Directory::generateDirectory()
{
  if (!pdsComm_.isSerial())
  {
    int procID = pdsComm_.procID();

    delete data_;
    data_ = new DirectoryData(pdsComm_);

    DirectoryData::NodeContainer nodes;
    for (CktNodeList::const_iterator iterCN = topology_.getOrderedNodeList().begin(), endCN = topology_.getOrderedNodeList().end(); iterCN != endCN; ++iterCN )
      if( (*iterCN)->get_IsOwned() && (*iterCN)->get_gID() != -1 )
      {
        DirectoryData::NodePtr new_node(
                new ParNode( NodeID( (*iterCN)->get_id(), (*iterCN)->get_gID() ),
                                     true, procID ) );
        nodes[(*iterCN)->get_nodeID()] = new_node;
      }
    
    data_->directory.addEntries( nodes );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Directory::getProcs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/19/01
//-----------------------------------------------------------------------------
std::vector<NodeID> Directory::getProcs( const std::vector<NodeID> & idVec,
                                         std::vector<int> & procVec )
{
  int size = idVec.size();
  procVec.resize( size, -1 );

  std::vector<NodeID> badIDs;

  if (!pdsComm_.isSerial())
  {
    DirectoryData::NodeDir::DataMap nodes;

    std::vector<NodeID> ids( idVec );
    badIDs = data_->directory.getEntries( ids, nodes );

    for( int i = 0; i < size; ++i )
    {
      if ( std::find( badIDs.begin(), badIDs.end(), idVec[i] ) == badIDs.end() )
        procVec[i] = nodes[idVec[i]]->proc();
    }
  }
  else
  {
    std::fill( procVec.begin(), procVec.end(), pdsComm_.procID() );
  }

  return badIDs;
}

//-----------------------------------------------------------------------------
// Function      : Directory::getSolnGIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/19/01
//-----------------------------------------------------------------------------
std::vector<NodeID> Directory::getSolnGIDs( const std::vector<NodeID> & idVec,
                                            std::vector< std::vector<int> > & gidVec,
                                            std::vector<int> & procVec )
{
  gidVec.resize( idVec.size() );

  // This method will find bad NodeIDs that weren't caught before in parallel execution
  // NOTE:  Serial execution will not find bad NodeIDs, they are caught below.
  std::vector<NodeID> badSolnNodes = getProcs( idVec, procVec );

  int local_bSN_size = badSolnNodes.size(), global_bSN_size = 0;
  pdsComm_.sumAll(&local_bSN_size, &global_bSN_size, 1); 

  if (global_bSN_size==0)
  {
    if (!pdsComm_.isSerial())
    {
      typedef Xyce::Parallel::Migrate< NodeID, std::vector<int> > SGMigrate;
      SGMigrate migrator(pdsComm_);

      std::vector<int> sortedProcVec( procVec );
      std::vector<NodeID> ids( idVec );
      SortContainer2( sortedProcVec, ids );

      std::vector<NodeID> inIDs;
      migrator( sortedProcVec, ids, inIDs );

      SGMigrate::DataMap outGIDs;
      for( unsigned int i = 0; i < inIDs.size(); ++i )
      {
        RCP< std::vector<int> > gids( rcp( new std::vector<int>() ) );
        const CktNode * cnp = topology_.findCktNode( inIDs[i] );
        gids->assign( cnp->get_SolnVarGIDList().begin(), cnp->get_SolnVarGIDList().end() );
        outGIDs[inIDs[i]] = gids;
      }

      SGMigrate::DataMap inGIDs;
      migrator.rvs( sortedProcVec, inIDs, outGIDs, inGIDs );

      for( unsigned int i = 0; i < idVec.size(); ++i )
        gidVec[i] = *(inGIDs[idVec[i]]);
    }
    else
    {
      for( unsigned int i = 0; i < idVec.size(); ++i )
      {
        const CktNode * cnp = topology_.findCktNode( idVec[i] );
        if( cnp )
        {
          gidVec[i].assign( cnp->get_SolnVarGIDList().begin(),
                          cnp->get_SolnVarGIDList().end() );
        }
        else
        {
          badSolnNodes.push_back( idVec[i] );
          gidVec[i].push_back( -1 );  // Put a dummy number as the GID value, error will be thrown.
        }
      }
    }
  }

  else
  {
    // In parallel, the bad solution node may be reported on a different processor by the parallel
    // directory. So, it is necessary to communicate all the bad solution nodes to create the correct 
    // error message.  First set the badSolnNodes to the local badSolnNodes and add off-processor
    // contributions.
    int numProc = pdsComm_.numProc();
    int procID = pdsComm_.procID();
    std::vector<int> local_numBadNodes( numProc, 0 ), numBadNodes( numProc, 0 );
    local_numBadNodes[ procID ] = local_bSN_size;   
    pdsComm_.sumAll(&local_numBadNodes[0], &numBadNodes[0], numProc );

    // Count the bytes for packing the strings for the bad nodes.
    int byteCount = 0;
    for ( std::vector< Xyce::NodeID >::iterator it = badSolnNodes.begin(); it != badSolnNodes.end(); ++it )
      byteCount += Xyce::packedByteCount( *it );

    // Loop over the processors and broadcast any bad solution nodes.
    for (int proc = 0; proc < numProc; ++proc)
    {
      // Broadcast the buffer size for this processor.
      int bsize=0;
      if (proc == procID) { bsize = byteCount; }
      pdsComm_.bcast( &bsize, 1, proc );
        
      if (bsize)
      {
        // Create buffer.
        int pos = 0;
        char * badNodeBuffer = new char[bsize];
        
        if ( proc == procID )
        {
          // Pack the bad node buffer
          for (int i=0; i<local_bSN_size; ++i)
            Xyce::pack(badSolnNodes[i], badNodeBuffer, bsize, pos, &pdsComm_);
        }

        // Broadcast packed buffer.
        pdsComm_.bcast( badNodeBuffer, bsize, proc );

        if ( proc != procID )
        {
          // Extract bad nodes and append them to the badSolnNodes list
          NodeID tmpNode;

          // Unpack each node and append to the badSolnNodes vector
          for (int i=0; i<numBadNodes[proc]; ++i)
          {
            Xyce::unpack(tmpNode, badNodeBuffer, bsize, pos, &pdsComm_);
            badSolnNodes.push_back( tmpNode );
          }
        }          
  
        // Clean up.
        delete [] badNodeBuffer;
      }
    } 
  }

  return badSolnNodes;
}

} // namespace Topo
} // namespace Xyce

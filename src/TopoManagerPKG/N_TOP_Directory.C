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

#include <N_PDS_Manager.h>
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

  DirectoryData( N_PDS_Comm & comm )
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
bool Directory::getProcs( const std::vector<NodeID> & idVec,
                                std::vector<int> & procVec )
{
  int size = idVec.size();
  procVec.resize( size );

  if (!pdsComm_.isSerial())
  {
    DirectoryData::NodeDir::DataMap nodes;

    std::vector<NodeID> ids( idVec );
    data_->directory.getEntries( ids, nodes );

    for( int i = 0; i < size; ++i )
      procVec[i] = nodes[idVec[i]]->proc();
  }
  else
  {
    std::fill( procVec.begin(), procVec.end(), pdsComm_.procID() );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Directory::getSolnGIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/19/01
//-----------------------------------------------------------------------------
bool Directory::getSolnGIDs( const std::vector<NodeID> & idVec,
                                   std::vector< std::vector<int> > & gidVec,
                                   std::vector<int> & procVec )
{
  gidVec.resize( idVec.size() );
  getProcs( idVec, procVec );

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
      if( !cnp )
      {
        Report::DevelFatal() << "Directory node not found: " << idVec[i].first + "\n";
      }
      gidVec[i].assign( cnp->get_SolnVarGIDList().begin(),
                        cnp->get_SolnVarGIDList().end() );
    }
  }

  return true;
}

} // namespace Topo
} // namespace Xyce

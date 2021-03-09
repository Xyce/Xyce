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

//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <numeric>
#include <iterator>

#include <N_ERH_ErrorMgr.h>
#include <N_PDS_Comm.h>
#include <N_PDS_Directory.h>
#include <N_PDS_GlobalAccessor.h>
#include <N_PDS_Manager.h>
#include <N_PDS_Migrate.h>
#include <N_PDS_Node.h>
#include <N_PDS_ParMap.h>
#include <N_TOP_CktGraph.h>
#include <N_TOP_CktNode.h>
#include <N_TOP_CktNode_Dev.h>
#include <N_TOP_ParLSUtil.h>
#include <N_TOP_Topology.h>
#include <N_IO_HangingResistor.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Misc.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Stats.h>
#include <N_UTL_fwd.h>

#include <Epetra_Util.h>
#include <Teuchos_Utils.hpp>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : ParLSUtil::ParLSUtil
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter  SNL, Parallel Computational Sciences
// Creation Date : 5/25/05
//-----------------------------------------------------------------------------
ParLSUtil::ParLSUtil(
  Topology &            topology,
  Parallel::Manager &   pds_manager)
  : Linear::QueryUtil(),
    topology_(topology),
    pdsManager_(pds_manager),
    numGlobalNodes_(0),
    numLocalNodes_(0),
    baseNodeGID_(0),
    numGlobalRows_(0),
    numLocalRows_(0),
    numExternRows_(0),
    numGlobalExternRows_(0),
    baseRowGID_(0),
    numGlobalStateVars_(0),
    numLocalStateVars_(0),
    baseStateVarGID_(0),
    numGlobalStoreVars_(0),
    numLocalStoreVars_(0),
    baseStoreVarGID_(0),
    numGlobalLeadCurrentVars_(0),
    numLocalLeadCurrentVars_(0),
    baseLeadCurrentVarGID_(0),
    numGlobalNZs_(0),
    numLocalNZs_(0)
{
}

//-----------------------------------------------------------------------------
// Function      : ParLSUtil::cleanRowLists
// Purpose       : reduce memory consumption of query utility
// Special Notes : the builder should call this after all maps and graphs
//                 have been made
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 1/25/19
//-----------------------------------------------------------------------------
void ParLSUtil::cleanRowLists()
{
  nodeList_GID_.clear();
  nodeList_ExternGID_.clear();

  // rowList_GID_.clear(); keep this for separated object analysis
  rowList_ExternGID_.clear(); 

  rowList_StateGID_.clear();
  rowList_StoreGID_.clear();
  rowList_LeadCurrentGID_.clear();

  rowList_NumNZs_.clear();
  for( int i = 0; i < rowList_ColList_.size(); ++i )
    rowList_ColList_[i].clear();
  rowList_ColList_.clear();

  isClean_ = true; // don't print out arrays, they have been cleared.
}

//-----------------------------------------------------------------------------
// Function      : ParLSUtil::operator<<
// Purpose       : generate utility with reference to topology
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/15/00
//-----------------------------------------------------------------------------
std::ostream & operator<< (std::ostream & os, const ParLSUtil & tlsu )
{
  os << "Parallel Topology LS Utility" << std::endl;
  os << "-------------------" << std::endl;
  os << "ProcID: " << tlsu.pdsManager_.getPDSComm()->procID() << std::endl;
  os << "Num Global Nodes: " << tlsu.numGlobalNodes_ << std::endl;
  os << "Num Local Nodes: " << tlsu.numLocalNodes_ << std::endl;
  os << "Num Global Rows: " << tlsu.numGlobalRows_ << std::endl;
  os << "Num Local Rows: " << tlsu.numLocalRows_ << std::endl;
  os << "Num Extern Rows: " << tlsu.numExternRows_ << std::endl;
  os << "Num Global Extern Rows: " << tlsu.numGlobalExternRows_ << std::endl;
  os << "Row GID Base: " << tlsu.baseRowGID_ << std::endl;

  os << "Num Global State Vars: " << tlsu.numGlobalStateVars_ << std::endl;
  os << "Num Local State Vars: " << tlsu.numLocalStateVars_ << std::endl;
  os << "State Var GID Base: " << tlsu.baseStateVarGID_ << std::endl;

  os << "Num Global Store Vars: " << tlsu.numGlobalStoreVars_ << std::endl;
  os << "Num Local Store Vars: " << tlsu.numLocalStoreVars_ << std::endl;
  os << "Store Var GID Base: " << tlsu.baseStoreVarGID_ << std::endl;

  os << "Num Global LeadCurrent Vars: " << tlsu.numGlobalLeadCurrentVars_ << std::endl;
  os << "Num Local LeadCurrent Vars: " << tlsu.numLocalLeadCurrentVars_ << std::endl;
  os << "LeadCurrent Var GID Base: " << tlsu.baseLeadCurrentVarGID_ << std::endl;

  os << "Num Global NZs: " << tlsu.numGlobalNZs_ << std::endl;
  os << "Num Local NZs: " << tlsu.numLocalNZs_ << std::endl;

  os << "GID Array: ";
  for( int i = 0; i < tlsu.numLocalRows_; ++i )
    os << tlsu.rowList_GID_[i] << " ";
  os << std::endl;

  if (!tlsu.isClean_)
  {
    os << "Extern GID Array: ";
    for( unsigned int i = 0; i < tlsu.rowList_ExternGID_.size(); ++i )
      os << tlsu.rowList_ExternGID_[i].first << " " <<
  	  tlsu.rowList_ExternGID_[i].second << "   ";
    os << std::endl;

    os << "Node Array: ";
    for( int i = 0; i < tlsu.numLocalNodes_; ++i )
      os << tlsu.nodeList_GID_[i] << " ";
    os << std::endl;

    os << "Extern Node Array: ";
    for( unsigned int i = 0; i < tlsu.nodeList_ExternGID_.size(); ++i )
      os << tlsu.nodeList_ExternGID_[i].first << " " <<
  	  tlsu.nodeList_ExternGID_[i].second << "   ";
    os << std::endl;

    os << "State GID Array: ";
    for( int i = 0; i < tlsu.numLocalStateVars_; ++i )
      os << tlsu.rowList_StateGID_[i] << " ";
    os << std::endl;

    os << "Store GID Array: ";
    for( int i = 0; i < tlsu.numLocalStoreVars_; ++i )
      os << tlsu.rowList_StoreGID_[i] << " ";
    os << std::endl;

    os << "LeadCurrent GID Array: ";
    for( int i = 0; i < tlsu.numLocalLeadCurrentVars_; ++i )
      os << tlsu.rowList_LeadCurrentGID_[i] << " ";
    os << std::endl;

    os << "NZ Array: ";
    for( int i = 0; i < tlsu.numLocalRows_; ++i )
      os << tlsu.rowList_NumNZs_[i] << " ";
    os << std::endl;

    os << "Col Index Array: " << std::endl;
    for( int i = 0; i < tlsu.numLocalRows_; ++i )
    {
      os << tlsu.rowList_GID_[i] << ": ";
      for( std::vector<int>::const_iterator it_iL = tlsu.rowList_ColList_[i].begin();
	  it_iL != tlsu.rowList_ColList_[i].end(); ++it_iL )
        os << (*it_iL) << " ";
      os << std::endl;
    }
  }

  return os;
}

//-----------------------------------------------------------------------------
// Function      : ParLSUtil::setupRowCol
// Purpose       : Setup row/col data for linear solver including reorder
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/1/01
//-----------------------------------------------------------------------------
bool ParLSUtil::setupRowCol()
{
  topology_.generateOrderedNodeList();

  Parallel::Communicator & comm = *pdsManager_.getPDSComm();
  int procCnt = comm.numProc();
  int procID = comm.procID();

  //construct v-node and d-node directories
  typedef RCP< Xyce::Parallel::IndexNode > IndexNodePtr;

  typedef std::multimap< std::string, IndexNodePtr > VNodeContainer;
  typedef std::map< std::string, IndexNodePtr > INodeContainer;

  typedef Xyce::Parallel::Hash<std::string> StringHash;
  typedef Xyce::Parallel::Migrate<std::string,Xyce::Parallel::IndexNode> INMigrate;

  typedef Xyce::Parallel::Directory< std::string,
                                     Xyce::Parallel::IndexNode,
                                     StringHash,
                                     VNodeContainer,
                                     INMigrate >
    VNodeDir;

  StringHash SHobj( procCnt );
  INMigrate Mobj( comm );

  //data for directory registration
  INodeContainer VData;

  //vectors for retrieval
  std::vector<std::string> VNames;
  VNames.reserve( int(topology_.getOrderedNodeList().size()/2) );

  //loop over node list and setup IndexNodes for Vs, register with VDir 
  VNodeDir *VDir;

  VDir = new VNodeDir( Mobj, SHobj );

  //find all voltage nodes connected to a voltage source in case we
  // need to force these nodes to be on the same processor
  int dSize=0;
  unordered_set<std::string> Vsrc_Connected_Nodes;
  for (CktNodeList::const_iterator it_cnL = topology_.getOrderedNodeList().begin(), end_cnL = topology_.getOrderedNodeList().end(); it_cnL != end_cnL; ++it_cnL )
  {
    const NodeID& nodeID = (*it_cnL)->get_nodeID();
    const std::string & id = (*it_cnL)->get_id();
    int type = (*it_cnL)->type();

    if (DEBUG_TOPOLOGY)
      Xyce::dout() << "Considering nodes for: " << id << std::endl;

    if  (type == _DNODE)
    {
      std::string::size_type col = id.find_first_of(':');
      if ( id[col+1] == 'V' || id[col+1] == 'v'
           || id.substr(col+1,col+6) == "yiso2"
           || id.substr(col+1,col+6) == "YISO2"
           || id.substr(col+1,col+5) == "yext"
           || id.substr(col+1,col+5) == "YEXT" )
      {
        if (DEBUG_TOPOLOGY)
          Xyce::dout() << "Getting adjacent nodes for: " << id << std::endl;

        std::vector<NodeID> adj_ids;
        topology_.returnAdjIDs( nodeID, adj_ids );
        int adjSize = adj_ids.size();
        for( int i = 0; i < adjSize; ++i )
        {
          if (DEBUG_TOPOLOGY)
            Xyce::dout() << "adj_ids["<<i<<"] = " << adj_ids[i] << std::endl;

          if( Xyce::get_node_id(adj_ids[i]) != "0" )
          {
            Vsrc_Connected_Nodes.insert(adj_ids[i].first);
          }
        }
      }
    
      // Count the number of device nodes. 
      dSize++;
    }
    else if (type == _VNODE && id != "0")
    {
      IndexNodePtr inode( new Xyce::Parallel::IndexNode( -99, procID ) );
      VData[id] = inode;
      VNames.push_back( id );
    }
  }

  if (DEBUG_TOPOLOGY)
    Xyce::dout() << "Vsrc_Connected_Nodes:"<<std::endl;

  unordered_set<std::string>::iterator its = Vsrc_Connected_Nodes.begin();
  unordered_set<std::string>::iterator fts = Vsrc_Connected_Nodes.end();
  for( ; its != fts; ++its )
  {
    if (DEBUG_TOPOLOGY)
      Xyce::dout() << *its << std::endl;

    VData[ *its ]->gid = -98;
  }

  VDir->addEntries( VData );
 
  //v-node use multimap, pick owner
  //locally loop over container in VDir, pick an owner
  //random choice of owner

  int currGID, baseDNodeGID;
  VNodeContainer & VNodes = VDir->container();
  VNodeContainer::iterator iterVN = VNodes.begin();
  VNodeContainer::iterator endVN = VNodes.end();
  VNodeContainer OwnedVNodes;

  std::string id("");
  std::vector<int> procBalance( procCnt, 0 );

  while( iterVN != endVN )
  {
    int proc = -1;

    id = iterVN->first;
    if( iterVN->second->gid == -98 ) proc = iterVN->second->pid;

    IndexNodePtr inode( new Xyce::Parallel::IndexNode( -99, iterVN->second->pid ) );

    ++iterVN;
    if( (iterVN != endVN) && (id == iterVN->first) )
    {
      std::vector<int> intVec;
      intVec.push_back(inode->pid);
      while( (iterVN != endVN) && (iterVN->first == id) )
      {
        if( iterVN->second->gid == -98 ) proc = iterVN->second->pid;

        intVec.push_back( iterVN->second->pid );
        ++iterVN;
      }
      if( proc != -1 )
      {
        inode->pid = proc;
        if (DEBUG_TOPOLOGY)
          Xyce::dout() << "pNode: " << id << " " << proc << std::endl;
      }
      else
      {
        //Bug 1239:  Commenting out random shuffle to address non-deterministic
        //           parallel behavior due to shared voltage node assignment.
        //Gitlab issue #149: more balanced assignment of voltage nodes required
        int new_pid = *(intVec.begin());
        int pid_count = procBalance[ new_pid ];
        std::vector<int>::iterator iV = intVec.begin(); iV++;
        while ( iV != intVec.end() )
        {
          if ( procBalance[ *iV ] < pid_count )
          {
            new_pid = *iV;
            pid_count = procBalance[ new_pid ];
          }
          iV++;
        }
        inode->pid = new_pid;
        procBalance[ new_pid ]++;
      }
    }

    OwnedVNodes.insert( VNodeContainer::value_type( id, inode ) );
  }

  VNodes = OwnedVNodes;
  OwnedVNodes.clear();

  //gids for both
  //locally loop over containers in VDir and DDir and set GIDs
  int VSize = VNodes.size();
  double tmpVar1, tmpVar2;

  tmpVar1 = static_cast<double>(VSize);
  comm.scanSum( &tmpVar1, &tmpVar2, 1 );
  int baseVNodeGID = static_cast<int>(tmpVar2) - VSize;

  comm.sumAll( &tmpVar1, &tmpVar2, 1 );
  int globalVNodeCnt = static_cast<int>(tmpVar2);

  currGID = baseVNodeGID;
  iterVN = VNodes.begin();
  for( ; iterVN != endVN; ++iterVN, ++currGID )
    iterVN->second->gid = currGID;

  //compute base GID for device nodes to be used later
  tmpVar1 = static_cast<double>(dSize);
  comm.scanSum( &tmpVar1, &tmpVar2, 1 );
  baseDNodeGID = static_cast<int>(tmpVar2) - dSize + globalVNodeCnt;

  //push indices back
  //do gets on VDir to get ownership and GIDs
  VData.clear();

  VDir->getEntries( VNames, VData );
  delete VDir;
  
  currGID = baseDNodeGID;
  for (CktNodeList::const_iterator it_cnL = topology_.getOrderedNodeList().begin(), end_cnL = topology_.getOrderedNodeList().end(); it_cnL != end_cnL; ++it_cnL )
  {
    CktNode & cn = **it_cnL;
    const std::string & id = cn.get_id();
    int type = cn.type();

    if( id != "0" )
    {
      if( type == _VNODE )
      {
        IndexNodePtr inode = VData[id];
        cn.set_gID( inode->gid );
        cn.set_ProcNum( inode->pid );
        cn.set_IsOwned( procID == inode->pid );
      }
      else if( type == _DNODE )
      {
        cn.set_gID( currGID++ );
        cn.set_ProcNum( procID );
        cn.set_IsOwned( true );
      }
    }
    else
    {
      cn.set_gID( -1 );
      cn.set_ProcNum( -1 );
      cn.set_IsOwned( false );
    }
  }

  if (DEBUG_TOPOLOGY)
  {
    comm.barrier();

    for( int i = 0; i < procCnt; ++i )
    {
      if( i == procID )
      {
        Xyce::dout() << "<<<<<<<<<<<<<ORIGINAL INDEX NODES>>>>>>>>>>>>>>>>>>>>>> " << procID << "\n";
        CktNodeList::const_iterator itL2 = topology_.getOrderedNodeList().begin();
        CktNodeList::const_iterator endL2 = topology_.getOrderedNodeList().end();
        for( int i = 0 ; itL2 != endL2; ++itL2, ++i )
        {
          Xyce::dout() << "Proc: " << procID << "\t" << "Node: " << i << std::endl;
          Xyce::dout() << **itL2;
        }
        Xyce::dout() << "<<<<<<<<ORIGINAL INDEX NODES END>>>>>>>>>>>>>>>>>>>>>> " << procID << "\n";
      }
      comm.barrier();
    }

    comm.barrier();
  }

  setupNodeGIDs();

  if (DEBUG_TOPOLOGY)
  {
  comm.barrier();

  for( int i = 0; i < procCnt; ++i )
  {
    if( i == procID )
    {
      Xyce::dout() << "<<<<<<<<<<<<<REINDEXED NODES>>>>>>>>>>>>>>>>>>>>>> " << procID << "\n";
      CktNodeList::const_iterator itL2 = topology_.getOrderedNodeList().begin();
      CktNodeList::const_iterator endL2 = topology_.getOrderedNodeList().end();
      for( int i = 0 ; itL2 != endL2; ++itL2, ++i )
      {
        Xyce::dout() << "Proc: " << procID << "\t" << "Node: " << i << std::endl;
        Xyce::dout() << **itL2;
      }
      Xyce::dout() << "<<<<<<<<REINDEXED NODES END>>>>>>>>>>>>>>>>>>>>>> " << procID << "\n";
    }
    comm.barrier();
  }

  comm.barrier();
  }

  setupSolnAndStateGIDs();

  if (DEBUG_TOPOLOGY)
    Xyce::dout() << topology_ << std::endl;

  if( checkConnectivity_ )
  { 
    //Setup index to GID map in CktGraph
    topology_.regenerateGIDNodeMap();
    
    testVoltageNodeConnectivity_();
  }

  //resolve late dependencies for devices
  topology_.resolveDependentVars();

  extractAllGIDsFromTopology();
  
  generateRowColData();

  topology_.writeNetlistAddResistors( connToOneTermIDs_, noDCPathIDs_ );

  //Don't remove this!!!!!!!!
  return true;
}

//-----------------------------------------------------------------------------
// Function      : ParLSUtil::testVoltageNodeConnectivity
// Purpose       : testing of voltage node connectivity for problems
// Special Notes : initially, just warn if a voltage node has only 1 connection
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/31/03
//-----------------------------------------------------------------------------
bool ParLSUtil::testVoltageNodeConnectivity_()
{
  CktNodeList::const_iterator it_cnL = topology_.getOrderedNodeList().begin();
  CktNodeList::const_iterator end_cnL = topology_.getOrderedNodeList().end();

  Parallel::Communicator & comm = *(pdsManager_.getPDSComm());
  int procCnt = comm.numProc();
  int procID = comm.procID();
  int m, n;
  int proc;
  int i, j, k;
  int gid, num_gid;
  std::vector<int> gidList;
  std::map<int, std::vector<int> > gid_map;
  std::map<int, std::vector<int> >::iterator gm_i, gm_end;

  std::map< int, std::set<int> > comm_gid;
  std::map< int, std::set<int> >::iterator cg_i, cg_end;
  std::set<int>::iterator vm_i, vm_end;

  std::vector<std::vector<int> > buf;

  it_cnL = topology_.getOrderedNodeList().begin();
  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    if ((*it_cnL)->type() == _VNODE)
    {
      gid = (*it_cnL)->get_gID();
      // If this voltage node is not owned by this processor, insert it
      // in the comm_gid list.
      if (gid >= 0 && !(*it_cnL)->get_IsOwned())
      {
        proc = (*it_cnL)->get_ProcNum();
        comm_gid[proc].insert(gid);
      }
    }
  }

  int num_proc = comm_gid.size();
  int max_num_proc = 0;
  comm.maxAll (&num_proc, &max_num_proc, 1);

  int tmpSize = procCnt*max_num_proc;
  std::vector<int> pcomm    (tmpSize,0);
  std::vector<int> pcomm_all(tmpSize,0);
  std::vector<int> ncomm    (tmpSize,0);
  std::vector<int> ncomm_all(tmpSize,0);

  i = 0;
  cg_i   = comm_gid.begin();
  cg_end = comm_gid.end() ;
  for (; cg_i!=cg_end; ++cg_i)
  {
    pcomm[procID*max_num_proc+i] = (*cg_i).first+1;
    ncomm[procID*max_num_proc+i] = (*cg_i).second.size();
    i++;
  }
  comm.sumAll(&pcomm[0], &pcomm_all[0], tmpSize);
  comm.sumAll(&ncomm[0], &ncomm_all[0], tmpSize);

  std::vector<int> sendBuf;
  std::vector<int> src;
  int numMsg = 0;

  // Allocate receive buffers
  for (i=0 ; i<procCnt*max_num_proc ; ++i)
  {
    if (pcomm_all[i] == procID+1)
    {
      j = i/max_num_proc;
      k = ncomm_all[i];
      buf.resize(numMsg+1);
      buf[numMsg].resize(k);
      src.push_back(j);
      numMsg++;
    }
  }

  // Issue receives
  for (i=0 ; i<numMsg ; ++i)
  {
    comm.iRecv (&buf[i][0], buf[i].size(), src[i]);
  }

  // Issue sends
  cg_i=comm_gid.begin();
  for ( ; cg_i!=cg_end; ++cg_i)
  {
    i = (*cg_i).second.size();
    sendBuf.resize(i);
    j = 0;
    vm_i=(*cg_i).second.begin();
    vm_end=(*cg_i).second.end();
    for ( ; vm_i !=vm_end; ++vm_i)
    {
      sendBuf[j++] = *vm_i;
    }
    comm.send (&sendBuf[0], i, (*cg_i).first);
  }

  // Wait for messages to complete
  comm.waitAll();

  // Use received data
  for (i=0 ; i<numMsg ; ++i)
  {
    k = buf[i].size();
    for (m=0 ; m<k ; ++m)
      comm_gid[src[i]].insert(buf[i][m]);
  }

  int n_bufs = comm_gid.size();
  std::vector<int> buf_dest(n_bufs), buf_len(n_bufs);
  std::vector<int *> buf_in(n_bufs), buf_out(n_bufs);

  int buf_tot = 0;
  cg_i   = comm_gid.begin();
  cg_end = comm_gid.end();
  for ( ; cg_i!=cg_end; ++cg_i)
  {
    buf_tot += (*cg_i).second.size();
  }

  std::vector<int> actual_buf_in(buf_tot,0);
  std::vector<int> actual_buf_out(buf_tot,0);

  i = 0;
  m = 0;
  n = 0;
  cg_i   = comm_gid.begin();
  cg_end = comm_gid.end();
  for ( ; cg_i!=cg_end; ++cg_i)
  {
    j = (*cg_i).second.size();
    buf_dest[i] = (*cg_i).first;
    buf_len[i] = j;
    buf_in[i] = &actual_buf_in[m];
    buf_out[i] = &actual_buf_out[m];
    m += j;

    k = 0;
    vm_i   = (*cg_i).second.begin();
    vm_end = (*cg_i).second.end();
    for ( ; vm_i !=vm_end; ++vm_i)
    {
      buf_in[i][k] = 0;
      buf_out[i][k] = 0;
      gid_map[*vm_i].push_back(n+k);
      ++k;
    }
    n += k;
    ++i;
  }

  // Determine which local voltage nodes are only connected to 1 device terminal.
  it_cnL = topology_.getOrderedNodeList().begin();
  for( ; it_cnL != end_cnL; ++it_cnL )
  { 
    if ((*it_cnL)->type() == _VNODE)
    { 
      gid = (*it_cnL)->get_gID();
      if (gid >= 0)
      { 
        num_gid = topology_.numAdjNodes( gid );
        if (num_gid > 0 && !((*it_cnL)->get_IsOwned()))
        { 
          gm_i = gid_map.find(gid);
          if (gm_i != gid_map.end())
          { 
            actual_buf_out[(gm_i->second)[0]] = num_gid;
          }
        }
      }
    }
  }

  comm_boundaries (gid_map, actual_buf_in, actual_buf_out,
                   buf_len, buf_dest, buf_in, buf_out, 1);

  // The gid_pos is necessary because the GIDs in the CktNode objects
  // include BOTH device and voltage nodes.  This is an easier way to
  // use a contiguous container (std::vector) to perform the node checking
  // operations.
  int num_nodes = 0;
  unordered_map<int, int> gid_pos;
  it_cnL = topology_.getOrderedNodeList().begin();
  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    if ((*it_cnL)->type() == _VNODE)
    {
      gid = (*it_cnL)->get_gID();
      if (gid >= 0)
      {
        gid_pos[gid] = num_nodes++;
        num_gid = topology_.numAdjNodes( gid );
        if (num_gid <= 1 && (*it_cnL)->get_IsOwned())
        {
          gm_i = gid_map.find(gid);
          if (gm_i == gid_map.end() ||
              (gm_i != gid_map.end() &&
               num_gid + actual_buf_in[(gm_i->second)[0]] <= 1))
          {
            connToOneTermIDs_.insert( (*it_cnL)->get_id() );
          }
        }
      }
    }
  }

  unordered_map<int, int>::iterator gp_i;
  unordered_map<int, int>::iterator gp_end;

  // Collect the GIDs for the external (ghost) nodes for this processor.
  std::vector<int> ext_gid  (gid_map.size(),0);
  gm_i   = gid_map.begin();
  gm_end = gid_map.end();
  for (i=0 ; gm_i != gm_end; ++gm_i, ++i)
  {
    gp_i = gid_pos.find((*gm_i).first);
    if (gp_i != gid_pos.end())
    {
      ext_gid[i] = gp_i->second;
    }
    else
    {
      Report::DevelFatal0()
        << "ParLSUtil::testVoltageNodeConnectivity_: External GID not found internally";
    }
  }

  // Go through all the devices and color the voltage nodes by their lead groupings.
  // This will induce a coloring of the voltage nodes on this processor, although it is possible
  // to have voltage nodes with no colors if they are ghost nodes.
  it_cnL = topology_.getOrderedNodeList().begin();
  int num_colors = 0;
  std::vector<int> node_color(num_nodes,0);  // If node_color == 0, then the node has not been colored yet.
  unordered_set<int> grounded_colors;        // Collect nodes that are grounded while processing DNODEs.
  std::map<int, std::vector<int> > color_pos;

  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    if ((*it_cnL)->type() == _DNODE)
    {
      CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_cnL);

      // Get the list of connected leads from the device
      const std::vector<int> & lead_conn = cktNodeDevPtr->leadConnect();
      i = lead_conn.size();

      if (i > 0)
      {
        // Get the adjacent voltage GIDs to this device GID that include ground
        int gid = (*it_cnL)->get_gID();
        int max_lead_value = 0;
        gidList.clear();
        topology_.returnAdjGIDsWithGround(gid, gidList);

        // Go through the groupings for the external nodes determined by lead_conn.
        // if lead_conn[k]==0, then the node is always grounded, so set the node_color to -1.
        std::vector<int>::const_iterator lead_it = lead_conn.begin();
        std::vector<int>::const_iterator lead_end = lead_conn.end();
        std::vector<int>::const_iterator gid_it = gidList.begin();
        std::vector<int>::const_iterator gid_end = gidList.end();
        for ( ; lead_it != lead_end; lead_it++, gid_it++ )
        {
          if (*lead_it == 0)
          {
            if (*gid_it >= 0)
            {
              int pos = gid_pos[*gid_it];
              node_color[pos] = -1;
            }
            i--; // Don't need to account for this node, it's grounded.
          }
          else if (*lead_it > max_lead_value)
          {
              max_lead_value = *lead_it;
          }
        }
        for (j=1 ; j<=max_lead_value; ++j)
        {
          int gid_min = num_colors+1;
          std::vector<int> gid_pos_set;
          unordered_set<int> node_color_set;
          lead_it = lead_conn.begin();
          lead_end = lead_conn.end();
          gid_it = gidList.begin();
          gid_end = gidList.end();
          for ( ; lead_it != lead_end; lead_it++, gid_it++ )
          {
            if (*lead_it == j)
            {
              i--; // Decrement node counter, this voltage node is being processed.

              if (*gid_it >= 0)
              {
                int pos = gid_pos[*gid_it];
                int curr_node_color = node_color[pos];

                // If the node color is not 0, check if it should be the minimum node color.
                // This logic will allow the minimum node color to be -1 or any positive number.
                if (curr_node_color)
                {
                  if (curr_node_color < gid_min)
                    gid_min = curr_node_color;
                  if (curr_node_color != -1)
                    node_color_set.insert( curr_node_color );
                }
                gid_pos_set.push_back( pos );
              }
              else
              {
                gid_min = -1;
              }
            }
          }

          // For the case when none of the nodes in the lead group are colored.
          if (gid_min == (num_colors+1) && gid_pos_set.size() > 1)
            gid_min = ++num_colors;

          // One of the nodes in this lead group is grounded, ground the colors
          // attached to this node or create a new color and ground it.
          if (gid_min == -1 && gid_pos_set.size())
          {
            std::vector<int>::iterator pos_it = gid_pos_set.begin();
            std::vector<int>::iterator pos_end = gid_pos_set.end();
            bool add_color = false;
            int new_color = num_colors+1;
            for( ; pos_it != pos_end; pos_it++ )
            {
              if (!node_color[*pos_it])
              {
                node_color[*pos_it] = new_color;
                color_pos[new_color].push_back( *pos_it );
                add_color = true;
              }
            }
            if (add_color)
            { 
              ++num_colors;  // Increment num_colors if we used a new color here.
              node_color_set.insert( new_color );
            }
            unordered_set<int>::iterator color_it = node_color_set.begin();
            unordered_set<int>::iterator color_end = node_color_set.end();
            for ( ; color_it != color_end; color_it++ )
            {
              grounded_colors.insert( *color_it ); 
            }
          }

          // We have a positive gid_min to assign to the GIDs of this lead group.
          // If there are other node colors for other GIDs in the lead group, reassign them to gid_min.
          if (gid_min != -1 && gid_pos_set.size() > 1)
          {
            std::vector<int>::iterator pos_it = gid_pos_set.begin();
            std::vector<int>::iterator pos_end = gid_pos_set.end();
            for( ; pos_it != pos_end; pos_it++ )
            {
              node_color[*pos_it] = gid_min;
              color_pos[gid_min].push_back( *pos_it );            
            }
            // If there are other colors that were assigned to this new color, consolidate them. 
            if (node_color_set.size() > 1)
            {
              unordered_set<int>::iterator color_it = node_color_set.begin();
              unordered_set<int>::iterator color_end = node_color_set.end();
              for ( ; color_it != color_end; color_it++ )
              {
                if ( *color_it != gid_min )
                {
                  // Copy all the GIDs from the additional color into the set of GIDs from gid_min.
                  // Change the node_color to gid_min and then clear the set of GIDs for color_it.
                  std::vector<int>::iterator add_it = color_pos[*color_it].begin();
                  std::vector<int>::iterator add_end = color_pos[*color_it].end();
                  std::copy( add_it, add_end, std::back_inserter(color_pos[gid_min]) );
                  for ( ; add_it != add_end; add_it++ )
                  {
                    node_color[*add_it] = gid_min;
                  }
                  color_pos.erase( *color_it );

                  // If the color being consolidated was grounded, add gid_min to the list of grounded colors.
                  if ( grounded_colors.count( *color_it ) )
                  {
                    grounded_colors.insert( gid_min );
                    grounded_colors.erase( *color_it );
                  }
                } 
              }
            }               
          } 
          if (i == 0)
          {
            break;
          }
        }
        if (i != 0)
        {
          Report::DevelFatal0() << " Connectivity checker error: run Xyce with this diagnostic turned off, via:\n  .options topology CHECK_CONNECTIVITY=0 ";
        }
      } // if (i > 0)
    } // if ((*it_cnL)->type() == _DNODE)
  } // for( ; it_cnL != end_cnL; ++it_cnL )

  // Now process any color that has been determined to be grounded already.
  unordered_set<int>::iterator ground_it = grounded_colors.begin();
  unordered_set<int>::iterator ground_end = grounded_colors.end();
  for ( ; ground_it != ground_end; ground_it++ )
  {
    std::vector<int>::iterator add_it = color_pos[*ground_it].begin();
    std::vector<int>::iterator add_end = color_pos[*ground_it].end();
    for ( ; add_it != add_end; add_it++ )
    {
      node_color[*add_it] = -1;
    }
    color_pos.erase( *ground_it );
  }

  // Create a copy of the external node pos vector, then sort it.
  std::vector<int> ext_pos = ext_gid;
  std::sort( ext_pos.begin(), ext_pos.end() );

  // Now go through the current node_color list to categorize each remaining node by color, single, or ungrounded.
  // NOTE:  Singles are the result of ghost nodes, so they should all be in the ext_gid vector.
  std::map<int,int> ext_pos_color;
  std::set<int> ext_pos_single;
  std::vector<int> ungrounded_pos;
  std::vector<int>::iterator node_color_it = node_color.begin();
  std::vector<int>::iterator node_color_end = node_color.end();
  for ( int i=0 ; node_color_it != node_color_end; i++, node_color_it++ )
  {
    if ( *node_color_it != -1 )
    {
      if ( *node_color_it == 0 )
      {
        // If this node is not in the ext_gid vector, then this node is ungrounded.
        std::vector<int>::iterator found_it = std::find( ext_pos.begin(), ext_pos.end(), i );
        if (found_it != ext_pos.end())
          ext_pos_single.insert( i );
        else
          ungrounded_pos.push_back( i );
      }
      else
      {
        ext_pos_color[ i ] = *node_color_it;
      }
    } 
  }

  bool work_to_do = true;

  while (work_to_do)
  {
    // Send / receive boundary information
    // NOTE:  The external gids (input) are only updated if not all of the local nodes are grounded.
    gm_i   = gid_map.begin();
    gm_end = gid_map.end();
    for (i=0 ; gm_i != gm_end; ++gm_i, ++i)
    {
      actual_buf_out[(*gm_i).second[0]] = node_color[ext_gid[i]];
    }
    comm_boundaries
      (gid_map, actual_buf_in, actual_buf_out,
       buf_len, buf_dest, buf_in, buf_out, 2);

    work_to_do = false;

    // This processor still might have work to do if there are ungrounded colors or singles.
    if (color_pos.size() || ext_pos_single.size())
    {
      gm_i   = gid_map.begin();
      gm_end = gid_map.end();
      //  Check if some processor sent a NEW grounded node to this processor.
      for (i=0 ; gm_i != gm_end; ++gm_i, ++i)
      {
        int new_node_color = actual_buf_in[(*gm_i).second[0]];
        if ( (node_color[ext_gid[i]] > new_node_color) && (new_node_color == -1) )
        {
          work_to_do = true;
          std::set<int>::iterator sfound_it = ext_pos_single.find( ext_gid[i] );
          if (sfound_it != ext_pos_single.end())
          {
            node_color[ext_gid[i]] = -1;
            ext_pos_single.erase( sfound_it );
          }
          else
          {
            // If this GID is part of a color, then set the color to -1.
            std::map<int,int>::iterator found_it = ext_pos_color.find(ext_gid[i]);
            if (found_it != ext_pos_color.end())
            {
              std::vector<int>::iterator add_it = color_pos[found_it->second].begin();
              std::vector<int>::iterator add_end = color_pos[found_it->second].end();
              for ( ; add_it != add_end; add_it++ )
              {
                node_color[*add_it] = -1;
              }
              color_pos.erase( found_it->second );
              ext_pos_color.erase( found_it );
            }
          }      
        }
      }
    } 

    // Check if every processor is finished. 
    j = work_to_do ? 1 : 0;
    comm.sumAll(&j, &k, 1);
    work_to_do = (k > 0);
  }

  // If any colors (groups), single nodes, or ungrounded nodes have entries, then there are ungrounded nodes.
  // This will be indicated by any node_color that is not set to -1. 
  // NOTE:  Now that coloring is being used, it is possible to discriminate if there is a group of nodes that
  // require a DC path to ground.  That would reduce the number of resistors that would need to be added.
  if (color_pos.size() || ext_pos_single.size() || ungrounded_pos.size())
  {
    it_cnL = topology_.getOrderedNodeList().begin();
    for( ; it_cnL != end_cnL; ++it_cnL )
    {
      if ((*it_cnL)->type() == _VNODE)
      {
        gid = (*it_cnL)->get_gID();
        if (gid >= 0 && node_color[gid_pos[gid]] >= 0)
        {
          noDCPathIDs_.insert( (*it_cnL)->get_id() );  
        }                                            
      }                                              
    }                                                
  }

  // Output the warnings through the topology object.
  topology_.outputTopoWarnings(connToOneTermIDs_, noDCPathIDs_);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ParLSUtil::comm_boundaries
// Purpose       : communicate boundary node data for topological checks
// Special Notes : mode=1 is sum, mode=2 is min
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/05
//-----------------------------------------------------------------------------
void ParLSUtil::comm_boundaries (std::map<int, std::vector<int> > & gid_map,
                                 std::vector<int> & actual_buf_in, std::vector<int> & actual_buf_out,
                                 std::vector<int> & buf_len, std::vector<int> & buf_dest,
                                 std::vector<int *> & buf_in, std::vector<int *> & buf_out, int mode)

{
  Parallel::Communicator & comm = *(pdsManager_.getPDSComm());
  unsigned int i;
  unsigned int n_bufs = buf_len.size();
  std::map< int, std::map<int, bool> >::iterator cg_i;
  std::map<int, std::vector<int> >::iterator g_i = gid_map.begin() ;
  std::map<int, std::vector<int> >::iterator g_end = gid_map.end();
  for ( ; g_i != g_end; ++g_i)
  {
    if ((*g_i).second.size() > 1)
    {
      for (i=1 ; i<(*g_i).second.size() ; ++i)
        actual_buf_out[(*g_i).second[i]] = actual_buf_out[(*g_i).second[0]];
    }
  }

  for (i = 0 ; i < n_bufs ; ++i)
  {
    comm.iRecv (buf_in[i], buf_len[i], buf_dest[i]);
  }
  for (i = 0 ; i < n_bufs ; ++i)
  {
    comm.send (buf_out[i], buf_len[i], buf_dest[i]);
  }
  comm.waitAll();

  g_i = gid_map.begin();
  g_end = gid_map.end();
  for ( ; g_i != g_end; ++g_i)
  {
    if ((*g_i).second.size() > 1)
    {
      for (i=1 ; i<(*g_i).second.size() ; ++i)
      {
        if (mode == 1)
          actual_buf_in[(*g_i).second[0]] += actual_buf_in[(*g_i).second[i]];
        else if (mode == 2)
        {
          if (actual_buf_in[(*g_i).second[i]] < actual_buf_in[(*g_i).second[0]])
            actual_buf_in[(*g_i).second[0]] = actual_buf_in[(*g_i).second[i]];
        }
      }
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ParLSUtil::setupNodeGIDs
// Purpose       : Generate Ordering and Var GIDs.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/22/01
//-----------------------------------------------------------------------------
bool ParLSUtil::setupNodeGIDs()
{
  //get lists of owned and boundary/ghost node GIDs
  topology_.generateOrderedNodeList();

  extractNodeGIDs();

  numLocalNodes_ = nodeList_GID_.size();

  //calculate base GID for this processors nodes, lex. ordering
  double tmpVar1, tmpVar2;
  tmpVar1 = numLocalNodes_;
  pdsManager_.getPDSComm()->scanSum( &tmpVar1, &tmpVar2, 1 );
  baseNodeGID_ = static_cast<int>(tmpVar2 - numLocalNodes_);
  pdsManager_.getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalNodes_ = static_cast<int>(tmpVar2);

  //Setup temp global accessors
  //Use to get new node GIDs for boundary/ghost nodes
  //---------------------------
  Parallel::GlobalAccessor * nodeGlobalAccessorPtr = pdsManager_.createGlobalAccessor();
  nodeGlobalAccessorPtr->registerExternGIDVector( nodeList_ExternGID_ );
  nodeGlobalAccessorPtr->generateMigrationPlan();

  std::map<int,int> nodeGIDMap;
  std::vector<int>::const_iterator it_nL = nodeList_GID_.begin();
  for( int iSV = 0; it_nL != nodeList_GID_.end(); ++it_nL, ++iSV )
    nodeGIDMap[ *it_nL ] = iSV + baseNodeGID_;

  std::map<int,int> externGIDMap;
  nodeGlobalAccessorPtr->migrateIntArray( nodeGIDMap, externGIDMap );
  delete nodeGlobalAccessorPtr;

  //loop over nodes and reset all GIDs to new lex. ordering
  for (CktNodeList::const_iterator it_cnL = topology_.getOrderedNodeList().begin(), end_cnL = topology_.getOrderedNodeList().end(); it_cnL != end_cnL; ++it_cnL )
  {
    std::map<int,int>::const_iterator gotGID = nodeGIDMap.find( (*it_cnL)->get_gID() );
    if( gotGID != nodeGIDMap.end() )
    {
      (*it_cnL)->set_gID( gotGID->second );

      if (DEBUG_TOPOLOGY)
        Xyce::dout() << "Node: " << (*it_cnL)->get_id() << " " << (*it_cnL)->get_gID() << std::endl;

    }
    else
    {
      std::map<int,int>::const_iterator gotExtGID = externGIDMap.find( (*it_cnL)->get_gID() );

      if( gotExtGID != externGIDMap.end() )
      {
        (*it_cnL)->set_gID( gotExtGID->second );
      }
      else if( (*it_cnL)->get_gID() != -1 )
      {
        Report::DevelFatal0() << "P" <<  Teuchos::Utils::toString(pdsManager_.getPDSComm()->procID())
                              << ": Node: " << (*it_cnL)->get_id() << ", global index ("
                              <<  Teuchos::Utils::toString( (*it_cnL)->get_gID() ) << ") is NOT found!";
      }
    }
  }

  //Reset global accessor for new lex. indexing
  topology_.generateOrderedNodeList();

  extractNodeGIDs();

  return true;
}


//-----------------------------------------------------------------------------
// Function      : ParLSUtil::setupSolnAndStateGIDs
// Purpose       : Generate Ordering and Var GIDs.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/2/01
//-----------------------------------------------------------------------------
bool ParLSUtil::setupSolnAndStateGIDs()
{
  //get soln and state var counts from all owned v-nodes and d-nodes
  std::vector<int> rowCountVec(numLocalNodes_);
  std::vector<int> stateCountVec(numLocalNodes_,0);
  std::vector<int> storeCountVec(numLocalNodes_,0);
  std::vector<int> leadCurrentCountVec(numLocalNodes_,0);

  double tmpVar1, tmpVar2;
  int Loc = 0;
  numLocalRows_ = 0;
  numLocalStateVars_ = 0;
  numLocalStoreVars_ = 0;
  numLocalLeadCurrentVars_ = 0;
  for (CktNodeList::const_iterator it_cnL = topology_.getOrderedNodeList().begin(), end_cnL = topology_.getOrderedNodeList().end(); it_cnL != end_cnL; ++it_cnL )
    if( (*it_cnL)->get_IsOwned() && (*it_cnL)->get_gID() != -1 )
    {
      rowCountVec[Loc] = (*it_cnL)->solnVarCount();
      numLocalRows_ += rowCountVec[Loc];

      if( (*it_cnL)->type() == _DNODE )
      {
         CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_cnL);
         stateCountVec[Loc] = cktNodeDevPtr->stateVarCount();
         numLocalStateVars_ += stateCountVec[Loc];
         storeCountVec[Loc] = cktNodeDevPtr->storeVarCount();
         numLocalStoreVars_ += storeCountVec[Loc];
         leadCurrentCountVec[Loc] = cktNodeDevPtr->branchDataVarCount();
         numLocalLeadCurrentVars_ += leadCurrentCountVec[Loc];
      }
      ++Loc;
    }

  //calculate base soln and state GIDs for this processor
  tmpVar1 = numLocalRows_;
  pdsManager_.getPDSComm()->scanSum( &tmpVar1, &tmpVar2, 1 );
  baseRowGID_ = static_cast<int>(tmpVar2 - numLocalRows_);
  pdsManager_.getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalRows_ = static_cast<int>(tmpVar2);

  tmpVar1 = numLocalStateVars_;
  pdsManager_.getPDSComm()->scanSum( &tmpVar1, &tmpVar2, 1 );
  baseStateVarGID_ = static_cast<int>(tmpVar2 - numLocalStateVars_);
  pdsManager_.getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalStateVars_ = static_cast<int>(tmpVar2);

  tmpVar1 = numLocalStoreVars_;
  pdsManager_.getPDSComm()->scanSum( &tmpVar1, &tmpVar2, 1 );
  baseStoreVarGID_ = static_cast<int>(tmpVar2 - numLocalStoreVars_);
  pdsManager_.getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalStoreVars_ = static_cast<int>(tmpVar2);

  tmpVar1 = numLocalLeadCurrentVars_;
  pdsManager_.getPDSComm()->scanSum( &tmpVar1, &tmpVar2, 1 );
  baseLeadCurrentVarGID_ = static_cast<int>(tmpVar2 - numLocalLeadCurrentVars_);
  pdsManager_.getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalLeadCurrentVars_ = static_cast<int>(tmpVar2);

  //Use Global Accessors to get GIDs for boundary/ghost nodes
  std::map< int,std::vector<int> > rowGIDMap, stateGIDMap;
  std::map< int,std::vector<int> > externRowGIDMap, externStateGIDMap;
  std::map< int,std::vector<int> > storeGIDMap;
  std::map< int,std::vector<int> > externStoreGIDMap;
  std::map< int,std::vector<int> > leadCurrentGIDMap;
  std::map< int,std::vector<int> > externLeadCurrentGIDMap;

  int currRowLoc = baseRowGID_;
  int currStateLoc = baseStateVarGID_;
  int currStoreLoc = baseStoreVarGID_;
  int currLeadCurrentLoc = baseLeadCurrentVarGID_;

  for( int i = 0; i < numLocalNodes_; ++i )
  {
    std::vector<int>& rowGIDvec = rowGIDMap[ nodeList_GID_[i] ];
    rowGIDvec.resize(rowCountVec[i]);
    iota( rowGIDvec.begin(), rowGIDvec.end(), currRowLoc );
    currRowLoc += rowCountVec[i];

    std::vector<int>& stateGIDvec = stateGIDMap[ nodeList_GID_[i] ];
    stateGIDvec.resize(stateCountVec[i]);
    iota( stateGIDvec.begin(), stateGIDvec.end(), currStateLoc );
    currStateLoc += stateCountVec[i];

    std::vector<int>& storeGIDvec = storeGIDMap[ nodeList_GID_[i] ];
    storeGIDvec.resize(storeCountVec[i]);
    iota( storeGIDvec.begin(), storeGIDvec.end(), currStoreLoc );
    currStoreLoc += storeCountVec[i];
   
    std::vector<int>& leadCurrGIDvec = leadCurrentGIDMap[ nodeList_GID_[i] ]; 
    leadCurrGIDvec.resize(leadCurrentCountVec[i]);
    iota( leadCurrGIDvec.begin(), leadCurrGIDvec.end(), currLeadCurrentLoc );
    currLeadCurrentLoc += leadCurrentCountVec[i];
  }

  Parallel::GlobalAccessor * nodeGlobalAccessorPtr = pdsManager_.createGlobalAccessor();
  nodeGlobalAccessorPtr->registerExternGIDVector( nodeList_ExternGID_ );
  nodeGlobalAccessorPtr->generateMigrationPlan();

  nodeGlobalAccessorPtr->migrateIntVecs( rowGIDMap, externRowGIDMap );
  nodeGlobalAccessorPtr->migrateIntVecs( stateGIDMap, externStateGIDMap );
  nodeGlobalAccessorPtr->migrateIntVecs( storeGIDMap, externStoreGIDMap );
  nodeGlobalAccessorPtr->migrateIntVecs( leadCurrentGIDMap, externLeadCurrentGIDMap );
  delete nodeGlobalAccessorPtr;

  //calculate number of boundary/ghost soln and state vars
  std::map< int,std::vector<int> >::iterator iterIVM = externRowGIDMap.begin();
  std::map< int,std::vector<int> >::iterator endIVM = externRowGIDMap.end();
  numExternRows_ = 0;
  for( ; iterIVM != endIVM; ++iterIVM )
    numExternRows_ += iterIVM->second.size();
  tmpVar1 = numExternRows_;
  pdsManager_.getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalExternRows_ = static_cast<int>(tmpVar2);

  //loop over nodes and assign soln and state GIDs
  for (CktNodeList::const_iterator it_cnL = topology_.getOrderedNodeList().begin(), end_cnL = topology_.getOrderedNodeList().end(); it_cnL != end_cnL; ++it_cnL )
  {
    if( rowGIDMap.count( (*it_cnL)->get_gID() ) )
    {
      (*it_cnL)->set_SolnVarGIDList( rowGIDMap[ (*it_cnL)->get_gID() ] );
      if( (*it_cnL)->type() == _DNODE )
      {
        CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_cnL);
        cktNodeDevPtr->set_StateVarGIDList( stateGIDMap[ (*it_cnL)->get_gID() ] );
        cktNodeDevPtr->set_StoreVarGIDList( storeGIDMap[ (*it_cnL)->get_gID() ] );
        cktNodeDevPtr->set_LeadCurrentVarGIDList( leadCurrentGIDMap[ (*it_cnL)->get_gID() ] );
      }
    }
    else if( externRowGIDMap.count( (*it_cnL)->get_gID() ) )
    {
      (*it_cnL)->set_SolnVarGIDList( externRowGIDMap[ (*it_cnL)->get_gID() ] );
    }
    else if( (*it_cnL)->get_gID() == -1 )
      (*it_cnL)->set_SolnVarGIDList( std::vector<int>(1,-1) );
    else
    {
      Report::DevelFatal0() << "P" <<  Teuchos::Utils::toString(pdsManager_.getPDSComm()->procID())
                            << ": Node: " << (*it_cnL)->get_id() << ", global index ("
                            << Teuchos::Utils::toString( (*it_cnL)->get_gID() ) << ") is NOT found!";
    }
  }

  topology_.registerGIDs();

  return true;

}

//-----------------------------------------------------------------------------
// Function      : ParLSUtil::generateRowColData
// Purpose       : Generate row/col data for linear solver.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/26/00
//-----------------------------------------------------------------------------
void ParLSUtil::generateRowColData()
{
  int procID = pdsManager_.getPDSComm()->procID();

  //--- add in dep soln var stuff
  if( !topology_.getDepSolnGIDMap().empty() )
  {
    std::map<int,int> tmpMap;
    std::map<int,int>::const_iterator iterIIM = topology_.getDepSolnGIDMap().begin();
    std::map<int,int>::const_iterator endIIM = topology_.getDepSolnGIDMap().end();
    for( ; iterIIM != endIIM; ++iterIIM )
      if( iterIIM->second != procID ) tmpMap.insert( *iterIIM );

    for( unsigned int i = 0; i < rowList_ExternGID_.size(); ++i )
      if( rowList_ExternGID_[i].second  != procID )
        tmpMap.insert( rowList_ExternGID_[i] );

    rowList_ExternGID_.resize( tmpMap.size() );
    iterIIM = tmpMap.begin();
    for( unsigned int i = 0; i < tmpMap.size(); ++iterIIM, ++i )
      rowList_ExternGID_[i] = *iterIIM;
  }

  //--- set numLocalRows_, numLocalStateVars_, and resize lists
  numLocalRows_ = rowList_GID_.size();
  numExternRows_ = rowList_ExternGID_.size();
  numLocalStateVars_ = rowList_StateGID_.size();
  numLocalStoreVars_ = rowList_StoreGID_.size();
  numLocalLeadCurrentVars_ = rowList_LeadCurrentGID_.size();

  double tmpVar1, tmpVar2;

  tmpVar1 = numLocalRows_;
  pdsManager_.getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalRows_ = static_cast<int>(tmpVar2);

  tmpVar1 = numExternRows_;
  pdsManager_.getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalExternRows_ = static_cast<int>(tmpVar2);

  tmpVar1 = numLocalStateVars_;
  pdsManager_.getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalStateVars_ = static_cast<int>(tmpVar2);

  tmpVar1 = numLocalStoreVars_;
  pdsManager_.getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalStoreVars_ = static_cast<int>(tmpVar2);

  tmpVar1 = numLocalLeadCurrentVars_;
  pdsManager_.getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalLeadCurrentVars_ = static_cast<int>(tmpVar2);

  //--- block extern gids by processor to support aztec requirements
  if( numExternRows_ )
  {
    std::vector<int> externGIDs( numExternRows_ );
    std::vector<int> PIDs( numExternRows_ );
    for( int i = 0; i < numExternRows_; ++i )
    {
      externGIDs[i] = rowList_ExternGID_[i].first;
      PIDs[i] = rowList_ExternGID_[i].second;
    }

    Epetra_Util Util;
    int ** listPtr = new int *[1];
    listPtr[0] = &externGIDs[0];
    Util.Sort( true, numExternRows_, &PIDs[0], 0, 0, 1, listPtr );
    delete [] listPtr;

    for( int i = 0; i < numExternRows_; ++i )
      rowList_ExternGID_[i] = std::pair<int,int>( externGIDs[i], PIDs[i] );
  }

  if (DEBUG_TOPOLOGY)
    Xyce::dout() << "Parallel Topology Util Vals: " << std::endl
                 << numGlobalRows_ << " " << numLocalRows_ << std::endl
                 << numGlobalStateVars_ << " " << numLocalStateVars_ << std::endl
                 << numGlobalStoreVars_ << " " << numLocalStoreVars_ << std::endl
                 << numGlobalLeadCurrentVars_ << " " << numLocalLeadCurrentVars_ << std::endl;

  int numRows = numLocalRows_ + numExternRows_ + 1;

  //Build global to local map for speed
  std::map<int,int> GtoL_Map;
  for( int i = 0; i < numLocalRows_; ++i )
    GtoL_Map[ rowList_GID_[i] ] = i;
  for( int i = 0; i < numExternRows_; ++i )
    GtoL_Map[ rowList_ExternGID_[i].first ] = i + numLocalRows_;
  GtoL_Map[ -1 ] = numLocalRows_ + numExternRows_;

  rowList_ColList_.resize( numRows );
  rowList_NumNZs_.resize( numRows );

  int rcCnt = rowList_ColList_.size();
  for( int i = 0; i < rcCnt; ++i ) rowList_ColList_[i].clear();

  //    to compile a list of col indices for each row
  for (CktNodeList::const_iterator it_cnL = topology_.getOrderedNodeList().begin(), end_cnL = topology_.getOrderedNodeList().end(); it_cnL != end_cnL; ++it_cnL )
  {
    if( (*it_cnL)->type() == _DNODE )
    {
      CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_cnL);

      const std::vector< std::vector<int> > & stamp = cktNodeDevPtr->jacobianStamp();
      const std::vector<int> & intGIDs = cktNodeDevPtr->get_SolnVarGIDList();
      const std::vector<int> & extGIDs = cktNodeDevPtr->get_ExtSolnVarGIDList();
      const std::vector<int> & depGIDs = cktNodeDevPtr->get_DepSolnGIDJacVec();
      std::vector<int> gids( intGIDs.size() + extGIDs.size() + depGIDs.size() );
      copy( extGIDs.begin(), extGIDs.end(), gids.begin() );
      copy( intGIDs.begin(), intGIDs.end(), gids.begin() + extGIDs.size() );
      copy( depGIDs.begin(), depGIDs.end(), gids.begin() + extGIDs.size() + intGIDs.size() );


      for( unsigned int i = 0; i < stamp.size(); ++i )
      {
        int length = stamp[i].size();
        for( int j = 0; j < length; ++j )
        {
          if (stamp[i][j] < (int)gids.size())
            rowList_ColList_[ GtoL_Map[ gids[i] ] ].push_back( gids[ stamp[i][j] ] );
        }
      }
    }
  }


  if (DEBUG_TOPOLOGY)
    Xyce::dout() << Xyce::section_divider << std::endl
                 << "Row List of NZ Cols extracted" << std::endl
                 << Xyce::section_divider << std::endl;

  numLocalNZs_ = 0;

#ifdef Xyce_NOX_LOCA_ARTIFICIAL_HOMOTOPY_SUPPORT
    //add in diagonal for homotopy support
    for( int i = 0; i < numLocalRows_; ++i )
      rowList_ColList_[i].push_back( rowList_GID_[i] );
    for( int i = 0; i < numExternRows_; ++i )
      rowList_ColList_[i+numLocalRows_].push_back( rowList_ExternGID_[i].first );
#endif

  //--- sort list of col indices and get rid of redundancies
  //    generate num of NZs data
  for( int i = 0; i < numRows; ++i )
  {
    std::sort(rowList_ColList_[i].begin(), rowList_ColList_[i].end());
    rowList_ColList_[i].erase(std::unique(rowList_ColList_[i].begin(), rowList_ColList_[i].end()), rowList_ColList_[i].end());

    rowList_NumNZs_[i] = rowList_ColList_[i].size();

    numLocalNZs_ += rowList_NumNZs_[i];
  }

  if (DEBUG_TOPOLOGY)
    Xyce::dout() << Xyce::section_divider << std::endl
                 << "Row List of NZ Cols cleaned up" << std::endl
                 << Xyce::section_divider << std::endl
                 << *this;

}

//-----------------------------------------------------------------------------
// Function      : ParLSUtil::extractNodeGIDs
// Purpose       : Fill all GID vectors using ordered node list.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/26/00
//-----------------------------------------------------------------------------
void ParLSUtil::extractNodeGIDs()
{
  nodeList_GID_.clear();
  nodeList_ExternGID_.clear();

  CktNodeList::const_iterator it = topology_.getOrderedNodeList().begin();
  CktNodeList::const_iterator end = topology_.getOrderedNodeList().end();

  for ( ; it != end; ++it )
  {
    if ( (*it)->get_gID() != -1 )
    {
      if( (*it)->get_IsOwned() )
      {
        nodeList_GID_.push_back( (*it)->get_gID() );
      }
      else
      {
        nodeList_ExternGID_.push_back( std::pair<int,int>( (*it)->get_gID(), (*it)->get_ProcNum() ) );
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : ParLSUtil::extractAllGIDsFromTopology
// Purpose       : Fill all GID vectors using ordered node list.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/26/00
//-----------------------------------------------------------------------------
void ParLSUtil::extractAllGIDsFromTopology()
{
  // Clear all vectors.
  rowList_GID_.clear();
  rowList_ExternGID_.clear();
  rowList_StateGID_.clear();
  rowList_StoreGID_.clear();
  rowList_LeadCurrentGID_.clear();
  vnodeGIDVector_.clear();
  vsrcGIDVector_.clear();
  nonlinGIDVector_.clear();

  CktNodeList::const_iterator it = topology_.getOrderedNodeList().begin();
  CktNodeList::const_iterator end = topology_.getOrderedNodeList().end();
  for ( ; it != end; ++it )
  {
    // Own this node, get solution, state, and store information.
    if( (*it)->get_IsOwned() && ( (*it)->get_gID() != -1 ) )
    {
      // Get solution variable GIDs.
      rowList_GID_.insert(rowList_GID_.end(),
          (*it)->get_SolnVarGIDList().begin(),
          (*it)->get_SolnVarGIDList().end());


      // Collect device nodes information.
      if ( (*it)->type() == _DNODE )
      {
        CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it);

        // Get state variable GIDs.
        if (cktNodeDevPtr->stateVarCount())
        {
          rowList_StateGID_.insert(rowList_StateGID_.end(),
              cktNodeDevPtr->get_StateVarGIDList().begin(),
              cktNodeDevPtr->get_StateVarGIDList().end());
        }

        // Get store variable GIDs.
        if (cktNodeDevPtr->storeVarCount())
        {
          rowList_StoreGID_.insert(rowList_StoreGID_.end(),
              cktNodeDevPtr->get_StoreVarGIDList().begin(),
              cktNodeDevPtr->get_StoreVarGIDList().end() );
        }

        // Get lead current variable GIDs.
        if (cktNodeDevPtr->branchDataVarCount())
        {
          rowList_LeadCurrentGID_.insert( rowList_LeadCurrentGID_.end(),
              cktNodeDevPtr->get_LeadCurrentVarGIDList().begin(),
              cktNodeDevPtr->get_LeadCurrentVarGIDList().end() );
        }

        const std::string & id = (*it)->get_id();
        std::string::size_type col = id.find_first_of(':');

        if ( id[col+1] == 'V' || id[col+1] == 'v' )
        {
          vsrcGIDVector_.insert( vsrcGIDVector_.end(),
              cktNodeDevPtr->get_SolnVarGIDList().begin(),
              cktNodeDevPtr->get_SolnVarGIDList().end() );

          vsrcGIDVector_.insert( vsrcGIDVector_.end(),
              cktNodeDevPtr->get_ExtSolnVarGIDList().begin(),
              cktNodeDevPtr->get_ExtSolnVarGIDList().end() );
        }

        if (!(cktNodeDevPtr->deviceInstance()->isLinearDevice()))
        {
          nonlinGIDVector_.insert( nonlinGIDVector_.end(),
              cktNodeDevPtr->get_SolnVarGIDList().begin(),
              cktNodeDevPtr->get_SolnVarGIDList().end() );

          nonlinGIDVector_.insert( nonlinGIDVector_.end(),
              cktNodeDevPtr->get_ExtSolnVarGIDList().begin(),
              cktNodeDevPtr->get_ExtSolnVarGIDList().end() );
        }
      }
   
      // Collect voltage nodes.
      if ( (*it)->type() == _VNODE )
      {
        vnodeGIDVector_.push_back( *((*it)->get_SolnVarGIDList().begin()) );
      }
    }

    // Do not own this node. Still get solution information.
    if( !( (*it)->get_IsOwned() ) && ( (*it)->get_gID() != -1 ) )
    {
      // Get solution variable GIDs.
      std::vector<int>::const_iterator it_svL = (*it)->get_SolnVarGIDList().begin();
      std::vector<int>::const_iterator it_svL_end = (*it)->get_SolnVarGIDList().end();
      for (; it_svL != it_svL_end; ++it_svL )
      {
        rowList_ExternGID_.push_back( std::pair<int,int>( *it_svL, (*it)->get_ProcNum() ) );
      }
    }
  } 

  // Make sure the nonlinGIDVector_ has unique entries.
  // First sort, then erase duplicates, then remove ground node GID.
  std::sort( nonlinGIDVector_.begin(), nonlinGIDVector_.end() );
  nonlinGIDVector_.erase(std::unique(nonlinGIDVector_.begin(), nonlinGIDVector_.end() ), nonlinGIDVector_.end() );
  if ( nonlinGIDVector_.size() > 0 )
  {
    if ( nonlinGIDVector_.front() == -1 )
      nonlinGIDVector_.erase( nonlinGIDVector_.begin() );
  }

}

} // namespace Topo
} // namespace Xyce

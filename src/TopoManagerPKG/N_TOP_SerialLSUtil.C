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
#include <N_PDS_Manager.h>
#include <N_PDS_Migrate.h>
#include <N_PDS_Node.h>
#include <N_PDS_ParMap.h>
#include <N_TOP_CktGraph.h>
#include <N_TOP_CktNode.h>
#include <N_TOP_CktNode_Dev.h>
#include <N_TOP_SerialLSUtil.h>
#include <N_TOP_Topology.h>
#include <N_IO_HangingResistor.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Misc.h>

#include <Teuchos_Utils.hpp>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : SerialLSUtil::SerialLSUtil
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter  SNL, Parallel Computational Sciences
// Creation Date : 5/25/05
//-----------------------------------------------------------------------------
SerialLSUtil::SerialLSUtil(
  Topology &            topology,
  Parallel::Manager &   pds_manager)
  : Linear::QueryUtil(),
    topology_(topology),
    pdsManager_(pds_manager),
    numGlobalNodes_(0),
    numGlobalRows_(0),
    numGlobalStateVars_(0),
    numGlobalStoreVars_(0),
    numGlobalLeadCurrentVars_(0),
    numGlobalNZs_(0)
{
}

//-----------------------------------------------------------------------------
// Function      : SerialLSUtil::cleanRowLists
// Purpose       : reduce memory consumption of query utility
// Special Notes : the builder should call this after all maps and graphs
//                 have been made
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 1/25/19
//-----------------------------------------------------------------------------
void SerialLSUtil::cleanRowLists()
{
  // rowList_GID_.clear(); keep this one for separated object analysis
  rowList_StateGID_.clear();
  rowList_StoreGID_.clear();
  rowList_LeadCurrentGID_.clear();
  rowList_NumNZs_.clear();
  for( int i = 0; i < numGlobalRows_; ++i )
    rowList_ColList_[i].clear();

  isClean_ = true; // don't print out arrays, they have been cleared.
}

//-----------------------------------------------------------------------------
// Function      : SerialLSUtil::operator<<
// Purpose       : generate utility with reference to topology
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/15/00
//-----------------------------------------------------------------------------
std::ostream & operator<< (std::ostream & os, const SerialLSUtil & tlsu )
{
  os << "Serial Topology LS Utility" << std::endl;
  os << "-------------------" << std::endl;
  os << "Num Global Nodes: " << tlsu.numGlobalNodes_ << std::endl;
  os << "Num Global Rows: " << tlsu.numGlobalRows_ << std::endl;

  os << "Num Global State Vars: " << tlsu.numGlobalStateVars_ << std::endl;

  os << "Num Global Store Vars: " << tlsu.numGlobalStoreVars_ << std::endl;

  os << "Num Global LeadCurrent Vars: " << tlsu.numGlobalLeadCurrentVars_ << std::endl;

  os << "Num Global NZs: " << tlsu.numGlobalNZs_ << std::endl;

  os << "GID Array: ";
  for( int i = 0; i < tlsu.numGlobalRows_; ++i )
    os << tlsu.rowList_GID_[i] << " ";
  os << std::endl;

  if (!tlsu.isClean_)
  {
    os << "State GID Array: ";
    for( int i = 0; i < tlsu.numGlobalStateVars_; ++i )
      os << tlsu.rowList_StateGID_[i] << " ";
    os << std::endl;

    os << "Store GID Array: ";
    for( int i = 0; i < tlsu.numGlobalStoreVars_; ++i )
      os << tlsu.rowList_StoreGID_[i] << " ";
    os << std::endl;

    os << "LeadCurrent GID Array: ";
    for( int i = 0; i < tlsu.numGlobalLeadCurrentVars_; ++i )
      os << tlsu.rowList_LeadCurrentGID_[i] << " ";
    os << std::endl;

    os << "NZ Array: ";
    for( int i = 0; i < tlsu.numGlobalRows_; ++i )
      os << tlsu.rowList_NumNZs_[i] << " ";
    os << std::endl;

    os << "Col Index Array: " << std::endl;
    for( int i = 0; i < tlsu.numGlobalRows_; ++i )
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
// Function      : SerialLSUtil::setupRowCol
// Purpose       : Setup row/col data for linear solver including reorder
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/1/01
//-----------------------------------------------------------------------------
bool SerialLSUtil::setupRowCol()
{

  setupNodeGIDs();

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
// Function      : SerialLSUtil::testVoltageNodeConnectivity
// Purpose       : testing of voltage node connectivity for problems
// Special Notes : initially, just warn if a voltage node has only 1 connection
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/31/03
//-----------------------------------------------------------------------------
bool SerialLSUtil::testVoltageNodeConnectivity_()
{
  CktNodeList::const_iterator it_cnL = topology_.getOrderedNodeList().begin();
  CktNodeList::const_iterator end_cnL = topology_.getOrderedNodeList().end();

  int i, j, k;
  int gid, num_gid;
  int num_nodes = 0;
  std::vector<int> gidList;
  std::vector<int> a,b;
  unordered_map<int, int> gid_pos;

  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    if ((*it_cnL)->type() == _VNODE)
    {
      gid = (*it_cnL)->get_gID();
      if (gid >= 0)
      {
        num_gid = topology_.numAdjNodes( gid );
        if (num_gid <= 1)
        {
          connToOneTermIDs_.insert( (*it_cnL)->get_id() );
        }
        gid_pos[gid] = num_nodes++; 
      }
    }
  }

  // Go through all the devices and color the voltage nodes by their lead groupings.
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
  // All that will remain is ungrounded colors or singles.
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

  it_cnL = topology_.getOrderedNodeList().begin();

  // Collect any nodes that do not have DC path to ground; indicated by any node_color NOT set to -1.
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

  // Output the warnings through the topology object.
  topology_.outputTopoWarnings(connToOneTermIDs_, noDCPathIDs_);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : SerialLSUtil::setupNodeGIDs
// Purpose       : Generate Ordering and Var GIDs.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/22/01
//-----------------------------------------------------------------------------
bool SerialLSUtil::setupNodeGIDs()
{
  //get lists of owned and boundary/ghost node GIDs
  topology_.generateOrderedNodeList();

  //loop over nodes and reset all GIDs to new lex. ordering
  int newGID = 0; 
  for (CktNodeList::const_iterator it_cnL = topology_.getOrderedNodeList().begin(), end_cnL = topology_.getOrderedNodeList().end(); it_cnL != end_cnL; ++it_cnL )
  {
    if( (*it_cnL)->get_id() != "0" )
    {
      (*it_cnL)->set_gID( newGID );
    
      // Increment count.
      newGID++;
    }
    else
    {
      (*it_cnL)->set_gID( -1 );
    }
  }

  numGlobalNodes_ = newGID;

  return true;
}


//-----------------------------------------------------------------------------
// Function      : SerialLSUtil::setupSolnAndStateGIDs
// Purpose       : Generate Ordering and Var GIDs.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/2/01
//-----------------------------------------------------------------------------
bool SerialLSUtil::setupSolnAndStateGIDs()
{
  //get soln and state var counts from all owned v-nodes and d-nodes
  std::vector<int> rowCountVec(numGlobalNodes_+1);
  std::vector<int> stateCountVec(numGlobalNodes_+1);
  std::vector<int> storeCountVec(numGlobalNodes_+1);
  std::vector<int> leadCurrentCountVec(numGlobalNodes_+1);

  int Loc = 0;
  rowCountVec[Loc] = 0;
  stateCountVec[Loc] = 0;
  storeCountVec[Loc] = 0;
  leadCurrentCountVec[Loc] = 0;

  for (CktNodeList::const_iterator it_cnL = topology_.getOrderedNodeList().begin(), end_cnL = topology_.getOrderedNodeList().end(); it_cnL != end_cnL; ++it_cnL )
  {
    if( (*it_cnL)->get_IsOwned() && (*it_cnL)->get_gID() != -1 )
    {
      rowCountVec[Loc+1] = rowCountVec[Loc] + (*it_cnL)->solnVarCount();
      if( (*it_cnL)->type() == _DNODE )
      {
        CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_cnL);
        stateCountVec[Loc+1] = stateCountVec[Loc] + cktNodeDevPtr->stateVarCount();
        storeCountVec[Loc+1] = storeCountVec[Loc] + cktNodeDevPtr->storeVarCount();
        leadCurrentCountVec[Loc+1] = leadCurrentCountVec[Loc] + cktNodeDevPtr->branchDataVarCount();
      }
      else
      {
        stateCountVec[Loc+1] = stateCountVec[Loc];
        storeCountVec[Loc+1] = storeCountVec[Loc];
        leadCurrentCountVec[Loc+1] = leadCurrentCountVec[Loc];
      }
      ++Loc;
    }
  }

  numGlobalRows_ = rowCountVec[Loc];
  numGlobalStateVars_ = stateCountVec[Loc];
  numGlobalStoreVars_ = storeCountVec[Loc];
  numGlobalLeadCurrentVars_ = leadCurrentCountVec[Loc];

  //loop over nodes and assign soln and state GIDs
  Loc = 0;
  for (CktNodeList::iterator it_cnL = topology_.getOrderedNodeList().begin(), end_cnL = topology_.getOrderedNodeList().end(); it_cnL != end_cnL; ++it_cnL )
  {
    if( (*it_cnL)->get_IsOwned() )
    {
      if ( (*it_cnL)->get_gID() != -1 )
      {
        // Set solution GIDs
        int count = rowCountVec[Loc+1] - rowCountVec[Loc];
        if ( count )
        {
          std::vector<int>& rowGIDs = (*it_cnL)->get_SolnVarGIDList();
          rowGIDs.resize( count );
          iota( rowGIDs.begin(), rowGIDs.end(), rowCountVec[Loc] );
        }

        if( (*it_cnL)->type() == _DNODE )
        {
          CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it_cnL);

          // Set state GIDs
          count = stateCountVec[Loc+1] - stateCountVec[Loc];
          if ( count )
          {
            std::vector<int>& stateGIDs = cktNodeDevPtr->get_StateVarGIDList();
            stateGIDs.resize( count );
            iota( stateGIDs.begin(), stateGIDs.end(), stateCountVec[Loc] );
          }
          // Set store GIDs
          count = storeCountVec[Loc+1] - storeCountVec[Loc];
          if ( count )
          {
            std::vector<int>& storeGIDs = cktNodeDevPtr->get_StoreVarGIDList();
            storeGIDs.resize( count );
            iota( storeGIDs.begin(), storeGIDs.end(), storeCountVec[Loc] );
          }
          // Set lead current GIDs
          count = leadCurrentCountVec[Loc+1] - leadCurrentCountVec[Loc];
          if ( count )
          {
            std::vector<int>& leadCurrentGIDs = cktNodeDevPtr->get_LeadCurrentVarGIDList();
            leadCurrentGIDs.resize( count );
            iota( leadCurrentGIDs.begin(), leadCurrentGIDs.end(), leadCurrentCountVec[Loc] );
          }
        }
        Loc++;
      }
      else
      {
        (*it_cnL)->set_SolnVarGIDList( std::vector<int>(1,-1) );
      }
    }
    else
    {
      Report::DevelFatal0() << "Node: " << (*it_cnL)->get_id() << ", global index ("
                            << Teuchos::Utils::toString( (*it_cnL)->get_gID() )
                            << ") is NOT found!";
    }
  }

  topology_.registerGIDs();

  return true;

}

//-----------------------------------------------------------------------------
// Function      : SerialLSUtil::generateRowColData
// Purpose       : Generate row/col data for linear solver.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/26/00
//-----------------------------------------------------------------------------
void SerialLSUtil::generateRowColData()
{
  //--- set numGlobalRows_, numGlobalStateVars_, and resize lists
  numGlobalRows_ = rowList_GID_.size();
  numGlobalStateVars_ = rowList_StateGID_.size();
  numGlobalStoreVars_ = rowList_StoreGID_.size();
  numGlobalLeadCurrentVars_ = rowList_LeadCurrentGID_.size();

  if (DEBUG_TOPOLOGY)
    Xyce::dout() << "Serial Topology Util Vals: " << std::endl
                 << numGlobalRows_ << " " << numGlobalRows_ << std::endl
                 << numGlobalStateVars_ << " " << numGlobalStateVars_ << std::endl
                 << numGlobalStoreVars_ << " " << numGlobalStoreVars_ << std::endl
                 << numGlobalLeadCurrentVars_ << " " << numGlobalLeadCurrentVars_ << std::endl;

  int numRows = numGlobalRows_ + 1;

  //Build global to local map for speed
  unordered_map<int,int> GtoL_Map;
  for( int i = 0; i < numGlobalRows_; ++i )
    GtoL_Map[ rowList_GID_[i] ] = i;
  GtoL_Map[ -1 ] = numGlobalRows_;

  rowList_ColList_.resize( numRows );
  rowList_NumNZs_.resize( numRows );

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
          {
            rowList_ColList_[ GtoL_Map[ gids[i] ] ].push_back( gids[ stamp[i][j] ] );
          }
        }
      }
    }
  }


  if (DEBUG_TOPOLOGY)
    Xyce::dout() << Xyce::section_divider << std::endl
                 << "Row List of NZ Cols extracted" << std::endl
                 << Xyce::section_divider << std::endl;

  numGlobalNZs_ = 0;

#ifdef Xyce_NOX_LOCA_ARTIFICIAL_HOMOTOPY_SUPPORT
    //add in diagonal for homotopy support
    for( int i = 0; i < numGlobalRows_; ++i )
      rowList_ColList_[i].push_back( rowList_GID_[i] );
#endif

  //--- sort list of col indices and get rid of redundancies
  //    generate num of NZs data
  for( int i = 0; i < numRows; ++i )
  {
    std::sort(rowList_ColList_[i].begin(), rowList_ColList_[i].end());
    rowList_ColList_[i].erase(std::unique(rowList_ColList_[i].begin(), rowList_ColList_[i].end()), rowList_ColList_[i].end());

    rowList_NumNZs_[i] = rowList_ColList_[i].size();

    numGlobalNZs_ += rowList_NumNZs_[i];
  }

  if (DEBUG_TOPOLOGY)
    Xyce::dout() << Xyce::section_divider << std::endl
                 << "Row List of NZ Cols cleaned up" << std::endl
                 << Xyce::section_divider << std::endl
                 << *this;

}

//-----------------------------------------------------------------------------
// Function      : SerialLSUtil::extractAllGIDsFromTopology
// Purpose       : Fill all GID vectors using ordered node list.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/26/00
//-----------------------------------------------------------------------------
void SerialLSUtil::extractAllGIDsFromTopology()
{
  // Clear all vectors.
  rowList_GID_.clear();
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
      // Get solution variable GIDs [for both voltage and device nodes]
      rowList_GID_.insert(rowList_GID_.end(),
          (*it)->get_SolnVarGIDList().begin(),
          (*it)->get_SolnVarGIDList().end());

      if ( (*it)->type() == _DNODE )
      {
        CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*it);
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
      }

      // Collect voltage nodes.
      if ( (*it)->type() == _VNODE )
      {
        vnodeGIDVector_.push_back( *((*it)->get_SolnVarGIDList().begin()) );
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

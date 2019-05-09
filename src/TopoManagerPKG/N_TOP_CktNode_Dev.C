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

#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceMgr.h>
#include <N_ERH_ErrorMgr.h>
#include <N_TOP_CktGraph.h>
#include <N_TOP_CktNode_Dev.h>
#include <N_TOP_Topology.h>
#include <N_UTL_Expression.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Interface_Enum_Types.h>
#include <N_UTL_Misc.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::~CktNode_Dev
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/22/03
//-----------------------------------------------------------------------------
CktNode_Dev::~CktNode_Dev()
{ 
  if (instanceBlock_) delete instanceBlock_; 
}

void
CktNode_Dev::loadNodeSymbols(
  Topology &            topology) const
{
  Util::SymbolTable &node_symbols = topology.getNodeSymbols();
  CktGraph &main_graph = topology.getMainGraph();

  deviceInstance_->loadNodeSymbols(node_symbols);

  const std::string &id = get_id();
  std::string::size_type col = id.find_first_of(':');

  if (id[col + 1] == 'V' || id[col + 1] == 'v' )
  {
    std::vector<NodeID> adj_ids;
    main_graph.returnAdjIDs(get_nodeID(), adj_ids);
    for (std::vector<NodeID>::const_iterator it_nodeid = adj_ids.begin(), end_nodeid = adj_ids.end(); it_nodeid != end_nodeid; ++it_nodeid)
      if ((*it_nodeid).first != "0" )
        node_symbols[Util::VSRC_SYMBOL][(*it_nodeid).first] = 0;
  }
}

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::instantiate
// Purpose       : instantiate device with dev pkg
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/22/03
//-----------------------------------------------------------------------------
bool CktNode_Dev::instantiate()
{
  if (deviceInstance_)
    return false;
  else
  {
    deviceInstance_ = deviceManager_->addDeviceInstance(*instanceBlock_);

    // for certain invalid device names, deviceInstance_ is set to the 
    // null pointer by the previous call.  Exit immediately if that happens, 
    // since Xyce will likely core dump later
    if (deviceInstance_ == 0)
    {
      Report::DevelFatal().in("CktNode_Dev::instantiate") 
         << "null Device Instance pointer";
    }

    // Delete the instanceBlock, it is not used after instantiation
    delete instanceBlock_; instanceBlock_ = 0; 
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::getDevState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
Device::DeviceState *
CktNode_Dev::getDevState()
{
  return deviceInstance_->getInternalState();
}

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::setDevState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
bool
CktNode_Dev::setDevState(
  const Device::DeviceState &   state)
{
  return deviceInstance_->setInternalState( state );
}

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::solnVarCount
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/26/01
//-----------------------------------------------------------------------------
int CktNode_Dev::solnVarCount() const
{
  return deviceInstance_->getNumIntVars();
}

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::stateVarCount
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/26/01
//-----------------------------------------------------------------------------
int CktNode_Dev::stateVarCount() const
{
  return deviceInstance_->getNumStateVars();
}

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::storeVarCount
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
int CktNode_Dev::storeVarCount() const
{
  return deviceInstance_->getNumStoreVars();
}

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::branchDataVarCount
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
int CktNode_Dev::branchDataVarCount() const
{
  return deviceInstance_->getNumBranchDataVars();
}
//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::leadConnect
// Purpose       : Find which leads have connections to each other
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 6/20/05
//-----------------------------------------------------------------------------
const std::vector<int> & CktNode_Dev::leadConnect() const
{
  return deviceInstance_->getDevConMap();
}

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::registerLIDswithDev
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/18/02
//-----------------------------------------------------------------------------
void
CktNode_Dev::registerLIDswithDev(
  const std::vector<int> &      intLIDVec,
  const std::vector<int> &      extLIDVec )
{
  deviceInstance_->registerLIDs( intLIDVec, extLIDVec );
}

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::registerStateLIDswithDev
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/18/02
//-----------------------------------------------------------------------------
void
CktNode_Dev::registerStateLIDswithDev(
  const std::vector<int> &      stateLIDVec )
{
  deviceInstance_->registerStateLIDs( stateLIDVec );
}

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::registerStoreLIDswithDev
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void
CktNode_Dev::registerStoreLIDswithDev(
  const std::vector<int> &      storeLIDVec)
{
  deviceInstance_->registerStoreLIDs( storeLIDVec );
}

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::registerLeadCurrentLIDswithDev
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void
CktNode_Dev::registerLeadCurrentLIDswithDev(
  const std::vector<int> &      leadCurrentLIDVec)
{
  deviceInstance_->registerBranchDataLIDs( leadCurrentLIDVec );
}


//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::registerDepLIDswithDev
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void
CktNode_Dev::registerDepLIDswithDev(
  const std::vector< std::vector<int> > &       depLIDVec )
{
  deviceInstance_->registerDepSolnLIDs( depLIDVec );
}

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::getDepSolnVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/05/01
//-----------------------------------------------------------------------------
void
CktNode_Dev::getDepSolnVars(
  std::vector<NodeID> &         dsVars )
{
  // Get dependent solution variables from the device and the Depend structure to determine which node type to look for in the graph.
  dsVars.clear();

  const std::vector<std::string> &solution_variables = deviceInstance_->getDepSolnVars();
  for (std::vector<std::string>::const_iterator it_sV = solution_variables.begin(), it_sV_end = solution_variables.end(); it_sV != it_sV_end; ++it_sV )
  {
    bool found = false;
    const std::vector<Device::Depend> &depParams = deviceInstance_->getDependentParams();
    for (std::vector<Device::Depend>::const_iterator it_sdV = depParams.begin(), it_sdV_end = depParams.end(); it_sdV != it_sdV_end; ++it_sdV )
    {
      int type = (it_sdV->expr)->get_type( *it_sV );
      if (type == XEXP_NODE)
      {
        dsVars.push_back( NodeID( *it_sV, _VNODE ) );
        found = true;
        break;
      }
      if (type == XEXP_INSTANCE || type == XEXP_LEAD)
      {
        dsVars.push_back( NodeID( *it_sV, _DNODE ) );
        found = true;
        break;
      }
    }
    // Punt if the node has not been found, hope that the graph only
    // has one node associated with this name.
    if (!found)
    {
      dsVars.push_back(NodeID( *it_sV, -1 ));
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::registerDepSolnGIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
void CktNode_Dev::registerDepSolnGIDs(
  const std::vector< std::vector<int> > & dsGIDs )
{
  if (DEBUG_TOPOLOGY)
  {
    Xyce::dout() << "reg secondary gids: " << get_id() << std::endl;

    for( unsigned int i = 0; i < dsGIDs.size(); ++i )
    {
      for( unsigned int j = 0; j < dsGIDs[i].size(); ++j )
        Xyce::dout() << " " << dsGIDs[i][j];
      Xyce::dout() << std::endl;
    }

    Xyce::dout() << std::endl;
  }

  deviceInstance_->registerDepSolnGIDs( dsGIDs );
}

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::put
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
std::ostream &
CktNode_Dev::put(std::ostream& os) const
{
  os << "CN_Dev: " << get_id() << std::endl;
  os << "   GID: " << gID_ << "  Proc: " << procNum_ << std::endl;
  os << "   Owned: " << isOwned_ << std::endl;
  os << "   Soln Var GID List: ";
  int count=0;
  std::vector<int>::const_iterator it_iL = solnVarGIDList_.begin();
  std::vector<int>::const_iterator it_iL_end = solnVarGIDList_.end();
  for( ; it_iL != it_iL_end; ++it_iL )
  {
    os << *it_iL << "  ";
    if (count >= 12) {os << std::endl;count=0;}
    else ++count;
  }
  os << std::endl;

  os << "   Ext Soln Var GID List: ";
  count=0;
  it_iL = extSolnVarGIDList_.begin();
  it_iL_end = extSolnVarGIDList_.end();
  for( ; it_iL != it_iL_end; ++it_iL )
  {
    os << *it_iL << "  ";
    if (count >= 12) {os << std::endl;count=0;}
    else ++count;
  }
  os << std::endl;

  os << "   State Var GID List: ";
  it_iL = stateVarGIDList_.begin();
  it_iL_end = stateVarGIDList_.end();
  for( ; it_iL != it_iL_end; ++it_iL )
  {
    os << *it_iL << "  ";
  }
  os << std::endl;

  os << "   Store Var GID List: ";
  it_iL = storeVarGIDList_.begin();
  it_iL_end = storeVarGIDList_.end();
  for( ; it_iL != it_iL_end; ++it_iL )
  {
    os << *it_iL << "  ";
  }
  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/20/02
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & CktNode_Dev::jacobianStamp() const
{
  return deviceInstance_->jacobianStamp();
}

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::registerJacLIDswithDev
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/23/02
//-----------------------------------------------------------------------------
void CktNode_Dev::registerJacLIDswithDev( const std::vector< std::vector<int> > & jacLIDVec )
{
  deviceInstance_->registerJacLIDs( jacLIDVec );
}

//-----------------------------------------------------------------------------
// Function      : CktNode_Dev::varTypeList
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/23/02
//-----------------------------------------------------------------------------
void CktNode_Dev::varTypeList( std::vector<char> & variableTypes_ ) const
{
  std::vector<char> typeList;
  deviceInstance_->varTypes(typeList);
  if (typeList.empty())
  {
    variableTypes_.insert(variableTypes_.end(), deviceInstance_->getNumIntVars(), 'V');
  }
  else
  {
    variableTypes_.insert(variableTypes_.end(), typeList.begin(), typeList.end());
  }
}

} // namespace Topo
} // namespace Xyce

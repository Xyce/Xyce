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

#ifndef N_TOP_CktNode_Dev_h
#define N_TOP_CktNode_Dev_h 1

#include <N_DEV_fwd.h>

#include <N_TOP_CktNode.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceBlock.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : CktNode_Dev
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class CktNode_Dev : public CktNode
{
public:
  CktNode_Dev(
    const Device::InstanceBlock &instance_block,
    Device::DeviceMgr &          device_manager)
    : CktNode(instance_block.getInstanceName().getEncodedName()),
      deviceInstance_(0),
      deviceManager_(&device_manager),
      instanceBlock_(0)
  {
    Xyce::set_node_type(nodeID_,_DNODE);
    instanceBlock_ = new Device::InstanceBlock(instance_block);
  }

  // Destructor
  virtual ~CktNode_Dev();

private:
  CktNode_Dev(const CktNode_Dev &);
  CktNode_Dev & operator=(const CktNode_Dev &);

public:

  virtual void loadNodeSymbols(Topology &topology) const;

  bool instantiated() const
  {
    return deviceInstance_!=0;
  }
  
  bool instantiate();

  const Device::InstanceBlock *deviceInstanceBlock()
  {
    return instanceBlock_;
  }

  const Device::DeviceInstance *deviceInstance()
  {
    return deviceInstance_;
  }

  // Get's the device state object.
  Device::DeviceState * getDevState();

  // Set's the device state object.
  bool setDevState(const Device::DeviceState & state);

  void registerLIDswithDev( const std::vector<int> & intLIDVec,
                            const std::vector<int> & extLIDVec );
  void registerStateLIDswithDev( const std::vector<int> & stateLIDVec );
  void registerStoreLIDswithDev( const std::vector<int> & storeLIDVec );
  void registerLeadCurrentLIDswithDev( const std::vector<int> & leadCurrentLIDVec);

  void registerDepLIDswithDev( const std::vector< std::vector<int> > & depLIDVec );

  // Setup secondary dependencies.
  void getDepSolnVars( std::vector< NodeID >& dsVars );
  void registerDepSolnGIDs(const std::vector< std::vector< int > > & dsGIDs);

  const std::vector<int> & get_ExtSolnVarGIDList() const {
    return extSolnVarGIDList_;
  }
 
  std::vector<int> & get_ExtSolnVarGIDList() {
    return extSolnVarGIDList_;
  }

  void set_ExtSolnVarGIDList(const std::vector<int> & svGIDList) {
    extSolnVarGIDList_ = svGIDList;
  }

  const std::vector<int> & get_StateVarGIDList() const {
    return stateVarGIDList_;
  }

  std::vector<int> & get_StateVarGIDList() {
    return stateVarGIDList_;
  }

  void set_StateVarGIDList(const std::vector<int> & svGIDList) {
    stateVarGIDList_ = svGIDList;
  }

  const std::vector<int> & get_StoreVarGIDList() const {
    return storeVarGIDList_;
  }

  std::vector<int> & get_StoreVarGIDList() {
    return storeVarGIDList_;
  }

  void set_StoreVarGIDList(const std::vector<int> & svGIDList) {
    storeVarGIDList_ = svGIDList;
  }

  const std::vector<int> & get_LeadCurrentVarGIDList() const {
    return leadCurrentVarGIDList_;
  }

  std::vector<int> & get_LeadCurrentVarGIDList() {
    return leadCurrentVarGIDList_;
  }

  void set_LeadCurrentVarGIDList(const std::vector<int> & lcvGIDList) {
    leadCurrentVarGIDList_ = lcvGIDList;
  }

  void get_DepSolnGIDVec(std::vector< std::vector<int> > & dsGIDs) const {
    dsGIDs = depSolnGIDVec_;
  }
  void set_DepSolnGIDVec(const std::vector< std::vector<int> > & dsGIDs) {
    depSolnGIDVec_ = dsGIDs;
  }

  int solnVarCount() const;
  int stateVarCount() const;
  int storeVarCount() const;
  int branchDataVarCount() const;

  int extSolnVarCount() const {
    return extSolnVarGIDList_.size();
  }

  int depSolnVarCount() const {
    return depSolnGIDVec_.size();
  }
 
  const std::vector<int> & leadConnect() const;

  const std::vector< std::vector<int> > & jacobianStamp() const;
  const std::vector<int> & get_DepSolnGIDJacVec() { return deviceInstance_->getDepSolnGIDVec(); }

  void registerJacLIDswithDev( const std::vector< std::vector<int> > & jacLIDVec );

  void varTypeList( std::vector<char> & varTypeVec ) const;

private:

  // Pointer to a device instance.
  Device::DeviceInstance *      deviceInstance_;
  Device::DeviceMgr *           deviceManager_;
  Device::InstanceBlock *       instanceBlock_;

  std::vector<int>              extSolnVarGIDList_;
  std::vector<int>              stateVarGIDList_;               /// State variable global id list.
  std::vector<int>              storeVarGIDList_;
  std::vector<int>              leadCurrentVarGIDList_;

  std::vector< std::vector<int> > depSolnGIDVec_;

public:

    std::ostream & put(std::ostream & os) const;

};

} // namespace Topo
} // namespace Xyce

#endif

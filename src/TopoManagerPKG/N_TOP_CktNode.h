//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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

#ifndef N_TOP_CktNode_h
#define N_TOP_CktNode_h

#include <iosfwd>
#include <string>
#include <list>
#include <vector>
#include <map>

#include <N_TOP_fwd.h>
#include <N_UTL_Misc.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : CktNode
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class CktNode
{
  friend std::ostream & operator << (std::ostream & os, const CktNode & cn);

public:
  CktNode(
    const std::string & nodeID,
    const int & globalID = 0)
    : gID_(globalID),
      procNum_(0),
      isOwned_(true),
      solnVarGIDList_()
  {
    Xyce::set_node_id(nodeID_,nodeID); 
  }

  // Destructor
  virtual ~CktNode() { }

private:
  CktNode(const CktNode &);
  CktNode &operator=(const CktNode &);

public:
  bool operator == (const CktNode & right) const { return nodeID_ == right.nodeID_; }
  bool operator != (const CktNode & right) const { return nodeID_ != right.nodeID_; }

  //-------- get and set methods for attributes

  // Get circuit NodeID.
  const NodeID& get_nodeID() const { return nodeID_; }

  // Get circuit node id string.
  const std::string & get_id() const { return Xyce::get_node_id(nodeID_); }

  // Get circuit node global id.
  const int & get_gID() const { return gID_; }

  // Set circuit node global id.
  void set_gID(const int & globalID) { gID_ = globalID; }

  int type() const { return Xyce::get_node_type(nodeID_); }

  const std::vector<int> & get_SolnVarGIDList() const {
    return solnVarGIDList_;
  }

  std::vector<int> & get_SolnVarGIDList() {
    return solnVarGIDList_;
  }

  void set_SolnVarGIDList(const std::vector<int> & svGIDVec) {
    solnVarGIDList_ = svGIDVec;
  }

  virtual int solnVarCount() const {
    return solnVarGIDList_.size();
  }

  // Get the processor number.
  const int & get_ProcNum() const { return procNum_; }

  // Set the processor number.
  void set_ProcNum(const int & pNum) { procNum_ = pNum; }

  // Get the processor ownership flag.
  const bool & get_IsOwned() const { return isOwned_; }

  // Set the processor ownership flag.
  void set_IsOwned(const bool & owned) { isOwned_ = owned; }

  //------- Added for use with the outputFileName function
  virtual void loadNodeSymbols(Topology &topology) const = 0;

  virtual void varTypeList( std::vector<char> & varTypeVec ) const = 0;

  virtual std::ostream & put(std::ostream & os) const = 0;

protected:
  NodeID                  nodeID_;        /// Node id.
  int                     gID_;           /// Global node id.
  int                     procNum_;       /// Processor number.
  bool                    isOwned_;       /// Processor ownership flag.

  std::vector<int>        solnVarGIDList_;/// Solution variable global id list.
};

} // namespace Topo
} // namespace Xyce

#endif

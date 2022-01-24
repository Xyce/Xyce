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

#ifndef N_TOP_CktNode_V_h
#define N_TOP_CktNode_V_h 1

#include <N_TOP_CktNode.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : CktNode_V
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class CktNode_V : public CktNode
{
public:
  // Constructor
  CktNode_V(const std::string & nodeID, const int & globalID = 0)
    : CktNode(nodeID, globalID)
  {
    Xyce::set_node_type(nodeID_,_VNODE);
  }

private:
  // Copy constructor
  CktNode_V(const CktNode_V & right);
  CktNode_V &operator=(const CktNode_V & right);

public:
  // Destructor
  virtual ~CktNode_V() {}

  int solnVarCount() const { return 1; }

  virtual void loadNodeSymbols(Topology &topology) const;

  virtual void varTypeList( std::vector<char> & varTypeVec ) const;
  
  std::ostream & put(std::ostream & os) const;

};

} // namespace Topo
} // namespace Xyce

#endif

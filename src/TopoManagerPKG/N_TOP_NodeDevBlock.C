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
// Creation Date  : 3/3/01
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_TOP_NodeDevBlock.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : NodeDevBlock::clear
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems modeling
// Creation Date : 2/8/2010
//-----------------------------------------------------------------------------
void NodeDevBlock::clear()
{
  devBlock_.clear();
  nodeList_.clear();
}


//-----------------------------------------------------------------------------
// Function      : NodeDevBlock::addNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 03/16/06
//-----------------------------------------------------------------------------
void NodeDevBlock::addNode(const std::string & p)
{
  nodeList_.push_back(p);
}

//-----------------------------------------------------------------------------
// Function      : NodeDevBlock::get_NodeList
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
const std::vector<std::string> & NodeDevBlock::get_NodeList() const
{
  return nodeList_;
}

//-----------------------------------------------------------------------------
// Function      : NodeDevBlock::get_NodeList
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
std::vector<std::string> & NodeDevBlock::get_NodeList() 
{
  return nodeList_;
}

//-----------------------------------------------------------------------------
//// Function      : NodeDevBlock::set_NodeList
//// Purpose       :
//// Special Notes :
//// Scope         : public
//// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
//// Creation Date : 5/18/00
////-----------------------------------------------------------------------------
void NodeDevBlock::set_NodeList(const std::vector<std::string> & nList)
{
  nodeList_.assign(nList.begin(), nList.end());
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/3/01
//-----------------------------------------------------------------------------
std::ostream & operator<< ( std::ostream & os, const NodeDevBlock & ndb )
{
  os << "NodeDevBlock: " << ndb.devBlock_.getInstanceName().getEncodedName() << std::endl;
  os << " Connected Nodes: ";
  std::vector<std::string>::const_iterator it_tpL = ndb.get_NodeList().begin();
  std::vector<std::string>::const_iterator it_tpL_end = ndb.get_NodeList().end();
  for( ; it_tpL != it_tpL_end; ++it_tpL)
    os << "   " << *it_tpL;
  os << std::endl;

  if( ndb.devBlock_.getInstanceName() != "" ) os << ndb.devBlock_ << std::endl;
  return os << std::endl;
}

} // namespace Topo

} // namespace Xyce

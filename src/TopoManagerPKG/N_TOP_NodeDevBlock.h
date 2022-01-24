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
// Purpose        : Block of Node and Dev data for passing between pkgs
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
//-----------------------------------------------------------------------------

#ifndef N_TOP_NodeDevBlock_h
#define N_TOP_NodeDevBlock_h 1

#include <iosfwd>

#include <N_DEV_DeviceBlock.h>
#include <N_UTL_Misc.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : NodeDevBlock
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/3/01
//-----------------------------------------------------------------------------
class NodeDevBlock
{
  friend std::ostream &operator<<(std::ostream & os, const NodeDevBlock & ndb);

public:
  NodeDevBlock()
  {}

  NodeDevBlock(
    const Device::InstanceBlock &       ib)
  : devBlock_(ib)
  {}

  NodeDevBlock(const NodeDevBlock & right)
  : devBlock_(right.devBlock_)
  {}

  ~NodeDevBlock()
  {}

private:
  // Assignment operator (autogen)
  NodeDevBlock & operator = (const NodeDevBlock & right);

public:
  // clear data for reuse of this object
  void clear();

  void addNode(const std::string&);
  std::vector<std::string> & get_NodeList();
  const std::vector<std::string> & get_NodeList() const;
  void set_NodeList(const std::vector<std::string> & nList);

  Device::InstanceBlock & getDevBlock() { return devBlock_; }
  const Device::InstanceBlock & getDevBlock() const { return devBlock_; }

protected:
  std::vector<std::string>       nodeList_;
  Device::InstanceBlock         devBlock_;
};

} // namespace Topo
} // namespace Xyce

#endif

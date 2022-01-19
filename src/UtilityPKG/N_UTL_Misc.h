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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 5/22/00
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_Misc_h
#define Xyce_N_UTL_Misc_h

#include <string>
#include <iostream>
#include <utility>
#include <vector>

#include <N_UTL_Pack.h>

namespace Xyce {
namespace Util {

std::ostream &word_wrap(std::ostream &os, const std::string &s, std::string::size_type line_length, const std::string &prefix, const std::string &prefix_first_line);

} // namespace Util

//-----------------------------------------------------------------------------
// Class         : NodeID 
// Purpose       : Container of node name and node type
// Special Notes :
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/27/2019
//-----------------------------------------------------------------------------
class NodeID : public std::pair< std::string, int >
{
  friend class Pack<NodeID>;
public:
  NodeID()
    : std::pair<std::string, int>()
  {}

  NodeID(const std::string &id, int type)
    : std::pair<std::string, int>( id, type )
  { /*std::cout << "Creating NodeID with node string: " << id << std::endl;*/ }
};

// Non-member helper functions
inline const std::string& get_node_id(const NodeID& nodeID)
{
  return nodeID.first;
}

inline void set_node_id(NodeID& nodeID, const std::string& id)
{
  nodeID.first = id;
 // std::cout << "Setting NodeID to node string: " << id << std::endl;
}

inline const int& get_node_type(const NodeID& nodeID)
{
  return nodeID.second;
}

inline void set_node_type(NodeID& nodeID, const int type)
{
  nodeID.second = type;
}

std::ostream& operator<< (std::ostream &os, const NodeID& n);


//-----------------------------------------------------------------------------
// Class         : NodeName
// Purpose       : Container of integers defining unique node name
// Special Notes :
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/27/2019 
//-----------------------------------------------------------------------------
class NodeName : public std::vector< int >
{
  public:

  NodeName()
    : std::vector<int>()
  {}
};

//-----------------------------------------------------------------------------
// Class         : StringName
// Purpose       : database of strings that are concatenated to make node names
// Special Notes :
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/27/2019 
//-----------------------------------------------------------------------------
class StringName
{
public:
  typedef std::vector<std::string> StringNameVector;

  /// Returns the filename vector
  ///
  /// The filename vector contains all the netlist filenames.
  ///
  /// @return const reference to the filename vector.
  static const StringNameVector &getStringNameVector();

  /// Returns the index with the string name
  ///
  /// @param name         string name
  ///
  /// @return index of the string in StringNameVector
  static int getStringNumber(const std::string& name);
};


} // namespace Xyce

namespace std {

template<>
struct equal_to<Xyce::NodeID> : public std::binary_function<Xyce::NodeID, Xyce::NodeID, bool>
{
  bool operator()(const Xyce::NodeID &lhs, const Xyce::NodeID &rhs) const
  {
    equal_to<std::string> x0;
    equal_to<int> x1;

    return x1(Xyce::get_node_type(lhs), Xyce::get_node_type(rhs)) && x0(Xyce::get_node_id(lhs), Xyce::get_node_id(rhs));
  }
};

template<>
struct hash<Xyce::NodeID> : public std::unary_function<Xyce::NodeID, size_t>
{
  size_t operator()(const Xyce::NodeID &node_id) const
  {
    hash<std::string> x0;
    hash<int> x1;

    return x0(Xyce::get_node_id(node_id)) ^ x1(Xyce::get_node_type(node_id));
  }
};

} // std namespace

#endif

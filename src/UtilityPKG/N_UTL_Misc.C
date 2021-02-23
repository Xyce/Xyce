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
// Purpose       :
//
// Special Notes :
//
// Creator       : Lon Waters, SNL
//
// Creation Date : 4/30/01
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>

#include <N_UTL_Misc.h>
#include <N_PDS_Comm.h>

#ifdef Xyce_PARALLEL_MPI
 #include <mpi.h>
#endif

namespace Xyce {

namespace {

std::string::const_iterator
find_next_char(
  std::string::const_iterator	p,
  std::string::const_iterator	end,
  char				c)
{
  while (p != end && *p != c)
    p++;
  return p;
}

std::string::const_iterator
find_next_not_char(
  std::string::const_iterator	p,
  std::string::const_iterator	end,
  char				c)
{
  while (p != end && *p == c)
    p++;
  return p;
}

inline std::string::const_iterator find_next_space(std::string::const_iterator p, std::string::const_iterator end) {
  return find_next_char(p, end, ' ');
}

inline std::string::const_iterator find_next_endl(std::string::const_iterator p, std::string::const_iterator end) {
  return find_next_char(p, end, '\n');
}

inline std::string::const_iterator find_next_nonspace(std::string::const_iterator p, std::string::const_iterator end) {
  return find_next_not_char(p, end, ' ');
}

} // namespace <null>

namespace Util {

std::ostream &
word_wrap(
  std::ostream &                os,
  const std::string &	        s,
  std::string::size_type        line_length,
  const std::string &	        prefix,
  const std::string &	        prefix_first_line)
{
  const std::string *u = &prefix_first_line;

  std::string::const_iterator p0, p1, p2, p3;
  p0 = p1 = p2 = s.begin();

  while (p2 != s.end() ) {

    // skip preceeding whitespace
    p1 = find_next_nonspace(p0, s.end());
    p3 = find_next_endl(p0, s.end());
    p2 = p1 = find_next_space(p1, s.end());
    do { // skip words
      p1 = find_next_nonspace(p1, s.end());
      p1 = find_next_space(p1, s.end());
      if (p3 < p1) {
	p2 = p3;
	break;
      }
      if (p1 - p0 + u->size() > line_length) // hit word end past line_length
	break;
      p2 = p1;
    } while (p2 != s.end());

    os << *u << std::string(p0, p2) << "\n";

    // if (p2 == p3) // If you want an embedded newline to mean
    //   u = &prefix_first_line;
    // else
    u = &prefix;

    p0 = p2 + 1;
  }

  return os;
}

} // namespace Util

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/11/01
//-----------------------------------------------------------------------------
std::ostream& operator<< (std::ostream &os, const NodeID& n)
{
  return os << "( " << Xyce::get_node_id(n) << " , " << Xyce::get_node_type(n) << " )";
}

//-----------------------------------------------------------------------------
// Function      : Node::packedByteCount
// Purpose       : Constructor
// Scope         : public
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/01
//-----------------------------------------------------------------------------
template<>
int Pack<NodeID>::packedByteCount(
  const NodeID &    node)
{
  int count = 0;

  count += sizeof(int); //Node name length
  count += Xyce::get_node_id(node).length(); //Node name

  count += sizeof(int); //Node type

  return count;
}

//-----------------------------------------------------------------------------
// Function      : Node::pack
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/11/01
//-----------------------------------------------------------------------------
template<>
void
Pack<NodeID>::pack(
  const NodeID &           node,
  char *                   buf,
  int                      bsize,
  int &                    pos,
  Parallel::Communicator * comm ) 
{
  int length = Xyce::get_node_id(node).length();
  comm->pack( &length, 1, buf, bsize, pos );

  comm->pack( Xyce::get_node_id(node).c_str(), length, buf, bsize, pos );

  length = Xyce::get_node_type(node);
  comm->pack( &length, 1, buf, bsize, pos );
}

//-----------------------------------------------------------------------------
// Function      : Node::unpack
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/11/01
//-----------------------------------------------------------------------------
template<>
void
Pack<NodeID>::unpack(
  NodeID &                  node,
  char *                    buf,
  int                       bsize,
  int &                     pos,
  Parallel::Communicator *  comm )
{
  int length;
  comm->unpack( buf, bsize, pos, &length, 1 );

  Xyce::set_node_id( node, std::string( (buf+pos), length ) );
  pos += length;

  comm->unpack( buf, bsize, pos, &length, 1 );
  Xyce::set_node_type( node, length );
}


// This is a simple struct to hold the string name vector.
struct StringData
{
  StringName::StringNameVector    stringNameVector_;             ///<Filename registration

  ~StringData() {}
};

//-----------------------------------------------------------------------------
// Function      : getStringData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/27/2019
//-----------------------------------------------------------------------------
/// returns the configuration data singleton.
///
/// @return reference to the configuration data singleton.
///
StringData &getStringData()
{
  static StringData data_;

  return data_;
}

//-----------------------------------------------------------------------------
// Function      : StringName::getStringNameVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/27/2019
//-----------------------------------------------------------------------------
const StringName::StringNameVector &
StringName::getStringNameVector()
{
  return getStringData().stringNameVector_;
}

//-----------------------------------------------------------------------------
// Function      : StringName::getStringNumber
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/27/2019
//-----------------------------------------------------------------------------
int
StringName::getStringNumber(const std::string& name)
{ 
  StringNameVector::iterator it;
  StringNameVector& stringVector = getStringData().stringNameVector_;
  it = std::find(stringVector.begin(), stringVector.end(), name);
  if ( it != stringVector.end() )
  { 
    return std::distance( stringVector.begin(), it );
  }
  else
  { 
    stringVector.push_back( name );
    return (stringVector.size()-1);
  }
}

} // namespace Xyce

//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Creation Date  : 07/05/01
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef N_TOP_Directory_h
#define N_TOP_Directory_h 1

#include <string>
#include <vector>
#include <iosfwd>

#include <N_TOP_fwd.h>
#include <N_PDS_fwd.h>
#include <N_UTL_Misc.h>

#include <N_TOP_Node.h>
#include <N_TOP_ParNode.h>

namespace Xyce {
namespace Topo {

struct DirectoryData;

//-----------------------------------------------------------------------------
// Class         : Directory
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/05/01
//-----------------------------------------------------------------------------
class Directory
{

public:

  // Constructor
  Directory( Topology &topology, Parallel::Communicator &pds_comm)
  : topology_(topology),
    pdsComm_(pds_comm),
    data_(0)
  { }

  // Destructor
  ~Directory();

  bool generateDirectory();

  std::vector<NodeID> getProcs( const std::vector<NodeID> & idVec,
                                std::vector<int> & procVec );

  std::vector<NodeID> getSolnGIDs( const std::vector<NodeID> & idVec,
                                   std::vector< std::vector<int> > & gidVec,
                                   std::vector<int> & procVec );

private:

  // Copy constructor (private).
  Directory(const Directory & right);

  // Assignment operator (private).
  Directory & operator=(const Directory & right);

  // Equality operator.
  bool operator==(const Directory & right) const;

  // Non-equality operator.
  bool operator!=(const Directory & right) const;

private:

  // Pointer to the topology manager.
  Topology &                    topology_;

  // Pointer to the PDS manager.
  Parallel::Communicator &      pdsComm_;

  // Pimpl
  DirectoryData *               data_;
};

} // namespace Topo
} // namespace Xyce

#endif

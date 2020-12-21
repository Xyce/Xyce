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

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/12/02
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef N_TOP_Indexor_h
#define N_TOP_Indexor_h 1

#include <vector>
#include <map>

#include <N_UTL_Misc.h>

#include <N_PDS_fwd.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : Indexor
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/12/02
//-----------------------------------------------------------------------------
class Indexor
{

public:

  // Constructor
  Indexor( Parallel::Manager &pds)
  : pdsMgr_(pds),
    accelMatrixIndex_(false)
  { }

  // Destructor
  ~Indexor() {}

  bool globalToLocal( int graph_id, std::vector<int> & ids );
  bool localToGlobal( int graph_id, std::vector<int> & ids );

  bool setupAcceleratedMatrixIndexing( int graph_id );
  bool deleteAcceleratedMatrixIndexing();

  bool matrixGlobalToLocal( int graph_id, const std::vector<int> & gids, std::vector< std::vector<int> > & stamp );

private:

  // Copy constructor (private).
  Indexor(const Indexor & right);
  // Assignment operator (private).
  Indexor & operator=(const Indexor & right);
  // Equality operator.
  bool operator==(const Indexor & right) const;
  // Non-equality operator.
  bool operator!=(const Indexor & right) const;

private:

  // Pointer to the PDS manager.
  Parallel::Manager &           pdsMgr_;
  bool                          accelMatrixIndex_;
  std::vector< std::map<int,int> > matrixIndexMap_;

};

} // namespace Topo
} // namespace Xyce

#endif

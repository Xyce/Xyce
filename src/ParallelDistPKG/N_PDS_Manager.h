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
// Purpose        : Manager for the parallel load-balance and distribution tools.
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/08/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_PDS_Manager_h
#define Xyce_N_PDS_Manager_h

#include <vector>
#include <string>
#include <map>

#include <N_PDS_fwd.h>
#include <N_LAS_fwd.h>

#ifdef Xyce_PARALLEL_MPI
#include <mpi.h>
#endif

namespace Xyce {
namespace Parallel {

/// Local to global and inverse maps.
enum {
  NODE,                                 ///< 

  SOLUTION,                             ///< Solution variables
  SOLUTION_OVERLAP,                     ///< Solution variables on processor boundary
  SOLUTION_OVERLAP_GND,                 ///< Solution variables on processor boundary with ground
  STATE,                                ///< State variables
  STORE,                                ///< Store variables
  LEADCURRENT,                          ///< Lead current variables
  JACOBIAN,                             ///< Jacobian values
  JACOBIAN_OVERLAP,                     ///< Jacobian values on processor boundary

  LINEAR_SOLUTION,                      ///< Only solution values for linear devices
  LINEAR_JACOBIAN,                      ///< Only Jacobian rows for linear devices, linear map
  GLOBAL_LINEAR_JACOBIAN,               ///< Only Jacobian rows for linear devices, global map
  NONLINEAR_SOLUTION,                   ///< Only solution values for nonlinear devices
  NONLINEAR_JACOBIAN,                   ///< Only Jacobian rows for nonlinear devices, nonlinear map
  GLOBAL_NONLINEAR_JACOBIAN,            ///< Only Jacobian rows for nonlinear devices, global map
  LIN_NONLIN_JACOBIAN,                  ///< Only Jacobian rows for linear nodes connected to nonlinear, linear map
  GLOBAL_LIN_NONLIN_JACOBIAN,           ///< Only Jacobian rows for linear nodes connected to nonlinear, global map
  NONLIN_LIN_JACOBIAN,                  ///< Only Jacobian rows for nonlinear nodes connected to linear, nonlinear map
  GLOBAL_NONLIN_LIN_JACOBIAN,           ///< Only Jacobian rows for nonlinear nodes connected to linear, global map

  MAP_COUNT                             ///< Number of maps
};

//-----------------------------------------------------------------------------
// Class         : Manager
// Purpose       : Manager for the parallel load-balance and distribution tools.
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Compuational Sciences
// Creation Date : 03/08/00
//-----------------------------------------------------------------------------
class Manager
{
public:
  Manager(int iargs, char **cargs, Xyce::Parallel::Machine comm);
  ~Manager();

private:
  Manager(const Manager &);
  Manager &operator=(const Manager &);

public:
  // Return ptr to common Comm object owned by Mgr.
  Communicator * getPDSComm() const
  {
    return pdsComm_;
  }

  ParMap *getParallelMap(int id) {
    return parMaps_[id];
  }

  // Gets a global accessor object which is associated with a parallel map.
  GlobalAccessor * getGlobalAccessor(int id)
  {
    return globalAccessors_[id];
  }

    // Gets a matrix graph object.
  const Linear::Graph * getMatrixGraph(int id) const
  {
    return matrixGraphs_[id];
  }

  bool addParallelMap(int id, ParMap * map);

  bool linkParallelMap(int new_id, int link_id);

  bool deleteParallelMap(int id);

  // Add a global accessor object associated with a parallel map.
  GlobalAccessor *addGlobalAccessor(int id);

  // Delete a global accessor object associated with a parallel map.
  bool deleteGlobalAccessor(int id);

  // Method which greats a global accessor object.
  GlobalAccessor * createGlobalAccessor();

  // Add a matrix graph object
  bool addMatrixGraph(
    int                         id,
    Linear::Graph *             graph );

  bool linkMatrixGraph(int new_id, int link_id);

  // Delete a matrix graph object
  bool deleteMatrixGraph(int id);

private:
  Communicator *                pdsComm_;               ///< Pointer to the PDS_Comm object.
  ParMap *                      parMaps_[MAP_COUNT];
  GlobalAccessor *              globalAccessors_[MAP_COUNT];
  Linear::Graph *               matrixGraphs_[MAP_COUNT];

  std::map<int,int>             linkedMapsGraphs_;
};

} // namespace Parallel
} // namespace Xyce

#endif

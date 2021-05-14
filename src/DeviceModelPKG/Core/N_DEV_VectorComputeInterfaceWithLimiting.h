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

//-----------------------------------------------------------------------------
//
// Purpose        : Abstract interface for general "vector loading"
//                  callback classes with limiting.
//
//
// Special Notes  : This can be used to implement a generic call-back interface
//                  for external simulators.  The GeneralExternal device
//                  will use this to call back to the external simulator to
//                  get the DAE vector and matrix elements associated with
//                  the device. This version enables limiting.
//
// Creator        : Paul Kuberry, SNL
//
// Creation Date  : 2/19/2021
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef N_DEV_VectorComputeInterfaceWithLimiting_H
#define N_DEV_VectorComputeInterfaceWithLimiting_H
#include "N_DEV_VectorComputeInterface.h"

namespace Xyce {
namespace Device {

class VectorComputeInterfaceWithLimiting : public VectorComputeInterface
{
public:
  VectorComputeInterfaceWithLimiting(){};
  virtual ~VectorComputeInterfaceWithLimiting(){};

  /// Call back function for time domain
  ///
  /// For VectorComputeInterfaceWithLimiting, they get their 
  /// computeXyceVectors from computeXyceVectorsWithLimiting with 
  /// arguments populated that will not be used.
  ///
  /// This enables the user to only define computeXyceVectorsWithLimiting
  /// and provide limiting + no limiting functionality.
  ///
  virtual bool computeXyceVectors(std::vector<double> & solutionVars,
                                  double time,
                                  std::vector<double>&F,
                                  std::vector<double>&Q,
                                  std::vector<double>&B,
                                  std::vector<std::vector<double> > &dFdx,
                                  std::vector<std::vector<double> > &dQdx) { return false; }


  /// Call back function for time domain
  ///
  /// Xyce will call this function at every newton iteration of every
  /// time step, and the external code must fill in F, Q, B,
  /// dFdx, dQdx, dFdxdVp, and dQdxdQp as appropriate.
  ///
  /// This must always be implemented.  
  ///
  /// @param[in] solutionVar   Solution variables
  /// @param[in] flagSolutionVars   Flag solution variables
  /// @param[in/out] storeVars   Store variables
  /// @param[in] time current value of time
  /// @param[in] deviceOptions   Options for device instance
  /// @param[in] solverState   Solver state
  /// @param[in/out] origFlag   Flag as to whether limited variable was modified
  /// @param[out] F   F vector
  /// @param[out] Q   Q vector
  /// @param[out] B   B vector
  /// @param[out] dFdx   Derivative of F with respect to solution
  /// @param[out] dQdx   Derivative of Q with respect to solution
  /// @param[out] dFdxdVp   Derivative of F with respect to limited solution
  /// @param[out] dQdxdVp   Derivative of Q with respect to limited solution
  ///
  virtual bool computeXyceVectorsWithLimiting(std::vector<double> & flagSolutionVars,
                                  std::vector<std::vector<double> > & solutionVars,
                                  std::vector<std::vector<double> > & storeVars,
                                  std::vector<std::vector<double> > & stateVars,
                                  const DeviceOptions & deviceOptions,
                                  const SolverState & solverState,
                                  bool & origFlag,
                                  std::vector<double>&F,
                                  std::vector<double>&Q,
                                  std::vector<double>&B,
                                  std::vector<std::vector<double> > &dFdx,
                                  std::vector<std::vector<double> > &dQdx,
                                  std::vector<double>&dFdxdVp,
                                  std::vector<double>&dQdxdVp) = 0; 

};

} // namespace Device
} // namespace Xyce
#endif

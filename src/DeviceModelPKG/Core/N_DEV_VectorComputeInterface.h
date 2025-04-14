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
// Purpose        : Abstract interface for general "vector loading"
//                  callback classes.
//
//
// Special Notes  : This can be used to implement a generic call-back interface
//                  for external simulators.  The GeneralExternal device
//                  will use this to call back to the external simulator to
//                  get the DAE vector and matrix elements associated with
//                  the device.
//
// Creator        : Tom Russo, SNL, Electrical Models and Simulation
//
// Creation Date  : 2/28/2017
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef N_DEV_VectorComputeInterface_H
#define N_DEV_VectorComputeInterface_H
#include <string>
#include <vector>
#include <complex>

namespace Xyce {
namespace Device {

/// Abstract interface class for allowing the GeneralExternal device
/// to obtain DAE vectors and matrices (e.g. from an external simulator
/// or otherwise)
///
/// Xyce solves the DAE
///   \f$
///     F(X) + \frac{dQ(X)}{dt} - B(t) = 0
///   \f$
/// where \f$X\f$ is the solution vector and \f$t\f$ is the time.
///
/// The equation is set up using modified nodal analysis, and as such
/// the total of all current contributions into nodes must sum to zero.
/// Additional equations may be needed for things like voltage sources.
/// Consult the Xyce Math Document for more details.
///
/// Each device is responsible for loading the components of F, Q, and
/// B.  By implementing this class and assigning an instance of it it
/// as the vector loader for a General External device (YGENEXT), an
/// external program can insert the results of an arbitrary
/// computation into the circuit equations.
///
/// Provisions are made for coupling to both time and frequency domain.
///
class VectorComputeInterface
{
public:
  VectorComputeInterface(){};
  virtual ~VectorComputeInterface(){};

  /// Call back function for time domain
  ///
  /// Xyce will call this function at every newton iteration of every
  /// time step, and the external code must fill in F, Q, B,
  /// dFdx, and dQdx as appropriate.
  ///
  /// This must always be implemented.  If you are implementing a
  /// *STRICTLY* frequency domain device, then simply implement this
  /// function to do nothing but return FALSE.  That will cause any
  /// attempt to use your model in time domain to fail.
  ///
  ///
  /// @param[in] sV   Solution variables
  /// @param[in] time current value of time
  /// @param[out] F   F vector
  /// @param[out] Q   Q vector
  /// @param[out] B   B vector
  /// @param[out] dFdx  Derivative of F with respect to solution
  /// @param[out] dQdx  Derivative of Q with respect to solution
  ///
  virtual bool computeXyceVectors(std::vector<double> & solutionVars,
                                  double time,
                                  std::vector<double>&F,
                                  std::vector<double>&Q,
                                  std::vector<double>&B,
                                  std::vector<std::vector<double> > &dFdx,
                                  std::vector<std::vector<double> > &dQdx) = 0;

  /// Call back function for frequency domain
  ///
  /// Xyce will call this function for every newton iteration of the nonlinear solve
  /// for every frequency, and the external code must fill in F, B,
  /// and dFdx as appropriate.
  ///
  /// This need not be implemented unless your model has a separate frequency
  /// domain evaluation.
  ///
  ///
  /// @param[in] sV   Solution variables
  /// @param[in] frequency current value of frequency
  /// @param[out] F   F vector
  /// @param[out] B   B vector
  /// @param[out] dFdx  Derivative of F with respect to solution
  ///
  /// @note While it is intended that this function could be used to
  /// implement device models that are nonlinear in the frequency
  /// domain, in the current implementation of the harmonic balance
  /// loader, it is explicitly assumed that all frequency domain
  /// models are linear.  Because of that upstream assumption, in
  /// reality this function is only called once per frequency at the
  /// beginning of a harmonic computation.  It is then assumed that
  /// the dFdx matrices returned are constant for the duration of the
  /// simulation.  Use of nonlinear frequency domain models is not yet
  /// supported at higher levels, even though the design of this class
  /// should support it when the upper level support is added.
  virtual bool computeXyceFDVectors(std::vector<std::complex<double> > & solutionVars,
                                   double frequency,
                                   std::vector<std::complex<double> >&F,
                                   std::vector<std::complex<double> >&B,
                                   std::vector<std::vector<std::complex<double> > > &dFdx)
  {
    return true;
  };

  /// query function for frequency domain loader
  ///
  /// If you implement computXyceFDVectors to return frequency domain information
  /// you must reimplement this to return "true"
  ///
  /// Xyce will call computeXyceVectors with time-domain solution
  /// vectors and Fourier transform the results when doing harmonic
  /// balance analysis unless this function returns true.  If true,
  /// Xyce will instead call the computeXyceFDVectors method when
  /// doing harmonic balance analysis, and will only call
  /// computeXyceVectors for time domain and small-signal AC analysis.
  ///
  virtual bool haveFDLoads() { return false;};

};

} // namespace Device
} // namespace Xyce
#endif

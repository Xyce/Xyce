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
// Purpose        : Neuron classes.
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 08/11/10
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MembraneHH_h
#define Xyce_N_DEV_MembraneHH_h

#include <N_DEV_MembraneModel.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : MembraneHH
// Purpose       : This is class defines a membrane with Hodgkin-Huxley equations
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
class MembraneHH : public MembraneModel
{
public:
  MembraneHH(const SolverState & ss1, double cMem, double gMem, double vRest, double eK, double gK, double eNa, double gNa);
  ~MembraneHH() {}

  void setJacStamp( int numExtVars, int segmentNumber, int vOffset, std::vector< std::vector< int > > & segmentJacStamp );
  void loadDAEQVector( int segmentNumber, std::vector< int > & lidIndexVector, Linear::Vector * solnVecPtr, Linear::Vector * daeQVecPtr, double segArea);
  void loadDAEFVector( int segmentNumber, std::vector< int > & lidIndexVector, Linear::Vector * solnVecPtr, Linear::Vector * daeFVecPtr, double segArea);
  void loadDAEdQdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector, std::vector< std::vector< int > > & jacobianOffsets, Linear::Vector * solnVecPtr, Linear::Matrix * dQdxMatPtr, double segArea);
  void loadDAEdFdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector, std::vector< std::vector< int > > & jacobianOffsets, Linear::Vector * solnVecPtr, Linear::Matrix * dFdxMatPtr, double segArea);

  double cMem_;     // membrane capacitance
  double gMem_;     // membrane conductance
  double vRest_;    // membrane rest voltage
  double eK_;
  double gK_;
  double eNa_;
  double gNa_;
};

} // namespace Device
} // namespace Xyce

#endif

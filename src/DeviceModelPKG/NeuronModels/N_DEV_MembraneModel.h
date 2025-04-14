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

#ifndef Xyce_N_DEV_MembraneModel_h
#define Xyce_N_DEV_MembraneModel_h

#include <vector>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_LAS_fwd.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : N_DEV_MembraneModel
// Purpose       : This is class is a virtual base class to define the base
//                 interface for membrane modeling equations.  Classes
//                 that derive from this will define the equations that
//                 need to be solved for various ion currents
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
class MembraneModel
{
public:
  MembraneModel(const SolverState & ss1)
    : numIndependentVars_(0),
      numExternalVars_(2),
      solState(ss1)
  {}
    
  ~MembraneModel() {}

  int numVars() { return numIndependentVars_; }

  virtual void setJacStamp( int numExtVars, int segmentNumber, int vOffset, std::vector< std::vector< int > > & segmentJacStamp ) {}
  virtual void loadDAEQVector( int segmentNumber, std::vector< int > & lidIndexVector, Linear::Vector * solnVecPtr, Linear::Vector * daeQVecPtr, double segArea) {}
  virtual void loadDAEFVector( int segmentNumber, std::vector< int > & lidIndexVector, Linear::Vector * solnVecPtr, Linear::Vector * daeFVecPtr, double segArea) {}
  virtual void loadDAEdQdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector, std::vector< std::vector< int > > & jacobianOffsets, Linear::Vector * solnVecPtr, Linear::Matrix * dQdxMatPtr, double segArea) {}
  virtual void loadDAEdFdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector, std::vector< std::vector< int > > & jacobianOffsets, Linear::Vector * solnVecPtr, Linear::Matrix * dFdxMatPtr, double segArea) {}

  int numIndependentVars_;

  // these are modeling parameters common to most model types
  // this is some data duplication between the devcie model the membrane model
  // need to think of a cleaner way to do this
  double cMem_;
  double gMem_;
  double vRest_;

  const int numExternalVars_;  // always assume that the owning cable equation has two external vars (a Vin and Vout)
  // if this assumption needs to be changed we'll have this constant to show where we made such an assumption.

  const SolverState & solState;  // this is here incase a model needs to be aware of simulator data like dcopflag or currtime.
};

} // namespace Device
} // namespace Xyce

#endif


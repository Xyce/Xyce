//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose        : Container class for linear system Jacobian, RHS, Soln Vec,
//                  Error Vec, and Creators N_LAS_QueryUtil and N_PDS_ParMap
//                  are registered to minimize input for matrix and vector
//                  creation
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/10/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef  Xyce_LAS_System_h
#define  Xyce_LAS_System_h

#include <list>

#include <N_ANP_fwd.h>
#include <N_LAS_fwd.h>

// to eliminate RCP warnings, putting N_PDS_Manager header here.
#include <N_PDS_Manager.h>

// ---------- Forward Declarations ----------

class Epetra_MapColoring;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : System
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
class System
{

public:
  System();

  ~System();

  // Registration methods for necessary utilities

  // Registers the Parallel Distribution Services (PDS) manager (i.e., sets the
  // pointer.
  bool registerPDSManager(N_PDS_Manager * PDS_Manager)
  { return ( pdsMgr_ = PDS_Manager ); }

  // Registers the Query utility
  bool registerQueryUtil(QueryUtil * LAS_QUtil)
  { return ( lasQueryUtil_ = LAS_QUtil ); }

  // Registers the LAS Builder object
  bool registerBuilder(Builder * builder)
  { return ( lasBuilder_ = builder ); }

  // Register pointer to device mask
  bool registerDeviceMaskVector(Vector * dmPP)
  { return (deviceMaskVectorPtr_ = dmPP); }

  // Generate the solution and IC coloring vectors using the builder.  
  bool generateSolnICColoring();

  // Accessor methods for linear system objects

  // Get method for the Jacobian matrix
  Matrix *  getJacobianMatrix() { return jacobianMatrixPtr_; }

  // Get method for the residual (RHS) vector
  Vector *  getRHSVector() { return rhsVectorPtr_; }

  // Get method for the Newton update vector
  Vector *  getNewtonVector() { return newtonVectorPtr_; }

  // Get method for the linear problem consisting of the Jacobian matrix,
  // Newton update, and RHS vector.
  Problem * getLinearProblem() { return lasProblemPtr_; }

  // Get method for the limiting (jdxp) vector
  Vector *  getJDXPVector() { return jdxpVectorPtr_; }

  // dFdxdvp vector, for new-DAE voltage limiting.
  Vector * getdFdxdVpVector () { return  dFdxdVpVectorPtr_; }

  // dQdxdvp vector, for new-DAE voltage limiting.
  Vector * getdQdxdVpVector () { return  dQdxdVpVectorPtr_; }

  // Get method for flag solution vector
  Vector * getFlagSolVector() { return flagSolVectorPtr_; }

  // Get method for the device mask vector
  Vector *  getDeviceMaskVector() { return deviceMaskVectorPtr_; }

  // Registers the Query utility
  QueryUtil * getQueryUtil()
  {
    return lasQueryUtil_;
  }

  // Create residual (RHS) and Jacobian
  bool initializeSystem(); // bool block_analysis_flag);

  // Builder access
  Builder & builder() { return *lasBuilder_; }
  const Builder & builder() const { return *lasBuilder_; }
 
  bool updateExternValsSolnVector(MultiVector * solnVector);

private:

  N_PDS_Manager   *  pdsMgr_;

  QueryUtil *  lasQueryUtil_;
  Builder   *  lasBuilder_;
  Problem   *  lasProblemPtr_;

  Matrix    *  jacobianMatrixPtr_;
  Vector    *  rhsVectorPtr_;
  Vector    *  newtonVectorPtr_;
  Vector    *  jdxpVectorPtr_;

  Vector    *  dFdxdVpVectorPtr_;
  Vector    *  dQdxdVpVectorPtr_;

  Vector    *  flagSolVectorPtr_;

  Vector    *  deviceMaskVectorPtr_;

  Epetra_MapColoring * solnColoringPtr_;

  Epetra_MapColoring * initialConditionColoringPtr_;
};

} // namespace Linear
} // namespace Xyce

#endif // Xyce_LAS_System_h


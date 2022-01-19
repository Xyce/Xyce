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

//-----------------------------------------------------------------------------
//
// Purpose        : Container class for linear system
//                  Jacobian, RHS, Soln Vec, Error Vec, and Creators
//                  N_LAS_QueryUtil and N_PDS_ParMap are registered to
//                  minimize input for matrix and vector creation
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

#include <Xyce_config.h>

#include <iostream>

#include <N_ANP_AnalysisManager.h>
#include <N_LAS_Builder.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_QueryUtil.h>
#include <N_LAS_System.h>
#include <N_LAS_SystemHelpers.h>
#include <N_LAS_Problem.h>
#include <N_LAS_Vector.h>
#include <N_PDS_GlobalAccessor.h>
#include <N_PDS_Manager.h>
#include <N_PDS_ParMap.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : System::System
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Oct  6 07:33:03 2014
//-----------------------------------------------------------------------------
System::System()
  : pdsMgr_(0),
    lasBuilder_(0),
    lasProblemPtr_(0),
    jacobianMatrixPtr_(0),
    rhsVectorPtr_(0),
    newtonVectorPtr_(0),
    jdxpVectorPtr_(0),
    dFdxdVpVectorPtr_(0),
    dQdxdVpVectorPtr_(0),
    flagSolVectorPtr_(0),
    deviceMaskVectorPtr_(0)
  {}

//-----------------------------------------------------------------------------
// Function      : System::~System
// Purpose       : destructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
System::~System()
{
  delete jacobianMatrixPtr_;
  delete jdxpVectorPtr_;
  delete rhsVectorPtr_;
  delete newtonVectorPtr_;
  delete lasProblemPtr_,
  delete dFdxdVpVectorPtr_;
  delete dQdxdVpVectorPtr_;
  delete flagSolVectorPtr_;
}

//-----------------------------------------------------------------------------
// Function      : System::numGlobalRows()
// Purpose       : Return the number of global rows from the builder.
// Special Notes :
// Scope         : Public
// Creator       : Heidi K. Thornquist, SNL
// Creation Date : 9/10/2019
//-----------------------------------------------------------------------------
// Return number of global rows as provided by the builder.
int System::numGlobalRows() const
{
  return lasBuilder_->getSolutionMap()->numGlobalEntities();
}

//-----------------------------------------------------------------------------
// Function      : System::initializeSystem
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
bool System::initializeSystem()
{
  bool bSuccess = true;

  rhsVectorPtr_ = lasBuilder_->createVector();
  newtonVectorPtr_ = lasBuilder_->createVector();
  jdxpVectorPtr_ = lasBuilder_->createVector();
  jacobianMatrixPtr_ = lasBuilder_->createMatrix();

  // If this is a matrix free analysis, there is no need to create a linear problem here.
  if (jacobianMatrixPtr_ != 0)
  { 
    lasProblemPtr_ = Xyce::Linear::createProblem( jacobianMatrixPtr_, newtonVectorPtr_, rhsVectorPtr_ );
  }

  // these are needed for new-DAE:
  dFdxdVpVectorPtr_ = lasBuilder_->createVector();
  dQdxdVpVectorPtr_ = lasBuilder_->createVector();

  flagSolVectorPtr_ = lasBuilder_->createVector();

  return bSuccess;
}

//-----------------------------------------------------------------------------
// Function      : System::updateExternValsSolnVector
// Purpose       : updates off proc values of soln vectors
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/18/00
//-----------------------------------------------------------------------------
bool System::updateExternValsSolnVector(MultiVector * solnVector)
{
#ifdef Xyce_PARALLEL_MPI
  Parallel::GlobalAccessor * Accessor = pdsMgr_->getGlobalAccessor( Parallel::SOLUTION );
  if (Accessor)
    Accessor->migrateMultiVector(solnVector);
#endif

  return true;
}

} // namespace Linear
} // namespace Xyce

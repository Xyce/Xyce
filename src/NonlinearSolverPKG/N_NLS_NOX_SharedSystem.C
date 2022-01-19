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

//-------------------------------------------------------------------------
//
// Purpose        : Interface to Xyce vectors for NOX.
//
// Special Notes  :
//
// Creator        : Tammy Kolda, NLS, 8950
//
// Creation Date  : 01/31/02
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include "N_NLS_NOX_SharedSystem.h"
#include "N_NLS_NOX_Interface.h"
#include "N_LAS_Vector.h"
#include "N_LAS_Matrix.h"
#include "N_LAS_System.h"
#include "N_LAS_Builder.h"
#include "N_ERH_ErrorMgr.h"

// ----------   NOX Includes   ----------

// ---------- Namespaces ---------------
namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

//-----------------------------------------------------------------------------
// Function      : SharedSystem::SharedSystem
// Purpose       : Constructor. Creates a shared system containing
//                 the soln vector, the previous solution vector,
//                 the RHS vector, the Newton vector, the Jacobian
//                 matrix, and a reference to the interface used to
//                 call the evaluation functions.
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
SharedSystem::SharedSystem(Linear::Vector& soln,
			   Linear::Vector& f,
			   Linear::Matrix& jacobian,
			   Linear::Vector& newton,
			   Linear::Vector& gradient,
			   Linear::System& lasSys,
			   Interface& interface) :

  xyceSolnPtr_(0),
  xyceFPtr_(0),
  xyceJacobianPtr_(0),
  xyceNewtonPtr_(0),
  xyceGradientPtr_(0),
  xyceLasSysPtr_(0),
  xyceInterfacePtr_(0),
  matrixFreeFlag_(false),
  ownerOfJacobian_(0),
  ownerOfStateVectors_(0)
{
  reset(soln, f, jacobian, newton, gradient, lasSys, interface);
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::~SharedSystem
// Purpose       : Destructor.
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
SharedSystem::~SharedSystem()
{
  delete xyceSolnPtr_;
  delete xyceFPtr_;
  delete xyceNewtonPtr_;
  delete xyceGradientPtr_;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::reset
// Purpose       : reset the Xyce fill objects - pointers may have changed!
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void SharedSystem::reset(Linear::Vector& x,
			 Linear::Vector& f,
			 Linear::Matrix& jacobian,
			 Linear::Vector& newton,
			 Linear::Vector& gradient,
			 Linear::System& lasSys,
			 Interface& interface)
{
  // Clear out old views
  delete xyceSolnPtr_;
  delete xyceFPtr_;
  delete xyceNewtonPtr_;
  delete xyceGradientPtr_;

  xyceJacobianPtr_ = &jacobian;
  xyceLasSysPtr_ = &lasSys;
  xyceInterfacePtr_ = &interface;
  matrixFreeFlag_ = xyceInterfacePtr_->getMatrixFreeFlag();

  // Create views of the data used for fills in xyce
  xyceSolnPtr_ = new Vector(x, lasSys);
  xyceFPtr_ = new Vector(f, lasSys);
  xyceNewtonPtr_ = new Vector(newton, lasSys);
  xyceGradientPtr_ = new Vector(gradient, lasSys);
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::computeF
// Purpose       : Compute the F corresponding to the current
//                 primary solution vector. Makes the primary
//                 solution vector owner in to the owner of the F.
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::computeF(const Vector& solution, Vector& F,
			    const Group* grp)
{
  ownerOfStateVectors_ = grp;

  *xyceSolnPtr_ = solution;
  bool status = xyceInterfacePtr_->computeF();

  if (status == false) {
    Report::DevelFatal0().in("SharedSystem::computeF")
      << "compute F failed!";
  }

  F = *xyceFPtr_;
  return status;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::computeJacobian
// Purpose       : Compute the Jacobian corresponding to the current
//                 primary solution vector. 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::computeJacobian(Group* grp)
{
  ownerOfJacobian_ = grp;

  // If this is matrix free, there is no Jacobian to compute.
  if (matrixFreeFlag_)
  {
    return true;
  }

  *xyceSolnPtr_ = grp->getX();

  if (!areStateVectors(grp)) {
#ifdef Xyce_VERBOSE_NOX
    if (1) { //RPP: Need to add priting utilities to group ctor
      dout() << "Warning: SharedSystem::computeJacobian() - State "
	   << "Vectors are not valid wrt solution!" << std::endl;
      dout() << "Calling computeResidual to fix this!" << std::endl;
    }
#endif
    // RPP: This line is not needed since we now call the group
    //ownerOfStateVectors_ = grp;
    
    NOX::Abstract::Group::ReturnType status = grp->computeF();

    if (status != NOX::Abstract::Group::Ok) {
      Report::DevelFatal0().in("SharedSystem::computeJacobian")
        << "compute F failed!";
    }

  }

  bool status = xyceInterfacePtr_->computeJacobian();

  if (status == false) {
      Report::DevelFatal0().in("SharedSystem::computeJacobian")
        << "SharedSystem::computeJacobian() - compute Jac failed!";
  }

  return status;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::computeNewton
// Purpose       : Compute the Newton corresponding to the current
//                 primary solution vector. 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::computeNewton(const Vector& F, Vector& Newton,
				 Teuchos::ParameterList& params)
{
  *xyceFPtr_ = F;
  // Zero out the Newton vector
  xyceNewtonPtr_->scale(0.0);
  bool status = xyceInterfacePtr_->computeNewton(params);
  Newton = *xyceNewtonPtr_;

  return status;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::applyJacobian
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::applyJacobian(const Vector& input, Vector& result) const
{
  if (!matrixFreeFlag_) {
    bool NoTranspose = false;
    xyceJacobianPtr_->matvec(NoTranspose, input.getNativeVectorRef(), result.getNativeVectorRef());
  } else {
    bool status = xyceInterfacePtr_->applyJacobian(input.getNativeVectorRef(), result.getNativeVectorRef());
    if (status == false) {
      Report::DevelFatal0().in("SharedSystem::applyJacobian")
        << "apply Jac failed!";
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::applyJacobianTranspose
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::applyJacobianTranspose(const Vector& input, Vector& result) const
{
  if (!matrixFreeFlag_) {
    bool Transpose = true;
    xyceJacobianPtr_->matvec(Transpose, input.getNativeVectorRef(), result.getNativeVectorRef());
  } else {
    Report::DevelFatal0().in("SharedSystem::applyJacobianTranspose")
      << "Not Supported for Matrix Free Loads!";
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::computeDfDpMulti	
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool SharedSystem::computeDfDpMulti	
  (const std::vector< int > & paramIDs, 
   NOX::Abstract::MultiVector & dfdp, 
   bool isValidF)
{
  bool status = xyceInterfacePtr_->computeDfDpMulti	
                   (paramIDs, dfdp, isValidF);
    
  return status;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::getSolutionVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
Vector& SharedSystem::getSolutionVector()
{
  return *xyceSolnPtr_;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::getJacobian
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
const Linear::Matrix& SharedSystem::getJacobian() const
{
  return *xyceJacobianPtr_;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::getJacobian
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
Linear::Matrix& SharedSystem::getJacobian(const Group* grp)
{
  ownerOfJacobian_ = grp;
  return *xyceJacobianPtr_;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::getStateVectors
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void SharedSystem::getStateVectors(const Group* grp)
{
  ownerOfStateVectors_ = grp;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::getLasSystem
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
Linear::System* SharedSystem::getLasSystem()
{
  return xyceLasSysPtr_;
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::cloneSolutionVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
Vector* SharedSystem::cloneSolutionVector() const
{
  Vector* tmpVectorPtr = 0;
  tmpVectorPtr = 
    dynamic_cast<Vector*>(xyceSolnPtr_->clone(NOX::DeepCopy).release().get());

  if (tmpVectorPtr == 0) {
    Report::DevelFatal0().in("SharedSystem::cloneSolutionVector")
      << " dynamic cast/ memory allocation failure!";
  }

  return (tmpVectorPtr);
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem:: getNewtonVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
const Vector & SharedSystem::getNewtonVector() const
{
  return *xyceNewtonPtr_;                                                       
}                                                                               

//-----------------------------------------------------------------------------
// Function      : SharedSystem::debugOutput1
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void SharedSystem::debugOutput1
   (Linear::Matrix & jacobian, Linear::Vector & rhs)
{
  xyceInterfacePtr_->debugOutput1(jacobian, rhs);
}

//-----------------------------------------------------------------------------
// Function      : SharedSystem::debugOutput3
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void SharedSystem::debugOutput3 
   (Linear::Vector & dxVector, Linear::Vector & xVector)
{
  xyceInterfacePtr_->debugOutput3(dxVector, xVector);
}

}}}

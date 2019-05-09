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

//-------------------------------------------------------------------------
//
// Purpose        : Interface to let multiple N_NLS::NOX::Group's
//                  share a single system of RHS Vector, Jacobian
//                  matrix, Newton vector, and gradient vector.
//                  Closely related to the NOX::Epetra::SharedSystem
//                  class.
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

#ifndef Xyce_N_NLS_NOX_SharedSystem_h
#define Xyce_N_NLS_NOX_SharedSystem_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

// ----------   NOX Includes   ----------

#include "N_NLS_fwd.h"
#include "N_NLS_NOX_Vector.h"
#include "N_NLS_NOX_Group.h"

#include <vector>

class Ifpack_IlukGraph;
class Ifpack_CrsRiluk;

//-----------------------------------------------------------------------------
// Class         : N_NLS::NOX::SharedSystem
//
// Purpose       :
//
//      Interface to let multiple N_NLS::NOX::Group's share the
//      vectors and matrices in the Xyce nonlinear solver.
//
//      Closely related conceptually to the
//      NOX::Epetra::SharedJacobian class.
//
// Creator       : Tammy Kolda, SNL, 8950
//
// Creation Date : 1/31/02
//-----------------------------------------------------------------------------

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {
class SharedSystem {

public:

  SharedSystem(Xyce::Linear::Vector& x,
	       Xyce::Linear::Vector& f,
	       Xyce::Linear::Matrix& jacobian,
	       Xyce::Linear::Vector& newton,
	       Xyce::Linear::Vector& gradient,
	       Xyce::Linear::System& lasSys,
	       Interface& interface);

  ~SharedSystem();


  void reset(Xyce::Linear::Vector& x,
	     Xyce::Linear::Vector& f,
	     Xyce::Linear::Matrix& jacobian,
	     Xyce::Linear::Vector& newton,
	     Xyce::Linear::Vector& gradient,
	     Xyce::Linear::System& lasSys,
	     Interface& interface);

  //---------------------------------------------------------------------------
  // Function      : isJacobianOwner
  // Purpose       : Verify that the group pointed to by grp is owner of the
  //                 Jacobian matrix.
  //---------------------------------------------------------------------------
  inline bool isJacobianOwner(const Group* grp) const
  {
    return (grp == ownerOfJacobian_);
  };

  //---------------------------------------------------------------------------
  // Function      : areStateVectors
  // Purpose       : To compute a Jacobian, the state vectors must be 
  //                 updated with respect to the solution in the group.  
  //                 However, the state vectors are updated ONLY during 
  //                 calls to compute the residual.  This method checks 
  //                 to see if the state vectors still correspond to this 
  //                 group.  Returns true if state vectors are correct.
  //---------------------------------------------------------------------------
  inline bool areStateVectors(const Group* grp) const
  {
    return (grp == ownerOfStateVectors_);
  };

  bool computeF(const Vector& solution, Vector& F, const Group* grp);

  bool computeJacobian(Group* grp);

  bool computeNewton(const Vector& F, Vector& Newton,
		     Teuchos::ParameterList& params);

  bool computeGradient(const Vector& F, Vector& Gradient);

  bool applyJacobian(const Vector& input, Vector& result) const;

  bool applyJacobianTranspose(const Vector& input, Vector& result) const;

  bool computeDfDpMulti (const std::vector< int > & paramIDs, 
                         NOX::Abstract::MultiVector & dfdp, 
                         bool isValidF);

  Vector& getSolutionVector();

  // Take ownership of const Jacobian.
  const Xyce::Linear::Matrix& getJacobian() const;

  // Take ownership of Jacobian and get a reference to it.
  Xyce::Linear::Matrix& getJacobian(const Group* grp);

  // Take ownership of the state vectors.
  void getStateVectors(const Group* grp);

  // Get a pointer to the Xyce::Linear::System object
  Xyce::Linear::System* getLasSystem();

  // Use for debugging (corresponding to the ones in N_NLS_NonLinearSolver).
  void debugOutput1 (Xyce::Linear::Matrix & jacobian, Xyce::Linear::Vector & rhs);
  void debugOutput3 (Xyce::Linear::Vector & dxVector, Xyce::Linear::Vector & xVector);

  // Preconditioning objects for the Group::applyRightPreconditioning method
  bool computePreconditioner();
  bool deletePreconditioner();
  bool applyRightPreconditioning(bool useTranspose, 
				 Teuchos::ParameterList& params,
				 const Vector& input, 
				 Vector& result);

  // This is used to construct vectors in the group.  We used to
  // clone a time integrator vector but when the DC Op point fails,
  // somewhere (I have no idea where) it decides to delete vectors
  // before the Nonlinear solver can delete theirs.  This causes
  // a seg fault.
  // 
  Vector* cloneSolutionVector() const;

  // Take ownership of const newton vector
  const Vector & getNewtonVector() const;

  void printSoln(std::ostream &os) {xyceSolnPtr_->print(os);}
  void printRes(std::ostream &os) {xyceFPtr_->print(os);}

private:

  // Views of xyce objects used in the fill
  Vector* xyceSolnPtr_;                     // Solution vector
  Vector* xyceFPtr_;                        // Residual Vector
  Xyce::Linear::Matrix* xyceJacobianPtr_;   // Jacobian matrix
  Vector* xyceNewtonPtr_;                   // Newton Vector
  Vector* xyceGradientPtr_;                 // gradient Vector
  Xyce::Linear::System* xyceLasSysPtr_;     // LAS System
  Interface* xyceInterfacePtr_;  // Nonlinear Solver Interface

  // Flag for Matrix Free Loads tscoffe/tmei 07/29/08
  bool matrixFreeFlag_;

  const Group* ownerOfJacobian_;
  const Group* ownerOfStateVectors_;

  // Ifpack preconditioning objects for applyRightPreconditioning method
  mutable Ifpack_IlukGraph* ifpackGraphPtr_;
  mutable Ifpack_CrsRiluk* ifpackPreconditionerPtr_;

}; // class SharedSystem
}}} // namespace N_NLS_NOX

#endif // Xyce_N_NLS_NOX_SharedSystem_h


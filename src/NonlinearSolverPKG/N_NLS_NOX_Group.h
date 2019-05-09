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
// Purpose        : Interface to Xyce for NOX groups.
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

#ifndef Xyce_N_NLS_NOX_Group_h
#define Xyce_N_NLS_NOX_Group_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

// ----------   NOX Includes   ----------

#include "NOX_Abstract_Group.H"
#include "Teuchos_RCP.hpp"
#include "Xyce_config.h"

// ---------- Forward Declarations ----------
namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {
  class Vector;
  class SharedSystem;
}}}
namespace NOX {
  namespace Abstract {
    class Vector;
    class Group;
  }
  namespace Parameter {
    class List;
  }
}
class Ifpack_IlukGraph;
class Ifpack_CrsRiluk;

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

//-----------------------------------------------------------------------------
// Class         : N_NLS::NOX::Group
//
// Purpose       :
//
//      NOX Group Interface for Xyce
//
// Creator       : Tammy Kolda, SNL, 8950
//
// Creation Date : 2/1/02
//-----------------------------------------------------------------------------

class Group : public virtual NOX::Abstract::Group {

public:

  Group(SharedSystem& s);
  Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

  ~Group();

  NOX::Abstract::Group& operator=(const Group& source);
  NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);

  void setX(const Vector& input);
  void setX(const NOX::Abstract::Vector& input);

  void computeX(const Group& grp, const Vector& d, double step);
  void computeX(const NOX::Abstract::Group& grp, const NOX::Abstract::Vector& d, double step);

  NOX::Abstract::Group::ReturnType computeF();
  NOX::Abstract::Group::ReturnType computeJacobian();
  NOX::Abstract::Group::ReturnType computeGradient();
  NOX::Abstract::Group::ReturnType computeNewton(Teuchos::ParameterList& params);
  NOX::Abstract::Group::ReturnType applyJacobian(const Vector& input, Vector& result) const;
  NOX::Abstract::Group::ReturnType applyJacobian(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const;

  NOX::Abstract::Group::ReturnType applyJacobianTranspose(const Vector& input, Vector& result) const;
  NOX::Abstract::Group::ReturnType applyJacobianTranspose(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const;
  NOX::Abstract::Group::ReturnType 
    applyJacobianInverse(Teuchos::ParameterList& params, 
				const Vector& input, Vector& result) const;
  NOX::Abstract::Group::ReturnType 
    applyJacobianInverse(Teuchos::ParameterList& params, 
				const NOX::Abstract::Vector& input,
				NOX::Abstract::Vector& result) const;
  
  NOX::Abstract::Group::ReturnType 
    applyRightPreconditioning(bool useTranspose,
				     Teuchos::ParameterList& params,
				     const Vector& input, 
				     Vector& result) const;
  NOX::Abstract::Group::ReturnType 
    applyRightPreconditioning(bool useTranspose,
				     Teuchos::ParameterList& params,
				     const NOX::Abstract::Vector& input, 
				     NOX::Abstract::Vector& result) const;

  bool isF() const;
  bool isJacobian() const;
  bool isGradient() const;
  bool isNewton() const;
  bool linearSolverStatus () const;
  void setLinearSolverStatus (bool status) { linearStatus_ = status; }

  const NOX::Abstract::Vector& getX() const;
  const NOX::Abstract::Vector& getF() const;
  double getNormF() const;
  const NOX::Abstract::Vector& getGradient() const;
  const NOX::Abstract::Vector& getNewton() const;

  //! Return RCP to solution vector.
  Teuchos::RCP< const NOX::Abstract::Vector > getXPtr() const;
  
  //! Return RCP to F(x)
  Teuchos::RCP< const NOX::Abstract::Vector > getFPtr() const;
  
  //! Return RCP to gradient.
  Teuchos::RCP< const NOX::Abstract::Vector > getGradientPtr() const;
  
  //! Return RCP to Newton direction.
  Teuchos::RCP< const NOX::Abstract::Vector > getNewtonPtr() const;

  Teuchos::RCP<NOX::Abstract::Group> 
    clone(NOX::CopyType type = NOX::DeepCopy) const;

protected:

  // resets the isValid flags to false
  void resetIsValid_();

  // Throws an error
  void throwError(std::string method, std::string message) const;


protected:

  // Reference to the shared Newton system
  SharedSystem* sharedSystemPtr_;

  // NOX Vectors for storing values
  Teuchos::RCP<Vector> xVecPtr_;
  Vector& xVec_;
  Teuchos::RCP<Vector> fVecPtr_;
  Vector& fVec_;
  Teuchos::RCP<Vector> newtonVecPtr_;
  Teuchos::RCP<Vector> gradVecPtr_;

  // Booleans for tracking whether or not these values have been
  // computed for the currect x.
  bool isValidF_;
  bool isValidJacobian_;
  bool isValidGradient_;
  bool isValidNewton_;
  mutable bool isValidPreconditioner_;
  mutable bool linearStatus_;

  // Value of the 2-Norm of F
  double normF_;

  // Flag to determine if the solver should refactor the 
  // preconditioner (iterative) or Jacobian (direct).
  mutable bool haveSolverFactors_;

}; // class SharedSystem

}}} // namespace N_NLS_NOX

#endif // Xyce_N_NLS_NOX_SharedSystem_h


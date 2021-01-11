//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
// Purpose        : This is collection of non-member functions that help
//                  in the construction of block linear systems, like those
//                  found in AC or HB analysis.
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical Systems Modeling
//
// Creation Date  : 06/22/11
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef  Xyce_LAS_EPETRAHELPERS_H
#define  Xyce_LAS_EPETRAHELPERS_H

// ---------- Standard Includes ----------
#include <vector>
#include <string>
// ----------   Xyce Includes   ----------

#include <N_LAS_fwd.h>
#include <Teuchos_RCP.hpp>

#include <Epetra_CrsMatrix.h>
#include <Epetra_Operator.h>
#include <Epetra_LinearProblem.h>

// ---------- Forward Declarations ----------

namespace Xyce {
namespace Linear {

// Helper function for writing Epetra_LinearProblem to file
void writeToFile(const Epetra_LinearProblem& problem, std::string prefix, 
                 int file_number, bool write_map);


//-----------------------------------------------------------------------------
// Class         : EpetraMatrixAccess
// Purpose       : Class defining epetraObj() methods for matrix classes.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
class EpetraMatrixAccess
{
  public:

  // Empty constructor
  EpetraMatrixAccess() {}
  virtual ~EpetraMatrixAccess() {}

  virtual Epetra_CrsMatrix & epetraObj() = 0;
  virtual const Epetra_CrsMatrix & epetraObj() const = 0;
};

//-----------------------------------------------------------------------------
// Class         : EpetraTransOp
// Purpose       : Class for swapping an operator and its transpose.
//               : This is needed for performing adjoint solves using Belos.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 2/17/2015
//-----------------------------------------------------------------------------
class EpetraTransOp : public virtual Epetra_Operator
  {
  public:
    /// Basic constructor.
    ///
    /// \param Op [in/out] The operator to wrap.  EpetraTransOp's
    ///   Apply() method will invoke this operator's ApplyInverse()
    ///   method, and vice versa.
    EpetraTransOp (const Teuchos::RCP<Epetra_Operator>& Op);

    //! Virtual destructor, for memory safety of derived classes.
    virtual ~EpetraTransOp() {};

    //! Whether the operator knows how to apply its transpose.
    bool HasApplyTranspose() const {
      return true;
    }

    //! Apply the transpose of the operator to x, putting the result in y.
    int
    Apply (const Epetra_MultiVector &X,
           Epetra_MultiVector &Y) const;

    //! Apply the inverse of the transposed operator to x, putting the result in y.
    int
    ApplyInverse (const Epetra_MultiVector &X,
                  Epetra_MultiVector &Y) const;

    //! Return a human-readable string describing the operator.
    const char* Label() const {
      return "Epetra_Operator applying A^{T} as A";
    }

    //! Return the current UseTranspose setting.
    bool UseTranspose() const {
      return (!Epetra_Op->UseTranspose ());
    }

    int SetUseTranspose (bool UseTranspose_in) {
      return Epetra_Op->SetUseTranspose (!UseTranspose_in);
    }

    /// \brief Return true if this object can provide an approximate inf-norm.
    ///
    /// If this method returns false, then the \c NormInf() method
    /// should not be used.
    bool HasNormInf () const {
      return Epetra_Op->HasNormInf ();
    }

    /// \brief Return the infinity norm of the global matrix.
    ///
    /// The returned value only makes sense if HasNormInf() == true.
    double NormInf() const  {
      return Epetra_Op->NormInf ();
    }

    //! Return the Epetra_Comm communicator associated with this operator.
    const Epetra_Comm& Comm() const { return Epetra_Op->Comm(); };

    //! Return the Epetra_Map object representing the domain of this operator.
    const Epetra_Map& OperatorDomainMap() const {
      return Epetra_Op->OperatorDomainMap();
    }

    //! Return the Epetra_Map object representing the range of this operator.
    const Epetra_Map& OperatorRangeMap() const {
      return Epetra_Op->OperatorRangeMap();
    }

  private:
    //! The underlying operator that this EpetraTransOp instance wraps.
    Teuchos::RCP<Epetra_Operator> Epetra_Op;
};

} // namespace Linear
} // namespace Xyce

#endif

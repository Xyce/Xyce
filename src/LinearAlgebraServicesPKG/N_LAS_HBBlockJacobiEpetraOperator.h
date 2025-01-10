//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : This is an HB specific class that derives off of
// Epetra_Operator and supports the matrix-free block Jacobi preconditioner
// that we need for HB.  It takes a pointer to the HB loader so that it can 
// correctly do the apply function.
//
// Creator        : Heidi Thornquist, 1437
//
// Creation Date  : 11/12/08
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_LAS_HBBlockJacobiEpetraOperator_h
#define Xyce_N_LAS_HBBlockJacobiEpetraOperator_h

// ---------- Standard Includes ----------

#include <vector>

#include <Teuchos_RCP.hpp>
#include <Epetra_Operator.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_FilteredMatrix.h>
#include <N_LAS_fwd.h>
#include <N_LOA_HBLoader.h>

class Epetra_LinearProblem;
class Epetra_MultiVector;
class Epetra_Comm;
class Epetra_Map;
class Epetra_Import;
class Amesos_BaseSolver;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : HBBlockJacobiEpetraOperator
// Purpose       : Matrix Free Epetra Operator concrete class
// Special Notes : 
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------

class HBBlockJacobiEpetraOperator : virtual public Epetra_Operator
{

public:

  HBBlockJacobiEpetraOperator();

  virtual ~HBBlockJacobiEpetraOperator();

  void initialize(
      std::vector<Teuchos::RCP<Epetra_LinearProblem> >& epetraProblems,
      const std::vector<Teuchos::RCP<Amesos_BaseSolver> >& amesosSolvers,
      const std::vector<Teuchos::RCP<FilteredMatrix> >& diffCMatrix,
      const std::vector<Teuchos::RCP<FilteredMatrix> >& diffGMatrix,
      const Teuchos::RCP<Loader::HBLoader>& hbLoader,
      const Teuchos::RCP<HBBuilder>& hbBuilder,
      const std::vector<double>& freqs,
      const std::pair<int,int>& localRange,
      const bool hbOsc
      );

    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
	does not support transpose use, this method should return a value of -1.
      
    \Param[in]
	   UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

    \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does not support transpose.
    */
  int SetUseTranspose(bool UseTranspose);

    //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
    /*
    \param[in]
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param[out]
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
    */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  int Apply(const MultiVector& X, MultiVector& Y) const;

    //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    /*
    \param[in]
	   X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param[out]
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method must 
              support the case where X and Y are the same object.
    */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  int ApplyInverse(const MultiVector& X, MultiVector& Y) const;
  int ApplyBlockJacobi(const MultiVector& X, MultiVector& Y) const;
  int ApplyCorrection(const MultiVector& X, MultiVector& Y) const;

    //! Returns the infinity norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

       \warning This method must not be called unless HasNormInf() returns true.
    */ 
  double NormInf() const;

  //! Returns a character string describing the operator
  const char * Label() const;

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const;

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const;

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & Comm() const;

  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & OperatorDomainMap() const;

  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & OperatorRangeMap() const;

private:
  bool isInitialized_, isCorrected_, hbOsc_;
  int N_, M_;
  int numAugRows_;
  std::pair<int,int> myN_;
  std::vector<double> freqs_;
  std::vector<Teuchos::RCP<Epetra_LinearProblem> > epetraProblems_;
  std::vector<Teuchos::RCP<Amesos_BaseSolver> > amesosSolvers_;
  std::vector<Teuchos::RCP<FilteredMatrix> > diffCMatrix_, diffGMatrix_;
  Teuchos::RCP<Loader::HBLoader> hbLoader_;
  Teuchos::RCP<HBBuilder> hbBuilder_;

  // Serialized objects for transferring RHS and solution vectors.
  std::vector<Teuchos::RCP<Epetra_Map> > serialEpetraMap_;
  mutable std::vector<Teuchos::RCP<Epetra_Import> > serialImporter_;
  mutable Teuchos::RCP<Epetra_MultiVector> serialX_, serialY_;
};

// Non-member constructor
Teuchos::RCP<HBBlockJacobiEpetraOperator> blockJacobiOperator(
    std::vector<Teuchos::RCP<Epetra_LinearProblem> >& epetraProblems,
    const std::vector<Teuchos::RCP<Amesos_BaseSolver> >& amesosSolvers,
    const std::vector<Teuchos::RCP<FilteredMatrix> >& diffCMatrix,
    const std::vector<Teuchos::RCP<FilteredMatrix> >& diffGMatrix,
    const Teuchos::RCP<Loader::HBLoader>& hbLoader,
    const Teuchos::RCP<HBBuilder>& hbBuilder,
    const std::vector<double> & freqs,
    const std::pair<int,int>& localRange,
    const bool hbOsc
    );

} // namespace Linear
} // namespace Xyce

#endif // Xyce_N_LAS_HBBlockJacobiEpetraOperator_h


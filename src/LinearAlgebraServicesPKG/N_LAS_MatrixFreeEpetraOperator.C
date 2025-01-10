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

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Creator        : Todd Coffey, 1414
//
// Creation Date  : 09/04/08
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Operator.h>
#include <N_LAS_EpetraMultiVector.h>
#include <N_LAS_MatrixFreeEpetraOperator.h>

#include <N_PDS_Comm.h>
#include <N_PDS_EpetraHelpers.h>
#include <N_PDS_EpetraParMap.h>
#include <N_LAS_Vector.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : matrixFreeEpetraOperator
// Purpose       : non-member constructor
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
RCP<MatrixFreeEpetraOperator> matrixFreeEpetraOperator(
    RCP<Operator> linearOperator,
    RCP<const Parallel::ParMap> solutionMap
    )
{
  RCP<MatrixFreeEpetraOperator> epetraOperator =
    rcp(new MatrixFreeEpetraOperator);
  epetraOperator->initialize(linearOperator,
      solutionMap
      );
  return epetraOperator;
}


//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::MatrixFreeEpetraOperator
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
MatrixFreeEpetraOperator::MatrixFreeEpetraOperator()
{
  isInitialized_ = false;
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::MatrixFreeEpetraOperator
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
MatrixFreeEpetraOperator::~MatrixFreeEpetraOperator()
{
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::initialize
// Purpose       : Initialization
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
void MatrixFreeEpetraOperator::initialize(
      RCP<Operator> linearOperator,
      RCP<const Parallel::ParMap> solutionMap
    )
{
  linearOperatorRCPtr_ = linearOperator;
  solutionMap_ = Teuchos::rcp_dynamic_cast<const Parallel::EpetraParMap>(solutionMap);
  isInitialized_ = true;
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::SetUseTranspose
// Purpose       : Define if transpose Apply and ApplyInverse is to be used.
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
int MatrixFreeEpetraOperator::SetUseTranspose(bool UseTranspose)
{
  // This is not supported for the HB load layers.
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::Apply
// Purpose       : Apply matrix free operator with Epetra_MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
int MatrixFreeEpetraOperator::Apply(
  const Epetra_MultiVector& X,
  Epetra_MultiVector& Y
  ) const
{
  // Convert these to Linear::MultiVectors and call the other Apply

  // COPY the multi-vector data into new objects on the stack.
  Epetra_MultiVector* Xcopy = new Epetra_MultiVector(X); // This gets deleted by the Linear::MultiVector below
  Epetra_MultiVector* Ycopy = new Epetra_MultiVector(Y); // This gets deleted by the Linear::MultiVector below
  Linear::EpetraMultiVector las_X(Xcopy, true); // this co-ops the Epetra_MultiVector and uses (and owns) its memory
  Linear::EpetraMultiVector las_Y(Ycopy, true); // this co-ops the Epetra_MultiVector and uses (and owns) its memory
  int status = Apply(las_X,las_Y);
  // COPY the Ycopy data back into Y
  Y = las_Y.epetraObj();
  return(status);
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::Apply
// Purpose       : Apply matrix free operator with Linear::MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
int MatrixFreeEpetraOperator::Apply(
  const Linear::MultiVector& X,
  Linear::MultiVector& Y
  ) const
{
  if (!isInitialized_)
  {
    Report::DevelFatal0().in("MatrixFreeEpetraOperator::Apply")
      << "I'm not initialized!";
  }

  int status = linearOperatorRCPtr_->apply(X,Y);

  return status;
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::ApplyInverse
// Purpose       : Apply inverse of matrix free operator with Epetra_MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
int MatrixFreeEpetraOperator::ApplyInverse(
  const Epetra_MultiVector& X,
  Epetra_MultiVector& Y
  ) const
{
  Report::DevelFatal0().in("MatrixFreeEpetraOperator::ApplyInverse")
    << "is not supported!";
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::ApplyInverse
// Purpose       : Apply inverse of matrix free operator with Linear::MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
int MatrixFreeEpetraOperator::ApplyInverse(
  const Linear::MultiVector& X,
  Linear::MultiVector& Y
  ) const
{
  Report::DevelFatal0().in("MatrixFreeEpetraOperator::ApplyInverse")
    << "is not supported!";
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::NormInf
// Purpose       : Norm Inf of matrix
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
double MatrixFreeEpetraOperator::NormInf() const
{
  Report::DevelFatal0().in("MatrixFreeEpetraOperator::NormInf")
    << "is not supported!";
  return -1.0;
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::Label
// Purpose       : Label for operator
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
const char * MatrixFreeEpetraOperator::Label() const
{
  return "Matrix Free Harmonic Balance Epetra Operator";
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::UseTranspose
// Purpose       : Query for useTranspose setting
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
bool MatrixFreeEpetraOperator::UseTranspose() const
{
  // Use Transpose is not supported, so always return false.
  return false;
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::HasNormInf
// Purpose       : Query for normInf support
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
bool MatrixFreeEpetraOperator::HasNormInf() const
{
  // Norm Inf is not supported, so always return false.
  return false;
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::Comm
// Purpose       : Return Epetra_Comm object
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
const Epetra_Comm & MatrixFreeEpetraOperator::Comm() const
{
  if (!isInitialized_)
  {
    Report::DevelFatal0().in("MatrixFreeEpetraOperator::Comm")
      << "I'm not initialized!";
  }
  return(*Parallel::getEpetraComm(&solutionMap_->pdsComm()));
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::OperatorDomainMap
// Purpose       : Return Epetra_Map corresponding to domain of operator
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
const Epetra_Map & MatrixFreeEpetraOperator::OperatorDomainMap() const
{
  if (!isInitialized_)
  {
    Report::DevelFatal0().in("MatrixFreeEpetraOperator::OperatorDomainMap")
      << "I'm not initialized!";
  }
  return(*(solutionMap_->petraMap()));
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::OperatorRangeMap
// Purpose       : Return Epetra_Map corresponding to range of operator
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
const Epetra_Map & MatrixFreeEpetraOperator::OperatorRangeMap() const
{
  if (!isInitialized_)
  {
    Report::DevelFatal0().in("MatrixFreeEpetraOperator::OperatorRangeMap")
      << "I'm not initialized!";
  }
  return(*(solutionMap_->petraMap()));
}

} // namespace Linear
} // namespace Xyce

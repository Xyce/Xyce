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
// Purpose        : interface to preconditioner
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 09/27/07
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_Preconditioner_h
#define Xyce_N_LAS_Preconditioner_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_UTL_fwd.h>
#include <N_LAS_fwd.h>

// ---------- Trilinos Includes ----------

#include <Teuchos_RCP.hpp>
#include <Teuchos_Describable.hpp>

class Epetra_Operator;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : Preconditioner
// Purpose       : interface to preconditioner
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
class Preconditioner : public Teuchos::Describable
{

public:
  // Constructors
  Preconditioner() {}

  // Destructor
  virtual ~Preconditioner() {}

  // Set the preconditioner options
  virtual bool setOptions(const Util::OptionBlock & OB) = 0;
  virtual bool setDefaultOptions() = 0;

  // Set individual preconditioner options
  virtual bool setParam( const Util::Param & param ) = 0;

  // Set the matrix pattern for the preconditioner
  virtual bool initGraph( const Teuchos::RCP<Problem> & problem ) = 0;

  // Set the matrix values for the preconditioner
  virtual bool initValues( const Teuchos::RCP<Problem> & problem ) = 0;

  // Compute the preconditioner using the current matrix values.
  virtual bool compute() = 0;

  // Apply the preconditioner; y = M*x.
  virtual int apply( MultiVector & x, MultiVector & y ) = 0;

  // Return the preconditioner as an Epetra_Operator object.
  virtual Teuchos::RCP<Epetra_Operator> epetraObj() = 0;

private:

  // No copying
  Preconditioner(const Preconditioner & right);
  Preconditioner & operator=(const Preconditioner & right);

  // No comparison
  bool operator==(const Preconditioner & right) const;
  bool operator!=(const Preconditioner & right) const;

};

} // namespace Linear
} // namespace Xyce

#endif

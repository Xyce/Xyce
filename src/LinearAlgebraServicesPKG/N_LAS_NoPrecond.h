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

#ifndef Xyce_N_LAS_NoPrecond_h
#define Xyce_N_LAS_NoPrecond_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_fwd.h>
#include <N_LAS_Preconditioner.h>

// ---------- Trilinos Includes ----------

#include <Teuchos_RCP.hpp>
#include <Teuchos_Describable.hpp>

// ----------  Fwd Declares     ----------

class Epetra_Operator;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : NoPrecond
// Purpose       : interface to preconditioner
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
class NoPrecond : public Preconditioner
{

public:
  // Constructors
  NoPrecond() {}

  // Destructor
  virtual ~NoPrecond() {}

  // Set the preconditioner options
  virtual bool setOptions(const Util::OptionBlock & OB) { return true; }
  virtual bool setDefaultOptions() { return true; }

  // Set individual preconditioner options
  virtual bool setParam( const Util::Param & param ) { return true; }

  // Set the matrix pattern for the preconditioner
  virtual bool initGraph( const Teuchos::RCP<Problem> & problem ) { return true; }

  // Set the matrix values for the preconditioner
  virtual bool initValues( const Teuchos::RCP<Problem> & problem ) { return true; }

  // Compute the preconditioner using the current matrix values.
  virtual bool compute() { return true; }

  // Apply the preconditioner; y = M*x.
  virtual int apply( MultiVector & x, MultiVector & y ) { return 0; }

  // Return the preconditioner as an Epetra_Operator object.
  virtual Teuchos::RCP<Epetra_Operator> epetraObj() { return Teuchos::null; }

private:

  // No copying
  NoPrecond(const NoPrecond & right);
  NoPrecond & operator=(const NoPrecond & right);

  // No comparison
  bool operator==(const NoPrecond & right) const;
  bool operator!=(const NoPrecond & right) const;

};

} // namespace Linear
} // namespace Xyce

#endif // Xyce_N_LAS_NoPrecond_h

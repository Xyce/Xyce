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
// Purpose        : Implementation file for the Abstract interface to the
//                  vector types (RDP, RSP, CDP or CSP).
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Computational Sciences
//
// Creation Date  : 10/13/00
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_LAS_Vector.h>
#include <N_UTL_FeatureTest.h>

// ---------  Other Includes  -----------

#include <Epetra_Vector.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : Vector:::Vector
// Purpose       : constructor 
// Special Notes : 
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
Vector::Vector( Epetra_Vector * origV, bool isOwned )
: MultiVector( dynamic_cast<Epetra_MultiVector *>(origV), isOwned )
{
}

//-----------------------------------------------------------------------------
// Function      : Vector:::Vector
// Purpose       : constructor
// Special Notes : 
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
Vector::Vector( Epetra_Vector * overlapV, const Epetra_BlockMap& parMap, bool isOwned )
: MultiVector( dynamic_cast<Epetra_MultiVector *>(overlapV), parMap, isOwned )
{
}

//-----------------------------------------------------------------------------
// Function      : Vector::dotProduct
// Purpose       : Returns the dot product of "this" vector and another.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 02/08/17
//-----------------------------------------------------------------------------
double Vector::dotProduct( const Vector & y ) const
{
  double result = 0.0;
  int PetraError = aMultiVector_->Dot(*(y.aMultiVector_), &result);

  if (DEBUG_LINEAR)
    processError( "Vector::dotProduct - ", PetraError );

  return result;
}


} // namespace Linear
} // namespace Xyce

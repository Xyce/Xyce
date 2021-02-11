//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Function      : Vector::Vector
// Purpose       : constructor 
// Special Notes : 
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
Vector::Vector( Epetra_Vector * origV, bool isOwned )
: MultiVector( dynamic_cast<Epetra_MultiVector *>(origV), isOwned ),
  groundNode_(0.0)
{
}

//-----------------------------------------------------------------------------
// Function      : Vector::Vector
// Purpose       : constructor
// Special Notes : 
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
Vector::Vector( Epetra_Vector * overlapV, const Epetra_BlockMap& parMap, bool isOwned )
: MultiVector( dynamic_cast<Epetra_MultiVector *>(overlapV), parMap, isOwned ),
  groundNode_(0.0)
{
}

//-----------------------------------------------------------------------------
// Function      : clone
// Purpose       : vector clone function 
// Special Notes : clones shape, not values
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
Vector* Vector::cloneVector() const
{
  Vector* new_vec = 0;
  if ( pmap() )
  {
    if ( pmap() == omap() )
      new_vec = new Vector( *pmap() );
    else
      new_vec = new Vector( *pmap(), *omap() );
  }
  else
  {
    // We don't have a map, so perform a cloneCopy
    new_vec = new Vector( *this );
  }
  return new_vec;
}
  
//-----------------------------------------------------------------------------
// Function      : cloneCopy
// Purpose       : vector clone function 
// Special Notes : clones shape and values
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 04/09/03
//-----------------------------------------------------------------------------
Vector* Vector::cloneCopyVector() const
{
  return new Vector( *this );
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
  int PetraError = epetraObj().Dot(y.epetraObj(), &result);

  return result;
}


} // namespace Linear
} // namespace Xyce

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
// Purpose        : Specification file for the Abstract interface to the
//                  vector types (RDP, RSP, CDP or CSP).
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Computational Sciences
//
// Creation Date  : 10/13/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_Vector_h
#define Xyce_N_LAS_Vector_h

#include <string>
#include <map>

#include <N_LAS_MultiVector.h>

#include <Epetra_Vector.h>

class N_PDS_ParMap;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : Vector
// Purpose       : Provides an abstract interface for the vector type (RDP,
//                 RSP, CDP or CSP).
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Compuational Sciences
// Creation Date : 10/13/00
//-----------------------------------------------------------------------------
class Vector : public MultiVector
{

public:

  // Constructors to map to Petra constructors.
  Vector(N_PDS_ParMap & map)
  : MultiVector(map, 1)
  {}

  Vector( N_PDS_ParMap & map, N_PDS_ParMap & ol_map )
  : MultiVector( map, ol_map, 1 )
  {}

  //Copy constructor
  Vector( const Vector & right )
  : MultiVector(right)
  {}

  // Constructor that wraps an Epetra multivector inside a Linear::MultiVector.
  // This is used in the nonlinear solver and linear solver interface.
  Vector( Epetra_Vector * origV, bool isOwned);

  // Constructor takes the oMultiVector and generates the aMultiVector
  Vector( Epetra_Vector * overlapMV, const Epetra_BlockMap& parMap, bool isOwned = true );

  // Destructor
  virtual ~Vector() {}

  // Operation: operator []
  double & operator[] (int index)
  {
    if (index >= 0)
      return (*oMultiVector_)[0][index];
    else
      return groundNode_;
  }

  // Operation: operator []
  const double & operator[] (int index) const
  {
    if (index >= 0)
      return (*oMultiVector_)[0][index];
    else
      return groundNode_;
  }

  // Dot product with another vector.
  double dotProduct(const Vector & y) const;
};

} // namespace Linear
} // namespace Xyce

#endif


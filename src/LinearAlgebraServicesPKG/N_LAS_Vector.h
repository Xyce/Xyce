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

#include <N_LAS_MultiVector.h>
#include <N_PDS_fwd.h>

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

  // Default constructor
  Vector() {}

  // Destructor
  virtual ~Vector() {}

  virtual Vector& operator=( const Vector& right ) = 0;
  MultiVector & operator=(const MultiVector & right)
  { Vector::operator=( right ); return *this; }

  // Clone operation:
  virtual Vector* cloneVector() const = 0;

  // Clone operation:
  virtual Vector* cloneCopyVector() const = 0;

  // Operation: operator []
  virtual double & operator[] (int index) = 0;

  // Operation: operator []
  virtual const double & operator[] (int index) const = 0;

  // Dot product with another vector.
  virtual double dotProduct(const Vector & y) const = 0;

};

} // namespace Linear
} // namespace Xyce

#endif


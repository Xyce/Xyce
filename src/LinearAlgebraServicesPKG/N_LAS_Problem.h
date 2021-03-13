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

//-----------------------------------------------------------------------------
//
// Purpose        : interface to linear problem
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/17/04
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_Problem_h
#define Xyce_N_LAS_Problem_h

#include <N_LAS_fwd.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : Problem
// Purpose       : interface to linear problem
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/17/04
//-----------------------------------------------------------------------------
class Problem
{
public:

  //Default constructor.
  Problem()
  : A_(0), Op_(0), x_(0), b_(0), matrixFreeFlag_(false)
  {}

  //Constructors
  Problem( Operator* Op, MultiVector* x, MultiVector* b );
  Problem( Matrix* A, MultiVector* x, MultiVector* b );

  //Copy constructor
  Problem( const Problem& prob );

  //Destructor
  virtual ~Problem() {}

  // Access linear problem components
  Matrix* getMatrix () { return A_; }
  Operator* getOp () { return Op_; }
  MultiVector* getRHS() { return b_; }
  MultiVector* getLHS() { return x_; }

  bool matrixFree() const { return(matrixFreeFlag_); }

protected:

  Matrix* A_;
  Operator* Op_;
  MultiVector* x_;
  MultiVector* b_;

  bool matrixFreeFlag_;
};

} // namespace Linear
} // namespace Xyce

#endif

//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Purpose        : Linear algebra utilities.
//
// Special Notes  :
//
// Creator        : Heidi K. Thornquist, SNL, Computational Sciences
//
// Creation Date  : 8/9/12
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_Util_h
#define Xyce_N_LAS_Util_h

//-----------------------------------------------------------------------------
// Class         : N_LAS_Util
// Purpose       : Provides some basic utilities for linear algebra computations.
// Special Notes :
// Creator       : Heidi K. Thornquist, SNL, Computational Sciences
// Creation Date : 8/9/12
//-----------------------------------------------------------------------------
namespace Xyce {
namespace Linear {

  // y = alpha*A*x + beta*y
  void crsAxpy( int Nrows, double alpha, double* Aval, int* ArowPtr, int* AcolInd, double* xval, double beta, double* yval )
  {
    int nnz = ArowPtr[Nrows];

    // y *= beta
    for (int i=0; i<Nrows; i++)
      yval[i] *= beta;

    for (int i=0; i<Nrows; i++)
    {
      double sum=0.0;
      for (int j=ArowPtr[i]; j<ArowPtr[i+1]; j++)
      {
        sum += Aval[j]*xval[AcolInd[j]];
      }
      yval[i] += alpha*sum;
    }
  }

} // namespace Linear
} // namespace Xyce

#endif

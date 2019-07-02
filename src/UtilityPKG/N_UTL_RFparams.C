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

#include <N_UTL_RFparams.h>
#include <Teuchos_SerialDenseSolver.hpp>

namespace Xyce{
namespace Util {

//-----------------------------------------------------------------------------
// Function      : Util::ytos
// Purpose       : convert a Y Matrix into an S matrix with same Z0
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
void ytos(const Teuchos::SerialDenseMatrix<int, std::complex<double> >& y,
          Teuchos::SerialDenseMatrix<int, std::complex<double> >& s,
          const std::vector<double> & Z0sVec )
{
  Teuchos::SerialDenseMatrix<int, std::complex<double> > eye(y.numRows(), y.numCols());
  Teuchos::SerialDenseMatrix<int, std::complex<double> > GrefInv(y.numRows(), y.numCols());
  Teuchos::SerialDenseMatrix<int, std::complex<double> > Zref(y.numRows(), y.numCols());
  Teuchos::SerialDenseMatrix<int, std::complex<double> > Gref(y.numRows(), y.numCols());

  double z0;

  for (int r = 0; r < y.numRows(); r++)
  {
    z0 = Z0sVec[r] ;
    for (int c = 0; c < y.numCols(); c++)
    {
      eye(r, c) = r == c ? 1.0 : 0.0;
      Zref(r, c) = r == c ? z0 : 0.0;
      Gref(r, c) = r == c ? 1.0 / sqrt(z0) : 0.0;
      GrefInv(r, c) = r == c ? sqrt(z0) : 0.0;
    }
  }

  s.putScalar(0.0);

  Teuchos::SerialDenseMatrix<int, std::complex<double> > t1(y.numRows(), y.numCols());

  t1.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Zref, y, 1.0);
  Teuchos::SerialDenseMatrix<int, std::complex<double> > t2(eye);
  t2 += t1; // t2 = 1 + Zref * Ymat
  Teuchos::SerialDenseMatrix<int, std::complex<double> > t3(eye);
  t3 -= t1; // t3 = 1 - Zref * Ymat

  s.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Gref, t3, 1.0);
  // Sout = Gref * (1 - Zref * Ymat)
  // Invert t2
  Teuchos::SerialDenseSolver<int, std::complex<double> > ftSolver;
  ftSolver.setMatrix(Teuchos::rcp(&t2, false));
  ftSolver.invert();
  // t2 = (1 + Zref * Ymat)^-1
  t3.putScalar(0.0);
  t3.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, s, t2, 1.0);
  // t3 = Gref * (1 - Zref * Ymat)*(1 + Zref * Ymat)^-1
  s.putScalar(0.0);
  s.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, t3, GrefInv, 1.0);

  return;
}

//-----------------------------------------------------------------------------
// Function      : Util::ytoz
// Purpose       : convert a Y Matrix into a Z matrix
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/02/2019
//-----------------------------------------------------------------------------
void ytoz(const Teuchos::SerialDenseMatrix<int, std::complex<double> >& y,
          Teuchos::SerialDenseMatrix<int, std::complex<double> >& z)
{
  z = y;
  Teuchos::SerialDenseSolver<int, std::complex<double> > ftSolver;
  ftSolver.setMatrix(Teuchos::rcp(&z, false));
  ftSolver.invert();

  return;
}

//-----------------------------------------------------------------------------
// Function      : Util::ztoy
// Purpose       : convert a Z Matrix into a Y matrix
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/02/2019
//-----------------------------------------------------------------------------
void ztoy(const Teuchos::SerialDenseMatrix<int, std::complex<double> >& z,
          Teuchos::SerialDenseMatrix<int, std::complex<double> >& y)
{
  y = z;
  Teuchos::SerialDenseSolver<int, std::complex<double> > ftSolver;
  ftSolver.setMatrix(Teuchos::rcp(&y, false));
  ftSolver.invert();

  return;
}

} // namespace Util
} // namespace Xyce

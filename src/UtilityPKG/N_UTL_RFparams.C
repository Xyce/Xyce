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

#include <N_ERH_Message.h>
#include <N_UTL_RFparams.h>
#include <Teuchos_SerialDenseSolver.hpp>

namespace Xyce{
namespace Util {
//-----------------------------------------------------------------------------
// Function      : Util::ytos
// Purpose       : convert a Y Matrix into a S matrix
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
void ytos(const Teuchos::SerialDenseMatrix<int, std::complex<double> >& y,
          Teuchos::SerialDenseMatrix<int, std::complex<double> >& s,
          const std::vector<double> & Z0sVec )
{
  // Input y matrix must be square.  Z0sVec must have the correct length.
  if ( (y.numRows() != y.numCols()) || (y.numRows() != Z0sVec.size()) )
  {
    Report::DevelFatal().in("Util::ytos") << "Invalid dimensions or size for input Y matrix or Z0 vector";
  }

  Teuchos::SerialDenseMatrix<int, std::complex<double> > idenMat(y.numRows(), y.numCols());
  Teuchos::SerialDenseMatrix<int, std::complex<double> > ZrsMat(y.numRows(), y.numCols());

  double z0;

  for (int i = 0; i < y.numRows(); i++)
  {
    z0 = Z0sVec[i] ;

    for (int j = 0; j < y.numCols(); j++)
    {
      idenMat(i, j) = (i == j) ? 1.0 : 0.0;
      ZrsMat(i, j) = (i == j) ? sqrt(z0) : 0.0;
    }
  }

//  s.putScalar(0.0);

  Teuchos::SerialDenseMatrix<int, std::complex<double> > ZrsYZrs(y.numRows(), y.numCols()), tmpMat(y.numRows(), y.numCols()), tmpMat2(idenMat);

  tmpMat.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, ZrsMat, y, 0.0); 
  ZrsYZrs.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, tmpMat, ZrsMat, 0.0);

  tmpMat = idenMat;
  tmpMat -= ZrsYZrs;
  tmpMat2 += ZrsYZrs;

  Teuchos::SerialDenseSolver<int, std::complex<double> > denseSolver;
  denseSolver.setMatrix( Teuchos::rcp( &tmpMat2, false ) );
  denseSolver.invert(); 

  s.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, tmpMat, tmpMat2, 0.0);

  return;
}



//-----------------------------------------------------------------------------
// Function      : Util::stoy
// Purpose       : convert a S Matrix into a Y matrix
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
void stoy(const Teuchos::SerialDenseMatrix<int, std::complex<double> >& s,
          Teuchos::SerialDenseMatrix<int, std::complex<double> >& y,
          const std::vector<double> & Z0sVec )
{
  // Input s matrix must be square.  Z0sVec must have the correct length.
  if ( (s.numRows() != s.numCols()) || (s.numRows() != Z0sVec.size()) )
  {
    Report::DevelFatal().in("Util::stoy") << "Invalid dimensions or size for input S matrix or Z0 vector";
  }

  Teuchos::SerialDenseMatrix<int, std::complex<double> > idenMat(s.numRows(), s.numCols());
  Teuchos::SerialDenseMatrix<int, std::complex<double> > YrsMat(s.numRows(), s.numCols());

  double z0;

  for (int i = 0; i < s.numRows(); i++)
  {
    z0 = Z0sVec[i] ;

    for (int j = 0; j < s.numCols(); j++)
    {
      idenMat(i, j) = (i == j) ? 1.0 : 0.0;
      YrsMat(i, j) = (i == j) ? 1.0/sqrt(z0) : 0.0;
    }
  }

//  s.putScalar(0.0);

  Teuchos::SerialDenseMatrix<int, std::complex<double> > yMat(s.numRows(), s.numCols()), tmpMat(idenMat), tmpMat2(idenMat);

  tmpMat -= s;
  tmpMat2 += s;

  Teuchos::SerialDenseSolver<int, std::complex<double> > denseSolver;
  denseSolver.setMatrix( Teuchos::rcp( &tmpMat2, false ) );
  denseSolver.invert(); 

  yMat.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, tmpMat, tmpMat2, 0.0);

  tmpMat.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0,  YrsMat,  yMat, 0.0);

  y.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, tmpMat, YrsMat, 0.0);


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
  // Input y matrix must be square.
  if ( y.numRows() != y.numCols() )
  {
    Report::DevelFatal().in("Util::ytoz") << "Invalid dimensions for input Y matrix";
  }

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
  // Input z matrix must be square.
  if ( z.numRows() != z.numCols() )
  {
    Report::DevelFatal().in("Util::ztoy") << "Invalid dimensions for input Z matrix";
  }

  y = z;
  Teuchos::SerialDenseSolver<int, std::complex<double> > ftSolver;
  ftSolver.setMatrix(Teuchos::rcp(&y, false));
  ftSolver.invert();

  return;
}

} // namespace Util
} // namespace Xyce

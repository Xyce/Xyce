//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
// Purpose        : Utilities related to .LIN analyses (e.g., conversions
//                  between S-parameters, Y-parameters, Z-parmeters and
//                  H-parameters).
// Special Notes  :
// Creator        : Pete Sholander Sandia National Laboratories, 1355
// Creation Date  : 2019/07/01
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_RFparams_h
#define Xyce_N_UTL_RFparams_h

#include <Teuchos_SerialDenseMatrix.hpp>

namespace Xyce {
namespace Util {

void ytos(const Teuchos::SerialDenseMatrix<int, std::complex<double> > &y,
          Teuchos::SerialDenseMatrix<int, std::complex<double> > &s,
          const std::vector<double> & Z0sVec);

void stoy(const Teuchos::SerialDenseMatrix<int, std::complex<double> > &s,
          Teuchos::SerialDenseMatrix<int, std::complex<double> > &y,
          const std::vector<double> & Z0sVec);

void ytoz(const Teuchos::SerialDenseMatrix<int, std::complex<double> > &y,
          Teuchos::SerialDenseMatrix<int, std::complex<double> > &z);

void ztoy(const Teuchos::SerialDenseMatrix<int, std::complex<double> > &z,
          Teuchos::SerialDenseMatrix<int, std::complex<double> > &y);

} // namespace Util
} // namespace Xyce

#endif // Xyce_N_UTL_RFparams_h

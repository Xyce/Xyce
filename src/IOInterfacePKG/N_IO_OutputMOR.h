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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Baur
//
// Creation Date  :
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputMOR_h
#define Xyce_N_IO_OutputMOR_h

#include <complex>

#include <N_PDS_fwd.h>

#include <Teuchos_SerialDenseMatrix.hpp>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// MOR
class OutputMOR
{
public:
  OutputMOR(const std::string &netlist_filename);

  virtual ~OutputMOR();

private:
  OutputMOR(const OutputMOR &);
  OutputMOR &operator=(const OutputMOR &);

public:
  void output(Parallel::Machine comm, bool origSystem, double freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H);
  void reset();

private:
  std::string           netlistFilename_;
  std::ostream *        os_;
};

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputMOR_h

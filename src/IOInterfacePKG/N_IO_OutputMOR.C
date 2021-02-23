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
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <fstream>

#include <N_IO_OutputMOR.h>
#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : OutputMOR
// Purpose       : Output class for OutputMOR runs
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : OutputMOR::OutputMOR
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
OutputMOR::OutputMOR(const std::string &netlist_filename)
  : netlistFilename_(netlist_filename),
    os_(0)
{}

//-----------------------------------------------------------------------------
// Function      : OutputMOR::~OutputMOR
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
OutputMOR::~OutputMOR()
{
  delete os_;
}

//-----------------------------------------------------------------------------
// Function      : OutputMOR::outputMORTF
// Purpose       : .PRINT output for mor runs for original and reduced system
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, Ting Mei
// Creation Date : 5/25/12
//-----------------------------------------------------------------------------
void OutputMOR::output(Parallel::Machine comm, bool orig, double freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> > &H)
{
  if (Parallel::rank(comm) == 0)
  {
    if (!os_) //Setup Output Stream and Print Out Header
    {
      std::string filename = netlistFilename_ + (orig ? ".Orig" : ".Red") + ".FD.prn";

      os_ = new std::ofstream(filename.c_str());

      (*os_).setf(std::ios::scientific);
      (*os_).precision(16);
      (*os_).setf(std::ios::left, std::ios::adjustfield);

      (*os_) << std::setw(22) << "Frequency";
      for (int i = 0; i < H.numRows(); ++i)
        for (int j = 0; j < H.numRows(); ++j)
        {
          std::ostringstream os;
          os << "Re(H(" << i << ", " << j << "))";
          (*os_) << " " << std::setw(22) << os.str();

          os.str("");
          os << "Im(H(" << i << ", " << j << "))";
          (*os_) << " " << std::setw(22) << os.str();
        }
      (*os_) << std::endl;
    }

    (*os_) << freq;

    for (int i = 0; i < H.numRows(); ++i)
      for (int j = 0; j < H.numRows(); ++j)
        (*os_) << " " << std::setw(22) << H(i, j).real() << " " << std::setw(22) << H(i, j).imag();

    (*os_) << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMOR::doResetOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void OutputMOR::reset()
{
  delete os_;
  os_ = 0;
}

} // namespace IO
} // namespace Xyce

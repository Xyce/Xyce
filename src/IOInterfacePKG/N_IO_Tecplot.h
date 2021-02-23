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

#ifndef Xyce_N_IO_Tecplot_h
#define Xyce_N_IO_Tecplot_h

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : getTecplotTimeDateStamp
// Purpose       : Get current date and time and format for .PRINT output
// Special Notes : tecplot version of getTimeDateStamp.
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date : 6/14/2013
//-----------------------------------------------------------------------------
std::string getTecplotTimeDateStamp();

void tecplotTimeHeader(std::ostream &os, bool print_title, const std::string title, const Util::Op::OpList &op_list, const OutputMgr &output_manager);
void tecplotFreqHeader(std::ostream &os, bool print_title, const std::string title, const Util::Op::OpList &op_list, const OutputMgr &output_manager);

} // namespace Outputter
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_Tecplot_h

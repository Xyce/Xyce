//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Purpose        : Functions to output parameter definitions
//
// Special Notes  :
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_OutputPars_h
#define Xyce_N_DEV_OutputPars_h

#include <iosfwd>
#include <string>
#include <map>

#include <N_DEV_fwd.h>
#include <N_DEV_Pars.h>
#include <N_DEV_Print.h>

namespace Xyce {
namespace Device {

std::ostream &laTexDevice(std::ostream &os, const std::string &name, const int level, const int type, const std::string &device_description, const ParametricData<void> &parameters, const OutputMode::Mode format);

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_OutputPars_h

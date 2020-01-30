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

//-----------------------------------------------------------------------------
//
// Purpose        : Declares the DistributionTool class.  Distribution tool
//                  buffers and distributes circuit blocks (and related data
//                  such as option blocks, metadata, etc) for/during parsing.
//
// Special Notes  :
//
// Creator        : Eric Rankin, SNL
//
// Creation Date  : 03/12/2003
//
//
//
//
//-----------------------------------------------------------------------------


#ifndef Xyce_N_IO_DistributionTool_h
#define Xyce_N_IO_DistributionTool_h

#include <string>
#include <vector>
#include <map>
#include <list>

#include <N_IO_fwd.h>
#include <N_PDS_fwd.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class          : DistributionTool
// Purpose        : Buffers and distributes circuit blocks (and related data
//                  such as option blocks, metadata, etc) for/during parsing.
//-----------------------------------------------------------------------------
class DistributionTool
{
public:

  virtual ~DistributionTool() {}

  // send options, metatdata, and context to all procs
  virtual bool broadcastGlobalData() = 0;

  // Distribute devices using hierarchical context object.
  virtual void distributeDevices() = 0;

  // Return any additional option blocks after device distribution.
  virtual std::list<Util::OptionBlock>& getAdditionalOptions() = 0;

};

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_DistributionTool_h

//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Purpose        : This file contains the details of the dope info class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  :
//
//
//
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ----------  Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_DEV_DevicePDEInstance.h>
#include <N_UTL_Expression.h>
#include <N_DEV_SpecieSource.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<SpecieSource>::ParametricData()
{
  addPar("NAME", "none", &SpecieSource::name);
}

ParametricData<SpecieSource> &SpecieSource::getParametricData() {
  static ParametricData<SpecieSource> parMap;

  return parMap;
}

// ----------------------------------------------------------------------------
// Function      : SpecieSource::SpecieSource
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
// ----------------------------------------------------------------------------
SpecieSource::SpecieSource ()
  : CompositeParam (getParametricData()),
    name("V0")
{}

// ----------------------------------------------------------------------------
// Function      : SpecieSource::processParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date :
// ----------------------------------------------------------------------------
bool SpecieSource::processParam
(Param & ndParam, std::string & param, DevicePDEInstance & di)
{
  bool bsuccess = true;

  return bsuccess;
}

// ----------------------------------------------------------------------------
// Function      : SpecieSource::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date :
// ----------------------------------------------------------------------------
void SpecieSource::processParams()
{}

} // namespace Device
} // namespace Xyce

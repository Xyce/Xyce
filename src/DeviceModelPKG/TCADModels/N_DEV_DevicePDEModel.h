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

//-----------------------------------------------------------------------------
//
// Purpose        : This file contains the PDE device model base class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/15/01
//
//
//
//
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_DevicePDEModel_h
#define Xyce_N_DEV_DevicePDEModel_h

// ---------- Standard Includes ----------

// ------------- Xyce Includes ------------
#include <N_DEV_DeviceModel.h>

// ---------- Forward Declarations ----------
namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : DevicePDEModel
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------
class DevicePDEModel: public DeviceModel
{
public:
  DevicePDEModel(
     const ModelBlock &          model_block,
     ParametricData<void> &      parametric_data,
     const FactoryBlock &        factory_block)
    : DeviceModel(model_block, parametric_data, factory_block),
      dopeInfoMap(),
      regionMap(),
      nodeMap()
  {}

  virtual ~DevicePDEModel () {};

private:
  DevicePDEModel(const DevicePDEModel &);
  DevicePDEModel &operator=(const DevicePDEModel &);

public:
  std::map<std::string, DopeInfo *> dopeInfoMap;
  std::map<std::string, CompositeParam *> regionMap;
  std::map<std::string, CompositeParam *> nodeMap;
};

} // namespace Device
} // namespace Xyce

#endif


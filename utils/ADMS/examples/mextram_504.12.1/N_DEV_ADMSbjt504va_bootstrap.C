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
// Purpose        : 
// A device must be registered with the Xyce device subsystem in order
// for models and classes of the devie to be created.  This class
// implements a call to the device registration function on
// construction and the static object is created upon shareable object
// load.
//
// Special Notes  : 
//
// Creator        : Tom Russo
//
// Creation Date  : 1/11/2012
//
//-------------------------------------------------------------------------
#include <N_DEV_ADMSbjt504va.h>

struct Bootstrap 
{
  Bootstrap() 
  {
    Xyce::Device::ADMSbjt504va::registerDevice();
  }
};

Bootstrap s_bootstrap;




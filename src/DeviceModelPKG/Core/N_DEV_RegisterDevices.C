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
// Creator        : David Baur
//
// Creation Date  : 3/15/2013
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <cstring>

#ifdef HAVE_DLFCN_H
#include <dlfcn.h>
#endif

#include <N_DEV_Message.h>

#include <N_DEV_RegisterDevices.h>
#include <N_DEV_RegisterOpenDevices.h>
#include <N_DEV_RegisterIBISDevices.h>
#include <N_DEV_RegisterTCADDevices.h>

#ifdef Xyce_NEURON_MODELS
#include <N_DEV_RegisterNeuronDevices.h>
#endif

#include <N_DEV_RegisterExternalDevices.h>

#ifdef Xyce_ADMS_MODELS
#include <N_DEV_RegisterADMSDevices.h>
#endif

#ifdef Xyce_NONFREE_MODELS
#include <N_DEV_RegisterNonFreeDevices.h>
#endif

#ifdef Xyce_RAD_MODELS
#include <N_DEV_RegisterSandiaDevices.h>
#endif

#ifdef Xyce_MODSPEC_MODELS
#include <N_DEV_RegisterModSpecDevices.h>
#endif
namespace Xyce {
namespace Device {

void 
registerDevices(const DeviceCountMap& deviceMap,
                const std::set<int>& levelSet, bool includeMI)
{
  registerOpenDevices(deviceMap, levelSet, includeMI);

#ifdef Xyce_NEURON_MODELS
  registerNeuronDevices(deviceMap, levelSet);
#endif

#ifdef Xyce_ADMS_MODELS
  registerADMSDevices(deviceMap,levelSet);
#endif

  registerTCADDevices(deviceMap, levelSet);
  registerIBISDevices();
  registerExternalDevices();

#ifdef Xyce_RAD_MODELS
  registerSandiaDevices(deviceMap, levelSet);
#endif

#ifdef Xyce_NONFREE_MODELS
  registerNonFreeDevices(deviceMap, levelSet);
#endif
#ifdef Xyce_MODSPEC_MODELS
  registerModSpecDevices();
#endif
}

void registerDL(const char *so_path, const char *function_key = 0);

void
registerPlugin(const char *name) 
{
  registerDL(name);
}

typedef void (*dl_register_t)();

void registerDL(const char *so_path, const char *function_key) {
#ifdef HAVE_DLFCN_H
  void *dl = dlopen(so_path, RTLD_NOW);
  if (!dl) {
    const char *error = dlerror();
    Report::UserError0() << "Failed to load plugin " << error;
  }
  else {
    if (function_key) {
      std::string s = std::strlen(function_key) ? function_key : "dl_register";

      dl_register_t f = (dl_register_t) dlsym(dl, s.c_str());
      if (!f) {
        f = (dl_register_t) dlsym(dl, s.c_str());
      }

      if (f) {
        (*f)();
      }
      else {
        if (std::strlen(function_key)) {
          Report::UserError0() << "Executing dynamic library " << so_path << " function " << s << "()";
        }
      }
    }
  }
#endif // HAVE_DLFCN_H
}

} // namespace Device
} // namespace Xyce

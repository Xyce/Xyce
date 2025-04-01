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

//-----------------------------------------------------------------------------
// Filename       : $N_DAK_DakotaController.h$
//
// Purpose        : This class defines a controller that Xyce can use to
//                  start a Dakota based analysis.  The object here is to
//                  allow Xyce to read a netlist with some Dakota commands
//                  embeded within it and then pass off control to Dakota
//                  to organize and run the Xyce simulations it needs.
//
// Special Notes  :
//
// Creator        : Richard L. Schiek
//
// Creation Date  : 09/07/2005
//
//-----------------------------------------------------------------------------

#ifndef N_DAK_DakotaController_h
#define N_DAK_DakotaController_h

#include <string>
#include <vector>
#include <list>
#include <iosfwd>

#ifdef Xyce_Dakota // Dakota headers
#define DISABLE_DAKOTA_CONFIG_H   1 // Don't need dakota's dakota_config.h 
#ifdef Xyce_PARALLEL_MPI
#define DAKOTA_HAVE_MPI
#endif
#include <ParallelLibrary.hpp>
#include <ProblemDescDB.hpp>
#include <LibraryEnvironment.hpp>
#undef DISABLE_DAKOTA_CONFIG_H  
#undef DAKOTA_HAVE_MPI
#else
namespace Dakota {
class LibraryEnvironment
{};
}
#endif // Xyce_Dakota

#include <N_CIR_fwd.h>

namespace Xyce {
namespace Dakota {

class Controller
{
public:
    Controller(int argc, char *argv[]);
    ~Controller();

    // get the Dakota problem ready to run
    bool initializeDakotaEnvironment();
    bool constructApplicationInterface();
    bool executeDakotaStrategy();
    void retrieveDakotaResults();
    
    // start the Dakota controlled simulations
    bool run();

private:
    int iargsReduced;             // argument line without -dakota <filename>
    char ** cargsReduced;         //
    
    std::string dakotaInputFilename_;   // dakota input file 
    std::string outputFilename_;
    std::string writeRestartFilename_;

  ::Dakota::LibraryEnvironment *  libraryEnvironment_;          ///< dakota objects we'll need
};

} // namespace Dakota
} // namespace Xyce

#endif


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

//-----------------------------------------------------------------------------
// Filename       : Xyce.C
//
// Purpose        : front end for standalone Xyce executable
//
// Special Notes  :
//
// Creator        : Eric Rankin
//
// Creation Date  : 01/28/04
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


#include <string>

#ifdef HAVE_LINUX_EXCEPTIONS
#include <fenv.h>
#endif

#include <N_CIR_Xyce.h>
#include <N_DAK_DakotaController.h>
#include <N_ERH_ErrorMgr.h>

// Function to be called if memory runs out:
void _new_handler(void)
{
  Xyce::Report::UserFatal0() << "OUT OF MEMORY (error in 'new')";
}


//
//-----------------------------------------------------------------------------
// Function      : checkForDakotaFlag()
// Purpose       : This function scans the argument line to determine if this
//                 is a Dakota controlled run. (i.e. has -dakota on arg list)
// Special Notes :
// Scope         :
// Creator       : Rich Schiek
// Creation Date : 10/14/2008
//-----------------------------------------------------------------------------
inline
bool checkForDakotaFlag(const int argc, const char * const argv[])
{
  for (int i = 0; i < argc; i++)
  {
    std::string currentArg(argv[i]);

    if (currentArg == "-dakota")
      return true;
  }
  return false;
}

//
//-----------------------------------------------------------------------------
// Function      : main
// Purpose       : front end for standalone Xyce executable
// Special Notes :
// Scope         :
// Creator       : Eric Rankin
// Creation Date : 01/28/2004
//-----------------------------------------------------------------------------
int main( int argc, char *argv[] )
{
  // Set divide by zero, and invalid operation handling on linux
#ifdef HAVE_LINUX_EXCEPTIONS
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif

  // Set out of memory detection on all systems
  std::set_new_handler (&_new_handler);

  bool bsuccess = false;

  // this will be a Xyce simulation where Dakota creates and
  // then runs Xyce. So pass control to it
  if (checkForDakotaFlag(argc, argv))
  {
    Xyce::Dakota::Controller dakotaController( argc, argv );

    bsuccess = dakotaController.run();
    if ( ! bsuccess)
      Xyce::Circuit::Xyce_exit(1,false);
  }
  else
  {
    Xyce::Circuit::Simulator xyce;

    bsuccess = xyce.run(argc, argv);
    if ( ! bsuccess)
      Xyce::Circuit::Xyce_exit(1,false);
  }

  // This should only be reached when exiting after a successful run.
  // The code must never exit with a non-zero exit code in parallel in
  // any manner other than through MPI_Abort(), because it is undefined
  // behavior in the MPI standard, and in fact can trip a race condition
  // in some versions on OpenMPI on some platforms.
  return 0;
}

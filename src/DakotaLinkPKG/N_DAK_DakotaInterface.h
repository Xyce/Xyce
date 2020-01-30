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
// Filename       : $N_DAK_DakotaInterface.h$
//
// Purpose        : This class defines an interface between the optimization
//                  and analysis routines in Dakota and the circuit routines
//                  in Xyce.  Since this class derives from a Dakota class,
//                  Dakota can work through it to access Xyce.
//
// Special Notes  :
//
//
// Creator        : Richard L. Schiek
//
// Creation Date  : 09/07/2005
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DAK_DakotaInterface_h
#define Xyce_N_DAK_DakotaInterface_h 1

#include <string>
#include <vector>

#ifdef Xyce_Dakota // Dakota headers
#define DISABLE_DAKOTA_CONFIG_H   1 // Don't need dakota's dakota_config.h
#ifdef Xyce_PARALLEL_MPI
#define DAKOTA_HAVE_MPI
#endif
#include <DirectApplicInterface.hpp>
#include <ProblemDescDB.hpp>
#undef DISABLE_DAKOTA_CONFIG_H     // Don't need dakota's dakota_config.h
#undef DAKOTA_HAVE_MPI
#else
namespace Dakota {
class ProblemDescDB
{};

class DirectApplicInterface
{
public:
  DirectApplicInterface(const ProblemDescDB &)
  {}
};

class String
{};
}
#endif // Xyce_Dakota

#include <N_CIR_fwd.h>

// standard includes

namespace Xyce {
namespace Dakota {

class Interface: public ::Dakota::DirectApplicInterface
{
public:
  Interface( const ::Dakota::ProblemDescDB& problem_db );
  ~Interface();

  void setArguments( int iargsIn, char * cargsIn [] );

protected:
  int derived_map_ac( const ::Dakota::String& ac_name );

private:
  // these are the arguments that Xyce would have received from the command line
  int xyceIargs;
  char ** xyceCargs;

  // This is how we pass variable names and values into Xyce
  std::vector< std::pair< std::string, std::string > > variableSubVec;

  // Vector to hold the endpoint simulation values when the user has requested
  // that as the response function.
  std::vector< double > simulationDataValues;

  void copyCargs( const int originalIargs, char ** const originalCargs, char ** & copyCargs );
  void deleteCargs( const int len, char ** & cargs );
};

} // namespace Dakota
} // namespace Xyce

#endif


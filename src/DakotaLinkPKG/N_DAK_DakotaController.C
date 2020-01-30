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
// Filename       : $N_DAK_DakotaController.C$
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

#include <Xyce_config.h>

#include <iostream>
#include <string>
#include <cstring>

#ifdef Xyce_Dakota // Dakota headers
#define DISABLE_DAKOTA_CONFIG_H   1 // Don't need dakota's dakota_config.h
#ifdef Xyce_PARALLEL_MPI
#define DAKOTA_HAVE_MPI
#endif
#include <DakotaInterface.hpp>
#include <DakotaModel.hpp>
#include <DakotaResponse.hpp>
#include <DakotaVariables.hpp>
#undef DISABLE_DAKOTA_CONFIG_H
#undef DAKOTA_HAVE_MPI
#endif // Xyce_Dakota

#include <N_DAK_DakotaController.h>
#include <N_DAK_DakotaInterface.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_FeatureTest.h>

// global variable needed by Dakota for output.
// int Dakota::write_precision = 10;

namespace Xyce {
namespace Dakota {

Controller::Controller(
  int           argc,
  char *        argv[]):
  outputFilename_( "outputDakota.txt" ),
  writeRestartFilename_("restartOutDakota.txt")
{
  if (DEBUG_DAKOTA)
    Xyce::dout() << "In Controller::Controller()" << std::endl;

  // save the argc and argv passed in to the constructor.
  // however, we want to remove the -dakota <filename> part so that 
  // we can pass this to another N_CIR_Xyce object and not have that 
  // object think we need to run Dakota.  In the copying process,
  // save the dakota input file name.
  
  iargsReduced = argc - 2;  // we'll take all the args except "-dakota" and <filename>
  if( iargsReduced < 1 )
  {
    // remaingin arguments aren't enough to get at least a netlist.
    Report::UserError() << "Not enough arguments passed in for a dakota run (i.e. runxyce <netlist> -dakota <dakota input file>";
  }
  
  cargsReduced = new char * [iargsReduced];
  
  bool dakotaFileFound = false;
  for (int i=0,k=0; i<argc;++i)
  {
    // this seems odd, but it was how IO::CmdParse did this as well.
    if (argv[i] == NULL) 
    {
      cargsReduced[k] = NULL;
      continue;
    }

    std::string tmpString(argv[i]);
    if ( tmpString == "-dakota" )
    {
      // skip the copy of this 
      i++;
      // save the next item
      std::string filename(argv[i]);
      dakotaInputFilename_ = filename;
      dakotaFileFound = true;
    }
    else
    {
      //everything else gets copied
      int size = tmpString.size()+2;
      cargsReduced[k] = new char[size];
      for (int j=0;j<size;++j)
      {
        cargsReduced[k][j] = 0;
      }
      std::strcpy(cargsReduced[k], tmpString.c_str());
      k++;
    }
  }
  
  if(!dakotaFileFound)
  {
    // didn't find a dakota input file name so throw an error
    Report::UserError() << "Could not find the dakota input filename on the command line. (i.e. -dakota <filename>)";
  }
 
  // Initialize Dakota environment objects. 
  initializeDakotaEnvironment();

  // Construct application interface, passing in Xyce's interface objects to Dakota.
  constructApplicationInterface();
}

//
// destructor
//
Controller::~Controller()
{
  if (DEBUG_DAKOTA)
    Xyce::dout() << "In Controller::~Controller()" << std::endl;

  // need to delete the cargsReduced we created
  for (int i=0; i<iargsReduced; i++)
    delete [] cargsReduced[i];

  delete cargsReduced;
  delete libraryEnvironment_;
}


bool Controller::run()
{
  bool result = true;
  result = result && executeDakotaStrategy();
  retrieveDakotaResults();
  return true;
}

//
// set up Dakota components
//
bool Controller::initializeDakotaEnvironment()
{
  bool result = true;

#ifdef Xyce_Dakota
  if (DEBUG_DAKOTA)
    Xyce::dout() << "In Controller::initializeDakotaEnvironment()" << std::endl;

  ::Dakota::ProgramOptions opts;
  opts.input_file(dakotaInputFilename_);
  opts.output_file(outputFilename_);
  opts.write_restart_file(writeRestartFilename_);

  // Defaults constructs the MPIManager, which assumes COMM_WORLD
  libraryEnvironment_ = new ::Dakota::LibraryEnvironment(opts);

  if(!libraryEnvironment_ )
  {
    Report::DevelFatal().in("Controller::initializeDakotaEnvironment") << "Could not create Dakota::LibraryEnvironment";
  }

#ifdef Xyce_Dakota_Parallel_Debug
  numCommunicators = static_cast< int >( libraryEnvironment_->parallel_library().analysis_intra_communicators().size() );
  Xyce::dout() << "In Controller::initializeDakotaEnvironment() numCommunicators = " << numCommunicators << std::endl;
  for( int i=0; i<numCommunicators; i++ )
  {
    Xyce::dout() << "MPI_Comm is " << (libraryEnvironment_->parallel_library()).analysis_intra_communicators())[i] << std::endl;
  }
#endif

#endif // Xyce_Dakota
  return result;
}


bool Controller::constructApplicationInterface()
{
  bool result = true;

#ifdef Xyce_Dakota
  if (DEBUG_DAKOTA)
    Xyce::dout() << "In Controller::constructApplicationInterface()" << std::endl;

  Interface* aDakotaInterface = new Interface( libraryEnvironment_->problem_description_db() );

  // Set arguments, plug in interface, Dakota will manage object.
  aDakotaInterface->setArguments(iargsReduced,cargsReduced);
  libraryEnvironment_->plugin_interface( "", "direct", "", aDakotaInterface );

#ifdef Xyce_Dakota_Parallel_Debug
  numCommunicators = static_cast< int >( libraryEnvironment_->parallel_library().analysis_intra_communicators().size() );
  Xyce::dout() << "In Controller::initializeDakotaEnvironment() numCommunicators = " << numCommunicators << std::endl;
  for( int i=0; i<numCommunicators; i++ )
  {
    Xyce::dout() << "MPI_Comm is " << (libraryEnvironment_->parallel_library()).analysis_intra_communicators())[i] << std::endl;
  }
#endif   

#endif // Xyce_Dakota

return result;
}

bool Controller::executeDakotaStrategy()
{
#ifdef Xyce_Dakota  
  if (DEBUG_DAKOTA)
    Xyce::dout() << "In Controller::executeDakotaStrategy()" << std::endl;

  libraryEnvironment_->execute();

#endif // Xyce_Dakota

  return true;
}

void Controller::retrieveDakotaResults()
{
#ifdef Xyce_Dakota
  if (DEBUG_DAKOTA)
    Xyce::dout() << "In Controller::retrieveDakotaResults()" << std::endl;

  const ::Dakota::Variables & vars = libraryEnvironment_->variables_results();
  const ::Dakota::Response & resp = libraryEnvironment_->response_results();
#endif // Xyce_Dakota
}

} // namespace Dakota
} // namespace Xyce

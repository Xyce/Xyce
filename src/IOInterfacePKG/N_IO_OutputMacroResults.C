//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose        : Output Manager
//
// Special Notes  :
//
// Creator        : Dave Baur
//
// Creation Date  : 07/25/2014
//
//
//
//
//-------------------------------------------------------------------------

#include <sstream>
#include <fstream>
#include <Xyce_config.h>

#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>

#include <N_IO_MeasureManager.h>
#include <N_IO_FourierMgr.h>
#include <N_UTL_ExtendedString.h>

#include <Teuchos_oblackholestream.hpp>

namespace Xyce {
namespace IO {

typedef std::vector<std::pair<std::string, std::string> > StringPairVector;

//-----------------------------------------------------------------------------
// Function      : outputMacroResults
// Purpose       : if any post-process analysis was specified(like objective or measure)
//                 output results.  Called after simulation is over.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL Electrical and Microsystem Modeling
// Creation Date : 11/05/08
//-----------------------------------------------------------------------------
void outputMacroResults(
  Parallel::Machine             comm,
  const Measure::Manager &      measure_manager,
  FourierMgr &                  fourier_manager,
  std::string                   netlist_filename,
  const StringPairVector &      response_functions,
  std::string                   outputResponseFilename,
  const int                     step_number,
  const double                  endSimTime)
{
  /// The use of this black hole stream was replaced with calls to Xyce::lout()
  // This is a null output stream that helps ensure a function that needs to be called
  // on every processor in parallel only outputs on one processor.
  //Teuchos::oblackholestream outputBHS;
  std::ofstream outputFileStream;

  // Output the Measure results only if .measure is being performed on any variables.
  if (measure_manager.isMeasureActive())
  {
    if (Parallel::rank(comm) == 0) 
    {
      // this is a hack until this is function is refactored.
      // If this is a ".step" analysis then the latest measure
      // results have already been written to a file.  So the following code 
      // would just overwrite the file with the same results.  Thus, we check
      // if step_number==0.  If so then it wasn't a step analysis and we need to 
      // write out the results to a file (if we are on proc 0).
      if ( step_number == 0 )
      {
        // Do output to .mt (or .ms or .ma) file and stdout, based on MEASPRINT 
        // option from .OPTIONS MEASURE line.  The default is to both .mt files 
        // and stdout.
        if (measure_manager.isMeasGlobalPrintEnabled())
	{
          measure_manager.outputResultsToMTFile(step_number);
        }
        if (measure_manager.isMeasGlobalVerbosePrintEnabled())
	{
          measure_manager.outputVerboseResults( Xyce::lout(), endSimTime );
        }
      }
    }
    else 
    {
      // Other procs will write to the log file when running in parallel
      // Not sure if we really need to do the log output on other proces.
      // Also enable .MEASURE output to log files, based on MEASPRINT option.
      if (measure_manager.isMeasGlobalPrintEnabled())
      {  
        measure_manager.outputResults( Xyce::lout() );
      }
    }
  }

  // Output the Fourier results to file if Fourier analysis is being performed.
  // Make sure the function gets called on all processors, but only one outputs it.
  if (fourier_manager.isFourierActive())
  {
    if (Parallel::rank(comm) == 0)
    {
      std::string filename = netlist_filename + ".four";
      outputFileStream.open( filename.c_str() );
      fourier_manager.outputResults(outputFileStream);
      outputFileStream.close();
    }
    else {
      fourier_manager.outputResults( Xyce::lout() );
    }
  }

  // if the response list is not empty, try to dump those results to a file
  // a big limitation here is that all responses must be measure functions.
  // need to make this more flexible to include all solution vars and
  // objectives and results.
  if (!response_functions.empty())
  {
    std::ofstream responseOFS(outputResponseFilename.c_str());

    for (StringPairVector::const_iterator currRespItr = response_functions.begin(), endRespItr = response_functions.end();
         currRespItr != endRespItr; ++currRespItr)
    {
      double respvalue = -1.0;
      // need to parse name from XXX_X:name So find last ":" and extract
      // remaining string as value that dakota is looking for.
      std::string tempName = currRespItr->first;
      std::string::size_type beginningOfVarName = tempName.find_last_of(":");

      if (beginningOfVarName != std::string::npos)
      {
        int numChars = (currRespItr->first).length() - beginningOfVarName;
        tempName.assign(currRespItr->first, beginningOfVarName + 1, numChars);
      }
      //Xyce::dout() << " looking for " << currRespItr->first << std::endl;
      ExtendedString es(tempName);
      bool found = measure_manager.getMeasureValue(es.toUpper(), respvalue);
      responseOFS << respvalue   << "   " << tempName << std::endl;
    }
    responseOFS.close();
  }
}

} // namespace IO
} // namespace Xyce

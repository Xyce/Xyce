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

//-----------------------------------------------------------------------------
// Filename       : $N_DAK_DakotaInterface.C$
//
// Purpose        : This class defines an interface between the optimization
//                  and analysis routines in Dakota and the circuit routines
//                  in Xyce.  The object here is allow Xyce to read a netlist
//                  with some Dakota commands embeded within the netlist and
//                  pass off control to Dakota to organize and run as many
//                  Xyce simulations as are needed for Dakota's analysis.
//
// Special Notes  :
//
//
// Creator        : Richard L. Schiek
//
// Creation Date  : 09/07/2005
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <sstream>
#include <iostream>

#include <N_ANP_AnalysisManager.h>
#include <N_CIR_Xyce.h>
#include <N_DAK_DakotaInterface.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_OutputMgr.h>
#include <N_LAS_System.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>

namespace Xyce {
namespace Dakota {

Interface::Interface( const ::Dakota::ProblemDescDB& problem_db)
  : ::Dakota::DirectApplicInterface( problem_db )
{
  if (DEBUG_DAKOTA)
    Xyce::dout() << "Interface::Interface() " << std::endl;
}

Interface::~Interface()
{
  if (DEBUG_DAKOTA)
    Xyce::dout() << "Interface::~Interface() " << std::endl;

  deleteCargs( xyceIargs, xyceCargs );
}

//
//
//
void Interface::setArguments( int iargsIn, char * cargsIn [] )
{
  // copy the iargs and cargs
  xyceIargs = iargsIn;
  copyCargs( xyceIargs, cargsIn, xyceCargs );

  if (DEBUG_DAKOTA)
  {
    Xyce::dout() << "Interface::setArguments" << std::endl;
    for(int i=0; i<xyceIargs; i++)
      Xyce::dout() << "arg[" << i << "] = \"" << xyceCargs[i] << "\"" << std::endl;
  }
}

//
// routines called directly by Dakota
//

int Interface::derived_map_ac( const ::Dakota::String& ac_name)
{
#ifdef Xyce_Dakota
  if (DEBUG_DAKOTA)
  {
    Xyce::dout() << "Interface::derived_map_ac() ac_name = \"" << ac_name << "\"" << std::endl
                 << "numACV = " << numACV << " numADIV = " << numADIV << " numADRV = " << numADRV << std::endl
                 << "xCLabels.size() = " << xCLabels.size()
                 << "  xDILabels.size() = " << xDILabels.size()
                 << "  xDRLabels.size() = " << xDRLabels.size() << std::endl;
  }

  //
  // loop over all the continuous and discrete variables labels and values
  //
  // make sure the list of variable substitutions is clear
  variableSubVec.clear();

  // now fill the variable substitution list with fresh values
  // pack in the real values
  for (int i = 0; i < numACV; i++)
  {
    std::ostringstream oss;
    oss << std::setprecision(16) << xC[i];
    variableSubVec.push_back(std::make_pair(xCLabels[i], oss.str()));
  }

  // now get discrete integer variables
  for(int i=0; i < numADIV; i++)
  {
    std::ostringstream oss;
    oss << xDI[i];
    variableSubVec.push_back( std::make_pair(xDILabels[i], oss.str()));
  }

  // now get discrete real variables
  for(int i=0; i < numADRV; i++)
  {
    std::ostringstream oss;
    oss << xDR[i];
    variableSubVec.push_back(std::make_pair(xDRLabels[i], oss.str()));
  }

  if (DEBUG_DAKOTA)
  {
    Xyce::dout() << "Interface::derived_map_ac() variableSubVec" << std::endl;
    for (std::vector< std::pair<std::string,std::string> >::iterator it = variableSubVec.begin(), end = variableSubVec.end(); it != end; ++it)
      Xyce::dout() << "\"" << (*it).first << "\", \"" << (*it).second << "\"" << std::endl;
  }

  // create the Xyce object and pass it the list of variables and arguments
  Circuit::Simulator xyce;

  xyce.setNetlistParameters( variableSubVec );
  xyce.initialize( xyceIargs, xyceCargs );

  //
  // look at the response functions requested
  //

  if (DEBUG_DAKOTA)
    Xyce::dout() << "Interface::derived_map_ac() Requested response functions: "
                 << numFns << std::endl;

  // Check that the Response variables have been registered with .measure or .objective.
  //

  // Size the endpoint simulation array and check that these function labels are valid.
  simulationDataValues.resize( numFns );
  for (std::vector<std::string>::const_iterator it = fnLabels.begin(), end = fnLabels.end(); it != end; ++it)
  {
    bool success = xyce.checkResponseVar(*it);
    if (!success)
      Report::UserError() << "Response variable " << *it << " is not defined in the netlist";
  }

  // Create a suffix to distinguish this run's output from all others
  std::ostringstream runSuffix;
  runSuffix << ".run" << evaluation_id();
  xyce.setOutputFileSuffix(runSuffix.str());
  
  bool simResult = xyce.runSimulation();

  if (DEBUG_DAKOTA)
    Xyce::dout() << "Interface::derived_map_ac() Xyce iteration done." << std::endl;

  //
  // collect results from response functions
  //

  for (int i = 0; i < numFns; ++i)
  {
    bool found = xyce.obtainResponse(fnLabels[i], fnVals[i]);
    if (!found)
    {
      // if the simulation failed, override the result by returning a nan
      // so the optimizer isn't confused by a zero result when there wasn't one
      const char * taqp=NULL;
      fnVals[i] = nan( taqp );
    }
  }

  //
  // now clean up the Xyce object
  //
  xyce.finalize();

#endif // Xyce_Dakota
  
  return 0;
}

void Interface::copyCargs( const int originalIargs, char ** const originalCargs, char ** & copyCargs )
{
  copyCargs = new char * [originalIargs];

  for (int i=0; i<originalIargs;++i)
  {
    // this seems odd, but it was how IO::CmdParse did this as well.
    if (originalCargs[i] == NULL)
    {
      copyCargs[i] = NULL;
      continue;
    }

    std::string tmpString(originalCargs[i]);
    int size = tmpString.size()+2;
    copyCargs[i] = new char[size];

    for (int j=0;j<size;++j)
    {
      copyCargs[i][j] = 0;
    }
    strcpy(copyCargs[i], tmpString.c_str());
  }
}

void Interface::deleteCargs( const int len, char ** & cargs )
{
  for(int i=0; i<len; i++)
  {
    if( cargs[i] )
    {
      delete [] cargs[i];
      cargs[i]=NULL;
    }
  }
  if( cargs )
  {
    delete cargs;
    cargs=NULL;
  }
}

} // namespace Dakota
} // namespace Xyce

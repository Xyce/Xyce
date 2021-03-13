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
// Purpose       : This file contains functions to call Xyce as a library
//                 and test the functionality of the Xyce library interface.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 02/14/2008
//
//-----------------------------------------------------------------------------

//
// To make this a better test for specific bugs in the bug repository,
// I'll have the code check for specific functionality and emit an error
// only if that funciton fails.  The bugs checked are:
//
// BUG 1466  -- simulateUntil() fails to return upon ADC output change
// BUG 1499  -- Use raw name for ADC and DAC
// BUG 1500  -- Breakpoint needed at end of ramp up/down DACC
// BUG 1628  -- Initial ADC values needed from getTimeVoltagePairs()
//
// Not checked:
// BUG 1467  -- ADC / DAC are not compatible with newdae
//
// From here I can't check if Xyce is running in new dae mode
// however, the netlist used for this test ADC_DACtest.cir
// does not force the time integrator to use the old dae load.



// ---------- Standard Includes ----------
#include <iostream>
#include <ctime>

// ----------   Xyce Includes   ----------
#include <Xyce_config.h>
#include <N_CIR_Xyce.h>
#include <N_DEV_Device.h>
#include <N_DEV_Algorithm.h>
#include <N_DEV_ADC.h>


int main( int argc, char * argv[] )
{
  // flags for the bug's we're checking
  bool pass_bug_1466 = false;
  bool pass_bug_1499 = false;
  bool pass_bug_1500 = false;
  bool pass_bug_1628 = false;

  // this is just the amount of time ahead of the current time where
  // we should place the next breakpoint.  Making this larger makes
  // this test faster because the simulation is paused fewer times.
  const double deltaTimeForPause = 1.0e-5;

  // make a Xyce object
  N_CIR_Xyce xyce;

  // Initialize Xyce.
  if ( ! xyce.initialize(argc, argv) )
  {
    std::cerr << "Failed in N_CIR_Xyce::initialize" << std::endl;
    exit(-1);
  }

  std::vector< std::string > dacNames;
  if (1) {
    // Get the names of all DACs in the analog circuit.
    if ( ! xyce.getDACDeviceNames(dacNames) )
    {
      std::cerr << "Failed to get the DAC device names" << std::endl;
    }
  }

  if (1) {
    // now check that the names start with "Y" this idicates that we have
    // the raw device names.
    pass_bug_1499 = true;
    
    for(std::vector<std::string>::iterator nameIter = dacNames.begin(), end = dacNames.end(); nameIter != end; nameIter++ )
    {
      if (nameIter->find("Y") == std::string::npos )
      {
        pass_bug_1499 = false;
      }
    }
  }

  std::map<std::string, std::map<std::string,double> > ADCParamsMap_;
  if (1) {
    // Setting up ADC's
    // set up the ADC's
    // map of <ADCNAME, <parametername,parametervalue>>
    xyce.getADCMap(ADCParamsMap_);
    std::cout << "getACDMAP pass" << std::endl;
  }

  if (1) {
    // Setting up ADC width
    // must construct a map of ADC name/bit vector widths
    std::map<std::string,int> ADCWidthMap;
    std::map<std::string, std::map<std::string,double> > ::iterator currentADC = ADCParamsMap_.begin();
    std::map<std::string, std::map<std::string,double> > ::iterator endADC = ADCParamsMap_.end();
    const int defaultWidth = 8;
    while( currentADC != endADC )
    {
      ADCWidthMap[currentADC->first] = defaultWidth;
      currentADC++;
    }
    xyce.setADCWidths(ADCWidthMap);
  }

  std::map< std::string, std::vector< std::pair<double,double> > > xyceTimeVoltageUpdateMap;
  std::map< std::string, std::vector< std::pair<double,double> >* > timeVoltageUpdateMap;
  double dacValue = 1.0;
  if (0) {
    // Set up time voltage pairs
    xyce.getTimeVoltagePairs(xyceTimeVoltageUpdateMap);
    std::map< std::string, std::vector< std::pair<double,double> > >::iterator tvCurrentPair = xyceTimeVoltageUpdateMap.begin();
    std::map< std::string, std::vector< std::pair<double,double> > >::iterator tvEndPair = xyceTimeVoltageUpdateMap.end();

    //
    // for debugging print out the set of time voltage pairs
    //
    //   Xyce::dout() << "xyceTimeVoltageUpdateMap = ";
    //   for ( ; tvCurrentPair != tvEndPair; ++tvCurrentPair )
    //   {
    //     Xyce::dout() << "\"" << tvCurrentPair->first << "\": ";
    //     std::vector< std::pair<double,double> >::iterator currentPair = (tvCurrentPair->second).begin();
    //     std::vector< std::pair<double,double> >::iterator endPair = (tvCurrentPair->second).end();
    //     for( ;currentPair != endPair; ++currentPair)
    //     {
    //       Xyce::dout() << "( " << currentPair->first << ", " << currentPair->second << " ) ";
    //     }
    //   }
    //   Xyce::dout() << std::endl;

    // For each DAC being simulated...
    for(std::vector<std::string>::iterator nameIter = dacNames.begin(), end = dacNames.end(); nameIter != end; nameIter++ )
    {
      // update this dac
      std::vector< std::pair<double,double> >* timeVoltageUpdatesPtr = new std::vector< std::pair<double,double> >;
      timeVoltageUpdatesPtr->push_back( std::make_pair( deltaTimeForPause, dacValue ) );
      timeVoltageUpdateMap[*nameIter] = timeVoltageUpdatesPtr;
    }

    //
    // for debugging print out the new set of time voltage pairs
    //
    //   Xyce::dout() << "New time voltage pairs " << std::endl;
    //   std::map< std::string, std::vector< std::pair<double,double> >* >::iterator ntvCurrentPair = timeVoltageUpdateMap.begin();
    //   std::map< std::string, std::vector< std::pair<double,double> >* >::iterator ntvEndPair = timeVoltageUpdateMap.end();
    //
    //   Xyce::dout() << "timeVoltageUpdateMap = ";
    //   for ( ; ntvCurrentPair != ntvEndPair; ++ntvCurrentPair )
    //   {
    //     Xyce::dout() << "\"" << ntvCurrentPair->first << "\": ";
    //     std::vector< std::pair<double,double> >::iterator currentPair = ntvCurrentPair->second->begin();
    //     std::vector< std::pair<double,double> >::iterator endPair = ntvCurrentPair->second->end();
    //     for( ;currentPair != endPair; ++currentPair)
    //     {
    //       Xyce::dout() << "( " << currentPair->first << ", " << currentPair->second << " ) ";
    //     }
    //   }
    //   Xyce::dout() << std::endl;

    // Update the time-voltage pairs in Xyce
    if (timeVoltageUpdateMap.size() != 0)
    {
      //Xyce::dout() << "Calling updateTimeVoltagePairs " << std::endl;
      xyce.updateTimeVoltagePairs(timeVoltageUpdateMap);
    }
  }

  if (0) {
    // simulate the circuit
    bool simStatus=false;
    double actualTime = 0.0;

    double requstedTime = 1.0e-4;  // there is still a bug in Xyce if this is called with 0.0 and later resumed. RLS
    double nextDacPauseTime = deltaTimeForPause;
    do
    {

      //Xyce::dout() << "Calling simulateUntil( " << requstedTime << ", " << actualTime << ")" << std::endl;
      simStatus = xyce.simulateUntil(requstedTime, actualTime);
      //Xyce::dout() << "Simulation Paused at " << actualTime << std::endl;
      if (!simStatus )
      {
        Xyce::dout() << "Simulation exited early. " << std::endl;
        break;
      }

      if (!pass_bug_1628 )
      {
        // Need to get the time voltage pairs and see if we get the dc op value too
        xyceTimeVoltageUpdateMap.clear();
        xyce.getTimeVoltagePairs(xyceTimeVoltageUpdateMap);

        //Xyce::dout() << "BUG 1628 Check: xyceTimeVoltageUpdateMap = ";
        for (std::map< std::string, std::vector< std::pair<double,double> > >::iterator tvCurrentPair = xyceTimeVoltageUpdateMap.begin(), end = xyceTimeVoltageUpdateMap.end(); tvCurrentPair != end; ++tvCurrentPair )
        {
          //Xyce::dout() << "\"" << tvCurrentPair->first << "\": ";
          for(std::vector< std::pair<double,double> >::iterator currentPair = (tvCurrentPair->second).begin(), end = (tvCurrentPair->second).end(); currentPair != end; ++currentPair)
          {
            //Xyce::dout() << "( " << currentPair->first << ", " << currentPair->second << " ) ";
            if (currentPair->first == 0.0 )
            {
              // we have a dc op value so but 1628 is considered a pass
              pass_bug_1628 = true;
            }
          }
        }
      }

      if (actualTime > nextDacPauseTime )
      {
        //Xyce::dout() << "actualTime > nextDacPauseTime " << actualTime << " > " << nextDacPauseTime << std::endl;
        if (!pass_bug_1500 )
        {
          // On bug 1500 the ADC/DACs were not setting breakpoints that Xyce could use
          // to pause the simulation.  If we get this this part of the test run, then
          // the ADC/DACs are setting breakpoints, so this bug is considered passed
          pass_bug_1500 = true;
        }

        // try and add a new state change on the DAC
        if (dacValue > 0.0 )
        {
          dacValue = 0.0;
        }
        else
        {
          dacValue = 1.0;
          if (!pass_bug_1466 )
          {
            // On bug 1466 a change in voltage of a DAC was not setting a breakpiont
            // if we make it here then the default of DAC value of 1.0 set before this
            // loop has been set to zero which caused a break point to be set and brought
            // back here to puch the DAC back to 1.0.  Thus, the dac's are functioning
            // properly and bug 1466 is passing
            pass_bug_1466 = true;
          }

        }
        for(std::vector<std::string>::iterator nameIter = dacNames.begin(), end = dacNames.end(); nameIter != end; nameIter++ )
        {
          // update this dac
          std::vector< std::pair<double,double> >* timeVoltageUpdatesPtr = new std::vector< std::pair<double,double> >;
          timeVoltageUpdatesPtr->push_back( std::make_pair( actualTime+deltaTimeForPause, dacValue ) );
          timeVoltageUpdateMap[*nameIter] = timeVoltageUpdatesPtr;
        }

        // Update the time-voltage pairs.
        if (timeVoltageUpdateMap.size() != 0)
        {

          //Xyce::dout() << "Calling updateTimeVoltagePairs " << std::endl;
          xyce.updateTimeVoltagePairs(timeVoltageUpdateMap);
        }
        nextDacPauseTime += deltaTimeForPause;
      }

      requstedTime = 1.0e-4;
    } while ( actualTime < requstedTime );
  }

  if ( ! xyce.finalize() )
  {
    std::cerr << "Failed to N_CIR_Xyce::finalize\n";
    //    return;
    exit(1);
  }

  return 0;
}

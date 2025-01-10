//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Creation Date : 11/11/2021
//
//-----------------------------------------------------------------------------



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
  int retValue = 0;

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
  // Get the names of all DACs in the analog circuit.
  if ( ! xyce.getDACDeviceNames(dacNames) )
  {
    std::cerr << "Failed to get the DAC device names" << std::endl;
    retValue = -1;
  }
  

  std::map<std::string, std::map<std::string,double> > ADCParamsMap_;
  {
    // Setting up ADC's
    // set up the ADC's
    // map of <ADCNAME, <parametername,parametervalue>>
    if(xyce.getADCMap(ADCParamsMap_))
    {
      std::cout << "getACDMAP pass" << std::endl;
    }
    else
    {
      // getADCMap returns false when there are no ADC's in the circuit.
      // if that occurs in this case then it's a failure.
      std::cout << "getACDMAP FAIL" << std::endl;
      retValue = -1;
    }
  }
  
  // haven't done the DC op yet, so the max size will be 1.
  int tvpMaxSize = -1;
  if(xyce.getTimeVoltagePairsSz(tvpMaxSize))
  {  
    if( tvpMaxSize == 1)
    {
      std::cout << "getTimeVoltagePairsSz pass" << std::endl;
    }
    else
    {
      std::cout << "getTimeVoltagePairsSz FAIL " << tvpMaxSize << std::endl;
      retValue = -1;
    }
  }
  else
  {
    // getTimeVoltagePairsSz() return false if there are no ADC's in the circiut
    // that would be an error in this case so the test should indicate a fail.
    std::cout << "No ADC's in circuit os getTimeVoltagePairsSz FAIL " << tvpMaxSize << std::endl;
    retValue = -1;
  }
  
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
  if(xyce.setADCWidths(ADCWidthMap))
  {
    std::cout << "setADCWidths pass" << std::endl;
  }
  else
  {
    std::cout << "setADCWidths FAIL" << std::endl;
    retValue = -1;
  }
  
  bool simStatus=false;
  double actualTime = 0.0;
  double requstedTime = 1.0e-11; 
  simStatus = xyce.simulateUntil(requstedTime, actualTime);
  std::cout << "Did first call to simulate " << std::endl;
  std::cout << "Simulation Paused at " << actualTime << std::endl;
  std::map< std::string, std::vector< std::pair<double,double> > > xyceTimeVoltageUpdateMap;
  std::map< std::string, std::vector< std::pair<double,double> >* > timeVoltageUpdateMap;
  double dacValue = 1.0;
  
  // Set up time voltage pairs
  xyce.getTimeVoltagePairs(xyceTimeVoltageUpdateMap);

  // For each DAC being simulated...
  for(auto nameIter = dacNames.begin(); nameIter != dacNames.end(); nameIter++ )
  {
    // update this dac
    std::vector< std::pair<double,double> >* timeVoltageUpdatesPtr = new std::vector< std::pair<double,double> >;
    timeVoltageUpdatesPtr->push_back( std::make_pair( deltaTimeForPause, dacValue ) );
    timeVoltageUpdateMap[*nameIter] = timeVoltageUpdatesPtr;
  }

  // Update the time-voltage pairs in Xyce
  if (timeVoltageUpdateMap.size() != 0)
  {
    xyce.updateTimeVoltagePairs(timeVoltageUpdateMap);
  }

  // simulate the circuit
  requstedTime = 1.0e-5;  // there is still a bug in Xyce if this is called with 0.0 and later resumed.
  double finalSimTime = 1.0e-4;
  std::cout << "Entering loop with actualtime " << actualTime << " and requested time " << requstedTime << std::endl;
  double nextDacPauseTime = deltaTimeForPause;
  do
  {
    Xyce::dout() << "Calling simulateUntil( " << requstedTime << ", " << actualTime << ")" << std::endl;
    simStatus = xyce.simulateUntil(requstedTime, actualTime);
    Xyce::dout() << "Simulation Paused at " << actualTime << std::endl;
    if (!simStatus )
    {
      Xyce::dout() << "Simulation exited early. " << std::endl;
      break;
    }
    
    if(xyce.getTimeVoltagePairsSz(tvpMaxSize))
    { 
      // compare the size returned in the getTimeVoltagePairsSz() call to the actual size
      xyce.getTimeVoltagePairs(xyceTimeVoltageUpdateMap);    
      // during the transient simulation, the size of the time voltage pairs vector 
      // will depend on the complexity of the signal, the ADC width and how much simulation 
      // time has elapsed.  So compare it to the size of the actual map.
      int maxFound=0;
      for( auto tvCurrentPair=xyceTimeVoltageUpdateMap.begin(); tvCurrentPair!=xyceTimeVoltageUpdateMap.end(); tvCurrentPair++ )
      {
        if( tvCurrentPair->second.size() > maxFound )
        {
          maxFound = tvCurrentPair->second.size();
        }
      }
      if( tvpMaxSize == maxFound )
      {
        std::cout << "getTimeVoltagePairsSz pass" << std::endl;
      }
      else
      {
        std::cout << "getTimeVoltagePairsSz FAIL " << tvpMaxSize << " != " << maxFound << std::endl;
        retValue = -1;
      }
    }
    else
    {
      // getTimeVoltagePairsSz() return false if there are no ADC's in the circiut
      // that would be an error in this case so the test should indicate a fail.
      std::cout << "No ADC's in circuit os getTimeVoltagePairsSz FAIL " << tvpMaxSize << std::endl;
      retValue = -1;
    }
    if (actualTime > nextDacPauseTime )
    {
      // try and add a new state change on the DAC
      if (dacValue > 0.0 )
      {
        dacValue = 0.0;
      }
      else
      {
        dacValue = 1.0;
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
        xyce.updateTimeVoltagePairs(timeVoltageUpdateMap);
      }
      nextDacPauseTime += deltaTimeForPause;
    }

    requstedTime += 1.0e-5;
  } while ( actualTime < finalSimTime );


  if ( ! xyce.finalize() )
  {
    std::cerr << "Failed to N_CIR_Xyce::finalize" << std::endl;
    // return error code.
    exit(-1);
  }

  return retValue;
}

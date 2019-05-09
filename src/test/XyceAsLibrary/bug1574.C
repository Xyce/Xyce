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

//-----------------------------------------------------------------------------
// Purpose       : This file contains functions to call Xyce as a library
//                 and test the functionality of the Xyce library interface.
// Special Notes :
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/15/200901/15/2009
//
//-----------------------------------------------------------------------------

// ---------- Standard Includes ----------
#include <iostream>
#include <ctime>

// ----------   Xyce Includes   ----------
#include <N_CIR_Xyce.h>
#include <N_UTL_Math.h>

//-----------------------------------------------------------------------------
// Function      : main
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/02/09
//-----------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
  
  // make a Xyce object
  N_CIR_Xyce xyce;
  
  // Initialize Xyce.
  if ( ! xyce.initialize(argc, argv) )
  {
    std::cerr << "Failed in N_CIR_Xyce::initialize" << std::endl;
    //    return;
    exit(-1);
  }
 
#if 0 
  // Set up the ADC's
  // map of <ADCNAME, <parametername,parametervalue>>
  std::map<string, std::map<string,double> > ADCParamsMap_;
  xyce.getADCMap(ADCParamsMap_);
  
  // Setting up ADC width
  // must construct a map of ADC name/bit vector widths
  map<string,int> ADCWidthMap;
  std::map<string, std::map<string,double> > ::iterator currentADC = ADCParamsMap_.begin();
  std::map<string, std::map<string,double> > ::iterator endADC = ADCParamsMap_.end();
  const int defaultWidth = 8;
  while( currentADC != endADC )
  {
    ADCWidthMap[currentADC->first] = defaultWidth;
    currentADC++;
  }
  xyce.setADCWidths(ADCWidthMap);
#endif
  
  // Set up time voltage pairs
  map< string, vector< pair<double,double> > > xyceTimeVoltageUpdateMap;
  xyce.getTimeVoltagePairs(xyceTimeVoltageUpdateMap);
  map< string, vector< pair<double,double> > >::iterator tvCurrentPair = xyceTimeVoltageUpdateMap.begin();
  map< string, vector< pair<double,double> > >::iterator tvEndPair = xyceTimeVoltageUpdateMap.end();
  
  map< string, vector< pair<double,double> >* > timeVoltageUpdateMap;

  // Get the names of all DACs in the analog circuit.
  vector< string > dacNames;
  if ( ! xyce.getDACDeviceNames(dacNames) )
  {
    std::cerr << "Failed to get the DAC device names" << std::endl;
    //    return;
    exit(-1);
  }

  cout.width(14);cout.precision(8);cout.setf(ios::scientific);
  
  // simulate the circuit
  bool stepStatus=false;
  double oldTime = 0.0;
  double actualTime = 0.0;
  double finalTime = xyce.getFinalTime ();
    
    //0.25e-3;  // there is still a bug in Xyce if this is called with 0.0 and later resumed. RLS

  double initialTime = 0.0;
  double dt = (finalTime-initialTime)/100.0;
  double curr_computed_dt = dt;

  // For each DAC being simulated...
  vector<string>::iterator nameIter = dacNames.begin();
  vector<string>::iterator endIter = dacNames.end();

  // this is just the amount of time ahead of the current time where 
  // we should place the next breakpoint.  Making this larger makes 
  // this test faster because the simulation is paused fewer times.
  double dacValue = 1.0;
  double time=0.0;
  double xycetime=0.0;
  for( ; nameIter != endIter; nameIter++ )
  {
    vector< pair<double,double> >* timeVoltageUpdatesPtr(0);
    timeVoltageUpdatesPtr = new vector< pair<double,double> >;

    // update this dac
    for (int i=0;i<10;++i)
    {
      double scale = time/finalTime;
      dacValue = sin(M_PI*scale);

      timeVoltageUpdatesPtr->push_back( make_pair( time, dacValue ) );
      cout.width(14);cout.precision(8);cout.setf(ios::scientific);
      cout << "bug1574:134:" << *nameIter<< ": time,value="<<time<<","<<dacValue<<endl;

      time += 10*dt;
    }
    cout << "size of the vector: " << timeVoltageUpdatesPtr->size() << endl;
    timeVoltageUpdateMap[*nameIter] = timeVoltageUpdatesPtr;
  }
  cout << "double-checking the time,value pairs:" << endl;
  nameIter = dacNames.begin();
  for( ; nameIter != endIter; nameIter++ )
  {
    vector< pair<double,double> >* tVUPtr(0);

    tVUPtr = timeVoltageUpdateMap[*nameIter];
    for (int i=0;i<tVUPtr->size();++i)
    {
      double time = (*tVUPtr)[i].first;
      double value = (*tVUPtr)[i].second;
      cout << "bug1574:151:" << *nameIter << ": time,value="<<time<<","<<value<<endl;
    }
  }

  bool bs1 = true;
  cout.width(14);cout.precision(8);cout.setf(ios::scientific);
  cout << "Calling dcop provisionalStep.  actualTime = " << actualTime <<endl;
  cout << "curr_computed_dt = " << curr_computed_dt << endl;
  stepStatus = xyce.provisionalStep(0.0, curr_computed_dt, timeVoltageUpdateMap);

  //(double maxTimeStep, 
   //double &timeStep, 
   //map< string, vector< pair<double,double> > > & timeVoltageUpdateMap)

  xyce.acceptProvisionalStep();
  actualTime = xyce.getTime();

  double maxdt = 1.0e-7;

  do 
  {
    cout.width(14);cout.precision(8);cout.setf(ios::scientific);
    cout << "Calling dcop provisionalStep.  actualTime = " << actualTime <<endl;

    stepStatus = xyce.provisionalStep(maxdt, curr_computed_dt, timeVoltageUpdateMap);

    cout.width(14);cout.precision(8);cout.setf(ios::scientific);
    cout << "actual dt = "<<curr_computed_dt<<endl;

    if (curr_computed_dt > maxdt)
    {
      xyce.rejectProvisionalStep();
    }
    else
    {
      xyce.acceptProvisionalStep();
    }

    actualTime = xyce.getTime ();

  } while ( actualTime < finalTime );
    
  if ( ! xyce.finalize() )
  {
    cerr << "Failed to N_CIR_Xyce::finalize\n";
    exit(1);
  }
  
  exit(0);
}

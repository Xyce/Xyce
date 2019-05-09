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

#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <set>

using namespace std;

#include "Xyce_config.h"
#include "N_CIR_MixedSignalSimulator.h"

inline char *strdup(const char *s) 
{
  return strcpy((char *) malloc(strlen(s) + 1), s);
}

int main()
{
  string netlist("c7i.cir");
  Xyce::Circuit::MixedSignalSimulator *xycePtr;
  vector< string > dacNames;
  vector< string >::iterator d_i, d_end;
  map< string, vector< pair<double,double> >* > updates;
  map< string, vector< pair<double,double> > > dummy;
  double time, step_time, actual_step_time;
  int i;
  bool accept;

  cout << "Starting stand alone test of dcop in Xyce for mixed signal simulation" << endl;

  xycePtr = new Xyce::Circuit::MixedSignalSimulator();
  if (xycePtr)
  {
    cout << "Xyce instance successfully created" << endl;
  }

  int argc = 2;
  char* argv[3];
  argv[0] = strdup("Xyce");
  argv[1] = strdup(netlist.c_str());
  argv[2] = 0;

  if ( ! xycePtr->initialize(argc, argv) )
  {
    cout << "Failed to initialize Xyce for netlist " << netlist << endl;
    return 2;
  }
  free(argv[0]);
  free(argv[1]);

  // Get the names of all DACs in the analog circuit.
  dacNames.clear();
  if ( ! xycePtr->getDACDeviceNames(dacNames) )
  {
    cout << "Failed to get the DAC device names for netlist " << netlist << endl;
    return 3;
  }

  vector< pair<double,double> > inp_update;

  inp_update.push_back(pair<double,double>(0.0, 3.3));

  updates[string("XINPUT:YDAC!DA_INVD0")] = &inp_update;

  if (updates.size() != dacNames.size())
  {
    cout << "updates size = " << updates.size() << " and number of DACs is " << dacNames.size() << endl;
  }
  d_i = dacNames.begin();
  d_end = dacNames.end();
  for ( ; d_i != d_end ; ++d_i )
  {
    cout << "Checking DAC named: '" << *d_i << "'" << endl;
    if (updates.find(*d_i) == updates.end())
      cout << "DAC " << *d_i << " is not initialized" << endl;
    else
      cout << "DAC " << *d_i << " is initialized" << endl;
  }

  int i1, i2;
  double transitionTime = 0;
  double transitionVal = 3.3;

  xycePtr->updateTimeVoltagePairs (updates);
  cout << "updateTimeVoltagePairs called" << endl;
  step_time = 1e-6;
  actual_step_time = 0;
  bool stat;
  dummy.clear();
  time = xycePtr->getTime();
  stat = xycePtr->provisionalStep(step_time, actual_step_time, dummy);
  xycePtr->acceptProvisionalStep();
  time = 0;
  i = 0;
  while (time < 3e-9)
  {
    dummy.clear();
    cout << "Preparing to call provisionalStep at time = " << time << endl;
    stat = xycePtr->provisionalStep(step_time, actual_step_time, dummy);
    cout << "From provisional step, stat = " << stat << ", actual step time = " << actual_step_time << endl;
    if (stat != 1)
    {
      cout << "status from provisionalStep = fail!" << endl;
      break;
    }

    i1 = static_cast<int>(time/1e-9);
    i2 = static_cast<int>((time+actual_step_time)/1e-9);
    if (time > transitionTime && i1 != i2)
    {
      cout << "Rejected" << endl;
      transitionTime += 1e-9;
      transitionVal = 3.3 - transitionVal;
      inp_update.clear();
      inp_update.push_back(pair<double,double>(transitionTime, transitionVal));
      xycePtr->updateTimeVoltagePairs (updates);
      cout << "updateTimeVoltagePairs called for transition to " << transitionVal <<
              " at time = " << transitionTime << endl;
      xycePtr->rejectProvisionalStep();
    }
    else
    {
      xycePtr->acceptProvisionalStep();
      cout << "Accepted" << endl;
      time += actual_step_time;
    }
  }

  cout << "Simulation Complete" << endl;

  delete xycePtr;

  return 0;
}


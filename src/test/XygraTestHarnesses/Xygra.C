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


#include <N_CIR_Xygra.h>
#include <N_ERH_ErrorMgr.h>
#ifdef HAVE_LINUX_EXCEPTIONS
#include <fenv.h>
#endif

//#define MixedResistiveSourceCoil
#define SourceCoil
//#define ResistiveCoil

// Function to be called if memory runs out:
void _new_handler (void)
{
  Xyce::Report::UserFatal0() << "OUT OF MEMORY (error in 'new')";
}

//-----------------------------------------------------------------------------
// Function      : main
// Purpose       : front end for standalone Xyce executable
// Special Notes :
// Scope         :
// Creator       : Eric Rankin
// Creation Date : 01/28/2004
//-----------------------------------------------------------------------------
int main( int iargs, char *cargs[] )
{
  // Set divide by zero, and invalid operation handling on linux
#ifdef HAVE_LINUX_EXCEPTIONS
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif
  // Set out of memory detection on all systems
  set_new_handler (&_new_handler);

  N_CIR_Xygra xygra;

  bool bsuccess = xygra.initialize(iargs, cargs);
  vector <string> deviceNames;
  vector <double> vN;

  if (bsuccess)
    bsuccess = xygra.getDeviceNames("YXYGRA",deviceNames);

  if (bsuccess)
  {
    vector<vector<int> >coilWindings;
    vector<vector<string> >coilNames;
    coilWindings.resize(deviceNames.size());
    coilNames.resize(deviceNames.size());
    for (int i=0; i < deviceNames.size(); ++i)
    {
      xygra.xygraGetCoilWindings(deviceNames[i],coilWindings[i]);
      xygra.xygraGetCoilNames(deviceNames[i],coilNames[i]);
      cout << " Xygra device " << deviceNames[i] << " has "
           << coilWindings[i].size() << " coils " << endl;
      for (int j=0; j<coilWindings[i].size(); j++)
      {
        cout << "    coil " << j << " is named " << coilNames[i][j] << " and has " << coilWindings[i][j]
             << "windings" << endl;
      }
    }

//        bsuccess = xygra.runSimulation();
    double completedTime, timeStep;
    completedTime = 0.0;
    timeStep = 1e-2;
    bool opComplete = false;

    while (!(xygra.simulationComplete()) && bsuccess)
    {
      cout << "Simulation incomplete, completedTime is " << completedTime
           << "." << endl;

#ifdef MixedResistiveSourceCoil
      // This test will have 2 Xygra devices, the resistive coils one, and
      // one that's just a single coil with a single winding that will
      // just be a source (this will be our acceleration)
      for (int i = 0; i<deviceNames.size(); ++i)
      {
        if (deviceNames[i] == "YXYGRA!ACC")
        {
          if (completedTime == 0)
          {
            vector<double> sV(1,0.0);
            xygra.xygraSetSources(deviceNames[i],sV,completedTime);
          }
          // sign of current negative so that voltage at the positive
          // node winds up positive.
          double currentValue =
            -.001*sin(2*3.14159265358979*10.0*(completedTime+timeStep));
          vector<double> sV(1,currentValue);
          xygra.xygraSetSources(deviceNames[i],sV,completedTime+timeStep);
        }
        else
        {
          int numNodes = xygra.xygraGetNumNodes(deviceNames[i]);
          int numWindings = xygra.xygraGetNumWindings(deviceNames[i]);
          cout << " " << i << " " << deviceNames[i] << " has " << numNodes << " nodes and " << numWindings << "windings."<< endl;

          // We're now going to kludge this by faking every winding out as a 1K
          // resistor.
          double G = 1.0/1000.0;

          vector<vector<double> > kM;

          kM.resize(numWindings);
          if (completedTime == 0)
          {
            for (int winding=0; winding<numWindings; ++winding)
            {
              kM[winding].resize(numWindings,0.0);
              kM[winding][winding] = G;  // only set diagonal
            }
            xygra.xygraSetK(deviceNames[i],kM,completedTime);
          }
          G = 1/(1000.0+completedTime+timeStep);
          for (int winding=0; winding<numWindings; ++winding)
          {
            kM[winding].resize(numWindings,0.0);
            kM[winding][winding] = G;  // only set diagonal
          }
          xygra.xygraSetK(deviceNames[i],kM,completedTime+timeStep);

        }
      }
#endif

#ifdef ResistiveCoil
      for (int i = 0; i<deviceNames.size(); ++i)
      {
        int numNodes = xygra.xygraGetNumNodes(deviceNames[i]);
        int numWindings = xygra.xygraGetNumWindings(deviceNames[i]);
        cout << " " << i << " " << deviceNames[i] << " has " << numNodes << " nodes and " << numWindings << "windings."<< endl;

        // We're now going to kludge this by faking every winding out as a 1K
        // resistor.
        double G = 1.0/1000.0;

        vector<vector<double> > kM;

        kM.resize(numWindings);
        for (int winding=0; winding<numWindings; ++winding)
        {
          kM[winding].resize(numWindings,0.0);
          kM[winding][winding] = G;  // only set diagonal
        }
        xygra.xygraSetK(deviceNames[i],kM);
      }
#endif

#ifdef SourceCoil
      // We'll set each winding to a current source to simulate a rough
      // sinusoidal voltage with one period over the entire run
      for (int i = 0; i<deviceNames.size(); ++i)
      {
        cout << " setting sources for device " << deviceNames[i] << endl;
        int numWindings = xygra.xygraGetNumWindings(deviceNames[i]);

        if (completedTime == 0)
        {
          // set the t=0 version first
          vector<double> sV(numWindings,0.0);
          cout << " setting sources for t=0 " << endl;
          xygra.xygraSetSources(deviceNames[i],sV,completedTime);
        }

        double currentValue =  .001*sin(2*3.14159265358979*10.0*(completedTime+timeStep));
        vector<double> sV(numWindings,currentValue);
        cout << " setting sources for t= "<< completedTime+timeStep << endl;
        xygra.xygraSetSources(deviceNames[i],sV,completedTime+timeStep);
      }
#endif

      {
        cout << "Calling simulateUntil with requested time " << completedTime+timeStep << endl;

        bsuccess = xygra.simulateUntil(completedTime + timeStep, completedTime);

        cout << "Simulated to " << completedTime << endl;
        for (int i=0; i<deviceNames.size(); i++)
        {
          int offset=0;
          xygra.xygraGetVoltages(deviceNames[i], vN);
          cout << " Nodal voltages for device " << deviceNames[i] << endl;
          for (int coil=0; coil<coilWindings[i].size(); coil++)
          {
            cout << "   Coil " << coil << ":" << endl;
            for (int node=0; node<coilWindings[i][coil]+1;node++)
            {
              cout << "    node " << node << " voltage " << vN[offset++]
                   << endl;
            }
          }
        }
      }
    }
  }

  return bsuccess ? 0 : -1;
}


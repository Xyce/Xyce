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
#include <iostream>
#include <fstream>

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
  // Set out of memory detection on all systems
  set_new_handler (&_new_handler);
 
  N_CIR_Xygra * XycePtr = new N_CIR_Xygra();

  bool bsuccess = XycePtr->initialize(iargs, cargs);
  vector <string> deviceNames;
  vector <double> vN;

  if (bsuccess) 
    bsuccess = XycePtr->getDeviceNames("YXYGRA",deviceNames);

  ifstream infile("s-and-k.out",ifstream::in);

  if (bsuccess)
    bsuccess=infile.good();

  if (bsuccess)
  {
    vector<vector<int> >coilWindings;
    vector<vector<string> >coilNames;
    vector<int> numNodes;
    vector<int> numWindings;
    coilWindings.resize(deviceNames.size());
    coilNames.resize(deviceNames.size());
    numNodes.resize(deviceNames.size());
    numWindings.resize(deviceNames.size());
    for (int i=0; i < deviceNames.size(); ++i)
    {
      XycePtr->xygraGetCoilWindings(deviceNames[i],coilWindings[i]);
      XycePtr->xygraGetCoilNames(deviceNames[i],coilNames[i]);
      numNodes[i]=XycePtr->xygraGetNumNodes(deviceNames[i]);
      numWindings[i]=XycePtr->xygraGetNumWindings(deviceNames[i]);
      
      cout << " Xygra device " << deviceNames[i] << " has " 
           << coilWindings[i].size() << " coils " << endl;
      for (int j=0; j<coilWindings[i].size(); j++)
      {
        cout << "    coil " << j << " is named " << coilNames[i][j] << " and has " << coilWindings[i][j] 
             << "windings" << endl;
      }
      cout << "     for a total of " << numNodes[i] << " nodes and " << numWindings[i] << " windings." << endl;
    }

    if (deviceNames.size() != 1)
    {
      cerr << " Sorry, this test is designed to work only with one Xygra device. " << endl;
    }
    else
    {
      
      double completedTime, timeStep;
      completedTime = 0.0;
      bool opComplete = false;

      while (!(XycePtr->simulationComplete()) && bsuccess)
      {
        cout << "Simulation incomplete, completedTime is " << completedTime 
             <<"." << endl;


        // We will read each line of the input file sequentially to get
        // the target final time, the s vector and K matrix for that target
        // time.
        double targetTime;
        vector<double> sV;
        vector<vector<double> > kM;
        sV.resize(numWindings[0]); // we only have 1 device for sure
        kM.resize(numWindings[0]);

        infile >> targetTime;
        if (infile.eof()) break; // won't know this till after we try reading
/********
        for (int i=0; i<numWindings[0]; i++)
          infile >> sV[i];
*/
        for (int i=0; i<numWindings[0]; i++)
        {
          kM[i].resize(numWindings[0]);
          for (int j=0; j<numWindings[0]; j++)
            infile >> kM[i][j];
        }
        for (int i=0; i<numWindings[0]; i++)
          infile >> sV[i];

        // Dump the vector and matrices:

        cout << "S vector for time " << targetTime << ":" << endl;
        for (int i=0; i<numWindings[0]; i++)
          cout << sV[i] << " " ;
        cout << endl;

        cout << "K matrix for time " << targetTime << ":" << endl;
        cout << "--------" << endl;
        for (int i=0; i<numWindings[0]; i++)
        {
          for (int j=0; j<numWindings[0]; j++)
            cout << kM[i][j] << " " ;
          cout << endl; 
        }
        cout << "--------" << endl;

        // By leaving off time argument, we don't interpolate.
        XycePtr->xygraSetSources(deviceNames[0],sV);
        XycePtr->xygraSetK(deviceNames[0],kM);

        if (targetTime > 0.0)
        {
          cout << "Calling simulateUntil with requested time " << targetTime << endl;

          bsuccess = XycePtr->simulateUntil(targetTime,completedTime);
          cout << "Simulated to " << completedTime << endl;
          for (int i=0; i<deviceNames.size(); i++)
          {
            int offset=0;
            XycePtr->xygraGetVoltages(deviceNames[i], vN);
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
  }
  delete XycePtr;

  (bsuccess) ? exit(0) : exit(-1);
}


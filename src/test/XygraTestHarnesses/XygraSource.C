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

#include <Xyce_config.h>

#include <N_CIR_Xygra.h>
#include <N_UTL_fwd.h>
#include <N_ERH_ErrorMgr.h>

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
  std::set_new_handler (&_new_handler);

  N_CIR_Xygra * XycePtr = new N_CIR_Xygra();

  bool bsuccess = XycePtr->initialize(iargs, cargs);
  std::vector <std::string> deviceNames;
  std::vector <double> vN;

  if (bsuccess)
    bsuccess = XycePtr->getDeviceNames("YXYGRA", deviceNames);

  Xyce::lout() << "Starting Xyce as library test" << std::endl;

  if (bsuccess)
  {
    std::vector<std::vector<int> >coilWindings;
    std::vector<std::vector<std::string> >coilNames;
    coilWindings.resize(deviceNames.size());
    coilNames.resize(deviceNames.size());
    for (int i=0; i < deviceNames.size(); ++i)
    {
      XycePtr->xygraGetCoilWindings(deviceNames[i], coilWindings[i]);
      XycePtr->xygraGetCoilNames(deviceNames[i], coilNames[i]);
      Xyce::dout() << " Xygra device " << deviceNames[i] << " has "
           << coilWindings[i].size() << " coils " << std::endl;
      for (int j=0; j<coilWindings[i].size(); j++)
      {
        Xyce::dout() << "    coil " << j << " is named " << coilNames[i][j] << " and has " << coilWindings[i][j]
             << " windings" << std::endl;
      }
    }

    double completedTime, timeStep;
    completedTime = 0.0;
    timeStep = 1e-3;
    bool opComplete = false;

    while (!(XycePtr->simulationComplete()) && bsuccess)
    {

      // We'll set each winding to a current source to simulate a rough
      // sinusoidal voltage with one period over the entire run
      for (int i = 0; i<deviceNames.size(); ++i)
      {
        int numWindings = XycePtr->xygraGetNumWindings(deviceNames[i]);

        if (completedTime == 0)
        {
          // set the t=0 version first
          std::vector<double> sV(numWindings,0.0);
          XycePtr->xygraSetSources(deviceNames[i],sV,completedTime);
        }

        double currentValue =  .001*sin(2*3.14159265358979*10.0*(completedTime+timeStep));
        std::vector<double> sV(numWindings,currentValue);
        XycePtr->xygraSetSources(deviceNames[i],sV,completedTime+timeStep);
      }

      bsuccess = XycePtr->simulateUntil(completedTime+timeStep,completedTime);
    }

    if (bsuccess)
      bsuccess = XycePtr->finalize();
  }

  delete XycePtr;

  (bsuccess) ? exit(0) : exit(-1);
}


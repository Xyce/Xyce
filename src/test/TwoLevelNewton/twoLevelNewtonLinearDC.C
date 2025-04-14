//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Filename       : twoLevelNewton.C
//
// Purpose        : test program for two level Newton algorithm
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 12/5/2018
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <baseNewton.h>
#include <twoLevelNewtonLinear.h>

//-----------------------------------------------------------------------------
bool topLevelNewton::runDC( int iargs, char *cargs[])
{
  bool initJctFlag = true;

  double dV000 = 0.1;
  int numSteps = 101;

  Xyce::Device::ExternalSimulationData  extSimData;

  int linearSystemSize=5; 
  allocateLinearSystem(linearSystemSize);
  setupNonInnerMatrix();

  bool bsuccess = false;
  if (simulator_)
  {
    setupSimDataDCOP(extSimData);

    for (int i_DCOP_Step=0;i_DCOP_Step<numSteps;++i_DCOP_Step)
    {
      Xyce::lout() 
        << "DCOP step (outer solve) # " 
        << i_DCOP_Step 
        << "  vconnect000 = " 
        << vconnect000 
        << std::endl;

      if (i_DCOP_Step==0) { initJctFlag = true; }
      else                { initJctFlag = false; }

      resetSimDataDC(extSimData);
      bsuccess = simulator_->startTimeStep(extSimData);
      bsuccess = newtonSolve();

      if (bsuccess)
      {
        // process "successful" time step
        simulator_->stepSuccess(Xyce::Analysis::TWO_LEVEL_MODE_DC_SWEEP);
        simulator_->endTimeStep(extSimData);
      }
      else
      {
        Xyce::lout() << "DC step failed to converge. Exiting" << std::endl;
        exit(0);
      }

      vconnect000 += dV000;
    }

    bsuccess = simulator_->finishSolvers();
    bsuccess = simulator_->finalize();
  }

  (bsuccess) ? exit(0) : exit(-1);
}

//-----------------------------------------------------------------------------
int main( int iargs, char *cargs[] )
{
  Xyce::lout() << "Starting Xyce as library test" << std::endl;

  // Set divide by zero, and invalid operation handling on linux
  // Set out of memory detection on all systems
  std::set_new_handler (&_new_handler);

  std::string netlist = "resInner.cir";
  Xyce::lout() << "Setting netlist to: " << netlist << std::endl;
  topLevelNewton circuitCalculation(netlist);
  bool bsuccess = circuitCalculation.runDC(iargs,cargs);

  return 0;
}


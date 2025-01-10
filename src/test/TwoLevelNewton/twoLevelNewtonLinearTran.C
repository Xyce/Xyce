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
// Filename       : twoLevelNewtonLinearTran.C
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
//
// This is kind of the "main" function, which contains the outer time loop.
//
//-----------------------------------------------------------------------------
bool topLevelNewton::runTran( int iargs, char *cargs[])
{
  Xyce::Device::ExternalSimulationData  extSimData;

  // This "endTimeStep" call is really just to get values for a bunch
  // of things like "finalTime".
  simulator_->endTimeStep(extSimData);

  int linearSystemSize=5; 
  allocateLinearSystem(linearSystemSize);
  setupNonInnerMatrix();

  bool bsuccess = false;
  if (simulator_)
  {
    //---------------------------------------------------------------------------
    // TRANOP (initial steady-state) solve:
    {
      vconnect000 = 1.0;
      vconnect001 = 2.0;

      Xyce::lout() << "Begin DCOP step (outer solve)"  << std::endl;
      initJctFlag = true;
      setupSimDataTRANOP(extSimData);
      bsuccess = simulator_->startTimeStep(extSimData);
      bsuccess = newtonSolve();

      if (bsuccess)
      {
        // process "successful" time step
        simulator_->stepSuccess(Xyce::Analysis::TWO_LEVEL_MODE_TRANSIENT_DCOP);
        simulator_->endTimeStep(extSimData);
      }
      else
      {
        Xyce::lout() << "Time step failed to converge. Exiting" << std::endl;
        exit(0);
      }
    }

    // sine source vars
    double V0 = 1.0;
    double VA = 2.0;
    double FREQ = 10.0;
    double mpi = 2.0*std::asin(1.0);

    // set up a max time step
    double max_dt = extSimData.finalTime/500.0;

    //---------------------------------------------------------------------------
    // transient calculation. (post DCOP)
    extSimData.beginIntegrationFlag=true;
    int iTimeStep=0;
    bool finished = false;

    // time loop
    while(!finished)
    {
      // before resetting the sim data, check against max_dt.
      if ( extSimData.nextTimeStep > max_dt)
      {
        extSimData.nextTime -= extSimData.nextTimeStep;
        extSimData.nextTime += max_dt;
        extSimData.nextTimeStep = max_dt;
      }

      resetSimData(extSimData);

      Xyce::lout() 
        << "Begin Time step (outer solve) # " 
        << iTimeStep+1 
        << " current,step = " 
        << extSimData.nextTime
        << "," 
        << extSimData.nextTimeStep
        << std::endl;

      initJctFlag = false;

      // update the sinewave source:
      double time = (extSimData.nextTime);
      vconnect000 = V0 + VA * std::sin(2.0*mpi*(FREQ*time));

      bsuccess = simulator_->startTimeStep(extSimData);
      bsuccess = newtonSolve();

      if (bsuccess)
      {
        // process "successful" time step
        simulator_->stepSuccess(Xyce::Analysis::TWO_LEVEL_MODE_TRANSIENT); 
      }
      else
      {
        Xyce::lout() << "Time step failed to converge. Exiting" << std::endl;
        exit(0);
      }

      simulator_->endTimeStep(extSimData);

      if (extSimData.current_time >= extSimData.finalTime) { finished = true; }

      Xyce::lout() << extSimData;

      extSimData.beginIntegrationFlag=false;
      ++iTimeStep;
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

  std::string netlist = "resInnerTran.cir";
  Xyce::lout() << "---------------\nRunning the 2-level calculation:\n---------------" << netlist << std::endl;
  topLevelNewton circuitCalculation(netlist);
  bool bsuccess = circuitCalculation.runTran(iargs,cargs);

  return 0;
}


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
//
// Purpose        : This is a container class for solver information that is
//                  used during coupled simulation runs.  The outer code
//                  will populate this struct and pass it into a xyce
//                  inner solve using the simulateStep method on N_CIR_Xyce.
//
//
// Special Notes  : ERK, 2018:  As currently written, this class contains both 
//                  "inputs" which go from the top level solver (i.e. Charon) 
//                  down to the inner solver (i.e. Xyce), as well as "outputs"
//                  which go the other way.  
//
//                  The "inputs" get applied in the SecondLevelSimulator::startTimeStep function.
//                  The "outputs" get extracted in the SecondLevelSimulator::endTimeStep function.
//
//                  The input and output variables have (rightfully) been 
//                  kept separate, even though some of them refer to the same quantity.
//
//                  It might make sense to put them in separate structures later.
//
// Creator        : Roger P. Pawlowski, SNL, Applied Math and Applications
//
// Creation Date  : 07/16/2009
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_ExternalSimulationData_h
#define Xyce_N_DEV_ExternalSimulationData_h

#include <map>
#include <string>
#include <vector>

namespace Xyce {
namespace Device {

struct ExternalSimulationData
{
  // "inputs", set in the SecondLevelSimulator::startTimeStep function:
  bool is_transient;
  double current_time;
  double final_time;
  double current_time_step_size;
  double previous_time_step_size;
  int time_step_number;

  bool forceOrder;
  int imposedTimeIntegrationOrder;

  bool forceBeginningIntegration;
  bool imposedBeginningIntegration;

  // "outputs" set in the SecondLevelSimulator::endTimeStep function:
  int    currentOrder;
  int    numberOfSteps;
  int    usedOrder;
  int    nscsco;

  double pdt;
  double nextTimeStep;
  double currTimeStep;
  double currentTime;
  double nextTime;
  double finalTime;
  double startingTimeStep;
  double bpTol;
  bool   dcopFlag;
  bool   acopFlag;
  bool   inputOPFlag;
  bool   tranopFlag;
  bool   transientFlag;
  bool   dcsweepFlag;
  int    timeStepNumber;
  bool   initTranFlag;
  bool   beginIntegrationFlag;
  int    doubleDCOPStep;
  bool   doubleDCOPEnabled;
  int    timeIntMode;
  bool   sweepSourceResetFlag;
};


// inline functions
//-----------------------------------------------------------------------------
// Function      : SpecieSource::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
inline std::ostream & operator<<(std::ostream & os, const ExternalSimulationData & esd)
{
  // "inputs", set in the SecondLevelSimulator::startTimeStep function:
  os << "ExternalSimulationData.  Input variables:" <<std::endl;
  os << "is_transient = ";
  if (esd.is_transient) { os << "true" <<std::endl;} 
  else                  { os << "false" <<std::endl;} 

  os << "current_time = " << esd.current_time << std::endl;
  os << "final_time = " << esd.final_time << std::endl;
  os << "current_time_step_size = " << esd.current_time_step_size << std::endl;
  os << "previous_time_step_size = " << esd.previous_time_step_size << std::endl;
  os << "time_step_number = " << esd.time_step_number << std::endl;

  os << "forceOrder = ";
  if (esd.forceOrder) { os << "true" << std::endl; }
  else                { os << "false" << std::endl; }
  os << "imposedTimeIntegrationOrder = " << esd.imposedTimeIntegrationOrder << std::endl;

  os << "forceBeginningIntegration = ";
  if (esd.forceBeginningIntegration) { os << "true" << std::endl; }
  else                               { os << "false" << std::endl; }
  os << "imposedBeginningIntegration = " << esd.imposedBeginningIntegration << std::endl;

  // "outputs" set in the SecondLevelSimulator::endTimeStep function:
  os << "ExternalSimulationData.  Output variables:" <<std::endl;
  os << "currentOrder = " << esd.currentOrder << std::endl;
  os << "numberOfSteps = " << esd.numberOfSteps << std::endl;
  os << "usedOrder = " << esd.usedOrder << std::endl;
  os << "nscsco = " << esd.nscsco << std::endl;

  os << "pdt = " << esd.pdt << std::endl;
  os << "nextTimeStep = " << esd.nextTimeStep << std::endl;
  os << "currTimeStep = " << esd.currTimeStep << std::endl;
  os << "currentTime = " << esd.currentTime << std::endl;
  os << "nextTime = " << esd.nextTime << std::endl;
  os << "finalTime = " << esd.finalTime << std::endl;
  os << "startingTimeStep = " << esd.startingTimeStep << std::endl;
  os << "bpTol = " << esd.bpTol << std::endl;

  os << "timeStepNumber = " << esd.timeStepNumber << std::endl;
  os << "doubleDCOPStep = " << esd.doubleDCOPStep << std::endl;
  os << "timeIntMode = " << esd.timeIntMode << std::endl;


  os << "dcopFlag = ";
  if(esd.dcopFlag) { os << "true" << std::endl; }
  else         { os << "true" << std::endl; }

  os << "acopFlag = ";
  if(esd.acopFlag) { os << "true" << std::endl; }
  else         { os << "true" << std::endl; }

  os << "inputOPFlag = ";
  if(esd.inputOPFlag) { os << "true" << std::endl; }
  else         { os << "true" << std::endl; }

  os << "tranopFlag = ";
  if(esd.tranopFlag) { os << "true" << std::endl; }
  else         { os << "true" << std::endl; }

  os << "transientFlag = ";
  if(esd.transientFlag) { os << "true" << std::endl; }
  else         { os << "true" << std::endl; }

  os << "dcsweepFlag = ";
  if(esd.dcsweepFlag) { os << "true" << std::endl; }
  else         { os << "true" << std::endl; }

  os << "initTranFlag = ";
  if(esd.initTranFlag) { os << "true" << std::endl; }
  else         { os << "true" << std::endl; }

  os << "beginIntegrationFlag = ";
  if(esd.beginIntegrationFlag) { os << "true" << std::endl; }
  else         { os << "true" << std::endl; }

  os << "doubleDCOPEnabled = ";
  if(esd.doubleDCOPEnabled) { os << "true" << std::endl; }
  else         { os << "true" << std::endl; }

  os << "sweepSourceResetFlag = ";
  if(esd.sweepSourceResetFlag) { os << "true" << std::endl; }
  else         { os << "true" << std::endl; }

  os << std::endl;
  return os;
}

} // namespace Device
} // namespace Xyce

#endif


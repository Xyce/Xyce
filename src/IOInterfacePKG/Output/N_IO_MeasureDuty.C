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
// Purpose       : Transient analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_MeasureDuty.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : Duty::Duty()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
Duty::Duty(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  TranStats(measureMgr, measureBlock),
  totalOnTime_(0.0),
  totalAveragingWindow_(0.0),
  inOnState_(false)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();
}


//-----------------------------------------------------------------------------
// Function      : Duty::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void Duty::reset() 
{
  resetBase();
  totalAveragingWindow_=0.0;
  totalOnTime_=0.0;
}

//-----------------------------------------------------------------------------
// Function      : Duty::updateMeasureVars()
// Purpose       : Updates the onTime
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 05/4/2020
//-----------------------------------------------------------------------------
void Duty::updateMeasureVars(double circuitTime, double signalVal)
{
  if( ( (signalVal + minval_) > onValue_ ) || ( inOnState_ && ( (signalVal + minval_) > offValue_ )) )
  {
    // add to On duty time
    totalOnTime_ += (circuitTime - lastTimeValue_);
    inOnState_=true;
  }
  else
  {
    inOnState_=false;
  }

  totalAveragingWindow_ += (circuitTime - lastTimeValue_);

  return;
}

//-----------------------------------------------------------------------------
// Function      : Duty::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double Duty::getMeasureResult()
{
  if( initialized_ )
  {
    calculationResult_ =  totalOnTime_ / totalAveragingWindow_;
  }
  return calculationResult_;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

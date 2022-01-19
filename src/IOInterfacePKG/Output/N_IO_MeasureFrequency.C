//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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

#include <N_IO_MeasureFrequency.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : Frequency::Frequency()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
Frequency::Frequency(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  TranStats(measureMgr, measureBlock),
  offToOnCount_(0),
  onToOffCount_(0),
  totalAveragingWindow_(0.0)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();
}


//-----------------------------------------------------------------------------
// Function      : Frequency::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void Frequency::reset() 
{
  resetBase();
  totalAveragingWindow_=0.0;
  offToOnCount_=0.0;
  onToOffCount_=0.0;
}

//-----------------------------------------------------------------------------
// Function      : Frequency::updateMeasureVars()
// Purpose       : Updates the onTime
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 04/28/2020
//-----------------------------------------------------------------------------
void Frequency::updateMeasureVars(double circuitTime, double signalVal)
{
  totalAveragingWindow_ += (circuitTime - lastTimeValue_);

  // did we just start the "on" segment of a cycle?
  // if so count it and mark that we're in "on"
  if( (signalVal + minval_ ) > onValue_ )
  {
    // did we just cross into a new cycle?
    if( lastSignalValue_ < onValue_ )
    {
      offToOnCount_++;
    }
  }

  if( (signalVal - minval_) < offValue_ )
  {
    // did we just cross into a new cycle?
    if( lastSignalValue_ > offValue_ )
    {
      onToOffCount_++;
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : Frequency::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double Frequency::getMeasureResult()
{
  if( initialized_ )
  {
    calculationResult_ =  (0.5 *(onToOffCount_ + offToOnCount_))/totalAveragingWindow_;
  }
  return calculationResult_;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

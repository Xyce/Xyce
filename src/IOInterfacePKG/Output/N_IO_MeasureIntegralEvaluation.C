//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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

#include <N_IO_MeasureIntegralEvaluation.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : IntegralEvaluation::IntegralEvaluation()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
IntegralEvaluation::IntegralEvaluation(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Stats(measureMgr, measureBlock),
  integralValue_(0.0)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();
}


//-----------------------------------------------------------------------------
// Function      : IntegralEvaluation::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void IntegralEvaluation::reset() 
{
  resetStats();
  integralValue_ = 0.0;
}

//-----------------------------------------------------------------------------
// Function      : IntegralEvaluation::setMeasureVarsForNewWindow()
// Purpose       : Called when entering a new RFC window for TRAN measures
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 04/28/2020
//-----------------------------------------------------------------------------
void IntegralEvaluation::setMeasureVarsForNewWindow()
{
  integralValue_ = 0.0;
  initialized_ = false;

  return;
}

//-----------------------------------------------------------------------------
// Function      : IntegralEvaluation::updateMeasureVars()
// Purpose       : Updates the integral value
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 or NOISE measures, it is frequency.  For DC measures, it
//                 is the value of the first variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 04/28/2020
//-----------------------------------------------------------------------------
void IntegralEvaluation::updateMeasureVars(double indepVarVal, double signalVal)
{
  if (fromGiven_ && toGiven_)
  {
    if (mode_ == "DC")
    {
      // For DC mode, we must account for both ascending and descending FROM-TO windows and
      // the "direction" (increasing/decreasing) of first swept variable on .DC line.
      if ((from_ <= to_) && dcSweepAscending_)
      {
        integralValue_ += 0.5 * (indepVarVal - lastIndepVarValue_) * (signalVal + lastSignalValue_);
      }
      else if ((from_ >= to_) && !dcSweepAscending_)
      {
        integralValue_ += 0.5 * (indepVarVal - lastIndepVarValue_) * (signalVal + lastSignalValue_);
      }
      else
      {
        integralValue_ -= 0.5 * (indepVarVal - lastIndepVarValue_) * (signalVal + lastSignalValue_);
      }
    }
    else
    {
      // handles AC, NOISE and TRAN cases
      integralValue_ += 0.5 * (indepVarVal - lastIndepVarValue_) * (signalVal + lastSignalValue_);
    }
  }
  else
  {
    // For the other three cases (FROM given, or TO given, or neither given, the
    // (indepVarVal - lastIndepVarValue_) term accounts for the sweep (integration) direction
    integralValue_ += 0.5 * (indepVarVal - lastIndepVarValue_) * (signalVal + lastSignalValue_);
  }

  return;
}


//-----------------------------------------------------------------------------
// Function      : IntegralEvaluation::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double IntegralEvaluation::getMeasureResult()
{
  if( initialized_ )
  {
    calculationResult_ =  integralValue_;
  }
  return calculationResult_;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

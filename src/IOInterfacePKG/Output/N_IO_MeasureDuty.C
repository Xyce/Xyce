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
  Base(measureMgr, measureBlock),
  totalOnTime_(0.0),
  totalAveragingWindow_(0.0),
  lastTimeValue_(0.0),
  inOnState_(false)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();
}

//-----------------------------------------------------------------------------
// Function      : Duty::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void Duty::prepareOutputVariables() 
{
  // this measurement should have only one dependent variable.
  // Error out if it doesn't
  numOutVars_ = outputVars_.size();
  if ( numOutVars_ > 1 )
  {
    std::string msg = "Too many dependent variables for DUTY measure, \"" + name_ + "\"";
    Report::UserError0() << msg;
  }

  outVarValues_.resize( numOutVars_, 0.0 );
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
// Function      : Duty::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void Duty::updateTran(
  Parallel::Machine comm,
  const double circuitTime,
  const Linear::Vector *solnVec,
  const Linear::Vector *stateVec,
  const Linear::Vector *storeVec,
  const Linear::Vector *lead_current_vector,
  const Linear::Vector *junction_voltage_vector,
  const Linear::Vector *lead_current_dqdt_vector)
{
  if( !calculationDone_ && withinTimeWindow( circuitTime ) )
  {
    // update our outVarValues_ vector
    updateOutputVars(comm, outVarValues_, circuitTime,
      solnVec, stateVec, storeVec, 0, lead_current_vector,
      junction_voltage_vector, lead_current_dqdt_vector, 0);

    // Need to set lastOutputValue_ variable to the current signal value
    // at the first time-step within the measurement window. (That
    // window is set by the TO-FROM and TD qualifiers if present.)  This is 
    // needed so that the RISE/FALL/CROSS count is not incremented at time=0, if
    // the measured waveform has a DC offset at time=0      
    if (!firstStepInMeasureWindow_)
    {
      lastOutputValue_ = outVarValues_[0]; 
      firstStepInMeasureWindow_ = true;
    }

    // rfcLevel_ has a default to 0.0 if RFC_LEVEL qualifier is not specified
    if( withinRiseFallCrossWindow( outVarValues_[0], rfcLevel_ ) )
    {
      // If LAST was specified then this is done
      // each time a new RFC window is entered.
      if( newRiseFallCrossWindowforLast() )
      {
        totalOnTime_ = 0.0;
        totalAveragingWindow_ = 0.0;
        initialized_ = false;
        firstStepInRfcWindow_ = false;
      }

      // record the start and end times of the RFC window
      if( !firstStepInRfcWindow_  )
      {
        firstStepInRfcWindow_ = true;
        rfcWindowFound_ = true;
        rfcWindowStartTime_ = circuitTime;
      }
      rfcWindowEndTime_ = circuitTime;

      if( initialized_  )
      {
        if( ( (outVarValues_[0] + minval_) > onValue_ ) || ( inOnState_ && ( (outVarValues_[0] + minval_) > offValue_ )) )
        {
          // add to On duty time
          totalOnTime_ += (circuitTime - lastTimeValue_);
          inOnState_=true;
        }
        else
        {
          inOnState_=false;
        }

        // Calculating an average of outVarValues_[0];
        totalAveragingWindow_ += (circuitTime - lastTimeValue_);
      }

      lastTimeValue_ = circuitTime;
      initialized_=true;
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : Duty::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void Duty::updateDC(
  Parallel::Machine comm,
  const std::vector<Analysis::SweepParam> & dcParamsVec,
  const Linear::Vector *solnVec,
  const Linear::Vector *stateVec,
  const Linear::Vector *storeVec,
  const Linear::Vector *lead_current_vector,
  const Linear::Vector *junction_voltage_vector,
  const Linear::Vector *lead_current_dqdt_vector)
{

}

//-----------------------------------------------------------------------------
// Function      : Duty::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical Models & Simulation
// Creation Date : 8/7/2019
//-----------------------------------------------------------------------------
void Duty::updateAC(
  Parallel::Machine comm,
  const double frequency,
  const Linear::Vector *solnVec,
  const Linear::Vector *imaginaryVec,
  const Util::Op::RFparamsData *RFparams)
{

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

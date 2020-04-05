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
  Base(measureMgr, measureBlock),
  integralValue_(0.0),
  lastIndepVarValue_(0.0),
  lastSignalValue_(0.0)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();
}

//-----------------------------------------------------------------------------
// Function      : IntegralEvaluation::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void IntegralEvaluation::prepareOutputVariables()
{
  outVarValues_.resize( outputVars_.size(), 0.0 );
  
  // this measurement should have only one dependent variable.
  if (outVarValues_.size() > 1 )
    Xyce::Report::UserError0() << "Too many dependent variables for INTEG measure, \"" << name_ << "\"";
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
  resetBase();
  integralValue_ = 0.0;
}

//-----------------------------------------------------------------------------
// Function      : IntegralEvaluation::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void IntegralEvaluation::updateTran(
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
      junction_voltage_vector, lead_current_dqdt_vector, 0 );

    // Need to set lastOutputValue_ variable to the current signal value
    // at the first time-step within the measurement window.  (That
    // window is set by the TO-FROM and TD qualifiers if present.)  This is 
    // needed so that the RISE/FALL/CROSS count is not incremented at time=0, if
    // the measured waveform has a DC offset at time=0    
    if (!firstStepInMeasureWindow_)
    {
      lastOutputValue_ = outVarValues_[0]; 
      firstStepInMeasureWindow_ = true;
    }

    if( withinRiseFallCrossWindow( outVarValues_[0], rfcLevel_ ) )
    {
      // If LAST was specified then this is done
      // each time a new RFC window is entered.
      if( newRiseFallCrossWindowforLast() )
      {
        integralValue_ = 0.0;
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
        // update integral of outVarValues_[0];
        integralValue_ += 0.5 * (circuitTime - lastIndepVarValue_) * (outVarValues_[0] + lastSignalValue_);
      }

      lastIndepVarValue_ = circuitTime;
      lastSignalValue_ = outVarValues_[0];
      initialized_=true;
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : IntegralEvaluation::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 3/25/2020
//-----------------------------------------------------------------------------
void IntegralEvaluation::updateDC(
  Parallel::Machine comm,
  const std::vector<Analysis::SweepParam> & dcParamsVec,
  const Linear::Vector *solnVec,
  const Linear::Vector *stateVec,
  const Linear::Vector *storeVec,
  const Linear::Vector *lead_current_vector,
  const Linear::Vector *junction_voltage_vector,
  const Linear::Vector *lead_current_dqdt_vector)
{
  // The dcParamsVec will be empty if the netlist has a .OP statement without a .DC statement.
  // In that case, a DC MEASURE will be reported as FAILED.
  if ( dcParamsVec.size() > 0 )
  {
    double dcSweepVal = dcParamsVec[0].currentVal;

    // Used in descriptive output to stdout. Store name and first/last values of
    // first variable found in the DC sweep vector
    sweepVar_= dcParamsVec[0].name;
    if (!firstSweepValueFound_)
    {
        startSweepValue_ = dcSweepVal;
        firstSweepValueFound_ = true;
    }
    endSweepValue_ = dcSweepVal;

    if( !calculationDone_ && withinDCsweepFromToWindow( dcSweepVal ) )
    {
      outVarValues_[0] = getOutputValue(comm, outputVars_[0],
                                        solnVec, stateVec, storeVec, 0,
                                        lead_current_vector,
                                        junction_voltage_vector,
                                        lead_current_dqdt_vector, 0);

      // Used in descriptive output to stdout. These are the first/last values
      // within the measurement window.
      if (!firstStepInMeasureWindow_)
      {
        startACDCmeasureWindow_ = dcSweepVal;
        firstStepInMeasureWindow_ = true;
      }
      endACDCmeasureWindow_ = dcSweepVal;

      if( initialized_  )
      {
        // update integral of outVarValues_[0].  Account for both ascending and descending
        // FROM-TO windows and "direction" (increasing/decreasing) of first swept variable on .DC line.
        if ( ((from_ <= to_) && (dcSweepVal > startSweepValue_)) || ((from_ >= to_) && (dcSweepVal < startSweepValue_)) )
          integralValue_ += 0.5 * (dcSweepVal - lastIndepVarValue_) * (outVarValues_[0] + lastSignalValue_);
        else
          integralValue_ -= 0.5 * (dcSweepVal - lastIndepVarValue_) * (outVarValues_[0] + lastSignalValue_);
      }

      lastIndepVarValue_ = dcSweepVal;
      lastSignalValue_ = outVarValues_[0];
      initialized_=true;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : IntegralEvaluation::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, 3/25/2020
// Creation Date : 3/25/2020
//-----------------------------------------------------------------------------
void IntegralEvaluation::updateAC(
  Parallel::Machine comm,
  const double frequency,
  const Linear::Vector *solnVec,
  const Linear::Vector *imaginaryVec,
  const Util::Op::RFparamsData *RFparams)
{
  // Used in descriptive output to stdout. Store first/last frequency values
  if (!firstSweepValueFound_)
  {
    startSweepValue_ = frequency;
    firstSweepValueFound_ = true;
  }
  endSweepValue_ = frequency;

  if( !calculationDone_ && withinFreqWindow( frequency ) )
  {
    // update our outVarValues_ vector
    updateOutputVars(comm, outVarValues_, frequency, solnVec, 0, 0,
                     imaginaryVec, 0, 0, 0, RFparams);

    // Used in descriptive output to stdout. These are the first/last values
    // within the measurement window.
    if (!firstStepInMeasureWindow_)
    {
      startACDCmeasureWindow_ = frequency;
      firstStepInMeasureWindow_ = true;
    }
    endACDCmeasureWindow_ = frequency;

    if( initialized_  )
    {
      // update integral of outVarValues_[0];
      integralValue_ += 0.5 * (frequency - lastIndepVarValue_) * (outVarValues_[0] + lastSignalValue_);
    }

    lastIndepVarValue_ = frequency;
    lastSignalValue_ = outVarValues_[0];
    initialized_=true;
  }
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

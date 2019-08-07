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
// Purpose       : Transient analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_MeasureMax.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : Max::Max()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
Max::Max(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Base(measureMgr, measureBlock),
  maximumValue_(0.0)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();
}

//-----------------------------------------------------------------------------
// Function      : Max::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void Max::prepareOutputVariables()
{
  // this measurement should have only one dependent variable.
  // Error for now if it doesn't
  numOutVars_ = outputVars_.size();

  if ( numOutVars_ > 1 )
  {
    std::string msg = "Too many dependent variables for statistical measure, \"" + name_ + "\" Exiting.";
    Report::UserFatal() << msg;
  }

  outVarValues_.resize( numOutVars_, 0.0 );
}


//-----------------------------------------------------------------------------
// Function      : Max::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void Max::reset() 
{
  resetBase();
  maximumValue_ = 0.0;
}

//-----------------------------------------------------------------------------
// Function      : Max::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void Max::updateTran(
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
    // we're in the time window, now we need to calculate the value of this
    // measure and see if it triggers any specified rise, fall, cross windowing.
    double tempResult = 0.0;

    // update our outVarValues_ vector
    updateOutputVars(comm, outVarValues_, circuitTime,
      solnVec, stateVec, storeVec, 0, lead_current_vector,
      junction_voltage_vector, lead_current_dqdt_vector, 0);

    // Need to set lastOutputValue_ variable to the current signal value
    // at the first time-step within the measurement window  (That
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
      // Processing needed on the first time step in the
      // RFC window.  If LAST was specified then this is done
      // each time a new RFC window is entered.
      if( !initialized_  || newRiseFallCrossWindowforLast() )
      {
        maximumValue_ = outVarValues_[0];
        calculationInstant_ = circuitTime;
        initialized_ = true;
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

      // calculation of the maximum value
      if( maximumValue_ < outVarValues_[0] )
      {
        maximumValue_ = outVarValues_[0];
        calculationInstant_ = circuitTime;
      }
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : Max::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 4/9/2017
//-----------------------------------------------------------------------------
void Max::updateDC(
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

      if ( !initialized_ )
      {
        maximumValue_ = outVarValues_[0];
        calculationInstant_ = dcSweepVal;
        initialized_ = true;
      }
      else if ( maximumValue_ < outVarValues_[0] )
      {
        maximumValue_ = outVarValues_[0];
        calculationInstant_ = dcSweepVal;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Max::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical Models & Simulation
// Creation Date : 1/29/2019
//-----------------------------------------------------------------------------
void Max::updateAC(
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

    if( !initialized_  )
    {
      maximumValue_ =  outVarValues_[0]; 
      calculationInstant_ = frequency;
      initialized_ = true;
    }
    else if( maximumValue_ <  outVarValues_[0])
    {
      // calculation of the maximum value
      maximumValue_ = outVarValues_[0];
      calculationInstant_ = frequency;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Max::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double Max::getMeasureResult()
{
  if( initialized_ )
  {
    calculationResult_ =  maximumValue_;
  }
  return calculationResult_;
}

//-----------------------------------------------------------------------------
// Function      : Max::printMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 2/09/2015
//-----------------------------------------------------------------------------
std::ostream& Max::printMeasureResult(std::ostream& os, bool printVerbose)
{
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os << std::scientific << std::setprecision(precision_);

  if (!printVerbose)
  {
    if ( !initialized_ && measureMgr_.isMeasFailGiven() && measureMgr_.getMeasFail() )
    {
      // output FAILED to .mt file if .OPTIONS MEASURE MEASFAIL=1 is given in the
      // netlist and this is a failed measure.
      os << name_ << " = FAILED" << std::endl;
    }
    else if ( (measureOutputOption_ == "TIME") || (measureOutputOption_ == "FREQ")
         || (measureOutputOption_ == "SV") )
    {
      // output the time (or frequency or value of the first variable in 
      // the DC sweep vector) when the maximum value occurs.
      os << name_ << " = " << calculationInstant_ << std::endl;
    }
    else
    {
      // output the maximum value
      os << name_ << " = " << this->getMeasureResult() << std::endl;
    }
  }
  else
  {
    if (initialized_)
    {
      os << name_ << " = " << this->getMeasureResult() ;

      // modeStr is "time" for TRAN mode, "freq" for AC mode and 
      // "<sweep variable> value" for DC mode.
      std::string modeStr = setModeStringForMeasureResultText();     
      os << " at " << modeStr << " = " << calculationInstant_ << std::endl;
    }
    else
    { 
      os << name_ << " = FAILED" << std::endl;
    }
  } 

  return os;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

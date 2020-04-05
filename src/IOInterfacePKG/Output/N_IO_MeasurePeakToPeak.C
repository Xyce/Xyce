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

#include <N_IO_MeasurePeakToPeak.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : PeakToPeak::PeakToPeak()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
PeakToPeak::PeakToPeak(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Base(measureMgr, measureBlock),
  maximumValue_(0.0),
  minimumValue_(0.0),
  maximumInstant_(0.0),
  minimumInstant_(0.0)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();
}

//-----------------------------------------------------------------------------
// Function      : PeakToPeak::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void PeakToPeak::prepareOutputVariables() 
{
  // this measurement should have only one dependent variable.
  // Error out if it doesn't
  numOutVars_ = outputVars_.size();

  if ( numOutVars_ > 1 )
  {
    std::string msg = "Too many dependent variables for PP measure, \"" + name_ + "\"";
    Report::UserError0() << msg;
  }

  outVarValues_.resize( numOutVars_, 0.0 );
}


//-----------------------------------------------------------------------------
// Function      : PeakToPeak::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void PeakToPeak::reset() 
{
  resetBase();
  maximumValue_ = 0.0;
  minimumValue_ = 0.0;
  maximumInstant_ = 0.0;
  minimumInstant_ = 0.0;
}

//-----------------------------------------------------------------------------
// Function      : PeakToPeak::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void PeakToPeak::updateTran(
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
    // at the first time-step within the measurement window   (That
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
        minimumValue_ = outVarValues_[0];
        maximumInstant_ = circuitTime;
        minimumInstant_ = circuitTime;
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

      // update the maximum and minimum values.  The Peak-to-Peak value
      // is calculated in the PeakToPeak::getMeasureResult() function.
      if( maximumValue_ < outVarValues_[0] )
      {
        maximumValue_ = outVarValues_[0];
        maximumInstant_ = circuitTime;
      }
      if( minimumValue_ > outVarValues_[0] )
      {
        minimumValue_ = outVarValues_[0];
        minimumInstant_ = circuitTime;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : PeakToPeak::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 4/9/2017
//-----------------------------------------------------------------------------
void PeakToPeak::updateDC(
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
        minimumValue_ = outVarValues_[0];
        maximumInstant_ = dcSweepVal;
        minimumInstant_ = dcSweepVal;
        initialized_ = true;
      }

      // update the maximum and minimum values.  The Peak-to-Peak value
      // is calculated in the PeakToPeak::getMeasureResult() function.
      if ( maximumValue_ < outVarValues_[0] )
      {
        maximumValue_ = outVarValues_[0];
        maximumInstant_ = dcSweepVal;
      }
      if ( minimumValue_ > outVarValues_[0] )
      {
        minimumValue_ = outVarValues_[0];
        minimumInstant_ = dcSweepVal;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : PeakToPeak::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical Models & Simulation
// Creation Date : 1/29/2019
//-----------------------------------------------------------------------------
void PeakToPeak::updateAC(
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
                     imaginaryVec, 0, 0, 0, RFparams );

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
      maximumValue_ = outVarValues_[0];
      minimumValue_ = outVarValues_[0];
      maximumInstant_ = frequency;
      minimumInstant_ = frequency;
      initialized_ = true;
    }

    // update the maximum and minimum values.  The Peak-to-Peak value
    // is calculated in the PeakToPeak::getMeasureResult() function.
    if( maximumValue_ < outVarValues_[0] )
    {
      maximumValue_ = outVarValues_[0];
      maximumInstant_ = frequency;
    }

    if( minimumValue_ > outVarValues_[0] )
    {
      minimumValue_ = outVarValues_[0];
      minimumInstant_ = frequency;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : PeakToPeak::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double PeakToPeak::getMeasureResult()
{
  if( initialized_ )
  {
    calculationResult_ =  maximumValue_ - minimumValue_;
  }
  return calculationResult_;
}

//-----------------------------------------------------------------------------
// Function      : PeakToPeak::printMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 2/09/2015
//-----------------------------------------------------------------------------
std::ostream& PeakToPeak::printMeasureResult(std::ostream& os, bool printVerbose)
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
    else 
    {
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
      os << " with max at " << modeStr << "= " << maximumInstant_ << 
      " and min at " << modeStr << "= " << minimumInstant_<< std::endl;
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

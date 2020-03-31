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

#include <N_IO_MeasureRMS.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : RMS::RMS()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
RMS::RMS(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Base(measureMgr, measureBlock),
  integralXsq_(0.0),
  lastIndepVarValue_(0.0),
  lastSignalValue_(0.0),
  totalIntegrationWindow_(0.0)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();
}

//-----------------------------------------------------------------------------
// Function      : RMS::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void RMS::prepareOutputVariables() 
{
  // this measurement should have only one dependent variable.
  // Error for now if it doesn't
  numOutVars_ = outputVars_.size();

  if ( numOutVars_ > 1 )
  {
    std::string msg = "Too many dependent variables for RMS measure, \"" + name_ + "\"";
    Report::UserError0() << msg;
  }

  outVarValues_.resize( numOutVars_, 0.0 );
}

//-----------------------------------------------------------------------------
// Function      : RMS::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void RMS::reset() 
{
  resetBase();
  integralXsq_ = 0.0;
  totalIntegrationWindow_ = 0.0;
}


//-----------------------------------------------------------------------------
// Function      : RMS::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void RMS::updateTran(
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
      // If LAST was specified then this is done
      // each time a new RFC window is entered.
      if( newRiseFallCrossWindowforLast() )
      {
        integralXsq_ = 0.0;
        totalIntegrationWindow_ = 0.0;
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
        integralXsq_ += 0.5* (circuitTime - lastIndepVarValue_) * (outVarValues_[0]*outVarValues_[0] + lastSignalValue_*lastSignalValue_);
        totalIntegrationWindow_ += (circuitTime - lastIndepVarValue_);
      }

      lastIndepVarValue_ = circuitTime;
      lastSignalValue_ = outVarValues_[0];
      initialized_=true;
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : RMS::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 3/25/2020
//-----------------------------------------------------------------------------
void RMS::updateDC(
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

      if ( initialized_ )
      {
        integralXsq_ += 0.5* (dcSweepVal - lastIndepVarValue_) * (outVarValues_[0]*outVarValues_[0] + lastSignalValue_*lastSignalValue_);
        totalIntegrationWindow_ += (dcSweepVal - lastIndepVarValue_);
      }

      lastIndepVarValue_ = dcSweepVal;
      lastSignalValue_ = outVarValues_[0];
      initialized_=true;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : RMS::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 3/25/2020
//-----------------------------------------------------------------------------
void RMS::updateAC(
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

    if ( initialized_ )
    {
      integralXsq_ += 0.5* (frequency - lastIndepVarValue_) * (outVarValues_[0]*outVarValues_[0] + lastSignalValue_*lastSignalValue_);
      totalIntegrationWindow_ += (frequency - lastIndepVarValue_);
    }

    lastIndepVarValue_ = frequency;
    lastSignalValue_ = outVarValues_[0];
    initialized_=true;
  }
}


//-----------------------------------------------------------------------------
// Function      : RMS::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double RMS::getMeasureResult()
{
  if( initialized_ )
  {
    if (abs(totalIntegrationWindow_) > 0)
      calculationResult_ =  sqrt(integralXsq_ / totalIntegrationWindow_);
    else
    {
      calculationResult_ = calculationDefaultVal_;
      initialized_=false;
    }

  }
  return calculationResult_;
}

//-----------------------------------------------------------------------------
// Function      : RMS::printMeasureWindow
// Purpose       : prints information related to measure window
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/25/2020
//-----------------------------------------------------------------------------
std::ostream& RMS::printMeasureWindow(std::ostream& os, const double indepVarValue)
{
  // Pathological case of FROM=TO within an otherwise valid FROM-TO window.
  // This a failed measure, but the FROM-TO window should be printed correctly.
  if ( (fromGiven_ || toGiven_) && (from_==to_) && firstSweepValueFound_ &&
       ((mode_ == "AC") || (mode_ == "DC")) )
  {
    basic_ios_all_saver<std::ostream::char_type> save(os);
    os << std::scientific << std::setprecision(precision_);
    std::string modeStr = setModeStringForMeasureWindowText();
    os << "Measure Start " << modeStr << "= " << startACDCmeasureWindow_
       << "\tMeasure End " << modeStr << "= " << endACDCmeasureWindow_ << std::endl;
  }
  else
  {
    Base::printMeasureWindow(os,indepVarValue);
  }

  return os;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

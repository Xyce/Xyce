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
// Purpose       :  Find the extrema (max, min or peak-to-peak value)
//                  of a simulation variable
//
// Special Notes : This class contains the functions that are common to
//                 the Max, Min and PeakToPeak classes.  It sits between
//                 those classes and the Base class.
//
// Creator       : Pete Sholander, SNL
// Creation Date : 04/28/2020
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_MeasureExtrema.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : Extrema::Extrema()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 04/28/2020
//-----------------------------------------------------------------------------
Extrema::Extrema(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Base(measureMgr, measureBlock)
{
  // Negative rise/fall/cross values < -1 are not supported for Extrema measures
  if (riseGiven_ && rise_ < -1)
  {
    rise_ = -1;
    Report::UserWarning0() << "RISE value for " << type_ << " measure " << name_
                           << " set to LAST (-1)";
  }
  else if (fallGiven_ && fall_ < -1)
  {
    fall_ = -1;
    Report::UserWarning0() << "FALL value for " << type_ << " measure " << name_
                           << " set to LAST (-1)";
  }
  else if (crossGiven_ && cross_ < -1)
  {
    cross_ = -1;
    Report::UserWarning0() << "CROSS value for " << type_ << " measure " << name_
                           << " set to LAST (-1)";
  }
}

//-----------------------------------------------------------------------------
// Function      : Extrema::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void Extrema::prepareOutputVariables()
{
  // this measurement should have only one dependent variable.
  // Error out if it doesn't
  numOutVars_ = outputVars_.size();

  if ( numOutVars_ > 1 )
  {
    std::string msg = "Too many dependent variables for " + type_ +  " measure, \"" + name_ + "\"";
    Report::UserError0() << msg;
  }

  outVarValues_.resize( numOutVars_, 0.0 );
}

//-----------------------------------------------------------------------------
// Function      : Extrema::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 4/28/2020
//-----------------------------------------------------------------------------
void Extrema::updateTran(
  Parallel::Machine comm,
  double circuitTime,
  double endSimTime,
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
      junction_voltage_vector, lead_current_dqdt_vector, 0, 0, 0, 0);

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
        setMeasureVarsForNewWindow(circuitTime, outVarValues_[0]);
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

      // calculation of the extrema
      updateMeasureVars(circuitTime, outVarValues_[0]);
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : Extrema::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 4/28/2020
//-----------------------------------------------------------------------------
void Extrema::updateDC(
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
    double dcSweepVal = getDCSweepVal(dcParamsVec);
    if (dcParamsVec[0].stepVal < 0)
      dcSweepAscending_=false;

    // Used in descriptive output to stdout. Store name of first variable found in
    // the DC sweep vector.
    sweepVar_ = getDCSweepVarName(dcParamsVec);
    firstSweepValueFound_ = true;

    if( !calculationDone_ && withinDCsweepFromToWindow(dcSweepVal) )
    {
      outVarValues_[0] = getOutputValue(comm, outputVars_[0],
                                        solnVec, stateVec, storeVec, 0,
                                        lead_current_vector,
                                        junction_voltage_vector,
                                        lead_current_dqdt_vector, 0, 0, 0, 0);

      if ( !initialized_ )
        setMeasureVarsForNewWindow(dcSweepVal, outVarValues_[0]);
      else
        updateMeasureVars(dcSweepVal, outVarValues_[0]);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Extrema::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical Models & Simulation
// Creation Date : 1/29/2019
//-----------------------------------------------------------------------------
void Extrema::updateAC(
  Parallel::Machine comm,
  double frequency,
  double fStart,
  double fStop,
  const Linear::Vector *solnVec,
  const Linear::Vector *imaginaryVec,
  const Util::Op::RFparamsData *RFparams)
{
  // Used for descriptive output to stdout.
  firstSweepValueFound_ = true;

  if( !calculationDone_ && withinFreqWindow( frequency ) )
  {
    // update our outVarValues_ vector 
    updateOutputVars(comm, outVarValues_, frequency, solnVec, 0, 0,
                     imaginaryVec, 0, 0, 0, 0, 0, 0, RFparams);

    if ( !initialized_ )
      setMeasureVarsForNewWindow(frequency, outVarValues_[0]);
    else
      updateMeasureVars(frequency, outVarValues_[0]);
  }
}

//-----------------------------------------------------------------------------
// Function      : Extrema::updateNoise()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical Models & Simulation
// Creation Date : 5/11/2020
//-----------------------------------------------------------------------------
void Extrema::updateNoise(
  Parallel::Machine comm,
  double frequency,
  double fStart,
  double fStop,
  const Linear::Vector *solnVec,
  const Linear::Vector *imaginaryVec,
  double totalOutputNoiseDens,
  double totalInputNoiseDens,
  const std::vector<Xyce::Analysis::NoiseData*> *noiseDataVec)
{
  // Used for descriptive output to stdout
  firstSweepValueFound_ = true;

  if( !calculationDone_ && withinFreqWindow( frequency ) )
  {
    // update our outVarValues_ vector
    updateOutputVars(comm, outVarValues_, frequency, solnVec, 0, 0,
                     imaginaryVec, 0, 0, 0,
                     totalOutputNoiseDens, totalInputNoiseDens, noiseDataVec, 0);

    if ( !initialized_ )
      setMeasureVarsForNewWindow(frequency, outVarValues_[0]);
    else
      updateMeasureVars(frequency, outVarValues_[0]);
  }
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

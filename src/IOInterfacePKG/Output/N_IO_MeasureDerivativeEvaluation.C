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

#include <N_IO_MeasureDerivativeEvaluation.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::DerivativeEvaluation()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
DerivativeEvaluation::DerivativeEvaluation(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Base(measureMgr, measureBlock),
  lastIndepVarValue_(0.0),
  lastDepVarValue_(0.0),
  lastTargValue_(0.0),
  numPointsFound_(0),
  whenIdx_(1)  
{
  // Note: the default value for whenIdx_ is 1.  This is for the default WHEN syntax,
  // which is WHEN v(a)=<val>
  if (atGiven_)
  {
    // the whenIdx_ is not actually used in the control logic in updateTran() for
    // the AT= syntax.  However, for AT= case, the output variable is in position 0 
    // of the outputVars_ vector.  Setting whenIdx_ = 0 prevents some invalid read 
    // accesses later in the updateTran() function when AT is used.
    whenIdx_=0;
  }  

  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void DerivativeEvaluation::prepareOutputVariables() 
{
  // If AT keyword is given then numOutVars should only have one entry
  numOutVars_ = outputVars_.size();

  if ( (numOutVars_ > 1) && atGiven_ )
  {
    std::string msg = "Too many dependent variables for DERIV measure, \"" + name_ + "\"";
    Report::UserError0() << msg;
  }

  outVarValues_.resize( numOutVars_, 0.0 );
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void DerivativeEvaluation::reset() 
{
  resetBase();
  lastIndepVarValue_=0.0;
  lastDepVarValue_=0.0;
  lastOutputVarValue_=0.0;
  lastTargValue_=0.0;
  numPointsFound_=0;
}


//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void DerivativeEvaluation::updateTran(
  Parallel::Machine comm,
  const double circuitTime,
  const Linear::Vector *solnVec,
  const Linear::Vector *stateVec,
  const Linear::Vector *storeVec,
  const Linear::Vector *lead_current_vector,
  const Linear::Vector *junction_voltage_vector,
  const Linear::Vector *lead_current_dqdt_vector)
{
  if (!calculationDone_ && withinTimeWindow( circuitTime ))
  {
    // update our outVarValues_ vector
    updateOutputVars(comm, outVarValues_, circuitTime,
      solnVec, stateVec, storeVec, 0, lead_current_vector,
      junction_voltage_vector, lead_current_dqdt_vector, 0, 0, 0, 0);

    if( !initialized_ )
    {
      setMeasureState(circuitTime);
    }

    if (atGiven_)
    { 
      // check and see if last point and this point bound the target point 
      double backDiff    = lastIndepVarValue_ - at_;
      double forwardDiff = circuitTime - at_;
    
      // if we bound the time target then either
      //  (backDiff < 0) && (forwardDiff > 0)  
      //   OR
      //  (backDiff > 0) && (forwardDiff < 0) 
      // or more simply sgn( backDiff ) = - sgn( forwardDiff )
      //
      // Also test for equality, to within the minval_ tolerance, as with the WHEN syntax.
      // Only do this equality test if the simulation time is > 0, since the
      // three-point difference approx. for the derivative needs a "previous point". 
      if( ((backDiff < 0.0) && (forwardDiff > 0.0)) || ((backDiff > 0.0) && (forwardDiff < 0.0)) ||
	  (((abs(backDiff) < minval_) || (abs(forwardDiff) < minval_)) && circuitTime > 0) )
      {
        updateCalculationResult(circuitTime);
        calculationDone_ = true;
      }
    }
    else if ( whenGiven_ && outputValueTargetGiven_ || (numOutVars_ == 3) )
    {
      double targVal= updateTargVal();

      // Need to set lastOutputValue_ variable to the current signal value
      // at the first time-step within the measurement window.  (That
      // window is set by the TO-FROM and TD qualifiers if present.)  This is 
      // needed so that the RISE/FALL/CROSS count is not incremented at time=0, if
      // the measured waveform has a DC offset at time=0    
      if (!firstStepInMeasureWindow_)
      {
        lastOutputValue_ = outVarValues_[whenIdx_]; 
        firstStepInMeasureWindow_ = true;
      }

      // for the WHEN qualifier, the rfc level used is either the target value of the WHEN
      // clause, or the value of the RFC_LEVEL qualifier if one is specified.
      if ( withinRiseFallCrossWindow( outVarValues_[whenIdx_], (rfcLevelGiven_ ? rfcLevel_ : targVal) ) )
      {
        // If LAST was specified then this is done
        // each time a new RFC window is entered.
        if( newRiseFallCrossWindowforLast() )
        {
          resultFound_ = false;
          calculationResult_= calculationDefaultVal_;
          firstStepInRfcWindow_ = false;
	  //Xyce::dout() << "found new rfc window at time= " << circuitTime << std::endl;
        }
        // only record the start time of the RFC window for the DERIVATIVE Measure.  The end
        // time is not recorded because the measure becomes inactive once the WHEN or AT clause
        // is satisfied.
        if( !firstStepInRfcWindow_  )
        {
          firstStepInRfcWindow_ = true;
          rfcWindowFound_ = true;
          rfcWindowStartTime_ = circuitTime;
        }

        // process WHEN qualifer
        if ( isWHENcondition(circuitTime, targVal) )
	{
          calculationInstant_ = interpolateCalculationInstant(circuitTime, targVal);
          updateCalculationResult(circuitTime);
          calculationDone_ = !measureLastRFC_;
          resultFound_=true;
        }
      }
    }
  }

  updateMeasureState(circuitTime);
}


//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 3/25/2020
//-----------------------------------------------------------------------------
void DerivativeEvaluation::updateDC(
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

    // Used in descriptive output to stdout. Store name of the first variable found in
    // the DC sweep vector.
    sweepVar_= getDCSweepVarName(dcParamsVec);
    firstSweepValueFound_ = true;
    ++numPointsFound_;

    // update our outVarValues_ vector
    updateOutputVars(comm, outVarValues_, dcSweepVal,
      solnVec, stateVec, storeVec, 0, lead_current_vector,
      junction_voltage_vector, lead_current_dqdt_vector, 0, 0, 0,0 );

    // Consider all intervals since the WHEN clause uses the interpolated calculation
    // instant.  We check that the interpolated WHEN time is within the measurement
    // window each time one is found.
    if (!calculationDone_)
    {
      if (atGiven_ && (numPointsFound_ > 1) && withinDCsweepFromToWindow(at_) )
      {
        // process AT qualifer. The AT value must be within the measurement window.
        if (isATforACDCNoise(dcSweepVal))
        {
          updateCalculationResult(dcSweepVal);
          calculationDone_ = true;
        }
      }
      else if ( (numPointsFound_ > 1) && (outputValueTargetGiven_ || (numOutVars_ == 3)) )
      {
        double targVal=updateTargVal();

        // process WHEN qualifer
        if ( isWHENcondition(dcSweepVal, targVal) )
        {
          // use same time interpolation algorithm as FIND-WHEN measure
          double whenTime = interpolateCalculationInstant(dcSweepVal, targVal);
          if ( withinDCsweepFromToWindow(whenTime) )
	  {
            calculationInstant_ = whenTime;
            updateCalculationResult(dcSweepVal);
            calculationDone_ = true;
            resultFound_=true;
          }
        }
      }
    }

    initialized_=true;
    updateMeasureState(dcSweepVal);
  }
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical Models & Simulation
// Creation Date : 3/25/2020
//-----------------------------------------------------------------------------
void DerivativeEvaluation::updateAC(
  Parallel::Machine comm,
  const double frequency,
  const Linear::Vector *solnVec,
  const Linear::Vector *imaginaryVec,
  const Util::Op::RFparamsData *RFparams)
{
  // Used in descriptive output to stdout
  firstSweepValueFound_ = true;
  ++numPointsFound_;

  // update our outVarValues_ vector
  updateOutputVars(comm, outVarValues_, frequency, solnVec, 0, 0,
                   imaginaryVec, 0, 0, 0, 0, 0, 0, RFparams);

  // Consider all intervals since the WHEN clause uses the interpolated calculation
  // instant.  We check that the interpolated WHEN time is within the measurement
  // window each time one is found.
  if (!calculationDone_)
  {
    if (atGiven_ && (numPointsFound_ > 1) && withinFreqWindow(at_))
    {
      // Process AT qualifer.  The AT value must be within the measurement window.
      if (isATforACDCNoise(frequency))
      {
        updateCalculationResult(frequency);
        calculationDone_ = true;
      }
    }
    else if (whenGiven_ && (numPointsFound_ > 1) && (outputValueTargetGiven_ || (numOutVars_ == 3)) )
    {
      double targVal=updateTargVal();

      // Process WHEN qualifer, using the same time interpolation algorithm as FIND-WHEN measure.
      // The measure will succeed if the interpolated WHEN time is inside of the measurement window.
      if (isWHENcondition(frequency, targVal))
      {
        double whenTime = interpolateCalculationInstant(frequency, targVal);
        if ( withinFreqWindow(whenTime) )
	{
          calculationInstant_ = whenTime;
          updateCalculationResult(frequency);
          calculationDone_ = true;
          resultFound_=true;
        }
      }
    }
  }

  initialized_ = true;
  updateMeasureState(frequency);
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::updateNoise()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical Models & Simulation
// Creation Date : 5/12/2010
//-----------------------------------------------------------------------------
void DerivativeEvaluation::updateNoise(
  Parallel::Machine comm,
  const double frequency,
  const Linear::Vector *solnVec,
  const Linear::Vector *imaginaryVec,
  const double totalOutputNoiseDens,
  const double totalInputNoiseDens,
  const std::vector<Xyce::Analysis::NoiseData*> *noiseDataVec)
{
  // Used in descriptive output to stdout
  firstSweepValueFound_ = true;
  ++numPointsFound_;

  // update our outVarValues_ vector
  updateOutputVars(comm, outVarValues_, frequency, solnVec, 0, 0,
                   imaginaryVec, 0, 0, 0,
                   totalOutputNoiseDens, totalInputNoiseDens, noiseDataVec, 0);

  // Consider all intervals since the WHEN clause uses the interpolated calculation
  // instant.  We check that the interpolated WHEN time is within the measurement
  // window each time one is found.
  if (!calculationDone_)
  {
    if (atGiven_ && (numPointsFound_ > 1) && withinFreqWindow(at_))
    {
      // Process AT qualifer.  The AT value must be within the measurement window.
      if (isATforACDCNoise(frequency))
      {
        updateCalculationResult(frequency);
        calculationDone_ = true;
      }
    }
    else if (whenGiven_ && (numPointsFound_ > 1) && (outputValueTargetGiven_ || (numOutVars_ == 3)) )
    {
      double targVal=updateTargVal();

      // Process WHEN qualifer, using the same time interpolation algorithm as FIND-WHEN measure.
      // The measure will succeed if the interpolated WHEN time is inside of the measurement window.
      if ( isWHENcondition(frequency, targVal) )
      {
        // use same time interpolation algorithm as FIND-WHEN measure
        double whenTime = interpolateCalculationInstant(frequency, targVal);
        if ( withinFreqWindow(whenTime) )
	{
          calculationInstant_ = whenTime;
          updateCalculationResult(frequency);
          calculationDone_ = true;
          resultFound_=true;
        }
      }
    }
  }

  initialized_ = true;
  updateMeasureState(frequency);
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::updateCalculationResult()
// Purpose       :
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 05/22/2020
//-----------------------------------------------------------------------------
void DerivativeEvaluation::updateCalculationResult(const double indepVarVal)
{
  // asymmetrical 3-point approximation for first derivative.
  calculationResult_ = (outVarValues_[0] - lastOutputVarValue_) / (indepVarVal - lastIndepVarValue_);

  return;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::isATforACDCNoise
// Purpose       : Evaluates if the AT condition is true for AC, DC and NOISE modes.
// Special Notes : For AC and NOISE measures, the independent variable is
//                 frequency.  For DC measures, it is the value of the first
//                 variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 05/21/2020
//-----------------------------------------------------------------------------
bool DerivativeEvaluation::isATforACDCNoise(const double indepVarVal)
{
  // check and see if last point and this point bound the target point
  double backDiff    = lastIndepVarValue_ - at_;
  double forwardDiff = indepVarVal - at_;

  // if we bound the frequency target then either
  //  (backDiff < 0) && (forwardDiff > 0)
  //   OR
  //  (backDiff > 0) && (forwardDiff < 0)
  // or more simply sgn( backDiff ) = - sgn( forwardDiff )
  //
  // Also test for equality, to within the minval_ tolerance, as with the WHEN syntax.
  return ( ((backDiff < 0.0) && (forwardDiff > 0.0)) || ((backDiff > 0.0) && (forwardDiff < 0.0)) ||
	   (((abs(backDiff) < minval_) || (abs(forwardDiff) < minval_))) );
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::isWHENcondition
// Purpose       : Evaluates if the WHEN condition is true for all modes
// Special Notes : For AC and NOISE measures, the independent variable is
//                 frequency.  For DC measures, it is the value of the first
//                 variable in the DC sweep vector.  For TRAN measures, it
//                 is the circuit time.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 05/21/2020
//-----------------------------------------------------------------------------
bool DerivativeEvaluation::isWHENcondition(const double indepVarVal, const double targVal)
{
  bool whenFound=false;

  if (!resultFound_)
  {
    // this is the simple case where Xyce output a value within tolerance
    // of the target value
    if( fabs(outVarValues_[whenIdx_] - targVal) < minval_ )
    {
      whenFound=true;
    }
    else
    {
      // check and see if last point and this point bound the target point
      double backDiff    = lastDepVarValue_ - targVal;
      double forwardDiff = outVarValues_[whenIdx_] - targVal;

      // if we bound the target then either
      //  (backDiff < 0) && (forwardDiff > 0)
      //   OR
      //  (backDiff > 0) && (forwardDiff < 0)
      // or more simply sgn( backDiff ) = - sgn( forwardDiff )
      if( ((backDiff < 0.0) && (forwardDiff > 0.0)) || ((backDiff > 0.0) && (forwardDiff < 0.0)) ||
          (backDiff == 0.0) )
      {
        whenFound = true;
      }
    }
  }

  return whenFound;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::setMeasureState()
// Purpose       : initializes the past values of the independent, dependent
//                 and measured variables, as well as the past target level.
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 05/05/2020
//-----------------------------------------------------------------------------
void DerivativeEvaluation::setMeasureState(const double indepVarVal)
{
  // assigned last dependent and independent var to current value of the independent
  // varible and outVarValue_[whenIdx_].  While we can't interpolate on this step, it
  // ensures that the initial history is something realistic
  lastIndepVarValue_=indepVarVal;
  lastDepVarValue_=outVarValues_[whenIdx_];
  lastOutputVarValue_=outVarValues_[0];
  if (outputValueTargetGiven_)
    lastTargValue_ = outputValueTarget_;
  else
    lastTargValue_ = outVarValues_[whenIdx_+1];

  initialized_=true;

  return;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::updateMeasureState()
// Purpose       : updates the past values of the independent, dependent
//                 and measured variables, as well as the past target level.
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 05/05/2020
//-----------------------------------------------------------------------------
void DerivativeEvaluation::updateMeasureState(const double indepVarVal)
{
  lastIndepVarValue_ = indepVarVal;
  lastDepVarValue_ = outVarValues_[whenIdx_];
  lastOutputVarValue_=outVarValues_[0];
  if (outputValueTargetGiven_)
    lastTargValue_ = outputValueTarget_;
  else
    lastTargValue_ = outVarValues_[whenIdx_+1];

  return;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::updateTargVal
// Purpose       : updates the target value for the WHEN clause
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 05/05/2020
//-----------------------------------------------------------------------------
double DerivativeEvaluation::updateTargVal()
{
  double targVal = 0.0;

  if( outputValueTargetGiven_ )
  {
    // This is the form WHEN v(a)=fixed value
    targVal = outputValueTarget_;
  }
  else
  {
    // This is the form WHEN v(a)= potentially changing value, such as v(a)=v(b)
    // in that case v(b) is in outVarValues_[whenIdx_+1]
    targVal = outVarValues_[whenIdx_+1];
  }

  return targVal;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double DerivativeEvaluation::getMeasureResult()
{
  return calculationResult_;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::printMeasureResult()
// Purpose       : used to print the measurement result to an output stream
//                 object, which is typically the mt0, ma0 or ms0 file
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 2/22/2015
//-----------------------------------------------------------------------------
std::ostream& DerivativeEvaluation::printMeasureResult(std::ostream& os)
{
    basic_ios_all_saver<std::ostream::char_type> save(os);
    os << std::scientific << std::setprecision(precision_);

    if ( !calculationDone_ && measureMgr_.isMeasFailGiven() && measureMgr_.getMeasFail() )
    {
      // output FAILED to .mt file if .OPTIONS MEASURE MEASFAIL=1 is given in the
      // netlist and this is a failed measure.
      os << name_ << " = FAILED" << std::endl;
    }
    else
    {
      os << name_ << " = " << this->getMeasureResult() << std::endl;
    }

    return os;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::printVerboseMeasureResult()
// Purpose       : used to print the "verbose" (more descriptive) measurement
//                 result to an output stream object, which is typically stdout
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 2/22/2015
//-----------------------------------------------------------------------------
std::ostream& DerivativeEvaluation::printVerboseMeasureResult(std::ostream& os)
{
    basic_ios_all_saver<std::ostream::char_type> save(os);
    os << std::scientific << std::setprecision(precision_);

    if (calculationDone_ || ( measureLastRFC_ && resultFound_ ) )
    {
      os << name_ << " = " << this->getMeasureResult() ;
      if (atGiven_)
      {
        os << " for AT = " << at_;
      }
      else
      {
        // modeStr is "time" for TRAN mode, "freq" for AC and NOISE modes, and
        // "<sweep variable> value" for DC mode.
        std::string modeStr = setModeStringForMeasureResultText();
        os << " at " << modeStr << " = " << calculationInstant_;
      }
    }
    else
    { 
      os << name_ << " = FAILED";
      if (atGiven_)
      { 
        os << " for AT = " << at_;
      }
    }
    os << std::endl;

    return os;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::printMeasureWindow
// Purpose       : prints information related to measure window
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 02/5/2015
//-----------------------------------------------------------------------------
std::ostream& DerivativeEvaluation::printMeasureWindow(std::ostream& os, const double endSimTime,
				                       const double startSweepVal, const double endSweepVal)
{
  //no op
  if (!atGiven_) 
  {
    Base::printMeasureWindow(os, endSimTime, startSweepVal, endSweepVal);
  }
  else
  {
    // no op if AT keyword was given
  }

  return os;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::printRFCWindow()
// Purpose       : print informaiton about the start time of the RISE or FALL
//                 window, if a valid one was found.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 09/21/2015
//-----------------------------------------------------------------------------
std::ostream& DerivativeEvaluation::printRFCWindow(std::ostream& os)
{  
  // no op, for any measure that supports WHEN     
  
  return os;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::interpolateCalculationInstant
// Purpose       : Interpolate the time for when the measure is satisifed.
//                 This accounts for case of WHEN V(1)=V(2) where both
//                 variables may be changing.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/25/2020
//-----------------------------------------------------------------------------
double DerivativeEvaluation::interpolateCalculationInstant(double currIndepVarValue, double targVal)
{
  // Calculate slopes and y-intercepts of the two lines, to get them into
  // canonical y=ax+c and y=bx+d form.  If the WHEN clause is of the form
  // WHEN V(1)=V(2) then the line with (a,c) is the value of V(1), which is the
  // "dependent variable".  The line with (b,d) is the value of V(2), which
  // is the "target value".
  double a = (outVarValues_[whenIdx_] - lastDepVarValue_)/(currIndepVarValue - lastIndepVarValue_);
  double b = (targVal - lastTargValue_)/(currIndepVarValue - lastIndepVarValue_);
  double c = outVarValues_[whenIdx_] - a*currIndepVarValue;
  double d = targVal - b*currIndepVarValue;

  // This is the algebra for when the time, frequency or DC sweep value when the two non-parallel
  // lines associated with WHEN clause intersect.
  double calcInstant = (d-c)/(a-b);

  return calcInstant;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

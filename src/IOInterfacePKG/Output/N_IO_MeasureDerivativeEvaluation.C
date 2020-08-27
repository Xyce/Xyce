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
// Function      : DerivativeEvaluationBase::DerivativeEvaluationBase()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
DerivativeEvaluationBase::DerivativeEvaluationBase(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Base(measureMgr, measureBlock),
  lastIndepVarValue_(0.0),
  lastDepVarValue_(0.0),
  lastTargValue_(0.0),
  startDCMeasureWindow_(0.0),
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

  // hard code this measure type to ignore the RFC_LEVEL qualifier
  if (rfcLevelGiven_)
  {
    rfcLevelGiven_=false;
    Xyce::Report::UserWarning0() << "RFC_LEVEL will be ignored for measure " << name_;
  }
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationBase::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void DerivativeEvaluationBase::prepareOutputVariables() 
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
// Function      : DerivativeEvaluationBase::resetDerivativeBase()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void DerivativeEvaluationBase::resetDerivativeEvaluationBase()
{
  resetBase();
  lastIndepVarValue_=0.0;
  lastDepVarValue_=0.0;
  lastOutputVarValue_=0.0;
  lastTargValue_=0.0;
  startDCMeasureWindow_=0;
  numPointsFound_=0;
  calculationResultVec_.clear();
  calculationInstantVec_.clear();
}


//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationBase::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void DerivativeEvaluationBase::updateTran(
  Parallel::Machine comm,
  const double circuitTime,
  const double endSimTime,
  const Linear::Vector *solnVec,
  const Linear::Vector *stateVec,
  const Linear::Vector *storeVec,
  const Linear::Vector *lead_current_vector,
  const Linear::Vector *junction_voltage_vector,
  const Linear::Vector *lead_current_dqdt_vector)
{
  numPointsFound_++;

  // update our outVarValues_ vector
  updateOutputVars(comm, outVarValues_, circuitTime,
    solnVec, stateVec, storeVec, 0, lead_current_vector,
    junction_voltage_vector, lead_current_dqdt_vector, 0, 0, 0, 0);

  if (numPointsFound_ == 1)
    setMeasureState(circuitTime);

  if ( !calculationDone_ && !isInvalidTimeWindow(endSimTime) )
  {
    initialized_ = true;

    if (atGiven_ && (numPointsFound_ > 1) && withinTimeWindow(at_))
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
        calculationResult_=getDerivativeValue(circuitTime);
        calculationDone_ = true;
        resultFound_ = true;
      }
    }
    else if ( whenGiven_ && (numPointsFound_ > 1) )
    {
      double targVal= updateTargVal();

      // The measure will succeed if the interpolated WHEN time is inside of the
      // measurement window, and the RISE/FALL/CROSS condition is met.
      if (isWHENcondition(circuitTime, targVal))
      {
        double whenTime = interpolateCalculationInstant(circuitTime, targVal);
        if (withinTimeWindow(whenTime))
	{
          updateRFCcountForWhen();
          if (withinRFCWindowForWhen())
	  {
            updateCalculationInstant(whenTime);
            updateCalculationResult(circuitTime);
            calculationDone_ = !measureLastRFC_;
            resultFound_=true;
          }
        }
      }
    }
  }

  updateMeasureState(circuitTime);
}


//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationBase::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 3/25/2020
//-----------------------------------------------------------------------------
void DerivativeEvaluationBase::updateDC(
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

    if (!firstStepInMeasureWindow_)
    {
      startDCMeasureWindow_ = dcSweepVal;
      firstStepInMeasureWindow_ = true;
    }

    // The second part of this conditional is needed to deal with multiple sweep
    // variables.  We need to reset the last value variables, each time the first
    // sweep variable rolls over.
    if( (numPointsFound_ == 1) || (initialized_ && dcSweepVal == startDCMeasureWindow_) )
    {
        setMeasureState(dcSweepVal);
    }

    // Consider all intervals since the WHEN clause uses the interpolated calculation
    // instant.  We check that the interpolated WHEN value is within the measurement
    // window when each one is found.
    if ( !calculationDone_ && !isInvalidDCsweepWindow(dcParamsVec[0].startVal, dcParamsVec[0].stopVal) )
    {
      initialized_=true;
      if (atGiven_ && (numPointsFound_ > 1) && withinDCsweepFromToWindow(at_) )
      {
        // process AT qualifer. The AT value must be within the measurement window.
        if (isATforACDCNoise(dcSweepVal))
        {
          calculationResult_=getDerivativeValue(dcSweepVal);
          calculationDone_ = true;
          resultFound_ = true;
        }
      }
      else if (whenGiven_ && (numPointsFound_ > 1))
      {
        double targVal=updateTargVal();

        // The measure will succeed if the interpolated WHEN time is inside of the
        // measurement window, and the RISE/FALL/CROSS condition is then met.
        if ( isWHENcondition(dcSweepVal, targVal) )
        {
          // use same time interpolation algorithm as FIND-WHEN measure
          double whenSweepVal = interpolateCalculationInstant(dcSweepVal, targVal);
          if ( withinDCsweepFromToWindow(whenSweepVal) )
	  {
            updateRFCcountForWhen();
            if (withinRFCWindowForWhen())
	    {
              updateCalculationInstant(whenSweepVal);
              updateCalculationResult(dcSweepVal);
              calculationDone_ = !measureLastRFC_;
              resultFound_=true;
            }
          }
        }
      }
    }

    updateMeasureState(dcSweepVal);
  }
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationBase::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical Models & Simulation
// Creation Date : 3/25/2020
//-----------------------------------------------------------------------------
void DerivativeEvaluationBase::updateAC(
  Parallel::Machine comm,
  const double frequency,
  const double fStart,
  const double fStop,
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

  if (numPointsFound_ == 1)
    setMeasureState(frequency);

  // Consider all intervals since the WHEN clause uses the interpolated calculation
  // instant.  We check that the interpolated WHEN time is within the measurement
  // window each time one is found.
  if ( !calculationDone_ && !isInvalidFreqWindow(fStart,fStop) )
  {
    initialized_ = true;
    if (atGiven_ && (numPointsFound_ > 1) && withinFreqWindow(at_))
    {
      // Process AT qualifer.  The AT value must be within the measurement window.
      if (isATforACDCNoise(frequency))
      {
        calculationResult_=getDerivativeValue(frequency);
        calculationDone_ = true;
        resultFound_ = true;
      }
    }
    else if (whenGiven_ && (numPointsFound_ > 1))
    {
      double targVal=updateTargVal();

      // Process WHEN qualifer, using the same time interpolation algorithm as FIND-WHEN measure.
      // The measure will succeed if the interpolated WHEN frequency is inside of the
      //  measurement window, and the RISE/FALL/CROSS condition is then met.
      if (isWHENcondition(frequency, targVal))
      {
        double whenFreq = interpolateCalculationInstant(frequency, targVal);
        if ( withinFreqWindow(whenFreq) )
	{
          updateRFCcountForWhen();
          if (withinRFCWindowForWhen())
	  {
            updateCalculationInstant(whenFreq);
            updateCalculationResult(frequency);
            calculationDone_ = !measureLastRFC_;
            resultFound_=true;
          }
        }
      }
    }
  }

  updateMeasureState(frequency);
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationBase::updateNoise()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical Models & Simulation
// Creation Date : 5/12/2010
//-----------------------------------------------------------------------------
void DerivativeEvaluationBase::updateNoise(
  Parallel::Machine comm,
  const double frequency,
  const double fStart,
  const double fStop,
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

  if (numPointsFound_ == 1)
    setMeasureState(frequency);

  // Consider all intervals since the WHEN clause uses the interpolated calculation
  // instant.  We check that the interpolated WHEN time is within the measurement
  // window each time one is found.
  if ( !calculationDone_ && !isInvalidFreqWindow(fStart,fStop) )
  {
    initialized_ = true;
    if (atGiven_ && (numPointsFound_ > 1) && withinFreqWindow(at_))
    {
      // Process AT qualifer.  The AT value must be within the measurement window.
      if (isATforACDCNoise(frequency))
      {
        calculationResult_=getDerivativeValue(frequency);
        calculationDone_ = true;
        resultFound_ = true;
      }
    }
    else if (whenGiven_ && (numPointsFound_ > 1))
    {
      double targVal=updateTargVal();

      // Process WHEN qualifer, using the same time interpolation algorithm as FIND-WHEN measure.
      // The measure will succeed if the interpolated WHEN frequency is inside of the
      // measurement window, and the RISE/FALL/CROSS condition is then met.
      if ( isWHENcondition(frequency, targVal) )
      {
        // use same time interpolation algorithm as FIND-WHEN measure
        double whenFreq = interpolateCalculationInstant(frequency, targVal);
        if ( withinFreqWindow(whenFreq) )
	{
           updateRFCcountForWhen();
          if (withinRFCWindowForWhen())
	  {
            updateCalculationInstant(whenFreq);
            updateCalculationResult(frequency);
            calculationDone_ = !measureLastRFC_;
            resultFound_=true;
          }
        }
      }
    }
  }

  updateMeasureState(frequency);
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationBase::isATforACDCNoise
// Purpose       : Evaluates if the AT condition is true for AC, DC and NOISE modes.
// Special Notes : For AC and NOISE measures, the independent variable is
//                 frequency.  For DC measures, it is the value of the first
//                 variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 05/21/2020
//-----------------------------------------------------------------------------
bool DerivativeEvaluationBase::isATforACDCNoise(const double indepVarVal)
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
// Function      : DerivativeEvaluationBase::isWHENcondition
// Purpose       : Evaluates if the WHEN condition is true for all modes
// Special Notes : For AC and NOISE measures, the independent variable is
//                 frequency.  For DC measures, it is the value of the first
//                 variable in the DC sweep vector.  For TRAN measures, it
//                 is the circuit time.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 05/21/2020
//-----------------------------------------------------------------------------
bool DerivativeEvaluationBase::isWHENcondition(const double indepVarVal, const double targVal)
{
  bool whenFound=false;

  //if (!resultFound_)
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
// Function      : DerivativeEvaluation::updateRFCcountForWhen()
// Purpose       : Updates the rise, fall and cross counts.
// Special Notes : This function is normally only called for a WHEN time that
//                 falls within the FROM-TO window.  So, the WHEN time is always
//                 a valid cross.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/21/2020
//-----------------------------------------------------------------------------
void DerivativeEvaluationBase::updateRFCcountForWhen()
{
  actualCross_++;
  if (outVarValues_[whenIdx_] > lastDepVarValue_)
    actualRise_++;
  else
    actualFall_++;

  return;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationBase::withinRFCWindowForWhen()
// Purpose       : Determine if the current WHEN time falls within the specified
//                 RFC value.
// Special Notes : Assumes that a WHEN measure without RISE, FALL or CROSS given
//                 on the .MEASURE line defaults to CROSS given.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/21/2020
//-----------------------------------------------------------------------------
bool DerivativeEvaluationBase::withinRFCWindowForWhen()
{
  bool retVal=false;

  if (riseGiven_ && (outVarValues_[whenIdx_] > lastDepVarValue_) && (actualRise_ >= rise_))
    retVal=true;
  else if (fallGiven_ && (outVarValues_[whenIdx_] < lastDepVarValue_) && (actualFall_ >= fall_))
    retVal=true;
  else if ( !(riseGiven_ || fallGiven_) && (actualCross_ >= cross_) )
    retVal=true;

  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationBase::setMeasureState()
// Purpose       : initializes the past values of the independent, dependent
//                 and measured variables, as well as the past target level.
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 05/05/2020
//-----------------------------------------------------------------------------
void DerivativeEvaluationBase::setMeasureState(const double indepVarVal)
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

  return;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationBase::updateMeasureState()
// Purpose       : updates the past values of the independent, dependent
//                 and measured variables, as well as the past target level.
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 05/05/2020
//-----------------------------------------------------------------------------
void DerivativeEvaluationBase::updateMeasureState(const double indepVarVal)
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
// Function      : DerivativeEvaluationBase::updateTargVal
// Purpose       : updates the target value for the WHEN clause
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 05/05/2020
//-----------------------------------------------------------------------------
double DerivativeEvaluationBase::updateTargVal()
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
// Function      : DerivativeEvaluationBase::printMeasureWindow
// Purpose       : prints information related to measure window
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 02/5/2015
//-----------------------------------------------------------------------------
std::ostream& DerivativeEvaluationBase::printMeasureWindow(std::ostream& os, const double endSimTime,
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
// Function      : DerivativeEvaluationBase::printRFCWindow()
// Purpose       : print informaiton about the start time of the RISE or FALL
//                 window, if a valid one was found.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 09/21/2015
//-----------------------------------------------------------------------------
std::ostream& DerivativeEvaluationBase::printRFCWindow(std::ostream& os)
{
  // no op, for any measure that supports WHEN

  return os;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationBase::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double DerivativeEvaluationBase::getMeasureResult()
{
  return calculationResult_;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationBase::interpolateCalculationInstant
// Purpose       : Interpolate the time for when the measure is satisifed.
//                 This accounts for case of WHEN V(1)=V(2) where both
//                 variables may be changing.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/25/2020
//-----------------------------------------------------------------------------
double DerivativeEvaluationBase::interpolateCalculationInstant(double currIndepVarValue, double targVal)
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

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationBase::getDerivativeValue
// Purpose       :
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 05/22/2020
//-----------------------------------------------------------------------------
double DerivativeEvaluationBase::getDerivativeValue(const double currIndepVarVal)
{
  // asymmetrical 3-point approximation for first derivative.
  return (outVarValues_[0] - lastOutputVarValue_) / (currIndepVarVal - lastIndepVarValue_);
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::DerivativeEvaluation()
// Purpose       : Constructor
// Special Notes : Non-continuous version of DERIV measure, that returns only
//                 only one measure value.  This class is invoked by
//                 .MEASURE TRAN, .MEASURE AC, .MEASURE NOISE and .MEASURE DC.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/21/2020
//-----------------------------------------------------------------------------
DerivativeEvaluation::DerivativeEvaluation(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  DerivativeEvaluationBase(measureMgr, measureBlock),
  RFC_(0)
{
  if (riseGiven_)
    RFC_=rise_;
  else if (fallGiven_)
    RFC_=fall_;
  else if (crossGiven_)
    RFC_=cross_;
  else
  {
    // default case when RISE, FALL or CROSS is not explicitly given on .MEASURE line
    crossGiven_=true;
    cross_=0;
    RFC_=0;
  }
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/21/2020
//-----------------------------------------------------------------------------
void DerivativeEvaluation::reset()
{
  resetDerivativeEvaluationBase();
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::updateCalculationResult
// Purpose       : Updates the vector that holds the measure values.  This
//                 vector may hold multiple values if the RISE, FALL or CROSS
//                 value is <0.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/20/2020
//-----------------------------------------------------------------------------
void DerivativeEvaluation::updateCalculationResult(double currentIndepVarVal)
{
  double val = getDerivativeValue(currentIndepVarVal);
  if (RFC_ >= 0)
  {
    // store the value once the requested rise (or fall or cross) number has been found
    // has been found
    if ( (riseGiven_ && (actualRise_>= rise_)) || (fallGiven_ && (actualFall_>= fall_)) ||
         (crossGiven_ && (actualCross_>= cross_)) )
    {
      calculationResultVec_.push_back(val);
      calculationResult_ = val;
    }
  }
  else
  {
    // For negative values, store at most the requested number of values (RFC_). If the
    // size of this vector is less than abs(RFC_) then the measure is "failed".  Otherwise
    // the current measure value is in calculationResultVec_[0].
    calculationResultVec_.push_back(val);
    if (calculationResultVec_.size() > abs(RFC_))
      calculationResultVec_.erase(calculationResultVec_.begin());

    if (calculationResultVec_.size() == abs(RFC_))
      calculationResult_ = calculationResultVec_[0];
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::updateCalculationInstant
// Purpose       : Updates the vector that holds the times (or frequencies or
//                 DC sweep values) when the measure was satisified.  This
//                 vector may hold multiple values if the RISE, FALL or CROSS
//                 value is <0.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/20/2020
//-----------------------------------------------------------------------------
void DerivativeEvaluation::updateCalculationInstant(double val)
{
  if (RFC_ >= 0)
  {
    // store the value once the requested rise (or fall or cross) number has been found
    // has been found
    if ( (riseGiven_ && (actualRise_>= rise_)) || (fallGiven_ && (actualFall_>= fall_)) ||
         (crossGiven_ && (actualCross_>= cross_)) )
    {
      calculationInstantVec_.push_back(val);
      calculationInstant_ = val;
    }
  }
  else
  {
    // For negative values, store at most the requested number of values (RFC_). If the
    // size of this vector is less than abs(RFC_) then the measure is "failed".  Otherwise
    // the current measure instant is in calculationInstantVec_[0].
    calculationInstantVec_.push_back(val);
    if (calculationInstantVec_.size() > abs(RFC_))
      calculationInstantVec_.erase(calculationInstantVec_.begin());

    if (calculationInstantVec_.size() == abs(RFC_))
      calculationInstant_ = calculationInstantVec_[0];
  }

  return;
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

  if (atGiven_ && resultFound_)
  {
    os << name_ << " = " << this->getMeasureResult() << std::endl;
  }
  else if (whenGiven_ && (((RFC_ >= 0) && (calculationResultVec_.size() > 0)) ||
			  ((RFC_ < 0) && (calculationResultVec_.size() == abs(RFC_)))) )
  {
    // For non-negative RFC values, the calculationResultVec_ will be non-empty for a successful
    // FIND-WHEN or WHEN meaure.  For negative RFC values, the calculationResultVec_ will
    // have the correct number of values in it.
    os << name_ << " = " << this->getMeasureResult() << std::endl;
  }
  else
  {
    if (measureMgr_.isMeasFailGiven() && measureMgr_.getMeasFail() )
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

  if (atGiven_ && resultFound_)
  {
    os << name_ << " = " << this->getMeasureResult() << " for AT = " << at_;
  }
  else if (whenGiven_ && (((RFC_ >= 0) && (calculationResultVec_.size() > 0)) ||
			  ((RFC_ < 0) && (calculationResultVec_.size() == abs(RFC_)))) )
  {
    // modeStr is "time" for TRAN or TRAN_CONT mode, "freq" for AC or AC_CONT mode
    // and "<sweep variable> value" for DC or DC_CONT mode.
    std::string modeStr = setModeStringForMeasureResultText();
    os << name_ << " = " << calculationResultVec_[0]
       << " at " << modeStr << " = " << calculationInstantVec_[0];
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
// Function      : DerivativeEvaluationCont::DerivativeEvaluationCont()
// Purpose       : Constructor
// Special Notes : Continuous version of DERIV measure that can return multiple
//                 measure values if RISE, FALL or CROSS is not specified.
//                 This class is invoked by .MEASURE TRAN_CONT, .MEASURE AC_CONT,
//                 .MEASURE NOISE_CONT and .MEASURE DC_CONT.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/21/2020
//-----------------------------------------------------------------------------
DerivativeEvaluationCont::DerivativeEvaluationCont(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  DerivativeEvaluationBase(measureMgr, measureBlock),
  contCross_(0),
  contRise_(0),
  contFall_(0),
  contRFC_(0)
{
  // these settings will find all times at which the WHEN clause is satisfied,
  // and may return multiple values.
  measureLastRFC_=true;

  if (riseGiven_)
  {
    contRise_=rise_;
    contRFC_=contRise_;
    rise_=-1;
  }
  else if (fallGiven_)
  {
    contFall_=fall_;
    contRFC_=contFall_;
    fall_=-1;
  }
  else if (crossGiven_)
  {
    contCross_=cross_;
    contRFC_=contCross_;
    cross_=-1;
  }
  else
  {
    // default case when RISE, FALL or CROSS is not explicitly given on .MEASURE line
    contCross_=1;
    contRFC_=contCross_;
    cross_=-1;
    crossGiven_=true;
  }
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationCont::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/03/2020
//-----------------------------------------------------------------------------
void DerivativeEvaluationCont::reset()
{
  resetDerivativeEvaluationBase();
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationCont::updateCalculationResult
// Purpose       : Updates the vector that holds the measure values.  This
//                 vector may hold multiple values if the RISE, FALL or CROSS
//                 value is <0.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/20/2020
//-----------------------------------------------------------------------------
void DerivativeEvaluationCont::updateCalculationResult(double currentIndepVarVal)
{
  double val = getDerivativeValue(currentIndepVarVal);

  if (contRFC_ >= 0)
  {
    // store all of the values once the requested number of rises (or falls or crosses)
    // has been found
    if ( (riseGiven_ && (actualRise_>= contRise_)) || (fallGiven_ && (actualFall_>= contFall_)) ||
         (crossGiven_ && (actualCross_>= contCross_)) )
    {
      calculationResultVec_.push_back(val);
      calculationResult_ = val;
    }
  }
  else
  {
    // For negative values, store at most the requested number of values (contRFC_). If the
    // size of this vector is less than abs(contRFC_) then the measure is "failed".  Otherwise
    // the current measure value is in calculationResultVec_[0].
    calculationResultVec_.push_back(val);
    if (calculationResultVec_.size() > abs(contRFC_))
      calculationResultVec_.erase(calculationResultVec_.begin());

    if (calculationResultVec_.size() == abs(contRFC_))
      calculationResult_ = calculationResultVec_[0];
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationCont::updateCalculationInstant
// Purpose       : Updates the vector that holds the times (or frequencies or
//                 DC sweep values) when the measure was satisified.  This
//                 vector may hold multiple values if the RISE, FALL or CROSS
//                 value is <0.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/20/2020
//-----------------------------------------------------------------------------
void DerivativeEvaluationCont::updateCalculationInstant(double val)
{
  if (contRFC_ >= 0)
  {
    // store all of the values once the requested number of rises (or falls or crosses)
    // has been found
    if ( (riseGiven_ && (actualRise_>= contRise_)) || (fallGiven_ && (actualFall_>= contFall_)) ||
         (crossGiven_ && (actualCross_>= contCross_)) )
    {
      calculationInstantVec_.push_back(val);
    }
  }
  else
  {
    // For negative values, store at most the requested number of values (contRFC_). If the
    // size of this vector is less than abs(contRFC_) then the measure is "failed".  Otherwise
    // the current instant value is in calculationInstantVec_[0].
    calculationInstantVec_.push_back(val);
    if (calculationInstantVec_.size() > abs(contRFC_))
      calculationInstantVec_.erase(calculationInstantVec_.begin());
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationCont::printMeasureResult()
// Purpose       : used to print the measurement result to an output stream
//                 object, which is typically the mt0, ma0 or ms0 file
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 8/21/2020
//-----------------------------------------------------------------------------
std::ostream& DerivativeEvaluationCont::printMeasureResult(std::ostream& os)
{
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os << std::scientific << std::setprecision(precision_);

  if (atGiven_ && resultFound_)
  {
    os << name_ << " = " << this->getMeasureResult() << std::endl;
  }
  else if (whenGiven_ && (((contRFC_ >= 0) && (calculationResultVec_.size() > 0)) ||
			  ((contRFC_ < 0) && (calculationResultVec_.size() == abs(contRFC_)))) )
  {
    // For non-negative RFC values, the calculationResultVec_ will be non-empty for a successful
    // FIND-WHEN or WHEN meaure.  For negative RFC values, the calculationResultVec_ will
    // have the correct number of values in it.
    if (contRFC_ >=0)
    {
      for (int i=0; i<calculationResultVec_.size(); i++)
        os << name_ << " = " << calculationResultVec_[i] << std::endl;
    }
    else
    {
      os << name_ << " = " << calculationResultVec_[0] << std::endl;
    }
  }
  else
  {
    if (measureMgr_.isMeasFailGiven() && measureMgr_.getMeasFail() )
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

  return os;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationCont::printVerboseMeasureResult()
// Purpose       : used to print the "verbose" (more descriptive) measurement
//                 result to an output stream object, which is typically stdout
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 8/21/2020
//-----------------------------------------------------------------------------
std::ostream& DerivativeEvaluationCont::printVerboseMeasureResult(std::ostream& os)
{
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os << std::scientific << std::setprecision(precision_);

  if (atGiven_ && resultFound_)
  {
    os << name_ << " = " << this->getMeasureResult() << " for AT = " << at_ << std::endl;
  }
  else if (whenGiven_ && (((contRFC_ >= 0) && (calculationResultVec_.size() > 0)) ||
			  ((contRFC_ < 0) && (calculationResultVec_.size() == abs(contRFC_)))) )
  {
    // modeStr is "time" for TRAN mode, "freq" for AC mode and
    // "<sweep variable> value" for DC mode.
    std::string modeStr = setModeStringForMeasureResultText();
    if (contRFC_ >=0)
    {
      for (int i=0; i<calculationResultVec_.size(); i++)
      {
         os << name_ << " = " << calculationResultVec_[i]
            << " at " << modeStr << " = " << calculationInstantVec_[i] << std::endl;
      }
    }
    else
    {
      os << name_ << " = " << calculationResultVec_[0]
         << " at " << modeStr << " = " << calculationInstantVec_[0] << std::endl;
    }
  }
  else
  {
    os << name_ << " = FAILED";
    if (atGiven_)
    {
      os << " for AT = " << at_;
    }
    os << std::endl;
  }

  return os;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

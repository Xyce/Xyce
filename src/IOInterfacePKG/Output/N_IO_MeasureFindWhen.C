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

#include <N_IO_MeasureFindWhen.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : FindWhenBase::FindWhenBase()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
FindWhenBase::FindWhenBase(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Base(measureMgr, measureBlock),
  lastIndepVarValue_(0.0),
  lastDepVarValue_(0.0),
  lastOutputVarValue_(0.0),
  lastTargValue_(0.0),
  startDCMeasureWindow_(0.0),
  numPointsFound_(0),
  whenIdx_(0)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();

  // Set the references into the outputVarValues_ array.
  // The default values of whenIdx_=0 works for WHEN measures.  For a FIND-WHEN
  // measure, the variable for the FIND clause is in outVarValues_[0], while the
  // variable for the WHEN clause is in outVarValues_[1].
  if (findGiven_ && !atGiven_)
  {
    whenIdx_ = 1;
  }

  // hard code this measure type to ignore the RFC_LEVEL qualifier
  if (rfcLevelGiven_)
  {
    rfcLevelGiven_=false;
    Xyce::Report::UserWarning0() << "RFC_LEVEL will be ignored for measure " << name_;
  }
}

//-----------------------------------------------------------------------------
// Function      : FindWhenBase::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void FindWhenBase::prepareOutputVariables()
{
  // This measurement can involve up to three solution variables.
  // However, if the AT keyword is given then numOutVars should only have one entry
  numOutVars_ = outputVars_.size();

  if ( (numOutVars_ > 1) && atGiven_ )
  {
    std::string msg = "Too many dependent variables for FIND-AT measure, \"" + name_ + "\"";
    Report::UserError0() << msg;
  }

  outVarValues_.resize( numOutVars_, 0.0 );
}

//-----------------------------------------------------------------------------
// Function      : FindWhenBase::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void FindWhenBase::resetFindWhenBase()
{
  resetBase();
  lastIndepVarValue_=0.0;
  lastDepVarValue_=0.0;
  lastOutputVarValue_=0.0;
  lastTargValue_=0.0;
  startDCMeasureWindow_=0.0;
  numPointsFound_=0;
  calculationResultVec_.clear();
  calculationInstantVec_.clear();
}


//-----------------------------------------------------------------------------
// Function      : FindWhenBase::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void FindWhenBase::updateTran(
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

  if( !calculationDone_ && !isInvalidTimeWindow(endSimTime) )
  {
    initialized_ = true;

    if (atGiven_ && withinTimeWindow(at_))
    {
      if (isATcondition(circuitTime))
      {
        calculationResult_= outVarValues_[0] - (circuitTime - at_)*
	        ( (outVarValues_[0] - lastOutputVarValue_)/(circuitTime - lastIndepVarValue_) );
        calculationDone_ = true;
        resultFound_ = true;
      }
    }
    else if (type_ == "WHEN")
    {
       double targVal=updateTargVal();

      // The measure will succeed if the interpolated WHEN time is inside of the
      // measurement window, and the RISE/FALL/CROSS condition is then met.
      if (isWHENcondition(circuitTime, targVal))
      {
        double whenTime;
        if (numPointsFound_ == 1)
	  whenTime = circuitTime;
        else
          whenTime = interpolateCalculationInstant(circuitTime, targVal);

        if (withinTimeWindow(whenTime))
	{
          updateRFCcountForWhen();
          if (withinRFCWindowForWhen())
            updateMeasureVars(circuitTime, targVal, whenTime);
        }
      }
    }
  }

  updateMeasureState(circuitTime);
}

//-----------------------------------------------------------------------------
// Function      : FindWhenBase::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void FindWhenBase::updateDC(
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
    sweepVar_= getDCSweepVarName(dcParamsVec);
    firstSweepValueFound_ = true;
    numPointsFound_++;

    // update our outVarValues_ vector
    updateOutputVars(comm, outVarValues_, dcSweepVal,
      solnVec, stateVec, storeVec, 0, lead_current_vector,
      junction_voltage_vector, lead_current_dqdt_vector, 0, 0, 0, 0);

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

    if( !calculationDone_ && !isInvalidDCsweepWindow(dcParamsVec[0].startVal, dcParamsVec[0].stopVal) )
    {
      initialized_=true;

      if (atGiven_ && withinDCsweepFromToWindow(at_))
      {
        // Process AT qualifer.  The AT value must be within the measurement window.
        if (isATcondition(dcSweepVal))
        {
          calculationResult_= outVarValues_[0] - (dcSweepVal - at_)*
	        ( (outVarValues_[0] - lastOutputVarValue_)/(dcSweepVal - lastIndepVarValue_) );
          calculationDone_ = true;
          resultFound_ = true;
        }
      }
      else if (type_ == "WHEN")
      {
        double targVal=updateTargVal();

        // The measure will succeed if the interpolated WHEN value is inside of the
	// measurement window, and the RISE/FALL/CROSS condition is then met.
        if (isWHENcondition(dcSweepVal, targVal))
        {
          double whenSweepVal;
          if (numPointsFound_ == 1)
	    whenSweepVal = dcSweepVal;
          else
            whenSweepVal = interpolateCalculationInstant(dcSweepVal, targVal);

          if (withinDCsweepFromToWindow(whenSweepVal))
	  {
            updateRFCcountForWhen();
            if (withinRFCWindowForWhen())
              updateMeasureVars(dcSweepVal, targVal, whenSweepVal);
          }
        }
      }
    }

    updateMeasureState(dcSweepVal);
  }
}

//-----------------------------------------------------------------------------
// Function      : FindWhenBase::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/21/2020
//-----------------------------------------------------------------------------
void FindWhenBase::updateAC(
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
  numPointsFound_++;

  // update our outVarValues_ vector
  updateOutputVars(comm, outVarValues_, frequency, solnVec, 0, 0,
                   imaginaryVec, 0, 0, 0, 0, 0, 0, RFparams);

  if (numPointsFound_ == 1)
    setMeasureState(frequency);

  if( !calculationDone_ && !isInvalidFreqWindow(fStart,fStop) )
  {
    initialized_ = true;
    if (atGiven_ && withinFreqWindow(at_))
    {
      // Process AT qualifer.  The AT value must be within the measurement window.
      if (isATcondition(frequency))
      {
        calculationResult_= outVarValues_[0] - (frequency - at_)*
	        ( (outVarValues_[0] - lastOutputVarValue_)/(frequency - lastIndepVarValue_) );
        calculationDone_ = true;
        resultFound_ = true;
      }
    }
    else if (type_ == "WHEN")
    {
        double targVal=updateTargVal();

        // The measure will succeed if the interpolated WHEN frequency is inside of the
	// measurement window, and the RISE/FALL/CROSS condition is then met.
        if (isWHENcondition(frequency, targVal))
        {
          double whenFreq;
          if (numPointsFound_ == 1)
	    whenFreq = frequency;
          else
            whenFreq = interpolateCalculationInstant(frequency, targVal);

          if (withinFreqWindow(whenFreq))
	  {
            updateRFCcountForWhen();
            if (withinRFCWindowForWhen())
	      updateMeasureVars(frequency, targVal, whenFreq);
          }
      }
    }
  }

  updateMeasureState(frequency);
}

//-----------------------------------------------------------------------------
// Function      : FindWhenBase::updateNoise()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 6/03/2020
//-----------------------------------------------------------------------------
void FindWhenBase::updateNoise(
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
  numPointsFound_++;

  // update our outVarValues_ vector
  updateOutputVars(comm, outVarValues_, frequency, solnVec, 0, 0,
                   imaginaryVec, 0, 0, 0,
                   totalOutputNoiseDens, totalInputNoiseDens, noiseDataVec, 0);

  if (numPointsFound_ == 1)
    setMeasureState(frequency);

  if( !calculationDone_ && !isInvalidFreqWindow(fStart,fStop) )
  {
    initialized_ = true;

    if (atGiven_ && withinFreqWindow(at_))
    {
      // Process AT qualifer.  The AT value must be within the measurement window.
      if (isATcondition(frequency))
      {
        calculationResult_= outVarValues_[0] - (frequency - at_)*
	        ( (outVarValues_[0] - lastOutputVarValue_)/(frequency - lastIndepVarValue_) );
        calculationDone_ = true;
        resultFound_ = true;
      }
    }
    else if (type_ == "WHEN")
    {
      {
        double targVal=updateTargVal();

        // The measure will succeed if the interpolated WHEN frequency is inside of the
	// measurement window, and the RISE/FALL/CROSS condition is then met.
        if (isWHENcondition(frequency, targVal))
        {
          double whenFreq;
          if (numPointsFound_ == 1)
	    whenFreq = frequency;
          else
            whenFreq = interpolateCalculationInstant(frequency, targVal);

          if (withinFreqWindow(whenFreq))
	  {
            updateRFCcountForWhen();
            if (withinRFCWindowForWhen())
              updateMeasureVars(frequency, targVal, whenFreq);
          }
        }
      }
    }
  }

  updateMeasureState(frequency);
}

//-----------------------------------------------------------------------------
// Function      : FindWhenBase::isATcondition
// Purpose       : Evaluates if the AT condition is true for all measure modes.
// Special Notes : For AC and NOISE measures, the independent variable is
//                 frequency.  For DC measures, it is the value of the first
//                 variable in the DC sweep vector.  For TRAN, it is circuit
//                 time.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 05/21/2020
//-----------------------------------------------------------------------------
bool FindWhenBase::isATcondition(const double indepVarVal)
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
// Function      : FindWhenBase::isWHENcondition
// Purpose       : Evaluates if the WHEN condition is true for all modes
// Special Notes : For AC and NOISE measures, the independent variable is
//                 frequency.  For DC measures, it is the value of the first
//                 variable in the DC sweep vector.  For TRAN measures, it
//                 is the circuit time.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 05/21/2020
//-----------------------------------------------------------------------------
bool FindWhenBase::isWHENcondition(const double indepVarVal, const double targVal)
{
  bool whenFound=false;

  if (outVarValues_[whenIdx_] == lastDepVarValue_)
  {
    // no cross can occur for a constant signal value.
    return false;
  }
  else if (numPointsFound_ > 1)
  {
    if (fabs(outVarValues_[whenIdx_] - targVal) < minval_)
    {
      // this is the simple case where Xyce output a value within tolerance
      // of the target value
      whenFound=true;
    }
    else
    {
      // check and see if last point and this point bound the target point
      double backDiff    = lastDepVarValue_ - lastTargValue_;
      double forwardDiff = outVarValues_[whenIdx_] - targVal;

      // if we bound the target then either
      //  (backDiff < 0) && (forwardDiff > 0)
      //   OR
      //  (backDiff > 0) && (forwardDiff < 0)
      // or more simply sgn( backDiff ) = - sgn( forwardDiff )
      if( ((backDiff < 0.0) && (forwardDiff > 0.0)) || ((backDiff > 0.0) && (forwardDiff < 0.0)) )
      {
        whenFound = true;
      }
    }
  }

  return whenFound;
}

//-----------------------------------------------------------------------------
// Function      : FindWhenBase::updateMeasureVars()
// Purpose       : Updates the calculation result and calculation instant vectors.
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it
//                 is the value of the first variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/21/2020
//-----------------------------------------------------------------------------
void FindWhenBase::updateMeasureVars(const double currIndepVarVal, const double targVal,
                                     const double whenInstant)
{
  updateCalculationInstant(whenInstant);
  if (findGiven_)
  {
    if (numPointsFound_ == 1)
      updateCalculationResult(outVarValues_[0]);
    else
      updateCalculationResult(interpolateFindValue(currIndepVarVal, targVal, whenInstant));
    }
  else
  {
    updateCalculationResult(whenInstant);
  }

  calculationDone_ = !measureLastRFC_;
  // resultFound_ is used to help control the descriptive output (to stdout) for a FIND-AT
  //  measure.  If it is false, the measure shows FAILED in stdout.
  resultFound_ = true;

  return;
}

//-----------------------------------------------------------------------------
// Function      : FindWhenBase::interpolateCalculationInstant
// Purpose       : Interpolate the time for when the measure is satisifed.
//                 This accounts for case of WHEN V(1)=V(2) where both
//                 variables may be changing.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/04/2020
//-----------------------------------------------------------------------------
double FindWhenBase::interpolateCalculationInstant(double currIndepVarValue, double targVal)
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

  double calcInstant;
  if (a==b && d==c)
  {
    // pathological case of two lines being identical
    calcInstant = currIndepVarValue;
  }
  else
  {
    // This is the algebra for the time, frequency or DC sweep value when the two non-parallel
    // lines associated with WHEN clause intersect.
    calcInstant = (d-c)/(a-b);
  }

  return calcInstant;
}

//-----------------------------------------------------------------------------
// Function      : FindWhenBase::interpolateFindValue
// Purpose       : Interpolate the find value (for a FIND-WHEN measure) based
//                 on the previously determined WHEN time (or frequency or
//                 DC sweep value).
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/18/2020
//-----------------------------------------------------------------------------
double FindWhenBase::interpolateFindValue(double currIndepVarValue, double targVal,
                                          double whenInstant)
{
  double findVal;
  if (fabs(outVarValues_[whenIdx_] - targVal) < minval_)
    findVal = outVarValues_[0];
  else
    findVal = outVarValues_[0] - (currIndepVarValue - whenInstant)*
	     ( (outVarValues_[0] - lastOutputVarValue_)/(currIndepVarValue - lastIndepVarValue_) );

  return findVal;
}

//-----------------------------------------------------------------------------
// Function      : FindWhenBase::updateRFCcountForWhen()
// Purpose       : Updates the rise, fall and cross counts.
// Special Notes : This function is normally only called for a WHEN time that
//                 falls within the FROM-TO window.  So, the WHEN time is always
//                 a valid cross.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/06/2020
//-----------------------------------------------------------------------------
void FindWhenBase::updateRFCcountForWhen()
{
  actualCross_++;
  if (outVarValues_[whenIdx_] > lastDepVarValue_)
    actualRise_++;
  else
    actualFall_++;

  return;
}

//-----------------------------------------------------------------------------
// Function      : FindWhenBase::withinRFCWindowForWhen()
// Purpose       : Determine if the current WHEN time falls within the specified
//                 RFC value.
// Special Notes : Assumes that a WHEN measure without RISE, FALL or CROSS given
//                 on the .MEASURE line defaults to CROSS given.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/06/2020
//-----------------------------------------------------------------------------
bool FindWhenBase::withinRFCWindowForWhen()
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
// Function      : FindWhenBase::setMeasureState()
// Purpose       : initializes the past values of the independent, dependent
//                 and measured variables, as well as the past target level.
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/04/2020
//-----------------------------------------------------------------------------
void FindWhenBase::setMeasureState(const double indepVarVal)
{
  // assigned last dependent and independent var to current value of the independent
  // varible and outVarValue_[whenIdx_].  While we can't interpolate on this step, it
  // ensures that the initial history is something realistic.
  lastOutputValue_ = outVarValues_[0]; // used for RFC counting
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
// Function      : FindWhenBase::updateMeasureState()
// Purpose       : updates the past values of the independent, dependent
//                 and measured variables, as well as the past target level.
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/04/2020
//-----------------------------------------------------------------------------
void FindWhenBase::updateMeasureState(const double indepVarVal)
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
// Function      : FindWhenBase::updateTargVal
// Purpose       : updates the target value for the WHEN clause
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/04/2020
//-----------------------------------------------------------------------------
double FindWhenBase::updateTargVal()
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
// Function      : FindWhenBase::printMeasureWindow
// Purpose       : prints information related to measure window
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 04/8/2020
//-----------------------------------------------------------------------------
std::ostream& FindWhenBase::printMeasureWindow(std::ostream& os, const double endSimTime,
				               const double startSweepVal, const double endSweepVal)
{
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
// Function      : FindWhenBase::checkMeasureLine
// Purpose       : check .MEASURE line for errors that will cause cause dumps
//               : later
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 09/01/2015
//-----------------------------------------------------------------------------
bool FindWhenBase::checkMeasureLine()
{
  bool bsuccess = true;
  // incorrect number of dependent solution variables will cause core dumps in
  // updateTran() function
  if (numDepSolVars_ <= 1 && findGiven_ && !atGiven_)
  {
    // FIND-WHEN
    bsuccess = false;
    Report::UserError0() << name_ << " has incomplete MEASURE line";
  }
  else if (numDepSolVars_ < 1 && !findGiven_)
  {
    // WHEN
    bsuccess = false;
    Report::UserError0() << name_ << " has incomplete MEASURE line";
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : FindWhenBase::printRFCWindow()
// Purpose       : print informaiton about the start time of the RISE or FALL
//                 window, if a valid one was found.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 09/21/2015
//-----------------------------------------------------------------------------
std::ostream& FindWhenBase::printRFCWindow(std::ostream& os)
{
  // no op, for any measure that supports WHEN

  return os;
}

//-----------------------------------------------------------------------------
// Function      : FindWhen::FindWhen()
// Purpose       : Constructor
// Special Notes : Non-continuous version of FindWhen, that returns only
//                 only one measure value.  This class is invoked by
//                 .MEASURE TRAN, .MEASURE AC and .MEASURE DC.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/03/2020
//-----------------------------------------------------------------------------
FindWhen::FindWhen(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  FindWhenBase(measureMgr, measureBlock),
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
// Function      : FindWhen::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/03/2020
//-----------------------------------------------------------------------------
void FindWhen::reset()
{
  resetFindWhenBase();
}

//-----------------------------------------------------------------------------
// Function      : FindWhen::updateCalculationResult
// Purpose       : Updates the vector that holds the measure values.  This
//                 vector may hold multiple values if the RISE, FALL or CROSS
//                 value is <0.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/20/2020
//-----------------------------------------------------------------------------
void FindWhen::updateCalculationResult(double val)
{
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
// Function      : FindWhen::updateCalculationInstant
// Purpose       : Updates the vector that holds the times (or frequencies or
//                 DC sweep values) when the measure was satisified.  This
//                 vector may hold multiple values if the RISE, FALL or CROSS
//                 value is <0.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/20/2020
//-----------------------------------------------------------------------------
void FindWhen::updateCalculationInstant(double val)
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
// Function      : FindWhen::printMeasureResult()
// Purpose       : used to print the measurement result to an output stream
//                 object, which is typically the mt0, ma0 or ms0 file
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 2/22/2015
//-----------------------------------------------------------------------------
std::ostream& FindWhen::printMeasureResult(std::ostream& os)
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
// Function      : FindWhen::printVerboseMeasureResult()
// Purpose       : used to print the "verbose" (more descriptive) measurement
//                 result to an output stream object, which is typically stdout
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 2/22/2015
//-----------------------------------------------------------------------------
std::ostream& FindWhen::printVerboseMeasureResult(std::ostream& os)
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
    os << name_ << " = " << calculationResultVec_[0];
    if (findGiven_)
      os << " at " << modeStr << " = " << calculationInstantVec_[0];
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
// Function      : FindWhenCont::FindWhenCont()
// Purpose       : Constructor
// Special Notes : Continuous version of FindWhen that can return multiple
//                 measure values if RISE, FALL or CROSS is not specified.
//                 This class is invoked by .MEASURE TRAN_CONT, .MEASURE AC_CONT
//                 and .MEASURE DC_CONT.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/03/2020
//-----------------------------------------------------------------------------
FindWhenCont::FindWhenCont(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  FindWhenBase(measureMgr, measureBlock),
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
  }
  else if (fallGiven_)
  {
    contFall_=fall_;
    contRFC_=contFall_;
  }
  else if (crossGiven_)
  {
    contCross_=cross_;
    contRFC_=contCross_;
  }
  else
  {
    //  // default case when RISE, FALL or CROSS is not explicitly given on .MEASURE line
    cross_=0;
    contCross_=cross_;
    contRFC_=contCross_;
    crossGiven_=true;
  }
}

//-----------------------------------------------------------------------------
// Function      : FindWhenCont::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/03/2020
//-----------------------------------------------------------------------------
void FindWhenCont::reset()
{
  resetFindWhenBase();
}

//-----------------------------------------------------------------------------
// Function      : FindWhenCont::updateCalculationResult
// Purpose       : Updates the vector that holds the measure values.  This
//                 vector may hold multiple values if the RISE, FALL or CROSS
//                 value is <0.
// Special Notes : For compatibility with other measure types, the calculationResult_
//                 variable is also updated with the current measure value.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/03/2020
//-----------------------------------------------------------------------------
void FindWhenCont::updateCalculationResult(double val)
{
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
// Function      : FindWhenCont::updateCalculationInstant
// Purpose       : Updates the vector that holds the times (or frequencies or
//                 DC sweep values) when the measure was satisified.  This
//                 vector may hold multiple values if the RISE, FALL or CROSS
//                 value is <0.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/031/2020
//-----------------------------------------------------------------------------
void FindWhenCont::updateCalculationInstant(double val)
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
// Function      : FindWhenCont::printMeasureResult()
// Purpose       : used to print the measurement result to an output stream
//                 object, which is typically the mt0, ma0 or ms0 file
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/03/2020
//-----------------------------------------------------------------------------
std::ostream& FindWhenCont::printMeasureResult(std::ostream& os)
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
// Function      : FindWhenCont::printVerboseMeasureResult()
// Purpose       : used to print the "verbose" (more descriptive) measurement
//                 result to an output stream object, which is typically stdout
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/16/2020
//-----------------------------------------------------------------------------
std::ostream& FindWhenCont::printVerboseMeasureResult(std::ostream& os)
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
         os << name_ << " = " << calculationResultVec_[i];
         if (findGiven_)
           os << " at " << modeStr << " = " << calculationInstantVec_[i];
         os << std::endl;
      }
    }
    else
    {
      os << name_ << " = " << calculationResultVec_[0];
      if (findGiven_)
        os << " at " << modeStr << " = " << calculationInstantVec_[0];
      os << std::endl;
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

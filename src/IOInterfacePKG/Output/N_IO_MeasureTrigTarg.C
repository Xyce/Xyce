//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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

#include <N_IO_MeasureTrigTarg.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::TrigTargBase()
// Purpose       : Base class for TrigTarg measures
// Special Notes : This is the new class, that was implemented for Xyce 7.5
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 9/14/2021
//-----------------------------------------------------------------------------
TrigTargBase::TrigTargBase(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Base(measureMgr, measureBlock),
  trigResult_(0.0),
  targResult_(0.0),
  trigResultFound_(false),
  targResultFound_(false),
  trigRFC_(0),
  targRFC_(0),
  actualTrigRise_(0),
  actualTrigFall_(0),
  actualTrigCross_(0),
  measureTrigLastRFC_(false),
  actualTargRise_(0),
  actualTargFall_(0),
  actualTargCross_(0),
  measureTargLastRFC_(0),
  startDCMeasureWindow_(0.0),
  numPointsFound_(0),
  targIdx_(1),
  lastIndepVarValue_(0.0),
  lastTrigOutputValue_(0.0),
  lastTargOutputValue_(0.0),
  lastTrigValue_(0.0),
  lastTargValue_(0.0),
  trigCalculationDone_(false),
  targCalculationDone_(false)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();

  // Set the references into the outputVarValues_ array.
  // The default values of targIdx_=1 works if the TRIG clause uses the v(a)=<val> syntax.
  // If the TRIG clause uses the AT=<time> syntax then the variable for the TARG clause is 
  // in outVarValues_[0].  If the TRIG clause uses the v(a)=v(b) syntax then the variable 
  // for the TARG clause is in outVarValues_[2].
  if (trigATgiven_) 
  { 
    // this handles the AT=<time> syntax in the TRIG clause
    targIdx_ = 0;
  } 
  else if ( !trigATgiven_ && !trigOutputValueTargetGiven_ )
  {
    // this handles the v(a)=v(b) syntax in the TRIG clause
    targIdx_ = 2;
  }

  if (trigRiseGiven_)
    trigRFC_=trigRise_;
  else if (trigFallGiven_)
    trigRFC_=trigFall_;
  else if (trigCrossGiven_)
    trigRFC_=trigCross_;
  else
  {
    // default case when RISE, FALL or CROSS is not explicitly given on .MEASURE line
    trigCrossGiven_=true;
    trigCross_=0;
    trigRFC_=0;
  }

  if (targRiseGiven_)
    targRFC_=targRise_;
  else if (targFallGiven_)
    targRFC_=targFall_;
  else if (targCrossGiven_)
    targRFC_=trigCross_;
  else
  {
    // default case when RISE, FALL or CROSS is not explicitly given on .MEASURE line
    targCrossGiven_=true;
    targCross_=0;
    targRFC_=0;
  }

  // If TD is only given for the TRIG clause then use that TD value for
  // both TRIG and TARG.  The converse is not true.
  if (trigTDgiven_ && !targTDgiven_)
  {
    targTD_ = trigTD_;
    targTDgiven_ = true;
  }

  // enforce precedence of AT over TD, for each clause.  This must come after the
  // preceeding conditional about using the TRIG TD for both, if only TRIG TD is given.
  if (trigATgiven_)
    trigTDgiven_ = false;
  if (targATgiven_)
    targTDgiven_ = false;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void TrigTargBase::prepareOutputVariables()
{
  // This measurement can involve up to three solution variables.
  // However, if the AT keyword is given then numOutVars should only have one entry
  numOutVars_ = outputVars_.size();
  outVarValues_.resize( numOutVars_, 0.0 );
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : protected
// Creator       : Pete Sholander, SNL
// Creation Date : 9/14/2021
//-----------------------------------------------------------------------------
void TrigTargBase::resetTrigTargBase()
{
  resetBase();
  lastIndepVarValue_=0.0;
  lastIndepVarValue_=0.0;
  lastTrigOutputValue_=0.0;
  lastTargOutputValue_=0.0;
  lastTrigValue_=0.0;
  lastTargValue_=0.0;
  trigResult_=0.0;
  targResult_=0.0;
  trigResultFound_=0.0;
  targResultFound_=false;
  trigCalculationDone_=false;
  targCalculationDone_=false;
  actualTrigRise_=0;
  actualTrigFall_=0;
  actualTrigCross_=0;
  actualTargRise_=0;
  actualTargFall_=0;
  actualTargCross_=0;
  startDCMeasureWindow_=0.0;
  numPointsFound_=0;
  trigResultVec_.clear();
  targResultVec_.clear();
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 9/14/2021
//-----------------------------------------------------------------------------
void TrigTargBase::updateTran(
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
  numPointsFound_++;

  // update our outVarValues_ vector
  updateOutputVars(comm, outVarValues_, circuitTime,
    solnVec, stateVec, storeVec, 0, lead_current_vector,
    junction_voltage_vector, lead_current_dqdt_vector, 0, 0, 0, 0);

  if (numPointsFound_ == 1)
    setMeasureState(circuitTime);

  // handle invalid AT times
  if (trigATgiven_ && isInvalidAT(trigAT_, 0, endSimTime))
    trigCalculationDone_ = true;
  if (targATgiven_ && isInvalidAT(targAT_, 0, endSimTime))
    targCalculationDone_ = true;

  if( !trigCalculationDone_ && !(trigTDgiven_&& (trigTD_ > endSimTime)) )
  {
    initialized_ = true;
    if ( trigATgiven_ )
    {
      if (circuitTime - minval_ >= trigAT_)
      {
        // Process AT qualifer.  The AT value must be within the measurement window.
        updateTrigVarsForAT();
      }
    }
    else
    {
      double trigTargetValue = updateTrigTargetVal();

      if ( isWHENcondition(circuitTime, outVarValues_[0], lastTrigOutputValue_, trigTargetValue, lastTrigValue_) )
      {
        double whenTime;
        if (numPointsFound_ == 1)
	  whenTime = circuitTime;
        else
	  whenTime = interpolateCalculationInstant(circuitTime, outVarValues_[0], lastTrigOutputValue_,
                                                   trigTargetValue, lastTrigValue_);

        if (withinTrigTDwindow(whenTime))
	{
          updateTrigRFCcount();
          if (withinTrigRFCWindow())
            updateTrigVarsForWhen(whenTime);
        }
      }

      lastTrigValue_ = trigTargetValue;
    }
  }

  if( !targCalculationDone_ && !(targTDgiven_&& (targTD_ > endSimTime)) )
  {
    initialized_ = true;
    if ( targATgiven_)
    {
      if (circuitTime - minval_ >= targAT_)
      {
        // Process AT qualifer.  The AT value must be within the measurement window.
        updateTargVarsForAT();
      }
    }
    else
    {
      double targTargetValue = updateTargTargetVal();

      if ( isWHENcondition(circuitTime, outVarValues_[targIdx_], lastTargOutputValue_, targTargetValue, lastTargValue_) )
      {
        double whenTime;
        if (numPointsFound_ == 1)
	   whenTime = circuitTime;
        else
	  whenTime = interpolateCalculationInstant(circuitTime, outVarValues_[targIdx_], lastTargOutputValue_,
                                                   targTargetValue, lastTargValue_);

        if (withinTargTDwindow(whenTime))
	{
          updateTargRFCcount();
          if (withinTargRFCWindow())
            updateTargVarsForWhen(whenTime);
        }
      }

      lastTargValue_ = targTargetValue;
    }
  }

  updateMeasureState(circuitTime);
}  

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 9/14/2021
//-----------------------------------------------------------------------------
void TrigTargBase::updateDC(
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

    // handle invalid AT sweep values
    if (trigATgiven_ && isInvalidAT(trigAT_, dcParamsVec[0].startVal, dcParamsVec[0].stopVal))
      trigCalculationDone_ = true;
    if (targATgiven_ && isInvalidAT(targAT_, dcParamsVec[0].startVal, dcParamsVec[0].stopVal))
      targCalculationDone_ = true;

    if( !trigCalculationDone_ && !(trigTDgiven_&& ( (dcSweepAscending_ && trigTD_ > dcParamsVec[0].stopVal) ||
						     (!dcSweepAscending_ && trigTD_ > dcParamsVec[0].startVal))) )
    {
      initialized_ = true;
      if ( trigATgiven_)
      {
        if ( (dcSweepAscending_ && (dcSweepVal - minval_ >= trigAT_)) ||
             (!dcSweepAscending_ && (dcSweepVal - minval_ <= trigAT_)) )
        {
          // Process AT qualifer.  The AT value must be within the measurement window.
          updateTrigVarsForAT();
        }
      }
      else
      {
        double trigTargetValue = updateTrigTargetVal();

        if ( isWHENcondition(dcSweepVal, outVarValues_[0], lastTrigOutputValue_, trigTargetValue, lastTrigValue_) )
        {
          double whenSweepVal;
          if (numPointsFound_ == 1)
	    whenSweepVal = dcSweepVal;
          else
	    whenSweepVal = interpolateCalculationInstant(dcSweepVal, outVarValues_[0], lastTrigOutputValue_,
                                                   trigTargetValue, lastTrigValue_);

          if (withinTrigTDwindow(whenSweepVal))
	  {
            updateTrigRFCcount();
            if (withinTrigRFCWindow())
              updateTrigVarsForWhen(whenSweepVal);
          }
        }

        lastTrigValue_ = trigTargetValue;
      }
    }

    if( !targCalculationDone_ && !( targTDgiven_ && ( (dcSweepAscending_ && targTD_ > dcParamsVec[0].stopVal) ||
				                      (!dcSweepAscending_ && targTD_ > dcParamsVec[0].startVal))) )
    {
      initialized_ = true;
      if ( targATgiven_)
      {
        if ( (dcSweepAscending_ && (dcSweepVal - minval_ >= targAT_)) ||
	     (!dcSweepAscending_ && (dcSweepVal - minval_ <= targAT_)) )
        {
          // Process AT qualifer.  The AT value must be within the measurement window.
          updateTargVarsForAT();
        }
      }
      else
      {
        double targTargetValue = updateTargTargetVal();

        if ( isWHENcondition(dcSweepVal, outVarValues_[targIdx_], lastTargOutputValue_, targTargetValue, lastTargValue_) )
        {
          double whenSweepVal;
          if (numPointsFound_ == 1)
	    whenSweepVal = dcSweepVal;
          else
	    whenSweepVal = interpolateCalculationInstant(dcSweepVal, outVarValues_[targIdx_], lastTargOutputValue_,
                                                   targTargetValue, lastTargValue_);

          if (withinTargTDwindow(whenSweepVal))
	  {
            updateTargRFCcount();
            if (withinTargRFCWindow())
              updateTargVarsForWhen(whenSweepVal);
          }
        }

        lastTargValue_ = targTargetValue;
      }
    }

    updateMeasureState(dcSweepVal);
  }
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 9/14/2021
//-----------------------------------------------------------------------------
void TrigTargBase::updateAC(
  Parallel::Machine comm,
  double frequency,
  double fStart,
  double fStop,
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

  // handle invalid AT frequencies
  if (trigATgiven_ && isInvalidAT(trigAT_, fStart, fStop))
    trigCalculationDone_ = true;
  if (targATgiven_ && isInvalidAT(targAT_, fStart, fStop))
    targCalculationDone_ = true;

  if( !trigCalculationDone_ && !(trigTDgiven_&& (trigTD_ > fStop)) )
  {
    initialized_ = true;
    if ( trigATgiven_)
    {
      if (frequency - minval_ >= trigAT_)
      {
        // Process AT qualifer.  The AT value must be within the measurement window.
        updateTrigVarsForAT();
      }
    }
    else
    {
      double trigTargetValue = updateTrigTargetVal();

      if ( isWHENcondition(frequency, outVarValues_[0], lastTrigOutputValue_, trigTargetValue, lastTrigValue_) )
      {
        double whenFreq;
        if (numPointsFound_ == 1)
	  whenFreq = frequency;
        else
	  whenFreq = interpolateCalculationInstant(frequency, outVarValues_[0], lastTrigOutputValue_,
                                                   trigTargetValue, lastTrigValue_);

        if (withinTrigTDwindow(whenFreq))
	{
          updateTrigRFCcount();
          if (withinTrigRFCWindow())
            updateTrigVarsForWhen(whenFreq);
        }
      }

      lastTrigValue_ = trigTargetValue;
    }
  }

  if( !targCalculationDone_ && !(targTDgiven_&& (targTD_ > fStop)) )
  {
    initialized_ = true;
    if ( targATgiven_)
    {
      if (frequency - minval_ >= targAT_)
      {
        // Process AT qualifer.  The AT value must be within the measurement window.
        updateTargVarsForAT();
      }
    }
    else
    {
      double targTargetValue = updateTargTargetVal();

      if ( isWHENcondition(frequency, outVarValues_[targIdx_], lastTargOutputValue_, targTargetValue, lastTargValue_) )
      {
        double whenFreq;
        if (numPointsFound_ == 1)
	   whenFreq = frequency;
        else
	  whenFreq = interpolateCalculationInstant(frequency, outVarValues_[targIdx_], lastTargOutputValue_,
                                                   targTargetValue, lastTargValue_);

        if (withinTargTDwindow(whenFreq))
	{
          updateTargRFCcount();
          if (withinTargRFCWindow())
            updateTargVarsForWhen(whenFreq);
        }
      }

      lastTargValue_ = targTargetValue;
    }
  }

  updateMeasureState(frequency);
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::updateNoise()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 9/14/2021
//-----------------------------------------------------------------------------
void TrigTargBase::updateNoise(
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

  // Used in descriptive output to stdout
  firstSweepValueFound_ = true;
  numPointsFound_++;

  // update our outVarValues_ vector
  updateOutputVars(comm, outVarValues_, frequency, solnVec, 0, 0,
                   imaginaryVec, 0, 0, 0,
                   totalOutputNoiseDens, totalInputNoiseDens, noiseDataVec, 0);

  if (numPointsFound_ == 1)
    setMeasureState(frequency);

  // handle invalid AT frequencies
  if (trigATgiven_ && isInvalidAT(trigAT_, fStart, fStop))
    trigCalculationDone_ = true;
  if (targATgiven_ && isInvalidAT(targAT_, fStart, fStop))
    targCalculationDone_ = true;

  if( !trigCalculationDone_ && !(trigTDgiven_&& (trigTD_ > fStop)) )
  {
    initialized_ = true;
    if ( trigATgiven_)
    {
      if (frequency - minval_ >= trigAT_)
      {
        // Process AT qualifer.  The AT value must be within the measurement window.
        updateTrigVarsForAT();
      }
    }
    else
    {
      double trigTargetValue = updateTrigTargetVal();

      if ( isWHENcondition(frequency, outVarValues_[0], lastTrigOutputValue_, trigTargetValue, lastTrigValue_) )
      {
        double whenFreq;
        if (numPointsFound_ == 1)
	  whenFreq = frequency;
        else
	  whenFreq = interpolateCalculationInstant(frequency, outVarValues_[0], lastTrigOutputValue_,
                                                   trigTargetValue, lastTrigValue_);

        if (withinTrigTDwindow(whenFreq))
	{
          updateTrigRFCcount();
          if (withinTrigRFCWindow())
            updateTrigVarsForWhen(whenFreq);
        }
      }

      lastTrigValue_ = trigTargetValue;
    }
  }

  if( !targCalculationDone_ && !(targTDgiven_&& (targTD_ > fStop)) )
  {
    initialized_ = true;
    if ( targATgiven_)
    {
      if (frequency - minval_ >= targAT_)
      {
        // Process AT qualifer.  The AT value must be within the measurement window.
        updateTargVarsForAT();
      }
    }
    else
    {
      double targTargetValue = updateTargTargetVal();

      if ( isWHENcondition(frequency, outVarValues_[targIdx_], lastTargOutputValue_, targTargetValue, lastTargValue_) )
      {
        double whenFreq;
        if (numPointsFound_ == 1)
	   whenFreq = frequency;
        else
	  whenFreq = interpolateCalculationInstant(frequency, outVarValues_[targIdx_], lastTargOutputValue_,
                                                   targTargetValue, lastTargValue_);

        if (withinTargTDwindow(whenFreq))
	{
          updateTargRFCcount();
          if (withinTargRFCWindow())
            updateTargVarsForWhen(whenFreq);
        }
      }

      lastTargValue_ = targTargetValue;
    }
  }

  updateMeasureState(frequency);
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::updateTrigVarsForAT
// Purpose       : Updates the variable associated with the TRIG calculation
// Special Notes : For AC and NOISE measures, the independent variable is
//                 frequency.  For DC measures, it is the value of the first
//                 variable in the DC sweep vector.  For TRAN, it is circuit
//                 time.
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 09/15/2021
//-----------------------------------------------------------------------------
void TrigTargBase::updateTrigVarsForAT()
{
  trigResult_ = trigAT_;
  trigCalculationDone_ = true;
  trigResultFound_ = true;
  
  return;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::updateTargVarsForAT
// Purpose       : Updates the variable associated with the TARG calculation
// Special Notes : For AC and NOISE measures, the independent variable is
//                 frequency.  For DC measures, it is the value of the first
//                 variable in the DC sweep vector.  For TRAN, it is circuit
//                 time.
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 09/15/2021
//-----------------------------------------------------------------------------
void TrigTargBase::updateTargVarsForAT()
{
  targResult_ = targAT_;
  targCalculationDone_ = true;
  targResultFound_ = true;
  
  return;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::isWHENcondition
// Purpose       : Evaluates if the WHEN condition is true for all modes
// Special Notes : For AC and NOISE measures, the independent variable is
//                 frequency.  For DC measures, it is the value of the first
//                 variable in the DC sweep vector.  For TRAN measures, it
//                 is the circuit time.
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
bool TrigTargBase::isWHENcondition(
  double indepVarValue,
  double depVarValue,
  double lastDepVarValue, 
  double targetValue,
  double lastTargetValue) const
{
  bool whenFound=false;

  if (depVarValue == lastDepVarValue)
  {
    // no cross can occur for a constant signal value.
    return false;
  }
  else if (numPointsFound_ > 1)
  {
    if (fabs(depVarValue - targetValue) < minval_)
    {
      // this is the simple case where Xyce output a value within tolerance
      // of the target value
      whenFound=true;
    }
    else
    {
      // check and see if last point and this point bound the target point
      double backDiff    = lastDepVarValue - lastTargetValue;
      double forwardDiff = depVarValue - targetValue;

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
// Function      : TrigTargBase::updateTrigVarsForWhen()
// Purpose       : Updates the variables associated with the TRIG calculations.
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it
//                 is the value of the first variable in the DC sweep vector.
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
void TrigTargBase::updateTrigVarsForWhen(double whenInstant)
{
  updateTrigResult(whenInstant);
  trigCalculationDone_ = !measureTrigLastRFC_;

  // trigResultFound_ is used to help control the descriptive output (to stdout)
  trigResultFound_ = true;

  return;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::updateTargVarsForWhen()
// Purpose       : Updates the variables associated with the TARG calculations.
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it
//                 is the value of the first variable in the DC sweep vector.
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
void TrigTargBase::updateTargVarsForWhen(double whenInstant)
{
  updateTargResult(whenInstant);
  targCalculationDone_ = !measureTargLastRFC_;

  // targResultFound_ is used to help control the descriptive output (to stdout)
  targResultFound_ = true;

  return;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::interpolateCalculationInstant
// Purpose       : Interpolate the time for when the TRIG or TARG condition
//                 measure is satisifed.
//                 This accounts for case of TRIG V(1)=V(2) where both
//                 variables may be changing.
// Special Notes :
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
double TrigTargBase::interpolateCalculationInstant(
  double currIndepVarValue,
  double currDepVarValue,
  double lastDepVarValue,
  double targetVal,
  double lastTargetVal) const
{
  // Calculate slopes and y-intercepts of the two lines, to get them into
  // canonical y=ax+c and y=bx+d form.  If the WHEN clause is of the form
  // WHEN V(1)=V(2) then the line with (a,c) is the value of V(1), which is the
  // "dependent variable".  The line with (b,d) is the value of V(2), which
  // is the "target value".
  double a = (currDepVarValue - lastDepVarValue)/(currIndepVarValue - lastIndepVarValue_);
  double b = (targetVal - lastTargetVal)/(currIndepVarValue - lastIndepVarValue_);
  double c = currDepVarValue - a*currIndepVarValue;
  double d = targetVal - b*currIndepVarValue;

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
// Function      : TrigTargBase::isInvalidAT
// Purpose       : Is the specified AT value outside of the simulation bounds         
// Special Notes : This function only works for AC, NOISE and TRAN.
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 10/12/2021
//-----------------------------------------------------------------------------
bool TrigTargBase::isInvalidAT(double atVal, double startVal, double stopVal) const
{
  if ( (fabs(atVal - startVal) < fabs(to_*minval_)) || (fabs(atVal - stopVal) < fabs(to_*minval_)) ) 
    return false;
  else if ( (startVal <= stopVal) && (atVal < startVal || atVal > stopVal) )
    return true;
  else if  ( (startVal > stopVal) && (atVal > startVal || atVal < stopVal) )
    return true;

  return false;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::updateTrigRFCcount()
// Purpose       : Updates the rise, fall and cross counts for the TRIG clause.
// Special Notes : This function is normally only called for a TRIG time that
//                 falls within the TRIG TD window.  So, the TRIG time is always
//                 a valid cross.
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
void TrigTargBase::updateTrigRFCcount()
{
  actualTrigCross_++;
  if (outVarValues_[0] > lastTrigOutputValue_)
    actualTrigRise_++;
  else
    actualTrigFall_++;

  return;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::updateTargRFCcount()
// Purpose       : Updates the rise, fall and cross counts for the TARG clause.
// Special Notes : This function is normally only called for a TARG time that
//                 falls within the TARG TD window
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
void TrigTargBase::updateTargRFCcount()
{
  actualTargCross_++;
  if (outVarValues_[targIdx_] > lastTargOutputValue_)
    actualTargRise_++;
  else
    actualTargFall_++;

  return;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::withinTrigRFCWindow()
// Purpose       : Determine if the current trig time falls within the specified
//                 RFC value.
// Special Notes : Assumes that a TRIG measure without RISE, FALL or CROSS given
//                 on the .MEASURE line defaults to CROSS given.
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 9/14/2021
//-----------------------------------------------------------------------------
bool TrigTargBase::withinTrigRFCWindow() const
{
  bool retVal=false;

  if (trigRiseGiven_ && (outVarValues_[0] > lastTrigOutputValue_) && (actualTrigRise_ >= trigRise_))
    retVal=true;
  else if (trigFallGiven_ && (outVarValues_[0] < lastTrigOutputValue_) && (actualTrigFall_ >= trigFall_))
    retVal=true;
  else if ( !(trigRiseGiven_ || trigFallGiven_) && (actualTrigCross_ >= trigCross_) )
    retVal=true;

  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::withinTargRFCWindow()
// Purpose       : Determine if the current targ time falls within the specified
//                 RFC value.
// Special Notes : Assumes that a TARG measure without RISE, FALL or CROSS given
//                 on the .MEASURE line defaults to CROSS given.
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
bool TrigTargBase::withinTargRFCWindow() const
{
  bool retVal=false;

  if (targRiseGiven_ && (outVarValues_[targIdx_] > lastTargOutputValue_) && (actualTargRise_ >= targRise_))
    retVal=true;
  else if (targFallGiven_ && (outVarValues_[targIdx_] < lastTargOutputValue_) && (actualTargFall_ >= targFall_))
    retVal=true;
  else if ( !(targRiseGiven_ || targFallGiven_) && (actualTargCross_ >= targCross_) )
    retVal=true;

  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::withinTrigTDwindow
// Purpose       : Checks if the independent variable is within the TRIG TD window 
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.   For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 10/12/2021
//-----------------------------------------------------------------------------
bool TrigTargBase::withinTrigTDwindow(double indepVarVal) const
{
  return !trigTDgiven_ || (trigTDgiven_ && (indepVarVal > trigTD_*(1-minval_)));
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::withinTargTDwindow
// Purpose       : Checks if the independent variable is within the TARG TD window 
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.   For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 10/12/2021
//-----------------------------------------------------------------------------
bool TrigTargBase::withinTargTDwindow(double indepVarVal) const
{
  return !targTDgiven_ || (targTDgiven_ && (indepVarVal > targTD_*(1-minval_)));
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::setMeasureState()
// Purpose       : initializes the past values of the independent, dependent
//                 and measured variables, as well as the past target level.
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
void TrigTargBase::setMeasureState(double indepVarVal)
{
  // assigned last dependent and independent var to current value of the independent
  // varible and outVarValue_[whenIdx_].  While we can't interpolate on this step, it
  // ensures that the initial history is something realistic.
  lastIndepVarValue_ = indepVarVal;
  if (!trigATgiven_)
    lastTrigOutputValue_ = outVarValues_[0];
  if (!targATgiven_)
    lastTargOutputValue_=outVarValues_[targIdx_];

  return;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::updateMeasureState()
// Purpose       : updates the past values of the independent, dependent
//                 and measured variables, as well as the past target level.
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
void TrigTargBase::updateMeasureState(double indepVarVal)
{
  lastIndepVarValue_ = indepVarVal;
  if (!trigATgiven_)
    lastTrigOutputValue_ = outVarValues_[0];
  if (!targATgiven_)
    lastTargOutputValue_ = outVarValues_[targIdx_];

  return;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::updateTrigTargetVal_
// Purpose       : updates the target value for the TRIG clause
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
double TrigTargBase::updateTrigTargetVal()
{
  double trigTargetVal = 0.0;

  if( trigOutputValueTargetGiven_ )
  {
    // This is the form TRIG v(a)=fixed value
    trigTargetVal = trigOutputValueTarget_;
  }
  else
  {
    // This is the form TRIG v(a)= potentially changing value, such as v(a)=v(b)
    // in that case v(b) is in outVarValues_[1]
    trigTargetVal = outVarValues_[1];
  }

  return trigTargetVal;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::updateTargTargetVal_
// Purpose       : updates the target value for the TRIG clause
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
double TrigTargBase::updateTargTargetVal()
{
  double targTargetVal = 0.0;

  if( targOutputValueTargetGiven_ )
  {
    // This is the form TARG v(a)=fixed value
    targTargetVal = targOutputValueTarget_;
  }
  else
  {
    // This is the form TARG v(a)= potentially changing value, such as v(a)=v(b)
    // in that case v(b) is in outVarValues_[targIdx+1]
    targTargetVal = outVarValues_[targIdx_+1];
  }

  return targTargetVal;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 9/14/2021
//-----------------------------------------------------------------------------
double TrigTargBase::getMeasureResult()
{
  // if both the TRIG and TARG times were not found then the measure value is default value (e.g., -1)
  if ( trigResultFound_ && targResultFound_ )
    calculationResult_ = targResult_ - trigResult_;

  return calculationResult_;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::printMeasureWindow()
// Purpose       : print information about the measurement window
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
std::ostream& TrigTargBase::printMeasureWindow(std::ostream& os, double endSimTime,
				               double startSweepVal, double endSweepVal) const
{
  // no op, for this measure type

  return os;
}
//-----------------------------------------------------------------------------
// Function      : TrigTargBase::printRFCWindow()
// Purpose       : print information about the start time of the RISE or FALL
//                 window, if a valid one was found.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
std::ostream& TrigTargBase::printRFCWindow(std::ostream& os) const
{
  // no op, for this measure type

  return os;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::checkMeasureLine
// Purpose       : check .MEASURE line for errors that will cause cause
//                 segfaults later
// Special Notes :
// Scope         : protected
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
bool TrigTargBase::checkMeasureLine() const
{
  bool bsuccess = true;
  int expectedNumDepSolVars=0;

  if ( !trigATgiven_ && trigOutputValueTargetGiven_ )
    expectedNumDepSolVars++;
  else if ( !trigATgiven_ && !trigOutputValueTargetGiven_ )
    expectedNumDepSolVars +=2;

  if ( !targATgiven_ && targOutputValueTargetGiven_ )
    expectedNumDepSolVars++;
  else if ( !targATgiven_ && !targOutputValueTargetGiven_ )
    expectedNumDepSolVars +=2;

  if ( ((numDepSolVars_ == 0) && !(trigATgiven_ && targATgiven_)) || (numDepSolVars_ > 4) ||
       (expectedNumDepSolVars != numDepSolVars_) )
       
  {
    bsuccess = false;
    Report::UserError0() << name_ << " has invalid MEASURE line";
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::printMeasureWarningsForAT()
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
void TrigTargBase::printMeasureWarningsForAT(double endSimTime) const
{
  // no op, for this measure type
}

//-----------------------------------------------------------------------------
// Function      : TrigTargBase::printMeasureWarnings()
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
void TrigTargBase::printMeasureWarnings(double endSimTime, double startSweepVal,
                                        double endSweepVal) const
{
  // no op, for this measure type
}


//-----------------------------------------------------------------------------
// Function      : TrigTarg::TrigTarg()
// Purpose       : Constructor
// Special Notes : Non-continuous version of TrigTarg, that returns only
//                 only one measure value.  This class is invoked by
//                 .MEASURE TRAN, .MEASURE AC and .MEASURE DC.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 9/14/2021
//-----------------------------------------------------------------------------
TrigTarg::TrigTarg(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  TrigTargBase(measureMgr, measureBlock)
{
  // RISE, FALL or CROSS values < -1 are not supported yet
  if ((trigRiseGiven_ && trigRise_ < -1) || (trigFallGiven_ && trigFall_ < -1) || (trigCrossGiven_ && trigCross_ < -1) ||
      (targRiseGiven_ && targRise_ < -1) || (targFallGiven_ && targFall_ < -1) || (targCrossGiven_ && targCross_ < -1))
  {
    Report::UserError0() << " RISE, FALL or CROSS values < -1 not supported for measure " << name_
                         << " for AC, DC, NOISE or TRAN measures";
  }

  measureTrigLastRFC_ = ((trigRiseGiven_ && trigRise_ < 0) || (trigFallGiven_ && trigFall_ < 0) || 
	                 (trigCrossGiven_ && trigCross_ < 0)) ? true : false;

  measureTargLastRFC_ = ((targRiseGiven_ && targRise_ < 0) || (targFallGiven_ && targFall_ < 0) || 
	                 (targCrossGiven_ && targCross_ < 0)) ? true : false;
}

//-----------------------------------------------------------------------------
// Function      : TrigTarg::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 9/14/2021
//-----------------------------------------------------------------------------
void TrigTarg::reset()
{
  resetTrigTargBase();
}

//-----------------------------------------------------------------------------
// Function      : TrigTarg::updateTrigResult
// Purpose       : Updates the vector that holds the measure values.  This
//                 vector may hold multiple values if the RISE, FALL or CROSS
//                 value is < 0.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
void TrigTarg::updateTrigResult(double val)
{
  if (trigRFC_ >= 0)
  {
    // store the value once the requested rise (or fall or cross) number has been found
    // has been found
    if ( (trigRiseGiven_ && (actualTrigRise_>= trigRise_)) || (trigFallGiven_ && (actualTrigFall_>= trigFall_)) ||
         (trigCrossGiven_ && (actualTrigCross_>= trigCross_)) )
    {
      trigResultVec_.push_back(val);
      trigResult_ = val;
    }
  }
  else
  {
    // For negative values, store at most the requested number of values (trigRFC_). If the
    // size of this vector is less than abs(trigRFC_) then the measure is "failed".  Otherwise
    // the current measure value is in trigResultVec_[0].
    trigResultVec_.push_back(val);
    if (trigResultVec_.size() > abs(trigRFC_))
      trigResultVec_.erase(trigResultVec_.begin());

    if (trigResultVec_.size() == abs(trigRFC_))
      trigResult_ = trigResultVec_[0];
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : TrigTarg::updateTargResult
// Purpose       : Updates the vector that holds the measure values.  This
//                 vector may hold multiple values if the RISE, FALL or CROSS
//                 value is < 0.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
void TrigTarg::updateTargResult(double val)
{
  if (targRFC_ >= 0)
  {
    // store the value once the requested rise (or fall or cross) number has been found
    // has been found
    if ( (targRiseGiven_ && (actualTargRise_>= targRise_)) || (targFallGiven_ && (actualTargFall_>= targFall_)) ||
         (targCrossGiven_ && (actualTargCross_>= targCross_)) )
    {
      targResultVec_.push_back(val);
      targResult_ = val;
    }
  }
  else
  {
    // For negative values, store at most the requested number of values (targRFC_). If the
    // size of this vector is less than abs(targRFC_) then the measure is "failed".  Otherwise
    // the current measure value is in targResultVec_[0].
    targResultVec_.push_back(val);
    if (targResultVec_.size() > abs(targRFC_))
      targResultVec_.erase(targResultVec_.begin());

    if (targResultVec_.size() == abs(targRFC_))
      targResult_ = targResultVec_[0];
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : TrigTarg::printMeasureResult()
// Purpose       : used to print the measurement result to an output stream
//                 object, which is typically the mt0, ma0 or ms0 file
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 9/14/2021
//-----------------------------------------------------------------------------
std::ostream& TrigTarg::printMeasureResult(std::ostream& os)
{
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os << std::scientific << std::setprecision(precision_);

  if (trigResultFound_ && targResultFound_)
  {
    os << name_ << " = " << this->getMeasureResult() << std::endl;
  }
  else
  {
    if ( measureMgr_.getMeasFail() )
    {
      // output FAILED to .mt file if .OPTIONS MEASURE MEASFAIL=1 and this is a failed measure.
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
// Function      : TrigTarg::printVerboseMeasureResult()
// Purpose       : used to print the "verbose" (more descriptive) measurement
//                 result to an output stream object, which is typically stdout
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 2/22/2015
//-----------------------------------------------------------------------------
std::ostream& TrigTarg::printVerboseMeasureResult(std::ostream& os)
{
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os << std::scientific << std::setprecision(precision_);

  if (trigResultFound_ && targResultFound_)
  {
    os << name_ << " = " << this->getMeasureResult()
       << " with targ = " << targResult_ << " and trig = " << trigResult_;
  }
  else if (!trigResultFound_ && targResultFound_)
    os << name_ << " = FAILED with targ = "  << targResult_ << " and trig = not found";
  else if (trigResultFound_ && !targResultFound_)
    os << name_ << " = FAILED with targ = not found and trig = " << trigResult_;
  else
    os << name_ << " = FAILED with targ = not found and trig = not found";
  
  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargCont::TrigTargCont()
// Purpose       : Constructor
// Special Notes : Continuous version of TrigTarg that can return multiple
//                 measure values if RISE, FALL or CROSS is not specified.
//                 This class is invoked by .MEASURE TRAN_CONT, .MEASURE AC_CONT
//                 and .MEASURE DC_CONT.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
TrigTargCont::TrigTargCont(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  TrigTargBase(measureMgr, measureBlock)
{
  if ((trigRiseGiven_ && trigRise_ < 0) || (trigFallGiven_ && trigFall_ < 0) || (trigCrossGiven_ && trigCross_ < 0) ||
      (targRiseGiven_ && targRise_ < 0) || (targFallGiven_ && targFall_ < 0) || (targCrossGiven_ && targCross_ < 0))
  {
    Report::UserError0() << " RISE, FALL or CROSS values < 0 not supported for measure " << name_
                         << " for AC_CONT, DC_CONT, NOISE_CONT or TRAN_CONT measures";
  }

  // these settings will find all times at which the WHEN clauses are satisfied,
  // and may return multiple values.
  measureTrigLastRFC_=true;
  measureTargLastRFC_=true;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargCont::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-----------------------------------------------------------------------------
void TrigTargCont::reset()
{
  resetTrigTargBase();
}

//-----------------------------------------------------------------------------
// Function      : TrigTargCont::updateTrigResult
// Purpose       : Updates the vector that holds the measure values.  This
//                 vector may hold multiple values if the RISE, FALL or CROSS
//                 value is < 0.
// Special Notes : For compatibility with other measure types, the calculationResult_
//                 variable is also updated with the current measure value.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 9/14/2021
//-----------------------------------------------------------------------------
void TrigTargCont::updateTrigResult(double val)
{
  if (trigRFC_ >= 0)
  {
    // store all of the values once the requested number of rises (or falls or crosses)
    // has been found
    if ( (trigRiseGiven_ && (actualTrigRise_>= trigRise_)) || (trigFallGiven_ && (actualTrigFall_>= trigFall_)) ||
         (trigCrossGiven_ && (actualTrigCross_>= trigCross_)) )
    {
      trigResultVec_.push_back(val);
      trigResult_ = val;
    }
  }
  else
  {
    // For negative values, store at most the requested number of values (trigRFC_). If the
    // size of this vector is less than abs(trigRFC_) then the measure is "failed".  Otherwise
    // the current measure value is in trigResultVec_[0].
    trigResultVec_.push_back(val);
    if (trigResultVec_.size() > abs(trigRFC_))
      trigResultVec_.erase(trigResultVec_.begin());

    if (trigResultVec_.size() == abs(trigRFC_))
      trigResult_ = trigResultVec_[0];
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargCont::updateTargResult
// Purpose       : Updates the vector that holds the measure values.  This
//                 vector may hold multiple values if the RISE, FALL or CROSS
//                 value is < 0.
// Special Notes : For compatibility with other measure types, the calculationResult_
//                 variable is also updated with the current measure value.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 9/21/2021
//-----------------------------------------------------------------------------
void TrigTargCont::updateTargResult(double val)
{
  if (targRFC_ >= 0)
  {
    // store all of the values once the requested number of rises (or falls or crosses)
    // has been found
    if ( (targRiseGiven_ && (actualTargRise_>= targRise_)) || (targFallGiven_ && (actualTargFall_>= targFall_)) ||
         (targCrossGiven_ && (actualTargCross_>= targCross_)) )
    {
      targResultVec_.push_back(val);
      targResult_ = val;
    }
  }
  else
  {
    // For negative values, store at most the requested number of values (targRFC_). If the
    // size of this vector is less than abs(targRFC_) then the measure is "failed".  Otherwise
    // the current measure value is in targResultVec_[0].
    targResultVec_.push_back(val);
    if (targResultVec_.size() > abs(targRFC_))
      targResultVec_.erase(targResultVec_.begin());

    if (targResultVec_.size() == abs(targRFC_))
      targResult_ = targResultVec_[0];
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargCont::printMeasureResult()
// Purpose       : used to print the measurement result to an output stream
//                 object, which is typically the mt0, ma0 or ms0 file
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 9/14/2021
//-----------------------------------------------------------------------------
std::ostream& TrigTargCont::printMeasureResult(std::ostream& os)
{
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os << std::scientific << std::setprecision(precision_);

  const int colWidth =20;

  if (trigResultFound_ && targResultFound_)
  {
    if (!trigATgiven_ && !targATgiven_)
    {
      const int numTrig = trigResultVec_.size();
      const int numTarg = targResultVec_.size();
      const int numResults = std::min(numTrig, numTarg);

      for (int i=0; i<numResults; i++)
      {
        os << name_ << " = " << targResultVec_[i] - trigResultVec_[i] << std::setw(colWidth) << " "
           << "targ = " << targResultVec_[i] << std::setw(colWidth) << " " << "trig = " << trigResultVec_[i]  << std::endl;
      }
    }
    else if (trigATgiven_ && !targATgiven_)
    {
      os << name_ << " = " << targResultVec_[0] - trigResult_ << std::setw(colWidth) << " "
         << "targ = " << targResultVec_[0] << std::setw(colWidth) << " " << "trig = " << trigResult_  << std::endl;
    }
    else if (!trigATgiven_ && targATgiven_)
    {
      os << name_ << " = " << targResult_ - trigResultVec_[0] << std::setw(colWidth) << " "
         << "targ = " << targResult_ << std::setw(colWidth) << " " << "trig = " << trigResultVec_[0]  << std::endl;
    }
    else
    { 
      os << name_ << " = " << targResult_ - trigResult_ << std::setw(colWidth) << " "
         << "targ = " << targResult_ << std::setw(colWidth) << " "<< "trig = " << trigResult_  << std::endl;
    }
  }
  else
  {
    if ( measureMgr_.getMeasFail() )
    {
      // output FAILED to .mt file if .OPTIONS MEASURE MEASFAIL=1 and this is a failed measure.
      os << name_ << " = FAILED";
    }
    else
    {
      os << name_ << " = " << this->getMeasureResult();
    }

    if (!trigResultFound_ && targResultFound_)
    {
      os << std::setw(colWidth) << " " << "targ = "  << targResultVec_[0]
         << std::setw(colWidth) << " " << "trig = not found" << std::endl;
    }
    else if (trigResultFound_ && !targResultFound_)
    {  
      os << std::setw(colWidth) << " " << "targ = not found" 
         << std::setw(colWidth) << " " << "trig = " << trigResultVec_[0] << std::endl;
    }
    else
    {
      os << std::setw(colWidth) << " " << "targ = not found"
         << std::setw(colWidth) << " " << "trig = not found" << std::endl;
    }
  }

  return os;
}

//-----------------------------------------------------------------------------
// Function      : TrigTargCont::printVerboseMeasureResult()
// Purpose       : used to print the "verbose" (more descriptive) measurement
//                 result to an output stream object, which is typically stdout
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 9/14/2021
//-----------------------------------------------------------------------------
std::ostream& TrigTargCont::printVerboseMeasureResult(std::ostream& os)
{
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os << std::scientific << std::setprecision(precision_);

  if (trigResultFound_ && targResultFound_)
  {
    if (!trigATgiven_ && !targATgiven_)
    {
      const int numTrig = trigResultVec_.size();
      const int numTarg = targResultVec_.size();
      const int numResults = std::min(numTrig, numTarg);

      for (int i=0; i<numResults; i++)
      {
        os << name_ << " = " << targResultVec_[i] - trigResultVec_[i]
           << " with targ = " << targResultVec_[i] << " and trig = " << trigResultVec_[i]  << std::endl;
      }
    }
    else if (trigATgiven_ && !targATgiven_)
    {
      os << name_ << " = " << targResultVec_[0] - trigResult_
         << " with targ = " << targResultVec_[0] << " and trig = " << trigResult_  << std::endl;
    }
    else if (!trigATgiven_ && targATgiven_)
    {
      os << name_ << " = " << targResult_ - trigResultVec_[0]
         << " with targ = " << targResult_ << " and trig = " << trigResultVec_[0]  << std::endl;
    }
    else
    { 
      os << name_ << " = " << targResult_ - trigResult_
         << " with targ = " << targResult_ << " and trig = " << trigResult_  << std::endl;
    }
  }
  else if (!trigResultFound_ && targResultFound_)
  {
    os << name_ << " = FAILED with targ = "  << targResultVec_[0] << " and trig = not found" << std::endl;
  }
  else if (trigResultFound_ && !targResultFound_)
  {  
    os << name_ << " = FAILED with targ = not found and trig = " << trigResultVec_[0] << std::endl;
  }
  else
  {
    os << name_ << " = FAILED with targ = not found and trig = not found" << std::endl;
  }

  return os;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

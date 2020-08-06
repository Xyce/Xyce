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
  doneIfFound_(true),
  lastIndepVarValue_(0.0),
  lastDepVarValue_(0.0),
  lastOutputVarValue_(0.0),
  lastTargValue_(0.0),
  startDCMeasureWindow_(0.0),
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
    // at the first time-step within the measurement window.  (That
    // window is set by the TO-FROM and TD qualifiers if present.)  This is 
    // needed so that the RISE/FALL/CROSS count is not incremented at time=0, if
    // the measured waveform has a DC offset at time=0    
    if (!firstStepInMeasureWindow_)
    {
      lastOutputValue_ = outVarValues_[0]; 
      firstStepInMeasureWindow_ = true;
    }

    if( !initialized_ )
    {
      // assigned last independent and dependent var to current time and outVarValue_[whenIdx_] 
      // while we can't interpolate on this step, it ensures that the initial history is
      // something realistic.
      lastIndepVarValue_=circuitTime;
      lastDepVarValue_=outVarValues_[whenIdx_];
      lastOutputVarValue_=outVarValues_[0];
      if (outputValueTargetGiven_)
        lastTargValue_ = outputValueTarget_;
      else
        lastTargValue_ = outVarValues_[whenIdx_+1];
      initialized_=true;
    }  

    double targVal=0.0;
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
    
    // for the FIND-WHEN and WHEN measures, the rfc level used is either the target value
    // of the WHEN clause, or the value of the RFC_LEVEL qualifier if one is specified.
    if (atGiven_)
    {
      // check and see if last point and this point bound the target point
      double backDiff    = lastIndepVarValue_ - at_;
      double forwardDiff = circuitTime - at_;

      // if we bound the frequency target then either
      //  (backDiff < 0) && (forwardDiff > 0)
      //   OR
      //  (backDiff > 0) && (forwardDiff < 0)
      // or more simply sgn( backDiff ) = - sgn( forwardDiff )
      //
      // Also test for equality, to within the minval_ tolerance, as with the WHEN syntax.
      if ( ((backDiff < 0.0) && (forwardDiff > 0.0)) || ((backDiff > 0.0) && (forwardDiff < 0.0)) ||
	     (((abs(backDiff) < minval_) || (abs(forwardDiff) < minval_))) )
      {
        calculationResult_= outVarValues_[0] - (circuitTime - at_)*
	        ( (outVarValues_[0] - lastOutputVarValue_)/(circuitTime - lastIndepVarValue_) );
        calculationDone_ = true;
        resultFound_ = true;
      }
    }
    else if( (type_ == "WHEN") &&  withinRiseFallCrossWindow(outVarValues_[whenIdx_],targVal) )
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
      // only record the start time of the RFC window for the FIND-WHEN Measure
      if( !firstStepInRfcWindow_  )
      {
        firstStepInRfcWindow_ = true;
        rfcWindowFound_ = true;
        rfcWindowStartTime_ = circuitTime;
      }

      if (!resultFound_)
      {
        // this is the simple case where Xyce output a value within a tolerance 
        // of the target value 
        if( fabs(outVarValues_[whenIdx_] - targVal) < minval_ )
        {
          updateCalculationInstant(circuitTime);
          if (findGiven_)
	  {
            updateCalculationResult(outVarValues_[0]);
	  }
          else
          {
            updateCalculationResult(circuitTime);
          }
          calculationDone_ = measureLastRFC_ ? false : doneIfFound_;
          // resultFound_ is used to control the descriptive output (to stdout) for a FIND-WHEN 
          //  measure.  If it is false, the measure shows FAILED in stdout.  This is needed for
          // compatibility with the LAST keyword, since FIND-WHEN does not use the intialized_ flag.
          resultFound_ = true;
        }
        else
        {
          // check and see if the lines defined by the current and previous values in the WHEN
          // clause indicate that the two lines, defined by those four values, have crossed.
          double prevDiff    = lastDepVarValue_ - lastTargValue_;
          double currentDiff = outVarValues_[whenIdx_] - targVal;

          // if the lines intersected then either
          //  (prevDiff < 0) && (currentDiff > 0)
          //   OR
          //  (prevDiff > 0) && (currentDiff < 0)
          // or more simply sgn( prevDiff ) = - sgn( currentDiff )
          if( ((prevDiff < 0.0) && (currentDiff > 0.0)) || ((prevDiff > 0.0) && (currentDiff < 0.0)) )
          {
            // Set the calculationInstant_ and calculationResult_ via interpolation
            interpolateResults(circuitTime, targVal);

            calculationDone_ = measureLastRFC_ ? false : doneIfFound_;
            // resultFound_ is used to control the descriptive output (to stdout) for a FIND-WHEN 
            //  measure.  If it is false, the measure shows FAILED in stdout.  This is needed for
            // compatibility with the LAST keyword, since FIND-WHEN does not use the intialized_ flag.
            resultFound_ = true;
          }
        }
      }
    }
  }

  // remember the last points in case we need to interpolate to the time when v(a)=x.  
  // lastDepVarValue_ is used to interpolate the time at which the measurement occurs.
  // lastOutputVarValue_ is used to interpolate the output value for a FIND-WHEN measure.
  lastIndepVarValue_=circuitTime;
  lastDepVarValue_=outVarValues_[whenIdx_]; 
  lastOutputVarValue_=outVarValues_[0];
  if (outputValueTargetGiven_)
    lastTargValue_ = outputValueTarget_;
  else
    lastTargValue_ = outVarValues_[whenIdx_+1];
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

    if( !calculationDone_ && withinDCsweepFromToWindow(dcSweepVal) )
    {
      outVarValues_[0] = getOutputValue(comm, outputVars_[0],
                                        solnVec, stateVec, storeVec, 0,
                                        lead_current_vector,
                                        junction_voltage_vector,
                                        lead_current_dqdt_vector, 0, 0, 0, 0);

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
      if( !initialized_ || (initialized_ && dcSweepVal == startDCMeasureWindow_) )
      {
        // Assigned last independent and dependent var to dcSweepVal and outVarValue_[whenIdx_]
        // While we can't interpolate on this step, it ensures that the initial history is
        // something realistic.
        lastOutputValue_ = outVarValues_[0]; // used for RFC counting
        lastIndepVarValue_=dcSweepVal;
        lastDepVarValue_=outVarValues_[whenIdx_];
        lastOutputVarValue_=outVarValues_[0];
        if (outputValueTargetGiven_)
          lastTargValue_ = outputValueTarget_;
        else
          lastTargValue_ = outVarValues_[whenIdx_+1];
        initialized_=true;
      }

      double targVal=0.0;
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

      if (atGiven_)
      {
        // check and see if last point and this point bound the target point
        double backDiff    = lastIndepVarValue_ - at_;
        double forwardDiff = dcSweepVal - at_;

        // if we bound the frequency target then either
        //  (backDiff < 0) && (forwardDiff > 0)
        //   OR
        //  (backDiff > 0) && (forwardDiff < 0)
        // or more simply sgn( backDiff ) = - sgn( forwardDiff )
        //
        // Also test for equality, to within the minval_ tolerance, as with the WHEN syntax.
        if ( ((backDiff < 0.0) && (forwardDiff > 0.0)) || ((backDiff > 0.0) && (forwardDiff < 0.0)) ||
	     (((abs(backDiff) < minval_) || (abs(forwardDiff) < minval_))) )
        {
          calculationResult_= outVarValues_[0] - (dcSweepVal - at_)*
	        ( (outVarValues_[0] - lastOutputVarValue_)/(dcSweepVal - lastIndepVarValue_) );
          calculationDone_ = true;
          resultFound_ = true;
        }
      }
      else if ( (type_ == "WHEN") && withinRiseFallCrossWindow(outVarValues_[whenIdx_], targVal) )
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

        if (!resultFound_)
        {
          // this is the simple case where Xyce output a value within a tolerance
          // of the target value
          if( fabs(outVarValues_[whenIdx_] - targVal) < minval_ )
          {
            updateCalculationInstant(dcSweepVal);
            if (findGiven_)
	    {
              updateCalculationResult(outVarValues_[0]);
	    }
            else
            {
              updateCalculationResult(dcSweepVal);
            }
            calculationDone_ = measureLastRFC_ ? false : doneIfFound_;
            // resultFound_ is used to control the descriptive output (to stdout) for a FIND-WHEN
            //  measure.  If it is false, the measure shows FAILED in stdout.
            resultFound_ = true;
          }
          else
          {
            // check and see if the lines defined by the current and previous values in the WHEN
            // clause indicate that the two lines, defined by those four values, have crossed.
            double prevDiff    = lastDepVarValue_ - lastTargValue_;
            double currentDiff = outVarValues_[whenIdx_] - targVal;

            // if the lines intersected then either
            //  (prevDiff < 0) && (currentDiff > 0)
            //   OR
            //  (prevDiff > 0) && (currentDiff < 0)
            // or more simply sgn( prevDiff ) = - sgn( currentDiff )
            if( ((prevDiff < 0.0) && (currentDiff > 0.0)) || ((prevDiff > 0.0) && (currentDiff < 0.0)) )
            {
              // Set the calculationInstant_ and calculationResult_ via interpolation
              interpolateResults(dcSweepVal, targVal);

              calculationDone_ = measureLastRFC_ ? false : doneIfFound_;
              // resultFound_ is used to control the descriptive output (to stdout) for a FIND-WHEN
              //  measure.  If it is false, the measure shows FAILED in stdout.
              resultFound_ = true;
            }
          }
        }
      }
    }

    // remember the last points in case we need to interpolate to the sweepVal when v(a)=x.
    // lastDepVarValue_ is used to interpolate the sweep value at which the measurement occurs.
    // lastOutputVarValue_ is used to interpolate the output value for a FIND-WHEN measure.
    lastIndepVarValue_=dcSweepVal;
    lastDepVarValue_=outVarValues_[whenIdx_];
    lastOutputVarValue_=outVarValues_[0];
    if (outputValueTargetGiven_)
      lastTargValue_ = outputValueTarget_;
    else
      lastTargValue_ = outVarValues_[whenIdx_+1];
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

  if( !calculationDone_ && withinFreqWindow(frequency) )
  {
    // update our outVarValues_ vector
    updateOutputVars(comm, outVarValues_, frequency, solnVec, 0, 0,
                     imaginaryVec, 0, 0, 0, 0, 0, 0, RFparams);

    if( !initialized_ )
    {
      // Assigned last independent and dependent var to frequency and outVarValue_[whenIdx_]
      // While we can't interpolate on this step, it ensures that the initial history is
      // something realistic.
      lastOutputValue_ = outVarValues_[0]; // used for RFC counting
      lastIndepVarValue_=frequency;
      lastDepVarValue_=outVarValues_[whenIdx_];
      lastOutputVarValue_=outVarValues_[0];
      if (outputValueTargetGiven_)
        lastTargValue_ = outputValueTarget_;
      else
        lastTargValue_ = outVarValues_[whenIdx_+1];
      initialized_=true;
    }

    double targVal=0.0;
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

    if (atGiven_)
    {
      // check and see if last point and this point bound the target point
      double backDiff    = lastIndepVarValue_ - at_;
      double forwardDiff = frequency - at_;

      // if we bound the frequency target then either
      //  (backDiff < 0) && (forwardDiff > 0)
      //   OR
      //  (backDiff > 0) && (forwardDiff < 0)
      // or more simply sgn( backDiff ) = - sgn( forwardDiff )
      //
      // Also test for equality, to within the minval_ tolerance, as with the WHEN syntax.
      if ( ((backDiff < 0.0) && (forwardDiff > 0.0)) || ((backDiff > 0.0) && (forwardDiff < 0.0)) ||
	     (((abs(backDiff) < minval_) || (abs(forwardDiff) < minval_))) )
      {
        calculationResult_= outVarValues_[0] - (frequency - at_)*
	        ( (outVarValues_[0] - lastOutputVarValue_)/(frequency - lastIndepVarValue_) );
        calculationDone_ = true;
        resultFound_ = true;
      }
    }
    else if ( (type_ == "WHEN") && withinRiseFallCrossWindow(outVarValues_[whenIdx_],targVal) )
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

      if (!resultFound_)
      {
        // this is the simple case where Xyce output a value within a tolerance
        // of the target value
        if( fabs(outVarValues_[whenIdx_] - targVal) < minval_ )
        {
          updateCalculationInstant(frequency);
          if (findGiven_)
	  {
            updateCalculationResult(outVarValues_[0]);
	  }
          else
          {
            updateCalculationResult(frequency);
          }
          calculationDone_ = measureLastRFC_ ? false : doneIfFound_;
          // resultFound_ is used to control the descriptive output (to stdout) for a FIND-WHEN
          //  measure.  If it is false, the measure shows FAILED in stdout.
          resultFound_ = true;
        }
        else
        {
          // check and see if the lines defined by the current and previous values in the WHEN
          // clause indicate that the two lines, defined by those four values, have crossed.
          double prevDiff    = lastDepVarValue_ - lastTargValue_;
          double currentDiff = outVarValues_[whenIdx_] - targVal;

          // if the lines intersected then either
          //  (prevDiff < 0) && (currentDiff > 0)
          //   OR
          //  (prevDiff > 0) && (currentDiff < 0)
          // or more simply sgn( prevDiff ) = - sgn( currentDiff )
          if( ((prevDiff < 0.0) && (currentDiff > 0.0)) || ((prevDiff > 0.0) && (currentDiff < 0.0)) )
          {
            // Set the calculationInstant_ and calculationResult_ via interpolation
            interpolateResults(frequency, targVal);

            calculationDone_ = measureLastRFC_ ? false: doneIfFound_;
            // resultFound_ is used to control the descriptive output (to stdout) for a FIND-WHEN
            //  measure.  If it is false, the measure shows FAILED in stdout.
            resultFound_ = true;
          }
        }
      }
    }
  }

  // remember the last points in case we need to interpolate to the frequency when v(a)=x.
  // lastDepVarValue_ is used to interpolate the frequency at which the measurement occurs.
  // lastOutputVarValue_ is used to interpolate the output value for a FIND-WHEN measure.
  lastIndepVarValue_=frequency;
  lastDepVarValue_=outVarValues_[whenIdx_];
  lastOutputVarValue_=outVarValues_[0];
  if (outputValueTargetGiven_)
    lastTargValue_ = outputValueTarget_;
  else
    lastTargValue_ = outVarValues_[whenIdx_+1];
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

  if( !calculationDone_ && withinFreqWindow(frequency) )
  {
    // update our outVarValues_ vector
    updateOutputVars(comm, outVarValues_, frequency, solnVec, 0, 0,
                     imaginaryVec, 0, 0, 0,
                     totalOutputNoiseDens, totalInputNoiseDens, noiseDataVec, 0);

    if( !initialized_ )
    {
      // Assigned last independent and dependent var to frequency and outVarValue_[whenIdx_]
      // While we can't interpolate on this step, it ensures that the initial history is
      // something realistic.
      lastOutputValue_ = outVarValues_[0]; // used for RFC counting
      lastIndepVarValue_=frequency;
      lastDepVarValue_=outVarValues_[whenIdx_];
      lastOutputVarValue_=outVarValues_[0];
      if (outputValueTargetGiven_)
        lastTargValue_ = outputValueTarget_;
      else
        lastTargValue_ = outVarValues_[whenIdx_+1];
      initialized_=true;
    }

    double targVal=0.0;
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

    if (atGiven_)
    {
      // check and see if last point and this point bound the target point
      double backDiff    = lastIndepVarValue_ - at_;
      double forwardDiff = frequency - at_;

      // if we bound the frequency target then either
      //  (backDiff < 0) && (forwardDiff > 0)
      //   OR
      //  (backDiff > 0) && (forwardDiff < 0)
      // or more simply sgn( backDiff ) = - sgn( forwardDiff )
      //
      // Also test for equality, to within the minval_ tolerance, as with the WHEN syntax.
      if ( ((backDiff < 0.0) && (forwardDiff > 0.0)) || ((backDiff > 0.0) && (forwardDiff < 0.0)) ||
	     (((abs(backDiff) < minval_) || (abs(forwardDiff) < minval_))) )
      {
        calculationResult_= outVarValues_[0] - (frequency - at_)*
	        ( (outVarValues_[0] - lastOutputVarValue_)/(frequency - lastIndepVarValue_) );
        calculationDone_ = true;
        resultFound_ = true;
      }
    }
    else if ( (type_ == "WHEN") && withinRiseFallCrossWindow(outVarValues_[whenIdx_],targVal) )
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

      if (!resultFound_)
      {
        // this is the simple case where Xyce output a value within a tolerance
        // of the target value
        if( fabs(outVarValues_[whenIdx_] - targVal) < minval_ )
        {
          updateCalculationInstant(frequency);
          if (findGiven_)
	  {
            updateCalculationResult(outVarValues_[0]);
	  }
          else
          {
            updateCalculationResult(frequency);
          }
          calculationDone_ = measureLastRFC_ ? false : doneIfFound_;
          // resultFound_ is used to control the descriptive output (to stdout) for a FIND-WHEN
          //  measure.  If it is false, the measure shows FAILED in stdout.
          resultFound_ = true;
        }
        else
        {
          // check and see if the lines defined by the current and previous values in the WHEN
          // clause indicate that the two lines, defined by those four values, have crossed.
          double prevDiff    = lastDepVarValue_ - lastTargValue_;
          double currentDiff = outVarValues_[whenIdx_] - targVal;

          // if the lines intersected then either
          //  (prevDiff < 0) && (currentDiff > 0)
          //   OR
          //  (prevDiff > 0) && (currentDiff < 0)
          // or more simply sgn( prevDiff ) = - sgn( currentDiff )
          if( ((prevDiff < 0.0) && (currentDiff > 0.0)) || ((prevDiff > 0.0) && (currentDiff < 0.0)) )
          {
            // Set the calculationInstant_ and calculationResult_ via interpolation
            interpolateResults(frequency, targVal);

            calculationDone_ = measureLastRFC_ ? false : doneIfFound_;
            // resultFound_ is used to control the descriptive output (to stdout) for a FIND-WHEN
            //  measure.  If it is false, the measure shows FAILED in stdout.
            resultFound_ = true;
          }
        }
      }
    }
  }

  // remember the last points in case we need to interpolate to the frequency when v(a)=x.
  // lastDepVarValue_ is used to interpolate the frequency at which the measurement occurs.
  // lastOutputVarValue_ is used to interpolate the output value for a FIND-WHEN measure.
  lastIndepVarValue_=frequency;
  lastDepVarValue_=outVarValues_[whenIdx_];
  lastOutputVarValue_=outVarValues_[0];
  if (outputValueTargetGiven_)
    lastTargValue_ = outputValueTarget_;
  else
    lastTargValue_ = outVarValues_[whenIdx_+1];
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
// Function      : FindWhenBase::interpolateResults()
// Purpose       : Interpolate the calculationInstance_ and calculationResult_
//                 values for when the measure is satisifed.  This accounts
//                 for case of WHEN V(1)=V(2) where both variables may be
//                 changing.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 02/13/2020
//-----------------------------------------------------------------------------
void FindWhenBase::interpolateResults(double currIndepVarValue, double targVal)
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

  // This is the algebra for when the two non-parallel lines associated with
  // the WHEN clause intersect.
  double ci = (d-c)/(a-b);
  updateCalculationInstant(ci);
  if (DEBUG_IO)
  {
    double targValAtCalcInstant = a*(d-c)/(a-b) + c;
    Xyce::dout() << "Target value at calculation instant for measure " << name_ << " is:"
                 << targValAtCalcInstant << std::endl;
  }

  // For a FIND measure, we need to interpolate the value of the variable in
  // the find clause and set the measure value to that interpolated value.
  if (findGiven_)
  {
    double tempVar = outVarValues_[0] - (currIndepVarValue - ci)*
	     ( (outVarValues_[0] - lastOutputVarValue_)/(currIndepVarValue - lastIndepVarValue_) );
    updateCalculationResult(tempVar);
  }
  else
    updateCalculationResult(ci);

  return;
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
  FindWhenBase(measureMgr, measureBlock)
{
  doneIfFound_=true;
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

    if ( !resultFound_ && measureMgr_.isMeasFailGiven() && measureMgr_.getMeasFail() )
    //if ( !calculationDone_ && measureMgr_.isMeasFailGiven() && measureMgr_.getMeasFail() )
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

    if (resultFound_)
    //if (calculationDone_ || ( measureLastRFC_ && resultFound_ ) )
    {
      os << name_ << " = " << this->getMeasureResult() ;
      if (atGiven_)
      {
        os << " for AT = " << at_;
      }
      else if (findGiven_)
      {
        // modeStr is "time" for TRAN or TRAN_CONT mode, "freq" for AC or AC_CONT mode
        // and "<sweep variable> value" for DC or DC_CONT mode.
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
  if (atGiven_)
    doneIfFound_=true;  // only return one value
  else
  {
    // these settings will find all times at which the WHEN clause is satisfied,
    // and may return multiple values.
    doneIfFound_=true;
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
      contCross_=1;
      contRFC_=contCross_;
      cross_=-1;
      crossGiven_=true;
    }
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
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
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
    // For negative values, store at most the requested number of values (rfcCont_). If the
    // size of this vector is less than abs(rfcCont_) then the measure is "failed".  Otherwise
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
// Purpose       :
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
    // For negative values, store at most the requested number of values (rfcCont_). If the
    // size of this vector is less than abs(rfcCont_) then the measure is "failed".  Otherwise
    // the current measure value is in calculationInstantVec_[0].
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

  // For non-negative RFC values, the calculationResultVec_ will be empty for a failed.
  // meaure.  For negative RFC values, the calculationResultVec_ will have too
  // few values for a failed measure.
  if ( (atGiven_ && !resultFound_) || (calculationResultVec_.size() == 0) ||
       ((contRFC_ < 0) && (abs(contRFC_) != calculationResultVec_.size())) )
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
  else
  {
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
  else if ( ((contRFC_ >= 0) && (calculationResultVec_.size() > 0)) ||
            ((contRFC_ < 0) && (calculationInstantVec_.size() == abs(contRFC_))) )
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

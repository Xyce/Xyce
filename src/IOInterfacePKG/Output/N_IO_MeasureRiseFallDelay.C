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

#include <N_IO_MeasureRiseFallDelay.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : RiseFallDelay::RiseFallDelay()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
RiseFallDelay::RiseFallDelay(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Base(measureMgr, measureBlock),
  trigVariableLengthHistoryNeeded_( false ),
  targVariableLengthHistoryNeeded_( false ),
  trigMax_(0.0),
  targMax_(0.0),
  trigResultIndex_(0),
  targResultIndex_(0),
  timeForTrig_(0.0),
  timeForTarg_(0.0),
  trigMaxChanged_(false),
  targMaxChanged_(false),
  timeForTrigFound_(false),
  timeForTargFound_(false),
  trigOutputValueTargetChanged_(false),
  targOutputValueTargetChanged_(false),
  prevOutputVar0_(0.0),
  prevOutputVar1_(0.0),
  targIdx_(1), // default syntax is TRIG v(a)=<val> TARG v(b)=<val>
  actualTrigRise_(0), // rise/fall/cross can be specified separately for TRIG and TARG
  actualTrigFall_(0),
  actualTrigCross_(0),
  isTrigRising_(false),
  isTrigFalling_(false),
  lastTrigOutputValue_(0.0),
  newTrigRiseWindow_(false),
  newTrigFallWindow_(false),
  newTrigCrossWindow_(false),
  actualTargRise_(0),
  actualTargFall_(0),
  actualTargCross_(0),
  isTargRising_(false),
  isTargFalling_(false),
  lastTargOutputValue_(0.0),
  newTargRiseWindow_(false),
  newTargFallWindow_(false),
  newTargCrossWindow_(false)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();

  // check if this measure will need to record any history.  It will need history
  // info if trig_frac_max or targ_frac_max are given.  (These are cases where the
  // trigger / target are a fraction of the TRIG/TARG maximum.  Thus a history will
  // be needed to find the extrema and than use it find the fraction of that extrema
  if( trigFracMaxGiven_ && !atGiven_ )
  {
    trigVariableLengthHistoryNeeded_ = true;
  }
  else
  {
    // we will have a fixed length history to bracket the 
    // trigger.
    trigIndepVarHistory_.resize(2);
    trigVarHistory_.resize(2);     
  }
  
  if( targFracMaxGiven_ )
  {
    targVariableLengthHistoryNeeded_ = true;
  }
  else
  {
    // we will have a fixed length history to bracket the 
    // target.
    targIndepVarHistory_.resize(2);
    targetVarHistory_.resize(2);   
  }

  // Set the references into the outputVarValues_ array.
  // The default values of targIdx_=1 works if the TRIG clause uses the v(a)=<val> syntax.
  // If the TRIG clause uses the AT=<time> syntax then the variable for the TARG clause is 
  // in outVarValues_[0].  If the TRIG clause uses the v(a)=v(b) syntax then the variable 
  // for the TARG clause is in outVarValues_[2].
  if (atGiven_) 
  { 
    // this handles the AT=<time> syntax in the TRIG clause
    targIdx_ = 0;
  } 
  else if (!atGiven_ && !(trigOutputValueTargetGiven_ || trigFracMaxGiven_) )
  {
    // this handles the v(a)=v(b) syntax in the TRIG clause
    targIdx_ = 2;
  }
}

//-----------------------------------------------------------------------------
// Function      : RiseFallDelay::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void RiseFallDelay::prepareOutputVariables() 
{
  // this measurement can involve up to three solution variables
  numOutVars_ = outputVars_.size();
  outVarValues_.resize( numOutVars_, 0.0 );
}

//-----------------------------------------------------------------------------
// Function      : RiseFallDelay::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void RiseFallDelay::reset() 
{
  resetBase();
  timeForTrigFound_=false;
  timeForTargFound_=false;
  if( trigVariableLengthHistoryNeeded_ )
  {
    trigIndepVarHistory_.clear();
    trigVarHistory_.clear();
  }
  if( targVariableLengthHistoryNeeded_ )
  {  
    targIndepVarHistory_.clear();
    targetVarHistory_.clear();
  }
  trigMax_ = 0.0;
  targMax_ = 0.0;
  trigResultIndex_ = 0;
  targResultIndex_ = 0;
  timeForTrig_ = 0.0;
  timeForTarg_ = 0.0;
  trigMaxChanged_ = false;
  targMaxChanged_ = false;
  timeForTrigFound_ = false;
  timeForTargFound_ = false;
  trigOutputValueTargetChanged_ = false;
  targOutputValueTargetChanged_ = false;
  prevOutputVar0_=0.0;
  prevOutputVar1_=0.0;

  actualTrigRise_=0;
  actualTrigFall_=0;
  actualTrigCross_=0;
  isTrigRising_=false;
  isTrigFalling_=false;
  lastTrigOutputValue_=0.0;
  newTrigRiseWindow_=false;
  newTrigFallWindow_=false;
  newTrigCrossWindow_=false;

  actualTargRise_=0;
  actualTargFall_=0;
  actualTargCross_=0;
  isTargRising_=false;
  isTargFalling_=false;
  lastTargOutputValue_=0.0;
  newTargRiseWindow_=false;
  newTargFallWindow_=false;
  newTargCrossWindow_=false;
}

//-----------------------------------------------------------------------------
// Function      : RiseFallDelay::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void RiseFallDelay::updateTran(
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
  if( !calculationDone_ && withinTimeWindow( circuitTime ) )
  {
    updateOutputVars(comm, outVarValues_, circuitTime,
      solnVec, stateVec, storeVec, 0, lead_current_vector,
      junction_voltage_vector, lead_current_dqdt_vector, 0, 0, 0, 0);

    // Need to set lastTrigOutputValue_ and lastTargOutputValue_ variables to the 
    // current signal values at the first time-step within the measurement window.
    // This is needed so that the RISE/FALL/CROSS count is not incremented at time=0, 
    // if the measured waveform has a DC offset at time=0    This will need to change
    // once separate TD (time delay) variables are allowed for the TRIG and TARG clauses.   
    if (!firstStepInMeasureWindow_)
    {
      if (!atGiven_) 
      { 
        lastTrigOutputValue_ = outVarValues_[0]; 
      }
      lastTargOutputValue_ = outVarValues_[targIdx_]; 
      firstStepInMeasureWindow_ = true;
    }

    // outVarValues has our TRIG and TARG values so, first store them in
    // our history array if FRAC_MAX was specified for TRIG and/or TARG
    // One memory savings is that We won't store the value if it hasn't changed
    // by less than minval_   Note: there is only one recordData flag so that
    // we only keep one history of the independent variable (time).
    bool recordData = false;

    if( trigVariableLengthHistoryNeeded_ )
    {
      int trigSize = trigVarHistory_.size();
      if( numOutVars_ > 0 )
      {
        if( trigSize == 0 )
        {
          // must record the first point
          recordData = true;
        }
        else if(fabs(trigVarHistory_[trigSize-1] - outVarValues_[0]) > minval_ )
        {
          //  trig changed enough that we need to record it.
          recordData = true;
        }
      }
    }
    else
    {
      // just store the current value and rotate the recent history variables 
      trigIndepVarHistory_[0] = trigIndepVarHistory_[1];
      trigVarHistory_[0] = trigVarHistory_[1];
      trigIndepVarHistory_[1] = circuitTime;
      trigVarHistory_[1] = outVarValues_[0];
    }

    if( targVariableLengthHistoryNeeded_ )
    {
      int targSize = targetVarHistory_.size();
      if( numOutVars_ > 0 )
      {
        if( targSize == 0 )
        {
          // must record the first point
          recordData = true;
        }
        else if( (numOutVars_ > targIdx_) && (fabs(targetVarHistory_[targSize-1] - outVarValues_[targIdx_]) > minval_ ) )
        {
          // either trig or targ changed enough that we need to record it.
          recordData = true;
        }
      }
    }
    else
    {
      // just store the current value and rotate the recent history varibles 
      targIndepVarHistory_[0] = targIndepVarHistory_[1];
      targetVarHistory_[0] = targetVarHistory_[1];
      targIndepVarHistory_[1] = circuitTime;
      targetVarHistory_[1] = outVarValues_[targIdx_];
    }
    
    // this block records data for the use case where we need a variable length
    // history of the signals (trigger and target signals) because FRAC_MAX was specified.
    if ( recordData )
    {
      // earlier checks indicated that we needed to record history
      // so do that for both trig and targ here.  We record both at each independent var
      // point so that we don't have to have two independent var histories.
      if( numOutVars_ >0 && withinTrigRiseFallCrossWindow() && !atGiven_)
      {
        trigIndepVarHistory_.push_back( circuitTime );
        trigVarHistory_.push_back( outVarValues_[0] );
        // update trigMax_ if needed.  This was changed to make frac_max work with 
        // RISE/FALL/CROSS windows for TRIG/TARG. We don't start recording history until
        // we are within the RISE/FALL/CROSS window.  So, the max value (after a FALL is 
        // sensed) may be the value at the prior time-step.
        if( outVarValues_[0] > trigMax_ || (prevOutputVar0_ > trigMax_) )
        {
          trigMax_ =  outVarValues_[0] >= prevOutputVar0_ ? outVarValues_[0] : prevOutputVar0_;
          trigMaxChanged_ = true;
        }
      }
      else if (atGiven_)
      {
        trigIndepVarHistory_.push_back( circuitTime );
      }

      if ( withinTargRiseFallCrossWindow() )  
      //if( numOutVars_ >1 && withinTargRiseFallCrossWindow() )
      {
        targIndepVarHistory_.push_back( circuitTime );
        targetVarHistory_.push_back( outVarValues_[targIdx_] );
        // update targMax_ if needed.  This was changed to make frac_max work with 
        // RISE/FALL/CROSS windows for TRIG/TARG. We don't start recording history until
        // we are within the RISE/FALL/CROSS window.  So, the max value (after a FALL is 
        // sensed) may be the value at the prior time-step.
        if( (outVarValues_[targIdx_] > targMax_) || (prevOutputVar1_ > targMax_) )
        {
          targMax_ =  outVarValues_[targIdx_] >= prevOutputVar1_ ? outVarValues_[targIdx_] : prevOutputVar1_;
          targMaxChanged_ = true;
        }
      }
    }
    
    // if TRIG uses the AT=<val> syntax then the trigOutputValueTarget_ is fixed.
    // Otherwise, in cases where trigFracMax_ is NOT given, then
    // trigOutputValueTarget_ will be fixed and is set 
    // during base class construction.  If it IS given then we may
    // need to re-calculate it here.
    // We also need to handle the v(a)=v(b) syntax in the TRIG clause
    if (!atGiven_)
    {
      if( trigFracMaxGiven_ && trigMaxChanged_ )
      {
        trigOutputValueTarget_ =  trigFracMax_ * trigMax_;
        trigOutputValueTargetChanged_ = true;
        trigMaxChanged_ = false;
      }
      else if ( !trigFracMaxGiven_ && !trigOutputValueTargetGiven_ && !atGiven_ )
      {
        // handles v(a)=v(b) in the TRIG clause
        trigOutputValueTarget_ = outVarValues_[1];
      }  
    }  

    // Rise/Fall/Cross values can be specified separately for TRIG and TARG clauses
    // Rise/Fall/Cross is not used if the TRIG clause uses the AT=<time> syntax
    // Also if FRAC_MAX is specified that Rise/Fall/Cross uses "absolute rise/fall
    // sensing with a cross level of 0.  If FRAC_MAX is not specified then the current
    // OutputTargetLevel is used as the RFC_LEVEL for rise, fall and cross.
    if (!atGiven_)
    {
      updateTrigRiseFallCrossCounts( outVarValues_[0], (trigFracMaxGiven_ ? 0 : trigOutputValueTarget_ ) );
      // If LAST was specified for TRIG then this is done
      // each time a new TRIG RFC window is entered.
      if( newTrigRiseFallCrossWindowforLast() )
      {
        timeForTrigFound_ = false;
        timeForTargFound_ = false;
	timeForTrig_ = 0;
        timeForTarg_ = 0;
	//Xyce::dout() << "found new TRIG rfc window at time= " << circuitTime << std::endl;
      }
    }
    updateTargRiseFallCrossCounts( outVarValues_[targIdx_], (targFracMaxGiven_ ? 0 : targOutputValueTarget_ ) );
    // If LAST was specified for TARG then this is done
    // each time a new TARG RFC window is entered.
    if( newTargRiseFallCrossWindowforLast() )
    {
      timeForTargFound_ = false;
      timeForTarg_ = 0;
      //Xyce::dout() << "found new TARG rfc window at time= " << circuitTime << std::endl;
    }

    // find time for TRIG
    if (atGiven_)
    {
      // at_ is a time value from the AT=<val> part of the TRIG clause
      if( (circuitTime >= at_) && !(timeForTrigFound_) )
      {
        timeForTrigFound_ = true;
        timeForTrig_ = circuitTime;
      }
    }
    else if( ( !timeForTrigFound_ || trigOutputValueTargetChanged_ ) && withinTrigRiseFallCrossWindow() )
    {
      // this block is used if the TRIG clause hasn't
      // the trigger yet or the trigger value has changed 
      int trigSize = trigVarHistory_.size();
      for( int i=trigResultIndex_;i<(trigSize-1);i++ )
      {
        // goal is to find trigFracMax_ * trigMax_
        double difference = trigVarHistory_[i] - trigOutputValueTarget_;
        double nextDifference = trigVarHistory_[i+1] - trigOutputValueTarget_;

        bool diffNeg = ( difference < 0 );
        bool nextDiffNeg = ( nextDifference < 0 );
        if( diffNeg != nextDiffNeg )
        {
          // crossed the target value
          // interpolate to get a better estimate of the targ time. 
          if( fabs(trigVarHistory_[i+1] - trigVarHistory_[i]) < minval_ )
          {
            // avoid potentially dividing by zero 
            timeForTrig_ = trigIndepVarHistory_[i];
            if (timeForTarg_ < timeForTrig_) 
            { 
              // need to re-find the TARG time, if the TRIG time has been re-set because of FRAC_MAX
              timeForTargFound_ = false;
            }
          }
          else
          {
            timeForTrig_ = (trigIndepVarHistory_[i+1] - trigIndepVarHistory_[i]) * ((trigOutputValueTarget_- trigVarHistory_[i]) / (trigVarHistory_[i+1] - trigVarHistory_[i]) ) + trigIndepVarHistory_[i];
            timeForTrigFound_=true;
            if (timeForTarg_ < timeForTrig_) 
            { 
              // need to re-find the TARG time, if the TRIG time has been re-set because of FRAC_MAX
              timeForTargFound_ = false;
            }
          }
          trigOutputValueTargetChanged_=false;
          trigResultIndex_=i;
          break;
        }
        else if ( (fabs(difference) < minval_) && (fabs(nextDifference) >= minval_) )
	{
          // equality to within MINVAL tolerance
          timeForTrig_ = trigIndepVarHistory_[i];
          timeForTrigFound_ = true;
          if (timeForTarg_ < timeForTrig_) 
          { 
              // need to re-find the TARG time, if the TRIG time has been re-set because of FRAC_MAX
              timeForTargFound_ = false;
          }
	}
      }
    }
  
    // For the TARG clause, in cases where targFracMax_ is NOT given, then
    // targOutputValueTarget_ will be fixed and is set 
    // during base class construction.  If it IS given then we may
    // need to re-calculate it here.
    // We also need to handle the v(a)=v(b) syntax in the TARG clause
    if( targFracMaxGiven_ && targMaxChanged_ )
    {
      targOutputValueTarget_ = targFracMax_ * targMax_;
      targOutputValueTargetChanged_ = true;
      targMaxChanged_ = false;
    }
    else if ( !targFracMaxGiven_ && !targOutputValueTargetGiven_ )
    {
      // handles v(a)=v(b) in the TARG clause
      targOutputValueTarget_ = outVarValues_[targIdx_+1];
    } 
  
    // find time for TARG.  We also only do this calculation if the TRIG time has been found
    // Also check that that the existing TARG time is still greater than the TRIG time, since the
    // TRIG time may have changed based on the frac_max specified for the TRIG  value 
    if ( timeForTrigFound_ && (!timeForTargFound_ || targOutputValueTargetChanged_ || (timeForTarg_ < timeForTrig_)) 
         && withinTargRiseFallCrossWindow() )
    {
      int targSize = targetVarHistory_.size();
      for( int i=targResultIndex_;i<(targSize-1); i++ )
      {
        double difference = targetVarHistory_[i] - targOutputValueTarget_;
        double nextDifference = targetVarHistory_[i+1] - targOutputValueTarget_;
        bool diffNeg = ( difference < 0 );
        bool nextDiffNeg = ( nextDifference < 0 );
        if( ( diffNeg != nextDiffNeg ) && targIndepVarHistory_[i+1] > timeForTrig_ )
        {
          // crossed the target value
          // interpolate to get a better estimate of the targ time. 
          if( fabs( targetVarHistory_[i+1] - targetVarHistory_[i]) < minval_ )
          {
            // avoid potentially dividing by zero 
            timeForTarg_ = targIndepVarHistory_[i];
          }
          else
          {
             timeForTarg_ = (targIndepVarHistory_[i+1] - targIndepVarHistory_[i]) * ((targOutputValueTarget_- targetVarHistory_[i]) /(targetVarHistory_[i+1] - targetVarHistory_[i])) + targIndepVarHistory_[i];
            timeForTargFound_=true;
          }
          targOutputValueTargetChanged_=false;
          targResultIndex_=i;
          break;
        }
        else if ( (fabs(difference) < minval_) && (fabs(nextDifference) >= minval_) 
                  && targIndepVarHistory_[i+1] > timeForTrig_ )
	{
          // equality to within MINVAL tolerance
          timeForTarg_ = targIndepVarHistory_[i];
          timeForTargFound_ = true;
	}
      }
    }

    // prune old history that isn't needed.
    // don't prune history every timestep.  Just when the trig/tarResultIndex_ is over
    // some max (set to 1000 here)
    const int pruningHistoryThreshold=1000;
    if( trigResultIndex_> pruningHistoryThreshold )
    {
      std::vector<double>::iterator pruneIttr = trigVarHistory_.begin();
      pruneIttr += trigResultIndex_;
      trigVarHistory_.erase( trigVarHistory_.begin(), pruneIttr );
      pruneIttr = trigIndepVarHistory_.begin();
      pruneIttr += trigResultIndex_;
      trigIndepVarHistory_.erase( trigIndepVarHistory_.begin(), pruneIttr );
      trigResultIndex_ = 0;
    }

    if( targResultIndex_ > pruningHistoryThreshold )
    {
      std::vector<double>::iterator pruneIttr = targetVarHistory_.begin();
      pruneIttr += targResultIndex_;
      targetVarHistory_.erase( targetVarHistory_.begin(), pruneIttr );
      pruneIttr = targIndepVarHistory_.begin();
      pruneIttr += targResultIndex_;
      targIndepVarHistory_.erase( targIndepVarHistory_.begin(), pruneIttr );
      targResultIndex_ = 0;
    }

    // added to make frac_max work with RISE/FALL/CROSS windows for TRIG/TARG.
    // Need because we don't start recording history until we are within the
    // RISE/FALL/CROSS window.  So, the max value (after a FALL is sensed) may
    // be the value at the prior time-step
    if ( numOutVars_ > 0 )
    {
      if ( !atGiven_ ) 
      {
        // outVarValues_[0] is associated with TRIG if the v(a)=<val> syntax is used
        // in the TRIG clause 
        prevOutputVar0_ = outVarValues_[0]; 
      }
      prevOutputVar1_ = outVarValues_[targIdx_];
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : RiseFallDelay::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void RiseFallDelay::updateDC(
  Parallel::Machine comm,
  const std::vector<Analysis::SweepParam> & dcParamsVec,
  const Linear::Vector *solnVec,
  const Linear::Vector *stateVec,
  const Linear::Vector *storeVec,
  const Linear::Vector *lead_current_vector,
  const Linear::Vector *junction_voltage_vector,
  const Linear::Vector *lead_current_dqdt_vector)
{

}

//-----------------------------------------------------------------------------
// Function      : RiseFallDelay::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical Models & Simulation
// Creation Date : 8/7/2019
//-----------------------------------------------------------------------------
void RiseFallDelay::updateAC(
  Parallel::Machine comm,
  const double frequency,
  const double fStart,
  const double fStop,
  const Linear::Vector *solnVec,
  const Linear::Vector *imaginaryVec,
  const Util::Op::RFparamsData *RFparams)
{

}

//-----------------------------------------------------------------------------
// Function      : RiseFallDelay::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/27/2012
//-----------------------------------------------------------------------------
double RiseFallDelay::getMeasureResult()
{
  // if both the TRIG and TARG times were not found then the measure value is default value (e.g., -1)
  if ( timeForTrigFound_ && timeForTargFound_ )
  {
    calculationResult_=timeForTarg_ - timeForTrig_;
  }
  return calculationResult_;
}

//-----------------------------------------------------------------------------
// Function      : RiseFallDelay::printMeasureResult()
// Purpose       : used to print the measurement result to an output stream
//                 object, which is typically the mt0, ma0 or ms0 file
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 2/09/2015
//-----------------------------------------------------------------------------
std::ostream& RiseFallDelay::printMeasureResult(std::ostream& os)
{
    basic_ios_all_saver<std::ostream::char_type> save(os);
    os << std::scientific << std::setprecision(precision_);

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

    return os;
}

//-----------------------------------------------------------------------------
// Function      : RiseFallDelay::printVerboseMeasureResult()
// Purpose       : used to print the "verbose" (more descriptive) measurement
//                 result to an output stream object, which is typically stdout
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 2/22/2015
//-----------------------------------------------------------------------------
std::ostream& RiseFallDelay::printVerboseMeasureResult(std::ostream& os)
{
    basic_ios_all_saver<std::ostream::char_type> save(os);
    os << std::scientific << std::setprecision(precision_);

    if ( timeForTrigFound_ && timeForTargFound_ )
    {
      os << name_ << " = " << this->getMeasureResult() ;
    }
    else
    { 
      os << name_ << " = FAILED";
    }     
    os << " with trig time= " << timeForTrig_ << " and targ time= " << timeForTarg_ << std::endl;

    return os;
}

//-----------------------------------------------------------------------------
// Function      : RiseFallDelay::withinTrigRiseFallCrossWindow
// Purpose       : Checks if current value is within measurement window for
//               : for the TRIG keyword
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 02/11/2015
//-----------------------------------------------------------------------------
bool RiseFallDelay::withinTrigRiseFallCrossWindow()
{
  bool retVal = true;
  if( trigRiseGiven_ || trigFallGiven_ || trigCrossGiven_ )
  {
    retVal = false;
    if( trigRiseGiven_ && ((trigRise_ < 0) || (trigRise_ == actualTrigRise_)))
    {
      retVal=true;
    }
    else if( trigFallGiven_ && ((trigFall_ < 0) || (trigFall_ == actualTrigFall_)))
    {
      retVal=true;
    }
    else if( trigCrossGiven_ && ((trigCross_< 0) || (trigCross_ == actualTrigCross_)))
    {
      retVal=true;
    }
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : RiseFallDelay::withinTargRiseFallCrossWindow
// Purpose       : Checks if current value is within measurement window for
//               : for the TARG keyword
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 02/11/2015
//-----------------------------------------------------------------------------
bool RiseFallDelay::withinTargRiseFallCrossWindow()
{
  bool retVal = true;
  if( targRiseGiven_ || targFallGiven_ || targCrossGiven_ )
  {
    retVal = false;

    if( targRiseGiven_ && ((targRise_ < 0) || (targRise_ == actualTargRise_)))
    {
      retVal=true;
    }
    else if( targFallGiven_ && ((targFall_ < 0) || (targFall_ == actualTargFall_)))
    {
      retVal=true;
    }
    else if( targCrossGiven_ && ((targCross_< 0) || (targCross_ == actualTargCross_)))
    {
      retVal=true;
    }
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::newTrigRiseFallCrossWindowForLast
// Purpose       : Determines if this is the start of a new
//                 Rise, Fall or Cross window.  It is used
//                 if the LAST keyword was specified for a TRIG clause.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 04/7/2016
//-----------------------------------------------------------------------------
bool RiseFallDelay::newTrigRiseFallCrossWindowforLast()
{
  bool retVal=false;

  // the variables newTrigRiseWindow_, newTrigFallWindow_ and newTrigCrossWindow_
  // were set in Base::withinTrigRiseFallCrossWindow().
  if ( (trigRise_ < 0 && trigRiseGiven_ && newTrigRiseWindow_) || 
       (trigFall_ < 0 && trigFallGiven_ && newTrigFallWindow_) || 
       (trigCross_ < 0 && trigCrossGiven_ && newTrigCrossWindow_) ) 
  {
    retVal=true;
  }

  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::newTargRiseFallCrossWindowForLast
// Purpose       : Determines if this is the start of a new
//                 Rise, Fall or Cross window.  It is used
//                 if the LAST keyword was specified for a TARG clause.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 04/7/2016
//-----------------------------------------------------------------------------
bool RiseFallDelay::newTargRiseFallCrossWindowforLast()
{
  bool retVal=false;

  // the variables newTargRiseWindow_, newTargFallWindow_ and newTargCrossWindow_
  // were set in Base::withinTargRiseFallCrossWindow().
  if ( (targRise_ < 0 && targRiseGiven_ && newTargRiseWindow_) || 
       (targFall_ < 0 && targFallGiven_ && newTargFallWindow_) || 
       (targCross_ < 0 && targCrossGiven_ && newTargCrossWindow_) ) 
  {
    retVal=true;
  }

  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : RiseFallDelay::checkMeasureLine
// Purpose       : check .MEASURE line for errors that will cause cause dumps
//               : later
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 09/01/2015
//-----------------------------------------------------------------------------
bool RiseFallDelay::checkMeasureLine()
{
  bool bsuccess = true;
  // incorrect number of dependent solution variables will cause core dumps in
  // updateTran() function
  if (numDepSolVars_ <= 1 && !atGiven_)
  {
    // TRIG clause uses v(a)=<val> syntax
    bsuccess = false;
    Report::UserError0() << name_ << " has incomplete MEASURE line";
  }
  else if (numDepSolVars_ < 1 && atGiven_)
  {
    // TRIG clause uses AT=<time> syntax
    bsuccess = false;
    Report::UserError0() << name_ << " has incomplete MEASURE line";
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : RiseFallDelay::updateTrigTargRiseFallCrossCounts
// Purpose       : Updates the actual number of rise, fall and crosses for both TRIG
//                 and TARG.  This function is not currently used.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 09/6/2015
//-----------------------------------------------------------------------------
void RiseFallDelay::updateTrigTargRiseFallCrossCounts( const double measureVal, const double crossVal,
		    const bool riseGiven, const bool fallGiven, const bool crossGiven, 
		    bool& isRising, bool& isFalling, int& actualRise, int& actualFall,
                    int& actualCross, double& lastOutputValue )
{
  if( riseGiven || fallGiven || crossGiven )
  {
    // first check if we need to adjust rise/fall/cross counts.
    if( (measureVal > lastOutputValue_) && !isRising )
    {
      // we've started a rise
      isRising= true;
      isFalling = false;
      actualRise++;
    }
    if( (measureVal < lastOutputValue) && !isFalling )
    {
      // we've started a fall
      isRising = false;
      isFalling = true;
      actualFall++;
    }
    if( (((measureVal-crossVal) < 0.0) && ((lastOutputValue-crossVal) > 0.0)) 
     || (((measureVal-crossVal) > 0.0) && ((lastOutputValue-crossVal) < 0.0)) )
    {
      // we've crossed measureVal-crossVal == 0 
      actualCross++;
    }

    lastOutputValue=measureVal;
  }
}

//-----------------------------------------------------------------------------
// Function      : RiseFallDelay::updateTrigRiseFallCrossCounts
// Purpose       : Updates the actual number of rise, fall and crosses for TRIG clause
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 02/11/2015
//-----------------------------------------------------------------------------
void RiseFallDelay::updateTrigRiseFallCrossCounts( double measureVal, double crossVal )
{
  // used to enable the LAST keyword for TRIG.  Reset to false, each time 
  // this function is called.
  newTrigRiseWindow_=false;
  newTrigFallWindow_=false;
  newTrigCrossWindow_=false;

  if( trigRiseGiven_ || trigFallGiven_ || trigCrossGiven_ )
  {
    // check if we need to adjust rise/fall/cross counts.
    if (trigFracMaxGiven_)
    {
      // sense absolute rises and falls, if using FRAC_MAX for TRIG clause
      if( (measureVal > lastTrigOutputValue_) && !isTrigRising_ )
      {
        // we've started a rise
        isTrigRising_= true;
        isTrigFalling_ = false;
        actualTrigRise_++;
      }
      if( (measureVal < lastTrigOutputValue_) && !isTrigFalling_ )
      {
        // we've started a fall
        isTrigRising_ = false;
        isTrigFalling_ = true;
        actualTrigFall_++;
      }
    }
    else
    {
      // use level-crossing approach otherwise, for HSpice compatibility. Also for 
      // HSpice compatibility, a rising (or falling) waveform that equals the cross 
      // level is considered to a rise (or fall).
      if ( ((measureVal-crossVal) >= 0.0) && ((lastTrigOutputValue_-crossVal) < 0.0) )
      {
        actualTrigRise_++;
        newTrigRiseWindow_=true;
      }
      else if( ((measureVal-crossVal) <= 0.0) && ((lastTrigOutputValue_-crossVal) > 0.0) )
      {   
        actualTrigFall_++;
        newTrigFallWindow_=true;
      }
    }

    // CROSS qualifier always uses level-crossing approach.  For HSpice compatibility, 
    // a rising (or falling) waveform that equals the cross level is considered to be
    // a cross.
    if( (((measureVal-crossVal) <= 0.0) && ((lastTrigOutputValue_-crossVal) > 0.0)) 
     || (((measureVal-crossVal) >= 0.0) && ((lastTrigOutputValue_-crossVal) < 0.0)) )
    {
      // we've crossed measureVal-crossVal == 0 
      actualTrigCross_++;
      if (!trigFracMaxGiven_) 
      {
        newTrigCrossWindow_=true;
      }
    }

    lastTrigOutputValue_=measureVal;
  }
}

//-----------------------------------------------------------------------------
// Function      : RiseFallDelay::updateTargRiseFallCrossCounts
// Purpose       : Updates the actual number of rise, fall and crosses
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 02/11/2015
//-----------------------------------------------------------------------------
void RiseFallDelay::updateTargRiseFallCrossCounts( double measureVal, double crossVal )
{
  // used to enable the LAST keyword for TARG.  Reset to false, each time 
  // this function is called.
  newTargRiseWindow_=false;
  newTargFallWindow_=false;
  newTargCrossWindow_=false;

  if( targRiseGiven_ || targFallGiven_ || targCrossGiven_)
  {
    // check if we need to adjust rise/fall/cross counts.
    if (targFracMaxGiven_)
    {
      // sense absolute rises and falls, if using FRAC_MAX for TRIG clause
      if( (measureVal > lastTargOutputValue_) && !isTargRising_ )
      {
        // we've started a rise
        isTargRising_= true;
        isTargFalling_ = false;
        actualTargRise_++;
      }
      if( (measureVal < lastTargOutputValue_) && !isTargFalling_ )
      {
        // we've started a fall
        isTargRising_ = false;
        isTargFalling_ = true;
        actualTargFall_++;
      }
    }
    else
    {
      // use level-crossing approach otherwise, for HSpice compatibility. Also for 
      // HSpice compatibility, a rising (or falling) waveform that equals the cross 
      // level is considered to a rise (or fall).
      if ( ((measureVal-crossVal) >= 0.0) && ((lastTargOutputValue_-crossVal) < 0.0) )
      {
        actualTargRise_++;
        newTargRiseWindow_=true;
      }
      else if( ((measureVal-crossVal) <= 0.0) && ((lastTargOutputValue_-crossVal) > 0.0) )
      {   
        actualTargFall_++;
        newTargFallWindow_=true;
      }
    }

    // CROSS qualifier always uses level-crossing approach.  For HSpice compatibility, 
    // a rising (or falling) waveform that equals the cross level is considered to be
    // a cross.
    if( (((measureVal-crossVal) <= 0.0) && ((lastTargOutputValue_-crossVal) > 0.0)) 
     || (((measureVal-crossVal) >= 0.0) && ((lastTargOutputValue_-crossVal) < 0.0)) )
    {
      // we've crossed measureVal-crossVal == 0 
      actualTargCross_++;
      if (!targFracMaxGiven_) 
      {
        newTargCrossWindow_=true;
      }
    }

    lastTargOutputValue_=measureVal;
  }
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

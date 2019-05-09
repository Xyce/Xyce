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

#include <N_IO_MeasureFindWhen.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : FindWhen::FindWhen()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
FindWhen::FindWhen(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Base(measureMgr, measureBlock),
  lastIndepVarValue_(0.0),
  lastDepVarValue_(0.0),
  lastOutputVarValue_(0.0),
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
  if (findGiven_) 
  { 
    whenIdx_ = 1;
  }
}

//-----------------------------------------------------------------------------
// Function      : FindWhen::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void FindWhen::prepareOutputVariables() 
{
  // this measurement can involve up to three solution variables
  numOutVars_ = outputVars_.size();
  outVarValues_.resize( numOutVars_, 0.0 );
}

//-----------------------------------------------------------------------------
// Function      : FindWhen::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void FindWhen::reset() 
{
  resetBase();
  lastIndepVarValue_=0.0;
  lastDepVarValue_=0.0;
  lastOutputVarValue_=0.0;
}


//-----------------------------------------------------------------------------
// Function      : FindWhen::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void FindWhen::updateTran(Parallel::Machine comm, const double circuitTime, const Linear::Vector *solnVec, const Linear::Vector *stateVec, const Linear::Vector *storeVec, const Linear::Vector *lead_current_vector, const Linear::Vector *junction_voltage_vector, const Linear::Vector *lead_current_dqdt_vector)
{
  if( !calculationDone_ && withinTimeWindow( circuitTime ) )
  {
    // we're in the time window, now we need to calculate the value of this
    // measure and see if it triggers any specified rise, fall, cross windowing.
    double tempResult = 0.0;

    // update our outVarValues_ vector
    updateOutputVars(comm, outVarValues_, circuitTime, solnVec, stateVec, storeVec, 0, lead_current_vector, junction_voltage_vector, lead_current_dqdt_vector );

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
      // assigned last dependent and independent var to current time and outVarValue_[whenIdx_] 
      // while we can't interpolate on this step, it ensures that the initial history is
      // something realistic.
      lastIndepVarValue_=circuitTime;
      lastDepVarValue_=outVarValues_[whenIdx_];
      lastOutputVarValue_=outVarValues_[0]; 
      initialized_=true;
    }  

    double targVal=0.0;
    bool doneIfFound=false;
    if( outputValueTargetGiven_ )
    {
      // This is the form WHEN v(a)=fixed value
      targVal = outputValueTarget_;
      doneIfFound = true;
    }
    else
    {
      // This is the form WHEN v(a)= potentially changing value, such as v(a)=v(b)
      // in that case v(b) is in outVarValues_[whenIdx_+1]
      targVal = outVarValues_[whenIdx_+1];
      // since we can't determine if the calculation is done at this point
      // we don't set calculationDone_ = true;
      //doneIfFound = false;
      // The doneIfFound usage was changed for Xyce 6.4, so that the measure returns the 
      // value (or time) at the first time that the WHEN clause is satisfied.
      doneIfFound=true;
    }
    
    // for the FIND-WHEN and WHEN measures, the rfc level used is either the target value
    // of the WHEN clause, or the value of the RFC_LEVEL qualifier if one is specified.
    if( (type_ == "WHEN") &&  withinRiseFallCrossWindow( outVarValues_[whenIdx_], 
                                     (rfcLevelGiven_ ? rfcLevel_ : targVal) ) )
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
          calculationInstant_ = circuitTime;
          if (findGiven_)
	  {
            calculationResult_ = outVarValues_[0];
	  }
          else
          {
            calculationResult_ = circuitTime;
          }
          calculationDone_ = measureLastRFC_ ? false : doneIfFound;
          // resultFound_ is used to control the descriptive output (to stdout) for a FIND-WHEN 
          //  measure.  If it is false, the measure shows FAILED in stdout.  This is needed for
          // compatibility with the LAST keyword, since FIND-WHEN does not use the intialized_ flag.
          resultFound_ = true;
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
          if( ((backDiff < 0.0) && (forwardDiff > 0.0)) || ((backDiff > 0.0) && (forwardDiff < 0.0)) )
          {
            // bound the solution so interpolate to find the target time (or frequency etc)
            calculationInstant_ = circuitTime - ( ((circuitTime - lastIndepVarValue_)/(outVarValues_[whenIdx_]-lastDepVarValue_)) * (outVarValues_[whenIdx_]-targVal) );
            if (findGiven_)
	    {
              // interpolate the value for the variable in the FIND clause, for a FIND-WHEN
              // measure, based on the interpolated measureTime variable
              calculationResult_= outVarValues_[0] - (circuitTime - calculationInstant_)*
	        ( (outVarValues_[0] - lastOutputVarValue_)/(circuitTime - lastIndepVarValue_) ); 
	    }
	    else
	    {
              // return the interpolated time if the measure is WHEN, rather than FIND-WHEN
              calculationResult_ = calculationInstant_;
	    }
            calculationDone_ = measureLastRFC_ ? false : doneIfFound;
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
}


//-----------------------------------------------------------------------------
// Function      : FindWhen::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void FindWhen::updateDC(Parallel::Machine comm, const std::vector<Analysis::SweepParam> & dcParamsVec, const Linear::Vector *solnVec, const Linear::Vector *stateVec, const Linear::Vector *storeVec, const Linear::Vector *lead_current_vector, const Linear::Vector *junction_voltage_vector, const Linear::Vector *lead_current_dqdt_vector)
{}

//-----------------------------------------------------------------------------
// Function      : FindWhen::printMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 2/22/2015
//-----------------------------------------------------------------------------
std::ostream& FindWhen::printMeasureResult(std::ostream& os, bool printVerbose)
{
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os << std::scientific << std::setprecision(precision_);
  if (!printVerbose)
  {
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
  }
  else
  {
    if (calculationDone_ || ( measureLastRFC_ && resultFound_ ) )
    {
      os << name_ << " = " << this->getMeasureResult() ;
      if (findGiven_)
      {   
        os << " at time = " << calculationInstant_;
      }
    }
    else
    { 
      os << name_ << " = FAILED";
    }
    os << std::endl;
  } 

  return os;
}

//-----------------------------------------------------------------------------
// Function      : FindWhen::checkMeasureLine
// Purpose       : check .MEASURE line for errors that will cause cause dumps
//               : later
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 09/01/2015
//-----------------------------------------------------------------------------
bool FindWhen::checkMeasureLine()
{
  bool bsuccess = true;
  // incorrect number of dependent solution variables will cause core dumps in
  // updateTran() function
  if (numDepSolVars_ <= 1 && findGiven_)
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
// Function      : FindWhen::printRFCWindow()
// Purpose       : print informaiton about the start time of the RISE or FALL
//                 window, if a valid one was found.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 09/21/2015
//-----------------------------------------------------------------------------
std::ostream& FindWhen::printRFCWindow(std::ostream& os)
{
  // no op, for any measure that supports WHEN 
 
  return os;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

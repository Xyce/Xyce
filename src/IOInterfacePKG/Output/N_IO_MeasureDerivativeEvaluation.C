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
  firstTimeValue_(0.0),
  lastTimeValue_(0.0),
  firstSignalValue_(0.0),
  lastSignalValue_(0.0),
  lastIndepVarValue_(0.0),
  lastDepVarValue_(0.0),
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
  firstTimeValue_ = 0.0;
  firstSignalValue_ = 0.0;
  lastIndepVarValue_=0.0;
  lastDepVarValue_=0.0;
  lastOutputVarValue_=0.0;
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
    // update our outVarValues_ vector
    updateOutputVars(comm, outVarValues_, circuitTime,
      solnVec, stateVec, storeVec, 0, lead_current_vector,
      junction_voltage_vector, lead_current_dqdt_vector, 0);

    if( !initialized_ )
    {
      // assigned last dependent and independent var to current time and outVarValue_[whenIdx_] 
      // while we can't interpolate on this step, it ensures that the initial history is
      // something realistic
      lastIndepVarValue_=circuitTime;
      lastDepVarValue_=outVarValues_[whenIdx_];
      lastOutputVarValue_=outVarValues_[0]; 
      initialized_=true;
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
        // The interpolated value at time = atVal_ is calculated for debugging purposes.  It is not
        // currently used in the measure calculation. 
        //double xVal = lastDepVarValue_ + (outVarValues_[0] - lastDepVarValue_)* 
	//            ((at_ - lastIndepVarValue_)/(circuitTime - lastIndepVarValue_));

        // asymmetrical 3-point approximation for first derivative.  
        calculationResult_ = (outVarValues_[0] - lastOutputVarValue_) / (circuitTime - lastIndepVarValue_); 
        calculationDone_ = true;
      }
    }
    else if ( ( outputValueTargetGiven_ || (numOutVars_ == 3) ) && !calculationDone_ 
	      && withinTimeWindow( circuitTime ) )
    {
      double targVal=0.0;
      bool doneIfFound=false;
  
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
        // since we can't determine if the calculation is done at this piont
        // we don't set calculationDone_ = true;
        //doneIfFound = false;
        // The doneIfFound usage was changed for Xyce 6.4, so that the measure returns the 
        // derivative at the first time that the WHEN clause is satisfied.
        doneIfFound = true;
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
 
        if (!resultFound_)
	{
          // this is the simple case where Xyce output a value within tolerance 
          // of the target value 
          if( fabs(outVarValues_[whenIdx_] - targVal) < minval_ )
          {
            calculationResult_ = circuitTime;
            calculationDone_ = measureLastRFC_ ? false : doneIfFound;
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
              // The interpolated value at time = atVal_ is calculated for debugging purposes.  It is not
              // currently used in the measure calculation. 
              //double xVal = lastDepVarValue_ + (outVarValues_[0] - lastDepVarValue_)* 
	      //              ((at_ - lastIndepVarValue_)/(circuitTime - lastIndepVarValue_));

              // use same time interpolation algorithm as FIND-WHEN measure
              calculationInstant_ = circuitTime - ( ((circuitTime - lastIndepVarValue_)/(outVarValues_[whenIdx_]-lastDepVarValue_)) * (outVarValues_[whenIdx_]-targVal) );

              // asymmetrical 3-point approximation for first derivative.  
              calculationResult_ = (outVarValues_[0] - lastOutputVarValue_) / (circuitTime - lastIndepVarValue_);         
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
    
    lastIndepVarValue_ = circuitTime;
    lastDepVarValue_ = outVarValues_[whenIdx_]; 
    lastOutputVarValue_=outVarValues_[0]; 
}


//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
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

}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical Models & Simulation
// Creation Date : 8/7/2019
//-----------------------------------------------------------------------------
void DerivativeEvaluation::updateAC(
  Parallel::Machine comm,
  const double frequency,
  const Linear::Vector *solnVec,
  const Linear::Vector *imaginaryVec,
  const Util::Op::RFparamsData *RFparams)
{

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
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 2/22/2015
//-----------------------------------------------------------------------------
std::ostream& DerivativeEvaluation::printMeasureResult(std::ostream& os, bool printVerbose)
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
      if (atGiven_)
      { 
        os << " for AT = " << at_;
      }
      else
      {    
        os << " at time = " << calculationInstant_;
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
  } 

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
std::ostream& DerivativeEvaluation::printMeasureWindow(std::ostream& os, const double endSimTime)
{
  //no op
  if (!atGiven_) 
  { 
    Base::printMeasureWindow(os,endSimTime); 
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

} // namespace Measure
} // namespace IO
} // namespace Xyce

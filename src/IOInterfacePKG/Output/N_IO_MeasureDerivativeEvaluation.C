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
  WhenAT(measureMgr, measureBlock, 1)
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
  numOutVars_ = outputVars_.size();
  outVarValues_.resize( numOutVars_, 0.0 );
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

  if ( !calculationDone_ && !isInvalidTimeWindow(endSimTime) )
  {
    initialized_ = true;

    if (atGiven_ && (numPointsFound_ > 1) && withinTimeWindow(at_))
    {
      if (isATcondition(circuitTime))
      {
        updateMeasureVarsForAT(circuitTime);
      }
    }
    else if ( whenGiven_ && (numPointsFound_ > 1) )
    {
      double targVal= getTargVal();

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
            updateMeasureVarsForWhen(circuitTime, whenTime);
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
        if (isATcondition(dcSweepVal))
        {
          updateMeasureVarsForAT(dcSweepVal);
        }
      }
      else if (whenGiven_ && (numPointsFound_ > 1))
      {
        double targVal=getTargVal();

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
              updateMeasureVarsForWhen(dcSweepVal, whenSweepVal);
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
// Creator       : Pete Sholander, SNL
// Creation Date : 3/25/2020
//-----------------------------------------------------------------------------
void DerivativeEvaluationBase::updateAC(
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
      if (isATcondition(frequency))
      {
        updateMeasureVarsForAT(frequency);
      }
    }
    else if (whenGiven_ && (numPointsFound_ > 1))
    {
      double targVal=getTargVal();

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
            updateMeasureVarsForWhen(frequency, whenFreq);
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
// Creator       : Pete Sholander, SNL
// Creation Date : 5/12/2010
//-----------------------------------------------------------------------------
void DerivativeEvaluationBase::updateNoise(
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
      if (isATcondition(frequency))
      {
        updateMeasureVarsForAT(frequency);
      }
    }
    else if (whenGiven_ && (numPointsFound_ > 1))
    {
      double targVal=getTargVal();

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
            updateMeasureVarsForWhen(frequency, whenFreq);
          }
        }
      }
    }
  }

  updateMeasureState(frequency);
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationBase::checkMeasureLine
// Purpose       : check .MEASURE line for errors that will cause cause dumps
//               : later
// Special Notes :
// Scope         : protected
// Creator       : Pete Sholander, SNL
// Creation Date : 10/20/2021
//-----------------------------------------------------------------------------
bool DerivativeEvaluationBase::checkMeasureLine() const
{
  bool bsuccess = true;
  // incorrect number of dependent solution variables will cause core dumps in
  // updateTran() function
  if (numDepSolVars_ == 0)
  {
    // this is wrong for any measure
    bsuccess = false;
  }
  else if ( (atGiven_ && numDepSolVars_ != 1) || 
            (!atGiven_ && outputValueTargetGiven_ && numDepSolVars_ != 2) ||
            (!atGiven_ && !outputValueTargetGiven_ && numDepSolVars_ != 3) )
  {
    // DERIV-AT
    // DERIV-WHEN with a fixed target
    // DERIV-WHEN with a variable target
    bsuccess = false;
  }

  if (!bsuccess)
     Report::UserError0() << name_ << " has invalid MEASURE line";

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationBase::updateMeasureVarsForAT
// Purpose       : Updates the calculation result, and associated flags, if
//                 the AT condition has been met.
// Special Notes : For AC and NOISE measures, the independent variable is
//                 frequency.  For DC measures, it is the value of the first
//                 variable in the DC sweep vector.  For TRAN, it is circuit
//                 time.
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 09/15/2021
//-----------------------------------------------------------------------------
void DerivativeEvaluationBase::updateMeasureVarsForAT(double currIndepVarVal)
{
  calculationResult_=getDerivativeValue(currIndepVarVal);
  calculationDone_ = true;
  resultFound_ = true;

  return;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationBase::updateMeasureVarsForWhen()
// Purpose       : Updates the calculation result and calculation instant vectors,
//                 and the associated flags, if the WHEN condition has been met.
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it
//                 is the value of the first variable in the DC sweep vector.
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 09/16/2021
//-----------------------------------------------------------------------------
void DerivativeEvaluationBase::updateMeasureVarsForWhen(double currIndepVarVal,
                                                        double whenInstance)
{
  updateCalculationInstant(whenInstance);
  updateCalculationResult(getDerivativeValue(currIndepVarVal));
  calculationDone_ = !measureLastRFC_;
  resultFound_=true;

  return;
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationBase::getDerivativeValue
// Purpose       :
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : protected
// Creator       : Pete Sholander, SNL
// Creation Date : 05/22/2020
//-----------------------------------------------------------------------------
double DerivativeEvaluationBase::getDerivativeValue(double currIndepVarVal) const
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
  DerivativeEvaluationBase(measureMgr, measureBlock)
{}

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
  resetWhenAT();
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::printMeasureResult()
// Purpose       : used to print the measurement result to an output stream
//                 object, which is typically the mt0, ma0 or ms0 file
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/21/2020
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
// Function      : DerivativeEvaluation::printVerboseMeasureResult()
// Purpose       : used to print the "verbose" (more descriptive) measurement
//                 result to an output stream object, which is typically stdout
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/21/2020
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
  DerivativeEvaluationBase(measureMgr, measureBlock)
{
  // this setting will find all times at which the WHEN clause is satisfied,
  // and may return multiple values.
  measureLastRFC_=true;
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
  resetWhenAT();
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluationCont::printMeasureResult()
// Purpose       : used to print the measurement result to an output stream
//                 object, which is typically the mt0, ma0 or ms0 file
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
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
  else if (whenGiven_ && (((RFC_ >= 0) && (calculationResultVec_.size() > 0)) ||
			  ((RFC_ < 0) && (calculationResultVec_.size() == abs(RFC_)))) )
  {
    // For non-negative RFC values, the calculationResultVec_ will be non-empty for a successful
    // FIND-WHEN or WHEN meaure.  For negative RFC values, the calculationResultVec_ will
    // have the correct number of values in it.
    if (RFC_ >=0)
    {
      for (size_t i=0; i<calculationResultVec_.size(); i++)
        os << name_ << " = " << calculationResultVec_[i] << std::endl;
    }
    else
    {
      os << name_ << " = " << calculationResultVec_[0] << std::endl;
    }
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
// Function      : DerivativeEvaluationCont::printVerboseMeasureResult()
// Purpose       : used to print the "verbose" (more descriptive) measurement
//                 result to an output stream object, which is typically stdout
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
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
  else if (whenGiven_ && (((RFC_ >= 0) && (calculationResultVec_.size() > 0)) ||
			  ((RFC_ < 0) && (calculationResultVec_.size() == abs(RFC_)))) )
  {
    // modeStr is "time" for TRAN mode, "freq" for AC mode and
    // "<sweep variable> value" for DC mode.
    std::string modeStr = setModeStringForMeasureResultText();
    if (RFC_ >=0)
    {
      for (size_t i=0; i<calculationResultVec_.size(); i++)
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

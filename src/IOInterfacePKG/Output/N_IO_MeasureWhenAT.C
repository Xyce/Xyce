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

#include <N_IO_MeasureBase.h>
#include <N_IO_MeasureWhenAT.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : WhenAT::WhenAT()
// Purpose       : Base class for measures (DERIV-AT, DERIV-WHEN, FIND-AT,
//                 FIND-WHEN and WHEN) that use an AT or WHEN clause
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 9/27/2021
//-----------------------------------------------------------------------------
  WhenAT::WhenAT(const Manager &measureMgr, const Util::OptionBlock & measureBlock, int whenIdx):
  Base(measureMgr, measureBlock),
  whenIdx_(whenIdx),
  RFC_(0),
  lastIndepVarValue_(0.0),
  lastDepVarValue_(0.0),
  lastOutputVarValue_(0.0),
  lastTargValue_(0.0),
  startDCMeasureWindow_(0.0),
  numPointsFound_(0)
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
// Function      : WhenAT::resetWhenAT
// Purpose       : Called when restarting a measure function.
// Scope         : protected
// Creator       : Pete Sholander, SNL
// Creation Date : 9/27/2021
//-----------------------------------------------------------------------------
void WhenAT::resetWhenAT()
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
// Function      : WhenAT::interpolateCalculationInstant
// Purpose       : Interpolate the time for when the measure is satisifed.
//                 This accounts for case of WHEN V(1)=V(2) where both
//                 variables may be changing.
// Special Notes :
// Scope         : protected
// Creator       : Pete Sholander, SNL
// Creation Date : 08/04/2020
//-----------------------------------------------------------------------------
double WhenAT::interpolateCalculationInstant(double currIndepVarValue, double targVal) const
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
// Function      : WhenAT::isATcondition
// Purpose       : Evaluates if the AT condition is true for all measure modes.
// Special Notes : For AC and NOISE measures, the independent variable is
//                 frequency.  For DC measures, it is the value of the first
//                 variable in the DC sweep vector.  For TRAN, it is circuit
//                 time.
// Scope         : protected
// Creator       : Pete Sholander, SNL
// Creation Date : 05/21/2020
//-----------------------------------------------------------------------------
bool WhenAT::isATcondition(double indepVarVal) const
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
	   (((fabs(backDiff) < minval_) || (fabs(forwardDiff) < minval_))) );
}

//-----------------------------------------------------------------------------
// Function      : WhenAT::isWHENcondition
// Purpose       : Evaluates if the WHEN condition is true for all modes
// Special Notes : For AC and NOISE measures, the independent variable is
//                 frequency.  For DC measures, it is the value of the first
//                 variable in the DC sweep vector.  For TRAN measures, it
//                 is the circuit time.
// Scope         : protected
// Creator       : Pete Sholander, SNL
// Creation Date : 05/21/2020
//-----------------------------------------------------------------------------
bool WhenAT::isWHENcondition(double indepVarVal, double targVal) const
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
// Function      : WhenAT::updateRFCcountForWhen
// Purpose       : Updates the rise, fall and cross counts.
// Special Notes : This function is normally only called for a WHEN time that
//                 falls within the FROM-TO window.  So, the WHEN time is always
//                 a valid cross.
// Scope         : protected
// Creator       : Pete Sholander, SNL
// Creation Date : 08/06/2020
//-----------------------------------------------------------------------------
void WhenAT::updateRFCcountForWhen()
{
  actualCross_++;
  if (outVarValues_[whenIdx_] > lastDepVarValue_)
    actualRise_++;
  else
    actualFall_++;

  return;
}

//-----------------------------------------------------------------------------
// Function      : WhenAT::withinRFCWindowForWhen
// Purpose       : Determine if the current WHEN time falls within the specified
//                 RFC value.
// Special Notes : Assumes that a WHEN measure without RISE, FALL or CROSS given
//                 on the .MEASURE line defaults to CROSS given.
// Scope         : protected
// Creator       : Pete Sholander, SNL
// Creation Date : 08/06/2020
//-----------------------------------------------------------------------------
bool WhenAT::withinRFCWindowForWhen() const
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
// Function      : WhenAT::setMeasureState
// Purpose       : initializes the past values of the independent, dependent
//                 and measured variables, as well as the past target level.
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : protected
// Creator       : Pete Sholander, SNL
// Creation Date : 08/04/2020
//-----------------------------------------------------------------------------
void WhenAT::setMeasureState(double indepVarVal)
{
  // assigned last dependent and independent var to current value of the independent
  // varible and outVarValue_[whenIdx_].  While we can't interpolate on this step, it
  // ensures that the initial history is something realistic.
  lastIndepVarValue_=indepVarVal;
  lastDepVarValue_=outVarValues_[whenIdx_];
  lastOutputVarValue_=outVarValues_[0];
  updateLastTargVal();

  return;
}

//-----------------------------------------------------------------------------
// Function      : WhenAT::updateMeasureState
// Purpose       : updates the past values of the independent, dependent
//                 and measured variables, as well as the past target level.
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : protected
// Creator       : Pete Sholander, SNL
// Creation Date : 08/04/2020
//-----------------------------------------------------------------------------
void WhenAT::updateMeasureState(double indepVarVal)
{
  lastIndepVarValue_ = indepVarVal;
  lastDepVarValue_ = outVarValues_[whenIdx_];
  lastOutputVarValue_=outVarValues_[0];
  updateLastTargVal();

  return;
}

//-----------------------------------------------------------------------------
// Function      : WhenAT::getTargVal
// Purpose       : updates the target value for the WHEN clause
// Special Notes :
// Scope         : protected
// Creator       : Pete Sholander, SNL
// Creation Date : 08/04/2020
//-----------------------------------------------------------------------------
double WhenAT::getTargVal() const
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
// Function      : WhenAT::updateLastTargVal
// Purpose       : updates the last target value for the WHEN clause
// Special Notes :
// Scope         : protected
// Creator       : Pete Sholander, SNL
// Creation Date : 09/16/2020
//-----------------------------------------------------------------------------
void WhenAT::updateLastTargVal()
{
  if (outputValueTargetGiven_)
    lastTargValue_ = outputValueTarget_;
  else
    lastTargValue_ = outVarValues_[whenIdx_+1];

  return;
}

//-----------------------------------------------------------------------------
// Function      : WhenAT::updateCalculationResult
// Purpose       : Updates the vector that holds the measure values.  This
//                 vector may hold multiple values if the RISE, FALL or CROSS
//                 value is < 0.
// Special Notes : For compatibility with other measure types, the calculationResult_
//                 variable is also updated with the current measure value.
// Scope         : protected
// Creator       : Pete Sholander, SNL
// Creation Date : 08/03/2020
//-----------------------------------------------------------------------------
void WhenAT::updateCalculationResult(double val)
{
  if (RFC_ >= 0)
  {
    // store all of the values once the requested number of rises (or falls or crosses)
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
// Function      : WhenAT::updateCalculationInstant
// Purpose       : Updates the vector that holds the times (or frequencies or
//                 DC sweep values) when the measure was satisified.  This
//                 vector may hold multiple values if the RISE, FALL or CROSS
//                 value is <0.
// Special Notes :
// Scope         : protected
// Creator       : Pete Sholander, SNL
// Creation Date : 08/20/2020
//-----------------------------------------------------------------------------
void WhenAT::updateCalculationInstant(double val)
{
  if (RFC_ >= 0)
  {
    // store all of the values once the requested number of rises (or falls or crosses)
    // has been found
    if ( (riseGiven_ && (actualRise_>= rise_)) || (fallGiven_ && (actualFall_>= fall_)) ||
         (crossGiven_ && (actualCross_>= cross_)) )
    {
      calculationInstantVec_.push_back(val);
    }
  }
  else
  {
    // For negative values, store at most the requested number of values (RFC_). If the
    // size of this vector is less than abs(RFC_) then the measure is "failed".  Otherwise
    // the current instant value is in calculationInstantVec_[0].
    calculationInstantVec_.push_back(val);
    if (calculationInstantVec_.size() > abs(RFC_))
      calculationInstantVec_.erase(calculationInstantVec_.begin());
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : WhenAT::printMeasureWindow
// Purpose       : prints information related to measure window
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 04/8/2020
//-----------------------------------------------------------------------------
std::ostream& WhenAT::printMeasureWindow(std::ostream& os, double endSimTime,
				               double startSweepVal, double endSweepVal) const
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
// Function      : WhenAT::printRFCWindow()
// Purpose       : print informaiton about the start time of the RISE or FALL
//                 window, if a valid one was found.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 09/21/2015
//-----------------------------------------------------------------------------
std::ostream& WhenAT::printRFCWindow(std::ostream& os) const
{
  // no op, for any measure that supports WHEN

  return os;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

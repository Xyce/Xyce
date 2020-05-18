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

#include <N_IO_MeasureRMS.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : RMS::RMS()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
RMS::RMS(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Stats(measureMgr, measureBlock),
  integralXsq_(0.0),
  totalIntegrationWindow_(0.0)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();
}


//-----------------------------------------------------------------------------
// Function      : RMS::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void RMS::reset() 
{
  resetStats();
  integralXsq_ = 0.0;
  totalIntegrationWindow_ = 0.0;
}


//-----------------------------------------------------------------------------
// Function      : RMS::setMeasureVarsForNewWindow()
// Purpose       : Called when entering a new RFC window for TRAN measures
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 04/28/2020
//-----------------------------------------------------------------------------
void RMS::setMeasureVarsForNewWindow()
{
  integralXsq_ = 0.0;
  totalIntegrationWindow_ = 0.0;
  initialized_ = false;

  return;
}


//-----------------------------------------------------------------------------
// Function      : RMS::updateMeasureVars()
// Purpose       : Updates the RMS calculation
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 measures, it is frequency.  For DC measures, it is the value
//                 of the first variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 04/28/2020
//-----------------------------------------------------------------------------
void RMS::updateMeasureVars(const double indepVarVal, const double signalVal)
{
  integralXsq_ += 0.5* (indepVarVal - lastIndepVarValue_) * (signalVal*signalVal + lastSignalValue_*lastSignalValue_);
  totalIntegrationWindow_ += (indepVarVal - lastIndepVarValue_);

  return;
}


//-----------------------------------------------------------------------------
// Function      : RMS::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double RMS::getMeasureResult()
{
  if ( initialized_ && (numPointsFound_ > 1) )
      calculationResult_ =  sqrt(integralXsq_ / totalIntegrationWindow_);

  return calculationResult_;
}

//-----------------------------------------------------------------------------
// Function      : RMS::printMeasureResult()
// Purpose       : used to print the measurement result to an output stream
//                 object, which is typically the mt0, ma0 or ms0 file
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 5/13/2020
//-----------------------------------------------------------------------------
std::ostream& RMS::printMeasureResult(std::ostream& os)
{
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os << std::scientific << std::setprecision(precision_);

  if ( (!initialized_ || (numPointsFound_ < 2)) && measureMgr_.isMeasFailGiven() && measureMgr_.getMeasFail() )
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
// Function      : RMS::printVerboseMeasureResult()
// Purpose       : used to print the measurement result to an output stream
//                 object, which is typically the mt0, ma0 or ms0 file
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 5/13/2020
//-----------------------------------------------------------------------------
std::ostream& RMS::printVerboseMeasureResult(std::ostream& os)
{
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os << std::scientific << std::setprecision(precision_);

  if ( initialized_ && (numPointsFound_ > 1) )
  {
    os << name_ << " = " << this->getMeasureResult() << std::endl;
  }
  else
  {
    os << name_ << " = FAILED" << std::endl;
  }

  return os;
}

//-----------------------------------------------------------------------------
// Function      : RMS::printMeasureWindow
// Purpose       : prints information related to measure window
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/25/2020
//-----------------------------------------------------------------------------
std::ostream& RMS::printMeasureWindow(std::ostream& os, const double indepVarValue)
{
  // Pathological case of FROM=TO within an otherwise valid FROM-TO window.
  // This a failed measure, but the FROM-TO window should be printed correctly.
  if ( (fromGiven_ || toGiven_) && (from_==to_) && firstSweepValueFound_ &&
       ((mode_ == "AC") || (mode_ == "DC")) )
  {
    basic_ios_all_saver<std::ostream::char_type> save(os);
    os << std::scientific << std::setprecision(precision_);
    std::string modeStr = setModeStringForMeasureWindowText();
    os << "Measure Start " << modeStr << "= " << startACDCmeasureWindow_
       << "\tMeasure End " << modeStr << "= " << endACDCmeasureWindow_ << std::endl;
  }
  else
  {
    Base::printMeasureWindow(os,indepVarValue);
  }

  return os;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

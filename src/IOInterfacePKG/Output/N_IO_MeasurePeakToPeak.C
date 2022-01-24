//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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

#include <N_IO_MeasurePeakToPeak.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : PeakToPeak::PeakToPeak()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
PeakToPeak::PeakToPeak(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Extrema(measureMgr, measureBlock),
  maximumValue_(0.0),
  minimumValue_(0.0),
  maximumInstant_(0.0),
  minimumInstant_(0.0)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();
}


//-----------------------------------------------------------------------------
// Function      : PeakToPeak::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void PeakToPeak::reset() 
{
  resetBase();
  maximumValue_ = 0.0;
  minimumValue_ = 0.0;
  maximumInstant_ = 0.0;
  minimumInstant_ = 0.0;
}


//-----------------------------------------------------------------------------
// Function      : PeakToPeak::setMeasureVarsForNewWindow()
// Purpose       : Called when entering a new RFC window for TRAN measures, or
//                 a new FROM-TO window for AC, DC and NOISE measures.
// Special Notes : For TRAN measure, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 04/28/2020
//-----------------------------------------------------------------------------
void PeakToPeak::setMeasureVarsForNewWindow(double indepVarVal, double depVarVal)
{
  maximumValue_ = depVarVal;
  minimumValue_ = depVarVal;
  maximumInstant_ = indepVarVal;
  minimumInstant_ = indepVarVal;
  initialized_ = true;

  return;
}

//-----------------------------------------------------------------------------
// Function      : PeakToPeak::updateMeasureVars()
// Purpose       : Updates the maximum and minimum values, and their
//                 calculation instants
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 04/28/2020
//-----------------------------------------------------------------------------
void PeakToPeak::updateMeasureVars(double indepVarVal, double depVarVal)
{
  // update the maximum and minimum values.  The Peak-to-Peak value
  // is calculated in the PeakToPeak::getMeasureResult() function.
  if( maximumValue_ < depVarVal )
  {
    maximumValue_ = depVarVal;
    maximumInstant_ = indepVarVal;
  }

  if( minimumValue_ > depVarVal )
  {
    minimumValue_ = depVarVal;
    minimumInstant_ = indepVarVal;
  }

  return;
}


//-----------------------------------------------------------------------------
// Function      : PeakToPeak::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double PeakToPeak::getMeasureResult()
{
  if( initialized_ )
  {
    calculationResult_ =  maximumValue_ - minimumValue_;
  }
  return calculationResult_;
}

//-----------------------------------------------------------------------------
// Function      : PeakToPeak::printMeasureResult()
// Purpose       : used to print the measurement result to an output stream
//                 object, which is typically the mt0, ma0 or ms0 file
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 2/09/2015
//-----------------------------------------------------------------------------
std::ostream& PeakToPeak::printMeasureResult(std::ostream& os)
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
// Function      : PeakToPeak::printVerboseMeasureResult()
// Purpose       : used to print the "verbose" (more descriptive) measurement
//                 result to an output stream object, which is typically stdout
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 2/22/2015
//-----------------------------------------------------------------------------
std::ostream& PeakToPeak::printVerboseMeasureResult(std::ostream& os)
{
    basic_ios_all_saver<std::ostream::char_type> save(os);
    os << std::scientific << std::setprecision(precision_);

    if (initialized_)
    {
      os << name_ << " = " << this->getMeasureResult() ;
      // modeStr is "time" for TRAN mode, "freq" for AC and NOISE modes, and
      // "<sweep variable> value" for DC mode.
      std::string modeStr = setModeStringForMeasureResultText();          
      os << " with max at " << modeStr << "= " << maximumInstant_ << 
      " and min at " << modeStr << "= " << minimumInstant_<< std::endl;
    }
    else
    { 
      os << name_ << " = FAILED" << std::endl;
    }

    return os;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

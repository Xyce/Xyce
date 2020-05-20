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

#include <N_IO_MeasureMin.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : Min::Min()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
Min::Min(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Extrema(measureMgr, measureBlock),
  minimumValue_(0.0)

{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();
}


//-----------------------------------------------------------------------------
// Function      : Min::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void Min::reset() 
{
  resetExtrema();
  minimumValue_ = 0.0;
}

//-----------------------------------------------------------------------------
// Function      : Min::setMeasureVarsForNewWindow()
// Purpose       : Called when entering a new RFC window for TRAN measures, or
//                 a new FROM-TO window for AC, DC and NOISE measures.
// Special Notes : For TRAN measure, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 04/28/2020
//-----------------------------------------------------------------------------
void Min::setMeasureVarsForNewWindow(const double indepVarVal, const double depVarVal)
{
  minimumValue_ = depVarVal;
  calculationInstant_ = indepVarVal;
  initialized_ = true;

  return;
}

//-----------------------------------------------------------------------------
// Function      : Min::updateMeasureVars()
// Purpose       : Updates the minimum value and calculation instant
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 04/28/2020
//-----------------------------------------------------------------------------
void Min::updateMeasureVars(const double indepVarVal, const double depVarVal)
{
  // calculation of the minimum value
  if( minimumValue_ > depVarVal )
  {
    minimumValue_ = depVarVal;
    calculationInstant_ = indepVarVal;
  }

  return;
}


//-----------------------------------------------------------------------------
// Function      : Min::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double Min::getMeasureResult()
{
  if( initialized_ )
  {
    calculationResult_ =  minimumValue_;
  }
  return calculationResult_;
}

//-----------------------------------------------------------------------------
// Function      : Min::printMeasureResult()
// Purpose       : used to print the measurement result to an output stream
//                 object, which is typically the mt0, ma0 or ms0 file
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 2/09/2015
//-----------------------------------------------------------------------------
std::ostream& Min::printMeasureResult(std::ostream& os)
{
    basic_ios_all_saver<std::ostream::char_type> save(os);
    os << std::scientific << std::setprecision(precision_);

    if ( !initialized_ && measureMgr_.isMeasFailGiven() && measureMgr_.getMeasFail() )
    {
      // output FAILED to .mt file if .OPTIONS MEASURE MEASFAIL=1 is given in the
      // netlist and this is a failed measure.
      os << name_ << " = FAILED" << std::endl;
    }
    else if ( (measureOutputOption_ == "TIME") || (measureOutputOption_ == "FREQ")
         || (measureOutputOption_ == "SV") )
    {
      // output the time (or frequency or value of the first variable in 
      // the DC sweep vector) when the minimum value occurs
      os << name_ << " = " << calculationInstant_ << std::endl;
    }
    else
    {
      // output the minimum value
      os << name_ << " = " << this->getMeasureResult() << std::endl;
    }

    return os;
}

//-----------------------------------------------------------------------------
// Function      : Min::printVerboseMeasureResult()
// Purpose       : used to print the "verbose" (more descriptive) measurement
//                 result to an output stream object, which is typically stdout
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 2/22/2015
//-----------------------------------------------------------------------------
std::ostream& Min::printVerboseMeasureResult(std::ostream& os)
{
    basic_ios_all_saver<std::ostream::char_type> save(os);
    os << std::scientific << std::setprecision(precision_);

    if (initialized_)
    {
      os << name_ << " = " << this->getMeasureResult() ;

      // modeStr is "time" for TRAN mode, "freq" for AC and NOISE modes, and
      // "<sweep variable> value" for DC mode.
      std::string modeStr = setModeStringForMeasureResultText();          
      os << " at " << modeStr << " = " << calculationInstant_ << std::endl;
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

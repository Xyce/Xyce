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

#include <N_IO_MeasureOffTime.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : OffTime::OffTime()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
OffTime::OffTime(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Base(measureMgr, measureBlock),
  totalOffTime_(0.0),
  numberOfCycles_(0),
  lastTimeValue_(0.0),
  lastSignalValue_(0.0)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();
}

//-----------------------------------------------------------------------------
// Function      : OffTime::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void OffTime::prepareOutputVariables() 
{
// this measurement should have only one dependent variable.
  // Error for now if it doesn't
  numOutVars_ = outputVars_.size();

  if ( numOutVars_ > 1 )
  {
    std::string msg = "Too many dependent variables for statistical measure, \"" + name_ + "\" Exiting.";
    Report::UserFatal() << msg;
  }

  outVarValues_.resize( numOutVars_ );
}

//-----------------------------------------------------------------------------
// Function      : OffTime::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void OffTime::reset() 
{
  resetBase();
  totalOffTime_ = 0.0;
  numberOfCycles_ = 0.0;
}


//-----------------------------------------------------------------------------
// Function      : OffTime::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void OffTime::updateTran(Parallel::Machine comm, const double circuitTime, const Linear::Vector *solnVec, const Linear::Vector *stateVec, const Linear::Vector *storeVec, const Linear::Vector *lead_current_vector, const Linear::Vector *junction_voltage_vector, const Linear::Vector *lead_current_dqdt_vector)
{
  if( !calculationDone_ && withinTimeWindow( circuitTime ) )
  {
    // we're in the time window, now we need to calculate the value of this
    // measure and see if it triggers any specified rise, fall, cross windowing.

    // update our outVarValues_ vector
    updateOutputVars(comm, outVarValues_, circuitTime, solnVec, stateVec, storeVec, 0, lead_current_vector, junction_voltage_vector, lead_current_dqdt_vector );

      if( initialized_  )
      {
        if( (outVarValues_[0] - minval_) <= offValue_ )
        {
          // add to Off duty time
          totalOffTime_ += (circuitTime - lastTimeValue_);

          // did we just cross into a new cycle?
          if( lastSignalValue_ > offValue_ )
          {
            numberOfCycles_++;
          }
        }
      }

      lastTimeValue_ = circuitTime;
      lastSignalValue_ = outVarValues_[0];
      initialized_=true;
  }
}


//-----------------------------------------------------------------------------
// Function      : OffTime::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void OffTime::updateDC(Parallel::Machine comm, const std::vector<Analysis::SweepParam> & dcParamsVec, const Linear::Vector *solnVec, const Linear::Vector *stateVec, const Linear::Vector *storeVec, const Linear::Vector *lead_current_vector, const Linear::Vector *junction_voltage_vector, const Linear::Vector *lead_current_dqdt_vector)
{

}


//-----------------------------------------------------------------------------
// Function      : OffTime::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double OffTime::getMeasureResult()
{
  // This will return the default value (e.g., -1) if the measure is never
  // initialized.  If the time window is valid but the OFF condition is
  // never met then this will return 0.0
  if( initialized_ )
  {
    calculationResult_ =  totalOffTime_;
    if( initialized_ && (numberOfCycles_ > 0))
    {
      calculationResult_ /=  numberOfCycles_;
    }
  }
  return calculationResult_;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

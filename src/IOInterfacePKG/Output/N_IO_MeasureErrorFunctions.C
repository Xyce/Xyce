//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Purpose       : Implement ERR, ERR1 and ERR2 measure types
// Special Notes : The Err1 class implements both the ERR and ERR1 measure types
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_MeasureErrorFunctions.h>
#include <N_DEV_DeviceSupport.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : ErrorFunctions::ErrorFunctions()
// Purpose       : Class for functions common to Err1 and Err2 classes
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-----------------------------------------------------------------------------
ErrorFunctions::ErrorFunctions(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Base(measureMgr, measureBlock)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();
}

//-----------------------------------------------------------------------------
// Function      : ErrorFunctions::prepareOutputVariables()
// Purpose       : Makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-----------------------------------------------------------------------------
void ErrorFunctions::prepareOutputVariables()
{
  numOutVars_ = outputVars_.size();
  outVarValues_.resize( numOutVars_, 0.0 );
}

//-----------------------------------------------------------------------------
// Function      : ErrorFunctions::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-----------------------------------------------------------------------------
void ErrorFunctions::updateTran(
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
  if( !calculationDone_ && withinTimeWindow( circuitTime ) )
  {
    // update our outVarValues_ vector
    updateOutputVars(comm, outVarValues_, circuitTime,
      solnVec, stateVec, storeVec, 0, lead_current_vector,
      junction_voltage_vector, lead_current_dqdt_vector, 0, 0, 0, 0);

    initialized_ = true;
    if ( withinYLimits_(outVarValues_[0]) )
      updateErrVars(outVarValues_[0], outVarValues_[1]);
  }
}

//-----------------------------------------------------------------------------
// Function      : ErrorFunctions::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-----------------------------------------------------------------------------
void ErrorFunctions::updateDC(
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
      // update our outVarValues_ vector
      updateOutputVars(comm, outVarValues_, dcSweepVal,
        solnVec, stateVec, storeVec, 0, lead_current_vector,
        junction_voltage_vector, lead_current_dqdt_vector, 0, 0, 0, 0);

      initialized_ = true;
      if ( withinYLimits_(outVarValues_[0]) )
        updateErrVars(outVarValues_[0], outVarValues_[1]);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : ErrorFunctions::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-----------------------------------------------------------------------------
void ErrorFunctions::updateAC(
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

  if( !calculationDone_ && withinFreqWindow( frequency ) )
  {
    // update our outVarValues_ vector
    updateOutputVars(comm, outVarValues_, frequency, solnVec, 0, 0,
                     imaginaryVec, 0, 0, 0, 0, 0, 0, RFparams);

    initialized_ = true;
    if ( withinYLimits_(outVarValues_[0]) )
      updateErrVars(outVarValues_[0], outVarValues_[1]);
  }
}

//-----------------------------------------------------------------------------
// Function      : ErrorFunctions::updateNoise()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical Models & Simulation
// Creation Date : 5/18/2010
//-----------------------------------------------------------------------------
void ErrorFunctions::updateNoise(
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

  if( !calculationDone_ && withinFreqWindow( frequency ) )
  {
    // update our outVarValues_ vector
    updateOutputVars(comm, outVarValues_, frequency, solnVec, 0, 0,
                     imaginaryVec, 0, 0, 0,
                     totalOutputNoiseDens, totalInputNoiseDens, noiseDataVec, 0);

    initialized_ = true;
    if ( withinYLimits_(outVarValues_[0]) )
      updateErrVars(outVarValues_[0], outVarValues_[1]);
  }
}

//-----------------------------------------------------------------------------
// Function      : ErrorFunctions::checkMeasureLine
// Purpose       : check .MEASURE line for errors that will cause cause dumps
//               : later
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-----------------------------------------------------------------------------
bool ErrorFunctions::checkMeasureLine() const
{
  bool bsuccess = true;
  // incorrect number of dependent solution variables will cause core dumps in
  // updateTran() function
  if (numDepSolVars_ !=2)
  {
    bsuccess = false;
    Report::UserError0() << name_ << " has incomplete MEASURE line";
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : ErrorFunctions::withinYLimits_
// Purpose       : use YMAX and YMIN to filter out "too large" and "too small"
//                 values from measure result
// Special Notes :
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-----------------------------------------------------------------------------
bool ErrorFunctions::withinYLimits_(double mVal) const
{
  return (fabs(mVal) <= ymax_) && (fabs(mVal) >= ymin_);
}

//-----------------------------------------------------------------------------
// Function      : Err1::Err1()
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-----------------------------------------------------------------------------
Err1::Err1(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  ErrorFunctions(measureMgr, measureBlock),
  err1SqSum_(0.0),
  numPts_(0)
{}

//-----------------------------------------------------------------------------
// Function      : Err1::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-----------------------------------------------------------------------------
void Err1::reset()
{
  resetBase();
  err1SqSum_ = 0.0;
  numPts_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : Err1::updateErrVars()
// Purpose       : Update the variables used to calculate the measure result
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-----------------------------------------------------------------------------
void Err1::updateErrVars(double mVal, double cVal)
{
  ++numPts_;
  double denom = std::max(fabs(mVal), minval_);
  err1SqSum_ += ((mVal - cVal)/denom) * ((mVal - cVal)/denom);

  return;
}

//-----------------------------------------------------------------------------
// Function      : Err1::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-----------------------------------------------------------------------------
double Err1::getMeasureResult()
{
  if( initialized_ )
  {
    calculationResult_ =  sqrt(err1SqSum_ / numPts_);
  }
  return calculationResult_;
}

//-----------------------------------------------------------------------------
// Function      : Err2::Err2()
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-----------------------------------------------------------------------------
Err2::Err2(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  ErrorFunctions(measureMgr, measureBlock),
  err2Sum_(0.0),
  numPts_(0)
{}

//-----------------------------------------------------------------------------
// Function      : Err2::reset()
// Purpose       : Called when restarting a measure function.  Resets any state.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-----------------------------------------------------------------------------
void Err2::reset()
{
  resetBase();
  err2Sum_ = 0.0;
  numPts_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : Err2::updateErrVars()
// Purpose       : Update the variables used to calculate the measure result
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-----------------------------------------------------------------------------
void Err2::updateErrVars(double mVal, double cVal)
{
  ++numPts_;
  double denom = std::max(fabs(mVal), minval_);
  err2Sum_ += fabs((mVal - cVal)/denom);
 
  return;
}

//-----------------------------------------------------------------------------
// Function      : Err2::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-----------------------------------------------------------------------------
double Err2::getMeasureResult()
{
  if( initialized_ )
  {
    calculationResult_ =  err2Sum_ / numPts_;
  }
  return calculationResult_;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

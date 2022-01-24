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
// Purpose       : Measure statistics of a simulation variable over
//                 an interval
// Special Notes : This class contains the functions that are common to
//                 the Average, IntegralEvaluation and RMS classes.  It
//                 sits between those classes and the Base class.
//
// Creator       : Pete Sholander, SNL
// Creation Date : 04/28/2020
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_MeasureStats.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : Stats::Stats()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 4/28/2020
//-----------------------------------------------------------------------------
Stats::Stats(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Base(measureMgr, measureBlock),
  lastIndepVarValue_(0.0),
  lastSignalValue_(0.0),
  numPointsFound_(0)
{}

//-----------------------------------------------------------------------------
// Function      : Stats::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void Stats::prepareOutputVariables()
{
  // this measurement should have only one dependent variable.
  // Error out if it doesn't.
  numOutVars_ = outputVars_.size();

  if ( numOutVars_ > 1 )
  {
    std::string msg = "Too many dependent variables for " + type_ + " measure, \"" + name_ + "\"";
    Report::UserError0() << msg;
  }
  outVarValues_.resize( numOutVars_, 0.0 );
}


//-----------------------------------------------------------------------------
// Function      : Stats::resetStats()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void Stats::resetStats()
{
  resetBase();
  numPointsFound_=0;
}

//-----------------------------------------------------------------------------
// Function      : Stats::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 4/28/2020
//-----------------------------------------------------------------------------
void Stats::updateTran(
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
      junction_voltage_vector, lead_current_dqdt_vector, 0, 0, 0 ,0);

    if( initialized_ )
    {
      updateMeasureVars(circuitTime, outVarValues_[0]);
    }

    updateMeasureState_(circuitTime, outVarValues_[0]);
  }
}

//-----------------------------------------------------------------------------
// Function      : Stats::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/21/2020
//-----------------------------------------------------------------------------
void Stats::updateDC(
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
      outVarValues_[0] = getOutputValue(comm, outputVars_[0],
                                        solnVec, stateVec, storeVec, 0,
                                        lead_current_vector,
                                        junction_voltage_vector,
                                        lead_current_dqdt_vector, 0, 0 , 0, 0);

      if ( initialized_ )
        updateMeasureVars(dcSweepVal, outVarValues_[0]);

      updateMeasureState_(dcSweepVal, outVarValues_[0]);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Stats::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/21/2020
//-----------------------------------------------------------------------------
void Stats::updateAC(
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

    if ( initialized_ )
      updateMeasureVars(frequency, outVarValues_[0]);

    updateMeasureState_(frequency, outVarValues_[0]);
  }
}

//-----------------------------------------------------------------------------
// Function      : Stats::updateNoise()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 5/11/2020
//-----------------------------------------------------------------------------
void Stats::updateNoise(
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

    if ( initialized_ )
      updateMeasureVars(frequency, outVarValues_[0]);

    updateMeasureState_(frequency, outVarValues_[0]);
  }
}

//-----------------------------------------------------------------------------
// Function      : Stats::updateMeasureState()
// Purpose       : Updates the past values of the independent and dependent
//                 variables.
// Special Notes : For TRAN measures, the independent variable is time.  For AC
//                 and NOISE measures, it is frequency.  For DC measures, it is
//                 the value of the first variable in the DC sweep vector.
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 04/28/2020
//-----------------------------------------------------------------------------
void Stats::updateMeasureState_(double indepVarVal, double depVarVal)
{
  lastIndepVarValue_ = indepVarVal;
  lastSignalValue_ = depVarVal;
  initialized_=true;
  ++numPointsFound_;

  return;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

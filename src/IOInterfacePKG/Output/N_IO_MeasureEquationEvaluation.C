//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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

#include <N_IO_MeasureEquationEvaluation.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : EquationEvaluation::EquationEvaluation()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
EquationEvaluation::EquationEvaluation(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Base(measureMgr, measureBlock)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;
  type_ = "EQN";  // change "PARAM" to "EQN"

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();

  // change FFT mode EQN/PARAM measure to be a TRAN mode one
  if (mode_ == "FFT")
    mode_ = "TRAN";
}

//-----------------------------------------------------------------------------
// Function      : EquationEvaluation::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void EquationEvaluation::prepareOutputVariables() 
{
  // this measurement should have only one dependent variable.
  // Error out if it doesn't
  numOutVars_ = outputVars_.size();
  
  if ( numOutVars_ > 1 )
  {
    std::string msg = "Too many dependent variables for EQN measure, \"" + name_ + "\"";
    Report::UserError0() << msg;
  }

  outVarValues_.resize( numOutVars_, 0.0 );
 
}

//-----------------------------------------------------------------------------
// Function      : EquationEvaluation::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void EquationEvaluation::reset() 
{
  resetBase();
}


//-----------------------------------------------------------------------------
// Function      : EquationEvaluation::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void EquationEvaluation:: updateTran(
  Parallel::Machine comm,
  const double circuitTime,
  const double endSimTime,
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
    for( int i=0; i< numOutVars_; i++ )
    {
      outVarValues_[i] = getOutputValue(comm, outputVars_[i], solnVec, stateVec, storeVec, 0, lead_current_vector, junction_voltage_vector, lead_current_dqdt_vector, 0, 0, 0, 0);
    }

    // not intuitive, but the output of this measure is just the outputVars_ operator evaluated 
    // within the FromToWindow.  At this time there shouldn't be more than one outVarValues_ so just
    // take the first element.
    initialized_ = true;
    calculationResult_=outVarValues_[0];
  }
}

//-----------------------------------------------------------------------------
// Function      : EquationEvaluation::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void EquationEvaluation::updateDC(
  Parallel::Machine comm,
  const std::vector<Analysis::SweepParam> & dcParamsVec,
  const Linear::Vector *solnVec,
  const Linear::Vector *stateVec,
  const Linear::Vector *storeVec,
  const Linear::Vector *lead_current_vector,
  const Linear::Vector *junction_voltage_vector,
  const Linear::Vector *lead_current_dqdt_vector)
{
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
      for( int i=0; i< numOutVars_; i++ )
      {
        outVarValues_[i] = getOutputValue(comm, outputVars_[i],
                                          solnVec, stateVec, storeVec, 0,
                                          lead_current_vector,
                                          junction_voltage_vector,
                                          lead_current_dqdt_vector, 0, 0, 0, 0);
      }
      // not intuitive, but the output of this measure is just the outputVars_ operator evaluated 
      // within the FromToWindow.  At this time there shouldn't be more than one outVarValues_ so just
      // take the first element.
      initialized_ = true;
      calculationResult_=outVarValues_[0];
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : EquationEvaluation::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2014
//-----------------------------------------------------------------------------
void EquationEvaluation::updateAC(
  Parallel::Machine comm,
  const double frequency,
  const double fStart,
  const double fStop,
  const Linear::Vector *solnVec,
  const Linear::Vector *imaginaryVec,
  const Util::Op::RFparamsData *RFparams)
{
  // Used in descriptive output to stdout
  firstSweepValueFound_ = true;

  if( !calculationDone_ && withinFreqWindow( frequency ) )
  {
    // update our outVarValues_ vector
    for( int i=0; i< numOutVars_; i++ )
    {
      outVarValues_[i] = getOutputValue(comm, outputVars_[i], solnVec, 0, 0,
                                        imaginaryVec, 0, 0, 0, 0, 0, 0, RFparams);
    }
    // not intuitive, but the output of this measure is just the outputVars_ operator evaluated 
    // within the FromToWindow.  At this time there shouldn't be more than one outVarValues_ so just
    // take the first element.
    initialized_ = true;
    calculationResult_=outVarValues_[0];
  }
}

//-----------------------------------------------------------------------------
// Function      : EquationEvaluation::updateNoise()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 5/11/2020
//-----------------------------------------------------------------------------
void EquationEvaluation::updateNoise(
  Parallel::Machine comm,
  const double frequency,
  const double fStart,
  const double fStop,
  const Linear::Vector *solnVec,
  const Linear::Vector *imaginaryVec,
  const double totalOutputNoiseDens,
  const double totalInputNoiseDens,
  const std::vector<Xyce::Analysis::NoiseData*> *noiseDataVec)
{
  // Used in descriptive output to stdout
  firstSweepValueFound_ = true;

  if( !calculationDone_ && withinFreqWindow( frequency ) )
  {
    // update our outVarValues_ vector
    for( int i=0; i< numOutVars_; i++ )
    {
      outVarValues_[i] = getOutputValue(comm, outputVars_[i], solnVec, 0, 0,
                                        imaginaryVec, 0, 0, 0, totalOutputNoiseDens, totalInputNoiseDens, noiseDataVec, 0);
    }
    // not intuitive, but the output of this measure is just the outputVars_ operator evaluated 
    // within the FromToWindow.  At this time there shouldn't be more than one outVarValues_ so just
    // take the first element.
    initialized_ = true;
    calculationResult_=outVarValues_[0];
  }
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

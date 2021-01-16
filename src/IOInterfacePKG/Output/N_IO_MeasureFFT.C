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
// Purpose       :  Support various FFT-based measures
//
// Special Notes :
//
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_MeasureFFT.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : FFT::FFT()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-----------------------------------------------------------------------------
FFT::FFT(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  Base(measureMgr, measureBlock),
  fftAnalysisPtr_(0)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();
}

//-----------------------------------------------------------------------------
// Function      : FFT::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void FFT::prepareOutputVariables()
{
  // this measurement should have only one dependent variable.
  // Error out if it doesn't
  numOutVars_ = outputVars_.size();

  if ( numOutVars_ > 1 )
  {
    std::string msg = "Too many dependent variables for " + type_ +  " measure, \"" + name_ + "\"";
    Report::UserError0() << msg;
  }

  outVarValues_.resize( numOutVars_, 0.0 );
}


//-----------------------------------------------------------------------------
// Function      : FFT::resetFFT()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void FFT::resetFFT()
{
  resetBase();
}

//-----------------------------------------------------------------------------
// Function      : FFT::fixupFFTMeasure
// Purpose       : Pass info from associated FFTAnalysis object to this measure
//                 object.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/18/2021
//-----------------------------------------------------------------------------
void FFT::fixupFFTMeasure(FFTAnalysis* fftAnalysisPtr)
{
  initialized_ = true;
  fftAnalysisPtr_ = fftAnalysisPtr;
}

//-----------------------------------------------------------------------------
// Function      : FFT::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/18/2021
//-----------------------------------------------------------------------------
void FFT::updateTran(
  Parallel::Machine comm,
  const double circuitTime,
  const Linear::Vector *solnVec,
  const Linear::Vector *stateVec,
  const Linear::Vector *storeVec,
  const Linear::Vector *lead_current_vector,
  const Linear::Vector *junction_voltage_vector,
  const Linear::Vector *lead_current_dqdt_vector)
{
  //  No Op for this measure mode, although this function is actually called
}


//-----------------------------------------------------------------------------
// Function      : FFT::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/18/2021
//-----------------------------------------------------------------------------
void FFT::updateDC(
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
// Function      : FFT::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/18/2021
//-----------------------------------------------------------------------------
void FFT::updateAC(
  Parallel::Machine comm,
  const double frequency,
  const Linear::Vector *solnVec,
  const Linear::Vector *imaginaryVec,
  const Util::Op::RFparamsData *RFparams)
{

}

//-----------------------------------------------------------------------------
// Function      : FFT::updateNoise()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical Models & Simulation
// Creation Date : 1/18/2021
//-----------------------------------------------------------------------------
void FFT::updateNoise(
  Parallel::Machine comm,
  const double frequency,
  const Linear::Vector *solnVec,
  const Linear::Vector *imaginaryVec,
  const double totalOutputNoiseDens,
  const double totalInputNoiseDens,
  const std::vector<Xyce::Analysis::NoiseData*> *noiseDataVec)
{

}

//-----------------------------------------------------------------------------
// Function      : FFTFind::FFTFind()
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-----------------------------------------------------------------------------
FFTFind::FFTFind(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  FFT(measureMgr, measureBlock)
{}

//-----------------------------------------------------------------------------
// Function      : FFTFind::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-----------------------------------------------------------------------------
void FFTFind::reset()
{
  resetFFT();
}

//-----------------------------------------------------------------------------
// Function      : FFTFind::getMeasureResult()
// Purpose       : Return the FFT magnitude at a specified frequency.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/18/2021
//-----------------------------------------------------------------------------
double FFTFind::getMeasureResult()
{
  if( initialized_ )
  {
    calculationResult_ = 1;
  }
  return calculationResult_;
}

//-----------------------------------------------------------------------------
// Function      : ENOB::ENOB()
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-----------------------------------------------------------------------------
ENOB::ENOB(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  FFT(measureMgr, measureBlock)
{}

//-----------------------------------------------------------------------------
// Function      : ENOB::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-----------------------------------------------------------------------------
void ENOB::reset()
{
  resetFFT();
}

//-----------------------------------------------------------------------------
// Function      : ENOB::getMeasureResult()
// Purpose       : Return Effective Number Of Bits (ENOB).  The units are "bits".
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/18/2021
//-----------------------------------------------------------------------------
double ENOB::getMeasureResult()
{
  if( initialized_ )
  {
    calculationResult_ = fftAnalysisPtr_->getENOB();
  }
  return calculationResult_;
}

//-----------------------------------------------------------------------------
// Function      : SFDR::SFDR()
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-----------------------------------------------------------------------------
SFDR::SFDR(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  FFT(measureMgr, measureBlock)
{}

//-----------------------------------------------------------------------------
// Function      : SFDR::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-----------------------------------------------------------------------------
void SFDR::reset()
{
  resetFFT();
}

//-----------------------------------------------------------------------------
// Function      : SFDR::getMeasureResult()
// Purpose       : Return Spurious Free Dynamic Range (SFDR) in dB
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/18/2021
//-----------------------------------------------------------------------------
double SFDR::getMeasureResult()
{
  if( initialized_ )
  {
    calculationResult_ = fftAnalysisPtr_->getSFDR();
  }
  return calculationResult_;
}

//-----------------------------------------------------------------------------
// Function      : SNDR::SNDR()
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-----------------------------------------------------------------------------
SNDR::SNDR(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  FFT(measureMgr, measureBlock)
{}

//-----------------------------------------------------------------------------
// Function      : SNDR::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-----------------------------------------------------------------------------
void SNDR::reset()
{
  resetFFT();
}

//-----------------------------------------------------------------------------
// Function      : SNDR::getMeasureResult()
// Purpose       : Return Signal to Noise-plus-Distortion Ratio (SNDR) in dB
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/18/2021
//-----------------------------------------------------------------------------
double SNDR::getMeasureResult()
{
  if( initialized_ )
  {
    calculationResult_ = fftAnalysisPtr_->getSNDR();
  }
  return calculationResult_;
}

//-----------------------------------------------------------------------------
// Function      : THD::THD()
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-----------------------------------------------------------------------------
THD::THD(const Manager &measureMgr, const Util::OptionBlock & measureBlock):
  FFT(measureMgr, measureBlock)
{}

//-----------------------------------------------------------------------------
// Function      : THD::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-----------------------------------------------------------------------------
void THD::reset()
{
  resetFFT();
}

//-----------------------------------------------------------------------------
// Function      : THD::getMeasureResult()
// Purpose       : Return Total Harmonic Distortion (THD) in dB
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/18/2021
//-----------------------------------------------------------------------------
double THD::getMeasureResult()
{
  if( initialized_ )
  {
    calculationResult_ = 20*log10(fftAnalysisPtr_->getTHD());
  }
  return calculationResult_;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

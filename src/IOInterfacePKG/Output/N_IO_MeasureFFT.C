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
  fftAnalysisPtr_(0),
  np_(0),
  atRounded_(0)
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
  if (isOpTypeAllowed())
  {
    initialized_ = true;
    fftAnalysisPtr_ = fftAnalysisPtr;
    np_ = fftAnalysisPtr_->getNP();

    if (findGiven_ && atGiven_)
      atRounded_ = std::round(at_/fftAnalysisPtr_->getFreq());
  }
}

//-----------------------------------------------------------------------------
// Function      : FFT::isOpTypeAllowed
// Purpose       : Determine if the specified operator type is allowed for a given
//                 FFT measure type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/26/2021
//-----------------------------------------------------------------------------
bool FFT::isOpTypeAllowed()
{
  bool bsuccess = true;

  std::string measureVarName = outputVars_[0]->getName();
  size_t parenIdx = measureVarName.find_first_of('(');
  if ((measureVarName[0] != '{') && (parenIdx != 1))
  {
    bsuccess = false;
    Report::UserError0() << "Complex operators such as " << measureVarName.substr(0,parenIdx)
                         <<  " not allowed for output variable for " << type_ <<  " measure "
			 << name_ << " for FFT measure mode";
  }

  return bsuccess;
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
  FFT(measureMgr, measureBlock),
  opType_("M")
{
  if (!findGiven_ || !atGiven_)
  {
    Xyce::Report::UserError0() << "Only FIND-AT supported for FFT measure " << name_;
  }
}

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
// Function      : FFTFind::isOpTypeAllowed
// Purpose       : Determine if the specified operator type is allowed for the
//                 FFTFind measure type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/26/2021
//-----------------------------------------------------------------------------
bool FFTFind::isOpTypeAllowed()
{
  bool bsuccess = true;

  std::string measureVarName = outputVars_[0]->getName();
  size_t parenIdx = measureVarName.find_first_of('(');
  if (measureVarName[0] == '{')
  {
    bsuccess = false;
    Report::UserError0() << "Expressions not allowed for output variable for FIND measure "
			 << name_ << " for FFT measure mode";
  }
  else if ( !((measureVarName[0] == 'V') || (measureVarName[0] == 'I')) )
  {
    bsuccess = false;
    Report::UserError0() << "Only V and I operators allowed for output variable for FIND measure "
			 << name_ << " for FFT measure mode";
  }
  else if (parenIdx > 1)
  {
    // get OpType for VR, IR, etc.
    opType_ = measureVarName.substr(1,parenIdx-1);
  }

  return bsuccess;
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
  if (initialized_ && fftAnalysisPtr_->isCalculated() && (atRounded_ >= 0) && (atRounded_ <= np_/2))
  {
    if (opType_ == "R")
      calculationResult_ = fftAnalysisPtr_->getFFTCoeffRealVal(atRounded_);
    else if (opType_ == "I")
      calculationResult_ = fftAnalysisPtr_->getFFTCoeffImagVal(atRounded_);
    else if (opType_ == "M")
      calculationResult_ = fftAnalysisPtr_->getMagVal(atRounded_);
    else if (opType_ == "P")
      calculationResult_ = fftAnalysisPtr_->getPhaseVal(atRounded_);
    else if (opType_ == "DB")
      calculationResult_ = 20*log10(fftAnalysisPtr_->getMagVal(atRounded_));
  }

  return calculationResult_;
}

//-----------------------------------------------------------------------------
// Function      : FFTFind::printVerboseMeasureResult()
// Purpose       : used to print the "verbose" (more descriptive) measurement
//                 result to an output stream object, which is typically stdout
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/21/2021
//-----------------------------------------------------------------------------
std::ostream& FFTFind::printVerboseMeasureResult(std::ostream& os)
{
    basic_ios_all_saver<std::ostream::char_type> save(os);
    os << std::scientific << std::setprecision(precision_);

    if (initialized_ && fftAnalysisPtr_->isCalculated())
    {
      os << name_ << " = " << this->getMeasureResult()
         << " at " <<  atRounded_ << " Hz (rounded from " << at_ << " Hz)" << std::endl;
    }
    else
    {
      os << name_ << " = FAILED" << std::endl;
    }

    return os;
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
  if( initialized_ && fftAnalysisPtr_->isCalculated() )
  {
    double thd=0;

    if (!nbHarmGiven_ || ((2*nbHarm_ >= np_) || (nbHarm_ <=0)) )
    {
      // use THD value calculated with all of the harmonics
      thd = fftAnalysisPtr_->getTHD();
    }
    else if (nbHarm_ == 1)
    {
      // to match HSPICE, rather than returning an inf
      thd = fftAnalysisPtr_->getNoiseFloor();
    }
    else
    {
      // need to recalculate the THD for the requested number of harmonics
      std::vector<double> mag = fftAnalysisPtr_->getMagVec();
      for (int i=2; i<=nbHarm_; i++)
        thd += mag[i]*mag[i];
      thd = sqrt(thd)/mag[1];
    }

    // return measure result in dB
    calculationResult_ = 20*log10(thd);
  }

  return calculationResult_;
}

//-----------------------------------------------------------------------------
// Function      : THD::printVerboseMeasureResult()
// Purpose       : used to print the "verbose" (more descriptive) measurement
//                 result to an output stream object, which is typically stdout
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/21/2021
//-----------------------------------------------------------------------------
std::ostream& THD::printVerboseMeasureResult(std::ostream& os)
{
    basic_ios_all_saver<std::ostream::char_type> save(os);
    os << std::scientific << std::setprecision(precision_);

    if (initialized_ && fftAnalysisPtr_->isCalculated())
    {
      os << name_ << " = " << this->getMeasureResult();
      if (nbHarmGiven_ && (nbHarm_ >=1))
      {
        int maxHarm = np_/2;
        os << " up to the harmonic: " << std::min(nbHarm_, maxHarm);
      }
      os << std::endl;
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

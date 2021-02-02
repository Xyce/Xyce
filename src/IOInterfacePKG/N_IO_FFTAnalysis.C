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
// Purpose       : This file contains the functions to perform an FFT analysis
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_DEV_DeviceMgr.h>
#include <N_ERH_ErrorMgr.h>

#include <N_IO_CircuitBlock.h>
#include <N_IO_FFTMgr.h>
#include <N_IO_FFTAnalysis.h>
#include <N_IO_NetlistImportTool.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>

#include <N_UTL_BreakPoint.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Interpolators.h>
#include <N_UTL_Math.h>
#include <N_UTL_Op.h>
#include <N_UTL_OpBuilder.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_SaveIOSState.h>

#include <Teuchos_ScalarTraits.hpp>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : FFTAnalysis::FFTAnalysis
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
FFTAnalysis::FFTAnalysis(const Util::OptionBlock & fftBlock )
  : startTime_(0.0),
    stopTime_(0.0),
    np_(1024),
    format_("NORM"),
    windowType_("RECT"),
    alpha_(3.0),
    freq_(0.0),
    fmin_(0.0),
    fmax_(0.0),
    startTimeGiven_(false),
    stopTimeGiven_(false),
    freqGiven_(false),
    fminGiven_(false),
    fmaxGiven_(false),
    calculated_(false),
    outputVarName_(""),
    fft_accurate_(1),
    fftout_(0),
    sampleIdx_(0),
    sampleTimeTol_(1e-20),
    noiseFloor_(1e-10),
    maxMag_(0.0),
    thd_(0.0),
    sndr_(0.0),
    enob_(0.0),
    sfdr_(0.0),
    sfdrIndex_(0)
{
  // based on what's in the option block passed in, we create the needed fft instance
  Util::ParamList variableList;   // Used to help register lead current requests with device manager.

  for (Util::ParamList::const_iterator currentParamIt = fftBlock.begin(); currentParamIt != const_cast<Util::OptionBlock &>(fftBlock).end(); ++currentParamIt)
  {
    std::string name = "";
    std::string tag = currentParamIt->tag();

    if ( tag == "FREQ" )
    {
      freq_ = currentParamIt->getImmutableValue<double>();
      freqGiven_ = true;

      // need to add error checking for the FREQ qualifier.
    }
    else if ( (tag == "V") || (tag[0] == 'I') || (tag == "P") || (tag == "W") || (tag == "N") )
    {
      // tag[0] is used for I because branch currents for transistors can have two
      // characters.  An example is IS for the M-Device.
      int nodes = currentParamIt->getImmutableValue<int>();
      Util::Param aParam;
      aParam.set( tag, nodes );
      outputVarName_ += tag;
      outputVarName_ += "(";

      // add to list of variables.  This will be used later in netlist parsing
      // to enable lead currents in the device manager.
      variableList.push_back(aParam);

      numDepSolVars_++;
      depSolVarIterVector_.push_back(aParam);

      // handle both V(a) and V(a,b) syntaxes
      for( int i=0; i<nodes; i++ )
      {
        currentParamIt++;
        aParam.set( (*currentParamIt).tag(), (*currentParamIt).getImmutableValue<double>() );
        depSolVarIterVector_.push_back( aParam );
        outputVarName_ += currentParamIt->tag();
        if (i != nodes-1 && nodes > 1) outputVarName_ += ",";
      }
      outputVarName_ += ")";
    }
    else if ( Xyce::Util::hasExpressionTag(tag) )
    {
      numDepSolVars_++;
      depSolVarIterVector_.push_back(*currentParamIt);

      outputVarName_ = tag;

      // add to list of variables.  This will be used later in netlist parsing
      // to enable lead currents in the device manager.
      Util::Param aParam;
      aParam.set( currentParamIt->tag(), currentParamIt->tag() );
      variableList.push_back(aParam);

    }
    else if ( (tag == "START") || (tag == "FROM") )
    {
      startTime_ = currentParamIt->getImmutableValue<double>();
      startTimeGiven_ = true;
      if ( startTime_ < 0 )
      {
        startTime_ = 0.0;
        Report::UserWarning0() << "Invalid start time on .FFT line, reset to transient start time of 0";
      }
    }
    else if ( (tag == "STOP") || (tag == "TO") )
    {
      stopTime_ = currentParamIt->getImmutableValue<double>();
      stopTimeGiven_ = true;
    }
    else if ( tag == "FORMAT")
    {
      ExtendedString tmpStr(currentParamIt->stringValue());
      format_ = tmpStr.toUpper();
      if ( !( (format_ == "NORM") || (format_ == "UNORM") ) )
      {
	Report::UserError0() << "Invalid FORMAT type " << format_ << " on .FFT line";
      }
    }
    else if ( tag == "NP" )
    {
      np_ = currentParamIt->getImmutableValue<int>();

      if (np_ <= 0)
      {
        Report::UserError0() << "NP value on .FFT line should be a power of 2, and >=4";
      }
      else if ((np_ >=1) && (np_ < 4))
      {
        np_ = 4;
	Report::UserWarning0() << "NP value on .FFT line reset to minimum value of 4";
      }
      else if ( (np_ > 0) && !( (np_ & (np_-1)) == 0) )
      {
        // Handle cases where NP is not equal to a power of two.  Values at midpoint round up.
        // Example is 48 rounds up to 64.  47 rounds down to 32.  (Note: This matches what HSPICE
        // does, and not what the manual says its does.  The manual says to "round up".
        int floorNP = std::floor(log2(np_));
        int ceilNP = std::ceil(log2(np_));
        int npExp = floorNP + std::round( (np_-std::pow(2,floorNP))/(std::pow(2,ceilNP)-std::pow(2,floorNP)) );
        np_ = std::pow(2, npExp);
        Report::UserWarning0() << "NP value on .FFT line not equal to power of 2.  Setting to " << np_;
      }
    }
    else if ( tag == "WINDOW")
    {
      ExtendedString tmpStr(currentParamIt->stringValue());
      windowType_ = tmpStr.toUpper();
      if ( !( (windowType_ == "RECT") || (windowType_ == "BART") || (windowType_ == "HANN") ||
              (windowType_ == "HAMM") || (windowType_ == "BLACK") || (windowType_ == "HARRIS") ||
              (windowType_ == "GAUSS") || (windowType_ == "KAISER") ) )
      {
	Report::UserError0() << "Invalid WINDOW type " << windowType_ << " on .FFT line";
      }
    }
    else if (tag == "ALFA")
    {
      // Parameter used in GAUSS and KAISER windows.  It is bounded by 1.0 and 20.0
      alpha_ = currentParamIt->getImmutableValue<double>();
      alpha_ = std::max(alpha_,1.0);
      alpha_ = std::min(alpha_,20.0);
    }

    else if (tag == "FMIN")
    {
      fmin_ = currentParamIt->getImmutableValue<double>();
      fminGiven_ = true;

      // need to add error checking for the FMIN qualifier.
    }
    else if (tag == "FMAX")
    {
      fmax_ = currentParamIt->getImmutableValue<double>();
      fmaxGiven_ = true;

      // need to add error checking for the FMAX qualifier.
    }
  }

  // Create FFT
  ftInData_.resize( np_);
  ftOutData_.resize( np_ + 2 );
  iftInData_.resize( np_ + 2 );
  iftOutData_.resize( np_ );

  ftInterface_ = Teuchos::rcp( new N_UTL_FFTInterface<std::vector<double> >( np_ ) );
  ftInterface_->registerVectors( ftInData_, &ftOutData_, iftInData_, &iftOutData_ );
}

//-----------------------------------------------------------------------------
// Function      : FFTAnalysis::FFTAnalysis
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
FFTAnalysis::~FFTAnalysis()
{
  for (Util::Op::OpList::iterator it = outputVars_.begin(); it != outputVars_.end(); ++it)
    delete *it;
}

//-----------------------------------------------------------------------------
// Function      : FFTAnalysis::fixupFFTParameters
// Purpose       : This function sets private variables that could be not set in
//                 the constructor, primarily because the end simulation time was
//                 not available at that point.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
void FFTAnalysis::fixupFFTParameters(Parallel::Machine comm,
  const Util::Op::BuilderManager &op_builder_manager,
  const double endSimTime,
  TimeIntg::StepErrorControl & sec,
  const int fft_accurate,
  const bool fftout)
{
  // set these to match the values stored in the FFTMgr class.
  fft_accurate_ = fft_accurate;
  fftout_ = fftout;

  // set some defaults, and additional error checking on in put values from the FFT line
  if (!stopTimeGiven_)
  {
    stopTime_ = endSimTime;
  }
  else if ( (stopTime_ > endSimTime) || (stopTime_ <= startTime_) )
  {
    stopTime_ = endSimTime;
    Report::UserWarning0() << "Invalid stop time on .FFT line, reset to transient stop time of " << endSimTime;
  }

  if (!freqGiven_)
    freq_ = 1.0/(stopTime_ - startTime_);

  if (!fminGiven_)
    fmin_ = 1.0/(stopTime_ - startTime_);

  if (!fmaxGiven_)
    fmax_ = 0.5*np_*fmin_;

  // set up vectors for the sample times and the sampled/interpolated data values.
  sampleTimes_.resize(np_,0.0);
  sampleValues_.resize(np_,0.0);

  mag_.resize(np_+1,0.0);
  phase_.resize(np_+1,0.0);
  fftRealCoeffs_.resize(np_+1,0.0);
  fftImagCoeffs_.resize(np_+1,0.0);

  // Compute new, equally spaced time points.
  double step = (stopTime_ - startTime_)/np_;
  if (startTimeGiven_)
  {
    sampleTimes_[0] = startTime_;
    sec.setBreakPoint (sampleTimes_[0]);
  }

  for (int i = 1; i < np_; i++)
  {
    sampleTimes_[i] = sampleTimes_[i-1] + step;
    if (fft_accurate_ == 1)
      sec.setBreakPoint (sampleTimes_[i]);
  }

  if(!(depSolVarIterVector_.empty())) // no point in calling this if depSolVarIterVector is empty
  {
    Util::Op::makeOps(comm, op_builder_manager, NetlistLocation(), depSolVarIterVector_.begin(), depSolVarIterVector_.end(), std::back_inserter(outputVars_));
  }
}

//-----------------------------------------------------------------------------
// Function      : FFTAnalysis::updateFFTData
// Purpose       : Called during the simulation to update the vectors with the
//                 time points and output variable values.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
void FFTAnalysis::updateFFTData(Parallel::Machine comm, const double circuitTime, const Linear::Vector *solnVec,
  const Linear::Vector *stateVec, const Linear::Vector * storeVec,
  const Linear::Vector *lead_current_vector, const Linear::Vector *junction_voltage_vector,
  const Linear::Vector *lead_current_dqdt_vector)
{
  // Save the time values, if interpolation is used.
  if (outputVars_.size() && (fft_accurate_ == 0))
    time_.push_back(circuitTime);

  int vecIndex = 0;
  for (std::vector<Util::Op::Operator *>::const_iterator it = outputVars_.begin(); it != outputVars_.end(); ++it)
  {
    double retVal = getValue(comm, *(*it), Util::Op::OpData(vecIndex, solnVec, 0, stateVec, storeVec, 0,
                             lead_current_vector, 0, junction_voltage_vector, 0)).real();
    if (fft_accurate_ == 0)
    {
      // save all of the output var values, if interpolation is used.
      outputVarsValues_.push_back(retVal);
    }
    else
    {
      // only save the values at the breakpoints, set at the specified sample times.
      if ( (sampleIdx_ < np_) && (abs( circuitTime - sampleTimes_[sampleIdx_]) <= sampleTimeTol_) )
      {
        sampleValues_[sampleIdx_] = retVal;
        sampleIdx_++;
      }
    }

    vecIndex++;
  }
}

//-----------------------------------------------------------------------------
// Function      : FFTAnalysis::outputResults
// Purpose       : Output results for this FFT Analysis at end of simulation
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
void FFTAnalysis::outputResults(std::ostream& outputStream)
{
  // Only calculate something if a .fft line was encountered and transient data was collected.
  int numOutVars = outputVars_.size();

  if ( (numOutVars>0) && !calculated_ )
  {
    // Calculate the fft coefficients for the given output variable.
    if (!time_.empty() && (fft_accurate_ == 0))
      interpolateData_();

    if (DEBUG_IO)
    {
      Xyce::dout() << std::endl << "Sample times and sampled/interpolated data values for FFT of " << outputVarName_ <<
                                 " are:" << std::endl;
      for (int i=0; i<sampleTimes_.size(); i++)
      {
        Xyce::dout() << "  " << sampleTimes_[i] << " , " << sampleValues_[i] << std::endl;
      }
      Xyce::dout() << std::endl;
    }

    applyWindowFunction_();
    calculateFFT_();
    calculated_ = true;
  }

  // Output the information to the outputStream
  printResult_( outputStream );
}

//-----------------------------------------------------------------------------
// Function      : FFTAnalysis::interpolateData_()
// Purpose       : evaluates interpolating polynomial at equidistant time pts
// Special Notes : In the final version of this class, this function will only
//                 be called if fft_accurate_ = 0.
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
bool FFTAnalysis::interpolateData_()
{
  if (time_.size() > 0)
  {
    Util::akima<double> interp;
    interp.init( time_, outputVarsValues_ );
    for (unsigned int i=0; i < np_; i++)
    {
      interp.eval( time_, outputVarsValues_, sampleTimes_[i], sampleValues_[i] );
     }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : FFTAnalysis::applyWindowFunction_()
// Purpose       : applies specified Windowing function to the interpolated
//                 data values.
// Special Notes :
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
bool FFTAnalysis::applyWindowFunction_()
{
  if (windowType_ == "RECT")
  {
    for (int i=0; i< np_; i++)
      ftInData_[i] = sampleValues_[i];
  }
  else if (windowType_ == "BART")
  {

    double weight;
    for (int i=0; i< np_; i++)
    {
      if (i < 0.5*(np_-1))
	weight = 2.0*i/(np_-1);
      else
        weight = 2.0 - 2.0*i/(np_-1);

      ftInData_[i] = sampleValues_[i]*weight;
    }
  }
  else if (windowType_ == "HANN")
  {
    for (int i=0; i< np_; i++)
      ftInData_[i] = sampleValues_[i]*sin(M_PI*i/(np_-1))*sin(M_PI*i/(np_-1));
  }
  else if (windowType_ == "HAMM")
  {
    for (int i=0; i< np_; i++)
      ftInData_[i] = sampleValues_[i] * (0.54 - 0.46*cos(2*M_PI*i/(np_-1)));
  }
  else if (windowType_ == "BLACK")
  {
    // "-67 dB Three-Term Blackman-Harris" window.  See SAND2017-4042.
    for (int i=0; i< np_; i++)
      ftInData_[i] = sampleValues_[i] *(0.42323 - 0.49755*cos(2*M_PI*i/(np_-1)) + 0.07922*cos(4*M_PI*i/(np_-1)));
  }
  else if (windowType_ == "HARRIS")
  {
    // "-92 dB Four-Term Blackman-Harris" window. See SAND2017-4042.
    for (int i=0; i< np_; i++)
      ftInData_[i] = sampleValues_[i] *(0.35875 - 0.48829*cos(2*M_PI*i/(np_-1)) + 0.14128*cos(4*M_PI*i/(np_-1))
					- 0.01168*cos(6*M_PI*i/(np_-1)));
  }
  else if ((windowType_=="GAUSS") || (windowType_=="KAISER"))
  {
    Report::UserWarning0() << "GAUSS and KAISER windows not supported yet. Defaulting to RECT";
  }

  if (DEBUG_IO)
  {
    Xyce::dout() << std::endl << "Sample times and windowed data values for FFT of " << outputVarName_
		              << " are:" << std::endl;
    for (int i=0; i<sampleTimes_.size(); i++)
    {
      Xyce::dout() << "  " << sampleTimes_[i] << " , " << ftInData_[i] << std::endl;
    }
    Xyce::dout() << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : FFTAnalysis::calculateFFT_()
// Purpose       : performs FFT analysis
// Special Notes :
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
void FFTAnalysis::calculateFFT_()
{
  ftInterface_->calculateFFT();

  fftRealCoeffs_[0] = ftOutData_[0]/np_;
  fftImagCoeffs_[0] = ftOutData_[1]/np_;
  double tmpVal;
  double noisePlusDist=0.0;
  double convRadDeg = 180.0/M_PI;

  // handle DC component
  mag_[0] = abs(fftRealCoeffs_[0]);
  phase_[0] = convRadDeg * atan2(fftImagCoeffs_[0], fftRealCoeffs_[0]);
  maxMag_ = mag_[0];

  // calculate mag, phase, total harmonic distortion (THD) and spurious free dynamic range (SFDR)
  for (int i=1;i<=np_/2; i++)
  {
    fftRealCoeffs_[i] = 2*ftOutData_[2*i]/np_;
    fftImagCoeffs_[i] = 2*ftOutData_[2*i+1]/np_;
    tmpVal = fftRealCoeffs_[i]*fftRealCoeffs_[i] + fftImagCoeffs_[i]*fftImagCoeffs_[i];
    mag_[i] = sqrt(tmpVal);
    if (mag_[i] > maxMag_)
      maxMag_ = mag_[i];
    if (fftout_ && (i>1))
      harmonicList_.push_back(std::make_pair(i,mag_[i])); // only make this vector, if it will be output

    // small magnitudes have phase set to 0
    if (mag_[i] > noiseFloor_)
      phase_[i] = convRadDeg * atan2(fftImagCoeffs_[i], fftRealCoeffs_[i]);

    if (i > 1)
    {
      noisePlusDist += tmpVal;
      thd_ += mag_[i]*mag_[i];
      if (mag_[i] > sfdr_)
      {
        sfdr_ = mag_[i];
        sfdrIndex_ = i;
      }
    }
  }

  // these metrics are output later if fftout_ =1
  sfdr_ = 20*log10(mag_[1]/sfdr_);                         // units are dB
  sndr_ = 20*log10(mag_[1] / sqrt(noisePlusDist)); // units are db
  enob_ = (sndr_ - 1.76)/6.02;                             // units are bits
  // don't take 20*log10() for THD, since both the actual value and dB value are output later.
  thd_ = sqrt(thd_)/mag_[1];

  // only sort the harmonicList_, if it will be output
  if (fftout_)
    std::sort(harmonicList_.begin(), harmonicList_.end(), fftMagCompFunc);

  if (DEBUG_IO)
  {
    std::string signStr;
    Xyce::dout() << "Coeffs from FFT interface and HSPICE-style coeffs for " << outputVarName_ << " are:" << std::endl;
    for (int i=0; i <= np_/2; i++)
    {
      signStr = (ftOutData_[2*i+1] < 0.0) ? "" : "+" ;
      Xyce::dout() << ftOutData_[2*i] << " " << signStr << ftOutData_[2*i+1] << "i     "
                   << fftRealCoeffs_[i] << " " << signStr << fftImagCoeffs_[i] << "i" <<std::endl;
    }
    Xyce::dout() << std::endl;
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : FFTAnalysis::printResult_( std::ostream& os )
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
std::ostream& FFTAnalysis::printResult_( std::ostream& os )
{
  basic_ios_all_saver<std::ostream::char_type> save(os);

  if (calculated_)
  {
    int colWidth1=12, colWidth2 = 16, precision = 6;
    std::string magString1, magString2;
    double normalization;

    // account for whether the output should be normalized or unnormalized values
    if (format_ == "NORM")
    {
      normalization = maxMag_;
      magString1 = "Norm. Mag (db)";
      magString2 = "Norm. Mag";
    }
    else
    {
      normalization = 1.0;
      magString1 = "Mag (db)";
      magString2 = "Mag";
    }

    os << "FFT analysis for " << outputVarName_ << ":" << std::endl
       << "  Window: " << windowType_ << std::scientific << std::setprecision(precision)
       << ", First Harmonic: " << freq_ <<", Start Freq: " << fmin_ << ", Stop Freq: " << fmax_ << std::endl;

    os << "DC component " << "   " << magString2 << "= " << mag_[0]/normalization
       << "   " << "Phase= " << phase_[0] << std::endl;

    os << std::setw(colWidth1) << "Index" << std::setw(colWidth2) << "Frequency"
       << std::setw(colWidth2) << std::setw(colWidth2) << magString2
       << std::setw(colWidth2) << "Phase" << std::endl;

    for (int i=1; i<=np_/2; i++)
    {
      os << std::setw(colWidth1) << i << std::setw(colWidth2) << i*freq_
         << std::setw(colWidth2) << mag_[i]/normalization
         << std::setw(colWidth2) << phase_[i] << std::endl;
    }

    if (fftout_)
    {
      os << std::endl
         << std::setw(colWidth1) << "THD = " << 20*log10(thd_) << " dB ( " << thd_ << " )" << std::endl
         << std::setw(colWidth1) << "SNDR = " << sndr_ << " dB" << std::endl
         << std::setw(colWidth1) << "ENOB = " << enob_ << " bit" << std::endl
         << std::setw(colWidth1) << "SFDR = " << sfdr_  << " dB at frequency " << sfdrIndex_*freq_ << std::endl;
    }

    os << std::endl;
  }

  return os;
}

} // namespace IO
} // namespace Xyce

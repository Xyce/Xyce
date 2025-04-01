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

//----------------------------------------------------------------
//
// Purpose       : This file is the class that implements an FFT analysis,
//                 as specified on a .FFT line.
//
// Special Notes :
//
// Creator       : Pete Sholander, SNL
//
// Creation Date : 1/3/2021
//---------------------------------------------------------------

#ifndef  Xyce_N_IO_FFTAnalysis_H
#define Xyce_N_IO_FFTAnalysis_H

#include <list>
#include <string>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_RCP.hpp>

#include <N_IO_fwd.h>
#include <N_IO_FFTMgr.h>

#include <N_LAS_Vector.h>
#include <N_PDS_fwd.h>
#include <N_TIA_StepErrorControl.h>

#include <N_UTL_fwd.h>
#include <N_UTL_FFTInterface.hpp>
#include <N_UTL_Param.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : FFTAnalysis
// Purpose       : Implements an FFT analysis, as specified on a .FFT line.
// Special Notes :
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 1/4/2021
//-----------------------------------------------------------------------------
class FFTAnalysis
{
public:
  FFTAnalysis(const Util::OptionBlock & fftBlock);

  // Destructor
  ~FFTAnalysis();

  // used to reset the FFTanalysis object at the start of each .STEP loop
  void reset();

  enum WindowType {RECT, BART, HANN, HAMM, BLACK, HARRIS, GAUSS, KAISER};

  // Return true if FFT analysis is being performed on any variables.
  bool isFFTActive() const { return !time_.empty(); }

  void fixupFFTParameters(Parallel::Machine comm,
                          const Util::Op::BuilderManager &op_builder_manager,
                          double endSimTime,
                          TimeIntg::StepErrorControl & sec,
                          bool fft_accurate,
                          bool fftout,
                          int fft_mode);

  void addSampleTimeBreakpoints();

  // Called during the simulation to update the fft objects held by this class
  void updateFFTData(Parallel::Machine comm,
                     double circuitTime,
                     const Linear::Vector *solnVec,
                     const Linear::Vector *stateVec,
                     const Linear::Vector * storeVec,
                     const Linear::Vector *lead_current_vector,
                     const Linear::Vector *junction_voltage_vector,
                     const Linear::Vector *lead_current_dqdt_vector);

  void outputResults( std::ostream& outputStream );
  void outputVerboseResults(std::ostream& outputStream);

  const Util::ParamList & getDepSolVarIterVector() const
  {
    return depSolVarIterVector_;
  }

  const std::string& getOutputVarName() const
  {
    return outputVarName_;
  }

  // getters used by .MEASURE FFT objects
  int getNP() const {return np_;}
  bool isCalculated() const {return calculated_;}

  // these functions are used by Measure FFT
  double getFundamentalFreq() const { return fundFreq_;}
  double getFirstHarmIdx() const { return fhIdx_;}
  double getNoiseFloor() const {return noiseFloor_;}
  const std::vector<double>& getMagVec() const {return mag_;}
  double getFFTCoeffRealVal(int index) const { return fftRealCoeffs_[index];}
  double getFFTCoeffImagVal(int index) const { return fftImagCoeffs_[index];}
  double getMagVal(int index) const { return mag_[index];}
  double getPhaseVal(int index) const { return phase_[index];}
  double calculateENOBforMeasFFT(int binSize) const;
  double calculateSNDRforMeasFFT(int binSize) const;
  double calculateSFDRforMeasFFT(int fminIndex, int fmaxIndex, bool fminGivn, int binSize) const;

  // these functions are used by both FFTAnalysis and Measure FFT
  double calculateSNR(int fmaxIndex) const;
  double calculateTHD(int fmaxIndex) const;
  double convertValuetoDB(double val) const;

private:
  void calculateResults_();
  bool interpolateData_();
  bool applyWindowFunction_();

  void calculateFFT_();
  void calculateSFDR_();
  int getLowerLimforSFDR_(int fminIndex, int fmaxIndex, bool fminGivn) const;
  void calculateSNDRandENOB_();

  std::ostream& printResult_( std::ostream& os );
  std::ostream& printVerboseResult_( std::ostream& os );
  std::ostream& printResultHeader_( std::ostream& os );
  std::ostream& printMetrics_( std::ostream& os );

  // used to sort a vector (in descending order) based on the double value in the pairs
  static bool fftMagCompFunc(const std::pair<int,double>& a, const std::pair<int,double>& b)
  {
    return a.second > b.second;
  }

private:
  TimeIntg::StepErrorControl* secPtr_; // ptr to step error control.  Used if FFT_ACCURATE=true
  double startTime_, stopTime_;
  int np_;
  std::string format_;
  bool formatGiven_;
  std::string windowType_;
  double alpha_;
  double fundFreq_; // fundamental frequency
  double freq_, fmin_, fmax_;  // values from FREQ (first harmonic), FMIN, FMAX qualifiers
  int fhIdx_, fminIdx_, fmaxIdx_; // values rounded to nearest harmonic index
  bool startTimeGiven_;
  bool stopTimeGiven_;
  bool freqGiven_;
  bool fminGiven_;
  bool fmaxGiven_;
  bool calculated_;
  std::string outputVarName_;
  bool fft_accurate_;
  bool fftout_;
  int fft_mode_;
  int sampleIdx_;
  double noiseFloor_;
  double maxMag_;
  double normalization_;  // used for outputting in NORM vs. UNORM format
  double thd_;
  double sndr_;
  double enob_;
  double snr_;
  double sfdr_;
  int sfdrIndex_;
  int colWidth1_;         // used for output formatting 
  int colWidth2_; 
  int precision_;
  std::vector<double> mag_;
  std::vector<double> phase_;
  std::vector<std::pair<int,double>> harmonicList_;

  int numDepSolVars_;

  // Xyce FFT interface
  Teuchos::RCP<N_UTL_FFTInterface<std::vector<double> > > ftInterface_;
  // these are FFT coefficients as reported by the Xyce FFT interface.
  std::vector<double> ftInData_, ftOutData_, iftInData_, iftOutData_;

  // the real and imaginary parts of the FFT coeffs, translated into HSPICE format
  std::vector<double> fftRealCoeffs_;
  std::vector<double> fftImagCoeffs_;

  std::vector<double> time_;
  std::vector<double> outputVarValues_;
  Util::ParamList depSolVarIterVector_;
  Util::Op::OpList outputVars_;
  std::vector<double> sampleTimes_, sampleValues_;
};

} // namespace IO
} // namespace Xyce

#endif  // Xyce_N_IO_FFTAnalysis_H

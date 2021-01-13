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

  enum WindowType {RECT, BART, HANN, HAMM, BLACK, HARRIS, GAUSS, KAISER};

  // Return true if FFT analysis is being performed on any variables.
  bool isFFTActive() const { return !time_.empty(); }

  void fixupFFTParameters(Parallel::Machine comm,
                          const Util::Op::BuilderManager &op_builder_manager,
                          const double endSimTime,
                          TimeIntg::StepErrorControl & sec,
                          const int fft_accurate,
                          const bool fftout);

  // Called during the simulation to update the fft objects held by this class
  void updateFFTData(Parallel::Machine comm,
                     const double circuitTime,
                     const Linear::Vector *solnVec,
                     const Linear::Vector *stateVec,
                     const Linear::Vector * storeVec,
                     const Linear::Vector *lead_current_vector,
                     const Linear::Vector *junction_voltage_vector,
                     const Linear::Vector *lead_current_dqdt_vector);

  void outputResults( std::ostream& outputStream );

  const Util::ParamList & getDepSolVarIterVector() const
  {
    return depSolVarIterVector_;
  }

private:
  bool interpolateData_();
  bool applyWindowFunction_();

  void calculateFFT_();

  std::ostream& printResult_( std::ostream& os );

private:
  double startTime_, stopTime_;
  int np_;
  std:: string format_;
  std::string windowType_;
  double alpha_;
  double freq_;
  double fmin_;
  double fmax_;
  bool stopTimeGiven_;
  bool freqGiven_;
  bool fminGiven_;
  bool fmaxGiven_;
  bool calculated_;
  std::string outputVarName_;
  int fft_accurate_;
  int fftout_;
  int sampleIdx_;
  double sampleTimeTol_;
  double thd_;
  double sndr_;
  double enob_;
  double sfdr_;
  int sfdrIndex_;
  std::vector<double> mag_;
  std::vector<double> normMag_;
  std::vector<double> phase_;

  int numDepSolVars_;

  Teuchos::RCP<N_UTL_FFTInterface<std::vector<double> > > ftInterface_;
  std::vector<double> ftInData_, ftOutData_, iftInData_, iftOutData_;

  std::vector<double> time_;
  std::vector<double> outputVarsValues_;
  Util::ParamList depSolVarIterVector_;
  Util::Op::OpList outputVars_;
  std::vector<double> sampleTimes_, sampleValues_;
};

} // namespace IO
} // namespace Xyce

#endif  // Xyce_N_IO_FFTAnalysis_H

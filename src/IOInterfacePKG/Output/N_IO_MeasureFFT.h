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
//
// Purpose        : Support various FFT-based measures.
//
// Special Notes  :
//
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 01/18/2021
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureFFT_h
#define Xyce_N_IO_MeasureFFT_h

#include <N_IO_MeasureBase.h>
#include <N_IO_FFTAnalysis.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-------------------------------------------------------------------------
// Class         : FFT
// Purpose       : Support various FFT-based measures
// Special Notes : This is the base class.
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-------------------------------------------------------------------------
class FFT : public Base
{
public:
  FFT(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~FFT() {};

  void prepareOutputVariables();
  void resetFFT();

  void updateTran(
    Parallel::Machine comm,
    const double circuitTime,
    const Linear::Vector *solnVec,
    const Linear::Vector *stateVec,
    const Linear::Vector *storeVec,
    const Linear::Vector *lead_current_vector,
    const Linear::Vector *junction_voltage_vector,
    const Linear::Vector *lead_current_dqdt_vector);

  void updateDC(
    Parallel::Machine comm,
    const std::vector<Analysis::SweepParam> & dcParamsVec,
    const Linear::Vector *solnVec,
    const Linear::Vector *stateVec,
    const Linear::Vector *storeVec,
    const Linear::Vector *lead_current_vector,
    const Linear::Vector *junction_voltage_vector,
    const Linear::Vector *lead_current_dqdt_vector);

  void updateAC(
    Parallel::Machine comm,
    const double frequency,
    const Linear::Vector *solnVec,
    const Linear::Vector *imaginaryVec,
    const Util::Op::RFparamsData *RFparams);

  void updateNoise(
    Parallel::Machine comm,
    const double frequency,
    const Linear::Vector *solnVec,
    const Linear::Vector *imaginaryVec,
    const double totalOutputNoiseDens,
    const double totalInputNoiseDens,
    const std::vector<Xyce::Analysis::NoiseData*> *noiseDataVec);

  virtual double getMeasureResult()=0;

  // used to pass info from associated FFTAnalysis object to this measure object
  void fixupFFTMeasure(FFTAnalysis* fftAnalysisPtr);

protected:
  FFTAnalysis* fftAnalysisPtr_;
  int np_;

private:
  int numOutVars_;
  std::vector<double> outVarValues_;
};

//-------------------------------------------------------------------------
// Class         : FFTFind
// Purpose       : Implement FIND measure for .MEASURE FFT
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-------------------------------------------------------------------------
class FFTFind : public FFT
{
public:
  FFTFind(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~FFTFind() {};

public:
  void reset();
  double getMeasureResult();

};

//-------------------------------------------------------------------------
// Class         : ENOB
// Purpose       : Implement ENOB measure for .MEASURE FFT
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-------------------------------------------------------------------------
class ENOB : public FFT
{
public:
  ENOB(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~ENOB() {};

public:
  void reset();
  double getMeasureResult();

};

//-------------------------------------------------------------------------
// Class         : SFDR
// Purpose       : Implement SFDR measure for .MEASURE FFT
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-------------------------------------------------------------------------
class SFDR : public FFT
{
public:
  SFDR(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~SFDR() {};

public:
  void reset();
  double getMeasureResult();

};

//-------------------------------------------------------------------------
// Class         : SNDR
// Purpose       : Implement SNDR measure for .MEASURE FFT
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-------------------------------------------------------------------------
class SNDR : public FFT
{
public:
  SNDR(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~SNDR() {};

public:
  void reset();
  double getMeasureResult();

};

//-------------------------------------------------------------------------
// Class         : THD
// Purpose       : Implement THD measure for .MEASURE FFT
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-------------------------------------------------------------------------
class THD : public FFT
{
public:
  THD(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~THD() {};

public:
  void reset();
  double getMeasureResult();
  std::ostream& printVerboseMeasureResult(std::ostream& os);
};

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif

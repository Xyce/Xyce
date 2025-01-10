//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : The outputter class for .PRINT HB_FD for FORMAT=STD
//
// Special Notes  :
//
// Creator        : Pete Sholander, SNL, Electrical Models and Simulation
//
// Creation Date  : 6/25/2018
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputterHBFreqPrn_h
#define Xyce_N_IO_OutputterHBFreqPrn_h

#include <N_IO_OutputterLocal.h>

namespace Xyce {
namespace IO {
namespace Outputter {

class HBFreqPrn : public HBInterface
{
public:
  HBFreqPrn(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &freq_print_parameters);

  virtual ~HBFreqPrn();

private:
  HBFreqPrn(const HBFreqPrn &);
  HBFreqPrn &operator=(const HBFreqPrn &);

public:
  virtual void doSetAnalysisMode(Analysis::Mode analysis_mode)
  {
    freqPrintParameters_.analysisMode_ = analysis_mode;
  }

  virtual void doOutputHB_FD(
    Parallel::Machine           comm,
    const std::vector<double> & freqPoints,
    const Linear::BlockVector & freqDomainSolutionVecReal,
    const Linear::BlockVector & freqDomainSolutionVecImaginary,
    const Linear::BlockVector & freqDomainLeadCurrentVecReal,
    const Linear::BlockVector & freqDomainLeadCurrentVecImaginary,
    const Linear::BlockVector & freqDomainJunctionVoltageVecReal,
    const Linear::BlockVector & freqDomainJunctionVoltageVecImaginary);

  virtual void doFinishOutput();

  virtual void doStartStep(int current_step, int number_of_step);

  virtual void doResetIndex();

  virtual void doSteppingComplete();

private:
  OutputMgr &           outputManager_;
  PrintParameters       freqPrintParameters_;
  std::string           freqFilename_;
  std::ostream *        fos_;
  int                   index_;
  int                   currentStep_;
  int                   numberOfSteps_;

  Util::Op::OpList      freqOpList_;
};


} // namespace Outputter
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputterHBFreqPrn_h

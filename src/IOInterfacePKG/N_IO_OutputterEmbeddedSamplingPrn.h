//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose        : Outputter class for prn files for Embedded Sampling info
//
// Special Notes  :
//
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 7/26/2019
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputterEmbeddedSamplingPrn_h
#define Xyce_N_IO_OutputterEmbeddedSamplingPrn_h

#include <N_IO_OutputterLocal.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Embedded Sampling outputters

class EmbeddedSamplingPrn : public Interface
{
public:
  EmbeddedSamplingPrn(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters);

  virtual ~EmbeddedSamplingPrn();

private:
  EmbeddedSamplingPrn(const EmbeddedSamplingPrn &);
  EmbeddedSamplingPrn &operator=(const EmbeddedSamplingPrn &);

public:

  virtual void doSetAnalysisMode(Analysis::Mode analysis_mode)
  {
    printParameters_.analysisMode_ = analysis_mode;
  }

  virtual void doFinishOutput();

  virtual void doStartStep(int step, int max_step);

  virtual void doResetIndex();

  virtual void doSteppingComplete();

  virtual void doOutputEmbeddedSampling(
    Parallel::Machine           comm,
    bool                        regressionPCEenable,
    bool                        projectionPCEenable,
    int                         numSamples,
    const std::vector<std::string> & regressionPCEcoeffs,
    const std::vector<std::string> & projectionPCEcoeffs,
    const std::vector<Xyce::Analysis::UQ::outputFunctionData*> & outFuncDataVec);

private:
  void EmbeddedSamplingHeader();

private:
  OutputMgr &                   outputManager_;
  PrintParameters               printParameters_;
  std::string                   outFilename_;
  std::ostream *                os_;
  int                           printCount_;
  int                           index_;
  int                           currentStep_;
  int                           numberOfSteps_;

  Util::Op::OpList              opList_;
};

} // namespace Outputter
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputterEmbeddedSamplingPrn_h

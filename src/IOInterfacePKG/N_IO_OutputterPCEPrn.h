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

//-----------------------------------------------------------------------------
//
// Purpose        : Outputter class for prn files for PCE info
//
// Special Notes  :
//
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 9/3/2019
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputterPCEPrn_h
#define Xyce_N_IO_OutputterPCEPrn_h

#include <N_IO_OutputterLocal.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// PCE outputters

class PCEPrn : public Interface
{
public:
  PCEPrn(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters);

  virtual ~PCEPrn();

private:
  PCEPrn(const PCEPrn &);
  PCEPrn &operator=(const PCEPrn &);

public:

  virtual void doSetAnalysisMode(Analysis::Mode analysis_mode)
  {
    printParameters_.analysisMode_ = analysis_mode;
  }

  virtual void doFinishOutput();

  virtual void doStartStep(int step, int max_step);

  virtual void doResetIndex();

  virtual void doSteppingComplete();

  virtual void doOutputPCE(
    Parallel::Machine           comm,
    int                         numQuadPoints,
    const std::vector<Xyce::Analysis::UQ::outputFunctionData*> & outFuncDataVec);

private:
  void PCEHeader();

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

#endif // Xyce_N_IO_OutputterPCEPrn_h

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
// Purpose        : The outputter class for .PRINT HB_TD for FORMAT=STD
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

#ifndef Xyce_N_IO_OutputterHBTimePrn_h
#define Xyce_N_IO_OutputterHBTimePrn_h

#include <N_IO_OutputterLocal.h>

namespace Xyce {
namespace IO {
namespace Outputter {

class HBTimePrn : public HBInterface
{
public:
  HBTimePrn(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &time_print_parameters);

  virtual ~HBTimePrn();

private:
  HBTimePrn(const HBTimePrn &);
  HBTimePrn &operator=(const HBTimePrn &);

public:
  virtual void doSetAnalysisMode(Analysis::Mode analysis_mode)
  {
    timePrintParameters_.analysisMode_ = analysis_mode;
  }

  virtual void doOutputHB_TD(
    Parallel::Machine           comm,
    const std::vector<double> & timePoints,
    const Linear::BlockVector & timeDomainSolutionVec,
    const Linear::BlockVector & timeDomainLeadCurrentVec,
    const Linear::BlockVector & timeDomainJunctionVoltageVec);

  virtual void doFinishOutput();

  virtual void doStartStep(int current_step, int number_of_step);

  virtual void doResetIndex();

  virtual void doSteppingComplete();

private:
  OutputMgr &           outputManager_;
  PrintParameters       timePrintParameters_;
  std::string           timeFilename_;
  std::ostream *        tos_;
  int                   index_;
  int                   currentStep_;
  int                   numberOfSteps_;

  Util::Op::OpList      timeOpList_;
};


} // namespace Outputter
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputterHBTimePrn_h

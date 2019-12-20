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
// Purpose        : Outputter class for csv files for HB_IC info
//
// Special Notes  :
//
// Creator        : Pete Sholander
//
// Creation Date  : 11/19/2019
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputterHBICCSV_h
#define Xyce_N_IO_OutputterHBICCSV_h

#include <N_IO_OutputterLocal.h>

namespace Xyce {
namespace IO {
namespace Outputter {


class HBICCSV : public TimeInterface
{
public:
  HBICCSV(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters);

  virtual ~HBICCSV();

private:
  HBICCSV(const HBICCSV &);
  HBICCSV &operator=(const HBICCSV &);

public:

  virtual void doSetAnalysisMode(Analysis::Mode analysis_mode)
  {
    printParameters_.analysisMode_ = analysis_mode;
  }

  virtual void doOutputTime(
    Parallel::Machine           comm,
    const Linear::Vector &        solution_vector,
    const Linear::Vector &        state_vector,
    const Linear::Vector &        store_vector,
    const Linear::Vector &  lead_current_vector,
    const Linear::Vector &  junction_voltage_vector) ;

  virtual void doFinishOutput();

  virtual void doStartStep(int current_step, int number_of_step);

  virtual void doResetIndex();

  virtual void doSteppingComplete();

  virtual void reopenTmpFile();

  virtual void copyTmpFileToOutputFile();

private:
  OutputMgr &           outputManager_;
  PrintParameters       printParameters_;
  std::string           outFilename_;
  std::ostream *        os_;  // used for normal output file
  std::ostream *        tmpOutStream_; // used for tmp output file, generated when
                                       // transient-assisted HB is used with .STEP
  int                   index_;
  int                   currentStep_;
  int                   numberOfSteps_;

  Util::Op::OpList      opList_;
};

} // namespace Outputter
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputterHBICCSV_h

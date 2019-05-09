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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Baur
//
// Creation Date  :
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputterTimeRaw_h
#define Xyce_N_IO_OutputterTimeRaw_h

#include <N_IO_OutputterLocal.h>

namespace Xyce {
namespace IO {
namespace Outputter {

class TimeRaw : public TimeInterface
{
public:
  TimeRaw(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters);

  virtual ~TimeRaw();

private:
  TimeRaw(const TimeRaw &);
  TimeRaw &operator=(const TimeRaw &);

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
  const Linear::Vector &  junction_voltage_vector);

  virtual void doFinishOutput();

private:
  void timeHeader();

private:
  OutputMgr &           outputManager_;
  PrintParameters       printParameters_;
  std::string           outFilename_;
  int                   numPoints_;
  long                  numPointsPos_;
  std::ostream *        os_;
  bool                  outputRAWTitleAndDate_;

  Util::Op::OpList      opList_;
};

} // namespace Outputter

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputterTimeRaw_h

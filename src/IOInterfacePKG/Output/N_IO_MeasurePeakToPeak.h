//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Purpose        : Measure statistics of a simulation variable
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 03/10/2009
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasurePeakToPeak_h
#define Xyce_N_IO_MeasurePeakToPeak_h

#include <N_IO_MeasureBase.h>
#include <N_IO_MeasureExtrema.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-------------------------------------------------------------------------
// Class         : PeakToPeak
// Purpose       : Measure the peak-to-peak value of a simulation variable
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-------------------------------------------------------------------------
class PeakToPeak : public Extrema
{
public:
  PeakToPeak(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~PeakToPeak() {};

  void reset();

  double getMeasureResult();
  std::ostream& printMeasureResult(std::ostream& os);
  std::ostream& printVerboseMeasureResult(std::ostream& os);

  void setMeasureVarsForNewWindow(double indepVarVal, double depVarVal);
  void updateMeasureVars(double indepVarVal, double depVarVal);

private:
  std::string type_;
  double maximumValue_;
  double minimumValue_;

  // maximumInstant_ will be a time value for TRAN mode, a frequency value for AC mode 
  // and a Sweep Variable value for DC mode.
  double maximumInstant_;
  double minimumInstant_;
};

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif

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
// Purpose        : Find the extrema (max, min or peak-to-peak value)
//                  of a simulation variable
//
// Special Notes  : This class contains the common functions for the Max, Min and
//                  PeakToPeak classes, and sits between those classes and the Base
//                  class. In general, it applies to measures that are based on the
//                  current value of the simulation variable.  For TRAN measures,
//                  the FROM, TO, TD, RISE, FALL and CROSS qualifiers can be used
//                  to limit the measurement window.  For AC and DC measures, the
//                  FROM and TO qualifiers can be used to limit the measurement
//                  window.
//
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 04/28/2020
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureExtrema_h
#define Xyce_N_IO_MeasureExtrema_h

#include <N_IO_MeasureBase.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-------------------------------------------------------------------------
// Class         : Extrema
// Purpose       : Find the extrema (max, min or peak-to-peak value)
//                 of a simulation variable
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 04/28/2010
//-------------------------------------------------------------------------
class Extrema : public Base
{
public:
  Extrema(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~Extrema() {};

  void prepareOutputVariables();
  void resetExtrema();

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

  virtual double getMeasureResult()=0;
  virtual std::ostream& printMeasureResult(std::ostream& os, bool printVerbose=false)=0;

  virtual void setMeasureVarsForNewWindow(const double indepVarVal, const double depVarVal)=0;
  virtual void updateMeasureVars(const double indepVarVal, const double depVarVal)=0;

private:
  int numOutVars_;
  std::vector<double> outVarValues_;
};

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif

//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Purpose        : Measure statistics of a simulation variable over
//                  an interval
//
// Special Notes  : This class contains the common functions for the Average,
//                  IntegralEvaluation and RMS classes, and sits between those
//                  classes and the Base class.  In general, it applies to
//                  measures that are calculated based on a function of the current
//                  and previous signal value, as calculated over the entire measurement
//                  window.  For TRAN measures, the FROM, TO, TD, RISE, FALL and
//                  CROSS qualifiers can be used to limit the measurement window.
//                  For AC and DC measures, the FROM and TO qualifiers can be used
//                  to limit the measurement window.

//
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 04/28/2020
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureStats_h
#define Xyce_N_IO_MeasureStats_h

#include <N_IO_MeasureBase.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-------------------------------------------------------------------------
// Class         : Stats
// Purpose       : Measure statistics of a simulation variable over
//                 an interval
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 04/29/2020
//-------------------------------------------------------------------------
class Stats : public Base
{
public:
  Stats(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~Stats() {};

  void prepareOutputVariables();
  void resetStats();

  void updateTran(
    Parallel::Machine comm,
    double circuitTime,
    double endSimTime,
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
    double frequency,
    double fStart,
    double fStop,
    const Linear::Vector *solnVec,
    const Linear::Vector *imaginaryVec,
    const Util::Op::RFparamsData *RFparams);

  void updateNoise(
    Parallel::Machine comm,
    double frequency,
    double fStart,
    double fStop,
    const Linear::Vector *solnVec,
    const Linear::Vector *imaginaryVec,
    double totalOutputNoiseDens,
    double totalInputNoiseDens,
    const std::vector<Xyce::Analysis::NoiseData*> *noiseDataVec);

  virtual double getMeasureResult()=0;

protected:
  // this function is only used for TRAN mode
  virtual void setMeasureVarsForNewWindow()=0;

  // this function is defined for AC, DC and TRAN modes
  virtual void updateMeasureVars(double indepVarVal, double depVarVal)=0;

  double lastIndepVarValue_;
  double lastSignalValue_;
  int numPointsFound_;

private:
  void updateMeasureState_(double indepVarVal, double depVarVal);

  int numOutVars_;
  std::vector<double> outVarValues_;
};

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif

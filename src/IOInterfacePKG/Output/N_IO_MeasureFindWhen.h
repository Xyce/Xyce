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
// Purpose        : Find the time when a variable hits a target value (WHEN measure),
//                  or find the value of a variable at the time that the WHEN
//                  clause is satisfied (FIND-WHEN measure).
//
// Special Notes  :
//
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 09/27/2021
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureFindWhen_h
#define Xyce_N_IO_MeasureFindWhen_h

#include <N_IO_MeasureWhenAT.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-------------------------------------------------------------------------
// Class         : FindWhenBase
// Purpose       : Find time when a variable hits a target value, or find
//                 the value of a variable at the time that the WHEN clause is
//                 satisfied. This class contains the variables and functions
//                 that are common to both the continous and non-continuous
//                 versions of this measure.  The non-continuous version
//                 (e.g., .MEASURE TRAN) will only return one measure value.
//                 The continuous version (e.g., .MEASURE TRAN_CONT) may
//                 return multiple values.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 08/03/2020
//-------------------------------------------------------------------------
class FindWhenBase : public WhenAT
{
public:
  FindWhenBase(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~FindWhenBase() {};

  void prepareOutputVariables();

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

protected:
  bool checkMeasureLine() const;

private:
  void updateMeasureVarsForAT(double currIndepVarVal);
  void updateMeasureVarsForWhen(double currIndepVarVal, double targVal, double whenInstant);
  double interpolateFindValue(double currIndepVarValue, double targVal, double whenTime) const;
};

//-------------------------------------------------------------------------
// Class         : FindWhen
// Purpose       : Non-continuous version of FindWhen, that returns only
//                 only one measure value.  This class is invoked by
//                 .MEASURE AC, .MEASURE DC, .MEASURE NOISE and .MEASURE TRAN.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 08/03/2020
//-------------------------------------------------------------------------
class FindWhen : public FindWhenBase
{
public:
  FindWhen(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~FindWhen() {};

  void reset();

  std::ostream& printMeasureResult(std::ostream& os);
  std::ostream& printVerboseMeasureResult(std::ostream& os);
};

//-------------------------------------------------------------------------
// Class         : FindWhenCont
// Purpose       : Continuous version of FindWhen that can return multiple
//                 measure values if RISE, FALL or CROSS is not specified.
//                 This class is invoked by .MEASURE AC_CONT, .MEASURE DC_CONT,
//                 .MEASURE NOISE_CONT and .MEASURE TRAN_CONT.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 08/03/2020
//-------------------------------------------------------------------------
class FindWhenCont : public FindWhenBase
{
public:
  FindWhenCont(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~FindWhenCont() {};

  void reset();

  std::ostream& printMeasureResult(std::ostream& os);
  std::ostream& printVerboseMeasureResult(std::ostream& os);
};

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif

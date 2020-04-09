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
// Purpose        : Implement the ERR, ERR1 and ERR2 measure types
//
// Special Notes  : The ERR and ERR1 measure types are both implemented by
//                  the Err1 class
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 03/08/2020
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureErrorFunctions_h
#define Xyce_N_IO_MeasureErrorFunctions_h

#include <N_IO_MeasureBase.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-------------------------------------------------------------------------
// Class         : ErrorFunctions
// Purpose       : Implement the ERRx measures
// Special Notes : This class contains the functions that are common to
//                 the Err1 and Err2 measure classes
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-------------------------------------------------------------------------
class ErrorFunctions : public Base
{
public:
  ErrorFunctions(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~ErrorFunctions() {};

  void prepareOutputVariables();
  bool checkMeasureLine();
  void resetErrorFunctions();

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

  virtual double getMeasureResult(){return 0.0;};
  virtual void updateErrVars(double mVal, double cVal){};

  bool withinYLimits(double mVal); // enforce YMAX and YMIN

private:
  std::string type_;
  int numOutVars_;
  std::vector<double> outVarValues_;
};

//-------------------------------------------------------------------------
// Class         : Err1
// Purpose       : Implement ERR1 measure
// Special Notes : This class also implements the ERR measure
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-------------------------------------------------------------------------
class Err1 : public ErrorFunctions
{
public:
  Err1(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~Err1() {};

public:
  void reset();
  void updateErrVars(double mVal, double cVal);
  double getMeasureResult();

private:
  double err1SqSum_;
  int numPts_;
};

//-------------------------------------------------------------------------
// Class         : Err2
// Purpose       : Implement ERR2 measure
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 03/08/2020
//-------------------------------------------------------------------------
class Err2 : public ErrorFunctions
{
public:
  Err2(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~Err2() {};

public:
  void reset();
  void updateErrVars(double mVal, double cVal);
  double getMeasureResult();

private:
  double err2Sum_;
  int numPts_;
};

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif

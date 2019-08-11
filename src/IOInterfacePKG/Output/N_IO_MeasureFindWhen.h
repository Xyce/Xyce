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
// Purpose        : Find time when a variable hits a target value
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 03/10/2009
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureFindWhen_h
#define Xyce_N_IO_MeasureFindWhen_h

#include <N_IO_MeasureBase.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-------------------------------------------------------------------------
// Class         : FindWhen
// Purpose       : Find time when a variable hits a target value
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-------------------------------------------------------------------------
class FindWhen : public Base
{
public:
  FindWhen(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~FindWhen() {};

  void prepareOutputVariables();
  bool checkMeasureLine();
  void reset();

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

  std::ostream& printMeasureResult(std::ostream& os, bool printVerbose=false);
  std::ostream& printRFCWindow(std::ostream& os);

private:
  int numOutVars_;
  std::vector<double> outVarValues_;

  // these are used to interpolate to the target independent var (time, frequency etc)  when 
  // the simulation has reported values that bound the target value
  double lastIndepVarValue_;
  double lastDepVarValue_;
  double lastOutputVarValue_;

  // This variable controls what is tested against in the WHEN clause.  
  // It refers to an index in the outputVarValues_ array.  Its value 
  // depends on whether the measure is a FIND-WHEN or a WHEN measure.
  // For a WHEN measure, whenIdx_ = 0 which is the default.  For a FIND-WHEN 
  // measure, whenIdx_ = 1  For the syntax WHEN v(a)=v(b), the value of v(a) 
  // is in outputVarValues_[whenIdx_] and the value of v(b) is in 
  // outputVarValues_[whenIdx_+1]
  int whenIdx_;

};

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif

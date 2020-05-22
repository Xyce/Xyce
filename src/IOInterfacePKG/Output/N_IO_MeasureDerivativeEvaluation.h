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
// Purpose        : Measure the derivative of an output variable
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 03/10/2009
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureDerivativeEvaluation_h
#define Xyce_N_IO_MeasureDerivativeEvaluation_h

#include <N_IO_MeasureBase.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-------------------------------------------------------------------------
// Class         : DerivativeEvaluation
// Purpose       : Measure the derivative of an output variable
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-------------------------------------------------------------------------
class DerivativeEvaluation : public Base
{
public:
  DerivativeEvaluation(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~DerivativeEvaluation() {};

  void prepareOutputVariables();
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

  void updateNoise(
    Parallel::Machine comm,
    const double frequency,
    const Linear::Vector *solnVec,
    const Linear::Vector *imaginaryVec,
    const double totalOutputNoiseDens,
    const double totalInputNoiseDens,
    const std::vector<Xyce::Analysis::NoiseData*> *noiseDataVec);

  double getMeasureResult();
  std::ostream& printMeasureWindow(std::ostream& os, const double endSimTime);
  std::ostream& printMeasureResult(std::ostream& os);
  std::ostream& printVerboseMeasureResult(std::ostream& os);
  std::ostream& printRFCWindow(std::ostream& os);

  void updateCalculationResult(const double indepVarVal);
  void processATforACDCNoise(const double indepVarVal);
  void processWHENforACDCNoise(const double indepVarVal);

  void setMeasureState(const double indepVarVal);
  void updateMeasureState(const double indepVarVal);
  void updateTargVal(double& targVal);

  // used for debugging purpose
  void interpolateCalculationInstant(double currIndepVarValue, double targVal);

private:
  std::string type_;
  int numOutVars_;
  std::vector<double> outVarValues_;

  double lastIndepVarValue_;
  double lastDepVarValue_;
  double lastOutputVarValue_;
  double lastTargValue_;
  int    numPointsFound_;  // AC and DC calculations need at least two points to be valid

  // This variable controls what is tested against in the WHEN clause.  
  // It refers to an index in the outputVarValues_ array.  For the syntax
  // WHEN v(a)=<val> the value of v(a) is in outputVarValues_[whenIdx_]
  // For the syntax WHEN v(a)=v(b) the value of v(a) is in outputVarValues_[whenIdx_]
  // and the value of v(b) is in outputVarValues_[whenIdx_+1]
  int whenIdx_;
};

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif

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
// Class         : DerivativeEvaluationBase
// Purpose       : Measure the derivative of an output variable.  This class
//                 contains the variables and functions that are common to
//                 both the continous and non-continuous versions of this
//                 measure.  The non-continuous version (e.g., .MEASURE TRAN)
//                 will only return one measure value.  The continuous version
//                 (e.g., .MEASURE TRAN_CONT) may return multiple values.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-------------------------------------------------------------------------
class DerivativeEvaluationBase : public Base
{
public:
  DerivativeEvaluationBase(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~DerivativeEvaluationBase() {};

  void prepareOutputVariables();
  void resetDerivativeEvaluationBase();
  void updateTran(
    Parallel::Machine comm,
    const double circuitTime,
    const double endSimTime,
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
    const double fStart,
    const double fStop,
    const Linear::Vector *solnVec,
    const Linear::Vector *imaginaryVec,
    const Util::Op::RFparamsData *RFparams);

  void updateNoise(
    Parallel::Machine comm,
    const double frequency,
    const double fStart,
    const double fStop,
    const Linear::Vector *solnVec,
    const Linear::Vector *imaginaryVec,
    const double totalOutputNoiseDens,
    const double totalInputNoiseDens,
    const std::vector<Xyce::Analysis::NoiseData*> *noiseDataVec);

  double getMeasureResult();
  std::ostream& printMeasureWindow(std::ostream& os, const double endSimTime,
				   const double startSweepVal, const double endSweepVal);
  std::ostream& printRFCWindow(std::ostream& os);

  virtual void updateCalculationResult(double val)=0;
  virtual void updateCalculationInstant(double val)=0;

  bool isATforACDCNoise(const double indepVarVal);
  bool isWHENcondition(const double indepVarVal, const double targVal);

  void setMeasureState(const double indepVarVal);
  void updateMeasureState(const double indepVarVal);
  double updateTargVal();
  double getDerivativeValue(const double currIndepVarValue);

protected:
  std::vector<double> calculationResultVec_;
  std::vector<double> calculationInstantVec_;

private:
  void updateMeasureVars(const double currIndepVarVal, const double targVal, const double whenInstant);
  double interpolateCalculationInstant(double currIndepVarValue, double targVal);

  void updateRFCcountForWhen();
  bool withinRFCWindowForWhen();

  std::string type_;
  int numOutVars_;
  std::vector<double> outVarValues_;

  double lastIndepVarValue_;
  double lastDepVarValue_;
  double lastOutputVarValue_;
  double lastTargValue_;
  double startDCMeasureWindow_;
  int    numPointsFound_;  // AC and DC calculations need at least two points to be valid

  // This variable controls what is tested against in the WHEN clause.  
  // It refers to an index in the outputVarValues_ array.  For the syntax
  // WHEN v(a)=<val> the value of v(a) is in outputVarValues_[whenIdx_]
  // For the syntax WHEN v(a)=v(b) the value of v(a) is in outputVarValues_[whenIdx_]
  // and the value of v(b) is in outputVarValues_[whenIdx_+1]
  int whenIdx_;
};

//-------------------------------------------------------------------------
// Class         : DerivativeEvaluation
// Purpose       : Non-continuous version of DERIV measure, that returns only
//                 only one measure value.  This class is invoked by
//                 .MEASURE AC, .MEASURE DC, .MEASURE NOISE and .MEASURE TRAN.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 08/21/2020
//-------------------------------------------------------------------------
class DerivativeEvaluation : public DerivativeEvaluationBase
{
public:
  DerivativeEvaluation(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~DerivativeEvaluation() {};

  void reset();

  void updateCalculationResult(double val);
  void updateCalculationInstant(double val);

  std::ostream& printMeasureResult(std::ostream& os);
  std::ostream& printVerboseMeasureResult(std::ostream& os);

private:
  // This will be set to the value of whichever qualifer (RISE, FALL or CROSS) is being
  // used by the measure.
  int RFC_;
};

//-------------------------------------------------------------------------
// Class         : DerivativeEvaluationCont
// Purpose       : Continuous version of DERIV measure that can return multiple
//                 measure values if RISE, FALL or CROSS is not specified.
//                 This class is invoked by .MEASURE AC_CONT, .MEASURE DC_CONT,
//                 .MEASURE NOISE_CONT and .MEASURE TRAN_CONT.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 08/21/2020
//-------------------------------------------------------------------------
class DerivativeEvaluationCont : public DerivativeEvaluationBase
{
public:
  DerivativeEvaluationCont(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~DerivativeEvaluationCont() {};

  void reset();
  std::ostream& printMeasureResult(std::ostream& os);
  std::ostream& printVerboseMeasureResult(std::ostream& os);

  void updateCalculationResult(double val);
  void updateCalculationInstant(double val);

private:
  int contCross_;
  int contRise_;
  int contFall_;
  int contRFC_;
};

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif

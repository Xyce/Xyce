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

#include <N_IO_MeasureWhenAT.h>

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
// Creator       : Pete Sholander, SNL
// Creation Date : 08/21/2020
//-------------------------------------------------------------------------
class DerivativeEvaluationBase : public WhenAT
{
public:
  DerivativeEvaluationBase(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~DerivativeEvaluationBase() {};

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
  double getDerivativeValue(double currIndepVarValue) const;

private:
  void updateMeasureVarsForAT(double currIndepVarVal);
  void updateMeasureVarsForWhen(double currIndepVarVal, double whenInstant);
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

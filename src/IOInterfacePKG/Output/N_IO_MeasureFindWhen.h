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
// Purpose        : Find the time when a variable hits a target value (WHEN measure),
//                  or find the value of a variable at the time that the WHEN
//                  clause is satisfied (FIND-WHEN measure).
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
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-------------------------------------------------------------------------
class FindWhenBase : public Base
{
public:
  FindWhenBase(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~FindWhenBase() {};

  void prepareOutputVariables();
  bool checkMeasureLine();
  void resetFindWhenBase();

  virtual void updateCalculationResult(double val)=0;
  virtual void updateCalculationInstant(double val)=0;

  bool isATcondition(const double indepVarVal);
  bool isWHENcondition(const double indepVarVal, const double targVal);

  void setMeasureState(const double indepVarVal);
  void updateMeasureState(const double indepVarVal);
  double updateTargVal();

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

  std::ostream& printMeasureWindow(std::ostream& os, const double endSimTime,
				   const double startSweepVal, const double endSweepVal);
  std::ostream& printRFCWindow(std::ostream& os);

protected:
  bool doneIfFound_;
  std::vector<double> calculationResultVec_;
  std::vector<double> calculationInstantVec_;

private:
  double interpolateCalculationInstant(double currIndepVarValue, double targVal);
  double interpolateFindValue(double currIndepVarValue, double targVal, double whenTime);

  void updateRFCcountForWhen();
  bool withinRFCWindowForWhen();

  int numOutVars_;
  std::vector<double> outVarValues_;

  // These are used to interpolate the independent variable (time, frequency or
  // dcSweepVal) value when the simulation has reported values that bound the target value.
  // If the WHEN clause is of the form WHEN V(1)=V(2) then the previous V(1) value
  // is in the lastDepVarValue_ variable. The previous V(2) value is in the
  // lastTargValue_ variable.
  double lastIndepVarValue_;
  double lastDepVarValue_;
  double lastOutputVarValue_;
  double lastTargValue_;

  double startDCMeasureWindow_;
  int numPointsFound_;

  // This variable controls what is tested against in the WHEN clause.  
  // It refers to an index in the outputVarValues_ array.  Its value 
  // depends on whether the measure is a FIND-WHEN or a WHEN measure.
  // For a WHEN measure, whenIdx_ = 0 which is the default.  For a FIND-WHEN 
  // measure, whenIdx_ = 1  For the syntax WHEN v(a)=v(b), the value of v(a) 
  // is in outputVarValues_[whenIdx_] and the value of v(b) is in 
  // outputVarValues_[whenIdx_+1]
  int whenIdx_;

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

  void updateCalculationResult(double val)
  {
    calculationResult_ = val;
  }

  void updateCalculationInstant(double val)
  {
    calculationInstant_ = val;
  }

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
// Creation Date : 08/03/20209
//-------------------------------------------------------------------------
class FindWhenCont : public FindWhenBase
{
public:
  FindWhenCont(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~FindWhenCont() {};

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

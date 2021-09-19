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
// Purpose        : Implement TRIG-TARG measure for AC, DC, NOISE and TRAN
//                  measure modes.
// Special Notes  : This is the new class, that was implemented for Xyce 7.5
//
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 09/14/2021
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureTrigTarg_h
#define Xyce_N_IO_MeasureTrigTarg_h

#include <N_IO_MeasureBase.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-------------------------------------------------------------------------
// Class         : TrigTargBase
// Purpose       : This class contains the variables and functions
//                 that are common to both the continous and non-continuous
//                 versions of this measure.  The non-continuous version
//                 (e.g., .MEASURE TRAN) will only return one measure value.
//                 The continuous version (e.g., .MEASURE TRAN_CONT) may
//                 return multiple values.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-------------------------------------------------------------------------
class TrigTargBase : public Base
{
public:
  TrigTargBase(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~TrigTargBase() {};

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

  double getMeasureResult();

  std::ostream& printRFCWindow(std::ostream& os) const;
  virtual void printMeasureWarnings(double endSimTime, double startSweepVal,
                                      double endSweepVal) const;
  virtual void printMeasureWarningsForAT(double endSimTime) const;

protected:
  void resetTrigTargBase();
  bool checkMeasureLine() const;

  virtual void updateTrigResult(double val)=0;
  virtual void updateTargResult(double val)=0;

  std::vector<double> trigResultVec_;
  std::vector<double> targResultVec_;
  double trigResult_;
  double targResult_;
  bool trigResultFound_;
  bool targResultFound_;

  // This will be set to the value of whichever qualifer (RISE, FALL or CROSS) is being
  //  used by the measure.
  int trigRFC_;
  int targRFC_;

  // added to allow RISE, FALL and CROSS to work separately for
  // TRIG and TARG variables
  int actualTrigRise_;
  int actualTrigFall_;
  int actualTrigCross_;
  bool measureTrigLastRFC_;

  int actualTargRise_;
  int actualTargFall_;
  int actualTargCross_;
  bool measureTargLastRFC_;

private:
  void setMeasureState(double indepVarVal);
  void updateMeasureState(double indepVarVal);

  double updateTrigTargetVal();
  double updateTargTargetVal();

  bool isWHENcondition(
    double indepVarValue, 
    double depVarValue,
    double lastDepVarValue,
    double targetValue,
    double lastTargetValue) const;

  void updateTrigVarsForAT();
  void updateTargVarsForAT();
  void updateTrigVarsForWhen(double whenInstant);
  void updateTargVarsForWhen(double whenInstant);

  double interpolateCalculationInstant(
    double currIndepVarValue,
    double currDepVarValue,
    double lastDepVarValue,
    double targetVal,
    double lastTargetVal) const;

  bool isInvalidAT(double atVal, double startVal, double stopVal) const;

  void updateTrigRFCcount();
  void updateTargRFCcount();

  bool withinTrigRFCWindow() const;
  bool withinTargRFCWindow() const;

  bool withinTrigTDwindow(double indepVarVal) const;
  bool withinTargTDwindow(double indepVarVal) const;

  bool withinTrigTDwindowForDC(double sweepVal) const;
  bool withinTargTDwindowForDC(double sweepVal) const;

  int numOutVars_;
  std::vector<double> outVarValues_;

  double startDCMeasureWindow_;
  int numPointsFound_;

  // This variable controls what is tested against in the TARG clause.  
  // It refers to an index in the outputVarValues_ array.  Its value 
  // depends on whether the TRIG clause uses the AT=<time> or v(a)=<val>
  // syntax. The default values of targIdx_=1 works if the TRIG clause uses 
  // the v(a)=<val> syntax.  If the TRIG clause uses the AT=<time> syntax then 
  // the variable for the TARG clause is in outVarValues_[0].  If the TRIG clause 
  // uses the v(a)=v(b) syntax then the variable for the TARG clause is in 
  // outVarValues_[2].
  int targIdx_;

  // These are used to interpolate the independent variable (time, frequency or
  // dcSweepVal) value when the simulation has reported values that bound the target value.
  double lastIndepVarValue_;
  double lastTrigOutputValue_;
  double lastTargOutputValue_;
  double lastTrigValue_;
  double lastTargValue_;

  bool trigCalculationDone_;
  bool targCalculationDone_;
};

//-------------------------------------------------------------------------
// Class         : TrigTarg
// Purpose       : Non-continuous version of TrigTarg, that returns only
//                 only one measure value.  This class is invoked by
//                 .MEASURE AC, .MEASURE DC, .MEASURE NOISE and .MEASURE TRAN.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-------------------------------------------------------------------------
class TrigTarg : public TrigTargBase
{
public:
  TrigTarg(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~TrigTarg() {};

  void reset();

  void updateTrigResult(double val);
  void updateTargResult(double val);

  std::ostream& printMeasureResult(std::ostream& os);
  std::ostream& printVerboseMeasureResult(std::ostream& os);
};

//-------------------------------------------------------------------------
// Class         : TrigTargCont
// Purpose       : Continuous version of TrigTarg that can return multiple
//                 measure values if RISE, FALL or CROSS is not specified.
//                 This class is invoked by .MEASURE AC_CONT, .MEASURE DC_CONT,
//                 .MEASURE NOISE_CONT and .MEASURE TRAN_CONT.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 09/14/2021
//-------------------------------------------------------------------------
class TrigTargCont : public TrigTargBase
{
public:
  TrigTargCont(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~TrigTargCont() {};

  void reset();
  
  void updateTrigResult(double val);
  void updateTargResult(double val);

  std::ostream& printMeasureResult(std::ostream& os);
  std::ostream& printVerboseMeasureResult(std::ostream& os);
};

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif

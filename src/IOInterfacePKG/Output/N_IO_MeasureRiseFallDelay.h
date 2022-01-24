//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Purpose        : Measure rise/fall delay times for TRIG-TARG measures
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 03/10/2009
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureRiseFallDelay_h
#define Xyce_N_IO_MeasureRiseFallDelay_h

#include <N_IO_MeasureBase.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-------------------------------------------------------------------------
// Class         : RiseFallDelay
// Purpose       : Measure Rise/fall/delay times
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-------------------------------------------------------------------------
class RiseFallDelay : public Base
{
public:
  RiseFallDelay(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
  ~RiseFallDelay() {};

  void prepareOutputVariables();
  bool checkMeasureLine() const; 
  void reset();

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

  double getMeasureResult();
  std::ostream& printMeasureResult(std::ostream& os);
  std::ostream& printVerboseMeasureResult(std::ostream& os);
  bool withinTrigRiseFallCrossWindow();
  bool withinTargRiseFallCrossWindow(); 
  bool newTrigRiseFallCrossWindowforLast();
  bool newTargRiseFallCrossWindowforLast();
  void updateTrigRiseFallCrossCounts(double, double);
  void updateTargRiseFallCrossCounts(double, double);
  void updateTrigTargRiseFallCrossCounts(double, double, bool, bool, bool,
                                         bool&, bool&, int& , int&, int&, double&);

private:
  bool trigVariableLengthHistoryNeeded_;
  bool targVariableLengthHistoryNeeded_;
  double trigMax_;
  double targMax_;
  int trigResultIndex_;
  int targResultIndex_;
  double timeForTrig_;
  double timeForTarg_;
  bool trigMaxChanged_;
  bool targMaxChanged_;
  bool timeForTrigFound_;
  bool timeForTargFound_;
  bool trigOutputValueTargetChanged_;
  bool targOutputValueTargetChanged_;
  int numOutVars_;

  // added to make frac_max work with RISE/FALL/CROSS windows for TRIG/TARG
  double prevOutputVar0_, prevOutputVar1_;
  std::vector<double> outVarValues_;
  // these are vectors to store history information.
  // we need two independent var vectors because we will
  // trim the vectors dynamically to keep down on memory use
  std::vector<double> trigIndepVarHistory_; // usually time
  std::vector<double> trigVarHistory_;      // the trigger history
  std::vector<double> targIndepVarHistory_; // usually time
  std::vector<double> targetVarHistory_;    // the target history

  // This variable controls what is tested against in the TARG clause.  
  // It refers to an index in the outputVarValues_ array.  Its value 
  // depends on whether the TRIG clause uses the AT=<time> or v(a)=<val>
  // syntax. The default values of targIdx_=1 works if the TRIG clause uses 
  // the v(a)=<val> syntax.  If the TRIG clause uses the AT=<time> syntax then 
  // the variable for the TARG clause is in outVarValues_[0].  If the TRIG clause 
  // uses the v(a)=v(b) syntax then the variable for the TARG clause is in 
  // outVarValues_[2].
  int targIdx_;

  // added to allow RISE, FALL and CROSS to work separately for
  // TRIG and TARG variables
  int actualTrigRise_;
  int actualTrigFall_;
  int actualTrigCross_;
  bool isTrigRising_;
  bool isTrigFalling_;
  double lastTrigOutputValue_;  
  bool newTrigRiseWindow_;
  bool newTrigFallWindow_;
  bool newTrigCrossWindow_;

  int actualTargRise_;
  int actualTargFall_;
  int actualTargCross_;
  bool isTargRising_;
  bool isTargFalling_;
  double lastTargOutputValue_;
  bool newTargRiseWindow_;
  bool newTargFallWindow_;
  bool newTargCrossWindow_;
};

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif

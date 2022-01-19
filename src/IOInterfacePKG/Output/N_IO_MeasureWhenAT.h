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
// Purpose        : Base class functions that are common to DERIV-AT, DERIV-WHEN
//                  FIND-AT, FIND-WHEN and WHEN measures
// Special Notes  :
//
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 09/27/2021
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureWhenAT_h
#define Xyce_N_IO_MeasureWhenAT_h

#include <N_IO_MeasureBase.h>

namespace Xyce {
namespace IO {
namespace Measure {

class WhenAT : public Base
{
public:
  WhenAT(const Manager &measureMgr, const Util::OptionBlock & measureBlock, int whenIdx);
  ~WhenAT() {};

  std::ostream& printMeasureWindow(std::ostream& os, double endSimTime,
			          double startSweepVal, double endSweepVal) const;
  std::ostream& printRFCWindow(std::ostream& os) const;

protected:
  void resetWhenAT();

  void updateCalculationResult(double val);
  void updateCalculationInstant(double val);
  std::vector<double> calculationResultVec_;
  std::vector<double> calculationInstantVec_;

  void setMeasureState(double indepVarVal);
  void updateMeasureState(double indepVarVal);

  double getTargVal() const;
  void updateLastTargVal();

  bool isATcondition(double indepVarVal) const;
  bool isWHENcondition(double indepVarVal, double targVal) const;

  double interpolateCalculationInstant(double currIndepVarValue, double targVal) const;
 
  void updateRFCcountForWhen();
  bool withinRFCWindowForWhen() const;

  int numOutVars_;
  std::vector<double> outVarValues_;

  // This variable controls what is tested against in the WHEN clause.  
  // It refers to an index in the outputVarValues_ array.  Its value 
  // depends on whether the measure is a FIND-WHEN or a WHEN measure.
  // For a WHEN measure, whenIdx_ = 0 which is the default.  For a FIND-WHEN 
  // or DERIV-WHEN measure, whenIdx_ = 1  For the syntax WHEN v(a)=v(b),
  // the value of v(a) is in outputVarValues_[whenIdx_] and the value of
  // v(b) is in outputVarValues_[whenIdx_+1]
  int whenIdx_;

  /// This will be set to the value of whichever qualifer (RISE, FALL or CROSS) is being
  // used by the measure.
  int RFC_;

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
};

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif

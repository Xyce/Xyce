//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : Measure statistics of a simulation variable
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 03/10/2009
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureFourier_h
#define Xyce_N_IO_MeasureFourier_h

#include <N_IO_MeasureBase.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-------------------------------------------------------------------------
// Class         : Fourier
// Purpose       : Calculate Fourier transform of a simulation variable
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 06/05/2013
//-------------------------------------------------------------------------
class Fourier : public Base
{
  public:
  Fourier(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
    ~Fourier() {};

    void prepareOutputVariables();
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

    std::ostream& printMeasureResult(std::ostream& os);
    std::ostream& printVerboseMeasureResult(std::ostream& os);
    void printMeasureWarningsForAT(double endSimTime) const;
    std::ostream& printMeasureWindow(std::ostream& os, double endSimTime,
				     double startSweepVal, double endSweepVal) const;

  private:
    void calculateMeasureResult_();
    void getLastPeriod_();
    bool interpolateData_();
    void calculateFT_();
    std::string type_;
    int numOutVars_, prdStart_; 
    std::vector<double> outVarValues_, time_, newTime_, newValues_, mag_, phase_, nmag_, nphase_, freq_;
    double period_, lastPrdStart_, thd_;
    bool calculated_;
    bool lastPeriodFound_;
};

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif

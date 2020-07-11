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
// Purpose        : Implement ERROR measure, which is different from the
//                  ERR, ERR1 and ERR2 measures
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 03/10/2009
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureError_h
#define Xyce_N_IO_MeasureError_h

#include <N_IO_MeasureBase.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-------------------------------------------------------------------------
// Class         : Error
// Purpose       : Implement ERROR measure
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-------------------------------------------------------------------------
class Error : public Base
{
  public:
  Error(const Manager &measureMgr, const Util::OptionBlock & measureBlock );
    ~Error() {};
    
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

    double getMeasureResult();
    std::ostream& printMeasureResult(std::ostream& os);
    std::ostream& printVerboseMeasureResult(std::ostream& os);

    // used to print message about measurement time/frequency window, etc.
    std::ostream& printMeasureWindow(std::ostream& os, const double endSimTime,
				     const double startSweepVal, const double endSweepVal);

  private:
    // these are used to hold the variable names, and data from the external file
    std::vector< std::string > varNames_;
    std::vector< double > indepVarValues_;
    std::vector< double > dataValues_;    

    int numOutVars_;
    std::vector<double> outVarValues_;

    // results from the simulation to compare to a column in dataValues.
    // simulationIndepVarVals_ is the "independent variable", or "X-axis.", in the
    // simulation data.  It is time or frequency for TRAN or AC modes.  It is the 
    // first variable in the DC Sweep Vector, for DC mode.  simulationDataVals_ 
    // is then the "dependent variable", or Y-axis, in the simulation data.
    std::vector<double> simulationIndepVarVals_, simulationDataVals_;       
};

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif

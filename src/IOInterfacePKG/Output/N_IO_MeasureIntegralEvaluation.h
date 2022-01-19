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
// Purpose        : Measure an integral of an output variable
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 03/10/2009
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureIntegralEvaluation_h
#define Xyce_N_IO_MeasureIntegralEvaluation_h

#include <N_IO_MeasureBase.h>
#include <N_IO_MeasureStats.h>

namespace Xyce {
namespace IO {
namespace Measure {


//-------------------------------------------------------------------------
// Class         : IntegralEvaluation
// Purpose       : Measure an integral of an output variable
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-------------------------------------------------------------------------
class IntegralEvaluation : public Stats
{
  public:
  IntegralEvaluation(const Manager &measureMgr, const Util::OptionBlock & measureBlock);
    ~IntegralEvaluation()
    {}

    void reset();

    double getMeasureResult();

    void setMeasureVarsForNewWindow();
    void updateMeasureVars(double indepVarVal, double signalVal);

  private:
    std::string type_;
    double integralValue_;
};

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif

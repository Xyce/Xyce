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

//-------------------------------------------------------------------------
//
// Purpose        : Forward declarations
//
// Special Notes  : Forward declaring everything as a class breaks if the implementation of the type changes (like during
//                  templatization)
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_fwd_h
#define Xyce_N_NLS_fwd_h

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Enum          : AnalysisMode
// Description   : These anlysis modes influence the choice of parameters for
//                 the nonlinear solver. The mode is set by the
//                 setAnalysisMode() function in the Manager class.
//-----------------------------------------------------------------------------
enum AnalysisMode
{
  DC_OP,
  DC_SWEEP,
  DC_NLPOISSON,
  TRANSIENT,
  HB_MODE,
  NUM_MODES
};

class Manager;
class ParamMgr;
class ReturnCodes;
class NonLinInfo;
class NLParams;
class NonLinearSolver;
class ConductanceExtractor;
class Sensitivity;
class TwoLevelNewton;
class ConstraintBT;
class DampedNewton;

namespace N_NLS_NOX {
  class SharedSystem;
  class Group;
  class Vector;
  class AugmentLinSys;
  class Interface;
}

namespace N_NLS_LOCA {
  class Group;
}

} // namespace Nonlinear
} // namespace Xyce

namespace LOCA {
  class GlobalData;
  class Stepper;
  namespace MultiContinuation {
    class AbstractGroup;
  }
  namespace StatusTest {
    class Wrapper;
  }
}

namespace NOX {
  namespace Abstract {
    class Vector;
    class Group;
  }
  namespace Parameter {
    class List;
  }
  namespace Status {
    class Test;
  }
}

#endif // Xyce_N_NLS_fwd_h

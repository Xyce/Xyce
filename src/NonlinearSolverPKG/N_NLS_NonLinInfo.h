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
// Purpose       : 
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 2/11/07
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NON_LIN_INFO_H
#define Xyce_N_NLS_NON_LIN_INFO_H

// ---------- Standard Declarations ----------
#include <N_NLS_TwoLevelEnum.h>

// ---------- Forward Declarations ----------

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Class         : NonLinInfo
// Purpose       : 
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 2/11/07
//-----------------------------------------------------------------------------
class NonLinInfo
{
public:
  NonLinInfo():
  newtonIter(0),
  twoLevelNewtonCouplingMode (FULL_PROBLEM),
  locaFlag(false),
  continuationStep(0),
  firstContinuationParam(false),
  firstSolveComplete(false)
  {};

  virtual ~NonLinInfo() {};

  int newtonIter;
  TwoLevelNewtonMode twoLevelNewtonCouplingMode;

  // LOCA/homotopy related stuff.
  bool locaFlag;
  int continuationStep;
  bool firstContinuationParam;
  bool firstSolveComplete;


};

} // namespace Nonlinear
} // namespace Xyce

#endif // Xyce_N_NLS_NON_LIN_INFO_H



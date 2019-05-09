//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
//
// Creator       : Eric Keiter
//
// Creation Date : 1/29/07
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_TIA_TwoLevelError.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace TimeIntg {

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : "<<" operator for two level error class.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/17/05
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const TwoLevelError & tlerror)
{
  os.width(20);os.precision(12);os.setf(std::ios::scientific);
  os << "\n-----------------------------------------" << std::endl;
  os << "\tTwoLevelError:\n";
  os << "\t    innerSize:\t" << tlerror.innerSize << std::endl;
  os << "\t    xErrorSum:\t" << tlerror.xErrorSum << std::endl;
  os << "\t    qErrorSum:\t" << tlerror.qErrorSum << std::endl;
  os << "\t xErrorSum_m1:\t" << tlerror.xErrorSum_m1 << std::endl;
  os << "\t xErrorSum_p1:\t" << tlerror.xErrorSum_p1 << std::endl;
  os << "\t q1HistorySum:\t" << tlerror.q1HistorySum << std::endl;
  os << Xyce::section_divider << std::endl;
  os << std::endl;

  return os;
}

} // namespace TimeIntg
} // namespace Xyce

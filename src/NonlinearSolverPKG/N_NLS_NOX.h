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

//-------------------------------------------------------------------------
//
// Purpose        : Error-handling for N_NLS_NOX
//
// Special Notes  :
//
// Creator        : Tammy Kolda, NLS, 8950
//
// Creation Date  : 01/31/02
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NOX_h
#define Xyce_N_NLS_NOX_h

// ---------- Standard Includes ----------

#include <string>

// ----------   Xyce Includes   ----------

// ----------   NOX Includes   ----------

// ---------- Forward Declarations ----------

// ---------- Namespace Declarations ----------

// ---------- Function Declarations ----------

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

  //---------------------------------------------------------------------------
  // Function      : error
  // Purpose       : This is a wrapper for
  //                 Xyce::Report::UserError0() <<  msg.
  //                 It throws "N_NLS_NOX Error" after calling report.
  //---------------------------------------------------------------------------
  void error(const std::string &msg);

  //---------------------------------------------------------------------------
  // Function      : warning
  // Purpose       : This is simply a wrapper for
  //                 Xyce::Report::UserWarning0() <<  msg.
  //---------------------------------------------------------------------------
  void warning(const std::string &msg);

  //---------------------------------------------------------------------------
  // Function      : info
  // Purpose       : This is simply a wrapper for
  //                 lout() << msg)
  //---------------------------------------------------------------------------
  void info(const std::string &msg);

}} // namespace N_NLS_NOX
} // namespace Xyce

#endif


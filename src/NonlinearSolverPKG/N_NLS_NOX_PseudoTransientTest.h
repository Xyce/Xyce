//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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

#ifndef Xyce_N_NLS_NOX_PseudoTransientTest_h
#define Xyce_N_NLS_NOX_PseudoTransientTest_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

// ----------   NOX Includes   ----------

#include <N_NLS_fwd.h>

#include "NOX_StatusTest_Generic.H"	// base class

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

//-----------------------------------------------------------------------------
// Class         : PseudoTransientTest
//
// Purpose       : a NOX-compatiable status test 
//                 for pseudo transient continuation.
//
//      
//
// Creator       : Roger Pawlowski
//
// Creation Date : 2/1/02
//-----------------------------------------------------------------------------

class PseudoTransientTest : public NOX::StatusTest::Generic {

public:

  //---------------------------------------------------------------------------
  // Function      : PseudoTransientTest (constructor)
  //
  // Purpose       : Constructs a NOX-compatiable status test 
  //                 for pseudo transient continuation.
  //
  //---------------------------------------------------------------------------
    PseudoTransientTest(double maxStepSize, double minNormF);

  //---------------------------------------------------------------------------
  // Function      : Destructor
  //---------------------------------------------------------------------------
  ~PseudoTransientTest();

  //---------------------------------------------------------------------------
  // Purpose       : Test stopping criterion given the current
  //                 nonlinear problem
  //---------------------------------------------------------------------------
  NOX::StatusTest::StatusType 
    checkStatus(const NOX::Solver::Generic& problem,
		NOX::StatusTest::CheckType checkType);
  //---------------------------------------------------------------------------
  // Purpose       : Test stopping criterion given the current
  //                 nonlinear problem
  //---------------------------------------------------------------------------
  NOX::StatusTest::StatusType getStatus() const { return status_; };

  //---------------------------------------------------------------------------
  // Purpose       : Output formatted description of stopping test to
  //                 output stream.
  //---------------------------------------------------------------------------
  std::ostream& print(std::ostream& stream, int indent = 0) const;

private:

  NOX::StatusTest::StatusType status_;
  double maxStepSize_;
  double currentStepSize_;
  double minNormF_;
  double currentNormF_;
  
}; 
}}} // namespace N_NLS_NOX

#endif


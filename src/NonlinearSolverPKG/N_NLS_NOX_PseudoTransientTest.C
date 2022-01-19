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

//-------------------------------------------------------------------------
//
// Purpose        : Status test.
//
// Special Notes  :
//
// Creator        : Roger Pawlowski, SNL 9233
//
// Creation Date  : 04/15/03
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_fwd.h>

#include "N_NLS_NOX_PseudoTransientTest.h"
#include "NOX.H"
#include "N_NLS_NOX_PseudoTransientSolver.h"

// ----------   Namespaces   ----------

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

// ----------   Code   ----------

PseudoTransientTest::PseudoTransientTest(double maxStepSize,
					 double minNormF) :
  status_(NOX::StatusTest::Unconverged),
  maxStepSize_(maxStepSize),
  currentStepSize_(0.0),
  minNormF_(minNormF),
  currentNormF_(1.0)
{
}

PseudoTransientTest::~PseudoTransientTest()
{
}

NOX::StatusTest::StatusType PseudoTransientTest::
checkStatus(const NOX::Solver::Generic& problem,
	    NOX::StatusTest::CheckType checkType)
{
  status_ = NOX::StatusTest::Unconverged;

  const PseudoTransientBased* psts = 0;
  psts = dynamic_cast<const PseudoTransientBased*>(&problem);

  if (psts == 0) {
    dout() << "NOX::StatusTest::PseudoTransientTest::checkStatus - failed dynamic_cast solver to PseudoTransientBased!" << std::endl;
    throw "NOX Error";
  }
  currentStepSize_ = psts->getPseudoTransientStepSize();

  currentNormF_ = problem.getSolutionGroup().getNormF();

  if ((currentStepSize_ >= maxStepSize_) && (currentNormF_ < minNormF_))
    status_ = NOX::StatusTest::Converged;

  return status_;
}

std::ostream& PseudoTransientTest::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << status_;
  stream << "Pseudo Transient Step Size = " 
	 << NOX::Utils::sciformat(currentStepSize_, 5) 
	 << " >= " << NOX::Utils::sciformat(maxStepSize_, 5);
  stream << std::endl;

  for (int j = 0; j < indent; ++j)
    stream << ' ';
  stream << status_;
  stream << "Pseudo Transient Residual Reduction = " 
	 << NOX::Utils::sciformat(currentNormF_, 5) 
	 << " < " << NOX::Utils::sciformat(minNormF_, 5);
  stream << std::endl;

  return stream;
}

}}}

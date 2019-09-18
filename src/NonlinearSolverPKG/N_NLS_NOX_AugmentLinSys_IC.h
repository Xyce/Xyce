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

//-------------------------------------------------------------------------
//
// Purpose        : Concrete class for augmenting the Jacobian for
//                  .IC simulations (initial condition)
//
// Special Notes  :
//
// Creator        : Eric R. Keiter
//
// Creation Date  : 09/15/07
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NOX_AugmentLinSys_IC_h
#define Xyce_N_NLS_NOX_AugmentLinSys_IC_h

#include "N_NLS_NOX_AugmentLinSys.h"
#include "N_PDS_ParMap.h"

#include <N_UTL_fwd.h>
#include <N_IO_InitialConditions.h>

#include <vector>

//-----------------------------------------------------------------------------
// Class         : N_NLS::NOX::AugmentLinSysIC
// Purpose       : Handles matrix augmentation to support .IC statements.
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems modeling.
// Creation Date : 09/15/2007
//-----------------------------------------------------------------------------
namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

class AugmentLinSysIC : public AugmentLinSys {

public:
  //! Ctor.
  AugmentLinSysIC(Xyce::IO::InitialConditionsData::NodeNamePairMap & op_in,
                  const std::vector<int>& colors,
                  Xyce::Linear::Vector* cloneVector);

  //! Dtor.
  virtual ~AugmentLinSysIC();

  void setProgressVariable(double dummy) {return;}

  void augmentResidual(const Xyce::Linear::Vector * solution,
		       Xyce::Linear::Vector * residual_vector);

  void augmentJacobian(Xyce::Linear::Matrix * jacobian);

 private:

  //! map of specified variables
  Xyce::IO::InitialConditionsData::NodeNamePairMap & op_;

  //! Color 0 are the voltage unknowns.
  const std::vector<int> & colors_;

  //! Temporary vector used to store diagonal.
  Xyce::Linear::Vector* tmp_vector_ptr_;

};

}}}

#endif


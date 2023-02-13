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

//-------------------------------------------------------------------------
//
// Purpose        : Concrete class for augmenting the Jacobian for
//                  pseudo-transient solves.
//
// Special Notes  :
//
// Creator        : Roger Pawlowski, NLS, 9233
//
// Creation Date  : 3/6/2006
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NOX_AugmentLinSys_PseudoTransient_h
#define Xyce_N_NLS_NOX_AugmentLinSys_PseudoTransient_h

//-----------------------------------------------------------------------------
// Class         : N_NLS::NOX::
//
// Purpose       :
//
//      NOX Group Interface for Xyce
//
// Creator       : Roger Pawlowski, SNL, 9233
//
// Creation Date : 3/6/2006
//-----------------------------------------------------------------------------

#include "N_NLS_NOX_AugmentLinSys.h"          
#include <vector>

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

  class AugmentLinSysPseudoTransient : public AugmentLinSys {

public:

  //! Ctor.
  AugmentLinSysPseudoTransient(const std::vector<int>& colors,
			       Xyce::Linear::Vector* cloneVector,
			       bool useVoltageScaleFactor=false,
			       double voltageScaleFactor=1.0);

  //! Dtor.
  virtual ~AugmentLinSysPseudoTransient();

  void setProgressVariable(double time_step_size);

  void augmentResidual(const Xyce::Linear::Vector * solution,
		       Xyce::Linear::Vector * residual_vector);
  
  void augmentJacobian(Xyce::Linear::Matrix * jacobian);

 private:

  bool use_voltage_scale_factor_;

  double voltage_scale_factor_;

  double time_step_size_;

  const std::vector<int>& colors_;

  Xyce::Linear::Vector* tmp_vector_ptr_;

};
 
}}} 

#endif 


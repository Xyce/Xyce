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
// Purpose        : Strategy pattern for allowing algorithms to augment
//                  the Jacobian and residual.
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

#ifndef Xyce_N_NLS_NOX_AugmentLinSys_h
#define Xyce_N_NLS_NOX_AugmentLinSys_h

#include <N_LAS_fwd.h>

//-----------------------------------------------------------------------------
// Class         : N_NLS_NOX::AugmentLinSys
// Purpose       :
// Creator       : Roger Pawlowski, SNL, 9233
// Creation Date : 3/6/2006
//-----------------------------------------------------------------------------

/*! \brief Pure virtual class to augment a linear system.

    Strategy pattern for allowing algorithms to augment the Jacobian
    and residual.  This class is used to provide an algorithm for
    augmenting the linear system for various solution techniques.
    Homotopy, pseudo-transient, and "gmin" stepping are examples.  In
    each case, the Jacobian's diagonal is changed/added to.  Some
    algorithms also augment the residual.
  
*/

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

class AugmentLinSys {

public:

  //! Ctor.
  AugmentLinSys(){};

  //! Dtor.
  virtual ~AugmentLinSys() {};

  //! Set the progress variable (time step size for pseudo transient).
  virtual void setProgressVariable(double value) = 0;

  //! Augments the Residual.
  virtual void augmentResidual(const Xyce::Linear::Vector * solution,
			       Xyce::Linear::Vector * residual_vector) = 0;
  
  //! Augments the Jacobian.
  virtual void augmentJacobian(Xyce::Linear::Matrix * jacobian) = 0;
  
};
 
}}} 

#endif 


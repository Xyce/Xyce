//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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

#ifndef Xyce_N_NLS_NOX_AugmentLinSys_GStepping_h
#define Xyce_N_NLS_NOX_AugmentLinSys_GStepping_h


#include <vector>
#include "N_NLS_NOX_AugmentLinSys.h"          

//-----------------------------------------------------------------------------
// Class         : N_NLS::NOX::GStepping
// Purpose       : Handles matrix augmentation to support GMIN stepping.
// Creator       : Roger Pawlowski, SNL, 9233
// Creation Date : 3/6/2006
//-----------------------------------------------------------------------------
namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {
  
  class GStepping : public AugmentLinSys {
    
  public:

    enum NodeListType {
      NLT_VoltageNodes,
      NLT_AllVoltageUnknowns
    };

  public:
    
    //! Ctor for the voltage nodes as a GID list.
    GStepping(NodeListType node_list_type,
              const std::vector<int>& vnodeVec,
	      Xyce::Linear::Vector* cloneVector,
	      double endValue,
              double residCond=0);
    
    //! Dtor.
    virtual ~GStepping();
    
    void setProgressVariable(double time_step_size);
    
    void augmentResidual(const Xyce::Linear::Vector * solution,
			 Xyce::Linear::Vector * residual_vector);
    
    void augmentJacobian(Xyce::Linear::Matrix * jacobian);

    inline void setResidualConductance(double c) {residualConductance_=c;};

  private:

    //! Type of list we are using.
    NodeListType node_list_type_;

    //! Conductance.
    double conductance_;
    
    //! List of voltage unknowns or node GIDs.
    const std::vector<int>& vnodeVec_;
    
    //! Temporary vector used to store diagonal.
    Xyce::Linear::Vector* tmp_vector_ptr_;
    
    //! low end of the exponential term.
    double scaled_end_value_;

    //! residual value of the conductance.  Should almost always be zero
    double residualConductance_;

  };
  
}}} 

#endif 


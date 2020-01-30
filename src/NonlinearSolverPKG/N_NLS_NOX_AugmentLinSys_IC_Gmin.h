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
// Purpose        : Concrete class for augmenting the Jacobian for
//                  .IC simulations (initial condition) with gmin stepping.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter
//
// Creation Date  : 04/29/12
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NOX_AugmentLinSys_IC_Gmin_h
#define Xyce_N_NLS_NOX_AugmentLinSys_IC_Gmin_h

#include <vector>

#include "N_NLS_NOX_AugmentLinSys.h"
#include <N_IO_InitialConditions.h>

#include <N_UTL_fwd.h>

//-----------------------------------------------------------------------------
// Class         : N_NLS::NOX::AugmentLinSysIC_Gmin
// Purpose       : Handles matrix augmentation to support .IC statements, if
//                 gmin stepping is also being applied.
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems modeling.
// Creation Date : 04/29/2012
//-----------------------------------------------------------------------------
namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

class AugmentLinSysIC_Gmin : public AugmentLinSys {

  public:
    enum NodeListType {
      NLT_VoltageNodes,
      NLT_AllVoltageUnknowns
    };

  public:
    //! Ctor.
    AugmentLinSysIC_Gmin
    ( NodeListType node_list_type,
      Xyce::IO::InitialConditionsData::NodeNamePairMap & op_in,
      const std::vector<int>& ic_colors,
      const std::vector<int>& vnodeVec,
      Xyce::Linear::Vector* cloneVector,
      double scaledEndValue,
      double resCond);

    //! Dtor.
    virtual ~AugmentLinSysIC_Gmin();

    void setProgressVariable(double dummy);

    void augmentResidual(const Xyce::Linear::Vector * solution,
             Xyce::Linear::Vector * residual_vector);

    void augmentJacobian(Xyce::Linear::Matrix * jacobian);

  private:

    //! Type of list we are using.
    NodeListType node_list_type_;

    //! Conductance.
    double conductance_;

    //! low end of the exponential term.
    double scaled_end_value_;

    //! residual value of the conductance.  Should almost always be zero
    double residualConductance_;

    //! map of specified variables
    Xyce::IO::InitialConditionsData::NodeNamePairMap & op_;

    //! Color 0 are the voltage unknowns.
    //! For the IC color map, the voltage nodes attached to
    //! independent voltage sources are not included.
    const std::vector<int>& ic_colors_;
    const std::vector<int>& vnodeVec_;

    //! Temporary vectors used to store diagonal.
    Xyce::Linear::Vector* vecptr1_;
    Xyce::Linear::Vector* vecptr2_;

};

}}}

#endif


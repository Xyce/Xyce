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
// Purpose        : Algorithm for augmenting the Jacobian for pseudo
//                  transient solves using vnode conductance.
//
// Special Notes  :
//
// Creator        : Roger Pawlowski, SNL 9233
//
// Creation Date  : 03/07/06
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include "N_LAS_Vector.h"
#include "N_LAS_Matrix.h"
#include "Epetra_MapColoring.h"
#include "N_NLS_NOX_AugmentLinSys_GStepping.h"

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

//-----------------------------------------------------------------------------
// Function      : GStepping::GStepping
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date : 
//-----------------------------------------------------------------------------
GStepping::GStepping(
        const std::vector<int>& vnodeGIDVec,
				Linear::Vector* cloneVector,
				double scaledEndValue,
        double resCond) :
  node_list_type_(NLT_VoltageNodes),
  vnodeGIDVec_(vnodeGIDVec),
  tmp_vector_ptr_(0),
  scaled_end_value_(scaledEndValue),
  residualConductance_(resCond)
{
  tmp_vector_ptr_ = new Linear::Vector(*cloneVector);
}

//-----------------------------------------------------------------------------
// Function      : GStepping::GStepping
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date : 
//-----------------------------------------------------------------------------
GStepping::
GStepping(const Teuchos::RCP<Epetra_MapColoring>& color_map,
	  Linear::Vector* cloneVector,
	  double scaledEndValue,
    double resCond) :
  node_list_type_(NLT_AllVoltageUnknowns),
  scaled_end_value_(scaledEndValue),
  residualConductance_(resCond)
{
  color_map_ = color_map;
  tmp_vector_ptr_ = new Linear::Vector(*cloneVector);
}

//-----------------------------------------------------------------------------
// Function      : GStepping::~GStepping
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date : 
//-----------------------------------------------------------------------------
GStepping::~GStepping()
{
  delete tmp_vector_ptr_;
}


//-----------------------------------------------------------------------------
// Function      : GStepping::setProgressVariable
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date : 
//-----------------------------------------------------------------------------
void GStepping::setProgressVariable(double conductance)
{
  // Exponential Continuation (con param goes from +4 -> -log10(endValue))
  conductance_ = pow(10.0, conductance) - pow(10.0, scaled_end_value_) + residualConductance_;
}

//-----------------------------------------------------------------------------
// Function      : GStepping::augmentResidual
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date : 
//-----------------------------------------------------------------------------
void GStepping::augmentResidual(const Linear::Vector * solution,
					   Linear::Vector * residualVector)
{
  if (node_list_type_ == NLT_VoltageNodes) 
  {
    std::vector<int>::const_iterator i = vnodeGIDVec_.begin();
    std::vector<int>::const_iterator stop = vnodeGIDVec_.end();
    for ( ; i < stop; ++i) 
    {
      double value = conductance_ * 
        (const_cast<Linear::Vector*>(solution))->getElementByGlobalIndex(*i);

      residualVector->sumElementByGlobalIndex(*i, value);
    }
  }
  else 
  {
    for (std::size_t i = 0; i <  tmp_vector_ptr_->localLength(); ++i) 
    {
      if ( (*color_map_)[i] == 0)
      {
        (*residualVector)[i] += conductance_ * (const_cast<Linear::Vector&>(*solution))[i]; 
      }
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : GStepping::augmentJacobian
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date : 
//-----------------------------------------------------------------------------
void GStepping::augmentJacobian(Linear::Matrix * jacobian)
{
  jacobian->getDiagonal(*tmp_vector_ptr_);
  
  if (node_list_type_ == NLT_VoltageNodes) 
  {
    std::vector<int>::const_iterator i = vnodeGIDVec_.begin();
    std::vector<int>::const_iterator stop = vnodeGIDVec_.end();
    for ( ; i < stop; ++i) 
    {
      tmp_vector_ptr_->sumElementByGlobalIndex(*i, conductance_);
    }
  }
  else 
  {
    for (std::size_t i = 0; i <  tmp_vector_ptr_->localLength(); ++i) 
    {
      if ( (*color_map_)[i] == 0)
      {
        (*tmp_vector_ptr_)[i] += conductance_; 
      }
    }
  }

  jacobian->replaceDiagonal(*tmp_vector_ptr_);
}

}}}

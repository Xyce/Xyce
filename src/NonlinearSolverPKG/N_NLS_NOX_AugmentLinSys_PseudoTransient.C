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
// Purpose        : Algorithm for augmenting the Jacobian for pseudo
//                  transient solves.
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
#include "N_NLS_NOX_AugmentLinSys_PseudoTransient.h"

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

//-----------------------------------------------------------------------------
// Function      : AugmentLinSysPseudoTransient::AugmentLinSysPseudoTransient
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date :
//-----------------------------------------------------------------------------
AugmentLinSysPseudoTransient::AugmentLinSysPseudoTransient
                            (const std::vector<int>& colors,
			     Linear::Vector* cloneVector,
			     bool useVoltageScaleFactor,
			     double voltageScaleFactor)
  : use_voltage_scale_factor_(useVoltageScaleFactor),
    voltage_scale_factor_(voltageScaleFactor),
    colors_(colors),
    tmp_vector_ptr_(0)
{
  tmp_vector_ptr_ = new Linear::Vector(*cloneVector);
}

//-----------------------------------------------------------------------------
// Function      : AugmentLinSysPseudoTransient::~AugmentLinSysPseudoTransient
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date :
//-----------------------------------------------------------------------------
AugmentLinSysPseudoTransient::~AugmentLinSysPseudoTransient()
{
  delete tmp_vector_ptr_;
}


//-----------------------------------------------------------------------------
// Function      : AugmentLinSysPseudoTransient::setProgressVariable
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date :
//-----------------------------------------------------------------------------
void AugmentLinSysPseudoTransient::setProgressVariable
  (double time_step_size)
{
  time_step_size_ = time_step_size;
}

//-----------------------------------------------------------------------------
// Function      : AugmentLinSysPseudoTransient::augmentResidual
// Purpose       :
// Special Notes : no-op for pseudo-transient.
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date :
//-----------------------------------------------------------------------------
void AugmentLinSysPseudoTransient::augmentResidual
  (const Linear::Vector * solution, Linear::Vector * residual_vector)
{
  // Nothing to do for pseudo transient!!
}

//-----------------------------------------------------------------------------
// Function      : AugmentLinSysPseudoTransient::augmentJacobian
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date :
//-----------------------------------------------------------------------------
void AugmentLinSysPseudoTransient::augmentJacobian
  (Linear::Matrix * jacobian)
{
  //cout << "Augmenting Jacobian for Pseudo Transient" << endl;
  //cout << "Pseudo Trans Step Size = " << pseudoTransientTimeStep_ << endl;
  
  //jacobian->print();

  //jacobian->scale(conParamValue);

  jacobian->getDiagonal(*tmp_vector_ptr_);
  
  //tmp_vector_ptr_->print();
    
  double value = 1.0 / time_step_size_;
    
  //cout << "Pseudo Transient Time Step Size = " << time_step_size_ << endl;

  if (!use_voltage_scale_factor_) 
  {
    tmp_vector_ptr_->addScalar(value);
  }
  else 
  {
    for (std::size_t i = 0; i <  tmp_vector_ptr_->localLength(); ++i) 
    {
      if ( colors_[i] == 0 )
      {
        (*tmp_vector_ptr_)[i] += value * voltage_scale_factor_; 
      }
      else
      {
        (*tmp_vector_ptr_)[i] += value;
      }
    }
    //RPP Might need to export local values for tmp_vector_ptr_ here
    //for parallel.
  }

  jacobian->replaceDiagonal(*tmp_vector_ptr_);

  //jacobian->print();
   
  //cout << "Time step size = " << time_step_size_ << endl;

  //if (use_voltage_scale_factor_)
  //cout << "Voltage scale factor = " << voltage_scale_factor_ << endl;
}

}}}

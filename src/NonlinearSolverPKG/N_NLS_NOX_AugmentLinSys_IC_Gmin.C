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
// Purpose        : Algorithm for augmenting the Jacobian for .IC_Gmin
//                  operating point solves.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 4/29/2012
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
#include "N_ERH_ErrorMgr.h"
#include "N_NLS_NOX_AugmentLinSys_IC_Gmin.h"

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

//-----------------------------------------------------------------------------
// Function      : AugmentLinSysIC_Gmin::AugmentLinSysIC_Gmin
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/29/2012
//-----------------------------------------------------------------------------
AugmentLinSysIC_Gmin::AugmentLinSysIC_Gmin(
      NodeListType node_list_type,
      Xyce::IO::InitialConditionsData::NodeNamePairMap & op_in,
      const std::vector<int>& ic_colors,
      const std::vector<int>& vnodeVec,
      Xyce::Linear::Vector* cloneVector,
      double scaledEndValue,
      double resCond)
    :
    node_list_type_(node_list_type),
    scaled_end_value_(scaledEndValue),
    residualConductance_(resCond),
    op_       (op_in),
    ic_colors_(ic_colors),
    vnodeVec_(vnodeVec),
    vecptr1_(0),
    vecptr2_(0)
{
  vecptr1_ = cloneVector->cloneVector();
  vecptr2_ = cloneVector->cloneVector();
}

//-----------------------------------------------------------------------------
// Function      : AugmentLinSysIC_Gmin::~AugmentLinSysIC_Gmin
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/29/2012
//-----------------------------------------------------------------------------
AugmentLinSysIC_Gmin::~AugmentLinSysIC_Gmin()
{
  delete vecptr1_;
  delete vecptr2_;
}

//-----------------------------------------------------------------------------
// Function      : IC_Gmin::setProgressVariable
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date :
//-----------------------------------------------------------------------------
void AugmentLinSysIC_Gmin::setProgressVariable(double conductance)
{
  // Direct continuation of conductance (con param goes from 1.0e4 -> 0.0
  //conductance_ = conductance;

  // Exponential Continuation (con param goes from +4 -> -log10(endValue))
  conductance_ = pow(10.0, conductance) - pow(10.0, scaled_end_value_) + residualConductance_;
}

//-----------------------------------------------------------------------------
// Function      : AugmentLinSysIC_Gmin::augmentResidual
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/29/2012
//-----------------------------------------------------------------------------
void AugmentLinSysIC_Gmin::augmentResidual
   (const Linear::Vector * solution, Linear::Vector * residual_vector)
{
#ifdef Xyce_DEBUG_IC_Gmin
  dout() << "Inside AugmentLinSysIC_Gmin::augmentResidual:"  << std::endl;
#endif

  // GMIN portion
   if (node_list_type_ == NLT_VoltageNodes)
  {
    std::vector<int>::const_iterator i = vnodeVec_.begin();
    std::vector<int>::const_iterator stop = vnodeVec_.end();
    for ( ; i < stop; ++i)
    {
      double value = conductance_ *
        (const_cast<Linear::Vector*>(solution))->getElementByGlobalIndex(*i);

      residual_vector->sumElementByGlobalIndex(*i, value);
    }
  }
  else
  {
    for (std::size_t i = 0; i <  vecptr1_->localLength(); ++i)
    {
      if ( vnodeVec_[i] == 0 )
      {
        (*residual_vector)[i] += conductance_ * (const_cast<Linear::Vector&>(*solution))[i];
      }
    }
  }

  // IC portion
  IO::InitialConditionsData::NodeNamePairMap::iterator op_i = op_.begin();
  IO::InitialConditionsData::NodeNamePairMap::iterator op_end = op_.end();
  for ( ; op_i != op_end ; ++op_i)
  {
    int row = (*op_i).second.first;
    int global_row(row);

    if ( ic_colors_[row] == 0 )
    {
      (*residual_vector)[row]  = 0.0;
    }
  }



  return;
}

//-----------------------------------------------------------------------------
// Function      : AugmentLinSysIC_Gmin::augmentJacobian
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/29/2012
//-----------------------------------------------------------------------------
void AugmentLinSysIC_Gmin::augmentJacobian(Linear::Matrix * jacobian)
{
#ifdef Xyce_DEBUG_IC_Gmin
  dout() << "Inside AugmentLinSysIC_Gmin::augmentJacobian:"  << std::endl;
#endif

  // GMIN portion
  jacobian->getDiagonal(*vecptr1_);
  if (node_list_type_ == NLT_VoltageNodes)
  {
    std::vector<int>::const_iterator i = vnodeVec_.begin();
    std::vector<int>::const_iterator stop = vnodeVec_.end();
    for ( ; i < stop; ++i)
    {
      vecptr1_->sumElementByGlobalIndex(*i, conductance_);
    }
  }
  else
  {
    for (std::size_t i = 0; i <  vecptr1_->localLength(); ++i)
    {
      if ( vnodeVec_[i] == 0 )
      {
        (*vecptr1_)[i] += conductance_;
      }
    }
  }
  jacobian->replaceDiagonal(*vecptr1_);

  // IC portion
  std::vector<int> col;
  std::vector<double> val;
  IO::InitialConditionsData::NodeNamePairMap::iterator op_i = op_.begin();
  IO::InitialConditionsData::NodeNamePairMap::iterator op_end = op_.end();

  jacobian->getDiagonal(*vecptr2_);

  for ( ; op_i != op_end ; ++op_i)
  {
    int row = (*op_i).second.first;
    int rowLen(0);
    int numEntries(0);

    if ( ic_colors_[row] == 0)
    {
      rowLen = jacobian->getLocalRowLength(row);

      col.resize(rowLen,0);
      val.resize(rowLen,0.0);
      jacobian->getLocalRowCopy(row, rowLen, numEntries, &val[0], &col[0]);

      // zero out the entire row.
      for (int i=0;i<val.size();++i) val[i] = 0.0;

      jacobian->putLocalRow(row, rowLen, &val[0], &col[0]);

      // set the diagonal to 1.0.
      (*vecptr2_)[row] = 1.0;
    }
  }

  jacobian->replaceDiagonal(*vecptr2_);

  return;
}

}}}

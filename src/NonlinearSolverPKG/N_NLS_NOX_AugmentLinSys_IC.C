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
// Purpose        : Algorithm for augmenting the Jacobian for .IC
//                  operating point solves.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 09/15/07
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_fwd.h>

// ----------   Xyce Includes   ----------

#include "N_LAS_Vector.h"
#include "N_LAS_Matrix.h"
#include "N_ERH_ErrorMgr.h"
#include "N_NLS_NOX_AugmentLinSys_IC.h"

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

//-----------------------------------------------------------------------------
// Function      : AugmentLinSysIC::AugmentLinSysIC
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/15/07
//-----------------------------------------------------------------------------
AugmentLinSysIC::AugmentLinSysIC
  (IO::InitialConditionsData::NodeNamePairMap & op_in,
  const std::vector<int>& colors,
  Linear::Vector* cloneVector
  )
  : op_       (op_in),
    colors_   (colors),
    tmp_vector_ptr_(0)
{
  tmp_vector_ptr_ = cloneVector->clone();
}

//-----------------------------------------------------------------------------
// Function      : AugmentLinSysIC::~AugmentLinSysIC
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/15/07
//-----------------------------------------------------------------------------
AugmentLinSysIC::~AugmentLinSysIC()
{
  delete tmp_vector_ptr_;
}

//-----------------------------------------------------------------------------
// Function      : AugmentLinSysIC::augmentResidual
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/15/07
//-----------------------------------------------------------------------------
void AugmentLinSysIC::augmentResidual
   (const Linear::Vector * solution, Linear::Vector * residual_vector)
{
#ifdef Xyce_DEBUG_IC
  dout() << "Inside AugmentLinSysIC::augmentResidual:"  << std::endl;
#endif

  IO::InitialConditionsData::NodeNamePairMap::iterator op_i = op_.begin();
  IO::InitialConditionsData::NodeNamePairMap::iterator op_end = op_.end();
  for ( ; op_i != op_end ; ++op_i)
  {
    int row = (*op_i).second.first;

    if ( colors_[row] == 0)
    {
      (*residual_vector)[row]  = 0.0;
    }
  }
  return;
}

//-----------------------------------------------------------------------------
// Function      : AugmentLinSysIC::augmentJacobian
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/15/07
//-----------------------------------------------------------------------------
void AugmentLinSysIC::augmentJacobian(Linear::Matrix * jacobian)
{
#ifdef Xyce_DEBUG_IC
  dout() << "Inside AugmentLinSysIC::augmentJacobian:"  << std::endl;
#endif

  std::vector<int> col;
  std::vector<double> val;
  IO::InitialConditionsData::NodeNamePairMap::iterator op_i = op_.begin();
  IO::InitialConditionsData::NodeNamePairMap::iterator op_end = op_.end();

  jacobian->getDiagonal(*tmp_vector_ptr_);

  for ( ; op_i != op_end ; ++op_i)
  {
    int row = (*op_i).second.first;
    int rowLen(0);
    int numEntries(0);

    if ( colors_[row] == 0)
    {
      rowLen = jacobian->getLocalRowLength(row);

      col.resize(rowLen,0);
      val.resize(rowLen,0.0);
      jacobian->getLocalRowCopy(row, rowLen, numEntries, &val[0], &col[0]);

      // zero out the entire row.
      for (int i=0;i<val.size();++i) val[i] = 0.0;

      jacobian->putLocalRow(row, rowLen, &val[0], &col[0]);

      // set the diagonal to 1.0.
      (*tmp_vector_ptr_)[row] = 1.0;
    }
  }

  jacobian->replaceDiagonal(*tmp_vector_ptr_);

  return;
}

}}}

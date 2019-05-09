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
//                  transient solves.
//
// Special Notes  :
//
// Creator        : Dave Shirley, PSSI
//
// Creation Date  : 05/08/06
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
#include "N_NLS_NOX_AugmentLinSys_OPStart.h"

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

//-----------------------------------------------------------------------------
// Function      : AugmentLinSysOPStart::AugmentLinSysOPStart
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date :
//-----------------------------------------------------------------------------
AugmentLinSysOPStart::AugmentLinSysOPStart(
  IO::InitialConditionsData::NodeNamePairMap & op_in,
  const NodeNameMap & allNodes_in, N_PDS_Comm * pdsCommPtr
  ) : op_      (op_in),
      allNodes_(allNodes_in),
      skipSet  (false),
      residualPtr_  (0),
      solutionPtr_  (0),
      rSize_   (0),
      pdsCommPtr_ (pdsCommPtr)
{
}

//-----------------------------------------------------------------------------
// Function      : AugmentLinSysOPStart::AugmentLinSysOPStart
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date :
//-----------------------------------------------------------------------------
AugmentLinSysOPStart::~AugmentLinSysOPStart()
{
}

//-----------------------------------------------------------------------------
// Function      : AugmentLinSysOPStart::augmentResidual
// Purpose       : augments the residual to support DCOP restart.
//
// Special Notes : erkeite:  This function doesn't seem to have much purpose.
//                 It saves pointers to the residual and solution vectors, so
//                 that it can use them in the augmentJacobian function, below.
//
//                 erkeite: The operations that need to be done on the
//                 residual depend on the jacobian structure.  In
//                 particular, they depend on what is happening with
//                 voltage sources.  Naively applying the 1 on the diagonal
//                 and 0 on the residual will result in singular matrices
//                 if the nodes in question are already being set by a
//                 voltage source.
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date :
//-----------------------------------------------------------------------------
void AugmentLinSysOPStart::augmentResidual
   (const Linear::Vector * solution, Linear::Vector * residual_vector)
{
  residualPtr_ = residual_vector;
  solutionPtr_ = solution;

#ifdef Xyce_PARALLEL_MPI
  pmap_ = residualPtr_->pmap();
#endif

  return;
}

//-----------------------------------------------------------------------------
// Function      : AugmentLinSysOPStart::augmentJacobian
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date :
//-----------------------------------------------------------------------------
void AugmentLinSysOPStart::augmentJacobian(Linear::Matrix * jacobian)
{
  IO::InitialConditionsData::NodeNamePairMap::iterator op_i;
  IO::InitialConditionsData::NodeNamePairMap::iterator op_end = op_.end();
  int i, row, rowLen, global_row, numRows;
  std::vector<int> col;
  std::vector<double> val;
  std::vector<double> resTmp;
  bool diag;
  int GID;
  std::map<int,std::string> rowOut;

  // erkeite: Determine which rows have diagonal elements already.
  // A lot of this code is needed to avoid causing the matrix to
  // be singular, which is easy to do when dealing with rows/cols
  // associated with Vsrc devices.
  numRows = jacobian->getLocalNumRows();
  if (!skipSet)
  {
    for (row=0 ; row<numRows ; ++row)
    {
#ifdef Xyce_PARALLEL_MPI
      global_row = pmap_->localToGlobalIndex(row);
#else
      global_row = row;
#endif
      rowLen = jacobian->getRowLength(global_row);
      //rowLen = jacobian->getLocalRowLength(row);
      col.resize(rowLen);
      val.resize(rowLen);
      jacobian->getRowCopy(global_row, rowLen, rowLen, &val[0], &col[0]);
      diag = false;
      for (i=0 ; i<rowLen ; ++i)
      {
        GID = col[i];
        if (GID == global_row)
        {
          if (val[i] != 0)
          {
            diag = true;
          }
        }
      }
      if (!diag)
      {
        skipLID.insert(row);
        skipGID.insert(global_row);
      }
    }
#ifdef Xyce_PARALLEL_MPI
    int numG = skipGID.size();
    int numG_tot;
    pdsCommPtr_->sumAll(&numG, &numG_tot, 1);
    int myPos;
    pdsCommPtr_->scanSum(&numG, &myPos, 1);
    myPos -= numG;
    std::vector<int> buf(numG_tot,0);
    std::vector<int> buf_g(numG_tot,0);
    std::set<int>::iterator skip_i=skipGID.begin();
    std::set<int>::iterator skip_end=skipGID.end();
    for ( ; skip_i != skip_end ; ++skip_i)
      buf[myPos++] = *skip_i;
    pdsCommPtr_->sumAll (&buf[0], &buf_g[0], numG_tot);
    skipGID.clear();
    for (i=0 ; i < numG_tot ; i++)
      skipGID.insert(buf_g[i]);
#endif
    skipSet = true;
  }

  op_i = op_.begin();
  std::set<int>::iterator skipLIDEnd = skipLID.end();
  std::set<int>::iterator skipGIDEnd = skipGID.end();
  for ( ; op_i != op_end ; ++op_i)
  {
    row = (*op_i).second.first;
    if (skipLID.find(row) == skipLIDEnd)
    {
#ifdef Xyce_PARALLEL_MPI
      global_row = pmap_->localToGlobalIndex(row);
#else
      global_row = row;
#endif
      rowLen = jacobian->getRowLength(global_row);
      col.resize(rowLen);
      val.resize(rowLen);
      jacobian->getRowCopy(global_row, rowLen, rowLen, &val[0], &col[0]);  // from the A-matrix
      // zero out all residual elements that are not in the skip list.
      Linear::Vector & residual = (*residualPtr_);
      residual[row] = 0;
      for (i=0 ; i<rowLen ; i++)
      {
        GID = col[i];
        // check if this column is in the skip list.  It will be in the skip list if
        // the diagonal element, (GID,GID) doensn't exist.
        if (skipGID.find(GID) == skipGIDEnd)
        {
          if (GID == global_row) // in other words, if this is the diagonal, then set to "1".
          {
            val[i] = 1;
          }
          else // otherwise, set all non-diagonal matrix entries to zero.
          {
            val[i] = 0;
          }
        }
        else // if (GID,GID) does not exist, then residual element is modified
        {
          residual[row] += val[i] *
            (const_cast<Linear::Vector*>(solutionPtr_))->getElementByGlobalIndex(GID);
        }
      }
      jacobian->putGlobalRow(global_row, rowLen, &val[0], &col[0]);
    }
  }

  return;
}

}}}

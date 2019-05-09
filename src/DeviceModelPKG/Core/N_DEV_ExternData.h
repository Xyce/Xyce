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

//-----------------------------------------------------------------------------
//
// Purpose        : This class is a container class for holding pointers,
//                  references, etc., to external linear solver data
//                  structures.  Occasionally, it will hold pointers
//                  to other things, but that is not the primary intention.
//
//                  In general, stuff that goes into this class should
//                  be stuff needed by more than one device instance type.
//
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_ExternData_h
#define Xyce_N_DEV_ExternData_h

#include <map>
#include <vector>

#include <N_DEV_fwd.h>
#include <N_LAS_fwd.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : ExternData
// Purpose       : Container for references, pointers, etc. to data structures
//                 outside the device package.  Mostly these are linear
//                 algebra objects like the Jacobian.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/25/02
//-----------------------------------------------------------------------------
class ExternData
{
public:
  ExternData()
    : lasSysPtr(0),

      JdxpVectorPtr(0),

      tmpdIdXPtr(0),
      tmpdQdXPtr(0),

      daeQVectorPtr(0),
      daeFVectorPtr(0),
      daeBVectorPtr(0),

      dFdxdVpVectorPtr(0),
      dQdxdVpVectorPtr(0),
      dQdxMatrixPtr(0),
      dFdxMatrixPtr(0),

      currSolVectorPtr(0),
      nextSolVectorPtr(0),
      lastSolVectorPtr(0),

      currStaVectorPtr(0),
      nextStaVectorPtr(0),
      lastStaVectorPtr(0),

      currStoVectorPtr(0),
      nextStoVectorPtr(0),
      
      nextLeadCurrFCompPtr(0),
      nextLeadCurrQCompPtr(0),
      nextJunctionVCompPtr(0),

      flagSolVectorPtr(0),

      nextStaDerivVectorPtr(0),

      deviceErrorWeightMask_(0),

      daeQVectorRawPtr(0),
      daeFVectorRawPtr(0),
      daeBVectorRawPtr(0),
      dFdxdVpVectorRawPtr(0),
      dQdxdVpVectorRawPtr(0),
      nextSolVectorRawPtr(0),
      currSolVectorRawPtr(0),
      lastSolVectorRawPtr(0),
      nextStaVectorRawPtr(0),
      currStaVectorRawPtr(0),
      lastStaVectorRawPtr(0),
      nextStoVectorRawPtr(0),
      currStoVectorRawPtr(0),
      nextStaDerivVectorRawPtr(0),
      initializeAllFlag(false)
  {}

  Linear::System * lasSysPtr;

  Linear::Vector * JdxpVectorPtr;

  Linear::Vector * tmpdIdXPtr;
  Linear::Vector * tmpdQdXPtr;

  // DAE formulation vectors
  Linear::Vector * daeQVectorPtr;
  Linear::Vector * daeFVectorPtr;
  Linear::Vector * daeBVectorPtr;

  Linear::Vector *  dFdxdVpVectorPtr;
  Linear::Vector *  dQdxdVpVectorPtr;

  // DAE formulation matrices
  Linear::Matrix * dQdxMatrixPtr;
  Linear::Matrix * dFdxMatrixPtr;

  Linear::Vector * currSolVectorPtr;
  Linear::Vector * nextSolVectorPtr;
  Linear::Vector * lastSolVectorPtr;

  Linear::Vector * currStaVectorPtr;
  Linear::Vector * nextStaVectorPtr;
  Linear::Vector * lastStaVectorPtr;

  Linear::Vector * currStoVectorPtr;
  Linear::Vector * nextStoVectorPtr;
  
  Linear::Vector * nextLeadCurrFCompPtr;
  Linear::Vector * nextLeadCurrQCompPtr;
  Linear::Vector * nextJunctionVCompPtr;
  
  Linear::Vector * flagSolVectorPtr;

  Linear::Vector  * nextStaDerivVectorPtr;

  Linear::Vector  * deviceErrorWeightMask_;

  // raw pointers (to internal vector data):
  double * daeQVectorRawPtr;
  double * daeFVectorRawPtr;
  double * daeBVectorRawPtr;
  double * dFdxdVpVectorRawPtr;
  double * dQdxdVpVectorRawPtr;

  double * nextSolVectorRawPtr;
  double * currSolVectorRawPtr;
  double * lastSolVectorRawPtr;

  double * nextStaVectorRawPtr;
  double * currStaVectorRawPtr;
  double * lastStaVectorRawPtr;

  double * nextStoVectorRawPtr;
  double * currStoVectorRawPtr;
  
  double * nextLeadCurrFCompRawPtr;
  double * nextLeadCurrQCompRawPtr;
  double * nextJunctionVCompRawPtr;

  double * nextStaDerivVectorRawPtr;

  bool initializeAllFlag;
};

} // namespace Device
} // namespace Xyce

#endif


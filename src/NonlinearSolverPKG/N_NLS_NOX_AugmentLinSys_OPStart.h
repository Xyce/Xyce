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
// Purpose        : Concrete class for augmenting the Jacobian for
//                  pseudo-transient solves.
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


#ifndef Xyce_N_NLS_NOX_AugmentLinSys_OPStart_h
#define Xyce_N_NLS_NOX_AugmentLinSys_OPStart_h

#include <set>

#include <N_UTL_fwd.h>

#include "N_NLS_NOX_AugmentLinSys.h"
#include "N_PDS_ParMap.h"
#include "N_PDS_Comm.h"

#include <N_IO_InitialConditions.h>

class Epetra_MapColoring;

//-----------------------------------------------------------------------------
// Class         : N_NLS::NOX::AugmentLinSysOPStart
// Purpose       :
// Creator       : Roger Pawlowski, SNL, 9233
// Creation Date :
//-----------------------------------------------------------------------------
namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

class AugmentLinSysOPStart : public AugmentLinSys {

public:
    AugmentLinSysOPStart(Xyce::IO::InitialConditionsData::NodeNamePairMap &, const Xyce::NodeNameMap &, N_PDS_Comm *);

  //! Dtor.
  virtual ~AugmentLinSysOPStart();

  void setProgressVariable(double dummy) {return;}

  void augmentResidual(const Xyce::Linear::Vector * solution, Xyce::Linear::Vector * residual_vector);

  void augmentJacobian(Xyce::Linear::Matrix * jacobian);

 private:

  //! map of specified variables
    Xyce::IO::InitialConditionsData::NodeNamePairMap & op_;
  const Xyce::NodeNameMap & allNodes_;

  bool skipSet;
  std::set<int> skipLID;
  std::set<int> skipGID;

  Xyce::Linear::Vector* residualPtr_;
  const Xyce::Linear::Vector* solutionPtr_;
  int rSize_;
  N_PDS_ParMap * pmap_;
  N_PDS_Comm * pdsCommPtr_;
};

}}}

#endif


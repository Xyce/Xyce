//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : Constructs Composite Transforms to be applied to
//                  LinearProblems
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 4/2/03
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_TransformTool_h
#define Xyce_N_LAS_TransformTool_h

#include <N_UTL_fwd.h>

#include <Teuchos_RCP.hpp>

#include <EpetraExt_Transform_Composite.h>

#include <Epetra_LinearProblem.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : Transform
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 4/2/03
//-----------------------------------------------------------------------------
struct Transform
: public EpetraExt::Transform_Composite<Epetra_LinearProblem>
{
};

//-----------------------------------------------------------------------------
// Class         : TransformTool
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 4/2/03
//-----------------------------------------------------------------------------
struct TransformTool
{
  // Construction of Composite Transform
  Teuchos::RCP<Transform> operator()( const Util::OptionBlock & options );
};

} // namespace Linear
} // namespace Xyce

#endif


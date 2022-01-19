//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Purpose        : Preconditioner Factory for Trilinos
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 10/01/07
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_TrilinosPrecondFactory_h
#define Xyce_N_LAS_TrilinosPrecondFactory_h

// ---------- Standard Includes ----------

#include <string>

#include <N_UTL_fwd.h>
#include <N_LAS_fwd.h>

// ----------   Xyce Includes   ----------

#include <N_LAS_PrecondFactory.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : PrecondFactory
// Purpose       :
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 10/01/07
//-----------------------------------------------------------------------------
class TrilinosPrecondFactory : public PrecondFactory
{
public:
  // Basic Constructor, sets preconditioner factory options.
  TrilinosPrecondFactory( const Util::OptionBlock & OB );

  // Destructor
  virtual ~TrilinosPrecondFactory() {}

  // Creates a new preconditioner (matrix based).
  Teuchos::RCP<Preconditioner> create( const Teuchos::RCP<Problem> & problem ) const;

  // Creates a new preconditioner (matrix free)
  // NOTE:  This creation type is not supported by the Trilinos factory at this time.
  Teuchos::RCP<Preconditioner> create( const Teuchos::RCP<System> & lasSystem ) const
  {
    Xyce::Report::DevelFatal0().in("TrilinosPrecondFactory::create()") << "using System is not supported!";
    return Teuchos::null;
  }

private:

  std::string precType_;
  Teuchos::RCP<const Util::OptionBlock> OB_;

  // Default constructor and copy constructor.
  TrilinosPrecondFactory();
  TrilinosPrecondFactory( const TrilinosPrecondFactory& pf );
};

} // namespace Linear
} // namespace Xyce

#endif // Xyce_N_LAS_TrilinosPrecondFactory_h

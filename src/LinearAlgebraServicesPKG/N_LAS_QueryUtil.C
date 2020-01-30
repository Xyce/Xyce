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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/26/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_LAS_QueryUtil.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : QueryUtil::QueryUtil
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter  SNL, Parallel Computational Sciences
// Creation Date : 5/25/05
//-----------------------------------------------------------------------------
QueryUtil::QueryUtil()
  : checkConnectivity_(true),
    supernode_(false),
    isClean_(false),
#ifdef Xyce_TEST_SOLN_VAR_MAP
    namesFile_(true)
#else
    namesFile_(false)
#endif
{
}

//-----------------------------------------------------------------------------
// Function      : QueryUtil::registerOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/05
//-----------------------------------------------------------------------------
bool QueryUtil::registerOptions(const Util::OptionBlock & OB)
{
  Util::ParamList::const_iterator it_tpL;
  Util::ParamList::const_iterator first = OB.begin();
  Util::ParamList::const_iterator last  = OB.end();

  for (it_tpL = first; it_tpL != last; ++it_tpL)
  {
    if (it_tpL->uTag()=="CHECK_CONNECTIVITY")
    {
      checkConnectivity_ = static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if(it_tpL->uTag()=="SUPERNODE")
    {
      supernode_ = static_cast<bool>(it_tpL->getImmutableValue<bool>());
    }
    else if(it_tpL->uTag()=="OUTPUTNAMESFILE")
    {
      namesFile_ = static_cast<bool>(it_tpL->getImmutableValue<bool>());
    }
  }

  return true;
}

void
QueryUtil::populateMetadata(
  IO::PkgOptionsMgr &   options_manager)
{
  Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("TOPOLOGY");

  parameters.insert(Util::ParamMap::value_type("REPLICATED_CKT", Util::Param("REPLICATED_CKT", 1)));
  parameters.insert(Util::ParamMap::value_type("CHECK_CONNECTIVITY", Util::Param("CHECK_CONNECTIVITY", 0)));
  parameters.insert(Util::ParamMap::value_type("SUPERNODE", Util::Param("SUPERNODE", false)));
  parameters.insert(Util::ParamMap::value_type("OUTPUTNAMESFILE", Util::Param("OUTPUTNAMESFILE", false)));
}

//-----------------------------------------------------------------------------
// Function      : QueryUtil::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool registerPkgOptionsMgr(QueryUtil &linear_solve_utility, IO::PkgOptionsMgr &options_manager)
{
  options_manager.addOptionsProcessor( "TIMEINT", IO::createRegistrationOptions(linear_solve_utility, &QueryUtil::registerTimeOptions));
  options_manager.addOptionsProcessor( "TOPOLOGY", IO::createRegistrationOptions(linear_solve_utility, &QueryUtil::registerOptions));

  return true;
}

} // namespace Linear
} // namespace Xyce

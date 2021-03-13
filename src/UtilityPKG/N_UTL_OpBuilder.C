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

//-----------------------------------------------------------------------------
//
// Purpose        : Provide tools for accessing output data in parallel or
//                  serial
//
// Special Notes  :
//
// Creator        : David Baur
//
// Creation Date  : 11/15/2013
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_UTL_OpBuilder.h>
#include <N_UTL_Param.h>

namespace Xyce {
namespace Util {
namespace Op {

//-----------------------------------------------------------------------------
// Function      : Builder::createOp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:56:45 2014
//-----------------------------------------------------------------------------
Operator * Builder::createOp(ParamList::const_iterator & it) const
{
  ParamList::const_iterator it2 = it;

  Operator *new_op = makeOp(it2);
  if (new_op)
    it = it2;

  return new_op;
}

//-----------------------------------------------------------------------------
// Function      : BuilderManager::createOp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 10:00:36 2014
//-----------------------------------------------------------------------------
Operator * BuilderManager::createOp(ParamList::const_iterator & it) const
{
  int index=0;
  for (BuilderVector::const_iterator it1 = opBuilderVector_.begin(), end1 = opBuilderVector_.end(); it1 != end1; ++it1, ++index) 
  {
#if 0
    if (index==1)
    {
      std::cout << "checkpoint for " << (*it).tag() << std::endl;
    }
#endif
    Operator *new_op = (*it1)->createOp(it);
    if (new_op)
    {
#if 0
      std::cout << "BuilderManager::createOp. index="<<index<<" succeeded in creating an op for " << (*it).tag() <<std::endl;
#endif
      return new_op;
    }
#if 0
    else
    {
      std::cout << "BuilderManager::createOp. index="<<index<<" failed    in creating an op for " << (*it).tag() <<std::endl;
    }
#endif
  }
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : BuilderManager::createOp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 10:00:36 2014
//-----------------------------------------------------------------------------
Operator *
BuilderManager::createOp(
  const std::string &   name) const
{
  ParamList param_list;
  param_list.push_back(Param(name, ""));
  ParamList::const_iterator it = param_list.begin();

  return createOp(it);
}

//-----------------------------------------------------------------------------
// Function      : BuilderManager::findCreateFunction
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 10:01:47 2014
//-----------------------------------------------------------------------------
CreateFunction
BuilderManager::findCreateFunction(
  Identifier    id) const
{
  CreateMap::const_iterator it = opCreateMap_.find(id);
  if (it != opCreateMap_.end())
    return (*it).second;

  return 0;
}

} // namespace Op
} // namespace Util
} // namespace Xyce

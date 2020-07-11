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

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 10/xx/2019
//
//
//
//
//-----------------------------------------------------------------------------

#include <iostream>
#include <unordered_map>
#include <string>

#include <xyceExpressionGroup.h>
#include <ast.h>
#include <newExpression.h>

namespace Xyce {
namespace Util {

//-------------------------------------------------------------------------------
// Function      : xyceExpressionGroup::xyceExpressionGroup
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 11/11/2019
//-------------------------------------------------------------------------------
xyceExpressionGroup::xyceExpressionGroup () :
  time_(0.0), temp_(0.0), VT_(0.0), freq_(0.0), dt_(0.0), alpha_(0.0)
{

}

//-------------------------------------------------------------------------------
// Function      : xyceExpressionGroup::getSolutionVal
// Purpose       : This was added to attempt to support the old API.  I hate it.
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : ???
//-------------------------------------------------------------------------------
bool xyceExpressionGroup::getSolutionVal(const std::string & nodeName, double & retval )
{
  bool success=true;
  retval = 0.0;
  std::string tmp = nodeName;
  Xyce::Util::toUpper(tmp);

  std::vector<std::string>::iterator it = std::find(names_.begin(), names_.end(), tmp);
  if (it != names_.end())
  {
    int index = it - names_.begin();
    int size = dvals_.size();
    if (index<size)
    {
      retval = dvals_[index];
    }

    std::cout << "Solution variable " << nodeName << " found by the xyceExpresionGroup! value = " << retval << std::endl;
  }
  else // not found
  {
    std::cout << "ERROR.  Solution variable " << nodeName << " not found by the xyceExpresionGroup!" << std::endl;
  }

  return success; // FIX THIS
}

//-------------------------------------------------------------------------------
// Function      : xyceExpressionGroup::getGlobalParameterVal
// Purpose       : This was added to attempt to support the old API.  I hate it.
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : ???
//-------------------------------------------------------------------------------
bool xyceExpressionGroup::getGlobalParameterVal (const std::string & nodeName, double & retval )
{
  bool success=true;
  retval = 0.0;
  std::string tmp = nodeName;
  Xyce::Util::toUpper(tmp);

  std::vector<std::string>::iterator it = std::find(names_.begin(), names_.end(), tmp);
  if (it != names_.end())
  {
    int index = it - names_.begin();
    retval = dvals_[index];
  }

  return success; // FIX THIS
}

//-------------------------------------------------------------------------------
// Function      : xyceExpressionGroup::getFunction
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 12/28/2019
//-------------------------------------------------------------------------------
bool xyceExpressionGroup::getFunction
  (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
{
#if 1
  std::cout << "xyceExpressionGroup::getFunction name = " << name <<std::endl;
#endif
  bool retval=true;

  std::string upperName = name;
  Xyce::Util::toUpper(upperName);

  if (functions_.find(upperName) != functions_.end()) { exp = functions_[upperName]; }
  else { retval = false; }

  return retval;
}

//-------------------------------------------------------------------------------
// Function      : xyceExpressionGroup::addFunction
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 12/28/2019
//-------------------------------------------------------------------------------
void xyceExpressionGroup::addFunction (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
{
#if 1
  std::cout << "xyceExpressionGroup::addFunction name = " << name <<std::endl;
#endif
  std::string upperName = name;
  Xyce::Util::toUpper(upperName);
  functions_[upperName] = exp;
}

//-------------------------------------------------------------------------------
// Function      : xyceExpressionGroup::addGlobalParam
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/4/2020
//-------------------------------------------------------------------------------
void xyceExpressionGroup::addGlobalParam (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
{
#if 1
  std::cout << "xyceExpressionGroup::addFunction name = " << name <<std::endl;
#endif
  std::string upperName = name;
  Xyce::Util::toUpper(upperName);
  globalParams_[upperName] = exp;
}

//-------------------------------------------------------------------------------
// Function      : xyceExpressionGroup::addParam
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/4/2020
//-------------------------------------------------------------------------------
void xyceExpressionGroup::addParam (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
{
#if 1
  std::cout << "xyceExpressionGroup::addFunction name = " << name <<std::endl;
#endif
  std::string upperName = name;
  Xyce::Util::toUpper(upperName);
  params_[upperName] = exp;
}

//-------------------------------------------------------------------------------
// Function      : xyceExpressionGroup::getParam
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 12/28/2019
//-------------------------------------------------------------------------------
bool xyceExpressionGroup::getParam
  (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
{
  bool retval=false;
  return retval;
}

//-------------------------------------------------------------------------------
// Function      : xyceExpressionGroup::getGlobalParam
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 12/28/2019
//-------------------------------------------------------------------------------
bool xyceExpressionGroup::getGlobalParam
  (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
{
  bool retval=false;
  return retval;
}

}
}

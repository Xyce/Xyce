
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

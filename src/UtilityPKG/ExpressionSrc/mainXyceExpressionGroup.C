
#include <iostream>
#include <unordered_map>
#include <string>

#include <mainXyceExpressionGroup.h>
#include <ast.h>
#include <newExpression.h>

namespace Xyce {
namespace Util {

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getSolutionVal
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 3/20/2020
//-------------------------------------------------------------------------------
bool mainXyceExpressionGroup::getSolutionVal(const std::string & nodeName, double & retval )
{
  bool success=true;
  retval = 0.0;
  std::string tmp = nodeName;
  Xyce::Util::toUpper(tmp);

  return success; // FIX THIS
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getFunction
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 3/20/2020
//-------------------------------------------------------------------------------
bool mainXyceExpressionGroup::getFunction
  (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
{
#if 1
  std::cout << "mainXyceExpressionGroup::getFunction name = " << name <<std::endl;
#endif
  bool retval=true;

  std::string upperName = name;
  Xyce::Util::toUpper(upperName);


  Util::ParamMap::const_iterator param_it = functions_.find(upperName);
  if (param_it != functions_.end())
  {
    const Util::Param & param = (*param_it).second;

    if ( param.getType()==Xyce::Util::EXPR )
    {
      Expression & expression = const_cast<Expression &>(param.getValue<Expression>());
      exp  = expression.newExpPtr_; 
    }
    else if (param.getType()== Xyce::Util::DBLE || param.getType() == Xyce::Util::STR) // ERK.  I want this option obsolete
    {
      double val = param.getImmutableValue<double>();
    }
  }
  else 
  { 
    retval = false; 
  }

  return retval;
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getParam
// Purpose       : 
// Special Notes :  
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 3/20/2020
//-------------------------------------------------------------------------------
bool mainXyceExpressionGroup::getParam
  (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
{
#if 1
  std::cout << "mainXyceExpressionGroup::getParam name = " << name <<std::endl;
#endif
  bool retval=true;

  std::string upperName = name;
  Xyce::Util::toUpper(upperName);


  Util::ParamMap::const_iterator param_it = params_.find(upperName);
  if (param_it != params_.end())
  {
    const Util::Param & param = (*param_it).second;

    if ( param.getType()==Xyce::Util::EXPR )
    {
      Expression & expression = const_cast<Expression &>(param.getValue<Expression>());
      exp  = expression.newExpPtr_; 
    }
    else if (param.getType()== Xyce::Util::DBLE || param.getType() == Xyce::Util::STR) // ERK.  I want this option obsolete
    {
      double val = param.getImmutableValue<double>();
    }
  }
  else 
  { 
    retval = false; 
  }

  return retval;
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getGlobalParam
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 3/20/2020
//-------------------------------------------------------------------------------
bool mainXyceExpressionGroup::getGlobalParam
  (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> & exp)
{
#if 1
  std::cout << "mainXyceExpressionGroup::getParam name = " << name <<std::endl;
#endif
  bool retval=true;

  std::string upperName = name;
  Xyce::Util::toUpper(upperName);

  Util::ParamMap::const_iterator param_it = globalParams_.find(upperName);
  if (param_it != globalParams_.end())
  {
    const Util::Param & param = (*param_it).second;

    if ( param.getType()==Xyce::Util::EXPR )
    {
      Expression & expression = const_cast<Expression &>(param.getValue<Expression>());
      exp  = expression.newExpPtr_; 
    }
    else if (param.getType()== Xyce::Util::DBLE || param.getType() == Xyce::Util::STR) // ERK.  I want this option obsolete
    {
      double val = param.getImmutableValue<double>();
    }
  }
  else 
  { 
    retval = false; 
  }

  return retval;
}

}
}
